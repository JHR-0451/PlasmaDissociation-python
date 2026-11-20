import numpy as np
import json
import pubchempy as pcp
import init
from scipy.spatial.distance import pdist, squareform
import os

def create_geom(n,nmod,T,modes,m,mom_num):
    Ax = modes[:,0]
    Ay = modes[:,1]
    Az = modes[:,2]
    Ax = Ax.reshape(n, nmod, order = 'F')
    Ay = Ay.reshape(n, nmod, order = 'F')
    Az = Az.reshape(n, nmod, order = 'F')
    rn = np.random.randn(nmod, mom_num)  # Use np.random.randn for standard normal distribution
    T=T*0.0000031668
    # Initialize arrays for random
    Meff = np.zeros(nmod)
    rv = np.zeros((nmod, mom_num))
    for i in range(nmod):
        for j in range(n):
            Meff[i] = Meff[i]+np.sum(((Ax[j, i]**2) + (Ay[j, i]**2) + (Az[j, i]**2)) * m[j])
        rv[i, :] = rn[i, :] * np.sqrt(2 * T / Meff[i])
    # Calculate the velocity by applying it through the tranformation matrix of normal modes.
    Vx = np.dot(Ax, rv)
    Vy = np.dot(Ay, rv)
    Vz = np.dot(Az, rv)
    Px = np.zeros((n,mom_num))
    Py = np.zeros((n,mom_num))
    Pz = np.zeros((n,mom_num))
    for i in range(n):
        Px[i,:] = Vx[i,:]*m[i]
        Py[i,:] = Vy[i,:]*m[i]
        Pz[i,:] = Vz[i,:]*m[i]
    
    return Px, Py, Pz


def create_wigner(n, nmod, T, modes, freqs_cm, m, mom_num, rng):
    """Thermal (harmonic) Wigner sampling of positions AND momenta, per normal mode.

    All atomic units (hbar = 1): masses in m_e (init.setup_masses), coordinates in bohr.
    For a mode with Cartesian displacement vector A (any normalisation), q is the coordinate along A:
    x = x_eq + A q,  Meff = sum_j m_j |A_j|^2,  KE = Meff qdot^2 / 2,  PE = Meff w^2 q^2 / 2.
    Thermal Wigner function of a harmonic oscillator = Gaussian in q and qdot with
        var(q)    = coth(w / 2kT) / (2 Meff w)
        var(qdot) = coth(w / 2kT) * w / (2 Meff)
    -> mean KE = mean PE = (w/4) coth(w/2kT); mean E = (w/2) coth(w/2kT); T = 0 gives zero-point sampling.
    Returns dX, dY, dZ (bohr displacements from the optimised geometry) and Px, Py, Pz, each (n, mom_num).
    """
    Ax = modes[:, 0].reshape(n, nmod, order='F')
    Ay = modes[:, 1].reshape(n, nmod, order='F')
    Az = modes[:, 2].reshape(n, nmod, order='F')
    w = np.asarray(freqs_cm, dtype=float) / 219474.6313632     # cm^-1 -> Hartree
    if np.any(w <= 0):
        raise ValueError(f"non-positive frequency in Wigner sampling: {freqs_cm}")
    kT = T * 0.0000031668                                        # same K -> Hartree factor as create_geom
    coth = 1.0 / np.tanh(w / (2 * kT)) if kT > 0 else np.ones_like(w)
    Meff = np.sum(m[:, None] * (Ax**2 + Ay**2 + Az**2), axis=0)  # (nmod,)
    sq = np.sqrt(coth / (2 * Meff * w))
    sv = np.sqrt(coth * w / (2 * Meff))
    q = rng.standard_normal((nmod, mom_num)) * sq[:, None]
    v = rng.standard_normal((nmod, mom_num)) * sv[:, None]
    dX, dY, dZ = Ax @ q, Ay @ q, Az @ q
    Px, Py, Pz = m[:, None] * (Ax @ v), m[:, None] * (Ay @ v), m[:, None] * (Az @ v)
    return dX, dY, dZ, Px, Py, Pz


def wigner_checks(dX, dY, dZ, Px, Py, Pz, m, T, freqs_cm, output_file="wigner_checks.txt"):
    """Write sample-mean KE and its expected value; harmonic PE needs the Hessian, so only KE is checked here."""
    w = np.asarray(freqs_cm, dtype=float) / 219474.6313632
    kT = T * 0.0000031668
    coth = 1.0 / np.tanh(w / (2 * kT)) if kT > 0 else np.ones_like(w)
    ke = 0.5 * np.sum((Px**2 + Py**2 + Pz**2) / m[:, None], axis=0)
    disp = np.sqrt(dX**2 + dY**2 + dZ**2)
    with open(output_file, "w") as f:
        f.write(f"Thermal Wigner sampling, T = {T} K, {Px.shape[1]} samples\n")
        f.write("Frequencies (cm^-1): " + " ".join(f"{x:.2f}" for x in freqs_cm) + "\n")
        f.write(f"Expected mean KE (sum w/4 coth): {np.sum(w / 4 * coth):.6f} Ha\n")
        f.write(f"Sample mean KE: {ke.mean():.6f} Ha (std {ke.std():.6f})\n")
        f.write(f"Expected mean vib. energy (sum w/2 coth): {np.sum(w / 2 * coth):.6f} Ha\n")
        f.write("RMS displacement per atom (bohr): " + " ".join(f"{x:.4f}" for x in np.sqrt((disp**2).mean(1))) + "\n")


def get_geometry_file():

    symbols = []
    coordinates = []
    with open("initial_guess.txt", "r") as file:
        for line in file:
            parts = line.split()
            if len(parts) == 4:  # Ensuring correct format
                symbols.append(parts[0])  # Atomic symbol
                x, y, z = map(float, parts[1:])
                coordinates.append([x * 1.88973, y * 1.88973, z * 1.88973])  # Convert to Bohr


    return symbols, coordinates

import numpy as np



def momenta_checks(Px, Py, Pz, symbols, temperature, output_file="momentum_checks.txt"):
    """
    Calculate the kinetic energy for each repeat, compare it to the expected value,
    and calculate the average magnitude of momentum for each atom across all repeats.
    Results are saved to a file.

    Parameters:
        Px, Py, Pz (np.ndarray): Momenta in x, y, z directions (shape: [natoms, repeats]).
        symbols (list): List of atomic symbols (e.g., ['H', 'C', 'O']).
        temperature (float): Temperature in Kelvin.
        output_file (str): File to save the results.
    """
    natoms = Px.shape[0]  # Number of atoms
    repeats = Px.shape[1]  # Number of repeats

    # Get the masses of the atoms from the ATOMIC_MASSES dictionary
    masses = init.setup_masses(symbols)

    # Initialize variables to store kinetic energy for each repeat
    kinetic_energies = np.zeros(repeats)

    # Loop through each repeat
    for j in range(repeats):
        file_kinetic_energy = 0.0  # Initialize kinetic energy for the current repeat
        for i in range(natoms):
            # Calculate the magnitude of the momentum vector
            magnitude = np.sqrt(Px[i, j]**2 + Py[i, j]**2 + Pz[i, j]**2)

            # Calculate the kinetic energy for this atom
            kinetic_energy = 0.5 * (magnitude**2 / masses[i])

            # Add the kinetic energy to the total for this repeat
            file_kinetic_energy += kinetic_energy

        # Store the total kinetic energy for this repeat
        kinetic_energies[j] = file_kinetic_energy

    # Calculate the average kinetic energy across all repeats
    average_kinetic_energy = np.mean(kinetic_energies)

    # Calculate the expected kinetic energy based on the equation
    # Expected KE = (3N - 6) * 0.5 * kT, where k = 3.16881e-6 (atomic units)
    expected_ke = (3 * natoms - 6) * 3.16881e-6 * temperature

    # Calculate the average magnitude of momentum for each atom across all repeats
    avg_momentum_magnitude = np.zeros(natoms)
    for i in range(natoms):
        total_momentum_magnitude = 0.0
        for j in range(repeats):
            total_momentum_magnitude += np.sqrt(Px[i, j]**2 + Py[i, j]**2 + Pz[i, j]**2)
        avg_momentum_magnitude[i] = total_momentum_magnitude / repeats

    # Write results to a file
    with open(output_file, "w") as file:
        file.write("Kinetic Energy Checks:\n")
        file.write(f"Expected Kinetic Energy: {expected_ke:.6f}\n")
        file.write(f"Average Kinetic Energy: {average_kinetic_energy:.6f}\n")
        file.write("Average Momentum Magnitude per Atom:\n")
        for i in range(natoms):
            file.write(f"Atom {i+1} ({symbols[i]}): {avg_momentum_magnitude[i]:.6f}\n")
        for j in range(repeats):
            file.write(f"Repeat {j+1}:\n")
            file.write(f"  Kinetic Energy: {kinetic_energies[j]:.6f}\n")
            if np.isclose(kinetic_energies[j], expected_ke, rtol=0.1*expected_ke):
                file.write("  Kinetic Energy check passed.\n")
            else:
                file.write("  Kinetic Energy check failed.\n")
        file.write("\n")






def generate_bondarr(symbols, coordinates):

    BOND_THRESHOLDS = {
        ("C", "C"): [(2.7401, "C=C"), (3.25, "C-C")],  # 1.45 Å and 1.70 Å * 1.88973 for C=C and C-C
        ("C", "F"): [(2.8346, "C-F")],  # 1.50 Å * 1.88973
        ("C", "H"): [(2.1732, "C-H")],  # 1.15 Å * 1.88973
        ("C", "O"): [(2.75, "C-O")],  # 1.43 Å * 1.88973
        ("O", "H"): [(2.0, "O-H")],  # 1.00 Å * 1.88973
        ("H", "O"): [(2.0, "O-H")],  # 1.00 Å * 1.88973
    }   

    natoms = len(symbols)

    # Compute pairwise distances
    distances = squareform(pdist(coordinates))  # Already in Bohr

    with open("../results/bondarr.txt", "w") as file:
        for i in range(natoms):
            for j in range(i + 1, natoms):
                elem1, elem2 = symbols[i], symbols[j]
                d = distances[i, j]  # Distance between atom i and j in Bohr
                
                # Sort elements alphabetically to match dictionary keys
                key = tuple(sorted((elem1, elem2)))
        
                # Check if bond exists in dictionary
                if key in BOND_THRESHOLDS:

                    # Loop through possible bond types for this pair
                    for max_dist, bond_type in BOND_THRESHOLDS[key]:
                        if d <= max_dist:
                            # Write bond type to file
                            file.write(f"{i+1}-{j+1}:{bond_type}\n")
                            
    
    print("bondarr.txt successfully written!")


def get_geometry_pubchem(molecule_name):
    search = pcp.get_compounds(molecule_name, 'name', record_type='3d')
    molecule = search[0]
    bonds = molecule.bonds
    atoms = molecule.atoms
    
    # Create a mapping from atom IDs to symbols
    atom_id_to_symbol = {atom.aid: atom.element for atom in atoms}
    symbols = [atom.element for atom in atoms]
    coordinates = [[atom.x*1.88973, atom.y*1.88973, atom.z*1.88973] for atom in atoms]
    with open("../results/bondarr.txt", "w") as file:
        for bond in bonds:
            atom1 = bond.aid1
            atom2 = bond.aid2
            order = bond.order
            
            # Determine the bond type based on the order
            if order == 1:
                bond_type = "-"
            elif order == 2:
                bond_type = "="
            elif order == 3:
                bond_type = "#"
            else:
                bond_type = "-"  # Default to single bond if unknown order
            
            # Get the symbols for the atoms
            atom1_symbol = atom_id_to_symbol.get(atom1, "Unknown")
            atom2_symbol = atom_id_to_symbol.get(atom2, "Unknown")
            
            # Write to file in the format 1-3:C-H
            file.write(f"{atom1}-{atom2}:{atom1_symbol}{bond_type}{atom2_symbol}\n")

    return symbols,coordinates

def organise_modes(modes,atoms):
    numeric_modes = np.zeros((len(modes)*atoms, 3))
    cnt=0
    for i in range(len(modes)):
        for j in range(atoms):
            numeric_modes[cnt,0]= float(modes[i]['displacement'][j][0])
            numeric_modes[cnt,1]= float(modes[i]['displacement'][j][1])
            numeric_modes[cnt,2]= float(modes[i]['displacement'][j][2])
            cnt+=1
   
    return numeric_modes

if __name__ == "__main__":
    with open('../inputs.json') as f:
        inputs=json.load(f)
    if inputs["Molecule_data"]["Geom_flg"] == "PubChem":
        symbols,coordinates = get_geometry_pubchem(inputs["Molecule_data"]["Molecule"])
    elif inputs["Molecule_data"]["Geom_flg"] == "NIST":
        from nist_geometry import get_geometry_nist
        symbols, coordinates = get_geometry_nist(inputs["Molecule_data"])
        generate_bondarr(symbols, coordinates)
    elif inputs["Molecule_data"]["Geom_flg"] == "Initial":
        symbols,coordinates = get_geometry_file()
        generate_bondarr(symbols, coordinates)
    if inputs["Elec_structure"]["method"] == "QChem":
        import qchem as qc
        opt_coords, modes = qc.initial_conditions(symbols,coordinates,inputs['HPC']['cores'], inputs['Elec_structure']['Basis'],
    opt_sym_ignore=inputs["Elec_structure"].get("Opt_sym_ignore", True),
    opt_sym_tol=inputs["Elec_structure"].get("Opt_sym_tol", 5),
    opt_geom_symflag=inputs["Elec_structure"].get("Opt_geom_symflag", True)
    )
    elif inputs["Elec_structure"]["method"] == "PySCF":
        import elec_pyscf as pyscf
        pyscf.initial_conditions(symbols,coordinates)
    
   
    modes = modes['modes']
    modes_raw = modes                     # keep the parsed dicts (frequencies) for Wigner sampling
    num_modes = len(modes)
    natoms = int((num_modes+6)/3)
    modes=organise_modes(modes,natoms)
    masses = init.setup_masses(symbols)
    # Initial-condition sampling. "Init_cond" absent or "classical" = original behaviour (unchanged code path).
    init_cond = inputs["Molecule_data"].get("Init_cond", "classical").lower()
    seed = inputs["Molecule_data"].get("Seed")              # optional int; absent = unseeded, as before
    T_in = inputs["Molecule_data"]["Temp"]
    nrep = inputs["HPC"]["repeats"]
    if init_cond == "classical":
        if seed is not None:
            np.random.seed(int(seed))
        Px, Py, Pz = create_geom(natoms,num_modes,T_in,modes,masses,nrep)
        momenta_checks(Px,Py,Pz,symbols,T_in)
        rep_coords = [opt_coords] * nrep
    elif init_cond == "wigner":
        freqs = [float(md['frequency']) for md in modes_raw]
        rng = np.random.default_rng(None if seed is None else int(seed))
        dX, dY, dZ, Px, Py, Pz = create_wigner(natoms,num_modes,T_in,modes,freqs,masses,nrep,rng)
        wigner_checks(dX, dY, dZ, Px, Py, Pz, masses, T_in, freqs)
        eq_lines = [ln.split() for ln in opt_coords.strip().splitlines()]
        eq_sym = [p[0] for p in eq_lines]
        eq_xyz = np.array([[float(x) for x in p[1:4]] for p in eq_lines])
        rep_coords = []
        for j in range(nrep):
            xyz = eq_xyz + np.stack([dX[:, j], dY[:, j], dZ[:, j]], axis=1)
            rep_coords.append("".join(f"{s} {float(x)!r} {float(y)!r} {float(z)!r} \n" for s, (x, y, z) in zip(eq_sym, xyz)))
    else:
        raise ValueError(f"Unknown Molecule_data.Init_cond: {init_cond!r} (use 'classical' or 'wigner')")
    for j in range(nrep):
        with open('../repetitions/rep-'+str(j+1)+'/Geometry', 'w') as file:
            file.write(rep_coords[j])
            file.write("momentum\n")
            # Write Px, Py, and Pz for each atom on the same line
            for atom in range(natoms):
                # Access the Px, Py, and Pz values using the corresponding indices
                px_value = Px[atom, j]
                py_value = Py[atom, j]
                pz_value = Pz[atom, j]
                file.write(f'{px_value}  {py_value}  {pz_value}\n')
        
