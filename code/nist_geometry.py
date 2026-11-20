"""Read experimental Cartesian geometries from NIST CCCBDB (no API required).

Uses only the Python standard library. Coordinates returned to setup are in Bohr,
using the existing project's conversion factor to keep this change isolated.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import re
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ANGSTROM_TO_BOHR = 1.88973
_ELEMENTS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split())


class NISTGeometryError(ValueError):
    """Missing, ambiguous, or invalid experimental geometry; never fall back."""


class _Tables(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables = []
        self.stack = []
        self.heading = []
        self.in_heading = False

    def handle_starttag(self, tag, attrs):
        if tag == 'h1':
            self.in_heading = True
        if tag == 'table':
            self.stack.append({'rows': [], 'row': None, 'cell': None})
        elif self.stack:
            table = self.stack[-1]
            if tag == 'tr':
                self._finish_row(table)
                table['row'] = []
            elif tag in ('td', 'th'):
                self._finish_cell(table)
                if table['row'] is None:
                    table['row'] = []
                table['cell'] = []

    @staticmethod
    def _finish_cell(table):
        if table['cell'] is not None:
            table['row'].append(' '.join(''.join(table['cell']).split()))
            table['cell'] = None

    def _finish_row(self, table):
        self._finish_cell(table)
        if table['row'] is not None:
            table['rows'].append(table['row'])
            table['row'] = None

    def handle_endtag(self, tag):
        if tag == 'h1':
            self.in_heading = False
        if not self.stack:
            return
        table = self.stack[-1]
        if tag in ('td', 'th'):
            self._finish_cell(table)
        elif tag == 'tr':
            self._finish_row(table)
        elif tag == 'table':
            self._finish_row(table)
            self.tables.append(self.stack.pop()['rows'])

    def handle_data(self, data):
        if self.in_heading:
            self.heading.append(data)
        if self.stack and self.stack[-1]['cell'] is not None:
            self.stack[-1]['cell'].append(data)


def validate_nist_request(cas_number, charge=0):
    """Return hyphen-free CAS, checking its checksum and workflow support."""
    digits = str(cas_number).strip().replace('-', '')
    if not re.fullmatch(r'\d{5,10}', digits):
        raise NISTGeometryError('NIST_CAS must be an explicit CAS number, e.g. 75-73-0.')
    checksum = sum(i * int(d) for i, d in enumerate(reversed(digits[:-1]), 1)) % 10
    if checksum != int(digits[-1]):
        raise NISTGeometryError('NIST_CAS has an invalid checksum.')
    # The existing initial_conditions workflow prepares neutral singlet molecules.
    if type(charge) is not int or charge != 0:
        raise NISTGeometryError('This geometry-source integration supports NIST_charge=0 only; ion initialisation needs separate code changes.')
    return digits


def _formula_counts(formula):
    if not isinstance(formula, str) or not re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', formula):
        raise NISTGeometryError('For NIST, Molecule must be a plain molecular formula, e.g. CF4; isomers are selected by NIST_CAS.')
    counts = Counter()
    for element, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
        if element not in _ELEMENTS or (count and int(count) < 1):
            raise NISTGeometryError('Invalid element or atom count in molecular formula.')
        counts[element] += int(count or 1)
    return counts


def parse_nist_geometry(html, expected_formula):
    """Extract exactly one Atom/x(Å)/y(Å)/z(Å) table and validate identity."""
    expected = _formula_counts(expected_formula)
    parser = _Tables()
    parser.feed(html)
    heading = ' '.join(''.join(parser.heading).split())
    match = re.search(r'Experimental data for\s+([^\s(]+)', heading)
    if not match or _formula_counts(match.group(1)) != expected:
        raise NISTGeometryError('NIST page heading does not match the requested molecular formula.')
    candidates = []
    for rows in parser.tables:
        if rows and [re.sub(r'\s+', '', c).lower() for c in rows[0]] == ['atom', 'x(å)', 'y(å)', 'z(å)']:
            candidates.append(rows[1:])
    if len(candidates) != 1:
        raise NISTGeometryError('Expected exactly one experimental Cartesian table in Å; none or multiple found.')
    symbols, coordinates, labels = [], [], []
    for row in candidates[0]:
        if len(row) != 4:
            raise NISTGeometryError('Malformed Cartesian row in NIST response.')
        atom = re.fullmatch(r'([A-Z][a-z]?)(\d+)', row[0])
        if not atom or atom.group(1) not in _ELEMENTS or row[0] in labels:
            raise NISTGeometryError('Invalid or duplicate NIST Cartesian atom label.')
        try:
            xyz = [float(v.replace('−', '-')) for v in row[1:]]
        except ValueError as exc:
            raise NISTGeometryError('Non-numeric NIST Cartesian coordinate.') from exc
        if not all(math.isfinite(v) for v in xyz):
            raise NISTGeometryError('Non-finite NIST Cartesian coordinate.')
        symbols.append(atom.group(1))
        coordinates.append(xyz)
        labels.append(row[0])
    if Counter(symbols) != expected:
        raise NISTGeometryError('NIST Cartesian atom counts do not match Molecule.')
    for i in range(len(coordinates)):
        for j in range(i):
            if math.dist(coordinates[i], coordinates[j]) < 1e-6:
                raise NISTGeometryError('Coincident atoms in NIST geometry.')
    return symbols, coordinates, heading


def get_geometry_nist(molecule_data, output_dir='.'):
    """Download and archive experimental coordinates; return symbols and Bohr."""
    cas = validate_nist_request(molecule_data.get('NIST_CAS', ''), molecule_data.get('NIST_charge', 0))
    formula = molecule_data['Molecule']
    _formula_counts(formula)
    url = 'https://cccbdb.nist.gov/exp2x.asp?' + urlencode({'casno': cas, 'charge': 0})
    request = Request(url, headers={
        'User-Agent': 'Mozilla/5.0 (compatible; PlasmaDissociation-python NIST geometry reader)',
        'Accept': 'text/html',
    })
    try:
        with urlopen(request, timeout=30) as response:
            if response.status != 200:
                raise NISTGeometryError(f'NIST returned HTTP {response.status}.')
            raw = response.read(5_000_001)
            encoding = response.headers.get_content_charset() or 'utf-8'
    except (HTTPError, URLError, TimeoutError, OSError) as exc:
        raise NISTGeometryError(f'Cannot retrieve NIST geometry from {url}: {exc}. No PubChem fallback was used.') from exc
    if len(raw) > 5_000_000:
        raise NISTGeometryError('NIST response exceeds the expected page size.')
    try:
        html = raw.decode(encoding)
    except (UnicodeError, LookupError) as exc:
        raise NISTGeometryError('Cannot decode NIST HTML response.') from exc
    symbols, angstrom, heading = parse_nist_geometry(html, formula)
    bohr = [[value * ANGSTROM_TO_BOHR for value in xyz] for xyz in angstrom]
    provenance = {
        'source': 'NIST CCCBDB experimental geometry', 'url': url,
        'retrieved_utc': datetime.now(timezone.utc).isoformat(),
        'html_sha256': hashlib.sha256(raw).hexdigest(),
        'cas_number_digits': cas, 'charge': 0, 'molecular_formula': formula,
        'page_heading': heading, 'symbols': symbols,
        'source_units': 'angstrom', 'coordinates_angstrom': angstrom,
        'output_units': 'bohr', 'angstrom_to_bohr': ANGSTROM_TO_BOHR,
        'coordinates_bohr': bohr,
    }
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'nist_source.html').write_bytes(raw)
    (folder / 'nist_geometry.json').write_text(json.dumps(provenance, indent=2) + '\n')
    (folder / 'nist_geometry_angstrom.xyz').write_text(
        str(len(symbols)) + '\nNIST CCCBDB experimental Cartesian geometry; angstrom\n' +
        ''.join(f'{symbol} {xyz[0]:.10f} {xyz[1]:.10f} {xyz[2]:.10f}\n' for symbol, xyz in zip(symbols, angstrom))
    )
    print(f'NIST geometry: {formula}, CAS {cas}, {len(symbols)} atoms; saved source and coordinate provenance.')
    return symbols, bohr


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description='Fetch and validate NIST geometry only; does not run Q-Chem or submit jobs.')
    cli.add_argument('--cas', required=True)
    cli.add_argument('--formula', required=True)
    cli.add_argument('--output-dir', required=True)
    args = cli.parse_args()
    atoms, xyz = get_geometry_nist({'NIST_CAS': args.cas, 'NIST_charge': 0, 'Molecule': args.formula}, args.output_dir)
    for atom, row in zip(atoms, xyz):
        print(atom, *(f'{v:.10f}' for v in row), '(Bohr)')
