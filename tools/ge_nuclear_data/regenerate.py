"""GE-NUCLEAR-DATA-001: reproducible exact-decimal supplement extraction.

Usage: python tools/ge_nuclear_data/regenerate.py --supplement PATH
Requires development-only SymPy 1.14.0. Does not evaluate arbitrary Mathematica.
"""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re

import sympy as s

from afh_symbolic import element, y

SOURCE_SHA256 = '943a403e216633840887bbfcbec55c56e0c5813b86335b49868acb89e74a49b9'
CONVENTION_ID = 'anand2014_prc89_065501_half_isospin_v1'
ISOTOPES = (70, 72, 73, 74, 76)
# Source order in ResponseNuclear, lines 988–995. Mixed source operators are
# deliberately stored in their SOURCE order until the explicit transpose below.
PAIRS = (
    ('M', 'M', 'M'),
    ('Sigma_double_prime', 'Sigma_double_prime', 'Sigma_double_prime'),
    ('Sigma_prime', 'Sigma_prime', 'Sigma_prime'),
    ('Phi_double_prime', 'Phi_double_prime', 'Phi_double_prime'),
    ('Phi_tilde_prime', 'Phi_tilde_prime', 'Phi_tilde_prime'),
    ('Delta', 'Delta', 'Delta'),
    ('Phi_double_prime_M', 'M', 'Phi_double_prime'),
    ('Delta_Sigma_prime', 'Sigma_prime', 'Delta'),
)
OPERATORS = tuple(pair[0] for pair in PAIRS[:6])


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_source(path):
    if sha(path) != SOURCE_SHA256:
        raise ValueError('supplement SHA-256 does not match the pinned primary source')
    text = Path(path).read_text()
    tables = {A: [] for A in ISOTOPES}
    pattern = re.compile(r'^pfGe(70|72|73|74|76)\[(\d+)\]=\{[^=]+='
        r'(-?\d+\.\d+),(\d+),(\d+),(\d+),(\d+),(\d+),(\d+)\};$')
    for number, line in enumerate(text.splitlines(), 1):
        match = pattern.fullmatch(line)
        if match:
            A, index, psi, *quantum = match.groups()
            A, index = int(A), int(index)
            if index != len(tables[A])+1:
                raise ValueError('duplicate/missing/out-of-order Ge row')
            tables[A].append(dict(row=index, source_line=number, psi_decimal=psi,
                                  quantum=list(map(int, quantum))))
    for A in ISOTOPES:
        declared = int(re.search(rf'pfGe{A}length=(\d+);', text).group(1))
        if len(tables[A]) != declared or declared != (83 if A == 73 else 15):
            raise ValueError('incomplete Ge density matrix')
    spins = {A: s.Rational(re.search(rf'\{{32,{A},(\d+(?:/\d+)?)\}}', text).group(1))
             for A in ISOTOPES}
    return tables, spins


def amplitudes(A, records, spin):
    """Scaled exp(y)*sqrt(4pi) p/n reduced amplitudes, source FF lines 963–977."""
    T = s.Rational(A, 2)-32
    result = {}
    for operator in OPERATORS:
        multipoles = defaultdict(lambda: [s.S.Zero, s.S.Zero])
        for row in records:
            K, tau, Np, jp, N, j = row['quantum']
            if K > 2*spin:
                raise ValueError('density rank exceeds isotope angular momentum')
            sp = element(operator, Np, jp, N, j, K)
            if sp == 0:
                continue
            density = s.Rational(row['psi_decimal'])
            # Source's p/n form-factor calculation; no full-sum Wilson input.
            if tau == 0:
                p = n = s.sqrt(2)/s.sqrt(2*T+1)/2
            elif tau == 1:
                p = -s.sqrt(6*T)/s.sqrt((2*T+1)*(T+1))/2
                n = -p
            else:
                raise ValueError('unknown source isospin')
            multipoles[K][0] += density*p*sp
            multipoles[K][1] += density*n*sp
        result[operator] = {K: tuple(s.simplify(s.expand(a)) for a in pair)
                            for K, pair in sorted(multipoles.items())}
    return result


def coefficients(poly):
    poly = s.Poly(s.expand(poly), y)
    if any(term.has(s.Float) for term in poly.all_coeffs()):
        raise ValueError('inexact arithmetic leaked into symbolic extraction')
    degree = max(0, poly.degree())
    return ['0' if poly.nth(k) == 0 else str(s.N(poly.nth(k), 40))
            for k in range(degree+1)]


def export_isotope(A, rows, spin):
    amp = amplitudes(A, rows, spin)
    responses, source_pn = [], []
    transform = ((1, 1), (1, -1))  # cp=c0+c1, cn=c0-c1
    for source_index, (channel, first, second) in enumerate(PAIRS, 1):
        ranks = sorted(set(amp[first]) & set(amp[second]))
        pn = [[s.expand(sum(amp[first][K][t]*amp[second][K][u]
                            for K in ranks)/(2*spin+1))
               for u in (0, 1)] for t in (0, 1)]
        canonical = [[s.expand((2*spin+1)/(4*s.pi)*sum(
            transform[t][p]*pn[p][n]*transform[u][n]
            for p in (0, 1) for n in (0, 1))) for u in (0, 1)] for t in (0, 1)]
        if source_index in (7, 8):
            canonical = list(zip(*canonical))
        for t in (0, 1):
            for u in (0, 1):
                source_pn.append(dict(channel=channel, nucleon=t, nucleon_prime=u,
                    source_operator_order=[first, second], coefficients=coefficients(pn[t][u])))
                polynomial = canonical[t][u]
                responses.append(dict(channel=channel, tau=t, tau_prime=u,
                    source_response_index=source_index,
                    source_isospin_indices=[u, t] if source_index in (7, 8) else [t, u],
                    coefficients=coefficients(polynomial), exact_zero=bool(polynomial == 0),
                    verification_status='symbolic_source_reconstruction'))
    M0 = amp['M'][0]
    closure = [s.N(v.subs(y, 0)/s.sqrt(2*spin+1), 40) for v in M0]
    return {
        'schema_version': 1, 'dataset_version': 'afh2014_ge_obdme_v1',
        'isotope': {'symbol': 'Ge', 'Z': 32, 'A': A, 'spin': str(spin),
                    'mass_GeV': None, 'mass_policy': 'caller-supplied bare-nucleus mass with provenance'},
        'primary_source': {'doi': '10.1103/PhysRevC.89.065501',
            'supplement_sha256': SOURCE_SHA256, 'density_table': f'pfGe{A}',
            'density_row_count': len(rows), 'density_source_lines': [rows[0]['source_line'], rows[-1]['source_line']],
            'spin_source': 'IsotopicSpins, supplement line 947',
            'upstream': 'Fitzpatrick et al., JCAP 02 (2013) 004, arXiv:1203.3542v3; source comment names v2',
            'acquisition': 'user-supplied attachment 2026-10-19; original APS download date/account not supplied'},
        'supporting_materials': {
            'readme_install.txt': {'sha256': 'cc844dc64f5b0cd255d2e8a21e219281634ee1e7d064439211c5eeb99e145648',
                                   'role': 'primary supplement installation instructions, not numerical input'},
            'GE_DMFF_SUPPLEMENT_VS_ELASTIC_COMPARISON.md': {'sha256': '230582e778f9653c68a22cc10e7f57461f2d4ea6125ab846a86b8fd81a2ba18f',
                                   'role': 'user-supplied audit, not primary nuclear-structure data'}},
        'model': {'label': 'GCN2850', 'label_evidence': 'supplied audit and published form-factor reproduction; not named in supplement'},
        'convention_id': CONVENTION_ID,
        'representation': {'formula': 'exp(-2*y)*sum(coefficients[k]*y**k)', 'coefficient_order': 'ascending',
            'export_significant_digits': 40, 'density_precision': 'written source decimals retained verbatim as exact rationals; original truncated decimal precision is not improved',
            'runtime_precision': 'float64', 'units': 'dimensionless raw W; no external q^2/m_N^2 or nuclear spin factor absorbed'},
        'momentum': {'variable': 'y', 'definition': 'y=q2_GeV2*b_GeV_inv^2/4',
            'default_b_fm_prescription': 'sqrt(41.467/(45*A^(-1/3)-25*A^(-2/3)))',
            'prescription_source': 'supplement SetIsotope line 1479',
            'source_fm_in_GeV_inv': '5.0677',
            'source_prescription_fm_conversion': '5.0677 GeV^-1/fm, exactly the source rounded unit; caller must explicitly select this prescription or supply another sourced b'},
        'domain': {'mathematical_domain': [0, None], 'physically_validated_range': None,
                   'validated_range_policy': 'no source-backed numerical cutoff was supplied; not a claim of physical validity at arbitrary y'},
        'normalization': {'from_source_W': '4*W_source; transpose ordered isospin indices for source indices 7 and 8',
            'pn_to_canonical': 'W=(2J+1)/(4pi) S F_source S^T, S=[[1,1],[1,-1]]; transpose mixed channels',
            'particle_number_rescaling': False,
            'effective_protons_q0': str(closure[0]), 'effective_neutrons_q0': str(closure[1]),
            'effective_A_q0': str(closure[0]+closure[1])},
        'responses': responses, 'source_pn_polynomials': source_pn,
        'density_rows': rows,
    }


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=True)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supplement', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[2]/'src/wimp_legend/targets/data/ge_afh2014_v1')
    args = parser.parse_args()
    tables, spins = parse_source(args.supplement)
    args.output.mkdir(parents=True, exist_ok=True)
    files = {}
    for A in ISOTOPES:
        data = export_isotope(A, tables[A], spins[A])
        file = args.output/f'ge{A}.json'
        dump(file, data)
        files[str(A)] = {'file': file.name, 'sha256': sha(file)}
        print(A, data['normalization'], flush=True)
    dump(args.output/'manifest.json', {
        'schema_version': 1, 'dataset_version': 'afh2014_ge_obdme_v1',
        'supplement_sha256': SOURCE_SHA256, 'convention_id': CONVENTION_ID,
        'regeneration': {'sympy': s.__version__, 'exact_decimal_arithmetic': True,
            'scripts_sha256': {p.name: sha(p) for p in (Path(__file__), Path(__file__).with_name('afh_symbolic.py'))}},
        'isotopes': files,
    })


if __name__ == '__main__':
    main()
