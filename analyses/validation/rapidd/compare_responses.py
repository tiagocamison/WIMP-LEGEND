"""Run the pinned, unmodified RAPIDD C form factors; no vendored coefficients.

This compares shared nuclear inputs, NOT independent nuclear-structure theory.
No production WIMP-LEGEND defaults or responses are changed.
"""
import argparse
import ctypes
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile

PIN = 'c36b290abfc8f68dc4e53445bc18c185f4b466b0'
CHANNELS = {
    'M': ('M', 'M'), 'Sigma_prime': ('Sp', 'Sp'),
    'Sigma_double_prime': ('Spp', 'Spp'),
    'Phi_double_prime': ('Ppp', 'Ppp'), 'Delta': ('D', 'D'),
    # Source p/n matrices are stored in the source's operator order,
    # before our canonical mixed-channel transposition.
    'Phi_double_prime_M': ('M', 'Ppp'),
    'Delta_Sigma_prime': ('Sp', 'D'),
}
PHI_TERMS = (('Pp', 'Pp'), ('Pt', 'Pt'))
PHI_MAPPING = 'RAPIDD Pp/Pp + Pt/Pt -> AFH Phi_tilde_prime (diagnostic)'
ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / 'src/wimp_legend/targets/data/ge_afh2014_v1'


def polynomial(coefficients, y):
    value = 0.
    for coefficient in reversed(coefficients):
        value = value*y + float(coefficient)
    return math.exp(-2*y)*value


def run(checkout, *, include_phi_candidate=False):
    checkout = Path(checkout).resolve()
    revision = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
    if revision != PIN:
        raise ValueError(f'Expected RAPIDD {PIN}; got {revision}')
    if subprocess.check_output(['git', '-C', str(checkout), 'status', '--porcelain'], text=True).strip():
        raise ValueError('RAPIDD checkout must be clean (including untracked files)')
    source = checkout/'lib/source/effFormFact.c'
    text = source.read_text()
    b2 = [float(x) for x in re.search(r'static const double B2\[199\] = \{([^}]+)', text)[1].split(',')]
    records = []
    with tempfile.TemporaryDirectory() as temporary:
        library = Path(temporary)/'responses.so'
        subprocess.run(['gcc', '-shared', '-fPIC', '-O2', str(source), '-lm', '-o', str(library)], check=True)
        response = ctypes.CDLL(str(library)).FF_pn_Fitz
        response.argtypes = [ctypes.c_int] + [ctypes.c_char_p]*4 + [ctypes.c_double]
        response.restype = ctypes.c_double
        for A in (70, 72, 73, 74, 76):
            data = json.loads((DATA/f'ge{A}.json').read_text())
            for row in data['source_pn_polynomials']:
                channel = row['channel']
                if channel in CHANNELS:
                    pairs = (CHANNELS[channel],)
                elif channel == 'Phi_tilde_prime' and include_phi_candidate:
                    pairs = PHI_TERMS
                else:
                    continue  # Preserve the scope of the first-stage artifact.
                n1, n2 = row['nucleon'], row['nucleon_prime']
                for y in (0., .001, .01, .1):
                    # Compare at SAME y, not the same recoil energy with differing
                    # nuclear masses/oscillator prescriptions hidden inside C.
                    energy = 2*y/(b2[A-1]*(32*.938+(A-32)*.940))*1e6
                    actual = sum(
                        response(A, first.encode(), second.encode(),
                                 b'p' if n1 == 0 else b'n', b'p' if n2 == 0 else b'n', energy)
                        for first, second in pairs
                    )
                    expected = polynomial(row['coefficients'], y)
                    records.append(dict(A=A, channel=channel, n1=n1, n2=n2,
                                        y=y, rapidd_ER_keV=energy, rapidd=actual,
                                        afh_source=expected, difference=actual-expected,
                                        relative_difference=None if expected == 0 else actual/expected-1))
    report = dict(schema_version=1, rapidd_commit=revision,
                rapidd_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                afh_manifest_sha256=hashlib.sha256((DATA/'manifest.json').read_bytes()).hexdigest(),
                comparison='source-order p/n F=e^(-2y)P(y), matched y; no canonical W rescaling',
                unsupported=['Phi_tilde_prime'],
                interpretation='Input comparison only; NOT a complete rate/detector/statistics validation',
                records=records)
    if include_phi_candidate:
        report['unsupported'] = []
        report['phi_mapping'] = PHI_MAPPING
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rapidd_checkout', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--include-phi-candidate', action='store_true')
    args = parser.parse_args()
    result = run(args.rapidd_checkout, include_phi_candidate=args.include_phi_candidate)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(f"Saved {len(result['records'])} comparisons to {args.output}")
