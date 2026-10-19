"""Replay the pinned source in Mathics3, without using the Python operator port.

Development only: SymPy 1.14.0, Mathics3 10.0.1. Supply the original attachment.
The single-particle reference uses the source operator expressions directly.
The source NineJ sum is shortened from 50 to 5: its SixJ triangle (1,l,g)
requires g<=l+1<=5 for every Ge orbit. No nonzero summands are removed.
Mathics uses machine arithmetic here; this is not a commercial Mathematica
certification or an independent nuclear-structure calculation.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
from mathics.core.load_builtin import import_and_load_builtins

from regenerate import parse_source, PAIRS, SOURCE_SHA256


NAMES = dict(M='MJ', Sigma_prime='SigmaPJ', Sigma_double_prime='SigmaPPJ',
             Delta='DeltaJ', Phi_double_prime='PhiPPJ', Phi_tilde_prime='PhiTPJ')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supplement', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows, spins = parse_source(args.supplement)
    import_and_load_builtins()
    from mathics.session import MathicsSession
    session = MathicsSession()
    session.evaluate('Get['+json.dumps(str(args.supplement.resolve()))+']')
    session.evaluate('dmformfactor`Private`cutsum=5;')
    coordinates = ('1/4', '1')
    cache, reference = {}, []
    for A, records in rows.items():
        T, J = A/2-32, float(spins[A])
        amplitudes = [{op: defaultdict(lambda: np.zeros(2)) for op in NAMES}
                      for _ in coordinates]
        for row in records:
            K, tau, Np, jp, N, j = row['quantum']
            p = math.sqrt(2)/math.sqrt(2*T+1)/2 if tau == 0 else -math.sqrt(6*T/((2*T+1)*(T+1)))/2
            nucleons = np.array([p, p if tau == 0 else -p])*float(row['psi_decimal'])
            for op, name in NAMES.items():
                for index, coordinate in enumerate(coordinates):
                    key = (op, Np, jp, N, j, K, coordinate)
                    if key not in cache:
                        expression = (f'N[Exp[{coordinate}]*Sqrt[4 Pi]*dmformfactor`Private`{name}'
                                      f'[{coordinate},{{{Np},{jp}/2}},{{{N},{j}/2}},{K}],17]')
                        value = session.evaluate(expression).to_python()
                        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
                            raise RuntimeError(f'unresolved source expression: {expression}: {value}')
                        cache[key] = value
                    amplitudes[index][op][K] += cache[key]*nucleons
        S = np.array([[1, 1], [1, -1]])
        for index, coordinate in enumerate(coordinates):
            values = []
            for source_index, (channel, left, right) in enumerate(PAIRS, 1):
                pn = sum((np.outer(amplitudes[index][left][K], amplitudes[index][right][K])
                          for K in amplitudes[index][left].keys() & amplitudes[index][right].keys()), np.zeros((2, 2)))/(2*J+1)
                W = (2*J+1)/(4*math.pi)*S@pn@S.T*math.exp(-2*float(Fraction(coordinate)))
                if source_index in (7, 8):
                    W = W.T
                values.extend(dict(channel=channel, tau=t, tau_prime=u, value=float(W[t, u]))
                              for t in (0, 1) for u in (0, 1))
            reference.append(dict(A=A, y=0.25 if index == 0 else 1.0, responses=values))
        print(f'Ge{A}: {len(cache)} cached source matrix elements', flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(
        source_sha256=SOURCE_SHA256, interpreter='Mathics3 10.0.1',
        scope='original single-particle expressions evaluated at two positive y; Python p/n assembly and explicit canonical map',
        limitations='shares SymPy angular backend; not commercial Mathematica or independent shell-model validation',
        cutsum=5, triangle_bound='g<=l+1<=5',
        single_particle=[dict(operator=k[0], quantum=list(k[1:6]), y=0.25 if k[6]=='1/4' else 1.0, value=float(v)) for k,v in sorted(cache.items())],
        datasets=reference), indent=2)+'\n')


if __name__ == '__main__':
    main()
