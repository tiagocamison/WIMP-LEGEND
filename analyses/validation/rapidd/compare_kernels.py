"""Shared-input comparison to RAPIDD's compiled C form factors.

Diagnostic only. Pointwise constant W fixtures are never production datasets.
The raw-velocity result follows difRateGen.c's unshifted halo-moment assembly;
substituting v_perp in that expression is a separately labelled diagnostic.
"""
import argparse
import ctypes
from dataclasses import replace
from itertools import product
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from compare_responses import CHANNELS, PHI_MAPPING, PHI_TERMS, PIN
from wimp_legend.interactions import CANONICAL_OPERATORS, WilsonCoefficients, RESPONSE_CHANNELS
from wimp_legend.rates import NREFTRateConfig, transition_probability
from wimp_legend.targets import (Isotope, SourceReference, MomentumVariable,
    PolynomialResponse, ResponseKey, NuclearResponseDataset,
    load_germanium_afh2014, germanium_afh2014_momentum)


HISTORICAL_OPERATORS = (1, 8, 11)
GE_ISOTOPES = (70, 72, 73, 74, 76)
COUPLINGS = ((1., 1.), (1., 0.), (1., -.7), (0., 1.))
C_KM_S = 299792.458


def shared_fixture(ff, isotope, energy, source, *, map_phi):
    transform = np.array([[1., 1.], [1., -1.]])
    base = load_germanium_afh2014(isotope, momentum=germanium_afh2014_momentum(isotope.A))
    responses = {}
    for channel in RESPONSE_CHANNELS:
        matrix = np.zeros((2, 2))
        pairs = (CHANNELS[channel],) if channel in CHANNELS else ()
        if channel == 'Phi_tilde_prime' and map_phi:
            pairs = PHI_TERMS
        for a, n1 in enumerate((b'p', b'n')):
            for b, n2 in enumerate((b'p', b'n')):
                matrix[a, b] = sum(
                    ff(isotope.A, first.encode(), second.encode(), n1, n2, energy)
                    for first, second in pairs
                )
        if channel in ('Phi_double_prime_M', 'Delta_Sigma_prime'):
            matrix = matrix.T
        matrix = (2*isotope.spin+1)/(4*math.pi)*transform@matrix@transform.T
        for t, u in product((0, 1), repeat=2):
            responses[ResponseKey(channel, t, u)] = PolynomialResponse(
                (matrix[t, u],), 0., None, source
            )
    metadata = replace(base.metadata, source=source, momentum=MomentumVariable('q2_GeV2'))
    return NuclearResponseDataset('diagnostic-only', isotope, metadata, responses)


def relative_difference(actual, reference):
    return None if reference == 0 else actual/reference-1


def run(checkout, *, operators=HISTORICAL_OPERATORS, isotopes=GE_ISOTOPES,
        energies_keV=(1., 10., 30.), masses_GeV=(40.,), speeds_km_s=(700.,),
        couplings=COUPLINGS):
    """Run diagonal kernels on kinematically allowed points of an editable grid."""
    checkout = Path(checkout).resolve()
    revision = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
    if revision != PIN or subprocess.check_output(['git', '-C', str(checkout), 'status', '--porcelain'], text=True).strip():
        raise ValueError('Require clean pinned RAPIDD checkout')
    operators, isotopes = tuple(operators), tuple(isotopes)
    energies_keV, masses_GeV, speeds_km_s = map(tuple, (energies_keV, masses_GeV, speeds_km_s))
    couplings = tuple(tuple(pair) for pair in couplings)
    if not operators or not set(operators) <= set(CANONICAL_OPERATORS):
        raise ValueError('Select operators from the canonical basis')
    if not isotopes or not set(isotopes) <= set(GE_ISOTOPES):
        raise ValueError('Select stable Ge isotope mass numbers')
    for name, values in (('energies_keV', energies_keV), ('masses_GeV', masses_GeV),
                         ('speeds_km_s', speeds_km_s)):
        if not values or any(not math.isfinite(v) or v <= 0 for v in values):
            raise ValueError(f'{name} must contain finite positive values')
    if not couplings or any(len(pair) != 2 or not all(map(math.isfinite, pair)) for pair in couplings):
        raise ValueError('Supply finite (cp, cn) coupling pairs')
    if max(speeds_km_s) >= C_KM_S:
        raise ValueError('Speeds must be less than c')
    map_phi = bool(set(operators) & {12, 13})
    source = SourceReference(
        'RAPIDD '+PIN, 'effFormFact.c, comparison-only pointwise fixture', 'implementation',
        notes=PHI_MAPPING if map_phi else 'Phi_tilde_prime unused for selected operators',
    )
    result = []
    skipped = []
    with tempfile.TemporaryDirectory() as temporary:
        library = Path(temporary)/'responses.so'
        subprocess.run(['gcc', '-shared', '-fPIC', '-O2', str(checkout/'lib/source/effFormFact.c'), '-lm', '-o', str(library)], check=True)
        lib = ctypes.CDLL(str(library))
        ff = lib.FF_pn_Fitz
        ff.argtypes = [ctypes.c_int]+[ctypes.c_char_p]*4+[ctypes.c_double]
        ff.restype = ctypes.c_double
        for name in ('FormFact_v0', 'FormFact_v2'):
            function = getattr(lib, name)
            function.argtypes = [ctypes.c_char_p]+[ctypes.c_int]*4+[ctypes.c_char_p]*2+[ctypes.c_double]*3
            function.restype = ctypes.c_double
        for A in isotopes:
            mass = 32*.938+(A-32)*.940  # explicit RAPIDD form-factor mass
            spin = 4.5 if A == 73 else 0.
            isotope = Isotope('Ge', A, 32, mass, spin, {k:source for k in ('identity','mass_GeV','spin')})
            for energy in energies_keV:
                q2 = 2*mass*energy*1e-6
                fixture = shared_fixture(ff, isotope, energy, source, map_phi=map_phi)
                for mchi, speed in product(masses_GeV, speeds_km_s):
                    mu = mchi*mass/(mchi+mass)
                    vmin2 = q2/(4*mu**2)
                    v2 = (speed/C_KM_S)**2
                    if v2 < vmin2:
                        skipped.append(dict(A=A, ER_keV=energy, mchi_GeV=mchi, v_km_s=speed))
                        continue
                    for op, (cp, cn) in product(operators, couplings):
                        coefficients = WilsonCoefficients.from_proton_neutron({op:(cp,cn)})
                        ours = transition_probability(fixture,coefficients,q2_GeV2=q2,V=v2-vmin2,j_chi=.5,config=NREFTRateConfig(.938))
                        terms = []
                        for name in ('FormFact_v0','FormFact_v2'):
                            terms.append(sum(ca*cb*getattr(lib,name)(b'pn_Fitz',A,32,op,op,na,nb,energy,mchi,.5) for na,ca in ((b'p',cp),(b'n',cn)) for nb,cb in ((b'p',cp),(b'n',cn))))
                        raw = terms[0]+v2*terms[1]
                        shifted = terms[0]+(v2-vmin2)*terms[1]
                        result.append(dict(
                            A=A, operator=op, ER_keV=energy, mchi_GeV=mchi, cp=cp, cn=cn,
                            v_km_s=speed, ours=ours, rapidd_raw_velocity_kernel=raw,
                            raw_relative_difference=relative_difference(raw, ours),
                            shifted_velocity_diagnostic=shifted,
                            shifted_relative_difference=relative_difference(shifted, ours),
                        ))
    report = dict(schema_version=1,rapidd_commit=revision,m_N_reference_GeV=.938,j_chi=.5,description='Shared C nuclear inputs, dimensionful pn coefficients; not full rates. Raw uses v^2; diagnostic shifted uses v_perp^2. No upstream or production changes.',records=result)
    if map_phi:
        report['phi_mapping'] = PHI_MAPPING
    if skipped:
        report['skipped_kinematic_points'] = skipped
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rapidd_checkout',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--all-operators', action='store_true')
    args = parser.parse_args()
    operators = CANONICAL_OPERATORS if args.all_operators else HISTORICAL_OPERATORS
    result = run(args.rapidd_checkout, operators=operators)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for op in operators:
        rows=[r for r in result['records'] if r['operator']==op]
        maxima = [max((abs(r[key]) for r in rows if r[key] is not None), default=None)
                  for key in ('raw_relative_difference', 'shifted_relative_difference')]
        zeros = sum(r['ours'] == r['rapidd_raw_velocity_kernel'] == r['shifted_velocity_diagnostic'] == 0
                    for r in rows)
        print(f'O{op}: raw max {maxima[0]}, shifted max {maxima[1]}, common zeros {zeros}')
