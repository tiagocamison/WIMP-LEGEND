"""Shared-input O1/O8/O11 comparison to RAPIDD's compiled C form factors.

Diagnostic only. Pointwise constant W fixtures are never production datasets.
The raw-velocity result follows difRateGen.c's unshifted halo-moment assembly;
substituting v_perp in that expression is a separately labelled diagnostic.
"""
import argparse
import ctypes
from dataclasses import replace
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from compare_responses import CHANNELS, PIN
from wimp_legend.interactions import WilsonCoefficients, RESPONSE_CHANNELS
from wimp_legend.rates import NREFTRateConfig, transition_probability
from wimp_legend.targets import (Isotope, SourceReference, MomentumVariable,
    PolynomialResponse, ResponseKey, NuclearResponseDataset,
    load_germanium_afh2014, germanium_afh2014_momentum)


def run(checkout):
    checkout = Path(checkout).resolve()
    revision = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
    if revision != PIN or subprocess.check_output(['git', '-C', str(checkout), 'status', '--porcelain'], text=True).strip():
        raise ValueError('Require clean pinned RAPIDD checkout')
    source = SourceReference('RAPIDD '+PIN, 'effFormFact.c, comparison-only pointwise fixture', 'implementation', notes='No production nuclear inputs replaced; missing Phi_tilde_prime unused by O1/O8/O11')
    result = []
    transform = np.array([[1., 1.], [1., -1.]])
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
        for A in (70, 72, 73, 74, 76):
            mass = 32*.938+(A-32)*.940  # explicit RAPIDD form-factor mass
            spin = 4.5 if A == 73 else 0.
            isotope = Isotope('Ge', A, 32, mass, spin, {k:source for k in ('identity','mass_GeV','spin')})
            base = load_germanium_afh2014(isotope, momentum=germanium_afh2014_momentum(A))
            for energy in (1., 10., 30.):
                q2 = 2*mass*energy*1e-6
                responses = {}
                for channel in RESPONSE_CHANNELS:
                    matrix = np.zeros((2, 2))
                    if channel in CHANNELS:
                        first, second = CHANNELS[channel]
                        for a, n1 in enumerate((b'p', b'n')):
                            for b, n2 in enumerate((b'p', b'n')):
                                matrix[a,b] = ff(A, first.encode(), second.encode(), n1, n2, energy)
                        if channel in ('Phi_double_prime_M','Delta_Sigma_prime'):
                            matrix = matrix.T
                        matrix = (2*spin+1)/(4*math.pi)*transform@matrix@transform.T
                    for t in (0,1):
                        for u in (0,1):
                            responses[ResponseKey(channel,t,u)] = PolynomialResponse((matrix[t,u],),0.,None,source)
                fixture = NuclearResponseDataset('diagnostic-only', isotope, replace(base.metadata,source=source,momentum=MomentumVariable('q2_GeV2')), responses)
                mchi = 40.
                mu = mchi*mass/(mchi+mass)
                vmin2 = q2/(4*mu**2)
                v2 = (700/299792.458)**2
                for op in (1,8,11):
                    for cp, cn in ((1.,1.),(1.,0.),(1.,-.7),(0.,1.)):
                        coefficients = WilsonCoefficients.from_proton_neutron({op:(cp,cn)})
                        ours = transition_probability(fixture,coefficients,q2_GeV2=q2,V=v2-vmin2,j_chi=.5,config=NREFTRateConfig(.938))
                        terms = []
                        for name in ('FormFact_v0','FormFact_v2'):
                            terms.append(sum(ca*cb*getattr(lib,name)(b'pn_Fitz',A,32,op,op,na,nb,energy,mchi,.5) for na,ca in ((b'p',cp),(b'n',cn)) for nb,cb in ((b'p',cp),(b'n',cn))))
                        raw = terms[0]+v2*terms[1]
                        shifted = terms[0]+(v2-vmin2)*terms[1]
                        result.append(dict(A=A,operator=op,ER_keV=energy,mchi_GeV=mchi,cp=cp,cn=cn,v_km_s=700.,ours=ours,rapidd_raw_velocity_kernel=raw,raw_relative_difference=raw/ours-1,shifted_velocity_diagnostic=shifted,shifted_relative_difference=shifted/ours-1))
    return dict(schema_version=1,rapidd_commit=revision,m_N_reference_GeV=.938,j_chi=.5,description='Shared C nuclear inputs, dimensionful pn coefficients; not full rates. Raw uses v^2; diagnostic shifted uses v_perp^2. No upstream or production changes.',records=result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rapidd_checkout',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = run(args.rapidd_checkout)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for op in (1,8,11):
        rows=[r for r in result['records'] if r['operator']==op]
        print(op,'raw max',max(abs(r['raw_relative_difference']) for r in rows),'shifted max',max(abs(r['shifted_relative_difference']) for r in rows))
