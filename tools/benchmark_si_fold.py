"""Matched cold detector-fold timing for generic and factorized SI, no new model."""
from pathlib import Path
from copy import deepcopy
import argparse,json,sys,tempfile,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from threadpoolctl import threadpool_limits
from analyses.legend_poster import load_config,resolve_profile
from analyses.legend_poster.models import target
from analyses.legend_poster.analysis import pure_signal_counts,observed_edges
from analyses.legend_poster.storage import source_identity,NumericalCache,file_hash


def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    c=resolve_profile(load_config(),'standard')
    isotope=next(i for i,_ in target(c,'natural') if i.spin>0)
    masses=c['signal']['masses_GeV'];mass=masses[len(masses)//2]
    identity=source_identity();measurements={};arrays={}
    for mode in ('generic','factorized_shm'):
        inputs=deepcopy(c);inputs['signal']['si_rate_engine']=mode
        values=[]
        for repeat in range(2):
            with tempfile.TemporaryDirectory(prefix='si-fold-') as directory:
                cache=NumericalCache(directory,identity['physics_fingerprint'])
                start=time.perf_counter()
                counts,check=pure_signal_counts(inputs,isotope,mass,1,c['detector']['default'],observed_edges(c),cache)
                elapsed=time.perf_counter()-start
                values.append(dict(seconds=elapsed,cache_hits=cache.hits,cache_misses=cache.misses,convergence=check))
                print(mode,'cold folded signal',repeat+1,elapsed,'seconds',flush=True)
                arrays[mode]=counts
        measurements[mode]=values
    reference=arrays['generic'];actual=arrays['factorized_shm']
    np.testing.assert_allclose(actual,reference,rtol=3e-8,atol=max(reference)*1e-13)
    data=dict(source=identity,benchmark_script_sha256=file_hash(__file__),configuration=c,
        isotope_A=isotope.A,mass_GeV=mass,
        method='Identical standard-profile pure_signal_counts; empty separate cache each repeat; includes construction, truth rates, interpolation refinements and detector integration.',
        measurements=measurements,
        median_speedup=float(np.median([r['seconds'] for r in measurements['generic']])/
                             np.median([r['seconds'] for r in measurements['factorized_shm']])),
        total_relative_difference=float(abs(actual.sum()/reference.sum()-1)))
    (out/'fold_benchmark.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print('Cold folded speedup',data['median_speedup'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    with threadpool_limits(limits=load_config()['numerical']['blas_threads']):run(args.output)
