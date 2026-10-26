"""Reproducible SI speed/memory and frozen-generic folded-limit comparison.

Run from the repository root: python tools/benchmark_si_fast.py --output DIR
All physics and quadrature inputs come from the saved poster master config.
The immutable published standard snapshot is the independent generic fold oracle.
"""
from pathlib import Path
from copy import deepcopy
import argparse,json,platform,sys,time,tracemalloc,resource
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from threadpoolctl import threadpool_limits
from analyses.legend_poster import load_config,resolve_profile,load_run
from analyses.legend_poster.models import target,isotope_weights
from analyses.legend_poster.analysis import signal_model,pure_signal_counts,observed_edges,count_limit
from analyses.legend_poster.storage import NumericalCache,source_identity,file_hash


def benchmark(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    c=resolve_profile(load_config(),'standard');fast=deepcopy(c)
    fast['signal']['si_rate_engine']='factorized_shm'
    old=load_run(ROOT/'analyses/legend_poster/validation/standard')
    prior=deepcopy(old['manifest']['resolved_configuration'])
    prior['signal']['si_rate_engine']='generic'
    assert prior==c,'Reference physics changed: do not compare unlike configurations'
    identity=source_identity();cache=NumericalCache(out/'cache',identity['physics_fingerprint'])
    bins=observed_edges(c);pure={};comparison=[]
    t=time.perf_counter()
    for mass in c['signal']['masses_GeV']:
        for isotope,fraction in target(c,'natural'):
            values,check=pure_signal_counts(fast,isotope,mass,1,c['detector']['default'],bins,cache)
            key=f'pure_Ge{isotope.A}__O1__m{mass:g}__benchmark'
            reference=old['signals'][key]
            np.testing.assert_allclose(values,reference,rtol=3e-8,atol=np.max(abs(reference))*1e-13)
            pure[isotope.A,mass]=values
            active=reference>max(reference)*1e-12
            comparison.append(dict(A=isotope.A,mass_GeV=mass,
                total_relative=float(abs(values.sum()/reference.sum()-1)),
                significant_bin_relative=float(np.max(abs(values[active]/reference[active]-1))),
                convergence=check))
    folded_seconds=time.perf_counter()-t
    print('Validated',len(comparison),'pure-isotope folds in',folded_seconds,'seconds;',cache.hits,'cache hits;',cache.misses,'misses',flush=True)
    limits=[]
    for row in old['limits']:
        if row['operator']!=1:continue
        scenario=row['scenario'];record=c['background']['scenarios'][scenario]
        mixed=sum((w['bulk_mass_fraction']*pure[a,row['mass_GeV']]
                   for a,w in isotope_weights(c,record['target']).items()),np.zeros(len(bins)-1))
        key=scenario+('__historical__Total' if row['background']=='historical' else '__Total')
        detector=c['detector']['presets'][record.get('detector',c['detector']['default'])]
        _,actual=count_limit(mixed,old['backgrounds'][key],bins,row['exposure_kg_day'],
                             detector['threshold_keV_obs'],c,row['method'])
        sigma=actual['mu_upper']*c['signal']['sigma_p_ref_cm2']
        np.testing.assert_allclose(sigma,row['limit_value'],rtol=3e-8,atol=0)
        assert actual['roi_lower_keV']==row['roi_lower_keV'] and actual['roi_upper_keV']==row['roi_upper_keV']
        limits.append(dict(scenario=scenario,background=row['background'],mass_GeV=row['mass_GeV'],
             exposure_kg_day=row['exposure_kg_day'],relative=float(abs(sigma/row['limit_value']-1))))
    print('Validated',len(limits),'physical SI limits',flush=True)
    # Same prepared model construction and same vectors for both implementations.
    masses=np.geomspace(min(c['signal']['masses_GeV']),max(c['signal']['masses_GeV']),12)
    prepared={}
    for mode,cfg in [('generic',c),('factorized_shm',fast)]:
        prepared[mode]=[(signal_model(cfg,i,float(m),1)) for m in masses for i,_ in target(c,'natural')]
    timings={};outputs={}
    for mode,models in prepared.items():
        samples=[]
        for repeat in range(2):
            t=time.perf_counter()
            values=[model(endpoint*np.linspace(0,.99,24)) for model,endpoint in models]
            samples.append(time.perf_counter()-t)
            print(mode,'raw scan',repeat+1,':',samples[-1],'seconds',flush=True)
        timings[mode]=samples;outputs[mode]=values
    for a,b in zip(outputs['generic'],outputs['factorized_shm']):
        np.testing.assert_allclose(a,b,rtol=3e-9,atol=np.max(abs(a))*2e-13)
    memory={}
    for mode,models in prepared.items():
        # Explicit tracked Python/NumPy allocation peak, not a native-RSS delta.
        model,endpoint=models[2]
        tracemalloc.start();model(endpoint*np.linspace(0,.99,24))
        current,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
        memory[mode]={'tracked_peak_bytes':peak,'tracked_current_bytes':current}
    result=dict(source=identity,configuration=c,benchmark_script_sha256=file_hash(__file__),
        reference_commit=old['manifest']['source']['git_commit'],
        reference_manifest_sha256=file_hash(old['path']/'manifest.json'),
        comparison_tolerances=dict(rate_rtol=3e-9,rate_atol_fraction=2e-13,fold_rtol=3e-8,fold_atol_fraction=1e-13,limit_rtol=3e-8),
        folded_comparison=comparison,limit_comparison=limits,fast_folded_35_signals_seconds=folded_seconds,
        folded_cache=dict(hits=cache.hits,misses=cache.misses),
        timing=dict(masses=12,isotopes=5,recoil_energies=24,repeats=2,
                    construction_included=False,seconds=timings,
                    median_speedup=float(np.median(timings['generic'])/np.median(timings['factorized_shm']))),
        memory=dict(method='tracemalloc per 24-point call; Python/NumPy tracked allocations, excludes untracked native memory',
                    measurements=memory,whole_benchmark_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        python=platform.python_version(),numpy=np.__version__)
    (out/'benchmark.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Passed',len(comparison),'pure-isotope folded comparisons and',len(limits),'SI limit comparisons')
    print(json.dumps(result['timing'],indent=2));print(json.dumps(result['memory'],indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    with threadpool_limits(limits=load_config()['numerical']['blas_threads']):benchmark(args.output)
