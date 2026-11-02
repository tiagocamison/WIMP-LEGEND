"""Cold/warm complete O1 poster analysis, matched optional detector backends.

Example: python tools/benchmark_folding.py --output /tmp/fold --masses 24
Uses unchanged standard integration settings, all configured scenarios/exposures,
all five isotopes and historical comparisons. Each backend has a fresh process
and independent initially empty cache. No plotted/interpolated count comparison.
"""
from pathlib import Path
import argparse,cProfile,json,os,pstats,resource,subprocess,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from analyses.legend_poster import load_config,resolve_profile,run_analysis,load_run
from analyses.legend_poster.storage import source_identity,file_hash


def child(args):
    from threadpoolctl import threadpool_limits
    c=resolve_profile(load_config(),'standard')
    c['signal']['si_rate_engine']='factorized_shm';c['signal']['operators']=[1]
    if args.masses:c['signal']['masses_GeV']=np.geomspace(3.,100.,args.masses).tolist()
    c['numerical']['folding_backend']=args.backend
    out=args.output/args.backend;out.mkdir(parents=True,exist_ok=True)
    cache=out/'cache'
    if cache.exists() and any(cache.iterdir()):raise ValueError('Cold benchmark requires an empty cache; choose a new output directory')
    measure=[]
    with threadpool_limits(limits=c['numerical']['blas_threads']):
        for label in ('cold','warm'):
            profiler=cProfile.Profile() if args.profile_analysis and label=='cold' else None
            if profiler:profiler.enable()
            start=time.perf_counter();run=run_analysis(c,cache_dir=cache,output_dir=out/'runs');seconds=time.perf_counter()-start
            if profiler:
                profiler.disable();profiler.dump_stats(out/'analysis.prof')
                with (out/'profile.txt').open('w') as f:pstats.Stats(profiler,stream=f).sort_stats('cumulative').print_stats(50)
            measure.append(dict(kind=label,seconds=seconds,cache=run['manifest']['cache']))
            print(args.backend,label,seconds,'seconds',flush=True)
    data=dict(source=source_identity(),script_sha256=file_hash(__file__),configuration=c,
        measurements=measure,run_path=str(run['path'].resolve()),
        peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        compilation_seconds=0.,profiling_enabled=args.profile_analysis,
        note='Independent process peak RSS includes imports, cold and warm runs. Timed calls exclude imports. Profiling adds overhead if explicitly enabled; use separate output for profiling.')
    (out/'timing.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def run(args):
    args.output.mkdir(parents=True,exist_ok=True)
    records={}
    for backend in ('reference','vectorized'):
        command=[sys.executable,str(Path(__file__).resolve()),'--child','--backend',backend,
            '--output',str(args.output.resolve()),'--masses',str(args.masses)]
        if args.profile_analysis:command+=['--profile-analysis']
        start=time.perf_counter()
        subprocess.run(command,cwd=ROOT,check=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
        record=json.loads((args.output/backend/'timing.json').read_text())
        record['subprocess_seconds_including_imports']=time.perf_counter()-start;records[backend]=record
    a,b=[load_run(records[k]['run_path']) for k in ('reference','vectorized')]
    maximum={}
    for kind in ('signals','backgrounds'):
        assert a[kind].keys()==b[kind].keys()
        diffs=[]
        for key,ref in a[kind].items():
            value=b[kind][key]
            np.testing.assert_allclose(value,ref,rtol=3e-8,atol=max(np.max(abs(ref)),1e-300)*1e-13)
            active=abs(ref)>max(np.max(abs(ref)),1e-300)*1e-13
            diffs.append(float(np.max(abs(value[active]/ref[active]-1),initial=0)))
        maximum[kind+'_max_relative']=max(diffs)
    assert len(a['limits'])==len(b['limits'])
    diffs=[]
    for ref,value in zip(a['limits'],b['limits']):
        for key in ('scenario','mass_GeV','operator','exposure_kg_day','method','background','status','roi_lower_keV','roi_upper_keV'):
            assert value[key]==ref[key],(key,value[key],ref[key])
        for key in ('limit_value','mu_upper','reference_signal_events','background_events'):
            np.testing.assert_allclose(value[key],ref[key],rtol=3e-8,atol=0)
        diffs.append(abs(value['limit_value']/ref['limit_value']-1))
    maximum['limit_max_relative']=max(diffs)
    assert all(x['passed'] for r in (a,b) for x in r['convergence'])
    result=dict(records=records,agreement=maximum,limits=len(a['limits']),
        speedups={label:records['reference']['measurements'][i]['seconds']/records['vectorized']['measurements'][i]['seconds'] for i,label in enumerate(('cold','warm'))},
        acceptance='Same bins/ROI, configurations, rate engine and quadrature; count/limit rtol=3e-8; no changed convergence criteria.',
        run_manifests={k:file_hash(Path(v['run_path'])/'manifest.json') for k,v in records.items()})
    (args.output/'comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(speedups=result['speedups'],agreement=maximum,limits=result['limits'])),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--masses',type=int,default=24,help='0 retains standard seven masses; otherwise logspace 3..100 GeV')
    parser.add_argument('--profile-analysis',action='store_true');parser.add_argument('--child',action='store_true')
    parser.add_argument('--backend',choices=('reference','vectorized'),default='reference')
    args=parser.parse_args()
    if args.masses<0:parser.error('masses must be >=0')
    (child if args.child else run)(args)
