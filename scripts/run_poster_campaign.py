"""Run the dedicated SI poster campaign without keeping Jupyter active.

Run from a merged WIMP-LEGEND checkout:
    python scripts/run_poster_campaign.py
    python scripts/run_poster_campaign.py --reuse-source --plot
"""
from pathlib import Path
import argparse
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analyses.legend_poster import load_config, resolve_profile, run_analysis, load_run
from analyses.legend_poster.campaign_scripts.poster_campaign import build_limits, load_spec


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reuse-source', action='store_true', help='Do not run detector physics; load saved source index')
    ap.add_argument('--plot', action='store_true', help='Also generate five poster figures')
    ap.add_argument('--ddp-root', type=Path, help='Path to local DirectDetectionPlots checkout (or use env)')
    ap.add_argument('--comparison-mode', choices=['ddp','frozen','none'], default='ddp',
                    help='Experimental limits: current DirectDetectionPlots, historical snapshots, or no limits')
    args = ap.parse_args()
    campaign = ROOT/'analyses/legend_poster/campaigns/poster_final'
    config_path = campaign/'config.json'
    spec_path = ROOT/'analyses/legend_poster/poster_campaign_spec.json'
    index_path = campaign/'source_run.json'
    config = load_config(config_path)
    spec = load_spec(spec_path)
    resolved = resolve_profile(config, 'standard')
    if (resolved['signal']['operators'] != [1] or
        resolved['analysis']['scenarios'] != [spec['input_scenario']] or
        resolved['signal'].get('si_rate_engine') != 'factorized_shm' or
        resolved['numerical'].get('folding_backend') != 'vectorized'):
        raise ValueError('Dedicated SI campaign configuration/optimizations not enabled')
    if args.reuse_source:
        state=json.loads(index_path.read_text(encoding='utf-8'))
        run=load_run(ROOT/state['run_relative'])
        if run['manifest']['run_id'] != state['run_id']:
            raise ValueError('Saved source index mismatch')
    else:
        t=time.perf_counter()
        run=run_analysis(config,profile='standard',
            cache_dir=ROOT/'analyses/legend_poster/cache',output_dir=campaign/'runs')
        campaign.mkdir(parents=True,exist_ok=True)
        index_path.write_text(json.dumps({'run_relative':str(run['path'].resolve().relative_to(ROOT)),
                    'run_id':run['manifest']['run_id']},indent=2)+'\n',encoding='utf-8')
        print(f'Numerical stage: {time.perf_counter()-t:.1f} seconds;'
              f" cache hits/misses {run['manifest']['cache']['hits']}/{run['manifest']['cache']['misses']}")
    t=time.perf_counter()
    rows,manifest=build_limits(run['path'],spec_path,campaign)
    print(f'Statistical stage: {time.perf_counter()-t:.1f} seconds;'
          f" {manifest['n_masses']} masses and {len(rows)} combinations")
    print('Saved:', campaign/'poster_limits.csv')
    if args.plot:
        from analyses.legend_poster.campaign_scripts.poster_figures import plot_five
        dest=plot_five(campaign,run['path'],spec_path,ROOT,ddp_root=args.ddp_root,
                       comparison_mode=args.comparison_mode)
        print('Figures:', dest)


if __name__=='__main__':
    main()
