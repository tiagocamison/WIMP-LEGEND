"""Run/replot the saved master configuration without a notebook."""
import argparse
from pathlib import Path
from .configuration import load_config,DEFAULT_CONFIG
from .analysis import run_analysis,load_run
from .plotting import plot_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=DEFAULT_CONFIG)
    parser.add_argument('--profile',default=None)
    parser.add_argument('--cache',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--plot-only',type=Path,help='Replot saved numerical outputs; no physics recomputation')
    args=parser.parse_args()
    config=load_config(args.config)
    run=load_run(args.plot_only) if args.plot_only else run_analysis(config,profile=args.profile,
        cache_dir=args.cache,output_dir=args.output)
    files=plot_run(run,plotting=config['plotting'])
    print('Numerical run:',run['path'])
    print('Cache:',run['manifest']['cache']['hits'],'hits;',run['manifest']['cache']['misses'],'misses')
    print('Saved',len(files),'figures and plot_manifest.json')


if __name__=='__main__':main()
