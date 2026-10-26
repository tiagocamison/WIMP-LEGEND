"""Reproduce historical backgrounds and plot BG-BENCHMARK-001 configurations.

Run from any directory: python examples/background_benchmark.py --output DIR
The original notebook is read, never edited. No likelihood or sensitivity run.
"""
import argparse,csv,hashlib,json,sys
from dataclasses import asdict
from pathlib import Path
from scipy.integrate import cumulative_trapezoid
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools/background_benchmark'))
from legacy import load_legacy
from folding import quadrature_counts
from wimp_legend.backgrounds.benchmark import BenchmarkConfig,build_benchmark
from wimp_legend.detector import ElectronRecoilResponse,DetectorResponse,GaussianResolution,LindhardQuenching,ConstantEfficiency
from wimp_legend.energy import EnergyDomain
from wimp_legend._analysis import scalar


def run(configuration,output,notebook):
    config=json.loads(Path(configuration).read_text());p=config['plotting'];output=Path(output);output.mkdir(parents=True,exist_ok=True)
    for field in ('min_energy_keV','threshold_keV','noise_sigma_keV','fano_factor','pair_energy_keV'):
        scalar(p[field],field)
    for field in ('bin_width_keV','exposure_kg_day','truth_step_keV','gaussian_tail_sigma'):
        scalar(p[field],field,positive=True)
    if not p['min_energy_keV']<=p['threshold_keV']<p['max_energy_keV']:
        raise ValueError('require plot minimum <= analysis threshold < maximum')
    legacy,notebook_sha=load_legacy(notebook)
    if notebook_sha!=config['legacy_notebook_sha256']:raise ValueError('legacy notebook checksum differs from locked configuration')
    sigma=lambda e:np.sqrt(p['noise_sigma_keV']**2+p['fano_factor']*p['pair_energy_keV']*np.asarray(e))
    efficiency=ConstantEfficiency(p['constant_detection_efficiency'])
    er=ElectronRecoilResponse(GaussianResolution(sigma),observed_efficiency=efficiency)
    nr=DetectorResponse(LindhardQuenching(k=p['lindhard_k'],Z=32),GaussianResolution(sigma),observed_efficiency=efficiency)
    # Threshold is an analysis boundary, not a second application of efficiency.
    bins=np.unique(np.r_[np.arange(p['min_energy_keV'],p['max_energy_keV'],p['bin_width_keV']),p['threshold_keV'],p['max_energy_keV']])
    bins=bins[(bins>=p['min_energy_keV'])&(bins<=p['max_energy_keV'])]
    centers=(bins[:-1]+bins[1:])/2;exposure=p['exposure_kg_day']
    curves={};records=[];summary={}
    for name,scenario in config['scenarios'].items():
        print('Folding',name,flush=True)
        model=build_benchmark(BenchmarkConfig(**scenario),er_response=er,nr_response=nr)
        def bound(component):
            upper=component.spectrum_energy_range_keV[1]
            upper_mean=upper if component.energy_domain is EnergyDomain.ER else nr.quenching.mean_energy_keV(upper)
            return float(sigma(upper_mean))
        rates={};counts={}
        for component,value in model.components.items():
            counts[component]=quadrature_counts(value,bins,exposure,step_keV=p['truth_step_keV'],
                order=p['quadrature_order'],tail_sigma=p['gaussian_tail_sigma'],sigma_bound=bound)
            rates[component]=counts[component]/exposure/np.diff(bins)*365.25
        rates['Total']=sum(rates.values());curves[name]=rates
        summary[name]=dict(assumptions=asdict(model.config),physical_provenance=json.loads(model.provenance),
            analysis_window_keV=[p['threshold_keV'],p['max_energy_keV']],
            expected_counts_in_analysis_window={k:float(v[bins[:-1]>=p['threshold_keV']].sum()) for k,v in counts.items()})
        for component,rate in rates.items():
            for lo,hi,value in zip(bins[:-1],bins[1:],rate):records.append([name,component,lo,hi,value/365.25])
    old={}
    for name,target,argon in [('natural_AAr','natGe','AAr'),('enriched_UAr','enrGe','UAr')]:
        dense=np.unique(np.r_[np.linspace(bins[0],bins[-1],10001),bins])
        values=legacy(dense,target=target,argon=argon)
        old[name]={}
        for component,density in values.items():
            integral=cumulative_trapezoid(density,dense,initial=0)
            old[name][component]=np.diff(np.interp(bins,dense,integral))/np.diff(bins)
        for component,rate in old[name].items():
            for lo,hi,value in zip(bins[:-1],bins[1:],rate):records.append(['legacy_'+name,component,lo,hi,float(value)/365.25])
    aliases={'Flat_Continuum':'flat_observed','Tritium':'H3','39Ar':'Ar39','68Ge':'Ge68',
             '65Zn':'Zn65','55Fe':'Fe55','CEvNS':'solar_CEvNS','nu_e_ES':'solar_ES'}
    colors={'H3':'purple','Ge68':'goldenrod','Ga68':'darkorange','Zn65':'royalblue',
            'Fe55':'forestgreen','Ar39':'darkturquoise','solar_CEvNS':'yellowgreen',
            'solar_ES':'seagreen','flat_observed':'firebrick','Total':'black'}
    labels={'H3':r'$^3$H','Ge68':r'$^{68}$Ge','Ga68':r'$^{68}$Ga','Zn65':r'$^{65}$Zn',
            'Fe55':r'$^{55}$Fe','Ar39':r'$^{39}$Ar','solar_CEvNS':r'solar CE$\nu$NS',
            'solar_ES':r'solar $\nu$-e ES','flat_observed':'flat (observed)','Total':'Total'}
    def plot(ax,values,title):
        for name,rate in values.items():
            canonical=aliases.get(name,name)
            ax.stairs(rate,bins,label=labels[canonical],linewidth=1.8 if name=='Total' else 1.,color=colors[canonical])
        ax.set(xlim=(bins[0],bins[-1]),yscale='log',ylim=(1e-5,1e4),title=title,xlabel='Observed energy [keV]')
        ax.axvline(p['threshold_keV'],ls=':',color='gray',label='analysis threshold')
        ax.grid(alpha=.2);ax.legend(fontsize=7,ncol=2)
    fig,axes=plt.subplots(2,2,figsize=(14,10),sharex=True)
    plot(axes[0,0],old['natural_AAr'],'Historical notebook: natural Ge + AAr')
    plot(axes[0,1],old['enriched_UAr'],'Historical notebook: enriched Ge + UAr')
    plot(axes[1,0],curves['natural_AAr'],'Improved benchmark: natural Ge + AAr')
    plot(axes[1,1],curves['enriched_UAr'],'Improved benchmark: enriched Ge + UAr')
    for row in axes:
        for ax in row:ax.set_ylabel('Rate [events / (kg year keV observed)]')
    fig.suptitle('BG-BENCHMARK-001 — illustrative benchmark, not a validated LEGEND prediction')
    fig.tight_layout();fig.savefig(output/'background_comparison.png',dpi=160);fig.savefig(output/'background_comparison.pdf');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(14,5.7))
    plot(axs[0],curves['legend_poster'],'Improved LEGEND poster benchmark')
    plot(axs[1],curves['legend_poster'],'K/L/M+ capture and threshold region')
    axs[1].set_xlim(bins[0],2.)
    for ax in axs:ax.set_ylabel('Rate [events / (kg year keV observed)]')
    fig.suptitle('BENCHMARK: activity history, LAr transfer, veto and response are assumptions',fontsize=10)
    fig.tight_layout();fig.savefig(output/'legend_background_benchmark.png',dpi=180);fig.savefig(output/'legend_background_benchmark.pdf');plt.close(fig)
    with (output/'component_rates.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['scenario','component','bin_lower_keV_obs','bin_upper_keV_obs','events_per_kg_day_keV_obs']);w.writerows(records)
    metadata=dict(status=config['status'],configuration_sha256=hashlib.sha256(Path(configuration).read_bytes()).hexdigest(),legacy_notebook_sha256=notebook_sha,
        quadrature=dict(order=p['quadrature_order'],truth_step_keV=p['truth_step_keV'],tail_sigma=p['gaussian_tail_sigma'],description='bin probabilities from existing responses; exact point-line folding; fixed GL continua; explicit Gaussian tail cutoff; legacy observed spectra bin-averaged on10001-point grid without refolding'),scenarios=summary)
    (output/'benchmark_run.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return metadata


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=ROOT/'examples/bg_benchmark_config.json')
    parser.add_argument('--output',type=Path,default=ROOT/'benchmark-output')
    parser.add_argument('--notebook',type=Path,default=ROOT/'WIMP_O1_spectrum.ipynb')
    args=parser.parse_args();run(args.config,args.output,args.notebook)
