"""Publication-style figures from saved numerical outputs, independent of physics."""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from .storage import file_hash,digest
from .analysis import load_run

DATA=Path(__file__).parent/'data/si_limits'


def comparison_curves(selected):
    manifest=json.loads((DATA/'manifest.json').read_text())
    lookup={r['id']:r for r in manifest['curves']}
    if set(selected)-set(lookup):raise ValueError('unknown comparison dataset')
    result=[]
    for name in selected:
        record=lookup[name];path=DATA/record['file']
        if file_hash(path)!=record['sha256']:raise ValueError('comparison source checksum mismatch')
        xy=np.loadtxt(path)
        if xy.ndim!=2 or xy.shape[1]!=2 or np.any(xy<=0) or np.any(~np.isfinite(xy)):
            raise ValueError('comparison curve must contain positive mass/cross-section pairs')
        # Preserve source vertex order; projections/contours need not be sorted.
        result.append((record,xy))
    return manifest,result


def limit_series(run,operator):
    """The same source of plotted ordinates used by tests/notebooks; no recalculation."""
    groups={}
    for row in run['limits']:
        if row['operator']!=operator:continue
        key=(row['scenario'],row['background'],row['exposure_kg_day'],row['method'])
        groups.setdefault(key,[]).append(row)
    return {key:sorted(values,key=lambda r:r['mass_GeV']) for key,values in groups.items()}


def plot_run(run,*,plotting=None,output_dir=None):
    run=load_run(run['path'] if isinstance(run,dict) else run)
    c=run['manifest']['resolved_configuration'];style=plotting or c['plotting']
    out=Path(output_dir or run['path']/'figures');out.mkdir(parents=True,exist_ok=True)
    bins=run['backgrounds']['observed_edges_keV'];width=np.diff(bins)
    year=c['analysis']['days_per_year'];files=[]
    compare,curves=comparison_curves(style['comparisons'])
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.9,
                         'xtick.direction':'in','ytick.direction':'in','savefig.bbox':'tight'}):
        for scenario in c['analysis']['scenarios']:
            fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
            for ax,limits in zip(axes,(style['energy_range_keV'],style['threshold_zoom_keV'])):
                prefix=scenario+'__'
                colors=iter(plt.get_cmap('tab10').colors)
                for key,values in run['backgrounds'].items():
                    if not key.startswith(prefix) or '__historical__' in key:continue
                    label=key[len(prefix):]
                    ax.stairs(values/width*year,bins,label=label,
                              color='black' if label=='Total' else next(colors),
                              linewidth=1.8 if label=='Total' else .9)
                preset=c['background']['scenarios'][scenario].get('detector',c['detector']['default'])
                ax.axvline(c['detector']['presets'][preset]['threshold_keV_obs'],ls=':',color='gray')
                if not np.any(run['backgrounds'][scenario+'__Total']>0):ax.set_ylim(1e-12,1)
                ax.set(xlim=limits,yscale='log',xlabel='Observed energy [keV]',
                       ylabel='Events / (kg year keV observed)')
                if not np.any(run['backgrounds'][scenario+'__Total']>0):
                    ax.set_ylim(1e-12,1);ax.text(.05,.5,'No active background',transform=ax.transAxes)
                else:
                    total=run['backgrounds'][scenario+'__Total']/width*year
                    visible=(bins[:-1]<limits[1])&(bins[1:]>limits[0])&(total>0)
                    if np.any(visible):
                        ax.set_ylim(float(total[visible].max())*1e-6,float(total[visible].max())*10)
                ax.grid(alpha=.2);ax.legend(fontsize=7,ncol=2)
            fig.suptitle(style['title']+' — '+scenario+'\n'+c['execution']['resolved_label'])
            name='background_'+scenario
            for ext in ('pdf','png'):
                path=out/(name+'.'+ext);fig.savefig(path,dpi=style['dpi']);files.append(path)
            plt.close(fig)
        for operator in c['signal']['operators']:
            fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
            if operator==1:
                for record,xy in curves:
                    projected=record['kind']=='projection'
                    label=record['id'].replace('_',' ') + (' projection' if projected else '')
                    label+=f" ({record['confidence_level']:.0%} CL)"
                    ax.plot(xy[:,0],xy[:,1],ls='--' if projected else '-',alpha=.7,lw=.9,label=label)
            selected_exposures={c['analysis']['exposures_kg_day'][i] for i in style['limit_exposure_indices']}
            for key,rows in limit_series(run,operator).items():
                if key[2] not in selected_exposures:continue
                x=np.array([r['mass_GeV'] for r in rows]);y=np.array([r['limit_value'] for r in rows])
                valid=np.isfinite(y)&(y>0)
                scenario,bg,exposure,method=key
                label=f'{scenario}, {exposure/year:g} kg yr, {bg}, {method}'
                ax.plot(x[valid],y[valid],marker='o',ms=3,lw=1.2,
                        ls=':' if bg=='historical' else '-',label=label)
            ax.set(xscale='log',yscale='log',xlabel=r'$m_\chi$ [GeV]',xlim=style['mass_range_GeV'])
            if operator==1:
                ax.set(ylabel=r'$\sigma_{p,\mathrm{upper}}$ [cm$^2$]',ylim=style['sigma_range_cm2'])
            else:
                ax.set_ylabel(r'$|c_{'+str(operator)+r'}^p|_\mathrm{upper}$ [GeV$^{-2}$]')
                if not any(np.isfinite(r['limit_value']) and r['limit_value']>0 for r in run['limits'] if r['operator']==operator):
                    ax.set_ylim(1e-12,1);ax.text(.05,.5,'No finite sensitivity',transform=ax.transAxes)
            ax.set_title(style['title']+f' — O{operator}\n'+c['execution']['resolved_label'])
            ax.text(.02,.02,f"Workbench CL={c['analysis']['confidence_level']:.0%}; jχ={c['signal']['j_chi']}; cn/cp={c['signal']['fn_over_fp']}",transform=ax.transAxes,fontsize=7)
            ax.grid(alpha=.2);ax.legend(fontsize=6.5,ncol=2,loc='best')
            fig.text(.02,-.04,'Provisional Ge/response/activity assumptions; external curves are frozen snapshots, not guaranteed current in 2026.',fontsize=7)
            name=f'limit_O{operator}'
            for ext in ('pdf','png'):
                path=out/(name+'.'+ext);fig.savefig(path,dpi=style['dpi']);files.append(path)
            plt.close(fig)
        if 1 in c['signal']['operators']:
            fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
            for key,rows in limit_series(run,1).items():
                scenario,bg,exposure,method=key
                if scenario!=style['exposure_comparison_scenario'] or bg!='revised':continue
                x=np.array([r['mass_GeV'] for r in rows]);y=np.array([r['limit_value'] for r in rows])
                valid=np.isfinite(y)&(y>0)
                ax.plot(x[valid],y[valid],marker='o',label=f'{exposure/year:g} kg yr; {method}')
            ax.set(xscale='log',yscale='log',xlabel='WIMP mass [GeV]',ylabel=r'$\sigma_{p,\mathrm{upper}}$ [cm$^2$]',
                   title='Exposure comparison — '+style['exposure_comparison_scenario'],
                   xlim=style['mass_range_GeV'],ylim=style['sigma_range_cm2'])
            ax.grid(alpha=.2)
            if ax.lines:ax.legend(fontsize=8)
            else:ax.text(.05,.5,'Exposure-plot scenario is not selected',transform=ax.transAxes)
            for ext in ('pdf','png'):
                path=out/('exposure_comparison.'+ext);fig.savefig(path,dpi=style['dpi']);files.append(path)
            plt.close(fig)
    metadata=dict(numerical_run_id=run['manifest']['run_id'],numerical_manifest_sha256=file_hash(run['path']/'manifest.json'),
        plotting=style,comparison_provenance=compare,
        figure_sha256={p.name:file_hash(p) for p in files})
    (out/'plot_manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return files
