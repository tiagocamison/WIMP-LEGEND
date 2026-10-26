"""Focused config, mixture, count/statistical and notebook workflow checks."""
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
pytest.importorskip('threadpoolctl')
pytest.importorskip('matplotlib')
from analyses.legend_poster import load_config,save_config,validate_config,resolve_profile,run_analysis,load_run
from analyses.legend_poster.analysis import observed_edges,roi_slice,count_limit,signal_model,fold_component,convergence
from analyses.legend_poster.models import target,isotope_weights,responses,build_background
from analyses.legend_poster.storage import NumericalCache,digest,file_hash
from analyses.legend_poster.folding import LocalGaussianBins
from analyses.legend_poster.plotting import comparison_curves,limit_series,plot_run
from wimp_legend.backgrounds import TruthBackground,MonoenergeticLine
from wimp_legend.energy import EnergyDomain
from wimp_legend.si import si_proton_coefficient_GeV_minus2,si_proton_cross_section_cm2
from wimp_legend.statistics import PoissonCLsCounting


@pytest.fixture(scope='module')
def config():return load_config()


def test_configuration_owned_round_trip(config,tmp_path):
    original=deepcopy(config);p=save_config(config,tmp_path/'config.json')
    assert load_config(p)==original
    owned=validate_config(config);owned['signal']['masses_GeV'].append(999)
    assert config==original
    # A failed save must not replace an existing valid configuration.
    broken=deepcopy(config);broken['signal']['j_chi']=.1
    before=p.read_bytes()
    with pytest.raises(ValueError):save_config(broken,p)
    assert p.read_bytes()==before


@pytest.mark.parametrize('path,value',[
 (('signal','convention_id'),'incompatible'),
 (('signal','operators'),[1,True]),
 (('signal','m_N_reference_GeV'),0),
 (('signal','halo','v0_km_s'),True),
 (('signal','masses_GeV'),[]),
 (('target','presets','natural','fraction_kind'),'unspecified'),
 (('target','presets','natural','fractions','76'),.5),
 (('target','isotopes','76','mass_GeV'),76*.9315),
 (('detector','presets','benchmark','efficiency'),1.01),
 (('analysis','confidence_level'),True),
 (('analysis','roi','signal_fraction'),True),
 (('analysis','methods'),['approximate_unlabelled']),
 (('numerical','gaussian_tail_sigma'),1),
 (('plotting','comparisons'),['AllLimits-2024']),
])
def test_invalid_configuration_rejected(config,path,value):
    c=deepcopy(config);node=c
    for key in path[:-1]:node=node[key]
    node[path[-1]]=value
    with pytest.raises(ValueError):validate_config(c)


def test_profiles_only_change_explicit_saved_choices(config):
    smoke=resolve_profile(config,'smoke')
    assert smoke['signal']['masses_GeV']==config['execution']['profiles']['smoke']['masses_GeV']
    assert smoke['target']==config['target'] and smoke['background']==config['background']
    assert smoke['execution']['resolved_label'].startswith('SMOKE')
    with pytest.raises(ValueError):resolve_profile(config,'unknown')


def test_isotope_weights_bulk_and_nuclear_semantics(config):
    for name in config['target']['presets']:
        w=isotope_weights(config,name);comp=target(config,name)
        assert sum(v['number_fraction'] for v in w.values())==pytest.approx(1,abs=1e-14)
        assert sum(v['bulk_mass_fraction'] for v in w.values())==pytest.approx(1,abs=1e-14)
        for a,row in w.items():
            assert row['bulk_mass_fraction']*row['pure_targets_per_kg']==pytest.approx(row['targets_per_kg'],rel=2e-15)
        mass_version=comp.to_mass_fractions()
        c=deepcopy(config);c['target']['presets'][name]['fraction_kind']='mass'
        c['target']['presets'][name]['fractions']={str(i.A):f for i,f in mass_version}
        back=isotope_weights(c,name)
        for a in w:assert back[a]['targets_per_kg']==pytest.approx(w[a]['targets_per_kg'],rel=2e-15)
    assert config['target']['presets']['enriched_provisional']['status']=='hypothetical'
    assert len([w for w in isotope_weights(config,'natural').values() if w['number_fraction']>0])==5
    assert len([w for w in isotope_weights(config,'enriched_provisional').values() if w['number_fraction']>0])==2


def test_nuclear_masses_keep_electron_correction_and_source(config):
    for isotope,_ in target(config,'natural'):
        row=config['target']['isotopes'][str(isotope.A)];conversion=config['target']['mass_conversion']
        assert isotope.mass_GeV<row['atomic_mass_u']*conversion['u_GeV']
        assert isotope.mass_GeV!=pytest.approx(isotope.A*.9315,rel=1e-5)
        assert 'NIST' in isotope.sources['mass_GeV'].reference


def test_bin_edges_fine_and_threshold_not_duplicated(config):
    c=deepcopy(config);c['detector']['presets']['benchmark']['threshold_keV_obs']=.125
    e=observed_edges(c)
    assert len(e)==len(np.unique(e)) and np.diff(e).min()>=.005-1e-12
    assert np.count_nonzero(e==.125)==1
    assert np.diff(e)[(e[:-1]>.2)&(e[1:]<1.9)]==pytest.approx(.01)
    c['numerical']['bin_segments'][0][2]=.005
    fine=observed_edges(c);assert len(fine)>len(e)
    s=np.ones(len(e)-1);roi,_=roi_slice(s,e,.12500000000000002,.95)
    assert e[roi.start]==.125


def test_background_mixture_and_domain_routing(config):
    m=build_background(config,'natural_AAr')
    assert len(m.components['solar_CEvNS'].components)==5
    nt=sum(v['targets_per_kg'] for v in isotope_weights(config,'natural').values())
    assert m.components['solar_ES'].spectrum.targets_per_kg==pytest.approx(nt,rel=2e-15)
    assert m.components['H3'].energy_domain is EnergyDomain.ER
    assert all(c.energy_domain is EnergyDomain.NR for c in m.components['solar_CEvNS'].components)
    c=deepcopy(config);c['background']['scenarios']['legend_poster']['components']['flat_observed']={'response':'benchmark'}
    with pytest.raises(ValueError):validate_config(c)
    c=deepcopy(config);c['detector']['presets']['sharper']=deepcopy(c['detector']['presets']['benchmark'])
    c['detector']['presets']['sharper']['noise_sigma_keV']/=2
    c['background']['scenarios']['legend_poster']['components']['Fe55']={'response':'sharper'}
    m=build_background(c,'legend_poster')
    assert m.components['Fe55'].response.resolution.sigma_at(0)<m.components['H3'].response.resolution.sigma_at(0)


def test_pruned_probabilities_delegate_core_and_preserve_lines(config):
    er,nr=responses(config)
    bins=np.linspace(-.1,4,500);energy=np.array([0,.1,.8,1.5,3])
    for response in (er,nr):
        proxy=LocalGaussianBins(response,12.)
        np.testing.assert_allclose(proxy.bin_probabilities(energy,bins),
            response.bin_probabilities(energy,bins),rtol=1e-12,atol=1e-30)
        np.testing.assert_array_equal(proxy.expected_line_counts(.13,2.,bins,3.),
                                      response.expected_line_counts(.13,2.,bins,3.))


def test_integrated_bins_recomputed_not_density_interpolation(config):
    er,_=responses(config);bins=np.array([.03,.125,.3,.8,1.5])
    source=TruthBackground(EnergyDomain.ER,er,lambda e:np.exp(-np.asarray(e)),(0.,2.),
        lines=(MonoenergeticLine(.16,.7,EnergyDomain.ER),))
    n=config['numerical'];coarse=fold_component(source,bins,n)
    fine_bins=np.sort(np.r_[bins,(bins[:-1]+bins[1:])/2]);fine=fold_component(source,fine_bins,n)
    np.testing.assert_allclose(coarse,fine.reshape(-1,2).sum(axis=1),rtol=5e-10,atol=1e-13)
    # The observed-flat component bypasses both detector kernels/efficiencies.
    model=build_background(config,'legend_poster')
    flat=model.components['flat_observed'];expected=flat.rate_per_kg_day_keV_obs*np.diff(bins)
    np.testing.assert_array_equal(fold_component(flat,bins,n),expected)


@pytest.mark.parametrize('operator',[1,8,11])
def test_si_and_single_operator_conventions(config,operator):
    isotope=next(i for i,_ in target(config,'natural') if i.A==76)
    model,_=signal_model(resolve_profile(config,'smoke'),isotope,6.,operator)
    coeff=model.coefficients if operator==1 else model.keywords['coefficients']
    assert coeff.convention_id==config['signal']['convention_id']
    assert all(coeff.get(op,t)==0 for op in (1,8,11) if op!=operator for t in (0,1))
    cp,cn=coeff.to_proton_neutron()[operator]
    assert cn==pytest.approx(cp*config['signal']['fn_over_fp'])
    if operator==1:
        assert si_proton_cross_section_cm2(cp,m_chi_GeV=6,m_proton_GeV=config['signal']['m_proton_GeV'])==pytest.approx(config['signal']['sigma_p_ref_cm2'],rel=3e-15)
        assert model.cross_section_limit_cm2(3)==3*config['signal']['sigma_p_ref_cm2']
    else:assert cp==config['signal']['coefficient_ref_GeV_minus2']


def test_count_roi_exact_cls_and_exposure(config):
    bins=np.array([0,.1,.2,.3,.4]);signal=np.array([9.,4.,3.,1.]);background=np.array([99.,.1,.2,.3])
    one,record=count_limit(signal,background,bins,2.,.1,config,'poisson_cls_counting')
    direct=PoissonCLsCounting().expected_upper_limit([16.],[1.2],confidence_level=.9)
    assert one.limit.mu_upper==pytest.approx(direct.mu_upper,rel=2e-14)
    assert record['roi_lower_keV']==.1 and record['retained_signal_fraction']==1
    two,_=count_limit(signal,background,bins,4.,.1,config,'poisson_cls_counting')
    np.testing.assert_array_equal(two.signal_counts,2*np.asarray(one.signal_counts))
    np.testing.assert_array_equal(two.background_counts,2*np.asarray(one.background_counts))
    binned,record=count_limit(signal,background,bins,2.,.1,config,'binned_asimov')
    assert binned.limit.kind=='expected_asimov' and 'asymptotic' in binned.limit.approximation


def test_cache_integrity_and_owned_arrays(tmp_path):
    cache=NumericalCache(tmp_path,'source-A');calls=[]
    def compute():calls.append(1);return {'x':np.array([1.,2.])}
    a=cache.arrays('test',{'parameter':1},compute);a['x'][0]=999
    b=cache.arrays('test',{'parameter':1},compute)
    np.testing.assert_array_equal(b['x'],[1,2]);assert len(calls)==1
    cache.arrays('test',{'parameter':2},compute);assert len(calls)==2
    file=next(tmp_path.glob('*.npz'));file.write_bytes(file.read_bytes()+b'broken')
    with pytest.raises(ValueError,match='integrity'):
        for parameter in (1,2):cache.arrays('test',{'parameter':parameter},compute)


def test_comparison_snapshots_verified_individually(config):
    manifest,curves=comparison_curves(config['plotting']['comparisons'])
    assert len(curves)==7 and len(manifest['commit'])==40
    assert {r['kind'] for r,_ in curves}=={'projection','observed_limit'}
    for record,xy in curves:
        raw=(ROOT/'analyses/legend_poster/data/si_limits'/record['file']).read_bytes()
        assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==record['upstream_git_blob_sha']
        assert record['confidence_level']==.9 and '2026' in record['currency']
        assert xy.shape[1]==2


@pytest.fixture(scope='module')
def completed_run(config,tmp_path_factory):
    root=tmp_path_factory.mktemp('poster-calculation')
    # Use the production cache if supplied: cross-process reuse is part of the
    # notebook workflow, but its source/config checks remain active.
    cache=Path(os.environ.get('WORKBENCH_TEST_CACHE',root/'cache')).resolve()
    r=run_analysis(config,profile='smoke',cache_dir=cache,output_dir=root/'runs')
    return r,cache,root


def test_real_run_mixtures_limits_and_manifests(completed_run,config):
    run,_,_=completed_run
    assert all(c['passed'] for c in run['convergence'])
    manifest=run['manifest'];assert 'physics_files' in manifest['source']
    assert manifest['resolved_configuration']['execution']['resolved_profile']=='smoke'
    assert set(manifest['nuclear_masses_GeV'])=={'70','72','73','74','76'}
    for scenario,record in config['background']['scenarios'].items():
        for mass in (6.,20.):
            for op in (1,8,11):
                mixed=run['signals'][f'{scenario}__O{op}__m{mass:g}']
                expected=np.zeros_like(mixed)
                for a,row in isotope_weights(config,record['target']).items():
                    if row['bulk_mass_fraction']:
                        expected+=row['bulk_mass_fraction']*run['signals'][f'pure_Ge{a}__O{op}__m{mass:g}__benchmark']
                np.testing.assert_array_equal(mixed,expected)
    si=[r for r in run['limits'] if r['operator']==1 and r['background']=='revised']
    assert len(set(round(r['limit_value']/config['signal']['sigma_p_ref_cm2'],10) for r in si))>2
    for row in run['limits']:
        expected=row['mu_upper']*config['signal']['sigma_p_ref_cm2'] if row['operator']==1 else np.sqrt(row['mu_upper'])*config['signal']['coefficient_ref_GeV_minus2']
        assert row['limit_value']==pytest.approx(expected,rel=1e-14)
        assert row['limit_unit']==('cm^2' if row['operator']==1 else 'GeV^-2')
    reloaded=load_run(run['path']);assert reloaded['limits']==run['limits']


@pytest.mark.parametrize('operator',[1,8,11])
def test_physical_reference_rescaling(config,operator):
    c=resolve_profile(config,'smoke');c['signal']['fn_over_fp']=.4
    isotope=next(i for i,_ in target(c,'natural') if i.A==76)
    reference,_=signal_model(c,isotope,6.,operator)
    changed=deepcopy(c)
    if operator==1:
        changed['signal']['sigma_p_ref_cm2']*=3;factor=3
    else:
        changed['signal']['coefficient_ref_GeV_minus2']*=2;factor=4
    scaled,_=signal_model(changed,isotope,6.,operator)
    energy=np.array([.1,.3,.8])
    original=np.asarray(reference(energy))
    assert np.all(original>0)
    np.testing.assert_allclose(scaled(energy),factor*original,rtol=3e-13,atol=0)


def test_physical_fine_binning_and_mass_dependence(completed_run,config):
    from analyses.legend_poster.analysis import pure_signal_counts
    from analyses.legend_poster.storage import source_identity
    run,cache_dir,_=completed_run
    c=resolve_profile(config,'smoke');c['numerical']['bin_segments'][0][2]=.005
    fine_bins=observed_edges(c);coarse_bins=run['signals']['observed_edges_keV']
    isotope=next(i for i,_ in target(c,'natural') if i.A==76)
    cache=NumericalCache(cache_dir,source_identity()['physics_fingerprint'])
    fine,check=pure_signal_counts(c,isotope,6.,1,'benchmark',fine_bins,cache)
    indices=np.searchsorted(fine_bins,coarse_bins)
    np.testing.assert_array_equal(fine_bins[indices],coarse_bins)
    rebinned=np.array([fine[a:b].sum() for a,b in zip(indices[:-1],indices[1:])])
    np.testing.assert_allclose(rebinned,run['signals']['pure_Ge76__O1__m6__benchmark'],rtol=3e-5,atol=1e-20)
    assert check['passed']
    for op in (1,8,11):
        rows=[r for r in run['limits'] if r['operator']==op and r['scenario']=='legend_poster'
              and r['exposure_kg_day']==config['analysis']['exposures_kg_day'][-1]]
        assert {r['mass_GeV'] for r in rows}=={6,20}
        # Sensitivity may have a broad plateau: mass dependence does not imply
        # a minimum one-percent change between arbitrary chosen masses.
        assert not np.isclose(rows[0]['limit_value']/rows[1]['limit_value'],1,rtol=1e-8,atol=0)


def test_exposure_and_styling_reuse_physics_cache(completed_run,config,monkeypatch):
    import analyses.legend_poster.analysis as module
    from wimp_legend.si import SIO1Signal
    run,cache,root=completed_run
    def forbidden(*args,**kwargs):raise AssertionError('expensive rate evaluation repeated')
    monkeypatch.setattr(SIO1Signal,'__call__',forbidden)
    monkeypatch.setattr(module,'differential_rate_per_kg_day_keV',forbidden)
    c=deepcopy(config);c['plotting']['title']='Changed styling only'
    c['analysis']['exposures_kg_day']=[2*x for x in c['analysis']['exposures_kg_day']]
    warm=run_analysis(c,profile='smoke',cache_dir=cache,output_dir=root/'other-runs')
    assert warm['manifest']['cache']['misses']==0
    for key,values in run['signals'].items():np.testing.assert_array_equal(values,warm['signals'][key])
    for r in run['limits']:
        w=next(v for v in warm['limits'] if all(v[k]==r[k] for k in ('scenario','background','mass_GeV','operator','method')) and v['exposure_kg_day']==2*r['exposure_kg_day'])
        assert w['reference_signal_events']==pytest.approx(2*r['reference_signal_events'],rel=1e-14)
        assert w['background_events']==pytest.approx(2*r['background_events'],rel=1e-14)


def test_notebook_plots_use_saved_outputs(completed_run,config,tmp_path):
    run,_,_=completed_run
    files=plot_run(run,output_dir=tmp_path)
    assert any(p.suffix=='.pdf' for p in files)
    record=json.loads((tmp_path/'plot_manifest.json').read_text())
    assert record['numerical_manifest_sha256']==file_hash(run['path']/'manifest.json')
    assert all(file_hash(tmp_path/name)==value for name,value in record['figure_sha256'].items())
    for key,rows in limit_series(run,8).items():
        assert all(r['operator']==8 and r['limit_unit']=='GeV^-2' for r in rows)
    notebook=json.loads((ROOT/'notebooks/poster_workbench.ipynb').read_text())
    code='\n'.join(''.join(c.get('source',[])) for c in notebook['cells'] if c['cell_type']=='code')
    assert 'run_analysis(' in code and 'plot_run(' in code and 'pilot_sensitivity_results.csv' not in code


def test_historical_assets_untouched(config):
    p=ROOT/config['provenance']['pilot_notebook']
    assert file_hash(p)==config['provenance']['pilot_sha256']
    base='f62e060ddedac75596a17e77af4053184dc07192'
    original=subprocess.check_output(['git','show',base+':notebooks/poster_workbench.ipynb'],cwd=ROOT)
    assert p.read_bytes()==original
    for file in (ROOT/'notebooks/poster_workbench_assets').iterdir():
        if file.is_file():assert file.read_bytes()==subprocess.check_output(['git','show',base+':'+str(file.relative_to(ROOT))],cwd=ROOT)


def test_notebooks_execute_clean_process(completed_run,config,tmp_path,monkeypatch):
    nbformat=pytest.importorskip('nbformat');nbclient=pytest.importorskip('nbclient')
    run,cache,root=completed_run
    p=save_config(config,tmp_path/'config.json')
    before=p.read_bytes()
    monkeypatch.setenv('WIMP_POSTER_CONFIG',str(p));monkeypatch.setenv('WIMP_POSTER_PROFILE','smoke')
    monkeypatch.setenv('WIMP_POSTER_CACHE',str(cache));monkeypatch.setenv('WIMP_POSTER_RUNS',str(tmp_path/'notebook-runs'))
    # Managed CI forbids Jupyter's NIC discovery and both TCP/IPC ZMQ sockets.
    # Execute every actual cell using a fresh IPython subprocess per notebook,
    # retaining rich outputs and failing on any cell error. No physics is mocked.
    runner = r"""
import sys
import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output
notebook=nbformat.read(sys.argv[1],as_version=4)
shell=InteractiveShell.instance()
for index,cell in enumerate(notebook.cells):
    if cell.cell_type!='code':continue
    with capture_output() as captured:
        result=shell.run_cell(cell.source,store_history=True)
    if result.error_before_exec or result.error_in_exec:
        raise RuntimeError('Notebook cell %d failed: %s\n%s' % (index,result.error_before_exec or result.error_in_exec,captured.stderr))
    cell.execution_count=shell.execution_count-1
    cell.outputs=[]
    for name,value in [('stdout',captured.stdout),('stderr',captured.stderr)]:
        if value:cell.outputs.append(nbformat.v4.new_output('stream',name=name,text=value))
    for output in captured.outputs:
        cell.outputs.append(nbformat.v4.new_output('display_data',data=output.data,metadata=output.metadata))
nbformat.validate(notebook)
nbformat.write(notebook,sys.argv[2])
"""
    for name in ('poster_setup.ipynb','poster_workbench.ipynb'):
        result=subprocess.run([sys.executable,'-c',runner,str(ROOT/'notebooks'/name),str(tmp_path/name)],
            cwd=ROOT/'notebooks',capture_output=True,text=True,timeout=600)
        assert result.returncode==0,result.stdout+'\n'+result.stderr
        notebook=nbformat.read(tmp_path/name,as_version=4)
        assert all(output.get('output_type')!='error' for cell in notebook.cells for output in cell.get('outputs',[]))
        assert all(cell.execution_count is not None for cell in notebook.cells if cell.cell_type=='code')
    assert p.read_bytes()==before
