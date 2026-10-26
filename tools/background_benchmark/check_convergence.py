"""Optional reproducible numerical refinement check for the plotting fold."""
from pathlib import Path
import json,sys
import numpy as np
import argparse
repo=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description='Check benchmark GL folding refinement; not a physical validation')
parser.add_argument('--config',type=Path,default=repo/'examples/bg_benchmark_config.json')
parser.add_argument('--scenario',default='legend_poster')
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
sys.path[:0]=[str(repo/'src'),str(repo/'tools/background_benchmark')]
from folding import quadrature_counts
from wimp_legend.backgrounds import BenchmarkConfig,build_benchmark
from wimp_legend.detector import ElectronRecoilResponse,DetectorResponse,GaussianResolution,ConstantEfficiency,LindhardQuenching
from wimp_legend.energy import EnergyDomain
cfg=json.loads(args.config.read_text());p=cfg['plotting']
sigma=lambda e:np.sqrt(p['noise_sigma_keV']**2+p['fano_factor']*p['pair_energy_keV']*np.asarray(e))
er=ElectronRecoilResponse(GaussianResolution(sigma),observed_efficiency=ConstantEfficiency(p['constant_detection_efficiency']))
nr=DetectorResponse(LindhardQuenching(32,p['lindhard_k']),GaussianResolution(sigma),observed_efficiency=ConstantEfficiency(p['constant_detection_efficiency']))
model=build_benchmark(BenchmarkConfig(**cfg['scenarios'][args.scenario]),er_response=er,nr_response=nr)
bins=np.unique(np.r_[np.arange(p['min_energy_keV'],p['max_energy_keV'],p['bin_width_keV']),p['threshold_keV'],p['max_energy_keV']])
result={};all_reference={}
for key,component in model.components.items():
 def bound(c):
  e=c.spectrum_energy_range_keV[1]
  return float(sigma(e if c.energy_domain is EnergyDomain.ER else nr.quenching.mean_energy_keV(e)))
 coarse=quadrature_counts(component,bins,p['exposure_kg_day'],step_keV=p['truth_step_keV'],order=p['quadrature_order'],tail_sigma=p['gaussian_tail_sigma'],sigma_bound=bound)
 fine=quadrature_counts(component,bins,p['exposure_kg_day'],step_keV=p['truth_step_keV']/2,order=10,tail_sigma=p['gaussian_tail_sigma'],sigma_bound=bound)
 mask=bins[:-1]>=p['threshold_keV'];den=max(float(fine[mask].sum()),1e-30)
 active=fine>fine.max()*1e-6
 result[key]=dict(coarse_events_in_window=float(coarse[mask].sum()),refined_events_in_window=float(fine[mask].sum()),relative_integrated_difference=float(abs(coarse[mask].sum()-fine[mask].sum())/den),max_absolute_bin_difference=float(np.max(abs(coarse-fine))),max_relative_bin_difference_above_one_millionth_peak=float(np.max(abs((coarse-fine)[active]/fine[active]))) if np.any(active) else 0.)
 all_reference[key]=fine.tolist()
 print(key,result[key],flush=True)
 assert result[key]['relative_integrated_difference']<1e-3,key
out=args.output
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(dict(coarse=dict(step_keV=p['truth_step_keV'],order=p['quadrature_order']),reference=dict(step_keV=p['truth_step_keV']/2,order=10),results=result,reference_bin_counts=all_reference),indent=2)+'\n')
