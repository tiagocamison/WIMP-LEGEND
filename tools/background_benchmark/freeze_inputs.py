"""Freeze portable source spectra from the documented raw inputs.

--raw-dir contains downloaded pp.txt / b8.txt and Ar39_run/Ar39/*.bs.
No spectrum is normalized on an analysis grid. See docs/BG_BENCHMARK.md.
The Ar39 input is a minimal transcription of ENSDF, not a raw NNDC download.
"""
import argparse,hashlib,json,re,shutil
from pathlib import Path
import numpy as np


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(p):
    data=[]
    for line in p.read_text().splitlines():
        tokens=line.split()
        if len(tokens)>=2:
            try:data.append([float(x) for x in tokens])
            except ValueError:pass
    return data

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[2]/'src/wimp_legend/backgrounds/data/benchmark_v1')
    a=parser.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=True)
    nuclear=Path(__file__).resolve().parents[2]/'src/wimp_legend/backgrounds/data/benchmark_v1/nuclear.json'
    if nuclear.resolve()!=(out/'nuclear.json').resolve():shutil.copyfile(nuclear,out/'nuclear.json')
    pp=np.array(rows(a.raw_dir/'pp.txt')).reshape(-1,2);pp=pp[np.argsort(pp[:,0])];pp=np.vstack(([0,0],pp))
    b8=np.array(rows(a.raw_dir/'b8.txt'))[:,:2]
    solar=dict(version='B16-GS98_pp_B8_benchmark_v1',flux_source='Vinyoles et al., ApJ 835,202(2017), DOI:10.3847/1538-4357/835/2/202, Table6',
        coverage='pp and B8 only; Be7/pep/CNO/hep omitted, not a full solar background',components=[])
    refs={'pp':('pp.txt',5.98e10,.006,'Bahcall, PRC56,3391(1997), TableIV','https://www.sns.ias.edu/~jnb/SNdata/Export/PPenergyspectrum/ppenergytab'),
          'B8':('b8.txt',5.46e6,.12,'Bahcall et al., PRC54,411(1996), central spectrum','https://www.sns.ias.edu/~jnb/SNdata/Export/B8spectrum/b8spectrum.txt')}
    for name,data in [('pp',pp),('B8',b8)]:
        file,flux,error,ref,url=refs[name]
        norm=float(np.sum(np.diff(data[:,0])*(data[:-1,1]+data[1:,1])/2))
        solar['components'].append(dict(name=name,total_flux_cm2_s=flux,relative_model_uncertainty=error,shape_reference=ref,url=url,
            raw_sha256=sha(a.raw_dir/file),source_table_integral=norm,normalization='piecewise-linear full-support integral, explicitly divided by source_table_integral',
            energies_MeV=data[:,0].tolist(),density_per_MeV=(data[:,1]/norm).tolist()))
        shutil.copyfile(a.raw_dir/file,out/('source_'+file))
    (out/'solar.json').write_text(json.dumps(solar,indent=2)+'\n')
    source=a.raw_dir/'Ar39_run/Ar39/beta-_Ar39_trans0_myEstep.bs';data=np.array(rows(source))[:,:2]
    if data[0,0]!=0 or data[-1,0]!=565 or data[-1,1]!=0:raise ValueError('Ar39 support is incomplete')
    np.savetxt(out/'ar39.csv',data,delimiter=',',header='T_keV,unnormalized_density_per_keV',comments='',fmt='%.12g')
    shutil.copyfile(a.raw_dir/'Ar39_run/Ar39/beta-_Ar39_trans0.bs',out/'source_ar39_default.bs')
    shutil.copyfile(a.raw_dir/'Ar39_run/Ar39.rpt',out/'source_ar39.rpt')
    shutil.copyfile(source,out/'source_ar39.bs');shutil.copyfile(a.raw_dir/'Ar39_minimal_ENSDF.txt',out/'Ar39_minimal_ENSDF.txt')
    ar=dict(tool='IAEA-NSDDNetwork/BetaShape 2.4 (06/2024)',tool_reference='X. Mougeot, Appl.Radiat.Isot.201(2023)111018, DOI:10.1016/j.apradiso.2023.111018',
        repository_commit=(a.raw_dir/'betashape_commit.txt').read_text().strip(),archive_sha256=sha(a.raw_dir/'BetaShape_v2.4.tar.bz2'),
        download='https://github.com/IAEA-NSDDNetwork/BetaShape',input_kind='minimal ENSDF transcription, not original downloaded evaluator file',
        input_sha256=sha(a.raw_dir/'Ar39_minimal_ENSDF.txt'),
        nuclear_reference='ENSDF 39K, Jun Chen, NDS149,1(2018); dataset 39Ar beta- decay (268Y), Q=565(5)keV, branch1, 7/2- -> 3/2+, 1U',
        source_url='https://www.nndc.bnl.gov/nudat3/getdecaydataset.jsp?dsid=39ar+bM+decay+%28268+y%29&nucleus=39K',
        command='betashape Ar39.txt -csv myEstep=0.1',corrections='default numerical screening, exchange, overlap and radiative; 1st-forbidden-unique shape; no experimental shape factor',
        original_output_sha256=sha(source),interpolation='BetaShape requested 0.1keV output; original default calculation grid1keV; runtime linear interpolation, full-support renormalization',
        uncertainty='endpoint 5keV source uncertainty; BetaShape point errors retained in source_ar39.bs, not propagated in benchmark',
        mean_energy_reference_keV=218.,normalization='runtime explicit full-support trapezoidal normalization; not analysis-window normalization')
    (out/'ar39_provenance.json').write_text(json.dumps(ar,indent=2)+'\n')
    manifest=dict(version='BG-BENCHMARK-001-v1',files={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='manifest.json'})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

if __name__=='__main__':main()
