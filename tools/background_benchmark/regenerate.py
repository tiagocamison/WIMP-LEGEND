"""Optional Linux/x86-64 input regeneration; no external tool at production load.

Download the hash-pinned official BetaShape2.4 binary package and Bahcall source
spectra into --work. Run the supplied benchmark, then the archived minimal ENSDF
input. Export portable data with freeze_inputs.py. No online auto-update occurs.
"""
import argparse,hashlib,json,os,platform,re,shutil,subprocess,sys,tarfile
from pathlib import Path
from urllib.request import urlopen

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'src/wimp_legend/backgrounds/data/benchmark_v1'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if platform.system()!='Linux' or platform.machine() not in ('x86_64','AMD64'):
        raise RuntimeError('reference regeneration uses the official Linux x86-64 binary; production data loading is portable')
    work=args.work.resolve();work.mkdir(parents=True,exist_ok=True)
    provenance=json.loads((DATA/'ar39_provenance.json').read_text())
    solar=json.loads((DATA/'solar.json').read_text())
    def download(url,path,expected):
        if not path.exists():
            with urlopen(url,timeout=120) as response:path.write_bytes(response.read())
        if digest(path)!=expected:raise ValueError(f'source checksum mismatch: {path.name}')
    revision=provenance['repository_commit']
    url=f'https://raw.githubusercontent.com/IAEA-NSDDNetwork/BetaShape/{revision}/packages/BetaShape_v2.4_Ubuntu_20.04.2_LTS_static.tar.bz2'
    archive=work/'BetaShape_v2.4.tar.bz2';download(url,archive,provenance['archive_sha256'])
    (work/'betashape_commit.txt').write_text(revision+'\n')
    for component,name in zip(solar['components'],('pp.txt','b8.txt')):
        download(component['url'],work/name,component['raw_sha256'])
    binary=work/'betashape/BetaShape_v2.4_Ubuntu_20.04.2_LTS_static'
    with tarfile.open(archive) as package:
        package.extractall(work/'betashape',filter='data')
        for member in package.getmembers():
            if member.isfile():
                original=package.extractfile(member).read()
                if hashlib.sha256(original).hexdigest()!=digest(work/'betashape'/member.name):
                    raise ValueError(f'extracted official package differs: {member.name}')
    env=os.environ.copy();env['BSINSTALL']=str(binary);env['PATH']=str(binary)+os.pathsep+env['PATH']
    with (work/'official_benchmark.log').open('w') as log:
        subprocess.run([str(binary/'mybench')],cwd=binary,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    log=(work/'official_benchmark.log').read_text(errors='replace')
    if any(word in log for word in ('Permission denied','Segmentation fault','not found')) or re.search(r'^\d+(?:,\d+)?[acd]\d+',log,re.M):
        raise RuntimeError('official benchmark has errors or reference differences; inspect official_benchmark.log')
    run=work/'Ar39_run';run.mkdir(exist_ok=True)
    shutil.copyfile(DATA/'Ar39_minimal_ENSDF.txt',work/'Ar39_minimal_ENSDF.txt')
    shutil.copyfile(DATA/'Ar39_minimal_ENSDF.txt',run/'Ar39.txt')
    with (work/'Ar39_betashape.log').open('w') as log:
        subprocess.run([str(binary/'betashape'),'Ar39.txt','-csv','myEstep=0.1'],cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    subprocess.run([sys.executable,str(Path(__file__).with_name('freeze_inputs.py')),
                    '--raw-dir',str(work),'--output',str(args.output.resolve())],check=True)


if __name__=='__main__':main()
