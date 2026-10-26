"""Content-keyed numerical storage; no pickle or executable cached objects."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import platform
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_identity():
    import wimp_legend
    if Path(wimp_legend.__file__).resolve().parent != (ROOT/'src/wimp_legend').resolve():
        raise ValueError('Workbench must use this checkout: install with python -m pip install -e .')
    # Fingerprint actual files, not only HEAD: dirty checkouts cannot reuse a
    # stale physics cache after an uncommitted implementation/data modification.
    paths=list((ROOT/'src/wimp_legend').rglob('*.py'))
    paths += [p for p in (ROOT/'src/wimp_legend').rglob('data/**/*') if p.is_file()]
    paths += [p for p in Path(__file__).parent.glob('*.py') if p.name!='plotting.py']
    paths += [ROOT/'tools/background_benchmark'/name for name in ('folding.py','legacy.py')]
    hashes={str(p.relative_to(ROOT)):file_hash(p) for p in sorted(set(paths))}
    def git(*args):
        try:return subprocess.check_output(['git',*args],cwd=ROOT,text=True).strip()
        except (OSError,subprocess.CalledProcessError):return 'unavailable'
    import scipy
    environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    return dict(git_commit=git('rev-parse','HEAD'),git_status=git('status','--porcelain'),
        physics_files=hashes,physics_fingerprint=digest(dict(files=hashes,environment=environment)),
        environment=environment)


class NumericalCache:
    def __init__(self,path,fingerprint):
        self.path=Path(path);self.path.mkdir(parents=True,exist_ok=True)
        self.fingerprint=fingerprint;self.hits=0;self.misses=0;self.used=[]

    def arrays(self,kind,inputs,compute):
        key=digest(dict(kind=kind,inputs=inputs,physics=self.fingerprint))
        file=self.path/(key+'.npz');meta=self.path/(key+'.json')
        if file.exists() and meta.exists():
            record=json.loads(meta.read_text())
            if record['key']!=key or record['sha256']!=file_hash(file):
                raise ValueError('numerical cache integrity failure: '+key)
            with np.load(file,allow_pickle=False) as data:result={k:data[k].copy() for k in data.files}
            self.hits+=1
        else:
            result={k:np.asarray(v) for k,v in compute().items()}
            if any(v.dtype.kind not in 'iufb' or np.any(~np.isfinite(v)) for v in result.values()):
                raise ValueError('cache arrays must be finite numeric data')
            temporary=file.with_suffix('.tmp.npz');np.savez_compressed(temporary,**result);temporary.replace(file)
            meta.write_text(json.dumps(dict(key=key,kind=kind,inputs=inputs,sha256=file_hash(file)),indent=2)+'\n')
            self.misses+=1
        self.used.append(key)
        return result
