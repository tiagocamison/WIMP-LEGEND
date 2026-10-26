"""Read-only execution of selected pure historical notebook definitions.

Only the named function definitions are loaded, never imports, cells, plotting,
statistics or stored outputs. A local NumPy proxy maps removed np.trapz to its
identical trapezoidal API; the notebook and NumPy module are not modified.
"""
import ast,hashlib,json
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn

NAMES=('get_detector_resolution','gaussian_peak','generate_EC_peaks',
       'beta_spectrum_tritium','beta_spectrum_Ar39','compute_neutrino_backgrounds','compute_background')


def load_legacy(notebook):
    path=Path(notebook);book=json.loads(path.read_text());functions={}
    for cell in book['cells']:
        if cell['cell_type']!='code':continue
        source=''.join(cell.get('source',[]))
        if not any('def '+name+'(' in source for name in NAMES):continue
        # Historical plotting blocks contain invalid Python escape sequences;
        # discard that unexecuted block before parsing, without editing source.
        source=source.split('if __name__')[0]
        try:tree=ast.parse(source)
        except SyntaxError:continue
        for node in tree.body:
            if isinstance(node,ast.FunctionDef) and node.name in NAMES:functions[node.name]=node
    if set(functions)!=set(NAMES):raise ValueError('notebook lacks expected historical background definitions')
    class NumpyCompatibility:
        def __getattr__(self,name):
            if name=='trapz':return getattr(np,'trapezoid',None) or np.trapz
            return getattr(np,name)
    scope={'np':NumpyCompatibility(),'spherical_jn':spherical_jn}
    tree=ast.Module(body=[functions[n] for n in NAMES],type_ignores=[])
    exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),scope)
    return scope['compute_background'],hashlib.sha256(path.read_bytes()).hexdigest()
