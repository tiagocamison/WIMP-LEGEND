"""Poster-only SI campaign postprocessing. Does not recompute particle or detector physics.

Drop-in helper for WIMP-LEGEND/analyses/legend_poster on poster-speedup after merge.
Uses verified unit-exposure folded spectra and the existing exact Poisson CLs API.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import csv
import hashlib
import json
import math

import numpy as np

from ..analysis import count_limit, load_run

POSTER_SPEC_VERSION = 1


def _sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def _write_json(path, payload):
    path = Path(path)
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def load_spec(path):
    spec = json.loads(Path(path).read_text(encoding='utf-8'))
    if spec.get('schema_version') != POSTER_SPEC_VERSION:
        raise ValueError('unsupported poster campaign spec version')
    if spec.get('background_modes') != {'none': [], 'neutrinos': ['solar_CEvNS', 'solar_ES'], 'full': ['Total']}:
        raise ValueError('Unexpected background component selection. Review the benchmark model first.')
    if spec.get('thresholds_keVee') != [0.1, 0.15, 1.0]:
        raise ValueError('Expected precisely 0.1, 0.15 and 1.0 keVee thresholds')
    if spec.get('exposures_tonne_year') != [0.1, 1.0, 10.0]:
        raise ValueError('Expected precisely 0.1, 1 and 10 tonne-years')
    if spec.get('statistics_method') != 'poisson_cls_counting':
        raise ValueError('Unsupported statistical construction')
    if spec.get('input_scenario') != 'legend_poster':
        raise ValueError('Unexpected benchmark scenario')
    return spec


def _validated_inputs(run, spec):
    """Reject configurations for which the cheap threshold/background scan is invalid."""
    config = run['manifest']['resolved_configuration']
    if config['signal']['operators'] != [1] or config['analysis']['scenarios'] != [spec['input_scenario']]:
        raise ValueError('Source run must contain ONLY O1 and legend_poster (dedicated campaign).')
    if config['analysis']['methods'] != [spec['statistics_method']]:
        raise ValueError('Unexpected statistical method in source run')
    if config['signal'].get('si_rate_engine') != 'factorized_shm':
        raise ValueError('Fast SI backend was not enabled in source run')
    if config['numerical'].get('folding_backend') != 'vectorized':
        raise ValueError('Vectorized folding backend was not enabled in source run')
    scenario = spec['input_scenario']
    record = config['background']['scenarios'][scenario]
    detector = config['detector']['presets'][record.get('detector', config['detector']['default'])]
    if not np.isclose(detector['threshold_keV_obs'], min(spec['thresholds_keVee']), atol=1e-12, rtol=0):
        raise ValueError('Source detector threshold must equal lowest campaign threshold')
    # This workbench detector uses a constant observed efficiency; the threshold
    # here is only the analysis ROI cut, not an efficiency-turn-on parameter.
    if config['background']['historical']['enabled']:
        raise ValueError('Historical comparison must be disabled in this dedicated campaign')
    s_bins = np.asarray(run['signals']['observed_edges_keV'], dtype=float)
    b_bins = np.asarray(run['backgrounds']['observed_edges_keV'], dtype=float)
    if not np.array_equal(s_bins, b_bins):
        raise ValueError('Signal and background energy bins differ')
    if not np.all(np.diff(s_bins) > 0):
        raise ValueError('Invalid energy bins')
    for threshold in spec['thresholds_keVee']:
        if not np.any(np.isclose(s_bins, threshold, atol=1e-10, rtol=0)):
            raise ValueError(f'{threshold:g} keVee is not an observed-energy bin edge: refold onto common bins')
    if not math.isclose(config['analysis']['confidence_level'], .90, rel_tol=0, abs_tol=1e-12):
        raise ValueError('Unexpected CL')
    return config, s_bins


def build_limits(run_directory, spec_path, campaign_dir):
    """Compute all 27 SI combinations using existing rate folding and counting CLs.

    The only numerical work performed here is statistical CLs/ROI evaluation.
    Outputs contain 27*N_mass rows and explicit provenance.
    """
    run = load_run(run_directory)  # verifies numerical-run checksums
    spec = load_spec(spec_path)
    config, edges = _validated_inputs(run, spec)
    scenario = spec['input_scenario']
    prefix = scenario + '__'
    components = {}
    for name in ('solar_CEvNS', 'solar_ES', 'Total'):
        key = prefix + name
        if key not in run['backgrounds']:
            raise KeyError(f'Background component missing: {key}')
        components[name] = np.asarray(run['backgrounds'][key], dtype=float)
    backgrounds = {
        'none': np.zeros(len(edges) - 1, dtype=float),
        'neutrinos': components['solar_CEvNS'] + components['solar_ES'],
        'full': components['Total'],
    }
    for name, values in backgrounds.items():
        if values.shape != (len(edges) - 1,) or np.any(~np.isfinite(values)) or np.any(values < -1e-14):
            raise ValueError(f'Invalid {name} background counts')
    if np.any(backgrounds['neutrinos'] > backgrounds['full'] + 1e-10 * max(1, backgrounds['full'].max())):
        raise ValueError('Neutrino background exceeds full background; investigate components')
    mass_grid = list(map(float, config['signal']['masses_GeV']))
    sigma_ref = float(config['signal']['sigma_p_ref_cm2'])
    days = float(config['analysis']['days_per_year'])
    rows = []
    for mass in mass_grid:
        key = f'{scenario}__O1__m{mass:g}'
        if key not in run['signals']:
            raise KeyError(f'Missing source signal: {key}')
        signal = np.asarray(run['signals'][key], dtype=float)
        if signal.shape != (len(edges)-1,) or np.any(~np.isfinite(signal)) or np.any(signal < 0):
            raise ValueError(f'Invalid source signal at {mass:g} GeV')
        for mode, background in backgrounds.items():
            for threshold in spec['thresholds_keVee']:
                for exposure_tyr in spec['exposures_tonne_year']:
                    exposure_kg_day = exposure_tyr * 1000.0 * days
                    _, result = count_limit(signal, background, edges, exposure_kg_day,
                                            threshold, config, spec['statistics_method'])
                    mu = float(result['mu_upper'])
                    if mu < 0 or math.isnan(mu):
                        raise ValueError('Invalid statistical limit')
                    limit = mu * sigma_ref if math.isfinite(mu) else math.inf
                    rows.append(dict(mass_GeV=mass, operator=1, background=mode,
                         threshold_keVee=threshold, exposure_tonne_year=exposure_tyr,
                         exposure_kg_day=exposure_kg_day, limit_cm2=limit,
                         mu_upper=mu, roi_lower_keVee=result['roi_lower_keV'],
                         roi_upper_keVee=result['roi_upper_keV'],
                         retained_signal_fraction=result['retained_signal_fraction'],
                         reference_signal_events=result['reference_signal_events'],
                         background_events=result['background_events'],
                         status=result['status'], confidence_level=config['analysis']['confidence_level'],
                         method=spec['statistics_method']))
    expected = len(mass_grid) * 3 * 3 * 3
    assert len(rows) == expected
    out = Path(campaign_dir)
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out/'poster_limits.csv'
    with csv_path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = dict(schema_version=1, description='Poster SI-only statistical grid from one folded source run',
       source_run_directory=str(Path(run_directory).resolve()),
       source_run_id=run['manifest']['run_id'], source_numerical_id=run['manifest']['numerical_id'],
       source_git_commit=run['manifest']['source']['git_commit'],
       source_fingerprint=run['manifest']['source']['physics_fingerprint'],
       source_manifest_sha256=_sha256(Path(run_directory)/'manifest.json'),
       source_signal_sha256=_sha256(Path(run_directory)/'signals.npz'),
       source_background_sha256=_sha256(Path(run_directory)/'backgrounds.npz'),
       campaign_spec_sha256=_sha256(spec_path),
       assumptions=['Thresholds are post-folding observed-energy analysis cuts, not changes in detector efficiency.',
                    'Exposure scales fixed per-kg-day spectra; cosmogenic time dependence is not re-evaluated.',
                    'No/background-neutrinos-only/full share identical Ge target and detector response.',
                    'Neutrinos include solar_CEvNS and solar_ES from the configured benchmark.',
                    'Backgrounds treated as perfectly known; exact median Poisson CLs counting.',
                    'Provisional backgrounds, enriched composition, quenching and efficiency.'],
       threshold_keVee=spec['thresholds_keVee'], exposure_tonne_year=spec['exposures_tonne_year'],
       background_modes=list(backgrounds), n_masses=len(mass_grid), n_limits=len(rows),
       signal_reference_sigma_p_cm2=sigma_ref, limit_csv_sha256=_sha256(csv_path))
    _write_json(out/'poster_manifest.json', manifest)
    return rows, manifest


def load_campaign(campaign_dir):
    base = Path(campaign_dir)
    manifest = json.loads((base/'poster_manifest.json').read_text(encoding='utf-8'))
    if _sha256(base/'poster_limits.csv') != manifest['limit_csv_sha256']:
        raise ValueError('Poster CSV checksum mismatch')
    with (base/'poster_limits.csv').open(encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        for field in ('mass_GeV','threshold_keVee','exposure_tonne_year','exposure_kg_day',
                      'limit_cm2','mu_upper','reference_signal_events','background_events',
                      'roi_lower_keVee','roi_upper_keVee','retained_signal_fraction','confidence_level'):
            row[field] = float(row[field])
        row['operator'] = int(row['operator'])
    if len(rows) != manifest['n_limits']:
        raise ValueError('Incomplete poster limit table')
    return rows, manifest


def select_curve(rows, background, threshold_keVee, exposure_tonne_year):
    series = sorted((x for x in rows if x['background'] == background and
                     abs(x['threshold_keVee']-threshold_keVee) < 1e-10 and
                     abs(x['exposure_tonne_year']-exposure_tonne_year) < 1e-10),
                    key=lambda x: x['mass_GeV'])
    if not series:
        raise ValueError('Curve absent from campaign')
    x = np.array([x['mass_GeV'] for x in series])
    y = np.array([x['limit_cm2'] for x in series])
    y[~np.isfinite(y) | (y <= 0)] = np.nan  # gaps denote no sensitivity
    return x, y
