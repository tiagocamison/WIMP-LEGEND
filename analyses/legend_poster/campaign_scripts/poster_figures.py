"""Plot-only LEGEND SI campaign figures using original DirectDetectionPlots styling.

All LEGEND curves come from verified campaign data. External experiment/projection
curves are rendered with DirectDetectionPlots.PlotFuncs.PlotBound, exactly as in
DM-Nucleon_LEGEND.ipynb (not digitized or replaced). No physics is recomputed.
"""
from __future__ import annotations

from contextlib import nullcontext
from pathlib import Path
import hashlib
import importlib.util
import os

import matplotlib.pyplot as plt
import numpy as np

from ..analysis import load_run
from .poster_campaign import load_campaign, load_spec, select_curve


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def locate_ddp(repo_root, ddp_root=None):
    """Resolve a full DirectDetectionPlots checkout, never an invented data path."""
    root = Path(repo_root).resolve()
    candidates = ([Path(ddp_root)] if ddp_root is not None else [
        Path(os.environ['DIRECT_DETECTION_PLOTS_ROOT'])
        if os.environ.get('DIRECT_DETECTION_PLOTS_ROOT') else None,
        root / 'DirectDetectionPlots',
        root.parent / 'DirectDetectionPlots',
        Path.cwd() / 'DirectDetectionPlots',
        Path.cwd(),
    ])
    for p in candidates:
        if p is None:
            continue
        p = p.expanduser().resolve()
        if p.name == 'notebooks' and (p / 'PlotFuncs.py').is_file():
            p = p.parent
        if (p / 'notebooks/PlotFuncs.py').is_file() and \
           (p / 'notebooks/sty.mplstyle').is_file() and \
           (p / 'data/limits/Nucleon/SI/Projections/CDEX.txt').is_file():
            return p
    raise FileNotFoundError(
        'DirectDetectionPlots checkout not found. In poster_plots_campaign.ipynb set '
        'DDP_ROOT = Path(r"<path-to-DirectDetectionPlots>") (the directory containing '
        'notebooks/PlotFuncs.py, notebooks/sty.mplstyle and data/), or set '
        'DIRECT_DETECTION_PLOTS_ROOT. No older comparison curves are substituted.'
    )


def _ddp_plotbound(ddp):
    """Import the user's original PlotFuncs.py without changing working directory."""
    p = ddp / 'notebooks/PlotFuncs.py'
    module_spec = importlib.util.spec_from_file_location('legend_poster_ddp_plotfuncs', p)
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    return module.PlotBound, module.col_alpha


def _draw_ddp_limits(ax, ddp, font=17, *, show_fog=True):
    """Reproduce the user's low-mass SI backgrounds and projection styling.

    Colors, datasets, labels, positions, and projection dashed lines follow
    cells 5 and 6 of DM-Nucleon_LEGEND(1).ipynb. DDP data are NEVER substituted.
    """
    bound, col_alpha = _ddp_plotbound(ddp)
    base = ddp / 'data/limits/Nucleon/SI'
    required = {
        'LZ 2026': base / 'LZ_2026_Combined.txt',
        'XENONnT 2026': base / 'XENONnT_2026_Combined.txt',
        'PandaX 2025 low mass': base / 'PandaX_2025_low_mass.txt',
        'SuperCDMS projection': base / 'Projections/SuperCDMS.txt',
        'CDEX-50 projection': base / 'Projections/CDEX.txt',
    }
    if show_fog:
        required['Xe neutrino fog'] = ddp / 'data/fog/Nucleon/SI/nuFog_n2_Xe.txt'
    missing = [f'{k}: {p}' for k, p in required.items() if not p.is_file()]
    if missing:
        raise FileNotFoundError('Missing required DirectDetectionPlots data:\n' + '\n'.join(missing))

    plt.sca(ax)  # PlotFuncs.PlotBound uses plt.text to place inline annotations.
    # Same datasets/shaded style as the supplied notebook; labels are annotations,
    # not legend entries, avoiding duplicate external-curve legend entries.
    bound(ax, str(required['LZ 2026']), edgecolor='darkgreen', textcolor='darkgreen',
          facecolor=[0., .62, .38], facealpha=.25, label='LZ', rotation=-40,
          label_pos=[7., 1.1e-46], fontsize=font, lw=2, path_effects=None,
          text_path_effects=None, zorder=0.5)
    bound(ax, str(required['XENONnT 2026']), edgecolor='teal', textcolor='teal',
          facecolor=col_alpha('teal', .5), facealpha=.25, label='XENONnT',
          rotation=-30, label_pos=[5., 2.e-46], fontsize=font, lw=2,
          path_effects=None, text_path_effects=None, zorder=.4)
    bound(ax, str(required['PandaX 2025 low mass']), edgecolor='navy',
          textcolor='navy', facecolor=col_alpha('navy',.5), facealpha=.25,
          zorder=-1.1, lw=2, path_effects=None)
    if show_fog:
        bound(ax, str(required['Xe neutrino fog']), y2=0., label='Neutrino fog\n(Xe)',
              label_pos=[1.5, 1.5e-46], fontsize=font, rotation=0,
              facecolor='#a7b9d6', lw=0, textcolor='#859bbc',
              facealpha=.5, text_path_effects=None, zorder=-2)
    bound(ax, str(required['SuperCDMS projection']), edgecolor='orange',
          textcolor='orange', facealpha=0, label='SuperCDMS', linestyle='--',
          rotation=-15, fontsize=font, label_pos=[1.5, 6.e-45], end_x=-10,
          path_effects=None, text_path_effects=None, lw=2.4, zorder=2)
    bound(ax, str(required['CDEX-50 projection']), edgecolor='royalblue',
          textcolor='royalblue', facealpha=0, label='CDEX', linestyle='--',
          rotation=-45, fontsize=font, label_pos=[3., 5.e-44],
          path_effects=None, text_path_effects=None, lw=2.4, zorder=2)


def _draw_frozen_limits(ax, repo_root):
    """Explicit older-snapshot option for environments lacking DDP (not default)."""
    folder = Path(repo_root) / 'analyses/legend_poster/data/si_limits'
    files = [
        ('SuperCDMS SNOLAB (historical projection)', 'SuperCDMS_SNOLAB.txt', 'orange', '--'),
        ('CDMSlite (historical)', 'CDMSLite.txt', 'royalblue', '--'),
        ('LZ 2024 (historical)', 'LZ_2024.txt', 'darkgreen', '-'),
        ('XENONnT 2025 (historical)', 'XENONnT_2025.txt', 'teal', '-'),
    ]
    for name, filename, color, linestyle in files:
        path = folder / filename
        if not path.is_file():
            raise FileNotFoundError(f'Frozen comparison file unavailable: {path}')
        data = np.loadtxt(path, comments='#')
        if data.ndim != 2 or data.shape[1] != 2 or np.any(~np.isfinite(data)) or np.any(data <= 0):
            raise ValueError(f'Invalid comparison dataset: {path}')
        ax.plot(data[:, 0], data[:, 1], color=color, ls=linestyle, lw=1.8,
                alpha=.8, label=name, zorder=1)


def _limit_axes(ax, *, xlim, ylim, label_size, tick_size):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(r'Dark matter mass [GeV/$c^2$]', fontsize=label_size)
    ax.set_ylabel(r'SI DM-nucleon cross section [cm$^2$]', fontsize=label_size)
    ax.tick_params(which='major', labelsize=tick_size)
    ax.grid(False)


def _draw_curve(ax, rows, bg, threshold, exposure, *, label, color, lw=2.8):
    m, limit = select_curve(rows, bg, threshold, exposure)
    m, limit = np.asarray(m, dtype=float), np.asarray(limit, dtype=float)
    good = np.isfinite(m) & np.isfinite(limit) & (m > 0) & (limit > 0)
    # Convert invalid points to NaN instead of connecting across an inaccessible
    # mass range. Never interpolate/alter a physical sensitivity curve.
    yy = np.where(good, limit, np.nan)
    ax.plot(m, yy, color=color, linewidth=lw, label=label, zorder=10,
            solid_capstyle='round')


def _save(fig, folder, name, dpi):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(folder / f'{name}.{ext}', dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def plot_five(campaign_dir, run_dir, spec_path, repo_root, *,
              ddp_root=None, xlim=(1., 10.), ylim=(1e-49, 1e-42), dpi=220,
              figure_dir=None, comparison_mode='ddp', figure_size=(11.5, 7.5),
              label_size=28, tick_size=17, annotation_size=16,
              legend_size=13, show_fog=True, show_titles=False):
    """Export five PDF/PNG figures; no physics or statistical recomputation.

    comparison_mode: 'ddp' (2026 datasets + shaded PlotBound), 'frozen'
    (explicit older snapshots), or 'none'. The external root is validated before
    saving any figure. The source arrays and campaign table are checksum-checked.
    """
    rows, manifest = load_campaign(campaign_dir)
    spec = load_spec(spec_path)
    run = load_run(run_dir)
    if manifest['source_run_id'] != run['manifest']['run_id']:
        raise ValueError('Campaign numerical run does not match plotting source')
    if manifest['source_signal_sha256'] != _sha(Path(run_dir) / 'signals.npz') or \
       manifest['source_background_sha256'] != _sha(Path(run_dir) / 'backgrounds.npz'):
        raise ValueError('Campaign/source checksum mismatch')
    if manifest['campaign_spec_sha256'] != _sha(spec_path):
        raise ValueError('Campaign specification changed since limits were computed')
    if comparison_mode not in {'ddp', 'frozen', 'none'}:
        raise ValueError("comparison_mode must be 'ddp', 'frozen', or 'none'")
    ddp = locate_ddp(repo_root, ddp_root) if comparison_mode == 'ddp' else None
    if ddp is not None:
        # Check EVERY dataset before generating figures rather than silently
        # falling back to another year's limits.
        data_base = ddp / 'data/limits/Nucleon/SI'
        required = [data_base / name for name in (
            'LZ_2026_Combined.txt', 'XENONnT_2026_Combined.txt',
            'PandaX_2025_low_mass.txt', 'Projections/SuperCDMS.txt',
            'Projections/CDEX.txt')]
        if show_fog:
            required += [ddp / 'data/fog/Nucleon/SI/nuFog_n2_Xe.txt']
        missing = [p for p in required if not p.is_file()]
        if missing:
            raise FileNotFoundError('Missing comparison files:\n' + '\n'.join(map(str, missing)))
    out = Path(figure_dir or Path(campaign_dir) / 'figures')
    c = run['manifest']['resolved_configuration']
    scenario = spec['input_scenario']
    bins = run['backgrounds']['observed_edges_keV']
    widths = np.diff(bins)
    centers = (bins[:-1] + bins[1:]) / 2
    year = c['analysis']['days_per_year']
    sigma_ref = c['signal']['sigma_p_ref_cm2']

    style = ddp / 'notebooks/sty.mplstyle' if ddp is not None else None
    context = plt.style.context(str(style)) if style is not None else nullcontext()
    with context, plt.rc_context({
        'text.usetex': False, 'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
        'axes.unicode_minus': False,
        'legend.frameon': True,
    }):
        # 1. Observed spectra (kept as bin-averaged densities from saved arrays).
        fig, ax = plt.subplots(figsize=figure_size)
        vals = run['backgrounds']
        ec = sum((np.asarray(vals[f'{scenario}__{name}'])
                  for name in ('Ge68', 'Ga68', 'Zn65', 'Fe55')), np.zeros_like(widths))
        groups = [
            ('EC lines (Ge, Ga, Zn, Fe)', ec, '#ad6b9c'),
            ('Tritium', vals[f'{scenario}__H3'], '#d0a048'),
            ('Argon-39', vals[f'{scenario}__Ar39'], '#4c9c88'),
            ('Compton continuum', vals[f'{scenario}__flat_observed'], '#777777'),
            ('Solar CEvNS', vals[f'{scenario}__solar_CEvNS'], '#3876aa'),
            ('Solar ES', vals[f'{scenario}__solar_ES'], '#7da3cd'),
        ]
        for label, per_bin, color in groups:
            d = np.asarray(per_bin) / widths * year
            if np.any(d > 0):
                ax.plot(centers, np.where(d > 0, d, np.nan), color=color,
                        linewidth=1.8, label=label)
        full = vals[f'{scenario}__Total'] / widths * year
        ax.plot(centers, np.where(full > 0, full, np.nan), lw=2.7,
                color='black', label='Total background', zorder=8)
        for color, mass in zip(('#cf534f','#ab232b','#761727'), spec['spectrum_masses_GeV']):
            key = f'{scenario}__O1__m{mass:g}'
            if key not in run['signals']:
                raise KeyError(f'Missing required WIMP spectrum: {key}')
            d = run['signals'][key] / widths * year * spec['spectrum_sigma_p_cm2'] / sigma_ref
            ax.plot(centers, np.where(d > 0, d, np.nan), color=color, lw=2.4,
                    ls='--', label=f'SI WIMP, {mass:g} GeV')
        ax.axvline(.1, color='gray', lw=1.5, ls=':')
        ax.set(xlim=(.04,2.), yscale='log')
        ax.set_ylim(1e-7, 1e3)
        ax.set_xlabel('Observed energy [keVee]', fontsize=label_size)
        ax.set_ylabel('Events / (kg·yr·keVee)', fontsize=label_size)
        ax.tick_params(which='major', labelsize=tick_size)
        ax.legend(loc='upper right', fontsize=max(9,legend_size-1), ncol=2)
        ax.text(.03,.03, f'Provisional Ge benchmark; WIMP $\\sigma_p={spec["spectrum_sigma_p_cm2"]:.0e}$ cm$^2$',
                fontsize=11, transform=ax.transAxes)
        if show_titles:
            ax.set_title('Energy spectra')
        _save(fig,out,'01_energy_spectra',dpi)

        # 2--5. Same external experimental boundaries in every sensitivity plot.
        specs = [
            ('02_neutrino_best_case', 'Neutrino-only benchmark | 100 eVee, 10 t·yr',
             [('neutrinos', .1, 10., 'LEGEND: neutrinos only, 10 t·yr', '#881d26')],
             '100 eVee · enriched Ge'),
            ('03_background_comparison', 'Background assumptions | 100 eVee, 10 t·yr',
             [('none', .1, 10., 'No background', '#404040'),
              ('neutrinos', .1, 10., 'Solar neutrinos only', '#dc8a23'),
              ('full', .1, 10., 'Full background', '#891e2a')],
             '100 eVee · 10 t·yr · enriched Ge'),
            ('04_threshold_comparison', 'Threshold assumptions | 10 t·yr, full background',
             [('full', .1, 10., '100 eVee', '#d46f65'),
              ('full', .15, 10., '150 eVee', '#c42c39'),
              ('full', 1., 10., '1 keVee', '#781825')],
             '10 t·yr · enriched Ge · full background'),
            ('05_exposure_comparison', 'Exposure assumptions | 150 eVee, full background',
             [('full', .15, .1, '0.1 t·yr', '#d46f65'),
              ('full', .15, 1., '1 t·yr', '#c42c39'),
              ('full', .15, 10., '10 t·yr', '#781825')],
             '150 eVee · enriched Ge · full background'),
        ]
        for filename, title, curves, legend_title in specs:
            fig, ax = plt.subplots(figsize=figure_size)
            _limit_axes(ax,xlim=xlim,ylim=ylim,label_size=label_size,tick_size=tick_size)
            if comparison_mode == 'ddp':
                _draw_ddp_limits(ax,ddp,annotation_size,show_fog=show_fog)
            elif comparison_mode == 'frozen':
                _draw_frozen_limits(ax,repo_root)
            for bg, threshold, exp, label, color in curves:
                _draw_curve(ax,rows,bg,threshold,exp,label=label,color=color)
            if show_titles:
                ax.set_title(title)
            handles, labels = ax.get_legend_handles_labels()
            # DDP labels are inline; only LEGEND runs appear in the legend.
            if comparison_mode == 'ddp':
                handles = handles[-len(curves):]
                labels = labels[-len(curves):]
            ax.legend(handles, labels, title=legend_title,
                      title_fontsize=max(9,legend_size-1), fontsize=legend_size,
                      loc='upper right', framealpha=.93)
            _save(fig,out,filename,dpi)
    return out
