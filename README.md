# WIMP-LEGEND

Modular dark-matter recoil and poster-analysis package for germanium detectors.
`main` contains the consolidated post-poster development baseline. Poster
scenarios remain benchmarks, not validated LEGEND background predictions.

- [Branch consolidation and historical snapshots](docs/POSTER_CONSOLIDATION_2026-11-18.md)
- [RAPIDD comparison and independent tester handoff](docs/RAPIDD_COMPARISON.md)
- [Frozen NREFT conventions](docs/NREFT_CONVENTIONS.md)
- [Germanium nuclear-data provenance](docs/GE_NUCLEAR_DATA.md)
- [Poster analysis configuration](analyses/legend_poster/config.json)
- [Notebooks](notebooks/)

Install with `python -m pip install -e '.[test,workbench]'`.
Historical notebooks, datasets, outputs and poster plots are retained as
references; changing plot styling does not validate a physics assumption.
