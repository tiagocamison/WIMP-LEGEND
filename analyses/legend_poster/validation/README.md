# Validation evidence and provisional standard run

`existing-suite.log` records the complete pre-existing 889-test suite.
`workbench-suite.log` records the final 37 focused tests, including both actual
notebooks. These are non-overlapping collections (926 distinct tests), not one
926-test invocation. No pre-existing test was edited. `notebook_execution.json`
records the actual notebook source hashes, executed cell counts and output types.
Notebook execution used fresh IPython subprocesses because managed CI blocked
Jupyter network-interface discovery and ZMQ services.

`standard/` is a **provisional benchmark calculation**, not a physical LEGEND
prediction. It is generated from the master configuration's standard profile
using implementation commit `a819570f0f0669ab7b286c3dc77fd36a05cf77c4`.
It includes the complete resolved input and all code/data hashes in its manifest.
The original pilot assets are elsewhere and were not used as revised signals.

Arrays contain integrated events/(kg day) in their stored nonuniform observed
bins. The CSV has SI limits in cm^2 and O8/O11 coefficient limits in GeV^-2;
it also records counts, exposures, ROI edges, retained fractions and the exact
statistical construction. Infinities, if present, mean no finite sensitivity,
not zero cross section. `summary.json` gives convergence and representative
limits; it is derived from those files, not a second physical configuration.

Reproduce from the repository root in the documented editable environment:

```bash
python -m analyses.legend_poster --profile standard
python -m analyses.legend_poster --plot-only analyses/legend_poster/validation/standard
```

The second command reads the included numerical outputs and current plot styles;
it performs no physics integration. Original figure hashes are in the saved plot
manifest; regenerated PDFs can have different creation timestamps. To preserve
the distributed snapshot, pass a separate directory directly to `plot_run`:

```python
from analyses.legend_poster import load_run
from analyses.legend_poster.plotting import plot_run
run = load_run('analyses/legend_poster/validation/standard')
plot_run(run, output_dir='/tmp/wimp-poster-figures')
```

Inspect the manifest's assumptions before using any figure publicly. The
enrichment, detector response, activities and LAr transfer/veto remain provisional;
the seven-mass grid is intended for experimentation and should be refined for a
final poster. External comparison curves are frozen upstream snapshots with
their own physical/statistical assumptions, not a certified latest-2026 inventory.
