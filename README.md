# APLE-LM

[中文说明](README_zh.md) · [MIT license](LICENSE) · [Citation](CITATION.cff)

MATLAB implementation of near-field multiuser localization with an extremely
large array and limited RF chains. This release reproduces the eight selected
thesis figures, using the audited implementation associated with:

Boyu Teng, Xiaojun Yuan, Rui Wang, Ying-Chang Liang, and Xinming Huang,
“Near-Field Multiuser Localization Based on Extremely Large Antenna Array With
Limited RF Chains,” IEEE Transactions on Wireless Communications, 2025.
[DOI: 10.1109/TWC.2025.3578399](https://doi.org/10.1109/TWC.2025.3578399)

The reference figures include genuine low-grid and low-RF failures. This is an
audited reproduction with documented numerical and evaluation changes; it is
not a verbatim archive of the original plotting code.

Version 0.2 renumbers the main figures in journal order and separates array
grouping from distance. See [the old/new mapping](docs/FIGURES.md). Use semantic
figure names (`grid`, `snr`, `cdf`, `aperture`, `rf`, `range`, `coverage`,
`iterations`) in scripts that should remain independent of figure numbering.

## Requirements

- MATLAB R2023b was used for validation. R2021a or newer is the compatibility
  target; other releases have not been tested.
- Base MATLAB and its JVM are sufficient for serial execution and plotting.
- Parallel Computing Toolbox is optional (`UseParallel=true`).
- Chinese plots need an installed Chinese font; English is the default.
- No Python, historical MAT files, other source repository, or external dataset
  is needed by users of this package.

Extract the complete folder, set it as MATLAB's current folder, and run:

```matlab
setup_aple;
plot_results('reference');                    % all eight figures, no simulation
plot_results('reference','Language','zh');    % Chinese axes and legends
```

PDF and PNG exports are written to `runs/reference_plots/`. The accepted
Chinese PDFs are already in `reference/figures/`. Plots have no overall titles;
experiment captions and interpretation are in [FIGURES.md](docs/FIGURES.md).

## Run simulations

```matlab
% Quick end-to-end check: one seed for the default array, all six SNRs,
% all four algorithms, position and channel errors (M2 and M3).
out = run_experiments('M2','smoke');
plot_results(out,'Figures',{'M2','M3'});

% Selected cases for every family, one seed each; not statistical evidence.
run_experiments('all','smoke');

% Full protocol: 1,000 seeds for every case, all original iteration budgets.
out = run_experiments('all','full');
plot_results(out);

% Optional parallel execution, bounded batches, automatic resume.
out = run_experiments('M5','full','UseParallel',true,'Workers',2,'BatchSize',25);
plot_results(out,'Figures',{'M5'});
```

Repeating the same call resumes unfinished jobs and skips completed jobs.
Use a **new OutputDir when changing seeds, numerical source, or MATLAB release**;
the runner checks a portable SHA-256 identity before resuming. Selection of
cases, batch size, or worker count can change without changing experiment
identity. Keep the same profile/seeds when extending a partial run.

Full reproduction is a substantial Monte Carlo computation: 181,000 trial/case
jobs, including the 143,000-job spatial scan. It does not run automatically when
plotting the bundled results. `smoke` changes only case/trial selection, not the
algorithm's iteration limits. Partial runs cannot be mistaken for a full result:
metadata reports completion and the plotter rejects incomplete/failed datasets
unless `AllowPartial=true` is explicitly supplied for inspection.

See [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the exact protocol, case
selection, bounded runs, and the historical separation rule.

## Files and storage

| Folder/file | Purpose |
|---|---|
| `run_experiments.m`, `plot_results.m`, `setup_aple.m` | Public entry points |
| `src/` | 46 unchanged numerical dependency files |
| `experiments/` | Seed generation, configurations, compact checkpoints, statistics |
| `plotting/` | Common M1–M8 figure styling, English/Chinese labels |
| `reference/` | CSV curve/CDF values, metadata, eight accepted PDF figures |
| `docs/` | Protocol, statistics, implementation changes, storage, validation |
| `examples/` | Small runnable scripts |
| `tests/` | API, statistics, configuration, resume/recovery and lock tests |
| `runs/` | Generated files; ignored by Git, not distributed |

**No experiment MAT file is bundled.** At runtime, each family stores a small
numeric checkpoint and one recovery copy. Dictionaries, observations, channels,
positions, Hessians, and solver traces stay in memory and are discarded after
their scalar errors are computed. There are no per-trial files, saved dictionary
caches, copied source trees, or growing solver logs. Details and recovery steps:
[STORAGE.md](docs/STORAGE.md).

## Metrics and interpretation

M2, M3, M6, M8, and M7 use complete samples. M1, M4, and M5 accuracy curves
use the declared **3 m trial-RMS threshold**. M1/M5 also show full-denominator
exceedance probabilities without interval bars. Full RMSE, retained counts,
and failure counts remain in `summary.csv`. No estimate is clipped to 3 m.
The common threshold is not tuned separately by method. BCRB is the figure label
for the implemented bound; see the precise convention in the protocol.

The five M4 grouping configurations are retained. M5 contains three APLE-LM
curves and two BCRBs. M8 displays the first ten of twenty actual updates in two
vertically stacked panels. M7 uses the original 13 × 11 scan cells without
interpolation.

[NUMERICAL_NOTES.md](docs/NUMERICAL_NOTES.md) distinguishes inherited numerical
corrections, evaluation changes, equivalent computation optimizations, and this
release's storage changes. In particular, the historical `2/Ns` user-separation
rule is retained to reproduce existing data. Ns = 9, 20, and 25 results must not
be described as using a common fixed separation threshold. No new low-RF cure or
end-to-end speed advantage over the baselines is claimed by this release.

## Tests and license

```matlab
setup_aple;
results = runtests('tests/test_release.m');
assert(all([results.Passed]));
```

Validation scope is recorded in [VALIDATION.md](docs/VALIDATION.md). Numerical
source provenance and hashes are in [NUMERICAL_SOURCES.json](docs/NUMERICAL_SOURCES.json).
Licensed under MIT, copyright © 2026 **Boyu Teng**. Original author comments are
preserved; see [AUTHORS.md](AUTHORS.md).
