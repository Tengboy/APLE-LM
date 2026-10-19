# Reproduction protocol

The bundled reference is the accepted 1,000-trial thesis result set. The code
regenerates inputs from seeds; reference CSVs are used only by the plotter.

## Common physical and numerical configuration

| Parameter | Value |
|---|---|
| Frequency / wavelength / element spacing | 6 GHz / 0.05 m / 0.025 m |
| Default array | 45 × 45 elements; M = 3 subarrays per axis, Ns = 15 elements per subarray axis |
| Default RF chains | 160 |
| Users | 3, except the single-user spatial scan |
| Angular sector / radial interval | 120 degrees / 5–10 m unless varied |
| Array normal | First Cartesian axis; elements lie in the other two axes |
| Angular grid parameter / radial layers | 60 / two endpoint radii unless varied |
| Initialization limits | APLE-LM 20, ACM 20, ES-GA 200, P-SIGW 0 |
| Coarse single-user / joint budgets | 20 / 200 |
| Fine budget / outer budget | 20 / 200 |
| Initialization / coarse / fine tolerances | 1e-3 / 1e-4 / 1e-5 |
| Damping | 0.7, historical floor 0.01 in the experiment adapter |
| Random stream | `mt19937ar`, `NormalTransform='Ziggurat'` |
| Full trial seeds | 1:1000 |
| Computational threads per numerical call | 1; caller setting restored on return |

M is the number of subarrays **per axis**: there are M² blocks. Ns is the number
of elements **per axis in each block**: the physical array has (M Ns)² elements.
The nominal grid parameter controls direction-cosine spacing `2/grid`; it is
not the total number of retained 3-D grid points. The circular angular mask and
two radii determine the actual dictionary size. `aple_cases` is the canonical
configuration listing.

P-SIGW's zero initialization limit disables only the initializer's local
refinement. Grid selection and the subsequent joint P-SIGW iteration remain.
The historical P-SIGW adapter can execute 201 iterations under its `> 200`
exit test; this has not been silently changed.

## Experiments and exact random streams

| Family / figures | Cases | Operating points | W seed |
|---|---:|---|---|
| `snr` / M2, M3 | 1 | SNR = 5:3:20 dB; all four methods | 0 |
| `grid` / M1 | 12 | M = 3, 5; Ns = 15; grid = 15:15:90; SNR = 15 dB | 0 |
| `aperture` / M4 | 5 | (M,Ns) = (3,15),(4,15),(5,15),(3,20),(3,25); SNR = 5:3:20 | 0 |
| `range` / M6 | 3 | 5–10, 10–15, 15–20 m; SNR = 5:3:20 | 0 |
| `rf` / M5 | 15 | (3,15),(5,9),(4,15); RF = 60:30:180; SNR = 15 dB | 1:5 for the five RF values |
| `iterations` / M8 | 2 | 5–10 and 10–15 m; SNR = 14 dB; 20 true updates | 0 |
| `coverage` / M7 | 143 | angle = 0:5:60 degrees; radius = 5:15 m; SNR = 15 dB; K = 1 | 0 |

The grid family compares APLE-LM, ES-GA, and P-SIGW. The other families use
APLE-LM only, except `snr`. Array grouping experiments use grid = 4 Ns.
Coverage is `ndgrid(angle,radius)`, with angle varying fastest, and position
`radius * [cosd(angle); sind(angle); 0]`.

For SNR sweeps, all six noise vectors are drawn **before** estimation. M8
regenerates this same six-SNR sequence and selects its fourth observation
(14 dB); drawing a fresh single 14 dB noise vector would not reproduce it.
Grid, RF, and coverage experiments use a single 15 dB draw. W has its own local
random stream, independent of trial seeds. Baselines reorder W columns and
element positions together to use the same physical observation as APLE-LM.
Different grouping configurations do not share a guaranteed identical physical W.

User generation samples polar angle and radius uniformly, not uniformly in
volume or solid angle. It retains the historical separation `2/Ns` in direction
cosines and the OR rule between the two axes. Consequently, Ns=9/20/25 use
different user populations from Ns=15. A common fixed-separation experiment
would be a new experiment, not a relabeling of these reference results.

## Metrics

For each trial, a one-to-one assignment minimizes the sum of squared user
position errors. Let `e² = mean_k ||p_hat_k - p_k||²` after this assignment.

- Full position RMSE: `sqrt(mean(e²))`.
- Conditional RMSE: `sqrt(mean(e²(e² <= 3²)))`, using finite valid samples.
  An empty accepted set gives NaN, not zero. Errors are selected, not clipped.
- Exceedance: count of `e² > 3²` or nonfinite/invalid estimates divided by
  attempted trials. M1/M5 show probability curves without interval bars. Pointwise 95% Wilson
  bounds remain in the summary CSV for reference.
- CDF: empirical distribution of `sqrt(e²)` with attempted trials as denominator.
  Invalid estimates leave missing probability mass at finite error; they are not
  silently discarded. The bundled reference has no nonfinite estimates.
- Any invalid estimate makes the corresponding **full** RMSE NaN. Its attempt
  and failure remain recorded. Numerically finite extreme errors remain valid
  observations, and they still contribute to full RMSE and exceedance.
- Channel NMSE is calculated for RF-combined per-user channels. Normalize each
  user's squared error by that user's true channel energy, average users and
  trials, then take `10*log10`. Channel assignment uses the position assignment.

The plot label is **BCRB**, following the paper/accepted figures. The implemented
quantity averages bounds obtained from each trial's Fisher information without
an explicit prior-information term. Position bounds are divided by the number
of users before averaging and square rooting; channel bounds use the same
per-user energy normalization as channel NMSE. BCRB curves use all samples;
they are references, not lower bounds derived for the selected <=3 m subset.

M4 also uses the common 3 m selection. Its `Included`, `Excluded`, full RMSE,
and exceedance columns remain available even though a separate probability
panel is not shown. No method-specific cleaning rule is applied.

## Running only what is needed

```matlab
% Case order is visible without building a dictionary:
d = aple_cases('rf');
struct2table(d.cases)

% Bound the work done in this call, then repeat to continue.
out = fullfile(pwd,'runs','my_rf_run');
run_experiments('M5','full','OutputDir',out,'Trials',1000, ...
                'BatchSize',25,'MaxBatches',4);
run_experiments('M5','full','OutputDir',out,'Trials',1000);

% A small selected-case run. This is deliberately incomplete.
out = fullfile(pwd,'runs','rf_check');
run_experiments('M5','full','OutputDir',out,'Seeds',[1,2], ...
                'Cases',[1,6,11],'BatchSize',1);
plot_results(out,'Figures',{'M5'},'AllowPartial',true);
```

`Seeds` overrides `Trials`. `Cases` applies to every selected family, so prefer
one family when selecting explicit indices. Profile `smoke` uses one seed and
representative cases (see `d.smoke_cases`); `full` uses every case and 1,000 seeds
unless overridden. Repeating a command skips committed jobs, including recorded
numerical failures; changing estimator behavior or retrying with changed code
requires a separate output folder. Plotting with `AllowPartial=true` is for
inspection and does not impute missing cells or count unfinished jobs as failures.

Parallelism operates over a bounded batch of independent trial seeds. The
numerical calls use one BLAS/computational thread per process, matching the
historical experiments; this applies to the client and every worker. The
default is serial; two workers are a conservative starting point. A worker
receives one active system/dictionary, so RAM use grows with worker count.
An already running pool must have the requested size. No per-worker/per-trial
result files are created. MATLAB may maintain its own parallel-pool job files
outside this experiment runner; those are controlled by the MATLAB cluster
profile, not stored as experiment outputs here.

Numerical BLAS differences, finite-difference Hessians, line searches, and
ill-conditioned stress cases can change individual trial outputs across
hardware, MATLAB releases, or thread settings. Fixed seeds specify input
generation; they do not guarantee bitwise identical floating-point estimates.
