# Validation scope

Version 0.2.1 is checked on Windows with MATLAB R2023b. Independent package
paths are added after `restoredefaultpath`; historical experiment folders are
not required. No full formal simulation suite was rerun for this layout revision.

| Current display and package check | Result |
|---|---|
| Existing public unit/regression suite | 8 passed, 0 failed |
| Original curve/image data | All 58 objects retained; maximum scaled CSV round-trip difference 3.915e-14 |
| Added grid bounds | Two horizontal BCRBs: 0.148783633769351 m and 0.0537883867232814 m |
| Grid / grouping / range | Grid has 2 panels; grouping and range are separate figures |
| Legends and titles | All 10 main-figure legends have opaque white backgrounds and borders; no overall plot titles |
| RF accuracy | 3 algorithm curves and 2 BCRBs retained |
| Iteration plot | 11.5 × 11.5 cm; 2 panels, 2 legends, iterations 1–10 |
| Source and reference data | 46 numerical files preserved byte for byte; summary and CDF CSV hashes unchanged |
| Extracted standalone package | All 8 reference figures redrawn from release paths; axes data match |
| Bundled reference | 306 summary rows, 8,008 CDF rows, 8 PDF figures; no MAT files |

`validation.json` separates this display audit from the historical v0.1 numerical
and orchestration evidence. Historical checks include 11 paired seed-1 cases,
86 position-MSE entries, bounds and 24 channel-NMSE entries matching exactly;
system regeneration; current/previous checkpoint recovery; a two-trial
interrupted/resumed run; skipping committed jobs; configuration mismatch and
partial-result rejection; and serial versus two-worker execution. The compact
two-trial checkpoint fixture used 7,259 bytes including its previous generation.

Historical simulations used one computational thread per process. The default
multi-threaded BLAS policy can change the rank-deficient M=5, Ns=9, RF=60 case.
Numerical entry points scope the historical single-thread policy and restore
the caller's setting; the existing regression suite checks both caller settings.
This is a reproducibility setting, not an algorithm repair.

The complete reference derives from the previously finished 181,000 trial jobs.
Selected-case regression does not establish bitwise portability on every future
machine or MATLAB release. No new full baseline speed benchmark, fixed-separation
experiment, or low-RF regularization experiment is included.

The 2026-10-05 display revision removes Wilson error bars from M1 and M5.
All 22 curves in those figures retain exactly the same X/Y data, colors,
line styles and markers; all 9 probability curves are ordinary marked lines.
Axis limits and opaque bordered legends are preserved. Reference CSV and all
46 numerical source hashes are unchanged. This style revision was checked by
redrawing both affected figures in Chinese and English, without estimator runs.
