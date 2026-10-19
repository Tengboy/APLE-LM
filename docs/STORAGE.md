# Storage and recovery

## What is saved

Each requested experiment family has one `<family>.mat` containing a single
`state` struct. It contains:

- Canonical case definitions, seeds, MATLAB release, schema and source signature.
- Per-trial position MSE and BCRB scalars; `snr` additionally stores channel NMSE
  and channel-bound scalars. M8 stores 20 scalar errors per trial.
- `done` and `failed` masks, generation count, uint16 iteration counts, single
  precision elapsed seconds, and at most ten failure examples of 300 characters.

The package does not save true/estimated positions, gains, observations,
physical channels, dictionaries, W matrices, Hessians, detailed trajectories,
or copies of source files. All these are regenerated or computed in RAM. One
system is active at a time; all coverage cells reuse the same dictionary.
Elapsed seconds exclude shared dictionary preparation, so this field is not an
end-to-end runtime comparison.

Full default numeric state for all eight families is about 6.3 MiB before MAT
compression and small metadata overhead. A current and previous state therefore
need about 13 MiB before compression; completed-run CSVs and exported figures
are additional. This is an allocation estimate, not a claimed measurement of
a fresh 181,000-job full run. Larger custom seed lists scale storage linearly.

## Commit and resume

A batch defaults to 25 trial/case jobs. After a complete batch:

1. Write `<family>.mat.pending.mat` and load it back to verify its content.
2. Copy a readable current checkpoint to `<family>.mat.previous.mat`.
3. Replace the current checkpoint with the verified pending file.

There are at most three checkpoint files for a family during a write and two
after normal completion. A corrupt current checkpoint is never copied over a
valid previous checkpoint. Loading tries the current file, then the previous
file with a visible recovery warning. An abandoned pending file alone is not
treated as committed work. If neither complete checkpoint is usable, the runner
stops rather than silently resetting progress.

Repeat the same experiment call after interruption. Only committed `done` jobs
are skipped; the interrupted uncommitted batch is recomputed. Recovery from a
corrupt latest checkpoint rolls back to the prior committed generation.
Choose a smaller batch for more frequent commits; choose a larger batch to
reduce write frequency. There is no promise of filesystem-level durability
through every possible OS/power failure, so retain separate backups for valuable
long runs if desired.

The signature includes experiment definitions, ordered seeds, MATLAB release,
numerical/experiment source files, and the runner. It excludes absolute paths
and plot styling, allowing the entire folder to move without invalidating a run.
It does not hash every external runtime library or record all CPU/BLAS settings.
Do not mix modified numerical code into an old output directory.

## Output lock

An atomic directory creation at `<output>/.aple-lock` prevents two runners from
writing the same output folder concurrently. Plotting a checkpoint directory
acquires the same lock while deriving and reading its tables, so a running
simulation cannot race a plot-time CSV export. Normal completion/errors release
the lock. A forcibly killed process may leave it behind.

If MATLAB reports `APLE:RunLocked`, first verify that no MATLAB process or worker
is still running an experiment in that output folder. Only then remove the
empty `.aple-lock` directory and repeat the original command. Do not delete
checkpoint files to bypass a lock or identity error. Independent output
directories can be run separately if memory permits.

## Tables and publishing

`summary.csv`, `cdf.csv`, and `metadata.json` are derived from checkpoints at the
end of a call or before plotting a run folder. They can be regenerated; checkpoint
files remain the resume authority. CSV exports do not contain the full original
trial outputs, so the bundled reference can reproduce the selected figures but
cannot re-evaluate arbitrary future metrics or thresholds for every family.
Fresh runs retain scalar MSE arrays, allowing new threshold analyses without
rerunning estimators.

`reference/` contains only lightweight CSVs, JSON, and eight PDF figures.
Generated `runs/` data are Git-ignored and excluded from the release ZIP. Custom
output directories outside `runs/` should also be excluded before publishing.
`SaveFig=false` is the default; avoid shipping editable FIG files unnecessarily.
