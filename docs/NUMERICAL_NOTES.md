# Implementation notes and changes

This release preserves the audited numerical kernels, rather than reverting to
the first journal plotting scripts. It adds independent experiment execution and
compact storage around those kernels. `NUMERICAL_SOURCES.json` records the
original relative source locations and SHA-256 of each unchanged copied file.

The subarray coarse estimate, freely fitted block gains, geometry-coupled fine
estimate, Gaussian information fusion, and grid-dictionary precomputation are
existing method components. They are not new contributions introduced by this
packaging work.

## Inherited numerical corrections

1. **Hessian in the actual optimization coordinates.** Parameters are
   `q=(u,v,r)` and `p=r*[sqrt(1-u²-v²);u;v]`. The audited implementation centrally
   differences the checked analytic gradient, refits nuisance gains at perturbed
   points, and symmetrizes the result. This includes the second-order coordinate
   chain term omitted by a simple `J' * H_p * J` conversion away from a stationary
   point. It also avoids inconsistent amplitude exponents and mixed product
   derivatives in the old handwritten Hessians. `hessian_from_gradient.m` limits
   perturbations to the physical domain. Symmetrization does not make an
   indefinite Hessian positive definite. This change is not an algebraically
   equivalent speed optimization.
2. **Mean and curvature at the same point.** Joint coarse routines recalculate
   curvature at the final returned parameters. In fine routines, `cal_Hess=1`
   evaluates curvature without taking another discarded position update. This
   prevents fusing a mean with curvature evaluated at a different point.
3. **Candidate feasibility and descent checks.** Line searches reject angles
   outside `u²+v²<1` and nonpositive distances. Newton candidates must be real,
   finite, feasible, and decrease the objective. Zero-direction outputs are
   initialized. The historical sign-direction and strict-decrease search remain;
   this is not a wholesale replacement with a new standard Armijo scheme.
4. **Initializer output synchronization.** The initializer returns Cartesian
   positions reconstructed from its final parameters, including when it finishes
   after an angular step. It no longer needs true positions just for allocation.
5. **General user count.** Gain phase generation uses `source_num` rather than a
   hard-coded three-user size. K=3 keeps its original random draw count; K=1 is
   used for M7.
6. **Stopping and numerical status.** The APLE outer loop respects its 200-update
   limit and checks final outputs for finite real positions. Curvature diagnostics
   do not alter the estimate. P-SIGW retains its separate historical stopping
   convention. A large finite error is not converted into a failed computation.

## Inherited evaluation corrections

7. **One-to-one user matching.** `match_user_positions.m` minimizes total squared
   position error using a toolbox-free assignment algorithm. Repeated nearest
   neighbors could previously hide a distant estimated user by assigning another
   estimate more than once. This is an evaluation correction, not a localization
   algorithm improvement. The same position assignment orders channel estimates.
8. **Per-user channel-bound normalization.** Both channel error and its bound are
   normalized per user's true RF-combined channel energy before averaging users.
   A ratio of total variances to total powers is generally a different quantity.
   The implemented channel is `W*h`, not the uncombined physical channel `h`.

## Inherited equivalent computation optimizations

9. **Factor fixed noise covariance once.** `prepare_grid_dictionary.m` uses a
   Cholesky factor of `W*W'` across grid points instead of repeatedly factoring
   the same matrix. Linear solves preserve the weighted projections.
10. **Preallocation and batched operations.** Page-wise matrix multiplication
    constructs block responses; vectorized grid scores and preallocated arrays
    reduce repeated allocation. Complex conjugation and weighting are retained.
11. **Reuse fixed responses.** Response and weighted matrices that are fixed in a
    gain/fine stage are reused. The implemented schedule computes coarse position
    once, then iterates gain/fine updates; it does not rerun coarse localization
    every outer iteration.
12. **Cache numerical-FIM base values.** Forward differences reuse the unperturbed
    function value across parameter columns. Calls requesting only the position
    information can skip the channel Jacobian. Step size and parameter ordering
    remain unchanged.
13. **Propagate only required covariance diagonals.** Channel variance uses the
    diagonal of `T*C*T'` without allocating the full propagated covariance.
    Position-only bound computation solves only required columns. These reduce
    local work/storage, but do not establish faster end-to-end localization than
    a baseline with an offline dictionary.

## Experiment adaptations retained in this release

14. **Same physical observation for baselines.** Full-array baselines reorder W
    columns and element positions together. Matching the random seed alone would
    not guarantee identical physical combining weights under different element
    storage orderings.
15. **Method-specific initialization budgets.** APLE/ACM/ES-GA/P-SIGW use
    20/20/200/0 initialization refinement limits. This preserves the corrected
    baseline protocol. P-SIGW's subsequent joint update remains enabled.
16. **Explicit physical/sampling protocol.** Default RF is 160, and all distance
    sweeps use 160 RF with 1,000 trials. The historical `2/Ns` separation and OR
    rule remain. The even-array reference-element convention is retained.
17. **Independent random streams and completion status.** Each seed regenerates
    positions, gains and observations, with a separate W stream. Input generation
    consumes the full SNR sequence before estimation. Completed, missing and
    numerically failed jobs remain distinguishable. The release retains the
    historical one-computational-thread policy in each numerical call, including
    workers, and restores the caller's thread setting afterward.
18. **Twenty true updates for M8.** `fig11_trace.m` records actual updates at every
    iteration, with snapshot curvature computations not fed back into the
    estimator. The first ten are plotted. There is no running minimum or repeated
    stopped endpoint.
19. **Accuracy plus reliability.** A common 3 m trial-RMS threshold defines the
    conditional curves. Complete-sample exceedance, Wilson intervals, raw RMSE
    columns and the key full CDF remain available. Statistics do not modify
    estimator outputs. Wilson bounds are retained as table fields; no interval
    bars are drawn in the exceedance panels.

## Changes introduced by the open-source packaging

The release replaces historical parent-directory/snapshot dependencies and
per-trial artifacts with seed-based generation and a uniform `run_experiments`
entry point. It holds one active dictionary per process, batches independent
seeds, verifies compact current/previous checkpoints, rejects incompatible
resume identities, and derives CSV results through one statistics routine.
It shares the SNR trial between position and channel evaluation, stores only
scalar metrics, and provides English/Chinese plotting from either fresh runs or
bundled CSVs. This changes storage and orchestration, not the localization
objectives, iteration budgets, or mathematical treatment of bad estimates.

## Remaining limits

Gaussian fusion is performed in direction-cosine/range coordinates, and the
implemented coarse stage runs only once. BCRB uses no explicit prior-information
term. These are retained implementation conventions to be distinguished from
more general formulations in the article.

For M=5, K=3, RF=60, the free coarse block-gain model has 75 columns but only
60 compressed measurements. Its rank deficiency and the resulting numerical
sensitivity are retained. There is no added geometric regularization, SVD
replacement, distance upper bound, positive-definite projection, or fallback
fusion in this release. Final Gaussian fusion does not enforce all the
intermediate solver's distance-domain checks. Moving to fewer blocks is a
different configuration, not a repair of the fixed failing configuration.

The grid parameter measures initialization resolution; dictionaries can be
prepared offline for a fixed W. Coarse refinement still evaluates `W*a(p)`.
This package therefore makes no unsupported claim that fewer initialization
grid points imply proportionally lower total online runtime.
