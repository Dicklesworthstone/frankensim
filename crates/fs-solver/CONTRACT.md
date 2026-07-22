# fs-solver CONTRACT

## Purpose and layer

Layer: **currently L3 FLUX** because the historical p-MG and Stokes
fixtures import `fs-feec` (L3); the other dependencies are `fs-la`,
`fs-sparse`, and `fs-spectral` (L1) plus `fs-tilelang` and `fs-math`
(L0). The generic nonlinear, block-operator, verification/admission,
and flexible-Krylov modules deliberately contain no `fs-feec` types
and form the lower-layer solver spine. Moving the package metadata to
L1 requires extracting the FEEC-specific adapters into an L3 crate;
this contract does not claim that extraction is complete. The solver
stack (plan §8.9) binds matrix-free methods to resumability,
cancellation boundaries, deterministic execution, adjoint support,
and error transparency (residual histories and structured stall
diagnoses, never timeout mysteries).

## Public types and semantics

- `LinearOp` — matrix-free operator trait: `n`, `apply`, and
  `apply_transpose` (defaulting to `apply` — correct ONLY for
  symmetric operators). `CsrOp::symmetric/general` adapts assembled
  matrices (general materializes the transpose once,
  deterministically).
  **NO-CLAIM: the crate does not verify operator symmetry.** There is
  no symmetry probe anywhere in fs-solver, so nothing detects a client
  `LinearOp` that inherits the default `apply_transpose` on a
  nonsymmetric operator. The transposed-solve battery pins the crate's
  OWN adapters — `CsrOp::general` and one nonsymmetric fs-opdsl
  convection–diffusion fixture that does override the method — and it
  cannot cover operators it never sees. An implementor who inherits the
  default on a nonsymmetric operator gets `GmresState::run(.., transposed
  = true)` solving `A·λ = rhs` instead of `Aᵀ·λ = rhs`; it converges
  normally and reports a small residual, but that residual is measured
  against the operator supplied, not against its transpose, so it is no
  evidence that the adjoint system was solved. Callers who need the
  guarantee must run their own `⟨Av, w⟩ = ⟨v, Aᵀw⟩` probe.
- `RectLinearOp` and `BlockOperator<N>` (`BlockOperator2` and
  `BlockOperator3` aliases) — dimension-checked rectangular block
  algebra with fixed-order primal and transpose application.
  `SquareBlock` explicitly adapts a `LinearOp` without a blanket
  implementation that makes method names ambiguous; `ZeroBlock`
  represents structural zeros, but an aggregate zero-dimensional
  block partition is refused as `BlockError::Empty` rather than
  becoming a vacuous solver operator. `RealEquivalentComplexOp`
  maps `A + iB` to the real block form `[[A,-B],[B,A]]`.
  `BlockSchur2` is an injected block-LDU preconditioner whose caller
  supplies the two diagonal solves and chooses the Schur sign.
- `CgState` — resumable preconditioned CG (SPD; fs-sparse `Precond`
  slot). `MinresState` — resumable MINRES (symmetric indefinite;
  Paige–Saunders with explicit two-level Givens memory in the state).
  `PminresState` — resumable preconditioned MINRES for symmetric
  indefinite systems with an SPD split preconditioner; the Stokes
  block-preconditioner battery is its first consumer. `GmresState` —
  resumable restarted GMRES(m) (general operators;
  `transposed = true` runs Aᵀ through the SAME machinery and
  preconditioner slot). Resume granularity: CG/MINRES per ITERATION,
  P-MINRES per ITERATION, GMRES per RESTART CYCLE (the Arnoldi basis
  is deliberately not checkpointed mid-cycle) — split runs at those
  boundaries are BITWISE-equal to straight runs, `clone()` is the
  checkpoint.
- `FgmresState` — resumable restarted flexible GMRES. The
  preconditioner is selected by logical iteration, so varying
  preconditioners are explicit and checkpoint/replay at restart
  boundaries remains bitwise deterministic.
- `LinearSystemVerifier` → `VerifiedLinearSystem` →
  `admit_linear_solver` — injected verification grammar and local
  decision receipt. CG requires symmetric positive-definite evidence
  and compatible/trivial nullspace evidence; MINRES requires verified
  symmetric-indefinite evidence plus no or fixed-SPD preconditioning;
  plain GMRES rejects variable preconditioners; FGMRES admits them.
  Every method refuses unresolved nullspace or unverified source
  compatibility. Refusals are structured and never silently widened.
- `NonlinearProblem`, `NewtonKrylovConfig`, and `NewtonKrylovState` —
  resumable inexact Newton-Krylov with Eisenstat-Walker forcing,
  FGMRES inner solves, Armijo line search or trust-region
  actual/predicted acceptance, and per-iteration telemetry. Outer
  checkpoints are `clone()` boundaries; construction failures and
  stalls are typed in `NewtonError`/`NewtonStallDiagnosis`.
- `spectral_service` — re-export of the one workspace `fs-spectral`
  authority. This is an ownership seam, not a claim that the generic
  eigensolver service is already implemented.
- `SolveReport { iters, rel_residual, converged, history, diagnosis }`
  with `StallDiagnosis::{Plateau, BudgetExhausted, Breakdown}`:
  Plateau = < 5% relative progress over the last 50-iteration window
  (calibrated: CG's recursive residual on inconsistent singular
  systems DIVERGES rather than plateauing — see `diag_probe.rs`).
  **`rel_residual` is not one quantity, and `converged` is not a
  Euclidean claim.** The field carries whatever the producing solver
  measures, and `converged` means only "that measure met `tol`":
  GMRES(m) and FGMRES recompute the true `‖b − Ax‖₂/‖b‖₂`; CG reports
  the RECURSIVELY propagated residual (never recomputed as `b − Ax`,
  so it drifts, and it drifts low); MINRES reports the Givens `|η|/‖b‖₂`
  estimate; and P-MINRES reports `|η|/‖b‖_M` — an estimate in the
  PRECONDITIONER's M-norm, which for an ill-conditioned SPD `M` is
  unrelated in magnitude to the Euclidean residual (with
  `M = diag(1, 1e-12)`, `r = (0,1)` against `b = (1,0)` gives
  `‖r‖_M/‖b‖_M = 1e-6` while `‖r‖₂/‖b‖₂ = 1`). `ResidualClaim::
  {TrueEuclidean, RecursiveEstimate, MNormEstimate}` names the three,
  and every solver state exposes `residual_claim()`; a driver that
  needs a Euclidean statement from CG/MINRES/P-MINRES must recompute
  `b − Ax` itself. The provenance is NOT carried on `SolveReport`
  itself — that struct is constructed by struct literal outside this
  crate, so the field could not be added without breaking those
  callers.
- `PMultigrid` — matrix-free p-MG V-cycle as a `Precond`: order
  hierarchy r → r/2 → … → 1 over fs-feec `TensorSpace`s; prolongation
  is EXACT INJECTION (the hierarchical Lobatto basis nests, so the
  Galerkin coarse operator IS the coarse-order operator — nothing
  assembled except r = 1); the r = 1 coarse level is assembled
  (interior Kronecker CSR) and solved by SA-AMG-preconditioned CG at a
  REQUESTED 1e-13 within a 2000-iteration cap (a loosely-solved coarse
  level makes the V-cycle a VARYING preconditioner and demonstrably
  breaks plain CG). "Near-exact" is the request, not a proven property:
  the receipt is retained on the preconditioner and readable through
  `coarse_solves_converged()` / `worst_coarse_rel_residual()` /
  `coarse_solves()`, and `reset_coarse_evidence()` scopes it to one
  outer solve.
  SMOOTHING (bead x08j): Chebyshev (band [λmax/16, λmax], power-method
  λmax of the PRECONDITIONED operator) accelerating PU-symmetrized
  vertex-patch additive Schwarz (all dofs of the ≤ 8 elements around
  each interior vertex; exact patch inverses by FAST DIAGONALIZATION —
  the patch operator is exactly the Kronecker sum of 1D windows, so
  per-axis generalized eigenproblems of size ≤ 2r+1 replace any
  (2r+1)³ dense factorization, eigendata shared across the ≤ 3 trim
  signatures) PLUS the exact r = 1 coarse term (Pavarino combination).
  Requires m ≥ 2 (asserted). MEASURED design ledger: PU weighting is
  load-bearing (without it counts jump 8 → 13 when the per-axis window
  multiplicity first hits 3 at m = 4); eigendata may be shared across
  same-signature vertices but window OFFSETS may not (sharing the
  representative's indices left cells uncovered at m ≥ 5: 29 iters at
  m = 5, outright failure at m = 6 — caught by the beyond-acceptance
  spot-checks, fixed, gated).
- `MaskedTensorOp` — the homogeneous-Dirichlet high-order Poisson
  apply as a `LinearOp`.
- `PminresState` (bead avuw) — resumable PRECONDITIONED MINRES with an
  SPD preconditioner (Lanczos in the M-inner product); `rel_residual`
  estimates the M-norm residual (cross-checked against TRUE residuals
  in the battery); resume is bitwise like the other Krylov states.
- `stokes::{StokesSystem, StokesOp, StokesBlockDiag}` (bead avuw) —
  the tensor Stokes saddle fixture Q_r³/P_{r−1}^disc with the
  Silvester–Wathen block-diagonal preconditioner (p-MG per velocity
  component + DIAGONAL pressure-mass inverse); constant-pressure null
  handled by projection. MEASURED rejection on record: the full-tensor
  Q_{r−1}^disc pressure (the bead's literal reading) is not uniformly
  inf-sup stable — counts grew 44 → 101 → 137 across m = 2..4 at
  r = 2 (the Q1/Q0 checkerboard family); the total-degree subset
  flattens them to 26..61.
- `StokesSystem`, `StokesOp`, `StokesBlockDiag` — unit-cube FEEC
  tensor Stokes fixture over Q_r^3/P_{r-1}^disc with projected
  constant-pressure null mode, block operator [[A, B^T], [B, 0]], and
  SPD block-diagonal inverse preconditioner (p-MG per velocity
  component plus diagonal pressure-mass inverse).
- `mixed::{CsrF32, mixed_cg_refine, MixedReport}` — f32 INNER CG
  (storage AND arithmetic: half the memory traffic per iteration)
  under f64 iterative refinement with TRUE f64 residuals and scaled
  corrections; f64-grade accuracy whenever κ(A) ≪ 1/ε_f32; stalls set
  `escalate` in the report (drivers fall back to the f64 path —
  reported, never silently absorbed). Deterministic (bitwise repeat
  tested).
- `dot`/`norm2` — deterministic fixed-shape reductions (fs-tilelang
  combiner; shape depends on length only).

## Invariants

- Checkpoint = `clone()`; split runs bitwise-equal at the stated
  granularity (tested at multiple cut points per method). Newton
  checkpoints are completed outer attempts; FGMRES checkpoints are
  restart boundaries.
- All inner products flow through the fixed-shape reduction — no
  thread- or tier-dependent bit patterns anywhere.
- The V-cycle preconditioner is symmetric by construction (identical
  pre/post smoothing, Galerkin-consistent transfers). Its FIXITY —
  the property plain CG actually needs — additionally requires the
  r = 1 coarse solve to be near-exact, and that is REQUESTED (1e-13
  within a 2000-iteration cap), not established: an application that
  exits on the cap makes the V-cycle an inexact, application-dependent
  operator. `PMultigrid::coarse_solves_converged()` (with
  `worst_coarse_rel_residual()` and `coarse_solves()`) is the retained
  receipt; a driver that does not consult it is asserting fixity
  without evidence. No seeded fault currently makes the cap fire.
- Transposed solves share every piece of primal infrastructure.
- The G3 bounded-integer nonsymmetric fixture satisfies the exact linear
  adjoint/finite-difference identity through `CsrOp::apply` and
  `CsrOp::apply_transpose`. This is not a claim about arbitrary nonlinear
  objectives or floating-roundoff regimes.
- Block row/column partitions and their nonzero aggregate dimension
  are validated once, then traversed in stable index order for both
  primal and transpose applications.
- A solver admission decision can only be constructed from a coherent
  local verifier finding; the decision function cannot silently
  upgrade missing or contradictory evidence into a stronger solver
  class. Execution authorization remains a higher-layer concern.

## Error model

Structured panics on dimension mismatches (programmer errors).
Non-finite legacy Krylov iteration quantities surface as
`StallDiagnosis::Breakdown`; Newton construction uses `NewtonError`
and in-flight attempts use `NewtonStallDiagnosis::NonFinite`, never UB
or silent NaN propagation.
Un-converged solves return reports with history + diagnosis.
Block construction returns `BlockError` for inconsistent partitions.
Linear verification/admission returns structured findings naming the
failed evidence rule. Nonlinear solves return typed dimension,
non-finite, inner-solve, line-search, trust-region, and budget stalls.

## Determinism class

Bit-deterministic across runs and (single-threaded v1) trivially
across thread counts; cross-ISA golden FNV-64 over CG/MINRES/GMRES
solutions and a pMG-preconditioned solve:
`0xbc00_5985_1f9c_4a8a`, recorded on Apple M4 Pro, verified on
Threadripper (x86_64).
Block traversal, verifier findings, FGMRES logical preconditioner
selection, and Newton globalization decisions use fixed iteration
orders. The new battery checks bitwise split/replay for FGMRES and
Newton at their documented checkpoint boundaries.

## Cancellation behavior

Iteration-granular: every state is complete between `run` calls, so
drivers interrupt by not continuing — request → drain (finish the
current iteration) → finalize (state is the checkpoint). fs-exec Cx
wiring lands with the drivers (workspace discipline).
Newton drains its active inner solve/globalization decision before an
outer checkpoint; FGMRES drains the active restart cycle. Neither new
state yet accepts a `Cx` directly.

## Unsafe boundary

None. `unsafe_code = "deny"`.

## Feature flags

None.

## Conformance tests

`tests/solver_battery.rs` (11 cases): CG vs dense-LU reference on the
FEEC Poisson system + bitwise resume at 3 cut points + G5 repeat;
MINRES on a genuinely indefinite shifted Poisson vs LU, with the |η|
residual estimate cross-checked against the TRUE residual, + bitwise
resume; GMRES on nonsymmetric convection–diffusion (fs-opdsl fixture)
vs LU + TRANSPOSED solve verified (Aᵀy = c residual ≤ 1e−10,
convergence within 2× of primal — adjoint readiness) + cycle-boundary
bitwise resume; structured diagnoses (BudgetExhausted on a
short-budget solve; Plateau on the classic restarted-GMRES stagnation
fixture — a cyclic shift matrix where every cycle reproduces x = 0);
mixed precision: f32-inner refinement reaching 5e−15 true relative
residual, matching plain f64 CG to 1e−9 and bitwise-repeatable;
p-MG ladder gates: converges everywhere with iteration counts inside
a fixed envelope (≤ 80) while identity-preconditioned counts blow up
(hard corner m = 4, r = 4: 9 vs 1192 — a 132× advantage), solutions
matching identity-CG to 1e−7; deterministic dot; golden hash (bumped
at x08j with justification: the smoother change is semantic).
`tests/metamorphic.rs`: relation `general-csr-linear-adjoint-difference`
exercises both production actions of `fs-solver::CsrOp::general` through the
shared, jointly shrinkable adjoint/finite-difference harness (seed
`0x2ACE_0004`, 384 cases, exact finite-bit tolerance). The fixture is one fixed
nonsymmetric 3x3 integer-valued CSR matrix with objective `[1, -2, 3]`; base
coordinates and a nonzero direction lie in `[-8, 8]`, and the jointly
shrinkable integer step lies in `[-3, 3]`. It checks that the primal finite
difference equals the `apply_transpose` directional derivative and that the
adjoint action is stable across the transform. This code-first declaration is
proof-pending until the owning-crate batch suite passes. The broader GMRES
transpose/LU pin above remains independent.
`tests/stokes_battery.rs` (5 cases): FEEC tensor Stokes fixture with
P-MINRES + blockdiag(p-MG, pressure mass) agrees with a dense
constant-pressure-pinned LU reference on m = 2, r = 2; velocity is
divergence-free to solver tolerance; mesh/order ladder iteration
envelope maxes at 61 iterations in the current battery; P-MINRES
resume is bitwise at cuts 1/7/23; Stokes golden hash
`0x5754_3908_cb41_7281`. The dense reference pins an actual
cell-constant pressure coefficient; pinning an arbitrary high-order
pressure mode leaves the nullspace intact and was rejected by the
fresh-eyes audit.
`tests/ladder_probe.rs` (bead x08j acceptance): FIXED smoothing degree
3 across both ladders — order ladder m = 3, r = 2..6 iters
[8, 8, 8, 9, 11] (max/min ≤ 1.5), mesh ladder r = 3, m = 2..5 iters
[2, 8, 8, 9] (max/min ≤ 1.5 over the nontrivial m ≥ 3; the m = 2
single-patch case is an EXACT solve and is gated as not-slower rather
than rewarding the trivial minimum), plus the m = 6 window-sharing
spot-check. `tests/diag_probe.rs`: the diagnosis-calibration
regression (singular-system CG diverges — must never read Plateau).
`tests/nonlinear_block_battery.rs` (G0/G1/G3/G5): 2x2 and 3x3 block
primal/transpose equivalence, real-equivalent complex algebra, an
exact manufactured Schur saddle solve, verifier/admission refusals,
non-collinear iteration-varying diagonal FGMRES preconditioners that
exercise stored `z_j` directions plus split replay, and Newton solves with
line-search/trust-region globalization, an explicit rejected-search
receipt, and outer split replay.

## No-claim boundaries

- p-growth of iteration counts: RESOLVED (bead x08j) — the
  vertex-patch Schwarz smoother holds counts flat across r = 2..6 and
  m = 3..6 (gated in `ladder_probe.rs`). Historical dead ends stay on
  record: element-matrix (EBE) and single-element Dirichlet-block
  Schwarz were measured WEAKER than Jacobi-Chebyshev (50 and 41 iters
  vs 8); the working design needed vertex-centered tensor windows +
  fast diagonalization + PU symmetrization + the Pavarino coarse term.
  Anisotropic/stretched meshes and non-tensor patches remain
  out-of-scope (the fast-diagonalization patch inverse requires the
  Kronecker structure).
- FGMRES now supports varying preconditioners at restart-boundary
  checkpoint granularity. It does not serialize an in-flight Arnoldi
  basis. Plain CG with a varying preconditioner remains known-broken
  (observed), and plain GMRES admission refuses that configuration.
- Saddle-point block preconditioners are implemented only for the
  unit-cube tensor Stokes fixture over Q_r^3/P_{r-1}^disc with
  pressure-mass Schur approximation. General fs-feec simplicial
  vector families, non-tensor mixed spaces, Navier-Stokes convection,
  pressure-robust production discretizations, and broader Schur
  approximations remain follow-up scope. The mixed-precision path
  measures its bandwidth SPEEDUP on the perf lane (fz2.4
  release-profile), not here; no ladder policy engine yet (fs-la's
  `Ladder` decides for factorizations; a Krylov-side policy joins when
  consumers exist).
- No dense output of Krylov bases, no eigenvalue estimation service,
  no threading (fs-exec drivers own parallelism), no G4 cancellation
  storms yet (needs the Cx wiring).
- The G3 CSR adopter proves only an exact relation on one fixed 3x3 matrix and
  bounded integer-valued vectors/steps whose binary64 arithmetic is exact. It
  does not establish adjoint accuracy under general floating-point roundoff,
  arbitrary sparse operators or dimensions, nonlinear finite differences,
  Krylov convergence, or derivative correctness for solver clients. Its
  declaration is not evidence of a passing batch run.
- The verifier is an injected evidence grammar, not an automatic
  symmetry/definiteness/nullspace certifier. Its receipt retains only
  dimension, verifier identity, and findings: it does not content-bind
  an opaque operator/RHS, and raw Krylov states remain callable. A
  ledgered caller must bind and consume the decision if it needs an
  authorization capability. Newton convergence is empirical residual
  evidence, not an interval enclosure or a proof of uniqueness. The
  generic `fs-spectral` eigensolver service remains pending even though
  its authority namespace is re-exported here.
- FGMRES checkpoints assume the exact same operator, RHS, and logical
  preconditioner policy on resume; Newton checkpoints assume the exact
  same problem implementation and parameters. These opaque input
  identities are caller/ledger invariants and are not authenticated by
  the lower-layer state itself.
- This package remains labeled L3 until its FEEC-specific p-MG/Stokes
  adapters are extracted; the dependency-light modules alone do not
  establish a completed architectural layer move.
- The p-MG is specialized to the unit-cube tensor Poisson fixtures;
  general-operator p-MG arrives with the fs-opdsl matrix-free atom
  integration.
