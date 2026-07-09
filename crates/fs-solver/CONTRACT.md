# fs-solver CONTRACT

## Purpose and layer

Layer: **L3 FLUX** (deps: fs-feec L3, fs-la/fs-sparse L1,
fs-tilelang/fs-math L0). The solver stack (plan §8.9): matrix-free
Krylov methods and p-multigrid bound by the four workspace contract
obligations from day one — resumable, cancellable, deterministic,
adjoint-equipped — with error transparency (residual histories and
structured stall diagnoses, never timeout mysteries).

## Public types and semantics

- `LinearOp` — matrix-free operator trait: `n`, `apply`, and
  `apply_transpose` (defaulting to `apply` — correct ONLY for
  symmetric operators; the transposed-solve battery catches
  violations). `CsrOp::symmetric/general` adapts assembled matrices
  (general materializes the transpose once, deterministically).
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
- `SolveReport { iters, rel_residual, converged, history, diagnosis }`
  with `StallDiagnosis::{Plateau, BudgetExhausted, Breakdown}`:
  Plateau = < 5% relative progress over the last 50-iteration window
  (calibrated: CG's recursive residual on inconsistent singular
  systems DIVERGES rather than plateauing — see `diag_probe.rs`).
- `PMultigrid` — matrix-free p-MG V-cycle as a `Precond`: order
  hierarchy r → r/2 → … → 1 over fs-feec `TensorSpace`s; prolongation
  is EXACT INJECTION (the hierarchical Lobatto basis nests, so the
  Galerkin coarse operator IS the coarse-order operator — nothing
  assembled except r = 1); the r = 1 coarse level is assembled
  (interior Kronecker CSR) and solved near-exactly by
  SA-AMG-preconditioned CG (a loosely-solved coarse level makes the
  V-cycle a VARYING preconditioner and demonstrably breaks plain CG).
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
  granularity (tested at multiple cut points per method).
- All inner products flow through the fixed-shape reduction — no
  thread- or tier-dependent bit patterns anywhere.
- The V-cycle preconditioner is symmetric (identical pre/post
  smoothing, Galerkin-consistent transfers, near-exact coarse).
- Transposed solves share every piece of primal infrastructure.

## Error model

Structured panics on dimension mismatches (programmer errors).
Non-finite iteration quantities surface as
`StallDiagnosis::Breakdown`, never UB or silent NaN propagation.
Un-converged solves return reports with history + diagnosis.

## Determinism class

Bit-deterministic across runs and (single-threaded v1) trivially
across thread counts; cross-ISA golden FNV-64 over CG/MINRES/GMRES
solutions and a pMG-preconditioned solve:
`0xbc00_5985_1f9c_4a8a`, recorded on Apple M4 Pro, verified on
Threadripper (x86_64).

## Cancellation behavior

Iteration-granular: every state is complete between `run` calls, so
drivers interrupt by not continuing — request → drain (finish the
current iteration) → finalize (state is the checkpoint). fs-exec Cx
wiring lands with the drivers (workspace discipline).

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
- Flexible-GMRES (varying preconditioners) is follow-up scope; plain
  CG with a varying preconditioner is known-broken (observed) — use
  GMRES or fix the preconditioner.
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
- The p-MG is specialized to the unit-cube tensor Poisson fixtures;
  general-operator p-MG arrives with the fs-opdsl matrix-free atom
  integration.
