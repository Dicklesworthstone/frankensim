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
  `GmresState` — resumable restarted GMRES(m) (general operators;
  `transposed = true` runs Aᵀ through the SAME machinery and
  preconditioner slot). Resume granularity: CG/MINRES per ITERATION,
  GMRES per RESTART CYCLE (the Arnoldi basis is deliberately not
  checkpointed mid-cycle) — split runs at those boundaries are
  BITWISE-equal to straight runs, `clone()` is the checkpoint.
- `SolveReport { iters, rel_residual, converged, history, diagnosis }`
  with `StallDiagnosis::{Plateau, BudgetExhausted, Breakdown}`:
  Plateau = < 5% relative progress over the last 50-iteration window
  (calibrated: CG's recursive residual on inconsistent singular
  systems DIVERGES rather than plateauing — see `diag_probe.rs`).
- `PMultigrid` — matrix-free p-MG V-cycle as a `Precond`: order
  hierarchy r → r/2 → … → 1 over fs-feec `TensorSpace`s; prolongation
  is EXACT INJECTION (the hierarchical Lobatto basis nests, so the
  Galerkin coarse operator IS the coarse-order operator — nothing
  assembled except r = 1); matrix-free Chebyshev smoothing on the
  Jacobi-scaled operator (band [λmax/16, λmax], hardened
  fixed-iteration power method for λmax); the r = 1 coarse level is
  assembled (interior Kronecker CSR) and solved near-exactly by
  SA-AMG-preconditioned CG (a loosely-solved coarse level makes the
  V-cycle a VARYING preconditioner and demonstrably breaks plain CG).
- `MaskedTensorOp` — the homogeneous-Dirichlet high-order Poisson
  apply as a `LinearOp`.
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
(≥ 10× advantage at the hard corner: measured 1192 vs ≈70 at
m = 4, r = 4), solutions matching identity-CG to 1e−7; deterministic
dot; golden hash. `tests/diag_probe.rs`: the diagnosis-calibration
regression (singular-system CG diverges — must never read Plateau).

## No-claim boundaries

- p-MG iteration counts grow MILDLY with r (hierarchical-injection
  CBS angle; measured, envelope-gated, not hidden). True
  p-independence needs an overlapping-patch Schwarz smoother —
  recorded follow-up; element-matrix (EBE) and single-element
  Dirichlet-block variants were implemented and measured WEAKER than
  Chebyshev here, so the naive versions are a dead end on this
  evidence.
- Preconditioned MINRES (needs an SPD split preconditioner) and
  flexible-GMRES (varying preconditioners) are follow-up scope; plain
  CG with a varying preconditioner is known-broken (observed) — use
  GMRES or fix the preconditioner.
- Saddle-point block preconditioners (Stokes-class, Schur
  approximations over the fs-feec mixed spaces) are split scope (the
  vector-MMS machinery they need is itself gated on the simplicial
  vector families). The mixed-precision path measures its bandwidth
  SPEEDUP on the perf lane (fz2.4 release-profile), not here; no
  ladder policy engine yet (fs-la's `Ladder` decides for
  factorizations; a Krylov-side policy joins when consumers exist).
- No dense output of Krylov bases, no eigenvalue estimation service,
  no threading (fs-exec drivers own parallelism), no G4 cancellation
  storms yet (needs the Cx wiring).
- The p-MG is specialized to the unit-cube tensor Poisson fixtures;
  general-operator p-MG arrives with the fs-opdsl matrix-free atom
  integration.
