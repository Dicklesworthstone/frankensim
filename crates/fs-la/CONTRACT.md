# CONTRACT: fs-la

> Status: PARTIAL — the GEMM, FACTORIZATION, MIXED-PRECISION, and
> EIGENSOLVER sections below are in force; batched small-dense remains
> skeleton scope.

## Purpose and layer
Dense linear algebra: GEMM, batched small dense, factorizations, eigensolvers. Layer: L1.

## Public types and semantics
- `gemm::{gemm_f64, gemm_f32, gemm_mixed}` — C = α·A·B + β·C on row-major
  contiguous slices (BLAS-shape signatures). β = 0 OVERWRITES C (NaN and
  garbage in C are ignored — the uninitialized-output convention).
  `gemm_mixed` is f32 STORAGE with f64 ACCUMULATION (exact widening; the
  bandwidth-vs-accuracy mode the plan uses throughout).
- f64 path: BLIS-style NC→KC→MC blocking with A/B panel packing and an
  MR×NR register-tiled microkernel (safe Rust, fused mul_add). f32/mixed
  paths share the loop order and KC chunking, unpacked in v1.
- `factor::{cholesky, lu, qr, tsqr_r, svd_jacobi}` + `FactorError` —
  dense factorizations. Failure is DATA: `NotSpd{index}` /
  `Singular{index}` typed diagnostics, never panics for data conditions.
  `Cholesky::solve`, `Lu::{solve, solve_transpose, condition_1}` (Hager
  1-norm estimate), `Qr::{apply_q, apply_qt, solve_ls}` are the
  refinement/consumer hooks. LU pivot tie-break: LOWEST index under equal
  magnitude (P2). `tsqr_r` computes the sign-canonicalized (non-negative
  diagonal) R over a binary combine tree whose shape is a pure function
  of (m, row_block, n). `svd_jacobi` is one-sided cyclic Jacobi (thin
  U·Σ·Vᵀ, σ descending, deterministic order and tie-breaks).
- `mixed::{solve_adaptive, ResidualTarget, Ladder, RefineReport}` —
  precision-ladder solves with iterative refinement AS POLICY:
  f32-factor/f64-refine → f64-direct → f64-factor/dd-refine. `backward`
  targets the normwise backward error; `forward` (optional) targets the
  relative forward error — the distinction is load-bearing: any stable
  f64 solve is backward-accurate to ~eps, but beating κ·eps FORWARD
  requires the dd (extended-residual) rung, and that is exactly when it
  engages. Every solve returns a `RefineReport` (rung, steps, achieved,
  escalations, condition estimate, correction-ratio forward estimate,
  full residual trajectory) — a precision decision is EVIDENCE, never a
  silent downgrade. Run-dry is reported honestly (converged = false,
  best-achieved recorded), never panicked.
- `eigen::{jacobi_eigh, lanczos_run/LanczosState,
  lobpcg_run/LobpcgState, EigenPair}` — symmetric eigensolvers.
  `jacobi_eigh`: dense cyclic Jacobi (ascending eigenvalues, orthonormal
  columns, lowest-index tie-break). Lanczos (full reorthogonalization)
  and LOBPCG (X−X_prev conjugate direction, optional preconditioner) are
  MATRIX-FREE (operator = `Fn(&[f64], &mut [f64])` closure — fs-sparse
  SpMV or stencils plug in without format coupling) and RESUMABLE
  (checkpoint = state clone; split runs are BITWISE equal to straight
  runs, tested). Every `EigenPair` carries the TRUE operator residual
  ‖A·v − λ·v‖₂ (recomputed, not the internal estimate). Start vectors
  use fs-math STRICT sin — platform libm here caused a real cross-ISA
  golden-hash divergence, caught by the trj pipeline and fixed; this is
  now a contract rule: NO platform libm in any path that feeds solver
  state.
- `eigen_complex::{eig, det_complex, EigFailure}` — general complex
  nonsymmetric eigenvalues: Givens Hessenberg reduction + explicitly
  shifted QR (Wilkinson shifts with total_cmp tie-breaks, exceptional
  shifts every 12 iterations, standard deflation, closed-form 2×2
  blocks). Canonical (re, im)-sorted output. Convergence exhaustion is
  a typed `EigFailure`, never a wrong answer. Tested: companion roots
  of unity, rotation conjugate pairs, agreement with `jacobi_eigh` on
  embedded symmetric matrices (1e-10), trace/determinant identities
  against an independent Gaussian-elimination `det_complex` oracle,
  Hermitian spectral reality, T6 colleague-matrix roots.
- `VERSION` for provenance stamping.

## Invariants
1. Packing is ARITHMETIC-NEUTRAL: the packed/blocked f64 path is bitwise
   identical to a same-order naive loop (tested across a 9-shape sweep
   including k=0, m=1, tails in every dimension, tall-skinny, wide).
2. Row (m) and column (n) tiling are bit-neutral; KC chunking is PART OF
   THE BIT CONTRACT (per-chunk register partials fold into C in chunk
   order — changing KC legitimately changes bits and requires a golden
   bump with justification). Submatrix consistency tested.
3. `gemm_mixed` output is bitwise equal to the f64 computation on
   exactly-widened inputs (tested).
4. (A·B)ᵀ = Bᵀ·Aᵀ within 1e-13 relative (order differs; not bitwise).
5. Factorization residuals (tested): ‖A−LLᵀ‖/‖A‖ ≤ n·1e-14 on SPD;
   LU solve round-trips at 1e-9 on random; A = QR reconstruction at
   1e-12 with Q orthogonal to 1e-13; TSQR R equals direct QR's
   canonicalized R (1e-10) for ANY tree shape and satisfies the Gram
   identity; SVD reconstructs to 1e-13 with U, V orthogonal to 1e-13
   (Hilbert-8 spectral condition lands in the known ~1.5e10 band).
6. Factorizations are bit-deterministic given the blocking constants
   (fixed loop orders; GEMM's KC contract inherited; TSQR tree fixed).
7. The mixed-precision LADDER DECISION is deterministic: fixed thresholds
   (κ·eps32·16 < 1 admits the f32 rung; κ·eps64·16 < forward-target
   admits working-precision rungs), fixed stall rule (two consecutive
   steps without halving), deterministic condition estimator — same
   input, same ladder, same trajectory bits (tested). f32 singularity
   escalates automatically (tested with an f32-collapsing matrix).
8. dd-refinement demonstrably converges to the exact solution of the
   STORED problem: at κ = 1e10 it beats the direct solve's forward error
   by ≥100× against a past-convergence reference (tested). Note the
   no-claim: it cannot recover accuracy already lost when b was rounded
   to f64 — ground truth is A⁻¹·fl(b), not the user's pre-rounding
   intent.

## Error model
Slice-length/shape mismatches panic with structured messages (programmer
errors). DATA conditions in factorizations return `FactorError` with the
offending index: non-SPD pivots, exactly-singular columns. LU `growth`
exposes the pivot-growth statistic for ledgering.

## Determinism class
GEMM: bit-deterministic CROSS-ISA by construction (fixed loop order,
fused mul_add, no threading in v1). Evidence: FNV-64 golden hash over a
48×36×300 α-scaled product = `0x1d7a_a3c6_b631_7ef0`, recorded on
aarch64-apple, required to match on x86-64 in the test suite.
Factorizations: same class; golden hash over Cholesky L + LU solve +
TSQR R + SVD σ = `0x181f_8f95_82d6_87ed`, verified identical on both
reference ISAs. Mixed precision: ladder decisions + solutions + reports
hashed over a κ ∈ {1e3, 1e7, 1e11} battery = `0x8e09_2d4a_ff1b_5028`
(bumped once: test fixture hardened from platform `powf` to strict
exp/ln). Eigensolvers: Lanczos + LOBPCG outputs on the 1D Laplacian =
`0x87da_0cb3_2344_b097` (bumped once: start vectors moved from platform
sin to strict sin after an ACTUAL cross-ISA divergence — the golden
discipline catching its first real bug). All verified identical on both
reference ISAs.

## Cancellation behavior
All future hot paths poll cancellation at tile boundaries (Decalogue P7).
No compute paths exist yet.

## Unsafe boundary
None. `unsafe_code` is denied workspace-wide; any future capsule must be
registered per docs/CONVENTIONS.md and ship a SAFETY.md.

## Feature flags
None. Frontier features use `frontier-*`, moonshots `moonshot-*`, default off.

## Conformance tests
In-crate GEMM suite: bitwise same-order oracle across shape sweep, β/α
edge semantics (β=0 NaN overwrite, α=0, k=0, empty m/n), transpose
identity, submatrix consistency, mixed == widened-f64 bitwise, f32
tolerance battery, determinism + golden hash. tests/conformance.rs
placeholder remains for the shared-harness migration.

## No-claim boundaries
- **No performance claims yet**: v1 microkernel is safe auto-vectorized
  Rust with fixed pre-autotuner blocking. The ≥75%-of-peak roofline
  target, arch-specific fs-simd capsule microkernels, autotuned blocking,
  CCD-aware fs-exec parallel tiling, and f32/mixed packing belong to the
  recorded perf follow-up bead (gated on the autotuner).
- No transposed-operand or strided (non-contiguous) input forms yet.
- Factorization v1 is single-threaded; fs-exec tile-parallel
  panel/update drivers and arena packing are recorded follow-up scope.
  The compact-WY trailing update is applied reflector-sequentially (the
  fused WY GEMM form joins the perf lane).
- `tsqr_r` returns R only; implicit-Q tree factors (for applying Qᵀ in
  parallel TSQR) join the fs-exec driver work.
- `condition_1` is an estimate (typically within a small factor; a lower
  bound in theory) — not a certified bound (fs-ivl owns certified claims).
- Jacobi SVD targets small/medium n (O(n²·m) per sweep); no blocked
  driver yet. Batched small-dense remains skeleton scope.
- `eigen_complex::eig` runs UNBALANCED explicit shifted QR (O(n²) per
  sweep; implicit bulge-chasing and a balancing pass are recorded perf/
  robustness refinements) and returns eigenvalues only (eigenvectors
  via inverse iteration are follow-up). Canonical output ordering is
  roundoff-sensitive for near-tied real parts — consumers comparing
  spectra should match as SETS (the battery does).
- Lanczos v1 uses FULL reorthogonalization (selective ω-recurrence is a
  recorded refinement); LOBPCG has no deflation/soft-locking yet and
  identity preconditioning by default. Eigenvector adjoints are recorded
  follow-up (dλ/dp = vᵀ(∂A/∂p)v composes caller-side with fs-ad).
- `condition_estimate` in `RefineReport` is a Hager-style ESTIMATE (not
  a certified bound); the ladder thresholds are engineering headroom,
  not proofs. Componentwise (per-entry) backward targets, Krylov
  inner/outer precision splits, and ledger event emission are recorded
  follow-up scope (fs-obs schema wiring + consumers do not exist yet).
