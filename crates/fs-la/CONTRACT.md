# CONTRACT: fs-la

> Status: PARTIAL — the GEMM, FACTORIZATION, MIXED-PRECISION,
> EIGENSOLVER, and BATCHED SMALL-DENSE sections below are in force.

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
- `gemm_f64_parallel_with_pool(..., pool, gate)` is the production
  cancellation-aware MC/NC engine. It retains shared-B NC/KC orchestration and
  dispatches disjoint M bands as `TileKernel` work through the caller-owned,
  reusable `TilePool`; every band receives the pool's real `Cx` (gate, budget,
  stream identity, scoped arena, and mode). It returns `GemmRunReport` with one
  fs-exec `RunReport` per NC/KC panel on a full transactional commit, or
  `GemmRunError` after draining all scoped workers; an error leaves
  caller-visible C bitwise unchanged. `gemm_f64_parallel_with_cancel` is the
  convenience wrapper that constructs an unpinned host pool and preserves the
  same structured cancellation/executor error surface.
  `gemm_f64_parallel_with_pool_budgeted` additionally requires a
  `GemmMemoryEnvelope`. Its checked plan covers transactional C staging, shared
  B packing, the reusable M-band table, fs-la's ordered panel-receipt vector,
  and one exact `ArenaPool::reservation_bytes_for_slice` A-panel reservation
  per active `min(pool.workers(), m_bands)` worker. No-product calls plan only
  C staging. Root capacities are fallibly pre-reserved before dispatch; every
  tile receives a finite cost quota equal to one arena reservation; typed arena
  allocation refusal drains siblings and returns `MemoryRefused` rather than a
  panic. `GemmMemoryReport` records the requested plan, conservative logical
  reservation high-water (arena attempts count when entered), and refused
  component bytes. Arithmetic uses checked u128 accounting; unrepresentable
  layouts or totals fail closed. Legacy wrappers explicitly use the unbounded
  envelope while retaining the same finite per-tile arena quota.
  `gemm_tuning_is_effective` is the producer-owned routing query used to
  avoid publishing tune evidence for single-thread, small-M, and no-op
  calls. `gemm_execution_tier`, `GEMM_IMPLEMENTATION_VERSION`, the executing
  `TilePool::placement_identity()`, and `GEMM_BUILD_FINGERPRINT` are tune/replay
  identity material. The execution tier is operation-specific: AVX-512 host
  capability maps to the current AVX2/FMA GEMM capsule only when both AVX2 and
  FMA are available; otherwise the scalar fallback is reported. The build fingerprint
  is BLAKE3 over compiler version plus the watched/resolved compiler and active
  wrapper executable bytes, Cargo profile/codegen inputs, target, explicit Rust
  flags, workspace manifests, and a sorted content snapshot of
  the fs-la/fs-simd/fs-exec/fs-alloc/fs-substrate/fs-blake3/fs-obs Rust source
  closure, the constellation lock, and the actual asupersync, proc-macro, and
  Franken evidence/decision/kernel source closure, its unconditional
  non-source compiler inputs, and its Git identity,
  plus an optional explicit `FRANKENSIM_GEMM_CODEGEN_ID` salt; rows from another
  generated-code identity cannot be adopted. Build systems that inject codegen
  settings or code generators outside that bounded closure must supply the
  salt; a missing required environment, source, or workspace identity input is
  a build error rather than a generic fallback. Relevant dirty asupersync
  source is content-addressed directly, so its identity remains exact instead
  of pretending that Git HEAD describes the working tree. This is not a claim
  to hash arbitrary external tools or other path-dependency sources outside
  the named closure.
- Dependency-graph evidence (bead fz2.6, codegen schema v2): source bytes do
  not pin dependency CODEGEN — Cargo feature unification can compile the
  normal-dependency closure (asupersync and its registry closure) under a
  different unified feature set without moving any other fingerprint input.
  Every build therefore binds exactly one evidence class into the
  fingerprint, exposed as `GEMM_GRAPH_EVIDENCE`: `receipt:<blake3>` when
  tooling exported `FRANKENSIM_DEPGRAPH_RECEIPT` (canonical resolved
  package+feature receipt from `cargo run -p xtask -- depgraph-receipt --
  <same selection as the build>`, verifiable post-build with `--verify`;
  ambiguous resolutions refuse a receipt), or `salt:<value>` from the
  workspace `.cargo/config.toml` `FRANKENSIM_DEPGRAPH_SALT` — an explicit,
  operator-visible declaration that all such builds form ONE equivalence
  class. Neither present → the build fails closed. Citable perf lanes must
  use receipts (see `perf-baselines/README.md`); durable tune rows can only
  cross binaries inside one evidence class because the class payload is part
  of the fingerprint. Receipt precision is `cargo tree` resolution
  (dev-dependencies participate in unification, a conservative
  over-approximation of plain builds; the invocation-exact `-Z unit-graph`
  no longer exists in Cargo).
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
- `rand_nla::{range_finder, rsvd, nystrom_psd, sketch_ls, hutchinson,
  hutch_pp, RangeReport, TraceReport}` — randomized NLA [F].
  REPLAYABILITY BY CONSTRUCTION: all randomness flows from keyed Philox
  streams (stable kernel-id registry in-module), so every "stochastic"
  estimate is a pure function of its seed — bitwise reproducible and
  cross-ISA (fs-rand's distributions are themselves bit-deterministic).
  Results carry evidence (probe counts, posterior error estimates,
  variance estimates — the Evidence<T> integration point). Tested:
  rank-r errors within budgeted tails on fast/slow/gap spectra with
  estimate coverage; RSVD leading singular values to 1e-6; Nyström PSD
  reconstruction to 1e-6; sketch-and-precondition LS matches direct QR
  to 1e-8 with fast CG convergence; Hutch++ MSE decisively below
  Hutchinson at matched probe budgets (measured, not cited).
- `VERSION` for provenance stamping.

- `batched::{BatchMat, BatchVec}` — batches of n dense k×k matrices /
  k-vectors in entry-plane SoA layout (one 128-byte-aligned plane per
  entry across the batch; `fs_soa::FieldBuf` storage, plane stride
  padded to 16 f64). Gather/scatter to row-major AoS; plane accessors
  are the SIMD surface (lanes run across ELEMENTS, never within one
  small matrix).
- `batched::batch_gemm` — per matrix C = α·A·B + β·C (β = 0 overwrites,
  the house convention); α = 0 scales/overwrites C without reading A or B,
  matching the core GEMM contract; size classes {4, 6, 8, 12, 16, 24, 32, 48}
  dispatch to monomorphized kernels, other k run the same code
  generically.
- The batched-f64 bit surface is `fs-la:batched-f64-bits=1`; its
  `0x0377_a8c9_5992_aee9` golden is registered in
  `golden-couplings.json` and covers reduction order and batch-membership
  invariance. Alpha-zero no-read behavior is locked by a poisoned-operand
  regression because the frozen nonzero-alpha fixture intentionally does
  not move.
- `batched::{batch_det, batch_inv}` — closed forms for k ≤ 4 (Jacobian
  hot path); exactly-zero determinants flagged `Singular`, batch
  continues (flagged outputs unspecified).
- `batched::{batch_cholesky, batch_solve_lower, batch_solve_upper,
  batch_cholesky_solve}` — non-positive pivots flagged
  `NotSpd { index: diagonal step }` per member with the pivot replaced
  by 1.0 (finite continuation, no NaN storm); flag list authoritative.
- `batched::{batch_lu, BatchLu}` — partial pivoting per matrix
  (strictly-greater comparison = lowest row among maximal |pivot|,
  `factor::lu`'s tie-break); compact LU + per-step permutations +
  `Singular` flags; `BatchLu::solve` applies P, L, U in place.
- `batched::{batch_eigh3_values, batch_eigh}` — symmetric 3×3
  eigenvalues in closed form (deterministic trigonometric method
  through strict kernels, ascending); general small k (the 6×6 path)
  per-matrix through `eigen::jacobi_eigh` (values + orthogonal
  vectors).
- Per-element failures reuse `factor::FactorError` as
  `Vec<(matrix_index, FactorError)>` — empty means every member
  succeeded.

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
4. Cancellation is transactional: a cancelled parallel GEMM may have
   completed private microtiles, but C remains bitwise unchanged after the
   request is observed and every scoped worker drains. A successful call
   crosses one documented non-cancellable final commit boundary.
5. (A·B)ᵀ = Bᵀ·Aᵀ within 1e-13 relative (order differs; not bitwise).
6. Factorization residuals (tested): ‖A−LLᵀ‖/‖A‖ ≤ n·1e-14 on SPD;
   LU solve round-trips at 1e-9 on random; A = QR reconstruction at
   1e-12 with Q orthogonal to 1e-13; TSQR R equals direct QR's
   canonicalized R (1e-10) for ANY tree shape and satisfies the Gram
   identity; SVD reconstructs to 1e-13 with U, V orthogonal to 1e-13
   (Hilbert-8 spectral condition lands in the known ~1.5e10 band).
7. Factorizations are bit-deterministic given the blocking constants
   (fixed loop orders; GEMM's KC contract inherited; TSQR tree fixed).
8. The mixed-precision LADDER DECISION is deterministic: fixed thresholds
   (κ·eps32·16 < 1 admits the f32 rung; κ·eps64·16 < forward-target
   admits working-precision rungs), fixed stall rule (two consecutive
   steps without halving), deterministic condition estimator — same
   input, same ladder, same trajectory bits (tested). f32 singularity
   escalates automatically (tested with an f32-collapsing matrix).
9. dd-refinement demonstrably converges to the exact solution of the
   STORED problem: at κ = 1e10 it beats the direct solve's forward error
   by ≥100× against a past-convergence reference (tested). Note the
   no-claim: it cannot recover accuracy already lost when b was rounded
   to f64 — ground truth is A⁻¹·fl(b), not the user's pre-rounding
   intent.

## Error model
Slice-length/shape mismatches panic with structured messages (programmer
errors). DATA conditions in factorizations return `FactorError` with the
offending index: non-SPD pivots, exactly-singular columns. LU `growth`
exposes the pivot-growth statistic for ledgering. Pool GEMM returns
`GemmRunError::MemoryRefused` when its checked plan exceeds the caller envelope
or a fallible reservation is declined, retaining drained progress and memory
accounting; `MemoryPlanOverflow` refuses before allocation when the plan cannot
be represented. Cancellation and executor failures likewise retain the full
`GemmRunReport`, and caller-visible `C` is unchanged on every error path. The
failure-only reports are boxed so the hot `Result` representation remains
small; successful runs still return their report directly.

## Determinism class
GEMM: bit-deterministic CROSS-ISA by construction (fixed per-element loop
order and fused mul_add; MC/NC work assignment is disjoint and bit-neutral
across thread counts). Evidence: FNV-64 golden hash over a
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
reference ISAs. Batched small-dense: bit-deterministic by CONSTRUCTION
and additionally BATCH-MEMBERSHIP INVARIANT — a matrix factored in a
batch of N is bitwise-equal to the same matrix alone in a batch of 1
(each member's arithmetic is independent, all loops fixed-order
mul_add; tested per op family). Golden hash over GEMM + Cholesky solve
+ LU + det/inv + eigh3 = `0x0377_a8c9_5992_aee9`, verified identical on
both reference ISAs.
Randomized NLA: keyed-Philox replayable (every estimate is a pure
function of its seed); FNV-64 golden over RSVD σ + posterior estimate +
Hutchinson trace = `0xeef1_0550_7daf_c0d5`, verified identical on both
reference ISAs (aarch64-apple M4, x86-64 trj) and in BOTH debug and
release on each (trj:/data/tmp/rn_verify2/run_{release,debug}.log,
2026-07-09). Bumped once with cause: the original fixture
built its spectrum with `f64::powi`, whose rounding is
optimization-level-dependent (1-ulp debug/release divergence from
exponent 4 up), so the sentinel bits depended on build mode
(0x3e92_8bac_8cf9_fd48 debug vs 0xf3dc_b63b_e63f_8ab9 release — the
latter was briefly and incorrectly pinned). The fixture now uses a
fixed-order product chain. Contract rule, sibling of the platform-libm
rule above: NO `f64::powi` with a variable or >3 exponent in any path
that feeds golden bits (tracked workspace-wide in bead
frankensim-powi-build-mode-determinism-4xnt).

## Cancellation behavior
`gemm_f64_parallel_with_pool` implements request → drain → finalize under an
`fs_exec::CancelGate`. It polls while initializing/staging output and
pack storage (every 4096 elements), before each packed A/B micro-panel (at
most MR×KC / NR×KC copies), and before every MR×NR×KC compute tile (at most
8224 fused multiply-adds including alpha/write-back). A request stops
TilePool acquisition; each NC/KC panel uses its deterministic ordinal as the
`run_declared` identity, and that call joins every worker and drops every
tile-scoped arena before returning, which is the drain barrier. All writes
before that barrier target private staging, so cancellation leaves C
unchanged. The current TilePool worker lifetime uses joined
`std::thread::scope`, not yet an asupersync child scope; fs-exec documents that
precise L0 no-claim.
After the final poll, copying the
complete staged result into exclusively borrowed C is the non-cancellable
finalize step; a later request applies to the next operation.

Batched kernels remain bounded synchronous loops; chunking a large batch
to tile quanta (and Cx poll points between chunks) is the fs-exec driver's
job.

## Unsafe boundary
None. `unsafe_code` is denied workspace-wide; any future capsule must be
registered per docs/CONVENTIONS.md and ship a SAFETY.md.

## Feature flags
None. Frontier features use `frontier-*`, moonshots `moonshot-*`, default off.

## Conformance tests
In-crate GEMM suite: bitwise same-order oracle across shape sweep, β/α
edge semantics (β=0 NaN overwrite, α=0, k=0, empty m/n), transpose
identity, submatrix consistency, mixed == widened-f64 bitwise, f32
tolerance battery, determinism + golden hash; G4 deterministic
mid-dispatch cancellation injection proves bounded polling, real TilePool
traversal receipts, worker drain, arena quiescence, and unchanged C; G5 proves
the success path is bitwise the serial contract across 1, 2, host-parallelism,
and advisory-pinned pool configurations.
tests/conformance.rs placeholder remains for the shared-harness migration.
Batched battery (tests/batched_battery.rs): GEMM vs scalar oracle
across size classes + β-accumulate path; batch-membership bitwise
invariance for Cholesky and pivoted LU; L·Lᵀ reconstruction and solve
for k ∈ {4, 6, 12, 24, 48}; NotSpd/Singular members flagged while 8/9
healthy neighbors stay bit-correct; PA = LU residual + solve; det
closed forms cross-checked against LU diagonal-product-with-parity;
inv·A = I; eigh3 closed form vs per-matrix Jacobi as sorted sets plus
trace identity and eigenvector residuals, degenerate (isotropic and
diagonal) fixtures; 128-byte plane alignment; cross-ISA golden hash.

## Perf-lane evidence (bead xdgf, measured)
- Release, macos-aarch64 (Apple M4 Pro, Mac16,11), single thread,
  fs-roofline `MachineAxes::probe()` peak 51.7 GFLOP/s (register-file-
  sized probe): gemm_f64 through the fs-simd NEON 8×4 capsule
  microkernel reaches 43.9–45.1 GFLOP/s at n ∈ {256, 512, 1024}
  square = 85–87% of measured peak — the ≥75% roofline gate PASSES
  (`tests/perf_lane.rs`, best-of-3, 2mnk flop model). n = 128: 38.4
  (0.74 — blocking overheads at small n, reported not gated).
- The capsule (`fs_simd::ops().mk8x4_f64`, NEON `vfmaq_laneq`) is
  BITWISE-identical to the scalar twin per element (same k-ascending
  fused order), so the GEMM golden 0x1d7a_a3c6_b631_7ef0 is
  tier-invariant — verified by the in-crate golden test on both the
  capsule and scalar paths (aarch64 + the twin), and by fs-simd's
  equivalence battery (kc ∈ 0..17 ∪ {256}, special values, nonzero
  starting accumulators).
- The second-ISA (x86-64/AVX) capsule + attainment row are ARMED
  PENDING x86 hardware (fleet ARM-only by census; scalar twin
  dispatches there meanwhile, correct but ungated).

## Perf-lane evidence (bead 9ekv, measured, TARGET NOT MET)
- Batched small-dense (tests/batched_perf_lane.rs, M4 Pro, roofline
  intensity model 24k² B/elem, 2k³ flops/elem, ~50 MB working sets):
  after the 9ekv slice-1 optimizations (m-chunking, fs-simd 4×4-tile
  capsule, power-of-two stride padding, plane-vectorized LU updates —
  all bitwise-neutral, golden 0x0377_a8c9_5992_aee9 unchanged),
  batch_gemm reaches 5.5–15.3 GFLOP/s = 10–29% of each binding roof;
  batch_lu reaches 22–32% of its binding roof. The plan's ≥60% target is
  explicitly compute-peak-relative and is NOT met on this machine: rows report
  both `attainment` (binding roof) and `target_attainment` (compute peak), and
  only the latter decides `target_met`. The asserted lane gate is an
  anti-collapse floor of 0.08 against the binding roof, not the target; an
  environment-invalid row fails the evidence lane rather than yielding a green
  no-claim run. Successor design notes remain in bead 9ekv.

## No-claim boundaries
- Cargo does not expose the resolved feature set of transitive dependencies to
  `fs-la`'s build script. The canonical FrankenSim graph fixes fs-exec's
  asupersync features through its manifest, but a build that feature-unifies an
  additional asupersync capability must set `FRANKENSIM_GEMM_CODEGEN_ID`; tune
  reuse across such an unsalted non-canonical graph is not claimed. A future
  resolved-unit-graph receipt should replace this operator obligation.
- Compiler and wrapper executable bytes are watched and content-addressed, but
  a wrapper that changes generated code through mutable configuration outside
  the declared Cargo/Rust environment must set `FRANKENSIM_GEMM_CODEGEN_ID`.
  Tune reuse across an unsalted non-transparent wrapper configuration is not
  claimed.
- `GemmMemoryReport` is deterministic logical reservation accounting, not RSS
  or allocator-overhead measurement. It excludes TilePool slots, deques,
  victim tables, worker stacks, and each `fs_exec::RunReport`'s dynamic
  `cancel_latencies_ns` and `tiles_by_worker` vectors;
  the shared generic executor lease is
  `frankensim-epic-substrate-wf9.16`. Arena high-water counts an allocation
  attempt before calling the allocator so concurrent live reservations cannot
  be underreported; a refused attempt can therefore make the logical peak
  conservatively exceed physically committed bytes.
- Batched small-dense throughput is NOT claimed roofline-competitive
  yet (see 9ekv evidence above): in-kernel AoS repacking or per-matrix
  packed-GEMM routing for large k, the autotuned interleave, and the
  f32/mixed batched variants remain open 9ekv scope.
- Perf scope still open after xdgf slice 1 (recorded successors):
  autotuned MC/NC/microkernel-shape sweep (KC retune = golden bump),
  CCD-aware fs-exec parallel tiling (fz2.2 lane), f32/mixed packed
  paths + capsules, AVX-512 microkernel, nightly perf-regression
  history through fs-roofline::regress.
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
  driver yet.
- Batched small-dense makes NO roofline-competitive claim yet. The current
  packed 4×4 f64/f32 capsules are bitwise-equivalent execution paths, but the
  ≥60%-of-compute-peak acceptance, autotuned interleave width, and a
  lane-vectorized pivoted LU remain open. Batched QR/SVD and eigh3 closed-form
  EIGENVECTORS (currently via per-matrix Jacobi) are recorded follow-up scope.
- `batched_f32` (9ekv scope e): `BatchMatF32` + `batch_gemm_f32` (fused
  f32 chain) + `batch_gemm_mixed` (f32 storage, EXACT widen, fused f64
  chain, exactly ONE f32 rounding per output — the intended substrate
  for a future LBM moment path; no production `fs-lbm` consumer yet).
  Bit-deterministic and membership-invariant by the same construction as
  the f64 path; own golden `0x5600_7cfe_6a6d_1f9a` (registered against
  `fs-la:batched-f32-bits=2`; v2 fixes the promised α = 0 no-read path),
  verified identical debug+release on
  aarch64. NO performance claim: v1 is plain plane sweeps (no MBLK
  chunking, no packed tiles, no capsules) — the perf treatment joins the
  9ekv lane; cross-ISA golden row pending the next x86 run.
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
- Randomized-NLA error estimates are PROBABILISTIC indicators with
  empirically validated coverage, not certified bounds (fs-ivl owns
  certificates); SRTT sketches (via fs-fft), an LSQR driver, and
  e-process stopping integration are recorded follow-up scope.
- `condition_estimate` in `RefineReport` is a Hager-style ESTIMATE (not
  a certified bound); the ladder thresholds are engineering headroom,
  not proofs. Componentwise (per-entry) backward targets, Krylov
  inner/outer precision splits, and ledger event emission are recorded
  follow-up scope (fs-obs schema wiring + consumers do not exist yet).
