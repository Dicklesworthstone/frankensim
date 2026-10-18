# CONTRACT: fs-modal

## Purpose and layer
The first-class vibration eigenproblem facility: K·φ = λ·M·φ for symmetric
K and symmetric positive-definite M, with per-mode certified intervals and
inertia-certified spectrum slicing. Layer: **L1 BEDROCK** (peers: fs-la,
fs-sparse, fs-math — the same tier fs-cheb occupies). Consumers: the
plates/shells, viscoelastic-damping, and vibroacoustic beads of the
musical-acoustics program, and any structural-dynamics caller.
Bead frankensim-fsim-vibration-eig-jw6yq.

## Public types and semantics
- `slice_window(k, m, (low, high], opts) -> SliceReport` — the certified
  sparse front door. Certifies the exact in-window eigenvalue count as
  `neg(K − high·M) − neg(K − low·M)` via `fs_sparse::direct` Sylvester
  inertia, then harvests exactly that many eigenpairs by shift-invert
  Lanczos in the M-inner product (full reorthogonalization, deterministic
  integer-LCG start vectors, RESTART-WITH-DEFLATION so degenerate clusters
  are recovered), and certifies each value with the M⁻¹-norm residual
  bound. One symbolic analysis serves every shift (union pattern).
  `SliceOptions.mass_diagonal_equilibration` is opt-in (default `false`):
  the sparse solve uses congruent coordinates `x = D y`, where
  `D_ii = 1/sqrt(M_ii)`. It preserves the pencil eigenvalues and inertia;
  returned modes and residual certificates are computed against original
  `K` and `M`. Unrepresentable scales or transformed entries refuse with
  `FS-MODAL-EQUILIBRATION` rather than changing pivot tolerances.
- `ModePair { lambda, phi, residual, interval }` — `φᵀMφ = 1` enforced;
  `interval = λ̂ ± ‖Kφ − λ̂Mφ‖_{M⁻¹}` contains at least one TRUE pencil
  eigenvalue (SPD-M reduction argument; see crate docs).
- `SliceReport { window, below_low, below_high, expected, modes, stats }` —
  `modes.len() == expected` is an invariant of every returned report;
  `below_*` are the raw inertia counts so downstream slicing consumers get
  certified evidence, not a boolean. `SliceStats` records the shift used,
  factorization count, Lanczos iterations, deflation restarts, and the
  factor's nnz(L)/peak-bytes/delayed-pivot numbers.
- `eigh_gen_dense(k, m, n) -> Vec<ModePair>` — dense strategy: fs-la
  Cholesky reduction (C = L⁻¹KL⁻ᵀ), cyclic Jacobi, back-transform,
  mass-normalization, the same M⁻¹-norm certificate. All n modes ascending.
- `shift_invert_modes_mfree(n, σ, want, budget, apply_m, solve_shifted)` —
  the matrix-free core (caller brings the shifted inverse): returns
  mass-normalized Ritz pairs nearest σ WITHOUT certified intervals or count
  certification (no M-solve, no factorization — recorded no-claim).
- `quadratic_radii(report, k, m) -> Vec<f64>` — per-mode Kato–Temple radii
  `‖R‖²/δ` for a certified slice. The radii are additive: the report's
  linear `interval`s and every digest over them are unchanged. The gap `δ`
  comes from the slice's exact in-window count:
  - Kahan's matching uses the full-block Frobenius residual;
  - clusters are runs of overlapping widened intervals;
  - the neighbouring clusters and the `[low, high)` window edges bound `δ`.

  A multi-mode cluster is tightened only with M-Gram defect
  ≤ `QUADRATIC_GRAM_DEFECT_MAX` (Mathias 1998). Its measured Galerkin defect
  `e = ‖QᵀKQ − diag θ‖_F` shrinks the gap to `δ − e`, and `e`, `2η|θ|` and a
  4-ulp rounding term are added. Otherwise, or whenever a widened interval
  leaves the window, the linear residual is kept.
- `quadratic_eigenvalues(m, c, k, n) -> Vec<C64>` — all 2n eigenvalues of
  `(λ²M + λC + K)φ = 0` via companion linearization `[[0, I], [−M⁻¹K,
  −M⁻¹C]]` and the fs-la dense complex QR path. Values only.
- `ModalError` — typed refusals with stable `FS-MODAL-*` display codes:
  dimension mismatch, mass-not-SPD, factorization refusal (wrapping the
  fs-sparse `FS-SPARSE-DIRECT-*` code and the shift), invalid window,
  window-unresolved (the count-certificate mutation gate), unrepresentable
  mass-coordinate equilibration, quadratic-eig failure.

## Invariants
1. A returned `SliceReport` satisfies `modes.len() == expected ==
   below_high − below_low`; a harvest that cannot meet the certified count
   REFUSES (`WindowUnresolved`) — never a silently short list.
2. Every returned `phi` satisfies `|φᵀMφ − 1| ≤` numerical roundoff
   (renormalized explicitly after harvest; tested at 1e-10).
3. Every `interval` is `λ̂ ± ‖Kφ − λ̂Mφ‖_{M⁻¹}` with the residual computed
   EXPLICITLY from K, M, φ (never from the Lanczos recurrence estimate).
   Analytic fixture values fall inside their intervals (tested on chains
   and the squared-Laplacian plate identity).
4. Metamorphic: slicing commutes with spectral shift (eig(K + cM, M) =
   eig(K, M) + c) and stiffness scaling (eig(sK, M) = s·eig(K, M)) within
   certificate widths (tested).
5. Repeat runs are bitwise identical: deterministic start vectors, fixed
   reorthogonalization order, `total_cmp` sorts (tested).

## Error model
Structural misuse (wrong slice lengths in the dense/quadratic paths) panics
with structured messages — programmer errors. Everything else refuses
through `ModalError`. Window endpoints that graze an eigenvalue make
(K − σM) exactly singular; the direct solver's refusal surfaces as
`Factor { shift, .. }` and the caller re-picks the endpoint — refusal beats
a fabricated count. Interior shifts retry three deterministic golden-ratio
fallbacks before giving up.

## Determinism class
Bit-deterministic across repeat runs on one host: sequential code, integer
LCG start vectors (no platform libm feeds solver state), index-ordered
reorthogonalization, `total_cmp` ordering. The pre-existing AMD digest
`tests/cross_isa_golden.rs::GOLDEN_HASH = 0x11dced94c7f67115` records the
aarch64/x86-64 debug comparison from 2026-08-23; IRP .67 registers that hash
without changing it, and requires fresh serial/parallel debug/release proof.
Underlying factorization determinism is inherited
from `fs_sparse::direct` (bitwise, tested there).

## Cancellation behavior
`slice_window` retains its context-free API and arithmetic. `slice_windows`
takes an `fs_exec::Cx`: it polls around factorization/analysis calls, at
Lanczos restarts, during certification and before returning reports. Native
windows use the existing executor bridge carrying the original cancellation
gate into scoped workers; the executor drains before returning. A cancelled
operation returns no partial multi-window report. Factorization, symbolic
analysis and an individual Lanczos sweep remain indivisible cancellation
units, so bounded-time cancellation within these calls is not claimed.

## Parallel windows and resource admission (IRP M1.5-W1)
`WindowPlan` requires ordered disjoint/adjacent windows, thread count, window
count cap, numerical-payload byte cap and relative certificate half-width.
`SliceOptions.ldlt.max_factor_bytes` must be explicit for parallel work.
Optional borrowed mesh information is owned once in `WindowPlan.coordinates`;
only ND consumes it. One mass analysis/factor and one shifted-pencil symbolic
analysis serve every window. Every distinct endpoint is factored once and its
inertia reused, so adjacent counts telescope exactly. Signed-zero endpoints
share one count. Endpoints coinciding with an eigenvalue continue to refuse.
Each window uses the same serial slicing body, seeds and accumulation order;
results are collected in window order independently of worker scheduling.
Wasm executes windows sequentially through the same body without spawning.

`MultiSliceReport` returns all slices, symbolic nnz(L), and timing/work/memory
accounting. `maximum_relative_halfwidth` is enforced against every explicit
residual; an unresolved or over-wide mode refuses the whole call. The memory
cap admits retained CSR/symbolic/mass payload, concurrent numeric factor caps,
Krylov scratch, scales and output vectors. It excludes allocator overhead,
thread stacks and transient symbolic-construction storage; reported peak is
an admitted conservative numerical payload, not measured process RSS. No
no-swap or hard whole-process memory guarantee is claimed.
The release-only tests contain the full 104544-DOF coupled P1 elastic block
and 1000000-DOF seven-point grid, not reduced stand-ins. Both must produce
nonempty counted windows with all relative half-widths <=1e-8 inside their
declared payload budget. P1/TET10 fill is separately reported. These native
scale/timing, 1/2/7-worker identity and wasm gates are batch-pending.

## Unsafe boundary
None. `#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags
None. The dev-only FrankenScipy oracle lives in `[dev-dependencies]`
(layer-exempt) and never enters the runtime graph.

## Conformance tests
In-crate: spring-mass chain window with certified count, interval
containment, and mass-normalization; squared-Laplacian plate identity
window spanning a DEGENERATE pair (asserts the deflation restart fired);
dense-vs-analytic and dense-vs-sparse cross-strategy agreement; shift
invariance and scaling covariance metamorphics; mass-not-SPD refusal;
starved-budget `WindowUnresolved` (the count-certificate mutation gate:
a harvester that misses modes cannot return success); invalid and empty
windows; matrix-free core vs assembled path; quadratic light-damping roots
vs exact per-mode formulas; bitwise rerun determinism.
`tests/modal_casebook.rs`: a 100,489-DoF grid pencil sliced over two
windows (one spanning degenerate pairs), counts certified against BOTH
inertia and the analytic spectrum, JSON-line evidence rows (shift,
factorizations, iterations, restarts, nnz(L), peak front bytes, delayed
pivots, max residual, wall time); plus the FrankenScipy `eigvalsh` oracle
cross-check on the dense strategy (deterministic seeded fixture, 1e-9
agreement gate).

## No-claim boundaries
- SPD M is REQUIRED; semidefinite or indefinite mass (free-free with a
  singular lumped mass, constraint-coupled formulations) is refused, not
  approximated. Buckling-style pencils (K, −K_G) are not this API.
- The matrix-free path certifies nothing: no inertia counts, no intervals
  (it has no factorization and no M-solve). Only the assembled
  `slice_window`/`slice_windows` paths carry certificates.
- The certificate bounds |λ_true − λ̂| for AT LEAST ONE true eigenvalue per
  interval; it does not by itself prove interval-eigenvalue pairing when
  intervals overlap, and it says nothing about eigenVECTOR accuracy
  (Davis–Kahan-style shape bounds are not claimed).
- `quadratic_eigenvalues` returns values only (the fs-la complex path has
  no eigenvectors), uses dense linearization only, and makes no claim of
  structure-preserving scaling for badly conditioned (M, C, K) — the
  viscoelastic bead's light-damping regime is the tested envelope.
- No LOBPCG/AMG tier for problems too large to factor: `slice_window`
  requires the sparse LDLᵀ to fit. That tier joins a follow-up bead when a
  consumer presents a pencil the direct factorization cannot hold; the
  matrix-free core is the seam it will plug into.
- Windows are half-open (low, high]: an eigenvalue EXACTLY at an endpoint
  makes the endpoint factorization singular and refuses (tested behavior in
  the fs-sparse casebook); callers nudge endpoints.
- Performance: the 100k-DoF casebook records wall time and memory evidence
  for its fixture; no throughput claim beyond it.
