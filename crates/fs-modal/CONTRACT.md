# CONTRACT: fs-modal

## Purpose and layer
The first-class vibration eigenproblem facility: K·φ = λ·M·φ for symmetric
K and symmetric positive-definite M, with per-mode certified intervals and
inertia-certified spectrum slicing. Layer: **L1 BEDROCK** (peers: fs-la,
fs-sparse, fs-math — the same tier fs-cheb occupies). Consumers: the
plates/shells, viscoelastic-damping, and vibroacoustic beads of the
musical-acoustics program, and any structural-dynamics caller.
Bead frankensim-fsim-vibration-eig-jw6yq. Since IRP M1.1-W1a the crate also owns the time-domain realization of
mass-normalized damped modes (`time`, moved from fs-couple under owner decision OD10; it needs only `std` and `fs-math`).

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
- `recentred_dense_radii(pairs, k, m) -> (Vec<ModePair>, Vec<f64>)` —
  for a COMPLETE dense solve (every mode present), recentre each λ on its
  own vector's Rayleigh quotient in the original pencil, then return the
  [`quadratic_radii`] with the window certified by dimension. An f64
  eigenvector of an ill-conditioned pencil keeps a linear residual near
  `ε·λ_max`; the quadratic radius bounds the eigenvalue itself.
  `eigh_gen_dense` is unchanged.
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

### `time` (IRP M1.1-W1a; moved from `fs-couple::modal_acoustic_time`)

`ModalAcousticTimeModel` is the generic fixed-rate realization of independent,
mass-normalized structural modes. It consumes generalized force in
`N / sqrt(kg)`, advances viscously damped modal coordinates by the exact
zero-order-hold transition in underdamped, critically damped, and overdamped
regimes, and emits physical observer pressure in pascals. Acoustic observation
uses a caller-supplied complex pressure-per-modal-velocity transfer at each
natural frequency under `exp(-i omega t)`; no material names, digital gains,
or mastering values occur in this layer. Energy/work diagnostics and
transactional state/pressure budgets accompany every sample.
The damped transition prepares its held-force convolution with a scaled entire
series and transition doubling, then reuses scalar gains in the sample loop.
It does not subtract a potentially huge static equilibrium, snap nearby damping
ratios to critical damping, or subtract nearly equal overdamped poles. Separate
increments and retentions preserve both small attacks and long decay tails.
Focused G1 tests cover those regimes and a one-second, 440-cycle high-Q decay;
cached full/half steps and direct steps share the same coefficient evaluator.
`held_force_for_port_velocity` solves the exact-ZOH scalar endpoint constraint
using free modal motion plus unit-port-force mobility, with conjugate force
and velocity projections.
The force mobility is evaluated from a zero state, rather than subtracting
unit-forced and free motions about the accepted state. This keeps a finite
physical force response from disappearing beneath a large background velocity;
solving the constraint still does not advance or mutate the state. It refuses invalid/nonpositive mobility and
nonfinite inputs/results without changing state. The ordinary step still
admits state and energy budgets. The bowed runtime's set-valued contact uses
the same free-motion plus unit-force endpoint response, cached per bow station
(see the Bow bullet under the string paths).
G1 compares an independent oscillator solution
and demonstrates the old pinning error; G3 checks mixed-sign multimode ports
in all damping regimes and full/half/subsample durations. This enforces only
endpoint velocity: continuous no-slip motion, exact within-step event times,
and a passive contact discretization are not claimed.
Each frame additionally retains mode-local viscous work-minus-energy differences
and their scale-aware roundoff allowances. Negative loss beyond a mode's allowance
refuses the entire step, even when another mode could mask it in the aggregate.
The aggregate work/energy calculation is unchanged; summing the mode-local
differences can differ by floating-point roundoff. Neither an undamped mode's
energy drift nor a negative value within its allowance is physical heating.
The default initial state is zero displacement and velocity. A caller whose
window begins after a load has already settled may explicitly initialize the
exact static compliance `q = F/omega^2`, `qdot = 0`; using that operation to
erase a real load-on transient is outside the contract.
`observer_pressure_with_transfers_about_static_equilibrium` additionally
removes the exact held-force compliance `F/omega^2` before applying the
complex transfer's displacement-quadrature term. This prevents a static load
from being misreported as DC acoustic pressure while leaving free modal motion
unchanged; the caller must provide the same generalized ZOH force that drove
the structural transition.

The acoustic realization is explicitly narrow-band: one complex transfer value
per mode is exact only at that mode's natural frequency. Broadband radiation,
propagation delay, and feedback impedance require the separately tracked stable
rational-fitting path; this module does not claim them or passivity.

Ownership and visibility: the module needs only `std` and `fs_math` (`det` transcendentals), so it sits at L1 below every
consumer. Mode cap `MAX_TIME_DOMAIN_ACOUSTIC_MODES = 4096`; `try_new` refuses 0 modes and more than 4096 with
`InvalidInput { what: "mode count must be in 1..=4096" }`. `advance_exact_zoh` and
`ModalAcousticTimeModel::advance_exact_zoh_step` are public (they were crate-private in fs-couple) because the bow contact,
the compact radiator and fs-couple's broadband observer take single exact steps; `ModalAcousticWorkspace::frame()` exposes the
last published diagnostics. The opt-in `broadband` observer is not part of this crate; it stays in
`fs_couple::modal_acoustic_time::broadband`. Determinism: one-host bitwise; no cross-ISA claim is made for the music stack.

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
The certified-damping paths refuse with `DampingInput`,
`DampingNotUnderdamped`, `DampingClusterTooLarge` (naming the cluster) and
`DampingBlockSpansWindow`; they never fall back to a default budget.
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

## Certified modal damping (IRP M0.2-W5)
`certify_modal_damping(omega, C, budget)` takes mass-normalized modes with
undamped frequencies `omega` [rad/s] and a dense symmetric modal damping
matrix `C` [1/s]. It decides which modes may keep one energy-quotient decay
rate (`ModalDampingBlock::Diagonal`, `ζ_i = C_ii/(2ω_i)`) and which must be
integrated as a full block (`Full`).

The certificate is Rouché's theorem on the Schur complement of
`Q(λ) = λ²I + λC + Ω²`. The derivation is in `src/damping.rs`.
- A diagonal radius is second order in the off-diagonal coupling.
- A block's poles come from a dense solve of its energy-coordinate matrix
  and are certified against the complement (Gohberg–Sigal).
- Modes whose radius breaks the budget merge with their strongest coupling
  partner, union–find style.
- Certified discs lie in the open upper half-plane and are pairwise
  disjoint. Their counts therefore sum to n: every pole is enclosed exactly
  once.
- `ModalDampingBudget` has no default. Its decay and cents tolerances must be
  finite and positive, and the cluster cap must lie in
  `1..=MAX_DAMPING_CLUSTER` (64, the connection ceiling).

`certify_modal_damping_windowed` never forms C.
- It reads `C_jk` only for `|ω_j − ω_k| < window`. `damping_window_neighbours`
  returns exactly those sets.
- It bounds each row's unread entries with Bessel's inequality from a
  caller-supplied column-norm bound. For FE callers that bound is
  `√(xᵀM⁻¹x)`, `x = C_FE φ̂_j`. The difference of squares puts a floor of
  about `√ε‖C e_j‖` under every far bound.
- `certify_modal_damping_windowed_far_rows` takes the far-row bound itself.
  FE callers form it as the mass norm of the residual
  `M⁻¹x_j − Σ_{k∈near∪{j}} φ̂_k C_kj`, with a floor of about `ε‖C e_j‖`.
  `residual_far_row_bounds` does this for a lumped mass. It charges the
  measured defect of `ΦᵀMΦ = I`, the caller's bound on the rounding of
  `x_j = C_FE φ̂_j`, and first-order `γ_k` allowances. The plate owner's
  windowed projection in fs-couple calls it.
- `FarCouplingBound` reports the bounds used.

`certify_modal_damping_owned_diagonal(omega, C, zeta_owner, budget)` is for a
caller whose diagonal decay rates have their own owner, such as a plate's
energy-quotient image. It certifies on C as above. Each `Diagonal` block then
reports `zeta_owner[i]` bit for bit, with its closed-form pole, and its radius
grows by the centre shift so the exact pole stays enclosed. An owner ζ that
breaks the budget or moves a disc into another refuses by name.

Comparisons carry `ROUCHE_SAFETY = 1e-9` relative; this is not interval
arithmetic.
- Singleton centres are closed-form roots. Their rounding bound `2εν` enters
  every disc.
- Block roots are seeded by fs-la's dense complex QR, which can be ~1e-12
  relative off. They are refined by inverse iteration and the two-sided
  Rayleigh functional, and carry 8× their first-order forward-error
  estimate.

The certifier adds no projection: plates pass
`ModalPlateViscosity::matrix_per_s()`.

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
`tests/certified_modal_damping.rs` (W5) tests the certificate against
refined dense truth. That truth is `quadratic_eigenvalues`, refined by
inverse iteration and the two-sided Rayleigh functional; the raw values move
up to ~1e-12 relative.
- **G0, 4096 seeded random PSD systems** (n 4..48, ζ 1e-4..5e-2, coupling
  strength 1e-4..1, 25 % with a planted 1e-4 pair): zero enclosure
  violations, with the maximum and median enclosure ratios reported.
- **G0 Rouché order:** radius order 2 ± 0.2 under off-diagonal halving.
- **G1:** a clustered pair forms a full block within budget, and forcing the
  diagonal breaks the decay budget (reported).
- **G0 cap and inputs:** the 64/65 cluster cap, and refusal of every budget,
  frequency and matrix input.
- **G0 Bessel:** the far bound dominates the true far norm and falls to its
  √ε floor on a complete space.
- **G0 windowed:** windowed certificates enclose dense truth.

The plate's own matrix is tested in
`fs-couple::thin_plate::material_loss::ports::tensor_tests`.
`tests/modal_casebook.rs`: a 100,489-DoF grid pencil sliced over two
windows (one spanning degenerate pairs), counts certified against BOTH
inertia and the analytic spectrum, JSON-line evidence rows (shift,
factorizations, iterations, restarts, nnz(L), peak front bytes, delayed
pivots, max residual, wall time); plus the FrankenScipy `eigvalsh` oracle
cross-check on the dense strategy (deterministic seeded fixture, 1e-9
agreement gate).

## No-claim boundaries
- Certified damping covers lightly damped (underdamped) modes only. A real
  pole refuses, and so does a non-PSD diagonal. The radii are floating-point
  certificates with a stated margin, not interval arithmetic. A block pole's
  allowance is a first-order a-posteriori estimate, not an enclosure, so an
  ill-conditioned (near-defective) block can break it.
- Only material damping is certified. Radiation and thermoelastic losses
  remain separate owners. Prony/frequency-dependent `C(λ)` branches are not
  implemented here.
- The Bessel far bound has an inherent √ε relative floor (a difference of
  squares), so a "vanishing" far norm is reported as ~1e-7, not 1e-10. The
  windowed path's per-mode cost is O(n) scans plus local Perron bounds. No
  wall-time claim is made for 10⁴ modes, which is untested.
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

## Physical transfer between changing modal bases

`transfer::transfer_modal_state` projects displacement and canonical momentum
using the actual source shapes and one fixed SPD reference mass matrix.
Equilibrium offsets are included before projection. Mode order, eigenvector
sign and rotations within a degenerate subspace do not identify state.
The proposal reports original incremental modal storage, full target-pencil
storage, retained storage and independent displacement/kinetic/elastic
projection residuals. Explicit orthogonality, eigen-residual, energy, work,
state-loss and allocation limits gate publication; inputs are never mutated.

Parameter work is the signed change of represented incremental storage before
projection. It is not the complete prestressed body's thermodynamic energy.
Discarded state is neither heat nor physical damping. This stable,
positive-frequency slice requires the same coordinate space and exactly
symmetric numeric operators; changed mesh/support, rigid modes and
constitutive memory require their own transfer.

`transfer::transfer_modal_state_with_mass_change` additionally admits distinct
SPD source/target masses on that same physical space. It preserves displacement
and canonical momentum before truncation: `p = M_old v_old`,
`v_new = M_new^-1 p`. Source shapes are checked in the old mass inner product;
target shapes and displacement projection use the new mass. The target mass
factorization is reused for its solve, whose independent relative residual has
an explicit caller-supplied limit. Target kinetic storage and discarded velocity
energy use the new mass. The signed kinetic-storage change is parameter work,
as is the stiffness/equilibrium contribution. Holding velocity or renormalizing
the sound amplitude would violate this contract.

This is an instantaneous Hamiltonian parameter update. It does not model mass
accretion/ejection, momentum flux, a moving coordinate system, added-fluid
coupling or new constraints; those need their own source/transfer laws. No live
material consumer is inferred from this generic transfer capability. Existing
fixed-mass tests cover nonidentity mass, mode signs/order, degenerate rotation,
equilibrium shifts, independent energy accounts, explicit truncation and
refused/cancelled inputs. Four additional G1/G3 tests cover analytic scalar
inertia changes in both directions, a dense changed mass with rotated degenerate
modes, target-mass truncation energy, and malformed/cancelled updates.
