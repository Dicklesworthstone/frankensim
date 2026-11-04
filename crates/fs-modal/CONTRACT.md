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
  `(λ²M + λC + K)φ = 0`. Exact power-of-two coordinate and time scaling of
  the original coefficients precedes Cholesky mass whitening and the fs-la
  dimensionless companion QR path. No `M⁻¹K` with dimensional squared
  frequencies is formed. All nonproportional and gyroscopic damping entries
  remain coupled. Exactly diagonal original pencils instead use analytic
  scalar quadratics with product-of-roots recovery, retaining slow overdamped
  poles, complex conjugates, free roots and unstable poles. There is no
  tolerance-based decoupling. Values only; canonical `(real, imaginary)` order.
  `quadratic_eigenvalues_with_cancel` shares this exact owner, adding an explicit
  per-window QR sweep cap and cancellation through native complex QR. The
  original entry point retains its 120-sweep cap and numerical arithmetic.
- `quadratic_modes::PreparedQuadraticModes` retains the complete original
  M/C/K and prepares binary coordinate/rate scaling plus a mass Cholesky once.
  `refine` uses the shared fs-la bordered Newton owner and analytic quadratic
  derivative to resolve one nearby simple pole and its complex right shape.
  Published displacement is in original coordinates with `phi^H M phi = 1`;
  velocity is the actual `s phi`. Full mutual inertia, nonproportional damping,
  signed gyroscopic terms and active/unstable coefficients are preserved.
  Physical recovery is independently gated by original M/C/K row contractions,
  with the declared reference rate resolving zero-pole normalization. Explicit
  bytes/work budgets include owner caches and the solver, and cancellation
  publishes no partial mode. This is local refinement, not a complete spectrum,
  eigenvalue enclosure or forward shape-error bound. Right modes alone do not
  provide general force participation/residues or a diagonal forced-response
  realization; those require left-mode information. Unresolved repeated or
  defective poles refuse through the existing Newton owner. Seven native
  `tests/complex_modes.rs` regressions cover analytic nonproportional phase,
  gyroscopic handedness, full inertia, unit modal mass, coordinate/rate scaling,
  active growth, rigid and defective poles, budgets, cancellation and replay.
- `quadratic_modes::PreparedQuadraticModes::prepare_response` rechecks one
  resolved simple pole against the original source rows and factors the full
  bordered linearization with `fs-la::pole_response::PreparedPole` once.
  `PreparedQuadraticResponse::forcing_residues` accepts complete batches of
  analytic forces evaluated at that pole and linear displacement/velocity
  sensors in original coordinates. It returns physical Laplace displacement,
  velocity and sensor residues plus participation in the canonical unit-mass
  right-mode gauge. Left/right nonnormal coupling comes from the full solve;
  no Hermitian right-mode force projection or proportional-damping shortcut is
  used. Binary force scaling and the physical spectral chain factor preserve
  the original mass, rate and coordinate units. Each RHS reuses the factor and
  audits the original bordered rows. Separate numerical/physical byte and work
  ceilings admit complete batches, and cancellation publishes no partial batch.
  These are one simple pole's residues, not full harmonic amplitudes, powers,
  observability certificates or a complete modal realization. Seven native
  `tests/quadratic_residues.rs` cases cover analytic scalar and full 2x2
  resolvents, nonnormal zero-right-projection excitation, nonproportional phase,
  signed gyroscopic coupling, mixed units and extreme rates, complex sensors,
  simple zero poles, active growth, refusal gates and deterministic factor reuse.
- `quadratic_modes::PreparedQuadraticModes::prepare_sensitivities` borrows the
  original source pencil and prepares one existing bordered response factor.
  `PreparedQuadraticSensitivities::parameter_derivatives` accepts batches of
  complete signed physical `dM/dtheta`, `dC/dtheta`, `dK/dtheta`, with symmetric
  mass directions and full nonsymmetric damping/stiffness directions. It
  returns original-coordinate pole, displacement and velocity derivatives.
  Each analytic source action uses the shared fs-la parameter-derivative owner;
  `response()` reuses that same factor for force residues. Binary direction,
  coordinate and rate scaling retain the declared physical parameter units.
  Recovery differentiates `phi^H M phi=1`, including the moving mass, and holds
  the **baseline** original phase-reference component real. That coordinate
  stays fixed across a batch; no differentiability of a moving largest-component
  convention is claimed. Published values are independently gated by the
  differentiated original M/C/K rows and the unit-mass derivative constraint.
  Separate physical work/payload envelopes and the shared RHS ceiling admit
  complete batches; cancellation or refusal changes no model or prepared factor.
  These are local simple-pole tangents, not finite-step predictions, a mode
  tracker across crossings, shape-error enclosures or repeated-mode derivatives.
  Seven native `tests/quadratic_sensitivities.rs` cases cover analytic scalar
  tuning and mass normalization, nonnormal and signed gyroscopic motion, full
  mutual mass/nonproportional directions, independent neighboring pencils,
  extreme coordinate/rate/parameter units, zero poles, refusal and exact reuse.
- `ModalError` — typed refusals with stable `FS-MODAL-*` display codes:
  dimension mismatch, mass-not-SPD, factorization refusal (wrapping the
  fs-sparse `FS-SPARSE-DIRECT-*` code and the shift), invalid window,
  window-unresolved (the count-certificate mutation gate), unrepresentable
  mass-coordinate equilibration, quadratic-eig failure.

### `harmonic`: full forced response and resonance

`harmonic::PreparedHarmonicResponse` prepares the complete retained mechanical
pencil `K - omega^2 M - i omega C` under `exp(-i omega t)`. The source mass,
stiffness and viscous damping remain fully coupled; a mode count or a list of
free-decay poles does not substitute for the forced port response. `try_new`
accepts real exactly symmetric row-major matrices. `try_modal` supplies the
mass-normalized diagonal stiffness convenience while retaining full damping;
factored frequency detuning preserves narrow peaks and exact nominal resonance.

An unloaded `try_modal` pencil with exactly zero off-diagonal damping uses
compact diagonal storage. Every nonzero coupling, however small, and every
port-loaded pencil retains the dense path. Layout selection uses the source
damping, including at DC where its off-diagonals affect the frequency tangent.
The complete caller-owned damping matrix still receives finite/symmetric
admission and requires `O(n*n)` inspection. Subsequent compact assembly and
prepared response/tangent evaluations use `O(n)` work and owned storage.
No mode, damping term or physical coupling is approximated by this dispatch.

`try_diagonal_modal` accepts two length-n vectors for an explicitly independent
modal model: natural frequencies [rad/s] and viscous damping coefficients
`c` [1/s] in `q'' + c*q' + omega_i^2*q = F`. These coefficients are damping
rates, not damping ratios. Its authored input and preparation both use `O(n)`
space/work. It shares modal frequency validation, M/K assembly, factored
detuning and the compact physical solver with `try_modal`; it does not project
or infer off-diagonal damping or a port load. Coordinate/payload admission
precedes coefficient allocation, and assembly polls cancellation per mode.
Finite nonnegative rates, representable squared frequencies, zero-frequency
coordinates and singular-response refusal retain the same physical rules.
The same `272*n` payload allowance covers its `68*n` preparation peak; retained
storage and single-response peaks remain `44*n` and `92*n` bytes. Existing
frequency-derivative payload/work controls apply unchanged.

Exact binary coordinate/time scaling reuses the quadratic-pencil owner. Native
fs-la complex LU is cached per frequency for multiple force/receiver queries;
exactly diagonal pencils use scalar division. The response returns original
coordinate displacement and velocity phasors, componentwise backward residual,
and independent cycle-averaged input and viscous powers for peak amplitudes.
The forcing convention is `v = -i omega q`, `P_in = Re(F^H v)/2` and
`P_loss = Re(v^H C v)/2`. Source units and conjugate force/velocity projection
remain the caller's responsibility. Physical radiation loss appears here only
if its loading is actually included in the supplied pencil.

Mass must pass positive-definite admission. Exactly diagonal damping and
stiffness receive nonnegative-coefficient checks. Coupled matrices pass a
Jacobi PSD check after their own diagonal congruence, so a small indefinite
subspace cannot hide behind a larger unrelated mode. The numerical allowance is
`64 n epsilon`; this is not an interval passivity proof. Source coefficients
are never clipped, diagonalized or shifted. Zero damping and stiffness
nullspaces remain intact; singular lossless resonances and unconstrained DC
free motion refuse instead of manufacturing finite amplitudes.

Caller controls bound retained coordinates, peak numerical-array payload and
relative solve residual. Dense owners admit `144*n*n + 128*n` bytes; compact
modal owners admit `272*n` bytes without a dense eigensolver work limit.
The compact owner's numerical arrays occupy `44*n` bytes; its preparation
and single-response peaks are `76*n` and `92*n` bytes respectively. Caller
input, allocator overhead and thread stacks are excluded. Preparation/solves
poll cancellation between dense
admission stages, per assembly/diagnostic row and inside LU. Cholesky/Jacobi
admission is one indivisible bounded-size call. Inputs remain immutable and no
partial response escapes a failure. Repeat arithmetic is deterministic on one
host; no new cross-ISA claim is made. Backward error is not a near-resonance
forward-error bound. Modal truncation, nonlinear response, frequency-dependent
material validity and startup transients are outside this linear response.
G1/G3 tests cover analytic resonance phase/linewidth, full nonproportional
cluster damping, antiresonant port cancellation, a full two-coordinate pencil,
reciprocity, unit-coordinate changes over 100 orders, physical power balance,
free motion, singular limits, PSD refusal, budgets, cancellation and replay.
Compact-path regressions include a 257-mode bank under its linear payload
limit and a tiny physical coupling that requires dense admission and produces
nonzero motion of an unforced coordinate.
Explicit-vector regressions author 4,097 modes without a dense input, check
closed-form motion and power, compare response/frequency tangents with the
matrix-input API, and cover malformed, nonfinite, unrepresentable, singular,
resource-limited and cancelled requests.

`try_with_port_impedance` and `try_modal_with_port_impedance` attach the complete
frequency-dependent physical load `Z(omega)` through a real force map B. Flow is
`B^T v`, opposing effort is `Z flow`, and the dynamic operator receives
`-i omega B Z B^T` BEFORE factorization. No inverse of the bare structure or of
the load is formed: physical loss/reactance can resolve bare resonances, and
singular load channels remain admissible. Projection uses `O(n*p*(n+p))` work,
exact binary port/coordinate scaling and the original shared LU/scalar solve.
The complete complex Hermitian part `(Z + Z^H)/2` passes the same numerical
PSD owner through its real representation; mutual, reactive, nonreciprocal and
exact lossless channels are retained without clipping or a damping floor.

`HarmonicPortControls` separately bounds ports, projection work and additional
numerical payload. The response's `port_load` reports original physical flow,
effort and `Re(flow^H effort)/2`. The overall balance subtracts that work and
original structural viscosity exactly once. Net boundary work is not an
inferred split into heat, storage or radiation. Finite velocity impedance
contributes zero at DC; a singular static impedance needs its owner's physical
limit, not an implicit regularization. Source frequency support, causality,
broadband passivity, physical attachment units and fit uncertainty remain
upstream responsibilities. Pointwise PSD is not a broadband certificate.
G1/G3 tests cover analytic added-mass/stiffness resonance shifts, bare-pole
regularization, full complex mutual/nonreciprocal response, power-conjugate
unit changes, rank-deficient and lossless loads, tiny physical damping, complex
Hermitian refusal, physical work, budgets and cancellation. Existing unloaded
harmonic regressions retain their original behavior and `port_load = None`.

`solve_with_frequency_derivative` reuses that same prepared factor for the exact
continuous frequency tangent, including authored `dF/domega` and the complete
physical `dZ/domega`. Loaded owners require an explicit finite impedance-slope
matrix; a zero matrix declares a locally constant load. The method applies
`2 omega M + i C + i B Z B^T + i omega B Z_omega B^T` without constructing a
second dense operator. Full mutual/nonproportional entries and physical port
units remain. Independent binary normalization of the tangent RHS preserves
the original force and high-Q response scales. It returns `q_omega`,
`v_omega = -i(q + omega q_omega)`, physical port flow/effort slopes, and
independently differentiated input, original viscous and absorbed boundary
work. A second row audit reconstructs the differentiated equation from the
published physical motion with uncancelled source bounds.

`HarmonicFrequencyControls` admits added payload `8*n*p + 512*n + 384*p` bytes
and dense work score `64*(n*n+n*p+p*p)+256*(n+p)` before allocation. An unloaded
compact modal owner admits work score `320*n`, with the same `512*n` added
payload bound. Its current tangent/audit arrays peak at `208*n` additional
bytes. Both paths retain the admitted scaled mass/damping and reuse the factor.
No additional factor, finite difference, material law or peak search is created.
M/C/K, geometry and the port map remain fixed; load causality, constitutive
support and the supplied slope are the physical owner's responsibility.
Frequency derivatives of a net load impedance do not identify its heat/storage
split. G1/G3 tests in `tests/harmonic_frequency.rs` cover high-Q and DC limits,
causal inertance/compliance and Maxwell loads, full mutual/nonreciprocal loads,
complex source slopes, exact antiresonances, zero forcing, binary coordinate
and time changes, independent work, explicit load-slope refusal, resources,
cancellation and deterministic reuse. Compact-path tests cover a 129-mode
bank under the linear work budget, agreement with the dense pencil, and a
nonzero DC tangent from tiny off-diagonal source damping.

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

`damping_rates_per_s()` is the authoritative physical viscous coefficient
`c` [1/s] for every admitted coordinate, including ratio-created banks.
`try_new_with_damping_rates` takes complete finite nonnegative rates in
`q''+c*q'+omega^2*q=F`; for a dimensional mass and drag coefficient this is
`c=b/m`. It admits exactly zero stiffness with positive drag and zero direct
pressure transfer. Full/half/arbitrary held-force steps, port mobility solves,
linear ramps, harmonic forces, simultaneous tones and affine tone envelopes
all use that same admitted law and the existing independent positive heat
owners. No artificial frequency, fitted decay floor or duplicate loss is added.

For an explicit-rate bank, `modes().damping_ratio` is diagnostic metadata:
`c/(2*omega)` when positive and representable, exact zero only when `c=0`, and
positive infinity when a positive ratio is undefined or would overflow or
underflow. This makes a missed legacy consumer fail visibly instead of treating
a damped free coordinate as lossless. Generator assembly, coefficient equality,
power evaluation and metadata copying must use `damping_rates_per_s()`, not
reconstruct rates from those diagnostics. The bank's arbitrary-duration
`advance_exact_zoh_duration_step` previews its real admitted law; passing a
diagnostic mode to the standalone ratio-based step does not copy that law.

`with_damping_rates` replaces complete rates without changing state, energy,
frequencies, observer transfers, clock or budgets. Changed rates receive a new
shared mechanics kernel and invalidate all preparation bound to the old one;
identical rates share existing caches. `with_damping_ratios` explicitly restores
the ratio law and cannot retain drag at zero frequency. Both admission routes
require representable physical rates. Legacy ratio-created banks retain their
established transition and heat operation order. `tests/damping_rates.rs`
checks free decay across drag scales, retained tiny heat, independent work,
all force laws, mixed free/soft/elastic modes, metadata, mobility, cache identity
and transactional refusals. Work and heat allowances remain floating-point
estimates rather than certified enclosures.

`ModalAcousticTimeModel::apply_impulse[_into]` applies already resolved
generalized impulses in `N s / sqrt(kg)` at fixed modal displacement:
`v+ = v- + J`, independent of the nominal sample period, stiffness and viscous
damping. `apply_port_impulse[_into]` projects a physical `N s` impulse with
the same reciprocal mode-shape weights that reconstruct the port velocity.
Signed work is the impulse dotted with mean pre/post velocity; it is computed
without subtracting stored mechanical energies. Endpoint means, products and
the signed reduction retain double-double tails so cancelling modal work does
not erase small energy transfers. A nonzero impulse swallowed by velocity
rounding or a nonzero physical projection underflowing to zero refuses the
whole jump. No time elapses and all viscous heat entries are exactly zero.
Free coordinates use the same mass normalization.
The `_into` methods reuse the existing workspace and perform one bounded scalar
pass without allocation or transcendental functions. All existing state,
energy, pressure and work-balance admission gates remain transactional. These
methods consume a resolved impulse; they infer neither restitution nor contact
duration, support reactions, omitted modes or broadband impact radiation.
Finite-duration excitation remains a force-trajectory problem. G1/G3 regressions
in `tests/modal_impulse.rs` check clock/damping independence, physical mass and
reciprocal port work, braking, subdivision/reversal, retained small work and
whole-bank refusal, including cancellation tails and instantaneous pressure
limits that finite-time decay could otherwise hide.

The damped transition prepares its held-force convolution with a scaled entire
series and transition doubling, then reuses scalar gains in the sample loop.
It does not subtract a potentially huge static equilibrium, snap nearby damping
ratios to critical damping, or subtract nearly equal overdamped poles. Separate
increments and retentions preserve both small attacks and long decay tails.
Focused G1 tests cover those regimes and a one-second, 440-cycle high-Q decay;
cached full/half steps and direct steps share the same coefficient evaluator.

The undamped held-force cache also retains the entire gains
`sinc(omega*h)` and `(1-cos(omega*h))/(omega*h)^2`. If squaring frequency loses
its normal range, or a static equilibrium/velocity quotient overflows, motion
uses these gains with scaled force impulse instead of dividing by stiffness.
This preserves finite near-rigid attacks and their signed physical work without
introducing a frequency floor. The corresponding exceptional spring-force and
stored-energy evaluation multiplies `omega*q` first, retaining finite restoring
acceleration and energy-budget admission even when `omega^2` underflows.
Ordinary normal-range transitions and energy retain their existing operation
order. `tests/soft_mode_time.rs` covers the free limit, cached/direct and
subdivision agreement, finite attacks with unrepresentable static equilibrium,
retained restoring force across damping regimes, and transactional energy caps.

`time::LinearForceTransition` adds a prepared exact first-order-hold path for
forces specified at both interval endpoints. It shares the entire-series and
transition-doubling owner, adding the analytically integrated force moment.
Preparation depends only on frequency, damping and duration; advancing a state
allocates nothing and uses scalar arithmetic only. Undamped, critical,
overdamped and explicit free coordinates use one continuous formula without
division by stiffness or pole separation. Scaled forcing and increments retain
soft-mode attacks and small motions; separate retentions preserve decay tails.
Inputs and unrepresentable scaled responses refuse without state mutation.
`advance_with_energy` returns the same motion plus independently integrated
force work and positive viscous heat. Exact first/second force moments compute
the ramp's work without subtracting absolute positions or dividing by a tiny
duration. The heat owner extends its square-root Gramian to velocity, scaled
initial force residual and constant force slope. The same resolved eight-point
quadrature and QR doubling retain all ramp-induced damping; no stored-energy
subtraction, loss clipping or pole decomposition enters the result. Zero loss
returns exact zero heat. `try_new_with_damping_rate` accepts the physical
mass-normalized viscous coefficient `c=b/m` [1/s] in `q''+c*q'+omega^2*q=F`,
including positive drag at exactly zero stiffness. The same entire transition
and positive heat Gramian retain the free body's decay distance, terminal
velocity, ramp work and dissipation without inventing stiffness or forming
the singular ratio `c/(2*omega)`. The ratio constructor retains its zero-drag
free-coordinate semantics. Warm calls allocate
nothing and evaluate no transcendental functions.

`ModalAcousticTimeModel::prepare_linear_forcing` cold-prepares immutable ramp
caches shared by model clones. `step_linear_into` uses the existing workspace,
publishes actual ramp work and per-mode heat, and enforces the same transactional
state/pressure/energy and independent balance gates as held-force stepping.
It refuses missing preparation instead of allocating in the callback.
`advance_exact_linear_step` previews the identical cached motion for coupled
observers without advancing accepted state. The allocating `step_linear`
convenience prepares and returns an owned frame. All model ramp intervals use
the nominal clock; standalone transitions support arbitrary positive duration.

The input is the piecewise-linear force itself; interpolation error in measured
forcing, finite-arithmetic phase error over very long intervals and a certified
roundoff enclosure are not claimed away. G1 compares independent 100-digit
augmented matrix exponentials, affine force balance, the free cubic limit and
resonant growth under a sampled linear forcing interpolant. G3 checks temporal
subdivision and replay, including both sides of critical damping and a strongly
overdamped slow pole. Independent integration of the physical power equations,
prescribed affine motion, heat below energy-subtraction roundoff, signed work
at large free-coordinate offsets, and additive ramp subdivision cover the new
energy accounts. Runtime regressions cover cold preparation, exact preview,
shared clone caches and transactional refusal. Existing ZOH stepping and energy
arithmetic are unchanged.

`time::HarmonicForceTransition` advances the same scalar oscillator under a
continuous load `Fc cos(Omega t) + Fs sin(Omega t)`. It shares the homogeneous
transition and prepares the exact force convolution by an entire series and
block composition. No division by detuning or pole separation occurs: exact
undamped resonance retains secular growth, and free/critical/overdamped modes
use the same formula. Its `try_new_with_damping_rate` constructor also admits
explicit viscous drag at zero stiffness under a tone or DC load, with the same
mass-normalized `c=b/m` units as the ramp owner. A separately evaluated cosine-minus-constant remainder
retains tiny resolved velocities about static equilibrium near zero forcing
frequency. Warm motion, phase advancement and zero-state harmonic responses
allocate nothing and call no transcendentals. Multiple harmonics superpose
after advancing the homogeneous state once; the public coupled-time reexport
exposes the same owner. Invalid or unrepresentable trials return no motion.
G1 checks analytic secular resonance, continuous complex amplitude and phase,
the free-coordinate and DC limits, and independent 115-digit dimensional
matrix exponentials; G3 checks Fourier superposition, temporal subdivision
and bitwise checkpoint replay. This is prescribed linear forcing, not a
nonlinear harmonic-balance solver. It supplies no work/heat or broadband
radiation account, and exact endpoints do not resolve unsampled observations.
Long oscillatory intervals and repeated cached phase rotations accumulate
floating-point error; no amplitude normalization is applied.

`time::PreparedHarmonicForce` adds the same independently integrated signed
work and positive viscous heat to a standalone scalar tone transition. Its
ratio and explicit-rate constructors share the existing motion owner and
`HarmonicIntegrals`; a free viscous coordinate therefore retains physical heat
and force work under DC or continuous harmonic forcing. It returns both
work/heat roundoff allowances, preserves cosine/sine and transient/forced
cross terms, and never derives heat from endpoint-energy subtraction. Warm
advancement allocates nothing and evaluates no transcendentals. The smaller
motion-only cache remains unchanged. `PreparedModalHarmonic` uses this same
prepared scalar type internally, retaining its existing admission and
transactional balance gates. Scalar preparation does not infer model budgets
or broadband radiation, and separate-tone heat cannot be summed to substitute
for the simultaneous-tone owner. Analytic free-drag power integrals, tiny
positive heat, negative DC source work, subdivision and bitwise bank agreement
are exercised in `tests/harmonic_energy.rs`.

`time::PreparedModalHarmonic` connects that motion owner to the actual shared
`ModalAcousticTimeModel::step_harmonic_into` runtime. Every retained mode has
its own cosine/sine force at one prescribed frequency, including DC. The
ordinary frame reports independently integrated input work `integral F*v dt`
and nonnegative viscous heat `integral 2*zeta*omega*v^2 dt`; harmonic work is
not approximated by endpoint force times displacement. Preparation integrates
on a resolved short interval with eight-point Gaussian quadrature, then
composes signed work forms and positive square-root heat factors. Residual-force
coordinates retain near-static motion and cosine/sine interference; heat is
never computed by subtracting endpoint mechanical energies. Arithmetic
allowances are numerical estimates, not rigorous interval certificates.
Stepping allocates nothing, calls no transcendentals, and shares existing
modal/aggregate state, energy and balance gates. Rejected steps commit no
modal state. Preparation is bound to the model's immutable modal kernel and
may be reused by its clones; checkpoint replay also retains caller-owned force
phase. The model's nominal Nyquist guard admits the prescribed frequency.
Nonzero legacy natural-frequency pressure transfers refuse: driven motion and
startup transients need a separate physical broadband observer, so this
mechanical frame reports zero pressure. Simultaneous distinct-frequency power
accounts are outside this API; summing separate single-tone heats would omit
cross terms. G1/G3 tests cover independent continuous-power quadrature,
resonance, critical/free/DC limits, high Q, finite-interval interference,
multimode subdivision/replay, tiny positive heat beneath a large lossless
energy background and transactional refusals. Existing ZOH gates and exact
motion regressions remain shared.

`time::PreparedModalTones` extends the shared physical step to simultaneous
continuous tones through `ModalAcousticTimeModel::step_tones_into`. Frequency
ratios need not be integers. Every mode receives a caller-ordered cosine/sine
row; the full force array is mode-major. One homogeneous state plus the shared
scalar convolutions supplies motion. Exact duplicate frequencies are coalesced
with compensated force sums; distinct nearby frequencies retain their beats.
The power owner reuses the resolved scalar flow and QR rotations. Applying
the same rotations to velocity and force rows retains their cross integral as
paired factors: opposing contributions cancel before products, so small beat
work survives without subtracting positive quadratic forms. Squared velocity
factors give nonnegative heat, including every cross-frequency term. Neither
account comes from endpoint energy subtraction or separate single-tone heats.

`ModalToneBudget` bounds input tones (hard maximum 32), owned coefficient and
preparation-scratch payload, and explicit preparation/step algebraic work
scores before allocation. These scores are not wall-clock guarantees. Warm
steps and mode previews allocate nothing and evaluate no transcendentals;
the existing single-tone fast path remains available. Model binding, nominal
Nyquist admission, transactional physical gates, externally checkpointed phase
and separate broadband observation retain the single-tone contract. Power
allowances remain numerical estimates; extreme cancellation and long repeated
phase rotations retain finite-arithmetic limits. G1/G3 tests compare independent
continuous-power integration across DC/inharmonic/critical/free/high-Q cases,
close destructive beats, exact duplicate cancellation, resonance, multimode
subdivision/replay, tiny near-DC heat, the active 32nd tone and budget refusals.

`LinearForceTransition::pullback` returns `LinearForcePullback` cotangents for
initial displacement/velocity and both endpoint forces, with frequency, damping
and duration fixed. It transposes the same cached flow and force moments;
near-identity corrections reuse the independent decrements, while long decays
reuse the homogeneous retentions. It allocates nothing, evaluates no
transcendentals and does not differentiate preparation iterations. Inputs and
scaled results must remain finite, with no mutation or partial gradient on
refusal. A reverse trajectory carries both state cotangents and adds both
adjacent contributions at every shared force sample. These are derivatives of
the caller's scalar objective, not physical state values. Parameter/time
sensitivities and derivatives of nonsmooth peak selection are outside scope.
G3 checks the transpose identity across free/under/critical/overdamped regimes
and tiny/finite steps, then compares all initial-state and shared-force
trajectory derivatives with independent finite differences. G0 covers invalid
seeds, overflow refusal, zero cotangents and deterministic retry.

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
Each frame retains independent mode-local integrals of physical viscous power
`2*zeta*omega*v(t)^2`, with scale-aware contraction roundoff allowances. Their
sum supplies the aggregate heat. Undamped and free modes report exactly zero;
high-Q losses are not recovered by subtracting large mechanical energies.
Per-mode and aggregate force-work/energy balance checks remain separate and
refuse disagreement beyond their floating-point allowances before state commit.
The held force's work is taken over the displacement increment the transition
itself adds to the stored displacement (free, undamped and short damped steps
carry it; the long damped step forms `q - q0`), never over `next - state`,
whose rounding the force would multiply; the allowances are unchanged.

A positive square-root Gramian is prepared beside each exact full/half-step
transition. The same scalar transition owner resolves a short interval before
composing its viscous integral across longer or stiff intervals. QR preserves
positivity without subtractive determinants or clipping. The eight-point
quadrature is confined to a dimensionless generator norm at most one half,
where its remainder is below floating-point rounding. Warm stepping reuses the
prepared factors without allocations; arbitrary durations prepare the same
integral in fixed-size local storage. Existing state-transition arithmetic and
the narrow-band acoustic observer remain unchanged. Six focused regressions
cover independent heat/work integration, high-Q cancellation, very soft forced
modes, long/stiff composition, exact zero-loss cases, and cached/subsample
consistency.
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

### Affine continuous tone envelopes

`PreparedModalTones::try_new_affine` and
`ModalAcousticTimeModel::step_tones_affine_into` admit
`(C0+(C1-C0)*t/h)*cos(Omega*t) + (S0+(S1-S0)*t/h)*sin(Omega*t)`.
Both endpoints use the same interval-start carrier basis and original
mass-normalized force units. Rotate the accepted endpoint pair before the
next interval. The caller retains carrier phase through zero amplitude;
neither a physical step nor a refusal changes caller-owned phase.

The original scalar entire convolution gains two normalized endpoint-delta
columns; no detuning denominator or force/time-rate representation is added.
Motion retains one homogeneous state per mode. Cross-tone work and positive
viscous heat use the original paired force/velocity QR owner with four
rotating source coordinates per distinct tone. The residual flow includes
the envelope derivative as well as carrier rotation. Exactly equal source
frequencies coalesce each endpoint in original order before differencing;
nearby tones remain independent and retain their interference.

The constant two-coordinate preparation and step remain available unchanged.
Mixing preparation and step families refuses. Affine preparation remains
bound to the original immutable modal kernel and requires zero legacy
natural-frequency observation maps. Original physical-state, energy and
independent work/heat balance gates apply before accepted state publication.
Cold byte/work admission includes the larger complete source blocks and
their power factors; no tone or mode is removed. Warm propagation reuses
cached coefficients without allocation or transcendental evaluation.
Finite-arithmetic and Gaussian-composition allowances remain numerical
estimates, not interval certificates or instrument-validity claims.

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
`DampingNotUnderdamped`, `DampingClusterTooLarge` (naming the cluster),
`DampingAccuracyUnresolved` (naming the block, pole and numerical radius), and
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
- Every returned pole uncertainty meets both requested accuracy limits,
  including isolated rows and full-system blocks with no complement.
  Retaining all coupling does not waive finite-precision accuracy: an
  unresolved pole refuses with `DampingAccuracyUnresolved` instead of
  returning a certificate whose own reported error exceeds its budget.

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
An exactly undamped independent channel preserves its exact pole `i*omega`
and zero radius, including through this owner-diagonal path; it needs no
artificial damping floor. Repeated independent frequencies retain their
analytic poles and multiplicities in one block without a dense eigensolve.
A block with exactly zero coupling to its complement needs no perturbation
allowance from unrelated modes.

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
- **G0 accuracy:** isolated rows and full-system fallback blocks enforce
  both error limits; a stricter-than-resolved request refuses deterministically.
- **G1 lossless limit:** owner-diagonal material channels with zero loss
  retain exact imaginary poles and zero radii alongside damped channels.
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
- `quadratic_eigenvalues` returns values only (no eigenvectors or pole
  certificates). Its matrix path removes arbitrary coordinate and time-unit
  scales, retaining nonproportional loss and gyroscopic coupling; analytic
  G1/G3 regressions cover rates `1e±200`, coefficient scales `1e±250`, and
  signed coordinate scales `1e±100`. This does not cure intrinsic mass
  ill-conditioning, defective poles or extreme separation within one coupled
  spectrum. Exactly diagonal pencils retain the scalar strong-overdamping
  path. Scaling refuses non-finite or erased nonzero coefficients and
  non-finite or entirely underflowed rescaled poles. The general QR result is
  not a certificate that every slow pole was resolved; the analytic diagonal
  path additionally refuses an erased nonzero scalar root.
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

## Rigid and free-motion state transfer

`transfer::transfer_modal_state_with_rigid_modes` explicitly admits zero-frequency
motion alongside oscillatory modes on the same physical coordinate space. The
complete target rigid subspace must be supplied with exactly zero eigenvalues,
M-orthonormal shapes and a declared relative K-null residual bound. The existing
sparse LDLT owner qualifies the original stiffness on its complement through
`[K, M Z; Z^T M, 0]`, requiring n positive and r negative pivots in the bordered
inertia. Congruence scaling applies only to this qualification matrix; physical
stiffness, mode eigenvalues and storage are never shifted or regularized. Missing
null directions, false rigid modes outside their residual bound, and unstable
elastic complements refuse. Border entries and original numeric factor payloads
have explicit independent budgets.

The original full-mass canonical-momentum transfer and projection energy/work
accounts are reused. Rigid coordinates contribute kinetic but zero elastic
storage. The existing positive-frequency APIs retain their admission. Three
additional native regressions cover translation/rotation and degenerate rigid
mixing, a fully free target, changed consistent inertia and reversible work,
independent physical momentum, explicit elastic truncation, instability,
incomplete nullspaces, factor/border budgets and cancellation. This remains an
instantaneous fixed-frame update; it does not infer support impulses, finite
rotations, remeshing or constitutive-history transfer.


## Phase-continuous prescribed frequency modulation

`time::harmonic::ToneModulationStep` prepares one immutable, allocation-free
proposal for a force with linearly varying instantaneous frequency and signed
gain. It integrates the desired source phase and supplies four endpoint
coefficients to the existing affine tone propagator, all in the same
interval-start carrier basis. `advance_force_phase` advances the unscaled
source through zero gain, release and re-attack only when the caller accepts
the complete physical/observation step. No natural mode, material parameter,
mechanical state or receiver memory is retuned or reset.

The fixed-carrier affine reconstruction has a dimensionless complex-envelope
error bound from `h^2 sup |z''|/8`; callers multiply it by their original force
coefficient norm and admit an explicit physical error allowance. Frequency,
gain and interval inputs and all resulting phases/bounds must be finite and
representable. Source users retain the original clock's frequency guard and
bandwidth obligations. This is exact propagation of a reconstructed force,
with second-order approximation to the prescribed chirp; the estimate does
not enclose floating-point roundoff or certify bandlimiting. String retuning
still requires the physical tension/geometry owner.

`tests/tone_modulation.rs` checks signed-gain/DC phase and envelope bounds,
zero-detuning parity, 1,000 silent phase advances followed by re-attack,
invalid/overflow refusals, and a 32/64/128 interval chirp crossing resonance
against an independently integrated physical ODE. Displacement, velocity,
signed input work and positive heat each converge at second order while the
existing reconstructed-force energy account closes.
