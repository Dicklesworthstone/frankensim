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
  window-unresolved (the count-certificate mutation gate), half-width
  exceeded (a parallel plan's converged mode over its certificate cap),
  unrepresentable mass-coordinate equilibration, quadratic-eig failure.

### `relaxation_modes`: complete Maxwell memory spectra

`MaxwellDampingBranch { strain, stiffness, relaxation_time_s }` describes one
original scalar strain channel `e_j = b_j^T q`, with n finite strain entries,
positive stiffness `k_j` and positive relaxation time `tau_j`. Omit branches
with zero stiffness. `maxwell_eigenvalues` and `maxwell_eigenvalues_with_cancel`
retain the full dense mechanical M/C/K and every branch through
`M q'' + C q' + K q + sum(b_j s_j) = 0` and
`s_j' = k_j b_j^T q' - s_j/tau_j`. M must be symmetric positive definite;
nonproportional and signed gyroscopic C entries remain coupled. K is the
**relaxed** stiffness: the branches dynamically supply their unrelaxed
stiffness, which must not also be added to K. Eliminating memory gives
`C(lambda) = C + sum(k_j tau_j b_j b_j^T / (1 + lambda tau_j))`.

Both entry points return all `2*n + branches.len()` state poles, including
internal material relaxation, in physical [1/s] under `exp(lambda*t)`, sorted
by real then imaginary part. Energy normalization `z_j = s_j/sqrt(k_j)` makes
velocity/memory coupling antisymmetric using the same rounded coefficient
on both sides. Branch storage is `s_j^2/(2*k_j)` and physical heat rate is
`s_j^2/(k_j*tau_j)`; gyroscopic coupling does no work. The symmetric part of
the supplied C is not certified positive semidefinite, and active source
coefficients are retained.

Exact binary coordinate/time scaling and Cholesky mass whitening feed the
existing fs-la dense complex QR owner. With no branches, the original
quadratic owner preserves its analytic scalar/diagonal paths. A free scalar
coordinate with one Maxwell branch factors out its exact zero pole and uses
the same robust quadratic root owner on
`m*tau*lambda^2 + (m+C*tau)*lambda + C+k_b*b^2*tau`.
Binary-factored source products and signed DD coefficient sums preserve weak
relaxation and active C before common scaling; unrepresentable coefficients
or poles refuse. A passive scalar resonator with one Maxwell branch uses a
bounded contraction when `C/(m*omega) <= 1` and `k_b*b^2/K <= 1/8`. It solves
the relaxation-pole shift directly and recovers the conjugate pair from its
sum/product identities, preserving tiny mechanical decay without trace
subtraction. The full cubic is retained; the bounds select a stable numerical
route, not an approximation of the material law. Binary-factored dimensional
products retain physical decay even when its dimensionless shift underflows;
unrepresentable physical poles refuse. G1/G3 tests in
`tests/relaxation_modes.rs` cover manufactured finite-loss roots, relative weak
decay down to 1e-300, separated relaxation rates and changes of units.

For two or more arms on a passive scalar resonator, `C/(m*omega) <= 1` and
`sum(k_j*b_j^2/K) <= 1/16` admit a bounded fixed point of the **original**
rational pencil. The dimensionless pair variables `d=decay/omega` and
`P=(decay^2+frequency^2)/omega^2` remain in `[0,0.55] x [1,1.1]`; the map's
infinity-norm derivative is below 0.4. Positive factored dimensional loss
sums recover the pair without trace subtraction or erasing loss when an
intermediate dimensionless ratio underflows. Each distinct active source
rate contributes a bracketed real material pole; repeated-time stress
combinations and zero-strain arms retain their exact dark poles. No branch
is collapsed or discarded. Brackets use the original factored source
coefficients and DD summation, remain short relative to each source rate,
and retain separated rates. Colliding rounded reciprocals of distinct
source times or an unresolved contraction/bracket select the dense fallback;
unrepresentable physical poles refuse. Work is quadratic in arm count with
linear storage and bounded 64-step iterations. Inline G0/G1/G3 scalar tests
cover complete manufactured finite-loss spectra, original characteristic
residuals, 200-decade rate separation, relative decay through 1e-300, exact
multiplicities, coordinate scaling, underflow recovery and cancellation.
The public-owner regression uses one QR sweep to guard this route's wiring.

Other scalar regimes and coupled mechanical systems retain the dense solver.
The default QR cap is 120 sweeps per active window;
`_with_cancel` accepts an explicit cap and cancellation, returning no partial
spectrum on refusal. Invalid source data and failed mass factors also refuse.
Exact source structure can bound zero multiplicity for triangular invertible
K and certain free scalar systems; excessive exact zeros from QR refuse.
Dense work is cubic in state count; callers bound their source model size.
These are numerical pole values, not eigenvectors, certified eigenvalue
enclosures, Rouché damping certificates, a full-system passivity certificate
or a time-domain realization. Dense QR does not guarantee relative accuracy
of very slow poles or tiny real parts of oscillatory poles across widely
separated rates. The zero-count guards detect only some such failures and
do not certify the remaining decay rates.

`PreparedMaxwellModes` retains the original full M/C/relaxed-K matrices and
Maxwell strain channels for repeated local mechanical mode refinement. Its
`refine` method supplies the complete rational operator
`T(s)=s²M+sC+K+sum(k_j*s*tau_j/(1+s*tau_j)*b_j*b_j^T)` and analytic derivative
to the existing fs-la bordered Newton owner. Binary congruence/rate scaling
and a cached mass Cholesky preserve the source units. `MaxwellMode` publishes
canonical physical displacement with `q^H M q=1`, velocity `s*q`, and every
original branch stress. Independent original-unit mechanical and branch-history
residuals gate publication using factored source products and DD accumulation.
Shared `PoleOptions` ceilings include cached branch storage and preparation/
refinement work; cancellation returns no partial mode.

This API refines one nearby simple mechanical right mode. Every rational
material singularity `s=-1/tau_j` is excluded, and nonzero mechanical
displacement is required. Pure material dark modes therefore remain outside
this API; the complete state-spectrum owner retains their pole values.
Multiple/defective modes can refuse. No arbitrary-seed convergence, spectral
completeness, relative weak-decay accuracy, eigenvalue/shape enclosure or
time-domain realization is claimed. Eight focused
actual-source native checks cover analytic material modes and stresses, coupled
gyroscopic shapes, full mutual mass, extreme binary source scaling, a nearby
weak decay, active-viscosity growth, and admission/resource/cancellation
refusals. These checks do not establish a full workspace Cargo/DSR gate.

`PreparedMaxwellModes::refine_residue` resolves a mechanical right mode and its
bilinear transpose mode through the same bordered Newton owner. The transpose
source equations are independently checked at the right pole. Original-unit
factored products and DD sums retain `ell^T T'(s) q`; authored lower bounds on
equilibrated derivative alignment and source cancellation reject unresolved
normalization. Original resource ceilings admit both solves and retained
responses together. `force_residue` applies `q*(ell^T f)/(ell^T T' q)` to an
original complex force and returns displacement, velocity and every branch
stress residue with linear work/storage and no further eigensolve.
`pole_derivative` applies the implicit simple-pole derivative to full M/C/K,
branch stiffness, log relaxation time and strain directions without
differentiating solver iterations. Gyroscopic and active coefficients remain.
These are local Laurent coefficients and infinitesimal derivatives, without a
complete forced response, complete expansion, spectral enclosure, time
realization or finite-parameter-step bound. Material singularities, dark modes
and unresolved multiplicity still refuse. Eight new focused source checks cover
independent partial fractions/adjugate response, scaling, weak decay,
normalization, cancellation and analytic/finite-difference derivatives; the
eight existing mode tests also pass.

### `harmonic::resonance`: actual driven peaks, notches and linewidths

`refine_transfer_extremum` accepts the actual complex transfer `H(omega)` and
its complete frequency derivative. `HarmonicFrequencyResponse` supplies fixed
real displacement/velocity receiver projections, retaining signed modal
interference and the full existing forced-response tangent. The stationary
condition is `Re(conj(H) H_omega) = 0`; strictly oriented endpoint slopes
distinguish a local maximum from a minimum. Source phase, off-diagonal
damping and frequency-dependent loading are not replaced by pole locations.

`refine_half_power_bandwidth` refines separately bracketed rising/falling
crossings of half the same peak's `|H|^2`. Its reported quality factor is the
actual peak frequency divided by that measured width, not an inferred modal
damping ratio or a statement about total absorbed mechanical power.
Normalized signs avoid overflowing/underflowing amplitude squares. Explicit
absolute/relative frequency resolution and one callback-evaluation ceiling
bound each query; refinement also stops at irreducible floating-point spacing,
visible in the returned brackets. Final half-power crossing brackets must be
strictly separated: unresolved positive linewidth refuses instead of reporting
Q from touching flanks. Exactly stationary samples are reused while their
strict sign brackets shrink, avoiding repeated physical solves without changing
the result or relaxing the resolution. Cancellation or a physical evaluator's
refusal publishes no partial result and preserves the original error.

Complex receiver transfer rows are first-class observations. Fixed complex
weights retain coherent phase between coordinates; optional per-radian weight
derivatives enter the response tangent through the complete chain rule
`r' x + r x'`. This supports phase-weighted sensor combinations, propagation
delays and linear receiver/filter laws without changing the mechanical solve.
The receiver law, its frequency support and its derivative remain caller-owned.
The original real-receiver methods and arithmetic remain unchanged.

These are local numerical queries. The caller isolates one continuous
extremum/crossing per bracket and retains the same physical source, receiver
and model. Opposite signs do not certify uniqueness, spectral completeness,
the nearest crossing, source forward error or absence of an unresolved pole.
Zero-slope endpoints and flat/unresolved interiors refuse. No damping floor,
Lorentzian fit or new eigensolver is introduced. G1/G3 tests cover
displacement/velocity/pole-frequency differences, high-Q linewidths, exact
two-mode antiresonance, full nonproportional damping, extreme response units,
singular evaluator refusal, budgets and cancellation.

`survey_transfer_extrema` and its callback counterpart
`survey_transfer_extrema_with_cancel` find both maxima and minima resolved by
an explicit nonnegative increasing frequency grid. They retain exact transfer
zeros, stationary grid samples, signed modal interference and the source's
complete frequency tangent. Every candidate uses the existing local refinement
owner. All input samples are cached, so shared endpoints and stationary grid
points do not repeat a physical solve; the grid and every refinement share one
evaluation cap. A separate extremum cap refuses truncation before refinement.
Cancellation or any physical/refinement refusal publishes no partial list.
The grid must separate the features the caller needs: an unresolved peak/notch
pair may lie between equal-sign samples, and opposite signs alone do not prove
uniqueness. Band endpoints are excluded, a silent observation returns no extrema,
and no linewidth or free-mode completeness is inferred. Storage is linear in
the input grid and number of returned extrema. `tests/harmonic_resonance_survey.rs`
checks a coupled damped two-mass response with both resonances and its exact
antiresonance, cached samples, and atomic budget/source/cancellation refusals.

### `harmonic`: bounded residual correction at cancellation notches

A coupled LU or compact solve that exceeds the caller's componentwise backward-error
limit may take at most two residual corrections with the same admitted factor.
Each correction uses the original dynamic operator and binary-normalized
double-double residual products, then rechecks the original residual gate before
publication. This recovers weak physical responses near cancellation notches
without changing damping, the source, the operator or the acceptance threshold.
Scalar-only owners retain their direct complex divisions and existing linear
work caps.

Correction scratch is two length-n complex vectors, within the existing payload
allowances. The coupled frequency, curvature and parameter work scores reserve
three times their documented baseline score for the bounded extra solves and
audits. Scalar-only owners retain their original work scores, including the
linear diagonal-bank limits. Cancellation or failure after the final pass still
returns no response and leaves the prepared owner reusable. Focused tests cover
the passive two-mass notch, its frequency tangent, coherent source phase and
extreme force units. The weak-coordinate forward comparison follows the authored
componentwise residual gate; that gate is not a blanket forward-error bound.

### `harmonic`: native frequency/parameter jets for resonance design

`solve_with_frequency_parameter_derivative` supplies the physical displacement and velocity jets needed by the shared peak/notch sensitivity query. It retains the full signed M/C/K direction and the authored F, F_omega, F_omegaomega, F_theta and F_omega_theta source data. The mixed equation is `A*q_omega_theta = F_omega_theta - A_omega_theta*q - A_theta*q_omega - A_omega*q_theta`. The original factor, binary coordinate/time scales and original-equation residual owner are reused; numerical scaling is not differentiated and no perturbed model, modal matching or new factor is introduced. Velocity retains `-i*(q_theta+omega*q_omega_theta)`, including its DC term.

`solve_with_loaded_frequency_parameter_derivative` extends the same owner to an attached frequency-dispersive impedance. Both complete Z_omega and Z_omegaomega arrays are mandatory. The shared normalized frequency action retains `i*B*Z*B^T*q_theta + i*omega*B*Z_omega*B^T*q_theta`, including mutual and nonreciprocal ports and compact factors at physically regularized bare poles. B and the entire Z(omega) law are held fixed with respect to theta: this entry does not claim changing attachment geometry, impedance-parameter derivatives, moving coordinates or frequency-dependent M/C/K directions. Explicit zero derivative arrays declare a constant load law; unsupported or malformed data refuses before solving.

`HarmonicFrequencyParameterResponse` retains the original frequency-curvature and parameter responses, mixed motion and its residual diagnostic. Borrowed displacement/velocity field-jet adapters and fixed real/complex transfer projections connect these directly to the shared resonance query; complete receiver jets can include authored observation changes. Added numerical payload is `8*n*p+3072*n+1536*p+128*k`; the baseline work score is `1024*n*n+2048*n+512*(n*p+p*p+k*k)+2048*(p+k)`, with `512*n*n*p` added for compact loading. Coupled LU/compact owners reserve three times that score for bounded residual corrections; direct scalar owners retain the baseline. The factor and borrowed arrays are excluded.

Focused native tests compare independent scalar quotient and Maxwell cubic derivatives, actual added mass/stiffness/damping, full gyroscopic and nonproportional directions, extreme coordinate units, DC, and compact mutual/nonreciprocal loading. End-to-end tests differentiate the distinct displacement and velocity peaks and the actual two-mass cancellation notch. These are local derivatives of the authored physical model, with no modal truncation, broadband passivity, finite-parameter-step or global peak selection guarantee.

### `harmonic`: changing attachments and dispersive material laws

`solve_with_frequency_parameter_and_port_derivative` extends the native mixed
jet to complete authored B_theta, Z_theta and Z_omega_theta, alongside explicit
Z_omega and Z_omegaomega. B and B_theta are frequency independent; original
coordinate/port dimensions and the physical coordinate basis remain fixed.
M/C/K directions retain the existing frequency-independent contract. The
derivative of a passive material need not itself be passive.

With L=B*Z*B^T, the additional mixed source is
`-L_theta*v_omega-L_omega_theta*v`. Both attachment occurrences are retained:
`L_theta=B_theta*Z*B^T+B*Z_theta*B^T+B*Z*B_theta^T` and
`L_omega_theta=B_theta*Z_omega*B^T+B*Z_omega_theta*B^T+B*Z_omega*B_theta^T`.
The existing complete `-A_omega*q_theta` action is included once. Binary-normalized
factored contractions take O(n*p+p*p) work and O(p) scratch, reusing the admitted
factor without dense load derivatives, perturbed models or finite differences.
A changing load adds 384*p numerical bytes and
`512*(n*p+p*p)+512*p` to the native baseline work score; coupled LU/compact
owners reserve three times the complete score. All authored directions are
admitted before the first solve/allocation.

`HarmonicPortFrequencyParameterResponse` publishes original-unit flow and effort
mixed derivatives, including `B_theta^T*v_omega` and all four terms of
`(Z*u)_omega_theta`. Its signed boundary-work derivative retains both distinct
cross contractions for mutual/nonreciprocal ports. The parent response also
reports independently differentiated input work, full symmetric-C viscous work
and their balance against port work. Symmetrizing C and C_theta before contraction
prevents changing gyroscopic coupling from creating direct heat. Normalized
double-double products retain finite cancellation before overflow; boundary
work is not reclassified as constitutive heat or storage. Fixed-load and unloaded
mixed entry points receive the same complete physical work accounting.

Focused regressions compare an independently assembled dense mass/damper/spring
pencil and the original Maxwell-memory rational quotient, with changing
attachments, material relaxation, all source jets, mutual/nonreciprocal loads,
DC and resonance. Port-unit changes by 2^+/-400 preserve mechanical motion and
work; individually overflowing reactive mixed-work terms cancel to finite zero.
End-to-end queries differentiate actual peak locations, both asymmetric
half-power crossings, bandwidth and Q for attachment, load mass, loss and
stiffness directions. Exact resource limits, invalid data, cancellation and
reusable-owner refusal are checked. A separate admission regression verifies
that coupled first-frequency solves actually reserve the documented three-pass
work allowance while scalar owners retain their direct baseline.

### `harmonic::resonance`: physical port-flow and effort observations

Native frequency responses expose fixed real/complex port-flow and port-effort
transfer rows, plus complete authored receiver frequency tangents. Native mixed
responses expose borrowed flow/effort field jets and fixed real/complex parameter
projections. These retain the physical load owner's original conjugate variables,
port ordering and units; coherent weights act as linear transfer rows without
conjugation. An absent load or incomplete physical jet refuses explicitly.

The effort tangent includes the complete dispersive Z_omega action. Changing
attachments and material laws remain in all five field derivatives, so the
existing peak/notch and actual half-power sensitivity queries apply directly to
transmitted-force or flow observations. Authored receiver jets retain changes in
the observation map, including the B_theta*effort term when observing a generalized
load force. Observation does not modify the mechanical load or power account.

First-derivative and complete parameter projections use the same existing
binary-normalized double-double product owner. Individually overflowing coherent
contributions may cancel before physical recovery. This avoids premature
overflow and subnormal rounding within the shared accumulator's working
precision. Projection takes linear work, constant scratch and no allocation; derivative shape and finiteness checks precede contraction.

Independent tests distinguish a damper-loaded oscillator's transmitted-force
peak from its displacement peak and verify both actual half-power crossings,
bandwidth, Q and their physical loss derivatives. A dispersive Maxwell load
with B(theta)=exp(theta)*B and Z(omega,theta)=exp(-2*theta)*Z(omega) checks that
changing coherent receiver weights preserve physical force, velocity and
resonance under a mere port-coordinate change. Complex receiver phase tangents,
missing jets, invalid weights and extreme coherent cancellations are covered.

### `harmonic::resonance`: signed physical-power extrema and sensitivities

`refine_power_extremum` and its complete parameter-derivative entry root the
actual signed P_omega through the existing numerical query owner. Native
frequency, curvature and mixed response adapters expose original input power,
structural viscous power and net load-port boundary work. Units and signs are
preserved; absent or incomplete load-port observations refuse.

Refinement evaluates only P and P_omega. One final curvature or mixed jet is
reserved inside the same evaluation cap, reproduces both sampled values within
32 machine epsilons without a unity floor, and supplies nonzero P_omegaomega
with the strict maximum/minimum orientation. This rejects a divergent
positive/negative power pole's reversed derivative without fabricating a finite
extremum. Exact stationary samples remain cached by the shared query owner.
Zero and negative stationary powers are admissible; the query never squares
power or constructs an artificial amplitude.

The local frequency derivative is `-P_omega_theta/P_omegaomega`. Total observed
power sensitivity retains `P_theta+P_omega*omega_theta` at the actual returned
frequency, using the existing normalized signed-work accumulator before
physical recovery. All source, attachment and constitutive directions belong
to the same authored jet. No new factor, perturbed solve, dense derivative
operator or numerical root engine is introduced.

Native regressions independently verify loaded mass/damper/spring source,
viscous and absorbed-power peak frequencies, heights and their complete physical
sensitivities. Signed/zero extrema, power units spanning 2^+/-700, exact stationary
sample reuse, finite-resolution total derivatives, sampled/unsampled pole
refusals, degenerate curvature, inconsistent jets, total budgets and final-call
cancellation are checked. These APIs resolve individual simple extrema of one
continuous physical branch; they do not infer power linewidths, global feature
completeness, evaluator forward-error bounds or constitutive heat/storage splits.

### `harmonic::resonance`: actual physical-power bandwidth and Q

`refine_power_bandwidth` and its parameter-derivative entry resolve both actual
half-height crossings of one positive signed-power maximum. The numerical
level is P-P_peak/2, normalized before subtraction; the observed work channel,
material law, source and units remain those of the peak. The caller separately
isolates each continuous flank. Neither a symmetric/Lorentzian line shape nor
an amplitude-square proxy is introduced.

The existing shared sign-search owner retains strict crossing orientation,
original evaluator errors and cancellation. Both final crossing brackets must
be disjoint, and both original power slopes must be nonzero and oriented toward
the same peak. A shared peak endpoint is evaluated once. One cap covers both
flanks, and the value-only entry requires no final curvature/parameter solve.
The returned Q is the actual peak frequency divided by measured width.

Parameter queries reserve exactly two final edge jets before refinement and
include them in the same total evaluation ceiling. Each reproduces its original
P/P_omega sample without a unity floor. Edge continuation uses
`omega_theta=(P_peak_total_theta/2-P_theta)/P_omega`, retaining motion of the
physical peak height. Width and Q then retain both edge and peak motion.
Existing binary-normalized double-double contractions form all differences and
ratios; no perturbed models, finite differences or new numerical root owner
are introduced. Only P, P_omega and P_theta are consumed at an edge.

Independent native tests verify input, viscous and port-power crossings and Q
against original loaded-oscillator formulas, including complete mechanics,
attachment, material and source directions. A Maxwell-memory model is checked
against its original cleared real cubic, an exact sqrt(2) power maximum and
independent 90-digit asymmetric half-height roots; six physical directions
include material relaxation. Its flow-amplitude slope is nonzero at the power
maximum. Source gain/phase invariance, power scaling through 2^+/-700 and paired
power/frequency scaling through 2^+/-600 exercise the normalized ratios.
Unresolved narrow widths, malformed or inconsistent data, exact evaluation
budgets and cancellation refuse atomically. This remains local continuation
of one selected positive peak and its isolated crossings; it does not infer
a free-pole decay rate, constitutive loss factor or global feature completeness.

### `harmonic::resonance`: coupled signed-power surveys

`survey_power_extrema` and its complete parameter-derivative entry find the
signed maxima and minima resolved by an authored increasing frequency grid.
They cache each original P/P_omega grid sample once. A shared generic grid
scanner also serves the existing complex-transfer survey; strict slope
reversals and their intervening stationary samples identify candidates without
double counting. Band endpoints lack two-sided evidence and are excluded.

Every candidate retains the existing numerical sign-search owner and original
power resolver. One final curvature or complete mixed jet per candidate is
reserved after discovery and before refinement. The global evaluation ceiling
includes the grid, all refinements and every final jet. A separate positive
extremum cap refuses truncation before refinement. Original evaluator errors,
cancellation, unresolved stationary intervals, inconsistent final observations
and wrongly oriented/degenerate curvature publish no partial collection.

`PowerParameterSurvey` reports each original extremum's local frequency and
total observed-power sensitivities, retaining the full source, material and
attachment direction. Only first derivatives are requested during refinement;
the complete jet is evaluated once at each resolved extremum. Signed and zero
powers retain their meaning. Temporary candidates store index pairs into one
ordered sample cache, so storage is O(grid+extrema), with constant refinement
scratch. No new numerical root engine, physical factor or dense derivative
operator is introduced.

Native regressions compare two coupled absorbed-power peaks and the intervening
minimum against an independent normal-coordinate rational response, including
a complete physical scaling direction. Signed/zero grid extrema, 2^+/-700 power
units, exact sample caching and total callback limits, candidate caps, source
errors, final-jet/cancellation refusals, degenerate curvature and poles are
checked. A deliberately insufficient grid demonstrates the stated resolution
boundary: these surveys do not certify global spectral completeness, assign
free modes, or differentiate peak switching. Net boundary work keeps its
authored physical interpretation.

### `harmonic::resonance`: complete observed-transfer sensitivity surveys

`survey_transfer_extrema_with_parameter_derivative` returns every simple
complex-response peak/notch resolved by an authored increasing frequency grid,
together with its complete local physical design derivative. It retains the
same source, receiver, material and attachment parameter across both callbacks.
The shared grid owner caches original H/H_omega observations once, identifies
strict sign reversals and retains stationary grid samples between oriented
neighbors. It reserves one final complete parameter jet per candidate before
refinement, within one ceiling covering the grid, all refinements and all jets.

Each candidate uses the original transfer refinement and parameter resolver.
Finite maxima retain their existing pole-resolution guard, and every final
jet reproduces the original observation and passes signed-curvature and
representability checks. A simple exact transfer zero retains its finite
notch-frequency and total complex derivatives; log-amplitude and phase
derivatives remain explicitly absent there. Full jets are evaluated only at
resolved extrema, never at every search trial.

`TransferParameterSurvey` reports frequency-ordered original extrema with
their frequency, total complex, log-amplitude and phase derivatives and existing
conditioning diagnostics. Each entry counts newly performed local first calls
plus its one final jet; shared grid calls belong to the survey total. Storage
is O(grid+extrema), with no dense jet cache, additional physical factor,
perturbed models or second numerical root engine. Extremum/evaluation caps,
original source errors and cancellation publish no partial list.

Native regressions compare both peaks and an exact omega=2 antiresonance of
a coupled two-mass pencil against independent original-quotient and polynomial
derivatives. Physical spring changes, complete frequency dilation and coherent
source gain/phase directions retain each feature's correct continuation.
Exact cancellation has zero observed amplitude and absent log/phase tangents,
while its notch-frequency derivative remains available. Cached grid reuse,
extreme coherent units, global evaluation reservations and atomic refusal
exercise the complete survey boundary.

The grid must separate every feature the caller needs: equal-sign endpoints
can hide a peak/notch pair, and endpoints lack two-sided evidence. These are
local derivatives of the selected observed features, not modal matching,
finite-parameter-step prediction, continuation through degeneracy or a global
spectral completeness certificate.

### `harmonic::resonance`: physical peak, notch and linewidth sensitivity

`refine_transfer_extremum_with_parameter_derivative` and its callback entry retain the original observed response and first-derivative sign search, then request one complete H/H_omega/H_omegaomega/H_theta/H_omega_theta jet at the resolved frequency. The implicit derivative of `Re(conj(H)*H_omega)=0` differentiates either a strict peak or a strict notch, including a simple exact transfer zero. Total complex response, log-magnitude and phase derivatives follow the moving extremum; log-magnitude and phase are explicitly absent at an exact zero. Signed curvature and stationarity diagnostics expose the remaining local conditioning and search resolution. Zero or wrongly oriented curvature refuses without an arbitrary curvature floor.

`refine_half_power_bandwidth_with_parameter_derivative` differentiates the two separately bracketed physical half-height crossings, their actual width and the observed response Q. Its threshold includes motion of the peak height, so changing source gain does not manufacture a linewidth change. The one or two final derivative calls are reserved before each search and included in its evaluation ceiling. Their H/H_omega samples must reproduce the search observations within 32 machine epsilons relative to their largest quadrature, without a unity floor. Full jets are never requested at every search trial. Original evaluator errors and cancellation remain atomic.

`HarmonicTransferFieldParameterJet` and `HarmonicReceiverParameterJet` project complete coherent physical observations, retaining all source/receiver frequency, parameter and mixed product-rule terms. Real and complex fixed receivers use the same owner. Binary-normalized double-double contractions keep signed products scaled until cancellation; observation uses constant scratch and linear work without another physical solve. These derivatives are local smooth-branch responses, with no finite-step, peak-switching, forward-error or global completeness claim. Focused callback and receiver tests cover mechanical M/C/K peak shifts, exact two-mass cancellation, physical propagation delay, asymmetric widths/Q, gain invariance, extreme units, coherent overflow cancellation and bounded refusal semantics.

### `harmonic`: full forced response and resonance

`harmonic::PreparedHarmonicInterference` prepares simultaneous coherent forces
on one original real M/C/K pencil. Each distinct frequency uses the existing
harmonic response owner; exactly repeated frequencies combine forces in DD
before their solve. Nearby tones retain their separate frequencies and beats.
The full coupled mass, stiffness and velocity operators remain intact, including
nonproportional damping and signed gyroscopic terms. No diagonal modal projection
or artificial loss is introduced. The source clock uses `exp(-i omega t)` and
peak force amplitudes; DC contributes its actual real static displacement.

The prepared observer applies M, K and D=(C+C^T)/2 to the original phasors once.
An observation through `fs-fft::harmonic_window::HarmonicWindow` then integrates
kinetic/elastic energy, signed input power, viscous power and the independently
formed storage rate, retaining all cross-frequency and sum-frequency products.
`F^T v - v^T D v - (v^T M a + v^T K q)` is reported as a balance residual.
A uniform window's integrated storage rate is its endpoint energy difference;
a triangular detector instead measures the weighted local rate. Powers from
separately cycle-averaged tones are never substituted for the coherent result.

Preparation admits at most 128 authored tones plus coordinate, peak-payload,
response-residual and product-work ceilings. Exactly equal frequencies share one
solve. Retained arrays use 144*n scalar bytes per distinct tone plus headers;
the peak budget additionally covers one existing dense response preparation and
bounded scratch. Warm observation allocates nothing and requires six scalar
harmonic products per coordinate and ordered frequency pair, with one shared
analytic kernel per pair and no matrix applications or solves. DD contractions
preserve cancellation residuals until final output, without claiming relative
accuracy at arbitrary nulls. Cancellation publishes no partial observation.

This is the particular harmonic motion of the supplied finite-dimensional
pencil, not a startup/free-decay solution, stochastic coherence inference,
nonlinear mixing model, modal truncation certificate or broadband impedance.
Singular frequencies retain the existing refusal. G1/G3 native tests in
`tests/harmonic_interference.rs` compare against independent scalar/2x2 inverse
solutions and real-signal quadrature; exercise coupled energy, gyroscopic work,
short shutters, beating, DC, exact duplicate cancellation and source-clock
shifts; and verify count, payload, product and singular-response refusals.

`PreparedHarmonicInterference::try_with_port_impedance` and the corresponding
ensemble constructor share the original loaded-response owner. A borrowed
`HarmonicInterferenceLoad` supplies the original real force map B and exactly
one `HarmonicImpedanceSample` per distinct source frequency, in strictly
increasing order. The full mutual, reactive impedance is retained; there is no
interpolation or inferred broadband causality certificate. One loaded
factorization per frequency serves all mean/factor right-hand sides. The shared
window observer independently integrates `flow^T effort`, including all
within-bank cross-frequency products, as `mean_boundary_power_w`. This signed
net boundary work can be negative when a passive load returns stored energy;
it is not automatically heat or radiation. The balance subtracts structural
viscous power, boundary work and independently differentiated structural
energy. Warm work becomes `(6*n + ports)` products per ordered tone pair;
admission also covers 32*ports retained flow/effort bytes per bank-frequency
and one existing loaded preparation. Native regressions use an independent
passive mutual mass/spring/resistor network and exact finite real/proper
ensembles integrated over uniform and triangular windows.

`harmonic::PreparedHarmonicResponse` prepares the complete retained mechanical
pencil `K - omega^2 M - i omega C` under `exp(-i omega t)`. The source mass,
stiffness and viscous damping remain fully coupled; a mode count or a list of
free-decay poles does not substitute for the forced port response. `try_new`
accepts real row-major matrices with exactly symmetric mass and stiffness.
The velocity operator may contain symmetric viscous damping and skew
gyroscopic coupling. `try_modal` supplies the mass-normalized diagonal
stiffness convenience while retaining both velocity contributions;
factored frequency detuning preserves narrow peaks and exact nominal resonance.

An unloaded `try_modal` pencil with exactly zero off-diagonal damping uses
compact diagonal storage. Every nonzero coupling, however small, and every
port-loaded pencil retains the dense path. Layout selection uses the source
damping, including at DC where its off-diagonals affect the frequency tangent.
The complete caller-owned velocity matrix still receives finite-entry
admission and requires `O(n*n)` inspection. Its symmetric viscous part must
pass the diagonal or dense PSD admission below. Subsequent compact assembly and
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
`P_loss = Re(v^H D v)/2`, where `D=(C+C^T)/2`. Source units and conjugate
force/velocity projection remain the caller's responsibility. Physical radiation
loss appears here only if its loading is actually included in the supplied pencil.

Mass must pass positive-definite admission. The symmetric viscous part and
stiffness receive nonnegative-coefficient checks when exactly diagonal.
Coupled symmetric parts pass a Jacobi PSD check after their own diagonal
congruence, so a small indefinite subspace cannot hide behind a larger
unrelated mode. The numerical allowance is
`64 n epsilon`; this is not an interval passivity proof. Source coefficients
are never clipped, diagonalized or shifted. Zero damping and stiffness
nullspaces remain intact; singular lossless resonances and unconstrained DC
free motion refuse instead of manufacturing finite amplitudes.

Gyroscopic coupling is retained with its original sign in the full dynamic
operator and its analytic frequency tangent; it can split circular resonances
and break reciprocity without introducing dissipation. PSD admission and work
diagnostics take the symmetric part before accumulation, so large gyroscopic
exchange cannot manufacture heat or hide an active viscous direction. This
uses no extra dense matrix or factorization. Pure skew operators with zero
viscous diagonals are admissible; finite source coefficients remain required.
G1/G3 regressions cover both signed circular branches through their split
resonances, lossless nonreciprocal transfer, analytic motion/heat frequency
tangents and refusal of indefinite symmetric parts.

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

`try_diagonal_modal_with_port_impedance` combines vector-authored independent
modes with the same complete physical load and Hermitian-PSD admission. It
keeps the structural operator diagonal and stores scaled B and BZ factors.
No mode or coupling is projected away; reactive, mutual, nonreciprocal,
rank-deficient and lossless impedances remain. The bounded complex LU acts on
ports plus retained pivot modes, not the entire modal bank.

For `A = D + L R^T`, `R = B_hat` and `L = -i*mu*B_hat*Z_hat`, every exact
zero diagonal of D is retained. Near-pole rows also enter the pivot block when
`|D_i| <= min(64*epsilon/residual_limit, 0.25) * h_i`, where
`h_i = mu * sum_a |(B_hat Z_hat)_ia| * sum_j |B_hat_ja|`.
This is an algebraic pivot heuristic, not a damping floor, modal truncation or
forward-error certificate. Eliminated modes N and retained modes P give the
bordered matrix `[I+R_N^T D_N^-1 L_N, -R_P^T; L_P, D_P]`. Its solution and
`q_N = D_N^-1(F_N-L_N*u)` retain all original equations. Neither Z nor a zero
bare diagonal is inverted. More exact bare zeros than ports necessarily
refuse as singular; oversized retained cores refuse their resource bounds.
Independent binary row/column equilibration changes no physical coefficient.
Each port RHS contraction receives its own exponent before the actual
equilibrated core rows are normalized together, so an uncoupled huge modal
response cannot suppress a tiny retained force.

Storage is `O(n*p+k*k)`, preparation `O(n*p*p+k*k*k)`, and a factor solve
`O(n*p+k*k)`, with `k = p + retained_pivot_modes`. Audits deliberately stream
the actual loaded matrix entries in `O(n*n*p)` work: the componentwise
denominator remains `|A||q|+|F|`, not a looser factorized envelope. The same
streamed coefficients enter the independently reconstructed frequency-tangent
audit. Existing original-coordinate flow, effort and physical work owners are
reused. At DC finite velocity loading vanishes from the solve but remains in
the frequency derivative.

Compact load work admission is `n*n*p + 2*n*p*p + 8*n*p + k*k*k`.
Additional load payload admission is
`64*n*p + 256*p*p + 128*p + 144*k*k + 256*k + W_LU(k,1)`, where W is the
native blocked-LU workspace query. The first terms cover physical/scaled port
copies and modal index arrays; the core allowance covers assembly, compensated
sums, equilibration, factor and solve scratch. Core dimension also passes the
existing fs-la bounded dense-work admission. These caps precede core
allocation/factorization. Arithmetic, cancellation and original-row residual
refusals publish no partial response.

`tests/harmonic_compact_ports.rs` covers a 1,025-mode loaded bare-pole closed
form, the full DC load tangent, dense parity for repeated/near bare poles,
mutual nonreciprocal, rank-deficient and lossless loads, physical power, and a
mixed-scale force fixture with a huge uncoupled response. Genuine unsupported
poles, core resource limits and cancellation refuse explicitly.

`solve_with_parameter_derivative` reuses the same prepared factor for signed
physical M/C/K/F directions at fixed frequency, coordinate basis, attachment
map and load impedance. It consumes the existing `QuadraticParameterDirection`
in original row-major coordinates: mass and stiffness directions are symmetric
but need not be positive; the full velocity direction can change both viscous
and gyroscopic coupling. A direction can introduce coupling into a diagonal
base. The differentiated equation is
`A q_theta = F_theta - K_theta q + omega^2 M_theta q - C_theta v`.
No perturbed factorization, eigenvector derivative or numerical-scaling
derivative is used. Binary-normalized products and DD RHS accumulation retain
physical units without first overflowing `omega^2` or matrix/state products.

The output includes original-coordinate motion derivatives and independently
differentiated input, viscous and fixed-load boundary work. In particular,
the viscous derivative includes `v^H D_theta v / 2`, while a changed purely
skew gyroscope contributes no direct heat. A published-coordinate row audit
uses the uncancelled component envelope as well as the ordinary tangent-solve
backward check. Neither is a forward-error certificate at a narrow resonance.
Signed power products retain binary exponents through cancellation before
recovering physical units, including lossless source-phase rotations and
motion in a damping direction's nullspace.

`harmonic_mean_power_derivative_with_cancel` exposes that existing physical
work contraction for original-state consumers. It contracts
`Re(flow_theta^H*effort + flow^H*effort_theta)/2` across equal finite port slices
using the same binary-normalized DD products, retaining mutual cancellation
before rounding a total or restoring physical units. It allocates nothing,
uses constant scratch and linear work, polls before/during/after contraction,
and refuses unrepresentable final values. It neither infers dissipative heat
nor supplies an exact-summation or forward-error guarantee at arbitrary nulls.
The accompanying focused cases in `harmonic_parameter` cover reactive phase
changes, overflowing products with finite mutual cancellation, original port
unit changes, a retained smallest-subnormal result, invalid inputs and late
cancellation. Those new cases remain unexecuted under the dependency resource
blocker.

Additional payload admission is `512*n + 384*p`; work admission is
`256*(n*n+n*p+p*p)+256*(n+p+k)`, with `64*(n*n*p+k*k)` added for a compact
loaded factor of dimension k. Direction/source inputs and the original factor
are excluded. No derivative matrix or second factor is allocated; dense
directions require quadratic work even on a compact base. Cancellation,
malformed directions, unrepresentable arithmetic and exhausted budgets
publish no partial result or changed owner.

Tests in `tests/harmonic_parameter.rs` compare scalar resonance slopes with
closed forms, keep a finite nonzero sensitivity at an exact dark response,
retain signed gyroscopic coupling, compare coupled loaded motion/work against
independent physical perturbations, check coordinate-unit changes and typed
refusals. These are derivatives of the supplied finite-dimensional mechanical
pencil, not constitutive identification, changing physical-basis derivatives,
derivatives of B or Z, material/geometry qualification, or truncation bounds.

`solve_with_parameter_and_port_derivative` additionally accepts the authored
`HarmonicPortParameterDirection`: a complete real n*p derivative of B and a
complete complex p*p derivative of Z in the original conjugate port units.
For `y=B^T v` and `e=Z y`, the differentiated mechanical RHS includes
`-B_theta e - B (Z_theta y + Z B_theta^T v)`. Both occurrences of the physical
attachment map remain, including newly coupled ports and nonreciprocal loads.
Published derivatives are `y_theta=B^T v_theta+B_theta^T v` and
`e_theta=Z y_theta+Z_theta y`; their signed boundary-work derivative is evaluated
independently. Direction matrices need not themselves be passive or reciprocal.

`solve_with_total_parameter_derivative` uses the same engine and also accepts
a finite signed `omega_theta`. Its authored M/C/K/F and optional B/Z directions
are total derivatives along the parameter, including constitutive dispersion.
The explicit frequency action is
`omega_theta*(2*omega*M*q+i*C*q+i*B*Z*B^T*q)`, and velocity is
`-i*(omega*q_theta+omega_theta*q)`, retaining the DC kinematic term. None for the
port direction fixes B and Z but keeps their explicit frequency action. This
does not replace a physical stiffness, damping or impedance slope with a frozen
material law. Existing fixed-frequency methods delegate with `omega_theta=0`;
the existing specialized frequency owner retains its separate fast modal path.

Load actions factor through port vectors in O(n*p+p*p) work and O(p) added
memory, without constructing an n*n load derivative or another factor. Binary
exponents survive intermediate flow/effort contractions; their uncancelled
envelopes enter the original-coordinate tangent audit. A supplied port
direction, or changing frequency with an attached load, adds `192*p` payload
bytes and work score `256*(n*p+p*p)+256*p`. Nonzero `omega_theta` also adds
`192*n*n+256*n` work. These allowances augment the parameter bounds above.
The coordinate basis, port count and physical state interpretation remain fixed.
`tests/harmonic_parameter_ports.rs` covers total scalar/DC formulas, independent
coupled nonreciprocal perturbations, compact loaded bare poles, extreme port
units and port reparameterization, cancelling huge attachment products,
frequency-owner agreement and transactional refusals.

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
A nonzero-frequency compact loaded owner admits work score
`64*(n*n*p+n*p+p*p+k*k)+256*(n+p+k)`, including the streamed audits and bordered
solves. Its added payload remains `8*n*p + 512*n + 384*p`; the retained core
belongs to the prepared owner's load budget. DC uses the diagonal structural
work score plus the existing port-tangent terms.
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

Frequency tangents use the same binary-normalized double-double physical-work
contraction as parameter tangents. Source work, full viscous work and boundary
work retain their signed products until cancellation, before recovering f64
units. No extra factorization or numerical arrays are required. Focused G1/G3
regressions retain zero work for a large purely reactive response with a moving
source phase, and the zero power slope of a finite driven resonance whose
individual differentiated products exceed f64. This protects representable
physical totals; it does not certify relative accuracy at arbitrary nulls.

`solve_with_frequency_curvature` extends that same immutable factor to the
complete second frequency derivative. It solves
`A*q_ww = F_ww - A_ww*q - 2*A_w*q_w`, where
`-A_ww = 2*M + 2*i*B*Z_w*B^T + i*w*B*Z_ww*B^T`. Both physical load
derivatives are mandatory for loaded owners, including explicit zero matrices
for locally zero derivatives. Source slope and curvature are independently
authored. The original M/C/K and conjugate port map remain fixed. No finite
frequency perturbation, additional factor, diagonal approximation or assumed
Lorentzian line shape enters the calculation.

First and second pencil actions share one factored numerical kernel, retaining
full mutual load and signed gyroscopic/nonproportional coupling. Each response
order has its own binary amplitude scale; the second derivative restores two
physical frequency-unit factors. Published motion obeys
`v_ww = -i*(2*q_w+w*q_ww)`, including DC. The independent second-equation audit
reconstructs the published physical curvature and retains separate uncancelled
source/action envelopes. Its residual is a backward check, not a near-pole
forward-accuracy certificate. Original port flow and effort include all three
second-order product-rule terms. Signed source, symmetric viscous and boundary
power curvatures are contracted independently with the existing binary-scaled
double-double owner; `harmonic_mean_power_curvature_with_cancel` exposes that
same contraction to other power-conjugate field owners. This does not infer a
load's material heat or stored energy from impedance samples.

`HarmonicFrequencyCurvatureControls` admits total additional payload
`8*n*p+1536*n+768*p` and work `256*(n*n+n*p+p*p)+1024*(n+p)`. Unloaded
diagonal banks retain linear work `1280*n`; compact loaded factors substitute
`n*n*p` for `n*n` and add `256*k*k+1024*k`. The prepared factor and caller
arrays are excluded, and refusals publish no partial result. Focused
`harmonic_curvature` regressions specify independent scalar quotient rules,
causal spring/inertance and Maxwell dispersion, DC gyroscopic coupling,
extreme finite high-Q power, dense/compact agreement at a loaded bare pole,
linear modal-bank budgets and transactional cancellation/reuse. Pure numerical
kernel tests additionally exercise the original factored actions and physical
work without requiring an executor or substituting a numerical solver.

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

`energy_j(states)` is the total modal mechanical energy of any state vector in
the model's modes, by the per-mode law and compensated sum a step uses for its
own totals, so after a step it equals that step's `total_modal_energy_j` bit
for bit. It admits nothing against the budgets and changes nothing; a wrong
state count is `StateCountMismatch`. It lets an owner outside this crate
(the composition-kernel Component in fs-nlmodal) measure stored energy from a
state alone.

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
`held_force_mean_velocity_relation` also exposes the exact same-interval mean
velocity and its affine response to an additional constant force. It uses the
prepared entire gains and the initial force imbalance rather than subtracting
endpoint displacements. Power-conjugate sampled loads can therefore close
their direct reaction against resolved motion, including small motion on a
large static equilibrium, without another transition or per-sample allocation.
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

`HarmonicForceTransition::prepare_sensitivities` adds analytic transient
partials with respect to independent mass-normalized stiffness `k = omega^2`,
physical viscous rate `c`, forcing angular frequency `Omega`, and common initial
source phase `phi`. The mechanics derivatives hold the forcing clock and phase
fixed; source derivatives hold the mechanics and interval-start state fixed.
Natural-frequency `omega` and damping-ratio changes require the k/c chain rules
documented for held-ramp forcing below; they are independent of driving-frequency
changes.

The k/c path differentiates the original continuous-tone convolution and its
exact doubling, sharing the homogeneous-flow derivative owner with held-ramp
forcing. The original near-DC force imbalance preserves static balance, and
independent homogeneous retentions preserve late damping-sensitivity tails.
For `F(t) = Fc*cos(Omega*t+phi) + Fs*sin(Omega*t+phi)` at `phi=0`, the source
phase derivative is `Fs*cos(Omega*t)-Fc*sin(Omega*t)`. Its response reuses the
original constant-tone owner with rotated quadratures `(Fs,-Fc)`. The frequency
derivative multiplies that same source by physical time `t`; it reuses the
existing entire affine-envelope convolution with the original unwrapped
`Omega*h`. No detuning denominator, finite difference or substituted
steady-state response is introduced. Startup, exact resonant secular growth,
free coordinates and critical/overdamped motion retain the original dynamics.

The original motion is returned unchanged. Cold preparation caches the added
affine gains; warm `advance` and `pullback` allocate nothing and call no
transcendental functions. The pullback returns initial-state, cosine/sine-force,
k/c, forcing-frequency and initial-phase cotangents. Initial phase means a
common additional rotation of the prescribed force quadratures at the interval
start, independent of `Omega`. At DC, a sine quadrature contributes no primal
force but can have a nonzero clock or phase derivative; it is retained.

`HarmonicForceTransition::pullback_force_phase` supplies the transpose and clock
derivative of the SAME source rotation as `advance_force_phase`. A reverse
trajectory must combine its source-force cotangents with the motion pullback
and sum both forcing-frequency contributions at every interval. The rotation's
frequency contribution is `h*(next_sine*bar_next_cosine -
next_cosine*bar_next_sine)`. A trajectory's common initial-phase derivative is
formed once from its original force covector, `Fs*bar_Fc-Fc*bar_Fs`; additionally
summing every local phase partial would double-count that dependence.

These are fixed-duration, fixed-coordinate and fixed-mass-normalization
partials. Material/basis changes, work/heat and nonsmooth observations retain
their separate derivatives. Finite arithmetic can lose accuracy; unrepresentable
results refuse without changing motion. G1/G3 checks in
`tests/harmonic_sensitivities.rs` retain the independent mechanical variational
ODEs, secular/DC limits, static balance, late critical decay and reverse sweeps.
`tests/harmonic_clock_sensitivity.rs` independently integrates the physical
frequency/phase variational equations across damping regimes, checks exact DC
and resonant limits, and compares an unequal-interval state-and-source reverse
sweep with the complete interval and independent ODE. Phase-transpose and
finite-input refusal checks cover reusable clock propagation.

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

`PreparedHarmonicForce::prepare_power_sensitivities` cold-prepares a borrowed
`PreparedHarmonicPowerSensitivity`. Its `gradients` returns separate
`HarmonicPowerGradient` entries for signed input work and viscous dissipation,
with physical initial-state, cosine/sine force, stiffness k, damping rate c,
forcing-frequency and common-phase derivatives. The original forward cache,
heat square root and independently integrated work form remain unchanged.
The derivative follows the original resolved Gaussian quadrature and exact
block composition, including every term in the derivative of `Q + E^T Q E`.
The balanced force-imbalance chart retains its fixed conditioning scale during
differentiation; its physical k chain includes the initial-coordinate change.
Heat differentiates `c*v^2` directly, so its right derivative at zero damping
retains the integral of squared velocity even though the heat itself is zero.
Neither power derivative is recovered from an endpoint-energy difference.

Queries allocate nothing and evaluate no transcendental functions. Reverse
trajectories add local power covectors to motion cotangents and propagate the
shared source phase through `pullback_force_phase`. Duration, changing basis
and mass normalization remain fixed. This is one coherent constant-envelope
tone: independently differentiated single-tone heats do not replace simultaneous
tone cross terms, and affine-envelope power derivatives are not supplied here.
The derivatives are finite-arithmetic estimates, not rigorous error enclosures;
nonfinite or unrepresentable queries refuse without mutating the preparation.
G1/G3 checks in `tests/harmonic_power_sensitivities.rs` compare all physical
partials with independently integrated variational ODEs, verify zero-loss/DC
limits, differentiated energy balance, unequal-step source/state reverse
chains, unit rescaling and reusable failure behavior.

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

`PreparedModalTones::prepare_sensitivities` cold-prepares optional exact
constant-tone or affine-envelope motion derivatives for every mode and distinct
nominal frequency.
The returned `PreparedModalToneSensitivity` borrows the unchanged bank and
reuses the original scalar harmonic sensitivity caches. `advance_mode` returns
the original coherent endpoint, independent physical k/c partials, and each
original source's forcing-frequency/initial-phase partials. `pullback_mode`
returns initial-state, k/c and original source-force/clock cotangents for an
endpoint objective. Fixed-capacity results expose only the complete original
caller order through `source_partials` or `source_cotangents`; warm queries
allocate nothing and evaluate no transcendental functions.

The shared scalar mechanical contraction receives one compensated aggregate
held-force imbalance and the summed source derivative gains. It differentiates
one physical initial state, preserving a nearly balanced DC load's tiny damping
derivative instead of subtracting separately differentiated trajectories.
Duplicate frequencies coalesce their primal forces and shared k/c response,
while each original input tone retains its own frequency and phase parameters.
Its clock partial therefore uses its own original force quadratures, even when
a duplicate cancels the nominal motion completely. Force cotangents map back
to every original source; nearby distinct tones retain their beat dynamics.

`PreparedModalTones::pullback_force_phase` delegates the original scalar phase
rotation's transpose and frequency derivative for the selected original tone.
A reverse trajectory combines source-state and mechanical cotangents at every
interval. Shared physical parameters require their authored chain rule across
tones and modes: a common fundamental with `Omega_j = n_j*Omega_0` weights each
source-frequency cotangent by `n_j`, and a shared initial harmonic phase uses
the corresponding phase factors. No frequency-ratio or parameter identity is
inferred from coincident nominal frequencies.

`ModalToneSensitivityBudget` separately caps the additional cache/scratch
payload, scalar preparation work and one mode's forward/reverse query work.
Admission precedes allocation; cancellation is polled during admission, between
scalar preparations and query contributions, and before publication. A refusal
returns no partial derivative and changes no bank/model state. Pure mode queries
leave physical state and observation admission to the caller, as do the original
mode previews. Constant and affine preparation/query laws must agree. Work/heat,
changing duration, basis/mass normalization and nonsmooth objective derivatives
remain outside these motion partials; the separate scalar power API above owns
constant-envelope work and heat gradients.

G1/G3 checks in `tests/tone_sensitivities.rs` compare independent physical
variational ODEs across free, resonant, critical, overdamped and near-frequency
motion; retain nonzero individual gradients at destructive duplicate
cancellation; resolve tiny near-DC equilibrium damping derivatives relatively;
and compare a multimode unequal-interval source/state reverse sweep with the
complete interval and independent ODEs. Transpose, budget, affine-law and
cancellation refusals cover reusable bank queries.

For affine banks, `advance_mode_affine` retains the original coherent endpoint
and returns k/c plus caller-ordered source-frequency/common-phase partials.
`pullback_mode_affine` returns initial-state, both envelope-endpoint pairs,
physical k/c and independent source-clock/phase cotangents. Both endpoint pairs
use the interval-start carrier basis. Exact duplicates share nominal motion
and mechanical tangents but retain independent source columns, including
complete destructive cancellation.

The affine owner differentiates the existing entire linear-envelope action
and its doubling. Source-frequency differentiation integrates the required
quadratic envelope. Separately doubled complementary Bernstein moments retain
small start-envelope and clock responses after long decay without subtracting
nearly equal constant/end moments. Aggregated k/c forcing still uses one held
force imbalance and one homogeneous state. Source-clock queries contract the
original endpoint pairs directly, so individually overflowing endpoint
differences need not erase an admitted coalesced cancellation. Affine storage
and preparation/query work scores are admitted separately; constant banks
retain their smaller coefficient path. Warm queries allocate nothing and
evaluate no transcendentals. Reverse trajectories sum adjacent knot
contributions, transpose carrier-basis rotations and include the accumulated
source-phase clock chain. Refusal and cancellation return no partial result.

Six G1/G3 regressions in `tests/affine_tone_sensitivities.rs` cover independent
variational ODEs, the complete endpoint transpose, constant-envelope reduction,
duplicate cancellation, tiny near-DC damping, long-decay clock moments,
unequal-step multimode reverse trajectories, and budget/law/cancellation reuse.

`LinearForceTransition::pullback` returns `LinearForcePullback` cotangents for
initial displacement/velocity and both endpoint forces, with frequency, damping
and duration fixed. It transposes the same cached flow and force moments;
near-identity corrections reuse the independent decrements, while long decays
reuse the homogeneous retentions. It allocates nothing, evaluates no
transcendentals and does not differentiate preparation iterations. Inputs and
scaled results must remain finite, with no mutation or partial gradient on
refusal. A reverse trajectory carries both state cotangents and adds both
adjacent contributions at every shared force sample. These are derivatives of
the caller's scalar objective, not physical state values. This fixed-coefficient
pullback does not supply parameter/time or nonsmooth peak-selection derivatives.
G3 checks the transpose identity across free/under/critical/overdamped regimes
and tiny/finite steps, then compares all initial-state and shared-force
trajectory derivatives with independent finite differences. G0 covers invalid
seeds, overflow refusal, zero cotangents and deterministic retry.

`LinearForceTransition::prepare_sensitivities` additionally prepares analytic
physical-coefficient sensitivities for transient damping calibration and inverse
design. `PreparedLinearForceSensitivity` reuses the original forward step and
returns derivatives with respect to independent mass-normalized stiffness
`k = omega^2` [1/s^2] and viscous rate `c` [1/s]. The small-interval entire flow
moments and their exact doubling identities are differentiated during scalar
preparation; no finite-difference trajectories, pole separation, fitted decay
law or special critical-damping band is introduced. Independent homogeneous
retentions preserve late free and held-force damping-sensitivity tails; the
force-imbalance form keeps exactly balanced static states independent of damping.
Runtime `advance` and `pullback` allocate nothing and evaluate no transcendental
functions. The latter combines both coefficient cotangents with the original
state and force-sample transpose.

These are partials at fixed duration, source coordinates, initial state and
forcing. A reverse trajectory accumulates coefficient contributions from every
step and both contributions at shared force samples. At fixed `c`, frequency
sensitivity is `2*omega*dJ/dk`; at fixed ratio with `c = 2*zeta*omega`, it is
`2*omega*dJ/dk + 2*zeta*dJ/dc`, and the ratio sensitivity is
`2*omega*dJ/dc`. Zero stiffness and zero damping have their one-sided physical
limits without a numerical floor. Material, mass-normalization, basis, time,
heat/work and nonsmooth-observation derivatives remain outside this helper.
Unrepresentable sensitivities refuse without mutating prepared or physical state.
G1 regressions compare an independent high-precision augmented variational
matrix exponential through free, critical, near-critical, high-Q and stiff-drag
regimes, late critical tails and exact static balance. G3 checks subdivision of
one continuous ramp, complete transient reverse sweeps and the ratio chain rule;
G0 checks refusal and deterministic replay in
`tests/linear_forcing_sensitivity.rs`.

`ModalAcousticTimeModel::prepare_linear_forcing_sensitivities` binds these
derivatives to the bank's original nominal-clock ramp cache. Preparation is
shared by immutable mechanics clones. Physical rate/ratio changes replace the
kernel and require fresh preparation, even when two different free-coordinate
rates have the same infinite diagnostic ratio. Unchanged laws and observation
edits retain the valid cache. `advance_exact_linear_step_with_sensitivities`
previews the original motion and coefficient partials; `pullback_exact_linear_step`
supplies the complete state/force/coefficient transpose without allocation,
transcendentals or accepted-state mutation. Units and fixed-input chain rules
are the scalar owner's; no dimensional mass or changing basis is inferred.
`tests/linear_forcing_bank_sensitivity.rs` checks mixed free/soft/elastic modes,
legacy arithmetic, cache lifetime and clone isolation, then recovers a shared
dimensional dashpot coefficient from actual bank transients with unequal masses
through the critical-damping point, including the physical `dc/dC` chain rule.

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

`unit_force_response_into` returns that unit-force response for every mode at
any duration: the displacement and velocity each mode reaches from rest under
a unit held generalized force, through the transition `step_duration_into`
takes for that duration (cached for the full and half sample, prepared
otherwise, explicit damping rates included). One private function forms the
response for it, for the held-force solve's mobility and for
`probe_linear_map`'s force column, so the gains it returns are the numbers the
solve inverts. A step is affine in its held force, so these are its gains in
that force from any state, to the rounding of the response itself; fs-nlmodal's
`ModalStructure` builds its port map's tangent from them. No budget is
checked: a unit force is a direction, not a load. A nonpositive or nonfinite
duration and a slice not of one state per mode are refused, and no state
changes. `tests/unit_force_response.rs` pins the response to the solve's
mobility to the bit at the full and half sample and an off-clock duration, for
a bank of ratios and one of explicit rates, and to the step's own change under
a held force.
Each frame retains independent mode-local integrals of physical viscous power
`2*zeta*omega*v(t)^2`, with scale-aware contraction roundoff allowances. Their
sum supplies the aggregate heat. Undamped and free modes report exactly zero;
high-Q losses are not recovered by subtracting large mechanical energies.
Per-mode and aggregate force-work/energy balance checks remain separate and
refuse disagreement beyond their floating-point allowances before state commit.
Bank energies, heat, allowances and ordinary signed work use fixed-order
compensated reductions. Pressure accumulates velocity and displacement
quadratures separately, retaining quiet coherent contributions between
cancelling modes or quadratures. Restored states and static equilibria use the
same energy reduction for budget admission. State transitions, double-double
impulse work, physical ceilings and balance-tolerance formulas are unchanged.
`tests/modal_bank_reduction.rs` covers a 4096-mode high-Q ringdown, retained
signed work and pressure, and energy/pressure refusals that naive summation
could conceal; this is floating-point compensation, not a certified enclosure.
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
Every ZOH frame (`step[_into]` and `step_duration[_into]`) uses this same
equilibrium-relative pressure before applying the acoustic ceiling. A settled
load therefore consumes no pressure headroom, while load-on, release and
oscillation about that load retain their signed amplitudes. Zero transfers and
free coordinates need no compliance division; near-rigid modes can evaluate
the finite quadrature `omega*q - F/omega` without representing `F/omega^2`.
This baseline is specific to held-force observations; ramp, harmonic and
impulse observation laws are unchanged. Inline G0/G1 regressions cover static
silence across damping regimes and step APIs, analytic oscillation and signed
load transients, pressure refusals and free/soft-coordinate limits.

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
residual; an unresolved or over-wide mode refuses the whole call. A window whose
iteration budget ends short of its certified count refuses as
`WindowUnresolved`; a converged mode whose residual / |lambda| exceeds the cap
(or is not finite) refuses as `HalfwidthExceeded`, carrying its window, Ritz
value, relative half-width and the cap, so the two causes never share a code
(`tests/parallel_windows.rs` gives both on one pencil). The memory
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

### Atomic installation into the live modal runtime

`time::ModalAcousticTimeModel::transfer_basis` connects that physical transfer
to the accepted oscillator state. `ModalBasisChange` binds both original
mass-normalized bases to their runtime spectra (only square-root roundoff is
allowed). The supplied replacement is an admitted at-rest bank on the same
sample clock. The original rigid/mass-transfer owner projects the current live
state, and the replacement runtime rechecks displacement, velocity and energy
budgets before its state, damping, observers and prepared transitions become
visible together. Mode count, ordering and degenerate eigenvector orientation
may change. Refusal or cancellation leaves the complete old bank unchanged.

`ModalRuntimeTransfer` returns the original projection and signed parameter
work plus independently contracted old/new runtime energies and the residual
`H_new-H_old-W_parameter+E_discarded`. The original energy tolerance gates this
additional check. An instantaneous stiffness change can supply real parametric
excitation; its energy is parameter work, never hidden force work, heat or an
amplitude normalization. The caller owns event timing and the physical source
of that work. No continuous coefficient interpolation is inferred.

This cold update allocates and qualifies original sparse operators; the target's
existing exact transition cache is reused after commit. External force maps,
receiver filters and material memory remain their owners' responsibility in an
enclosing coupled transaction. The shared `fs-couple::modal_acoustic_time`
re-export exposes the same API. Five focused G1/G3 regressions in
`tests/modal_transfer.rs` cover changed-inertia momentum and subsequent damped
motion, degenerate basis rotation, truncation/free release, analytically exact
stiffness-switching parametric growth with a complete work account, and atomic
refusal. Full crate execution still requires the canonical workspace siblings.


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

### Mixed coherent and statistical mechanical interference

`harmonic::PreparedHarmonicInterferenceEnsemble` retains an explicitly authored
deterministic mean plus real and jointly proper complex force factors through
the same original coupled M/C/K pencil. Each factor spans its own complete
coordinates and frequencies. Real factors have zero mean and identity real
second moments. Proper factors have zero mean, identity conjugate second
moments and zero pseudo-covariance; real/proper cross moments are zero.
Physical coefficients include square-root statistical weights. There is no
inferred independence law, Gaussian assumption or normalization of source PSDs.

One existing harmonic response owner is prepared for each distinct frequency
globally and solves all mean/factor right-hand sides at that frequency.
Exact duplicate forces combine in DD within their own factor before solving;
opposite forces in separate factors retain separate physical variance.
Empty factors preserve their original indices, and a pure stochastic input
needs no fabricated deterministic tone. Completely empty sources refuse.

The coherent and statistical observers share DD contractions of
original-coordinate M*v, K*q, D*v and force, with D=(C+C^T)/2, plus the
physical flow/effort work when loaded. Real factors
retain sum and difference frequencies; proper factors retain the complete
within-factor difference-frequency interference. Expected kinetic/elastic
energy, input work, viscous heat, independent storage rate and power balance
use the existing normalized uniform/triangular observation windows. The
gyroscopic part changes response without being assigned dissipative heat.

Count, payload and product-work admission precede factorization. At most 64
random factors and 4096 authored tones are admitted, with the original 128-tone
ceiling per bank. Ordering records and bank headers are included in the
payload alongside 144*n retained scalar bytes per bank-frequency and one dense
response preparation. Warm observation allocates nothing and performs
6*n*sum_b(L_b^2) scalar products, with L_b the distinct frequencies in bank b.
Different factors need no cross contractions, and no warm linear solve or
matrix application occurs. These are work/storage bounds, not timing claims.

This reports ensemble second moments, not variances of energy or other
fourth-order statistics, random realizations, startup motion, inferred source
coherence, or a general near-null accuracy certificate. Native regressions use
independent real trajectories and finite mixed ensembles, exact duplicate
cancellation versus independent factors, proper beats and DC, real/unitary
factor-basis changes, coherent parity, and explicit resource/cancellation gates.

### Prescribed harmonic supports on the original pencil

`harmonic::PreparedSupportedHarmonicResponse` retains full M/C/K and a strictly
increasing original-coordinate support map. Under `exp(-i omega t)`, it solves
`A_ff q_f = F_f - A_fc q_c`, where `A = K-omega^2 M-i omega C`. Cross inertia,
viscous coupling and signed gyroscopic entries survive support elimination.
The original source passes the existing binary-scaled SPD/PSD admission; only
the free block is factored. Rigid full-body stiffness modes therefore remain
compatible with supported DC response. Empty and all-prescribed maps work;
singular free blocks refuse without invented mass, damping or penalty springs.

Full displacement/velocity and `R_c=(Aq-F)_c` are returned in the source chart.
Reactions are forces BY supports ON the structure; externally applied forces
at supported coordinates are subtracted once. Applied and support cycle work
and full symmetric-C viscous power are independently contracted in DD on the
shared scaled pencil. Gyroscopic coupling creates no heat, support power may
be signed, and cancellation in the work balance is not clipped. The original
free-row residual is checked after the solve. Repeated amplitudes and phases
reuse the same factor without changing accepted source state.

The complete source is bounded by the existing dense work cap and by
`320*n*n+512*n` numerical bytes, including full admission, retained full/free
operators, free LU and one solve. Allocation metadata/caller inputs are
excluded. Cancellation is polled during assembly, support elimination, solves
and observation. This is a finite-dimensional particular harmonic solution,
not a startup transient, actuator model, spatial/material certificate or a
forward-error guarantee near a pole. `tests/harmonic_supports.rs` independently
checks consistent cross-mass base excitation, relative support phases,
gyroscopic work, applied support-force redistribution, stationary parity,
DC, all-prescribed motion, singular free resonance and admission/cancellation.

### Physical receiver interference and averaged-reading covariance

Both prepared harmonic interference owners expose `prepare_receivers` for
signed real rows in the original pencil coordinates. A row observes
displacement, velocity, or acceleration; physical point/area/differential
pickups and linear strain/stress recovery use the caller's original map and
units. Receiver interpolation and frequency-dependent sensing dynamics are
not inferred. Projection uses DD before any temporal products, preserving
spatial cancellation. Acceleration differentiates projected velocity once,
without forming an unnecessary omega-squared intermediate.

`PreparedHarmonicReceivers::observe_pair` gives signed receiver cross moments
or self mean squares over the actual uniform/triangular observation window.
It reports deterministic means, coherent products, stochastic products, and
their total. Real factors retain sum/difference terms; proper factors retain
within-factor beats. Stochastic products accumulate directly, so a large
coherent bias does not erase small physical variance through subtraction.
Distinct declared factors are never combined coherently.

The separate `averaged_covariance` is the ensemble covariance of the two
window-averaged readings. Each factor's window amplitude is summed before
multiplication. It differs from the average instantaneous covariance: a
finite exposure can suppress a carrier in its averaged reading while the
mechanical mean-square motion remains positive. No stationary spectrum,
Gaussian law, fourth-order statistic, random realization, or near-null
relative accuracy certificate is inferred.

Count/payload/projection/product budgets are checked before projection
allocation. At most 128 receiver rows are retained, with DD real/imaginary
amplitudes per bank-frequency-row. Warm pair observation allocates nothing and
does no solve or original-coordinate contraction; requesting selected pairs
does not construct a dense receiver covariance matrix. Source response
residuals remain available and signed rounding is never clipped.

### Gain and phase control of physical receiver interference

`harmonic::PreparedHarmonicReceiverControl::try_with_shared_factors` binds
already solved receiver projections of a bias and complex source-control
bases. Each gain scales its source's complete authored tone bank. This
supports independent actuator/frequency bases and shared stochastic drivers
without solving or retuning the pencil during evaluation. Receiver order,
physical units, pencil, and source clock remain the caller's common model.

The constructor explicitly declares a shared factor basis: random bank k is
the same real or proper-complex variable in every input. Counts and laws must
match. Independent sources occupy different factor slots with explicit zero
columns elsewhere; unknown dependence is never treated as independence.
The deterministic mean is bank zero. These are authored second moments,
without a Gaussian law or fourth-order uncertainty claim.

`evaluate_into` reports the weighted finite-window expected mean-square
motion and its real Cartesian derivatives with respect to each complex gain.
For g=A*exp(i*phi), the phase derivative is
`-g.im*gradient.re + g.re*gradient.im`. The objective includes coherent and
stochastic contributions, all real sum-frequency terms, and proper-factor
beats under the original uniform/triangular window. It is the average
instantaneous vibration, distinct from the mean square of an averaged reading.

Exact frequencies merge across source bases before products. DD amplitudes
are recomposed first, then the derivative contracts twice the residual field
with each source response. This preserves quiet residuals under destructive
interference that can disappear when independently squared source powers are
subtracted. Half of the symmetric tone-pair kernels suffice; flat per-bank
arrays retain responses and reusable mixed/gradient scratch.

At most 32 source gains are admitted. Tone count, payload and conservative
product work are checked before allocation. Warm evaluation allocates nothing
and performs no mechanical solve or original-coordinate contraction. On
refusal the caller's gradient is unchanged; subsequent calls rebuild private
scratch. Source backward residuals remain available. This does not impose
actuator constraints, establish a causal broadband control filter, add startup
motion, or certify relative error near a null.

### Optimal source gains with an authored actuator-effort law

`PreparedHarmonicReceiverControl::prepare_optimum` caches the real gain-space
quadratic for one physical window and an explicitly supplied full SPD effort
metric W. The real coordinates are [Re(g0), Im(g0), Re(g1), Im(g1), ...].
DD contractions of the original response bases retain cross-source,
cross-frequency and shared real/proper factor terms. The bias/source and
source/source terms are assembled directly, not obtained by subtracting large
objectives or finite-differencing gradients. Receiver weights and source clocks
remain those of the original owner.

A query minimizes expected window-averaged vibration plus lambda*x^T*W*x for
an authored strictly positive lambda. The effort metric and coefficient carry
the caller's physical calibration; they are not inferred mechanical input
power, material loss, or hidden numerical regularization. W can retain
cross-actuator/quadrature coupling. An optional strictly positive total effort
cap uses a nonnegative multiplier mu and the gain system
(Q+(lambda+mu)*W)x=-b. Existing cancellable Cholesky factors operate on
diagonally equilibrated real systems. No additional damping, altered source
pencil, pseudoinverse threshold or automatic ridge is inserted.

Without a binding effort cap, one gain-system factorization suffices. A
binding cap uses the metric dual-norm estimate to bracket its monotone
multiplier, then retains a feasible endpoint throughout bounded bisection.
A query reports no approximate optimum when its factorization, effort
resolution or original-gradient stationarity requirements are unmet. Different
effort weights and caps reuse the physical quadratic.

Published vibration is recomputed by the original receiver-control owner
after DD field recomposition. Original-metric effort and componentwise
stationarity are checked independently of evaluating Q*x+b. The result
reports complex gains, coherent/stochastic vibration, effort, the penalized
objective, effort multiplier and numerical residual. These are floating
discrete-model optimality diagnostics, not interval or near-null forward-error
certificates. The specialized solve covers this fixed harmonic model and one
total quadratic effort cap; general nonlinear/per-actuator control constraints
remain authored optimization problems outside this owner.

Preparation admits source-quadratic product work and conservative additional
payload, including the original blocked factor's packing scratch. Queries
bound gain factorizations and poll cancellation in the original numerical
kernels. The borrowed controller freezes the source model for the query family;
the physical model is not modified by preparation or a refused query.

## Explicit physical-rate simultaneous tones

`PreparedModalTones::try_new_with_damping_rates` and its affine-envelope
variant admit complete physical oscillator banks for
`q'' + c*q' + omega^2*q = sum_j F_j(t)`. The supplied nonnegative rate is
the coefficient c [1/s], not c/2 or a damping ratio. Zero-stiffness coordinates
may have positive drag. No sampling clock, acoustic observer, Nyquist bound
or acoustic mode cap is inferred. All modes and tones remain under one
aggregate storage/preparation/step budget; widened integer counts prevent
capacity arithmetic from wrapping before coefficient allocation.

Pure `advance_mode_with_energy` and `advance_mode_affine_with_energy`
return motion, independent signed work and positive viscous heat through the
existing entire-convolution and paired-QR power owner. Exact duplicate tones
coalesce coherently; distinct close frequencies retain beats and all power
cross terms. Affine endpoint quadratures share the interval-start carrier
basis. Neither accepted state nor external source phase changes in a preview.

Existing model-based preparations retain their original ratio/rate arithmetic,
Nyquist/observer admission and immutable-kernel identity. An explicit bank is
unbound and cannot be passed to a model-mutating step. The controlled physical
constructor polls before admission, between cost rows and prepared modes/tones,
and before return; one bounded power preparation remains atomic. Caller-owned
physical state/energy acceptance and observation resolution remain separate;
these endpoint previews are not peak finders or interval certificates.

`PreparedModalTones::preparation_storage_bound_bytes()` exposes the SAME
payload-plus-peak-scratch formula used by native admission, with its original
mode count, distinct/input tone counts and envelope kind. It is a conservative
bound, not an allocator measurement, and retains the same exclusions for
caller workspace and model-owned state. Physical consumers can reserve an
old bank while preparing its replacement without duplicating the native
coefficient-size formula or understating simultaneous storage.

`tests/physical_tone_banks.rs` covers duplicate cancellation, near-tone beats,
affine work/heat against independent ODE integration across damping regimes,
analytic rigid drag, high-frequency admission, preserved model binding,
aggregate budgets and cooperative preparation cancellation.


## Physical port loading on prescribed harmonic supports

`PreparedSupportedHarmonicResponse::try_with_port_impedance` assembles the complete
original `K - omega^2 M - i omega C - i omega B Z B^T` before prescribed-support
elimination. The original real force map also supplies `B^T v`, including every
moving support row; mutual, reactive and nonreciprocal impedance entries remain.
Only the free principal block is factored, allowing physical load to resolve an
otherwise singular bare free resonance. Empty and all-prescribed maps remain
valid; finite velocity impedance vanishes at DC.

The original port owner admits Hermitian-PSD load data and separately bounds its
payload and projection work. It observes full flow, opposing effort and absorbed
boundary power independently. Reactions retain the loaded original equations;
applied work plus signed support work balances intrinsic viscous heat plus port
absorption. Boundary work is never assigned from the balance or reclassified as
material heat. Structural payload, binary scaling, source residual, factor reuse
and cancellation retain the existing supported-solve contract. The load's owner
must supply frequency validity and any radiation, heat or storage interpretation.
Regression cases cover nonreciprocal mutual loading, moving-support port flow,
free-resonance regularization by a physical load, static/empty/all-support limits,
independent work accounting, refusals and cancellation.


### Coherent tone work and heat control gradients

`PreparedModalTones::prepare_power_gradients` admits a borrowed
`PreparedModalTonePowerGradient` under separate fixed-array scratch and
single-mode query-work caps. No new numerical coefficients are prepared.
`gradients_mode` and `gradients_mode_affine` return separate input-work and
viscous-heat derivatives for initial q/v, every original source amplitude or
both envelope endpoints, and independent common source phase. The original
paired velocity/force QR factors are contracted before products: work
differentiates `(R*x)^T*(T*x)` and heat differentiates `||R*x||^2` with their
original physical scales. All coherent cross terms and the aggregate force
imbalance chart remain. No independent-tone heat sum or endpoint-energy
subtraction is used.

The physical input transpose is applied per factor row before accumulating
endpoint covectors. Exact duplicates map back to every original source, and
affine endpoints are coalesced before differences, admitting representable
cancellation even when individual deltas overflow. Both endpoint pairs use
the interval-start carrier basis. Reverse trajectories add power and motion
covectors, sum shared knot contributions and transpose carrier rotations.
Mechanics k/c, source frequencies, duration and basis/mass normalization stay
fixed on this API. Warm queries are O(n^2), allocate nothing, evaluate no
transcendentals and poll cancellation between factor rows and source mappings.
Refusal publishes no partial result. G1/G3 tests in `tone_power_gradients.rs`
cover independent physical variational ODEs, duplicate/near-tone controls,
complete energy transpose, constant-envelope reduction, the active 32nd source,
unequal-step multimode reverse chains, units, budgets and reusable refusals.
Finite-arithmetic gradients remain estimates, not rigorous enclosures.

### Coherent source-frequency work and heat gradients

`PreparedModalTones::prepare_power_frequency_gradients` adds separately admitted
source-clock work/heat derivatives for constant and affine-envelope banks. It
borrows the original bank. Constant-envelope bases prepare existing affine
TonePower factors per mode, under conservative affine byte/work admission and
separate affine-layout query caps. Their exact direction is
`dF_j/dOmega_j = t*(S_j*cos(Omega_j*t)-C_j*sin(Omega_j*t))`; contracting the full
affine endpoint power covector with `h*(S_j,-C_j)` includes direct source work and
changed coherent velocity together. Original input order is retained at duplicate
frequencies, including opposite nonzero independent detuning derivatives when
nominal forces cancel. Frequencies sharing an authored fundamental require their
explicit parameter-chain weights.

Local derivatives hold initial state, interval-start quadratures, mechanics,
duration and basis/mass fixed. Multi-interval reverse sweeps also include motion
and source-phase advancement clock chains. Affine bases use the complementary
quadratic power companions described below. Warm
queries allocate nothing, evaluate no transcendentals and poll between factor
rows/source contractions; cancellation or nonfinite results publish no partial
derivative. These are operation properties, not measured timing claims.

Five G1/G3 tests in `tone_power_frequency.rs` cover independent power variational
ODEs, analytic duplicate/DC detuning, relative tiny near-DC heat derivatives,
scalar/32nd-source agreement, unequal-step multimode clock reversal, unit rescaling
and resource/law/refusal reuse. The six existing coherent power-control tests also
pass on this source. Gradients remain finite-arithmetic estimates, not enclosures.

## Complete physical power parameter jets

`harmonic_mean_power_parameter_jet_with_cancel` contracts original
power-conjugate peak phasors into P, P_omega, P_omegaomega, P_theta and
P_omega_theta. It includes all source/flow/effort product rules, including
both distinct mixed cross terms. Binary-normalized double-double products
are accumulated before recovering each physical result. Work is O(n), with
constant scratch, no allocation and a cancellation checkpoint per coordinate.
Shapes must be complete, nonempty and finite; refusal returns no partial jet.
Signed boundary work does not identify constitutive heat or stored energy.
G1/G3 tests in `harmonic_physical_power_parameter_jet.rs` check independent
resistive heat with a changing reactive/source law and cancelling products
above the ordinary floating-point range.

### Forced response on the general quadratic source pencil

`quadratic_modes::forced` consumes the original admitted
`PreparedQuadraticModes` and factors its same compensated polynomial at
s=-i*omega. It retains SPD M, signed/nonsymmetric C and K, and the original
binary coordinate/time scaling. Physical forces use the dual scale; repeated
load cases reuse one native cancellable complex LU. DC gives local static
compliance, while true singular resonance refuses without a damping floor.

Original separate force actions independently gate each published
displacement. Shared `evaluate_harmonic_balance` reports signed applied,
ideal-reaction, inertial, velocity-operator and stiffness powers. Active and
circulatory work is not classified as heat; a particular response does not
imply decay of unstable transients. Residuals are backward checks, not
forward-response bounds near resonance. Explicit factor/work/payload gates
and cancellation are transactional. The existing passive harmonic API retains
its admission policy.

Shared compensated polynomial arithmetic prevents a separately rounded
squared frequency from erasing original binary detuning. Six full-owner
regressions and two arithmetic tests cover static and damped response,
active/signed operators, nonsymmetric and gyroscopic work, represented narrow
detuning, extreme force scaling and atomic refusal. These new tests are
source-reviewed but unexecuted while the development executor is disconnected.

### Coherent stiffness and viscous-rate work/heat gradients

`PreparedModalTones::prepare_power_mechanical_gradients` separately admits
coherent finite-interval work and viscous-heat derivatives with respect to
independent mass-normalized stiffness `k=omega^2` and physical viscous rate c.
Both constant and affine envelopes retain every source interference term,
including exactly coalesced duplicates. The original resolved scalar flows and
analytic Frechet/envelope directions feed joint velocity/force QR; both
mechanical tangents follow the same orthogonal rotations without differentiating
QR pivots. Keeping force-only rows preserves nonzero work derivatives at a
motionless forced equilibrium. The fixed conditioning chart contributes its
original initial-force-imbalance derivative `-q*h/scale` before each physical
contraction. `dQ/dc` includes `integral(v^2 dt)`, including at c=0.

Cold preparation admits 8*n^2 doubles per mode plus peak scratch under explicit
work/storage caps. Warm queries use allocation/transcendental-free O(n^2)
contractions with cancellation between rows and no partial publication. Initial
physical state, interval-start source amplitudes/endpoints, source clocks,
duration and basis/mass normalization stay fixed. Trajectory reversals also
propagate state cotangents and sum motion and local-power contributions; authored
natural-frequency/damping-ratio/material chains remain explicit. Six G1/G3 tests
in `tone_power_mechanical.rs` cover independent variational ODEs in all damping
regimes, differentiated independent energy balance, analytic static DC and
zero-loss limits, tiny near-DC heat, duplicate cancellation, scalar/32nd-source
agreement, unequal affine-step reversal, physical units and budget/cancellation
reuse. Derivatives remain finite-arithmetic estimates, not enclosures.

### Affine-envelope source-clock work and heat

`prepare_power_frequency_gradients` now admits original affine tone banks,
with `gradients_mode_affine` retaining both original interval-start carrier
endpoint pairs. The exact direction is
`h*[u*(1-u)*(S0,-C0)+u^2*(S1,-C1)]`, with `u=t/h`. The original paired
power-factor owner prepares degree-two companions in this complementary clock
basis. Existing scalar endpoint moments supply the flow; structured source
blocks and the same QR/doubling owner compose long intervals. No separate ODE
solver, numerical source perturbation, independent-tone heat sum or endpoint
energy subtraction is introduced.

The nominal affine envelope uses equal complementary coefficients. Its clock
transpose contracts the original endpoints per factor row, avoiding individual
endpoint differences and a subtraction of late-time complementary covectors.
Exactly coincident sources retain separate detuning controls even when their
nominal envelopes cancel. The original force-imbalance chart, initial transients,
source phase and cross-frequency products remain. Multi-interval derivatives
still require the motion/state and advancing source-phase chains.

The quadratic layout has n=2+6K, up to 194 coordinates for 32 distinct tones.
Separate cold storage/work and fixed-array warm scratch/work admission precede
allocation. Warm work is bounded by `32*n^2+128*n*input_tones`; queries allocate
nothing, evaluate no transcendentals and preserve reuse after cancellation or
refusal. Four new tests check independent physical variational ODEs and constant
reduction, cancelling duplicates with overflowing individual endpoint deltas,
weak complementary work after a long decay, all 32 distinct frequencies,
physical-unit conversion and resource/cancellation boundaries. Together with
mechanical and directly affected original tests, 26 native modal tests pass.
These are finite-arithmetic local derivatives, not rigorous error enclosures or
measured performance guarantees.

### Symmetric quadratic physical-work curvature

`harmonic_quadratic_power_curvature_with_cancel` evaluates original
`Re(v^H*sym(C)*v)/2` and both frequency derivatives for a fixed finite
row-major C. It retains all product-rule terms and permits signed directions;
observation does not admit a constitutive law. Scalar matrices also observe
independently owned material resistance or spring-storage coefficients.

Upper-triangle traversal combines each off-diagonal coefficient pair before
field contraction, preserving the sum's low word. Exact gyroscopic pairs
cancel before reaching physical-work accumulators. Subnormal coefficients are
never halved before multiplication; overflowing same-sign coefficient sums
remain separate normalized products. Neither C*v nor twice a large coefficient
is formed. The existing binary-normalized double-double work owner provides
O(n^2) work, constant scratch, no allocation, and row cancellation checkpoints.
Malformed, nonfinite, unrepresentable or cancelled observations return no
partial result. Focused tests cover complex physical curvature, tiny symmetric
work beside huge gyroscopic entries, subnormal/large coefficient products and
coefficient remainders surviving leading quadratic cancellation.

### Complete physical quadratic-work parameter jets

`harmonic_quadratic_power_parameter_jet_with_cancel` observes
`Re(v^H*sym(C)*v)/2` through value, first/second frequency, parameter and mixed
frequency/parameter order. Its complete original field jet and optional C_theta
retain every product-rule term, including direct matrix changes and both mixed
field/matrix contributions. C and C_theta are frequency independent; signed
directions are observations, not newly admitted positive constitutive laws.

The same upper-triangle coefficient and normalized contraction primitives serve
the curvature-only API. Exact skew pairs cancel before physical accumulation;
coefficient low words, subnormal factors and overflowing same-sign sums retain
the earlier range guarantees. A skew base matrix can have a nonzero symmetric
direction, which contributes physical work. No C*v array, heap allocation,
factorization or unused parameter computation is added to curvature-only calls.
Work remains O(n^2), with constant scratch, input admission and row cancellation.

Focused native regressions compare a complex moving-field/matrix polynomial,
extreme coefficient and gyroscopic scales, direct and mixed loss directions,
malformed inputs and cancellation/replay; the existing curvature tests remain.

### First-order quadratic physical work for resonance refinement

`harmonic_quadratic_power_derivative_with_cancel` observes original
`Re(v^H sym(C) v)/2` and its complete frequency slope without requiring a
curvature field or evaluating second-order work. It shares the higher-jet
owner's original symmetric coefficient pairs and normalized double-double
products, preserving gyroscopic cancellation, weak loss and subnormal
coefficients. Signed forms retain their physical sign. First-order results
remain available when an unused second derivative would overflow. The
allocation-free traversal uses constant scratch and O(n^2) work, with bounded
row cancellation checks and no partial result after refusal.
