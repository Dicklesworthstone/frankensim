# fs-atmo — CONTRACT

Bead: frankensim-wf-root-guzez.4.5 (E3.3a, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md §5.4
(ROUND 6 steady state). Evidence: air-state-v1.json (E1.8) under the frozen
registry (E1.7). Conventions: frame-conventions-v1; heights are altitude h
above the aerodynamic ground plane.

## Purpose and layer

L2. The atmosphere foundation: FlatSiteLogLaw mean wind, wall-compatible
solenoidal modal turbulence with exact analytic derivatives, and the
provenance-bound `sample_air_state` API. E3.3b adds the exact-discrete OU
amplitude evolution (sequential checkpointed state) and a fit of selected
one-point stresses from the Mann-class spectral tensor; E3.3c the optional
FetchAdjustedMassConsistent mode.

## Public types and semantics

- `FlatSiteLogLaw { scenario_effective_z0, displacement_height,
  reference_height, reference_speed }` — ONE scenario-level z₀; `speed(h)`
  = (u*/κ)·ln((h−d)/z₀) above the sublayer, 0 inside; analytic `dspeed_dh`.
- `TurbulenceField::build(seed, n_modes, sigma, length, u_adv)` — philox
  counter-addressed modes (StreamKey{seed, ATMO_KERNEL, tile = mode});
  `sample(x, y, h, tick)` → velocity + EXACT analytic gradient. Per mode:
  ψ_h ∝ sin(k_h·h)·cos θ, ψ_v ∝ cos(k_h·h)·cos θ, u = ∇×ψ, with the
  frozen-phase clock θ = k_x(x − U_adv·t) + k_y·y + φ.
- `AirScenario` / `DEC17_AIR` — ρ/μ/T/p with provenance (E1.8 derivations).
- `Atmosphere::sample_air_state(x, y, h, tick)` → `AirState` whose
  `dynamic_pressure_pa()` and `reynolds(chord)` derive from the SAME
  provenance-bound state as the velocity.
- `Atmosphere::sample_flow` and `EvolvingAtmosphere::sample_flow` expose
  the exact original mean-plus-turbulence velocity and gradient without
  imposing constant scenario thermodynamics. L3/L4 material compositions
  attach their own coherent gas/profile; no upward dependency is introduced.
  Original air-state sampling delegates to this same flow evaluation.
- `evolving::EvolvingTurbulenceField::at_anchor(basis, seed, times)` —
  opt-in space–time field with one declared OU correlation time per spatial
  mode. Two independent OU quadratures per mode multiply the vector potential
  before its curl; both have stationary variance 1/2. This preserves the
  basis realization's complete phase-averaged stress tensor and makes its
  conditional one-point statistics horizontally homogeneous. Per-mode
  covariance decays exponentially while its phase advects with the declared
  mean speed. `sample` reads the current signed tick; `sample_batch` fills
  caller-owned storage without allocations or random draws. `advance_to`
  alone evolves the sequential state. A self-contained checkpoint retains
  basis, times, and both OU paths; resume is bitwise.
- `Atmosphere::into_evolving` / `evolving::EvolvingAtmosphere::sample_air_state`
  — the same mean law and provenance-bound air properties with the evolved
  vector velocity and analytic gradient, so dynamic pressure and Reynolds
  use that same instantaneous field.
- `fit::{fit_amplitudes, build_fitted_field}` fits (u², v², w², uw) at a
  declared reference height, with a positive-semidefinite amplitude covariance
  and explicit finite-basis residuals. The fit retains its seed, Mann target
  and target shear stress. Building a realization checks that binding and
  reuses the retained stress instead of repeating the target quadrature.
  Γ = 0 fits the variances with zero cross-correlation, without dividing by a
  zero shear-stress target. Nonfinite heights/covariances/advection refuse.
- `mann::stress_integrals` supplies those four target stresses over the declared
  band `0.01 <= |k|L <= 1000`. A deterministic 64×48×24 Gauss–Legendre rule
  integrates dimensionless log radius, polar angle and symmetry-reduced azimuth.
  Eddy lifetime is evaluated once per radius; angular trigonometry is reused.
  Compensated accumulation retains shear stress, and reflection symmetry makes
  `uv` and `vw` exactly zero. The isotropic limit uses its analytic angular
  integral. Dimensional scaling is `alpha_eps23 * L^(2/3)`, avoiding artificial
  overflow from dimensional wavevector powers. No omitted-band energy is added.
- `mann_field::{MannSpectralModel, HomogeneousMannField}` provides a separate
  homogeneous prescribed velocity field from the full Mann spectral tensor.
  Each finite-band wavevector has two transverse polarizations whose complete
  covariance equals its weighted tensor, including cross-components. Four
  independent unit-variance OU paths supply both polarizations' sine/cosine
  quadratures. Geometry and amplitude seeds are separate; vector advection,
  decorrelation times, mode count and both kL cutoffs are explicit. Initialization
  performs the tensor algebra; velocity/gradient samples allocate no storage.
  `covariance` exposes the exact represented spatial or space-time tensor,
  including finite-band and quadrature effects. Checkpoints retain model inputs
  and all four paths and rebuild the deterministic model on resume. The
  original homogeneous observation has no wall condition; z = 0 is an ordinary
  plane. The explicit plane-blocked observation below adds its own boundary
  treatment without changing the retained model or quadratures.
- `mann_field::fixed_step::FixedStepMannField` starts from a retained homogeneous
  snapshot and evolves its four quadratures at an explicitly supplied physical
  time step. Each mode uses the exact OU transition
  `rho = exp(-dt/tau)`, `sigma = sqrt(-expm1(-2dt/tau))`, following Gillespie
  (1996), DOI `10.1103/PhysRevE.54.2084`. Its distinct counter-addressed stream
  is bound to the authored step size and seed; changing the time step does not
  preserve the same realization or the original 120 Hz path. The original
  integer-tick owner and its checkpoint semantics remain intact.
  Sampling reuses the same spectral tensor, transverse polarizations and
  spatial filters, with advection applied once at the retained phase origin
  plus elapsed physical time. The immutable model is shared; accepted/trial
  amplitudes and transition coefficients are prepared once. Transition work
  and actual numerical capacities are admitted explicitly, and cancelled or
  refused advances retain the accepted state. Checkpoint/resume and temporal
  state copying bind the complete model, seed, clock and correlation times.
  The nominal numerical payload is 88 bytes per mode, excluding immutable
  spectral geometry and fixed metadata; shared coefficient storage is counted
  conservatively per owner. Authored OU decorrelation is a synthetic temporal
  model, not resolved nonlinear eddy turnover, and is not temporally bandlimited.

### Balanced finite-mode Mann synthesis

`MannSpectralModel::build_with_rule` admits the explicit
`MannQuadratureRule::BalancedEnergyV1` geometry rule. The original `build`
constructor selects `LegacyLogRadiusV1` and preserves its deterministic
geometry and samples. The balanced rule uses an even count of 2 through 256
modes and the authored finite kL band. It improves how a limited mode budget
resolves the full anisotropic tensor, including spatial covariance.

Preparation evaluates the original Mann tensor on a fixed 32×32×16 proposal
grid in log radius, vertical direction cosine and quarter azimuth. Each
diagonal component contributes equal proposal probability; a 5% uniform
mixture retains strictly positive support. Those normalizations affect only
where modes are sampled. The final positive integration weight includes the
exact conditional-CDF Jacobian and the physical radial volume factor; every
mode retains the full original tensor at its wavevector. No component is
whitened, fitted or rescaled to a desired realized covariance, and no energy
outside the authored band is added.

Each quarter-domain draw supplies a pair reflected in ky, while the real
sine/cosine representation supplies conjugate wavevectors. The two partners
retain distinct mode indices and independent paths for all four OU
quadratures. Their represented equal-time, one-point ensemble `uv` and `vw`
covariances cancel exactly. That cancellation is not asserted at arbitrary
separations or at nonzero lag with unequal authored partner lifetimes.

The proposal uses 16,384 tensor evaluations and less than 1 MiB of temporary
numerical storage once during model construction. It is then discarded.
Retained mode storage and the O(modes) velocity/gradient sampling kernel are
unchanged. Both temporal owners, bridges and existing linear observation
operators consume the resulting model. Homogeneous checkpoints explicitly
retain the rule, fixed-step checkpoints retain the full model, and prepared
spatial filters reject a different rule even when other model inputs match.

`mann_field::importance::tests` independently checks the positive coordinate
Jacobian, all nine tensor components per mode, physical scaling, independent
paired amplitudes over 4096 realizations, and both clocks' replay and filter
binding. Held-out geometry seeds 32 through 63 reduce full-tensor normalized
RMS stress error at Γ = 0, 1, 3.9, 7 and 10 with both 64 and 256 modes.
Independent refined continuum integrals also check full-tensor RMS error at
two oblique separations for Γ = 0, 3.9 and 10. Each component error is divided
by `sqrt(R_ii(0) R_jj(0))` before averaging its square over all nine entries
and held-out seeds. These checks establish the
declared cases; they do not guarantee every component, seed, finite band or
separation improves, or establish agreement with measured atmospheric data.

### Full-tensor Mann observations above a plane

`mann_field::wall::MannPlaneWall` declares a finite plane `z = height_m` in
the authored Mann frame. Both temporal owners expose `sample_wall_blocked`
and `covariance_wall_blocked`; the spectral model also exposes equal-time
covariance. Sampling preserves the original accepted quadratures and adds the
harmonic, irrotational correction of Hunt & Graham (1978),
<https://doi.org/10.1017/S0022112078000130>. With horizontal wavenumber magnitude
`q` and height `h`, a modal amplitude `a` receives
`a_z * (i*kx/q, i*ky/q, -1) * exp(i*k.wall_point - i*k.U*t - q*h)`.

The correction enforces zero normal velocity at the plane, preserves
solenoidality and decays aloft. All nine derivatives are analytic. Covariance
transforms both original polarizations and all cross-components at both
observation points, retaining the original OU decay and phase advection.
Both coordinates are required because normal homogeneity is lost. The common
first-point phase is removed before covariance evaluation. Near-wall
differences use half-angle identities and `expm1`; exact vertical modes have
a separate zero-correction branch, and an exactly vanished correction does
not evaluate an irrelevant distant-wall phase.

Only finite points in the exact upper half-space are admitted. The same
retained modal arrays are read without new allocations, draws or clock
changes. This is an outer-flow slip-plane construction, without a viscous
no-slip layer, ground-generated turbulence or a calibrated near-ground
spectrum. Homogeneous sinc filters alone do not supply the evanescent
correction's spatial average. Six focused Rust regressions accompany the
point observation. After execution was restored,
all six passed in the native debug build, including independent tensor
covariance, all nine gradients and extreme far-field recovery.

`mann_field::wall::filter::MannWallBoxFilter` supplies exact uniform averages
over boxes, parallelepipeds and collapsed faces/lines/points. Both temporal
owners expose `sample_wall_blocked_filtered` and
`covariance_wall_blocked_filtered`; the spectral model exposes the equal-time
covariance. The original homogeneous gain is combined with the analytic
evanescent gain `product(sinh(z_j)/z_j)`, where
`z_j = -q*a_j.z + i*(kx*a_j.x + ky*a_j.y)`. Both observations of the full
tensor receive their own height-dependent response; energy is not rescaled.
The whole footprint must remain above the plane, without clipping.

Prepared correction gains are normalized at the footprint's lowest height,
preserving deep wall-touching averages even when the center correction has
underflowed. Compensated increments preserve tiny near-wall differences.
Numerical gain storage is zero for points, 8 bytes/mode for tangential
footprints and 32 bytes/mode with vertical extent; actual capacity is admitted.
The source-bound observation allocates nothing and draws no randomness when
sampled. Zero extent retains the original point arithmetic. Original physical
clocks, spectra and OU paths are unchanged. This averages velocity and its
analytic gradient, not nonlinear momentum/energy fluxes or unresolved physics.

Six native filter regression groups passed against independent spatial and
two-observer covariance quadrature, gradient differences, solenoidality,
wall parity, both clocks/replay, deep/tiny footprints and admission. The
complete affected Mann module suite passed 28 tests. These establish the
represented finite-field identities, not a DNS or empirical calibration.

## Invariants

- Wall-compatible solenoidality: u is a curl — ∇·u ≡ 0 analytically; the
  analytic-gradient trace cancels to machine precision (V-04a relative bound).
- Wall-compatible parity: sin(k_h·0) = 0 makes the vertical component vanish
  at h = 0 IDENTICALLY (bitwise zero, battery-asserted).
- Pointwise z₀ insertion is FORBIDDEN while solenoidal claims stand
  (Round-2): the law carries exactly one scenario-level effective z₀.
- The analytic gradient is term-wise exact (central differences converge to
  it at measured order 2).
- Mode draws are a pure function of (seed, kernel, tile): independent of
  mode count and of other fields built in between for the original wall basis.
- Homogeneous Mann modes instead use an explicit finite-band quadrature whose
  geometry depends on its mode count. Every velocity polarization is transverse
  to its wavevector, so its instantaneous analytic divergence vanishes without
  a wall projection. Unresolved energy is omitted, never rescaled away.

- `ou::OuMode::from_correlation_time(sigma_st, tau)` — exact-discrete
  ρ = e^(−Δt/τ), σ_innov = σ_st·√(1−ρ²) (the stationary law is a FIXED
  POINT of the transition, battery-proven).
- `ou::StationaryOuPath` — StationaryOuPathV1: exact stationary draw at
  the FIXED anchor tick −3840; philox counter-addressed innovations
  (pure function of (seed, OU_KERNEL, mode, tick)); sequential
  `advance_to` (never backwards, span-capped), `checkpoint`/`resume`
  with BITWISE reconstruction equivalence; every preroll window is a
  suffix of ONE physical path.

## Error model

Typed `Refusal { code, message, ranked_repairs }`. Codes:
`non-finite-input`, `z0-outside-domain`, `displacement-invalid`,
`reference-height-invalid`, `reference-speed-invalid`,
`mode-count-invalid`, `turbulence-params-invalid`, `below-surface-query`;
OU layer: `ou-params-invalid`, `ou-advance-backwards`,
`ou-advance-span-exceeded`, `ou-checkpoint-invalid`.
Caps tested at cap AND cap+1.
Evolving-field admission also reports `evolving-params-invalid`,
`evolving-checkpoint-invalid`, and `sample-count-mismatch`. The two
quadratures must share a clock and use distinct prescribed stream keys;
invalid advance requests and invalid batches leave state/output unchanged.

## Determinism class

Deterministic: det:: transcendentals, philox counter-addressed draws, no
global state; V-04a golden pinned under
`org.frankensim.fs-atmo.v04a-golden.v1` (golden-bump protocol).

## Cancellation behavior

Synchronous pure functions; nothing to cancel.

## Unsafe boundary

Workspace `deny(unsafe_code)`; no unsafe.

## Feature flags

None.

## Conformance tests

`tests/v04a_battery.rs` — V-04a EXECUTED: relative divergence < 1e-12 over
288 probe/tick combinations; bitwise wall parity; FD-vs-analytic gradient
Richardson order ≈ 2 with finest-step error < 1e-6; air-state consistency
(Re/q hand-recomputed; log law recovers the reference point; E1.8 rho with
provenance string); seed determinism bit-identical + counter-addressed
partition law; refusal caps; log-law identities (U(d+z₀) = 0 exactly,
monotone, instrument-height ordering); pinned golden. JSONL receipts.

`tests/ou_battery.rs` — E3.3b-i / V-17 core: exact-discrete relations vs
hand closed forms + the stationary fixed point; no-spin-up window
variances; the suffix property (pause/resume vs straight advance,
bitwise); checkpoint reconstruction equivalence (bitwise) + checkpoint
refusals; the 39-minute (284k-tick) recurrence battery in checkpointed
chunks; caps at cap AND cap+1; pinned OU golden.

`evolving::tests` — G0/G3/G5: all nine stationary stress components against
an independent curl calculation (2e-14 absolute tolerance), horizontal
homogeneity, advected OU autocovariance and stationary variance over 4096
independent paths (0.07 normalized tolerance), exact wall parity, analytic
divergence below 1e-14 on the fixture, gradient finite-difference order
2 ± 0.05, read-only query/batch equivalence, bitwise checkpoint/resume, and
velocity/gradient/air-property consistency.

`tests/mann_battery.rs` checks target integration against independently converged
high-resolution references at Γ = 0, 1, 3.9, 7 and the admitted cap 10, with
2e-6 relative tolerance on all four nonzero stresses. An independent converging
Simpson energy integral checks the isotropic limit. Exact tensor symmetries and
dimensional rescaling, including otherwise ill-conditioned wavevector scales,
are covered. These numerical checks concern the declared finite band; they do
not establish a continuum spectral tail or empirical atmospheric calibration.

`tests/fit_battery.rs` additionally checks the finite Γ = 0 limit, seed/target
binding, input refusal, and the nonzero-shear fit golden. The Mann-integral and
fit goldens change with the corrected quadrature target. The
existing per-truncation and bootstrap checks concern the fitted one-point
stresses, not full spectral/coherence agreement.

`mann_field::tests` checks all nine modal tensor components for Γ = 0, 3.9
and 10, independent isotropic radial-energy integration, sheared two-point
covariance at oblique separations over 4096 independent amplitude seeds,
advected OU covariance over 2048 independent seeds, analytic gradients and
divergence, absence of an invented wall condition, exact checkpoint/resume,
and model/state refusals. These validate the represented finite-band model;
they do not certify continuum quadrature convergence or measured wind spectra.

## No-claim boundaries

- The modal amplitude decay is an ANALYTIC-CONSTRUCTION placeholder: no
  statistical claim of any kind (spectra, coherence, TI) — those are
  E3.3b's Mann-class fit under V-04b1/V-04b2 with its own artifacts.
- `TurbulenceField` retains static per-realization amplitudes and pure mean
  advection. `EvolvingTurbulenceField` supplies the declared finite-mode OU
  temporal model; its correlation times are explicit model inputs. Its
  conditional stationarity and covariance tests do not validate atmospheric
  spectra, intermittency, nonlinear energy transfer, or the full Mann target.
  The existing Flyer lifecycle's reduced gust forcing is not changed by
  opting into the new fs-atmo API elsewhere.
- Neutral stability only (declared, per air-state-v1); other classes are
  future additive modes.
- No measured/historical gust trace or modeled thermals. Compact gust events
  use only the explicit kinematic potential construction below.
- `DEC17_AIR` is the ensemble MEAN of the pre-registered distributions;
  H-case members draw their own scenario per the frozen registry.

## Exact native-step OU velocity observations

`fixed_step::MannStepAverage` is an opt-in, prepared observation of one complete
physical interval of the existing advected `FixedStepMannField`. The next endpoint
and the interval-mean quadratures are sampled from their joint Gaussian law. The
endpoint retains the original counter stream and multiply/add order; a disjoint
`0x4d464156` stream supplies the independent integral bridge. Linear interpolation
between endpoint samples is not substituted for that bridge.

Preparation admits the original immutable spectrum, clock, phase origin and RNG
identity, plus an explicit numerical workspace ceiling (104 additional bytes per
mode, excluding shared source definitions and fixed metadata). Each joint advance
uses eight normal draws per mode, reuses two amplitude banks, polls cancellation,
and publishes endpoint and observation together only after all admission passes.
Small decay/phase coefficients use analytic divided-difference and covariance
series; scaled complex arithmetic admits large finite phases without squaring
them. Nonrepresentable decay, phase, coefficient or bridge scale is refused.

The observation supports the original point, oriented box, harmonic slip-plane,
and complete wall/box kernels, with their analytic spatial gradients. Those static
linear operators act on the same averaged quadratures; querying draws no noise and
allocates no amplitude storage. Existing endpoint checkpoint format and replay
remain unchanged, and determine all future means. The observation workspace is
not a new accepted flow owner and need not be checkpointed to replay future steps.

This is the exact linear time average of the represented finite-band OU velocity
model, conditional on its authored correlation times. It does not average
nonlinear momentum or energy fluxes, supply unresolved stress/kinetic energy,
validate physical eddy turnover, or make OU forcing temporally band-limited.

Verification: the native `mann_field::` tests include an independent continuous
stochastic-convolution quadrature oracle,
small-decay and rotating/extreme-phase limits, 4096-realization joint moments,
unchanged endpoint/restart bits, every publication cancellation boundary, and
all four spatial observations with divergence and gradient checks.

`fixed_step::TranslatingMannStepAverage` observes the same OU source along
`x(t) = x_end - V * (t_end - t)`, so encounter frequencies use `k·(U_adv - V)`.
Its samples take the probe/footprint position at the interval end and return
velocity components in the original frame; they do not subtract probe speed.
Point and fixed-orientation footprint averages support three-dimensional uniform
translation. Wall-blocked observations require `V_z = 0`: normal motion changes
the evanescent amplitude during the interval and is explicitly refused. The
distinct type cannot be passed as a stationary inlet observation. Every sample
from one owner shares its one translation velocity and integral bridge; separate
velocities do not establish a joint continuous-path covariance.

`advance_one_with_translating_average` reuses the original endpoint transaction,
work budget, storage and cancellation path. Ordinary endpoint bits and restart
behavior remain unchanged. Focused regressions exercise a translating finite
sensor against the exact Galilean-frame stationary observation, all spatial
observations, wall admission, replay, and a complete wavelength traversal versus
an observation travelling with the eddies in the nearly frozen limit.

The older 120 Hz OU parameter constructor uses `-expm1(-2 dt/tau)` when
`dt/tau < sqrt(epsilon)` and avoids overflowing `120 * tau`, preserving finite
nonzero innovation variance at long correlation times. Ordinary-path arithmetic
and its existing golden remain unchanged. At decay below f64 resolution, `rho`
can still round to one; this correction preserves the representable short-step
diffusion and does not claim to resolve mean reversion over that enormous time.

## Replayable conditional OU observations between fixed steps

`fixed_step::MannOuBridge` prepares the accepted state of a `FixedStepMannField`
and that owner's exact next endpoint without advancing it. The right endpoint
uses the original fixed-step counter and multiply/add order, so observing either
endpoint reproduces the ordinary owner bit for bit. Arbitrary admitted f64
physical times are observed through a canonical recursively bisected Gaussian
OU bridge. Shared ancestors provide the correct conditional cross-time covariance;
reordering observations, rejected CFD stages and smaller-step retries do not
change the represented path. There is no endpoint interpolation and no independent
per-query noise.

For left and right OU durations `dl, dr` and correlation time `tau`, define
`a = exp(-dl/tau)`, `b = exp(-dr/tau)`, `qL = -expm1(-2dl/tau)`,
`qR = -expm1(-2dr/tau)` and `qT = -expm1(-2(dl+dr)/tau)`. The conditional
midpoint mean has weights `a*qR/qT` and `b*qL/qT`, and its variance is
`(qL/qT)*qR`. Dividing before multiplication avoids premature variance underflow;
small positive innovations survive even when an interior decay rounds to one.
Rounded physical endpoint spacing is mapped affinely onto the original nominal
fixed dt. Spatial synthesis retains the requested physical time exactly. Claims
of temporal exactness concern this declared finite-mode Gaussian model and its
nominal clock, subject to floating-point rounding.

The independent `0x4d425247` kernel addresses each canonical midpoint using its
finite f64 magnitude bits, sign, mode and quadrature, with the retained interval
step mixed injectively into the bridge seed. Every normal consumes a disjoint
pair of counters within one interval. Existing endpoint and interval-average
streams remain unchanged. A source checkpoint plus the bridge seed and depth
admission reproduces subsequent interval observations; frames are observations,
not additional accepted source owners.

Preparation costs `4 * mode_count` normal draws. Each interior observation costs
exactly `4 * mode_count * tree_depth`; endpoint observations require none.
`required_normal_draws` permits the caller to charge stages and retries against a
single outer budget. The explicit maximum depth is at most 128. Out-of-interval,
nonfinite, nonrepresentable and deeper queries are refused without snapping.
One bounded ancestry and one coefficient table are reused across all modes;
consecutive equal correlation times reuse coefficients. Successful frame updates
and point, box, harmonic slip-plane and wall/box samples allocate no numerical
storage and use the original complete Mann tensor and analytic spatial gradients.

The immutable interval and reusable frame each admit nominal numerical storage
of 88 bytes per mode and validate actual capacities. This conservatively counts
the shared original temporal arrays once per owner, excludes immutable spectral
geometry and fixed metadata, and separately uses at most 128 stack ancestry nodes
and 128 coefficient triples. Frame updates validate the original definition and
bridge seed, poll cancellation within 64 normal draws and 16 coefficient nodes,
and publish amplitudes and time together only after the final poll. A refused or
cancelled update preserves its previously accepted frame and the original source.
Frames can be reused across successive intervals of the same source and seed.

Verification includes original endpoint and filtered-synthesis bit equality,
non-dyadic query-order/retry replay, interval reuse, an independent Gaussian
Schur-complement oracle for conditional means, variances and cross-time covariance
over 65,536 scalar realizations, slow/fast and rounded-clock limits, signed-time
addressing, depth/work/storage refusal and mid-operation/final cancellation.
These observations do not validate atmospheric intermittency, nonlinear energy
transfer or unresolved stresses, and OU forcing remains temporally unbounded in
frequency. The pre-existing interval-average API retains its independently keyed
integral bridge; this API does not establish a joint realization with those means.

## Exact represented spatial momentum and subfilter stress

`mann_field::filter::moments::MannBoxMoments` prepares an instantaneous quadratic
observation of the existing homogeneous Mann realization. Both temporal owners
expose `sample_filtered_moments`, returning the original box-averaged velocity
and gradient plus the complete central stress `average(u_i u_j) - mean(u_i)
mean(u_j)`. The result also exposes convective momentum flux per unit density and
the subfilter kinetic energy. These quantities retain represented fluctuations
that disappear when a coupling squares the filtered velocity alone.

For the same centered footprint gain `G(k) = product(sinc(k dot half_edge))`,
the cosine/cosine central coefficient for modes `m,n` is
`(G(k_m-k_n)+G(k_m+k_n))/2 - G(k_m)G(k_n)` and the sine/sine coefficient is
`(G(k_m-k_n)-G(k_m+k_n))/2`. Every cross-mode term is retained. This is the
spatial stress of one accepted realization, without an isotropic approximation
or an ensemble replacement. Boxes, skew parallelepipeds and collapsed faces and
lines use the same observation; an exact point returns exactly zero stress.

Preparation integrates small modal phases as central uniform moments and combines
the independent box axes with angle-addition identities. It retains the quadratic
sine and quartic cosine covariances without subtracting rounded unit gains. Under
footprint refinement the stress therefore recovers
`sum_edge (grad(u_i) dot edge)(grad(u_j) dot edge)/3`, including anisotropic and
cross-component stresses. The series is confined to preparation: stored pair
counts, storage budgets and the sampling arithmetic/cost remain unchanged.

Preparation admits `n(n+1)/2` unordered pairs and their storage explicitly: two
f64 coefficients per pair plus the original eight-byte gain per mode. Sampling
allocates no heap storage, reuses the accepted OU amplitudes and polls
cancellation throughout the quadratic work. It requires at most 12 KiB of stack
scratch. Prepared geometry remains bound to the complete original spectral
model. Neither sampling nor cancellation advances or mutates either time owner.

The stress is positive semidefinite in exact arithmetic; small floating-point
eigenvalues are not clipped. Cross-mode cancellation and underflow can still
limit accuracy when the physical stress itself is nearly singular or below the
representable range. Tests cover the analytic refinement limit down to 1e-10
footprints, the distinct quadratic/quartic modal limits, independent quadrature
across the series boundary, and edge reordering/orientation invariance.
This observation does not add energy outside the
represented finite spectral band, provide a nonlinear time average, or include
wall-blocking corrections, density fluctuations, pressure or a separate mean
flow. The existing linear filters retain their original O(n) cost.

## Density-stratified fetch mean wind

`fetch_mass_consistent::FetchAdjustedAnelastic` binds the existing fetch log law
to an immutable, caller-authored `ExponentialDensityProfile` on an explicit finite
height interval. It returns local physical density, velocity and `rho * velocity`
mass flux together. For inverse scale height `k >= 0`, its vertical correction
uses `(exp(k * (h - d - z0)) - 1) / k`, evaluated through `expm1` with the exact
constant-density limit. Thus `div(rho * velocity) = 0` on each smooth piece while
`div(velocity) = k * w` generally differs from zero. Sampling allocates no storage
and performs no quadrature; the reference speed scale is prepared once.

The original constant-density API, zero-velocity roughness cutoff and fetch clamps
remain unchanged. Cutoff and fetch endpoints are nonsmooth; endpoints retain the
interior fetch derivative convention. Nonfinite sources, out-of-domain heights,
unrepresentable density, velocity or mass flux refuse. This is an estimated
kinematic mean-flow model, without inferred gas properties, pressure/temperature
closure, turbulent density fluctuations or a momentum/energy solution. The
existing fetch battery checks independent finite-difference mass continuity for
both roughness-gradient signs, nonzero velocity divergence, constant/weak-density
limits, impermeability, density-scale invariance and domain/overflow refusals.

`Atmosphere::sample_anelastic_flow` and
`EvolvingAtmosphere::sample_anelastic_flow` additionally compose the original
horizontal log-law mean with the original frozen or OU-evolving modal
turbulence under that same explicitly prescribed density profile. If `v` is
the original solenoidal turbulence, the new perturbation is
`u' = (rho_ref/rho) v`. Thus `rho*u'` is solenoidal and the complete field
satisfies `div(rho*u) = 0`; the unmodified horizontal mean independently
satisfies that constraint. The exact gradient adds `k*u'_i` to each vertical
derivative after rescaling. Wall impermeability survives. Reference-height
velocity amplitudes and the constant-density field are retained exactly;
away from that height the velocity covariance scales with `(rho_ref/rho)^2`.

`AnelasticFlowSample` returns the actual background density, complete velocity
and gradient, and complete mass flux together. It does not reuse the constant
`AirScenario` thermodynamic values or assert an equation of state for the
authored profile. Querying allocates nothing and consumes no random draws or
clock advances. The existing fetch correction is not applied again. This is
an opt-in kinematic realization of the background mass constraint, as stated
by Bannon et al. (2006), <https://doi.org/10.1175/MWR3228.1>; it supplies no
stratified momentum/energy evolution, buoyancy, pressure response, density
fluctuations or empirical turbulence calibration. The dedicated anelastic
battery checks full frozen/evolving mass-flux divergence by independent
differences, all nine velocity derivatives, exact uniform-density recovery,
unchanged mean, wall condition, reference amplitude, density scaling, OU
noninterference/replay and atomic query refusals.

The full-tensor `HomogeneousMannField` and `FixedStepMannField` owners expose
`sample_anelastic` and `sample_wall_blocked_anelastic` under the same immutable
density profile. Their +z coordinate is the profile height, with a shared height
origin for wall observations. The original solenoidal perturbation is rescaled
once, after any harmonic wall blocking; the exact gradient includes the density
product rule. Uniform-density observations preserve the original bits, and
normal mass flux at the wall remains zero. No mean flow is added. These are
point observations: the original homogeneous velocity covariance is multiplied
by both local density-scale factors, and existing spatial-filter/moment APIs
are not reinterpreted as anelastic averages. Sampling preserves both clocks,
quadratures and checkpoints without allocation. The same fixed-background
kinematic limitations apply. The anelastic battery covers both clocks, independent
mass-flux differences, all gradients, wall flux, scaling, replay and refusals.

## Exact represented wall-blocked momentum and subfilter stress

`mann_field::wall::filter::moments::MannWallBoxMoments` extends the same
instantaneous quadratic observation to the harmonic slip-plane field. Both
temporal owners expose `sample_wall_blocked_filtered_moments`, returning the
existing exact filtered velocity and gradient together with the complete
spatial central stress, convective momentum flux per unit density and subfilter
kinetic energy. Every bulk/bulk, bulk/wall and wall/wall cross-mode product is
retained. The whole footprint must remain above the authored plane; skew
volumes, faces, lines and exact points share this observation.

The Fourier and evanescent wall functions have analytically integrated complex
exponential products on each independent uniform footprint coordinate. Wall
functions are normalized at the footprint's lowest point, so neither growing
exponentials nor underflow of a center correction destroys a deep wall-touching
average. Small exponents use centered polynomial moments directly; independent
axes are combined in central form. Tangential footprints combine each mode's
bulk and wall vector coefficients before quadratic accumulation, preserving
exact normal impermeability on the plane and tiny normal stress above it.
The original linear filter's mean and analytic-gradient arithmetic is retained.

Preparation accepts explicit retained-byte and unordered-function-pair budgets.
A nonpoint tangential footprint needs `n` functions, and a footprint with vertical
extent needs `2n`; `N(N+1)/2` pairs each retain four f64 covariance coefficients.
Storage also includes the existing linear filter's gains (8 or 32 bytes per
mode respectively). An exact point uses zero numerical storage and returns
exactly zero central stress. Preparation and sampling each use at most 24 KiB of
bounded modal stack scratch, in addition to small fixed intermediates. Sampling
allocates no heap storage, polls cancellation throughout the quadratic work,
draws no randomness and changes neither accepted amplitudes nor either clock.

Tests compare independent point quadrature and positive semidefinite directional
energies on wall-touching lines, faces and skew volumes; verify conservation of
momentum flux under footprint subdivision; and recover the analytic gradient
limit down to 1e-10 footprint scales. Manufactured-mode regressions retain a
deep wall layer after its center amplitude underflows and a near-wall quartic
normal stress around 2e-58. Both temporal owners, replay, geometry orientation,
exact point control, source binding, storage/work refusal and cancellation are
covered. All eight new tests and the complete current 62-test fs-atmo unit suite
passed using the pinned compiler and actual native dependencies. This integrates
the represented outer slip-plane field; it adds no unresolved spectral energy,
no-slip viscous layer, wall-generated turbulence, density fluctuations, pressure
or separate mean flow. Positive semidefiniteness holds in exact arithmetic;
roundoff eigenvalues are not clipped.

## Budgeted affine observations of represented Mann turbulence

`mann_field::affine::MannAffineField` is an immutable affine snapshot of one
accepted realization at one authored anchor. `HomogeneousMannField`,
`FixedStepMannField` and `MannOuBridgeFrame` expose `affine_observation`;
preparation uses the original value/analytic gradient, actual retained
quadratures and O(modes) work, without allocations or random draws. Each
subsequent sample costs O(1) and returns explicit Euclidean velocity and
Frobenius gradient truncation estimates. Homogeneous modes contribute
`A |k|^2`, where `A = sqrt(|a|^2 + |b|^2)`. A wall correction contributes
`R q^2` to velocity curvature and `sqrt(2) R q^2` to gradient curvature, with
`R = hypot(a_z,b_z)` and `q = |k_horizontal|`; horizontal phase and vertical
exponential derivatives are both included. Taylor remainder estimates are
`H_v r^2/2` and `H_g r`. These are analytical bounds evaluated in ordinary
floating-point arithmetic, not interval certificates or bounds on roundoff,
unresolved spectrum, time integration or trajectory error.

A wall snapshot requires its anchor exactly on the declared plane; the affine
field then retains no penetration there and the source's zero divergence up
to roundoff. Queries below the plane, invalid budgets, excessive spatial
remainder, or nonrepresentable nonzero bounds refuse. No automatic
reanchoring, query-order-dependent patches or silent exact fallback occurs.
Snapshots must be rebuilt after source time/amplitudes change; original exact
observations remain unchanged.

Five focused regression groups cover exact-anchor equivalence, realization
error estimates, wall-envelope derivatives and slip, reordered bridge stages,
representability, cancellation and zero fields. The actual atmospheric crate
and its unit-test modules passed a pinned compiler metadata check using real
dependencies; these new tests were not linked or executed.

## Compact advecting turbulent gust packets

`mann_field::wall::gust::MannWallGust::new(wall, center_m, half_width_m,
reference_time_s, strength)` defines a compact envelope on the original
wall-blocked Mann realization. Both `HomogeneousMannField::sample_wall_gust`
and `FixedStepMannField::sample_wall_gust` return its velocity and all nine
analytic spatial derivatives. Coordinates use the original upward Mann frame.
The center moves with the source model's uniform advection at its accepted
physical time. The same spectral polarizations, quadratures and OU histories
drive every point; observation introduces no new random process, clock,
checkpoint format or numerical allocation. Preserve the immutable packet
definition alongside the original source checkpoint for replay.

For wall-blocked velocity `v`, the fixed gauge is
`A = (integral_wall^z v_y ds, -integral_wall^z v_x ds, 0)`, with `curl(A)=v`.
The compact potential envelope is
`e = strength * product((1-s_i²)^3)` for every `|s_i|<1`, zero otherwise,
where `s_i` is the advected center-relative coordinate divided by its half
width. The returned velocity is `curl(e*A) = e*v + grad(e) cross A`.
Its gradient includes the full product rule, envelope Hessian and analytic
derivatives of `A`. Scaling only velocity would create spurious divergence
`grad(e) dot v`; the curl correction removes it. The tangential potential
vanishes at the wall, preserving exact zero normal velocity and its two
tangential derivatives even when the packet intersects the plane.

Vertical Fourier integrals use stable analytic sinc/expm1 expressions,
including zero vertical or horizontal wavenumber and thin near-wall limits.
The C2 potential envelope makes velocity and its first spatial derivatives
continuous across the compact support boundary. Sampling outside support
returns exact zeros before any spectral work. Inside support, cost is
O(modes), without quadrature or per-query allocation. Invalid geometry,
nonrepresentable reciprocal-square widths, below-plane queries, or
nonrepresentable evolved positions, velocities or gradients refuse without
changing the source. Strength is a nonnegative potential multiplier, not a
peak gust velocity or prescribed RMS. At the packet center only the velocity
is strength times the original wall velocity; the gradient still contains
the envelope curvature contribution.

The potential construction follows the divergence constraint used in
inhomogeneous synthetic turbulence (Yu & Bai, JCP 256 (2014), 234–253,
doi:10.1016/j.jcp.2013.08.055). This explicit gauge and compact envelope define
a prescribed inhomogeneous gust, not that paper's calibrated construction.
They change the stationary Mann spectrum/covariance; existing homogeneous
covariance and filter APIs are not reinterpreted as packet statistics. There
is no momentum/pressure solve, buoyancy, resolved transition or measured-gust
calibration. Native wind consumers can use the returned velocity/gradient
through their existing callbacks; no legacy flight forcing is changed.

Eight native regression groups cover a curl oracle built from independent
vertical integration of the public wall velocity, all nine gradients by
spatial differences, closed-surface mass flux without using analytic
gradients, exact wall/support zeros and edge smoothness, degenerate Fourier
limits, source advection and both OU checkpoint histories, strength scaling,
and invalid/unrepresentable queries. These test the prescribed kinematic
construction and its integration with the existing clocks.
