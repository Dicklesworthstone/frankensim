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

`tests/fit_battery.rs` additionally checks the finite Γ = 0 limit, seed/target
binding, input refusal, and the unchanged nonzero-shear fit golden. The
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
- No gust events, no thermals, no site heterogeneity (E3.3c territory).
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
