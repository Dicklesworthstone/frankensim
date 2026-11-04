# fs-aeroac — Contract

## Purpose and layer

Layer L3. Low-Mach aeroacoustic source models (bead
frankensim-fsim-aeroacoustic-sources-9ok02, IN PROGRESS — this crate
is the first slice): the acoustic-analogy half of the honest hybrid
over fs-lbm's incompressible base flows. That flow-derived family is a
SHAPE/SCALING authority; ABSOLUTE SPL IS NEVER CLAIMED (the crate's
scope law, pinned as data and by test).

### Physical three-dimensional compact volume sources

`monopole::CompactMonopoleObserver` is a separate physical-source observation:
`p'(r,t+r/c) = rho V''(t)/(4 pi r)`, using one original fluid-acoustic snapshot.
It returns pressure in pascals and its causal arrival time. Radius-dependent
compactness, observer distance and a conservative source-surface pressure/bulk-
modulus ratio must satisfy explicit limits. The declared maximum frequency must
cover the emitted source's nonlinear harmonics; a drive frequency alone is not
that guarantee. Nonfinite or unrepresentable values refuse without an output.
Preparation and observation allocate no memory.

The caller supplies physical volume acceleration. Spherical bubble dynamics can
provide it through `4 pi R (R R'' + 2 R'^2)`. Incompressible-flow pressure is not
silently reinterpreted as acoustic pressure. This free-field, compact, one-way
linear model supplies neither boundary scattering nor shock propagation,
attenuation, source radiation reaction or measured SPL qualification. Four G1/G3
tests check independent harmonic pressure, causal delay, radiative spreading,
material scaling and domain/representability refusals.

## Public types and semantics

- `bessel::{j0, j1, y0, y1, hankel0_outgoing, hankel1_outgoing}` —
  deterministic Bessel/Hankel functions built EXCLUSIVELY from
  fs_math primitives (det transcendentals, Dd double-double series
  accumulation below the x = 18 crossover, recurrence-GENERATED
  asymptotic sums above it — deliberately no transcribed coefficient
  tables). `e^{-i omega t}` convention: the outgoing cylindrical wave
  is `H^(1) = J + iY`. Y at non-positive arguments is NaN (real
  domain).
- `curle2d::dipole_pressure` — 2D frequency-domain Curle compact
  dipole over the outgoing Hankel Green's function
  `G = (i/4) H0^(1)(kr)`: `p = (ik/4) H1^(1)(kr) (rhat . F)` per unit
  span. The 2D-vs-3D Green's-function trap from the bead's polish
  round is structurally excluded (frequency domain, Hankel kernel);
  every output embeds [`SCOPE_STATEMENT`].
- `curle2d::PreparedDipoleObserver` prepares the actual outgoing Hankel
  propagation and analytic pressure-gradient transfer for multiple compact line
  sources in one `FluidAcousticState`. Complex peak forces retain one explicit
  frequency and common time origin. Pressure and particle velocity are summed
  before computing `p conj(u)/2`; active intensity includes interference terms
  and may have signed components, while reactive intensity is kept separately.
  The source layout need not be compact, but each constituent source still is.
  Preparation checks the caller's source-count budget; application is O(N),
  allocation-free, with six complex cached coefficients per source and
  double-double reductions. `MaterialJetLabiumRun::prepare_coherent_observer`
  binds the exact source medium and rejects frequencies at/above its Nyquist
  limit. Independent records must be physically phase-aligned by the caller.
  No scattering, mutual loading, absorption, or
  absolute-SPL validation is inferred. G0/G1/G3 `tests/interference.rs` checks
  cancellation, quadrature, path phase, near-field momentum, standing-wave
  intensity, rigid motion, dynamic-range residuals, and bounded admission;
  `tests/material.rs` exercises the real material-bound flow consumer.
- `PreparedDipoleObserver::evaluate_coherence_modes` propagates explicit source
  coherence factors: `F = B xi`, `E[xi xi*] = I`. The source second moments
  `B B*` are positive semidefinite by construction, without a dense source
  covariance matrix, phase synthesis or a post-hoc positivity repair. Each
  factor coherently combines all source positions and x/y forces before forming
  quadratic observations; only independent factors add powers. Returned
  `DipoleFieldStatistics` contains `E[|p|²]/2`, the Hermitian velocity spectral
  covariance `E[u_i conj(u_j)]/2`, and signed active/reactive intensity
  `E[p conj(u)]/2`. These are spectral-line powers, not PSD/Hz: one-sided PSD
  factors for a bin width df must be scaled by sqrt(2 df) into peak-force units.
  Broadband bin addition assumes stationary resolved frequencies, not arbitrary
  finite-time signals. Application is allocation-free O(N r), with explicit
  `CoherenceBudget` limits on factor count r and source-factor products N r,
  and double-double reductions both within and between factors.
  `evaluate_independent_sources` gives the equivalent O(N) path for mutually
  uncorrelated sources with coherent x/y force components, without N sparse
  N-vector factors. Independence and any supplied correlation are caller-declared
  physical assumptions; this API does not estimate a cross-spectrum or validate
  a turbulence model. G0/G1/G3 cases in `tests/interference.rs` compare the
  factor propagation with independent dense cross-spectral algebra (including
  complex correlation), coherent nulls versus independent power, standing waves,
  unitary factor changes, common-clock shifts, cancellation residuals and work
  limits. The material-run test exercises the independent-source observation
  through its actual source-medium binding. The existing 2D and one-way source
  authority boundaries remain in every statistical result.
- `curle3d::CompactForceObserver` prepares a stationary 3D point-force observer
  from one coherent `FluidAcousticState` and physical source/observer positions.
  For actual force ON the fluid [N], full pressure is
  `n.F(t-r/c)/(4 pi r²) + n.Fdot(t-r/c)/(4 pi c r)`; its far-field method
  explicitly omits induction. There is no extra density factor: material
  density belongs in source mechanics, while material sound speed remains in
  propagation. Forces must be 3D SI totals, never force per unit span. The
  prepared operator exposes `r/c`; callers evaluate retarded source history or
  use their existing delay filter exactly once. Evaluation is constant-work,
  allocation-free and uses fixed-order reductions. `tests/curle3d.rs` checks
  the independently differentiated retarded Green function, radial laws,
  signed directivity, rigid motion, density independence and refusals. This is
  an ideal free-space transfer for prescribed forces, not validation of a
  source model, finite-size compactness, enclosure scattering or measured SPL.
  `field_at_retarded` also accepts the original physical force impulse [N s]
  at that same retarded time, and returns all three particle-velocity components
  plus instantaneous signed `p u` energy flux. Its `1/r³` impulse and `1/r²`
  force tensors retain transverse motion even at a pressure node; the `1/r`
  term alone has the plane-wave relation. The caller owns the integration
  constant and complete source history; no block-local reset or inferred
  impulse is permitted. Cached coefficients and double-double reductions make
  evaluation constant-work and allocation-free, and keep small residual
  amplitudes when Cartesian forces cancel. Unrepresentable velocity/flux
  refuses without narrowing the pressure-only API. G1/G3 tests compare the
  full harmonic owner and cycle-mean intensity, the source-free momentum law,
  pressure-null transverse motion, signed amplitude and quadratic flux scaling,
  inverse-density response and nonfinite histories. The same compact linear
  one-way model limits remain; no source dynamics or radiation reaction is added.
- `curle3d::moving::MovingCompactObserver` prepares eight real pressure
  coefficients for physical Q/Qdot and force/force-rate at one emission event,
  retaining Doppler amplitude, source acceleration and near-field induction.
  Its constant-velocity trajectory owner analytically inverts retarded time
  under an explicit global subsonic margin, clearance and measured light-cone
  residual limit; supplied instantaneous kinematics establish no global motion
  guarantee. The original fluid snapshot and source clock are preserved, with
  no extra density on force. Evaluation is constant-work, allocation-free and
  coherently sums using double-double arithmetic. Source-history coverage and
  antialiasing belong to the consumer. This is prescribed compact 3D pressure
  in a quiescent lossless medium, without finite-size validation, scattering,
  absorption, feedback, transient force particle velocity or SPL qualification.
- `curle3d::moving::flow::UniformFlowCompactObserver` observes fixed laboratory
  sources and microphones in an explicitly supplied uniform subsonic fluid
  velocity. It reuses the original moving-source pressure coefficients after a
  Galilean transformation, retaining convective volume-source pressure, force
  direction and induction. The original fluid density and sound speed are
  unchanged. Direction-dependent travel time solves `|x-y-U T| = c T`; its
  downstream root remains finite near Mach one. Transient and harmonic queries
  use the same coefficients and original source clock, with the delay applied
  once. Both endpoints are stationary in the laboratory, so no spurious Doppler
  frequency shift or receiver amplitude multiplier is introduced. Preparation
  and observation allocate nothing and perform bounded constant work. Explicit
  Mach, geometric clearance and represented-time residual limits refuse invalid
  observations. `tests/uniform_flow.rs` checks the independent convected Green
  function, stationary-owner parity, axial/crosswind delay, steady injection,
  transient/harmonic agreement, rotation, material scaling and refusals. This
  is pressure in a homogeneous, lossless, uniform flow, without shear-layer
  refraction, boundaries, dispersion, particle velocity, energy-flux inference,
  source feedback, finite-source validation or measured-SPL qualification.
- `curle3d::moving::trajectory::PreparedSourceTrajectory` binds a finite C2
  quintic position history to the original microphone and fluid. Its position,
  velocity and acceleration come from the same accepted-knot interpolant.
  Outward `fs_ivl` Bernstein hulls admit speed and continuous receiver clearance
  over complete spans, with bounded subdivision, storage and cancellation;
  inconclusive bounds refuse. The resulting global subsonic bound makes the
  arrival map monotone and supplies a conservative Doppler-rate ceiling.
  Allocation-free safeguarded retarded queries reuse the original pressure law,
  preserve the requested receiver clock, and check the light-cone residual at
  the actual public f64 source time. Missing trajectory history never
  extrapolates. Bounds apply to this interpolant, not the source integrator's
  error or physical model; finite-source and receiver bandwidth admission remain
  separate. Existing uniform-motion queries retain their analytic fast path.
- `curle3d::moving::trajectory::PreparedKinematicTrajectory` reuses the same
  C2 quintic compiler and outward speed admission for independent source and
  microphone paths. `PreparedMovingReceiverSource` borrows two such paths in
  one original fluid state and admits a finite receiver window, continuous
  retarded clearance and source-history coverage. Bernstein span-pair bounds
  carry both physical time intervals, exclude causally unreachable events and
  subdivide under explicit pair/node/depth caps; inconclusive bounds refuse.
  Queries allocate nothing, preserve the exact receiver timestamp and measure
  residual at the actual represented emission time. Pressure is the original
  Eulerian source field at the moving microphone, without relative-velocity
  substitution or a receiver-Doppler amplitude multiplier. Separate metadata
  reports `(1-n.v_receiver/c)/(1-n.v_source/c)` and a conservative whole-window
  ceiling for bandwidth admission. Trial root positions may have zero distance;
  the physical clearance applies at admitted retarded events. Bounds concern
  prescribed interpolation, not mechanical error; microphone aerodynamic
  self-noise, directivity, finite-source validity and geometric-modulation
  bandwidth are not inferred. Source strengths remain on their original clock.
- `curle3d::PreparedCompactSourceObserver` prepares actual 3D source positions
  at one frequency and checked material state. `CompactSourcePhasor` explicitly
  separates volume injection Q [m³/s] from force ON the fluid F [N]. In
  `exp(-i omega t)` convention it applies the outgoing spherical Green function
  `G = exp(ikr)/(4 pi r)` with `p = -i omega rho G Q - F.grad(G)` and
  `u = grad(p)/(i omega rho)`. The full gradient includes angular derivatives:
  a dipole pressure null may still have tangential particle motion. Sources and
  mechanisms sum coherently before signed active/reactive intensity is formed.
  Sixteen cached complex coefficients per source give allocation-free O(N)
  application, with explicit preparation budget and double-double reductions.
  `CompactForceObserver::prepare_harmonic` preserves the same original source
  and observer coordinates, material and retarded phase. Force pressure does
  not acquire an extra density factor; monopole pressure does depend on density.
  `evaluate_coherence_modes` accepts explicit peak-unit factors of the joint
  Q/F source second moment, under `curle3d::CoherenceBudget` limits; its O(N r)
  application forms pressure variance, the full Hermitian 3D velocity covariance
  and active/reactive intensity. `evaluate_independent_sources` is the O(N)
  specialization, retaining coherence among Q/F components of each source.
  Neither path invents correlation or treats a spectral-line power as PSD/Hz.
  No periodic-flow force or N/m line source is automatically made an isolated
  3D radiator. These are prescribed compact sources with no finite-size
  validation, source motion, scattering, absorption or mutual loading.
  `tests/compact_interference.rs` checks both existing retarded owners,
  independently differentiated pressure, monopole/dipole spherical power,
  pressure-null velocity, separated-source and mixed-mechanism interference,
  standing waves, complex partial coherence against dense algebra, material
  density scaling, dynamic-range residuals and bounded admission.
- `PreparedCompactSourceObserver::with_plane` adds an `IdealAcousticPlane` to
  the original lossless or molecular-relaxation owner. The resulting
  `PreparedHalfSpaceObserver` solves one infinite rigid (`u.n=0`) or ideal
  pressure-release (`p=0`) plane by the exact image Green function. For inward
  unit normal n the image volume flow is `s Q` and image force is
  `s (I-2nn^T) F`, with s=+1/-1 respectively. It reflects dipole orientation as
  well as position. Both pressure and every near-field velocity component
  retain their direct/reflected path phase before intensity or variance is
  formed. Even mutually independent physical sources each interfere with
  their own reflected field; an image is never an independent noise source.
  The normal selects the fluid half-space; sources lie strictly inside it,
  receivers may lie on the plane, and unrepresentable image geometry refuses.
  Preparation bounds both paths per source and precomputes force reflection
  in the cached transfer. Applications are allocation-free O(N), or O(N r)
  for r explicitly authored coherence factors; the coherence work budget
  counts both paths. Cancellable preparation/application polls per source or
  factor and publishes no partial result. `evaluate_harmonic_line` feeds the
  existing finite-window and ensemble observations without a second delay.
  The original material, equilibrium and phase-delay distinction, source
  geometry, clock, and force/volume units survive. This is an ideal half-space
  boundary solution for prescribed point sources, not finite-source validation,
  a finite-impedance ground approximation, wall-edge diffraction, a room solver,
  source feedback or measured-SPL qualification. See Allen & Berkley, JASA 65
  (1979), 943–950, doi:10.1121/1.382599, for the image construction; the force
  reflection follows by differentiating its Green function in source position.
  `tests/compact_interference/plane.rs` checks both exact boundary conditions
  for all four source components, analytic reflection notches and delays,
  near-field momentum, rigid coordinate changes, complex partial coherence,
  molecular absorption/dispersion, finite-window beats, quiet cancellation
  residuals, cancellation and work/domain refusals. The `plane_interference`
  example prints physical pressure and energy flux through this same owner.
- `PreparedCompactSourceObserver::with_orthogonal_planes` extends that original
  transfer to one, two or three perpendicular ideal planes, including mixed
  rigid/pressure-release conditions. The full finite group has 2^m paths per
  source, including every multiple reflection. Preparation checks that physical
  sources and receiver share all selected half-spaces, refuses nonperpendicular
  normals (dot product above 64 machine epsilons), and bounds the full path count.
  It never snaps geometry. It folds all sixteen coefficients with compensated
  sums into each original source transfer, retaining O(N) allocation-free
  application and original source coherence, gas, geometry and source clock.
  Only O(N) image workspace is retained during preparation; the final cache does
  not grow by 2^m. Direct-path delay queries remain direct-path quantities, not
  an invented delay for the combined field. This is the finite ideal corner
  Green function, with no general-wedge or enclosed-room reverberation claim.
  `tests/compact_interference/corner.rs` checks independent image pressure and
  momentum, all rigid/release wall conditions, edges/vertices, rotated frames,
  source statistics, molecular propagation and bounded/cancelled admission.
- `observation::PreparedHarmonicField` observes deterministic lines on a common
  receiver clock over an explicit rectangular interval. The existing 2D/3D
  prepared observers export `evaluate_harmonic_line` with the propagation delay
  already applied. Exactly equal frequencies coalesce using double-double sums;
  close frequencies remain distinct. Pressure mean, mean square and physical
  `<p(t) u(t)>` retain both difference-frequency beats and sum-frequency terms,
  including subperiod windows. The shared `fs-fft::harmonic_window` integral is
  reused per line pair across all velocity components. Preparation and each
  allocation-free observation are bounded by explicit line/pair budgets. This
  is floating-point deterministic-signal integration, not inferred stochastic
  coherence, a certified enclosure or a new source/SPL authority. Severe
  cancellation between distinct frequencies remains roundoff-limited; negative
  computed pressure squares refuse rather than being clamped. G0/G1/G3 cases
  in `tests/harmonic_observation.rs` check a subperiod analytic solution,
  coherent nulls and small residuals, beats against independent real-signal
  quadrature, clock shifts, actual 2D/3D source-consumer paths and admission.
- `observation::PreparedProperHarmonicField` adds finite-window second moments
  for explicitly authored joint source/frequency factors. For factor fields
  `A_r(t) = sum_j A_rj exp(-i omega_rj t)`, the real field is
  `Re sum_r A_r(t) xi_r`, with declared moments `E[xi_r]=0`,
  `E[xi_r conj(xi_s)]=delta_rs` and `E[xi_r xi_s]=0`. These assumptions give
  positive-semidefinite covariance by construction; they do not require a
  Gaussian distribution, imply stationarity, or infer unknown source phases.
  Properness removes sum-frequency terms even in subperiod windows. Authored
  cross-frequency coherence still produces difference-frequency beating;
  only distinct uncorrelated factors add second moments without cross terms.
  Same-frequency source contributions coalesce before squaring within each
  factor. The original receiver, frame, source clock, material propagation and
  source-model authority remain the caller's responsibility, with propagation
  supplied by the existing observer owners. Peak-unit coefficients use
  covariance `2 df S` for a one-sided PSD bin of width df, not unscaled PSD/Hz.
  `WindowedAcousticStatistics` returns the ensemble/window mean pressure square,
  full real symmetric velocity covariance and signed physical `E[p(t)u(t)]`
  energy flux. These are averages of instantaneous products, not products of
  window averages, complex spectral covariances, reactive intensities or
  uncertainty estimates. A separately observed deterministic mean field adds
  its pressure square and intensity to the zero-mean fluctuation contributions;
  nonzero pseudocovariance is outside this model. Preparation checks total
  modes, supplied lines and sum of squared per-mode line counts before any
  allocation. Observation allocates nothing, reuses one proper temporal kernel
  per frequency pair, and accumulates every quadratic contribution in
  double-double arithmetic until the complete bank is reduced. Physical
  frequency differences are formed in Hz before radian conversion, preserving
  a slow beat in the presence of a much larger common carrier. Negative
  pressure/velocity variances refuse without clamping; arbitrary nulls have no
  relative-accuracy guarantee. Callback cancellation returns typed refusal
  without a partial result and leaves the immutable bank reusable. G0/G1/G3
  cases in `tests/proper_harmonic_observation.rs` check subperiod normalization,
  complex partial correlation against independently synthesized four-phase
  ensembles, beating versus declared independence, pressure-null motion,
  dynamic-range residuals, unitary factors, physical clock and common-carrier
  shifts, sub-ulp beat energy, bounded admission and cancellation/reuse.
- `PreparedCompactSourceObserver::new_molecular` retains an original
  `fs_material::relaxation::RelaxationState` and propagates its relaxing
  compression law through the same coherent 3D source operator. The original
  molecular owner supplies `k² = rho omega²/K(omega)`; pressure, radial/angular
  gradients, particle velocity and every spectral second moment use that same
  complex `k`. The momentum denominator stays `i omega rho`, so this is
  inviscid momentum with molecular bulk absorption and dispersion. Shear
  viscosity, thermal conduction and wall losses are explicitly absent; no
  plane-wave attenuation coefficient is reinterpreted as a bulk modulus.
  Source and wave-number accessors retain the actual material and frequency.
  `source_delay_s` remains the equilibrium `r/c` reference; the new
  `source_phase_delay_s` exposes `r Re(k)/omega`, which is neither group delay
  nor a broadband signal-front time. Application retains the same O(N),
  allocation-free cached transfer and O(N r) coherence paths, with a distinct
  model scope carried by coherent and statistical outputs. Original material
  errors survive through `MaterialAeroacError`; nonfinite or underflowed paths
  refuse instead of silently declaring physical silence. G1/G3 additions in
  `tests/compact_interference.rs` recover the empty-reservoir lossless limit,
  compare a one-reservoir gas with independent caloric elimination, check a
  three-rung pressure-gradient ladder and the exact volumetric absorption /
  spherical-flux balance, and retain coherence and admission behavior. These
  are numerical model checks, not experimental gas or source qualification.
- `bickley::bickley_rayleigh_mode` — Rayleigh-equation shooting
  solver (RK4 + complex secant, half-line with symmetry boundary
  conditions and the exact `phi' = -alpha phi` far-field decay) for
  the Bickley jet `U = sech^2(y)`: the inviscid instability oracle
  the fs-lbm jet validation chain runs against. Non-convergence
  REFUSES with the residual disclosed; no partial eigenvalue.
- `bickley::rayleigh_residual_closed_form` — the self-verification
  surface for the analytic pins.
- `jetlab::run_jet_labium` — the edge-tone base-flow fixture: an
  fs-lbm slot jet (smoothed top-hat) impinging on a two-cell splitter
  plate (DISCLOSED simplification of a wedge — the sharp edge, not
  the wedge angle, drives the oscillation class), in a PERIODIC-x
  domain closed by a per-row-profile FRINGE layer (the spectral-DNS
  fringe method: the sponge re-conditions outflow to the authored
  inflow profile, and doubles as the measured acoustic absorber).
  Records the plate force series via momentum exchange; returns
  per-run diagnostics (max Mach, plate-vs-fringe plane fluxes,
  Reynolds) and embeds the scope statement.
  `jetlab::dipole_spectrum_line` radiates a caller-FFT'd force line
  through the 2D Curle dipole.
- `jetlab::transverse_force_peak` — the shared Hann-periodogram peak
  instrument (mean removal, guard-band skip, prominence vs median
  admitted power, unwindowed RMS floor value); reports the peak BIN
  so the +-1/(2 bin) quantization travels with every Strouhal claim.
- Lattice destabilization is a TYPED refusal: a blow-up inside
  fs-lbm's stepping (its density assert panics on the first bad
  cell, so no pre-scan can precede it) is contained at the
  settle/ramp/record boundary via `catch_unwind` and returned as a
  named `NonFinite` refusal; the rig state is discarded on that
  path. Executed trigger: plain-BGK destabilizes at cell Reynolds
  `Re/delta` in (36, 48) on this rig family (measured at the fixed
  rig speed u = 0.08 across delta 6 and 12; the Mach/u dependence of
  the boundary was not probed and is not claimed).
- `jetlab::run_adiabatic_ramp` — the hysteresis-following protocol
  for the MULTI-STABLE rig: after one base lock-in, Reynolds is
  swept quasi-statically (per-step linear relaxation-time ramp at
  FIXED jet speed, so Mach, fringe profile, nozzle slit, and seed
  are all unchanged) up to `reynolds_end` and back, with a spectral
  rung measured after every transition. Returns `2 rungs - 1`
  direction-tagged measurements; adiabaticity (ramp steps >> the
  oscillation period) is the caller's declared obligation. Refuses
  degenerate spans and a top-of-ramp tau at or below the 0.5005
  stability floor.
- `SCOPE_STATEMENT` — the no-absolute-SPL / 2D-to-3D-span-correction
  law as data (the marketing-mutation guard asserts it).
- `regime` — bead l011o: pinned CentralMoment 2D tonal-lock rows
  (`PINNED_2D_CENTRAL_MOMENT_TONAL`), the named inverse-cascade mechanism,
  `admit_broadband_spectrum` (refuses values outside `[0, 1]` and at or below
  `TONAL_FLATNESS_CEILING`), `measure_spectral_flatness` /
  `classify_spectrum` (geometric/arithmetic power-mean over all supplied bins;
  zero bins give zero flatness in a nonzero band, negative/nonfinite power and
  a silent band refuse, and no noise floor is inserted; no table is
  minted), and `evaluate_slot_jet_3d_operator` (D3Q19 CentralMoment
  smoke). `SlotJet3dOperatorSmoke::incorporate_measured_spectrum`
  raises `broadband_demonstrated` only from a real spectrum that
  clears the floor. The smoke alone does not mint a 3-D broadband
  table.
- `slot_jet_3d` — bead frankensim-music-v8-root-3ez8g.10.1: the 3-D
  slot-jet rig on the D3Q19 central-moment operator, the executed
  follow-up lane to the recorded 2-D tonal lock. Periodic domain +
  right-side per-y-row fringe (`Fringe3`, the Sponge2 port; its 3-D
  reflection coefficient is MEASURED by the test battery with a
  bounce-back-wall control), binary-slit nozzle wall (the executed
  2-D lesson kept), two-cell splitter plate across the full span,
  momentum-exchange plate force, `run_slot_jet_3d` /
  `classify_rung` / `SlotJet3dRung::to_jsonl` receipts with the FFT
  bin disclosure and the `FORCE_RMS_AMPLITUDE_FLOOR` amplitude law.
  Chunked checkpoints fail closed on unknown format versions, files beyond the
  configured size envelope, inconsistent settle/record/force counters, any
  truncated or trailing state bytes, non-finite force/population values, and
  non-positive fluid density before resumed stepping can mutate the rig. Their
  identity covers the full run configuration, including settle and record
  horizons; terminal checkpoints are retained for idempotent replay.
  The Re actuator is the central-moment second-order rate at fixed
  `u_jet`; independent rungs, not adiabatic continuation. Both
  outcomes are honest wins: a broadband rung above the tonal
  ceiling or a quantified tonal refusal. Driver:
  `src/bin/slot_jet_3d_sweep.rs` +
  `scripts/e2e/music/slot_jet_3d_sweep.sh`. No experimental,
  video-backed, or absolute-level claim survives this lane; the
  crate scope statement applies to every spectral quantity.
- `jetcard::{JetCard, JetCardClaim, CardAuthority, mint_tonal_interim_card,
  mint_broadband_card, mint_refusal_boundary_card}` — minted jet
  cards (bead 3ez8g.10.2), the typed boundary object between this
  lab lane and the audio-rate performance lane under the downstream
  "no card -> refuse" menu law. Authority law enforced at minting
  and revalidated on parse: `X-Struct` requires a residual issued
  against lab data, `X-Est` otherwise; a broadband claim requires
  16-band noise content; every card embeds the scope statement and
  refuses queries outside its demonstrated Reynolds range. The
  serialization is a pinned-field-order v1 JSON document with a
  strict fail-closed parser (foreign schema/claim kinds, tampered
  scope, trailing bytes all refuse by name); `content_hash` is
  FNV-1a over those canonical bytes — a deterministic correlation
  identity, NOT a cryptographic or authenticity anchor (the
  workspace `ProvenanceHash` caveat applies verbatim). The tonal
  interim card is minted from the recorded stage-I staging point
  (St 0.03662 vs Brown 0.03554, +3.0% inside the ±6% bin; ramp
  hysteresis and 0.0366/0.0458 multi-stability recorded) with the
  exact staging-rig FNV fingerprint as provenance.
- `forced_jet::{ForcedJetResolution, ForcedJet, JetStepResponse,
  bickley_step_response, FRONT_MARGIN}` — beads pexfp.135/.138/.321,
  ruling 7376: the causal response of a planar Bickley jet to a uniform
  transverse acoustic velocity at a baffled flue exit, integrated in
  time from the linearised vorticity equation (Thom wall at the exit:
  the jet leaves undisplaced and parallel; the exact discrete harmonic
  exterior beyond the lateral band; a downstream wall kept
  `FRONT_MARGIN` beyond the disturbance front). Units: `b`, `U0`,
  `b / U0`, `Re = U0 b / nu`. `bickley_step_response` records the
  centreline displacement at one station for a unit velocity step from
  rest over the declared window and holds it afterwards; `transfer(w)`
  is that held response's exact Fourier transform and `hat_weights`
  the audio-rate convolution weights in seconds for a mouth velocity
  linear between samples (`g_0` multiplies the current sample, the
  weights end one sample after the window and sum to the held gain).
  `leading_hat_weights` returns their count and only the first `count`
  of them, bit for bit, for a slow jet whose caller reads only its
  leading lags.
  `interval_hat_weights` integrates this SAME interpolant's derivative
  against both physical-time endpoint hats on a supplied convective interval.
  It returns newer/older weights in seconds, leaves the held tail inactive,
  and refuses invalid ages, advances or durations. At uniform spacing,
  adjacent interval weights recompose the existing hat weights. The consumer
  owns the accepted clock history; this method adds no variable-Re flow law.
  Refuses non-power-of-two or short boxes (`length < window + 8`), an
  explicit viscous number `4 dt / (Re dx^2)` above 2, stations outside
  the box, and non-finite results.

## Invariants

### Shared material integration

`material::MaterialJetLabium` binds geometry and physical speed/dx/dt inputs
to one `FluidTransport + FluidAcoustics` source and the existing
`fs-lbm::material::MaterialScaling` owner. Gas coherence and each liquid or
real-gas owner's original state/phase admission remain active. Acoustic and
transport density must match exactly before compilation. The existing jet-labium solver
receives material-derived viscosity; its actual surface force is converted
to N/m with `rho * dx^3 / dt^2`. Physical frequency and the same fluid's sound
speed supply `k = 2 pi f / c` to the existing Curle observer. No independent
rho, viscosity, Reynolds number or acoustic speed can override the binding.
Pure, humid, mixed and thermochemical gases, nonideal gas states and acoustic
liquids use the same path. Original source owners remain borrowed and available
through `material()`/`scaling().material()`. The frozen `acoustic()` snapshot
supplies both physical Mach and propagation; observers never re-query a source
or reconstruct liquid/real-gas properties through an ideal-gas surrogate.
`MaterialJetLabiumRun::observe_spectral_line` transforms its own accepted
force history with `fs-fft`, mean removal and the existing symmetric Hann
window. Coherent-gain normalization preserves resolved sinusoid amplitude;
conjugation makes the phase convention `exp(-i omega t)` relative to the first
recorded sample. Physical bin frequency and same-medium pressure are returned
together. An independent direct DFT checks force amplitude and phase.

Both material-run types provide `prepare_force_spectrum(max_samples)` to retain
the two force-component FFTs once. The sample budget is checked before O(N)
storage and O(N log N) transformation. Cached interior-bin queries do not repeat
an FFT or allocate its workspace; they preserve the existing mean, symmetric
Hann, coherent-gain and phase arithmetic. The 2-D `JetForceSpectrum` observes
N/m sources in its original fluid. `SlotJet3dForceSpectrum` retains actual N,
the original periodic span and the same-fluid acoustic wavenumber. DC and
Nyquist bins remain refused, and no isolated 3-D radiator is inferred.

For a known fundamental that falls between FFT bins, both runs also expose
`fit_force_harmonics(&fs_fft::harmonic::HarmonicFitPlan)`. The shared L1 owner
jointly fits DC and every declared sine/cosine pair through one prepared QR;
the two physical force axes reuse that factorization and their work buffers.
The plan must match the original sample count and timestep exactly. Its memory
budget and reciprocal-condition gate remain explicit, with original typed
numerical refusals preserved as `AeroacError::HarmonicFit`. The adapter adds one
reused N-sample axis buffer to the fitter's scratch/output storage. Returned
peak phasors use exp(-i omega t), phase zero at the first recorded sample, and
retain N/m or N without a span conversion. `JetForceHarmonics::observe_harmonic`
uses the fitted physical frequency in the original Curle transfer;
`SlotJet3dForceHarmonics::force_harmonic` retains periodic-cell source semantics.
The fitted DC and per-axis sample residual RMS are exposed separately. This is
known-frequency estimation, not pitch detection, a stationarity certificate,
statistical uncertainty, broadband PSD, recovery of already aliased samples,
or source/radiation validation. Unexplained sampled content remains in the
residual; unmodeled tones and frequency error can also bias fitted components.
G0/G1/G3 tests check independent Hann DFT coefficients, weak
off-bin even harmonics, phase/units, same-medium radiation, clock mismatches and
bounded admission; the actual 2-D/3-D flow tests exercise both consumer routes.

For independent or inharmonic lines, `fit_force_tones(&ToneFitPlan)` uses the
same two-axis buffers and L1 QR/scaling owner. Both physical run types accept
strictly increasing caller-supplied frequencies and preserve their exact
sample clock. `JetForceTones::observe_tone(index, ...)` applies the original
same-fluid 2-D Curle transfer at that declared frequency;
`SlotJet3dForceTones::force_tone(index)` returns the integrated N and acoustic
wavenumber without converting the periodic span into an isolated radiator.
Indices are zero-based in the plan's exposed `frequencies_hz()` order. There
is no nearest-bin, integer-order or tolerance-based frequency substitution.
The result retains joint DC, complex peak phase and per-axis residual RMS.
These are known-tone estimates, not source separation from unknown pitches,
an inharmonicity estimator, or a replacement for a broadband spectrum. G0/G1/G3
checks resolve a weak line beside a stronger tone within one FFT bin, retain
inharmonic radiation frequencies and periodic-span force units, compare the
same frequencies with the harmonic API, and preserve typed nonfinite/clock
refusals. No new force waveform is generated by either analysis path.

`fit_force_harmonics_refined(seed_hz, orders, frequency_budget)` now estimates
one local fundamental jointly from both recorded force components. The shared
L1 vector refiner preserves their relative physical amplitudes, phase at sample
zero and independent DC terms; constant axes add no frequency information.
The original timestep and frozen fluid are retained. The 2-D result continues
through the existing same-medium Curle transfer in N/m at the refined harmonic
frequencies; the periodic 3-D result retains integrated N and the actual span,
without inventing a missing force component or isolated-source radiation.
The explicit total workspace budget includes the transposed force record and
returned force/order arrays. This is a bounded stationary local fit, not global
pitch detection, a source-compactness test or an alias correction. Diagnostics
report evaluations and remaining local correction, not frequency uncertainty.
Focused checks cover a missing fundamental, weak cross-axis harmonics,
physical phase and radiation, periodic-span units and total-budget refusals.

`fit_force_tones_refined(seeds_hz, bounds_hz, frequency_budget)` extends both
material-bound runs to several locally unknown, independent frequencies. Each
tone retains a positive, ordered, nonoverlapping sub-Nyquist interval; no
harmonic-spacing assumption is imposed on beating or inharmonic forces. One
L1 variable-projection fit uses all recorded components with their original
physical relative strengths. The frequencies, peak phasors and DC are solved
jointly, with phase at the original first sample. The existing force-tone
observers use those inferred frequencies and the original fluid for propagation.
The 3-D result keeps its integrated N and periodic span. The complete storage
budget includes transposed samples and assembled force arrays, while returned
frequency and correction arrays are moved without duplication. Diagnostics
expose candidate count, local corrections and joint frequency conditioning;
none is an uncertainty bound. Missing or inseparable requested tones refuse,
as do exhausted work and active bounds without convergence. This does not
discover the number of sources, establish stationarity, or undo sample aliasing.

The physical Mach envelope is strictly below 0.3 and is checked separately
from numerical lattice Mach. Grid allocation is budgeted before construction.
An observed speed-envelope violation refuses output; diagnostics sample the
existing two planes every 16 steps and do not prove the full-field bound.
Observer frequencies at or above the force record's Nyquist limit refuse.
The medium is homogeneous and frozen during a run. Compressible feedback,
reaction, cavitation, phase transitions, relaxation/absorption, global flow validity and quantitative
experimental accuracy are not established by this adapter; the crate's
two-dimensional shape/scaling scope remains unchanged.

`tests/material.rs` verifies the actual lattice-to-SI force/frequency mapping,
same-state radiation, independent physical/numerical Mach admission, and
data-only pressure/temperature/mixture changes through executed flow and
observer outputs (G0/G3). These short transient tests establish integration,
not a converged tonal spectrum or absolute sound level.
The original IAPWS IF97 water EOS and industrial transport laws also drive
executed jet histories, SI forces and same-state radiation at two temperatures.
Independent force-unit and wavenumber checks cover this liquid path; incompatible
transport/acoustic densities refuse before either spatial rig runs.

`material::MaterialSlotJet3d` applies the same material/scaling owner to the
existing D3Q19 central-moment slot-jet rig. Its second-order relaxation rate is
derived from the fluid's actual kinematic viscosity; only higher-order numerical
relaxation remains an independent control. Three-dimensional momentum exchange
converts with `rho * dx^4 / dt^2` into N across one explicitly retained periodic
span. Both material adapters share one Hann/coherent-gain/phase FFT extraction.
The physical line frequency and acoustic wavenumber use the same timestep and
fluid as the flow. A three-dimensional periodic-cell force is not passed to the
two-dimensional Curle Green function or relabeled as microphone pressure.

The 3-D initial-speed admission includes the rig's 1.5 maximum spanwise seed
factor. Final-plane speed diagnostics use all three velocity components
(`momentum / density`), and plane mass flux is momentum itself. These are
final-plane samples, not a global speed bound or maximum over the recording.
`tests/material_slot_jet_3d.rs` executes pressure, temperature and composition
substitutions through the original solver, independently checks SI force and
DFT amplitude/phase, and checks necessary material/Mach/allocation refusals.
An original IF97 liquid also reaches the native D3Q19 solver, physical force
conversion, Mach report and spectral wavenumber without replacing its EOS.
This frozen homogeneous low-Mach adapter supplies no compressible feedback,
chemical evolution, acoustic radiation, turbulence-convergence or SPL claim.

1. EXACT Wronskian: `J1 Y0 - J0 Y1 = 2/(pi x)` holds within 5e-13
   over a log grid 1e-3..3e3 including the series/asymptotic
   crossover from both sides. Honest limits (review-measured against
   a 120-digit reference): the identity is blind to COMMON phase
   errors, and pointwise accuracy beyond the certified grid degrades
   as ~x*eps from f64 phase representation (~2e-11 at x = 1e6) —
   stated in the module doc; intended kr ranges sit inside the
   certified band.
2. Independent-path derivatives: central-difference `J0'` matches
   `-J1` and `Y0'` matches `-Y1` at FD accuracy (1e-8) across both
   regimes.
3. Cross-implementation oracle: fsci-special's cephes-heritage
   j0/j1/y0/y1 agree within 5e-12 of envelope (fsci's plain-f64
   Y-series carry ~1e-12 cancellation near their own crossover —
   measured; the Wronskian is the precision arbiter).
4. Small-argument limits exact (J0(0) = 1, J1 ~ x/2, bounded Y0 log
   constant); far-field `|H0| -> sqrt(2/(pi x))` within the O(x^-2)
   envelope correction.
5. Curle dipole physics: far-field decay exponent is -1/2
   (CYLINDRICAL spreading; the 3D-Green's-swap mutation reads -1 and
   fails the gate), the perpendicular directivity null is EXACT, the
   45-degree amplitude matches cos(theta) within 1e-3, and the scope
   statement rides on every output.
6. Bickley analytic pins are SELF-VERIFIED per run: `phi = sech^2` at
   `(alpha, c) = (2, 2/3)` and `phi = sech tanh` at `(1, 2/3)` drive
   the closed-form Rayleigh residual to machine zero (< 1e-14) at
   every probe point, and the residual is proven LIVE (perturbing
   alpha or c by 0.01-0.05 moves it above 1e-3). No literature value
   is transcribed on trust.
7. The shooting solver reproduces the pins from the unstable side
   (alpha = 1.95 sinuous: c within 0.02 of 2/3 with small positive
   Im c; alpha = 0.95 varicose: within 0.03), shows positive
   mid-band growth exceeding the near-neutral growth, and the
   eigenvalue is grid-converged (2048 vs 4096 RK4 steps within
   1e-8).
8. VERIFICATION CHAIN link 1 (dev-dep composition with fs-lbm): the
   D2Q9 Bickley jet's measured linear growth rate CONVERGES toward
   the inviscid Rayleigh oracle under simultaneous grid + Reynolds
   refinement — measured -21.4% (b = 8 lu, Re 240) closing to -12.0%
   (b = 12 lu, Re 432), both biased LOW as finite-Re physics demands,
   with clean-exponential window consistency, low-Mach and linearity
   diagnostics green per run, and gates authored just outside the
   measured values (a fixed loose envelope would hide a broken
   convergence direction). Executed lesson in the test doc: the free
   jet's own viscous diffusion (~nu/b^2 per step) decays the base
   flow and drags the measured rate down over long runs — fit early,
   at high Re.
9. Jet-labium fixture: THE VACUOUS-OSCILLATION TRAP (executed and
   now gated): the mirror-symmetric rig, unseeded, preserves
   symmetry to roundoff and its force spectrum shows high-prominence
   structure in ~1e-15 amplified machine noise — prominence-style
   ratios cannot detect it, so every amplitude-bearing test now
   REQUIRES the deterministic sinuous seed (`seed_amplitude`) and a
   force-RMS floor (1e-6 lattice). Seeded, the default configuration
   (Re 240, edge at 3 slot heights) saturates into a REAL limit
   cycle (Fy rms >> floor, essentially pure tone), max Mach 0.07,
   flux imbalance under 1% (the outlet-reflection pathology read
   6%). Geometry/regime refusals typed; bitwise deterministic.
10. EDGE-TONE STAGING (ignored heavy test, executed on record;
    seed-provenance trio in the test doc): with a nozzle wall at the
    jet root and the 0.005 sinuous seed, the SATURATED rig at
    Re 144, h/delta = 10 locks stage-I St = 0.03662 vs Brown's
    (1937) 0.03554 — +3.0%, within the record's +-6% bin
    quantization and INSIDE the published spread (Vaik/Paal exp
    0.03723, CFD 0.03497; two fetched sources). The unseeded run
    selects the SAME frequency in amplified roundoff
    (frequency-selection only); a 0.02 seed lands a NEIGHBORING
    locked state (St 0.0458) — edge-tone multi-stability, recorded.
    STRUCTURAL FINDING: without the nozzle the rig locks to the free
    jet's own mode (St 0.46, 12x the ladder) — the jet-root
    receptivity edge closes Brown's loop. SLIT-SHOULDER FIX
    (executed): with a nozzle the fringe target is the BINARY slit
    profile matching the wall opening (a smooth target slammed into
    the shoulder every wrap and fed a slit-lip mode at St_delta 2-3
    that blocked lattice refinement); the fine-lattice regression
    test (delta = 7.5, geometric similarity, recorded St 0.0458 =
    1.29x Brown) gates the LADDER BAND [0.7, 1.4] x Brown.
    MULTI-STABILITY (executed): attractor selection varies with seed
    and resolution (recorded states 0.0366 and 0.0458), so strict
    cross-resolution convergence of the selected attractor is not
    claimable from fixed-parameter runs; the adiabatic ramp
    (invariant 13) is the follow-up instrument.
13. ADIABATIC RAMP PROTOCOL (ignored heavy test, executed on record
    at TWO ramp rates — 15 and 3 stage-I periods per transition —
    with qualitative agreement): branch-following over
    Re 144 -> 264 -> 144 at h/delta = 10 stays inside the stage-I
    NEIGHBORHOOD [0.6, 1.7] x Brown on every rung of both legs
    (executed extremes 0.717/1.546; the slit-lip >5x, free-jet >12x,
    and stage-II ~2.3x signatures never appear), wanders between the
    neighboring 0.0366/0.0458 locked families (mode competition),
    and returns to Re 144 on a DIFFERENT attractor than it started
    from (bin 9 vs bin 7 on the committed protocol) — REAL,
    bitwise-reproducible hysteresis. NO stage-II lock occurs by
    Re 264 despite the literature onset ~220-250 for a wedge: a
    recorded NO-CLAIM boundary of this two-cell-plate periodic-
    fringe 2D rig, not a validated staging result.
11. NOISE TABLES (the product deliverable): `noisetable` fits
    Strouhal-band power-DENSITY shapes (record-length-independent —
    a band-SUM convention broke the synth round trip by the log-band
    width ratio, executed) over a velocity sweep with per-entry
    regime gates and real-amplitude floors; the JSON export embeds
    the scope statement (marketing-mutation guard) and the MEASURED
    source-force mean-square exponent as data (the saturated tonal limit cycle is NOT
    monotone in u near mode switches — multi-stability, reported not
    prescribed); the demo synth consumes the table and its output's
    band densities match within 6 dB (measured worst 2.5 dB).
    Band-limited force RMS integrates the one-sided periodogram BEFORE
    dividing each band's power by its bin count. It uses `N * sum(w²)`
    normalization for the symmetric Hann window; interior bins have their
    conjugate partners, while Nyquist is counted once and DC is excluded.
    This is a windowed finite-record estimate with leakage, not an isolated
    tone's coherent-gain correction or a radiated-power estimate. Relative
    density shapes still normalize their strongest band to 0 dB. The
    dimensionless shape-only synth filters uniform [-0.5, 0.5] input noise
    through the normalized inverse FFT, without a duration-dependent gain
    or implicit conversion to force or microphone pressure.
12. Refusals typed by name (non-positive wavenumber, coincident
   observer/source, non-finite inputs, bad solver parameters,
   non-convergence with residual); Hankel functions outside their
   domain are FULLY NaN (both components — a half-valid complex once
   leaked, review-caught); the Rayleigh solver's result is
   documented GUESS-DEPENDENT (any eigenvalue of the symmetry class
   zeroes the mismatch); determinism bitwise.
14. FORCED JET (`tests/forced_jet.rs`, executed on trj lane q): at
   station 11 b and Re 600, the step kernel's transfer equals
   independent harmonic forcing from rest within 0.1 % at w 0.5. At
   w 1-2, holding the response after window 90 cuts a slowly decaying
   tail. That costs 2.5-3.9 % (1.2 % at w 3; the test holds w 2 to
   4 %). At w 2, window 150 halves it, while the harmonic side moves
   1e-3, and spacing and time step leave it at 3.2-3.8 %;
   the kernel converges at second order on the ladder dx 0.2/0.1/0.05
   (window 90 b/U0; change 0.1 -> 0.05 is 0.07-0.11 % at w 0.03-0.16,
   0.7-1.9 % at w 0.5-1, held gain 14.8625 -> 14.8663 b/U0) and is bit
   identical at 1 and 5 workers; doubling the box (102.4 -> 204.8 b)
   or narrowing the band (8 -> 6 b) moves it at most 3.3e-3 / 1.9e-3
   at w >= 0.06. Inviscid, the downstream wave carries the published
   sinuous spatial mode at w 0.8 (Sabatini & Bailly 2015,
   1.449709 - 0.134109i): wavenumber within 0.1 %, growth rate
   converging from below with the critical-layer resolution (4.4 % low
   at dx 0.1 with window 125; a dx 0.05 prototype is 1.2 % low). The
   hat weights for a recorder-like jet (22.5 m/s, 0.4 mm, 44.1 kHz)
   are 73 causal taps whose current-sample weight is the early
   `s ~ t` response (0.5 % from half a period) and whose sum is the
   held gain to 1e-9. SnowySummit's causality admission: the transfer
   sampled on [0, 10] b/U0 and inverted puts less energy and L1 before
   time zero than any of the probe's pure-delay controls (5, 10, 15,
   20 b) on the same grid: energy 1.1e-5 against 4.2e-4 and L1 2.1e-3
   against 0.15 untapered, 3.0e-9 against 2.6e-6 and 6.8e-5 against
   0.015 under the 0.25 Gaussian taper (prototype kernels hold the same
   level over windows 30-150 b/U0). Windows stay at
   or below ~150 b/U0: beyond ~200 the downstream packet's growth
   brings FFT round-off into the global Poisson solve and every
   station blows up at once; from 90 to 200 the transfer moves at most
   3 % and 4 degrees at w >= 0.04.

## Error model

Typed `AeroacError` (non-finite, invalid-parameter, cancelled,
not-converged with disclosed residual, and the preserved shared harmonic-fit
refusal). No silent degradation, no fabricated eigenvalues.

## Determinism class

Deterministic: fs_math::det transcendentals + Dd arithmetic only, no
platform libm, no RNG. `forced_jet` splits its rows across
`std::thread::scope` workers; every row of both transforms is computed
whole by one worker, so the result is bit-identical for every worker
count (tested at 1 and 5).

## Cancellation behavior

Synchronous, milliseconds-class. No `Cx` integration (workspace
`frankensim-ccmn`). Exception: `bickley_step_response` is
seconds-to-minutes class (window 90 b/U0 at dx 0.1: 100-175 s on a
loaded 8-worker trj; dx 0.2: 52 s) and is not interruptible; callers
poll their gate before and after it.

`PreparedProperHarmonicField::observe_window_cancellable` accepts an explicit
caller gate and polls before work, at each factor and frequency pair, and before
publication. Cancellation returns `AeroacError::Cancelled` with no partial
statistics; the prepared bank is unchanged and can be observed again.
Preparation is separately bounded and synchronous. This callback does not
provide `Cx` ownership or implicitly alter other synchronous APIs.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/bickley_lbm.rs` (1): the LBM-vs-Rayleigh growth-rate
convergence fixture (invariant 8).

`tests/edgetone_staging.rs` (4, ignored heavy): stage-I Strouhal vs
Brown + Vaik/Paal published values; fine-lattice slit-fix regression
(invariant 10); adiabatic ramp hysteresis protocol (invariant 13);
turbulent-regime probe (stability-boundary bracket + tonal-window +
typed-destabilization pins).

`tests/noisetable.rs` (3): sweep + export + synth round trip;
typed refusals; duration-independent dimensionless synthesis amplitude.
Inline `noisetable` G1/G3 tests independently check windowed Parseval with
Nyquist and two-tone RMS under changed duration, source scale and band mapping.

`tests/jetlab.rs` (5): edge-tone oscillation + diagnostics +
radiation with scope; typed refusals; bitwise determinism; ramp
structure + rung grid + bitwise determinism; ramp + peak-instrument
typed refusals.

`tests/jetcard.rs` (13): tonal interim card mints the recorded
stage-I point with X-Struct residual + narrow claim + validity
refusals; exact round trip + stable content hash; broadband/refusal
minting laws (no broadband rung -> refuse; out-of-regime rung ->
refuse; boundary needs >= 2 admitted all-tonal rungs); the
no-residual X-Est cap and the forged-authority refusal; parser
refusals on foreign schema, laundered scope, truncation, and
trailing bytes; momentum-thickness tanh-edge law + refusals; pinned
schema constant; the sweep-receipt reader (`SlotJet3dRung::from_jsonl`,
`parse_sweep_receipts`): bitwise round trip incl. a non-finite flux
imbalance, pinned-order/schema/trailing-byte/non-finite refusals,
header-first and foreign-schema file refusals, and every archived
`tests/receipts/slot-jet-3d-re-sweep*.jsonl` parsing; the
temporal-Nyquist law (`peak_at_nyquist_edge`, `NYQUIST_EDGE_BINS`):
a rung whose raw peak sits within 8 bins of the record's Nyquist edge
is the bounce-back parity artifact (executed at rate 1.60: bin 4095 of
4096) and no minter admits it; the parity filter
(`parity_filtered_force_series`, `classify_rung_parity_filtered`):
non-overlapping pair averaging annihilates the period-2 line exactly,
keeps a tone's bin index and Strouhal (record duration unchanged; the
pipeline Strouhal is halved), refuses odd/short records.

`tests/forced_jet.rs` (6, full suite minutes-class): the published spatial mode
downstream; step kernel against harmonic forcing; the dx 0.2/0.1/0.05
ladder and worker-count bits; box and band invariance; hat weights and
the causality admission against a pure delay (invariant 14); local interval
weights against uniform hats at three clocks, tail support and invalid clocks.

`tests/aeroac.rs` (11): Wronskian identity; derivative cross-checks;
fsci-special oracle; small-argument limits; Hankel far-field
amplitude; Curle spreading + directivity + scope guard; Curle typed
refusals; Bickley pin self-verification + falsifiers; shooting solver
vs pins + grid convergence; Bickley typed refusals; bitwise
determinism.

## No-claim boundaries (the bead's remaining scope — OPEN)

- Edge-tone staging BEYOND stage I: only the stage-I point
  (Re 144, h/delta = 10) is validated. Hysteresis and stage-I-band
  confinement under adiabatic ramping ARE now measured (invariant
  13), but no stage-II/III LOCK exists on this rig by Re 264 —
  stage-ladder claims above stage I remain refused.
- Grid-convergence of source spectra across refinement levels: not
  implemented (single-resolution v1).
- Noise tables catalog THIS rig's tonal limit cycle in 2D; they are
  not turbulent flute-noise spectra. QUANTIFIED (executed probe,
  delta = 12 lu): the plain-BGK stable window ends at cell Reynolds
  `Re/delta` in (36, 48), and the spectrum is still an essentially
  pure tone (flatness ~1e-15) at the last stable rung. The fs-lbm
  `CollisionModel2::CentralMoment` operator (bead 3zkcr) removes the
  stability wall — the same rig runs through Re 2304 (Re_cell 192),
  now selectable via `JetLabiumConfig::collision` — but every
  executed rung REMAINS an essentially pure tone (flatness ~1e-18;
  free-jet-scale St ~ 0.467 at Re >= 1152). Bead l011o records that
  landscape as [`crate::regime::PINNED_2D_CENTRAL_MOMENT_TONAL`] and
  [`crate::regime::two_d_broadband_refusal`]: cataloging those rows as
  broadband flute noise is a typed refusal. The D3Q19 central-moment
  operator is live ([`crate::regime::evaluate_slot_jet_3d_operator`])
  but no 3D Re-sweep has been executed; a demonstrated 3D broadband
  regime remains a no-claim. The 2D inverse-cascade mechanism is the
  recorded physics boundary.
- Jet cards: the tonal interim card claims EDGE-TONE CLASS ONLY —
  it is not flute broadband, carries no band noise content, no
  absolute SPL, and no oscillation-amplitude authority; consuming
  it as a noise source is a type error by construction. The
  broadband and refusal-boundary minters are law + schema today:
  no card of either kind has been minted from executed 3-D sweep
  receipts (the .10.1 sweep is the pending producer). The card
  content hash is FNV-1a: deterministic correlation identity only,
  never authenticity.
- Quadrupole (volume) sources: dipoles dominate at low Mach over
  rigid surfaces; the Lighthill volume term is out of scope v1.
- Only orders 0 and 1 of the Bessel functions are provided (all the
  2D monopole/dipole machinery needs); higher orders are not
  implemented.
- The Rayleigh oracle is INVISCID (no critical-layer viscous
  correction); exactly-neutral real-c modes are ill-posed for
  shooting and are validated by approach from the unstable side.
- Forced jet: the base flow is parallel (its own viscous spreading
  over the window is neglected; the disturbance is viscous); 2-D (no
  span ends); linear (saturation belongs to the consumer); the exit is
  a baffle wall (no channel profile development, no chamfer); the
  labium does not react back on the jet, which runs on freely past
  the station; the response is held after the window, so content
  below w ~ 2 pi / window is the held limit, not computed, and the
  slowly decaying tail the hold cuts costs 2.5-3.9 % at w 1-2 at
  window 90 (invariant 14); the
  spatial growth rate is low by the critical-layer resolution (4.4 %
  at dx 0.1, invariant 14).

## Complete-sphere compact-array radiation power

`curle3d::PreparedCompactRadiatedPower` integrates the cycle-mean outward power
of prescribed peak volume-flow sources Q [m^3/s] and forces ON the fluid F [N].
All sources share one frequency, coordinate frame and emission clock, using the
original `FluidAcousticState`. With Z = rho c, its authority is

`P = k^2/(8 pi Z) <|sum_j (Z Q_j + n.F_j) exp(-i k n.s_j)|^2>_sphere`.

The angular average uses analytic scalar, vector and tensor spherical moments,
retaining separated monopole/monopole, force/force and mixed Q/F interference.
There is no angular quadrature in preparation or evaluation. Each individual
source must admit the compact-source model; the complete array need not be
compact. The `from_observer` adapter retains the original material, frequency
and source geometry. It explicitly refuses molecular observers: attenuation
changes the outward flux with radius, so a lossless integral cannot represent
their emitted power.

Exactly colocated sources are field-summed before quadratic products, in their
original input order. Normalization, field sums and contractions use
double-double arithmetic. A globally field-summed colocated baseline is combined
with separate `j0(x)-1` and `j1(x)/x-1/3` corrections; small source-separation
power is not obtained by subtracting two rounded order-one powers. This improves
weak residuals but is not a certified relative enclosure at every null.
Nonfinite or negative numerical power refuses instead of being clamped.

For N original sources and M distinct positions, preparation is O(N^2), stored
data are O(N + M^2), and coherent evaluation is O(N M + M^2), allocation-free and
without trigonometry. Group sums are recomputed to keep application immutable
and avoid caller-owned scratch. The independent-source path takes O(N).
Explicit coherence factors obey E[xi_r conj(xi_s)] = delta_rs: sources and Q/F
mechanisms interfere within a factor, while factors add powers. This
cycle-average observable needs no proper-complex or Gaussian assumption.
Factors use peak line units; a one-sided PSD bin of width df needs covariance
2 df S, not an unscaled PSD factor.

Source count, N^2, factor count and N^2 times factor count are checked before
their admitted work or allocation. Preparation and evaluation poll cancellation
and publish no partial operator or power. This is a homogeneous, stationary,
lossless free-space result for prescribed sources; it does not add scattering,
mutual loading, moving sources, finite-source validation or measured-SPL
qualification.

Nine authored regression tests cover SI normalization and observer binding,
separated and mixed-source phase signs, full pressure-velocity flux through
enclosing spheres, weak colocated residuals, tiny monopole and parallel/transverse
force pairs, explicit coherence factors and unitary mixing, rigid transformations,
resource limits, cancellation and molecular-model refusal. Rust execution remains
unverified in this session because the native process runner is unavailable.
Independent JavaScript angular and full finite-radius flux checks agreed to
1.33e-15 relative error for the documented three-source fixture; those checks
validate the reference equations, not Rust compilation.

## Mixed deterministic, real and proper acoustic source ensembles

`observation::PreparedHarmonicEnsemble` combines an optional deterministic
mean with explicit real or proper-complex fluctuation factors on one original
receiver/frame/source clock. A real factor has E[xi]=0 and E[xi^2]=1, retaining
its authored waveform phase and both sum/difference-frequency products. A proper
factor has E[z]=E[z^2]=0 and E[|z|^2]=1, so its sum-frequency terms vanish.
All factors have the declared joint identity covariance, zero cross-factor
moments and zero mixed real/proper moments. No Gaussian law, unknown-phase
independence, stationarity or PSD-bin normalization is inferred.

Real factors can represent phase-locked source-amplitude variability and
nonzero complex pseudocovariance. A real covariance square root of the stacked
real/imaginary source coefficients describes general positive-semidefinite
nonproper fluctuations by construction. Proper factors can coexist without
doubling their rank. The supplied physical peak field coefficients already
carry the square-root weights and propagation phases.

The deterministic, proper-only and mixed owners now share one field-moment
implementation. Coalesced equal-frequency fields retain their DD real/imaginary
components through every later contraction; they are not rounded back to C64
before products. Hz-to-radian conversion and deterministic rotations use the
shared extended kernel, preserving high-carrier beat frequencies and small
cosine decrements. Means integrate the retained field against an exact unit
DC field through the same kernel. These measures improve cancellation but do
not establish relative accuracy at arbitrary nulls or certify transcendental
errors.

`WindowedAcousticEnsemble` reports ensemble/window means, uncentered pressure
mean square, the full real symmetric velocity second-moment tensor and signed
mean p*u [W/m^2]. Deterministic motion remains in the tensor; it is not centered
fluctuation covariance or a detector-confidence interval. The deterministic
owner also exposes full moments through `observe_window_moments` and its
cancellable counterpart.

For fixed spatial dimension, admission covers total factors, supplied lines
and the sum of per-bank squared line counts before allocation. Preparation
takes O(sum N_r^2), stores O(sum N_r) fields, and polls within grouping work.
Observation takes O(sum M_r^2) for the distinct frequency counts M_r, allocates
nothing and reuses each temporal kernel across the field components. Both
legacy owners now have cancellable preparation; deterministic queries and full
moment queries also have cancellation callbacks. Refusal publishes no partial
operator/result and does not mutate an accepted owner.

Seven authored G0/G1/G3/G4 regressions cover an exact eight-state real/proper
ensemble against direct real-time integration, short-window nonproper power
and signed flux, real factor-basis invariance and proper-as-real decomposition,
high-carrier Hz beats, weak coherent residuals, physical clock/frame covariance,
resource boundaries, cancellation/replay and darkness. The exact shared kernel,
updated acoustic owners, seven new tests and eight prior proper-observation
tests pass one focused native Rust type analysis (31 test definitions including
16 shared-kernel tests). This is not execution or a full-crate build; both
remain blocked by exhausted storage. An independent analytic versus eight-state
real-time reference agrees within 3.22e-15 maximum absolute error for the
nine reported fixture moments.

## Time-dependent pressure/velocity product spectra

The original `PreparedHarmonicField`, `PreparedProperHarmonicField` and
`PreparedHarmonicEnsemble` owners expose `product_spectrum` and cancellable
preparation. `AcousticFieldComponent` selects pressure or a velocity axis
in the original Cartesian frame. The returned `AcousticProductSpectrum`
exposes its components and shared analytic spectrum, evaluates it at a
physical source-clock time, and observes either normalized detector profile.

All adapters traverse their retained coalesced DD fields and route bilinear
terms through the one fs-fft product-spectrum engine. No parallel temporal
implementation or new propagation model is introduced. Deterministic means
and real factors retain sum and difference tones; proper factors retain
difference tones. The prepared mixed owner's declared joint second moments
govern which factor contributions add, without inferring independence from
unknown phase or assuming stationarity or Gaussian statistics.

Pressure/pressure is in Pa^2, velocity/velocity in m^2/s^2 and either ordering
of pressure/velocity in W/m^2. The last quantity is signed physical p*u,
including reversible storage exchange and broadband beats; it is not a
local heating rate or a pressure-only plane-wave estimate. Ensemble results
are expected quadratic waveforms, not random realizations or fluctuation
noise spectra.

Before allocation, admission counts the sum of per-bank squared distinct
line counts and the real/proper generated-term counts with checked arithmetic.
The lazy contraction iterator avoids a second O(M^2) input array. The caller's
workspace limit and the shared fixed term ceiling apply. Invalid axes,
resource refusal and cancellation publish nothing and leave original owners
reusable. Waveform and window queries allocate nothing.

Four integration tests use the original public owners to check signed
standing-wave flux and broadband p/u products against direct real fields,
an exact eight-state mixed ensemble and its phase-locked subcycle distinction,
Uniform/Triangular means against all three original moment owners, and
retained weak coherent fields with budgets, cancellation and replay.
A focused native Rust harness compiled and executed the exact shared/acoustic
modules with the real fs_math library and only a nominal AeroacError boundary:
all 31 tests passed, including the five new shared-spectrum tests and four
original-owner acoustic integration tests. Independent standard-library
reference calculations agree with the eight-state time-domain products
within 6.67e-16 maximum absolute error; Uniform and Triangular weighted
real-time quadrature agree within 1.32e-14 and 2.69e-14. This scoped execution
and these references are not a full-crate build or repository quality gate.

## Coherent volume sources in uniform subsonic flow

`curle3d::array::convected::PreparedConvectedMonopoleObserver` extends the
prescribed volume-source lane to a constant laboratory-frame fluid velocity.
The original `FluidAcousticState` supplies density and adiabatic sound speed;
all sources and the observer remain fixed. The source injects volume at the
mean velocity with no extra force or entropy. Force-dipole wakes, nonuniform
or vortical mean flows, boundaries, absorption, source compactness validation
and reciprocal source loading remain outside this owner.

The original `moving::flow::UniformFlowCompactObserver` supplies every source
path's pressure coefficient and travel time through its Galilean moving-source
law. The array adds the exact convected Green potential's particle velocity
through its spatial gradient, using the SAME retained travel phase once.
The original downstream rationalization and explicit subsonic Mach margin
remain authoritative. Source phases are summed before quadratic observables, and
the original shared complex accumulation and source-factor statistics retain
coherent cancellation and complex partial correlation. A harmonic-line export
feeds the existing finite-window owner without adding another delay.

Pressure-work intensity `<p u>` is distinct from total quadratic perturbation
transport `J=<p u>+U<E>`, with
`E=p²/(2 rho c²)+rho |u|²/2`. The latter satisfies
`partial_t E + div J = 0` away from sources for the declared uniform-flow
linearized Euler system; it does not claim the change in total mean-flow energy.
Both coherent and ensemble outputs retain this distinction. Preparation is
O(N), retaining four complex coefficients per source; field applications are
allocation-free O(N), and authored factor applications O(N rank). Source-count,
factor-work, distance, finite-range and cancellation checks remain explicit.

Nine G0/G1/G3/G4 regressions cover original-owner pressure/delay parity,
analytic quiescent and axial/transverse
fields, finite-difference Euler/continuity and energy divergence, coherent nulls
and partial coherence, material scaling, coordinate transforms, flow-reversal
reciprocity, nearly sonic admitted travel times, work refusal and cancellation.
They pass a focused native exact-source leaf build using the production
`FluidAcousticState` implementation, production complex statistics and complete
`fs-math` crate, plus the original Galilean pressure and moving-source owners,
with nominal error carriers only. This is not a full-crate
Cargo/DSR gate or experimental source validation.

## Measured force coherence through the original line-dipole observer

`curle2d::PreparedDipoleObserver::evaluate_welch_bin` propagates a joint
recorded force cross spectrum through the original near-field pressure and
particle-velocity transfers. Channels are ordered
`[Fx_0,Fy_0,Fx_1,Fy_1,...]` in N/m, with one source frame and physical clock.
The prepared frequency must exactly match the selected estimator bin.
Each retained segment vector is summed coherently across force axes and
source paths before its second moments enter the ensemble average. Measured
cross-source and cross-axis phase, partial correlation and cancellation
therefore survive; unknown phase is never replaced with a synthesized one.

The shared one-sided /Hz factors are multiplied by `sqrt(2*df)` once, so
the existing peak-phasor factor of one half produces this bin's pressure
power, velocity covariance and signed active/reactive intensity. Results
are integrated bin quantities, not densities per hertz or window-corrected
tone amplitudes. In particular, a periodic-Hann resolved tone contributes
two thirds of its energy to its central bin. DC and Nyquist are refused.
Adding bins assumes a stationary resolved spectrum; finite-time
cross-frequency beats and statistical confidence require their own evidence.

`MaterialJetLabiumRun::prepare_force_cross_spectrum` observes the accepted
run's actual Cartesian force history without copying or transposing it and
retains its original timestep and `FluidAcousticState`. The resulting
`JetForceCrossSpectrum::observe_bin` reuses that same medium and the original
single-source observer. It cannot silently substitute nominal air, scalar
force magnitude, a new clock or inferred absolute-SPL authority.

Propagation is allocation-free O(sources * retained segments) after
preparation. Explicit source-factor work limits and cancellation guard the
observation. The source compactness, two-dimensional free-space model,
homogeneous lossless medium and original absolute-level boundaries remain.

Four focused native tests exercise the actual FFT/Curle/material snapshot
algorithms: independently authored complex forces against all field moments;
rectangular and Hann bin normalization; equal, opposite and partially
correlated source histories; a pressure antinode with zero energy flux; and
shape/frequency/work refusal, cancellation and deterministic replay. The
manufactured accepted-run adapter regression additionally checks original
medium/clock provenance and the public path, but requires the missing
constellation dependencies to run. The native probe uses only a nominal
outer error carrier and is not a full-crate integration gate.
