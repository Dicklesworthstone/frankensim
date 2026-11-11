## Joint recorded volume/force sources from resolved flow

`curle3d::array::recorded` binds original source positions, ordering, physical
clock and acoustic reference to one streaming Welch factorization. Each source
retains real `[Q,Fx,Fy,Fz]` in m³/s and N. Cross-source, cross-axis and
volume/force coherence survive into the existing complete pressure/vector-
velocity operator before bin-integrated moments and signed intensity are formed.
Lossless, original molecular and existing ideal-image propagation remain
available without a second window, conjugation, attenuation or normalization.

`PreparedPermeableSurface3::begin_recording` includes its current accepted
mass/momentum sample once; `record_sample` appends subsequent original accepted
samples. A changed reference, source order, geometry or off-clock time refuses
before acceptance. There is no implicit adaptive-time resampling, missing
prehistory, per-face power addition, source derivative or additional CFD step.
The selected stationary FW-H source scope remains: exterior volume sources
may still be needed, and enclosed sources must not be counted a second time.

Geometry is included in the explicit transform memory cap. Accepted sampling
allocates no adapter buffers; the original segment buffer and spectral factors
replace retention of the complete source history or dense source covariance.
Focused regressions cover phased source interference, a pressure null with
nonzero velocity, exact joint batch/stream agreement, cancellation rollback and
the existing accepted-CFD versus independent-DFT field comparison. The new
source and test module pass pinned rustfmt; native tests remain unexecuted after
a resource-limited real-source build stopped in asupersync before fs-couple.

## Permeable control-surface sound from accepted compressible flow

`curle3d::quadrupole::flux::surface::PreparedPermeableSurface3` extracts a
closed stationary Cartesian FW-H surface from the original frozen-gas or
evolving-mixture flow owner. Each original face supplies volume source
`Q=rho*(u.n)*dA/rho0` and force source
`F=[(p-p0)*n+rho*u*(u.n)-tau.n]*dA`. Euler constructors omit tau; molecular
constructors retain the actual configured molecular/bulk and selected SGS
traction. Conserved momentum and the original EOS provide the physical
advective and pressure terms, with each adjacent cell's flux formed before
centered face interpolation. The viscous term comes directly from the original
conservative face, without cell-stress averaging or another gradient estimate.
The FW-H convention follows Lockard and Casper (2005), equations (1)-(6):
<https://ntrs.nasa.gov/api/citations/20050182105/downloads/20050182105.pdf>.

The complete face bank retains physical midpoint positions, areas, outward
orientation, source order and proper stationary frame rotation. Native
`tau[:,axis]` is read in the positive grid-axis direction, rotated once, and
subtracted with the outward sign and full physical face area. Subtraction
remains in double-double arithmetic until the final force component, so a
small residual is not erased by first rounding a large inviscid loading.
Exact original conserved momentum supplies the unchanged mass source.
The absent-stress branch preserves the original Euler arithmetic.

Actual accepted-time `[Q,Fx,Fy,Fz]` histories feed the existing coherent
compact-field and Welch owners. All faces and channels retain their phases
before pressure, vector velocity or intensity is formed. Preparation checks
the explicit ambient rho0/p0/c0 against the original EOS; a mixture requires
its original species-ordered ambient partial-density vector. No reference
composition is guessed, and no separate propagation speed is fitted.

Surface reconstruction uses bounded O(faces) work and reusable adapter
storage. Molecular entry points additionally evaluate the original full-grid
stress once, reusing its allocated native caches. All channels and the native
clock publish atomically, so failed/cancelled sampling preserves the preceding
accepted observation. The retained `includes_viscous_stress` mode refuses
sampling through the other source lane. No accepted flow state, closure
history, interval report or boundary impulse changes. The selected box
requires an actual cell on both sides of every face; ghosts and periodic
wrapping are not invented. Source, geometry, workspace, medium and strictly
increasing clock admission remain explicit.

Scope is stationary geometry and a quiescent homogeneous exterior reference.
Face reconstruction and surface quadrature require spatial refinement. The
FW-H exterior volume term remains necessary wherever exterior stresses are
significant; enclosed volume stresses must not be added twice. No moving
surface, missing SGS closure, wake correction or measured source qualification
is inferred. Existing tests specify physical density/momentum/rotation,
independent monopole recovery and exterior-wave cancellation, original mixture
calorics, transactional refusals and coherent history observation against a
direct DFT. Added molecular cases cover all-six-face viscous signs/rotation/
areas, constant-stress closed-surface cancellation, exact zero-law recovery,
a 1 N residual after cancellation of 1e16 N loads, and original species with
later accepted molecular sampling.

## Full moving volume, force and stress arrays

`PreparedMovingSourceArray::multipole_field_at` composes the original moving
compact and stress-source owners at one physical microphone event. Every
source/image path solves its causal root once; that same emission time and
original position, velocity, acceleration and jerk prepare both mechanisms.
The caller supplies their original physical strengths and the distinct,
spatially weighted force/stress past-worldline integrals. Missing prior forcing
or an initial acoustic field is not reconstructed from a scalar impulse.

All source mechanisms and direct/ideal-plane/wedge paths remain in compensated
pressure and Cartesian particle-velocity sums before rounding or signed `p*u`.
The original image parity and vector transformation apply once to the complete
field. Microphone motion changes the causal observation event; it adds no
acoustic velocity or pressure multiplier. Existing uniform-flow preparations
retain their fluid-rest-frame observation and exclude mean-flow energy advection.

Queries allocate no additional storage and preserve original aggregate source,
path and retarded-iteration limits. History callbacks retain their own
clearance, integration accuracy and work limits. Native compactness, lossless
propagation and ideal-boundary assumptions remain; this is not finite-body
scattering or feedback into the producing structure/flow. Five regressions in
`tests/moving_source_array/multipole.rs` check independent stationary Green
derivatives, accelerating/wedge composition, coherent cancellation with weak
stress or surviving velocity, original Q/F recovery and complete-history refusal.

## Broadband coherent mixed sources at a dispersive impedance wall

`PreparedImpedanceWallCompactSourceObserver::evaluate_welch_bin` propagates
original synchronized `[Q,Fx,Fy,Fz]` factors. Its
`evaluate_welch_bin_with_quadrupoles` adds physical
`[xx,yy,zz,xy,xz,yz]` stress blocks from the SAME record. Both banks require
identical fluid, oriented plane, impedance, receiver, exact frequency and
reflected/total field. Every source, tensor component and direct/returned path
enters one compensated pressure/vector-velocity sum before second moments.
The native peak-bin owner applies `sqrt(2 df)` once. Geometry and channel order
of raw mixed records remain caller-bound; no correlation or duplicate
permeable-surface volume source is inferred.

`wall::recorded::observe_welch_band` evaluates the caller's original wall law
once at each selected interior-bin frequency. Only impedance may vary across
the band. It retains one frequency's native transfer banks at a time, shares
a total quadrature cap across all sources/frequencies, and checks aggregate
bin/segment/source work before propagation. One compensated accumulator forms
pressure power, full velocity covariance and signed active/reactive intensity.
Results add resolved stationary-bin moments; they do not supply deterministic
cross-frequency phase, finite-window beats or a transient waveform. Source and
wall variation inside a bin, and spectral leakage, require frequency refinement.
DC, even-record Nyquist, changed domains, incomplete channel layouts, exhausted
work and cancellation refuse without publishing a partial band.

The existing original-recorded-stress observers and their geometry/clock
binding remain available. `tests/impedance_wall_mixed_recorded.rs` adds
independent spherical image derivatives, mixed pressure cancellation with
surviving velocity, finite wall dispersion against an explicit windowed DFT,
original recorded geometry/clock, and combined admission/refusal checks.
The original `tests/impedance_wall_recorded.rs` stress-only suite is retained.

## Moving-receiver history composition seams

`PreparedKinematicTrajectory::next_source_knot_after` exposes only the next
original accepted geometry knot so the existing history renderer can integrate
each original polynomial piece without rebuilding motion.
`PreparedReflectedMovingReceiverSource::force_history_path_at` exposes the
actual virtual receiver at one arrival event, pressure parity and the same
parity-weighted transpose image Jacobian used by native full-field reflection.
A past-force integral must use this current spatial point throughout its source
history; an integral from a different image or receiver event is invalid.

`PreparedMovingForceFieldSample::evaluate_dd` evaluates the existing local
Q/F field with a DD past-force integral and retains DD pressure and velocity
for coherent downstream source/image accumulation. No new propagation law,
root, strength interpolation, homogeneous initial field or per-source intensity
is introduced. The accepted-history owner still supplies causal coverage,
past-worldline clearance and bounded quadrature. The original f64 field APIs
continue through the same native evaluation.

# fs-aeroac — Contract

## Spherical reflection at a local impedance plane

`curle3d::refraction::full_wave::wall` propagates prescribed compact volume and
force sources through the outgoing angular spectrum of an infinite locally reacting
wall. Its outward normal points from fluid into wall, `p=Z*u_n`, with specific
impedance in Pa s/m and exp(-i omega t). The original lossless fluid supplies
density and sound speed. Frequency-dependent wall laws must be evaluated by
their original material/mechanical owner at the requested frequency.

The reflection is `(Z*zeta-rho*omega)/(Z*zeta+rho*omega)` for every propagating
and evanescent component. The same integral supplies pressure and all particle
velocity components, including induction and coherent direct/reflected
interference. Normal-incidence matching therefore does not erase a spherical
near-field reflection. Exact rigid and pressure-release walls use analytic
images without numerical quadrature.

Finite passive impedances retain evanescent amplification, including a lossy
compliant surface mode. Its analytic maximum enters the omitted-tail bound;
splitting at that maximum prevents a narrow mode from escaping panel resolution.
Fractional coefficient-variation controls resolve high gain without requiring
work proportional to the peak amplitude; the requested field-error tolerances
still govern quadrature acceptance.
Pure compliance with zero resistance and positive reactance retains its outgoing
surface wave by exact limiting absorption. For normalized `z=iX`, `a=k/X`, its
evanescent coefficient is `1+2a/(t-a-i0)`. Each complete pressure or velocity
derivative channel is integrated with analytic principal-value subtraction and
the positive `i*pi` residue, retaining the same source and receiver derivatives.
An analytic derivative supplies the removable-singularity value when a sampled
coordinate rounds to the pole. The tail begins at `T>=2a`, where `|R|<=3`.
The four cached pole/derivative samples count against the original work budget;
absolute tolerance controls the regularized integral because the analytic pole
term can cancel it. Small reactance can require a large cutoff and bounded
preparation may refuse rather than claim uniform efficiency in that limit.
No resistance floor or reflection clipping is introduced. Nonlocal/porous tangential motion,
finite edges, multiple walls, mean flow, bulk attenuation and source feedback
are outside this model.

Preparation shares the existing evaluation/depth budgets, compensated adaptive
Gauss--Kronrod integrator, cancellation polls and Bessel accuracy limitations.
Quadrature estimates and the separate floating-point analytic tail bound cover
the reflected transfer, not geometry or roundoff; no interval certificate is
claimed. Application reuses four cached complex coefficients without allocation
or transcendental work. A coincident source/receiver may query the finite
reflected field but cannot query the singular total field. G1/G3 tests in
`impedance_wall_reflection` check independent analytic ideal images, the Robin
boundary law and wall absorption, lossy surface-mode response, reciprocity,
independent exponential-integral references for narrow passive and exactly lossless
surface modes within fixed work budgets, positive-resistance continuity on and off
axis, three quadrature tolerance rungs, cancellation and work/physical refusals.

`PreparedImpedanceWallForce` differentiates the original source position in this
same Green field for a physical complex force ON fluid [N]. Seven analytic
channels retain pressure, the full velocity tensor, opposite mixed-derivative
signs and transverse motion on the symmetry axis. Exact rigid/release image
forces reverse their normal component. Finite walls share the monopole's
coefficient, surface-mode split, lossless pole prescription, resolution and tail owner, with the extra
derivative powers included in force-specific error bounds. Preparation caches
reflected and coherent total 4-by-3 Cartesian matrices; application performs
compensated linear reductions without allocation, integration or numerical
differentiation. `impedance_wall_force` tests rotated analytic images, all
pressure/velocity source derivatives including pure compliance, Robin pressure/velocity and passive wall
absorption, independent narrow-surface-mode references within a fixed work
budget, finite reflected self-fields, original-clock harmonic lines and refusals.

`PreparedImpedanceWallQuadrupole` carries a full original symmetric Lighthill
volume-stress moment [N m] through the same passive-wall Green field. Its six
physical components `(xx,yy,zz,xy,xz,yz)` retain trace and both shear contractions.
Four pressure Hessian and ten velocity third-derivative channels preserve
transverse motion at pressure nulls. Exact rigid/release limits reflect both
tensor indices with `S M S`, `S=I-2nn^T`; finite walls reuse the original spectral
coefficient, lossy surface-mode split, lossless pole prescription and evanescent amplification bounds with
second/third derivative powers and Frobenius shear weights. Cached reflected
and coherent total 4-by-6 transfers apply without allocation or integration;
harmonic lines retain the original source clock. A receiver may lie on the wall,
or coincide with the source for reflected-only evaluation. Native diagnostics
separate estimated quadrature error from analytic tail bounds per N m.
`impedance_wall_quadrupole` checks independently differentiated spherical image
fields, every source derivative of the wall force field, the rotated Robin law
and passive energy absorption, and independent exponential-integral derivatives
for a narrow lossy surface mode including axial pressure-null shear velocity.
This remains a prescribed compact moment, without finite-source resolution,
unresolved subgrid sound, source loading, finite edges or bulk propagation loss.

`PreparedImpedanceWallQuadrupoleObserver` extends the same native wall kernel
to original volume-stress cell assemblies, with reflected-only and coherent
total views. One cached physical 4-by-6 map per cell feeds the existing generic
stress-array reducer: tensor entries, trace and both shear contractions remain
coherent, and explicitly independent cells still retain their own direct/returned
cross terms. Caller-supplied low-rank factors retain inter-cell phase before
forming pressure power, full velocity covariance and signed active/reactive
intensity. Original-clock harmonic lines add no delay. A shared preparation
budget bounds all cell quadrature, retaining native per-cell diagnostics;
warm application allocates nothing and performs no spectral integration.
Reflected self-fields and wall receivers are admitted; singular total self-fields,
invalid source counts, exhausted budgets and cancellation refuse atomically.
`impedance_wall_quadrupole_array` exercises native-cell superposition, finite-wall
Robin absorption, destructive and partial coherence, ideal-wall path cross terms,
aggregate work, cancellation and self-field admission.

`RecordedQuadrupoleSpectrum` now binds its original stress cells and exact
interior-bin frequency to reflected-only or coherent total impedance-wall
observation. The prepared wall array delegates each original Welch factor to
the unchanged generic stress reducer, retaining full pressure/velocity moments
and cell, tensor and direct/reflected cross terms. The native peak-bin owner
applies sqrt(2 df) once; no extra density, cell volume, phase origin, propagation
delay or FFT is introduced. Recording-aware evaluation rejects changed positions
or order and enforces the original frequency and combined cell/segment work.
The caller supplies the original material impedance at each requested frequency;
a prepared response cannot be applied at another bin. Results are bin-integrated
moments with the original window, aliasing and finite-record limitations, without
waveform continuation, cross-frequency phase or inferred source independence.
`impedance_wall_recorded` specifies known-tone normalization and full moments,
ideal-image interference, finite-wall Robin absorption and physical binding,
work and cancellation refusals.

`PreparedImpedanceWallCompactSourceObserver` prepares complete physical Q/F
source assemblies for reflected-only or coherent total observation, including
receivers on the wall. One 4-by-4 Cartesian transfer per original source feeds
the existing compensated array/statistics owner. Source and mechanism phases
survive until pressure and vector velocity are formed; explicitly independent
sources still retain their own Q/F and direct/reflected cross terms. Caller
supplied low-rank coherence factors propagate the full pressure/velocity second
moments without inventing independence from unknown phase. Harmonic lines keep
the original source clock for finite-window observation. Preparation shares one
aggregate Q/F evaluation cap across all sources and retains separate native
unit-amplitude diagnostics; application performs no allocation or quadrature.
`impedance_wall_array` compares rotated ideal-wall assemblies against analytic
Green-function images and dense source covariance algebra, checks finite-wall
Robin pressure/velocity and absorption for coherent and partially coherent
sources, and exercises aggregate work, cancellation and direct-self-field limits.

## Exact physical records and sampled harmonic receivers

`observation::AngularHarmonicFieldLine` and `PreparedHarmonicField::new_angular`
retain source-authored rad/s without a round trip through f64 Hz. The original
angular values govern coalescing, samples, finite-window p/u moments, product
spectra and output-Nyquist admission. Existing `HarmonicFieldLine` constructors
retain their original mathematical Hz conversion and observation semantics.
Both paths share amplitude accumulation, work budgets and cancellation. The
`fs-couple` port-harmonic bridge uses this angular path; its long-epoch regression
checks original-clock pressure, intensity, quadratic spectra and calibrated PCM.

Material-bound 2D and periodic 3D jets retain every admitted recording step,
including odd and prime record lengths. Their prepared force spectra and the
shared transverse periodogram use fs-fft's cached GeneralFft. The physical grid
remains k/(N*dt); only the internal Bluestein convolution is padded. Even-N
Nyquist is excluded and odd-N floor(N/2) is retained. Hann windowing, DC removal,
peak-force normalization, original material and N/m or periodic-cell N units
keep their existing owners. Sample-count admission bounds work before numeric
allocation; the derived FFT payload is checked explicitly. The 2D diagnostic
average counts every actual 16-step sample, including a partial final group.

JetForceTones::prepare_observer_field and JetForceHarmonics::prepare_observer_field
propagate each actual fitted frequency through the original near-field Curle
pressure/velocity operator. This produces one coherent PreparedHarmonicField<2>
for analytic windows or sampled pressure, with phase at original force sample
zero. DC force and unexplained residuals are excluded from acoustic tones.
The original frozen homogeneous lossless 2D source and shape/scaling limits
remain; this adds no isolated 3D radiation or absolute-SPL authority.

PreparedHarmonicField samples physical pascals through HarmonicSampleClock at
epoch + integer_index/sample_rate. DD clock arithmetic retains long counters
and sub-ULP sample intervals. Shared Hz phase conversion retains the low 2*pi
part, preventing rounded-constant drift at large source epochs; authored rad/s
semantics remain separate. Coalesced amplitudes and coherent reduction remain
in DD through the final scalar. Every line must be strictly below output
Nyquist, even at a cancellation null. Sample/line work is bounded before output,
and sampling allocates nothing. Absolute integer indexing fixes the arithmetic
independently of callback partitions and avoids accumulated oscillator drift.

These are observations of settled prescribed harmonics. No attack, source
evolution, missing phase/coherence, extrapolation of a nonstationary force
record, or recovery of already aliased measurements is inferred. Direct DFT,
physical source-to-pressure, independent real-waveform, analytic-window,
large-epoch, cancellation and output-admission tests cover these paths.



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
- `PreparedMovingCompactSample::evaluate_monopole_field` observes the full
  instantaneous field of physical volume injection on that SAME retarded
  event. The original potential `phi=-[Q/(4 pi R kappa)]ret` supplies pressure
  and all three particle-velocity components, including acceleration,
  transverse motion and `1/R²` induction. Q and Qdot keep the source clock;
  the microphone's observed Doppler rate is not an input multiplier. Signed
  intensity is `p u`, with no plane-wave substitution or phasor half-factor.
  Evaluation is allocation-free, bounded constant work with DD reductions.
  Nonfinite/unrepresentable velocity or flux refuses without narrowing the
  existing pressure API. Force-source velocity still requires additional
  history/initial conditions and is not supplied by this volume-only method.
  Existing source/receiver trajectory queries return the same prepared type.
  Ideal reflected paths transform its gradient by the transpose image map
  into physical receiver axes while pressure retains the original force frame.
  Coherent paths must sum p and u before forming total intensity. In a shared
  uniform-flow frame, u is the acoustic perturbation; `p u` excludes mean-flow
  advection. G1/G3 cases in `moving_compact`, `moving_trajectory`,
  `moving_receiver` and `moving_reflection` compare independently differentiated
  retarded potentials, analytic induction/Doppler limits, signed/material
  scaling, ideal-wall conditions and coordinate covariance. These tests do
  not qualify source dynamics, finite size, nonlinear sound or measured SPL.
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
  `UniformFlowReceiver` binds multiple source paths to one admitted fluid,
  laboratory velocity, microphone and clearance policy. Exact emission-endpoint
  preparation preserves cached relative geometry at large clock epochs. Its
  native Eulerian pressure view carries the actual emission/arrival timestamps
  and residual; the retained fluid-frame `1/kappa` is not the laboratory
  observed source-time rate, which is exactly one for these fixed endpoints.
  `UniformFlowFrame` also maps accepted laboratory events into one shared
  fluid-rest frame: position uses the explicit reference event, velocity is
  relative to the fluid, and acceleration/time/source strengths are unchanged.
  `PreparedKinematicTrajectory::in_uniform_flow` transforms each original knot
  before the existing quintic compilation and whole-span speed admission, with
  no intermediate knot array or added runtime allocation. Both paths must bind
  the identical frame and original fluid. The existing moving-receiver causal
  owner and `fs-couple` pressure renderer then retain their original root,
  clearance, bandwidth, finite-history and cancellation policies. Explicit
  laboratory kinematic observation uses the inverse map. The general pair's
  whole-window rate bound remains conservative; instantaneous clock metadata
  uses both fluid-relative velocities without multiplying pressure by it.
  `tests/moving_receiver/uniform_flow.rs` checks moving/accelerating source and
  receiver pressure against independently differentiated convected Green
  potentials, actual source-time rate, bitwise zero-flow parity, large reference
  events, relative-Mach admission and frame/domain/resource refusals.
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
  The source contractions preserve both scalar-product and summation residuals
  through the coherent reduction, then round the combined physical field.
  This also applies to the native 2D dipole, convected monopole, volume-stress
  and recorded-source paths; a nearly cancelling pair must reproduce its small
  residual source, not a rounded-product noise floor. Signed statistical cross
  moments preserve their product residuals across independent factors. These
  arithmetic guarantees do not improve the cached Green kernel's geometric or
  physical-model accuracy, and no allocation or transcendental is added to the
  O(source count) application.
  Harmonic monopole/dipole and volume-stress preparation now share an extended
  source-to-receiver path: authored binary64 coordinates are subtracted, scaled
  and normed in DD, and `Re(k)*r` stays extended through trigonometric reduction.
  Short path differences therefore survive long carrier paths before the
  original C64 field coefficients are cached. This applies to the same native
  lossless/molecular kernels and their ideal images; source geometry and the
  material's original wave number are not altered or recalibrated. Distance
  metadata and cached coefficients remain binary64, and arbitrary near-null
  relative accuracy is not certified. `tests/compact_interference/path_phase.rs`
  checks an ultrasonic spatial null, a sub-ULP separation on a long path, and
  coherent mixed volume/stress cancellation against factored analytic Green
  derivatives. Repeated field/coherence evaluation incurs no added work.
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
- `PreparedCompactSourceObserver::with_wedge` adds the exact finite dihedral
  group for two infinite ideal walls opening by pi/m, with integer m in 2..=128.
  Equal rigid or pressure-release walls admit every supported order; mixed
  conditions require even m because the boundary signs must close after m
  double reflections. Geometry is admitted within 64 machine epsilons in the
  normal dot product and relative cross-product magnitude, without snapping.
  The complete 2m images carry both affine source positions and their full
  Cartesian force transforms, with source and receiver inside the physical
  wedge. Both wall conditions include reactive near-field velocity and the
  common edge. Preparation bounds all 2mN paths, uses O(N) scratch and folds
  sixteen coefficients per original source with compensated sums; runtime
  remains allocation-free O(N), including independent-source statistics that
  retain every source's own image interference. Original material relaxation,
  direct-path delays, source clock and statistical source count survive.
  General angles, odd mixed wedges, finite impedance and room reverberation
  require other Green functions and are not approximated by truncated images.
  `tests/compact_interference/wedge.rs` checks an independent polar image oracle,
  near-field momentum, all four source components on both walls and the edge,
  pi/2 parity with the corner owner, rigid frame and wall-order covariance,
  molecular dispersion, source coherence, cancellation and admission budgets.
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
  The rung verdict (peak, prominence, flatness) reads the admitted jet
  band only: periodogram bins from `RUNG_DC_GUARD_BINS` (4) up to
  `RUNG_BAND_TOP_STROUHAL` (St_h 1). A record whose band holds fewer
  than `RUNG_MIN_BAND_BINS` (128) bins refuses by name, so a pure tone
  cannot read broadband for want of bins. Verdicts written before this
  band law used an n/8 guard and flatness to Nyquist; those rows are
  not re-read.
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
- `forced_jet::{BaseSpreading, bickley_step_response_with_base,
  ForcedJet::with_base, JetStepResponse::station_width_ratio}` — MB
  97/125/216, bead pexfp.135.1: an opt-in base development,
  quasi-parallel. At each `x` the base is `U = a(x) sech^2(y / beta(x))`,
  with `beta = b(x) / b0` and the centre speed `a = beta^(-1/2)`, so the
  momentum flux is conserved. The deviation keeps the linearised equations
  with that local `U` and `U''`. The base's own transverse velocity and every
  streamwise derivative of the base are dropped (the quasi-parallel
  truncation). The centreline deviation moves at the local centre speed. The
  exit is the `beta = 1` profile, so a measured exit's receptivity is
  unchanged. Past the station the base is held at the station's.
  `LaminarBickley`: `beta = (1 + 6 x / Re)^(2/3)`, the centreline momentum
  balance with `u_c^2 b` conserved, virtual origin `Re / 6` half-widths
  upstream. `Linear { width_slope }`: `beta = 1 + width_slope x`.
  `Parallel` is the default, and its arithmetic is
  `bickley_step_response_with_exit`'s, bit for bit. Refuses a laminar
  Bickley base at a non-finite Reynolds number, a width slope that is not
  positive and finite, a hold station that is not positive and finite, and a
  band narrower than 4 of the held half-widths. NOT CLAIMED: the base's
  transverse velocity and streamwise derivatives, and the jet past the
  station.
- `forced_jet::{bickley_step_span, JetStepSpan}` — MB 23643, a jet length
  that moves during a note: ONE integration over stations `[near, far]`,
  with the base held past `far`, records every grid node the cubic reads
  from one node before `near` to two past `far`.
  `JetStepSpan::response_at(X)` rebuilds the step response at any `X` in
  the span by `displacement_at`'s own weights and the same trapezoid
  integrals. A span of one station is the single-station response bit for
  bit. A parallel base never reads the station, so every station of a
  parallel span is its own integration bit for bit. A spread base is the
  single-station response bit for bit at `far`; nearer stations read a
  base that spreads out to `far` (stated, small). Refuses stations that are
  not finite, positive and ordered, and a read outside the span.

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

`material::MaterialSlotJet3d` applies the same material/scaling owner to the
existing D3Q19 central-moment slot-jet rig. Its second-order relaxation rate is
derived from the gas's actual kinematic viscosity; only higher-order numerical
relaxation remains an independent control. Three-dimensional momentum exchange
converts with `rho * dx^4 / dt^2` into N across one explicitly retained periodic
span. Both material adapters share one Hann/coherent-gain/phase FFT extraction.
The physical line frequency and acoustic wavenumber use the same timestep and
gas as the flow. A three-dimensional periodic-cell force is not passed to the
two-dimensional Curle Green function or relabeled as microphone pressure.

The 3-D initial-speed admission includes the rig's 1.5 maximum spanwise seed
factor. Final-plane speed diagnostics use all three velocity components
(`momentum / density`), and plane mass flux is momentum itself. These are
final-plane samples, not a global speed bound or maximum over the recording.
`tests/material_slot_jet_3d.rs` executes pressure, temperature and composition
substitutions through the original solver, independently checks SI force and
DFT amplitude/phase, and checks necessary material/Mach/allocation refusals.
This frozen homogeneous low-Mach adapter supplies no compressible feedback,
chemical evolution, acoustic radiation, turbulence-convergence or SPL claim.

1. EXACT Wronskian: `J1 Y0 - J0 Y1 = 2/(pi x)` holds within 5e-13
   over a log grid 1e-3..3e3 including the series/asymptotic
   crossover from both sides. At large arguments all four J/Y owners
   retain their pi/4 or 3pi/4 shift by angle addition after one shared
   reduction of x; the outgoing H1 no longer combines differently
   rounded J1/Y1 phases. The amplitude uses `sqrt((2/pi)/x)`, retaining
   finite nonzero waves even where `pi*x` would overflow. Independent
   256-bit MPFR 4.2.1 references check actual J0/Y0/J1/Y1 phases at
   exact f64 inputs 1e4, 1e6, 1e12, 1e16, 1e100 and f64::MAX within
   8e-15 of the asymptotic amplitude envelope. Separate normalized
   Wronskian/cylindrical-energy-flux and modulus regressions cover the
   asymptotic crossover through f64::MAX. The Wronskian alone is blind
   to a common phase error. These are sampled accuracy checks, not
   uniform relative-error bounds at zeros or bounds on input kr error;
   downstream quadrature budgets and applicability limits are unchanged.
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
verdict reads the halved record at sample period 2, so its admitted
band and Strouhal are the raw record's), refuses odd/short records.

`tests/forced_jet.rs` (6, full suite minutes-class): the published spatial mode
downstream; step kernel against harmonic forcing; the dx 0.2/0.1/0.05
ladder and worker-count bits; box and band invariance; hat weights and
the causality admission against a pure delay (invariant 14); local interval
weights against uniform hats at three clocks, tail support and invalid clocks.

`tests/forced_jet_spreading.rs` (7): the parallel default unchanged; a
spread base that never widens returns the parallel jet bit for bit; the laws'
widths are their closed forms, and Bickley's satisfies its centreline momentum
balance; a laminar Bickley base at Re 1e8 tends to the parallel jet
(|dG|/|G| 2.990e-8); the Fletcher slope at the station moves |G(0.3)| by
0.1376; the measured exit stays in phase at the exit; laws, holds and
bands outside their domain are refused.

`tests/forced_jet_span.rs` (4): a span of one station is the single-station
response bit for bit (values, integrals, metadata) under the parallel and the
laminar Bickley bases and both exits (the linear base is not run); every station of a parallel span is its own integration bit for bit;
a spread span is exact at its far end and within 1e-2 of the single station
at its near end (printed); stations outside their domain are refused.

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
- Supplied symmetric 3D Lighthill volume-stress moments propagate through
  `curle3d::quadrupole`; `quadrupole::lbm` and optional `quadrupole::flux`
  extract selected native 3D flow moments into the joint recorder.
  Automatic source-region completeness,
  unresolved SGS completion and inclusion in the existing surface-only rigs
  remain out of scope; this capability alone does not qualify a jet-noise source.
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

## Directional interference and outgoing wave momentum

`PreparedCompactRadiatedPower::with_momentum` retains that original free-space
owner and prepares `PreparedCompactRadiatedMomentum`. For its normalized
amplitude A(n), it evaluates both P = <|A|^2> and the vector first moment
M = <n |A|^2>. The outgoing wave momentum flux is M/c in newtons, using the
same original material sound speed. M/P, when P is nonzero, is the
power-weighted mean direction; the exact model obeys |M| <= P.

A colocated volume/force source has M = (2/3) Re(a f*) in the original
sqrt(watt) normalization. Its relative Q/F phase can reverse directional
emission while leaving total sphere power unchanged. Separated-source
monopole, mixed Q/F and force/force interference uses the original lower
angular moments plus the rank-three moment, from the Legendre integral
[DLMF 10.54.2](https://dlmf.nist.gov/10.54.E2) and spherical-Bessel series
[DLMF 10.53.1](https://dlmf.nist.gov/10.53.E1). A small-argument series retains
j3 without upward-recurrence cancellation; j2/x is retained directly at large
separation. The original coalesced DD baseline and separate small remainders
remain intact, including directional residuals below order-one f64 power ulps.
This does not certify relative accuracy at every interference null.

Preparation adds two DD coefficients per distinct-position pair under the
original source/pair limits. Repeated coherent application forms power and
M in one traversal, allocating nothing and evaluating no transcendental
functions. Explicit independent sources and low-rank coherence factors retain
each source's Q/F phase; independent factors add both power and momentum.
Factor units and peak-phasor averaging remain those of the original power
owner. Cancellation publishes no partial result.

This is an outgoing-wave observation of prescribed Q/F sources in a
homogeneous stationary lossless free space. It is not an applied mechanical
recoil force or a coupled source/mean-flow momentum balance. Ideal-wall
owners explicitly refuse: dividing full-image scalar power by image-group
order does not give the vector momentum in the physical wall domain.
Molecular-model refusal remains inherited from the original power adapter.

Seven G1/G3/G4 regressions in `tests/compact_momentum.rs` cover analytic
cardioid phase steering, independent angular quadrature of separated mixed
and force-only sources, opposite monopole phase steering, weak cancellation,
the power/momentum bound, coordinate and material covariance, independent
and coherently mixed source factors, budgets, ideal-wall refusal and
cancellation. Native execution is required before claiming these Rust
regressions have passed.

## Directional radiation of complete physical multipoles

`power::PreparedMultipoleRadiatedPower::with_momentum` prepares
`PreparedMultipoleRadiatedMomentum` from the original Q/F/Lighthill-stress
owner. It returns the same power P, first angular moment M [W] and outgoing
wave momentum flux M/c [N] as the original compact momentum result, now
retaining all Q/F/T interference. The physical far-pressure source amplitude
is `rho c Q + n.F - i k n.T.n`; T is the original symmetric integrated
Lighthill stress [N m], with physical shear entries and its complete trace.
The trace is absorbed coherently into the scalar amplitude before quadratic
products, leaving a traceless tensor for the analytic contractions.

In normalized sqrt(watt) amplitudes (a, f, T), a colocated source has
`M = (2/3) Re(a f*) + (4/15) Re(T f*)`. A force/stress pair can reverse
directional emission through its physical relative phase without changing
sphere power. Separated scalar/stress, force/stress and stress/stress terms
use angular moments through rank five. The new rank-five identity is
`<n^5 exp(-i x n.e)> = -i [(j3/x^2) sym(delta delta e)
- (j4/x) sym(delta e e e) + j5 e^5]`, with 15 and 10 distinct terms
in the symmetrized sums. It follows from the same DLMF Legendre integral
and spherical-Bessel series cited above.

The compact and multipole owners share the rank-three coefficient evaluation
and momentum-result validation. Small arguments use fixed series; quotient
identities retain j3/x^2 even when j3 itself underflows. At large arguments,
j2/x, j3/x and j4/x are kept directly to avoid cancelling adjacent Bessel
orders. A DD globally coalesced baseline plus separate separation corrections
retains weak source and trace cancellation. This is not a certified relative
enclosure at arbitrary nulls; nonfinite DD components or negative numerical
power refuse instead of being silently repaired.

Preparation retains the original geometry, material, source count and pair
budget, adding three DD coefficients per distinct-position pair. Coherent
application computes P and M in one O(N M + M^2) traversal, allocation-free
and without transcendental calls. Independent sources take O(N) while keeping
each source's intrinsic Q/F/T coherence. Explicit source factors add both
power and signed momentum, with the original factor count and pair-product
budgets. All paths poll cancellation and publish no partial result.

Only the original homogeneous stationary lossless free-space sphere is
admitted. Ideal-wall owners explicitly refuse the vector extension: scalar
image-group division does not supply physical-domain momentum or wall recoil.
This observation does not add applied mechanical force, source feedback,
mean-flow momentum, unresolved turbulence, finite-source validity or
measured-SPL authority. MR35 and MR55 remain broader open physical programs.

Seven G1/G3/G4 regression groups in `tests/multipole_momentum.rs` cover
physical phase and SI normalization, independent direct angular integration
of all source mixtures, tiny-separation and weak tensor cancellation,
trace cancellation and the original compact limit, frame/material covariance,
explicit source coherence, power/momentum bounds, budgets, cancellation and
wall-domain refusal. Native execution is still required before claiming
these Rust tests passed.

## Finite-window real multipole radiation

`power::PreparedMultipoleRadiatedEnergy` binds original
`PreparedMultipoleRadiatedPower` owners and physical Q/F/stress banks through
`MultipoleRadiationLine`. Every bank must retain the same original
`FluidAcousticState`, free-space domain, Cartesian frame and source clock.
A supplied bank establishes the medium even for darkness; an absent bank
has no implicit medium. Unknown phases, independence, stochastic covariance
and PSD/bin normalization are not inferred.

The explicit sphere center defines x = center + R n and retarded sphere time
tau = t - R/c. In the radiation-zone limit, the original normalized amplitude
A gives instantaneous outward power
`P(tau) = 2 <[Re sum_l (-i A_l(n)) exp(-i omega_l tau)]^2>`.
This is the limiting outward pressure/velocity flux over enclosing spheres,
not local near-field pressure squared or a finite-radius storage balance.
Only uniform `HarmonicWindow` intervals are admitted as physical integrals.
Results retain mean power, energy [J], the directional first energy moment [J],
and outgoing wave momentum [N s], using the original c.

Both difference- and sum-frequency products survive. Their spatial wave
vectors are respectively
`k_a (s_a-center) - k_b (s_b-center)` and
`k_a (s_a-center) + k_b (s_b-center)`.
The native angular moments through rank five integrate every Q/F/T cross term;
`fs_fft::harmonic_window` supplies the real temporal Gram product, including
short-window sum terms and frequency beats on the authored Hz clock.
Joint translation of source positions and sphere center preserves the result.
Moving sources while holding the observing center fixed changes directional
arrival times and can change a finite window's energy, even for a monopole.

Exactly equal frequencies and positions coalesce physical source terms before
normalization. A globally coalesced baseline contracts directly with the
temporal kernel; complex spatial separation corrections stay separate in DD.
The cache retains conjugate and unconjugated angular products without rounding
each source-pair result to binary64. These measures preserve weak interference,
but do not certify relative accuracy at every null. Unrepresentable phases,
coefficients or output, and negative numerical power, refuse without clamping.

For L supplied banks, P original source terms and F distinct frequencies,
preparation costs O(L^2 + P^2), retained data are O(L + F^2), and each window
takes O(F^2) work with no allocation or source/angle traversal. Admission
checks counts, L^2 + P^2 and a conservative peak workspace bound before
allocation, including preparation buffers, the frequency-pair cache and
rank-five scratch. Every path polls cancellation and publishes no partial
operator or observation.

The seven regression groups in `tests/multipole_window.rs` compare independent
space/time integration of real Q/F/T fields, short-window monopoles, physical
off-center arrival times, destructive beats and large source epochs, interval
additivity and the original cycle-power/momentum limit, weak coalesced sources,
joint translation, exact budgets, common-medium admission and cancellation.
Native Rust execution remains required before claiming those tests passed.
This is prescribed finite-spectrum observation, not source turn-on dynamics,
a source-work debit, mechanical impulse, mean-flow/attenuating propagation,
ideal-wall momentum or reciprocal structural loading. MR35/MR56 remain open.

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
or vortical mean flows, finite-impedance boundaries, absorption, source
compactness validation and reciprocal source loading remain outside this
owner. The ideal reflection adapters are specified below.

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

## Integrated acoustic power in uniform subsonic flow

`curle3d::PreparedConvectedMonopolePower` integrates the total cycle-average
quadratic perturbation-energy transport `<p u> + U <E>` of fixed prescribed
volume sources Q [m^3/s], retaining the original material, mean velocity,
frequency and common source clock. Its source law and convected Green function
are those of `PreparedConvectedMonopoleObserver`, including the original
[Wang (2020), equations 3.28–3.30](https://pure.tue.nl/ws/portalfiles/portal/149260007/20200511_Z.Wang.pdf)
reference. Volume enters at the mean velocity without extra force or entropy.
This observable is distinct from pressure work alone or mean-flow energy change.

For M = U/c, beta^2 = 1-|M|^2, e parallel to M, and
A = beta I + (beta^2-beta) e e^T, each source pair uses
`d = A^-1(s_j-s_l)`, `x = k |d|`, `a = M.d/|d|`, and the exact angular
kernel `exp(i a x) [j0(x) + i a j1(x)]`. The common power factor is
`rho omega^2 / (8 pi c beta^4)`. The scalar/vector sphere integrals follow
[DLMF 10.54.2](https://dlmf.nist.gov/10.54.E2). Thus separated-source phases and
flow direction affect interference; even a single fixed volume source has the
Mach factor `1/(1-|M|^2)^2` relative to the same source in a stationary fluid.

Exactly coincident source sums, the coincident baseline and pair contractions
use double-double arithmetic. The separate `j0-1` and `cos(a x)-1` decrements
preserve subwavelength opposed-source radiation without subtracting rounded
order-one powers. This is not a certified relative enclosure at every null;
negative or nonfinite numerical power refuses. For P direct/image paths and G
distinct positions, preparation takes O(P^2) work and stores O(P + G^2) data.
Immutable evaluation is allocation-free O(P G + G^2) without trigonometry or
angular quadrature.
Explicit source-coherence factors add powers between factors and preserve
interference within each factor; they use peak spectral-line units, not PSD/Hz.
Source/pair/factor work bounds and cancellation precede result publication.

The `from_observer` adapter preserves the original source geometry and any
already admitted set of one to three perpendicular infinite ideal rigid or
pressure-release planes, or an ideal pi/m wedge. Mean flow is tangent to every
plane, hence parallel to the edge for a wedge. The complete
signed image group remains coherent within each original physical source and
each statistical factor. For reflection R, the total energy flux transforms
as `J(Rx)=R J(x)`, so whole-image power divided by the group order is power into
the physical fluid sector: `2^n` for n perpendicular planes and `2m` for a
pi/m wedge. The wedge adapter reuses the original observer's `WedgeImages`
affine maps and signs, including translated planes with different anchors.
The same admission allows m in 2..=128; mixed rigid/release walls require even
m, while equal wall conditions admit every supported order. The original
observer owns geometry and tangency admission; no new boundary approximation
or source feedback is introduced.
Source/pair/factor budgets include every image path before allocation, and
the cancellable adapter publishes only complete groups. Input vectors still
name only physical sources. The boundary scope also survives an empty source
array. Preparation retains the original analytic pair kernel and evaluation
allocates no additional storage.

The model assumes homogeneous, subsonic, lossless flow and individually compact
prescribed sources. It adds no force-source vortical wake, nonuniform flow,
absorption, wall layers, finite impedance, general diffraction, source feedback,
finite-source validation or measured-SPL authority.
`tests/convected_power.rs` compares integrated power with the original full
pressure/velocity/energy flux on two enclosing spheres and checks the stationary
limit, Mach scaling, weak directional pairs, coherence, rigid transformations,
source-owner binding, resource admission and cancellation.
`tests/convected_power_boundaries.rs` checks independent one/two-wall interference
formulas, original near-field energy flux integrated over physical hemispheres
and right-angle sectors at two radii, the stationary image-array limit including
three planes, original-source coherence, full image budgets and cancellation.
`tests/convected_power_wedges.rs` adds independent single-source chord sums,
physical 60-degree and 45-degree sector flux at two radii, stationary signed
image-array parity, rotated/translated wall anchors, complete 2m path budgets
through m=128, cancellation and original-source statistical coherence.
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

## Ideal acoustic reflection in uniform tangential flow

`PreparedConvectedMonopoleObserver::with_orthogonal_planes` admits one to
three perpendicular infinite ideal rigid or pressure-release planes only when
the original physical mean flow is tangent to every wall. Tangency is checked
within 64 machine epsilons in normalized directions; no flow is projected.
The complete 2^m direct/image paths retain the original material, flow,
frequency and clearance limits. Preparation bounds every path and folds four
pressure/velocity coefficients per physical source, preserving O(N) runtime
and source-factor work counts of N. Each physical source remains coherent
with all its images, including under independent-source statistics.

`with_wedge` and its cancellable counterpart extend the same convected
potential to ideal wedges of opening `pi/m`, `2 <= m <= 128`, with flow
parallel to the edge (or zero). The shared `WedgeImages` geometry supplies
the complete `2m` dihedral group. Reflection preserves both distance and
the mean-flow displacement projection when the flow is tangent to both
faces, so this is an exact image solution for the admitted geometry.
Both walls may be rigid or pressure release; mixed conditions require even
`m` so their image signs close consistently around the edge. Other angles
and normal flow are refused. Preparation costs O(mN) work and O(N) scratch;
the same four coefficients per physical source preserve O(N) warm evaluation,
source coherence, material and direct-path travel-time metadata. The retained
`reflection_wedge_order` and every output scope identify the admitted model.

The same uniform-flow perturbation energy flux `<p u> + U <E>` is formed from
the folded field. Both ideal boundary conditions give zero wall-normal flux
for tangential flow. Shear flow, finite impedance, general wedge diffraction,
enclosed rooms and source feedback remain outside this owner. G1/G3/G4
reflection regressions check independent convected Green-function derivatives
and polar image locations through `m=128`, both flow directions, rotated and
translated faces, zero normal energy flux, zero-flow parity with the original
compact owner, the right-angle corner limit, physical-source coherence and
clock, and domain/budget/cancellation refusals with deterministic replay.

## Measured three-dimensional compact-source interference

`PreparedCompactSourceObserver`, `PreparedHalfSpaceObserver`, and
`PreparedConvectedMonopoleObserver` now accept interior bins from the shared
`WelchCrossSpectrum`. Compact channels are source-major `[Q, Fx, Fy, Fz]` in
m³/s and N; uniform-flow channels contain one physical Q per original source.
The shared admitted-bin view supplies exact-frequency checking and checked
`sqrt(2 df)` conversion from one-sided density factors to peak bin factors.

Original pressure/velocity transfers propagate every measured factor
coherently before forming second moments, preserving cross-source, cross-axis,
volume/force, direct/image, molecular-dispersion and convective interference.
Images never become independent recorded source coordinates. Both direct and
image paths count toward work admission in the separate half-space owner;
folded orthogonal-plane paths keep the original prepared transfer cost.
Queries allocate nothing, enforce source-factor work budgets and poll
cancellation at source boundaries. Uniform-flow pressure work and transported
perturbation energy remain separately reported.

Results integrate one resolved bin. Welch leakage and resolution remain;
cross-frequency covariance, overlapping-segment independence, confidence
intervals and additional source-model authority are not inferred. Six native
public-consumer regressions pass for complex partial coherence, mixed-mechanism
cancellation, center/adjacent Hann bins, rigid/release reflection, uniform-flow
energy, refusal and replay. The nine original convected regressions also pass.
These runs use the actual production compact/image/flow/statistics algorithms
and original fluid snapshot; the unavailable molecular constructor is omitted
from the leaf and no replacement material behavior is supplied. The seventh
public regression covers original molecular relaxation and reflection but
remains unrun without the constellation. This is not full workspace validation.

### Spatial coherence and physical microphone-array readouts

`observation::array::PreparedHarmonicArray<D>` binds borrowed original 2D/3D
receiver ensembles to one caller-authored source clock and unique ordered
physical fluctuation coordinates. Equal factor counts or laws alone do not
establish correlation. Allocation-free bounded pair queries retain complete
p/p, p/u, u/p and u/u finite-window second moments. Covariance of time-averaged
signals is computed separately and directly from real/proper fluctuation
factors, without subtracting a large deterministic mean.

`prepare_readout` coherently projects actual pressure and velocity through
ideal first-order microphones using each receiver's original fluid impedance.
Signed gains and causal harmonic delays apply once; projection, phase and
exact-frequency grouping retain double-double coefficients before products.
Prepared observations delegate to the existing harmonic ensemble/window owner.
Joint processed-channel statistics require matching clock/factor identities
and explicit pair-work admission. No propagation law is duplicated.

Cross-location p/u products are not local energy flux. Processed channels are
pressure-equivalent signals, not new spatial receivers, loading ports,
radiation estimates, measured microphone calibration or transient pulse
filters. Harmonic continuation and the original source/propagation limits
remain explicit. Five focused native public tests pass for independent
time/ensemble quadrature, beamforming/directional nulls, velocity response at
a pressure null, weak covariance beside large coherent terms, real/proper
coordinates, resource admission and cancellation. They compile the production
observation modules, shared harmonic/math algorithms and original fluid-state
constructor/methods; only root error carriers are nominal. The numerical leaf
also compiles with warnings denied. Full workspace integration remains unrun.

## Same-step plate momentum with explicit jet SGS

The 3D slot-jet records the native splitter-plate momentum exchange from the
actual post-collision populations before fringe forcing. Canonical plate-link
selection excludes the nozzle. Ordinary and chunked recording share this step;
checkpoint version 3 refuses histories containing the old post-stream force
proxy, rather than combining two force definitions in a resumed spectrum.

`run_slot_jet_3d_with_sgs` and `MaterialSlotJet3d::run_with_sgs` explicitly
consume the existing Smagorinsky, WALE or Vreman boundary closure during
settling and recording, using one reusable workspace and same-step wall loads.
Model coefficients, filter width and viscosity-ratio limits remain authored.
Zero coefficient preserves the corrected ordinary trajectory. Recorded shear
dissipation means and eddy-viscosity maxima retain their native lattice units.
Material density, molecular viscosity and periodic-cell force conversion
retain their original owners. This adds no random forcing, unresolved acoustic
source, broadband-regime validation or absolute-SPL authority.

`slot_jet_3d/tests.rs` covers independently assembled post-collision impulse,
fluid-plus-wall momentum balance, zero-closure bits, nonzero closure loads,
chunk partition equivalence, source/closure refusals and material SI scaling.

## Retarded reflection of finite moving compact sources

`curle3d::moving::trajectory::{IdealReflectionObserver,PreparedReflectedSourceTrajectory}` carries the complete 2^m image group of one through three perpendicular ideal rigid/pressure-release walls into moving-source pressure. The original microphone, medium and walls share a stationary fluid-rest frame. Virtual receivers implement `G(x,H y)=G(H^-1 x,y)`, retaining the original source position, velocity, acceleration, force and force rate without cloned trajectories or a second interpolation. Every path has its own retarded emission time and signed pressure coefficients; source strengths remain on their original clock.

A shared observer owns the complete bounded image geometry. Each borrowed source owner retains one arrival interval per image, reuses the original quintic and safeguarded root solver, and requires complete common receiver coverage. Outward Bernstein subdivision admits source confinement within every physical half-space and clearance from all virtual receivers; work/depth limits and cancellation remain explicit. Queries allocate nothing, preserve the represented-time residual policy, and never drop unsupported image histories. Path zero retains original direct-path arithmetic.

G1/G3 regressions compare accelerating source images with independently differentiated retarded Green potentials, stationary harmonic corner interference, each physical delay, rigid/release boundary conditions, coordinate covariance, whole-span wall crossing, budgets and cancellation. This is prescribed compact-source pressure, without finite-impedance walls, general wedges, enclosed rooms, moving receivers/walls, mean flow, absorption, radiation reaction or inferred source-model/experimental authority.

### Coherent harmonic source realizations and physical waveform output

`observation::array::PreparedHarmonicArray::realize` applies explicit caller-supplied real/proper-complex coefficients to the original mean/factor fields. Each named physical coefficient is shared across every receiver and frequency and stays fixed throughout a realization. Individual draws are finite but need not have unit magnitude; there is no normalization, extra sqrt(2), conjugation, RNG, inferred phase independence or per-sample redraw. Original source-clock/factor identities, receiver geometry and material snapshots remain attached. Preparation is bounded by channels, factors, lines and the squared per-channel combined line count; cancellation/refusal publishes no partial realization.

Mean and factor coefficients combine coherently in DD before sampling. `PreparedHarmonicArrayRealization::realize_readout` reuses the same draws with an existing matching microphone/beamforming readout; the original actual-p/u projection, material impedance, signed gains and causal delays are retained once. `PreparedHarmonicField::sample_field_into` and array `sample_fields_into` emit simultaneous pressure [Pa] and actual particle velocity [m/s] on the original integer/DD sample clock. Array output is sample-major and preserves receiver order; allocation-free absolute-index sampling is bit-stable under block partitioning. These are point samples of settled harmonic continuation, not finite-window means, transient attacks, a complete stochastic process or PCM normalization.

Shape/work/clock and all physical frequencies are admitted before output writes. Every line remains strictly below output Nyquist even after a zero draw, source cancellation or microphone null. Runtime numeric failure/cancellation may leave a prefix in caller output, which must be discarded. Focused tests recover spatial ensemble moments from eight realized waveform histories, verify physical directional/readout delays with nonunit complex draws, preserve weak pressure/velocity cancellation, refuse hidden out-of-band lines, and retain phase across sample indices beyond 2^53.

## Moving-microphone ideal-wall reflection on original trajectories

`MovingReflectionObserver` binds one original `PreparedKinematicTrajectory` microphone to one through three perpendicular ideal rigid/pressure-release walls. `PreparedReflectedMovingReceiverSource` retains all 2^m paths for an original source. Each borrowed image view applies the inverse wall transform to the receiver's position and derivatives at the physical observation time, while the original source polynomial, Q/F history and force frame remain unchanged. The existing moving-pair solver resolves every path's retarded emission time. Pressure remains the Eulerian source field; the path-specific observed source-clock derivative accounts for both motions and is never multiplied into pressure.

The existing continuous relative-velocity Doppler bound now transforms each receiver image's velocity hull, together with a conservative bound on the stored reflection operator's norm. Direct paths preserve their existing arithmetic. Original source and microphone paths receive whole-span half-space admission; the source remains strictly inside and exactly represented tangential receiver motion can lie on an ideal wall. The original causal Bernstein pair bounds admit every image's clearance. Pair and node limits are aggregate across all images, root caps remain per path, finite source coverage is required at both receiver endpoints, and cancellation publishes no partial owner or pressure sum.

The bounded image owner borrows histories without copying or recompiling them, and queries allocate nothing. Ten original physics/audio regressions cover independent accelerating image-source Green potentials, both moving clocks, reflected Doppler bandwidth, exact release-wall pressure cancellation, the fixed-microphone limit, coordinate covariance, complete histories, work limits and transactional rendering. Together with 30 existing controls their prior focused original-source build passed. That prior run retained narrow material/FFT/error-wiring facades and was not full-workspace or encoder integration validation. The model remains prescribed compact-source pressure in a homogeneous lossless fluid with stationary ideal walls; moving walls, finite impedance, arbitrary wedge angles, enclosed rooms, absorption and source feedback are not added.

`MovingReflectionObserver::in_uniform_flow` and `wedge_in_uniform_flow` also
admit stationary laboratory walls when the declared uniform mean flow is tangent
to every wall. The laboratory wall points are translated to the original
trajectory frame once, preserving stored normals and complete image groups.
Tangency is checked against normalized flow direction, so a small nonzero
normal flow cannot masquerade as a stationary fluid-frame boundary. Original
source and microphone paths must share the exact same `UniformFlowFrame`.
Retarded roots, whole-path clearance, reflected Doppler bounds and source
strength histories use the existing fluid-frame owners without recompilation
or extra query allocations. Pressure is the convected Eulerian field sampled at
the moving microphone. Tests compare accelerating orthogonal and wedge images,
monopoles and dipoles against independently differentiated laboratory convected
Green potentials, verify observed clock rates, and refuse incompatible frames,
normal flow and incomplete work. This adds no sheared flow or moving-wall model.

`IdealReflectionObserver::new_wedge` and `MovingReflectionObserver::new_wedge`
now share the harmonic owner's ideal pi/m dihedral admission, m in 2..=128.
All 2m paths are retained; equal rigid/release conditions admit every supported
order, and mixed conditions require even order. There is no geometry snapping
or reflection-order truncation. Original perpendicular constructors preserve
their transform arithmetic. The shared group prepares each inverse affine
receiver map once, so moving queries perform fixed work per path rather than
repeating m reflections or constructing image trajectories. An anchored map
retains physical translation and force/velocity orientation. Interval bounds
use the stored coefficients, and the outward norm bound
`sqrt(1 + ||A^T A - I||_infinity)` accounts for their rounding in Doppler admission.
The cached inverse Jacobian's transpose also carries monopole particle velocity
back to physical receiver axes, including wall parity. The same transform reaches
stationary reflected fields and both retarded-query and emission-endpoint paths;
signed pressure-work flux follows the transformed pressure and velocity together.

Preparation retains O(m) image metadata and borrowed pair views. Per-source
path caps, aggregate pair/node work, original source coverage and cancellation
apply to every image. The existing `fs-couple` reflected-pressure renderers
consume these owners unchanged, charging all 2m paths before output.
G1/G3 tests compare accelerating source/receiver paths against independent
polar source images (including dipoles and observed-clock rates), recover the
harmonic wedge limit and physical frame covariance, and refuse incompatible
wall signs, angles, incomplete budgets and cancellation. These finite ideal
wedge solutions do not add arbitrary-angle diffraction, finite impedance,
enclosed-room reverberation, moving walls or experimental qualification.

## Distributed three-dimensional splitter forces and coherent observation

`slot_jet_3d::distributed` adds an opt-in x/z partition of the existing splitter's
actual halfway-reflection links. `run_slot_jet_3d_patches` retains all three force
components ON THE FLUID, measured from the same ordinary or explicitly selected
SGS collision that streams each recorded step, before fringe forcing. The
original central-moment trajectory, seed, molecular viscosity and diagnostics
remain authoritative. Each patch bounds both plate faces, edges and diagonal
links; periodic-seam coordinates are represented around their selected solid
cell. Patch count, force-record bytes and native probe bytes are explicitly
admitted before the associated allocations.

`MaterialSlotJet3d::run_patches` binds those forces to the original material and
physical clock, using its existing `rho * dx^4 / dt^2` force scale and `dx` geometry
scale. It retains Fx, Fy and Fz without span/area normalization. The returned
record's `estimate_welch` feeds sample-major data directly to the original joint
Welch owner as source-major `[0, Fx, Fy, Fz]` channels; no full-record transpose or
monopole-zero history is allocated. Sample count derives from the private force
buffer. The resulting spectrum retains its material, positions and radii.

`observe_bin` applies the existing compact-source-array pressure/velocity Green
transfer to each synchronized Welch factor before calculating second moments.
It preserves measured cross-patch and cross-axis phase and coherence, including
radiation from separated opposing forces whose aggregate vanishes. Every patch
must satisfy caller-authored `k*radius` and radius/range limits at that observer
and frequency. Output is bin-integrated rather than PSD/Hz. No acoustic relative-
error bound follows from those geometric admissions, especially near destructive
interference.

This is resolved stationary-wall dipole radiation from an authored distribution
within one periodic flow cell, in the same homogeneous stationary lossless
material and free space. It does not qualify an isolated physical radiator,
replicate a periodic array, model finite-body scattering/source feedback, resolve
volume quadrupoles, or synthesize stochastic SGS sound. Existing aggregate
records and checkpoint version 3 are unchanged. Refused/cancelled runs expose no
partial force record; native patch exchanges publish only with accepted fluid.

Five new aeroac unit groups cover canonical x/z partition and periodic geometry,
record/probe budgets, native aggregate force sign, original material SI scaling,
and propagation/coherence of opposing spatial forces including Fz. Together with
three native LBM patch-load groups, these targets typecheck through full native
fs-aeroac production and unit sources after concurrent boundary changes were
merged. These runtime tests and turbulent acoustic benchmarks were not executed
in this session; no measured-spectrum or SPL qualification is claimed.

## Radiated power in original ideal reflecting domains

`PreparedCompactRadiatedPower::from_observer` retains the original folded
perpendicular-plane or pi/m-wedge image group. `from_half_space` binds the
separate direct/reflected owner without rebuilding a field observer. Both use
the original physical source positions, material snapshot, wall signs and
Cartesian force-reflection matrices. The complete image pressure transforms by
the wall's sign, so its intensity is invariant under the finite reflection
group. The existing analytic full-sphere integral divided by that group's order
is therefore the power into the physical fluid solid angle. Free-space
preparation retains its original arithmetic and full-sphere meaning.

For one monopole a height h above a plane, the result is
`P/P_free = 1 + s*sinc(2kh)`, where s is +1 for a rigid wall and -1 for pressure
release. Normal and tangential dipoles retain their different reflection
parities; volume/force cross terms remain. Existing DD moment decrements retain
weak pressure-release radiation without subtracting two rounded unit powers.
No clipping or positivity repair is applied.

Every input/coherence-factor entry still names one physical source. Images of
that source stay coherent even when different physical sources are declared
independent. Source/path metadata and scopes distinguish the physical domain.
Preparation budgets count all P direct/image paths and P² work before image
allocation; factor budgets include that same expanded pair count. Coherent
application reuses cached spherical moments without allocation or trigonometry.
Independent-source work is O(N g²) for N physical sources and group order g,
without allocating sparse factors; the original free-space O(N) path remains.
Cancellation/refusal publishes no partial owner or power.

G1/G3 cases in `tests/compact_power/reflection.rs` cover analytic monopole and
both dipole parities, tiny release nulls, original physical hemisphere flux,
rotated/translated geometry, physical wedge flux and corner equivalence,
source-factor identity, expanded budgets and cancellation/replay. The model
remains prescribed compact sources in a stationary homogeneous lossless fluid
with ideal infinite walls. Molecular absorption still refuses; finite wall
impedance, general diffraction, source feedback and experimental SPL authority
are not inferred.

## Stationary planar acoustic refraction

`curle3d::refraction::PreparedRefractedMonopole` observes one prescribed physical
volume source across an infinite plane between two original `FluidAcousticState`
snapshots. A bounded monotone Fermat solve locates the crossing. Sound speeds set
the two ray directions and travel time; both densities and speeds set the oblique
pressure transmission and normal-power partition. Refracted spreading uses both
principal ray-tube widths, not the sum of the ray lengths. The source-clock
exp(-i omega t) transfer is `-i omega rho1 Tp exp(i omega T)/(4 pi D)`.

The API requires explicit minimum k*h on each side, minimum direction cosine on
both legs, crossing half-width tolerance and root-iteration cap. These are
caller-selected high-frequency/grazing controls, not certified bounds on the
geometrical-acoustics error. Sources and receivers must be strictly on their
declared sides. Cancellation, exhaustion and unrepresentable fields refuse
without returning a partial owner. Preparation and application allocate nothing;
application reuses four cached complex coefficients without transcendental work.

Output uses the existing coherent-field and harmonic-line types, including
integration with the original finite-window observation owner. Velocity is the
leading-order transmitted `p n2/Z2`; reactive near-field induction is absent.
The model returns the selected transmitted geometrical arrival only. It does
not add head/lateral waves, evanescent transmission, critical diffraction,
reflected source feedback, flow, absorption or experimental source authority.
G1/G3 tests in `tests/refraction.rs` check Snell fields, normal energy partition,
equal-speed spherical pressure, reciprocal source/receiver exchange, independent
stationary-phase Hessians, rotation/scaling, physical clock integration and
explicit admission/work refusals.

## Full-wave transmission of a physical spherical source

`curle3d::refraction::full_wave::PreparedFullWaveMonopole` prepares the transmitted
pressure and all three particle-velocity components of one prescribed volume
source across the original stationary fluid plane. The outgoing Sommerfeld
spectrum retains propagating, critical and evanescent components. Pressure and
velocity derive from the same Green function and original densities/speeds;
near-field velocity and reactive intensity are retained. Equal original sound
speeds use the exact scaled spherical field, including receiving-density
velocity scaling. Different but nearby speeds still use distinct spectra.
The scalar identity and outgoing exp(-i omega t) branch follow Chew, Purdue
ECE 604 lecture 35, equations (35.1.11)--(35.1.14):
<https://engineering.purdue.edu/wcchew/ece604s20/Lecture%20Notes/Lect35.pdf>.

Branch-regularizing substitutions split the spectrum at both original fluid
wavenumbers. A bounded adaptive Gauss--Kronrod rule integrates pressure/radial/
normal velocity together, with DD coherent reductions. Explicit evaluation and
depth limits bound temporary storage and work; cancellation or exhaustion
publishes no owner. Monotone phase and attenuation variation must be resolved
before embedded error estimates can pass. An analytic omitted-tail bound is
reported separately from estimated finite quadrature error. Neither includes
geometry, Bessel or floating-point error; this is not an interval certificate.
Diagnostics report maximum q*R against the existing Bessel owner's documented
accuracy range. There is no implicit high-frequency fallback or spectral loss.

The cached four-coefficient transfer applies original source-clock amplitudes
without allocation or transcendental work and feeds existing coherent-field and
harmonic-window observers. The model covers strictly positive source/receiver
clearances on opposite sides of one infinite lossless plane. It does not model
flow, absorption, finite surfaces, layers or source feedback. Force sources
have their own complete tensor transfer described below.
`tests/full_wave_refraction.rs` covers independent high-precision direct-spectrum
references, important evanescent contributions, exact reactive equal-speed
fields, reciprocity, pressure-gradient velocity, ray-tube asymptotics, geometry,
refinement, original observation clocks and explicit work/cancellation refusals.

## Full-wave refraction of physical force dipoles

`curle3d::refraction::full_wave::PreparedFullWaveForce` transports an actual
three-dimensional force **on the incident fluid** [N] across the same original
fluid plane. This supplies the transmitted field of Curle/structural compact
forces without treating them as volume sources or multiplying by density a
second time. All force components retain their complex source-clock phase.

Pressure uses the analytic source derivative of the transmitted Green
function; velocity uses its complete mixed source/receiver derivative.
Source-normal factors use the incident vertical wavenumber, receiver-normal
factors the receiving wavenumber. The two mixed tensor entries remain distinct
when speeds differ. Tangential velocity remains present at dipole pressure
nulls, including both transverse channels on the normal axis. Equal original
speeds use exact derivatives of the scaled spherical Green function, including
the receiving-density velocity correction; nearly equal speeds still use
their distinct spectra.

Seven complex cylindrical coefficients share the original outgoing branch
transformations, bounded adaptive Gauss--Kronrod nodes, DD reductions and
phase/attenuation resolution controls. The original monopole remains a
three-channel specialization. `FullWaveForceLimits` controls the pressure
vector's Euclidean norm [Pa/N] and velocity tensor's Frobenius norm [m/s/N],
so multiplying their estimated errors by the force norm covers any prescribed
orientation. Relative transfer tolerances do not promise relative field error
at a directional null. Diagnostics report each radial/normal pressure and each
tt/ss/tn/nt/nn velocity error separately. Analytic degree-one/degree-two spectral
tail envelopes retain near-interface gradients; their logarithmic evaluation
prevents premature underflow. Strongly attenuated physical integrands also use
logarithmic scaling when ordinary intermediates cannot represent the field.
Quadrature remains estimated, not certified; Bessel, geometry, coordinate
rotation and floating-point errors remain outside both reported quantities.

Preparation stores twelve Cartesian complex coefficients. Applying a force or
feeding `HarmonicFieldLine<3>` is allocation-free with fixed-order DD sums and
no trigonometry, quadrature or numerical differentiation. Explicit budgets and
cancellation return no partial owner. The source remains prescribed, compact,
stationary and three-dimensional; finite source size, moving/absorbing fluids,
finite or layered interfaces, mutual fluid loading and measured-SPL authority
remain outside this model.

`tests/full_wave_force_refraction.rs` pins all seven channels against independent
35/50-digit direct-spectrum references, the homogeneous density-jump limit,
source-derivative signs, receiving-fluid momentum, pressure-null transverse
motion, rotation covariance, tensor-norm refinement, physical observation clocks
and work/cancellation refusals. A balanced-half-decay unit regression covers
representable force gradients after raw evanescent exponential underflow.

## Reflected and coherent incident-side full-wave observations

`curle3d::refraction::full_wave::reflection::PreparedReflectedFullWaveMonopole`
places both original source and observer in the incident fluid at the existing
`PlanarFluidInterface`. It retains the complete spherical reflection spectrum,
including both critical circles, total internal reflection, evanescent and
reactive contributions. `evaluate_reflected` returns only the returning field;
`evaluate_total` coherently adds the original homogeneous direct field before
forming pressure/velocity intensity. Original phases and interference cross
terms are retained. Both views expose native harmonic lines without another
propagation delay. A coincident source/receiver admits a finite reflected field
but refuses the singular direct/total field; it does not infer source feedback.

With outgoing `zeta_j=sqrt(k_j^2-q^2)` and incident-side clearances a,b, the
pressure reflection factor is `R=(rho2*zeta1-rho1*zeta2)/(rho2*zeta1+rho1*zeta2)`.
The shared Weyl weight is `q*dq*R*exp(i*zeta1*(a+b))/zeta1`; its pressure scale
is `omega*rho1/(4 pi)`, while radial and normal velocity follow analytic
receiving derivatives. Returning normal velocity retains its minus sign.
Equal original sound speeds use an exact density-jump spherical image,
including induction; identical fluids return exact zero reflected field.

The original critical-branch map and bounded Gauss--Kronrod engine are shared
with transmission. The incident-grazing quotient is cancelled analytically
before evaluation. Log-scaled admittance products and spectral weights avoid
premature overflow/underflow. Phase and attenuation resolution are supplemented
by a bound on within-panel reflection-factor variation, preserving narrow
critical transitions when original speeds nearly coincide. Analytic omitted
tail bounds use `|R|<=1` and the complete reflected clearance `a+b`.

`FullWaveMonopoleLimits` and diagnostics control the reflected transfer.
The coherent total inherits its absolute numerical error estimates but has no
relative-accuracy guarantee at direct/reflected cancellation. Quadrature errors
remain estimates; roundoff, geometry and Bessel errors remain outside those
estimates and the analytic tail bounds. Preparation is bounded/cancellable;
application uses cached coefficients without allocation or transcendental work.
This is a stationary, lossless, infinite two-fluid plane with prescribed compact
sources; absorption, finite edges, layers, moving fluids, loading and measured
SPL authority are outside the model.

`tests/full_wave_reflection.rs` covers independently refined direct-q spectra
for both critical orderings and nearly equal speeds; exact image/interference
fields; interface pressure and normal-velocity continuity against the original
transmission owner; momentum, reciprocity, rotation, finite reflected self-field,
native observation clocks and bounded refusals. Private regressions exercise
critical endpoints, raw-exponential underflow and narrow reflection transitions.

### Reflected physical force fields

`PreparedReflectedFullWaveForce` uses analytic source and receiver derivatives
of that same reflected Green function for an original force ON fluid [N].
All seven cylindrical pressure/velocity channels share the branch maps,
reflection coefficient, resolution controls and adaptive nodes. The mixed
velocity coefficients have opposite signs, and the normal-normal coefficient
retains the returning-wave minus sign. At equal speeds this is the exact image
force `R*(F-2*n*(n dot F))`, including transverse/induction motion. The original
homogeneous force owner supplies the direct contribution in the actual source
geometry, and direct/reflected amplitudes are added before intensity is formed.

Original `FullWaveForceLimits` separately bound reflected pressure-vector and
velocity-tensor norm estimates per newton. Degree-one/two tail envelopes retain
all physical derivatives and the incident-fluid density. The total field
inherits those absolute estimates but no relative guarantee at an interference
null. The two cached 4-by-3 Cartesian matrices serve arbitrary complex force
vectors with fixed-order reductions and no allocation or further integration.
Reflected and total harmonic methods use the original source clock. Coincident
source/observer positions retain reflected fields and refuse singular total
queries; this does not introduce source loading or feedback.

`tests/full_wave_reflected_force.rs` covers all seven coefficients against three
independent refined spectra, rotated image-force signs, coherent total intensity,
source derivatives, incident-fluid momentum, transverse motion on axis, finite
reflected self fields, original clocks and bounded input/work refusals.

## Full-wave compact-source arrays and source coherence

`curle3d::refraction::full_wave::array::PreparedFullWaveCompactSourceObserver`
prepares the original volume and force transfers for multiple physical source
locations behind one stationary fluid plane. Each source retains the native
`CompactSourcePhasor` [Q,Fx,Fy,Fz] and its sixteen Cartesian pressure/velocity
coefficients. Sources lie strictly in the incident fluid, the observer strictly
in the receiving fluid. All propagation phase is already in the transfer;
there is no invented single ray delay for the full-wave spectrum.

`FullWaveCompactSourceLimits` bounds source allocation and total integrand
evaluations across every Q and F preparation, in addition to each kernel's
own numerical and work limits. Exact equal-speed kernels consume no spectral
evaluations. Exhaustion, invalid input or cancellation publishes no partial
array. An empty array is validated silence. Per-source diagnostics keep Q and
F units distinct: multiply monopole errors by |Q|, pressure-vector and
velocity-tensor errors by ||F||, and sum across sources for an absolute field
estimate. Relative source tolerances do not guarantee relative array error at
an interference null. All original numerical and physical scope limits remain.

Application reuses the cached 4-by-4 matrices in O(N), without allocation,
quadrature or transcendental functions, with fixed-order DD sums. Coherent
evaluation combines all source and Q/F amplitudes before intensity; native
`HarmonicFieldLine<3>` output feeds the existing finite-window observer.
Explicit independent sources and bounded low-rank coherence factors reuse
the original `StatisticsSum` owner for pressure power, full velocity covariance,
and signed active/reactive intensity. Unknown phase does not imply independence.
Factors retain peak-source units; one-sided source PSD in a bin of width df
requires factors of 2 df S. Results are bin/line moments, not PSD/Hz and not
an assertion of finite-window cross-frequency incoherence.

`tests/full_wave_refraction_array.rs` covers homogeneous density-jump reduction
with mixed Q/F sources, actual unequal-speed kernels, coherent cancellation,
independent versus factored source statistics, all p/u moments, native harmonic
lines, per-source diagnostics, shared work exhaustion, empty arrays and input
or cancellation refusals.

### Reflected and coherent total arrays in the incident fluid

The same prepared array owner adds `new_reflected` and `new_total`; the original
`new` remains transmitted-only. `FullWaveCompactFieldKind` and `scope()` retain
which physical field every coherent, independent-source, and low-rank
statistical result contains. Original sources and the observer are strictly in
the incident fluid for both new constructors. The reflection owners supply
each original Q/F column, including complete particle velocity and evanescence.

For total arrays, each source's direct and returned complex amplitudes are
combined before source statistics. Even explicitly independent physical sources
retain their own direct/reflected interference; an image is not an independent
noise source. Original aggregate work and allocation budgets still cover every
Q/F preparation. Evaluation reuses the same fixed-order DD and StatisticsSum
paths without allocation, new quadrature, or added propagation delay.

Reflected arrays admit finite source-coincident observations; total arrays
refuse a singular direct self-field at preparation even if its eventual source
amplitude would be zero. Numerical diagnostics control the reflected transfer
and do not promise relative accuracy at total-field interference nulls.
`tests/full_wave_reflected_array.rs` covers exact density-jump images with mixed
Q/F sources, all covariance/intensity components, independent and factored
coherence, unequal-speed composition, shared work exhaustion, singular versus
finite self-fields, harmonic clocks, cancellation and invalid inputs.

## Full-wave reflection from a finite fluid layer

`curle3d::refraction::full_wave::layer::reflection::PreparedLayeredReflectedFullWaveMonopole`
observes the original incident-side source and receiver at one
`PlanarFluidLayer`. Its returning Green function retains every coherent
internal reflection, all three critical circles, evanescence, and the full
pressure/particle-velocity field. The original prescribed-source clock is
retained; no ray delay, truncated image train, artificial damping, or source
feedback is introduced. `evaluate_reflected` and `evaluate_total` retain
separate cached responses, and native harmonic-line methods carry those same
amplitudes. A finite reflected self-field is admitted; its singular direct
total field is refused.

The original transmission denominator supplies the effective reflection
coefficient through a signed numerator from the same continuity equations.
Shared critical-point maps, density normalization, scaled trigonometric
functions, and denominator-variation guards resolve internal resonances without
forming growing transfer matrices. Analytic incident-grazing cancellation and
logarithmic numerator/denominator scaling preserve representable fields. The
tail starts beyond every critical circle, where `|R|<=1`, and uses the exterior
return clearance `a+b`; layer thickness is not incorrectly added to that decay.
This bound is not asserted for intermediate evanescent regions.

Original monopole limits and layer diagnostics control the reflected transfer.
Total observations inherit absolute reflected error estimates without a
relative guarantee at interference nulls. A zero-thickness layer or one matched
to its receiving half-space delegates to the original interface owner exactly.
Applying prepared responses uses fixed work without allocation or quadrature.
The original slow-layer guided-pole refusal, lossless/static/infinite-plane
scope, and absence of interval or measured-SPL authority remain.

`tests/full_wave_layer_reflection.rs` covers independently refined 35/50-digit
direct-q references, zero-thickness and matched-interface reductions, a shifted
original interface, coherent equal-speed image trains and total intensity,
momentum, reciprocity, original clocks, finite self-fields, and bounded input,
work and cancellation refusals. Private tests cover every grazing branch,
reflection zeros and physically scaled fields beyond raw exponential underflow.

### Original forces reflected by the finite fluid layer

`PreparedLayeredReflectedFullWaveForce` in the same module extends that original
returning Green function to a prescribed three-component force ON fluid [N].
Its two pressure-vector and five velocity-tensor coefficients reuse the exact
analytic derivative expressions of the single-interface reflection owner;
there is no finite-difference force or pressure/impedance velocity substitute.
The opposite mixed signs and negative normal-normal term preserve the original
normal-force image reversal through every coherent internal return.

Two cached Cartesian 4-by-3 transfers serve reflected and coherent total fields
at the original source clock. The native layer source/receiver restrictions,
shared quadrature and resonance rules, exact first-interface reductions, and
finite reflected versus singular total self-field distinction are unchanged.
Force limits bound pressure-vector and velocity-tensor norms per newton; layer
diagnostics retain those units and the reflected-only relative-error scope.
Application uses fixed-order DD accumulation without allocation, quadrature,
finite differences or transcendental work. No source feedback is inferred.

`tests/full_wave_layer_reflected_force.rs` pins all seven channels against the
independently refined direct-q oracle, compares rotated equal-speed image-force
trains and coherent total intensity, verifies analytic source derivatives of
the monopole, retains on-axis transverse motion, and covers matched-fluid,
clock, self-field, invalid-input, work-limit and cancellation behavior.

### Prepared compact-source arrays around the original layer

`curle3d::refraction::full_wave::layer::array::PreparedLayeredFullWaveCompactSourceObserver`
composes the native layered Q/F owners for multiple original source positions.
`new`, `new_reflected`, and `new_total` respectively select transmission beyond
the second face, the returning incident-side field, or coherent incident-side
direct plus return. The original `PlanarFluidLayer` supplies all three material
snapshots and both face clearances without reconstructing or renormalizing its
plane. No observation inside the layer or guided-mode solver is inferred.

Each source retains one cached Cartesian 4-by-4 transfer and separate native
monopole/force diagnostics. `FullWaveCompactSourceLimits` bounds allocation,
per-mechanism integration and aggregate Q/F evaluations across the complete
array. Only an actual native interface reduction with equal exterior speeds
has zero quadrature work; equal speeds alone do not eliminate a finite density
layer. Empty arrays are validated silence. Total construction refuses any
singular direct source/observer coincidence, while reflected construction
retains finite returning self-fields.

Interface and layer arrays share the same private borrowed `ArrayField` and
`Transfer` implementation: fixed-order DD accumulation, coherent source phases,
independent-source statistics, bounded low-rank factors, full velocity covariance,
and signed active/reactive intensity. Application adds no allocation, quadrature,
transcendental work or dynamic dispatch. Within each original source, direct
and returned amplitudes interfere before any independence assumption across
different sources. Field-kind metadata and result scopes distinguish all three
observations; reflected numerical diagnostics retain their absolute-only
meaning near a total-field or array interference null.

`tests/full_wave_layer_array.rs` covers collapsed/homogeneous and invisible-face
limits, actual finite-layer Q/F composition and all moments, path/source coherence,
global work exhaustion, equal-speed density-layer work, finite versus singular
self-fields, original clocks, explicit domains, empty arrays and cancellation.

`PreparedLayeredFullWaveCompactSourceObserver::evaluate_welch_bin` now
propagates original source-major Q/F records through its native layer maps.
The companion `evaluate_welch_bin_with_quadrupoles` consumes one synchronized
record containing compact `[Q,Fx,Fy,Fz]` blocks followed by physical stress
`[xx,yy,zz,xy,xz,yz]` blocks. Both banks must share all three original fluid
states, first-face point and oriented normal, thickness, receiver, exact bin
frequency and transmitted/reflected/incident-total field. Interior stress
observations cannot be paired with these exterior compact-source views.
One common DD pressure/velocity sum per original factor retains Q/F/stress,
cell, tensor, direct/returned and internal-layer interference before moments.
The original peak-bin conversion occurs once. The existing compact factor loop
is shared unchanged with the fluid-interface observer; warm application adds
no FFT, quadrature, covariance expansion, delay or waveform reconstruction.
Combined source/segment work and cancellation are checked before reduction.
Original channel ordering, disjoint-source responsibility, window/aliasing
limitations and the nontrapping lossless-layer scope remain explicit.
`full_wave_layer_mixed_welch` specifies native finite-layer moments in all
three views, absent/homogeneous limits, measured phase cancellation with
surviving velocity, and material/frame/thickness/clock/work refusals.

## Full spherical-wave transmission through a finite fluid layer

`curle3d::refraction::full_wave::layer::PreparedLayeredFullWaveMonopole`
extends the prescribed original volume source to two parallel faces enclosing
one homogeneous lossless fluid layer. `PlanarFluidLayer` retains all three
original `FluidAcousticState` snapshots, a first-face point, a unit receiving
normal and the finite nonnegative thickness. The source and receiver lie
strictly outside opposite faces. The Green function sums every coherent
internal return and retains propagating, critical and evanescent incident
components, including frustrated transmission through an evanescent layer.
Pressure and both cylindrical velocity components share that Green function;
the returned Cartesian field retains active and signed reactive intensity.

With outgoing `zeta_j = sqrt(k_j^2-q^2)`, `Y_j=zeta_j/rho_j`, and original
exterior clearances a and b, the transmitted source spectrum is
`p/Q = omega/(2 pi) integral q J0(qR) exp(i zeta_1 a+i zeta_3 b)/L dq`, where
`L=(Y_1+Y_3)cos(zeta_2 h)-i h(zeta_2^2/rho_2+rho_2 Y_1 Y_3)sinc(zeta_2 h)`.
This follows pressure and normal-velocity continuity, as in Mei, MIT Wave
Propagation, chapter 3, equations (1.18)--(1.25):
<https://web.mit.edu/1.138j/www/material/chap-3.pdf>.
The sinc form is regular at internal grazing. Multiplying the denominator by
`exp(-Im(zeta_2) h)` removes growing hyperbolic intermediates; the matching
physical decay remains attached to the exterior spectrum. A common vertical
wavenumber factor is cancelled analytically when all original speeds agree.
Density scaling and logarithmic weighting preserve representable attenuated
fields. Unrepresentable original coefficient ratios or distinct wavenumbers
rounded to the same value refuse explicitly.

The shared bounded Gauss--Kronrod engine splits at every distinct original
fluid wavenumber. It requires phase/attenuation resolution and a conservative
bound on within-panel denominator variation before accepting any embedded
quadrature estimate. A narrow layer resonance cannot be accepted merely
because the embedded samples both miss it. The all-evanescent omitted tail has
an analytic envelope with total clearance `a+h+b`, preserving layer decay;
pressure and velocity tails are reported separately in physical units. These
are numerical estimates plus analytic spectral-tail bounds, not interval
certificates: geometry, Bessel, floating-point and coordinate errors remain
outside the estimates. There is no implicit ray fallback or added damping.

A nonzero layer must satisfy `c_layer >= min(c_incident,c_receiving)`.
Slower-than-both layers may have real guided poles and require a separate
pole/residue or physical-loss treatment; this owner refuses that class.
Zero thickness delegates the original two-fluid interface exactly, and three
identical fluid snapshots delegate its exact spherical field. Preparation
retains explicit work/depth limits and cancellation; failure publishes no
partial owner. The cached four-coefficient transfer applies complex peak Q
or produces a native harmonic line without allocation, quadrature, or a second
propagation delay. Absorption, multiple layers, finite lateral edges, moving
fluids and source feedback are outside this model.

`tests/full_wave_layer_refraction.rs` covers independent 35/50-digit direct-q
references for unequal exteriors, an evanescent barrier and a narrow coherent
resonance; an independent spherical-image series with equal speeds and unequal
densities; exact and vanishing-thickness limits; receiving-fluid momentum,
reciprocity, rotation, refinement, observation clocks, work and cancellation.
Private regressions check internal grazing, large evanescent layer depth and a
resonance that passes the phase-span test but fails denominator resolution.

### Physical force sources through the finite layer

`PreparedLayeredFullWaveForce` applies the same coherent finite-layer Green
function to a prescribed original compact force ON the incident fluid, in
newtons. The source derivative is `i/(omega rho_incident) grad_source(p/Q)`;
analytic receiving derivatives retain all five cylindrical velocity-tensor
entries, including transverse motion at a pressure null. The seven pressure
and velocity channels share the parent's critical branches, physical spectral
weight, resonance resolution and adaptive nodes. No finite-difference dipole,
single ray, truncated return series or `p/Z` velocity substitutes for the tensor.

The original `FullWaveForceLimits` separately control the pressure-vector norm
and velocity-tensor Frobenius norm per newton. Derivative-aware analytic
evanescent tails retain the layer thickness and original density factors;
reported quadrature errors remain estimates rather than interval certificates.
The prepared 4-by-3 Cartesian matrix applies arbitrary complex force vectors
without allocation, quadrature or transcendental evaluation and exposes native
harmonic lines with the existing clock. Zero thickness and identical fluids
delegate the original single-interface force owner. The parent's exterior-only,
lossless, nontrapping and prescribed-source scope remains unchanged.

`tests/full_wave_layer_force.rs` checks all seven tensor channels against
independently refined direct-q spectra, a coherent spherical-image force series,
source derivatives and receiving-fluid momentum. It also checks exact collapsed
limits, axial transverse motion, harmonic observation, zero/nonfinite forces,
work exhaustion, cancellation and original geometry/domain refusals.

## Resolved volume-stress quadrupoles

`curle3d::quadrupole::PreparedQuadrupoleObserver` propagates prescribed symmetric
Lighthill stress moments in the original `FluidAcousticState`. Each cell supplies
`Q_ij = integral_cell T_ij dV` in N m, ordered xx, yy, zz, xy, xz, yz; shear
entries are physical tensor entries, not pre-doubled. The complete trace and
cross stresses are retained. Source cells keep their own phase and position;
the cell quadrature is compact, not the entire distributed source region.

For the existing peak-phasor convention exp(-i omega t), the free-space Green
function is `G = exp(ikr)/(4 pi r)`. Analytic Hessians give `p = Q_ij d_i d_j G`;
the corresponding third derivatives give all three components of
`u = grad(p)/(i rho omega)`. This includes transverse and reactive near fields,
including velocity at pressure-null directions. The isotropic identity
`Q = q I => p = -k^2 q G` retains the physically relevant trace. Source physics
follows Lighthill (1952), DOI 10.1098/rspa.1952.0060; a resolved-flow application
is Sarkar & Hussaini, NASA-CR-191543 (1993),
<https://ntrs.nasa.gov/citations/19940017407>.

`resolved_reynolds_moment_n_m` computes instantaneous `rho u_i u_j dV` from
actual density, volume and velocity in the propagation frame. The original
stress history must be integrated or Fourier transformed on its original
clock; squaring velocity phasors at the same frequency is not a source model.
The helper infers no mean subtraction, detrending, pressure/density correction,
viscous stress or SGS completion. Full prescribed Lighthill moments can instead
be supplied directly. Surface-force and volume-source decompositions remain
the caller's physical responsibility; no source is added automatically.

Preparation caches 24 complex coefficients per cell. Application is O(N),
allocation-free and contains no transcendental work; fixed-order double-double
reductions retain coherent cancellation. The separate observer leaves existing
monopole/dipole transfer storage unchanged. `evaluate_harmonic_line` feeds the
original pressure/vector-velocity finite-window observer; equal-frequency
lines combine before quadratic statistics. Explicit coherence modes use the
existing bounded statistics reduction without inferring independence or random
phase. Source, geometry, finite-value, work and cancellation admission publish
no partial prepared observer or field.

`PreparedQuadrupoleObserver::with_orthogonal_planes` folds the complete 2^m
image group of one to three perpendicular ideal walls into those same 24
coefficients. `with_wedge` uses the existing 2m dihedral group for opening pi/m;
equal rigid/release walls are admitted for m in 2..=128 and mixed walls require
even m. Each original source stays strictly in every physical half-space;
receivers may lie on a wall. The original direct path, source positions,
material, frequency and source-count semantics are retained.

For a source image H and pressure character s, the physical stress image is
`s H Q H^T`, including trace and all cross stresses. Preparation equivalently
uses the original tensor at virtual receiver `H^-1 x`, then pulls the complete
velocity back by `(H^-1)^T`. This transpose matters for wedge rotations.
Images remain coherent within every physical source and every recorded factor;
pressure and particle velocity are summed before energy products. Preparation
is bounded by physical sources times all image paths, checks cancellation per
path/source, and uses O(N) scratch. Evaluation retains allocation-free O(N)
work with no tensor transformation or transcendental work.

G1/G3/G4 `quadrupole_reflection` tests compare with independently transformed
source tensors at oblique walls, check each tensor basis against rigid/release
conditions in a mixed three-wall corner, and verify rotated wedge faces,
harmonic/statistical readouts and geometry/work/cancellation refusals.
These ideal infinite boundaries do not model finite impedance, edge or general
wedge diffraction, closed-room reverberation, absorption, motion or source
feedback, and do not validate unresolved turbulent sound or measured SPL.

`evaluate_welch_bin` consumes the original synchronized `WelchCrossSpectrum`
factors with six source-major N m channels per prepared cell, ordered xx, yy,
zz, xy, xz, yz. Every cell and tensor component is summed coherently within
one recorded factor before p/u second moments are formed. The original
`sqrt(2 df)` peak conversion is applied once, yielding bin-integrated moments
rather than PSD/Hz or fitted tone amplitudes. Native window/detrending, leakage,
finite-record and source limitations remain; overlapping segments imply no
independence, and distinct bins imply no cross-frequency coherence. Frequency,
channel count and cell/factor work are admitted before allocation-free cached
propagation, with cancellation at least every 256 cells and between factors.
The G1/G3/G4 `quadrupole_welch` fixtures cover independently prescribed tensor
phasors and all p/u moments, Hann bin power, instantaneous velocity-to-stress
frequency doubling, cross-cell interference, pressure-null transverse motion,
and shape/frequency/work/cancellation refusals.

`quadrupole::recorded::QuadrupoleWelchStream` binds the fixed cell positions and
six physical stress channels to an explicit epoch, dt and timestamp-error
budget. It consumes instantaneous accepted stress arrays through the original
`WelchStream::push_with`; no velocity-spectrum multiplication, new estimator,
source independence, resampling or implicit decimation is introduced. Source
geometry/order, timestamp and record-length checks precede the atomic Welch
push. Timestamps derive from the epoch and integer index using DD arithmetic;
the admitted range must resolve adjacent samples even at a large epoch. The
caller still establishes that the supplied moments belong to that time and
that the source sampling resolves its bandwidth.

Preparation charges retained positions plus the original transform's entire
workspace/factor bank to the stated memory budget. Accepted samples allocate
nothing and retain only the native segment buffers, not full flow histories.
`RecordedQuadrupoleSpectrum` moves the completed factors without a copy and
binds observer preparation/evaluation to the same source positions and actual
interior bin. The original unused-tail, leakage, detrending, statistical and
point-quadrature limits remain; there is no deterministic waveform or
cross-frequency phase reconstruction. Four G0/G1/G4/G5 fixtures cover exact
streaming/batch factor agreement with overlap/tails, changed-source and clock
refusals, cancelled-transform retry, bound observer use, memory/shape admission
and timestamp-resolution refusal.

Kernel regressions cover units, the isotropic Helmholtz identity, independent
Green-function Hessian refinement, all-component linear momentum, rotation
covariance, far-field directivity/impedance, pressure-null transverse velocity,
and the actual quadratic source's DC/doubled-frequency content. Observer
fixtures cover coherent cancellation, explicit independent modes, native
harmonic-field bridging and admission/cancellation refusals. The ordinary
constructor is homogeneous stationary lossless point quadrature: cell-size
error requires spatial refinement. There is no moving-source correction,
mean-flow propagation, arbitrary-boundary scattering, feedback, unresolved
SGS sound or measured-SPL qualification.

`PreparedQuadrupoleObserver::new_molecular` and its cancellable variant retain
the original `RelaxationState` and its `k²=rho omega²/K(omega)`. The SAME complex
wavenumber supplies the Hessian, third derivatives and every ideal-wall/wedge
image. Pressure attenuation alone would give incorrect reactive velocity and
energy transport; momentum remains `grad(p)=i rho omega u`. Supplied N m stress
moments receive no additional density or volume factor. Classical shear and
conduction losses remain separate, and no attenuation filter is added. The
existing real-wavenumber kernel arithmetic is retained. Reference `r/c_eq`,
actual direct-path phase `Re(k) r/omega`, and broadband arrival are distinct;
neither exposed delay is applied again during harmonic evaluation.

`RecordedQuadrupoleSpectrum::prepare_molecular_observer` binds those material
coefficients to the actual interior Welch bin and original source layout.
Coherent, low-rank, recorded and finite-window readouts use the original cached
coefficients with unchanged O(N), allocation-free evaluation. Material refusals
remain typed separately from geometry, work and cancellation refusals. Focused
G1/G3/G4 `quadrupole_molecular` tests cover an independently eliminated caloric
reservoir and Green Hessian, refined momentum gradients, isotropic trace and
shell-integrated caloric dissipation, the empty-reservoir limit, molecular
tensor images and wedge boundary conditions, and recorded-source coherence.
The existing mixed compact/stress Welch observer admits molecular propagation
only when both owners retain the identical original relaxation state and
wavenumber, in addition to matching equilibrium medium and boundary geometry.
It retains cross-mechanism interference in free space and ideal image domains;
different relaxation rates cannot share a combined observation merely because
equilibrium density and sound speed agree.

## Native D3Q19 resolved-volume source acquisition

`curle3d::quadrupole::lbm::PreparedLbmReynoldsVolume3` observes the original
`BoundaryGrid3` through the original `MaterialScaling` density, dx and dt.
The native fallible `macroscopic` accessor returns density and Guo half-force
velocity in one population traversal; existing `velocity` retains the same
accepted-state arithmetic. Each full fluid cell contributes instantaneous
`rho u_i u_j dx^3` in xx, yy, zz, xy, xz, yz order. Density is the actual local
lattice density times its material reference, not an assumed constant one.

The caller authors a stationary proper frame, physical center of cell zero,
and half-open source box. Positions use the original on-site indices times dx
without an inferred half-cell shift. Source order is x-fastest within y/z,
skipping only prepared voxel solids. Every sample checks dimensions and all
selected-cell occupancy; it does not claim object identity or topology outside
that box. The full-cell quadrature has no fill fractions, cut-cell volumes,
periodic images, moving-frame correction or automatic source-extent completion.
No mean velocity, SGS stress, viscous stress or thermodynamic correction is
inferred. Spatial source and observation refinement remain necessary.

Explicit cell/work and backing-buffer byte budgets precede allocation. Sampling
reuses accepted/trial moment arrays and publishes all components plus a strictly
increasing physical timestamp together. Data, topology or cancellation refusal
retains the accepted sample and grid. `record_sample` forwards only accepted
moments/time into the original geometry-bound joint Welch owner and requires
exact equality of the original lattice and recording timesteps. Recorder
refusal leaves that source sample available for retry and advances no record.
The grid exposes no physical clock: the caller still binds accepted solver
states to their real times; the recorder owns uniform cadence admission.

Five G0/G1/G4 fixtures cover SI/frame and half-force tensor values, a real native
force step, fixed fluid selection, memory/work/frame admission, late malformed
cell/topology and cancellation refusal, and native endpoint-forced sinusoidal
flow through source acquisition, joint recording and the coherent observer.
The end-to-end oracle checks the actual quadratic stress's doubled frequency
and distributed coherent pressure/vector-velocity power. These fixtures add no
measured-SPL, developed-turbulence or unresolved-noise qualification.

## Resolved compressible volume-stress sources

With `compressible-volume`, `curle3d::quadrupole::flux::PreparedCompressibleVolume3`
extracts instantaneous sources from accepted `fs-flux` `Euler3d`,
`MixtureEuler3d`, or their original spatial `MolecularEuler3d` owners. One
prepared layout retains the complete original EOS, including the mixture
species axis and caloric references, the native Cartesian grid/boundaries, and
a nonempty increasing selection of cell indices. Native cell centers and
volumes are in SI; an explicit stationary proper rotation and origin map them
to the propagation frame. Velocity is absolute in that frame; no mean flow,
reference temperature or alternate gamma is inferred.

The Euler constructors retain the resolved inviscid Lighthill moments
`Q_ij = [rho u_i u_j + ((p-p0)-c0^2(rho-rho0)) delta_ij] dV`.
The molecular constructors `from_molecular_frozen` and
`from_molecular_mixture` additionally subtract `(R tau R^T) dV`, using the
original configured native stress, including authored bulk viscosity,
source-dependent shear viscosity and any explicitly selected SGS closure.
The six entries are xx, yy, zz, xy, xz, yz, without an engineering-shear factor.
The retained `FluidAcousticState` supplies rho0/c0 and the caller supplies the
positive absolute p0 of that same ambient medium. That medium can differ from
the source gas; p0 is not guessed from its composition. Every sample decodes
pressure through its original source owner. Double-double arithmetic retains
the full Reynolds, pressure-density and rotated viscous contributions together
until the final molecular moment is rounded; an exactly zero stress preserves
the original Euler arithmetic. The trace, signed isotropic
stress and physical cross stresses are preserved. The full tensor convention
follows NASA report 19770026218, equations (8)-(9):
<https://ntrs.nasa.gov/api/citations/19770026218/downloads/19770026218.pdf>.

Molecular sampling borrows one original `StressField3d` evaluation at the
accepted endpoint. The native operator owns all caloric inversions, actual
transport data, velocity gradients, wall/reservoir conditions and configured
closure contributions. Opposing conservative face tractions are averaged to
cell centers and off-diagonal pairs are symmetrized. This reconstruction is
exact for affine stress and second order for smooth resolved interiors;
material jumps, walls and modeled SGS retain their native limitations.
No acoustic-side constitutive law or duplicate gradient field is introduced.
This is an instantaneous stress, not the preceding time step's integrated
load. Original accepted flow, clock, interval reports, boundary loads and
Lagrangian closure history remain unchanged; reusable constitutive caches may
change. Mixture observations pass the supplied original dissipation controls.

All constructors publish the initial accepted source sample. Later extraction
requires the same full source/grid binding and a strictly later native accepted
time. The retained `includes_viscous_stress` selection prevents using an Euler
sampling method on a molecular sampler, or a molecular method on an Euler
sampler. All selected moments are staged before a final cancellation check and
buffer swap; refusal changes neither the accepted source sample nor its time.
There is no caller-supplied replacement clock. `record_sample` forwards
accepted positions, moments and time to the original `QuadrupoleWelchStream`;
that owner refuses a nonuniform cadence rather than resampling. Source
positions/order remain bound through propagation.

Prepared numerical payload comprises indices, positions and two six-moment
buffers (nominally 128 bytes/cell on 64-bit hosts); actual vector capacities are
admitted against the caller's bound. Original immutable EOS metadata and
existing native transport storage are separate. Sampling allocates no adapter
field buffers. Source-owned constitutive and cancellation work remains in the
native operator. Existing kernel and owner regressions cover inviscid
cancellation, original source/clock/rotation, and an accepted Euler history
through coherent recording against a direct DFT. Added molecular owner cases
specify independent full affine shear/bulk stress, positive-axis face signs,
rigid-motion/dilatation/zero-law limits, smooth periodic shear convergence,
full tensor rotation, evolving species and final-cancellation atomicity.

Euler entry points remain explicitly inviscid. Molecular entry points include
the original configured resolved stress, without inferring missing SGS noise,
applied body-force dipoles or missing source-volume contributions. Products of
cell averages do not reconstruct unresolved products. Periodic fields supply
only the selected physical box, with no implicit image sum. Cell quadrature
and source-region truncation require their own refinement; wall scattering,
noise calibration and measured-SPL qualification remain absent.

## Coherent recorded force and volume-stress sound

`PreparedCompactSourceObserver::evaluate_welch_bin_with_quadrupoles` propagates
one joint record containing 4N compact channels [Q, Fx, Fy, Fz] followed by 6M
volume-stress channels [xx, yy, zz, xy, xz, yz]. Both source groups enter the
same pressure and three-component velocity sums inside each original Welch
factor, before the existing second-moment reduction. Recorded cross-mechanism
interference is retained; separate mechanism powers are not substituted for
their coherent sum. The original peak-equivalent bin normalization applies once.

Both prepared observers must have identical acoustic state, receiver and exact
frequency. Homogeneous stationary lossless or molecularly relaxing free space
and matching ideal perpendicular-plane or pi/m-wedge domains are admitted.
Molecular observers must retain the same original relaxation state and complex
wavenumber; matching equilibrium density and sound speed alone is insufficient.
Original plane order, geometry, boundary conditions and wedge order must match
exactly; mismatches refuse.
The caller must provide disjoint physical mechanisms; a volume source already
represented by a permeable acoustic surface must not be counted twice.

Application reuses both prepared transfer banks and the shared quadrupole factor
accumulator without allocating source vectors, another spectrum or a dense
covariance. Checked combined channel and (N+M)*segment work bounds precede the
bounded-cancellation reduction. Three G1/G3/G4 fixtures cover independently
specified mixed-source field moments, destructive/constructive/decorrelated
cross-mechanism sound including surviving pressure-null velocity, and
medium/receiver/frequency/path/shape/work/cancellation refusals. Window leakage,
sampling and source-volume limits remain; this is bin-integrated statistics,
not deterministic cross-frequency phase, inferred independence or calibrated
turbulence sound.

## Full-wave observations and physical sources inside fluid layers

`curle3d::refraction::full_wave::layer::interior::PreparedLayerInteriorFullWaveMonopole`
observes an original exterior volume source strictly inside the finite layer.
The original pressure/normal-velocity transfer equations retain both internal
traveling directions, every coherent return, evanescence and reactive velocity.
The existing critical branch maps, bounded quadrature and denominator-resolution
criterion remain shared with exterior transmission. Scaled cosine/sinc
numerators cancel their growing factors before physical attenuation is applied;
near-face fields do not disappear when transmission across the whole layer
underflows. Independent positive all-evanescent envelopes use the actual
exterior clearance plus interior depth, with separate pressure and velocity
bounds. An identical layer/receiving fluid delegates the original interface.

`interior::source::PreparedInteriorSourceFullWaveMonopole` supplies the reciprocal
physical volume source inside the layer, observed in its incident exterior.
Pressure reciprocity follows the original variable-density acoustic operator;
exterior velocity uses the exterior-coordinate derivatives and actual receiving
density. It is not copied from the interior receiver field. Its tail bounds
include those exterior derivatives. Explicitly reversing the same layer's
faces and exterior states observes the other half-space. Both APIs cache four
complex transfer coefficients, apply physical Q without allocations, and feed
native harmonic observations on the original source clock without another delay.

Original stationary, lossless, infinite-plane and nontrapping-layer admission
remains. There is no source loading, absorption, flow, finite edge or guided-pole
solver. Quadrature estimates and analytic
spectral-tail bounds are not interval certificates or experimental accuracy.
Focused G1/G3 tests cover original spherical limits, independent two-family
spherical images, both-face pressure/normal-velocity continuity, receiving-fluid
momentum, reciprocal pressure, rotation, near-face evanescent survival, native
clock identity and transactional domain/work/cancellation refusals.

### Exterior physical forces observed inside the layer

`interior::PreparedInteriorFullWaveForce` applies prescribed original forces ON
the incident fluid, in newtons, to strict interior observations. It shares the
existing monopole's outgoing Green spectrum, all internal returns, scaled
pressure and normal-velocity numerators, critical branches and resonance guards.
The original `+i/(omega rho_incident) grad_source` derivative supplies both
pressure coefficients and all five cylindrical velocity-tensor coefficients.
The mixed normal/tangential entries retain different standing-wave factors;
they need not coincide or have opposite signs. Transverse and reactive motion,
including velocity at a pressure null, remain present in the Cartesian field.

All seven channels use the same adaptive samples. Separate derivative-aware
all-evanescent bounds at clearance `a+s` retain the actual incident and layer
densities and bound the pressure-vector and velocity-tensor norms per newton.
The original work/depth controls and cancellation remain authoritative. A
cached 4-by-3 transfer applies arbitrary complex forces with fixed double-double
accumulation, without allocation, quadrature or another propagation delay.
An invisible second face delegates the original single-interface force owner.
Original face distances remain double-double through oscillatory phase and
attenuation in the shared kernel, for both interior receivers and reciprocal
interior volume sources. The existing large-evanescence half-space limit is
preserved. Forces originating inside the layer are outside this force API.

`tests/full_wave_layer_interior_force.rs` covers all seven coefficients against
independent 35/50-digit direct-q references, a rotated two-family image-force
series with the original density gain and reflected normal force, and all
three original source derivatives. Exact limits, axial transverse motion,
native clocks, finite forces, strict domains, work exhaustion and cancellation
are also covered. A private phase regression retains a second-face distance
smaller than one ulp of its leading component. Numerical diagnostics remain
estimates plus analytic spectral-tail bounds, not interval certificates.

### Physical forces originating inside the layer

`interior::source::PreparedInteriorSourceFullWaveForce` radiates a prescribed
original force ON the layer fluid into the incident exterior. The source lies
strictly between the original faces and the observer strictly before the first
face. Explicitly reversing the same original layer observes the other exterior;
the force remains a Cartesian vector in its actual source fluid. The reciprocal
interior Green kernel supplies every coherent return, critical branch, scaled
numerator and resonance check without another boundary solve.

For source depth s, the original pressure numerator satisfies
`dA/ds = i rho_layer B`. Thus the normal source-force pressure uses `-WB`,
while tangential source-force pressure uses `i q WA/rho_layer`, both divided
by `2 pi`. Receiving derivatives use the actual exterior density and outward
propagation sign. The physical receiver radial basis is opposite the swapped
interior kernel's tangent. Both mixed velocity coefficients keep their distinct
standing-wave factors. The Cartesian force/velocity transfer is the transpose
of the swapped exterior-force/interior-receiver transfer; the pressure-force
vector is the negative swapped monopole velocity. Pressure reciprocity alone
does not justify copying either velocity vector or a force's source density.

Seven channels share adaptive samples and derivative-aware all-evanescent tails
at the actual exterior clearance plus source depth. The pressure bound retains
`max(rho_exterior/rho_layer,1)` and the velocity bound retains both fluid
densities. Numerical limits and cancellation remain explicit. Preparation
caches a 4-by-3 Cartesian response; force evaluation and native harmonic output
use fixed DD accumulation without allocation, quadrature or an extra delay.
An invisible second face delegates the original reversed single-interface
force owner, including its exact equal-speed and homogeneous limits.

`tests/full_wave_layer_interior_source_force.rs` uses the independently refined
35/50-digit reciprocal spectra for all seven channels, a rotated two-family
force-image series into both original exteriors, Cartesian reciprocity,
original source derivatives and exterior-fluid momentum. Exact reductions,
axial pressure-null motion, native clocks, strict domains, finite forces, work
exhaustion and cancellation are covered. Interior receivers, source loading,
guided poles, absorption, flow and finite lateral edges remain outside this
owner; quadrature estimates and tail bounds are not interval certificates.

### Coherent exterior source arrays observed inside the layer

`interior::array::PreparedInteriorFullWaveCompactSourceObserver` composes the
native interior monopole and force owners for original incident-exterior source
positions. It retains the original layer, strict interior receiver and physical
source clock. One cached 4-by-4 Cartesian transfer per source maps the original
Q [m³/s] and three force [N] channels into actual pressure and vector velocity.
All forward/backward propagation and internal returns stay coherent within
each original source. Applying source amplitudes uses the existing fixed-order
double-double array reduction, without allocation or additional quadrature.

The same transfer bank supports coherent fields and native harmonic lines,
explicit independent-source statistics, and caller-supplied low-rank coherence
factors. Phase is retained across Q/F channels, sources and internal paths
before pressure/velocity moments, active intensity and signed reactive intensity
are formed. Independence is never inferred from an unknown phase. These are
peak-phasor line powers, not spectral densities or an added propagation clock.

Per-source Q/F tolerances retain their own physical units, while one shared
evaluation cap bounds all preparations. Source and observation counts, strict
face clearances, numerical limits and cancellation refuse partial publication;
empty arrays enforce the same interior receiver domain. An invisible second
face with equal original exterior speeds uses the exact zero-quadrature path.
Equal speeds alone do not erase unequal-density internal reflections. Native
source diagnostics remain available; relative per-source tolerances do not
guarantee relative array accuracy at an interference null. Interior sources,
face traces, zero-thickness interiors and guided poles remain outside this API.

`tests/full_wave_layer_interior_array.rs` checks exact original reductions,
actual layer Q/F composition and all pressure/velocity moments, coherent versus
independent source powers, shared work exhaustion, equal-speed density-layer
work, original clocks, empty arrays and domain/cancellation refusal.

### Coherent original source arrays inside the layer

`interior::source::array::PreparedInteriorSourceFullWaveCompactSourceObserver`
propagates prescribed Q/F sources strictly inside the original fluid layer to
an observer strictly in its incident exterior. Reversing that same physical
layer explicitly observes the other exterior. Native interior-source monopole
and force owners supply one cached 4-by-4 Cartesian transfer per original
location. Q remains peak volume flow [m³/s]; force remains peak force ON the
actual layer fluid [N]. The receiver keeps its original exterior density,
pressure/vector-velocity units and physical harmonic clock.

The shared borrowed array evaluator supplies coherent fields, native harmonic
lines, independent-source statistics and explicit low-rank coherence factors.
All internal histories and Q/F phase stay within each original physical source
before statistical reduction. Pressure/velocity moments and signed active and
reactive flux retain source interference, including cancellation at pressure
nulls. Applying a prepared transfer uses fixed-order DD accumulation, without
allocation or further spectral integration; factor evaluation remains linear
in the number of sources times coherence factors.

One aggregate work cap bounds the native Q and force preparations. Per-source
diagnostics retain their different physical units and do not imply relative
array accuracy at interference nulls. Empty arrays still require positive layer
thickness, a strict original DD exterior observer, positive frequency and valid
numerical controls. Equal speeds use zero quadrature only when the second face
is also physically absent; two unequal-density interfaces retain their returns.
Source/observer geometry, source counts, work and cancellation admission publish
no partial array. Interior receivers, zero-thickness source layers, face sources,
guided poles, absorption, moving fluids and source loading remain outside scope.

`tests/full_wave_layer_interior_source_array.rs` covers original homogeneous and
reversed-interface limits, actual finite-layer Q/F composition, all coherent and
partial-coherence moments, aggregate work exhaustion, source nulls, clocks,
strict DD domains, empty arrays and cancellation.

## Causal transient Lighthill volume-stress observation

`curle3d::quadrupole::transient::PreparedTransientQuadrupoleObserver` observes
stationary symmetric stress moments in the original homogeneous lossless
`FluidAcousticState`. Pressure retains the complete retarded stress,
first-derivative and second-derivative terms. Vector particle velocity also
retains the causal stress integral and transverse induction field. Tensor
trace and physical off-diagonal entries are preserved. Source pressure and
velocity are summed coherently before forming the signed instantaneous flux
`p u`; neither a harmonic half-product nor a plane-wave intensity shortcut
is used.

Preparation caches each source's geometry, material coefficients and travel
time under an explicit source-count budget. Evaluation allocates nothing,
uses fixed-order double-double accumulation and checks cancellation. The
receiver-time callback observes each source at its own emission time, applying
propagation delay exactly once. Caller-owned source histories retain their
original derivatives, integral, interpolation and bandwidth; callback/block
boundaries do not reset physical history. Invalid inputs, unresolved timestamps
and nonfinite outputs refuse without publishing a partial field.

`evaluate_at_clock_dd` and its cancellable counterpart retain the receiver and
emission timestamps in double-double seconds through the source-history
callback. They accept the existing `HarmonicSampleClock::time_at` clock and
preserve short source-dependent delays under a common large epoch translation.
Callers retain DD timestamps until subtracting their own history epoch or
evaluating phase; no phase origin, source alignment or interpolation is inferred.
The compatible `evaluate_at_time` API rounds callback timestamps to f64: even
when a rounded delay is positive, its phase accuracy can be insufficient for
destructive interference. G1/G3 DD-clock regressions exercise full tensor fields
and an analytic 16 kHz two-source pressure null at a 2^30-second epoch. G4 checks
normalization, nonfinite/overflow/unresolved times and bounded cancellation.

The same transient observer admits `with_orthogonal_planes` and `with_wedge`
using the original ideal-wall geometry and finite image groups. Every image
keeps its independent retarded time and complete tensor/vector transform;
the harmonic folded coefficient is never treated as a single transient delay.
`source_count`, `source_position_m` and `source_delay_s` retain physical-source
and direct-path meanings. `path_count`, `path_source_index`, `path_delay_s` and
`source_delay_bounds_s` expose all image paths for complete history-window
admission. The callback reads the same original physical source at each path's
emission time; explicit retarded input contains one history per path. Images
apply `s H Q H^T` equivalently through virtual receiver `H^-1 x`, pulling
velocity back by `(H^-1)^T`, including nonsymmetric wedge rotations. Preparation
and storage are O(paths) under an explicit source-path cap; repeated observation
is allocation-free O(paths) with bounded cancellation and no transcendental work.
The complete pressure/velocity sum precedes instantaneous flux, retaining
interference between direct and reflected arrivals. No image is silently cut
off, and callback boundaries never advance or reset the physical source history.

G1/G3 regressions compare harmonic field/mean-flux parity, independently
differentiated retarded Green potentials, all three momentum equations, causal
pulse arrivals, transverse velocity at pressure nulls and coherent interference.
G1/G4 `transient_quadrupole_reflection` controls independently transform oblique
source images, distinguish direct/echo onset, enforce all six tensor components
at mixed rigid/release walls, check wedge harmonic and momentum parity, and
exercise source identity, complete delay bounds, budgets and cancellation.
These ideal infinite walls/wedges add no source-model validation, finite-cell
quadrature bound, moving quadrupole, absorption, finite wall impedance, general
diffraction, closed-room reverberation, reciprocal loading or measured-SPL claim.

## Sources and receivers inside the original refracting layer

`curle3d::refraction::full_wave::layer::interior::internal::PreparedLayerInternalFullWaveMonopole`
prepares complex pressure and all three particle-velocity components with both
the original volume-flow source and receiver strictly inside a finite fluid
layer. The direct spherical field uses the original layer-fluid owner; the
returned field retains all four coherent image families and the original
physical source clock. Reflected-only evaluation is finite at a coincident
source and receiver; total-field evaluation explicitly refuses that singularity.

The two returning spectral families reuse the native layer denominator, branch
maps, bounded quadrature and scaled interior functions. Each uses its actual
positive source-plus-receiver return clearance for the evanescent cutoff.
Separate cutoffs avoid forcing long-distance returns onto a near-face spectral
tail. Double-double depths and rationalized admittance differences preserve
original geometry and weak material contrast. Absolute error allocations and
tail bounds add across families even at destructive interference; preparation
shares one evaluation cap and cancellation contract. Warm evaluation applies
cached coefficients without allocation or further quadrature.

An absent second face delegates to the original reversed single-interface
owner, including the exact homogeneous limit. Native harmonic-line output
retains pressure, vector velocity, and signed active/reactive intensity.
The scope is stationary, lossless, infinite planar, nontrapping layers and
prescribed monopoles. Face sources, guided poles, absorption, moving media,
source loading and a finite total self-field are not supplied.

`tests/full_wave_layer_internal.rs` adds independent spherical-image references,
deep evanescent and equal-depth cases, reciprocity, all receiving momentum
components, continuity-equation divergence, both-face pressure/velocity
continuity, rotated geometry, finite reflected self-fields, symmetry,
aggregate work exhaustion, cancellation and invalid-input refusals.

## Physical force sources and receivers inside the original fluid layer

`interior::internal::PreparedLayerInternalFullWaveForce` extends the original
internal-layer Green function to prescribed three-dimensional forces ON the
layer fluid in newtons. It caches boundary-return pressure and all Cartesian
particle-velocity components, plus the coherent sum with the native direct
force field when source and receiver are distinct. Reflected self-fields remain
finite; total self-fields explicitly refuse. The original exp(-i omega t)
source clock and harmonic-line intensity conventions remain unchanged.

Source derivatives use `+i grad_source/(omega*rho_layer)`. Both returning
families reuse the monopole owner's scaled pressure and normal-velocity
numerators, original layer denominator, critical maps and DD clearances.
Receiving momentum retains all seven cylindrical force channels, including the
two distinct mixed velocity entries and their source/receiver normal signs.
Each sample computes the common spectral basis once; warm evaluation applies
cached Cartesian coefficients without quadrature or allocation.

The two active families divide the absolute tolerances, use zero relative
tolerance and share one evaluation budget. Per-face analytic tails retain the
actual positive source-plus-receiver return length and the additional source
and receiver derivative powers. Summed diagnostics bound omitted returns and
estimate reflected-transfer quadrature; they do not certify relative accuracy
at a coherent total-field null. An absent second face delegates to the native
reversed single-interface force owner, including its exact equal-speed and
homogeneous limits. Cancellation and exhaustion publish no partial result.

`tests/full_wave_layer_internal_force.rs` covers all seven force channels against
independent direct-q reference values refined at 35 and 50 decimal digits,
rotated four-family spherical-force images, original source derivatives,
receiving momentum, velocity reciprocity, coherent direct/total separation,
self-fields, invisible-face limits and domain/work/cancellation refusals.
The existing scalar implementation and its regressions remain in place.
The scope remains prescribed sources in stationary, lossless, infinite planar,
nontrapping fluid layers; guided poles, absorption, source loading and finite
edges are not supplied. Native Rust execution is not claimed by this extension.

## Coherent compact-source arrays inside the original fluid layer

`interior::internal::array::PreparedLayerInternalFullWaveCompactSourceObserver`
prepares original volume-flow and three-dimensional force sources, with every
source and the receiver strictly inside the same finite layer. `new` retains
the coherent direct-plus-returning field; `new_reflected` retains only the
boundary returns and admits the finite reflected self-field. Total-field
preparation rejects a coincident source before its spectral work. Both
constructors validate the original DD receiver domain even for an empty array.

The adapter caches one Cartesian 4-by-4 transfer per physical source from the
native internal monopole and force owners. It reuses the common fixed-order DD
array evaluation, harmonic clock, explicit low-rank coherence factors and
independent-source statistics. Q/F amplitudes and every direct/reflected path
remain coherent within each original source. Pressure, all velocity components,
and signed active/reactive intensity use the actual field; no extra propagation
delay or plane-wave intensity approximation is introduced. Warm application
allocates nothing and repeats no quadrature.

The original per-mechanism controls and one aggregate evaluation cap cover all
sources, Q/F preparations and active returning families. Original unit-Q and
unit-force diagnostics remain separate. Exact absent-second-face/equal-speed
reductions consume no spectral evaluations; equal sound speeds alone do not
remove unequal-density interfaces. Cancellation or exhaustion publishes no
partial array. Per-source estimates give no relative array-error guarantee at
interference nulls, and independent statistics require caller-supplied source
independence. Coherence factors describe line powers, not automatically PSD/Hz.

`tests/full_wave_layer_internal_array.rs` covers direct/return Q/F composition,
coherent and factorized statistics, source cancellation, original clocks,
exact reductions, strict domains, finite reflected self-fields and aggregate
preparation/coherence budgets. The same stationary, lossless, infinite planar,
nontrapping-layer and prescribed-source scope applies; source loading, guided
poles, absorption, moving fluids and finite edges remain outside this API.

## Full-wave refraction of original volume-stress sources

`curle3d::refraction::full_wave::quadrupole::PreparedFullWaveQuadrupole`
transmits the existing symmetric `QuadrupolePhasor` [N m] across an original
planar interface between two stationary homogeneous lossless fluids. Analytic
source Hessians and receiving derivatives retain all six physical stress
entries, their trace and shear, full complex pressure and vector velocity.
No receiving plane-wave impedance or refracted-direction approximation replaces
the near field. Source and receiver normal derivatives retain their distinct
original fluid wavenumbers.

Fourteen normalized cylindrical channels reuse the native spectral branch maps,
bounded adaptive integration and cancellation. Square-root-two shear weights
make pressure/velocity diagnostics complete Frobenius transfer norms. A stable
axis continuation and logarithmically balanced products preserve physical
derivatives near zero lateral range and after raw evanescent decay underflows.
The shared positive tail envelope now admits cubic derivative powers explicitly.
Exactly equal original speeds use the native homogeneous quadrupole owner,
including the correct pressure and velocity factors for a density jump.

Preparation caches a Cartesian 4-by-6 transfer; warm field and harmonic-line
evaluation allocates nothing and performs no quadrature or additional delay.
Each cell remains a prescribed compact stress moment. Source-size error,
unresolved SGS stresses, mean flow, absorption, finite boundaries, loading and
measured-SPL qualification are outside scope. Quadrature estimates and analytic
spectral tails are separate from geometry/Bessel/floating-point errors.

Six native tests cover homogeneous/density-jump limits, independent derivatives
of force sources, trace and axis shear, receiving momentum/divergence, full
tensor rotation, density scaling, normalized errors, admission/cancellation and
third derivatives after raw-decay underflow. A shared-tail test independently
checks the cubic envelope by integration by parts.

## Refracted volume-cell arrays and original recorded stress spectra

`full_wave::quadrupole::array::PreparedFullWaveQuadrupoleObserver` prepares
the native transmitted 4-by-6 stress transfer at every original cell position.
One source-count and aggregate quadrature budget bounds the complete array;
per-source tensor diagnostics remain in their physical units. Coherent fields,
explicitly independent cells and low-rank coherence factors preserve all six
physical stress entries and the caller's correlation assumptions. Pressure and
vector velocity are reduced before signed intensity or second moments.

`RecordedQuadrupoleSpectrum::prepare_refracted_observer` binds the existing
record's geometry and frequency to these transfers.
`evaluate_refracted_bin` refuses reordered or relocated cells before applying
their original joint Welch factors. The existing peak-bin owner applies the
one-sided sqrt(2 df) normalization exactly once. No extra density, cell volume,
phase origin, ray delay or independent-cell approximation is introduced.

Warm observation uses fixed-order double-double reductions without allocation,
FFT, quadrature or dense source covariance. Work is linear in cells times
recorded factors. Cancellation and invalid source/frequency/work admission
publish no partial field or statistics. These outputs are bin-integrated
moments, not a reconstructed waveform or a physical-validation claim for the
underlying flow/source model. Source-window/leakage/spatial-quadrature and
single-interface propagation limits remain explicit.

Four G1/G3/G4 tests compare native tensor sums and flux, coherent cancellation
against independent/partial-coherence statistics, recorded Reynolds-stress
transmission against the independent equal-speed pressure/velocity factors,
and geometry/frequency/aggregate-budget/cancellation refusals.

## Reflected and coherent total original quadrupole fields

`PreparedReflectedFullWaveQuadrupole` extends the original stress source to
incident-side observations at one infinite planar fluid interface. Source and
receiver remain strictly in the original incident fluid. The source is the
original six physical symmetric stress moments in N m on its existing
`exp(-i omega t)` clock. Reflected-only and coherent direct-plus-reflected
queries expose complex pressure, all three particle-velocity components,
signed active/reactive intensity and the original harmonic-field line type.

The original reflected Green kernel supplies critical-circle branch maps,
grazing cancellation and logarithmically scaled spectral weights. Analytic
source Hessians and receiving derivatives retain physical shear, tensor trace,
near-field, reactive and evanescent contributions. The returning normal
velocity has the opposite propagation sign to transmission. The published
transmission owner's fourteen-channel tensor reconstruction and ten velocity
diagnostic entries are shared without changing its public API or its coherent
array/recorded-source consumers. Pressure and velocity retain their respective
degree-two and degree-three analytic omitted-tail envelopes.

Exactly equal original sound speeds use the existing homogeneous quadrupole
observer with image stress `S M S`, where `S = I - 2 n n^T`, and reflection
factor `(rho2-rho1)/(rho2+rho1)`. The image displacement is formed with the
original double-double geometric separation. Matched original fluids have
exactly zero reflected transfer. The direct native contribution and reflected
transfer are cached separately; total pressure and velocity are added before
intensity is formed, retaining their interference cross terms. Warm source
application uses fixed 4-by-6 reductions without allocation, quadrature or
transcendental evaluation.

A coincident prescribed source and receiver can query the finite reflected
field; a total-field query refuses its singular direct contribution. Native
cancellation, material domains and bounded work remain explicit. Numerical
relative tolerances control the reflected tensor transfer, not a coherent
total-field cancellation; its absolute estimates carry to the total field.
The owner does not add source loading, unresolved stress generation,
absorption, mean flow, finite edges, multilayer reflection or an interval
certificate.

Focused regressions specify exact rotated image tensors and density limits,
independent full-wave reference coefficients, native source/receiver
derivative relations, interface pressure/normal-velocity continuity, coherent
energy interference and self-field/domain/work refusals. Independent direct
transverse-wavenumber integrations with J2/J3 angular moments were executed at
35 and 50 decimal digits for both sound-speed orders; all fourteen complex
coefficients agreed at all 30 printed significant digits. Pinned
formatter/parser checks passed; the native Rust regressions were added but
were not executed in the shared editing environment.


## Full volume-stress transmission through a coherent finite fluid layer

`full_wave::layer::PreparedLayeredFullWaveQuadrupole` propagates the original
six physical symmetric stress moments [N m] from the incident fluid to a
receiver beyond the layer. Analytic source Hessians and receiving derivatives
retain trace, both shear contractions, complex pressure and every velocity
component, including transverse, near, critical and evanescent contributions.
Source-normal and receiving-normal derivatives use their respective original
fluid wavenumbers; pressure divided by a plane-wave impedance is not substituted.

Fourteen normalized cylindrical channels share the existing resonant layer
Green kernel, branch regularization, logarithmic scaling, cancellation and
bounded adaptive quadrature. The original layer denominator retains all coherent
internal returns. Quadratic pressure and cubic velocity tail envelopes use the
physical source/layer/receiver clearance and original density factors. Diagnostics
retain separate quadrature estimates and analytic omitted-spectrum bounds in
Frobenius norms per unit stress; neither is an interval or total floating-point
certificate. An absent layer or three identical fluids delegate to the original
single-interface owner with its exact reductions.

Preparation caches one Cartesian 4-by-6 transfer. Warm field and harmonic-line
observations perform fixed-order DD reductions without allocation, quadrature,
transcendental calls or another propagation delay. Signed active/reactive
intensity follows coherent pressure/velocity contraction. The scope remains
prescribed compact sources in stationary, lossless, infinite planar nontrapping
fluids. Source loading, unresolved SGS stresses, finite edges, absorption, moving
media and guided-pole regimes are not supplied.

Six integration regressions specify collapsed/identical reductions, independent
convergent spherical images for every stress component, native force-source
finite-difference references in both exterior speed orders, Helmholtz trace and
on-axis shear, receiving momentum/mass continuity, density scaling, physical
error norms and domain/work/cancellation refusals. One unit regression checks
the cubic tail by an independent integration-by-parts recurrence. Pinned Rust
formatting/parser checks passed; the native regressions were authored but not
executed in the shared editing environment.


## Incident-side stress arrays and their original joint recordings

`PreparedFullWaveQuadrupoleObserver::new_reflected` and `new_total` extend the
existing array and recorded-spectrum path to the original incident fluid.
`FullWaveQuadrupoleField` identifies transmitted, returned-only or coherent
incident-total transfers explicitly. Every cell reuses the native reflected
quadrupole owner. Preparation caches its complete physical 4-by-6 transfer;
all original shear, image orientation, near-field and evanescent conventions
remain with that owner. Warm coherent, independent-cell, low-rank and Welch
observations use the same allocation-free fixed-order DD reductions.

Direct and returning pressure/velocity interfere before intensity or statistics
are formed, including within each original recorded factor. The recording adds
`prepare_reflected_observer` and `prepare_total_observer`; the existing
`evaluate_refracted_bin` consumes any admitted field kind while retaining
original source positions/order, clock, six-channel stress normalization and
one-sided bin factors. No extra density, cell volume, phase delay, FFT or dense
source covariance is introduced.

Incident-side receivers are validated even for empty arrays. Returned-only
observations admit finite reflected self-fields; total observations refuse
coincident prescribed sources. One aggregate quadrature budget covers all cells,
with original per-source diagnostics and transactional cancellation. For total
fields, reflected-transfer absolute errors carry through the exact direct sum;
relative accuracy at interference nulls is not asserted. The existing stationary,
lossless, infinite-plane, prescribed-source limitations remain in force.

Three focused G1/G3/G4 regression groups specify native returned tensors and
coherent direct/return intensity, an independent direct-plus-image four-cell
recording against original two-cell total observation, matched-fluid silence,
geometry/frequency preservation, finite self-fields, strict fluid sides and
aggregate-work/cancellation refusals. Pinned Rust formatting/parser and bounded
source/API review passed; native Rust test execution remains unperformed in the
shared editing environment.

## Joint recorded compact and volume-stress refraction

`PreparedFullWaveCompactSourceObserver::evaluate_welch_bin` now propagates the
original synchronized source-major `[Q,Fx,Fy,Fz]` channels through its cached
pressure and full vector-velocity transfers. The companion
`evaluate_welch_bin_with_quadrupoles` consumes one joint original record with
those compact channels followed by the prepared physical stress blocks
`[xx,yy,zz,xy,xz,yz]`. Volume flow is m³/s, force on fluid is N, and integrated
stress moment is N m; shear entries are not pre-doubled.

Both banks must share original fluid states, oriented interface, receiver,
exact frequency and transmitted, reflected-only or incident-total field kind.
Each original Welch factor enters one common DD pressure/velocity reduction
before second moments are formed. This retains Q/F/stress, cell/component and
direct/reflected interference, including velocity at pressure cancellation.
The native peak-bin owner applies its normalization once; outputs are
bin-integrated moments, not PSD/Hz or deterministic cross-frequency phases.

Application is allocation-free O(total sources × original segments), reuses
prepared transfers and checks the combined source/factor budget before reduction.
No FFT, quadrature, source covariance or independent source spectrum is rebuilt.
Original recording/channel geometry, disjoint-source responsibility,
finite-record limitations and the stationary lossless infinite-plane model
remain. Unknown correlation is not replaced by independence.

Three focused G1/G3/G4 regression groups specify known synchronized phasors
for all three fields and both equal/unequal sound speeds; analytic compact/stress
pressure cancellation with surviving velocity and measured correlation changes;
and source-domain, clock, shape, aggregate-work and cancellation refusals.
Pinned Rust formatting/parser checks and bounded source/API review passed.
Native Cargo tests were authored but not executed in the shared editing environment.

## Finite-layer reflected and incident-total compact quadrupoles

`curle3d::refraction::full_wave::layer::reflection::PreparedLayeredReflectedFullWaveQuadrupole` propagates the original six symmetric stress moments [N m] to pressure and all three particle-velocity components in the incident fluid. The native layer reflection denominator retains every coherent internal return, critical and evanescent contributions, and the original physical shear and trace conventions. Source Hessians and receiving derivatives share the original fourteen returning angular channels and logarithmically scaled products; this owner introduces no product finite differences or second wave solver.

Preparation caches reflected and, where finite, coherent direct-plus-return Cartesian transfers. Warm application uses fixed-order DD contraction before signed active/reactive intensity. A coincident source and receiver has a finite returned field but refuses the singular total point-source field. A vanished layer or absent second face delegates exactly to the native single-interface owner. Four pressure and ten velocity diagnostics retain bounded native quadrature, quadratic/cubic analytic tails, cancellation, and work limits. Total-field error estimates are absolute reflected-field estimates and give no relative guarantee at interference nulls. The existing stationary, lossless, infinite-plane, nontrapping-layer and prescribed-source applicability limits remain in force.

## External compact quadrupoles observed inside a finite layer

`full_wave::layer::interior::PreparedInteriorFullWaveQuadrupole` maps the original incident-exterior symmetric stress [N m] to pressure and vector velocity at a receiver strictly inside the finite layer. The original source Hessian uses the incident wavenumber. The native interior pressure numerator and independent normal-velocity numerator retain both internal propagation directions, including normal motion at pressure nodes; tangential momentum uses the actual layer density. This is not a traveling-wave impedance reconstruction.

The original fourteen tensor channels share critical maps, resonance resolution, scaled numerators, quadrature, cancellation and work bounds. Quadratic/cubic analytic tail envelopes use the physical source-depth-plus-receiver-depth clearance and the incident/layer density bounds. An absent rear face delegates to the native single-interface quadrupole owner. Preparation caches a 4-by-6 Cartesian map; warm evaluation uses the shared fixed-order DD stress contraction, retaining the original harmonic clock and coherent signed flux. Sources inside the layer remain outside this owner's domain, as do trapping, lossy, moving or finite-edge layers. Quadrature estimates and analytic tails are not interval certificates or measured-SPL qualification.

## Finite-layer quadrupole arrays and original recorded stresses

`full_wave::layer::{LayerQuadrupoleField, PreparedLayeredFullWaveQuadrupoleObserver}` prepares bounded arrays of original exterior stress moments for transmitted, reflected-only, coherent incident-total, or interior receivers. The explicit field kind admits the actual receiver fluid even for an empty array. Native per-cell owners retain complete tensor transfers and diagnostics; a single aggregate work cap bounds the entire preparation. Exactly reduced interfaces preserve their zero-work path when applicable. Equal sound speeds alone never erase density-layer internal returns. Reflected point-source self-fields remain finite; incident-total self-fields refuse.

Coherent phasors, independent cells, low-rank coherence factors, and joint Welch bins reuse the original stress-array DD reducers. Cell, tensor, and internal-path cross terms are retained within each coherent factor before pressure, velocity, and signed active/reactive statistics are accumulated. Warm reductions introduce no quadrature or covariance-matrix expansion. The original single-interface reducer and mixed-source `accumulate_welch_factor` entry point retain their arithmetic and normalization.

`RecordedQuadrupoleSpectrum::prepare_layered_observer` and `evaluate_layered_bin` carry original recorded geometry, source order, harmonic clock, physical N m stress units, and native Welch peak-bin conversion into all four layer views. Frequency/bin and geometry mismatches refuse. No added density, volume, phase reset, propagation delay, FFT, or inferred independence is introduced. Outputs retain bin-integrated power and the native stationary-source and layer-domain limits. Per-cell absolute estimates do not imply relative error control at an array interference null.

## Original compact quadrupoles inside a layer radiating outward

`full_wave::layer::interior::source::PreparedInteriorSourceFullWaveQuadrupole` maps an original symmetric stress moment [N m] strictly inside a finite fluid layer to pressure and vector velocity in its incident exterior. The opposite exterior uses the existing explicit reversal of the same layer's outside fluids, first-face point and normal; the source tensor stays in its original Cartesian frame and source clock.

The original reciprocal interior Green kernel supplies independent pressure/normal weights. Mixed normal/tangential source derivatives use A_s = i rho_layer B; two normal source derivatives use A_ss = -zeta_layer² A. Exterior receiver momentum uses the exterior wavenumber and density. Both internal directions, every coherent return, physical trace/shear, near and evanescent contributions survive the shared fourteen-channel Cartesian reconstruction. Native branches, denominator, resonance checks and logarithmic attenuation recovery remain in one solver.

Preparation caches a complete 4-by-6 transfer with the existing four-pressure/ten-velocity diagnostics. Quadratic/cubic analytic tail envelopes use actual exterior clearance plus source depth and original source/receiver density bounds. An absent second face delegates to the reversed native single-interface owner. Warm applications use the shared fixed-order DD contraction without allocation, quadrature, product finite differences or another delay. Strict source/receiver domains, nontrapping lossless infinite-plane scope, cancellation and finite work limits remain enforced. Interior receivers and source loading are outside this owner.

Six focused physical regression groups specify exact interface/homogeneous reductions, rotated two-family image sums in both exterior orientations, independent native-force source derivatives for both speed orders, source-layer Helmholtz trace and axis shear, exterior momentum/divergence/density scaling, tail norms and domain/work/cancellation refusals. Formatter checks passed; native numerical execution was unavailable in the editing environment.

## Original inside-layer stress arrays and recorded exterior observations

`layer::interior::source::quadrupole::array::PreparedInteriorSourceFullWaveQuadrupoleObserver` prepares complete native 4-by-6 pressure/velocity maps for bounded arrays of original inside-layer stress cells. The observer lies in the incident exterior; explicit reversal of the original layer admits the other exterior while retaining source positions, Cartesian moments and phase. Positive layer thickness, original receiver domain, frequency and numerical controls are admitted even for an empty array. Native source owners admit every cell. A single aggregate quadrature cap bounds preparation; exact source-fluid interface reductions retain zero-work evaluation when applicable, and equal speeds alone do not erase layer reflections.

Coherent fields, independent cells, supplied coherence modes and joint Welch bins reuse the same generic fixed-order DD stress-array reducers. Every tensor and internal-path cross term stays coherent within its original cell or supplied factor before statistics are formed. Warm reductions require no quadrature, alternative covariance expansion or changed shear convention. Diagnostics retain their native per-N-m pressure and velocity units and absolute-error scope at interference nulls.

`RecordedQuadrupoleSpectrum::prepare_interior_source_layered_observer` and `evaluate_interior_source_layered_bin` bind the recording's actual inside-layer source order and positions, original sample clock, stress units, window and joint spectrum to exterior propagation. Native peak-bin normalization is applied once. Changed geometry/order, wrong bin/frequency and unsupported endpoints refuse. No additional density, volume, phase origin, propagation delay or inferred independence is introduced.

Four focused regression groups specify both exterior orientations against native single-source fields and flux, coherent cancellation versus independent-source power, empty-array/domain/exact-reduction/work/cancellation behavior, and original recorded clock/geometry/bin power against the native homogeneous observer. Formatter checks passed; native numerical execution was unavailable in the editing environment.


## Moving compact force fields with explicit past forcing

`curle3d::moving::force` extends the original moving compact pressure operator
with full Cartesian particle velocity for mixed volume and force sources.
`MovingCompactObserver::prepare_force_emission` retains the original pressure
coefficients and adds constant-work boundary coefficients only for this opt-in
field owner. The original pressure and monopole preparation costs are unchanged.

The force velocity is the Cartesian gradient of the divergence of
`integral F(s)/|x-y(s)| ds / (4*pi*rho)`, with the upper limit equal to the
retarded source time. Its local boundary terms retain source velocity,
acceleration, transverse force, induction and radiation. The caller supplies
the unscaled past-worldline integral
`J = integral (3*n*(n.F)-F)/R^3 ds` [N s/m³] at the SAME receiver position;
`force_history_integrand` evaluates this kernel without a Doppler multiplier.
A zero prefix declares zero prior force contribution. Unknown initial fields
or a prefix evaluated at a different microphone position cannot be inferred
from the local force sample.

Applications preserve the original pressure law, sum force and volume velocity
before forming signed instantaneous `p*u`, and allocate nothing. Continuous
clearance over the complete retained force history, integration accuracy,
initial forcing, source bandwidth and causal root coverage remain the history
owner's obligations. Pointwise quadrature checks do not certify inter-node
clearance. Nonfinite or unrepresentable coefficients, history, fields or flux
refuse the full field without changing pressure-only admission.

G1/G3 `tests/moving_force_field.rs` checks the original stationary-force limit,
coherent mixed volume/force fields, an independent retarded vector-potential
Hessian for accelerating approach/recession/transverse motion at three
refinement steps, rotation covariance and admission refusals. This remains
prescribed linear compact-source acoustics in a homogeneous lossless medium,
without finite-source, viscous, scattering, feedback or measured-SPL claims.

### Coherent moving-force fields through ideal reflections

The fixed and moving reflected trajectory owners expose `force_field_at` for
the complete orthogonal or admitted pi/m image group. Every path retains its
own original causal source event. Full Q/F pressure and particle velocity,
including induction and the past-force term, are summed in double-double
before the signed instantaneous pressure-work flux `p*u` is formed. A pressure
node therefore does not erase particle motion, and image powers are not added
independently. Velocity follows parity times the transpose virtual-receiver
Jacobian into physical fluid-rest axes; microphone velocity is never added.

The history callback receives the path index, its actual virtual observer and
its retarded source time. It supplies the original `J` at that fixed virtual
position. A direct-path history or one evaluated at a previous receiver event
cannot be substituted. Prehistory, continuous separation over the entire past
worldline and integration accuracy retain the original force owner's explicit
obligations. Existing causal, finite-history, image-group and cancellation
admission remains authoritative; no failed path publishes a partial field.
Preparation reuses each native retarded pressure event, and field queries
allocate nothing. In admitted tangential mean flow, `p*u` remains the pressure
work flux and does not include mean-flow advection of acoustic energy.

G1/G3 `tests/moving_reflected_force.rs` compares all eight corner boundary
combinations with independent stationary Green-function derivatives and source
images; distinguishes coherent flux from summed path flux; checks complete
rotated wedge forces at a moving microphone; and verifies rigid/release wall
conditions and nonzero pressure-node motion for an accelerating source.
Missing or nonfinite image history, incomplete time coverage and cancellation
refuse the complete result.

### Coherent source arrays at moving microphones

`trajectory::PreparedMovingSourceArray` observes the complete instantaneous
Q/F field of independently moving original sources at one original microphone.
The direct constructor borrows admitted moving pairs; the reflected constructor
borrows complete orthogonal or pi/m image groups. Every pair must reference the
same original receiver object, and reflected sources must share the same image
owner. Matching coordinates at one event or equal path counts are insufficient.
The common query window is the intersection of native source windows. Checked
aggregate source, source-path and causal-root refinement budgets are admitted
before field queries; native knot searches and history-callback work retain
their original explicit bounds.

Callbacks preserve input source order, exact image identity, each path's own
retarded source clock, and its current virtual observer. Original pressure and
Cartesian velocity contributions remain in DD across BOTH source and image
sums, before one signed instantaneous `p*u`. Weak force fields therefore
survive destructive cancellation of stronger volume fields even when separately
rounded public source fields would lose them. Receiver motion supplies no
amplitude multiplier. Construction and queries allocate no owner buffers,
clone no trajectories, and neither interpolate nor advance source mechanics.
Past forcing, continuous worldline clearance and integral accuracy remain with
the history owner; a previous receiver position's force integral cannot be
substituted. This is prescribed coherent propagation, without an inferred noise
law, audio resampling or mean-flow energy-advection model.

Four focused public tests compare three accelerating sources and a moving
microphone, with direct and mixed ideal-wall paths, against independently
root-solved Green-potential derivatives under refinement; retain a 2^-60 force
signal after volume-source cancellation; enforce common receiver/boundary,
coverage and aggregate work; and refuse late missing/nonfinite histories or
cancellation with exact replay. Those four groups and the four prior reflected
full-field groups pass in the focused native source build. The run uses the
original moving/trajectory, interval, material-state and image algorithms, with
unrelated database owners and unused root error variants omitted; it is not a
full-crate or repository gate.

## Uniformly moving transient quadrupole pressure

`MovingCompactObserver::prepare_uniform_quadrupole` evaluates the complete
spatial Hessian of `Q_ij(tau)/(4*pi*R*kappa)`, `kappa = 1-n.v_source/c`,
using the original analytic retarded clock and Mach/clearance/residual limits.
Eighteen fixed coefficients retain near, induction and radiation pressure,
the physical stress trace and doubled symmetric shear entries. Q retains
physical N m units; Qdot/Qddot are its source-clock derivatives. All use fixed
medium-frame axes; translation adds no rotation or density factor. Preparation and DD contraction
allocate nothing. The optional moving-receiver clock ratio is
`(1-n.v_receiver/c)/(1-n.v_source/c)`, without a receiver pressure multiplier.

Scope remains uniform subsonic compact-source translation in a homogeneous
quiescent lossless medium. Source history and derivative consistency are
caller-owned. Pressure does not require the stress integral; moving quadrupole
particle velocity, acceleration, reflected paths and finite-cell quadrature
are not inferred. `tests/moving_quadrupole.rs` specifies stationary radial
orders, independent retarded Green-Hessian refinement, physical clock ratios,
rotation covariance and domain refusal.

## Accelerating transient quadrupole pressure

`MovingCompactObserver::prepare_emission_quadrupole` and
`PreparedSourceTrajectory::prepare_quadrupole_receiver` extend the same
18-coefficient stress-pressure evaluator to accelerating subsonic motion.
They differentiate the original retarded potential, rather than attaching a
uniform-motion amplitude to a curved delay. With `D=R*kappa`,
`C=1-|v/c|²+r.a/c²`, their symmetric spatial derivatives are
`tau_ij=[C*n_i*n_j-kappa*(n_i*M_j+M_i*n_j)-kappa²*delta_ij]/(c*R*kappa³)`
and `D_ij=-c*C*tau_ij+(n_i*a_j+a_i*n_j)/(c²*kappa)
-(r.jerk-3*v.a)*n_i*n_j/(c³*kappa²)`. These retain acceleration and jerk in
the Green amplitude as well as the original Q/Qdot/Qddot source-clock terms.
Uniform translation preserves the existing analytic preparation arithmetic.

The trajectory entry point reuses the original admitted quintic, safeguarded
retarded inverse, residual tolerance and continuous Mach/clearance bounds.
`kinematics_with_jerk_at` analytically differentiates its original acceleration
controls in fixed medium-frame axes. It adds no stored trajectory arrays or
new integration law. Original C2 knots need not have continuous jerk: internal
knots use the left polynomial, matching the root owner's endpoint convention;
the first knot uses the right polynomial. Pressure at such a knot uses this
one-sided convention and does not assert a unique classical Hessian there.

The represented polynomial, source history and their convergence remain
caller-owned. A bound on retarded source-clock rate does not certify bandwidth
of geometric modulation or of a jerk jump. This is prescribed compact stress
in the original homogeneous lossless fluid rest frame; moving quadrupole
velocity fields, finite-cell quadrature, material absorption and source
feedback are not supplied. Preparation and evaluation allocate nothing.

The same pressure owner now prepares original moving microphones through
`PreparedMovingReceiverSource::prepare_quadrupole_receiver` and every ideal
image through `PreparedReflectedSourceTrajectory::prepare_quadrupole_path`
and `PreparedReflectedMovingReceiverSource::prepare_quadrupole_path`. These
reuse the native causal roots, complete wall/wedge groups, source and receiver
frames, clearance admission and root caps. Independent source paths share
the SAME analytic quintic-jerk evaluator and one-sided C2-knot convention.

At each fixed arrival event the spatial Green-potential Hessian uses source
velocity, acceleration and jerk only. Microphone motion sets its observation
position and retarded clock rate; it adds neither a pressure multiplier nor
receiver-acceleration terms to that Hessian. Each image reads the original
six-entry Q/Qdot/Qddot at its own source time. Evaluating at the inverse-image
receiver implements the full rank-two transformation, equivalent to image
stress `s H Q H^T`, with pressure character s applied once. Trace, physical
shear entries and nonsymmetric wedge rotations are retained. The optional
clock-rate query takes physical receiver velocity in the common fluid frame;
its ray is pulled back without pressure parity.

`PreparedMovingQuadrupoleSample::pressure_dd` retains tensor/radial
contraction through later source/image sums; `pressure_pa` rounds that same
result. This is a pressure capability, not a moving quadrupole particle-velocity
or flux reconstruction. The five G1/G3/G4 cross-crate regressions in
`fs-couple/tests/moving_reflected_quadrupole.rs` exercise independently
bisected image-source spatial Hessians under accelerating source/microphone
motion, six-image rigid/release wedges, weak coherent shear residuals,
complete image support/work/Doppler admission and transactional source reads.

G1/G3 `tests/accelerating_quadrupole.rs` checks stationary/uniform recovery,
independent bisection and spatial-potential Hessian refinement with acceleration
and jerk, exact polynomial jets and C2 knot sides, constant-stress jerk
radiation, rotation covariance and domain/cancellation refusal. Test execution
and full workspace validation are separate from these specified regressions.

## Recorded stress histories in transient quadrupole observation

`curle3d::quadrupole::transient::history::PreparedQuadrupoleHistory` connects
original finite, physical six-component stress records [N m] to the existing
full transient pressure/vector-velocity propagator. It retains source positions
and order, uses the original `QuadrupoleSampleClock` timestamp admission, and
requires explicit first/last stress rates and the stress integral through the
first sample. All six components use one clamped C2 cubic reconstruction;
value, first/second derivatives and integral are evaluations of that same
piecewise polynomial. The initial integral persists through every interval
and query order, preserving the near-field velocity history.

One common positive tridiagonal factor and a reused six-channel DD workspace
prepare all sources in O(sources*samples) work. The original moments, scaled
curvatures and cumulative DD integrals are retained under explicit source,
sample and peak backing-buffer byte limits. Local coefficients are reconstructed
in constant work; warm history and acoustic observations allocate nothing.
Queries retain DD time through endpoint-anchored clock subtraction; the stated
inclusive arrivals map to the original source endpoints without a timing
tolerance. The original per-query strict-retardation admission remains.
Prepared reciprocals remove division from warm interpolation. There is no
finite-difference derivative in the observation loop or dense source solve.

`prepare_observer` and `bind_observer` attach only matching physical source
layouts. Binding intersects the history coverage of every direct and ideal
reflected path. An uncovered source/receiver time, changed source order,
incompatible echo-delay spread, malformed sample, budget overflow or
cancellation returns no partial owner or field. Images query the same physical
history at their own delays; pressure and vector velocity sum before p*u.

The interpolation is offline and uses the complete supplied record. It does
not certify source bandwidth, antialiasing, positivity, derivative error or
turbulence-source validity, and does not invent prehistory, silence, periodic
continuation or a post-record tail. Its source-clock and interpolation errors
require independent refinement. The original history-bound observer retains
homogeneous stationary lossless propagation; the molecular observer below
supplies a separate, explicitly initially resting dispersive observation path.
`tests/quadrupole_history.rs` specifies cubic reproduction with nonzero initial
integrals, original-clock refinement, full reflected fields, three momentum
equations, DD epoch translation and coverage/ownership/cancellation admission.

## Causal molecular transient stress propagation

`curle3d::quadrupole::transient::molecular::PreparedMolecularQuadrupoleObserver`
retains the original `RelaxationState`, full physical stress histories and every
direct, ideal-wall or wedge path. `RelaxationState::laplace_response` analytically
continues the ORIGINAL caloric reservoir sum into finite Re(s)>=0; s=-i omega
matches this exterior owner's exp(-i omega t) convention. Its stable excess
propagation rate preserves weak loss without subtracting nearly equal
wavenumbers. No fitted spectrum, effective bulk viscosity or classical shear
and conduction attenuation is added.

The observer extracts each exact frozen-speed r/c_inf wavefront and prepares
four bounded radial symbols. They act on the same Q, Qdot, Qddot and causal
integral used by the native full pressure/vector-velocity kernel, including
transverse near fields, tensor image transformations and signed interference.
Instantaneous p*u is formed after coherent path/source sums. The empty-reservoir
limit uses the complete original lossless field without FFT propagation error.

An initially resting call explicitly declares zero prehistory and unexcited
exterior/material state. Nonzero authored initial Q, Qdot or integral refuses;
a finite initial Qddot jump is admitted without a smooth-start order claim.
Native recorded histories retain authored endpoint data separately from
spline-solve roundoff. Geometry, source order, clock and the full finite output
horizon are checked; later source and received tails remain unspecified.

BDF2 convolution quadrature samples the analytic material response at
right-half-plane Laplace points. Contour inversion prepares causal weights,
then a separate zero-padded FFT computes linear convolution. The first N
weights suffice for the first N outputs without truncating a physical tail.
Time-step error, contour alias/roundoff, original interpolation/bandwidth and
cell-size errors require separate refinement. With transform length L, work is
O(paths*L*log L + paths*N); numerical storage/work and cancellation are explicitly
bounded. This is a finite-prefix observation cache, not a streaming predictor,
coupled radiation reaction or automatic error certificate.

Six focused regression groups specify independently eliminated molecular Green
derivatives and harmonic signs, lossless full-field/large-clock parity, BDF2
refinement against independent Bromwich inversion, coherent ideal boundaries
and cancellation, original recorded binding and no circular future wrap, and
explicit onset, ownership, clock, budget and cancellation refusals. Selective
formatter/diff checks and a separate developer numerical oracle check passed;
native Rust execution was unavailable.

## Causal molecular volume, force and stress composition

`curle3d::molecular::PreparedMolecularSourceObserver` composes original physical
volume injection Q [m^3/s], force ON the fluid F [N], and symmetric stress
moments [N m] in the existing molecular transient cache. `TransientAcousticSource`
reuses `CompactSourceStrength`, retains the original force impulse [N s], and
carries the native stress value, rate, acceleration and integral. Source positions,
order and every direct, ideal-wall or wedge delay remain those of the original
path owner.

For G=exp(-kappa*r)/(4*pi*r), volume and force satisfy
`p=rho*s*Q*G-F.grad(G)` and `u=-Q*grad(G)+Hessian(G)*F/(rho*s)`, with
kappa=s*sqrt(rho/K(s)) from the original reservoirs. Their F0/F1/F2 projections
join stress's F0..F3 projections in DD before the shared filters. Physical
pressure and all Cartesian velocity components, including transverse near
fields, combine across mechanisms and paths before instantaneous p*u. Force
receives no additional density factor. Empty-reservoir propagation recovers
the native complete lossless sources.

One original kernel bank, causal DD clock, contour and padded linear-convolution
engine serve the whole assembly. Only six signed volume/force scales per path
are added, charged to the same workspace-byte limit. FFT work bounds and
cancellation remain shared; no additional material sampling or convolution
bank is introduced. `stress_observer()` exposes the same cache for the existing
stress-only and recorded-stress APIs without rebuilding.

Calling initially-resting observation declares zero prehistory and unexcited
exterior/material state. At the epoch, Q, F, force impulse, stress, stress rate
and stress integral must vanish. Finite Qdot, Fdot and stress acceleration jumps
are accepted without a smooth-start order claim. Original derivative/integral
consistency and bandwidth remain caller-owned; later source/received tails are
unspecified. The model remains stationary, homogeneous molecular bulk relaxation
with inviscid momentum; classical shear/conduction losses, motion, finite-impedance
scattering and reciprocal source loading are not supplied.

Six G1/G3/G4 regression groups specify independent complete Green fields and
harmonic signs, native lossless/source-cache recovery on a large original clock,
full-field refinement against independent Bromwich inversion, ideal boundaries
and cross-mechanism interference, contour/future-pulse causality, and onset,
storage, clock and cancellation refusals. Selective formatting/diff checks,
bounded independent source review and a separate numerical fixture calculation
passed; native Rust execution remains unverified.

## Coherent FW-H surface and exterior-volume recordings

`curle3d::array::recorded::joint::JointSourceWelchStream` binds disjoint compact
and Lighthill-stress source layouts, one acoustic reference and one physical
clock to a SINGLE native Welch stream. Source-major `[Q,Fx,Fy,Fz]` blocks
(m3/s, N) precede `[xx,yy,zz,xy,xz,yz]` integrated stress blocks (N m), with
physical undoubled shear. The original mixed observer coherently sums both
mechanisms' pressure and full vector velocity before forming second moments.
Cross-source, cross-axis and cross-mechanism terms survive, with native phase
and peak-bin scaling applied once. Existing lossless, molecular and matching
ideal-image paths remain available. Finite-record factors are not deterministic
cross-frequency waveform phases or an assumption of independent segments.

`curle3d::quadrupole::flux::surface::exterior::ExteriorSurfaceWelchStream`
privately owns the original surface and volume extractors. Frozen-gas and
original mixture entry points read ONE borrowed accepted flow state per sample;
matching timestamps on separately supplied states cannot substitute for that
identity. Constructors include the current sample exactly once. All selected
volume cells must lie outside the original half-open control box. Face-adjacent
exterior cells remain valid; enclosed cells refuse instead of being radiated
a second time. Original gas/species, grid, boundaries, frame, acoustic reference,
source order and uniform recording clock remain bound. No CFD step, source
kernel, area/volume multiplier or source derivative is added.

The `from_molecular_frozen`/`from_molecular_mixture` constructors and matching
sample methods evaluate ONE original `StressField3d` per physical frame.
The same accepted primitives, exact conserved momentum and face stresses supply
both surface `-tau.n dA` and exterior-volume `-tau dV`. Source identity,
geometry, immutable stress mode and native/recording clocks are admitted before
the shared evaluation. This preserves matching viscous terms across the FW-H
partition without paying for a second full-grid constitutive/gradient pass.
Native transport controls, walls and selected SGS remain in their original
owner, and no accepted flow or closure history is advanced.


Every unsuccessful sample restores both accepted source histories and the
joint recording. A failed second extractor restores only the already committed
surface; a rejected spectral push restores both caches. The original Welch
transaction retains its prior factors, and no fallible action follows a
successful push. Warm observation adds no adapter allocation. Three explicit
byte caps bound both native caches and the joint recording; geometry is charged
to its recording allowance, source counts bound local extraction, and one native
Welch work budget covers all `4*S + 6*V` channels. The original CFD owner retains
the storage and constitutive work of a full-grid molecular evaluation. No dense covariance or full history is
constructed.

`RecordedExteriorSurfaceSpectrum` retains the original control box and frame.
Its preparation/evaluation checks use the actual frame inverse and reject
inside, on-surface and roundoff-near-surface receivers. The original mixed
propagator additionally requires identical receiver, frequency, medium,
molecular law and ideal boundary domains. The source scope remains stationary,
resolved and finite in its selected exterior extent. Euler entry points retain
inviscid stress; molecular entry points include the original configured viscous
and SGS stress. Missing SGS sound, exterior volume, moving geometry and source
feedback are not inferred. Spatial, temporal and source-model qualification
remain necessary.

Analytical regressions specify compact/stress interference with zero pressure
but nonzero particle velocity, compact-only parity under zero exterior stress,
and propagation/binding refusals. Native-source regressions use actual frozen
and ten-species mixture histories, joint batch factors/full fields, exterior
geometry and rotated receivers, source/clock refusal, both rollback paths and
exact retry. Pinned formatting and independent source/API/physics review passed.
These original Rust regressions were authored without a complete dependency
build in their introduction session. Added molecular groups specify exact
single-evaluation acquisition through an instrumented original transport law,
coherent factors, source-mode refusal, disjoint exterior selection, both
rollback paths, exact zero-law recovery and changing-composition source
agreement. Their execution status is reported with the associated change.

## Coherent recorded sound-power spectra

`RecordedJointSourceSpectrum::prepare_radiated_power` binds the original two
source layouts, material, clock and Welch factors to the existing analytic
`PreparedMultipoleRadiatedPower`. Each factor retains every volume/force/stress
cross term before spatial power integration. Native peak-bin scaling is applied
once. Results give complete-bin watts and one-sided W/Hz; repeated prepared
evaluation reuses one source workspace without allocation or angular quadrature.
`radiated_power_band` sums explicitly selected complete interior bins in DD,
admits the whole source-pair/factor/bin workload before allocation, and retains
one spatial operator at a time. DC and even-length Nyquist refuse; the odd last
bin remains physical. Cancellation publishes no partial band and permits retry.

The native `RecordedExteriorSurfaceSpectrum` delegates both operations to its
same privately retained joint record. The scope is stationary homogeneous
lossless free-space radiation, with the original disjoint-source, compact-cell,
finite exterior-volume, sampling, window and leakage limits. Band summation
does not reconstruct cross-frequency phases, missing cells, absorption,
scattering, source feedback or statistical confidence.

Five focused integration tests passed through a native exact-source leaf using
the original FFT/GeneralFft/Welch, acquisition/clock, material and analytic power
owners: coherent Q/trace cancellation, separated mixed-source sphere flux,
overlapping-window derivative Parseval, odd/even endpoint treatment, and bounded
cancellation/retry. SIMD dispatch used the original scalar kernels. Unrelated
CFD, molecular and receiver adapters were omitted from this leaf; full-workspace
and native CFD extraction integration are not claimed executed.

## Coherent accepted-flow histories for transient sound

`ExteriorSurfaceHistoryRecorder` retains the accepted real FW-H surface Q/F
and disjoint exterior Lighthill volume T histories from one immutable native
Euler or mixture-flow borrow. The original gas, composition, grid, ambient
reference, source order and uniform clock remain authoritative. Sampling uses
admitted buffers; a refused extraction or recording restores the preceding
accepted source caches and prefix. Flow is never advanced or resampled here.

`AcousticSourceHistoryRecorder::prepare_clamped_cubic` explicitly selects an
offline C2 reconstruction with authored endpoint rates and force/stress integral
constants. A source-band declaration must lie strictly below the original
Nyquist frequency; it does not certify bandwidth or derivative accuracy. One
private const-generic implementation preserves the original stress-history API
and its scaled tridiagonal solve, analytic derivatives and DD primitives.
Preparation borrows the completed record and can be retried after cancellation.

The prepared history binds the same source positions, clock and equilibrium
reference to the existing molecular Q/F/T propagator. Explicit initially-resting
observation retains original frozen-speed fronts, material memory, coherent
pressure and full particle velocity. The exact authored epoch supplies onset
conditions; incompatible zero-past data are refused. RecordedExteriorSurfaceHistory
retains the FW-H exterior receiver guard for both binding and observation.

Source reconstruction, finite exterior-volume selection, inviscid stationary
source assumptions and propagation convergence remain separate obligations.
This finite recorded prefix does not reconstruct Welch phases, infer unrecorded
prehistory, supply a later material tail or claim a calibrated live audio stream.

## Exact straight-edge plane-wave diffraction

`diffraction::HalfPlaneDiffraction` solves Sommerfeld's ideal half-plane
problem in the local physical frame `x >= 0, y = 0`, with its edge along z.
It binds the original `FluidAcousticState`, positive frequency and a unit
incident propagation direction; both rigid and pressure-release sheets are
available. Edge-parallel phase is separated analytically, so skew incidence
retains the physical transverse wavelength. Signed `y = +0` and `y = -0`
select the two distinct sheet faces; the open half-plane remains continuous.

The coherent total solution includes incident, reflected and diffracted
fields through the illuminated, transition and shadow regions. Analytic
Cartesian Euler velocity retains the `1/sqrt(r)` edge term. Pressure and
velocity feed the existing `CoherentCompactField` and `HarmonicFieldLine<3>`
observations without a visibility mask, second propagation delay or pressure
squared substitute for signed energy flux. The exact edge and exactly grazing
incidence refuse. This is an infinite, zero-thickness, stationary ideal sheet;
finite barriers, impedance, point-source incidence, arbitrary wedges, losses,
flow across the edge and source loading remain separate models.

The implementation follows Sommerfeld A3--A4 in
[the NIST half-plane study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4844525/).
Its small-argument Fresnel series uses DD accumulation; the large-argument
incomplete-gamma continued fraction ([DLMF 8.9.2](https://dlmf.nist.gov/8.9.E2))
retains the scaled tail directly. Shared deterministic phase reduction and
the analytic common edge phase avoid cancellation of two rounded large
phases. Work is allocation-free and bounded by two 256-term transition
evaluations, with cancellation at most eight iterations apart. Numerical
convergence is floating-point evidence, not an interval or null-relative bound.

`diffraction::tests` checks independent Fresnel/Faddeeva landmarks and full
Sommerfeld reference values, both faces' exact boundary conditions, analytic
velocity against independent pressure derivatives, the full Helmholtz residual,
illuminated/shadow/transition asymptotic limits, skew incidence, material density
scaling, observation-line parity and domain/cancellation refusals (G1/G2/G3).

### Exact diffraction with edge-parallel mean flow

`HalfPlaneDiffraction::with_edge_parallel_flow` and
`PreparedHalfPlaneSpectrum::with_edge_parallel_flow` extend the same exact
Sommerfeld owner to a stationary sheet in uniform subsonic flow along its edge.
For each incident wave normal n, the laboratory frequency omega obeys
`omega = c*k + U*kz`, with `kz = k*n_edge`. The intrinsic frequency is
`Omega = omega/(1 + U*n_edge/c)`. Translational symmetry retains kz in every
diffracted term, so the material derivative is exactly `-i*Omega`: both
pressure and physical particle velocity come from the original separated
Helmholtz/Euler solution. This is a direct reduction of the uniform-flow
linearized Euler equations, rather than a general-flow approximation.

The incident direction is a wave normal; laboratory group velocity is
`c*n + U`. All observation lines keep the original laboratory frequency.
World-space preparation uses each component's own spatial wave number when
moving its source phase to the edge, then retains coherent interference and
rotates the complete velocity into world coordinates. Different directions
at one laboratory frequency generally have different intrinsic frequencies.

Existing active intensity remains pressure work `<p*u>`.
`evaluate_energy_transport` reuses the original quadratic-energy owner to
return `<E>` and `<p*u> + U*<E>` after coherent assembly. It therefore retains
the advection of interference energy even when pressure-work divergence alone
is nonzero. Preparation remains O(N); observation is allocation-free O(N),
with the original cancellation checks and component bounds.

Zero flow preserves the original solution. The admitted speed is strictly
subsonic and parallel to the edge; both rigid and pressure-release faces remain
available. Cross-edge flow, shear, wakes/Kutta conditions, finite barriers,
impedance, losses and source feedback are not supplied. The new
`diffraction::flow_tests` exercise the original pressure/momentum equations,
both sheet faces, zero-flow equality, laboratory clocks, placed coherent
spectra, and energy-flux divergence (G1/G3). Independent Fresnel quadrature and
finite-difference equation checks were executed during implementation; the
native regression file and full crate integration still require a Cargo run
because the execution host became unavailable.

### Placed straight edges and coherent angular spectra

`diffraction::PreparedHalfPlaneSpectrum` reuses that Sommerfeld owner for every
supplied propagating plane-wave component. `HalfPlaneFrame` places the edge in
world coordinates with explicit origin, edge direction and sheet normal. The
finite nonzero axes are normalized; nonorthogonality beyond 64 epsilons refuses.
The incident complex peak pressures [Pa] are referenced to one explicit world
point. DD displacement/phase products move them to the edge origin exactly once,
so translating the boundary never silently resets the source phase.

All components' complex pressure and actual Euler velocity accumulate in DD.
The complete velocity is rotated into the world frame before active/reactive
intensity is formed; coherent interference cross terms and quiet residuals are
retained. `evaluate_harmonic_line` feeds the original finite-window and sampled
observers without a second delay or assumed pressure-to-velocity conversion.
Exact surface traces use `HalfPlaneReceiver::SheetFace` and a named face; world
points resolving exactly on the sheet refuse instead of inferring signed zero
after a rotation. Nearby points are never snapped onto the boundary.

Preparation and storage are O(N); each receiver is allocation-free O(N), with
the original transition cancellation checks and a caller component budget plus
a hard cap of 4096. There is no mode cutoff, implicit quadrature normalization
or angular interpolation. A continuous incident spectrum requires its own
angular-refinement evidence. The operator accepts only propagating components
whose resolved normal direction exceeds 64 epsilons in magnitude, refusing
unresolved grazing incidence after world projection. It adds no evanescent incidence, source reconstruction, finite
barrier, impedance edge or cross-edge/sheared-flow model. The angular-spectrum
convention follows [JASA 148, EL333 (2020), §2.1](https://pmc.ncbi.nlm.nih.gov/articles/PMC7567611/).

`diffraction::spectrum::tests` checks world rotations and phase-reference changes,
bitwise whole-world translation, coherent power cross terms, a 1e-20 quiet
component between cancelling unit fields, both rotated sheet traces, world
pressure-gradient velocity and original observation-line parity. A smooth
propagating angular spectrum converges over 16/32/64 components against a 1024
component reference; the finest combined pressure/impedance-scaled-velocity
absolute error is below 3e-5 Pa for that fixture (G1/G3). Geometry, budget,
edge/ambiguous-face and cancellation refusals are checked separately.

## Complete causal returns between opposing ideal walls

`curle3d::slab::PreparedSlabSourceObserver` observes the original stationary
volume, force and symmetric-stress histories between two infinite parallel
rigid or pressure-release planes. Both equal and mixed boundary conditions are
supported. The same source projections used by the molecular propagator retain
full near/far pressure, vector particle velocity, force impulse and stress
integral; this lossless owner introduces no molecular filter or other damping.
Image parity and vector/tensor orientation are applied to cached geometry.
Physical source indices remain unchanged and each image queries its own
retarded time on the original double-double clock.

An explicit initially-resting epoch and finite duration T make the infinite
image series finite without a bounce-order approximation. For slab width L,
the retained periods satisfy N > cT/(2L), including a rounding guard. Omitted
images have path length at least 2NL and therefore cannot arrive during the
declared window. Insufficient path budgets refuse the complete preparation;
later receiver times refuse rather than returning a truncated reverberant
field. Queries check compatible zero-onset data, observe no source before its
wavefront, and sum coherent fields before forming instantaneous p*u. Preparation
is O(paths) with bounded cancellation polls; evaluation allocates nothing and
uses cached coefficients with O(paths) work and no transcendental evaluation.

`tests/transient_slab_reflection.rs` contains G1/G3/G4 checks against an
independent, larger physical-source image sum, first and repeated echoes,
both ideal wall conditions for full force/stress orientation, rotated frames,
large epochs, window extension, coherent cancellation, and preparation/query
refusals. These are numerical checks of prescribed compact sources in a
stationary homogeneous lossless slab. Finite wall impedance, finite edges,
arbitrary rooms, moving sources/walls, radiation reaction, source bandwidth
qualification and measured-SPL validation are not supplied.

## Finite-distance spherical source diffraction by the original straight edge

`diffraction::PreparedHalfPlaneMonopole` binds one original volume-flow
source, placed ideal half-plane, receiver, material and positive frequency.
It resolves the full three-dimensional Biot--Tolstoy/Sommerfeld Green function
for rigid or pressure-release sheets, including spherical spreading, reactive
near fields, geometrical shadow transitions and either explicit sheet face.
The stationary-fluid source clock is retained: `p/Q=-i*omega*rho*G` and `u/Q=-grad(G)`
under exp(-i*omega*t). Velocity comes from differentiating the same field;
there is no plane-wave pressure/impedance replacement or independent ray gain.

The original impulse solution is Ouis (2002), Journal of Sound and Vibration
252(2), equations 3--6, DOI 10.1006/jsvi.2000.3553:
<https://faculty.kfupm.edu.sa/are/djamel/Publications/hard_half_plane_useful.pdf>.
An integration-by-parts form isolates the geometrical-boundary cancellation.
The remaining contour rotates by pi/4 without crossing its poles or square-root
cut. Dimensionless substitutions expose both shadow-transition widths and map
the complete infinite contour into finite numerical regions. Pressure and all
three analytic velocity components share each quadrature node.

Preparation reuses the existing adaptive Gauss--Kronrod engine, double-double
coherent sums and logarithmic physical scaling from spherical-wave refraction.
It respects the original positive pressure/vector-velocity tolerances, at most
1,048,576 evaluations, depth at most 48, fallible panel allocation and cancellation.
Tail envelopes control phase-resolution refinement only; they do not remove
a spectral region. Embedded error estimates exclude geometry and floating-point
errors and are not interval certificates. Quiet coherent cancellation is tested
against the complete physical field, including its geometrical contribution.
Diagnostics retain the work, component errors and minimum source-edge-receiver
path. Preparation publishes no partial transfer on refusal.

The immutable prepared observation caches four complex coefficients.
Changing the prescribed source phasor then requires four complex products,
with no allocation, quadrature or transcendental evaluation. The original
material, source, sheet frame, boundary and receiver selection remain queryable.
The existing `HarmonicFieldLine<3>` and `CoherentCompactField` outputs allow
coherent composition before intensity or time-domain observation. Occluded
image coincidence is finite; a visible source singularity, edge evaluation,
source on the closed sheet, or ambiguous sheet-face receiver refuses.

The default entry point uses a stationary homogeneous lossless fluid. Its
uniform edge-parallel-flow extension is described below. Neither entry point
adds a finite-barrier, impedance coating, arbitrary-wedge, moving-source,
absorption or reciprocal source-loading approximation.

Eight focused native tests compile the exact new implementation, original
plane-wave/frame owner, actual shared quadrature and fs-math, with the original
numerical material snapshot and narrowed error carriers. They check independent
static and free-space limits, analytic Euler velocity and Helmholtz residual,
reciprocity, both exact sheet conditions, occluded image coincidence, smooth
shadow-boundary gradients, convergence to the original plane-wave owner,
placed-source phase/velocity covariance, physical-scale preservation after raw
contour-decay underflow, and bounded refusals. This is focused source-slice
validation with the pinned compiler, not a full fs-aeroac/workspace Cargo result.

## Positive short-window pressure and velocity moments

Short uniform observations use the shared `fs-fft` positive temporal factor before forming pressure or velocity products: 3–32 real waveform tones when the full carrier span is admitted, or 3–64 proper-complex tones when their difference band is admitted, with at most 17 modes. Real modes retain sum frequencies and the true constant-mode mean; proper modes fold the peak-unit half into products before overflow. Every pressure/velocity diagonal and signed cross moment uses the same factor. Supplied Hz/rad/s clocks and source statistics remain unchanged.

Fixed scratch preserves allocation-free observation and the existing supplied-line/squared-line admission. Unsupported profiles, spans and counts retain the original pair kernels. Mode truncation and DD arithmetic do not certify arbitrary null-relative error. `tests/harmonic_observation_short_band.rs` exercises independent tiny-null waveform references, real/proper mixtures, native angular input, signed flux, cancellation and finite large physical moments.

## Short-window cross-receiver interference

Eligible uniform pair observations use one shared union-frequency temporal factor, with zeros for absent receiver frequencies, original order for self-pairs and canonical order for cross-pairs. Pressure/pressure, pressure/velocity, velocity/pressure and velocity/velocity moments retain signed products and source-factor identity. Real time-average covariance comes from the constant mode; proper time-average covariance integrates the original DD carrier fields before scaled-half outer products.

Queries remain allocation-free, with at most 64 analytic or 32 real union tones and 17 modes. Unsupported unions retain the original kernels. DD cancellation and finite modal truncation do not provide an arbitrary-null-relative or global array PSD certificate; different receiver pairs can select different temporal truncations. `tests/harmonic_array_short_band.rs` exercises self/local consistency, unequal-support phase and reciprocity, independent-factor zero covariance, physical-carrier averaging and finite peak normalization.

## Finite source diffraction in uniform edge-parallel flow

`PreparedHalfPlaneMonopole::with_edge_parallel_flow` extends the original
point-source owner to uniform subsonic mean fluid motion along the placed
straight edge. Source and sheet remain stationary; the original laboratory
frequency, material snapshot, source amplitude, and world geometry are retained.
The prescribed volume source enters continuity as Q*delta with no separate
force or vorticity injection.

For M=U/c and beta²=1-M², the exact potential Green function is
`G_U=exp(-i*k*M*dz/beta²)*G_0(x',s';k/beta)/beta`, with z'=z/beta.
The original ideal sheet is unchanged by this transform because the flow is
tangent to its faces and edge. The potential gives
`p/Q=rho*(-i*omega+U*partial_z)*G_U` and `u/Q=-grad(G_U)`.
Thus pressure and axial velocity mix during the Lorentz transform; substituting
one effective sound speed or one plane-wave intrinsic frequency is insufficient.
The governing transform and free-space fundamental solution are Proposition 2.1
and (2.2) of <https://www.waves.kit.edu/downloads/CRC1173_Preprint_2022-16.pdf>.

The existing stationary contour, both shadow-transition regions and analytic
gradients are reused. Every quadrature sample is transformed into actual
world-space pressure and velocity BEFORE the original error test. Conservative
back-propagated tolerances control only phase refinement of its tail.
Double-double contractions retain cancellation between pressure and axial
velocity contributions. The caller's evaluation, depth, allocation and
cancellation bounds remain authoritative. Diagnostics retain physical geometric
path length, not stretched coordinates. Zero flow uses the unchanged original
kernel arithmetic and reproduces its prepared transfer exactly.

`evaluate_energy_transport` uses the original coherent energy owner to report
mean perturbation energy and laboratory flux `<p*u>+U*<E>`; `evaluate`
continues to report signed active/reactive pressure work. Preparation adds no
larger solve or sampled edge spectrum, and subsequent source amplitudes still
use four cached complex products without allocation or quadrature.

Six focused G1/G3 flow tests and all eight original stationary tests passed
with the pinned July 6 compiler. Checks include the independent convected
spherical Green function, momentum/continuity/convected-Helmholtz residuals,
divergence-free mean energy transport away from the source, both sheet faces,
flow-reversed reciprocity, distant-source convergence in pressure AND velocity,
world/phase covariance, zero-flow parity and bounded refusals. This is an
exact-source numerical leaf using the original fs-math, field, quadrature,
frame and numerical material snapshot with narrowed error carriers; full
Cargo/workspace validation was not executed in the storage-constrained session.

Nonuniform/sheared or cross-edge flow, source motion, supersonic propagation,
finite geometry, thermal/viscous bulk attenuation, and source feedback remain
outside this admitted flow class.

## Physical finite-force diffraction by the original straight edge

`diffraction::PreparedHalfPlaneForce` resolves the complete Cartesian force
transfer through the original spherical half-plane Green function. The
prescribed source is a force ON the fluid [N], on its original harmonic clock:
`p=F.grad_source(G)` and `u=grad_receiver(p)/(i*omega*rho)` under
`exp(-i*omega*t)`. Its three pressure and nine velocity coefficients are
analytic derivatives of the same Biot--Tolstoy contour. Source displacement,
plane-wave impedance, a directional ray gain or an independently fitted
velocity never replaces that mixed receiver/source derivative. This extends
the finite-volume source owner using the same material, placed geometry,
ideal rigid/release boundaries, explicit sheet traces and physical field
carriers.

All twelve physical transfer channels share each original Gauss--Kronrod
node and deterministic DD accumulation. Force tolerances reuse the native
Euclidean pressure-transfer and Frobenius velocity-transfer limits in Pa/N
and m/s/N. Diagnostics retain Cartesian component error estimates, original
source-edge-receiver distance, evaluations, depth and panel count. These
estimates do not include geometry or floating-point error and are not relative
error certificates at a directional/coherent force null.

The mixed angular derivative retains finite integrated weight at arbitrarily
narrow geometrical transitions. Finite contour regions use logarithmic
coordinates between disparate transition/decay widths, so neither side of
that contribution can hide between initial quadrature nodes. At exactly
q=0 the visible geometrical Hessian has half weight, matching the continuous
sum with the one-sided diffraction integral. Force-specific tail envelopes
control phase refinement only; the complete infinite contour remains in
the integration, with original bounded depth/work and cancellation.

The immutable observer retains twelve complex coefficients. Changing the
original complex force vector requires twelve products and DD reductions,
without allocation, integration or transcendental evaluation. Output includes
pressure, full particle velocity, signed active/reactive intensity and the
existing perturbation-energy carrier. The tensor and force vector use world
Cartesian axes; both source and receiver transformations retain the original
phase and physical units.

Six new G1/G3 regressions check independent free spherical-force limits,
source derivatives of the original monopole in light and shadow, Euler and
continuity/energy balance, both exact sheet faces and occluded image coincidence,
shadow-transition continuity through offsets from 1e-7 to 1e-13 m, reciprocal
mixed tensors, source-clock/frame covariance, and physical/resource/cancellation
refusals. Together with the fourteen original stationary/flow point-source
tests, twenty tests passed on the pinned nightly-2026-07-06 native numerical
slice. It compiles the exact production derivative/point/flow kernels, original
frame/plane-wave and adaptive integration owners, actual fs-math and original
material numerical snapshot. Narrowed error carriers and excluded unrelated
source/material wrappers mean this is not full-crate, workspace or DSR proof.

This force owner is stationary and homogeneous/lossless. It does not inherit
the monopole's edge-flow model: a body force in moving fluid also has convected
vortical response. Finite barriers, wall impedance/absorption, moving force
sources, source feedback and measured acoustics remain outside this model.

## Physical volume-stress diffraction by the original straight edge

`diffraction::PreparedHalfPlaneQuadrupole` binds the original symmetric
Lighthill stress moment [N m], acoustic state, source position, placed edge,
receiver and positive frequency. It completes the ideal half-plane's physical
volume/force/stress source ladder. The original `QuadrupolePhasor` entries
xx, yy, zz, xy, xz, yz retain their trace, physical shear and emission clock.
The source Hessian supplies pressure and its receiving derivative supplies
all three Euler velocity components. No source-volume or density multiplier
is added to an already integrated stress. All source contributions and tensor
entries remain coherent before active/reactive intensity or energy is formed.

The original Biot--Tolstoy Green integral, geometry admission, DD propagation
phase and adaptive Gauss--Kronrod owner remain shared with the finite volume
and force propagators. Third derivatives require an additional exact
regularization at geometrical shadow boundaries: subtract
A*s/(1+s²/L²), A=(d/dL)[exp(i*k*L)/L]/L, from the integration-by-parts
integrand and integrate that term analytically. Its logarithmic primitive
cancels the visible Green term's q|q| contribution. The remaining contour
integrand vanishes cubically at its origin, so third derivatives are regular
at q=0 as well as on either one-sided approach. Small divided exponentials
and the visible Green remainder use analytic continuations; no displaced
source stencil, high-frequency ray approximation or numerical source
regularization changes the physical problem. The entire rotated contour is
integrated, including its algebraic subtraction tail.

A fixed, allocation-free 40-scalar derivative algebra retains exactly two
symmetric source derivatives and one receiving derivative. Geometry
derivatives are prepared in the original world axes and dimensionless edge
length; the physical length/momentum scales are applied once. The 24 shared
quadrature channels weight shear by sqrt(2), so the pressure and velocity
acceptance norms are their complete Frobenius norms. The retained 4-by-6
physical map instead includes the required factor two for physical shear.
Warm evaluation uses 24 complex products and fixed-order DD sums with no
allocation, new quadrature, transcendental evaluation or phase reset.

Diagnostics retain total evaluations, maximum depth, retained panels, shortest
edge path and absolute errors for every physical Cartesian transfer column.
For an authored moment, sum the column errors times the corresponding moment
magnitudes. These are embedded quadrature estimates, excluding geometry,
material and floating-point errors; they do not establish relative accuracy
at a tensor-pattern or coherent null. Explicit resource limits and bounded
cancellation publish no partial prepared observer. The original point and
force numerical kernels are unchanged.

Seven new physical regressions cover all six entries against independent
spherical Hessian/third-derivative fields, original force-source finite
differences in light/shadow and at q=0, the isotropic Helmholtz trace,
momentum/compression/energy flux, both exact sheet faces, occluded image
coincidence, shadow offsets through 1e-13 m, full world-tensor rotation,
original source phase, zero moments and invalid/cancelled/unresolved work.
All seven and the twenty existing stationary/flow point and force checks
passed in the focused native Rust slice with the pinned July 2026 compiler.
The slice compiled the complete new implementation, original quadrature,
geometry, diffraction and energy numerical code, and actual fs-math.
It narrowed unrelated crate/error/material-construction wrappers; full Cargo
and the repository DSR gate were unavailable and are not claimed.

The model remains a stationary prescribed compact stress in a homogeneous
lossless fluid around one infinite, ideal straight-edged sheet. Finite source
size and stress acquisition require their original source evidence. Finite
barriers, absorbing faces, nonuniform or mean flow, source feedback,
unresolved SGS stresses and measured-SPL qualification are outside this owner.

## Profile-aware observer and array interference

Local and array field projections now dispatch actual Uniform/Triangular
normalized temporal factors in the original Hz or rad/s source domain.
Pressure, intensity and supported cross-channel products retain coherent
field interference before products. Proper-factor averaged readout covariance
continues to integrate the original carrier separately from the relative-band
second moments. Existing union-spectrum bounds and fallback remain explicit.

Three focused actual-helper regressions pass with warnings denied, including
triangular dark leakage, real/proper statistics, retained source clocks and
unequal-support array correlations with independent integration. Original
ComplexSum, frequency, kernel and averaged-amplitude functions were executed.
The full aeroac source/observer pipeline and workspace suite were not run.

## Causal finite-source half-plane diffraction

`diffraction::PreparedHalfPlaneTransient` prepares the actual real
Biot--Tolstoy half-plane kernel for the original
`curle3d::moving::TimedCompactSource` volume-flow clock. The stationary
source and receiver share the existing placed `HalfPlaneFrame`,
`HalfPlaneReceiver`, ideal rigid/pressure-release boundary and original
`FluidAcousticState`. The original harmonic point owner and this time owner
share one private geometric lowering, including the regular shadow-boundary
coordinate and exact sheet-face semantics.

The owner accepts original Q and dQ/dt at uniform physical knots. They define
a C1 cubic Hermite source history; the first knot must declare Q=dQ/dt=0 and
every original force/force-rate component must be exactly zero. No force
input is silently dropped. The real retarded integral is analytically
differentiated to obtain pressure and all three particle-velocity components.
An integration-by-parts form remains regular at both shadow transitions;
time-basis knots and narrow geometric transition widths split the adaptive
integration. The eight real tap channels act on Q and dt*dQ/dt, retaining
their actual physical units. They are not sampled harmonic transfers or a
frequency-grid propagator.

Preparation retains every causal tap required by the declared observation
horizon (1..=65536 samples, including the zero start). No old source sample
or physical edge tail is discarded within that window. Longer records
refuse. Explicitly appended zero source knots observe the post-emission
tail within the same horizon; there is no inferred terminal time or claim
of indefinitely running finite memory. The complete record uses cached
zero-padded linear convolution through the original `fs-fft::GeneralFft`.
Its work is O(H log H) for the prepared horizon H, with eight cached
coefficient spectra. `observe_at` evaluates one accepted index directly
with double-double accumulation and no allocation. Complete-record
observation allocates its declared scratch and result vectors.

The source clock is checked against the original start plus i*dt with
double-double arithmetic; adjacent timestamps must remain representable
and strictly ordered. Causal-support comparisons use the prepared relative
index/dt clock, so rounding a large absolute epoch cannot erase a physical
first arrival. The mathematical convolution uses no later source knot. Known empty source prefixes retain exact zero fields up to the first
possible physical arrival, including a delayed emission prefix; this
does not truncate later edge memory or overwrite physical cancellation
nodes. Remaining FFT and floating-point roundoff is not covered by the
quadrature estimates.

Returned fields contain real pressure, full velocity, instantaneous p*u
without a phasor factor of one half, and the original perturbation energy
p^2/(2K)+rho*|u|^2/2. Coherent fields must be summed before forming
intensity or energy. Direct observations accumulate per-tap embedded
quadrature estimates against the actual source magnitudes. Batch
observations use the conservative sum of tap estimates times record-wide
maxima of |Q| and |dt*dQ/dt|. These estimates exclude source interpolation,
geometry, FFT and other floating-point error; source/time refinement is
an independent accuracy axis.

Preparation admits the original per-tap quadrature limits, a total
integrand-evaluation ceiling of at most 16777216, checked FFT extents and
an explicit conservative numerical-payload envelope. Caller-owned source
records and allocator metadata are outside that envelope. Cancellation
polls surround original synchronous FFT operations and occur throughout
coefficient preparation and physical contractions. Invalid geometry,
clock, material-independent source scope, work exhaustion, nonfinite
physical fields and cancellation publish no partial prepared owner or
record.

Seven new native regression groups cover:
- independent spherical pressure, full velocity and exact causal support
  for cubic histories, including delayed emission, zero records and large
  absolute epochs whose timestamp differences differ from the prepared dt;
- an independently integrated, untransformed real Biot--Tolstoy reference
  in light and shadow for both ideal boundary conditions;
- temporal Hermite refinement for an independently prescribed quartic
  history (observed relative errors 1.9174e-6, 1.4698e-7, 2.7287e-8);
- both exact sheet-face conditions and narrow shadow transitions;
- independently differenced momentum, continuity and instantaneous
  energy-flux equations;
- placed-frame rotation, original clock offset, linearity, independence
  from later source samples, and a nonzero post-emission edge tail;
- source, clock, horizon, payload, quadrature-work and cancellation
  refusals with exact retry.

The final focused native run passes these seven groups plus all 27
existing point/flow/force/quadrupole diffraction regressions. It uses the
actual existing diffraction/quadrature owners, original fs-math, original
one-dimensional FFT algorithms and original scalar SIMD leaves. The
narrow build omits unrelated executor/material-receipt integration and
does not establish full Cargo/workspace compilation, ISA speed, or a
measured performance ratio.

Physical scope remains stationary, linear, homogeneous and lossless
acoustics around one infinite ideal half-plane with an open-fluid volume
source. This does not supply moving-source or convected time-domain
diffraction, finite barriers, impedance edges, thermal/molecular
absorption, stress-quadrupole time histories, structural feedback or an
experimental SPL certificate. Those require their own original physical
owners and evidence.

## Two-tone local and cross-receiver projections

The bounded local and union-spectrum array factors now include two real or
proper source frequencies in native Hz or rad/s. Original-carrier averaged
covariance, signed cross moments, temporal profiles and one-tone fallback are
preserved. New physical null and cross-receiver checks join the existing profile
tests; the combined duct/aeroac actual-helper execution passes eight checks with
warnings denied. Full source/observer pipeline integration remains unrun.

## Causal edge diffraction with original mean-flow transport

`PreparedHalfPlaneTransient::with_edge_parallel_flow` extends the preceding
causal owner to a uniform strictly subsonic mean fluid velocity parallel to
the infinite edge. The stationary scope of `new` remains unchanged; zero
flow retains its original tap and quadrature arithmetic. Source and receiver
stay fixed in the laboratory. Original Q and dQ/dt, absolute epoch, relative
sample clock, original density and adiabatic modulus remain authoritative.

The actual causal Lorentz image stretches the axial separation by 1/beta,
where beta^2=1-(U/c)^2. A physical path has delay
(sqrt(beta^2*r_perp^2+dz^2)-(U/c)*dz)/(beta^2*c). The downstream expression is
rationalized before evaluation. Around the first edge arrival, extra delay
uses sqrt(1+w^2)-1=w^2/(sqrt(1+w^2)+1), avoiding subtraction of a large
stretched flight time and nearly equal clock advance. Source-basis support
boundaries are computed from that same physical arrival.

The image couples pressure and axial particle velocity, using the same
physical transform as the existing harmonic point-flow owner. Every
integrand sample and exact direct/edge term is transformed before the
original physical pressure/vector-velocity quadrature norms are tested.
The internal Lorentz operator speed does not replace the original material
state or become a shifted source frequency. World-frame rotations retain
the actual signed mean-fluid velocity along the authored edge.

The same eight real kernels, complete declared-window memory, native FFT
linear convolution, direct observations and explicit total-work/storage
limits remain in use. No source samples or physical post-emission edge
tail are dropped. Flow does not authorize force injection, source motion,
cross-edge/sheared convection, finite barriers, impedance edges or loss.

`HalfPlaneTransientField::instantaneous_energy_flux_w_m2` now reports
laboratory perturbation transport p*u+U*E in addition to pressure work p*u.
The existing original `curle3d::energy::instantaneous` owner supplies energy
and both fluxes for stationary and flowing fields. This is full quadratic
perturbation energy, not exclusively radiated sound power. A moving control
surface would require its own relative flux. The original material
coefficients are retained and coherent fields must be assembled before
the quadratic reduction.

Seven new G1/G3 regression groups verify the unchanged zero-flow taps;
independent real convected free-space fields and arrival at Mach numbers
-0.75, 0.4, 0.85 and 0.99; the real unintegrated-by-parts Biot--Tolstoy
potential with independent space/time differentiation; temporal source
refinement; linearized-Euler momentum, continuity and advected-energy
balance; both exact sheet-face conditions and shadow transitions; placed
world frames, original clocks, source cessation and future-source
independence; and physical/budget/cancellation refusals with exact retry.
On the declared quartic-history fixture the relative field errors refine
as 1.7043e-6, 1.6904e-7 and 1.6016e-8.

The final focused native run passes these seven groups and all seven
existing causal diffraction groups (14 total). It compiles the actual
causal operators, original geometry/frame and Gauss--Kronrod numerical
owners, original acoustic material snapshot and energy reduction, actual
fs-math and original FFT/scalar-SIMD algorithms. Unrelated material
receipts and executor/dependency construction are outside this narrow
build. This is not full Cargo/workspace or DSR success, a near-sonic
relative-error certificate, or a measured performance benchmark.

## Causal physical-force diffraction and retained impulse

`diffraction::PreparedHalfPlaneForceTransient` propagates original stationary
force-on-fluid records through the same real Biot--Tolstoy half-plane Green
function as the volume-source owner. It evaluates
`p = grad_source G[F]` and
`u = -grad_receiver grad_source G[integral(F dt)] / rho`.
Analytic derivatives precede the radial substitution; the shared original
mixed source/receiver geometry also serves the harmonic physical-force owner.
No displaced source pair, harmonic reconstruction, fitted attenuation or
pressure-to-impedance velocity approximation is introduced.

The original uniform source clock carries all three force and force-rate
components. C1 cubic Hermite interpolation has an analytically integrated
primitive: value bases retain their unit impulse beyond support, while rate
bases retain zero total impulse. The complete finite-horizon response includes
post-emission near-field velocity and the causal edge tail. All original world
source axes and receiver velocity axes are rotated independently. Rigid and
pressure-release sheet faces, exact shadow boundaries and arbitrarily narrow
resolved shadow transitions retain the original half-weighted mixed terms.

Preparation caches 24 real lag coefficients and corresponding original FFT
spectra. Record observation uses six input transforms and zero-padded linear
convolution, with no circular wrap or tail cutoff. One-knot observation is
allocation-free. Both paths use signed DD reductions and the original coherent
perturbation-energy owner. Exact pre-arrival zeros are retained; elsewhere the
FFT retains its ordinary floating-point error. Explicit per-tap/total work,
finite-history and peak-payload limits cover all force directions and the
original FFT workspace. Cancellation and refusal publish no partial owner or
record. Input must start quiescent; volume injection and incompatible clocks
refuse rather than being silently discarded.

Seven new physical groups pass alongside all fourteen original stationary/
convected transient groups and six harmonic-force groups. The independent
free-space pulse and retained-impulse reference agrees within 4.469e-14 in the
direct observation; FFT error normalized by each complete record component is
3.432e-15, including pressure zeros. Independent unintegrated real Green
derivatives agree within 2.140e-11 in the sampled light/shadow cases. Additional
tests cover time refinement, all force directions, both sheet faces, narrow
transitions, placed frames, causal future independence, resources/cancellation
and original momentum, continuity and energy equations. The largest sampled
finite-difference PDE defect was 1.238e-7.

The focused native build executes unchanged production numerical owners with
the original gas material, one-dimensional FFT, mixed geometry, adaptive
quadrature and energy arithmetic; unrelated modules/dependency surfaces are
omitted. This is not a full Cargo/workspace test or experimental validation.
Quadrature estimates exclude interpolation, geometry and floating-point error;
time refinement remains required. This force owner is a stationary homogeneous
lossless-fluid model without source motion, mean flow, finite wall impedance,
finite barriers, source feedback or measured-SPL qualification.
