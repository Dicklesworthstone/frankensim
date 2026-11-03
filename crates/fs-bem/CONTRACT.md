# fs-bem CONTRACT

## Exterior observation FIRs and CQ (IRP M2.5-W6)

`observation_transfers::BemObservation` binds the existing dense Burton–Miller
owner to every supplied signed boundary-velocity trace and explicit exterior
receiver. `SourceId` distinguishes structural-coordinate velocity from aperture
volume flow; the latter integrates to one SI volume flow using the original
panel areas. Sources above 256 are partitioned into complete batches, without
truncation. Every source/receiver retains pressure and all three Cartesian
pressure-gradient components. These are the same centroid Green rows and
analytic normal derivatives used by the Helmholtz owner; the near-field owner
supplies exact closed-surface/exterior/clearance admission. This route does not
claim adaptively integrated near-field quadrature accuracy.

`bake_observation_transfers` uses `fs-fft` on four to eight fixed doubling DFT
grids (32..65536). Fourier data use exp(-iwt) and are conjugated exactly once
at the DFT boundary. Propagation has an exact integer prefix; its fractional
part remains in the taps. The caller explicitly declares sampling rate, a
pass band, a C4 transition to zero strictly before Nyquist, reconstruction
latency, source/receiver identities, geometric earliest arrivals, producer
provenance/uncertainty, and an entry budget. The taper is exactly one in the
pass band. Output time includes the declared additional filter latency. No
passivity enforcement, allpass approximation or rational fit is applied.

Every DFT rung records worst-quantity tail energy after N/4, wrap-around energy
in the second half, and the maximum per-quantity relative L2 DTFT change on a
fixed 127-point withheld pass-band grid. A quantity is the pressure [Pa] or the
whole Cartesian gradient vector [Pa/m] of one source/receiver pair; Pa and Pa/m
are never mixed in a norm. Per-component fractions are frame-dependent: an
on-axis receiver's zero-by-symmetry components hold only roundoff.
Only the first half is retained, after measuring the full untrimmed inverse.
The final tail must meet the caller's <=1e-10 tolerance, wrap must be <1e-10,
and inter-rung DTFT change must be strictly below the caller's <=1e-4 tolerance.
Exhaustion returns `NotConverged`, with no partial bake or clipping repair.
Zero-frequency pressure/gradient for a static prescribed velocity are exactly
zero from q=-rho*s*v and the exterior zero-pressure condition at infinity.

`fs-vfit::observation::{FirKernel,ObservationFir}` is the portable data/replay
owner. Receiver engines can consume the FIRs without depending on fs-bem.
M2.6-W1 must retain the explicit p/gradient ordering, SI source-coordinate
units, pass band, geometric delays, latency and uncertainty; no microphone or
atmosphere defaults are supplied here. Actual M2.6 receiver integration remains
with that consumer bead.

`helmholtz::laplace` analytically continues the SAME point kernels, self terms,
static row closures and dense Burton–Miller assembly to k=i*s/c, Re(s)>=0;
q=-rho*s*v. The real-f64 kernel/self arithmetic is preserved, including the
iterative owner's inherited kernels. The legacy CQ adapter does not select
the routed triangle owner; the separate complex triangle assembly is used by
the molecular exterior below. Dense/resolution/batch caps still apply at |k|. A known
propagation delay may be extracted inside each Green exponential, bounded by
the actual surface clearance/c, to avoid exp(-s*tau)*exp(s*tau) under/overflow.
No radiated-power or passivity claim is made at a complex frequency.

`cq_time_response` is a BDF2 reference, with delta(zeta)=3/2-2*zeta+zeta^2/2
and an explicitly declared damped contour. It evaluates the original
Laplace-domain operator, not a fit. Its uncorrected second-order claim requires
smooth onset and zero prehistory. See Lubich, Numer. Math. 52 (1988), and
Hassell–Sayas, sections 2–3, <https://arxiv.org/abs/1407.0345>. CQ contour error,
time error, FIR truncation and BEM spatial error remain separate.

All cold operations are deterministic on one ISA. Cancellation is checked
before/between frequency/contour solves; a dense solve drains before the next
poll, with no mid-LU latency guarantee. Frequency storage, taps and input
history have separate 4194304-entry caps; source/output caps are 4096/256.
Failed FIR replay is transactional and leaves its shared histories unchanged.

`ProducerUncertainty` preserves exact analytic, provenance-bearing relative
estimate, or explicitly unresolved spatial error. BEM refuses ExactAnalytic.
An unresolved BEM bake only qualifies realization of that discrete operator;
it does not establish continuum or instrument accuracy. The strict 1e-6
sphere realization test uses the independent closed-form operator and time
response. Sphere/lathed-bell FIR-versus-CQ tests use actual identical BEM
panels, both source kinds and p/gradient, with four time rungs and fs-mms's
second-order gate. These do not upgrade the upstream BEM mesh qualification.
No real-instrument full-band claim is made; central native results are pending.

Focused targets: `tests/observation_transfers.rs` (G0/G1/G2/G3/G5),
`helmholtz::real_kernel_bits` and `helmholtz::laplace::tests` (G0), plus
`fs-vfit::observation::tests` (G0). Receipts are emitted by these tests; this
bead has no render-producing default or separate e2e framework.

## Material molecular exterior observation

`observation_transfers::MolecularBemObservation` binds the original
`helmholtz::molecular::MolecularExterior` to the same signed structural and
unit-volume-flow source basis, finite receivers and bounded Fourier owner.
One evaluated `RelaxationState` supplies the complex wavenumber to both the
selected boundary solve and the receiver Green representation. Neumann
data remain `i rho omega v`; `impedance_row` projects those same solved pressures
against the original receiving trace and panel areas. Thus an omitted-load
receipt and a microphone cannot silently use different material compression.

The default `TriangleCbie` formulation requires the whole band below the
geometry's conservative CBIE resonance-exclusion bound, using equilibrium sound
speed as a conservative bound on wavenumber for positive reservoirs.
`new_with_formulation` also accepts explicit `CentroidBurtonMiller`, retaining
the original complex centroid/disc combined equation throughout the entire
band, including low frequencies. That arm crosses the CBIE exclusion bound
without switching to or from triangle quadrature. Its low-frequency resistance
can be spatially biased; no passivity or spatial-accuracy certificate follows.
The original assembly reserves respectively 512 or 8 kernel terms per panel
pair. Source batches, dense storage, declared assembly work, wavelength and
receiver-clearance admission remain bounded by the original owners. Material,
boundary, LU and Fourier cancellation refusals are preserved. The FIR prefix
uses geometric clearance divided by frozen sound speed; sampled fields already
contain the complete dispersive propagation phase, without a second delay.

Operator identity includes the original gas, frozen heat capacity and every
evaluated species/reservoir capacity and relaxation time, as well as boundary
geometry, signed traces, fixed formulation and declared band. Different molecular
rates cannot reuse a bake merely because equilibrium density and sound speed match.
This is relaxing compression with inviscid momentum, not exterior shear or
thermal diffusion, and does not extend the reciprocal Galerkin reaction owner.
The existing centroid receiver representation and numerical FIR convergence
carry no new spatial or measured-instrument accuracy claim. Focused regressions
in `observation_transfers/molecular/tests.rs` compare reaction, distance phase,
pressure changes, empty-reservoir recovery, causal prefixes and admission.
An independent pulsating-sphere case spans low frequency and the first
interior Dirichlet resonance under one combined operator; formulation changes
also invalidate a cached transfer while retaining the same material and arrival.

`MolecularBemObservation::prepare_harmonic` binds those same complex-wavenumber
boundary and receiver fields to the shared coherent and partially coherent
observation bank. Particle velocity follows `grad(p)/(i rho omega)` from the
original inviscid momentum law. Intensity and microphone cross spectra retain
source interference, dispersive phase and reactive velocity without another
absorption correction or boundary solve during evaluation. Ideal first-order
microphones retain their real equilibrium `rho*c` pressure-equivalent
calibration; that reference is not the frequency-dependent wave impedance, so
unity plane-wave sensitivity in a dispersive medium is not implied. Material
band, cancellation and unresolved spatial-uncertainty boundaries remain.
G1/G3 tests check independently differentiated spherical waves, radial
transported-power attenuation, mixed structural/volume-flow interference,
microphone cross-phase and empty-reservoir recovery.

## Partially coherent harmonic observations and microphone arrays

`BemObservation`, `ReciprocalBemObservation`, `ScatteringObservation`,
`LocalBoundaryObservation` and `IterativeBoundaryObservation` expose
`prepare_harmonic` at one positive angular
frequency. This caches each existing owner's complete source-to-receiver
pressure and Cartesian particle-velocity transfers. The scattering owner
combines incident and scattered amplitudes before any quadratic observation;
the local-wall owner uses its physical, unwindowed constitutive response and
retains its material-band refusals. These are the same physical operators used
by the existing observation pipeline, with no new propagation/load model.
The matrix-free adapter retains all original-equation source solve reports
alongside the prepared field bank and introduces no dense fallback or extra
GMRES solve during coherence evaluation.

`PreparedHarmonicObservation` accepts coherent peak source phasors, explicit
mode-major factors B of the source second moment C = B B*, or diagonal peak
second moments for independent sources. Source coordinates retain their own
declared SI units, including mixed structural-mode and aperture-flow bases.
Output contains pressure mean square, the full Hermitian three-axis velocity
second moment, and signed complex intensity E[p conj(u)]/2. The negative-time
phasor convention fixes the reactive sign. Pressure power is never substituted
for transported intensity in a standing wave or reactive near field.

`microphones` prepares the existing physical first-order pressure/velocity
projection at each receiver, then forms the complete receiver cross-spectrum
C_ij = E[m_i conj(m_j)]/2 in Pa². Cross-receiver phase and pressure/velocity
interference remain available for arrays and beamforming. Independent sources
can still produce correlated receivers. Factors imply positive-semidefinite
statistics up to arithmetic roundoff; no eigenvalue repair or coherence guess
is applied. Zero rank explicitly means silence.

Preparation owns allocation and the original bounded BEM solve. Applications
reuse transfers, compensated signed reductions and microphone workspaces with
no allocation. Field factor work is O(sources × receivers × rank), independent
field work O(sources × receivers), and array factor work
O(rank × (sources × receivers + receivers²)). Explicit product budgets, a
4096-factor cap and the existing 64-receiver cap bound work before contraction;
cancellation and nonfinite statistics refuse. Discard caller output on error.
These are spectral-line averages, not PSD/Hz, a source-coherence estimator, an
in-loop load, or a spatial/model accuracy certificate: producer uncertainty
remains explicitly unresolved. G0/G1/G3/G5 tests cover a standing wave, dense
covariance and beamformer references, exact cancellation, diagonal factors,
weak signed residuals, and independently solved finite-source diffraction.

## Purpose and layer

Layer: L3 (FLUX). Laplace BEM panel methods (plan §8.3 [F], bead
tfz.20): potential-flow screening for exterior aerodynamics — the
ornithoid flagship's wide-search stage. INVISCID HONESTY LABELS apply
everywhere: this is screening, not a viscous truth source.

## Public types and semantics

- `Medium::from(&GasState)` copies the caller's density and sound speed bits
  directly. It evaluates no default atmosphere and adds no source authority;
  NASA-9/source receipts remain with the resolved material state. G3 tests pin
  both copied values, including values distinct from the air convenience.

- `panel3d`: validated `SpherePanels` (centroid/normal/area panelization of
  fs-rep-mesh icospheres, with optional exact oriented triangle vertices);
  `dense_matrix` — the collocation Neumann
  operator with the outside-limit jump −σ/2 on the diagonal and
  centroid-monopole off-diagonal rows (screening-grade; measured
  convergence is the gate); `fmm_matvec` — the SAME operator through
  three fs-fmm gradient-kernel passes dotted with target normals;
  `fmm_transpose_matvec` — the adjoint operator through the same FMM
  kernels, with panel-area placement and gradient antisymmetry tested
  against the dense transpose;
  `solve_exterior` — GMRES over the FMM matvec, returning a converged
  `ExteriorSolution` with its full `SolveReport` or an
  `ExteriorSolveError::NotConverged` carrying the last iterate and report;
  every returned report is a `ResidualClaim::TrueEuclidean` one — GMRES
  recomputes `‖b − Ax‖₂/‖b‖₂`, and the zero-onset-flow short circuit is
  built through `SolveReport::from_claim` with an EXACTLY zero residual
  rather than a bare struct literal, so `report.euclidean_rel_residual()`
  is `Some` on both paths (bead `frankensim-extreal-program-f85xj.2.24`);
  a fallible operator refusal is distinct and preserves its first `BemError`
  plus the iterate/report rather than laundering the refusal into convergence;
  `surface_velocity`.
- `panel2d`: `Airfoil2d` + `naca4_symmetric`; Hess–Smith `solve` —
  constant sources per panel plus one shared vortex density, the KUTTA
  row closing the system (equal tangential speeds leaving the two
  trailing-edge panels; circulation DETERMINED, not assumed); lift by
  PRESSURE INTEGRATION of the enforced surface field (the Γ-accounting
  shortcut was measurably wrong bookkeeping and is gone);
  `Cl = L'/(q∞ c)` uses the distance from trailing-edge node zero to the
  farthest leading-edge node as its reference chord. Both chord endpoints must
  be present. Pressure coefficients, lift, and the adjoint are invariant under
  uniform scaling and joint rigid rotation/translation of section and flow;
  `dcl_dalpha_adjoint` — one transposed solve for the solution
  sensitivity plus solve-free finite-difference output partials, FD-gated
  and explicitly not claimed as an exact symbolic derivative. This historical
  adjoint path and its arithmetic remain unchanged.
  `solve_with_lift_slope -> LiftSolution2d { solution, dcl_dalpha }` is a joint
  steady/tangent API: one geometry-only Hess–Smith LU solves both the primal
  RHS and its analytic angle derivative. The tangential velocity and its
  tangent share the same panel integrals; `Cp' = -2*vt*vt'`, and the lift
  derivative includes rotation of the force projection as well as pressure
  changes. It adds O(n^2) tangent work and O(n) storage to the ordinary solve,
  with no perturbed solves or finite-difference output partials. The steady
  `sources`, `gamma`, `vt` and `cl` retain the ordinary solve's bits.
  Non-finite tangent solutions and output derivatives return typed errors.
  This is the derivative of the admitted discrete inviscid model at fixed
  geometry, not a viscous, stall or shape derivative.
  The constant-panel integrals carry a battery-pinned lesson: the normal
  component is (θ₂−θ₁)/2π — the reversed order self-cancels a closed
  sheet's field (caught by the single-panel-vs-quadrature and
  uniform-sheet probes). `solve_naca0012_prestall` is the narrower G2
  validation entry point: it constructs the unit-chord NACA 0012 section and
  refuses non-finite angles or `|alpha| > 10 degrees` before allocation so an
  inviscid solve cannot be mistaken for a stall prediction.
- `wake2d`: `WakeSim` — impulsive-start free wake; Kelvin-conserving
  trailing-edge shedding, regularized point-vortex convection, the
  quasi-steady bound circulation relaxing against wake downwash;
  ledgered traces. The global freestream direction is resolved relative to the
  actual chord. Wake positions are in the local unit-chord frame (LE = (0, 0),
  TE = (1, 0)), circulation is scaled by U∞c, and time by c/U∞. The
  quarter-chord bound vortex and three-quarter-chord control point therefore
  move and scale with the supplied section, preserving its nondimensional lift
  history. A cached panel lift slope supplies a bounded harmonic response to
  local wake downwash: `Γ_target = Γ_steady cos(δα) + (dCl/dα)/2 sin(δα)`.
  The response retains camber lift at zero incidence and the section's actual
  zero-lift angle without dividing by `sin(alpha)`. This is a local
  attached-flow, quasi-steady screening closure; it is not a resolved unsteady
  pressure, stall, or separated-flow model. Construction obtains steady lift
  and its analytic slope from the joint API using one dense factorization;
  the prior separate steady/adjoint path required three factorizations and
  cubic finite-difference output work. There are no panel solves or
  sensitivity calculations inside wake steps. Analytic tangents remove the
  finite-difference noise from transformed and near-zero-incidence wake
  histories. The battery checks independent five-point derivatives for
  symmetric/cambered sections at signed angles, bitwise joint/ordinary
  steady equality, and rigid-frame/scale invariance. The wake-history
  similarity tolerance is tightened from 2e-7 to 2e-10.

### `helmholtz` (bead frankensim-fsim-helmholtz-bem-k1ryv)

- `solve_radiation(surface, k, medium, velocity, formulation)` — exterior
  Neumann radiation on a closed panel surface: complex surface pressure,
  radiated power, panels-per-wavelength diagnostic. Time convention
  `e^{-i omega t}`, `G = e^{ikr}/(4 pi r)`; with the module's HBIE row
  sign, the Burton-Miller coupling is `alpha = +i/k`. It is pinned by
  measurement, not by convention label: with the exact static rows it tracks
  plain-CBIE accuracy and convergence away from interior resonances, while
  `-i/k` is 3-5x less accurate for `ka >= 2`
  (`burton_miller_sign_tracks_cbie_accuracy_off_resonance`).
- `solve_radiation_batch(surface, k, medium, velocity_fields, formulation)` —
  the same solve semantics for a nonempty, bounded batch (at most 256 fields)
  at one frequency, preserving input order while sharing one validated dense
  matrix, LU factorization, and condition diagnostic across every right-hand
  side.
- `PreparedCbieGeometry` — optional immutable triangle preparation for repeated
  plain-CBIE batches. Exact static single/double layers occupy `16 * n^2`
  logical bytes, admitted against a caller-supplied byte budget only after the
  highest planned wavenumber passes the ordinary dense/resolution caps.
  Reservation is fallible; a cancellation callback is polled before allocation
  and between target rows. The borrowed surface cannot be cross-wired. Every
  `solve_batch` retains the original parameter/field/resolution checks and
  builds and factors a fresh frequency-dependent operator with unchanged
  dynamic quadrature. Existing one-shot entry points do not allocate this cache.
- `Formulation::{PlainCbie, BurtonMiller, BurtonMillerWrongAlphaSign}` —
  measured roles: PlainCbie is the accurate non-resonant arm (1.7-3.4%
  on the pulsating sphere across ka in [0.05, 5]); BurtonMiller is the
  resonance-safe production arm for ka >= 0.5 (1.1-5.2%; 1.7% at the
  ka = pi fictitious frequency where PlainCbie hits 50%). The wrong-sign
  arm exists so the mutation fails loudly (measured 3.3x worse at
  resonance).
- `radiation_impedance_matrix` — one factorization, n unit-velocity
  solves; feeds the vibroacoustic-coupling bead.
- `far_field` — centroid-rule sampled directivity amplitudes (monopole
  uniform to 0.00%, dipole cos-theta correlation 1.0000 measured).
  `far_field_integrated` integrates the retained triangles' spatial phase
  with a 4x4 Gauss rule; centroid-only surfaces retain the original rule.
  Both return `F` with the same `p -> F e^{ikr}/r` convention.
- `FarFieldTable` / `tabulate_far_field` — solver-neutral
  `(ω, direction) → F` samples with `p → F e^{ikr}/r`. Downstream
  observers (fs-couple) consume the table; they do not take a
  production dependency on this crate.
- `exterior_pressure_at_points` — deterministic batched finite-distance
  centroid pressure via `SUM (dG/dn_y p - G i omega rho v) A`; mismatched source
  surface or medium and invalid/non-exterior points are refused.
  `exterior_pressure_at_points_integrated` shares this admission and integrates
  retained triangles: exact static terms plus 8x8 regular remainders near a
  panel, or 4x4 complete Green kernels for well-separated points. The latter
  avoids cancellation in static edge formulas at large receiver distances.
  Centroid-only surfaces retain the original rule. Percussion's default far
  and finite receivers explicitly select these integrated entry points;
  existing centroid consumers and the adaptive `near_field` owner are unchanged.
- `directivity_sh_table` / `DirectivityTable` — far field projected
  onto orthonormal complex spherical harmonics (Condon–Shortley,
  `Y_{l,-m} = (-1)^m conj(Y_lm)`) up to `l_max <= MAX_SH_DEGREE = 64`
  by independently refined Gauss–Legendre x uniform-phi quadrature.
  Retained SH degree does not select a source bandlimit: the former minimum
  grid aliased an omitted x-dipole into a 100%-captured monopole at degree zero.
  `DirectivityProjectionOptions` declares total direction work and relative
  tolerance; two successive refinements must agree in both retained
  coefficients and full angular power before any table is published.
  Default tolerance is 1e-10, with at most 262144 source directions over the
  entire refinement. `DirectivityNotConverged` retains the spent budget and
  last changes when another complete grid cannot fit. Explicit-option entry
  points poll cancellation between latitude rows and source batches of at most
  256 directions; an opaque evaluator owns its in-call latency.
  Scale-normalized compensated power and coefficient sums preserve diagnostic
  ratios for tiny/large finite sources. Separable longitude Fourier sums and
  latitude projection avoid a dense SH design matrix and reduce accumulation
  from O(directions * degree^2) to O(directions * degree + latitudes * degree^2).
  Truncation is the reported `captured_fraction`, and
  `evaluate`/`power_by_degree`/`coefficient` serve the rendering and
  diagnostics consumers. Measured at ka = 1 on 320 panels: monopole
  l = 0 fraction and captured fraction > 0.999, dipole (l=1, m=0)
  fraction > 0.99, reconstruction vs directly sampled far field
  < 1e-3 relative L2 at held-out directions.
  The focused `sh_projection_*` tests additionally pin omitted polar and
  azimuthal multipoles, complex phase/normalization, source scaling and
  rotations, coherent displaced-source interference, and work/cancellation
  refusals. These are G1/G3 numerical identities and G4 bounded-work checks,
  not real-instrument qualification or a continuum error certificate.
- `radiation_efficiency` — `sigma = W / (1/2 rho c INT |v_n|^2 dS)`
  per solved velocity pattern (the per-mode radiated-power ratio).
  Oracles: monopole `(ka)^2/(1+(ka)^2)` within 5% at ka in
  {0.5, 1, 2} (PlainCbie arm), dipole `(ka)^4/(4+(ka)^4)` within 10%
  at ka = 1.
- Per-solve diagnostics on `RadiationSolution`:
  `condition_lower_bound` (a probe-based LOWER bound on the 1-norm
  condition number — a warning signal, not a certificate) and
  `dense_cap_utilization` (panel count over the dense cap, the
  headroom the aperture pilot asked for).
  `radiated_power_roundoff_interval` is an outward-rounded enclosure of the
  exact panel-power dot product over the already-computed pressure, velocity,
  and area values. A negative rounded reduction is mapped to neutral zero only
  when that interval straddles zero; an interval wholly below zero retains the
  negative value. This enclosure does not bound BEM solve, quadrature,
  discretization, or model error. Measured: plain CBIE's bound
  spikes ~21x over its own off-resonance value at the DISCRETE first
  interior resonance (ka ~ 3.20 on the 320-panel sphere — shifted from
  the continuum pi by the flat-panel discretization) while
  Burton-Miller's moves ~1.4x.
- `helmholtz_casebook` binary — the e2e casebook lane: mesh ->
  Burton-Miller solves -> impedance matrix -> spherical-harmonic table
  -> radiation efficiency, one JSON line per stage with FNV-1a-64
  result hashes. Every field is deterministic except `elapsed_ms`
  (wall-clock run evidence); `estimated_dense_bytes` is a deterministic
  3 x 16 n^2 model, not a measured RSS.
- `baffled_piston_face_impedance` — mean face `p/v` in an infinite baffle,
  converted by area from the sole `fs_phs::baffled_piston_impedance` owner
  (Bessel–Struve closed form, acoustic `p/U`). The former bem, phs and
  thin_plate ring copies were algebraically one algorithm; one test-only
  historical rule remains. G3 `historical_ring_and_overlap_integral_match_piston`
  records its 8/16 coarse-error reversal, requires decreasing errors on
  16/32/64 rings and final complex relative error < 0.2%, and reports the
  late observed order via fs-mms. Equivalent-disc singular near-field
  replacement allows O(h) error and coarse cancellation; the ladder does
  not establish an error floor or a smooth-quadrature order. An independent
  disk-overlap reduction of the Rayleigh integral uses smooth composite
  Simpson integration on 64/128/256/512 intervals, requires fourth order
  (3.8–4.2) and final relative error < 1e-10. Both acceptance bands reject
  explicitly planted 1% errors. Revised native results remain batch-pending.
  No unbaffled-pipe or nonuniform-face equivalence is claimed.
- `HelmholtzError` — stable `FS-BEM-HELM-*` refusals: bad parameter,
  too-coarse (< 6 panels/wavelength), dense work cap (8192 panels),
  shape mismatch, singular; optional geometry preparation additionally reports
  cache-byte budget, allocation, and cancellation refusals.

Wavelength resolution uses the longest retained triangle edge, including on
graded and highly anisotropic meshes. A tiny panel area cannot admit a long
edge spanning an unresolved acoustic phase. `SpherePanels::resolution_length`
and `OpenPanels::resolution_length` cache this geometric scale once; dense,
matrix-free, complex/molecular, penetrable-fluid and reciprocal closed/open
surface paths share it. Complex waves retain the existing conservative `|k|`
criterion, resolving attenuation as well as phase. Centroid-only surfaces have
no recoverable aspect ratio and retain `sqrt(maximum area)` as an explicitly
shape-regular estimate. Six panels per wavelength remains a necessary floor,
not a bound on solution error or a replacement for refinement.

### `helmholtz::iterative` (IRP M2.5-W5c)

The optional `IterativeOptions::directional` (M2.5-W5b) selects the directional
FmmC64 constructors through this same operator/state API. Its explicit order
and leaf capacity become the effective values, and all directional parameters
participate in checkpoint cross-wiring checks. `None` preserves the original
W5a route and kD>10 refusal. Dynamic layers and the RHS use the directional
owner; zero-k closures reuse its low-frequency W5a path, retaining the same
centroid quadrature. Each route keeps its own declared FMM approximation error.
`bem_directional_matches_dense_above_w5a_regime` exercises a slender authored
flared surface at k=8/m and root kD>10, with a three-order ladder, true GMRES
residual, and finest pressure error <=10*GMRES tolerance against dense LU.
Directional statistics travel in the existing resource receipt. This does not
claim a new full-band bell bake or relax resolution/solver/admission gates.

`IterativeRadiation::prepare/start/advance/finish` and
`solve_radiation_iterative` solve the existing constant-panel Burton-Miller
system through complex fs-fmm double/mixed matvecs and fs-solver `GmresState`.
The stacked real representation applies the complex operator once per matvec.
This path retains the surface cap of 81,920 panels and the resolution check,
without the dense LU cap of 8,192. Dense arithmetic and panel-bit goldens are
unchanged. Plain CBIE has its own operator on this owner; see below.

The actual dense Burton-Miller arm uses centroid off-diagonal kernels and
equivalent-disc self terms, not the triangle quadrature used by Plain CBIE.
The iterative path reuses those exact parent self terms and the same point
kernels in adjacent leaves; it introduces no second singular quadrature.
Static double/mixed row closures use the same effective order. W5a shares its
traversal with the dynamic operator; directional zero-k closures select W5a's
low-frequency partition. The static constant-density identities close under
the computed zero-k operator. Static-sum and cross-route FMM approximation
remain part of the operator error, independently tested against dense LU.

That centroid rule is first order. One-point weights near the collocation
point carry O(h) errors in S, D and N; under exact integration plain CBIE
converges at h^2 while Burton-Miller stays at h^1 (pexfp.111). The declared
`PanelRule::TriangleNear { max_pairs }` (`helmholtz_triangle_rule.rs`)
replaces every `triangle_is_near` pair as follows:
- S and D: exact static layers plus 8x8 remainders, as plain CBIE;
- N: the exact constant-density finite part,
  `(n_x.n)[k^2 S_T + SUM_e m_e.I_e] - SUM_e (n_x.m_e)(n.I_e)` with
  `I_e = INT_e grad G`;
- self terms: the flat triangle's;
- both static closures: recomputed under the same rule;
- D': keeps its centroid kernel and point closure. Exact flat-facet D' is
  itself first order, since `D'_0[1] + 1/2 = O(h)` on facets.

Routes and measurements:
- Dense: `solve_radiation_with_rule` / `solve_radiation_batch_with_rule`.
- Iterative: `IterativeOptions.panel_rule`, through the one `NearRows`
  mechanism; `impedance_sphere_triangle_rule_iterative_matches_dense` holds
  iterative to dense within 1.2e-6.
- `impedance_sphere_triangle_rule_ladder` measures dense orders (pulsating /
  scattering) of 2.07 / 2.03 at ka 1 and 1.91 / 2.01 at ka 3, against 1.4-1.7
  for the centroid rows.
- `PanelRule::Centroid` is the unchanged default, bit for bit.

TriangleNear's far pairs keep one-point weights. They leave an O(h) remnant
that alpha = i/k amplifies at low ka. `PanelRule::Triangle { max_pairs }`
keeps the same exact near zone and applies Radon's seven-point degree-5
rule to S, D and N on every other pair. D' keeps its centroid kernel.
- Dense, ka 0.05 dipole radiation resistance
  (`impedance_sphere_dipole_resistance_far_rule_ladder`). Relative error
  over rungs 2-4:
  - centroid: 33, 20, 10.6 (the resistance at 34x to 11.6x exact);
  - TriangleNear: 0.32, 2.24, 1.78 (stalls);
  - Triangle: 0.69, 0.37, 0.15.
- Dense orders under Triangle, pulsating / scattering: 2.02 / 2.01 at
  ka 1, and 2.01 / 2.02 at ka 3.
- Iterative (`BmCloudLayers`): the far field is far-field-only W5a passes
  over the Radon cloud that plain CBIE also uses. The near rows hold the
  dense rule's entry minus the node sums the pass covers.
  `impedance_sphere_triangle_far_cloud_matches_dense` holds it within
  1.3e-6 of the dense Triangle (tolerance 1e-5) at rungs 2-4, ka 0.05 and
  pi.
- Directional options are refused, because there is no far-only
  directional pass; that regime uses TriangleNear.

Centroid stays the default.

Plain CBIE (`Formulation::PlainCbie`, IRP .105) keeps the dense owner's
triangle quadrature; it is not replaced by centroid kernels.
- The far field is a far-field-only W5a double layer (single layer for the
  RHS) over a cloud of Radon's seven-point degree-5 rule on every triangle.
  Its first node is the centroid, so the collocation points are cloud points.
  `leaf_capacity` counts panels; the cloud capacity is seven times it.
- For every target, the near table stores the dense owner's `(S, D)` entry
  (`assemble_dense_with_static`'s own dispatch: Duffy self term, 8x8 near and
  4x4 far rules on the analytic static layers) for every panel with a node in
  the target leaf's direct-interaction leaves, minus the exact point sums of
  that panel's nodes the far pass covers. Entry plus far pass is the dense
  entry to FMM accuracy. A panel with every node separated enters through its
  seven-point sum alone.
- `CBIE_MIN_LEAF_RADII` (4) refuses leaves narrower than four panel radii,
  where the seven-point sum would stand in for the dense rule at more than
  1.1e-6 of the single layer (7.5e-8 at 6 radii, 1.3e-8 at 8).
- `NearFieldBlockJacobi` factors same-leaf blocks of the stored table, which
  are exact dense entries there. `max_near_entries` also charges the table,
  two complex entries per stored pair; `m2l_cache_bytes` goes to the one
  far-field double layer. Workers and the cache change no bit.
- Retained triangles are required. The directional route and Robin
  boundaries refuse: plain CBIE rows lie below the interior Dirichlet bound,
  where root kD is within W5a.

The near table is this owner's one sparse near-field mechanism, keyed by
`NearOperator` (S, D, N; D' always stays centroid). A crate-internal
`NearRows` rule plugs into Burton-Miller through `prepare_with_near`.
- Its rows correct the centroid operator, P2P near field included: A gains
  `D + alpha N`, B gains `S`, and the rule's `(dA_ii, dB_ii)` join both
  diagonals and the block-Jacobi blocks.
- Row lengths are known before allocation, and `max_near_entries` charges
  pairs times operators. Rows must be finite, strictly ascending and
  worker-independent, or the build refuses.
- `prepare` is `prepare_with_near(None)`, bitwise, and plain CBIE refuses a
  rule.
- `burton_miller_near_rows_enter_as_declared` plants a three-operator rule
  and checks both sides against the declared sums within 1e-14, a converged
  preconditioned solve, the inclusive budget, and the CBIE refusal.

`helmholtz::iterative::cbie_tests` checks every pair of a 1280-panel
icosphere at k = 1.3 against the assembled dense matrix: 354080 pairs with no
separated node are bitwise the dense entries; 48192 corrected pairs rebuild it
within 3.5e-16; 1236128 separated pairs' seven-point sums are within 4.2e-8
(single) and 3.6e-8 (double, times the panel radius) of the single layer (ts1,
6a478871b). `tests/iterative_cbie.rs` solves 5120-panel bodies at fmm order 7
and GMRES 1e-10 against dense plain CBIE: sphere 6.5e-7 (k = 0.05, depth 2),
8.3e-7 (k = 1) and 9.0e-7 (k = 2, depth 3), ellipsoid 1, 0.7, 0.5 m 7.6e-7
(k = 1), all within 1e-5, in 9 to 11 iterations. Its order ladder gives
1.8e-3, 4.4e-5 and 8.3e-7 at p = 3, 5 and 7. Workers 1, 5 and 16 and the row
cache are bitwise (ts1, 88e7dd212).

`IterativeOptions` explicitly declares GMRES tolerance, restart (1..256),
iteration cap, FMM order/leaf capacity, and near-field storage budget.
`Preconditioner::NearFieldBlockJacobi` factors exact point-panel blocks within
individual direct-interaction leaves, split at 1..128 panels in stable order.
Its total complex entry count is checked before block allocation. `None`
provides the unpreconditioned comparison. Right preconditioning preserves the
true Euclidean residual of the original FMM equation rather than replacing it
with a left-preconditioned norm. No SPAI, Calderon or transpose-solve claim is
made. FMM's existing kD<=10, translation, coefficient and near-pair caps remain
in force, with typed refusals and no dense fallback.

`IterativeOptions::near_field_cache_bytes` caches the double and mixed
operators' adjacent-leaf pair values once per prepared operator, half the
budget each. The default is 512 MiB.
- Without the cache, every matvec re-evaluates those values.
- The cache changes no bits.
- An operator whose cache exceeds its half stays uncached, and so does the
  directional route.
- `ResourceReceipt::near_cache_bytes` reports what was cached.
- `bem_near_field_cache_is_bitwise` checks the default against a zero budget
  and a too-small one.

`IterativeOptions::workers` (default 1) is set on every W5a engine the owner
builds (`FmmC64::set_workers`). `IterativeOptions::m2l_cache_bytes` (default
0, off) caches the double and mixed operators' M2L kernel rows, half the
budget each, like the near-field cache.
- Neither changes a bit.
- `declared_fmm_workers_leave_the_iterative_solve_bitwise` pins the pressures,
  the iteration count and the residual claim across workers and the row cache.
- On ts1, a 5120-panel solve takes 173.4 s on one uncached worker and 15.7 s
  on 16 workers with the row cache.

Checkpoints clone the `GmresState` only at completed restart boundaries.
Cancellation is polled before each matvec and between its two synchronous FMM
passes. The current FMM pass drains before cancellation returns; no mid-pass
latency claim is made. A latched error suppresses further operator calls while
the current Arnoldi cycle unwinds; that cycle is discarded and the preceding
checkpoint remains resumable. A deterministic matvec-index predicate permits
race-free cancellation tests. Checkpoints reject a different prepared geometry
or numerical configuration. An exhausted solve returns `NotConverged` with
last pressure, full report and receipt; it never publishes an unconverged
`RadiationSolution`.

`ResourceReceipt` records preparation/RHS/solve wall time, actual matvec and
iteration counts, FMM order/stats, stored block/indices bytes, and Linux
process-wide VmHWM (None on unsupported hosts). These measurements are not
portable speed or memory guarantees. The shared solution's
`condition_lower_bound=1` is the universal bound, with no dense probe claimed;
`dense_cap_utilization=N/8192` can exceed one. Radiated power reuses the parent
reduction and its existing roundoff-only enclosure.

`tests/iterative.rs` provides G3 dense comparisons through s=4 (5,120 panels)
and an authored lathed bell, including ka=pi; a live FMM-order ladder; G2 s=6
(81,920 panels), ka<=5, orders 2/3/4 with residual <=1e-6 and finest-order
analytic mean-pressure error <=1e-3; G1 four-rung fixed-panels/wavelength
preconditioned/unpreconditioned counts; G0/G5 bitwise restart replay; G4
cancellation; and storage/input/panel/iteration refusal boundaries. JSONL
records numerical and process evidence. Coarse order rows report observed
errors even when they miss the final qualification band. All numerical and
resource qualifications remain pending central native verification until the
bead's run receipt exists. In particular a small Krylov residual does not bound
FMM, panel or model error, and order 4 is not an automatic accuracy certificate.

### `helmholtz::aperture` (IRP M2.5-W3a)

`aperture_impedance(surface, apertures, omegas, medium, route, gate)` returns
one row per frequency of the self-plus-mutual aperture matrix
`Z[a][b] = SUM_i A_i shape_a(i) p_b(i)` [N s/m] under `exp(-i w t)`: the
generalized force on aperture `a` per unit generalized velocity of aperture
`b`. The body force is `-Z v`, and a mass load has a negative imaginary part.
An `Aperture` is a real per-panel outward velocity shape plus a declared
reference point and unit normal. The `declaration` hash binds the surface,
shapes, reference planes, band, medium and route.

There are two routes, and neither is symmetrized, clipped or loaded:
- `Collocation` uses Burton-Miller. It solves one dense batch up to
  `MAX_DENSE_PANELS`. Above that cap it uses the FMM/GMRES owner when
  `IterativeOptions` are given, and refuses otherwise.
- `Galerkin` uses the reciprocal P1/P0 velocity ports, within their restricted
  band.

Each row reports two measured quantities:
- `reciprocity_defect = ||Z - Z^T||_F / ||Z||_F`;
- `passivity_margin`, the smallest eigenvalue of `(Z + Z^H)/2` over
  `||Z||_F`.

`volume_gains[a] = SUM_i A_i shape_a(i)` [m³ per m/s] converts generalized
velocity to volume flow. `acoustic(row)` returns the acoustic matrix
`Z[a][b] / (q_a q_b)` [Pa s/m³] and refuses a zero gain. For indicator shapes
it maps volume flows to aperture-mean pressures. This is the coupled load that
fs-duct `CoupledTermination` takes.

`baffled_piston_array_impedance` covers coplanar baffled circular pistons.
- Self terms are the closed-form face impedance times the piston area.
- Mutual terms are the Rayleigh double integral on a polar Gauss-Legendre ×
  midpoint disc ladder. A term is accepted at the first two rungs that agree
  within `PISTON_MUTUAL_TOLERANCE = 1e-9`. Overlapping, touching or unresolved
  pistons are refused.

The tests check the following:
- **G2, two pistons:** the mutual term matches the exact King
  (Hankel-transform) representation, with its tail bound, within 1e-4 for
  kd in [0.5, 20] and ka <= 2. Measured agreement is 1e-12 to 1e-10.
  Pritchard's series itself is not evaluated.
- **G0, reciprocity and passivity:** pistons and Galerkin meet both at
  roundoff. The partitioned-sphere sum is compared to the pulsating sphere
  for both routes.
- **G1, ladder:** a closed equal-area-faceted tube (two wall holes plus the
  mouth) runs on a 4-rung ladder of 360 / 1488 / 5856 / 23232 panels. The
  fourth rung routes past the dense cap.
  - It converges monotonically at a consistent observed order of about
    h^1.25: successive Frobenius steps 2.1e-2 / 8.6e-3 / 3.6e-3, ratios 2.4
    (trj b187). The piston velocity jumps at every aperture edge.
  - The bead's 1e-3 G1 target is NOT met at these rungs.
  - Collocation and Galerkin also converge together on a partitioned sphere.
- **Duct-network wiring:** `aperture_matrix_terminates_a_duct_network`
  terminates the coarse tube's plane-wave bore at its two holes and mouth
  through one coupled load. It checks `Re Z_in >= 0` on a 121-point
  300-3000 Hz scan. It also checks that the mutual terms move at least one
  golden-refined input-impedance peak by more than 0.01 cent against the
  self-only diagonal. Measured on the 360-panel matrix: -2.33 cent at
  1040 Hz and +2.70 cent at 2109 Hz (yto y2). This is a wiring check, not an
  accuracy claim.

No-claim: aperture shapes are panel unions. The collocation route is
reciprocal only to discretization error. A passive row is not a causal
positive-real time model; that is M2.5-W3b's realization gate.

### `plate_aperture` (IRP pexfp.366)

`plate_aperture_conductance(raster, options)` returns the potential-flow
conductance `K = Q / delta-phi` [m] of a raster of open square cells cut
straight through a rigid plate `0 <= z <= t` in an infinite baffle.
Rayleigh's thin circular hole gives `K = 2a`. Each raster cell is split
`m x m`.
- Each half-space is exact through the Neumann kernel `1/(2 pi r)`, by P0
  collocation with exact cell integrals, applied as an FFT convolution.
- The duct is cell-centred finite volumes on the voxel prism, in
  `round(t / h)` layers.
- The coupled system `(I + Lambda H2) u = Lambda e` is solved by CG in the
  `H2` inner product. A thin plate solves `2 H w = 1`.
- Refusals: invalid inputs, the cell, FFT-grid and voxel budgets, and a CG
  miss, including an outer residual above ten times the tolerance once
  recomputed.
- Results are bitwise identical at every worker count and in any listed
  cell order.

`K` converges at about first order in `h` (thin plate 0.94-0.95, thick
1.2), so a caller ladders `m` and extrapolates. Gates (tests/plate_aperture.rs):

| case | measured | gate |
|---|---|---|
| Rayleigh disc, cells of `R/rpx`, `rpx` 4 to 64 | error 0.18 % at 64 | order >= 0.8, finest <= 0.25 % |
| 42 x 3.1 mm ellipse at the CT's 0.5 mm, `m` 1 to 8 | limit +0.40 % of `2a AGM(1, b/a)` | within 0.5 % |
| flanged duct, `t = 4a` | end correction 0.8189 vs 0.8216 | within 1 % |

- Rotations, a mirror and a translation of the raster agree to 1e-9.
- The 3-D finite-volume Schur complement agrees with the duct map.

`cannone_f_hole_conductance_ladder` measures the CT f-holes of il Cannone
as a thin plate, each hole alone and both together, against the sum of
their equivalent ellipses, and gates only convergence.
`cannone_f_hole_conductance_table_reproduces` solves the specimen's
`top-plate-fhole-conductance.tsv` again (`examples/fhole_conductance.rs`:
both holes, and each alone, through the measured pair bore and through its
own bore) and requires every printed digit. The h8 body table's neck must
carry the pair row.

No-claim: the plate is flat and infinite, and each side is a half-space,
not the arched top, the finite body or its closed cavity. The flow is
inviscid and incompressible. A raster stands in for the outline, and on a
slender ellipse its staircase raises `K` by 0.4 %.

## Invariants

### Solid-mode radiation loading (IRP M2.5-W4a)

`modal_loading` consolidates the original Helmholtz source batches and complete
power-consistent Rayleigh loading. Closed bodies require retained outward closed
triangle geometry through the existing `near_field::Geometry` admission. Planar
cells carry physical square metres and every signed mode trace. Both routes
retain the work-conjugate projection `Z_rs = integral phi_r p_s dS`; no structural
mode or cross term is discarded, and a BEM basis above 256 fields is partitioned
deterministically into complete batches. The caller declares the output-entry
budget for the full frequency sweep and a provenance hash. The bake declaration
hash includes geometry, all traces, route, medium, frequencies and budget.

For `k * source_diameter <= 0.5`, power-consistent Rayleigh resistance uses
compensated Cartesian source moments and a positive angular Gram rule. This is
the same point-area sinc operator, evaluated without subtracting almost equal
monopole powers for a weak dipole, quadrupole or higher cancelling field.
Degree-16 outgoing amplitudes have an absolute exponential remainder below
`2.2e-20 * sum(abs(area * shape))`; the 17-point polar / 34-point azimuth rule
integrates their degree-32 products exactly in exact arithmetic. Cartesian
moments are projected once into compensated direction/degree coefficients;
each further frequency contracts 17 coefficients instead of 153 monomials per
direction and field. This geometry cache uses `153 * 17 * fields` doubles
(2.54 MiB at the 128-field cap), replacing the Cartesian moment cache. Angular
nodes are also reused. Both geometry projections and recurring frequency work
are included in the existing modal ceiling. The equivalent-disc reactance
remains bit-identical.
G1 tests compare orders 1, 2 and 4 with independent positive sine-power angular
integrals down to 0.001 Hz; G3 compares the overlapping sinc branch and translated,
rotated and reordered sources. No eigenvalue clipping or loss floor is applied.
This absolute approximation is not a relative guarantee for arbitrary dark modes,
nor does it cure cancellation in a downstream ill-conditioned matrix quadratic form.

`PreparedSurface::coherent_radiated_power` supplies a direct outgoing-power
observation that avoids that downstream matrix cancellation. For compact sources,
it contracts complex modal velocities into each cached directional coefficient
BEFORE summing Taylor degrees and squaring. Larger sources combine modal velocities
at each cell and integrate the resulting field with positive angular weights.
The opposite azimuths are paired as `|even|² + |odd|²`; simply halving the angular
domain would be wrong for complex velocities. Double-double contractions retain
weak residuals, without negative-power clipping or an artificial loss floor.

The noncompact rule selects polar order n and 2n azimuth points from
`epsilon = (kD)^n/n! <= 1e-18`, capped at n=128, and refuses if unresolved or over
the existing scalar-work ceiling. In exact arithmetic it integrates products of
degree-(n-1) plane-wave polynomials exactly. With
`A = sum_cells,modes |area*shape| |velocity|` and
`s = rho*omega²/(4*pi*c)`, the absolute angular truncation bound is
`s*A²*(2*epsilon+epsilon²)` for the compact polynomial, twice this for the
noncompact exact-exponential rule. The reported bound excludes floating-point,
geometry, finite-cell quadrature and material errors; it is not a relative
accuracy or passivity certificate for arbitrarily dark sources. A single cell
uses its exact point-source limit. The observation retains the six-cell-per-
wavelength admission and does not alter the loading or time realization.
Cancellation is polled at entry, during source contraction and each angular
direction; the existing compact geometry-cache initialization remains bounded
but noninterruptible. G1/G3 regressions include separate-coordinate order-1/2/4
cancellation down to 0.001 Hz, complex directional phases, angular-branch
agreement, 128 retained modes, source/basis/frame changes and exact nulls.

Rows use `exp(-i omega t)`: resistance is `Re Z`, added mass is `-Im Z/omega`.
The unchanged matrix must have relative reciprocity defect at most
`DENSE_RECIPROCITY_BOUND` = 1e-12 (an iterative row's bound is below) and
minimum resistance eigenvalue at least `-1e-12 ||Z||_F`, with a resolved full
eigendecomposition residual. Failed numerical data refuse instead of being
symmetrized, eigenvalue-clipped or supplied with extra damping. This is a
sampled floating-point diagnostic, not an interval or positive-real certificate.
Frequency lookup admits exact declared samples only; it does not extrapolate.
Cancellation is checked before geometry work, each frequency/source batch, and
publication. Original dense factorization still has no mid-call cancellation.

`tests/modal_loading.rs` contains G0 passivity and cap boundaries, G5 replay,
257-field batch equivalence, an independent Rayleigh double sum of the declared
point/disc quadrature, a three-rung six-mode plate-efficiency ladder against an
independent angular sine-transform integral (2%), and the strict pulsating and
translating sphere closed forms (1e-3). The sphere reference is Sutherland
WR68-12, NASA 19690015869, Appendix B equations B-9/B-18 and Fig.8, re-derived
under this time convention. The Wallace paper (JASA 51, 1972, 946-952,
doi:10.1121/1.1912943) full-text figure/equation digits were not available; the
independent integral is not claimed as a digitized-paper comparison.

Central verification found negative resistance at `ka=0.05` on all three
sphere meshes (2208/3968/6240 panels). The assembly now also enforces the
static double-layer/weighted-adjoint compatibility described below; native
verification of that correction is pending. The strict sphere tests retain
their passivity refusal and 1e-3 threshold. Sampled PSD does not
provide the M2.5-W3b causal passive MIMO realization, a solid-mesh-to-acoustic
mortar transfer, or the required source-backed free-aluminium-bar loaded-pole
comparison. The existing vibroacoustic casebook consumes planar bake rows and
checks loaded frequency-response force/power balance; it does not satisfy the
bar's 0.1-cent frequency/2% decay requirement. Full M2.5-W4a remains open for
these numerical and composition gates. No piano example or panel-bits golden
was changed.

`modal_radiation_loading_routed` is the M2.5-W5 route past the dense cap.
- At or below `MAX_DENSE_PANELS` it is the dense bake under the options'
  declared panel rule; the other options are unused. Under the centroid rule
  (or no options) that is the dense bake bitwise, declaration identity
  included. A triangle rule is `solve_radiation_batch_with_rule`'s bake, with
  the rule bound into the declaration. Plain-CBIE rows ignore every rule.
- Above it, a closed body's rows go to `IterativeRadiation` with the declared
  options, each in its own formulation (the route's or `GeometryPolicy`'s).
  Every field must converge to the GMRES tolerance or the bake refuses with
  `Iterative`. The declaration identity binds the numerical options, not
  workers, caches or the panel rule's pair budget; the rule itself is bound.
  Closure and orientation are admitted by
  `near_field::components`, the owner `near_field::Geometry` uses below the
  cap; that receiver admission itself stops at the dense cap.
  `GeometryPolicy::new` admits surfaces up to `MAX_SURFACE_PANELS`, since its
  box bound holds at any size; its `solve_batch` keeps the dense owner's cap.
- Reciprocity of an iterative row (the .105 G2 ruling, 2026-10-23):
  - FMM truncation is not reciprocal, and the GMRES residual cannot see it:
    the residual is measured against the same truncated operator.
  - On a unit icosphere at ka 0.05 (plain CBIE), a 1e4-tighter GMRES
    tolerance moves the defect 2%, while +2 FMM orders cut it 30-67x (ts1,
    2026-10-23).
  - So each row also applies order p + 2 and p + 4 twins of its operator to
    the converged fields, one operator at a time, under the declared budgets.
    eps_fmm is the largest `||A_{p+2} x - A_p x|| / ||A_p x||`, and
    eps_fmm_next is the same one step up.
  - The bound is `max(10 eps_gmres, eps_fmm)`, where eps_gmres is the largest
    achieved Euclidean GMRES residual.
  - The p + 2 estimate is Richardson-style, valid while the FMM error decays
    geometrically. A row with `eps_fmm_next > 0.3 eps_fmm` refuses with
    `UnresolvedFmmEstimate`.
  - fs-fmm admits orders up to 12, so a row declared above order 8 refuses.
    The twins obey the row's declared budgets, including
    `IterativeOptions::max_coefficient_slots` (default fs-fmm's v1
    `MAX_COMPLEX_COEFFICIENT_SLOTS`), a budget that changes no bits. Twins at
    p + 2 and p + 4 on the 20480- and 81224-panel spheres need more than v1.
  - The row's receipt, `IterativeErrorScale`, records eps_gmres, the three
    orders, both estimates and each twin's budget, applications and wall time
    (`TwinWork`). `reciprocity_bound` records the bound. Nothing is
    symmetrised, and the twins change no bit of Z.
- The passivity gate is unchanged. `NonReciprocal` (with its bound and scale)
  and `NonPassive` carry the refused matrix as a diagnostic, never as a row.
- `iterative_rows_are_reciprocal_to_their_measured_error_scale` (unit test)
  forces 1280- and 20480-panel icospheres onto the iterative owner at FMM
  orders 4/6/8 and 6/8:
  - every row passes its bound;
  - the defect and eps_fmm each shrink at least 10x per +2 orders;
  - the same row with a swapped projection refuses by more than 3 decades;
  - an unresolved estimate refuses;
  - a direct row keeps 1e-12.
- `routed_modal_loading_is_the_dense_bake_below_the_cap` checks:
  - the bitwise dense equality for a plain-CBIE and a Burton-Miller row;
  - the unrouted W5c refusal;
  - a live closure refusal;
  - a worker-bitwise, option-bound declaration at 20480 panels;
  - that declaration's two-mode row admitted with its full receipt;
  - the routed rows' Z, printed bit for bit, equal to the same bakes before
    the twins (ts1, 2026-10-23).

G3, loaded modes of a real free aluminium bar (`fs-couple/tests/loaded_free_bar.rs`):
- Fixture: beam 40 of NASA TM X-2852 (specimen `aluminum-beams-howell1974`)
  with the 2024-T3 matdb constants, a declared temper substitution, on the
  TET10 solid. Faces of the same Kuhn grid carry the vertex average of u.n,
  and `GeometryPolicy` selects plain CBIE.
- `loaded_modes_match_coupled_solve` compares two routes for every elastic
  mode below 1.6 kHz, on [12,2,1] and [16,4,1] grids:
  - one routed single-mode bake per mode, which is reciprocal by
    construction;
  - the coupled FE+BEM system in the complete in-vacuo basis of the same FE
    model, rigid modes included.
  They agree within 4.1e-4 cents and 1.5% in decay (gates 0.1 cent and 2%).
  The decay difference is the fluid coupling to other modes that a
  single-mode bake leaves out. Multi-mode bakes stay refused: the complete
  matrix's collocation reciprocity defect is 1e-2, against the 1e-12 gate.
- `loaded_free_bar_radiation_decay_ladder` puts the fundamental on a fixed
  [96,16,2] TET10 structure and nests BEM panels at 496/1760/6592. In vacuo
  it gives 371.4 Hz against the measured 369 Hz; air lowers it by 0.80 Hz.
  The radiation decay rate 5.69/5.94/6.05e-5 1/s changes by 4.1% then 1.8%
  (order about 1.2). Its log decrement 1.6e-7 is 5e-4 of the measured
  air-minus-vacuum 3.2e-4, which also holds viscous loss: a one-sided bound
  (ts1, a2e2fb443).

1. G0 Gauss identity: the assembled Neumann operator applied to ones
   gives −1 at every centroid within discretization tolerance
   (bem-001) — sign conventions cannot drift silently.
2. Sphere analytic (G2): mean surface-speed error vs 1.5·U·sinθ
   < 0.03 at 1280 panels and decreasing under refinement (bem-002).
3. The FMM path IS the dense operator: matvec and transpose relative
   deviations are < 1e-4 at order 6; GMRES(FMM) reproduces the
   dense-LU solution to < 1e-3 with iterations ledgered (bem-003).
4. Hess–Smith: lift slope within 5% of the thickness-corrected
   2π(1+0.77t) and above thin-airfoil 2π; stagnation Cp = 1 within 5%;
   Kutta row satisfied to roundoff; adjoint dCl/dα matches central FD
   to 1e-6 (bem-004).
5. Free wake: Wagner-like start (first/steady in [0.3, 0.7]),
   asymptote within [0.9, 1.05] of the pressure-derived screening
   circulation scale,
   coarse-grained monotone growth (early lumped-starting-vortex dips
   are ledgered, not hidden), Kelvin circulation bookkeeping, bounded
   stable roll-up, and bitwise determinism of the complete wake/history/trace
   state (bem-005).
6. NACA 0012 (G2): a least-squares lift slope over `|alpha| <= 8 degrees`
   remains from zero to 20% above the independent NASA TM-4074 table-I slope
   at Mach 0.15, Reynolds number 5.97 million, and free transition; odd
   symmetry is retained and the validation API refuses `|alpha| > 10 degrees`
   (bem-006). The one-sided band records the report's observed inviscid-theory
   overprediction and is an honesty envelope, not viscous parity.

## Error model

Public constructors and numerical entry points return typed `BemError` values
for malformed/non-finite geometry, mismatched vectors, singular systems,
invalid tolerances, zero trace stride, and explicit dense/FMM/transient work
envelopes. `solve_exterior` never publishes an unconverged iterate as ordinary
success. Airfoil, sphere-panel, and wake state storage is read-only after
validated construction. Physical honesty: every battery verdict carries the
`inviscid-screening` model label; no viscous claims anywhere.
The NACA 0012 validation boundary is checked before geometry allocation and is
inclusive at ten degrees. The underlying generic `panel2d::solve` stays
available for explicitly inviscid mathematical screening outside that evidence
envelope.
`BemError::AllocationFailed` covers explicitly reserved BEM geometry, dense,
wake, and exactly sized trace buffers. The separately documented process-level
allocator no-claim still applies inside fs-fmm passes, and fs-solver's current
GMRES state allocation remains infallible after the bounded BEM admission step.

## Determinism class

Bit-deterministic across runs on a platform (dense LU, deterministic
FMM underneath, fixed shedding/convection order).

## Cancellation behavior

Wake state is cloneable and callers can chunk at fallible `step` boundaries.
`PreparedCbieGeometry::new_with_cancel` additionally polls a caller callback
before allocation and between rows of its frequency-independent static table;
no partially prepared object is returned on cancellation.
Dense panel assembly/LU and each FMM/GMRES call do **not** currently accept a
`Cx`, poll cancellation, or expose mid-call resume state. Cross-crate Cx/resume
integration is tracked separately under `frankensim-ccmn`; no cancellation
latency claim is made here.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`src/panel3d.rs` unit tests: the private `LinearOp::apply_transpose`
wrapper matches the dense transpose, and invalid `SpherePanels` vector
shapes are rejected before FMM math. `tests/battery.rs`:
bem-001 Gauss identity; bem-002 sphere analytic; bem-003 FMM-vs-dense
matvec, transpose + GMRES; bem-004 Hess–Smith slope band, Cp sanity,
Kutta, adjoint gate; bem-005 impulsive-start free wake; bem-006 NACA 0012
pre-stall lift slope against NASA TM-4074 table I plus the ten-degree refusal;
invalid-input/work/trace refusal; unconverged exterior-solve refusal with
retained report. The source report is pinned by URL and SHA-256 in the battery;
NASA marks it as U.S. Government work with public use permitted.

### Helmholtz invariants

1. Kernel formulas are pinned by central-finite-difference tests of G.
   Plain-CBIE triangle surfaces use exact static single/double-layer integrals
   for every interaction, preserving `D_0[1] = -1/2` at a panel centroid.
   Fixed Gauss integration covers the regular Helmholtz remainders: 8x8 near
   a panel and 4x4 for far interactions. A far-panel centroid approximation
   loses spatial moments that can dominate weak bending radiation; replacing
   only its static part is insufficient. The self and near-face single layers
   share the same static integral. Independent
   rectangular-panel single-layer and solid-angle formulas test thin opposing
   faces down to a 1 micrometre gap and both orientations. Gauss-Duffy
   integration handles the self remainders; scale covariance and the analytic
   pulsating sphere pin this path. A closed 1 mm plate boundary additionally
   tests the static row identity and positive radiation from sixteen bending
   fields at 40 and 1640 Hz, without a power floor or relaxed power gate.
   Far remainders now require sixteen kernel evaluations per pair; dense
   factorization and batch counts are unchanged. These collocation tests are
   not a general passivity or convergence certificate. Burton-Miller retains
   the established equivalent-disc self and
   centroid off-diagonal scheme because its hypersingular finite-part
   operator needs a different Galerkin-quality rule.
2. The hypersingular static self entry uses the exact closed-surface
   identity `N_0[1] = 0` as a discrete row sum, so the same point-panel
   quadrature appears on both sides and its error cancels. Burton-Miller also
   enforces `(D_0+I/2)1=0` by its static self coefficient, shared with `D'_0`
   to retain the area-weighted adjoint. Consequently `A^T(D'_0+I/2)=0`
   for the panel-area vector `A`, compatible with `A^T N_0=0`.
   Leaving the double-layer self zero while closing only `N_0` violates this
   compatibility and introduces a net static flux amplified by `|alpha| = 1/k`.
   No pressure, resistance, spectral clipping or empirical parameter enters
   the correction. G0 checks the actual assembled static blocks on a sphere
   and an ellipsoid; the pulsating-sphere test now asserts the same 8%/positive
   power envelope at `ka=0.05` and `0.2` as at the higher frequencies.
3. Pulsating-sphere impedance within the authored per-arm envelopes
   above; radiated power positive for mesh-resolved velocity fields;
   area-weighted impedance-matrix reciprocity to 0.5% measured.
4. Repeat solves are bitwise identical, including batched solves, their
   equivalence to independent single-field solves, the condition diagnostic,
   and spherical-harmonic coefficient tables.
5. The spherical-harmonic basis is self-checked: Gauss–Legendre weight
   and moment identities to 1e-13 and normalized associated-Legendre
   orthonormality to 1e-12 under the same quadrature the projection
   uses. The phase convention is pinned by a genuinely non-axisymmetric
   falsifier: an x-axis dipole must satisfy `a_(1,-1) = -a_(1,1)`
   (measured defect 4e-16) and reconstruct at held-out directions —
   a one-sided Condon–Shortley or conjugation slip cannot survive it.
6. The condition diagnostic discriminates the resonant arm: CBIE's
   bound inflates across the fictitious-frequency band while
   Burton-Miller's stays flat (measured 20.8x vs 1.4x; asserted > 2x
   separation with the peak found by a coarse-then-refined scan).
7. Pulsating-sphere finite pressure matches within 8% at `ka = 1` and
   converges to the direct far field (2% remainder).
8. Integrated observer regressions compare the triangle phase to an analytic
   rectangular-aperture Fourier integral and the default percussion finite
   microphone to independently adaptive Green rows for a thin-plate bending
   field. The same trace exposes the material difference from the preserved
   centroid observers. Fixed observer quadrature is not an error certificate.

## No-claim boundaries

- 3D LIFTING surfaces (Kutta strips, wake SHEETS) and the fs-vpm
  pairing for flapping gaits — the 2D shedding loop ships; 3D is the
  flagship successor.
- Closed-form Helmholtz panel-integral far fields. Explicit integrated
  observers use fixed Gaussian rules; legacy observers retain centroid rules.
  Neither constitutes a bound on pressure error for arbitrary boundary traces.
- Induced-drag decomposition and force/moment beyond lift (Cp
  machinery exists; the Trefftz-plane analysis is successor scope).
- Elastostatic BEM (staged later per the bead, noted not promised).
- XFOIL-class viscous corrections (never claimed — screening only).
- Post-stall or separation behavior for NACA 0012. `bem-006` validates only the
  pre-stall lift-slope envelope; it makes no drag, transition, boundary-layer,
  maximum-lift, or stall-onset claim.
- FMM-accelerated 2D wake convection. The shipped path is a direct all-pairs
  screening kernel with an explicit 1,024-vortex / 1,048,576-pair per-step
  admission ceiling.
- Helmholtz (`helmholtz` module): the original BurtonMiller centroid scheme
  produced negative low-frequency resistance (measured -9.4 vs +1.03 Pa s/m
  at ka = 0.05 on 320 panels). The static compatibility correction above
  has not yet passed native verification; no expanded low-frequency accuracy
  claim follows from its algebraic identity. General singular finite-part
  triangle quadrature and the strict modal sphere ladder remain open gates.
  Passivity is claimed only for mesh-resolved velocity fields: a
  quadrupole's true radiated power at ka = 1 sits below the noise floor
  on the 80-panel fixture (measured, recorded). `condition_lower_bound`
  is a LOWER bound from six deterministic probe solves: a large value
  is a reliable warning, but a small value certifies nothing (a rigorous
  Hager/Higham estimator needs adjoint solves and is the recorded
  follow-up). `DirectivityTable` truncation is estimated via
  `captured_fraction` only after independent angular quadrature refinement.
  Finite grid agreement is not a rigorous bound for an arbitrary evaluator
  and cannot exclude adversarial unresolved angular structure. The table inherits every
  accuracy boundary of the underlying solve arm. No FMM acceleration
  (dense cap 8192 panels), no scattering/incident-field path,
  no half-space or impedance boundary conditions, no Bessel-backed
  piston closed form (small-ka series only until the duct bead's special
  functions land). Finite-point admission is a non-certifying discrete
  solid-angle guard. The legacy finite observer retains centroid-panel error;
  its explicit integrated variant still has fixed-quadrature and boundary
  discretization error. Use the separate `near_field` owner when adaptive
  integration and retained quadrature diagnostics are required.

## radiation_bake (bead zolja)

- `bake_radiation_load(surface, driven, omegas, medium, l_max,
  source_id) -> RadiationLoadBake`: an OFFLINE frequency sweep over
  `solve_radiation` producing an area-averaged `Z_L(omega)` table
  (specific and acoustic normalizations both retained) plus
  per-frequency SH directivity for a driven mouth on a closed exterior
  mesh. Formulation policy: plain CBIE below mouth `ka = 0.5` (the
  recorded Burton-Miller low-`ka` resistance boundary), Burton-Miller
  above. Every row keeps the guardrail receipts
  (panels-per-wavelength, condition lower bound, radiated power,
  passivity margin); the sweep REFUSES past the mesh's
  panels-per-wavelength bound rather than extrapolating.
- `lathe_profile(profile, circumferential) -> LathedBell`: a closed
  bell fixture (flare wall + throat cap rigid, mouth cap = the driven
  piston, panel order suffix-masked). Caps are subdivided into
  CONCENTRIC ANNULI at the wall panel scale — an executed lesson: a
  single-fan 60 mm mouth cap made every Burton-Miller row's `Re Z`
  negative by a near-constant offset (the centroid-quadrature
  resistance artifact scales with panel size); annular caps restored
  positive, smooth margins across the band.
- Conformance: zb-001 (the bake reproduces the pulsating sphere inside
  the solver's own tested bands: CBIE < 4%, Burton-Miller < 8%),
  zb-002 (sweep refusals), zb-003 (divergence-theorem closure and
  outward orientation of the lathe against the frustum-stack volume).
- No-claim: the bake certifies nothing beyond the solver's own bands;
  table consumers own interpolation and out-of-support refusal
  (fs-duct `TabulatedLoad`); no rational realization is performed here
  (the fs-vfit certificate-class lane belongs to the runtime consumer).

### Declaration-bound scalar refinement (IRP M2.5-W1)

`radiation_bake::declared::BakeDeclaration` validates and hashes source geometry
identity, the full requested frequency grid and band endpoints, relative
operator tolerance, physical reference surface, mode content, SH degree and
medium before solving. Private fields prevent later mutation under the same
identity. Source and mesh identities are distinct. A `PanelSource` regenerates
each `PanelLadder` rung directly from its immutable source geometry; the trait's
source identity/count/geometry contract is trusted, not independently certified.
Generated driven panels must lie on the declared mouth plane and face its
outward normal. Whole-surface body modes explicitly refuse duct-load export.

`bake_declared` calls the unchanged `bake_radiation_load` owner at every
frequency/rung. Complete requires at least three solved rungs and two successive
complex impedance and SH-vector relative L2 changes within the declared
tolerance. This measures the full retained complex SH vector through the
declared degree; it is not a relative guarantee for every quiet coefficient or
a bound on omitted SH degrees. It is finite-frequency-grid refinement evidence,
not a uniform-in-frequency interpolation error enclosure or physical validation.
Negative resistance beyond the same relative tolerance refuses. No panel
arithmetic, formulation threshold or legacy bake result changes here.
`bake_declared` stays dense: a rung past `MAX_DENSE_PANELS` is
`Blocked { route: "M2.5-W5" }`. `bake_declared_routed(.., Some(&options), ..)`
is that route: such a rung is solved by `bake_radiation_load_iterative`
(the matrix-free FMM + GMRES owner, under the caller's explicit tolerance,
FMM order and storage budgets), while rungs within the cap stay dense. Each
`RefinementRow` and receipt line records `route` ("dense" or "iterative").
An iterative refusal is `BakeError::IterativeRoute { at_hz, what }`, and an
unconverged solve is never a row. The same ladder, tolerance and three-rung
rule decide completeness.

The exact next source panel count is checked before allocating its surface or
dense matrix. More than 8192 panels returns `Blocked` with frequency, needed
panels and route `M2.5-W5`. If three solved rungs are impossible within the
declared ladder's dense budget, it blocks before wasting earlier solves.
Exhausted below-cap ladders return `NotConverged` with retained rung diagnostics.
Neither status exports a full-band load. `partial_band` explicitly exposes only
the contiguous proven prefix, with a `band_refusal` and the declaration hash;
an expected-hash mismatch refuses. Consumer tables must keep that exact support
and their `OutOfTable` behavior. Rung receipts retain source area, estimated
minimum panels, actual panels/wavelength, changes, status and measured time;
actual runner host/command are supplied when rendering receipts.

Cancellation is polled before/between existing solves. The inherited dense LU
still has no mid-call cancellation latency guarantee. Higher azimuthal/radial
mode requests explicitly name dependency `M2.1-W3c`; no scalar fallback is used.
Realistic full-band bells remain blocked until a suitable W5 or independently
qualified axisymmetric route exists.

`tests/bake_band.rs` contains strict 1e-3 pulsating-sphere and finite-wall
pipe/Levine-Schwinger qualification tests. The pipe has an outer wall, inward
bore wall, annular lip and backed throat piston, with explicit t/a and panel
ladders. The exact LS load is transported to that actual driven reference plane
by the existing lossless fs-duct owner. These strict tests may expose the current
dense route's limit and must not have their tolerances relaxed. Separate 0.25
declared-tolerance tests exercise software completion, bitwise replay and
partial-support refusal only; they confer no 1e-3 numerical/physical claim.
Native execution remains pending. No converged real instrument bake is claimed.

### Open pressure-jump screens and Robin boundaries (IRP M2.5-W7)

`helmholtz::bc` substitutes the normalized admittance `beta=rho*c*Y` into
**both** existing Burton-Miller rows. For outward normals and `exp(-i omega t)`,
passive wall flow is `v=v_drive-Y*p`, so `q=i*omega*rho*v_drive-i*k*beta*p` and
`(A+i*k*B*beta)*p=B*q_drive`. A local coefficient may vary with position and
frequency; a `BoundaryOperator` applies an explicitly supplied nonlocal law
(including a qualified tangential BL operator). It is not approximated by a
diagonal. Its fingerprint binds checkpoints, not scientific authority. The
caller owns the nonlocal operator's units, discretization and weighted-work
passivity. Re(beta)<0 refuses unless declared active; active results retain
that declaration and cannot establish passive acceptance.

The rigid/all-zero-beta path calls the existing dense or iterative owner with
unchanged pressure/velocity arithmetic. Closed Robin geometry requires retained
outward, conforming closed triangles; topology admission reuses `near_field`.
The dense row pair and LU are the existing owner; `.109` supplies FMM/GMRES,
preconditioning, original-equation residuals, cancellation/checkpoints and work
receipts. `solve_radiation_bc_auto` selects it above 8192 panels. Preconditioning
continues to use the rigid near blocks; it does not replace the Robin operator
or its residual. Plane-wave scattering supplies the incident Robin datum to
the same equation; wall absorption uses the total incident+scattered pressure,
while the outgoing field remains the scattered trace. Reported surface work is
not an independent far-field accuracy certificate. The inherited dense solve
has no mid-factorization cancellation guarantee.

`helmholtz::bc::incident::IncidentField::Sources` extends the finite-source
owner above to coherent mixtures of plane waves, volume monopoles and force dipoles, using the same
analytic kernels and velocity convention. `IncidentTraces` samples that field
or accepts pressure/normal-velocity traces from a caller-owned Helmholtz field.
`ScatteringSolution` retains incident pressure and normal velocity separately
from the outgoing scattered solution, with explicit total boundary accessors.
Absorption uses total wall pressure, including interference; its passivity
tolerance scales with apparent wall power. A field containing a plane wave
has no outgoing total `1/r` far field and the existing observer refuses it.

The analytic source list is bounded at 64. Volume and force positions share finite
receivers' component-wise winding and point-to-triangle clearance admission,
including the original iterative surface range; interior, touching and
unresolved source positions refuse. A positive source clearance is required,
but it is not a near-source discretization error bound. Analytic sampling is
fixed-order and binds traces immutably to geometry, frequency and medium.
External traces must represent a field regular in the bodies; their producer
owns that PDE/transfer obligation. This lane is uniform, stationary, lossless
acoustics; it does not extend molecular-relaxation or convected-wave physics.

`PreparedBoundary` factors the same dense Robin operator once and reuses its
matrix/LU across separate calls and source batches, with the inherited
256-field batch cap. Rigid/all-zero admittance uses the original radiation
operator and batch worker scheduling unchanged. `scatter_field` retains the
source-bound finite receiver observer; `scatter` consumes bound traces.
Local and nonlocal wall laws are retained; a changed nonlocal fingerprint
refuses reuse. `solve_incident_bc_iterative` and
`incident::solve_scattering_iterative` send the same coherent traces through
the existing bounded FMM/GMRES owner. Its cancellation, residual, resolution
and work admissions remain in force. The finite receiver observer retains
its independent dense panel cap. Dense preparation retains its original
cancellation limitations.

`tests/incident_scattering.rs` covers analytic monopole Euler velocity and
conserved outward power (G0/G2), coherent cancellation/superposition and total
Robin absorption (G0/G3), agreement with the original plane-wave solve,
rigid/all-zero-beta bit preservation and prepared-batch replay (G5), and
source/trace admission. An independent spherical-Bessel addition-theorem
oracle supplies a three-rung finite-source Robin-sphere far-field ladder
(G1/G2), with an explicit 3% finest-rung gate, not a 1e-3 qualification.

`reciprocal_radiation::open_screen` admits indexed oriented open two-manifolds,
with every component having a free edge. It refuses closed components,
duplicate faces, inconsistent/nonmanifold edges, disconnected vertex links and
degenerate geometry. This is not a global self-intersection certificate. The
unknown is the P1 pressure jump `[p]`, zero at free-edge vertices, in the
surface-curl weak hypersingular formulation. Geometry and graded-edge level are
supplied as data; the solver never fabricates a closed thin body or silently
changes the mesh. A quadratic edge-graded disc mesh is test data. The observed
sqrt-edge convergence order is reported with `fs_mms`, not assumed smooth.

Analytic Laplace triangle moments, Duffy rules and the outgoing angular Gram
identity share the existing reciprocal owner. The dynamic remainder uses the
existing real/damped `WaveNumber` seam. The existing RHP `solve_laplace` uses the
full complex weak kernel and makes no real-frequency power claim. Far fields
use the jump's double-layer potential; receivers use the existing adaptive P1
near-field integration. Both screen sides are represented without coincident
opposite faces. A screen's optional sheet law is `v=v_drive-Y*[p]`. Two
independently coated faces would require an additional derivative-jump trace
and are outside this single-jump model.

Open P1 assembly has explicit triangle order (4/8), angular refinement (0..2),
edge-level, unknown, retained-cache and scalar-work ceilings. Its dense panel
and unknown cap is 8192. **The closed constant-panel `.109` operator does not
implement an over-cap open P1 hypersingular route**; open requests above this
cap refuse. Logical cache/work ceilings are not RSS or timing guarantees. The
open owner polls cancellation through preparation, assembly, LU rows, RHS and
residuals; arbitrary caller-supplied nonlocal applications remain indivisible.
There is no repair of negative resistance, clipping or tolerance promotion.

`OpenOptions.workers` runs four steps on declared threads: the static and
dynamic facet-pair moments, the outgoing Gram form, and the LU (through
`fs_la::eigen_complex::with_lu_workers`).
- The results are bitwise those of one worker. Pair values are computed in
  row blocks and scattered in the serial order. Each Gram entry takes its
  directions in the serial order.
- One worker keeps the serial loops and their per-pair poll points. More
  workers poll once per block of `4 * workers` facet rows.
- Zero workers refuse. `open_screen_workers_are_bitwise` checks 2, 3 and 7
  workers at a real and a damped frequency.

`PreparedScreen::solve_many` solves several real-frequency drives against one
assembly and LU, and is bitwise that many separate `solve` calls
(`open_screen_many.rs`, workers 1, 2, 3 and 7, rigid and local admittance).
Each drive passes the single solve's backward-residual and passive-work
checks; one failure, a drive count above the declared ceiling, or a
cancellation (also polled between back-substitutions) refuses the whole call.

`PreparedScreen::solve_velocity_ports` adds complete structural self and mutual
loading on this same impermeable screen. Each real source trace retains all
three normal-velocity corner values of every oriented facet. Curved/folded
facets may therefore use their own normal projections of one shared structural
nodal vector; no averaged normal or free-edge velocity clamp is introduced.
The exact P1 triangle mass pairing supplies both the weak acoustic drive and
its mechanical work transpose. With `exp(-i*w*t)`, the unchanged matrix is
`Z_rs = integral b_r [p_s] dS`, the fluid force is `-Z*rate`, and source power
is `0.5 Re(rate^H Z rate)`. The jump already represents both wetted faces in
one fluid, without a factor two or an invented closing wall.

One original assembly/LU serves every source. Declared whole-matrix, retained
jump, source-count and solve/projection-work caps precede batch allocation;
all original geometry, resolution and source-residual gates remain in force.
The result retains the original unit-source `JumpSolution` columns. Coherent
near pressure, far field and radiated power combine those jumps before applying
the original P1 observer once, retaining every interference term. Failed or
cancelled batches publish no partial matrix. This is real-frequency response;
it supplies neither a causal loading realization nor complex loaded poles.

Three focused tests authored in `open_screen_many.rs` compare work with
independent degree-two triangle quadrature, coherent fields with an original
combined-drive solve, folded-facet traces with an independently projected
equivalent weak load, and whole-batch caps/cancellation/replay. These tests have
not been executed in the current workspace with its missing pinned dependencies.

`PreparedScreen::complex_operator` supplies a separate, bounded holomorphic
continuation of the **same discrete real-frequency operator** for loaded-pole
pencils. It retains the analytic static terms and cosine Green remainder, then
continues the outgoing angular Gram form through symmetric antipodal products
`(c_-,i*c_+,j + c_+,i*c_-,j)/2`, evaluated at one complex wavenumber without
conjugation. Real-axis matrix arithmetic is the original assembly order; both
`H` and its analytic derivative `dH/dk` are returned in the original active P1
pressure ordering. Their units are metres and square metres. Zero belongs to
this chart; there is no closed-body `1/k` coupling or pressure elimination.

`velocity_port_pairing` exposes the original exact real P1 mass pairing `G`,
mode-major, and `pressure_dof_map` exposes the original zero-free-edge pressure
map. With `exp(-i*w*t)` and modal displacement `x`, the augmented shell rows are
`(K-w^2 I-i*w*C)*x + G*p = 0` and
`H(w/c)*p - rho*w^2*G^T*x = 0`. Thus the source and work transpose keep the
same sign and coordinate units as the real-frequency forced response. No
inverse of `H` or fitted impedance is needed to define the full pole pencil.

`ScreenComplexLimits` admits both matrices, bounded parallel pair scratch,
angular/coefficient payloads and derivative assembly work before allocation.
The original frequency band and preparation work limit still apply. Whole
contours can be preflighted with `admit_complex_domain`; every subsequent
sample checks its own domain as well. Both signs of the imaginary part are
admitted only when `|Im(k)|*2R` fits the declared positive growth limit, at most
8. `R` covers actual cached samples and vertices about a fixed geometry center,
so translating a screen cannot create exponentially large individual factors.
Original deterministic pair order and block cancellation remain in force;
failure returns no partial operator. These are logical payload/work bounds,
not process RSS or accuracy certificates.

This finite operator pole chart does not certify positive-real behavior over
the whole Laplace half-plane, supply a production CQ reaction, or alter
`solve_laplace`. Loaded poles still require the mesh/quadrature/angular ladder
and the downstream physical residual checks. Two tests of the exact production
scalar kernel passed with the native compiler and production `fs_math`: analytic
derivatives in independent complex directions, zero limits, acoustic reality
symmetry and translation-invariant antipodal products. Three additional matrix
tests cover original real-axis bits and work projection, full analytic matrix
derivatives, domain/resource refusal and cancellation/replay; these are authored
but have not run with the unavailable pinned package dependencies.

### Molecular gas loading and receivers on two-sided open screens

`PreparedScreen::solve_molecular_velocity_ports` binds the original evaluated
`MolecularExterior` to the existing P1 pressure-jump screen. The original
material supplies complex k at a real physical frequency; the weak drive is
`H(k)*jump = i*rho*omega*G^T*rate`. The real physical `rho*omega` is preserved,
and the same exact P1 G supplies every self/mutual resisting effort. Both faces
already belong to this one jump, with no duplicated wet area or closing wall.

`complex_matrix` is a values-only specialization of the original centered
cosine/antipodal continuation. `ScreenComplexValueLimits` counts one H; the
unchanged `ScreenComplexLimits` counts H and its derivative. The new path forms
no derivative matrix or derivative pair moments. One H/LU serves every source
column. Whole-batch work includes assembly, factorization, solves, residuals
and full work projection. Source and solve payloads, original geometry limits,
material band, spatial resolution and centered exponential-growth chart all
remain explicit. These are logical admission bounds, not timing/RSS claims.

The retained `MolecularScreenVelocityPortSolution` exposes the complete matrix,
original jump columns and signed coherent boundary work. That work includes
exterior molecular heat and outgoing energy; it is not distant-sphere flux.
No repair or whole-half-plane passivity/causality certificate is introduced.
The medium is stationary with inviscid momentum and molecular compression;
viscothermal boundary layers, independent face coatings and cavity backing
remain outside this owner.

`prepare_receivers` borrows that exact immutable solution and integrates
reusable complex P1 pressure/gradient rows. It keeps exact retained triangles,
both screen sides, independent componentwise quadrature gates, clearance and
row/work/depth limits. Applying actual complex rates combines the original
jumps before observing pressure and Euler velocity, retaining interference
without another solve or kernel integration. Error estimates cover integration
of fixed solved jumps, not boundary mesh/material error or relative accuracy
near coherent cancellation. Six focused source/receiver regressions cover
native lossless recovery, independent P1 work, physical molecular forcing,
absorbing P1 field/gradient references, coherence and admission/cancellation.
They are authored and syntax checked; native execution remains pending the
existing missing dependency and unavailable build-runner blockers.

`MolecularExterior::prepare_incident` retains one evaluated molecular wave and
the original coherent source declaration. Its reusable analytic sampler and
`with_field` avoid repeated material evaluation at quadrature points or across
incident cases. Physical plane pressure, monopole volume flow and force-dipole
strengths use the original shared kernels and real `rho*omega`; whole point
batches have explicit source-point work admission and transactional refusal.

`PreparedScreen::solve_molecular_illumination` adds reciprocal sound excitation
of the same two-sided screen. One original H/LU serves every mechanical source
column and every blocked incident case. The blocked drive is the weak integral
of `-incident_velocity.dot(normal)` at the original triangle quadrature points.
It is not a centroid or free-edge-clamped velocity. Incoming analytic pressure
is continuous across this zero-thickness surface, so the blocked scattered
jump alone supplies force `-G*jump` on structural coordinates. Finite emitters
may lie on either side but must clear every retained triangle face and edge.
Declared quadrature order and source clearance do not certify incident-load
spatial accuracy; that remains a separate source/mesh/quadrature refinement.

The result retains original blocked jumps, forces and immutable illuminations.
Its prepared receivers cache analytic incident fields and reuse one P1 row
bank for the coherent blocked-plus-motion jump; incident pressure and velocity
are added exactly once. Budgets include whole-batch boundary, sampling,
projection and receiver payloads, with bounded cancellation and no partial
publication. Signed total interface work under illumination is not outgoing
scattered power or gas absorption alone. Five focused sampler/screen regressions
cover independent source pressure/gradient, native blocked-screen recovery,
independent weak drive/work and total-field integration, interference and
resource/source-clearance/cancellation boundaries. They are authored and syntax
checked; native execution remains blocked by the existing unavailable runner
and missing sibling dependency.

`reciprocal_radiation::pipe_end` (IRP .100, RULING (2)) measures the open
end of a thin unflanged pipe on this owner, in Levine and Schwinger's own
geometry.
- The pipe is a zero-thickness "cup": a cylinder closed by a hemisphere, so
  its only free edge is the mouth rim. The rim is quadratically graded.
- Two source rings (declared radius and points; the tests use eight points
  on the node ring) with strengths 1 and `-exp(-i k delta)` cancel their
  backward plane waves, so the forward amplitude cannot vanish. On the node
  ring they do not excite the first radial mode. The screen is rigid.
- Everything behind the stations changes only the incident amplitude, never
  `R = B/A`.
- The total pressure is read on the node ring `r = 0.6276 a` of the first
  radial mode, one station per receiver request. It is fitted to a standing
  wave referred to the rim plane.
- Each row carries `R`, `l/a` from the fit and from the interior zero
  crossing, the fit residual and condition, `k * span` and the station
  pressures.

`tests/pipe_end.rs` measures; it has no accuracy gate, which .100's owner
declares. The dense cap binds on facets (8192), so the rungs are 36, 46 and
58 vertices per ring (3114, 5060 and 8062 facets, ratio about 1.27; the
order estimate is weak at that ratio).
- Against `fs_duct::levine_schwinger` at the finest rung, with node-ring
  sources, for ka 0.05 to 3.8: `|R - R_LS|` is 5.9e-5 to 8.6e-4 and
  `|l/a - LS|` is 3.3e-5 to 5.1e-4. The observed order is about 2.3.
- Axial sources (`source_radius` 0) excite the first radial mode. It is
  barely evanescent at ka 3.8 (`exp(-0.49 z/a)`), reaches the mouth and
  converts, so `R` misses by 3.6-3.7e-2 on every rung. The node-ring source
  takes the miss to 1.0e-3 on the same rung: the prediction was 3x, the
  measurement 36x (`open_pipe_end_node_ring_sources_at_the_radial_cut_on`).
- The cup's finite outer wall leaves a length dependence at ka 1-3 that is
  not yet resolved (.100's gate, ruling 6396). That is why the 1e-3 claim is
  not gated here.
- The interior zero crossing is biased: 1.3e-2 at ka 0.05, and meaningless
  for `k * span >= 0.45`.
- Planted negatives:
  - the conjugate-convention fit flips `l/a` exactly and misses `R` by 1.1-1.9;
  - receivers on the axis at ka 3.6 raise the fit residual 260x (0.25
    against 9.6e-4).
- Cup lengths 4.6 a and 5.6 a move `R` by 1.3e-3 at ka 1 and 9.6e-5 at
  ka 3, within the rung's own discretization error.

`tests/open_and_impedance.rs` authors strict 1e-3 thin-disc mass/complex
impedance/scattering and Robin-sphere pulsating/plane-wave gates, plus G0
passivity/topology/cap tests, G1 edge ladder, G3 fictitious-frequency condition
checks and G5 rigid-bit/repeat checks. The disc oracle is the independent
80..140-digit, M=12/18/24/30 evaluation in `tests/disc_reference.py` of
Mellow & Karkkainen, JASA 118 (2005), 1311–1325,
doi:10.1121/1.2000828, eqs.47–54, b=a; it conjugates the paper's time convention.
No recording or plot was fitted. Native execution is pending: these are gates,
not passing claims. Existing panel-bits goldens are unchanged; there is no
schema migration, new dependency, instrument qualification or render claim.


## Finite-source diffraction by closed bodies

`helmholtz::bc::incident` illuminates the original closed-boundary operator
with a plane wave or a coherent array of at most 64 ideal volume/force sources. A
monopole is specified by its physical position and outward complex volume
flow [m³/s]; the declared `exp(-i omega t)` convention gives
`p=-i*omega*rho*Q*exp(ikr)/(4*pi*r)` and retains the reactive `1/r²` particle
velocity. The scattered drive is `-v_incident-Y*p_incident`, applied to both
original Robin/Burton–Miller rows. Multiple disconnected closed obstacles
share one solve and therefore retain mutual scattering and diffraction.

`ForceDipole` supplies Cartesian peak force ON the fluid [N], not N/m or N/s.
It uses `p=-F.grad(G)` and `u=-Hessian(G)F/(i*omega*rho)`, with the same
radial Green derivatives as the original boundary kernels. Pressure has no
extra density multiplier. The full angular derivative retains tangential
velocity at a pressure-directivity null. `IncidentField::ForceDipoles` and
mixed `IncidentSource::ForceDipole` feed the original trace, dense/batch and
FMM/GMRES scattering owners; no source-pair approximation is used in production.
Iterative source clearance/winding and boundary sampling poll the existing gate.

The source positions pass the existing retained-triangle clearance and
component-wise winding admission. Interior sources, boundary collisions,
empty/over-cap source arrays and unresolved source/receiver coincidences
refuse; point sources are not silently moved. `ScatteredField` borrows the
immutable source array and body. Its finite receiver result adds the analytic
incident field to the original adaptive scattered-field representation. Its
total far field exists for volume/force illumination; a plane wave is not
misrepresented as an outgoing `1/r` amplitude.

`solve_scattering_batch` and `solve_radiation_bc_batch` assemble and factor
once for all admitted right-hand sides, including local and nonlocal Robin
laws. The rigid branch delegates to the existing rigid batch. The original
dense panel/resolution ceilings and 256-field batch cap apply; the dense
assembly/LU cancellation limitation remains. Near-field discrepancy estimates
exclude BEM discretization, ideal-source approximation and material error.
This is homogeneous, lossless, stationary scalar acoustics, not an interior
source, transmission, deformable obstacle, atmospheric-refraction or measured
HRTF model.

`solve_scattering_iterative` sends the identical admitted incident boundary
data through the existing FMM/GMRES Robin owner and returns its original
convergence report and resource receipt. It avoids dense factorization; source
geometry supports the existing iterative panel range while adaptive receiver
geometry retains its 8192-panel ceiling. Cancellation is
checked before source admission and around the existing iterative solve;
source winding and individual FMM passes remain bounded indivisible work.
The inline dense-versus-iterative Robin-source check uses an 80-panel direct
near tree and a 1e-8 pressure agreement band, plus immediate cancellation.
It verifies adapter equivalence, not a new FMM truncation-accuracy claim.

`ScatteredField::total_pressure_velocity_at` observes coherent incident plus
scattered pressure and particle velocity through the original adaptive vector
receiver. `total_first_order_at` applies the existing ideal first-order
microphone responses to that total field, retaining pressure-equivalent units.
The reported discrepancies cover receiver quadrature only; they do not bound
source, material or BEM discretization error. Receivers are ideal and have no
calibration or backreaction. The original 64-receiver and quadrature-work caps
remain, and a focused regression checks the Euler gradient by displaced total
pressure along with omni, cardioid and figure-eight responses.

Inline tests in `helmholtz_incident.rs` cover G0 Euler-gradient consistency
(2e-9 relative) and outward low-frequency volume flux, exterior/singularity
refusals, G5 bitwise batch/individual and historical plane-wave equivalence,
G1/G2 finite-source total-pressure refinement against the independent rigid
sphere series (20/80/320 panels, 2% finest complex-field band), and G3 source
and receiver reciprocity (1e-3 relative). The sphere oracle uses the spherical
addition theorem (NIST DLMF §10.60.1–2) and the rigid condition
`a_l=-j_l'(ka)/h_l^(1)'(ka)`; it does not fit BEM output.

### Finite-source diffraction in the audio observation pipeline

`observation_transfers::ScatteringObservation` binds prescribed point-volume-
flow and directed force sources, fixed closed rigid obstacles and physical receivers to the
existing causal `ObservationTransferBake` / `ObservationFir` pipeline.
`PointSourceSpec` and the original volume-only constructor remain available.
`new_mixed` accepts `CompactSourceSpec::VolumeFlow` and `ForceDipole`; the latter
has an explicit finite unit direction and a scalar force signal [N]. Separate
co-located directions can represent a full Cartesian force with coherent phases.
`SourceId::PointVolumeFlow` [m³/s] and `SourceId::PointForce` [N] cannot be used
as surface-velocity traces. Existing source-kind tags remain unchanged; force
has a distinct tag. The operator identity includes source kinds and directions.

Each frequency shares one boundary assembly/factorization across sources and
one adaptive pressure/velocity Green-row preparation across receivers. Admitted
immutable source geometry is reused across frequencies through the original
prepared-boundary and incident-trace owners. The
analytic direct field and the solved scattered field are added coherently.
The source-major, receiver-major output is `(p, dp/dx, dp/dy, dp/dz)` per
declared source coordinate (m³/s or N);
the gradient uses the original Euler relation `grad(p)=i*omega*rho*v`.
The transfer preserves phase cancellation in acoustic shadows, without a
binary visibility mask or an empirical shadow gain.

Source/receiver exterior winding, actual triangle clearance, mutual separation,
unique identities, finite band and existing quadrature/work caps are admitted.
The straight-line source-to-receiver delay is a conservative earliest arrival,
even across an obstructed segment. The existing four-or-more-rung DFT, tail,
wrap-around and DTFT gates then admit the temporal realization, retaining its
explicit reconstruction delay and unresolved or supplied producer uncertainty.
An exact-analytic producer claim is refused. Cancellation is checked between
bounded stages and propagates as `TransferError::Cancelled`; dense LU and
adaptive receiver preparation retain their existing interruption limits.

A force source has nonzero static pressure and singular zero-frequency Euler
velocity; this positive-frequency wave owner has no admitted static solution.
Unwindowed force DC and FIR bakes explicitly refuse, including the generic
bake entry that otherwise supplies its own zero bin. Force FIR use requires
the existing explicit `LowerBandWindow` through `LocalBoundaryObservation` or
its iterative owner. The band-limited artifact preserves that window; its
zero DC is a declared observation filter, not a physical static-field claim.
Positive-frequency harmonic/coherence observations remain unwindowed and retain
mixed force/volume cross terms, SI pressure, velocity and complex intensity.

`tests/force_dipole.rs` defines independent monopole-pair limit and Euler-gradient
checks, pressure-null angular velocity, actual mixed-source diffraction and
density scaling, dense/iterative parity, source clearance and kind/direction
binding, DC/FIR refusal, and an actual lower-windowed force/flow FIR replay.
These test definitions require native execution before numerical pass claims.

Inline G0/G1/G3/G4 tests cover displaced-receiver pressure-gradient agreement
within 2e-5 relative, distinct source semantics, invalid geometry/band and
cross-wired declarations, cancellation before and inside a frequency solve,
and an actual four-rung BEM-to-FIR bake/replay with complex transfer error
below 2e-4. This qualifies only the declared discrete stationary observer:
it supplies no source loading, moving geometry, finite-aperture equivalence,
spatial convergence certificate, or measured instrument/head/room response.

### Material and compliant patches in the observation pipeline

`observation_transfers::LocalBoundary` assigns a finite table of passive laws
to physical boundary panels. `LocalWall::Mechanical` evaluates the original
`fs_phs::WallPin` specific impedance with its native exp(-i omega t) convention;
`LocalWall::PorousRigid` uses `PorousRigidBoundary` and its existing single
conversion from JCA(L). Both preserve complex phase and SI specific impedance,
without a panel-area factor. Each distinct law is evaluated once per frequency,
then expanded to its declared patches. The table admits at most 256 laws and
the original dense panel ceiling. The complete constitutive state, including
porous source support and the original fluid bulk modulus, binds its hash.
Porous saturating-fluid density and sound speed must match the exterior medium.

`LocalBoundaryObservation::for_scattering` extends the existing prescribed
point-flow observer with these walls; `for_surface` consumes the original
finite structural-velocity or unit-aperture-flow traces. The original prepared
Robin owner resolves all patches and closed components together, reusing one
boundary factorization and adaptive pressure/velocity receiver preparation per
frequency across all sources. Each panel with nonzero prescribed source drive
must be explicitly rigid: placing a passive shunt there would change that
source's actual prescribed velocity or volume flow. The surface adapter admits
the existing centroid owner and refuses routed/triangle-ruled owners rather
than silently changing their discretization. It produces stationary observer
transfers; it does not supply a new in-loop load or moving-frame constitutive law.

`physical_frequency_sample` retains unwindowed material physics and propagates
native material/mechanical errors. A declared band containing a lossless wall's
zero-impedance resonance refuses before sampling; no admittance is clipped.
Source cards with a positive lower frequency require an explicit
`LowerBandWindow`. Its C4 lower stop/pass transition is bound into the operator
hash and retained by `BoundaryObservationBake`, alongside the original upper
band transfer and boundary hash. `band_limited_frequency_sample` is exactly
zero in the declared lower stop region and does not query unsupported material
frequencies there; it still checks cancellation and the declared upper band.
The lower stop must cover every selected material's unsupported lower range,
and the upper band must remain within their common source support. This is an
explicit band-limited reconstruction, never a full-band material extrapolation.

The existing FIR owner retains its reconstruction latency, four-or-more-rung
ladder, 1e-10-or-stricter tail/wrap gates and 1e-4-or-stricter inter-rung DTFT
gate. Windowed DTFT checks span the actual unity passband. Every rung must
resolve at least four unity-passband bins and four bin intervals across each
transition, preventing an empty/narrow sampled band from appearing converged.
These are necessary resolution checks, not bounds for unsampled resonances.
Original full-band bake signatures and arithmetic remain unchanged. Focused
tests cover constitutive phase/units, state and patch identity, rigid parity,
total Robin pressure gradients, passive work, preserved finite-aperture flow,
source-support and cancellation refusals, and porous BEM-to-FIR playback.
Authored tests require native execution; source review and formatting alone
establish no passing numerical or experimental qualification.

### Iterative material diffraction and large-surface observation

`helmholtz::bc::PreparedIterativeBoundary` retains one original
`IterativeRadiation::prepare_bc` FMM/preconditioner for repeated prescribed or
incident source solves. `solve` and `scatter` use the original Robin drive,
actual boundary velocity and total-wall work finalizers. Every source must pass
the original GMRES residual gate. Native nonconvergence, FMM work/accuracy
refusals and cancellation diagnostics remain intact; a changed nonlocal wall
fingerprint refuses reuse. Resource receipts include the shared setup time,
so summing receipts does not measure disjoint preparation costs.

`observation_transfers::IterativeBoundaryObservation` uses that prepared owner
for point-flow/force sources and finite prescribed structural/aperture traces.
`for_mixed_scattering` retains the distinct SI source coordinates. One
frequency shares both the FMM/preconditioner and adaptive pressure/velocity
receiver rows across all source columns. Point-source topology/clearance is
admitted against immutable geometry once; frequency sampling reuses the
original incident trace sampler and adds the exact direct field coherently.
Finite traces retain the disjoint-drive/passive-patch rule and SI normalization.
`IterativeObservationSample` carries complete p/gradient columns, actual
per-source solver reports, resource receipts and total-wall absorption.

The explicit iterative constructors for the material table and receiver
geometry admit up to the existing `MAX_SURFACE_PANELS` (81,920), while the old
dense constructors keep their 8,192-panel limit. Wavelength resolution, closed
outward topology, componentwise exterior winding, triangle clearance and
constitutive support are unchanged. Retained finite source coefficients must
fit `MAX_OBSERVATION_ENTRIES` before copying. Receiver integration keeps its
existing depth, phase/separation and at-most-20-million kernel-evaluation
budgets; large geometry admission does not guarantee that a requested receiver
table fits those budgets. Iterative receiver admission polls cancellation in
topology, clearance and winding tiles; vector preparation polls at recursive
triangle tiles. Original FMM setup and individual passes retain their declared
bounded interruption limitations.

`IterativeObservationOptions` records the explicit FMM/GMRES, receiver and band
settings. Every solver/preconditioner/directional/rule setting binds the new
operator identity. The bake retains those settings and the original material
band metadata; all nonzero source solves must succeed before the existing FIR
owner can admit temporal convergence. Generic private error propagation does
not change the public dense FIR API or its arithmetic. GMRES residual does not
bound FMM approximation, receiver quadrature, spatial discretization or source
model error; those remain separate from FIR tail/wrap/DTFT convergence.

Focused regressions define a three-order FMM-to-dense complete-field ladder,
finite-aperture/passive-work agreement, native GMRES failure propagation into
FIR baking, a 20,480-panel constructor with original geometry/work refusals,
porous band-limited playback and changed-nonlocal-law refusal. Native execution
is required before claiming these numerical regressions pass. This extension
does not qualify measured player/chair geometry or close the full W9 milestone.

### Shared material diffraction, finite-source reaction and angular fields

`LocalBoundaryObservation::combine_surface_and_scattering` and its iterative
counterpart consume already-admitted finite-source and compact-source owners
to create one coherent bank. Finite structural velocity/aperture flow columns
precede compact volume-flow/force columns and keep their original SI units.
The exact surface object, medium, receiver identities/order, clearance,
material boundary, frequency limit, receiver quadrature and lower window must
agree; iterative solver settings must also agree. Duplicate names and combined
source/output extents exceeding the original caps refuse before any solve.

One boundary preparation and receiver quadrature table serve the whole bank
at each frequency. Prescribed surface velocity and incident Robin traces reach
the original solve owner; analytic compact pressure and full particle velocity
are added exactly once. Retained source dispatch carries that same distinction
through harmonic statistics, microphones, directivity and positive angular
power. Cross terms therefore survive between a signed wall mode, an aperture,
a volume injection and a force on the fluid, including destructive interference.
Original source clearance, passive-patch disjointness and cancellation remain
active. A force column refuses raw DC and an unwindowed FIR; a matching explicit
lower window is retained. A mixed compact basis still refuses finite velocity
self-impedance extraction: no point-source self-work port or automatic
structural loading law is invented.

Focused G1/G3 tests define an independently assembled combined Robin drive,
dense/iterative total p/u parity, signed-coordinate reversal, microphone nulls
with nonzero velocity/outgoing power, factorized partial coherence and exact
direct compact far-field terms. Identity, source clearance, cancellation and
window refusals use the original admission paths. These algebraic tests do not
establish mesh convergence or measured source fidelity.

`IterativeBoundaryObservation::solve_frequency` retains the original pressure
and actual normal velocity for every source at one positive physical frequency.
One shared FMM/preconditioner solves every original source column; the streaming
receiver path and retained path use the same solve traversal and finalizers.
The caller explicitly budgets both pressure and velocity coefficients before
any source solve, within `MAX_OBSERVATION_ENTRIES`. Native material support and
solver refusals propagate unchanged; DC and unsupported material frequencies
do not manufacture retained fields. This response is unwindowed even when its
observer declares a lower FIR reconstruction window.

`IterativeBoundaryResponse::observe_receivers` applies the original adaptive
vector observation owner without another boundary solve. Finite radiator
fields are the total outgoing field. Compact-source traces are the scattered
field; the exact original unit-flow monopole or unit-force dipole is added once
at receivers and exposed through typed `source_incident` for angular consumers.
The original `source_monopole` adapter remains volume-only and explicitly
refuses a force column. Immutable geometry, medium, source
coordinates, materials and numerical settings remain bound by the borrowed
observer and its operator identity. Individual original GMRES and resource
reports remain available; sharing field arrays does not upgrade their error
claims.

`velocity_port_impedance` projects those same finite-source pressure fields
against the original real velocity basis: `Z_ij = integral b_i p_j dA`.
Disconnected passive material bodies participate in the same solve, so their
multiple scattering changes both the self/mutual source reaction and observed
sound. Prescribed nonzero source support still requires rigid patches, keeping
each source coordinate's actual velocity unchanged. The matrix is row-major
(effort, rate), with no symmetrization, diagonal replacement or positivity
repair. Structural velocity coordinates retain force as their work-dual effort;
unit-volume-flow apertures retain work-equivalent pressure. Only an all-aperture
basis exposes `acoustic_impedance` in Pa s/m^3. Compact point sources have no finite
surface self-work port and refuse loading extraction.

`resisting_efforts` returns `Z * rates`; the fluid's reaction on the driven
structure has the opposite sign. `input_power_w` evaluates
`Re(rates^H Z rates)/2`, preserving coherent cross terms and phases. For passive
walls this is the sum of outgoing boundary work and wall absorption of the same
combined discrete fields, not the sum of individual source powers. Matrix
entries and all source-pair/panel products have separate explicit budgets;
projection polls cancellation in bounded panel tiles. Transpose reciprocity is
measured on the unrepaired matrix. The normalized smallest Hermitian-work
eigenvalue uses the original aperture/Jacobi diagnostic only within fs-la's
conservative unpolled-work allowance; `passivity_margin == None` explicitly
leaves that diagnostic unresolved for larger port counts while retaining the
complete physical matrix. A negative measured margin remains negative. These
coordinate-dependent diagnostics are not an in-loop reciprocal/passive rational
realization or an error bound for the continuum problem.

`far_field` returns total outgoing amplitudes under `p ~ F exp(ikr)/r` using
the original integrated triangle phase owner and exact direct volume monopoles
and force dipoles. Force amplitudes are `-i k (direction dot F)/(4 pi)` times
the source-position phase, with no density factor. Both direct and boundary
terms use the same canonical unit direction. Original source-coordinate units
remain attached to every column and drive.
`coherent_directivity` first combines complex-driven pressure and velocity
traces, then applies the original refined SH projection; all source phases and
interference remain present without another boundary solve. Explicit angular,
temporary-entry and total kernel-work budgets apply across every refinement
grid. Kernel integration polls between bounded panel tiles. SH captured fraction
and power by degree describe the requested angular representation; they do not
bound BEM/FMM or panel-integration error, and truncated coefficients alone do not
establish total radiated power.

`IterativeBoundaryResponse::prepare_radiation` retains each physical source's
complex outgoing angular column once. `PreparedRadiationObservation` separately
stores caller-requested directions and a declared finite ladder of the SAME
positive Gauss–Legendre/uniform-longitude quadrature used by SH projection.
`PreparedRadiationOptions` requires at least three grids, initial latitude
count at least eight, the original direction ceiling/tolerance, and explicit
aggregate kernel-work and source-entry storage caps. The ladder grows by
`n/2+1` with staggered longitudes; no preparation is returned on cancellation
or an exceeded cap. No new boundary solve occurs.

Repeated coherent, factorized partial-coherence and independent-source
contractions allocate nothing. The `observe_*` methods contract only requested
directions and require no angular convergence; `evaluate_*` also return the
converged full-sphere power/moments. The factor layout is mode-major/source-minor,
with peak second moment `C = B B*`; each entry retains its original coordinate
units (including m³/s, N or a declared surface velocity coordinate). Complex
source fields are summed with compensation BEFORE squaring, then independent
modes are added. Power is the positive angular sum of
`E[|F|²]/(2 rho c)`, in W; requested-direction statistics are W/sr. This is a
spectral-line time/ensemble average, not a density in W/Hz. No rounded source
Gram matrix, diagonal power replacement, sign clipping or wall-loss double
debit is used. First and second directional power moments retain their signed
cross terms; the trace of the second moment is total outgoing power.

Every actual drive/factorization must pass two successive grid-to-grid changes
in total power AND all first/second moments, normalized by the larger total
power. Quiet cancelling fields are checked on their own surviving power,
not against separate source self-powers. Exhausting the cached ladder returns
`AngularPowerNotConverged`; callers must prepare a finer bank before retrying.
Finite quadrature agreement remains an estimate, not an arbitrary-integrand,
continuum, BEM/FMM or panel-quadrature certificate. SH truncation does not enter
this full-field power integral. Runtime product/rank caps and bounded
cancellation reuse the harmonic-observation work contract.

`BemObservation::prepare_radiation_from_solutions` provides the same angular
bank for retained dense surface solutions. `prepare_harmonic_from_solutions`
likewise supplies the original centroid p/gradient receiver rows without
re-solving a source. Both require exact source order/velocity traces,
frequency and medium, with the original surface-fingerprint and finite-field
admission. Their source types and original operator identity are preserved;
these guards do not certify the accuracy of publicly supplied solution data.

Focused regressions define a two-body, two-finite-source material load and
receiver case, independent coherent surface-work balance, obstacle-removal
sensitivity, a three-order iterative-to-dense mutual-loading ladder, signed
structural-coordinate units, point-source total-field parity, budget refusals
and cancellation. These tests require native execution before a passing
numerical claim. The synthetic bodies do not qualify measured player/chair
models, open two-sided shells, a full-band compact aperture model or the full
W9 milestone.

The iterative observer also exposes its exact medium, panel count and declared
maximum angular frequency for downstream physical load admission.
`aperture_support_areas_m2` sums original panel areas on each nonzero finite
unit-flow trace without replacing its velocity profile or renormalizing flow.
`admit_loading_band_hz` checks the complete positive unwindowed interval against
native material support and the original resolution declaration, including an
unsampled lossless zero-impedance resonance. It shares the existing material
resonance check; the observation window still has no authority to extend a
physical loading domain.

## Reciprocal complete source loads and matched receiver transfers

`modal_loading::modal_radiation_loading_reciprocal` exposes the existing
symmetric P1/P0 Galerkin owner's complete work-conjugate velocity-port matrix
through the unchanged strict modal reciprocity and resistance admissions.
Independent localized structural and aperture columns remain intact. Reported
unknowns, retained cache, original source/projection work, angular refinement,
backward residual and flux balance describe the actual bounded discrete owner.
No collocation output is symmetrized and no failed admission is repaired.

`observation_transfers::ReciprocalBemObservation` observes the same retained
P1 pressure and P0 auxiliary density through the original adaptive near-field
receiver, including pressure gradients from that field's particle velocity.
It reuses the existing source/receiver declaration and causal FIR realization.
Cross-wired declarations, out-of-band requests, original geometry/resonance
exclusions, work limits and cancellation refuse. The existing continuum and
finite-band limitations remain; algebraic reciprocity does not establish
measured acoustic accuracy. Focused original-owner equality regressions are
authored; native execution remains pending.

## Source-owned molecular exterior reaction and receivers

`helmholtz::molecular::MolecularExterior` binds the original evaluated
`RelaxationState` and an explicit positive angular-frequency band to the existing
complex Helmholtz owner. The relaxing caloric response gives
`k² = rho*omega²/K(omega)`; the physical Neumann trace remains
`q = i*rho*omega*v` for exp(-i omega t). Molecular dispersion therefore changes
both the boundary reaction and exterior propagation. No equilibrium-speed
substitution into `rho*c*k`, attenuation-to-modulus conversion, duplicate
molecular kernel or downstream absorption filter is used.

One batch shares its boundary assembly and LU. The existing complex owner now
has a values-only specialization with two dense matrices instead of four,
omitting derivative evaluation/storage on identical quadrature. Original
derivative APIs and real-frequency paths retain their arithmetic. Resolution
admission uses |k|, including the absorption length. Explicit panel/matrix/work
budgets and original field caps remain. Cancellation is checked at admission,
each assembly/RHS row and between fields. The original cancellable LU polls
at elimination and substitution boundaries, retaining its declared-worker block
policy. Nothing is published on refusal.

`MolecularRadiation` retains immutable source/frequency/geometry identity and
uses the existing centroid Green representation for finite pressure and
Cartesian gradient. Receiver geometry/clearance checks are unchanged. Its raw
boundary work is neither clipped nor equated with distant radiated flux: it
includes energy ultimately absorbed in the exterior gas.

`solve_boundary_batch` also admits the existing local or nonlocal surface
admittance, including mixed rigid/absorbing patches. Its normalization remains
`beta = rho*c_eq*Y`; the physical surface condition is
`v = v_drive - Y p`. Substitution gives
`(A + i*(omega/c_eq)*B*beta)p = B*(i*rho*omega*v_drive)` on the original
complex-k boundary pair. Complex k is used for propagation and never replaces
omega/c_eq in this boundary conversion. Nonzero laws require the original
closed-triangle admission. Zero laws preserve the original rigid arithmetic.
Local laws add O(N²) work with no extra dense matrix; nonlocal laws preserve
their full image through N applications and O(N³) dense composition. Their
application work remains the law owner's responsibility; ComplexOptions still
bounds the original complex assembly. Cancellation reaches each composition
column/row, the original LU and each drive, without partial publication.

`MolecularBoundarySolution` retains the actual loaded velocity and reports
drive work, surface absorption and exterior boundary work separately.
Surface absorption is checked against its declared passive/active policy;
the exterior term includes gas absorption and outgoing flux. No term is
clipped or debited twice, and a caller-owned nonlocal law must retain its
fingerprint during the solve. G1/G3 tests use a three-mesh analytic Robin
sphere, empty-reservoir agreement with the original mixed-wall owner,
independent rank-one pressure feedback for a spatially nonlocal backing,
batch superposition and the drive/surface/exterior work balance.

`aperture_response` reuses original aperture admission, force projection,
reciprocity/Hermitian-part diagnostics and `ApertureImpedance::acoustic(0)`.
This matrix is directly compatible with the existing `CoupledTermination` and
molecular duct multiport closure, without another production dependency.
`observe_volume_flows` superposes the retained boundary traces before a single
receiver traversal; solved outward opening flows and their reaction see one
source gas. Generalized shapes with zero net flow keep their mechanical matrix
meaning but refuse conversion to plane-port flow observation.

Focused tests cover values/derivative-owner bit parity and half-storage budget,
a three-mesh absorbing-sphere impedance/receiver/gradient ladder, empty-source
parity and source-rate liveness, original molecular duct/full-aperture closure,
and band/geometry/work/cancellation refusal. Native execution results are
recorded with the implementation; the test definitions alone are not a claim
of full workspace or instrument qualification.

This harmonic model has inviscid momentum and molecular relaxing compression
in one spatially uniform, stationary exterior gas. It does not model exterior
shear viscosity, conductive or viscous boundary layers, mean flow, evolving gas
inventories, nonuniform exteriors, a rational time-domain load, or experimental
source accuracy. The original centroid/disc Burton–Miller and triangle-CBIE
accuracy limits remain. No passivity repair or BEM discretization-error bound
is introduced, and receiver near-field convergence still requires refinement.

### Adaptive molecular receivers on the original boundary traces

`MolecularRadiation::prepare_near_field` and `observe_near_field` extend the
original retained-triangle receiver to passive complex wavenumbers. The real
and molecular routes share quadrature rules, subdivision and componentwise
pressure/gradient discrepancy tests. Complex |k| resolves phase and attenuation
on each tile; physical velocity still uses the original real rho*omega.
Real-k kernel arithmetic and the existing P1 receiver consumer are retained.
Geometry winding/clearance, maximum depth, complete-array kernel work and
explicit retained-row bytes are admitted. Cancellation reaches integration
tiles and contractions; exhaustion returns no partial field or centroid fallback.

`PreparedMolecularNearField` binds original source ownership, surface,
frequency and receiver positions. One integrated row bank serves all source
columns. It returns paired pressure/particle velocity and propagated absolute
quadrature discrepancy estimates. These estimates cover integration of the
fixed solved traces only, excluding boundary discretization, material/source
error and cancellation-conditioned relative error. Receiver positions/output,
geometry scratch, allocator metadata and stacks remain outside the row-byte cap.

Coherent molecular scattering adds the original analytic incident field once
at the prepared positions. Its prepared-row query reuses quadrature across an
incident batch. Aperture-flow observation coherently combines the original
loaded traces before one adaptive contraction. `MolecularBemObservation::
with_near_field` selects the same receiver machinery for harmonic/FIR samples
and matched source responses; controls enter the operator hash. The source
bank retains unit-source discrepancy estimates and actual integration count,
and propagates estimates with absolute source-rate weights. Boundary loads,
wall/exterior work and material laws are unchanged. Adaptive source responses
reserve another `2*s*r` complex-entry equivalents for four real discrepancy
components per source/receiver, in addition to their separate row-storage cap.

G1 source/core regressions cover independent complex Green derivatives and
refined triangle integration near faces/edges, native real-k recovery, coherent
row reuse, original momentum, source/frequency/geometry binding, budgets and
late cancellation/replay. A source-bank case checks unchanged loading and
shared-row discrepancy propagation. They are authored and syntax checked;
native execution remains pending the dependency/resource blocker. Scalar
oracle checks qualify only their reference formulas, not production BEM passes.

### Material walls and matched molecular source response

`MolecularBemObservation::with_local_boundary` attaches the existing immutable
`LocalBoundary` mechanical and sourced porous laws to this same exterior.
Each law is evaluated once per frequency with its original phase, specific
impedance, physical band and equilibrium fluid. Every panel supporting a
prescribed structural velocity or aperture volume flow must remain explicitly
rigid; an absorbing shunt cannot silently change that source coordinate.
Original wall identity enters the observer hash. Reaction, pressure gradients,
harmonic intensity and microphone statistics consume the same Robin solutions.
Mechanical-wall FIR uses the existing observer. Porous laws with a positive
lower support remain in-band harmonic/reaction inputs: DC and this owner's FIR
bake refuse because its artifact has no lower-band window. Their original
equivalent-fluid pore model is retained; molecular pore relaxation is not added.

`prepare_response` returns `PreparedMolecularResponse` with the complete
work-conjugate source impedance and matched `PreparedHarmonicObservation`.
All reaction rows and receiver fields reuse one original solution bank, with
one assembly/LU per bounded source batch. The input, exterior and wall matrices
project pressure against the prescribed, actual, and absorbed velocities,
respectively. Separate direct projections preserve small wall losses and
coherent mutual terms. Arbitrary complex source rates yield resisting efforts
and signed mean powers through allocation-free O(sources²) queries, with no
additional BEM solve, clipping, symmetrization or passivity repair. An all-flow
source bank exposes its pressure/volume-flow impedance in the original source
order for a duct multiport closure; mixed mechanical coordinates refuse that
conversion. The reaction on a driven structure is opposite to resisting effort.

The same already-admitted source bank now retains its original pressure and
actual normal-velocity traces. `evaluate_boundary_coherent` reconstructs both
in original panel order with compensated sums and no allocation or additional
solve. Its `2*s*p` complex products remain within the source-bank reservation;
the retained velocity includes the original passive wall response. Consumers
can therefore contract incident/motion interference directly at the boundary.
Output shapes and rates are admitted before writes; a late cancellation or
overflow requires discarding the output buffers and leaves the bank reusable.

Additional explicit caps admit `3*s² + 2*s*p + 8*s*r` live complex entries
outside assembly and `3*s²*p` projection products, for s sources, p panels and
r receivers. Cancellation reaches source solves, projection rows/panel tiles
and receiver queries. Source/medium/band, budget and finite-value refusals are
preserved. These sampled loads are not a fitted passive time-domain model.
New G0/G1/G3 tests cover independent mechanical-wall pressure feedback and
resistive work, original porous source support, mixed source/receiver coherence,
three-mesh loaded-sphere reaction and receiver convergence, and admission.
Those new native regressions are authored but unexecuted in the authoring
environment: pinned sibling dependencies and shared storage/memory pressure
blocked compilation. Formatting does not qualify their numerical assertions.

### Coherent finite-source scattering in relaxing gas

`MolecularExterior::solve_scattering_batch` now illuminates the same original
closed-body Robin operator with coherent plane waves, nearby volume sources,
compact forces, or mixtures of them. The existing incident sampler shares
the original `WaveNumber` Green/derivative seam across real and complex k;
it retains compensated coherent sums and the real-frequency operation graph.
Source strengths remain physical outward flow [m³/s] and total force on fluid
[N]. For exp(-i omega t), `p=-i*rho*omega*Q*G-F.grad(G)` and
`u=grad(p)/(i*rho*omega)`, including angular force-source velocity at pressure
nodes. A dispersive plane wave therefore has `u=k*p*d/(rho*omega)`, not the
equilibrium plane-wave ratio. Its amplitude refers to pressure at the origin.

For a stationary surface, `v_total=-Y*p_total` gives the scattered drive
`-v_incident-Y*p_incident`. Local and complete nonlocal admittances use this
same equation and one complex assembly/LU for the complete incident batch.
Only outgoing scattered traces enter its Green representation; analytic
incident pressure and velocity are added coherently at each total receiver.
Wall absorption uses total pressure, including all interference terms, and
is checked against the original active/passive declaration. Signed total
outward boundary work is separate from the scattered exterior work; neither
is mislabeled as a far-field cross section in an absorbing medium.

The original matrix panel/storage/work admission is shared with assembly and
runs before source preparation. Additional source-point work is bounded for
the whole incident batch and every receiver query. Original closed-component
winding, source/receiver clearance, cardinality, finite-value and source-law
identity checks remain. Cancellation reaches source evaluations, boundary
preparation, the original assembly/LU, total-work reductions and receivers.
This uses the existing centroid finite receiver quadrature; it introduces no
near-field convergence certificate, moving-medium model or duplicate loss.

Six new regressions cover complex source gradients and physical momentum,
lossless-owner parity, coherent nonlocal wall work, early budget and geometry
refusal, preservation of a weak source between large opposing fields, and a
three-mesh Robin sphere with independent complex spherical-Bessel/Hankel
pressure, radial velocity and analytically integrated absorption. The sphere
series follows NIST DLMF 10.53.1, 10.51.2 and 10.60.7. Native execution remains
pending for the compilation blockers recorded above; the authored convergence
thresholds are not measured outcomes.

## Source-bound gas and liquid acoustic media

`helmholtz::Medium::from_material` consumes the shared checked
`fs_material::fluid_acoustics::FluidAcoustics` view. Gas and liquid density and
sound speed come from their original material state without another EOS or an
air default. The original owner retains source identities and query receipts;
`Medium` is the numerical radiation input. The existing raw `Medium` and
infallible gas copy remain low-level numerical entry points.

Focused G1/G3 regressions check coherent-gas bit preservation and stale-state
refusal, manufactured-liquid pulsating-sphere impedance against its analytic
value, and radiated-power scaling with `rho c` at fixed wavenumber. This does not
add liquid attenuation, cavitation or an experimental accuracy claim; existing
boundary-integral discretization and wavelength limits remain.

### Porous material reflection on closed Robin surfaces

`helmholtz::bc::PorousRigidBoundary` resolves the existing
`fs_material::porous::PorousFluidLaw` once per frequency and fluid state.
It retains the original skeleton source and sampled thermoviscous fluid,
and derives both the exterior medium and `beta=rho*c/conj(Zs)` from that
same response. The conjugation preserves the exp(-i omega t) pressure phase;
specific impedance is not divided by panel area. Dense radiation/scattering
and iterative scattering use the existing Robin owners and their work,
geometry, resolution and cancellation limits. The scalar beta is also exposed
for mixed rigid/coated patches. No constitutive kernel or absorber fit is copied.

The complex planar diagnostic uses `R=(cos(theta)-beta)/(cos(theta)+beta)`.
This is explicitly a locally reacting, planar rigid-backing approximation:
no lateral layer propagation, frame motion, curved-layer correction, diffuse
room model or experimental qualification is inferred. Exterior propagation
retains the existing lossless-medium model. `tests/porous_boundary.rs` checks
independent high-precision phase/absorption references, oblique normal-power
balance, shared-fluid substitution, original Robin solve equality, total-wall
work and invalid source support. Authored Rust tests require native execution;
the test definitions alone establish no passing numerical qualification.

`helmholtz::bc::PorousLayersBoundary` borrows a prepared
`fs_material::porous::layers::RigidBackedLayers` and carries its original
porous source laws, layer order and coherent exterior fluid through those same
dense radiation/scattering and iterative scattering routes. The adapter takes
the stack's **normal-incidence** normalized admittance and conjugates it for
the BEM phasor convention; all panel rows then use that fixed locally reacting
coefficient. The material owner's angle-dependent, extended-reaction stack
reflection remains a separate diagnostic, not a nonlocal BEM boundary model.
Empty stacks reduce to rigid walls. A finite pressure reflection at a pressure
node does not imply a representable finite Robin admittance: this conversion
explicitly refuses such a pole. The porous-boundary tests also cover the
single-layer reduction, porous/air-gap phase and wall work, dense/iterative
agreement and original-fluid propagation.


## Penetrable homogeneous-fluid interface

`helmholtz::transmission::FluidInterface` solves reflection, transmission and
diffraction through one bounded homogeneous fluid body, with continuous pressure
and outward normal velocity. Both the exterior and interior retain their original
`FluidAcoustics` owner and one coherent numerical snapshot at a shared angular
frequency. Each source is queried once during preparation; repeated incident
fields reuse the coupled factorization. The snapshots do not independently
confer measured or calibrated material authority.

The scope is stationary, inviscid, nondispersive, positive-density fluids, a
source-free bounded interior and one consistently outward closed triangle
component. Original closure and winding checks do not certify global mesh
nonintersection. Interior sources, disconnected/nested material regions, solid
shear, mean flow, porous-frame motion and bulk attenuation are outside this
interface. Both material wavelengths must satisfy the existing six-panel
resolution gate.

The direct Green identities use `exp(-i omega t)` and the original outgoing
`exp(i k r)` kernel. Pressure is common to both sides; `q_i/rho_i=q_o/rho_o`
enforces velocity continuity. Subtracting interior/exterior trace identities
yields the weakly singular direct transmission combination (Bruno and
Perez-Arancibia, *Proceedings of the Royal Society A* 473 (2017), 20170161,
equations 3.3–3.7, doi:10.1098/rspa.2017.0161). Unknowns `(p,q_o/k_o)` have
matching pressure units. Original off-diagonal constant-panel kernels,
equivalent-disc self integrals and static double-layer closure remain the
discretization.

The implementation solves scattered traces directly using `M-I`. Saved diagonal
contrasts and analytically cancelled outgoing-kernel differences retain weak
density/speed reflections. A midpoint sine identity removes the common
hypersingular static term before floating-point subtraction. Equivalent-disc
self contrasts use a factored series within the original wavelength admission
(`k*a <= sqrt(pi)/3`). Identical exterior/interior acoustic coefficients return
exact zero scattered traces and unchanged incident boundary data.

Preparation admits the doubled dense dimension and caller-declared matrix/vector
byte budget before quadratic allocation. Its peak includes the actual
worker-dependent LU descriptor workspace through the original algebra owner's
`lu_complex_blocked_workspace_bytes`, and is compared with the solve/output
peak. This is a logical payload budget excluding allocator/thread runtime
metadata and stacks, not a process-RSS limit. Assembly polls cancellation;
the original cancellable blocked complex LU and triangular solves remain the
factorization/solve owners. The inherited condition diagnostic is a probe-based
lower bound, not a conditioning certificate.

`IncidentTraces` and `scatter_field` retain the original regular exterior
plane-wave/coherent-monopole sampling and source-clearance admission. The
outgoing `RadiationSolution` feeds existing finite and far-field observers.
`FluidScatteredField::total_pressure_at` adds the same retained incident field.
Interior pressure uses `S_i q_i-D_i p`: negating both total traces reuses the
original adaptive triangle receiver, with exact winding, face/edge clearance
and its existing quadrature/depth/work controls. Its quadrature discrepancy
excludes boundary discretization and material-model error.

Signed total interface work is published separately. Its real part vanishes
for the continuum lossless, source-free interior; a discrete defect is not
reported as absorption or repaired by clipping. Imaginary work is the physical
reactive exchange. Regression controls compare complex fluid-sphere partial-wave
solutions on 20, 80 and 320 panels for density-only, speed-only and combined
contrasts, check actual interior/exterior receivers and finite-source phase,
and resolve the first-Born forward limit down to a `1e-12` speed contrast.

### Molecular gas transmission

`helmholtz::transmission::molecular::MolecularFluidInterface` extends the same
single-body direct transmission combination to two original `RelaxationState`
owners. Each material is evaluated once at the actual angular frequency;
its complex `relaxing_wavenumber` supplies propagation without an additional
Stokes–Kirchhoff term. Both original gas states remain borrowed and their mean
pressures must agree for the stationary interface. Density remains real.
Pressure and physical normal velocity are continuous; the real equilibrium
`k_ref=omega/c_eq,out` only scales `q_out/k_ref=i*rho_out*c_eq,out*v`.
Complex propagation never substitutes `k*c_eq` for physical `rho*omega`.

The real-k owner retains its original kernel callbacks and operation graph.
Both routes share dense assembly, payload admission and factorization.
Complex weak-contrast kernels use factored midpoint moment series; larger
contrasts use the two decaying endpoint exponentials and their polynomials,
avoiding hyperbolic overflow and cancellation between large midpoint squares.
Equivalent-disc self differences retain factored powers. Identical density and
complex wavenumber return exactly zero outgoing traces, including in an
absorbing medium. Repeated coherent plane, volume and force illuminations
reuse the same factorization and preserve source strengths and phase.

`MolecularTransmissionField` retains outgoing and total physical traces,
paired finite exterior/interior pressure and particle velocity, and signed
outward interface work. Its negative real part equals continuum molecular
loss in the source-free interior; numerical defects are never clipped into
absorption. Receivers use the original complex centroid Green derivative and
actual local `rho*omega`, with winding and triangle-clearance admission.
This does not inherit the real owner's adaptive receiver quadrature or confer
spatial accuracy. Both complex wavelengths, doubled dense storage (including
worker scratch), source/receiver product caps and cancellation are admitted.
The scope excludes multiple bodies, mean-flow/material transport, shear loss,
solid/frame dynamics and experimental qualification.

Four authored regressions cover weak/static/extreme-loss kernel limits,
empty-reservoir recovery of the original nondispersive solver, exact molecular
transparency, mixed sources and resource refusals, and three-mesh complex-sphere
pressure/velocity convergence. The independent sphere also compares surface
work to the radial integral of `omega*Im(1/K)*|p|²/2`. Native execution remains
blocked; formatting and an independent scalar oracle check do not establish
production numerical passes.

## Mutually reflecting penetrable fluid bodies

`helmholtz::transmission::multiple::FluidBodies` extends the original fluid
interface to multiple separated bodies, each retaining its original surface
and `FluidAcoustics` source. The common exterior is queried once and each
interior once. A complete coupled factorization and cached exterior Green
transfer are reused across coherent plane, volume, force and mixed incidents.
There is no independent-body superposition or reflection-order truncation.

For the original single-body contrast `C_i` and pressure-scaled outgoing traces
`s_i=(p_i,q_i/k_o)`, the exciting field is `inc_i + sum_j E_ij s_j`, with
`E_ij=[D_o,-k_o S_o; N_o/k_o,-D'_o]`. The coupled rows are therefore
`(I+C_i)s_i + sum_j C_i E_ij s_j = -C_i inc_i`. This reuses the original stable
contrast kernels and preserves weak scattering. Identical-fluid bodies have
the exact algebraic equation `s_i=0`: their zero outgoing columns are removed
from the coupling without changing other bodies' fields. They still receive
and transmit the complete exciting field. All-transparent arrays preserve the
original incident traces and return bitwise zero outgoing coefficients.

Local outgoing coefficients are distinct from physical scattered boundary
traces, which include every other body's contribution. `Solution` exposes
both explicitly. Its complete scattered work includes mutual interference;
per-body total interface work retains signed discretization defects. The
internal outgoing representation is private. Source-bound finite receivers
and bounded, cancellable retained-triangle far-field evaluation reuse the
original observers. Body-indexed interior receivers use each body's own
wavenumber, density, total traces and original winding/clearance checks.

Each body must be one outward closed, non-self-intersecting component. Positive
separation of retained-vertex axis-aligned boxes is conservatively required;
overlapping boxes refuse even when the actual surfaces might be disjoint.
Nested domains, touching bodies and self-intersection certification remain
outside this owner. Both wavelengths on every body pass the original resolution
gate. Peak logical numerical-payload admission includes dense coupled and
transfer matrices, local contrast blocks, original worker-dependent LU scratch,
vectors and copied surface arrays; allocator/runtime metadata and stacks are
excluded. This is not a process-RSS bound.

`tests/fluid_multiple.rs` checks original single-body agreement, exact and mixed
transparency, coherent-source identity, permutation and translation covariance,
mutual interference work, weak-contrast first-Born volume scaling, and the
geometry/source/receiver/resource/cancellation boundaries. An independent
coaxial multiple-sphere reference uses spherical Bessel/Hankel waves,
Gauss-Legendre translation projection and its own small linear solve. Separate
wave-order and 20/80/320-panel-per-body mesh ladders are required. The last mesh
has a declared 5% far-field fixture tolerance; this is not a general accuracy
certificate. Stationary inviscid nondispersive fluids, constant-panel centroid
interaction, equivalent-disc self terms, and source-free interiors retain the
original owner's limitations. Close-interaction accuracy requires refinement.

### Multiple bodies with original molecular absorption and dispersion

`helmholtz::transmission::molecular::multiple::{Body, MolecularFluidBodies,
MolecularMultipleField}` extends these coupled equations to distinct relaxing
gas bodies in one original relaxing exterior. Each `RelaxationState` is
retained and evaluated at the shared frequency; stationary interfaces require
identical mean pressure. Both complex |k| wavelengths constrain every mesh.
Pressure and physical normal velocity remain continuous. The real reference
`k_ref=omega/c_eq,out` scales the pressure coordinates; physical momentum
remains `q=i*rho*omega*v` with the original real density on each side.

The original real-k `FluidBodies` and this molecular owner share geometry and
payload admission, transfer, coupled assembly/LU and physical-trace recovery
through a private engine. The native real-k operation order and cached real
impedance remain unchanged. Molecular blocks use the existing stable complex
contrast differences and `E_ij=[D,-k_ref*S;N/k_ref,-D']`. One retained coupled
factorization includes every mutual reflection for repeated coherent plane,
volume-flow and force incidents. Exact transparent bodies emit zero outgoing
coefficients while receiving all neighbors. Their interior receivers evaluate
the original incident plus global outgoing field, so an artificial interface
does not change the field.

Complete scattered and total p/v remain distinct from local outgoing
coefficients. Per-body interface work and complete scattered work retain
coherent cross terms, reactive work and signed numerical defects. Roundoff
intervals bound panel-dot arithmetic only; absorbing-exterior surface work is
not power at infinity. Interior/exterior p/u receivers use original centroid
Green derivatives, actual rho*omega, winding and clearance admission. No
adaptive near-field or close-body continuum accuracy is inferred. Disjoint
vertex boxes, total panels, peak numerical payload, source/receiver products
and cancellation remain bounded. Receiver geometry scratch, allocator/runtime
metadata and thread stacks are outside the logical numerical-payload cap.

Five new G0/G1/G3 regressions cover empty-reservoir native and single-body
molecular parity, transparent-neighbor removal, coherent cross-work and body
ordering, resource/cancellation refusal, and an independent complex two-sphere
mode-order and 20/80/320-panel-per-body ladder for paired fields and interface
work. The continuum oracle also compares inward work with the original
caloric bulk heat. These definitions passed syntax/formatting checks; native
compilation and BEM test execution remain blocked in the authoring environment.

### Adaptive molecular transmission receivers

Single and multiple molecular transmission fields expose explicit adaptive
scattered, total and interior pressure/velocity queries. They reuse the shared
complex receiver rows and preserve every propagated discrepancy component.
Interior integration reverses the Green sign once with the original interior
rho*omega. Exterior multiple-body integration uses local outgoing coefficients;
using complete scattered traces there would count neighboring fields repeatedly.
A transparent body's interior admits its selected component and the combined
surface, then integrates one global outgoing row bank plus the incident field.
All-panel work/row caps therefore remain effective. Fully transparent fields
retain exact incident values and zero quadrature discrepancy.

Five G0/G1/G3 cases cover the native limit, independent refined integration of
fixed traces on both sides, single/multiple parity, transparent neighbors,
resource/source/geometry admission and cancellation/replay. These adaptive
queries improve integration of existing traces; they do not refine the boundary
solution or certify close-body interaction accuracy. Native tests are authored
but unexecuted in the current environment.

## Cancellation in coherent field evaluation

Analytic incident pressure and particle velocity, constant-panel finite and far
BEM receiver fields, spherical-harmonic reconstruction, and signed surface work use the same
componentwise Neumaier accumulator as prepared harmonic-source statistics.
The low part is retained until the complete contraction finishes, preserving a
weak field between much larger opposing contributions with constant auxiliary
storage. Sources remain complex amplitudes until the physical quadratic
observation; no source-power addition replaces interference.

G0/G1/G3 regressions exercise all analytic source kinds and mixed illuminations,
source ordering, a weak survivor through the prepared scattering owner, an
analytic Green-layer receiver, a spherical-harmonic dipole node, and signed
active/reactive surface work. Deterministic traversal is retained. Compensation
does not certify arbitrary permutation equality, resolve kernel/phase errors,
or bound a nearly cancelling field's relative error; the original spatial,
quadrature, conditioning and material limitations still apply.

## Finite-window observations of prepared BEM fields

`PreparedHarmonicExposure` consumes the existing
`PreparedHarmonicObservation` source-to-receiver pressure/particle-velocity
transfers at each authored positive angular frequency. The canonical owner
and explicit source clock define physical scope, including an empty/dark
observation. Every line must match its producer hash, ordered physical source
coordinates, exact receiver names/positions and original density/sound speed.
Matching receiver names alone cannot authorize mixing different physical
operators. The common emission clock is explicitly declared, not reconstructed
from unrelated phase data.

An optional deterministic mean and explicitly authored real/proper-complex
factor banks retain all source, incident/scattered and cross-frequency
interference. Every factor is zero mean. Real coordinates have E[xi^2]=1;
proper coordinates have E[|xi|^2]=1 and E[xi^2]=0. Distinct factors have zero
conjugate AND unconjugated cross moments, including mixed real/proper pairs.
Different factors and the mean add second moments; frequencies within a
factor do not become independent. Complex source coefficients on real factors
represent general nonproper covariance through a real square root of the
stacked real/imaginary source covariance. Factor amplitudes carry peak source
units and square-root weights. No Gaussian distribution, stationarity or
PSD-bin normalization is inferred. Existing constructors retain all-proper
semantics; explicit mixed statistics use new_with_statistics[_cancellable].

Preparation contracts retained p/u rows in double-double arithmetic, without
another BEM solve, then coalesces exactly equal frequencies in authored order.
Each query reuses the shared `fs_fft::harmonic_window` analytic rad/s kernels
across all receivers and components. Deterministic means integrate against a
unit DC field through the same DD kernel, preserving low components before
cross-frequency cancellation. Output includes pressure/velocity means,
uncentered pressure mean square, the full real symmetric Cartesian velocity
second moment and signed finite-window mean p u [W/m^2]. No plane-wave
pressure-to-intensity substitution, reactive-intensity surrogate, positive
clipping or heating attribution occurs.

The immutable prepared owner retains its original producer uncertainty; these
moments do not certify source fidelity, boundary discretization or quadrature
accuracy. Double-double arithmetic improves cancellation but supplies no
certified relative guarantee at arbitrary nulls. Nonfinite or negative
computed self moments refuse.

At most 4096 supplied lines and 256 factors are admitted. Explicit source-row,
receiver-frequency-pair and logical workspace budgets precede allocation.
Preparation takes O(receivers * sources * lines) transfer work; queries take
O(receivers * sum(distinct frequencies per bank squared)) fixed-size
four-component contractions. Preparation/query cancellation publishes no
partial result and leaves accepted owners reusable. The returned point array
borrows original source/receiver identity and the declared clock.

Eight authored G0/G1/G3/G4/G5 regressions cover two-frequency physical
incident-plus-scattered fields against independent simultaneous-source BEM
solves and real-time integration; exact finite proper-factor phase ensembles;
full velocity covariance and orthogonal factor mixing; negative subperiod
standing-wave energy flux; sub-ulp proper beats; weak coherent source/group
residuals; clock covariance; and original-owner, budget, cancellation and dark
field behavior. Source formatting and focused production-module Rust type analysis passed. Native numerical execution is
pending because the authoring environment's storage is exhausted; source
inspection and independent reference calculations are not a Rust test pass.

An independent JavaScript reconstruction of the two-frequency/four-component
proper-factor fixture agreed with direct Simpson integration over its exact
16-state phase ensemble to 3.96e-16 maximum absolute error. The same analytic
references give -0.001546695268142812 W/m^2 standing-wave flux and
1.378634354255046e-25 Pa^2 for the weak proper beat; these are reference
calculations, not executed Rust results.

## Joint uncertainty of averaged acoustic sensor readouts

`PreparedHarmonicExposure::observe_readout[_cancellable]` observes a linear
sensor that reports its normalized window-averaged pressure and Cartesian
particle velocity. It returns ensemble means for the original receiver array
and centered covariance blocks for explicitly requested ordered receiver pairs.
Pressure/pressure, both pressure/velocity blocks and the complete velocity
cross tensor retain their physical units. Reversing a receiver pair transposes
the block; correlations between different receivers can have either sign.

Each source factor is integrated coherently across ALL its frequencies
before its authored real/proper contraction is formed. Thus this result is
Cov(average(field at receiver A), average(field at receiver B)), rather than
average of the instantaneous field covariance. A full-period linear pressure
readout can have negligible variance while the acoustic mean-square pressure
remains nonzero. Pressure/velocity readout covariance has intensity units but
is not an energy flux. Deterministic fields contribute to the means only;
the centered covariance is formed directly from fluctuations without
subtracting a large squared deterministic mean.

The original source factors, transfers, geometry, medium and common clock are
retained. There is no additional BEM solve or propagation approximation.
Each complex integrated factor retains DD components using the shared unit-DC
real-product contraction for both its real and imaginary parts. Receiver
covariance then uses Re(A)*Re(B) for a real coordinate and
Re(A*conj(B))/2 for a proper coordinate, preserving the former's
pseudocovariance. Both contractions use the shared DD zero-frequency kernel.
The selected normalized window profile is honored: a triangular detector has
a sinc-squared linear-response envelope, and covariance multiplies integrated
amplitudes with the original factor law. Its physical field is weighted once
before a readout is squared.

These pair covariances support caller-authored real microphone combinations:
mean = w^T m and variance = w^T C w when all needed pairs are requested.
Neither detector electronics/noise, calibration, source fidelity, Gaussian
statistics nor a sampling-bandwidth model is inferred. Original producer
uncertainty remains unresolved.

Per-query work is O(receivers * retained lines + requested pairs * factors),
without a quadratic frequency-pair loop. Explicit product and logical
workspace budgets precede allocation; requested receiver indices and the
65,536-pair hard cap are admitted first. The existing source/frequency/receiver
preparation caps still apply. Cancellation or refusal returns no partial
observation and leaves the immutable prepared owner reusable.

Six authored regressions cover finite-ensemble real-time readout integration
and all covariance blocks; full-period averaging versus nonzero acoustic
energy; common-noise differential cancellation beside a large deterministic
mean; weak grouped factors and clock translation; pairs, budgets, cancellation
and dark scope; and triangular weighted integration. Native formatting and
focused Rust type analysis passed for the exact production modules and six
new tests with the shared 22-test temporal module. This is not a full BEM
crate build or executed Rust tests. An independent triangular four-phase
reference agreed with weighted real-time integration to 5.30e-14 maximum
absolute discrepancy.

## Delayed pressure-array channels from original BEM fields

`PreparedHarmonicPressureArray` prepares explicitly authored causal
delay-and-sum channels from `PreparedHarmonicExposure`. A channel is
`y_c(t)=sum_j weight_j*p_receiver_j(t-delay_j)`. Each tap names an ORIGINAL
physical receiver, a finite signed real dimensionless gain and a finite
nonnegative delay. Multiple taps on one receiver are allowed, tap order is
retained, and no gain normalization or steering direction is inferred.
Channel identities are unique and bounded. Empty taps explicitly mean a dark
channel; zero channels yield an empty observation.

Under the original exp(-i omega t) convention, delay uses
exp(+i omega delay). The original propagated pressure fields are summed with
these complex gains BEFORE any product. Coherence across receivers, sources,
frequencies and the deterministic mean/real/proper-factor decomposition survives.
Delay phases are prepared once per distinct frequency across all banks, with
DD products, near-unit cosine decrements and DD tap accumulation. Repeated
window queries require no BEM solve, spatial reconstruction or delay-phase
evaluation. Transcendental leaves remain deterministic f64; DD is not an
arbitrary-null relative error certificate.

Outputs retain original receiver geometry, physical source coordinates,
producer hash/uncertainty and source clock separately from the processed
channel definitions. Channel pressure is a signal-processing readout, not a
new spatial receiver or physical acoustic intensity. These authored harmonic
continuations do not infer a pulse onset, sampled fractional-delay filter,
microphone calibration, electronic noise or source bandwidth.

Every query returns per-channel E[average(y)] and E[average(y²)], plus the
full real symmetric centered Cov(average(y_a),average(y_b)) matrix. Mean-square
pressure and variance of an averaged pressure are distinct quantities.
Covariance integrates each complete factor before applying its original
real/proper contraction. Deterministic means never contaminate a small
covariance through subtraction of a large squared mean. The copied factor
laws survive after the source exposure is dropped. Real factors preserve
carrier-sum interference in channel mean squares and pseudocovariance in
averaged-channel covariance; proper factors retain their physical one-half
normalization. No independent-microphone, stationary-source or Gaussian
assumption is added.
Uniform and Triangular normalized detector sensitivities are supported.

With C channels, T total taps, K unique frequencies, B banks and L retained
lines, preparation uses O(K*T + K*B*log(L+1) + L*T) work. Pressure-moment
queries cost O(C*sum_b(F_b²) + C*L + B*C²). Logical storage is
O(C*L + T + C²), plus bounded metadata and frequency lists. Fixed ceilings
are 64 channels, 4096 taps and 128 bytes per channel identity. Explicit
preparation, query and logical-workspace budgets are checked before allocation;
the original exposure's source/frequency/factor bounds remain in force.
Cancellation/refusal publishes no partial channel bank or observation.

Six authored G0/G1/G3/G4 regressions cover delayed real-time integration over
an exact 16-state proper ensemble under both shutters; broadband steering and
directional cancellation; identity with original point/readout observers;
sub-ulp cosine decrements and weak tap survivors; common latency/channel
permutation; and metadata, budgets, cancellation, darkness and replay.
Native formatting and focused no-codegen Rust analysis passed for the exact
new production/test files, existing exposure/readout code and shared temporal
module (28 test definitions). This is not a full BEM crate build or executed
Rust tests. Independent Python delayed-signal quadrature agrees with the
analytical moments to 2.00e-14 maximum absolute discrepancy for Uniform and
2.10e-14 for Triangular across means, mean squares and readout covariance.

## Mixed real/proper BEM source statistics

`HarmonicExposureFactorStatistics::{Real, ProperComplex}` declares the
unit second moments of each fluctuation coordinate. The new explicit
statistics constructors require one law per supplied factor, including empty
factors, and copy that law into the original frequency bank and every prepared
delay-and-sum bank. Ordered factor-statistics iterators expose the retained
declaration. Legacy constructors still make every fluctuation factor proper.

For a factor carrying amplitudes A and B at angular frequencies a and b,
the windowed product contains
`Re(A*conj(B)*Phi(a-b) + kappa*A*B*Phi(a+b))/2`,
where Phi is the normalized negative-time window transform and kappa is one
for Real and zero for ProperComplex. Thus real uncertainty about a coherent
source retains carrier-sum interference; it must not be silently circularized.
This law is used consistently for pressure, all velocity components, signed
p/u flux, physical receiver readouts and delayed processed pressure channels.
The complete-cycle monochromatic observers continue to depend only on the
conjugate source covariance.

The shared temporal kernels and original DD transfer contractions are reused.
There is no additional source/receiver covariance matrix, propagation solve,
sampling, geometry approximation or frequency discretization. Arithmetic
complexity and existing hard work caps are unchanged. Logical heap admission
includes the enlarged bank headers that own the factor laws. A malformed
statistics count, budget refusal or cancellation publishes no partial result.

Six focused regressions cover the real-versus-proper subperiod p/u distinction,
an exact eight-state real/proper ensemble against direct delayed real-time
quadrature under Uniform and Triangular shutters, proper-to-two-real
quadrature equivalence, real orthogonal basis rotations and mixed-factor
permutations, physical clock-epoch covariance, and all-proper compatibility,
identity, workspace admission, darkness and cancellation. The numerical
fixture exercises the actual prepared-observer algebra with synthetic original
transfer rows; it does not claim to validate a boundary solve.

All six new mixed-statistics regressions compiled, linked and executed in a
native Rust executable with the real fs_math library: 6 passed, 0 failed.
The exact exposure/readout/array implementation and current shared window
source were included. The bounded harness supplies nominal outer
source/receiver/medium/hash declarations, the original prepared-observation
layout/admission and the three used error variants; it omits the unused
product-spectrum export and does not compile or solve the complete BEM crate.
The 22 shared window tests were filtered from this run. Native rustfmt and
focused Rust no-codegen analysis also passed.

## Oriented acoustic surface-power spectra

`PreparedHarmonicSurfacePower` contracts the original prepared BEM pressure
and Cartesian particle velocity on an explicitly authored receiver quadrature.
Each `HarmonicSurfaceQuadrature` entry names an original physical receiver
index and a finite oriented area vector [m²]. The expected signed power is
`P(t) = sum_q area_q dot E[p_q(t) u_q(t)]` [W]. Signed and zero vectors
remain explicit; repeated receivers support subdivision or separate patches.
No surface closure, outward orientation, normal interpolation, plane-wave
pressure-to-intensity substitution, dissipation or absorption is inferred.

The optional mean and each real/proper source factor retain their complete
frequency and source coherence. For every ordered frequency pair within a
bank, spatial DD contraction first forms the difference amplitude
`0.5 * sum_q p_a(q) * conj(area_q dot u_b(q))`.
Mean/real banks additionally contribute the sum amplitude
`0.5 * sum_q p_a(q) * (area_q dot u_b(q))`;
proper banks omit that sum under their original zero pseudocovariance.
DD angular sums/differences pass directly to the shared spectrum without
rounding the source clock through f64 Hz. The resulting observable retains DC,
beats, carrier-sum exchange and negative power. It is an expected power waveform,
not a random realization or a fourth-moment power-noise spectrum.

The owner exposes its spectrum, instantaneous power [W], normalized Uniform or
Triangular detector averages [W], and unweighted physical interval energy [J].
Interval energy uses the shared zero-decay causal integral; a DD source-clock
entry point preserves an epoch plus elapsed time. A shaped detector average
is not silently multiplied into physical energy. Original receiver/source
identity, medium owner, producer uncertainty, clock, quadrature and copied
factor statistics survive after the source exposure is dropped. Source,
spatial, quadrature and arbitrary-null relative accuracy remain uncertified.

For Q quadrature entries, F_b retained frequencies in each bank and T generated
terms, preparation costs O(Q*sum_b(F_b²) + T*log(T)); memory is
O(Q + T + factors). Spatial products are reduced before emitting terms, so
there is no receiver-by-frequency-pair spectral expansion. Explicit admission
bounds quadrature count, Q*sum_b(F_b²), generated terms before coalescing, and
peak logical workspace including both temporary terms and retained-tone
capacity. Existing physical fields are borrowed, with no new solve or dense
covariance. There is a hard 65,536-entry quadrature ceiling plus the shared
131,072-term ceiling and original fixed work cap. Empty quadrature yields
explicit zero flux with original physical identity. Cancellation/refusal
publishes no partial surface or observation.

Five focused G0/G1/G3 regressions exercise an exact eight-state mixed real/proper
ensemble against independently evaluated real-time surface flux, agreement
with original finite-window p/u observations, signed standing-wave reactive
exchange, orientation reversal, subdivision/permutation, copied provenance,
and explicit silent entries. A high-carrier proper beat and cancelling oriented
areas test weak flux/energy and frequency-pair storage independent of spatial
sample count. Admission, cancellation and empty/dark surfaces remain atomic.

Native Rust compilation, linking and execution passed all five new surface
regressions: 5 passed, 0 failed. The exact surface implementation/tests,
current exposure implementation and current shared window/spectrum/causal
implementation were exercised with the real fs_math library; the 43 shared
tests already passed separately and were filtered here. The bounded harness
supplies nominal outer identity/medium/hash/error declarations and the original
prepared-observation layout/admission, while omitting untouched array/readout
children. This verifies the physical observation algebra, not a complete BEM
crate build, boundary solve or repository quality gate. Native rustfmt passed.
