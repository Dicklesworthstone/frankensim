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
iterative owner's inherited kernels. Plain-CBIE triangle integration remains
real-frequency only. Dense/resolution/batch caps still apply at |k|. A known
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
  `dcl_dalpha_adjoint` — one transposed solve for the solution
  sensitivity plus solve-free finite-difference output partials, FD-gated
  and explicitly not claimed as an exact symbolic derivative. The
  constant-panel integrals carry a battery-pinned lesson: the normal
  component is (θ₂−θ₁)/2π — the reversed order self-cancels a closed
  sheet's field (caught by the single-panel-vs-quadrature and
  uniform-sheet probes). `solve_naca0012_prestall` is the narrower G2
  validation entry point: it constructs the unit-chord NACA 0012 section and
  refuses non-finite angles or `|alpha| > 10 degrees` before allocation so an
  inviscid solve cannot be mistaken for a stall prediction.
- `wake2d`: `WakeSim` — impulsive-start free wake; Kelvin-conserving
  trailing-edge shedding, regularized point-vortex convection, the
  quasi-steady bound circulation relaxing against wake downwash;
  ledgered traces.

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
in the body table's top thickness. It reports their ratio to the body
example's equivalent-ellipse neck and gates only convergence.

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
existing real/damped `WaveNumber` seam. Complex-frequency continuation uses the
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
