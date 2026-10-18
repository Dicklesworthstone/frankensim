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
  sign, the Burton-Miller coupling is `alpha = -i/k` (pinned by the
  interior-resonance contrast, not by convention label).
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
  by exact band-limited quadrature (Gauss–Legendre x uniform phi);
  truncation honesty is the reported `captured_fraction`, and
  `evaluate`/`power_by_degree`/`coefficient` serve the rendering and
  diagnostics consumers. Measured at ka = 1 on 320 panels: monopole
  l = 0 fraction and captured fraction > 0.999, dipole (l=1, m=0)
  fraction > 0.99, reconstruction vs directly sampled far field
  < 1e-3 relative L2 at held-out directions.
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
unchanged. Plain CBIE is explicitly refused by this owner.

The actual dense Burton-Miller arm uses centroid off-diagonal kernels and
equivalent-disc self terms, not the triangle quadrature used by Plain CBIE.
The iterative path reuses those exact parent self terms and the same point
kernels in adjacent leaves; it introduces no second singular quadrature.
Static double/mixed row closures use the same effective order. W5a shares its
traversal with the dynamic operator; directional zero-k closures select W5a's
low-frequency partition. The static constant-density identities close under
the computed zero-k operator. Static-sum and cross-route FMM approximation
remain part of the operator error, independently tested against dense LU.

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

Rows use `exp(-i omega t)`: resistance is `Re Z`, added mass is `-Im Z/omega`.
The unchanged matrix must have relative reciprocity defect at most 1e-12 and
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
   compatibility and introduces a net static flux amplified by `alpha=-i/k`.
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
  follow-up). `DirectivityTable` truncation error is reported via
  `captured_fraction`, never bounded a priori; the table inherits every
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
