# fs-bem CONTRACT

## Purpose and layer

Layer: L3 (FLUX). Laplace BEM panel methods (plan §8.3 [F], bead
tfz.20): potential-flow screening for exterior aerodynamics — the
ornithoid flagship's wide-search stage. INVISCID HONESTY LABELS apply
everywhere: this is screening, not a viscous truth source.

## Public types and semantics

- `panel3d`: validated `SpherePanels` (centroid/normal/area panelization of
  fs-rep-mesh icospheres); `dense_matrix` — the collocation Neumann
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
- `Formulation::{PlainCbie, BurtonMiller, BurtonMillerWrongAlphaSign}` —
  measured roles: PlainCbie is the accurate non-resonant arm (1.7-3.4%
  on the pulsating sphere across ka in [0.05, 5]); BurtonMiller is the
  resonance-safe production arm for ka >= 0.5 (1.1-5.2%; 1.7% at the
  ka = pi fictitious frequency where PlainCbie hits 50%). The wrong-sign
  arm exists so the mutation fails loudly (measured 3.3x worse at
  resonance).
- `radiation_impedance_matrix` — one factorization, n unit-velocity
  solves; feeds the vibroacoustic-coupling bead.
- `far_field` — sampled directivity amplitudes (monopole uniform to
  0.00%, dipole cos-theta correlation 1.0000 measured).
- `exterior_pressure_at_points` — deterministic batched finite-distance
  pressure via `SUM (dG/dn_y p - G i omega rho v) A`; mismatched source
  surface or medium and invalid/non-exterior points are refused.
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
- `baffled_piston_impedance` — Rayleigh-integral piston (half-space),
  validated against the Bessel-free small-ka series.
- `HelmholtzError` — stable `FS-BEM-HELM-*` refusals: bad parameter,
  too-coarse (< 6 panels/wavelength), dense work cap (8192 panels),
  shape mismatch, singular.

## Invariants

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

1. Kernel formulas are pinned by central-finite-difference tests of G;
   disc self terms by numerical quadrature of the regularized kernels.
2. The hypersingular static self entry uses the exact closed-surface
   identity `N_0[1] = 0` as a discrete row sum, so the same point-panel
   quadrature appears on both sides and its error cancels.
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

## No-claim boundaries

- 3D LIFTING surfaces (Kutta strips, wake SHEETS) and the fs-vpm
  pairing for flapping gaits — the 2D shedding loop ships; 3D is the
  flagship successor.
- Exact panel-integral far fields (centroid monopoles ship for
  off-diagonal rows; analytic quadrilateral/triangle integrals are
  follow-up under the same operator surface).
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
- Helmholtz (`helmholtz` module): BurtonMiller's radiation RESISTANCE
  below ka ~ 0.5 is NOT resolved by centroid quadrature (measured -9.4
  vs +1.03 Pa s/m at ka = 0.05 on 320 panels; recorded unasserted in
  the JSON evidence) — use PlainCbie below the first interior resonance;
  exact singular triangle quadrature is the recorded fix trigger.
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
  solid-angle guard and retains centroid-panel error.
