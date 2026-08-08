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
4. Repeat solves are bitwise identical.

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
  on the 80-panel fixture (measured, recorded). No FMM acceleration
  (dense cap 8192 panels), no spherical-harmonic directivity coefficient
  tables (sampled directivity only), no scattering/incident-field path,
  no half-space or impedance boundary conditions, no Bessel-backed
  piston closed form (small-ka series only until the duct bead's special
  functions land).
