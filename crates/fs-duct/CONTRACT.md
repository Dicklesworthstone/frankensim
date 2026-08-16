# fs-duct — Contract

## Purpose and layer

Layer L3 (FLUX). 1D viscothermal duct/horn acoustics
(bead frankensim-fsim-duct-acoustics-zdmm1): transfer matrices over
cylinder/cone segment chains with Zwikker–Kosten wide-tube losses,
low-`ka` radiation terminations, and input impedance. Generic by
doctrine: instrument bores, mufflers, HVAC runs and waveguides are the
same object, and every medium property derives from
`fs_material::gas::GasState` — no hardcoded air.

## Public types and semantics

- `Segment::{Cylinder, Cone}` — geometry in SI; a cone is a linear
  taper propagated with exact 1D spherical waves. Lossy cones
  cascade spherical substations at each slice's own mid-radius
  (lossless stays the one-shot `e^{±ikx}/x` 2-port).
- `LossModel::{Lossless, WideTube, AllRegime, Bessel}` — lossless,
  first-order ZK, piecewise wide/Poiseuille, and frequency-by-frequency
  Bessel Zwikker–Kosten `F(r_v)` (every shear number). A cone still
  uses the spherical `e^{\pm ikx}/x` basis.
- `Termination::{Closed, IdealOpen, UnflangedOpen, FlangedOpen}` —
  ideal limits plus the classic low-`ka` radiation fits (end
  corrections 0.6133a / 0.8216a) with a named `ka` ceiling.
- `compact_radiation_impedance` — the radiating-mouth load as a
  named primitive (same fit the TMM already uses).
- `absorbed_spherical_pressure` — free-field observer: spherical
  spreading times ISO 9613-1 molecular absorption. Humidity is an
  explicit `[0, 1]` argument. Do not add Stokes–Kirchhoff on top.
- `segment_wave` — complex `k` and characteristic impedance from the
  gas state: `k = k0 [1 + (1+i)/(sqrt2 rv)(1 + (gamma-1)/sqrt Pr)]`,
  `Zc` with the `(1 - (gamma-1)/sqrt Pr)` factor; shear number
  `rv = r sqrt(rho omega / mu)` reported, refused below 10.
- `input_impedance` / `impedance_sweep` / `impedance_peaks` —
  `Z_in = p/U` at the inlet with per-solve diagnostics (minimum shear
  number, mouth `ka`).
- `Segment::ToneHole` + `HoleState` + `tone_hole_shunt` — side
  branch (`[[1,0],[1/Z_h,1]]`). The chimney is a short cylinder:
  OPEN is that run plus a flanged mouth (Dalmont inner length
  on `b/a`; `0.8216 a` or the Rayleigh piston lives in the
  termination). CLOSED is the same run with a rigid cap. A
  compact chimney reprints lumped `L`/`C`; a long one carries
  its own quarter-wave. `tone_hole_shunt_wall` and
  `input_impedance_wall` put the same `WallPin` on that
  chimney cylinder. WideTube on the neck is Bessel so
  a narrow chimney does not raise the bore's `r_v` refusal
  or jump at the AllRegime `r_v = 10` cliff.
  An OPEN hole is still the T-junction
  `series(Z_s/2)·shunt(Z_h)·series(Z_s/2)` with Nederveen
  `t_s = −0.37 b²/a`. Hole radius must stay below the bore
  radius.
- `DuctError` — stable `FS-DUCT-*` refusals: bad parameter, too-narrow
  (wide-tube floor), radiation-`ka` ceiling, empty duct, singular.

Time convention `e^{-i omega t}` (matching `fs_bem::helmholtz`):
`Im k > 0` decays, mass-like reactance is negative imaginary, closed
pipe `Z_in = +i Zc cot(kL)`.

- `modal` module: the multimodal (m = 0 radial mode) image —
  `mm_input_impedance(duct, state, omega, loss, termination, n_modes,
  extra_slices) -> ModalResponse` with the N x N input impedance matrix,
  the plane (0,0) element, and per-mode local cutoffs at the throat.
  Modes `psi_n = J0(gamma_n rho/R)/|J0(gamma_n)|` with `gamma_n` the
  roots of `J1` found in-crate (double-double power series + bisection —
  no transcribed constants); per-mode `k_n = sqrt(k0^2 - (gamma_n/R)^2)`
  from the LOSSY plane wavenumber, evanescent modes kept decaying, never
  clipped; `Zc_n = Zc_0 k_0/k_n`. The chain recurses a REFLECTION matrix
  (bounded under evanescence) from the mouth backward; junctions couple
  mode sets through closed-form Lommel projection integrals; cones and
  flares STAIRCASE into short cylinders whose density is one arm of the
  convergence-ladder disclosure. Tone holes refuse in this image
  (`BadParameter`); the plane-wave image keeps them (mm-001..mm-006).
  `modal_characteristic_impedances` and
  `modal_reflection_from_impedance` are the public seams for runtime
  consumers (the fs-couple MM characteristic-line realizer), keeping
  the `Im k >= 0` branch rule and the `Zc_n = Zc_0 k_0/k_n`
  derivation single-sourced here.

- `TabulatedLoad` + `input_impedance_load` / `input_impedance_tabulated`
  / `modal::mm_input_impedance_tabulated`: the measured/baked mouth-load
  lane (bead zolja). The TABLE is the source of truth: admission
  enforces strictly increasing positive omegas, finite rows, and
  passivity (`Re Z >= -tol |Z|`); queries linearly interpolate and
  REFUSE outside the support (`FS-DUCT-OUT-OF-TABLE`) — nothing
  extrapolates. THE `ka` LIFT: the low-`ka` analytic fits keep their
  `MAX_RADIATION_KA` refusal exactly as before; only tabulated loads
  play above it, inside their own support. In the modal image the table
  loads the PLANE mode; higher modes keep the disclosed matched-mouth
  closure (tl-001..tl-004).

## Invariants

1. Segment 2-ports are built numerically from exact analytic basis
   solutions (plane/spherical waves with the transmission-line
   relation `U = S p'/(i k z_specific)`) — never transcribed matrices.
   Oracles: lossless closed/open cylinder matches `+i Zc cot(kL)` /
   `-i Zc tan(kL)` to 1e-9; the cone reproduces the cylinder in the
   zero-taper limit; the (p, U) determinant is exactly 1 (constant
   Wronskian, `r` proportional to apex distance); half-cone
   composition equals the whole cone to 1e-9.
2. A near-complete closed-open cone resonates at ALL harmonics
   `n c / (2 L_apex)` (the classic conical-bore result, within 3% on
   the committed fixture).
3. The lossless unflanged quarter-wave ladder lands on
   `(2n-1) c / (4 (L + 0.6133 a))` within 0.5%; the viscothermal
   ladder sits below it by the independently computed dispersion
   deficit (measured ratio 1.006/0.998/1.000 across three peaks).
4. Resonance Q of a closed viscothermal cylinder matches the
   closed-form Kirchhoff wall-loss `Q = k/(2 alpha)` within 3%
   (measured 0.6%) — the loss model validated through the full TMM
   machinery against an independent expression.
5. The thermal loss share for air is the `(gamma-1)/sqrt(Pr)` factor
   (measured 1.475x viscous): a dropped thermal term fails loudly.
6. `Re Z_in >= 0` across the band for radiating viscothermal ducts
   (sign pin); lossy peaks flatten below lossless ones (dispersion
   sign pin).
7. Ambient parameterization: lossless resonances scale with
   `sqrt(T_hot/T_cold)` against the independently computed constant.
8. Repeat evaluations are bitwise identical.
9. MEASURED-DATA VALIDATION: the Ernoult 2021 four-hole cylinder
   (Acta Acustica 5:47, CC-BY Table 1 geometry; measured curves
   published via openwind, GPLv3) — the five-fingering first-peak
   ladder reproduces the measured 283/332/449/619/770 Hz within an
   authored 30-cent envelope (Dalmont inner matching on `b/a`
   removed the old systematic -8..-20 cent flat bias) with the
   monotone fingering-ladder doctrine check; plus the exact-cascade
   tone-hole algebra pin
   (1e-12), open-raises/closed-perturbs contrasts, and hole refusals.
10. (Review round) The Zc/impedance correction is pinned by the
   INDEPENDENT sqrt(Z_series/Y_shunt) transmission-line route (complex
   square roots, physical-branch selection) with the review's three
   surviving eps_z mutants asserted OUTSIDE the second-order band;
   contracting cones match port reversal and the lumped
   cavity-compliance limit (9.5e-4); the flanged ladder lands on
   0.8216a; chained same-radius halves equal the whole through
   input_impedance to 1e-10 in both loss arms.

11. Multimodal gates (music bead 3ez8g.4.1): the N = 1 modal image
    degenerates to the scalar plane-wave path on stepped-cylinder
    chains within 1e-9 relative (measured 1e-10 class); modal
    wavenumbers are the analytic `sqrt(k^2 - (gamma_n/R)^2)` in both
    regimes; the junction projection matrix matches independent in-test
    quadrature to 5e-9; a sudden expansion's evanescent modes add a
    purely reactive, convergent impedance shift; and on a trumpet-like
    flare the plane-wave image misses the multimodal peak structure by
    >= 8 cents (measured 14.9) while the mode ladder is settled at the
    top (N=4 vs N=5 <= 1 cent, measured 0.18) — THE recorded trigger for
    the multimodal expansion, executed (mm-001..mm-005).

12. Tabulated-load gates (bead zolja): admission refuses short,
    unsorted, non-finite, and non-passive tables; interpolation is
    exact at nodes and linear between them; out-of-table refuses on
    both sides; the analytic `UnflangedOpen` still refuses at
    `ka > 1` while a table covering the same frequency plays
    (tl-001/tl-002); and an fs-bem-baked bell table moves a throat
    cylinder's dominant peak by hundreds of cents relative to the bare
    unflanged fit, on the record (tl-003: 1188 -> 852 Hz, -576 cents).
    The committed bake artifact `data/radiation/bell-fixture-zl.tsv`
    is schema-, resolution- (ppw >= 6), and passivity-gated on every
    run; minting is an explicit `--ignored` test.

## Error model

Typed `DuctError` refusals with stable codes; validity is refused by
name (narrow tubes, radiation `ka`), never silently degraded.
Diagnostics (`min_shear_number`, `mouth_ka`) ride on every response.

## Determinism class

Bit-deterministic across runs on a platform: sequential fixed-order
assembly, `fs_math::det` transcendentals, dense 2x2 complex solves via
`fs_la::eigen_complex::lu_complex`.

## Cancellation behavior

Pure short-running functions; sweeps are caller-chunkable per
frequency. No `Cx` integration (tracked with the workspace-wide
`frankensim-ccmn` effort).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

Inline `tests` module: lossless closed forms; quarter-wave end
correction + dispersion deficit; four cone oracles; Kirchhoff-Q;
passivity/flattening/thermal-share pins; hot-duct sqrt(T) law; named
refusals; bitwise determinism.

`src/modal.rs` `modal_tests`/`modal_flare_tests`, cases mm-001..mm-006 —
Bessel self-certification (J0' = -J1, root residuals 1e-14, mode
orthogonality and junction projections vs independent quadrature),
analytic cutoffs both regimes, plane-image degeneracy (stepped chains
1e-10 class; cone staircase converging with slice density), evanescent
junction inertance (reactive, two-step-net convergent — equal-N mode
matching has the classic oscillating truncation tail, disclosed),
refusals + bitwise repeats, and the flare convergence ladder + trigger
(peak tables logged per N and per staircase density; the first
measure-mode run caught a coarse-staircase artifact posing as a treble
mode shift — the trigger is therefore judged with BOTH arms on the fine
staircase).

`src/lib.rs` `tabulated_load_tests` + `bell_bake_artifact`, cases
tl-001..tl-004 and the committed-artifact gate — table admission /
interpolation / refusals, the `ka` lift (analytic refusal and tabulated
success at the same frequency), the fs-bem end-to-end bell bake with the
logged peak shift, and the modal plane-mode table lane. The bake driver
itself is certified in fs-bem (`radiation_bake`, zb-001..zb-003:
pulsating-sphere oracle inside the solver's own tested bands, sweep
refusals including too-coarse-never-extrapolate, and the
divergence-theorem closure/orientation pin on the lathe).

## No-claim boundaries

- Wide-tube first order refuses `rv < 10`. [`LossModel::AllRegime`]
  falls to Poiseuille there; [`LossModel::Bessel`] is the
  frequency-by-frequency Zwikker–Kosten `F(r_v)` at every shear
  number. Both store Helmholtz `k` (`e^{ikx}`), not telegraph `γ`.
- Straight smooth isothermal walls. A [`fs_phs::WallPin`] on
  `input_impedance_wall` adds the locally reacting shunt
  `Y' = 2π a slant / (r − iωσ + i K/ω)` to the gas `Y'`
  (`slant = √(1+(dr/dx)²)`, `1` on a cylinder; same pin
  as the ODE cell shunt). `input_impedance` is rigid.
  No roughness, no porous liners, no mean flow, no nonlinear
  (high amplitude) propagation.
- Lossy cones cascade spherical substations with `k, Zc` at each
  slice's own mid-radius (lossless without a wall stays the
  exact one-shot `e^{±ikx}/x` 2-port; a wall follows the
  local radius and slant, so that path slices too). The multimodal
  expansion for strongly flaring horns now EXISTS (`modal` module;
  trigger executed: >= 8 cents of plane-wave peak error on a
  trumpet-like flare). Its own v1 boundaries: m = 0 modes only
  (axisymmetric bores and sources), no tone holes, no mean flow;
  higher modes terminate into their own characteristic impedance at
  the mouth (a matched-mouth approximation, disclosed — the plane
  mode can now play a TABULATED bell load via
  `mm_input_impedance_tabulated` (bead zolja), while the analytic
  plane load keeps its `ka` refusals; the higher-mode mouth coupling
  remains the disclosed boundary); per-mode viscothermal loss is carried through the
  LOSSY plane wavenumber (`k_n^2 = k_0^2 - (gamma_n/R)^2`), so the
  wide-tube validity boundary binds per mode through the plane
  mode's shear number, and no independent higher-mode boundary-layer
  model is claimed; validation against measured instrument
  impedance remains the brass-gates bead's scope.
- Unflanged radiation is the low-`ka` Levine–Schwinger fit
  (ceiling 1.0, refused by name). A flanged mouth uses that fit
  below `ka = 1` and the Rayleigh baffled piston
  (`fs_phs::baffled_piston_impedance`, same half-space kernel as
  `fs_bem::helmholtz`) above it. A full exterior BEM of an
  unflanged pipe is still the recorded successor.
- Tone holes are the compact T-junction (Nederveen series + Dalmont
  inner + chimney wall law + `Y = σ Y_open + (1−σ) Y_closed` vent
  mix + mutual series `t_m ∝ e^{−s/a}` from neighboring open
  holes). The Ernoult dataset is the regression floor.
- The gas state inherits `fs_material::gas` boundaries (ideal gas,
  phase validity unchecked, calorically perfect gamma).
