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

`levine_schwinger::{LsQuadrature,LsSolution,levine_schwinger,
levine_schwinger_impedance}` adds the original thin-wall unflanged-pipe
Wiener–Hopf solution (Levine & Schwinger, Phys. Rev. 73, 383–406, p.383 and
Eq.V.16, DOI:10.1103/PhysRev.73.383). `LsQuadrature { nodes }` explicitly
selects 16..=1024 fixed Gauss–Legendre nodes. `LsSolution` retains ka, |R|,
ln|R|, l/a and the realized node count; `reflection()` and
`normalized_impedance()` use exp(-i omega t), `R=-|R| exp(2 i ka l/a)` and
`Z/Zc=(1+R)/(1-R)`. The physical impedance is in Pa s/m3, with
`Zc=rho*c/(pi*a²)`, referred to the geometric mouth plane.

`Termination::UnflangedThinWall { quadrature }` is additive and does not
change any existing termination/default/golden. `LS_MAX_KA` is a lazy value
computed and refined with the shared J1 root finder, not a transcribed
constant. ka at or above it refuses as `FS-DUCT-LS-CUTOFF`; nonfinite/nonpositive
ka or SI parameters and out-of-budget nodes refuse as `FS-DUCT-BAD-PARAMETER`.
The model is for an axisymmetrically excited semi-infinite rigid circular pipe
of zero wall thickness. Finite walls, bells, external obstacles, flow, nonlinear
edge loss, evanescent feedback from a nearby upstream discontinuity and
propagating higher modes need their own model.

The WH tail is over **0..infinity**, with denominator `x sqrt(x²+(ka)²)`.
Its quartic compactification integrates the entire domain; the propagating
integral uses `x=ka sin(theta)`, with `theta=(pi/2)(3u²-2u³)` smoothing both
angular endpoints. `atan2(J1,-Y1)` retains the continuous phase
through the Y1 zero. DLMF 10.8.1/10.31.1 series and 10.32.9 scaled K1 integral
avoid cancellation/overflow; the auxiliary Gaussian integral's u>12 tail is
below 1e-60 for x>=1. At x>=32, the 24-term DLMF 10.40.1/2 expansions provide
scaled I1*K1 without ever forming exp(+x). Quadrature refinement is numerical
evidence, not an interval certificate; accuracy is not uniform arbitrarily
close to the higher-mode cutoff.

The original article's numerical low-frequency constant 0.6133 is historical:
the original integral gives **0.61270103593...**, consistent with Rienstra,
Fundamentals of Duct Acoustics (2015), section 7.2, delta=0.6127:
<https://sjoerdr.win.tue.nl/papers/VKI_Rienstra.pdf>.
Rienstra & Hirschberg, An Introduction to Acoustics, section 7.9 p.214
(6 Sep 2026), also gives delta=0.6127; V.V. Tupov's 2019 recalculation,
AIP Conf. Proc.2195,020034, DOI:10.1063/1.5140134, reports the correction
in its abstract. `ls_historical_fit_retains_its_distinct_end_correction`
pins the historical load coefficient and its 0.00100654 relative discrepancy
from the exact l/a at ka=0.01 without changing that load.
Independent 35-digit mpmath integration gives l/a(0.01)=0.6126826892082684
and |R|(0.05)=0.9987549288169632. Thus the leading loss `(ka)²/2` at 0.05
has relative error 0.00394305, and is not a 1e-3 oracle at that finite ka.
The low-frequency test uses the exact-integral constant and ka=0.01/0.005/0.0025;
at ka=0.05 it also bounds the leading loss expansion's measured relative
remainder between 0.0039 and 0.0040, alongside the exact-reflection oracle.
The old UnflangedOpen fit and its frozen TMM sweep remain unchanged.

Model-form error of the legacy fit `[(ka)²/4 - i*0.6133*ka]` against the WH
normalized load, from independent mpmath Bessels and 256/512-node quadrature:

| ka | relative complex impedance error |
| --- | --- |
| 0.1 | 0.002556512 |
| 0.25 | 0.009532061 |
| 0.5 | 0.028737865 |
| 1.0 | 0.094032364 |

The first 1e-3 crossing is near **ka=0.0113471** (about 0.01135), so the
legacy validity ceiling of ka=1 does not imply a 1e-3 model-form error band.
These are independent mathematical reference calculations, not native build
qualification; `ls_low_ka_fit_model_form_band` recomputes the errors/crossing
on the production path. The L2 fs-phs MouthFlange retains its disclosed low-ka
fit and cannot call this L3 owner; no default consumer is migrated here.

- `Segment::{Cylinder, Cone}` — geometry in SI; a cone is a linear
  taper propagated with exact 1D spherical waves. Lossy cones
  cascade spherical substations at each slice's own mid-radius
  (lossless stays the one-shot `e^{±ikx}/x` 2-port).
- `LossModel::{Lossless, WideTube, Bessel}` — lossless, the retained
  first-order wide-tube benchmark, and exact frequency-by-frequency
  Bessel Zwikker–Kosten `F(r_v)`. The interpolated variant is removed;
  no Bessel refusal falls back to another material model. A cone still
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
  or jump at a narrow/wide material-model boundary.
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

- `FingeringTable` / `Fingering` (music bead 3ez8g.6.1 — the deferred
  fingering-tables slice consumed): typed chart-side fingering data
  over a duct template. The template's tone-hole GEOMETRY never
  changes; a fingering re-states only the hole states, and both the
  TMM and characteristic images consume the same table. Admission:
  template must carry holes, one state per hole per fingering, unique
  labels; an unknown label REFUSES by name (never a silent no-op).

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

- Breath vs hall (music bead 3ez8g.3.5): the emergent demonstration
  that ambient state is physics — one clarinet-scale bore, three
  moist-air states, first-impedance-peak pitch in cents. Cold dry
  hall (10 °C) vs warm saturated breath (37 °C): 94.7 cents; at fixed
  body temperature HUMIDITY ALONE moves the peak 17.5 cents, within
  0.2 cents of the analytic sound-speed-ratio prediction (17.3) — the
  executed duct and the mixture rule agree independently.

## Conformance tests

`src/levine_schwinger.rs::tests` (IRP M2.1-W1): G2 eight independent
35-digit WH reference rows at ka=0.01/0.1/0.5/1/2/3/3.5/3.8 (1e-3 relative),
plus independent fsci-special Bessels and tanh-sinh quadrature (1e-8 relative);
the corrected static/low-ka limits; G1 32/64/128/256-node refinement, final
pair below 1e-10 and informative coarse rung; G0 cutoff and node cap/cap+1
refusals, passivity and negative mass-like reactance; G3 measured legacy-fit
model-form band; G5 bitwise repeat. `tests/duct_conformance.rs::dt_006_levine_schwinger_thin_wall_termination`
checks the lossless 0.30 m x 8 mm cylinder against an independent closed-form
load transfer at 50 frequencies through ka=3.8 (1e-10 relative), with per-check
JSONL receipts. The existing cross-ISA sweep golden must remain unchanged.

`tests/duct_conformance.rs` (bead 3ez8g.13.2; the reimplementation
contract through the public surface, independent oracles): dt-001
the closed lossless cylinder matches `Z_in = +i Z0 cot(kL)` under
e^{-i omega t} to 1e-16 (the -i form belongs to e^{+i omega t});
dt-002 the vanishing-taper cone degenerates to the cylinder (8e-6)
and the ideal-open pipe peaks at c/4L (2e-4); dt-003 wide-tube Q
follows the sqrt-frequency law across the first two resonances
(0.3%); dt-004 refusals by name (narrow-tube shear floor,
non-physical parameters, empty duct); dt-005 the unflanged
(0.6133 a) and flanged (0.8216 a) end corrections move the
quarter-wave peak by the textbook effective-length ratios (3e-4).
`tests/ernoult_corpus.rs` carries the registered corpus fixture.
Inline unit modules stay where they are cheap; nothing was weakened.

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

- Wide-tube first order refuses `rv < 10`. `LossModel::Bessel` uses
  the shared fs-phs exact law on its admitted shear range [1e-6,1e6],
  with cancellation-free `1-F` below |ζ|=12. Invalid/out-of-domain data
  refuse; there is no interpolated narrow branch or silent zero-loss result.
  Both store Helmholtz `k` (`e^{ikx}`), not telegraph `γ`.
  G2 `bessel_small_rv_series_matches_poiseuille_limit` checks derived
  finite-rv R/L/C corrections at 1e-10 relative with explicit omitted-term
  bounds. The leading Poiseuille/isothermal limits alone differ by more
  than 1e-10 at rv=.1. `bessel_refuses_instead_of_falling_back` covers
  bad inputs; `narrow_cone_bessel_no_fallback` requires decreasing
  12/24/48-slice errors against a 96-slice A coefficient, order 1.5..3
  and final error <1%. The specified 0.5..2mm cone at 50Hz has rv approximately 2..10
  in room air, below the wide-tube floor rather than below one.
  Native qualification, before/after peak shifts and golden replay remain
  batch-pending. No continuous-profile/3-D convergence claim is made.
- Straight smooth isothermal walls. A [`fs_phs::WallPin`] on
  `input_impedance_wall` adds the locally reacting shunt
  `Y' = 2π a slant / (r − iωσ + i K/ω)` to the gas `Y'`
  (`slant = √(1+(dr/dx)²)`, `1` on a cylinder; same pin
  as the ODE cell shunt). `input_impedance` is rigid.
  No roughness, no porous liners, no mean flow (QUIESCENT medium,
  `M̄ = 0`; through-flow absent by construction), no nonlinear
  (high-amplitude) propagation. AMPLITUDE BOUNDARY: these claims hold
  for peak local acoustic Mach < 0.01; beyond that they shrink or
  escalate one rung per the falsifier in
  `docs/DECISION_MUSIC_MEAN_FLOW.md` (rung A; the `k/(1-M²)` hack is
  NEVER-listed there).
  slice's own mid-radius (lossless without a wall stays the
  exact one-shot `e^{±ikx}/x` 2-port; a wall follows the
  local radius and slant, so that path slices too). The multimodal
  expansion for strongly flaring horns now EXISTS (`modal` module;
  trigger executed: >= 8 cents of plane-wave peak error on a
  trumpet-like flare). Its own v1 boundaries: m = 0 modes only
  (axisymmetric bores and sources), no tone holes, no mean flow (the
  same quiescent/amplitude boundary as the plane path above; escalation
  per `docs/DECISION_MUSIC_MEAN_FLOW.md`);
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
  below `ka = 1` and the shared closed-form Bessel–Struve baffled piston
  (`fs_phs::baffled_piston_impedance`, same half-space kernel as
  `fs_bem::helmholtz`) above it. A full exterior BEM of an
  unflanged pipe is still the recorded successor.
- Tone holes are the compact T-junction (Nederveen series + Dalmont
  inner + chimney wall law + `Y = σ Y_open + (1−σ) Y_closed` vent
  mix + mutual series `t_m ∝ e^{−s/a}` from neighboring open
  holes). The Ernoult dataset is the regression floor.
- The gas state inherits `fs_material::gas` boundaries (ideal gas,
  phase validity unchecked, calorically perfect gamma).
