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

### Smooth plane-wave profiles (IRP M2.1-W3a)

`profile::RadiusProfile` admits SI data through `exponential(r0,area_flare,L)`
(`S=S0 exp(area_flare*x)`), `power(r0,x0,exponent,L)`
(`r=r0 ((x0+x)/x0)^exponent`) and `spline(stations,endpoint_slopes)`.
Spline stations remain accessible unchanged; None chooses natural ends and
Some([left,right]) supplies dr/dx for a clamped C2 cubic. Every cubic span
must have a strictly positive Bernstein control hull. This is conservative
positivity admission, not clipping or a certificate of specimen uncertainty.
The measured-data owner (e.g. fs-scenario's BoreProfile::stations()) retains
source records and uncertainties; this acoustic image consumes its SI values.

`Segment::Profile { profile: Arc<RadiusProfile>, steps: RiccatiSteps }` embeds
that data in the existing Duct chain. Segment is Clone, no longer Copy; no
second duct type or global profile registry exists. `Segment::axial_length`
also includes profiles. The profile carries no material defaults: propagation
takes the caller's GasState, frequency, loss model, mouth load and step axis.

`profile_matrix` integrates [p,Z_ref U] backward with fourth-order, two-node
Magnus and a traceless 2x2 exponential. Forward off-diagonal coefficients are
`i k Zc/Z_ref` and `i k Z_ref/Zc` under exp(-i omega t), from the existing
`segment_wave` at BOTH quadrature radii. Exact Bessel ZK therefore has one
constitutive owner. `profile_impedance` returns inlet p/U and mouth-pressure/
inlet-flow and mouth-flow/inlet-flow ratios; (+infinity,0) is a rigid load.
Transfer evolution crosses scalar impedance poles without dividing by p or U.
Extreme transfer conditioning/overflow still refuses; arbitrary long strongly
evanescent chains need the separate scattering-composition work item.

`RiccatiSteps` and `ConeSlices::{Auto,Fixed}` admit 1..=65536. Auto preserves
the old minimum-two, relative-radius rule whenever it was <=12, removes its
upper clamp and refuses above the named cap. `cone_matrix` returns the actual
slice count. Existing 12-slice narrow-cone comparison now requests Fixed(12)
explicitly; Auto is a different, live axis. Cylinder/tone-hole golden arithmetic
is unchanged. Profile wall pins and the multimodal profile image refuse by
name: this module converges to Webster's plane-wave equation, not to 3-D.

Determinism is fixed-order sequential det:: arithmetic. `profile_matrix_controlled`
polls a caller cancellation closure before each step and returns
FS-DUCT-CANCELLED with the completed count, publishing no partial transfer.
The ordinary profile/chain calls are synchronous, with no ambient Cx claim.

G2 profile tests compare 100-frequency exponential and order-3/2 Bessel-horn
impedances to independent analytic solutions at1e-9, including exponential
cutoff/evanescence. Constant cylinders (lossless/Bessel) and lossless cones
match existing exact bases at1e-12. Finite lossy-cone coefficients vary with
radius, so that image is checked against the continuous limit with the cone
slice ladder, not asserted equal to one frozen-radius cone. G1 fs-mms gates
Magnus order4, clamped-spline station order4 and spherical-slice order2 within
0.2; mandatory coarse rungs remain visible and finer slices reach1e-8.
`duct_conformance::dt_007_smooth_profile_chain` checks the existing chain with
a TabulatedLoad against independent ABCD impedance/peak positions at1e-9.
G0 cap/cap+1, invalid data, cancellation and G5 bitwise-repeat tests accompany
these claims. Native verification is required before promotion.

Coordinator ruling2858 corrects two incompatible literal bead criteria:
closed-form exponential and lossless conical two-ports satisfy1e-12; finite
lossy cones instead converge under the2±0.2 gate, with a bitwise check of
every midpoint material sample against the shared exact-ZK owner and of the
resulting64-slice cascade. The4/8/16/32/64 rungs remain, extended through4096
so the actual finest pair meets1e-8. No acceptance band is widened.

### Continuous material fields on the same profile propagator

`profile::gas_field::GasFieldDuct` binds the existing `Duct` geometry to
`fs_material::gas_field::AxialGasField`. The latter remains the only primitive
field/source owner. For a varying field, the acoustic binding evaluates the
NASA-9/material state at the actual two Gaussian quadrature points of each
Magnus step, then retains those admitted states for frequency sweeps. Derived
density, sound speed, viscosity, conductivity, gamma and Prandtl number are
never interpolated. Both constant-gas profiles and fields call the same local
generator, exact-ZK law and traceless Magnus exponential.

Cylinders, cones and smooth radius profiles compose in one pressure/volume-flow
chain. Authored geometry boundaries, primitive field stations and profile
spline knots split smooth intervals before integration; derivative jumps do
not pass through an allegedly fourth-order step. `RiccatiSteps(n)` requests n
steps per such interval in the varying-field image. Its explicit field axis
replaces the original uniform-profile integration axis, and the actual total
is available as `steps_used()`, capped at `MAX_RICCATI_STEPS` without clamping.
A uniform field performs no field discretization and delegates the original
duct and its original profile/cone axes bit-for-bit.

Each resolved step becomes a `compose::TwoPort`; the existing stable chain
owner produces impedance, outlet flow/pressure transfer, and direct input
reflection. A frequency sweep never needs a global growing ABCD product.
The exact inlet state supplies the port impedance; the default radiation
medium is the exact outlet state. An explicit exterior must be coherent and
at the same static pressure. All loads retain their existing validity checks.
Travel delay integrates local 1/c with the same quadrature. Cancellation during
material binding or block construction returns the completed-step count and
publishes no partial result.

The `fs-couple::RegionalDuct::from_continuous_gas_field` adapter sends this
operator into its existing characteristic-line and baffled-outlet paths;
`ReedBoreVoice::new_with_continuous_gas_field` drives the original reed loop
with it. The older midpoint `from_gas_field` route remains an explicit
piecewise-constant comparison rung. Scheduled `ReedPerformance` selects
`GasFieldResolution::Continuous(RiccatiSteps(n))` through
`from_bytes_with_field_resolution`; its existing schedule, massive reed,
terminal loading and baffled microphone share that binding. The existing
midpoint constructor and input records retain their behavior.
Nonuniform fields with tone holes or
deformable-wall pins refuse: no spatial chimney/mutual-hole/wall law is inferred.

Verification covers exact uniform bits, independent analytic linear-temperature
cylinder transfer and closed/open impedance (a convergent Frobenius form of
the J0/Y0 solution), fourth-order joint temperature/composition/radius
refinement with Bessel losses, local primitive-state evaluation, pressure and
mixture changes reaching physical transfer, and real rendered-pressure changes
above their step-refinement discrepancy. Source coefficients in these tests
are declared analytic fixtures, not experimental material validation. The
operator remains a frozen, rigid-wall, no-mean-flow, plane-wave model; it does
not generate breath fields, advect species, evolve wall temperature, model
condensation or qualify a complete instrument.

### Segment chains and terminations

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
  (ceiling 1.0, refused by name). A flanged mouth, a bell or an open
  tone hole, is ONE multimodal flanged-pipe model at every `ka`
  (`flanged_pipe_impedance`, Zorumski 1973). The aperture radiation
  matrix over `FLANGED_MODES = 16` axisymmetric duct modes comes from the
  King spectral integral, and evanescent modes are eliminated by Schur
  complement. One mode is exactly the closed-form piston
  (`fs_phs::baffled_piston_impedance`). Measured ladder at `ka = 1e-3`:
  end correction 0.8488 → 0.8330 → 0.8263 → 0.8236 → 0.8224 → 0.8221 a
  for 1/2/4/8/16/24 modes, toward the flanged-pipe 0.8216 a (Norris &
  Sheng 1989).
  - The previous fit-below / piston-above switch jumped by ~0.19 z0 at
    `ka = 1`, which no causal characteristic line can realize.
  - Refused above `FLANGED_MAX_KA = 20`.
  - The plane-mode TMM above a bore's first cross-mode cutoff stays the
    disclosed 1-D boundary.
  - A full exterior BEM of an unflanged pipe is still the recorded successor.
- Tone holes are the compact T-junction (Nederveen series + Dalmont
  inner + chimney wall law + `Y = σ Y_open + (1−σ) Y_closed` vent
  mix + mutual series `t_m ∝ e^{−s/a}` from neighboring open
  holes). The Ernoult dataset is the regression floor.
- The gas state inherits `fs_material::gas` boundaries (ideal gas,
  phase validity unchecked, calorically perfect gamma).

## Generic frequency-domain networks (IRP M2.1-W8)

`network::{DuctNetwork,NetworkNode,NetworkEdge,network_port_impedance}` solves
connected pressure/volume-flow graphs, including parallel paths and self-loops.
Every edge owns a complete `Duct` chain, oriented input->output; nodes are ideal
common-pressure junctions, explicit loads or labelled external ports. Currents
are positive into the network. `NetworkResponse::z` is row-major `p=Z U`, in
Pa s/m³ under exp(-i omega t). The port list follows caller node-array order.

Stable node/edge IDs are required. Assembly and LU use canonical integer indices
sorted by those identities, not a hash-map traversal; reindexing nodes changes
only the returned port permutation. Each edge consumes the existing chain owner
at two independent algebraic boundary columns (+/- a local rho*c/S reference).
This recovers its ABCD block without duplicating propagation/loss physics or
altering any old chain/default/golden. Negative reference columns are algebraic
probes, not active boundary loads in the solved network. The reciprocal
Wronskian det(M)=1 is checked, then used to stamp equal off-diagonals.

`JunctionCorrections::ChainLocal` retains existing compact-hole corrections.
`Disabled` retains each complete chimney and mouth, but removes its compact
inner/series/mutual junction corrections. Bare graph junctions are ideal;
corrections are never added implicitly. `JunctionMultiport` instead supplies
an explicitly identified, frequency-specific plane-port admittance, with exact
bindings to every incident edge end. Its terminal pressures remain distinct.
The supplied matrix must be reciprocal and have PSD real Hermitian part within
1e-12 normalized numerical admission. This is not evidence that a 3-D/modal
producer's geometry, discretization or evanescent-to-plane projection is valid.
No such producer or local junction correction is synthesized by this module.

Caps are 128 graph nodes, 256 edges, 32 labelled ports, 4096 aggregate segments
and 256 pressure unknowns after feature expansion. Structural admission runs
before dense allocation. Disconnected graphs, invalid indices, duplicate IDs,
zero axial-length edges and absent ports refuse. Edge normalized ABCD, edge
admittance and the grounded nodal solve have full infinity-norm condition
estimates; values above 1e10 (or nonfinite/singular blocks) return `Singular`.
The stable scattering successor M2.1-W9 is required outside this envelope.
Accepted edge estimates and the nodal estimate accompany the response. These
are floating-point solves, not interval accuracy certificates for every graph.

Loads are explicit closed, pressure-release, passive finite impedances,
tabulated loads, or an explicit radiation model/radius. Exact thin-wall LS is
consumed through its existing owner and explicit quadrature; legacy fits keep
their own refusals. `NetworkResponse::drive` / `driving_point` require exactly
one assignment for EVERY other external port, and the same `GasState` as the
network solve. They return pressures and currents as well as input impedance;
there is no omitted-port load, medium, wall or clock fallback. Uniform gas is
explicit per solve; axial-field edges and time-domain characteristic networks
are not claimed. Calls are synchronous and bounded; frequency sweeps can be
caller-chunked, with the same existing no-Cx boundary.

Unit tests in `src/network/tests.rs`: G2 T-branch parallel admittance, lossless
and Bessel (1e-12); G2 200-frequency Herschel–Quincke TL (1e-10 dB), independently
re-derived from Selamet, Dickey & Novak, JASA 96,3177 (1994), Eq.(8)
<https://mae.osu.edu/sites/default/files/2021-12/J11_0.pdf>. The bead's 1997
attribution differs from this primary paper. Finite TL notches are included;
an infinite-TL pole is checked as zero transmission within 1e-10, without a
finite dB clamp. G3 long-chimney ToneHole equivalence with corrections disabled
(1e-10), G0 reciprocity/PSD at 100 fs-rand frequencies (1e-12), cap/cap+1 and
invalid/singular cases, G3 exact node-reindexing and G5 bitwise repeat. Additional
tests cover self-loops, supplied-junction stamps and actual exact-LS consumption.
The existing `cross_isa_golden` chain is unchanged and must still pass centrally.
The new tests emit stderr JSON receipts. Numeric claims await native execution.

## Stable composition (IRP M2.1-W9)

`compose` derives local power-wave two-ports from the existing cylinder/cone
bases and constitutive `SegmentWave`. Forward and backward amplitudes are
referenced at opposite endpoints, so only `exp(i*k*L)` with Im(k)>=0 is
evaluated. Spherical endpoint slopes and area factors are retained. Lossy
cones keep their existing midpoint slice axis. Smooth profiles expose the
same local Magnus steps through a shared callback; the stable chain never
forms their global ABCD product. `profile_matrix` retains its explicit matrix
export and original operation order inside its representable envelope.

`compose_chain` offers load-to-input reflection recursion and an independently
ordered Redheffer cross-check. Mouth transfers accumulate local forward-wave
ratios, with log magnitudes and fixed-order unit phases to avoid intermediate
product underflow. A final unrepresentable attenuated tail rounds to zero;
this is not an interval enclosure of subnormal amplitudes. Feedback ratios may
exceed one in a standing wave or at an area change: that physical gain is not
clipped. Cancellation estimates for the internal denominators and input-flow
normalization must not exceed the declared bound, inclusively. Above it the
call returns `Singular`, with no partial/nonfinite physical response published.
The default bound is 1e12; block/mode caps are 131072/32. A strongly damped
matched line remains well-conditioned even when its ABCD chart is unusable.
These are numerical sensitivity estimates, not certified forward-error bounds.

The ordinary `input_impedance_core` first evaluates this stable rule, including
the pole guard. To preserve the existing cross-ISA/reed bit semantics, it may
then use the EXACT historical ABCD operation graph only when the unloaded
transfer-condition indicator `4/|S21|^2` is <=1e8 and that evaluation is finite.
`DuctResponse::composition` records `LegacyAbcd` or the stable route and
`condition_estimate` reports terminated feedback sensitivity. An unrepresentable
ABCD indicator is `None`, not a fabricated finite bound. Outside this numerical
compatibility envelope the raw product is not evaluated. No golden constant or
baseline is re-pinned; the compatibility feedback estimate must also be <=1e8.
Debug/release cross-ISA proof remains a native gate under
`docs/GOLDEN_POLICY.md`. The driving-point LineLadder is not changed here.

`ChainResult` includes the physical closed-cap pressure transfer. The historical
`DuctResponse::p_mouth_over_u_in==0` closed-mouth radiation sentinel is retained,
as required by its existing contract/tests. Both return zero closed-cap flow.
The positive real power-wave reference must match at connected ports. Supplied
feature scattering matrices own their passivity, units and evanescent-basis
qualification; algebraic composition does not mint that evidence. Generic
equal-mode Redheffer products use the existing complex LU; modal diagonal
reflection propagation calls the same bounded primitive with its exact prior
multiplication order. No second factorization or loss model is introduced.

Network assembly now consumes the finite admittance `(I-S)(I+S)^-1/R` from this
same owner, including long opaque edges. `network_chain_two_port` remains an
explicit finite, conditioned ABCD export; graph assembly no longer depends on
its envelope or reconstructs ABCD through opposite signed load probes. Existing
network nodal/admittance condition caps remain 1e10. Local compact corrections
remain explicit and disabled means no inner/series/mutual correction.

Controlled composition polls cancellation between blocks, profile/cone steps,
matrix rows and existing LU rows. Default APIs remain synchronous with an
uncancelled caller; supplied external matrix entries are not themselves a
convergence certificate. All new arithmetic uses deterministic shared owners;
there is no stochastic seed, material default or instrument-shaped solver.

`src/compose/tests.rs` contains G2 1e-12 closed-form/DD fixtures: the 10m, 1mm
Bessel tube at 20kHz (and 100/1000m overflow extensions), a declared 50-segment
18kHz staircase stop-band horn and a 24-closed-hole 15kHz lattice. Their first
native receipts report the old response defect/refusal and determinant digit
loss alongside the independent DD result. Further gates cover 1e-13 legacy
overlap across cylinder/cone/profile/hole/wall data, exact cylinder transfers,
the two rules, condition/block/mode cap boundaries, cancellation, multiport
unit-column solutions, opaque network admittance and G5 repeats. No previous
tolerance, assert, golden or line-ladder baseline is weakened. Native proof,
wasm compilation, causality and reed rerender remain pending centrally.
