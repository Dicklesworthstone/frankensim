# fs-duct — Contract

## Purpose and layer

Layer L3 (FLUX). 1D viscothermal duct/horn acoustics
(bead frankensim-fsim-duct-acoustics-zdmm1): transfer matrices over
cylinder/cone segment chains with Zwikker–Kosten wide-tube losses,
low-`ka` radiation terminations, and input impedance. Generic by
doctrine: instrument bores, mufflers, HVAC runs and waveguides are the
same object, and every medium property derives from
`fs_material::gas::GasState` — no hardcoded air.

### Independent exterior gas

`exterior::{input_impedance,input_impedance_and_holes,input_reflection,input_impedance_load,
input_impedance_load_and_holes}` use a separately admitted `GasState` for the
terminal radiation and every open/vented side-hole mouth. The finite chimney
keeps its local interior gas, viscothermal transport and material wall; its
pressure/volume-flow port receives the exterior radiation impedance directly.
Closed chimney compliance and internal junction/mutual corrections remain
internal. The existing low-ka vortex reference estimate converts exterior end
masses to interior equivalent lengths with `rho_exterior/rho_interior`.
Identical media retain the original arithmetic. Different hole media cannot
fall back to the old single-medium ABCD path.

The flanged aperture's radiation matrix and piston term use exterior density
and sound speed. Its eliminated duct modes retain the local density, sound
speed and axial wavenumber. `flanged_pipe_impedance_with_exterior` therefore
uses `(rho_i c_i)/(rho_e c_e) K_i/sqrt(K_i^2-alpha_n^2)` for the modal diagonal
in exterior-normalized coordinates. Both media retain the existing `ka` cap;
identical states reproduce the original flanged load. A dimensional two-mode
continuity oracle straddles the internal cutoff, while a double-double
T-junction cascade checks flanged terminal loading, input reflection and open/
vented chimney-mouth volume flow. The existing unflanged approximations retain
their prior scope.

`GasFieldDuct::chain_with_wall_exterior_controlled` carries the same explicit
exterior through continuously varying fields, including local chimneys.
Response and reflection use it together. An omitted exterior preserves the
historical local-gas hole convention; every explicit gas must be coherent and
at the common static pressure. These are frozen acoustic interfaces, with no
mean flow, species diffusion, condensation or varying-medium mutual-hole law.
G1 tests compare independent finite-line ABCD assembly and chimney/outlet flow
continuity; G0/G3 cover identical-gas bits, closed-hole isolation, pressure
scaling, invalid media, and a warm field feeding a complete hole cluster.
`GasFieldDuct::response_with_wall_and_holes` retains open/vented chimney flows
from their actual locations in that same stable cascade. Uniform and varying
fields share one chimney/series-branch preparation and flow-reconstruction
owner. `open_hole_sites` exposes their physical geometry and the original
quadrature travel times, including local chimney travel. This feeds the
`fs-couple` bell-and-hole observer without dropping sources when a continuous
field is selected. A separate warm-prefix/uniform-cluster cascade checks both
chimney entrance and mouth transfers; closed holes never radiate. Existing
frozen-field, mutual-hole and infinite-baffle observation limits still apply.

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
is unchanged. `profile_matrix_wall` and ordinary `input_impedance_wall` chains
apply the existing locally reacting `WallPin` at both quadrature points, with
the actual local radius and `sqrt(1+(dr/dx)^2)` surface-area factor. Per-area
mass/stiffness/resistance are explicitly supplied and constant over this wall
image; no material or shell parameters are inferred from the gas or profile.
`WallLaw` is the shared local specific-impedance seam for source-owned material
walls. `WallModel::{Rigid,Pin,Law}` keeps the existing constant-pin API and
arithmetic while `input_impedance_with_wall`, stable composition and smooth/
continuous-field methods pass a complete frequency law through every actual
cylinder, cone slice, profile quadrature point and admitted chimney radius.
Geometry supplies the same circumference and slant factor; gas and wall
admittance still enter the one original telegraph-to-wave conversion. The law
must return a passive, finite, nonzero impedance under exp(-i omega t) and
refuse unsupported material, frequency or geometry. General laws stay on the
stable cascade, never a constant-pin legacy reconstruction. The source-bound
thin-membrane/Maxwell adapter and original material receipts live in fs-couple;
fs-duct does not duplicate that constitutive law or infer shell physics.
The multimodal profile image remains separate: this module converges to
Webster's plane-wave equation, not to 3-D.

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
Tone-hole junctions retain their exact authored station and use its resolved
gas in the original junction and homogeneous-chimney law. Closed holes may
lie anywhere in the field. Multiple open/vented holes require a constant
primitive temperature/composition state throughout the span between them:
the original mutual-hole correction is retained in full. A single open hole
may lie in a gradient. A varying inter-hole span refuses, including equal
endpoint states with an intervening gradient; no effective medium or omitted
mutual term is inferred. The chimney is explicitly the local homogeneous
model, not transverse heat/species transport. Junction blocks consume the
same work budget and cancellation checkpoints as propagation blocks.
`chain_wall_controlled`, `response_wall` and `reflection_wall`
apply the SAME local wall-admittance owner at each retained gas/geometry
quadrature point. Regional responses, characteristic feedback and baffled
outlet observations all pass through this wall-aware path. Uniform fields
delegate the original wall-aware chain exactly. Missing or invalid wall data
refuses; it does not select rigid walls or air defaults.

Verification covers exact uniform bits, independent analytic linear-temperature
cylinder transfer and closed/open impedance (a convergent Frobenius form of
the J0/Y0 solution), fourth-order joint temperature/composition/radius
refinement with Bessel losses, local primitive-state evaluation, pressure and
mixture changes reaching physical transfer, and real rendered-pressure changes
above their step-refinement discrepancy. Source coefficients in these tests
are declared analytic fixtures, not experimental material validation.
Additional G1/G3 tests check an independent elastic-cylinder compliance
solution, a tapered wall against spherical-slice refinement, fourth-order
joint material-field/wall convergence, passive reflection, and actual feedback
and observer-pressure changes. These prove numerical model composition, not
experimental wall validation. The wall is a locally reacting supplied law;
axial shell modes, thermal evolution and nonlocal wall radiation are absent.
The operator remains a frozen, no-mean-flow, plane-wave model; it does
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
  It needs an `fs_material::iso9613::TerrestrialAir` admission of the
  state's own composition. Anything else is refused by
  `DuctError::NonTerrestrialAir`, which names the species.
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
  `mm_input_impedance_load` supplies an explicit finite plane-mode impedance
  with nonnegative resistance to that same recurrence. It retains every
  interior modal coupling and the matched higher-mode mouth closure. Runtime
  consumers use it for an independently declared DC limit without extending a
  positive-frequency table. It establishes neither causality of an arbitrary
  supplied impedance law nor experimental validity of a radiation bake.

- `modal::mm_profile_impedance(profile, state, omega, loss, load, modes,
  steps) -> ModalResponse` (IRP M2.1-W3b): CONTINUOUS multimodal m = 0
  propagation on a smooth `RadiusProfile`, after Pagneux, Amir &
  Kergomard (1996, Part I). It integrates the coupled modal equations
  `P' = (R'/R) A P + Z' U`, `U' = Y' P - (R'/R) A^T U` from the mouth to
  the inlet with no staircase. The coupling is closed form:
  `A_mn = 2 gamma_n^2 s_m s_n / (gamma_n^2 - gamma_m^2)`, `s_n = sgn J0(gamma_n)`.
  It is derived in the module doc from the Bessel ODE and orthogonality,
  with no quadrature.
  - Each fixed step is a fourth-order two-node Gauss-Magnus exponential
    of the `2N x 2N` generator, applied to the scaled impedance as a
    Moebius map, so evanescent modes relax onto their decaying branch
    without explicit Riccati marching.
  - `MAX_EVANESCENT_STEP = 8` refuses a step whose `max_n |Im k_n| |h|`
    would put the O(1) propagating entries below accuracy. The error is
    "refine Riccati steps"; it never silently degrades.
  - `ModeCount` is admitted up to `certified_mode_cap()`, which is
    computed as 1 plus the number of `J1` roots below
    `BESSEL_CERTIFIED_X = 64` (21 modes). `RiccatiSteps` uses the W3a caps.
  - `MouthLoad` selects the plane-mode load: `Analytic` (the scalar
    termination with its refusals), `Tabulated`, or `Matched` (Pagneux's
    infinite-pipe end). Higher modes stay matched, the staircase's
    disclosed closure.
  - Losses follow the staircase rule: `segment_wave` for the plane mode,
    `k_n` from the lossy `k_0`, `Zc_n = Zc_0 k_0/k_n`. A uniform section
    reproduces both exactly.
  - `ModalResponse` gains `steps_used` and, for the continuous image,
    `tail_bound_estimate = |Z00(N) - Z00(N/2)|`.
  - `RadiusProfile::slope` supplies the analytic `R'` for every shape.
  - `bessel_j0` keeps its double-double series below 26, so every
    staircase argument keeps its bits. Above 26 it switches to the
    generated Hankel expansion.
  - The staircase image, its goldens and `MAX_MODES = 8` are unchanged;
    it is now the O(dx) comparison rung.
- `modal::mm_profile_field(profile, state, omega, loss, load, modes, steps,
  inlet_flow) -> ModalField` (IRP M2.1-W3b): the continuous field behind
  that impedance, driven by a prescribed inlet modal flow (a plane piston
  of volume flow `q` is `[q, 0, ..]`).
  - It keeps the scaled impedance at every step station and each step's
    inverse flow map. It then carries the inlet flow back toward the mouth
    (`u_j = Y_j^-1 u_{j+1}`, which contracts evanescent components) and
    sets `P = Z U` station by station.
  - `ModalField` holds the stations (ascending from the inlet), `P_n`, `U_n`
    and the axis pressure `sum_n P_n psi_n(0)`, `psi_n(0) = 1/|J0(gamma_n)|`.
  - It refuses an inlet flow that is not one finite value per retained
    mode, and a mouth with no finite plane-mode impedance. Memory is
    `O(steps N^2)`. The impedance path is bit-identical.

- `profiled_tone_hole_shunt` / `profiled_tone_hole_junction` (IRP M5.5-W2
  seam qualification): a tone hole whose chimney is a chain of duct sections
  (cylinders, cones, profiles) from the bore outward, for example an undercut
  cone under a cylinder. There is no new `Segment` variant.
  - The chimney is the plane-wave chain of its sections, terminated flanged
    open, rigid closed, or as the vent admittance blend.
  - The junction keeps the cylindrical hole's Dalmont inner length (it
    extends a cylindrical first section, otherwise a cylinder of the
    bore-side radius is prepended) and the Nederveen series length, at the
    bore-side radius. There is no mutual term.
  - One cylindrical section is `tone_hole_shunt` bit for bit.
  - Refuses: no sections, a tone-hole section, an invalid section, or a
    bore-side radius not below the bore.
  - `PROFILED_CHIMNEY_LIMITATION` states the 1-D model's omission of the
    junction's evanescent near field. Comparisons through it carry that
    limitation, so they never authorize omitting a real undercut; the 3-D
    multiport (M2.4-W2) is the acceptance route.
  - Tested in `fs-couple/tests/simplification_equivalence.rs`
    (`profiled_chimney_reduces_to_cylindrical_tone_hole`): the shunt is bit
    identical, a constant-radius profile agrees to 8.2e-14, and the 2-port
    agrees with the `Segment::ToneHole` chain to 1.9e-15. Refusals are tested
    there too.

- `dg` (IRP M2.2-W3a): the lossless variable-area plane-wave system
  `(A/ρc²) p_t + q_x = 0`, `(ρ/A) q_t + p_x = 0` in the time domain, by
  nodal DG on Legendre–Gauss–Lobatto elements of order `p`.
  - `DgDuct::new(duct, gas, DgMesh { elements, order }, inlet, outlet)`
    meshes `Cylinder`, `Cone` and `Profile` segments (at least one element
    each). Its mass matrices integrate the exact `A(x)` by Gauss–Legendre
    rules doubled until they repeat to 64 ulps.
  - `DgDuct::new_with_gas_field` consumes the same `AxialGasField` as the
    frequency-domain propagator. The material owner reevaluates primitive
    temperature/composition at every mass quadrature point, giving the actual
    `A/(rho*c*c)` and `rho/A` weights. Derived gas coefficients are never
    interpolated. Both matched boundaries use their own exact endpoint gas.
    Geometry boundaries and nonuniform-field stations are exact element
    boundaries; the requested total element count must cover those spans.
    A uniform field delegates the original constructor bit-for-bit. The
    complete primitive field and source cards remain available through
    `gas_field()`; `gas()` is the inlet reference state, with explicit
    `inlet_gas()`/`outlet_gas()` accessors. Unequal field/geometry lengths,
    insufficient elements and material/source/phase refusals return before
    publishing an operator. Temperature and composition are prescribed and
    frozen: no mean flow, heat/species evolution or viscothermal loss is added.
  - `DgDuct::with_viscothermal(poles, band_hz, maximum_relative_error)` adds
    both viscous-flow and thermal-pressure memory from `fs_phs::ZkDiffusive`.
    This attachment admits identical circular cylinders with uniform material
    gas only; taper, profile and nonuniform gas-field requests refuse explicitly.
    Density, viscosity, gamma and Prandtl number come from the bound material.
    Constant normalized wall coefficients commute with the exact DG mass maps.
    Every positive arm adds a negative-semidefinite coupling to its retained
    storage coordinate. Poiseuille drag and both complex Bessel wall operators
    are retained without adding the inviscid L/C twice. The declared pole band
    is a separate error axis from spatial and exponential convergence.
    `total_dofs()` includes every memory state; the same 2048-DOF cap applies.
    The augmented physical operator changes the state identity, so pre-attachment
    states cannot silently enter it. Initialization starts wall memory at zero.
  - `DgDuct::with_spatial_viscothermal(DgWallResolution)` extends the SAME
    memory/state/exponential owner to the original cone/profile geometry and
    primitive gas field. Positive Gauss points reevaluate the source gas and
    original radius; derived transport is never interpolated or taken from
    the inlet. Each arm couples through `sqrt(w*C0)*Mp^-1/2*phi` or
    `sqrt(w*L0)*Mq^-1/2*phi`, preserving the full DG mass maps.
    Wall quadrature (order+1 through 64), mesh and pole count are independent
    convergence axes. The reported band errors are sampled constitutive
    errors at those points, not a spatial or broadband accuracy certificate.
    A single exclusive memory enum prevents attaching either wall route twice.
    `wall_stored_energy_j` retains wall energy; `wall_dissipated_power_w`
    independently evaluates nonnegative physical drag/relaxation squares for
    BOTH routes. The cumulative ledger still infers loss rather than separately
    integrating it. Uniform exact-mass and spatial routes are compared against
    each other, exact-Bessel cylinder impedance, Poiseuille drag and the closed
    volume isothermal limit. Native execution of the new spatial tests is pending.
  - `DgDuct::with_molecular_relaxation(RelaxationState)` consumes the shared
    material's evaluated internal heat capacities and relaxation times. It
    couples molecular and circular thermal diffusion in one pressure law,
    `kappa(s)=1/(rho R T) [1-R sum_n w_n s/(s cp(s)+d_n)]`, where
    `cp(s)=cp_inf+sum_j C_j/(1+s tau_j)` and the original positive Bessel
    measure supplies `w_n,d_n`. It REPLACES the old thermal block and retains
    the old shear block; adding independent thermal and molecular shunts would
    count compressibility twice. Frozen pressure mass and positive thermal/
    molecular loss squares advance in the same augmented exponential.
    Original equilibrium gas/source metadata and physical p/q observations
    survive. Empty/zero-capacity reservoirs preserve the old operator bitwise.
    The attachment accepts uniform gas with either uniform cylindrical wall
    memory or the original spatial wall quadrature on cones and radius
    profiles, or a bare variable-area adiabatic duct. Each spatial point keeps
    its original local Bessel measure and normalized pressure projection;
    surviving viscous coordinates are explicitly reindexed before the joint
    thermal/molecular coordinates are appended. No radius average is used.
    A spatial attachment adds `elements * quadrature_points * poles * J`
    states for J active molecular reservoirs. That evaluated-state API refuses
    nonuniform gas; the pointwise-law API below handles the original field.
    Evaluated-only attachments lack a source law for live gas changes. The
    source-backed uniform transfer is described below.
    Compliant and Maxwell material walls can subsequently use the frozen
    pressure map. Independent molecular loss is separate from wall loss in
    harmonic responses and their coupled power accounts. The existing dense
    cap includes every new radial/reservoir coordinate.
    The frozen mass rescaling preserves the original pressure mass, not a
    quadrature-lumped replacement. At finite axial quadrature, the effective
    pressure operator is `s [I + sum_q (B_inf*kappa_q(s)-1) v_q v_q^T]`.
    The quadrature projectors sum to identity only when that rule resolves
    the original area-weighted mass. Consequently the exact isothermal or
    equilibrium storage limit is not claimed at unconverged quadrature on a
    tapered or curved profile. Axial quadrature, mesh and radial pole counts
    remain independent refinement axes. This is the local-cylinder Webster
    wall model; spherical-cone and three-dimensional corrections are absent.
    `with_molecular_relaxation_field[_controlled](RelaxationAbsorption)`
    resolves that SAME material law through the original `AxialGasField` at
    every original thermal quadrature point and adaptive pressure-mass point.
    A varying field requires the spatial wall owner. The frozen pressure mass
    is integrated on the original element basis, retaining all source/geometry
    station splits. For `T=Mp_inf^(-1/2) Mp_eq^(1/2)`, retained pressure
    couplings use the congruence `T L T^T`, and ports use `T B`; the roots are
    not assumed to commute. Every old thermal functional becomes
    `sqrt(C_inf(x)/C_eq(x)) T v_eq`, with its own local cp and reservoir rates.
    Shear histories, original matched boundary impedances and physical p/q
    remain intact. Subsequent compliant walls consume the frozen mass map.
    Uniform fields delegate the original path bitwise; an all-zero law keeps
    the original operator. Non-relaxing regions retain their thermal arms.
    `molecular_relaxation_field_dofs[_controlled]` counts the exact complete
    state before dense allocation; inactive zero-capacity reservoirs do not
    consume states. Count admission validates field stations and thermal
    points, while construction also checks adaptive mass points. At finite
    quadrature the physical pressure map is the exact frozen mass plus the
    quadrature of `A(x)*(kappa(s,x)-1/B_inf(x))`; its isothermal limit still
    needs axial refinement. No temperature/composition average, background
    transport or live molecular-history transfer is substituted. Five native
    regressions cover local constitutive Schur elimination, quadrature
    refinement, original port/wall/loss work, uniform/zero parity and local
    reservoir count/domain refusals; their native execution remains pending.
    The original equilibrium Bessel band/error is retained; it is NOT a new
    accuracy certificate for the combined constitutive approximation. No bulk
    SK term, axial heat conduction, evolving chemistry or finite wall thermal
    feedback is implied. Four focused regressions cover shared adiabatic
    compliance, combined isothermal storage, physical loss, time/spatial
    refinement and rate-data pressure changes. Spatial regressions additionally
    compare the actual augmented Schur complement with the independently
    integrated local heat equation, refine a curved-profile quadrature, compare
    the uniform and spatial coupled pressure/loss observations, and check
    remapped shear loss, side ports and exact-time halving. Native execution
    of these new molecular tests is pending.
  - States carry a reproducible identity of the assembled operator, geometry,
    mass-coordinate transforms, boundaries and typed ports. `step`, fallible
    `sync_nodal`, and fallible `l2_pressure_error` refuse incompatible states
    before mutation, including same-sized states from a different medium.
    Independently rebuilt identical operators remain compatible. Changing
    `dt` on the same operator reanchors the lossless modal cache at the current
    physical state; it does not reinterpret elapsed steps with the new `dt`.
  - Interior faces use the alternating flux (`q̂` from the left, `p̂` from the
    right). It is energy-conserving and optimal (`p + 1`) for every `p`;
    `K_q = -K_pᵀ` holds exactly through LGL summation by parts.
  - Boundaries: `Closed`, `IdealOpen`, `Matched` (upwind characteristic)
    and `Port(PortId)`. A port holds an inflow `q` over a step and returns
    the exact impulse `ΔX = ∫ p dt`, so its work is exactly `e ΔX`.
  - `propagator(dt, tol)` is the exact exponential, with no RK, trapezoid,
    implicit-midpoint or generalized-alpha step.
    - Without wall memory and with lossless boundaries the energy-coordinate operator is skew. It is
      diagonalized by the SVD of its coupling block, from the symmetric
      Jordan–Wielandt matrix (`householder_eigh`). The certificate is the
      larger of the relative singular residual and the basis orthogonality
      defects.
    - `Matched` and viscothermal memory are dissipative and non-normal. They use a scaled Taylor
      series with squaring, certified by a rigorous truncation bound carried
      through the squarings plus an a-posteriori halving estimate.
    - A certificate above `tol` refuses (`FS-DUCT-PROPAGATOR-CERTIFICATE`).
    - There is no Krylov path: more than `MAX_DG_PROPAGATOR_DOFS = 2048` DOFs
      refuse in `new` and `propagator`.
  - `DgPropagator::step` returns each port's `ΔX` and an `EnergyLedgerRow`
    whose stored energy is that of the reconstructed physical state.
    - For the lossless image and lossless boundaries `E(t) - E(0) - W` is a conservation check.
    - For `Matched` or wall memory the absorbed energy is inferred from the deficit
      (unclamped), so `balance_checked` is false.
  - Origin: QuietLotus's uncommitted 2026-10-16 work, reworked (see the
    commit for what changed).

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

  `TabulatedMouth` optionally binds that table to an explicit axial mouth
  reference, radius, resolved exterior gas and finite nonnegative real DC
  impedance. `validate_mouth` requires geometry agreement within 1 nm, exact
  gas agreement and every stored impedance row to be passive at a fixed
  1e-12 relative tolerance, regardless of the source import tolerance.
  This is caller-declared metadata, not validation of the source measurement
  or bake. The root `MouthLoad` re-export is the same selection type used by
  the modal authority. Existing frequency-domain calls remain compatible;
  scalar time-domain consumers require the physical reference and full band.

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

`src/dg/tests.rs` (IRP M2.2-W3a), one JSON line per rung:
- `dg_lgl_nodes_exact` (G0): exact to degree `2p - 1` (1e-14), and
  summation by parts.
- `dg_mass_matrices_integrate_the_exact_area` (G0): a unit field mapped
  through `to_energy_coords` gives `1ᵀ M 1` per cone element, which matches
  the closed-form `∫ A/(ρc²)` and `∫ ρ/A` to 1e-13; measured 2.7e-15.
- `dg_cylinder_ir_matches_dalembert` and `dg_cone_ir_matches_spherical_exact`
  (G1/G2): Gaussian pulses at rest travel 0.75 L and reflect once.
  - Exact solutions use images: even at a closed end and odd at an open end
    for the cylinder; `r p` with odd images for the open/open cone.
  - The h-ladder, K = 24..96, is gated at `p + 1 ± 0.2` through fs-mms.
    Measured orders: p = 2 gives 3.01 and 3.06, p = 3 gives 4.02 and 4.03,
    p = 4 gives 5.00 and 5.02 (cylinder, then cone). K = 16 was
    pre-asymptotic.
  - The p-ladder at K = 16 at least halves the error per order.
  - A linear `Profile` reproduces the `Cone` to 2.4e-14.
- `dg_energy_ledger_lossless` (G0):
  - Over 1e5 steps on closed/closed and closed/open, physical-energy drift
    is 1.3e-15 and 1.6e-15 (1e-12 allowed). The ledger agrees with the
    independently computed nodal energy to 6e-16.
  - Two ports balance their exact `e ΔX` work to 2.7e-13.
  - A matched end is passive and absorbs the pulse.
- `dg_impedance_matches_tmm_lossless` (G2): the first 6 peaks of a port/open
  cylinder against the Lossless TMM poles, within 0.1 c.
  - The DG impedance's own poles (port-coupled `σ_k / 2π`) match to 7.2e-8 c.
  - FFT of a 1 s windowed impulse response matches to 0.0065 c. The window
    biases peaks as `σ²`: 0.25 s gives 0.099 c.
- `dg_port_coupling_dt_ladder` (G1): a step-averaged held effort converges
  at order 2 ± 0.2 in Δt; measured 2.05.
- `dg_propagator_certificate` (G0): measured 1.5e-15 (skew) and 1.3e-14
  (matched) against 1e-13; an unachievable tolerance refuses by name.
- `dg_refusals` (G0): order and element caps at cap and cap + 1, the DOF
  cap, `A <= 0`, tone holes, an empty duct, shared port ids, `dt` and
  `tol`, and unknown, non-finite or foreign port inputs.
- `dg_bitwise_repeat` (G5).

`src/modal/riccati/tests.rs` (IRP M2.1-W3b), one JSON line per rung:
- `mm_bessel_certified_to_gamma16` (G0): `bessel_j0_certified`, `J1` and
  the `J0' = -J1` identity match an independent spectral quadrature to
  1e-13 on [0, 64]; measured about 1.7e-15. It reports, without gating,
  the legacy `bessel_j0`'s 1.4e-7 at 26.
- `mm_mode_cap_refusal` (G0): the computed cap is 21; cap and cap+1;
  steps 0 and cap+1; the evanescent-step guard at its boundary.
- `profile_slope_matches_central_differences` (G0).
- `mm_riccati_n1_equals_plane_riccati` (G2): N = 1 matches the W3a
  plane Riccati to 1e-10.
- `mm_riccati_cylinder_decouples` (G2): closed-form `i Zc cot kL` and
  matched `Zc_n` to 1e-12; off-diagonals exactly zero.
- `mm_riccati_matches_extrapolated_staircase` (G1/G3): the staircase is
  O(dx), gated at order 1 +- 0.2 on slice multipliers 4/8/16/32. Its
  first-order Richardson limit agrees with the continuous result to
  1e-4 on a cone and a 400-cone exponential polygon; measured 2.5e-5
  and 1.3e-5. Multipliers 1/2/4/8 were pre-asymptotic (order 1.38).
- `mm_riccati_step_ladder` (G1): order 4 +- 0.2 via fs-mms; observed
  4.00.
- `mm_riccati_mode_ladder` (G1): uses the Amir, Pagneux & Kergomard 1997
  Table I trombone bell. `|dZ00|` shrinks monotonically for N = 1 to 16.
  Liveness: N = 1 vs N = 8 moves its peaks by -15.5, -36.6 and -61.0 c.
- `mm_riccati_pagneux_cases` (G2): the exact one-mode exponential horn
  (Part I App. E) to 1e-9, above and below cutoff; measured 9e-14.
  Also the Part I 45-degree cone at 7 and 17 modes: reciprocal, with a
  shrinking tail.
  - It also checks the one open figure by the same authors: Amir, Pagneux
    and Kergomard 1995 (WIT Computational Acoustics, Fig. 1). That is `|p|`
    on the axis of a 0.5 to 4 m exponential horn at 0.3 and 0.9 f_c, read
    at 20 points under OD17 with half-widths 0.005 and 0.015.
  - The reads follow the paper's ANSYS line, a converged reference. So the
    comparison uses our converged value: |p| at the admitted mode cap (21).
    The allowance is the half-width plus the declared algebraic mode tail
    (`mode_tail`, shared with `mm_riccati_mode_ladder`). That tail is fitted
    per point on the N = 1..16 increments, carried to the cap and multiplied
    by the safety factor 2.
  - The agreement is reported, not gated. Measured: 18 of 20 reads agree.
    - z = 4 m at 0.9 f_c misses by 0.0003 (0.0168 against 0.0150 + 0.0015).
    - z = 4.5 m at 0.9 f_c has no declared tail: its fitted order is 0.26,
      against the model's own order > 0.5 condition. The on-axis series
      alternates between even and odd N near the mouth.
  - Gated: N = 1 misses both piston reads by more than the half-width (it
    misses 19 of 20 reads), and step refinement moves every compared value
    by less than 1e-3 of the half-width (at most 3.3e-9).
- `mm_field_cylinder_modes_closed_form` (G2): per-mode outgoing and
  decaying waves on a cylinder, with axis weights checked against quadrature,
  to 1e-12; measured 9.4e-14.
- `mm_field_n1_equals_exact_webster` (G2): the one-mode field of the
  exponential horn station by station to 1e-9; measured 8.5e-13.
- `mm_field_power_flux_conserved` (G0): the lossless axial power flux is
  constant over every station to 1e-10, with propagating (8 kHz) and
  evanescent (600 Hz) higher modes; measured 4.9e-14. Also checked: the inlet
  `P = Z U_in` and the matched mouth `P_n = Zc_n U_n`.
- `mm_field_refusals` (G0).
- `mm_riccati_reciprocity_and_power` (G0): `Z = Z^T` and lossless power
  balance to 1e-10.
- `mm_riccati_bitwise_repeat` (G5).

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

- Time-domain DG is frozen, no-mean-flow plane-wave propagation. Spatial wall
  memory uses the local cylindrical ZK/Webster approximation on a taper;
  it is not the spherical-cone correction or a full 3-D wall solution.
  - Finite wall-temperature feedback, independent integrated wall dissipation
    and the full impedance-peak gates of M2.2-W3b remain open.
    There are no side branches or radiation (W3c), moving geometry (W3d) or
    valve dynamics (W3e). Tone-hole segments refuse.
  - The matched boundary's certificate includes an a-posteriori estimate,
    not an enclosure, and its absorbed energy is inferred, not checked.
  - Propagators are dense, `O(N³)` to build, and capped at 2048 DOFs; there
    is no Krylov path.
  - Port coupling holds the effort for a step, which is first order in Δt
    for point-sampled inputs and second order for step averages.

- Continuous multimodal image (M2.1-W3b): published Pagneux benchmark
  results exist only as figures, so no 1 % comparison against printed
  values is claimed. The one open figure (WIT 1995 Fig. 1) agrees at 18
  of 20 reads, at read error (1 to 5 %) plus the declared mode tail. That is
  reported, not gated. Otherwise agreement is shown against exact
  solutions, the independent staircase limit and conservation laws.
  - Mode convergence is algebraic. On the 45-degree cone, 7 vs 17 modes
    differ by 3.1 % in `Z00`, so `tail_bound_estimate` is an indicator,
    not a bound.
  - Higher modes terminate matched (no radiation matrix), per-mode
    boundary-layer loss is approximated as in the staircase, and holes
    and C0 profiles are out of scope.
  - Accuracy above `MAX_EVANESCENT_STEP` is refused, not degraded.

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
  No roughness or distributed porous liners; the network's separate
  `PorousRigid` termination is a normal-incidence boundary. No mean flow (QUIESCENT medium,
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
  - Refused above `FLANGED_MAX_KA = 24`. This is the fixed-node
    quadrature's design limit, not a validity edge of the model (the
    multimodal load is exact for every `ka`). It was 20 until
    `flanged_design_is_measured_across_and_beyond_the_cap` measured the
    design: the quadrature `Z_00` against the closed-form piston is 3.5e-14
    at `ka = 20`, at most 6.5e-14 through 26, 2.9e-13 at 28, and loses the
    1e-12 gate at 30 (2.0e-11; 3.6e-8 at 34, 1.5e-5 at 38). The 16-mode
    against 24-mode load differs by 1.1e-3 to 1.4e-3 over `ka` 20 to 28, flat
    in `ka`, so mode truncation does not worsen with the lift. A mode exactly
    at cutoff drops out of the Schur complement (`ζ = ∞`, the continuous
    limit), so the new cutoff inside the band (`α_7 = 22.76`) is handled the
    way `α_1` to `α_6` were. Values at `ka <= 20` are bit-identical to before.
  - The plane-mode TMM above a bore's first cross-mode cutoff stays the
    disclosed 1-D boundary.
  - A full exterior BEM of an unflanged pipe is still the recorded successor.
- Tone holes are the compact T-junction (Nederveen series + Dalmont
  inner + chimney wall law + `Y = σ Y_open + (1−σ) Y_closed` vent
  mix + mutual series `t_m ∝ e^{−s/a}` from neighboring open
  holes). That mutual series term is a bore-junction correction, not an
  exterior radiation matrix. The Ernoult dataset is the regression floor.
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
through its stable local-wave scattering composition and physical admittance
stamp. If an edge's pressure-to-flow chart or the resulting nodal solve is
singular or exceeds its condition allowance, the same composed edge scattering
operators are retained in the complete graph solve. The equations
`(I-S)p - (I+S) R U = 0` and nodal volume-flow balance preserve finite loaded
half-wave resonances and antiresonances without an artificial loss, frequency
shift or geometry subdivision. No opposite signed load probes or global
growing ABCD product are needed. Propagation, loss and local chimney physics
have one owner; this alternate algebraic elimination does not change them.

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
zero axial-length edges and absent ports refuse. The retained-wave route adds
two inward-flow unknowns per edge, capped at 768 total unknowns before allocation.
It uses pressure-valued currents `R_edge U` and row equilibration; the ordinary
nodal route retains its original arithmetic and smaller factorization. Both
complete solves have full infinity-norm condition estimates and refuse above
1e10 or on nonfinite/singular operators. In particular a genuine undamped
closed-network pole still refuses. `wave_constraint_edges` records the selected
representation (zero for the nodal route). Accepted edge estimates report
composition and admittance conditioning on the nodal route, and composition
conditioning alone when currents are retained. The graph condition then refers
to the dimensionless row-equilibrated constraints. These are floating-point
solves, not interval accuracy certificates for every graph. The returned
open-circuit N-port impedance itself must exist: a load attached later through
`drive` cannot regularize an already singular open-circuit network query.

G1/G3 tests `network_retains_scattering_at_half_wave_and_low_frequency`,
`network_half_wave_subdivision_preserves_loaded_response`,
`network_half_wave_antiresonance_is_not_a_closed_cavity_pole` and
`network_half_wave_multiport_retains_reciprocal_mutual_loading` compare analytic
matched loads, complex terminated cylinders, zero-pressure antiresonances and
reciprocal rank-one N-port impedances. They retain the exact queried frequency;
the subdivision comparison includes detuning on both sides of the edge poles.

Loads are explicit closed, pressure-release, passive finite impedances,
tabulated loads, a material-backed rigid porous termination, or an explicit radiation model/radius. Exact thin-wall LS is
consumed through its existing owner and explicit quadrature; legacy fits keep
their own refusals. `NetworkResponse::drive` / `driving_point` require exactly
one assignment for EVERY other external port, and the same `GasState` as the
network solve. They return pressures and currents as well as input impedance;
there is no omitted-port load, medium, wall or clock fallback. The uniform API
keeps its single explicit gas; the regional API below binds each edge separately.
Continuous axial fields within an edge and time-domain characteristic networks
are not claimed here. Calls are synchronous and bounded; frequency sweeps can be
caller-chunked, with the same existing no-Cx boundary.

### Directly loaded driving-point response

`network_drive`, `network_drive_coupled` and their `_regional` variants impose
unit flow at one active port and an explicit scalar or coupled load at every
other port. They preserve the existing open-circuit query followed by `drive`
operation order whenever that solve succeeds. Only `Singular` selects the
complete loaded scattering system; invalid inputs retain their original refusals.

The loaded KCL rows retain external currents, while boundary rows impose
`p + Z U = 0`, with `U` into the graph. Closed ports impose `U=0`. The load matrix
is never inverted, so pressure-release boundaries and singular passive mutual
loads remain valid. At most 31 external currents add to the retained edge
currents, capped at 799 total unknowns before the dense allocation. The same
row equilibration and full infinity-norm condition ceiling of 1e10 apply; a
true pole of the fully loaded network still refuses. Labels, active-port
partition, frequency, passivity, edge gas records and exterior admission remain
mandatory. This solves the supplied load problem without requiring its
open-circuit N-port impedance to exist.

All four `fs-couple::NetworkDrivingPoint` response methods use these entry
points, including their causal line, tract and radiation consumers. The original
N-port query continues to require finite open-circuit `Z`.
`network/tests/loaded.rs` defines regressions for half-wave scalar and mutual
loads, pressure release, a rank-one passive load with a pressure-release
combination, mean power balance, regional gas binding, bit parity with the
original admitted route, and invalid load refusals. The driving-point adapter
has a matching half-wave regression.

### Regional material and exterior binding

`NetworkMedia` binds a complete original `GasRecord` to every exact edge id,
plus one independently specified homogeneous exterior. Mixture composition,
NASA-9 query states and source-pack receipts remain in those records and in
the resulting `NetworkResponse::media()`. Immutable shared ownership avoids
copying all source records for each frequency. Duplicate, omitted and unknown
edge ids, incoherent media or unequal static pressures refuse before assembly.
No temperature, composition or transport coefficient is inferred or averaged.

`network_port_impedance_regional` uses the same nodal assembly and dense solve
as `network_port_impedance`. Each edge's own rho, c, gamma, viscosity and thermal
transport feed the existing wave/viscothermal owner. Embedded chimneys retain
their edge's interior gas while their exposed mouths receive the declared
exterior. All graph radiation and porous loads use that same exterior. A flanged
mouth also retains its actual incident interior gas in the original higher-mode
Schur load. Multiple unlike incident gases at one mouth do not define one modal
port plane and refuse that load; an explicit finite/coupled boundary remains
available. Exact thin-wall Levine–Schwinger radiation retains its single-fluid
domain and refuses an interior/exterior contrast rather than reinterpreting it.
`drive_bound` and `drive_coupled_bound` retain it for deferred labelled-port
loads. Pressure and volume flow are continuous at the original graph junctions;
there is no added reflection law or material-dependent output gain.

This is a stationary, pressure-matched gas field. Mixing, diffusion, mean flow,
moving interfaces and heat/species evolution are absent. Each edge is spatially
uniform; callers may explicitly subdivide geometry into edges. Cross-edge
exterior interaction still requires the existing full `CoupledTermination`
matrix from its physical producer; a regional binding does not supply one.
Existing chain-local compact-hole and radiation validity limits remain in force.

`tests/network_media.rs` defines seven focused regressions: uniform bit parity
and original-record retention; independent serial/parallel branch solutions
with N2/O2, He and CO2; separate interior/exterior and deferred load equivalence;
embedded chimney equivalence; reciprocity, incoming power and lossless pressure
scaling; ambiguous or unsupported modal mouth media; and invalid
identity/pressure/coherence refusal. These definitions and
source checks are not native numerical or experimental validation.

`NetworkLoad::PorousRigid` uses the existing `PorousGasLaw` and JCA(L)
material owner in both graph-boundary assembly and labelled-port termination.
The exact solve gas refreshes drag, inertia and thermal relaxation while pore
geometry and intrinsic permeability remain fixed. Its normal-incidence surface
impedance is conjugated from `exp(+i omega t)` to the network convention and
divided by exposed area to obtain `p/U`. Invalid geometry, incoherent gas,
source-band departure and unrepresentable conversions refuse. The existing
response gas check also prevents using an old network with a new atmosphere.
Independent cylinder solutions and `0.5 Re(p conj(U))` power balances cover
both routes, lossless and Bessel duct walls, and data-only gas/pore changes.
This uniform-pressure termination does not establish distributed lining,
moving porous frames, liquid saturation or transient thermal deposition.

`NetworkResponse::drive_coupled` also takes one `CoupledTermination`: a
row-major acoustic impedance matrix [Pa s/m³] over named external ports, with
`p = -Z U` for currents into the network. It carries the response's exact
`omega` (bitwise) and a nonempty `source_id`, for example an fs-bem aperture
matrix declaration. Its ports are always unknowns and are disjoint from the
SISO `loads`; together they terminate every non-driven port exactly once. The
block enters the port Schur complement whole, `A = Z_live + Z_L`. Admission
refuses non-finite or non-square matrices, a frequency mismatch, and a
Hermitian part with an eigenvalue below `-1e-12 max|Z|`. Reciprocity is not
required: a collocation radiation matrix is reciprocal only to its
discretization error. `drive` is the no-coupled case with unchanged arithmetic.

Compact chains retain independent exterior self loads for embedded tone-hole
mouths. The automatic six-iteration scalar correction has been removed: it
assumed a flanged downstream bell before the actual terminal was supplied and
clipped signed mutual work/reactance. It did not solve the declared coupled
boundary problem. Existing compact bore-junction corrections remain intact;
axial hole spacing alone does not declare three-dimensional mouth geometry.

Exterior interaction requires explicit mouth ports and a complete
`CoupledTermination` from the physical radiation owner, such as the shared
fs-bem aperture matrix consumed by fs-couple. Passivity constrains the whole
Hermitian matrix; a particular opening's mutual work may be negative. The
standalone `tests/compact_boundaries.rs` regressions check local linear stamps
under different actual terminals and an analytic passive matrix with negative
mutual resistance and opposite-signed even-mode reactance. They are numerical
test definitions, not experimental validation of an instrument or aperture.

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
`network_coupled_termination_even_mode_and_diagonal_limit` drives a symmetric
Y-branch. A coupled `[[z,m],[m,z]]` load equals the SISO even-mode load `z+m`
within 1e-12 and differs from the uncoupled answer by more than 1e-3. The
diagonal-matrix limit and a mixed one-port coupled case are bitwise equal to
SISO `drive`, and each admission rule refuses.
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

Network assembly consumes the finite admittance `(I-S)(I+S)^-1/R` from this
same owner, including long opaque edges, when that pressure-to-flow chart is
admissible. It retains the original scattering constraints at component poles
as described above. `network_chain_two_port` remains an
explicit finite, conditioned ABCD export; graph assembly no longer depends on
its envelope or reconstructs ABCD through opposite signed load probes. Existing
network complete-solve condition caps remain 1e10. Local compact corrections
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


## Continuous DG response and retained-state composition

`DgDuct::continuous_operator()` borrows the original complete energy-coordinate
generator and original typed inflow columns: `x' = L*x + B*U`, pressure `B^T*x`,
and storage `x^T*x/2`. Gas, geometry, endpoint laws, side ports, spatial mass
maps, every viscothermal arm and compliant/Maxwell wall coordinate are retained.
The view carries the original operator identity. It permits joint assembly at
an isolated duct impedance pole without dividing by that impedance.

For lossless boundaries, the stored generator is exactly `[0,-G;G^T,0]`, the
operator already advanced by `SkewModal`. The independently mass-normalized
`H` must first agree with `-G^T` within `64*n_local^2*epsilon` of their common
maximum norm (`n_local <= 9`). Larger or nonfinite transpose defects refuse;
floating-point normalization roundoff is not represented as physical damping.

`DgDuct::harmonic_response` solves `(-i*omega*I-L)*x = B*U` with the existing
cancellable complex LU. Phasors are peak amplitudes in `exp(-i*omega*t)`.
`DgContinuous::observe_harmonic` observes a complete supplied original state,
including one solved jointly with another physical owner, and checks those
same equations using the actual supplied port flows. Undeclared or duplicate
input ports refuse; omitted declared ports receive zero flow. The admitted
viscothermal frequency band is enforced. Accepted time histories are untouched.

Responses include physical nodal pressure/axial flow, conjugate port pressure,
mean gas/wall storage, original positive wall-loss squares and matched-boundary
absorption. The latter is `(|p|^2/Zc + Zc*|q|^2)/4` at each matched endpoint,
using its actual local gas and area. Mean port work is `Re(p*conj(U))/2`.
The independently computed power residual is diagnostic. The infinity-norm
backward residual is gated by the caller's explicit tolerance; it is not a
conditioning or forward-error certificate near an undamped resonance.

Before allocation, the standalone solve charges `n^2*(n+ports+8)` work and
observation charges `n^2*(ports+8)`, with the existing complete-state DG cap.
Cancellation, singular LU, invalid budgets/tolerance, nonfinite solutions or
observations, and failed original-equation residual checks return no response.
This is the frozen linear DG spatial model, with its original plane-wave,
constitutive-band and geometric validity limits. It is not the finite-step
held-flow transfer or a replacement for spatial, wall-pole or coupling-clock
convergence. External loads enter only through an explicitly composed owner.

## DG interior point-flow junctions

`DgDuct::with_side_ports` attaches at most 32 explicitly identified interior
mesh faces. The downstream pressure trace is the original alternating-flux
trace. Its held inflow changes the local through-flow jump by that inflow and
uses the same pressure integral for exact port work. Only conjugate port
columns and their operator identity change: the acoustic, constitutive wall,
energy-coordinate and exponential owners remain unchanged. Existing states
cannot be reinterpreted after changing this port basis. An empty attachment
preserves the original owner and identity.

`side_port_source` returns the actual face station, local bore area and complete
original gas-field evaluation; a bore area is never presented as a hole area.
Finite chimneys, aperture geometry, junction corrections, moving pads and
exterior radiation require explicit coupled owners. All unknown loads must
consume the complete `multiport_response` together to retain their same-step
interaction; independent scalar boundary solves are not equivalent. This is a
zero-volume plane-wave junction, not a resolved three-dimensional tone hole.

Focused regressions cover pressure/work conjugacy, integrated discrete
continuity, source identity and local material fields, wall attachment order,
and equal-pressure/opposite-work composition with a second original DG duct.
Native execution and instrument qualification remain pending.

## Reciprocal compliant walls in the DG time domain

`DgDuct::with_compliant_wall` appends the original `WallPin` radial displacement
and velocity to the same energy-normalized DG operator and dense exponential.
For circumference P, area A and gas bulk modulus B, continuity receives -P v
while the wall receives pressure through sigma v_t = p - K w - r v.
The reciprocal terms cancel in the complete energy rate; wall storage is
integral P(sigma v^2 + K w^2)/2 and its independent loss is integral P r v^2.
`compliant_wall_state` reads physical motion and storage; the existing
`wall_dissipated_power_w` includes both material and gas-wall loss once.

This cold attachment requires uniform gas and identical circular segments.
Viscothermal memory must be attached first and is retained unchanged. New
states start with the wall at rest; old states refuse by operator identity.
All boundary/side ports see wall reaction through their original pressure
columns. A matched end remains the gas characteristic resistance, not an
assertion of reflectionless termination for an elastic tube. Source bandwidth,
small-displacement applicability and refinement remain caller responsibilities.
Axial shell waves, thermal evolution and external wall radiation are not
supplied. The dense ledger still labels dissipation inferred from its
deficit; the independent instantaneous loss is not an integrated heat proof.

G0/G1 regressions compare actual pressure, wall displacement and velocity to
the independent sealed-tube oscillator, check conservative port work and total
storage, and reconstruct the combined gas/material dissipation from x^T L x.

`with_compliant_wall_relaxation` consumes the original `GeneralizedMaxwell`
spectrum using the same equilibrium thin-wall reduction, K_j=(K_inf/E_inf)E_j.
Each positive arm adds radial reaction pressure z_j with z_j'=K_j v-z_j/tau_j,
work-conjugate reaction -z_j, storage integral P z_j^2/(2K_j), and loss integral
P z_j^2/(K_j tau_j). Every arm advances together with gas and radial motion in
the existing exact-ZOH exponential; no midpoint or separately lagged material
step is introduced. Exact zero arms add no state and preserve the original
operator identity. Positive coefficients that cannot be represented and total
states beyond the original dense cap refuse. New memories start relaxed.

`compliant_wall_state` reports each arm's physical radial pressure and its
spring energy; complete wall storage and loss include all arms exactly once.
Original source provenance and frequency applicability belong to the material
adapter. Frozen geometry/state and uniform-gas restrictions remain in force.
An independent decaying sealed-gas eigenmotion exercises pressure, wall motion,
branch pressure, complete energy and physical loss. Further regressions cover
unchanged original gas/port blocks, combined passivity, zero-arm equivalence,
foreign states and invalid spectra. These new tests await native execution.

## Nonuniform material gas coupled to the same compliant DG wall

The original prescribed `AxialGasField` now composes with elastic, Kelvin and
source-Maxwell walls on identical circular segments. The existing pointwise
material evaluations still own each acoustic pressure mass
`Mp = integral A/(rho c^2) phi phi^T dx`. Wall inertia and stiffness instead use
the exact geometric polynomial Gram matrix `G = integral phi phi^T dx`.
For circumference P and surface density sigma, the shared energy operator uses
`C = sqrt(P/sigma) Mp^(-1/2) G^(1/2)`: pressure receives `-C` and wall velocity
receives `+C^T`. These matrices need not commute. No inlet, nodal or averaged
bulk modulus is substituted for the field. Physical wall velocity is recovered
by `G^(-1/2)/sqrt(P sigma)`; displacement and every original Maxwell reaction
use the corresponding spring-energy maps. Gas and material memories remain
separate, and all evolve in the original exact DG exponential.

The uniform gas and uniform-field routes retain their original operator and
identity. Existing gas diffusion, boundary/side pressure ports, Maxwell terms,
energy coordinates and dense state caps remain authoritative. The controlled
wall constructors support bounded cold cancellation; no interrupted candidate
is returned. Frozen gas/material state, uniform cylindrical geometry, local
thin-wall applicability and caller-resolved source bandwidth still apply.
This does not supply evolving gas composition, tapered shell dynamics or a
separately integrated heat receipt; the existing dense-ledger limitations remain.

Three exact-production coordinate tests passed natively: independent two-region
sealed-pocket oscillators, noncommuting mass roots with physical force/storage,
and uniform limits plus invalid/cancelled preparation. Five full DG/source-pack
regressions cover pointwise weak continuity, Maxwell force/storage/loss, retained
gas memory and port work, original uniform identity and caps, actual composition
and solid-source changes in propagated pressure/wall motion, and source refusals.
Those full crate tests await native execution; leaf evidence does not certify
the complete material-query/DG exponential integration.

`with_compliant_wall_regions` extends the same operator to one independent
elastic/Kelvin--Voigt `WallPin` per original element. Material boundaries align
with element boundaries; each region retains its own surface mass, stiffness
and resistance. The existing gas/wall coupling is rescaled by the local wall
mass, with opposite transpose entries retained exactly. Physical displacement,
velocity and loss use those same local coefficients. Gas mass matrices,
viscothermal memories, boundary/side ports and the time integrator are retained.
Identical regional coefficients take the original homogeneous path bitwise.
`compliant_wall_at_element` reports actual base assignments; `compliant_wall`
returns no homogeneous pin for differing base coefficients. Equal base pins do
not imply equal hereditary spectra. Preparation is cancellable and old states
refuse by operator identity. G1 tests reconstruct each local physical
equation and independently integrate energy/loss for uniform and spatial gas;
G3 tests exercise original pressure-port evolution and data-only region swaps.
`with_compliant_wall_relaxation_regions` retains each element's independently
declared generalized-Maxwell spectrum beside its original base pin. Every
positive arm allocates only that element's nodal state, with its own stiffness
and relaxation rate. The same reciprocal gas/wall operator supplies its force,
storage and loss; no common spectrum or additional constitutive stepper appears.
Homogeneous assignments retain the original full-field ordering and operator
identity. Regional pressure rows follow element then arm order and are exactly
zero outside their support. Native source, pressure, force, storage and port
regressions cover differing spectra and mixed hereditary/Kelvin regions.
The local membrane model still omits axial shell waves and stress transfer at
material joints. These attachment APIs initialize cold mechanics; no internal
history is silently discarded. Live elastic/Kelvin--Voigt and explicitly mapped
Maxwell transfers use the separate continuing-state operations below.

## Live material updates on the native gas/wall state

`DgDuct::rebind_compliant_wall_regions` changes fixed-geometry wall coefficients
on an accepted trajectory. It preserves gas wave and viscothermal coordinates
exactly, radial displacement, and conjugate wall momentum. The original native
reciprocal operator is rebuilt with each element's new coefficients. The actual
spring/kinetic storage change enters `cumulative_parameter_work_j` once, apart
from port work, initial energy and inferred dissipation. The parameter-changing
physical owner must settle that signed work; this operation does not infer heat
or a mass-transfer process. A finite explicit work budget bounds each update.

Original state/operator identity, passive coefficients and cancellation are
checked before publication; old propagators refuse after an operator change.
The elastic/Kelvin--Voigt API refuses retained Maxwell histories; callers use
`rebind_compliant_wall_relaxation_regions` with a complete per-element
`DgWallRelaxationTransfer` mapping for hereditary updates. The source owner must
establish each arm's physical identity. The native transfer preserves every
elastic extension `z_j/K_j` and hence viscous offset `w-z_j/K_j`, including when
uniform and regional storage order differ. It requires one-to-one coverage of
every original and replacement positive history on its original element.
Missing/duplicate histories, activated/extinguished arms and unsupported maps
refuse. No projection or reset substitutes for correspondence.

Each normalized memory coordinate scales by `sqrt(K_new/K_old)`. Its actual
spring-storage change joins the existing parameter-work account exactly once.
A relaxation-time-only change owes zero instantaneous work and changes future
native dissipation. The same regional assembler and propagator continue the
trajectory. Native tests check analytical post-change wall
motion/pressure, physical quadrature of parameter work, retained gas-memory
coordinates, one-time energy accounting and transactional refusal. This does not
remap gas, mesh or expanding geometry, or infer physical arm identity. Five
focused native/source regressions add analytic Maxwell continuation, physical
storage quadrature, source-temperature response and nonzero-history rollback.
Rust parsing and source review passed; these new native regressions have not
been executed in the implementation environment.

## Live uniform gas changes on a continuing DG trajectory

`DgDuct::rebind_uniform_gas` and its source-field counterpart retain acoustic
compression `p/B` and axial momentum `rho*q/A` through a prescribed,
zero-duration uniform material intervention at fixed geometry. The original
energy coordinates scale by `sqrt(B_new/B_old)` and
`sqrt(rho_old/rho_new)`. Native geometric mass maps, the wave operator,
matched boundaries and pressure/volume-flow ports receive the same gas
coefficients; physical pressure and flow change consistently.

Wall displacement, velocity, elastic/Kelvin response, regional Maxwell
histories, clock, initial energy and prior dissipation remain intact. The
actual acoustic storage change enters signed parameter work once. Physical
wall observation maps are updated consistently with the gas normalization.
This work belongs to the owner changing the background material; it is not
heat, loss or an inferred mass/species injection.

Source-backed updates retain the complete original replacement field catalog,
queries, receipts and phase policy, including catalog species whose local
fraction is zero. The complete field must cover the original length and be
spatially uniform. Unsupported nonuniform gas and native viscothermal gas
closures refuse; their histories are neither averaged nor reset. Regional
Maxwell wall histories remain supported.

All owner/state, work-budget, arithmetic and cancellation checks precede
publication. A new native propagator is required; the old operator identity
cannot continue the changed state. This interface supplies no background
thermal/species dynamics, mean flow, static pressure traction or nonuniform
temporal-interface model.

Three native regressions cover an independent one-element wave continuation,
signed work and subsequent ledger closure, regional wall memory and physical
acceleration, port/matched-boundary equivalence to fresh assembly, and atomic
refusals. Source-wall integration tests live in fs-couple. Parsing and source
review do not establish executed native or empirical validation.


### Live uniform gas changes with native viscothermal history

`DgDuct::rebind_uniform_viscothermal_gas` and
`rebind_uniform_viscothermal_gas_field` extend live material intervention to
the original uniform circular `ZkDiffusive` owner. Fixed radius, mesh, exact
angular band, pole count and dimensionless Jacobi sections identify the
physical arms. This is not arbitrary state-space matching. The replacement
gas independently passes the original band's complex-response and real-loss
admission under the supplied error budget. Spatial/nonuniform memories,
changed support and creation or removal of thermal arms remain refused.

Each viscous arm preserves momentum and each thermal arm preserves compression
volume. Their physical inertance/compliance are respectively
`(rho/A)*(gain/rate)` and `(A/B)*(gain/rate)`. Main acoustic compression
and momentum retain their original invariants. The same native stamping
owner rebuilds mass matrices, reciprocal wave/diffusion operators, matched
impedances and physical port maps. Structural coordinates, regional wall
data and Maxwell histories remain attached; only the reciprocal pressure
couplings and corresponding observations rescale.

The actual changed wave plus diffusion-memory storage enters cumulative
parameter work once. Viscosity/conductivity-only changes with unchanged
density, bulk modulus and gamma produce zero instantaneous work while
changing subsequent losses. The transfer never labels coefficient work as
heat or infers mass injection, mean-pressure traction or a thermochemical
background balance. Its external owner owes those exchanges. Native callers
must rebuild the propagator after an accepted operator change; stale operators
refuse. All numerical/source/band/work/cancellation failures retain the
original duct and state.

Three focused native regressions specify independent closed-form viscous
continuation, every arm's momentum/compression invariant, independent
instantaneous dissipation, cold-operator/port equivalence, regional wall
history retention, transport-only changes and atomic late refusals.
Bounded independent source/API/physical-account review found no concrete
defect. The reconstructed final Rust sources/tests received whitespace and
source checks but could not be parsed, compiled or numerically executed
after the executor became unavailable.


### Live regional viscothermal gas-state transfer

`DgDuct::rebind_spatial_viscothermal_gas_field` updates the existing spatial
viscothermal owner at the original element boundaries, basis, Gauss points,
radii and native Jacobi-arm support. The complete replacement `AxialGasField`
feeds the original gas-mass, local Bessel, matched-boundary and pressure-port
assemblies. Every new quadrature state retains the original band/error admission.
The uniform exact-mass transfer API keeps its separate existing rules.

The main conserved quantities are the weak coordinates `Mp*p` and `Mq*q`.
For `x=M^(1/2)u`, the transfer is
`x_new=M_new^(-1/2) M_old^(1/2) x_old`; the roots need not commute.
This does not claim exact pointwise conservation under unresolved spatial
variation. Each fixed physical diffusion section separately preserves its
viscous momentum or thermal compression volume. Gas-independent Jacobi weights
cancel from the old/new scalar arm storage-root ratio. Main and ALL diffusion
storage changes enter one signed external parameter-work account, without
subtracting unchanged large structural storage or advancing the physical clock.

Original structural displacement, velocity and regional Maxwell coordinates
remain unchanged. Reciprocal pressure/wall coupling is transformed with the
same pressure map and its opposite transpose; physical wall observation retains
its original geometric map. Rebuilding the gas must not reinterpret those
unchanged wall coordinates with the new gas mass root. The state keeps its
initial energy and prior work/loss history, and receives the new operator
identity and current physical pressure/flow.

Changed element subdivision, quadrature/pole/thermal support, invalid source
states, unresolved band, stale state, excessive/unrepresentable parameter work
and cancellation refuse before publication. The source owner prescribing the
background update owes the recorded work. This is not a mean-flow,
thermochemical heat/species evolution, pressure-prestress, remeshing or history
projection model.

Two focused native regression definitions independently integrate weak
compression/momentum and main energy at 32 Gauss points; check every physical
arm's stored quantity and energy, noncommuting roots, current-pressure wall
force, dissipation and continuing motion; and exercise uniform/spatial
observation retention and late refusal rollback. Changed source and test files
passed Rust syntax parsing and source/API review in the implementation
environment. Native compilation and execution were unavailable (no Rust/Cargo,
DSR or RCH); these definitions are not executed or experimental evidence.
Focused native command:
`cargo test -p fs-duct --lib dg::viscothermal::spatial::update::tests`.

## Experimental nonlinear Euler propagation

Feature `nonlinear-euler` exposes `nonlinear::NonlinearDuct`, an inviscid
finite-volume propagation image of the existing `Duct`. It consumes the
original source-owning `fs_gas::dynamics::FrozenGas` and delegates every
accepted step to the shared `Euler1d` solver. Density, momentum, total energy,
temperature-dependent heat capacities, shock heating and bidirectional
reflections are evolved together. No separate instrument-specific gas law is
introduced.

Each original segment boundary is a mesh face. Cylinder and cone cell volumes
are integrated analytically. Profile areas use the original radius evaluator;
four-point Gaussian quadrature is split at original spline knots and refined
within a 256-subdivision cap until repeated integrals agree to roundoff.
This is a bounded numerical estimate for nonpolynomial profiles, not a
certified enclosure. Physical face areas and integrated volumes feed the
shared pressure-area balance and local CFL calculation.

Construction requires explicit source gas, per-segment mesh counts, physical
boundaries, and initial absolute density/velocity/pressure. Every call to
`step(Time, StepControl, cancelled)` accepts at most one CFL-limited step;
refusal/cancellation leaves that step's original state and time unchanged.
Construction and snapshot allocations are fallible. The gas solver reuses
its numerical workspace.

`sample` and `samples` expose pressure, temperature, density, sound speed,
axial velocity, volume flow and mass flow at original geometric cell centres.
These local diagnostics are not a substitute for accepted boundary fluxes or
the gas solver's extensive conservation report.

Empty geometry, invalid counts, tone holes, and abrupt radius steps refuse.
A side branch or abrupt junction needs its own conservative interface law.
This image does not inherit linear DG's exact exponential or high-order
certificate, viscothermal memory, radiation-impedance boundaries, wall
compliance, multidimensional modes or an entropy-stability theorem. It
contributes to `frankensim-irp-v5-instrument-realism-pexfp.91` without closing
that full nonlinear DG/physical-loss program.

Eight focused integration tests cover exact conical volumes, exponential and
spline profile integrals, tapered rest, Sod propagation, finite-amplitude wall
reflection with reversed volume flow, evolving source-dependent heat capacity,
source identity retention, unsupported topology and atomic cancellation.
Focused native checks compile the complete gas solver and the actual duct
geometry source; a full `fs-duct`/DSR run remains a separate verification gate.

Run the source-backed, finite-amplitude tapered-horn consumer with:

```bash
cargo run -p fs-duct --features nonlinear-euler --example nonlinear_horn -- 160
```

It emits solved gas observables as CSV and conservation diagnostics on stderr.
The example's nitrogen source association is explicitly declared; it is not
an experimentally validated horn or shock-tube prediction.

Focused native execution of this consumer at 32 and 160 cells reached 0.0024 s
in 89 and 435 accepted steps respectively. Both stayed within the declared
200--1000 K source domain. The 160-cell run's mass residual was zero and its
physical total-energy balance residual was -1.421085e-14 J. These checks used
the complete wall-enabled gas source and exact production duct geometry;
they are execution/conservation evidence, not a spatial-convergence study.

## Full multimodal flanged-mouth radiation

`flanged::flanged_modal_impedance` exposes the existing finite-basis Zorumski
aperture matrix as row-major `P=ZU` in the signed `J0/|J0|` radial basis,
with impedance in Pa s/m³. Signed congruence and one aperture-area conversion
retain mutual work; peak-phasor mean radiation power is `Re(U^H P)/2`.
`modal::{mm_input_impedance_flanged,mm_profile_impedance_flanged,
mm_profile_field_flanged}` load every retained radial mode through the
original staircase and continuous propagators. The continuous boundary does
not require infinite characteristic impedance at a transverse cutoff.

Interior propagation and exterior radiation retain separate coherent gases
at equal static pressure. The exterior ka ceiling, flanged mode ceiling and
original axial-step/refinement limits remain mandatory. The continuous
half-order indicator rebuilds the radiation matrix at that smaller order.
One retained mode is the baffled piston; the legacy scalar `FlangedOpen`
continues to use its condensed higher-mode correction. This finite-mode,
infinite-baffle axisymmetric operator does not add independent higher-mode
wall-loss physics or establish arbitrary-geometry radiation accuracy.

`modal/flanged_tests.rs` checks signed mutual power against an independent
Rayleigh far-field integral, physical cylinder transfer through both original
propagators, full field work conservation, cutoff continuity and a 1/2/4-mode
resonance-frequency/linewidth ladder. Mode and axial refinement remain
separate convergence axes; a single retained order is not continuum truth.

## Finite-amplitude manifold consumer

`NonlinearDuct::into_solver` transfers the exact original geometric mesh,
source and accepted conservative field into the shared gas-network owner,
without resampling. `nonlinear_manifold` connects a cylindrical feeder and
cylindrical/conical daughters through a finite rigid 30 mL plenum, retaining
the existing NASA-9 nitrogen source and reflecting exterior ends.

Run with:

```bash
cargo run -p fs-duct --features nonlinear-euler --example nonlinear_manifold -- 64
```

The example emits branch fields and compact-volume thermodynamics as CSV,
plus common-step and conservation diagnostics. At 64 cells per branch it
executes 527 steps to 4 ms with maximum stage CFL 0.3627042, no halvings, and
zero reported aggregate mass and physical-energy residuals. Its two daughter
snapshot peak overpressures are 12.108 kPa and 7.575 kPa; temperatures range
from 295.107 to 388.735 K. These are actual execution observations, not a
resolution-converged or experimentally validated manifold prediction.

The existing eight nonlinear duct tests and the 160-cell horn example also
pass against the network-enabled gas source. Checks compile the complete
gas implementation and exact production duct geometry source slice; the
full workspace/DSR verification boundary remains explicit. The plenum is
an acoustically compact storage assumption and introduces no empirical
junction-loss or multidimensional-mixing claim.
## Harmonic molecular and wall propagation

`molecular::MolecularDuct` binds existing cylinder/cone/profile geometry to the
original `RelaxationState`. It consumes `caloric_response` once per frequency
and combines molecular heat capacity with the same positive circular thermal
diffusion spectrum used by the DG owner. Shear diffusion retains `ZkDiffusive`;
there is no added ISO, bulk-viscosity or Stokes--Kirchhoff loss. Slip/adiabatic
walls recover the material's exact relaxing wave; empty reservoirs recover
the existing lossless or finite-spectrum Bessel law. One isothermal compliance
remains in the low-frequency wall limit.

Distinct radius spectra are prepared once and reused because the admitted
gas is homogeneous. Existing stable two-port composition solves inlet
impedance and outlet pressure/flow, including strongly attenuating lines.
Cone slices retain the local-radius spherical approximation and profiles
retain their fourth-order Magnus axis. The explicit load owns exterior
radiation; no exterior medium is inferred. Molecular tone-hole chimneys use
explicit mouth loads as described below.
No mean-flow, evolving-composition or multidimensional claim is added.
Source-owned deformable walls are supported by the explicit binding below. Wall band admission is the original sampled equilibrium Bessel
comparison, not a bound for the full coupled molecular response.

`molecular/tests.rs` covers closed-form molecular waves, the original DG
operator and its power balance, empty-reservoir and isothermal limits,
pressure/rate changes, split invariance, attenuation and spatial refinement,
plus budget/band/cancellation refusal. The `molecular_duct` example accepts
the existing gas and relaxation grammars and reports solved impedance,
through-flow and input/output power across six frequencies. Its terrestrial-air
invocation composes existing engineering gas and relaxation models; it does
not establish experimental calibration.

### Molecular finite tone-hole chimneys

`MolecularDuct` now admits open, closed and partially vented compact tone holes.
The original inner, series and mutual junction corrections retain their
existing owners. Finite chimneys use the same original gas reservoirs and
cached radial wall spectra as the bore. Every open/vented mouth needs an
explicit finite passive pressure/volume-flow load, in authored open-hole order;
closed holes need no exterior load. `response_with_holes` and its controlled
variant return physical junction/mouth pressure and flow per unit inlet flow,
plus independent branch-input, exterior-load and dissipated work coefficients.
The existing partial-vent admittance weights scale flow and work once.

Decaying homogeneous cylinder coordinates use expm1 and sine/cosine identities
to preserve low-frequency compliance, quarter-wave pressure nodes and strongly
attenuating spans without growing exponentials or impedance-pole division.
The original bounded scattering owner reconstructs incoming waves with its
logarithmic transmission products. Prepared chains validate the existing
composition options before assembly; cancellation and complete block caps
apply to each junction and original propagation block. No loss is clipped.

Independent tests cover finite-line open/closed/vent responses, integrated
molecular volume work, compact gas-spring limits, three-rung wall refinement,
long attenuating chimneys, source sensitivity and explicit-load admission.
The scalar response API retains compact junctions, plane-wave chimneys and
independent exterior mouth loads. Partial venting retains the existing
admittance-mixture model; resolved pad motion and exterior radiation data are
not inferred. The field and coupled-matrix APIs extend this path as below. Original diffusion accuracy and geometry limits
remain separate from molecular-rate data quality.

### Molecular propagation through primitive gas fields

`MolecularDuct::from_gas_field` binds the original `AxialGasField` and declared
`RelaxationAbsorption` to the existing geometry/field Magnus quadrature. Local
temperature and composition are resolved through the thermochemical material
owner at each original sample and exact tone-hole station. Density, transport,
heat capacities and rates are not independently interpolated. The field covers
the exact duct length at common static pressure; geometry boundaries, source
stations and profile knots remain integration boundaries. Uniform fields
delegate the original homogeneous constructor and preserve its arithmetic.
`state()` is the exact inlet state; `outlet_state`, `gas_field` and
`field_relaxation_law` expose the retained binding.

Frequency preparation evaluates each distinct caloric state once. Exact state
and gas/radius caches use ordered maps. One `ZkDiffusiveBasis` supplies the
dimensionless radial spectra for the nonuniform field; each local binding keeps
the original complex/loss Bessel band admission. Wall poles and axial steps are
separate refinement axes. Chimneys use their local gas, including density in
the compact series correction. A varying primitive field between multiple open
holes still refuses: the existing internal mutual-junction law is homogeneous.

Focused tests compare with independent mass/momentum integration and a
modified-Bessel radial heat/momentum oracle, exercise three-rung axial and wall
ladders, and check source sensitivity, uniform identity, local chimney work,
cache reuse, coverage, caps and cancellation. Scope is frozen, no-mean-flow,
rigid-wall plane-wave acoustics; no gas transport, condensation, evolving
reservoir transfer, resolved pads or three-dimensional qualification is added.

### Molecular multi-aperture exterior feedback

`MolecularDuct::prepare_exterior_ports` exposes authored open/vented mouths,
followed by the terminal/bell, as one frequency-prepared scattering operator.
The inlet is port zero. `MolecularPortScattering::response` closes the exterior
ports with the existing frequency-specific `network::CoupledTermination` in
physical pressure/outward-volume-flow coordinates: p = Z_rad U. Self and mutual
entries participate in the loaded solve. The matrix is never inverted,
symmetrized or clipped; singular, rank-deficient and zero passive loads are
admissible. A rigid infinite terminal remains on the scalar response route.
The original network owner checks the complete Hermitian part. Exact frequency,
complete unique labels and an explicit source identity bind the input. Matrix
label permutations preserve physics and observations retain prepared order.

Original decaying bore/chimney blocks and homogeneous sealed-stub coordinates
remain authoritative. Consecutive bore spans are reduced before aperture
composition, so the dense closure has at most 33 unknowns (inlet plus 32 exterior
ports), independently of axial and diffusion resolution. Existing block and
condition limits and bounded cancellation apply. Component scattering retains
reciprocity; closed-system reciprocity also depends on the supplied exterior
matrix. BEM collocation reciprocity error is preserved rather than repaired.
Physical vent flow carries sigma once; recovery of a legacy full-radius scalar
mouth load Z requires Z/sigma on the corresponding physical matrix diagonal.
This retains the original compact vent-admittance mixture.

The solution reports all mouth pressures/flows, total exterior work and raw
input-minus-exterior dissipation, including bore/chimney molecular and wall
heat. Individual mouths can receive work under mutual coupling; only summed
exterior work must be passive. Focused tests cover independent finite-line
pressure/flow closure, scalar/vent recovery, reciprocal lossless power balance,
sealed-stub poles, opaque chimneys, rank-deficient loads and admission limits.
The `molecular_exterior` example consumes explicit gas/relaxation declarations
and a physical three-mouth impedance matrix file and reports solved pressures,
flows and powers. It adds no default exterior medium, BEM kernel, nonlinear
oscillation loop, radiation fit, resolved pads or time-domain realization.

### Molecular gas with source-owned compliant walls

`MolecularDuct::{new_with_compliant_wall,from_gas_field_with_compliant_wall}`
bind an `Arc<dyn WallLaw + Send + Sync>` beside the existing molecular gas and
optional radial heat/shear spectra. The shared `WallModel` constitutive map
adds `2*pi*r*sqrt(1+(dr/dx)^2)/Z_surface` to the physical gas volume admittance.
Molecular relaxation, thermal diffusion and displacement therefore act in one
mass/momentum operator; wall displacement is not another gas thermal shunt.
The original rigid constructors retain their wave arithmetic and cone route.
The explicit compliant route makes `ConeSlices` live even without radial
diffusion. Profiles and primitive gas fields use their original quadrature,
including the actual surface slopes. Radius/material/slope response caching
cannot alias different physical wall areas.

Every existing scalar, finite-chimney, flow/work and coupled-exterior response
uses that same immutable binding. Closed stubs retain wall loading, vent flows
retain the original single opening weight, and exterior self/mutual matrices
still act on physical outward volume flow. Passive finite nonzero impedance
and source applicability are admitted by the original wall-law boundary at
each queried frequency and radius; cancellation and composition caps remain.
Input work minus exterior work includes gas loss and work delivered to the
mechanical wall. A wall law may include radiation, so this balance is not
exclusively gas heat. Local wall reaction does not resolve axial shell modes,
nonlocal wall coupling, large displacement or thermoelastic heat exchange.

Focused tests compare a cylinder to independently assembled Debye plus radial
Bessel heat/momentum equations and a spatial work integral; refine mechanical
cone slices without thermal diffusion; compare a varying primitive field and
profile against independent physical-equation RK4; preserve uniform-field
bits; and compare finite chimney flows and scalar/mutual exterior closures.
Source-domain/passivity refusals, work caps, cancellation and repeated-radius
response caching remain observable.

## Continuing uniform molecular gas backgrounds

`DgDuct::rebind_uniform_molecular_gas_field` continues the original molecular/thermal coordinates through a prescribed change of uniform background gas. Source-backed field attachment retains its complete `RelaxationAbsorption` declaration, including the uniform specialization; `molecular_relaxation_law()` exposes that declaration. The existing `molecular_relaxation_field()` still identifies the nonuniform pointwise route. Evaluated-only attachments cannot invent a missing source law.

Both fields must remain uniform, cover the exact original geometry, and retain identical complete species fractions and evaluated reservoir support. The original law reevaluates every material coefficient. Native exact-mass circular diffusion, positive spatial wall quadrature on fixed variable geometry, and adiabatic molecular gas are supported. Every spatial point, original radius, radial section count, mass basis and angular admission band survives; new local Bessel models must meet the same complex-response and real-loss gates. Inactive species retain their original source declaration and the material owner's original omission semantics.

The intervention preserves weak acoustic compression, axial/viscous momentum, and each thermal/molecular mode's physical entropy perturbation. With `t_n=y_n+a_n P_n x_p`, compression is proportional to `x_p - sum a_n P_n^T y_n`, so pressure includes the history-dependent caloric correction. Thermal and reservoir entropy densities scale with `sqrt(rho cp_inf/T)` and `sqrt(rho C_j/T)`. This is an explicit prescribed intervention without injected acoustic entropy or compression, not a background heat, inventory or species transport solve. The external material owner owes the signed change in actual quadratic storage as parameter work.

Structural displacement, velocity, viscosity and Maxwell coordinates retain their existing owner and histories. Accepted sample count, initial storage, port work and dissipated energy survive. The new operator and state publish atomically after finite/work/cancellation checks; a changed operator invalidates its old propagator. Unchanged material preserves coordinates and operator bitwise; an empty molecular declaration delegates to the original nonmolecular transfer. Nonuniform background refresh, species-fraction changes and changed reservoir support refuse.

Focused native tests under `dg::molecular::update::tests` cover independent physical entropy/compression/momentum reconstruction, a conical spatial-wall weak equation-of-state check, nonzero hereditary wall continuation and energy balance, round trips, exact no-op/empty-law parity, and refusal/cancellation rollback.
