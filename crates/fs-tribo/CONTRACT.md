# CONTRACT: fs-tribo

> Status: `frankensim-ext-tribo-dry-baseline-tgbj` baseline. This is not an
> admitted material-data result or an Euler-disc validation result.

## Purpose and layer

`fs-tribo` is an L3 dry-contact constitutive leaf. It provides typed dry
friction, elastic Hertz reference forms, scalar rolling/contour loss laws,
explicit dissipated-work accounting, a bounded flash-temperature candidate,
and caller-owned Archard wear state. It does not depend on contact-solver
internals: the consuming solver owns contact complementarity, tangent
reactions, and residual embedding.

## Public types and semantics

- `InterfaceSystemRef` names ordered surfaces, history, medium, and a nonblank
  caller `source_id`. `InputAuthority` is a caller-declared ceiling
  (`CallerDeclared`, `SyntheticFixture`, or `Estimated`), never an admission or
  verification receipt. Every constitutive response returns this provenance.
- `DryInterfaceSystemCard` is a dependency-independent, immutable
  InterfaceSystemCard-style input. It combines that ordered interface/history
  identity with one `FrictionLaw` and closed SI temperature, nominal-pressure,
  and slip-speed ranges. Its `query` refuses all out-of-domain states and a
  mismatch between the declared slip-speed magnitude and the supplied tangent
  vector. It is not an alias for, or a substitute for, the pending seed-data
  `fs-matdb::InterfaceSystemCard` query/receipt path.
- `ContactFrame` normalizes a finite contact normal. `TangentialSlip` accepts
  only a finite velocity with no material normal component; it refuses rather
  than silently projects a closing/separating velocity into friction.
- `FrictionLaw` provides Coulomb, velocity-dependent, and Stribeck rungs. At
  zero slip `evaluate` reports static capacity only; it never invents a stick
  reaction. `kinetic_coefficient` exposes `μ(|v|)` for a strictly positive
  speed. `regularized_traction_1d` is the declared ODE embedding: a linear
  ramp through the origin to the kinetic value at scale `v_reg`, then the kinetic rung,
  with force opposing `v_rel`. That map is a solver regularization (the 1-D
  analogue of the fs-contact smooth adapter), not complementarity. A bow, a
  brake pad, and a fault are the same law.
- `TangentialLaw` is the stateful, solver-independent pHS constitutive seam.
  `FrictionLaw` and `DryInterfaceSystemCard` implement it with an implicit static
  reaction followed by deterministic sliding bisection, never an ODE ramp.
  Input motion is the released driver-minus-body increment and the explicit
  isotropic body compliance gives `du = du_free - C F`. Returned traction is
  effort entering the interface, conjugate to the **solved** increment, so
  `Delta H + D - F dot du = 0`. Dry laws have `H = 0`, export `D >= 0` through
  `DissipationStep`/`WorkLedger`, and retain accumulated slip plus the optional
  caller thermal boundary. They do not evolve temperature or a third-body layer.
  `TangentialSolveConfig` requires explicit force/residual/iteration budgets
  (1 through 128 bisections); unresolved or unbracketed roots refuse. Candidate
  state is immutable until the consuming solver publishes it. The card impl
  requires temperature and pressure without defaults, and checks the solved
  speed against its closed applicability interval. Anisotropic body response,
  thermal transport and pHS network wrappers belong to their respective owners.
- `layer_contact::ThirdBodyLayerContact` implements the same seam with retained
  elastic and plastic shear. The midpoint radial return solves the full
  isotropic body load line in two tangent axes. Reversible energy is
  `H = k |e|²/2`; physical dissipation is `Y |dp| + c |dp|²/dt`, with no
  numerical energy deficit relabeled as heat. `Sticking` means no plastic flow;
  elastic relative displacement may remain nonzero. Every step requires
  positive compressive load, the original frozen temperature, the caller's
  explicit elastic-slip bound and force budget. Opening and foreign-law
  checkpoints refuse without changing state. The numerical law has caller
  authority; the separate fs-material source owner binds actual layer geometry
  and the distinct plastic-shear-viscosity claim. Thermal softening, layer
  transfer/wear, recontact and measured rosin behavior are outside this rung.
- `HertzSpherePlane` and `HertzCylinderPlane` are G1 elastic closed forms,
  expressed in SI. They use caller-supplied radius/modulus and retain the
  caller ceiling; a mathematical closed form does not validate its inputs.
- `ConstantRollingMoment` uses `M = -sign(omega) N a` (N m) and
  `P = N a |omega|` (W). `ConstantContourForce` uses
  `F = -sign(v) F_c` (N) and `P = F_c |v|` (W). They are separate mechanisms;
  callers must not combine them without independent source support.
- `HeatPartition`, `DissipationStep`, `WorkLedger`, and `WearState` have
  private invariant fields and checked getters. `WorkLedger` retains total
  dissipated work and every explicit heat/other channel, and commits all four
  candidate totals only after their closure checks pass. `ArchardLaw` computes
  `dV = k F ds/H` in m3 and commits only a validated candidate state.
- `flash_temperature_candidate` is a typed model-form branch. With complete
  caller thermal conductivity/diffusivity data for every nonzero heat share,
  it reports the uniform-flux semi-infinite candidate
  `2 q'' sqrt(alpha * t / pi) / k`, with `t = traverse_length/slip_speed`.
  Missing properties return a typed `Unknown`, rather than guessed material
  data or a partial temperature result.
- `surface_excitation` admits named, uniformly sampled finite or periodic
  surface-height traces and filters them by the exact one-dimensional marginal
  of an elliptic Hertz pressure footprint. It analytically integrates the
  piecewise-linear trace and its path derivative, then applies only the
  explicitly bounded contact tangent `delta_F = k_n delta_h`.
  `SelfAffinePeriodicProfileSpectrum` can construct a deterministic,
  random-phase, band-limited one-dimensional profile whose spatial PSD follows
  `Phi(q) proportional to q^-(2 H + 1)` and whose Fourier coefficients are
  normalized to an explicit RMS height. Its seed selects spatial phase only;
  it creates no temporal/audio spectrum or apparatus-specific frequency.

### Partial-slip return-map rung

`partial_slip` is a solver-independent, two-region Jenkins/Iwan
*return-map* rung (`two-region-iwan-frame-transport-v3`). `NormalPatchView`, `PartialSlipInterface`, and
`GeneralizedWorkOwnership` retain the caller's named normal-patch, ordered
interface/history, and work-owner inputs; their authority remains explicitly
caller-declared. `PartialSlipLaw::advance` has no `fs-contact` dependency and
cannot solve normal contact, evolve patch geometry, or admit a material card.
Two parallel regions retain independent elastic shear/twist and kinetic-branch
history. Their stiffness fractions are `1-h` and `h`, their Coulomb-capacity
fractions are `a(1-h)` and `1-a(1-h)`, and `0<a,h<1`. Holding the displacement
preserves preload without relaxation or heat; reversal first unloads elastically.
One yielding region means partial slip, both means gross slide. The mixed
translation/spin return satisfies maximum-dissipation normality in the elastic
compliance metric; pure translation or spin uses a closed-form return.
This finite two-region approximation is not a resolved Cattaneo--Mindlin
pressure/traction field or a measured slipping-area fraction.

Work is independently computed as `-F dot du - M dtheta`. The balance separates
reversible storage, physical plastic-slip heat, and the backward-Euler spring
loss `sum(k_i |Delta e_i|^2 / 2)`. Numerical loss is never thermalized. Its
refinement claim applies to smooth fixed-branch histories; unresolved
static-to-kinetic jumps require event resolution and have no such claim.
The endpoint wrench is unchanged by energy bookkeeping. The work is a discrete
endpoint contraction, not an exact continuous-time integral or midpoint work.
Purely elastic closed cycles generate no physical heat, while their numerical
loss decreases at first order as the interval is refined.
Retained shear is transported before applying the next relative displacement.
The declared geometric convention is the minimum rotation between successive
contact normals, followed by the supplied common material spin about the new
normal. The final tangent-gauge change only re-expresses coordinates. A single
normalized plane rotation acts on both regions, preserving spring lengths;
normal transport therefore introduces neither projection damping nor heat.
Relative spin alone changes stored twist. Common spin is a distinct mandatory
kinematic input, `(omega_a + omega_b) dot normal / 2` for two rigid surfaces;
it cannot be inferred from an arbitrary tangent gauge. Normal tilt and spin
are composed as an interval rotation, with the same time-resolution limitation
as the surrounding constitutive step. Nearly antipodal normal jumps refuse
because the minimum-rotation path is unresolved; supply intermediate geometry.
Checkpoints contain both regional strains, kinetic branches, and the last
accepted frame, plus all law, patch and interface identities. A nonempty
frame-free or aggregate-strain checkpoint is insufficient. Equality of those
inputs is not external physical validation.
Rolling-deformation loss remains a separate, zero channel in this rung.

### Rolling and contour-deformation loss rungs

`rolling_loss` provides three immutable, separately sourced candidates:
the `Leine-style Coulomb contour-force rung` (`F_c = -sign(v_c)c_cN`), a
viscous contour alternative (`F_v = -c_vv_c`, with `c_v` in N s/m), and a
finite-patch hysteretic rolling moment (`M_h = -sign(omega_r)N l_h eta_h`).
The last card declares its loss length `l_h` [m] and bounded loss factor
`eta_h`; it is a model-form candidate, not a resolved deformation field.
Every response retains a named patch/load/curvature receipt, ordered
interface/history, temperature/frequency range, source-card authority, and a
dedicated generalized-work key. Principal curvatures are retained directly;
an equivalent-radius substitute requires explicit caller approximation
authority. The rungs have no blend/preference operation, no numerical error
certificate, and no torsional microslip or gross-slide contribution. Pure spin
therefore produces no torque in this slice; a separately sourced torsional
material-hysteresis rung would need a disjoint work channel.

## Invariants

- All public construction and evaluation paths refuse missing identity,
  non-dry medium, negative/non-finite input, and every non-finite derived
  candidate from finite extreme inputs.
- Dry tangential and scalar resistance oppose their declared relative rates.
  Elastic layers can return stored work; total contact work equals storage
  change plus finite, non-negative physical dissipation.
- Heat shares are finite, non-negative, and sum to one. A dissipation step's
  channels are finite, non-negative, and close to its total before a ledger can
  mutate. All ledger channels and wear candidate totals are checked before any
  assignment.
- Card state ranges are closed and finite; absolute temperature has a positive
  lower bound. A direct `FrictionLaw` remains a low-level constitutive rung,
  while state-domain enforcement occurs at `DryInterfaceSystemCard::query`.
- A flash result is either complete for both receiving surfaces or an explicit
  `Unknown`; it never silently assigns missing thermal properties.
- Partial-slip refuses nonfinite or nonpositive patch and law inputs, malformed
  tangent frames/kinematics/work ownership, and checkpoint mismatches. It
  computes only its declared lumped return-map state and never creates a
  contact-solver or material-admission receipt.
- Each rolling-loss card refuses non-dry interfaces, malformed patch/source
  receipts, out-of-range temperature/frequency, nonfinite derived work, and a
  wrong work channel. A rolling key conservatively refuses sharing the same
  patch/interval pair with partial-slip because that module does not expose its
  coordinate labels for a finer disjointness proof.
- Rolling/contour candidates are passive: endpoint body power is non-positive,
  heat is non-negative, reversible storage is zero, and the retained work
  identity reconstructs heat over the caller interval.
- No mutable material table, hidden state, or unordered material-pair lookup
  exists. Ordered interface and history identities remain explicit.
- Surface excitation refuses missing texture/source identity, uncovered finite
  footprints, excessive periodic repetition/work, non-finite kinematics, a
  height outside the caller's tangent-linearization fraction, and a
  perturbation that would open the unilateral contact. Both surface sources
  and authority ceilings remain separate in the receipt.
- `PeriodicHarmonicSurface` admits either a measured Fourier decomposition or
  an explicitly declared hypothetical periodic material track. Signed
  sine/cosine height coefficients are keyed by positive integer spatial cycle,
  canonicalized by cycle count, and sampled only with at least eight points per
  shortest period. Duplicate cycles and underresolved grids refuse. An empty
  spectrum is an exact smooth track. No material name, acoustic frequency,
  stochastic realization, seam correction, or authority upgrade is inferred.
- `SelfAffinePeriodicProfileSpectrum` refuses nonpositive/nonfinite RMS height,
  Hurst exponents outside `0 < H < 1`, inverted or zero cycle bands, and
  harmonic bands whose shortest wavelength cannot satisfy the eight-samples-
  per-period rule within the one-million-sample trace ceiling, as well as
  unbounded harmonic counts. Every retained cycle owns a counter-based phase
  stream keyed by the explicit seed and cycle, so extending the band preserves
  shared-cycle phases. Phase changes preserve the declared Fourier power, and
  the continuous periodic profile retains the requested RMS height to
  floating-point roundoff.

## Quasistatic Reynolds film

`resistive_film` owns the solver-independent incompressible thin-film pressure
law shared with the compatibility surface `fs-flux::resistive_film`. At each
configuration it solves the bounded pressure island `L(h) p = A B v` and applies
the reciprocal resisting force `B^T A p`. The analytic tangent includes both
gap-dependent conductance and independent mechanical effort changes. Cell areas,
port weights and Poiseuille conductances are never independently normalized.
This law carries no clock or pressure history; the mechanical caller owns
stepping, cancellation between calls and integrated energy accounting.

The fixed limits are 64 pressure cells, 256 channels and 256 mechanical ports.
Nonpositive channel apertures close exactly; collapsed volume cells, trapped
pressure pockets, pressure/gap limit violations and unresolved factorizations
refuse without modifying output slices. No leakage, gap floor, pressure clipping
or diagonal repair is introduced. Compressibility, thermal/slip effects, fluid
inertia, evolving topology, fluid radiation and elastohydrodynamic closure are
outside this law. The existing nine pressure-law tests run through the fs-flux
compatibility surface, while fs-couple tests exercise reciprocal mechanical
reaction, joint dissipation and atomic retry in its existing time integrator.

## Error model

`TriboError` is a total refusal surface for missing identity, non-dry media,
invalid/overflowing physical input, non-finite vectors, normal slip, malformed
or out-of-domain applicability ranges, invalid partitions, and forged/invalid
dissipation states. A refusal makes no partial state change to `WorkLedger` or
`WearState`. Missing flash thermal data is a typed `Unknown`, not an error or a
candidate.

## Determinism class

Pure scalar operations are deterministic for the same input and ISA. Every
solver-state transcendental is routed through `fs_math::det` (bead
3ez8g.7.3): norms use `det::hypot`, which is overflow/underflow-safe BY
CONSTRUCTION (max-scaling, so the ratio never exceeds 1) — the original
platform-`hypot` rationale is preserved, and the kernel is bit-identical on
every conforming target; the Stribeck decay uses `det::exp`, tangent-frame
rotation `det::sin`/`det::cos`, and the Hertz `x^{3/2}` is the IEEE-exact
`x*sqrt(x)`. CROSS-ISA BIT-STABILITY IS CLAIMED for these friction paths and
enforced by a single committed golden hash executed on both aarch64 and
x86-64 (`cross_isa_golden_hash_over_the_routed_friction_paths`), plus the
workspace `check-libm` lint (fs-tribo is a doctrine crate). The recorded
counter-argument: platform `hypot` is near-correctly-rounded while
`det::hypot` spends up to 2 ULP — a cost orders below Stribeck-parameter
uncertainty, accepted for reproducibility. Composed images (e.g. the bowed
string) remain capped by their OTHER operands until audited.

## Cancellation behavior

This bounded scalar leaf has no asynchronous work and no cancellation scope.
Callers scheduling many contact points own cancellation between calls.

## Unsafe boundary

None. The crate forbids unsafe code.

## Feature flags

`benchmark-laws` (default off) exposes `fs_tribo::benchmark`, the M4.2-W4
benchmark tangential laws. They are also compiled under `cfg(test)`. Only
dev-dependencies may enable the feature;
`tests/benchmark_laws.rs::g0_benchmark_laws_not_reachable_from_render_paths`
scans every workspace manifest and fails if a `[dependencies]` table enables
it.

### Benchmark tangential laws (`benchmark-laws`)

These laws exist for named comparisons, such as the X-Bw2 bow-law bakeoff.
They are never render laws.

`RateState` is Dieterich–Ruina rate-and-state friction.
- Friction: `mu = mu0 + a ln(V/V0) + b ln(V0 theta/Dc)`.
- `Evolution::Aging` is `d theta/dt = 1 - V theta/Dc`.
- `Evolution::Slip` is `d theta/dt = -(V theta/Dc) ln(V theta/Dc)`.
- Admission: `a > 0`, `b >= 0`, positive `V0`/`Dc`, and a declared `speed_range_mps`. Outside that range a step returns `OutsideApplicability`. A zero-speed step refuses because `ln V` is undefined.
- Over one held-effort interval, both state laws integrate exactly in closed form using `expm1`.
- Friction uses the state at the start of the interval. It is then strictly increasing in `V`, so the compliant root is unique.
- `critical_stiffness_n_per_m(N) = N (b - a)/Dc`.

`LuGre` (Canudas de Wit et al. 1995) is one-dimensional.
- Bristle law: `dz/dt = v - sigma0 |v| z/g(v)`, where `g(v) = Fc + (Fs - Fc) e^{-(v/vs)^2}`. It is integrated exactly at held `v`.
- Force: `F = sigma0 z1 + sigma1 (z1 - z0)/dt + sigma2 v`.
- `sigma0 z^2/2` is reversible storage, and the remainder of the port work is heat.
- A step whose dissipation would be negative refuses, rather than assuming the `sigma1` passivity condition.
- A two-dimensional slip increment refuses.

Both laws implement `TangentialLaw`:
- The compliant branch bisects with the caller's `TangentialSolveConfig` budget.
- A rate-and-state root that needs a speed outside the admitted range returns `OutsideApplicability`.
- Parameter caps are checked at `cap` and at the next representable value above it.

`SmithWoodhouse` is the plastic thermal model of rosin from Smith & Woodhouse
(2000), J. Mech. Phys. Solids 48:1633-1681. Equation numbers are the paper's.
- Friction (Eq. 21): `mu N = A k_y(T) sgn v`, set-valued in sticking, with `A = pi a^2` proportional to `N`.
- Temperature: the heat balance of Eq. (6) per unit area. Conduction into the stationary and the moving solids uses the discrete half-space Green's function of Eqs. (10)-(14). The moving solid sees only the history within slip `8a/(3 pi)` (Eq. 16 at constant speed), with the cut pulse weighted exactly.
- `k_y(T)` is calibrated from the steady-sliding fit (Eq. 2) through Eq. (20) (section 5.3). `SmithWoodhouse::new` certifies with an interval bound that the calibration temperature increases with speed, and refuses otherwise.
- Time stepping: backward Euler in the losses, with friction from the step-start temperature, over a declared uniform step. The paper uses a forward difference.
- Refusals: a different step length, zero load, a temperature above the calibrated range, a history beyond the declared capacity (at most `MAX_HISTORY_STEPS`), and non-positive parameters.

## Conformance tests

Inline and external G0/G1/G3 tests cover identity/dry/tangent refusals;
ordered-card applicability and slip-speed consistency; analytic rigid-block
equilibrium; independent numerical Hertz values and force-pressure-radius
cross-relations; reversal/scaling/provenance; resistance sign; partition and
ledger-channel closure; typed flash-data insufficiency and power scaling;
forged negative/non-finite work rejection; and rollback/deterministic replay
for work and wear. Surface-excitation tests additionally cover periodic-spectrum
order invariance, exact sample-grid realization, duplicate-cycle refusal, and
the eight-samples-per-shortest-period gate. They also cover deterministic
self-affine replay, exact RMS normalization, the declared PSD slope, seed-only
phase changes, and malformed statistic/band refusal. Test coefficients are
explicitly synthetic fixtures.
The separately owned partial-slip tests cover its scalar return-map admission,
reversal, checkpoint replay, and receipt-mutation refusals.
`tests/tangential_law.rs` adds G0 DEF-06 continuity at `v_reg*(1 +/- 1e-12)`
with a `1e-12*mu_s*N` band for all three rungs; 10,000 Philox-seeded steps
per rung with non-negative heat and per-step closure within
`1e-12*max(abs(F dot du), f64::EPSILON J)`; G1 incline acceleration within
`1e-12 m/s2`; closed-domain next-up refusal; force/iteration-budget refusal;
and G5 checkpoint replay. The actual modal owner is compared bitwise in
`fs-couple`'s `friction::scalar::trait_tests`, avoiding an upward dependency.
These are named verification targets; their execution is required before
claiming numerical qualification.
Rolling-loss tests cover zero load/rate/factor, force and moment reversal,
load scaling and SI dimensions, pure roll/spin/mixed motion, passivity with an
independent work integral, temperature/frequency refusal, work-key overlap,
source-card/patch replay mutation, and a no-target-field source scan.
`tests/benchmark_laws.rs` (feature `benchmark-laws`) covers:
- **G1 rate-and-state steady state:** the `mu_ss = mu0 + (a - b) ln(V/V0)` log law at relative 1e-10, for both evolution laws across 1e-8..1e-2 m/s.
- **G1 velocity step:** the 1/e transient decays over `Dc` of slip, within 1 % at 12.5/25/50 steps per `Dc`.
- **G1 Ruina spring-slider:** stiffness sweeps through the compliant `advance` path bracket `k_c = N (b - a)/Dc` between `0.98 k_c` and `1.02 k_c`. They compare oscillation amplitude early against late, over 300 `Dc` of slip.
- **G1 LuGre:** the steady-state Stribeck curve at relative 1e-10, and presliding stiffness `sigma0` at relative 1e-8.
- **G0 LuGre bristle map:** passivity over Philox-seeded histories (12,800 intervals).
- **G1 LuGre drive-spring:** stick-slip exists and is periodic.
- **G0:** parameter cap and cap+1.
- **G1 Smith & Woodhouse Green's function:** the weights match Eq. (13), its telescoped partial sums and Eq. (15).
- **G1 Smith & Woodhouse calibration:** Eq. (20) temperatures and Eq. (2) friction at relative 1e-12, over 1e-4..0.3 m/s for four rig geometries.
- **G1 Smith & Woodhouse steady sliding:** friction approaches Eq. (2) as `t^{-1/2}`, the stationary half-space tail. One Richardson step leaves under 10 % of the 4 s deviation.
- **G1 Smith & Woodhouse sinusoidal sliding (Fig. 18):** clockwise loops confined to low speed, which close at the reversal as `O(h)`. The three-rung step ladder has order in [0.8, 1.2]. The temperature levels are printed and not gated: they do not reproduce Fig. 18(a) or the 75-115 K stated for a 0.25 mm radius.
- **G0 Smith & Woodhouse:** caps, refusals and the set-valued stick/slip branch.
- **G0 reachability scan.**

## No-claim boundaries

The rate-and-state and LuGre benchmarks carry synthetic parameters. The Smith
& Woodhouse tests use that paper's published rosin, perspex and steady-sliding
inputs, for comparison only. No benchmark claims a material admission, bow
data, or equivalence to any render law. A benchmark result names a comparison and does
not admit a law.

This crate does not mint material admission, calibration, the upstream
`fs-matdb::InterfaceSystemCard` query receipt, a flash-temperature error bound
or thermal-port solution, a roughness measurement or statistical roughness
model, normal adhesion or adhesion-induced contact-area evolution, general three-dimensional plasticity, finite-patch partial
slip as a resolved contact field, general lubrication/EHL, contact geometry evolution,
wear geometry updates, stop time, Euler-disc ranking, one-millimetre optimum,
sound pressure, loudness, or experimental/video correspondence. The surface
trace leaf consumes a caller-supplied one-dimensional profile and a local
contact tangent; it does not replace a two-dimensional roughness map, nonlinear
re-solution of contact, structural response, or acoustic radiation. Its flash candidate is a declared
semi-infinite uniform-flux model-form estimate, not a temperature measurement
or validation result. A caller may record an estimated or synthetic result,
but must not promote it based on this crate alone.

The self-affine constructor is likewise a geometry realization, not evidence
that a specimen is self-affine or that its RMS height, Hurst exponent,
wavelength cutoffs, or phase realization are correct. Those inputs require
profilometry or an explicitly weak estimate. The one-dimensional profile does
not establish an isotropic two-dimensional areal PSD or resolve asperity-scale
partial contact.

The rolling and contour rungs likewise do not establish deformation-loss
dominance, coupon calibration validity, normal-contact admission, an error
bound, or any apparatus-specific selection between their rival model forms.

## Explicit semi-infinite heat allocation

`HeatPartition::semi_infinite_surfaces` allocates an interface-generated heat
pulse in proportion to `e = k / sqrt(alpha) = sqrt(k rho cp)`. It uses the
existing checked surface properties and heat partition. The normalized ratio
does not overflow when the sum of two finite effusivities would. This is the
equal-exposure, homogeneous isotropic half-space approximation: the two
incremental surface-temperature rises agree and the two fractions sum to one.
The caller owns contact resistance, moving-source effects, finite thickness,
material state and applicability. No third-body heat share is inferred.

G1/G3 coverage in `tests/dry_baseline.rs` checks an independent effusivity
ratio, equal incremental temperatures through the existing flash-temperature
formula, ordered-surface reversal and extreme finite normalization. The
semi-infinite boundary assumptions and finite-size limitations follow A. Jain,
*International Journal of Heat and Mass Transfer* 202 (2023), 123721,
<https://doi.org/10.1016/j.ijheatmasstransfer.2022.123721>. This is numerical
verification of a selected model, not experimental material qualification.

## Adhesive junction shear

`adhesive::AdhesiveFrictionLaw` retains shear resistance of an attached junction
at zero or tensile net normal force. Its static and sliding capacities are
`tau * A + mu * N_rep`, with independently supplied real contact area, junction
shear strengths, and the nonnegative repulsive contribution to normal force.
The signed net adhesive force must not be substituted for `N_rep`. Setting the
junction strengths to zero recovers the existing pressure law; setting its
coefficients to zero selects constant junction shear strength. Calibrations
must avoid counting the same junction resistance in both terms.

The leaf borrows the original interface provenance, allocates nothing during
evaluation, opposes sliding, and reports nonnegative frictional power. The
scalar regularization ramps the complete capacity continuously near zero slip.
`tests/adhesive.rs` checks finite zero-load friction, detachment, pressure-law
recovery, extensive scaling, vector power, reversal, regularization, input
refusals, and subnormal velocities. Contact geometry, partial slip, shear-induced
peeling, bond aging and measured strength admission remain the caller's concern.
No strength is inferred from adhesion energy or bulk material names. The
area-proportional term follows the measured single-asperity relation in Carpick
et al., *Langmuir* 12 (1996), <https://doi.org/10.1021/la9509007>; the optional
pressure term is an explicit uncoupled model choice, not a universal fit.

## Changing-coefficient partial-slip continuation

`PartialSlipLaw::advance_with_coefficients` continues the same contact/source model with changed Coulomb coefficients or translational/torsional stiffness. It retains both regions' elastic strains and kinetic branches and runs the original return map; regional stiffness/capacity fractions and the torsional capacity factor remain fixed. Its independent parameter-work port supplies `C = psi(new coefficients, old history) - psi(old coefficients, old history)` before motion. Returned storage change refers to the actual old law, so `delta(psi) + physical slip heat + numerical loss = wrench work + C`. Softening can export parameter work; reduced yield strength can relax stored strains even at zero imposed motion. Ordinary restore remains coefficient-exact.

The coupled owner must accept the coefficient driver and contact history together and balance the signed parameter port. This first-order parameter-then-motion split supplies no temperature-dependent material data or caloric closure: temperature-dependent spring free energy requires mechanical internal energy `U = psi - T d(psi)/dT`, so parameter work alone is not heat. Numerical dissipation also remains nonthermal. G0/G3 tests in `tests/partial_slip_coefficients.rs` cover signed stiffness work, independent mixed shear/spin energy, retained preload and normal-load changes, zero-motion relaxation, unchanged-law identity, and strict source/decomposition refusal.
