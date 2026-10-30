# fs-material — CONTRACT

The constitutive-law and authority-separated identifiability kernel: materials
as mathematical objects — calibration domains, CONSISTENT tangent operators,
thermodynamic guardrails, hysteresis, uncertainty — plus source-closed,
multi-case questions whose execution and assessment identities stay distinct.

Ambition tags: elastic/hyperelastic/J2/RC-fiber laws + calibration [S];
multi-case identifiability authority schema [F]. Ogden and tension-only
HGO/GOH are [S] constitutive laws; specimen/tissue calibration is separate.

## Source-bound axial modulus references

`axial_reference::ResolvedAxialModulusReference` resolves density and an
explicitly named tensile/Young modulus for tension or compression modulus for
compression from one immutable material
card and one complete typed source point. It retains both original usage
receipts and uncertainties. Tension and compression are separate modes;
dynamic or flexural modulus, strength and pointwise stress cannot substitute
for an axial modulus. A source-reported Young modulus does not establish a
tensile-test protocol, directional equivalence or a measured linear range. No Poisson ratio or three-dimensional isotropy is inferred.

This is an authored, zero-origin **linear modulus reference**. A reported
compression modulus fitted over a finite strain interval is not reclassified
as a measured initial tangent. The caller declares a strain ceiling, original
length and area. The existing generalized-Maxwell equilibrium-spring limit
supplies nominal stress and ramp work; the coupon returns `E*A/L`, mass,
force, displacement and **model loading work**.

Loading must increase monotonically in the selected direction. Holds,
unloading, reversal, excess strain and nonrepresentable values refuse before
history changes. When the selected claim declares a crosshead-speed axis, the
caller must identify it and the displacement/time ramp must match. When source
speed is unstated, ramp duration remains an authored model choice. A query's
unconstrained extra axis cannot acquire source-rate authority.

The model does not predict relaxation, creep, cyclic hysteresis, heat generation,
recoverable material energy, lateral contraction or failure. Sampled strain,
rate and history laws cannot be silently frozen into this reference.

Focused checks cover independent axial stiffness/force/mass/work references,
separate loading modes, source/protocol refusals and unchanged state after
rejected steps:

```sh
cargo test -p fs-material axial_reference
```

The `material_axial_reference` example reads a canonical source pack and derives
axis types from its selected source properties. It accepts explicit geometry,
loading mode, source coordinates and up to 64 monotonic displacement ramps.
It validates the entire requested path before emitting a successful response.
The material seed README includes complete PEEK tension and ROHACELL compression
recipes with independent force, mass and model-work expectations.


## Source-bound compression observation references

`compression_reference::ResolvedCompressionReference` binds a selected nominal
engineering-compression stress curve and a strain-independent reference density
from one immutable material card. Original typed queries, claim selections,
uncertainties and exact knot receipts survive. Source interpolation is the
default requirement. The explicit `PiecewiseLinearBetweenObservations` choice
authors a bounded linear reference between discrete observations without
changing their `TabulatedOnly` source policy or issuing an interpolated source
receipt. All original knots must remain supported; no origin or missing segment
is fabricated. Scalars, strengths, moduli, true-stress keys, sampled density and
unowned rate dependence refuse.

A uniform specimen supplies original length and area. It starts at the first
source compression and keeps its original mass. Strictly increasing shortening
produces nominal force, one-sided segment stiffness and piecewise-integrated
loading work **relative to that first state**. Known crosshead speed is enforced;
unknown rate leaves ramp duration authored. Rejected steps preserve history.
No unloading, holds, creep, recovery, cyclic hysteresis, lateral contraction,
failure, dissipated heat or recoverable-energy claim is made. The model does
not identify a source test's geometry with the caller's geometry.

`material_compression_reference` loads canonical packs and validates the entire
requested path before printing successful output. Its explicit
`--interpolate-observations` flag exposes the model-authoring choice.
`tests/compression_reference.rs` checks the source and numerical boundaries;
the catalog's compression-reference tests exercise four actual foam grades
and reinforced PEEK through the production compiler and this owner.

## Purpose and layer

Layer **L3** (FLUX support). Runtime deps: `std`, fs-ad (dual-number
energy derivatives), fs-blake3 (typed content identities), fs-evidence
(model cards, Evidence, and V&V artifacts), fs-matdb (immutable material and
constitutive-model cards), fs-qty, fs-math.
Consumers: fs-solid elasticity (tfz.13), fiber-section beams, lattice
homogenization, the P2 milestone.

## Public types and semantics

- `wood_moisture::ResolvedWoodMoisture` binds the original 12% dry-basis
  density and independent green-to-ovendry radial, tangential and volumetric
  shrinkages, plus the source-declared fiber-saturation approximation. Exact
  keys, typed source queries, card identity and all five receipts survive.
  Equation (4-9) gives dimension ratios relative to 12% moisture; independent
  volume data and water mass then give density while conserving ovendry mass.
  The admitted moisture interval intersects the coefficient domains with the
  below-fiber-saturation range. No longitudinal dimension or deformation
  tensor is inferred from independent population averages. Moisture is
  prescribed; this does not solve drying, diffusion, hysteresis or stress.

- `section_temperature::thermoelastic::ConcentricThermoelasticLaw` derives
  layer force, Helmholtz storage, entropy and thermal curvature from the same
  source-bound axial potential. Each layer retains its own temperature, card,
  original annular geometry and integrated thermal-expansion path. Scalar
  reference density and expansion, plus scalar or globally affine `E(T)`, are
  admitted; every original modulus knot is checked. Continuous material support
  and all original property selections survive reevaluation, including auxiliary
  electrical/thermal properties. For elastic strain `e`, the per-reference-length
  potential is `psi=E*A*e^2/2`, entropy is
  `s=A*(E*alpha*e-E'*e^2/2)` and curvature is
  `psi_TT=A*(E*alpha^2-2*E'*alpha*e)`. Consumers add `psi+T*s` to the original
  stress-free caloric energy and must admit positive combined heat capacity,
  `C_caloric-T*psi_TT`, on their actual displacement and temperature domain.
  This shared law supplies constitutive derivatives; continuing motion, heat
  transfer and atomic publication belong to the native coupling owners.
  It does not infer zero expansion, average materials, evolve mass, model slip,
  or certify a radial/finite-strain/phase-changing constitutive law. The focused
  G0/G1 tests cover unequal layer temperatures, opposing stresses, source pins,
  negative expansion and independent free-energy derivative identities. Those
  tests were added but not executed natively in the connector-only session;
  syntax parsing and independent derivative arithmetic were checked.
- `section::winding::CloseWoundSection` retains a sourced circular core and a
  distinct round-wire winding. Exact helical length gives winding line mass;
  the close-coiled spring law gives its axial rigidity with core and winding
  constrained to share end displacement while slipping along the span.
  `WoundCoreSection` retains either the core's sourced Kelvin–Voigt bending
  viscosity or its independent per-layer Prony histories, with common-state
  and intersected core/wrap source-band admission. The consumer projects core
  bending arms using total winding-plus-core modal mass; the wrapping material
  never inherits a bending spectrum or becomes an extra bonded bending layer.
  Geometric closest-turn separation rejects interpenetration. Thin wire,
  small helix angle and long wavelengths are explicit approximation domains.
  The consumer propagates these coefficients through the existing string/body
  mechanics, gas drag and physical pressure. G1 compares mass against resolved
  helix quadrature and stiffness against spring strain energy; G3 changes only
  winding material data. This sliding/core-dominated approximation supplies no
  winding friction/loss, torsion, thermal mismatch or resolved contact, and
  does not claim acoustic equivalence to a measured or resolved specimen.

- `section::wound::HelicalWrapping` is the mass-only interim for wraps outside
  the close-coiled admission domain (thick wire/coil ratios, several wraps, as
  on piano bass strings). Each wrap lies on the diameter beneath it; its line
  mass is the wire density times area times the exact helix length per axial
  length. Neighbouring turns must clear one wire diameter at their exact
  closest approach (any helix angle); `close_wound_pitch_m` solves the
  tightest admissible pitch. Wraps resolve at the core's query coordinates.
  They carry no tension and add no EA, EI, loss, friction or wrap prestress:
  the core section owns those. G0 tests check the helix integral, the second
  wrap's radius, brute-force turn clearance and every refusal.
- `section::strength::AxialTensileCapacity` binds a caller-named tensile-limit
  claim [Pa] (yield, proof or ultimate strength; never chosen, factored or
  converted to an allowable here) from every bonded layer. A common axial
  strain gives layer stress `E_i N / EA`; the capacity is
  `min_i(limit_i EA / E_i)` and that layer governs. A layer without the claim
  refuses as missing data, never as unbreakable. Only axial-force stress is
  checked: residual drawing/winding stress, thermal eigenstrain, vibration
  bending stress and fatigue are absent.
- `section::loss_angle::ConcentricLossAngleSection` binds a named dimensionless
  loss tangent `Im(E)/Re(E)` from every bonded layer: the complex bending
  stiffness is `EI + iH`, `H = sum(E_i I_i tan(delta_i))`, with no averaged
  layer loss. A quality factor or decrement is never converted into it. A
  modal consumer uses the equal-energy-per-cycle rate `H k^4/(mu omega)` at
  each actual mode, inside the intersected source band. Missing, negative,
  wrong-quantity, state-varying and band-less claims refuse; explicit zero is
  lossless. Thermoelastic relaxation from thermal data is in `section::ted`.
- `section::ted::ConcentricThermoelasticBending` computes Zener thermoelastic
  bending loss of a homogeneous solid circular member from its one material
  state: `Delta = E alpha^2 T/(rho c_p)` with T from the state's single typed
  absolute-temperature axis, and `visco::circular_ted`'s exact series over the
  m = 1 Neumann thermal modes, `tau_n = rho c_p R^2/(k xi_n^2)` with
  `J1'(xi_n) = 0` and weights `8/(xi_n^2 (xi_n^2 - 1))` (32 modes; omitted
  weight below 2e-6). `xi_n` comes from a periodic-trapezoid `J1'` and is
  checked against fs-math's certified J1; the weights against independent
  projection integrals. Layered, hollow, untyped-temperature and incomplete
  inputs refuse. It is a lower bound on bending loss: intrinsic anelastic
  excess is the consumer's measured addend or a reported unknown.

- `section::relaxation::ConcentricPronySection` and
  `section::torsion::ConcentricPronyTorsionalSection` compose independent
  generalized-Maxwell bending/shear spectra on the existing bonded concentric
  geometry. Both use one `LayerRelaxationProperties` coefficient resolver;
  the bending-specific spelling names the same source-role type. Each layer
  selects its equilibrium modulus and zero or more modulus/time pairs.
  An explicit empty spectrum is elastic; an omitted layer or absent claim
  refuses. Equilibrium must match the actual channel's Young/shear modulus.
  Each arm retains `delta_E * I` for bending or `delta_G * J`, `J=2I`,
  for circular torsion. Time constants and layer identities are never averaged.
  The original torsional Kelvin viscosity remains a separate parallel dashpot;
  it is not counted again in the Prony spectrum. Geometry, mass, equilibrium
  stiffness and all source receipts survive without authority promotion.
  Density and every coefficient must be constant at the resolved temperature.
  Each relaxation time (or equilibrium modulus for an empty spectrum) supplies
  frequency support, intersected with every other active coefficient domain.
  At most 128 arms per selected section channel are admitted. Missing, ambiguous,
  wrong-quantity, incompatible-band, negative-increment, nonpositive-time and
  unrepresentable-stiffness inputs refuse atomically. G1/G3 tests check independent
  complex bending/torsional stiffness, unchanged parallel viscous loss and
  homogeneous subdivision; G0 covers empty layers and invalid source roles.
  This is frozen-state small-strain bonded bending/circular torsion, not helical
  winding, axial creep, interlayer slip, heating or experimental validation.
- `state_point` is the fail-closed bridge from immutable `fs-matdb`
  `MaterialCard`s to executable scalar property bundles. Consumers declare
  property names, dimensions, numerical domains, the complete physical query
  point, and an explicit claim-selection rule. Resolution is atomic, retains
  each `PropertyUsageReceipt`, canonicalizes requirement order, and binds the
  card, state point, selected values, consumer domains, and receipts into one
  identity. It never infers properties from chemistry names or substitutes a
  representative material.
  `ScalarPropertyRequirement::try_with_quantity` carries an exact `QuantitySpec`
  through both bulk and interface queries, including caller-pinned selection.
  The answer retains the same quantity kind and value form. Dimension-only
  requirements remain explicit and cannot consume semantic claims. Bundles
  containing semantic requirements use identity domain v2 and bind every
  requirement's canonical quantity token; all-dimensional bundles retain the
  exact v1 preimage. Named convenience resolvers still declare their existing
  dimensional property requirements; axis descriptors travel exactly as
  declared by the caller, without inference from names.
  `ScalarPropertyRequirement::try_with_key` also retains the existing matdb
  hardness test context. Both material and interface resolution forward it
  under all three selection policies; a claim pin cannot fill in an omitted
  apparatus/protocol/specimen requirement. The resolved requirement exposes
  that context, and selected claim/usage-receipt identities bind it into the
  existing bundle identity. The focused synthetic G0 regression proves this
  adapter behavior, not measured hardness accuracy or standards compliance.
  Distinct complete keys may share a name. Duplicate keys refuse, and same-name
  pins must identify distinct claims matching the required contexts. Bundles
  order same-name answers by selected claim identity; reversed requirement/pin
  order cannot alter the result. `property(name)` returns `None` for an ambiguous
  name; `property_by_key` requires the complete quantity/test/tensor context.
  The elastic resolver uses that exact lookup for all 36 matrix coefficients.
  `ResolvedScalarProperty::is_constant_at_fixed_temperature` retains the exact
  selected claim's independent axes and admits either validity-wide scalar constants
  or T-only curves evaluated at positive absolute temperature. Frequency curves,
  multiaxis surfaces, other state axes and temperature-difference coordinates refuse. Legacy T
  remains kelvin; typed coordinates must be dimensional temperature or static
  absolute temperature. The original interpolation decision, evidence and
  identities remain unchanged. This classifies a declared source model;
  consumers must still check every operating frequency against its validity.
  G0 tests cover axis semantics, identical values from differently pinned
  curves, and unchanged source receipts. It does not establish measured accuracy
  or admit a changing-temperature constitutive history.
  `ResolvedScalarProperty::angular_frequency_band_rad_s` reads one unambiguous
  frequency coordinate from the exact query and converts its selected-claim
  applicability to rad/s. Static cyclic frequency is multiplied by 2*pi;
  static angular frequency and legacy `omega` retain their values. Semantic
  descriptors take precedence over axis names. Multiple frequency coordinates,
  incompatible `omega` semantics, negative or overflowing bounds refuse;
  inverse-time dimensions alone do not identify an unrelated axis as frequency.
  Absent source bounds remain unconstrained. This nominal unit conversion does
  not modify receipts, extend source support or freeze a sampled frequency curve.
  G0 tests cover conversions, source preservation, ambiguity and invalid bounds.
- `ResolvedMaterialStatePoint::query` reconstructs its complete physical query;
  `at_query_point` refreshes every retained property from the exact original
  card and selected claims. Temperature, pressure, moisture, shear rate and
  other declared coordinate values can change together. All coordinate names
  and quantity conventions, complete property keys (including same-name tensor
  and test contexts), and consumer admissibility remain fixed. A dependent
  property leaving source support refuses the whole candidate instead of
  switching sources or leaving a stale value. `ResolvedLiquidState` uses this
  same refresh and recalculates all derived transport coefficients with its
  original dynamic/kinematic viscosity basis; section and liquid temperature
  laws delegate to it. Focused G0 regressions cover simultaneous multiaxis
  updates, complete-key preservation, typed-coordinate refusals, source pins
  and atomic failure. These are endpoint evaluations, not history integration,
  support along an intervening path, or measured multiaxis validation.
  `ResolvedInterfaceStatePoint::{query, at_query_point}` uses the same refresh
  implementation for ordered surface systems. Every property, complete key,
  selected claim, coordinate convention and consumer domain is retained; any
  unsupported property refuses the whole candidate. A different surface card,
  preparation/history identity or axis schema cannot be substituted by refresh.
  `ResolvedElectricalInterface` and `ResolvedLiquidWettingState` can adopt an
  already resolved multiphysics interface bundle without querying again. Their
  electrical-resistance and equilibrium-angle admission remains active even
  when the generic bundle requested only finite values. Both adapters and
  `ResolvedLinearInterface` rebuild their derived response through this shared
  endpoint refresh, retaining all other interface properties and original source
  receipts. G0/G3 `tests/interface_refresh.rs` exercises simultaneous temperature
  and pressure changes across electrical, mechanical and wetting projections,
  physical coefficient reductions, deterministic refresh, source pinning and
  atomic refusal. This supplies no wear/hysteresis integration, contact geometry
  update, path-support proof or empirical validation; the caller still owns
  evolution and the accepted-state transaction.
  The shared `Surface2D` payload evaluates both source coordinates together
  through the original matdb evaluator. `ResolvedScalarProperty::source_axes`
  retains every declared dependence even at an exact vertex; a surface cannot
  masquerade as a constant or temperature-only coefficient. Sparse source holes
  remain refusals during complete-bundle refresh. Fixed-context temperature
  paths use the canonical `ClaimSet::query_line_typed` seam to retain every
  crossed source-cell breakpoint and prove continuous support. Caloric,
  thermal-expansion and Seebeck consumers integrate those same source slices
  with their existing physics owners. Endpoint support alone never authorizes
  integration across a source hole or nonzero discrete-only span.
- `conductor::resolve_uniform_conductor` binds exactly one sourced
  `electrical_resistivity` claim in ohm-metre dimensions to caller-owned
  positive uniform length and cross-sectional area, retaining the resolved
  material bundle and receipt with `R = rho_e L / A`. It refuses invalid
  geometry, source/unit/domain/selection failures, and unrepresentable derived
  resistance. `joule_power_w` evaluates finite signed-current `I^2 R` only;
  it is not a circuit, contact-resistance, current-distribution, temperature,
  or energy-balance model. The adapter is synchronous, fixed-cost arithmetic
  and has no cancellation claim.
  `OhmicDrive` distinguishes ideal current and voltage control;
  `current_a_for_drive` resolves `I` or `V/R`, refusing nonfinite inputs and
  unrepresentable nonzero currents. The same `joule_power_w` law then supplies
  nonnegative source dissipation. No dynamic circuit or source impedance is inferred.
  `ElectricalPropertyBasis` also admits explicit scalar DC conductivity through
  `resolve_conductor_state_with_basis` and `resolve_uniform_conductor_with_basis`.
  The source is evaluated at the query point before taking its reciprocal;
  reciprocal knot interpolation is not substituted. The original property,
  uncertainty and receipts remain intact under `source_property()`; derived
  resistivity is never inserted as a measured claim. Existing resistivity query
  entrypoints keep their requirement unchanged. `try_from_material` can adopt a
  complete multiphysics bundle with either one source basis; two bases require
  an explicit independent-variable choice and are not averaged. Consequently,
  concentric sections can mix conductivity- and resistivity-sourced layers in
  their existing parallel-current and layer-specific Joule-heat calculation.
  This scalar reciprocal law does not invert tensors or admit AC response.
  `resolve_conductor_state_with_property` and
  `resolve_uniform_conductor_with_property` bind a caller-selected exact
  `PropertyKey` to that same explicit resistivity/conductivity basis. They do
  not normalize names or infer aliases. The key must retain the basis's exact
  dimensional quantity; the complete original bundle, conditions and receipts
  remain authoritative. `source_key()` survives bulk-to-uniform lowering and
  material refresh. Planar source-name iterators borrow the actual bound names,
  and coupled refresh guards retain the selected key even when two properties
  in the same bundle have equal values. The canonical entrypoints retain their
  original requirements and refusal behavior.
- `field_storage::{ResolvedCapacitor, ResolvedInductor}` adopt complete material
  state bundles and derive `C = epsilon A / d` and `L = mu N^2 A / l` from
  explicit uniform geometry. Absolute SI electric permittivity [F/m] and
  magnetic permeability [H/m] are required; relative coefficients cannot be
  silently substituted. Original source claims and state receipts are retained.
  Charge and flux-linkage energies use `q^2/(2C)` and `lambda^2/(2L)`.
  Only linear isotropic, nondispersive, lossless uniform-field reductions are
  admitted. Curves and surfaces may be resolved at fixed temperature and typed
  pressure, composition, moisture or concentration. The shared fixed-external
  source policy checks every value axis and applicability bound; dimensional
  equality does not turn stress into pressure or strain into composition.
  Frequency-, field-, strain- or stress-dependent/restricted sources require
  another storage model. Changing a resolved coefficient is a new snapshot;
  live evolution must account its parameter work and any field entropy.
  Invalid geometry, ambiguous sources and unrepresentable coefficients/energy
  refuse. These reductions introduce no circuit solver, hysteresis or loss law.
- `thermoelectric::ResolvedSeebeckSpan` resolves one source-owned scalar, linear-inside curve
  or fixed-context two-coordinate surface slice of
  `absolute_seebeck_coefficient` [V/K] at both ends of a fixed-context temperature
  span and at its actual bulk temperature, preserving original receipts and the
  same selected claim. Joule-heated bulk states cannot escape source support
  merely because the junctions remain cooler. Its
  zero-magnetic-field convention is `J = sigma(-grad(V)-S grad(T))`,
  `q = -k grad(T)+S T J`. Exact piecewise-linear integration supplies `S_left`
  and `S_right`, the integral means on a-to-bulk and bulk-to-b. The voltage drop
  is `integral(S dT)`; heat into the three stores is `[-I S_left T_a,
  I (S_left-S_right) T_bulk, I S_right T_b]`. These cross terms exchange equal
  and opposite electrical/thermal work and zero net entropy analytically.
  The bulk term is an interior-node Thomson allocation. Endpoint edge averages
  contain the remaining spatial allocation, so these are discrete nodal heat
  transfers, not exact continuum endpoint Peltier fluxes or a whole-leg Thomson
  source. Spatial refinement is required to resolve that distribution. A
  constant scalar retains the original zero-bulk-heat Kelvin law. Isothermal
  edges use the actual sourced S; signed integral cancellation is permitted.
  Source-cell breakpoints are retained before integration; curves on a fixed
  non-temperature axis contribute their evaluated constant coefficient.
  Relative thermopower, nonzero tabulated-only spans, source holes and
  extrapolation refuse. The entire temperature path, including the actual
  bulk state, must remain supported by the SAME selected claim.
  Original conductor, Fourier and caloric owners supply all
  dissipative transport and storage. This module supplies no independent circuit
  or thermal integrator and no tensor, magnetic-field or phase-change model.
- `liquid::resolve_liquid_state` resolves positive density, specific heat,
  conductivity and dynamic viscosity at one caller-supplied source condition.
  It retains the complete four-property bundle and derives kinematic viscosity,
  thermal diffusivity and Prandtl number, refusing non-finite or non-positive
  results. Liquid phase, formulation, concentration and operating support come
  from the selected card and query; the adapter adds no EOS, concentration
  conversion, pressure-independence or phase-stability evidence. Arithmetic
  after bounded resolution is synchronous and fixed-cost. Source-driven heat
  and flow comparisons remain numerical checks of the declared properties,
  not experimental qualification. The separate
  `liquid::resolve_liquid_state_from_kinematic` entry point requires a
  kinematic-viscosity claim instead and derives dynamic viscosity as `rho*nu`
  at the same query point. It retains the original kinematic receipt and
  refuses unrepresentable conversion; neither API silently substitutes the
  other viscosity basis or certifies consistency of independent source tables.
- `resolve_isotropic_elastic_state_point` requests only density, Young's
  modulus, and Poisson ratio for linear vibration/elasticity. The distinct
  `resolve_isotropic_solid_state_point` additionally requires yield stress for
  contact/plastic admission; consumers do not need to fabricate an unrelated
  yield datum merely to compute modes.
  `resolve_isotropic_thermoelastic_state_point` atomically resolves density,
  Young's modulus, Poisson ratio, linear expansion, specific heat, and thermal
  conductivity at an explicit positive `T` in kelvin. It retains all six
  property-use receipts and their original uncertainty and supplies the existing
  Zener law. Missing properties and out-of-domain states refuse; no density or
  chemistry-name preset is selected. This is the homogeneous isotropic,
  single thermal relaxation mode approximation. Thermal conductivity may be
  explicitly dimensional (legacy data) or the exact static ThermalConductivity
  schema; a semantic claim of any other kind/form refuses. The NASA TN D-6448
  seed's v3 manifest declares its source-to-runtime names, thermal kind and
  absolute-temperature axis. Typed T coordinates must be absolute temperature;
  a temperature difference cannot supply the Zener law's reference temperature.
  The compiler-to-state integration test uses that actual source at 300 K,
  retaining its Unstated uncertainty and point-only temperature validity.
  Resolved material/interface bundles expose `axis_quantities()` from their
  member receipts. Existing bundle hashes transitively bind every descriptor
  through these receipt identities without duplicating the quantity system.
  This is not an anisotropic law or a
  resolved thermomechanical energy balance. Its synthetic temperature curves
  and independent Zener checks establish numerical behavior, not experimental
  validation.
- `resolve_isotropic_thermoelastic_state_point_with_requirements` resolves the
  same mandatory six-property bundle together with explicitly requested auxiliary
  properties. One card, state point, exact pin plan and resolved identity cover
  the complete bundle; auxiliary values cannot replace or ambiguate required
  thermal/elastic names. The existing 64-property cap is checked before cloning
  auxiliary requirements. The original convenience resolver uses an empty
  auxiliary set and retains its prior requirements and identities.
  The material-driven plate compiler consumes this API to bind thermal loss and
  independently selected bending viscosity to the same physical state. G0/G1/G3
  integration tests in fs-couple cover required-name shadowing, count limits,
  exact claim pins in the presence of conflicting heat-capacity claims, typed
  conductivity, temperature-dependent input curves and actual emitted pressure.
- `resolve_isotropic_thermal_expansion_state_point` independently resolves the
  instantaneous isotropic linear expansion coefficient, including its exact
  state-point validity and property-use receipt. Negative expansion is
  admissible; material names imply neither its sign nor its value. This result
  is a coefficient, not a total strain: a coupled driver must integrate it over
  the admitted temperature path before solid mechanics use.
  `integrate_isotropic_thermal_expansion` performs that integration exactly for
  scalar, piecewise-linear curve and fixed-context surface claims. It requires
  absolute temperature (including legacy kelvin `T`), keeps every other state
  axis fixed, proves continuous support through the canonical source line,
  and pins the complete path to one selected claim. Crossed source-cell
  breakpoints preserve changes in interpolation slope. The complete path and
  signed `integral(alpha(T) dT)` result are identity-bound.
  `integrate_orthotropic_thermal_expansion` uses that same integration owner
  for the three explicit `linear_thermal_expansion_coefficient_1/2/3`
  properties in inverse-kelvin dimensions. Each principal coefficient retains
  its own reference claim through the complete fixed-context path. The result
  retains all three source records and signed principal strains; negative
  coefficients, cooling and fixed-context surface slices use the same rules.
  Missing directions, ambiguous claims, wrong dimensions and support gaps
  refuse the whole result. These are additive small strains in one fixed
  principal frame. The consumer must supply its proper rotation into the mesh;
  neither a general six-component coefficient law nor changing principal axes
  are inferred. The focused regression source is
  `tests/orthotropic_thermal_expansion.rs`.
  `resolve_mean_thermal_expansion_strain` separately resolves source-declared
  interval means as engineering strain `mean_alpha * (T_end - T_reference)`.
  Both positive absolute-temperature coordinates must be typed and constrained
  by the selected source validity. The derived result retains the parent
  evidence and binds its reference/endpoint roles in its identity; it refuses
  nonfinite results and nonpositive length ratios. Callers must explicitly
  select a mean coefficient: equal dimensions do not establish that convention.
  This endpoint consequence supplies no continuous alpha law, density law,
  thermal stress, or propagated statistical uncertainty. The stored steel
  source integration in `xtask/tests/matdb_pack_cli.rs` covers the five observed
  intervals, receipt preservation, replay and unsupported-state refusals.
  `resolve_orthotropic_elastic_state_point`
  atomically resolves density plus the nine principal-axis engineering
  constants and refuses a non-positive-definite compliance. Spatial material
  orientation remains an independently identity-bound field owned by the
  geometry/solid consumer; no chemistry name implies isotropy or a grain axis.
  `resolve_orthotropic_elastic_state_point_from_reference_ratios` derives E/G
  as the exact binary64 products of receipted E_L and the source's ratios.
  L/R/T map to axes1/2/3, with shears in G_LR,G_RT,G_TL order. The caller
  declares Full3d or an ordered material plane and the MajorRatios Poisson
  rule. Full3d returns the existing complete law; a plane result has only its
  two Young moduli, major Poisson value and shear modulus, never fabricated
  out-of-plane constants. Missing required ratios/density refuse by name.
  Independently measured minor Poisson ratios for requested pairs, when
  present, retain receipts and signed relative reciprocity residuals; they
  are never averaged into the major set. Unrepresentable residuals refuse.
  `ReferenceRatioDerivation` uses domain
  `org.frankensim.fs-material.orthotropic-reference-ratio.v1` to bind source
  bundle, set, axes, rule, strain limit and residuals. Canonical-key states
  retain their existing bundle identity and have no ratio derivation.
  Derived products are not inserted as fake source claims. FPL E_L remains
  the quasi-static reference modulus with its original moisture validity and
  Unstated uncertainty; no dynamic correction, moisture law, spatial grain
  orientation or propagated covariance is inferred. This is bounded exact
  algebra after ordinary source resolution, with no iterative cancellation
  claim. G0 tests cover compliance roundtrip2e-15 (8ulp), product/receipt bits,
  reciprocity, plane support and refusals; the source-to-consumer lane checks
  all23 FPL packs (12 Full3d,18 Plane12) and bitwise plate section/mode controls.
  Native qualification is recorded separately; these tests are not specimen
  validation. No material pack schema or source row is changed.
  `resolve_elastic_tensor_state_point` resolves density and all 36 explicitly
  addressed coefficients from one card at one supported point. The requested
  keys must share a source tensor identity, basis/frame, notation, order and
  symmetry, and address each cell exactly; missing entries and implicit zeros
  refuse. The existing scalar requirements carry the complete descriptor through
  all bulk/interface selection policies and retain each usage receipt.
  `ElasticTensorStatePoint` exposes the nominal source matrix and its evidence
  bundle. Full-law symmetry/stability and spatial rotation are admitted by
  `fs-solid::linear3::TetElasticMaterial::from_resolved_elastic_tensor`.
  Marginal uncertainties remain upstream; this adapter does not assume
  independence or produce a certified rotated uncertainty bound.
  `resolve_strain_tensor_state_point` resolves all six explicit symmetric strain
  coordinates through the same scalar requirement/selection path. Exact keys
  must share notation, order, frame and source tensor identity and address each
  position once. Negative and zero strain components are admitted; missing
  components never imply zero. `StrainTensorStatePoint` retains the six nominal
  values and all selected evidence at the query point. The solid thermal adapter
  can consume this bundle with a caller-declared integrated stress-free strain
  interpretation. Neither point support nor a source tensor name proves that
  temperature integration was performed or that a complete path is supported.
  `resolve_stress_tensor_state_point` likewise resolves six explicit symmetric
  Cauchy stress coordinates from one material card and supported query point.
  Complete scalar requirements retain the stress descriptor through all bulk
  and interface selection policies; equal pressure dimensions do not alias
  stiffness. `StressTensorStatePoint` exposes pascal values in the declared
  notation/order/frame, source tensor identity and the selected evidence bundle.
  Missing entries, mixed tensors/bases and context-free requirements refuse.
  Signed stress is admitted; engineering stress shear is physical sigma_ij.
  Source compiler v8 and normalized-pack v7 feed this path. The synthetic
  Pa/kPa import test compares it with the existing rotated isotropic law's
  constrained thermal response. It adds no solved field, initial-stress load,
  inferred covariance, finite-deformation stress or experimental validation.
  `resolve_joint_stress_tensor_state_point` resolves the same six nominal
  stress coordinates and their joint receipt from one immutable material-card
  pack, point and selection policy. Each nominal member receipt must match the
  joint answer. Private `JointStressTensorStatePoint` fields bind both under
  `joint-stress-state-point.v1`; a caller cannot substitute unrelated covariance.
  The source-coordinate covariance or exact unknown-correlation reason is
  retained. Nominal values are not asserted to be distribution means, marginal
  widths are not converted to covariance, and no orientation uncertainty is
  inferred. The existing two joint-query policies are supported. The fixed-frame
  numerical transform is owned by fs-solid, with no new pack schema or dependency.
  `resolve_joint_strain_tensor_state_point` binds that nominal state to the
  existing joint query on the same `NormalizedMaterialCardPack`, component
  requirements, physical point and selection policy. It compares all member
  receipts and retains the source joint receipt, including covariance in source
  component order or the exact unknown-correlation reason. Its private
  `JointStrainTensorStatePoint` fields prevent injection of an unrelated matrix.
  Domain `joint-strain-state-point.v1` binds nominal state and joint receipt.
  The existing single-claim and observation-backed joint policies are supported;
  this adapter adds no pinned joint policy or marginal-width-to-variance inference.
  The two-material reduced modulus on the contact-capable isotropic rung uses
  the standard ordered Hertz half-space compliance and binds both resolved
  state identities. Out-of-range temperature or phase-dependent data refuse
  through the material query; this module does not silently extend a solid
  card through melting or select a new phase/model.
- `phase::EquilibriumEnthalpyPhaseCurve` is the generic phase-state ingress:
  source-provided monotone specific-enthalpy knots resolve temperature, bulk
  density, and solid/liquid mass fractions, including an isothermal
  latent-heat plateau. It contains no material-name presets and never
  extrapolates. `try_single_phase` admits only explicitly all-solid or
  all-liquid curves with strictly increasing temperature; it refuses a
  latent plateau or mixed fractions disguised as a one-phase curve. It is not
  a heat-transfer, deformation, remeshing, free-surface, acoustic, or optical
  solver.
  `temperature_derivative_at_specific_enthalpy` returns the admitted chart's
  piecewise-linear dT/dh, exactly zero on a latent plateau. Interior knots use
  the right segment; the upper endpoint uses the left. This is an explicit
  generalized derivative for semismooth transport, not a smoothed heat capacity
  or a classical derivative at a corner. Nonfinite/out-of-domain inputs and
  unrepresentable positive slopes refuse. Signed-zero enthalpy endpoints resolve
  identically in the state and derivative queries. The focused public derivative
  tests include independent segment differences, latent plateaus, knot sides,
  single-phase charts, domain errors and representability boundaries.
  `try_from_heat_capacity` derives that same one-phase chart from positive,
  strictly increasing source `Cp(T)` knots, taking its enthalpy reference at
  the first source temperature. Linear `Cp` segments are integrated exactly
  and subdivided under a caller-provided real-arithmetic chord-error bound;
  all source endpoints remain knots. The bound excludes source uncertainty,
  floating-point error and the linear-in-temperature density approximation.
  Generated charts cap at 4096 knots before allocation. This is bounded
  numerical construction, not experimental validation or a mechanical model.
  `ResolvedHeatCapacityCurve::try_new` admits that chart directly from one
  immutable material card using explicit Cp/density keys and complete typed
  endpoint queries. Both properties require continuous support over the same
  fixed-context temperature span through `ClaimSet::query_line_typed`.
  The union of their curve knots and crossed surface-cell breakpoints is
  retained (at most 4096 distinct temperatures), including density-only knots. Ordered
  Cp/density receipts remain accessible and are bound into the curve identity.
  Wrong dimensions, changing non-temperature conditions, unsupported spans and
  discrete-only data refuse. The reference enthalpy and single-phase choice
  remain caller assumptions; this ingress does not infer phase stability,
  phase-transition endpoints, source accuracy or material compatibility.
  `UniformEnthalpyStepInput` and `UniformEnthalpyStep<R>` carry a proposed
  uniform-body heat-transfer step between physical owners without a reverse
  solver dependency. They include specimen mass/geometry, initial chart/state,
  duration, internal heat, actual boundary heat, a distinct solve-residual
  budget and the concrete solver report. They are data, not an admission or
  transport solver. The transport must justify lumping; the coupled owner must
  verify chart identity, energy balance and regime before publishing, and the
  callback must not debit an external source before that commit.
- `resolve_interface_state_point` applies the same atomic query contract to an
  ordered `InterfaceSystemCard`. Its identity therefore includes both surface
  material states and texture frames, medium, third body, environment,
  preparation/history, complete query point, claim-selection policy, selected
  values, and property-use receipts. Interface properties never migrate onto
  either bulk material merely because the two chemistry names match.
- `ResolvedLinearInterface` lowers explicit stiffness/viscosity-per-area claims
  from that same resolved interface into passive bilateral `k` [N/m] and `c`
  [N s/m] over a declared positive area. It retains the complete ordered source,
  selected property names and area; exact receiving material states and actual
  angular frequencies can be checked without replacing receipts. Both claims
  require explicit intersected frequency support and fixed-temperature constant
  coefficients. Missing loss, frequency-dependent curves, two zero coefficients,
  negative values, overflow and disappearing nonzero products refuse. Geometry,
  work-conjugate force maps, rest displacement and stepping remain the consumer's
  existing owners; this is shared material lowering, not another joint solver.

- `SmallStrainLaw` trait (Voigt 6-space, TENSOR shear components):
  `stress` / `tangent` / `update_state` / `admissibility` / `card`.
  **The tangent contract**: `tangent` is the exact derivative of the
  ALGORITHMIC stress update at the same committed state — FD-gated for
  every law in conformance (merge-gate discipline, same as adjoints).
- `IsotropicElastic` (E, ν → Lamé), `OrthotropicElastic` (engineering
  constants; construction REFUSES thermodynamically inadmissible Poisson
  sets via compliance positive-definiteness minors).
- `J2Plasticity`: radial-return mapping with linear isotropic or explicitly
  split isotropic/Prager kinematic hardening
  and the Simo–Hughes algorithmic moduli
  `C = κ I⊗I + 2μθ I_dev − 2μθ̄ n̂⊗n̂`.
- `Hyperelastic` (`NeoHookean`, `MooneyRivlin`, decoupled `Ogden`,
  `HolzapfelGasserOgden`): stored energies written
  ONCE generic over the fs-ad `Real` scalar; `piola` is the exact dual
  gradient, `tangent` the exact nested-dual Hessian (9×9). `det F ≤ 0`
  and nonfinite states refuse structurally. Ogden uses the original
  `mu/alpha` convention: infinitesimal shear `mu = sum(mu_p alpha_p)/2`.
  At most eight finite terms with `mu_p alpha_p > 0` are admitted, including
  negative-mu/negative-alpha terms. The shared principal-stretch helper uses
  the existing symmetric Jacobi owner; a spectral trace lift carries the
  exact value/gradient/Hessian through repeated roots, using the analytic
  divided-difference limit without perturbing stretches or differentiating
  eigenvectors. Derivatives above second order are not claimed.
  HGO admits at most four explicitly unit reference directions, `c10 > 0`,
  `k1 >= 0`, `k2 > 0`, dispersion in `[0,1/3]`, and an explicit full or
  isochoric structural invariant. Its energy per family is
  `k1/(2 k2) (exp(k2 <E>_+²)-1)`. At zero dispersion `I4 <= 1` contributes
  exactly zero; with dispersion, tension refers to the distributed
  structural strain `E`, preserving C1 continuity and direction independence
  at dispersion 1/3. The inactive-side Hessian is returned at E=0; a
  two-sided Hessian there is not claimed. A single aligned family supplies
  `2 k1 (lambda²-1) exp(k2 (lambda²-1)²) lambda²` incompressible Cauchy
  stress; two identical families supply twice this value.
  All four variants enforce all principal stretches in `[1/limit,limit]` in
  stress and tangent evaluation, with eight-EPSILON eigensolver boundary
  roundoff; limit*(1+1e-12) refuses. The limit must be finite and at least one,
  including after mutation of the public field. Generic-scalar energy
  evaluation retains its documented admissible-F precondition. The Ogden and
  HGO variants offer quadratic or log-squared
  bulk energies and separate distortional/bulk accessors. Full-invariant
  HGO's fibre accessor intentionally retains volume dependence.
- `hyper::cards` instantiates strict, dimensioned `ConstitutiveModelCard`s
  for `neo-hookean-log`, `mooney-rivlin-decoupled`, `ogden-decoupled` and
  `hgo-dispersed` v1 and registers the real
  stress/tangent/energy factories in `graph::LawRegistry`. Graph inputs are
  H=F-I, outputs P, so the zero-input admission probe is the reference state.
  Moduli require Pa, exponents/directions/counts/selectors are dimensionless;
  missing/extra parameters, fractional counts, versions and undeclared
  invariant selectors refuse. FSMODPK v1 remains unchanged.
  Neo-Hookean supplies `mu` and Lamé `lambda` in Pa for the native
  `mu/2 (I1-3) - mu ln(J) + lambda/2 ln(J)^2` energy. Mooney–Rivlin
  supplies `c10`, `c01` and `kappa` in Pa for its native isochoric-invariant
  energy plus `kappa/2 (J-1)^2`. Both require dimensionless `stretch_limit`;
  neither accepts an unused bulk selector. Native admission requires finite
  coefficients and representable reference shear/longitudinal stiffness.
  All four laws traverse the existing model-pack decode and registry test;
  the isotropic laws and Ogden declare smooth nodes, HGO piecewise smooth.
  `hyper::cards::ResolvedHyperelastic` retains the exact source model card and
  admitted `QueryPoint` beside its native law. `try_from_model_card_at` requires
  the typed source domain to contain that fixed query; bounded conditions
  support static temperature, pressure, composition, moisture and concentration
  quantities, plus legacy kelvin `T`. The context-free constructor admits only
  unconstrained domains. Evolving kinematic calibration axes refuse until an
  actual coordinate owner exists. `HyperelasticMaterial` lets spatial consumers
  borrow this same energy/Piola/tangent owner. Native parameter and actual
  principal-stretch checks remain unchanged; the wrapper adds no thermal
  evolution or uncertainty propagation.
  The offline `xtask matdb-pack` `constitutive-v1` profile accepts explicit
  model/law/state/initialization/validity declarations and dimensioned
  parameter/validity rows, normalizes through fs-qty, retains each transform
  receipt, and self-decodes the existing runtime model pack. Its compiler
  tests transport the Ogden/HGO laws into the real registry. A generic profile does
  not certify any law's physical parameters or infer a tissue identity.
  Polyconvexity is parameter-derived: Ogden sufficient exponents are
  alpha>=3/2 (F power perspective) or alpha<=-3 (cof F perspective), with
  positive coefficients and convex bulk. Full-invariant HGO has convex
  positive-part structural energy. Other sets return `None`; no numerical
  spot check upgrades an unknown result to a proof. Log-squared bulk is
  certified only on declared J<=e, where its second derivative is nonnegative.
- `active_tissue::ActiveTissue` composes a passive `Hyperelastic` with a
  prescribed `ActiveFibreStress` along a reference fibre `a0` (unit; for an
  HGO passive card it must be one of its fibre families). The active Cauchy
  stress is `T_max a f_L(lambda_f) m (x) m` with `lambda_f = |F a0|` and the
  activation `a` in `[0, 1]`; its first Piola image is
  `(J / lambda_f^2) T_max a f_L (F a0) (x) a0`. The factor `f_L` is a Gaussian,
  a parabola or a tabulated piecewise-linear family (at most 64 knots). The
  law refuses a fibre stretch outside its declared inclusive domain (eight ULPs
  of roundoff at the bound) and an activation outside `[0, 1]`
  (`MaterialError::State`), and a non-positive or non-finite `T_max`, a non-unit
  fibre or an invalid family (`MaterialError::Parameters`). Zero activation
  returns the passive stress and tangent bit for bit and does not consult the
  domain. The active stress is a PRESCRIBED, non-conservative input: the tangent
  is the passive (symmetric) tangent plus the exact forward-mode derivative of
  the active stress, which is not symmetric (`tangent_symmetric: false`), and
  its work `P_a : Fdot` is booked as an external power input, never as stored
  energy. Activation is a performance control that no law input, card or pack
  sets here.
- `Uniaxial` trait + the RC flagship pair: `MenegottoPintoSteel`
  (R0/a1/a2 curvature degradation, Bauschinger via branch-state
  asymptote intersection) and `ManderConcrete` (confined envelope
  `f = f′cc·x·r/(r−1+xʳ)`, elastic unload/reload lines with residual
  strain, zero tension).
- `calibrate_bilinear`: segmented least squares recovering (E, σ_y, H)
  from monotonic data with standard-error envelopes and RMS residual.
- `evidence_stress`: wraps any law's stress in `Evidence` with the card
  attached and `in_domain` FLAGGING calibration-domain exit.
- `tensor`: Voigt helpers (deviator, contraction with shear doubling,
  von Mises, Rodrigues rotations) used by the objectivity gates.
  Shared `StressTensorBasis` / `StressTensorNotation` descriptors distinguish
  physical/engineering stress shear (both `sigma_ij`) from Mandel stress shear.
  They are re-exported from fs-matdb for fs-solid's recovered thermal-stress
  outputs; they do not reuse engineering-strain factor two or imply Piola stress
  or material/frame calibration. Source-pack stress transport uses the separate
  StressTensorComponent context described above.

### `gas` — first-principles ambient gas state (zdmm1 slice 1)

`thermochemical_gas::packs::from_packs` joins explicit verified species,
NASA-9 and pure-transport material-card packs. Names supply no properties.
The model and each resolved state retain the three pack hashes per canonical
component, alongside existing thermodynamic/transport receipts. Missing NASA-9
or transport data refuse. `gas_species::input::read_gas_record` requires an
already-loaded registry for thermochemical input and returns the full resolved
state; the legacy two-argument reader preserves its old bits and refuses that
record without a registry. A named `mixture` record keeps its species and
mole fractions beside its state (`GasRecord::Mixture(state, composition)`;
a mass record is converted by the species' own molar masses), so a consumer
that needs the composition reads it and never infers it from the state.
Phase admission is explicit: frozen gas makes no
condensation claim; the liquid-water option applies the existing Buck limit
only inside its stated temperature window. No real-gas or molecular-relaxation
accuracy is supplied by this adapter. Literal NASA/JANAF disagreements remain
test failures after accounting separately for each printed half-unit. Relative
bands are unchanged; cp also satisfies its relative-only gate. Rows absent from
the source are recorded as absent, never interpolated into alleged facts.

Frozen gas parcels and nozzles use the thermochemical owner's direct entropy
differences instead of subtracting absolute reference entropies. The relative
density/entropy closure targets a requested change from the retained initial
state and preserves the original bracket, cancellation and source checks.
Reversible volume changes target zero entropy change. Midpoint advancement
reports its actual numerical entropy change without asserting exact adiabatic
reversibility. Iterative entropy and nozzle trial states use allocation-free
caloric properties; accepted states retain complete source receipts.

The converging `nozzle_flow_to` closure accepts vacuum back pressure without
querying an EOS at zero pressure. Sonic and kinetic-energy residuals use the
thermochemical owner's direct enthalpy differences; common formation offsets
cannot erase available nozzle energy.

`nozzle::rocket` adds ideal converging-diverging propulsion using the same
frozen NASA-9 composition. A bounded sonic solve determines mass flow from
chamber state and throat area, then a bounded supersonic temperature solve
closes exit mass flux using source enthalpy and entropy differences. Trial
caloric evaluations allocate no material receipts; accepted throat and exit
states retain full source and transport resolution. Explicit mass, energy and
entropy residuals are checked before publication. Outputs include exit Mach,
momentum and pressure thrust, total thrust, effective exhaust velocity,
characteristic velocity and specific impulse. The latter uses standard
gravity only as a unit convention. The off-design extension below admits
source-supported unchoked discharge, ideal internal normal shocks and attached
external compression in addition to matched/underexpanded and vacuum flow.
It does not infer viscous separation, combustion equilibrium, real-gas
validity, nozzle losses, adjoints or resolved multidimensional flow.

`nozzle::blowdown` removes finite mass from a frozen, well-mixed, adiabatic
fixed-volume reservoir on the original density/entropy path: `d(mu)=h dm`
implies `ds=0` for the remaining gas even with temperature-dependent heat
capacity. It exports actual `U0-U1` once, retaining formation conventions and
every species. `rocket_blowdown` closes `dm=dt*mdot(m0-dm/2)` by bounded
implicit midpoint iteration and reports midpoint thrust impulse. Initial,
midpoint and endpoint nozzle/source domains and the donor-fraction limit are
checked transactionally. Smooth time accuracy is second order, as checked
against analytic perfect-gas mass, pressure and impulse refinement; a mass
fraction limit does not certify time error. No wall heat, moving tank,
combustion, resolved chamber/jet dynamics or double-counted exhaust kinetic
energy is supplied by this lane.

Source-backed gas temperature and pressure support comes from the selected
NASA-9 regions and every pure-transport claim, not from the legacy USSA/Sutherland
air constructor's 50..2000 K and 100 bar bounds. `GasState::validate_coherence`
checks positive finite coordinates and the common EOS/caloric/transport/acoustic
identities; it does not assert source coverage. Energy and entropy inversion
accept positive finite temperature brackets and retain their original source
checks, residual tolerances, cancellation and atomic updates. Consequently
source-supported hot gases can heat, compress and supply the same density,
viscosity, conductivity and sound speed to existing consumers without changing
their physics code. Legacy `GasSpec` constructors keep their original limits.
G1/G3 manufactured tests cover heating through 2000 K, variable-cp reversible
compression, data-only composition changes and unchanged state after leaving
transport support. The NASA-9 family covers species-specific high-temperature
regions (NASA/TP-2002-211556, <https://ntrs.nasa.gov/citations/20020085330>);
the active source bounds still govern each calculation. This does not admit
condensed phases, real-gas behavior, equilibrium chemistry or plasma physics,
and source support alone is not an experimental-accuracy claim.

`gas_species::input::read_isochoric_parcel` retains the same parsed
thermochemical model, typed conditions and canonical source associations in
the existing `IsochoricGasParcel` owner. Its mass is rho*V from the declared
volume and initial state; later signed heat evolves the original NASA-9 energy,
EOS and transport rather than a frozen initial cp. `input::packs::read_source_parcel`
uses the same bounded canonical pack loader and complete-record checks as the
snapshot reader, and the parcel owns its sources after caller bytes are dropped.
Borrowed prefix input is unchanged on every refusal and retains a caller-owned
suffix on success; the complete pack reader rejects trailing fields.

Finite-parcel ingress explicitly requires `FrozenGas`: an initial liquid-water
saturation check cannot certify later cooling, so `LiquidWaterLimit` is refused
instead of silently weakened. The caller remains responsible for a valid
frozen-gas phase domain during evolution. No condensation, reaction, spatial
flow or experimental accuracy is supplied by this input adapter. The existing
snapshot grammar, source selection and phase checks are shared without changing
their numerical evaluation. Focused G1/G3/G4 regressions in
`fs-couple/tests/gas_source_packet.rs` check source retention after input drop,
nonlinear caloric heating and pressure in both composition bases, exact snapshot
parity, and transactional source/volume/phase/suffix refusals.

`gas_field::AxialGasField` stores 2..=4096 supplied axial stations with
temperature and canonical mole composition at one explicit common pressure.
It interpolates only primitive coordinates; NASA-9 thermodynamics, EOS, Wilke
viscosity and WMS conductivity are reevaluated by the existing owner at every
sample. Uniform fields and original station values preserve direct-evaluation
bits. No extrapolation, invented transport data or independent interpolation of
rho/c/cp/mu is allowed. `with_composition` retains the active species' original
laws, transport selections and pack identities; explicit zeros omit inactive
species without fabricating trace concentrations. The full source model remains
on the field. Only generated interpolation roundoff receives a bounded simplex
correction; authored fractions remain subject to the existing exact-sum rule.
The existing liquid-water phase policy is explicit; supersaturation refuses
instead of silently removing water. The field is prescribed, not a heat/species
transport or condensation solver, and source-station measurement uncertainty
remains with its caller. Construction and sampling are synchronous bounded
operations; callers can chunk sampling. Focused tests cover uniform bit identity,
fresh primitive evaluation, zero-concentration ingress, source retention and
station/sample/domain/phase refusals. Native verification is recorded separately.

The five-species acoustic fixture binds NASA/TP-2002-211556 Appendix D to
primary zero-density transport correlations: Lemmon/Jacobsen (N2/O2/Ar),
Laesecke/Muzny and Huber (CO2), and IAPWS R12-08/R15-11 (H2O). The stored
linear curves use 256 uniform intervals plus declared query knots; tests
measure 64/128/256 rungs and require relative representation error <= 1e-4.
Source uncertainty remains separate, without invented confidence intervals.
Water's 250 K dilute-vapor continuation has no quantified sub-triple accuracy
certificate. These curves do not certify finite-density or mixture errors.
The SRM1720 fixture uses the pack's information values, not certified complete
air. Cramer1993 remains source-unverified: the 108-state humid grid computes
NASA-9 sound speeds but emits NO-DATA and fails until the primary comparison
and its stated uncertainty can be verified. No secondary transcription is
used to qualify the requested 5e-4 gate.

`relaxation::RelaxationAbsorption` (M0.1-W2) is the molecular-relaxation law
for any composition. Each species declares its relaxing internal reservoirs:
vibrational modes with Planck-Einstein heat capacity, and classical rotation.
They enter one exact complex `c_v(omega)`, and `k = omega/c(omega) + i alpha_SK`.
Stokes-Kirchhoff is counted once, there is no separate bulk viscosity, and the
strengths are not linearized.

`RelaxationState::caloric_response` exposes that same bounded, allocation-free
reservoir sum as complex mass-specific `cv`, `cp = cv + R`, and bulk modulus
`K = p cp/cv`, including the equilibrium DC limit. Under `exp(-i omega t)`,
`Im(K) <= 0`; the cycle-mean compressional work is
`-omega Im(K) |strain|^2 / 2` per volume for peak phasors. Duct wall thermodynamics
and coupled cavity pressure use this constitutive owner directly. It adds no
shear/thermal propagation attenuation, boundary-layer loss or second bulk
viscosity. Single-reservoir real-arithmetic rational response, DC/frozen limits,
the original wavenumber relation and the zero-reservoir lossless limit are
checked in `tests/relaxation_absorption.rs`.

`RelaxationState::bulk_relaxation` realizes that same finite caloric law as
equilibrium bulk storage plus positive Maxwell compression arms. It uses the
existing symmetric Jacobi owner on the diagonal-plus-rank-one caloric matrix,
checks eigenpairs and static/frozen limits, and retains every original source
reservoir. Equal molecular times combine exactly; no fitting, pole truncation,
propagation loss or new gas-specific constitutive law is introduced. The work
is bounded by the original 32-reservoir cap. Single-reservoir analytic
elimination, multi-reservoir complex response, equal-time/zero-strength limits
and pressure/rate scaling are checked in `relaxation::bulk::tests`.

Relaxation times are data: Landau-Teller collision-partner rates with a linear
mole-fraction rule, a partner-independent power law, or ISO 9613-1's `f_rO`
and `f_rN` (`iso9613_relaxation_frequencies_hz`, bit-identical to the ISO
formula's). The ISO laws are a declared terrestrial-air variant: dry-basis O2
and N2 must be within 1 % of 0.209 and 0.781, with CO2 at most 1e-3. Any other
mixture is refused, naming the species (the same admission as
`iso9613::TerrestrialAir`). `RelaxationAbsorption::iso9613_terrestrial_air()`
is that data set: the O2 and N2 vibrations at ISO's theta and rates, with Ar,
CO2 and H2O declared without reservoirs (ISO's own assumption). Rotation is
not a reservoir in it, so the rotational bulk loss inside ISO's classical
constant is left out.
`RelaxationAbsorption::evaluate_iso9613_rates_held(medium)` evaluates the
same reservoirs for N2 and O2 in any proportion (with Ar, H2O, and CO2 at or
below ISO's ceiling) with ISO's relaxation frequencies at the medium's own T, p
and water fraction HELD: a declared assumption, since the partner rates of any
other proportion are not sourced. At ISO's own fractions it is ISO's rate
law. Any other species present refuses by name.

The law also refuses:
- a present species without data;
- a present collision partner without a rate;
- T outside a rate's validity;
- more than 32 reservoirs;
- relaxing heat capacity above translation.

G0 results (`tests/relaxation_absorption.rs`):
- zero reservoirs equal Stokes-Kirchhoff bitwise;
- the single-process `alpha lambda` peak equals its closed form (2.6e-15);
- `alpha >= 0` and `c0 <= c(omega) <= c_inf` at 10^4 random states;
- the Kramers-Kronig transform of `alpha` matches the dispersion (1.6e-9).

Against ISO 9613-1's own vibrational terms, with ISO's theta and rates:
- dry air differs by 2.0e-3 (253 K) to 1.0e-2 (323 K). This is the standard's
  small-strength, equilibrium-heat-capacity linearization.
- humid air additionally differs by about -x_H2O, the O2/N2 dilution that the
  standard omits.

The 1e-3 G1/G2 gates against ISO are therefore not met and remain open.
NO-CLAIM: no sourced H2O, CO2 or rotational rate data is bundled, so the law
refuses those species rather than defaulting them.

The universal ambient-medium primitive of the musical-acoustics
doctrine: `GasState::try_new(&GasSpec, T, p)` derives the COMPLETE
tuple (density, sound speed, dynamic viscosity, thermal conductivity,
gamma, specific gas constant, cp, Prandtl, characteristic impedance)
bottoms-up — ideal gas + Sutherland + a declared conductivity model
(USSA-1976 empirical air fit or the kinetic-theory Eucken relation).
`GasSpec::dry_air_ussa1976()` carries the public-domain U.S. Standard
Atmosphere 1976 constants and reproduces its printed sea-level values
(rho 1.2250, c 340.294, mu 1.7894e-5) to printed precision; evaluated
at 20 C it IS the pair the acoustic crates' `air()` conveniences
hardcode (asserted in the fs-couple casebook).
`GasState::stokes_kirchhoff_absorption(ω)` is the classical
viscothermal `α ~ ω²` law from `(μ, κ, γ, c_p)`.
`GasState::iso9613_absorption(air, h_r, ω)` is ISO 9613-1
(molecular + classical) with humidity as an explicit `[0, 1]`
argument; do not add Stokes on top. It needs `air`, an
`iso9613::TerrestrialAir` admission bound to this state. The admission is
made from the gas's explicit composition (`TerrestrialAir::admit`, the
relaxation law's own policy) or with the legacy moist-air state
(`TerrestrialAir::moist_air`, U.S. Standard Atmosphere 1976 dry air). Any
other composition is refused by `Iso9613Refusal`, which names the species:
a foreign species (He, SF6, ...), O2 or N2 off ISO's dry fraction, or CO2
above 1e-3. `iso9613::iso9613_terms` exposes ISO's own decomposition, its
classical coefficient and each relaxing species' strength and `f_r`.
Consumers that need ISO's dispersion (its exact Debye pairs) build it from
there, not from a second transcription. `iso9613_relative_humidity` inverts
ISO's humidity-to-water conversion. Pr is DERIVED, never an
input — the emergent ~0.71 for air is the built-in kinetic-theory
falsifier, and the Eucken-vs-USSA divergence (measured 0.1%..12% over
200..600 K, monotone-worsening) is asserted as the relation's honest
boundary. Validity window T in [50, 2000] K, p in (0, 1e7] Pa refuses
outside; the calorically-perfect gamma degrades above ~600 K
(vibrational excitation, order 1-2% there), densities near the 100 bar
ceiling carry a 1-3% real-gas error at 300 K (worse when cold),
and the USSA/Sutherland fits are EXTRAPOLATED beyond their published
bands (~187-288 K atmosphere tables; ~170-1900 K Sutherland-air) —
all documented as accuracy boundaries so extreme-regime studies
evaluate with honest error language. NO-CLAIM: phase validity is not
checked — an ideal-gas model cannot know about condensation, so a
returned state is the gas-phase extrapolation and never evidence the
gas phase exists at that (T, p); callers own the phase check until
vapor-pressure/melting taxonomy rows land. The Eucken 5/4 coefficient
is pinned ABSOLUTELY by the monatomic identity Pr = 2/3 (exact,
mu-independent), alongside the relative Eucken-vs-USSA divergence
envelope.

`resolve_sutherland_gas_state` connects a material card to this same evaluator.
It resolves five exact parameter keys at a shared, explicitly typed temperature
and pressure, retaining the complete source bundle and receipts. Reference
viscosity is converted back to Sutherland beta; it is not mistaken for viscosity
at the operating temperature. Reference temperature is absolute, Sutherland S
is an interval, and conductivity model selection is explicit. The returned
`ResolvedGasState` exposes the selected model separately; its parameter-bundle
identity does not claim to bind that model choice. Missing, ambiguous, mistyped,
out-of-source-domain or invalid parameters refuse atomically. The
`air-dry-ussa1976` seed defines a narrower 273.15–313.15 K / 80–110 kPa dry
engineering application range, with its 1976 composition and model limits.
This is source-backed model evaluation, not real-gas or humidity qualification.
The G1/G3 `g1_g3_sourced_dry_air_reaches_acoustic_loss` integration case in
`xtask/tests/matdb_pack_cli.rs` exercises real source compilation, persistent
store reopening, discovery, this adapter and the existing cylinder air-loss
law at four T/P states, with independent references, receipts, replay and
source-sensitivity/refusal controls. This evidence is numerical, not empirical.

`moist_air::resolve_moist_air_state` resolves the carrier and dilute water-vapor
cards separately at the same typed temperature and total mixture pressure,
then uses the shared moist-gas evaluator. Both parameter bundles and all ten
usage receipts remain accessible; relative humidity, vapor mole fraction and
vapor mass fraction are distinct outputs. Explicit component flags and normal
source applicability prevent silently substituting another composition.
Neither component identity alone identifies the mixture. USSA carrier
conductivity, Eucken vapor conductivity and Wilke/Wassiljewa–Mason–Saxena
mixing are explicit fixed model choices. The component evaluation does not
assert stable pure steam at the mixture's total pressure. This is bounded
source-backed engineering model evaluation, not condensation or a real-gas EOS.

### `visco` — viscoelastic damping tiers (bead ybc75)

Four tiers, one runtime path. `RayleighDamping` (tier 0: αM + βK
coefficients and the two-point fit); `FractionalZener` (fitting
canonical, 4 params, `Real`-generic so exact AD parameter tangents share
the value's code path); `GeneralizedMaxwell` (runtime Prony series with
the unconditionally stable exact-exponential recursion for
piecewise-linear strain and a per-step work/stored/dissipated energy
ledger; `matching_loss` builds one branch that hits an authored
`η(ω₀)`); `lower_to_prony` (deterministic relaxation-time ladder + NNLS,
emitting a `LoweredModel` whose certificate is the MEASURED supremum
relative modulus error over a dense in-band grid and whose evaluation
REFUSES outside the certified band, `FS-MAT-VISCO-OUT-OF-BAND`);
`ThermoelasticZener` (closed-form beam/plate thermoelastic loss,
  plus handbook `aluminum` / `structural_steel` constructors);
`fit_fractional_zener` (deterministic Levenberg–Marquardt over an
admissibility-preserving reparameterization, exact AD Jacobian, and a
caller-AUTHORED residual budget `tol_rel`: a stationary point above the
budget refuses with `FS-MAT-VISCO-FIT-DIVERGED` instead of returning a
model the data cannot support). The e2e evidence lane is the
`visco_casebook` binary: synthesize (deterministic FNV-keyed ±0.1%
perturbation, no RNG state) → fit → lower → single-DOF velocity-Verlet
ring-down → re-measure η from the decay envelope, one JSON line per
stage, bit-identical across reruns except `elapsed_ms`.

`visco::tensor` (`TensorProny`, `PronyBranch`, `TensorPronyState`, `Mandel6`, IRP M0.2-W2a) extends the runtime Prony model to passive anisotropic 6×6 Mandel tensors (`[xx, yy, zz, √2 xy, √2 yz, √2 zx]`). Admission requires a symmetric positive-definite equilibrium stiffness `C_eq ≻ 0` (diagonally scaled Cholesky + deterministic Jacobi eigendecomposition via `fs_la::eigen::jacobi_eigh`), an optional positive-semidefinite Kelvin–Voigt tensor `kv_pa_s ⪰ 0`, at most 64 relaxation branches `(C_j, τ_j)` with `C_j ⪰ 0` (eigenvalues verified non-negative without clipping) and `τ_j > 0`, and a finite applicability band `band_rad_s`. Each branch is factored once at admission as `C_j = L_j L_jᵀ`. Constructors include general `TensorProny::new`, shared-eigenbasis `TensorProny::kelvin(c_eq, per_mode)`, and `TensorProny::isotropic(k, g)` (whose pure shear and pure volumetric histories reduce bit-for-bit to scalar `GeneralizedMaxwell::step`). Frame rotation `TensorProny::rotated` uses the single crate-wide Mandel Bond rotation owner `tensor::rotate_mandel` (`C → Q_6 C Q_6ᵀ`, shared with `fs-solid::linear3`). `TensorProny::plane_stress_complex` computes the frequency-dependent complex Schur complement `(Q′(ω), Q″(ω))` in any ordered `MaterialPlane` and checks passivity (`Re Q ≻ 0`, `Im Q ⪰ 0`) at `ω`; it does not claim a frequency-independent branch-wise plane-stress Prony series (lowering a plane-stress Schur target to reduced-plane Prony branches is M0.2-W2b's responsibility). Time stepping `TensorProny::step` advances `q_j ∈ ℝ⁶` componentwise via `prony_ramp` and integrates dashpot heat by polarization as a sum of squares `A_j ‖L_jᵀ (q_{j,0} + h_j δ)‖² + R_j ‖L_jᵀ δ‖² ≥ 0` with no clamp, closing `work = Δψ + dissipated` to `≤ 1e-12` relative. Intrinsic material loss only is represented here; thermoelastic, contact, support, and acoustic radiation losses are separate owners and are not folded into `TensorProny`.

`visco::tensor::lower` (`TensorComplexModulus`, `KelvinFractional`, `ComplexEngineering`, `PlaneStressTarget`, `TensorLoweringTolerance`, `TensorCertificate`, `LoweredTensorModel`, `lower_tensor_to_prony`, IRP M0.2-W2b) lowers causal frequency-domain anisotropic complex-modulus targets `C*(ω) = C′(ω) + i C″(ω)` to passive `TensorProny` laws over `[f_lo_hz, f_hi_hz]`: 1×1 scalar `FractionalZener` targets delegate to `lower_to_prony` bit-identically; shared-eigenbasis `KelvinFractional` targets solve Lawson–Hanson NNLS (`visco::nnls`) per Kelvin mode on the shared log-spaced `τ` ladder; general non-commuting targets (`ComplexEngineering`, `PlaneStressTarget`) solve PSD-constrained weighted least squares in a simultaneously diagonalized congruence basis at the band midpoint with `jacobi_eigh` PSD cone projection inside the iteration loop (`ViscoError::FitDiverged` on cap exhaustion, no post-fit clipping). `TensorCertificate` certifies `grid_sup_storage + lipschitz_term` and `grid_sup_loss + lipschitz_term` on a 512-point log-spaced verification grid using generalized eigenvalue relative quadratic-form errors `|λ(C_fit − C, C)|` plus the inter-sample Lipschitz enclosure `L · Δln(ω) / 2`, refusing with `ViscoError::FitTolerance` when either bound exceeds `TensorLoweringTolerance`.

## Invariants

1. **Tangent consistency (the merge gate)**: every law's tangent matches
   central FD of its own stress to ≤1e−5 relative across elastic branch,
   plastic branch, cyclic states, and 9×9 hyperelastic components.
2. **Frame indifference**: isotropic small-strain σ(QεQᵀ) = Qσ(ε)Qᵀ;
   hyperelastic P(QF) = Q·P(F) — randomized rotation battery.
3. **Return-map consistency**: after every J2 update, the yield function
   at the returned stress satisfies f ≤ tolerance; dissipation increments
   σ:Δεₚ ≥ 0 (associative flow), total cycle dissipation > 0.
4. **Hysteresis fixture behavior**: M-P virgin curve approaches the b·E₀
   asymptote, tangents E₀/b·E₀ at the extremes, reverse branches soften
   below the elastic line (Bauschinger), symmetric cycles dissipate;
   Mander peaks exactly at (ε_cc, f′cc) with slope 0, softens post-peak,
   unloads to the residual strain, reloads rejoining the envelope.
5. **Calibration round-trip**: synthetic bilinear data recovers E within
   1%, H within 5%, σ_y within 2%, truth inside the fitted envelope.
6. **Inadmissible parameters refuse at construction** (ν bounds,
   compliance definiteness, Ec > Esec, b ∈ [0,1), positive yield).

## Error model

`MaterialError`: `Parameters`, `State` (e.g. det F ≤ 0), `Calibration`.
Out-of-calibration-domain USE is not an error — it is flagged through
`Evidence.model.in_domain` so upstream policy decides.

## Determinism class

**D0**: pure f64 arithmetic with fs_math::det transcendentals; no
iteration counts depend on ambient state (radial return is closed-form).

## Cancellation behavior

All updates are closed-form, allocation-free; P7 by boundedness.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

- Moist air (music bead 3ez8g.3.5): `GasState::try_new_moist_air`
  derives the ideal dry-air + water-vapor mixture from relative
  humidity — breath vs hall as physics, not a detune knob. Both
  mixture rules are EXACT for an ideal mixture (molar mass linear in
  `x_w`; the isochoric-heat identity `1/(γ−1) = Σ x_i/(γ_i−1)`); the
  RH → vapor-fraction step is the Buck 1996 liquid-water fit
  (`e_s = 611.21 exp((18.678 − t/234.5)t/(257.14 + t))` Pa, quoted
  −20..+50 °C, refused outside when RH > 0). Water-vapor spec
  provenance: M = 18.01528e-3 kg/mol, γ = 33.590/(33.590−8.31432)
  from NIST-JANAF H-064 cp(H2O g, 298.15 K) = 33.590 J/(mol K).
  The COMSOL/White steam reference μ = 1.12e-5 Pa s at 350 K and
  S = 1064 K gives β = 2.4186073348106757e-6. These correct the
  earlier rounded/misattributed Cp and β. The dry limit RH = 0 is
  the SAME CODE PATH — bitwise, so no golden anywhere moves.

- Cane packs (music bead 3ez8g.3.6): the reed card's material
  population, licensing-first — three CC-BY primaries in
  data/matdb/seed-v1 (reed-BLANK density + E_L with the radial cortex
  gradient named, Karczewski 2025 MDPI; reed-HEEL E_L/E_T dry AND
  water-soaked, Ukshini & Dirckx 2020 MDPI; whole-wall DMA storage
  modulus + logarithmic decrement, Weidenfeller 2018 SciELO). Blank
  and heel measure DIFFERENT zones of one plant (cortex-weighted vs
  inner) — recorded separately, never averaged; moisture conditions
  recorded per row (unstated conditions carry a
  moisture_content_known 0 flag). NO licensable shear modulus or E_R
  exists (Obataya/Taillard/Kawasaki are publisher-copyrighted) — the
  absence is asserted, not papered over.

- Tissue packs (music bead 3ez8g.3.7): the vocal-fold/lip population
  for the reduce lab and glottis islands — SIX CC-BY packs where a
  refusal was expected: human vertical stiffness gradient
  (3.2/5.5/7.0 kPa indentation with full protocol), porcine dynamic
  modulus across the phonation band (postmortem time typed per row),
  the felid shear rheometry loss-factor find (tan delta 0.33/0.15 —
  the only licensable native-tissue damping), the lone in-vivo lip
  stiffness (45 kPa WITH the hard-teeth-in-series caveat), the
  canonical two-mass MODEL CARD behind a register firewall
  (model_-prefixed names, declared not-constitutive), and the
  STATED-not-measured soft-tissue density register (1.04 g/cm3 as
  stated; the measured Perlman-Titze chain is unlicensable and no
  pack claims it). Method-scale disjointness is a first-class fact:
  Pa-scale matrix shear vs kPa-scale indentation registers are
  recorded separately, never averaged.

- Porous-absorber facility (music bead 3ez8g.5.4): `porous` module
  — JCA(L) equivalent fluid (five transport parameters + optional
  Lafarge thermal permeability; `None` degrades to Champoux–Allard
  with the degeneracy EXPLICIT: k0' = φΛ'²/8) yielding effective
  density/bulk, characteristic impedance/wavenumber, rigid-backed
  surface impedance, and normal-incidence absorption. Cards are
  licensing-first: Nguyen 2024 measured melamine characterization
  (CC-BY, uncertainties quoted in provenance) and the FOAM 02
  Basotect/Pinta identified rows (CC-BY, thick-sample fits only —
  thin-sample rows pin parameters at fit bounds and are refused as
  card material). ORACLE: the committed FOAM 02 measured mean
  spectra (machine-readable CC-BY data; no figure digitization) —
  the model reproduces measured absorption to worst 0.041 / mean
  0.018 (Basotect 60 mm) and 0.029 / 0.006 (Pinta 50 mm) inside
  authored 0.06/0.03 bands. `porous::gas::PorousGasLaw` retains the
  original card, provenance, frequency support and explicitly supplied
  characterization-gas viscosity. Intrinsic permeability
  `k0 = mu_ref / sigma_ref` stays fixed; current `GasState` viscosity
  supplies `sigma = mu / k0`. Density, sound speed, viscosity, gamma,
  pressure and Prandtl number come from that same coherent gas snapshot,
  then the original JCA(L) evaluator. `fs-duct::network::NetworkLoad::PorousRigid`
  consumes the rigid-backed response as a normal-incidence boundary,
  with the exact phasor and exposed-area conversion. Source/geometry
  refusals remain explicit. No porous constitutive model is duplicated.
  The skeleton remains rigid and gas-filled; gas substitution does not
  establish liquid saturation, condensation, swelling, moving-frame waves,
  distributed lining or transient deposition of absorbed acoustic heat.
  Focused gas-source and duct-cylinder tests exercise data-only changes of
  temperature, pressure, gas mixture and pore permeability, including
  independent acoustic power balance through both network load paths.

## Conformance tests

`tests/active_tissue.rs` (suite `irp-s0-materials`, bead M0.3-W3): zero
activation bitwise equal to the passive law (1000 random F, stress and
tangent), the isometric Cauchy stress `T_max a f_L` against closed forms of
three force-length fixtures on a 50-point stretch grid (1e-12, orthogonal
components exactly zero), the traction-free free-contraction root at prescribed
J = 1 against an analytic balance (1e-10, three activations), the active-work
ledger `W_total = Delta Psi_passive + W_active` over a time-varying activation
(1e-10), objectivity, the finite-difference tangent of the active stress and of
the total, and the refusals above at the bound and just beyond it.

`tests/conformance.rs` (JSON verdicts, suite `fs-material/conformance`):
mt-001 FD tangent gate (every law, cyclic states); mt-002 objectivity;
mt-003 return-map consistency + dissipation; mt-004 hysteresis fixtures;
mt-005 calibration round-trip (+ degenerate refusals); mt-006 rank-one
convexity sampling for NH/MR, constructor refusals, card completeness,
Evidence domain flagging.

Visco (inline in `src/visco.rs`, 14 tests): closed forms over six
decades, α→1 classical-Zener reduction, Prony recursion vs direct
convolution at 1e-12, dissipation-ledger closure with sign-flip
mutation, ≤8-term plateau lowering with out-of-band refusal and
dropped-term mutation, modal injection vs exact quadratic roots, AD
gradients vs finite differences, bitwise-deterministic lowering,
dt→0 elastic-jump/small-dt expm1 regression, LM fit ground-truth
recovery from an offset start with negated-loss mutation refusal,
degenerate-input refusals, and the η ≥ 0 / monotone-relaxation property
sweep. `visco_casebook` is the end-to-end lane (fit noise floor,
certified lowering, ring-down η within an authored, measured gate).

`tests/prony_tensor.rs` (9 tests, suite `irp-s0` / `prony_tensor`):
PD/PSD/τ admission refusals, 64 random multiaxial energy-closure and
non-negative dissipation histories, bitwise isotropic reduction to
scalar `GeneralizedMaxwell`, Kelvin constructor vs general assembly,
SO(3) Bond rotation spectrum/dissipation invariance, Kramers–Kronig
consistency, Fourier-quadrature relaxation vs complex modulus,
plane-stress complex Schur complement passivity and elastic limit,
and out-of-band refusal.

`tests/prony_tensor_lowering.rs` (8 tests, suite `irp-s0` / `prony_tensor_lowering`):
1×1 scalar bitwise reproduction of `lower_to_prony`, KelvinFractional
shared-eigenbasis per-mode lowering within storage/loss certificates,
orthotropic `ComplexEngineering` non-commuting lowering with PSD
branches, `PlaneStressTarget` Schur-complement lowering, 65,536-point
independent dense-grid enclosure check, branch-removal mutation breaking
the certificate, term-count/band/tolerance/divergence refusals, and
bitwise determinism across repeated runs.

`src/gas.rs` moist battery — the two-source saturation rule (Buck vs
the independent Alduchov–Eskridge Magnus constants, ≤0.5% across the
window; IAPWS absolute pins 611.657 Pa at the triple point and
2339.2 Pa at 20 °C, measured 611.65/2338.3); the OIML R 111 / CIPM
conventional-density anchor (ρ = 1.2 kg/m³ DEFINED at 20 °C, 101325
Pa, 50% RH; measured 1.19885 — an independent metrology-community
check on the whole mixture chain, and humid air comes out LIGHTER);
the bitwise dry limit; mixture rules re-derived in the test from raw
literals (the anti-tautology discipline: never compare two quantities
off the same GasState); RH-monotone c/ρ with the saturated-20 °C
uplift in the published 0.8–1.8 m/s class (measured 1.263, the
Wong-1986/Cramer-1993 ~0.35% class); refusals (RH ∉ [0,1], Buck
window, the x_w ≤ 0.15 transport ceiling) and bitwise repeats.

`tests/cane_pack.rs` g2 gate — positivity + schema hygiene (every
scalar carries uncertainty + validity rows), fiber anisotropy in the
source's own band (E_L/E_T = 10 dry, wider wet), the moisture law
(water softens both axes, transverse more), cross-source RANGE OVERLAP
between the heel's static band and the whole-wall DMA range, the
blank-vs-heel DISJOINT-zones assertion with the gradient named in both
observations, the derived loss factor delta/pi in the woody class with
conversion ownership stated, and the named-absence assertion that no
shear row exists without a source.

`tests/tissue_pack.rs` g2 gate — positivity + schema hygiene, the
human vertical gradient (the source's own finding), method-scale
disjointness with named method dependence in every observation,
postmortem-time separation, loss tangents in the soft-tissue class
with tiger above lion, the lip teeth-caveat presence, the model-card
firewall (model_ prefix + declared register), and the
stated-density assertions (no measured tissue-density row may exist
without a licensable source).

`src/porous.rs` `porous_tests` (pa-001..pa-006): both measured-
spectrum oracles; analytic degeneracies (thin layer -> 0, rigid-
backed compliance limit Zs -> K/(jωd) at 7e-5, air-like card ->
free-medium Zc within the disclosed 0.07% constant mismatch);
refusals by name (parameters, out-of-band frequency, negative
thickness); the composed lined-tube demonstration (a simulated
ISO-10534-2 two-microphone rig recovers the direct absorption to
2e-16, and the melamine-lined tube's fundamental decays at
T60 = 0.13 s where the rigid tube never decays); card validation +
bitwise determinism.

## No-claim boundaries

- The porous facility claims LININGS (surface impedance/absorption
  of a rigid-framed layer), never the hammer's compression stack;
  no elastic-frame Biot waves (equivalent fluid only — frame-wave
  escalation is a named successor); no room acoustics; no oblique
  incidence in v1 (normal incidence + locally-reacting consumption).
  Pinta Plano Polar is polyester-fleece FELT-CLASS nonwoven, not
  wool piano felt — the wool-felt characterization is a recorded
  vein.


- **No intrinsic tissue calibration is inferred from the law family.**
  `tests/hyperelastic_ogden_hgo.rs` gates Ogden-1972's published three-term
  analytic nominal curves under prescribed J=1 with lateral-pressure
  elimination (uniaxial 1..7.5, equibiaxial 1..4.5, pure shear 1..5, 1e-3)
  and HGO one-/two-family incompressible Cauchy oracles (1e-6). Separate
  traction-free finite-bulk roots measure the O(1/kappa) modelling residual:
  finite kappa/mu=1e4 is not identically incompressible. The literal 1e-3
  finite-bulk Ogden endpoint band and 1e-6 HGO band at kappa/c10=1e5 are
  contradicted by the exact specified laws; tests retain those counterexamples
  and their derived residuals rather than changing the law or widening gates.
  Treloar observations retrieved/digitized by Bucchi & Hearn (2013, Appendix A)
  are diagnostic only: source uncertainty, redistribution licence and fit band
  are unstated. The M8-W1 ledger rows `ogden-1972-treloar-three-term-parameters`,
  `treloar-1944-points-bucchi-hearn-2013` (fixture-only) and
  `hgo-2000-fig14-media-parameters` record the provenance. The original Ogden
  (1972) paper was not opened: its coefficients were read from an open reproduction
  (Twizell and Ogden 1981, Table I, no unit printed). Literature coefficients are
  not a calibrated tissue card.
- **The active fibre stress is isometric, prescribed and uncalibrated.** There is
  no force-velocity factor, no fatigue or length history, and no activation
  dynamics. The force-length families in `tests/active_tissue.rs` are test
  fixtures, not muscle data; no vocal-fold or lip active curve is gated, and no
  `T_max` is sourced here (none was found for orbicularis oris). The G2 gate
  that reproduces a published force-length or prestress curve with its source's
  own parameters is not shipped yet.
- **Mander cyclic rules are the simplified elastic-unload variant**
  (declared on the card): no cyclic stiffness degradation beyond the
  residual-strain rule, no tension stiffening. Full Mander cyclic rules
  are follow-up work for the fiber-beam bead.
- **J2 cyclic response is bounded**: the opt-in linear Prager component
  translates the yield surface, but supplies no nonlinear backstress saturation,
  calibrated ratcheting, damage accumulation or fatigue-life prediction.
- **No damage/softening 3D laws yet**; failure envelopes are declarative.
- **Rank-one convexity is SAMPLED**, a necessary-condition spot test —
  not a polyconvexity proof (interval-certified convexity belongs to the
  certifier tier).
- **Calibration v0 is bilinear segmented LSQ**; nonlinear (Mander/M-P
  parameter) fitting and posterior envelopes land with fs-io CSV
  catalogs and the UQ stack.
- **Homogenization hooks**: homogenized laws register as ordinary
  `ModelCard`-carrying laws; the unit-cell pipeline itself is the
  lattice bead's scope.
- **Visco lowering certificates are measured, not analytic**: the
  `LoweredModel` bound is a dense-grid supremum on the stated band —
  no global proof, and nothing is claimed outside the band (evaluation
  there refuses). The casebook's ring-down comparison is a NUMERICAL
  consistency check of the fit/lower/step chain against its own
  declared physics; it calibrates no real material, and its η gate is
  authored from an observed run with headroom, not derived from an
  integrator error bound. The fractional→Prony path serves 3D FEM via
  fs-time's damping seam; a spatially assembled viscoelastic FEM
  operator is NOT part of this crate's claim.

## Authority-separated multi-case identifiability schema (bead I10.1)

`identifiability/authoritative.rs` is the current public I10.1 contract. It
models a campaign, not a single coupon, because complementary specimens,
protocols, environments, and observation operators can jointly break a gauge
that no individual case resolves. The schema preserves ambitious continuous,
discrete, mixed, and stratified gauge/theorem targets while refusing to report
their truth before an evidence-bearing assessment exists.

The authority chain is deliberately monotone:

1. `IdentifiabilityProblemDocument` is a canonical but **unresolved** physical
   and statistical question. Raw decode returns only this type.
2. `AdmittedIdentifiabilityProblem::resolve_and_admit` requires concrete
   Context-of-Use/material/model/case artifacts and an exact closed
   `SourceResolutionSet` for every other source. It then mints `ProblemId` and a
   separate `SourceAdmissionId`.
3. `IdentifiabilityExecutionPlan` adds coordinates, parameter actions,
   algorithms, analyzer/build/derivative semantics, seeds, budgets,
   numerical controls, arithmetic policy, initialization, stopping, and
   determinism. Every confirmatory proposition is preregistered as an exact
   `ClaimRequest` together with a sourced dimensionless error metric,
   nondimensionalization, and maximum certified error. The plan also requires
   an exact locally verified authority closure for every transitive execution
   source, and replay requires the caller to present that same closure. It mints
   `ExecutionId` without rewriting `ProblemId`.
4. `IdentifiabilityAssessment` adds product-typed propositions and explicit
   claimed-established/claimed-refuted/claimed-inconclusive/not-assessed
   conclusions with receipts. Its proposition registry must equal the exact
   preregistered execution registry; matching only coarse claim axes cannot
   substitute a new proposition. Its caller-held source authority must agree
   with both problem and execution authority on every overlapping key. It mints
   `AssessmentId` without rewriting the problem or execution. These
   deliberately named claims are not promoted to theorem tokens until a
   method-specific verifier exists.

### Current public schema

The authoritative API and all four independently domain-separated identity
stages are schema **v3**. Decoders require the exact stage version; v1/v2 bytes
are not migrated or reinterpreted implicitly.

- `SourceRef` binds a `SourceKey`, semantic `SourceKind`, expected content hash,
  exact fs-blake3 digest domain, and positive source-contract version. A hash is
  not issuer authentication. Opaque resolution must supply bytes that reproduce
  the digest while retaining the digest domain, preventing cross-domain replay.
  `SourceVerification` records the narrower fact actually established:
  `TypedArtifact`, `HashPreimage { byte_len }`, or `Unverified`; it never
  promotes byte equality into provenance or issuer trust. Synchronous generic
  hash-preimage verification refuses inputs larger than the 4 MiB canonical
  transport bound; larger sources require a future typed or streaming
  artifact-store verifier rather than an unbounded allocation disguised as
  admission.
  Typed Context/V&V/material/model artifacts additionally require their actual
  domain and schema version. `AuthorityDisposition` distinguishes byte-content
  verification, an unvalidated external trust-policy receipt, and an explicit
  unresolved state; unresolved sources cannot mint `ProblemId`.
- `StudyParameter` separates decision purpose (`Estimand`, `Nuisance`,
  `Hyperparameter`, `CalibrationControl`) from inferential treatment
  (`Estimated`, `Profiled`, `Marginalized`, source-bound `Conditioned`, or a
  DAG-checked `Derived` value). Physical quantity/domain/prior, concrete owner
  payload, population/case/specimen/field/hierarchical scope, and honest
  connectivity coverage remain coordinate-free. Case applicability is explicit:
  material-lot, field, and hierarchical scopes carry nonempty exact case sets;
  specimen scope names its one case; and only `Global` implies every campaign
  case without enumerating it.
  `ParameterPrior::Uniform` denotes a continuous uniform measure and therefore
  requires finite, positive-width support contained in the physical domain.
  Singleton support, including signed-zero aliases, refuses rather than silently
  manufacturing atomic mass. A future discrete/Dirac prior family may represent
  atoms explicitly without weakening this continuous-measure invariant.
- `JointConstraint` supports dimension-checked affine, simplex, ordered,
  external-manifold, and stochastic-coupling domains. This prevents an
  apparently valid Cartesian product of marginal domains from silently
  admitting physically impossible parameter combinations.
- `AdmissibleDomainWitness` supplies one finite value for every parameter and
  proves that the *same global point* satisfies every locally evaluable
  marginal, conditioned, affine, simplex, and ordered predicate. Derived
  parameters or opaque external/stochastic constraints require exactly one
  `OpaqueDomainMembershipClaim`. Its admitted
  `AdmissibleDomainCertificate` source has an exact kind/domain/version, and
  its binding commits to the witness, all parameter semantics, the complete
  constraint conjunction, and the exact `SourceRef`s for opaque joint
  constraints. Derived-definition source keys remain bound through the
  parameter semantics and the enclosing problem source registry. This closes
  the former gap where locally evaluable constraints could each look feasible
  while their intersection was empty. The opaque certificate remains a
  source-bound claim, not a theorem minted by fs-material.
- `StudyCaseDocument` binds case purpose, exact initial state, specimen,
  protocol, forward model, prospective or retrospective data intent,
  observations, and discrepancy semantics. `CasePhysicsSources` closes every
  embedded frame, geometry, process, preparation, load, environment, time-grid,
  and explicit-initial-state digest to an exact source kind, role-specific
  domain, contract version, and hash. `ObservationKey { case, channel }`
  prevents local channel-name aliasing across cases.
- `StudyObservation` carries exact QoI/unit/`QuantitySpec`, frame, graph port,
  observation/aggregation/sensor sources, instrument, acquisition channel,
  clock, versions, marginal noise, missingness, saturation, and prospective or
  exact raw-row semantics. At source admission, each QoI/unit is closed against
  the concrete Context of Use;
  every retrospective row is closed against the re-derived `DataLineage` and
  must match its manifest QoI, instrument, acquisition channel, and clock.
  Every `ObservationSourceRef` additionally binds the exact dataset-byte hash,
  locator domain and positive locator-contract version, immutable locator
  hash, and extraction-receipt hash. Distinct manifest IDs cannot alias one
  receipt-independent raw locator identity: changing only an extraction
  receipt cannot manufacture a second observation or cross a holdout
  partition. Receipts still move the complete typed source and artifact
  identities. Every row's dataset hash must equal the experiment authenticity
  source-byte hash. The exact instrument must occur
  exactly once in the current experiment roster, its calibration-certificate
  hash must equal the sensor source, and the observation/protocol clock must
  occur in the experiment clock topology. A blind-falsification case also needs
  an exact `BlindReleaseReceipt` accepted by
  `CalibrationSplit::blind_selection`; its authority receipt moves the source
  admission identity rather than the physical question.
- `JointNoiseModel` keeps marginals and dependence separate. Dense correlation
  requires the exact composite channel set and marginals with finite standard
  deviation; bounded/unknown noise is never assigned a fictitious Gaussian
  scale. `Independent` is not a default: it requires an exact source-bound
  `Assumption`. External kernels and honest unknown dependence remain
  representable.
- `ObservationSharingGroup` closes intentional within-case reuse of one raw row
  by naming every repeated row, its exact consumer-channel set, one
  joint-likelihood source, and a justification. The row coverage is exact:
  omitted repeats and declarations for rows that are not actually repeated both
  refuse.
- `DataReusePolicy` represents an explicit pairwise-disjoint raw-data policy or
  explicit sharing groups; there is no inferred reuse default. Admission
  compares exact raw-byte digests, manifest digests, and immutable row-source
  identities, so distinct experiment wrappers cannot hide reuse. Intentional
  reuse requires a non-overlapping `DataSharingGroup`, exact case membership,
  joint-likelihood source, and human justification. Every within- and cross-case
  sharing declaration must name the one global `DenseCorrelation` or
  `ExternalKernel` likelihood; sharing under `Independent` or `Unknown`
  dependence refuses.
- `DistributionFunctional` names location, log-scale, cross-channel
  correlation, missingness-logit, and censoring-logit targets. Derivative units
  are derived from functional and physical-parameter quantities; callers cannot
  supply contradictory units. `InfluenceRepresentation` supports direct,
  state-mediated, composite-DAG, and externally defined mathematical routes;
  computational derivative providers live only in the execution plan. No route
  carries a nonzero/proof receipt. Every influence endpoint must lie in a case
  where its parameter is active, and an influence claim's case scope must
  contain every endpoint and parameter reached through the full transitive
  composite-input DAG, not merely the final functional's endpoint.
- `GaugeDeclaration` represents continuous, discrete, mixed, stratified, and
  explicitly suspected gauges plus quotient, local-section, slice, retained,
  unresolved, or continuous-reduction-with-discrete-residual handling. A
  `GaugeSlice` binds its exact joint constraint, fixed codimension or
  stratum-indexed codimension-profile source, and an exact `GaugeSection`
  transversality/coverage source. Its constraint support must be a nonempty
  subset of the gauge members; internally understood equality/manifold
  codimensions must agree with the declaration, while inequality/order geometry
  requires the external coverage declaration. Mixed gauges can name an exact
  reduction of only their continuous component plus the induced residual
  discrete action without conflating that state with a full quotient.
  `GaugeCompositionDeclaration` is an exact hyperedge over a simultaneously
  active gauge system. It binds the complete member set, an authority-owned
  composition law, `IndependentProduct` or `Generated` semantics, effective
  algebra/orbit geometry, status, and validity cells. Every simultaneously
  active assumed system needs exactly one matching assumed hyperedge; missing,
  partial, duplicate, or scope-incompatible composition refuses instead of
  inheriting iteration order. Independent products must preserve the member
  invariants, while generated actions may declare semidirect, ordered, bracket,
  or richer interaction geometry within the union of coordinates their members
  can move. This schema is deliberately fertile ground for future groupoid,
  Lie-algebroid, sheaf, and cohomological composition theorems; the declared law
  is their substitution-resistant input, not their proof.
  Claim compatibility is fail-closed: a
  raw `Unique` fiber cannot intersect any declared gauge; `DiscreteOrbit`
  requires a purely discrete gauge containing every claimed parameter;
  `MixedOrbit` requires a mixed gauge whose continuous component has an explicit
  reduction and whose residual discrete action is source-bound;
  `QuotientUnique` requires that containment plus a non-suspected gauge and an
  explicit full quotient or full-coverage slice (a partial mixed reduction is
  insufficient); and stratified fibers require a typed
  stratification source. These checks prevent internally contradictory claim
  shapes without asserting that the declared gauge or fiber theorem is true.
- `ParameterExecutionAction` covers every physical parameter exactly once and
  must agree with its declared treatment. Built-in coordinates are checked as
  dimensionally valid full-domain bijections. Numerical rank/singular-value
  controls carry an exact `Nondimensionalization` source. All Five Explicits,
  claim requests, numerical scaling, and rank/conditioning/arithmetic choices
  are execution-identity inputs.
- `TypedIdentifiabilityClaim` is a product of information regime, local/global
  extent, fiber geometry/cardinality, mathematical quantifier, scalar/algebraic
  domain, subject, and campaign/case/stratum scope. Complex extensions,
  mixed/stratified domains, derived functionals, and stratified fibers carry
  exact admitted source definitions. Finite, orbit, quotient, positive-
  dimensional, and stratified fibers remain representable independently of
  local/global extent. These axes are not collapsed into an ordinal evidence
  color.
  Every applicable information/extent/genericity request is conjunctive rather
  than precedence-selected. Cross-axis validation refuses finite-data claims
  over prospective cases, parameters inactive throughout the requested cases,
  transitive influence endpoints outside the requested scope, and gauge-fiber
  mismatch.
  A posterior request additionally carries a `ProbabilityMeasure` joint-prior
  source and requires a declared distribution prior for every free inferential
  parameter active in any requested case. Quantifier sources are narrow:
  `QuantifierRealization` for point claims, `ReferenceMeasure` or
  `ProbabilityMeasure` for almost-everywhere claims, `QuantifierDomain` for
  universal claims, and `ProbabilityMeasure` for probability thresholds.
- `ClaimRequest` binds the complete `TypedIdentifiabilityClaim` proposition and
  a `DimensionlessErrorPolicy` into execution identity. The policy carries
  exact `DimensionlessErrorMetric` and `Nondimensionalization` sources plus a
  finite nonnegative maximum certified error. `requested_axes()` is only a
  derived planner projection; it is not independent authority or identity
  state. A positive/refuting `ClaimAssessment` must reproduce the exact
  preregistered proposition, analyzer, metric, and nondimensionalization, carry
  a separately byte-verified evidence receipt, and report a finite nonnegative
  certified error bound no larger than the preregistered maximum. Inconclusive
  and not-assessed states retain explicit reasons.
- Four owner-local identity declarations register the problem,
  source-admission, execution, and assessment surfaces independently. Each has
  its own version/domain/magic, exhaustive no-`..` field classifier, semantic
  bindings, dependency edge, mutation evidence, and coupling surface.
  `ArtifactHeader.id` remains part of exact execution/assessment transport as a
  ledger label, but is deliberately excluded from their scientific identity
  projections; all other Five Explicits remain identity-bearing.

### Current invariants

1. **No authority from transport:** problem decode re-runs structural admission,
   rejects duplicate/trailing/noncanonical encodings, and still returns only an
   unresolved document. Execution and assessment decode both require their
   exact caller-held, locally verified source sets; serialized verification
   markers cannot mint authority. Problem identity exists only after exact
   source resolution.
2. **Identity non-interference:** coordinates, algorithms, seeds, budgets,
   exact claim requests/error policies, and builds move only `ExecutionId`;
   conclusions and receipts move only `AssessmentId`; trust-policy receipts
   move `SourceAdmissionId` but not the physical `ProblemId`. A ledger
   artifact-label mutation changes exact transport but does not change
   `ExecutionId` or `AssessmentId`; execution source-authority changes do
   change `ExecutionId`.
3. **Multi-case closure:** all parameter scopes, observation endpoints, claim
   scopes, observation- and data-sharing groups, discrepancy roles,
   constraints, influences, and gauges reference the exact canonical campaign
   graph. Non-global ambiguous parameter scopes carry exact case applicability;
   influence endpoints outside parameter applicability refuse. Unreachable
   source registry entries refuse instead of perturbing identity.
4. **Model/source closure:** concrete material/model membership, parameter
   roster/quantity/nominal bounds, initial-state policy, state/protocol version,
   Context-of-Use QoI/unit, experiment QoI membership, instrument/calibration
   binding, clock-topology membership, exact split partition, blind-row source
   binding, blind-release split/hash/commitment authority, and cross-case
   raw-row reuse are re-derived from concrete artifacts rather than trusted from
   problem bytes. Shared concrete keys and shared blind releases must agree
   exactly instead of taking last-writer-wins state.
5. **No hidden defaults:** prospective/retrospective data, missingness,
   discrepancy, dependence, prior absence/not-applicability, connectivity
   absence, gauge handling, and every claim conclusion are explicit variants.
6. **Dimensional closure:** physical parameters retain exact `QuantitySpec`;
   affine joint constraints check coefficient-times-parameter dimensions;
   simplex/ordered constraints require exact common quantities; influence
   derivatives are derived; numerical controls and claim error bounds bind
   exact nondimensionalization semantics.
7. **Global domain closure:** one canonical witness covers every parameter and
   the full joint-constraint conjunction. Built-in predicates are evaluated at
   that point; opaque predicates require the exact full-conjunction membership
   binding, so individually plausible but jointly contradictory domains cannot
   mint a problem identity.
8. **Dependence closure:** independence is explicit and source-bound. Repeated
   row consumption and cross-case provenance overlap require exact sharing
   declarations, and every sharing declaration must agree with the one global
   joint-likelihood source.
9. **Confirmatory-claim closure:** execution owns a bounded nonempty registry of
   exact propositions and sourced dimensionless error policies. Assessment must
   reproduce that registry exactly; its metric and nondimensionalization must
   equal the request, and its certified error cannot exceed the preregistered
   maximum. Posterior claims bind a joint probability-measure source and require
   distribution priors for every free parameter whose applicability overlaps
   the claimed cases. Quantifier roles accept only their narrow source kinds.
10. **Graph closure:** derived parameters and composite influences must be DAGs;
   every free parameter has a declared route or an explicit no-connectivity
   reason; modeled-discrepancy parameter use must remain inside the parameter's
   exact case applicability; influence-claim support closes recursively over
   every composite input; self-correlation and incompatible
   missingness/censoring functionals refuse.
11. **Bounded deterministic transport:** maps/sets and commutative affine terms,
   symmetric correlation endpoints, sharing groups, and dense-correlation
   permutations are canonicalized; counts and text are bounded; floats retain
   exact bits with signed-zero normalization; schemas reject stale/future
   versions. The current public authority schema is v3; v1/v2 transports fail
   closed rather than being reinterpreted under the expanded domain-witness,
   case-source, row-lineage, sharing, scope, and exact-claim semantics. Generic
   synchronous source-preimage verification is bounded by the same 4 MiB
   ceiling.

### Current verification and observability

`tests/identifiability_authority.rs` and
`tests/identifiability_retrospective.rs` emit deterministic JSON diagnostics and
cover canonical problem/execution/assessment round trips, multi-case
ordering and composite endpoints, source/dataset closure, constructive
admissible-domain fixtures, source-authority separation, identity
non-interference, disconnected parameters, derived and composite cycles,
unit-invalid constraints, dense-correlation marginal rules, accidental and
declared within-/cross-case raw reuse, explicit generated and independent gauge
composition hyperedges, typed slice support/codimension/coverage, pure-discrete
versus reduced-mixed residual orbits, modeled-discrepancy case scope,
transitive composite-claim closure, treatment/action coverage, product-typed
claims, real ExperimentArtifact/CalibrationSplit
admission, partition leakage, experiment-QoI/instrument/channel/clock closure,
blind-release absence/replay/commitment/purpose/conflict cases, transitive
problem/execution/assessment authority agreement, and adversarial
transport/version inputs. `fs-evidence/tests/vv.rs` independently exercises
typed row-source injectivity, dataset anchoring, and codec/identity mutations.
Identity tests also exercise stage-domain/magic separation, caller-order
canonicalization, artifact-label exclusion, and exact stage version refusal.

The v3 verification gate additionally requires adversarial mutations for exact
claim substitution, metric/nondimensionalization substitution, certified error
above the preregistered maximum, broad or mismatched quantifier kinds, missing
joint-prior authority, a missing prior on any overlap-active free parameter,
ambiguous-scope case-set drift, and influence endpoints outside parameter or
claim applicability. Omitting any of these mutation classes is a verification
gap, not evidence that the corresponding invariant is vacuous.

The schema/codec path is deterministic D0, has bounded transport and collection
ingress, performs no I/O or parallel work, and contains no `unsafe`. Its final
aggregate structural cap is defense in depth, not yet a complete pre-allocation
work ledger for every nested public input. Downstream symbolic, numerical,
algebraic, differential, sheaf/cohomology, and theorem-discovery analyzers are
not covered by that bounded-work exemption: they must be cancellable
asupersync tile programs and finalize receipts only after drain.

### Current no-claim boundaries

- Structural admission proves closure of the declared question, not
  identifiability, observability, a nonzero sensitivity, laboratory validity, or
  model adequacy.
- A `GaugeSlice` coverage source and a residual discrete-action source are exact,
  authority-checked semantic inputs, not internal proofs of transversality,
  global orbit coverage, induced-action well-definedness, freeness, or quotient
  regularity. They make those ambitious theorem obligations preregistrable and
  substitution-resistant; a future sealed gauge verifier must discharge them.
- The admissible-domain witness is a constructive finite-f64 point, not a proof
  of robust interior, positive measure, connectedness, smooth-manifold
  regularity, or absence of disconnected or singular feasible components.
  Built-in predicates are evaluated locally and exactly under the schema's f64
  semantics. For derived/external/stochastic predicates, v3 verifies the
  *structural binding* between the witness, full conjunction, and one admitted
  certificate source; it does not parse that source or prove the external
  membership assertion true. Future [F/M] certifiers may add rational/interval/
  Taylor witnesses, independently checked decision procedures, constructive
  infeasibility certificates, and sheaf-local witnesses with explicit
  obstruction cocycles.
- Content verification proves byte equality only. `ExternalTrustReceipt`
  retains a receipt issued under an external policy that fs-material does not
  authenticate; it does not prove the issuer, experiment, or scientific
  proposition correct.
- Generic `SourceResolution::verify` is deliberately synchronous and bounded:
  opaque preimages above 4 MiB refuse. v3 does not ship a generic typed or
  streaming artifact-store verifier for arbitrary oversized opaque sources, a
  Merkle-range proof protocol, or a resumable verification receipt. Callers
  must not truncate a large artifact and treat a prefix digest as source
  admission; the existing concrete typed artifact admission paths do not imply
  that generic capability.
- Exact `CasePhysicsSources` close embedded digests to source
  kind/domain/version/hash and admitted bytes; they do not yet parse every
  geometry, process, loading, environment, or time-grid source into a common
  authority-grade semantic schema. Byte closure therefore cannot substitute
  for physical-law or laboratory validation.
- `InfluenceDeclaration` is reachability semantics. Only a separately admitted
  assessment may establish/refute nonzero influence, rank, genericity,
  globality, practical resolution, or a gauge theorem.
- A source-bound zero-discrepancy assumption is still an assumption. Missing or
  uncharacterized discrepancy never means zero.
- `DeclaredSyntheticSelfModel` additionally binds the declared producer and
  the complete forward-model `SourceRef` through a domain-separated production
  preimage. That caller-computable consistency relation is not evidence that
  the experiment actually ran that executable or model. Theorem promotion
  requires a future signed execution receipt binding the experiment hash,
  executable/build identity, producer, and exact forward-model source.
- Declared `ArtifactHeader` time and memory budgets are identity-bearing
  execution policy, but v3 does not yet estimate a document's working set or
  admit it against `memory_bytes`. Fixed transport, per-collection, and
  aggregate structural ceilings remain independent safety limits. A future
  deterministic estimator should charge nested allocation and semantic-work
  units before construction and canonical decoding, then bind that admission
  receipt into execution identity.
- Custom `Debug` implementations on direct authority/capability carriers and
  stage aggregates reduce accidental disclosure and keep aggregate output
  bounded. They are not an authorization boundary: explicit getters expose
  exact audit material, some public nested value and error diagnostics retain
  ordinary `Debug`, and v3 does not claim panic-path-wide redaction. A future
  authority-graph redaction certifier should cover every reachable public leaf,
  error path, and maximum formatted size without removing exact audit access.
- The gauge vocabulary intentionally leaves room for powerful new quotient,
  slice, stratification, groupoid, sheaf, and cohomological theorems. The v3
  fiber/gauge compatibility checks eliminate contradictory declarations; they
  do not establish freeness, properness, orbit completeness, quotient
  regularity, slice coverage, stabilizer strata, or a global fiber theorem.
  Future evidence can promote genuinely new theorems without shrinking the
  ambition of the representable claim space.
- The physical problem identity is coordinate-free but not magically invariant
  under arbitrary changes of units, likelihood, prior measure, model family,
  protocol, discretization, or data reuse. A future proved equivalence may bind
  multiple `ProblemId`s in a new theorem receipt rather than silently changing
  v3 identity semantics.
- Retrospective v3 re-derives the experiment/split lineage directly. It does not
  yet require an `AdmittedVvCase`, so full validation-plan, physical-referent,
  solution-verification, and assumptions-ledger closure is a promotion blocker,
  not a current claim. Direct experiment instrument-roster/calibration-hash and
  clock-membership checks are enforced, but synchronized-clock method/skew
  sufficiency and the complete validation-plan topology still require the full
  V&V-case authority chain. Likewise, `UnitId`/`QuantitySpec` agreement is
  declared and identity-bound but cannot be independently checked until the
  Context-of-Use schema carries dimensional unit semantics.
- Raw rows now bind dataset bytes, a versioned locator, extraction receipt, QoI,
  instrument, acquisition channel, clock, and admitted partition membership.
  Locator injectivity deliberately ignores the extraction-receipt hash, so
  reissuing or relabelling a receipt cannot duplicate one immutable raw row;
  the receipt remains artifact-identity-bearing provenance.
  They still do not carry a typed per-sample value, unit, timestamp, spatial
  coordinate/support, uncertainty contribution, parser transformation trace,
  or authenticated extraction-receipt issuer. That is stronger than a row label
  but remains weaker than a complete authority-grade measurement table.
- Reuse closure prevents an observation row or cross-case dataset from being
  silently counted as independent. It does not prove that the declared joint
  likelihood is statistically adequate, identifiable, correctly normalized,
  or robust to misspecification; `Independent` likewise remains an admitted
  assumption rather than a discovered fact.
- Material-lot, field, and hierarchical parameters now carry exact case sets,
  and applicability/endpoint checks enforce structural overlap. Those checks do
  not prove that a named lot was physically used in a specimen, that a field
  support source denotes the intended region, or that a hierarchy source has
  the claimed population semantics. Likewise, narrow quantifier source kinds
  do not by themselves prove that a realization belongs to the claimed scope,
  that a measure is normalized, or that its support equals the declared domain.
- Posterior claims bind a joint `ProbabilityMeasure` source and require
  individual distribution priors on every overlap-active free parameter. v3
  does not independently prove that the joint measure is normalized, that its
  marginals equal those parameter priors, or that its dependence/support match
  the claimed campaign scope.
- `ClaimedEstablished` and `ClaimedRefuted` are substitution-resistant content
  claims, not theorem receipts. A future sealed verifier must bind the full
  claim digest, problem/source-admission/execution identities, method/build,
  exact consumed partitioned row sources, numerical-error policy, issuer,
  checker, and trust policy before exposing an unprefixed theorem verdict.
- Execution v3 preregisters the complete proposition and sourced dimensionless
  error policy, and assessment enforces exact substitution resistance plus the
  declared certified-error ceiling. Source admission still proves only the
  metric/nondimensionalization/receipt bytes: it does not establish that those
  bytes define a mathematically valid dimensionless metric, that the analyzer
  computed the bound correctly, or that the receipt's issuer is trustworthy.
  A sealed method-specific verifier must independently check those obligations
  before the claimed bound can support a theorem token.

## Historical single-case identifiability prototype (non-authoritative)

The following description documents the initial single-case draft for design
archaeology. Its implementation and `StudySpecId`/`PhysicalStudyId` wrappers
live only in repository history; its obsolete prototype integration test remains
explicitly compile-disabled. None of these historical materials are the current
I10.1 public contract or can mint current authority.

The removed `identifiability.rs` prototype owned the admitted *subject* of
later structural, local, generic, global, and practical identifiability
analysis. Its job was to make an inverse problem closed, replayable, and
impossible to silently widen; it did not itself prove an identifiability
theorem. The authoritative multi-case schema above now owns that role.

Ambition tags: closed law/experiment schema and canonical identities [S];
downstream symbolic, numerical, algebraic, and sheaf-theoretic evidence [F/M]
remains external and is carried only through explicit receipts.

### Historical prototype types and semantics

- `MaterialModelBinding::from_cards` binds the complete immutable
  `MaterialCard`, its exact member `ConstitutiveModelCard`, the narrow canonical
  parameter block, law/state/card schema versions, parameter roster, and an
  explicitly supplied constitutive-graph content binding.
- `ParameterSpec` separates canonical physical quantity/domain/prior from the
  optimizer `ParameterCoordinate`; it also records semantic owner, population
  scope, target/nuisance/fixed class, and honest structural-observability state.
- `InitialStateBinding`, `SpecimenBinding`, `FrameBinding`, and
  `ProtocolBinding` bind initialization, geometry/process/preparation,
  handedness/orientation, load/environment/time paths, experiment clock, and
  refinement semantics.
- `DataLineage::from_vv` accepts concrete `ExperimentArtifact` and
  `CalibrationSplit` values rather than arbitrary references. It rebinds their
  exact canonical hashes, raw manifest/source/custody identities,
  preregistration, blind commitment, partition counts, parser/preprocessing,
  split-grouping identity, and V&V schema version. Only preregistered
  calibration row IDs and their immutable source identities are exposed for
  estimation; validation and blind row identities remain capability-restricted.
- `ObservationSpec` binds the measured quantity, frame, model graph node/port,
  observation operator and version, aggregation, sensor/channel model,
  calibration certificate, transfer/filter/support, clock/delay/anti-aliasing,
  marginal noise, missingness, saturation, protocol/refinement, and exact raw
  calibration rows.
- `ObservationPath` states where a parameter can enter the observation
  distribution (mean, variance, covariance, censoring, missingness, or hidden
  state). Its status distinguishes declared connectivity, symbolic nonzero,
  numerical witness, proven zero, and unresolved paths.
- `NoiseDependence` represents the dimensionless correlation `R` in
  `Sigma = D R D^T`; each channel's dimensional marginal scale remains in its
  `NoiseModel`. Channel order is explicit on input and canonicalized on
  admission.
- `GaugeClass` binds group-action, quotient, slice/inverse, stabilizer-strata,
  and evidence artifacts without conflating a gauge declaration with a proof.
- `DiscrepancyModel` keeps `NoModel`, evidence-backed `Zero`, and `Modeled`
  discrepancy distinct. Every observation needs exactly one row.
- `IdentifiabilityEvidence` has five independent fields—structural, local,
  generic, global, and practical. There is deliberately no ordering or
  automatic promotion among them.
- `IdentifiabilityStudySpec` is the opaque, canonicalized closed graph.
  `AdmittedIdentifiabilityStudy` retains two complete `StudyIdentityReceipt`
  preimages:
  - `StudySpecId` binds every field including the header artifact ID and exact
    optimizer coordinates; it is the replay identity.
  - `PhysicalStudyId` omits only the wire-only header ID and an already
    validated built-in bijective coordinate chart. Priors live in canonical
    physical coordinates, so a prior/support/version change still moves this
    identity.

### Identifiability-schema invariants

1. **Exact model membership:** a model card must be an exact content member of
   the material card. Every constitutive-owned parameter must match the model
   roster and full `QuantitySpec` dimensions; undeclared constitutive roles
   refuse.
2. **Version closure:** schema, V&V, material-card, law, state, experiment
   protocol, refinement, observation-operator, sensor, parser, prior, and
   discrepancy versions are identity-bearing. The versions duplicated in the
   Five-Explicits header must agree exactly; no migration/default is inferred.
3. **State closure:** zero state is accepted only for a model card that permits
   it. A model requiring declared state needs a nonzero exact artifact at the
   exact state-schema version.
4. **Coordinate closure:** physical and optimizer domains are finite. Identity,
   affine, and positive-log built-ins must be nonsingular, dimensionally valid,
   and map the full optimizer interval exactly onto the canonical physical
   interval. Estimated domains are nondegenerate. Signed zero is canonicalized.
5. **Observation closure:** every target/nuisance candidate has at least one
   non-proven-zero declared path, or is explicitly retained as unidentifiable
   with a nonzero witness. This is a connectivity invariant only.
6. **Nuisance and data closure:** each nuisance binds the exact admitted split;
   observation operators may consume only the split's calibration rows.
   Validation and blind rows cannot enter pre-release estimation.
7. **Unit-safe likelihood closure:** derivative `QuantitySpec` dimensions equal
   output divided by input; correlation is finite, positive semidefinite,
   exactly unit-diagonal, and covers the exact observation-channel set.
8. **Gauge closure:** gauge members are distinct, known, nonfixed parameters;
   continuous dimension is bounded by member count; action/quotient/slice and
   stabilizer artifacts are all nonzero.
9. **Explicit discrepancy closure:** missing discrepancy is never interpreted
   as zero. Modeled discrepancy requires versioned content, confounding policy,
   and a separate evidence state.
10. **Canonical transport:** exact transport is length-framed, bounded to 4 MiB,
    domain-separated, stale/future-version rejecting, duplicate rejecting, and
    trailing-byte rejecting. Decode re-runs all admission rules and then requires
    byte-for-byte re-encoding, so noncanonical order and signed-zero aliases do
    not mint alternative preimages.
11. **Identity adjudication:** receipts retain the complete canonical bytes,
    schema version, and bounded row count rather than discarding the preimage
    after hashing.

### Error, determinism, cancellation, and unsafe boundaries

`IdentifiabilityError` names malformed text/numerics, zero identities,
cardinality/duplicate/reference failures, exact version/state/nuisance/gauge/
covariance failures, V&V/material refusals, and byte offsets for canonical
transport failures. Admission never returns a partial receipt.

Schema validation and canonical encoding are deterministic D0 programs over
bounded vectors and `BTreeMap`/`BTreeSet` order. They perform no parallel work,
I/O, waiting, or unbounded iteration and therefore satisfy P7 cancellation by
boundedness. The later symbolic/numerical identifiability analyzers are not
covered by this exemption: they must run as cancellable asupersync tile programs
and emit evidence artifacts only after drain/finalize. This module contains
zero `unsafe` and has no feature flags.

### Archived compile-disabled test inventory

`tests/identifiability.rs` is guarded by `#![cfg(any())]` and therefore does not
compile or run in the current conformance gate. It is retained only as design
archaeology; the list below records what that historical prototype suite was
intended to cover. None of these archived cases is evidence for the
authoritative multi-case schema unless an enabled current test above proves the
same obligation:

- G0 exact decode/re-admit fixed points and retained exact/physical preimages;
- G0 canonical invariance under parameter, observation, path, discrepancy, and
  covariance channel-order permutations;
- G3 affine-coordinate quotient invariance while replay identity moves;
- physical prior-version identity sensitivity;
- disconnected estimates versus explicit witnessed unidentifiability;
- nuisance split mismatch and validation/blind-row leakage refusal;
- independent state/protocol/refinement/clock version refusals;
- correlation normalization and exact channel-set closure;
- valid and dangling gauge declarations;
- missing/no-model/zero discrepancy separation;
- independent structural/local/generic/global/practical evidence fields;
- stale/future schema, trailing bytes, truncation, count bombs, and
  noncanonical unit ordering; and
- signed-zero canonicalization plus NaN, reversed-domain, singular-affine, and
  invalid-log refusal.

### Identifiability no-claim boundaries

- **Schema admission is not an identifiability result.** A connected path is not
  a nonzero sensitivity; a numerical witness is not structural rank; local is
  not global; generic is not uniform; practical is not structural. Evidence
  status changes only when a separately admitted analyzer receipt says so.
- **Hashes bind bytes but do not authenticate issuers.** Neither a content hash,
  the current V&V `authenticated` metadata, nor a boolean/current calibration
  flag proves laboratory authority. Capability/issuer trust and signature
  verification remain external admission policy.
- **The constitutive graph binding is caller supplied.** `ConstitutiveGraph`
  does not yet own a canonical authority-grade semantic identity, so this module
  cannot recompute or authenticate that digest. Downstream theorem promotion
  must wait for a graph-owned identity/cross-check or state this limitation.
- **Physical quotienting is intentionally narrow.** It covers only the built-in
  coordinate charts validated here. It does not quotient gauge orbits, custom
  transforms, prior pushforwards, discretizations, frames, units, likelihoods,
  or discrepancy families. Stronger equivalence may absolutely be proved by
  future theorem-producing machinery, but it must mint a new evidence-bearing
  identity/receipt rather than silently broadening this retired prototype
  identity.
- **Gauge artifacts are declarations until proved.** A nonzero action/quotient/
  slice hash does not establish freeness, properness, orbit completeness,
  singular-stratum coverage, or a global slice theorem.
- **The historical correlation surface is cross-channel only.** Temporal
  kernels, parameter-
  dependent covariance, richer censoring/dropout likelihoods, repeated-trial
  hierarchies, derived-parameter graphs, bounded-logit/simplex/manifold charts,
  and, in this historical prototype, blind-release consumption belong to
  follow-on versioned schemas. The authoritative multi-case successor above now
  enforces blind-release admission without retroactively changing this retired
  prototype identity.
- **No physical-law validation is inferred.** Exact cards, protocols, sensor
  models, discrepancy artifacts, and evidence receipts can still encode a
  scientifically wrong claim. Their independent Gauntlet and authority checks
  remain mandatory.

## ConstitutiveGraph and law-node protocol (bead kagp)

Matter is a typed constitutive graph, not a bag of scalars. `graph.rs`
owns the seven-role decomposition as an executable protocol:

- `NodeRole` — the seven roles; TopologyBalance and BulkStorage are
  DECLARABLE but execution-refused (fs-feec/fs-rep-mesh own them); bulk
  transport, reversible coupling, interface, reaction/source, and
  internal memory are executable.
- `NodeDeclaration` / `LawNode` — every node declares ports (name +
  `fs_qty::Dims` + `TimeParity`), state slots + schema version, a
  calibration `ValidityDomain`, a differentiability class, an
  `EnergyBehavior` (including the EXPLICIT `Empirical` no-claim), and
  whether it claims a consistent tangent and/or a free energy ψ.
  `admit_node` refuses incomplete declarations and probes every claim
  (tangent shape, ψ presence for storage-claiming nodes), naming node,
  law, and failed obligation in each typed `GraphError`.
- Thermodynamic gates (test/audit surface): `check_consistent_tangent`
  (analytic vs central-FD, per entry); `check_free_energy_consistency`
  (outputs are the conjugate forces ∂ψ/∂inputs AND the tangent — the
  Hessian of ψ — is symmetric: Maxwell reciprocity);
  `check_psd_symmetric_part` (second law for force→flux blocks via
  Sylvester on the symmetric part); `check_onsager_casimir`
  (`L[i][j] = εᵢεⱼ L[j][i]` from declared port parities: even–even
  symmetric, mixed-parity antisymmetric).
- `LawRegistry` — implementations keyed by the immutable fs-matdb
  `(LawId, LawVersion)`; instantiation validates the card, checks the
  built node's identity and state-schema agreement, and admits it.
  fs-material CONSUMES card metadata, never redefines it.
- `AggregateStateSchema` — the runtime-state codec when laws coexist:
  exact round trip; version (layout-sensitive FNV fold), length, and
  count drift all refuse.
- `ConstitutiveGraph` — admitted nodes composed by typed edges (dims
  must match EXACTLY; one driver per input port), executed as ONE
  deterministic single pass in topological order (insertion-order tie
  breaks). Cycles refuse: implicit coupling belongs to the solver loop
  wrapped around the graph, never a hidden fixed point inside it.
  Execution audits declared-dissipative nodes for non-negative reported
  rates and totals the dissipation.

`ConstitutiveGraph::linearize` uses this same value evaluator and propagates
directional derivatives through every declared consistent tangent, without
assembling a global Jacobian or advancing committed history during a probe.
`execute_cancellable` and `linearize` poll before and after each law callback;
individual law work must remain bounded. Arity, finite input/output/state,
finite tangent actions and dissipation are admitted before returning any
candidate. External directions cannot drive an internally connected input.
The four `graph_linearized` regressions compare a composed history-dependent
chain against analytic and independent central-difference derivatives, and
exercise sparse directions, cancellation and malformed-law refusal.

GRAPH NO-CLAIMS: single-pass execution is not equilibrium; the
dissipation audit checks REPORTED rates (a law that misreports is
caught only by its own gates/fixtures); free-energy and reciprocity
gates run at caller-chosen points and prove nothing globally;
objectivity/frame-indifference remains per-law fixture scope, not a
graph-level proof.
- Moist air transport uses the Wilke viscosity mixture rule and the
  Wassiljewa–Mason–Saxena conductivity rule with viscosity-based mixing
  factors. Ideal molar heat capacities determine mixture Cp/Cv; Prandtl
  number is recomputed from the resulting Cp, viscosity and conductivity.
  These dilute-gas rules, vapor Sutherland extrapolation and Eucken
  conductivity are engineering approximations without a measured error bound.
  The existing x_w = 0.15 ceiling remains; the saturation ENHANCEMENT
  FACTOR (~0.5% at 1 atm) is neglected; RH > 1
  is refused (a condensing state the ideal-gas mixture cannot
  represent — the phase no-claim made structural at this constructor's
  input; the plain try_new window no-claim stands). CO2 fraction is
  not modeled (the USSA dry-air M already carries the standard CO2).

### Material-card J2 plasticity

`plastic::source::resolve_j2_plasticity` resolves Young modulus, Poisson ratio,
yield stress and an explicitly selected hardening datum at one common physical
state. `ResolvedJ2Plasticity::try_from_material` retains an existing complete
multiphysics bundle, including additional properties and original source receipts.
Its `law()` supplies the existing J2 radial return and consistent tangent directly
to native small-strain solvers, including `fs-solid::SmallStrainTetProblem`.
Scalar, curve and sampled-surface coefficients may depend on typed fixed
pressure, composition, moisture or concentration, alongside fixed absolute
temperature. The shared external-condition admission rejects mechanical,
frequency and rate coordinates that would require an evolving physical owner.

The original hardening bases are the plastic-strain modulus H or a
monotonic uniaxial total-strain tangent Et. For the latter, H = Et E/(E-Et) is
derived after resolving both coefficients at the same state; 0 <= Et < E is
required. Perfect plasticity needs an explicit zero source. The original Et
claim and uncertainty remain intact, and no derived measured-H claim is created.
Source, quantity, domain, parameter and representability errors refuse.
Coefficients remain frozen within one load path; temperature re-resolution
starts a different law, not an energy-accounted transfer of plastic history.
The declared strain bound is limited to 0.05 and remains a caller obligation.
Creep, rate effects, softening, phase changes and finite strain are outside this law.

### Source-backed mixed J2 hardening and cyclic reversal

`J2Plasticity::new_mixed` adds an explicitly supplied nonnegative Prager
modulus `H_k` to the original associative return map. The backstress is
`a=(2/3) H_k dev(eps_p)`, reconstructed from the same retained plastic strain;
there is no second history vector or reset on load reversal. The yield function
is `q(stress-a) - sigma_y - H_iso alpha`. Its radial return and algorithmic
tangent use the same effective trial direction and denominator
`3 mu + H_iso + H_k`. `yield_function_with_state` accepts physical stress and
subtracts this backstress; the original `yield_function(stress, alpha)` is the
centered effective-stress diagnostic. Existing `new` and `H_k=0` retain the
original isotropic numerical path and model metadata.

`HardeningPropertyBasis::MixedPlasticStrainModuli` requires separate
`isotropic_hardening_modulus` and `kinematic_hardening_modulus` source claims
in Pa, with the convention above. Each retains its source identity and actual
state query. No monotonic tangent measurement identifies or silently supplies
the split. The original isotropic choices remain isotropic even if a shared
bundle contains a kinematic datum. The exact requirements are returned as a
vector because the mixed reduction has five mechanical roles.

Recoverable free energy includes `H_k dev(eps_p):dev(eps_p)/3`, reported as
`kinematic_hardening_j_m3` separately from isotropic hardening storage.
The corresponding quadratic increment belongs to backward-Euler numerical
loss. Physical dissipation remains `sigma_y delta_alpha`; selected native
linear overstress adds `eta delta_alpha^2/dt` once and includes `H_k` in the
same return denominator. Neither recoverable backstress storage nor numerical
loss becomes plastic heat.

The original `ThermoplasticSource` accepts this same explicit split, but `H_k`
must be constant at every source-line knot across its whole caloric chart.
It then adds no omitted temperature derivative to mechanical entropy. Existing
affine E/H_iso, thermal expansion, yield softening, actual finite caloric
feedback and sourced overstress retain their original owners and accounts.
Temperature-dependent H_k refuses, even if its initial and endpoint values
coincide. Native 3-D small-strain solids consume the resolved law unchanged.

Authored regressions compare an independent reverse-yield limit, distinct
cyclic response at identical monotonic tangent, nonproportional consistent
tangents, separate recoverable/physical/numerical energy, the zero-H_k limit,
source refusal, a native P1 cyclic coupon and finite caloric heating.
They do not establish executed test results or empirical cyclic calibration.
This advances MR22/MR62 without claiming nonlinear cyclic hardening, fatigue
counting, life prediction, finite-strain plasticity or broad bead completion.

### Isothermal liquid momentum and free-surface source bundle

`liquid::surface::resolve_liquid_surface_state` requires sourced density and an
explicit dynamic/kinematic viscosity basis. `LiquidCapillarity::Sourced`
additionally requires positive surface tension at the same complete query;
`Neglected` is an explicit approximation, not a zero-valued material claim.
An existing shared material bundle can be adopted without losing extra receipts.
Kinematic conversion uses the same resolved density. Source/domain/unit/selection
errors and unsupported dynamic coefficients refuse the complete reduction.

The actual consumer is `fs-lbm::material::free_surface`. No heat-capacity or
conductivity tuple is needed for this isothermal momentum model. Source conditions
and surface context remain frozen; the adapter supplies no phase diagram, ambient
gas, wetting or purity default, runtime material-state transfer, or evolving
surface-energy account. Scalar constants and supported temperature coefficients
do not authorize frequency- or rate-dependent Newtonian flow.

`gas_field::steady::SteadyAxialGas::solve` lowers its source packs, axial volumes,
supplies and walls through `to_network` into the shared
`network::SteadyGasNetwork`. That owner also admits branched and inlet-fed
recirculating graphs, with at most 64 cells and 256 interior/boundary connections.
The existing fs-thermal `Advected` observer owns both donor-cell enthalpy and
stationary species transport; the original thermal timeline owns wall heat,
combined CFL and energy updates. Cyclic species balance uses the existing
fs-la prepared LU with bounded work and a numerical conditioning allowance;
unanchored closed loops refuse instead of inventing a composition or leakage.
`GasEnthalpyLaw` derives h(T), T(h,p), and a conservative thermal slope bound from
the original NASA-9 source. Caloric laws may differ by stationary composition
only after exact full species/fraction and common enthalpy-reference checks.
Every inlet carries its own source-mixture enthalpy. Initial EOS cell masses
are pseudo-time relaxation weights, not a physical variable-density history.
Only converged, phase-admitted common-pressure states are published, retaining
their thermochemical and material transport sources for existing flow, thermal
and acoustic consumers. Each final power gate bounds every donor/outlet caloric
closure flux separately so inverse errors cannot conceal an inlet energy error.
The axial adapter additionally checks its reconstructed inlet face, then emits
the same `AxialGasField` used by existing ducts and reeds. Balances, work, flow
and bulk-cell low-Mach admission survive this lowering. Downstream faces carry
first-order mixed-cell states; field interpolation and acoustic sampling still
need independent spatial refinement. The scope excludes pressure-drop and
branch-throat velocity prediction, transient variable-density flow, inferred
wall conductances, diffusion, reaction and condensation. Source, phase, work,
convergence and publication failures return no partial solution. This does not
close the complete M0.1-W3 qualification.

### Equilibrium substrate wetting on the liquid surface owner

`liquid::surface::resolve_liquid_wetting_state` resolves only
`equilibrium-contact-angle` [radians through liquid] from an ordered
(substrate, liquid) interface. The original interface identity, preparation
frames, source receipt, history and named environment survive. The caller
assigns that environment; this response does not resolve an ambient gas EOS.
The angle lies strictly inside (0, pi), and the same frozen-coefficient rules
as the momentum/surface bundle refuse dynamically supported data.
Missing equilibrium data never falls back to advancing/receding endpoints.
`require_liquid_state` checks the declared liquid phase/process/revision and
complete query/quantity descriptors; SurfaceSpec does not carry a bulk-card
content hash. Liquid/gas sigma remains solely in `ResolvedLiquidSurfaceState`.
The sparse free-surface consumer assigns this response to one planar domain
wall. No hysteresis evolution, moving-contact-line regularization, Marangoni
traction or measured-accuracy claim follows from that assignment.

### Incremental J2 work and physical heat

`plastic::energy` audits the original radial-return stress and history update
under the explicitly energetic interpretation of linear isotropic hardening:
elastic energy is `0.5 e:C:e`, hardening storage is `0.5 H alpha^2`, and
physical dissipation is `sigma_y delta_alpha`. Endpoint stress work also
contains the backward-Euler loss
`0.5 delta_e:C:delta_e + 0.5 H (delta_alpha)^2`; that numerical term is
reported separately and is never a physical heat source.

The mixed extension above adds its separately reported Prager storage and
quadratic increment to these accounts. Callers must retain the same elastic
coefficients and both hardening moduli across the increment.
Changing only the yield offset can release stored energy at fixed total strain.
The owner checks coefficients, tensor strain bounds, finite deviatoric plastic
history and its accumulated-strain bound, endpoint yield consistency, and an
independent absolute work residual. Proposals return immutable native stress,
history and energy accounts. These gates do not establish a measured heat
partition for an arbitrary empirical hardening curve or supply creep, phase
change, thermal expansion, finite-strain work, or parameter-work accounting.

The regression source covers independently derived pure-shear work and storage,
elastic unloading without physical heat, reversal with accumulated history,
yield reduction at fixed strain, and invalid-history refusal.

### Energetic J2 caloric feedback

`plastic::thermal::CaloricPoint` binds a virgin native J2 history to a
source-resolved, strictly increasing solid caloric chart. Yield stress and heat
capacity may vary with absolute temperature; Young's modulus E and native
energetic isotropic hardening H may be globally affine in temperature. Poisson
ratio, any selected Prager modulus, density and reference thermal geometry
remain fixed. All property roles share
one material card and source context. Each selected claim supplies one
continuously supported canonical temperature line over the full chart. Every
crossed curve knot and surface-cell edge participates in affine-slope, fixed
coefficient and positive yield/Cp admission; equal endpoints cannot hide an
interior derivative jump. Mechanically bounded strain, stress or rate
coefficients cannot be frozen while their histories evolve. Renamed temperature
coordinates require explicit absolute-temperature semantics. Variable energetic
H requires the explicit native plastic-strain-modulus basis. Uniaxial tangent
conversion retains its constant-E/Et path; variable E with exactly zero Et is
also admitted, but a generally rational derived H(T) is not treated as affine.

This reduction declares energetic isotropic hardening: `H alpha^2/2` is stored
energy; optional Prager storage is separate and also recoverable.
`sigma_y(T1) delta_alpha` is the physical plastic heat supplied to the
existing `UniformEnthalpyStep` transport owner; selected expansion adds the
signed reversible exchange below. Elastic storage, hardening storage,
boundary heat and backward-Euler numerical loss remain distinct. Every thermal
proposal starts from the same accepted mechanical/caloric history; the actual
returned temperature re-evaluates the native J2 increment before separate heat
and combined-work gates accept the candidate. A rejected or cancelled proposal
changes neither history nor time. Optional sourced isotropic thermal
expansion is supported as described below. The model excludes nonlinear or
piecewise-affine E/H slopes, variable Poisson ratio, phase change and
finite-deformation volume work; source receipts and residual gates do not
establish empirical accuracy.

### Sourced thermoelastic J2 expansion

`SourceQuery::thermal_expansion` explicitly selects either zero eigenstrain or
the canonical instantaneous linear expansion coefficient on the same card.
The selected finite scalar, curve or fixed-context surface supplies its complete
pinned temperature line, including every source knot and crossed surface cell.
Negative coefficients and zero crossings are admitted. The existing shared
piecewise-linear primitive computes `f(T)=integral(alpha dT)` without resampling;
`ThermalJ2` passes `eps-f(T)I` to the original J2 stress, consistent tangent and
history update. The reference temperature stays fixed at specimen birth.
Renamed temperature coordinates retain their explicit absolute-temperature
semantics and the original fixed nonthermal conditions.

With fixed Poisson ratio and globally affine E(T)/H(T), free energy remains
`psi=e:C(T):e/2+H(T) p²/2`. Mechanical entropy is
`s=alpha(T) tr(sigma)-(E'/E) psi_elastic-H' p²/2`; internal energy adds `T s`.
The native J2 energy account still evaluates both histories with the new law.
`coefficient_storage_change_j_m3` exposes the missing frozen-history coefficient
change `Cpsi=psi(T_new,e_old,p_old)-psi(T_old,e_old,p_old)`; the increment's
`stored_change_j_m3()` includes it exactly once. Signed reversible caloric
transfer is `tr(sigma_new) integral_Told^Tnew(alpha dT)-Cpsi-delta(T s)`, using
the direct step integral and actual pinned endpoint coefficients. With fixed
E/H and constant alpha, the original `-T_old alpha delta(tr(sigma))` arithmetic
is retained. Physical plastic heat is still `sigma_y delta_p`; numerical loss
remains separate. This closes the actual two-temperature endpoint work balance
without replacing the native return map, history, tangent or time refinement.

Source Cp and actual caloric segment slopes must stay positive after
thermomechanical coupling. Let A=max|alpha|, A'=max|alpha'| and K_max be the
maximum bulk modulus over the complete source span. The sufficient gate is
`rho c_min > T_max [9K_max A² + 3sqrt(3)K_max strain_limit A' + 6sqrt(3)|K'| strain_limit A]`.
The native energy owner bounds the actual elastic tensor by the tensor-Voigt
norm (diagonal squares plus twice shear squares), giving
`|tr(sigma)| <= 3sqrt(3)K_max strain_limit` and
`|tr(C' e)| <= 3sqrt(3)|K'| strain_limit` without assuming exact zero plastic
trace. These bound the signed terms in
`c_eff=rho c_caloric+T[alpha' tr(sigma)+2alpha tr(C' e)-9K alpha²]`, including
both sides of expansion-curve knots. Affine E/H have no second derivatives;
no unsupported plastic-history bound is needed for H. Fixed E/H retain the
previous capacity arithmetic. Failed admission publishes no material or
thermal candidate. Reference mass, geometry and Poisson ratio remain fixed;
phase change, derivative jumps in E/H and finite volume work remain outside
this model.

Focused synthetic tests cover expansion-knot integration and cooling, signed
and zero expansion, affine-modulus curves and typed source surfaces, independent
free-energy/entropy identities, data-only changes to adiabatic temperature and
native 3D loaded shape/plasticity/heat, continuation with coefficient work, and
atomic effective-capacity refusal. These fixtures are not measured material
validation; native execution of the new modulus regressions remains pending
in the configured build environment.

### Structural source bindings for thermal transport

`ThermoplasticSource::thermal_material_point` exposes the exact original card,
complete fixed context and declared temperature axis. `pin_thermal_transport`
selects a nonnegative static property at the actual accepted caloric state,
applies the same mechanical-coordinate/frequency guard as the six native
thermoplastic roles, and returns its original claim pin. The consuming
transport owner admits the property quantity. No material state, caloric law,
plastic history or transport formula is duplicated. `ThermalMaterialPoint`
admits the same typed fixed external conditions for conductivity source axes.
Its explicitly coupled additional axes keep their original physical owners;
an additional temperature cannot be frozen or introduced as a nonthermal axis.


### Declared mushy-phase volume and mobility

`phase::mushy::MushyPhaseLaw` binds the original equilibrium h chart to explicit
constant constituent densities, a permeability scale and a positive minimum
liquid-volume fraction. `PhaseVolumeFractions` converts mass fractions using
constituent specific volumes and checks the chart's bulk density. The declared
Carman–Kozeny closure supplies inverse permeability from liquid volume fraction;
it is neither a universal melting law nor a measured-data assertion. The source
state, original chart and model identity remain available. Fully liquid material
has zero added drag; states below the support floor refuse rather than inventing
solid locking. Equal constant phase densities are separately admitted for
fixed-volume consumers. Separate constituent velocities, elastic skeleton
response, pressure/composition evolution and moving boundaries remain outside
this closure.

### Source surfaces in hereditary thermal materials

`thermovisco::source::ResolvedThermoMaxwell` lowers the canonical typed
temperature line from scalar, curve and sampled-surface claims. Every crossed
source-cell edge, original role key, selected claim and fixed material context
is retained. Reference spring and relaxation coefficients must remain exactly
constant across the complete line; the declared modulus/time/free-strain laws
provide the supported temperature response. Source gaps and discrete-only
intervals refuse rather than supplying interpolated coverage.

Hyperelastic, hereditary and J2 source bindings share a single
fixed-external-coordinate policy. Typed pressure, composition, moisture and
concentration may describe fixed external conditions. Evolving mechanical
coordinates and a second temperature coordinate cannot be frozen into a thermal
law. Frequency-dependent values refuse; thermovisco frequency applicability
alone must retain its static/relaxed endpoint. Caloric feedback admits a globally
affine common modulus multiplier and free thermal strain across the complete
retained knot union. It uses mechanical entropy `s=alpha*sigma-S_prime*Xi`
from the free-energy derivative and accounts for
internal energy separately from spring free energy. Physical dashpot loss and
signed reversible exchange, including the native thermal-eigenstrain work,
enter the original caloric transport owner; numerical loss remains separate.
Expansion requires positive total heat capacity over the complete chart and
admitted strain/history range, including instantaneous branch stiffness and
the mixed expansion/modulus derivative term. The native branch-strain envelope
is enforced on expansion-bearing endpoints. Changing reference density, phase,
or either affine derivative still requires additional terms and refuses.

### Shared phase-dependent conductivity

`phase::conductivity::LiquidMassFractionConductivity` owns the declared positive
multiplier linear in liquid mass fraction. It preserves the existing FEM law's
nearer-endpoint arithmetic and derivative. This explicit phenomenological model
scales a separately admitted base k(T); it does not infer a volume-mixture rule,
solid/liquid material measurements, or additional source authority. The original
fs-conduction facade and live material-flow thermal coupling consume this same
law. The immutable equilibrium h chart remains the phase-state owner.

### Shared linear magnetic permeability admission

`field_storage::linear_magnetic_permeability` exposes the existing absolute
scalar permeability admission used by `ResolvedInductor`. The linear planar
`fs-em` field owner calls this same function, so circuit storage and spatial
fields share dimensions, source support and linear lossless model restrictions.
Fixed temperature and explicitly typed pressure/composition/moisture/concentration
source snapshots are allowed, including source curves and surfaces. Both source
axes and applicability domains use the public
`ResolvedScalarProperty::is_static_at_fixed_external_conditions` classification,
backed by the same fixed-external policy as mechanical and thermal laws. Field/frequency/strain/stress dependence or bounds still require
a richer owner and refuse here. The existing circuit rebind path owns explicit
parameter work; fixed-storage winding and moving-gap owners retain their stricter
source-evolution guards. Focused field-storage tests cover pressure/composition
surfaces, typed external curves, original claim retention and dynamic-role refusal.
This adds no independent constitutive formula or magnetic heating law.

## Reversible thermal field storage

field_storage::thermal supplies source-bound reversible thermal field response for fixed-geometry linear C(T)=epsilon(T) A/d and L(T)=mu(T) N^2 A/l. The exact original material card, selected claims, typed absolute temperature and fixed external context remain attached to each refreshed native element. The existing continuous source-path owner checks the full requested interval; the shared affine-temperature derivative admission is reused by thermoplastic and field-storage consumers. Globally affine positive coefficient paths (including constant and fixed-context surface slices) supply field Helmholtz energy psi=z^2/(2K), entropy s=psi Kprime/K, internal energy U=psi+T s, and the nonpositive fixed-coordinate heat-capacity contribution -2T psi(Kprime/K)^2. A conservative full-interval bound supports the coupled owner admission of positive total caloric plus field capacity under an enforced coordinate limit. This material law does not advance a circuit or thermal state. Coupled owners must account reversible field-caloric exchange separately from Joule heat and add U, not psi, to total internal energy. Derivative jumps, field/frequency-dependent coefficients, second thermal source axes, hysteresis, phase transitions, moving geometry and changing nonthermal context refuse or require another owner. Synthetic G1 tests check potential derivatives, electric/magnetic signs, native geometry, source retention, surface paths, coefficient swaps and validity/refusal boundaries; they are not experimental material validation.

## Spatial affine magnetic thermal material

`field_storage::thermal::magnetic::ThermalMagneticMaterial` pins a complete original material snapshot, its immutable card and selected claims, an explicit absolute-temperature axis, fixed remaining context, and a positive globally affine permeability path. It shares source-path admission, Helmholtz derivatives and capacity arithmetic with lumped C(T)/L(T). For actual cell volume V and gauge-invariant |B| it returns psi=V B²/(2mu), entropy=psi mu_prime/mu, U=psi+T entropy and fixed-B heat capacity=−2T psi(mu_prime/mu)². The full-range conservative capacity bound requires the caller to enforce its B envelope and add all field contributions to the caloric capacity. This is a linear reversible fixed-geometry rung; kinked/nonaffine paths, hysteresis, saturation, transition physics and experimental validation are not inferred.

## Source-integrated directional hygroscopic strain

`state_point::integrate_orthotropic_hygroscopic_expansion` integrates three
explicit dimensionless `linear_hygroscopic_expansion_coefficient_{1,2,3}`
claims with respect to dry-basis moisture (kg water/kg dry solid). Coefficients
are engineering-strain derivatives relative to the original stress-free length.
The typed `hygroscopic_strain_reference_moisture_content` claim must equal the
reference query throughout the admitted path. Temperature is explicitly absolute
and fixed; all other query coordinates and schemas remain identical. Selected
claims are pinned and the original material-line owner admits continuous support
and exact piecewise-linear integration, including surfaces and drying. Missing
axes/directions, ambiguous sources, datum changes and domain gaps refuse.
This is a reversible additive small-strain constitutive path, not a finite-strain,
hysteretic, sorption-energy or stress-assisted transport law. Source and datum
identities accompany the complete endpoint evaluations. Tests use authored
analytic coefficients, not experimental wood validation.

## Source-backed permanent-magnet recoil

`field_storage::permanent::ResolvedRecoilMaterial` binds the reversible isotropic law `B = mu_recoil H + Br n` at one complete fixed material state. Required source properties are `recoil_magnetic_permeability` [H/m], `remanent_flux_density` [T], signed `recoil_minimum_parallel_field_strength` and `recoil_maximum_parallel_field_strength` [A/m], and `recoil_maximum_transverse_field_strength` [A/m]. The parallel interval must contain zero; all values retain their original units, selected claims, card identity, query coordinates and receipts. Generic permeability, coercivity and energy-product-derived model coefficients cannot substitute for recoil data. All coefficients and limits must be static at the selected external conditions; field/frequency dependence cannot be frozen into this linear law.

`at_query_point` refreshes the complete source-pinned bundle, including Br, slope and operating limits. It creates a new fixed-state snapshot; it does not supply live thermal parameter work or irreversible demagnetization. Specimen orientation belongs to the field owner. Computed parallel/transverse H and signed margins are checked against every source limit; those margins are model-admission estimates, not certified bounds on the numerical field. Remanence-only catalog entries remain insufficient input.

## Source-backed reversible nonlinear magnetics

`field_storage::nonlinear::ReversibleMagneticMaterial` consumes the explicit
`reversible_magnetic_field_strength` source property [A/m] along a named,
dimensioned flux-density coordinate [T]. The original matdb line query pins
one source across the whole interval from zero to the declared maximum,
including fixed-context slices of two-axis source surfaces. Every breakpoint
retains its canonical answer and receipt. Field-dependent secant permeability,
coercivity and material names are not alternative ingress paths.

Admission requires H(0)=0 and finite strictly positive segment slopes. Other
source axes and applicability bounds must be fixed external thermodynamic
conditions. Discrete-only spans, gaps, source switching, nonmonotone response
and extrapolation refuse. Runtime H evaluates the compiled canonical source
path without rebuilding material queries in each field/solver iteration.
Binary search selects the original segment in O(log knots), with no allocation
or independent permeability fit. The same piecewise-affine path supplies H/B,
dH/dB, exact segment-integrated energy `psi(B)=integral_0^B H(b) db` and
coenergy `integral_0^H B(h) dh`. Both integrals are accumulated directly so
small positive coenergy survives cancellation in `B H-psi`.
`flux_density_at_field_strength` inverts this same curve and its original
supported range; it does not introduce a separately interpolated inverse law.
The tangent at zero uses the first segment; an interior knot uses its right
derivative. Bregman loss visits only the original segments crossed by the two
field magnitudes. `bregman_density_3d_j_m3` extends the same original-slope
radial/angular decomposition to spatial fields, including rotation and crossed
saturation segments, without subtracting large stored energies. The planar
entry point keeps its reduction order. This loss is numerical dissipation,
never hysteretic or Joule heat. Six direct convexity tests cover both dimensions,
independent energy identities, rotation covariance and tiny angular losses.
These are floating-point model evaluations, not uncertainty
propagation or enclosures. Focused material regressions compare compiled
curves and fixed-temperature surface slices to canonical source queries, check
inverse and Legendre derivatives, and preserve domain refusals.

The law is reversible, isotropic and without remanence. A hysteresis loop or
initial magnetization branch does not establish that physical model. Static
temperature/material substitution creates a new law; it does not supply
magnetic heat, parameter work, hysteresis or irreversible demagnetization.
G1/G3 regressions in `fs-em/tests/nonlinear_magnetic.rs` exercise analytical
energy/tangents, material-only solved saturation, mixed-material oblique fields,
the linear limit, source-domain refusals and cancellation.

## Reversible longitudinal piezoelectric material

`piezoelectric::LongitudinalPiezoelectric` resolves density, directional
`piezoelectric_compliance_33_at_constant_field` [1/Pa], signed
`piezoelectric_strain_charge_33` [m/V = C/N], and absolute
`piezoelectric_permittivity_33_at_constant_stress` [F/m] from one original
material state. The common local direction is material axis 3, transverse
stress is free, and the electric field is longitudinal. Scalar/curve/surface
queries retain their complete selected claims and fixed external conditions;
dynamic field, stress, strain or frequency support cannot be frozen.

The strain-charge convention is `strain=s_E stress+d E`,
`D=d stress+epsilon_T E`. Coupled stability requires both `s_E>0` and
`epsilon_axial=epsilon_T-d*d/s_E>0`; positive individual diagonal coefficients
alone do not suffice. The single energy is
`w(strain,D)=strain^2/(2s_E)+(D-(d/s_E)strain)^2/(2epsilon_axial)`.
Its derivatives supply both stress and electric field. The reduced
`epsilon_axial` fixes axial strain with transverse stress still free; it is
not asserted to be the fully clamped three-dimensional permittivity.

This is a reversible, fixed-state, small-signal one-dimensional law. It does
not infer off-axis tensors, arbitrary tensor rotation, saturation, depoling,
hysteresis, leakage, dielectric heating or thermal coefficient work. Source
selection preserves evidence, not measured predictive accuracy. G1/G3 tests in
`tests/piezoelectric.rs` independently check strain-charge equations, energy
derivatives, open/short stiffness, signed/zero coupling, SI-sized synthetic
coefficients, and whole-block/source-domain refusal. Actual coupled motion,
charge and native circuit composition belong to `fs-couple::piezoelectric`.

## Explicit liquid source admission for evolving shear rate

`ResolvedLiquidSurfaceState::try_for_shear_rate_evolution` retains the original
liquid card, viscosity basis and selected claims while explicitly identifying
the temperature and physical shear-rate axes that its consumer must evolve.
It introduces no viscosity law, guessed rate or surface-energy model. The
existing frozen source constructor continues to refuse rate-dependent claims;
ordinary frozen flow ingress also refuses the evolving marker. Only an owner
that solves those physical coordinates may use this bundle. The moving-liquid
enthalpy/rheology owner in `fs-lbm` provides that integration. Focused admission
tests cover dynamic/kinematic viscosity bases, source retention and refusal of
unowned axes; caloric response with shear-dependent storage remains unsupported.

## Source-backed planar orthotropic electrical response

`conductor::planar` resolves two independent principal conductivity or
resistivity properties and an explicit fixed material-to-geometry frame.

Scalar and planar tensor snapshots may explicitly retain one solver-owned
liquid-mass-fraction coordinate using `with_liquid_mass_fraction_axis` and
`try_from_phase_material`. It must exist as a dimensional scalar in `[0,1]`.
The existing static source policy applies to every other coordinate. Refresh
preserves the declared axis, source basis and tensor frame; it does not evolve
phase or infer a solid/liquid mixture law. Coupled consumers must bind that
axis to their actual caloric liquid mass fraction. Unowned response dependence
remains refused.

Resistivity is inverted in its principal basis after source evaluation; the
existing tensor rotation computes `sigma = Q diag(sigma1,sigma2) Q^T`.
The complete material snapshot, source basis, frame identities and orientation
remain bound to the result. Missing directions never become scalar defaults.
The sealed common response also admits the existing scalar owner and a mixed
scalar/orthotropic region enum. Both feed the same native `fs-em` field solve
and its `fs-couple` temperature-feedback driver. This is a decoupled physical
plane with symmetric static Ohmic response; Hall/AC effects, inferred normal
transport and changing principal directions require additional physical laws.
Program anchor: MR43 (`frankensim-material-reality-ui2if.7.1`), remains open.

## Finite ideal-gas internal-energy capacity for coupled exchange

`ThermochemicalGasModel::maximum_isochoric_temperature_slope` bounds `dT/du`
over an explicit source-supported temperature bracket using the original NASA-9
coefficient interval bounds with bounded deterministic temperature subdivision
and each species' `cv = cp - R/M`. Every accepted subinterval proves positive
`cv`; an unresolved or unstable source range refuses. The bound is a timestep
aid; the existing source model still owns actual `u(T)` and inversion.

`ThermochemicalGasModel::with_caloric_continuity` applies the explicit
thermochemical reference policy to each retained species while preserving
composition, transport models, and source-pack identities. `GasEnthalpyLaw`
and isochoric capacity bounds may cross NASA-9 region boundaries only for
species carrying that admitted policy. Raw models retain the single-region
requirement. Each intersected source region receives its own outward-rounded
Cp/Cv bound, including the upper chart selected when the bracket ends exactly
at an internal seam. The minimum admitted capacity bounds the complete bracket;
reconciliation alone does not establish caloric stability. A fixed 13-entry
stack limits subdivision to 8,191 interval evaluations per intersected region,
with at most 16 regions per species.

For corrected models, transport enthalpy-reference identity includes the
complete retained species catalog and every explicit caloric policy,
independently of the current mixture fractions. Mixtures of the same corrected
catalog remain compatible across composition changes; raw and corrected
references, or different correction policies, remain distinguishable. Source
thermochemistry, energy inversion, thermal entropy changes, and reaction
thermochemistry therefore use the same declared caloric reference.

`IsochoricGasParcel::model` exposes that same immutable composition/source model
for coupled film queries. No constant heat capacity, enthalpy substitution or
second gas caloric law is introduced. Regressions cover source-supported
multi-region nitrogen heating, inversion, entropy and pressure response,
unchanged raw admission, and refusal of nonpositive Cv at a selected seam.

## Biased incremental longitudinal magnetostriction

`magnetostrictive::LongitudinalMagnetostrictive` resolves density, incremental
axial compliance at fixed H, the signed reciprocal piezomagnetic d33, and
incremental absolute permeability at fixed axial stress from one original
material bundle. The coefficient property names explicitly distinguish these
incremental, direction-3, stress/field conventions from an ordinary magnetic
permeability or a global magnetostriction curve. No grade name selects data.

The original query must contain an explicitly typed, finite nonzero H-bias
coordinate in A/m. This owner holds that bias and every other external condition
fixed. The shared owned-coordinate/fixed-external admission retains typed
temperature, pressure and composition and refuses unowned dynamic fields,
stress, frequency and rate dependence. Bias- or temperature-dependent source
curves select a tangent at initialization; they are not silently refreshed
along the trajectory. A scalar requires `ConstantWithinValidity`.

The reversible perturbation equations are `epsilon=s_H*sigma+d*H` and
`B=d*sigma+mu_T*H`. All four variables are increments about the already biased,
stress-free prepared reference. `s_H>0` and `mu_axial=mu_T-d*d/s_H>0` admit the
whole coupled block. With independent epsilon and B, one energy
`w=epsilon^2/(2*s_H)+(B-(d/s_H)*epsilon)^2/(2*mu_axial)` gives both stress and H.
This is incremental storage, not the absolute energy of the DC-biased device;
bias preparation and its source work remain outside the model. The constitutive
convention follows Deng and Dapino, *Magnetic flux biasing of magnetostrictive
sensors*, Smart Materials and Structures 26 (2017) 055027, Eq. (1),
<https://doi.org/10.1088/1361-665X/aa688b>. No measured constants are imported.

This is a local small-signal, longitudinal, reversible tangent with free
transverse stress. It supplies neither a full tensor, a stress-preload law,
magnetization/domain history, hysteresis, saturation, heating nor parameter work.
The dynamic consumer must declare a strain/field increment envelope; avoiding
zero-bias crossing does not prove that a measured law is linear in that envelope.
Synthetic material tests compare an independent symmetric-block inverse and
energy, work-conjugate derivatives, signed/bias/temperature data changes,
coupled instability and source-domain refusal. Motion and native circuit
composition belong to `fs-couple::magnetostrictive`.

## Source-defined linear J2 overstress and physical creep heat

`plastic::rate::LinearOverstressJ2` extends the original radial-return owner
with the explicitly declared associative Perzyna n=1 law
`alpha_dot = max(q - sigma_y - H alpha, 0)/eta`. The native tensor convention,
history, trial stress and consistent algorithmic tangent remain shared. For an
actual physical interval dt, backward Euler gives
`delta_alpha = max(f_trial,0)/(3 mu + H_iso + H_k + eta/dt)`, where `H_k=0`
for the original isotropic law and q is formed from stress minus backstress.
Positive endpoint overstress
is a legitimate committed state; the consistency gate checks
`q - sigma_y - H alpha - eta delta_alpha/dt`. Every new interval starts from
the accepted internal variables. Native viscosity and duration must be finite,
positive and representable together with the return denominator.

The energetic-hardening interpretation remains explicit:
`psi=e:C:e/2 + H_iso alpha^2/2 + H_k dev(eps_p):dev(eps_p)/3`. Physical dissipation is
`sigma_y delta_alpha + eta delta_alpha^2/dt`; the viscosity is not added to
hardening storage. The original backward-Euler numerical loss is separate and
never deposited as heat. Fixed-strain stress relaxation can therefore warm the
finite thermal store while reducing actual elastic storage. The original
rate-independent law and arithmetic remain the default when no viscosity is
selected.

`ThermoplasticSource::with_linear_overstress` selects
`plastic_overstress_viscosity` [Pa s] from the same original mechanical/caloric
material card and complete fixed external context. This role is distinct from
fluid shear viscosity. The original claim is pinned; its continuous source line
must stay finite and positive across the complete caloric interval. Piecewise
temperature dependence is admitted; unowned mechanical, rate, history and
frequency source coordinates refuse. Every thermal trial and the accepted
endpoint refresh this claim together with the other material coefficients.
Viscosity affects kinetics and dissipation, contributes no stored free energy,
and creates no additional mechanical entropy term. Existing affine E/H and
integrated expansion accounts remain intact.

`CaloricPoint` and the spatial thermal owner supply actual physical dt.
`resolve_expanded_enthalpy_for_step` and `advance_with_energy_for_step` expose
the same interval-aware source law; their duration-free counterparts refuse a
source that explicitly selected overstress. Coupling proposals and thermal
substeps cannot accumulate extra plastic time. This is small-strain linear
overstress with energetic linear hardening, not a nonlinear rate exponent,
finite-strain creep, damage, phase transformation or empirical validation.

Synthetic `tests/plastic_rate.rs` compares independent radial-return/storage/
heat expressions, a fixed-strain analytic stress-relaxation solution with
first-order refinement, all discrete tangent entries and eta/dt scaling.
The spatial consumer tests actual loaded-solid creep, finite-store heating and
endpoint temperature feedback in `fs-couple`. This advances MR22
(`frankensim-material-reality-ui2if.3.4`); its broader material/history scope
remains open.

## Constant dielectric storage admission for native charge relaxation

`field_storage::constant_electric_permittivity` admits the existing absolute
`electric_permittivity` requirement only when its selected source declares a
constant law over its admitted validity domain. A coincidentally unchanged
sample from a temperature or phase curve does not establish this constitutive
contract. Units, positivity, complete source receipt and original query
semantics remain owned by the common material resolver.

The native leaky-dielectric trajectory uses this admission to change
conductivity at the actual thermal state without silently introducing
permittivity parameter work or reversible thermal entropy. It retains the
complete common source bundle; changing material identity, source claims,
coordinate meanings, tensor frame or phase-axis binding on an existing
trajectory is refused. Variable dielectric storage needs its own work/entropy
law and is outside this entry point.

## Composition replacement at a shared heat and chemical energy endpoint

`ClosedGasParcel::replace_model_with_heat` extends the original conservative
composition boundary with one signed heat increment. At unchanged mass and
volume it reconstructs the new source-backed thermochemical model directly at
`U_end = U_start + Q` through the existing density/internal-energy inversion.
It does not require an intermediate adiabatic composition state or an
intermediate heat-only state to lie inside the material support. The original
conserved-energy replacement delegates with zero heat and keeps its existing
meaning.

`GasCompositionStep` now reports supplied heat, actual stored-energy change and
their residual. Formation energies remain part of the original species internal
energies; callers must not add their release again as an independent heat term.
All endpoint EOS, caloric, transport and acoustic properties come from the same
accepted composition, temperature and pressure. Nonfinite or unrepresentable
energy increments, source/domain failures, failed inversion and cancellation
publish no parcel. This primitive supplies a material state transition, not
reaction kinetics, species transport or open-system mass/chemical work.

## Composition, volume and mechanical work at one caloric endpoint

`ClosedGasParcel::replace_model_with_heat_work_volume` reconstructs the new
source-backed material model at unchanged mass, specified final volume and
`U_end = U_start + Q - W`. Positive W is work delivered by the gas. The original
density/internal-energy inversion resolves temperature, pressure, transport and
acoustics together, without an intermediate composition, heat-only or volume-only
state. The existing fixed-volume heat and conserved-energy methods delegate with
zero work and retain their meaning.

`GasCompositionStep` records both volumes, supplied heat, delivered work,
actual storage change and `delta U + W - Q`. The mechanical consumer must
independently close the supplied work against its actual pressure/volume
trajectory; this material primitive supplies no work quadrature or entropy
production law. Formation references remain inside species internal energies.
Nonfinite/lost exchanges, unsupported endpoints, failed energy inversion and
cancellation preserve the original parcel.

G1/G4 regression cases in the original composition module independently check
an ideal-gas isothermal expansion with equal supplied heat and delivered work,
and failed/cancelled moving-volume reconstruction. Parser checks completed;
native Rust execution was unavailable because DSR/RCH are absent and the shared
filesystem is full.


### Reference-aware hyperelastic material interface

`hyper::cards::HyperelasticMaterial` exposes `energy_density`, `first_piola`
and `consistent_tangent` with defaults delegating to the original native law.
A reference-deformation adapter must transform all three consistently and
retain the underlying elastic-domain owner. Spatial consumers can preserve the
complete material binding instead of stripping it to bare coefficients.
The public `state_point::is_fixed_external_coordinate` classifies one coordinate;
it does not replace typed-query, source-axis, validity or path admission.


### Explicit planar dielectric source laws

`field_storage::planar` admits explicit scalar or fixed-principal-axis
in-plane permittivity. Orthotropic sources provide independently selected
`electric_permittivity_1` and `electric_permittivity_2` in absolute F/m,
with the existing named `PlanarElectricalFrame`. The shared storage source
gate retains units, selected claims and applicability; the original tensor
rotation supplies `epsilon = Q diag(epsilon1,epsilon2) Q^T`. No conductivity
orientation or omitted normal response is inferred. The physical plane must
be decoupled from the omitted normal direction.

Snapshot resolution and trajectory admission remain distinct. Retained
charge owners require every principal coefficient to be explicitly constant
throughout source validity. Equal samples do not prove constancy. Changing
dielectric coefficients or axes during a trajectory still requires an
independently accounted parameter-work/thermal-entropy model; this increment
does not supply one. All source coefficients remain in the original complete
material bundle.


### Complete liquid bundles through temperature and shear evolution

`ResolvedLiquidState::try_from_material` adopts an already resolved material
bundle after checking the original four transport requirements and explicit
dynamic/kinematic viscosity basis. Additional coupled properties remain
attached. `LiquidTemperatureLaw::try_from_material` also checks the original
immutable card and temperature coordinate. Its existing exact-source refresh
therefore retains all selected claims, query meanings and receipts through
later temperatures, including an active surface-tension claim.

`ResolvedLiquidSurfaceState::try_for_shear_rate_evolution_with_capillarity`
keeps the scenario's explicit capillary choice in that same complete bundle
while preserving the original temperature/shear ownership. The old entry
continues to declare neglected capillarity. Admission of an evolving bundle
does not authorize freezing an unsupported coefficient; the receiving
runtime must own or refuse each active response. The thermal free-surface
consumer currently requires source-declared validity-wide constant sigma.
No surface-energy, Marangoni or wetting dynamics are inferred by adoption.


### Relative entropy on the admitted enthalpy/phase chart

`EquilibriumEnthalpyPhaseCurve::specific_entropy_change_j_kg_k(h0,h1)`
integrates `dh/T(h)` on the SAME bounded piecewise-linear caloric chart.
Sensible intervals use their logarithmic integral, evaluated with log1p near
isothermal limits; latent plateaus contribute `delta_h/T_transition`.
Fixed-order compensated accumulation, original endpoint/source-domain
admission and explicit representability checks apply. Reversing endpoints
reverses the sign. No apparent plateau heat capacity, additional latent energy,
absolute entropy reference, chemical mixing or variable-pressure work is
invented. Its thermodynamic scope is the chart's fixed pressure/composition
convention. Coupled owners retain responsibility for boundary entropy and any
additional state/work coordinates.

Two focused regression definitions cover independent sensible-plus-latent
integrals, subinterval additivity, reversal/domain refusal and a nearly
isothermal sensible limit. Source syntax parsing and independent scalar
arithmetic passed; native Rust tests were unavailable in the implementation
environment. This is not a continuous transport or entropy-production proof.


### Source-backed planar dielectric thermodynamics

`ThermalPlanarDielectricMaterial` pins one original complete material bundle,
selected claims, nonthermal context and explicit scalar or fixed orthotropic
dielectric axes. Positive globally affine principal epsilon(T) must be supported
throughout the declared absolute-temperature interval, including interior source
knots. It reuses the original affine-path, Helmholtz derivative and capacity
helpers, and re-evaluates the full source bundle at each accepted temperature.

At fixed displacement D and volume, the law supplies
`psi=V D.epsilon(T)^-1.D/2`,
`S=-partial_T psi=V E.epsilon_T.E/2`,
`U=psi+TS` and the nonpositive fixed-D field heat capacity. The field contribution
is relative to the original zero-field caloric chart, not another heat store.
A conservative capacity bound applies over the explicit T and |D| envelope.
A thermal consumer must add ALL field contributions assigned to each caloric
receiver and keep their combined heat capacity positive.

Source continuity compares original card, exact claims, full query semantics
and actual values; a requery's stronger pinned selection receipt is not confused
with a change in physical material. Axes and composition/pressure context stay
fixed. No dispersion, nonlinear polarization, phase-dependent entropy or
electrostriction is inferred. The focused source-derivative and coupled tests
are in the fs-couple dielectric thermal test child. Final source syntax parsing,
bounded source review and independent scalar arithmetic passed; native Rust
tests were unavailable, so no executed native or empirical qualification is
asserted.


### Source-bound thermal longitudinal piezoelectric potential

`ThermalLongitudinalPiezoelectric` pins the complete original longitudinal
material bundle, card, source claims, axis quantities and fixed nonthermal
coordinates. The owned absolute-temperature interval admits affine s_E(T),
signed d(T), epsilon_T(T), constant reference density, and positive definiteness
of the full coupled block A=[[s_E,d],[d,epsilon_T]] throughout. In natural
coordinates x=[strain,D], y=A^-1*x, the same native response supplies
Psi=V*x^T*y/2, S=V*y^T*A'*y/2, U=Psi+T*S, and
C_x=-T*V*(A'*y)^T*A^-1*(A'*y).

Signed/zero coupling is valid; signed-zero reference checks use the admitted
source-path scale without replacing actual source values. A normalized
full-block eigenvalue/derivative bound supplies a conservative minimum field
heat capacity over the caller's strain/D envelope. This is an admitted model
bound, not a numerical certificate. Common affine-source and capacity routines
are shared with the existing field-storage thermodynamics; positive-coefficient
source admission retains its original behavior. No separate constitutive engine,
pyroelectric offset, expansion, depoling, hysteresis, phase transition or
dissipative ceramic law is inferred.

Focused synthetic regressions in fs-couple compare the complete analytic
potentials, thermal derivatives, signed/zero coupling and actual coupled
circuit/motion/caloric evolution. Source syntax and independent scalar references
were checked; native compilation/tests were not run in the editing environment.


### Final nozzle admission after enclosing heat or mechanical work

`GasTransportStep::admit_endpoint_nozzle_pressures` reapplies the existing
active-nozzle pressure-order check to a complete endpoint parcel map after an
enclosing physical owner supplies heat or mechanical work. Original parcel IDs
and actual positive transfers remain authoritative; the adapter samples no new
flow, spends no inventory and changes no accepted parcel. It checks pressure
ordering only and does not extend the original transport-stage entropy account
to subsequent heat/work. The multi-chamber DG integration exercises a late
heated-reservoir reversal and complete rollback. Source/API review and syntax
parsing completed; native execution remains unverified in this environment.


### Source-backed reversible isotropic thermoelastic storage

`elastic::thermal::AffineIsotropicThermoelastic` retains the original material
card, selected elastic and expansion claims, full query semantics and declared
positive temperature interval. It admits positive affine Young's modulus,
constant Poisson ratio in (-1, 1/2), and constant instantaneous linear expansion.
Dimensionally equal semantic coefficient kinds are not interchangeable with
the native dimension-only elastic schema. External coordinates and selected
claims stay pinned through every refresh.

The original `IsotropicElastic` supplies stress. Explicit conversions between
its tensor-shear Voigt convention and the solid owner's Mandel convention
preserve work. At fixed total strain, the adapter evaluates elastic Helmholtz
energy Psi, entropy -dPsi/dT, internal energy Psi+T*S and fixed-total-strain
capacity -T*d2Psi/dT2. Total and elastic strain components must stay inside the
declared small-strain envelope. The stress-free reference is the retained
expansion source, not an inferred density or material-name temperature.

The full-interval conservative lower capacity bound is intended to be added
to the actual caloric inventory by a coupled physical consumer. This adapter
does not create heat storage, infer current density, evolve temperature, solve
mechanical equilibrium, supply plasticity or claim finite-deformation validity.
A consumer must separately retain its reference mass, admit the actual caloric
chart and check positive combined capacity.

Two focused regression definitions cover independent isotropic invariants,
pure shear work, fixed-total-strain thermal derivatives, signed modulus slopes
and expansion, exact relaxed/constant-source limits, capacity bounds and
source/schema/domain refusals. All three changed Rust files passed syntax
parsing. Independent finite-difference reference calculations and bounded
physics/API reviews completed. Native compilation, Cargo tests and rustfmt
remain unexecuted in this environment.

## Source-defined bilinear cohesive interfaces

`cohesive::source` projects an original ordered `InterfaceSystemCard` or an
already resolved multiphysics interface bundle into the existing
`BilinearCohesive` law. Explicit `cohesive_stiffness` [Pa/m],
`cohesive_strength` [Pa], `cohesive_fracture_energy` [J/m²], and
`cohesive_compression_stiffness` [Pa/m] are required. The numerical maximum
separation is supplied separately and must include complete failure. No bulk
material name, absent property, or geometric proximity creates an adhesive bond.

The full bundle, exact source selections, ordered surfaces and texture frames
remain accessible. The canonical fixed-external-condition policy applies to
both value dependence and applicability: evolving gap, damage, stress, rate or
frequency coordinates cannot be frozen into this rate-independent law. All
external conditions remain fixed during the accepted damage history; there is
no implicit coefficient/history refresh or inference of fracture heat.

The existing common-mode cohesive law still owns traction, tangents,
irreversibility, compression after failure and fracture-energy accounting.
This is not a calibrated mixed-mode adhesion model. The focused
`cohesive_source` regressions cover shared electrical/mechanical source bundles,
fixed-temperature substitution, independent fracture work and source/unit/domain
refusals; structural source-to-force tests live in `fs-solid`.

## Explicit independent isotropic elastic parameters

`state_point::resolve_isotropic_elastic_state_point_from_parameters` admits a
caller-selected Young/Poisson, Young/shear, or bulk/shear pair with density.
Arbitrary source property names and complete typed query conditions are retained.
The scalar role requirements reject stress, pressure, hardness and tensor-component
aliases; this adapter does not infer isotropy from a material name or directional
measurements. Positive, finite and numerically representable isotropic stiffness
is required before publishing the existing `IsotropicElasticStatePoint`.

The original `resolved()` bundle contains only the selected source properties.
`derivation()` retains their roles and nominal algebraic consequences separately;
no derived constant is inserted as a measured claim. `identity()` binds the
derivation and is used by the solid element, thermal-strain and reduction ingress.
Canonical Young/Poisson states retain their previous identities. Source marginal
uncertainty remains attached; joint covariance is not assumed and uncertainty is
not propagated. Additional unselected measurements are not reconciled automatically.

Focused regressions cover independent constitutive stress, source conditions,
invalid pairs, conflicting claim pins, unit rescaling, and the actual tetrahedral
mass/strain-energy and thermal-expansion identity path. This is a bounded MR11
parameter-conversion slice, not complete derived-property or material validation.

### Reversible adhesion driving equilibrium liquid wetting

`liquid::surface::adhesion::resolve_young_dupre_wetting` derives the existing
`ResolvedLiquidWettingState` from an explicit ordered interface claim
`equilibrium-work-of-adhesion` [J/m2] and the selected liquid/gas surface
tension [N/m]. The explicitly chosen ideal smooth, homogeneous, rigid-wall
model obeys `W = sigma (1 + cos(theta))`; the half-angle evaluation remains
well conditioned near zero angle. The model requires `0 < W < 2 sigma`.
Complete wetting, film formation, complete dewetting, roughness and hysteresis
require other physical owners and are refused rather than clipped.

The angle retains both original source receipts and is identified as a
Young--Dupre prediction, never inserted as a measured contact-angle claim.
The boundary checks the exact liquid card, selected tension receipt, ordered
surface roles, query point and axis semantics before the existing native LBM
wall operator consumes it. Refresh requires the original interface and the
original selected tension source at the same new physical point; no source
substitution or partial refresh is accepted. This is a frozen endpoint, not a
thermal surface-energy evolution. Directly sourced contact angles retain their
existing resolution and refresh path.

`YoungDupreLaw` also supplies the signed dry-substrate-reference surface free
energy `sigma A_lv + (sigma-W) A_sl` and unbalanced Young force per contact-line
length. These fixed-size evaluations allocate no memory. Reversible equilibrium
work is distinct from fracture toughness, peel dissipation and advancing or
receding hysteresis measurements; none are silently substituted. No droplet
pull-off force is inferred. The relation and energy convention follow
Nagy, Scientific Reports 14, 29991 (2024), equations (1)-(3),
<https://doi.org/10.1038/s41598-024-81710-6>.

Numerical regressions check the Young relation, force direction, independent
fixed-volume spherical-cap energy minimization, area scaling and domain
refusals. Native source/wall tests exercise distinct adhesion-driven capillary
responses, mass conservation, exact paired-source refresh, substitution refusal
and the unchanged directly sourced angle path.

## Mode-dependent cohesive fracture

`BilinearCohesive::new_mixed_mode` selects separate opening/shear stiffnesses,
peak strengths and fracture energies. Quadratic nominal-traction interaction
localizes the first damage-initiation point on each straight jump segment. The
undamaged shear-energy fraction at that point selects the Benzeggagh--Kenane
toughness. This implementation explicitly freezes that onset mix and envelope
for the point's subsequent history. A shared maximum energy norm then drives
irreversible scalar degradation, with undamaged compression after failure.
Pure-mode and proportional mixed-mode work equal their selected fracture energy;
nonproportional post-initiation paths retain exact analytic work/storage/loss
accounting and cannot reselect toughness or heal. The onset-mix freeze is a
declared model extension, not experimental qualification of arbitrary histories.

The consistent algorithmic tangent includes the sensitivity of the initiation
crossing and selected toughness when a trial first initiates damage. It can be
nonsymmetric; the existing Newton/FGMRES owner consumes it without an SPD
assumption. Evaluation allocates nothing; initiation and timed damage/failure
events use analytic quadratic crossings, split at the normal compression switch.
Radius-scaled coefficients and a rationalized outward root retain initiation
fractions below `2^-64`, preventing a long increment from changing the physical
strength threshold or causing premature rupture. Retained-envelope steps
use closed-form arithmetic. Equal mode
parameters recover the original common-mode response. The scalar law getters
continue to describe pure opening; each point exposes its selected mixed envelope.

`CohesivePoint::trial_with_mean_traction` additionally returns the analytic
three-component traction average along that same straight jump path. It splits
retained-history elasticity, softening, compression and complete failure using
the original initiation/envelope history. Affine-vector/norm integrals preserve
transverse traction under nonproportional loading; a short-interval binomial
evaluation retains small tractions near failure. The mean's dot product with the
full jump change must reproduce the original independent path work within the
supplied J/m2 gate. No scalar work projection, quadrature, new damage law or
allocation is introduced. Endpoint traction and endpoint tangents remain separate
and do not differentiate the path average. Regressions compare vector means with
independent traction quadrature, shear rotations, full rupture, reversals and
straight-path subdivision.

`resolve_mixed_mode_bilinear_cohesive` requires eight explicit ordered-interface
claims: normal/shear stiffness [Pa/m], normal/shear strength [Pa], mode-I/mode-II
fracture energy [J/m2], dimensionless BK exponent, and compression stiffness. No
missing mode is inferred. The existing `ResolvedBilinearCohesive` source bundle
and solid/P1/thermal point paths remain the owners; fixed external source
conditions and all original source-role refusals apply. The existing thermal
joint is normal-only and therefore consumes the pure-opening envelope through
the same exact point work/storage interface.

Quadratic initiation and BK interaction reference Camanho and Davila (2002),
NASA/TM-2002-211737, <https://ntrs.nasa.gov/citations/20020053651>. Tests check
pure-mode and proportional mixed fracture work, equal-mode recovery,
initiation-sensitive and retained-history tangents, nonproportional work,
straight-path subdivision, compression/shear, source-unit/domain refusals and
actual native-equilibrium load histories. Temperature/rate-dependent toughness,
fatigue, evolving mode-mix calibration and fracture heat partition remain separate.


## Solid-vapor constitutive laws

`sublimation` supplies the shared `fs_thermal::sublimation::SublimationLaw`
interface. `WaterIceSublimation` evaluates the IAPWS ice-Ih pressure curve on
50..=273.16 K; `CarbonDioxideSublimation` evaluates the NIST solid-CO2 Antoine
fit on 154.26..=195.89 K, with explicit bar-to-Pa conversion. Named curves
refuse extrapolation. Pressure derivatives are analytic. Their temperature-
dependent latent heat and derivative follow Clausius-Clapeyron with ideal
vapor and negligible solid specific volume; they are not independently
measured calorimetry or a real-fluid equation of state. Original sources and
approximations are documented beside each law. Density is explicitly supplied
by the caller and held constant; no porosity, accommodation or geometry is
inferred from the material name.

`ClausiusClapeyronSublimation` admits caller-calibrated reference pressure,
reference temperature, constant latent heat, molar mass, density and a bounded
validity interval. Log-pressure evaluation avoids avoidable overflow and keeps
the exact reference anchor. All scalar laws allocate nothing and use the shared
`fs_math` deterministic transcendental kernels. Ice saturation pressure is also
available independently of density for humidity consumers. The six focused
regressions check published pressure anchors, analytic derivatives, units,
latent-heat behavior, explicit domains and representability refusals. These
laws supply constitutive data; phase inventory, transport and energy evolution
remain with their physical consumers.


## Explicit ice-relative humidity

`gas::HumidityReference` and `try_new_moist_air_with_reference` select liquid
water or ice as the relative-humidity pressure reference. Ice uses the same
bounded IAPWS curve as the sublimation law. Existing moist-air entry points
retain their liquid-water reference, including below freezing; zero humidity
retains the original dry-gas domain. Actual vapor mole fraction flows through
the original ideal-mixture, Wilke and WMS evaluator, so acoustic and transport
consumers see the same physical composition. The existing RH [0,1] and dilute
water-mole-fraction ceiling remain in force.

`resolve_moist_air_state_with_reference` retains the reference alongside both
original component source bundles. Ice pressure does not relax either card's
applicability or qualify cold transport fits. These APIs prescribe a gas state;
they do not evolve condensation or change a downstream acoustic standard's
humidity convention. Seven native shared-mixture regressions passed, including
four new reference/domain cases. The source-resolution wrapper was syntax
checked and reviewed, but was outside that native numerical closure.

### Scalar magnetic hysteresis and uniform-core constitutive coupling

`field_storage::hysteresis::ScalarHysteresis` is an opt-in, fixed-condition
energy-based play model with explicit caller-supplied branch susceptibilities,
pinning fields and saturation polarizations. Each branch owns retained internal
polarization. The immutable law and accepted state stay separate; trial field or
flux updates do not mutate accepted memory or silently replace model parameters.
This supplies remanence, reversal/minor-loop memory and saturation independently
of the existing reversible B-H sources.

For B = mu_b H + sum(z), stored energy is
W = mu_b H²/2 + sum(z²/(2 chi)); irreversible heat is sum(Hc |delta z|).
The exact scalar constrained minimizer gives pinning/sliding/saturation behavior.
For each declared monotone segment, magnetic work is the analytical integral
H dB = delta W + Q, including crossed corners. Every field reversal must be
supplied; endpoint trapezoids do not replace this work. Positive background
permeability makes the B-driven inverse strictly monotone. Safeguarded Newton
uses the active-set tangent, has bounded work, and allocates state only for the
final candidate. Tangent values at kinks select the pinned/saturated side.

`UniformHystereticCore` exposes current- and linkage-driven candidates through
lambda=N A B and I=l H/N, including signed winding orientation. Core energy and
pinning heat scale by the actual A l volume; reported mean current is conjugate
to linkage work. Circuit trials must all use the same accepted magnetic history
until the complete coupled step accepts. This is a constitutive adapter, not an
alternate circuit time integrator or automatic spatial-field hookup.

The branch potential is quadratic with explicit saturation bounds, following
the scalar dry-friction energy construction discussed in Prigozhin et al.,
IEEE Trans. Magn. 52 (2016), doi:10.1109/TMAG.2016.2599143.
It does not reproduce the paper's general vector model or confer measured
material authority on supplied coefficients. Temperature, stress and parameters
remain fixed. Eddy/excess losses, dynamic/vector hysteresis, demagnetization
damage and thermal parameter work require separate models. Never add a loss map
for the same resolved hysteretic heat. Eleven focused tests cover analytic
segments/loop area, minor memory, monotone subdivision, inverse histories,
tangents, reversible limits, refusal and winding energy/heat conjugacy.


### Isotropic vector magnetic hysteresis

`field_storage::hysteresis::vector::VectorHysteresis` extends the original
fixed-condition branch parameters to three-dimensional polarization memory.
Each branch minimizes |z|²/(2chi) - H.z + hc|z-z_old| inside its saturation
ball. Exact interior shrinkage and a bounded scalar multiplier solve give the
constrained vector minimizer. Its consistent active-set dz/dH is symmetric
positive semidefinite; the positive background permeability makes the complete
dB/dH positive definite for the bounded, backtracked B-driven inverse.
All trials start from the same accepted memory. Branch-state storage allocates
only for the accepted field candidate, not for inverse trial responses.

Reports distinguish physical pinning heat sum(hc|delta z|), stored-energy
change, endpoint H_new.delta B, and nonnegative numerical dissipation. The
numerical term is not heat; endpoint work is not an exact arbitrary rotating
path integral. Resolve direction changes and refine the history. Original
source/field/flux domains, constitutive and work residuals, bounded iterations
and cancellation gate acceptance without modifying the previous state.

Eight physical tests cover the scalar state/heat limit, independent spherical
minimization, rotation covariance including tangents, an analytic rotating-loop
loss limit and history refinement, finite-difference tangents, B-driven
histories, reversible limits and atomic refusal. Eleven original scalar tests
remain passing. This is an explicit fitted, fixed-temperature/stress
quadratic-potential specialization of Prigozhin et al. (2016), not experimental
material validation, texture anisotropy, dynamic/eddy/excess loss,
demagnetization damage or automatic field/circuit/thermal integration.

### Compressible porous leakage and selective membrane permeation

`porous::permeability` separates intrinsic permeability [m²] from gas viscosity,
dense-membrane molar permeability [mol/(m·s·Pa)] and thickness-normalized
permeance [mol/(m²·s·Pa)]. The porous owner compiles orthotropic normal
permeability and series layer resistance once. Its isothermal ideal-gas Darcy
law conserves mass with distinct face volume rates, includes an explicitly
gas-model-bound per-layer Klinkenberg factors when calibrated, and refuses unsupported
temperature, absolute pressure, speed or material conditions. A zero normal
permeability seals the path exactly. The acoustic gas-substitution owner uses the
same intrinsic-permeability type. This is dry rigid-skeleton quasi-steady flow;
Knudsen/transition diffusion, pore storage, wet multiphase flow, pore heat storage,
Forchheimer inertia and evolving deformation are not inferred.

`thermochemical_gas::energy::parcel::transfer::species::network::transport::step_candidate_with_paths` accepts `GasTransportPaths` containing nozzles,
Maxwell–Stefan diffusion, porous barriers and dense membranes. All rates sample
the same original inventory and enter one conservative species/partial-enthalpy
exchange. Gross outgoing species allowances are shared across every mechanism;
incoming gas cannot replenish those allowances or traverse another edge during
the same explicit step. Final endpoints must preserve active total-pressure or
species partial-pressure order and pass the existing extensive entropy gate.
The original two-mechanism API delegates with empty porous/membrane slices.

The membrane owner uses source-selected, constant-at-fixed-temperature molar
permeability and actual species partial pressures, allowing counterpermeation.
It checks the selected claim at both gas boundaries, retains explicit zero
permeability for blocked species, and refuses missing species, stale temperature
aliases, wrong units and unsupported pressure-dependent sorption. The closure
assumes a steady negligible-storage, near-isothermal dense film with dilute Henry
sorption; competitive sorption, plasticization, swelling, transient breakthrough,
nonideal fugacity and thermal transpiration need additional constitutive laws.
Authored mixture/humidity applicability remains a source requirement.

Analytic regressions cover compressible and layered rates, anisotropy, gas
substitution, slip calibration, selective flux and reversal, common species
withdrawal, mass/enthalpy conservation, equilibrium overshoot and cancellation.
These are constitutive and numerical checks, not validation for an unnamed
porous solid or polymer. The existing measured acoustic fixtures remain active.


## External species streams through the original gas energy owner

`thermochemical_gas::energy::parcel::transfer::species::boundary` accepts
integrated species mass and complete stream energy on the original gas
thermochemical reference. `ClosedGasParcel::exchange_species_boundary` reuses
the existing species-inventory, composition, EOS, transport and relative-energy
inversion paths to return one immutable endpoint with extensive balances.
Carrier inventories, volume and external query coordinates remain intact.
Exact individual-species depletion and reentry use a retained complete source
catalog; no trace gas or replacement caloric law is inserted. Active caloric,
transport, source-pack and water-identity substitutions refuse. Candidates
borrow the immutable catalog to avoid copying source cards on each trial.

Original-owner queries expose species mass, partial pressure and partial
specific enthalpy. Opposing streams may carry energy with zero net species
mass. The external phase must debit the same signed energy on a compatible
reference; latent heat alone does not determine NASA-9 vapor enthalpy. This
primitive supplies neither an interface rate nor automatic thermal-reference
alignment. Six focused native regressions passed for influx/deposition,
depletion/reentry, reference shifts, counterflow, ordering and atomic refusal.

### Thermal entropy adapters

The shared ThermalLaw adapter for EquilibriumEnthalpyPhaseCurve delegates to
the original exact piecewise dh/T integral, including latent plateaus once.
GasEnthalpyLaw supplies the same optional entropy interface by reconstructing
both enthalpy endpoints with its original NASA/EOS law and evaluating the
cancellation-safe frozen-mixture entropy change at the retained pressure.
Source domains, bounded inversion and frozen-composition restrictions remain
in force. These adapters enable shared thermal entropy audits without a second
caloric law. Independent sensible/latent and variable-cp gas regressions cover
storage, endpoint reversal, no-ops and refused source-domain exits.


## Molecular sublimation with evolving gas energy

`sublimation::gas::advance_surface_with_gas` couples one original solid caloric
law to the original `ClosedGasParcel` and its retained complete source catalog.
Both temperatures, species inventories, pressure and gas properties evolve.
The explicit species binding supplies a fixed solid enthalpy reference offset;
source molar mass must agree, and h_v(T)-(h_s+offset)=L(T) must hold within the
caller's declared absolute tolerance at both accepted surface endpoints.
Incompatible caloric temperature dependence refuses rather than being hidden
by a new energy datum.

The shared molecular law supplies gross incoming and outgoing streams. Each
stream carries the original species h(T)-R_specific*T/2, with separate surface
and gas temperatures. Counterflow can therefore transfer heat at zero net mass.
The original gas boundary owner reconstructs the endpoint once. The candidate
checks actual combined mass and energy and retains signed recession; neither
input changes until the surrounding consumer accepts both endpoints together.

This is a bounded first-order explicit interval. Both gross withdrawals must
fit their initial inventories and the caller's maximum fraction, even when
net mass is zero; timestep refinement resolves thermal and pressure feedback.
Unsupported domains, incompatible references, unrepresentable increments,
exhaustion and cancellation refuse atomically. The local Maxwellian/sticking
assumptions exclude a gas film, carrier thermal accommodation, external heat,
resolved flow and geometry changes. Multi-surface shared-chamber evolution is
separate and retains its thermostatted contract.

Six new native coupled regressions passed with the original thermochemistry,
EOS, transport and gas-energy implementation: dry-carrier sublimation,
deposition, zero-net-mass heat transfer, formation-reference invariance,
timestep convergence and atomic refusal. The six existing species-boundary
regressions also passed in the same actual-source native executable.


### Source-derived two-surface thermal Hertz compliance

`elastic::thermal::contact::ThermalHertzPair` keeps both original isotropic solid cards, complete query contexts and exact source selections. Each body has its own actual temperature. Continuous source lines admit positive affine Young's moduli and fixed Poisson ratios/yield limits; unsupported derivative jumps, hidden evolving axes and extrapolation refuse. The original `IsotropicSolidStatePoint::reduced_modulus_against` computes the harmonic compliance. Its full two-temperature gradient and Hessian are analytical; the effective modulus is not approximated as affine. Immutable source data are shared between candidate checkpoints.

For a fixed-geometry native Hertz potential `psi=G(delta)*E*`, the derivatives support surface entropy `s_i=-G*dE*/dT_i`, internal energy `U=psi+sum(T_i*s_i)` and reciprocal thermal capacity. Returned source states retain their original identities and receipts. This is a constitutive capability, not a contact or thermal integrator; coupled owners must retain physical heat/work, total heat capacity, source domains and mechanical state. No geometry expansion, changing reference mass, phase change, plasticity, measured calibration or full Euler thermal closure is inferred.

Focused native G1/G3 tests compare the modulus and all first/second derivatives with independent rational compliance and finite differences, check thermodynamic internal-energy derivatives, ordered surface exchange, the constant-modulus limit, and original-card/domain/non-affine-source refusal. Synthetic source fixtures establish implementation behavior only.

### Physical cohesive viscosity

`BilinearCohesive::with_viscosity` selects explicit Kelvin--Voigt traction
`(1-d) eta_i gdot_i` with independent normal/shear interface viscosities in Pa s/m.
`trial_over_time` requires positive physical time and a declared jump-rate limit;
untimed trials refuse this law. The original intrinsic common-mode or frozen-onset
BK fracture owner still owns damage and storage. Normal viscosity acts in both
opening and compression while bonds survive, and all viscous traction vanishes
after complete failure; the separate compression penalty remains.

Endpoint total traction and its consistent, possibly nonsymmetric tangent are
separate from the original equilibrium account. Exact piecewise analytic
constant-rate path integration reports nonnegative viscous loss independently
from intrinsic fracture. Energy uses that integral, not endpoint traction times
jump. The short-interval binomial evaluation has below-roundoff truncation in
its admitted range. Longer ramps evaluate the logarithmic integral directly,
preserving nonzero loss at large failure/onset ratios without cancellation of
the leading terms. G0/G3 regressions compare full-rupture viscous work with its
independent closed form and retain the same work under reversal-path subdivision.
`trial_over_time_with_mean_traction` adds the same exactly integrated surviving
bond fraction times `eta_i gdot_i` to the equilibrium vector mean. Its three-mode
work equals intrinsic work plus the separately reported physical viscous loss;
it reuses the timed owner's integral rather than a second damage evolution.
No stored dashpot state, healing, fracture heat partition,
evolving-temperature refresh or universal adhesive calibration is inferred.

The source adapter queries intrinsic coefficients and explicit interface
viscosities in the same ordered original bundle, retains receipts, and rejects
bulk viscosity units and unowned evolving axes. Intrinsic fracture work must
exclude viscous losses to avoid double counting measured finite-rate work.


### Heterogeneous calibrated slip in layered gas barriers

`porous::permeability::DarcyGasBarrier` admits a different calibrated Klinkenberg
factor in each layer. One bounded scalar solve conserves mass through the
original ordered layers; factors are never averaged into one apparent material.
The exact common-factor fast path remains constant-time and ordinary flow
queries remain allocation-free. The shared support intersection and actual gas
model/composition admission still apply to every layer. The continuum linear
resistance adapter refuses any nonzero slip.

`mass_flow_profile[_for_model]` additionally returns the absolute layer-interface
pressures, iteration count and normalized pressure-closure residual. Sealed and
zero-area paths do not invent interior pressures. Reversing both the physical
stack and its reservoirs reverses the flow; reversing only the reservoirs can
change its magnitude when the layer factors differ. Equilibrium conductance
uses the exact small-signal series resistance. This is rigid, isothermal,
constant-composition slip flow, not a Knudsen or thermal-transpiration model.

Six new native regressions and all 19 existing porous/acoustic regressions pass,
including independent per-layer mass balance, manufactured interfaces, stack
orientation, subdivision, adjacent-float pressure differences, exact closed
limits, calibration refusal, extreme slip contrast and bounded 64-layer work.


## Joint nonisothermal multi-surface sublimation

`sublimation::gas::chamber::advance_chamber_with_gas` advances heterogeneous
surfaces against one original gas snapshot. Each surface retains its explicit
identity, species, phase law, caloric callback and fixed energy-reference offset.
Every participating species is evaluated once at the original gas state.
Canonical surface/species accumulation precedes exactly one original gas
boundary reconstruction; input ordering cannot change the accepted result.

Gross incoming streams share each species' OLD inventory allowance. Gas emitted
by another surface cannot fund an additional withdrawal in the same explicit
interval. Each solid separately bounds gross emission. `InventoryLimit` returns
a conservative retry duration instead of clipping streams. Common-reference
energy aggregation preserves heat transfer during counterflow without repeatedly
adding large opposing formation energies. All solids and gas pass actual total
and per-species mass balance, original caloric compatibility and energy checks
before one immutable joint candidate is returned.

The existing single-surface API shares the same physical preparation. Fixed
area, local Maxwellian vapor, original source domains, specular nonsticking
reflection and pure-solid restrictions remain. This explicit model does not
supply gas transport, mixed phases, extinction topology or external heating.
Six new native regressions cover warm/cold redistribution, different volatile
species, exact single-surface equivalence, order independence, shared gross
withdrawal, energy-reference invariance and atomic refusal. All 18 actual-source
joint/scalar/species-boundary tests passed with the original gas reconstruction.


### Adaptive closed-chamber histories

`sublimation::gas::chamber::advance::advance_chamber_adaptive` compares one joint
explicit interval with two half intervals and accepts the finer conservative
history without extrapolation. Local absolute/relative tolerances cover every
solid mass and temperature, gas temperature and EACH retained species mass,
including trace vapor. Stable surface ordering makes states and retry decisions
independent of caller ordering on the same ISA. Typed inventory-limit refusals
reduce the whole interval; invalid source/caloric data and cancellation remain
immediate failures. Trial, minimum-step and maximum-step controls bound work.

The candidate retains only accepted endpoints and cumulative molecular streams.
Actual full-history total and per-species mass, combined energy, and each side's
stream/storage energy account must meet the original absolute balance budgets;
those tolerances are not multiplied by step count. Gas/solid inputs remain
unchanged on a late refusal. Areas and source bindings stay fixed throughout
this driver; thermal events use their separate owner. Step doubling supplies a
local time-error estimate, not a rigorous global accuracy or physical-validation
certificate. Six focused adaptive tests passed, including convergence toward a
2048-step reference, shared-inventory retries, trace-species refinement, exact
order invariance, work limits and late cancellation. The combined native
adaptive/joint/scalar/boundary suite passed 24 tests with no failures.


### Coupled molecular, Knudsen and viscous transport in porous gas paths

`thermochemical_gas::energy::parcel::transfer::species::diffusion::dusty_gas`
implements the isothermal ideal-gas dusty-gas equations in the stationary-solid
frame. One homogeneous dry pore path couples all selected binary molecular
frictions, explicit Knudsen wall frictions and the original no-slip Darcy owner.
There is no zero-molar or zero-mass-flux constraint: unequal counterdiffusion and
uphill individual-species flux can follow from the coupled molecular friction.
A separate Darcy/nozzle/diffusion edge must not count that same pore area again.

The immutable medium retains the original intrinsic-permeability barrier and
explicit pore radius, gas-filled porosity, transport tortuosity and provenance.
Effective molecular diffusivity is porosity/tortuosity times the source value;
Knudsen diffusivity uses the actual thermochemical species molar mass. The
transport tortuosity is the declared coefficient in that equation, not an
inferred squared path ratio or the acoustic high-frequency tortuosity. Added
Klinkenberg slip and heterogeneous pore layers refuse to prevent overlapping
rarefaction models or an invented unresolved interface composition.

The face composition is concentration weighted at the declared mean temperature.
Original selected binary source laws are sampled at the face and admitted at
both endpoints, including pressure-dependent diffusivity. The shared material
owner supplies mixture viscosity and thermochemistry. A symmetric positive
friction system uses the existing dense Cholesky owner, retains every species
row and rechecks the original equations after total flux assembly. Work is
bounded to 32 species; loops poll cancellation. Reports retain source states,
individual Knudsen/diffusive/viscous rates and nonnegative local dissipation.
Nonrepresentable nonzero forces refuse instead of becoming false equilibrium.

`network::transport::GasTransportPaths::dusty_gas` samples every face from the
unchanged network snapshot. Its signed species transfers share one gross
withdrawal allowance, original donor partial enthalpies, material reconstruction
and final conserved-species entropy admission with every other mechanism.
Counterflow does not replenish an outgoing allowance. Local dissipation is a
diagnostic and is never added again as gas heat. Zero area queries no inactive
sources and returns no invented face state; zero intrinsic Darcy permeability
suppresses only viscous flow when explicit conducting pore diffusivities remain.

This is a bounded central face approximation, requiring spatial and physical
time refinement. It does not resolve thermal transpiration, Soret effects,
adsorption/surface diffusion, water saturation, changing pore geometry,
pore storage, pore-size distributions, inertial Forchheimer flow or a nonlinear
steady composition profile through a finite heterogeneous laminate. Source and
geometry applicability must be supplied for the actual material.

Before an executor restart, a native executable of the actual owners passed 89
distinct focused tests, including nine new dusty-gas regressions: pure-gas
Knudsen plus compressible Darcy limits, Graham counterdiffusion, an independent
ternary uphill-flux reference, concentration weighting, source/units/closed-path
refusals, absent-side species, finite mass/energy conservation, shared nozzle
withdrawal and constitutive admission. An additional extreme-range underflow
regression and the merge onto the concurrently improved entropy API were
source-reviewed after recovery; full native rerun and workspace validation
remain unverified because the restart removed the dependency caches.


### Runnable evolving-gas ice cold trap

`cargo run -p fs-material --example nonisothermal_ice_cold_trap` advances two
10-micrometer ice films at 260/252 K through one 100 mL water-vapor chamber,
initially 256 K and 130 Pa, for 2 ms. The original NASA-9 H2O coefficients and
IAPWS ice pressure/latent law remain the thermodynamic owners. The example
explicitly derives h_s(T)=h_v(T)-L(T)-datum on 250--270 K; this compatible reduced
ice chart is not independently measured ice calorimetry. Original dilute-vapor
transport remains within its stated sub-triple continuation. Fixed film area,
917 kg/m3 density, accommodation one and a well-mixed local vapor are explicit
scenario assumptions; no chamber-design or spatial-flow validation is claimed.

The actual native example passed both default and tenfold tighter local-error
runs. Warm loss/cold deposition were 26.567540/23.881098 micrograms by default,
with gas temperature 255.821736 K; tighter results were 26.564590/23.878148
micrograms and 255.821858 K. Both runs assert whole-history original-datum mass
and energy conservation and report gross counterflux, adaptive work and residuals.
Run with argument `0.1` to tighten local temporal error targets tenfold.


### Evolving gas on the original multi-cell heat timeline

`sublimation::gas::timeline::advance_deposition_with_gas` composes original
thermal transport with simultaneous molecular exchange against one evolving
gas. Every surface retains its original cell identity, phase/caloric source,
species binding and fixed datum offset. Original heat observers check each
selected pure-solid enthalpy branch after EVERY accepted thermal substep and
birth. One joint phase interval then updates the actual masses and enthalpies
through the thermal owner's transactional boundary.

Start-of-interval exposure and exact birth/contact times remain authoritative.
The original thermal work ceiling applies to the whole requested advance.
Empty pre-birth periods jump to their next event without renewing empty work.
Actual combined gas/solid storage is compared with incoming births and external
heat, using stable initial/birth-relative energy changes and explicit datum
conversions. Late source/domain/balance failure or cancellation discards the
complete candidate, including elapsed time and consumed events.

This fixed-area heat/phase composition is first order; refine the thermal
max_dt_s to resolve splitting and explicit molecular exchange. The standalone
adaptive chamber driver does not certify this operator-splitting error. The
caller still owns pure-solid admission, supported source domains and topology;
melting, extinction, moving geometry and gas flow are not inferred. Ten new
actual-source native integration tests passed for physical heat/gas feedback,
order, births/datum shifts, exposure events, work limits, solid-branch refusal,
conservation and rollback.

## Sourced sharp solid/liquid melting

`phase::melting::ResolvedMeltingCurve` joins two continuously supported source
Cp/density branches and independently sourced melting temperature/fusion
enthalpy into the original `EquilibriumEnthalpyPhaseCurve`. The caller explicitly
declares one sharp equilibrium transition. Solid support must end, liquid
support must start, and the transition query must lie exactly at the sourced
melting temperature. Every query carries an explicitly typed absolute
temperature. Pressure, composition and other external context, including quantity
descriptors, must be identical throughout; an optional typed categorical source
selector is zero for solid, one for liquid, and absent from the transition query.
The constructor refuses missing density/Cp/fusion data, disconnected or discrete
support, incompatible context, nonpositive fusion energy, unrepresentable energy
increments and combined charts exceeding the original generated-knot limit.

The existing heat-capacity owner integrates sensible heat on both branches. One
fusion increment separates the solid/liquid endpoint enthalpies at the exact
same temperature; independently tabulated assigned enthalpy is not also added.
The resulting chart, its states and relative entropy remain consumable by the
original enthalpy owners. Density retains branch-specific values and the original
additive-specific-volume coexistence rule; specimen mass stays invariant and the
consumer must account for the resulting volume change. Every selected Cp,
density, melting-temperature and fusion-energy receipt is retained and bound
into the curve identity. The declared sensible interpolation allowance is not a
source-uncertainty claim. This admission does not infer a phase diagram, nucleation,
kinetics, supercooling, hysteresis, moving geometry or a fixed-volume transport
closure.


### Prescribed wet-pore gas response

`porous::wet` supplies explicitly calibrated gas relative permeability and pore
diffusivity at a prescribed liquid saturation, total porosity and gas-pore radius.
The nonwetting van Genuchten--Mualem relation uses supplied residual liquid/gas
saturations and shape parameter. Diffusion independently selects either the
Millington--Quirk 1961 relation or an explicit transport tortuosity. These are
separate mobilities: residual Darcy immobility does not invent a diffusion
connectivity threshold. Full liquid saturation returns exact zeros and no
fabricated finite gas tortuosity.

Immutable samples retain the complete calibration through shared storage.
Sampling uses bounded deterministic scalar work with no sample allocation,
retains positive transport near closure, and refuses unrepresentable positive
coefficients. Domain bounds apply before the declared endpoint extensions.
Porosity and radius evolution are prescribed; intrinsic permeability, liquid
mass, capillary pressure, mechanical work and hysteresis are not inferred.
Nine focused native tests passed with the original MaterialError and fs_math,
including independent analytic limits, exact cube-root diffusion fixtures,
geometry scaling, monotonicity, source-domain admission and numerical limits.


## Lower standard atmosphere for flight consumers

`gas::atmosphere::StandardAtmosphere1976` derives the seven hydrostatic layer
bases from the U.S. Standard Atmosphere 1976 surface state and table 4 lapse
rates. Geometric and geopotential altitude have separate query methods; the
source's effective-radius mapping and gravity agree with the pressure
integration. Each accepted sample uses the original `GasState` dry-air owner
for density, sound speed and transport. Layer-base pressure is derived once
instead of using independently rounded values. Successful construction and
queries use fixed arrays without heap allocation. Near layer boundaries,
`ln_1p` retains small increments. Repeatability is same-platform/toolchain;
no cross-ISA bit-equality claim is made for the host logarithm.

The admitted branch is 0..80000 m geometric altitude, where the source retains
constant molecular mass. Above that boundary, USSA table 8 requires additional
composition variation; this implementation refuses extrapolation. The profile
is a dry standard reference, not measured weather, wind, humidity, terrain or
a rotating/inertial Earth model. Its effective hydrostatic gravity must not be
silently substituted for a different flight-frame gravity law. Focused tests
cover source layer values, the common gas state, layer continuity, the altitude
crosswalk, hydrostatic balance and upper-domain refusal.


### Saturation in the shared gas transport owners

`PorousGasLayer::with_wet_state` retains intrinsic permeability, original axes,
source records and gas support while applying the sampled gas relative mobility.
`DarcyGasBarrier::with_wet_states` recompiles each original series layer once;
repeated saturation updates replace their preceding factors. A flooded layer
closes the pressure-driven path exactly. Existing Darcy consumers receive this
behavior through the same barrier. Dry-gas Klinkenberg calibrations refuse wet
use rather than silently extending gas-specific slip authority.

`DustyGasMedium` distinguishes reference pore geometry from its effective wet
geometry. The molecular and Knudsen terms use the independently sampled absolute
diffusivity factor; the Knudsen term also uses the current gas-pore radius.
Relative Darcy mobility is not multiplied into diffusion. Full liquid saturation
closes all gas transport before inactive gas or binary sources are queried.
The complete wet calibration remains attached to the reported medium.

These are prescribed instantaneous constitutive updates. They do not move liquid,
resize a retained gas inventory, account for displaced-gas work, infer pore
connectivity or evolve a hysteresis branch. Those operations require a coupled
storage/geometry owner. The active gas law retains its original domain and
finite-step mass, energy and entropy admission.

Eight added integration regressions cover independently derived pure/binary
flux scaling, residual-gas diffusion, exact flooding, repeated wet updates,
radius dependence, retained source identity/slip refusal and anisotropic series
layers. Both modified production files and all new tests passed Rustfmt parsing.
The nine standalone wet-law tests passed natively. The complete selected-source
Rust test typecheck also passed after restoring the real gas dependencies.
Native execution of the eight added integration tests remains blocked by full
executor storage; it is not claimed here.


### Transient heterogeneous pore-body gas storage

`network::transport::pore_domain` retains finite reservoirs and an ordered set
of actual gas inventories. Each pore body's gas volume is its authored gas-filled
porosity times bulk area times thickness. Initialization uses the original EOS;
provided inventories must already match that volume. Independently authored
homogeneous throats retain their own length, area, material and dusty-gas closure.
No averaged heterogeneous interface law or continuous laminate profile is inferred.

Each accepted substep samples every throat from one unchanged gas snapshot and
uses the existing joint species/enthalpy reconstruction. Incoming gas cannot fund
another withdrawal within that snapshot. Typed inventory and caloric-solve
refusals trigger bounded timestep halving; invalid sources and domains propagate.
A failed late substep or cancellation discards the complete requested window.
Whole-window mass, individual species, energy and entropy are rechecked from the
actual initial and final inventories using the original absolute budgets.

Reports aggregate directional species and donor-enthalpy throughput, pore uptake,
work counts and whole-window balances. Only the last substep's detailed source
report is retained. The total omitted throat void volume is measured against
stored pore-body volume and must fit an explicitly declared geometric bound.
This first-order fixed-volume model needs temporal refinement. Moving porosity,
solid heat, sorption and a resolved nonlinear throat profile are separate owners.

Five focused regressions cover delayed forwarding, first-order convergence
against the independent three-volume diffusion eigenmodes, heterogeneous
storage/throat coefficients, bounded adaptive withdrawal and atomic refusal.
New files passed Rustfmt parsing and the selected actual-source Rust test
typecheck completed successfully. Native execution remains blocked by full
executor storage. These five tests are not claimed as run.


### Finite gas/solid adsorption with binding-energy storage

`adsorption` binds one retained gas species to a homogeneous Langmuir site
population, an original capacity/affinity/kinetic material source and two original
single-phase caloric charts. Capacity is per dry skeleton mass, affinity has
inverse-pressure units, and the bound adsorbate has an explicit enthalpy-datum
offset into the original NASA gas reference. The reference affinity fixes its
entropy datum. Temperature-dependent affinity follows the same gas/bound caloric
integrals and satisfies the corresponding van't Hoff relation; adsorption heat
is their enthalpy difference and is not added again as an external heat source.

Frozen original partial pressure and surface temperature give an analytic site
update and both integrated gross streams. Gross adsorption is limited by the
original gas inventory; total site turnover is bounded separately. Newly adsorbed
molecules may desorb within this exact frozen site solution. Original gas bulk
enthalpy travels with each direction. One finite gas reconstruction and one
bracketed, safeguarded-Newton solid caloric inversion produce a joint immutable
candidate. Both temperatures evolve; the original source domains, thermal-change
bounds, final affinity, mass/energy budgets and combined gas/site entropy must
all admit the endpoint. Cancellation or any late refusal leaves inputs intact.

The source adapter rejects stale temperature/pressure aliases and unowned
conditions. Adsorbate caloric, transport and pack identities remain bound.
Stable relative caloric operations cancel arbitrary reference energies before
thermal inversion and entropy admission. The shared
`ResolvedHeatCapacityCurve::specific_entropy_change` now integrates the original
chart's chord Cp/T directly between temperatures, preserving tiny thermal changes
even when a large enthalpy datum would round their absolute h values together.

This rigid, near-isothermal local kinetic closure uses one noncompetitive site
population and a declared bulk-enthalpy stream convention. It does not supply
free-molecular impingement, surface diffusion, liquid displacement, swelling,
hysteresis or solid heat conduction. Time refinement is required for the finite
bath's pressure/temperature feedback. It does not impose a hidden thermostat.

Eleven adsorption regressions and one new shared caloric regression cover exact
frozen kinetics and turnover, independent affinity/entropy derivatives, finite
uptake/desorption energy and entropy, equilibrium turnover, reference invariance,
piecewise caloric inversion and atomic/source refusals. All pass actual-source
Rust test typechecking together with the wet/DGM/pore-domain tests. Native
execution of these new coupled tests remains blocked by full executor storage.


## Variable-caloric off-design converging-diverging nozzle

`RocketNozzleRegime` distinguishes `UnchokedSubsonic`, `ChokedSupersonic`,
`ChokedExternalCompression` and `ChokedNormalShock`. Subsonic exit pressure at
the actual exit/throat area ratio determines unchoking; sonic throat pressure
alone is not the back-pressure threshold. A valid near-stagnation subsonic
solution does not require source support at a hypothetical sonic station.

`ClosedGasParcel::normal_shock_flow` treats its unchanged source parcel as the
upstream static state. Actual mass and momentum invariants determine the
velocity/temperature relation, and a bounded variable-caloric energy solve
factors out the trivial unchanged upstream root. The accepted shock reports
original source-resolved upstream/downstream states, actual entropy production
and conservation residuals. Internal nozzle shocks are located by section
area, followed by the same downstream isentrope and subsonic pressure recovery;
axial shock location requires a separately supplied nozzle contour.

`RocketNozzleFlow` retains throat Mach, regime, optional shock report and an
optional back-pressure residual for subsonic discharge. Its existing entropy
energy residual now subtracts the separately reported physical shock entropy
rise. Total axial thrust remains `mdot*u_exit+(p_exit-p_ambient)*A_exit` with
pressure thrust counted once, including at shock-on-exit transitions. Inner
caloric stations allocate no material receipts; accepted physical stations
retain the original complete property/transport evidence. Vacuum and low back
pressure skip the unnecessary subsonic-threshold solve. Existing blowdown and
vehicle consumers use the same entry point and original accepted impulse.

The model is quasi-steady, quasi-one-dimensional, inviscid, adiabatic frozen
ideal gas on an ordinary branch with one sonic crossing. External compression
is the attached ideal-flow branch, not a claim that a real nozzle will remain
attached. No viscous separation/side-load prediction, oblique plume solution,
shock train, transient shock tracking or reaction relaxation is supplied.
Equal/higher ambient pressure still refuses forward propulsion flow.

Native actual-source checks: 38 nozzle-family tests passed, including perfect-
gas shock oracles from Mach 1.0001 through 5, variable-cp conservation and
entropy, independent internal shock areas, unchoking thresholds, pressure-thrust
continuity, restricted caloric/transport domains, cancellation and actual finite
blowdown through the shock regime. The lower standard atmosphere's three native
tests also passed. This is a scoped native build, not a full workspace result.


## Equilibrium chamber to source-preserving frozen products

`ThermochemicalGasModel::equilibrium_chamber_hp` and
`equilibrium_chamber_hp_enthalpy` bind the original neutral ideal-gas
composition owner to the existing material parcel owner. The material bank
must cover every chemistry species and retain exactly equal complete NASA-9
models, including conventions, coefficient cards, provenance and caloric
corrections. Feed amounts are independent of the bank's catalog fractions.
Chamber pressure and all other material query axes come from the typed context.

A common inlet temperature or an explicit sum of original stream enthalpies
sets the HP target. Formation enthalpy enters once. The accepted composition
uses the existing bounded generated-simplex closure, reports its maximum
roundoff correction, and refuses eligible positive traces that cannot be
represented by the frozen material owner. `with_composition` preserves active
thermodynamic models, transport selections and source-pack identities.

The final source-bearing material state is resolved once. Equilibrium mass
and its actual density determine volume, and the original `ClosedGasParcel`
is assembled without a duplicate transport evaluation. Publication checks
actual parcel mass against both feed and equilibrium amounts and actual
parcel enthalpy against both the requested target and the equilibrium result.
Original chemistry receipts, mass/energy residuals and frozen parcel are
retained together; source mismatch, transport-domain refusal, unrepresentable
trace, failed balance or cancellation publishes nothing.

The resulting storage, acoustic, caloric and nozzle composition is explicitly
frozen. The HP closure is an initial thermodynamic chamber state, not a finite-
rate combustor or a continuously fed engine. No reaction heat is added after
freezing; equilibrium nozzle relaxation, kinetics, ions, condensed phases,
real-gas behavior and flame stability remain separate physical models.

Native actual-source verification passed all 152 material thermochemical-gas
family tests, including eight new chamber/pack integration tests, the existing
storage/mixing/transfer family and off-design nozzle tests. The new cases
include an independent exothermic HP root followed by the actual nozzle,
different inlet temperatures, complete/exact source-bank admission, structural
absence versus trace underflow, original admitted-pack decode and receipt
retention, final source refusal and cancellation. Independent equation/source
review found no conservation, reference or trace-freezing blocker.

`cargo run -p fs-material --example equilibrium_rocket -- 300 3000000 100000 2`
runs a manufactured variable-cp A/B source example through the actual HP
chamber and frozen nozzle. Its native execution produced 902.699280458 K chamber
temperature with 3.47e-18 kg mass residual and 3.64e-12 J enthalpy residual.
These values validate that authored integration fixture, not real propellant
or engine performance. The full repository Cargo workspace was not run.

### Sealed constant-volume equilibrium and conserved parcel chemistry

`energy::equilibrium::volume` connects the original TV/UV chemistry owner to
the same source-bearing gas material. `equilibrium_chamber_uv` takes a feed
temperature; `equilibrium_chamber_uv_energy` takes explicit total internal
energy. Both hold the declared volume and original feed mass, allowing molecular
count, pressure and temperature to change. Original formation energy is included
once; constant-pressure enthalpy and frozen heat-capacity formulas are not used
to solve the reacting endpoint.

The common product bank must contain exactly the chemistry's complete NASA-9
models. The existing composition-freezing boundary retains active transport
cards, source packs and external query axes, and refuses unrepresentable positive
traces. Final material pressure is reclosed at the exact conserved mass/volume;
the adjustment must fit the declared chemical volume/mass tolerances. Actual
mass, density, target energy, freezing energy and pressure discrepancies remain
in the result. Transport and acoustic properties are resolved at this final
temperature/pressure/composition, with no stale feed coefficients.

`ClosedGasParcel::equilibrate_adiabatic` performs the same operation atomically
on an existing rigid parcel. Its mass and volume retain their exact original
bits, and its actual stored internal energy supplies the target. The complete
bank must retain the active feed's exact caloric/transport models and source-pack
identities; it may supply initially absent products through the declared element
constraints. Source refusal, failed balances or cancellation preserve the whole
original parcel. The resulting storage and acoustic composition is frozen;
the reacting equilibrium response is returned separately. This is a bounded
thermodynamic relaxation endpoint, not a time integrator, ignition criterion,
finite-rate combustor, detonation model or experimental validation.

Focused G0/G1/G3 regressions cover an independent exothermic isomer UV root and
pressure rise, transport/acoustic consistency, reacting versus frozen heat
capacity, exact parcel mass/volume preservation, extensive scaling, source
identity, impossible energy brackets and late cancellation.



### Generated amounts and finite-rate gas reactors

`ThermochemicalGasModel::with_mole_amounts` accepts nonnegative amounts in the
bank's canonical species order and reuses the original exact-simplex and
source-preserving composition owners. Largest-amount scaling avoids overflow
while normalizing extensive inventories. Exact zeros may disappear; a positive
amount that cannot retain a positive mole fraction refuses. The resulting
immutable model retains active NASA-9, transport, water and source-pack identity.
This constructor does not change parcel mass or energy.

The L3 `fs-gas` optional `material` feature supplies the actual finite-rate
consumer through `advance_closed_parcel_adiabatic`. It integrates a declared
reversible elementary mechanism for a specified elapsed time, then uses the
existing material composition/energy replacement to publish one common gas
state. Original stored internal energy is the target; no extra reaction heat
is applied. Mass and volume retain their original bits. Source identity,
actual material balances, zero-duration identity and final cancellation are
checked before publication. Kinetics stays in L3; this crate gains no time
integrator dependency. The endpoint's acoustic composition remains frozen.

### Stable simultaneous species and additional-energy exchange

`species::network::exchange_candidate_with_heat` combines original species
streams with antisymmetric signed additional energy before one reconstruction
per affected finite gas. Heat-only nodes use the same owner. Additional energy
is separate from the existing transported partial enthalpy and may supply a
molecular energy correction; prescribing it does not establish a physical rate
law or entropy guarantee. Canonical endpoint/amount ordering and compensated
summation retain deterministic aggregation. A shared input bound covers mass
and energy records, and original gross species-withdrawal limits remain intact.

At each original node temperature T0, the thermal target uses incoming original
species enthalpy increments plus net species mass times R_i*T0, then additional
energy. The target composition is reconstructed once through the original
NASA temperature-increment owner. This cancels conserved species reference
constants before arithmetic, preserving counterflow heat when absolute formation
energies are enormous. Imported species must support T0 and the declared common
caloric bracket; source extrapolation and missing identities refuse atomically.

The local energy residual is the original-law sensible storage minus that
thermal target. The global residual sums these local closures after the paired
species/energy ledger cancels analytically. Actual total and per-species inventory
rounding is independently admitted against the intended mass changes. Displayed
absolute energy/enthalpy accounts retain their source references and are not
subtracted to judge a small caloric residual. The existing no-additional-energy
entry point delegates to this same stable owner.

Seven new regressions cover heat-only nodes, combined mass/heat, canonical order,
rollback, extreme reference shifts, zero-net-mass counterflow heat and changing
mixtures with distinct species datums. Final files passed Rustfmt parsing and
independent source/equation review. A final test typecheck was blocked before
source analysis when the executor's genuine core dependency archives disappeared;
these seven tests are not claimed as compiled or executed.


### Reciprocal free-molecular species and thermal transport

`species::effusion` samples both original Maxwellian molecular beams through a
qualified fixed thin aperture. Species i has directional molar rate
A*tau_i*p_i/sqrt(2*pi*M_i*R*T). Reciprocal transmissions are explicitly authored
or bound to original selected dimensionless constant source claims. The same
probability applies in both directions and is independent of molecular energy
and internal state. Gas caloric/transport/pack identities and active water
binding are checked directly without rebuilding a calibration composition.

Crossing energy is original partial h_i(T)-R*T/(2*M_i). The monatomic limit is
2RT per mole. Polyatomic use explicitly assumes equilibrated internal modes
independent of translational velocity. Net energy uses original relative NASA
enthalpy increments for paired counterstreams, so thermal transfer survives
large arbitrary formation-energy references. Scaled binary products retain
representable small positive rates; true underflow/overflow refuses. Species
and caller-input counts are bounded before sorting.

The authored qualification covers all mixtures/subsets of the retained original
gas catalog over its complete T/p support. Aperture diameter/minimum mean free
path <=0.1 and wall thickness/diameter <=0.01 are explicit approximation gates,
not experimental validation or rigorous error bounds. Reservoirs are stationary,
well mixed and locally Maxwellian; interface collisions, wall accommodation,
adsorption, energy filtering, bulk drift and moving-wall work are excluded.
A long accommodating pore requires its own wall-energy/transport closure.
Zero exposed area queries no inactive sources. A selected zero transmission
still admits its actual source domain at both reservoir states.

`network::transport::GasTransportPaths::effusion` retains both directional
species withdrawals and their molecular energy correction. The correction uses
the actual rounded integrated beam masses and a stable temperature difference.
It joins every other mechanism's original-inventory allowance and the ONE
shared species/energy reconstruction. No pressure-order gate is applied to
thermal transpiration. The original final entropy admission remains in force.
In the pure-gas limit p_left/sqrt(T_left)=p_right/sqrt(T_right) gives zero net
mass while unequal temperatures still exchange energy. Separate apertures may
be composed with other paths; the same exposed area must not be counted twice.

Nine face regressions and five finite-network regressions cover independent
monatomic/polyatomic energy, mass dependence, thermal balance, reversal, large
reference shifts, source/geometry/domain limits, tiny-rate representability,
finite analytical temperatures, uphill-pressure mass transfer, shared nozzle
withdrawal and atomic closure/cancellation. Final files passed actual Rustfmt
parsing and independent equation/source/API review. The final combined test
typecheck was blocked by missing genuine core dependency archives; these
fourteen new tests are not claimed as compiled or executed.

The earlier wet/DGM/pore/adsorption selection passed test typechecking at commit
642daa253c4136a65787ad22eae82c20c7344438, and nine standalone wet-law tests ran
successfully. Full-workspace validation and execution of the expanded coupled
suite, including subsequent network/effusion changes, remain unverified.

## Source-bound Prony histories and exact stress control

`visco::stress_control::StressControlledMaxwell` extends the existing
`GeneralizedMaxwell`/ `PronyState` owner to exact piecewise-linear stress
ramps, holds and zero-load recovery. The symmetric positive-definite transformed
branch system is diagonalized once through the existing Jacobi eigensolver.
Work integrates stress times strain rate; dashpot heat is integrated
independently. State admission checks the constitutive stress, finite values,
positive storage and numerical energy balance before mutation. At most 64
strictly positive branches are admitted, or an empty elastic spectrum.
Instantaneous stress jumps have zero viscous heat.

`visco::source::ResolvedPronyReference` consumes dimensioned scalar
`prony-branch-count`, `prony-equilibrium-modulus`, and indexed
`prony-modulus-N` / `prony-relaxation-time-N` claims under the original typed
query and selection policy. It verifies and retains every source receipt and
unstated uncertainty. Optional density supplies mass only when present.
Sampled curves, missing arms, wrong dimensions, invalid spectra, unowned query
coordinates and mismatched bound original geometry refuse.

The specimen starts relaxed in incremental coordinates about the fixed source
operating point. Its nominal stress/engineering strain and work use original
area/length/volume. Strain and force controls share one persistent history and
allow holds, unloading and reversal. Positive duration and authored strain,
absolute strain-rate and elapsed-time limits are mandatory. Stress-controlled
whole-ramp bounds are conservative floating-point interval sums, not certified
outward-rounded enclosures; a large step can require subdivision. Rejected
steps preserve the entire previous state and identity. Accepted responses expose
cumulative and incremental model work and dissipation, current storage, energy
residual and a history identity. Unloading can return work; dissipation remains
nonnegative. A closed strain loop can retain branch storage.

The runnable `material_viscoelastic_reference` example requires explicit
`--time-domain-reference`, bounds input/step counts and emits output only after
the complete history is admitted. Source-fitted coefficients do not independently
validate the simulated history. No unknown preload stress/energy, nonlinear
crushing, densification, permanent set, damage, acoustic absorption or thermal
evolution is inferred. The Deng foam's dynamic calibration and the printed TPU
model-versus-measured cycle-work discrepancy remain explicit in the seed catalog.
Focused G0/G1/G3 tests cover analytic creep/recovery, independent work/heat
quadrature, subdivision, source transport, mixed controls, identity and atomic
refusals; the existing strain-ramp energy regressions remain applicable.

## Prepared thermal transport over source temperature paths

`thermal::transport::PreparedThermalTransport` binds density, caller-selected
specific heat and directional scalar conductivity from one immutable material
card to a positive, explicitly typed absolute-temperature interval. Independent
claim selections and every fixed non-temperature condition remain explicit.
Canonical `ClaimSet::query_line_typed` supplies complete support and original
source knots/receipts; the prepared union grid never smooths across a breakpoint
or silently changes the source model. The union is bounded by three canonical
line budgets. Every selected source explicitly bounds temperature; value-axis
dimensions must match the typed temperature or typed fixed static coordinates.
Second source temperatures and nonzero frequency/rate coordinates refuse.

Runtime evaluation uses one binary search, no allocation, no hashing and no
database query. It interpolates rho, c and k separately before computing rho*c,
k/(rho*c) and sqrt(k)*sqrt(rho*c). The latter avoids unnecessary intermediate
product overflow. Preparation validates source support and knot values;
unrepresentable interior results still refuse without publishing a sample or
mutating the model. Splitting source segments on the union grid can change
roundoff relative to a direct query, not the underlying piecewise-affine law.

The identity binds the original card and ordered source receipts, including fixed
conditions and selection. Original uncertainty remains on the retained answers;
derived values do not acquire an invented measured claim or uncertainty bound.
No isotropy, Cp/Cv conversion, caloric integration, phase change, evolving
moisture, covariance or coupled mechanics is inferred. Re-prepare after changing
the card, fixed conditions or selected claims. Focused integration tests cover
numerical evaluation, source selection, typed refusal and the real NASA 316 pack.

## Unified gas enthalpy transport

Bulk pair transfer, bulk snapshot networks and selective pair exchange share
the original species-network caloric reconstruction. Conserved species
formation datums cancel before sensible enthalpy/internal-energy increments
are evaluated; ideal-gas pressure-flow work remains explicit. Each affected
volume is reconstructed once from its simultaneous original-inventory fluxes.
This reference invariance applies to nonreacting transport with unchanged
species source conventions; it does not remove reaction formation energy.
Displayed extensive changes retain their original reference convention, while
conservation uses local datum-cancelled energy residuals. Actual total and
per-species inventory rounding remain separately bounded. Bulk-edge and public
selective-network work limits, source conflicts, gross donor allowances,
cancellation and atomic publication remain enforced.

The sequential nozzle network retains its stable-ID Lie splitting. Its energy
audit telescopes the accepted original-caloric link residuals, since each staged
endpoint becomes the next link's initial state; actual final total and species
masses are independently checked. Per-parcel energy reports accumulate the
same accepted storage increments. The global residual is not reconstructed by
subtracting large rounded formation-energy reports. This is a balance of the
chosen first-order split trajectory, not an entropy or temporal-error bound.

## Source-relative external gas and adsorption energy

`ClosedGasParcel::exchange_species_boundary_relative` accepts each species'
complete integrated energy as `dm_i*h_i(T_original) + Q_relative_i`. The external
owner supplies `Q_relative_i` directly from its original law; it is never
recovered by subtracting rounded absolute stream totals. The fixed-volume gas
endpoint solves the original final-mixture caloric increment against
`sum(dm_i*R_i*T_original + Q_relative_i)`. Pressure-flow work is included once,
without freezing heat capacity or discarding reaction energy. This API covers
nonreacting species exchange with explicitly compatible source conventions.

The complete retained catalog supports exact individual-species depletion and
reentry. Total vacuum, unavailable or mismatched sources, unrepresentable mass
changes and unsupported caloric states refuse atomically. Actual total mass,
actual species changes versus requested transfers, and reconstructed inventory
versus intended inventory are separately bounded. Zero-net-mass thermal
exchange is permitted. Absolute reporting fields retain the original energy
datum; their rounded subtraction is not the authoritative thermal residual.
The legacy absolute-energy API remains available and cannot recover precision
already lost when a caller constructed its inputs. The partial-enthalpy-change
query evaluates the original NASA-9 source; it does not certify source-region
continuity or invent a phase/latent reference convention.

Single-site and competitive adsorption pass their original-law desorption
enthalpy correction through this relative boundary. Finite solid caloric
storage, binding heat, gross site/species inventories, equilibrium-crossing
checks and combined entropy remain under their existing owners. The combined
energy audit sums the gas and solid thermal residuals before large formation
constants enter display reports. Source alignment and independently admitted
binding enthalpy remain required; a reference shift cannot manufacture source
precision or validate an empirical adsorption model.

## Source-relative molecular phase exchange and histories

Solid/vapor surface exchange evaluates original-source enthalpy differences
before formation constants are rounded into stream totals. The existing
Maxwellian molecular energy remains `h_i(T) - R_i*T/2`; its kinetic correction
and reflected heat enter the same gas reconstruction exactly once. The shared
solid update uses the admitted sublimation latent anchor together with the
actual source-alignment residual. A compatibility tolerance does not turn the
solid/vapor enthalpy gap into an exact identity or erase its integrated energy.
Mixed solid/liquid storage contributes only bulk enthalpy minus pure-solid
enthalpy; fusion already present in the original chart is never added again.
The solid update checks the caloric change actually represented by its final
stored enthalpy and independently bounds its actual mass-change roundoff.

Specular, accommodating, equilibrium-phase and joint chamber callers all pass
source-relative species energy to the original finite gas owner. Reflected
carrier heat retains its species identity and zero net mass. Gross incoming
and outgoing inventory limits, original source domains, latent alignment and
explicit exposed-solid area remain required. Reference invariance is limited
by the precision of those independently admitted sources and stored states;
an extreme datum cannot manufacture accurate latent or binding source data.

Chamber and thermal-timeline histories accumulate only accepted original
caloric residuals and storage increments. Rejected trials and discarded coarse
solutions contribute nothing. Combined, gas-stream and solid-stream defects
remain subject to the declared whole-history energy tolerance; actual total
and species mass checks remain independent. Thermal heat and incoming birth
energy retain their original owners and explicit datum offsets. Absolute
storage/stream reports retain reference energies, while their rounded
subtraction is not used to reconstruct small heat or to accept the trajectory.
These changes preserve the existing explicit kinetics, splitting, source
support and molecular assumptions; they do not introduce a new evaporation,
reaction, accommodation or real-gas model.

## Directional passive plane-stress viscosity

`PlaneStressViscosity` admits positive-semidefinite directional loss while
retaining all six source coefficients, units, identity and frequency support.
`factor_passive_resultant` normalizes the three-component work form and reuses
the tensor Prony PSD owner, including its `64*EPSILON*norm` zero-eigenvalue
roundoff allowance, with reconstruction checked at `1e-12`. Rank-deficient
stress and power use the retained Gram factor; no positive loss floor is
introduced. The established positive-definite path retains its arithmetic.
The native material tests cover axial and coupled rank-one nullspaces, rotated
work, thickness scaling, coefficient scales `1e-200..1e200`, and active-tensor
refusals. This admits declared reduced viscosity; it does not infer a bulk
constitutive law from absent data.

## Direct volumetric transport and free-expansion density

`thermal::transport::volumetric::ResolvedVolumetricThermalTransport` resolves
positive volumetric heat capacity and scalar static conductivity from one card
at an exact typed condition, retaining original answers, uncertainty and receipts.
It shares static-temperature/source-axis admission and factorized effusivity
arithmetic with prepared transport. It derives k/C_vol and sqrt(k)*sqrt(C_vol)
without inventing density, mass-specific heat or independent statistical evidence.
Refresh pins the original claims and requires the original card and coordinate
schema; a failed update leaves the old snapshot intact. A point observation does
not acquire interval support or become a caloric integral.

`thermal::transport::expansion::PreparedFreeExpansionDensity` admits a positive
reference density and continuously supported engineering free length changes.
The caller explicitly declares isotropy or fixed orthogonal principal directions,
constant mass and unconstrained free expansion. Each length ratio is normalized
at the density reference; the finite volume ratio is their product and current
density is reference density divided by that product. The model retains original
source lines and density evidence, fixed context, basis/frame, selections and
nonzero source datum residuals. CTE, infinitesimal tensor strain and unknown
temperature support cannot substitute for the declared source roles.

The union grid includes the reference temperature. Evaluation uses one binary
search and no allocation, hashing or source query. It derives stretches, volume
and density from interpolated engineering strains, rather than interpolating the
derived density. Positive stretches, volume and density must remain representable;
an extreme interior refusal publishes no sample and leaves the path unchanged.
There is no inferred mass transfer, constrained strain, phase law, evolving frame,
cross-source specimen match or propagated uncertainty.


## Saturated-water boiling and condensation

`phase::vaporization::IsobaricWaterVaporization` binds the original
`fs-thermochem::water` IAPWS source to a fixed-pressure liquid/vapor enthalpy
plateau. It evolves actual vapor mass fraction, additive specific volume,
bulk density, internal energy and entropy; vapor is not represented as liquid.
The fixed total mass may expand against an external pressure load. The
returned pressure is the original source evaluation within the caller's
pressure-inversion tolerance, not a patched target pressure.

Finite heat increments independently audit Q - Delta H and
Q - Delta U - p Delta V against one explicit absolute joule budget.
Pressure work, internal energy and entropy changes use the actual represented
vapor-mass change and the original source's direct phase gaps. Exact requested
heat to either pure endpoint selects that source endpoint; excess heat still
refuses rather than extrapolating into subcooling or superheating. Critical
coalescence, unrepresentable increments and cancellation also refuse.

The existing shared-thermal `ThermalLaw` consumer uses T(h)=T_sat and the exact
plateau entropy integral Delta h/T_sat. Its geometry is reference heat-exchange
geometry, not a rigid-vessel constraint or a dynamically updated fluid mesh.
Real timeline tests cover heating/cooling with ambient loss, finite contact
between a steam condenser and a lower-pressure boiler, positive heat-transfer
entropy production, complete phase depletion/reentry and atomic refusal.
`cargo run -p fs-material --example water_boiling` runs a heating/cooling cycle
with phase inventory, volume, expansion work and independent energy checks.

### Adiabatic pressure-drop flashing

`WaterVaporizationState::flash_to_pressure` resolves the equilibrium endpoint
of a steady adiabatic throttle with no shaft work and negligible kinetic and
potential energy changes. Specific enthalpy is retained bitwise. The original
IAPWS source at the lower pressure determines actual liquid/vapor inventory,
volume, internal energy and entropy. Independent source Delta u + Delta(pv)
and signed entropy-generation checks use explicit caller budgets; p_out Delta v
is not substituted for the inlet/outlet flow-work difference.

The report retains the actual pressure drop, original pressure-inversion
residual, flow-work change and signed entropy generation. Pressure increases,
critical collapse and enthalpies requiring subcooling or superheating refuse;
no off-saturation source or valve flow-rate model is inferred. The returned
state's `law()` reuses its exact retained chart for further heating or cooling.
Tests include vapor creation from saturated 450 K liquid throttled to
101325 Pa, a heating/flash/condensation chain, conservation and source-domain,
accuracy and cancellation refusals.


## Original sources at native shear and electrical operating coordinates

`phase::transport::resolve_phase_transport_state_with_shear_rate` extends the
existing phase-first material resolution with one explicitly typed inverse-time
coordinate owned by the native flow. Temperature and phase still come from the
original enthalpy chart. Liquid viscosity and conductivity may depend on actual
shear; density and caloric storage may not acquire an unowned rate dependence.
Solid endpoints retain the existing omission of liquid-only properties.

`conductor::response::ScalarConductorResponse` retains one complete pinned scalar
conductor bundle and explicitly declares at most one electric-field magnitude
and one current-density magnitude coordinate. Native spatial consumers query
the original source at actual |E| and |J|, including all original requirements,
claim pins, coordinate quantities and external context. Conductivity and
resistivity remain distinct source bases. Every original electrical response
coordinate must be owned; frequency, strain and other unowned responses refuse.

The shared multiple-owned-coordinate predicate preserves the existing static
external-coordinate rules; it provides no evolution law itself. The electrical
response adapter is an isotropic, passive, memoryless rung, with no tensor, Hall,
hysteresis or field-storage claim. Ordinary phase metadata is retained; combined
phase-fraction/E/J source admission is not established merely by that metadata.
Actual field/current closure and finite-heat ownership belong to their original
native consumers, not this material adapter.

## Source-bound complete tensor relaxation

`visco_tensor::source::TensorPronyProperties` binds the original resolved
material state to a complete 3D Maxwell law. Equilibrium stiffness comes from
the canonical isotropic or full orthotropic elastic properties. Each source
branch names all 21 symmetric upper-triangle Mandel entries and its relaxation
time; zero entries remain explicit data. A positive small-strain limit is also
a selected source claim. Shared names resolve once, incompatible quantity roles
refuse, and the original material state retains all selected usage receipts.

Every coefficient must be finite and constant at the admitted external state.
The law retains the intersection of all selected frequency domains, including
equilibrium, strain, time and branch data. The original tensor constructor owns
PSD admission and its roundoff policy. Missing data, frequency-sampled elastic
coefficients offered as equilibrium springs, invalid complete elasticity,
nonpositive times or strain limits, and more than 64 source branches refuse.
This binding does not condense plane stress, invent a full tensor from in-plane
coefficients, evolve temperature/moisture, or establish measured-material
authority. Five focused source tests cover full coefficients, shear convention,
orthotropy, shared names, source bands and refusals; native execution of these
new tests remains unverified in the resource-constrained authoring environment.

## Continuous volumetric storage and expansion-aware transport

`thermal::transport::volumetric::PreparedVolumetricThermalTransport` prepares
positive, continuously supported volumetric capacity C_vol(T) and directional
scalar k(T) on one fixed typed source context. It retains both original source
lines, selections, receipts, uncertainties and the immutable card identity.
Every source breakpoint enters the union grid. Point evaluation uses one binary
search, interpolates the independent inputs, then derives k/C_vol and effusivity;
no density, mass-specific heat or covariance is invented.

`storage_change_j_m3(T0,T1)` integrates the piecewise-linear C_vol chart directly
over crossed segments. `temperature_after_storage_change(T0,delta_storage)`
inverts that same monotone relation for heating and cooling using bounded knot
traversal and a scaled quadratic solve. Both use the caller's fixed volume and
caloric convention. They do not automatically represent internal energy,
enthalpy, pressure work, changing volume, latent heat or phase change. Numerical
exactness is for the selected affine model, subject to roundoff, not a measured
interpolation-error claim. Nonrepresentable updates and unsupported temperatures
refuse without mutation. A representable local inverse remains available when
the unused whole-interval integral overflows. Successful runtime calls allocate
nothing and perform no material queries.

Thermal path admission distinguishes an explicitly typed, fixed
`TemperatureDifference` test condition from a second absolute temperature.
Unknown or extra absolute-temperature response coordinates remain unsupported.

`thermal::transport::expanding::PreparedExpandingThermalTransport` composes the
original `PreparedFreeExpansionDensity` with selected Cp(T) and k(T), requiring
the same card, exact endpoints, axis and fixed context. It keeps the original
density model and its volume/reference assumptions intact. At runtime it
evaluates that model and the Cp/k grid, then derives transport; it never replaces
free-expansion density with a linear endpoint fit or a new measured claim.
Preparation checks density-only knots as well as thermal knots. Evaluation uses
at most two binary searches without allocation, hashing, queries or mutation.

Focused integration tests are `material_expanding_transport`,
`material_volumetric_curves` and `fsri_plastic_thermal_profiles` in xtask.
They cover analytic piecewise storage/inverse identities, source and numerical
refusals, original NIST expansion/Cp/k with an explicitly assumed reference
density, and real FSRI foam/plastic heating and cooling at fixed volume.

## Caloric source provenance survives the original phase-curve handoff

`EquilibriumEnthalpyPhaseCurve::source_receipts` retains the exact original
source usages from both `ResolvedHeatCapacityCurve` and
`ResolvedMeltingCurve`. Cp, density, transition temperature and latent-heat
selection/context remain available when a thermal consumer receives only the
plain curve. Existing wrapper `receipts()` methods delegate to this same
inventory. Source identity hashing and caloric interpolation are unchanged:
the receipts already contributed to the sourced chart identity.

This closes a loss of source ownership when another subsystem evolves a
coordinate held fixed while the caloric temperature line was constructed.
Consumers must inspect the actual selected claims before treating that chart as
independent of the evolving coordinate. Explicitly authored h(T,phase) charts
carry no resolved-source receipts and retain their declared law; an empty
receipt list is not empirical material validation. Tests retain Cp/density and
melting receipts through the original curve and its clone.

## IF97 water heating, boiling and superheat

`phase::water::IsobaricWater` gives one original-source, constant-pressure
enthalpy chart spanning compressed liquid, liquid/vapor coexistence and
superheated steam. It consumes only `fs-thermochem::water::if97`; original h, u,
s and v share that source's reference and are never mixed with SR1-86.
Temperature bounds must bracket saturation, remain within 273.16–1073.15 K and
admit the Region 1/2 coexistence join only through 623.15 K. This is not an ice,
Region 3, supercritical, metastable, kinetic phase-change or rigid-vessel model.

Original-source p,h and p,s inversions use safeguarded Newton/bisection, caller
residual and iteration limits, and cancellation. Coexistence resolves quality
analytically at the fixed saturation temperature. Actual source endpoints and
residuals are returned without shifting h or s. Single-phase states report no
vapor quality; mass fractions, density and specific volume remain available in
every phase. Temperature alone at saturation refuses as ambiguous.

Construction proves positive Cp over both complete sensible branches using the
source's bounded interval arithmetic. The retained global inverse slope admits
stable use of the existing `fs-thermal` timeline across latent plateaus; it is
not a sampled or locally guessed heat-capacity bound. Chart construction is
separate from repeated state queries.

The original IF97 regional fits have a finite saturation matching error:
delta_s=(h_v-h_l)/Ts-(s_v-s_l). Admission requires an explicit bound on its
magnitude. Raw source entropy and this defect remain visible. The thermal law
uses the continuous pressure-specific primitive S*=s+chi*delta_s, where chi is
zero in liquid, vapor quality in coexistence and one in vapor. It integrates
dh/T on this isobar, and must never be substituted for absolute source entropy
in a pressure-changing device. The correction is not physical entropy
generation or an added heat source.

`advance_heat` retains fixed mass and reports both phase inventories, physical
volume change, delta_h, delta_u and outward p*delta_v work. Original-source
stable increments and actual rounded endpoint differences independently satisfy
the caller's unchanged absolute joule budget. Nonzero heat or phase inventory
lost to representability refuses. The shared thermal geometry remains reference
heat-exchange geometry, while physical expansion is explicitly reported.

Six native material tests exercise p,h/p,s source round trips; complete
heating/condensation and phase depletion; explicit entropy-seam accounting;
domain/work/cancellation refusals; the real shared thermal history through all
three regimes; and finite steam/liquid thermal contact with conserved total
enthalpy and positive modeled entropy production. The runnable
`cargo run -p fs-material --example water_heating` heats and cools 50 g from
300 K through boiling to 700 K and back, independently reconciling every
accepted thermal interval with material calorimetry. This establishes a bounded
liquid/vapor consumer for bead `frankensim-7sga6.2.7`, not its complete ice and
all-regime scope.

## IF97 adiabatic pressure devices

`phase::water::flow` supplies steady throttles, liquid pumps, vapor compressors
and steam expanders on the same original IF97 charts. All quantities are per
kilogram of throughput. Pressure directions are explicit; kinetic/potential
changes, transient machine dynamics and flow-rate laws are outside this owner.

A throttle resolves h_out=h_in at a lower pressure, admitting compressed liquid,
wet or superheated outlets where the prepared source chart supports them.
A shaft device first resolves the original isentropic outlet from
s_out=s_in. Positive shaft work is outward: an expander gives
w=eta*(h_in-h_isentropic), while a compressor/pump gives
w=(h_in-h_isentropic)/eta. Efficiency is an explicit input in (0,1], with no
invented empirical default. The actual outlet follows from h_out=h_in-w.
For eta=1 the admitted isentropic source state is reused directly, avoiding an
unnecessary second inverse and cancellation in reconstructed target enthalpy.

Each result independently admits delta_h+w and delta_u+delta(p*v)+w against
the caller's unchanged numerical energy budget. Flow work is p_out*v_out minus
p_in*v_in, not p_out*delta_v, and is not added again to shaft work. Cross-pressure
entropy always uses original source s. Its signed change and each chart's
saturation matching defect are reported separately; the fixed-pressure thermal
entropy correction is never used as a global reference. Source-domain failures,
lost nonzero work or efficiency loss, entropy decrease outside the caller's
roundoff allowance, budget failures and cancellation publish no partial result.

Five focused native device tests pass alongside the six material tests. They
cover liquid/steam throttles, pump work against independent integration of
v*dp along a raw-source isentrope, superheated compression, wet turbine outlets,
a complete Rankine cycle, and physical/domain/cancellation refusals.
`cargo run -p fs-material --example water_rankine` exercises a 10 kPa
condenser, a 4 MPa boiler at 773.15 K, pump efficiency 0.80 and turbine efficiency
0.85. It reports states and all four process balances, including a separate
internal-energy/flow-work cycle check. Its approximately 31.22% efficiency is
the result for those illustrative inputs, not a calibrated plant prediction.
Heat-exchanger reservoir entropy and pressure losses are not supplied.

## Shared tensor strain admission in all material-history interfaces

`TensorProny::validate_mandel_strain` checks the original finite small-strain
bound in the Frobenius norm, with Mandel shear normalization. The original
Mandel `step`, tensor-shear `history().advance` and projected strain observer
use this one admission rule before constitutive history can change. Source
limits therefore survive a change of consumer or coordinate representation.
The focused regression checks shear-stress parity and identical refusal without
mutation in both time-domain interfaces. This admits supplied endpoints only;
it does not enclose strain extrema between time samples or quadrature points.

## Volumetric storage on the shared thermal timeline

`thermal::transport::volumetric::timeline::FixedVolumeThermalLaw` borrows the
existing prepared volumetric source path and implements the original
`fs_thermal::enthalpy::ThermalLaw`. Its explicit reference mass/volume ratio
rho0 and zero-storage datum give h(T) = integral C_vol(T) dT / rho0. A compatible
cell therefore stores m*h = V0*integral C_vol(T) dT. rho0 is an authored scenario
basis, not a new measured density, temperature-dependent density or source Cp
claim. Mass and volume remain fixed. Callers must select a compatible caloric
constraint; this does not infer pressure work, phase change or moisture exchange.

Preparation retains the original source identity, includes the mass basis and
datum in the adapter identity, checks representable storage bounds and computes
an upward-rounded global temperature-slope bound from the minimum affine
capacity. Runtime storage/inversion reuse the original exact segment integrals.
`admit_cell` checks actual mass versus reference volume and supported storage;
its 64-epsilon comparison permits arithmetic roundoff only.

`FixedVolumeConductionLink::try_aligned` binds an existing ordered interface
to the two original borrowed caloric laws. Each evaluation verifies the active
link geometry, actual cell mass/volume bases and original material-law objects.
It evaluates each cell's sourced k at its current temperature, then uses the
existing series/contact conductance owner. Distances, transfer area, contact
resistance and alignment are explicit caller inputs. This is a uniform-cell,
cell-centre conductivity approximation, not an exact nonuniform half-cell
integral. Return the value through the existing implicit
`StepCoefficients::conductances_w_k`; the original timeline owns history,
sources, integration and local/global energy acceptance. Evaluation allocates
no source tables, queries no card and computes no content hashes.

Four focused G1/G3 tests in `xtask/tests/material_volumetric_timeline.rs`
cover reference-basis/datum invariance, an independently solved nonlinear
two-cell balance, original-law/geometry/domain refusals, and actual compiled
FSRI foam and density-free FRP heating on the existing timeline. The explicit
authored mass basis cancels from those fixed-volume temperature results.
An entropy interpretation is available only through the explicit opt-in below.
No experimental qualification or full spatial FEM result is added.

## Sensible entropy from the retained volumetric path (2026-10-28)

`PreparedVolumetricThermalTransport::sensible_entropy_change_j_m3_k`
integrates the retained piecewise-affine capacity as
`integral C_vol(T)/T dT` in J/(m3 K). Calling it explicitly selects a
reversible sensible-caloric interpretation at fixed work/composition
coordinates and the source reference-volume convention. An apparent
experimental capacity does not automatically become thermodynamic Cv.
No absolute entropy datum, latent heat, pressure work, chemistry, moisture
exchange or hysteresis is supplied.

Each crossed segment uses positive endpoint weights for its analytic
logarithmic integral. A bounded small-interval series avoids cancellation;
compensated accumulation retains the original source breaks. Reversal is
antisymmetric. Unsupported temperatures and unrepresentable arithmetic
refuse without changing the source path. Runtime traversal is bounded by
the original knots, allocates nothing, queries no card and hashes nothing.
The platform logarithm remains subject to floating-point roundoff; no
cross-ISA transcendental or measurement-uncertainty guarantee is added.

`FixedVolumeThermalLaw::with_sensible_entropy` opts the SAME borrowed source
law into the existing `ThermalLaw::specific_entropy_change_j_kg_k` account.
Default bindings still report entropy unavailable. Opt-in changes the
binding identity exactly once, retaining the source receipts, authored
reference mass/volume ratio and storage datum. The caller must justify
the capacity's caloric/work constraint. In exact arithmetic,
`m * delta_s = V0 * integral C_vol(T)/T dT`; the authored reference mass basis
cancels.

The adapter uses the actual stored enthalpy increment and a capacity-weighted
mean reciprocal temperature. Distinct enthalpy coordinates that round to
one temperature retain their nonzero `delta_h/T` increment. Existing
storage/inversion admission remains active: even finite endpoint
coordinates may refuse if their difference or another required intermediate
is unrepresentable. No broad all-finite-input success guarantee is made.

The original `fs_thermal` explicit and backward-Euler owners continue to
own heat flux, state acceptance, reservoir entropy, conductive generation,
unclassified source watts and temporal defects. The adapter neither
reclassifies a source as boundary heat nor counts numerical entropy as
material production. Three additional tests in
`xtask/tests/material_volumetric_timeline.rs` cover independent analytic
integrals, distinct source breaks, adjacent-temperature/enthalpy increments,
reference-basis/datum invariance, default refusal and opt-in identity, an
isolated implicit conduction balance, and actual compiled FSRI PP/mattress
foam heating. All seven timeline tests pass, including the prior four
energy/storage regressions.

## Integrated conductivity layers on the volumetric timeline (2026-10-29)

`FixedVolumeConductionLink::into_integrated_layers` explicitly selects the
existing quasi-steady, one-dimensional layer model. The returned
`FixedVolumeLayerConduction` copies each original conductivity line once,
preserving its knots and temperature support. The borrowed caloric paths retain
their source receipts, fixed reference mass basis, datum and entropy selection.
It integrates k(T) through both legs using the shared `TwoLayerConduction`
kernel; the kernel's public API and numerical equations remain unchanged.

Every evaluation admits the actual active link, ordered endpoints, transfer
area, original material-law objects, mass/volume bases and stored-energy
coordinates. Both conductivity charts must cover the full endpoint-temperature
interval. Positive lengths and area, nonnegative area-specific contact
resistance, an iteration budget and numerical tolerances are explicit inputs.
The caller declares homogeneous, internally source-free legs and alignment of
the selected scalar k with heat flow. This does not establish isotropy or source
qualification for unmeasured conditions.

The response carries one signed heat transfer, a secant conductance, both
contact temperatures, endpoint derivatives and continuity diagnostics. Feed
its conductance into the original implicit `StepCoefficients` callback on every
candidate state. The secant conductance is distinct from the returned endpoint
Newton derivatives. The existing timeline owns storage, acceptance, sources,
energy conservation and optional sensible entropy. No layer storage or second
time integration is introduced. Successful runtime evaluation allocates no
memory, performs no database query and computes no content hash; work remains
bounded by the retained knots and the existing solver's iteration budget.

Three additional tests in `xtask/tests/material_volumetric_timeline.rs` cover
independent affine-conductivity integrals and derivatives, contact between two
different materials, source/geometry/domain refusals, cancellation and budget
exhaustion, and compiled FSRI polypropylene on the original implicit energy and
entropy account. All ten tests in that file pass in the focused native build.
These are bounded constitutive and timeline results; spatial convergence,
experimental qualification, phase change and dynamic moisture transport remain
outside this model.


## Exact equal-time material-history coordinates

`visco_tensor::coordinates::RelaxationCoordinates` groups displacement-driven
history rows only by bit-identical positive relaxation time. For each original
factor `A`, the existing Householder QR owner supplies `A = Q R`; complete thin
columns, including rank-deficient and zero directions, remain. With relaxed
`z = A q` or unrelaxed `z = 0` initialization, `z = Q y` is invariant and
`y' = (R q - y) / tau` preserves the original force, storage and dissipated
heat. Nearby poles are not merged and no mechanical mode or material coefficient
is fitted. Families no larger than the displacement space keep original rows.

Column scaling, orthogonality and source-entry reconstruction checks refuse
unrepresentable transforms. Explicit state/work budgets and cancellation cause
refusal rather than truncation. Immutable `RelaxationMemoryMap` slices share one
complete inverse map and recover every original material-point history for its
existing constitutive and source-domain checks. Independently authored original
histories outside this invariant subspace cannot be silently projected into it.
This is a coordinate representation of the declared original material dynamics,
not permission to lower spatial, modal or temporal convergence requirements.

## Joint phase and electrical source coordinates

`ScalarConductorResponse::new_with_liquid_mass_fraction_axis` binds the original enthalpy-derived liquid mass fraction and explicitly owned actual E/J magnitudes together. Complete source requirements, selected claims, typed coordinates and support survive phase, temperature and field refresh. Public static phase admission and generic material refresh remain strict. The existing source data supplies its admitted law; no universal phase mixture rule or higher-dimensional table is invented. Focused source regressions cover retained roles/pins, actual-field queries, phase support and unowned-coordinate refusal.

## Free dimensional change without a density requirement (2026-10-29)

`thermal::transport::expansion::PreparedFreeExpansion` consumes original finite
engineering length-change sources with an explicit reference temperature. The
new `FreeExpansionGeometryQuery` selects either an authored isotropic
interpretation or three fixed orthogonal principal directions. No density or
mass input is needed to compute the supported dimensional change.

At temperature T, each principal stretch is
`(1 + epsilon_i(T)) / (1 + epsilon_i(T_ref))`; the relative volume is their
product. Source breakpoints, nonzero fitted datum residuals, complete typed
query context and original receipts remain intact. The reference is retained as
an exact union knot. Only the independent source strains interpolate; stretches
and volume are derived from them. The identity includes the reference, frame,
original card and source receipts. Evaluation uses one binary search and
allocates nothing, queries no card and computes no hash.

`PreparedFreeExpansionDensity::try_from_expansion` optionally attaches a density
answer at that same reference and computes `rho_ref / J` under the caller's
constant-mass assumption. It requires the original immutable card and full
reference context, and reuses the prepared strain lines. The existing density
constructor, public sample accessors and v1 density identity byte order are
preserved. Existing expanding thermal transport continues to consume this
composition without separately interpolating density products.

Four additional tests in `xtask/tests/material_expansion_density.rs` exercise
independent anisotropic geometry, reference-dependent identities, optional
density composition, source/domain/frame/length refusals, and compiled NIST 304
dimensional change without an inserted density claim. All eight tests in that
suite and three existing expanding-transport tests pass in the focused native
build. Thermal coefficients retain their distinct dimensions and are refused
by this finite-strain API; mean and differential coefficients are not inferred
from one another. Constrained mechanics, changing directions, phase transitions
and uncertainty propagation remain separate model requirements.

## Rigid water vessels with evolving pressure

`phase::water::isochoric::IsochoricWater` closes equilibrium water at fixed
specific volume with internal energy as the storage coordinate. Original IF97
Regions 1, 2 and 4 supply compressed liquid, liquid/vapor coexistence and
superheated vapor over a declared 273.16–623.15 K window, with declared pressure
bounds through 100 MPa. Saturation quality follows the actual phase volumes.
Off saturation, a safeguarded pressure inverse uses the original analytic
volume derivative. A second safeguarded inverse resolves u(v,T); its derivative
is the actual isochoric heat capacity, including phase conversion and the
original Region 4 pressure derivative. It does not substitute the approximate
Clapeyron identity for that independently fitted source equation.

Every heat step keeps mass and vessel volume fixed, admits Q=delta_U and
Q=delta_H-V*delta_p independently, and reports zero boundary work. The pressure
part of enthalpy is not an additional heat source or wall work. Original-source
endpoint volumes, volume residuals, phase inventories, entropy and saturation
Gibbs mismatch remain visible; caloric values are not patched to fit budgets.
Nonzero heat or phase mass lost to representation refuses. Controls bound
volume/energy residuals, pressure/temperature iterations and source windows;
cancellation or refusal leaves the accepted state unchanged.

The temperature window is a bounded inverse search interval, not a certified
global Cv enclosure. Unsupported intermediate source queries refuse. This
owner does not claim Region 3, ice, mass flow, wall compliance, phase kinetics
or a spatially resolved vessel.

Six focused native tests pass: original p/v/u round trips; analytic two-phase
Cv against independent finite differences; sealed evaporation/condensation;
liquid pressurization and the cold-water density anomaly; heat subdivision and
small increments; domain, work, energy and cancellation refusals. The runnable
`cargo run -p fs-material --example water_rigid_vessel` heats 50 g in 5 L from
300 to 600 K and back. Pressure rises from approximately 3.54 kPa to 2.61 MPa,
with maximum complete-history energy discrepancy 5.27e-8 J and zero wall work.
These are numerical balance checks, not source-fit uncertainty estimates.
This advances the bounded phase-transition consumer in
`frankensim-7sga6.2.7`.

## Dense-water rigid-vessel calorimetry

`phase::water::dense::DenseIsochore` connects the original IF97 Region 3 source
to fixed-mass, fixed-volume heat evolution. The constructor retains one
`DenseIsochoricInterval`: original regional admission, stable derivatives and
positive Cv are established over the entire requested temperature path, with
explicit work limits and cancellation. Its enclosure and dT/du bound are reused
across state queries and heat steps; no proof is repeated inside the caloric
inverse.

A safeguarded Newton/bisection inverse uses original u(rho,T) and Cv to resolve
a heat-driven endpoint. Pressure follows the same Helmholtz equation. Original
stable increments separately provide delta_u, delta_h, delta_p and delta_s.
Accepted heat must satisfy both Q=delta_U and Q=delta_H-V*delta_p using those
increments and independently using the actual rounded source endpoints.
V*delta_p is the pressure contribution to enthalpy, while rigid-wall work is
exactly zero. No heat, entropy, caloric reference or nonzero increment is silently
clipped or reconciled. Unrepresentable heat, failed energy admission, unsupported
source states, work exhaustion and cancellation return no partial state.

This is a single-phase Region 3 equilibrium owner. It does not automatically
stitch across the Region 1/2/3 boundaries, provide high-temperature saturated
pairs, resolve a moving interface, or model vessel compliance and strength.
The separate `isochoric` owner retains the lower-temperature liquid/vapor
coexistence model. Neither substitutes its internal-energy coordinate for the
existing enthalpy timeline's declared caloric convention.

Four focused native tests pass: source u/T round trips and retained Cv bounds;
high-pressure heating/cooling with separate pressure enthalpy and zero wall
work; independent Cv/Cv-over-T integration plus subdivision and small heat;
and domain, work, energy, chart-identity and late-cancellation refusals.
The complete original water module compiles with both new volume owners, and
both vessel examples plus the existing Rankine example execute successfully.

`cargo run -p fs-material --example water_dense_vessel` heats 0.5 kg in 1 L
from 650 to approximately 750 K and cools it back. Pressure rises from about
25.584 to 78.310 MPa. Heating supplies 144.903 kJ to internal energy; enthalpy
increases by 197.629 kJ, including 52.726 kJ of V*delta_p and zero wall work.
The whole-path preparation uses 113 interval evaluations for its 649–751 K
source window, reused over 40 heat steps. Maximum complete-history numerical
energy discrepancy is 1.47e-6 J. Source uncertainty, reservoir entropy and a
container-strength prediction are outside that numerical balance claim.

## Compositional sensible heat storage (2026-10-29)

`thermal::transport::storage::PreparedHeatStorage` prepares one or more explicit
constituents at a shared temperature. `HeatStorageAmount::MassKg` multiplies a
source specific capacity; `ReferenceVolumeM3` multiplies native volumetric
capacity. Neither requires conductivity or an inferred density. Each component
retains its immutable card, exact claim selection, typed fixed conditions,
original source line, uncertainty and physical amount. The caller declares
compatible caloric/work constraints and the shared-temperature approximation.

Preparation intersects the selected temperature supports and retains every
source breakpoint inside that intersection. Constituent IDs must be unique;
sorting them before compensated summation and hashing makes input permutation
irrelevant. Capacity amounts must be positive, finite and dimensionally matched.
Preparation admits at most 64 constituents and 262,144 union knots. Runtime
capacity needs one binary search; energy, inverse temperature and optional
sensible entropy use the shared positive-affine caloric kernel. Successful
evaluation allocates nothing and performs no source query, content hashing or
numerical root iteration. Integrals traverse relative temperature intervals to
retain small changes; finite representability limits can still refuse.

`HeatStorageThermalLaw` binds this total storage to the existing thermal
timeline's specific-energy coordinate with an explicit total mass and datum.
Mass-only assemblies must match their known total mass; mixed assemblies require
additional positive mass for their volume-based constituents. Call `admit_cell`
before inserting a birth event, since the generic `ThermalLaw` evaluation cannot
inspect cell mass. The timeline continues to own power, evolution, cancellation
and energy acceptance. Sensible entropy is opt-in and changes the law identity.

Four focused assembly tests exercise independent mass/volume integrals,
inversion, entropy, amount scaling, source choice, permutation, domain refusals
and actual thermal timelines. Compiled MDF plus acrylic/mineral countertop works
below MDF's conductivity support; beryllium heats from its Cp source alone.
The three existing volumetric-curve tests and ten existing volumetric-timeline
tests also pass with the extracted shared caloric kernel. This addition supplies
no conductivity homogenization, phase or chemical mixing law, thermal-expansion
work, moisture evolution, absolute entropy, propagated uncertainty or automatic
justification of a uniform-temperature body.

## Dense isobaric heating, boiling and supercritical water

`phase::water::pressure::DenseIsobaricWater` binds the original Region 3
pressure/density solver to a fixed-pressure enthalpy owner. Explicit temperature
windows stay inside Region 3. Subcritical windows may retain the original
high-temperature Region 3/4 saturation pair; supercritical states have no
liquid/vapor quality or invented phase inventories. Original h, u, s and phase
volumes remain unshifted. Vapor quality is a mass fraction, and mixture volume
uses the constituent specific volumes.

Temperature, specific-enthalpy and raw-specific-entropy queries share the
retained chart. Caloric inverses use safeguarded Newton/bisection with original
Cp or Cp/T, explicit residual tolerances, source work counts, and 1..=128
outer iterations. Density inversion has its own bounded controls. Authored
caloric endpoints and latent endpoints are retained exactly. This search
window is not a global Cp certificate: unresolved interior source queries,
including near-critical resolution limits, can refuse. No `ThermalLaw`
stability bound or automatic Region 1/2/3 temperature-path stitching is claimed.

Fixed-mass heat deposition evolves H, with the independent physical balance
Q=Delta U+p_requested*Delta V. The requested mechanical pressure and actual
Helmholtz endpoint pressures remain distinct. Reports expose source pressure
residuals, endpoint h-u-p_requested*v defects, expansion work, actual phase
masses and volumes, and raw source entropy. The isobar-local entropy primitive
uses the explicit saturation defect delta_s=(h_v-h_l)/Ts-(s_v-s_l);
its correction is not an absolute entropy for pressure-changing equipment.

Stable bivariate source differences and latent lever-rule changes avoid
subtracting large reference energies for small heating steps. Acceptance gates
both stable and independently represented endpoint versions of Delta H-Q and
Delta U+p*Delta V-Q. Nonzero heat or phase mass lost to representation,
unsupported states, chart mismatch, exhausted work and early/late cancellation
return no partial state or source debit. A numerical residual budget is not a
physical source uncertainty or a license to replace source pressure.

Six focused native tests pass: 25 MPa dense-liquid/supercritical h and raw-s
round trips; 20 MPa high-temperature boiling and condensation with distinct
mass/volume fractions and one latent account; independent Cp and Cp/T integration
through the pseudocritical peak; subdivided heating/cooling, expansion work and
small/zero heat; coarse pressure closure rejected by the physical joule budget;
and domain, ambiguity, work, seam, chart, representability and cancellation
refusals. No transport coefficient, spatial boiling, moving interface, vessel
mechanics or reservoir entropy model is supplied by this equilibrium owner.

## Pressure equipment across original water source regions

`phase::water::flow::across` connects existing Region 1/2/4 Gibbs charts and the
new Region 3 Helmholtz pressure chart. Typed state/chart enums retain the exact
source owner and reuse the existing pressure-device controls and balance
report. Steady adiabatic throttles preserve original h; compressors and expanders
resolve original s at the outlet pressure, then apply the caller's explicit
isentropic efficiency. No pressure-local entropy correction is used across
pressures, and no caloric datum is shifted to manufacture regional agreement.

Every accepted device independently closes Delta h+W and
Delta u+Delta(p*v)+W, with outward shaft work positive and flow work included
once. Original h-u-p_requested*v defects at the inlet, outlet and isentropic
reference are individually limited by the energy budget so pressure errors
cannot cancel invisibly. Reports retain source entropy production, numerical
entropy residual, source-junction defects, density-source evaluation counts and
caloric iterations. Reversible devices reuse the isentropic endpoint. Numerical,
chart-domain, direction, efficiency, work and cancellation refusals return no
partial device state. The production dispatch and solves allocate nothing.

Four focused tests pass: a dense-liquid throttle reaching wet low-pressure
steam; supercritical expansion with raw entropy and declared efficiency;
hot-liquid compression across the Gibbs/dense interface with reversible return;
and invalid direction/efficiency, pressure closure, unsupported outlet, work and
cancellation refusals. The combined production water library also passes the six
new dense-calorimetry tests and all five existing pressure-device/Rankine tests:
15 integrated native tests in total.

`cargo run -p fs-material --example water_high_pressure_steam` runs a partial
steady steam train: 20 MPa water heats from 630 to 645 K through its
638.895911546 K saturation plateau, expands through an explicitly illustrative
85% efficient turbine to 10 kPa, and condenses to saturated liquid. Boiler heat
is 851.285931773 kJ/kg, including 584.286587172 kJ/kg of latent enthalpy.
Shaft output is 788.761576280 kJ/kg; the wet turbine outlet has vapor quality
0.659460798. Maximum independent stage energy residual is 1.08e-8 J/kg and
the complete open-train residual is 2.10e-9 J/kg on the native validation run.
These are numerical balances, not source or equipment measurement accuracy.

The inlet already contains substantial thermal energy. Feed-pump/preheater
return and reservoir temperatures are not modeled, so this example supplies
neither a closed-cycle efficiency nor a reservoir-entropy claim. Equipment
endpoint connections do not imply automatic temperature-path stitching of
Region 1/2/3, transient machine dynamics, pressure-drop/flow-rate laws, phase
kinetics, or resolved two-phase flow.

## Ordered material layer stacks (2026-10-29)

`thermal::transport::stack::PreparedLayerStackConduction` binds 1–64 physical
layers to their original scalar conductivity sources, thicknesses, common area
and explicit area-specific contact resistances. Each layer retains its stable
ID, immutable card identity, complete typed source context, selected receipts
and uncertainty. Physical order is preserved and contributes to identity;
repeated materials may have distinct layer IDs. No common-grid resampling or
effective-conductivity claim replaces the individual source curves.

The shared `thermal::conduction::LayerStackConduction` kernel solves one
conserved flux with exact piecewise-affine conductivity integrals. A signed
analytic Kirchhoff inverse traverses original breakpoints; safeguarded
Newton/bisection closes the final layer without a dense interface solve.
Explicit contact jumps determine the layer-face temperatures. Every source
chart must cover the entire endpoint-temperature interval, and trial steps
that cross the cold boundary tighten the bracket without extrapolation.

The response exposes all ordered face temperatures, signed power, secant
conductance, both physical endpoint derivatives, and separate layer-heat and
contact-temperature residuals. Each layer's reconstructed power must meet the
requested heat tolerance. Contact representation allowances remain distinct
from that physical acceptance test. Equal-temperature and single-layer cases
have direct solutions. Successful evaluation allocates nothing, performs no
source queries or hashes, and runs checkpoints per layer and iteration under
the caller's finite iteration budget. Cancelled or failed work publishes no
response and mutates no prepared input.

Feed the secant conductance into the existing thermal owner's candidate-state
coefficient callback; endpoint derivatives serve a separate Newton assembly.
The original timeline continues to own storage, time integration, sources,
energy acceptance and optional sensible entropy. The caller declares fixed
planar geometry, source-free interiors and alignment of scalar conductivity
with heat flow. Unknown source orientation remains unknown. No internal layer
storage, radiation, lateral transport, phase law or inferred contact law is
supplied by this steady constitutive relation.

Four focused tests in `xtask/tests/material_layer_stack.rs` pass against the
production modules: signed inversion and analytic resistance, independent
nonlinear four-layer interfaces and derivatives, refusal/cancellation/budget
behavior, and actual compiled FSRI gypsum/mineral-wool/plywood sources. The
source integration independently constructs 0.1 W conduction and exercises
50 J transfer on the original thermal timeline with energy and sensible-entropy
closure. Existing two-layer parity and the 64-layer limit are included.

## Regional liquid/dense-water/steam caloric paths

`phase::water::regional::RegionalWater` prepares one fixed-pressure original
IF97 path across Regions 1, 3 and 2, including Region 4 liquid/vapor coexistence
when present. Its admitted temperature window lies within 273.16..=1073.15 K
and its pressure within 0<p<=100 MPa; the exact critical singular state remains
unsupported. Only source pieces, saturation and junctions intersecting the
requested window are prepared. Cold compressed feedwater, high-temperature
boiling and superheated or supercritical steam can consequently share one
caloric owner. Low-pressure vapor-only windows and exact saturation endpoint
windows retain their original source admission.

Preparation retains bounded fixed arrays of original source pieces and events,
with no heap allocation or coefficient/reference changes. Each single-phase
regional junction carries the original signed h, u, s and v differences and
temperature/pressure residuals. Callers separately budget enthalpy, entropy,
volume and boundary-temperature matching. A source matching defect is not
latent heat, a new phase, an entropy-production model, or a global equation of
state reconciliation. Original IF97 fits can produce positive gaps or negative
overlaps in h and s independently.

Temperature, h and raw-s queries share the prepared path. Caloric inverses
use original Cp or Cp/T in safeguarded Newton/bisection, bounded density
solves and a total outer-iteration budget shared across candidate pieces.
A generic caloric query refuses an unrepresented gap or multiple original
roots. Fixed-mass heat continuation retains the initial state and chooses the
first admissible root in the heating/cooling direction; it still refuses a
gap rather than inventing storage at a regional boundary. Exact zero heat
retains the state without an inverse. Numerical work and residuals remain
observable. Gibbs evaluations contribute outer iterations; the reported
source-evaluation counter counts the nested dense pressure/stiffness work.

Heat deposition evolves H and independently audits Q=Delta U+p_requested*Delta V.
Stable original-source differences are accumulated separately as sensible,
physical latent, and regional matching contributions in both enthalpy and
entropy. Reports retain phase mass and volume, endpoint pressure-closure
defects, stable differences and independent represented-endpoint residuals.
Supercritical states have no vapor quality or fictitious phase inventory, even
when their source routing reaches Region 2. Invalid budgets, exhausted work,
unrepresentable nonzero heat or phase mass, unsupported/ambiguous states,
chart mismatch and early/late cancellation produce no partial result.

Six focused native tests cover 20/25 MPa heating and cooling across both source
junctions; independent h/s gap and overlap behavior; Cp and Cp/T quadrature;
window-scoped preparation, low-pressure vapor and exact latent endpoints;
subdivision, tiny heat/quality and zero-heat identity; and source, matching,
work, energy and cancellation refusal. This equilibrium pressure owner does
not supply a whole-path Cp stability certificate for `ThermalLaw`, ice,
resolved phase kinetics, spatial two-phase flow, or vessel mechanics.

## Complete high-pressure water cycles through the regional chart

`phase::water::flow::across` accepts `RegionalWater` and
`RegionalWaterState` through additional typed enum variants. The existing
throttle, compressor and expander equations retain original h and raw s,
reuse the regional inverses, propagate typed domain/gap/ambiguity refusal,
and keep the established independent enthalpy and internal-energy/flow-work
audits. Regional preparation and matching budgets remain owned by the chart.
The dispatch introduces no heap allocation, caloric datum shift, duplicate
flow work or alternative device efficiency convention.

`cargo run -p fs-material --example water_regional_rankine` runs two complete
four-device Rankine cycles: a 10 kPa saturated-liquid condenser outlet enters
an explicitly 80% efficient feed pump, a single regional boiler path heats to
823.15 K at 20 MPa or 25 MPa, an explicitly 85% efficient turbine returns to
10 kPa, and the condenser returns exactly to the original saturated-liquid
state. Pump discharge is cold Region 1 water; each boiler traverses Region 3
and reaches Region 2 steam. The 20 MPa boiler includes physical vaporization;
the 25 MPa boiler has zero latent heat and no invented supercritical quality.

The example reports every stage, both original regional junctions, sensible/
latent/matching heat contributions, phase quality, net shaft work, thermal
efficiency, fluid entropy changes and independent energy residuals.
It uses illustrative efficiencies, no pressure-drop model, and no assigned
thermal reservoirs; it does not claim measured machinery performance or total
reservoir entropy production. The earlier hot-inlet steam-train example remains
a partial train and retains its stated limits.

Three integrated regressions exercise both complete cycles and cold-feed
region crossings, recover the existing 4 MPa Gibbs Rankine reference without
changing its h/s datum or device behavior, and check regional outlet-domain and
cancellation refusal. Combined with the six regional caloric tests and the
15 existing dense/caloric/pressure-device tests, the original-source native
water integration has 24 passing tests. Both water examples compile and run.
The full workspace/DSR gate remains unavailable in this checkout because the
required sibling constellation repositories and DSR executable are absent.

## Coherent thermal-fluid transport view

`fluid_transport::ThermalFluidTransport` extends the existing transport view
with its original absolute temperature and optional isobaric volumetric
expansion. Density, Cp, conductivity and viscosity remain the same immutable
checked tuple used by mechanical consumers. Original ideal-gas owners derive
beta = 1/T at fixed composition after their EOS/caloric/acoustic coherence gate.
Source-resolved liquids require one explicitly typed absolute-temperature
coordinate and retain unknown expansion; a missing coefficient never means
zero expansion. An explicit authored-state constructor can declare its own
temperature and optional beta, without acquiring source-data authority.

This is one evaluated state, not a temperature-dependent fluid law or evolving
reservoir. Consumers must keep source/domain and frozen-property assumptions,
and must verify agreement when binding mechanical and thermal snapshots.
Three focused unit definitions cover coherent gas data, authored expansion and
typed liquid temperature. Runtime qualification is reported separately.

## Source-bound Ohmic constitutive graph

`conductor::graph::MaterialOhmicNode` exposes the original scalar conductor as a public constitutive-law node with temperature-change and signed-current inputs, and signed-voltage and nonnegative-Joule-power outputs. It preserves the selected electrical property key and basis, original material bundle and claims, and continuously supported temperature slice, including fixed-context sampled surfaces. Conductivity is queried before conversion to resistance. Its piecewise-consistent tangent differentiates the native resistance and Joule-power equations; interior knots use the right interval and the final endpoint uses the left interval. Construction admits positive, representable native transport across every electrical source knot and refuses unsupported spans, unowned response coordinates, and phase-fraction evolution. The default uses fixed geometry and memoryless scalar DC transport; sourced free expansion can be bound as described below. Thermal storage, time, contacts, skin effects and phase evolution remain with their existing physical owners.

### Source-defined expansion in electrical and thermal feedback

`PreparedFreeExpansion::at_temperature_with_tangent` derives normalized
principal stretches, their exact piecewise temperature derivatives and the
volume derivative from original source segments. It uses the right interval
at interior knots and the left interval at the final endpoint, with one binary
search and no runtime allocation or source query for this geometric evaluation.

`MaterialOhmicNode::with_free_expansion` requires the exact original material
card, reference query, full temperature interval and an explicit current axis.
It derives L=L0*s_parallel and A=A0*s_transverse1*s_transverse2 at the same
temperature as the original electrical property. The graph tangent composes
electrical resistivity/conductivity and all three geometric derivatives.
If an expansion property already belongs to the original conductor bundle,
its selected claim must remain identical. Physical cancellation of material
and geometric temperature derivatives is admitted.

The existing graph/enthalpy solver thus closes strain data, geometry,
resistance/current, Joule heat and fixed-mass temperature on one endpoint.
Five material regressions cover normalized tangents, all current axes, both
electrical source bases, derivative cancellation, fixed geometry and source
pin/refusal behavior. One coupled regression compares three strain tables over
two continuing steps against an independent finite-capacity quadratic and
terminal-work balance. Both focused test targets pass `cargo check`; execution
remains pending. This is unconstrained, quasistatic, stress-free expansion.
Source density claims, reference thermal geometry and fixed-mass caloric
storage retain their original owners. No mechanical work, moving thermal
boundary areas or constrained deformation is inferred from the new binding.

## Spatial heat storage inside sourced layers (2026-10-29)

`thermal::transport::stack::transient::PreparedLayeredThermalBody` prepares
1–64 planar material layers with up to 512 spatial control volumes. Each layer
has an original conductivity query, a separately selected capacity query,
explicit whole-layer mass and a positive uniform subdivision count. Mass-based
capacity must name that whole mass; native volumetric capacity must name the
area-times-thickness volume. Preparation divides the amount among cells and
retains one per-cell storage curve per layer, shared by that layer's immutable
caloric law. Spatial refinement does not duplicate source preparations.

Capacity and conductivity retain the same immutable material card and agree
on shared typed fixed conditions. The caller still declares compatible source
populations and caloric/work constraints; common card identity does not prove
matched specimens. All capacity and conductivity paths share a positive-width
temperature interval containing the explicit common energy datum. Original
receipts, amounts, physical ordering, geometry and subdivision counts contribute
to prepared identity. No density law is inferred from the required mass input.

`BoundLayeredThermalBody` attaches the exact original borrowed storage laws and
checked cell/material/link IDs. Its birth events create independent cell
enthalpies, then contacts. Each internal face uses the original two-layer
conduction kernel across the adjacent half-cell distances. Contact resistance
applies only at actual material interfaces; subdivisions inside a layer have
perfect contact. The response retains each edge's heat flow, interface
temperatures, endpoint derivatives and continuity diagnostics.

Call `step_coefficients` on every candidate of the existing backward-Euler
residual owner. It checks the actual law objects, cell mass and geometry, source
layer identities and link orientation/area before returning sparse conductance
replacements. The original timeline owns energy, sources, evolution, optional
sensible entropy and acceptance. This adapter supplies no second time integrator
or nonlinear solver. Initial conductances are only initial values. Successful
edge kernels allocate nothing; candidate temperature, response and coefficient
containers use O(cells) storage. Checkpoints and explicit edge iteration budgets
bound work; cancellation/refusal changes no physical state.

The model holds geometry, mass, composition and source conditions fixed. Only
the two axial faces are represented; lateral surfaces are insulated. Boundary
heat sources/reservoirs use the existing owner's events. Each cell has one
temperature; spatial refinement and time-step refinement remain independent.
Phase change, moisture transport, changing cure, expansion work and experimental
qualification are not supplied by this sensible finite-volume assembly.

All four tests in `xtask/tests/material_layered_transient.rs` pass against the
production modules. An independently derived insulated cosine mode exercises
2/4/8/16 cells and 1/2/4/8 backward-Euler steps with separate spatial and temporal
convergence. Manufactured nonlinear heterogeneous layers check half-cell and
contact physics, exact caloric integrals, finite-difference tangents, energy
and entropy. A six-cell FSRI gypsum/mineral-wool/plywood wall retains real source
receipts, and a material substitution changes its transient from the same initial
temperature field. Binding/domain/budget/cancellation refusals preserve the
original timeline.

## Pressure-controlled ice/liquid freezing and thawing

`phase::water::freezing::FreezingWater` combines original IAPWS Ice Ih,
the released pressure-dependent Ih melting curve, and original IF97 Region 1
liquid properties. This is an ambient-pressure equilibrium owner: the melting
temperature must be in 273.15..=273.16 K, where the two source charts overlap.
Pressure cannot exceed p_melt(273.15 K), and original liquid source admission
also applies near the triple point and boiling boundary. The immutable
temperature window brackets melting, with ice no colder than 50 K and a
stable Region 1 liquid upper endpoint. A latent-only window is supported.
High-pressure freezing, supercooled liquid, other ice polymorphs, vapor,
nucleation and resolved interfaces are outside this owner.

The prepared pair retains original h, u, s, specific volumes and the actual
melting inverse diagnostics. Both Gibbs mismatch and the pressure-local fusion
entropy defect L/Tm-(s_liquid-s_ice) have explicit independent admission budgets.
Neither property reference is shifted. The defect is a source-consistency
diagnostic, not modeled entropy production.

Temperature states and explicit liquid mass fractions share one source owner.
Enthalpy inversion uses physical lever-rule storage on the fusion plateau and
safeguarded Newton/bisection with analytic Cp on each sensible branch. Exact
source knots and complete phase depletion are retained; temperature alone at
melting is ambiguous. Liquid/ice mass fractions and volume fractions remain
distinct because the phases have different densities.

Fixed-mass heat steps evolve H and independently close Delta U+p*Delta V.
Stable original-source increments and direct endpoint subtraction must both
fit the unchanged caller joule budget. Reports separate sensible and fusion
enthalpy while retaining actual phase masses, physical expansion/contraction,
boundary work, raw source entropy and the pressure-local integral of dh/T.
Latent enthalpy is counted once. The external pressure load permits expansion;
this is not a sealed rigid-vessel pressure model.

Six native tests pass: released normal-pressure melting and original-source
identity/Clapeyron consistency; complete freezing/thawing and expansion;
independent Cp/CpT quadrature; subdivided, latent-only, small and zero heat;
source overlap, matching and work refusal; and chart identity, numerical
representation and early/late cancellation. The tests compile the exact
original production sources against the new ice thermochemistry library.

## Freezing on the shared conservative thermal timeline

`phase::water::freezing::thermal::FreezingWaterThermalLaw` connects the same
ice/liquid enthalpy owner to `fs_thermal::enthalpy::ThermalLaw`. Preparation
obtains outward-rounded whole-interval Cp bounds independently for original
Ice Ih and IF97 Region 1. The maximum reciprocal lower Cp bounds dT/dh over
the entire immutable caloric path; the physical fusion plateau has exactly
zero slope. A latent-only chart prepares no sensible-branch certificate.
The retained bounds are reused across thermal steps and checkpoints.

The adapter uses the existing temperature inversion and stable pressure-local
integral of dh/T. It does not introduce a second heat solver, apparent-Cp
regularization, a sampled stability minimum, or a new energy reference.
Unresolved source bounds and cancellation yield no prepared law. Runtime
domain errors propagate through the existing transactional thermal owner.

`cargo run -p fs-material --example water_freeze_thaw` evolves a uniform
50 g sample at 101325 Pa through cold and warm reservoir contact. The original
sources give Tm=273.152519079769 K and fusion enthalpy 333426.514501674 J/kg.
An explicitly supplied 5 W/K contact conducts to 253.15 K until 300 s and then
293.15 K until 650 s. Temperature and phase inventory emerge from conserved
enthalpy: the sample crosses the fusion plateau, becomes cold ice, thaws, and
returns to warm liquid. The example separately reports liquid mass and volume
fractions, physical free volume, Delta U, pressure work, heat, sample/reservoir
entropy and the temporal entropy defect.

The source Gibbs mismatch is 0.016931408 J/kg and the explicit pressure-local
fusion entropy correction is 6.198518076e-5 J/(kg K); neither is hidden as
physical entropy generation. The native run's maximum independent energy
residual is 6.84e-10 J. This is arithmetic/solver closure, not experimental
prediction accuracy. The example assumes internal equilibrium and supplied
contact conductance. Reference thermal geometry remains fixed while material
volume is reported separately; no resolved front, nucleation, freezing-induced
stress, moving contact area or rigid-vessel behavior is inferred.

Ten integrated native tests pass: the six freezing-core tests plus four
thermal tests for complete branch Cp bounds and latent-only preparation,
reservoir-driven freeze/thaw with timestep-refined entropy defect, isolated
ice/water contact with conserved heat and positive thermal entropy, and
cancellation/domain atomicity. The complete original fs-thermal library,
the freezing material namespace and the executable example build and run
without numerical stubs. Full workspace/DSR verification remains unavailable
because the required sibling constellation checkouts and DSR runner are absent.

## Sparse implicit tangents for layered thermal bodies (2026-10-29)

`thermal::transport::stack::transient::implicit::LayeredThermalStep` evaluates
the existing `BackwardEulerStep` residual and attaches a derivative of that
same material-backed candidate. Its unknowns are cell energy increments in
the original owner's CellId order. For `R = delta_E - dt Q(T(E))`, the retained
operator is `J = I - dt (dQ/dT) diag(1/C)`, using original per-cell capacity and
physical endpoint derivatives of integrated conductivity/contact flux.
Secant conductance is not substituted for a Newton derivative.

The operator admits exactly its bound body cells and links. Authored constant
cell powers, convection on the original exposed areas, and reservoir exchange
remain with the original heat-rate owner. Its boundary slopes include reservoir
coverage and fixed convection/conductance. Pending event admission and target
left-limit semantics remain those of `BackwardEulerStep`; larger coupled graphs
continue to use the existing coefficient adapter and their complete Jacobian.

Only three diagonals are retained. Tangent assembly and matrix actions require
O(cells) storage/work after material flux evaluation. Forward, transpose and
Jacobi inverse-diagonal actions allocate nothing and validate arithmetic and
cancellation before publishing their bounded output. Unequal capacities and
nonlinear conductivity generally make the operator nonsymmetric; the transpose
uses the actual opposite diagonals. Every tangent is bound to its exact trial,
with original temperatures, capacities and link responses available for clients.

The existing thermal trial still owns residuals, optional entropy and independent
local/total joule acceptance. This addition supplies no new nonlinear solver,
time integrator, trajectory adjoint, condition-number certificate or automatic
time-step accuracy claim. The existing `fs-la` prepared LU is used only by the
focused integration test's bounded small-system Newton iteration.

All four tests in `xtask/tests/material_layered_implicit.rs` pass. Independent
constant coefficients include sources, convection and reservoir coverage;
nonlinear finite differences verify the full energy-coordinate matrix and the
transpose inner-product identity. Invalid/cancelled matrix actions preserve
outputs. An actual six-cell FSRI wall converges at a 3,000-second implicit step
and passes original energy/entropy acceptance; material substitution changes
the endpoint. This establishes algebraic convergence at that step, separately
from time-discretization accuracy. Scope, source, event, identity, overflow and
cancellation refusals preserve the initial timeline.


## Source-backed spatial water freezing

`phase::water::freezing::transport` binds the unchanged original freezing chart
to independent liquid IAPWS SR6-08(2011) and empirical freshwater-ice Yen (1981)
conductivity. This composition deliberately admits standard atmospheric
pressure (101325 Pa) only, with the complete caloric temperature window inside
the transport source domains. Absent phases are not extrapolated. Melting
transport values are prepared once; sensible-phase queries retain original
source admission. Partially frozen cells require the caller to select series
or parallel layers. These directional closures use physical phase VOLUME
fractions from distinct source densities, never mass fractions or an inferred
isotropic mixture. The ice correlation supplies no pressure response or
experimental uncertainty envelope. Transport does not alter h, u, s or fusion
enthalpy.

`freezing::spatial::FreezingSlab` assembles a stationary one-dimensional plane
reference mesh on the existing `fs-thermal` timeline and geometric conduction
observer. Initial source density fixes conserved reference mass. Both faces
may have explicit temperature reservoirs with finite areal resistance;
otherwise they are insulated, as are all lateral faces. Scheduled reservoir
temperature changes use the original event owner. State-dependent internal and
half-cell boundary conductances are refreshed on every actual CFL-limited
substep. The latent plateau remains the exact original h-to-phase relation;
no apparent-Cp regularization or separate latent source is introduced.

Accepted advances independently check original-source Delta H = Q and
Delta U + p*Delta V = Q against the supplied joule tolerance. The original
thermal and reservoir entropy accounts accompany the separate temporal
entropy defect. Cell profiles expose actual enthalpy, temperature, phase mass
and volume fractions; inventory reports distinct phase masses and free volumes.
The reported equivalent frozen REFERENCE length is thickness times frozen mass
fraction, not a mechanically resolved interface or physical ice thickness.
Reference heat-transfer geometry remains frozen while free expansion and
pressure work are reported separately. This conduction model does not solve
convection, stress, rigid-vessel pressure, nucleation, supercooling, front
kinetics or a deforming boundary. It has no unconditional spatial/time accuracy
claim. Source-backed mesh/time refinement and the shared owner's independent
Stefan reference are numerical implementation evidence, not experimental
validation of every freezing specimen.

The runnable `water_freezing_front` example computes an advancing phase field
inside a 10 mm water layer and reports the resolved profile, phase inventories,
heat, internal energy, pressure work and entropy. Focused tests also cover
adiabatic equilibrium, a bath reversal, source/chart binding, energy balance,
work limits and cancellation without altering accepted state.

## Native thermoelectric constitutive graph

`thermoelectric::graph::MaterialSeebeckNode` binds the original absolute thermopower source to the public constitutive graph through three temperature-change inputs and signed current. Every value is supplied by `ResolvedSeebeckSpan`, retaining both terminal Peltier transfers, the original bulk Thomson exchange, and electrical supply. The selected claim, query quantities, fixed external context, and continuously supported temperature slice remain pinned. Its consistent tangent integrates piecewise source derivatives, including reversed edges and the isothermal limit, without differencing nearly equal thermopowers. The node declares reversible transfer with no local storage and explicitly zero dissipation; graph execution and equation-compiler execution retain that declaration, and the compiler codec preserves it. Resistance, Joule heat, Fourier conduction, caloric storage, physical time, and acceptance of complete coupled balances remain with their existing owners. No phase evolution, magnetic transport, relative-to-reference thermopower substitution, or measured authority beyond the original source is inferred.

## Radial material conduction and finite thermal feedback

`thermal::conduction::LayerGeometry` admits planar, concentric cylindrical and
spherical source-free conduction. It lowers the exact metric
`A_contact * integral(ds/A(s))` to equivalent planar distances, retaining
physical radii and contact area. The original `TwoLayerConduction` integrates
each source k(T), solves contact temperatures and supplies endpoint derivatives.
The contact jump uses the actual area; thin-radius differences use stable
logarithm/difference forms. No source value is inferred from a material name.

`PreparedLayerConduction::new_geometric` binds those metrics to the original
active contact, source-qualified pure-phase caloric objects and per-side
coverage. `LayerConductionObserver` supplies their integrated conductances to
the existing explicit enthalpy/heat-rate/CFL owner; `resolve` supplies the same
implicit endpoint residual. One signed internal transfer changes both finite
stores, with no external heat or duplicate storage. Source/domain/event and
cancellation refusal discard the full proposed advance. Source selection and
geometric metrics are prepared once; flux evaluation performs no database query.

The six `radial_conduction` regressions cover independent cylindrical/spherical
resistance, nonlinear conductivity integrals/tangents, thin shells and the
exterior limit, data-only rho/Cp/k changes between finite solid/liquid cells,
native time refinement/energy/entropy, implicit endpoint parity and refusal
with replay. Their integration target and the existing `layer_conduction` and
`layer_conduction_network` targets compile in the scoped Cargo check; runtime
execution is not established by that check.

This is stationary radial conduction at fixed mass/reference geometry with
quiescent liquid and explicit symmetry assumptions. The legs have no independent
thermal storage; transient spatial resolution and coefficient-freezing error
require refinement. No convection, fluid pressure/motion, confinement correction,
phase-transition propagation or experimental accuracy is claimed.

### Shared electrical geometry for native thermal networks

`MaterialOhmicNode::resolve_from_state` lowers an already resolved complete
electrical/thermal bundle through the same prepared L/A and optional free
expansion used by graph evaluation. It retains the card, full source
requirements/keys/claims, quantity schema and fixed conditions; only the
declared temperature may change. No second property query occurs. Original
source-domain admission belongs to the canonical resolved bundle, and geometry
remains bounded by the prepared interval. A focused regression checks parity
with ordinary node resolution, bundle preservation and foreign/out-of-range
refusal. The `conductor_graph` integration target compiles; this check does not
establish runtime execution of that regression.


## Implicit spatial freezing with original source tangents

`phase::water::freezing::transport::FreezingConductivityState` now retains
analytic local dT/dh and dk/dh with its original material/transport response.
Sensible branches use the reciprocal original Ice Ih or IF97 Cp and the
original conductivity temperature derivative. Fusion retains exactly zero
dT/dh. Its conductivity derivative follows the physical liquid VOLUME
fraction, dphi/dh = vi*vl/(v_mix^2*L), through the caller's series or parallel
closure. Both exact fusion endpoints select the plateau-side generalized
derivative. No apparent heat capacity, smoothed transition, finite-difference
runtime tangent, extra latent source or shifted energy reference is used.

`FreezingWaterImplicitLaw` borrows the exact original caloric and transport
objects, verifies chart equality, and provides one fused source query to the
shared thermal chain adapter. Bounds and all source-domain refusals remain
those of the original immutable material chart. The extra scalar tangents are
computed during the existing query and add no allocation or inverse query to
its getters.

`freezing::spatial::FreezingSlab::advance_implicit_to` uses the original
backward-Euler residual and the reusable `ImplicitConductionWorkspace`.
Its endpoint conductances include the actual phase closure and face half-cell
resistance. The same slab can continue with explicit or implicit transport.
Both paths share independent original-source Delta H=Q and Delta U+p*Delta V=Q
checks and retain original thermal/reservoir entropy. The nested implicit
report repeats the physical energy account for comparison; callers must not
add it a second time. Cancellation or failed nonlinear/work/energy admission
leaves the original slab intact, and numerical scratch remains reusable.

The existing `water_freezing_front` example accepts `--implicit`,
`--cells N` and `--dt SECONDS`. For example,
`cargo run -p fs-material --example water_freezing_front -- --implicit --cells 96 --dt 0.5`
computes the 120-second 10 mm slab with 244 accepted intervals and a maximum
independent energy residual about 1.1e-10 J in the native run. This is algebraic
closure, not physical or time-discretization accuracy. At 120 seconds the
frozen-reference equivalent length is 5.456292577 mm; halving the maximum
step gives 5.449258361 mm. These are phase-inventory lengths on the fixed
reference mesh, not mechanically resolved ice thicknesses.

All 24 focused freezing-namespace tests pass against the actual updated
thermochemistry and thermal libraries. New tests compare analytic source and
phase-volume tangents with independent differences, verify exact borrowed
bindings, resolve an actual source-backed front beyond the initial explicit
CFL limit, compare three implicit timesteps with a fine explicit reference,
reverse the bath at an exact authored time, reproduce chunked continuation
bit-for-bit, and reuse scratch after cancellation/work refusal. Refinement
errors remain separate from the much smaller joule acceptance tolerance.

All previous physical limits remain: atmospheric source/transport overlap,
declared phase morphology, fixed reference geometry and external pressure.
This addition supplies no convection, stress, nucleation, supercooling,
rigid-vessel coupling or moving surface. Implicit interval retries do not
provide automatic temporal error control or unconditional nonlinear
convergence. Full workspace/DSR qualification remains unavailable in the
current checkout; the native validation compiles the original module sources
against actual production dependency artifacts, without numerical stand-ins.

## Source-defined expanding thermal surfaces

`thermal::transport::expansion::surface::PreparedExpansionSurface` binds explicit
reference areas grouped by principal normal to one original
`PreparedFreeExpansion`. At each admitted temperature it derives
`A_i = A_i,ref * s_j * s_k` and the exact temperature tangent. Another consumer
can reuse an already evaluated original tangent only when the full prepared
path identity matches. Area evaluation performs no runtime allocation, hashing
or material query. Original card/claim references, fixed query conditions,
reference normalization and the declared principal/isotropic basis remain
authoritative. Zero selected faces remain zero; electrical L/A does not infer
a thermal surface or boundary coverage.

Three `expansion_surface` regressions cover normalized anisotropy and one-sided
source knots, exact tangents and selected/zero patches, legitimate derivative
cancellation, foreign paths and domain/area refusal. This owner admits aligned
faces under diagonal unconstrained free expansion. Arbitrary curved normals,
shear, constrained stress and moving interbody contact coverage require their
own geometry owners. Thermal boundary laws and disjoint coverage admission
remain with the consumer; surface geometry alone adds no heat or storage.

## Sealed dense-water heat exchange on the shared thermal timeline

`phase::water::dense::thermal::DenseIsochoricThermalLaw` retains the original
prepared IF97 Region 3 isochore. Its declared storage basis is internal energy,
with original T(u), whole-interval positive-Cv slope bound, local analytic
dT/du = 1/Cv, and original-source entropy increments. Birth admission requires
m = rho V within a 64-epsilon arithmetic allowance and refuses unsupported u.
Mass advection and the existing sublimation owners refuse this storage basis.

`DenseWaterThermalNetwork` couples sealed rigid volumes through the original
thermal contacts, reservoir boundaries, events, timestep owner and immutable
candidate publication. All volumes are born at time zero. Stable source
differences and represented endpoints independently check Delta U = Q and
Delta H - V Delta p = Q. Moving-boundary work is exactly zero; V Delta p is
an enthalpy contribution, not wall work. Thermal-coordinate entropy and original
source entropy are reported separately, exposing finite inverse representation
effects without discarding tiny representable heat.

Four new integration tests, together with four existing dense-source tests, cover
source Cv/tangents and tiny-heat entropy, mass/volume admission and actual
advection refusal, finite-body equilibration, exact bath events, continuation,
cancellation, budgets and source-domain failure. The bath transient converges
against independent quadrature of t = (m/G) integral Cv(T)/(T_bath - T) dT:
at 20.695034947 s, temperature errors for maximum steps 0.4/0.2/0.1 s are
0.241386644/0.120461412/0.060267601 K.

`water_dense_vessel --bath` computes finite-conductance heating and cooling,
pressure and reservoir entropy; its original prescribed-heat mode remains
available. The example's 0.5 kg, 1-litre vessel with authored 25 W/K conductance
reaches 732.426872 K and 68.727035 MPa at 100 s, with Q = 120810.286679 J,
Delta H = 163953.619467 J and V Delta p = 43143.332788 J. Across 800 steps its
maximum history energy residual is 2.6663e-7 J.

This is a spatially uniform network on one fixed-density chart. Conductances are
authored; dense-water conductivity, container strength, wall storage/compliance,
internal convection, multiphase behavior and regional stitching are not supplied.
Native validation compiles the actual dense namespace and complete thermochemistry
and thermal owners; it does not constitute a full workspace gate.

## Adaptive spatial freezing with actual phase-inventory accuracy

`FreezingWaterImplicitLaw::time_observables` obtains temperature and liquid
mass fraction from one original caloric inverse. `FreezingSlab::advance_adaptive_to`
uses the shared implicit accuracy controller and retains the existing independent
source Delta H, Delta U + p Delta V and entropy accounts. No apparent heat
capacity, duplicate fusion source, extrapolated accepted state or synthetic
temperature-to-phase relation is added.

`FreezingSlabAdaptiveAudit` separates the physical audit from the controller's
local h/T/phase estimates and all attempted nonlinear work. Source and controller
failures remain transactional, including cancellation after the final original-
source audit. A returned `suggested_dt_s` explicitly seeds continuation.
Scheduled face-temperature changes retain their exact event times and actual
reservoir entropy temperatures.

Two additional source tests bring the actual freezing namespace to 26 passing
tests. A 24-cell, 12-second comparison against independent explicit transport at
a 0.002-second cap reduces mass-weighted absolute enthalpy-field error from
1.637753637 to 0.226923967 J when local h/T/liquid-fraction allowances tighten
by 100. Maximum temperature discrepancy falls from 0.046597616 to 0.006387667 K,
maximum liquid-fraction discrepancy from 0.011288234 to 0.001564972, and bath
heat discrepancy from 1.498449183 to 0.207880577 J. Accepted comparison pairs
increase from 141 to 1382: tighter accuracy costs additional original solves.
The second case tests an exact cold-to-warm bath event, continuation at that
event with the returned proposal, two-bath entropy and final-audit cancellation.

`water_freezing_front --adaptive` exposes dimensioned enthalpy and temperature
allowances, a liquid-mass-fraction allowance, a relative change allowance, and
minimum/maximum comparison spans. The 24-cell, 120-second default adaptive run
accepts 762 pairs (1524 physical half steps) from 817 comparisons. Its maximum
accepted error ratio is 0.925807011 and maximum independent calorimetric residual
is 1.1102e-9 J. Computed frozen reference-equivalent depth is 5.225752655 mm.

These are local temporal estimates and computed reference-geometry results,
not global or spatial accuracy certificates. Spatial refinement, source ranges,
the explicit mixed-cell conductivity closure, and the fixed-pressure/no-flow
limits of the original slab remain independent requirements.
