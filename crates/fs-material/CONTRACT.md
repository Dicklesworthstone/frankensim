# fs-material — CONTRACT

The constitutive-law and authority-separated identifiability kernel: materials
as mathematical objects — calibration domains, CONSISTENT tangent operators,
thermodynamic guardrails, hysteresis, uncertainty — plus source-closed,
multi-case questions whose execution and assessment identities stay distinct.

Ambition tags: elastic/hyperelastic/J2/RC-fiber laws + calibration [S];
multi-case identifiability authority schema [F]. Ogden and tension-only
HGO/GOH are [S] constitutive laws; specimen/tissue calibration is separate.

## Purpose and layer

Layer **L3** (FLUX support). Runtime deps: `std`, fs-ad (dual-number
energy derivatives), fs-blake3 (typed content identities), fs-evidence
(model cards, Evidence, and V&V artifacts), fs-matdb (immutable material and
constitutive-model cards), fs-qty, fs-math.
Consumers: fs-solid elasticity (tfz.13), fiber-section beams, lattice
homogenization, the P2 milestone.

## Public types and semantics


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
- `J2Plasticity`: radial-return mapping with linear isotropic hardening
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
  New variants enforce all principal stretches in `[1/limit,limit]` with
  eight-EPSILON eigensolver boundary roundoff; limit*(1+1e-12) refuses.
  Legacy NH/MR retain their existing card-only stretch limit (an explicit
  remaining admission gap). Both new variants offer quadratic or log-squared
  bulk energies and separate distortional/bulk accessors. Full-invariant
  HGO's fibre accessor intentionally retains volume dependence.
- `hyper::cards` instantiates strict, dimensioned `ConstitutiveModelCard`s
  for `ogden-decoupled` v1 and `hgo-dispersed` v1 and registers the real
  stress/tangent/energy factories in `graph::LawRegistry`. Graph inputs are
  H=F-I, outputs P, so the zero-input admission probe is the reference state.
  Moduli require Pa, exponents/directions/counts/selectors are dimensionless;
  missing/extra parameters, fractional counts, versions and undeclared
  invariant selectors refuse. FSMODPK v1 remains unchanged.
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
  tests transport both laws into the real registry. A generic profile does
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
- **J2 has no Bauschinger effect** (isotropic hardening only; kinematic/
  mixed hardening is a follow-up law, the card says so).
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

The selected hardening basis is either the plastic-strain modulus H or a
monotonic uniaxial total-strain tangent Et. For the latter, H = Et E/(E-Et) is
derived after resolving both coefficients at the same state; 0 <= Et < E is
required. Perfect plasticity needs an explicit zero source. The original Et
claim and uncertainty remain intact, and no derived measured-H claim is created.
Source, quantity, domain, parameter and representability errors refuse.
Coefficients remain frozen within one load path; temperature re-resolution
starts a different law, not an energy-accounted transfer of plastic history.
The declared strain bound is limited to 0.05 and remains a caller obligation.
Creep, rate effects, softening, phase changes and finite strain are outside this law.

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

Callers must retain the same elastic coefficients and H across the increment.
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
energetic hardening H may be globally affine in temperature. Poisson ratio,
density and reference thermal geometry remain fixed. All property roles share
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
energy. `sigma_y(T1) delta_alpha` is the physical plastic heat supplied to the
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
and extrapolation refuse. Runtime H uses the canonical pinned query. The same
piecewise-affine path supplies H/B, dH/dB, exact segment-integrated energy
`psi(B)=integral_0^B H(b) db` and coenergy `B H-psi`. The tangent at zero uses
the first segment; an interior knot uses its right derivative. These are
floating-point model evaluations, not uncertainty propagation or enclosures.

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
over an explicit single-region temperature bracket using the original NASA-9
coefficient interval bound and each species' `cv = cp - R/M`. The bound is a
timestep aid; the existing source model still owns actual `u(T)` and inversion.
`IsochoricGasParcel::model` exposes that same immutable composition/source model
for coupled film queries. No constant heat capacity, enthalpy substitution or
second gas caloric law is introduced. A focused test compares the bound with
actual mixture heat capacities and refuses a bracket outside source support.

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
`delta_alpha = max(f_trial,0)/(3 mu + H + eta/dt)`. Positive endpoint overstress
is a legitimate committed state; the consistency gate checks
`q - sigma_y - H alpha - eta delta_alpha/dt`. Every new interval starts from
the accepted internal variables. Native viscosity and duration must be finite,
positive and representable together with the return denominator.

The energetic-hardening interpretation remains explicit:
`psi=e:C:e/2 + H alpha^2/2`. Physical dissipation is
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
overstress with energetic isotropic hardening, not a nonlinear rate exponent,
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
