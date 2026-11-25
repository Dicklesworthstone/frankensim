# CONTRACT: fs-multiphase

## Purpose and layer

`fs-multiphase` is an **L3 (`FLUX`)** crate providing deterministic, constitutive
two-phase wall heat- and mass-transfer closures:

- Full five-regime Nukiyama pool boiling (`boiling::BoilingCurveV1`): natural
  convection, Rohsenow nucleate boiling multiplied by exact implicit
  Schlunder/Thome-Shakir binary mixture suppression (driven by
  `fs-thermochem::solution::nrtl` `Delta T_bp = T_dew - T_bub`), Lienhard-Dhir /
  Zuber critical heat flux (`C_CHF = 0.149` flat plate default), safeguarded
  `C^1`-smooth monotone cubic Hermite transition boiling, and Berenson/Bromley +
  Henry film boiling / Leidenfrost minimum (distinguishing saturated pool from
  subcooled/impinging dynamic Leidenfrost).
- Nusselt laminar, wavy-laminar, and turbulent filmwise condensation on vertical
  plates, inclined plates, and horizontal cylinders (`condensation::CondensationModelV1`)
  coupled with Colburn-Hougen non-condensable gas boundary-layer resistance.
- Unified single-source phase-change surface interface
  (`BoilingCondensationInterfaceV1`) owning wall phase-change strictly when
  `T_w > T_sat(P)` (boiling) or `T_w < T_dew(P, w_nc)` (condensation) with zero
  double-counting against sub-boiling convection (`fs-convection`).
- Stirred multiphase liquid/vapor/oil/emulsion vessel component
  (`vessel::MultiphaseVesselComponent`) implementing
  `fs_ports::component::Component` for stockpots, saucepans, skillets, pressure
  cookers, and commercial deep-fryer vats with dryout regularization, NRTL wine
  deglazing & internal fond dissolution, Carreau/Herschel-Bulkley sauce rheology,
  Bikerman boil-over foam, butter water-flash to *beurre noisette* browning,
  sub-tube fryer cold-zone stratification, thermal safety bounds, and
  atom-balanced triglyceride hydrolysis (`triglyceride + 3 H2O -> 3 ffa + glycerol`).
- Closed-loop Péneloux-translated Peng--Robinson vapor-compression refrigeration
  cycle components (`cycle::RefrigerationCompressorComponent`,
  `cycle::CondenserEvaporatorHXComponent`, `cycle::ExpansionThrottleComponent`,
  and `cycle::IceHarvesterComponent`) across `R-290`, `R-134a`, `R-404A`
  (3-component blend), `R-448A` (5-component zeotropic blend), and subcritical /
  transcritical `R-744` CO2 gas-coolers on `ENTHALPY_REF_IIR_273K` (`iir-273.15K`).

Layer dependencies stay strictly within L0/L1/L3/UTIL (`fs-math`, `fs-exec`,
`fs-qty`, `fs-evidence`, `fs-matdb`, `fs-thermochem`, `fs-ports`), with zero
external runtime dependencies and zero edits to `fs-couple` or `fs-material`.

## Public types and semantics

- `boiling`:
  - `BoilingFluidKindV1` (`PureWater`, `EthanolWaterNrtl`, `CustomFluid`,
    `TriglycerideOil` which explicitly refuses Nukiyama boiling).
  - `BoilingFluidPropertiesV1` (saturated liquid/vapor densities, transport,
    surface tension, latent heat, boiling range `Delta T_bp`, and
    `janaf-298.15K-formation` specific enthalpies).
  - `BoilingSurfaceCardV1` (Rohsenow `C_sf`, Prandtl exponent `n`,
    Lienhard-Dhir `C_CHF`, contact angle, wall thermal properties for Henry
    effusivity shift, emissivity, characteristic length; loadable from
    `fs-matdb::ConstitutiveModelCard` with law id `"rohsenow-boiling-surface"` v1).
  - `LeidenfrostModeV1` (`SaturatedPool` vs `SubcooledImpinging`).
  - `BoilingCurveV1`, `BoilingEvaluationV1`, `BoilingRegimeV1`,
    `BoilingRegimeBoundariesV1`, `BoilingStepAdviceV1`, `BoilingErrorV1`.
- `condensation`:
  - `CondensationGeometryV1` (`VerticalPlate`, `InclinedPlate`, `HorizontalCylinder`).
  - `CondensationFluidPropertiesV1`, `NonCondensableGasV1`, `CondensationModelV1`,
    `CondensationEvaluationV1`, `CondensationRegimeV1`, `CondensationErrorV1`.
- Root (`lib.rs`):
  - `BoilingCondensationInterfaceV1`, `SurfacePhaseChangeEvaluationV1`,
    `SurfacePhaseChangeModeV1`, plus re-exports of the primary boiling,
    condensation, and vessel types and `JANAF_FORMATION_DATUM_ID_V1`.
- `vessel`:
  - `MultiphaseVesselComponent` implementing `fs_ports::component::Component`
    with 7 ports (`wall_heat`, `immersed_food_heat` field width `M >= 1`,
    `liquid_stream` on `kitchen_liquid_v1`, `vapor_stream` on `kitchen_gas_v1`,
    `stir_rpm`, `p_vessel_pa`, and 12-wide `probes`) and literal
    `StateIdentityCharterV2` (`owner = "fs-multiphase::MultiphaseVessel"`,
    `state_family = "stirred-vle-v1"`, `codec_version = 1`).
  - `MultiphaseVesselConfig`, `MultiphaseVesselState`, `VesselInitialStateV1`,
    `OilCardV1` (loadable from `kitchen-cooking-oil` v1 cards),
    `OpenVentSpecV1`, `VesselSafetyFlagsV1`, `VesselErrorV1`.
- `cycle`:
  - `RefrigerantModelV1` built directly on `fs-thermochem::eos::peng_robinson`
    with Péneloux molar volume shifts (`new_with_volume_shifts`),
    `pure_saturation_at_temperature`, and `two_phase_flash_pt` across all five
    `KM0.3` commercial refrigerants (`R-290`, `R-134a`, `R-404A` 3-component,
    `R-448A` 5-component zeotropic blend, and `R-744` CO2 subcritical +
    transcritical gas-cooler) anchored on `ENTHALPY_REF_IIR_273K`
    (`iir-273.15K`, `h_l(0 C) = 200 kJ/kg`, `s_l(0 C) = 1 kJ/(kg K)`), plus
    `from_card` loader for `peng-robinson-refrigerant-v1` cards.
  - Four reusable `fs_ports::component::Component` physics owners, each with a
    self-contained literal `StateIdentityCharterV2` and deterministic bit-exact
    word codec:
    - `RefrigerationCompressorComponent` (`owner = "fs-multiphase::RefrigerationCompressor"`,
      `state_family = "pr-compressor-v1"`, `ElectricalVoltageCurrent` + `refrigerant_v1`
      + `shell_heat` + `fouling_factor`).
    - `CondenserEvaporatorHXComponent` (`owner = "fs-multiphase::CondenserEvaporatorHX"`,
      `state_family = "pr-hx-frost-v1"`, `CondenserGasCooler` and `FrostingEvaporator`
      modes with dynamic frost thickness ODE and `drain_stream` meltwater removal).
    - `ExpansionThrottleComponent` (`owner = "fs-multiphase::ExpansionThrottle"`,
      `state_family = "txv-isenthalpic-v1"`, exact isenthalpic `h_out == h_in`,
      Péneloux liquid density orifice flow, superheat-feedback modulation).
    - `IceHarvesterComponent` (`owner = "fs-multiphase::IceHarvester"`,
      `state_family = "pr-ice-harvester-v1"`, air-cooled commercial ice machine
      with potable `water_feed_stream`, Ice Ih bin, and kitchen `room_heat` rejection).

## Invariants

1. **No hardcoded sea-level pressure or boiling point:** Every boiling,
   condensation, and vessel query takes an explicit caller-supplied pressure
   (`P_atm(z)` or regulated vessel pressure `p_vessel_pa`).
2. **`C^1` smoothness across all five Nukiyama regimes:** Heat flux `q''(T_w)`
   and its analytical derivative `dq''/dT_w` are continuous across every regime
   knot (`T_sat`, ONB blend, CHF cap, transition spline, Leidenfrost valley,
   film boiling), with exact analytical `dq''/dT_w` matching central finite
   differences to `< 1e-8`.
3. **Strict non-boiling lipid refusal (Round-2 Physics #3):** Routing pure
   triglyceride cooking oils into Nukiyama boiling unconditionally returns
   `Err(BoilingErrorV1::NonBoilingLipidPhase { .. })`.
4. **Unified `janaf-298.15K-formation` enthalpy datum (Round-2 Physics #1):**
   All water/steam liquid and vapor specific enthalpies and all
   `MultiphaseVesselComponent` species enthalpies use single-accounting total
   enthalpy on `janaf-298.15K-formation`, guaranteeing First-Law energy closure
   `< 1e-11` across boiling, wine deglazing, butter flashing, oil hydrolysis,
   cold-zone exchange, and boil-over.
5. **Single-source ownership handoff:** Boiling returns zero phase-change flux
   for `T_w <= T_sat`; condensation returns zero flux for `T_w >= T_dew`.
6. **Non-condensable monotonicity:** Increasing bulk non-condensable air mass
   fraction strictly decreases the Colburn-Hougen interface temperature `T_i`
   and effective condensation heat-transfer coefficient.
7. **Dryout regularization & non-negativity:** `f_wet = tanh(m_liq / (rho_l * A_base * delta_min))`
   clamps evaporation strictly to available volatile mass and emits
   `Proposal::Shorten` on dryout approach; stored mass and carried inventories
   never become negative or `NaN`.
8. **Exact internal fond ownership & hydrolysis stoichiometry:** Adhered pan
   fond (`m_fond_adhered_kg`) is stored inside `MultiphaseVesselState` and
   included in `stored_mass_kg()` and `stored_carried_quantity("solute_fond")`
   before and after dissolution (`< 1e-11` closure). Triglyceride hydrolysis
   (`C57H104O6 + 3 H2O -> 3 C18H34O2 + C3H8O3`) closes `ElementalMatrix`
   (`A * (m_dot / M) == 0`) and `sum(m_dot) == 0` under `audit_carried_reactions`.
9. **Closed-loop vapor-compression First & Second Laws (`G0`):** Every cycle
   evaluation across `R-290`, `R-134a`, `R-404A`, `R-448A`, and subcritical /
   transcritical `R-744` satisfies `|Q_cond - (Q_evap + W_comp - Q_shell)| < 1e-11`
   relative, `COP < COP_Carnot` strictly, `S_gen >= 0`, exact isenthalpic TXV
   throttling (`h_out == h_in`), and explicit near-critical uncertainty
   disclosure (`18%` band for `R-744` above `0.95 * T_c`).
10. **Conservative frosting & meltwater drainage:** Evaporator frost growth
    monotonically degrades `UA_evap(delta_frost)` and defrost drains liquid
    meltwater strictly OUT through `drain_stream` (`kitchen_liquid_v1`) with
    zero re-evaporation into box air and `< 1e-11` mass and enthalpy closure.

## Error model

- `BoilingErrorV1`: typed refusals for `NonBoilingLipidPhase`, `InvalidPressure`,
  `InvalidTemperature`, `InvalidParameter`, `InvalidComposition`,
  `InvalidRegimeOrdering`, `Thermochem`, and `Unrepresentable`.
- `CondensationErrorV1`: typed refusals for `InvalidPressure`,
  `InvalidTemperature`, `InvalidGeometry`, `InvalidNonCondensableFraction`,
  `InvalidParameter`, `Thermochem`, `ConvergenceFailed`, and `Unrepresentable`.
- `VesselErrorV1`: typed refusals for `InvalidConfig`, `InvalidScalar`,
  `InvalidPressure`, `InvalidCard`, `InconsistentStreamMassFlow`,
  `DepletedInventory`, `EmptyVesselThermalSingularity`, `Thermochem`,
  `PortContract`, `SolverConvergence`, and `Cancelled`.
- `CycleErrorV1`: typed refusals for `InvalidConfig`, `InvalidScalar`,
  `InvalidPressure`, `InvalidTemperature`, `InvalidCard`,
  `ThermodynamicInversion`, `InconsistentStreamMassFlow`, `DepletedInventory`,
  `PengRobinson`, `Saturation`, `Flash`, `ConvergenceFailed`, `PortContract`,
  and `Cancelled`.
- No panics or silent `NaN`/inf propagation on any invalid input.

## Determinism class

`cross-isa` (deterministic binary64): all transcendental evaluations (`exp`,
`expm1`, `ln`, `pow`, `sin`, `sqrt`, `tanh`) route strictly through
`fs_math::det` with zero raw platform `libm` calls, zero `f64::powi`, and zero
hash-ordered iteration. Repeated evaluations and `Component::propose` steps are
bit-identical.

## Cancellation behavior

All boiling and condensation evaluations are stateless, allocation-free (or
bounded fixed-iteration scalar solves of at most 64-128 iterations) leaf
constitutive calls that complete in microseconds. `MultiphaseVesselComponent::propose`
polls `cx.checkpoint()?` at step entry and executes bounded `O(1)` analytical
enthalpy inversions plus one bounded scalar VLE root solve with zero locks or
background threads.

## Unsafe boundary

`#![forbid(unsafe_code)]` at crate root (`src/lib.rs`). Zero `unsafe` blocks
anywhere in `fs-multiphase`.

## Feature flags

None in v1 (all boiling, condensation, and vessel laws compile in the default
crate configuration).

## Conformance tests

- Unit tests in `src/boiling.rs`, `src/condensation.rs`, and `src/vessel.rs`:
  pure-water CHF and Berenson/Henry limits, Schluender/Thome zero-range identity,
  analytical vs central-FD `dq''/dT_w`, Nusselt laminar/wavy/turbulent branches,
  Colburn-Hougen air-fraction monotonicity, exact hydrolysis molar-mass closure,
  quadratic enthalpy inversion round-trip, durable codec round-trip, 3-5x
  stirring HTC scaling, and cancellation/pressure refusals.
- Integration suite in `tests/nukiyama_and_nusselt.rs`:
  - Incropera 1 atm water Nukiyama sweep (`q_CHF ~ 1.26 MW/m^2 +- 5%`,
    `T_min ~ 195-220 C` saturated pool vs `250-300+ C` dynamic skillet, SS304
    vs Cu `C_sf` shifts) and Denver 83.5 kPa boiling-point drop (`~94.6 C`).
  - Wine-deglazing (12-14% ABV ethanol-water via `KM0.1` NRTL) demonstrating
    30-60% Schluender/Thome nucleate boiling suppression vs ideal pure Rohsenow.
  - Analytical `dq''/dT_w` central-FD verification to `< 1e-8` across all five
    Nukiyama regimes and `C^1` Hermite transitions.
  - Coupled implicit-midpoint pan-quench transient sweeping Leidenfrost ->
    Transition -> CHF -> Nucleate with zero 2-cycling and `< 1e-11` energy
    conservation.
  - Exact laminar vertical-plate Nusselt oracle to `< 1e-10`, inclined plate &
    horizontal cylinder checks, and monotonic Colburn-Hougen air degradation.
  - Planted negatives: triglyceride oil refusal, non-positive pressures /
    temperatures / `C_sf`, out-of-range air mass fractions, and `NaN`/inf inputs.
- Integration suite in `tests/vessel_physics.rs`:
  - `fs_ports::component::conformance::check_component` + durable
    `StateIdentityCharterV2` checkpoint round-trip.
  - Emergent USDA wine-simmer ethanol retention driven by `KM0.1` NRTL + `KM2.1`
    boiling flux (zero hardcoded decay tables) + internal fond dissolution with
    `< 1e-11` mass and energy conservation.
  - Bikerman pasta-pot boil-over (`boil_over_flag == 1.0`) and collapse on
    power drop.
  - Unsalted & salted butter water-flash plateau near `100 C` -> `> 120 C`
    milk-solid *beurre noisette* browning vs dry clarified butter/ghee and canola.
  - 40-hour commercial deep fryer: sub-tube cold-zone thermal stratification,
    FFA/TPM/smoke-point/viscosity aging curves, bubble-agitation HTC boost,
    flashpoint/autoignition safety flags, exact `audit_carried_reactions` pass,
    and planted stoichiometry-fault refusal.
  - Dryout-to-zero stress test (`Proposal::Shorten`, zero `NaN`, zero negative mass).
  - Denver `83.4 kPa` altitude boiling shift, `12 psig` pressure cooker
    elevation, sauce reduction BPE + Carreau shear-thinning + scorching risk,
    `kitchen-cooking-oil` matdb card loading, and planted invalid-input refusals.

- Integration suite in `tests/vapor_compression_and_ice.rs`:
  - `G0` First Law (`< 1e-11`) and Second Law (`COP < COP_Carnot`, `S_gen >= 0`)
    across all 5 refrigerants (`R-290`, `R-134a`, `R-404A`, `R-448A` capturing
    ~5.5-6.7 K evaporator glide, and `R-744` subcritical + supercritical gas-cooler
    at 95 F and 115 F with declared 18% near-critical density band).
  - `fs_ports::component::conformance::check_component` and bit-exact durable
    `encode_state`/`decode_state` round-trips for all 4 components.
  - Evaporator frosting-and-defrost transient verifying monotonic `UA_evap`
    degradation, conservative `drain_stream` meltwater removal, temporary box-air
    warming spike, and clean capacity restoration.
  - Planted `+25%` condenser fouling factor raising head pressure and compressor
    current/power.
  - Air-cooled `IceHarvesterComponent` 70 F vs 95 F harvest derating strictly in
    25-40%, continuous 1.5-3.5 kW room waste heat rejection, and `< 1e-11` mass
    conservation.
  - `peng-robinson-refrigerant-v1` card loading and planted negative refusals.

## No-claim boundaries

- Pool-boiling correlations apply to horizontal upward-facing heated surfaces
  under steady or quasi-steady local thermal equilibrium; flow boiling inside
  channels, subcooled nucleate micro-layer collapse dynamics, and near-critical
  boiling (`P > 0.90 P_c`) are outside this crate's claim.
- Pure triglyceride cooking oils are explicitly refused by `BoilingCurveV1`
  rather than assigned synthetic boiling curves; thermal cracking, smoking,
  hydrolysis, TPM polymerization, and flash/autoignition safety boundaries are
  owned by `MultiphaseVesselComponent`.
- Condensation assumes filmwise condensation on smooth surfaces; dropwise
  condensation on superhydrophobic coatings and non-condensable fogging inside
  the bulk gas phase are outside scope.
- `MultiphaseVesselComponent` treats the active liquid bath as well-mixed
  (plus an optional lumped sub-tube fryer cold zone); 3D internal velocity
  fields, boiling-bubble coalescence CFD, and oil-fire combustion flame plumes
  above autoignition are outside scope.
- Peng--Robinson cubic EOS with constant Péneloux shifts carries a declared
  near-critical density error band (up to 18% for `R-744` above `0.95 * T_c`);
  Span--Wagner multi-parameter critical-point scaling, two-phase maldistribution
  across parallel evaporator circuits, and lubricant oil-return miscibility
  dynamics are outside v1 claim.
