//! Original IF97 water plus a manufactured NASA-9 variable-Cv gas.
//! The high wet-mixture transport bounds below isolate the homogeneous
//! closure's equations; they are NOT empirical water/air transport data.
use super::super::ConductingVaporBubble;
use super::super::tests::{
    RADIUS, TEMPERATURE, equilibrium_pressure as pure_pressure, model as pure_model, numerics,
};
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    required_nasa9_parameter_dims,
};

const MOLAR_MASS: f64 = 0.028;
const GAS_PRESSURE: f64 = 30_000.0;

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn source() -> Provenance {
    Provenance {
        source: "manufactured conserved gas-vapor caloric identity".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn inventory(partial_pressure: f64, formation: f64, transport_maximum: f64) -> ConservedBubbleGas {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
    coefficients[3] = 0.0005;
    coefficients[7] = formation;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let species = SpeciesId::new("conserved-gas-vapor-fixture").unwrap();
    let thermodynamics = Nasa9StandardStateModelV1::new(
        species.clone(),
        MolarMass::new(MOLAR_MASS),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("synthetic").unwrap(),
        )
        .unwrap(),
        vec![
            Nasa9RegionV1::from_card(ConstitutiveModelCard {
                law: LawId(NASA9_LAW_ID_V1.into()),
                law_version: NASA9_LAW_VERSION_V1,
                parameters,
                state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
                sources: vec![],
                provenance: source(),
            })
            .unwrap(),
        ],
    )
    .unwrap();
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 2e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.03),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (2000.0, 2.0 * value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        transport_maximum,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 1e-8, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: source(),
            })
            .unwrap();
    }
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![thermodynamics],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![SpeciesTransportMaterial {
            species: species.clone(),
            material: MaterialCard::assemble(
                MaterialStateId {
                    chemistry: "synthetic bubble gas".into(),
                    phase: "gas".into(),
                    process: "test".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        None,
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            TEMPERATURE,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            quantity(QuantityKind::Pressure),
            partial_pressure,
        )
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    let volume = (4.0 * PI / 3.0) * RADIUS * RADIUS * RADIUS;
    let parcel = ClosedGasParcel::new(model, point, volume).unwrap();
    ConservedBubbleGas::new(
        parcel,
        GasEnergySolve {
            minimum_temperature_k: TEMPERATURE - 2.0,
            maximum_temperature_k: TEMPERATURE + 2.0,
            energy_tolerance_j_kg: 1e-8,
            maximum_iterations: 128,
        },
    )
    .unwrap()
}

fn limits() -> GasVaporLimits {
    GasVaporLimits {
        transport: GasVaporTransportBounds {
            minimum_temperature_k: TEMPERATURE - 2.0,
            maximum_temperature_k: TEMPERATURE + 2.0,
            minimum_total_pressure_pa: 1000.0,
            maximum_total_pressure_pa: 1e7,
            minimum_vapor_mass_fraction: 0.0,
            maximum_vapor_mass_fraction: 1.0,
            minimum_thermal_conductivity_w_m_k: 10.0,
            minimum_binary_diffusivity_m2_s: 0.01,
        },
        maximum_total_density_ratio: 0.02,
        maximum_mixture_mach_number: 0.05,
        maximum_pressure_relaxation_ratio: 0.05,
        maximum_internal_temperature_difference_k: 0.1,
        maximum_relative_composition_gradient: 0.05,
        maximum_thermal_response_ratio: 0.05,
        maximum_diffusion_response_ratio: 0.05,
        maximum_imposed_rate_per_s: 0.0,
    }
}

fn loaded(offset: f64, velocity: f64, formation: f64) -> ConductingVaporBubble {
    pure_model(offset, velocity)
        .with_noncondensable_gas(inventory(GAS_PRESSURE, formation, 2000.0), limits())
        .unwrap()
}

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn gas_supported_equilibrium_and_loaded_species_momentum_first_law_are_coherent() {
    let pressure = pure_pressure() + GAS_PRESSURE;
    let stationary = loaded(0.0, 0.0, 0.0);
    let equilibrium = stationary.observe(pressure, TEMPERATURE).unwrap();
    near(equilibrium.bubble.liquid_acceleration_m_s2, 0.0, 0.0, 1e-8);
    assert_eq!(equilibrium.bubble.vapor_mass_rate_kg_s, 0.0);
    assert_eq!(equilibrium.bubble.temperature_rate_k_s, 0.0);
    assert!(equilibrium.gas.unwrap().mass_kg > 0.0);
    for offset in [-0.001, 0.001] {
        for velocity in [-0.005, 0.005] {
            let bubble = loaded(offset, velocity, 0.0);
            let frame = bubble.observe(pressure, TEMPERATURE + offset).unwrap();
            let o = frame.bubble;
            let g = frame.gas.unwrap();
            let water = bubble.model().source().state();
            let area = 4.0 * PI * RADIUS * RADIUS;
            let rho_l = water.liquid_density_kg_m3();
            let rho_v = o.vapor.vapor_density_kg_m3();
            near(
                o.bubble_pressure_pa,
                o.vapor.pressure_pa() + g.partial_pressure_pa,
                1e-15,
                0.0,
            );
            near(
                o.wall_velocity_m_s,
                velocity + o.vapor_mass_flux_kg_m2_s / rho_l,
                1e-14,
                1e-18,
            );
            let vapor_energy_rate = o.vapor.vapor_internal_energy_j_kg() * o.vapor_mass_rate_kg_s
                + o.vapor_mass_kg
                    * o.vapor
                        .vapor_internal_energy_temperature_derivative_j_kg_k()
                    * o.temperature_rate_k_s;
            near(
                vapor_energy_rate + g.energy_rate_w,
                o.heat_power_w + o.vapor.liquid_enthalpy_j_kg() * o.vapor_mass_rate_kg_s
                    - o.bubble_pressure_pa * o.volume_rate_m3_s,
                3e-12,
                1e-16,
            );
            let j = o.vapor_mass_flux_kg_m2_s;
            near(
                o.vapor_recoil_pressure_pa,
                j * j * (1.0 / rho_v - 1.0 / rho_l),
                3e-15,
                1e-24,
            );
            // Inserting rho_mix in this species recoil would lose the dry
            // component's zero interfacial flux and change momentum transfer.
            assert!(
                (o.vapor_recoil_pressure_pa
                    - j * j * (1.0 / g.mixture_density_kg_m3 - 1.0 / rho_l))
                    .abs()
                    > 1e-12
            );
            let momentum_left = RADIUS * o.liquid_acceleration_m_s2
                + 2.0 * o.wall_velocity_m_s * velocity
                - 0.5 * velocity * velocity;
            let momentum_right = (o.bubble_pressure_pa
                - pressure
                - 2.0 * water.surface_tension_n_m() / RADIUS
                - 4.0 * water.dynamic_viscosity_pa_s() * velocity / RADIUS
                + o.vapor_recoil_pressure_pa)
                / rho_l;
            near(momentum_left, momentum_right, 5e-13, 1e-15);
            near(o.liquid_volume_flow_m3_s, area * velocity, 2e-15, 0.0);
            near(
                o.liquid_volume_flow_rate_m3_s2,
                area * (o.liquid_acceleration_m_s2 + 2.0 * o.wall_velocity_m_s * velocity / RADIUS),
                2e-15,
                1e-24,
            );
            near(
                g.heat_power_w,
                g.energy_rate_w + g.pressure_work_power_w,
                0.0,
                0.0,
            );
            assert_eq!(o.heat_power_w, frame.liquid.interface_heat_power_w);
        }
    }
}

#[test]
fn accepted_mixture_conserves_inventory_and_original_gas_energy_under_formation_shifts() {
    let mut original = loaded(0.001, 0.0, 0.0);
    let mut shifted = loaded(0.001, 0.0, 1e12);
    let reference_mass = original.gas().unwrap().mass_kg();
    let gas_constant = original.gas().unwrap().gas_constant_j_kg_k();
    let pressure = pure_pressure() + GAS_PRESSURE;
    let controls = numerics();
    let mut steps = Vec::new();
    for bubble in [&mut original, &mut shifted] {
        let step = bubble
            .advance_controlled(
                2e-6,
                &controls,
                &|_| pressure,
                &|_| TEMPERATURE + 0.001,
                &mut || false,
            )
            .unwrap();
        let temperature = bubble.state().temperature_k();
        // Independent analytic integral for the authored variable-Cv NASA law.
        let expected = reference_mass
            * gas_constant
            * (temperature - TEMPERATURE)
            * (2.5 + 0.00025 * (temperature + TEMPERATURE));
        near(
            step.bubble.noncondensable_energy_change_j,
            expected,
            5e-13,
            1e-23,
        );
        near(
            step.integrated_noncondensable_energy_j,
            expected,
            0.0,
            controls.bubble.energy_balance_tolerance_j,
        );
        near(
            step.bubble.noncondensable_energy_change_j + step.noncondensable_pressure_work_j,
            step.noncondensable_heat_j,
            0.0,
            controls.bubble.energy_balance_tolerance_j,
        );
        near(
            step.bubble.contents_energy_change_j,
            step.bubble.vapor_energy_change_j + expected,
            5e-13,
            1e-23,
        );
        assert_eq!(
            bubble.gas().unwrap().mass_kg().to_bits(),
            reference_mass.to_bits()
        );
        assert_eq!(
            step.gas.unwrap().mass_kg.to_bits(),
            reference_mass.to_bits()
        );
        assert_eq!(
            bubble
                .gas()
                .unwrap()
                .reference_parcel()
                .state()
                .state()
                .temperature,
            TEMPERATURE
        );
        assert!(
            step.noncondensable_energy_residual_j.abs()
                <= controls.bubble.energy_balance_tolerance_j
        );
        assert!(step.total_energy_residual_j.abs() <= controls.total_energy_balance_tolerance_j);
        assert!(step.liquid_energy_residual_j.abs() <= controls.liquid_energy_balance_tolerance_j);
        steps.push(step);
    }
    assert_eq!(original.state(), shifted.state());
    assert_eq!(
        original.liquid_temperatures_k(),
        shifted.liquid_temperatures_k()
    );
    assert_eq!(steps[0], steps[1]);
}

#[test]
fn dry_inventory_cushions_a_compression_and_vanishing_inventory_recovers_the_pure_limit() {
    let with_gas = loaded(0.0, -0.01, 0.0);
    let without_gas = pure_model(0.0, -0.01);
    let pressure = pure_pressure() + GAS_PRESSURE;
    let loaded_frame = with_gas.observe(pressure, TEMPERATURE).unwrap();
    let pure_frame = without_gas.observe(pressure, TEMPERATURE).unwrap();
    assert!(
        loaded_frame.bubble.liquid_acceleration_m_s2 > pure_frame.bubble.liquid_acceleration_m_s2
    );
    assert!(loaded_frame.gas.unwrap().pressure_work_power_w < 0.0);
    let trace = pure_model(0.0, 0.0)
        .with_noncondensable_gas(inventory(0.01, 0.0, 2000.0), limits())
        .unwrap();
    let trace_frame = trace.observe(pure_pressure() + 0.01, TEMPERATURE).unwrap();
    near(
        trace_frame.gas.unwrap().mixture_sound_speed_m_s,
        trace_frame.bubble.vapor.vapor_sound_speed_m_s(),
        1e-6,
        0.0,
    );
    let mut pure = pure_model(0.001, 0.0);
    let step = pure
        .advance_controlled(
            1e-6,
            &numerics(),
            &|_| pure_pressure(),
            &|_| TEMPERATURE + 0.001,
            &mut || false,
        )
        .unwrap();
    assert!(step.gas.is_none());
    assert_eq!(step.integrated_noncondensable_energy_j, 0.0);
    assert_eq!(step.noncondensable_pressure_work_j, 0.0);
    assert_eq!(step.noncondensable_heat_j, 0.0);
    assert_eq!(step.bubble.noncondensable_energy_change_j, 0.0);
    assert_eq!(
        step.bubble.contents_energy_change_j,
        step.bubble.vapor_energy_change_j
    );
    assert_eq!(
        step.bubble.observation.bubble_pressure_pa,
        step.bubble.observation.vapor.pressure_pa()
    );
}

#[test]
fn wet_transport_domains_and_failed_mixing_estimates_refuse_without_changing_state() {
    let pressure = pure_pressure() + GAS_PRESSURE;
    for selection in 0..4 {
        let mut configured = limits();
        match selection {
            0 => configured.transport.minimum_thermal_conductivity_w_m_k = 1e-12,
            1 => configured.transport.minimum_binary_diffusivity_m2_s = 1e-12,
            2 => configured.maximum_internal_temperature_difference_k = 1e-12,
            _ => configured.maximum_imposed_rate_per_s = 1e10,
        }
        let mut bubble = pure_model(0.001, 0.0)
            .with_noncondensable_gas(inventory(GAS_PRESSURE, 0.0, 2000.0), configured)
            .unwrap();
        let original_state = bubble.state();
        let original_profile = bubble.liquid_temperatures_k().to_vec();
        let original_mass = bubble.gas().unwrap().mass_kg();
        assert!(
            bubble
                .advance_controlled(
                    1e-6,
                    &numerics(),
                    &|_| pressure,
                    &|_| TEMPERATURE + 0.001,
                    &mut || false
                )
                .is_err()
        );
        assert_eq!(bubble.state(), original_state);
        assert_eq!(bubble.liquid_temperatures_k(), original_profile);
        assert_eq!(bubble.gas().unwrap().mass_kg(), original_mass);
    }
    let mut restricted = limits();
    restricted.transport.maximum_total_pressure_pa = 100_000.0;
    assert!(matches!(
        pure_model(0.0, 0.0)
            .with_noncondensable_gas(inventory(GAS_PRESSURE, 0.0, 2000.0), restricted),
        Err(ConductingVaporBubbleError::OutsideDomain { .. })
    ));
    let mut restricted = limits();
    restricted.transport.minimum_vapor_mass_fraction = 0.99;
    assert!(matches!(
        pure_model(0.0, 0.0)
            .with_noncondensable_gas(inventory(GAS_PRESSURE, 0.0, 2000.0), restricted),
        Err(ConductingVaporBubbleError::OutsideDomain { .. })
    ));
}

#[test]
fn cancellation_dropped_candidate_and_original_endpoint_source_refusal_are_atomic() {
    let pressure = pure_pressure() + GAS_PRESSURE;
    let mut bubble = loaded(0.001, 0.0, 0.0);
    let original_state = bubble.state();
    let original_profile = bubble.liquid_temperatures_k().to_vec();
    let original_inventory = bubble.gas().unwrap().clone();
    let mut complete_polls = 0;
    {
        let proposal = bubble
            .prepare_advance_controlled(
                1e-6,
                &numerics(),
                &|_| pressure,
                &|_| TEMPERATURE + 0.001,
                &mut || {
                    complete_polls += 1;
                    false
                },
            )
            .unwrap();
        assert_ne!(proposal.step().bubble.observation.state, original_state);
    }
    assert_eq!(bubble.state(), original_state);
    assert_eq!(bubble.liquid_temperatures_k(), original_profile);
    assert_eq!(bubble.gas().unwrap(), &original_inventory);
    let mut polls = 0;
    assert!(matches!(
        bubble.advance_controlled(
            1e-6,
            &numerics(),
            &|_| pressure,
            &|_| TEMPERATURE + 0.001,
            &mut || {
                polls += 1;
                polls >= complete_polls
            }
        ),
        Err(ConductingVaporBubbleError::Cancelled)
    ));
    assert_eq!(bubble.state(), original_state);
    assert_eq!(bubble.liquid_temperatures_k(), original_profile);
    assert_eq!(bubble.gas().unwrap(), &original_inventory);
    let mut refused = pure_model(0.001, 0.0)
        .with_noncondensable_gas(inventory(GAS_PRESSURE, 0.0, TEMPERATURE + 1e-8), limits())
        .unwrap();
    let before = refused.state();
    assert!(matches!(
        refused.advance_controlled(
            1e-6,
            &numerics(),
            &|_| pressure,
            &|_| TEMPERATURE + 0.001,
            &mut || false
        ),
        Err(ConductingVaporBubbleError::Gas(_))
    ));
    assert_eq!(refused.state(), before);
    assert_eq!(refused.liquid_temperatures_k(), original_profile);
}

#[test]
fn inventory_rebinding_and_thermal_geometry_mismatch_are_rejected() {
    let bubble = loaded(0.0, 0.0, 0.0);
    assert!(matches!(
        bubble.with_noncondensable_gas(inventory(GAS_PRESSURE, 0.0, 2000.0), limits()),
        Err(ConductingVaporBubbleError::InvalidInput(_))
    ));
    let mut advanced = pure_model(0.001, 0.0);
    advanced
        .advance_controlled(
            1e-6,
            &numerics(),
            &|_| pure_pressure(),
            &|_| TEMPERATURE + 0.001,
            &mut || false,
        )
        .unwrap();
    assert!(matches!(
        advanced.with_noncondensable_gas(inventory(GAS_PRESSURE, 0.0, 2000.0), limits()),
        Err(ConductingVaporBubbleError::InvalidInput(_))
    ));
}
