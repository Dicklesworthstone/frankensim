//! G1/G3/G4: independent caloric/entropy identities and native equation limits.
use super::*;
use crate::cavitation::{BubbleGas, KellerMiksisBubble, RayleighPlessetBubble};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
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

const RADIUS: f64 = 1e-3;
const TEMPERATURE: f64 = 400.0;
const REFERENCE: f64 = 100_000.0;
const MOLAR_MASS: f64 = 0.028;
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn source() -> Provenance {
    Provenance {
        source: "synthetic caloric bubble identity".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn liquid() -> CavitationState {
    CavitationState::try_new(293.15, 1000.0, 0.001, 0.072, 2000.0, 1500.0).unwrap()
}
fn pg() -> f64 {
    REFERENCE - liquid().saturation_pressure_pa() + 2.0 * liquid().surface_tension_n_m() / RADIUS
}
fn gas(slope: f64, formation: f64, transport_maximum: f64, water: bool) -> ClosedGasParcel {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
    coefficients[3] = slope;
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
    let species = SpeciesId::new("caloric-bubble-fixture").unwrap();
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
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7),
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
        water.then_some(&species),
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            TEMPERATURE,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), pg())
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    ClosedGasParcel::new(model, point, volume(RADIUS).unwrap()).unwrap()
}
fn limits() -> CaloricBubbleLimits {
    CaloricBubbleLimits {
        maximum_gas_liquid_density_ratio: 0.02,
        maximum_vapor_pressure_fraction: 0.1,
        maximum_gas_wall_mach_number: 0.05,
    }
}
fn equation(compressible: bool) -> CaloricBubbleEquation {
    let bubble = BubbleLimits {
        minimum_radius_m: 0.3 * RADIUS,
        maximum_radius_m: 2.0 * RADIUS,
        maximum_wall_mach_number: 0.05,
        maximum_vapor_density_ratio: 0.01,
    };
    if compressible {
        CaloricBubbleEquation::KellerMiksis(KellerMiksisLimits {
            bubble,
            maximum_pressure_ratio: 0.01,
        })
    } else {
        CaloricBubbleEquation::RayleighPlesset(bubble)
    }
}
fn solve() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 220.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}
fn numerics() -> CaloricBubbleNumerics {
    CaloricBubbleNumerics {
        radial: BubbleNumerics {
            relative_tolerance: 1e-9,
            radius_absolute_tolerance_m: 1e-13,
            velocity_absolute_tolerance_m_s: 1e-8,
            initial_step_s: 1e-7,
            minimum_step_s: 1e-13,
            maximum_step_s: 1e-6,
            maximum_attempts: 50_000,
        },
        energy_absolute_tolerance_j: 1e-13,
        energy_balance_tolerance_j: 1e-12,
    }
}
fn drive(time: f64) -> BubblePressureSample {
    let omega = 2.0 * core::f64::consts::PI * 2000.0;
    BubblePressureSample {
        absolute_pressure_pa: REFERENCE + 40_000.0 * (omega * time).sin(),
        time_derivative_pa_s: 40_000.0 * omega * (omega * time).cos(),
    }
}
fn bubble(compressible: bool, slope: f64, formation: f64) -> CaloricBubble {
    CaloricBubble::new(
        liquid(),
        gas(slope, formation, 2000.0, false),
        equation(compressible),
        BubbleState::try_new(0.0, RADIUS, 0.0).unwrap(),
        drive(0.0),
        limits(),
        solve(),
    )
    .unwrap()
}
fn no_heat(_: f64, _: &CaloricGasObservation) -> Result<f64, CaloricBubbleError> {
    Ok(0.0)
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} (tolerance {tolerance:e})"
    );
}

#[test]
fn constant_cv_adiabatic_dynamics_match_both_original_polytropic_equations() {
    for compressible in [false, true] {
        let mut coupled = bubble(compressible, 0.0, 0.0);
        let mut original = coupled.state();
        let gas = BubbleGas::try_new(RADIUS, pg(), 1.4).unwrap();
        let controls = numerics();
        let end = 1.2e-4;
        let step = coupled
            .advance_controlled(end, &controls, &drive, &no_heat, &mut || false)
            .unwrap();
        let native = match coupled.equation() {
            CaloricBubbleEquation::RayleighPlesset(limits) => {
                RayleighPlessetBubble::try_new(liquid(), gas, limits)
                    .unwrap()
                    .advance(&mut original, end, &controls.radial, &|t| {
                        drive(t).absolute_pressure_pa
                    })
                    .unwrap()
            }
            CaloricBubbleEquation::KellerMiksis(limits) => {
                KellerMiksisBubble::try_new(liquid(), gas, limits)
                    .unwrap()
                    .advance(&mut original, end, &controls.radial, &drive)
                    .unwrap()
            }
        };
        near(coupled.state().radius_m(), original.radius_m(), 2e-11);
        near(
            coupled.state().wall_velocity_m_s(),
            original.wall_velocity_m_s(),
            5e-7,
        );
        near(
            step.observation.bubble.gas_partial_pressure_pa,
            native.observation.gas_partial_pressure_pa,
            0.02,
        );
        let expected_t = TEMPERATURE * (RADIUS / coupled.state().radius_m()).powf(1.2);
        near(coupled.gas().state().state().temperature, expected_t, 2e-5);
        assert_eq!(step.supplied_heat_j, 0.0);
        assert!(step.gas_energy_residual_j.abs() <= controls.energy_balance_tolerance_j);
    }
}

#[test]
fn variable_cv_entropy_and_source_integral_are_independent_of_formation_offsets() {
    let slope = 0.003;
    let mut plain = bubble(true, slope, 0.0);
    let mut shifted = bubble(true, slope, 1e14);
    let initial = plain.gas().clone();
    let mass = initial.mass_kg();
    let controls = numerics();
    let a = plain
        .advance_controlled(1.2e-4, &controls, &drive, &no_heat, &mut || false)
        .unwrap();
    let b = shifted
        .advance_controlled(1.2e-4, &controls, &drive, &no_heat, &mut || false)
        .unwrap();
    assert_eq!(plain.state(), shifted.state());
    assert_eq!(a, b);
    assert_eq!(plain.gas().mass_kg(), mass);
    assert_eq!(plain.gas().model(), initial.model());
    assert_eq!(
        plain.gas().point().axes().get("specimen-lot"),
        initial.point().axes().get("specimen-lot")
    );
    let t = plain.gas().state().state().temperature;
    // Cv/R = 2.5 + b T. Integrating ds/R = Cv/R dT/T + dV/V
    // gives this exact adiabatic identity without consulting a native stepper.
    let entropy = 2.5 * (t / TEMPERATURE).ln()
        + slope * (t - TEMPERATURE)
        + (plain.gas().volume_m3() / initial.volume_m3()).ln();
    near(entropy, 0.0, 2e-8);
    let gas_constant = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS;
    let expected = mass
        * gas_constant
        * (2.5 * (t - TEMPERATURE) + 0.5 * slope * (t - TEMPERATURE) * (t + TEMPERATURE));
    near(a.gas_energy_change_j, expected, 2e-15);
    near(
        initial
            .stored_energy_change_to(plain.gas(), || false)
            .unwrap(),
        expected,
        2e-15,
    );
    let mut constant = bubble(true, 0.0, 0.0);
    constant
        .advance_controlled(1.2e-4, &controls, &drive, &no_heat, &mut || false)
        .unwrap();
    assert!((constant.state().radius_m() - plain.state().radius_m()).abs() > 1e-6);
}

#[test]
fn signed_heat_changes_motion_and_closes_original_calorics_and_pressure_work() {
    let controls = numerics();
    let mut hot = bubble(true, 0.001, 0.0);
    let mut cold = hot.clone();
    let initial = hot.gas().clone();
    let held = |_: f64| BubblePressureSample {
        absolute_pressure_pa: REFERENCE,
        time_derivative_pa_s: 0.0,
    };
    let heating = hot
        .advance_controlled(5e-5, &controls, &held, &|_, _| Ok(0.5), &mut || false)
        .unwrap();
    let cooling = cold
        .advance_controlled(5e-5, &controls, &held, &|_, _| Ok(-0.5), &mut || false)
        .unwrap();
    assert!(hot.gas().state().state().temperature > TEMPERATURE);
    assert!(cold.gas().state().state().temperature < TEMPERATURE);
    assert!(hot.state().radius_m() > cold.state().radius_m());
    near(heating.supplied_heat_j, 0.5 * 5e-5, 1e-17);
    near(cooling.supplied_heat_j, -0.5 * 5e-5, 1e-17);
    near(
        initial
            .stored_energy_change_to(hot.gas(), || false)
            .unwrap()
            + heating.work_by_gas_j
            - heating.supplied_heat_j,
        0.0,
        controls.energy_balance_tolerance_j,
    );
    // KM's analytic pressure derivative includes heat even at zero R'.
    let at_rest = bubble(true, 0.001, 0.0);
    let warm = at_rest.observe(held(0.0), 0.5).unwrap();
    let neutral = at_rest.observe(held(0.0), 0.0).unwrap();
    let ratio = at_rest.gas().state().state().specific_gas_constant / warm.gas.cv_j_kg_k;
    let expected = (RADIUS / liquid().liquid_sound_speed_m_s()) * ratio / warm.gas.volume_m3 * 0.5
        / liquid().liquid_density_kg_m3()
        / (RADIUS
            + 4.0 * liquid().dynamic_viscosity_pa_s()
                / liquid().liquid_density_kg_m3()
                / liquid().liquid_sound_speed_m_s());
    near(
        warm.bubble.wall_acceleration_m_s2 - neutral.bubble.wall_acceleration_m_s2,
        expected,
        1e-10,
    );
}

#[test]
fn proposal_drop_cancellation_work_limit_and_endpoint_source_failure_leave_all_state_unchanged() {
    let mut owner = bubble(true, 0.001, 0.0);
    let before = owner.clone();
    let proposal = owner
        .prepare_advance_controlled(1e-5, &numerics(), &drive, &no_heat, &mut || false)
        .unwrap();
    assert_eq!(proposal.step().observation.bubble.state.time_s(), 1e-5);
    drop(proposal);
    assert_eq!(owner, before);
    // A maximum step of 1 us forces accepted prefixes before this late poll.
    let last_heat_time = std::cell::Cell::new(0.0);
    let error = owner
        .advance_controlled(
            1e-4,
            &numerics(),
            &drive,
            &|time, _| {
                last_heat_time.set(time);
                Ok(0.0)
            },
            &mut || last_heat_time.get() > 5e-6,
        )
        .unwrap_err();
    assert_eq!(error, CaloricBubbleError::Cancelled);
    assert_eq!(owner, before);
    let mut limited = numerics();
    limited.radial.maximum_attempts = 1;
    assert!(matches!(
        owner.advance_controlled(1e-4, &limited, &drive, &no_heat, &mut || false),
        Err(CaloricBubbleError::Dynamics(
            CavitationError::StepBudgetExceeded { .. }
        ))
    ));
    assert_eq!(owner, before);
    let mut narrow = CaloricBubble::new(
        liquid(),
        gas(0.001, 0.0, TEMPERATURE + 0.001, false),
        equation(true),
        before.state(),
        drive(0.0),
        limits(),
        solve(),
    )
    .unwrap();
    let snapshot = narrow.clone();
    assert!(matches!(
        narrow.advance_controlled(5e-5, &numerics(), &drive, &|_, _| Ok(0.5), &mut || false),
        Err(CaloricBubbleError::Gas(GasEnergyError::Material(_)))
    ));
    assert_eq!(narrow, snapshot);
}

#[test]
fn noncondensable_material_gas_inertia_mach_vapor_and_energy_scales_are_admitted() {
    let base = bubble(true, 0.0, 0.0);
    assert!(
        CaloricBubble::new(
            liquid(),
            gas(0.0, 0.0, 2000.0, true),
            equation(true),
            base.state(),
            drive(0.0),
            limits(),
            solve()
        )
        .is_err()
    );
    for limit in [
        CaloricBubbleLimits {
            maximum_gas_liquid_density_ratio: 1e-6,
            ..limits()
        },
        CaloricBubbleLimits {
            maximum_vapor_pressure_fraction: 0.001,
            ..limits()
        },
    ] {
        assert!(matches!(
            CaloricBubble::new(
                liquid(),
                base.gas().clone(),
                equation(true),
                base.state(),
                drive(0.0),
                limit,
                solve()
            ),
            Err(CaloricBubbleError::OutsideModel { .. })
        ));
    }
    let fast = BubbleState::try_new(0.0, RADIUS, 30.0).unwrap();
    assert!(matches!(
        CaloricBubble::new(
            liquid(),
            base.gas().clone(),
            equation(true),
            fast,
            drive(0.0),
            limits(),
            solve()
        ),
        Err(CaloricBubbleError::OutsideModel {
            quantity: "gas wall Mach number",
            ..
        })
    ));
    let mut owner = base.clone();
    let mut invalid = numerics();
    invalid.energy_absolute_tolerance_j =
        owner.gas().mass_kg() * solve().energy_tolerance_j_kg * 0.5;
    assert!(matches!(
        owner.advance_controlled(1e-5, &invalid, &drive, &no_heat, &mut || false),
        Err(CaloricBubbleError::InvalidInput(_))
    ));
    assert_eq!(owner, base);
}
