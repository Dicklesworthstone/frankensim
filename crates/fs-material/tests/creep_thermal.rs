//! Source-bound zero-expansion thermocreep; coefficients are manufactured, not measured.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::creep::source::*;
use fs_material::creep::thermal::{
    FeedbackConfig, FeedbackError, InitialState, ZeroExpansionCaloricPoint,
};
use fs_material::creep::{CreepState, CreepTolerances};
use fs_material::phase::{
    HeatCapacitySourceQuery, PhaseStateError, ResolvedHeatCapacityCurve, SolidLiquidPhase,
    UniformEnthalpyStep,
};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use std::cell::Cell;

const MASS: f64 = 1e-3;
fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, value), (500.0, value)],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
            500.0,
        ),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured thermocreep fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn card(cp: f64, varying_density: bool) -> MaterialCard {
    let mut properties = vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1.2e9),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        scalar(CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, 2e6),
        scalar(CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, 1e-4),
        scalar(CREEP_STRESS_EXPONENT_PROPERTY, Dims([0; 6]), 3.0),
        scalar(
            CREEP_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            30_000.0,
        ),
        scalar(
            CREEP_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            400.0,
        ),
        scalar("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0]), cp),
        scalar("density", Dims([-3, 1, 0, 0, 0, 0]), 1000.0),
    ];
    if varying_density {
        properties[8].value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 1000.0), (400.0, 900.0), (500.0, 1000.0)],
            dims: Dims([-3, 1, 0, 0, 0, 0]),
        };
        properties[8].interpolation = InterpolationPolicy::LinearInside;
    }
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermocreep solid".into(),
            phase: "solid".into(),
            process: "manufactured reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn source(card: &MaterialCard) -> NortonCreepSource<'_> {
    NortonCreepSource::new(
        card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.05,
        },
        || false,
    )
    .unwrap()
}
fn caloric(card: &MaterialCard, phase: SolidLiquidPhase) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
            density: &PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0])),
            temperature_axis: "T",
            lower: &point(300.0),
            upper: &point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        phase,
        0.0,
        0.001,
    )
    .unwrap()
}
fn initial(caloric: &ResolvedHeatCapacityCurve) -> InitialState {
    InitialState {
        strain: [0.0; 6],
        history: CreepState {
            creep_strain: [0.0; 6],
            equivalent_creep_strain: 0.0,
        },
        specific_enthalpy_j_kg: caloric
            .state_at_temperature(400.0)
            .unwrap()
            .specific_enthalpy_j_kg(),
        mass_kg: MASS,
        surface_area_m2: 0.001,
        time_s: 12.0,
    }
}
fn config() -> FeedbackConfig {
    FeedbackConfig {
        max_iterations: 100,
        temperature_tolerance_k: 1e-9,
        energy_tolerance_j: 1e-9,
        stress_tolerance_pa: 1e-5,
        relaxation: 0.75,
    }
}
fn strain() -> [f64; 6] {
    [0.02, -0.01, -0.01, 0.0, 0.0, 0.0]
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

#[test]
fn g1_physical_creep_heat_accelerates_flow_and_excludes_numerical_loss() {
    let card = card(5.0, false);
    let caloric = caloric(&card, SolidLiquidPhase::Solid);
    let old =
        ZeroExpansionCaloricPoint::new(&card, source(&card), &caloric, initial(&caloric), || false)
            .unwrap();
    let next = old
        .advance_adiabatic(strain(), 0.1, config(), || false)
        .unwrap();
    let frozen = old
        .source()
        .at_temperature(400.0, 0.1)
        .unwrap()
        .law()
        .advance_with_energy(
            old.strain(),
            old.history(),
            &strain(),
            CreepTolerances {
                stress_pa: 1e-5,
                energy_j_m3: 0.001,
            },
        )
        .unwrap();
    assert!(next.next.caloric_state().temperature_k() > 420.0);
    assert!(next.next.history().equivalent_creep_strain > frozen.state.equivalent_creep_strain);
    assert!(next.numerical_loss_j > 0.0);
    let thermal = MASS
        * (next.next.caloric_state().specific_enthalpy_j_kg()
            - old.caloric_state().specific_enthalpy_j_kg());
    close(thermal, next.physical_heat_j, config().energy_tolerance_j);
    assert!((thermal - next.physical_heat_j - next.numerical_loss_j).abs() > 1e-3);
    assert!(next.combined_residual_j.abs() < config().energy_tolerance_j);
    assert!(next.iterations > 1);
    assert_eq!(old.time_s(), 12.0);
    assert_eq!(old.history().equivalent_creep_strain, 0.0);
    close(next.next.time_s(), 12.1, 0.0);
}

#[test]
fn g3_caloric_capacity_changes_temperature_and_same_native_creep_history() {
    let low_card = card(5.0, false);
    let high_card = card(100.0, false);
    let low_caloric = caloric(&low_card, SolidLiquidPhase::Solid);
    let high_caloric = caloric(&high_card, SolidLiquidPhase::Solid);
    let low = ZeroExpansionCaloricPoint::new(
        &low_card,
        source(&low_card),
        &low_caloric,
        initial(&low_caloric),
        || false,
    )
    .unwrap()
    .advance_adiabatic(strain(), 0.1, config(), || false)
    .unwrap();
    let high = ZeroExpansionCaloricPoint::new(
        &high_card,
        source(&high_card),
        &high_caloric,
        initial(&high_caloric),
        || false,
    )
    .unwrap()
    .advance_adiabatic(strain(), 0.1, config(), || false)
    .unwrap();
    assert!(low.next.caloric_state().temperature_k() > high.next.caloric_state().temperature_k());
    assert!(
        low.next.history().equivalent_creep_strain > high.next.history().equivalent_creep_strain
    );
    for step in [&low, &high] {
        let bound = step
            .next
            .source()
            .at_temperature(step.next.caloric_state().temperature_k(), 0.1)
            .unwrap();
        let expected = bound
            .law()
            .advance_with_energy(
                &[0.0; 6],
                &initial(&low_caloric).history,
                &strain(),
                CreepTolerances {
                    stress_pa: 1e-5,
                    energy_j_m3: 0.001,
                },
            )
            .unwrap();
        assert_eq!(&expected.state, step.next.history());
        assert_eq!(&expected.stress_pa, step.next.stress_pa());
    }
}

#[test]
fn g0_external_transport_retains_reports_and_trials_never_advance_old_history() {
    struct Report(usize);
    let card = card(5.0, false);
    let caloric = caloric(&card, SolidLiquidPhase::Solid);
    let old =
        ZeroExpansionCaloricPoint::new(&card, source(&card), &caloric, initial(&caloric), || false)
            .unwrap();
    let mut calls = 0;
    let next = old
        .advance(
            strain(),
            0.1,
            config(),
            |input| {
                calls += 1;
                assert_eq!(input.initial, *old.caloric_state());
                assert_eq!(input.duration_s, 0.1);
                assert_eq!(input.mass_kg, MASS);
                let external_heat_j = 0.025;
                let state = input.curve.advance_specific_energy(
                    input.initial,
                    (input.internal_heat_j + external_heat_j) / input.mass_kg,
                )?;
                Ok::<_, PhaseStateError>(UniformEnthalpyStep {
                    state,
                    external_heat_j,
                    energy_residual_tolerance_j: 1e-10,
                    report: Report(calls),
                })
            },
            || false,
        )
        .unwrap();
    assert_eq!(next.transport.report.0, calls);
    assert_eq!(next.iterations, calls);
    let thermal = MASS
        * (next.next.caloric_state().specific_enthalpy_j_kg()
            - old.caloric_state().specific_enthalpy_j_kg());
    close(
        thermal,
        next.physical_heat_j + 0.025,
        config().energy_tolerance_j,
    );
    assert_eq!(old.history().equivalent_creep_strain, 0.0);
    assert_eq!(old.time_s(), 12.0);
}

#[test]
fn g0_cancel_failure_and_exhaustion_publish_no_aged_state() {
    let card = card(5.0, false);
    let caloric = caloric(&card, SolidLiquidPhase::Solid);
    let old =
        ZeroExpansionCaloricPoint::new(&card, source(&card), &caloric, initial(&caloric), || false)
            .unwrap();
    let preserved = (
        old.strain().to_owned(),
        old.history().clone(),
        *old.caloric_state(),
        old.time_s(),
    );
    assert!(matches!(
        old.advance_adiabatic(strain(), 0.1, config(), || true),
        Err(FeedbackError::Cancelled)
    ));
    assert!(matches!(
        old.advance_adiabatic(
            strain(),
            0.1,
            FeedbackConfig {
                max_iterations: 1,
                ..config()
            },
            || false
        ),
        Err(FeedbackError::NotConverged { .. })
    ));
    let calls = Cell::new(0);
    let cancelled = old.advance(
        strain(),
        0.1,
        config(),
        |input| {
            calls.set(calls.get() + 1);
            Ok::<_, PhaseStateError>(UniformEnthalpyStep {
                state: input.curve.advance_specific_energy(
                    input.initial,
                    input.internal_heat_j / input.mass_kg,
                )?,
                external_heat_j: 0.0,
                energy_residual_tolerance_j: 1e-10,
                report: (),
            })
        },
        || calls.get() > 0,
    );
    assert!(matches!(cancelled, Err(FeedbackError::Cancelled)));
    let bad = old.advance(
        strain(),
        0.1,
        config(),
        |input| {
            Ok::<_, PhaseStateError>(UniformEnthalpyStep {
                state: input.initial,
                external_heat_j: 0.0,
                energy_residual_tolerance_j: 1e-10,
                report: (),
            })
        },
        || false,
    );
    assert!(matches!(bad, Err(FeedbackError::Invalid(_))));
    assert_eq!(
        (
            old.strain().to_owned(),
            old.history().clone(),
            *old.caloric_state(),
            old.time_s()
        ),
        preserved
    );
}

#[test]
fn g0_foreign_phase_and_variable_density_are_refused_before_feedback() {
    let original = card(5.0, false);
    let different = card(10.0, false);
    let foreign = caloric(&different, SolidLiquidPhase::Solid);
    assert!(
        ZeroExpansionCaloricPoint::new(
            &original,
            source(&original),
            &foreign,
            initial(&foreign),
            || false
        )
        .is_err()
    );
    let liquid = caloric(&original, SolidLiquidPhase::Liquid);
    assert!(
        ZeroExpansionCaloricPoint::new(
            &original,
            source(&original),
            &liquid,
            initial(&liquid),
            || false
        )
        .is_err()
    );
    let varying = card(5.0, true);
    let varying_caloric = caloric(&varying, SolidLiquidPhase::Solid);
    assert!(
        ZeroExpansionCaloricPoint::new(
            &varying,
            source(&varying),
            &varying_caloric,
            initial(&varying_caloric),
            || false
        )
        .is_err()
    );
}
