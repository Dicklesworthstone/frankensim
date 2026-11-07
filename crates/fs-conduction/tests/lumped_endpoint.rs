//! Native implicit source feedback retains the original enthalpy/phase owner.

use std::cell::Cell;

use fs_blake3::ContentHash;
use fs_conduction::ConductionError;
use fs_conduction::lumped::{
    LUMPED_ENTHALPY_BISECTION_ITERATIONS, LumpedThermalEnvironment, LumpedThermalTransport,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};

fn knot(h: f64, t: f64, liquid: f64) -> EnthalpyPhaseKnot {
    EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: liquid,
        bulk_density_kg_m3: 1000.0,
    }
}

fn sensible_curve(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x62; 32]),
        SolidLiquidPhase::Solid,
        vec![
            knot(datum, 300.0, 0.0),
            knot(datum + 100_000.0, 1300.0, 0.0),
        ],
    )
    .unwrap()
}

fn environment() -> LumpedThermalEnvironment {
    LumpedThermalEnvironment {
        temperature_k: 300.0,
        radiation_temperature_k: 300.0,
        convection_w_per_m2_k: 0.0,
        transport: LumpedThermalTransport::try_declared(1000.0, 0.0).unwrap(),
        enthalpy_tolerance_j_kg: 1e-9,
        maximum_thermal_residual_j: 1e-9,
    }
}

fn input(
    curve: &EquilibriumEnthalpyPhaseCurve,
    h: f64,
    duration: f64,
    heat: f64,
) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: 2.0,
        volume_m3: 0.002,
        surface_area_m2: 0.1,
        internal_heat_j: heat,
        duration_s: duration,
    }
}

#[test]
fn constant_source_adds_to_prepared_heat_once_with_bounded_work() {
    let curve = sensible_curve(0.0);
    let endpoint = environment()
        .prepare_endpoint(input(&curve, 0.0, 2.0, 30.0))
        .unwrap();
    let calls = Cell::new(0);
    let step = endpoint
        .solve_with_internal_power(
            0.0,
            |_, _| {
                calls.set(calls.get() + 1);
                Ok(50.0)
            },
            || Ok(()),
        )
        .unwrap();
    // Two kg at 100 J/(kg K): 30 J baseline + 100 J additional heat.
    assert!((step.state.temperature_k() - 300.65).abs() < 1e-10);
    assert!((step.endpoint.stored_energy_change_j - 130.0).abs() <= 1e-9);
    assert_eq!(step.additional_internal_heat_j, 100.0);
    assert!(step.endpoint.residual_j.abs() <= 1e-9);
    assert_eq!(step.source_evaluations, calls.get());
    assert!(step.source_evaluations <= LUMPED_ENTHALPY_BISECTION_ITERATIONS + 2);

    let cooling = environment()
        .prepare_endpoint(input(&curve, 10_000.0, 2.0, 30.0))
        .unwrap()
        .solve_with_internal_power(0.0, |_, _| Ok(-50.0), || Ok(()))
        .unwrap();
    assert!((cooling.state.temperature_k() - 399.65).abs() < 1e-10);
    assert_eq!(cooling.additional_internal_heat_j, -100.0);
    assert!((cooling.endpoint.stored_energy_change_j + 70.0).abs() <= 1e-9);
}

#[test]
fn temperature_feedback_uses_implicit_endpoint_and_explicit_step_bound() {
    let curve = sensible_curve(0.0);
    let endpoint = environment()
        .prepare_endpoint(input(&curve, 0.0, 10.0, 0.0))
        .unwrap();
    let feedback = 0.2 / 200.0; // dP/dE, not dP/dT or dP/dh.
    let limit = endpoint
        .maximum_step_s_for_internal_power(feedback, || Ok(()))
        .unwrap();
    assert!(limit < 500.0 && limit > 499.999_999_999);
    let extreme = endpoint
        .maximum_step_s_for_internal_power(f64::MAX, || Ok(()))
        .unwrap();
    assert!(extreme.is_finite() && extreme > 0.0);
    assert!(
        endpoint
            .maximum_step_s_for_internal_power(0.0, || Ok(()))
            .unwrap()
            .is_infinite()
    );
    let step = endpoint
        .solve_with_internal_power(feedback, |_, t| Ok(20.0 + 0.2 * (t - 300.0)), || Ok(()))
        .unwrap();
    let expected_energy = 20.0 * 10.0 / (1.0 - feedback * 10.0);
    assert!((step.endpoint.stored_energy_change_j - expected_energy).abs() < 1e-8);
    assert!((step.additional_internal_heat_j - expected_energy).abs() < 1e-8);
    assert!(
        endpoint
            .solve_with_internal_power(
                0.051,
                |_, _| panic!("source must not run beyond its admitted step"),
                || Ok(())
            )
            .is_err()
    );
}

#[test]
fn driven_source_shares_native_convection_and_radiation_balance() {
    let curve = sensible_curve(0.0);
    let mut bath = environment();
    bath.convection_w_per_m2_k = 5.0;
    bath.transport = LumpedThermalTransport::try_declared(1000.0, 0.8).unwrap();
    let endpoint = bath
        .prepare_endpoint(input(&curve, 50_000.0, 10.0, 100.0))
        .unwrap();
    // Choose an independent analytic endpoint at 900 K from an 800 K start.
    let radiation_out = 0.8 * 5.670_374_419e-8 * 0.1 * (900.0_f64.powi(4) - 300.0_f64.powi(4));
    let power = 200.0 * 100.0 / 10.0 + 5.0 * 0.1 * 600.0 + radiation_out - 10.0;
    let step = endpoint
        .solve_with_internal_power(0.0, |_, _| Ok(power), || Ok(()))
        .unwrap();
    assert!((step.state.temperature_k() - 900.0).abs() < 1e-9);
    assert!((step.endpoint.convection_heat_j + 3000.0).abs() < 1e-8);
    assert!((step.endpoint.radiation_heat_j + 10.0 * radiation_out).abs() < 1e-8);
    assert!(step.endpoint.residual_j.abs() <= bath.maximum_thermal_residual_j);
    let sampled = endpoint
        .evaluate_specific_enthalpy(
            step.state.specific_enthalpy_j_kg(),
            bath.temperature_k,
            100.0 + step.additional_internal_heat_j,
        )
        .unwrap();
    assert_eq!(sampled.convection_heat_j, step.endpoint.convection_heat_j);
    assert_eq!(sampled.radiation_heat_j, step.endpoint.radiation_heat_j);
    assert!((sampled.residual_j - step.endpoint.residual_j).abs() < 1e-10);
}

#[test]
fn source_retains_enthalpy_coordinate_inside_latent_plateau() {
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x63; 32]),
        vec![
            knot(0.0, 300.0, 0.0),
            knot(20_000.0, 500.0, 0.0),
            knot(70_000.0, 500.0, 1.0),
            knot(90_000.0, 700.0, 1.0),
        ],
    )
    .unwrap();
    let endpoint = environment()
        .prepare_endpoint(input(&curve, 20_000.0, 1.0, 0.0))
        .unwrap();
    let step = endpoint
        .solve_with_internal_power(
            0.1 / 2.0,
            |h, _| Ok(10_000.0 + 0.1 * (h - 20_000.0)),
            || Ok(()),
        )
        .unwrap();
    let expected_h = 20_000.0 + 10_000.0 / (2.0 - 0.1);
    assert!((step.state.specific_enthalpy_j_kg() - expected_h).abs() < 1e-8);
    assert_eq!(step.state.temperature_k(), 500.0);
    assert!((step.state.liquid_mass_fraction() - (expected_h - 20_000.0) / 50_000.0).abs() < 1e-12);
}

#[test]
fn exact_equilibrium_is_preserved_and_cancelled_proposals_publish_nothing() {
    let curve = sensible_curve(0.0);
    let original = input(&curve, 10_000.0, 1.0, 0.0);
    let endpoint = environment().prepare_endpoint(original).unwrap();
    let equilibrium = endpoint
        .solve_with_internal_power(0.0, |_, _| Ok(0.0), || Ok(()))
        .unwrap();
    assert_eq!(equilibrium.state, original.initial);
    assert_eq!(equilibrium.endpoint.residual_j, 0.0);
    assert_eq!(equilibrium.source_evaluations, 1);
    let checks = Cell::new(0);
    let calls = Cell::new(0);
    let cancelled = endpoint.solve_with_internal_power(
        0.0,
        |_, _| {
            calls.set(calls.get() + 1);
            Ok(50.0)
        },
        || {
            checks.set(checks.get() + 1);
            if checks.get() == 4 {
                Err(ConductionError::Cancelled {
                    stage: "driven-source-test",
                    at: calls.get(),
                })
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(cancelled, Err(ConductionError::Cancelled { .. })));
    assert_eq!(calls.get(), 1);
    assert_eq!(
        endpoint
            .solve_with_internal_power(0.0, |_, _| Ok(0.0), || Ok(()))
            .unwrap(),
        equilibrium
    );
}

#[test]
fn either_residual_budget_refuses_unrepresentable_energy() {
    let curve = sensible_curve(1e16);
    // Adjacent enthalpy values differ by 2 J/kg, hence 4 J for this body.
    // A 2 J source cannot be represented by any endpoint on this datum.
    for (h_tolerance, energy_tolerance) in [(10.0, 1.0), (0.1, 10.0)] {
        let mut bath = environment();
        bath.enthalpy_tolerance_j_kg = h_tolerance;
        bath.maximum_thermal_residual_j = energy_tolerance;
        let endpoint = bath
            .prepare_endpoint(input(&curve, 1e16, 1.0, 0.0))
            .unwrap();
        assert!(
            endpoint
                .solve_with_internal_power(0.0, |_, _| Ok(2.0), || Ok(()))
                .is_err()
        );
    }
}

#[test]
fn source_and_feedback_nonfinite_values_are_refused() {
    let curve = sensible_curve(0.0);
    let endpoint = environment()
        .prepare_endpoint(input(&curve, 0.0, 1.0, 0.0))
        .unwrap();
    for bound in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            endpoint
                .solve_with_internal_power(
                    bound,
                    |_, _| panic!("invalid bound must be rejected first"),
                    || Ok(())
                )
                .is_err()
        );
    }
    assert!(
        endpoint
            .solve_with_internal_power(0.0, |_, _| Ok(f64::NAN), || Ok(()))
            .is_err()
    );
    assert!(environment()
        .prepare_endpoint(input(&curve, 0.0, f64::MAX, f64::MIN_POSITIVE))
        .is_err());
}
