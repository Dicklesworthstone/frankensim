//! Whole-owner energy/entropy checks on synthetic source-backed material laws.

use super::*;
use crate::cavity_phs::piston::tests::{controls, drive, gas};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_ports::PortedModule;
use std::convert::Infallible;

fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-100.0, 100.0],
        source_tolerance_j: 1e-10,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
    }
}
fn sensible() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        fs_blake3::hash_bytes(b"synthetic entropy contact cp=500"),
        SolidLiquidPhase::Solid,
        [200.0, 1000.0]
            .into_iter()
            .map(|temperature_k| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 500.0 * (temperature_k - 200.0),
                temperature_k,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 4000.0,
            })
            .collect(),
    )
    .unwrap()
}
fn latent(transition: f64) -> EquilibriumEnthalpyPhaseCurve {
    let knot = |h, t, f| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: 4000.0,
    };
    EquilibriumEnthalpyPhaseCurve::try_new(
        fs_blake3::hash_bytes(format!("synthetic latent contact T={transition}").as_bytes()),
        vec![
            knot(0.0, 200.0, 0.0),
            knot(100000.0, transition, 0.0),
            knot(160000.0, transition, 1.0),
            knot(400000.0, 900.0, 1.0),
        ],
    )
    .unwrap()
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: 0.001,
        volume_m3: 0.001 / 4000.0,
        surface_area_m2: 0.02,
        internal_heat_j: 0.0,
        duration_s: drive(0.0).duration_s,
    }
}
fn link(
    g: f64,
) -> impl FnMut(
    &ResolvedThermochemicalGasState,
    EquilibriumPhaseState,
) -> Result<ConductanceSample<(EquilibriumPhaseState, f64)>, Infallible> {
    move |gas, body| {
        Ok(ConductanceSample {
            conductance_w_k: g,
            source: (body, gas.state().temperature),
        })
    }
}
fn fixed(slope: f64) -> GasPiston {
    GasPiston::new(gas(0.75, slope), 1.0, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap()
}
fn logarithmic_mean(a: f64, b: f64) -> f64 {
    if a == b {
        a
    } else {
        (b - a) / ((b - a) / a).ln_1p()
    }
}
fn check_account<S>(candidate: &DiscreteGradientThermalCandidate<S>) {
    let step = candidate.step();
    assert!(step.source_residual_j.abs() <= exchange().source_tolerance_j);
    assert!(step.thermal_energy_residual_j.abs() <= exchange().energy_tolerance_j);
    assert!(step.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
    assert!(step.contact_entropy_production_j_k >= 0.0);
    assert!(step.entropy_residual_j_k.abs() < 5e-10);
    assert!(step.entropy_tolerance_j_k.is_finite() && step.entropy_tolerance_j_k > 0.0);
    assert!(step.entropy_residual_j_k.abs() <= step.entropy_tolerance_j_k);
    assert_eq!(candidate.thermal().external_heat_j, -step.exchanged_heat_j);
    assert!(candidate.piston().audit_last_window().is_ok());
}

#[test]
fn fixed_sensible_contact_matches_original_entropy_and_logarithmic_means_in_both_directions() {
    let old = fixed(0.0);
    let curve = sensible();
    let cg = old.gas().mass_kg()
        * old
            .gas()
            .state()
            .thermodynamics()
            .properties()
            .mass_specific_cv()
            .value();
    for initial_body_temperature in [300.0, 650.0] {
        let thermal = body(&curve, 500.0 * (initial_body_temperature - 200.0));
        let candidate = old
            .thermal_discrete_gradient_candidate(
                drive(0.0),
                thermal,
                controls(),
                exchange(),
                link(4000.0),
                || false,
            )
            .unwrap();
        check_account(&candidate);
        let report = candidate.step();
        let final_gas_temperature = candidate.piston().gas().state().state().temperature;
        let final_body_temperature = candidate.thermal().state.temperature_k();
        let tg = logarithmic_mean(400.0, final_gas_temperature);
        let tb = logarithmic_mean(initial_body_temperature, final_body_temperature);
        assert!((report.piston.thermodynamic_mean.state().temperature - tg).abs() < 1e-9);
        assert!((report.thermal_mean_temperature_k - tb).abs() < 1e-9);
        assert!(
            (report.thermal_mean_temperature_k - report.thermal_midpoint.temperature_k()).abs()
                > 0.01
        );
        let entropy = cg * (final_gas_temperature / 400.0).ln()
            + 0.5 * (final_body_temperature / initial_body_temperature).ln();
        assert!((report.total_entropy_change_j_k - entropy).abs() < 1e-12);
        assert!((entropy - report.contact_entropy_production_j_k).abs() < 5e-10);
        assert!(report.total_entropy_change_j_k > 0.0);
        assert_eq!(report.prescribed_source_entropy_j_k, 0.0);
        assert_eq!(candidate.thermal().report.0, report.thermal_midpoint);
        assert_eq!(
            candidate.thermal().report.1,
            report.piston.thermodynamic_mean.state().temperature
        );
        assert_eq!(candidate.piston().gas().volume_m3(), old.gas().volume_m3());
        assert_eq!(report.piston.support_work_j, 0.0);
        assert!(report.exchanged_heat_j * (initial_body_temperature - 400.0) > 0.0);
    }
}

#[test]
fn moving_variable_capacity_gas_and_finite_body_close_energy_and_accumulated_entropy() {
    let curve = sensible();
    for velocity in [-2.0, 2.0] {
        let mut piston = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, velocity).unwrap();
        let initial = piston.clone();
        let mut thermal = body(&curve, 225000.0);
        let h0 = thermal.initial.specific_enthalpy_j_kg();
        let mut produced = 0.0;
        let mut mechanical_input = 0.0;
        for _ in 0..12 {
            let candidate = piston
                .thermal_discrete_gradient_candidate(
                    drive(0.0),
                    thermal,
                    controls(),
                    exchange(),
                    link(2000.0),
                    || false,
                )
                .unwrap();
            check_account(&candidate);
            let (next, accepted, report) = candidate.into_parts();
            produced += report.contact_entropy_production_j_k;
            mechanical_input += report.piston.external_force_work_j
                - report.piston.exterior_pressure_work_j
                + report.piston.support_work_j;
            assert!(report.piston.displacement_m != 0.0);
            assert_eq!(accepted.state.phase_curve_identity(), curve.identity());
            thermal.initial = accepted.state;
            piston = next;
        }
        let gas_entropy = initial.gas().mass_kg()
            * initial
                .gas()
                .model()
                .specific_entropy_change_j_kg_k(initial.gas().point(), piston.gas().point())
                .unwrap();
        let body_entropy = thermal.mass_kg
            * curve
                .specific_entropy_change_j_kg_k(h0, thermal.initial.specific_enthalpy_j_kg())
                .unwrap();
        assert!((gas_entropy + body_entropy - produced).abs() < 2e-9);
        assert!(gas_entropy + body_entropy > 0.0);
        let gas_energy = initial.gas().mass_kg()
            * initial
                .gas()
                .model()
                .thermodynamics()
                .specific_internal_energy_change(
                    fs_qty::Temperature::new(initial.gas().state().state().temperature),
                    fs_qty::Temperature::new(piston.gas().state().state().temperature),
                )
                .unwrap()
                .value();
        let kinetic = 0.5
            * piston.mass_kg()
            * (piston.velocity_m_s().powi(2) - initial.velocity_m_s().powi(2));
        let body_energy = thermal.mass_kg * (thermal.initial.specific_enthalpy_j_kg() - h0);
        assert!(
            (gas_energy + kinetic + body_energy - mechanical_input).abs()
                < 12.0 * exchange().energy_tolerance_j
        );
        assert_eq!(piston.gas().mass_kg(), initial.gas().mass_kg());
        assert_eq!(piston.gas().point().axes().get("specimen-lot"), Some(&42.0));
    }
}

#[test]
fn melting_and_freezing_keep_latent_entropy_and_phase_dependent_conductance() {
    let old = fixed(0.001);
    for transition in [350.0, 500.0] {
        let curve = latent(transition);
        let thermal = body(&curve, 130000.0);
        let candidate = old
            .thermal_discrete_gradient_candidate(
                drive(0.0),
                thermal,
                controls(),
                exchange(),
                |gas, phase| {
                    Ok::<_, Infallible>(ConductanceSample {
                        conductance_w_k: 2000.0 * (1.0 + phase.liquid_mass_fraction()),
                        source: (phase, gas.state().temperature),
                    })
                },
                || false,
            )
            .unwrap();
        check_account(&candidate);
        let report = candidate.step();
        assert_eq!(candidate.thermal().state.temperature_k(), transition);
        assert!((report.thermal_mean_temperature_k - transition).abs() < 1e-12);
        assert!(
            (report.thermal_entropy_change_j_k - report.thermal_energy_change_j / transition).abs()
                < 1e-13
        );
        assert!(
            (candidate.thermal().state.liquid_mass_fraction()
                - thermal.initial.liquid_mass_fraction())
                * (400.0 - transition)
                > 0.0
        );
        assert_eq!(
            report.conductance_w_k,
            2000.0 * (1.0 + report.thermal_midpoint.liquid_mass_fraction())
        );
        assert_eq!(candidate.thermal().report.0, report.thermal_midpoint);
    }
}

#[test]
fn prescribed_sources_cross_the_latent_interval_without_claiming_reservoir_entropy() {
    let old = fixed(0.001);
    let curve = latent(350.0);
    let thermal = UniformEnthalpyStepInput {
        internal_heat_j: 100.0,
        ..body(&curve, 90000.0)
    };
    let candidate = old
        .thermal_discrete_gradient_candidate(
            drive(20.0),
            thermal,
            controls(),
            exchange(),
            link(4000.0),
            || false,
        )
        .unwrap();
    check_account(&candidate);
    let report = candidate.step();
    let h1 = candidate.thermal().state.specific_enthalpy_j_kg();
    assert!(h1 > 160000.0);
    let cp_solid = 100000.0 / 150.0;
    let cp_liquid = 240000.0 / 550.0;
    let expected_body_entropy = thermal.mass_kg
        * (cp_solid * (350.0 / thermal.initial.temperature_k()).ln()
            + 60000.0 / 350.0
            + cp_liquid * (candidate.thermal().state.temperature_k() / 350.0).ln());
    assert!((report.thermal_entropy_change_j_k - expected_body_entropy).abs() < 1e-12);
    let expected_source_entropy = 20.0 / report.piston.thermodynamic_mean.state().temperature
        + 100.0 / report.thermal_mean_temperature_k;
    assert_eq!(
        report.prescribed_source_entropy_j_k,
        expected_source_entropy
    );
    assert!((report.thermal_energy_change_j + report.exchanged_heat_j - 100.0).abs() < 1e-7);
    assert!((report.piston.gas.supplied_heat_j - report.exchanged_heat_j - 20.0).abs() < 1e-12);
    // The independently heated body can physically overtake the gas. The closed
    // pair's equilibrium-order admission must not reject this driven case.
    assert!(thermal.initial.temperature_k() < old.gas().state().state().temperature);
    assert!(
        candidate.thermal().state.temperature_k()
            > candidate.piston().gas().state().state().temperature
    );

    let curve = sensible();
    let cooling = UniformEnthalpyStepInput {
        internal_heat_j: -10.0,
        ..body(&curve, 225000.0)
    };
    let candidate = old
        .thermal_discrete_gradient_candidate(
            drive(-5.0),
            cooling,
            controls(),
            exchange(),
            link(1.0),
            || false,
        )
        .unwrap();
    check_account(&candidate);
    assert!(candidate.step().total_entropy_change_j_k < 0.0);
    assert!(candidate.step().prescribed_source_entropy_j_k < 0.0);
    assert!(candidate.step().contact_entropy_production_j_k > 0.0);

    // Externally maintained contact: the body source replaces exactly the
    // exported heat while the gas source removes it. Conduction still produces
    // entropy even though the material states can remain stationary.
    let maintained = UniformEnthalpyStepInput {
        internal_heat_j: 10.0,
        ..body(&curve, 150000.0)
    };
    let candidate = old
        .thermal_discrete_gradient_candidate(
            drive(-10.0),
            maintained,
            controls(),
            exchange(),
            link(2000.0),
            || false,
        )
        .unwrap();
    check_account(&candidate);
    assert!((candidate.step().exchanged_heat_j - 10.0).abs() < 1e-8);
    assert!((candidate.piston().gas().state().state().temperature - 400.0).abs() < 1e-8);
    assert!((candidate.thermal().state.temperature_k() - 500.0).abs() < 1e-8);
    assert!((candidate.step().contact_entropy_production_j_k - 0.005).abs() < 1e-10);
    assert!((candidate.step().prescribed_source_entropy_j_k + 0.005).abs() < 1e-10);
}

#[test]
fn equilibrium_is_stationary_and_insulation_preserves_independent_source_updates() {
    let old = fixed(0.001);
    let curve = sensible();
    let candidate = old
        .thermal_discrete_gradient_candidate(
            drive(0.0),
            body(&curve, 100000.0),
            controls(),
            exchange(),
            link(4000.0),
            || false,
        )
        .unwrap();
    assert_eq!(candidate.piston(), &old);
    assert_eq!(candidate.step().exchanged_heat_j, 0.0);
    assert_eq!(candidate.step().total_entropy_change_j_k, 0.0);
    assert_eq!(candidate.step().exchange_trials, 1);

    let thermal = UniformEnthalpyStepInput {
        internal_heat_j: 25.0,
        ..body(&curve, 150000.0)
    };
    let candidate = old
        .thermal_discrete_gradient_candidate(
            drive(20.0),
            thermal,
            controls(),
            exchange(),
            link(0.0),
            || false,
        )
        .unwrap();
    check_account(&candidate);
    let mut reference = old.clone();
    let gas_step = reference
        .step_discrete_gradient(drive(20.0), controls(), || false)
        .unwrap();
    assert_eq!(candidate.piston(), &reference);
    assert_eq!(candidate.step().piston, gas_step);
    assert_eq!(
        candidate.thermal().state,
        curve.state_at_specific_enthalpy(175000.0).unwrap()
    );
    assert_eq!(candidate.step().exchanged_heat_j, 0.0);
    assert_eq!(candidate.step().contact_entropy_production_j_k, 0.0);
}

#[test]
fn stiff_closed_fixed_contact_requests_time_refinement_instead_of_crossing_equilibrium() {
    let old = fixed(0.0);
    let curve = sensible();
    let thermal = body(&curve, 225000.0);
    let result = old.thermal_discrete_gradient_candidate(
        drive(0.0),
        thermal,
        controls(),
        exchange(),
        link(30000.0),
        || false,
    );
    assert!(
        matches!(result, Err(ThermalPistonError::Invalid(message)) if message.contains("equilibrium"))
    );
    let input = PistonInput {
        duration_s: 0.1 * drive(0.0).duration_s,
        ..drive(0.0)
    };
    let candidate = old
        .thermal_discrete_gradient_candidate(
            input,
            UniformEnthalpyStepInput {
                duration_s: input.duration_s,
                ..thermal
            },
            controls(),
            exchange(),
            link(30000.0),
            || false,
        )
        .unwrap();
    check_account(&candidate);
    assert!(
        candidate.thermal().state.temperature_k()
            > candidate.piston().gas().state().state().temperature
    );
    assert_eq!(old, fixed(0.0));
    assert_eq!(
        thermal.initial,
        curve.state_at_specific_enthalpy(225000.0).unwrap()
    );
}

#[test]
fn source_refusal_and_late_cancellation_preserve_both_physical_owners_and_accepted_ports() {
    let mut old = fixed(0.001);
    old.step_discrete_gradient(drive(0.0), controls(), || false)
        .unwrap();
    let saved = old.clone();
    let transfers = old.last_window_transfers().map(<[_]>::to_vec);
    let duration = old.last_window_duration_s();
    let curve = sensible();
    let thermal = body(&curve, 225000.0);
    let failure = |_: &ResolvedThermochemicalGasState,
                   _: EquilibriumPhaseState|
     -> Result<ConductanceSample<()>, &str> { Err("missing contact source") };
    assert!(matches!(
        old.thermal_discrete_gradient_candidate(
            drive(0.0),
            thermal,
            controls(),
            exchange(),
            failure,
            || false,
        ),
        Err(ThermalPistonError::Interface("missing contact source"))
    ));
    assert!(
        old.thermal_discrete_gradient_candidate(
            drive(0.0),
            thermal,
            controls(),
            exchange(),
            link(-1.0),
            || false,
        )
        .is_err()
    );
    assert!(matches!(
        old.thermal_discrete_gradient_candidate(
            drive(0.0),
            thermal,
            controls(),
            HeatExchangeSolve {
                heat_bracket_j: [-1e-6, 1e-6],
                ..exchange()
            },
            link(4000.0),
            || false,
        ),
        Err(ThermalPistonError::OutsideBracket { .. })
    ));
    let other = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        fs_blake3::hash_bytes(b"another physical source"),
        SolidLiquidPhase::Solid,
        curve.knots().to_vec(),
    )
    .unwrap();
    assert!(
        old.thermal_discrete_gradient_candidate(
            drive(0.0),
            UniformEnthalpyStepInput {
                curve: &other,
                ..thermal
            },
            controls(),
            exchange(),
            link(4000.0),
            || false,
        )
        .is_err()
    );
    let mut calls = 0;
    old.thermal_discrete_gradient_candidate(
        drive(0.0),
        thermal,
        controls(),
        exchange(),
        link(4000.0),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        old.thermal_discrete_gradient_candidate(
            drive(0.0),
            thermal,
            controls(),
            exchange(),
            link(4000.0),
            || {
                seen += 1;
                seen == calls
            },
        ),
        Err(ThermalPistonError::Cancelled)
    ));
    assert!(matches!(
        old.thermal_discrete_gradient_candidate(
            drive(0.0),
            thermal,
            controls(),
            exchange(),
            link(4000.0),
            || true,
        ),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(old, saved);
    assert_eq!(old.last_window_duration_s(), duration);
    assert_eq!(old.last_window_transfers().map(<[_]>::to_vec), transfers);
    assert_eq!(
        thermal.initial,
        curve.state_at_specific_enthalpy(225000.0).unwrap()
    );
}

mod real_gas;
