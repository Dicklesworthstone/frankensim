use super::*;
use fs_couple::friction_heat::plastic_wear::adaptive::*;

fn velocity(speed: f64) -> Vec<[f64; 2]> {
    let mut values = vec![[0.0; 2]; 25];
    values[12][0] = speed;
    values
}

fn march_controls() -> PlasticSlidingMarchControls {
    PlasticSlidingMarchControls {
        initial_step_s: 0.1,
        minimum_step_s: 0.0001,
        maximum_step_s: 0.1,
        maximum_attempts: 256,
        maximum_accepted_substeps: 256,
        tolerance: PlasticSlidingMarchTolerance {
            cell_energy_j: 1.0,
            interface_temperature_k: 100.0,
            pressure_pa: 1e8,
            length_m: 1.0,
            work_j: 1.0,
            relative: 0.0,
        },
        energy_tolerance_j: 1e-9,
    }
}

#[test]
fn accepted_half_steps_are_exact_native_states_with_original_payloads_and_work() {
    with_cx(|cx, _| {
        let data = fixture(cx);
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind(&data, &thermal, &route, loaded(&data, cx), limits(), cx)
            .with_plastic_work_heat(PlasticWorkHeatPartition::new([0.3, 0.2]).unwrap());
        let saved = contact.state().clone();
        let saved_cells = thermal.cells().clone();
        let result = contact
            .march_prescribed_velocity_candidate(
                &thermal,
                0.1,
                &velocity(0.05),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                march_controls(),
                T0,
                transport,
                cx,
            )
            .unwrap();
        let first = contact
            .step_candidate(
                &thermal,
                0.05,
                &motion(0.05 * 0.05),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx,
            )
            .unwrap();
        let second = first
            .contact()
            .step_candidate(
                first.thermal(),
                0.1,
                &motion(0.05 * 0.05),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                first.report().heat.interface_temperature_k,
                transport,
                cx,
            )
            .unwrap();
        let report = result.report();
        assert_eq!(report.attempts, 1);
        assert_eq!(report.rejected_attempts, 0);
        assert_eq!(report.steps.len(), 2);
        assert_eq!(result.contact().state(), second.contact().state());
        assert_eq!(result.thermal().cells(), second.thermal().cells());
        assert_eq!(report.steps[0].report.source, first.report().source);
        assert_eq!(report.steps[1].report.source, second.report().source);
        assert!(report.thermal_trials > first.report().trials + second.report().trials);
        close(
            report.frictional_work_j,
            first.report().frictional_work_j + second.report().frictional_work_j,
            0.0,
            1e-18,
        );
        assert!(report.normal_plastic_work_j > 0.0);
        close(
            report.normal_plastic_heat_j[0],
            0.3 * report.normal_plastic_work_j,
            1e-14,
            1e-18,
        );
        close(
            report.normal_plastic_heat_j[1],
            0.2 * report.normal_plastic_work_j,
            1e-14,
            1e-18,
        );
        close(
            report.nonthermal_plastic_work_j,
            0.5 * report.normal_plastic_work_j,
            1e-14,
            1e-18,
        );
        close(
            report.thermal_energy_change_j,
            report.frictional_work_j
                + report.normal_plastic_heat_j.iter().sum::<f64>()
                + report.other_heat_j,
            0.0,
            1e-9,
        );
        close(
            report.contact_elastic_energy_change_j
                + report.thermal_energy_change_j
                + report.nonthermal_plastic_work_j,
            report.normal_work_j
                + report.temperature_parameter_work_j
                + report.wear_parameter_work_j
                + report.frictional_work_j
                + report.other_heat_j,
            0.0,
            1e-9,
        );
        close(
            result.contact().state().tangential_states()[12].slip_m()[0],
            0.005,
            0.0,
            1e-16,
        );
        assert_eq!(result.contact().time_s(), 0.1);
        assert_eq!(result.thermal().time_s(), 0.1);
        assert_eq!(contact.state(), &saved);
        assert_eq!(thermal.cells(), &saved_cells);
    });
}

// Independent continuous scalar reduction for the symmetric synthetic fixture:
// 2*g*(Ti-Tb)=A*H(Tb)*mu(Ti)*v, 2*C*dTb/dt=that friction power.
// H controls normal saturation, but no plastic work is converted into heat here.
// RK4 is used only for this independent smooth reference, not product state.
fn continuous_temperature(speed: f64, duration: f64) -> f64 {
    let capacity = RHO * AREA * LENGTH * CP;
    let g = CONDUCTIVITY * AREA / (0.5 * LENGTH);
    let rate = |theta: f64| {
        let force_velocity = AREA * 1e6 * (1.0 - 0.01 * theta) * speed;
        let interface =
            (2.0 * g * theta + 0.4 * force_velocity) / (2.0 * g - 0.01 * force_velocity);
        g * (interface - theta) / capacity
    };
    let dt = duration / 4096.0;
    let mut theta = 0.0;
    for _ in 0..4096 {
        let a = rate(theta);
        let b = rate(theta + 0.5 * dt * a);
        let c = rate(theta + 0.5 * dt * b);
        let d = rate(theta + dt * c);
        theta += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    T0 + theta
}

#[test]
fn physical_error_refinement_converges_to_independent_caloric_hardness_reference() {
    with_cx(|cx, _| {
        let data = fixture(cx);
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind(&data, &thermal, &route, loaded(&data, cx), limits(), cx);
        let run = |settings| {
            contact
                .march_prescribed_velocity_candidate(
                    &thermal,
                    0.1,
                    &velocity(0.05),
                    PlasticSlidingHeatBoundary::HeldApproach,
                    controls(),
                    settings,
                    T0,
                    transport,
                    cx,
                )
                .unwrap()
        };
        let coarse = run(march_controls());
        let mut tight = march_controls();
        tight.tolerance = PlasticSlidingMarchTolerance {
            cell_energy_j: 1e-9,
            interface_temperature_k: 0.003,
            pressure_pa: 30.0,
            length_m: 5e-10,
            work_j: 1e-9,
            relative: 1e-3,
        };
        let refined = run(tight);
        let reference = continuous_temperature(0.05, 0.1);
        let coarse_error = (coarse.thermal().temperature_k(CELLS[0]).unwrap() - reference).abs();
        let refined_error = (refined.thermal().temperature_k(CELLS[0]).unwrap() - reference).abs();
        assert!(refined.report().rejected_attempts > 0);
        assert!(refined.report().steps.len() > coarse.report().steps.len());
        assert!(refined.report().maximum_accepted_normalized_defect <= 1.0);
        assert!(
            refined_error < 0.5 * coarse_error,
            "coarse={coarse_error}, refined={refined_error}"
        );
        close(
            refined.thermal().temperature_k(CELLS[0]).unwrap(),
            reference,
            0.0,
            0.015,
        );
        // Body A's separately sourced hardness equals the local saturated
        // pressure, so its native Archard integral is exactly k*A*v*t.
        close(
            refined.contact().state().wear_state().volumes_m3()[0],
            1e-6 * AREA * 0.05 * 0.1,
            2e-6,
            1e-24,
        );
        close(
            refined.contact().state().tangential_states()[12].slip_m()[0],
            0.005,
            0.0,
            1e-15,
        );
        assert!(refined.report().combined_energy_residual_j.abs() <= tight.energy_tolerance_j);
    });
}

#[test]
fn late_budget_source_and_cancellation_refusals_publish_no_partial_span() {
    with_cx(|cx, gate| {
        let data = fixture(cx);
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind(&data, &thermal, &route, loaded(&data, cx), limits(), cx);
        let saved = contact.state().clone();
        let cells = thermal.cells().clone();
        let mut settings = march_controls();
        settings.initial_step_s = 0.025;
        settings.maximum_step_s = 0.025;
        settings.maximum_attempts = 1;
        let result = contact.march_prescribed_velocity_candidate(
            &thermal,
            0.1,
            &velocity(0.05),
            PlasticSlidingHeatBoundary::HeldApproach,
            controls(),
            settings,
            T0,
            transport,
            cx,
        );
        assert!(matches!(
            result,
            Err(PlasticSlidingMarchError::AttemptBudget)
        ));
        settings.maximum_attempts = 16;
        let seen_time = Cell::new(0.0_f64);
        let result = contact.march_prescribed_velocity_candidate(
            &thermal,
            0.1,
            &velocity(0.05),
            PlasticSlidingHeatBoundary::HeldApproach,
            controls(),
            settings,
            T0,
            |trial, temperature| {
                seen_time.set(seen_time.get().max(trial.time_s()));
                if trial.time_s() > 0.04 {
                    Err(DepositionError::InvalidInput(
                        "synthetic missing transport domain",
                    ))
                } else {
                    transport(trial, temperature)
                }
            },
            cx,
        );
        assert!(matches!(
            result,
            Err(PlasticSlidingMarchError::Step(
                PlasticSlidingHeatError::Source(_)
            ))
        ));
        assert!(seen_time.get() > 0.04);
        assert_eq!(contact.state(), &saved);
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(contact.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        gate.request();
        assert!(
            contact
                .march_prescribed_velocity_candidate(
                    &thermal,
                    0.1,
                    &velocity(0.05),
                    PlasticSlidingHeatBoundary::HeldApproach,
                    controls(),
                    settings,
                    T0,
                    transport,
                    cx,
                )
                .is_err()
        );
        assert_eq!(contact.state(), &saved);
        assert_eq!(thermal.cells(), &cells);
    });
}

#[test]
fn whole_span_energy_gate_is_independent_of_lax_native_solve_tolerances() {
    with_cx(|cx, _| {
        let data = fixture(cx);
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind(&data, &thermal, &route, loaded(&data, cx), limits(), cx);
        let saved = contact.state().clone();
        let mut native = controls();
        native.newton.absolute_tolerance = 1.0;
        native.newton.relative_tolerance = 0.0;
        native.heat = FrictionHeatTolerances {
            cell_j: 1.0,
            thermal_j: 1.0,
            interface_j: 1.0,
            total_j: 1.0,
        };
        let mut settings = march_controls();
        settings.energy_tolerance_j = 1e-12;
        let result = contact.march_prescribed_velocity_candidate(
            &thermal,
            0.1,
            &velocity(0.05),
            PlasticSlidingHeatBoundary::HeldApproach,
            native,
            settings,
            T0,
            transport,
            cx,
        );
        assert!(matches!(
            result,
            Err(PlasticSlidingMarchError::Energy { .. })
        ));
        assert_eq!(contact.state(), &saved);
        assert_eq!(contact.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
