use super::*;
use fs_couple::friction_heat::plastic_wear::adaptive::{
    PlasticSlidingMarchControls, PlasticSlidingMarchError, PlasticSlidingMarchTolerance,
};
use fs_couple::friction_heat::plastic_wear::depth::adaptive::*;

fn velocity(speed: f64) -> Vec<[f64; 2]> {
    let mut values = vec![[0.0; 2]; 25];
    values[12][0] = speed;
    values
}

fn march_controls(duration: f64) -> DepthSlidingMarchControls {
    DepthSlidingMarchControls {
        march: PlasticSlidingMarchControls {
            initial_step_s: duration,
            minimum_step_s: 1e-5,
            maximum_step_s: duration,
            maximum_attempts: 512,
            maximum_accepted_substeps: 1024,
            tolerance: PlasticSlidingMarchTolerance {
                cell_energy_j: 1.0,
                interface_temperature_k: 100.0,
                pressure_pa: 1e8,
                length_m: 1.0,
                work_j: 1.0,
                relative: 0.0,
            },
            energy_tolerance_j: 1e-9,
        },
        body_temperature_k: 100.0,
        hardness_pa: 1e8,
    }
}

#[test]
fn g1_depth_march_retains_native_fine_history_payloads_and_all_three_work_accounts() {
    with_cx(|cx, _| {
        let data = depth_fixture(320.0, cx);
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind_depth(
            &data,
            &thermal,
            &route,
            depth_loaded(&data, cx),
            limits(),
            cx,
        )
        .with_plastic_work_heat(PlasticWorkHeatPartition::new([0.3, 0.2]).unwrap());
        let saved = contact.source().clone();
        let saved_cells = thermal.cells().clone();
        let duration = 0.13;
        let speed = 0.05;
        let result = contact
            .march_prescribed_velocity_candidate(
                &thermal,
                duration,
                &velocity(speed),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                march_controls(duration),
                T0,
                transport,
                cx,
            )
            .unwrap();
        let first = contact
            .step_candidate(
                &thermal,
                0.5 * duration,
                &motion(speed * (0.5 * duration)),
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
                duration,
                &motion(speed * (duration - first.contact().time_s())),
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
        assert_eq!(result.contact().source(), second.contact().source());
        assert_eq!(result.thermal().cells(), second.thermal().cells());
        assert_eq!(report.steps[0].report.source, first.report().source);
        assert_eq!(report.steps[1].report.source, second.report().source);
        assert_eq!(report.steps[0].report.mechanics, first.report().mechanics);
        assert_eq!(report.steps[1].report.mechanics, second.report().mechanics);
        assert!(report.thermal_trials > first.report().trials + second.report().trials);
        let mut elastic = 0.0;
        let mut plastic = 0.0;
        let mut normal_work = 0.0;
        let mut cap_parameter = 0.0;
        for native in [&first, &second] {
            let mechanics = &native.report().mechanics;
            assert_eq!(mechanics.wear_normals().len(), 2);
            for r in [
                mechanics.temperature_normal(),
                &mechanics.wear_normals()[0].normal,
                &mechanics.wear_normals()[1].normal,
            ] {
                elastic += r.elastic_energy_change_j;
                plastic += r.response.plastic_work_j;
                normal_work += r.normal_work_j;
            }
            cap_parameter += mechanics.wear_normals()[1].normal.parameter_work_j;
        }
        assert!(cap_parameter.abs() > 0.0);
        close(
            report.contact_elastic_energy_change_j,
            elastic,
            1e-13,
            1e-18,
        );
        close(report.normal_plastic_work_j, plastic, 1e-13, 1e-18);
        close(report.normal_work_j, normal_work, 1e-13, 1e-18);
        close(report.normal_plastic_heat_j[0], 0.3 * plastic, 1e-13, 1e-18);
        close(report.normal_plastic_heat_j[1], 0.2 * plastic, 1e-13, 1e-18);
        close(
            report.nonthermal_plastic_work_j,
            0.5 * plastic,
            1e-13,
            1e-18,
        );
        close(
            report.frictional_work_j,
            first.report().frictional_work_j + second.report().frictional_work_j,
            0.0,
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
            result.contact().source().tangential_states()[12].slip_m()[0],
            speed * duration,
            0.0,
            1e-16,
        );
        assert_eq!(result.contact().time_s(), duration);
        assert_eq!(result.thermal().time_s(), duration);
        assert_eq!(contact.source(), &saved);
        assert_eq!(thermal.cells(), &saved_cells);
    });
}

// Independent continuous reduction, away from any active-set change. The
// softer graded body stays saturated: p=H_A(T,z_A). Each finite receiver obeys
// C*dT/dt=g*(Ti-T), the massless interface balances mu(Ti)*p*A*v, and each
// body's dz/dt uses its OWN H(T,z). This reference advances only three scalars;
// it is not used to create a product state or to emulate the native split.
fn continuous_reference(speed: f64, duration: f64) -> [f64; 3] {
    let capacity = RHO * AREA * LENGTH * CP;
    let g = CONDUCTIVITY * AREA / (0.5 * LENGTH);
    let rate = |[theta, za, zb]: [f64; 3]| {
        let ha = 1e6 * (1.0 - 0.01 * theta) + GRADIENT[0] * za;
        let hb = 2e6 + GRADIENT[1] * zb;
        let force_velocity = AREA * ha * speed;
        let interface =
            (2.0 * g * theta + 0.4 * force_velocity) / (2.0 * g - 0.01 * force_velocity);
        let coefficient = 1.0 + 0.02 * interface;
        [
            g * (interface - theta) / capacity,
            WEAR_K[0] * coefficient * speed,
            WEAR_K[1] * coefficient * ha * speed / hb,
        ]
    };
    let dt = duration / 4096.0;
    let mut y = [0.0; 3];
    for _ in 0..4096 {
        let a = rate(y);
        let b = rate(std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = rate(std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = rate(std::array::from_fn(|i| y[i] + dt * c[i]));
        y = std::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    [T0 + y[0], y[1], y[2]]
}

#[test]
fn g1_depth_refinement_converges_to_independent_live_hardness_heat_and_two_body_wear() {
    with_cx(|cx, _| {
        let data = depth_fixture(320.0, cx);
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind_depth(
            &data,
            &thermal,
            &route,
            depth_loaded(&data, cx),
            limits(),
            cx,
        );
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
        let coarse = run(march_controls(0.1));
        let mut tight = march_controls(0.1);
        tight.march.tolerance = PlasticSlidingMarchTolerance {
            cell_energy_j: 1e-9,
            interface_temperature_k: 0.003,
            pressure_pa: 30.0,
            length_m: 5e-11,
            work_j: 1e-9,
            relative: 1e-3,
        };
        tight.body_temperature_k = 0.003;
        tight.hardness_pa = 30.0;
        let fine = run(tight);
        assert!(fine.report().rejected_attempts > 0);
        assert!(fine.report().steps.len() > coarse.report().steps.len());
        assert!(fine.report().maximum_accepted_normalized_defect <= 1.0);
        let reference = continuous_reference(0.05, 0.1);
        let coarse_values = [
            coarse.thermal().temperature_k(CELLS[0]).unwrap(),
            coarse.contact().source().wear_state().recession_m()[0][12],
            coarse.contact().source().wear_state().recession_m()[1][12],
        ];
        let fine_values = [
            fine.thermal().temperature_k(CELLS[0]).unwrap(),
            fine.contact().source().wear_state().recession_m()[0][12],
            fine.contact().source().wear_state().recession_m()[1][12],
        ];
        for i in 0..3 {
            let coarse_error = (coarse_values[i] - reference[i]).abs();
            let fine_error = (fine_values[i] - reference[i]).abs();
            assert!(
                fine_error < 0.5 * coarse_error,
                "component {i}: coarse={coarse_error:e}, fine={fine_error:e}"
            );
        }
        close(fine_values[0], reference[0], 0.0, 0.015);
        for side in 0..2 {
            close(fine_values[side + 1], reference[side + 1], 0.003, 1e-15);
            close(
                fine.contact().source().wear_state().volumes_m3()[side],
                AREA * fine_values[side + 1],
                1e-13,
                1e-27,
            );
        }
        let expected_h = 1e6 * (1.0 - 0.01 * (reference[0] - T0)) + GRADIENT[0] * reference[1];
        close(
            fine.contact()
                .source()
                .wear_state()
                .contact_state()
                .pressure_pa()[12],
            expected_h,
            5e-4,
            1.0,
        );
        assert!(fine.report().combined_energy_residual_j.abs() <= tight.march.energy_tolerance_j);
    });
}

#[test]
fn g2_body_temperature_and_original_hardness_have_independent_error_gates() {
    with_cx(|cx, _| {
        let data = depth_fixture(320.0, cx);
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind_depth(
            &data,
            &thermal,
            &route,
            depth_loaded(&data, cx),
            limits(),
            cx,
        );
        let coarse = contact
            .step_candidate(
                &thermal,
                0.1,
                &motion(0.005),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx,
            )
            .unwrap();
        let fine = contact
            .march_prescribed_velocity_candidate(
                &thermal,
                0.1,
                &velocity(0.05),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                march_controls(0.1),
                T0,
                transport,
                cx,
            )
            .unwrap();
        let temperature_error = CELLS
            .into_iter()
            .map(|id| {
                (coarse.thermal().temperature_k(id).unwrap()
                    - fine.thermal().temperature_k(id).unwrap())
                .abs()
            })
            .fold(0.0_f64, f64::max);
        let mut hardness_error = 0.0_f64;
        for side in 0..2 {
            for i in 0..25 {
                let a = coarse.contact().source().wear_state();
                let b = fine.contact().source().wear_state();
                hardness_error = hardness_error.max(
                    (a.hardness_profiles()[side]
                        .hardness_pa_at(a.recession_m()[side][i])
                        .unwrap()
                        - b.hardness_profiles()[side]
                            .hardness_pa_at(b.recession_m()[side][i])
                            .unwrap())
                    .abs(),
                );
                hardness_error = hardness_error.max(
                    (a.model().config().hardness_pa[i] - b.model().config().hardness_pa[i]).abs(),
                );
            }
        }
        assert!(temperature_error > 0.0);
        assert!(hardness_error > 0.0);
        for hardness_only in [false, true] {
            let mut settings = march_controls(0.1);
            settings.march.minimum_step_s = 0.1;
            if hardness_only {
                settings.hardness_pa = 0.5 * hardness_error;
            } else {
                settings.body_temperature_k = 0.5 * temperature_error;
            }
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
            let Err(PlasticSlidingMarchError::Unresolved { defect, .. }) = result else {
                panic!("a separate physical scale did not refuse the unresolved interval");
            };
            assert!(defect.common.normalized < 1.0);
            close(defect.body_temperature_k, temperature_error, 1e-8, 1e-10);
            close(defect.hardness_pa, hardness_error, 1e-8, 1e-5);
            close(defect.normalized, 2.0, 1e-8, 1e-8);
        }
    });
}

#[test]
fn g3_depth_march_late_budget_source_and_cancellation_refuse_the_complete_span() {
    with_cx(|cx, gate| {
        let data = depth_fixture(320.0, cx);
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind_depth(
            &data,
            &thermal,
            &route,
            depth_loaded(&data, cx),
            limits(),
            cx,
        );
        let saved = contact.source().clone();
        let cells = thermal.cells().clone();
        let mut settings = march_controls(0.025);
        settings.march.maximum_attempts = 1;
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
        settings.march.maximum_attempts = 16;
        settings.march.maximum_accepted_substeps = 2;
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
            Err(PlasticSlidingMarchError::AcceptedSubstepBudget)
        ));
        settings.march.maximum_accepted_substeps = 32;
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
        assert_eq!(contact.source(), &saved);
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
        assert_eq!(contact.source(), &saved);
        assert_eq!(thermal.cells(), &cells);
    });
}
