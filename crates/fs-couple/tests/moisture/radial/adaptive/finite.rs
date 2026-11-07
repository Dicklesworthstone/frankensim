//! Independent finite-gas decay and actual unequal radial inventories.

use super::*;
use crate::atmosphere::{
    atmosphere as gas, carrier_mass, controls as finite_controls, sorption as sorption_law,
};
use fs_couple::moisture::atmosphere::FiniteMoistureAtmosphere;

fn check_account<Step>(
    air: &FiniteMoistureAtmosphere,
    original: &MoistureRadial,
    next_air: &FiniteMoistureAtmosphere,
    next: &MoistureRadial,
    report: &MoistureAdvanceReport<Step, 2>,
) {
    let tolerance = 1e-10;
    near(
        next_air.vapor_mass_kg().unwrap() + next.water_mass_kg(),
        air.vapor_mass_kg().unwrap() + original.water_mass_kg(),
        tolerance,
    );
    near(
        next_air.gas().mass_kg() + next.wet_mass_kg(),
        air.gas().mass_kg() + original.wet_mass_kg(),
        tolerance,
    );
    near(carrier_mass(next_air), carrier_mass(air), 1e-12);
    near(report.stored_water_kg[0], original.water_mass_kg(), 1e-13);
    near(report.stored_water_kg[1], next.water_mass_kg(), 1e-13);
    assert_eq!(
        report.reservoir_water_kg.unwrap(),
        [
            air.vapor_mass_kg().unwrap(),
            next_air.vapor_mass_kg().unwrap()
        ]
    );
    assert_eq!(
        report.activity.unwrap(),
        [
            air.relative_humidity().unwrap(),
            next_air.relative_humidity().unwrap()
        ]
    );
    assert_eq!(report.boundary_water_kg[0], 0.0);
    near(
        report.boundary_water_kg[1],
        report.water_change_kg,
        tolerance,
    );
    near(report.transport_balance_residual_kg, 0.0, tolerance);
    near(report.water_balance_residual_kg, 0.0, tolerance);
    assert!(report.maximum_carrier_mass_residual_kg.unwrap().abs() <= tolerance);
    assert!(report.total_mass_residual_kg.unwrap().abs() <= tolerance);
    assert!(report.maximum_accepted_error_ratio <= 1.0);
    assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
}

fn unequal_shell_exact(capacity: f64, initial_rh: f64, time: f64) -> ([f64; 2], f64) {
    let (mass, _) = two_shell_parameters(RadialMoistureSymmetry::Sphere);
    let slope = 0.25;
    let conductance = 500.0 * 0.001 * 4.0 * core::f64::consts::PI / (1.0 / 0.0275 - 1.0 / 0.0675);
    let boundary = 500.0 * 0.001 * 4.0 * core::f64::consts::PI / (1.0 / 0.0675 - 1.0 / 0.1);
    let mean = slope * (capacity * initial_rh + mass[0] * 0.22 + mass[1] * 0.08)
        / (capacity + slope * (mass[0] + mass[1]));
    let a = conductance / mass[0];
    let b = conductance / mass[1];
    let c = boundary / mass[1];
    // Eliminate RH through C*RH + m0*M0 + m1*M1 = constant. The two
    // remaining physical diffusion modes have the exact matrix exponential.
    let matrix = [
        [-a, a],
        [
            b - c * slope * mass[0] / capacity,
            -b - c * (1.0 + slope * mass[1] / capacity),
        ],
    ];
    let trace = matrix[0][0] + matrix[1][1];
    let determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    let discriminant = (trace * trace - 4.0 * determinant).sqrt();
    let rates = [0.5 * (trace + discriminant), 0.5 * (trace - discriminant)];
    let deviation = [0.22 - mean, 0.08 - mean];
    let derivative: [f64; 2] =
        core::array::from_fn(|i| matrix[i][0] * deviation[0] + matrix[i][1] * deviation[1]);
    let exact: [f64; 2] = core::array::from_fn(|i| {
        mean + ((derivative[i] - rates[1] * deviation[i]) * (rates[0] * time).exp()
            + (rates[0] * deviation[i] - derivative[i]) * (rates[1] * time).exp())
            / (rates[0] - rates[1])
    });
    let humidity =
        initial_rh - (mass[0] * (exact[0] - 0.22) + mass[1] * (exact[1] - 0.08)) / capacity;
    (exact, humidity)
}

#[test]
fn finite_radial_adaptive_refines_two_unequal_shells_and_the_real_gas_inventory() {
    with_cx(|cx| {
        let original = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.02, 0.035, 0.1],
            vec![0.22, 0.08],
            [0.001; 2],
            cx,
        );
        let air = gas(0.8, 1.0);
        let original_air = air.clone();
        let rh = air.relative_humidity().unwrap();
        let capacity = air.vapor_mass_kg().unwrap() / rh;
        let (exact, exact_rh) = unequal_shell_exact(capacity, rh, 0.1);
        let faces = [None, Some(sorption_law(0.25))];
        for second_order in [false, true] {
            let mut errors = Vec::new();
            let mut humidity_errors = Vec::new();
            let mut work = Vec::new();
            for tolerance in [1e-3, 1e-5] {
                let time = controls(0.1, tolerance);
                let (next_air, next, evaluations) = if second_order {
                    let (next_air, next, report) = air
                        .advance_radial_sdirk2_adaptive_to(
                            &original,
                            0.1,
                            &faces,
                            config(),
                            nonlinear(),
                            finite_controls(),
                            time,
                            cx,
                        )
                        .unwrap();
                    check_account(&air, &original, &next_air, &next, &report);
                    assert!(report.rejected_intervals > 0);
                    (next_air, next, report.evaluations)
                } else {
                    let (next_air, next, report) = air
                        .advance_radial_adaptive_to(
                            &original,
                            0.1,
                            &faces,
                            config(),
                            nonlinear(),
                            finite_controls(),
                            time,
                            cx,
                        )
                        .unwrap();
                    check_account(&air, &original, &next_air, &next, &report);
                    assert!(report.rejected_intervals > 0);
                    (next_air, next, report.evaluations)
                };
                assert_eq!(next.time_s(), 0.1);
                errors.push(
                    next.moisture()
                        .iter()
                        .zip(exact)
                        .map(|(&actual, expected)| (actual - expected).abs())
                        .fold(0.0_f64, f64::max),
                );
                humidity_errors.push((next_air.relative_humidity().unwrap() - exact_rh).abs());
                work.push(evaluations);
            }
            assert!(errors[1] < 0.2 * errors[0], "{errors:?}");
            assert!(
                humidity_errors[1] < 0.2 * humidity_errors[0],
                "{humidity_errors:?}"
            );
            assert!(errors[1] < if second_order { 5e-7 } else { 3e-5 });
            assert!(humidity_errors[1] < if second_order { 3e-5 } else { 5e-4 });
            assert!(work[1] > work[0]);
        }
        assert_eq!(air, original_air);
        assert_eq!(original.moisture(), &[0.22, 0.08]);
        assert_eq!(original.time_s(), 0.0);
    });
}

fn one_shell_decay(z: f64, second_order: bool) -> f64 {
    if second_order {
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        (1.0 + (2.0 * gamma - 1.0) * z) / (1.0 + gamma * z).powi(2)
    } else {
        1.0 / (1.0 + z)
    }
}

#[test]
fn finite_radial_adaptive_observes_content_humidity_and_physical_water() {
    with_cx(|cx| {
        let radius: f64 = 0.05;
        let diffusion = 0.001;
        let slope = 0.25;
        let initial = 0.1;
        let original = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![initial],
            [diffusion; 2],
            cx,
        );
        let air = gas(0.8, 1.0);
        let rh = air.relative_humidity().unwrap();
        let capacity = air.vapor_mass_kg().unwrap() / rh;
        let dry = 500.0 * 4.0 * core::f64::consts::PI * radius.powi(3) / 3.0;
        let b = 1.0 + slope * dry / capacity;
        let equilibrium = (slope * rh + slope * dry * initial / capacity) / b;
        let lambda = 3.0 * diffusion / radius.powi(2) * b;
        let target = 0.4 / lambda;
        let faces = [None, Some(sorption_law(slope))];
        for second_order in [false, true] {
            let coarse = equilibrium
                + (initial - equilibrium) * one_shell_decay(lambda * target, second_order);
            let fine = equilibrium
                + (initial - equilibrium)
                    * one_shell_decay(lambda * target / 2.0, second_order).powi(2);
            let difference = (coarse - fine).abs();
            for (coordinate, difference) in
                [difference, dry * difference / capacity, dry * difference]
                    .into_iter()
                    .enumerate()
            {
                let mut time = controls(target, 1.0);
                time.maximum_step_s = 2.0 * target;
                let factor = if second_order { 1.0 / 3.0 } else { 1.0 };
                let tolerance = factor * difference / 0.2;
                match coordinate {
                    0 => time.moisture_absolute_tolerance = tolerance,
                    1 => time.activity_absolute_tolerance = tolerance,
                    _ => time.reservoir_water_absolute_tolerance_kg = tolerance,
                }
                let inspect = |next_air: &FiniteMoistureAtmosphere, next: &MoistureRadial| {
                    near(next.moisture()[0], fine, 3e-11);
                    near(
                        next_air.relative_humidity().unwrap(),
                        rh - dry * (fine - initial) / capacity,
                        3e-9,
                    );
                };
                if second_order {
                    let (next_air, next, report) = air
                        .advance_radial_sdirk2_adaptive_to(
                            &original,
                            target,
                            &faces,
                            config(),
                            nonlinear(),
                            finite_controls(),
                            time,
                            cx,
                        )
                        .unwrap();
                    inspect(&next_air, &next);
                    check_account(&air, &original, &next_air, &next, &report);
                    near(report.maximum_accepted_error_ratio, 0.2, 3e-7);
                    assert_eq!(report.accepted_intervals, 1);
                    assert_eq!(report.evaluations, 6);
                    assert_eq!(report.suggested_step_s, target);
                    assert_eq!(report.rejected_intervals, 0);
                } else {
                    let (next_air, next, report) = air
                        .advance_radial_adaptive_to(
                            &original,
                            target,
                            &faces,
                            config(),
                            nonlinear(),
                            finite_controls(),
                            time,
                            cx,
                        )
                        .unwrap();
                    inspect(&next_air, &next);
                    check_account(&air, &original, &next_air, &next, &report);
                    near(report.maximum_accepted_error_ratio, 0.2, 3e-7);
                    assert_eq!(report.accepted_intervals, 1);
                    assert_eq!(report.evaluations, 3);
                    assert_eq!(report.suggested_step_s, 2.0 * target);
                    assert_eq!(report.rejected_intervals, 0);
                }
            }
        }
    });
}

#[test]
fn finite_radial_adaptive_work_limit_publishes_neither_accepted_prefix() {
    with_cx(|cx| {
        let original = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.05],
            vec![0.1],
            [0.001; 2],
            cx,
        );
        let air = gas(0.8, 1.0);
        let original_air = air.clone();
        let faces = [None, Some(sorption_law(0.25))];
        let step = 1.0 / 256.0;
        for second_order in [false, true] {
            let budget = if second_order { 6 } else { 3 };
            let mut time = controls(step, 1.0);
            time.maximum_evaluations = budget;
            if second_order {
                assert!(matches!(
                    air.advance_radial_sdirk2_adaptive_to(
                        &original, 2.0 * step, &faces, config(), nonlinear(), finite_controls(), time, cx,
                    ),
                    Err(MoistureAdvanceError::WorkLimit { time_s, accepted_intervals: 1, evaluations })
                        if time_s == step && evaluations == budget
                ));
            } else {
                assert!(matches!(
                    air.advance_radial_adaptive_to(
                        &original, 2.0 * step, &faces, config(), nonlinear(), finite_controls(), time, cx,
                    ),
                    Err(MoistureAdvanceError::WorkLimit { time_s, accepted_intervals: 1, evaluations })
                        if time_s == step && evaluations == budget
                ));
            }
            assert_eq!(air, original_air);
            assert_eq!(original.moisture(), &[0.1]);
            assert_eq!(original.time_s(), 0.0);
        }
        let mut time = controls(step, 1.0);
        time.maximum_evaluations = 12;
        let (next_air, next, report) = air
            .advance_radial_sdirk2_adaptive_to(
                &original,
                2.0 * step,
                &faces,
                config(),
                nonlinear(),
                finite_controls(),
                time,
                cx,
            )
            .unwrap();
        assert_eq!(report.accepted_intervals, 2);
        assert_eq!(report.accepted_substeps, 4);
        assert_eq!(report.evaluations, 12);
        assert!(next.moisture()[0] > original.moisture()[0]);
        assert!(next_air.relative_humidity().unwrap() < air.relative_humidity().unwrap());
        check_account(&air, &original, &next_air, &next, &report);
        assert_eq!(air, original_air);
    });
}

#[test]
fn finite_radial_adaptive_cancellation_preserves_both_physical_inventories() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let original = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.05],
            vec![0.1],
            [0.001; 2],
            &cx,
        );
        let air = gas(0.8, 1.0);
        let original_air = air.clone();
        let faces = [None, Some(sorption_law(0.25))];
        gate.request();
        assert!(matches!(
            air.advance_radial_adaptive_to(
                &original,
                1.0,
                &faces,
                config(),
                nonlinear(),
                finite_controls(),
                controls(1.0, 1e-5),
                &cx,
            ),
            Err(MoistureAdvanceError::Cancelled)
        ));
        assert!(matches!(
            air.advance_radial_sdirk2_adaptive_to(
                &original,
                1.0,
                &faces,
                config(),
                nonlinear(),
                finite_controls(),
                controls(1.0, 1e-5),
                &cx,
            ),
            Err(MoistureAdvanceError::Cancelled)
        ));
        assert_eq!(air, original_air);
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}
