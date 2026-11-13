//! Coupled continuous and discrete two-inventory references for adaptive SDIRK.
//! The isotherm, diffusion and gas packs are authored verification inputs.
use super::*;
use fs_couple::moisture::atmosphere::AtmosphereMoistureError;

fn decay(z: f64) -> f64 {
    let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
    (1.0 - (1.0 - 2.0 * gamma) * z) / (1.0 + gamma * z).powi(2)
}

#[test]
fn g1_finite_sdirk_adaptive_controls_moisture_humidity_and_vapor() {
    with_cx(|cx| {
        let initial = 0.1;
        let slope = 0.25;
        let original = grid([0.1; 3], [1; 3], 0.01, vec![initial], cx);
        let air = atmosphere::atmosphere(0.8, 1.0);
        let original_air = air.clone();
        let rh = air.relative_humidity().unwrap();
        let capacity = air.vapor_mass_kg().unwrap() / rh;
        let dry = original.dry_mass_kg();
        let carrier = atmosphere::carrier_mass(&air);
        let faces = std::array::from_fn(|_| Some(atmosphere::sorption(slope)));
        // Eliminate vapor with C*RH + dry*M = constant. Six half-cell
        // faces give k=12D/L^2; the closed two-store mode has rate k*b.
        let b = 1.0 + slope * dry / capacity;
        let equilibrium = (slope * rh + slope * dry * initial / capacity) / b;
        let lambda = 12.0 * 0.01 / 0.1_f64.powi(2) * b;
        let target = 0.4 / lambda;
        let coarse = equilibrium + (initial - equilibrium) * decay(lambda * target);
        let fine = equilibrium + (initial - equilibrium) * decay(lambda * target / 2.0).powi(2);
        let difference = (coarse - fine).abs();
        let differences = [difference, dry * difference / capacity, dry * difference];

        // Make each physical coordinate independently dominate the estimate;
        // an omitted RH or vapor observable cannot pass through a content-only test.
        for (coordinate, difference) in differences.into_iter().enumerate() {
            let mut time = time_controls(target, 1.0);
            time.maximum_step_s = 2.0 * target;
            let tolerance = difference / (3.0 * 0.2);
            match coordinate {
                0 => time.moisture_absolute_tolerance = tolerance,
                1 => time.activity_absolute_tolerance = tolerance,
                _ => time.reservoir_water_absolute_tolerance_kg = tolerance,
            }
            let (next_air, next, report) = air
                .advance_cartesian_sdirk2_adaptive_to(
                    &original,
                    target,
                    &faces,
                    config(),
                    nonlinear(),
                    atmosphere::controls(),
                    time,
                    cx,
                )
                .unwrap();
            near(next.moisture()[0], fine, 2e-11);
            near(
                next_air.relative_humidity().unwrap(),
                rh - dry * (fine - initial) / capacity,
                2e-9,
            );
            near(report.maximum_accepted_error_ratio, 0.2, 2e-7);
            assert_eq!(report.comparisons, 1);
            assert_eq!(report.accepted_intervals, 1);
            assert_eq!(report.accepted_substeps, 2);
            assert_eq!(report.evaluations, 6);
            assert_eq!(report.rejected_intervals, 0);
            // A ratio of .2 must not use the first-order .25 growth threshold.
            assert_eq!(report.suggested_step_s, target);
            assert_eq!(next.time_s(), target);
            assert_eq!(report.interval_s, [0.0, target]);
            assert_eq!(
                report.last_substep.transport.endpoint.interval_s,
                [target / 2.0, target]
            );
            assert_eq!(
                report.last_substep.transport.first_stage.interval_s[0],
                target / 2.0
            );
            let water = dry * (fine - initial);
            for face_water in report.boundary_water_kg {
                near(face_water, water / 6.0, 1e-11);
            }
            near(report.water_change_kg, water, 1e-11);
            near(
                report.water_balance_residual_kg,
                0.0,
                time.water_balance_tolerance_kg,
            );
            near(
                report.transport_balance_residual_kg,
                0.0,
                time.water_balance_tolerance_kg,
            );
            near(
                report.reservoir_water_kg.unwrap()[1],
                next_air.vapor_mass_kg().unwrap(),
                0.0,
            );
            near(
                report.activity.unwrap()[1],
                next_air.relative_humidity().unwrap(),
                0.0,
            );
            near(atmosphere::carrier_mass(&next_air), carrier, 1e-12);
            assert!(
                report.maximum_carrier_mass_residual_kg.unwrap() <= time.water_balance_tolerance_kg
            );
            assert!(
                report.total_mass_residual_kg.unwrap().abs() <= time.water_balance_tolerance_kg
            );
            if coordinate == 2 {
                time.reservoir_water_absolute_tolerance_kg *= 2.0;
                let (_, _, growing) = air
                    .advance_cartesian_sdirk2_adaptive_to(
                        &original,
                        target,
                        &faces,
                        config(),
                        nonlinear(),
                        atmosphere::controls(),
                        time,
                        cx,
                    )
                    .unwrap();
                near(growing.maximum_accepted_error_ratio, 0.1, 1e-7);
                assert_eq!(growing.suggested_step_s, 2.0 * target);
            }
        }
        assert_eq!(air, original_air);
        assert_eq!(original.moisture(), &[initial]);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn g4_finite_sdirk_adaptive_recovers_stiff_stage_range_without_clipping() {
    with_cx(|cx| {
        let initial = 0.2;
        let slope = 0.25;
        let target = 0.5;
        let original = grid([0.01; 3], [1; 3], 1e-4, vec![initial], cx);
        let air = atmosphere::atmosphere(0.05, 1.0);
        let original_air = air.clone();
        let rh = air.relative_humidity().unwrap();
        let capacity = air.vapor_mass_kg().unwrap() / rh;
        let dry = original.dry_mass_kg();
        let b = 1.0 + slope * dry / capacity;
        let equilibrium = (slope * rh + slope * dry * initial / capacity) / b;
        let lambda = 12.0 * 1e-4 / 0.01_f64.powi(2) * b;
        let expected = equilibrium + (initial - equilibrium) * (-lambda * target).exp();
        let faces = std::array::from_fn(|_| Some(atmosphere::sorption(slope)));
        let mut time = time_controls(target, 1e-5);
        let (next_air, next, report) = air
            .advance_cartesian_sdirk2_adaptive_to(
                &original,
                target,
                &faces,
                config(),
                nonlinear(),
                atmosphere::controls(),
                time,
                cx,
            )
            .unwrap();
        // The initial span has k*dt=6 and its second stage is inadmissible.
        // A positive endpoint follows the continuous coupled decay after retries.
        near(next.moisture()[0], expected, 1e-4);
        assert!(next.moisture()[0] > 0.0);
        assert!(next_air.relative_humidity().unwrap() > rh);
        assert_eq!(next.time_s(), target);
        assert!(report.rejected_intervals > 0);
        assert!(report.evaluations < 6 * report.comparisons);
        assert_eq!(report.evaluations % 2, 0);
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        near(
            next_air.vapor_mass_kg().unwrap() + next.water_mass_kg(),
            air.vapor_mass_kg().unwrap() + original.water_mass_kg(),
            time.water_balance_tolerance_kg,
        );
        near(
            report.boundary_water_kg.iter().sum(),
            report.water_change_kg,
            time.water_balance_tolerance_kg,
        );
        near(
            report.water_change_kg,
            dry * (next.moisture()[0] - initial),
            1e-13,
        );
        time.minimum_step_s = target;
        assert!(matches!(
            air.advance_cartesian_sdirk2_adaptive_to(
                &original,
                target,
                &faces,
                config(),
                nonlinear(),
                atmosphere::controls(),
                time,
                cx,
            ),
            Err(MoistureAdvanceError::Step {
                evaluations: 2,
                error: AtmosphereMoistureError::Moisture(MoistureError::StageRange { .. }),
                ..
            })
        ));
        assert_eq!(air, original_air);
        assert_eq!(original.moisture(), &[initial]);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn g4_finite_sdirk_adaptive_reserves_stages_and_discards_private_prefixes() {
    with_cx(|cx| {
        let initial = 0.1;
        let slope = 0.25;
        let step = 1.0 / 1024.0;
        let original = grid([0.1; 3], [1; 3], 0.01, vec![initial], cx);
        let air = atmosphere::atmosphere(0.8, 1.0);
        let original_air = air.clone();
        let faces = std::array::from_fn(|_| Some(atmosphere::sorption(slope)));
        let mut time = time_controls(step, 1.0);
        for (budget, evaluated, accepted, failed_at) in
            [(5, 4, 0, 0.5 * step), (11, 10, 1, 1.5 * step)]
        {
            time.maximum_evaluations = budget;
            assert!(matches!(
                air.advance_cartesian_sdirk2_adaptive_to(
                    &original, 2.0 * step, &faces, config(), nonlinear(),
                    atmosphere::controls(), time, cx,
                ),
                Err(MoistureAdvanceError::WorkLimit {
                    evaluations, accepted_intervals, time_s,
                }) if evaluations == evaluated && accepted_intervals == accepted && time_s == failed_at
            ));
        }
        time.maximum_evaluations = 12;
        let (next_air, next, report) = air
            .advance_cartesian_sdirk2_adaptive_to(
                &original,
                2.0 * step,
                &faces,
                config(),
                nonlinear(),
                atmosphere::controls(),
                time,
                cx,
            )
            .unwrap();
        let rh = air.relative_humidity().unwrap();
        let capacity = air.vapor_mass_kg().unwrap() / rh;
        let dry = original.dry_mass_kg();
        let b = 1.0 + slope * dry / capacity;
        let equilibrium = (slope * rh + slope * dry * initial / capacity) / b;
        let lambda = 12.0 * 0.01 / 0.1_f64.powi(2) * b;
        // Four accepted half intervals evolve the SAME finite reservoir.
        let expected = equilibrium + (initial - equilibrium) * decay(lambda * step / 2.0).powi(4);
        near(next.moisture()[0], expected, 2e-11);
        near(
            next_air.relative_humidity().unwrap(),
            rh - dry * (expected - initial) / capacity,
            2e-9,
        );
        assert_eq!(report.evaluations, 12);
        assert_eq!(report.accepted_intervals, 2);
        assert_eq!(report.accepted_substeps, 4);
        assert_eq!(report.rejected_intervals, 0);
        assert_eq!(next.time_s(), 2.0 * step);
        near(
            report.boundary_water_kg.iter().sum(),
            report.water_change_kg,
            time.water_balance_tolerance_kg,
        );
        assert_eq!(air, original_air);
        assert_eq!(original.moisture(), &[initial]);
        assert_eq!(original.time_s(), 0.0);
    });
}
