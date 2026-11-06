//! Independent coupled modes and opposed local transfers under shared SDIRK.
use super::*;

fn decay(z: f64) -> f64 {
    let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
    (1.0 - (1.0 - 2.0 * gamma) * z) / (1.0 + gamma * z).powi(2)
}

#[test]
fn g1_shared_sdirk_adaptive_resolves_coupled_mean_and_difference_modes() {
    with_cx(|cx| {
        let target = 2.0;
        let slope = 0.25;
        let air = atmosphere::atmosphere(0.8, 1.0);
        let initial_air = air.clone();
        let a = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.05], cx);
        let b = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.20], cx);
        let faces = [
            Some(atmosphere::sorption(slope)),
            Some(atmosphere::sorption(slope)),
            None,
            None,
            None,
            None,
        ];
        let participants = [spatial(8, &b, &faces), spatial(3, &a, &faces)];
        let rh = air.relative_humidity().unwrap();
        let capacity = air.vapor_mass_kg().unwrap() / rh;
        let dry = a.dry_mass_kg() + b.dry_mass_kg();
        let mean0 = 0.125;
        let k = 0.04;
        let scale = 1.0 + slope * dry / capacity;
        let equilibrium = (slope * rh + slope * dry * mean0 / capacity) / scale;
        // The mean couples to the ONE vapor store; the difference has no
        // reservoir debit. These are independent continuous-time eigenmodes.
        let mean = equilibrium + (mean0 - equilibrium) * (-k * scale * target).exp();
        let difference = 0.075 * (-k * target).exp();
        let expected = [mean - difference, mean + difference];
        let loose = air
            .advance_shared_cartesian_sdirk2_adaptive_to(
                &participants,
                target,
                2,
                atmosphere::controls(),
                time_controls(target, 2e-4),
                cx,
            )
            .unwrap();
        let time = time_controls(target, 2e-6);
        let fine = air
            .advance_shared_cartesian_sdirk2_adaptive_to(
                &participants,
                target,
                2,
                atmosphere::controls(),
                time,
                cx,
            )
            .unwrap();
        let error = |candidate: &SharedCartesianAtmosphereMoistureAdvanceCandidate<
            MoistureCartesianSdirk2Report,
        >| {
            candidate
                .solids()
                .iter()
                .zip(expected)
                .map(|(p, m)| (p.moisture.moisture()[0] - m).abs())
                .fold(0.0_f64, f64::max)
        };
        assert!(error(&fine) < 0.15 * error(&loose));
        assert!(error(&fine) < 1e-6);
        near(
            fine.reservoir().relative_humidity().unwrap(),
            rh - dry * (mean - mean0) / capacity,
            2e-5,
        );
        let report = fine.report();
        assert!(report.summary.rejected_intervals > 0);
        assert!(report.summary.evaluations > loose.report().summary.evaluations);
        assert_eq!(report.summary.evaluations, 6 * report.summary.comparisons);
        assert_eq!(
            report.summary.accepted_substeps,
            2 * report.summary.accepted_intervals
        );
        assert_eq!(report.summary.interval_s, [0.0, target]);
        assert_eq!(report.summary.last_substep.interval_s[1], target);
        assert_eq!(
            report.specimens.iter().map(|p| p.id).collect::<Vec<_>>(),
            vec![3, 8]
        );
        for (i, (field, account)) in fine.solids().iter().zip(&report.specimens).enumerate() {
            assert_eq!(field.id, account.id);
            assert_eq!(field.moisture.time_s(), target);
            assert_eq!(field.transport.endpoint.interval_s[1], target);
            assert!(field.transport.first_stage.interval_s[1] < target);
            assert_eq!(&account.boundary_water_kg[2..], &[0.0; 4]);
            assert_eq!(account.water_change_kg.is_sign_positive(), i == 0);
            let change =
                field.moisture.dry_mass_kg() * (field.moisture.moisture()[0] - [0.05, 0.20][i]);
            near(account.water_change_kg, change, 1e-13);
            near(
                account.boundary_water_kg.iter().sum(),
                change,
                time.water_balance_tolerance_kg,
            );
            near(
                field.transport.endpoint.transport_end[0].query().axes()["moisture_content"],
                field.moisture.moisture()[0],
                0.0,
            );
            for source in field.transport.endpoint.surface_sorption.iter().flatten() {
                near(
                    source.query().axes()["relative_humidity"],
                    fine.reservoir().relative_humidity().unwrap(),
                    1e-10,
                );
            }
        }
        let before = air.vapor_mass_kg().unwrap() + a.water_mass_kg() + b.water_mass_kg();
        let after = fine.reservoir().vapor_mass_kg().unwrap()
            + fine
                .solids()
                .iter()
                .map(|p| p.moisture.water_mass_kg())
                .sum::<f64>();
        near(after, before, time.water_balance_tolerance_kg);
        near(
            report.summary.water_balance_residual_kg,
            0.0,
            time.water_balance_tolerance_kg,
        );
        near(
            atmosphere::carrier_mass(fine.reservoir()),
            atmosphere::carrier_mass(&air),
            1e-12,
        );
        assert!(
            report.summary.maximum_carrier_mass_residual_kg.unwrap()
                <= time.water_balance_tolerance_kg
        );
        assert!(
            report.summary.total_mass_residual_kg.unwrap().abs() <= time.water_balance_tolerance_kg
        );

        let reversed = air
            .advance_shared_cartesian_sdirk2_adaptive_to(
                &[spatial(3, &a, &faces), spatial(8, &b, &faces)],
                target,
                2,
                atmosphere::controls(),
                time,
                cx,
            )
            .unwrap();
        assert_eq!(fine.reservoir(), reversed.reservoir());
        assert_eq!(
            report.summary.evaluations,
            reversed.report().summary.evaluations
        );
        for ((left, right), (a, b)) in fine
            .solids()
            .iter()
            .zip(reversed.solids())
            .zip(report.specimens.iter().zip(&reversed.report().specimens))
        {
            assert_eq!(left.id, right.id);
            assert_eq!(left.moisture.moisture(), right.moisture.moisture());
            assert_eq!(a.boundary_water_kg, b.boundary_water_kg);
        }
        assert_eq!(air, initial_air);
        assert_eq!(a.moisture(), &[0.05]);
        assert_eq!(b.moisture(), &[0.20]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}

#[test]
fn g1_shared_sdirk_adaptive_keeps_opposed_specimen_errors_when_vapor_is_stationary() {
    with_cx(|cx| {
        let target = 10.0;
        let air = atmosphere::atmosphere(0.8, 1.0);
        let a = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.1], cx);
        let b = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.3], cx);
        let faces = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
            None,
            None,
            None,
            None,
        ];
        // Equal dry masses and opposite deviations from M_eq=.2 leave vapor
        // unchanged, but each field has a nonzero local truncation error.
        let coarse = 0.1 * decay(0.04 * target);
        let fine = 0.1 * decay(0.02 * target).powi(2);
        let mut time = time_controls(target, (coarse - fine).abs() / (3.0 * 0.2));
        time.activity_absolute_tolerance = 1e-8;
        time.reservoir_water_absolute_tolerance_kg = 1e-10;
        time.maximum_step_s = 2.0 * target;
        let endpoint = air
            .advance_shared_cartesian_sdirk2_adaptive_to(
                &[spatial(7, &a, &faces), spatial(2, &b, &faces)],
                target,
                2,
                atmosphere::controls(),
                time,
                cx,
            )
            .unwrap();
        let report = endpoint.report();
        assert_eq!(report.summary.evaluations, 6);
        assert_eq!(report.summary.accepted_intervals, 1);
        assert_eq!(report.summary.accepted_substeps, 2);
        assert_eq!(report.summary.rejected_intervals, 0);
        near(report.summary.maximum_accepted_error_ratio, 0.2, 2e-7);
        assert_eq!(report.summary.suggested_step_s, target);
        near(
            endpoint.reservoir().relative_humidity().unwrap(),
            air.relative_humidity().unwrap(),
            1e-12,
        );
        near(
            endpoint.reservoir().vapor_mass_kg().unwrap(),
            air.vapor_mass_kg().unwrap(),
            1e-13,
        );
        near(
            report.summary.water_change_kg,
            0.0,
            time.water_balance_tolerance_kg,
        );
        for (i, (field, account)) in endpoint.solids().iter().zip(&report.specimens).enumerate() {
            let sign = if i == 0 { 1.0 } else { -1.0 };
            near(field.moisture.moisture()[0], 0.2 + sign * fine, 2e-12);
            near(account.maximum_accepted_error_ratio, 0.2, 2e-7);
            let water = -sign * field.moisture.dry_mass_kg() * (0.1 - fine);
            near(account.water_change_kg, water, 2e-12);
            assert!(water.abs() > 0.01);
            near(account.boundary_water_kg[0], water / 2.0, 2e-12);
            near(account.boundary_water_kg[1], water / 2.0, 2e-12);
            assert_eq!(&account.boundary_water_kg[2..], &[0.0; 4]);
        }
        near(
            report.specimens.iter().map(|p| p.water_change_kg).sum(),
            0.0,
            1e-12,
        );
        assert_eq!(a.moisture(), &[0.1]);
        assert_eq!(b.moisture(), &[0.3]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}

#[test]
fn g4_shared_sdirk_adaptive_reserves_complete_stages_and_discards_all_private_fields() {
    with_cx(|cx| {
        let air = atmosphere::atmosphere(0.8, 1.0);
        let original_air = air.clone();
        let a = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.05], cx);
        let b = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.20], cx);
        let faces = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
            None,
            None,
            None,
            None,
        ];
        let participants = [spatial(8, &b, &faces), spatial(3, &a, &faces)];
        let mut time = time_controls(1.0, 1.0);
        for (budget, evaluated, accepted, failed_at) in [(5, 4, 0, 0.5), (11, 10, 1, 1.5)] {
            time.maximum_evaluations = budget;
            assert!(matches!(
                air.advance_shared_cartesian_sdirk2_adaptive_to(
                    &participants, 2.0, 2, atmosphere::controls(), time, cx,
                ),
                Err(MoistureAdvanceError::WorkLimit {
                    evaluations, accepted_intervals, time_s,
                }) if evaluations == evaluated && accepted_intervals == accepted && time_s == failed_at
            ));
        }
        time.maximum_evaluations = 12;
        let complete = air
            .advance_shared_cartesian_sdirk2_adaptive_to(
                &participants,
                2.0,
                2,
                atmosphere::controls(),
                time,
                cx,
            )
            .unwrap();
        assert_eq!(complete.report().summary.evaluations, 12);
        assert_eq!(complete.report().summary.accepted_intervals, 2);
        assert_eq!(complete.report().summary.accepted_substeps, 4);
        assert_eq!(complete.report().summary.rejected_intervals, 0);
        assert_eq!(
            complete.solids().iter().map(|p| p.id).collect::<Vec<_>>(),
            vec![3, 8]
        );
        assert!(complete.solids().iter().all(|p| p.moisture.time_s() == 2.0));
        assert!(complete.solids()[0].moisture.moisture()[0] > 0.05);
        assert!(complete.solids()[1].moisture.moisture()[0] < 0.20);
        assert!(matches!(
            air.advance_shared_cartesian_sdirk2_adaptive_to(
                &participants,
                2.0,
                1,
                atmosphere::controls(),
                time,
                cx,
            ),
            Err(MoistureAdvanceError::Invalid(_))
        ));
        assert_eq!(air, original_air);
        assert_eq!(a.moisture(), &[0.05]);
        assert_eq!(b.moisture(), &[0.20]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}
