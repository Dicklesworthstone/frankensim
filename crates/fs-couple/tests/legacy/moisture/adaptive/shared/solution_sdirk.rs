//! G1 independent specimen errors and finite solution activity under shared SDIRK.
use super::*;

#[test]
fn g1_shared_solution_sdirk_adaptive_retains_opposed_fields_and_full_stage_costs() {
    with_cx(|cx| {
        let target = 10.0;
        let volume = 1e-3;
        let bath = solution::reservoir(volume, vec![0.0]);
        let original_solution = bath.solution().clone();
        let a = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.1], cx);
        let b = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.3], cx);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
            None,
            None,
            None,
            None,
        ];
        let participants = [spatial(7, &a, &faces), spatial(2, &b, &faces)];
        // The mean and pure-solvent inventory remain stationary. Each genuine
        // one-cell specimen still has its own nonzero mode with k=4D/Lx².
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        let decay = |z: f64| (1.0 - (1.0 - 2.0 * gamma) * z) / (1.0 + gamma * z).powi(2);
        let coarse = 0.1 * decay(0.04 * target);
        let fine = 0.1 * decay(0.02 * target).powi(2);
        let mut time = time_controls(target, (coarse - fine).abs() / (3.0 * 0.2));
        time.maximum_step_s = 2.0 * target;
        time.activity_absolute_tolerance = 1e-8;
        time.reservoir_water_absolute_tolerance_kg = 1e-10;
        let endpoint = bath
            .advance_shared_cartesian_sdirk2_adaptive_to(
                &participants,
                target,
                2,
                solution::controls(volume),
                time,
                cx,
            )
            .unwrap();
        let report = endpoint.report();
        assert_eq!(report.summary.interval_s, [0.0, target]);
        assert_eq!(report.summary.comparisons, 1);
        assert_eq!(report.summary.evaluations, 6);
        assert_eq!(report.summary.accepted_intervals, 1);
        assert_eq!(report.summary.accepted_substeps, 2);
        assert_eq!(report.summary.rejected_intervals, 0);
        assert_eq!(report.summary.activity, Some([1.0; 2]));
        near(report.summary.maximum_accepted_error_ratio, 0.2, 2e-7);
        // Cubic local error admits this span but must not double at ratio 0.2.
        assert_eq!(report.summary.suggested_step_s, target);
        assert_eq!(
            report.specimens.iter().map(|p| p.id).collect::<Vec<_>>(),
            vec![2, 7]
        );
        for (i, (field, account)) in endpoint.solids().iter().zip(&report.specimens).enumerate() {
            let sign = if i == 0 { 1.0 } else { -1.0 };
            let water = -sign * field.moisture.dry_mass_kg() * (0.1 - fine);
            assert_eq!(field.id, account.id);
            assert_eq!(field.moisture.time_s(), target);
            assert_eq!(field.transport.endpoint.interval_s[1], target);
            near(field.moisture.moisture()[0], 0.2 + sign * fine, 2e-12);
            near(account.maximum_accepted_error_ratio, 0.2, 2e-7);
            near(account.water_change_kg, water, 2e-12);
            near(
                account.stored_water_kg[1] - account.stored_water_kg[0],
                water,
                2e-12,
            );
            assert!(account.water_change_kg.abs() > 0.01);
            for transfer in &account.boundary_water_kg[..2] {
                near(*transfer, 0.5 * water, 2e-12);
            }
            assert_eq!(&account.boundary_water_kg[2..], &[0.0; 4]);
        }
        near(report.summary.water_change_kg, 0.0, 1e-12);
        near(
            endpoint.reservoir().solvent_mass_kg(),
            bath.solvent_mass_kg(),
            1e-12,
        );
        near(
            endpoint.reservoir().solution().solvent_volume_m3(),
            volume,
            1e-15,
        );
        assert_eq!(endpoint.reservoir().solution().solute_moles(), &[0.0]);

        // One complete comparison costs six stages. A second coarse and
        // first half consume four more; the final half cannot fit into eleven.
        let mut short = time_controls(1.0, 1.0);
        short.maximum_evaluations = 11;
        assert!(matches!(
            bath.advance_shared_cartesian_sdirk2_adaptive_to(
                &participants,
                2.0,
                2,
                solution::controls(volume),
                short,
                cx,
            ),
            Err(MoistureAdvanceError::WorkLimit {
                evaluations: 10,
                accepted_intervals: 1,
                time_s,
            }) if time_s == 1.5
        ));
        assert_eq!(bath.solution(), &original_solution);
        assert_eq!(a.moisture(), &[0.1]);
        assert_eq!(b.moisture(), &[0.3]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}

#[test]
fn g1_shared_solution_sdirk_adaptive_refines_finite_brine_and_preserves_solute() {
    with_cx(|cx| {
        let volume = 1e-4;
        let initial_activity = 0.95_f64;
        let total_moles = -initial_activity.ln() * volume / solution::WATER_VOLUME;
        let moles = vec![0.6 * total_moles, 0.4 * total_moles];
        let bath = solution::reservoir(volume, moles.clone());
        let original_solution = bath.solution().clone();
        let a = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.05], cx);
        let b = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.20], cx);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
            None,
            None,
            None,
            None,
        ];
        let participants = [spatial(8, &b, &faces), spatial(3, &a, &faces)];
        let target = 10.0;
        let loose = bath
            .advance_shared_cartesian_sdirk2_adaptive_to(
                &participants,
                target,
                2,
                solution::controls(volume),
                time_controls(target, 1e-5),
                cx,
            )
            .unwrap();
        let time = time_controls(target, 1e-8);
        let fine = bath
            .advance_shared_cartesian_sdirk2_adaptive_to(
                &participants,
                target,
                2,
                solution::controls(volume),
                time,
                cx,
            )
            .unwrap();
        // Independent continuous reduction: total dry mass is 1 kg, water
        // density 1000 kg/m³, V(m)=1e-4-(m-0.125)/1000 and A=-ln(0.95)*1e-4.
        // The mean obeys m'=0.04*(0.2*exp(-A/V(m))-m).
        // RK4 at 4096/8192 steps agrees within 4e-16; inverse quadrature
        // integral(dm/m') independently returns 10 s within 1e-13.
        let mean = 0.145_977_933_390_675_77;
        let difference = 0.075 * (-0.04_f64 * target).exp();
        let expected = [mean - difference, mean + difference];
        let error = |candidate: &SharedCartesianSolutionMoistureAdvanceCandidate<
            MoistureCartesianSdirk2Report,
        >| {
            candidate
                .solids()
                .iter()
                .zip(expected)
                .map(|(p, m)| (p.moisture.moisture()[0] - m).abs())
                .fold(0.0_f64, f64::max)
        };
        assert!(error(&fine) < 0.05 * error(&loose));
        assert!(error(&fine) < 1e-7);
        let report = fine.report();
        assert_eq!(report.summary.interval_s, [0.0, target]);
        assert_eq!(report.summary.last_substep.interval_s[1], target);
        assert!(report.summary.rejected_intervals > 0);
        assert!(report.summary.evaluations > loose.report().summary.evaluations);
        assert_eq!(report.summary.evaluations, 6 * report.summary.comparisons);
        assert_eq!(
            report.summary.accepted_substeps,
            2 * report.summary.accepted_intervals
        );
        assert_eq!(
            report.specimens.iter().map(|p| p.id).collect::<Vec<_>>(),
            vec![3, 8]
        );
        let final_solution = fine.reservoir().solution();
        let activity = final_solution.water_activity().unwrap();
        assert!(activity < initial_activity - 1e-3);
        assert!(final_solution.solvent_volume_m3() < 0.95 * volume);
        assert_eq!(final_solution.solute_moles(), moles.as_slice());
        near(
            activity,
            (-solution::WATER_VOLUME * total_moles / final_solution.solvent_volume_m3()).exp(),
            2e-14,
        );
        near(activity, 0.937_151_719_100_423_5, 1e-7);
        let activity_account = report.summary.activity.unwrap();
        near(activity_account[0], initial_activity, 1e-14);
        near(activity_account[1], activity, 0.0);
        let reservoir_account = report.summary.reservoir_water_kg.unwrap();
        near(reservoir_account[0], bath.solvent_mass_kg(), 0.0);
        near(
            reservoir_account[1],
            fine.reservoir().solvent_mass_kg(),
            0.0,
        );
        let mut water_before = bath.solvent_mass_kg();
        let mut water_after = fine.reservoir().solvent_mass_kg();
        for (i, (field, account)) in fine.solids().iter().zip(&report.specimens).enumerate() {
            let before = [0.05, 0.20][i];
            let after = field.moisture.moisture()[0];
            let dry = field.moisture.dry_mass_kg();
            assert_eq!(field.id, account.id);
            assert_eq!(field.moisture.time_s(), target);
            assert_eq!(field.transport.endpoint.interval_s[1], target);
            assert_eq!(account.water_change_kg.is_sign_positive(), i == 0);
            assert!(account.maximum_accepted_error_ratio > 0.0);
            assert!(account.maximum_accepted_error_ratio <= 1.0);
            near(account.water_change_kg, dry * (after - before), 1e-14);
            near(
                account.boundary_water_kg.iter().sum(),
                account.water_change_kg,
                time.water_balance_tolerance_kg,
            );
            assert_eq!(&account.boundary_water_kg[2..], &[0.0; 4]);
            for source in field.transport.endpoint.surface_sorption.iter().flatten() {
                near(source.query().axes()["relative_humidity"], activity, 1e-10);
            }
            water_before += dry * before;
            water_after += dry * after;
        }
        near(water_after, water_before, time.water_balance_tolerance_kg);
        near(
            report.summary.water_balance_residual_kg,
            0.0,
            time.water_balance_tolerance_kg,
        );
        assert_eq!(bath.solution(), &original_solution);
        assert_eq!(a.moisture(), &[0.05]);
        assert_eq!(b.moisture(), &[0.20]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}
