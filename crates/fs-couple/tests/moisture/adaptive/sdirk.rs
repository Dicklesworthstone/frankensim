//! Independent decay and physical-work checks for the shared SDIRK controller.
use super::*;

fn decay(z: f64) -> f64 {
    let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
    (1.0 - (1.0 - 2.0 * gamma) * z) / (1.0 + gamma * z).powi(2)
}

#[test]
fn sdirk_adaptive_uses_second_order_error_and_combined_endpoint_water() {
    with_cx(|cx| {
        let original = slab(card([0.1; 2], 0.3), 1.0, vec![0.1], cx);
        // One cell, two prescribed half-cell faces: lambda=4D/L^2=.4.
        let coarse = 0.2 - 0.1 * decay(0.4);
        let fine = 0.2 - 0.1 * decay(0.2).powi(2);
        let mut controls = time_controls(1.0, (coarse - fine).abs() / (3.0 * 0.2));
        controls.maximum_step_s = 2.0;
        let (next, report) = original
            .advance_sdirk2_adaptive_to(1.0, &faces(0.2), config(), nonlinear(), controls, cx)
            .unwrap();
        near(next.moisture()[0], fine, 2e-13);
        near(report.maximum_accepted_error_ratio, 0.2, 2e-10);
        assert_eq!(report.comparisons, 1);
        assert_eq!(report.accepted_intervals, 1);
        assert_eq!(report.accepted_substeps, 2);
        assert_eq!(report.evaluations, 6);
        assert_eq!(report.rejected_intervals, 0);
        // Doubling a second-order step magnifies local error eightfold, so a
        // ratio of .2 must not trigger the first-order .25 growth threshold.
        assert_eq!(report.suggested_step_s, 1.0);
        assert_eq!(report.last_substep.endpoint.interval_s, [0.5, 1.0]);
        assert_eq!(report.last_substep.first_stage.interval_s[0], 0.5);
        let water = original.dry_mass_kg() * (fine - 0.1);
        for actual in report.boundary_water_kg {
            near(actual, water / 2.0, 2e-13);
        }
        near(report.water_change_kg, water, 2e-13);
        near(report.water_balance_residual_kg, 0.0, 2e-13);
        // A local ratio of .1 permits doubling under the cubic error model.
        controls.moisture_absolute_tolerance *= 2.0;
        let (_, growing) = original
            .advance_sdirk2_adaptive_to(1.0, &faces(0.2), config(), nonlinear(), controls, cx)
            .unwrap();
        assert_eq!(growing.suggested_step_s, 2.0);
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn sdirk_adaptive_refines_to_the_discrete_diffusion_semigroup() {
    with_cx(|cx| {
        let n = 8;
        let diffusivity = 0.1;
        let wave = |i: usize| (core::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos();
        let original = slab(
            card([diffusivity; 2], 0.3),
            1.0,
            (0..n).map(|i| 0.15 + 0.05 * wave(i)).collect(),
            cx,
        );
        let sealed = [MoistureFace::Sealed, MoistureFace::Sealed];
        let lambda = 4.0
            * diffusivity
            * (n * n) as f64
            * (core::f64::consts::PI / (2.0 * n as f64)).sin().powi(2);
        let error = |state: &MoistureSlab| {
            state
                .moisture()
                .iter()
                .enumerate()
                .map(|(i, &m)| (m - (0.15 + 0.05 * (-lambda).exp() * wave(i))).abs())
                .fold(0.0_f64, f64::max)
        };
        let (loose, loose_report) = original
            .advance_sdirk2_adaptive_to(
                1.0,
                &sealed,
                config(),
                nonlinear(),
                time_controls(1.0, 1e-4),
                cx,
            )
            .unwrap();
        let tight = time_controls(1.0, 1e-7);
        let (fine, report) = original
            .advance_sdirk2_adaptive_to(1.0, &sealed, config(), nonlinear(), tight, cx)
            .unwrap();
        assert!(error(&fine) < 0.1 * error(&loose));
        assert!(error(&fine) < 1e-5);
        assert!(report.evaluations > loose_report.evaluations);
        assert!(report.rejected_intervals > 0);
        assert_eq!(report.evaluations, 6 * report.comparisons);
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        assert!(report.maximum_accepted_error_ratio <= 1.0);
        assert_eq!(report.boundary_water_kg, [0.0; 2]);
        assert_eq!(report.interval_s, [0.0, 1.0]);
        assert_eq!(fine.time_s(), 1.0);
        near(
            report.water_change_kg,
            0.0,
            tight.water_balance_tolerance_kg,
        );
        let (replay, replayed) = original
            .advance_sdirk2_adaptive_to(1.0, &sealed, config(), nonlinear(), tight, cx)
            .unwrap();
        assert_eq!(replay.moisture(), fine.moisture());
        assert_eq!(replayed.evaluations, report.evaluations);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn sdirk_adaptive_shortens_stage_range_refusals_without_clipping() {
    with_cx(|cx| {
        let original = slab(card([0.01; 2], 0.3), 1.0, vec![0.2], cx);
        let mut controls = time_controls(100.0, 1e-5);
        let (next, report) = original
            .advance_sdirk2_adaptive_to(100.0, &faces(0.0), config(), nonlinear(), controls, cx)
            .unwrap();
        // z=4 on the first coarse probe exceeds the SDIRK monotonicity range.
        // Refinement recovers the independent continuous decay, without a BE
        // fallback or publishing the inadmissible stage.
        near(next.moisture()[0], 0.2 * (-4.0_f64).exp(), 1e-4);
        assert!(next.moisture()[0] > 0.0);
        assert_eq!(next.time_s(), 100.0);
        assert!(report.rejected_intervals > 0);
        assert!(report.evaluations < 6 * report.comparisons);
        assert_eq!(report.evaluations % 2, 0);
        near(
            report.boundary_water_kg.iter().sum(),
            original.dry_mass_kg() * (next.moisture()[0] - 0.2),
            controls.water_balance_tolerance_kg,
        );
        controls.minimum_step_s = 100.0;
        assert!(matches!(
            original.advance_sdirk2_adaptive_to(
                100.0,
                &faces(0.0),
                config(),
                nonlinear(),
                controls,
                cx,
            ),
            Err(MoistureAdvanceError::Step {
                evaluations: 2,
                error: MoistureError::StageRange { .. },
                ..
            })
        ));
        assert_eq!(original.moisture(), &[0.2]);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn sdirk_adaptive_reserves_complete_stage_work_and_discards_private_prefixes() {
    with_cx(|cx| {
        let original = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let mut controls = time_controls(1.0, 1e-3);
        controls.maximum_evaluations = 5;
        assert!(matches!(
            original.advance_sdirk2_adaptive_to(
                1.0,
                &faces(0.1),
                config(),
                nonlinear(),
                controls,
                cx,
            ),
            Err(MoistureAdvanceError::WorkLimit {
                evaluations: 4,
                accepted_intervals: 0,
                ..
            })
        ));
        controls.maximum_evaluations = 11;
        assert!(matches!(
            original.advance_sdirk2_adaptive_to(
                2.0, &faces(0.1), config(), nonlinear(), controls, cx,
            ),
            Err(MoistureAdvanceError::WorkLimit {
                evaluations: 10, accepted_intervals: 1, time_s,
            }) if time_s == 1.5
        ));
        controls.maximum_evaluations = 12;
        let (next, report) = original
            .advance_sdirk2_adaptive_to(2.0, &faces(0.1), config(), nonlinear(), controls, cx)
            .unwrap();
        assert_eq!(report.evaluations, 12);
        assert_eq!(report.accepted_intervals, 2);
        assert_eq!(report.accepted_substeps, 4);
        assert_eq!(next.time_s(), 2.0);
        assert_eq!(next.moisture(), original.moisture());
        // Unsupported source/domain values remain terminal, charged once.
        assert!(matches!(
            original.advance_sdirk2_adaptive_to(
                1.0,
                &faces(0.31),
                config(),
                nonlinear(),
                controls,
                cx,
            ),
            Err(MoistureAdvanceError::Step {
                evaluations: 2,
                error: MoistureError::Invalid(_),
                ..
            })
        ));
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}
