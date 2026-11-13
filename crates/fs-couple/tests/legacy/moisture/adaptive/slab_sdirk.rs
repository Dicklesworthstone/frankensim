//! G1/G3/G4: native slab accuracy against independent finite-vapor modes.
//! Coefficients are authored numerical fixtures, not measured wood validation.
use super::*;
use fs_couple::moisture::atmosphere::{DryBasisSorption, MoistureParticipant};

fn participant<'a>(
    id: u64,
    slab: &'a MoistureSlab,
    faces: &'a [Option<DryBasisSorption>; 2],
) -> MoistureParticipant<'a> {
    MoistureParticipant {
        id,
        slab,
        faces,
        transport: config(),
    }
}

fn decay(z: f64) -> f64 {
    let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
    (1.0 - (1.0 - 2.0 * gamma) * z) / (1.0 + gamma * z).powi(2)
}

#[test]
fn g1_finite_slab_sdirk_adaptive_resolves_coupled_uptake_and_drying() {
    with_cx(|cx| {
        let slope = 0.25;
        let faces = [
            Some(atmosphere::sorption(slope)),
            Some(atmosphere::sorption(slope)),
        ];
        for rh in [0.8, 0.2] {
            for kirchhoff in [false, true] {
                let air = atmosphere::atmosphere(rh, 1.0);
                let original_air = air.clone();
                let mut original = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
                if kirchhoff {
                    original = original.with_kirchhoff_diffusivity(nonlinear()).unwrap();
                }
                let capacity = air.vapor_mass_kg().unwrap() / air.relative_humidity().unwrap();
                let dry = original.dry_mass_kg();
                let total = air.vapor_mass_kg().unwrap() + original.water_mass_kg();
                let equilibrium = slope * total / (capacity + slope * dry);
                let lambda = 0.04 * (1.0 + slope * dry / capacity);
                let target = 0.8 / lambda;
                let exact = equilibrium + (0.1 - equilibrium) * (-lambda * target).exp();
                let (_, loose, loose_report) = air
                    .advance_sdirk2_adaptive_to(
                        &original,
                        target,
                        &faces,
                        config(),
                        nonlinear(),
                        atmosphere::controls(),
                        time_controls(target, 2e-4),
                        cx,
                    )
                    .unwrap();
                let time = time_controls(target, 2e-6);
                let (next_air, fine, report) = air
                    .advance_sdirk2_adaptive_to(
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
                let error = (fine.moisture()[0] - exact).abs();
                assert!(error < 0.2 * (loose.moisture()[0] - exact).abs());
                assert!(error < 1e-6, "{error}");
                near(
                    next_air.relative_humidity().unwrap(),
                    (total - dry * exact) / capacity,
                    2e-5,
                );
                assert_eq!(fine.time_s(), target);
                assert_eq!(report.interval_s, [0.0, target]);
                assert!(report.evaluations > loose_report.evaluations);
                assert!(report.rejected_intervals > 0);
                assert_eq!(report.evaluations, 6 * report.comparisons);
                assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
                assert!(report.maximum_accepted_error_ratio <= 1.0);
                let water = dry * (fine.moisture()[0] - 0.1);
                near(report.water_change_kg, water, 1e-13);
                near(
                    report.boundary_water_kg.iter().sum(),
                    water,
                    time.water_balance_tolerance_kg,
                );
                near(
                    next_air.vapor_mass_kg().unwrap() + fine.water_mass_kg(),
                    total,
                    time.water_balance_tolerance_kg,
                );
                near(
                    atmosphere::carrier_mass(&next_air),
                    atmosphere::carrier_mass(&air),
                    1e-12,
                );
                assert_eq!(report.last_substep.transport.endpoint.interval_s[1], target);
                assert_eq!(
                    report
                        .last_substep
                        .transport
                        .endpoint
                        .kirchhoff_source
                        .is_some(),
                    kirchhoff
                );
                assert!(
                    report.maximum_carrier_mass_residual_kg.unwrap()
                        <= time.water_balance_tolerance_kg
                );
                assert!(
                    report.total_mass_residual_kg.unwrap().abs() <= time.water_balance_tolerance_kg
                );
                assert_eq!(air, original_air);
                assert_eq!(original.moisture(), &[0.1]);
                assert_eq!(original.time_s(), 0.0);
            }
        }
    });
}

#[test]
fn g1_shared_slab_sdirk_adaptive_resolves_mean_and_difference_modes() {
    with_cx(|cx| {
        let target = 2.0;
        let slope = 0.25;
        let air = atmosphere::atmosphere(0.8, 1.0);
        let original_air = air.clone();
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.05], cx);
        let b = slab(card([0.01; 2], 0.3), 1.0, vec![0.20], cx)
            .with_kirchhoff_diffusivity(nonlinear())
            .unwrap();
        let faces = [
            Some(atmosphere::sorption(slope)),
            Some(atmosphere::sorption(slope)),
        ];
        let fields = [participant(8, &b, &faces), participant(3, &a, &faces)];
        let rh = air.relative_humidity().unwrap();
        let capacity = air.vapor_mass_kg().unwrap() / rh;
        let dry = a.dry_mass_kg() + b.dry_mass_kg();
        let equilibrium = slope * (capacity * rh + dry * 0.125) / (capacity + slope * dry);
        // The mean exchanges water with ONE gas store; equal/opposite field
        // deviations form a separate mode with zero net reservoir debit.
        let mean = equilibrium
            + (0.125 - equilibrium) * (-0.04 * (1.0 + slope * dry / capacity) * target).exp();
        let difference = 0.075 * (-0.04 * target).exp();
        let expected = [mean - difference, mean + difference];
        let loose = air
            .advance_shared_sdirk2_adaptive_to(
                &fields,
                target,
                2,
                nonlinear(),
                atmosphere::controls(),
                time_controls(target, 2e-4),
                cx,
            )
            .unwrap();
        let time = time_controls(target, 2e-6);
        let fine = air
            .advance_shared_sdirk2_adaptive_to(
                &fields,
                target,
                2,
                nonlinear(),
                atmosphere::controls(),
                time,
                cx,
            )
            .unwrap();
        let error = |candidate: &SharedAtmosphereMoistureAdvanceCandidate<MoistureSdirk2Report>| {
            candidate
                .solids()
                .iter()
                .zip(expected)
                .map(|(field, expected)| (field.moisture.moisture()[0] - expected).abs())
                .fold(0.0_f64, f64::max)
        };
        assert!(error(&fine) < 0.2 * error(&loose));
        assert!(error(&fine) < 1e-6);
        near(
            fine.reservoir().relative_humidity().unwrap(),
            rh - dry * (mean - 0.125) / capacity,
            2e-5,
        );
        let report = fine.report();
        assert!(report.summary.rejected_intervals > 0);
        assert!(report.summary.evaluations > loose.report().summary.evaluations);
        assert_eq!(report.summary.evaluations, 6 * report.summary.comparisons);
        assert_eq!(report.summary.interval_s, [0.0, target]);
        assert_eq!(
            report.specimens.iter().map(|s| s.id).collect::<Vec<_>>(),
            vec![3, 8]
        );
        for (i, (field, account)) in fine.solids().iter().zip(&report.specimens).enumerate() {
            assert_eq!(field.id, account.id);
            assert_eq!(field.moisture.time_s(), target);
            assert_eq!(field.transport.endpoint.interval_s[1], target);
            assert_eq!(account.water_change_kg.is_sign_positive(), i == 0);
            let water =
                field.moisture.dry_mass_kg() * (field.moisture.moisture()[0] - [0.05, 0.20][i]);
            near(account.water_change_kg, water, 1e-13);
            near(
                account.boundary_water_kg.iter().sum(),
                water,
                time.water_balance_tolerance_kg,
            );
            assert_eq!(field.transport.endpoint.kirchhoff_source.is_some(), i == 1);
        }
        near(
            fine.reservoir().vapor_mass_kg().unwrap()
                + fine
                    .solids()
                    .iter()
                    .map(|s| s.moisture.water_mass_kg())
                    .sum::<f64>(),
            air.vapor_mass_kg().unwrap() + a.water_mass_kg() + b.water_mass_kg(),
            time.water_balance_tolerance_kg,
        );
        near(
            report.summary.water_balance_residual_kg,
            0.0,
            time.water_balance_tolerance_kg,
        );
        let reversed = air
            .advance_shared_sdirk2_adaptive_to(
                &[participant(3, &a, &faces), participant(8, &b, &faces)],
                target,
                2,
                nonlinear(),
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
        for ((left, right), (left_account, right_account)) in fine
            .solids()
            .iter()
            .zip(reversed.solids())
            .zip(report.specimens.iter().zip(&reversed.report().specimens))
        {
            assert_eq!(left.moisture.moisture(), right.moisture.moisture());
            assert_eq!(
                left_account.boundary_water_kg,
                right_account.boundary_water_kg
            );
        }
        assert_eq!(air, original_air);
        assert_eq!(a.moisture(), &[0.05]);
        assert_eq!(b.moisture(), &[0.20]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}

#[test]
fn g1_shared_slab_sdirk_adaptive_keeps_local_error_when_vapor_is_stationary() {
    with_cx(|cx| {
        let target = 10.0;
        let air = atmosphere::atmosphere(0.8, 1.0);
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let b = slab(card([0.01; 2], 0.3), 1.0, vec![0.3], cx);
        let faces = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
        ];
        let coarse = 0.1 * decay(0.04 * target);
        let fine = 0.1 * decay(0.02 * target).powi(2);
        let mut time = time_controls(target, (coarse - fine).abs() / (3.0 * 0.2));
        time.activity_absolute_tolerance = 1e-8;
        time.reservoir_water_absolute_tolerance_kg = 1e-10;
        time.maximum_step_s = 2.0 * target;
        let accepted = air
            .advance_shared_sdirk2_adaptive_to(
                &[participant(7, &a, &faces), participant(2, &b, &faces)],
                target,
                2,
                nonlinear(),
                atmosphere::controls(),
                time,
                cx,
            )
            .unwrap();
        let report = accepted.report();
        assert_eq!(report.summary.evaluations, 6);
        assert_eq!(report.summary.accepted_intervals, 1);
        assert_eq!(report.summary.accepted_substeps, 2);
        assert_eq!(report.summary.rejected_intervals, 0);
        near(report.summary.maximum_accepted_error_ratio, 0.2, 2e-7);
        // A cubic growth rule doubles at ratio <= 1/8; this comparison must
        // retain its span despite a stationary gas and zero net uptake.
        assert_eq!(report.summary.suggested_step_s, target);
        assert_eq!(accepted.reservoir(), &air);
        near(
            report.summary.water_change_kg,
            0.0,
            time.water_balance_tolerance_kg,
        );
        for (i, (field, account)) in accepted.solids().iter().zip(&report.specimens).enumerate() {
            let sign = if i == 0 { 1.0 } else { -1.0 };
            near(field.moisture.moisture()[0], 0.2 + sign * fine, 2e-12);
            near(account.maximum_accepted_error_ratio, 0.2, 2e-7);
            let water = -sign * field.moisture.dry_mass_kg() * (0.1 - fine);
            assert!(water.abs() > 0.01);
            near(account.water_change_kg, water, 2e-12);
            for face in account.boundary_water_kg {
                near(face, water / 2.0, 2e-12);
            }
        }
        assert_eq!(a.moisture(), &[0.1]);
        assert_eq!(b.moisture(), &[0.3]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}

#[test]
fn g4_slab_sdirk_adaptive_reserves_stages_and_discards_complete_private_windows() {
    with_cx(|cx| {
        let air = atmosphere::atmosphere(0.8, 1.0);
        let original_air = air.clone();
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.05], cx);
        let b = slab(card([0.01; 2], 0.3), 1.0, vec![0.20], cx);
        let faces = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
        ];
        let fields = [participant(8, &b, &faces), participant(3, &a, &faces)];
        let mut time = time_controls(1.0, 1.0);
        for (budget, evaluated, accepted, failed_at) in [(5, 4, 0, 0.5), (11, 10, 1, 1.5)] {
            time.maximum_evaluations = budget;
            let shared = air
                .advance_shared_sdirk2_adaptive_to(
                    &fields,
                    2.0,
                    2,
                    nonlinear(),
                    atmosphere::controls(),
                    time,
                    cx,
                )
                .map(|_| ());
            let single = air
                .advance_sdirk2_adaptive_to(
                    &a,
                    2.0,
                    &faces,
                    config(),
                    nonlinear(),
                    atmosphere::controls(),
                    time,
                    cx,
                )
                .map(|_| ());
            for refusal in [shared, single] {
                assert!(matches!(refusal,
                    Err(MoistureAdvanceError::WorkLimit { evaluations, accepted_intervals, time_s })
                        if evaluations == evaluated && accepted_intervals == accepted && time_s == failed_at
                ));
            }
        }
        time.maximum_evaluations = 12;
        let complete = air
            .advance_shared_sdirk2_adaptive_to(
                &fields,
                2.0,
                2,
                nonlinear(),
                atmosphere::controls(),
                time,
                cx,
            )
            .unwrap();
        assert_eq!(complete.report().summary.evaluations, 12);
        assert_eq!(complete.report().summary.accepted_intervals, 2);
        assert_eq!(complete.report().summary.accepted_substeps, 4);
        assert_eq!(complete.report().summary.rejected_intervals, 0);
        assert!(complete.solids().iter().all(|s| s.moisture.time_s() == 2.0));
        assert!(matches!(
            air.advance_shared_sdirk2_adaptive_to(
                &fields,
                2.0,
                1,
                nonlinear(),
                atmosphere::controls(),
                time,
                cx,
            ),
            Err(MoistureAdvanceError::Invalid(_))
        ));
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
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
            assert!(matches!(
                air.advance_shared_sdirk2_adaptive_to(
                    &fields,
                    2.0,
                    2,
                    nonlinear(),
                    atmosphere::controls(),
                    time,
                    &cancelled,
                ),
                Err(MoistureAdvanceError::Cancelled)
            ));
            assert!(matches!(
                air.advance_sdirk2_adaptive_to(
                    &a,
                    2.0,
                    &faces,
                    config(),
                    nonlinear(),
                    atmosphere::controls(),
                    time,
                    &cancelled,
                ),
                Err(MoistureAdvanceError::Cancelled)
            ));
        });
        assert_eq!(air, original_air);
        assert_eq!(a.moisture(), &[0.05]);
        assert_eq!(b.moisture(), &[0.20]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}
