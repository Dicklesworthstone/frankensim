//! G1/G4 adaptive native slabs retain finite solution mass and local field errors.
//! Thermodynamic and transport coefficients are authored numerical fixtures.
use super::*;
use fs_couple::moisture::{
    atmosphere::DryBasisSorption,
    solution::{FiniteMoistureSolution, MoistureParticipant},
};
use fs_material::osmosis::{AqueousSolution, nonideal::VirialSolutionModel};
use std::sync::Arc;

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

fn rk4(initial: f64, duration: f64, steps: usize, rate: impl Fn(f64) -> f64) -> f64 {
    let h = duration / steps as f64;
    let mut value = initial;
    for _ in 0..steps {
        let a = rate(value);
        let b = rate(value + 0.5 * h * a);
        let c = rate(value + 0.5 * h * b);
        let d = rate(value + h * c);
        value += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    value
}

#[test]
fn g1_g4_solution_slab_sdirk_adaptive_controls_pure_solvent_mass_and_private_work() {
    with_cx(|cx| {
        let volume = 1e-4;
        let bath = solution::reservoir(volume, Vec::new());
        let original_bath = bath.solution().clone();
        let original = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let target = 10.0;
        let dry = original.dry_mass_kg();
        let coarse = 0.2 - 0.1 * decay(0.04 * target);
        let fine = 0.2 - 0.1 * decay(0.02 * target).powi(2);
        let exact = 0.2 - 0.1 * (-0.04 * target).exp();
        let mut time = time_controls(target, 1.0);
        time.activity_absolute_tolerance = 1e-15;
        time.reservoir_water_absolute_tolerance_kg = dry * (coarse - fine).abs() / (3.0 * 0.2);
        time.maximum_step_s = 2.0 * target;
        let (next_bath, next, report) = bath
            .advance_sdirk2_adaptive_to(
                &original,
                target,
                &faces,
                config(),
                nonlinear(),
                solution::controls(volume),
                time,
                cx,
            )
            .unwrap();
        // Activity is exactly one and the moisture tolerance is deliberately
        // loose. The independent solvent-mass error must set this ratio.
        near(report.maximum_accepted_error_ratio, 0.2, 2e-7);
        assert_eq!(report.evaluations, 6);
        assert_eq!(report.accepted_intervals, 1);
        assert_eq!(report.accepted_substeps, 2);
        assert_eq!(report.rejected_intervals, 0);
        assert_eq!(report.suggested_step_s, target);
        assert_eq!(report.activity, Some([1.0; 2]));
        near(next.moisture()[0], fine, 3e-13);
        let water = dry * (fine - 0.1);
        near(
            next_bath.solvent_mass_kg(),
            bath.solvent_mass_kg() - water,
            2e-13,
        );
        near(report.water_change_kg, water, 2e-13);
        for face in report.boundary_water_kg {
            near(face, water / 2.0, 2e-13);
        }
        near(
            next_bath.solvent_mass_kg() + next.water_mass_kg(),
            bath.solvent_mass_kg() + original.water_mass_kg(),
            time.water_balance_tolerance_kg,
        );
        assert_eq!(report.last_substep.transport.endpoint.interval_s[1], target);

        time.reservoir_water_absolute_tolerance_kg *= 0.01;
        let (_, tight, tight_report) = bath
            .advance_sdirk2_adaptive_to(
                &original,
                target,
                &faces,
                config(),
                nonlinear(),
                solution::controls(volume),
                time,
                cx,
            )
            .unwrap();
        assert!((tight.moisture()[0] - exact).abs() < 0.2 * (fine - exact).abs());
        assert!(tight_report.rejected_intervals > 0);
        assert!(tight_report.evaluations > report.evaluations);
        assert_eq!(tight.time_s(), target);

        let mut limited = time_controls(1.0, 1.0);
        for (budget, evaluated, accepted, failed_at) in [(5, 4, 0, 0.5), (11, 10, 1, 1.5)] {
            limited.maximum_evaluations = budget;
            let single = bath
                .advance_sdirk2_adaptive_to(
                    &original,
                    2.0,
                    &faces,
                    config(),
                    nonlinear(),
                    solution::controls(volume),
                    limited,
                    cx,
                )
                .map(|_| ());
            let shared = bath
                .advance_shared_sdirk2_adaptive_to(
                    &[participant(7, &original, &faces)],
                    2.0,
                    1,
                    nonlinear(),
                    solution::controls(volume),
                    limited,
                    cx,
                )
                .map(|_| ());
            for result in [single, shared] {
                assert!(matches!(result, Err(MoistureAdvanceError::WorkLimit {
                    evaluations, accepted_intervals, time_s,
                }) if evaluations == evaluated && accepted_intervals == accepted && time_s == failed_at));
            }
        }
        assert_eq!(bath.solution(), &original_bath);
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn g1_shared_solution_slab_sdirk_adaptive_refines_original_nonideal_mean_and_difference_modes() {
    with_cx(|cx| {
        let volume = 1e-4;
        let amount = 0.1;
        let target = 10.0;
        let slope = 0.25;
        let model = Arc::new(
            VirialSolutionModel::new(
                293.15,
                solution::WATER_VOLUME,
                1.0,
                vec![3000.0],
                vec![0.002],
            )
            .unwrap(),
        );
        let bath = FiniteMoistureSolution::new(
            AqueousSolution::nonideal(model.clone(), volume, vec![amount]).unwrap(),
            0.018,
        )
        .unwrap();
        let before_bath = bath.solution().clone();
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.05], cx);
        let b = slab(card([0.01; 2], 0.3), 1.0, vec![0.28], cx)
            .with_kirchhoff_diffusivity(nonlinear())
            .unwrap();
        let faces = [
            Some(atmosphere::sorption(slope)),
            Some(atmosphere::sorption(slope)),
        ];
        let fields = [participant(8, &b, &faces), participant(3, &a, &faces)];
        let dry = a.dry_mass_kg() + b.dry_mass_kg();
        let density = bath.solvent_mass_kg() / volume;
        let mean0 = 0.165;
        let activity = |v: f64| {
            let c = amount / v;
            (-solution::WATER_VOLUME * (c + 0.001 * c * c)).exp()
        };
        // Eliminating the ONE solvent store leaves a nonlinear scalar mean
        // equation; the independent difference still decays with rate 4D/L².
        let rate = |mean: f64| {
            let v = volume - dry * (mean - mean0) / density;
            0.04 * (slope * activity(v) - mean)
        };
        let mean = rk4(mean0, target, 2048, rate);
        near(mean, rk4(mean0, target, 4096, rate), 5e-13);
        let difference = 0.115 * (-0.04 * target).exp();
        let expected = [mean - difference, mean + difference];
        let expected_volume = volume - dry * (mean - mean0) / density;
        let loose = bath
            .advance_shared_sdirk2_adaptive_to(
                &fields,
                target,
                2,
                nonlinear(),
                solution::controls(volume),
                time_controls(target, 2e-4),
                cx,
            )
            .unwrap();
        let time = time_controls(target, 2e-7);
        let tight = bath
            .advance_shared_sdirk2_adaptive_to(
                &fields,
                target,
                2,
                nonlinear(),
                solution::controls(volume),
                time,
                cx,
            )
            .unwrap();
        let error = |candidate: &SharedSolutionMoistureAdvanceCandidate<MoistureSdirk2Report>| {
            candidate
                .solids()
                .iter()
                .zip(expected)
                .map(|(field, expected)| (field.moisture.moisture()[0] - expected).abs())
                .fold(0.0_f64, f64::max)
        };
        assert!(error(&tight) < 0.2 * error(&loose));
        assert!(error(&tight) < 1e-6);
        near(
            tight.reservoir().solution().solvent_volume_m3(),
            expected_volume,
            2e-9,
        );
        near(
            tight.reservoir().solution().water_activity().unwrap(),
            activity(expected_volume),
            2e-6,
        );
        assert_eq!(tight.reservoir().solution().solute_moles(), &[amount]);
        assert!(std::ptr::eq(
            tight.reservoir().solution().nonideal_model().unwrap(),
            model.as_ref()
        ));
        let report = tight.report();
        assert!(report.summary.rejected_intervals > 0);
        assert!(report.summary.evaluations > loose.report().summary.evaluations);
        assert_eq!(report.summary.evaluations, 6 * report.summary.comparisons);
        assert_eq!(
            report.summary.accepted_substeps,
            2 * report.summary.accepted_intervals
        );
        assert_eq!(report.summary.interval_s, [0.0, target]);
        assert_eq!(
            report.specimens.iter().map(|s| s.id).collect::<Vec<_>>(),
            vec![3, 8]
        );
        for (i, (field, account)) in tight.solids().iter().zip(&report.specimens).enumerate() {
            assert_eq!(field.id, account.id);
            assert_eq!(field.moisture.time_s(), target);
            assert_eq!(field.transport.endpoint.interval_s[1], target);
            assert_eq!(account.water_change_kg.is_sign_positive(), i == 0);
            let water =
                field.moisture.dry_mass_kg() * (field.moisture.moisture()[0] - [0.05, 0.28][i]);
            near(account.water_change_kg, water, 2e-13);
            near(
                account.boundary_water_kg.iter().sum(),
                water,
                time.water_balance_tolerance_kg,
            );
            assert_eq!(field.transport.endpoint.kirchhoff_source.is_some(), i == 1);
            for source in field.transport.endpoint.surface_sorption.iter().flatten() {
                near(
                    source.query().axes()["relative_humidity"],
                    tight.reservoir().solution().water_activity().unwrap(),
                    2e-12,
                );
            }
        }
        near(
            tight.reservoir().solvent_mass_kg()
                + tight
                    .solids()
                    .iter()
                    .map(|s| s.moisture.water_mass_kg())
                    .sum::<f64>(),
            bath.solvent_mass_kg() + a.water_mass_kg() + b.water_mass_kg(),
            time.water_balance_tolerance_kg,
        );
        near(
            report.summary.water_balance_residual_kg,
            0.0,
            time.water_balance_tolerance_kg,
        );
        let reversed = bath
            .advance_shared_sdirk2_adaptive_to(
                &[participant(3, &a, &faces), participant(8, &b, &faces)],
                target,
                2,
                nonlinear(),
                solution::controls(volume),
                time,
                cx,
            )
            .unwrap();
        assert_eq!(
            tight.reservoir().solution(),
            reversed.reservoir().solution()
        );
        assert_eq!(
            report.summary.evaluations,
            reversed.report().summary.evaluations
        );
        for (left, right) in tight.solids().iter().zip(reversed.solids()) {
            assert_eq!(left.moisture.moisture(), right.moisture.moisture());
        }
        assert_eq!(bath.solution(), &before_bath);
        assert_eq!(a.moisture(), &[0.05]);
        assert_eq!(b.moisture(), &[0.28]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}
