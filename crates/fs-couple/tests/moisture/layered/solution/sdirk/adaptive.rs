//! Whole-duration finite-bath accuracy on the independent material references.
use super::*;
use fs_couple::moisture::adaptive::{
    LayeredMoistureAdvanceReport, MoistureAdvanceError, MoistureTimeControls,
};

type Report = LayeredMoistureAdvanceReport<SolutionMoistureStep<MoistureLayeredSdirk2Report>>;

fn time(maximum: f64, tolerance: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: maximum,
        minimum_step_s: 1e-7,
        maximum_step_s: maximum,
        maximum_evaluations: 5000,
        moisture_absolute_tolerance: tolerance,
        activity_absolute_tolerance: tolerance,
        reservoir_water_absolute_tolerance_kg: 1.0,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-10,
    }
}

fn account(report: &Report, p: &Physics, initial: [f64; 2], next: &[f64], duration: f64) {
    assert_eq!(report.summary.interval_s, [0.0, duration]);
    assert_eq!(report.layers.len(), 2);
    assert_eq!(report.summary.boundary_water_kg[0], 0.0);
    assert_eq!(
        report.layers[0].boundary_water_kg[1],
        -report.layers[1].boundary_water_kg[0]
    );
    assert_eq!(
        report.summary.accepted_substeps,
        2 * report.summary.accepted_intervals
    );
    assert!(report.summary.maximum_accepted_error_ratio <= 1.0);
    assert_eq!(report.summary.evaluations % 2, 0);
    for i in 0..2 {
        let change = p.capacity[i] * (next[i] - initial[i]);
        near(report.layers[i].water_change_kg, change, 1e-10);
        near(
            report.layers[i].boundary_water_kg.iter().sum(),
            change,
            1e-10,
        );
    }
    near(report.summary.water_balance_residual_kg, 0.0, 1e-10);
    near(report.summary.transport_balance_residual_kg, 0.0, 1e-10);
    near(
        report.summary.water_change_kg,
        report.summary.boundary_water_kg[1],
        1e-10,
    );
}

#[test]
fn layered_finite_bath_accuracy_refines_sdirk_and_each_material_account() {
    with_cx(|cx| {
        let initial = [0.1, 0.2];
        let volume = 1e-4;
        let solute = 0.2;
        let duration = 8.0;
        for symmetry in [
            None,
            Some(RadialMoistureSymmetry::Sphere),
            Some(RadialMoistureSymmetry::Cylinder { length_m: 0.1 }),
        ] {
            let p = coefficients(symmetry);
            let exact = continuous(&p, initial, duration, volume, solute, 32768);
            let mut errors = Vec::new();
            let mut evaluations = Vec::new();
            for tolerance in [1e-3, 1e-6] {
                let bath = aqueous::reservoir(volume, vec![solute]);
                let (next_bath, activity, final_water, report) = if let Some(symmetry) = symmetry {
                    let old = radial_body(symmetry, initial, cx);
                    let (next_bath, next, report) = bath
                        .advance_layered_radial_sdirk2_adaptive_to(
                            &old,
                            duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            aqueous::controls(volume),
                            time(duration, tolerance),
                            cx,
                        )
                        .unwrap();
                    assert_eq!(old.time_s(), 0.0);
                    assert_eq!(next.time_s(), duration);
                    (
                        next_bath,
                        next.activity().to_vec(),
                        next.water_mass_kg(),
                        report,
                    )
                } else {
                    let old = slab_body(initial, cx);
                    let (next_bath, next, report) = bath
                        .advance_layered_sdirk2_adaptive_to(
                            &old,
                            duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            aqueous::controls(volume),
                            time(duration, tolerance),
                            cx,
                        )
                        .unwrap();
                    assert_eq!(old.time_s(), 0.0);
                    assert_eq!(next.time_s(), duration);
                    (
                        next_bath,
                        next.activity().to_vec(),
                        next.water_mass_kg(),
                        report,
                    )
                };
                account(&report, &p, initial, &activity, duration);
                let original_water: f64 = (0..2).map(|i| p.capacity[i] * initial[i]).sum();
                near(
                    next_bath.solvent_mass_kg() + final_water,
                    bath.solvent_mass_kg() + original_water,
                    1e-10,
                );
                assert_eq!(next_bath.solution().solute_moles(), &[solute]);
                assert_eq!(
                    report.summary.reservoir_water_kg.unwrap(),
                    [bath.solvent_mass_kg(), next_bath.solvent_mass_kg()]
                );
                errors.push(
                    (0..2)
                        .map(|i| (activity[i] - exact[i]).abs())
                        .fold(0.0_f64, f64::max),
                );
                evaluations.push(report.summary.evaluations);
            }
            assert!(errors[1] < 0.3 * errors[0], "{symmetry:?}: {errors:?}");
            assert!(
                evaluations[1] > evaluations[0],
                "{symmetry:?}: {evaluations:?}"
            );
        }
    });
}

#[test]
fn loose_layered_solution_accuracy_recovers_native_fine_path_without_stage_double_debit() {
    with_cx(|cx| {
        let initial = [0.1, 0.2];
        let volume = 1e-4;
        let duration = 2.0;
        for symmetry in [
            None,
            Some(RadialMoistureSymmetry::Sphere),
            Some(RadialMoistureSymmetry::Cylinder { length_m: 0.1 }),
        ] {
            let bath = aqueous::reservoir(volume, vec![0.2]);
            let (actual_bath, actual, report, expected_bath, expected, receipts) =
                if let Some(symmetry) = symmetry {
                    let old = radial_body(symmetry, initial, cx);
                    let (half_bath, half, first) = bath
                        .advance_layered_radial_sdirk2(
                            &old,
                            0.5 * duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            aqueous::controls(volume),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    let (expected_bath, expected, second) = half_bath
                        .advance_layered_radial_sdirk2(
                            &half,
                            0.5 * duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            aqueous::controls(volume),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    let (actual_bath, actual, report) = bath
                        .advance_layered_radial_sdirk2_adaptive_to(
                            &old,
                            duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            aqueous::controls(volume),
                            time(duration, 1.0),
                            cx,
                        )
                        .unwrap();
                    (
                        actual_bath,
                        actual.activity().to_vec(),
                        report,
                        expected_bath,
                        expected.activity().to_vec(),
                        [first, second],
                    )
                } else {
                    let old = slab_body(initial, cx);
                    let (half_bath, half, first) = bath
                        .advance_layered_sdirk2(
                            &old,
                            0.5 * duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            aqueous::controls(volume),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    let (expected_bath, expected, second) = half_bath
                        .advance_layered_sdirk2(
                            &half,
                            0.5 * duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            aqueous::controls(volume),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    let (actual_bath, actual, report) = bath
                        .advance_layered_sdirk2_adaptive_to(
                            &old,
                            duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            aqueous::controls(volume),
                            time(duration, 1.0),
                            cx,
                        )
                        .unwrap();
                    (
                        actual_bath,
                        actual.activity().to_vec(),
                        report,
                        expected_bath,
                        expected.activity().to_vec(),
                        [first, second],
                    )
                };
            assert_eq!(actual, expected);
            assert_eq!(actual_bath.solution(), expected_bath.solution());
            assert_eq!(report.summary.evaluations, 6);
            assert_eq!(report.summary.accepted_intervals, 1);
            let exterior: f64 = receipts
                .iter()
                .map(|r| r.transport.endpoint.boundary_water_kg[1])
                .sum();
            let interface: f64 = receipts
                .iter()
                .map(|r| r.transport.endpoint.interfaces[0].water_left_to_right_kg)
                .sum();
            near(report.summary.boundary_water_kg[1], exterior, 2e-14);
            near(report.layers[0].boundary_water_kg[1], -interface, 2e-14);
            near(report.layers[1].boundary_water_kg[0], interface, 2e-14);
            account(&report, &coefficients(symmetry), initial, &actual, duration);
        }
    });
}

#[test]
fn layered_sdirk_solution_work_exhaustion_discards_an_accepted_private_prefix() {
    with_cx(|cx| {
        let volume = 1e-4;
        let bath = aqueous::reservoir(volume, vec![0.2]);
        let saved = bath.solution().clone();
        let step = 1.0 / 64.0;
        let mut limited = time(step, 1.0);
        limited.maximum_evaluations = 6;
        let slab = slab_body([0.1, 0.2], cx);
        assert!(matches!(bath.advance_layered_sdirk2_adaptive_to(
            &slab,2.0*step,[false,true],config(),nonlinear(),aqueous::controls(volume),limited,cx,
        ),Err(MoistureAdvanceError::WorkLimit{time_s,accepted_intervals:1,evaluations:6}) if time_s==step));
        assert_eq!(slab.time_s(), 0.0);
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let radial = radial_body(symmetry, [0.1, 0.2], cx);
            let initial = radial.activity().to_vec();
            assert!(matches!(bath.advance_layered_radial_sdirk2_adaptive_to(
                &radial,2.0*step,[false,true],config(),nonlinear(),aqueous::controls(volume),limited,cx,
            ),Err(MoistureAdvanceError::WorkLimit{time_s,accepted_intervals:1,evaluations:6}) if time_s==step));
            assert_eq!(radial.time_s(), 0.0);
            assert_eq!(radial.activity(), initial);
        }
        assert_eq!(bath.solution(), &saved);
    });
}
