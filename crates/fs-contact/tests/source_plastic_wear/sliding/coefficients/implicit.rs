//! G1/G3 implicit endpoint wear against independent one-cell equations.

use super::*;
use fs_contact::normal_patch::plastic_source::wear::implicit::{
    ImplicitWearBoundary, ImplicitWearBoundaryReceipt, ImplicitWearControls, ImplicitWearError,
};
use fs_solver::NewtonKrylovConfig;

fn solve_controls() -> ImplicitWearControls {
    ImplicitWearControls {
        maximum_unknowns: 25,
        maximum_evaluations: 300,
        maximum_newton_iterations: 24,
        difference_step: 1e-6,
        pressure_tolerance_pa: 1e-3,
        recession_tolerance_m: 1e-13,
        relative_tolerance: 1e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 4,
            ..NewtonKrylovConfig::default()
        },
    }
}

fn slip(distance: f64) -> Vec<f64> {
    let mut values = vec![0.0; 25];
    values[12] = distance;
    values
}

fn self_compliance(source: &SourcePlasticArchardWear) -> f64 {
    // Independent square-cell integral of the half-space 1/(pi E* r) kernel.
    let width = source.contact().model().config().grid.cell_width_m;
    let modulus = source
        .contact()
        .normal_source()
        .elastic()
        .reduced_modulus_pa();
    4.0 * width * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * modulus)
}

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

#[test]
fn implicit_large_step_preserves_a_positive_contact_and_the_original_dent() {
    with_cx(|cx, _| {
        let source = fixture([0.01, 0.03], [1e6, 2e6], [1e-4; 2], false, cx);
        let c = self_compliance(&source);
        let dented = source
            .advance_normal(&source.initial_state(), 5e-5, cx)
            .unwrap()
            .candidate;
        let dent = dented.contact_state().contact_state().plastic_depth_m()[12];
        let p0 = 0.5e6;
        let before = source
            .advance_normal(&dented, dent + c * p0, cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let distance = 1e-3;
        let beta = (0.01 / 1e6 + 0.03 / 2e6) * distance;
        assert!(beta > c);
        let explicit = source
            .advance_wear_during(&before, &slip(distance), 1.0, cx)
            .unwrap();
        assert_eq!(
            explicit
                .candidate
                .contact_state()
                .contact_state()
                .pressure_pa()[12],
            0.0
        );
        let result = source
            .advance_wear_implicit(
                &before,
                &slip(distance),
                1.0,
                ImplicitWearBoundary::HeldApproach,
                before.contact_state().contact_state().pressure_pa(),
                solve_controls(),
                cx,
            )
            .unwrap();
        let next = result.wear.candidate.contact_state().contact_state();
        let expected = p0 / (1.0 + beta / c);
        near(next.pressure_pa()[12], expected, 1e-8, 1e-3);
        near(result.quadrature_pressure_pa[12], expected, 1e-8, 1e-3);
        for side in 0..2 {
            let depth = [0.01 / 1e6, 0.03 / 2e6][side] * expected * distance;
            near(
                result.wear.recession_increment_m[side][12],
                depth,
                1e-8,
                1e-13,
            );
            near(
                result.wear.wear[side].increment_m3,
                depth * 1e-8,
                1e-8,
                1e-21,
            );
            near(
                result.wear.candidate.volumes_m3()[side],
                depth * 1e-8,
                1e-8,
                1e-21,
            );
        }
        assert_eq!(
            next.plastic_depth_m(),
            before.contact_state().contact_state().plastic_depth_m()
        );
        assert_eq!(
            next.cumulative_plastic_work_j(),
            before
                .contact_state()
                .contact_state()
                .cumulative_plastic_work_j()
        );
        assert_eq!(result.wear.receipt.response.plastic_work_j, 0.0);
        near(
            result.wear.receipt.elastic_energy_change_j,
            0.5 * 1e-8 * c * (expected * expected - p0 * p0),
            1e-8,
            1e-16,
        );
        assert_eq!(result.wear.receipt.normal_work_j, 0.0);
        near(
            result.wear.receipt.parameter_work_j,
            result.wear.receipt.elastic_energy_change_j,
            1e-12,
            1e-18,
        );
        assert!(result.maximum_normalized_residual <= 1.0);
        assert!(result.nonlinear.as_ref().unwrap().converged);
        assert!(result.evaluations > 2);
        assert_eq!(before, saved);
        let replay = source
            .advance_wear_implicit(
                &before,
                &slip(distance),
                1.0,
                ImplicitWearBoundary::HeldApproach,
                before.contact_state().contact_state().pressure_pa(),
                solve_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result, replay);

        // Repeated native intervals retain both inventories and converge toward
        // the independently integrated constant-coefficient pressure decay.
        let mut refined = before.clone();
        for _ in 0..8 {
            refined = source
                .advance_wear_implicit(
                    &refined,
                    &slip(distance / 8.0),
                    1.0 / 8.0,
                    ImplicitWearBoundary::HeldApproach,
                    refined.contact_state().contact_state().pressure_pa(),
                    solve_controls(),
                    cx,
                )
                .unwrap()
                .wear
                .candidate;
        }
        let p8 = refined.contact_state().contact_state().pressure_pa()[12];
        near(p8, p0 / (1.0 + beta / (8.0 * c)).powi(8), 1e-8, 1e-3);
        let continuous = p0 * (-beta / c).exp();
        assert!((p8 - continuous).abs() < 0.5 * (expected - continuous).abs());
        let removed = refined.recession_m()[0][12] + refined.recession_m()[1][12];
        near(removed, c * (p0 - p8), 1e-8, 1e-13);
    });
}

#[test]
fn maintained_force_and_loading_spring_keep_native_work_and_separate_hardness() {
    with_cx(|cx, _| {
        let source = fixture([0.01, 0.03], [1e6, 2e6], [1e-4; 2], false, cx);
        let c = self_compliance(&source);
        let area = 1e-8;
        let p0 = 0.5e6;
        let distance = 5e-4;
        let beta = (0.01 / 1e6 + 0.03 / 2e6) * distance;
        let loaded = source
            .advance_normal_at_load(&source.initial_state(), area * p0, normal_controls(), cx)
            .unwrap();
        let before = loaded.candidate;
        let maintained = source
            .advance_wear_implicit(
                &before,
                &slip(distance),
                0.5,
                ImplicitWearBoundary::MaintainedLoad {
                    force_n: area * p0,
                    controls: normal_controls(),
                },
                before.contact_state().contact_state().pressure_pa(),
                solve_controls(),
                cx,
            )
            .unwrap();
        let next = maintained.wear.candidate.contact_state().contact_state();
        near(next.pressure_pa()[12], p0, 1e-9, 1e-4);
        near(
            next.approach_m() - before.contact_state().contact_state().approach_m(),
            beta * p0,
            1e-9,
            1e-14,
        );
        let work = area * p0 * beta * p0;
        near(maintained.wear.receipt.normal_work_j, work, 1e-9, 1e-17);
        near(maintained.wear.receipt.parameter_work_j, -work, 1e-9, 1e-17);
        near(
            maintained.wear.receipt.elastic_energy_change_j,
            0.0,
            0.0,
            1e-16,
        );
        assert!(matches!(
            maintained.boundary,
            ImplicitWearBoundaryReceipt::MaintainedLoad(_)
        ));

        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: area / c,
            position_m: 2.0 * c * p0,
        };
        let loaded = source
            .advance_normal_with_drive(&source.initial_state(), drive, normal_controls(), cx)
            .unwrap();
        let spring_energy = loaded.drive.spring_energy_j;
        let before = loaded.candidate;
        let driven = source
            .advance_wear_implicit(
                &before,
                &slip(distance),
                0.5,
                ImplicitWearBoundary::FixedDrive {
                    drive,
                    controls: normal_controls(),
                },
                before.contact_state().contact_state().pressure_pa(),
                solve_controls(),
                cx,
            )
            .unwrap();
        let expected = p0 / (1.0 + beta / (c + area / drive.stiffness_n_per_m));
        let next = driven.wear.candidate.contact_state().contact_state();
        near(next.pressure_pa()[12], expected, 1e-8, 1e-3);
        near(
            next.approach_m(),
            drive.position_m - area * expected / drive.stiffness_n_per_m,
            1e-9,
            1e-13,
        );
        let ImplicitWearBoundaryReceipt::FixedDrive(receipt) = &driven.boundary else {
            panic!("spring receipt");
        };
        let release =
            0.5 * (area * area / drive.stiffness_n_per_m) * (p0 * p0 - expected * expected);
        near(
            spring_energy - receipt.spring_energy_j,
            release,
            1e-8,
            1e-16,
        );
        near(driven.wear.receipt.normal_work_j, release, 1e-8, 1e-16);
        near(
            driven.wear.receipt.elastic_energy_change_j
                + driven.wear.receipt.response.plastic_work_j,
            driven.wear.receipt.normal_work_j + driven.wear.receipt.parameter_work_j,
            1e-12,
            1e-18,
        );
        for side in 0..2 {
            near(
                driven.wear.recession_increment_m[side][12],
                [0.01 / 1e6, 0.03 / 2e6][side] * expected * distance,
                1e-8,
                1e-13,
            );
        }
    });
}

#[test]
fn nonlinear_local_coefficients_follow_the_solved_endpoint_and_actual_speed() {
    with_cx(|cx, _| {
        let sliding = local_source([1e6, 2e6], false, false, cx);
        let source = sliding.wear_source();
        let c = self_compliance(source);
        let p0 = 0.5e6;
        let before = source
            .advance_normal(&source.initial_state(), c * p0, cx)
            .unwrap()
            .candidate;
        let distance = 1e-4;
        let duration = 0.1;
        let a = (0.01 / 1e6 + 0.03 / 2e6) * distance;
        // k_i(p,v)=k_i0*(0.5+p/1e6+v/0.002); here v=0.001.
        // Therefore (a/1e6)*p^2+(c+a)*p-c*p0=0.
        let b = c + a;
        let expected = 2.0 * c * p0 / (b + (b * b + 4.0 * (a / 1e6) * c * p0).sqrt());
        let result = source
            .advance_wear_implicit(
                &before,
                &slip(distance),
                duration,
                ImplicitWearBoundary::HeldApproach,
                before.contact_state().contact_state().pressure_pa(),
                solve_controls(),
                cx,
            )
            .unwrap();
        near(
            result
                .wear
                .candidate
                .contact_state()
                .contact_state()
                .pressure_pa()[12],
            expected,
            1e-8,
            1e-3,
        );
        let samples = result.wear.coefficient_samples.as_ref().unwrap();
        assert_eq!(samples.sources.len(), 1);
        assert_eq!(samples.cell_source_indices[12], Some(0));
        let query = &samples.sources[0];
        assert_eq!(
            query.state_coordinate("wear_p"),
            Some(result.quadrature_pressure_pa[12])
        );
        near(
            query.state_coordinate("wear_v").unwrap(),
            distance / duration,
            1e-15,
            1e-18,
        );
        for side in 0..2 {
            let k = coefficient(side, expected, distance / duration, false);
            let depth = k * expected * distance / [1e6, 2e6][side];
            near(
                result.wear.recession_increment_m[side][12],
                depth,
                1e-8,
                1e-13,
            );
            let key = &source.coefficient_keys()[side];
            assert_eq!(
                query
                    .property_by_key(key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                source
                    .coefficient_source()
                    .property_by_key(key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
        // A loose numerical root must not bypass the separately measured
        // physical wear equation, even if pressure's gate is deliberately lax.
        let mut lax = solve_controls();
        lax.newton.absolute_tolerance = 1.0;
        lax.newton.relative_tolerance = 0.0;
        lax.pressure_tolerance_pa = 1e9;
        let refused = source.advance_wear_implicit(
            &before,
            &slip(distance),
            duration,
            ImplicitWearBoundary::HeldApproach,
            before.contact_state().contact_state().pressure_pa(),
            lax,
            cx,
        );
        assert!(
            matches!(refused, Err(ImplicitWearError::Residual { recession_m, .. }) if recession_m > lax.recession_tolerance_m)
        );
    });
}

#[test]
fn zero_path_and_late_source_budget_or_cancellation_keep_the_accepted_history() {
    with_cx(|cx, gate| {
        let sliding = local_source([1e6, 2e6], false, false, cx);
        let source = sliding.wear_source();
        let before = source
            .advance_normal(&source.initial_state(), self_compliance(source) * 0.5e6, cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let mut zero_controls = solve_controls();
        zero_controls.maximum_unknowns = 0;
        zero_controls.maximum_evaluations = 1;
        let zero = source
            .advance_wear_implicit(
                &before,
                &[0.0; 25],
                0.1,
                ImplicitWearBoundary::HeldApproach,
                before.contact_state().contact_state().pressure_pa(),
                zero_controls,
                cx,
            )
            .unwrap();
        assert!(zero.nonlinear.is_none());
        assert_eq!(zero.evaluations, 1);
        assert_eq!(zero.wear.candidate.volumes_m3(), before.volumes_m3());
        assert_eq!(zero.wear.candidate.recession_m(), before.recession_m());
        assert_eq!(
            zero.wear
                .candidate
                .contact_state()
                .contact_state()
                .pressure_pa(),
            before.contact_state().contact_state().pressure_pa()
        );
        assert!(
            zero.wear
                .coefficient_samples
                .as_ref()
                .unwrap()
                .sources
                .is_empty()
        );
        let mut bounded = solve_controls();
        bounded.maximum_evaluations = 1;
        assert!(matches!(
            source.advance_wear_implicit(
                &before,
                &slip(1e-4),
                0.1,
                ImplicitWearBoundary::HeldApproach,
                before.contact_state().contact_state().pressure_pa(),
                bounded,
                cx,
            ),
            Err(ImplicitWearError::EvaluationBudget)
        ));
        assert!(matches!(
            source.advance_wear_implicit(
                &before,
                &slip(1e-4),
                0.001,
                ImplicitWearBoundary::HeldApproach,
                before.contact_state().contact_state().pressure_pa(),
                solve_controls(),
                cx,
            ),
            Err(ImplicitWearError::Step(SourcePlasticWearError::Material(_)))
        ));
        assert_eq!(before, saved);
        gate.request();
        assert!(
            source
                .advance_wear_implicit(
                    &before,
                    &slip(1e-4),
                    0.1,
                    ImplicitWearBoundary::HeldApproach,
                    before.contact_state().contact_state().pressure_pa(),
                    solve_controls(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
