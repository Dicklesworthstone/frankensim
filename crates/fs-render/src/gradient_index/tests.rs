use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};

fn glass(n: f64, sigma: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::try_constant(sigma).unwrap(),
        GlassProvenance::Custom,
    )
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(1.0, 0.0, 0.0),
    }
}
fn axis() -> Vec3 {
    Vec3::new(0.0, 0.0, 1.0)
}
fn bounds() -> Aabb {
    Aabb::new(Point3::new(-2.0, -0.5, -1.0), Point3::new(2.0, 1.0, 1.0))
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{a} != {b} (tol {tolerance})");
}

#[test]
fn g1_homogeneous_ray_optical_length_and_absorption() {
    let glass = glass(1.5, 0.7);
    let result = propagate(
        &HomogeneousIndexField::new(&glass),
        ray(),
        axis(),
        550.0,
        1.0,
        GradientIndexOptions::default(),
    )
    .unwrap();
    close(result.ray.origin.x, 1.0, 1.0e-14);
    assert_eq!(result.ray.origin.y, 0.0);
    assert_eq!(result.ray.dir, ray().dir);
    assert_eq!(result.transverse_axis, axis());
    close(result.optical_path_length_m, 1.5, 2.0e-14);
    close(result.group_optical_path_length_m.unwrap(), 1.5, 2.0e-14);
    close(result.optical_depth, 0.7, 2.0e-14);
    close(result.power_transmittance(), (-0.7_f64).exp(), 2.0e-14);
    assert_eq!(result.field_evaluations, 1);
    assert_eq!(result.accepted_steps, 1);
}

#[test]
fn g1_affine_index_matches_closed_form_and_snell_invariant() {
    let glass = glass(1.5, 0.7);
    let gradient = 0.8;
    let field = AffineIndexField::try_new(
        &glass,
        ray().origin,
        Vec3::new(0.0, gradient, 0.0),
        bounds(),
    )
    .unwrap();
    let length = 1.3;
    let result = propagate(
        &field,
        ray(),
        axis(),
        550.0,
        length,
        GradientIndexOptions {
            relative_tolerance: 1.0e-12,
            position_tolerance_m: 1.0e-12,
            direction_tolerance: 1.0e-12,
            optical_path_tolerance_m: 1.0e-12,
            ..GradientIndexOptions::default()
        },
    )
    .unwrap();
    let n = (1.5_f64 * 1.5 + (gradient * length).powi(2)).sqrt();
    let inverse_hyperbolic = (gradient * length / 1.5).asinh();
    close(
        result.ray.origin.x,
        1.5 / gradient * inverse_hyperbolic,
        2.0e-11,
    );
    close(result.ray.origin.y, (n - 1.5) / gradient, 2.0e-11);
    close(
        result.final_sample.phase_index() * result.ray.dir.x,
        1.5,
        2.0e-11,
    );
    close(result.ray.dir.y, gradient * length / n, 2.0e-11);
    close(
        result.optical_path_length_m,
        0.5 * (length * n + 1.5 * 1.5 / gradient * inverse_hyperbolic),
        2.0e-11,
    );
    close(result.optical_depth, 0.7 * length, 2.0e-13);
    assert_eq!(result.transverse_axis, axis());
    assert_eq!(
        result.field_evaluations,
        1 + 6 * (result.accepted_steps + result.rejected_steps)
    );
}

#[test]
fn g3_reversal_recovers_position_direction_and_parallel_frame() {
    let glass = glass(1.5, 0.02);
    let field = AffineIndexField::try_new(&glass, ray().origin, Vec3::new(0.0, 0.8, 0.0), bounds())
        .unwrap();
    let options = GradientIndexOptions {
        max_step_m: 0.02,
        ..GradientIndexOptions::default()
    };
    let forward = propagate(&field, ray(), Vec3::new(0.0, 1.0, 1.0), 550.0, 1.0, options).unwrap();
    let backward = propagate(
        &field,
        Ray {
            origin: forward.ray.origin,
            dir: forward.ray.dir.scale(-1.0),
        },
        forward.transverse_axis,
        550.0,
        1.0,
        options,
    )
    .unwrap();
    assert!(norm(backward.ray.origin.delta_from(ray().origin)) < 1.0e-10);
    assert!(norm(sub(backward.ray.dir, ray().dir.scale(-1.0))) < 1.0e-10);
    let original_axis = transverse(ray().dir, Vec3::new(0.0, 1.0, 1.0)).unwrap();
    assert!(norm(sub(backward.transverse_axis, original_axis)) < 1.0e-10);
    close(
        forward.optical_path_length_m,
        backward.optical_path_length_m,
        1.0e-10,
    );
    close(norm(forward.transverse_axis), 1.0, 1.0e-14);
    close(forward.transverse_axis.dot(forward.ray.dir), 0.0, 1.0e-14);
}

#[test]
fn g1_tighter_error_controls_refine_curvature() {
    let glass = glass(1.5, 0.0);
    let field = AffineIndexField::try_new(&glass, ray().origin, Vec3::new(0.0, 0.8, 0.0), bounds())
        .unwrap();
    let loose = GradientIndexOptions {
        position_tolerance_m: 1.0e-3,
        direction_tolerance: 1.0e-3,
        optical_path_tolerance_m: 1.0e-3,
        relative_tolerance: 1.0e-3,
        initial_step_m: 1.2,
        max_step_m: 1.2,
        ..GradientIndexOptions::default()
    };
    let fine = GradientIndexOptions {
        position_tolerance_m: 1.0e-12,
        direction_tolerance: 1.0e-12,
        optical_path_tolerance_m: 1.0e-12,
        relative_tolerance: 1.0e-12,
        ..loose
    };
    let a = propagate(&field, ray(), axis(), 550.0, 1.2, loose).unwrap();
    let b = propagate(&field, ray(), axis(), 550.0, 1.2, fine).unwrap();
    let exact_y = ((1.5_f64 * 1.5 + 0.96 * 0.96).sqrt() - 1.5) / 0.8;
    assert!((b.ray.origin.y - exact_y).abs() < (a.ray.origin.y - exact_y).abs() * 0.01);
    assert!(b.accepted_steps > a.accepted_steps);
    assert!(b.rejected_steps > 0);
    assert!(b.max_normalized_error_estimate <= 1.0);
}

#[test]
fn g1_phase_and_group_index_remain_distinct() {
    let glass = DielectricGlass::new(
        CauchyIor::try_new(1.4, 0.01, 0.0).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    );
    let field = HomogeneousIndexField::new(&glass);
    let result = propagate(
        &field,
        ray(),
        axis(),
        500.0,
        0.4,
        GradientIndexOptions::default(),
    )
    .unwrap();
    close(result.optical_path_length_m, 1.44 * 0.4, 1.0e-14);
    close(
        result.group_optical_path_length_m.unwrap(),
        1.52 * 0.4,
        1.0e-14,
    );
}

#[test]
fn g0_refuses_bad_samples_domains_options_and_frames() {
    assert_eq!(
        GradientIndexSample::try_new(0.8, axis(), 0.0, None),
        Err(GradientIndexError::InvalidSample)
    );
    assert_eq!(
        GradientIndexSample::try_new(1.5, axis(), -1.0, None),
        Err(GradientIndexError::InvalidSample)
    );
    let glass = glass(1.5, 0.0);
    let field = HomogeneousIndexField::new(&glass);
    assert_eq!(
        propagate(
            &field,
            ray(),
            ray().dir,
            550.0,
            0.1,
            GradientIndexOptions::default()
        ),
        Err(GradientIndexError::InvalidTransverseAxis)
    );
    assert_eq!(
        propagate(
            &field,
            ray(),
            axis(),
            550.0,
            -1.0,
            GradientIndexOptions::default()
        ),
        Err(GradientIndexError::InvalidOptions)
    );
    let bounded =
        AffineIndexField::try_new(&glass, ray().origin, Vec3::new(0.0, 0.8, 0.0), bounds())
            .unwrap();
    assert_eq!(
        propagate(
            &bounded,
            ray(),
            axis(),
            550.0,
            1.0,
            GradientIndexOptions {
                max_steps: 1,
                ..GradientIndexOptions::default()
            }
        ),
        Err(GradientIndexError::StepBudgetExceeded)
    );
    assert_eq!(
        bounded.sample(Point3::new(4.0, 0.0, 0.0), 550.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
}

#[test]
fn g4_cancellation_is_polled_even_on_rejected_trials() {
    let glass = glass(1.5, 0.0);
    let field = AffineIndexField::try_new(&glass, ray().origin, Vec3::new(0.0, 0.8, 0.0), bounds())
        .unwrap();
    let mut polls = 0;
    let error = propagate_with_control(
        &field,
        ray(),
        axis(),
        550.0,
        1.0,
        GradientIndexOptions {
            initial_step_m: 1.0,
            max_step_m: 1.0,
            direction_tolerance: 1.0e-15,
            relative_tolerance: 0.0,
            ..GradientIndexOptions::default()
        },
        || {
            polls += 1;
            polls <= 2
        },
    )
    .unwrap_err();
    assert_eq!(error, GradientIndexError::Cancelled);
    assert_eq!(polls, 3);
}

#[test]
fn g0_zero_distance_and_unavailable_group_are_explicit() {
    struct PhaseOnly;
    impl GradientIndexField for PhaseOnly {
        fn sample(&self, _: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            GradientIndexSample::try_new(1.4, Vec3::new(0.0, 0.0, 0.0), 1.0e-15, None)
        }
    }
    let zero = propagate(
        &PhaseOnly,
        ray(),
        axis(),
        550.0,
        0.0,
        GradientIndexOptions::default(),
    )
    .unwrap();
    assert_eq!(zero.ray, ray());
    assert_eq!(zero.accepted_steps, 0);
    assert_eq!(zero.field_evaluations, 1);
    assert_eq!(zero.power_transmittance(), 1.0);
    assert_eq!(zero.group_optical_path_length_m, None);
    let thin = propagate(
        &PhaseOnly,
        ray(),
        axis(),
        550.0,
        0.1,
        GradientIndexOptions::default(),
    )
    .unwrap();
    close(thin.absorbed_power_fraction(), 1.0e-16, 1.0e-29);
}
