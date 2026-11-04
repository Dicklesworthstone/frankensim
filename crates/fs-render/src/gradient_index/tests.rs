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

fn infrared_glass() -> DielectricGlass {
    use crate::dielectric::{DielectricIorSample, TabulatedIor};
    // An exactly linear phase law: n = 1.71875 - lambda_nm / 64000,
    // hence n_g = 1.71875. Absorption has a narrower, independent source band.
    let phase = TabulatedIor::try_new_in_band(
        std::array::from_fn(|i| {
            let wavelength = 8000.0 + 750.0 * i as f64;
            DielectricIorSample::try_new_in_band(wavelength, 1.71875 - wavelength / 64000.0)
                .unwrap()
        }),
        750.0,
    )
    .unwrap();
    DielectricGlass::new(
        phase,
        BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
            (9000.0, 0.2),
            (11000.0, 0.6),
            (13000.0, 1.0),
        ])
        .unwrap(),
        GlassProvenance::Custom,
    )
}

struct Unbounded;
impl event::BoundaryOracle for Unbounded {
    fn query(&mut self, _: Point3) -> Result<event::BoundaryQuery, GradientIndexError> {
        Ok(event::BoundaryQuery {
            clearance_m: f64::INFINITY,
            witness: None,
        })
    }
}

#[test]
fn g1_infrared_scalar_bundle_and_event_paths_keep_phase_group_and_absorption() {
    let glass = infrared_glass();
    let field = HomogeneousIndexField::new(&glass);
    let length = 0.4;
    for wavelength in [9000.0, 11000.0, 13000.0] {
        let n = 1.71875 - wavelength / 64000.0;
        let alpha = 0.2 + (wavelength - 9000.0) / 5000.0;
        let scalar = propagate(
            &field,
            ray(),
            axis(),
            wavelength,
            length,
            GradientIndexOptions::default(),
        )
        .unwrap();
        let bundle = bundle::propagate_bundle(
            &field,
            ray(),
            axis(),
            wavelength,
            length,
            bundle::BundleOptions::default(),
        )
        .unwrap();
        let event = event::propagate_to_boundary(
            &field,
            &mut Unbounded,
            ray(),
            axis(),
            wavelength,
            length,
            None,
            event::GradientEventOptions::default(),
        )
        .unwrap();
        assert_eq!(event.termination, event::BoundaryTermination::DistanceLimit);
        for result in [&scalar, bundle.propagation(), &event.propagation] {
            close(result.ray.origin.x, length, 2e-14);
            assert_eq!(result.ray.dir, ray().dir);
            close(result.optical_path_length_m, n * length, 2e-14);
            close(
                result.group_optical_path_length_m.unwrap(),
                1.71875 * length,
                2e-14,
            );
            close(result.optical_depth, alpha * length, 2e-14);
            close(result.power_transmittance(), (-alpha * length).exp(), 2e-14);
            assert_eq!(result.field_evaluations, 1);
        }
    }
}

#[test]
fn g1_infrared_affine_and_quadratic_fields_use_the_common_source_band() {
    let glass = infrared_glass();
    assert_eq!(glass.wavelength_bounds_nm().unwrap(), (9000.0, 13000.0));
    let field = AffineIndexField::try_new(&glass, ray().origin, Vec3::new(0.0, 0.2, 0.0), bounds())
        .unwrap();
    let length = 0.6;
    let wavelength = 11000.0;
    let n0 = 1.71875 - wavelength / 64000.0;
    let scalar = propagate(
        &field,
        ray(),
        axis(),
        wavelength,
        length,
        GradientIndexOptions::default(),
    )
    .unwrap();
    let n = (n0 * n0 + (0.2_f64 * length).powi(2)).sqrt();
    close(scalar.ray.origin.y, (n - n0) / 0.2, 2e-10);
    close(
        scalar.final_sample.phase_index() * scalar.ray.dir.x,
        n0,
        2e-10,
    );
    close(
        scalar.group_optical_path_length_m.unwrap() - scalar.optical_path_length_m,
        wavelength / 64000.0 * length,
        2e-10,
    );
    let lens = bundle::QuadraticIndexField::try_new(
        &glass,
        ray().origin,
        Vec3::new(0.0, 0.0, 0.0),
        bundle::IndexHessian::try_new([0.0, -n0, -n0], [0.0; 3]).unwrap(),
        Aabb::new(Point3::new(-1.0, -0.1, -0.1), Point3::new(1.0, 0.1, 0.1)),
    )
    .unwrap();
    let focused = bundle::propagate_bundle(
        &lens,
        ray(),
        axis(),
        wavelength,
        length,
        bundle::BundleOptions::default(),
    )
    .unwrap();
    for (column, launch) in focused
        .angular_derivatives()
        .into_iter()
        .zip(focused.launch_basis())
    {
        close(
            norm(sub(column.position_m_per_rad, launch.scale(length.sin()))),
            0.0,
            2e-9,
        );
        close(
            norm(sub(column.direction_per_rad, launch.scale(length.cos()))),
            0.0,
            2e-9,
        );
    }
    for wavelength in [8999.0, 13001.0, 550.0, 0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            propagate(
                &field,
                ray(),
                axis(),
                wavelength,
                length,
                GradientIndexOptions::default(),
            )
            .is_err()
        );
        assert!(
            bundle::propagate_bundle(
                &lens,
                ray(),
                axis(),
                wavelength,
                length,
                bundle::BundleOptions::default(),
            )
            .is_err()
        );
        assert!(
            event::propagate_to_boundary(
                &field,
                &mut Unbounded,
                ray(),
                axis(),
                wavelength,
                length,
                None,
                event::GradientEventOptions::default(),
            )
            .is_err()
        );
    }
    assert_eq!(
        propagate_with_control(
            &field,
            ray(),
            axis(),
            wavelength,
            length,
            GradientIndexOptions::default(),
            || false,
        ),
        Err(GradientIndexError::Cancelled)
    );
}

#[test]
fn g1_infrared_lte_emission_uses_the_same_absorbing_path_and_si_planck_law() {
    let glass = infrared_glass();
    let field = HomogeneousIndexField::new(&glass);
    let temperature = 420.0;
    let wavelength = 11000.0;
    let length = 0.4;
    let (path, emission) = thermal::propagate_emission_to_boundary_with_control(
        &field,
        &|_: Point3| Ok(temperature),
        1e-13,
        &mut Unbounded,
        ray(),
        axis(),
        wavelength,
        length,
        None,
        event::GradientEventOptions::default(),
        || true,
    )
    .unwrap();
    let lambda_m: f64 = wavelength * 1e-9;
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0_f64;
    let k = 1.380_649e-23;
    let blackbody = (2.0 * h * c * c / lambda_m.powi(5)) * 1e-9
        / (h * c / (lambda_m * k * temperature)).exp_m1();
    let expected = blackbody * -(-0.6_f64 * length).exp_m1();
    close(emission, expected, 2e-11 * expected);
    close(path.propagation.optical_depth, 0.6 * length, 2e-14);
}
