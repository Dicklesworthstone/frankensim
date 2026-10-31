use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};

fn glass(n: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}
fn z_ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn x_axis() -> Vec3 {
    Vec3::new(1.0, 0.0, 0.0)
}
fn tight() -> BundleOptions {
    BundleOptions {
        integration: GradientIndexOptions {
            position_tolerance_m: 1e-12,
            direction_tolerance: 1e-12,
            optical_path_tolerance_m: 1e-12,
            optical_depth_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            initial_step_m: 0.05,
            min_step_m: 1e-12,
            max_step_m: 0.3,
            max_steps: 20_000,
        },
        position_derivative_tolerance_m: 1e-12,
        direction_derivative_tolerance: 1e-12,
    }
}
fn lens(glass: &DielectricGlass, curvature: f64) -> QuadraticIndexField<'_> {
    QuadraticIndexField::try_new(
        glass,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 0.0),
        IndexHessian::try_new([curvature, curvature, 0.0], [0.0; 3]).unwrap(),
        Aabb::new(Point3::new(-0.2, -0.2, -0.1), Point3::new(0.2, 0.2, 4.0)),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn close_vec(actual: Vec3, expected: Vec3, tolerance: f64) {
    assert!(
        vec_norm(sub(actual, expected)) <= tolerance,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn g1_homogeneous_bundle_recovers_inverse_square_oblique_irradiance_and_metadata() {
    let material = glass(1.5, 0.3);
    let field = HomogeneousIndexField::new(&material);
    let result = propagate_bundle(
        &field,
        Ray {
            origin: Point3::new(3.0, -2.0, 1.0),
            dir: Vec3::new(0.0, 0.0, 5.0),
        },
        Vec3::new(4.0, 0.0, 0.0),
        550.0,
        2.0,
        tight(),
    )
    .unwrap();
    assert_eq!(result.initial_ray().origin, Point3::new(3.0, -2.0, 1.0));
    assert_eq!(result.initial_ray().dir, z_ray().dir);
    assert_eq!(result.wavelength_nm(), 550.0);
    assert_eq!(result.launch_basis(), [x_axis(), Vec3::new(0.0, 1.0, 0.0)]);
    assert_eq!(result.propagation().field_evaluations, 1);
    assert_eq!(result.propagation().optical_path_length_m, 3.0);
    assert_eq!(result.propagation().optical_depth, 0.6);
    for (column, axis) in result
        .angular_derivatives()
        .into_iter()
        .zip(result.launch_basis())
    {
        assert_eq!(column.position_m_per_rad, axis.scale(2.0));
        assert_eq!(column.direction_per_rad, axis);
    }
    let projection = result
        .project_at_endpoint(Vec3::new(0.0, 0.6, 0.8), ProjectionOptions::default())
        .unwrap();
    close(projection.area_per_solid_angle_m2(), 5.0, 2e-14);
    close(projection.condition_number(), 1.25, 2e-14);
    close(
        projection.point_source_irradiance(7.0).unwrap(),
        7.0 * (-0.6_f64).exp() * 0.8 / 4.0,
        2e-14,
    );
    close(
        projection.arc_length_derivatives_m_per_rad()[1],
        -1.5,
        2e-14,
    );
    close_vec(
        projection.angular_derivatives()[1].position_m_per_rad,
        Vec3::new(0.0, 2.0, -1.5),
        2e-14,
    );
}

#[test]
fn g1_quadratic_on_axis_bundle_matches_exact_focusing_and_defocusing() {
    // At the axial central ray n=1.5 and grad(n)=0. Exact first variations
    // obey X''=-X for H_xx=H_yy=-1.5, and X''=+X for the opposite sign.
    // They are not inferred from tracing finitely separated neighboring rays.
    let material = glass(1.5, 0.2);
    for curvature in [-1.5, 1.5] {
        let field = lens(&material, curvature);
        for distance in [0.4_f64, 1.2, 2.3] {
            let result =
                propagate_bundle(&field, z_ray(), x_axis(), 550.0, distance, tight()).unwrap();
            let (position, direction) = if curvature < 0.0 {
                (distance.sin(), distance.cos())
            } else {
                (distance.sinh(), distance.cosh())
            };
            for (column, axis) in result
                .angular_derivatives()
                .into_iter()
                .zip(result.launch_basis())
            {
                close_vec(column.position_m_per_rad, axis.scale(position), 2e-10);
                close_vec(column.direction_per_rad, axis.scale(direction), 2e-10);
                close(
                    column.direction_per_rad.dot(result.propagation().ray.dir),
                    0.0,
                    2e-13,
                );
            }
            let projection = result
                .project_at_endpoint(z_ray().dir, ProjectionOptions::default())
                .unwrap();
            close(
                projection.area_per_solid_angle_m2(),
                position * position,
                3e-9,
            );
            close(
                projection.point_source_irradiance(2.0).unwrap(),
                2.0 * (-0.2 * distance).exp() / (position * position),
                2e-9,
            );
            close(result.propagation().distance_m, distance, 0.0);
            close(result.propagation().optical_depth, 0.2 * distance, 2e-13);
        }
    }
}

#[test]
fn g1_affine_bundle_matches_independent_closed_form_derivatives() {
    let n = 1.5;
    let g = 0.8;
    let distance = 0.9;
    let material = glass(n, 0.0);
    let field = AffineIndexField::try_new(
        &material,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, g, 0.0),
        Aabb::new(Point3::new(-0.1, -0.5, -0.5), Point3::new(2.0, 1.0, 0.5)),
    )
    .unwrap();
    let result = propagate_bundle(
        &field,
        Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: x_axis(),
        },
        Vec3::new(0.0, 1.0, 0.0),
        550.0,
        distance,
        tight(),
    )
    .unwrap();
    // n*t_x is constant; d(n*t_y)/ds=g. Differentiating the resulting
    // analytic trajectory at launch alpha=beta=0 supplies both columns.
    let end_n = (n * n + g * g * distance * distance).sqrt();
    let first = result.angular_derivatives()[0];
    close_vec(
        first.position_m_per_rad,
        Vec3::new(n / g * (n / end_n - 1.0), n * distance / end_n, 0.0),
        2e-10,
    );
    close_vec(
        first.direction_per_rad,
        Vec3::new(
            -n * n * g * distance / end_n.powi(3),
            n.powi(3) / end_n.powi(3),
            0.0,
        ),
        2e-10,
    );
    let second = result.angular_derivatives()[1];
    close_vec(
        second.position_m_per_rad,
        Vec3::new(0.0, 0.0, n / g * (g * distance / n).asinh()),
        2e-10,
    );
    close_vec(
        second.direction_per_rad,
        Vec3::new(0.0, 0.0, n / end_n),
        2e-10,
    );
}

#[test]
fn g2_mixed_hessian_and_normalized_launch_match_neighboring_ray_falsifier() {
    let material = glass(1.6, 0.1);
    let field = QuadraticIndexField::try_new(
        &material,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.1, -0.04, 0.03),
        IndexHessian::try_new([-0.08, 0.05, 0.02], [0.03, -0.02, 0.01]).unwrap(),
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0)),
    )
    .unwrap();
    let ray = Ray {
        origin: Point3::new(0.01, 0.02, -0.03),
        dir: Vec3::new(2.0, 0.4, 0.2),
    };
    let result =
        propagate_bundle(&field, ray, Vec3::new(0.0, 0.0, 1.0), 550.0, 0.7, tight()).unwrap();
    let epsilon = 2e-4;
    for (column, axis) in result
        .angular_derivatives()
        .into_iter()
        .zip(result.launch_basis())
    {
        let trace = |sign: f64| {
            super::super::propagate(
                &field,
                Ray {
                    origin: ray.origin,
                    dir: sub(result.initial_ray().dir, axis.scale(-sign * epsilon)),
                },
                Vec3::new(0.0, 0.0, 1.0),
                550.0,
                0.7,
                tight().integration,
            )
            .unwrap()
        };
        let plus = trace(1.0);
        let minus = trace(-1.0);
        close_vec(
            column.position_m_per_rad,
            plus.ray
                .origin
                .delta_from(minus.ray.origin)
                .scale(0.5 / epsilon),
            1e-7,
        );
        close_vec(
            column.direction_per_rad,
            sub(plus.ray.dir, minus.ray.dir).scale(0.5 / epsilon),
            1e-7,
        );
    }
}

#[test]
fn g0_caustics_grazing_and_unresolved_footprints_refuse_irradiance() {
    let material = glass(1.5, 0.0);
    let focused = propagate_bundle(
        &lens(&material, -1.5),
        z_ray(),
        x_axis(),
        550.0,
        core::f64::consts::PI,
        tight(),
    )
    .unwrap();
    assert_eq!(
        focused.project_at_endpoint(z_ray().dir, ProjectionOptions::default()),
        Err(BundleError::SingularProjection)
    );
    let regular = propagate_bundle(
        &lens(&material, -1.5),
        z_ray(),
        x_axis(),
        550.0,
        1.0,
        tight(),
    )
    .unwrap();
    assert_eq!(
        regular.project_at_endpoint(x_axis(), ProjectionOptions::default()),
        Err(BundleError::GrazingProjection)
    );
    assert_eq!(
        regular.project_at_endpoint(
            z_ray().dir,
            ProjectionOptions {
                min_area_m2_per_sr: 1.0,
                ..ProjectionOptions::default()
            }
        ),
        Err(BundleError::SingularProjection)
    );
    let projected = regular
        .project_at_endpoint(z_ray().dir, ProjectionOptions::default())
        .unwrap();
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            projected.point_source_irradiance(bad),
            Err(BundleError::InvalidIntensity)
        );
    }
    let zero = propagate_bundle(
        &HomogeneousIndexField::new(&material),
        z_ray(),
        x_axis(),
        550.0,
        0.0,
        tight(),
    )
    .unwrap();
    assert_eq!(
        zero.project_at_endpoint(z_ray().dir, ProjectionOptions::default()),
        Err(BundleError::SingularProjection)
    );
}

#[test]
fn g0_footprint_and_intensity_scaling_do_not_lose_representable_results() {
    let material = glass(1.5, 800.0);
    let field = HomogeneousIndexField::new(&material);
    let result = propagate_bundle(&field, z_ray(), x_axis(), 550.0, 1.0, tight()).unwrap();
    assert_eq!(result.propagation().power_transmittance(), 0.0);
    let projection = result
        .project_at_endpoint(z_ray().dir, ProjectionOptions::default())
        .unwrap();
    let expected_half = 1e150 * (-400.0_f64).exp();
    let expected = expected_half * expected_half;
    let actual = projection.point_source_irradiance(1e300).unwrap();
    assert!(actual > 0.0);
    close(actual / expected, 1.0, 5e-13);
    let tiny = propagate_bundle(&field, z_ray(), x_axis(), 550.0, 1e-310, tight()).unwrap();
    assert_eq!(
        tiny.project_at_endpoint(z_ray().dir, ProjectionOptions::default()),
        Err(BundleError::SingularProjection)
    );
}

#[test]
fn g0_quadratic_admission_covers_interior_extrema_and_field_domain() {
    let material = glass(1.1, 0.0);
    let result = QuadraticIndexField::try_new(
        &material,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        IndexHessian::try_new([2.0, 0.0, 0.0], [0.0; 3]).unwrap(),
        Aabb::new(Point3::new(-1.0, -0.1, -0.1), Point3::new(0.0, 0.1, 0.1)),
    );
    // Every box corner has n=1.1, but the interior minimum is n=0.85.
    assert!(matches!(result, Err(BundleError::InvalidDifferential)));
    assert!(IndexHessian::try_new([f64::NAN, 0.0, 0.0], [0.0; 3]).is_err());
    let material = glass(1.5, 0.0);
    let field = lens(&material, -1.5);
    assert!(matches!(
        propagate_bundle(&field, z_ray(), x_axis(), 550.0, 5.0, tight()),
        Err(BundleError::Propagation(
            GradientIndexError::OutsideFieldDomain
        ))
    ));
}

#[test]
fn g3_derivative_error_controls_refine_an_on_axis_ray_with_zero_local_gradient() {
    let material = glass(1.5, 0.0);
    let field = lens(&material, -1.5);
    let loose = BundleOptions {
        integration: GradientIndexOptions {
            relative_tolerance: 1e-3,
            initial_step_m: 1.0,
            max_step_m: 1.0,
            ..GradientIndexOptions::default()
        },
        position_derivative_tolerance_m: 1e-3,
        direction_derivative_tolerance: 1e-3,
    };
    let coarse = propagate_bundle(&field, z_ray(), x_axis(), 550.0, 2.0, loose).unwrap();
    let fine = propagate_bundle(&field, z_ray(), x_axis(), 550.0, 2.0, tight()).unwrap();
    let expected = 2.0_f64.sin();
    assert!(
        (fine.angular_derivatives()[0].position_m_per_rad.x - expected).abs()
            < (coarse.angular_derivatives()[0].position_m_per_rad.x - expected).abs() * 0.001
    );
    assert!(fine.propagation().accepted_steps > coarse.propagation().accepted_steps);
    assert!(fine.propagation().max_normalized_error_estimate <= 1.0);
}

#[test]
fn g4_bundle_work_limit_and_rejected_step_cancellation_are_transactional() {
    use core::cell::Cell;
    struct Count<'a> {
        field: QuadraticIndexField<'a>,
        queries: Cell<u32>,
    }
    impl GradientIndexField for Count<'_> {
        fn sample(&self, p: Point3, w: f64) -> Result<GradientIndexSample, GradientIndexError> {
            self.field.sample(p, w)
        }
    }
    impl SecondOrderIndexField for Count<'_> {
        fn sample_second_order(&self, p: Point3, w: f64) -> Result<IndexJet, GradientIndexError> {
            self.queries.set(self.queries.get() + 1);
            self.field.sample_second_order(p, w)
        }
    }
    let material = glass(1.5, 0.0);
    let field = Count {
        field: lens(&material, -1.5),
        queries: Cell::new(0),
    };
    let options = BundleOptions {
        integration: GradientIndexOptions {
            initial_step_m: 1.0,
            max_step_m: 1.0,
            ..tight().integration
        },
        ..tight()
    };
    let mut polls = 0;
    assert_eq!(
        propagate_bundle_with_control(&field, z_ray(), x_axis(), 550.0, 2.0, options, || {
            polls += 1;
            polls <= 2
        }),
        Err(BundleError::Propagation(GradientIndexError::Cancelled))
    );
    assert_eq!(field.queries.get(), 7);
    assert_eq!(
        propagate_bundle(
            &field,
            z_ray(),
            x_axis(),
            550.0,
            2.0,
            BundleOptions {
                integration: GradientIndexOptions {
                    max_steps: 1,
                    ..options.integration
                },
                ..options
            }
        ),
        Err(BundleError::Propagation(
            GradientIndexError::StepBudgetExceeded
        ))
    );
}
