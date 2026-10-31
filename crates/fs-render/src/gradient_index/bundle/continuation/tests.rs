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
fn lens(material: &DielectricGlass, curvature: f64) -> QuadraticIndexField<'_> {
    QuadraticIndexField::try_new(
        material,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 0.0),
        IndexHessian::try_new([curvature, curvature, 0.0], [0.0; 3]).unwrap(),
        Aabb::new(Point3::new(-0.2, -0.2, -0.1), Point3::new(0.2, 0.2, 4.0)),
    )
    .unwrap()
}

#[test]
fn g1_uniform_split_retains_source_solid_angle_and_all_optical_history() {
    let material = glass(1.5, 0.3);
    let field = HomogeneousIndexField::new(&material);
    let first = propagate_bundle(&field, z_ray(), x_axis(), 527.0, 0.75, tight()).unwrap();
    let seed = first.continuation();
    assert_eq!(seed.at_endpoint(), first);
    let split = propagate_seeded_bundle(&field, &seed, 1.25, tight()).unwrap();
    let direct = propagate_bundle(&field, z_ray(), x_axis(), 527.0, 2.0, tight()).unwrap();
    assert_eq!(split.initial_ray(), direct.initial_ray());
    assert_eq!(split.launch_basis(), direct.launch_basis());
    assert_eq!(split.wavelength_nm(), 527.0);
    assert_eq!(split.interface_events(), 0);
    assert_eq!(split.angular_derivatives(), direct.angular_derivatives());
    assert_eq!(split.propagation().ray, direct.propagation().ray);
    assert_eq!(split.propagation().distance_m, 2.0);
    assert_eq!(split.propagation().optical_path_length_m, 3.0);
    assert_eq!(split.propagation().group_optical_path_length_m, Some(3.0));
    close(split.propagation().optical_depth, 0.6, 2e-16);
    assert_eq!(split.propagation().field_evaluations, 2);
    assert_eq!(split.propagation().accepted_steps, 2);
    close(
        split
            .project_at_endpoint(z_ray().dir, ProjectionOptions::default())
            .unwrap()
            .area_per_solid_angle_m2(),
        4.0,
        1e-14,
    );
}

#[test]
fn g1_curved_split_obeys_the_original_analytic_angular_flow() {
    let material = glass(1.5, 0.2);
    let field = lens(&material, -1.5);
    let first = propagate_bundle(&field, z_ray(), x_axis(), 550.0, 0.4, tight()).unwrap();
    let split = propagate_seeded_bundle(&field, &first.continuation(), 0.7, tight()).unwrap();
    let direct = propagate_bundle(&field, z_ray(), x_axis(), 550.0, 1.1, tight()).unwrap();
    for ((actual, reference), axis) in split
        .angular_derivatives()
        .into_iter()
        .zip(direct.angular_derivatives())
        .zip(split.launch_basis())
    {
        close_vec(actual.position_m_per_rad, axis.scale(1.1_f64.sin()), 2e-11);
        close_vec(actual.direction_per_rad, axis.scale(1.1_f64.cos()), 2e-11);
        close_vec(
            actual.position_m_per_rad,
            reference.position_m_per_rad,
            2e-11,
        );
        close_vec(actual.direction_per_rad, reference.direction_per_rad, 2e-11);
    }
    assert_eq!(split.initial_ray(), first.initial_ray());
    assert_eq!(split.launch_basis(), first.launch_basis());
    close(split.propagation().distance_m, 1.1, 1e-15);
    close(split.propagation().optical_path_length_m, 1.65, 3e-13);
    close(split.propagation().optical_depth, 0.22, 2e-14);
    assert!(split.propagation().accepted_steps > first.propagation().accepted_steps);
    assert!(split.propagation().field_evaluations > first.propagation().field_evaluations);
}

#[test]
fn g1_mapped_uniform_columns_retain_zero_distance_endpoint_and_history() {
    let incident = glass(1.2, 0.3);
    let target = glass(1.8, 0.6);
    let field = HomogeneousIndexField::new(&target);
    let first = propagate_bundle(
        &HomogeneousIndexField::new(&incident),
        z_ray(),
        x_axis(),
        550.0,
        0.75,
        tight(),
    )
    .unwrap();
    let outgoing = Ray {
        origin: first.propagation().ray.origin,
        dir: Vec3::new(0.0, 1.0, 0.0),
    };
    let columns = [
        AngularDerivative {
            position_m_per_rad: Vec3::new(0.9, 0.2, 0.1),
            direction_per_rad: Vec3::new(0.4, 0.0, 0.3),
        },
        AngularDerivative {
            position_m_per_rad: Vec3::new(-0.2, -0.1, 0.7),
            direction_per_rad: Vec3::new(-0.6, 0.0, 0.8),
        },
    ];
    // Exercise representation admission with arbitrary tangent columns. The
    // separate physical interface tests supply actual differentiated Snell maps.
    let sample = field.sample(outgoing.origin, 550.0).unwrap();
    let seed =
        BundleContinuation::from_interface(&first, outgoing, x_axis(), sample, columns).unwrap();
    let endpoint = seed.at_endpoint();
    assert_eq!(endpoint.propagation().ray, outgoing);
    assert_eq!(endpoint.propagation().final_sample, sample);
    assert_eq!(endpoint.angular_derivatives(), columns);
    assert_eq!(
        endpoint.propagation().distance_m,
        first.propagation().distance_m
    );
    assert_eq!(
        endpoint.propagation().optical_path_length_m,
        first.propagation().optical_path_length_m
    );
    assert_eq!(
        endpoint.propagation().optical_depth,
        first.propagation().optical_depth
    );
    assert_eq!(
        endpoint.propagation().field_evaluations,
        first.propagation().field_evaluations
    );
    assert_eq!(endpoint.interface_events(), 1);
    let result = propagate_seeded_bundle(&field, &seed, 1.25, tight()).unwrap();
    assert_eq!(
        result.propagation().ray.origin,
        Point3::new(0.0, 1.25, 0.75)
    );
    assert_eq!(result.initial_ray(), first.initial_ray());
    assert_eq!(result.launch_basis(), first.launch_basis());
    assert_eq!(result.interface_events(), 1);
    assert_eq!(
        result.propagation().initial_sample,
        first.propagation().initial_sample
    );
    for (actual, original) in result.angular_derivatives().into_iter().zip(columns) {
        close_vec(
            actual.position_m_per_rad,
            sub(
                original.position_m_per_rad,
                original.direction_per_rad.scale(-1.25),
            ),
            1e-15,
        );
        assert_eq!(actual.direction_per_rad, original.direction_per_rad);
    }
    close(result.propagation().optical_path_length_m, 3.15, 5e-16);
    close(
        result.propagation().group_optical_path_length_m.unwrap(),
        3.15,
        5e-16,
    );
    close(result.propagation().optical_depth, 0.975, 2e-16);
    assert_eq!(result.propagation().field_evaluations, 2);
}

#[test]
fn g1_arbitrary_seeded_columns_obey_exact_harmonic_and_hyperbolic_flow() {
    let material = glass(1.5, 0.2);
    let first = propagate_bundle(
        &HomogeneousIndexField::new(&material),
        z_ray(),
        x_axis(),
        550.0,
        0.7,
        tight(),
    )
    .unwrap();
    let columns = [
        AngularDerivative {
            position_m_per_rad: Vec3::new(0.8, 0.2, 0.1),
            direction_per_rad: Vec3::new(0.3, -0.4, 0.0),
        },
        AngularDerivative {
            position_m_per_rad: Vec3::new(-0.1, 0.7, -0.2),
            direction_per_rad: Vec3::new(0.5, 0.2, 0.0),
        },
    ];
    let distance = 1.2_f64;
    for curvature in [-1.5, 1.5] {
        let field = lens(&material, curvature);
        let seed = BundleContinuation::from_interface(
            &first,
            first.propagation().ray,
            x_axis(),
            field.sample(first.propagation().ray.origin, 550.0).unwrap(),
            columns,
        )
        .unwrap();
        let result = propagate_seeded_bundle(&field, &seed, distance, tight()).unwrap();
        let (c, s, sign) = if curvature < 0.0 {
            (distance.cos(), distance.sin(), -1.0)
        } else {
            (distance.cosh(), distance.sinh(), 1.0)
        };
        for (actual, initial) in result.angular_derivatives().into_iter().zip(columns) {
            let x = initial.position_m_per_rad;
            let v = initial.direction_per_rad;
            close_vec(
                actual.position_m_per_rad,
                Vec3::new(c * x.x + s * v.x, c * x.y + s * v.y, x.z),
                3e-11,
            );
            close_vec(
                actual.direction_per_rad,
                Vec3::new(sign * s * x.x + c * v.x, sign * s * x.y + c * v.y, 0.0),
                3e-11,
            );
        }
        close(
            result.propagation().optical_depth,
            0.2 * (0.7 + distance),
            2e-13,
        );
        assert_eq!(result.initial_ray(), first.initial_ray());
        assert_eq!(result.interface_events(), 1);
    }
}

#[test]
fn g0_continuation_refuses_relocation_non_tangent_columns_and_wrong_field() {
    let material = glass(1.5, 0.2);
    let field = HomogeneousIndexField::new(&material);
    let first = propagate_bundle(&field, z_ray(), x_axis(), 550.0, 0.5, tight()).unwrap();
    let ray = first.propagation().ray;
    let sample = first.propagation().final_sample;
    let columns = first.angular_derivatives();
    for invalid in [
        Ray {
            origin: ray.origin.offset(x_axis().scale(1e-10)),
            ..ray
        },
        Ray {
            dir: ray.dir.scale(1.01),
            ..ray
        },
    ] {
        assert_eq!(
            BundleContinuation::from_interface(&first, invalid, x_axis(), sample, columns),
            Err(BundleError::InvalidContinuation)
        );
    }
    for invalid in [
        AngularDerivative {
            direction_per_rad: ray.dir,
            ..columns[0]
        },
        AngularDerivative {
            position_m_per_rad: Vec3::new(f64::NAN, 0.0, 0.0),
            ..columns[0]
        },
    ] {
        assert_eq!(
            BundleContinuation::from_interface(
                &first,
                ray,
                x_axis(),
                sample,
                [invalid, columns[1]]
            ),
            Err(BundleError::InvalidContinuation)
        );
    }
    let seed = first.continuation();
    let changed_absorption = glass(1.5, 0.3);
    assert_eq!(
        propagate_seeded_bundle(
            &HomogeneousIndexField::new(&changed_absorption),
            &seed,
            0.5,
            tight()
        ),
        Err(BundleError::InitialSampleMismatch)
    );
    assert_eq!(seed.at_endpoint(), first);
}

#[test]
fn g4_seeded_work_refusal_and_initial_cancellation_leave_the_seed_reusable() {
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
    let first = propagate_bundle(&field.field, z_ray(), x_axis(), 550.0, 0.3, tight()).unwrap();
    let seed = first.continuation();
    assert_eq!(
        propagate_seeded_bundle_with_control(&field, &seed, 1.0, tight(), || false),
        Err(BundleError::Propagation(GradientIndexError::Cancelled))
    );
    assert_eq!(field.queries.get(), 0);
    let limited = BundleOptions {
        integration: GradientIndexOptions {
            max_steps: 1,
            ..tight().integration
        },
        ..tight()
    };
    assert_eq!(
        propagate_seeded_bundle(&field, &seed, 1.0, limited),
        Err(BundleError::Propagation(
            GradientIndexError::StepBudgetExceeded
        ))
    );
    assert_eq!(seed.at_endpoint(), first);
    let completed = propagate_seeded_bundle(&field, &seed, 1.0, tight()).unwrap();
    close(
        completed.angular_derivatives()[0].position_m_per_rad.x,
        1.3_f64.sin(),
        3e-11,
    );
}

#[test]
fn g0_group_delay_is_retained_only_when_all_positive_length_legs_support_it() {
    struct Constant(GradientIndexSample);
    impl GradientIndexField for Constant {
        fn sample(&self, _: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            Ok(self.0)
        }
        fn is_spatially_uniform(&self) -> bool {
            true
        }
    }
    impl SecondOrderIndexField for Constant {
        fn sample_second_order(&self, _: Point3, _: f64) -> Result<IndexJet, GradientIndexError> {
            Ok(IndexJet::new(self.0, IndexHessian::ZERO))
        }
    }
    let known = Constant(
        GradientIndexSample::try_new(1.5, Vec3::new(0.0, 0.0, 0.0), 0.1, Some(1.6)).unwrap(),
    );
    let unknown =
        Constant(GradientIndexSample::try_new(1.5, Vec3::new(0.0, 0.0, 0.0), 0.1, None).unwrap());
    let first = propagate_bundle(&known, z_ray(), x_axis(), 550.0, 0.5, tight()).unwrap();
    let seed = BundleContinuation::from_interface(
        &first,
        first.propagation().ray,
        x_axis(),
        unknown.0,
        first.angular_derivatives(),
    )
    .unwrap();
    let no_travel = propagate_seeded_bundle(&unknown, &seed, 0.0, tight()).unwrap();
    assert_eq!(
        no_travel.propagation().group_optical_path_length_m,
        Some(0.8)
    );
    let next = propagate_seeded_bundle(&unknown, &seed, 0.25, tight()).unwrap();
    assert_eq!(next.propagation().group_optical_path_length_m, None);
    let seed = BundleContinuation::from_interface(
        &next,
        next.propagation().ray,
        x_axis(),
        known.0,
        next.angular_derivatives(),
    )
    .unwrap();
    let last = propagate_seeded_bundle(&known, &seed, 0.25, tight()).unwrap();
    assert_eq!(last.propagation().group_optical_path_length_m, None);
    close(last.propagation().optical_path_length_m, 1.5, 1e-15);
    close(last.propagation().optical_depth, 0.1, 2e-17);
    assert_eq!(last.interface_events(), 2);
}
