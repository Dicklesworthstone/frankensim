use super::*;
use crate::gradient_index::event::{
    BoundaryQuery, BoundarySide, BoundaryWitness, DepartureNeighborhood, DepartureRequest,
    SignedBoundaryDistance, propagate_to_boundary,
};
use std::cell::Cell;

struct Lens {
    calls: Cell<u32>,
    scalar_calls: Cell<u32>,
    uniform: bool,
}
impl Lens {
    fn new(uniform: bool) -> Self {
        Self {
            calls: Cell::new(0),
            scalar_calls: Cell::new(0),
            uniform,
        }
    }
    fn optical(&self, p: Point3) -> Result<GradientIndexSample, GradientIndexError> {
        let k = if self.uniform { 0.0 } else { -1.5 };
        GradientIndexSample::try_new(
            1.5 + 0.5 * k * (p.x * p.x + p.y * p.y),
            Vec3::new(k * p.x, k * p.y, 0.0),
            0.2,
            Some(1.6),
        )
    }
}
impl GradientIndexField for Lens {
    fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        self.scalar_calls.set(self.scalar_calls.get() + 1);
        self.optical(p)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.uniform
    }
}
impl SecondOrderIndexField for Lens {
    fn sample_second_order(&self, p: Point3, _: f64) -> Result<IndexJet, GradientIndexError> {
        self.calls.set(self.calls.get() + 1);
        let k = if self.uniform { 0.0 } else { -1.5 };
        Ok(IndexJet::new(
            self.optical(p)?,
            IndexHessian::try_new([k, k, 0.0], [0.0; 3]).unwrap(),
        ))
    }
}

struct Plane(f64);
impl BoundaryOracle for Plane {
    fn query(&mut self, p: Point3) -> Result<BoundaryQuery, GradientIndexError> {
        let gap = p.z - self.0;
        Ok(BoundaryQuery {
            clearance_m: gap.abs().next_down().max(0.0),
            witness: Some(BoundaryWitness {
                boundary_id: 0,
                feature_id: 0,
                surface_point: Point3::new(p.x, p.y, self.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                distance_upper_m: gap.abs().next_up(),
                side: if gap < 0.0 {
                    BoundarySide::Negative
                } else {
                    BoundarySide::Positive
                },
                departure_radius_m: Some(1.0),
            }),
        })
    }
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn axis() -> Vec3 {
    Vec3::new(1.0, 0.0, 0.0)
}
fn launch(field: &impl SecondOrderIndexField, distance: f64) -> RayBundlePropagation {
    propagate_bundle(field, ray(), axis(), 550.0, distance, options().bundle()).unwrap()
}
fn options() -> BundleEventOptions {
    BundleEventOptions {
        boundary: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-11,
                optical_path_tolerance_m: 1e-11,
                optical_depth_tolerance: 1e-11,
                relative_tolerance: 1e-12,
                initial_step_m: 0.2,
                min_step_m: 1e-12,
                max_step_m: 0.4,
                max_steps: 20_000,
            },
            ..GradientEventOptions::default()
        },
        position_derivative_tolerance_m: 1e-11,
        direction_derivative_tolerance: 1e-11,
    }
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    assert!(vec_norm(sub(a, b)) <= tolerance, "{a:?} != {b:?}");
}

#[test]
fn g1_boundary_bundle_retains_original_source_and_analytic_grin_focusing() {
    let field = Lens::new(false);
    let first = launch(&field, 0.3);
    let before = field.calls.get();
    let result = propagate_seeded_bundle_to_boundary(
        &field,
        &mut Plane(1.2),
        &first.continuation(),
        2.0,
        None,
        options(),
    )
    .unwrap();
    let BoundaryTermination::Surface(witness) = result.termination else {
        panic!("expected contact")
    };
    assert!(witness.distance_upper_m <= options().boundary.surface_tolerance_m);
    assert!(result.bundle.propagation().ray.origin.z < 1.2);
    assert_eq!(result.bundle.initial_ray(), first.initial_ray());
    assert_eq!(result.bundle.launch_basis(), first.launch_basis());
    assert_eq!(
        result.bundle.propagation().initial_sample,
        first.propagation().initial_sample
    );
    let distance = result.bundle.propagation().distance_m;
    for (column, axis) in result
        .bundle
        .angular_derivatives()
        .into_iter()
        .zip(first.launch_basis())
    {
        close_vec(column.position_m_per_rad, axis.scale(distance.sin()), 1e-9);
        close_vec(column.direction_per_rad, axis.scale(distance.cos()), 1e-9);
    }
    assert!((result.bundle.propagation().optical_path_length_m - 1.5 * distance).abs() < 1e-12);
    assert!(
        (result
            .bundle
            .propagation()
            .group_optical_path_length_m
            .unwrap()
            - 1.6 * distance)
            .abs()
            < 1e-12
    );
    assert!((result.bundle.propagation().optical_depth - 0.2 * distance).abs() < 1e-12);
    assert_eq!(field.calls.get() - before, result.segment.field_evaluations);
    assert_eq!(field.scalar_calls.get(), 0, "no separate central-ray trace");
    assert_eq!(
        result.segment.field_evaluations,
        1 + 6 * (result.segment.accepted_steps + result.segment.rejected_steps)
    );
    assert_eq!(
        result.bundle.propagation().field_evaluations,
        first.propagation().field_evaluations + result.segment.field_evaluations
    );
}

#[test]
fn g1_uniform_boundary_bundle_uses_one_jet_and_the_same_optical_event_path() {
    let field = Lens::new(true);
    let first = launch(&field, 0.25);
    let before = field.calls.get();
    let result = propagate_seeded_bundle_to_boundary(
        &field,
        &mut Plane(1.0),
        &first.continuation(),
        2.0,
        None,
        options(),
    )
    .unwrap();
    let optical = propagate_to_boundary(
        &field,
        &mut Plane(1.0),
        first.propagation().ray,
        first.propagation().transverse_axis,
        550.0,
        2.0,
        None,
        options().boundary,
    )
    .unwrap();
    assert_eq!(result.segment, optical.propagation);
    assert_eq!(result.termination, optical.termination);
    assert_eq!(result.boundary_evaluations, optical.boundary_evaluations);
    assert_eq!(result.contact_resolution_m, optical.contact_resolution_m);
    assert_eq!(field.calls.get() - before, 1);
    assert_eq!(result.segment.field_evaluations, 1);
    for (column, axis) in result
        .bundle
        .angular_derivatives()
        .into_iter()
        .zip(first.launch_basis())
    {
        close_vec(
            column.position_m_per_rad,
            axis.scale(result.bundle.propagation().distance_m),
            1e-15,
        );
        assert_eq!(column.direction_per_rad, axis);
    }
}

#[test]
fn g1_boundary_steps_are_controlled_by_hessian_columns_even_on_a_straight_axis() {
    let field = Lens::new(false);
    let first = launch(&field, 0.0);
    let mut loose = options();
    loose.boundary.integration.relative_tolerance = 0.0;
    loose.boundary.integration.initial_step_m = 0.8;
    loose.boundary.integration.max_step_m = 0.8;
    loose.position_derivative_tolerance_m = 1e-3;
    loose.direction_derivative_tolerance = 1e-3;
    let mut tight = loose;
    tight.position_derivative_tolerance_m = 1e-12;
    tight.direction_derivative_tolerance = 1e-12;
    let run = |options| {
        propagate_seeded_bundle_to_boundary(
            &field,
            &mut Plane(1.2),
            &first.continuation(),
            2.0,
            None,
            options,
        )
        .unwrap()
    };
    let coarse = run(loose);
    let precise = run(tight);
    assert!(precise.segment.accepted_steps > coarse.segment.accepted_steps);
    assert!(precise.segment.rejected_steps > 0);
    let distance = precise.bundle.propagation().distance_m;
    close_vec(
        precise.bundle.angular_derivatives()[0].position_m_per_rad,
        axis().scale(distance.sin()),
        1e-10,
    );
}

struct SmoothDeparture {
    prepare: u32,
    signed: u32,
    failure: Option<GradientIndexError>,
}
impl BoundaryOracle for SmoothDeparture {
    fn query(&mut self, _: Point3) -> Result<BoundaryQuery, GradientIndexError> {
        panic!("must refuse inside the departure collar");
    }
    fn prepare_departure(
        &mut self,
        _: DepartureRequest,
    ) -> Result<DepartureNeighborhood, GradientIndexError> {
        self.prepare += 1;
        Ok(DepartureNeighborhood::Smooth {
            radius_m: 1.0,
            normal_axis: Vec3::new(0.0, 0.0, 1.0),
            max_normal_deviation: 0.0,
        })
    }
    fn signed_departure(
        &mut self,
        _: BoundaryWitness,
        p: Point3,
    ) -> Result<SignedBoundaryDistance, GradientIndexError> {
        self.signed += 1;
        if self.signed == 2
            && let Some(error) = self.failure
        {
            return Err(error);
        }
        Ok(SignedBoundaryDistance {
            lower_m: p.z.next_down(),
            upper_m: p.z.next_up(),
        })
    }
}

#[test]
fn g4_hessian_stages_preserve_signed_departure_errors_budgets_and_cancellation() {
    let field = Lens::new(false);
    let ray = Ray {
        origin: Point3::new(0.0, 0.0, -0.25e-7),
        ..ray()
    };
    let first = propagate_bundle(&field, ray, axis(), 550.0, 0.0, options().bundle()).unwrap();
    let departure = GradientDeparture {
        witness: BoundaryWitness {
            boundary_id: 0,
            feature_id: 0,
            surface_point: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            distance_upper_m: 0.25e-7,
            side: BoundarySide::Negative,
            departure_radius_m: None,
        },
        side: BoundarySide::Positive,
    };
    for failure in [
        Some(GradientIndexError::Cancelled),
        Some(GradientIndexError::BoundaryQueryBudgetExceeded),
        None,
    ] {
        let mut geometry = SmoothDeparture {
            prepare: 0,
            signed: 0,
            failure,
        };
        let mut controls = options();
        if failure.is_none() {
            controls.boundary.max_boundary_queries = 2;
        }
        let before = field.calls.get();
        let result = propagate_seeded_bundle_to_boundary(
            &field,
            &mut geometry,
            &first.continuation(),
            1.0,
            Some(departure),
            controls,
        );
        assert_eq!(
            result,
            Err(BundleError::Propagation(
                failure.unwrap_or(GradientIndexError::BoundaryQueryBudgetExceeded)
            ))
        );
        assert_eq!(
            field.calls.get() - before,
            1,
            "signed guard precedes the Hessian query"
        );
        assert_eq!(geometry.prepare, 1);
        assert_eq!(geometry.signed, if failure.is_some() { 2 } else { 1 });
    }
    let mut geometry = SmoothDeparture {
        prepare: 0,
        signed: 0,
        failure: None,
    };
    let before = field.calls.get();
    assert_eq!(
        propagate_seeded_bundle_to_boundary_with_control(
            &field,
            &mut geometry,
            &first.continuation(),
            1.0,
            Some(departure),
            options(),
            || false,
        ),
        Err(BundleError::Propagation(GradientIndexError::Cancelled))
    );
    assert_eq!(field.calls.get(), before);
    assert_eq!(geometry.prepare, 0);
    assert_eq!(geometry.signed, 0);
}

#[test]
fn g4_boundary_bundle_refuses_an_unmatched_next_material_without_retracing() {
    let field = Lens::new(false);
    let first = launch(&field, 0.2);
    struct Changed;
    impl GradientIndexField for Changed {
        fn sample(&self, _: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            GradientIndexSample::try_new(1.6, Vec3::new(0.0, 0.0, 0.0), 0.2, Some(1.6))
        }
    }
    impl SecondOrderIndexField for Changed {
        fn sample_second_order(
            &self,
            p: Point3,
            wavelength: f64,
        ) -> Result<IndexJet, GradientIndexError> {
            Ok(IndexJet::new(
                self.sample(p, wavelength)?,
                IndexHessian::ZERO,
            ))
        }
    }
    assert_eq!(
        propagate_seeded_bundle_to_boundary(
            &Changed,
            &mut Plane(1.0),
            &first.continuation(),
            1.0,
            None,
            options(),
        ),
        Err(BundleError::InitialSampleMismatch)
    );
}

#[test]
fn g4_rejected_internal_stages_never_query_a_hessian_outside_the_free_ball() {
    struct Steep<'a> {
        center: &'a Cell<Point3>,
        calls: Cell<u32>,
    }
    impl GradientIndexField for Steep<'_> {
        fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            let q = (100.0 * p.x).tanh();
            GradientIndexSample::try_new(
                1.5 + 0.4 * q,
                Vec3::new(40.0 * (1.0 - q * q), 0.0, 0.0),
                0.0,
                None,
            )
        }
    }
    impl SecondOrderIndexField for Steep<'_> {
        fn sample_second_order(
            &self,
            p: Point3,
            wavelength: f64,
        ) -> Result<IndexJet, GradientIndexError> {
            assert!(
                vec_norm(p.delta_from(self.center.get())) <= 0.1,
                "an inadmissible stage reached the authored Hessian field"
            );
            self.calls.set(self.calls.get() + 1);
            let q = (100.0 * p.x).tanh();
            Ok(IndexJet::new(
                self.sample(p, wavelength)?,
                IndexHessian::try_new([-8000.0 * q * (1.0 - q * q), 0.0, 0.0], [0.0; 3]).unwrap(),
            ))
        }
    }
    struct FreeBall<'a>(&'a Cell<Point3>);
    impl BoundaryOracle for FreeBall<'_> {
        fn query(&mut self, point: Point3) -> Result<BoundaryQuery, GradientIndexError> {
            self.0.set(point);
            Ok(BoundaryQuery {
                clearance_m: 0.1,
                witness: None,
            })
        }
    }
    let center = Cell::new(ray().origin);
    let field = Steep {
        center: &center,
        calls: Cell::new(0),
    };
    let first = launch(&field, 0.0);
    let mut controls = options();
    controls.boundary.integration.initial_step_m = 0.1;
    controls.boundary.integration.max_step_m = 0.1;
    let before = field.calls.get();
    let result = propagate_seeded_bundle_to_boundary(
        &field,
        &mut FreeBall(&center),
        &first.continuation(),
        0.2,
        None,
        controls,
    )
    .unwrap();
    assert_eq!(result.termination, BoundaryTermination::DistanceLimit);
    assert_eq!(field.calls.get() - before, result.segment.field_evaluations);
    assert!(result.segment.rejected_steps > 0);
    assert!(
        result.segment.field_evaluations
            < 1 + 6 * (result.segment.accepted_steps + result.segment.rejected_steps),
        "an out-of-ball trial must stop before spending all six Hessian queries"
    );
}
