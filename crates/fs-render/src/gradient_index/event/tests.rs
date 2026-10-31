use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};

struct Planes {
    positions: Vec<f64>,
}

impl BoundaryOracle for Planes {
    fn query(&mut self, point: Point3) -> Result<BoundaryQuery, GradientIndexError> {
        let Some((index, &x)) = self
            .positions
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| (point.x - **a).abs().total_cmp(&(point.x - **b).abs()))
        else {
            return Ok(BoundaryQuery {
                clearance_m: f64::INFINITY,
                witness: None,
            });
        };
        let gap = point.x - x;
        let distance = gap.abs();
        let departure_radius = self
            .positions
            .iter()
            .enumerate()
            .filter(|(other, _)| *other != index)
            .map(|(_, other)| 0.25 * (other - x).abs())
            .fold(1.0_f64, f64::min);
        Ok(BoundaryQuery {
            clearance_m: distance.next_down().max(0.0),
            witness: Some(BoundaryWitness {
                boundary_id: index,
                feature_id: index,
                surface_point: Point3::new(x, point.y, point.z),
                normal: Vec3::new(1.0, 0.0, 0.0),
                distance_upper_m: if distance == 0.0 {
                    0.0
                } else {
                    distance.next_up()
                },
                side: if gap > 0.0 {
                    BoundarySide::Positive
                } else if gap < 0.0 {
                    BoundarySide::Negative
                } else {
                    BoundarySide::OnSurface
                },
                departure_radius_m: Some(departure_radius),
            }),
        })
    }
}

fn glass() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::try_constant(0.3).unwrap(),
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

#[test]
fn g1_uniform_boundary_arrival_retains_integrated_position_and_power() {
    let glass = glass();
    let field = HomogeneousIndexField::new(&glass);
    let options = GradientEventOptions::default();
    let result = propagate_to_boundary(
        &field,
        &mut Planes {
            positions: vec![1.0],
        },
        ray(),
        axis(),
        550.0,
        4.0,
        None,
        options,
    )
    .unwrap();
    let BoundaryTermination::Surface(witness) = result.termination else {
        panic!("expected contact");
    };
    assert_eq!(witness.boundary_id, 0);
    assert!(witness.distance_upper_m <= options.surface_tolerance_m);
    assert!(result.propagation.ray.origin.x < 1.0);
    assert!(1.0 - result.propagation.ray.origin.x <= options.surface_tolerance_m);
    assert_eq!(
        result.propagation.distance_m,
        result.propagation.ray.origin.x
    );
    assert_eq!(
        result.propagation.optical_path_length_m,
        1.5 * result.propagation.distance_m
    );
    assert_eq!(
        result.propagation.optical_depth,
        0.3 * result.propagation.distance_m
    );
    assert_eq!(result.propagation.field_evaluations, 1);
    assert!(result.contact_resolution_m > 0.0);
}

#[test]
fn g1_curved_affine_ray_hits_plane_at_analytic_height() {
    let glass = glass();
    let field = AffineIndexField::try_new(
        &glass,
        ray().origin,
        Vec3::new(0.0, 0.8, 0.0),
        Aabb::new(Point3::new(-1.0, -0.5, -1.0), Point3::new(2.0, 1.0, 1.0)),
    )
    .unwrap();
    let result = propagate_to_boundary(
        &field,
        &mut Planes {
            positions: vec![0.9],
        },
        ray(),
        axis(),
        550.0,
        2.0,
        None,
        GradientEventOptions::default(),
    )
    .unwrap();
    assert!(matches!(
        result.termination,
        BoundaryTermination::Surface(_)
    ));
    let point = result.propagation.ray.origin;
    let expected_y = 1.5 / 0.8 * ((0.8 * point.x / 1.5).cosh() - 1.0);
    assert!((point.y - expected_y).abs() < 2.0e-9);
    assert!(
        (result.propagation.final_sample.phase_index() * result.propagation.ray.dir.x - 1.5).abs()
            < 2.0e-9
    );
    assert!(result.propagation.distance_m > point.x);
    assert!(
        (result.propagation.optical_depth - 0.3 * result.propagation.distance_m).abs() < 1.0e-13
    );
}

#[test]
fn g3_departure_does_not_skip_thin_opposite_surface() {
    let glass = glass();
    let field = HomogeneousIndexField::new(&glass);
    let mut planes = Planes {
        positions: vec![1.0, 1.00001],
    };
    let options = GradientEventOptions::default();
    let first = propagate_to_boundary(
        &field,
        &mut planes,
        ray(),
        axis(),
        550.0,
        2.0,
        None,
        options,
    )
    .unwrap();
    let BoundaryTermination::Surface(witness) = first.termination else {
        panic!("expected first face");
    };
    assert_eq!(witness.boundary_id, 0);
    let second = propagate_to_boundary(
        &field,
        &mut planes,
        first.propagation.ray,
        first.propagation.transverse_axis,
        550.0,
        0.1,
        Some(GradientDeparture {
            witness,
            side: BoundarySide::Positive,
        }),
        options,
    )
    .unwrap();
    let BoundaryTermination::Surface(next) = second.termination else {
        panic!("skipped thin face");
    };
    assert_eq!(next.boundary_id, 1);
    assert!(second.propagation.distance_m > 0.0 && second.propagation.distance_m < 1.01e-5);
    let total_distance = first.propagation.distance_m + second.propagation.distance_m;
    assert!((total_distance - second.propagation.ray.origin.x).abs() < 1.0e-14);
    assert!(
        (first.propagation.optical_depth + second.propagation.optical_depth - 0.3 * total_distance)
            .abs()
            < 1.0e-14
    );
}

#[test]
fn g3_reflection_departure_integrates_the_return_path() {
    let glass = glass();
    let field = HomogeneousIndexField::new(&glass);
    let mut planes = Planes {
        positions: vec![1.0],
    };
    let options = GradientEventOptions::default();
    let first = propagate_to_boundary(
        &field,
        &mut planes,
        ray(),
        axis(),
        550.0,
        2.0,
        None,
        options,
    )
    .unwrap();
    let BoundaryTermination::Surface(witness) = first.termination else {
        panic!("expected face");
    };
    let reflected = Ray {
        origin: first.propagation.ray.origin,
        dir: first.propagation.ray.dir.scale(-1.0),
    };
    let result = propagate_to_boundary(
        &field,
        &mut planes,
        reflected,
        axis(),
        550.0,
        0.5,
        Some(GradientDeparture {
            witness,
            side: BoundarySide::Negative,
        }),
        options,
    )
    .unwrap();
    assert_eq!(result.termination, BoundaryTermination::DistanceLimit);
    assert_eq!(result.propagation.distance_m, 0.5);
    assert!((result.propagation.ray.origin.x - (reflected.origin.x - 0.5)).abs() < 1.0e-14);
    assert_eq!(result.propagation.optical_depth, 0.15);
}

#[test]
fn g0_residual_only_or_grazing_contact_is_not_a_crossing() {
    struct ResidualOnly;
    impl BoundaryOracle for ResidualOnly {
        fn query(&mut self, _: Point3) -> Result<BoundaryQuery, GradientIndexError> {
            Ok(BoundaryQuery {
                clearance_m: 0.0,
                witness: None,
            })
        }
    }
    let glass = glass();
    let field = HomogeneousIndexField::new(&glass);
    assert_eq!(
        propagate_to_boundary(
            &field,
            &mut ResidualOnly,
            ray(),
            axis(),
            550.0,
            1.0,
            None,
            GradientEventOptions::default()
        ),
        Err(GradientIndexError::BoundaryResolutionLimit)
    );
    let grazing = Ray {
        origin: Point3::new(1.0 - 1.0e-8, 0.0, 0.0),
        dir: Vec3::new(0.0, 1.0, 0.0),
    };
    assert_eq!(
        propagate_to_boundary(
            &field,
            &mut Planes {
                positions: vec![1.0]
            },
            grazing,
            axis(),
            550.0,
            1.0,
            None,
            GradientEventOptions::default()
        ),
        Err(GradientIndexError::BoundaryResolutionLimit)
    );
}

// The plane x=1, truncated at y=0. An interior witness's nearest patch edge
// is exactly -y away; positive-y rays can miss it at shallow incidence.
struct HalfPlane;
impl BoundaryOracle for HalfPlane {
    fn query(&mut self, point: Point3) -> Result<BoundaryQuery, GradientIndexError> {
        let surface_point = Point3::new(1.0, point.y.min(0.0), point.z);
        let gap = point.x - 1.0;
        Ok(BoundaryQuery {
            clearance_m: gap.abs().next_down().max(0.0),
            witness: Some(BoundaryWitness {
                boundary_id: 0,
                feature_id: 0,
                surface_point,
                normal: Vec3::new(1.0, 0.0, 0.0),
                distance_upper_m: distance_upper(point, surface_point),
                side: if gap < 0.0 {
                    BoundarySide::Negative
                } else if gap > 0.0 {
                    BoundarySide::Positive
                } else {
                    BoundarySide::OnSurface
                },
                departure_radius_m: (point.y < 0.0).then(|| (-point.y).next_down()),
            }),
        })
    }
}

#[test]
fn g0_near_edge_tangent_miss_is_not_a_contact() {
    let glass = glass();
    let field = HomogeneousIndexField::new(&glass);
    let ray = Ray {
        origin: Point3::new(1.0 - 5.0e-8, -2.0e-6, 0.0),
        dir: Vec3::new(0.01, (1.0_f64 - 0.0001).sqrt(), 0.0),
    };
    let plane_arc = (1.0 - ray.origin.x) / ray.dir.x;
    assert!(ray.origin.y + plane_arc * ray.dir.y > 0.0);
    assert_eq!(
        propagate_to_boundary(
            &field,
            &mut HalfPlane,
            ray,
            axis(),
            550.0,
            0.01,
            None,
            GradientEventOptions::default(),
        ),
        Err(GradientIndexError::BoundaryResolutionLimit)
    );
}

#[test]
fn g1_low_angle_interior_contact_fits_the_planar_patch() {
    let glass = glass();
    let field = HomogeneousIndexField::new(&glass);
    let ray = Ray {
        origin: Point3::new(1.0 - 2.0e-6, -1.0e-3, 0.0),
        dir: Vec3::new(0.01, (1.0_f64 - 0.0001).sqrt(), 0.0),
    };
    let plane_arc = (1.0 - ray.origin.x) / ray.dir.x;
    assert!(ray.origin.y + plane_arc * ray.dir.y < 0.0);
    let result = propagate_to_boundary(
        &field,
        &mut HalfPlane,
        ray,
        axis(),
        550.0,
        0.01,
        None,
        GradientEventOptions::default(),
    )
    .unwrap();
    assert!(matches!(
        result.termination,
        BoundaryTermination::Surface(_)
    ));
    assert!(result.propagation.distance_m > 0.0);
    assert!(result.propagation.distance_m < plane_arc);
    assert!(plane_arc - result.propagation.distance_m <= result.contact_resolution_m * 1.00000001);
    assert!(result.propagation.ray.origin.y < 0.0);
}

#[test]
fn g0_departure_requires_resolved_patch_and_field_collar() {
    let glass = glass();
    let field = HomogeneousIndexField::new(&glass);
    let mut planes = Planes {
        positions: vec![1.0],
    };
    let first = propagate_to_boundary(
        &field,
        &mut planes,
        ray(),
        axis(),
        550.0,
        2.0,
        None,
        GradientEventOptions::default(),
    )
    .unwrap();
    let BoundaryTermination::Surface(mut witness) = first.termination else {
        panic!("expected face");
    };
    witness.departure_radius_m = None;
    assert_eq!(
        propagate_to_boundary(
            &field,
            &mut planes,
            first.propagation.ray,
            axis(),
            550.0,
            1.0,
            Some(GradientDeparture {
                witness,
                side: BoundarySide::Positive
            }),
            GradientEventOptions::default()
        ),
        Err(GradientIndexError::UnsupportedBoundaryDeparture)
    );
    witness.departure_radius_m = Some(1.0);
    struct StrictOutside;
    impl GradientIndexField for StrictOutside {
        fn sample(&self, point: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            if point.x < 1.0 {
                return Err(GradientIndexError::OutsideFieldDomain);
            }
            GradientIndexSample::try_new(1.0, Vec3::new(0.0, 0.0, 0.0), 0.0, Some(1.0))
        }
    }
    assert_eq!(
        propagate_to_boundary(
            &StrictOutside,
            &mut planes,
            first.propagation.ray,
            axis(),
            550.0,
            1.0,
            Some(GradientDeparture {
                witness,
                side: BoundarySide::Positive
            }),
            GradientEventOptions::default()
        ),
        Err(GradientIndexError::OutsideFieldDomain)
    );
}

#[test]
fn g3_empty_geometry_preserves_adaptive_propagation_bits() {
    let glass = glass();
    let field = AffineIndexField::try_new(
        &glass,
        ray().origin,
        Vec3::new(0.0, 0.8, 0.0),
        Aabb::new(Point3::new(-1.0, -0.5, -1.0), Point3::new(2.0, 1.0, 1.0)),
    )
    .unwrap();
    let options = GradientEventOptions::default();
    let reference = propagate(&field, ray(), axis(), 550.0, 1.0, options.integration).unwrap();
    let result = propagate_to_boundary(
        &field,
        &mut Planes { positions: vec![] },
        ray(),
        axis(),
        550.0,
        1.0,
        None,
        options,
    )
    .unwrap();
    assert_eq!(result.termination, BoundaryTermination::DistanceLimit);
    assert_eq!(result.propagation, reference);
}

#[test]
fn g4_query_budget_and_cancellation_stop_the_walk() {
    let glass = glass();
    let field = HomogeneousIndexField::new(&glass);
    let mut planes = Planes {
        positions: vec![1.0],
    };
    assert_eq!(
        propagate_to_boundary(
            &field,
            &mut planes,
            ray(),
            axis(),
            550.0,
            2.0,
            None,
            GradientEventOptions {
                max_boundary_queries: 1,
                ..GradientEventOptions::default()
            }
        ),
        Err(GradientIndexError::BoundaryQueryBudgetExceeded)
    );
    let mut checks = 0;
    assert_eq!(
        propagate_to_boundary_with_control(
            &field,
            &mut planes,
            ray(),
            axis(),
            550.0,
            2.0,
            None,
            GradientEventOptions::default(),
            || {
                checks += 1;
                checks < 4
            }
        ),
        Err(GradientIndexError::Cancelled)
    );
    assert_eq!(checks, 4);
}
