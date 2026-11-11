use super::super::{CurvedPolarizedScene, CurvedSceneError, CurvedTraceSettings, SceneMedium};
use super::*;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricEvent, DielectricGlass, DielectricSurface,
    GlassProvenance,
};
use crate::gradient_index::bundle::ProjectionOptions;
use crate::gradient_index::bundle::interface::{BundleInterfaceError, InterfaceOptions};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{AffineIndexField, GradientIndexOptions};
use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{
    Interaction, PathTermination, SceneError, SurfaceBinding,
};
use crate::polarization::dielectric::transport::bundle::PolarizedBundleError;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Point3, Vec3};

pub(super) fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

pub(super) fn glass(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

pub(super) fn sphere() -> GeometryInstance {
    GeometryInstance::try_new(
        1,
        hash_domain("automatic-bundle-test", b"native-unit-sphere"),
        SharedGeometry::chart(SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        }),
        RigidTransform::identity(),
    )
    .unwrap()
}

pub(super) fn plane(id: u64, center: Point3, radius: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("automatic-bundle-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [center.x - radius, center.y - radius, center.z],
                [center.x + radius, center.y - radius, center.z],
                [center.x, center.y + radius, center.z],
            ],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

pub(super) fn bindings(detector_axis: Vec3) -> Vec<SurfaceBinding> {
    vec![
        SurfaceBinding {
            object_id: 1,
            interaction: Interaction::Dielectric {
                outside: 1,
                inside: 0,
                surface: DielectricSurface::SMOOTH,
            },
        },
        SurfaceBinding {
            object_id: 9,
            interaction: Interaction::Detector {
                reference_axis: detector_axis,
            },
        },
    ]
}

pub(super) fn settings() -> CurvedBundleSettings {
    CurvedBundleSettings {
        trace: CurvedTraceSettings {
            max_surface_events: 8,
            max_distance_m: 8.0,
            max_total_steps: 16_384,
            max_query_work: 1_000_000,
            event: GradientEventOptions {
                integration: GradientIndexOptions {
                    position_tolerance_m: 1e-11,
                    direction_tolerance: 1e-12,
                    optical_path_tolerance_m: 1e-11,
                    optical_depth_tolerance: 1e-12,
                    relative_tolerance: 1e-11,
                    ..GradientIndexOptions::default()
                },
                surface_tolerance_m: 1e-8,
                ..GradientEventOptions::default()
            },
        },
        position_derivative_tolerance_m: 1e-10,
        direction_derivative_tolerance: 1e-11,
        interface: InterfaceOptions::default(),
    }
}

fn add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

pub(super) fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.dot(v).sqrt())
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

pub(super) fn launch_basis(direction: Vec3) -> [Vec3; 2] {
    let first = unit(add(Vec3::new(1.0, 0.0, 0.0), direction.scale(-direction.x)));
    [first, cross(direction, first)]
}

pub(super) fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e} (tolerance {tolerance:.3e})"
    );
}

pub(super) fn close_vec(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.x, expected.x, tolerance);
    close(actual.y, expected.y, tolerance);
    close(actual.z, expected.z, tolerance);
}

fn snell(direction: Vec3, outward: Vec3, ratio: f64) -> Vec3 {
    let normal = if direction.dot(outward) < 0.0 {
        outward
    } else {
        outward.scale(-1.0)
    };
    let ci = -direction.dot(normal);
    let ct = (1.0 - ratio * ratio * (1.0 - ci * ci)).sqrt();
    unit(add(direction.scale(ratio), normal.scale(ratio * ci - ct)))
}

/// Independent quadratic intersections and vector Snell, with no production
/// event walker, normal differential, ODE, or prescribed production leg lengths.
pub(super) fn ball_reference(ray: Ray, index: f64) -> (Ray, f64, f64, f64) {
    let origin = Vec3::new(ray.origin.x, ray.origin.y, ray.origin.z);
    let b = origin.dot(ray.dir);
    let before = -b - (b * b - origin.dot(origin) + 1.0).sqrt();
    let entry = ray.at(before);
    let entry_normal = unit(Vec3::new(entry.x, entry.y, entry.z));
    let inside = snell(ray.dir, entry_normal, 1.0 / index);
    let chord = -2.0 * inside.dot(Vec3::new(entry.x, entry.y, entry.z));
    let exit = entry.offset(inside.scale(chord));
    let outgoing = snell(inside, unit(Vec3::new(exit.x, exit.y, exit.z)), index);
    let after = (2.5 - exit.z) / outgoing.z;
    let ci = -ray.dir.dot(entry_normal);
    let ct = -inside.dot(entry_normal);
    let ts = 1.0 - ((ci - index * ct) / (ci + index * ct)).powi(2);
    let tp = 1.0 - ((index * ci - ct) / (index * ci + ct)).powi(2);
    (
        Ray {
            origin: exit.offset(outgoing.scale(after)),
            dir: outgoing,
        },
        before + after,
        chord,
        0.5 * (ts * ts + tp * tp),
    )
}

pub(super) fn reference_columns(ray: Ray, endpoint: impl Fn(Ray) -> Point3) -> [Vec3; 2] {
    let h = 1e-6;
    launch_basis(ray.dir).map(|axis| {
        let displaced = |angle| {
            endpoint(Ray {
                origin: ray.origin,
                dir: unit(add(ray.dir, axis.scale(angle))),
            })
        };
        displaced(h).delta_from(displaced(-h)).scale(0.5 / h)
    })
}

#[test]
fn g1_automatic_native_ball_lens_matches_dispersive_thick_lens_area_and_power() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![sphere(), plane(9, Point3::new(0.0, 0.0, 2.5), 0.05)])
                .unwrap();
        let lens = DielectricGlass::new(
            CauchyIor::try_new(1.45, 0.012, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.15).unwrap(),
            GlassProvenance::Custom,
        );
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            bindings(Vec3::new(0.0, 1.0, 0.0)),
            100_000,
        )
        .unwrap();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, -3.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        for wavelength in [400.0_f64, 550.0, 700.0] {
            let input = scene
                .launch_bundle(
                    ray,
                    Vec3::new(1.0, 0.0, 0.0),
                    wavelength,
                    Stokes::unpolarized(1.0),
                    1,
                )
                .unwrap();
            let mut branches = 0;
            let result = scene
                .trace_bundle_with_branches(cx, input, 1, settings(), |contact, context| {
                    assert_eq!(contact.object_id, 1);
                    assert!(!context.total_internal_reflection);
                    branches += 1;
                    DielectricEvent::Transmission
                })
                .unwrap();
            assert_eq!(branches, 2);
            assert_eq!(result.termination(), PathTermination::Detected);
            assert_eq!(result.surface_events(), 3);
            assert_eq!(result.medium_index(), 1);
            let path = result.bundle();
            let geometric = path.propagation();
            let propagation = geometric.propagation();
            assert_eq!(geometric.initial_ray(), ray);
            assert_eq!(geometric.interface_events(), 2);
            let index = 1.45 + 0.012 / (wavelength * 1e-3).powi(2);
            let group = 1.45 + 0.036 / (wavelength * 1e-3).powi(2);
            // Exact first variation of the finite-thickness ball lens:
            // 2 m air + two R=1 interfaces + 2 m glass + 1.5 m air.
            let b = 15.0 / index - 9.5;
            let normal = Vec3::new(0.0, 0.0, -1.0);
            let projection = geometric
                .project_at_endpoint(normal, ProjectionOptions::default())
                .unwrap();
            close(projection.area_per_solid_angle_m2(), b * b, 2e-6);
            let columns = projection.angular_derivatives();
            close_vec(columns[0].position_m_per_rad, Vec3::new(b, 0.0, 0.0), 2e-6);
            close_vec(columns[1].position_m_per_rad, Vec3::new(0.0, b, 0.0), 2e-6);
            close(propagation.ray.origin.z, 2.5, 2e-7);
            close(propagation.distance_m, 5.5, 3e-7);
            close(propagation.optical_path_length_m, 3.5 + 2.0 * index, 5e-7);
            close(
                propagation.group_optical_path_length_m.unwrap(),
                3.5 + 2.0 * group,
                5e-7,
            );
            close(propagation.optical_depth, 0.3, 1e-7);
            let transmission = 4.0 * index / (1.0 + index).powi(2);
            close(path.interface_log_power().exp(), transmission.powi(2), 1e-7);
            let irradiance = 3.0 * transmission.powi(2) * (-0.3_f64).exp() / (b * b);
            close(
                path.point_source_irradiance(normal, ProjectionOptions::default(), 3.0)
                    .unwrap(),
                irradiance,
                3e-5 * irradiance,
            );
            assert!(result.max_contact_resolution_m() > 0.0);
        }
    });
}

#[test]
fn g1_automatic_off_axis_ball_matches_independent_snell_footprint_and_detector_contact() {
    with_cx(|cx, _| {
        let ray = Ray {
            origin: Point3::new(0.12, 0.07, -3.0),
            dir: unit(Vec3::new(0.025, -0.015, 1.0)),
        };
        let (expected, outside, inside, transfer) = ball_reference(ray, 1.5);
        let expected_columns = reference_columns(ray, |r| ball_reference(r, 1.5).0.origin);
        let expected_area = cross(expected_columns[0], expected_columns[1]).z.abs();
        let geometry =
            InstanceScene::try_new(vec![sphere(), plane(9, expected.origin, 0.015)]).unwrap();
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let mut previous_bundle = None;
        // The detector reference controls its polarization readout. It must
        // neither rotate the original launch-angle chart nor erase the native
        // contact normal/identities used by irradiance and connection consumers.
        for detector_axis in [Vec3::new(0.0, 1.0, 0.0), Vec3::new(1.0, 0.0, 0.0)] {
            let scene = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                bindings(detector_axis),
                100_000,
            )
            .unwrap();
            let input = scene
                .launch_bundle(
                    ray,
                    launch_basis(ray.dir)[0],
                    550.0,
                    Stokes::unpolarized(1.0),
                    1,
                )
                .unwrap();
            let result = scene
                .trace_bundle_with_branches(cx, input, 1, settings(), |_, _| {
                    DielectricEvent::Transmission
                })
                .unwrap();
            assert_eq!(result.termination(), PathTermination::Detected);
            assert_eq!(result.detector_reference_axis(), Some(detector_axis));
            let contact = result.last_contact().unwrap();
            assert_eq!(contact.object_id, 9);
            assert_eq!(
                contact.geometry_identity,
                geometry.instances()[1].geometry_identity()
            );
            assert_eq!(
                contact.frame_identity,
                geometry.instances()[1].frame_identity()
            );
            close_vec(contact.witness.normal, Vec3::new(0.0, 0.0, 1.0), 1e-14);
            let path = result.bundle();
            let geometric = path.propagation();
            let propagation = geometric.propagation();
            close_vec(
                propagation.ray.origin.delta_from(expected.origin),
                Vec3::new(0.0, 0.0, 0.0),
                2e-6,
            );
            close_vec(propagation.ray.dir, expected.dir, 2e-6);
            let displacement = propagation
                .ray
                .origin
                .delta_from(contact.witness.surface_point);
            assert!(displacement.dot(displacement).sqrt() <= contact.resolution_m + 1e-13);
            let projection = geometric
                .project_at_endpoint(contact.witness.normal, ProjectionOptions::default())
                .unwrap();
            for (actual, expected) in projection
                .angular_derivatives()
                .into_iter()
                .zip(expected_columns)
            {
                close_vec(actual.position_m_per_rad, expected, 3e-6);
            }
            close(projection.area_per_solid_angle_m2(), expected_area, 3e-6);
            close(propagation.distance_m, outside + inside, 3e-6);
            close(
                propagation.optical_path_length_m,
                outside + 1.5 * inside,
                4e-6,
            );
            close(path.interface_log_power().exp(), transfer, 2e-6);
            let expected_irradiance = transfer * (-0.2 * inside).exp() / expected_area;
            close(
                path.point_source_irradiance(
                    contact.witness.normal,
                    ProjectionOptions::default(),
                    1.0,
                )
                .unwrap(),
                expected_irradiance,
                5e-5 * expected_irradiance,
            );
            if let Some(previous) = previous_bundle {
                assert_eq!(*path, previous);
            }
            previous_bundle = Some(*path);
        }
        // This is a genuinely finite receiver; the original straight ray misses.
        assert!(
            geometry.instances()[1]
                .intersect(cx, &ray, 8.0, 1e-10)
                .unwrap()
                .is_none()
        );
    });
}

/// Analytic n=1.4+0.2*x ray in the z parameter, followed by an independently
/// located sphere/plane crossing and Snell transmission to z=2. The formula
/// remains valid under both launch-angle perturbations used below.
pub(super) fn grin_reference(ray: Ray, curved: bool) -> (Ray, f64, f64, f64, f64) {
    let (n0, gradient) = (1.4_f64, 0.2_f64);
    let py = n0 * ray.dir.y;
    let pz = n0 * ray.dir.z;
    let perpendicular = py.hypot(pz);
    let u0 = (n0 * ray.dir.x / perpendicular).asinh();
    let at_z = |z: f64| {
        let u = u0 + gradient * z / pz;
        let n = perpendicular * u.cosh();
        let point = Point3::new((n - n0) / gradient, py * z / pz, z);
        let direction = Vec3::new(u.tanh(), py / n, pz / n);
        let distance = perpendicular * (u.sinh() - u0.sinh()) / gradient;
        let optical = perpendicular.powi(2) / (2.0 * gradient)
            * (u - u0 + u.sinh() * u.cosh() - u0.sinh() * u0.cosh());
        (point, direction, n, distance, optical)
    };
    let z = if curved {
        let (mut lower, mut upper) = (0.0, 1.0);
        for _ in 0..60 {
            let z = 0.5 * (lower + upper);
            let p = at_z(z).0;
            if p.x * p.x + p.y * p.y + p.z * p.z > 1.0 {
                upper = z;
            } else {
                lower = z;
            }
        }
        0.5 * (lower + upper)
    } else {
        1.0
    };
    let (hit, incoming, index, inside, optical) = at_z(z);
    let normal = if curved {
        unit(Vec3::new(hit.x, hit.y, hit.z))
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let outgoing = snell(incoming, normal, index);
    let after = (2.0 - hit.z) / outgoing.z;
    let ci = incoming.dot(normal);
    let ct = outgoing.dot(normal);
    let rs = ((index * ci - ct) / (index * ci + ct)).powi(2);
    let rp = ((ci - index * ct) / (ci + index * ct)).powi(2);
    (
        Ray {
            origin: hit.offset(outgoing.scale(after)),
            dir: outgoing,
        },
        inside,
        after,
        optical + after,
        1.0 - 0.5 * (rs + rp),
    )
}

#[test]
fn g1_automatic_grin_crossings_match_analytic_sphere_and_plane_footprints() {
    with_cx(|cx, _| {
        let material = glass(1.4, 0.3);
        let field = AffineIndexField::try_new(
            &material,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.2, 0.0, 0.0),
            Aabb {
                min: Point3::new(-1.5, -1.5, -2.0),
                max: Point3::new(1.5, 1.5, 3.0),
            },
        )
        .unwrap();
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Differentiable {
                field: &field,
                identity: ContentHash([39; 32]),
            },
            SceneMedium::Homogeneous(&air),
        ];
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        for curved in [false, true] {
            let (expected, inside, outside, optical, transfer) = grin_reference(ray, curved);
            let columns = reference_columns(ray, |r| grin_reference(r, curved).0.origin);
            let area = cross(columns[0], columns[1]).z.abs();
            let boundary = if curved {
                sphere()
            } else {
                plane(1, Point3::new(0.0, 0.0, 1.0), 5.0)
            };
            let geometry =
                InstanceScene::try_new(vec![boundary, plane(9, expected.origin, 0.02)]).unwrap();
            let scene = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                bindings(Vec3::new(0.0, 1.0, 0.0)),
                100_000,
            )
            .unwrap();
            let input = scene
                .launch_bundle(
                    ray,
                    Vec3::new(1.0, 0.0, 0.0),
                    550.0,
                    Stokes::unpolarized(1.0),
                    0,
                )
                .unwrap();
            let result = scene
                .trace_bundle_with_branches(cx, input, 0, settings(), |_, _| {
                    DielectricEvent::Transmission
                })
                .unwrap();
            assert_eq!(result.termination(), PathTermination::Detected);
            assert_eq!(result.surface_events(), 2);
            assert_eq!(result.medium_index(), 1);
            let path = result.bundle();
            let geometric = path.propagation();
            let propagation = geometric.propagation();
            close_vec(
                propagation.ray.origin.delta_from(expected.origin),
                Vec3::new(0.0, 0.0, 0.0),
                2e-6,
            );
            close_vec(propagation.ray.dir, expected.dir, 2e-6);
            close(propagation.distance_m, inside + outside, 3e-6);
            close(propagation.optical_path_length_m, optical, 4e-6);
            close(
                propagation.group_optical_path_length_m.unwrap(),
                optical,
                4e-6,
            );
            close(propagation.optical_depth, 0.3 * inside, 1e-6);
            close(path.interface_log_power().exp(), transfer, 1e-6);
            let normal = Vec3::new(0.0, 0.0, -1.0);
            let projection = geometric
                .project_at_endpoint(normal, ProjectionOptions::default())
                .unwrap();
            for (actual, expected) in projection.angular_derivatives().into_iter().zip(columns) {
                close_vec(actual.position_m_per_rad, expected, 5e-6);
            }
            close(projection.area_per_solid_angle_m2(), area, 2e-5);
            let irradiance = transfer * (-0.3 * inside).exp() / area;
            close(
                path.point_source_irradiance(normal, ProjectionOptions::default(), 1.0)
                    .unwrap(),
                irradiance,
                1e-5 * irradiance,
            );
            assert!(
                geometry.instances()[1]
                    .intersect(cx, &ray, 4.0, 1e-10)
                    .unwrap()
                    .is_none()
            );
            if !curved {
                // Erasing the explicit Hessian capability must refuse, even
                // when the underlying analytic field happens to be affine.
                let first_order_media = [
                    SceneMedium::GradientIndex {
                        field: &field,
                        identity: ContentHash([39; 32]),
                    },
                    SceneMedium::Homogeneous(&air),
                ];
                let first_order_scene = CurvedPolarizedScene::try_new(
                    cx,
                    &geometry,
                    &first_order_media,
                    bindings(Vec3::new(0.0, 1.0, 0.0)),
                    100_000,
                )
                .unwrap();
                assert!(matches!(
                    first_order_scene.launch_bundle(
                        ray,
                        Vec3::new(1.0, 0.0, 0.0),
                        550.0,
                        Stokes::unpolarized(1.0),
                        0,
                    ),
                    Err(CurvedBundleSceneError::MissingSecondOrderField { medium: 0 })
                ));
            }
        }
    });
}

#[test]
fn g1_automatic_native_tir_keeps_phase_and_refuses_forbidden_transmission() {
    with_cx(|cx, _| {
        let after = 0.2 / 0.28;
        let receiver = Point3::new(0.8 - 0.96 * after, 0.0, 0.8);
        let geometry = InstanceScene::try_new(vec![sphere(), plane(9, receiver, 0.04)]).unwrap();
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            bindings(Vec3::new(0.0, 1.0, 0.0)),
            100_000,
        )
        .unwrap();
        let ray = Ray {
            origin: Point3::new(0.8, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let input = scene
            .launch_bundle(
                ray,
                Vec3::new(0.0, 1.0, 0.0),
                550.0,
                Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(),
                0,
            )
            .unwrap();
        let result = scene
            .trace_bundle_with_branches(cx, input, 0, settings(), |_, context| {
                assert!(context.total_internal_reflection);
                DielectricEvent::Reflection
            })
            .unwrap();
        assert_eq!(result.termination(), PathTermination::Detected);
        assert_eq!(result.medium_index(), 0);
        assert_eq!(result.surface_events(), 2);
        let path = result.bundle();
        let propagation = path.propagation().propagation();
        close_vec(
            propagation.ray.origin.delta_from(receiver),
            Vec3::new(0.0, 0.0, 0.0),
            4e-6,
        );
        close_vec(propagation.ray.dir, Vec3::new(-0.96, 0.0, 0.28), 2e-6);
        close(propagation.distance_m, 0.6 + after, 5e-6);
        close(propagation.optical_path_length_m, 1.5 * (0.6 + after), 8e-6);
        close(propagation.optical_depth, 0.2 * (0.6 + after), 2e-6);
        close(path.interface_log_power(), 0.0, 1e-13);
        let gamma = (1.5_f64.powi(2) * 0.8_f64.powi(2) - 1.0).sqrt();
        let phase = 2.0 * ((1.5 * gamma).atan2(0.6) - gamma.atan2(1.5 * 0.6));
        close(path.interface_state().u, phase.cos(), 3e-6);
        close(path.interface_state().v, phase.sin(), 3e-6);
        assert!(matches!(
            scene.trace_bundle_with_branches(cx, input, 0, settings(), |_, _| {
                DielectricEvent::Transmission
            }),
            Err(CurvedBundleSceneError::Bundle(
                PolarizedBundleError::Interface(BundleInterfaceError::TotalInternalReflection)
            ))
        ));
    });
}

#[test]
fn g0_automatic_scene_preserves_occlusion_and_medium_topology() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            sphere(),
            plane(5, Point3::new(0.0, 0.0, 1.5), 0.2),
            plane(9, Point3::new(0.0, 0.0, 2.5), 0.1),
        ])
        .unwrap();
        let lens = glass(1.5, 0.0);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let mut surfaces = bindings(Vec3::new(0.0, 1.0, 0.0));
        surfaces.push(SurfaceBinding {
            object_id: 5,
            interaction: Interaction::Opaque,
        });
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, surfaces, 100_000).unwrap();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, -3.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let launch = |medium| {
            scene
                .launch_bundle(
                    ray,
                    Vec3::new(1.0, 0.0, 0.0),
                    550.0,
                    Stokes::unpolarized(1.0),
                    medium,
                )
                .unwrap()
        };
        let blocked = scene
            .trace_bundle_with_branches(cx, launch(1), 1, settings(), |_, _| {
                DielectricEvent::Transmission
            })
            .unwrap();
        assert_eq!(blocked.termination(), PathTermination::Opaque);
        assert_eq!(blocked.last_contact().unwrap().object_id, 5);
        assert_eq!(blocked.surface_events(), 3);
        assert_eq!(blocked.detector_reference_axis(), None);
        close(
            blocked.bundle().propagation().propagation().ray.origin.z,
            1.5,
            2e-7,
        );
        assert!(matches!(
            scene.trace_bundle_with_branches(cx, launch(0), 0, settings(), |_, _| panic!(
                "medium check must precede branch selection"
            )),
            Err(CurvedBundleSceneError::Scene(CurvedSceneError::Scene(
                SceneError::MediumMismatch {
                    object_id: 1,
                    expected: 1,
                    actual: 0
                }
            )))
        ));
        assert!(
            scene
                .trace_bundle_with_branches(cx, launch(1), 0, settings(), |_, _| {
                    DielectricEvent::Transmission
                })
                .is_err()
        );
    });
}

#[test]
fn g4_automatic_scene_accumulates_step_and_query_budgets_and_cancels_after_interface() {
    with_cx(|cx, gate| {
        let geometry =
            InstanceScene::try_new(vec![sphere(), plane(9, Point3::new(0.0, 0.0, 2.5), 0.1)])
                .unwrap();
        let lens = glass(1.5, 0.1);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            bindings(Vec3::new(0.0, 1.0, 0.0)),
            100_000,
        )
        .unwrap();
        let input = scene
            .launch_bundle(
                Ray {
                    origin: Point3::new(0.0, 0.0, -3.0),
                    dir: Vec3::new(0.0, 0.0, 1.0),
                },
                Vec3::new(1.0, 0.0, 0.0),
                550.0,
                Stokes::unpolarized(1.0),
                1,
            )
            .unwrap();
        let complete = scene
            .trace_bundle_with_branches(cx, input, 1, settings(), |_, _| {
                DielectricEvent::Transmission
            })
            .unwrap();
        let propagation = complete.bundle().propagation().propagation();
        let total_steps = propagation.accepted_steps + propagation.rejected_steps;
        assert!(total_steps > 3);
        let mut limited = settings();
        limited.trace.max_total_steps = total_steps - 1;
        let mut interfaces = 0;
        assert!(
            scene
                .trace_bundle_with_branches(cx, input, 1, limited, |_, _| {
                    interfaces += 1;
                    DielectricEvent::Transmission
                })
                .is_err()
        );
        assert_eq!(
            interfaces, 2,
            "the cumulative step limit must survive both interfaces"
        );
        limited = settings();
        limited.trace.max_query_work = complete.query_work() - 1;
        assert!(
            scene
                .trace_bundle_with_branches(cx, input, 1, limited, |_, _| {
                    DielectricEvent::Transmission
                })
                .is_err()
        );
        limited = settings();
        limited.trace.max_distance_m = 3.0;
        let bounded = scene
            .trace_bundle_with_branches(cx, input, 1, limited, |_, _| DielectricEvent::Transmission)
            .unwrap();
        assert_eq!(bounded.termination(), PathTermination::DistanceLimit);
        assert_eq!(bounded.surface_events(), 1);
        close(
            bounded.bundle().propagation().propagation().distance_m,
            3.0,
            1e-12,
        );
        assert!(bounded.bundle().interface_log_power().is_finite());
        limited = settings();
        limited.trace.max_surface_events = 2;
        let bounded = scene
            .trace_bundle_with_branches(cx, input, 1, limited, |_, _| DielectricEvent::Transmission)
            .unwrap();
        assert_eq!(bounded.termination(), PathTermination::SurfaceLimit);
        assert_eq!(bounded.surface_events(), 2);
        assert_eq!(bounded.bundle().propagation().interface_events(), 2);
        let cancelled = scene.trace_bundle_with_branches(cx, input, 1, settings(), |_, _| {
            gate.request();
            DielectricEvent::Transmission
        });
        assert!(matches!(
            cancelled,
            Err(CurvedBundleSceneError::Scene(CurvedSceneError::Scene(
                SceneError::Cancelled
            )))
        ));
        assert_eq!(input.propagation().propagation().distance_m, 0.0);
        assert_eq!(input.interface_state(), Stokes::unpolarized(1.0));
    });
}

/// One uniform medium and a small detector far up the z axis, which every
/// launch below misses, so each trace escapes.
fn far_detector_scene_parts() -> (InstanceScene, DielectricGlass) {
    let geometry =
        InstanceScene::try_new(vec![plane(9, Point3::new(0.0, 0.0, 50.0), 0.1)]).unwrap();
    (geometry, glass(1.0, 0.0))
}

fn far_detector_binding() -> Vec<SurfaceBinding> {
    vec![SurfaceBinding {
        object_id: 9,
        interaction: Interaction::Detector {
            reference_axis: Vec3::new(0.0, 1.0, 0.0),
        },
    }]
}

#[test]
fn g1_uniform_launches_renormalized_in_their_last_bit_trace_as_original_launches() {
    // In a uniform medium a zero-distance launch carries initial_ray().dir = unit(d)
    // and propagation().ray.dir = unit(unit(d)), and the two can differ in the last
    // bit. Before the bounded check, every such launch refused with
    // InvalidInitialBundle. A Mac Python replay of unit() predicts 14 of these 343
    // directions differ. The test counts them on the engine's own launches and
    // requires at least one, so the scan exercises the renormalized case.
    with_cx(|cx, _| {
        let (geometry, air) = far_detector_scene_parts();
        let media = [SceneMedium::Homogeneous(&air)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, far_detector_binding(), 100_000)
                .unwrap();
        let mut renormalized = 0;
        for i in -3..=3 {
            for j in -3..=3 {
                for k in -3..=3 {
                    let ray = Ray {
                        origin: Point3::new(0.0, 0.0, 0.0),
                        dir: Vec3::new(f64::from(i) + 0.1, f64::from(j) + 0.2, f64::from(k) + 0.3),
                    };
                    let launch = scene
                        .launch_bundle(
                            ray,
                            Vec3::new(1.0, 0.0, 0.0),
                            550.0,
                            Stokes::unpolarized(1.0),
                            0,
                        )
                        .unwrap();
                    let propagation = launch.propagation();
                    if propagation.initial_ray().dir != propagation.propagation().ray.dir {
                        renormalized += 1;
                    }
                    let path = scene
                        .trace_bundle_with_branches(cx, launch, 0, settings(), |_, _| {
                            DielectricEvent::Transmission
                        })
                        .unwrap_or_else(|error| panic!("direction {:?}: {error:?}", ray.dir));
                    assert_eq!(path.termination(), PathTermination::Escaped);
                    assert_eq!(path.surface_events(), 0);
                }
            }
        }
        assert!(
            renormalized > 0,
            "the scan must contain launches whose second normalization moved a bit"
        );
    });
}

#[test]
fn g0_continued_bundle_is_not_an_original_launch() {
    with_cx(|cx, _| {
        let (geometry, air) = far_detector_scene_parts();
        let media = [SceneMedium::Homogeneous(&air)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, far_detector_binding(), 100_000)
                .unwrap();
        let launch = scene
            .launch_bundle(
                Ray {
                    origin: Point3::new(0.0, 0.0, 0.0),
                    dir: Vec3::new(0.1, 0.2, 1.0),
                },
                Vec3::new(1.0, 0.0, 0.0),
                550.0,
                Stokes::unpolarized(1.0),
                0,
            )
            .unwrap();
        let continued = launch
            .continue_in(
                &crate::gradient_index::HomogeneousIndexField::new(&air),
                0.5,
                BundleOptions::default(),
            )
            .unwrap();
        assert_eq!(continued.propagation().interface_events(), 0);
        assert!(matches!(
            scene.trace_bundle_with_branches(cx, continued, 0, settings(), |_, _| {
                DielectricEvent::Transmission
            }),
            Err(CurvedBundleSceneError::InvalidInitialBundle)
        ));
        assert!(
            scene
                .trace_bundle_with_branches(cx, launch, 0, settings(), |_, _| {
                    DielectricEvent::Transmission
                })
                .is_ok()
        );
    });
}

#[test]
fn g0_renormalized_direction_admits_one_renormalization_and_refuses_seven_eps() {
    let bits =
        |x: u64, y: u64, z: u64| Vec3::new(f64::from_bits(x), f64::from_bits(y), f64::from_bits(z));
    // R-c printed these launches (initial, propagated): the molecular fixture ray
    // and the coating area-emission ray, each b = unit(a) exactly.
    let molecular = (
        bits(
            0xbfe3_c1aa_0689_dff8,
            0xbfbb_89c6_1262_166a,
            0x3fe8_efc1_a214_0d40,
        ),
        bits(
            0xbfe3_c1aa_0689_dff8,
            0xbfbb_89c6_1262_1669,
            0x3fe8_efc1_a214_0d3f,
        ),
    );
    let coating = (
        bits(0x3fe9_35fd_eba9_579c, 0, 0x3fe3_b547_bec7_7728),
        bits(0x3fe9_35fd_eba9_579d, 0, 0x3fe3_b547_bec7_7728),
    );
    for (a, b) in [molecular, coating] {
        assert!(renormalized_direction_matches(a, a));
        assert!(renormalized_direction_matches(a, b));
        assert!(renormalized_direction_matches(b, a));
    }
    let a = molecular.0;
    let moved = |relative: f64| Vec3::new(a.x, a.y + relative * f64::EPSILON * a.y.abs(), a.z);
    // The sum rounds by at most half an ulp of a.y, under eps/2 relative, so a
    // 5 eps move stays inside the 6 eps bound and a 7 eps move stays outside.
    assert!(renormalized_direction_matches(a, moved(5.0)));
    assert!(!renormalized_direction_matches(a, moved(7.0)));
    // A zero component must stay exactly zero, up to two subnormal units.
    let zero = coating.0;
    assert!(renormalized_direction_matches(
        zero,
        Vec3::new(zero.x, 2.0 * f64::from_bits(1), zero.z)
    ));
    assert!(!renormalized_direction_matches(
        zero,
        Vec3::new(zero.x, 1e-300, zero.z)
    ));
    assert!(!renormalized_direction_matches(
        a,
        Vec3::new(f64::NAN, a.y, a.z)
    ));
}
