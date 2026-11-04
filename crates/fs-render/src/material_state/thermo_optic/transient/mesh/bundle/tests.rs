use super::*;
use crate::dielectric::{BeerLambertAbsorption, DielectricGlass, GlassProvenance, SellmeierIor};
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::bundle::{AngularDerivative, BundleError, ProjectionOptions};
use crate::gradient_index::event::propagate_to_boundary;
use fs_blake3::ContentHash;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig, ScalarField,
    ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::thermo_optic::{SchottThermoOpticLaw, ThermoOpticCoefficients};
use fs_rep_mesh::TetComplex;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 96,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx);
    });
}
fn norm(v: Vec3) -> f64 {
    v.dot(v).sqrt()
}
fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / norm(v))
}
fn add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}
fn sub(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:.17e} != {b:.17e}; tolerance {t}");
}
fn vector_near(a: Vec3, b: Vec3, t: f64) {
    near(norm(sub(a, b)), 0.0, t);
}
fn point_near(a: Point3, b: Point3, t: f64) {
    vector_near(a.delta_from(b), Vec3::new(0.0, 0.0, 0.0), t);
}
fn reference(dispersive: bool) -> DielectricGlass {
    DielectricGlass::new(
        SellmeierIor::try_new_in_band(
            [1.25, 0.0, 0.0],
            [if dispersive { 0.04 } else { 0.0 }, 0.0, 0.0],
            (1000.0, 2000.0),
        )
        .unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}
fn law(nonlinear: bool) -> SchottThermoOpticLaw {
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [0.003, if nonlinear { 2.0e-5 } else { 0.0 }, 0.0],
            e: [0.0; 2],
            resonance_um: 0.0,
        },
        300.0,
        [200.0, 450.0],
        [1000.0, 2000.0],
        ContentHash([96; 32]),
    )
    .unwrap()
}
fn options() -> P1MeshBundleOptions {
    P1MeshBundleOptions {
        mesh: P1MeshOptions {
            boundary_id: 43,
            max_internal_faces: 128,
            event: GradientEventOptions {
                integration: GradientIndexOptions {
                    initial_step_m: 0.01,
                    max_step_m: 0.04,
                    min_step_m: 1.0e-13,
                    relative_tolerance: 2.0e-12,
                    position_tolerance_m: 2.0e-12,
                    direction_tolerance: 2.0e-12,
                    optical_path_tolerance_m: 2.0e-12,
                    optical_depth_tolerance: 2.0e-12,
                    ..GradientIndexOptions::default()
                },
                surface_tolerance_m: 1.0e-8,
                ..GradientEventOptions::default()
            },
        },
        position_derivative_tolerance_m: 2.0e-12,
        direction_derivative_tolerance: 2.0e-12,
    }
}
fn bipyramid(split: bool) -> ConductionMesh {
    let mut p = vec![
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 1.0],
    ];
    let tets = if split {
        p.extend([p[0], p[1], p[2]]);
        vec![[0, 1, 2, 3], [5, 6, 7, 4]]
    } else {
        vec![[0, 1, 2, 3], [0, 1, 2, 4]]
    };
    ConductionMesh::new(TetComplex::from_tets(p.len(), tets), p).unwrap()
}
fn column(reverse: bool) -> ConductionMesh {
    let mut p = Vec::new();
    for z in [-1.0, 0.0, 1.0] {
        p.extend([
            [-1.0, -1.0, z],
            [1.0, -1.0, z],
            [-1.0, 1.0, z],
            [1.0, 1.0, z],
        ]);
    }
    let mut tets = Vec::new();
    for base in [0, 4] {
        for a in [
            [1, 2, 4],
            [1, 4, 2],
            [2, 1, 4],
            [2, 4, 1],
            [4, 1, 2],
            [4, 2, 1],
        ] {
            let mut t = [base, base + a[0], base + a[0] + a[1], base + 7];
            if reverse {
                t.swap(0, 1);
            }
            tets.push(t);
        }
    }
    ConductionMesh::new(TetComplex::from_tets(p.len(), tets), p).unwrap()
}
fn locate(snapshot: P1ConductionTemperatureSnapshot<'_>, p: Point3) -> P1ConductionElement<'_> {
    // Fixture setup only; the production caller must supply its original cell.
    (0..snapshot.mesh().element_count())
        .map(|i| snapshot.element(i).unwrap())
        .find(|e| e.barycentric(p).is_ok())
        .unwrap()
}
fn trace(
    cx: &Cx<'_>,
    mesh: &P1ConductionMesh<'_>,
    ray: Ray,
    length: f64,
    opts: P1MeshBundleOptions,
    nonlinear: bool,
    dispersive: bool,
) -> Result<P1MeshBundlePropagation, P1MeshError> {
    let e = locate(mesh.snapshot, ray.origin);
    let glass = reference(dispersive);
    let law = law(nonlinear);
    let binding =
        MaterialThermoOpticField::try_new(&glass, &law, &e, ContentHash([97; 32])).unwrap();
    mesh.propagate_bundle(
        cx,
        &binding,
        ray,
        Vec3::new(1.0, 0.0, 0.0),
        1500.0,
        length,
        opts,
    )
}
fn primal(
    cx: &Cx<'_>,
    mesh: &P1ConductionMesh<'_>,
    ray: Ray,
    length: f64,
    opts: P1MeshBundleOptions,
    nonlinear: bool,
    dispersive: bool,
) -> P1MeshPropagation {
    let e = locate(mesh.snapshot, ray.origin);
    let glass = reference(dispersive);
    let law = law(nonlinear);
    let binding =
        MaterialThermoOpticField::try_new(&glass, &law, &e, ContentHash([97; 32])).unwrap();
    mesh.propagate(
        cx,
        &binding,
        ray,
        Vec3::new(1.0, 0.0, 0.0),
        1500.0,
        length,
        opts.mesh,
    )
    .unwrap()
}
fn launch() -> Ray {
    Ray {
        origin: Point3::new(-0.06, -0.04, -0.6),
        dir: unit(Vec3::new(0.18, 0.07, 1.0)),
    }
}
fn temperature(p: &[f64; 3]) -> f64 {
    300.0 + 40.0 * p[0] - 16.0 * p[1]
        + if p[2] < 0.0 {
            60.0 * p[2]
        } else {
            -50.0 * p[2]
        }
}

// Independent exact ray in n(r)=n0+g.(r-r0): d(nu)/ds=g. This
// oracle has no native RK stage, event collar or thermal reconstruction.
fn affine(ray: Ray, n0: f64, g: Vec3, s: f64) -> (Ray, f64) {
    let b = norm(g);
    if b == 0.0 {
        return (
            Ray {
                origin: ray.origin.offset(ray.dir.scale(s)),
                dir: ray.dir,
            },
            n0 * s,
        );
    }
    let axis = g.scale(1.0 / b);
    let q0 = n0 * ray.dir.dot(axis);
    let perpendicular = sub(ray.dir.scale(n0), axis.scale(q0));
    let a = norm(perpendicular);
    let q = q0 + b * s;
    let n = a.hypot(q);
    let integral = ((q / a).asinh() - (q0 / a).asinh()) / b;
    let along = s * (2.0 * q0 + b * s) / (n + n0);
    let point = ray
        .origin
        .offset(add(perpendicular.scale(integral), axis.scale(along)));
    let direction = add(perpendicular, axis.scale(q)).scale(1.0 / n);
    let primitive = |q: f64| 0.5 * (q * a.hypot(q) + a * a * (q / a).asinh());
    (
        Ray {
            origin: point,
            dir: direction,
        },
        (primitive(q) - primitive(q0)) / b,
    )
}
fn crossing(ray: Ray) -> (Ray, f64, f64) {
    let n0 = 1.5 + 0.00125 * (temperature(&[ray.origin.x, ray.origin.y, ray.origin.z]) - 300.0);
    let g = Vec3::new(0.05, -0.02, 0.075);
    let mut a = 0.0;
    let mut b = 1.0;
    for _ in 0..60 {
        let s = 0.5 * (a + b);
        if affine(ray, n0, g, s).0.origin.z < 0.0 {
            a = s;
        } else {
            b = s;
        }
    }
    let s = 0.5 * (a + b);
    let (mut hit, phase) = affine(ray, n0, g, s);
    hit.origin.z = 0.0;
    (hit, s, phase)
}
fn exact(ray: Ray, s: f64) -> (Ray, f64) {
    let (hit, at, phase) = crossing(ray);
    if s <= at {
        let n0 = 1.5 + 0.00125 * (temperature(&[ray.origin.x, ray.origin.y, ray.origin.z]) - 300.0);
        return affine(ray, n0, Vec3::new(0.05, -0.02, 0.075), s);
    }
    let n = 1.5 + 0.05 * hit.origin.x - 0.02 * hit.origin.y;
    let (end, tail) = affine(hit, n, Vec3::new(0.05, -0.02, -0.0625), s - at);
    (end, phase + tail)
}
fn finite_difference(
    ray: Ray,
    basis: [Vec3; 2],
    epsilon: f64,
    f: impl Fn(Ray) -> Ray,
) -> [AngularDerivative; 2] {
    basis.map(|axis| {
        let a = f(Ray {
            dir: unit(add(ray.dir, axis.scale(epsilon))),
            ..ray
        });
        let b = f(Ray {
            dir: unit(add(ray.dir, axis.scale(-epsilon))),
            ..ray
        });
        AngularDerivative {
            position_m_per_rad: a.origin.delta_from(b.origin).scale(0.5 / epsilon),
            direction_per_rad: sub(a.dir, b.dir).scale(0.5 / epsilon),
        }
    })
}
fn error(a: [AngularDerivative; 2], b: [AngularDerivative; 2]) -> f64 {
    a.into_iter()
        .zip(b)
        .map(|(a, b)| {
            norm(sub(a.position_m_per_rad, b.position_m_per_rad))
                .max(norm(sub(a.direction_per_rad, b.direction_per_rad)))
        })
        .fold(0.0, f64::max)
}
fn optics(result: &P1MeshBundlePropagation, expected: Ray, phase: f64, length: f64, tol: f64) {
    let p = result.bundle.propagation();
    point_near(p.ray.origin, expected.origin, tol);
    vector_near(p.ray.dir, expected.dir, tol);
    near(p.distance_m, length, 5.0e-15);
    near(p.optical_path_length_m, phase, tol);
    near(p.group_optical_path_length_m.unwrap(), phase, tol);
    near(p.optical_depth, 0.2 * length, 3.0e-14);
    near(p.power_transmittance(), (-0.2 * length).exp(), 3.0e-14);
    assert_eq!(result.bundle.interface_events(), 0);
}

#[test]
fn original_launch_family_and_footprint_cross_many_cells_in_both_orientations() {
    with_cx(|_, cx| {
        for reverse in [false, true] {
            for slope in [0.0, 18.0] {
                let mesh = column(reverse);
                let values: Vec<_> = mesh
                    .positions()
                    .iter()
                    .map(|p| 300.0 + slope * p[0])
                    .collect();
                let snapshot =
                    P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
                let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
                let ray = Ray {
                    origin: Point3::new(0.17, -0.23, -0.9),
                    dir: Vec3::new(0.0, 0.0, 1.0),
                };
                let length = 1.7;
                let result = trace(cx, &prepared, ray, length, options(), false, false).unwrap();
                let (expected, phase) = affine(
                    ray,
                    1.5 + 0.00125 * slope * ray.origin.x,
                    Vec3::new(0.00125 * slope, 0.0, 0.0),
                    length,
                );
                optics(&result, expected, phase, length, 4.0e-9);
                assert!(result.internal_faces >= 5);
                assert_eq!(result.bundle.initial_ray(), ray);
                let expected = finite_difference(ray, result.bundle.launch_basis(), 2.0e-5, |r| {
                    affine(
                        r,
                        1.5 + 0.00125 * slope * r.origin.x,
                        Vec3::new(0.00125 * slope, 0.0, 0.0),
                        length,
                    )
                    .0
                });
                assert!(error(result.bundle.angular_derivatives(), expected) < 3.0e-8);
                if slope == 0.0 {
                    let footprint = result
                        .bundle
                        .project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())
                        .unwrap();
                    near(
                        footprint.area_per_solid_angle_m2(),
                        length * length,
                        2.0e-12,
                    );
                    near(
                        footprint.point_source_irradiance(4.0).unwrap(),
                        4.0 * (-0.2 * length).exp() / (length * length),
                        2.0e-12,
                    );
                }
            }
        }
    });
}

#[test]
fn gradient_jump_matches_independent_equal_arc_rays_and_refines_three_times() {
    with_cx(|_, cx| {
        let mesh = bipyramid(false);
        let values: Vec<_> = mesh.positions().iter().map(temperature).collect();
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let ray = launch();
        let length = 1.0;
        let mut errors = Vec::new();
        for (step, tol) in [(0.2, 2.0e-5), (0.08, 2.0e-8), (0.03, 2.0e-11)] {
            let mut opts = options();
            let i = &mut opts.mesh.event.integration;
            i.initial_step_m = step;
            i.max_step_m = step;
            i.relative_tolerance = tol;
            i.position_tolerance_m = tol;
            i.direction_tolerance = tol;
            i.optical_path_tolerance_m = tol;
            i.optical_depth_tolerance = tol;
            opts.position_derivative_tolerance_m = tol;
            opts.direction_derivative_tolerance = tol;
            opts.mesh.event.surface_tolerance_m = (tol * 20.0).max(1.0e-8);
            let result = trace(cx, &prepared, ray, length, opts, false, false).unwrap();
            assert_eq!(result.internal_faces, 1);
            assert_eq!(result.element_index, 1);
            let expected = finite_difference(ray, result.bundle.launch_basis(), 1.0e-5, |r| {
                exact(r, length).0
            });
            let err = error(result.bundle.angular_derivatives(), expected);
            errors.push(err);
            let (central, phase) = exact(ray, length);
            optics(&result, central, phase, length, 2.0e-4_f64.max(20.0 * tol));
        }
        eprintln!("gradient-jump bundle errors: {errors:?}");
        assert!(errors[1] < errors[0] * 0.3, "{errors:?}");
        assert!(errors[2] < errors[1] * 0.3, "{errors:?}");
        assert!(errors[2] < 2.0e-8, "{errors:?}");
        // Equal-distance primal finite differences use the separate first-order
        // native engine, not the new bundle stage/coordinate implementation.
        let result = trace(cx, &prepared, ray, length, options(), false, false).unwrap();
        let expected = finite_difference(ray, result.bundle.launch_basis(), 2.0e-5, |r| {
            primal(cx, &prepared, r, length, options(), false, false)
                .propagation
                .ray
        });
        assert!(error(result.bundle.angular_derivatives(), expected) < 3.0e-7);
        let (central, phase) = exact(ray, length);
        optics(&result, central, phase, length, 5.0e-9);
    });
}

#[test]
fn equal_arc_columns_convert_on_both_sides_inside_the_finite_departure_collar() {
    with_cx(|_, cx| {
        let mesh = bipyramid(false);
        let values: Vec<_> = mesh.positions().iter().map(temperature).collect();
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let ray = launch();
        let e = snapshot.element(0).unwrap();
        let glass = reference(false);
        let law = law(false);
        let field =
            MaterialThermoOpticField::try_new(&glass, &law, &e, ContentHash([97; 32])).unwrap();
        let field = field.prepare(1500.0).unwrap();
        let mut opts = options();
        opts.mesh.event.surface_tolerance_m = 1.0e-4;
        let mut boundary = e
            .boundary_oracle(43, opts.mesh.event.surface_tolerance_m)
            .unwrap();
        let arrival = propagate_to_boundary(
            &field,
            &mut boundary,
            ray,
            Vec3::new(1.0, 0.0, 0.0),
            1500.0,
            1.0,
            None,
            opts.mesh.event,
        )
        .unwrap();
        let (_, at, _) = crossing(ray);
        let gap = at - arrival.propagation.distance_m;
        assert!(gap > 1.0e-6 && gap < 2.0e-4);
        let BoundaryTermination::Surface(witness) = arrival.termination else {
            panic!("the original finite face must be reached")
        };
        let next = prepared.neighbours[0][witness.feature_id].unwrap();
        let neighbour = snapshot.element(next.element).unwrap();
        let guarded = MeshBoundary::with_transition(
            neighbour
                .boundary_oracle(43, opts.mesh.event.surface_tolerance_m)
                .unwrap(),
            Transition::new(
                e,
                neighbour,
                next.local_face,
                arrival.propagation.ray.origin,
                witness,
            ),
        );
        // Admission concerns the actual integrated endpoint. Do not snap an
        // approximately reached ray to this exact mathematical plane merely
        // to create a test. Exercise the original exact-plane guard directly.
        assert_eq!(
            guarded.admit_bundle_endpoint(crossing(ray).0.origin),
            Err(GradientIndexError::BoundaryResolutionLimit),
        );
        for (sign, element) in [(-1.0, 0), (1.0, 1)] {
            let length = at + sign * 0.25 * gap;
            let result = trace(cx, &prepared, ray, length, opts, false, false).unwrap();
            assert_eq!(result.internal_faces, 1);
            assert_eq!(result.element_index, element);
            assert_eq!(result.termination, BoundaryTermination::DistanceLimit);
            let expected = finite_difference(ray, result.bundle.launch_basis(), 1.0e-6, |r| {
                exact(r, length).0
            });
            let err = error(result.bundle.angular_derivatives(), expected);
            assert!(err < 6.0e-8, "side {element}: {err}");
            let expected = finite_difference(ray, result.bundle.launch_basis(), 2.0e-6, |r| {
                primal(cx, &prepared, r, length, options(), false, false)
                    .propagation
                    .ray
            });
            assert!(error(result.bundle.angular_derivatives(), expected) < 2.0e-6);
            let (central, phase) = exact(ray, length);
            optics(&result, central, phase, length, 4.0e-10);
        }
    });
}

#[test]
fn nonlinear_material_hessian_and_distinct_group_path_survive_gradient_jumps() {
    with_cx(|_, cx| {
        let mesh = bipyramid(false);
        let values: Vec<_> = mesh.positions().iter().map(temperature).collect();
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let ray = launch();
        let length = 1.0;
        let result = trace(cx, &prepared, ray, length, options(), true, true).unwrap();
        let central = primal(cx, &prepared, ray, length, options(), true, true);
        point_near(
            result.bundle.propagation().ray.origin,
            central.propagation.ray.origin,
            5.0e-10,
        );
        vector_near(
            result.bundle.propagation().ray.dir,
            central.propagation.ray.dir,
            5.0e-10,
        );
        let p = result.bundle.propagation();
        near(
            p.optical_path_length_m,
            central.propagation.optical_path_length_m,
            5.0e-10,
        );
        near(
            p.group_optical_path_length_m.unwrap(),
            central.propagation.group_optical_path_length_m.unwrap(),
            5.0e-10,
        );
        assert!((p.optical_path_length_m - p.group_optical_path_length_m.unwrap()).abs() > 1.0e-3);
        near(p.optical_depth, 0.2 * length, 3.0e-14);
        let expected = finite_difference(ray, result.bundle.launch_basis(), 2.0e-5, |r| {
            primal(cx, &prepared, r, length, options(), true, true)
                .propagation
                .ray
        });
        assert!(error(result.bundle.angular_derivatives(), expected) < 3.0e-7);
    });
}

#[test]
fn accepted_backward_euler_snapshot_changes_the_native_bundle_and_footprint() {
    with_cx(|_, cx| {
        let mesh = bipyramid(false);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let time =
            BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(2.0).unwrap())
                .unwrap();
        let old = [280.0, 320.0, 300.0, 300.0, 300.0];
        let accepted = time
            .advance(
                cx,
                problem,
                None,
                &old,
                1.0,
                StepConfig {
                    linear: LinearConfig::default(),
                    energy_tolerance_j: 1.0e-10,
                },
            )
            .unwrap();
        for (&actual, expected) in accepted
            .temperature
            .iter()
            .zip([290.0, 310.0, 300.0, 300.0, 300.0])
        {
            near(actual, expected, 3.0e-11);
        }
        near(accepted.energy_residual_j, 0.0, 1.0e-10);
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, -0.5),
            dir: unit(Vec3::new(0.1, 0.0, 1.0)),
        };
        let length = 0.9;
        let mut bends = [0.0; 2];
        for (i, (values, amplitude)) in [(&old[..], 20.0), (&accepted.temperature[..], 10.0)]
            .into_iter()
            .enumerate()
        {
            let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, values).unwrap();
            let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
            let result = trace(cx, &prepared, ray, length, options(), false, false).unwrap();
            assert_eq!(result.internal_faces, 1);
            let g = Vec3::new(0.00125 * amplitude, 0.0, 0.0);
            let expected = finite_difference(ray, result.bundle.launch_basis(), 1.0e-5, |r| {
                affine(r, 1.5, g, length).0
            });
            assert!(error(result.bundle.angular_derivatives(), expected) < 3.0e-8);
            let (end, phase) = affine(ray, 1.5, g, length);
            optics(&result, end, phase, length, 3.0e-9);
            let footprint = result
                .bundle
                .project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())
                .unwrap();
            assert!(footprint.area_per_solid_angle_m2() > 0.8);
            bends[i] = result.bundle.propagation().ray.dir.x - ray.dir.x;
        }
        assert!(bends[0] > 0.01 && bends[1] > 0.49 * bends[0] && bends[1] < 0.51 * bends[0]);
    });
}

#[test]
fn contacts_cancellation_and_global_bundle_budgets_remain_explicit() {
    with_cx(|gate, cx| {
        let ray = launch();
        let mesh = bipyramid(false);
        let values: Vec<_> = mesh.positions().iter().map(temperature).collect();
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let result = trace(cx, &prepared, ray, 1.0, options(), false, false).unwrap();
        let mut opts = options();
        opts.mesh.max_internal_faces = 0;
        assert_eq!(
            trace(cx, &prepared, ray, 1.0, opts, false, false).unwrap_err(),
            P1MeshError::FaceBudgetExceeded
        );
        opts = options();
        let p = result.bundle.propagation();
        opts.mesh.event.integration.max_steps = p.accepted_steps + p.rejected_steps - 1;
        assert!(matches!(
            trace(cx, &prepared, ray, 1.0, opts, false, false),
            Err(P1MeshError::Propagation(
                GradientIndexError::StepBudgetExceeded
            )) | Err(P1MeshError::Bundle(BundleError::Propagation(
                GradientIndexError::StepBudgetExceeded
            )))
        ));
        opts = options();
        opts.mesh.event.max_boundary_queries = result.boundary_evaluations - 1;
        assert!(matches!(
            trace(cx, &prepared, ray, 1.0, opts, false, false),
            Err(P1MeshError::Propagation(
                GradientIndexError::BoundaryQueryBudgetExceeded
            )) | Err(P1MeshError::Bundle(BundleError::Propagation(
                GradientIndexError::BoundaryQueryBudgetExceeded
            )))
        ));
        let contact = bipyramid(true);
        let contact_values = vec![300.0; contact.vertex_count()];
        let contact_snapshot =
            P1ConductionTemperatureSnapshot::try_new(cx, &contact, &contact_values).unwrap();
        let contact_mesh = P1ConductionMesh::try_new(cx, contact_snapshot).unwrap();
        let stop = trace(cx, &contact_mesh, ray, 1.0, options(), false, false).unwrap();
        assert!(matches!(stop.termination, BoundaryTermination::Surface(_)));
        assert_eq!(stop.internal_faces, 0);
        gate.request();
        assert!(matches!(
            trace(cx, &prepared, ray, 1.0, options(), false, false),
            Err(P1MeshError::Propagation(GradientIndexError::Cancelled))
                | Err(P1MeshError::Bundle(BundleError::Propagation(
                    GradientIndexError::Cancelled
                )))
        ));
    });
}

#[test]
fn insufficient_transversality_and_foreign_snapshot_refuse() {
    with_cx(|_, cx| {
        let mesh = bipyramid(false);
        let values: Vec<_> = mesh.positions().iter().map(temperature).collect();
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let mut opts = options();
        opts.mesh.event.min_approach_cosine = 0.999;
        assert!(matches!(
            trace(cx, &prepared, launch(), 1.0, opts, false, false),
            Err(P1MeshError::Propagation(
                GradientIndexError::BoundaryResolutionLimit
            )) | Err(P1MeshError::Bundle(BundleError::Propagation(
                GradientIndexError::BoundaryResolutionLimit
            )))
        ));
        // Equal numerical values are not the accepted snapshot's identity.
        let copied = values.clone();
        let other = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &copied).unwrap();
        let element = other.element(0).unwrap();
        let glass = reference(false);
        let law = law(false);
        let binding =
            MaterialThermoOpticField::try_new(&glass, &law, &element, ContentHash([97; 32]))
                .unwrap();
        assert_eq!(
            prepared
                .propagate_bundle(
                    cx,
                    &binding,
                    launch(),
                    Vec3::new(1.0, 0.0, 0.0),
                    1500.0,
                    1.0,
                    options(),
                )
                .unwrap_err(),
            P1MeshError::ForeignSnapshot
        );
    });
}
