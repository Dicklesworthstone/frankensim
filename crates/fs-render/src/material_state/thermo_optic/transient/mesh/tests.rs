use super::*;
use crate::dielectric::{BeerLambertAbsorption, DielectricGlass, GlassProvenance, SellmeierIor};
use crate::gradient_index::GradientIndexOptions;
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
                seed: 94,
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
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance}"
    );
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn point_near(a: Point3, b: Point3, tolerance: f64) {
    vector_near(
        Vec3::new(a.x, a.y, a.z),
        Vec3::new(b.x, b.y, b.z),
        tolerance,
    );
}
fn reference() -> DielectricGlass {
    DielectricGlass::new(
        SellmeierIor::try_new_in_band([1.25, 0.0, 0.0], [0.0; 3], (1000.0, 2000.0)).unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}
fn law() -> SchottThermoOpticLaw {
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [0.003, 0.0, 0.0],
            e: [0.0; 2],
            resonance_um: 0.0,
        },
        300.0,
        [200.0, 450.0],
        [1000.0, 2000.0],
        ContentHash([94; 32]),
    )
    .unwrap()
}
fn options() -> P1MeshOptions {
    P1MeshOptions {
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
    }
}
fn column(reverse: bool) -> ConductionMesh {
    let mut positions = Vec::new();
    for z in [-1.0, 0.0, 1.0] {
        positions.extend([
            [-1.0, -1.0, z],
            [1.0, -1.0, z],
            [-1.0, 1.0, z],
            [1.0, 1.0, z],
        ]);
    }
    let mut tets = Vec::new();
    for base in [0, 4] {
        for axes in [
            [1, 2, 4],
            [1, 4, 2],
            [2, 1, 4],
            [2, 4, 1],
            [4, 1, 2],
            [4, 2, 1],
        ] {
            let mut tet = [base, base + axes[0], base + axes[0] + axes[1], base + 7];
            if reverse {
                tet.swap(0, 1);
            }
            tets.push(tet);
        }
    }
    ConductionMesh::new(TetComplex::from_tets(12, tets), positions).unwrap()
}
fn bipyramid(height: f64, split_contact: bool) -> ConductionMesh {
    let mut positions = vec![
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.0, height],
    ];
    let tets = if split_contact {
        positions.extend([[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0]]);
        vec![[0, 1, 2, 3], [5, 6, 7, 4]]
    } else {
        vec![[0, 1, 2, 3], [0, 1, 2, 4]]
    };
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}
fn launch() -> Ray {
    Ray {
        origin: Point3::new(0.17, -0.23, -0.9),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn locate(snapshot: P1ConductionTemperatureSnapshot<'_>, p: Point3) -> P1ConductionElement<'_> {
    // Test setup only. Production callers explicitly select their original
    // containing element; traversal performs no global point search.
    (0..snapshot.mesh().element_count())
        .map(|i| snapshot.element(i).unwrap())
        .find(|e| e.barycentric(p).is_ok())
        .unwrap()
}
fn trace(
    cx: &Cx<'_>,
    mesh: &P1ConductionMesh<'_>,
    values: &[f64],
    ray: Ray,
    length: f64,
    opts: P1MeshOptions,
) -> Result<P1MeshPropagation, P1MeshError> {
    let snapshot =
        P1ConductionTemperatureSnapshot::try_new(cx, mesh.snapshot.mesh(), values).unwrap();
    let element = locate(snapshot, ray.origin);
    let glass = reference();
    let law = law();
    let binding =
        MaterialThermoOpticField::try_new(&glass, &law, &element, ContentHash([95; 32])).unwrap();
    mesh.propagate(
        cx,
        &binding,
        ray,
        Vec3::new(1.0, 0.0, 0.0),
        1500.0,
        length,
        opts,
    )
}

// Independent exact affine-index ray: d(n u)/ds = g. Resolve the conserved
// momentum into components parallel/perpendicular to g and integrate in s.
fn affine(ray: Ray, n0: f64, g: Vec3, s: f64) -> (Ray, f64) {
    let b = g.dot(g).sqrt();
    let axis = Vec3::new(g.x / b, g.y / b, g.z / b);
    let parallel = n0 * ray.dir.dot(axis);
    let perpendicular = Vec3::new(
        n0 * ray.dir.x - parallel * axis.x,
        n0 * ray.dir.y - parallel * axis.y,
        n0 * ray.dir.z - parallel * axis.z,
    );
    let a = perpendicular.dot(perpendicular).sqrt();
    let q = parallel + b * s;
    let n = a.hypot(q);
    let integral = ((q / a).asinh() - (parallel / a).asinh()) / b;
    let parallel_distance = s * (2.0 * parallel + b * s) / (n + n0);
    let point = ray.origin.offset(Vec3::new(
        perpendicular.x * integral + axis.x * parallel_distance,
        perpendicular.y * integral + axis.y * parallel_distance,
        perpendicular.z * integral + axis.z * parallel_distance,
    ));
    let direction = Vec3::new(
        (perpendicular.x + q * axis.x) / n,
        (perpendicular.y + q * axis.y) / n,
        (perpendicular.z + q * axis.z) / n,
    );
    let primitive = |q: f64| 0.5 * (q * a.hypot(q) + a * a * (q / a).asinh());
    (
        Ray {
            origin: point,
            dir: direction,
        },
        (primitive(q) - primitive(parallel)) / b,
    )
}

#[test]
fn original_face_ids_carry_uniform_and_affine_rays_through_twelve_tetrahedra() {
    with_cx(|_, cx| {
        for reverse in [false, true] {
            let mesh = column(reverse);
            for gradient in [[0.0, 0.0], [12.0, -7.0]] {
                let values: Vec<_> = mesh
                    .positions()
                    .iter()
                    .map(|p| 300.0 + gradient[0] * p[0] + gradient[1] * p[1])
                    .collect();
                let snapshot =
                    P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
                let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
                let ray = launch();
                let length = 1.75;
                let result = trace(cx, &prepared, &values, ray, length, options()).unwrap();
                assert!(result.internal_faces >= 5);
                assert_eq!(result.termination, BoundaryTermination::DistanceLimit);
                assert_eq!(result.propagation.distance_m, length);
                let (expected, phase) = if gradient == [0.0, 0.0] {
                    (
                        Ray {
                            origin: ray.origin.offset(Vec3::new(0.0, 0.0, length)),
                            dir: ray.dir,
                        },
                        1.5 * length,
                    )
                } else {
                    let n0 =
                        1.5 + 0.00125 * (gradient[0] * ray.origin.x + gradient[1] * ray.origin.y);
                    affine(
                        ray,
                        n0,
                        Vec3::new(0.00125 * gradient[0], 0.00125 * gradient[1], 0.0),
                        length,
                    )
                };
                point_near(result.propagation.ray.origin, expected.origin, 5.0e-10);
                vector_near(result.propagation.ray.dir, expected.dir, 5.0e-10);
                near(result.propagation.optical_path_length_m, phase, 7.0e-10);
                near(
                    result.propagation.group_optical_path_length_m.unwrap(),
                    phase,
                    7.0e-10,
                );
                near(result.propagation.optical_depth, 0.2 * length, 3.0e-12);
                assert!(
                    snapshot
                        .element(result.element_index)
                        .unwrap()
                        .barycentric(result.propagation.ray.origin)
                        .is_ok()
                );
                assert_eq!(
                    result,
                    trace(cx, &prepared, &values, ray, length, options()).unwrap()
                );
            }
        }
    });
}

#[test]
fn continuous_temperature_gradient_jump_matches_two_independent_affine_solutions() {
    with_cx(|_, cx| {
        let mesh = column(false);
        let values: Vec<_> = mesh
            .positions()
            .iter()
            .map(|p| {
                300.0 + 8.0 * p[0] - 5.0 * p[1] + (if p[2] < 0.0 { 12.0 } else { -7.0 }) * p[2]
            })
            .collect();
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let ray = launch();
        let length = 1.7;
        let result = trace(cx, &prepared, &values, ray, length, options()).unwrap();
        let n0 = 1.5 + 0.00125 * (8.0 * ray.origin.x - 5.0 * ray.origin.y + 12.0 * ray.origin.z);
        let g0 = Vec3::new(0.01, -0.00625, 0.015);
        let (mut low, mut high) = (0.0, 1.5);
        for _ in 0..60 {
            let mid = 0.5 * (low + high);
            if affine(ray, n0, g0, mid).0.origin.z < 0.0 {
                low = mid
            } else {
                high = mid
            }
        }
        let crossing = 0.5 * (low + high);
        let (mut first, phase0) = affine(ray, n0, g0, crossing);
        first.origin.z = 0.0;
        let n1 = 1.5 + 0.00125 * (8.0 * first.origin.x - 5.0 * first.origin.y);
        let (expected, phase1) = affine(
            first,
            n1,
            Vec3::new(0.01, -0.00625, -0.00875),
            length - crossing,
        );
        point_near(result.propagation.ray.origin, expected.origin, 3.0e-9);
        vector_near(result.propagation.ray.dir, expected.dir, 3.0e-9);
        near(
            result.propagation.optical_path_length_m,
            phase0 + phase1,
            4.0e-9,
        );
        near(
            result.propagation.group_optical_path_length_m.unwrap(),
            phase0 + phase1,
            4.0e-9,
        );
        near(
            result.propagation.power_transmittance(),
            (-0.2 * length).exp(),
            3.0e-12,
        );
        assert!(result.internal_faces >= 5);
    });
}

#[test]
fn exact_distance_limits_inside_both_sides_of_a_transition_collar_lose_no_path() {
    with_cx(|_, cx| {
        let mesh = bipyramid(1.0, false);
        let values = [300.0; 5];
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let start = snapshot.element(0).unwrap();
        let glass = reference();
        let law = law();
        let binding =
            MaterialThermoOpticField::try_new(&glass, &law, &start, ContentHash([96; 32])).unwrap();
        let field = binding.prepare(1500.0).unwrap();
        let ray = Ray {
            origin: Point3::new(0.1, 0.0, -0.5),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let mut opts = options();
        opts.event.surface_tolerance_m = 1.0e-4;
        let mut oracle = start
            .boundary_oracle(43, opts.event.surface_tolerance_m)
            .unwrap();
        let first = propagate_to_boundary(
            &field,
            &mut oracle,
            ray,
            Vec3::new(1.0, 0.0, 0.0),
            1500.0,
            1.0,
            None,
            opts.event,
        )
        .unwrap();
        let gap = 0.5 - first.propagation.distance_m;
        assert!(gap > 0.0 && gap < opts.event.surface_tolerance_m);
        for (length, element) in [(0.5 - 0.25 * gap, 0), (0.5 + 0.25 * gap, 1)] {
            let result = prepared
                .propagate(
                    cx,
                    &binding,
                    ray,
                    Vec3::new(1.0, 0.0, 0.0),
                    1500.0,
                    length,
                    opts,
                )
                .unwrap();
            assert_eq!(result.termination, BoundaryTermination::DistanceLimit);
            assert_eq!(result.element_index, element);
            assert_eq!(result.internal_faces, 1);
            assert_eq!(result.propagation.distance_m, length);
            point_near(
                result.propagation.ray.origin,
                ray.origin.offset(Vec3::new(0.0, 0.0, length)),
                3.0e-14,
            );
            near(
                result.propagation.optical_path_length_m,
                1.5 * length,
                3.0e-14,
            );
            near(
                result.propagation.group_optical_path_length_m.unwrap(),
                1.5 * length,
                3.0e-14,
            );
            near(result.propagation.optical_depth, 0.2 * length, 3.0e-14);
        }
    });
}

#[test]
fn external_and_duplicated_contact_faces_stop_and_narrow_neighbours_refuse() {
    with_cx(|_, cx| {
        for split in [false, true] {
            let mesh = bipyramid(1.0, split);
            let values = vec![300.0; mesh.vertex_count()];
            let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
            let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
            let ray = Ray {
                origin: Point3::new(0.1, 0.0, -0.5),
                dir: Vec3::new(0.0, 0.0, 1.0),
            };
            let result = trace(cx, &prepared, &values, ray, 2.0, options()).unwrap();
            let BoundaryTermination::Surface(witness) = result.termination else {
                panic!("mesh boundary was crossed")
            };
            assert_eq!(witness.boundary_id, 43);
            assert!(
                mesh.boundary()
                    .iter()
                    .any(|face| face.face == witness.feature_id)
            );
            if split {
                assert_eq!(result.element_index, 0);
                assert_eq!(result.internal_faces, 0);
                near(result.propagation.distance_m, 0.5, 2.0e-8);
            } else {
                assert_eq!(result.element_index, 1);
                assert_eq!(result.internal_faces, 1);
                near(result.propagation.distance_m, 1.3, 4.0e-8);
            }
        }
        let mesh = bipyramid(1.0e-7, false);
        let values = [300.0; 5];
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let ray = Ray {
            origin: Point3::new(0.1, 0.0, -0.5),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let mut opts = options();
        opts.event.surface_tolerance_m = 1.0e-6;
        assert_eq!(
            trace(cx, &prepared, &values, ray, 0.6, opts).unwrap_err(),
            P1MeshError::Propagation(GradientIndexError::BoundaryResolutionLimit)
        );
    });
}

#[test]
fn original_transient_cooling_drives_the_complete_two_element_optical_path() {
    with_cx(|_, cx| {
        let mesh = bipyramid(1.0, false);
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
        // For the x-odd mode [-1,1,0,0,0], the two original tetrahedra
        // have total K eigenvalue 2/3 and lumped C=capacity/3=2/3.
        // Backward Euler at dt=1 halves the amplitude, conserving total heat.
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
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let length = 0.9;
        let mut bends = [0.0; 2];
        for (i, (values, amplitude)) in [(&old[..], 20.0), (&accepted.temperature[..], 10.0)]
            .into_iter()
            .enumerate()
        {
            let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, values).unwrap();
            let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
            let result = trace(cx, &prepared, values, ray, length, options()).unwrap();
            let (expected, phase) =
                affine(ray, 1.5, Vec3::new(0.00125 * amplitude, 0.0, 0.0), length);
            point_near(result.propagation.ray.origin, expected.origin, 3.0e-10);
            vector_near(result.propagation.ray.dir, expected.dir, 3.0e-10);
            near(result.propagation.optical_path_length_m, phase, 3.0e-10);
            near(
                result.propagation.group_optical_path_length_m.unwrap(),
                phase,
                3.0e-10,
            );
            near(
                result.propagation.power_transmittance(),
                (-0.2 * length).exp(),
                3.0e-12,
            );
            assert_eq!(result.internal_faces, 1);
            assert_eq!(result.element_index, 1);
            bends[i] = result.propagation.ray.dir.x;
        }
        assert!(bends[0] > 0.01 && bends[1] > 0.49 * bends[0] && bends[1] < 0.51 * bends[0]);
    });
}

#[test]
fn identity_orientation_cancellation_and_whole_path_work_budgets_are_retained() {
    with_cx(|gate, cx| {
        let mesh = column(false);
        let values = vec![300.0; mesh.vertex_count()];
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let other = values.clone();
        assert_eq!(
            trace(cx, &prepared, &other, launch(), 1.75, options()).unwrap_err(),
            P1MeshError::ForeignSnapshot
        );
        let result = trace(cx, &prepared, &values, launch(), 1.75, options()).unwrap();
        let mut opts = options();
        opts.max_internal_faces = 0;
        assert_eq!(
            trace(cx, &prepared, &values, launch(), 1.75, opts).unwrap_err(),
            P1MeshError::FaceBudgetExceeded
        );
        opts = options();
        opts.event.integration.max_steps =
            result.propagation.accepted_steps + result.propagation.rejected_steps - 1;
        assert_eq!(
            trace(cx, &prepared, &values, launch(), 1.75, opts).unwrap_err(),
            P1MeshError::Propagation(GradientIndexError::StepBudgetExceeded)
        );
        opts = options();
        opts.event.max_boundary_queries = result.boundary_evaluations - 1;
        assert_eq!(
            trace(cx, &prepared, &values, launch(), 1.75, opts).unwrap_err(),
            P1MeshError::Propagation(GradientIndexError::BoundaryQueryBudgetExceeded)
        );

        let overlap = bipyramid(-0.5, false);
        let temperatures = [300.0; 5];
        let invalid =
            P1ConductionTemperatureSnapshot::try_new(cx, &overlap, &temperatures).unwrap();
        assert!(matches!(
            P1ConductionMesh::try_new(cx, invalid),
            Err(P1MeshError::InvalidTopology { .. })
        ));
        gate.request();
        assert!(matches!(
            P1ConductionMesh::try_new(cx, snapshot),
            Err(P1MeshError::Propagation(GradientIndexError::Cancelled))
        ));
        let element = snapshot.element(0).unwrap();
        let glass = reference();
        let law = law();
        let binding =
            MaterialThermoOpticField::try_new(&glass, &law, &element, ContentHash([95; 32]))
                .unwrap();
        assert_eq!(
            prepared
                .propagate(
                    cx,
                    &binding,
                    launch(),
                    Vec3::new(1.0, 0.0, 0.0),
                    1500.0,
                    1.75,
                    options()
                )
                .unwrap_err(),
            P1MeshError::Propagation(GradientIndexError::Cancelled)
        );
    });
}
