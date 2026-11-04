use super::*;
use crate::dielectric::BeerLambertAbsorption;
use crate::gradient_index::{GradientIndexOptions, GradientIndexPropagation};
use crate::gradient_index::bundle::{BundleError, BundleOptions, RayBundlePropagation, propagate_bundle};
use crate::gradient_index::event::{BoundaryTermination, GradientEventOptions, propagate_to_boundary};
use crate::material_state::thermo_optic::{TemperatureGradientSample, TemperatureJet};
use crate::material_state::thermo_optic::transient::mesh::transition::{MeshBoundary, Transition};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 121, kernel_id: 1, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&gate, &cx);
    });
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}; tolerance {tolerance}");
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance); near(a.y, b.y, tolerance); near(a.z, b.z, tolerance);
}
fn point(x: [f64; 3]) -> Point3 { Point3::new(x[0], x[1], x[2]) }
fn unit(v: Vec3) -> Vec3 { v.scale(1.0 / v.dot(v).sqrt()) }
fn same_path(a: &GradientIndexPropagation, b: &GradientIndexPropagation, tolerance: f64) {
    vector_near(a.ray.origin.delta_from(b.ray.origin), Vec3::new(0.0, 0.0, 0.0), tolerance);
    vector_near(a.ray.dir, b.ray.dir, tolerance);
    near(a.distance_m, b.distance_m, 5e-13);
    near(a.optical_path_length_m, b.optical_path_length_m, tolerance);
    near(a.group_optical_path_length_m.unwrap(), b.group_optical_path_length_m.unwrap(), tolerance);
    near(a.optical_depth, 0.02 * a.distance_m, 5e-13);
    near(a.power_transmittance(), (-0.02 * a.distance_m).exp(), 5e-13);
}
fn options() -> P1MeshBundleOptions {
    P1MeshBundleOptions {
        mesh: P1MeshOptions {
            boundary_id: 121, max_internal_faces: 128,
            event: GradientEventOptions {
                integration: GradientIndexOptions {
                    initial_step_m: 0.01, max_step_m: 0.04, min_step_m: 1e-13,
                    relative_tolerance: 2e-12, position_tolerance_m: 2e-12,
                    direction_tolerance: 2e-12, optical_path_tolerance_m: 2e-12,
                    optical_depth_tolerance: 2e-12, ..GradientIndexOptions::default()
                },
                surface_tolerance_m: 1e-8, ..GradientEventOptions::default()
            },
        },
        position_derivative_tolerance_m: 2e-12,
        direction_derivative_tolerance: 2e-12,
    }
}
fn column(reverse: bool) -> ConductionMesh {
    let mut positions = Vec::new();
    for z in [-1.0, 0.0, 1.0] {
        positions.extend([[-1.0, -1.0, z], [1.0, -1.0, z], [-1.0, 1.0, z], [1.0, 1.0, z]]);
    }
    let mut tets = Vec::new();
    for base in [0, 4] {
        for a in [[1, 2, 4], [1, 4, 2], [2, 1, 4], [2, 4, 1], [4, 1, 2], [4, 2, 1]] {
            let mut tet = [base, base + a[0], base + a[0] + a[1], base + 7];
            if reverse { tet.swap(0, 1); }
            tets.push(tet);
        }
    }
    ConductionMesh::new(TetComplex::from_tets(12, tets), positions).unwrap()
}
fn bipyramid(contact: bool) -> ConductionMesh {
    let mut positions = vec![
        [-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0], [0.0, 0.0, 1.0],
    ];
    let tets = if contact {
        positions.extend([positions[0], positions[1], positions[2]]);
        vec![[0, 1, 2, 3], [5, 6, 7, 4]]
    } else { vec![[0, 1, 2, 3], [0, 1, 2, 4]] };
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}
fn source<'a>(
    cx: &Cx<'_>, mesh: &'a ConductionMesh, temperature: &'a [f64], pressure: &'a [f64],
) -> P1PressureTemperatureSnapshot<'a> {
    let thermal = P1ConductionTemperatureSnapshot::try_new(cx, mesh, temperature).unwrap();
    let pressure = P1PressureSnapshot::try_new(
        cx, mesh, NodalPressureValues::Pascals(pressure), ContentHash([121; 32]),
    ).unwrap();
    P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([122; 32]), pressure).unwrap()
}
fn locate(source: P1PressureTemperatureSnapshot<'_>, p: Point3) -> P1PressureTemperatureElement<'_> {
    // Fixture setup only; production traversal requires the original launch cell.
    (0..source.temperature_snapshot().mesh().element_count())
        .map(|i| source.element(i).unwrap())
        .find(|e| e.temperature_element().barycentric(p).is_ok()).unwrap()
}
fn trace<'a>(
    cx: &Cx<'_>, mesh: &P1ConductionMesh<'a>, source: (P1PressureTemperatureSnapshot<'a>, If97Region),
    ray: Ray, wavelength: f64, length: f64, options: P1MeshBundleOptions,
) -> Result<P1MeshBundlePropagation, PressureTemperatureMeshError> {
    let element = locate(source.0, ray.origin);
    let state = element.bind_water(source.1).unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    // The outer optical identity is deliberately distinct from the paired state.
    let binding = MaterialWaterRefractiveField::try_new(&state, &absorption, ContentHash([123; 32])).unwrap();
    mesh.propagate_pressure_temperature_water_bundle(
        cx, &binding, ray, Vec3::new(1.0, 0.0, 0.0), wavelength, length, options,
    )
}
fn primal<'a>(
    cx: &Cx<'_>, mesh: &P1ConductionMesh<'a>, source: (P1PressureTemperatureSnapshot<'a>, If97Region),
    ray: Ray, wavelength: f64, length: f64, options: P1MeshOptions,
) -> Result<P1MeshPropagation, PressureTemperatureMeshError> {
    let element = locate(source.0, ray.origin);
    let state = element.bind_water(source.1).unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    let binding = MaterialWaterRefractiveField::try_new(&state, &absorption, ContentHash([123; 32])).unwrap();
    mesh.propagate_pressure_temperature_water(
        cx, &binding, ray, Vec3::new(1.0, 0.0, 0.0), wavelength, length, options,
    )
}
struct AffineState { temperature: f64, pressure: f64, gt: Vec3, gp: Vec3 }
impl PressureTemperatureGradientField for AffineState {
    fn pressure_temperature_gradient(&self, p: Point3) -> Result<PressureTemperatureGradientSample, GradientIndexError> {
        let x = Vec3::new(p.x, p.y, p.z);
        Ok(PressureTemperatureGradientSample {
            temperature: TemperatureGradientSample {
                temperature_k: self.temperature + self.gt.dot(x), gradient_k_per_m: self.gt,
            },
            pressure_pa: self.pressure + self.gp.dot(x), pressure_gradient_pa_per_m: self.gp,
        })
    }
}
impl PressureTemperatureSecondOrderField for AffineState {
    fn pressure_temperature_jet(&self, p: Point3) -> Result<PressureTemperatureJet, GradientIndexError> {
        let sample = self.pressure_temperature_gradient(p)?;
        Ok(PressureTemperatureJet {
            temperature: TemperatureJet { sample: sample.temperature, diagonal_k_per_m2: [0.0; 3], off_diagonal_k_per_m2: [0.0; 3] },
            pressure_pa: sample.pressure_pa, pressure_gradient_pa_per_m: self.gp,
            pressure_diagonal_pa_per_m2: [0.0; 3], pressure_off_diagonal_pa_per_m2: [0.0; 3],
        })
    }
}
fn smooth(source: &AffineState, region: If97Region, ray: Ray, wavelength: f64, length: f64) -> RayBundlePropagation {
    let state = PressureTemperatureWaterField::try_new(source, region, ContentHash([124; 32])).unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    let binding = state.bind_refraction(&absorption).unwrap();
    let prepared = binding.prepare(wavelength).unwrap();
    propagate_bundle(&prepared, ray, Vec3::new(1.0, 0.0, 0.0), wavelength, length, BundleOptions {
        integration: options().mesh.event.integration,
        position_derivative_tolerance_m: options().position_derivative_tolerance_m,
        direction_derivative_tolerance: options().direction_derivative_tolerance,
    }).unwrap()
}

#[test]
fn g1_joint_nonparallel_fields_cross_original_cells_in_both_phases_and_orientations() {
    with_cx(|_, cx| {
        for reverse in [false, true] {
            let mesh = column(reverse);
            for (t0, p0, gp, region) in [
                (300.0, 30e6, Vec3::new(8e6, -6e6, 3e6), If97Region::Region1),
                (400.0, 8e4, Vec3::new(1e4, -0.8e4, 0.5e4), If97Region::Region2),
            ] {
                let original = AffineState { temperature: t0, pressure: p0, gt: Vec3::new(5.0, -3.0, 2.0), gp };
                let temperature: Vec<_> = mesh.positions().iter().map(|p| t0 + original.gt.dot(Vec3::new(p[0], p[1], p[2]))).collect();
                let pressure: Vec<_> = mesh.positions().iter().map(|p| p0 + gp.dot(Vec3::new(p[0], p[1], p[2]))).collect();
                let joint = source(cx, &mesh, &temperature, &pressure);
                let prepared = P1ConductionMesh::try_new(cx, joint.temperature_snapshot()).unwrap();
                let ray = Ray { origin: point([0.17, -0.23, -0.9]), dir: Vec3::new(0.0, 0.0, 1.0) };
                for wavelength in [589.0, 900.0] {
                    let result = trace(cx, &prepared, (joint, region), ray, wavelength, 1.7, options()).unwrap();
                    let central = primal(cx, &prepared, (joint, region), ray, wavelength, 1.7, options().mesh).unwrap();
                    let expected = smooth(&original, region, ray, wavelength, 1.7);
                    same_path(result.bundle.propagation(), expected.propagation(), 5e-9);
                    same_path(&central.propagation, expected.propagation(), 5e-9);
                    assert!(result.internal_faces >= 5 && central.internal_faces >= 5);
                    assert_eq!(result.bundle.interface_events(), 0);
                    for (a, b) in result.bundle.angular_derivatives().into_iter().zip(expected.angular_derivatives()) {
                        vector_near(a.position_m_per_rad, b.position_m_per_rad, 3e-8);
                        vector_near(a.direction_per_rad, b.direction_per_rad, 3e-8);
                    }
                    let p = result.bundle.propagation();
                    assert!((p.group_optical_path_length_m.unwrap() - p.optical_path_length_m).abs() > 1e-8);
                }
            }
        }
    });
}

#[test]
fn g1_pressure_only_gradient_jump_retains_native_bundle_variations() {
    with_cx(|_, cx| {
        let mesh = bipyramid(false);
        let temperature = [300.0; 5];
        let pressure: Vec<_> = mesh.positions().iter().map(|p| {
            30e6 + 4e6 * p[0] - 1e6 * p[1] + if p[2] < 0.0 { 8e6 * p[2] } else { -6e6 * p[2] }
        }).collect();
        let joint = source(cx, &mesh, &temperature, &pressure);
        let prepared = P1ConductionMesh::try_new(cx, joint.temperature_snapshot()).unwrap();
        let ray = Ray { origin: point([-0.06, -0.04, -0.6]), dir: unit(Vec3::new(0.18, 0.07, 1.0)) };
        let result = trace(cx, &prepared, (joint, If97Region::Region1), ray, 589.0, 1.0, options()).unwrap();
        let central = primal(cx, &prepared, (joint, If97Region::Region1), ray, 589.0, 1.0, options().mesh).unwrap();
        same_path(result.bundle.propagation(), &central.propagation, 5e-10);
        assert_eq!(result.internal_faces, 1);
        assert_eq!(result.element_index, 1);
        assert_eq!(result.bundle.interface_events(), 0);
        let epsilon = 2e-5;
        for (basis, actual) in result.bundle.launch_basis().into_iter().zip(result.bundle.angular_derivatives()) {
            let neighbour = |sign: f64| {
                let d = Vec3::new(ray.dir.x + sign * epsilon * basis.x,
                    ray.dir.y + sign * epsilon * basis.y, ray.dir.z + sign * epsilon * basis.z);
                primal(cx, &prepared, (joint, If97Region::Region1), Ray { dir: unit(d), ..ray },
                    589.0, 1.0, options().mesh).unwrap().propagation.ray
            };
            let (a, b) = (neighbour(1.0), neighbour(-1.0));
            vector_near(actual.position_m_per_rad, a.origin.delta_from(b.origin).scale(0.5 / epsilon), 5e-7);
            vector_near(actual.direction_per_rad,
                Vec3::new(a.dir.x - b.dir.x, a.dir.y - b.dir.y, a.dir.z - b.dir.z).scale(0.5 / epsilon), 5e-7);
        }
    });
}

#[test]
fn g1_departure_collar_queries_keep_original_pressure_temperature_pairs() {
    with_cx(|_, cx| {
        let mesh = bipyramid(false);
        let temperature: Vec<_> = mesh.positions().iter().map(|p| {
            355.0 + 5.0 * p[0] - 3.0 * p[1] + if p[2] < 0.0 { 15.0 * p[2] } else { -20.0 * p[2] }
        }).collect();
        let pressure: Vec<_> = mesh.positions().iter().map(|p| {
            30e6 + 4e6 * p[0] - 1e6 * p[1] + if p[2] < 0.0 { 8e6 * p[2] } else { -6e6 * p[2] }
        }).collect();
        let snapshot = source(cx, &mesh, &temperature, &pressure);
        let (previous, next) = (snapshot.element(0).unwrap(), snapshot.element(1).unwrap());
        let cursor = MeshTemperature { current: Cell::new(previous.temperature_element()), previous: Cell::new(None) };
        let field = MeshPressureTemperature::new(&cursor, previous);
        let lower = point([0.02, 0.01, -1e-4]);
        assert_eq!(field.sample_jet(lower).unwrap(), previous.sample_jet(lower).unwrap());
        cursor.previous.set(Some(previous.temperature_element()));
        cursor.current.set(next.temperature_element());
        assert!(field.pressure_gradient_jump(previous.temperature_element(), next.temperature_element()).unwrap());
        for z in [-1e-4, 1e-4, -2e-5, 2e-5, -1e-6, 1e-6] {
            let p = point([0.02, 0.01, z]);
            let original = if z < 0.0 { previous } else { next };
            // Compare both values AND both one-sided gradients from the same
            // original element. A mixed-side state cannot satisfy this equality.
            assert_eq!(field.sample_jet(p).unwrap(), original.sample_jet(p).unwrap());
        }
        assert!(matches!(field.sample_jet(point([2.0, 0.0, 0.0])), Err(GradientIndexError::OutsideFieldDomain)));
    });
}

#[test]
fn g3_pressure_only_jump_refuses_an_unresolved_face_bundle_endpoint() {
    with_cx(|_, cx| {
        for jump in [false, true] {
            let mesh = bipyramid(false);
            let temperature = [300.0; 5];
            let pressure = [30e6, 30e6, 30e6, 22e6, if jump { 22e6 } else { 38e6 }];
            let snapshot = source(cx, &mesh, &temperature, &pressure);
            let mesh_owner = P1ConductionMesh::try_new(cx, snapshot.temperature_snapshot()).unwrap();
            let element = snapshot.element(0).unwrap();
            let thermal = element.temperature_element();
            let state = element.bind_water(If97Region::Region1).unwrap();
            let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
            let binding = state.bind_refraction(&absorption).unwrap();
            let field = binding.prepare(589.0).unwrap();
            let ray = Ray { origin: point([0.0, 0.0, -0.5]), dir: Vec3::new(0.0, 0.0, 1.0) };
            let mut opts = options();
            opts.mesh.event.surface_tolerance_m = 1e-4;
            let mut boundary = thermal.boundary_oracle(opts.mesh.boundary_id, opts.mesh.event.surface_tolerance_m).unwrap();
            let arrival = propagate_to_boundary(&field, &mut boundary, ray, Vec3::new(1.0, 0.0, 0.0),
                589.0, 1.0, None, opts.mesh.event).unwrap();
            let BoundaryTermination::Surface(witness) = arrival.termination else { panic!("missing original face"); };
            let neighbour = mesh_owner.neighbours[0][witness.feature_id].unwrap();
            let next = snapshot.element(neighbour.element).unwrap().temperature_element();
            assert_eq!(thermal.gradient, next.gradient); // temperature alone sees no jump
            let cursor = MeshTemperature { current: Cell::new(thermal), previous: Cell::new(None) };
            let joint = MeshPressureTemperature::new(&cursor, element);
            let pressure_jump = joint.pressure_gradient_jump(thermal, next).unwrap();
            assert_eq!(pressure_jump, jump);
            let guarded = MeshBoundary::with_transition(
                next.boundary_oracle(opts.mesh.boundary_id, opts.mesh.event.surface_tolerance_m).unwrap(),
                Transition::new(thermal, next, neighbour.local_face, arrival.propagation.ray.origin, witness)
                    .with_additional_gradient_jump(pressure_jump),
            );
            // Exercise the original exact-plane guard directly. Do not snap a
            // numerically arrived ray onto a face to manufacture an endpoint.
            let admission = guarded.admit_bundle_endpoint(point([0.0, 0.0, 0.0]));
            if jump {
                assert_eq!(admission, Err(GradientIndexError::BoundaryResolutionLimit));
            } else {
                assert_eq!(admission, Ok(()));
            }
        }
    });
}

#[test]
fn g3_joint_mesh_preserves_source_identity_contacts_phase_and_work_limits() {
    with_cx(|gate, cx| {
        let mesh = bipyramid(false);
        let temperature = [300.0; 5];
        let pressure = [30e6; 5];
        let snapshot = source(cx, &mesh, &temperature, &pressure);
        let prepared = P1ConductionMesh::try_new(cx, snapshot.temperature_snapshot()).unwrap();
        let ray = Ray { origin: point([0.0, 0.0, -0.5]), dir: Vec3::new(0.0, 0.0, 1.0) };
        let element = snapshot.element(0).unwrap();
        let wrong = PressureTemperatureWaterField::try_new(&element, If97Region::Region1, ContentHash([125; 32])).unwrap();
        let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
        let binding = MaterialWaterRefractiveField::try_new(&wrong, &absorption, ContentHash([126; 32])).unwrap();
        assert!(matches!(prepared.propagate_pressure_temperature_water(
            cx, &binding, ray, Vec3::new(1.0, 0.0, 0.0), 589.0, 0.9, options().mesh,
        ), Err(PressureTemperatureMeshError::SourceIdentityMismatch)));
        assert!(matches!(prepared.propagate_pressure_temperature_water_bundle(
            cx, &binding, ray, Vec3::new(1.0, 0.0, 0.0), 589.0, 0.9, options(),
        ), Err(PressureTemperatureMeshError::SourceIdentityMismatch)));
        let copied = temperature.to_vec();
        let foreign = source(cx, &mesh, &copied, &pressure);
        assert!(matches!(trace(cx, &prepared, (foreign, If97Region::Region1), ray, 589.0, 0.9, options()),
            Err(PressureTemperatureMeshError::Mesh(P1MeshError::ForeignSnapshot))));
        let mut limited = options();
        limited.mesh.max_internal_faces = 0;
        assert!(matches!(trace(cx, &prepared, (snapshot, If97Region::Region1), ray, 589.0, 0.9, limited),
            Err(PressureTemperatureMeshError::Mesh(P1MeshError::FaceBudgetExceeded))));

        let contact = bipyramid(true);
        let contact_temperature = vec![300.0; contact.vertex_count()];
        let mut contact_pressure = vec![30e6; contact.vertex_count()];
        for vertex in [4, 5, 6, 7] { contact_pressure[vertex] = 40e6; }
        let contact_source = source(cx, &contact, &contact_temperature, &contact_pressure);
        let contact_mesh = P1ConductionMesh::try_new(cx, contact_source.temperature_snapshot()).unwrap();
        let stop = trace(cx, &contact_mesh, (contact_source, If97Region::Region1), ray, 589.0, 0.9, options()).unwrap();
        assert!(matches!(stop.termination, BoundaryTermination::Surface(_)));
        assert_eq!(stop.internal_faces, 0);
        assert_eq!(stop.bundle.interface_events(), 0);
        near(stop.bundle.propagation().distance_m, 0.5, 1e-8);

        let boiling_temperature = [373.15; 5];
        let boiling_pressure = [2e5, 2e5, 2e5, 2.5e5, 0.5e5];
        let boiling = source(cx, &mesh, &boiling_temperature, &boiling_pressure);
        let boiling_mesh = P1ConductionMesh::try_new(cx, boiling.temperature_snapshot()).unwrap();
        assert!(matches!(primal(cx, &boiling_mesh, (boiling, If97Region::Region1), ray, 589.0, 1.2, options().mesh),
            Err(PressureTemperatureMeshError::Mesh(P1MeshError::Propagation(GradientIndexError::OutsideFieldDomain)))));
        gate.request();
        assert!(matches!(trace(cx, &prepared, (snapshot, If97Region::Region1), ray, 589.0, 0.9, options()),
            Err(PressureTemperatureMeshError::Mesh(P1MeshError::Propagation(GradientIndexError::Cancelled)))
            | Err(PressureTemperatureMeshError::Mesh(P1MeshError::Bundle(BundleError::Propagation(GradientIndexError::Cancelled))))));
    });
}
