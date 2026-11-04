use super::*;
use crate::dielectric::BeerLambertAbsorption;
use crate::gradient_index::{GradientIndexOptions, GradientIndexPropagation};
use crate::gradient_index::bundle::{
    BundleError, BundleOptions, ProjectionOptions, RayBundlePropagation, propagate_bundle,
};
use crate::material_state::thermo_optic::{
    TemperatureGradientField, TemperatureSecondOrderField,
};
use fs_blake3::ContentHash;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig, ScalarField,
    ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_qty::Pressure;
use fs_rep_mesh::TetComplex;
use fs_thermochem::water::if97::If97Region;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            StreamKey { seed: 103, kernel_id: 1, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        );
        f(&gate, &cx);
    });
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}; tolerance {tolerance}");
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.dot(v).sqrt())
}
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
            boundary_id: 103,
            max_internal_faces: 128,
            event: GradientEventOptions {
                integration: GradientIndexOptions {
                    initial_step_m: 0.01,
                    max_step_m: 0.04,
                    min_step_m: 1e-13,
                    relative_tolerance: 2e-12,
                    position_tolerance_m: 2e-12,
                    direction_tolerance: 2e-12,
                    optical_path_tolerance_m: 2e-12,
                    optical_depth_tolerance: 2e-12,
                    ..GradientIndexOptions::default()
                },
                surface_tolerance_m: 1e-8,
                ..GradientEventOptions::default()
            },
        },
        position_derivative_tolerance_m: 2e-12,
        direction_derivative_tolerance: 2e-12,
    }
}
fn column(reverse: bool) -> ConductionMesh {
    let mut positions = Vec::new();
    for z in [-1.0, 0.0, 1.0] {
        positions.extend([
            [-1.0, -1.0, z], [1.0, -1.0, z], [-1.0, 1.0, z], [1.0, 1.0, z],
        ]);
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
    } else {
        vec![[0, 1, 2, 3], [0, 1, 2, 4]]
    };
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}
fn locate(snapshot: P1ConductionTemperatureSnapshot<'_>, p: Point3) -> P1ConductionElement<'_> {
    // Fixture setup only. Production callers select their original launch cell.
    (0..snapshot.mesh().element_count())
        .map(|i| snapshot.element(i).unwrap())
        .find(|e| e.barycentric(p).is_ok()).unwrap()
}
fn trace(
    cx: &Cx<'_>, mesh: &P1ConductionMesh<'_>, ray: Ray, length: f64,
    water: (f64, If97Region), wavelength: f64, options: P1MeshBundleOptions,
) -> Result<P1MeshBundlePropagation, WaterMeshError> {
    let element = locate(mesh.snapshot(), ray.origin);
    let state = IsobaricWaterField::try_new(
        &element, Pressure::new(water.0), water.1, ContentHash([104; 32]),
    ).unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    // Deliberately distinct inner-state and outer-optical snapshot bindings.
    let binding = MaterialWaterRefractiveField::try_new(&state, &absorption, ContentHash([105; 32])).unwrap();
    mesh.propagate_water_bundle(
        cx, &binding, ray, Vec3::new(1.0, 0.0, 0.0), wavelength, length, options,
    )
}
fn primal(
    cx: &Cx<'_>, mesh: &P1ConductionMesh<'_>, ray: Ray, length: f64,
    water: (f64, If97Region), wavelength: f64, options: P1MeshOptions,
) -> Result<P1MeshPropagation, WaterMeshError> {
    let element = locate(mesh.snapshot(), ray.origin);
    let state = IsobaricWaterField::try_new(
        &element, Pressure::new(water.0), water.1, ContentHash([104; 32]),
    ).unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    let binding = MaterialWaterRefractiveField::try_new(&state, &absorption, ContentHash([105; 32])).unwrap();
    mesh.propagate_water(
        cx, &binding, ray, Vec3::new(1.0, 0.0, 0.0), wavelength, length, options,
    )
}
struct AffineTemperature { reference: f64, gradient: Vec3 }
impl TemperatureField for AffineTemperature {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        Ok(self.reference + self.gradient.dot(Vec3::new(p.x, p.y, p.z)))
    }
}
impl TemperatureGradientField for AffineTemperature {
    fn temperature_gradient(&self, p: Point3) -> Result<TemperatureGradientSample, GradientIndexError> {
        Ok(TemperatureGradientSample {
            temperature_k: self.temperature_k(p)?, gradient_k_per_m: self.gradient,
        })
    }
}
impl TemperatureSecondOrderField for AffineTemperature {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        Ok(TemperatureJet {
            sample: self.temperature_gradient(p)?,
            diagonal_k_per_m2: [0.0; 3], off_diagonal_k_per_m2: [0.0; 3],
        })
    }
}
fn smooth(
    temperature: &AffineTemperature, ray: Ray, length: f64,
    water: (f64, If97Region), wavelength: f64,
) -> RayBundlePropagation {
    let state = IsobaricWaterField::try_new(
        temperature, Pressure::new(water.0), water.1, ContentHash([106; 32]),
    ).unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    let binding = state.bind_refraction(&absorption).unwrap();
    let prepared = binding.prepare(wavelength).unwrap();
    // A P1 temperature Hessian is zero, but the original water optical Hessian
    // is nonlinear and must remain present in the variational integration.
    assert!(prepared.optical_jet(ray.origin).unwrap().hessian().diagonal()[0].abs() > 1e-10);
    propagate_bundle(
        &prepared, ray, Vec3::new(1.0, 0.0, 0.0), wavelength, length,
        BundleOptions {
            integration: options().mesh.event.integration,
            position_derivative_tolerance_m: options().position_derivative_tolerance_m,
            direction_derivative_tolerance: options().direction_derivative_tolerance,
        },
    ).unwrap()
}

#[test]
fn g1_water_mesh_rays_and_bundles_retain_the_original_field_in_both_phases() {
    with_cx(|_, cx| {
        for reverse in [false, true] {
            let mesh = column(reverse);
            for (pressure, temperature, region) in [
                (3e6, 300.0, If97Region::Region1), (8e4, 400.0, If97Region::Region2),
            ] {
                let original = AffineTemperature { reference: temperature, gradient: Vec3::new(8.0, 0.0, 0.0) };
                let values: Vec<_> = mesh.positions().iter().map(|p| temperature + 8.0 * p[0]).collect();
                let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
                let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
                let ray = Ray { origin: Point3::new(0.17, -0.23, -0.9), dir: Vec3::new(0.0, 0.0, 1.0) };
                for wavelength in [589.0, 900.0] {
                    let result = trace(cx, &prepared, ray, 1.7, (pressure, region), wavelength, options()).unwrap();
                    let expected = smooth(&original, ray, 1.7, (pressure, region), wavelength);
                    same_path(result.bundle.propagation(), expected.propagation(), 5e-9);
                    let central = primal(cx, &prepared, ray, 1.7, (pressure, region), wavelength, options().mesh).unwrap();
                    same_path(&central.propagation, expected.propagation(), 5e-9);
                    assert!(result.internal_faces >= 5 && central.internal_faces >= 5);
                    assert_eq!(result.bundle.interface_events(), 0);
                    for (a, b) in result.bundle.angular_derivatives().into_iter().zip(expected.angular_derivatives()) {
                        vector_near(a.position_m_per_rad, b.position_m_per_rad, 3e-8);
                        vector_near(a.direction_per_rad, b.direction_per_rad, 3e-8);
                    }
                    let footprint = |bundle: &RayBundlePropagation| bundle
                        .project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())
                        .unwrap().area_per_solid_angle_m2();
                    near(footprint(&result.bundle), footprint(&expected), 8e-8);
                    let p = result.bundle.propagation();
                    assert!((p.optical_path_length_m - p.group_optical_path_length_m.unwrap()).abs() > 1e-8);
                }
            }
        }
    });
}

#[test]
fn g1_water_gradient_jump_bundle_matches_neighboring_original_rays() {
    with_cx(|_, cx| {
        let mesh = bipyramid(false);
        let values: Vec<_> = mesh.positions().iter().map(|p| {
            355.0 + 10.0 * p[0] - 4.0 * p[1] + if p[2] < 0.0 { 15.0 * p[2] } else { -20.0 * p[2] }
        }).collect();
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let ray = Ray { origin: Point3::new(-0.06, -0.04, -0.6), dir: unit(Vec3::new(0.18, 0.07, 1.0)) };
        let water = (3e6, If97Region::Region1);
        let result = trace(cx, &prepared, ray, 1.0, water, 589.0, options()).unwrap();
        let central = primal(cx, &prepared, ray, 1.0, water, 589.0, options().mesh).unwrap();
        same_path(result.bundle.propagation(), &central.propagation, 5e-10);
        assert_eq!(result.internal_faces, 1);
        assert_eq!(result.element_index, 1);
        assert_eq!(result.bundle.interface_events(), 0);
        let epsilon = 2e-5;
        for (axis, actual) in result.bundle.launch_basis().into_iter().zip(result.bundle.angular_derivatives()) {
            let neighbour = |sign: f64| {
                let d = Vec3::new(
                    ray.dir.x + sign * epsilon * axis.x,
                    ray.dir.y + sign * epsilon * axis.y,
                    ray.dir.z + sign * epsilon * axis.z,
                );
                primal(cx, &prepared, Ray { dir: unit(d), ..ray }, 1.0, water, 589.0, options().mesh)
                    .unwrap().propagation.ray
            };
            let (a, b) = (neighbour(1.0), neighbour(-1.0));
            vector_near(actual.position_m_per_rad, a.origin.delta_from(b.origin).scale(0.5 / epsilon), 5e-7);
            let direction = Vec3::new(a.dir.x - b.dir.x, a.dir.y - b.dir.y, a.dir.z - b.dir.z).scale(0.5 / epsilon);
            vector_near(actual.direction_per_rad, direction, 5e-7);
        }
    });
}

#[test]
fn g1_accepted_backward_euler_cooling_changes_water_bending_and_bundle_transport() {
    with_cx(|_, cx| {
        let mesh = bipyramid(false);
        let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material, element_materials: None, source: &source };
        let time = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(2.0).unwrap()).unwrap();
        let old = [280.0, 320.0, 300.0, 300.0, 300.0];
        let accepted = time.advance(
            cx, problem, None, &old, 1.0,
            StepConfig { linear: LinearConfig::default(), energy_tolerance_j: 1e-10 },
        ).unwrap();
        for (&actual, expected) in accepted.temperature.iter().zip([290.0, 310.0, 300.0, 300.0, 300.0]) {
            near(actual, expected, 3e-11);
        }
        near(accepted.energy_residual_j, 0.0, 1e-10);
        let ray = Ray { origin: Point3::new(0.0, 0.0, -0.5), dir: Vec3::new(0.0, 0.0, 1.0) };
        let mut bends = [0.0; 2];
        for (i, (values, slope)) in [(&old[..], 20.0), (&accepted.temperature[..], 10.0)].into_iter().enumerate() {
            let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, values).unwrap();
            let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
            assert!(core::ptr::eq(prepared.snapshot().nodal_temperature_k(), values));
            let result = trace(cx, &prepared, ray, 0.9, (3e6, If97Region::Region1), 589.0, options()).unwrap();
            let original = AffineTemperature { reference: 300.0, gradient: Vec3::new(slope, 0.0, 0.0) };
            let expected = smooth(&original, ray, 0.9, (3e6, If97Region::Region1), 589.0);
            same_path(result.bundle.propagation(), expected.propagation(), 5e-9);
            assert_eq!(result.internal_faces, 1);
            for (a, b) in result.bundle.angular_derivatives().into_iter().zip(expected.angular_derivatives()) {
                vector_near(a.position_m_per_rad, b.position_m_per_rad, 3e-8);
                vector_near(a.direction_per_rad, b.direction_per_rad, 3e-8);
            }
            bends[i] = result.bundle.propagation().ray.dir.x;
        }
        assert!(bends[0] < -1e-4);
        assert!(bends[1].abs() > 0.49 * bends[0].abs() && bends[1].abs() < 0.51 * bends[0].abs());
    });
}

#[test]
fn g3_water_mesh_retains_contact_phase_snapshot_and_work_boundaries() {
    with_cx(|gate, cx| {
        let mesh = bipyramid(false);
        let values = [300.0; 5];
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let ray = Ray { origin: Point3::new(0.0, 0.0, -0.5), dir: Vec3::new(0.0, 0.0, 1.0) };
        let water = (3e6, If97Region::Region1);
        let mut limited = options();
        limited.mesh.max_internal_faces = 0;
        assert!(matches!(trace(cx, &prepared, ray, 0.9, water, 589.0, limited),
            Err(WaterMeshError::Mesh(P1MeshError::FaceBudgetExceeded))));
        assert!(matches!(primal(cx, &prepared, ray, 0.9, water, 589.0, limited.mesh),
            Err(WaterMeshError::Mesh(P1MeshError::FaceBudgetExceeded))));

        let copied = values.to_vec();
        let other = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &copied).unwrap();
        let element = other.element(0).unwrap();
        let state = IsobaricWaterField::try_new(&element, Pressure::new(water.0), water.1, ContentHash([107; 32])).unwrap();
        let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
        let binding = state.bind_refraction(&absorption).unwrap();
        assert!(matches!(prepared.propagate_water_bundle(
            cx, &binding, ray, Vec3::new(1.0, 0.0, 0.0), 589.0, 0.9, options(),
        ), Err(WaterMeshError::Mesh(P1MeshError::ForeignSnapshot))));

        let contact = bipyramid(true);
        let contact_values = vec![300.0; contact.vertex_count()];
        let contact_snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &contact, &contact_values).unwrap();
        let contact_mesh = P1ConductionMesh::try_new(cx, contact_snapshot).unwrap();
        let stop = trace(cx, &contact_mesh, ray, 0.9, water, 589.0, options()).unwrap();
        assert!(matches!(stop.termination, BoundaryTermination::Surface(_)));
        assert_eq!(stop.internal_faces, 0);
        assert_eq!(stop.bundle.interface_events(), 0);
        near(stop.bundle.propagation().distance_m, 0.5, 1e-8);

        let crossing = [350.0, 350.0, 350.0, 340.0, 440.0];
        let boiling_snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &crossing).unwrap();
        let boiling_mesh = P1ConductionMesh::try_new(cx, boiling_snapshot).unwrap();
        assert!(matches!(primal(cx, &boiling_mesh, ray, 0.9, (1e5, If97Region::Region1), 589.0, options().mesh),
            Err(WaterMeshError::Mesh(P1MeshError::Propagation(GradientIndexError::OutsideFieldDomain)))));

        gate.request();
        assert!(matches!(trace(cx, &prepared, ray, 0.9, water, 589.0, options()),
            Err(WaterMeshError::Mesh(P1MeshError::Propagation(GradientIndexError::Cancelled)))
            | Err(WaterMeshError::Mesh(P1MeshError::Bundle(BundleError::Propagation(GradientIndexError::Cancelled))))));
    });
}
