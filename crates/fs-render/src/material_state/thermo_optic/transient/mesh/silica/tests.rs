use super::*;
use crate::dielectric::BeerLambertAbsorption;
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::bundle::{BundleOptions, propagate_bundle};
use crate::material_state::thermo_optic::TemperatureSecondOrderField;
use fs_blake3::ContentHash;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig, ScalarField,
    ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7980,
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
fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.dot(v).sqrt())
}
fn mesh(contact: bool) -> ConductionMesh {
    let mut positions = vec![
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 1.0],
    ];
    let tets = if contact {
        positions.extend([positions[0], positions[1], positions[2]]);
        vec![[0, 1, 2, 3], [5, 6, 7, 4]]
    } else {
        vec![[0, 1, 2, 3], [0, 1, 2, 4]]
    };
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}
fn options() -> P1MeshBundleOptions {
    P1MeshBundleOptions {
        mesh: P1MeshOptions {
            boundary_id: 7980,
            max_internal_faces: 8,
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
fn primal(
    cx: &Cx<'_>,
    mesh: &P1ConductionMesh<'_>,
    ray: Ray,
    length: f64,
    wavelength: f64,
    options: P1MeshOptions,
) -> Result<P1MeshPropagation, SilicaMeshError> {
    let element = mesh.snapshot().element(0).unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    let binding =
        MaterialSilicaRefractiveField::try_new(&element, &absorption, ContentHash([79; 32]))
            .unwrap();
    mesh.propagate_silica(
        cx,
        &binding,
        ray,
        Vec3::new(1.0, 0.0, 0.0),
        wavelength,
        length,
        options,
    )
}
fn bundle(
    cx: &Cx<'_>,
    mesh: &P1ConductionMesh<'_>,
    ray: Ray,
    length: f64,
    wavelength: f64,
) -> P1MeshBundlePropagation {
    let element = mesh.snapshot().element(0).unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    let binding =
        MaterialSilicaRefractiveField::try_new(&element, &absorption, ContentHash([79; 32]))
            .unwrap();
    mesh.propagate_silica_bundle(
        cx,
        &binding,
        ray,
        Vec3::new(1.0, 0.0, 0.0),
        wavelength,
        length,
        options(),
    )
    .unwrap()
}
struct Affine {
    base: f64,
    gradient: Vec3,
}
impl TemperatureField for Affine {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        Ok(self.base + self.gradient.dot(Vec3::new(p.x, p.y, p.z)))
    }
}
impl TemperatureGradientField for Affine {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        Ok(TemperatureGradientSample {
            temperature_k: self.temperature_k(p)?,
            gradient_k_per_m: self.gradient,
        })
    }
}
impl TemperatureSecondOrderField for Affine {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        Ok(TemperatureJet {
            sample: self.temperature_gradient(p)?,
            diagonal_k_per_m2: [0.0; 3],
            off_diagonal_k_per_m2: [0.0; 3],
        })
    }
}
fn same_path(a: &GradientIndexPropagation, b: &GradientIndexPropagation) {
    vector_near(
        a.ray.origin.delta_from(b.ray.origin),
        Vec3::new(0.0, 0.0, 0.0),
        5e-9,
    );
    vector_near(a.ray.dir, b.ray.dir, 5e-9);
    near(a.distance_m, b.distance_m, 5e-13);
    near(a.optical_path_length_m, b.optical_path_length_m, 5e-9);
    near(
        a.group_optical_path_length_m.unwrap(),
        b.group_optical_path_length_m.unwrap(),
        5e-9,
    );
    near(a.power_transmittance(), (-0.02 * a.distance_m).exp(), 5e-13);
}

#[test]
fn g1_accepted_cooling_changes_native_silica_bending_and_bundle_curvature() {
    with_cx(|_, cx| {
        let mesh = mesh(false);
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
        let old = [200.0, 240.0, 220.0, 220.0, 220.0];
        let accepted = time
            .advance(
                cx,
                problem,
                None,
                &old,
                1.0,
                StepConfig {
                    linear: LinearConfig::default(),
                    energy_tolerance_j: 1e-10,
                },
            )
            .unwrap();
        for (&actual, expected) in accepted
            .temperature
            .iter()
            .zip([210.0, 230.0, 220.0, 220.0, 220.0])
        {
            near(actual, expected, 3e-11);
        }
        near(accepted.energy_residual_j, 0.0, 1e-10);
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, -0.5),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
        let mut bends = [0.0; 2];
        for (i, (values, slope)) in [(&old[..], 20.0), (&accepted.temperature[..], 10.0)]
            .into_iter()
            .enumerate()
        {
            let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, values).unwrap();
            let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
            assert!(core::ptr::eq(
                prepared.snapshot().nodal_temperature_k(),
                values
            ));
            let affine = Affine {
                base: 220.0,
                gradient: Vec3::new(slope, 0.0, 0.0),
            };
            let original =
                MaterialSilicaRefractiveField::try_new(&affine, &absorption, ContentHash([80; 32]))
                    .unwrap();
            for wavelength in [550.0, 1550.0] {
                let optical = original.prepare(wavelength).unwrap();
                assert!(
                    optical
                        .optical_jet(ray.origin)
                        .unwrap()
                        .hessian()
                        .diagonal()[0]
                        .abs()
                        > 1e-12
                );
                let expected = propagate_bundle(
                    &optical,
                    ray,
                    Vec3::new(1.0, 0.0, 0.0),
                    wavelength,
                    0.9,
                    BundleOptions {
                        integration: options().mesh.event.integration,
                        position_derivative_tolerance_m: 2e-12,
                        direction_derivative_tolerance: 2e-12,
                    },
                )
                .unwrap();
                let actual = bundle(cx, &prepared, ray, 0.9, wavelength);
                let central = primal(cx, &prepared, ray, 0.9, wavelength, options().mesh).unwrap();
                same_path(actual.bundle.propagation(), expected.propagation());
                same_path(&central.propagation, expected.propagation());
                assert_eq!(actual.internal_faces, 1);
                assert_eq!(actual.bundle.interface_events(), 0);
                for (a, b) in actual
                    .bundle
                    .angular_derivatives()
                    .into_iter()
                    .zip(expected.angular_derivatives())
                {
                    vector_near(a.position_m_per_rad, b.position_m_per_rad, 3e-8);
                    vector_near(a.direction_per_rad, b.direction_per_rad, 3e-8);
                }
                bends[i] = actual.bundle.propagation().ray.dir.x;
            }
        }
        assert!(bends[0] > 1e-5);
        assert!(bends[1] > 0.49 * bends[0] && bends[1] < 0.51 * bends[0]);
    });
}

#[test]
fn g1_silica_gradient_jump_bundle_matches_actual_neighboring_rays() {
    with_cx(|_, cx| {
        let mesh = mesh(false);
        let values: Vec<_> = mesh
            .positions()
            .iter()
            .map(|p| {
                220.0 + 10.0 * p[0] - 4.0 * p[1]
                    + if p[2] < 0.0 {
                        15.0 * p[2]
                    } else {
                        -20.0 * p[2]
                    }
            })
            .collect();
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let ray = Ray {
            origin: Point3::new(-0.06, -0.04, -0.6),
            dir: unit(Vec3::new(0.18, 0.07, 1.0)),
        };
        let actual = bundle(cx, &prepared, ray, 1.0, 1550.0);
        same_path(
            actual.bundle.propagation(),
            &primal(cx, &prepared, ray, 1.0, 1550.0, options().mesh)
                .unwrap()
                .propagation,
        );
        assert_eq!(actual.internal_faces, 1);
        assert_eq!(actual.element_index, 1);
        assert_eq!(actual.bundle.interface_events(), 0);
        let epsilon = 2e-5;
        for (axis, derivative) in actual
            .bundle
            .launch_basis()
            .into_iter()
            .zip(actual.bundle.angular_derivatives())
        {
            let neighbour = |sign: f64| {
                primal(
                    cx,
                    &prepared,
                    Ray {
                        dir: unit(Vec3::new(
                            ray.dir.x + sign * epsilon * axis.x,
                            ray.dir.y + sign * epsilon * axis.y,
                            ray.dir.z + sign * epsilon * axis.z,
                        )),
                        ..ray
                    },
                    1.0,
                    1550.0,
                    options().mesh,
                )
                .unwrap()
                .propagation
                .ray
            };
            let (a, b) = (neighbour(1.0), neighbour(-1.0));
            vector_near(
                derivative.position_m_per_rad,
                a.origin.delta_from(b.origin).scale(0.5 / epsilon),
                5e-7,
            );
            vector_near(
                derivative.direction_per_rad,
                Vec3::new(a.dir.x - b.dir.x, a.dir.y - b.dir.y, a.dir.z - b.dir.z)
                    .scale(0.5 / epsilon),
                5e-7,
            );
        }
    });
}

#[test]
fn g3_silica_mesh_preserves_source_contact_snapshot_and_work_boundaries() {
    with_cx(|gate, cx| {
        let original_mesh = mesh(false);
        let values = [220.0; 5];
        let snapshot =
            P1ConductionTemperatureSnapshot::try_new(cx, &original_mesh, &values).unwrap();
        let prepared = P1ConductionMesh::try_new(cx, snapshot).unwrap();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, -0.5),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let mut limited = options().mesh;
        limited.max_internal_faces = 0;
        assert!(matches!(
            primal(cx, &prepared, ray, 0.9, 550.0, limited),
            Err(SilicaMeshError::Mesh(P1MeshError::FaceBudgetExceeded))
        ));
        let copy = values.to_vec();
        let foreign = P1ConductionTemperatureSnapshot::try_new(cx, &original_mesh, &copy).unwrap();
        let element = foreign.element(0).unwrap();
        let absorption = BeerLambertAbsorption::CLEAR;
        let binding =
            MaterialSilicaRefractiveField::try_new(&element, &absorption, ContentHash([81; 32]))
                .unwrap();
        assert!(matches!(
            prepared.propagate_silica_bundle(
                cx,
                &binding,
                ray,
                Vec3::new(1.0, 0.0, 0.0),
                550.0,
                0.9,
                options()
            ),
            Err(SilicaMeshError::Mesh(P1MeshError::ForeignSnapshot))
        ));
        let contact = mesh(true);
        let contact_values = vec![220.0; contact.vertex_count()];
        let contact_snapshot =
            P1ConductionTemperatureSnapshot::try_new(cx, &contact, &contact_values).unwrap();
        let contact_mesh = P1ConductionMesh::try_new(cx, contact_snapshot).unwrap();
        let stop = bundle(cx, &contact_mesh, ray, 0.9, 550.0);
        assert!(matches!(stop.termination, BoundaryTermination::Surface(_)));
        assert_eq!(stop.internal_faces, 0);
        near(stop.bundle.propagation().distance_m, 0.5, 1e-8);
        let crossing = [220.0, 220.0, 220.0, 210.0, 450.0];
        let out_of_band =
            P1ConductionTemperatureSnapshot::try_new(cx, &original_mesh, &crossing).unwrap();
        let out_of_band = P1ConductionMesh::try_new(cx, out_of_band).unwrap();
        assert!(matches!(
            primal(cx, &out_of_band, ray, 1.1, 550.0, options().mesh),
            Err(SilicaMeshError::Mesh(P1MeshError::Propagation(
                GradientIndexError::OutsideFieldDomain
            )))
        ));
        assert!(matches!(
            primal(cx, &prepared, ray, 0.9, 399.0, options().mesh),
            Err(SilicaMeshError::Optical(_))
        ));
        gate.request();
        assert!(matches!(
            primal(cx, &prepared, ray, 0.9, 550.0, options().mesh),
            Err(SilicaMeshError::Mesh(P1MeshError::Propagation(
                GradientIndexError::Cancelled
            )))
        ));
    });
}
