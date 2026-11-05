use super::*;
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::bundle::{BundleOptions, propagate_bundle};
use crate::gradient_index::thermal::TemperatureField;
use crate::material_state::thermo_optic::transient::mesh::water::pressure_temperature::{
    NodalPressureValues, P1PressureSnapshot, P1PressureTemperatureSnapshot,
};
use crate::material_state::thermo_optic::{
    TemperatureGradientSample, TemperatureJet, TemperatureSecondOrderField,
};
use fs_blake3::ContentHash;
use fs_conduction::ConductionMesh;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::gas::GasSpec;
use fs_material::gas_scattering::{MolecularRayleigh, MolecularRayleighMixture};
use fs_qty::semantic::{Composition, CompositionBasis};
use fs_rep_mesh::TetComplex;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 173,
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
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.dot(v).sqrt())
}
fn model() -> MolecularRayleighMixture {
    MolecularRayleighMixture::try_new(
        &[GasSpec::dry_air_ussa1976()],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        &[
            MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, 0.03, ContentHash([173; 32]))
                .unwrap(),
        ],
        None,
    )
    .unwrap()
}
fn bipyramid(contact: bool, reverse: bool) -> ConductionMesh {
    let mut positions = vec![
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 1.0],
    ];
    let mut tets = if contact {
        positions.extend([positions[0], positions[1], positions[2]]);
        vec![[0, 1, 2, 3], [5, 6, 7, 4]]
    } else {
        vec![[0, 1, 2, 3], [0, 1, 2, 4]]
    };
    if reverse {
        for tet in &mut tets {
            tet.swap(0, 1);
        }
    }
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}
fn snapshot<'a>(
    cx: &Cx<'_>,
    mesh: &'a ConductionMesh,
    temperature: &'a [f64],
    pressure: &'a [f64],
) -> P1PressureTemperatureSnapshot<'a> {
    let thermal = P1ConductionTemperatureSnapshot::try_new(cx, mesh, temperature).unwrap();
    let pressure = P1PressureSnapshot::try_new(
        cx,
        mesh,
        NodalPressureValues::Pascals(pressure),
        ContentHash([174; 32]),
    )
    .unwrap();
    P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([175; 32]), pressure).unwrap()
}
fn state<'a, 'mesh>(
    element: &'a P1PressureTemperatureElement<'mesh>,
    model: &'a MolecularRayleighMixture,
) -> PressureTemperatureGasField<'a, P1PressureTemperatureElement<'mesh>> {
    PressureTemperatureGasField::try_new(
        element,
        model,
        Temperature::new(300.0),
        Pressure::new(100_000.0),
        element.snapshot().content_identity(),
    )
    .unwrap()
}
fn options() -> P1MeshBundleOptions {
    P1MeshBundleOptions {
        mesh: P1MeshOptions {
            boundary_id: 173,
            max_internal_faces: 32,
            event: GradientEventOptions {
                integration: GradientIndexOptions {
                    initial_step_m: 0.01,
                    max_step_m: 0.05,
                    relative_tolerance: 2e-12,
                    position_tolerance_m: 2e-12,
                    direction_tolerance: 2e-12,
                    optical_path_tolerance_m: 2e-12,
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
fn ray() -> Ray {
    Ray {
        origin: Point3::new(-0.06, -0.04, -0.6),
        dir: unit(Vec3::new(0.18, 0.07, 1.0)),
    }
}
fn same_path(a: &GradientIndexPropagation, b: &GradientIndexPropagation, tolerance: f64) {
    vector_near(
        a.ray.origin.delta_from(b.ray.origin),
        Vec3::new(0.0, 0.0, 0.0),
        tolerance,
    );
    vector_near(a.ray.dir, b.ray.dir, tolerance);
    near(a.optical_path_length_m, b.optical_path_length_m, tolerance);
    near(
        a.group_optical_path_length_m.unwrap(),
        b.group_optical_path_length_m.unwrap(),
        tolerance,
    );
    assert_eq!(a.power_transmittance(), 1.0);
}

struct AffineTemperature {
    gradient: Vec3,
}
impl TemperatureField for AffineTemperature {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        Ok(300.0 + self.gradient.dot(Vec3::new(p.x, p.y, p.z)))
    }
}
impl TemperatureGradientField for AffineTemperature {
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
impl TemperatureSecondOrderField for AffineTemperature {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        Ok(TemperatureJet {
            sample: self.temperature_gradient(p)?,
            diagonal_k_per_m2: [0.0; 3],
            off_diagonal_k_per_m2: [0.0; 3],
        })
    }
}

#[test]
fn g1_native_mesh_gas_matches_unsplit_thermal_rays_and_angular_derivatives() {
    with_cx(|_, cx| {
        let model = model();
        let analytic = AffineTemperature {
            gradient: Vec3::new(12.0, -7.0, 18.0),
        };
        let original = IsobaricGasField::try_new(
            &analytic,
            &model,
            Pressure::new(100_000.0),
            Temperature::new(300.0),
            ContentHash([176; 32]),
        )
        .unwrap();
        let optics = original.bind_refraction().unwrap();
        for reverse in [false, true] {
            let geometry = bipyramid(false, reverse);
            let temperature: Vec<_> = geometry
                .positions()
                .iter()
                .map(|p| {
                    analytic
                        .temperature_k(Point3::new(p[0], p[1], p[2]))
                        .unwrap()
                })
                .collect();
            let pressure = vec![100_000.0; temperature.len()];
            let joint = snapshot(cx, &geometry, &temperature, &pressure);
            let mesh = P1ConductionMesh::try_new(cx, joint.temperature_snapshot()).unwrap();
            let element = joint.element(0).unwrap();
            let thermal_element = element.temperature_element();
            let gas = IsobaricGasField::try_new(
                &thermal_element,
                &model,
                Pressure::new(100_000.0),
                Temperature::new(300.0),
                ContentHash([176; 32]),
            )
            .unwrap();
            let actual = mesh
                .propagate_gas_bundle(
                    cx,
                    &gas,
                    ray(),
                    Vec3::new(1.0, 0.0, 0.0),
                    550.0,
                    1.0,
                    options(),
                )
                .unwrap();
            let expected = propagate_bundle(
                &optics,
                ray(),
                Vec3::new(1.0, 0.0, 0.0),
                550.0,
                1.0,
                BundleOptions {
                    integration: options().mesh.event.integration,
                    position_derivative_tolerance_m: options().position_derivative_tolerance_m,
                    direction_derivative_tolerance: options().direction_derivative_tolerance,
                },
            )
            .unwrap();
            same_path(actual.bundle.propagation(), expected.propagation(), 5e-10);
            for (a, b) in actual
                .bundle
                .angular_derivatives()
                .into_iter()
                .zip(expected.angular_derivatives())
            {
                vector_near(a.position_m_per_rad, b.position_m_per_rad, 5e-9);
                vector_near(a.direction_per_rad, b.direction_per_rad, 5e-9);
            }
            assert_eq!(actual.internal_faces, 1);
            assert_eq!(actual.bundle.interface_events(), 0);
            let primal = mesh
                .propagate_gas(
                    cx,
                    &gas,
                    ray(),
                    Vec3::new(1.0, 0.0, 0.0),
                    550.0,
                    1.0,
                    options().mesh,
                )
                .unwrap();
            same_path(&primal.propagation, expected.propagation(), 5e-10);
            let coupled = mesh
                .propagate_pressure_temperature_gas_bundle(
                    cx,
                    &state(&element, &model),
                    ray(),
                    Vec3::new(1.0, 0.0, 0.0),
                    550.0,
                    1.0,
                    options(),
                )
                .unwrap();
            same_path(coupled.bundle.propagation(), expected.propagation(), 5e-10);
        }
    });
}

#[test]
fn g1_pressure_only_jump_preserves_receiving_angular_family() {
    with_cx(|_, cx| {
        let geometry = bipyramid(false, false);
        let temperature = [300.0; 5];
        let pressure: Vec<_> = geometry
            .positions()
            .iter()
            .map(|p| {
                100_000.0 + 4_000.0 * p[0] - 2_000.0 * p[1]
                    + if p[2] < 0.0 {
                        15_000.0 * p[2]
                    } else {
                        95_000.0 * p[2]
                    }
            })
            .collect();
        let joint = snapshot(cx, &geometry, &temperature, &pressure);
        let mesh = P1ConductionMesh::try_new(cx, joint.temperature_snapshot()).unwrap();
        let element = joint.element(0).unwrap();
        let model = model();
        let gas = state(&element, &model);
        let result = mesh
            .propagate_pressure_temperature_gas_bundle(
                cx,
                &gas,
                ray(),
                Vec3::new(1.0, 0.0, 0.0),
                550.0,
                1.0,
                options(),
            )
            .unwrap();
        assert_eq!(result.internal_faces, 1);
        assert_eq!(result.element_index, 1);
        assert_eq!(result.bundle.interface_events(), 0);
        let epsilon = 2e-5;
        for (basis, actual) in result
            .bundle
            .launch_basis()
            .into_iter()
            .zip(result.bundle.angular_derivatives())
        {
            let neighbour = |sign: f64| {
                let d = ray().dir;
                let direction = unit(Vec3::new(
                    d.x + sign * epsilon * basis.x,
                    d.y + sign * epsilon * basis.y,
                    d.z + sign * epsilon * basis.z,
                ));
                mesh.propagate_pressure_temperature_gas(
                    cx,
                    &gas,
                    Ray {
                        dir: direction,
                        ..ray()
                    },
                    Vec3::new(1.0, 0.0, 0.0),
                    550.0,
                    1.0,
                    options().mesh,
                )
                .unwrap()
                .propagation
                .ray
            };
            let (a, b) = (neighbour(1.0), neighbour(-1.0));
            vector_near(
                actual.position_m_per_rad,
                a.origin.delta_from(b.origin).scale(0.5 / epsilon),
                5e-7,
            );
            vector_near(
                actual.direction_per_rad,
                Vec3::new(a.dir.x - b.dir.x, a.dir.y - b.dir.y, a.dir.z - b.dir.z)
                    .scale(0.5 / epsilon),
                5e-7,
            );
        }
    });
}

#[test]
fn g1_correlated_native_pressure_and_temperature_keep_constant_index_across_a_jump() {
    with_cx(|_, cx| {
        let geometry = bipyramid(false, false);
        let temperature: Vec<_> = geometry
            .positions()
            .iter()
            .map(|p| {
                300.0
                    + if p[2] < 0.0 {
                        20.0 * p[2]
                    } else {
                        120.0 * p[2]
                    }
            })
            .collect();
        let pressure: Vec<_> = temperature
            .iter()
            .map(|&t| (100_000.0 / 300.0) * t)
            .collect();
        let joint = snapshot(cx, &geometry, &temperature, &pressure);
        let mesh = P1ConductionMesh::try_new(cx, joint.temperature_snapshot()).unwrap();
        let element = joint.element(0).unwrap();
        let model = model();
        let gas = state(&element, &model);
        let result = mesh
            .propagate_pressure_temperature_gas_bundle(
                cx,
                &gas,
                ray(),
                Vec3::new(1.0, 0.0, 0.0),
                550.0,
                1.0,
                options(),
            )
            .unwrap();
        let propagation = result.bundle.propagation();
        vector_near(
            propagation.ray.origin.delta_from(ray().origin),
            ray().dir,
            5e-10,
        );
        vector_near(propagation.ray.dir, ray().dir, 5e-10);
        let n = model
            .evaluate(300.0, 100_000.0)
            .unwrap()
            .refractive_index(550.0)
            .unwrap()
            .phase_index;
        near(propagation.optical_path_length_m, n, 5e-10);
        near(propagation.group_optical_path_length_m.unwrap(), n, 5e-10);
        for (basis, derivative) in result
            .bundle
            .launch_basis()
            .into_iter()
            .zip(result.bundle.angular_derivatives())
        {
            vector_near(derivative.position_m_per_rad, basis, 5e-9);
            vector_near(derivative.direction_per_rad, basis, 5e-9);
        }
        assert_eq!(result.internal_faces, 1);
        assert_eq!(result.bundle.interface_events(), 0);
    });
}

#[test]
fn g3_gas_mesh_preserves_native_sources_contacts_work_and_cancellation() {
    with_cx(|gate, cx| {
        let geometry = bipyramid(false, false);
        let temperature = [300.0; 5];
        let pressure = [100_000.0; 5];
        let joint = snapshot(cx, &geometry, &temperature, &pressure);
        let mesh = P1ConductionMesh::try_new(cx, joint.temperature_snapshot()).unwrap();
        let element = joint.element(0).unwrap();
        let model = model();
        let gas = state(&element, &model);
        let axis = Vec3::new(1.0, 0.0, 0.0);
        let wrong = PressureTemperatureGasField::try_new(
            &element,
            &model,
            Temperature::new(300.0),
            Pressure::new(100_000.0),
            ContentHash([177; 32]),
        )
        .unwrap();
        assert!(matches!(
            mesh.propagate_pressure_temperature_gas(
                cx,
                &wrong,
                ray(),
                axis,
                550.0,
                1.0,
                options().mesh,
            ),
            Err(GasMeshError::SourceIdentityMismatch)
        ));
        let copied = temperature.to_vec();
        let foreign = snapshot(cx, &geometry, &copied, &pressure)
            .element(0)
            .unwrap();
        assert!(matches!(
            mesh.propagate_pressure_temperature_gas_bundle(
                cx,
                &state(&foreign, &model),
                ray(),
                axis,
                550.0,
                1.0,
                options(),
            ),
            Err(GasMeshError::Mesh(P1MeshError::ForeignSnapshot))
        ));
        let mut limited = options();
        limited.mesh.max_internal_faces = 0;
        assert!(matches!(
            mesh.propagate_pressure_temperature_gas_bundle(
                cx,
                &gas,
                ray(),
                axis,
                550.0,
                1.0,
                limited,
            ),
            Err(GasMeshError::Mesh(P1MeshError::FaceBudgetExceeded))
        ));
        let contact = bipyramid(true, false);
        let contact_temperature = vec![300.0; contact.vertex_count()];
        let mut contact_pressure = vec![100_000.0; contact.vertex_count()];
        for vertex in [4, 5, 6, 7] {
            contact_pressure[vertex] = 150_000.0;
        }
        let contact_joint = snapshot(cx, &contact, &contact_temperature, &contact_pressure);
        let contact_mesh =
            P1ConductionMesh::try_new(cx, contact_joint.temperature_snapshot()).unwrap();
        let contact_element = contact_joint.element(0).unwrap();
        let stop = contact_mesh
            .propagate_pressure_temperature_gas_bundle(
                cx,
                &state(&contact_element, &model),
                ray(),
                axis,
                550.0,
                1.0,
                options(),
            )
            .unwrap();
        assert!(matches!(stop.termination, BoundaryTermination::Surface(_)));
        assert_eq!(stop.internal_faces, 0);
        assert_eq!(stop.bundle.interface_events(), 0);
        gate.request();
        assert!(
            mesh.propagate_pressure_temperature_gas_bundle(
                cx,
                &gas,
                ray(),
                axis,
                550.0,
                1.0,
                options(),
            )
            .is_err()
        );
    });
}
