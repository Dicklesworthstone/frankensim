use super::*;
use crate::charts::Ray;
use crate::dielectric::{BeerLambertAbsorption, DielectricGlass, GlassProvenance, SellmeierIor};
use crate::gradient_index::bundle::{BundleOptions, propagate_bundle};
use crate::gradient_index::event::{
    BoundaryOracle, BoundarySide, BoundaryTermination, DepartureRequest, GradientDeparture,
    GradientEventOptions, propagate_to_boundary_with_control,
};
use crate::gradient_index::{GradientIndexOptions, propagate};
use crate::material_state::thermo_optic::MaterialThermoOpticField;
use fs_blake3::ContentHash;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{
    ConductionProblem, ConductivityModel, LinearConfig, ScalarField, ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::thermo_optic::{SchottThermoOpticLaw, ThermoOpticCoefficients};
use fs_rep_mesh::TetComplex;

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}
fn mesh(tet: [u32; 4]) -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![tet]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
    .unwrap()
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
fn integration() -> GradientIndexOptions {
    GradientIndexOptions {
        initial_step_m: 0.01,
        max_step_m: 0.04,
        relative_tolerance: 2.0e-12,
        position_tolerance_m: 2.0e-12,
        direction_tolerance: 2.0e-12,
        optical_path_tolerance_m: 2.0e-12,
        optical_depth_tolerance: 2.0e-12,
        ..GradientIndexOptions::default()
    }
}
fn reference() -> DielectricGlass {
    DielectricGlass::new(
        SellmeierIor::try_new_in_band([1.25, 0.0, 0.0], [0.0; 3], (1000.0, 2000.0)).unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}
fn law(d: [f64; 3]) -> SchottThermoOpticLaw {
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d,
            e: [0.0; 2],
            resonance_um: 0.0,
        },
        300.0,
        [200.0, 450.0],
        [1000.0, 2000.0],
        ContentHash([81; 32]),
    )
    .unwrap()
}

#[test]
fn original_affine_geometry_and_both_orientations_retain_exact_p1_state() {
    with_cx(|_, cx| {
        let positions = vec![
            [10.0, -3.0, 5.0],
            [11.0, -2.0, 5.0],
            [10.0, -1.0, 6.0],
            [9.0, -3.0, 7.0],
        ];
        let values: Vec<_> = positions
            .iter()
            .map(|p| 300.0 + 3.0 * p[0] - 5.0 * p[1] + 7.0 * p[2])
            .collect();
        for tet in [[0, 1, 2, 3], [0, 2, 1, 3]] {
            let mesh = ConductionMesh::new(TetComplex::from_tets(4, vec![tet]), positions.clone())
                .unwrap();
            let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
            assert!(core::ptr::eq(snapshot.mesh(), &mesh));
            assert!(core::ptr::eq(
                snapshot.nodal_temperature_k(),
                values.as_slice()
            ));
            let element = snapshot.element(0).unwrap();
            for (v, p) in positions.iter().enumerate() {
                assert_eq!(
                    element
                        .temperature_k(Point3::new(p[0], p[1], p[2]))
                        .unwrap(),
                    values[v]
                );
            }
            let p = Point3::new(10.0, -2.25, 5.75);
            let jet = element.sample_temperature_jet(p).unwrap();
            near(
                jet.sample.temperature_k,
                300.0 + 3.0 * p.x - 5.0 * p.y + 7.0 * p.z,
                1.0e-12,
            );
            vector_near(
                jet.sample.gradient_k_per_m,
                Vec3::new(3.0, -5.0, 7.0),
                2.0e-13,
            );
            assert_eq!(jet.diagonal_k_per_m2, [0.0; 3]);
            assert_eq!(jet.off_diagonal_k_per_m2, [0.0; 3]);
            // Every original face is outward despite the reversed stored tet.
            for face in element.faces_m() {
                let a = Point3::new(face[0][0], face[0][1], face[0][2]);
                let b = Point3::new(face[1][0], face[1][1], face[1][2]);
                let c = Point3::new(face[2][0], face[2][1], face[2][2]);
                let u = Vec3::new(b.x - a.x, b.y - a.y, b.z - a.z);
                let v = Vec3::new(c.x - a.x, c.y - a.y, c.z - a.z);
                let normal = Vec3::new(
                    u.y * v.z - u.z * v.y,
                    u.z * v.x - u.x * v.z,
                    u.x * v.y - u.y * v.x,
                );
                assert!(normal.dot(Vec3::new(p.x - a.x, p.y - a.y, p.z - a.z)) < 0.0);
            }
            assert_eq!(
                element.temperature_k(Point3::new(50.0, 0.0, 0.0)),
                Err(GradientIndexError::OutsideFieldDomain)
            );
        }
    });
}

#[test]
fn original_contact_temperatures_and_gradients_stay_on_their_own_side() {
    with_cx(|_, cx| {
        // Duplicated contact nodes carry two physical temperatures at one face.
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 6, 5, 7]]),
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, -1.0],
            ],
        )
        .unwrap();
        let values = [300.0, 310.0, 300.0, 305.0, 400.0, 380.0, 400.0, 403.0];
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let a = snapshot.element(0).unwrap();
        let b = snapshot.element(1).unwrap();
        let contact = Point3::new(0.2, 0.2, 0.0);
        near(a.temperature_k(contact).unwrap(), 302.0, 1.0e-13);
        near(b.temperature_k(contact).unwrap(), 396.0, 1.0e-13);
        vector_near(
            a.temperature_gradient(contact).unwrap().gradient_k_per_m,
            Vec3::new(10.0, 0.0, 5.0),
            1.0e-13,
        );
        vector_near(
            b.temperature_gradient(contact).unwrap().gradient_k_per_m,
            Vec3::new(-20.0, 0.0, -3.0),
            1.0e-13,
        );
        assert_eq!(
            a.temperature_k(Point3::new(0.2, 0.2, -0.01)),
            Err(GradientIndexError::OutsideFieldDomain)
        );
        assert_eq!(
            b.temperature_k(Point3::new(0.2, 0.2, 0.01)),
            Err(GradientIndexError::OutsideFieldDomain)
        );
    });
}

#[test]
fn state_admission_and_finite_face_oracle_refuse_unsupported_queries() {
    with_cx(|gate, cx| {
        let mesh = mesh([0, 1, 2, 3]);
        assert!(matches!(
            P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &[300.0; 3]),
            Err(TransientTemperatureFieldError::FieldLength { .. })
        ));
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                P1ConductionTemperatureSnapshot::try_new(
                    cx,
                    &mesh,
                    &[300.0, invalid, 300.0, 300.0]
                ),
                Err(TransientTemperatureFieldError::InvalidTemperature { vertex: 1 })
            ));
        }
        let values = [300.0; 4];
        let state = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        assert!(matches!(
            state.element(1),
            Err(TransientTemperatureFieldError::InvalidElement { element: 1 })
        ));
        let element = state.element(0).unwrap();
        assert!(!element.is_spatially_uniform());
        let mut oracle = element.boundary_oracle(17, 1.0e-8).unwrap();
        let centre = oracle.query(Point3::new(0.25, 0.25, 0.25)).unwrap();
        assert!(centre.clearance_m > 0.14 && centre.clearance_m <= 0.25 / 3.0f64.sqrt());
        assert!(centre.witness.is_none());
        let face = oracle
            .query(Point3::new(0.2, 0.2, 1.0e-10))
            .unwrap()
            .witness
            .unwrap();
        assert_eq!(face.boundary_id, 17);
        assert_eq!(face.feature_id, 3);
        assert_eq!(face.side, BoundarySide::Negative);
        vector_near(face.normal, Vec3::new(0.0, 0.0, -1.0), 1.0e-14);
        assert!(face.distance_upper_m <= 1.0e-8);
        // The radius supports finite-face arrival; it authorizes no external
        // continuation, even when the supplied contact is this oracle's own.
        assert!(face.departure_radius_m.unwrap() > 4.0e-8);
        assert_eq!(
            oracle.prepare_departure(DepartureRequest {
                contact: GradientDeparture {
                    witness: face,
                    side: BoundarySide::Negative
                },
                origin: Point3::new(0.2, 0.2, 1.0e-10),
                direction: Vec3::new(0.0, 0.0, 1.0),
                surface_tolerance_m: 1.0e-8,
                min_approach_cosine: 1.0e-5,
            }),
            Err(GradientIndexError::UnsupportedBoundaryDeparture)
        );
        assert_eq!(
            oracle.query(Point3::new(1.0e-11, 0.2, 1.0e-11)),
            Err(GradientIndexError::BoundaryResolutionLimit)
        );
        gate.request();
        assert!(matches!(
            P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values),
            Err(TransientTemperatureFieldError::Cancelled)
        ));
    });
}

#[test]
fn native_backward_euler_cooling_drives_bending_delay_and_absorption_to_first_face() {
    with_cx(|_, cx| {
        let mesh = mesh([0, 1, 2, 3]);
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
            BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(4.0).unwrap())
                .unwrap();
        let old = [300.0, 320.0, 280.0, 300.0];
        let accepted = time
            .advance(
                cx,
                problem,
                None,
                &old,
                1.0,
                StepConfig {
                    linear: LinearConfig::default(),
                    energy_tolerance_j: 1.0e-11,
                },
            )
            .unwrap();
        // Native lumped P1 C=capacity/24 I and K v=(k/6)v for
        // v=[0,1,-1,0]. Backward Euler halves this original thermal mode.
        for (actual, expected) in accepted
            .temperature
            .iter()
            .zip([300.0, 310.0, 290.0, 300.0])
        {
            near(*actual, expected, 2.0e-11);
        }
        near(accepted.stored_energy_change_j, 0.0, 1.0e-11);
        near(accepted.energy_residual_j, 0.0, 1.0e-11);
        assert_eq!(old, [300.0, 320.0, 280.0, 300.0]);
        let reference = reference();
        let law = law([0.003, 0.0, 0.0]);
        let mut bends = [0.0; 2];
        for (i, (values, amplitude)) in [(&old[..], 20.0), (&accepted.temperature[..], 10.0)]
            .into_iter()
            .enumerate()
        {
            let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, values).unwrap();
            let element = snapshot.element(0).unwrap();
            let field = MaterialThermoOpticField::try_new(
                &reference,
                &law,
                &element,
                ContentHash([91 + i as u8; 32]),
            )
            .unwrap();
            let prepared = field.prepare(1500.0).unwrap();
            let tolerance = 1.0e-9;
            let mut oracle = element.boundary_oracle(27, tolerance).unwrap();
            let ray = Ray {
                origin: Point3::new(0.1, 0.1, 0.1),
                dir: Vec3::new(0.0, 0.0, 1.0),
            };
            let path = propagate_to_boundary_with_control(
                &prepared,
                &mut oracle,
                ray,
                Vec3::new(1.0, 0.0, 0.0),
                1500.0,
                1.0,
                None,
                GradientEventOptions {
                    integration: integration(),
                    surface_tolerance_m: tolerance,
                    ..GradientEventOptions::default()
                },
                || cx.checkpoint().is_ok(),
            )
            .unwrap();
            let BoundaryTermination::Surface(contact) = path.termination else {
                panic!("original tetrahedron face was skipped")
            };
            assert_eq!(contact.boundary_id, 27);
            assert_eq!(contact.feature_id, 0);
            // Independent exact affine-index ray: d(n u)/ds=grad(n).
            // n0=1.5, grad(n)=5/12 D0 amplitude (1,-1,0).
            let n0 = 1.5;
            let b = 2.0f64.sqrt() * (5.0 / 12.0) * 0.003 * amplitude;
            let p = path.propagation;
            let s = p.distance_m;
            let n = (n0 * n0 + b * b * s * s).sqrt();
            let z = n0 / b * (b * s / n0).asinh();
            let transverse = b * s * s / (n + n0) / 2.0f64.sqrt();
            near(p.ray.origin.x, 0.1 + transverse, 4.0e-10);
            near(p.ray.origin.y, 0.1 - transverse, 4.0e-10);
            near(p.ray.origin.z, 0.1 + z, 4.0e-10);
            vector_near(
                p.ray.dir,
                Vec3::new(
                    b * s / n / 2.0f64.sqrt(),
                    -b * s / n / 2.0f64.sqrt(),
                    n0 / n,
                ),
                4.0e-10,
            );
            near(s, n0 / b * (b * 0.7 / n0).sinh(), 8.0e-9);
            let phase = 0.5 * (s * n + n0 * n0 / b * (b * s / n0).asinh());
            near(p.optical_path_length_m, phase, 4.0e-10);
            near(p.group_optical_path_length_m.unwrap(), phase, 4.0e-10);
            near(p.power_transmittance(), (-0.2 * s).exp(), 4.0e-11);
            assert!(p.field_evaluations > p.accepted_steps);
            assert!(element.temperature_k(p.ray.origin).is_ok());
            bends[i] = p.ray.dir.x;
        }
        assert!(bends[0] > 0.01 && bends[1] < 0.51 * bends[0] && bends[1] > 0.49 * bends[0]);
    });
}

#[test]
fn nonlinear_thermo_optic_hessian_drives_native_bundle_inside_exact_p1_element() {
    with_cx(|_, cx| {
        let mesh = mesh([0, 1, 2, 3]);
        let values = [300.0, 320.0, 280.0, 300.0];
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let element = snapshot.element(0).unwrap();
        let reference = reference();
        let law = law([0.003, 0.0006, 0.0]);
        let field =
            MaterialThermoOpticField::try_new(&reference, &law, &element, ContentHash([93; 32]))
                .unwrap();
        let prepared = field.prepare(1500.0).unwrap();
        let ray = Ray {
            origin: Point3::new(0.2, 0.2, 0.1),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        assert_eq!(
            element
                .temperature_jet(ray.origin)
                .unwrap()
                .diagonal_k_per_m2,
            [0.0; 3]
        );
        let hessian = prepared.optical_jet(ray.origin).unwrap().hessian();
        // n_TT=(n-1/n) D1=0.0005; grad(T)=(20,-20,0).
        for (actual, expected) in hessian.diagonal().into_iter().zip([0.2, 0.2, 0.0]) {
            near(actual, expected, 1.0e-14);
        }
        for (actual, expected) in hessian.off_diagonal().into_iter().zip([-0.2, 0.0, 0.0]) {
            near(actual, expected, 1.0e-14);
        }
        let length = 0.3;
        let bundle = propagate_bundle(
            &prepared,
            ray,
            Vec3::new(1.0, 0.0, 0.0),
            1500.0,
            length,
            BundleOptions {
                integration: integration(),
                position_derivative_tolerance_m: 2.0e-12,
                direction_derivative_tolerance: 2.0e-12,
            },
        )
        .unwrap();
        // Independently vary both original launch angles in the native primal.
        let epsilon = 2.0e-5;
        for (axis, derivative) in bundle
            .launch_basis()
            .into_iter()
            .zip(bundle.angular_derivatives())
        {
            let trace = |sign: f64| {
                let dir = Vec3::new(
                    sign * epsilon * axis.x,
                    sign * epsilon * axis.y,
                    1.0 + sign * epsilon * axis.z,
                );
                propagate(
                    &prepared,
                    Ray {
                        origin: ray.origin,
                        dir,
                    },
                    Vec3::new(1.0, 0.0, 0.0),
                    1500.0,
                    length,
                    integration(),
                )
                .unwrap()
            };
            let plus = trace(1.0);
            let minus = trace(-1.0);
            let scale = 0.5 / epsilon;
            vector_near(
                derivative.position_m_per_rad,
                Vec3::new(
                    (plus.ray.origin.x - minus.ray.origin.x) * scale,
                    (plus.ray.origin.y - minus.ray.origin.y) * scale,
                    (plus.ray.origin.z - minus.ray.origin.z) * scale,
                ),
                2.0e-9,
            );
            vector_near(
                derivative.direction_per_rad,
                Vec3::new(
                    (plus.ray.dir.x - minus.ray.dir.x) * scale,
                    (plus.ray.dir.y - minus.ray.dir.y) * scale,
                    (plus.ray.dir.z - minus.ray.dir.z) * scale,
                ),
                2.0e-9,
            );
        }
        assert!((bundle.angular_derivatives()[0].position_m_per_rad.x - length).abs() > 1.0e-4);
    });
}
