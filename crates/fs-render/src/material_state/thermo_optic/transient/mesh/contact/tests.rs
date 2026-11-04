use super::*;
use crate::dielectric::{BeerLambertAbsorption, DielectricGlass, GlassProvenance, SellmeierIor};
use fs_blake3::ContentHash;
use fs_conduction::ConductionMesh;
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
                seed: 113,
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
fn vector(a: Vec3, b: Vec3) {
    near(a.x, b.x, 3e-14);
    near(a.y, b.y, 3e-14);
    near(a.z, b.z, 3e-14);
}
fn glass(n: f64, absorption: f64) -> DielectricGlass {
    // Phase n at 1500 nm and distinct group index n + (n²-1)/(24 n).
    // The original Sellmeier pole remains outside the supported band.
    DielectricGlass::new(
        SellmeierIor::try_new_in_band(
            [0.96 * (n * n - 1.0), 0.0, 0.0],
            [0.09, 0.0, 0.0],
            (1000.0, 2000.0),
        )
        .unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
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
        ContentHash([113; 32]),
    )
    .unwrap()
}
fn pair(reverse: bool) -> ConductionMesh {
    let positions = vec![
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut tets = vec![[0, 1, 2, 3], [4, 5, 6, 7]];
    if reverse {
        for tet in &mut tets {
            tet.swap(0, 1);
        }
    }
    ConductionMesh::new(TetComplex::from_tets(8, tets), positions).unwrap()
}
fn fixture(
    reverse: bool,
    temperature: [f64; 2],
    indices: [f64; 2],
    absorption: [f64; 2],
    f: impl FnOnce(&CancelGate, &Cx<'_>, &P1ContactInterface<'_>),
) {
    with_cx(|gate, cx| {
        let mesh = pair(reverse);
        let values = [
            temperature[0],
            temperature[0],
            temperature[0],
            temperature[0],
            temperature[1],
            temperature[1],
            temperature[1],
            temperature[1],
        ];
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let a = snapshot.element(0).unwrap();
        let b = snapshot.element(1).unwrap();
        let ga = glass(indices[0], absorption[0]);
        let gb = glass(indices[1], absorption[1]);
        let law = law();
        let ba = MaterialThermoOpticField::try_new(&ga, &law, &a, ContentHash([114; 32])).unwrap();
        let bb = MaterialThermoOpticField::try_new(&gb, &law, &b, ContentHash([114; 32])).unwrap();
        let contact = P1ContactInterface::try_new(cx, [&ba, &bb], [3, 3]).unwrap();
        f(gate, cx, &contact);
    });
}
fn launch(sine: f64, from_top: bool) -> Ray {
    Ray {
        origin: Point3::new(0.0, -0.1, if from_top { 0.2 } else { -0.2 }),
        dir: Vec3::new(
            sine,
            0.0,
            (1.0 - sine * sine).sqrt() * if from_top { -1.0 } else { 1.0 },
        ),
    }
}
const AXIS: Vec3 = Vec3::new(0.0, 1.0, 0.0);

#[test]
fn g1_native_snell_split_preserves_each_materials_path_and_power() {
    for reverse in [false, true] {
        fixture(
            reverse,
            [300.0; 2],
            [1.2, 1.8],
            [0.2, 0.6],
            |_, cx, contact| {
                let ray = launch(0.6, false);
                let split = contact.split(cx, 0, ray, AXIS, 1500.0, 0.1).unwrap();
                let incident_length = 0.25;
                let incident = split.incident;
                near(incident.distance_m, incident_length, 3e-15);
                near(incident.ray.origin.z, 0.0, 3e-15);
                near(incident.ray.origin.x, 0.15, 3e-15);
                let r = split.reflected.unwrap();
                let t = split.transmitted.unwrap();
                vector(r.propagation().ray.dir, Vec3::new(0.6, 0.0, -0.8));
                vector(
                    t.propagation().ray.dir,
                    Vec3::new(0.4, 0.0, 0.84_f64.sqrt()),
                );
                near(
                    t.propagation().optical_path_length_m,
                    1.2 * incident_length + 1.8 * 0.1,
                    3e-15,
                );
                near(
                    t.propagation().group_optical_path_length_m.unwrap(),
                    (1.2 + (1.2 * 1.2 - 1.0) / (24.0 * 1.2)) * incident_length
                        + (1.8 + (1.8 * 1.8 - 1.0) / (24.0 * 1.8)) * 0.1,
                    3e-15,
                );
                near(
                    t.propagation().optical_depth,
                    0.2 * incident_length + 0.6 * 0.1,
                    3e-15,
                );
                near(
                    r.propagation().optical_depth,
                    0.2 * (incident_length + 0.1),
                    3e-15,
                );
                let ct = 0.84_f64.sqrt();
                let rs = (1.2 * 0.8 - 1.8 * ct) / (1.2 * 0.8 + 1.8 * ct);
                let rp = (1.8 * 0.8 - 1.2 * ct) / (1.8 * 0.8 + 1.2 * ct);
                let reflection = 0.5 * (rs * rs + rp * rp);
                near(split.fresnel.reflectance, reflection, 3e-15);
                near(
                    t.power_fraction(),
                    (1.0 - reflection) * (-0.11_f64).exp(),
                    3e-15,
                );
                near(r.power_fraction(), reflection * (-0.07_f64).exp(), 3e-15);
                assert_eq!((r.element_index(), t.element_index()), (0, 1));
                assert_eq!(t.arrival_face(), None);
                vector(t.propagation().transverse_axis, AXIS);
                // Arrival admission is bounded by the requested outgoing segment:
                // a distant edge hit cannot invalidate a valid shorter path.
                contact
                    .split(cx, 0, launch(0.0, false), AXIS, 1500.0, 0.01)
                    .unwrap();
            },
        );
    }
}

#[test]
fn g1_original_thermal_contact_jump_changes_refraction_without_inventing_absorption() {
    fixture(
        false,
        [320.0, 280.0],
        [1.5; 2],
        [0.0; 2],
        |_, cx, contact| {
            let split = contact
                .split(cx, 0, launch(0.4, false), AXIS, 1500.0, 0.1)
                .unwrap();
            let ni = 1.5 + 0.5 * (1.5 - 1.0 / 1.5) * 0.003 * 20.0;
            let nt = 1.5 - 0.5 * (1.5 - 1.0 / 1.5) * 0.003 * 20.0;
            let transmitted = split.transmitted.unwrap();
            near(transmitted.propagation().ray.dir.x, ni / nt * 0.4, 3e-15);
            near(
                transmitted.power_fraction() + split.reflected.unwrap().power_fraction(),
                1.0,
                3e-15,
            );
            near(transmitted.propagation().optical_depth, 0.0, 0.0);
            // Reverse the actual refracted direction in the same original state.
            let direction = transmitted.propagation().ray.dir.scale(-1.0);
            let reverse_ray = Ray {
                origin: split.incident.ray.origin.offset(direction.scale(-0.05)),
                dir: direction,
            };
            let reverse = contact
                .split(cx, 1, reverse_ray, AXIS, 1500.0, 0.05)
                .unwrap();
            vector(
                reverse.transmitted.unwrap().propagation().ray.dir,
                launch(0.4, false).dir.scale(-1.0),
            );
            near(
                reverse.fresnel.reflectance,
                split.fresnel.reflectance,
                4e-15,
            );
        },
    );
    fixture(false, [300.0; 2], [1.5; 2], [0.0; 2], |_, cx, contact| {
        let split = contact
            .split(cx, 0, launch(0.4, false), AXIS, 1500.0, 0.1)
            .unwrap();
        assert!(split.reflected.is_none());
        vector(
            split.transmitted.unwrap().propagation().ray.dir,
            launch(0.4, false).dir,
        );
        assert_eq!(split.fresnel.transmittance, 1.0);
    });
}

#[test]
fn g1_total_internal_reflection_keeps_the_original_cell_and_exact_power() {
    fixture(false, [300.0; 2], [1.0, 1.5], [0.0; 2], |_, cx, contact| {
        for (sine, tir) in [(0.6, false), (0.8, true)] {
            let split = contact
                .split(cx, 1, launch(sine, true), AXIS, 1500.0, 0.05)
                .unwrap();
            assert_eq!(split.fresnel.total_internal_reflection, tir);
            let reflected = split.reflected.unwrap();
            assert_eq!(reflected.element_index(), 1);
            assert!(reflected.propagation().ray.dir.z > 0.0);
            if tir {
                assert!(split.transmitted.is_none());
                assert_eq!(reflected.power_fraction(), 1.0);
            } else {
                near(
                    reflected.power_fraction() + split.transmitted.unwrap().power_fraction(),
                    1.0,
                    3e-15,
                );
            }
        }
    });
}

#[test]
fn g3_retained_face_arrivals_compose_without_losing_paths_or_changing_owners() {
    with_cx(|_, cx| {
        // Third cell shares the sloping face opposite local vertex 0 of cell 1.
        let mut positions = pair(false).positions().to_vec();
        positions.extend([
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [2.0, 1.0, 1.0],
        ]);
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(12, vec![[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]]),
            positions,
        )
        .unwrap();
        let values = [300.0; 12];
        let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &values).unwrap();
        let elements = [
            snapshot.element(0).unwrap(),
            snapshot.element(1).unwrap(),
            snapshot.element(2).unwrap(),
        ];
        let glasses = [glass(1.2, 0.2), glass(1.8, 0.3), glass(1.5, 0.4)];
        let law = law();
        let a = MaterialThermoOpticField::try_new(
            &glasses[0],
            &law,
            &elements[0],
            ContentHash([114; 32]),
        )
        .unwrap();
        let b = MaterialThermoOpticField::try_new(
            &glasses[1],
            &law,
            &elements[1],
            ContentHash([114; 32]),
        )
        .unwrap();
        let c = MaterialThermoOpticField::try_new(
            &glasses[2],
            &law,
            &elements[2],
            ContentHash([114; 32]),
        )
        .unwrap();
        let first = P1ContactInterface::try_new(cx, [&a, &b], [3, 3]).unwrap();
        let second = P1ContactInterface::try_new(cx, [&b, &c], [0, 3]).unwrap();
        let result = first
            .split(cx, 0, launch(0.6, false), AXIS, 1500.0, 10.0)
            .unwrap();
        let arrived = result.transmitted.unwrap();
        assert_eq!(arrived.arrival_face(), Some(0));
        let next = second.split_arrival(cx, &arrived, 0.01).unwrap();
        assert_eq!(next.incident, *arrived.propagation());
        let branch = next.transmitted.unwrap();
        near(
            branch.propagation().distance_m,
            arrived.propagation().distance_m + 0.01,
            3e-15,
        );
        near(
            branch.propagation().optical_path_length_m,
            arrived.propagation().optical_path_length_m + 0.015,
            3e-15,
        );
        // Independent s/p powers with rotation between the two incidence
        // planes. Multiplying scalar averages would erase the first contact's
        // polarization and give a different second-contact power.
        let transmission = |ni: f64, nt: f64, ci: f64| {
            let ct = (1.0 - (ni / nt).powi(2) * (1.0 - ci * ci)).sqrt();
            [
                4.0 * ni * nt * ci * ct / (ni * ci + nt * ct).powi(2),
                4.0 * ni * nt * ci * ct / (nt * ci + ni * ct).powi(2),
            ]
        };
        let t1 = transmission(1.2, 1.8, 0.8);
        let normal = Vec3::new(2.0, 1.0, 1.0).scale(1.0 / 6.0_f64.sqrt());
        let direction = arrived.propagation().ray.dir;
        let t2 = transmission(1.8, 1.5, direction.dot(normal));
        let s = cross(direction, normal);
        let c2 = arrived.frame().horizontal().dot(s).powi(2) / s.dot(s);
        let power = 0.5
            * ((t1[0] * c2 + t1[1] * (1.0 - c2)) * t2[0]
                + (t1[0] * (1.0 - c2) + t1[1] * c2) * t2[1])
            * (-arrived.propagation().optical_depth - 0.004).exp();
        near(branch.power_fraction(), power, 4e-15);
        let erased = arrived.power_fraction() * next.fresnel.transmittance * (-0.004_f64).exp();
        assert!((erased - power).abs() > 1e-5);
        assert!(arrived.stokes().degree() > 0.01);
        assert!(matches!(
            first.split_arrival(cx, &arrived, 0.1),
            Err(P1ContactError::ForeignArrival)
        ));
    });
}

#[test]
fn g1_native_tir_phase_and_reference_axis_covariance_survive_contact_transport() {
    fixture(false, [300.0; 2], [1.0, 1.5], [0.0; 2], |_, cx, contact| {
        let ray = launch(0.5_f64.sqrt(), true);
        let input = Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap();
        let a = contact
            .split_polarized(cx, 1, ray, AXIS, input, 1500.0, 0.05)
            .unwrap();
        assert!(a.transmitted.is_none());
        let a = a.reflected.unwrap();
        let out = a.stokes();
        near(out.i, 1.0, 3e-15);
        near(out.q, 0.0, 3e-15);
        near(out.u, 0.8, 3e-14);
        near(out.v, 0.6, 3e-14);
        near(out.degree(), 1.0, 3e-15);
        let b = contact
            .split_polarized(
                cx,
                1,
                ray,
                cross(ray.dir, AXIS),
                Stokes::try_new(1.0, 0.0, -1.0, 0.0).unwrap(),
                1500.0,
                0.05,
            )
            .unwrap()
            .reflected
            .unwrap();
        for (a, b) in [out.i, out.q, out.u, out.v].into_iter().zip([
            b.stokes().i,
            b.stokes().q,
            b.stokes().u,
            b.stokes().v,
        ]) {
            near(a, b, 3e-14);
        }
        vector(a.frame().horizontal(), b.frame().horizontal());
    });
}

#[test]
fn g3_nonuniform_unpaired_foreign_spectral_and_cancelled_inputs_refuse() {
    fixture(
        false,
        [300.0; 2],
        [1.2, 1.5],
        [0.0; 2],
        |gate, cx, contact| {
            assert!(
                contact
                    .split(cx, 0, launch(0.2, false), AXIS, 900.0, 0.1)
                    .is_err()
            );
            assert!(matches!(
                P1ContactInterface::try_new(
                    cx,
                    [contact.sides[0].binding, contact.sides[1].binding],
                    [0, 0]
                ),
                Err(P1ContactError::InvalidPair)
            ));
            let source = contact.sides[0].binding;
            let mesh = source.temperature_source().snapshot().mesh();
            let mut values = source
                .temperature_source()
                .snapshot()
                .nodal_temperature_k()
                .to_vec();
            values[3] = 301.0;
            let snapshot = P1ConductionTemperatureSnapshot::try_new(cx, mesh, &values).unwrap();
            let ea = snapshot.element(0).unwrap();
            let eb = snapshot.element(1).unwrap();
            let a = MaterialThermoOpticField::try_new(
                source.reference_glass(),
                source.law(),
                &ea,
                source.snapshot_identity(),
            )
            .unwrap();
            let b = MaterialThermoOpticField::try_new(
                source.reference_glass(),
                source.law(),
                &eb,
                source.snapshot_identity(),
            )
            .unwrap();
            assert!(matches!(
                P1ContactInterface::try_new(cx, [&a, &b], [3, 3]),
                Err(P1ContactError::NonuniformElement { element: 0 })
            ));
            assert!(matches!(
                P1ContactInterface::try_new(cx, [source, &b], [3, 3]),
                Err(P1ContactError::InvalidPair)
            ));
            let mut outside = launch(0.0, false);
            outside.origin.z = -2.0;
            assert!(matches!(
                contact.split(cx, 0, outside, AXIS, 1500.0, 0.1),
                Err(P1ContactError::InvalidLaunch)
            ));
            gate.request();
            assert!(matches!(
                contact.split(cx, 0, launch(0.2, false), AXIS, 1500.0, 0.1),
                Err(P1ContactError::Mesh(P1MeshError::Propagation(
                    GradientIndexError::Cancelled
                )))
            ));
        },
    );
}
