//! G1/G3 spatial sensor/actuator fixtures with synthetic SI material data.
//! Analytical constant-field solutions test the actual coupled native solve;
//! these fixtures are not experimental validation of a named ceramic.
use fs_blake3::ContentHash;
use fs_couple::piezoelectric::spatial::{
    ElectrodeCondition as Condition, PiezoelectricElectrode as Electrode,
    TetPiezoelectricControls as Controls, TetPiezoelectricError, TetPiezoelectricModel as Model,
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::piezoelectric::tensor::TensorPiezoelectric;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::tensor::{mandel_rotation, rotation};
use fs_solid::linear3::TetAssemblyBudget;

const FRAME: ContentHash = ContentHash([79; 32]);
const LENGTH: f64 = 1e-3;
const STIFFNESS: [f64; 6] = [11e10, 13e10, 10e10, 5e10, 6e10, 7e10];
const EPSILON: [f64; 3] = [2e-8, 1.5e-8, 1e-8];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 79,
                    kernel_id: 79,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}
fn material(e3: [f64; 6]) -> TensorPiezoelectric {
    let mut claims = ClaimSet::new();
    for requirement in TensorPiezoelectric::requirements().unwrap() {
        let name = requirement.name();
        let value = if name == "density" {
            7500.0
        } else {
            let index = name.rsplit('_').next().unwrap().as_bytes();
            let row = (index[0] - b'1') as usize;
            let col = (index[1] - b'1') as usize;
            if name.starts_with("piezoelectric_stiffness_e_mandel_") {
                if row == col { STIFFNESS[row] } else { 0.0 }
            } else if name.starts_with("piezoelectric_stress_charge_mandel_") {
                if row == 2 { e3[col] } else { 0.0 }
            } else {
                assert!(name.starts_with("piezoelectric_permittivity_s_"));
                if row == col { EPSILON[row] } else { 0.0 }
            }
        };
        let dims = requirement.dims();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored spatial piezoelectric SI fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic reciprocal ceramic".into(),
            phase: "solid".into(),
            process: "authored diagonal Mandel stress-charge law".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    TensorPiezoelectric::resolve(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap()
}
struct Cube {
    nodes: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    bottom: Vec<usize>,
    top: Vec<usize>,
    supports: Vec<(usize, f64)>,
}
fn cube(refinement: usize) -> Cube {
    let (mesh, nodes) = fs_feec::fixtures::kuhn_cube(refinement);
    let nodes: Vec<_> = nodes.into_iter().map(|p| p.map(|x| x * LENGTH)).collect();
    let bottom = nodes
        .iter()
        .enumerate()
        .filter_map(|(i, p)| (p[2] == 0.0).then_some(i))
        .collect();
    let top = nodes
        .iter()
        .enumerate()
        .filter_map(|(i, p)| (p[2] == LENGTH).then_some(i))
        .collect();
    let x = refinement * (refinement + 1).pow(2);
    let y = refinement * (refinement + 1);
    Cube {
        nodes,
        tets: mesh.tets.iter().map(|t| t.map(|n| n as usize)).collect(),
        bottom,
        top,
        supports: vec![
            (0, 0.0),
            (1, 0.0),
            (2, 0.0),
            (3 * x + 1, 0.0),
            (3 * x + 2, 0.0),
            (3 * y + 2, 0.0),
        ],
    }
}
fn budget() -> TetAssemblyBudget {
    TetAssemblyBudget {
        maximum_nodes: 64,
        maximum_tetrahedra: 192,
        maximum_free_dofs: 192,
        minimum_scaled_jacobian: 1e-8,
    }
}
fn controls() -> Controls {
    Controls {
        relative_tolerance: 1e-11,
        maximum_iterations: 1024,
        force_tolerance_n: 1e-6,
        charge_tolerance_c: 1e-18,
        energy_tolerance_j: 1e-14,
        maximum_strain_norm: 1e-3,
        maximum_electric_field_v_m: 2e6,
    }
}
fn electrodes(cube: &Cube, reference: f64, top: Condition) -> [Electrode<'_>; 2] {
    [
        Electrode {
            nodes: &cube.bottom,
            condition: Condition::Voltage(reference),
        },
        Electrode {
            nodes: &cube.top,
            condition: top,
        },
    ]
}
fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn apply<const R: usize, const C: usize>(m: &[[f64; C]; R], x: &[f64; C]) -> [f64; R] {
    core::array::from_fn(|i| m[i].iter().zip(x).map(|(a, b)| a * b).sum())
}
fn axial_traction(cube: &Cube, stress: f64) -> Vec<[f64; 3]> {
    // Uniform axial stress needs equal/opposite tractions on both end faces.
    // The six point supports remove rigid modes; they are not a clamped face
    // that could supply the missing distributed bottom traction.
    let mut force = vec![[0.0; 3]; cube.nodes.len()];
    for tet in &cube.tets {
        for omit in 0..4 {
            let face: Vec<_> = (0..4).filter(|&i| i != omit).map(|i| tet[i]).collect();
            let sign = if face.iter().all(|&i| cube.nodes[i][2] == LENGTH) {
                1.0
            } else if face.iter().all(|&i| cube.nodes[i][2] == 0.0) {
                -1.0
            } else {
                continue;
            };
            let a = cube.nodes[face[0]];
            let b = cube.nodes[face[1]];
            let c = cube.nodes[face[2]];
            let area = 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs();
            for node in face {
                force[node][2] += sign * stress * area / 3.0;
            }
        }
    }
    force
}

#[test]
fn g1_free_voltage_actuation_resolves_axial_transverse_and_shear_response() {
    with_cx(|cx, _| {
        let e = [-4.0, -3.0, 12.0, 0.0, 2.0, 5.0];
        let materials = [material(e)];
        let field = 100.0 / LENGTH;
        let strain: [f64; 6] = core::array::from_fn(|i| e[i] * field / STIFFNESS[i]);
        let epsilon_free = EPSILON[2] + (0..6).map(|i| e[i] * e[i] / STIFFNESS[i]).sum::<f64>();
        for refinement in [1, 2] {
            let cube = cube(refinement);
            let n = cube.nodes.len();
            let model =
                Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
            let solved = model
                .solve(
                    &cube.supports,
                    &electrodes(&cube, 0.0, Condition::Voltage(-100.0)),
                    &vec![[0.0; 3]; n],
                    &vec![0.0; n],
                    controls(),
                    cx,
                )
                .unwrap();
            for (point, u) in cube.nodes.iter().zip(&solved.displacement_m) {
                let expected = [
                    strain[0] * point[0] + 2.0_f64.sqrt() * strain[5] * point[2],
                    strain[1] * point[1] + 2.0_f64.sqrt() * strain[4] * point[2],
                    strain[2] * point[2],
                ];
                for axis in 0..3 {
                    near(u[axis], expected[axis], 2e-15, 3e-8);
                }
            }
            for cell in &solved.cells {
                for i in 0..6 {
                    near(cell.strain_mandel[i], strain[i], 2e-12, 3e-8);
                    near(cell.response.stress_mandel_pa[i], 0.0, 1e-2, 0.0);
                }
            }
            let charge = epsilon_free * field * LENGTH * LENGTH;
            near(solved.electrodes[0].charge_c, charge, 2e-18, 3e-8);
            near(solved.electrodes[1].charge_c, -charge, 2e-18, 3e-8);
            near(solved.stored_energy_j, 0.5 * 100.0 * charge, 1e-14, 3e-8);
            near(solved.mechanical_ramp_work_j, 0.0, 1e-14, 0.0);
            near(solved.energy_residual_j, 0.0, 1e-14, 0.0);
            // Nonzero total charge must recover the independent voltage oracle;
            // the open-circuit zero-charge case alone cannot catch a RHS sign bug.
            let charged = model
                .solve(
                    &cube.supports,
                    &electrodes(&cube, 0.0, Condition::Charge(-charge)),
                    &vec![[0.0; 3]; n],
                    &vec![0.0; n],
                    controls(),
                    cx,
                )
                .unwrap();
            near(charged.electrodes[1].voltage_v, -100.0, 2e-8, 3e-8);
            near(charged.stored_energy_j, solved.stored_energy_j, 1e-14, 3e-8);
            // The same prepared operators must solve a translated potential gauge.
            let shifted = model
                .solve(
                    &cube.supports,
                    &electrodes(&cube, 1e9, Condition::Voltage(1e9 - 100.0)),
                    &vec![[0.0; 3]; n],
                    &vec![0.0; n],
                    controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(shifted.displacement_m, solved.displacement_m);
            assert_eq!(shifted.stored_energy_j, solved.stored_energy_j);
            for (a, b) in shifted.voltage_v.iter().zip(&solved.voltage_v) {
                near(a - 1e9, *b, 2e-7, 0.0);
            }
        }
    })
}

#[test]
fn g1_loaded_open_electrodes_stiffen_the_solid_and_generate_voltage() {
    with_cx(|cx, _| {
        let cube = cube(2);
        let n = cube.nodes.len();
        let e = 10.0;
        let stress = 1e6;
        let materials = [material([0.0, 0.0, e, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let force = axial_traction(&cube, stress);
        let charge = vec![0.0; n];
        let short = model
            .solve(
                &cube.supports,
                &electrodes(&cube, 0.0, Condition::Voltage(0.0)),
                &force,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let open = model
            .solve(
                &cube.supports,
                &electrodes(&cube, 0.0, Condition::Charge(0.0)),
                &force,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let s_short = stress / STIFFNESS[2];
        let s_open = stress / (STIFFNESS[2] + e * e / EPSILON[2]);
        for &node in &cube.top {
            near(short.displacement_m[node][2], s_short * LENGTH, 2e-15, 3e-8);
            near(open.displacement_m[node][2], s_open * LENGTH, 2e-15, 3e-8);
            assert!(open.displacement_m[node][2] < short.displacement_m[node][2]);
        }
        near(
            open.electrodes[1].voltage_v,
            e * s_open * LENGTH / EPSILON[2],
            2e-8,
            3e-8,
        );
        near(open.electrodes[1].charge_c, 0.0, 1e-18, 0.0);
        near(
            short.electrodes[1].charge_c,
            -e * s_short * LENGTH * LENGTH,
            1e-18,
            3e-8,
        );
        for (solved, strain) in [(&short, s_short), (&open, s_open)] {
            near(
                solved.stored_energy_j,
                0.5 * stress * strain * LENGTH.powi(3),
                1e-14,
                3e-8,
            );
            near(solved.electrical_ramp_work_j, 0.0, 1e-14, 0.0);
            near(
                solved.mechanical_ramp_work_j,
                solved.stored_energy_j,
                1e-14,
                3e-8,
            );
        }
    })
}

#[test]
fn g3_rotating_the_actual_mesh_and_tensor_preserves_field_and_work() {
    with_cx(|cx, _| {
        let cube = cube(2);
        let n = cube.nodes.len();
        let q = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.63);
        let m = mandel_rotation(q);
        let frame = ContentHash([80; 32]);
        let materials = [material([-4.0, -3.0, 12.0, 0.0, 2.0, 5.0])];
        let rotated_materials = [materials[0].rotated(frame, q).unwrap()];
        let rotated_nodes: Vec<_> = cube.nodes.iter().map(|x| apply(&q, x)).collect();
        let strain = [2e-6, -1e-6, 3e-6, 0.0, 4e-7, -2e-7];
        let f = core::f64::consts::FRAC_1_SQRT_2;
        let affine = [
            [strain[0], strain[3] * f, strain[5] * f],
            [strain[3] * f, strain[1], strain[4] * f],
            [strain[5] * f, strain[4] * f, strain[2]],
        ];
        let u: Vec<_> = cube.nodes.iter().map(|x| apply(&affine, x)).collect();
        let prescribed: Vec<_> = u
            .iter()
            .enumerate()
            .flat_map(|(i, u)| u.iter().enumerate().map(move |(j, &v)| (3 * i + j, v)))
            .collect();
        let rotated_u: Vec<_> = u.iter().map(|u| apply(&q, u)).collect();
        let rotated_prescribed: Vec<_> = rotated_u
            .iter()
            .enumerate()
            .flat_map(|(i, u)| u.iter().enumerate().map(move |(j, &v)| (3 * i + j, v)))
            .collect();
        let electrodes = electrodes(&cube, 0.0, Condition::Voltage(-30.0));
        let force = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let base = Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx)
            .unwrap()
            .solve(&prescribed, &electrodes, &force, &charge, controls(), cx)
            .unwrap();
        let rotated = Model::assemble(
            &rotated_nodes,
            &cube.tets,
            &rotated_materials,
            frame,
            budget(),
            cx,
        )
        .unwrap()
        .solve(
            &rotated_prescribed,
            &electrodes,
            &force,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        for (a, b) in base.cells.iter().zip(&rotated.cells) {
            let stress = apply(&m, &a.response.stress_mandel_pa);
            let displacement = apply(&q, &a.response.displacement_c_m2);
            let field = apply(&q, &a.response.electric_field_v_m);
            for i in 0..6 {
                near(b.response.stress_mandel_pa[i], stress[i], 1e-4, 3e-8);
            }
            for i in 0..3 {
                near(
                    b.response.displacement_c_m2[i],
                    displacement[i],
                    1e-16,
                    3e-8,
                );
                near(b.response.electric_field_v_m[i], field[i], 1e-5, 3e-8);
            }
        }
        near(rotated.stored_energy_j, base.stored_energy_j, 1e-14, 3e-8);
        near(
            rotated.mechanical_ramp_work_j,
            base.mechanical_ramp_work_j,
            1e-14,
            3e-8,
        );
        near(
            rotated.electrical_ramp_work_j,
            base.electrical_ramp_work_j,
            1e-14,
            3e-8,
        );
    })
}

#[test]
fn g0_g4_boundary_frame_envelope_and_cancellation_refuse_without_publication() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 10.0, 0.0, 0.0, 0.0])];
        assert!(
            Model::assemble(
                &cube.nodes,
                &cube.tets,
                &materials,
                ContentHash([80; 32]),
                budget(),
                cx
            )
            .is_err()
        );
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let electrodes = electrodes(&cube, 0.0, Condition::Voltage(-100.0));
        let force = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let valid = model
            .solve(&cube.supports, &electrodes, &force, &charge, controls(), cx)
            .unwrap();
        assert!(
            model
                .solve(&[], &electrodes, &force, &charge, controls(), cx)
                .is_err()
        );
        let mut duplicate = cube.supports.clone();
        duplicate.push(cube.supports[0]);
        assert!(
            model
                .solve(&duplicate, &electrodes, &force, &charge, controls(), cx)
                .is_err()
        );
        let overlap = [
            electrodes[0],
            Electrode {
                nodes: &cube.bottom,
                condition: Condition::Voltage(1.0),
            },
        ];
        assert!(
            model
                .solve(&cube.supports, &overlap, &force, &charge, controls(), cx)
                .is_err()
        );
        let point = [
            Electrode {
                nodes: &cube.bottom[..1],
                condition: Condition::Voltage(0.0),
            },
            electrodes[1],
        ];
        assert!(
            model
                .solve(&cube.supports, &point, &force, &charge, controls(), cx)
                .is_err()
        );
        let floating = [
            Electrode {
                nodes: &cube.bottom,
                condition: Condition::Charge(0.0),
            },
            Electrode {
                nodes: &cube.top,
                condition: Condition::Charge(0.0),
            },
        ];
        assert!(
            model
                .solve(&cube.supports, &floating, &force, &charge, controls(), cx)
                .is_err()
        );
        let mut narrow = controls();
        narrow.maximum_electric_field_v_m = 100.0;
        assert!(
            model
                .solve(&cube.supports, &electrodes, &force, &charge, narrow, cx)
                .is_err()
        );
        assert_eq!(
            valid.displacement_m,
            model
                .solve(&cube.supports, &electrodes, &force, &charge, controls(), cx)
                .unwrap()
                .displacement_m
        );
        gate.request();
        assert!(matches!(
            model.solve(&cube.supports, &electrodes, &force, &charge, controls(), cx),
            Err(TetPiezoelectricError::Cancelled)
        ));
        assert!(matches!(
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx),
            Err(TetPiezoelectricError::Cancelled)
        ));
    })
}
