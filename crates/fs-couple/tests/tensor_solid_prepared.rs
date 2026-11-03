//! Original quadrature, consistent inertia, exact material memory and nodal work.
//! Synthetic source fixtures; no material calibration or mesh convergence claim.
use fs_couple::tensor_material::{
    InitialMemory, TensorMaterialBody, TensorMaterialBudget, prepared::PreparedTensorControls,
    solid::prepared::PreparedTensorSolidBody,
};
use fs_exec::{CancelGate, Cx};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_solid::linear3::{
    TetAssemblyBudget,
    tet10::{TetQuadrature, promote_p1_to_p2},
    visco::{
        SolidMeshRef, SolidProjectionBudget, TetViscoField, TetViscoMaterial, ViscoElements,
        ViscoStateBudget,
    },
};

const CORNERS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 0.2, 0.0],
    [0.0, 0.0, 0.2],
];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 62,
                    kernel_id: 44,
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

fn controls() -> PreparedTensorControls {
    PreparedTensorControls {
        maximum_coordinates: 30,
        maximum_material_points: 64,
        maximum_state_coordinates: 120,
        maximum_heat_groups: 2,
        maximum_preparation_bytes: 512 * 1024 * 1024,
        maximum_preparation_multiply_terms: 1_000_000_000_000,
        maximum_step_multiply_terms: 100_000_000,
        maximum_relative_operator_error: 1e-9,
        maximum_work_error_j: 1e-10,
        maximum_steps: 1024,
    }
}

fn material(kelvin_time: f64, fraction: f64) -> TetViscoMaterial {
    let diagonal = [120e6, 90e6, 120e6, 60e6, 50e6, 40e6];
    let equilibrium: Mandel6 =
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { diagonal[i] } else { 0.0 }));
    let branches = if fraction == 0.0 {
        vec![]
    } else {
        vec![PronyBranch::new(equilibrium.map(|r| r.map(|x| fraction * x)), 0.004).unwrap()]
    };
    TetViscoMaterial {
        law: TensorProny::new(
            equilibrium,
            if kelvin_time == 0.0 {
                None
            } else {
                Some(equilibrium.map(|r| r.map(|x| kelvin_time * x)))
            },
            branches,
            (0.0, 1e6),
        )
        .unwrap(),
        density_kg_m3: 1000.0,
        material_state_identity: [44; 32],
    }
}

fn projection_budget() -> SolidProjectionBudget {
    SolidProjectionBudget {
        max_coordinates: 30,
        max_basis_scalars: 900,
        max_material_points: 64,
        max_scalar_terms: 20_000_000,
    }
}

fn identity(n: usize) -> Vec<Vec<f64>> {
    (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect()
}

fn quadratic(matrix: &fs_sparse::Csr, values: &[f64]) -> f64 {
    let mut product = vec![0.0; values.len()];
    matrix.spmv(values, &mut product);
    0.5 * values.iter().zip(product).map(|(a, b)| a * b).sum::<f64>()
}

#[test]
fn g1_all_thirty_tet10_coordinates_retain_original_memory_heat_and_work() {
    with_cx(|cx, gate| {
        let (nodes, cells) = promote_p1_to_p2(&CORNERS, &[[0, 1, 2, 3]], None, cx).unwrap();
        let material = material(0.0001, 0.4);
        let source = ViscoElements::new(
            SolidMeshRef::Tet10 {
                nodes_m: &nodes,
                elements: &cells,
                quadrature: TetQuadrature {
                    stiffness_points_per_dir: 3,
                    mass_points_per_dir: 4,
                },
            },
            TetViscoField::Uniform(&material),
            &[],
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 1000,
            },
            cx,
        )
        .unwrap();
        let projection = source
            .project_basis(&identity(30), projection_budget(), cx)
            .unwrap();
        assert_eq!(projection.material_points().len(), 27);
        assert!(projection.mass_matrix_kg().iter().any(|m| *m < 0.0));
        // The unchanged owner retains 60 + 27*6 = 222 physical states and
        // requires a 506-coordinate exponential block. It must still refuse.
        assert!(
            TensorMaterialBody::from_mass_matrix(
                30,
                projection.mass_matrix_kg().to_vec(),
                projection.material_points().to_vec(),
                TensorMaterialBudget {
                    max_state_coordinates: 256,
                    max_material_points: 64
                },
                gate,
            )
            .is_err()
        );
        let dt = 1.0 / 48000.0;
        let body = PreparedTensorSolidBody::new(projection, dt, controls(), gate).unwrap();
        assert_eq!(body.material_body().coordinates(), 30);
        assert_eq!(body.material_body().state_dimension(), 90);
        assert_eq!(body.material_body().relaxation_times_s(), &[0.004]);
        let strain = 1e-5;
        let initial: Vec<_> = nodes
            .iter()
            .flat_map(|p| [strain * p[0], 0.0, 0.0])
            .collect();
        let mut history = body
            .history(&initial, &[0.0; 30], InitialMemory::Instantaneous, gate)
            .unwrap();
        let expected_energy = 0.5 * (0.2 * 0.2 / 6.0) * 120e6 * 1.4 * strain * strain;
        assert!((history.material_history().stored_energy_j() - expected_energy).abs() < 2e-16);
        let mut independent_work = 0.0;
        let mut heat = 0.0;
        let mut maximum_balance = 0.0_f64;
        for step in 0..480 {
            let before = history.nodal_displacement_m(gate).unwrap();
            let mut force = [[0.0; 3]; 10];
            if step < 24 {
                force[0][0] = -0.1;
                force[1][0] = 0.1;
            }
            let report = history.advance(&force, gate).unwrap();
            let after = history.nodal_displacement_m(gate).unwrap();
            independent_work += force
                .iter()
                .zip(before.iter().zip(&after))
                .map(|(f, (a, b))| {
                    (0..3)
                        .map(|axis| f[axis] * (b[axis] - a[axis]))
                        .sum::<f64>()
                })
                .sum::<f64>();
            assert_eq!(report.element_heat_j().len(), 1);
            assert!(report.element_heat_j()[0] >= 0.0);
            heat += report.element_heat_j()[0];
            let physical = history.material_history();
            let equilibrium = quadratic(&source.equilibrium().stiffness, physical.displacement_m());
            let kinetic = quadratic(&source.equilibrium().mass, physical.velocity_mps());
            let memory: f64 = body
                .projection()
                .material_points()
                .iter()
                .flat_map(|p| p.branches())
                .flat_map(|b| b.rows())
                .map(|row| {
                    let z = row
                        .iter()
                        .zip(physical.memory_displacement_m())
                        .map(|(a, y)| a * y)
                        .sum::<f64>();
                    0.5 * z * z
                })
                .sum();
            assert!((report.material.energy.equilibrium_after_j - equilibrium).abs() < 2e-15);
            assert!((report.material.energy.kinetic_after_j - kinetic).abs() < 2e-15);
            assert!((report.material.energy.memory_after_j - memory).abs() < 2e-15);
            let balance =
                equilibrium + kinetic + memory + heat - expected_energy - independent_work;
            maximum_balance = maximum_balance.max(balance.abs());
            assert!(
                balance.abs() < 1e-10,
                "original physical balance {balance:e}"
            );
            assert!(report.nodal_work_error_j.abs() < 1e-14);
        }
        assert!(heat > 0.5 * expected_energy);
        println!(
            "complete TET10: points=27, coordinates=30, original_states=222, shared_states=90, maximum_original_balance_j={maximum_balance:e}, heat_j={heat:e}"
        );
    });
}

#[test]
fn g1_original_elements_receive_their_own_kelvin_and_maxwell_heat() {
    with_cx(|cx, gate| {
        let mut nodes = CORNERS.to_vec();
        nodes.extend(CORNERS.map(|mut p| {
            p[0] += 2.0;
            p
        }));
        let cells = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let fixed: Vec<_> = (0..24).filter(|&i| i != 11 && i != 23).collect();
        let original = [material(0.0002, 0.0), material(0.0, 0.5)];
        let swapped = [original[1].clone(), original[0].clone()];
        let make_source = |materials| {
            ViscoElements::new(
                SolidMeshRef::Tet4 {
                    nodes_m: &nodes,
                    elements: &cells,
                },
                TetViscoField::PerElement(materials),
                &fixed,
                TetAssemblyBudget::standard(),
                ViscoStateBudget {
                    max_internal_variables: 1000,
                },
                cx,
            )
            .unwrap()
        };
        let a = make_source(&original);
        let b = make_source(&swapped);
        let ab = PreparedTensorSolidBody::new(
            a.project_basis(&identity(2), projection_budget(), cx)
                .unwrap(),
            1e-5,
            controls(),
            gate,
        )
        .unwrap();
        let bb = PreparedTensorSolidBody::new(
            b.project_basis(&identity(2), projection_budget(), cx)
                .unwrap(),
            1e-5,
            controls(),
            gate,
        )
        .unwrap();
        let mut ah = ab
            .history(
                &[1e-6, 1e-6],
                &[0.01, 0.01],
                InitialMemory::Instantaneous,
                gate,
            )
            .unwrap();
        let mut bh = bb
            .history(
                &[1e-6, 1e-6],
                &[0.01, 0.01],
                InitialMemory::Instantaneous,
                gate,
            )
            .unwrap();
        let mut force = [[0.0; 3]; 8];
        force[0] = [100.0; 3]; // Fixed rows do no work.
        let ar = ah.advance(&force, gate).unwrap();
        let br = bh.advance(&force, gate).unwrap();
        assert!(ar.element_heat_j().iter().all(|h| *h > 0.0));
        assert!((ar.element_heat_j()[0] - ar.element_heat_j()[1]).abs() > 1e-10);
        assert!((ar.element_heat_j()[0] - br.element_heat_j()[1]).abs() < 1e-12);
        assert!((ar.element_heat_j()[1] - br.element_heat_j()[0]).abs() < 1e-12);
        assert_eq!(ar.material.energy.supplied_work_j, 0.0);
        let before = ah.material_history().displacement_m().to_vec();
        let time = ah.material_history().time_s();
        let heat = ah.material_history().dissipated_heat_j();
        assert!(ah.advance(&force[..7], gate).is_err());
        assert_eq!(ah.material_history().displacement_m(), before.as_slice());
        assert_eq!(ah.material_history().time_s(), time);
        assert_eq!(ah.material_history().dissipated_heat_j(), heat);
    });
}
