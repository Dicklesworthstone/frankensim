//! G0/G1/G3: material-bound thermal strain drives the existing quadratic
//! elasticity operators and solver, with no thermal-only stiffness model.
#[allow(dead_code, unused_imports)]
#[path = "hex/support.rs"]
mod hex_support;
#[allow(dead_code)]
#[path = "tet10/support.rs"]
mod tet_support;

use fs_blake3::ContentHash;
use fs_solid::linear3::tet10::QuadraticTetProblem;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetElasticMaterial, TetMaterialField,
    TetThermalLoad, TetThermalStrainField, TetThermalStrainState,
};
use fs_sparse::{DirectOrdering, LdltOptions, SymbolicLdlt};
use tet_support::{block, material, promote, with_cx};

fn thermal(material: &TetElasticMaterial, strain: [f64; 6]) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        400.0,
        300.0,
        strain,
        material.material_state_identity,
        ContentHash([77; 32]),
        0.01,
    )
    .unwrap()
}

fn check_affine_patch(
    nodes: &[[f64; 3]],
    assembly: &TetElasticAssembly,
    load: &TetThermalLoad,
    strain: [f64; 6],
    material: &TetElasticMaterial,
) {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let displacement: Vec<_> = nodes
        .iter()
        .flat_map(|&[x, y, z]| {
            [
                strain[0] * x + s * strain[3] * y + s * strain[5] * z,
                s * strain[3] * x + strain[1] * y + s * strain[4] * z,
                s * strain[5] * x + s * strain[4] * y + strain[2] * z,
            ]
        })
        .collect();
    let mut applied = vec![0.0; displacement.len()];
    assembly.stiffness.spmv(&displacement, &mut applied);
    let scale = load
        .full_equivalent_force_n
        .iter()
        .map(|v| v.abs())
        .fold(0.0f64, f64::max);
    for (elastic, thermal) in applied.iter().zip(&load.full_equivalent_force_n) {
        assert!(
            (elastic - thermal).abs() < 1e-11 * scale,
            "K u = {elastic}, thermal force = {thermal}"
        );
    }
    let c = material.stiffness_mandel_pa();
    let density: f64 = (0..6)
        .map(|i| (0..6).map(|j| strain[i] * c[i][j] * strain[j]).sum::<f64>())
        .sum();
    let exact_work = assembly.element_volumes_m3.iter().sum::<f64>() * density;
    let work: f64 = displacement
        .iter()
        .zip(&load.full_equivalent_force_n)
        .map(|(u, f)| u * f)
        .sum();
    assert!((work - exact_work).abs() < 1e-12 * exact_work);
    // Internal eigenstrain cannot supply a net force or torque to the body.
    let mut resultant = [0.0; 3];
    let mut moment = [0.0; 3];
    for (p, f) in nodes
        .iter()
        .zip(load.full_equivalent_force_n.chunks_exact(3))
    {
        for i in 0..3 {
            resultant[i] += f[i];
        }
        moment[0] += p[1] * f[2] - p[2] * f[1];
        moment[1] += p[2] * f[0] - p[0] * f[2];
        moment[2] += p[0] * f[1] - p[1] * f[0];
    }
    assert!(
        resultant
            .iter()
            .chain(&moment)
            .all(|v| v.abs() < 1e-11 * scale)
    );
}

#[test]
fn g1_quadratic_thermal_patch_includes_shear_and_anisotropic_elasticity() {
    with_cx(|cx| {
        // Deliberately unequal shear/normal moduli catch Mandel-order mistakes.
        let mut c = [[0.0; 6]; 6];
        for (i, row) in c.iter_mut().enumerate() {
            row[i] = (i + 2) as f64 * 1e6;
        }
        c[0][1] = 0.2e6;
        c[1][0] = c[0][1];
        let law = TetElasticMaterial::try_new_mandel(1200.0, c, [58; 32]).unwrap();
        let strain = [0.001, -0.0004, 0.0002, 0.0003, -0.0001, 0.0005];
        let state = thermal(&law, strain);
        let (mut nodes, elements) = promote(&block([1, 1, 1], [1.0; 3]), cx);
        // Isoparametric geometry still represents an affine physical u exactly.
        for p in &mut nodes {
            *p = hex_support::distortion(*p);
        }
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        let load = problem
            .assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx)
            .unwrap();
        check_affine_patch(&nodes, &problem.assemble(cx).unwrap(), &load, strain, &law);
        assert_eq!(
            load,
            problem
                .assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx)
                .unwrap()
        );
        for family in hex_support::FAMILIES {
            let mut mesh = hex_support::mesh([1, 1, 1], [1.0; 3], family, cx);
            for p in &mut mesh.nodes_m {
                *p = hex_support::distortion(*p);
            }
            let problem = hex_support::problem(&mesh, &law, &[]);
            let load = problem
                .assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx)
                .unwrap();
            check_affine_patch(
                &mesh.nodes_m,
                &problem.assemble(cx).unwrap(),
                &load,
                strain,
                &law,
            );
            assert_eq!(
                load,
                problem
                    .assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx)
                    .unwrap()
            );
        }
    });
}

fn symmetry_supports(nodes: &[[f64; 3]]) -> Vec<usize> {
    nodes
        .iter()
        .enumerate()
        .flat_map(|(node, p)| {
            (0..3)
                .filter(move |&axis| p[axis] == 0.0)
                .map(move |axis| 3 * node + axis)
        })
        .collect()
}

fn check_free_expansion(nodes: &[[f64; 3]], assembly: TetElasticAssembly, load: TetThermalLoad) {
    assert_eq!(assembly.free_dofs, load.free_dofs);
    let symbolic = SymbolicLdlt::analyze(&assembly.stiffness, DirectOrdering::Amd).unwrap();
    let factor = symbolic
        .factor(&assembly.stiffness, &LdltOptions::default())
        .unwrap();
    let u = factor.solve(&load.reduced_force_n);
    for (&dof, value) in load.free_dofs.iter().zip(u) {
        let expected = 0.001 * nodes[dof / 3][dof % 3];
        assert!(
            (value - expected).abs() < 1e-12,
            "dof {dof}: {value} != {expected}"
        );
    }
}

#[test]
fn g1_existing_solver_recovers_thermal_expansion_on_all_quadratic_families() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let state = thermal(&law, [0.001, 0.001, 0.001, 0.0, 0.0, 0.0]);
        let (nodes, elements) = promote(&block([1, 1, 1], [1.0; 3]), cx);
        let fixed = symmetry_supports(&nodes);
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        check_free_expansion(
            &nodes,
            problem.assemble(cx).unwrap(),
            problem
                .assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx)
                .unwrap(),
        );
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1, 1, 1], [1.0; 3], family, cx);
            let fixed = symmetry_supports(&mesh.nodes_m);
            let problem = hex_support::problem(&mesh, &law, &fixed);
            check_free_expansion(
                &mesh.nodes_m,
                problem.assemble(cx).unwrap(),
                problem
                    .assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx)
                    .unwrap(),
            );
        }
    });
}

#[test]
fn g0_source_mismatch_shape_budget_and_constrained_loads_are_preserved() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let state = thermal(&law, [0.001, 0.001, 0.001, 0.0, 0.0, 0.0]);
        let foreign = thermal(
            &TetElasticMaterial::try_new(1.0, 2e6, 0.25, [99; 32]).unwrap(),
            [0.001, 0.001, 0.001, 0.0, 0.0, 0.0],
        );
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1, 1, 1], [1.0; 3], family, cx);
            let fixed: Vec<_> = (0..3 * mesh.nodes_m.len()).collect();
            let mut problem = hex_support::problem(&mesh, &law, &fixed);
            let load = problem
                .assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx)
                .unwrap();
            assert!(load.reduced_force_n.is_empty());
            assert!(load.free_dofs.is_empty());
            // Fully restrained expansion has sigma_xx=-E*epsilon/(1-2 nu).
            let face_force: f64 = mesh
                .nodes_m
                .iter()
                .enumerate()
                .filter(|(_, p)| p[0] == 1.0)
                .map(|(i, _)| load.full_equivalent_force_n[3 * i])
                .sum();
            assert!((face_force - 4000.0).abs() < 1e-9);
            assert!(matches!(
                problem.assemble_thermal_load(TetThermalStrainField::Uniform(&foreign), cx),
                Err(TetElasticError::ThermalElasticStateMismatch { element: 0 })
            ));
            assert!(matches!(
                problem.assemble_thermal_load(TetThermalStrainField::PerElement(&[]), cx),
                Err(TetElasticError::ThermalStrainCountMismatch { .. })
            ));
            problem.budget.maximum_nodes = mesh.nodes_m.len() - 1;
            assert!(matches!(
                problem.assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx),
                Err(TetElasticError::BudgetExceeded { .. })
            ));
        }
        let (nodes, elements) = promote(&block([1, 1, 1], [1.0; 3]), cx);
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        assert!(matches!(
            problem.assemble_thermal_load(TetThermalStrainField::Uniform(&foreign), cx),
            Err(TetElasticError::ThermalElasticStateMismatch { element: 0 })
        ));
        assert!(matches!(
            problem.assemble_thermal_load(TetThermalStrainField::PerElement(&[]), cx),
            Err(TetElasticError::ThermalStrainCountMismatch { .. })
        ));
    });
}

#[test]
fn g4_cancelled_request_never_returns_quadratic_thermal_forces() {
    let mesh =
        with_cx(|cx| hex_support::mesh([1; 3], [1.0; 3], hex_support::HexFamily::Lagrange27, cx));
    let (nodes, elements) = with_cx(|cx| promote(&block([1; 3], [1.0; 3]), cx));
    let law = material(2e6, 0.25);
    let state = thermal(&law, [0.001, 0.001, 0.001, 0.0, 0.0, 0.0]);
    let gate = fs_exec::CancelGate::new_clock_free();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 58,
                kernel_id: 58,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        assert!(matches!(
            problem.assemble_thermal_load(TetThermalStrainField::Uniform(&state), &cx),
            Err(TetElasticError::Cancelled)
        ));
        assert!(matches!(
            hex_support::problem(&mesh, &law, &[])
                .assemble_thermal_load(TetThermalStrainField::Uniform(&state), &cx),
            Err(TetElasticError::Cancelled)
        ));
    });
}
