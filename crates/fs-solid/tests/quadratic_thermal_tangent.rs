//! Quadratic thermal equilibrium reaches the existing signed modal pencil.
//! Analytic contractions distinguish physical stress from thermal load or means.
#[allow(dead_code, unused_imports)]
#[path = "hex/support.rs"]
mod hex_support;
#[allow(dead_code)]
#[path = "tet10/support.rs"]
mod tet_support;

use fs_blake3::ContentHash;
use fs_modal::eigh_gen_dense;
use fs_solid::linear3::tet10::QuadraticTetProblem;
use fs_solid::linear3::{
    QuadraticPrestressedAssembly, TetAssemblyBudget, TetElasticAssembly, TetElasticError,
    TetElasticMaterial, TetMaterialField, TetStaticSolveConfig, TetThermalStrainField,
    TetThermalStrainState,
};
use tet_support::{block, material, promote, quadratic_form, with_cx};

const FRAME: ContentHash = ContentHash([91; 32]);
const YOUNG: f64 = 2e6;
const STRAIN: f64 = 0.001;

fn thermal(law: &TetElasticMaterial, axial: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        300.0 + 10000.0 * axial,
        300.0,
        [axial, 0.0, 0.0, 0.0, 0.0, 0.0],
        law.material_state_identity,
        ContentHash([92; 32]),
        0.02,
    )
    .unwrap()
}

fn supports(nodes: &[[f64; 3]], both_ends: bool) -> Vec<usize> {
    nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 0.0 || both_ends && p[0] == 1.0)
        .flat_map(|(node, _)| [3 * node, 3 * node + 1, 3 * node + 2])
        .collect()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn eigenvalues(assembly: &TetElasticAssembly) -> Vec<f64> {
    eigh_gen_dense(
        &assembly.stiffness.to_dense(),
        &assembly.mass.to_dense(),
        assembly.free_dofs.len(),
    )
    .unwrap()
    .into_iter()
    .map(|m| m.lambda)
    .collect()
}

fn check_pencils(
    nodes: &[[f64; 3]],
    original: &TetElasticAssembly,
    pencils: &[QuadraticPrestressedAssembly; 3],
) {
    // With nu=0 and both x faces fixed, u=0 and sigma_xx=-E*epsilon_free.
    // All three families represent the admissible v_y=x(1-x) exactly.
    // Thus v^T K_g v = integral sigma_xx*(1-2x)^2 dV = sigma_xx/3.
    let probe: Vec<_> = original
        .free_dofs
        .iter()
        .map(|&dof| {
            let x = nodes[dof / 3][0];
            if dof % 3 == 1 { x * (1.0 - x) } else { 0.0 }
        })
        .collect();
    assert!(probe.iter().any(|&v| v != 0.0));
    for (axial, assembly) in [-STRAIN, 0.0, STRAIN].into_iter().zip(pencils) {
        let pencil = assembly.pencil();
        assert_eq!(pencil.free_dofs, original.free_dofs);
        assert_eq!(pencil.mass.to_dense(), original.mass.to_dense());
        assert_eq!(pencil.total_mass_kg, original.total_mass_kg);
        assert_eq!(pencil.element_volumes_m3, original.element_volumes_m3);
        assert_eq!(
            assembly.elastic_stiffness().to_dense(),
            original.stiffness.to_dense()
        );
        assert!(
            assembly
                .response()
                .displacement_m()
                .iter()
                .flatten()
                .all(|u| u.abs() < 1e-11)
        );
        for element in assembly.response().elements() {
            close(element.mean_stress_mandel_pa()[0], -YOUNG * axial);
            for &stress in &element.mean_stress_mandel_pa()[1..] {
                close(stress, 0.0);
            }
        }
        close(
            quadratic_form(assembly.geometric_stiffness(), &probe),
            -YOUNG * axial / 3.0,
        );
        close(
            quadratic_form(&pencil.stiffness, &probe),
            quadratic_form(&original.stiffness, &probe) - YOUNG * axial / 3.0,
        );
        for row in 0..pencil.free_dofs.len() {
            for col in 0..pencil.free_dofs.len() {
                assert_eq!(
                    assembly.geometric_stiffness().get(row, col).to_bits(),
                    assembly.geometric_stiffness().get(col, row).to_bits()
                );
            }
        }
    }
    let [cold, zero, hot] = pencils;
    assert!(
        zero.geometric_stiffness()
            .to_dense()
            .iter()
            .all(|&v| v == 0.0)
    );
    assert_eq!(
        zero.pencil().stiffness.to_dense(),
        original.stiffness.to_dense()
    );
    assert_ne!(cold.identity(), hot.identity());
    let (tension, unloaded, compression) = (
        eigenvalues(cold.pencil()),
        eigenvalues(zero.pencil()),
        eigenvalues(hot.pencil()),
    );
    assert_eq!(tension.len(), original.free_dofs.len());
    for ((&t, &u), &c) in tension.iter().zip(&unloaded).zip(&compression) {
        let roundoff = 1e-9 * u.abs().max(1.0);
        assert!(
            c <= u + roundoff && u <= t + roundoff,
            "compression {c:e}, unloaded {u:e}, tension {t:e}"
        );
    }
    assert!(
        compression
            .iter()
            .zip(&unloaded)
            .any(|(&c, &u)| c < u * (1.0 - 1e-7))
    );
    assert!(
        tension
            .iter()
            .zip(&unloaded)
            .any(|(&t, &u)| t > u * (1.0 + 1e-7))
    );
}

#[test]
fn g1_all_quadratic_families_retain_mass_and_receive_signed_thermal_mode_shifts() {
    with_cx(|cx| {
        let law = material(YOUNG, 0.0);
        let states = [-STRAIN, 0.0, STRAIN].map(|a| thermal(&law, a));
        let (nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
        let fixed = supports(&nodes, true);
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        let pencils = states.each_ref().map(|state| {
            problem
                .solve_thermal_tangent(
                    TetThermalStrainField::Uniform(state),
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap()
        });
        check_pencils(&nodes, &problem.assemble(cx).unwrap(), &pencils);
        assert_eq!(
            pencils[1].identity(),
            problem
                .solve_thermal_tangent(
                    TetThermalStrainField::Uniform(&states[1]),
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap()
                .identity()
        );
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let fixed = supports(&mesh.nodes_m, true);
            let problem = hex_support::problem(&mesh, &law, &fixed);
            let pencils = states.each_ref().map(|state| {
                problem
                    .solve_thermal_tangent(
                        TetThermalStrainField::Uniform(state),
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.02,
                        cx,
                    )
                    .unwrap()
            });
            check_pencils(&mesh.nodes_m, &problem.assemble(cx).unwrap(), &pencils);
        }
    });
}

// Use the original elasticity operator solely to specify an exactly representable
// manufactured dead load. The stress and geometric energy references below are
// independent continuum formulas, not another copy of the K_g contraction.
fn manufactured_forces(nodes: &[[f64; 3]], original: &TetElasticAssembly) -> Vec<(usize, f64)> {
    let displacement: Vec<_> = original
        .free_dofs
        .iter()
        .map(|&dof| {
            if dof % 3 == 0 {
                STRAIN * nodes[dof / 3][0].powi(2)
            } else {
                0.0
            }
        })
        .collect();
    let mut force = vec![0.0; displacement.len()];
    original.stiffness.spmv(&displacement, &mut force);
    original.free_dofs.iter().copied().zip(force).collect()
}

fn check_variable_stress(nodes: &[[f64; 3]], tangent: &QuadraticPrestressedAssembly) {
    // u_x=b*x^2, nu=0 => sigma_xx=2*E*b*x.
    // For v_y=x^2, integral sigma_xx*(2x)^2 dV = 2*E*b.
    // One HEX mean-stress substitution instead gives 4*E*b/3.
    for (p, u) in nodes.iter().zip(tangent.response().displacement_m()) {
        close(u[0], STRAIN * p[0].powi(2));
        close(u[1], 0.0);
        close(u[2], 0.0);
    }
    close(
        tangent.response().elastic_energy_j(),
        2.0 * YOUNG * STRAIN.powi(2) / 3.0,
    );
    let probe: Vec<_> = tangent
        .pencil()
        .free_dofs
        .iter()
        .map(|&dof| {
            if dof % 3 == 1 {
                nodes[dof / 3][0].powi(2)
            } else {
                0.0
            }
        })
        .collect();
    close(
        quadratic_form(tangent.geometric_stiffness(), &probe),
        2.0 * YOUNG * STRAIN,
    );
}

#[test]
fn g1_quadratic_mixed_load_uses_actual_spatial_stress_instead_of_element_means() {
    with_cx(|cx| {
        let law = material(YOUNG, 0.0);
        let state = thermal(&law, 0.0);
        let field = TetThermalStrainField::Uniform(&state);
        let (nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
        let fixed = supports(&nodes, false);
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        let original = problem.assemble(cx).unwrap();
        let forces = manufactured_forces(&nodes, &original);
        check_variable_stress(
            &nodes,
            &problem
                .solve_thermal_tangent_with_nodal_forces(
                    field,
                    &forces,
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap(),
        );
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let fixed = supports(&mesh.nodes_m, false);
            let problem = hex_support::problem(&mesh, &law, &fixed);
            let original = problem.assemble(cx).unwrap();
            let forces = manufactured_forces(&mesh.nodes_m, &original);
            check_variable_stress(
                &mesh.nodes_m,
                &problem
                    .solve_thermal_tangent_with_nodal_forces(
                        field,
                        &forces,
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.02,
                        cx,
                    )
                    .unwrap(),
            );
        }
    });
}

#[test]
fn g0_g4_no_vibration_dofs_or_cancellation_returns_no_quadratic_pencil() {
    let law = material(YOUNG, 0.0);
    let state = thermal(&law, STRAIN);
    let field = TetThermalStrainField::Uniform(&state);
    let (nodes, elements) = with_cx(|cx| promote(&block([1; 3], [1.0; 3]), cx));
    let meshes =
        with_cx(|cx| hex_support::FAMILIES.map(|f| hex_support::mesh([1; 3], [1.0; 3], f, cx)));
    let fixed: Vec<_> = (0..3 * nodes.len()).collect();
    let problem = QuadraticTetProblem {
        nodes_m: &nodes,
        elements: &elements,
        materials: TetMaterialField::Uniform(&law),
        fixed_dofs: &fixed,
        budget: TetAssemblyBudget::standard(),
        quadrature: tet_support::Q,
    };
    with_cx(|cx| {
        assert!(matches!(
            problem.solve_thermal_tangent(field, TetStaticSolveConfig::default(), FRAME, 0.02, cx,),
            Err(TetElasticError::NoFreeDofs)
        ));
        for mesh in &meshes {
            let fixed: Vec<_> = (0..3 * mesh.nodes_m.len()).collect();
            assert!(matches!(
                hex_support::problem(mesh, &law, &fixed).solve_thermal_tangent(
                    field,
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                ),
                Err(TetElasticError::NoFreeDofs)
            ));
        }
    });
    let gate = fs_exec::CancelGate::new_clock_free();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 91,
                kernel_id: 91,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        assert!(matches!(
            problem
                .solve_thermal_tangent(field, TetStaticSolveConfig::default(), FRAME, 0.02, &cx,),
            Err(TetElasticError::Cancelled)
        ));
        for mesh in &meshes {
            let fixed = supports(&mesh.nodes_m, true);
            assert!(matches!(
                hex_support::problem(mesh, &law, &fixed).solve_thermal_tangent(
                    field,
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    &cx,
                ),
                Err(TetElasticError::Cancelled)
            ));
        }
    });
}
