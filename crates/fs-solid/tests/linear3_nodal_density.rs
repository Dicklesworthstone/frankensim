//! Exact graded reference inertia through the native linear/modal pencil.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_modal::eigh_gen_dense;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticError, TetElasticMaterial, TetLinearElasticProblem,
    TetMaterialField,
};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const DENSITY: [f64; 4] = [1000.0, 1000.0, 1000.0, 5000.0];

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 7,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn material() -> TetElasticMaterial {
    TetElasticMaterial::try_new(2000.0, 12e9, 0.25, [7; 32]).unwrap()
}

fn problem<'a>(
    material: &'a TetElasticMaterial,
    fixed: &'a [usize],
) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard(),
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn density_location_changes_native_modes_at_identical_total_mass_and_stiffness() {
    with_cx(|cx, _| {
        let material = material();
        // Fix the base; all three free components have shape N_3.
        let p = problem(&material, &[0, 1, 2, 3, 4, 5, 6, 7, 8]);
        let uniform = p.assemble(cx).unwrap();
        let graded = p
            .assemble_with_nodal_reference_density(&DENSITY, cx)
            .unwrap();
        close(graded.total_mass_kg, uniform.total_mass_kg);
        assert_eq!(graded.stiffness.to_dense(), uniform.stiffness.to_dense());
        assert_eq!(graded.free_dofs, [9, 10, 11]);
        let k = graded.stiffness.to_dense();
        let m = graded.mass.to_dense();
        let modes = eigh_gen_dense(&k, &m, 3).unwrap();
        let reference = eigh_gen_dense(&k, &uniform.mass.to_dense(), 3).unwrap();
        // Integral rho*N_3² = 50 kg, versus 100/3 kg with the same mean rho.
        // Thus every eigenvalue must fall by exactly 2/3, despite equal mass.
        for (mode, before) in modes.iter().zip(reference) {
            close(mode.lambda / before.lambda, 2.0 / 3.0);
            close(
                mode.lambda.sqrt() / before.lambda.sqrt(),
                (2.0_f64 / 3.0).sqrt(),
            );
            assert!(mode.residual <= 1e-12 * mode.lambda);
            let mut momentum = vec![0.0; 3];
            graded.mass.spmv(&mode.phi, &mut momentum);
            close(mode.phi.iter().zip(momentum).map(|(v, p)| v * p).sum(), 1.0);
        }
    });
}

#[test]
fn assembly_retains_cubic_moments_and_uniform_bits_under_vertex_permutation() {
    with_cx(|cx, _| {
        let material = material();
        let mut p = problem(&material, &[]);
        let graded = p
            .assemble_with_nodal_reference_density(&DENSITY, cx)
            .unwrap();
        // Independent barycentric monomial integration: product(p_i!)/6!
        // on this unit tetrahedron, with total polynomial degree three.
        let factorial = [1.0, 1.0, 2.0, 6.0];
        for a in 0..4 {
            for b in 0..4 {
                let exact = DENSITY
                    .iter()
                    .enumerate()
                    .map(|(i, rho)| {
                        let mut powers = [0; 4];
                        powers[i] += 1;
                        powers[a] += 1;
                        powers[b] += 1;
                        rho * powers.into_iter().map(|p| factorial[p]).product::<f64>() / 720.0
                    })
                    .sum::<f64>();
                for c in 0..3 {
                    close(graded.mass.get(3 * a + c, 3 * b + c), exact);
                    close(graded.mass.get(3 * a + c, 3 * b + (c + 1) % 3), 0.0);
                }
            }
        }
        let permutations = [[3, 1, 0, 2]];
        p.tetrahedra = &permutations;
        let permuted = p
            .assemble_with_nodal_reference_density(&DENSITY, cx)
            .unwrap();
        for (a, b) in graded.mass.to_dense().iter().zip(permuted.mass.to_dense()) {
            close(*a, b);
        }
        let uniform = p.assemble(cx).unwrap();
        let nodal_uniform = p
            .assemble_with_nodal_reference_density(&[2000.0; 4], cx)
            .unwrap();
        assert_eq!(
            uniform.total_mass_kg.to_bits(),
            nodal_uniform.total_mass_kg.to_bits()
        );
        for (a, b) in uniform
            .mass
            .to_dense()
            .iter()
            .zip(nodal_uniform.mass.to_dense())
        {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(
            uniform.stiffness.to_dense(),
            nodal_uniform.stiffness.to_dense()
        );
    });
}

#[test]
fn invalid_density_budget_and_cancellation_refuse_before_publishing_a_pencil() {
    with_cx(|cx, gate| {
        let material = material();
        let mut p = problem(&material, &[]);
        for density in [
            vec![1.0; 3],
            vec![1.0, 2.0, 3.0, 0.0],
            vec![1.0, 2.0, 3.0, -1.0],
            vec![1.0, 2.0, 3.0, f64::NAN],
            vec![1.0, 2.0, 3.0, f64::INFINITY],
            vec![f64::from_bits(1); 4],
        ] {
            assert!(matches!(
                p.assemble_with_nodal_reference_density(&density, cx),
                Err(TetElasticError::InvalidMaterial { .. })
            ));
        }
        let large_nodes = NODES.map(|p| p.map(|x| 12.0 * x));
        p.nodes_m = &large_nodes;
        assert!(matches!(
            p.assemble_with_nodal_reference_density(&[f64::MAX; 4], cx),
            Err(TetElasticError::InvalidMaterial { .. })
        ));
        p.budget.maximum_nodes = 3;
        assert!(matches!(
            p.assemble_with_nodal_reference_density(&DENSITY, cx),
            Err(TetElasticError::BudgetExceeded { what: "nodes", .. })
        ));
        gate.request();
        assert!(matches!(
            p.assemble_with_nodal_reference_density(&DENSITY, cx),
            Err(TetElasticError::Cancelled)
        ));
    });
}
