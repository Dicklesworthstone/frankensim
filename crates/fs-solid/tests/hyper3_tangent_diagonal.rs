//! Native material preconditioner diagonals without global stiffness assembly.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::{HyperTetProblem, HyperTetSettings, HyperelasticLawField, TetAssemblyBudget};

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 17,
                    kernel_id: 97,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}

#[test]
fn local_diagonal_matches_original_heterogeneous_stiffness_and_unit_actions() {
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.1, 0.1, 0.0],
        [0.2, 1.0, 0.1],
        [0.0, 0.2, 1.2],
        [-0.8, 0.1, 0.0],
    ];
    let cells = [[0, 1, 2, 3], [0, 4, 3, 2]];
    let laws = [1.0, 3.0].map(|scale| {
        Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0 * scale,
                lambda: 4.0 * scale,
            },
            4.0,
        )
        .unwrap()
    });
    let problem = HyperTetProblem {
        nodes_m: &nodes,
        tetrahedra: &cells,
        material: HyperelasticLawField::PerElement(&laws),
        prescribed_m: &[(0, 0.03)],
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings::default(),
    };
    let displacement: Vec<_> = nodes
        .iter()
        .flat_map(|p| [0.03 + 0.2 * p[0] + 0.05 * p[1], -0.04 * p[1], 0.1 * p[2]])
        .collect();
    with_cx(|gate, cx| {
        let prepared = problem.prepare(cx).unwrap();
        let matrix = prepared.tangent_matrix(&displacement, cx).unwrap();
        let tangent = prepared
            .linearized_tangent(&displacement, 2 * 81 * 8, cx)
            .unwrap();
        let mut diagonal = [11.0; 15];
        let mut scratch = [f64::NAN; 15];
        tangent.diagonal(&mut diagonal, &mut scratch, cx).unwrap();
        for dof in 0..15 {
            let mut unit = [0.0; 15];
            unit[dof] = 1.0;
            let mut action = [0.0; 15];
            tangent.apply(&unit, &mut action, &mut scratch, cx).unwrap();
            for expected in [matrix.get(dof, dof), action[dof]] {
                assert!((diagonal[dof] - expected).abs() <= 3e-13 * (1.0 + expected.abs()));
            }
        }
        // This material-only observer retains prescribed rows unchanged.
        assert!(diagonal[0] > 0.0);
        let saved = diagonal;
        assert!(
            tangent
                .diagonal(&mut diagonal, &mut scratch[..14], cx)
                .is_err()
        );
        assert_eq!(diagonal, saved);
        gate.request();
        assert!(tangent.diagonal(&mut diagonal, &mut scratch, cx).is_err());
        assert_eq!(diagonal, saved);
    });
}
