//! G1/G3: exact shared-face aliases preserve modes without assigning pressures.

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_sparse::Coo;

fn with_cx(run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 193,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn close(a: f64, b: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= 1e-11 * a.abs().max(b.abs()).max(1.0),
        "{a} versus {b}"
    );
}

fn diagonal(values: &[f64]) -> Csr {
    let mut coo = Coo::new(values.len(), values.len());
    for (i, &value) in values.iter().enumerate() {
        coo.push(i, i, value);
    }
    coo.assemble()
}

fn pencil(free: &[usize], k: &[f64], m: &[f64], full_dofs: usize) -> FrozenBrittlePencil {
    FrozenBrittlePencil {
        stiffness: diagonal(k),
        mass: diagonal(m),
        free_dofs: free.to_vec(),
        full_dofs,
        reference_mass_kg: m.iter().sum(),
        element_damage: vec![],
    }
}

fn row(face_id: u64, nodes: [usize; 2], direction: [f64; 3]) -> BrittleModeConstraint {
    BrittleModeConstraint {
        row: BrittleFaceImpactRow {
            face_id,
            vertex: 0,
            nodes,
        },
        direction,
    }
}

#[test]
fn shared_faces_reversed_pairs_and_support_aliases_keep_the_physical_spectrum() {
    with_cx(|cx| {
        // Three scalar masses [2,3,5] kg with ground springs [3,5,7] N/m.
        // Closing u0=u1 leaves exact eigenvalues 7/5 and (3+5)/(2+3).
        let make = || pencil(&[2, 5, 8], &[3.0, 5.0, 7.0], &[2.0, 3.0, 5.0], 9);
        let normal = row(10, [0, 1], [0.0, 0.0, 1.0]);
        let original = vec![
            normal,
            row(11, [0, 1], [0.0, 0.0, 1.0]),
            row(12, [1, 0], [0.0, 0.0, -1.0]),
            row(13, [1, 0], [0.0, 0.0, 1.0]),
            row(14, [0, 1], [0.0, 0.0, -1.0]),
            // Its physical x row is wholly carried by the supports.
            row(16, [0, 2], [1.0, 0.0, 0.0]),
        ];
        let baseline =
            project_contact_pencil(make(), vec![normal], Default::default(), cx).unwrap();
        let result =
            project_contact_pencil(make(), original.clone(), Default::default(), cx).unwrap();
        assert_eq!(result.constraints(), &[normal]);
        assert_eq!(result.original_constraints(), original);
        assert!(baseline.individual_reactions_determined());
        assert!(!result.individual_reactions_determined());
        assert_eq!(result.dimension(), 2);
        assert_eq!(
            result.projected_stiffness_n_m(),
            baseline.projected_stiffness_n_m()
        );
        assert_eq!(result.projected_mass_kg(), baseline.projected_mass_kg());
        for (alias, sign) in result
            .constraint_aliases()
            .iter()
            .zip([1., 1., 1., -1., -1., 0.])
        {
            assert_eq!(
                alias.independent_constraint,
                if sign == 0. { None } else { Some(0) }
            );
            assert_eq!(alias.orientation, sign);
        }
        for ((mode, reference), expected) in
            result.modes().iter().zip(baseline.modes()).zip([1.4, 1.6])
        {
            close(mode.coordinates.lambda, expected);
            assert_eq!(
                mode.constraint_force_shape,
                reference.constraint_force_shape
            );
            let phi = &mode.displacement_shape;
            close(
                2.0 * phi[2] * phi[2] + 3.0 * phi[5] * phi[5] + 5.0 * phi[8] * phi[8],
                1.0,
            );
            for original in result.original_constraints() {
                let [a, b] = original.row.nodes;
                let jump: f64 = (0..3)
                    .map(|j| original.direction[j] * (phi[3 * b + j] - phi[3 * a + j]))
                    .sum();
                close(jump, 0.0);
            }
            // Only this aggregate is determined: J^T df = K phi-lambda M phi.
            close(
                mode.constraint_force_shape[0],
                (5.0 - 3.0 * expected) * phi[5],
            );
            assert_eq!(mode.constraint_force_shape.len(), 1);
        }
        let capped = FrozenBrittleContactControls {
            maximum_contact_rows: 1,
            ..Default::default()
        };
        assert!(
            project_contact_pencil(make(), original, capped, cx).is_err(),
            "incoming rows must fit the cap before any alias reduction"
        );
    });
}

#[test]
fn sticking_rows_alias_componentwise_and_empty_physical_space_retains_witnesses() {
    with_cx(|cx| {
        let make = || pencil(&[0, 1, 2, 3, 4, 5], &[3.0; 6], &[2.0; 6], 6);
        let directions = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let one_face: Vec<_> = directions
            .into_iter()
            .map(|direction| row(20, [0, 1], direction))
            .collect();
        let mut shared = one_face.clone();
        shared.extend(
            directions
                .into_iter()
                .map(|direction| row(21, [1, 0], direction.map(|v| -v))),
        );
        let result = project_contact_pencil(make(), shared, Default::default(), cx).unwrap();
        assert_eq!(result.constraints().len(), 3);
        assert_eq!(result.original_constraints().len(), 6);
        assert_eq!(result.dimension(), 3);
        for mode in result.modes() {
            close(mode.coordinates.lambda, 1.5);
        }
        for (index, alias) in result.constraint_aliases().iter().enumerate() {
            assert_eq!(alias.independent_constraint, Some(index % 3));
            assert_eq!(alias.orientation, 1.0);
        }
        // Node 1 is now supported: the one remaining displacement is closed.
        let closed = project_contact_pencil(
            pencil(&[2], &[3.0], &[2.0], 9),
            vec![
                row(30, [0, 1], [0., 0., 1.]),
                row(31, [1, 0], [0., 0., 1.]),
                row(32, [1, 2], [0., 0., 1.]),
                row(33, [0, 2], [0., 0., 1.]),
            ],
            Default::default(),
            cx,
        )
        .unwrap();
        assert_eq!(closed.dimension(), 0);
        assert!(closed.modes().is_empty());
        assert!(closed.basis().is_empty());
        assert_eq!(closed.original_constraints().len(), 4);
        assert_eq!(closed.constraint_aliases()[1].orientation, -1.0);
        assert_eq!(closed.constraint_aliases()[2].independent_constraint, None);
        // Different held partners give identical free covectors, but cannot
        // identify which supported node receives the original contact force.
        assert_eq!(
            closed.constraint_aliases()[3].independent_constraint,
            Some(0)
        );
        assert_eq!(closed.constraint_aliases()[3].orientation, 1.0);
        assert!(!closed.individual_reactions_determined());
    });
}

#[test]
fn near_duplicates_and_nonalias_dependencies_are_not_silently_truncated() {
    with_cx(|cx| {
        let nearly = project_contact_pencil(
            pencil(&[0, 1, 2, 3, 4, 5], &[3.0; 6], &[2.0; 6], 6),
            vec![
                row(40, [0, 1], [1., 0., 0.]),
                row(41, [0, 1], [1., 1e-14, 0.]),
            ],
            Default::default(),
            cx,
        );
        assert!(matches!(
            nearly,
            Err(Error::Invalid(
                "unresolved contact constraint rank after exact alias reduction"
            ))
        ));
        // J02 = J01 + J12 is a true dependency, but not a signed duplicate.
        let loop_rows = project_contact_pencil(
            pencil(&[2, 5, 8], &[3.0; 3], &[2.0; 3], 9),
            vec![
                row(50, [0, 1], [0., 0., 1.]),
                row(51, [1, 2], [0., 0., 1.]),
                row(52, [0, 2], [0., 0., 1.]),
            ],
            Default::default(),
            cx,
        );
        assert!(matches!(
            loop_rows,
            Err(Error::Invalid(
                "unresolved contact constraint rank after exact alias reduction"
            ))
        ));
    });
}
