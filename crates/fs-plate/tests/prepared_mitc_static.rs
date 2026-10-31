//! Repeated native static solves retain their physical support/operator owner.

use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, PreparedMitcStatic,
    ShearShellSection, ShellConstraint, ShellMass, assemble_mitc3p, solve_static,
};
use std::sync::Arc;

fn model(supported: bool) -> MitcShellModel {
    let mesh = ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
        vec![[0, 1, 3], [0, 3, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 4], 0.9).unwrap();
    let section = ShearShellSection::from_orthotropic_solid(
        10.0,
        10.0,
        0.3,
        10.0 / 2.6,
        10.0 / 2.6,
        10.0 / 2.6,
        0.1,
        2.0,
        0.0,
        5.0 / 6.0,
    )
    .unwrap();
    let constraints = if supported {
        vec![
            ShellConstraint {
                node: 0,
                mask: [true; 6],
                frame: ConstraintFrame::Global,
            },
            ShellConstraint {
                node: 3,
                mask: [false, true, false, false, false, false],
                frame: ConstraintFrame::Local([[0.8, 0.6, 0.0], [-0.6, 0.8, 0.0], [0.0, 0.0, 1.0]]),
            },
        ]
    } else {
        vec![]
    };
    assemble_mitc3p(
        &mesh,
        &directors,
        &vec![section; 2],
        &constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 4,
            max_triangles: 2,
        },
        1.0e-3,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_repeated_loads_recover_manufactured_world_motion_with_oblique_supports() {
    let original = Arc::new(model(true));
    let prepared = PreparedMitcStatic::try_new(Arc::clone(&original), || false).unwrap();
    assert!(Arc::ptr_eq(&original, &prepared.shared_model()));
    assert!(std::ptr::eq(original.as_ref(), prepared.model()));

    for (seed, amplitude) in [1.0, -0.25, 2.0, 0.0].into_iter().enumerate() {
        let reduced: Vec<_> = (0..original.pencil().free)
            .map(|i| amplitude * 1.0e-4 * (((i * 7 + seed * 11) % 19) as f64 - 9.0))
            .collect();
        let mut reduced_load = vec![0.0; reduced.len()];
        original.pencil().k.spmv(&reduced, &mut reduced_load);
        let world_load = original.expand(&reduced_load).unwrap();
        let expected = original.expand(&reduced).unwrap();
        let actual = prepared.solve(&world_load, || false).unwrap();
        let scale = expected.iter().map(|v| v.abs()).fold(1.0e-8_f64, f64::max);
        for (a, e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() <= 2.0e-8 * scale, "{a:e} != {e:e}");
        }
        // The borrowed compatibility API retains the same deterministic math.
        assert_eq!(actual, solve_static(&original, &world_load).unwrap());
    }
}

#[test]
fn g0_invalid_and_cancelled_candidates_leave_cached_factor_usable() {
    let prepared = PreparedMitcStatic::try_new(model(true), || false).unwrap();
    let mut load = vec![0.0; 24];
    load[8] = -0.001;
    load[20] = 0.0005;
    let expected = prepared.solve(&load, || false).unwrap();
    assert!(matches!(
        prepared.solve(&load[..23], || false),
        Err(PlateError::BadConstraint { .. })
    ));
    let mut invalid = load.clone();
    invalid[0] = f64::NAN;
    assert!(matches!(
        prepared.solve(&invalid, || false),
        Err(PlateError::BadConstraint { .. })
    ));
    assert_eq!(prepared.solve(&load, || true), Err(PlateError::Cancelled));
    let mut checkpoints = 0;
    assert_eq!(
        prepared.solve(&load, || {
            checkpoints += 1;
            checkpoints >= 3
        }),
        Err(PlateError::Cancelled)
    );
    assert_eq!(prepared.solve(&load, || false).unwrap(), expected);
}

#[test]
fn g0_preparation_preserves_free_body_refusal_and_cancellation() {
    assert!(PreparedMitcStatic::try_new(model(false), || false).is_err());
    assert!(matches!(
        PreparedMitcStatic::try_new(model(true), || true),
        Err(PlateError::Cancelled)
    ));
}
