use fs_la::{LuWorkspace, LuWorkspaceError};
use fs_la::factor::lu;

#[test]
fn prepared_lu_matches_blocked_reference_across_panel_boundary() {
    for n in [0, 1, 3, 31, 32, 33, 65] {
        let mut a = vec![0.0; n * n];
        for row in 0..n {
            for col in 0..n {
                a[row * n + col] = if row == col { 5.0 } else {
                    ((row * 17 + col * 11 + 3) % 13) as f64 / (20.0 * n as f64)
                };
            }
        }
        // Force a nontrivial permutation without losing diagonal dominance.
        if n > 1 { for col in 0..n { a.swap(col, (n - 1) * n + col); } }
        let rhs: Vec<_> = (0..n).map(|i| (i as f64 + 1.0) / 7.0).collect();
        let mut expected = rhs.clone();
        lu(&a, n).unwrap().solve(&mut expected);
        let mut workspace = LuWorkspace::new(n).unwrap();
        let mut result = vec![0.0; n];
        for _ in 0..4 {
            workspace.solve_into(&a, &rhs, &mut result).unwrap();
            for (actual, expected) in result.iter().zip(&expected) {
                assert!((actual - expected).abs() <= 2e-12 * (1.0 + expected.abs()));
            }
        }
    }
}

#[test]
fn prepared_lu_refuses_transactionally_and_recovers() {
    let mut workspace = LuWorkspace::new(2).unwrap();
    let mut result = [17.0, -19.0];
    assert_eq!(workspace.solve_into(&[1.0], &[1.0, 2.0], &mut result), Err(LuWorkspaceError::Dimension));
    assert!(matches!(workspace.solve_into(&[1.0, 2.0, 2.0, 4.0], &[1.0, 2.0], &mut result), Err(LuWorkspaceError::Singular { .. })));
    assert_eq!(workspace.solve_into(&[f64::NAN, 0.0, 0.0, 1.0], &[1.0, 2.0], &mut result), Err(LuWorkspaceError::NonFinite));
    assert_eq!(result, [17.0, -19.0]);
    workspace.solve_into(&[0.0, 2.0, 1.0, 0.0], &[6.0, 4.0], &mut result).unwrap();
    assert_eq!(result, [4.0, 3.0]);
}

#[test]
fn prepared_lu_extent_overflow_is_typed() {
    assert!(matches!(LuWorkspace::new(usize::MAX), Err(LuWorkspaceError::Capacity)));
}
