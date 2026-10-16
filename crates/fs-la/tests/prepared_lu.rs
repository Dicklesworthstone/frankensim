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

#[test]
fn retained_factors_preserve_pivots_across_many_rhs_and_replacements() {
    for n in [0,1,2,7,33] {
        let mut a=vec![0.0;n*n];
        for i in 0..n {for j in 0..n {
            a[i*n+j]=if i==j {4.0} else {((i*7+j*13)%11) as f64/100.0};
        }}
        // A sequence of row swaps, not only a single pivot at column zero.
        for i in 0..n/2 {for j in 0..n {a.swap(i*n+j,(n-1-i)*n+j);}}
        let mut cached=LuWorkspace::new(n).unwrap();cached.factor(&a).unwrap();
        let mut fresh=LuWorkspace::new(n).unwrap();
        for seed in 0..16 {
            let b:Vec<_>=(0..n).map(|i|((i+seed)*17%19) as f64-9.0).collect();
            let (mut x,mut y)=(vec![0.0;n],vec![0.0;n]);
            cached.solve_factored_into(&b,&mut x).unwrap();
            fresh.solve_into(&a,&b,&mut y).unwrap();assert_eq!(x,y);
            let mut reference=b.clone();lu(&a,n).unwrap().solve(&mut reference);
            for (x,y) in x.iter().zip(reference) {assert!((x-y).abs()<1e-12);}
        }
    }
}

#[test]
fn failed_factor_cannot_reuse_stale_matrix_but_bad_rhs_does_not_poison_factor() {
    let mut work=LuWorkspace::new(2).unwrap();let mut out=[19.0,-17.0];
    assert_eq!(work.solve_factored_into(&[1.0,2.0],&mut out),Err(LuWorkspaceError::Unfactored));
    let a=[0.0,2.0,1.0,0.0];work.factor(&a).unwrap();
    assert_eq!(work.solve_factored_into(&[f64::NAN,1.0],&mut out),Err(LuWorkspaceError::NonFinite));
    assert_eq!(out,[19.0,-17.0]);
    work.solve_factored_into(&[6.0,4.0],&mut out).unwrap();assert_eq!(out,[4.0,3.0]);
    for bad in [&[1.0][..],&[1.0,2.0,2.0,4.0],&[1.0,0.0,0.0,f64::INFINITY]] {
        work.factor(&a).unwrap();assert!(work.factor(bad).is_err());
        assert_eq!(work.solve_factored_into(&[6.0,4.0],&mut out),Err(LuWorkspaceError::Unfactored));
        assert_eq!(out,[4.0,3.0]);
    }
    work.factor(&a).unwrap();work.solve_factored_into(&[8.0,5.0],&mut out).unwrap();assert_eq!(out,[5.0,4.0]);
}
