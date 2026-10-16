//! Sparse inverse certificates against independent closed-form systems.
use fs_solver::goal::{GoalResidualError, GoalResidualLimits, enclose_goal_error};
use fs_solver::goal::inverse::sparse::{SparseInverseLimits, SparseInverseStatus, VerifiedSparseInverse};
use fs_sparse::{Coo, Csr};

fn limits() -> SparseInverseLimits {
    SparseInverseLimits { max_rows: 2_000, max_input_nonzeros: 50_000,
        max_entries: 100_000, max_updates: 200_000 }
}
fn residual_limits() -> GoalResidualLimits {
    GoalResidualLimits { max_rows: 2_000, max_nonzeros: 50_000 }
}
fn blocks(count: usize, scale: f64) -> Csr {
    // Each block is I + 3*ones(3,3): SPD, but its comparison matrix is not.
    let mut a = Coo::new(3*count, 3*count);
    for block in 0..count {
        for i in 0..3 { for j in 0..3 {
            a.push(3*block+i, 3*block+j, scale * if i==j { 4.0 } else { 3.0 });
        } }
    }
    a.assemble()
}
fn prepare(a: &Csr) -> VerifiedSparseInverse {
    VerifiedSparseInverse::prepare(a, limits(), || true).unwrap()
}

#[test]
fn sparse_nondominant_system_exceeds_dense_inverse_row_limit() {
    let a = blocks(300, 1.0); // 900 unknowns, without a 900-column inverse.
    let zero = vec![0.0; a.nrows()];
    let old = enclose_goal_error(&a, &zero, &zero, &zero, &zero,
        None, residual_limits(), || true).unwrap();
    assert!(old.inverse_infinity_upper().is_none());
    let checked = prepare(&a);
    assert_eq!(checked.status(), SparseInverseStatus::Enclosed);
    // (I+3 J)^-1 = I - 0.3 J, with row absolute sum exactly 13/10.
    assert!(checked.inverse_upper().unwrap() >= 1.3);
    assert!(checked.inverse_upper().unwrap() < 10.0);
    assert_eq!(checked.work().eliminated_rows, 900);
    assert!(checked.work().factor_entries < 3_000);
    assert!(checked.work().updates < 3_000);
    assert!(checked.work().peak_entries <= limits().max_entries);
}

#[test]
fn sparse_path_bound_uses_triangular_sweeps_not_a_product_of_global_norms() {
    let n = 513;
    let mut a = Coo::new(n,n);
    for i in 0..n {
        a.push(i,i,2.0);
        if i>0 { a.push(i,i-1,1.0); }
        if i+1<n { a.push(i,i+1,1.0); }
    }
    let checked = prepare(&a.assemble());
    // Sign similarity to tridiag(-1,2,-1) preserves |A^-1|. Its row sum is
    // i*(n+1-i)/2; n=513 gives an exactly representable maximum 33024.5.
    let exact = 33024.5;
    assert!(checked.inverse_upper().unwrap() >= exact);
    assert!(checked.inverse_upper().unwrap() < exact*1.01);
    assert_eq!(checked.work().updates, n-1);
}

#[test]
fn stored_asymmetry_is_bounded_and_matrix_substitution_refuses() {
    let a = Csr::from_parts(2,2,vec![0,2,4],vec![0,1,0,1],vec![3.0,1.125,1.0,2.0]);
    let checked = prepare(&a);
    let exact = 4.0/(6.0-1.125);
    assert!(checked.inverse_upper().unwrap() >= exact);
    assert!(checked.perturbation_upper().unwrap() >= 0.125);
    let zero = [0.0;2];
    checked.enclose_goal(&a,&zero,&zero,&zero,&zero,None,residual_limits(),||true).unwrap();
    let changed = Csr::from_parts(2,2,vec![0,2,4],vec![0,1,0,1],vec![3.0,1.25,1.0,2.0]);
    assert!(matches!(checked.enclose_goal(&changed,&zero,&zero,&zero,&zero,None,
        residual_limits(),||true), Err(GoalResidualError::MatrixMismatch)));
    // The unmirrored upper coefficient may not be silently dropped either.
    let upper = Csr::from_parts(2,2,vec![0,2,3],vec![0,1,1],vec![1.0,2.0,1.0]);
    let refused = prepare(&upper);
    assert_eq!(refused.status(),SparseInverseStatus::PerturbationNotEstablished);
    assert!(refused.inverse_upper().is_none());
}

#[test]
fn cached_goal_enclosure_retains_inexact_dual_error_and_replays() {
    let a = blocks(1,1.0);
    let checked = prepare(&a);
    let rhs = [19.0,20.0,21.0]; // exact solution [1,2,3]
    let g = [1.0,-1.0,0.5]; // exact goal 0.5
    let z = [0.0;3]; // deliberately unsolved dual, not exact authority
    let first = checked.enclose_goal(&a,&rhs,&z,&g,&z,None,residual_limits(),||true).unwrap();
    let error = first.goal_error().unwrap();
    assert!(error.lower() <= 0.5 && error.upper() >= 0.5);
    assert!(first.dual_error_upper().unwrap() > 0.0);
    assert_eq!(first, checked.enclose_goal(&a,&rhs,&z,&g,&z,None,residual_limits(),||true).unwrap());
    let work = checked.work();
    let solved = checked.enclose_goal(&a,&rhs,&[1.0,2.0,3.0],&g,&z,None,
        residual_limits(),||true).unwrap();
    assert!(solved.goal_error().unwrap().magnitude_upper() < 1e-10);
    assert_eq!(checked.work(),work);
}

#[test]
fn scaling_singular_pivots_and_finite_range_never_create_false_authority() {
    for exponent in [-450,-200,0,200,450] {
        let scale = 2.0_f64.powi(exponent);
        let checked = prepare(&blocks(1,scale));
        assert_eq!(checked.status(),SparseInverseStatus::Enclosed);
        assert!(checked.inverse_upper().unwrap()*scale >= 1.3);
        assert!(checked.inverse_upper().unwrap()*scale < 10.0);
    }
    let singular = Csr::from_parts(2,2,vec![0,2,4],vec![0,1,0,1],vec![1.0;4]);
    let checked = prepare(&singular);
    assert_eq!(checked.status(),SparseInverseStatus::PivotNotPositive);
    assert!(checked.inverse_upper().is_none());
    let tiny = Csr::from_parts(1,1,vec![0,1],vec![0],vec![f64::from_bits(1)]);
    assert_eq!(prepare(&tiny).status(),SparseInverseStatus::ArithmeticRange);
}

#[test]
fn sparse_fill_and_update_caps_stop_without_a_partial_inverse() {
    let a = blocks(4,1.0);
    for cap in [0,1,4] {
        let checked = VerifiedSparseInverse::prepare(&a,
            SparseInverseLimits { max_updates:cap, ..limits() },||true).unwrap();
        assert_eq!(checked.status(),SparseInverseStatus::UpdateLimit);
        assert!(checked.work().updates<=cap);
        assert!(checked.inverse_upper().is_none());
    }
    for cap in [0,1,10,36] {
        let checked = VerifiedSparseInverse::prepare(&a,
            SparseInverseLimits { max_entries:cap, ..limits() },||true).unwrap();
        assert_eq!(checked.status(),SparseInverseStatus::EntryLimit);
        assert!(checked.work().peak_entries<=cap);
        assert!(checked.inverse_upper().is_none());
    }
    assert!(matches!(VerifiedSparseInverse::prepare(&a,
        SparseInverseLimits { max_rows:1, ..limits() },||true),Err(GoalResidualError::Limit{..})));
}

#[test]
fn every_cancellation_boundary_refuses_and_exact_retry_is_unchanged() {
    let a = blocks(12,1.0);
    let original = a.clone();
    let mut calls = 0;
    let expected = VerifiedSparseInverse::prepare(&a,limits(),||{calls+=1;true}).unwrap();
    assert!(calls>3);
    for stop in 1..=calls {
        let mut seen=0;
        assert!(matches!(VerifiedSparseInverse::prepare(&a,limits(),||{seen+=1;seen!=stop}),
            Err(GoalResidualError::Cancelled)));
    }
    assert_eq!(a,original);
    let retry=prepare(&a);
    assert_eq!(expected.inverse_upper(),retry.inverse_upper());
    assert_eq!(expected.work(),retry.work());
    let zero=vec![0.0;a.nrows()];
    assert!(matches!(retry.enclose_goal(&a,&zero,&zero,&zero,&zero,None,residual_limits(),||false),
        Err(GoalResidualError::Cancelled)));
}
