//! Coupled bounds on a solid that no positive comparison scaling can certify.
//! A=I+11^T, A^-1=I-11^T/4; independent closed-form rank-one oracles below.
use fs_solver::goal::{GoalResidualError, GoalResidualLimits};
use fs_solver::goal::feedback::{
    FeedbackBoundStatus, FeedbackInverseMethod, FeedbackResidualLimits, FeedbackResidualReport,
    enclose_affine_feedback_error_with_inverse, enclose_affine_feedback_error_with_schur,
};
use fs_sparse::{Coo, Csr};

fn matrix(rows: &[&[f64]]) -> Csr {
    let mut coo = Coo::new(rows.len(), rows[0].len());
    for (i, row) in rows.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            if value != 0.0 { coo.push(i, j, value); }
        }
    }
    coo.assemble()
}
fn solid() -> Csr { matrix(&[&[2., 1., 1.], &[1., 2., 1.], &[1., 1., 2.]]) }
fn inverse() -> Vec<Vec<f64>> {
    vec![vec![0.75, -0.25, -0.25], vec![-0.25, 0.75, -0.25], vec![-0.25, -0.25, 0.75]]
}
fn limits() -> FeedbackResidualLimits {
    FeedbackResidualLimits {
        solid: GoalResidualLimits { max_rows: 16, max_nonzeros: 10_000 },
        max_ports: 4, max_transfer_nonzeros: 32, max_response_entries: 32,
        max_verification_entries: 10_000,
    }
}
fn run(c: f64, columns: &[Vec<f64>], response: &[Vec<f64>], cap: usize) -> Result<FeedbackResidualReport, GoalResidualError> {
    let mut lim = limits(); lim.max_verification_entries = cap;
    enclose_affine_feedback_error_with_inverse(
        &solid(), &[1., 2., 3.], &[0., 0., 0.],
        &matrix(&[&[1.], &[1.], &[1.]]), &matrix(&[&[c, c, c]]), &[0.25],
        Some(response), None, columns, lim, || true,
    )
}
fn oracle(report: &FeedbackResidualReport, c: f64, system_scale: f64) {
    // M=I+(1-c)11^T, load=(1,2,3)+1/4. No production solve/checker is used.
    let k = 1.0-c;
    let beta = k/(1.0+3.0*k);
    let expected: Vec<f64> = [1.25, 2.25, 3.25].iter().map(|b| b-beta*6.75).collect();
    let error = expected.iter().map(|x| x.abs()).fold(0.0, f64::max);
    assert_eq!(report.status(), FeedbackBoundStatus::Enclosed, "{report:?}");
    assert!(report.state_error_infinity_upper().unwrap() >= error);
    let norm = (1.0-beta).abs()+2.0*beta.abs();
    assert!(report.coupled_inverse_infinity_upper().unwrap()*system_scale >= norm);
    assert!(report.state_error_infinity_upper().unwrap() < 20.0);
}

#[test]
fn nondominant_spd_solid_uses_checked_inverse_for_complete_feedback() {
    let a = solid(); let b = matrix(&[&[1.], &[1.], &[1.]]);
    let c = matrix(&[&[0.5, 0.5, 0.5]]); let response = vec![vec![0.25; 3]];
    let old = enclose_affine_feedback_error_with_schur(
        &a, &[1., 2., 3.], &[0.; 3], &b, &c, &[0.25], Some(&response), None, limits(), || true,
    ).unwrap();
    assert_eq!(old.status(), FeedbackBoundStatus::SolidInverseUnavailable);
    let report = run(0.5, &inverse(), &response, 10_000).unwrap();
    oracle(&report, 0.5, 1.0);
    assert_eq!(report.inverse_method(), Some(FeedbackInverseMethod::StateContraction));
    assert!(report.solid_inverse_infinity_upper().unwrap() >= 1.25);
    assert_eq!(report.response_residual_infinity_upper().len(), 1);
}

#[test]
fn inverse_evidence_reaches_schur_without_assuming_state_contraction() {
    let report = run(-4.0, &inverse(), &[vec![0.25; 3]], 10_000).unwrap();
    oracle(&report, -4.0, 1.0);
    assert!(report.gain_infinity_upper().unwrap() >= 3.0);
    assert_eq!(report.inverse_method(), Some(FeedbackInverseMethod::PortSchurDominance));
    assert!(report.schur_inverse_infinity_upper().unwrap() >= 0.25);
}

#[test]
fn inaccurate_inverse_and_response_columns_cannot_mint_authority() {
    let zero_inverse = vec![vec![0.0; 3]; 3];
    let report = run(0.5, &zero_inverse, &[vec![0.25; 3]], 10_000).unwrap();
    assert_eq!(report.status(), FeedbackBoundStatus::SolidInverseUnavailable);
    assert!(report.state_error_infinity_upper().is_none());
    let mut doubled = inverse();
    for column in &mut doubled { for value in column { *value *= 2.0; } }
    assert!(run(0.5, &doubled, &[vec![0.25; 3]], 10_000).unwrap().state_error_infinity_upper().is_none());
    let bad_response = run(0.5, &inverse(), &[vec![0.0; 3]], 10_000).unwrap();
    assert!(bad_response.response_residual_infinity_upper()[0] >= 1.0);
    assert!(bad_response.state_error_infinity_upper().is_none());
    // Exact singular feedback: C H = 1, despite an invertible solid.
    let singular = enclose_affine_feedback_error_with_inverse(
        &solid(), &[1.; 3], &[0.; 3], &matrix(&[&[1.], &[1.], &[1.]]),
        &matrix(&[&[4., 0., 0.]]), &[0.], Some(&[vec![0.25; 3]]), None,
        &inverse(), limits(), || true,
    ).unwrap();
    assert!(singular.state_error_infinity_upper().is_none());
}

#[test]
fn inverse_and_schur_share_one_preflighted_verification_budget() {
    // Base=24, inverse validation=9, inverse verification=36. Schur adds 13.
    let columns = inverse(); let response = [vec![0.25; 3]];
    assert!(matches!(run(0.5, &columns, &response, 68), Err(GoalResidualError::Limit { field: "verification entries", .. })));
    oracle(&run(0.5, &columns, &response, 69).unwrap(), 0.5, 1.0);
    for cap in [69, 81] {
        let report = run(-4.0, &columns, &response, cap).unwrap();
        assert!(report.state_error_infinity_upper().is_none(), "{report:?}");
        assert_eq!(report.status(), FeedbackBoundStatus::ContractionNotEstablished);
    }
    oracle(&run(-4.0, &columns, &response, 82).unwrap(), -4.0, 1.0);
}

#[test]
fn malformed_proposals_refuse_and_dimension_limits_are_not_bypassed() {
    let response = [vec![0.25; 3]];
    assert!(matches!(run(0.5, &inverse()[..2], &response, 10_000), Err(GoalResidualError::Length { field: "inverse columns", .. })));
    let mut bad = inverse(); bad[1].pop();
    assert!(matches!(run(0.5, &bad, &response, 10_000), Err(GoalResidualError::Length { field: "inverse column", .. })));
    bad = inverse(); bad[1][1] = f64::NAN;
    assert!(matches!(run(0.5, &bad, &response, 10_000), Err(GoalResidualError::NonFinite { field: "inverse column", .. })));
}

#[test]
fn cancellation_in_each_phase_refuses_publication_and_retry_is_exact() {
    let a = solid(); let b = matrix(&[&[1.], &[1.], &[1.]]); let c = matrix(&[&[-4.; 3]]);
    let columns = inverse(); let saved = columns.clone(); let response = [vec![0.25; 3]];
    let check = |poll: &mut dyn FnMut() -> bool| enclose_affine_feedback_error_with_inverse(
        &a, &[1., 2., 3.], &[0.; 3], &b, &c, &[0.25], Some(&response), None,
        &columns, limits(), poll,
    );
    let mut calls = 0;
    let expected = check(&mut || { calls += 1; true }).unwrap();
    assert!(calls > 5);
    for fail_at in 1..=calls {
        let mut seen = 0;
        assert!(matches!(check(&mut || { seen += 1; seen < fail_at }), Err(GoalResidualError::Cancelled)));
    }
    assert_eq!(columns, saved);
    assert_eq!(check(&mut || true).unwrap(), expected);
}

#[test]
fn physical_system_scaling_keeps_a_real_coupled_bound() {
    for power in [-450, -100, 0, 100, 450] {
        let scale = 2.0_f64.powi(power);
        let a = matrix(&[&[2.*scale, scale, scale], &[scale, 2.*scale, scale], &[scale, scale, 2.*scale]]);
        let b = matrix(&[&[scale], &[scale], &[scale]]);
        let columns: Vec<Vec<_>> = inverse().iter().map(|col| col.iter().map(|x| x/scale).collect()).collect();
        let report = enclose_affine_feedback_error_with_inverse(
            &a, &[scale, 2.*scale, 3.*scale], &[0.; 3], &b,
            &matrix(&[&[0.5; 3]]), &[0.25], Some(&[vec![0.25; 3]]), None,
            &columns, limits(), || true,
        ).unwrap();
        oracle(&report, 0.5, scale);
    }
}
