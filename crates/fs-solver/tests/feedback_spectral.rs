//! The full coupled equation must retain response errors after sparse proof.
use fs_solver::goal::{GoalResidualLimits, GoalResidualError};
use fs_solver::goal::feedback::{FeedbackBoundStatus, FeedbackInverseMethod,
    FeedbackResidualLimits, enclose_affine_feedback_error_with_schur,
    enclose_affine_feedback_error_with_spectral};
use fs_solver::goal::inverse::spectral::{SpectralError, SpectralInverseLimits, SpectralStop};
use fs_sparse::{Coo, Csr};

fn solid(n: usize) -> Csr {
    let mut a = Coo::new(n, n);
    for i in 0..n {
        for j in (i / 3) * 3..(i / 3) * 3 + 3 {
            a.push(i, j, if i == j { 3.0 } else { 1.5 });
        }
    }
    a.assemble()
}
fn limits(n: usize) -> FeedbackResidualLimits {
    FeedbackResidualLimits { solid: GoalResidualLimits { max_rows: n, max_nonzeros: 3*n },
        max_ports: 1, max_transfer_nonzeros: 2*n, max_response_entries: n,
        max_verification_entries: 1_000_000 }
}
fn spectral(n: usize) -> SpectralInverseLimits {
    SpectralInverseLimits { system: limits(n).solid, max_storage_entries: 40*n, max_work_entries: 500_000, max_shift_attempts: 8 }
}
fn transfer(n: usize, feedback_scale: f64) -> (Csr, Csr) {
    let mut b = Coo::new(n, 1);
    let mut c = Coo::new(1, n);
    for i in 0..n { b.push(i, 0, 1.0); c.push(0, i, feedback_scale / n as f64); }
    (b.assemble(), c.assemble())
}

#[test]
fn full_coupled_error_is_enclosed_above_the_dense_inverse_limit() {
    let n = 300;
    let a = solid(n);
    let (b,c) = transfer(n, 1.0);
    let rhs = vec![5.0; n];
    let candidate = vec![0.0; n];
    let responses = [vec![1.0 / 6.0; n]];
    let old = enclose_affine_feedback_error_with_schur(
        &a,&rhs,&candidate,&b,&c,&[0.0],Some(&responses),None,limits(n),||true,
    ).unwrap();
    assert_eq!(old.status(), FeedbackBoundStatus::SolidInverseUnavailable);
    let result = enclose_affine_feedback_error_with_spectral(
        &a,&rhs,&candidate,&b,&c,&[0.0],Some(&responses),None,limits(n),spectral(n),||true,
    ).unwrap();
    assert_eq!(result.status(), FeedbackBoundStatus::Enclosed);
    assert_eq!(result.inverse_method(), Some(FeedbackInverseMethod::StateContraction));
    assert!(result.solid_spectral().is_some());
    // A*1=6*1, BC*1=1: the exact constant solution is 1 apart from the
    // tiny representation error in 1/n. A 1e-12 outward oracle band covers it.
    assert!(result.state_error_infinity_upper().unwrap() >= 1.0 - 1e-12);
    assert_eq!(result.response_residual_infinity_upper().len(),1);
}

#[test]
fn sparse_proof_does_not_erase_inaccurate_response_columns() {
    let n = 3;
    let a = solid(n);let(b,c)=transfer(n,1.0);
    let bad = [vec![0.0;n]];
    let checked = enclose_affine_feedback_error_with_spectral(
        &a,&[5.;3],&[0.;3],&b,&c,&[0.],Some(&bad),None,limits(n),spectral(n),||true,
    ).unwrap();
    assert!(checked.solid_spectral().is_some());
    assert!(checked.response_residual_infinity_upper()[0]>=1.0);
    assert_ne!(checked.status(),FeedbackBoundStatus::Enclosed);
    assert!(checked.state_error_infinity_upper().is_none());
}

#[test]
fn negative_feedback_uses_the_error_enclosed_schur_route() {
    let n = 3;
    let a = solid(n);let(b,c)=transfer(n,-12.0);
    let responses=[vec![1.0/6.0;n]];
    let checked=enclose_affine_feedback_error_with_spectral(
        &a,&[18.;3],&[0.;3],&b,&c,&[0.],Some(&responses),None,limits(n),spectral(n),||true,
    ).unwrap();
    assert_eq!(checked.status(),FeedbackBoundStatus::Enclosed);
    assert_eq!(checked.inverse_method(),Some(FeedbackInverseMethod::PortSchurDominance));
    assert!(checked.gain_infinity_upper().unwrap()>=2.0);
    assert!(checked.state_error_infinity_upper().unwrap()>=1.0);
}

#[test]
fn one_shared_budget_and_every_cancellation_point_are_respected() {
    let n=3;let a=solid(n);let(b,c)=transfer(n,1.0);let responses=[vec![1.0/6.0;n]];
    let invoke=|cap, cancelled: &mut dyn FnMut()->bool| {
        let mut lim=limits(n);lim.max_verification_entries=cap;
        enclose_affine_feedback_error_with_spectral(
            &a,&[5.;3],&[0.;3],&b,&c,&[0.],Some(&responses),None,lim,spectral(n),cancelled,
        )
    };
    let report=invoke(1_000_000,&mut ||true).unwrap();
    // Both ordinary assessments conservatively charge (n+nnz)*(p+1).
    let need=2*(n+a.nnz())*2+report.solid_spectral().unwrap().work_entries();
    assert_eq!(invoke(need,&mut ||true).unwrap().status(),FeedbackBoundStatus::Enclosed);
    let stopped=invoke(need-1,&mut ||true).unwrap();
    assert_eq!(stopped.solid_spectral().unwrap().stop(), SpectralStop::WorkLimit);
    assert!(stopped.state_error_infinity_upper().is_none());
    let mut polls=0;invoke(need,&mut ||{polls+=1;true}).unwrap();
    for at in 1..=polls {
        let mut seen=0;
        assert!(matches!(invoke(need,&mut ||{seen+=1;seen!=at}),
            Err(SpectralError::Residual(GoalResidualError::Cancelled))));
    }
}

#[test]
fn cached_spectral_authority_keeps_its_owned_operator_and_never_refactors() {
    use fs_solver::goal::feedback::enclose_affine_feedback_error_with_spectral_inverse;
    use fs_solver::goal::inverse::spectral::prepare_spectral_inverse;
    let n=3; let a=solid(n); let (b,c)=transfer(n,1.0); let responses=[vec![1.0/6.0;n]];
    let preparation=prepare_spectral_inverse(&a,1.0,spectral(n),||true).unwrap();
    let certificate=preparation.certificate.unwrap();
    let first=enclose_affine_feedback_error_with_spectral_inverse(
        &certificate,&[5.;3],&[0.;3],&b,&c,&[0.],Some(&responses),None,limits(n),||true,
    ).unwrap();
    assert_eq!(first.status(),FeedbackBoundStatus::Enclosed);
    assert_eq!(first.solid_spectral().unwrap().work_entries(),0);
    assert_eq!(first.solid_spectral().unwrap().shift_attempts(),0);
    let second=enclose_affine_feedback_error_with_spectral_inverse(
        &certificate,&[5.;3],&[0.9;3],&b,&c,&[0.],Some(&responses),None,limits(n),||true,
    ).unwrap();
    assert!(second.state_error_infinity_upper().unwrap()<first.state_error_infinity_upper().unwrap());
    assert_eq!(first.solid_inverse_infinity_upper(),second.solid_inverse_infinity_upper());
}

#[test]
fn existing_comparison_success_does_not_spend_spectral_preparation() {
    let mut a=Coo::new(1,1); a.push(0,0,2.); let a=a.assemble();
    let (b,c)=transfer(1,0.25);
    let ordinary=enclose_affine_feedback_error_with_schur(
        &a,&[1.],&[0.],&b,&c,&[0.],None,None,limits(1),||true,
    ).unwrap();
    let mut disabled=spectral(1); disabled.max_work_entries=0; disabled.max_storage_entries=0;
    let checked=enclose_affine_feedback_error_with_spectral(
        &a,&[1.],&[0.],&b,&c,&[0.],None,None,limits(1),disabled,||true,
    ).unwrap();
    assert_eq!(ordinary,checked);
    assert!(checked.solid_spectral().is_none());
}
