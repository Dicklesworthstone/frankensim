use super::*;
use fs_sparse::Coo;

fn limits() -> SpectralInverseLimits {
    SpectralInverseLimits {
        system: GoalResidualLimits { max_rows: 4096, max_nonzeros: 100_000 },
        max_storage_entries: 200_000, max_work_entries: 4_000_000,
        max_shift_attempts: 12,
    }
}
fn matrix(values: &[&[f64]]) -> Csr {
    let n = values.len();
    let mut a = Coo::new(n, n);
    for (i, row) in values.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() { a.push(i, j, v); }
    }
    a.assemble()
}
fn positive_off_diagonals(n: usize, scale: f64) -> Csr {
    assert_eq!(n % 3, 0);
    let mut a = Coo::new(n, n);
    for i in 0..n {
        for j in (i / 3 * 3)..(i / 3 * 3 + 3) {
            a.push(i, j, scale * if i == j { 2.0 } else { 1.5 });
        }
        // A chain of PSD graph-Laplacian contributions between the blocks.
        if i + 3 < n {
            a.push(i, i, scale * 0.125); a.push(i + 3, i + 3, scale * 0.125);
            a.push(i, i + 3, -scale * 0.125); a.push(i + 3, i, -scale * 0.125);
        }
    }
    a.assemble()
}
fn multiply(a: &Csr, x: &[f64]) -> Vec<f64> {
    (0..a.nrows()).map(|i| {
        a.row(i).0.iter().zip(a.row(i).1).fold(0.0, |sum, (&j, &v)| v.mul_add(x[j], sum))
    }).collect()
}

#[test]
fn sparse_proof_crosses_dense_inverse_cap_on_nondominant_spd_systems() {
    // Every 3x3 block is .5 I + 1.5 11^T; the connecting Laplacian is PSD.
    // Thus lambda_min >= .5, but the comparison matrix cannot be positive.
    let a = positive_off_diagonals(600, 1.0);
    let zeros = vec![0.0; a.nrows()];
    let ordinary = enclose_goal_error(&a, &zeros, &zeros, &zeros, &zeros, None,
        limits().system, || true).unwrap();
    assert!(ordinary.inverse_infinity_upper().is_none());
    let prepared = prepare_spectral_inverse(&a, 0.25, limits(), || true).unwrap();
    assert_eq!(prepared.stop, SpectralStop::Certified);
    assert!(prepared.peak_storage_entries < a.nrows() * a.nrows());
    let proof = prepared.certificate.unwrap();
    assert!(proof.coercivity_lower() > 0.24 && proof.coercivity_lower() <= 0.5);
    assert!(proof.defect_upper() < 1e-10);
    assert_eq!(proof.matrix(), &a);
    let exact = (0..a.nrows()).map(|i| (i % 7) as f64).collect::<Vec<_>>();
    let rhs = multiply(&a, &exact); // This dyadic fixture has exact small sums.
    let approximate = exact.iter().enumerate().map(|(i, &v)| v + if i % 2 == 0 { 0.125 } else { -0.25 }).collect::<Vec<_>>();
    let goal = (0..a.nrows()).map(|i| if i == 7 { 1.0 } else { 0.0 }).collect::<Vec<_>>();
    let report = proof.enclose_goal_error(&rhs, &approximate, &goal, &zeros, None,
        limits().system, || true).unwrap();
    let error = exact[7] - approximate[7];
    let bound = report.goal_error().unwrap();
    assert!(bound.lower() <= error && error <= bound.upper());
    assert!(report.dual_error_upper().unwrap() >= error.abs(), "zero dual is not exact");
}

#[test]
fn shifted_gram_checks_missing_modes_and_unstored_fill() {
    // A good-looking factor on the first mode says nothing about the hidden
    // negative second mode. The full original operator is the proof target.
    let indefinite = matrix(&[&[4.0, 0.0], &[0.0, -0.125]]);
    let hidden = certify_shifted_gram(&indefinite, 0.25, &[vec![(0, 1.9375)]], limits(), || true).unwrap();
    assert_eq!(hidden.stop, SpectralStop::NotEstablished);
    assert!(hidden.certificate.is_none());
    // Off-diagonal Gram products must be checked even when A omits them.
    let a = Csr::from_parts(2, 2, vec![0, 1, 2], vec![0, 1], vec![2.0, 2.0]);
    let wrong = certify_shifted_gram(&a, 1.0, &[vec![(0, 1.0), (1, 1.0)]], limits(), || true).unwrap();
    assert!(wrong.certificate.is_none());
    let right = certify_shifted_gram(&a, 1.0, &[vec![(0, 1.0)], vec![(1, 1.0)]], limits(), || true).unwrap();
    assert!(right.certificate.is_some());
}

#[test]
fn original_asymmetry_is_never_discarded_by_factor_proposal() {
    let a = matrix(&[&[4.0, 1.0 + 1e-8], &[1.0, 4.0]]);
    let prepared = prepare_spectral_inverse(&a, 1.0, limits(), || true).unwrap();
    let proof = prepared.certificate.unwrap();
    assert!(proof.defect_upper() >= 4e-9);
    assert!(proof.coercivity_lower() < 1.0);
    assert_eq!(proof.matrix().row(0).1[1].to_bits(), (1.0_f64 + 1e-8).to_bits());
    let hostile = matrix(&[&[2.0, 0.0], &[5.0, 2.0]]);
    let bogus = certify_shifted_gram(&hostile, 1.0,
        &[vec![(0, 1.0)], vec![(1, 1.0)]], limits(), || true).unwrap();
    assert!(bogus.certificate.is_none(), "the discarded upper/lower mismatch must refuse");
}

#[test]
fn proposal_errors_are_typed_and_do_not_produce_authority() {
    let a = matrix(&[&[2.0, 1.5], &[1.5, 2.0]]);
    for columns in [vec![vec![(0, 1.0), (0, 1.0)]], vec![vec![(1, 1.0), (0, 1.0)]],
        vec![vec![(2, 1.0)]], vec![vec![(0, f64::NAN)]]] {
        assert!(matches!(certify_shifted_gram(&a, 0.25, &columns, limits(), || true),
            Err(SpectralError::InvalidProposal(_))));
    }
    for shift in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(prepare_spectral_inverse(&a, shift, limits(), || true),
            Err(SpectralError::InvalidProposal(_))));
    }
    let mut small = limits(); small.system.max_rows = 1;
    assert!(matches!(prepare_spectral_inverse(&a, 0.25, small, || true),
        Err(SpectralError::Residual(GoalResidualError::Limit { .. }))));
}

#[test]
fn every_cancellation_boundary_prevents_publication() {
    let a = positive_off_diagonals(30, 1.0);
    let mut polls = 0;
    let complete = prepare_spectral_inverse(&a, 0.25, limits(), || { polls += 1; true }).unwrap();
    assert!(complete.certificate.is_some());
    for stop in 1..=polls {
        let mut visited = 0;
        let cancelled = prepare_spectral_inverse(&a, 0.25, limits(), || {
            visited += 1; visited < stop
        });
        assert!(matches!(cancelled,
            Err(SpectralError::Residual(GoalResidualError::Cancelled))), "stop {stop}");
    }
    let proof = complete.certificate.unwrap();
    let zero = vec![0.0; a.nrows()];
    assert!(matches!(proof.enclose_goal_error(&zero, &zero, &zero, &zero, None,
        limits().system, || false), Err(GoalResidualError::Cancelled)));
}

#[test]
fn budgets_cover_failed_shifts_and_exact_success_boundaries() {
    let a = positive_off_diagonals(30, 1.0);
    let baseline = prepare_spectral_inverse(&a, 4.0, limits(), || true).unwrap();
    assert!(baseline.certificate.is_some());
    assert!(baseline.shift_attempts > 1, "retry must be exercised");
    let mut exact = limits();
    exact.max_work_entries = baseline.work_entries;
    exact.max_storage_entries = baseline.peak_storage_entries;
    let repeated = prepare_spectral_inverse(&a, 4.0, exact, || true).unwrap();
    assert!(repeated.certificate.is_some());
    assert_eq!(repeated.work_entries, baseline.work_entries);
    exact.max_work_entries -= 1;
    let exhausted = prepare_spectral_inverse(&a, 4.0, exact, || true).unwrap();
    assert_eq!(exhausted.stop, SpectralStop::WorkLimit);
    assert_eq!(exhausted.work_entries, exact.max_work_entries);
    assert!(exhausted.certificate.is_none());
    exact.max_work_entries = baseline.work_entries;
    exact.max_storage_entries -= 1;
    let exhausted = prepare_spectral_inverse(&a, 4.0, exact, || true).unwrap();
    assert_eq!(exhausted.stop, SpectralStop::StorageLimit);
    assert!(exhausted.certificate.is_none());
}

#[test]
fn rescaling_and_replay_preserve_verified_operator_binding() {
    for power in [-400, 0, 400] {
        let scale = 2.0_f64.powi(power);
        let a = positive_off_diagonals(12, scale);
        let left = prepare_spectral_inverse(&a, 0.25 * scale, limits(), || true).unwrap();
        let right = prepare_spectral_inverse(&a, 0.25 * scale, limits(), || true).unwrap();
        assert_eq!(left.work_entries, right.work_entries);
        let l = left.certificate.unwrap(); let r = right.certificate.unwrap();
        assert_eq!(l.coercivity_lower().to_bits(), r.coercivity_lower().to_bits());
        assert_eq!(l.inverse_infinity_upper().to_bits(), r.inverse_infinity_upper().to_bits());
        assert!(l.coercivity_lower() <= 0.5 * scale);
        assert!(l.coercivity_lower() > 0.24 * scale);
        // A certificate owns its exact coefficients, not an ambient pointer.
        let snapshot = l.clone();
        assert_eq!(snapshot.matrix(), &a);
    }
}

#[test]
fn singular_and_indefinite_operators_do_not_gain_a_spectral_inverse() {
    for a in [matrix(&[&[1.0, -1.0], &[-1.0, 1.0]]),
        matrix(&[&[1.0, 2.0], &[2.0, 1.0]])] {
        let result = prepare_spectral_inverse(&a, 0.25, limits(), || true).unwrap();
        assert_eq!(result.stop, SpectralStop::NotEstablished);
        assert!(result.certificate.is_none());
    }
}
