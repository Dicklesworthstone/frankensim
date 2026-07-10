//! Battery for proof-carrying optimization (fs-sos). The load-bearing property
//! is ZERO FALSE CERTIFICATES: a sum-of-squares certificate is verified by a
//! polynomial identity, so a claimed bound above the true minimum fails, while a
//! verified bound is sound (holds for every sampled point).

use fs_sos::{
    Poly, SosCertificate, certify_quadratic, is_psd, lyapunov_certifies_stability, square,
};

#[test]
fn polynomial_arithmetic_is_correct() {
    // (1+x)(1-x) = 1 - x².
    let prod = Poly::new(vec![1.0, 1.0]).mul(&Poly::new(vec![1.0, -1.0]));
    assert_eq!(prod.coeffs(), &[1.0, 0.0, -1.0]);
    assert!((prod.eval(2.0) - (-3.0)).abs() < 1e-12);
    assert_eq!(prod.degree(), 2);
    // (x-2)² = x² - 4x + 4.
    assert_eq!(
        square(&Poly::new(vec![-2.0, 1.0])).coeffs(),
        &[4.0, -4.0, 1.0]
    );
}

#[test]
fn a_quadratic_certificate_is_exact_and_sound() {
    // p(x) = x² - 4x + 7, global min 3 at x = 2.
    let p = Poly::new(vec![7.0, -4.0, 1.0]);
    let cert = certify_quadratic(1.0, -4.0, 7.0).unwrap();
    assert!(cert.verify(&p, 1e-9));
    // Radius-scoped value claim (bead wa8i S2: the old tol-based global
    // claim was forgeable): sound and barely degraded for an honest cert.
    let bound = cert.certified_bound_on(&p, 1e3).unwrap();
    assert!(bound <= 3.0 + 1e-12 && bound > 3.0 - 1e-6, "bound {bound}");
    // the certified bound is a true lower bound everywhere.
    for k in -50..=50 {
        let x = f64::from(k) * 0.2;
        assert!(p.eval(x) >= 3.0 - 1e-9, "p({x}) = {} < 3", p.eval(x));
    }
}

#[test]
fn a_multi_square_certificate_proves_a_tight_bound() {
    // p(x) = x⁴ + 2x² + 1, min 1 at x = 0; p - 1 = (x²)² + (√2 x)².
    let p = Poly::new(vec![1.0, 0.0, 2.0, 0.0, 1.0]);
    let cert = SosCertificate {
        squares: vec![
            Poly::new(vec![0.0, 0.0, 1.0]),
            Poly::new(vec![0.0, 2.0_f64.sqrt()]),
        ],
        lower_bound: 1.0,
    };
    assert!(cert.verify(&p, 1e-9));
    let bound = cert.certified_bound_on(&p, 1e3).unwrap();
    assert!(bound <= 1.0 + 1e-12 && bound > 1.0 - 1e-6, "bound {bound}");
    for k in -50..=50 {
        let x = f64::from(k) * 0.2;
        assert!(p.eval(x) >= 1.0 - 1e-9);
    }
}

#[test]
fn there_are_zero_false_certificates() {
    let p = Poly::new(vec![7.0, -4.0, 1.0]); // min 3
    // a certificate claiming a bound ABOVE the true minimum does not verify.
    let liar = SosCertificate {
        squares: vec![Poly::new(vec![-2.0, 1.0])], // the correct square (x-2)
        lower_bound: 5.0,                          // but a false, too-high bound
    };
    assert!(!liar.verify(&p, 1e-9));
    // The liar cannot extract a false value claim from ANY door — and
    // the global door does something better than refuse: this liar's
    // residual is EXACTLY the constant −2 (dyadic arithmetic), so the
    // exact absorption CORRECTS the lie: p − 5 = (x−2)² − 2 ⇒ p ≥ 3,
    // the TRUE minimum. The claimed 5.0 is never repeated.
    let corrected = liar.certified_bound_global(&p).unwrap();
    assert!(
        corrected < liar.lower_bound && corrected <= 3.0 + 1e-12 && corrected > 3.0 - 1e-9,
        "the lie is corrected to the true bound, not repeated: {corrected}"
    );
    let liar_bound = liar.certified_bound_on(&p, 1e3).unwrap();
    assert!(
        liar_bound <= 3.0 + 1e-12,
        "radius door also never exceeds truth: {liar_bound}"
    );
    // a bogus square set does not verify either.
    let bogus = SosCertificate {
        squares: vec![Poly::new(vec![0.0, 1.0])], // x², missing the +1
        lower_bound: 0.0,
    };
    assert!(!bogus.verify(&Poly::new(vec![1.0, 0.0, 1.0]), 1e-9)); // p = x² + 1
}

#[test]
fn certify_quadratic_rejects_unbounded_forms() {
    assert!(certify_quadratic(-1.0, 0.0, 0.0).is_none()); // opens downward
    assert!(certify_quadratic(0.0, 1.0, 0.0).is_none()); // linear
    assert!(certify_quadratic(f64::NAN, 0.0, 0.0).is_none());
    assert!(certify_quadratic(1.0, f64::INFINITY, 0.0).is_none());
    assert!(certify_quadratic(f64::MIN_POSITIVE, f64::MAX, 0.0).is_none());
}

#[test]
fn psd_feasibility_is_decided_correctly() {
    assert!(is_psd(&[vec![2.0, 0.0], vec![0.0, 3.0]], 1e-9));
    assert!(is_psd(&[vec![2.0, 1.0], vec![1.0, 2.0]], 1e-9)); // eigenvalues 1, 3
    assert!(is_psd(&[vec![0.0, 0.0], vec![0.0, 0.0]], 1e-9)); // boundary
    assert!(!is_psd(&[vec![1.0, 2.0], vec![2.0, 1.0]], 1e-9)); // eigenvalues 3, -1
    assert!(!is_psd(&[vec![1.0, 0.0], vec![0.0, -0.5]], 1e-9));
}

#[test]
fn a_lyapunov_certificate_decides_linear_stability() {
    let identity = [[1.0, 0.0], [0.0, 1.0]];
    // stable: ẋ = diag(-1,-2) x, V = xᵀx certifies stability.
    assert!(lyapunov_certifies_stability(
        [[-1.0, 0.0], [0.0, -2.0]],
        identity
    ));
    // unstable: an eigenvalue at +1 -> no quadratic certificate with P = I.
    assert!(!lyapunov_certifies_stability(
        [[1.0, 0.0], [0.0, -1.0]],
        identity
    ));
}

#[test]
fn non_symmetric_input_cannot_forge_a_certificate() {
    // A = [[0,-1],[-1,0]] has eigenvalues ±1 → UNSTABLE. A quadratic form
    // xᵀPx depends only on the SYMMETRIC part (P+Pᵀ)/2, so a non-symmetric P
    // must never forge a certificate the symmetric part cannot support. Here
    // (P+Pᵀ)/2 = [[1, 1.5],[1.5, 1]] is INDEFINITE (eigenvalues 2.5, −0.5).
    let a = [[0.0, -1.0], [-1.0, 0.0]];
    let p = [[1.0, 0.0], [3.0, 1.0]];
    assert!(
        !lyapunov_certifies_stability(a, p),
        "forged a stability certificate for an UNSTABLE system via non-symmetric P"
    );
    // is_psd must judge the quadratic form (its symmetric part), which is
    // indefinite → NOT PSD.
    assert!(
        !is_psd(&[vec![1.0, 0.0], vec![3.0, 1.0]], 1e-9),
        "non-symmetric matrix falsely judged PSD (symmetric part is indefinite)"
    );
    // A genuinely PSD form given non-symmetrically must still read PSD
    // (symmetric part [[2,1],[1,2]], eigenvalues 1,3).
    assert!(is_psd(&[vec![2.0, 0.0], vec![2.0, 2.0]], 1e-9));
}

#[test]
fn certification_is_deterministic() {
    let p = Poly::new(vec![7.0, -4.0, 1.0]);
    let a = certify_quadratic(1.0, -4.0, 7.0).unwrap();
    let b = certify_quadratic(1.0, -4.0, 7.0).unwrap();
    assert_eq!(a.residual(&p).to_bits(), b.residual(&p).to_bits());
}

#[test]
fn s2_counterexample_cannot_forge_a_value_claim() {
    // Bead wa8i S2: squares = [x − 10⁶], lower_bound = 0,
    // p = (x − 10⁶)² − 1e−9 − 1e−9·x. The coefficient residual is
    // ~9.3e−10 ≤ 1e−9, so the IDENTITY check passes — but
    // p(10⁶) ≈ −9.8e−4 < 0, so a global "p ≥ 0" claim is a FORGERY.
    let q = Poly::new(vec![-1e6, 1.0]);
    let sq = fs_sos::square(&q);
    let c = sq.coeffs().to_vec();
    let p = Poly::new(vec![c[0] - 1e-9, c[1] - 1e-9, c[2]]);
    let cert = SosCertificate {
        squares: vec![q],
        lower_bound: 0.0,
    };
    // The identity check is honest about what it checks…
    assert!(
        cert.verify(&p, 1e-9),
        "the decomposition identity DOES hold to tol"
    );
    // …but no value claim is available: global requires exactly-zero
    // non-constant residual coefficients.
    assert_eq!(
        cert.certified_bound_global(&p),
        None,
        "a nonzero linear residual coefficient must block the global claim"
    );
    // The radius-scoped bound is a TRUE statement that correctly does
    // not certify p ≥ 0: at radius 2e6 the envelope exceeds |p(10⁶)|.
    let bound = cert.certified_bound_on(&p, 2e6).expect("finite radius");
    let p_at_root: f64 = p
        .coeffs()
        .iter()
        .rev()
        .fold(0.0, |acc, &ck| acc.mul_add(1e6, ck));
    assert!(
        p_at_root < 0.0,
        "the counterexample dips negative: {p_at_root:.3e}"
    );
    assert!(
        bound <= p_at_root,
        "the rigorous bound must sit at or below the true dip: {bound:.3e} vs {p_at_root:.3e}"
    );
}

#[test]
fn radius_bound_is_sound_and_useful_for_honest_certs() {
    // certify_quadratic's own certificates: the float residual is a few
    // ulps, so the radius bound degrades lower_bound only negligibly —
    // and stays SOUND (≤ the true minimum + nothing).
    let cert = certify_quadratic(1.0, -4.0, 7.0).expect("a > 0"); // min 3 at x=2
    let p = Poly::new(vec![7.0, -4.0, 1.0]);
    let bound = cert.certified_bound_on(&p, 1e3).expect("finite radius");
    assert!(
        bound <= 3.0 + 1e-12,
        "never above the true minimum: {bound}"
    );
    assert!(
        bound > 3.0 - 1e-6,
        "an honest cert must not be degraded into uselessness: {bound}"
    );
    // A LIAR certificate (claims a bound above the true minimum) stays
    // structurally unable to forge: its residual is large, so the
    // radius bound lands far BELOW its claim (the lie becomes visible
    // as degradation, never as a false claim).
    let liar = SosCertificate {
        squares: vec![Poly::new(vec![-2.0, 1.0])],
        lower_bound: 5.0,
    };
    let liar_bound = liar.certified_bound_on(&p, 1e3).expect("finite radius");
    assert!(
        liar_bound < 3.0,
        "the liar's rigorous bound must fall below the true minimum: {liar_bound}"
    );
}

#[test]
fn exact_zero_residual_certifies_globally() {
    // A certificate built from EXACT dyadic coefficients: p = (x − 1)²
    // + 2 with q = x − 1 (all coefficients and products exactly
    // representable) ⇒ the residual is identically zero and the GLOBAL
    // claim is a theorem.
    let q = Poly::new(vec![-1.0, 1.0]);
    let p = Poly::new(vec![3.0, -2.0, 1.0]); // (x−1)² + 2
    let cert = SosCertificate {
        squares: vec![q],
        lower_bound: 2.0,
    };
    let bound = cert
        .certified_bound_global(&p)
        .expect("exact dyadic cert must certify globally");
    assert!(
        (bound - 2.0).abs() < 1e-12,
        "the global bound is the exact lower bound: {bound}"
    );
}
