//! G0/G1/G4: analytic history rates use the original fallible recurrence.

use fs_cheb::{Cheb1, ChebError};

#[test]
fn checked_derivative_retains_native_coefficients_and_physical_chain_rule() {
    // On [2, 6], x = 4 + 2t. Coefficients describe x^3 exactly:
    // x^3 = 88 + 102 T1(t) + 24 T2(t) + 2 T3(t).
    let cubic = Cheb1::from_coeffs(2.0, 6.0, vec![176.0, 102.0, 24.0, 2.0]);
    let checked = cubic
        .try_differentiate_with_checkpoint(&mut || Ok(()))
        .unwrap();
    assert_eq!(checked.coeffs(), &[108.0, 48.0, 6.0]);
    assert_eq!(checked, cubic.differentiate());
    for x in [2.0, 2.5, 4.0, 5.5, 6.0] {
        assert!((checked.eval(x) - 3.0 * x * x).abs() < 1e-12);
    }

    // The reference derivative overflows but the physical derivative does not.
    let wide = Cheb1::from_coeffs(-1e308, 1e308, vec![0.0, 1e308]);
    let derivative = wide
        .try_differentiate_with_checkpoint(&mut || Ok(()))
        .unwrap();
    assert_eq!(derivative.coeffs(), &[2.0]);
    assert_eq!(derivative.eval(0.0), 1.0);
    assert_eq!(derivative, wide.differentiate());
}

#[test]
fn finite_source_with_unrepresentable_rate_returns_a_numerical_refusal() {
    let tiny = f64::from_bits(1);
    let source = Cheb1::from_coeffs(0.0, tiny, vec![0.0, 1.0]);
    assert!(matches!(
        source.try_differentiate_with_checkpoint(&mut || Ok(())),
        Err(ChebError::Numerical { .. })
    ));
}

#[test]
fn derivative_cancels_mid_recurrence_without_a_partial_object() {
    let source = Cheb1::from_coeffs(-2.0, 3.0, vec![2.0, -0.5, 0.25, -0.125, 0.0625]);
    let mut polls = 0;
    let result = source.try_differentiate_with_checkpoint(&mut || {
        polls += 1;
        if polls == 3 {
            Err(ChebError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert_eq!(result, Err(ChebError::Cancelled));
    assert_eq!(polls, 3);
    assert_eq!(source.coeffs(), &[2.0, -0.5, 0.25, -0.125, 0.0625]);
}
