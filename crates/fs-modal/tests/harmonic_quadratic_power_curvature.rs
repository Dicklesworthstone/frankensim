//! G1 original quadratic work and gyroscopic separation; G3 coefficient/field
//! scaling; G0 bounded shape, arithmetic and cancellation refusals.
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicResponseError, harmonic_quadratic_power_curvature_with_cancel as observe,
};

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-15 * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn original_complex_work_curvature_ignores_gyroscopic_exchange() {
    // D=[[2,1],[1,5]], C=D+[[0,7],[-7,0]]. Independent real/imaginary
    // quadratic polynomials contribute (P,P',P'')=(13,1,53)+(4.5,-9,12).
    let v = [C64::new(1., 2.), C64::new(2., -1.)];
    let first = [C64::new(3., -1.), C64::new(-1., 2.)];
    let second = [C64::new(-2., 1.), C64::new(4., 3.)];
    let with_skew = observe(&v, &first, &second, &[2., 8., -6., 5.], || false).unwrap();
    let symmetric = observe(&v, &first, &second, &[2., 1., 1., 5.], || false).unwrap();
    assert_eq!(with_skew, symmetric);
    near(with_skew.sample.average_power_w, 17.5);
    near(with_skew.sample.angular_derivative, -8.);
    near(with_skew.angular_second_derivative, 65.);

    // The enormous skew work would exceed f64 and DD's relative precision
    // beside the weak diagonal dissipation if accumulated before cancellation.
    let weak = 2.0_f64.powi(-900);
    let gyro = 2.0_f64.powi(1000);
    let field = 2.0_f64.powi(500);
    let v = [C64::from_re(field); 2];
    let first = [C64::from_re(3. * field); 2];
    let second = [C64::from_re(-2. * field); 2];
    let physical = observe(&v, &first, &second, &[weak, gyro, -gyro, 2. * weak], || {
        false
    })
    .unwrap();
    let base = 1.5 * 2.0_f64.powi(100);
    near(physical.sample.average_power_w, base);
    near(physical.sample.angular_derivative, 6. * base);
    near(physical.angular_second_derivative, 14. * base);
}

#[test]
fn original_coefficients_survive_subnormal_halves_and_overflowing_pair_sums() {
    let tiny = f64::from_bits(1);
    let field = 2.0_f64.powi(500);
    let v = [C64::from_re(field); 2];
    let zero = [C64::ZERO; 2];
    for (a, b, multiplier) in [(tiny, 0., 0.5), (tiny, 2. * tiny, 1.5)] {
        // Neither half of minsubnormal nor the average of the unequal pair
        // is representable, but their original quadratic work is normal.
        let result = observe(&v, &v, &zero, &[0., a, b, 0.], || false).unwrap();
        let expected = multiplier * 2.0_f64.powi(-74);
        near(result.sample.average_power_w, expected);
        near(result.sample.angular_derivative, 2. * expected);
        near(result.angular_second_derivative, 2. * expected);
    }

    let coefficient = 2.0_f64.powi(1023);
    assert!((coefficient + coefficient).is_infinite());
    let small = [C64::from_re(2.0_f64.powi(-510)); 2];
    let result = observe(
        &small,
        &small,
        &zero,
        &[0., coefficient, coefficient, 0.],
        || false,
    )
    .unwrap();
    near(result.sample.average_power_w, 8.);
    near(result.sample.angular_derivative, 16.);
    near(result.angular_second_derivative, 16.);

    // Signed coefficient directions are legitimate observations. The low word
    // of C01+C10 remains after the representable leading work cancels.
    let v = [C64::ONE, -C64::ONE];
    let low = 2.0_f64.powi(-60);
    let result = observe(&v, &v, &zero, &[1., 2., low, 1.], || false).unwrap();
    near(result.sample.average_power_w, -0.5 * low);
    near(result.sample.angular_derivative, -low);
    near(result.angular_second_derivative, -low);
}

#[test]
fn malformed_nonfinite_or_cancelled_work_publishes_no_partial_result() {
    let v = [C64::ONE, C64::new(0.5, 0.25)];
    let zero = [C64::ZERO; 2];
    let matrix = [2., 0.5, 0.5, 3.];
    let baseline = observe(&v, &zero, &zero, &matrix, || false).unwrap();
    assert!(observe(&[], &[], &[], &[], || false).is_err());
    assert!(observe(&v, &zero[..1], &zero, &matrix, || false).is_err());
    assert!(observe(&v, &zero, &zero, &matrix[..3], || false).is_err());
    assert!(observe(&v, &zero, &zero, &[2., f64::INFINITY, 0.5, 3.], || false).is_err());
    assert!(
        observe(
            &v,
            &[C64::new(f64::NAN, 0.), C64::ZERO],
            &zero,
            &matrix,
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        observe(&v, &zero, &zero, &matrix, || {
            polls += 1;
            polls == 3
        }),
        Err(HarmonicResponseError::Cancelled)
    ));
    assert_eq!(polls, 3);
    assert_eq!(
        baseline,
        observe(&v, &zero, &zero, &matrix, || false).unwrap()
    );
    let huge = [C64::from_re(2.0_f64.powi(900))];
    assert!(observe(&huge, &[C64::ZERO], &[C64::ZERO], &[1.], || false).is_err());
}
