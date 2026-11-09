//! G1 independent original coefficient/field product rules; G3 normalized
//! physical units and gyro separation; G0 complete-jet refusals and replay.
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{HarmonicPowerParameterJet, HarmonicTransferFieldParameterJet};
use fs_modal::harmonic::{
    HarmonicResponseError, harmonic_quadratic_power_curvature_with_cancel as curvature,
    harmonic_quadratic_power_parameter_jet_with_cancel as observe,
};

fn field<const N: usize>(orders: &[[C64; N]; 5]) -> HarmonicTransferFieldParameterJet<'_> {
    HarmonicTransferFieldParameterJet {
        values: &orders[0],
        angular_derivative: &orders[1],
        angular_second_derivative: &orders[2],
        parameter_derivative: &orders[3],
        angular_parameter_derivative: &orders[4],
    }
}
fn values(jet: &HarmonicPowerParameterJet) -> [f64; 5] {
    [
        jet.sample.average_power_w,
        jet.sample.angular_derivative,
        jet.angular_second_derivative,
        jet.parameter_derivative,
        jet.angular_parameter_derivative,
    ]
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-15 * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn complex_material_and_field_changes_retain_every_mixed_work_term() {
    // D=[[2,1],[1,5]], D_theta=[[3,-2],[-2,0.5]]. Original operators
    // additionally contain independent skew parts of strengths7 and11.
    let orders = [
        [C64::new(1., 2.), C64::new(2., -1.)],
        [C64::new(3., -1.), C64::new(-1., 2.)],
        [C64::new(-2., 1.), C64::new(4., 3.)],
        [C64::new(0.5, -1.), C64::new(-2., 0.75)],
        [C64::new(1.5, 0.25), C64::new(0.25, -0.5)],
    ];
    let matrix = [2., 8., -6., 5.];
    let direction = [3., 9., -13., 0.5];
    let result = observe(field(&orders), &matrix, Some(&direction), || false).unwrap();
    let symmetric = observe(
        field(&orders),
        &[2., 1., 1., 5.],
        Some(&[3., -2., -2., 0.5]),
        || false,
    )
    .unwrap();
    assert_eq!(result, symmetric);
    // Independent real/imaginary polynomial expansion gives the first3
    // orders. Theta is -25.25+8.75; mixed is 11+13.25-19.
    for (actual, expected) in values(&result)
        .into_iter()
        .zip([17.5, -8., 65., -16.5, 5.25])
    {
        near(actual, expected);
    }
    let fixed = observe(field(&orders), &matrix, None, || false).unwrap();
    near(fixed.parameter_derivative, -25.25);
    near(fixed.angular_parameter_derivative, 24.25);
    let ordinary = curvature(&orders[0], &orders[1], &orders[2], &matrix, || false).unwrap();
    assert_eq!(ordinary.sample, result.sample);
    assert_eq!(
        ordinary.angular_second_derivative,
        result.angular_second_derivative
    );
}

#[test]
fn coefficient_direction_retains_weak_work_subnormal_pairs_and_finite_extreme_products() {
    let weak = 2.0_f64.powi(-900);
    let gyro = 2.0_f64.powi(1000);
    let amplitude = 2.0_f64.powi(500);
    let orders = [
        [C64::from_re(amplitude); 2],
        [C64::from_re(3. * amplitude); 2],
        [C64::from_re(-2. * amplitude); 2],
        [C64::from_re(4. * amplitude); 2],
        [C64::from_re(5. * amplitude); 2],
    ];
    let matrix = [weak, gyro, -gyro, 2. * weak];
    let direction = [-weak, -gyro, gyro, 3. * weak];
    let result = observe(field(&orders), &matrix, Some(&direction), || false).unwrap();
    let unit = 2.0_f64.powi(100);
    for (actual, multiplier) in values(&result).into_iter().zip([1.5, 9., 21., 13., 57.]) {
        near(actual, multiplier * unit);
    }

    // A pure gyro primal still has an independently authored symmetric
    // material direction. Halving its minsubnormal coefficient would erase
    // the finite response. Changing fields here must not add a C_theta*v_theta.
    let pure_gyro = [0., gyro, -gyro, 0.];
    let tiny = f64::from_bits(1);
    for (a, b, multiplier) in [(tiny, 0., 0.5), (tiny, 2. * tiny, 1.5)] {
        let result = observe(field(&orders), &pure_gyro, Some(&[0., a, b, 0.]), || false).unwrap();
        let expected = multiplier * 2.0_f64.powi(-74);
        assert_eq!(&values(&result)[..3], &[0.; 3]);
        near(result.parameter_derivative, expected);
        near(result.angular_parameter_derivative, 6. * expected);
    }

    // C_theta01+C_theta10 overflows, but its field-weighted work is finite.
    let coefficient = 2.0_f64.powi(1023);
    let small = [C64::from_re(2.0_f64.powi(-510)); 2];
    let zero = [C64::ZERO; 2];
    let orders = [small, small, zero, zero, zero];
    let result = observe(
        field(&orders),
        &[0.; 4],
        Some(&[0., coefficient, coefficient, 0.]),
        || false,
    )
    .unwrap();
    near(result.parameter_derivative, 8.);
    near(result.angular_parameter_derivative, 16.);

    // Preserve the coefficient sum's low word through leading-work cancellation.
    let v = [C64::ONE, -C64::ONE];
    let orders = [v, v, zero, zero, zero];
    let low = 2.0_f64.powi(-60);
    let result = observe(field(&orders), &[0.; 4], Some(&[1., 2., low, 1.]), || false).unwrap();
    near(result.parameter_derivative, -0.5 * low);
    near(result.angular_parameter_derivative, -low);
}

#[test]
fn malformed_direction_and_unrepresentable_or_cancelled_jets_publish_nothing() {
    let values = [C64::ONE, C64::new(0.5, 0.25)];
    let zero = [C64::ZERO; 2];
    let orders = [values, zero, zero, zero, zero];
    let matrix = [2., 0.5, 0.5, 3.];
    let direction = [0.1, 0., 0., -0.2];
    let baseline = observe(field(&orders), &matrix, Some(&direction), || false).unwrap();
    assert!(observe(field(&orders), &matrix, Some(&direction[..3]), || false).is_err());
    assert!(
        observe(
            field(&orders),
            &matrix,
            Some(&[0., f64::INFINITY, 0., 0.]),
            || false
        )
        .is_err()
    );
    assert!(observe(field(&orders), &matrix[..3], Some(&direction), || false).is_err());
    let mut bad = field(&orders);
    bad.angular_parameter_derivative = &zero[..1];
    assert!(observe(bad, &matrix, None, || false).is_err());
    let mut nonfinite = orders;
    nonfinite[4][1] = C64::new(f64::NAN, 0.);
    assert!(observe(field(&nonfinite), &matrix, None, || false).is_err());
    let empty = [[C64::ZERO; 0]; 5];
    assert!(observe(field(&empty), &[], None, || false).is_err());
    let mut polls = 0;
    assert!(matches!(
        observe(field(&orders), &matrix, Some(&direction), || {
            polls += 1;
            polls == 3
        }),
        Err(HarmonicResponseError::Cancelled)
    ));
    assert_eq!(polls, 3);
    assert_eq!(
        baseline,
        observe(field(&orders), &matrix, Some(&direction), || false).unwrap()
    );

    // The unchanged matrix has zero work, but an unrepresentable physical
    // parameter direction refuses the complete result instead of returning
    // the otherwise finite frequency sub-jet.
    let huge = [
        [C64::from_re(2.0_f64.powi(900))],
        [C64::ZERO],
        [C64::ZERO],
        [C64::ZERO],
        [C64::ZERO],
    ];
    assert!(observe(field(&huge), &[0.], Some(&[1.]), || false).is_err());
}
