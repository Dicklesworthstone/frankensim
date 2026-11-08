use super::*;

#[test]
fn velocity_curvature_cancels_before_doubling_overflows() {
    let amplitude = 0.75 * f64::MAX;
    let first = C64::new(amplitude, -amplitude);
    assert!((2.0 * amplitude).is_infinite());
    assert_eq!(velocity_curvature(2.0, first, -first).unwrap(), C64::ZERO);
}

#[test]
fn lossless_source_phase_curvature_has_zero_work_despite_overflowing_products() {
    // M=1, K=2, C=0, omega=1 and F(omega)=F0*exp(i*omega).
    // The complete analytic jets are binary-exact multiples at this point.
    let amplitude = 2.0_f64.powi(600);
    let force = C64::new(amplitude, amplitude);
    let force_first = C64::new(-amplitude, amplitude);
    let force_second = -force;
    let q = force;
    let q_first = C64::new(2.0, 1.0) * force;
    let q_second = C64::new(9.0, 4.0) * force;
    let velocity = C64::new(q.im, -q.re);
    let first_sum = q + q_first;
    let velocity_first = C64::new(first_sum.im, -first_sum.re);
    let velocity_second = velocity_curvature(1.0, q_first, q_second).unwrap();
    assert!((force_first.re * velocity_first.re).is_infinite());
    assert!((force.re * velocity_second.re).is_infinite());
    let power = harmonic_mean_power_curvature_with_cancel(
        &[velocity],
        &[velocity_first],
        &[velocity_second],
        &[force],
        &[force_first],
        &[force_second],
        || false,
    )
    .unwrap();
    assert_eq!(power, 0.0);
}

#[test]
fn finite_high_q_resonance_work_retains_all_viscous_cross_terms() {
    // M=K=omega=1, constant F=2^350 and C=2^-100.
    // P''=-4*M^2*|F|^2/C^3=-2^1002 is finite. Squaring v' first is not.
    let force = C64::from_re(2.0_f64.powi(350));
    let damping = 2.0_f64.powi(-100);
    let velocity = C64::from_re(2.0_f64.powi(450));
    let velocity_first = C64::new(0.0, 2.0_f64.powi(551));
    let velocity_second = C64::new(-2.0_f64.powi(653), -2.0_f64.powi(551));
    assert!((velocity_first.im * velocity_first.im).is_infinite());
    let source = harmonic_mean_power_curvature_with_cancel(
        &[velocity],
        &[velocity_first],
        &[velocity_second],
        &[force],
        &[C64::ZERO],
        &[C64::ZERO],
        || false,
    )
    .unwrap();
    let mut viscous = PowerSum::default();
    viscous
        .half_inner(velocity_second, velocity, damping)
        .unwrap();
    viscous
        .half_inner(velocity_first, velocity_first, damping)
        .unwrap();
    viscous
        .half_inner(velocity_first, velocity_first, damping)
        .unwrap();
    viscous
        .half_inner(velocity, velocity_second, damping)
        .unwrap();
    assert_eq!(source, -2.0_f64.powi(1002));
    assert_eq!(viscous.value().unwrap(), source);
}

#[test]
fn cross_term_preserves_minimum_subnormal_work() {
    let tiny = f64::from_bits(1);
    let power = harmonic_mean_power_curvature_with_cancel(
        &[C64::ZERO],
        &[C64::ONE],
        &[C64::ZERO],
        &[C64::ZERO],
        &[C64::from_re(tiny)],
        &[C64::ZERO],
        || false,
    )
    .unwrap();
    assert_eq!(power.to_bits(), tiny.to_bits());
}

#[test]
fn work_refusal_and_final_cancellation_publish_no_partial_sum() {
    let z = [C64::ZERO];
    let mut calls = 0;
    let result = harmonic_mean_power_curvature_with_cancel(&z, &z, &z, &z, &z, &z, || {
        calls += 1;
        calls == 3
    });
    assert!(matches!(result, Err(HarmonicResponseError::Cancelled)));
    assert_eq!(calls, 3);
    assert!(
        harmonic_mean_power_curvature_with_cancel(
            &z,
            &z,
            &z,
            &z,
            &z,
            &[C64::new(f64::NAN, 0.0)],
            || false,
        )
        .is_err()
    );
    assert!(
        harmonic_mean_power_curvature_with_cancel(&[], &[], &[], &[], &[], &[], || false,).is_err()
    );
}
