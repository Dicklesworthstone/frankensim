use super::*;

fn close(actual: C64, expected: C64, relative: f64) {
    let scale = component_scale(actual).max(component_scale(expected));
    assert!(
        component_scale(actual - expected) <= relative * scale,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn high_q_frequency_slope_reuses_the_physical_factor_with_a_complex_drive_slope() {
    let gate = CancelGate::new_clock_free();
    let natural = core::f64::consts::TAU * 440.0;
    let damping = 2.0e-10 * natural;
    let force = C64::new(0.3, -0.2);
    let drive_slope = C64::new(0.0007, 0.0003);
    for omega in [
        natural,
        f64::from_bits(natural.to_bits() + 37),
        natural * 1.000_001,
    ] {
        let dynamic = C64::new(frequency_detuning(natural, omega), -omega * damping);
        let matrix = [dynamic];
        let factor = HarmonicFactor::new(&matrix, 1, &gate).unwrap();
        let (state, residual) = factor.solve(&[force], 1e-12, &gate).unwrap();
        let rhs = drive_slope - pole_derivative(omega, damping) * state[0];
        let (tangent, derivative_residual) = factor.solve(&[rhs], 1e-12, &gate).unwrap();
        // Independent quotient rule, with no differencing of a narrow peak.
        let expected =
            drive_slope / dynamic + (force / dynamic) * (C64::new(2.0 * omega, damping) / dynamic);
        close(tangent[0], expected, 2e-14);
        close(
            velocity_derivative(omega, state[0], tangent[0]),
            C64::new(0.0, -1.0) * (state[0] + tangent[0].scale(omega)),
            2e-14,
        );
        assert!(residual <= 1e-12 && derivative_residual <= 1e-12);
    }
}

#[test]
fn causal_maxwell_gain_tangent_keeps_small_and_large_relaxation_scales() {
    // At omega=r the analytic derivative is real, 1/(2*r). A squared
    // denominator would overflow at the large scale and underflow at the small.
    for rate in [1e-150, 2.0, 1e200] {
        close(
            relaxation_gain_derivative(rate, rate),
            C64::from_re(0.5 / rate),
            3e-14,
        );
    }
    let rate: f64 = 700.0;
    for omega in [3.0, 300.0, 6000.0] {
        let square = omega * omega + rate * rate;
        let expected = C64::new(
            2.0 * omega * rate * rate / square.powi(2),
            rate * (omega * omega - rate * rate) / square.powi(2),
        );
        close(relaxation_gain_derivative(omega, rate), expected, 3e-14);
        let denominator = C64::new(rate, -omega);
        let strain = C64::new(0.003, -0.002);
        let strain_derivative = C64::new(-0.0002, 0.0007);
        let gain = C64::new(0.0, -omega) / denominator;
        let memory_derivative = (C64::from_re(rate) / denominator) * strain_derivative
            - relaxation_gain_derivative(omega, rate) * strain;
        let mismatch_derivative =
            gain * strain_derivative + relaxation_gain_derivative(omega, rate) * strain;
        close(
            memory_derivative + mismatch_derivative,
            strain_derivative,
            3e-14,
        );
    }
}

#[test]
fn signed_rank_tangent_preserves_original_support_and_all_cross_terms() {
    let column = [0.5, 0.0, -1.25];
    let state = [C64::new(0.2, -0.3), C64::new(1.0, 2.0), C64::new(-0.1, 0.4)];
    for gain in [
        C64::new(0.0, -2.7),
        relaxation_gain_derivative(650.0, 500.0),
    ] {
        let mut rhs = [C64::ZERO; 4];
        subtract_rank_action(&mut rhs, &column, &state, gain);
        for i in 0..3 {
            let expected = (0..3).fold(C64::ZERO, |sum, j| {
                sum - gain.scale(column[i] * column[j]) * state[j]
            });
            if expected == C64::ZERO {
                assert_eq!(rhs[i], C64::ZERO);
            } else {
                close(rhs[i], expected, 3e-14);
            }
        }
        assert_eq!(rhs[3], C64::ZERO);
    }
}

#[test]
fn reciprocal_air_antiresonance_has_a_finite_nonzero_displacement_slope() {
    let gate = CancelGate::new_clock_free();
    let omega: f64 = 700.0;
    let coupling = 100.0;
    let source = C64::new(0.3, 0.1);
    let a = C64::new(frequency_detuning(500.0, omega), -omega * 20.0);
    let b = C64::new(0.0, -omega * coupling);
    let matrix = [a, b, -b, C64::ZERO];
    let factor = HarmonicFactor::new(&matrix, 2, &gate).unwrap();
    let (state, _) = factor.solve(&[source], 1e-12, &gate).unwrap();
    assert!(state[0].abs() < 1e-20);
    let mut rhs = [
        -pole_derivative(omega, 20.0) * state[0],
        -pole_derivative(omega, 0.0) * state[1],
    ];
    subtract_reciprocal_air_action(&mut rhs, &state, 0, 1, coupling);
    let (derivative, residual) = factor.solve(&rhs, 1e-12, &gate).unwrap();
    // q_mech = d_air*F/(a*d_air+b^2), so at d_air=0 only d_air' survives.
    let expected = source.scale(-2.0 * omega) / (b * b);
    close(derivative[0], expected, 3e-14);
    assert!(derivative[0].abs() > 0.0 && residual < 1e-12);
    // The omitted reciprocal term would corrupt this acoustic-state slope.
    let expected_air = source / b;
    let expected_air_slope = -a * expected / b - expected_air.scale(omega.recip());
    close(derivative[1], expected_air_slope, 3e-14);
}

#[test]
fn quadratic_storage_and_force_velocity_slopes_obey_the_product_rule() {
    let omega: f64 = 650.0;
    let natural: f64 = 700.0;
    let damping = 14.0;
    let force = C64::new(0.3, 0.1);
    let dynamic = C64::new(frequency_detuning(natural, omega), -omega * damping);
    let q = force / dynamic;
    let dq = -(pole_derivative(omega, damping) * q) / dynamic;
    let velocity = C64::new(0.0, -omega) * q;
    let dv = velocity_derivative(omega, q, dq);
    let supply_slope = 0.5 * (force * dv.conj()).re;
    let heat_slope = damping * half_norm_derivative(velocity, dv);
    assert!((heat_slope - supply_slope).abs() <= 3e-14 * heat_slope.abs());
    let stored_slope = 0.5
        * (half_norm_derivative(q.scale(natural), dq.scale(natural))
            + half_norm_derivative(velocity, dv));
    let independent =
        0.5 * ((natural * natural + omega * omega) * (q.conj() * dq).re + omega * q.norm_sq());
    assert!((stored_slope - independent).abs() <= 3e-14 * independent.abs());
}
