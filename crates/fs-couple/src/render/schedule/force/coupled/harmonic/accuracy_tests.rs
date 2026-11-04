//! G1/G3: original linear solves across resonant detuning and drive scales.
use super::*;
use fs_math::dd::Dd;

fn relative_close(actual: C64, expected: C64, tolerance: f64) {
    let scale = component_scale(expected);
    assert!(scale > 0.0);
    let error = C64::new(
        actual.re / scale - expected.re / scale,
        actual.im / scale - expected.im / scale,
    )
    .abs();
    assert!(
        error <= tolerance,
        "actual={actual:?}, expected={expected:?}, error={error}"
    );
}

#[test]
fn high_q_detuning_agrees_with_independent_extended_precision_squares() {
    for natural in [core::f64::consts::TAU * 440.0, 12_345.678_9] {
        for offset in [1_u64, 17, 113] {
            for drive in [
                f64::from_bits(natural.to_bits() - offset),
                f64::from_bits(natural.to_bits() + offset),
            ] {
                let w = Dd::from_f64(natural);
                let d = Dd::from_f64(drive);
                let expected = (w * w - d * d).to_f64();
                let actual = frequency_detuning(natural, drive);
                assert!(
                    (actual / expected - 1.0).abs() <= 3.0 * f64::EPSILON,
                    "natural={natural}, drive={drive}, actual={actual}, expected={expected}"
                );
            }
        }
        assert_eq!(frequency_detuning(natural, natural), 0.0);
        assert_eq!(frequency_detuning(natural, 0.0), natural * natural);
    }
}

#[test]
fn weak_resonant_drive_survives_stiff_coordinate_and_larger_force() {
    // The old F/operator_scale sends the first force to zero. Its residual
    // then hides beneath the stiff second row, while the missing response is
    // sixteen orders larger than the response that survives.
    let gate = CancelGate::new_clock_free();
    let diagonal = [C64::new(0.0, -2e-30), C64::new(3e6, -10.0)];
    let matrix = [diagonal[0], C64::ZERO, C64::ZERO, diagonal[1]];
    let force = [C64::from_re(1e-320), C64::from_re(1e-300)];
    let factor = HarmonicFactor::new(&matrix, 2, &gate).unwrap();
    let (response, residual) = factor.solve(&force, 1e-12, &gate).unwrap();
    assert!(residual < 1e-12);
    for i in 0..2 {
        relative_close(response[i], force[i] / diagonal[i], 2e-14);
    }
    assert!(response[0].abs() > 1e15 * response[1].abs());
}

#[test]
fn harmonic_restores_finite_response_when_scale_ratio_underflows_or_overflows() {
    let gate = CancelGate::new_clock_free();
    let matrix = [
        C64::new(0.0, -1e-35),
        C64::ZERO,
        C64::ZERO,
        C64::from_re(1e6),
    ];
    let force = [C64::from_re(2e-320), C64::ZERO];
    assert_eq!(force[0].re / 1e6, 0.0);
    let (response, _) = solve_harmonic(&matrix, 2, &force, 1e-12, &gate).unwrap();
    relative_close(response[0], force[0] / matrix[0], 2e-14);
    assert_eq!(response[1], C64::ZERO);

    // Both real/imaginary components are finite even though the naive real
    // scale ratio overflows. This probes the same physical amplitude restore.
    let matrix = [C64::new(1e-308, 1e-308)];
    let (response, _) = solve_harmonic(&matrix, 1, &[C64::from_re(2.0)], 1e-12, &gate).unwrap();
    relative_close(response[0], C64::new(1e308, -1e308), 2e-14);
}

#[test]
fn shared_harmonic_factor_preserves_signed_complex_amplitude_and_exact_silence() {
    let gate = CancelGate::new_clock_free();
    let matrix = [
        C64::new(0.002, -0.003),
        C64::new(-0.0004, 0.0001),
        C64::new(-0.0004, 0.0001),
        C64::new(3.0, -0.3),
    ];
    let base = [C64::new(0.3, -0.1), C64::new(-0.2, 0.05)];
    let factor = HarmonicFactor::new(&matrix, 2, &gate).unwrap();
    let (reference, _) = factor.solve(&base, 1e-12, &gate).unwrap();
    for amplitude in [-1e150, -1e-150, 1e-150, 1e150] {
        let force = base.map(|f| f.scale(amplitude));
        let (response, residual) = factor.solve(&force, 1e-12, &gate).unwrap();
        assert!(residual < 1e-12);
        for i in 0..2 {
            relative_close(response[i], reference[i].scale(amplitude), 2e-14);
        }
    }
    let (silence, residual) = factor.solve(&[C64::ZERO; 2], 1e-12, &gate).unwrap();
    assert_eq!(silence, [C64::ZERO; 2]);
    assert_eq!(residual, 0.0);
}

#[test]
fn normalization_preserves_original_finite_solutions_at_extreme_resonance() {
    // Normalizing this RHS to one would overflow the intermediate solution,
    // although the source force and the original physical response both fit.
    let gate = CancelGate::new_clock_free();
    let matrix = [C64::from_re(1e-320), C64::ZERO, C64::ZERO, C64::ONE];
    let force = [C64::from_re(1e-100), C64::ZERO];
    let (response, residual) = solve_harmonic(&matrix, 2, &force, 1e-12, &gate).unwrap();
    relative_close(response[0], force[0] / matrix[0], 2e-14);
    assert_eq!(response[1], C64::ZERO);
    assert!(residual < 1e-12);

    // A single global RHS normalization cannot retain this amplitude ratio.
    // The original equation can, and both coordinates must remain present.
    let identity = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    let force = [C64::from_re(1e300), C64::from_re(1e-300)];
    let (response, residual) = solve_harmonic(&identity, 2, &force, 1e-12, &gate).unwrap();
    assert_eq!(response, force);
    assert!(residual < 1e-12);
}

#[test]
fn harmonic_accuracy_keeps_singularity_cancellation_and_representability_refusals() {
    let gate = CancelGate::new_clock_free();
    assert!(HarmonicFactor::new(&[C64::ONE, C64::ONE, C64::ONE, C64::ONE], 2, &gate).is_err());
    let matrix = [C64::from_re(1e-10)];
    let factor = HarmonicFactor::new(&matrix, 1, &gate).unwrap();
    assert!(factor.solve(&[C64::from_re(1e308)], 1e-12, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        factor.solve(&[C64::ZERO], 1e-12, &cancelled),
        Err(ModalCouplingError::Cancelled)
    ));
    // Do not silently remove an operator coefficient or accept nonfinite
    // input. No threshold truncation is used.
    assert!(
        HarmonicFactor::new(
            &[
                C64::from_re(1e308),
                C64::from_re(1e-308),
                C64::ZERO,
                C64::ONE
            ],
            2,
            &gate
        )
        .is_err()
    );
    let matrix = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    let identity = HarmonicFactor::new(&matrix, 2, &gate).unwrap();
    assert!(
        identity
            .solve(&[C64::new(f64::NAN, 0.0), C64::ZERO], 1e-12, &gate)
            .is_err()
    );
}
