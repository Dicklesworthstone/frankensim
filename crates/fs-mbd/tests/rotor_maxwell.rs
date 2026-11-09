//! Physical oracles for retained hereditary bearing resonance.

use fs_math::c64::C64;
use fs_mbd::rotordynamics::{
    RotorBearingModel, RotorError,
    relaxation::{
        MAX_ROTOR_MAXWELL_BRANCHES, MaxwellDampingBranch, PoleOptions, ResonanceControls,
    },
    resonance::{RotorHarmonicDirection, RotorObservation},
};

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:?}, expected={expected:?}"
    );
}
fn scalar(actual: f64, expected: f64, tolerance: f64) {
    close(C64::from_re(actual), C64::from_re(expected), tolerance);
}
fn branch(strain: &[f64], stiffness: f64, tau: f64) -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        strain: strain.to_vec(),
        stiffness,
        relaxation_time_s: tau,
    }
}
fn weights(w: f64, k: f64, tau: f64) -> (C64, C64) {
    // Direct complex constitutive equation, independent of production's
    // bounded real/reciprocal elimination and original-row implementation.
    let denominator = C64::new(1.0, w * tau);
    (
        C64::new(0.0, k * w * tau) / denominator,
        C64::new(0.0, k * tau) / (denominator * denominator),
    )
}
fn options() -> PoleOptions {
    PoleOptions {
        spectral_scale: 3.0,
        step_tolerance: 2e-13,
        residual_tolerance: 2e-12,
        max_iterations: 40,
        max_dimension: 8,
        max_dense_bytes: 16 * 1024 * 1024,
        max_work: 100_000_000,
    }
}

#[test]
fn scalar_memory_matches_original_complex_constitutive_law_and_loss_peak() {
    let model = RotorBearingModel::new(1, &[2.0], &[0.3], &[10.0], &[0.0], || false).unwrap();
    let branches = [branch(&[1.25], 6.0, 0.4)];
    let point = model
        .at_speed_with_maxwell_memory(8.0, &branches, || false)
        .unwrap();
    let target = C64::new(0.7, -0.4);
    let mut loss = Vec::new();
    for w in [-10.0, -2.5, 0.0, 0.625, 2.5, 10.0] {
        let (h, dh) = weights(w, 6.0, 0.4);
        let d = C64::new(10.0 - 2.0 * w * w, 0.3 * w) + h.scale(1.25 * 1.25);
        let dd = C64::new(-4.0 * w, 0.3) + dh.scale(1.25 * 1.25);
        let force = [d * target];
        let result = point
            .harmonic_response_with_derivative(
                w,
                &force,
                RotorHarmonicDirection {
                    frequency_derivative: 1.0,
                    spin_derivative: 0.0,
                    force_derivative: &[C64::ZERO],
                },
                3e-12,
                || false,
            )
            .unwrap();
        let q = result.response.motion.displacement[0];
        let dq = -(dd * target) / d;
        let stress = h * target.scale(1.25);
        close(q, target, 3e-12);
        close(result.displacement_derivative[0], dq, 3e-12);
        close(
            result.velocity_derivative[0],
            C64::new(0.0, w) * dq + C64::new(0.0, 1.0) * target,
            3e-12,
        );
        close(result.response.branch_stress[0], stress, 3e-12);
        close(
            result.branch_stress_derivative[0],
            (dh * target + h * dq).scale(1.25),
            3e-12,
        );
        let heat = stress.norm_sq() / (2.0 * 6.0 * 0.4);
        scalar(result.response.branch_mean_heat_w[0], heat, 3e-12);
        scalar(result.response.branch_mean_work_w[0], heat, 3e-12);
        scalar(result.response.mean_power_residual_w, 0.0, 2e-10);
        assert!(result.response.motion.relative_residual < 3e-12);
        assert!(result.response.relaxation_relative_residual < 3e-12);
        assert!(result.derivative_relative_residual < 3e-12);
        assert!(result.relaxation_derivative_relative_residual < 3e-12);
        if w > 0.0 {
            loss.push((result.response.branch_stress[0] / target.scale(1.25)).im);
        }
    }
    // The constitutive loss modulus, not an inferred damping ratio, peaks at w*tau=1.
    scalar(loss[1], 3.0, 3e-12);
    assert!(loss[1] > loss[0] && loss[1] > loss[2]);
}

#[test]
fn signed_gyro_source_and_memory_slopes_follow_the_complete_circular_equation() {
    let model = RotorBearingModel::new(
        2,
        &[2.0, 0.0, 0.0, 2.0],
        &[0.4, 0.0, 0.0, 0.4],
        &[8.0, 0.0, 0.0, 8.0],
        &[0.0, -0.6, 0.6, 0.0],
        || false,
    )
    .unwrap();
    let branches = [branch(&[1.0, 0.0], 3.0, 0.3), branch(&[0.0, 1.0], 3.0, 0.3)];
    let (dw, dspin) = (0.7, -0.2);
    let (f, df) = (C64::new(1.2, -0.4), C64::new(0.1, -0.3));
    for spin in [-2.0, 2.0] {
        let point = model
            .at_speed_with_maxwell_memory(spin, &branches, || false)
            .unwrap();
        for handedness in [-1.0, 1.0] {
            let circular = C64::new(0.0, -handedness);
            for w in [-1.7, 0.5, 2.3] {
                let (h, dh) = weights(w, 3.0, 0.3);
                let d = C64::new(8.0 - 2.0 * w * w - handedness * spin * 0.6 * w, 0.4 * w) + h;
                let dd = C64::new(
                    -4.0 * w * dw - handedness * 0.6 * (spin * dw + w * dspin),
                    0.4 * dw,
                ) + dh.scale(dw);
                let q = f / d;
                let dq = (df - dd * q) / d;
                let result = point
                    .harmonic_response_with_derivative(
                        w,
                        &[f, circular * f],
                        RotorHarmonicDirection {
                            frequency_derivative: dw,
                            spin_derivative: dspin,
                            force_derivative: &[df, circular * df],
                        },
                        4e-12,
                        || false,
                    )
                    .unwrap();
                for (i, phase) in [C64::ONE, circular].into_iter().enumerate() {
                    close(result.response.motion.displacement[i], phase * q, 4e-12);
                    close(result.displacement_derivative[i], phase * dq, 4e-12);
                    close(
                        result.velocity_derivative[i],
                        phase * C64::new(0.0, 1.0) * (dq.scale(w) + q.scale(dw)),
                        4e-12,
                    );
                    close(
                        result.branch_stress_derivative[i],
                        phase * (h * dq + dh.scale(dw) * q),
                        4e-12,
                    );
                }
                scalar(result.response.motion.mean_gyroscopic_power_w, 0.0, 3e-12);
                scalar(result.response.mean_power_residual_w, 0.0, 4e-11);
            }
        }
    }
}

#[test]
fn hereditary_bearings_preserve_nonsymmetric_transfer_and_active_baseline_growth() {
    let model = RotorBearingModel::new(
        2,
        &[1.0, 0.0, 0.0, 1.0],
        &[0.2, 0.4, 0.0, -0.1],
        &[4.0, 3.0, 0.0, 9.0],
        &[0.0; 4],
        || false,
    )
    .unwrap();
    let branches = [branch(&[1.0, 0.0], 2.0, 0.5)];
    let point = model
        .at_speed_with_maxwell_memory(-5.0, &branches, || false)
        .unwrap();
    let w = 1.4;
    let (h, dh) = weights(w, 2.0, 0.5);
    let (a, d, upper) = (
        C64::new(4.0 - w * w, 0.2 * w) + h,
        C64::new(9.0 - w * w, -0.1 * w),
        C64::new(3.0, 0.4 * w),
    );
    let (da, dd, du) = (
        C64::new(-2.0 * w, 0.2) + dh,
        C64::new(-2.0 * w, -0.1),
        C64::new(0.0, 0.4),
    );
    let q = [-upper / (a * d), C64::ONE / d];
    let dq = [
        -du / (a * d) + upper * (da * d + a * dd) / (a * d * a * d),
        -dd / (d * d),
    ];
    let response = point
        .harmonic_response_with_derivative(
            w,
            &[C64::ZERO, C64::ONE],
            RotorHarmonicDirection {
                frequency_derivative: 1.0,
                spin_derivative: 0.0,
                force_derivative: &[C64::ZERO; 2],
            },
            3e-12,
            || false,
        )
        .unwrap();
    for i in 0..2 {
        close(response.response.motion.displacement[i], q[i], 3e-12);
        close(response.displacement_derivative[i], dq[i], 3e-12);
    }
    let v = &response.response.motion.velocity;
    let damping_work = 0.5
        * (((v[0].scale(0.2) + v[1].scale(0.4)).conj() * v[0]).re
            + (v[1].scale(-0.1).conj() * v[1]).re);
    scalar(
        response.response.motion.mean_velocity_force_power_w,
        damping_work,
        2e-12,
    );
    let poles = point.poles(|| false).unwrap();
    assert_eq!(poles.poles_per_s.len(), 5);
    let unstable = C64::new(0.05, (9.0_f64 - 0.0025).sqrt());
    assert!(
        poles
            .poles_per_s
            .iter()
            .any(|&s| (s - unstable).abs() < 2e-9)
    );
}

#[test]
fn full_memory_spectrum_keeps_equal_rate_dark_states_and_original_force_residue() {
    let model = RotorBearingModel::new(1, &[1.0], &[2.0], &[1.5], &[0.0], || false).unwrap();
    let branches = [
        branch(&[1.0], 0.75, 0.25),
        branch(&[1.0], 0.75, 0.25),
        branch(&[0.0], 7.0, 0.2),
    ];
    let point = model
        .at_speed_with_maxwell_memory(9.0, &branches, || false)
        .unwrap();
    // (s+4)(s^2+2s+1.5)+1.5s=(s+1)(s+2)(s+3).
    // The independent equal-rate stress difference adds -4; zero strain adds -5.
    let spectrum = point.poles(|| false).unwrap();
    assert_eq!(spectrum.memory_states, 3);
    assert_eq!(spectrum.poles_per_s.len(), 5);
    for (actual, expected) in spectrum
        .poles_per_s
        .iter()
        .zip([-5.0, -4.0, -3.0, -2.0, -1.0])
    {
        close(*actual, C64::from_re(expected), 2e-9);
    }
    let prepared = point.prepare_modes(options(), || false).unwrap();
    let mode = prepared
        .refine(C64::new(-1.98, 0.01), &[C64::ONE], || false)
        .unwrap();
    close(mode.laplace_frequency_per_s, C64::from_re(-2.0), 2e-11);
    close(mode.branch_stress[0], C64::from_re(-0.75), 2e-11);
    close(mode.branch_stress[1], C64::from_re(-0.75), 2e-11);
    close(mode.branch_stress[2], C64::ZERO, 2e-11);
    let prepared_residue = prepared
        .refine_residue(C64::from_re(-2.0), &[C64::ONE], &[C64::ONE], 1e-10, || {
            false
        })
        .unwrap();
    let residue = prepared_residue
        .force_residue(&[C64::ONE], || false)
        .unwrap();
    // T'(s0)=-1/2, hence residue=1/T'=-2, retaining each original stress.
    close(residue.displacement[0], C64::from_re(-2.0), 2e-11);
    close(residue.velocity[0], C64::from_re(4.0), 2e-11);
    close(residue.branch_stress[0], C64::from_re(1.5), 2e-11);
    close(residue.branch_stress[1], C64::from_re(1.5), 2e-11);
    assert!(
        prepared
            .refine(C64::from_re(-4.0), &[C64::ONE], || false)
            .is_err()
    );
}

#[test]
fn driven_memory_mobility_peak_and_half_height_follow_exact_rational_curve() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let branches = [branch(&[1.0], 1.0, 1.0)];
    let point = model
        .at_speed_with_maxwell_memory(0.0, &branches, || false)
        .unwrap();
    let result = point
        .driven_resonance(
            (0.05, 2.0),
            &[C64::ONE],
            RotorObservation::Velocity(&[1.0]),
            3e-12,
            ResonanceControls {
                absolute_frequency_tolerance_rad_s: 1e-11,
                relative_frequency_tolerance: 1e-12,
                maximum_evaluations: 512,
                measure_half_power_bandwidth: true,
            },
            || false,
        )
        .unwrap();
    // For y=w^2, inverse mobility magnitude squared is
    // (y^2-y+1)/(1+y); its minimum is at y=sqrt(3)-1.
    let peak = (3.0_f64.sqrt() - 1.0).sqrt();
    let minimum = 2.0 * 3.0_f64.sqrt() - 3.0;
    let level = 2.0 * minimum;
    let discriminant = ((1.0 + level) * (1.0 + level) - 4.0 * (1.0 - level)).sqrt();
    let low = ((1.0 + level - discriminant) / 2.0).sqrt();
    let high = ((1.0 + level + discriminant) / 2.0).sqrt();
    scalar(result.peak.angular_frequency_rad_s, peak, 3e-10);
    scalar(result.observation.norm_sq(), 1.0 / minimum, 3e-10);
    let width = result.half_power.unwrap();
    scalar(width.lower.angular_frequency_rad_s, low, 3e-10);
    scalar(width.upper.angular_frequency_rad_s, high, 3e-10);
    scalar(width.bandwidth_rad_s, high - low, 3e-10);
}

#[test]
fn empty_memory_reduces_to_native_rotor_and_refusals_are_atomic() {
    let model = RotorBearingModel::new(1, &[2.0], &[0.4], &[8.0], &[0.0], || false).unwrap();
    let native = model.at_speed(-3.0, || false).unwrap();
    let empty = model
        .at_speed_with_maxwell_memory(-3.0, &[], || false)
        .unwrap();
    let force = [C64::new(1.0, 0.4)];
    let direction = RotorHarmonicDirection {
        frequency_derivative: 0.8,
        spin_derivative: -0.2,
        force_derivative: &[C64::new(-0.1, 0.3)],
    };
    let a = native
        .harmonic_response_with_derivative(1.9, &force, direction, 2e-12, || false)
        .unwrap();
    let b = empty
        .harmonic_response_with_derivative(1.9, &force, direction, 2e-12, || false)
        .unwrap();
    close(
        a.response.displacement[0],
        b.response.motion.displacement[0],
        2e-12,
    );
    close(
        a.displacement_derivative[0],
        b.displacement_derivative[0],
        2e-12,
    );
    assert!(b.response.branch_stress.is_empty());
    assert_eq!(
        empty.poles(|| false).unwrap().poles_per_s,
        native.poles(|| false).unwrap().poles_per_s
    );
    for invalid in [
        branch(&[1.0], 0.0, 1.0),
        branch(&[1.0], 1.0, -1.0),
        branch(&[], 1.0, 1.0),
    ] {
        assert!(
            model
                .at_speed_with_maxwell_memory(0.0, &[invalid], || false)
                .is_err()
        );
    }
    let too_many = vec![branch(&[1.0], 1.0, 1.0); MAX_ROTOR_MAXWELL_BRANCHES + 1];
    assert!(
        model
            .at_speed_with_maxwell_memory(0.0, &too_many, || false)
            .is_err()
    );
    let memories = [branch(&[1.0], 1.0, 1.0)];
    let point = model
        .at_speed_with_maxwell_memory(0.0, &memories, || false)
        .unwrap();
    assert!(matches!(point.poles(|| true), Err(RotorError::Cancelled)));
    assert!(matches!(
        point.harmonic_response(1.0, &force, 2e-12, || true),
        Err(RotorError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        point.harmonic_response_with_derivative(1.0, &force, direction, 2e-12, || {
            polls += 1;
            polls >= 8
        }),
        Err(RotorError::Cancelled)
    ));
    let first = point
        .harmonic_response(1.0, &force, 2e-12, || false)
        .unwrap();
    let second = point
        .harmonic_response(1.0, &force, 2e-12, || false)
        .unwrap();
    assert_eq!(first, second);
    let free = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let free_point = free
        .at_speed_with_maxwell_memory(0.0, &memories, || false)
        .unwrap();
    assert!(
        free_point
            .harmonic_response(0.0, &[C64::ONE], 2e-12, || false)
            .is_err()
    );
    // omega*tau=2^1200 loses the positive loss quadrature below binary64.
    // Its heat can still be representable: refuse rather than silently turn
    // this material into a conservative spring in the forced-response solve.
    let extreme = RotorBearingModel::new(
        1,
        &[2.0_f64.powi(-1000)],
        &[0.0],
        &[2.0_f64.powi(200)],
        &[0.0],
        || false,
    )
    .unwrap();
    let slow = [branch(&[1.0], 1.0, 2.0_f64.powi(600))];
    let extreme_point = extreme
        .at_speed_with_maxwell_memory(0.0, &slow, || false)
        .unwrap();
    assert!(matches!(
        extreme_point.harmonic_response(2.0_f64.powi(600), &[C64::ONE], 2e-12, || false),
        Err(RotorError::Invalid("Maxwell loss quadrature underflowed"))
    ));
}
