use super::{DuffingBackboneError, DuffingSpectrum, MAX_DUFFING_HARMONICS, duffing_spectrum_exact};
use crate::{SosModalStorage, StressChannel};

fn close(got: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (got - expected).abs() <= absolute + relative * expected.abs(),
        "got {got:.17e}, expected {expected:.17e}"
    );
}

fn waveform(spectrum: &DuffingSpectrum, phase: f64) -> (f64, f64) {
    let mut displacement = 0.0;
    let mut velocity = 0.0;
    for (j, &amplitude) in spectrum.cosine_amplitudes.iter().enumerate() {
        let order = (2 * j + 1) as f64;
        let (sine, cosine) = (order * phase).sin_cos();
        displacement += amplitude * cosine;
        velocity -= amplitude * sine * order * spectrum.angular_frequency;
    }
    (displacement, velocity)
}

// Independent energy-separated quarter-period quadrature; no elliptic
// functions, nome or spectral coefficients enter the trajectory oracle.
fn energy_period(omega0: f64, beta: f64, amplitude: f64) -> f64 {
    let count = 16_384;
    let step = core::f64::consts::FRAC_PI_2 / f64::from(count);
    let mut sum = 0.0;
    for j in 0..=count {
        let cosine = (f64::from(j) * step).cos();
        let weight = if j == 0 || j == count {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight
            / (omega0 * omega0 + 0.5 * beta * amplitude * amplitude * (1.0 + cosine * cosine))
                .sqrt();
    }
    4.0 * step * sum / 3.0
}

fn rk4(x: f64, v: f64, step: f64, omega0: f64, beta: f64) -> (f64, f64) {
    let acceleration = |q: f64| -(omega0 * omega0 + beta * q * q) * q;
    let (q1, v1) = (v, acceleration(x));
    let (q2, v2) = (v + 0.5 * step * v1, acceleration(x + 0.5 * step * q1));
    let (q3, v3) = (v + 0.5 * step * v2, acceleration(x + 0.5 * step * q2));
    let (q4, v4) = (v + step * v3, acceleration(x + step * q3));
    (
        x + step * (q1 + 2.0 * q2 + 2.0 * q3 + q4) / 6.0,
        v + step * (v1 + 2.0 * v2 + 2.0 * v3 + v4) / 6.0,
    )
}

#[test]
fn g1_harmonics_waveforms_and_tail_match_independent_physical_trajectories() {
    for (omega0, beta) in [
        (1.0, -0.95),
        (1.0, -0.5),
        (1.0, 0.3),
        (1.0, 20.0),
        (0.0, 1.0),
    ] {
        let amplitude = 1.0;
        let resolved = duffing_spectrum_exact(omega0, beta, amplitude, 64).unwrap();
        let coarse = duffing_spectrum_exact(omega0, beta, amplitude, 2).unwrap();
        let period = energy_period(omega0, beta, amplitude);
        close(
            resolved.angular_frequency * period,
            core::f64::consts::TAU,
            2e-13,
            0.0,
        );
        let steps = 16_384;
        let step = period / f64::from(steps);
        let (mut x, mut v) = (amplitude, 0.0);
        let mut projected = [0.0; 4];
        for j in 0..steps {
            let phase = core::f64::consts::TAU * f64::from(j) / f64::from(steps);
            for (n, coefficient) in projected.iter_mut().enumerate() {
                *coefficient += 2.0 * x * (((2 * n + 1) as f64) * phase).cos() / f64::from(steps);
            }
            if j % 128 == 0 {
                let physical_phase = resolved.angular_frequency * step * f64::from(j);
                let (q, velocity) = waveform(&resolved, physical_phase);
                close(q, x, 0.0, 2e-10);
                close(velocity, v, 0.0, 2e-9);
                assert!(
                    (waveform(&coarse, physical_phase).0 - x).abs()
                        <= coarse.displacement_tail_bound + 2e-10
                );
            }
            (x, v) = rk4(x, v, step, omega0, beta);
        }
        close(x, amplitude, 0.0, 3e-10);
        close(v, 0.0, 0.0, 3e-9);
        for (got, expected) in resolved.cosine_amplitudes.iter().zip(projected) {
            close(*got, expected, 0.0, 3e-10);
        }
    }
}

#[test]
fn g1_linear_limit_and_weak_nonlinear_harmonic_signs() {
    let linear = duffing_spectrum_exact(7.0, 0.0, 0.4, 5).unwrap();
    assert_eq!(linear.cosine_amplitudes, [0.4, 0.0, 0.0, 0.0, 0.0]);
    assert_eq!(linear.displacement_tail_bound, 0.0);
    assert_eq!(
        duffing_spectrum_exact(7.0, 1.0, 0.0, 3)
            .unwrap()
            .cosine_amplitudes,
        [0.0; 3]
    );
    for ratio in [-1e-5, 1e-5] {
        let result = duffing_spectrum_exact(1.0, ratio, 1.0, 3).unwrap();
        // Lindstedt expansion at a PHYSICAL turning amplitude:
        // q/A = cos(theta) + ratio*(cos(3theta)-cos(theta))/32 + O(ratio^2).
        close(
            (result.cosine_amplitudes[0] - 1.0) / ratio,
            -1.0 / 32.0,
            0.0,
            4e-7,
        );
        close(result.cosine_amplitudes[1] / ratio, 1.0 / 32.0, 0.0, 4e-7);
        assert!(result.cosine_amplitudes[2] > 0.0);
    }
}

#[test]
fn g3_retained_order_never_renormalizes_physical_amplitudes() {
    for beta in [-0.999, -0.2, 0.2, 100.0] {
        let reference = duffing_spectrum_exact(1.0, beta, 1.0, 128).unwrap();
        let mut previous_tail = f64::INFINITY;
        for count in [1, 2, 4, 8, 16, 32] {
            let result = duffing_spectrum_exact(1.0, beta, 1.0, count).unwrap();
            assert_eq!(result.angular_frequency, reference.angular_frequency);
            assert_eq!(
                result.cosine_amplitudes,
                reference.cosine_amplitudes[..count]
            );
            assert!(result.displacement_tail_bound < previous_tail);
            previous_tail = result.displacement_tail_bound;
            for j in 0..257 {
                let phase = core::f64::consts::TAU * f64::from(j) / 257.0;
                assert!(
                    (waveform(&result, phase).0 - waveform(&reference, phase).0).abs()
                        <= result.displacement_tail_bound + 2e-14
                );
            }
        }
        if beta > 0.0 {
            assert!(reference.cosine_amplitudes[0] < 1.0);
            assert!(reference.cosine_amplitudes.iter().all(|x| *x >= 0.0));
        } else {
            assert!(reference.cosine_amplitudes[0] > 1.0);
            assert!(reference.cosine_amplitudes[1] < 0.0);
        }
    }
}

#[test]
fn g3_physical_units_and_subnormal_nome_preserve_harmonics() {
    for beta in [-0.75, 3.0] {
        let reference = duffing_spectrum_exact(1.0, beta, 1.0, 8).unwrap();
        for (time_exponent, amplitude_exponent) in [(400, 0), (-400, 0), (800, 400), (-800, -400)] {
            let rate = 2.0_f64.powi(time_exponent);
            let amplitude = 2.0_f64.powi(amplitude_exponent);
            let scaled_beta = beta * 2.0_f64.powi(2 * time_exponent - 2 * amplitude_exponent);
            let result = duffing_spectrum_exact(rate, scaled_beta, amplitude, 8).unwrap();
            close(
                result.angular_frequency / rate,
                reference.angular_frequency,
                3e-13,
                0.0,
            );
            for (&got, &expected) in result
                .cosine_amplitudes
                .iter()
                .zip(&reference.cosine_amplitudes)
            {
                close(got / amplitude, expected, 3e-12, 0.0);
            }
            close(
                result.displacement_tail_bound / amplitude,
                reference.displacement_tail_bound,
                3e-12,
                0.0,
            );
        }
    }
    let extreme = duffing_spectrum_exact(
        2.0_f64.powi(1000),
        2.0_f64.powi(-1000),
        2.0_f64.powi(800),
        2,
    )
    .unwrap();
    // beta*A^2/omega0^2 = 2^-1400; q underflows, but A*q does not.
    close(extreme.cosine_amplitudes[1], 2.0_f64.powi(-605), 5e-12, 0.0);
    assert!(extreme.displacement_tail_bound > 0.0);
}

#[test]
fn g0_softening_separatrix_and_parameter_refusals() {
    let t = 2.0_f64.powi(-53);
    let amplitude = 1.0 - t;
    let result = duffing_spectrum_exact(1.0, -(1.0 + 2.0 * t), amplitude, 512).unwrap();
    assert!(result.angular_frequency > 0.0);
    close(waveform(&result, 0.0).0, amplitude, 0.0, 3e-14);
    assert!(result.cosine_amplitudes[0] > amplitude);
    assert!(result.cosine_amplitudes[1] < 0.0);
    assert_eq!(
        duffing_spectrum_exact(1.0, -1.0, 1.0, 8),
        Err(DuffingBackboneError::NoPeriodicOrbit)
    );
    for count in [0, MAX_DUFFING_HARMONICS + 1, usize::MAX] {
        assert_eq!(
            duffing_spectrum_exact(1.0, 1.0, 1.0, count),
            Err(DuffingBackboneError::InvalidInput)
        );
    }
    assert_eq!(
        duffing_spectrum_exact(f64::NAN, 1.0, 1.0, 8),
        Err(DuffingBackboneError::InvalidInput)
    );
}

#[test]
fn g1_source_adapter_keeps_every_channel_and_refuses_coupled_axes() {
    let mut storage = SosModalStorage {
        omegas: vec![7.0, 11.0],
        channels: vec![
            StressChannel {
                coefficient: 2.0,
                coupling: vec![2.0, 0.0, 0.0, -3.0],
            },
            StressChannel {
                coefficient: 3.0,
                coupling: vec![-4.0, 0.0, 0.0, 5.0],
            },
        ],
    };
    for (mode, beta) in [(0, 56.0), (1, 93.0)] {
        let actual = storage.invariant_axis_spectrum(mode, 0.7, 8).unwrap();
        let expected = duffing_spectrum_exact(storage.omegas[mode], beta, 0.7, 8).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(
            actual.angular_frequency,
            storage.invariant_axis_backbone(mode, 0.7).unwrap()
        );
    }
    storage.channels[0].coupling[1] = 0.5;
    storage.channels[0].coupling[2] = 0.5;
    assert!(storage.invariant_axis_spectrum(0, 0.7, 8).is_err());
    assert!(storage.invariant_axis_spectrum(2, 0.7, 8).is_err());
}
