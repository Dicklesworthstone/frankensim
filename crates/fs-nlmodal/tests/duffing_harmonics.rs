//! Physical derivative spectra, original-ODE checks and retained-band tails.

use fs_nlmodal::{
    DuffingBackboneError, DuffingKinematicSpectrum, duffing_kinematic_spectrum_exact,
    duffing_spectrum_exact,
};

fn close(got: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (got - expected).abs() <= absolute + relative * expected.abs(),
        "got {got:.17e}, expected {expected:.17e}"
    );
}

fn state_at_phase(spectrum: &DuffingKinematicSpectrum, phase: f64) -> [f64; 3] {
    let mut result = [0.0; 3];
    for (j, &coefficient) in spectrum.displacement.cosine_amplitudes.iter().enumerate() {
        let (sine, cosine) = ((2 * j + 1) as f64 * phase).sin_cos();
        result[0] += coefficient * cosine;
        result[1] += spectrum.velocity_sine_amplitudes[j] * sine;
        result[2] += spectrum.acceleration_cosine_amplitudes[j] * cosine;
    }
    result
}

// Independent energy-separated period, without elliptic functions or returned
// spectral data. RK4 then integrates the original physical force directly.
fn energy_period(omega0: f64, beta: f64, amplitude: f64) -> f64 {
    let count = 16_384;
    let step = core::f64::consts::FRAC_PI_2 / f64::from(count);
    let mut integral = 0.0;
    for i in 0..=count {
        let sine = (f64::from(i) * step).sin();
        let weight = if i == 0 || i == count {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight
            / (omega0 * omega0 + 0.5 * beta * amplitude * amplitude * (1.0 + sine * sine)).sqrt();
    }
    4.0 * step * integral / 3.0
}

fn rk4_step(x: [f64; 2], dt: f64, omega0: f64, beta: f64) -> [f64; 2] {
    let force = |q: f64| -(omega0 * omega0 + beta * q * q) * q;
    let k1 = [x[1], force(x[0])];
    let k2 = [x[1] + 0.5 * dt * k1[1], force(x[0] + 0.5 * dt * k1[0])];
    let k3 = [x[1] + 0.5 * dt * k2[1], force(x[0] + 0.5 * dt * k2[0])];
    let k4 = [x[1] + dt * k3[1], force(x[0] + dt * k3[0])];
    core::array::from_fn(|i| x[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
}

#[test]
fn g1_physical_velocity_and_acceleration_match_the_original_ode() {
    for (omega0, beta, amplitude) in [
        (2.0, -3.96, 1.0),
        (2.0, -0.4, 1.0),
        (2.0, 0.0004, 1.0),
        (2.0, 4.0, 1.0),
        (2.0, 40_000.0, 1.0),
        (0.0, 3.0, 0.7),
    ] {
        let spectrum = duffing_kinematic_spectrum_exact(omega0, beta, amplitude, 32).unwrap();
        assert_eq!(
            spectrum.displacement,
            duffing_spectrum_exact(omega0, beta, amplitude, 32).unwrap()
        );
        let period = energy_period(omega0, beta, amplitude);
        let independent_frequency = core::f64::consts::TAU / period;
        let count = 32_768;
        let dt = period / f64::from(count);
        let mut x = [amplitude, 0.0];
        let mut measured_velocity = [0.0; 4];
        let mut measured_acceleration = [0.0; 4];
        for i in 0..count {
            let phase = core::f64::consts::TAU * f64::from(i) / f64::from(count);
            let acceleration = -(omega0 * omega0 + beta * x[0] * x[0]) * x[0];
            for j in 0..4 {
                let (sine, cosine) = ((2 * j + 1) as f64 * phase).sin_cos();
                measured_velocity[j] += x[1] * sine;
                measured_acceleration[j] += acceleration * cosine;
            }
            if i % 512 == 0 {
                let predicted = state_at_phase(&spectrum, phase);
                close(predicted[0] / amplitude, x[0] / amplitude, 0.0, 3e-10);
                close(
                    predicted[1] / (amplitude * independent_frequency),
                    x[1] / (amplitude * independent_frequency),
                    0.0,
                    3e-10,
                );
                close(
                    predicted[2] / (amplitude * independent_frequency.powi(2)),
                    acceleration / (amplitude * independent_frequency.powi(2)),
                    0.0,
                    3e-10,
                );
            }
            x = rk4_step(x, dt, omega0, beta);
        }
        for j in 0..4 {
            close(
                spectrum.velocity_sine_amplitudes[j] / (amplitude * independent_frequency),
                2.0 * measured_velocity[j] / f64::from(count) / (amplitude * independent_frequency),
                0.0,
                3e-10,
            );
            close(
                spectrum.acceleration_cosine_amplitudes[j]
                    / (amplitude * independent_frequency.powi(2)),
                2.0 * measured_acceleration[j]
                    / f64::from(count)
                    / (amplitude * independent_frequency.powi(2)),
                0.0,
                3e-10,
            );
        }
    }
}

#[test]
fn g1_derivative_tails_bound_the_omitted_physical_band_and_preserve_refinement() {
    for beta in [-0.999_999, -0.9, 1.0, 1e6] {
        let reference = duffing_kinematic_spectrum_exact(1.0, beta, 1.0, 512).unwrap();
        let mut previous = f64::INFINITY;
        for count in [1, 3, 8, 32] {
            let coarse = duffing_kinematic_spectrum_exact(1.0, beta, 1.0, count).unwrap();
            assert_eq!(
                coarse.velocity_sine_amplitudes,
                reference.velocity_sine_amplitudes[..count]
            );
            assert_eq!(
                coarse.acceleration_cosine_amplitudes,
                reference.acceleration_cosine_amplitudes[..count]
            );
            let omitted_velocity: f64 = reference.velocity_sine_amplitudes[count..]
                .iter()
                .map(|a| a.abs())
                .sum();
            let omitted_acceleration: f64 = reference.acceleration_cosine_amplitudes[count..]
                .iter()
                .map(|a| a.abs())
                .sum();
            assert!(omitted_velocity <= coarse.velocity_tail_bound * (1.0 + 3e-13));
            assert!(omitted_acceleration <= coarse.acceleration_tail_bound * (1.0 + 3e-13));
            assert!(coarse.acceleration_tail_bound < previous);
            previous = coarse.acceleration_tail_bound;
        }
    }
    // A displacement-only criterion badly understates derivative bandwidth.
    let one = duffing_kinematic_spectrum_exact(0.0, 1.0, 1.0, 1).unwrap();
    let omega = one.displacement.angular_frequency;
    assert!(one.velocity_tail_bound / omega > 3.0 * one.displacement.displacement_tail_bound);
    assert!(
        one.acceleration_tail_bound / omega.powi(2)
            > 9.0 * one.displacement.displacement_tail_bound
    );
}

#[test]
fn g3_near_separatrix_derivatives_keep_the_original_well_margin() {
    let t = 2.0_f64.powi(-53);
    let amplitude = 1.0 - t;
    let beta = -(1.0 + 2.0 * t);
    assert_eq!(1.0 + beta * (amplitude * amplitude), 0.0);
    let spectrum = duffing_kinematic_spectrum_exact(1.0, beta, amplitude, 512).unwrap();
    assert!(spectrum.displacement.cosine_amplitudes[0] > amplitude);
    assert!(spectrum.velocity_sine_amplitudes[0] < 0.0);
    assert!(spectrum.velocity_sine_amplitudes[1] > 0.0);
    assert!(spectrum.acceleration_tail_bound < 1e-20);
    for i in 0..256 {
        let phase = core::f64::consts::TAU * f64::from(i) / 256.0;
        let [q, v, acceleration] = state_at_phase(&spectrum, phase);
        close(acceleration, -(1.0 + beta * q * q) * q, 0.0, 3e-12);
        close(
            0.5 * v * v + 0.5 * q * q + 0.25 * beta * q.powi(4),
            0.5 * amplitude * amplitude + 0.25 * beta * amplitude.powi(4),
            0.0,
            3e-12,
        );
    }
    let coarse = duffing_kinematic_spectrum_exact(1.0, beta, amplitude, 8).unwrap();
    assert!(coarse.acceleration_tail_bound > 0.01);
}

#[test]
fn g3_derivatives_survive_a_rounded_zero_displacement_and_frequency_square_overflow() {
    for (rate_exp, amplitude_exp, beta_exp) in [(500, -100, 0), (700, -400, 600)] {
        let rate = 2.0_f64.powi(rate_exp);
        let amplitude = 2.0_f64.powi(amplitude_exp);
        let beta = 2.0_f64.powi(beta_exp);
        let spectrum = duffing_kinematic_spectrum_exact(rate, beta, amplitude, 2).unwrap();
        assert_eq!(spectrum.displacement.cosine_amplitudes[1], 0.0);
        // The Lindstedt remainder is below 2^-1200 relative in these cases.
        let velocity_exp = beta_exp + 3 * amplitude_exp - rate_exp - 5;
        let acceleration_exp = beta_exp + 3 * amplitude_exp - 5;
        close(
            spectrum.velocity_sine_amplitudes[1],
            -3.0 * 2.0_f64.powi(velocity_exp),
            6e-12,
            0.0,
        );
        close(
            spectrum.acceleration_cosine_amplitudes[1],
            -9.0 * 2.0_f64.powi(acceleration_exp),
            6e-12,
            0.0,
        );
        let coarse = duffing_kinematic_spectrum_exact(rate, beta, amplitude, 1).unwrap();
        close(
            coarse.acceleration_tail_bound,
            spectrum.acceleration_cosine_amplitudes[1].abs(),
            6e-12,
            0.0,
        );
        assert!(coarse.acceleration_tail_bound.is_normal());
    }
}

#[test]
fn g3_velocity_and_acceleration_follow_physical_unit_scaling() {
    for beta in [-0.75, 3.0] {
        let reference = duffing_kinematic_spectrum_exact(1.0, beta, 1.0, 16).unwrap();
        for (time_exp, amplitude_exp) in [(400, 0), (-400, 0), (200, 400), (-200, -400)] {
            let scaled = duffing_kinematic_spectrum_exact(
                2.0_f64.powi(time_exp),
                beta * 2.0_f64.powi(2 * (time_exp - amplitude_exp)),
                2.0_f64.powi(amplitude_exp),
                16,
            )
            .unwrap();
            let velocity_scale = 2.0_f64.powi(time_exp + amplitude_exp);
            let acceleration_scale = 2.0_f64.powi(2 * time_exp + amplitude_exp);
            for j in 0..16 {
                close(
                    scaled.velocity_sine_amplitudes[j] / velocity_scale,
                    reference.velocity_sine_amplitudes[j],
                    3e-12,
                    0.0,
                );
                close(
                    scaled.acceleration_cosine_amplitudes[j] / acceleration_scale,
                    reference.acceleration_cosine_amplitudes[j],
                    3e-12,
                    0.0,
                );
            }
            close(
                scaled.velocity_tail_bound / velocity_scale,
                reference.velocity_tail_bound,
                3e-12,
                0.0,
            );
            close(
                scaled.acceleration_tail_bound / acceleration_scale,
                reference.acceleration_tail_bound,
                3e-12,
                0.0,
            );
        }
    }
}

#[test]
fn g0_opt_in_derivative_refusal_keeps_the_displacement_query_available() {
    let rate = 2.0_f64.powi(600);
    assert!(duffing_spectrum_exact(rate, 0.0, 1.0, 8).is_ok());
    assert_eq!(
        duffing_kinematic_spectrum_exact(rate, 0.0, 1.0, 8),
        Err(DuffingBackboneError::NumericalRange)
    );
    let zero = duffing_kinematic_spectrum_exact(rate, 1.0, 0.0, 8).unwrap();
    assert!(
        zero.velocity_sine_amplitudes
            .iter()
            .chain(&zero.acceleration_cosine_amplitudes)
            .all(|&a| a == 0.0)
    );
    assert_eq!(zero.velocity_tail_bound, 0.0);
    assert_eq!(zero.acceleration_tail_bound, 0.0);
    assert_eq!(
        duffing_kinematic_spectrum_exact(1.0, -1.0, 1.0, 8),
        Err(DuffingBackboneError::NoPeriodicOrbit)
    );
    assert_eq!(
        duffing_kinematic_spectrum_exact(1.0, 1.0, 1.0, 0),
        Err(DuffingBackboneError::InvalidInput)
    );
}
