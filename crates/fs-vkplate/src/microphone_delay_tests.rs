//! Analytic propagation oracles: delay moments, complex gain, and actual PCM
//! interference. These tests qualify interpolation, not measured instruments.

use super::*;

fn response(start: usize, weights: &[f64], omega: f64) -> [f64; 2] {
    let mut result = [0.0; 2];
    for (i, &weight) in weights.iter().enumerate() {
        let phase = omega * (start + i) as f64;
        result[0] += weight * phase.cos();
        result[1] -= weight * phase.sin();
    }
    result
}

fn complex_error(actual: [f64; 2], expected: [f64; 2]) -> f64 {
    (actual[0] - expected[0]).hypot(actual[1] - expected[1])
}

#[test]
fn g0_fractional_propagation_preserves_time_moments_without_added_latency() {
    for delay in [1.125_f64, 1.5, 1.875, 12.125, 12.5, 12.875] {
        let (start, weights) = microphone_delay_taps(delay);
        for order in 0..=3 {
            let actual = weights
                .iter()
                .enumerate()
                .map(|(i, weight)| weight * ((start + i) as f64).powi(order))
                .sum::<f64>();
            let expected = delay.powi(order);
            assert!((actual - expected).abs() < 2e-14 * expected.max(1.0));
        }
    }
    for delay in [0.0, 1.0, 12.0] {
        assert_eq!(
            microphone_delay_taps(delay),
            (delay as usize, [1.0, 0.0, 0.0, 0.0])
        );
    }
    // Never wrap a negative centered-kernel start around the pressure ring.
    assert_eq!(microphone_delay_taps(0.25), (0, [0.75, 0.25, 0.0, 0.0]));
}

#[test]
fn g1_complex_propagation_error_is_reduced_across_fractional_paths() {
    // exp(-i omega tau) is the independent exact retarded-pressure oracle.
    // An amplitude-only test misses errors that move cancellation notches.
    for omega in [0.1 * PI, 0.2 * PI, 0.4 * PI] {
        for fraction in [0.125, 0.25, 0.5, 0.75, 0.875] {
            let delay = 12.0 + fraction;
            let phase = omega * delay;
            let exact = [phase.cos(), -phase.sin()];
            let (start, weights) = microphone_delay_taps(delay);
            let cubic_error = complex_error(response(start, &weights, omega), exact);
            let linear_error =
                complex_error(response(12, &[1.0 - fraction, fraction], omega), exact);
            assert!(
                cubic_error < 0.3 * linear_error,
                "omega={omega}, fraction={fraction}: cubic={cubic_error}, linear={linear_error}"
            );
        }
    }
}

#[test]
fn g3_export_retains_every_signed_coefficient_and_the_physical_flight_time() {
    let mut ins = Instrument::new(48_000.0);
    let mut body = Body::point("piston", 1.0);
    body.radiates = true;
    body.h = 0.01;
    let delay = 12.25;
    let distance = delay * ins.gas.sound_speed / ins.rate;
    body.origin[2] = -distance;
    ins.add(body);
    ins.mics.push(Mic::new([0.0; 3]));
    let image = ins.image();
    let taps = &image.mics[0].taps;
    assert_eq!(taps.len(), 2);
    assert!(taps[0].3 < 0.0 && taps[1].4 < 0.0);
    let gain = ins.gas.density * 0.01_f64.powi(2) * ins.rate.powi(2) / (2.0 * PI * distance);
    let total = taps.iter().map(|t| t.3 + t.4).sum::<f64>();
    let moment = taps
        .iter()
        .map(|t| t.2 as f64 * t.3 + (t.2 + 1) as f64 * t.4)
        .sum::<f64>();
    assert!((total / gain - 1.0).abs() < 1e-14);
    assert!((moment / total - delay).abs() < 1e-13);
    assert!(image.mics[0].ring_len > taps.last().unwrap().2 + 1);
}

#[test]
fn g3_live_pressure_improves_the_null_from_two_half_wavelength_separated_paths() {
    let mut ins = Instrument::new(48_000.0);
    // At 9.6 kHz, 2.5 samples of extra flight are exactly half a wavelength.
    // Different fractions expose the old filter's error; equal fractions
    // would cancel even with the old interpolation error.
    let delays = [10.25_f64, 12.75];
    let omega = 0.4 * PI;
    let distances = delays.map(|d| d * ins.gas.sound_speed / ins.rate);
    for distance in distances {
        let mut body = Body::point("driven piston", 1.0);
        body.radiates = true;
        body.h = 0.01;
        body.origin[2] = -distance;
        let body = ins.add(body);
        ins.drives.push(Drive {
            body,
            node: 0,
            knots: vec![(0.0, 0.0)],
        });
    }
    ins.mics.push(Mic::new([0.0; 3]));
    let gain = ins.gas.density * 0.01_f64.powi(2) / (2.0 * PI);
    let mut sample = [0.0];
    let mut cubic_energy = 0.0;
    let mut linear_energy = 0.0;
    for n in 0..1032 {
        for (drive, distance) in ins.drives.iter_mut().zip(distances) {
            // Equal pressure amplitudes after each path's 1/r spreading.
            drive.knots[0].1 = distance * (omega * n as f64).cos();
        }
        ins.step(&mut sample);
        if n < 32 {
            continue;
        }
        let exact = delays
            .iter()
            .map(|d| (omega * (n as f64 - d)).cos())
            .sum::<f64>();
        assert!(exact.abs() < 2e-12);
        let linear = delays
            .iter()
            .map(|d| {
                let whole = d.floor();
                let fraction = d - whole;
                (1.0 - fraction) * (omega * (n as f64 - whole)).cos()
                    + fraction * (omega * (n as f64 - whole - 1.0)).cos()
            })
            .sum::<f64>();
        cubic_energy += (sample[0] / gain).powi(2);
        linear_energy += linear * linear;
    }
    assert!(
        cubic_energy < 0.11 * linear_energy,
        "destructive-interference error: cubic={cubic_energy}, linear={linear_energy}"
    );
    assert!((cubic_energy / 1000.0).sqrt() < 0.014);
}
