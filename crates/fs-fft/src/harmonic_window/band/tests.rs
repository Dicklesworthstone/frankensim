use super::*;
use crate::harmonic_window::HarmonicWindowProfile;

mod physical;
mod profile;

fn field(real: f64, imaginary: f64) -> [Dd; 2] {
    [Dd::from_f64(real), Dd::from_f64(imaginary)]
}

fn close(actual: Dd, expected: f64, relative: f64) {
    assert!(
        (actual.to_f64() - expected).abs() <= expected.abs() * relative,
        "actual={actual:?}, expected={expected:e}"
    );
}

#[test]
fn g1_three_tone_null_retains_positive_fourth_order_beat_energy() {
    let frequencies = [5e14, 5e14 + 1.0, 5e14 + 2.0];
    let fields = [field(1.0, 0.0), field(-2.0, 0.0), field(1.0, 0.0)];
    for duration in [1e-6, 1e-9, 1e-12] {
        let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
        let band = window
            .prepare_uniform_band_hz(&frequencies)
            .unwrap()
            .unwrap();
        // |1-2*exp(-i*w*t)+exp(-2*i*w*t)|² = 16*sin(w*t/2)^4.
        // Its leading centered-window integral is pi^4*T^4/5; the next
        // relative term is O(T²), below 2e-12 throughout this fixture.
        let expected = core::f64::consts::PI.powi(4) * duration.powi(4) / 5.0;
        let power = band.mean_square_parts_dd(&fields).unwrap();
        assert!(power.to_f64() > 0.0);
        close(power, expected, 3e-11);
        assert_eq!(
            band.mean_cross_parts_dd(&fields, &fields).unwrap(),
            [power, Dd::ZERO]
        );
    }
}

#[test]
fn g1_band_factor_matches_independent_quadrature_at_admission_edge() {
    let frequencies = [5e14, 5e14 + 0.75, 5e14 + 2.0];
    let duration = 0.49 / core::f64::consts::PI;
    let start = 0.137;
    let fields = [field(0.3, -0.7), field(-1.1, 0.2), field(0.8, 0.9)];
    let other = [field(-0.2, 0.3), field(0.4, -0.5), field(0.6, 0.7)];
    let window = HarmonicWindow::new(start, duration).unwrap();
    let band = window
        .prepare_uniform_band_hz(&frequencies)
        .unwrap()
        .unwrap();
    assert_eq!(band.mode_count(), MAX_UNIFORM_HARMONIC_BAND_MODES);
    let evaluate = |fields: &[[Dd; 2]; 3], time: f64| {
        fields
            .iter()
            .zip(frequencies)
            .fold([0.0; 2], |mut sum, (amplitude, frequency)| {
                let phase = -core::f64::consts::TAU * (frequency - frequencies[0]) * time;
                let (s, c) = phase.sin_cos();
                sum[0] += amplitude[0].to_f64() * c - amplitude[1].to_f64() * s;
                sum[1] += amplitude[0].to_f64() * s + amplitude[1].to_f64() * c;
                sum
            })
    };
    let mut expected = [0.0; 3];
    let steps = 4096;
    for i in 0..=steps {
        let t = start + duration * f64::from(i) / f64::from(steps);
        let a = evaluate(&fields, t);
        let b = evaluate(&other, t);
        let weight = if i == 0 || i == steps {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        expected[0] += weight * (a[0] * a[0] + a[1] * a[1]);
        expected[1] += weight * (a[0] * b[0] + a[1] * b[1]);
        expected[2] += weight * (a[1] * b[0] - a[0] * b[1]);
    }
    expected = expected.map(|x| x / (3.0 * f64::from(steps)));
    close(
        band.mean_square_parts_dd(&fields).unwrap(),
        expected[0],
        2e-13,
    );
    let cross = band.mean_cross_parts_dd(&fields, &other).unwrap();
    close(cross[0], expected[1], 2e-13);
    close(cross[1], expected[2], 2e-13);
    let reverse = band.mean_cross_parts_dd(&other, &fields).unwrap();
    assert_eq!(cross, [reverse[0], -reverse[1]]);
    let mut a = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
    let mut b = a;
    band.project_parts_dd(&fields, &mut a[..band.mode_count()])
        .unwrap();
    band.project_parts_dd(&other, &mut b[..band.mode_count()])
        .unwrap();
    let mut modal_cross = [Dd::ZERO; 2];
    let mut modal_power = Dd::ZERO;
    for mode in 0..band.mode_count() {
        let product = band.mode_cross_parts_dd(mode, a[mode], b[mode]).unwrap();
        modal_cross[0] = modal_cross[0] + product[0];
        modal_cross[1] = modal_cross[1] + product[1];
        let norm = band.mode_cross_parts_dd(mode, a[mode], a[mode]).unwrap();
        assert_eq!(norm[1], Dd::ZERO);
        modal_power = modal_power + norm[0];
    }
    assert_eq!(cross, modal_cross);
    assert_eq!(band.mean_square_parts_dd(&fields).unwrap(), modal_power);
}

#[test]
fn g3_absolute_carrier_does_not_change_envelope_clock_or_work() {
    let window = HarmonicWindow::new(1e15 + 0.25, 0.001).unwrap();
    let low = [0.0, 0.25, 1.0];
    let high = low.map(|x| 5e14 + x);
    let a = window.prepare_uniform_band_hz(&low).unwrap().unwrap();
    let b = window.prepare_uniform_band_hz(&high).unwrap().unwrap();
    let fields = [field(0.3, 0.7), field(-1.0, 0.5), field(-0.7, 0.2)];
    assert_eq!(
        a.mean_square_parts_dd(&fields).unwrap(),
        b.mean_square_parts_dd(&fields).unwrap()
    );
    assert_eq!(a.power_work(), 4 * a.mode_count());
    assert_eq!(a.cross_work(), 7 * a.mode_count());
    assert_eq!(a.coherency_work(), 10 * a.mode_count());
}

#[test]
fn g0_bounded_projection_retains_large_nulls_and_weighted_subnormal_power() {
    let duration = 1e-9;
    let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
    let band = window
        .prepare_uniform_band_hz(&[5e14, 5e14 + 1.0, 5e14 + 2.0])
        .unwrap()
        .unwrap();
    let base = core::f64::consts::PI.powi(4) * duration.powi(4) / 5.0;
    for amplitude in [1e100, 1e154] {
        let fields = [
            field(amplitude, 0.0),
            field(-2.0 * amplitude, 0.0),
            field(amplitude, 0.0),
        ];
        close(
            band.mean_square_parts_dd(&fields).unwrap(),
            base * amplitude * amplitude,
            3e-12,
        );
    }
    let amplitude = 1e-144;
    let fields = [
        field(amplitude, 0.0),
        field(-2.0 * amplitude, 0.0),
        field(amplitude, 0.0),
    ];
    let expected = base * amplitude * amplitude;
    let actual = band.mean_square_parts_dd(&fields).unwrap().to_f64();
    assert!(actual > 0.0);
    assert!(
        (actual - expected).abs() <= f64::from_bits(1),
        "actual={actual:e}, expected={expected:e}"
    );

    let boundary = 0.5 / core::f64::consts::PI;
    let inside = HarmonicWindow::new(0.0, boundary.next_down()).unwrap();
    let outside = HarmonicWindow::new(0.0, boundary.next_up()).unwrap();
    assert_eq!(
        inside.uniform_band_mode_count_hz(&[0.0, 1.0, 2.0]).unwrap(),
        Some(17)
    );
    assert_eq!(
        outside
            .uniform_band_mode_count_hz(&[0.0, 1.0, 2.0])
            .unwrap(),
        None
    );
}

#[test]
fn g0_unsupported_bands_fall_back_and_bad_fields_refuse() {
    let window = HarmonicWindow::new(0.0, 0.1).unwrap();
    for frequencies in [&[1.0, 2.0][..], &[0.0, 1.0, 10.0][..], &[1.0; 65][..]] {
        assert!(
            window
                .prepare_uniform_band_hz(frequencies)
                .unwrap()
                .is_none()
        );
    }
    let shaped = HarmonicWindow::with_profile(0.0, 0.1, HarmonicWindowProfile::Triangular).unwrap();
    assert!(
        shaped
            .prepare_uniform_band_hz(&[1.0, 2.0, 3.0])
            .unwrap()
            .is_none()
    );
    assert_eq!(
        window.uniform_band_mode_count_hz(&[1.0, 2.0, f64::NAN]),
        Err(WindowError::InvalidFrequency)
    );
    let band = window
        .prepare_uniform_band_hz(&[1.0, 2.0, 3.0])
        .unwrap()
        .unwrap();
    assert_eq!(
        band.mean_square_parts_dd(&[field(1.0, 0.0); 2]),
        Err(WindowError::InvalidBandShape)
    );
    let mut fields = [field(1.0, 0.0); 3];
    fields[1][1] = Dd::from_f64(f64::NAN);
    assert_eq!(
        band.mean_square_parts_dd(&fields),
        Err(WindowError::NonFiniteAmplitude)
    );
    assert_eq!(
        band.mode_cross_parts_dd(band.mode_count(), field(1.0, 0.0), field(1.0, 0.0)),
        Err(WindowError::InvalidBandShape)
    );
    assert_eq!(
        band.mode_cross_parts_dd(0, field(f64::NAN, 0.0), field(1.0, 0.0)),
        Err(WindowError::NonFiniteAmplitude)
    );
}
