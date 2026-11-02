use super::*;

fn limits(n: usize) -> HarmonicProductSpectrumLimits {
    HarmonicProductSpectrumLimits {
        maximum_products: n,
        maximum_generated_terms: 2 * n,
        maximum_workspace_bytes: 1_048_576,
    }
}

fn input(wa: f64, a: [f64; 2], wb: f64, b: [f64; 2]) -> AngularHarmonicProductInput {
    AngularHarmonicProductInput {
        left_angular_frequency_rad_s: wa,
        right_angular_frequency_rad_s: wb,
        left_amplitude: a.map(Dd::from_f64),
        right_amplitude: b.map(Dd::from_f64),
        law: HarmonicProductLaw::Real,
    }
}

fn prepare(inputs: Vec<AngularHarmonicProductInput>) -> PreparedHarmonicProductSpectrum {
    let count = inputs.len();
    PreparedHarmonicProductSpectrum::from_angular_frequencies(inputs, limits(count)).unwrap()
}

fn close(a: Dd, b: f64, relative: f64, scale: f64) {
    assert!(
        (a.to_f64() - b).abs() <= relative * scale,
        "{a:?} != {b:e} (scale {scale:e})"
    );
}

fn simpson(mut f: impl FnMut(f64) -> f64, duration: f64) -> f64 {
    let n = 12_000;
    let mut total = 0.0;
    for i in 0..=n {
        let u = duration * f64::from(i) / f64::from(n);
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        total += weight * f(u);
    }
    total * duration / (3.0 * f64::from(n))
}

#[test]
fn dc_response_preserves_the_short_decay_decrement_and_stiff_limit() {
    let spectrum = prepare(vec![input(0.0, [2.0, 0.0], 0.0, [1.0, 0.0])]);
    let duration = 0.4;
    for rate in [0.0, 1e-18, 1e-6, 0.5, 8.0, 2000.0] {
        let result = spectrum.causal_integral(-0.25, duration, rate).unwrap();
        let expected = if rate == 0.0 {
            2.0 * duration
        } else {
            2.0 * -(-rate * duration).exp_m1() / rate
        };
        close(result.integral_dd(), 0.8, 2e-15, 1.0);
        close(result.weighted_integral_dd(), expected, 2e-15, expected);
        close(result.decay_dd(), (-rate * duration).exp(), 2e-15, 1.0);
        close(
            result.decay_dd() + result.one_minus_decay_dd(),
            1.0,
            2e-30,
            1.0,
        );
    }
    let alpha = Dd::from_f64(1e-18);
    let short = spectrum.causal_integral(0.0, 1.0, 1e-18).unwrap();
    assert!(short.decay_dd().lo < 0.0);
    assert_eq!(short.decay_dd().hi, 1.0);
    close(
        (short.one_minus_decay_dd() - alpha) / (alpha * alpha),
        -0.5,
        3e-14,
        1.0,
    );
}

#[test]
fn causal_products_match_independent_real_time_quadrature() {
    let left = [(2.1, [0.7, -0.2]), (7.4, [-0.3, 0.5])];
    let right = [(-1.4, [0.9, 0.4]), (3.6, [0.2, -0.7])];
    let mut products = Vec::new();
    for &(wa, a) in &left {
        for &(wb, b) in &right {
            products.push(input(wa, a, wb, b));
        }
    }
    products.push(input(0.0, [0.4, 0.0], 0.0, [1.0, 0.0]));
    let spectrum = prepare(products);
    let start = -0.17;
    let real_signal = |lines: &[(f64, [f64; 2])], t: f64| -> f64 {
        lines
            .iter()
            .map(|(w, a)| a[0] * (w * t).cos() + a[1] * (w * t).sin())
            .sum()
    };
    let product = |s: f64| real_signal(&left, start + s) * real_signal(&right, start + s) + 0.4;
    for duration in [1e-5, 0.2, 1.4] {
        for rate in [0.0, 1e-12, 0.3, 4.0] {
            let result = spectrum.causal_integral(start, duration, rate).unwrap();
            let plain = simpson(product, duration);
            let weighted = simpson(|s| (-rate * (duration - s)).exp() * product(s), duration);
            close(
                result.integral_dd(),
                plain,
                5e-12,
                duration.max(plain.abs()),
            );
            close(
                result.weighted_integral_dd(),
                weighted,
                5e-12,
                duration.max(weighted.abs()),
            );
            if rate == 0.0 {
                assert_eq!(result.integral_dd(), result.weighted_integral_dd());
            }
        }
    }
}

#[test]
fn accepted_window_partitions_obey_the_causal_semigroup() {
    let spectrum = prepare(vec![
        input(3.0, [0.3, -0.7], 5.5, [-0.9, 0.1]),
        input(0.0, [1.2, 0.0], 0.0, [1.0, 0.0]),
    ]);
    let start = 0.375;
    let first_duration = 0.125;
    let second_duration = 0.5;
    for rate in [0.0, 1e-15, 0.8, 30.0] {
        let all = spectrum
            .causal_integral(start, first_duration + second_duration, rate)
            .unwrap();
        let first = spectrum
            .causal_integral(start, first_duration, rate)
            .unwrap();
        let second = spectrum
            .causal_integral(start + first_duration, second_duration, rate)
            .unwrap();
        close(
            first.integral_dd() + second.integral_dd(),
            all.integral_dd().to_f64(),
            3e-14,
            1.0,
        );
        close(
            first.weighted_integral_dd() * second.decay_dd() + second.weighted_integral_dd(),
            all.weighted_integral_dd().to_f64(),
            3e-14,
            1.0,
        );
        close(
            first.decay_dd() * second.decay_dd(),
            all.decay_dd().to_f64(),
            3e-14,
            1.0,
        );
    }
}

#[test]
fn short_weak_sine_energy_survives_causal_dc_carrier_cancellation() {
    let omega = 3e-5;
    let duration = 1e-4;
    let spectrum = prepare(vec![input(omega, [0.0, 1.0], omega, [0.0, 1.0])]);
    for rate in [0.0, 1e-5, 1.0, 100.0] {
        let result = spectrum.causal_integral(0.0, duration, rate).unwrap();
        let expected = simpson(
            |s| (omega * s).sin().powi(2) * (-rate * (duration - s)).exp(),
            duration,
        );
        assert!(result.weighted_integral_dd().to_f64() > 0.0);
        close(result.weighted_integral_dd(), expected, 4e-12, expected);
    }
}

#[test]
fn native_clock_survives_unrepresentable_endpoint_and_large_carrier_beats() {
    let start = 2.0_f64.powi(50);
    let duration = 0.03125;
    assert_eq!(start + duration, start);
    let omega = 2.0_f64.powi(-40);
    let spectrum = prepare(vec![input(omega, [0.7, 0.3], 0.0, [1.0, 0.0])]);
    let result = spectrum.causal_integral(start, duration, 0.4).unwrap();
    // The phase at the start is exactly 1024 rad. Integrate in local s;
    // forming start+s in one f64 would erase the entire accepted interval.
    let expected = simpson(
        |s| {
            let phase = 1024.0 + omega * s;
            (0.7 * phase.cos() + 0.3 * phase.sin()) * (-0.4 * (duration - s)).exp()
        },
        duration,
    );
    close(result.weighted_integral_dd(), expected, 3e-14, duration);

    let carrier = 2.0_f64.powi(50);
    let mut beat = input(carrier + 0.25, [0.7, -0.2], carrier, [0.3, 0.8]);
    beat.law = HarmonicProductLaw::ProperComplex;
    let spectrum = prepare(vec![beat]);
    assert_eq!(spectrum.tones().len(), 1);
    assert_eq!(
        spectrum.tones()[0].angular_frequency_dd(),
        Dd::from_f64(0.25)
    );
    let result = spectrum.causal_integral(0.2, 0.7, 0.9).unwrap();
    let expected = simpson(
        |s| {
            let phase = 0.25 * (0.2 + s);
            let real = 0.5 * (0.7 * 0.3 + (-0.2) * 0.8);
            let imaginary = 0.5 * ((-0.2) * 0.3 - 0.7 * 0.8);
            (real * phase.cos() + imaginary * phase.sin()) * (-0.9 * (0.7 - s)).exp()
        },
        0.7,
    );
    close(result.weighted_integral_dd(), expected, 3e-14, 1.0);
}

#[test]
fn zero_duration_darkness_invalid_domains_and_cancellation_are_explicit() {
    let spectrum = prepare(vec![input(1.0, [1.0, 0.0], 0.0, [1.0, 0.0])]);
    let zero = spectrum.causal_integral(1e300, 0.0, 1e300).unwrap();
    assert_eq!(zero.integral_dd(), Dd::ZERO);
    assert_eq!(zero.weighted_integral_dd(), Dd::ZERO);
    assert_eq!(zero.decay_dd(), Dd::ONE);
    assert_eq!(zero.one_minus_decay_dd(), Dd::ZERO);
    for (start, duration, rate) in [
        (f64::NAN, 1.0, 1.0),
        (0.0, -1.0, 1.0),
        (0.0, 1.0, -1.0),
        (0.0, f64::INFINITY, 1.0),
        (0.0, 1.0, f64::INFINITY),
    ] {
        assert_eq!(
            spectrum.causal_integral(start, duration, rate).unwrap_err(),
            ProductSpectrumError::Window(WindowError::InvalidWindow),
        );
    }
    assert_eq!(
        spectrum
            .causal_integral_cancellable(0.0, 1.0, 0.2, || true)
            .unwrap_err(),
        ProductSpectrumError::Cancelled,
    );
    let baseline = spectrum.causal_integral(0.0, 1.0, 0.2).unwrap();
    let mut calls = 0;
    assert_eq!(
        spectrum
            .causal_integral_cancellable(0.0, 1.0, 0.2, || {
                calls += 1;
                calls == 2
            })
            .unwrap_err(),
        ProductSpectrumError::Cancelled,
    );
    assert_eq!(
        spectrum
            .causal_integral(0.0, 1.0, 0.2)
            .unwrap()
            .weighted_integral_dd(),
        baseline.weighted_integral_dd(),
    );
    let dark = prepare(vec![]);
    assert_eq!(
        dark.causal_integral(0.0, 1.0, 0.2)
            .unwrap()
            .weighted_integral_dd(),
        Dd::ZERO
    );
}
