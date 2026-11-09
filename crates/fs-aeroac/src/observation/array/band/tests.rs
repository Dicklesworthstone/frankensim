use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn field<const N: usize>(
    frequencies: [f64; N],
    amplitudes: [[f64; 2]; N],
    gain: f64,
) -> PreparedHarmonicField<2> {
    PreparedHarmonicField {
        lines: frequencies
            .into_iter()
            .zip(amplitudes)
            .map(|(frequency, a)| {
                let component = |g| ComplexSum {
                    re: Dd::from_f64(gain * g * a[0]),
                    im: Dd::from_f64(gain * g * a[1]),
                };
                PreparedLine {
                    frequency: PreparedFrequency::Angular(frequency),
                    pressure: component(1.0),
                    velocity: [component(-2.0), component(0.5)],
                }
            })
            .collect(),
    }
}

#[test]
fn triangular_cross_receiver_null_retains_negative_correlation_and_weighted_readout() {
    let amplitudes = [[1.0, 0.0], [-2.0, 0.0], [1.0, 0.0]];
    let a = field([1.0, 2.0, 3.0], amplitudes, 1.0);
    let b = field([1.0, 2.0, 3.0], amplitudes, -0.5);
    let duration = 1e-12_f64;
    let window =
        HarmonicWindow::with_profile(-duration / 2.0, duration, HarmonicWindowProfile::Triangular)
            .unwrap();
    let mut sums = CrossSums::new();
    let mut covariance = CrossSums::new();
    assert!(
        accumulate(
            &a,
            &b,
            window,
            FactorKind::Deterministic,
            &mut sums,
            Some(&mut covariance),
            &mut || false
        )
        .unwrap()
    );
    let square = -0.5 * duration.powi(4) / 240.0;
    for (actual, expected) in [
        sums.pp,
        sums.pu[0],
        sums.up[0],
        sums.uu[0][1],
        covariance.pp,
    ]
    .into_iter()
    .zip([
        square,
        -2.0 * square,
        -2.0 * square,
        -square,
        -0.5 * duration.powi(4) / 576.0,
    ]) {
        assert!(
            (actual.to_f64() - expected).abs() < 5e-10 * expected.abs(),
            "{actual:?} != {expected:e}"
        );
    }
    let mut reversed = CrossSums::new();
    assert!(
        accumulate(
            &b,
            &a,
            window,
            FactorKind::Deterministic,
            &mut reversed,
            None,
            &mut || false
        )
        .unwrap()
    );
    assert_eq!(sums.pp, reversed.pp);
    assert_eq!(sums.pu, reversed.up);
}

#[test]
fn shaped_union_correlation_and_proper_average_keep_the_original_carrier() {
    let a = field(
        [8.0, 9.0, 10.0],
        [[0.7, 0.2], [-0.4, 0.3], [0.1, -0.6]],
        1.0,
    );
    let b = field(
        [9.0, 10.0, 11.0],
        [[-0.2, 0.1], [0.3, 0.7], [0.6, -0.2]],
        -0.5,
    );
    let window =
        HarmonicWindow::with_profile(0.03, 0.08, HarmonicWindowProfile::Triangular).unwrap();
    let mut sums = CrossSums::new();
    let mut covariance = CrossSums::new();
    assert!(
        accumulate(
            &a,
            &b,
            window,
            FactorKind::Proper,
            &mut sums,
            Some(&mut covariance),
            &mut || false
        )
        .unwrap()
    );
    let mut means = [[0.0; 2]; 2];
    let mut product = 0.0;
    let steps = 2048;
    for k in 0..=steps {
        let fraction = k as f64 / steps as f64;
        let time = 0.03 + 0.08 * fraction;
        let weight = if k == 0 || k == steps {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        } * 2.0
            * (1.0 - (2.0 * fraction - 1.0).abs())
            / (3.0 * steps as f64);
        let mut fields = [[0.0; 2]; 2];
        for (r, source) in [&a, &b].into_iter().enumerate() {
            for line in &source.lines {
                let (s, c) = (line.frequency.native() * time).sin_cos();
                let re = line.pressure.re.to_f64();
                let im = line.pressure.im.to_f64();
                fields[r][0] += re * c + im * s;
                fields[r][1] += im * c - re * s;
            }
            for component in 0..2 {
                means[r][component] += weight * fields[r][component];
            }
        }
        product += weight * 0.5 * (fields[0][0] * fields[1][0] + fields[0][1] * fields[1][1]);
    }
    let averaged = 0.5 * (means[0][0] * means[1][0] + means[0][1] * means[1][1]);
    assert!((sums.pp.to_f64() - product).abs() < 2e-12);
    assert!((covariance.pp.to_f64() - averaged).abs() < 2e-12);
    assert!((sums.pu[0].to_f64() + 2.0 * product).abs() < 4e-12);
    let before = sums.pp;
    let broad = HarmonicWindow::with_profile(0.0, 10.0, HarmonicWindowProfile::Triangular).unwrap();
    assert!(
        !accumulate(
            &a,
            &b,
            broad,
            FactorKind::Proper,
            &mut sums,
            None,
            &mut || false
        )
        .unwrap()
    );
    assert_eq!(sums.pp, before);
}

#[test]
fn two_tone_cross_receiver_gram_retains_real_and_proper_signed_power() {
    let duration = 1e-12_f64;
    for proper in [false, true] {
        let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
        let a = field([carrier, carrier + 1.0], [[1.0, 0.0], [-1.0, 0.0]], 1.0);
        let b = field([carrier, carrier + 1.0], [[1.0, 0.0], [-1.0, 0.0]], -0.5);
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let triangle = profile == HarmonicWindowProfile::Triangular;
            let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
            let mut sums = CrossSums::new();
            let mut covariance = CrossSums::new();
            assert!(
                accumulate(
                    &a,
                    &b,
                    window,
                    if proper {
                        FactorKind::Proper
                    } else {
                        FactorKind::Deterministic
                    },
                    &mut sums,
                    Some(&mut covariance),
                    &mut || false
                )
                .unwrap()
            );
            let square = -0.5
                * if proper {
                    duration.powi(2) / if triangle { 48.0 } else { 24.0 }
                } else {
                    duration.powi(4) * if triangle { 3.0 / 320.0 } else { 9.0 / 320.0 }
                };
            for (actual, expected) in [sums.pp, sums.pu[0], sums.up[0], sums.uu[0][1]]
                .into_iter()
                .zip([square, -2.0 * square, -2.0 * square, -square])
            {
                assert!(
                    (actual.to_f64() - expected).abs() < 5e-10 * expected.abs(),
                    "{actual:?} != {expected:e}"
                );
            }
            if !proper {
                let mean = duration.powi(2) / if triangle { 16.0 } else { 8.0 };
                assert!((covariance.pp.to_f64() + 0.5 * mean * mean).abs() < 5e-10 * mean * mean);
            }
        }
    }
}
