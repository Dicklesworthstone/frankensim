use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn close(actual: Dd, expected: f64) {
    assert!(
        (actual.to_f64() - expected).abs() <= 5e-10 * expected.abs().max(1e-300),
        "actual={actual:?}, expected={expected:e}"
    );
}

#[test]
fn triangular_duct_null_keeps_real_means_and_signed_port_work() {
    let duration = 1e-12_f64;
    let window =
        HarmonicWindow::with_profile(-0.5 * duration, duration, HarmonicWindowProfile::Triangular)
            .unwrap();
    for proper in [false, true] {
        let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
        let frequencies = [carrier, carrier + 1.0, carrier + 2.0];
        let fields = |i: usize, _| {
            let p = Dd::from_f64([1.0, -2.0, 1.0][i]);
            [p, Dd::ZERO, Dd::from_f64(-2.0) * p, Dd::ZERO]
        };
        let mut sums = [[Dd::ZERO; 5]];
        assert!(
            accumulate(
                window,
                frequencies.into_iter(),
                fields,
                proper,
                !proper,
                &mut sums,
                &mut || false
            )
            .unwrap()
        );
        let mean = if proper {
            0.0
        } else {
            -duration.powi(2) / 24.0
        };
        let square = duration.powi(4) / if proper { 480.0 } else { 240.0 };
        for (actual, expected) in
            sums[0]
                .into_iter()
                .zip([mean, -2.0 * mean, square, 4.0 * square, -2.0 * square])
        {
            close(actual, expected);
        }
    }
}

#[test]
fn duct_profile_fallback_and_cancellation_leave_the_pair_owner_available() {
    let window =
        HarmonicWindow::with_profile(0.0, 10.0, HarmonicWindowProfile::Triangular).unwrap();
    let original = [[Dd::ONE; 5]];
    let mut sums = original;
    assert!(
        !accumulate(
            window,
            [1.0, 2.0, 3.0].into_iter(),
            |_, _| [Dd::ONE; 4],
            false,
            true,
            &mut sums,
            &mut || false
        )
        .unwrap()
    );
    assert_eq!(sums, original);
    let short =
        HarmonicWindow::with_profile(-0.05, 0.1, HarmonicWindowProfile::Triangular).unwrap();
    assert!(matches!(
        accumulate(
            short,
            [1.0, 2.0, 3.0].into_iter(),
            |_, _| [Dd::ONE; 4],
            false,
            true,
            &mut sums,
            &mut || true
        ),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(sums, original);
    let frequencies = [1.0; 33];
    assert!(
        !accumulate(
            short,
            frequencies.into_iter(),
            |_, _| [Dd::ONE; 4],
            false,
            true,
            &mut sums,
            &mut || false
        )
        .unwrap()
    );
    assert_eq!(sums, original);
}

#[test]
fn two_tone_duct_pressure_and_signed_flow_keep_real_and_proper_null_power() {
    let duration = 1e-12_f64;
    for proper in [false, true] {
        let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let triangle = profile == HarmonicWindowProfile::Triangular;
            let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
            let mut sums = [[Dd::ZERO; 5]];
            assert!(
                accumulate(
                    window,
                    [carrier, carrier + 1.0].into_iter(),
                    |i, _| {
                        let amplitude = Dd::from_f64(if i == 0 { 1.0 } else { -1.0 });
                        [
                            amplitude,
                            Dd::ZERO,
                            -Dd::from_f64(2.0) * amplitude,
                            Dd::ZERO,
                        ]
                    },
                    proper,
                    !proper,
                    &mut sums,
                    &mut || false
                )
                .unwrap()
            );
            let mean = if proper {
                0.0
            } else {
                duration.powi(2) / if triangle { 16.0 } else { 8.0 }
            };
            let square = if proper {
                duration.powi(2) / if triangle { 48.0 } else { 24.0 }
            } else {
                duration.powi(4) * if triangle { 3.0 / 320.0 } else { 9.0 / 320.0 }
            };
            for (actual, expected) in
                sums[0]
                    .into_iter()
                    .zip([mean, -2.0 * mean, square, 4.0 * square, -2.0 * square])
            {
                close(actual, expected);
            }
            let before = sums;
            assert!(
                !accumulate(
                    window,
                    [carrier].into_iter(),
                    |_, _| [Dd::ONE; 4],
                    proper,
                    !proper,
                    &mut sums,
                    &mut || false
                )
                .unwrap()
            );
            assert_eq!(sums, before);
        }
    }
}
