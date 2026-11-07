use super::super::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};
use super::*;
use core::f64::consts::TAU;
use std::cell::Cell;

fn limits() -> WelchChannelTransferBudget {
    WelchChannelTransferBudget {
        maximum_output_channels: 8,
        maximum_complex_products: 1_000_000,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn observe<const D: usize>(
    rows: &[[f64; D]],
    n: usize,
    window: WelchWindow,
    detrend: WelchDetrend,
) -> WelchCrossSpectrum {
    WelchCrossSpectrum::estimate_interleaved(
        rows,
        0.002,
        WelchConfig {
            segment_samples: n,
            hop_samples: n / 2,
            window,
            detrend,
        },
        WelchBudget {
            maximum_channels: 8,
            maximum_samples_per_channel: 1024,
            maximum_transforms: 1024,
            maximum_workspace_bytes: 1_000_000,
        },
        || false,
    )
    .unwrap()
}

fn close(a: C64, b: C64, scale: f64) {
    assert!(
        (a.re - b.re).abs().max((a.im - b.im).abs()) <= 2e-11 * scale.max(1e-12),
        "{a:?} != {b:?}, scale {scale:e}"
    );
}

fn mul(a: C64, b: C64) -> C64 {
    C64::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

fn conjugate(a: C64) -> C64 {
    C64::new(a.re, -a.im)
}

#[test]
fn real_calibration_matches_record_mixing_with_original_windows_and_clock() {
    let gain = [[1.25, -0.75, 0.2], [0.0, -1.0, 0.5]];
    let rows: Vec<[f64; 3]> = (0..117)
        .map(|j| {
            let t = j as f64;
            [
                4.0 + (0.21 * t).cos(),
                -0.5 + 0.3 * (0.21 * t + 0.7).sin(),
                0.1 * ((j * 7 % 11) as f64 - 5.0),
            ]
        })
        .collect();
    let mixed: Vec<[f64; 2]> = rows
        .iter()
        .map(|x| core::array::from_fn(|o| (0..3).map(|i| gain[o][i] * x[i]).sum()))
        .collect();
    for n in [15, 32] {
        for window in [WelchWindow::Rectangular, WelchWindow::PeriodicHann] {
            for detrend in [WelchDetrend::None, WelchDetrend::SegmentMean] {
                let input = observe(&rows, n, window, detrend);
                let expected = observe(&mixed, n, window, detrend);
                let output = input
                    .transfer_channels(
                        2,
                        |_, _, o, i| Ok(C64::new(gain[o][i], 0.0)),
                        limits(),
                        || false,
                    )
                    .unwrap();
                assert_eq!(output.config(), input.config());
                assert_eq!(output.timestep_s(), input.timestep_s());
                assert_eq!(output.bin_width_hz(), input.bin_width_hz());
                assert_eq!(output.segment_count(), input.segment_count());
                assert_eq!(output.used_sample_count(), input.used_sample_count());
                assert_eq!(
                    output.discarded_tail_samples(),
                    input.discarded_tail_samples()
                );
                assert_eq!(output.channel_count(), 2);
                for bin in 0..output.bin_count() {
                    for a in 0..2 {
                        for b in 0..2 {
                            let oracle = expected.cross_spectrum(bin, a, b).unwrap();
                            close(
                                output.cross_spectrum(bin, a, b).unwrap(),
                                oracle,
                                oracle.re.abs().max(oracle.im.abs()),
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn complex_maps_preserve_full_cross_spectra_against_dense_covariance_oracle() {
    let rows: Vec<[f64; 3]> = (0..96)
        .map(|j| {
            let x = j as f64;
            [(0.4 * x).cos(), (0.4 * x + 0.8).sin(), (0.73 * x).cos()]
        })
        .collect();
    let input = observe(&rows, 31, WelchWindow::PeriodicHann, WelchDetrend::None);
    let gain = |bin: usize, output: usize, channel: usize| {
        let real = (output as f64 + 1.0) * (channel as f64 - 0.3);
        let imaginary = if bin == 0 {
            0.0
        } else {
            0.1 * bin as f64 * (output as f64 - channel as f64)
        };
        C64::new(real, imaginary)
    };
    let output = input
        .transfer_channels(2, |b, _, o, i| Ok(gain(b, o, i)), limits(), || false)
        .unwrap();
    for bin in 0..input.bin_count() {
        for a in 0..2 {
            for b in 0..2 {
                // Independent dense H S H* contraction; production instead
                // propagates each raw factor and never constructs S.
                let mut oracle = C64::default();
                for i in 0..3 {
                    for j in 0..3 {
                        let term = mul(
                            mul(gain(bin, a, i), input.cross_spectrum(bin, i, j).unwrap()),
                            conjugate(gain(bin, b, j)),
                        );
                        oracle.re += term.re;
                        oracle.im += term.im;
                    }
                }
                close(output.cross_spectrum(bin, a, b).unwrap(), oracle, 1.0);
            }
        }
    }
    let top = output.bin_count() - 1;
    let physical = output
        .peak_bin_hz(top, output.frequency_hz(top).unwrap())
        .unwrap();
    let mut peak_power = 0.0;
    for factor in physical.factors() {
        let p = physical.peak_factor(factor[0]).unwrap();
        peak_power += 0.5 * (p.re * p.re + p.im * p.im);
    }
    let density = output.cross_spectrum(top, 0, 0).unwrap().re;
    close(
        C64::new(peak_power, 0.0),
        C64::new(density * output.bin_width_hz(), 0.0),
        1.0,
    );
}

#[test]
fn phase_aligned_delayed_channels_construct_and_cancel_including_odd_top_bin() {
    const DELAY: usize = 3;
    for (n, tone) in [(31usize, 15usize), (32, 5)] {
        let rows: Vec<[f64; 2]> = (0..3 * n)
            .map(|j| {
                let angle = TAU * tone as f64 / n as f64;
                [
                    (angle * j as f64).cos(),
                    (angle * (j as f64 - DELAY as f64)).cos(),
                ]
            })
            .collect();
        let input = observe(&rows, n, WelchWindow::Rectangular, WelchDetrend::None);
        let output = input
            .transfer_channels(
                2,
                |bin, frequency, output, input| {
                    if input == 0 {
                        return Ok(C64::new(1.0, 0.0));
                    }
                    let angle = -TAU * frequency * DELAY as f64 * 0.002;
                    let (s, c) = angle.sin_cos();
                    let sign = if output == 0 { 1.0 } else { -1.0 };
                    Ok(C64::new(
                        sign * c,
                        if bin == 0 || (n % 2 == 0 && bin == n / 2) {
                            0.0
                        } else {
                            sign * s
                        },
                    ))
                },
                limits(),
                || false,
            )
            .unwrap();
        let original = input.cross_spectrum(tone, 0, 0).unwrap().re;
        close(
            output.cross_spectrum(tone, 0, 0).unwrap(),
            C64::new(4.0 * original, 0.0),
            original,
        );
        assert!(
            output.cross_spectrum(tone, 1, 1).unwrap().re < 1e-22 * original,
            "the original relative phase must survive delay compensation"
        );
    }
}

fn factor_fixture(channels: usize, segments: usize) -> WelchCrossSpectrum {
    WelchCrossSpectrum {
        config: WelchConfig {
            segment_samples: 4,
            hop_samples: 4,
            window: WelchWindow::Rectangular,
            detrend: WelchDetrend::None,
        },
        timestep_s: 0.25,
        bin_width_hz: 1.0,
        channels,
        segments,
        input_samples: 4 * segments,
        used_samples: 4 * segments,
        factors: vec![C64::default(); 3 * channels * segments],
    }
}

#[test]
fn weak_fields_and_product_residuals_survive_before_any_power_is_formed() {
    let mut input = factor_fixture(3, 2);
    input.factors[6..12].copy_from_slice(&[
        C64::new(1e100, 0.0),
        C64::new(1.0, 0.0),
        C64::new(-1e100, 0.0),
        C64::new(0.0, 1e100),
        C64::new(0.0, 1.0),
        C64::new(0.0, -1e100),
    ]);
    let output = input
        .transfer_channels(1, |_, _, _, _| Ok(C64::new(1.0, 0.0)), limits(), || false)
        .unwrap();
    assert_eq!(output.factors[2], C64::new(1.0, 0.0));
    assert_eq!(output.factors[3], C64::new(0.0, 1.0));
    assert_eq!(output.cross_spectrum(1, 0, 0).unwrap(), C64::new(2.0, 0.0));

    let eps = 2.0_f64.powi(-27);
    input.factors[6..12].copy_from_slice(&[
        C64::new(1.0 - eps, 0.0),
        C64::new(1.0, 0.0),
        C64::default(),
        C64::new(0.0, 1.0 - eps),
        C64::new(0.0, 1.0),
        C64::default(),
    ]);
    let output = input
        .transfer_channels(
            1,
            |_, _, _, i| Ok(C64::new([1.0 + eps, -1.0, 0.0][i], 0.0)),
            limits(),
            || false,
        )
        .unwrap();
    let residual = -2.0_f64.powi(-54);
    assert_eq!(output.factors[2], C64::new(residual, 0.0));
    assert_eq!(output.factors[3], C64::new(0.0, residual));
}

#[test]
fn complete_transfer_is_preflighted_and_coefficients_are_cached_per_frequency() {
    let input = factor_fixture(3, 2);
    let needs = input.channel_transfer_requirements(2).unwrap();
    assert_eq!(needs.coefficient_values, 18);
    assert_eq!(needs.complex_products, 36);
    assert_eq!(needs.output_factors, 12);
    assert_eq!(needs.workspace_bytes, 18 * core::mem::size_of::<C64>());
    let exact = WelchChannelTransferBudget {
        maximum_output_channels: 2,
        maximum_complex_products: needs.complex_products,
        maximum_workspace_bytes: needs.workspace_bytes,
    };
    for denied in [
        WelchChannelTransferBudget {
            maximum_output_channels: 1,
            ..exact
        },
        WelchChannelTransferBudget {
            maximum_complex_products: needs.complex_products - 1,
            ..exact
        },
        WelchChannelTransferBudget {
            maximum_workspace_bytes: needs.workspace_bytes - 1,
            ..exact
        },
    ] {
        assert!(matches!(
            input.transfer_channels(
                2,
                |_, _, _, _| panic!("budget admission must precede coefficient callbacks"),
                denied,
                || false,
            ),
            Err(WelchError::BudgetExceeded(_))
        ));
    }
    let calls = Cell::new(0usize);
    input
        .transfer_channels(
            2,
            |_, _, _, _| {
                calls.set(calls.get() + 1);
                Ok(C64::new(1.0, 0.0))
            },
            exact,
            || false,
        )
        .unwrap();
    assert_eq!(calls.get(), needs.coefficient_values);
    assert!(input.channel_transfer_requirements(0).is_err());
    assert!(matches!(
        input.channel_transfer_requirements(usize::MAX),
        Err(WelchError::BudgetExceeded(_))
    ));
}

#[test]
fn invalid_endpoints_callbacks_and_mid_contraction_cancellation_are_transactional() {
    let mut input = factor_fixture(140, 5);
    input.factors.fill(C64::new(1.0, 0.0));
    let original = input.factors.clone();
    let polls = Cell::new(0usize);
    assert!(matches!(
        input.transfer_channels(
            8,
            |_, _, _, _| Ok(C64::new(0.25, 0.0)),
            limits(),
            || {
                polls.set(polls.get() + 1);
                polls.get() == 4
            },
        ),
        Err(WelchError::Cancelled)
    ));
    assert_eq!(input.factors, original);
    let first = input
        .transfer_channels(8, |_, _, _, _| Ok(C64::new(0.25, 0.0)), limits(), || false)
        .unwrap();
    let replay = input
        .transfer_channels(8, |_, _, _, _| Ok(C64::new(0.25, 0.0)), limits(), || false)
        .unwrap();
    assert_eq!(first.factors, replay.factors);
    assert!(first.factors.iter().all(|v| *v == C64::new(35.0, 0.0)));
    for bin in [0, 2] {
        assert!(matches!(
            input.transfer_channels(
                1,
                |b, _, _, _| Ok(C64::new(1.0, if b == bin { 0.1 } else { 0.0 })),
                limits(),
                || false,
            ),
            Err(WelchError::Invalid(_))
        ));
    }
    assert!(matches!(
        input.transfer_channels(
            1,
            |_, _, _, _| Ok(C64::new(f64::NAN, 0.0)),
            limits(),
            || false
        ),
        Err(WelchError::NonFinite)
    ));
    assert!(matches!(
        input.transfer_channels(
            1,
            |_, _, _, _| Ok(C64::new(f64::MAX, 0.0)),
            limits(),
            || false
        ),
        Err(WelchError::Unrepresentable)
    ));
    let error = WelchError::Invalid("physical transfer is unavailable at this frequency");
    let actual = input.transfer_channels(
        1,
        |b, _, _, _| {
            if b == 1 {
                Err(error.clone())
            } else {
                Ok(C64::new(1.0, 0.0))
            }
        },
        limits(),
        || false,
    );
    assert_eq!(actual.unwrap_err(), error);
    assert_eq!(input.factors, original);
}
