use super::*;
use core::f64::consts::TAU;

fn budget() -> WelchBudget {
    WelchBudget {
        maximum_channels: 8,
        maximum_samples_per_channel: 1024,
        maximum_transforms: 128,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn config(window: WelchWindow, detrend: WelchDetrend) -> WelchConfig {
    WelchConfig {
        segment_samples: 32,
        hop_samples: 16,
        window,
        detrend,
    }
}

fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 4e-12 * scale,
        "{a:e} != {b:e}, scale {scale:e}"
    );
}

// G1: independent direct DFT of windowed real records, with the CSD formed
// directly instead of using the production factor/cross-spectrum routines.
fn oracle(channels: &[&[f64]], dt: f64, cfg: WelchConfig, bin: usize, a: usize, b: usize) -> C64 {
    let n = cfg.segment_samples;
    let count = 1 + (channels[0].len() - n) / cfg.hop_samples;
    let mut total = C64::default();
    for segment in 0..count {
        let start = segment * cfg.hop_samples;
        let transform = |channel: usize| {
            let x = &channels[channel][start..start + n];
            let mean = if cfg.detrend == WelchDetrend::SegmentMean {
                x.iter().sum::<f64>() / n as f64
            } else {
                0.0
            };
            let mut out = C64::default();
            let mut energy = 0.0;
            for (j, value) in x.iter().enumerate() {
                let w = if cfg.window == WelchWindow::Rectangular {
                    1.0
                } else {
                    (1.0 - (TAU * j as f64 / n as f64).cos()) * 0.5
                };
                let angle = TAU * bin as f64 * j as f64 / n as f64;
                out.re += (value - mean) * w * angle.cos();
                out.im -= (value - mean) * w * angle.sin();
                energy += w * w;
            }
            (out, energy)
        };
        let (x, energy) = transform(a);
        let (y, _) = transform(b);
        let sided = if bin == 0 || bin == n / 2 { 1.0 } else { 2.0 };
        let gain = sided * dt / (count as f64 * energy);
        // conj(X_a) * X_b matches the physical exp(-i omega t) convention.
        total.re += gain * (x.re * y.re + x.im * y.im);
        total.im += gain * (x.re * y.im - x.im * y.re);
    }
    total
}

#[test]
fn g1_joint_welch_matches_direct_dft_for_complex_cross_phase() {
    let rows: Vec<[f64; 3]> = (0..85)
        .map(|j| {
            let t = j as f64;
            [
                0.8 + (0.41 * t).cos() + 0.2 * (0.77 * t).sin(),
                -0.2 + 0.4 * (0.41 * t + 0.63).cos() - 0.1 * (0.77 * t).sin(),
                0.05 * ((j * 7 % 11) as f64 - 5.0),
            ]
        })
        .collect();
    let columns: [Vec<f64>; 3] = core::array::from_fn(|c| rows.iter().map(|r| r[c]).collect());
    let views: [&[f64]; 3] = core::array::from_fn(|c| columns[c].as_slice());
    for window in [WelchWindow::Rectangular, WelchWindow::PeriodicHann] {
        for detrend in [WelchDetrend::None, WelchDetrend::SegmentMean] {
            let cfg = config(window, detrend);
            let spectrum =
                WelchCrossSpectrum::estimate(&views, 0.002, cfg, budget(), || false).unwrap();
            let interleaved =
                WelchCrossSpectrum::estimate_interleaved(&rows, 0.002, cfg, budget(), || false)
                    .unwrap();
            assert_eq!(spectrum.factors, interleaved.factors);
            assert_eq!(spectrum.segment_count(), 4);
            assert_eq!(spectrum.used_sample_count(), 80);
            assert_eq!(spectrum.discarded_tail_samples(), 5);
            for bin in 0..spectrum.bin_count() {
                for a in 0..3 {
                    for b in 0..3 {
                        let actual = spectrum.cross_spectrum(bin, a, b).unwrap();
                        let expected = oracle(&views, 0.002, cfg, bin, a, b);
                        close(actual.re, expected.re, 0.02);
                        close(actual.im, expected.im, 0.02);
                        let transpose = spectrum.cross_spectrum(bin, b, a).unwrap();
                        close(actual.re, transpose.re, 0.02);
                        close(actual.im, -transpose.im, 0.02);
                    }
                }
            }
        }
    }
}

#[test]
fn g1_parseval_keeps_dc_nyquist_and_the_window_energy_normalization() {
    let n = 32;
    let samples: Vec<f64> = (0..80)
        .map(|j| {
            2.0 + 0.7 * (TAU * 4.0 * j as f64 / n as f64).cos()
                + if j % 2 == 0 { 0.3 } else { -0.3 }
        })
        .collect();
    for window in [WelchWindow::Rectangular, WelchWindow::PeriodicHann] {
        for detrend in [WelchDetrend::None, WelchDetrend::SegmentMean] {
            let cfg = config(window, detrend);
            let s =
                WelchCrossSpectrum::estimate(&[&samples], 0.001, cfg, budget(), || false).unwrap();
            let integral: f64 = (0..s.bin_count())
                .map(|k| s.cross_spectrum(k, 0, 0).unwrap().re * s.bin_width_hz())
                .sum();
            let mut expected = 0.0;
            for r in 0..s.segment_count() {
                let x = &samples[r * cfg.hop_samples..r * cfg.hop_samples + n];
                let mean = if detrend == WelchDetrend::SegmentMean {
                    x.iter().sum::<f64>() / n as f64
                } else {
                    0.0
                };
                let mut square = 0.0;
                let mut weight = 0.0;
                for (j, value) in x.iter().enumerate() {
                    let w = if window == WelchWindow::Rectangular {
                        1.0
                    } else {
                        0.5 - 0.5 * (TAU * j as f64 / n as f64).cos()
                    };
                    square += ((value - mean) * w).powi(2);
                    weight += w * w;
                }
                expected += square / weight / s.segment_count() as f64;
            }
            close(integral, expected, expected);
            assert!(s.cross_spectrum(n / 2, 0, 0).unwrap().re > 0.0);
        }
    }
    let tone: Vec<f64> = (0..32)
        .map(|j| (TAU * 4.0 * j as f64 / 32.0).cos())
        .collect();
    let s = WelchCrossSpectrum::estimate(
        &[&tone],
        0.001,
        config(WelchWindow::PeriodicHann, WelchDetrend::None),
        budget(),
        || false,
    )
    .unwrap();
    close(
        s.cross_spectrum(4, 0, 0).unwrap().re * s.bin_width_hz(),
        1.0 / 3.0,
        1.0,
    );
    for k in [3, 5] {
        close(
            s.cross_spectrum(k, 0, 0).unwrap().re * s.bin_width_hz(),
            1.0 / 12.0,
            1.0,
        );
    }
}

#[test]
fn g3_measured_factors_preserve_reinforcement_cancellation_and_partial_coherence() {
    let a: Vec<f64> = (0..64)
        .map(|j| (TAU * 4.0 * j as f64 / 32.0).cos())
        .collect();
    let cfg = WelchConfig {
        hop_samples: 32,
        ..config(WelchWindow::Rectangular, WelchDetrend::None)
    };
    for (mode, expected) in [(0, 2.0), (1, 0.0), (2, 1.0)] {
        let b: Vec<f64> = a
            .iter()
            .enumerate()
            .map(|(j, x)| {
                if mode == 1 || (mode == 2 && j >= 32) {
                    -*x
                } else {
                    *x
                }
            })
            .collect();
        let s = WelchCrossSpectrum::estimate(&[&a, &b], 0.001, cfg, budget(), || false).unwrap();
        // A physical linear observer adds amplitudes inside each factor, then
        // adds powers between factors. No dense CSD or inferred source phase.
        let output_power: f64 = s
            .factors(4)
            .unwrap()
            .map(|f| (f[0].re + f[1].re).powi(2) + (f[0].im + f[1].im).powi(2))
            .sum::<f64>()
            * s.bin_width_hz();
        close(output_power, expected, 2.0);
        let cross = s.cross_spectrum(4, 0, 1).unwrap();
        if mode == 2 {
            close(cross.re, 0.0, 0.05);
        }
    }
}

#[test]
fn g3_clock_rotation_and_extreme_input_scaling_keep_physical_cross_spectra() {
    for scale in [1e-100, 1.0, 1e100] {
        let mut baseline = None;
        for phase in [0.0, 0.83] {
            let rows: Vec<[f64; 2]> = (0..64)
                .map(|j| {
                    let t = TAU * 4.0 * j as f64 / 32.0 + phase;
                    [scale * t.cos(), scale * t.sin()]
                })
                .collect();
            let s = WelchCrossSpectrum::estimate_interleaved(
                &rows,
                0.001,
                config(WelchWindow::Rectangular, WelchDetrend::None),
                budget(),
                || false,
            )
            .unwrap();
            let cross = s.cross_spectrum(4, 0, 1).unwrap();
            let expected = -0.5 * scale * scale / s.bin_width_hz();
            close(cross.im, expected, expected.abs());
            if let Some(previous) = baseline {
                close(cross.im, previous, expected.abs());
            }
            baseline = Some(cross.im);
        }
    }
    // Segment means are formed after normalization, so a large finite constant
    // does not overflow its intermediate sum or create spurious fluctuation.
    let large = [[1e308]; 32];
    let s = WelchCrossSpectrum::estimate_interleaved(
        &large,
        0.001,
        config(WelchWindow::PeriodicHann, WelchDetrend::SegmentMean),
        budget(),
        || false,
    )
    .unwrap();
    assert!(s.factors.iter().all(|v| *v == C64::default()));
}

#[test]
fn g0_admission_budgets_cancel_and_invalid_queries_do_not_publish_partial_spectra() {
    let samples = [1.0; 65];
    let cfg = config(WelchWindow::PeriodicHann, WelchDetrend::SegmentMean);
    assert!(matches!(
        WelchCrossSpectrum::estimate(&[&samples, &samples[..64]], 0.001, cfg, budget(), || false),
        Err(WelchError::Invalid(_))
    ));
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(WelchCrossSpectrum::estimate(&[&samples], dt, cfg, budget(), || false).is_err());
    }
    for bad_cfg in [
        WelchConfig {
            segment_samples: 31,
            ..cfg
        },
        WelchConfig {
            hop_samples: 0,
            ..cfg
        },
        WelchConfig {
            hop_samples: 33,
            ..cfg
        },
    ] {
        assert!(
            WelchCrossSpectrum::estimate(&[&samples], 0.001, bad_cfg, budget(), || false).is_err()
        );
    }
    for limits in [
        WelchBudget {
            maximum_channels: 0,
            ..budget()
        },
        WelchBudget {
            maximum_samples_per_channel: 64,
            ..budget()
        },
        WelchBudget {
            maximum_transforms: 2,
            ..budget()
        },
        WelchBudget {
            maximum_workspace_bytes: 0,
            ..budget()
        },
    ] {
        assert!(matches!(
            WelchCrossSpectrum::estimate(&[&samples], 0.001, cfg, limits, || false),
            Err(WelchError::BudgetExceeded(_))
        ));
    }
    let mut calls = 0;
    assert!(matches!(
        WelchCrossSpectrum::estimate(&[&samples], 0.001, cfg, budget(), || {
            calls += 1;
            calls == 7
        }),
        Err(WelchError::Cancelled)
    ));
    let mut bad = samples;
    bad[64] = f64::NAN; // Even an incomplete, explicitly discarded tail is checked.
    assert!(matches!(
        WelchCrossSpectrum::estimate(&[&bad], 0.001, cfg, budget(), || false),
        Err(WelchError::NonFinite)
    ));
    let s = WelchCrossSpectrum::estimate(&[&samples], 0.001, cfg, budget(), || false).unwrap();
    assert_eq!(s.discarded_tail_samples(), 1);
    assert!(s.factors(s.bin_count()).is_err());
    assert!(s.cross_spectrum(0, 1, 0).is_err());
}

#[test]
fn g1_rescaling_preserves_representable_bins_after_tiny_intermediate_products() {
    let cfg = WelchConfig {
        segment_samples: 4,
        hop_samples: 4,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::None,
    };
    for (large, dt) in [(1e308, 1e-300), (1e300, 1e-100)] {
        // Exact DFT: X[1] = -2i regardless of the equal large samples.
        // Multiplying its normalized coefficient by sqrt(dt) first would
        // underflow, even though the restored physical factor is finite.
        let samples = [[large], [1.0], [large], [-1.0]];
        let spectrum =
            WelchCrossSpectrum::estimate_interleaved(&samples, dt, cfg, budget(), || false)
                .unwrap();
        let factor = spectrum.factors(1).unwrap().next().unwrap()[0];
        let expected = (2.0 * dt).sqrt();
        close(factor.re, 0.0, expected);
        close(factor.im, expected, expected);
        close(
            spectrum.cross_spectrum(1, 0, 0).unwrap().re,
            2.0 * dt,
            2.0 * dt,
        );
        assert!(
            spectrum
                .factors
                .iter()
                .all(|v| v.re.is_finite() && v.im.is_finite())
        );
    }
}

#[test]
fn g0_normalization_underflow_refuses_to_erase_a_representable_source_bin() {
    let cfg = WelchConfig {
        segment_samples: 4,
        hop_samples: 4,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::None,
    };
    // The true bin-1 factor is i*sqrt(2)*1e-80 and its PSD is 2e-160.
    // Dividing the small samples by 1e300 cannot retain them in f64;
    // publishing a zero factor would silently invent perfect cancellation.
    let samples = [[1e300], [1e-30], [1e300], [-1e-30]];
    assert!(matches!(
        WelchCrossSpectrum::estimate_interleaved(&samples, 1e-100, cfg, budget(), || false),
        Err(WelchError::Unrepresentable)
    ));
}

#[test]
fn g1_detrending_preserves_representable_fluctuations_beside_dc_offsets() {
    let cfg = WelchConfig {
        segment_samples: 4,
        hop_samples: 4,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::SegmentMean,
    };
    let d = f64::EPSILON;
    for offset in [0.0, 1.0, 1.5, 1.75, -1.5] {
        // All four authored values and the centered +/-d signal are exact.
        // F64 normalization before detrending changes the 1.5-offset PSD by
        // +26.5625% and the 1.75-offset PSD by -23.4375%.
        let samples = [offset + d, offset - d, offset + d, offset - d];
        let spectrum =
            WelchCrossSpectrum::estimate(&[&samples], 1.0, cfg, budget(), || false).unwrap();
        let expected = 4.0 * d * d;
        close(
            spectrum.cross_spectrum(2, 0, 0).unwrap().re,
            expected,
            expected,
        );
    }
    // Weak quadrature signals on DIFFERENT offsets retain both their channel
    // amplitudes and relative phase. X=2d, Y=-2id in the ordinary DFT; this
    // owner's negative-time one-sided cross spectrum is -2id² per Hz.
    let x = [1.5 + d, 1.5, 1.5 - d, 1.5];
    let y = [1.75, 1.75 + d, 1.75, 1.75 - d];
    let spectrum = WelchCrossSpectrum::estimate(&[&x, &y], 1.0, cfg, budget(), || false).unwrap();
    let expected = 2.0 * d * d;
    let cross = spectrum.cross_spectrum(1, 0, 1).unwrap();
    close(cross.re, 0.0, expected);
    close(cross.im, -expected, expected);
    for channel in 0..2 {
        close(
            spectrum.cross_spectrum(1, channel, channel).unwrap().re,
            expected,
            expected,
        );
    }
}
