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
fn g1_borrowed_peak_bin_retains_phase_and_converts_bin_power_once() {
    let rows: Vec<[f64; 2]> = (0..64)
        .map(|j| {
            let phase = TAU * 3.0 * j as f64 / 32.0;
            [2.0 * phase.cos(), -3.0 * phase.sin()]
        })
        .collect();
    let cfg = WelchConfig {
        hop_samples: 32,
        ..config(WelchWindow::Rectangular, WelchDetrend::None)
    };
    let spectrum =
        WelchCrossSpectrum::estimate_interleaved(&rows, 0.003, cfg, budget(), || false).unwrap();
    let frequency = spectrum.frequency_hz(3).unwrap();
    let view = spectrum.peak_bin_hz(3, frequency).unwrap();
    assert_eq!(view.channel_count(), 2);
    assert_eq!(view.segment_count(), 2);
    assert_eq!(view.frequency_hz(), frequency);
    let original = spectrum.factors(3).unwrap().next().unwrap();
    assert!(core::ptr::eq(
        view.factors().next().unwrap().as_ptr(),
        original.as_ptr()
    ));
    let mut powers = [0.0; 2];
    for factor in view.factors() {
        let first = view.peak_factor(factor[0]).unwrap();
        let second = view.peak_factor(factor[1]).unwrap();
        close(first.re, core::f64::consts::SQRT_2, 2.0);
        close(first.im, 0.0, 2.0);
        close(second.re, 0.0, 3.0);
        close(second.im, -3.0 / core::f64::consts::SQRT_2, 3.0);
        powers[0] += first.norm_sq() * 0.5;
        powers[1] += second.norm_sq() * 0.5;
    }
    close(powers[0], 2.0, 2.0);
    close(powers[1], 4.5, 4.5);
}

#[test]
fn g0_peak_bin_admission_respects_exact_hz_and_original_angular_clock() {
    let spectrum = WelchCrossSpectrum::estimate(
        &[&[0.0; 32]],
        0.003,
        config(WelchWindow::Rectangular, WelchDetrend::None),
        budget(),
        || false,
    )
    .unwrap();
    let frequency = spectrum.frequency_hz(3).unwrap();
    let omega = TAU * frequency;
    // This actual bin loses one ULP if its authored angular clock is divided
    // back to Hz. Admission must compare clocks in the observer's units.
    assert_ne!((omega / TAU).to_bits(), frequency.to_bits());
    assert!(spectrum.peak_bin_rad_s(3, omega).is_ok());
    assert!(spectrum.peak_bin_hz(3, frequency).is_ok());
    for bin in [0, 16, 17, usize::MAX] {
        assert!(spectrum.peak_bin_hz(bin, frequency).is_err());
        assert!(spectrum.peak_bin_rad_s(bin, omega).is_err());
    }
    for bad in [
        0.0,
        -frequency,
        f64::NAN,
        f64::INFINITY,
        f64::from_bits(frequency.to_bits() + 1),
    ] {
        assert!(spectrum.peak_bin_hz(3, bad).is_err());
    }
    for bad in [
        0.0,
        -omega,
        f64::NAN,
        f64::INFINITY,
        f64::from_bits(omega.to_bits() + 1),
    ] {
        assert!(spectrum.peak_bin_rad_s(3, bad).is_err());
    }
}

#[test]
fn g0_peak_conversion_preserves_representable_extremes_and_refuses_range_loss() {
    let cfg = WelchConfig {
        segment_samples: 4,
        hop_samples: 4,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::None,
    };
    let zero_record = [0.0; 4];
    let at_clock =
        |dt| WelchCrossSpectrum::estimate(&[&zero_record], dt, cfg, budget(), || false).unwrap();
    let high = at_clock(2.0_f64.powi(-1001));
    let high = high.peak_bin_hz(1, high.frequency_hz(1).unwrap()).unwrap();
    let low = at_clock(2.0_f64.powi(999));
    let low = low.peak_bin_hz(1, low.frequency_hz(1).unwrap()).unwrap();
    let small = 2.0_f64.powi(-600);
    let large = 2.0_f64.powi(600);
    close(
        high.peak_factor(C64::new(small, -small)).unwrap().re,
        2.0_f64.powi(-100),
        2.0_f64.powi(-100),
    );
    close(
        low.peak_factor(C64::new(large, -large)).unwrap().im,
        -2.0_f64.powi(100),
        2.0_f64.powi(100),
    );
    assert_eq!(
        low.peak_factor(C64::new(2.0_f64.powi(-570), 0.0))
            .unwrap()
            .re,
        2.0_f64.powi(-1070)
    );
    assert!(matches!(
        high.peak_factor(C64::new(large, 0.0)),
        Err(WelchError::Unrepresentable)
    ));
    assert!(matches!(
        low.peak_factor(C64::new(0.0, small)),
        Err(WelchError::Unrepresentable)
    ));
    assert_eq!(low.peak_factor(C64::default()).unwrap(), C64::default());
    assert!(matches!(
        high.peak_factor(C64::new(f64::NAN, 0.0)),
        Err(WelchError::NonFinite)
    ));
    assert!(matches!(
        low.peak_factor(C64::new(0.0, f64::INFINITY)),
        Err(WelchError::NonFinite)
    ));
    let angular_overflow = at_clock(6.0e-309);
    let frequency = angular_overflow.frequency_hz(1).unwrap();
    assert!(angular_overflow.peak_bin_hz(1, frequency).is_ok());
    assert!(angular_overflow.peak_bin_rad_s(1, TAU * frequency).is_err());
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

fn same_stream_spectrum(actual: &WelchCrossSpectrum, expected: &WelchCrossSpectrum) {
    assert_eq!(actual.config(), expected.config());
    assert_eq!(
        actual.timestep_s().to_bits(),
        expected.timestep_s().to_bits()
    );
    assert_eq!(
        actual.bin_width_hz().to_bits(),
        expected.bin_width_hz().to_bits()
    );
    assert_eq!(actual.channel_count(), expected.channel_count());
    assert_eq!(actual.segment_count(), expected.segment_count());
    assert_eq!(actual.used_sample_count(), expected.used_sample_count());
    assert_eq!(
        actual.discarded_tail_samples(),
        expected.discarded_tail_samples()
    );
    assert_eq!(actual.factors.len(), expected.factors.len());
    for (a, b) in actual.factors.iter().zip(&expected.factors) {
        assert_eq!(
            (a.re.to_bits(), a.im.to_bits()),
            (b.re.to_bits(), b.im.to_bits())
        );
    }
}

#[test]
fn g3_stream_partitions_preserve_exact_batch_phase_overlap_windows_and_tails() {
    let rows: Vec<[f64; 3]> = (0..85)
        .map(|j| {
            let t = j as f64;
            [
                0.8 + (0.41 * t).cos(),
                -0.2 + (0.41 * t + 0.63).cos(),
                0.05 * ((j * 7 % 11) as f64 - 5.0),
            ]
        })
        .collect();
    let columns: [Vec<f64>; 3] = core::array::from_fn(|c| rows.iter().map(|r| r[c]).collect());
    let limits = WelchBudget {
        maximum_transforms: 256,
        ..budget()
    };
    for window in [WelchWindow::Rectangular, WelchWindow::PeriodicHann] {
        for detrend in [WelchDetrend::None, WelchDetrend::SegmentMean] {
            for hop_samples in [1, 5, 16, 32] {
                let cfg = WelchConfig {
                    hop_samples,
                    ..config(window, detrend)
                };
                let expected =
                    WelchCrossSpectrum::estimate_interleaved(&rows, 0.002, cfg, limits, || false)
                        .unwrap();
                for pattern in [
                    &[1][..],
                    &[7, 1, 23, 4][..],
                    &[31, 1, 32, 13, 8][..],
                    &[85][..],
                ] {
                    let mut stream =
                        WelchStream::new(3, rows.len(), 0.002, cfg, limits, || false).unwrap();
                    let (mut start, mut chunk) = (0, 0);
                    while start < rows.len() {
                        let end = (start + pattern[chunk % pattern.len()]).min(rows.len());
                        if chunk % 2 == 0 {
                            stream
                                .push_interleaved(&rows[start..end], || false)
                                .unwrap();
                        } else {
                            let slices: [&[f64]; 3] =
                                core::array::from_fn(|c| &columns[c][start..end]);
                            stream.push(&slices, || false).unwrap();
                        }
                        assert_eq!(stream.accepted_sample_count(), end);
                        assert_eq!(stream.remaining_sample_count(), rows.len() - end);
                        assert_eq!(
                            stream.completed_segment_count(),
                            if end < 32 {
                                0
                            } else {
                                1 + (end - 32) / hop_samples
                            }
                        );
                        start = end;
                        chunk += 1;
                    }
                    stream.push_interleaved::<3>(&[], || false).unwrap();
                    same_stream_spectrum(&stream.finish().unwrap(), &expected);
                }
            }
        }
    }
}

#[test]
fn g4_stream_cancelled_chunks_replay_with_different_data_without_stale_factors() {
    let cfg = WelchConfig {
        segment_samples: 8,
        hop_samples: 4,
        window: WelchWindow::PeriodicHann,
        detrend: WelchDetrend::None,
    };
    let original: Vec<[f64; 2]> = (0..64)
        .map(|j| [(0.31 * j as f64).cos(), (0.27 * j as f64).sin()])
        .collect();
    let accepted: Vec<_> = original[..17]
        .iter()
        .copied()
        .chain([[0.0; 2]; 47])
        .collect();
    let expected =
        WelchCrossSpectrum::estimate_interleaved(&accepted, 0.001, cfg, budget(), || false)
            .unwrap();
    let prepare = || {
        let mut stream = WelchStream::new(2, 64, 0.001, cfg, budget(), || false).unwrap();
        stream.push_interleaved(&original[..17], || false).unwrap();
        stream
    };
    let mut total_polls = 0;
    prepare()
        .push_interleaved(&original[17..41], || {
            total_polls += 1;
            false
        })
        .unwrap();
    assert!(total_polls > 35);
    // Before input admission, after multiple FFTs, and after replacement-tail
    // construction. The retry is dark where the cancelled chunk was nonzero.
    for stop in [1, 35, total_polls] {
        let mut stream = prepare();
        let mut polls = 0;
        assert!(matches!(
            stream.push_interleaved(&original[17..41], || {
                polls += 1;
                polls == stop
            }),
            Err(WelchError::Cancelled)
        ));
        assert_eq!(stream.accepted_sample_count(), 17);
        assert_eq!(stream.completed_segment_count(), 3);
        stream.push_interleaved(&accepted[17..], || false).unwrap();
        same_stream_spectrum(&stream.finish().unwrap(), &expected);
    }
    // The no-complete-window append path also rolls back a partially copied
    // chunk without overwriting the previous accepted prefix.
    let cfg = config(WelchWindow::Rectangular, WelchDetrend::SegmentMean);
    let accepted: Vec<_> = original[..5]
        .iter()
        .copied()
        .chain([[0.0; 2]; 27])
        .collect();
    let expected =
        WelchCrossSpectrum::estimate_interleaved(&accepted, 0.001, cfg, budget(), || false)
            .unwrap();
    let mut stream = WelchStream::new(2, 32, 0.001, cfg, budget(), || false).unwrap();
    stream.push_interleaved(&original[..5], || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        stream.push_interleaved(&original[5..22], || {
            polls += 1;
            polls == 11
        }),
        Err(WelchError::Cancelled)
    ));
    assert_eq!(stream.accepted_sample_count(), 5);
    stream.push_interleaved(&accepted[5..], || false).unwrap();
    same_stream_spectrum(&stream.finish().unwrap(), &expected);
}

#[test]
fn g0_stream_rejects_shape_overrun_and_nonfinite_discarded_tail_atomically() {
    let cfg = WelchConfig {
        segment_samples: 8,
        hop_samples: 4,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::None,
    };
    let rows = [[1.0, -1.0]; 19];
    let mut stream = WelchStream::new(2, rows.len(), 0.001, cfg, budget(), || false).unwrap();
    stream.push_interleaved(&rows[..16], || false).unwrap();
    assert!(matches!(
        stream.push_interleaved(&[[f64::NAN; 2]; 4], || false),
        Err(WelchError::Invalid(_))
    ));
    assert!(matches!(
        stream.push(&[&[1.0], &[]], || false),
        Err(WelchError::Invalid(_))
    ));
    assert!(matches!(
        stream.push_interleaved(&[[1.0; 3]], || false),
        Err(WelchError::Invalid(_))
    ));
    let bad_tail = [[1.0, -1.0], [f64::INFINITY, -1.0], [1.0, -1.0]];
    assert!(matches!(
        stream.push_interleaved(&bad_tail, || false),
        Err(WelchError::NonFinite)
    ));
    assert_eq!(stream.accepted_sample_count(), 16);
    assert_eq!(stream.completed_segment_count(), 3);
    stream.push_interleaved(&rows[16..], || false).unwrap();
    let actual = stream.finish().unwrap();
    assert_eq!(actual.used_sample_count(), 16);
    assert_eq!(actual.discarded_tail_samples(), 3);
    same_stream_spectrum(
        &actual,
        &WelchCrossSpectrum::estimate_interleaved(&rows, 0.001, cfg, budget(), || false).unwrap(),
    );
}

#[test]
fn g1_stream_keeps_batch_dynamic_range_and_transactional_numerical_refusal() {
    let cfg = WelchConfig {
        segment_samples: 4,
        hop_samples: 4,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::None,
    };
    for (large, dt) in [(1e308, 1e-300), (1e300, 1e-100)] {
        let values = [large, 1.0, large, -1.0];
        let expected =
            WelchCrossSpectrum::estimate(&[&values], dt, cfg, budget(), || false).unwrap();
        let mut stream = WelchStream::new(1, 4, dt, cfg, budget(), || false).unwrap();
        for chunk in values.chunks(1) {
            stream.push(&[chunk], || false).unwrap();
        }
        same_stream_spectrum(&stream.finish().unwrap(), &expected);
    }
    let mut stream = WelchStream::new(1, 4, 1e-100, cfg, budget(), || false).unwrap();
    stream.push(&[&[1e300]], || false).unwrap();
    assert!(matches!(
        stream.push(&[&[1e-30, 1e300, -1e-30]], || false),
        Err(WelchError::Unrepresentable)
    ));
    assert_eq!(stream.accepted_sample_count(), 1);
    assert_eq!(stream.completed_segment_count(), 0);
    stream.push(&[&[0.0, 1e300, 0.0]], || false).unwrap();
    let expected =
        WelchCrossSpectrum::estimate(&[&[1e300, 0.0, 1e300, 0.0]], 1e-100, cfg, budget(), || {
            false
        })
        .unwrap();
    same_stream_spectrum(&stream.finish().unwrap(), &expected);

    let cfg = WelchConfig {
        detrend: WelchDetrend::SegmentMean,
        ..cfg
    };
    let d = f64::EPSILON;
    let rows = [
        [1.5 + d, 1.75],
        [1.5, 1.75 + d],
        [1.5 - d, 1.75],
        [1.5, 1.75 - d],
    ];
    let expected =
        WelchCrossSpectrum::estimate_interleaved(&rows, 1.0, cfg, budget(), || false).unwrap();
    let mut stream = WelchStream::new(2, 4, 1.0, cfg, budget(), || false).unwrap();
    for row in &rows {
        stream
            .push_interleaved(core::slice::from_ref(row), || false)
            .unwrap();
    }
    same_stream_spectrum(&stream.finish().unwrap(), &expected);
}

#[test]
fn g0_stream_admits_full_record_resources_and_never_pads_missing_samples() {
    let cfg = config(WelchWindow::Rectangular, WelchDetrend::None);
    for limits in [
        WelchBudget {
            maximum_channels: 1,
            ..budget()
        },
        WelchBudget {
            maximum_samples_per_channel: 63,
            ..budget()
        },
        WelchBudget {
            maximum_transforms: 5,
            ..budget()
        },
        WelchBudget {
            maximum_workspace_bytes: 128,
            ..budget()
        },
    ] {
        assert!(matches!(
            WelchStream::new(2, 64, 0.001, cfg, limits, || false),
            Err(WelchError::BudgetExceeded(_))
        ));
    }
    assert!(matches!(
        WelchStream::new(2, 64, 0.001, cfg, budget(), || true),
        Err(WelchError::Cancelled)
    ));
    let mut stream = WelchStream::new(2, 64, 0.001, cfg, budget(), || false).unwrap();
    stream.push_interleaved(&[[1.0; 2]; 63], || false).unwrap();
    assert_eq!(stream.remaining_sample_count(), 1);
    assert!(matches!(stream.finish(), Err(WelchError::Invalid(_))));
}
