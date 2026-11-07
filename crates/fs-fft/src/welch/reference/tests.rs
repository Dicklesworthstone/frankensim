use super::*;

fn fit_budget() -> WelchReferenceFitBudget {
    WelchReferenceFitBudget {
        maximum_reference_channels: 16,
        maximum_target_channels: 16,
        maximum_work_units: 1 << 26,
        maximum_workspace_bytes: 1 << 24,
    }
}
fn apply_budget() -> WelchChannelTransferBudget {
    WelchChannelTransferBudget {
        maximum_output_channels: 16,
        maximum_complex_products: 1 << 24,
        maximum_workspace_bytes: 1 << 24,
    }
}
fn config(n: usize) -> WelchConfig {
    WelchConfig {
        segment_samples: n,
        hop_samples: n,
        window: super::super::WelchWindow::Rectangular,
        detrend: super::super::WelchDetrend::None,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}
fn complex_close(a: C64, b: C64, tolerance: f64) {
    close(a.re, b.re, tolerance);
    close(a.im, b.im, tolerance);
}

// These are realizable real-record DFT factors with one explicitly selected bin.
// Endpoint construction enforces the same real-signal symmetry as estimate().
fn at_bin(n: usize, bin: usize, rows: &[Vec<C64>]) -> WelchCrossSpectrum {
    let channels = rows[0].len();
    let segments = rows.len();
    let mut factors = vec![C64::default(); (n / 2 + 1) * channels * segments];
    for (s, row) in rows.iter().enumerate() {
        assert_eq!(row.len(), channels);
        for (c, &value) in row.iter().enumerate() {
            if bin == 0 || (n % 2 == 0 && bin == n / 2) {
                assert_eq!(value.im, 0.0);
            }
            factors[(bin * segments + s) * channels + c] = value;
        }
    }
    WelchCrossSpectrum {
        config: config(n),
        timestep_s: 0.125,
        bin_width_hz: 8.0 / n as f64,
        channels,
        segments,
        input_samples: n * segments,
        used_samples: n * segments,
        factors,
    }
}

fn observed(n: usize, bin: usize, amplitudes: &[Vec<C64>]) -> WelchCrossSpectrum {
    let channels = amplitudes[0].len();
    let mut samples = vec![vec![0.0; n * amplitudes.len()]; channels];
    for (r, row) in amplitudes.iter().enumerate() {
        for j in 0..n {
            let (s, c) = (core::f64::consts::TAU * bin as f64 * j as f64 / n as f64).sin_cos();
            for k in 0..channels {
                samples[k][r * n + j] = row[k].re * c + row[k].im * s;
            }
        }
    }
    let views: Vec<&[f64]> = samples.iter().map(Vec::as_slice).collect();
    WelchCrossSpectrum::estimate(
        &views,
        0.125,
        config(n),
        super::super::WelchBudget {
            maximum_channels: 16,
            maximum_samples_per_channel: 4096,
            maximum_transforms: 1024,
            maximum_workspace_bytes: 1 << 24,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn g1_real_records_learn_complex_gains_and_cancel_a_separate_record() {
    let (g, h) = (
        [C64::new(0.4, 0.7), C64::new(-0.3, 0.2)],
        [C64::new(-0.2, 0.4), C64::new(0.8, -0.1)],
    );
    let training: Vec<Vec<C64>> = (0..8)
        .map(|r| {
            let a = C64::new(if r % 4 < 2 { 1.0 } else { -1.0 }, 0.0);
            let b = C64::new(0.0, if r % 2 == 0 { 1.0 } else { -1.0 });
            let noise = C64::new(if r < 4 { 0.3 } else { -0.3 }, 0.0);
            vec![
                a,
                b,
                g[0].mul(a).add(g[1].mul(b)).add(noise),
                h[0].mul(a)
                    .add(h[1].mul(b))
                    .add(C64::new(0.0, 2.0).mul(noise)),
            ]
        })
        .collect();
    let record = observed(16, 2, &training);
    let model = record
        .fit_reference_transfer(&[2, 3], &[0, 1], 1e-12, fit_budget(), || false)
        .unwrap();
    assert_eq!(model.rank(2).unwrap(), 2);
    assert_eq!(model.training_segment_count(), 8);
    assert_eq!(model.target_channels(), &[2, 3]);
    assert_eq!(model.reference_channels(), &[0, 1]);
    for r in 0..2 {
        complex_close(model.gain(2, 0, r).unwrap(), g[r], 2e-14);
        complex_close(model.gain(2, 1, r).unwrap(), h[r], 2e-14);
    }
    let mut holdout = Vec::new();
    let mut noises = Vec::new();
    for r in 0..5 {
        let a = C64::new(0.2 + r as f64 * 0.13, -0.1 + r as f64 * 0.03);
        let b = C64::new(-0.7 + r as f64 * 0.09, 0.5 - r as f64 * 0.11);
        // These held-out residuals need not be orthogonal to the references.
        let e = C64::new(0.03 * (r + 1) as f64, -0.02 * (r + 2) as f64);
        let f = C64::new(0.0, 2.0).mul(e);
        holdout.push(vec![
            a,
            b,
            g[0].mul(a).add(g[1].mul(b)).add(e),
            h[0].mul(a).add(h[1].mul(b)).add(f),
        ]);
        noises.push(vec![e, f]);
    }
    let held = observed(16, 2, &holdout);
    let truth = observed(16, 2, &noises);
    let residual = model
        .residual_spectrum(&held, apply_budget(), || false)
        .unwrap();
    assert_eq!(residual.segment_count(), 5);
    assert_eq!(residual.config(), held.config());
    assert_eq!(residual.timestep_s(), held.timestep_s());
    assert_eq!(residual.bin_width_hz(), held.bin_width_hz());
    assert_eq!(residual.used_sample_count(), held.used_sample_count());
    for a in 0..2 {
        for b in 0..2 {
            complex_close(
                residual.cross_spectrum(2, a, b).unwrap(),
                truth.cross_spectrum(2, a, b).unwrap(),
                3e-14,
            );
        }
    }
    assert!(residual.cross_spectrum(2, 0, 1).unwrap().im < 0.0);
    let freq = residual.frequency_hz(2).unwrap();
    assert_eq!(residual.peak_bin_hz(2, freq).unwrap().frequency_hz(), freq);
}

#[test]
fn g3_dependent_zero_and_rescaled_references_preserve_the_residual() {
    let g = [C64::new(0.2, 0.7), C64::new(-0.4, 0.1)];
    let rows: Vec<Vec<C64>> = (0..8)
        .map(|r| {
            let a = C64::new(if r % 4 < 2 { 1.0 } else { -1.0 }, 0.0);
            let b = C64::new(0.0, if r % 2 == 0 { 1.0 } else { -1.0 });
            let noise = C64::new(if r < 4 { 0.125 } else { -0.125 }, 0.0);
            vec![
                a,
                b,
                a.scale(2.0),
                C64::default(),
                g[0].mul(a).add(g[1].mul(b)).add(noise),
            ]
        })
        .collect();
    let record = at_bin(8, 1, &rows);
    let model = record
        .fit_reference_transfer(&[4], &[0, 1, 2, 3], 1e-12, fit_budget(), || false)
        .unwrap();
    let reversed = record
        .fit_reference_transfer(&[4], &[2, 1, 0, 3], 1e-12, fit_budget(), || false)
        .unwrap();
    assert_eq!(model.rank(1).unwrap(), 2);
    assert_eq!(reversed.rank(1).unwrap(), 2);
    assert_eq!(model.rank(0).unwrap(), 0);
    let expected = model
        .residual_spectrum(&record, apply_budget(), || false)
        .unwrap();
    let reverse = reversed
        .residual_spectrum(&record, apply_budget(), || false)
        .unwrap();
    complex_close(
        expected.cross_spectrum(1, 0, 0).unwrap(),
        C64::new(0.125, 0.0),
        2e-14,
    );
    complex_close(
        reverse.cross_spectrum(1, 0, 0).unwrap(),
        expected.cross_spectrum(1, 0, 0).unwrap(),
        2e-14,
    );
    let scaled: Vec<Vec<C64>> = rows
        .iter()
        .map(|row| {
            let mut values = row.clone();
            values[0] = values[0].scale(1e-120);
            values[1] = values[1].scale(1e120);
            values
        })
        .collect();
    let record = at_bin(8, 1, &scaled);
    let units = record
        .fit_reference_transfer(&[4], &[0, 1], 1e-12, fit_budget(), || false)
        .unwrap();
    let actual = units
        .residual_spectrum(&record, apply_budget(), || false)
        .unwrap();
    complex_close(
        actual.cross_spectrum(1, 0, 0).unwrap(),
        expected.cross_spectrum(1, 0, 0).unwrap(),
        3e-14,
    );
    complex_close(units.gain(1, 0, 0).unwrap().scale(1e-120), g[0], 3e-14);
    complex_close(units.gain(1, 0, 1).unwrap().scale(1e120), g[1], 3e-14);
    let zero = record
        .fit_reference_transfer(&[4], &[3], 1e-12, fit_budget(), || false)
        .unwrap();
    assert_eq!(zero.rank(1).unwrap(), 0);
    let unchanged = zero
        .residual_spectrum(&record, apply_budget(), || false)
        .unwrap();
    for (a, b) in unchanged
        .factors(1)
        .unwrap()
        .zip(record.factors(1).unwrap())
    {
        assert_eq!(a[0], b[4]);
    }
}

#[test]
fn g1_factor_cancellation_preserves_power_lost_by_subtracting_spectra() {
    let epsilon = 2.0f64.powi(-50);
    let rows: Vec<Vec<C64>> = (0..4)
        .map(|r| {
            let reference = if r < 2 { 1.0 } else { -1.0 };
            let quiet = if r % 2 == 0 { epsilon } else { -epsilon };
            vec![C64::new(reference, 0.0), C64::new(reference + quiet, 0.0)]
        })
        .collect();
    let record = at_bin(8, 1, &rows);
    assert_eq!(
        record.cross_spectrum(1, 1, 1).unwrap().re - record.cross_spectrum(1, 0, 0).unwrap().re,
        0.0
    );
    let model = record
        .fit_reference_transfer(&[1], &[0], 1e-14, fit_budget(), || false)
        .unwrap();
    let residual = model
        .residual_spectrum(&record, apply_budget(), || false)
        .unwrap();
    let power = residual.cross_spectrum(1, 0, 0).unwrap().re;
    assert!(power > 0.0);
    assert!((power / (4.0 * epsilon * epsilon) - 1.0).abs() < 3e-14);
}

#[test]
fn g1_nearly_parallel_references_use_declared_rank_and_stable_qr() {
    let delta = 2.0f64.powi(-20);
    let epsilon = 2.0f64.powi(-30);
    let rows: Vec<Vec<C64>> = (0..4)
        .map(|r| {
            let u = if r % 2 == 0 { 1.0 } else { -1.0 };
            let w = if r < 2 { 1.0 } else { -1.0 };
            vec![
                C64::new(1.0, 0.0),
                C64::new(1.0 + delta * u, 0.0),
                C64::new(u + epsilon * w, 0.0),
            ]
        })
        .collect();
    let record = at_bin(8, 1, &rows);
    let model = record
        .fit_reference_transfer(&[2], &[0, 1], 1e-8, fit_budget(), || false)
        .unwrap();
    assert_eq!(model.rank(1).unwrap(), 2);
    close(model.gain(1, 0, 0).unwrap().re, -1.0 / delta, 2e-14);
    close(model.gain(1, 0, 1).unwrap().re, 1.0 / delta, 2e-14);
    let residual = model
        .residual_spectrum(&record, apply_budget(), || false)
        .unwrap();
    let expected = 4.0 * epsilon * epsilon;
    assert!((residual.cross_spectrum(1, 0, 0).unwrap().re / expected - 1.0).abs() < 1e-10);
    let truncated = record
        .fit_reference_transfer(&[2], &[0, 1], 1e-5, fit_budget(), || false)
        .unwrap();
    assert_eq!(truncated.rank(1).unwrap(), 1);
    assert!(
        truncated
            .residual_spectrum(&record, apply_budget(), || false)
            .unwrap()
            .cross_spectrum(1, 0, 0)
            .unwrap()
            .re
            > 3.9
    );
}

#[test]
fn g0_real_endpoints_and_odd_top_bin_keep_their_original_phase_laws() {
    for (n, bin, gain) in [
        (8, 0, C64::new(0.3, 0.0)),
        (8, 4, C64::new(-0.7, 0.0)),
        (9, 4, C64::new(0.2, 0.8)),
    ] {
        let rows: Vec<Vec<C64>> = [1.0, -1.0, 2.0, -2.0]
            .into_iter()
            .map(|a| {
                let x = C64::new(a, 0.0);
                vec![x, gain.mul(x)]
            })
            .collect();
        let record = at_bin(n, bin, &rows);
        let model = record
            .fit_reference_transfer(&[1], &[0], 1e-12, fit_budget(), || false)
            .unwrap();
        complex_close(model.gain(bin, 0, 0).unwrap(), gain, 2e-14);
        let residual = model
            .residual_spectrum(&record, apply_budget(), || false)
            .unwrap();
        assert!(residual.cross_spectrum(bin, 0, 0).unwrap().re < 1e-28);
        if n % 2 == 0 {
            assert_eq!(model.gain(bin, 0, 0).unwrap().im, 0.0);
            for row in residual.factors(bin).unwrap() {
                assert_eq!(row[0].im, 0.0);
            }
        }
    }
}

#[test]
fn g3_single_segment_fit_is_not_claimed_as_holdout_coherence() {
    let training = at_bin(8, 1, &[vec![C64::new(1.0, 0.0), C64::new(0.0, 1.0)]]);
    let model = training
        .fit_reference_transfer(&[1], &[0], 1e-12, fit_budget(), || false)
        .unwrap();
    assert_eq!(model.rank(1).unwrap(), 1);
    assert_eq!(
        model
            .residual_spectrum(&training, apply_budget(), || false)
            .unwrap()
            .cross_spectrum(1, 0, 0)
            .unwrap()
            .re,
        0.0
    );
    let held = at_bin(
        8,
        1,
        &[
            vec![C64::new(1.0, 0.0), C64::new(1.0, 0.0)],
            vec![C64::new(1.0, 0.0), C64::new(-1.0, 0.0)],
        ],
    );
    let residual = model
        .residual_spectrum(&held, apply_budget(), || false)
        .unwrap();
    close(residual.cross_spectrum(1, 0, 0).unwrap().re, 4.0, 1e-14);
}

#[test]
fn g4_bounds_invalid_shapes_cancellation_and_replay_are_transactional() {
    let rows: Vec<Vec<C64>> = (0..200)
        .map(|r| {
            let x = C64::new((r as f64 * 0.37).sin(), (r as f64 * 0.21).cos());
            let y = C64::new((r as f64 * 0.17).cos(), (r as f64 * 0.13).sin());
            vec![x, y, x.add(y), x.sub(y)]
        })
        .collect();
    let record = at_bin(8, 1, &rows);
    let original = record.factors.clone();
    let needs = record.reference_fit_requirements(2, 2).unwrap();
    let exact = WelchReferenceFitBudget {
        maximum_reference_channels: 2,
        maximum_target_channels: 2,
        maximum_work_units: needs.work_units,
        maximum_workspace_bytes: needs.workspace_bytes,
    };
    let model = record
        .fit_reference_transfer(&[2, 3], &[0, 1], 1e-12, exact, || false)
        .unwrap();
    for budget in [
        WelchReferenceFitBudget {
            maximum_work_units: needs.work_units - 1,
            ..exact
        },
        WelchReferenceFitBudget {
            maximum_workspace_bytes: needs.workspace_bytes - 1,
            ..exact
        },
        WelchReferenceFitBudget {
            maximum_reference_channels: 1,
            ..exact
        },
    ] {
        assert!(
            record
                .fit_reference_transfer(&[2, 3], &[0, 1], 1e-12, budget, || false)
                .is_err()
        );
    }
    for (targets, refs) in [
        (&[][..], &[0][..]),
        (&[2][..], &[][..]),
        (&[2, 2][..], &[0][..]),
        (&[2][..], &[0, 0][..]),
        (&[2][..], &[2][..]),
        (&[8][..], &[0][..]),
    ] {
        assert!(
            record
                .fit_reference_transfer(targets, refs, 1e-12, fit_budget(), || false)
                .is_err()
        );
    }
    for tolerance in [0.0, 1.0, -0.1, f64::NAN] {
        assert!(
            record
                .fit_reference_transfer(&[2], &[0], tolerance, fit_budget(), || false)
                .is_err()
        );
    }
    assert!(record.reference_fit_requirements(usize::MAX, 1).is_err());
    let mut polls = 0;
    assert!(matches!(
        record.fit_reference_transfer(&[2, 3], &[0, 1], 1e-12, exact, || {
            polls += 1;
            polls == 3
        }),
        Err(WelchError::Cancelled)
    ));
    assert_eq!(record.factors, original);
    let replay = record
        .fit_reference_transfer(&[2, 3], &[0, 1], 1e-12, exact, || false)
        .unwrap();
    assert_eq!(model.gains, replay.gains);
    assert_eq!(model.ranks, replay.ranks);
    let needs = model.residual_requirements(&record).unwrap();
    assert_eq!(needs.coefficient_values, 0);
    let exact = WelchChannelTransferBudget {
        maximum_output_channels: 2,
        maximum_complex_products: needs.complex_products,
        maximum_workspace_bytes: needs.workspace_bytes,
    };
    let residual = model.residual_spectrum(&record, exact, || false).unwrap();
    for budget in [
        WelchChannelTransferBudget {
            maximum_complex_products: needs.complex_products - 1,
            ..exact
        },
        WelchChannelTransferBudget {
            maximum_workspace_bytes: needs.workspace_bytes - 1,
            ..exact
        },
    ] {
        assert!(model.residual_spectrum(&record, budget, || false).is_err());
    }
    let mut polls = 0;
    assert!(matches!(
        model.residual_spectrum(&record, exact, || {
            polls += 1;
            polls == 3
        }),
        Err(WelchError::Cancelled)
    ));
    assert_eq!(model.gains, replay.gains);
    assert_eq!(
        model
            .residual_spectrum(&record, exact, || false)
            .unwrap()
            .factors,
        residual.factors
    );
    let mut wrong = record.clone();
    wrong.timestep_s *= 2.0;
    assert!(model.residual_spectrum(&wrong, exact, || false).is_err());
    wrong = record.clone();
    wrong.config.window = super::super::WelchWindow::PeriodicHann;
    assert!(model.residual_spectrum(&wrong, exact, || false).is_err());
    // Every segment can be representable while across-segment normalization
    // would erase a nonzero reference component. Refuse that range loss.
    let wide = at_bin(
        8,
        1,
        &[
            vec![C64::new(1e150, 0.0), C64::new(1.0, 0.0)],
            vec![C64::new(1e-200, 0.0), C64::new(-1.0, 0.0)],
        ],
    );
    assert!(matches!(
        wide.fit_reference_transfer(&[1], &[0], 1e-12, fit_budget(), || false),
        Err(WelchError::Unrepresentable)
    ));
    let absent = at_bin(
        8,
        1,
        &[
            vec![C64::default(), C64::new(1e150, 0.0)],
            vec![C64::default(), C64::new(1e-200, 0.0)],
        ],
    );
    let identity = absent
        .fit_reference_transfer(&[1], &[0], 1e-12, fit_budget(), || false)
        .unwrap();
    let unchanged = identity
        .residual_spectrum(&absent, apply_budget(), || false)
        .unwrap();
    for (residual, input) in unchanged
        .factors(1)
        .unwrap()
        .zip(absent.factors(1).unwrap())
    {
        assert_eq!(residual[0], input[1]);
    }
}
