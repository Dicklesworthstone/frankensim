//! Simultaneous physical components retain one frequency and their vector metric.

use fs_fft::C64;
use fs_fft::harmonic::{
    HarmonicFitBudget, HarmonicFitError, HarmonicFrequencyBudget, HarmonicFrequencyError,
    refine_harmonic_frequency, refine_harmonic_frequency_channels,
};

fn budget() -> HarmonicFrequencyBudget {
    HarmonicFrequencyBudget {
        frequency_bounds_hz: [215.0, 220.0],
        absolute_tolerance_hz: 1e-8,
        max_evaluations: 24,
        fit: HarmonicFitBudget {
            max_workspace_bytes: 8 * 1024 * 1024,
            min_reciprocal_condition: 1e-10,
        },
    }
}

fn record(frequency: f64, dc: f64, phasors: &[C64]) -> Vec<f64> {
    (0..2039)
        .map(|j| {
            dc + phasors
                .iter()
                .enumerate()
                .map(|(k, p)| {
                    let (s, c) = (core::f64::consts::TAU * frequency * (k + 1) as f64 * j as f64
                        / 8000.0)
                        .sin_cos();
                    p.re * c + p.im * s
                })
                .sum::<f64>()
        })
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn scalar_front_reuses_exactly_the_same_single_channel_numerics() {
    let data = record(217.31, 7.0, &[C64::new(0.7, -0.2), C64::new(0.001, 0.002)]);
    let scalar = refine_harmonic_frequency(&data, 1.0 / 8000.0, 217.8, &[1, 2], budget()).unwrap();
    let vector =
        refine_harmonic_frequency_channels(&[&data], 1.0 / 8000.0, 217.8, &[1, 2], budget())
            .unwrap();
    assert_eq!(scalar.fundamental_hz, vector.fundamental_hz);
    assert_eq!(scalar.fit, vector.fits[0]);
    assert_eq!(scalar.evaluations, vector.evaluations);
    assert_eq!(scalar.last_correction_hz, vector.last_correction_hz);
}

#[test]
fn physical_vector_rotation_retains_pitch_and_rotates_every_phasor() {
    let frequency = 217.31;
    let channels = [
        record(frequency, 0.2, &[C64::new(0.0, 0.0); 3]),
        record(
            frequency,
            -1.3,
            &[
                C64::new(0.0, 0.0),
                C64::new(0.7, -0.2),
                C64::new(1e-5, 2e-5),
            ],
        ),
        record(
            frequency,
            2.1,
            &[C64::new(0.0, 0.0), C64::new(0.2, 0.4), C64::new(-0.3, 0.1)],
        ),
    ];
    let refs: Vec<_> = channels.iter().map(Vec::as_slice).collect();
    let original =
        refine_harmonic_frequency_channels(&refs, 1.0 / 8000.0, 217.8, &[1, 2, 3], budget())
            .unwrap();
    assert_eq!(original.fits[0].dc, 0.2);
    assert_eq!(original.fits[0].phasors, [C64::new(0.0, 0.0); 3]);
    let (s, c) = 0.37_f64.sin_cos();
    let (u, v) = 0.61_f64.sin_cos();
    let rotation = [[c, -s * v, s * u], [s, c * v, -c * u], [0.0, u, v]];
    let rotated: Vec<Vec<f64>> = rotation
        .iter()
        .map(|row| {
            (0..channels[0].len())
                .map(|j| (0..3).map(|i| row[i] * channels[i][j]).sum())
                .collect()
        })
        .collect();
    let refs: Vec<_> = rotated.iter().map(Vec::as_slice).collect();
    let result =
        refine_harmonic_frequency_channels(&refs, 1.0 / 8000.0, 217.8, &[1, 2, 3], budget())
            .unwrap();
    close(result.fundamental_hz, frequency, 2e-8);
    close(result.fundamental_hz, original.fundamental_hz, 2e-10);
    for (i, row) in rotation.iter().enumerate() {
        close(
            result.fits[i].dc,
            (0..3).map(|k| row[k] * original.fits[k].dc).sum(),
            2e-10,
        );
        for h in 0..3 {
            close(
                result.fits[i].phasors[h].re,
                (0..3)
                    .map(|k| row[k] * original.fits[k].phasors[h].re)
                    .sum(),
                2e-9,
            );
            close(
                result.fits[i].phasors[h].im,
                (0..3)
                    .map(|k| row[k] * original.fits[k].phasors[h].im)
                    .sum(),
                2e-9,
            );
        }
        assert!(result.fits[i].residual_rms < 1e-9);
    }
}

#[test]
fn weak_conflicting_channel_does_not_receive_invented_equal_amplitude_weight() {
    // The one-frequency model is deliberately imperfect: its physical least
    // squares compromise must be dominated by the strong source, not by an
    // average of separate pitches or unit normalization of both channels.
    let strong = record(217.31, 3.0, &[C64::new(1.0, 0.0)]);
    let weak = record(217.91, -2.0, &[C64::new(0.001, 0.0)]);
    let mut b = budget();
    b.absolute_tolerance_hz = 1e-7;
    let result =
        refine_harmonic_frequency_channels(&[&strong, &weak], 1.0 / 8000.0, 217.5, &[1], b)
            .unwrap();
    close(result.fundamental_hz, 217.31, 2e-6);
    assert!(result.fits[1].residual_rms > 1e-4);
    assert!(result.fits[0].residual_rms < 1e-6);
}

#[test]
fn vector_fit_preserves_common_unit_scaling_and_constant_preloads() {
    let signal = record(217.31, 3.0, &[C64::new(0.7, -0.2), C64::new(0.001, 0.002)]);
    for scale in [1e-150, 1.0, 1e150] {
        let data: Vec<_> = signal.iter().map(|x| x * scale).collect();
        let constant = vec![9.0 * scale; data.len()];
        let result = refine_harmonic_frequency_channels(
            &[&constant, &data],
            1.0 / 8000.0,
            217.8,
            &[1, 2],
            budget(),
        )
        .unwrap();
        close(result.fundamental_hz, 217.31, 1e-8);
        close(result.fits[1].dc / scale, 3.0, 1e-10);
        close(result.fits[1].phasors[1].re / scale, 0.001, 1e-10);
        close(result.fits[1].phasors[1].im / scale, 0.002, 1e-10);
        assert_eq!(result.fits[0].dc, 9.0 * scale);
        assert_eq!(result.fits[0].residual_rms, 0.0);
    }
}

#[test]
fn all_channels_are_admitted_and_retained_storage_is_budgeted() {
    let data = record(217.31, 3.0, &[C64::new(1.0, 0.0)]);
    let mut bad = data.clone();
    bad[300] = f64::NAN;
    assert!(matches!(
        refine_harmonic_frequency_channels(&[&data, &bad], 1.0 / 8000.0, 217.8, &[1], budget()),
        Err(HarmonicFrequencyError::Fit(
            HarmonicFitError::NonFinite { .. }
        ))
    ));
    assert!(matches!(
        refine_harmonic_frequency_channels(
            &[&data, &data[..300]],
            1.0 / 8000.0,
            217.8,
            &[1],
            budget()
        ),
        Err(HarmonicFrequencyError::InvalidParameter { .. })
    ));
    assert!(matches!(
        refine_harmonic_frequency_channels(&[], 1.0 / 8000.0, 217.8, &[1], budget()),
        Err(HarmonicFrequencyError::InvalidParameter { .. })
    ));
    assert_eq!(
        refine_harmonic_frequency_channels(
            &[&[3.0; 31], &[7.0; 31]],
            1.0 / 8000.0,
            217.8,
            &[1],
            budget()
        ),
        Err(HarmonicFrequencyError::Unidentifiable)
    );
    let mut b = budget();
    b.fit.max_workspace_bytes = 0;
    let Err(HarmonicFrequencyError::Fit(HarmonicFitError::BudgetExceeded {
        required_bytes, ..
    })) = refine_harmonic_frequency_channels(&[&data], 1.0 / 8000.0, 217.8, &[1], b)
    else {
        panic!("storage refusal required")
    };
    b.fit.max_workspace_bytes = required_bytes;
    refine_harmonic_frequency_channels(&[&data], 1.0 / 8000.0, 217.8, &[1], b).unwrap();
    assert!(matches!(
        refine_harmonic_frequency_channels(&[&data, &data], 1.0 / 8000.0, 217.8, &[1], b),
        Err(HarmonicFrequencyError::Fit(
            HarmonicFitError::BudgetExceeded { .. }
        ))
    ));
}
