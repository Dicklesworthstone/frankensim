//! G0/G1/G3: independent force-spectrum and known-frequency projection checks.

use super::*;
use fs_fft::harmonic::{HarmonicFitBudget, HarmonicFitError};
use std::f64::consts::TAU;

fn budget() -> HarmonicFitBudget {
    HarmonicFitBudget {
        max_workspace_bytes: 4 * 1024 * 1024,
        min_reciprocal_condition: 1e-10,
    }
}

fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "actual={actual:e}, expected={expected:e}, tolerance={absolute:e}"
    );
}

#[test]
fn g1_prepared_force_spectrum_matches_independent_hann_dft() {
    let n = 128;
    let dt = 0.0001;
    let forces: Vec<_> = (0..n)
        .map(|i| {
            let theta = TAU * i as f64 / n as f64;
            [
                2.0 + 1.3 * (4.0 * theta).cos() - 0.7 * (9.0 * theta).sin(),
                -3.0 + 0.4 * (4.0 * theta).sin() + 0.03 * (31.0 * theta).cos(),
            ]
        })
        .collect();
    let spectrum = ForceSpectrum::new(&forces, dt, n).unwrap();
    for bin in [1, 4, 9, 31, 63] {
        let (actual, frequency) = spectrum.line(bin).unwrap();
        close(frequency, bin as f64 / (n as f64 * dt), 0.0);
        for axis in 0..2 {
            let mean = forces.iter().map(|f| f[axis]).sum::<f64>() / n as f64;
            let mut expected = C64::ZERO;
            for (i, force) in forces.iter().enumerate() {
                let theta = TAU * bin as f64 * i as f64 / n as f64;
                let hann = 0.5 - 0.5 * (TAU * i as f64 / (n - 1) as f64).cos();
                expected = expected
                    + C64::new(theta.cos(), theta.sin())
                        .scale((force[axis] - mean) * hann * 4.0 / (n - 1) as f64);
            }
            assert!((actual[axis] - expected).abs() < 3e-14);
        }
        assert_eq!(spectrum.line(bin).unwrap().0, actual);
    }
}

#[test]
fn g1_physical_force_spectrum_retains_prime_and_odd_window_frequencies() {
    for n in [65_usize, 96, 101] {
        let dt = 1.0 / 8000.0;
        let edge = (n - 1) / 2;
        let record: Vec<_> = (0..n)
            .map(|i| {
                let phase = TAU * i as f64 / n as f64;
                [
                    0.7 + (7.0 * phase).cos(),
                    -0.4 + (edge as f64 * phase).sin(),
                ]
            })
            .collect();
        let spectrum = ForceSpectrum::new(&record, dt, n).unwrap();
        assert_eq!(spectrum.sample_count, n);
        for bin in [1, 7, edge] {
            let (actual, frequency) = spectrum.line(bin).unwrap();
            assert_eq!(frequency, bin as f64 / (n as f64 * dt));
            for axis in 0..2 {
                let mean = record.iter().map(|f| f[axis]).sum::<f64>() / n as f64;
                let mut expected = C64::ZERO;
                for (i, row) in record.iter().enumerate() {
                    let angle = TAU * bin as f64 * i as f64 / n as f64;
                    let window = 0.5 - 0.5 * (TAU * i as f64 / (n - 1) as f64).cos();
                    expected = expected
                        + C64::new(angle.cos(), angle.sin())
                            .scale((row[axis] - mean) * window * 4.0 / (n - 1) as f64);
                }
                assert!((actual[axis] - expected).abs() < 2e-12);
            }
        }
        assert!(spectrum.line(n.div_ceil(2)).is_err());
    }
}

fn harmonic_record(n: usize, dt: f64, fundamental: f64) -> Vec<[f64; 2]> {
    (0..n)
        .map(|i| {
            let theta = TAU * fundamental * i as f64 * dt;
            [
                0.4 + 1.3 * theta.cos() - 0.2 * theta.sin()
                    + 0.00001 * (2.0 * theta).cos()
                    + 0.00002 * (2.0 * theta).sin()
                    - 0.15 * (3.0 * theta).cos(),
                -0.8 + 0.3 * theta.sin() - 0.07 * (2.0 * theta).cos() + 0.25 * (3.0 * theta).sin(),
            ]
        })
        .collect()
}

#[test]
fn g1_joint_force_fit_resolves_weak_off_bin_even_harmonics_and_radiation() {
    let n = 512;
    let dt = 1.0 / 8000.0;
    let f0 = 132.7;
    let plan = HarmonicFitPlan::new(n, dt, f0, &[1, 2, 3], budget()).unwrap();
    let record = harmonic_record(n, dt, f0);
    let (force, summary) = fit_force_harmonics(&record, dt, &plan).unwrap();
    let expected = [
        [C64::new(1.3, -0.2), C64::new(0.0, 0.3)],
        [C64::new(0.00001, 0.00002), C64::new(-0.07, 0.0)],
        [C64::new(-0.15, 0.0), C64::new(0.0, 0.25)],
    ];
    for (actual, expected) in force.iter().zip(expected) {
        for axis in 0..2 {
            assert!((actual[axis] - expected[axis]).abs() < 3e-12);
        }
    }
    close(summary[0].dc, 0.4, 3e-12);
    close(summary[1].dc, -0.8, 3e-12);
    for fit in &summary {
        assert!(fit.residual_rms < 3e-12);
    }
    let acoustic = FluidAcousticState::try_new(1.2, 340.0).unwrap();
    let fitted = JetForceHarmonics {
        acoustic,
        harmonics: plan.harmonics().to_vec(),
        fundamental_hz: f0,
        force_n_m: force,
        dc_force_n_m: [summary[0].dc, summary[1].dc],
        residual_rms_n_m: [summary[0].residual_rms, summary[1].residual_rms],
    };
    let observed = fitted.observe_harmonic(2, [0.4, 0.3], [0.0, 0.0]).unwrap();
    assert_eq!(observed.frequency_hz, 2.0 * f0);
    let reference =
        dipole_pressure(expected[1], TAU * 2.0 * f0 / 340.0, [0.4, 0.3], [0.0, 0.0]).unwrap();
    assert!((observed.field.pressure - reference.pressure).abs() < 1e-11);
    assert_eq!(observed.field.scope, SCOPE_STATEMENT);
    assert!(fitted.observe_harmonic(4, [0.4, 0.3], [0.0, 0.0]).is_err());
}

#[test]
fn g0_force_spectrum_budget_bins_and_original_fit_clock_refuse() {
    let forces = harmonic_record(64, 0.001, 31.7);
    assert!(ForceSpectrum::new(&forces, 0.001, 63).is_err());
    let odd = ForceSpectrum::new(&forces[..63], 0.001, 64).unwrap();
    assert!(odd.line(31).is_ok());
    assert!(odd.line(32).is_err());
    assert!(ForceSpectrum::new(&forces, 0.0, 64).is_err());
    let spectrum = ForceSpectrum::new(&forces, 0.001, 64).unwrap();
    for bin in [0, 32, usize::MAX] {
        assert!(spectrum.line(bin).is_err());
    }
    let plan = HarmonicFitPlan::new(64, 0.001, 31.7, &[1, 2, 3], budget()).unwrap();
    assert!(fit_force_harmonics(&forces[..63], 0.001, &plan).is_err());
    assert!(fit_force_harmonics(&forces, 0.002, &plan).is_err());
    let mut invalid = forces;
    invalid[5][1] = f64::NAN;
    assert!(matches!(
        fit_force_harmonics(&invalid, 0.001, &plan),
        Err(AeroacError::HarmonicFit(HarmonicFitError::NonFinite { .. }))
    ));
}

#[test]
fn g3_joint_force_fit_preserves_phase_and_physical_force_scaling() {
    let n = 256;
    let dt = 1.0 / 8000.0;
    let plan = HarmonicFitPlan::new(n, dt, 132.7, &[1, 2, 3], budget()).unwrap();
    let record = harmonic_record(n, dt, 132.7);
    let scale = 0.037;
    let scaled: Vec<_> = record.iter().map(|f| f.map(|x| x * scale)).collect();
    let (original, baseline) = fit_force_harmonics(&record, dt, &plan).unwrap();
    let (changed, summary) = fit_force_harmonics(&scaled, dt, &plan).unwrap();
    for (actual, expected) in changed.iter().zip(original) {
        for axis in 0..2 {
            assert!((actual[axis] - expected[axis].scale(scale)).abs() < 3e-14);
        }
    }
    for axis in 0..2 {
        close(summary[axis].dc, baseline[axis].dc * scale, 3e-14);
    }
}

#[test]
fn g1_noncommensurate_force_tones_preserve_weak_lines_and_true_frequency_radiation() {
    let n = 1024;
    let dt = 1.0 / 4000.0;
    // First two tones are closer than one FFT bin; the third is inharmonic.
    let frequencies = [132.7, 134.1, 267.2578];
    let amplitudes = [
        [C64::new(1.3, -0.2), C64::new(0.0, 0.3)],
        [C64::new(0.00001, 0.00002), C64::new(0.025, -0.07)],
        [C64::new(-0.15, 0.04), C64::new(-0.03, 0.25)],
    ];
    let dc = [0.4, -0.8];
    let record: Vec<_> = (0..n)
        .map(|i| {
            let mut force = dc;
            for (frequency, phasors) in frequencies.iter().zip(amplitudes) {
                let theta = TAU * frequency * i as f64 * dt;
                for axis in 0..2 {
                    force[axis] += phasors[axis].re * theta.cos() + phasors[axis].im * theta.sin();
                }
            }
            force
        })
        .collect();
    let plan = ToneFitPlan::new(n, dt, &frequencies, budget()).unwrap();
    let (force, summary) = fit_force_tones(&record, dt, &plan).unwrap();
    for (actual, expected) in force.iter().zip(amplitudes) {
        for axis in 0..2 {
            assert!((actual[axis] - expected[axis]).abs() < 3e-11);
        }
    }
    for axis in 0..2 {
        close(summary[axis].dc, dc[axis], 3e-12);
        assert!(summary[axis].residual_rms < 3e-12);
    }
    let acoustic = FluidAcousticState::try_new(1.2, 340.0).unwrap();
    let fitted = JetForceTones {
        acoustic,
        frequencies_hz: plan.frequencies_hz().to_vec(),
        force_n_m: force,
        dc_force_n_m: summary.map(|s| s.dc),
        residual_rms_n_m: summary.map(|s| s.residual_rms),
    };
    assert_eq!(fitted.frequencies_hz(), frequencies);
    let field_budget = crate::observation::HarmonicObservationBudget {
        max_lines: 3,
        max_line_pairs: 9,
    };
    let field = fitted
        .prepare_observer_field([0.4, 0.3], [0.0, 0.0], field_budget)
        .unwrap();
    let clock = crate::observation::HarmonicSampleClock::new(fs_math::dd::Dd::ZERO, 4000).unwrap();
    let mut pressure = [0.0; 73];
    field
        .sample_pressure_into(
            clock,
            19,
            &mut pressure,
            crate::observation::HarmonicSampleBudget {
                max_samples: 73,
                max_tone_samples: 219,
            },
        )
        .unwrap();
    for (offset, actual) in pressure.iter().enumerate() {
        let time = (19 + offset) as f64 / 4000.0;
        let expected: f64 = frequencies
            .iter()
            .zip(amplitudes)
            .map(|(&frequency, force)| {
                let p = dipole_pressure(force, TAU * frequency / 340.0, [0.4, 0.3], [0.0, 0.0])
                    .unwrap()
                    .pressure;
                let theta = TAU * frequency * time;
                p.re * theta.cos() + p.im * theta.sin()
            })
            .sum();
        assert!((actual - expected).abs() < 2e-10);
    }
    assert!(matches!(
        fitted.prepare_observer_field_cancellable([0.4, 0.3], [0.0, 0.0], field_budget, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    for index in 0..frequencies.len() {
        let observed = fitted.observe_tone(index, [0.4, 0.3], [0.0, 0.0]).unwrap();
        assert_eq!(observed.frequency_hz, frequencies[index]);
        let expected = dipole_pressure(
            amplitudes[index],
            TAU * frequencies[index] / 340.0,
            [0.4, 0.3],
            [0.0, 0.0],
        )
        .unwrap();
        assert!((observed.field.pressure - expected.pressure).abs() < 1e-10);
        assert_eq!(observed.field.scope, SCOPE_STATEMENT);
    }
    assert!(
        fitted
            .observe_tone(usize::MAX, [0.4, 0.3], [0.0, 0.0])
            .is_err()
    );
}

#[test]
fn g0_g3_tone_force_fit_preserves_harmonic_results_and_original_clock_refusals() {
    let n = 128;
    let dt = 1.0 / 8000.0;
    let f0 = 132.7;
    let record = harmonic_record(n, dt, f0);
    let harmonic = HarmonicFitPlan::new(n, dt, f0, &[1, 2, 3], budget()).unwrap();
    let tones = ToneFitPlan::new(n, dt, &[f0, 2.0 * f0, 3.0 * f0], budget()).unwrap();
    let (h, hs) = fit_force_harmonics(&record, dt, &harmonic).unwrap();
    let (t, ts) = fit_force_tones(&record, dt, &tones).unwrap();
    for (a, b) in h.iter().zip(t) {
        for axis in 0..2 {
            assert!((a[axis] - b[axis]).abs() < 3e-12);
        }
    }
    for axis in 0..2 {
        close(hs[axis].dc, ts[axis].dc, 3e-12);
    }
    assert!(fit_force_tones(&record[..n - 1], dt, &tones).is_err());
    assert!(fit_force_tones(&record, 2.0 * dt, &tones).is_err());
    let mut invalid = record;
    invalid[17][1] = f64::INFINITY;
    assert!(matches!(
        fit_force_tones(&invalid, dt, &tones),
        Err(AeroacError::HarmonicFit(HarmonicFitError::NonFinite { .. }))
    ));
}

#[test]
fn g0_g1_recorded_force_cross_spectrum_preserves_medium_clock_and_field_statistics() {
    use crate::curle2d::CoherenceBudget;
    use fs_fft::welch::{
        WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchError, WelchWindow,
    };

    let gas = GasState::try_new_moist_air(307.0, 83_000.0, 0.4).unwrap();
    let acoustic = gas.acoustic_state().unwrap();
    let dt = 1.0 / 6000.0;
    // A manufactured accepted history at the observation seam, as in the
    // refined-tone tests. This does not assert convergence of a jet simulation.
    let run = MaterialJetLabiumRun {
        material: &gas,
        acoustic,
        timestep_s: dt,
        force_series_n_m: harmonic_record(265, dt, 282.6),
        sampled_physical_mach_max: 0.0,
        lattice_diagnostics: JetLabiumDiagnostics {
            mach_max_lattice: 0.0,
            flux_plate_plane: 0.0,
            flux_fringe_plane: 0.0,
            reynolds: 0.0,
        },
    };
    let config = WelchConfig {
        segment_samples: 64,
        hop_samples: 32,
        window: WelchWindow::PeriodicHann,
        detrend: WelchDetrend::SegmentMean,
    };
    let work = WelchBudget {
        maximum_channels: 2,
        maximum_samples_per_channel: 265,
        maximum_transforms: 14,
        maximum_workspace_bytes: 1024 * 1024,
    };
    let prepared = run
        .prepare_force_cross_spectrum(config, work, || false)
        .unwrap();
    assert_eq!(prepared.acoustic(), acoustic);
    let spectrum = prepared.spectrum();
    assert_eq!(spectrum.timestep_s(), dt);
    assert_eq!(spectrum.config(), config);
    assert_eq!(spectrum.channel_count(), 2);
    assert_eq!(spectrum.segment_count(), 7);
    assert_eq!(spectrum.bin_count(), 33);
    assert_eq!(spectrum.used_sample_count(), 256);
    assert_eq!(spectrum.discarded_tail_samples(), 9);
    close(spectrum.bin_width_hz(), 6000.0 / 64.0, 0.0);

    // Independently submit the same physical history as separate x/y channels.
    // This catches channel/clock wiring mistakes in the interleaved adapter.
    let x: Vec<_> = run.force_series_n_m().iter().map(|row| row[0]).collect();
    let y: Vec<_> = run.force_series_n_m().iter().map(|row| row[1]).collect();
    let reference = WelchCrossSpectrum::estimate(&[&x, &y], dt, config, work, || false).unwrap();
    let receiver = [0.37, -0.24];
    let source = [-0.03, 0.06];
    let coherence = CoherenceBudget {
        max_modes: 7,
        max_source_mode_pairs: 7,
    };
    for bin in [3, 6, 9] {
        let actual = prepared
            .observe_bin(bin, receiver, source, coherence, || false)
            .unwrap();
        let observer = PreparedDipoleObserver::new(
            acoustic,
            reference.frequency_hz(bin).unwrap(),
            receiver,
            &[source],
            1,
        )
        .unwrap();
        let expected = observer
            .evaluate_welch_bin(&reference, bin, coherence, || false)
            .unwrap();
        assert!(actual.mean_square_pressure_pa2 > 0.0);
        assert_eq!(
            actual.mean_square_pressure_pa2,
            expected.mean_square_pressure_pa2
        );
        assert_eq!(
            actual.velocity_covariance_m2_s2,
            expected.velocity_covariance_m2_s2
        );
        assert_eq!(actual.active_intensity_w_m2, expected.active_intensity_w_m2);
        assert_eq!(
            actual.reactive_intensity_w_m2,
            expected.reactive_intensity_w_m2
        );
        assert_eq!(actual.scope, SCOPE_STATEMENT);
    }
    for bin in [0, spectrum.bin_count() - 1, usize::MAX] {
        assert!(matches!(
            prepared.observe_bin(bin, receiver, source, coherence, || false),
            Err(AeroacError::InvalidParameter { .. })
        ));
    }
    for limited in [
        CoherenceBudget {
            max_modes: 6,
            ..coherence
        },
        CoherenceBudget {
            max_source_mode_pairs: 6,
            ..coherence
        },
    ] {
        assert!(matches!(
            prepared.observe_bin(3, receiver, source, limited, || false),
            Err(AeroacError::InvalidParameter { .. })
        ));
    }
    assert!(matches!(
        run.prepare_force_cross_spectrum(
            config,
            WelchBudget {
                maximum_transforms: 13,
                ..work
            },
            || false,
        ),
        Err(WelchError::BudgetExceeded(_))
    ));
}
