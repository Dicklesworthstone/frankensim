use super::*;

fn field(lines: &[(f64, C64)]) -> PreparedHarmonicField<2> {
    PreparedHarmonicField::new(
        &lines
            .iter()
            .map(|&(frequency_hz, pressure_pa)| HarmonicFieldLine {
                frequency_hz,
                pressure_pa,
                particle_velocity_m_s: [C64::ZERO; 2],
            })
            .collect::<Vec<_>>(),
        HarmonicObservationBudget {
            max_lines: 16,
            max_line_pairs: 256,
        },
    )
    .unwrap()
}

fn budget(samples: usize) -> HarmonicSampleBudget {
    HarmonicSampleBudget {
        max_samples: samples,
        max_tone_samples: samples * 16,
    }
}

#[test]
fn g1_sampled_off_grid_beats_match_independent_real_waveform() {
    let tones = [
        (213.17, C64::new(0.7, -0.2)),
        (217.83, C64::new(-0.11, 0.3)),
        (431.29, C64::new(0.06, 0.04)),
    ];
    let field = field(&tones);
    let clock = HarmonicSampleClock::new(Dd::from_f64(-0.031), 8000).unwrap();
    let mut output = [0.0; 257];
    field
        .sample_pressure_into(clock, 79, &mut output, budget(257))
        .unwrap();
    for (index, actual) in output.iter().enumerate() {
        let t = -0.031 + (79 + index) as f64 / 8000.0;
        let expected: f64 = tones
            .iter()
            .map(|(f, a)| {
                let angle = core::f64::consts::TAU * f * t;
                a.re * angle.cos() + a.im * angle.sin()
            })
            .sum();
        assert!((actual - expected).abs() < 4e-14);
    }
}

#[test]
fn g3_callback_partition_and_long_integer_clock_preserve_samples() {
    let field = field(&[
        (239.317, C64::new(0.2, -0.9)),
        (240.1, C64::new(0.07, 0.04)),
    ]);
    let clock = HarmonicSampleClock::new(Dd::from_f64(0.137), 48_000).unwrap();
    let mut complete = [0.0; 99];
    let mut partitioned = [0.0; 99];
    field
        .sample_pressure_into(clock, 123, &mut complete, budget(99))
        .unwrap();
    let mut first = 123;
    for chunk in partitioned.chunks_mut(7) {
        field
            .sample_pressure_into(clock, first, chunk, budget(7))
            .unwrap();
        first += chunk.len() as u64;
    }
    assert_eq!(complete, partitioned);
    let long = (1_u64 << 53) + 7;
    let interval = clock.time_at(long + 1).unwrap() - clock.time_at(long).unwrap();
    assert!((interval.to_f64() - 1.0 / 48_000.0).abs() < 2e-20);
    let late = HarmonicSampleClock::new(Dd::from_f64(2_f64.powi(40)), 48_000).unwrap();
    assert_eq!(late.time_at(1).unwrap().hi, late.epoch_s().hi);
    assert!(late.time_at(1).unwrap().lo > 0.0);
}

#[test]
fn g3_integer_cycles_on_large_epochs_preserve_physical_hz_phase() {
    let tones = [
        (1024.0, C64::new(0.7, -0.2)),
        (1536.0, C64::new(-0.11, 0.3)),
        (1030.25, C64::new(0.06, 0.04)),
    ];
    let field = field(&tones);
    let duration = 37.0 / 8192.0;
    let reference = field
        .observe_window(HarmonicWindow::new(0.0, duration).unwrap())
        .unwrap();
    for epoch in [2.0_f64.powi(40), -2.0_f64.powi(40)] {
        // Every frequency completes an integer number of cycles at this
        // epoch. Reducing that integer analytically avoids using the same
        // large-angle Hz conversion as the implementation under test.
        for (frequency, _) in tones {
            assert_eq!((frequency * epoch).fract(), 0.0);
        }
        let clock = HarmonicSampleClock::new(Dd::from_f64(epoch), 8192).unwrap();
        let mut samples = [0.0; 65];
        field
            .sample_pressure_into(clock, 0, &mut samples, budget(65))
            .unwrap();
        for (index, actual) in samples.iter().enumerate() {
            let time = index as f64 / 8192.0;
            let expected: f64 = tones
                .iter()
                .map(|(frequency, amplitude)| {
                    let phase = core::f64::consts::TAU * frequency * time;
                    amplitude.re * phase.cos() + amplitude.im * phase.sin()
                })
                .sum();
            assert!((actual - expected).abs() < 2e-14);
        }
        // Analytic windows must use the same corrected physical Hz convention
        // as sampled output, retaining both coherent sum and beat terms.
        let observed = field
            .observe_window(HarmonicWindow::new(epoch, duration).unwrap())
            .unwrap();
        assert!((observed.mean_pressure_pa - reference.mean_pressure_pa).abs() < 2e-14);
        assert!(
            (observed.mean_square_pressure_pa2 - reference.mean_square_pressure_pa2).abs() < 2e-14
        );
    }
}

#[test]
fn g1_samples_agree_with_analytic_coherent_window_energy() {
    let field = field(&[(73.2, C64::new(1.0, -0.3)), (76.7, C64::new(-0.6, 0.7))]);
    let clock = HarmonicSampleClock::new(Dd::from_f64(0.041), 200_000).unwrap();
    let mut samples = vec![0.0; 4001];
    field
        .sample_pressure_into(clock, 0, &mut samples, budget(4001))
        .unwrap();
    let integral = samples
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let weight = if i == 0 || i == 4000 {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            weight * p * p
        })
        .sum::<f64>()
        / (3.0 * 4000.0);
    let exact = field
        .observe_window(HarmonicWindow::new(0.041, 0.02).unwrap())
        .unwrap();
    assert!((integral - exact.mean_square_pressure_pa2).abs() < 1e-10);
}

#[test]
fn g3_sampling_retains_weak_equal_frequency_cancellation_residual() {
    let field = field(&[
        (100.0, C64::new(1e16, 0.0)),
        (100.0, C64::new(1.0, 0.0)),
        (100.0, C64::new(-1e16, 0.0)),
    ]);
    assert_eq!(field.frequency_count(), 1);
    assert_eq!(field.pressure_at_clock_dd(Dd::ZERO).unwrap(), 1.0);
}

#[test]
fn g0_sampling_refuses_aliasing_clock_work_and_cancellation() {
    let f = field(&[(4000.0, C64::ONE)]);
    let clock = HarmonicSampleClock::new(Dd::ZERO, 8000).unwrap();
    let mut output = [123.0; 8];
    assert!(
        f.sample_pressure_into(clock, 0, &mut output, budget(8))
            .is_err()
    );
    assert_eq!(output, [123.0; 8]);
    let f = field(&[(100.0, C64::ONE)]);
    assert!(f.validate_sample_range(clock, u64::MAX, 1).is_err());
    assert!(HarmonicSampleClock::new(Dd::ZERO, 0).is_err());
    assert!(HarmonicSampleClock::new(Dd::from_f64(1e200), 8000).is_err());
    assert!(
        f.sample_pressure_into(clock, 0, &mut output, budget(7))
            .is_err()
    );
    assert!(matches!(
        f.sample_pressure_into_cancellable(clock, 0, &mut output, budget(8), || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(output, [123.0; 8]);
}
