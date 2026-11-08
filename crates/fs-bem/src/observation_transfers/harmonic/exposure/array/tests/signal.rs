use super::*;

const DRAWS: [HarmonicFactorDraw; 2] = [
    HarmonicFactorDraw::Real(-0.7),
    HarmonicFactorDraw::ProperComplex(C64::new(0.3, -1.1)),
];

fn realized(fixture: &Fixture) -> PreparedHarmonicRealization<'_> {
    let banks: [[HarmonicExposureLine<'_>; 2]; 3] = core::array::from_fn(|b| {
        [
            line(&fixture.owners[0], &fixture.amplitudes[b][..1]),
            line(&fixture.owners[1], &fixture.amplitudes[b][1..]),
        ]
    });
    PreparedHarmonicExposure::new_with_statistics(
        &fixture.owners[0],
        "directional-waveform-clock",
        &banks[0],
        &[&banks[1], &banks[2]],
        &[
            HarmonicExposureFactorStatistics::Real,
            HarmonicExposureFactorStatistics::ProperComplex,
        ],
        preparation(),
    )
    .unwrap()
    .realize(&DRAWS, realization_limits())
    .unwrap()
}

fn realization_limits() -> HarmonicRealizationLimits {
    HarmonicRealizationLimits {
        maximum_lines: 32,
        maximum_field_products: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn sampling(frames: usize) -> HarmonicPressureSampleLimits {
    HarmonicPressureSampleLimits {
        maximum_frames: frames,
        maximum_channel_line_samples: 6 * frames,
    }
}

fn directional_fixture() -> Fixture {
    let mut fixture = Fixture::new();
    for (frequency, owner) in fixture.owners.iter_mut().enumerate() {
        let z = owner.medium.density * owner.medium.sound_speed;
        for (receiver, field) in owner.transfers.iter_mut().enumerate() {
            for c in 0..3 {
                field[c + 1] = C64::new(
                    (0.3 + 0.2 * receiver as f64 - 0.1 * c as f64) / z,
                    (-0.6 + 0.4 * frequency as f64 + 0.3 * c as f64) / z,
                );
            }
        }
    }
    fixture
}

fn patterns() -> [FirstOrder; 2] {
    [
        FirstOrder::new(0.25, [0.6, 0.8, 0.0]).unwrap(),
        FirstOrder::new(0.7, [-2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0]).unwrap(),
    ]
}

// Independent real-time delayed p/u projection; no prepared channel phasors.
fn pressure(
    fixture: &Fixture,
    patterns: &[FirstOrder],
    channel: &HarmonicPressureChannel,
    time: f64,
) -> f64 {
    let mut value = 0.0;
    for tap in &channel.taps {
        for frequency in 0..2 {
            let owner = &fixture.owners[frequency];
            let coefficient = fixture.amplitudes[0][frequency]
                + fixture.amplitudes[1][frequency].scale(-0.7)
                + fixture.amplitudes[2][frequency] * C64::new(0.3, -1.1);
            let angle = owner.omega * (time - tap.delay_s);
            let (sine, cosine) = angle.sin_cos();
            let field = owner.transfers[tap.receiver].map(|p| {
                let amplitude = p * coefficient;
                amplitude.re * cosine + amplitude.im * sine
            });
            let pattern = patterns[tap.receiver];
            let axial: f64 = pattern
                .front_axis()
                .iter()
                .zip(&field[1..])
                .map(|(n, u)| n * u)
                .sum();
            value += tap.weight
                * (pattern.pressure_fraction() * field[0]
                    - (1.0 - pattern.pressure_fraction())
                        * owner.medium.density
                        * owner.medium.sound_speed
                        * axial);
        }
    }
    value
}

#[test]
fn directional_signal_uses_the_same_draw_actual_velocity_and_fractional_delays() {
    let fixture = directional_fixture();
    let realization = realized(&fixture);
    let patterns = patterns();
    let channels = channels();
    let signal =
        PreparedHarmonicPressureSignal::new(&realization, &channels, Some(&patterns), limits())
            .unwrap();
    assert!(core::ptr::eq(signal.realization(), &realization));
    assert_eq!(signal.realization().factor_draws(), &DRAWS);
    assert_eq!(signal.pressure_array().channels(), &channels);
    assert_eq!(
        signal.pressure_array().microphone_patterns(),
        Some(patterns.as_slice())
    );
    assert_eq!(signal.pressure_array().factor_statistics().count(), 0);
    let clock = HarmonicSampleClock::new(Dd::from_f64(0.017), 8192).unwrap();
    let mut output = vec![0.0; 3 * 33];
    signal
        .sample_pressure_into(clock, 0, &mut output, sampling(33))
        .unwrap();
    for (frame, values) in output.chunks_exact(3).enumerate() {
        for (channel, &value) in channels.iter().zip(values) {
            close(
                value,
                pressure(&fixture, &patterns, channel, 0.017 + frame as f64 / 8192.0),
                3e-13,
            );
        }
    }
    let mut blocks = vec![0.0; output.len()];
    let mut first = 0;
    for block in blocks.chunks_mut(3 * 7) {
        signal
            .sample_pressure_into(clock, first, block, sampling(33))
            .unwrap();
        first += (block.len() / 3) as u64;
    }
    assert_eq!(
        output.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        blocks.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
    );
}

#[test]
fn sampled_waveform_moments_match_the_retained_deterministic_array() {
    let fixture = directional_fixture();
    let realization = realized(&fixture);
    let signal =
        PreparedHarmonicPressureSignal::new(&realization, &channels(), Some(&patterns()), limits())
            .unwrap();
    let n = 1024usize;
    let start = 0.13;
    let clock = HarmonicSampleClock::new(Dd::from_f64(start), 4096).unwrap();
    let window = HarmonicWindow::new(start, n as f64 / 4096.0).unwrap();
    let analytic = signal.pressure_array().observe_window(window).unwrap();
    let mut output = vec![0.0; 3 * (n + 1)];
    signal
        .sample_pressure_into(clock, 0, &mut output, sampling(n + 1))
        .unwrap();
    let mut means = [0.0; 3];
    let mut squares = [0.0; 3];
    for (i, frame) in output.chunks_exact(3).enumerate() {
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for c in 0..3 {
            means[c] += weight * frame[c] / (3.0 * n as f64);
            squares[c] += weight * frame[c] * frame[c] / (3.0 * n as f64);
        }
    }
    for (c, point) in analytic.moments().iter().enumerate() {
        close(means[c], point.mean_pressure_pa, 3e-11);
        close(squares[c], point.mean_square_pressure_pa2, 3e-11);
    }
    assert!(analytic.averaged_covariance_pa2().iter().all(|&v| v == 0.0));
}

#[test]
fn directional_cancellation_retains_a_residual_lost_by_projecting_rounded_samples() {
    let mut reference = owner(7.0, &[C64::ONE]);
    reference.medium.density = 1.0;
    reference.medium.sound_speed = 1.0;
    reference.transfers[0][1] = C64::ONE;
    reference
        .sources
        .push(SourceId::StructuralMode("weak-pressure".into()));
    reference
        .transfers
        .push([C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]);
    let amplitudes = [C64::from_re(1e16), C64::ONE];
    let exposure = PreparedHarmonicExposure::new(
        &reference,
        "quiet-residual",
        &[line(&reference, &amplitudes)],
        &[],
        preparation(),
    )
    .unwrap();
    let realization = exposure.realize(&[], realization_limits()).unwrap();
    let clock = HarmonicSampleClock::new(Dd::ZERO, 1024).unwrap();
    let mut physical = [HarmonicFieldSample::default()];
    realization
        .sample_fields_into(
            clock,
            0,
            &mut physical,
            HarmonicSampleLimits {
                maximum_frames: 1,
                maximum_receiver_line_samples: 1,
            },
        )
        .unwrap();
    assert_eq!(
        physical[0].pressure_pa,
        physical[0].particle_velocity_m_s[0]
    );
    let signal = PreparedHarmonicPressureSignal::new(
        &realization,
        &[
            channel("directional", &[tap(0, 1.0, 0.0)]),
            channel("null", &[tap(0, 1.0, 0.0), tap(0, -1.0, 0.0)]),
        ],
        Some(&[FirstOrder::new(0.5, [1.0, 0.0, 0.0]).unwrap()]),
        limits(),
    )
    .unwrap();
    let mut output = [0.0; 2];
    signal
        .sample_pressure_into(clock, 0, &mut output, sampling(1))
        .unwrap();
    assert_eq!(output, [0.5, 0.0]);
}

#[test]
fn dark_and_empty_channels_keep_bandwidth_budgets_and_atomic_admission() {
    let fixture = directional_fixture();
    let realization = realized(&fixture);
    let definitions = [channel("dark", &[])];
    let signal =
        PreparedHarmonicPressureSignal::new(&realization, &definitions, None, limits()).unwrap();
    let mut output = [17.0; 4];
    let aliased = HarmonicSampleClock::new(Dd::ZERO, 1).unwrap();
    assert!(
        signal
            .sample_pressure_into(aliased, 0, &mut output, sampling(4))
            .is_err()
    );
    assert_eq!(output, [17.0; 4]);
    let clock = HarmonicSampleClock::new(Dd::ZERO, 4096).unwrap();
    assert!(
        signal
            .sample_pressure_into(
                clock,
                0,
                &mut output,
                HarmonicPressureSampleLimits {
                    maximum_channel_line_samples: 0,
                    ..sampling(4)
                }
            )
            .is_err()
    );
    assert!(
        signal
            .sample_pressure_into(clock, u64::MAX, &mut output, sampling(4))
            .is_err()
    );
    assert!(matches!(
        signal.sample_pressure_into_cancellable(clock, 0, &mut output, sampling(4), || true),
        Err(TransferError::Cancelled)
    ));
    assert_eq!(output, [17.0; 4]);
    signal
        .sample_pressure_into(clock, 0, &mut output, sampling(4))
        .unwrap();
    assert_eq!(output, [0.0; 4]);
    assert!(matches!(
        PreparedHarmonicPressureSignal::new_cancellable(
            &realization,
            &definitions,
            None,
            limits(),
            || true
        ),
        Err(TransferError::Cancelled)
    ));
    let empty = PreparedHarmonicPressureSignal::new(&realization, &[], None, limits()).unwrap();
    assert!(
        empty
            .sample_pressure_into(clock, 0, &mut output, sampling(4))
            .is_err()
    );
    empty
        .sample_pressure_into(clock, 0, &mut [], sampling(0))
        .unwrap();
    assert!(
        empty
            .sample_pressure_into(aliased, 0, &mut [], sampling(0))
            .is_err()
    );
}
