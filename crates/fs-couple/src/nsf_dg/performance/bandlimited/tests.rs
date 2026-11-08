use super::*;

fn policy() -> BandlimitedObservationPolicy {
    BandlimitedObservationPolicy {
        passband_hz: 6000.0,
        stopband_hz: 15000.0,
        phase_count: 128,
        max_preparation_terms: 1 << 24,
        max_history_values: 1 << 16,
        max_output_terms: 1 << 16,
    }
}

#[test]
fn pressure_terminals_and_distinct_tile_groups_keep_one_physical_phase() {
    // Original interval means of a sinusoid, not endpoint samples. The two
    // independently placed reservoir bases have different sizes and signs.
    // Their coherent scaling must survive the SAME filter before Rayleigh.
    let dt = 1.0 / 96000.0;
    let frequency = 2000.0;
    let omega = core::f64::consts::TAU * frequency;
    let mean_gain = (0.5 * omega * dt).sin() / (0.5 * omega * dt);
    for rate in [72000, 144000] {
        let mut observer =
            BandlimitedObservation::new(dt, rate, 256, 2, [2, 3], policy(), &mut || false).unwrap();
        let latency = observer.resampler().report().latency_s;
        let mut accepted = 0;
        let mut peak = 0.0_f64;
        for sample in 0..256 {
            let required = observer.required_source_intervals().unwrap();
            while accepted < required {
                let value = mean_gain * (omega * (accepted as f64 + 0.5) * dt).sin();
                observer
                    .push_source(
                        value,
                        &[2.0 * value, -value],
                        [
                            &[4.0 * value, 0.0][..],
                            &[-2.0 * value, 8.0 * value, value][..],
                        ],
                    )
                    .unwrap();
                accepted += 1;
            }
            let actual = observer.render_next().unwrap();
            for (channel, scale) in [1.0, 2.0, -1.0, 4.0, 0.0, -2.0, 8.0, 1.0]
                .into_iter()
                .enumerate()
            {
                assert_eq!(
                    actual[channel],
                    scale * actual[0],
                    "channel {channel}, output {sample}"
                );
            }
            let time = (sample as f64 + 0.5) / f64::from(rate);
            if time > 2.0 * latency + 2.0 * dt {
                let expected = mean_gain * (omega * (time - latency)).sin();
                assert!(
                    (actual[0] - expected).abs() < 2.0e-4,
                    "the declared physical mean/latency must survive: {actual:?}, expected {expected}"
                );
                peak = peak.max(actual[0].abs());
            }
        }
        assert!(
            peak > 0.9,
            "the phase check must reach settled nonzero sound"
        );
        assert_eq!(accepted, source_intervals_for(dt, rate, 256).unwrap());
        assert_eq!(observer.resampler().accepted_source_intervals(), accepted);
        assert_eq!(observer.resampler().rendered_samples(), 256);
        assert!(
            observer.required_source_intervals().is_err(),
            "no extra physical tail"
        );
    }
}

#[test]
fn original_channel_basis_is_checked_before_any_filter_history_is_accepted() {
    let mut observer =
        BandlimitedObservation::new(1.0 / 96000.0, 48000, 8, 2, [2, 1], policy(), &mut || false)
            .unwrap();
    for (terminals, reservoirs) in [
        (&[1.0][..], vec![&[3.0, 4.0][..], &[5.0][..]]),
        (&[1.0, 2.0][..], vec![&[3.0, 4.0][..]]),
        (&[1.0, 2.0][..], vec![&[3.0][..], &[4.0, 5.0][..]]),
        (&[1.0, 2.0][..], vec![&[3.0, 4.0][..], &[f64::NAN][..]]),
    ] {
        assert!(observer.push_source(7.0, terminals, reservoirs).is_err());
        assert_eq!(observer.resampler().accepted_source_intervals(), 0);
        assert_eq!(observer.resampler().rendered_samples(), 0);
    }
    observer
        .push_source(7.0, &[1.0, 2.0], [&[3.0, 4.0][..], &[5.0][..]])
        .unwrap();
    assert!(
        observer.render_next().is_err(),
        "the complete physical window is required"
    );
    assert_eq!(observer.resampler().accepted_source_intervals(), 1);
    assert_eq!(observer.resampler().rendered_samples(), 0);
}

#[test]
fn cold_transfer_keeps_exact_clocks_window_channels_and_zero_history() {
    let dt: f64 = 1.0 / 96000.0;
    let mut observer =
        BandlimitedObservation::new(dt, 72000, 16, 1, [2, 3], policy(), &mut || false).unwrap();
    observer
        .validate_cold_basis(dt, 72000, 16, 1, [2, 3])
        .unwrap();
    assert!(
        observer
            .validate_cold_basis(f64::from_bits(dt.to_bits() + 1), 72000, 16, 1, [2, 3])
            .is_err()
    );
    assert!(
        observer
            .validate_cold_basis(dt, 72000, 17, 1, [2, 3])
            .is_err()
    );
    assert!(
        observer
            .validate_cold_basis(dt, 72000, 16, 2, [2, 3])
            .is_err()
    );
    assert!(
        observer
            .validate_cold_basis(dt, 72000, 16, 1, [3, 2])
            .is_err()
    );
    observer
        .push_source(1.0, &[2.0], [&[3.0, 4.0][..], &[5.0, 6.0, 7.0][..]])
        .unwrap();
    assert!(
        observer
            .validate_cold_basis(dt, 72000, 16, 1, [2, 3])
            .is_err()
    );
}
