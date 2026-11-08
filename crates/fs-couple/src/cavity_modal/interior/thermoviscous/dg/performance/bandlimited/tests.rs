use super::*;

fn policy() -> BandlimitedObservationPolicy {
    BandlimitedObservationPolicy {
        passband_hz: 4000.0,
        stopband_hz: 12000.0,
        phase_count: 128,
        max_preparation_terms: 1 << 22,
        max_history_values: 1 << 16,
        max_output_terms: 1 << 16,
    }
}

#[test]
fn selected_pressure_and_every_terminal_share_original_mean_gain_and_delay() {
    let dt: f64 = 7.3e-6;
    let frequency = 1800.0;
    let omega = core::f64::consts::TAU * frequency;
    let mean_gain = (0.5 * omega * dt).sin() / (0.5 * omega * dt);
    for rate in [48000, 192000] {
        let count = 512;
        let sources = (count as f64 / (dt * f64::from(rate))).ceil() as usize;
        let gate = CancelGate::new_clock_free();
        let mut observer =
            BandlimitedDgObservation::new(dt, rate, count, sources, 2, policy(), &gate).unwrap();
        let report = observer.resampler().report();
        let mut accepted = 0;
        let mut flows = [0.0; 2];
        let mut peak = 0.0_f64;
        for sample in 0..count {
            let required = observer.required_source_intervals().unwrap();
            while accepted < required {
                let angle = omega * (accepted as f64 + 0.5) * dt;
                observer
                    .push_source(
                        mean_gain * angle.sin(),
                        &[
                            1e-4 * mean_gain * angle.cos(),
                            -2e-4 * mean_gain * angle.sin(),
                        ],
                    )
                    .unwrap();
                accepted += 1;
            }
            let pressure = observer.render_next(&mut flows).unwrap();
            let time = (sample as f64 + 0.5) / f64::from(rate);
            if time > 2.0 * report.latency_s + 2.0 * dt {
                let angle = omega * (time - report.latency_s);
                assert!((pressure - mean_gain * angle.sin()).abs() < 1e-5);
                assert!((flows[0] / 1e-4 - mean_gain * angle.cos()).abs() < 1e-5);
                assert!((flows[1] / -2e-4 - mean_gain * angle.sin()).abs() < 1e-5);
                peak = peak.max(pressure.abs());
            }
        }
        assert!(peak > 0.9, "must reach settled audible pressure");
        assert_eq!(accepted, sources);
        assert_eq!(observer.resampler().rendered_samples(), count);
        assert_eq!(
            observer.resampler().report().input_period_s.to_bits(),
            dt.to_bits()
        );
        assert!(observer.required_source_intervals().is_err());
        assert!(observer.push_source(0.0, &[0.0; 2]).is_err());
    }
}

#[test]
fn finite_control_count_and_original_channel_basis_are_admitted_exactly() {
    let dt = 1.0 / 96000.0;
    let gate = CancelGate::new_clock_free();
    for wrong in [28, 30] {
        assert!(BandlimitedDgObservation::new(dt, 44100, 13, wrong, 2, policy(), &gate).is_err());
    }
    let mut observer =
        BandlimitedDgObservation::new(dt, 44100, 13, 29, 2, policy(), &gate).unwrap();
    assert!(observer.push_source(1.0, &[2.0]).is_err());
    assert!(observer.push_source(1.0, &[2.0, f64::NAN]).is_err());
    assert_eq!(observer.resampler().accepted_source_intervals(), 0);
    let mut untouched = [8.0, 9.0];
    assert!(observer.render_next(&mut untouched).is_err());
    assert_eq!(untouched, [8.0, 9.0]);
    let required = observer.required_source_intervals().unwrap();
    for _ in 0..required {
        observer.push_source(1.0, &[2.0, 3.0]).unwrap();
    }
    assert!(observer.render_next(&mut [0.0]).is_err());
    assert_eq!(observer.resampler().rendered_samples(), 0);
    observer.render_next(&mut untouched).unwrap();
    assert_eq!(observer.resampler().rendered_samples(), 1);
}

#[test]
fn pressure_only_observation_and_cancellation_keep_their_existing_meanings() {
    let gate = CancelGate::new_clock_free();
    let mut observer =
        BandlimitedDgObservation::new(1.0 / 48000.0, 48000, 1, 1, 0, policy(), &gate).unwrap();
    observer.push_source(1.0, &[]).unwrap();
    assert!(observer.render_next(&mut []).unwrap().is_finite());
    assert_eq!(observer.resampler().report().channels, 1);
    gate.request();
    assert!(matches!(
        BandlimitedDgObservation::new(1.0 / 48000.0, 48000, 1, 1, 0, policy(), &gate),
        Err(RenderError::Voice(AcousticRealizeError::Cancelled))
    ));
}
