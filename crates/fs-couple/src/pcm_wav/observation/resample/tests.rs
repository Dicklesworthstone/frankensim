use super::*;

fn policy(dt: f64, rate: u32, phases: usize) -> BandlimitedObservationPolicy {
    let slow = dt.recip().min(f64::from(rate));
    BandlimitedObservationPolicy {
        passband_hz: 0.4 * slow,
        stopband_hz: 0.5 * slow,
        phase_count: phases,
        max_preparation_terms: 1 << 27,
        max_history_values: 1 << 22,
        max_output_terms: 1 << 20,
    }
}

fn sinc(x: f64) -> f64 {
    if x == 0.0 { 1.0 } else { x.sin() / x }
}

// Exact integral mean of an independently specified continuous sinusoid.
// This factor must remain in the answer: observation must not deconvolve it.
fn mean_tone(index: usize, dt: f64, frequency: f64, phase: f64) -> f64 {
    sinc(PI * frequency * dt) * (2.0 * PI * frequency * (index as f64 + 0.5) * dt + phase).sin()
}

fn render(
    dt: f64,
    rate: u32,
    count: usize,
    phases: usize,
    channels: usize,
    mut input: impl FnMut(usize, &mut [f64]),
    chunks: &[usize],
) -> (Vec<f64>, BandlimitedResampler) {
    let mut filter =
        BandlimitedResampler::new(dt, rate, channels, count as u64, policy(dt, rate, phases))
            .unwrap();
    let mut frame = vec![0.0; channels];
    let mut output = vec![0.0; count * channels];
    let mut rendered = 0;
    let mut block = 0;
    while rendered < count {
        let end = count.min(rendered + chunks[block % chunks.len()]);
        for n in rendered..end {
            while filter.accepted_source_intervals() < filter.required_source_intervals().unwrap() {
                input(filter.accepted_source_intervals(), &mut frame);
                filter.push_frame(&frame).unwrap();
            }
            filter
                .render_next(&mut output[n * channels..(n + 1) * channels])
                .unwrap();
        }
        rendered = end;
        block += 1;
    }
    assert_eq!(filter.rendered_samples(), count as u64);
    assert_eq!(
        filter.accepted_source_intervals(),
        filter.source_intervals()
    );
    (output, filter)
}

fn settled_start(report: BandlimitedObservationReport) -> usize {
    // Entire finite stencil lies after t=0, with a modest edge margin.
    ((2 * report.half_support + 3) as f64
        * report.input_period_s
        * f64::from(report.output_rate_hz))
    .ceil() as usize
        + 4
}

#[test]
fn fractional_rates_retain_interval_mean_gain_phase_and_si_channel_ratios() {
    for (dt, rate) in [
        (1.0_f64 / 44100.0, 48000),
        (1.0 / 48000.0, 44100),
        (7.3e-6, 48000),
        (1.0 / 32000.0, 192000),
    ] {
        let frequency = 0.395 * dt.recip().min(f64::from(rate));
        let (output, filter) = render(
            dt,
            rate,
            2048,
            1024,
            3,
            |i, frame| {
                frame[0] = mean_tone(i, dt, frequency, 0.0);
                frame[1] = 2.5e-5 * mean_tone(i, dt, frequency, PI / 2.0);
                frame[2] = -17.0;
            },
            &[1, 13, 256, 3],
        );
        let info = filter.report();
        assert_eq!(info.input_period_s.to_bits(), dt.to_bits());
        assert_eq!(info.latency_s, (info.half_support as f64 + 0.5) * dt);
        let mut maximum_error = 0.0_f64;
        for n in settled_start(info)..2048 {
            let time = (n as f64 + 0.5) / f64::from(rate) - info.latency_s;
            let angle = 2.0 * PI * frequency * time;
            let gain = sinc(PI * frequency * dt);
            maximum_error = maximum_error.max((output[3 * n] - gain * angle.sin()).abs());
            maximum_error =
                maximum_error.max((output[3 * n + 1] / 2.5e-5 - gain * angle.cos()).abs());
            assert!((output[3 * n + 2] + 17.0).abs() < 1e-12);
        }
        assert!(
            maximum_error < 7e-6,
            "dt={dt}, rate={rate}, error={maximum_error}"
        );
        // The original physical box-mean attenuation is material at this tone.
        assert!(sinc(PI * frequency * dt) < 0.98);
    }
}

// Independent definition of the legacy conservative projection: integrate
// piecewise-constant original source means over each output time interval.
fn overlap_tone(dt: f64, rate: u32, n: usize, frequency: f64) -> f64 {
    let begin = n as f64 / f64::from(rate);
    let end = (n + 1) as f64 / f64::from(rate);
    let first = (begin / dt).floor() as usize;
    let last = (end / dt).ceil() as usize;
    (first..last)
        .map(|i| {
            let length = end.min((i + 1) as f64 * dt) - begin.max(i as f64 * dt);
            length.max(0.0) * f64::from(rate) * mean_tone(i, dt, frequency, 0.0)
        })
        .sum()
}

#[test]
fn removes_audible_foldover_that_conservative_overlap_retains() {
    let dt = 1.0 / 96000.0;
    let rate = 44100;
    let frequency = 30000.0; // Original resolved tone folds to 14.1 kHz.
    let count = 4096;
    let (output, filter) = render(
        dt,
        rate,
        count,
        1024,
        1,
        |i, frame| frame[0] = mean_tone(i, dt, frequency, 0.0),
        &[73],
    );
    let start = settled_start(filter.report());
    let filtered_rms =
        (output[start..].iter().map(|x| x * x).sum::<f64>() / (count - start) as f64).sqrt();
    let overlap_rms = ((start..count)
        .map(|i| overlap_tone(dt, rate, i, frequency).powi(2))
        .sum::<f64>()
        / (count - start) as f64)
        .sqrt();
    println!(
        "30 kHz -> 14.1 kHz foldover: overlap RMS={overlap_rms:.9e}, bandlimited RMS={filtered_rms:.9e}"
    );
    assert!(overlap_rms > 0.1);
    assert!(filtered_rms < 4e-6);
    assert!(overlap_rms / filtered_rms > 25000.0);
}

#[test]
fn stop_edge_and_near_input_nyquist_remain_rejected_at_noninteger_ratios() {
    for (dt, rate) in [
        (1.0 / 48000.0, 44100),
        (7.3e-6, 48000),
        (1.0 / 96000.0, 44100),
    ] {
        let stop = policy(dt, rate, 1024).stopband_hz;
        for frequency in [stop, stop * 1.001, 0.499 * dt.recip()] {
            let (output, filter) = render(
                dt,
                rate,
                1024,
                1024,
                1,
                |i, frame| frame[0] = mean_tone(i, dt, frequency, 0.31),
                &[97],
            );
            let maximum = output[settled_start(filter.report())..]
                .iter()
                .fold(0.0_f64, |m, x| m.max(x.abs()));
            assert!(
                maximum < 7e-6,
                "dt={dt}, rate={rate}, f={frequency}, leakage={maximum}"
            );
        }
    }
}

#[test]
fn fractional_phase_resolution_converges_separately_from_the_window() {
    let dt = 1.0 / 44100.0;
    let rate = 48000;
    let frequency = 17000.0;
    let mut errors = Vec::new();
    for phases in [16, 64, 256] {
        let (output, filter) = render(
            dt,
            rate,
            1024,
            phases,
            1,
            |i, frame| frame[0] = mean_tone(i, dt, frequency, 0.0),
            &[23, 71],
        );
        let info = filter.report();
        assert_eq!(info.phase_spacing_input_samples, 1.0 / phases as f64);
        let error = (settled_start(info)..output.len())
            .map(|n| {
                let t = (n as f64 + 0.5) / f64::from(rate) - info.latency_s;
                (output[n] - sinc(PI * frequency * dt) * (2.0 * PI * frequency * t).sin()).abs()
            })
            .fold(0.0_f64, f64::max);
        errors.push(error);
    }
    assert!(errors[0] > 1e-3, "{errors:?}");
    assert!(errors[1] < errors[0] / 10.0, "{errors:?}");
    assert!(errors[2] < errors[1] / 10.0, "{errors:?}");
    assert!(errors[2] < 2e-5, "{errors:?}");
}

#[test]
fn zero_prehistory_and_original_finite_coverage_are_causal_without_tail() {
    let dt = 1.0 / 8000.0;
    let rate = 64000;
    let count = 129;
    let (output, mut filter) = render(
        dt,
        rate,
        count,
        64,
        1,
        |i, frame| frame[0] = if i == 0 { 1.0 } else { 0.0 },
        &[1],
    );
    // Source interval zero is only available at t=dt; the first eight output
    // centers precede it, although mechanics has accepted that whole interval.
    assert!(output[..8].iter().all(|x| *x == 0.0));
    assert!(output[8..].iter().any(|x| *x != 0.0));
    assert_eq!(filter.accepted_source_intervals(), 17);
    assert_eq!(filter.rendered_samples(), count as u64);
    let mut sentinel = [123.0];
    assert!(filter.render_next(&mut sentinel).is_err());
    assert!(filter.push_frame(&[0.0]).is_err());
    assert_eq!(sentinel, [123.0]);
    assert_eq!(filter.accepted_source_intervals(), 17);
}

#[test]
fn near_matching_clock_keeps_exactly_one_original_interval_per_output() {
    let rate = 48000;
    let nominal = 1.0 / f64::from(rate);
    for dt in [
        f64::from_bits(nominal.to_bits() - 1),
        nominal,
        f64::from_bits(nominal.to_bits() + 1),
    ] {
        let mut filter =
            BandlimitedResampler::new(dt, rate, 1, 1024, policy(dt, rate, 64)).unwrap();
        assert_eq!(filter.report().input_period_s.to_bits(), dt.to_bits());
        assert_eq!(filter.source_intervals(), 1024);
        for i in 0..1024 {
            assert_eq!(filter.required_source_intervals().unwrap(), i + 1);
            filter.push_frame(&[1.0]).unwrap();
            filter.render_next(&mut [0.0]).unwrap();
        }
        assert_eq!(filter.accepted_source_intervals(), 1024);
    }
}

#[test]
fn output_is_bit_identical_across_blocks_and_wide_physical_channel_arrays() {
    let dt = 1.0 / 61731.25;
    let rate = 48000;
    let channels = 257; // Native reservoir surface tiles are not WAV channels.
    let input = |i: usize, frame: &mut [f64]| {
        for (channel, x) in frame.iter_mut().enumerate() {
            *x = (channel as f64 - 128.0) * mean_tone(i, dt, 3123.0, 0.3);
        }
    };
    let (whole, a) = render(dt, rate, 257, 64, channels, input, &[257]);
    let (split, b) = render(dt, rate, 257, 64, channels, input, &[1, 7, 19, 33]);
    assert_eq!(whole, split);
    assert_eq!(a.accepted_source_intervals(), b.accepted_source_intervals());
    assert_eq!(a.report(), b.report());
    for frame in split.chunks_exact(channels) {
        assert_eq!(frame[128], 0.0);
        assert_eq!(frame[0], -frame[256]);
    }
}

#[test]
fn budgets_invalid_frames_and_incomplete_coverage_refuse_before_mutation() {
    let dt = 1.0 / 96000.0;
    let rate = 44100;
    let original = policy(dt, rate, 32);
    let reference = BandlimitedResampler::new(dt, rate, 2, 16, original).unwrap();
    let info = reference.report();
    let exact = BandlimitedObservationPolicy {
        max_preparation_terms: info.preparation_terms,
        max_output_terms: info.output_terms,
        max_history_values: info.history_values,
        ..original
    };
    assert!(BandlimitedResampler::new(dt, rate, 2, 16, exact).is_ok());
    for limited in [
        BandlimitedObservationPolicy {
            max_preparation_terms: info.preparation_terms - 1,
            ..exact
        },
        BandlimitedObservationPolicy {
            max_output_terms: info.output_terms - 1,
            ..exact
        },
        BandlimitedObservationPolicy {
            max_history_values: info.history_values - 1,
            ..exact
        },
    ] {
        assert!(BandlimitedResampler::new(dt, rate, 2, 16, limited).is_err());
    }
    for invalid in [
        BandlimitedObservationPolicy {
            passband_hz: 0.0,
            ..original
        },
        BandlimitedObservationPolicy {
            stopband_hz: 24000.0,
            ..original
        },
        BandlimitedObservationPolicy {
            passband_hz: original.stopband_hz,
            ..original
        },
        BandlimitedObservationPolicy {
            phase_count: 15,
            ..original
        },
        BandlimitedObservationPolicy {
            stopband_hz: f64::NAN,
            ..original
        },
    ] {
        assert!(BandlimitedResampler::new(dt, rate, 2, 16, invalid).is_err());
    }
    assert!(BandlimitedResampler::new(dt, rate, usize::MAX, 16, original).is_err());
    assert!(BandlimitedResampler::new(dt, rate, 2, 1 << 52, original).is_err());
    let mut filter = reference;
    let mut sentinel = [10.0, 20.0];
    assert!(filter.render_next(&mut sentinel).is_err());
    assert!(filter.push_frame(&[0.0]).is_err());
    assert!(filter.push_frame(&[0.0, f64::NAN]).is_err());
    assert_eq!(filter.accepted_source_intervals(), 0);
    assert_eq!(filter.rendered_samples(), 0);
    assert_eq!(sentinel, [10.0, 20.0]);
    let needed = filter.required_source_intervals().unwrap();
    for _ in 0..needed {
        filter.push_frame(&[2.0, 3.0]).unwrap();
    }
    assert!(filter.push_frame(&[2.0, 3.0]).is_err());
    assert!(filter.render_next(&mut [0.0]).is_err());
    assert_eq!(filter.rendered_samples(), 0);
    filter.render_next(&mut sentinel).unwrap();
    assert_eq!(filter.rendered_samples(), 1);
    assert_eq!(filter.accepted_source_intervals(), needed);
}

#[test]
fn cancellation_interrupts_cold_preparation_before_it_finishes() {
    let dt = 1.0 / 96000.0;
    let rate = 44100;
    let policy = policy(dt, rate, 1024);
    assert!(matches!(
        BandlimitedResampler::new_cancellable(dt, rate, 1, 100, policy, &mut || true),
        Err("bandlimited observation preparation cancelled")
    ));
    let mut polls = 0;
    let result = BandlimitedResampler::new_cancellable(dt, rate, 1, 100, policy, &mut || {
        polls += 1;
        polls == 4
    });
    assert!(matches!(
        result,
        Err("bandlimited observation preparation cancelled")
    ));
    assert_eq!(polls, 4);
}
