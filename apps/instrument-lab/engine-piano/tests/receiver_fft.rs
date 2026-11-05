//! The restored receiver changes only observer evaluation, not mechanics or latency.
use fs_live_piano::{Config, Counters, Description, LivePiano};
use std::time::{Duration, Instant};

fn discrete_counts(counters: Counters) -> [u64; 9] {
    // Exhaustive destructuring makes a newly added counter require a decision
    // here. Peak pressure is the receiver's floating-point result, not a count.
    let Counters {
        samples,
        note_ons,
        refused_note_ons,
        rebuilds,
        retirements,
        forced_retirements,
        contact_failures,
        clipped_samples,
        peak_pa: _,
        stretch_out_of_domain,
    } = counters;
    [
        samples,
        note_ons,
        refused_note_ons,
        rebuilds,
        retirements,
        forced_retirements,
        contact_failures,
        clipped_samples,
        stretch_out_of_domain,
    ]
}

fn event(piano: &mut LivePiano, block: usize) {
    match block {
        0 => {
            for key in [36, 55, 60, 67] {
                piano.note_on(key, 2.5).unwrap();
            }
        }
        16 => piano.set_sustain(0.8).unwrap(),
        32 => {
            for key in [36, 55, 60, 67] {
                piano.note_off(key);
            }
        }
        48 => {
            piano.note_on(79, 1.7).unwrap();
            piano.note_on(43, 2.2).unwrap();
        }
        64 => piano.set_sustain(0.0).unwrap(),
        80 => {
            piano.note_off(79);
            piano.note_off(43);
        }
        _ => {}
    }
}

#[test]
fn full_piano_fft_preserves_pressure_mechanics_and_clone() {
    let description = Description::embedded().unwrap();
    for (bank_hz, browser_tone) in [(0.0, false), (600.0, false), (600.0, true)] {
        let mut config = Config::live(48_000);
        config.silent_bank_hz = bank_hz;
        // Match the actual browser's string-loss and measured-hammer choices
        // as well as the native defaults; preserve both existing comparisons.
        config.rt0425_string_damping = browser_tone;
        config.felt_pack_hammers = browser_tone;
        config.partitioned_receiver = false;
        let mut direct = LivePiano::new(&description, config.clone()).unwrap();
        config.partitioned_receiver = true;
        let mut fast = LivePiano::new(&description, config).unwrap();
        assert!(!direct.uses_partitioned_receiver());
        assert!(fast.uses_partitioned_receiver());
        let (mut signal, mut error, mut peak_error) = (0.0_f64, 0.0_f64, 0.0_f64);
        let mut timings = [Duration::ZERO; 2];
        let mut worst = [Duration::ZERO; 2];
        let mut audio = [[0.0_f32; 128]; 2];
        // Interleave both engines in the same invocation and alternate which
        // renders first. Timings are measurements, not a CI performance gate.
        for block in 0..256 {
            event(&mut direct, block);
            event(&mut fast, block);
            for lane in if block % 2 == 0 { [0, 1] } else { [1, 0] } {
                let start = Instant::now();
                if lane == 0 {
                    direct.render(&mut audio[lane]);
                } else {
                    fast.render(&mut audio[lane]);
                }
                let elapsed = start.elapsed();
                timings[lane] += elapsed;
                worst[lane] = worst[lane].max(elapsed);
            }
            for (&a, &b) in audio[0].iter().zip(&audio[1]) {
                assert!(a.is_finite() && b.is_finite());
                let delta = f64::from(a) - f64::from(b);
                signal += f64::from(a).powi(2);
                error += delta.powi(2);
                peak_error = peak_error.max(delta.abs());
            }
            assert_eq!(
                format!("{:?}", direct.close_ledger()),
                format!("{:?}", fast.close_ledger())
            );
            assert_eq!(
                discrete_counts(direct.counters),
                discrete_counts(fast.counters)
            );
            // FFT arithmetic may change the last bits of physical pressure.
            // This is a much tighter bound than the separately checked f32 PCM.
            let reference_peak = direct.counters.peak_pa;
            let fft_peak = fast.counters.peak_pa;
            assert!(reference_peak.is_finite() && fft_peak.is_finite());
            assert!((reference_peak - fft_peak).abs() < 1e-10 * (1.0 + reference_peak));
            assert_eq!(direct.dynamic_modes(), fast.dynamic_modes());
            assert_eq!(direct.take_error(), None);
            assert_eq!(fast.take_error(), None);
            if block == 97 {
                direct = direct.clone();
                fast = fast.clone();
            }
        }
        assert!(signal > 1e-6, "silent comparison");
        assert!(peak_error < 1e-7, "peak PCM error {peak_error}");
        assert!(error / signal < 1e-12, "residual must stay below -120 dB");
        eprintln!(
            "receiver bank_hz={bank_hz} browser_tone={browser_tone} samples=32768 direct_ms={:.3} fft_ms={:.3} worst_direct_ms={:.3} worst_fft_ms={:.3} residual_db={} peak_error={peak_error}",
            timings[0].as_secs_f64() * 1000.0,
            timings[1].as_secs_f64() * 1000.0,
            worst[0].as_secs_f64() * 1000.0,
            worst[1].as_secs_f64() * 1000.0,
            10.0 * (error / signal).log10()
        );
    }
}

#[test]
fn active_coefficients_preserve_pcm_rebuilds_and_rollback() {
    let description = Description::embedded().unwrap();
    for (bank_hz, browser_tone) in [(0.0, false), (600.0, false), (600.0, true)] {
        let mut config = Config::live(48_000);
        config.silent_bank_hz = bank_hz;
        config.rt0425_string_damping = browser_tone;
        config.felt_pack_hammers = browser_tone;
        config.partitioned_receiver = true;
        config.active_coefficients = false;
        let mut original = LivePiano::new(&description, config.clone()).unwrap();
        config.active_coefficients = true;
        let mut packed = LivePiano::new(&description, config).unwrap();
        let mut timings = [Duration::ZERO; 2];
        let mut worst = [Duration::ZERO; 2];
        let mut audio = [[0.0_f32; 128]; 2];
        let mut heard = false;
        for block in 0..256 {
            for piano in [&mut original, &mut packed] {
                event(piano, block);
                match block {
                    97 => *piano = piano.clone(),
                    100 => {
                        piano.fail_next_rebuild();
                        piano.note_on(84, 2.2).unwrap();
                    }
                    101 => piano.note_on(84, 2.2).unwrap(),
                    128 => piano.set_quality(48, 3000.0, 2400).unwrap(),
                    132 => piano.note_on(42, 2.2).unwrap(),
                    144 => piano.set_quality(200, 6000.0, 4800).unwrap(),
                    148 => piano.note_on(72, 2.2).unwrap(),
                    160 => {
                        for key in [42, 72, 84] {
                            piano.note_off(key);
                        }
                    }
                    _ => {}
                }
            }
            for lane in if block % 2 == 0 { [0, 1] } else { [1, 0] } {
                let start = Instant::now();
                if lane == 0 {
                    original.render(&mut audio[lane]);
                } else {
                    packed.render(&mut audio[lane]);
                }
                let elapsed = start.elapsed();
                timings[lane] += elapsed;
                worst[lane] = worst[lane].max(elapsed);
            }
            for (&a, &b) in audio[0].iter().zip(&audio[1]) {
                assert!(a.is_finite() && b.is_finite());
                assert_eq!(a.to_bits(), b.to_bits(), "bank={bank_hz}, block={block}");
                heard |= a.abs() > 1e-4;
            }
            assert_eq!(
                format!("{:?}", original.close_ledger()),
                format!("{:?}", packed.close_ledger())
            );
            assert_eq!(
                format!("{:?}", original.counters),
                format!("{:?}", packed.counters)
            );
            assert_eq!(original.dynamic_modes(), packed.dynamic_modes());
            let error = original.take_error();
            assert_eq!(error, packed.take_error());
            assert_eq!(
                error.is_some(),
                block == 100,
                "failed rebuild must be exercised"
            );
        }
        assert!(heard, "silent comparison");
        eprintln!(
            "active coefficients bank_hz={bank_hz} browser_tone={browser_tone} exact_samples=32768 original_ms={:.3} packed_ms={:.3} worst_original_ms={:.3} worst_packed_ms={:.3}",
            timings[0].as_secs_f64() * 1000.0,
            timings[1].as_secs_f64() * 1000.0,
            worst[0].as_secs_f64() * 1000.0,
            worst[1].as_secs_f64() * 1000.0
        );
    }
}
