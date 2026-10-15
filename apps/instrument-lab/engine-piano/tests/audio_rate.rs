//! G0/G3 real-owner integration: no surrogate piano, new material or fitted sound.
use fs_live_piano::{Config, Description, LivePiano, PianoAudio};
use fs_live_piano::refmods::{board_geometry::SurfaceSample, geometry::Course,
    linear::BoardMode, steinway_scale};
use fs_couple::pcm_wav::decimate::Decimator;

fn fixture() -> (Description, Config) {
    let course = steinway_scale::courses().unwrap().into_iter().find(|c| c.midi == 69).unwrap();
    let mut bridge = [0.0; 88]; bridge[48] = 0.4;
    let board = vec![BoardMode { frequency_hz: 310.0, damping_ratio: 0.01,
        bridge, volume: 0.02 }];
    let desc = Description {
        courses: vec![Course { duplex_length_m: 0.0, ..course }], board,
        damping: None, surface: vec![SurfaceSample {
            position_m: [0.0, 0.0, 0.0], area_m2: 0.02, mode_shape: vec![1.0],
        }], board_provenance: "manufactured test board, source string and felt".into(),
        board_mass_kg: 1.0, board_area_m2: 0.02, exact_low_modes: 1,
    };
    let mut cfg = Config::live(48_000);
    cfg.max_partials = 16;
    cfg.course_mode_budget = 48;
    cfg.max_dynamic_modes = 64;
    cfg.partial_band_hz = 6_000.0;
    cfg.duplex = false;
    (desc, cfg)
}

#[test]
fn factor_one_retains_the_original_pressure_and_float_bits() {
    let (desc, cfg) = fixture();
    let mut old = LivePiano::new(&desc, cfg.clone()).unwrap();
    let mut new = PianoAudio::new(&desc, cfg, 1).unwrap();
    old.note_on(69, 1.0).unwrap(); new.note_on(69, 1.0).unwrap();
    let mut a = [0.0; 2048]; let mut b = [0.0; 2048];
    old.render(&mut a); new.render(&mut b).unwrap();
    assert!(a.iter().any(|x| *x != 0.0));
    assert!(a.iter().zip(b).all(|(x,y)| x.to_bits() == y.to_bits()));
    assert_eq!(old.counters.samples, new.frames());
    assert_eq!(old.counters.peak_pa.to_bits(), new.peak_pa().to_bits());
    assert_eq!(format!("{:?}", old.close_ledger()), format!("{:?}", new.close_ledger()));
    assert_eq!(new.decimation_delay_frames(), 0.0);
}

#[test]
fn every_high_rate_physical_pressure_sample_reaches_the_shared_filter() {
    let (desc, cfg) = fixture();
    for ratio in [2, 4] {
        let mut high_cfg = cfg.clone(); high_cfg.sample_rate *= ratio as u32;
        let mut direct = LivePiano::new(&desc, high_cfg).unwrap();
        let gas = direct.acoustic_gas();
        let mut output = PianoAudio::new_with_gas(&desc, cfg.clone(), ratio, gas).unwrap();
        let mut filter = Decimator::new(ratio, 1).unwrap();
        let mut packet = vec![0.0; ratio];
        let mut endpoint_diff = 0.0_f64;
        direct.note_on(69, 1.2).unwrap(); output.note_on(69, 1.2).unwrap();
        for frame in 0..2048 {
            if frame == 700 {
                direct.set_sustain(0.7).unwrap(); output.set_sustain(0.7).unwrap();
                direct.note_off(69); output.note_off(69);
            }
            for x in &mut packet { *x = direct.step(); }
            assert!(direct.take_error().is_none());
            let expected = filter.preview(&packet).unwrap()[0]; filter.commit();
            let actual = output.step_pressure().unwrap();
            assert_eq!(actual.to_bits(), expected.to_bits(), "ratio {ratio}, frame {frame}");
            endpoint_diff = endpoint_diff.max((actual-packet[ratio-1]).abs());
        }
        assert!(endpoint_diff > 1e-10, "endpoint-only downsampling was not distinguished");
        assert_eq!(output.frames(), 2048);
        assert_eq!(output.piano().counters.samples, 2048 * ratio as u64);
        assert_eq!(output.mechanical_rate(), 48_000 * ratio as u32);
        assert_eq!(output.piano().acoustic_gas().sound_speed.to_bits(), gas.sound_speed.to_bits());
        assert_eq!(direct.partials_hz(69), output.piano().partials_hz(69));
        assert_eq!(format!("{:?}", direct.close_ledger()), format!("{:?}", output.close_ledger()));
        let (travel, delay) = direct.microphone_delay_samples();
        assert_eq!(output.microphone_delay_frames(),
            (travel as f64 / ratio as f64, delay / ratio as f64 + filter.delay_output_frames()));
    }
}

#[test]
fn callbacks_and_mid_decay_clones_preserve_both_histories() {
    let (desc, cfg) = fixture();
    let mut start = PianoAudio::new(&desc, cfg, 4).unwrap();
    start.note_on(69, 0.9).unwrap();
    start.render(&mut [0.0; 777]).unwrap(); // clone a nonempty FIR and live felt history
    let mut whole = start.clone();
    let mut expected = vec![0.0; 2048]; whole.render(&mut expected).unwrap();
    for block in [1, 37, 128, 1024] {
        let mut split = start.clone();
        let mut actual = vec![0.0; expected.len()];
        for chunk in actual.chunks_mut(block) { split.render(chunk).unwrap(); }
        assert_eq!(actual, expected, "callback {block}");
        assert_eq!(split.frames(), whole.frames());
        assert_eq!(split.piano().counters.samples, whole.piano().counters.samples);
        assert_eq!(format!("{:?}", split.close_ledger()), format!("{:?}", whole.close_ledger()));
    }
}

#[test]
fn rate_and_output_band_admission_precede_physics_and_do_not_retune() {
    let (desc, cfg) = fixture();
    for ratio in [0, 3, 8, usize::MAX] {
        assert!(PianoAudio::new(&desc, cfg.clone(), ratio).is_err());
    }
    for rate in [0, 7999, 48_001, 192_000, u32::MAX] {
        let mut bad = cfg.clone(); bad.sample_rate = rate;
        assert!(PianoAudio::new(&desc, bad, 4).is_err());
    }
    let mut bad = cfg.clone(); bad.partial_band_hz = 23_000.0;
    assert!(PianoAudio::new(&desc, bad, 4).is_err()); // NOT admitted by 192 kHz Nyquist
    let mut audio = PianoAudio::new(&desc, cfg.clone(), 2).unwrap();
    let original = audio.piano().partials_hz(69);
    for band in [0.0, f64::NAN, f64::INFINITY, 23_000.0] {
        assert!(audio.set_quality(48, band, 64).is_err());
    }
    assert_eq!(audio.piano().partials_hz(69), original);
    assert_eq!(audio.piano().config().partial_band_hz, cfg.partial_band_hz);
    assert!(!audio.is_poisoned()); assert_eq!(audio.frames(), 0);
}

#[test]
fn partial_source_failure_never_restarts_a_missing_output_packet() {
    let (desc, mut cfg) = fixture(); cfg.max_dynamic_modes = 1;
    let mut audio = PianoAudio::new(&desc, cfg, 4).unwrap();
    audio.note_on(69, 1.0).unwrap(); // original engine queues admission until step
    assert!(audio.step_pressure().is_err());
    assert!(audio.is_poisoned()); assert_eq!(audio.frames(), 0);
    let accepted = audio.piano().counters.samples;
    assert_eq!(accepted, 1); // source advanced; no complete audio packet accepted
    assert!(audio.take_error().is_some()); assert!(audio.take_error().is_none());
    let mut out = [123.0; 20];
    assert!(audio.render(&mut out).is_err()); assert_eq!(out, [0.0; 20]);
    assert_eq!(audio.piano().counters.samples, accepted);
    assert!(audio.note_on(69, 1.0).is_err());
    assert!(audio.clone().step_pressure().is_err());
}

#[test]
fn browser_description_reports_host_clock_and_complete_numerical_latency() {
    let (desc, cfg) = fixture();
    for ratio in [1, 2, 4] {
        let audio = PianoAudio::new(&desc, cfg.clone(), ratio).unwrap();
        let info = fs_live_piano::describe_audio(&audio);
        assert!(info.contains("\"sampleRate\":48000,"));
        assert!(info.contains(&format!("\"mechanicalRateHz\":{},\"substeps\":{ratio},", 48000 * ratio)));
        let (_, delay) = audio.microphone_delay_frames();
        assert!(info.contains(&format!("\"microphoneFilterDelaySamples\":{delay},")));
        assert!(info.contains("\"partialBandHz\":6000,"));
        if ratio == 1 { assert_eq!(info, fs_live_piano::describe(audio.piano())); }
    }
}
