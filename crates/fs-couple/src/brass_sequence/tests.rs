//! Manufactured resistive paths, not a measured instrument or microphone.
use super::*;
use crate::mm_line::MmLoad;
use crate::pcm_wav::observation::DecimatedRenderer;
use fs_duct::{LossModel, Segment, TabulatedLoad, TabulatedMouth};
use fs_material::gas::GasSpec;
use fs_math::c64::C64;

struct Fixture {
    gas: GasState,
    ducts: Vec<Duct>,
    tables: Vec<TabulatedLoad>,
    line: MmLineConfig,
}
fn lip() -> LipIsland {
    LipIsland { mass_kg: 1.8e-4, stiffness_n_m: 350.0, damping_n_s_m: 0.05,
        width_m: 0.012, rest_gap_m: 5e-4, face_area_m2: 2e-5,
        provenance: "manufactured sequence regression".into() }
}
fn config() -> BrassSequenceConfig {
    BrassSequenceConfig { samples: 384, max_block_len: 384,
        initial_blowing_pressure_pa: 1500.0, initial_lip_tension: 1.0,
        initial_combo: 0 }
}
fn event(sample: u64, control: BrassControl) -> BrassEvent { BrassEvent { sample, control } }
fn valve(sample: u64, combo: usize, fade_samples: usize) -> BrassEvent {
    event(sample, BrassControl::SetValve { combo, fade_samples })
}
fn events() -> Vec<BrassEvent> {
    vec![event(0, BrassControl::SetBlowingPressure(400.0)),
        event(13, BrassControl::SetLipTension(1.1)), valve(31, 1, 11),
        event(31, BrassControl::SetBlowingPressure(600.0)),
        event(43, BrassControl::SetBlowingPressure(1000.0)), valve(127, 0, 0),
        event(128, BrassControl::SetLipTension(0.8)), valve(192, 1, 16),
        event(383, BrassControl::SetBlowingPressure(0.0))]
}
impl Fixture {
    fn new() -> Self {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
        let line = MmLineConfig { sample_rate_hz: 48000, n_modes: 2,
            extra_slices: 1, loss: LossModel::Lossless };
        let dt = 1.0 / f64::from(line.sample_rate_hz);
        let radius = 0.002;
        let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
        let ceiling = crate::driving_point::LineLadder::default().max_taps;
        let mut ducts = Vec::new();
        let mut tables = Vec::new();
        for (transit, resistance) in [(8.0, 0.7), (12.0, 1.3)] {
            let length = transit * gas.sound_speed * dt;
            ducts.push(Duct { segments: vec![Segment::Cylinder { radius, length }] });
            tables.push(TabulatedLoad::try_new(
                vec![core::f64::consts::TAU / (ceiling as f64 * dt), core::f64::consts::PI / dt],
                vec![C64::new(resistance * zc, 0.0); 2], "manufactured resistive path", 0.0,
            ).unwrap().with_mouth(TabulatedMouth { axial_position_m: length,
                radius_m: radius, gas, dc_impedance: resistance * zc }).unwrap());
        }
        Self { gas, ducts, tables, line }
    }
    fn sequence(&self, config: BrassSequenceConfig, events: Vec<BrassEvent>) -> BrassSequence {
        let loads: Vec<_> = self.tables.iter().map(MmLoad::Tabulated).collect();
        BrassSequence::new(&self.ducts, &["short", "long"], &self.gas,
            MmBankLoad::PerCombo(&loads), &self.line, lip(), 1.5e-6, config, events).unwrap()
    }
    fn manual(&self) -> BrassVoice {
        let loads: Vec<_> = self.tables.iter().map(MmLoad::Tabulated).collect();
        let mut voice = BrassVoice::new(&self.ducts, &["short", "long"], &self.gas,
            MmBankLoad::PerCombo(&loads), &self.line, lip(), 1.5e-6).unwrap();
        let c = config();
        for control in [BrassControl::SetBlowingPressure(c.initial_blowing_pressure_pa),
            BrassControl::SetLipTension(c.initial_lip_tension),
            BrassControl::SetValve { combo: c.initial_combo, fade_samples: 0 }]
        { voice.apply(control).unwrap(); }
        voice
    }
}

#[test]
fn events_and_outgoing_waves_match_original_voice_for_every_callback_partition() {
    let fixture = Fixture::new();
    let authored = events();
    let mut original = fixture.manual();
    let mut pressure = vec![0.0; 384];
    let mut waves = vec![0.0; 384];
    for sample in 0..384 {
        for event in authored.iter().filter(|event| event.sample == sample as u64) {
            original.apply(event.control).unwrap();
        }
        original.step_block_outgoing(&mut pressure[sample..sample+1], &mut waves[sample..sample+1]).unwrap();
    }
    assert!(pressure.iter().any(|p| p.abs() > 1.0));
    for block in [1, 17, 64, 384] {
        let mut sequence = fixture.sequence(config(), events());
        let mut actual = vec![0.0; 384];
        let mut outgoing = vec![0.0; 384];
        for (p, a) in actual.chunks_mut(block).zip(outgoing.chunks_mut(block)) {
            sequence.block_outgoing(p, a).unwrap();
        }
        assert_eq!(actual, pressure, "mouthpiece callback {block}");
        assert_eq!(outgoing, waves, "outgoing callback {block}");
        assert_eq!(sequence.samples_rendered(), 384);
        assert_eq!(sequence.remaining_samples(), 0);
        assert_eq!(sequence.next_event, authored.len());
        assert_eq!(sequence.voice().lift_log().len(), 3);
        for (a, b) in sequence.voice().lift_log().iter().zip(original.lift_log()) {
            assert_eq!((a.from, a.to, a.carried, a.lift), (b.from, b.to, b.carried, b.lift));
        }
        assert!(sequence.voice().receipts().iter().all(|r| r.fallback_bins == 0));
    }
    let mut plain = fixture.sequence(config(), events());
    let mut actual = vec![0.0; 384];
    plain.block(&mut actual).unwrap();
    assert_eq!(actual, pressure, "optional tracing cannot alter source physics");
}

#[test]
fn invalid_schedule_refuses_before_acoustic_realization() {
    let pressure = |sample| event(sample, BrassControl::SetBlowingPressure(1.0));
    for invalid in [vec![pressure(2), pressure(1)], vec![pressure(384)],
        vec![pressure(5), pressure(5)], vec![valve(10, 2, 0)],
        vec![valve(10, 1, 20), valve(20, 0, 0)], vec![valve(380, 1, 5)],
        vec![event(5, BrassControl::SetBlowingPressure(f64::NAN))],
        vec![event(5, BrassControl::SetLipTension(0.0))],
        vec![event(5, BrassControl::SetLipTension(f64::MAX))]]
    {
        assert!(validate_events(config(), 48000, 2, &lip(), &invalid).is_err());
    }
    // A fade may finish exactly at the next actual switch or the horizon.
    validate_events(config(), 48000, 2, &lip(),
        &[valve(10, 1, 20), valve(15, 1, 0), valve(30, 0, 354)]).unwrap();
    for changed in [BrassSequenceConfig { samples: 0, ..config() },
        BrassSequenceConfig { max_block_len: 0, ..config() },
        BrassSequenceConfig { initial_combo: 2, ..config() },
        BrassSequenceConfig { initial_blowing_pressure_pa: -1.0, ..config() }]
    { assert!(validate_events(changed, 48000, 2, &lip(), &[]).is_err()); }
    let cap = BrassSequenceConfig { samples: MAX_BRASS_EVENTS as u64 + 1, ..config() };
    let mut many: Vec<_> = (0..MAX_BRASS_EVENTS as u64).map(pressure).collect();
    validate_events(cap, 48000, 2, &lip(), &many).unwrap();
    many.push(pressure(MAX_BRASS_EVENTS as u64));
    assert!(validate_events(cap, 48000, 2, &lip(), &many).is_err());
    let invalid_duct = Duct { segments: Vec::new() };
    let gas = Fixture::new().gas;
    let line = MmLineConfig { sample_rate_hz: 48000, n_modes: 2,
        extra_slices: 1, loss: LossModel::Lossless };
    // This must be the schedule refusal, not the deliberately invalid bore.
    assert!(matches!(BrassSequence::new(&[invalid_duct], &["invalid"], &gas,
        MmLoad::Analytic(fs_duct::Termination::Closed), &line, lip(), 1.5e-6,
        config(), vec![pressure(384)]), Err(RenderError::Control { .. })));
}

#[test]
fn callback_refusals_leave_the_original_sequence_untouched() {
    let fixture = Fixture::new();
    let mut sequence = fixture.sequence(config(), events());
    assert!(sequence.block(&mut []).is_err());
    let mut too_long = vec![7.0; 385];
    assert!(sequence.block(&mut too_long).is_err());
    assert!(too_long.iter().all(|x| *x == 7.0));
    assert!(sequence.block_outgoing(&mut [7.0; 2], &mut [9.0; 1]).is_err());
    assert!(sequence.validate_sample_rate(24000).is_err());
    assert_eq!(sequence.samples_rendered(), 0);
    assert_eq!(sequence.next_event, 0);
    let mut actual = vec![0.0; 384];
    sequence.block(&mut actual).unwrap();
    let mut expected = vec![0.0; 384];
    fixture.sequence(config(), events()).block(&mut expected).unwrap();
    assert_eq!(actual, expected);
    assert!(sequence.block(&mut [0.0]).is_err());
    assert_eq!(sequence.samples_rendered(), 384);
}

#[test]
fn a_physical_failure_poisons_continuation_after_any_advanced_prefix() {
    let fixture = Fixture::new();
    // A finite control may still overflow the physical nonlinear solve. The
    // wrapper must not retry a partly advanced lip/line trajectory as healthy.
    let mut sequence = fixture.sequence(config(),
        vec![event(8, BrassControl::SetBlowingPressure(f64::MAX))]);
    assert!(sequence.block(&mut [0.0; 16]).is_err());
    assert!(matches!(sequence.block(&mut [0.0]), Err(RenderError::Poisoned)));
    assert!(matches!(sequence.validate_sample_count(1), Err(RenderError::Poisoned)));
    assert!(matches!(sequence.validate_sample_rate(48000), Err(RenderError::Poisoned)));
}

#[test]
fn existing_decimation_changes_only_the_observation_clock() {
    let fixture = Fixture::new();
    let render = |block| {
        let source = fixture.sequence(config(), events());
        let mut observer = DecimatedRenderer::new(source, 48000, 24000, 17).unwrap();
        let mut pressure = vec![0.0; 192];
        for chunk in pressure.chunks_mut(block) { observer.block(chunk).unwrap(); }
        assert_eq!(observer.samples_rendered(), 192);
        pressure
    };
    assert_eq!(render(1), render(17));
}
