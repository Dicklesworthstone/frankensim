//! Explicit silent-key scheduling through the shared physical audio clock.
use super::*;
use super::super::{audio::AudioStream, board_geometry::SurfaceSample,
    geometry::{self, Course}, linear::BoardMode};
use fs_bem::helmholtz::Medium;

#[test]
fn silent_controls_validate_keys_values_and_do_not_count_as_excitation() {
    let rows = format!("{HEADER}\n0,silent_key_down,64,0\n0,sostenuto,0,1\n200,note_on,60,2\n");
    let p = Performance::read(&rows, &[60, 64], 2000).unwrap();
    assert_eq!(p.events[0].control, Control::SilentKeyDown { key: 64 });
    assert_eq!(p.events[1].control, Control::Sostenuto(true));
    assert_eq!(p.single_excitation_key().unwrap(), 60);
    assert_eq!(p.cursor, 0);
    assert!(Performance::from_events(p.events.clone(), &[60, 64], 2000).is_ok());
    assert!(Performance::from_events(p.events.clone(), &[60], 2000).is_err());
    let only_silent = Performance::read(&format!("{HEADER}\n0,silent_key_down,64,0\n"),
        &[64], 2000).unwrap();
    assert!(only_silent.single_excitation_key().is_err());
    for row in ["0,silent_key_down,61,0", "0,silent_key_down,60,1",
        "0,silent_key_down,60,NaN", "0,note_on,60,0"] {
        assert!(Performance::read(&format!("{HEADER}\n{row}"), &[60, 64], 2000).is_err());
    }
}

fn stream() -> AudioStream {
    let c = geometry::demonstration_scale().unwrap()[39];
    let c = Course { unison: 1, duplex_length_m: 0.0, detune_cents: 0.0, ..c };
    let board = [BoardMode { frequency_hz: 180.0, damping_ratio: 0.01,
        bridge: [0.5; 88], volume: 0.1 }];
    let piano = Instrument::new(vec![c, Course { midi: 64, ..c }], &board,
        48_000, 4, 12, true).unwrap();
    let surface = [SurfaceSample { position_m: [0.0; 3], area_m2: 0.2,
        mode_shape: vec![board[0].volume/0.2] }];
    let score = Performance::read(
        "sample,event,key,value\n0,silent_key_down,64,0\n0,sostenuto,0,1\n\
         100,note_off,64,0\n200,note_on,60,2\n800,note_off,60,0\n\
         900,sustain,0,0.5\n1000,sostenuto,0,0\n1200,sustain,0,0\n",
        &[60, 64], 1800).unwrap();
    AudioStream::new(piano, score, Some(&surface), [0.0, 0.0, 1.0],
        Medium::air(), 1.0).unwrap()
}

#[test]
fn silent_capture_release_and_half_pedal_share_the_pressure_sample_clock() {
    let mut whole = stream();
    let mut split = stream();
    let mut expected = vec![0.0; 1800];
    whole.render_block(&mut expected).unwrap();
    let mut actual = vec![0.0; 1800];
    // Silent depression, capture and release consume time without supplying
    // hammer work; the scheduled strike occurs at the next sample, exactly 200.
    split.render_block(&mut actual[..200]).unwrap();
    assert_eq!(split.sample_position(), 200);
    assert_eq!(split.instrument().accounting.input_work_j, 0.0);
    assert!(actual[..200].iter().all(|v| *v == 0.0));
    split.render_block(&mut actual[200..201]).unwrap();
    assert!(split.instrument().accounting.input_work_j > 0.0);
    split.render_block(&mut []).unwrap();
    assert_eq!(split.sample_position(), 201);
    for chunk in actual[201..].chunks_mut(127) { split.render_block(chunk).unwrap(); }
    assert_eq!(actual, expected);
    assert!(actual.iter().any(|v| v.abs() > 1e-10));
    assert_eq!(whole.sample_position(), 1800);
    assert_eq!(split.sample_position(), 1800);
    let a = whole.instrument();
    let b = split.instrument();
    assert_eq!(a.bank.q, b.bank.q);
    assert_eq!(a.bank.v, b.bank.v);
    assert_eq!(a.accounting.input_work_j, b.accounting.input_work_j);
    assert_eq!(a.accounting.dissipated_j(), b.accounting.dissipated_j());
    assert!(a.accounting.damper_loss_j > 0.0);
    let defect = a.accounting.input_work_j-a.energy_j()-a.accounting.dissipated_j();
    assert!(defect.abs() < 1e-7, "unaccounted work {defect:e} J");
}
