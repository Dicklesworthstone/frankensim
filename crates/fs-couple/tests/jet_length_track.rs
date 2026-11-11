//! The jet-length track (`GestureTarget::JetLength`, MB 23643) on the committed D4 flute card at the
//! coarse forced-jet row with the Bickley base and the Fletcher swing limit, cut to 1024 output
//! samples of its opening ramp. A track held at the flue row's cut-up renders the card's bits; a
//! moving track moves the jet, the swing limit's lip and the sidecar's rate; a wall-cut window (the
//! D5 whistle's), the island's chain and a zero length refuse it by name.

use fs_couple::flue_loop::JET_LENGTH_NEEDS_LIP_JET;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;
use fs_scenario::gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue};

const FLUTE: &str =
    include_str!("../../../data/instruments/descriptions/flute/flute-d4.performance");
const WHISTLE: &str =
    include_str!("../../../data/instruments/descriptions/tin-whistle/tin-whistle-d5.performance");
const JET: &str = "jet 90 1e-6 3 0\n";
const ROWS: &str =
    "jet_response forced-bickley 0.25 40\njet_spreading bickley\njet_swing_limit fletcher-1976\n";
const OUTPUT: &str = "output 120000 5 1000000\n";
const SHORT: &str = "output 1024 5 1000000\n";

/// `text` cut to 1024 output samples, `rows` after its jet row, and its gesture's opening pressure
/// ramp with `lip` beside it.
fn card(text: &str, rows: &str, lip: Option<GestureTrack>) -> String {
    assert_eq!(text.matches(JET).count(), 1, "the jet row changed; update this test");
    assert_eq!(text.matches(OUTPUT).count(), 1, "the output row changed; update this test");
    let (head, gesture) = text.split_once("\ngesture\n").expect("a gesture section");
    let schedule = GestureSchedule::from_canonical_bytes(gesture.as_bytes()).expect("canonical");
    let mut blow = schedule.tracks()[0].clone();
    assert_eq!(blow.target, GestureTarget::BlowingPressure);
    blow.events.truncate(1);
    let tracks = core::iter::once(blow).chain(lip).collect();
    let schedule = GestureSchedule::try_new(schedule.control_rate_hz, tracks).unwrap();
    let head = head.replace(OUTPUT, SHORT).replace(JET, &format!("{JET}{rows}"));
    format!("{head}\ngesture\n{}", String::from_utf8(schedule.to_canonical_bytes()).unwrap())
}

fn parse(text: &str) -> Result<FlueFilePerformance, String> {
    FlueFilePerformance::from_bytes(text.as_bytes(), 256, &CancelGate::new_clock_free())
        .map_err(|e| format!("{e:?}"))
}

fn lip(initial: f64, events: Vec<GestureEvent>) -> GestureTrack {
    GestureTrack { id: "lip".into(), target: GestureTarget::JetLength,
        initial: GestureValue::LengthM(initial), events }
}

fn render(file: FlueFilePerformance) -> Vec<u64> {
    let mut pressure = vec![0.0; file.info().samples as usize];
    let mut renderer = file.into_renderer();
    for part in pressure.chunks_mut(256) {
        renderer.block(part).unwrap();
    }
    assert!(pressure.iter().all(|p| p.is_finite()));
    pressure.iter().map(|p| p.to_bits()).collect()
}

#[test]
fn a_track_held_at_the_flue_rows_cut_up_renders_the_cards_bits() {
    let plain = parse(&card(FLUTE, ROWS, None)).expect("parses");
    let cut_up = plain.geometry().cut_up_m;
    assert!(!plain.renderer().voice().jet_length_admitted());
    assert_eq!(plain.renderer().voice().forced_jet().unwrap().jet_length_span_m(), None);
    let held = parse(&card(FLUTE, ROWS, Some(lip(cut_up, vec![])))).expect("parses");
    assert!(held.renderer().voice().jet_length_admitted());
    assert_eq!(held.renderer().voice().forced_jet().unwrap().jet_length_span_m(), Some((cut_up, cut_up)));
    assert_eq!(held.renderer().jet_length_rate_m_s(), Some(0.0));
    let expected = render(plain);
    assert!(expected.iter().any(|p| f64::from_bits(*p) != 0.0));
    assert_eq!(render(held), expected);
}

#[test]
fn a_moving_track_moves_the_jet_and_the_swing_limits_lip() {
    let plain = parse(&card(FLUTE, ROWS, None)).expect("parses");
    let cut_up = plain.geometry().cut_up_m;
    let short = cut_up - 5e-4;
    // From the cut-up to 0.5 mm shorter over 10 ms from 5 ms on, then held.
    let moving = parse(&card(FLUTE, ROWS, Some(lip(cut_up, vec![GestureEvent {
        time_s: 5e-3, transition_s: 1e-2, value: GestureValue::LengthM(short) }])))).expect("parses");
    {
        let voice = moving.renderer().voice();
        let kernel = voice.forced_jet().unwrap();
        assert_eq!(kernel.jet_length_span_m(), Some((short, cut_up)));
        assert_eq!(kernel.jet_length_m(), cut_up, "read at the cut-up, inside the span");
        let limit = voice.jet_swing_limit().expect("limit");
        assert_eq!(limit.cut_up_m(), cut_up);
        let rate = moving.renderer().jet_length_rate_m_s().unwrap();
        assert!((rate - 5e-4 / 1e-2).abs() < 1e-9, "{rate:e}");
    }
    let mut pressure = vec![0.0; moving.info().samples as usize];
    let mut renderer = moving.into_renderer();
    for part in pressure.chunks_mut(256) {
        renderer.block(part).unwrap();
    }
    let voice = renderer.voice();
    let kernel = voice.forced_jet().unwrap();
    assert_eq!(kernel.jet_length_m(), short);
    let (beta, _) = kernel.station_spread().unwrap();
    let limit = voice.jet_swing_limit().unwrap();
    assert_eq!(limit.cut_up_m(), short);
    assert_eq!(limit.spread_width_m(), beta * kernel.half_width_m());
    assert!(pressure.iter().all(|p| p.is_finite()));
    let moved: Vec<u64> = pressure.iter().map(|p| p.to_bits()).collect();
    assert_ne!(moved, render(plain), "the jet length must reach the microphone");
}

#[test]
fn a_length_track_is_refused_by_name_off_the_forced_jet_on_an_embouchure() {
    for (text, rows, length, refusal) in [
        (WHISTLE, "", 0.0053, JET_LENGTH_NEEDS_LIP_JET),
        (FLUTE, "", 0.0075, "a jet-length track is refused by the island's convection chain"),
        (FLUTE, ROWS, 0.0, "jet length must be finite and positive"),
    ] {
        match parse(&card(text, rows, Some(lip(length, vec![])))) {
            Ok(_) => panic!("accepted a jet-length track: {refusal}"),
            Err(got) => assert!(got.contains(refusal), "refused as `{got}`, expected `{refusal}`"),
        }
    }
}
