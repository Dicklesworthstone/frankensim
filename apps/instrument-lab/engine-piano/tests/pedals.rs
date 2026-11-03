//! The live una corda and sostenuto pedals. Una corda against the offline reference engine
//! (same source files, one course on the bridge, as examples/validate.rs): the strikes agree
//! to roundoff. Sostenuto against its definition: a key the rod caught sounds exactly as a key
//! held down, and a key pressed after the rod rose is damped exactly as without it.
use fs_live_piano::refmods::{engine, geometry, steinway_scale};
use fs_live_piano::{Config, Description, LivePiano};

const RATE: u32 = 48_000;
const PARTIALS: usize = 256;

/// Relative bridge-velocity difference (max over the run) and the bridge velocity's peak.
fn una_corda_run(midi: u8, una_corda: bool) -> (f64, f64, Vec<f64>) {
    let scale = geometry::read_scale(fs_live_piano::kernel::SCALE_CSV).unwrap();
    let course = scale[usize::from(midi - 21)];
    let desc = Description::with_courses(vec![course], 800.0).unwrap();
    let material = steinway_scale::hammer_material(&course).unwrap();
    // 256 partials: the table reaches the partial band, so no cap binds. A binding cap ends the
    // live duplex segments where the speaking strings end, which the reference does not do.
    let mut reference = engine::Instrument::new_with_string_damping(vec![course], &desc.board, RATE, 1, PARTIALS, true,
        vec![material], Some(engine::ShankGeometry::published()), None, None, false).unwrap();
    let mut cfg = Config::live(RATE);
    cfg.board_band_hz = 800.0;
    cfg.max_partials = PARTIALS;
    cfg.course_mode_budget = usize::MAX / 4;
    cfg.static_attachment = false;
    cfg.damper_points = 1;
    let mut live = LivePiano::new(&desc, cfg).unwrap();
    reference.set_una_corda(una_corda);
    live.set_una_corda(una_corda);
    reference.note_on(midi, 2.5).unwrap();
    live.note_on(midi, 2.5).unwrap();
    live.apply_pending_notes();
    let row: Vec<f64> = reference.bank.strings[0].bridge.clone();
    let n = reference.bank.modes.len();
    let (mut peak, mut worst, mut trace) = (0.0_f64, 0.0_f64, Vec::new());
    for _ in 0..RATE / 5 {
        reference.step().unwrap();
        live.step();
        let vr: f64 = row.iter().zip(&reference.bank.v[n..]).map(|(g, v)| g * v).sum();
        let vl = live.bridge_velocity(midi).unwrap();
        peak = peak.max(vr.abs());
        worst = worst.max((vr - vl).abs());
        trace.push(vl);
    }
    assert!(live.take_error().is_none());
    (worst / peak, peak, trace)
}

#[test]
fn una_corda_strikes_match_the_reference_engine() {
    // A trichord keeps two strings, a bichord one, a single string its own.
    for (midi, changes) in [(57u8, true), (35, true), (25, false)] {
        let (rel, _, shifted) = una_corda_run(midi, true);
        assert!(rel < 1.0e-9, "key {midi}: live vs reference with una corda, {rel:e}");
        let (_, peak, normal) = una_corda_run(midi, false);
        let moved = shifted.iter().zip(&normal).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max) / peak;
        assert_eq!(moved > 1.0e-3, changes, "key {midi}: una corda moved the bridge by {moved:e}");
    }
}

fn piano() -> LivePiano {
    LivePiano::new(&Description::embedded().unwrap(), Config::live(RATE)).unwrap()
}
fn play(p: &mut LivePiano, seconds: f64) -> Vec<u32> {
    let mut out = vec![0.0_f32; (seconds * f64::from(RATE)) as usize];
    p.render(&mut out);
    out.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn a_key_caught_by_the_sostenuto_sounds_as_a_held_key() {
    let mut held = piano();
    held.note_on(57, 2.0).unwrap();
    let mut a = play(&mut held, 0.3);
    a.extend(play(&mut held, 1.0));

    let mut caught = piano();
    caught.note_on(57, 2.0).unwrap();
    let mut b = play(&mut caught, 0.3);
    caught.set_sostenuto(true);
    caught.note_off(57);
    b.extend(play(&mut caught, 1.0));
    assert!(a == b, "a latched key must ring exactly as a held one");

    let mut released = piano();
    released.note_on(57, 2.0).unwrap();
    let mut c = play(&mut released, 0.3);
    released.note_off(57);
    c.extend(play(&mut released, 1.0));
    let tail = |x: &[u32]| x[x.len() - 4800..].iter().map(|v| f64::from(f32::from_bits(*v)).powi(2)).sum::<f64>();
    assert!(tail(&a) > 100.0 * tail(&c), "the dampers must fall on a released key without the sostenuto");

    // The rod rising lets the dampers fall on the caught key.
    caught.set_sostenuto(false);
    assert!(play(&mut caught, 0.5) != play(&mut held, 0.5));
}

#[test]
fn the_sostenuto_catches_only_what_is_up_when_it_is_pressed() {
    let released = |prepare: &dyn Fn(&mut LivePiano)| {
        let mut p = piano();
        prepare(&mut p);
        p.note_on(57, 2.0).unwrap();
        let mut x = play(&mut p, 0.3);
        p.note_off(57);
        x.extend(play(&mut p, 1.0));
        x
    };
    let plain = released(&|_| {});
    // Pressed before the key went down: the key is not caught.
    assert!(released(&|p| p.set_sostenuto(true)) == plain);
    // Pressed while the sustain pedal is fully down: every damper is caught.
    let caught_by_sustain = released(&|p| {
        p.set_sustain(1.0).unwrap();
        p.set_sostenuto(true);
        p.set_sustain(0.0).unwrap();
    });
    let mut held = piano();
    held.note_on(57, 2.0).unwrap();
    let held = play(&mut held, 1.3);
    assert!(caught_by_sustain == held);
    // Half a sustain pedal raises no damper fully: nothing is caught.
    assert!(released(&|p| {
        p.set_sustain(0.5).unwrap();
        p.set_sostenuto(true);
        p.set_sustain(0.0).unwrap();
    }) == plain);
    // Every key up lifts the sostenuto too.
    let mut p = piano();
    p.set_sostenuto(true);
    p.set_una_corda(true);
    p.all_notes_off();
    assert!(!p.sostenuto() && !p.una_corda());
}

fn damper_description(keys: &[u8]) -> Description {
    let mut desc = Description::embedded().unwrap().truncated(4).unwrap();
    desc.courses.retain(|c| keys.contains(&c.midi));
    desc
}

fn damper_config(text: String) -> Config {
    let mut cfg = Config::live(RATE);
    cfg.max_partials = 16;
    cfg.course_mode_budget = 48;
    cfg.duplex = false;
    cfg.retire_ringing_j = 0.0;
    cfg.retire_damped_j = 0.0;
    cfg.damper_specification = Some(text);
    cfg
}

#[test]
fn supplied_free_and_pad_rows_control_release_and_sostenuto_across_the_old_break() {
    let desc = damper_description(&[57, 89]);
    let length = desc.courses[1].length_m;
    let cfg = damper_config(format!("frankensim-piano-dampers-v1\nfree,57\npad,89,{},{},0.4\n",
        0.31 * length, 0.39 * length));
    for midi in [57, 89] {
        let mut held = LivePiano::new(&desc, cfg.clone()).unwrap();
        held.note_on(midi, 2.0).unwrap();
        play(&mut held, 0.03);
        let mut released = held.clone();
        let mut caught = held.clone();
        let mut sustained = held.clone();
        released.note_off(midi);
        caught.set_sostenuto(true);
        caught.note_off(midi);
        sustained.set_sustain(1.0).unwrap();
        sustained.note_off(midi);
        let ringing = play(&mut held, 0.06);
        let stopped = play(&mut released, 0.06);
        assert!(play(&mut caught, 0.06) == ringing);
        assert!(play(&mut sustained, 0.06) == ringing);
        if midi == 57 {
            assert!(stopped == ringing, "an explicit free key has no release drag");
            assert_eq!(released.ledger.damper_loss_j, 0.0);
        } else {
            assert!(stopped != ringing, "the supplied pad above MIDI 88 must act");
            assert!(released.ledger.damper_loss_j > 0.0);
            caught.set_sostenuto(false);
            assert!(play(&mut caught, 0.03) != play(&mut held, 0.03));
            assert!(caught.ledger.damper_loss_j > 0.0);
        }
        for piano in [&mut held, &mut released, &mut caught, &mut sustained] {
            let ledger = piano.close_ledger();
            assert!(ledger.residual_j.abs() < 1e-8 * ledger.input_j);
            assert!(piano.take_error().is_none());
        }
    }
}

#[test]
fn authored_drag_uses_the_existing_squared_pedal_factor_and_pad_geometry() {
    let desc = damper_description(&[57]);
    let length = desc.courses[0].length_m;
    let make = |start: f64, end: f64, drag: f64| {
        let text = format!("frankensim-piano-dampers-v1\npad,57,{},{},{drag}\n", start * length, end * length);
        LivePiano::new(&desc, damper_config(text)).unwrap()
    };
    let mut half_pedal = make(0.31, 0.39, 0.4);
    let mut quarter_drag = make(0.31, 0.39, 0.1);
    let mut moved = make(0.21, 0.29, 0.1);
    for piano in [&mut half_pedal, &mut quarter_drag, &mut moved] { piano.note_on(57, 2.0).unwrap(); }
    let initial = play(&mut half_pedal, 0.03);
    assert!(play(&mut quarter_drag, 0.03) == initial);
    assert!(play(&mut moved, 0.03) == initial, "raised pads do not affect the strike");
    half_pedal.set_sustain(0.5).unwrap();
    for piano in [&mut half_pedal, &mut quarter_drag, &mut moved] { piano.note_off(57); }
    let a = play(&mut half_pedal, 0.04);
    assert!(play(&mut quarter_drag, 0.04) == a, "0.4 * (1 - 0.5)^2 is the same physical drag as 0.1");
    assert!(play(&mut moved, 0.04) != a, "actual pad stations must reach the modal projection");
    assert!(half_pedal.ledger.damper_loss_j > 0.0 && moved.ledger.damper_loss_j > 0.0);
    let mut copy = moved.clone();
    assert!(play(&mut copy, 0.02) == play(&mut moved, 0.02), "timing copies retain the supplied pads");
    for piano in [&mut half_pedal, &mut quarter_drag, &mut moved, &mut copy] {
        assert!(piano.take_error().is_none());
    }
}

#[test]
fn the_complete_damper_source_is_admitted_before_any_note() {
    let desc = damper_description(&[57, 89]);
    let header = "frankensim-piano-dampers-v1";
    let length = desc.courses[0].length_m;
    let start = 0.3 * length;
    let end = f64::from_bits(start.to_bits() + 1);
    for text in [String::new(), format!("{header}\nfree,57\n"),
        format!("{header}\nfree,57\nfree,57\nfree,89\n"),
        format!("{header}\nfree,57\nfree,89\nfree,60\n"),
        format!("{header}\npad,57,0.01,0.02,NaN\nfree,89\n"),
        format!("{header}\npad,57,0.01,{length},0.4\nfree,89\n"),
        format!("{header}\npad,57,{start},{end},0.4\nfree,89\n")] {
        assert!(LivePiano::new(&desc, damper_config(text)).is_err());
    }
    assert!(LivePiano::new(&desc, damper_config(format!("{header}\nfree,57\nfree,89\n"))).is_ok());
}
