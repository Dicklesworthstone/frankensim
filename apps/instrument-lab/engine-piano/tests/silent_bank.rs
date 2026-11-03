//! The silent-string bank (`Config::silent_bank_hz`): every partial up to the band of every
//! string that is not sounding moves, driven by the bridge.
//! - Against the same partials attached as dynamic strings (every other key held down
//!   silently, the partial band at the bank's band), both under the full sustain pedal so
//!   no pad touches either, an A3 strike sounds the same to within -50 dB (measured -57.7).
//! - The books close with the bank's explicit exchange booked, through a struck key handing
//!   its bank partials to the dynamic set, returning them when it retires, and a re-strike.
use fs_live_piano::{Config, Description, LivePiano};

const RATE: u32 = 48_000;

fn render(p: &mut LivePiano, seconds: f64) -> Vec<f32> {
    let mut out = vec![0.0_f32; (seconds * f64::from(RATE)) as usize];
    p.render(&mut out);
    out
}

#[test]
fn the_bank_sounds_as_the_same_partials_attached_as_strings() {
    let desc = Description::embedded().unwrap();
    let band = 300.0;
    let mut cfg = Config::live(RATE);
    cfg.max_courses = 88;
    cfg.max_dynamic_modes = 1 << 20;
    cfg.retire_ringing_j = 0.0;
    cfg.retire_damped_j = 0.0;

    let mut strings = LivePiano::new(&desc, cfg.clone()).unwrap();
    strings.set_sustain(1.0).unwrap();
    strings.set_quality(200, band, 1 << 20).unwrap();
    // Only keys with partials in the band: an attached string always keeps at least one.
    let held: Vec<u8> = (21..=108).filter(|&m| m != 57 && strings.partials_hz(m)[0] <= band).collect();
    for midi in held {
        strings.silent_key_down(midi).unwrap();
    }
    strings.apply_pending_notes();
    strings.set_quality(200, 0.45 * f64::from(RATE), 1 << 20).unwrap();
    strings.note_on(57, 2.0).unwrap();
    let a = render(&mut strings, 1.0);

    let mut cfg_bank = cfg;
    cfg_bank.silent_bank_hz = band;
    let mut bank = LivePiano::new(&desc, cfg_bank).unwrap();
    bank.set_sustain(1.0).unwrap();
    bank.note_on(57, 2.0).unwrap();
    let b = render(&mut bank, 1.0);

    let (mut num, mut den) = (0.0, 0.0);
    for (x, y) in a.iter().zip(&b) {
        num += f64::from(x - y).powi(2);
        den += f64::from(*x).powi(2);
    }
    let db = 10.0 * (num / den).log10();
    println!("bank vs attached strings: {db:.1} dB");
    assert!(db < -50.0, "the bank differs from the attached strings by {db:.1} dB");
    // And the bank changes the strike: without it the silent strings are rigid chords.
    let mut none = LivePiano::new(&desc, Config::live(RATE)).unwrap();
    none.set_sustain(1.0).unwrap();
    none.note_on(57, 2.0).unwrap();
    let c = render(&mut none, 1.0);
    let moved: f64 = c.iter().zip(&a).map(|(x, y)| f64::from(x - y).powi(2)).sum::<f64>() / den;
    assert!(10.0 * moved.log10() > -20.0, "the silent strings must matter here");
}

#[test]
fn the_books_close_with_the_bank_exchange_booked() {
    let mut cfg = Config::live(RATE);
    cfg.silent_bank_hz = 600.0;
    cfg.rt0425_string_damping = true;
    cfg.rt0425_hammer_stiffness = true;
    cfg.rt0425_hammer_dissipation = true;
    let mut p = LivePiano::new(&Description::embedded().unwrap(), cfg).unwrap();
    p.note_on(57, 3.0).unwrap();
    p.note_on(45, 2.0).unwrap();
    render(&mut p, 0.4);
    p.note_off(57);
    p.note_off(45);
    render(&mut p, 2.0);
    assert!(p.counters.retirements > 0, "the struck keys must have returned their partials to the bank");
    p.note_on(57, 2.0).unwrap();
    p.set_sustain(1.0).unwrap();
    render(&mut p, 0.5);
    let l = p.close_ledger();
    println!("{l:?}");
    assert!(l.residual_j.abs() < 1.0e-9 * l.input_j, "residual {:e} of {:e} J", l.residual_j, l.input_j);
    // Measured 5.5e-4 of the input.
    assert!(l.coupling_j.abs() < 2.0e-3 * l.input_j, "coupling {:e} of {:e} J", l.coupling_j, l.input_j);
    assert!(l.damper_loss_j > 0.0 && p.take_error().is_none());
}

#[test]
fn silent_strings_use_supplied_pads_above_the_estimated_break_and_the_same_pedal_factor() {
    let mut desc = Description::embedded().unwrap().truncated(4).unwrap();
    desc.courses.retain(|c| [57, 89].contains(&c.midi));
    let length = desc.courses[1].length_m;
    let make = |pad: Option<(f64, f64, f64)>, sustain: f64| {
        let row = pad.map_or_else(|| "free,89".into(), |(start, end, drag)|
            format!("pad,89,{},{},{drag}", start * length, end * length));
        let mut cfg = Config::live(RATE);
        cfg.max_partials = 16;
        cfg.course_mode_budget = 48;
        cfg.duplex = false;
        cfg.silent_bank_hz = 2_000.0;
        cfg.retire_ringing_j = 0.0;
        cfg.retire_damped_j = 0.0;
        cfg.damper_specification = Some(format!("frankensim-piano-dampers-v1\nfree,57\n{row}\n"));
        let mut piano = LivePiano::new(&desc, cfg).unwrap();
        piano.set_sustain(sustain).unwrap();
        piano.note_on(57, 2.0).unwrap();
        piano
    };
    let mut half_pedal = make(Some((0.31, 0.39, 0.4)), 0.5);
    let mut quarter_drag = make(Some((0.31, 0.39, 0.1)), 0.0);
    let mut moved = make(Some((0.21, 0.29, 0.1)), 0.0);
    let mut free = make(None, 0.0);
    let mut raised = make(Some((0.31, 0.39, 0.4)), 1.0);
    let a = render(&mut half_pedal, 0.08);
    assert!(render(&mut quarter_drag, 0.08) == a);
    let b = render(&mut moved, 0.08);
    let c = render(&mut free, 0.08);
    assert!(render(&mut raised, 0.08) == c, "a fully raised supplied pad acts as a free row");
    assert!(half_pedal.grant_for(89).is_none(), "the padded string must still be owned by the silent bank");
    assert!(half_pedal.ledger.damper_loss_j > 0.0, "an unplayed pad above MIDI 88 must dissipate bridge-driven motion");
    assert_eq!(half_pedal.ledger.damper_loss_j.to_bits(), quarter_drag.ledger.damper_loss_j.to_bits());
    assert!(moved.ledger.damper_loss_j > 0.0 && b != a, "the bank must use the supplied station geometry");
    assert_eq!(free.ledger.damper_loss_j, 0.0);
    assert_eq!(raised.ledger.damper_loss_j, 0.0);

    // The same pad remains attached when the sympathetic key moves from the
    // silent bank to the dynamic set and is caught, then released by the rod.
    let before = half_pedal.stored_energy_j();
    half_pedal.silent_key_down(89).unwrap();
    half_pedal.apply_pending_notes();
    assert!(half_pedal.grant_for(89).is_some());
    assert!((half_pedal.stored_energy_j() - before).abs() < 1e-10 * half_pedal.ledger.input_j);
    half_pedal.set_sostenuto(true);
    half_pedal.note_off(89);
    let loss = half_pedal.ledger.damper_loss_j;
    render(&mut half_pedal, 0.02);
    assert_eq!(half_pedal.ledger.damper_loss_j, loss);
    half_pedal.set_sostenuto(false);
    render(&mut half_pedal, 0.02);
    assert!(half_pedal.ledger.damper_loss_j > loss);
    for piano in [&mut half_pedal, &mut quarter_drag, &mut moved, &mut free, &mut raised] {
        let ledger = piano.close_ledger();
        assert!(ledger.residual_j.abs() < 1e-8 * ledger.input_j);
        assert!(piano.take_error().is_none());
    }
}

#[test]
fn a_supplied_pad_cannot_silently_overdamp_the_silent_bank_approximation() {
    let mut desc = Description::embedded().unwrap().truncated(4).unwrap();
    desc.courses.retain(|c| c.midi == 57);
    let length = desc.courses[0].length_m;
    let mut cfg = Config::live(RATE);
    cfg.max_partials = 16;
    cfg.damper_specification = Some(format!("frankensim-piano-dampers-v1\npad,57,{},{},1000000\n",
        0.31 * length, 0.39 * length));
    // The active pad's exact viscous flows admit this physical source. The
    // silent bank keeps its existing underdamped transition domain and refuses.
    assert!(LivePiano::new(&desc, cfg.clone()).is_ok());
    cfg.silent_bank_hz = 600.0;
    let error = LivePiano::new(&desc, cfg).err().unwrap();
    assert!(error.contains("silent-bank partial would be overdamped"), "{error}");
}
