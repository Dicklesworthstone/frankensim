//! The silent-string bank (`Config::silent_bank_hz`): every partial up to the band of every
//! string that is not sounding moves, driven by the bridge.
//! - Against the same partials attached as dynamic strings (every other key held down
//!   silently, the partial band at the bank's band), both under the full sustain pedal so
//!   no pad touches either, an A3 strike sounds the same to within -30 dB.
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
    assert!(db < -30.0, "the bank differs from the attached strings by {db:.1} dB");
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
    assert!(l.coupling_j.abs() < 1.0e-2 * l.input_j, "coupling {:e} of {:e} J", l.coupling_j, l.input_j);
    assert!(l.damper_loss_j > 0.0 && p.take_error().is_none());
}
