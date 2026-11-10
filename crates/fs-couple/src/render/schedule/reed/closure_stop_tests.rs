//! The card row `closure_stop on|off` (the reed step's rigid closure stop at zero opening, commit a857c6f17, unconditional until
//! this row), on a card whose reed is pressed shut through a very soft lay. G0: the default and `closure_stop on` render the same
//! bytes; the row refuses a bad word and a massless reed. G1: with the stop the opening never goes below zero; without it the
//! soft lay lets the reed through and the pressure differs. G2: a reed that stays open renders the same bytes with the stop on
//! or off. Synthetic numbers: a byte-identity and sign test of the switch, not an instrument.
use super::*;

const HEAD: &str = "frankensim-reed-performance-v1\naudio 48000 4800 20000\nambient 293.15 101325 0\n";
// Blown at 9 kPa, 1.5 times the reed's 6 kPa closing pressure.
const TAIL: &str = "listener_m 1\ntermination flanged\nsegments 1\ncylinder 0.0073 0.5\n\
    compile_limits 100000 1024\nschedule\nfrankensim-gesture-schedule-v1\n\
    control_rate_hz\t700\ntracks\t1\ntrack\tblow\tblowing-pressure\n\
    initial\tpressure-pa\t0e0\nevents\t1\nevent\t0e0\t1e-2\tpressure-pa\t9e3\n";
// A massive reed with a very soft lay (stiffness 1e3, exponent 2, onset 0.2 mm, no loss).
const SOFT: &str = "reed 0.0004 0.013 6000 3.8e-6 1400 0.077 1e3 2 2e-4 0";
const MASSLESS: &str = "reed 0.0004 0.013 6000 0 0 0.35";

fn card(reed_row: &str, rows: &str) -> String {
    format!("{HEAD}{reed_row}\n{rows}{TAIL}")
}

fn refusal(text: &str) -> &'static str {
    match ReedPerformance::from_bytes(text.as_bytes(), 64) {
        Err(ReedPerformanceError::Input { what, .. }) => what,
        other => panic!("expected an input refusal, got {:?}", other.err()),
    }
}

/// The pressure bits of the whole window and the lowest reed opening seen at a block end.
fn run(text: &str) -> (Vec<u64>, f64, bool) {
    let mut performance = ReedPerformance::from_bytes(text.as_bytes(), 64).expect("the card must load");
    let mut bits = Vec::with_capacity(4800);
    let mut lowest = f64::INFINITY;
    let mut out = vec![0.0; 64];
    let mut stop = true;
    for _ in 0..75 {
        performance.block(&mut out).expect("the short window must render");
        bits.extend(out.iter().map(|x| x.to_bits()));
        let RenderVoice::ReedBore(voice) = &performance.renderer.context().voices[0] else {
            panic!("expected the original reed voice");
        };
        lowest = lowest.min(voice.reed_y);
        stop = voice.closure_stop();
    }
    (bits, lowest, stop)
}

#[test]
fn g0_the_default_and_closure_stop_on_render_the_same_bytes() {
    let (absent, lowest_absent, stop_absent) = run(&card(SOFT, ""));
    let (on, lowest_on, stop_on) = run(&card(SOFT, "closure_stop on\n"));
    assert_eq!(absent, on, "an explicit `on` must be the default, bit for bit");
    assert_eq!(lowest_absent.to_bits(), lowest_on.to_bits());
    assert!(stop_absent && stop_on);
}

#[test]
fn g1_without_the_stop_the_soft_lay_lets_the_reed_through_and_the_pressure_differs() {
    let (on, lowest_on, _) = run(&card(SOFT, "closure_stop on\n"));
    let (off, lowest_off, stop_off) = run(&card(SOFT, "closure_stop off\n"));
    assert!(!stop_off, "the row must reach the voice");
    assert!(lowest_on >= 0.0, "with the stop the opening never goes below zero: {lowest_on:e}");
    assert!(lowest_off < 0.0, "without it the soft lay lets the reed through: {lowest_off:e}");
    assert_ne!(on, off, "the two steps must differ once the reed closes");
    assert!(off.iter().all(|bits| f64::from_bits(*bits).is_finite()));
}

#[test]
fn g2_a_reed_that_stays_open_renders_the_same_bytes_with_the_stop_on_or_off() {
    // Blown at 1.5 kPa, a quarter of the closing pressure and below the oscillation threshold: the reed settles at its
    // static opening and never reaches zero, so the stop never acts and `off` must be the same step as `on`, bit for bit.
    let low = |rows: &str| card(SOFT, rows).replace("pressure-pa\t9e3", "pressure-pa\t1.5e3");
    let (on, lowest_on, _) = run(&low("closure_stop on\n"));
    let (off, lowest_off, stop_off) = run(&low("closure_stop off\n"));
    assert!(!stop_off, "the row must reach the voice");
    assert!(lowest_on > 0.0 && lowest_off > 0.0, "the reed must stay open: {lowest_on:e} {lowest_off:e}");
    assert_eq!(on, off, "a step whose free solution stays at or above zero is the same step either way");
}

#[test]
fn g0_closure_stop_refuses_a_bad_word_and_a_massless_reed() {
    assert_eq!(refusal(&card(SOFT, "closure_stop maybe\n")), "closure_stop takes on or off");
    assert_eq!(refusal(&card(MASSLESS, "closure_stop off\n")), "closure_stop requires a massive reed");
    // A trailing field is refused by the row reader (its own message).
    let _ = refusal(&card(SOFT, "closure_stop off extra\n"));
}
