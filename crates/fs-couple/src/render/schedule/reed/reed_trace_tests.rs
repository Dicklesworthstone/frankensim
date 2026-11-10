//! The reed opening tap (MB 22212 item 1): an opt-in line `y v` per solver sample, the blade's opening [m] and velocity [m/s] as the closure stop and the lay see them.
//! It must change no audio sample, write one line per sample, and show the stop: with it on the opening never goes below zero; with `closure_stop off` on a very soft lay
//! pressed shut by a blow of 1.5 times the closing pressure it does. Synthetic numbers, as in `closure_stop_tests`.
use super::*;
use std::sync::{Arc, Mutex};

const HEAD: &str = "frankensim-reed-performance-v1\naudio 48000 4800 20000\nambient 293.15 101325 0\n";
const TAIL: &str = "listener_m 1\ntermination flanged\nsegments 1\ncylinder 0.0073 0.5\n\
    compile_limits 100000 1024\nschedule\nfrankensim-gesture-schedule-v1\n\
    control_rate_hz\t700\ntracks\t1\ntrack\tblow\tblowing-pressure\n\
    initial\tpressure-pa\t0e0\nevents\t1\nevent\t0e0\t1e-2\tpressure-pa\t9e3\n";
const SOFT: &str = "reed 0.0004 0.013 6000 3.8e-6 1400 0.077 1e3 2 2e-4 0";

fn card(rows: &str) -> String {
    format!("{HEAD}{SOFT}\n{rows}{TAIL}")
}

struct Tap(Arc<Mutex<Vec<u8>>>);

impl std::io::Write for Tap {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.lock().expect("tap buffer").extend_from_slice(buf);
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The pressure bits of the whole window, and the tap's lines when a tap was attached.
fn run(text: &str, tapped: bool) -> (Vec<u64>, Vec<(f64, f64)>) {
    let mut performance = ReedPerformance::from_bytes(text.as_bytes(), 64).expect("the card must load");
    let buffer = Arc::new(Mutex::new(Vec::new()));
    if tapped {
        performance.record_reed_opening(Box::new(Tap(buffer.clone()))).expect("a massive reed takes the tap");
    }
    let mut bits = Vec::with_capacity(4800);
    let mut out = vec![0.0; 64];
    for _ in 0..75 {
        performance.block(&mut out).expect("the short window must render");
        bits.extend(out.iter().map(|x| x.to_bits()));
    }
    let text = String::from_utf8(buffer.lock().expect("tap buffer").clone()).expect("the tap writes text");
    let lines = text
        .lines()
        .map(|line| {
            let mut words = line.split_ascii_whitespace().map(|w| w.parse::<f64>().expect("a number"));
            (words.next().expect("y"), words.next().expect("v"))
        })
        .collect();
    (bits, lines)
}

#[test]
fn the_opening_tap_changes_no_audio_and_writes_one_line_per_sample() {
    for rows in ["closure_stop on\n", "closure_stop off\n", ""] {
        let (plain, none) = run(&card(rows), false);
        let (tapped, lines) = run(&card(rows), true);
        assert_eq!(plain, tapped, "the tap must leave every sample of the audio as it was ({rows:?})");
        assert!(none.is_empty());
        assert_eq!(lines.len(), 4800, "one line per solver sample ({rows:?})");
        assert!(lines.iter().all(|(y, v)| y.is_finite() && v.is_finite()));
    }
}

#[test]
fn the_tap_shows_the_stop_holding_the_opening_at_zero_and_the_off_row_letting_it_through() {
    let lowest = |rows: &str| run(&card(rows), true).1.iter().map(|(y, _)| *y).fold(f64::INFINITY, f64::min);
    let (on, off) = (lowest("closure_stop on\n"), lowest("closure_stop off\n"));
    assert!(on >= 0.0, "with the stop the opening never goes below zero: {on:e}");
    assert!(off < 0.0, "without it the soft lay lets the reed through: {off:e}");
}

#[test]
fn the_tap_refuses_a_massless_reed_and_an_advanced_voice() {
    let massless = format!("{HEAD}reed 0.0004 0.013 6000 0 0 0.35\n{TAIL}");
    let mut p = ReedPerformance::from_bytes(massless.as_bytes(), 64).expect("a massless reed loads");
    assert!(p.record_reed_opening(Box::new(Tap(Arc::new(Mutex::new(Vec::new()))))).is_err());
    let mut q = ReedPerformance::from_bytes(card("").as_bytes(), 64).expect("the card loads");
    let mut out = vec![0.0; 64];
    q.block(&mut out).expect("one block");
    assert!(q.record_reed_opening(Box::new(Tap(Arc::new(Mutex::new(Vec::new()))))).is_err(), "an advanced voice refuses the tap");
}
