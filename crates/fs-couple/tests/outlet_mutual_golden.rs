//! The `outlet_mutual free-field` row keeps its bits: the two posted flute
//! cards that carry it (D4 and G4) and the recorder a1 card with the row added,
//! each cut to a 2048-sample window with its first gesture event only, render
//! the pressure bits recorded on the code before `outlet_mutual pritchard`
//! (FNV-1a 64 of the little-endian f64 bits, printed by this file's first
//! version on that code; the recorder card's were recorded again on 2026-11-09
//! at ca28955c9, when the card took the measured windway height); without the
//! row the same window renders other bits.
use fs_couple::flue_loop::network::outlet_mutual::OutletMutualForm;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;

const FLUTE_D4: &str = include_str!(
    "../../../data/instruments/descriptions/flute/flute-d4-pf661-5ft-unflanged-exact-held-flow-bore-and-chimney-loss-outlet-mutual.performance");
const FLUTE_G4: &str = include_str!(
    "../../../data/instruments/descriptions/flute/flute-g4-pf661-5ft-unflanged-exact-held-flow-bore-and-chimney-loss-outlet-mutual.performance");
const RECORDER_A1: &str = include_str!("../../../data/instruments/descriptions/recorder/recorder-a1.performance");
const LAST_OUTLET: &str = "outlet 21 -0.0863000 0\n";
const FREE: &str = "\noutlet_mutual free-field\n";
const WINDOW: usize = 2048;
/// Digests recorded before `pritchard` existed; recorder-a1's again on 2026-11-09 at
/// ca28955c9 for the card at the measured windway height (h 1.05428 mm).
const GOLDEN: [(&str, u64); 3] =
    [("flute-d4", 0x8027_0d04_9e43_fb48), ("flute-g4", 0x1acd_072d_973e_4078), ("recorder-a1", 0x480a_2b02_aa03_de52)];

fn gate() -> CancelGate { CancelGate::new_clock_free() }

/// The card with `WINDOW` output samples at full scale 1 and its first gesture event only.
fn cut(text: &str) -> String {
    let start = text.find("\noutput ").expect("an output row") + 1;
    let end = start + text[start..].find('\n').expect("output row end");
    let fields: Vec<&str> = text[start..end].split_ascii_whitespace().collect();
    assert_eq!(fields.len(), 4, "output row: {fields:?}");
    let text = format!("{}output {WINDOW} 1 {}{}", &text[..start], fields[3], &text[end..]);
    let at = text.find("\nevents\t").expect("an events row");
    let (head, tail) = text.split_at(at);
    let lines: Vec<&str> = tail[1..].split('\n').collect();
    assert!(lines[1].starts_with("event\t"), "{}", lines[1]);
    let newline = if tail.ends_with('\n') { "\n" } else { "" };
    format!("{head}\nevents\t1\n{}{newline}", lines[1])
}

fn digest(pressure: &[f64]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for p in pressure {
        for byte in p.to_bits().to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
    }
    h
}

fn card(name: &str) -> String {
    match name {
        "flute-d4" => cut(FLUTE_D4),
        "flute-g4" => cut(FLUTE_G4),
        "recorder-a1" => {
            assert!(RECORDER_A1.contains(LAST_OUTLET));
            cut(&RECORDER_A1.replace(LAST_OUTLET, &format!("{LAST_OUTLET}outlet_mutual free-field\n")))
        }
        _ => unreachable!("{name}"),
    }
}

/// The card's mutual form (None without the row) and its rendered window's digest.
fn render(text: &str) -> (Option<OutletMutualForm>, u64) {
    let file = FlueFilePerformance::from_bytes(text.as_bytes(), 256, &gate()).unwrap_or_else(|e| panic!("{e}"));
    let form = file.renderer().voice().network().unwrap().outlet_mutual().map(|m| m.record().form);
    let mut renderer = file.into_renderer();
    let mut pressure = vec![0.0_f64; WINDOW];
    for block in pressure.chunks_mut(256) { renderer.block(block).unwrap(); }
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-10));
    (form, digest(&pressure))
}

fn check(name: &str) {
    let want = GOLDEN.iter().find(|(n, _)| *n == name).expect("a golden digest").1;
    let text = card(name);
    assert_eq!(text.matches(FREE).count(), 1, "{name}");
    let (form, got) = render(&text);
    assert_eq!(form, Some(OutletMutualForm::FreeField));
    assert_eq!(got, want, "{name}: the free-field row's bits moved ({got:016x} against {want:016x})");
    // Planted negative: the same window without the row is not the golden one, so the digest sees the coupling.
    let (form, plain) = render(&text.replace(FREE, "\n"));
    assert_eq!(form, None);
    assert_ne!(plain, want, "{name}: the row did not reach the window");
}

#[test]
fn flute_d4_keeps_its_free_field_bits() { check("flute-d4"); }
#[test]
fn flute_g4_keeps_its_free_field_bits() { check("flute-g4"); }
#[test]
fn recorder_a1_keeps_its_free_field_bits() { check("recorder-a1"); }
