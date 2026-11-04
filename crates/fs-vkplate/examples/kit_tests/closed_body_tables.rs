//! The committed closed-body tables (`data/closed-body/`, baked by fs-bem's
//! `examples/drum_observer.rs`, each reproducible byte for byte by its `--verify`) attach to the
//! drums exactly as this example builds them for a render (`drum_cards::drum`, the card heads,
//! the renders' own stick speeds, the oversampling that picks): rate, gas, microphone, floor,
//! radius, depth and the kernels fs-bem reconstructs all bind, as `--observer` needs them to.
use super::*;
use fs_vkplate::closed_body::{ClosedBodyTable, RECIPROCITY_LIMIT};

#[test]
fn the_committed_tables_attach_to_the_kit_drums_as_rendered() {
    let tables: [(&str, kit::DrumSpec, bool, f64, &str); 5] = [
        ("kick", kit::kick(), true, 3.0, include_str!("../../data/closed-body/kick.tsv")),
        ("snare", kit::snare(), false, 4.5, include_str!("../../data/closed-body/snare.tsv")),
        ("snare-off", kit::snare_wires_off(), false, 4.5, include_str!("../../data/closed-body/snare-off.tsv")),
        ("tom-12", kit::tom(12.0), false, 3.5, include_str!("../../data/closed-body/tom-12.tsv")),
        ("tom-14", kit::tom(14.0), false, 3.5, include_str!("../../data/closed-body/tom-14.tsv")),
    ];
    for (name, spec, kick, speed, text) in tables {
        let table = ClosedBodyTable::parse(text).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(table.reciprocity <= RECIPROCITY_LIMIT, "{name}: reciprocity {}", table.reciprocity);
        let mut p = drum_cards::drum(&spec, 1.0, &[(0.005, speed)], kick, false).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(table.rate, p.instrument.rate, "{name}: the table's clock is the rendered piece's");
        kit::attach_closed_body(&mut p, &spec, text).unwrap_or_else(|e| panic!("{name}: {e}"));
        let observer = p.instrument.mics[0].closed_body().unwrap_or_else(|| panic!("{name}: no observer on the overhead"));
        eprintln!("{name}: {} kHz, FFT {}, latency {} samples; {}", table.rate / 1000.0, table.n_fft, observer.latency(), table.receipt());
    }
}

/// `--head-dipole`'s route: the table's length and the air modes' f_max reach the card build through `drum_cards::drum`
/// (the spec survives `card_drum` and `retune`), and the card's tension sets the law's frequency.
#[test]
fn the_head_dipole_survives_the_card_route() {
    let text = include_str!("../../data/closed-body/snare-off.tsv");
    let head_dipole = Some(kit::head_dipole_from_table(text, Some(8000.0)).unwrap());
    let spec = kit::DrumSpec { head_dipole, ..kit::snare_wires_off() };
    let resistance = |tune: f64| {
        let p = drum_cards::drum(&spec, tune, &[(0.005, 4.5)], false, false).unwrap();
        for n in p.notes.iter().filter(|n| n.starts_with("head dipole")) {
            eprintln!("tune {tune}: {n}");
        }
        p.instrument.rad_dipole.as_ref().expect("the spec asked for the dipole").resistance
    };
    let (card, tuned) = (resistance(1.0), resistance(1.5));
    assert!(card > 0.0);
    // Omega^4 at the two hybrids; the mirrored solve puts them at 152.40 and 183.58 Hz.
    assert!((tuned / card / (183.5781_f64 / 152.3994).powi(4) - 1.0).abs() < 1e-3, "{card} {tuned}");
    let without = drum_cards::drum(&kit::snare_wires_off(), 1.0, &[(0.005, 4.5)], false, false).unwrap();
    assert!(without.instrument.rad_dipole.is_none());
}
