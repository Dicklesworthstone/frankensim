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
