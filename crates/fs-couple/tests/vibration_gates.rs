//! Vibration claim gates (music bead `frankensim-music-v8-root-3ez8g.7.1`):
//! the truncation-honesty gate with its under-truncated falsifier, the
//! selector-threshold-as-data pin, and the gate-summary enumeration —
//! the review that promotes the string/plate/nonlinear registry rows.
//! The Gonzalez caveat travels with the KC/von-Kármán row: conservation
//! under discrete-gradient stepping is structurally blind to gradient
//! errors, so nonlinear evidence cites FD-gradient/trajectory oracles,
//! never energy conservation.

use fs_couple::piano_vertical::{HammerLaw, PedalState, PianoStringSpec, PianoVertical};

fn spec_with_modes(n_modes: usize) -> PianoStringSpec {
    PianoStringSpec {
        f0_hz: 220.0,
        b_inharmonicity: 3.5e-4,
        detune_cents: 0.0,
        n_modes,
        damping_ratio: 4.0e-4,
    }
}

fn top_mode_energy_share(n_modes: usize) -> f64 {
    let mut pv = PianoVertical::new(
        spec_with_modes(n_modes),
        0.0,
        HammerLaw::Felt,
        PedalState {
            sustain: true,
            una_corda: false,
        },
        0.03,
    )
    .expect("vertical");
    pv.strike(2.0);
    // Past the contact (a few ms) but early enough that damping has not
    // reshaped the distribution.
    for _ in 0..4_800usize {
        let _ = pv.step();
    }
    let e = pv.middle_string_modal_energies();
    let total: f64 = e.iter().sum();
    assert!(total > 0.0, "the strike must deposit energy");
    e[n_modes - 1] / total
}

#[test]
fn vg_001_truncation_gate_and_under_truncated_falsifier() {
    // THE TRUNCATION GATE: a modal row is honest only if its retained
    // series has CONVERGED — the top retained mode must carry a
    // negligible share of the string energy, so the (absent, disclosed)
    // modes above it would carry less still. The disclosed fixture
    // (N=12 at 220 Hz under the felt hammer) passes; the deliberately
    // under-truncated config (N=3: the top mode sits mid-band where the
    // felt still delivers) FAILS THE SAME GATE — an executed falsifier,
    // not a vacuous check.
    let share_full = top_mode_energy_share(12);
    let share_under = top_mode_energy_share(3);
    // Gate authored from measurement: the disclosed fixture measured
    // 6.2e-20 (the ~2 ms felt contact simply never reaches partial 12
    // at 2.7 kHz) and the falsifier 1.0e-1 — the gate sits decisively
    // between, ~17 orders above the pass side and 50x below the fail
    // side, so neither a fixture re-dimension nor solver noise can
    // flip it silently.
    let gate = 2.0e-3;
    assert!(
        share_full < gate,
        "the disclosed N=12 fixture must pass its own truncation gate \
         (top-mode share {share_full:.3e} vs gate {gate:.1e})"
    );
    assert!(
        share_under > 10.0 * gate,
        "the under-truncated N=3 falsifier must fail the gate DECISIVELY \
         (top-mode share {share_under:.3e} vs gate {gate:.1e})"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"vg-001-truncation\",\"verdict\":\"pass\",\
         \"share_n12\":{share_full:.3e},\"share_n3\":{share_under:.3e},\"gate\":{gate:.1e}}}"
    );
}

#[test]
fn vg_002_selector_threshold_is_self_verified_data() {
    // The linear<->nonlinear switch point is COMMITTED DATA with fixture
    // provenance, and this gate RE-DERIVES it from the KC glide law
    // (3/32)(EA/T0)(k pi A/L)^2 * (1200/ln 2) cents — the same law
    // fs-nlmodal pins analytically at 1e-12. A hand-edited threshold
    // fails here.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf();
    let tsv = std::fs::read_to_string(root.join("data/claims/vibration-selector-thresholds.tsv"))
        .expect("committed selector thresholds");
    let field = |name: &str| -> f64 {
        tsv.lines()
            .find_map(|l| l.strip_prefix(&format!("{name}\t")))
            .unwrap_or_else(|| panic!("missing field {name}"))
            .parse()
            .expect("numeric")
    };
    let (target_cents, ea, t0, len, k) = (
        field("target_glide_cents"),
        field("ea_n"),
        field("t0_n"),
        field("length_m"),
        field("mode_k"),
    );
    let committed = field("amplitude_threshold_m");
    let r_star = target_cents / (1200.0 / core::f64::consts::LN_2);
    let derived =
        (len / (k * core::f64::consts::PI)) * (r_star / ((3.0 / 32.0) * (ea / t0))).sqrt();
    let rel = ((committed - derived) / derived).abs();
    assert!(
        rel < 1.0e-8,
        "the committed threshold must re-derive from its own provenance \
         (committed {committed:.9e} vs derived {derived:.9e}, rel {rel:.2e})"
    );
    // Sanity: the threshold sits in the physically sensible band for a
    // 0.65 m nylon-like string — millimetres, not microns or metres.
    assert!(
        (1.0e-3..3.0e-2).contains(&committed),
        "threshold {committed:.3e} m out of the plausible band"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"vg-002-selector\",\"verdict\":\"pass\",\
         \"amplitude_threshold_m\":{committed:.6e},\"target_glide_cents\":{target_cents}}}"
    );
}

#[test]
fn vg_003_vibration_gate_summary_enumerates_every_row() {
    // The committed vibration gate summary must name every registry row
    // for the vibration fillings — an omitted row is a red test — and
    // the honesty statuses must MATCH the registry (a summary that says
    // green while the registry says ungated is the lie this gate
    // exists to catch).
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf();
    let registry = std::fs::read_to_string(root.join("instrument-claims.json")).expect("registry");
    let summary = std::fs::read_to_string(root.join("data/claims/vibration-gate-summary.tsv"))
        .expect("committed vibration gate summary");
    assert!(summary.starts_with("# frankensim-vibration-gate-summary-v1"));
    let fillings = [
        "string",
        "plate",
        "string-plate-nl",
        "bowed",
        "string-plate-cavity",
    ];
    let mut rows = 0usize;
    let mut cursor = 0usize;
    while let Some(hit) = registry[cursor..].find("\"filling\": \"") {
        let at = cursor + hit + "\"filling\": \"".len();
        let end = registry[at..].find('"').expect("end");
        let filling = &registry[at..at + end];
        cursor = at + end;
        if !fillings.contains(&filling) {
            continue;
        }
        let tail = &registry[cursor..];
        let grab = |tag: &str| -> String {
            let t = format!("\"{tag}\": \"");
            let s = tail.find(&t).expect(tag) + t.len();
            let e = tail[s..].find('"').expect("end");
            tail[s..s + e].to_string()
        };
        let (image, gate) = (grab("image"), grab("gate"));
        rows += 1;
        let line = summary
            .lines()
            .find(|l| l.starts_with(&format!("{filling}\t{image}\t")))
            .unwrap_or_else(|| panic!("summary omits {filling}/{image}"));
        assert!(
            line.split('\t').nth(3) == Some(gate.as_str()),
            "summary status for {filling}/{image} must match the registry ({gate})"
        );
    }
    assert!(rows >= 7, "the vibration menu: {rows} rows");
    // The named-missing honesty row: the plate 1-port stays ungated
    // until an EXECUTED fit artifact exists.
    assert!(
        summary.contains("vfit-driving-point\tpassive-one-port\tungated"),
        "the plate 1-port must stay honestly ungated"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"vg-003-gate-summary\",\"verdict\":\"pass\",\
         \"rows\":{rows}}}"
    );
}
