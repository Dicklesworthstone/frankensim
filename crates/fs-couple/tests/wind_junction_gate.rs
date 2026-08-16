//! Wind char-vs-TMM junction gate (music bead
//! `frankensim-music-v8-root-3ez8g.6.1`): the characteristic image's
//! realized reflectance must MATCH the TMM authority's `R(omega)` on a
//! holed bore, within an authored band — the menu's recorded failure
//! mode ("junctions drift from TMM R") made an executable gate, with the
//! FALSIFIER: a deliberately detuned hole model (the Ernoult-magnitude
//! 15 mm position error) must FAIL the same gate. Both images consume
//! the SAME typed `FingeringTable` (the chart-side fingering data —
//! never inline hole states).

use fs_couple::driving_point::characteristic_line;
use fs_duct::{
    Duct, Fingering, FingeringTable, HoleState, LossModel, Segment, Termination, input_impedance,
};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_vfit::discretize::reflectance;

const RATE: u32 = 48_000;

fn air20() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

/// The Ernoult-class template with `detune_m` added to hole 2's position
/// (0 = the honest geometry; 0.015 = the falsifier's data error).
fn template(detune_m: f64) -> Duct {
    let bore = 0.002f64;
    let hole = |r: f64, chimney: f64| Segment::ToneHole {
        hole_radius: r,
        chimney_height: chimney,
        bore_radius: bore,
        state: HoleState::Closed,
    };
    Duct {
        segments: vec![
            Segment::Cylinder {
                radius: bore,
                length: 0.100,
            },
            hole(0.0015, 0.0017),
            Segment::Cylinder {
                radius: bore,
                length: 0.030,
            },
            hole(0.00175, 0.0013),
            Segment::Cylinder {
                radius: bore,
                length: 0.050 + detune_m,
            },
            hole(0.00175, 0.0015),
            Segment::Cylinder {
                radius: bore,
                length: 0.060 - detune_m,
            },
            hole(0.00125, 0.0014),
            Segment::Cylinder {
                radius: bore,
                length: 0.0475,
            },
        ],
    }
}

fn table(detune_m: f64) -> FingeringTable {
    let closed = |states: [HoleState; 4]| states.to_vec();
    FingeringTable::try_new(
        template(detune_m),
        vec![
            Fingering {
                label: "xxxx".to_string(),
                holes: closed([HoleState::Closed; 4]),
            },
            Fingering {
                label: "xxxo".to_string(),
                holes: closed([
                    HoleState::Closed,
                    HoleState::Closed,
                    HoleState::Closed,
                    HoleState::Open,
                ]),
            },
            Fingering {
                label: "xxox".to_string(),
                holes: closed([
                    HoleState::Closed,
                    HoleState::Closed,
                    HoleState::Open,
                    HoleState::Closed,
                ]),
            },
        ],
    )
    .expect("table admits")
}

/// Worst |R_char - R_tmm| over the working band: the char image's
/// realized taps (from `char_duct`) DTFT'd back to the frequency domain
/// against the TMM reflectance of `tmm_duct`. The two ducts are the
/// SAME object in the gate; they differ only in the falsifier, which
/// models the char image's junction handling drifting from the
/// authority (an A-vs-B consistency gate is structurally blind to a
/// SHARED geometry error — the recorded metamorphic-blindness law — so
/// the falsifier perturbs one arm, not the shared input).
fn worst_reflectance_gap(char_duct: &Duct, tmm_duct: &Duct, gas: &GasState) -> f64 {
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * 0.002 * 0.002);
    let line = characteristic_line(
        char_duct,
        gas,
        Termination::UnflangedOpen,
        RATE,
        4096,
        zc,
        None,
    )
    .expect("line realizes");
    let taps = line_taps(&line);
    let dt = 1.0 / f64::from(RATE);
    let mut worst = 0.0f64;
    for i in 0..60 {
        let f = 200.0 + 900.0 * f64::from(i) / 59.0;
        let omega = core::f64::consts::TAU * f;
        let z = input_impedance(
            tmm_duct,
            gas,
            omega,
            LossModel::Bessel,
            Termination::UnflangedOpen,
        )
        .expect("tmm")
        .impedance;
        let r_tmm = reflectance(C64::new(z.re, z.im), zc);
        // The realized IR lives on the e^{+iwt} DFT side (conjugated at
        // realization); conjugate back for the comparison.
        let mut re = 0.0f64;
        let mut im = 0.0f64;
        for (k, &h) in taps.iter().enumerate() {
            let ph = omega * dt * k as f64;
            re += h * ph.cos();
            im -= h * ph.sin();
        }
        let r_char = C64::new(re, -im);
        worst = worst.max((r_char - r_tmm).abs());
    }
    worst
}

/// The line's FIR taps via the public history-free surface: the impulse
/// response is re-read by pushing a unit impulse through a CLONE.
fn line_taps(line: &fs_vfit::discretize::DelayedFilter) -> Vec<f64> {
    let mut probe = line.clone();
    let mut taps = Vec::with_capacity(4096);
    for k in 0..4096usize {
        let x = if k == 0 { 1.0 } else { 0.0 };
        taps.push(probe.push(x).expect("push"));
    }
    taps
}

#[test]
fn char_junctions_match_tmm_r_and_the_detuned_hole_fails() {
    let gas = air20();
    let honest = table(0.0);
    let detuned = table(0.015);
    // The gate: every fingering's char reflectance tracks the TMM R
    // within the authored band.
    let band = 0.08f64;
    let mut worst_honest = 0.0f64;
    for label in honest.labels() {
        let duct = honest.duct(label).expect("fingering");
        let gap = worst_reflectance_gap(&duct, &duct, &gas);
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"wind-junction-gate\",\"fingering\":\"{label}\",\
             \"worst_reflectance_gap\":{gap:.4}}}"
        );
        worst_honest = worst_honest.max(gap);
    }
    // THE FALSIFIER: the char image's effective junction geometry
    // drifting 15 mm (Ernoult magnitude) from the TMM authority must
    // exceed the band on an open-hole fingering. Perturbing ONE arm is
    // the point: perturbing the shared duct leaves both arms agreeing
    // (executed on the first run — 0.0185 vs band 0.08 — the
    // metamorphic-blindness law in action, now documented here).
    let duct_bad = detuned.duct("xxox").expect("fingering");
    let duct_ref = honest.duct("xxox").expect("fingering");
    let gap_bad = worst_reflectance_gap(&duct_bad, &duct_ref, &gas);
    // An unknown fingering refuses by name.
    let unknown = honest.duct("oooo").is_err();
    assert!(
        worst_honest < band,
        "honest fingerings must sit inside the band ({worst_honest:.4} >= {band})"
    );
    assert!(
        gap_bad > band,
        "the detuned falsifier must FAIL the gate ({gap_bad:.4} <= {band})"
    );
    assert!(unknown, "unknown fingering label must refuse");
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"wind-junction-gate\",\"verdict\":\"pass\",\
         \"worst_honest\":{worst_honest:.4},\"detuned_falsifier\":{gap_bad:.4},\
         \"band\":{band}}}"
    );
}

#[test]
fn wind_gate_summary_enumerates_every_registry_row() {
    // The committed wind gate summary (3ez8g.6.1 logging clause) must
    // name every wind-reed row — an omitted row is a silently-missing
    // gate; and the promoted rows must read green in BOTH files.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("workspace root")
        .to_path_buf();
    let registry = std::fs::read_to_string(root.join("instrument-claims.json")).expect("registry");
    let summary = std::fs::read_to_string(root.join("data/claims/wind-gate-summary.tsv"))
        .expect("committed wind gate summary");
    assert!(summary.starts_with("# frankensim-wind-gate-summary-v1"));
    let mut images = Vec::new();
    let mut cursor = 0usize;
    while let Some(hit) = registry[cursor..].find("\"filling\": \"wind-reed\"") {
        let at = cursor + hit;
        let tail = &registry[at..];
        let tag = "\"image\": \"";
        let img_at = tail.find(tag).expect("image") + tag.len();
        let img_end = tail[img_at..].find('"').expect("end");
        images.push(tail[img_at..img_at + img_end].to_string());
        cursor = at + tag.len();
    }
    assert!(images.len() >= 5, "the wind menu: {images:?}");
    for image in &images {
        assert!(summary.contains(image.as_str()), "summary omits {image}");
    }
    for green in [
        "tmm\tfingering-peak-cents\tgreen",
        "char-line\tjunction-match-tmm-r\tgreen",
        "char-line\tquarter-wave-lock\tgreen",
        "vfit-hold\theld-fingering-cents\tgreen",
        "phs-chain\ttd-oracle-agreement\tungated",
    ] {
        assert!(summary.contains(green), "summary must carry: {green}");
    }
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"wind-gate-summary\",\"verdict\":\"pass\",\
         \"rows\":{}}}",
        images.len()
    );
}
