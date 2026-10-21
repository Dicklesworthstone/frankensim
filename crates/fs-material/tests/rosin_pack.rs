//! Physicality gates for the violin-rosin matdb packs (bead `frankensim-irp-v5-instrument-realism-pexfp.192`, plan item M8-W2b, "rosin
//! (third-body data)"): bulk properties and the viscosity law of Smith and Woodhouse 2000, and two independent glass-transition values
//! (Woodhouse and Galluzzo 2025, CC BY 4.0; Woodhouse, Euphonics 9.6).
//!
//! The gates encode what cannot be wrong whichever rosin was measured: positivity, the printed "approximately four orders of magnitude"
//! of viscosity between 45 and 85 C against the printed exponential law, the qualitative softening temperatures lying inside the range the
//! law was fitted over and in order, the two sources agreeing on the violin rosin's transition (about 50 C and 49 C), the double-bass
//! rosin's transition lying below room temperature (the page says it is already above it at room temperature) and the violin rosin's
//! above it, and schema hygiene (every scalar has an uncertainty and a validity row).

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

fn load_pack(dir: &str, license: &str) -> BTreeMap<String, f64> {
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(
        manifest.contains(&format!("\nlicense\t{license}\n")),
        "{dir}: license {license}"
    );
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let (mut scalars, mut ids) = (BTreeMap::new(), Vec::new());
    let (mut unc, mut val) = (Vec::new(), Vec::new());
    for line in text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "scalar" => {
                let value: f64 = cols[4].parse().expect("scalar value");
                assert!(
                    value.is_finite() && value > 0.0,
                    "{dir}: {} must be positive",
                    cols[3]
                );
                scalars.insert(cols[3].to_string(), value);
                ids.push(cols[1].to_string());
            }
            "uncertainty" => unc.push(cols[1].to_string()),
            "validity" => val.push(cols[1].to_string()),
            "observation" => {}
            other => panic!("unknown row kind {other}"),
        }
    }
    for id in &ids {
        assert!(
            unc.contains(id),
            "{dir}: scalar {id} lacks an uncertainty row"
        );
        assert!(val.contains(id), "{dir}: scalar {id} lacks a validity row");
    }
    scalars
}

/// The printed Eq. (3): eta = prefactor x 10^((reference - T)/interval) Pa s, T in kelvin above the ambient temperature.
fn viscosity(p: &BTreeMap<String, f64>, kelvin: f64) -> f64 {
    let excess = kelvin - p["viscosity_law_ambient_temperature"];
    p["viscosity_law_prefactor"]
        * 10f64.powf(
            (p["viscosity_law_reference_excess_temperature"] - excess)
                / p["viscosity_law_decade_temperature_interval"],
        )
}

#[test]
fn rosin_packs_pass_physicality_gates() {
    let sw = load_pack(
        "rosin-smith-woodhouse-jmps2000",
        "Published-Scientific-Facts-Attribution-Retained",
    );
    let wg = load_pack(
        "rosin-violin-glass-transition-woodhouse-galluzzo-2025",
        "CC-BY-4.0",
    );
    let eu = load_pack(
        "rosin-glass-transition-woodhouse-euphonics-9-6",
        "Published-Scientific-Facts-Attribution-Retained",
    );

    // Bulk density of a resin acid solid: close to water, a little above.
    assert!(
        (1000.0..1100.0).contains(&sw["density"]),
        "density {}",
        sw["density"]
    );

    // The printed law gives the printed "approximately four orders of magnitude" between 45 and 85 C, and its stated anchors.
    let (t45, t85) = (318.15, 358.15);
    let decades = (viscosity(&sw, t45) / viscosity(&sw, t85)).log10();
    assert!(
        (decades - sw["viscosity_decades_lost_over_45_to_85_c_printed"]).abs() < 0.25,
        "{decades} decades"
    );
    let at70 = viscosity(&sw, 343.15);
    assert!((at70 - 3000.0).abs() < 1.0, "eta(70 C) = {at70}");

    // Softening statements are ordered and lie inside the fitted range.
    let (indent, sticky) = (
        sw["knife_indent_temperature_about"],
        sw["thick_sticky_liquid_temperature_about"],
    );
    assert!(indent < sticky);
    assert!(t45 <= indent && sticky <= t85, "{indent} {sticky}");

    // Two sources on the violin rosin's transition agree to a few kelvin; the bass rosin's lies below room temperature.
    let (tg_wg, tg_eu) = (
        wg["glass_transition_temperature_about"],
        eu["violin_rosin_glass_transition_temperature"],
    );
    assert!((tg_wg - tg_eu).abs() < 3.0, "{tg_wg} against {tg_eu}");
    let (room, tg_bass) = (293.15, eu["double_bass_rosin_glass_transition_temperature"]);
    assert!(
        tg_bass < room && room < tg_eu,
        "bass {tg_bass}, room {room}, violin {tg_eu}"
    );
    // The viscosity law's fitted range straddles the violin transition (Woodhouse and Galluzzo: properties vary rapidly across it).
    assert!(t45 < tg_eu && tg_eu < t85);
}
