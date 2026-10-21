//! Physicality gates for the instrument-string matdb packs (bead `frankensim-irp-v5-instrument-realism-pexfp.192`, plan item M8-W2b):
//! nylon, fluorocarbon and natural gut monofilament, from Woodhouse and Lynch-Aird 2019 (CC BY 4.0).
//!
//! The gates encode what cannot be wrong whichever string was measured: positivity, the paper's own statement that fluorocarbon is
//! denser than gut or nylon, the loss factor eta_E = E''/E' the paper quotes for each material (about 0.04 at low stress for all three,
//! falling to about 0.02 for nylon over the tested stress range), the fit's stress axis giving physical stresses at those loss factors,
//! and schema hygiene (every scalar carries an uncertainty row and at least one validity row).

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

struct Pack {
    scalars: BTreeMap<String, f64>,
    uncertainty_rows: Vec<String>,
    validity_rows: Vec<String>,
}

fn load_pack(dir: &str) -> Pack {
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(
        manifest.contains("CC-BY-4.0"),
        "{dir}: licensing-first CC BY"
    );
    assert!(
        manifest.contains("https://doi.org/10.3813/AAA.919333"),
        "{dir}: DOI required"
    );
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let mut pack = Pack {
        scalars: BTreeMap::new(),
        uncertainty_rows: Vec::new(),
        validity_rows: Vec::new(),
    };
    let mut scalar_ids = Vec::new();
    for line in text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "scalar" => {
                let value: f64 = cols[4].parse().expect("scalar value");
                pack.scalars.insert(cols[3].to_string(), value);
                scalar_ids.push(cols[1].to_string());
            }
            "uncertainty" => pack.uncertainty_rows.push(cols[1].to_string()),
            "validity" => pack.validity_rows.push(cols[1].to_string()),
            "observation" => {}
            other => panic!("unknown row kind {other}"),
        }
    }
    for id in &scalar_ids {
        assert!(
            pack.uncertainty_rows.iter().any(|u| u == id),
            "{dir}: scalar {id} lacks an uncertainty row"
        );
        assert!(
            pack.validity_rows.iter().any(|v| v == id),
            "{dir}: scalar {id} lacks a validity row"
        );
    }
    for (name, &value) in &pack.scalars {
        assert!(
            value.is_finite() && value > 0.0,
            "{dir}: {name} must be positive and finite"
        );
    }
    pack
}

/// Stress in GPa at which E''/E' reaches `eta`, for E' = a + b sigma (a in Pa, b dimensionless with sigma in GPa) and constant E''.
fn stress_gpa_at_loss_factor(p: &Pack, eta: f64) -> f64 {
    let a_gpa = p.scalars["young_modulus_storage_intercept"] / 1e9;
    let b = p.scalars["young_modulus_storage_per_gpa_stress"];
    let loss_gpa = p.scalars["young_modulus_loss"] / 1e9;
    (loss_gpa / eta - a_gpa) / b
}

#[test]
fn string_packs_pass_physicality_gates() {
    let nylon = load_pack("string-nylon-woodhouse-lynch-aird-2019");
    let fluoro = load_pack("string-fluorocarbon-woodhouse-lynch-aird-2019");
    let gut = load_pack("string-gut-woodhouse-lynch-aird-2019");

    // The paper: "fluorocarbon [has] a significantly higher density than natural gut or nylon".
    let (rn, rf, rg) = (
        nylon.scalars["density"],
        fluoro.scalars["density"],
        gut.scalars["density"],
    );
    assert!(rf > rg && rf > rn, "fluorocarbon must be the densest");
    assert!(
        rf / rg.max(rn) > 1.3,
        "and by a wide margin: {rf} vs {rg}, {rn}"
    );

    // Gut: E' and E'' fixed, and the quoted eta_E = 0.04 must be their ratio within the table's rounding.
    let eta_gut = gut.scalars["young_modulus_loss"] / gut.scalars["young_modulus_storage"];
    assert!(
        (eta_gut - gut.scalars["loss_factor"]).abs() < 5e-4,
        "gut eta_E {eta_gut} vs quoted {}",
        gut.scalars["loss_factor"]
    );

    // Nylon and fluorocarbon: the quoted low-stress eta_E = 0.04 is reached at a positive, physical stress of the straight-line fit.
    for (name, p) in [("nylon", &nylon), ("fluorocarbon", &fluoro)] {
        let s = stress_gpa_at_loss_factor(p, p.scalars["loss_factor_at_low_stress"]);
        assert!(
            (0.005..0.15).contains(&s),
            "{name}: eta_E = 0.04 at {s} GPa is outside a plausible low-stress band"
        );
    }

    // Nylon: eta_E falls from 0.04 to 0.02 "as stress increases over the range tested"; the fit must reach 0.02 above the 0.04 stress
    // and below a typical nylon breaking stress (0.5 GPa), and E' must rise with stress.
    let s_lo = stress_gpa_at_loss_factor(&nylon, nylon.scalars["loss_factor_at_low_stress"]);
    let s_hi = stress_gpa_at_loss_factor(
        &nylon,
        nylon.scalars["loss_factor_at_highest_tested_stress"],
    );
    assert!(
        s_hi > s_lo && s_hi < 0.5,
        "nylon stress span {s_lo}..{s_hi} GPa"
    );
    assert!(nylon.scalars["young_modulus_storage_per_gpa_stress"] > 0.0);

    // The kHz storage modulus of all three lies in the polymer band 3-20 GPa at a harp-string stress of 0.1 GPa (gut has none).
    let at = |p: &Pack, sigma_gpa: f64| {
        p.scalars["young_modulus_storage_intercept"]
            + p.scalars["young_modulus_storage_per_gpa_stress"] * sigma_gpa * 1e9
    };
    for (name, e) in [
        ("nylon", at(&nylon, 0.1)),
        ("fluorocarbon", at(&fluoro, 0.1)),
        ("gut", gut.scalars["young_modulus_storage"]),
    ] {
        assert!(
            (3e9..20e9).contains(&e),
            "{name}: E' {e} Pa outside the polymer band"
        );
    }
}
