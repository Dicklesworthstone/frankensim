//! Physicality gates for the harp-string bulk-density matdb packs (bead `frankensim-irp-v5-instrument-realism-pexfp.192`, plan item
//! M8-W2b): per-string unstretched diameter and bulk density of Bowbrand nylon and gut and Savarez fluorocarbon harp strings (Lynch-Aird
//! and Woodhouse 2018, CC BY 4.0) and the gut twist/finish study's ranges (2023, CC BY 4.0).
//!
//! The gates encode what cannot be wrong whichever string was weighed: positivity; the paper's own statement that the nylon strings split
//! into a lower- and a higher-density group (the thicker strings in the higher group) and that fluorocarbon is the densest and wound
//! fluorocarbon strings less dense than monofilament ones; the representative densities of Woodhouse and Lynch-Aird 2019 (the
//! `string-*-woodhouse-lynch-aird-2019` packs) lying within 3 % of the mean of the thinner strings of the same material here;
//! the 2023 gut density range containing every 2018 gut density; schema hygiene (every scalar has uncertainty and validity rows).

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

fn load_pack(dir: &str) -> BTreeMap<String, f64> {
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(manifest.contains("\nlicense\tCC-BY-4.0\n"), "{dir}: licensing-first CC BY");
    assert!(manifest.contains("https://doi.org/10."), "{dir}: DOI required");
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
                assert!(value.is_finite() && value > 0.0, "{dir}: {} must be positive", cols[3]);
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
        assert!(unc.contains(id), "{dir}: scalar {id} lacks an uncertainty row");
        assert!(val.contains(id), "{dir}: scalar {id} lacks a validity row");
    }
    scalars
}

#[test]
fn harp_string_density_packs_pass_physicality_gates() {
    let h = load_pack("string-harp-bulk-density-lynch-aird-woodhouse-2018");
    let g23 = load_pack("string-gut-harp-twist-lynch-aird-woodhouse-2023");
    let rho = |id: &str| h[&format!("bulk_density_{id}")];
    let dia = |id: &str| h[&format!("unstretched_diameter_{id}")];

    // Nylon: two density groups, the thicker strings in the higher group.
    let low = ["n8", "n14", "n5"];
    let high = ["n23a", "n29"];
    let low_max = low.iter().map(|i| rho(i)).fold(0.0, f64::max);
    let high_min = high.iter().map(|i| rho(i)).fold(f64::INFINITY, f64::min);
    assert!(high_min > low_max, "nylon density groups overlap: {low_max} vs {high_min}");
    let low_dia_max = low.iter().map(|i| dia(i)).fold(0.0, f64::max);
    let high_dia_min = high.iter().map(|i| dia(i)).fold(f64::INFINITY, f64::min);
    assert!(high_dia_min > low_dia_max, "the higher-density group is the thicker strings");

    // Fluorocarbon: monofilaments (C1, C3) denser than the wound strings (C5b, C7a); fluorocarbon denser than gut denser than nylon.
    assert!(rho("c1").min(rho("c3")) > rho("c5b").max(rho("c7a")));
    let gut = ["g2b", "g3b", "g3c", "g5b", "g5c"];
    let gut_min = gut.iter().map(|i| rho(i)).fold(f64::INFINITY, f64::min);
    let gut_max = gut.iter().map(|i| rho(i)).fold(0.0, f64::max);
    assert!(rho("c5b").min(rho("c7a")) > gut_max, "even the wound fluorocarbon exceeds gut");
    assert!(gut_min > high.iter().map(|i| rho(i)).fold(0.0, f64::max), "gut exceeds nylon");

    // Cross-pack: the 2019 representative densities against the thinner strings here.
    let p19 = |m: &str| load_pack(&format!("string-{m}-woodhouse-lynch-aird-2019"))["density"];
    let mean = |ids: &[&str]| ids.iter().map(|i| rho(i)).sum::<f64>() / ids.len() as f64;
    for (name, ids) in [("nylon", &low[..]), ("fluorocarbon", &["c1", "c3"][..]), ("gut", &gut[..])] {
        let (a, b) = (p19(name), mean(ids));
        assert!((a - b).abs() / b < 0.03, "{name}: 2019 density {a} vs thinner-string mean {b}");
    }

    // 2023: its density range contains every gut density of 2018; stretching leaves 97-98 % of the density.
    assert!(g23["bulk_density_unstretched_table_range_low"] <= gut_min);
    assert!(g23["bulk_density_unstretched_table_range_high"] >= gut_max);
    assert!(g23["stretched_density_fraction_of_unstretched_low"] < g23["stretched_density_fraction_of_unstretched_high"]);
    assert!(g23["stretched_density_fraction_of_unstretched_high"] < 1.0);
    assert!(g23["breaking_strain_typical_low"] < g23["breaking_strain_typical_high"]);
    assert!(g23["twist_angle_degrees_table_range_low"] < g23["twist_angle_degrees_table_range_high"]);
    assert!(g23["unstretched_diameter_table_range_low"] < g23["unstretched_diameter_table_range_high"]);
}
