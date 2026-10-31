//! Physicality gates for the horse bow-hair matdb pack (bead `frankensim-irp-v5-instrument-realism-pexfp.192`, plan item M8-W2b), from
//! Pitteroff and Woodhouse 1998, Acustica - acta acustica 84:543-562 (numerical facts with attribution; no figure or table is copied).
//!
//! The gates encode what cannot be wrong whichever hair was measured: positivity, the paper's own cross-check E_H = c_H^2 rho_H of about
//! 7 GPa, the dynamic modulus being the stiffer of the two measured at low frequency, the hair counts and ribbon dimensions ordering
//! violin < viola < cello < bass, and the mass per hair length (ribbon mass / (count x length)) of every bow lying inside the line density
//! of a hair of the printed density and printed diameter range. Schema hygiene: every scalar has an uncertainty and a validity row.

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
    assert!(manifest.contains("Published-Scientific-Facts-Attribution-Retained"));
    assert!(manifest.contains("Acustica - acta acustica 84:543-562"));
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

#[test]
fn bow_hair_pack_passes_physicality_gates() {
    let p = load_pack("bow-hair-horse-pitteroff-woodhouse-acustica1998");

    // The paper's second estimate: E_H = c_H^2 rho_H "of about 7 GPa", within the first estimate's printed +-0.5 GPa.
    let e_from_speed = p["longitudinal_wave_speed_quoted_from_schumacher"].powi(2) * p["density"];
    let e_static = p["young_modulus_static"];
    assert!(
        (e_from_speed - e_static).abs() <= p["young_modulus_static_printed_half_width"],
        "c^2 rho = {e_from_speed} against the static {e_static}"
    );

    // Dynamic (10-100 Hz) E_H is the stiffer one, and the two printed intervals touch.
    let e_dyn = p["young_modulus_dynamic"];
    assert!(e_dyn > e_static);
    assert!(
        e_dyn - p["young_modulus_dynamic_printed_half_width"]
            <= e_static + p["young_modulus_static_printed_half_width"]
    );
    assert!(p["loss_factor_printed_range_low"] < p["loss_factor_printed_range_high"]);
    assert!(
        p["exponential_decay_rate_printed_range_low"]
            < p["exponential_decay_rate_printed_range_high"]
    );
    assert!(p["hair_diameter_printed_range_low"] < p["hair_diameter_printed_range_high"]);
    assert!(p["hairs_in_contact_with_string_at_most"] <= p["violin_bow_hair_count_estimate"]);

    // Ribbon dimensions grow from violin to bass (Table I).
    let names = ["violin", "viola", "cello", "double_bass"];
    let get = |n: &str, what: &str| p[&format!("{n}_bow_{what}")];
    let mut prev = (0.0, 0.0, 0.0);
    for n in names {
        let count = get(n, "hair_count_estimate");
        let mass = get(n, "hair_ribbon_mass");
        let width = get(n, "ribbon_width");
        assert!(
            count > prev.0 && mass > prev.1 && width >= prev.2,
            "{n}: ordering"
        );
        prev = (count, mass, width);
    }

    // Mass per hair length of each bow lies in the line density of a hair of the printed density and printed diameter range
    // (the paper's ribbon mass is for a hank longer than the hair length, so the line density is an upper-leaning estimate).
    let line = |d: f64| p["density"] * std::f64::consts::PI * (d / 2.0).powi(2);
    let (lo, hi) = (
        line(p["hair_diameter_printed_range_low"]),
        line(p["hair_diameter_printed_range_high"]),
    );
    for n in names {
        let per_hair_per_m = get(n, "hair_ribbon_mass")
            / (get(n, "hair_count_estimate") * get(n, "hair_length_frog_to_tip"));
        assert!(
            (lo..=hi).contains(&per_hair_per_m),
            "{n}: {per_hair_per_m} kg/m per hair outside {lo}..{hi}"
        );
    }
}

#[test]
fn published_dashpot_reduction_retains_its_nominal_geometry_and_model_status() {
    let p = load_pack("bow-hair-horse-pitteroff-woodhouse-acustica1998");
    // P&W Part I pp. 546-548: the model uses 50 hairs of nominal diameter
    // 0.19 mm and length 0.65 m, driven at their midpoint. This conversion is
    // geometry normalization of the published D, not a fit to instrument audio.
    let factor = 50.0 * std::f64::consts::PI * 0.000095_f64.powi(2) / (0.65 * 0.5 * 0.5);
    let damping = factor * p["kelvin_voigt_young_viscosity_pw98_model"];
    assert!((damping - p["nominal_ribbon_damping_pw98_model"]).abs() < 1e-13);
    assert_eq!(p["nominal_ribbon_damping_pw98_model"], 10.0);
    let spring = factor * p["young_modulus_static"];
    assert!((spring / 60_000.0 - 1.0).abs() < 0.02, "the paper rounds S to 6e4 N/m");
    let text = std::fs::read_to_string(repo_root().join(
        "data/matdb/seed-v1/bow-hair-horse-pitteroff-woodhouse-acustica1998/properties.tsv")).unwrap();
    assert!(text.contains("MODEL DATA, NOT A MEASURED BROADBAND CONSTITUTIVE LAW"));
}
