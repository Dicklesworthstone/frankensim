//! Physicality gates for the clarinet-reed cane packs of lane P12 (plan item M8-W2a): Obataya 1999 (J. Wood Sci. 45:106,
//! dry and soaked cane, L and LT moduli and loss), Obataya 1996 (Wood Res. Tech. Notes 32:30, moisture and extractive
//! effects) and Taillard et al. 2014 (Acta Acust. united Ac. 100:555, Zener fits of 55 Vandoren V12 reeds). All three are
//! numerical facts with attribution; no figure or table is copied.
//!
//! The gates encode what cannot be wrong whichever specimen was measured, plus the sources' own cross-checks: positivity,
//! soaking softens both moduli, the shear loss exceeds the fibre loss, the anisotropy E/G lies in the band the source
//! prints, Table 5 against Table 6 of the 1996 report, the 1996 and 1999 untreated air-dry moduli overlapping, the
//! moisture curve of untreated cane (loss peak near 7 percent moisture content, modulus falling to the soaked plateau),
//! and each Zener fit's printed 4 kHz storage modulus reproduced by its printed E1, E2 and relaxation frequency. Schema
//! hygiene: every scalar has an uncertainty and a validity row.

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

/// Claim id -> (value, absolute SD when the row states one).
fn load_pack(dir: &str, citation: &str) -> BTreeMap<String, (f64, Option<f64>)> {
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(manifest.contains("Published-Scientific-Facts-Attribution-Retained"));
    assert!(manifest.contains(citation), "{dir}: citation");
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let mut values = BTreeMap::new();
    let (mut unc, mut val) = (BTreeMap::new(), Vec::new());
    for line in text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "scalar" => {
                let value: f64 = cols[4].parse().expect("scalar value");
                assert!(
                    value.is_finite() && value > 0.0,
                    "{dir}: {} must be positive",
                    cols[1]
                );
                values.insert(cols[1].to_string(), value);
            }
            "uncertainty" => {
                let sd = (cols[2] == "absolute").then(|| cols[3].parse::<f64>().expect("sd"));
                unc.insert(cols[1].to_string(), sd);
            }
            "validity" => val.push(cols[1].to_string()),
            "observation" => {}
            other => panic!("unknown row kind {other}"),
        }
    }
    values
        .into_iter()
        .map(|(id, v)| {
            let sd = *unc
                .get(&id)
                .unwrap_or_else(|| panic!("{dir}: scalar {id} lacks an uncertainty row"));
            assert!(val.contains(&id), "{dir}: scalar {id} lacks a validity row");
            (id, (v, sd))
        })
        .collect()
}

#[test]
fn obataya_1999_dry_and_soaked_cane_pass_physicality_gates() {
    let p = load_pack(
        "cane-arundo-wetdry-obataya-jws1999",
        "Journal of Wood Science 45:106-112",
    );
    let v = |id: &str| p[id].0;
    for treat in ["untreated", "extracted"] {
        let (dry, wet) = (format!("{treat}-airdry"), format!("{treat}-wet"));
        // Soaking softens the fibre and the LT shear modulus (Tables 2 and 3).
        assert!(
            v(&format!("e-l-{wet}")) < v(&format!("e-l-{dry}")),
            "{treat}: E_L"
        );
        assert!(
            v(&format!("g-lt-{wet}")) < v(&format!("g-lt-{dry}")),
            "{treat}: G_LT"
        );
        // Soaked specimens (moisture content 189 and 232 percent, Table 1) are denser than water; air-dry ones are not.
        assert!(v(&format!("density-{wet}")) > 1000.0);
        assert!((300.0..600.0).contains(&v(&format!("density-{dry}"))));
        for state in [&dry, &wet] {
            // The printed tan delta_s / tan delta_L ratios are 1.81-2.70: shear dissipates more than the fibre.
            assert!(v(&format!("eta-lt-{state}")) > v(&format!("eta-l-{state}")));
            // The printed E/G is 8.39-10.10; ratios of the printed means sit within 0.1 of it.
            let eg = v(&format!("e-l-{state}")) / v(&format!("g-lt-{state}"));
            assert!((8.29..10.2).contains(&eg), "{state}: E/G {eg}");
            assert!(v(&format!("eta-l-{state}")) < 0.05);
        }
    }
    // The water-soluble extractives double the air-dry fibre loss (the paper's central cane finding).
    assert!(v("eta-l-untreated-airdry") > 1.5 * v("eta-l-extracted-airdry"));
}

#[test]
fn obataya_1996_moisture_and_extractive_effects_are_consistent() {
    let p = load_pack(
        "cane-arundo-moisture-obataya-wrtn1996",
        "Wood Research and Technical Notes 32:30-65",
    );
    let v = |id: &str| p[id].0;
    let sd = |id: &str| p[id].1.expect("sd");
    // Table 5 against Table 6 (same 60 specimens): the relative effects of the ratio of means lie within one SD of the
    // printed mean of the per-specimen effects.
    let de = v("e-l-t5-untreated") / v("e-l-t5-extracted") - 1.0;
    assert!(
        (de - v("t6-storage-modulus")).abs() <= sd("t6-storage-modulus"),
        "dE' {de}"
    );
    let dt = v("eta-l-t5-untreated") / v("eta-l-t5-extracted") - 1.0;
    assert!(
        (dt - v("t6-loss-factor")).abs() <= sd("t6-loss-factor"),
        "dtan {dt}"
    );
    assert!(v("sg-t5-untreated") > v("sg-t5-extracted"));
    // Same laboratory, same cane, same method: the 1996 Table 5 and the 1999 untreated air-dry moduli overlap.
    let q = load_pack(
        "cane-arundo-wetdry-obataya-jws1999",
        "Journal of Wood Science 45:106-112",
    );
    let e99 = q["e-l-untreated-airdry"].0;
    assert!((v("e-l-t5-untreated") - e99).abs() <= sd("e-l-t5-untreated"));
    // Figs 25 and 27 (read error 0.1 GPa and 0.0002): untreated modulus falls from oven-dry to the soaked plateau, its loss
    // peaks near 7 percent moisture content, and above the fibre saturation point both treatments agree within the reads.
    let read_e = 0.1e9;
    assert!(v("e-l-fig-untreated-mc0") > v("e-l-fig-untreated-mc7p2"));
    assert!(v("e-l-fig-untreated-mc7p2") > v("e-l-fig-untreated-mc19p0"));
    assert!((v("e-l-fig-untreated-mc19p0") - v("e-l-fig-untreated-soaked")).abs() <= read_e);
    assert!((v("e-l-fig-untreated-soaked") - v("e-l-fig-extracted-soaked")).abs() <= 2.0 * read_e);
    let peak = v("eta-l-fig-untreated-mc7p1");
    for id in [
        "eta-l-fig-untreated-mc0",
        "eta-l-fig-untreated-mc19p0",
        "eta-l-fig-untreated-soaked",
    ] {
        assert!(peak > v(id) + 0.0004, "{id}");
    }
    assert!(v("eta-l-fig-extracted-mc5p3") < v("eta-l-fig-extracted-mc0"));
}

#[test]
fn taillard_2014_zener_fits_reproduce_their_printed_storage_moduli() {
    let p = load_pack(
        "cane-arundo-vandoren-v12-zener-taillard-aaua2014",
        "Acta Acustica united with Acustica 100(3):555-573",
    );
    let v = |id: &str| p[id].0;
    let sd = |id: &str| p[id].1.expect("sd");
    for fit in ["h4-a", "h5-a", "h4-b", "h5-b"] {
        for (m, fr, storage) in [
            ("zener-e", "zener-fr-e", "e-l-storage-4khz"),
            ("zener-g", "zener-fr-g", "g-lt-storage-4khz"),
        ] {
            let (m1, m2) = (v(&format!("{m}1-{fit}")), v(&format!("{m}2-{fit}")));
            let r = (4000.0 / v(&format!("{fr}-{fit}"))).powi(2);
            let model = m1 + m2 * r / (1.0 + r);
            let printed = v(&format!("{storage}-{fit}"));
            // The printed value is the mean over 55 reeds of a nonlinear function of each reed's parameters, so it
            // differs from the function of the means by a second-order term; measured here at 0.07-0.13 SD.
            let spread = sd(&format!("{storage}-{fit}"));
            assert!(
                (model - printed).abs() <= 0.2 * spread,
                "{fit} {storage}: model {model} printed {printed}"
            );
            assert!(m1 < printed && printed < m1 + m2, "{fit} {storage}");
        }
    }
    // The authors' own summary: drying (series B) lowers E_L and raises G_LT.
    for h in ["h4", "h5"] {
        assert!(v(&format!("e-l-storage-4khz-{h}-b")) < v(&format!("e-l-storage-4khz-{h}-a")));
        assert!(v(&format!("g-lt-storage-4khz-{h}-b")) > v(&format!("g-lt-storage-4khz-{h}-a")));
    }
}
