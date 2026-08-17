//! g2 physicality gates for the vocal-fold / lip tissue matdb packs
//! (music bead `frankensim-music-v8-root-3ez8g.3.7`) — the tissue card
//! population the reduce lab and the glottis islands read.
//!
//! The gates encode what CANNOT be wrong: positivity + schema hygiene;
//! the human vertical stiffness gradient (superior < medial <
//! inferior, the source's own finding); METHOD-SCALE DISJOINTNESS
//! (Pa-scale felid shear matrix vs kPa-scale human indentation vs
//! tens-of-kPa porcine dynamic — different methods measure different
//! things, recorded never averaged, and every observation must say
//! why); the loss tangents in the soft-tissue class; the MODEL-CARD
//! firewall (two-mass parameters must declare themselves NOT tissue
//! properties); and the stated-not-measured density register.

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
    units: BTreeMap<String, String>,
    observation_notes: String,
}

fn load_pack(dir: &str) -> Pack {
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(manifest.contains("CC-BY-4.0"), "{dir}: licensing-first");
    assert!(manifest.contains("https://doi.org/"), "{dir}: DOI required");
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let mut pack = Pack {
        scalars: BTreeMap::new(),
        units: BTreeMap::new(),
        observation_notes: String::new(),
    };
    let mut uncertainty_ids = Vec::new();
    let mut validity_ids = Vec::new();
    let mut scalar_ids = Vec::new();
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
                pack.scalars.insert(cols[3].to_string(), value);
                pack.units.insert(cols[3].to_string(), cols[5].to_string());
                scalar_ids.push(cols[1].to_string());
            }
            "uncertainty" => uncertainty_ids.push(cols[1].to_string()),
            "validity" => validity_ids.push(cols[1].to_string()),
            "observation" => pack.observation_notes.push_str(cols[4]),
            other => panic!("unknown row kind {other}"),
        }
    }
    for id in &scalar_ids {
        assert!(
            uncertainty_ids.iter().any(|u| u == id),
            "{dir}: scalar {id} lacks an uncertainty row"
        );
        assert!(
            validity_ids.iter().any(|v| v == id),
            "{dir}: scalar {id} lacks a validity row"
        );
    }
    pack
}

#[test]
fn g2_tissue_packs_pass_physicality_and_register_gates() {
    let human = load_pack("vocalfold-human-indentation-jvoice-pmc12180296");
    let porcine = load_pack("vocalfold-porcine-aspiration-mdpi-s21092923");
    let felid = load_pack("vocalfold-felid-rheometry-plos-pone0027029");
    let lip = load_pack("lip-perioral-invivo-mdpi-ma17153654");
    let model = load_pack("twomass-if72-standard-plos-pone0187486");
    let density = load_pack("tissue-density-stated-plos-pcbi1004907");

    // The human vertical stiffness gradient — the source's own finding.
    let sup = human.scalars["young_modulus_effective_superior"];
    let med = human.scalars["young_modulus_effective_medial"];
    let inf = human.scalars["young_modulus_effective_inferior"];
    assert!(
        sup < med && med < inf,
        "the vertical gradient must rise superior->inferior ({sup} < {med} < {inf} kPa)"
    );

    // METHOD-SCALE DISJOINTNESS: felid shear-matrix (Pa) sits orders
    // below human indentation (kPa) which sits below porcine dynamic
    // top (kPa) — different methods, and every observation must NAME
    // the method dependence so nobody averages across registers.
    let g_hi = felid.scalars["shear_storage_modulus_printed_range_high"];
    assert_eq!(
        felid.units["shear_storage_modulus_printed_range_high"],
        "Pa"
    );
    assert_eq!(human.units["young_modulus_effective_medial"], "kPa");
    assert!(
        g_hi / 1000.0 < sup,
        "felid shear (Pa scale) must sit below human indentation (kPa scale)"
    );
    let e_dyn_hi = porcine.scalars["young_modulus_dynamic_printed_range_high"];
    assert!(
        e_dyn_hi > inf,
        "porcine dynamic top above quasi-static human"
    );
    for (pack, name) in [(&human, "human"), (&porcine, "porcine"), (&felid, "felid")] {
        assert!(
            pack.observation_notes.contains("method")
                || pack.observation_notes.contains("never averaged"),
            "{name}: the observation must name the method dependence"
        );
    }

    // Postmortem time is load-bearing: the 1-day porcine range must
    // differ from the 4-hour range (recorded separately).
    assert!(
        porcine.scalars["young_modulus_dynamic_day_printed_range_high"]
            != porcine.scalars["young_modulus_dynamic_printed_range_high"]
    );

    // Loss tangents in the soft-tissue class, tiger above lion (the
    // source's roar finding).
    let td_tiger = felid.scalars["loss_tangent_tiger"];
    let td_lion = felid.scalars["loss_tangent_lion"];
    assert!((0.05..0.6).contains(&td_tiger) && (0.05..0.6).contains(&td_lion));
    assert!(td_tiger > td_lion);

    // The lip pack carries its structure-in-series caveat.
    assert!(
        lip.observation_notes.contains("teeth"),
        "the lip value's hard-teeth-in-series caveat is load-bearing"
    );
    assert!(
        lip.scalars["young_modulus_effective_upper_lip"]
            > lip.scalars["young_modulus_effective_cheek_left"]
    );

    // THE MODEL-CARD FIREWALL: the two-mass pack must declare itself
    // NOT tissue data, and its property names must be model_-prefixed
    // so a tissue consumer cannot pick them up by accident.
    assert!(
        model.observation_notes.contains("MODEL")
            && model.observation_notes.contains("not constitutive"),
        "the two-mass pack must declare its register"
    );
    for name in model.scalars.keys() {
        assert!(
            name.starts_with("model_"),
            "model-card property {name} must be model_-prefixed"
        );
    }
    // Sanity on the canonical set (S-H defaults).
    assert!((model.scalars["model_mass_lower"] - 0.125).abs() < 1e-12);
    assert!((model.scalars["model_stiffness_lower"] - 80.0).abs() < 1e-12);

    // The stated-not-measured density register: the value carries the
    // stated flag and the observation names the unlicensable measured
    // chain; no pack anywhere claims a MEASURED tissue density.
    assert!((density.scalars["density_stated_soft_tissue"] - 1040.0).abs() < 1e-12);
    assert!(
        density.observation_notes.contains("STATED")
            && density.observation_notes.contains("Perlman"),
        "the stated-density register must name the unlicensable measured chain"
    );
    for (pack, name) in [
        (&human, "human"),
        (&porcine, "porcine"),
        (&felid, "felid"),
        (&lip, "lip"),
    ] {
        assert!(
            !pack.scalars.keys().any(|k| k.starts_with("density")),
            "{name}: no measured tissue-density row may exist without a licensable source"
        );
    }
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"g2-tissue-packs\",\"verdict\":\"pass\",\
         \"human_gradient_kpa\":[{sup},{med},{inf}],\"tan_delta\":[{td_tiger},{td_lion}],\
         \"lip_kpa\":{},\"packs\":6}}",
        lip.scalars["young_modulus_effective_upper_lip"]
    );
}
