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
        porcine.scalars["young_modulus_dynamic_day_printed_range_high"].to_bits()
            != porcine.scalars["young_modulus_dynamic_printed_range_high"].to_bits()
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

// ---------------------------------------------------------------------------------------------------------------
// HUMAN layered vocal-fold packs of music bead `frankensim-irp-v5-instrument-realism-pexfp.194` (M8-W2d): the layered
// human lamina propria of Kelleher et al. 2013 and the human elastic set of Alipour and Vigmostad 2012. Both are
// publisher-copyright papers, so the packs hold numeric facts with attribution (licence token
// `Published-Scientific-Facts-Attribution-Retained`), never a table layout, and the loader below does not require
// CC-BY. The gates encode what the sources themselves state: every typed value equals the printed one; the
// per-specimen ratios recomputed from the typed rows reproduce the means the source prints; the layer ordering the
// authors report; E_L > E_T > G_L in every specimen; every observation binds the specimen state (age, sex, hours
// postmortem, storage) to its values; and the two independent human sources agree on the longitudinal low-strain
// scale within a factor 2. The older packs above stay untouched; their statement that no pack holds a measured
// density is superseded for measured use by `tf-rho` / `ff-rho` of the Alipour pack (a facts-only row).
// ---------------------------------------------------------------------------------------------------------------

struct HumanPack {
    manifest: String,
    /// claim id -> (value, unit, observation id)
    scalar: BTreeMap<String, (f64, String, String)>,
    /// claim id -> uncertainty kind and, for absolute kinds, the spread
    uncertainty: BTreeMap<String, (String, Option<f64>)>,
    /// claim id -> number of validity rows
    validity: BTreeMap<String, usize>,
    /// observation id -> specimen, method and caveats joined
    observation: BTreeMap<String, String>,
}

impl HumanPack {
    fn value(&self, claim: &str) -> f64 {
        self.scalar
            .get(claim)
            .unwrap_or_else(|| panic!("claim {claim} missing"))
            .0
    }
}

fn load_human_pack(dir: &str) -> HumanPack {
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(
        manifest.contains("Published-Scientific-Facts-Attribution-Retained"),
        "{dir}: facts-only licence token"
    );
    assert!(manifest.contains("https://doi.org/"), "{dir}: DOI required");
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let mut pack = HumanPack {
        manifest,
        scalar: BTreeMap::new(),
        uncertainty: BTreeMap::new(),
        validity: BTreeMap::new(),
        observation: BTreeMap::new(),
    };
    for line in text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "scalar" => {
                assert_eq!(cols.len(), 7, "{dir}: scalar arity");
                let value: f64 = cols[4].parse().expect("scalar value");
                assert!(
                    value.is_finite() && value > 0.0,
                    "{dir}: {} positive",
                    cols[1]
                );
                let previous = pack.scalar.insert(
                    cols[1].to_string(),
                    (value, cols[5].to_string(), cols[2].to_string()),
                );
                assert!(previous.is_none(), "{dir}: duplicate claim {}", cols[1]);
            }
            "uncertainty" => {
                assert_eq!(cols.len(), 7, "{dir}: uncertainty arity");
                let spread = match cols[2] {
                    "unstated" => None,
                    "absolute" => Some(cols[3].parse::<f64>().expect("spread")),
                    other => panic!("{dir}: unexpected uncertainty kind {other}"),
                };
                pack.uncertainty
                    .insert(cols[1].to_string(), (cols[2].to_string(), spread));
            }
            "validity" => {
                assert_eq!(cols.len(), 6, "{dir}: validity arity");
                *pack.validity.entry(cols[1].to_string()).or_insert(0) += 1;
            }
            "observation" => {
                assert_eq!(cols.len(), 5, "{dir}: observation arity");
                pack.observation.insert(
                    cols[1].to_string(),
                    format!("{} {} {}", cols[2], cols[3], cols[4]),
                );
            }
            other => panic!("unknown row kind {other}"),
        }
    }
    for (id, (_, _, obs)) in &pack.scalar {
        assert!(
            pack.uncertainty.contains_key(id),
            "{dir}: {id} lacks an uncertainty row"
        );
        assert!(
            pack.validity.contains_key(id),
            "{dir}: {id} lacks a validity row"
        );
        assert!(
            pack.observation.contains_key(obs),
            "{dir}: {id} names missing observation {obs}"
        );
    }
    pack
}

/// Kelleher et al. 2013 Table 2: (subject, layer, E_L, E_T, G_L [kPa], R2 of the contact fit).
const TABLE_2: [(&str, &str, f64, f64, f64, f64); 12] = [
    ("a", "ligament", 28.20, 2.53, 0.91, 1.000),
    ("a", "cover", 196.40, 5.18, 4.18, 0.997),
    ("b", "ligament", 29.95, 1.18, 0.89, 0.997),
    ("b", "cover", 579.20, 7.17, 6.01, 0.978),
    ("c", "ligament", 24.32, 1.48, 0.77, 0.997),
    ("c", "cover", 268.20, 8.25, 7.50, 0.976),
    ("d", "ligament", 25.82, 1.51, 0.78, 0.996),
    ("d", "cover", 302.25, 14.35, 6.42, 0.999),
    ("e", "ligament", 38.56, 5.06, 1.59, 0.999),
    ("e", "cover", 199.80, 4.80, 3.24, 0.976),
    ("f", "ligament", 14.11, 3.08, 1.03, 0.999),
    ("f", "cover", 253.85, 11.93, 7.62, 0.988),
];

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * b.abs().max(1.0)
}

#[test]
fn g2_human_layer_pack_pins_the_printed_table_and_the_printed_means() {
    let pack = load_human_pack("vocalfold-human-layers-kelleher-bmmb-pmc3745772");

    // Every typed value equals the printed one, with the specimen state bound to it.
    for (subject, layer, el, et, gl, r2) in TABLE_2 {
        let id = format!("{subject}-{layer}");
        assert!(near(pack.value(&format!("{id}-el")), el), "{id} E_L");
        assert!(near(pack.value(&format!("{id}-et")), et), "{id} E_T");
        assert!(near(pack.value(&format!("{id}-gl")), gl), "{id} G_L");
        assert!(near(pack.value(&format!("{id}-r2")), r2), "{id} R2");
        assert_eq!(pack.scalar[&format!("{id}-el")].1, "kPa");
        let state = &pack.observation[&id];
        assert!(
            state.contains("h postmortem") && state.contains("PBS") && state.contains(" y, "),
            "{id}: the observation must carry age, postmortem hours and storage"
        );
        assert!(
            state.contains("1 Hz") && state.contains("coupled"),
            "{id}: loading rate and the E_T/G_L coupling caveat must be carried"
        );
        // Per-value uncertainty is not printed by the source: it must say so.
        for part in ["el", "et", "gl", "r2"] {
            assert_eq!(pack.uncertainty[&format!("{id}-{part}")].0, "unstated");
        }
    }

    // The per-specimen ratios recomputed from the typed rows reproduce the means the source prints
    // (14 and 39 for E_L/E_T, 28 and 54 for E_L/G_L, rounded to integers in the text).
    fn mean_ratio(layer: &str, f: fn(f64, f64, f64) -> f64) -> f64 {
        let rows: Vec<f64> = TABLE_2
            .iter()
            .filter(|r| r.1 == layer)
            .map(|r| f(r.2, r.3, r.4))
            .collect();
        assert_eq!(rows.len(), 6);
        rows.iter().sum::<f64>() / 6.0
    }
    fn el_over_et(el: f64, et: f64, _gl: f64) -> f64 {
        el / et
    }
    fn el_over_gl(el: f64, _et: f64, gl: f64) -> f64 {
        el / gl
    }
    let cases: [(&str, &str, f64, f64, fn(f64, f64, f64) -> f64); 4] = [
        ("ligament", "e-l-over-e-t", 14.0, 8.0, el_over_et),
        ("cover", "e-l-over-e-t", 39.0, 22.0, el_over_et),
        ("ligament", "e-l-over-g-l", 28.0, 8.0, el_over_gl),
        ("cover", "e-l-over-g-l", 54.0, 23.0, el_over_gl),
    ];
    for (layer, what, printed, spread, f) in cases {
        let id = format!("ratio-{what}-{layer}");
        assert!(near(pack.value(&id), printed), "{id} typed mean");
        assert_eq!(pack.uncertainty[&id].1, Some(spread), "{id} typed spread");
        let recomputed = mean_ratio(layer, f);
        assert!(
            (recomputed - printed).abs() < 0.5,
            "{id}: recomputed mean {recomputed} against printed {printed}"
        );
    }

    // The layer ordering the authors report.
    let mut stiffer_ratio = Vec::new();
    for subject in ["a", "b", "c", "d", "e", "f"] {
        let lig = |part: &str| pack.value(&format!("{subject}-ligament-{part}"));
        let cov = |part: &str| pack.value(&format!("{subject}-cover-{part}"));
        let ratio = cov("el") / lig("el");
        assert!(
            (3.0..30.0).contains(&ratio),
            "{subject}: cover/ligament E_L {ratio}"
        );
        stiffer_ratio.push(ratio);
        for layer in ["ligament", "cover"] {
            let (el, et, gl) = (
                pack.value(&format!("{subject}-{layer}-el")),
                pack.value(&format!("{subject}-{layer}-et")),
                pack.value(&format!("{subject}-{layer}-gl")),
            );
            assert!(el > et && et > gl, "{subject} {layer}: E_L > E_T > G_L");
        }
        // Cover E_T and G_L exceed the ligament's in every subject except E, whose cover E_T is lower
        // (G_L of subject E is still higher): the exception the source names.
        assert!(cov("gl") > lig("gl"), "{subject}: cover G_L above ligament");
        if subject == "e" {
            assert!(cov("et") < lig("et"), "subject E: the named exception");
        } else {
            assert!(cov("et") > lig("et"), "{subject}: cover E_T above ligament");
        }
    }
    let geometric = (stiffer_ratio.iter().map(|r| r.ln()).sum::<f64>() / 6.0).exp();
    assert!(
        (5.0..20.0).contains(&geometric),
        "cover is stiffer than ligament in tension by roughly a factor 10 (geometric mean {geometric})"
    );
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"g2-human-layers\",\"verdict\":\"pass\",\
         \"cover_over_ligament_el_geometric_mean\":{geometric},\"packs\":1}}"
    );
}

#[test]
fn g2_human_elastic_set_pins_the_table_and_agrees_with_the_layered_pack() {
    let set = load_human_pack("vocalfold-human-elastic-set-alipour-jvoice-pmc3508138");
    let layers = load_human_pack("vocalfold-human-layers-kelleher-bmmb-pmc3745772");

    // Table 1: mean and spread, in kPa, kg/m3 (g/mL x 1000) and dimensionless.
    for (id, mean, spread, unit) in [
        ("tf-el", 30.0, 13.4, "kPa"),
        ("tf-et", 1.0, 0.5, "kPa"),
        ("ff-el", 13.9, 7.6, "kPa"),
        ("tf-gl", 0.67, 0.31, "kPa"),
        ("tf-nu", 0.57, 0.17, "1"),
        ("tf-rho", 1029.9, 5.1, "kg/m3"),
        ("ff-rho", 1019.3, 4.6, "kg/m3"),
    ] {
        assert!(near(set.value(id), mean), "{id} mean");
        assert_eq!(set.uncertainty[id].1, Some(spread), "{id} spread");
        assert_eq!(set.scalar[id].1, unit, "{id} unit");
    }

    // A measured density of human true and false fold tissue; true above false (the source's paired t-test
    // finding), both in the soft-tissue class.
    let (rho_t, rho_f) = (set.value("tf-rho"), set.value("ff-rho"));
    assert!(
        rho_t > rho_f && (1000.0..1100.0).contains(&rho_f) && (1000.0..1100.0).contains(&rho_t)
    );
    assert!(set.observation["tissue-density"].contains("MEASURED"));
    assert!(set.observation["tissue-density"].contains("Archimedes"));

    // Transverse stiffness far below longitudinal (the source: about 30 times), Poisson ratio near 0.5.
    let anisotropy = set.value("tf-el") / set.value("tf-et");
    assert!((20.0..40.0).contains(&anisotropy), "E_L/E_T {anisotropy}");
    assert!((0.4..0.7).contains(&set.value("tf-nu")));
    // The nonlinearity statement is a printed range, and the typical sample sits inside it.
    assert!(set.value("stiffening-low") < set.value("stiffening-high"));
    let typical = set.value("typical-high") / set.value("typical-low");
    assert!(
        (set.value("stiffening-low")..=set.value("stiffening-high")).contains(&typical),
        "typical sample slope ratio {typical} inside the printed 10-15"
    );

    // State is bound to the values: donors, storage, rate, and the register of the modulus.
    for id in [
        "tf-longitudinal",
        "tf-transverse",
        "ff-longitudinal",
        "tf-shear",
        "tf-poisson",
    ] {
        let state = &set.observation[id];
        assert!(
            state.contains("Geriatric") && state.contains("1 Hz") && state.contains("nonlinear"),
            "{id}: donors, rate and nonlinearity must be carried"
        );
    }

    // Two independent human sources, two methods (ergometer on lamina propria strips; tension with optical strain
    // on separated ligament): the longitudinal low-strain scale agrees within a factor 2.
    let ligament: Vec<f64> = TABLE_2
        .iter()
        .filter(|r| r.1 == "ligament")
        .map(|r| r.2)
        .collect();
    let ligament_mean = ligament.iter().sum::<f64>() / 6.0;
    let ratio = set.value("tf-el") / ligament_mean;
    assert!(
        (0.5..2.0).contains(&ratio),
        "Alipour {} vs ligament mean {ligament_mean}",
        set.value("tf-el")
    );
    assert!(near(layers.value("a-ligament-el"), 28.2));

    // The measured density sits within 2 percent of the stated 1040 kg/m3 constant of the older register pack.
    assert!((rho_t - 1040.0).abs() / 1040.0 < 0.02);
    assert!(set.manifest.contains("Alipour"));
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"g2-human-elastic-set\",\"verdict\":\"pass\",\
         \"rho_true_fold\":{rho_t},\"rho_false_fold\":{rho_f},\"el_over_et\":{anisotropy}}}"
    );
}
