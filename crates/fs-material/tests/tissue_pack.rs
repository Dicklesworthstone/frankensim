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
        lip.scalars["young_modulus_upper_lip_apparatus_effective"]
            > lip.scalars["young_modulus_cheek_left_apparatus_effective"]
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
        lip.scalars["young_modulus_upper_lip_apparatus_effective"]
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

/// Chhetri et al. 2011 Tables 1-3: (larynx, side, sample, location group, D mm, T mm (0 = not printed), h_max mm, E kPa, iqr kPa).
const INDENTATION: [(&str, &str, &str, &str, f64, f64, f64, f64, f64); 22] = [
    (
        "m49",
        "s1",
        "hemilarynx-medial-surface",
        "hemilarynx",
        0.5,
        0.0,
        0.8,
        7.4,
        1.6,
    ),
    (
        "m49",
        "s1",
        "cover-medial-surface",
        "medial-cover",
        0.5,
        2.0,
        0.5,
        5.4,
        2.0,
    ),
    (
        "m49",
        "s1",
        "cover-inferior-medial-surface",
        "inferior-medial-cover",
        0.5,
        3.0,
        0.75,
        7.9,
        4.6,
    ),
    (
        "m49",
        "s1",
        "ta-attached-to-cartilage",
        "ta-muscle",
        0.5,
        5.0,
        0.6,
        2.0,
        2.1,
    ),
    (
        "m49",
        "s1",
        "ta-attached-to-cartilage",
        "ta-muscle",
        2.0,
        5.0,
        0.8,
        1.5,
        0.8,
    ),
    (
        "m49",
        "s2",
        "hemilarynx-medial-surface",
        "hemilarynx",
        0.5,
        0.0,
        0.75,
        5.3,
        0.5,
    ),
    (
        "m49",
        "s2",
        "ta-attached-to-cartilage",
        "ta-muscle",
        0.5,
        5.0,
        0.75,
        1.3,
        2.4,
    ),
    (
        "m49",
        "s2",
        "ta-isolated",
        "ta-muscle",
        0.5,
        5.0,
        0.75,
        2.1,
        1.0,
    ),
    (
        "m49",
        "s2",
        "ta-isolated",
        "ta-muscle",
        1.0,
        5.0,
        0.75,
        1.6,
        1.0,
    ),
    (
        "f19",
        "s1",
        "hemilarynx",
        "hemilarynx",
        0.5,
        0.0,
        1.0,
        13.1,
        1.5,
    ),
    (
        "f19",
        "s1",
        "cover-superior-surface",
        "superior-cover",
        0.5,
        2.0,
        0.75,
        2.9,
        1.1,
    ),
    (
        "f19",
        "s1",
        "cover-medial-surface",
        "medial-cover",
        0.5,
        2.0,
        0.75,
        4.7,
        1.2,
    ),
    (
        "f19",
        "s2",
        "vocal-fold-medial-surface",
        "medial-cover",
        0.5,
        0.0,
        1.0,
        3.9,
        0.7,
    ),
    (
        "f19",
        "s2",
        "vocal-fold-inferior-medial-surface",
        "inferior-medial-cover",
        0.5,
        0.0,
        1.0,
        7.0,
        0.9,
    ),
    (
        "m30",
        "s1",
        "cover-superior-surface",
        "superior-cover",
        1.0,
        1.6,
        0.4,
        3.2,
        1.8,
    ),
    (
        "m30",
        "s1",
        "cover-medial-surface",
        "medial-cover",
        1.0,
        2.0,
        0.4,
        5.7,
        2.0,
    ),
    (
        "m30",
        "s1",
        "cover-medial-surface-lp-layer-up",
        "",
        1.0,
        2.0,
        0.4,
        2.7,
        2.3,
    ),
    (
        "m30",
        "s1",
        "ta-isolated",
        "ta-muscle",
        1.0,
        3.3,
        0.5,
        2.7,
        1.5,
    ),
    (
        "m30",
        "s2",
        "cover-superior-surface",
        "superior-cover",
        1.0,
        1.6,
        0.4,
        2.7,
        1.6,
    ),
    (
        "m30",
        "s2",
        "cover-medial-surface",
        "medial-cover",
        1.0,
        1.6,
        0.4,
        4.3,
        1.9,
    ),
    (
        "m30",
        "s2",
        "cover-medial-surface-lp-layer-up",
        "",
        1.0,
        1.6,
        0.4,
        4.4,
        1.9,
    ),
    (
        "m30",
        "s2",
        "ta-isolated",
        "ta-muscle",
        1.0,
        4.0,
        0.4,
        3.1,
        1.4,
    ),
];

#[test]
fn g2_human_indentation_pack_pins_the_rows_and_records_the_source_inconsistency() {
    let pack = load_human_pack("vocalfold-human-indentation-layers-chhetri-jvoice-pmc2891297");

    // Every typed row equals the printed one, with its indenter regime carried as validity.
    for (lx, side, sample, _group, d, t, h, e, iqr) in INDENTATION {
        let id = format!("{lx}-{side}-{sample}-d{d}");
        assert!(near(pack.value(&id), e), "{id} E");
        assert_eq!(pack.scalar[&id].1, "kPa");
        assert_eq!(pack.uncertainty[&id].1, Some(iqr), "{id} iqr");
        let expected_validity = if t > 0.0 { 5 } else { 3 };
        assert_eq!(
            pack.validity[&id], expected_validity,
            "{id}: h/D and D/T regime rows"
        );
        assert!(h > 0.0 && d > 0.0);
    }
    for lx in ["m49", "f19", "m30"] {
        let state = &pack.observation[&format!("apparatus-larynx-{lx}")];
        assert!(
            state.contains("postmortem") && state.contains("AN APPARATUS RESPONSE"),
            "{lx}: the apparatus-response register and the specimen state must be carried"
        );
    }

    // The printed location means (and N) are reproduced by the typed rows; the ranges too, except the one the source
    // gets wrong: the printed TA range 1.3-2.7 kPa does not contain the Table 3 value 3.1 kPa.
    let printed = [
        ("hemilarynx", 8.6, 5.3, 13.1, 3_usize),
        ("inferior-medial-cover", 7.5, 7.0, 7.9, 2),
        ("medial-cover", 4.8, 3.9, 5.7, 5),
        ("superior-cover", 2.9, 2.7, 3.2, 3),
        ("ta-muscle", 2.0, 1.3, 2.7, 7),
    ];
    for (group, mean_printed, lo, hi, n) in printed {
        let rows: Vec<f64> = INDENTATION
            .iter()
            .filter(|r| r.3 == group)
            .map(|r| r.7)
            .collect();
        assert_eq!(rows.len(), n, "{group}: N");
        let mean = rows.iter().sum::<f64>() / n as f64;
        assert!(
            (mean - mean_printed).abs() < 0.06,
            "{group}: recomputed {mean} vs printed {mean_printed}"
        );
        assert!(near(
            pack.value(&format!("summary-{group}-mean")),
            mean_printed
        ));
        assert!(near(pack.value(&format!("summary-{group}-low")), lo));
        assert!(near(pack.value(&format!("summary-{group}-high")), hi));
        let (row_lo, row_hi) = rows
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| {
                (a.min(v), b.max(v))
            });
        assert!(
            near(row_lo, lo),
            "{group}: printed low {lo} vs rows {row_lo}"
        );
        if group == "ta-muscle" {
            assert!(
                row_hi > hi,
                "the source's TA range excludes its own Table 3 value"
            );
            assert!(
                pack.observation["apparatus-printed-location-summary"].contains("does NOT contain")
            );
        } else {
            assert!(
                near(row_hi, hi),
                "{group}: printed high {hi} vs rows {row_hi}"
            );
        }
    }

    // The stiffness order the authors report: intact hemilarynx, then the isolated cover from inferior to superior,
    // then the isolated thyroarytenoid muscle.
    let mean = |g: &str| pack.value(&format!("summary-{g}-mean"));
    assert!(
        mean("hemilarynx") > mean("inferior-medial-cover")
            && mean("inferior-medial-cover") > mean("medial-cover")
            && mean("medial-cover") > mean("superior-cover")
            && mean("superior-cover") > mean("ta-muscle")
    );
    // The bead's marking rules: an apparatus response is stored under `*_apparatus_effective` with `apparatus-*`
    // observations, and layer thickness is geometry, never a material row.
    let text = std::fs::read_to_string(repo_root().join(
        "data/matdb/seed-v1/vocalfold-human-indentation-layers-chhetri-jvoice-pmc2891297/properties.tsv",
    ))
    .expect("properties");
    let mut scalars = 0;
    for line in text.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "scalar" => {
                scalars += 1;
                assert!(
                    cols[3].ends_with("_apparatus_effective"),
                    "{} must be an apparatus property",
                    cols[3]
                );
                assert!(
                    cols[2].starts_with("apparatus-"),
                    "{} must name an apparatus observation",
                    cols[2]
                );
                assert!(
                    !cols[3].contains("thickness"),
                    "{}: layer thickness is geometry, not material data",
                    cols[3]
                );
            }
            "observation" => assert!(
                cols[1].starts_with("apparatus-"),
                "observation id {}",
                cols[1]
            ),
            _ => {}
        }
    }
    assert_eq!(
        scalars,
        INDENTATION.len() + 15,
        "22 rows and 5 printed location summaries of three values"
    );

    // Two independent labs, two methods: the transverse stiffness of the cover (indentation of isolated cover,
    // effective) and the transverse E_T of the layered pack's cover specimens agree within a factor 3.
    let cover_et: Vec<f64> = TABLE_2
        .iter()
        .filter(|r| r.1 == "cover")
        .map(|r| r.3)
        .collect();
    let cover_et_mean = cover_et.iter().sum::<f64>() / 6.0;
    let ratio = mean("medial-cover") / cover_et_mean;
    assert!(
        (0.3..3.0).contains(&ratio),
        "Chhetri medial cover vs Kelleher cover E_T mean {cover_et_mean}"
    );
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"g2-human-indentation\",\"verdict\":\"pass\",\
         \"cover_vs_kelleher_et\":{ratio},\"ta_mean_kpa\":{}}}",
        mean("ta-muscle")
    );
}

// ---------------------------------------------------------------------------------------------------------------
// LIP / perioral packs of music bead `frankensim-irp-v5-instrument-realism-pexfp.194` (M8-W2d), the P5 half: the
// in vivo mechanical response of trombonists' lips (Newton 2009 thesis), the lower-lip aspiration modulus (Luboz,
// Promayon and Payan 2014) and the perioral mouth-corner stiffness (Chu, Barlow and colleagues 2009 and 2010).
// Every value here is an APPARATUS response (a whole embouchure, an aspiration device acting on a layered face, a
// mouth-corner stretch): the bead's rule is that such a response is stored under a `*_apparatus_effective` property
// name with an `apparatus-*` observation id and is used only through a forward model of that apparatus. The gates
// pin the typed values to the printed ones, refuse any value without that marking, and run the one cross-check the
// sources allow (the unit of the perioral regression variable is not printed; the OroSTIFF mean must lie inside the
// range of the regression over 0 to 24 mm of stretch, which it does only if the variable is a stretch in mm).
// ---------------------------------------------------------------------------------------------------------------

struct LipPack {
    manifest: String,
    /// claim id -> (property, value, unit, observation id)
    scalar: BTreeMap<String, (String, f64, String, String)>,
    /// claim id -> uncertainty kind
    uncertainty: BTreeMap<String, String>,
    /// claim id -> number of validity rows
    validity: BTreeMap<String, usize>,
    /// observation id -> specimen, method and caveats joined
    observation: BTreeMap<String, String>,
}

impl LipPack {
    fn value(&self, claim: &str) -> f64 {
        self.scalar
            .get(claim)
            .unwrap_or_else(|| panic!("claim {claim} missing"))
            .1
    }

    fn unit(&self, claim: &str) -> &str {
        &self
            .scalar
            .get(claim)
            .unwrap_or_else(|| panic!("claim {claim} missing"))
            .2
    }
}

fn load_lip_pack(dir: &str) -> LipPack {
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(
        manifest.contains("Published-Scientific-Facts-Attribution-Retained"),
        "{dir}: facts-only licence token"
    );
    assert!(manifest.contains("https://"), "{dir}: a locator is required");
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let mut pack = LipPack {
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
                assert!(value.is_finite(), "{dir}: {} finite", cols[1]);
                let previous = pack.scalar.insert(
                    cols[1].to_string(),
                    (
                        cols[3].to_string(),
                        value,
                        cols[5].to_string(),
                        cols[2].to_string(),
                    ),
                );
                assert!(previous.is_none(), "{dir}: duplicate claim {}", cols[1]);
            }
            "uncertainty" => {
                assert_eq!(cols.len(), 7, "{dir}: uncertainty arity");
                assert!(
                    cols[2] == "unstated" || cols[2] == "absolute",
                    "{dir}: unexpected uncertainty kind {}",
                    cols[2]
                );
                pack.uncertainty
                    .insert(cols[1].to_string(), cols[2].to_string());
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
    for (id, (_, _, _, obs)) in &pack.scalar {
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

const LIP_APPARATUS_PACKS: [&str; 3] = [
    "lip-resonance-human-video-newton-thesis2009",
    "lip-lower-aspiration-luboz-abme2014",
    "perioral-stiffness-chu-barlow-jslhr2009",
];

#[test]
fn g2_lip_packs_mark_every_value_as_an_apparatus_response() {
    let mut claims = 0usize;
    for dir in LIP_APPARATUS_PACKS {
        let pack = load_lip_pack(dir);
        assert!(!pack.scalar.is_empty(), "{dir}: no scalar rows");
        for (id, (property, _, _, obs)) in &pack.scalar {
            assert!(
                property.ends_with("_apparatus_effective"),
                "{dir}: {id} property {property} must end with _apparatus_effective"
            );
            assert!(
                obs.starts_with("apparatus-"),
                "{dir}: {id} observation {obs} must start with apparatus-"
            );
            claims += 1;
        }
        for (obs, text) in &pack.observation {
            assert!(
                obs.starts_with("apparatus-"),
                "{dir}: observation id {obs}"
            );
            assert!(
                text.contains("apparatus-effective"),
                "{dir}: {obs} must state that its values are apparatus responses"
            );
        }
    }
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"g2-lip-apparatus-marking\",\"verdict\":\"pass\",\
         \"packs\":{},\"claims\":{claims}}}",
        LIP_APPARATUS_PACKS.len()
    );
}

#[test]
fn g2_newton_human_lip_pack_pins_the_printed_values() {
    let pack = load_lip_pack("lip-resonance-human-video-newton-thesis2009");
    assert!(pack.manifest.contains("era.ed.ac.uk/handle/1842/3140"));
    // (claim id, printed value, unit): the thesis' own sentences (Sec. 4.4.1 and 4.4.2) and Table 4.2, human row.
    let printed: [(&str, f64, &str); 34] = [
        ("aj-bb1-f-out", 38.0, "Hz"),
        ("aj-bb1-f-out-pm", 3.0, "Hz"),
        ("aj-bb1-q-out-low", 1.5, "1"),
        ("aj-bb1-q-out-high", 2.4, "1"),
        ("aj-bb1-q-out-mean", 1.8, "1"),
        ("aj-bb1-q-in-low", 2.0, "1"),
        ("aj-bb1-q-in-high", 2.5, "1"),
        ("aj-bb2-q-out", 1.7, "1"),
        ("aj-bb2-q-in", 5.0, "1"),
        ("aj-bb2-f-in", 140.0, "Hz"),
        ("tj-bb2-f-out-low", 55.0, "Hz"),
        ("tj-bb2-f-out-high", 58.0, "Hz"),
        ("tj-bb2-q-out-low", 1.6, "1"),
        ("tj-bb2-q-out-high", 2.4, "1"),
        ("tj-bb2-f-second-low", 165.0, "Hz"),
        ("tj-bb2-f-second-high", 173.0, "Hz"),
        ("tj-bb2-q-second-low", 2.8, "1"),
        ("tj-bb2-q-second-high", 5.0, "1"),
        ("tj-f3-f-out", 100.0, "Hz"),
        ("tj-f3-q-out", 1.7, "1"),
        ("tj-f3-f-in", 200.0, "Hz"),
        ("tj-f3-q-in-low", 4.0, "1"),
        ("tj-f3-q-in-high", 5.0, "1"),
        ("tj-f3-f-played", 170.0, "Hz"),
        ("t42-fi-low", 32.0, "Hz"),
        ("t42-fi-high", 60.0, "Hz"),
        ("t42-qi-low", 1.5, "1"),
        ("t42-qi-high", 2.4, "1"),
        ("t42-fj-low", 70.0, "Hz"),
        ("t42-fj-high", 100.0, "Hz"),
        ("t42-qj-low", 1.2, "1"),
        ("t42-qj-high", 5.0, "1"),
        ("t42-f0-low", 50.0, "Hz"),
        ("t42-f0-high", 80.0, "Hz"),
    ];
    for (id, want, unit) in printed {
        assert!(
            near(pack.value(id), want),
            "{id}: typed {} against printed {want}",
            pack.value(id)
        );
        assert_eq!(pack.unit(id), unit, "{id}: unit");
    }
    // Every printed range is stored as a low and a high row, never averaged: low < high, and the printed mean sits inside its range.
    for stem in [
        "aj-bb1-q-out",
        "aj-bb1-q-in",
        "tj-bb2-f-out",
        "tj-bb2-q-out",
        "tj-bb2-f-second",
        "tj-bb2-q-second",
        "tj-f3-q-in",
        "t42-fi",
        "t42-qi",
        "t42-fj",
        "t42-qj",
        "t42-f0",
    ] {
        let (lo, hi) = (
            pack.value(&format!("{stem}-low")),
            pack.value(&format!("{stem}-high")),
        );
        assert!(lo < hi, "{stem}: low {lo} must be below high {hi}");
    }
    let (mean, lo, hi) = (
        pack.value("aj-bb1-q-out-mean"),
        pack.value("aj-bb1-q-out-low"),
        pack.value("aj-bb1-q-out-high"),
    );
    assert!(lo < mean && mean < hi, "printed mean Q {mean} inside {lo}..{hi}");
    // The thesis' own statement for these data: the playing frequency of the target lies above the outward resonance.
    assert!(pack.value("aj-bb1-f-out") < 58.0);
    assert!(pack.value("tj-bb2-f-out-high") < 116.0);
    assert!(pack.value("tj-f3-f-out") < pack.value("tj-f3-f-played"));
    // The pack must say that the Table 4.2 second-resonance range and the text disagree, and that nothing is averaged.
    let table = &pack.observation["apparatus-human-table-4-2"];
    assert!(table.contains("never averaged") && table.contains("does not match"));
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"g2-newton-human-lip\",\"verdict\":\"pass\",\
         \"claims\":{},\"outward_hz_aj_bb1\":{}}}",
        printed.len(),
        pack.value("aj-bb1-f-out")
    );
}

#[test]
fn g2_lower_lip_aspiration_and_perioral_stiffness_pin_the_printed_values() {
    let luboz = load_lip_pack("lip-lower-aspiration-luboz-abme2014");
    for (id, want, unit) in [
        ("lower-lip", 33.7, "kPa"),
        ("lower-lip-pm", 7.3, "kPa"),
        ("cheek", 31.0, "kPa"),
        ("cheek-pm", 4.6, "kPa"),
        ("cheekbone", 34.9, "kPa"),
        ("cheekbone-pm", 6.6, "kPa"),
        ("forehead", 17.3, "kPa"),
        ("forehead-pm", 4.1, "kPa"),
    ] {
        assert!(near(luboz.value(id), want), "{id}: typed {}", luboz.value(id));
        assert_eq!(luboz.unit(id), unit, "{id}: unit");
    }
    assert!(
        luboz.observation["apparatus-luboz-aspiration"]
            .contains("standard deviation or a standard error"),
        "the unnamed +- statistic must be recorded as unknown"
    );
    // The abstract's own ordering: the forehead is the softest site, and the lip sits within the +- of the cheeks.
    assert!(luboz.value("forehead") < luboz.value("cheek"));
    assert!((luboz.value("lower-lip") - luboz.value("cheek")).abs() < luboz.value("lower-lip-pm"));

    let chu = load_lip_pack("perioral-stiffness-chu-barlow-jslhr2009");
    // The printed polynomials, N/mm in the paper, typed exactly (x1000) in N/m with the span in mm.
    for (id, want, unit) in [
        ("m-const", 76.7, "N/m"),
        ("m-lin-abs", 3.0, "N/m/mm"),
        ("m-quad", 0.5, "N/m/mm/mm"),
        ("m-r2", 0.58, "1"),
        ("f-const", 94.5, "N/m"),
        ("f-lin-abs", 5.5, "N/m/mm"),
        ("f-quad", 0.5, "N/m/mm/mm"),
        ("f-r2", 0.77, "1"),
        ("orostiff-mean", 73.0, "N/m"),
    ] {
        assert!(near(chu.value(id), want), "{id}: typed {}", chu.value(id));
        assert_eq!(chu.unit(id), unit, "{id}: unit");
    }
    let male = |span: f64| {
        chu.value("m-const") - chu.value("m-lin-abs") * span + chu.value("m-quad") * span * span
    };
    let female = |span: f64| {
        chu.value("f-const") - chu.value("f-lin-abs") * span + chu.value("f-quad") * span * span
    };
    // Recomputed from the typed coefficients against the printed equation 0.0767 - 0.0030 S + 0.0005 S^2 (N/mm).
    for span in [0.0, 3.0, 12.0, 24.0] {
        let printed_n_per_mm = 0.0767 - 0.0030 * span + 0.0005 * span * span;
        assert!(
            (male(span) - 1000.0 * printed_n_per_mm).abs() < 1e-9,
            "male polynomial at {span}"
        );
        let printed_female = 0.0945 - 0.0055 * span + 0.0005 * span * span;
        assert!(
            (female(span) - 1000.0 * printed_female).abs() < 1e-9,
            "female polynomial at {span}"
        );
    }
    // The unit cross-check: the OroSTIFF mean over all trials lies inside the male regression over 0 to 24 mm of
    // stretch (its minimum is at S = 3 mm, 72.2 N/m; its value at 24 mm is 292.7 N/m). If the variable were an
    // absolute span of about 40 to 70 mm the polynomial would exceed 1 N/mm and this would fail.
    let (lo, hi) = (male(3.0), male(24.0));
    assert!(lo < male(0.0) && lo < male(12.0), "S = 3 is the minimum");
    assert!(
        (lo..=hi).contains(&chu.value("orostiff-mean")),
        "OroSTIFF mean {} outside the regression range {lo}..{hi}",
        chu.value("orostiff-mean")
    );
    assert!(
        chu.observation["apparatus-chu-2009-male"].contains("an inference"),
        "the unit reading of the span must be recorded as an inference"
    );
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"g2-lip-aspiration-and-perioral\",\"verdict\":\"pass\",\
         \"lower_lip_kpa\":{},\"male_n_per_m_at_24mm\":{hi}}}",
        luboz.value("lower-lip")
    );
}

#[test]
fn g2_human_cover_shear_pack_pins_the_fits_and_reproduces_the_crossover() {
    let pack = load_human_pack("vocalfold-human-cover-shear-chan-rodriguez-jasa-pmc2561277");

    // Table 2 of Chan and Rodriguez 2008, as printed: p, q and R2 of the power laws fitted to the mean of seven covers.
    for (id, value) in [
        ("storage-p", 17.246),
        ("storage-q", 0.6320),
        ("storage-r2", 0.9493),
        ("loss-p", 58.448),
        ("loss-q", 0.3243),
        ("loss-r2", 0.8363),
        ("viscosity-p", 10.509),
        ("viscosity-q", 0.7391),
        ("viscosity-r2", 0.9934),
        ("damping-ratio-100-250", 0.7383),
        ("crossover-frequency", 50.0),
    ] {
        assert!(near(pack.value(id), value), "{id}");
        assert_eq!(
            pack.uncertainty[id].0,
            if id == "damping-ratio-100-250" {
                "absolute"
            } else {
                "unstated"
            }
        );
    }
    assert_eq!(pack.uncertainty["damping-ratio-100-250"].1, Some(0.1211));
    let state = &pack.observation["cover-seven-subjects"];
    for needle in [
        "seven subjects",
        "SIMPLE-SHEAR",
        "1.0 percent shear strain",
        "MEAN of the seven",
        "total laryngectomy",
        "must not be extrapolated",
    ] {
        assert!(
            state.contains(needle),
            "the observation must carry: {needle}"
        );
    }

    let (pg, qg) = (pack.value("storage-p"), pack.value("storage-q"));
    let (pl, ql) = (pack.value("loss-p"), pack.value("loss-q"));
    let (pe, qe) = (pack.value("viscosity-p"), pack.value("viscosity-q"));
    let g_storage = |f: f64| pg * f.powf(qg);
    let g_loss = |f: f64| pl * f.powf(ql);

    // The crossover the authors report ("around 50 Hz") follows from the two printed power laws, and G'' is above G'
    // below it and under it above it.
    let crossover = (pl / pg).powf(1.0 / (qg - ql));
    assert!(
        (40.0..65.0).contains(&crossover),
        "crossover {crossover} Hz against the printed about 50 Hz"
    );
    assert!(g_loss(10.0) > g_storage(10.0) && g_loss(200.0) < g_storage(200.0));

    // The derived evaluations equal the printed laws at the phonatory frequencies (stored rounded to 1e-3 Pa).
    for (id, f, value) in [
        ("derived-gprime-100", 100.0, g_storage(100.0)),
        ("derived-gprime-250", 250.0, g_storage(250.0)),
        ("derived-gdoubleprime-100", 100.0, g_loss(100.0)),
        ("derived-gdoubleprime-250", 250.0, g_loss(250.0)),
    ] {
        assert!(
            (pack.value(id) - value).abs() <= 1e-3 * 0.5 + 1e-9 * value,
            "{id} at {f} Hz: {} against {value}",
            pack.value(id)
        );
        assert_eq!(pack.scalar[id].1, "Pa");
        assert!(pack.observation["derived-evaluations"].contains("COMPUTED"));
    }

    // Consistency between the three printed fits: eta' = G''/(2 pi f) within 30 percent at 100 Hz, and the printed
    // mean loss factor over 100-250 Hz within 10 percent of the mean of G''/G' from the fits.
    let eta_100 = pe * 100.0_f64.powf(-qe);
    let eta_from_loss = g_loss(100.0) / (2.0 * core::f64::consts::PI * 100.0);
    assert!(
        (eta_100 / eta_from_loss - 1.0).abs() < 0.3,
        "eta' {eta_100} against G''/omega {eta_from_loss}"
    );
    let samples = 1501_u32;
    let mean_loss_factor = (0..samples)
        .map(|i| {
            let f = 100.0 + 150.0 * f64::from(i) / f64::from(samples - 1);
            g_loss(f) / g_storage(f)
        })
        .sum::<f64>()
        / f64::from(samples);
    let printed = pack.value("damping-ratio-100-250");
    assert!(
        (mean_loss_factor / printed - 1.0).abs() < 0.1,
        "mean G''/G' {mean_loss_factor} against printed {printed}"
    );

    // Contrast with the felid register: the human cover's phonatory-range loss factor lies far above the large-felid
    // lamina propria tan delta of the older pack, and the source itself says it exceeds the 0.1-0.2 of fold models.
    let felid = load_pack("vocalfold-felid-rheometry-plos-pone0027029");
    assert!(printed > 2.0 * felid.scalars["loss_tangent_tiger"]);
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"g2-human-cover-shear\",\"verdict\":\"pass\",\
         \"crossover_hz\":{crossover},\"mean_loss_factor_100_250\":{mean_loss_factor}}}"
    );
}
