//! Physicality gates for the wound-string and multifilament-core matdb packs of the .192 winding hunt (bead
//! `frankensim-irp-v5-instrument-realism-pexfp.192`, plan item M8-W2b): four cello G strings (Lampis et al. 2024, CC BY 4.0), the cores of
//! two patents (Behringer US 3,610,084; Kureha US 5,427,008), a wound bass-guitar string (Bavu, Smith and Wolfe 2005) and the winding and
//! core densities of commercial guitar strings (Parker and Rohrbacher US 6,348,646).
//!
//! Each gate is an identity the source's own numbers must satisfy: wave impedance = sqrt(tension x linear density) and tension =
//! linear density x (2 L f)^2; a wound string far less stiff in bending than a solid steel rod of its diameter; 7 filaments of 205 denier
//! = the printed 1435 denier; the invention's core weight per metre = yarns x 300 denier; a twisted yarn's printed diameter at least that
//! of the same mass at solid density; the second torsional resonance near twice the first and the effective diameter implied by the
//! torsional stiffness and the torsional wave speed the size of a bass string; densities in the band of the named alloys. Schema hygiene
//! (every scalar has an uncertainty and a validity row) for every pack.

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
        let c: Vec<&str> = line.split('\t').collect();
        match c[0] {
            "scalar" => {
                let v: f64 = c[4].parse().expect("value");
                assert!(v.is_finite() && v > 0.0, "{dir}: {} must be positive", c[3]);
                scalars.insert(c[3].to_string(), v);
                ids.push(c[1].to_string());
            }
            "uncertainty" => unc.push(c[1].to_string()),
            "validity" => val.push(c[1].to_string()),
            "observation" => {}
            other => panic!("unknown row kind {other}"),
        }
    }
    for id in &ids {
        assert!(
            unc.contains(id) && val.contains(id),
            "{dir}: {id} needs uncertainty and validity rows"
        );
    }
    scalars
}

const PERMITTED_FACTS: &str = "Published-Scientific-Facts-Attribution-Retained";
const DENIER: f64 = 1e-3 / 9000.0;

#[test]
fn cello_g_strings_satisfy_the_string_identities() {
    let s = load_pack("cello-g-strings-whole-lampis-2024", "CC-BY-4.0");
    let l = s["playing_length"];
    for name in ["prelude", "helicore", "dominant", "kaplan"] {
        let (t, mu, z) = (
            s[&format!("tension_{name}")],
            s[&format!("linear_density_{name}")],
            s[&format!("wave_impedance_{name}")],
        );
        assert!(
            ((t * mu).sqrt() - z).abs() / z < 0.01,
            "{name}: Z = sqrt(T mu)"
        );
        // Nominal G2 = 98 Hz on the 0.70 m monochord.
        let t_nominal = mu * (2.0 * l * 98.0_f64).powi(2);
        assert!((t_nominal - t).abs() / t < 0.02, "{name}: T = mu (2 L f)^2");
        // A wound string is far less stiff than a solid steel rod of its diameter.
        let d = s[&format!("diameter_{name}")];
        let solid_steel = 2.1e11 * std::f64::consts::PI * (d / 2.0).powi(4) / 4.0;
        let b = s[&format!("bending_stiffness_{name}")];
        assert!(
            b < 0.1 * solid_steel,
            "{name}: B {b} against solid steel {solid_steel}"
        );
    }
    // The second paper's session of the same Dominant string.
    let (t, f0, mu, z) = (
        s["tension_dominant_acta_session"],
        s["fundamental_frequency_dominant_acta_session"],
        s["linear_density_dominant_acta_session"],
        s["wave_impedance_dominant_acta_session"],
    );
    assert!(
        (t / (2.0 * l * f0).powi(2) - mu).abs() / mu < 0.005,
        "l = T/(2 L f0)^2"
    );
    assert!(((mu * t).sqrt() - z).abs() / z < 0.005, "Z = sqrt(l T)");
    assert!(
        (t - s["tension_dominant"]).abs() / t < 0.005,
        "same string, two sessions"
    );
    assert!(
        (s["bending_stiffness_dominant_acta_session"] - s["bending_stiffness_dominant"]).abs()
            < 1e-5,
        "same B in both papers"
    );
}

#[test]
fn multifilament_core_patents_satisfy_their_identities() {
    let s = load_pack(
        "string-core-multifilament-patents-behringer-kureha",
        PERMITTED_FACTS,
    );
    // Behringer: 7 filaments of 205 denier are the printed 1435 denier.
    let total = s["behringer_violin_e_core_total_linear_density"];
    let sum = s["behringer_violin_e_core_filament_count"]
        * s["behringer_violin_e_core_filament_linear_density"];
    assert!(
        (total - sum).abs() / total < 1e-3,
        "7 x 205 denier = 1435 denier"
    );
    assert!((total - 1435.0 * DENIER).abs() / total < 1e-9);
    for key in [
        "behringer_core_breaking_elongation_range",
        "behringer_core_tensile_strength_range",
        "behringer_core_young_modulus_range",
    ] {
        assert!(
            s[&format!("{key}_low")] < s[&format!("{key}_high")],
            "{key} ordered"
        );
    }
    // Kureha: the invention's core weight per metre is its yarns x 300 denier; a twisted yarn is no thinner than its mass at solid density.
    let yarn = s["kureha_pvdf_multifilament_yarn_linear_density"];
    let solid_density = 1780.0; // specific gravity 1.78 against water at 1000 kg/m^3
    for (tag, st) in [("example1_8_yarns", "a"), ("example2_6_yarns", "g")] {
        let n = s[&format!("kureha_{tag}_yarn_count")];
        let dia = s[&format!("kureha_{tag}_core_diameter")];
        let d_solid = (4.0 * n * yarn / solid_density / std::f64::consts::PI).sqrt();
        assert!(
            dia > d_solid,
            "{tag}: printed {dia} against solid {d_solid}"
        );
        // The invention's A string was made on the 8-yarn core and its G string on the 6-yarn core (Table 1 of the patent).
        let core = s[&format!("kureha_table2_invention_pvdf_{st}_core_linear_density")];
        assert!(
            (core - n * yarn).abs() / core < 0.07,
            "{tag}: printed {core} against {}",
            n * yarn
        );
        assert!(
            s[&format!("kureha_{tag}_young_modulus")] >= s["kureha_claimed_young_modulus_minimum"],
            "{tag}: modulus above the claimed minimum"
        );
    }
    // Table 2: a finished string weighs more than its core, for all nine.
    for kind in ["gut_eudoxa", "invention_pvdf", "nylon_dominant"] {
        for st in ["a", "d", "g"] {
            let core = s[&format!("kureha_table2_{kind}_{st}_core_linear_density")];
            let whole = s[&format!("kureha_table2_{kind}_{st}_string_linear_density")];
            assert!(
                whole > core,
                "{kind} {st}: string {whole} against core {core}"
            );
        }
    }
}

#[test]
fn wound_bass_string_torsion_is_the_size_of_a_bass_string() {
    let s = load_pack(
        "string-wound-bass-guitar-bavu-smith-wolfe-2005",
        PERMITTED_FACTS,
    );
    let (l, mu, k) = (
        s["bridge_separation"],
        s["linear_density"],
        s["torsional_stiffness_in_situ"],
    );
    let mut prev = f64::INFINITY;
    for case in ["a", "b"] {
        let (f1, f2) = (
            s[&format!("first_torsional_resonance_case_{case}")],
            s[&format!("second_torsional_resonance_case_{case}")],
        );
        assert!(
            (1.9..2.1).contains(&(f2 / f1)),
            "case {case}: f2/f1 = {}",
            f2 / f1
        );
        assert!(
            f1 < prev,
            "the torsional frequency falls as the tension rises"
        );
        prev = f1;
        // Torsional wave speed c = 2 L f1 (fixed ends); polar moment per length = K / c^2; a uniform cylinder of mass per length mu
        // and that moment has diameter 2 sqrt(2 theta / mu).
        let c = 2.0 * l * f1;
        let d_eff = 2.0 * (2.0 * (k / (c * c)) / mu).sqrt();
        assert!(
            (2.0e-3..4.0e-3).contains(&d_eff),
            "case {case}: effective diameter {d_eff} m"
        );
    }
    assert!(s["second_torsional_resonance_case_b"] < s["second_torsional_resonance_case_a"]);
    assert!(s["transverse_fundamental_case_a"] < s["transverse_fundamental_case_b"]);
}

#[test]
fn guitar_string_winding_densities_are_in_the_alloy_bands() {
    let s = load_pack(
        "guitar-string-winding-densities-parker-rohrbacher-us6348646",
        PERMITTED_FACTS,
    );
    let pb = s["phosphor_bronze_winding_wire_density_ej16_e_string"];
    assert!((8700.0..9000.0).contains(&pb), "phosphor bronze {pb}");
    let (ni, pearse_w, pearse_c) = (
        s["nickel_plated_winding_density_super_slinky"],
        s["winding_density_john_pearse_medium_measured"],
        s["core_density_john_pearse_medium_measured"],
    );
    for (name, v) in [
        ("nickel plated", ni),
        ("pearse winding", pearse_w),
        ("pearse core", pearse_c),
    ] {
        assert!((7000.0..9200.0).contains(&v), "{name}: {v}");
    }
    assert!(
        ni < pb,
        "the nickel-plated winding as measured is lighter than the phosphor bronze"
    );
    assert!(
        s["phosphor_bronze_winding_wire_diameter_ej16_e_string"]
            > s["nickel_plated_winding_diameter_super_slinky"]
    );
}
