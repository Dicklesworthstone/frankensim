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

/// The Dominant's two loss rows are derived, not printed. The paper prints one loss number for the string, the quality factor averaged
/// over the first 30 modes of the pluck decays (1064); the rows read it as a constant loss factor, eta_F = eta_B = 1 / 1064. With the two
/// equal, the friction and bending terms of Woodhouse 2004 eq. (8) give zeta = (T eta_F + EI eta_B k^2) / (2 (T + EI k^2)) = eta / 2 for
/// every mode number k, so Q_k = 1 / eta = 1064 for all of them.
#[test]
fn the_dominant_loss_rows_read_the_printed_mean_q_as_a_constant_loss_factor() {
    let s = load_pack("cello-g-strings-whole-lampis-2024", "CC-BY-4.0");
    let q = s["mean_quality_factor_first_30_modes_dominant"];
    assert_eq!(q, 1064.0);
    let (eta_f, eta_b) = (s["loss_factor_tension_dominant"], s["loss_factor_bending_dominant"]);
    // Two distinct property names (the engine's selection needs two) holding the same number, 1 / Q.
    assert_eq!(eta_f, eta_b);
    assert!((eta_f * q - 1.0).abs() < 1e-12, "eta = {eta_f:e} is not 1 / {q}");
    // The law, evaluated from the pack's own tension and bending stiffness for 64 modes: zeta_k = eta / 2, Q_k = 1064.
    let (t, ei) = (s["tension_dominant"], s["bending_stiffness_dominant"]);
    for k in 1..=64 {
        let k2 = f64::from(k * k);
        let zeta = (t * eta_f + ei * eta_b * k2) / (2.0 * (t + ei * k2));
        assert!((1.0 / (2.0 * zeta) - q).abs() < 1e-9 * q, "mode {k}: Q = {}", 1.0 / (2.0 * zeta));
    }
    // Planted negatives: only the Dominant has the derived rows, and they sit under their own observation, which says so.
    for other in ["prelude", "helicore", "kaplan"] {
        assert!(!s.contains_key(&format!("loss_factor_tension_{other}")), "{other}");
    }
    let text = std::fs::read_to_string(
        repo_root().join("data/matdb/seed-v1/cello-g-strings-whole-lampis-2024/properties.tsv"),
    )
    .expect("properties");
    let observation = text
        .lines()
        .find(|line| line.starts_with("observation\tlampis-jasael-2024-loss-derived\t"))
        .expect("the derived rows have their own observation");
    assert!(observation.contains("DERIVED row, not a measurement"), "{observation}");
    let rows = text.lines().filter(|line| line.contains("\tlampis-jasael-2024-loss-derived\t")).count();
    assert_eq!(rows, 3, "the observation and its two scalars");
}

/// The Kaplan KS510 cello strings behind the posted C2, D3 and A3 notes: tension, mass per length, bending stiffness and the two loss
/// coefficients of Mansour, Woodhouse and Scavone (author manuscript, Table 1; Table 2.1 of Mansour's 2016 McGill thesis), and the outer
/// diameters of Percival 2013 Table 3.2 (cello-I, strings of unknown make, +-10 percent), which the Kaplan table does not print.
/// Each row is checked against the number printed in the table, then against the table's own identities: the printed impedance
/// Z0 = sqrt(T ms), T = ms (2 L f0)^2 on the 0.69 m effective length, and the printed inharmonicity B = pi^2 EI / (T L^2).
/// The loss rows are the cards' `loss-factors` rows: (eta_F, eta_B) = (1.2e-4, 0.047) C2, (2.3e-4, 0.125) D3, (2.2e-4, 0.114) A3.
#[test]
fn kaplan_ks510_rows_are_the_mansour_table_and_satisfy_its_identities() {
    let s = load_pack("cello-ks510-mansour-2016-percival-2013", PERMITTED_FACTS);
    let l = s["playing_length"];
    assert_eq!(l, 0.69, "effective open-string length of Table 1");
    // note, f0 [Hz], T [N], ms [kg/m], EI [N m^2], eta_F, eta_B, printed Z0 [kg/s], printed B, Percival diameter [m]
    let table = [
        ("a3", 220.0, 171.0, 1.85e-3, 3.26e-4, 2.2e-4, 0.114, 0.56, 39.5e-6, 0.68e-3),
        ("d3", 146.8, 135.9, 3.31e-3, 2.48e-4, 2.3e-4, 0.125, 0.67, 37.9e-6, 0.80e-3),
        ("c2", 65.4, 131.5, 16.14e-3, 6.20e-4, 1.2e-4, 0.047, 1.46, 97.8e-6, 1.40e-3),
    ];
    for (n, f0, t, ms, ei, eta_f, eta_b, z0, b, d) in table {
        for (stem, printed) in [
            ("tension", t),
            ("linear_density", ms),
            ("bending_stiffness", ei),
            ("loss_friction", eta_f),
            ("loss_bending", eta_b),
            ("diameter", d),
        ] {
            assert_eq!(s[&format!("{stem}_ks510_{n}")], printed, "{stem} {n}");
        }
        assert!(((t * ms).sqrt() / z0 - 1.0).abs() < 0.01, "{n}: Z0 = sqrt(T ms)");
        assert!((ms * (2.0 * l * f0).powi(2) / t - 1.0).abs() < 0.01, "{n}: T = ms (2 L f0)^2");
        let b_law = std::f64::consts::PI.powi(2) * ei / (t * l * l);
        assert!((b_law / b - 1.0).abs() < 0.01, "{n}: B = pi^2 EI / (T L^2) = {b_law:e}");
    }
    // Only the three notes of the posted cards: the Table 1 G2 column is not entered, and neither is eta_A.
    for absent in ["tension_ks510_g2", "loss_friction_ks510_g2", "loss_air_ks510_a3"] {
        assert!(!s.contains_key(absent), "{absent}");
    }
    // The two observations say what they are; the diameters' one says they are NOT the Kaplan strings.
    let text = std::fs::read_to_string(
        repo_root().join("data/matdb/seed-v1/cello-ks510-mansour-2016-percival-2013/properties.tsv"),
    )
    .expect("properties");
    let percival = text
        .lines()
        .find(|line| line.starts_with("observation\tpercival-2013-table32\t"))
        .expect("the diameters' observation");
    assert!(percival.contains("NOT the Kaplan strings"), "{percival}");
}

/// The Thomastik Spirocore Orchestra 4/4 medium double-bass strings behind the posted bass notes E1, A1, D2 and G2: the catalogue tensions at the 110.0 cm vibrating length
/// (E 34.0, A 34.0, D 32.0, G 30.0 kg, read as kilogram-force), the line masses derived from them by Mersenne's law at the equal-tempered open-string frequencies (A4 = 440 Hz),
/// the gauges of a player-compiled table (read as inches) and the two bounds on the unpublished core bending stiffness: a solid steel rod of core diameter 0.75 mm x sqrt(T / T_A)
/// and 7/81 of it for a close-packed bundle of seven free wires. Each row is checked against its derivation, and the derived mass against the gauge: the string's mean density
/// stays below solid steel's.
#[test]
fn spirocore_bass_rows_are_their_catalogue_derivations() {
    let s = load_pack("bass-strings-thomastik-spirocore-orchestra-catalogue-2019", PERMITTED_FACTS);
    let l = s["playing_length"];
    assert_eq!(l, 1.1, "the catalogue's 110.0 cm vibrating length");
    let near = |a: f64, b: f64| ((a - b) / b).abs() < 1e-12;
    let t_a = 34.0 * 9.80665;
    // note, catalogue tension [kg], semitones below A4, gauge [in]
    for (n, kg, semitones, inch) in [
        ("e1", 34.0, -41.0, 0.111),
        ("a1", 34.0, -36.0, 0.085),
        ("d2", 32.0, -31.0, 0.070),
        ("g2", 30.0, -26.0, 0.052),
    ] {
        let t = kg * 9.80665;
        let f = 440.0 * 2f64.powf(semitones / 12.0);
        let ms = t / (2.0 * l * f).powi(2);
        let d = inch * 0.0254;
        assert!(near(s[&format!("tension_bass_{n}")], t), "tension {n}");
        assert!(near(s[&format!("linear_density_bass_{n}")], ms), "Mersenne mass {n}");
        assert!(near(s[&format!("diameter_bass_{n}")], d), "gauge {n}");
        let core = 0.75e-3 * (t / t_a).sqrt();
        let hi = 200e9 * std::f64::consts::PI * core.powi(4) / 64.0;
        assert!(near(s[&format!("bending_stiffness_bass_{n}_hi")], hi), "upper bound {n}");
        assert!(near(s[&format!("bending_stiffness_bass_{n}_lo")], hi * 7.0 / 81.0), "lower bound {n}");
        let rho = ms / (std::f64::consts::PI * d * d / 4.0);
        assert!(rho < 7850.0, "{n}: the mean density {rho} kg/m3 must stay below steel's");
    }
}

/// The SURROGATE string loss of the posted surrogate bass pages: the Spirocore catalogue seed with one observation and two scalar rows added, the loss factors 1/900 (the translational-wave Q of typically 900
/// that Bavu, Smith and Wolfe (2005) report for a bass-guitar string, not a measurement of any Spirocore string). Every catalogue row is kept unchanged, the two new rows are 1/900 on both factors, and the
/// observation says in its own words that it is a surrogate.
#[test]
fn spirocore_surrogate_seed_is_the_catalogue_seed_plus_two_loss_rows_of_one_over_900() {
    let parent = load_pack("bass-strings-thomastik-spirocore-orchestra-catalogue-2019", PERMITTED_FACTS);
    let dir = "bass-strings-thomastik-spirocore-orchestra-surrogate-loss-bavu2005";
    let s = load_pack(dir, PERMITTED_FACTS);
    for (k, v) in &parent {
        assert_eq!(s.get(k), Some(v), "{k}: the surrogate seed keeps the catalogue row");
    }
    assert_eq!(s.len(), parent.len() + 2, "exactly two scalar rows are added");
    for k in ["loss_factor_tension_bass_surrogate", "loss_factor_bending_bass_surrogate"] {
        assert!(((s[k] - 1.0 / 900.0) * 900.0).abs() < 1e-12, "{k} is 1/900, got {}", s[k]);
    }
    let text = std::fs::read_to_string(repo_root().join(format!("data/matdb/seed-v1/{dir}/properties.tsv"))).expect("properties");
    assert!(text.contains("SURROGATE row, not a measurement of any Spirocore string"), "the observation names itself a surrogate");
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
