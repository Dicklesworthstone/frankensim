//! The shielded-string pack of Woodhouse (2003), Acta Acustica united with Acustica 89:355-368: Table I's "String"
//! column, the stationary surface of that paper's bowed-string thermal simulations (a metal string thermally shielded
//! by a transferred rosin layer). Galluzzo, Woodhouse and Mansour (2017) state they used it unchanged. The thermal
//! rosin row of `material_strings` can name it as its string-surface pack.
//!
//! The gates: schema hygiene (every scalar has an uncertainty and a validity row), the three values as printed, and
//! the printed Eq. (3) thermal parameter lambda0 = K/sqrt(D) = sqrt(K rho c): 2050 for this column, and 170 for the
//! Rosin column, whose density, specific heat and conductivity are those of the Smith and Woodhouse (2000) rosin pack.
//! Table I prints lambda0 to the tens, so each must agree within 5. The paper says the column is the geometric mean
//! of an unnamed typical steel and polymer, so each value must also lie strictly between the steel music-wire pack and
//! the perspex pack; no equality with those two packs is claimed.

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

fn read(dir: &str, file: &str) -> String {
    std::fs::read_to_string(repo_root().join(format!("data/matdb/seed-v1/{dir}/{file}")))
        .expect(file)
}

/// The v1 pack's scalars by property name, with its schema hygiene checked.
fn load_pack(dir: &str, license: &str) -> BTreeMap<String, f64> {
    let manifest = read(dir, "manifest.tsv");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(
        manifest.contains(&format!("\nlicense\t{license}\n")),
        "{dir}: license {license}"
    );
    let text = read(dir, "properties.tsv");
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

/// Scalar rows of any pack, by property name; other row kinds are skipped.
fn scalars_of(dir: &str) -> BTreeMap<String, f64> {
    read(dir, "properties.tsv")
        .lines()
        .filter_map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            (cols[0] == "scalar").then(|| (cols[3].to_string(), cols[4].parse().expect("value")))
        })
        .collect()
}

/// Eq. (3) of Woodhouse (2003): lambda0 = K/sqrt(D) = sqrt(K rho c).
fn lambda0(p: &BTreeMap<String, f64>) -> f64 {
    (p["thermal_conductivity"] * p["density"] * p["specific_heat_capacity"]).sqrt()
}

#[test]
fn shielded_string_pack_matches_woodhouse_2003_table_i() {
    let lic = "Published-Scientific-Facts-Attribution-Retained";
    let s = load_pack("string-shielded-woodhouse-acta2003", lic);
    let names: Vec<&str> = s.keys().map(String::as_str).collect();
    assert_eq!(names, ["density", "specific_heat_capacity", "thermal_conductivity"]);
    // As printed in Table I, column String.
    assert_eq!(s["density"], 3000.0);
    assert_eq!(s["specific_heat_capacity"], 1000.0);
    assert_eq!(s["thermal_conductivity"], 1.4);
    // Eq. (3), printed to the tens: String 2050, Rosin 170.
    let rosin = load_pack("rosin-smith-woodhouse-jmps2000", lic);
    for (column, ours, printed) in [("String", lambda0(&s), 2050.0), ("Rosin", lambda0(&rosin), 170.0)] {
        println!("{{\"table\":\"W2003 Table I\",\"column\":\"{column}\",\"lambda0\":{ours},\"printed\":{printed}}}");
        assert!((ours - printed).abs() <= 5.0, "{column}: lambda0 {ours} against the printed {printed}");
    }
    // A geometric mean of a steel and a polymer lies strictly between them.
    let steel = scalars_of("steel-music-wire-string-heptonstall-ligo-p1400229");
    let polymer = scalars_of("perspex-smith-woodhouse-jmps2000");
    for name in ["density", "specific_heat_capacity", "thermal_conductivity"] {
        let (lo, hi) = (steel[name].min(polymer[name]), steel[name].max(polymer[name]));
        assert!(lo < s[name] && s[name] < hi, "{name}: {} not strictly between {lo} and {hi}", s[name]);
    }
}
