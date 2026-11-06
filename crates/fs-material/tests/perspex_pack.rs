//! Physicality gates for the perspex pack of Smith and Woodhouse 2000 (Table 1), the bounding solid of
//! their steady-sliding rig. Woodhouse and Galluzzo's (2025) eq. (13) refits that rig's data, so the
//! `fs_tribo::rosin_thermal` calibration conducts into this solid.
//!
//! The gates encode what cannot be wrong for a clear cast PMMA: positivity, and schema hygiene
//! (every scalar has an uncertainty and a validity row). Density, conductivity, specific heat and
//! diffusivity must agree with an independent PMMA source: the PLEXIGLAS GS reference pack, a
//! manufacturer's typical values, to 5 % for density and specific heat, and 10 % for conductivity
//! and the diffusivity that carries it (0.2 against 0.19 W/m/K is 5.3 %).

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

#[test]
fn perspex_pack_passes_physicality_gates() {
    let sw = load_pack(
        "perspex-smith-woodhouse-jmps2000",
        "Published-Scientific-Facts-Attribution-Retained",
    );
    let gs = scalars_of("plexiglas-gs-clear-reference");
    let diffusivity = |p: &BTreeMap<String, f64>| {
        p["thermal_conductivity"] / (p["density"] * p["specific_heat_capacity"])
    };
    let checks = [
        ("density", sw["density"], gs["density"], 0.05),
        (
            "thermal_conductivity",
            sw["thermal_conductivity"],
            gs["thermal_conductivity"],
            0.10,
        ),
        (
            "specific_heat_capacity",
            sw["specific_heat_capacity"],
            gs["specific_heat_capacity"],
            0.05,
        ),
        ("diffusivity", diffusivity(&sw), diffusivity(&gs), 0.10),
    ];
    for (name, ours, reference, tol) in checks {
        let rel = (ours - reference).abs() / reference;
        println!(
            "{{\"pack\":\"perspex-smith-woodhouse-jmps2000\",\"check\":\"{name}\",\"value\":{ours},\"plexiglas_gs\":{reference},\"rel\":{rel},\"tol\":{tol}}}"
        );
        assert!(
            rel <= tol,
            "{name}: {ours} against PLEXIGLAS GS {reference}"
        );
    }
}
