//! Gates for the two winding-metal packs of the .192 winding hunt (bead `frankensim-irp-v5-instrument-realism-pexfp.192`): nickel from NBS
//! Circular 485 (1950) and swaged tungsten rod from NASA TN D-6794 (1972), both US Government works.
//!
//! Nickel: the printed shear modulus of annealed nickel against E / (2 (1 + nu)) from the printed Young's modulus and Poisson ratio (to
//! 8 percent), cast above commercial density, the order of the reported moduli (sand-cast < average < cold-drawn wire < Z-B wire).
//! Tungsten: the strain-gage modulus within 2 percent of the extensometer one, the conclusion's value between them, the modulus at
//! 1900 K between 0.5 and 0.7 of the room-temperature one, and the single-observation precision above the strain-gage uncertainty.

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
    assert!(manifest.contains("\nlicense\tWork-of-the-US-Government-Public-Use-Permitted\n"));
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    let (mut s, mut ids) = (BTreeMap::new(), Vec::new());
    let (mut unc, mut val) = (Vec::new(), Vec::new());
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        match c[0] {
            "scalar" => {
                let v: f64 = c[4].parse().expect("value");
                assert!(v.is_finite() && v > 0.0, "{} must be positive", c[3]);
                s.insert(c[3].to_string(), v);
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
            "{dir}: {id} uncertainty and validity rows"
        );
    }
    s
}

#[test]
fn nickel_pack_is_consistent() {
    let s = load_pack("nickel-nbs-circular485-1950");
    let (e, nu) = (
        s["young_modulus_nickel_average_about"],
        s["poisson_ratio_nickel"],
    );
    let g_from_e = e / (2.0 * (1.0 + nu));
    let g = s["shear_modulus_annealed_nickel_99_2_percent"];
    assert!(
        (g - g_from_e).abs() / g < 0.08,
        "G {g} against E/(2(1+nu)) {g_from_e}"
    );
    assert!(s["density_high_purity_cast_nickel"] > s["density_commercial_malleable_nickel_about"]);
    assert!((8800.0..8950.0).contains(&s["density_commercial_malleable_nickel_about"]));
    assert!(
        s["young_modulus_a_nickel_sand_cast"] < e
            && e < s["young_modulus_cold_drawn_nickel_wire"]
            && s["young_modulus_cold_drawn_nickel_wire"] < s["young_modulus_z_b_nickel_wire"]
    );
    assert!(s["shear_modulus_unannealed_nickel"] < g);
}

#[test]
fn tungsten_rod_pack_is_consistent() {
    let s = load_pack("tungsten-rod-nasa-tnd6794-1972");
    let (ext, gage, concl) = (
        s["young_modulus_tungsten_room_temperature_extensometer"],
        s["young_modulus_tungsten_room_temperature_strain_gages"],
        s["young_modulus_tungsten_room_temperature_conclusion"],
    );
    assert!(
        (gage - ext).abs() / ext < 0.02,
        "strain gages {gage} against extensometer {ext}"
    );
    assert!(
        ext <= concl && concl <= gage,
        "the conclusion's value lies between the two"
    );
    let ratio = s["young_modulus_tungsten_1900_k"] / concl;
    assert!((0.5..0.7).contains(&ratio), "E(1900 K) / E(room) = {ratio}");
    assert!(
        s["single_observation_precision_one_sigma_relative"]
            > s["strain_gage_uncertainty_relative_about"]
    );
}
