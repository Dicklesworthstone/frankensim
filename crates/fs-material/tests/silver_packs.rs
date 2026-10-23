//! Physicality gates for the silver string packs (bead `frankensim-irp-v5-instrument-realism-pexfp.369`: silver piano strings as a
//! data swap): hard-drawn silver from NBS Circular C447 (1943) Table 60 serial 3801 and Table 67, and the drawn-hard tensile strength
//! from Bureau of Standards Circular 101 (1921) Table 25 (US Government works, public domain).
//!
//! The gates encode what cannot be wrong whichever silver was tested: positivity; the printed kips/in2, kg/mm2, cal and W/cm values
//! reproduced from the SI rows at the stated conversion factors; drawing raising the modulus above the annealed element's; Poisson's
//! ratio in the metals' band; Dulong-Petit; and schema hygiene. The strength claim states no wire diameter (nothing sourced bounds the
//! diameter effect), the packs hold no yield, and the strength lies below the stress any concert-pitch string carries.

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

const MPA_PER_KSI: f64 = 6.894757293168;
const MPA_PER_KGF_MM2: f64 = 9.80665;

/// Scalars [SI] of one pack after its manifest and schema checks.
fn scalars(dir: &str, citation: &str) -> BTreeMap<String, f64> {
    let root = repo_root().join("data/matdb/seed-v1").join(dir);
    let manifest = std::fs::read_to_string(root.join("manifest.tsv")).expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v3\n"));
    assert!(manifest.contains(&format!("\npack_id\t{dir}\n")));
    assert!(manifest.contains("\nlicense\tWork-of-the-US-Government-Public-Use-Permitted\n"));
    assert!(manifest.contains(citation), "{dir}: {citation}");
    assert!(manifest.contains("\naxis\tT\tT\tabsolute-temperature\n"));
    // Neither pack states a wire diameter: no silver source resolves one.
    assert!(!manifest.contains("wire_diameter"), "{dir}");
    let text = std::fs::read_to_string(root.join("properties.tsv")).expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let mut s = BTreeMap::new();
    let (mut ids, mut unc, mut val) = (Vec::new(), Vec::new(), Vec::new());
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        match c[0] {
            "scalar" => {
                let v: f64 = c[4].parse().expect("value");
                assert!(v.is_finite() && v > 0.0, "{} must be positive", c[3]);
                let si = match c[5] {
                    "MPa" => v * 1e6,
                    "kg/m^3" | "1" | "K^-1" | "J/kg/K" | "W/m/K" => v,
                    unit => panic!("{}: unexpected unit {unit}", c[3]),
                };
                assert!(s.insert(c[3].to_string(), si).is_none(), "{} is duplicated", c[3]);
                assert_eq!(c[6], "constant", "{}", c[3]);
                ids.push(c[1].to_string());
            }
            "uncertainty" => unc.push(c[1].to_string()),
            "validity" => {
                assert_eq!(&c[2..], ["T", "293.15", "293.15", "K"], "{}: one admitted temperature", c[1]);
                val.push(c[1].to_string());
            }
            "observation" => {}
            other => panic!("unknown row kind {other}"),
        }
    }
    for id in &ids {
        assert!(unc.contains(id) && val.contains(id), "{id}: uncertainty and validity rows");
    }
    s
}

#[test]
fn silver_packs_pass_physicality_gates() {
    let s = scalars("silver-pure-nbs-c447", "nbscircular447.pdf");
    let names: Vec<&str> = s.keys().map(String::as_str).collect();
    assert_eq!(names, ["density", "linear_thermal_expansion_coefficient", "poisson_ratio",
        "specific_heat_capacity", "thermal_conductivity", "young_modulus"]);
    // Serial 3801 as printed: 10,900 kips/in2, specific gravity 10.43, Poisson 0.392.
    assert!((s["young_modulus"] / 1e6 / MPA_PER_KSI / 10_900.0 - 1.0).abs() < 5e-9);
    assert_eq!(s["density"], 10_430.0);
    assert_eq!(s["poisson_ratio"], 0.392);
    // Drawing raises the modulus above serial 3802's annealed 10,600 and Table 67's 10,300 kips/in2,
    // and the printed shear modulus 3,900 kips/in2 agrees with E/(2(1+nu)) to under 1 percent.
    assert!(s["young_modulus"] > 10_600.0 * MPA_PER_KSI * 1e6);
    let shear = s["young_modulus"] / (2.0 * (1.0 + s["poisson_ratio"]));
    assert!((shear / (3_900.0 * MPA_PER_KSI * 1e6) - 1.0).abs() < 0.01);
    // Lighter than lead (11340) and gold (19300), heavier than copper (8960).
    assert!((9_000.0..11_000.0).contains(&s["density"]));
    assert!((0.25..0.5).contains(&s["poisson_ratio"]));
    // Table 67's thermal columns at the stated factors: 18.9e-6 per deg C, 0.056 cal/g/deg C at
    // 4.184 J/cal, 4.08 W/cm/deg C.
    assert_eq!(s["linear_thermal_expansion_coefficient"], 18.9e-6);
    assert!((s["specific_heat_capacity"] / 4184.0 - 0.056).abs() < 1e-9);
    assert_eq!(s["thermal_conductivity"], 408.0);
    // Dulong-Petit with Table 67's atomic weight 107.880: 3R/M = 231.2 J/kg/K.
    let dulong_petit = 3.0 * 8.314_462_618 / 0.107_88;
    assert!((s["specific_heat_capacity"] / dulong_petit - 1.0).abs() < 0.05);

    let strength = scalars("silver-drawn-hard-uts-nbs-c101", "nbscircular101.pdf");
    let names: Vec<&str> = strength.keys().map(String::as_str).collect();
    // One tensile strength and nothing else: no yield of drawn silver is sourced.
    assert_eq!(names, ["ultimate_tensile_strength_drawn_hard"]);
    let uts = strength["ultimate_tensile_strength_drawn_hard"];
    // Printed 36.0 kg/mm2 = 51 200 lbs./in.2 (the same value rounded to three figures).
    assert!((uts / 1e6 / MPA_PER_KGF_MM2 / 36.0 - 1.0).abs() < 1e-12);
    assert!((uts / 1e6 / MPA_PER_KSI / 51.2 - 1.0).abs() < 2e-3);
    // Drawn silver is stronger than annealed rod (Table 60 serial 3798: 21.5 kips/in2) and far
    // weaker than steel music wire. A concert-pitch string at 4 L^2 f^2 rho, even A4 (440 Hz) at a
    // short 0.4 m speaking length, already carries more.
    assert!(uts > 21.5 * MPA_PER_KSI * 1e6 && uts < 1e9);
    assert!(4.0 * 0.4_f64.powi(2) * 440.0_f64.powi(2) * s["density"] > uts);
}
