//! Gates for the string-section form of the Heptonstall, Waller and Robertson music-wire data (bead
//! `frankensim-irp-v5-instrument-realism-pexfp.363`, U2): canonical names over the same LIGO P1400229 facts.
//!
//! The gates encode what cannot be wrong: every value shared with the original pack
//! `steel-music-wire-heptonstall-ligo-p1400229` is identical to it (density, untreated modulus, the 465 um untreated loss
//! tangent and its 8-3000 Hz band); every breaking stress lies in the band of cold-drawn high-carbon steel wire and
//! exceeds the stress a steel string must carry at concert pitch at the treble end of the scale; schema hygiene (every
//! scalar has an uncertainty and a validity row, every row is admitted at the same absolute temperature).

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

fn scalars(dir: &str) -> (BTreeMap<String, f64>, Vec<Vec<String>>) {
    let text = std::fs::read_to_string(repo_root().join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
        .expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let (mut s, mut validity) = (BTreeMap::new(), Vec::new());
    let (mut ids, mut unc, mut val) = (Vec::new(), Vec::new(), Vec::new());
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        match c[0] {
            "scalar" => {
                let v: f64 = c[4].parse().expect("value");
                assert!(v.is_finite() && v > 0.0, "{dir}: {} must be positive", c[3]);
                let si = match c[5] { "GPa" => v * 1e9, "MPa" => v * 1e6,
                    "Pa" | "kg/m^3" | "1" | "K^-1" | "J/kg/K" | "W/m/K" => v,
                    unit => panic!("{dir}: {} unit {unit}", c[3]) };
                assert!(s.insert(c[3].to_string(), si).is_none(), "{dir}: {} is duplicated", c[3]);
                ids.push(c[1].to_string());
            }
            "uncertainty" => unc.push(c[1].to_string()),
            "validity" => { val.push(c[1].to_string()); validity.push(c.iter().map(|x| x.to_string()).collect()); }
            "observation" => {}
            other => panic!("unknown row kind {other}"),
        }
    }
    for id in &ids {
        assert!(unc.contains(id) && val.contains(id), "{dir}: {id} needs uncertainty and validity rows");
    }
    (s, validity)
}

#[test]
fn steel_string_pack_matches_its_source_and_passes_physicality_gates() {
    let dir = "steel-music-wire-string-heptonstall-ligo-p1400229";
    let manifest = std::fs::read_to_string(repo_root().join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
        .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v3\n"));
    assert!(manifest.contains("\nlicense\tPublished-Scientific-Facts-Attribution-Retained\n"));
    assert!(manifest.contains("scientific facts") && manifest.contains("arXiv:1506.03778"));
    assert!(manifest.contains("\naxis\tT\tT\tabsolute-temperature\n"));
    assert!(manifest.contains("\naxis\tfrequency\tfrequency\tcyclic-frequency\n"));
    let (s, validity) = scalars(dir);
    let (original, _) = scalars("steel-music-wire-heptonstall-ligo-p1400229");
    assert_eq!(s["density"], original["density_literature_value_used"]);
    assert_eq!(s["young_modulus"], original["young_modulus_untreated"]);
    assert_eq!(s["bending_loss_tangent"], original["material_loss_angle_465um_untreated"]);
    for row in &validity {
        match row[2].as_str() {
            "T" => assert_eq!(&row[3..], ["293.15", "293.15", "K"], "{}", row[1]),
            "frequency" => {
                assert_eq!(row[1], "loss-tangent-465um");
                assert_eq!(&row[3..], ["8", "3000", "Hz"]);
            }
            axis => panic!("unexpected validity axis {axis}"),
        }
    }
    let breaking: Vec<(&String, &f64)> = s.iter().filter(|(k, _)| k.starts_with("breaking_stress_")).collect();
    assert_eq!(breaking.len(), 13);
    // A plain string at fixed pitch carries 4 L^2 f^2 rho whatever its diameter. At C8 on the
    // Model D (RT-0425: 0.052 m, about 4186 Hz) steel needs about 1.5 GPa; every measured
    // music wire is stronger, and none exceeds the cold-drawn high-carbon steel band.
    let concert_c8 = 4.0 * 0.052_f64.powi(2) * 4186.0_f64.powi(2) * s["density"];
    for (name, stress) in breaking {
        assert!(*stress > concert_c8 && *stress < 3.5e9, "{name}: {stress}");
    }
    assert!((190e9..220e9).contains(&s["young_modulus"]));
    assert!((7700.0..7900.0).contains(&s["density"]));
    // Table 3's thermal inputs reproduce the paper's own printed thermoelastic
    // corner for wire 1 (465 um): Table 4 lists 129 Hz, eq. 6 f_c = 2.16 k/(rho C d^2),
    // which is 4 xi_1^2/(2 pi) k/(rho C d^2) of the circular-section series.
    let xi = fs_material::visco::circular_ted::circular_thermal_modes()[0].0;
    let corner = 4.0 * xi * xi / (2.0 * std::f64::consts::PI) * s["thermal_conductivity"]
        / (s["density"] * s["specific_heat_capacity"] * 465e-6_f64.powi(2));
    assert!((corner / 129.0 - 1.0).abs() < 0.01, "{corner}");
    // Dulong-Petit for iron (55.845 g/mol): 3R/M = 446.7 J/kg/K; steel is within 10 %.
    assert!((s["specific_heat_capacity"] / (3.0 * 8.314_462_618 / 0.055_845) - 1.0).abs() < 0.10);
    assert!((1.0e-5..1.3e-5).contains(&s["linear_thermal_expansion_coefficient"]));
}
