//! Physicality gates for the pure-gold matdb pack (bead `frankensim-irp-v5-instrument-realism-pexfp.363`, U2: gold piano strings as a
//! data swap), from NBS Circular C447 (1943), Table 60 serials 3641-3645 and Table 67 (a US Government work, public domain).
//!
//! The gates encode what cannot be wrong whichever gold was tested: positivity; the printed kips/in2 and g/cm3 values reproduced from
//! the SI rows at the stated conversion factors; for one annealed rod, proportional limit below the 0.01 percent yield below the
//! tensile strength; cold work raising the wire's tensile strength above every annealed or cast value; an endurance limit below its
//! sheet's tensile strength; a heavy fcc metal's density; Poisson's ratio in the metals' band; schema hygiene (every scalar has an
//! uncertainty and a validity row, and every row is admitted at the same absolute temperature).

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

const MPA_PER_KSI: f64 = 6.894757293168;

#[test]
fn gold_pack_passes_physicality_gates() {
    let dir = "gold-pure-nbs-c447";
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v3\n"));
    assert!(manifest.contains("\nlicense\tWork-of-the-US-Government-Public-Use-Permitted\n"));
    assert!(manifest.contains("nbscircular447.pdf"));
    assert!(manifest.contains("\naxis\tT\tT\tabsolute-temperature\n"));
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let mut s: BTreeMap<String, f64> = BTreeMap::new();
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
    // The printed values, recovered from the SI rows at the stated factors.
    let printed = |name: &str| s[name] / 1e6 / MPA_PER_KSI;
    for (name, ksi) in [
        ("young_modulus", 10_300.0),
        ("proportional_limit_annealed_rod", 2.9),
        ("yield_strength_0p01pct_offset_annealed_rod", 8.1),
        ("ultimate_tensile_strength_annealed_rod", 17.9),
        ("shear_modulus_printed_range_low", 4_100.0),
        ("shear_modulus_printed_range_high", 5_600.0),
        ("ultimate_tensile_strength_cast_rod", 16.5),
        ("ultimate_tensile_strength_cold_rolled_wire", 31.4),
        ("ultimate_tensile_strength_annealed_sheet", 17.7),
        ("endurance_limit_1e7_cycles_annealed_sheet", 4.6),
    ] {
        assert!((printed(name) / ksi - 1.0).abs() < 5e-6, "{name}: {} ksi printed {ksi}", printed(name));
    }
    assert_eq!(s["density"], 19_300.0);
    // One annealed rod: proportional limit < 0.01 % yield < tensile strength.
    assert!(s["proportional_limit_annealed_rod"] < s["yield_strength_0p01pct_offset_annealed_rod"]);
    assert!(s["yield_strength_0p01pct_offset_annealed_rod"] < s["ultimate_tensile_strength_annealed_rod"]);
    // Cold work (60 % reduction) strengthens: the wire beats every annealed or cast value.
    for soft in ["ultimate_tensile_strength_annealed_rod", "ultimate_tensile_strength_cast_rod",
        "ultimate_tensile_strength_annealed_sheet"] {
        assert!(s["ultimate_tensile_strength_cold_rolled_wire"] > 1.5 * s[soft], "{soft}");
    }
    assert!(s["endurance_limit_1e7_cycles_annealed_sheet"] < s["ultimate_tensile_strength_annealed_sheet"]);
    assert!(s["shear_modulus_printed_range_low"] < s["shear_modulus_printed_range_high"]);
    // A heavy fcc metal, denser than silver (10492 kg/m3) and lead (11340 kg/m3).
    assert!((18_000.0..21_000.0).contains(&s["density"]));
    assert!((0.25..0.5).contains(&s["poisson_ratio"]));
    // Table 67's thermal columns, recovered at the stated factors: 14.2e-6 per deg C,
    // 0.031 cal/g/deg C at 4.184 J/cal, 2.96 W/cm/deg C.
    assert_eq!(s["linear_thermal_expansion_coefficient"], 14.2e-6);
    assert!((s["specific_heat_capacity"] / 4184.0 - 0.031).abs() < 1e-9);
    assert_eq!(s["thermal_conductivity"], 296.0);
    // Dulong-Petit: a solid element's heat capacity near room temperature is
    // close to 3R per mole; with Table 67's atomic weight 197.2, 3R/M = 126.5 J/kg/K.
    let dulong_petit = 3.0 * 8.314_462_618 / 0.1972;
    assert!((s["specific_heat_capacity"] / dulong_petit - 1.0).abs() < 0.05);
}
