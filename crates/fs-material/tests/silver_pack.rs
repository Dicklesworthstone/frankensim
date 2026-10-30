//! Physicality gates for the pure-silver matdb pack (bead `frankensim-irp-v5-instrument-realism-pexfp.192`, plan item M8-W2b: silver
//! windings), from Smith and Fickett 1995, J. Res. NIST 100(2):119-171 (a US Government work, public domain).
//!
//! The gates encode what cannot be wrong whichever silver was measured: positivity; the printed 0 K polycrystal Young's modulus lying
//! between the Reuss and Voigt bounds computed from the printed single-crystal constants C11, C12 and C44; silver stiffening on cooling;
//! the annealed moduli below the handbook value (annealed silver is the softer state); Poisson's ratio in the metals' band; density of a
//! heavy fcc metal; schema hygiene (every scalar has an uncertainty and a validity row).

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

#[test]
fn silver_pack_passes_physicality_gates() {
    let dir = "silver-pure-smith-fickett-nist1995";
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.contains("\nlicense\tWork-of-the-US-Government-Public-Use-Permitted\n"));
    assert!(manifest.contains("https://doi.org/10.6028/jres.100.012"));
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    let mut s: BTreeMap<String, f64> = BTreeMap::new();
    let (mut ids, mut unc, mut val) = (Vec::new(), Vec::new(), Vec::new());
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        match c[0] {
            "scalar" => {
                let v: f64 = c[4].parse().expect("value");
                assert!(v.is_finite() && v > 0.0, "{} must be positive", c[3]);
                s.insert(c[3].to_string(), v);
                ids.push(c[1].to_string());
            }
            "curve" => {
                let mut previous = f64::NEG_INFINITY;
                for knot in c[7].split(',') {
                    let (t, v) = knot.split_once(':').expect("temperature/value");
                    let t: f64 = t.parse().expect("temperature");
                    let v: f64 = v.parse().expect("value");
                    assert!(t.is_finite() && t > previous);
                    assert!(v.is_finite() && v > 0.0);
                    previous = t;
                }
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
            "{id}: uncertainty and validity rows"
        );
    }

    // The 0 K polycrystal modulus against the Voigt and Reuss bounds of the printed single-crystal constants (cubic).
    let (c11, c12, c44) = (
        s["single_crystal_c11_at_0_k"],
        s["single_crystal_c12_at_0_k"],
        s["single_crystal_c44_at_0_k"],
    );
    let bulk = (c11 + 2.0 * c12) / 3.0;
    let g_voigt = (c11 - c12 + 3.0 * c44) / 5.0;
    let g_reuss = 5.0 * (c11 - c12) * c44 / (4.0 * c44 + 3.0 * (c11 - c12));
    let young = |g: f64| 9.0 * bulk * g / (3.0 * bulk + g);
    let (e_reuss, e_voigt) = (young(g_reuss), young(g_voigt));
    let e0 = s["young_modulus_at_0_k"];
    assert!(e_reuss < e0 && e0 < e_voigt, "{e_reuss} < {e0} < {e_voigt}");

    let e293 = s["young_modulus"];
    assert!(e293 < e0, "silver stiffens on cooling");
    assert!(
        s["young_modulus_annealed_as_received"]
            < s["young_modulus_annealed_700c_after_5_percent_overstrain"]
    );
    assert!(s["young_modulus_annealed_700c_after_5_percent_overstrain"] < e293);
    assert!((0.3..0.4).contains(&s["poisson_ratio"]));
    assert!((10400.0..10700.0).contains(&s["density"]));
}
