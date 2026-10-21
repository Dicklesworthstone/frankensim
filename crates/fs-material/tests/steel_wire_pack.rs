//! Physicality gates for the steel music-wire matdb pack (bead `frankensim-irp-v5-instrument-realism-pexfp.192`, plan item M8-W2b: wire
//! cores of solid-steel and wound strings), from Heptonstall, Waller and Robertson (LIGO P1400229, arXiv:1506.03778): Young's modulus from
//! 50 tensile tests and the material loss angle of four wires from ringdowns in vacuum.
//!
//! The gates encode what cannot be wrong whichever wire was measured: positivity; the printed Young's moduli of untreated and cryogenically
//! treated wire consistent within their printed half-widths (the paper's own conclusion); the printed mean loss angle equal to the mean of
//! the four wires' values within its printed half-width; the loss angles in the band of a high-carbon steel (1e-5 to 1e-3); density of
//! steel; schema hygiene (every scalar has an uncertainty and a validity row).

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

#[test]
fn steel_wire_pack_passes_physicality_gates() {
    let dir = "steel-music-wire-heptonstall-ligo-p1400229";
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.contains("\nlicense\tPublished-Scientific-Facts-Attribution-Retained\n"));
    assert!(manifest.contains("arXiv:1506.03778"));
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

    // Untreated and cryogenically treated moduli agree within their printed half-widths ("both consistent").
    let (e_u, e_c) = (
        s["young_modulus_untreated"],
        s["young_modulus_cryo_treated"],
    );
    let tol = s["young_modulus_untreated_printed_half_width"]
        + s["young_modulus_cryo_treated_printed_half_width"];
    assert!((e_u - e_c).abs() <= tol, "{e_u} against {e_c}");
    assert!((1.9e11..2.3e11).contains(&e_u), "steel modulus band");
    assert!((7700.0..7900.0).contains(&s["density_literature_value_used"]));

    // The printed mean loss angle is the mean of the four wires' (assuming no cryogenic effect), within its half-width.
    let four = [
        "material_loss_angle_465um_untreated",
        "material_loss_angle_457um_cryo_treated",
        "material_loss_angle_197um_untreated",
        "material_loss_angle_197um_cryo_treated",
    ];
    let mean = four.iter().map(|k| s[*k]).sum::<f64>() / 4.0;
    let (printed, half) = (
        s["material_loss_angle_mean_of_all_wires_printed"],
        s["material_loss_angle_mean_of_all_wires_printed_half_width"],
    );
    assert!(
        (mean - printed).abs() <= half,
        "mean of four {mean} against printed {printed}"
    );
    for k in four {
        assert!((1e-5..1e-3).contains(&s[k]), "{k} outside the steel band");
    }
}
