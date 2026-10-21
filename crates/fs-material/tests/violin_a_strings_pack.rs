//! Gates for the six violin A strings of Pickering 1990 (CASJ 1(5):23): each printed diameter in millimetres is the printed diameter in
//! inches times 25.4 and each printed tension in newtons is the printed pounds times 4.44822 (both to 0.6 percent, the table's
//! rounding), the tension rises with the diameter inside each core type (heavy over light), and the stretched A string's linear mass
//! derived from T/(2 L f)^2 for the 328 mm bed at 440 Hz is that of a string of its diameter and a density between that of gut and that of
//! steel.

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

#[test]
fn violin_a_strings_pack_passes_unit_and_density_gates() {
    let dir = "violin-a-strings-by-core-pickering-casj1990";
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.contains("\nlicense\tPublished-Scientific-Facts-Attribution-Retained\n"));
    assert!(manifest.contains("Catgut Acoustical Society Journal 1(5, Series II)"));
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
    let keys = [
        "steel_nickel",
        "steel_aluminum",
        "perlon_aluminum_heavy",
        "perlon_aluminum_light",
        "gut_aluminum_heavy",
        "gut_aluminum_light",
    ];
    for k in keys {
        let (d, d_in) = (
            s[&format!("diameter_{k}")],
            s[&format!("diameter_printed_in_inches_{k}")],
        );
        assert!(
            (d - d_in * 0.0254).abs() / d < 0.006,
            "{k}: diameter {d} against {d_in} in"
        );
        let (t, t_lb) = (
            s[&format!("tension_{k}")],
            s[&format!("tension_printed_in_pounds_force_{k}")],
        );
        assert!(
            (t - t_lb * 4.44822).abs() / t < 0.006,
            "{k}: tension {t} against {t_lb} lbf"
        );
        // Linear mass from T = mu (2 L f)^2, L = 0.328 m, f = 440 Hz, and the density of a solid string of that mass and diameter.
        let mu = t / (2.0 * 0.328 * 440.0_f64).powi(2);
        let rho = mu / (std::f64::consts::PI * (d / 2.0).powi(2));
        assert!(
            (1000.0..9000.0).contains(&rho),
            "{k}: equivalent density {rho} kg/m^3"
        );
    }
    for core in ["perlon_aluminum", "gut_aluminum"] {
        let (h, l) = (format!("{core}_heavy"), format!("{core}_light"));
        assert!(s[&format!("tension_{h}")] > s[&format!("tension_{l}")]);
        assert!(s[&format!("diameter_{h}")] > s[&format!("diameter_{l}")]);
    }
}
