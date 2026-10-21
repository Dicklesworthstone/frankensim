//! Physicality gates for the tonewood matdb packs (bead `frankensim-irp-v5-instrument-realism-pexfp.193`, plan item M8-W2c):
//! sycamore maple from Sonderegger et al. 2013 (full orthotropic set) and Spycher et al. 2008 (L and R flexural moduli), Indian
//! rosewood from Sprossmann et al. 2017 (CC BY 4.0), the four single-specimen axial rows of Bremaud et al. 2012, and the Dalbergia
//! nigra DMA rows of de Barros and Politano 2019 (CC BY 4.0).
//!
//! The gates encode what the sources themselves state or what cannot be wrong for any wood: positivity; E_L > E_R > E_T for every
//! test mode; ultrasonic moduli above the static ones (Sonderegger p. 97); the paper's own L:R:T modulus ratio 9.8:1.7:1 and shear
//! ratios 4.5 and 2.9 (whose LR/LT labels the paper's text and Table 4 disagree on); reciprocity nu_ij/E_i = nu_ji/E_j, which only
//! holds after the paper's reversed Poisson subscripts are converted to the FPL order (the swapped reading fails by a factor of
//! three); positive definiteness of the Poisson set; Spycher's printed sound speeds equal to sqrt(E/rho) of the transcribed moduli;
//! E' = (E'/gamma) x gamma for the derived Bremaud rows; storage modulus rising and tan delta falling from 0.025 to 50 Hz in the DMA
//! rows. Schema hygiene: every scalar has an uncertainty and a validity row.

use std::collections::BTreeMap;

const FACTS: &str = "Published-Scientific-Facts-Attribution-Retained";

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

/// Claim id -> value, after asserting the manifest licence and schema hygiene.
fn load_pack(dir: &str, license: &str) -> BTreeMap<String, f64> {
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(
        manifest.contains(&format!("\nlicense\t{license}\n")),
        "{dir}: licence {license}"
    );
    assert!(manifest.contains("doi:10."), "{dir}: DOI required");
    if license == FACTS {
        assert!(
            manifest.contains("facts"),
            "{dir}: facts-only terms must say the pack keeps facts"
        );
    }
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let (mut scalars, mut unc, mut val) = (BTreeMap::new(), Vec::new(), Vec::new());
    for line in text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "scalar" => {
                let value: f64 = cols[4].parse().expect("scalar value");
                assert!(
                    value.is_finite() && value > 0.0,
                    "{dir}: {} must be positive",
                    cols[1]
                );
                assert!(
                    scalars.insert(cols[1].to_string(), value).is_none(),
                    "{dir}: duplicate claim {}",
                    cols[1]
                );
            }
            "uncertainty" => unc.push(cols[1].to_string()),
            "validity" => val.push(cols[1].to_string()),
            "observation" => {}
            other => panic!("unknown row kind {other}"),
        }
    }
    for id in scalars.keys() {
        assert!(
            unc.contains(id),
            "{dir}: scalar {id} lacks an uncertainty row"
        );
        assert!(val.contains(id), "{dir}: scalar {id} lacks a validity row");
    }
    scalars
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / (0.5 * (a + b))
}

#[test]
fn sycamore_maple_sonderegger_full_set_passes_physicality_gates() {
    let p = load_pack("maple-sycamore-sonderegger2013", FACTS);
    let e = |dir: &str, mode: &str| p[&format!("young-modulus-{dir}-{mode}")];

    for mode in ["tension", "compression", "ultrasonic"] {
        let (l, r, t) = (
            e("longitudinal", mode),
            e("radial", mode),
            e("tangential", mode),
        );
        assert!(l > r && r > t, "{mode}: E_L > E_R > E_T, got {l} {r} {t}");
    }
    // Ultrasonic moduli exceed the static ones in every direction (p. 97: "up to 2.5 times higher").
    for dir in ["longitudinal", "radial", "tangential"] {
        let us = e(dir, "ultrasonic");
        for mode in ["tension", "compression"] {
            let st = e(dir, mode);
            assert!(
                us > st && us <= 2.6 * st,
                "{dir}: ultrasonic {us} vs {mode} {st}"
            );
        }
    }
    // The paper's text: L:R:T = 9.8:1.7:1 from the ultrasonic moduli of Table 4.
    let t_us = e("tangential", "ultrasonic");
    assert!((e("longitudinal", "ultrasonic") / t_us - 9.8).abs() < 0.05);
    assert!((e("radial", "ultrasonic") / t_us - 1.7).abs() < 0.05);
    // Shear: G_RT is the smallest; the two other planes are 4.5x and 2.9x G_RT as the text says, whichever label is right.
    let (glr, glt, grt) = (
        p["shear-modulus-lr-ultrasonic"],
        p["shear-modulus-lt-ultrasonic"],
        p["shear-modulus-rt-ultrasonic"],
    );
    assert!(grt < glr.min(glt));
    let (hi, lo) = (glr.max(glt) / grt, glr.min(glt) / grt);
    assert!(
        (hi - 4.5).abs() < 0.15 && (lo - 2.9).abs() < 0.15,
        "{hi} {lo}"
    );

    // Reciprocity with the tension moduli: nu_RT/E_R = nu_TR/E_T and nu_LR/E_L = nu_RL/E_R.
    let nu = |k: &str| p[&format!("nu-{k}")];
    let (el, er, et) = (
        e("longitudinal", "tension"),
        e("radial", "tension"),
        e("tangential", "tension"),
    );
    let rt = rel(nu("rt-tension") / er, nu("tr-tension") / et);
    assert!(rt < 0.05, "R-T reciprocity off by {rt}");
    let lr = rel(nu("lr-tension") / el, nu("rl-tension") / er);
    assert!(lr < 0.15, "L-R reciprocity off by {lr}");
    // Falsifier: reading the printed subscripts in the FPL order (unconverted) breaks R-T reciprocity by about a factor of three.
    let swapped = (nu("rt-tension") / et) / (nu("tr-tension") / er);
    assert!(
        swapped > 2.5,
        "unconverted subscripts would satisfy reciprocity: {swapped}"
    );

    // Positive definiteness of the compression set and of the all-climate means (each product < 1 and the determinant > 0).
    for suffix in ["compression", "mean"] {
        let n = |k: &str| nu(&format!("{k}-{suffix}"));
        let (lr_, rl, rt_, tr, lt, tl) = (n("lr"), n("rl"), n("rt"), n("tr"), n("lt"), n("tl"));
        for prod in [lr_ * rl, rt_ * tr, lt * tl] {
            assert!(prod < 1.0, "{suffix}: {prod}");
        }
        let det = 1.0 - lr_ * rl - rt_ * tr - tl * lt - 2.0 * rl * tr * lt;
        assert!(det > 0.0, "{suffix}: Poisson determinant {det}");
    }
    assert!(p["density-ovendry"] < p["density-normal-climate"]);
}

#[test]
fn sycamore_maple_spycher_bars_pass_physicality_gates() {
    let p = load_pack("maple-sycamore-spycher2008", FACTS);
    // Table 3 prints the sound speed of each group (axial, radial), m/s; it must equal sqrt(E/rho) of the transcribed rows.
    for (g, c_ax, c_rad) in [
        ("a", 3533.0, 1348.0),
        ("b", 3894.0, 1662.0),
        ("control", 4180.0, 1703.0),
    ] {
        let rho = p[&format!("density-{g}")];
        let el = p[&format!("dynamic-young-modulus-longitudinal-{g}")] * 1.0e6;
        let er = p[&format!("dynamic-young-modulus-radial-{g}")] * 1.0e6;
        assert!(el > er, "{g}: E_L > E_R");
        assert!((500.0..700.0).contains(&rho), "{g}: density {rho}");
        assert!(rel((el / rho).sqrt(), c_ax) < 0.015, "{g}: axial speed");
        assert!(rel((er / rho).sqrt(), c_rad) < 0.015, "{g}: radial speed");
    }
    assert!((p["loss-factor-longitudinal"] - 0.043 / std::f64::consts::PI).abs() < 1.0e-6);
}

#[test]
fn indian_rosewood_sprossmann_bars_pass_physicality_gates() {
    let p = load_pack("rosewood-indian-sprossmann2017", "CC-BY-4.0");
    let (es, ed, g) = (
        p["young-modulus-longitudinal-static-bending"],
        p["dynamic-young-modulus-longitudinal"],
        p["dynamic-shear-modulus-lr"],
    );
    assert!(rel(es, ed) < 0.10, "static {es} and dynamic {ed} E_L agree");
    let anisotropy = ed / g;
    assert!((5.0..20.0).contains(&anisotropy), "E_L/G_LR {anisotropy}");
    assert!((0.001..0.05).contains(&p["loss-factor-longitudinal"]));
    assert!((600.0..1000.0).contains(&p["density"]));
}

#[test]
fn bremaud_axial_rows_pass_physicality_gates() {
    let mut rows = BTreeMap::new();
    for dir in [
        "rosewood-brazilian-bremaud2012",
        "rosewood-indian-bremaud2012",
        "maple-sycamore-bremaud2012",
        "cedar-western-red-bremaud2012",
    ] {
        let p = load_pack(dir, FACTS);
        let derived = p["specific-dynamic-modulus-longitudinal"] * p["specific-gravity"] * 1000.0;
        assert!(
            rel(p["dynamic-young-modulus-longitudinal"], derived) < 1.0e-9,
            "{dir}: E' = (E'/gamma) x gamma"
        );
        assert!((0.001..0.05).contains(&p["loss-factor-longitudinal"]));
        rows.insert(dir, p);
    }
    // The standard relation the paper tests: the low-specific-modulus maple carries the highest damping, the Dalbergia rows the lowest.
    let get = |d: &str, k: &str| rows[d][k];
    let maple = "maple-sycamore-bremaud2012";
    for d in rows.keys().filter(|d| **d != maple) {
        assert!(get(maple, "loss-factor-longitudinal") > get(d, "loss-factor-longitudinal"));
        assert!(
            get(maple, "specific-dynamic-modulus-longitudinal")
                < get(d, "specific-dynamic-modulus-longitudinal")
        );
    }
    // Same species, independent source: the single Bremaud D. latifolia specimen against Sprossmann's 20-bar mean.
    let spr = load_pack("rosewood-indian-sprossmann2017", "CC-BY-4.0");
    let one = get(
        "rosewood-indian-bremaud2012",
        "dynamic-young-modulus-longitudinal",
    );
    assert!(rel(one, spr["dynamic-young-modulus-longitudinal"]) < 0.15);
}

#[test]
fn dalbergia_nigra_dma_rows_pass_physicality_gates() {
    let p = load_pack("rosewood-brazilian-barros2019", "CC-BY-4.0");
    for s in ["1", "2"] {
        let k = |name: &str| p[&format!("{name}-{s}")];
        assert!(k("storage-modulus-longitudinal-50hz") > k("storage-modulus-longitudinal-0p025hz"));
        assert!(k("tan-delta-longitudinal-50hz") < k("tan-delta-longitudinal-0p025hz"));
        assert!((1000.0..1300.0).contains(&k("density")));
    }
}
