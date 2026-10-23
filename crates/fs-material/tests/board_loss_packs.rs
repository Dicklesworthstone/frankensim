//! Physicality gates for the soundboard loss seeds of bead
//! `frankensim-irp-v5-instrument-realism-pexfp.373`: the Sitka spruce
//! along-grain loss over moisture content (Qiu et al. 2026) and the Norway
//! spruce bar damping over frequency and moisture (Merhar and Moutou Pitti
//! 2025), both CC BY 4.0 (Forests); and the along-grain loss of five Pinus
//! species (Bremaud et al. 2012, facts only), the composition-proxy bracket of
//! the sugar pine ribs.
//!
//! The gates encode what cannot be wrong whichever wood was tested: the printed
//! values reproduced exactly; typed v3 axes; loss in wood's small-damping
//! regime; Qiu's loss rising with moisture content as the source states;
//! Merhar's 24 cells complete, each at one mode frequency and one equilibrium
//! moisture content, with the three modes' frequencies ordered and the higher
//! bending modes more damped than the first at every moisture content (the
//! rise the source reports, partly shear loss in these deep bars).

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).parent().and_then(std::path::Path::parent)
        .expect("root").to_path_buf()
}
fn read(dir: &str, file: &str) -> String {
    std::fs::read_to_string(repo_root().join(format!("data/matdb/seed-v1/{dir}/{file}"))).expect(file)
}
fn typed_manifest(dir: &str) {
    typed_manifest_licensed(dir, "CC-BY-4.0");
}
fn typed_manifest_licensed(dir: &str, license: &str) {
    let manifest = read(dir, "manifest.tsv");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v3\n"), "{dir}");
    assert!(manifest.contains(&format!("\npack_id\t{dir}\n")), "{dir}");
    assert!(manifest.contains(&format!("\nlicense\t{license}\n")), "{dir}");
    assert!(manifest.contains("\naxis\tmoisture_content\tmoisture_content\tdimensional\n"), "{dir}");
    assert!(manifest.contains("\naxis\tfrequency\tfrequency\tcyclic-frequency\n"), "{dir}");
}

#[test]
fn sitka_loss_moisture_curve_reproduces_the_printed_means_and_rises_with_moisture() {
    let dir = "spruce-sitka-loss-moisture-qiu2026";
    typed_manifest(dir);
    let text = read(dir, "properties.tsv");
    let rows: Vec<Vec<&str>> = text.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    let kinds: Vec<&str> = rows.iter().map(|r| r[0]).collect();
    assert_eq!(kinds, ["observation", "curve", "uncertainty", "validity", "validity"]);
    let curve = &rows[1];
    assert_eq!(&curve[3..7], ["loss_factor_longitudinal", "moisture_content", "%", "1"]);
    assert_eq!(curve[8], "linear");
    let knots: Vec<(f64, f64)> = curve[7].split(',').map(|k| {
        let (x, y) = k.split_once(':').unwrap();
        (x.parse().unwrap(), y.parse().unwrap())
    }).collect();
    // Text pp.11-12: the six group means at 2, 4, ..., 12 percent.
    let printed = [(2.0, 0.01012), (4.0, 0.01123), (6.0, 0.01217), (8.0, 0.01297), (10.0, 0.01366), (12.0, 0.01424)];
    assert_eq!(knots, printed);
    assert!(knots.windows(2).all(|w| w[1].1 > w[0].1), "loss rises with moisture content");
    assert!(knots.iter().all(|k| (0.004..0.02).contains(&k.1)), "wood's small-damping regime");
    assert_eq!(&rows[3][2..], ["moisture_content", "2", "12", "%"]);
    // The first mode's mean frequency fell from 68.34 Hz at 2 % to 62.36 Hz at 12 %.
    assert_eq!(&rows[4][2..], ["frequency", "62.36", "68.34", "Hz"]);
}

#[test]
fn norway_spruce_bar_damping_reproduces_table_2_and_rises_from_the_first_mode() {
    let dir = "spruce-norway-loss-frequency-moisture-merhar2025";
    typed_manifest(dir);
    let text = read(dir, "properties.tsv");
    // (plane, mode, emc) -> (zeta, frequency)
    let mut cells: BTreeMap<(String, u32, String), (f64, f64)> = BTreeMap::new();
    let mut zeta: BTreeMap<String, f64> = BTreeMap::new();
    let mut freq: BTreeMap<String, f64> = BTreeMap::new();
    let mut emc: BTreeMap<String, String> = BTreeMap::new();
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        match c[0] {
            "observation" => assert!(c[4].contains("shear") && c[4].contains("Norway spruce, not Sitka")),
            "scalar" => {
                assert_eq!((c[5], c[6]), ("1", "constant"));
                assert!(c[3] == "flexural_damping_ratio_lr_plane" || c[3] == "flexural_damping_ratio_lt_plane");
                zeta.insert(c[1].to_string(), c[4].parse().unwrap());
            }
            "uncertainty" => assert_eq!(&c[2..], ["unstated", "-", "-", "-", "-"]),
            "validity" if c[2] == "frequency" => { assert_eq!(c[3], c[4]); freq.insert(c[1].to_string(), c[3].parse().unwrap()); }
            "validity" if c[2] == "moisture_content" => { assert_eq!(c[3], c[4]); emc.insert(c[1].to_string(), c[3].to_string()); }
            other => panic!("unexpected row {other}"),
        }
    }
    for (id, z) in &zeta {
        let parts: Vec<&str> = id.split('-').collect(); // damping-ratio-<plane>-mode<k>-rh<rh>
        let mode: u32 = parts[3].trim_start_matches("mode").parse().unwrap();
        cells.insert((parts[2].to_string(), mode, emc[id].clone()), (*z, freq[id]));
    }
    assert_eq!(cells.len(), 24, "3 modes x 4 moisture contents x 2 planes");
    // Spot checks against Table 2 (untreated): LR mode 1 at 44 % RH, LT mode 3 at 88 % RH.
    assert_eq!(cells[&("lr".into(), 1, "9.9".into())], (0.0050, 2673.0));
    assert_eq!(cells[&("lt".into(), 3, "19.1".into())], (0.0097, 8291.0));
    for plane in ["lr", "lt"] {
        for mc in ["7.2", "9.9", "13.3", "19.1"] {
            let at = |mode: u32| cells[&(plane.to_string(), mode, mc.to_string())];
            assert!(at(1).1 < at(2).1 && at(2).1 < at(3).1, "{plane} {mc}: mode order");
            assert!(at(2).0 > at(1).0 && at(3).0 > at(1).0, "{plane} {mc}: higher modes more damped");
            for mode in 1..=3 { assert!((0.002..0.02).contains(&at(mode).0), "{plane} {mc} mode {mode}"); }
        }
    }
}

#[test]
fn pinus_along_grain_loss_reproduces_table_iii_and_is_flagged_as_a_composition_proxy() {
    let dir = "pine-loss-longitudinal-bremaud2012";
    typed_manifest_licensed(dir, "Published-Scientific-Facts-Attribution-Retained");
    let text = read(dir, "properties.tsv");
    // species -> (tan delta, spread, f_lo, f_hi, EMC percent)
    let mut cells: BTreeMap<String, (f64, f64, f64, f64, f64)> = BTreeMap::new();
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        let species = || c[1].rsplit("pinus-").next().unwrap().to_string();
        match c[0] {
            "observation" => {
                assert!(c[4].starts_with("composition-proxy: a Pinus species, NOT sugar pine (Pinus lambertiana)"), "{}", c[1]);
                assert!(c[3].contains("free-free beams, first mode") && c[3].contains("200 and 600 Hz"), "{}", c[1]);
            }
            "scalar" => {
                assert_eq!((c[3], c[5], c[6]), ("loss_factor_longitudinal", "1", "constant"));
                cells.entry(species()).or_default().0 = c[4].parse().unwrap();
            }
            "uncertainty" => {
                assert_eq!((c[2], c[4], c[5], c[6]), ("absolute", "1", "0.68", "1"));
                cells.entry(species()).or_default().1 = c[3].parse().unwrap();
            }
            "validity" if c[2] == "frequency" => {
                let cell = cells.entry(species()).or_default();
                (cell.2, cell.3) = (c[3].parse().unwrap(), c[4].parse().unwrap());
            }
            "validity" if c[2] == "moisture_content" => {
                assert_eq!((c[3], c[5]), (c[4], "%"), "one measured moisture content");
                cells.entry(species()).or_default().4 = c[3].parse().unwrap();
            }
            other => panic!("unexpected row {other}"),
        }
    }
    // Table III (tan delta, coefficient of variation) and thesis Tableau 10 (EMC), normal heartwood.
    let printed = [("caribaea", 0.0132, 0.03, 10.1), ("nigra-laricio", 0.0088, 0.06, 9.6), ("patula", 0.0066, 0.05, 10.6),
        ("pinaster", 0.0120, 0.10, 10.5), ("sylvestris", 0.0106, 0.13, 10.5)];
    assert_eq!(cells.len(), printed.len());
    for (species, tan_delta, cv, emc) in printed {
        let (value, spread, lo, hi, mc) = cells[species];
        assert_eq!((value, lo, hi, mc), (tan_delta, 200.0, 600.0, emc), "{species}");
        assert!((spread - cv * tan_delta).abs() < 1e-12, "{species}: the spread is cv x mean");
        assert!((0.004..0.02).contains(&value), "{species}: wood's small-damping regime");
    }
}
