//! Physicality gates for the whole-string matdb pack of 35 violin strings (bead `frankensim-irp-v5-instrument-realism-pexfp.192`, plan item
//! M8-W2b), from the table of Pickering 1985, J. Catgut Acoust. Soc. 44:6-8 (numerical facts with attribution; the table is not copied).
//!
//! The table's own columns give two identities, which are the gates that the transcription is right: the printed wave impedance is
//! sqrt(tension x linear mass) (to 0.3 %), and the printed tension at the open-string frequency f over the test bed's 329 mm is
//! linear mass x (2 L f)^2 of the relaxed string (to 3 %: the relaxed linear mass is not the stretched one). Also: maximum diameter not
//! below minimum, impedance and linear mass rising from E to G, every scalar with an uncertainty and a validity row.

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

/// Open-string frequency in Hz of each string (the tension in the table is at these frequencies).
fn open_hz(pitch: char) -> f64 {
    match pitch {
        'G' => 196.0,
        'D' => 293.7,
        'A' => 440.0,
        'E' => 659.3,
        _ => unreachable!(),
    }
}

#[test]
fn violin_strings_pack_passes_physicality_gates() {
    let dir = "violin-strings-whole-pickering-jcas1985";
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.contains("\nlicense\tPublished-Scientific-Facts-Attribution-Retained\n"));
    assert!(manifest.contains("Journal of the Catgut Acoustical Society 44:6-8"));
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    let mut s: BTreeMap<String, f64> = BTreeMap::new();
    let (mut ids, mut unc, mut val) = (Vec::new(), Vec::new(), Vec::new());
    let mut observations = Vec::new();
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        match c[0] {
            "scalar" => {
                s.insert(c[3].to_string(), c[4].parse().expect("value"));
                ids.push(c[1].to_string());
            }
            "uncertainty" => unc.push(c[1].to_string()),
            "validity" => val.push(c[1].to_string()),
            "observation" => observations.push(c[2].to_string()),
            other => panic!("unknown row kind {other}"),
        }
    }
    for id in &ids {
        assert!(
            unc.contains(id) && val.contains(id),
            "{id}: uncertainty and validity rows"
        );
    }
    assert_eq!(observations.len(), 35, "one observation per string");

    let mut by_pitch: BTreeMap<char, Vec<(f64, f64)>> = BTreeMap::new();
    for (i, obs) in observations.iter().enumerate() {
        let k = format!("s{:02}", i + 1);
        let pitch = obs
            .split(' ')
            .find(|w| ["G", "D", "A", "E"].contains(w))
            .expect("pitch letter")
            .chars()
            .next()
            .unwrap();
        let (t, z, mu) = (
            s[&format!("tension_{k}")],
            s[&format!("wave_impedance_{k}")],
            s[&format!("linear_density_{k}")],
        );
        assert!(
            s[&format!("diameter_max_{k}")] >= s[&format!("diameter_min_{k}")],
            "{k}: diameter"
        );
        let z_calc = (t * mu).sqrt();
        assert!(
            (z_calc - z).abs() / z < 0.003,
            "{k}: sqrt(T mu) = {z_calc} against the printed {z}"
        );
        let t_calc = mu * (2.0 * 0.329 * open_hz(pitch)).powi(2);
        assert!(
            (t_calc - t).abs() / t < 0.03,
            "{k}: mu (2 L f)^2 = {t_calc} against the printed {t}"
        );
        by_pitch.entry(pitch).or_default().push((z, mu));
    }
    // Impedance and linear mass rise from the E strings to the G strings (the means of each pitch).
    let mean = |p: char, f: fn(&(f64, f64)) -> f64| {
        let v = &by_pitch[&p];
        v.iter().map(f).sum::<f64>() / v.len() as f64
    };
    for f in [(|x: &(f64, f64)| x.0) as fn(&(f64, f64)) -> f64, |x| x.1] {
        assert!(
            mean('E', f) < mean('A', f)
                && mean('A', f) < mean('D', f)
                && mean('D', f) < mean('G', f)
        );
    }
    assert_eq!(by_pitch.values().map(Vec::len).sum::<usize>(), 35);
}
