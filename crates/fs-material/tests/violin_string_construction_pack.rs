//! Physicality gates for the violin-string construction pack (bead `frankensim-irp-v5-instrument-realism-pexfp.192`, plan item M8-W2b),
//! from Table 1 of Firth 1987, J. Catgut Acoust. Soc. 47:17 (numerical facts with attribution; the table is not copied): 16 Pirastro
//! violin strings with core, core-with-underwrap and finished diameters, gut twist pitch, floss counts and the metal wrap's wire diameter.
//!
//! The gates are the geometry of a wrapped string: the finished diameter is not below the core with its underwrap, which is not below
//! the core, and a single layer of the metal wire on the underwrap gives a diameter that the polished finished string undercuts by at
//! most 12 percent and never exceeds (OD <= underwrap + 2 wire diameters, OD >= 0.88 of it). Schema hygiene: every scalar has an
//! uncertainty and a validity row.

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

#[test]
fn violin_string_construction_pack_passes_geometry_gates() {
    let dir = "violin-strings-construction-firth-jcas1987";
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.contains("\nlicense\tPublished-Scientific-Facts-Attribution-Retained\n"));
    assert!(manifest.contains("Journal of the Catgut Acoustical Society 47:17-18"));
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    let mut s: BTreeMap<String, f64> = BTreeMap::new();
    let (mut ids, mut unc, mut val) = (Vec::new(), Vec::new(), Vec::new());
    let mut observations = 0;
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
            "observation" => observations += 1,
            other => panic!("unknown row kind {other}"),
        }
    }
    for id in &ids {
        assert!(
            unc.contains(id) && val.contains(id),
            "{id}: uncertainty and validity rows"
        );
    }
    assert_eq!(observations, 16, "one observation per string");

    let keys: Vec<String> = s
        .keys()
        .filter_map(|k| k.strip_prefix("outer_diameter_").map(str::to_string))
        .collect();
    assert_eq!(keys.len(), 16);
    let mut with_underwrap = 0;
    for k in &keys {
        let (od, core) = (
            s[&format!("outer_diameter_{k}")],
            s[&format!("core_diameter_{k}")],
        );
        let wire = s[&format!("wrap_wire_diameter_{k}")];
        assert!(od > core, "{k}: finished string thicker than its core");
        if let Some(cu) = s.get(&format!("core_with_underwrap_diameter_{k}")) {
            with_underwrap += 1;
            assert!(*cu >= core, "{k}: underwrap not thinner than the core");
            assert!(
                od >= *cu,
                "{k}: finished string not thinner than the underwrap"
            );
            let one_layer = cu + 2.0 * wire;
            assert!(
                od <= one_layer * 1.0001 && od >= 0.88 * one_layer,
                "{k}: OD {od} against one wire layer {one_layer}"
            );
        }
    }
    assert_eq!(with_underwrap, 14);
}
