//! E3.4-iii battery (bead wf-root-guzez.4.8.3): bilinear terrain queries
//! exact at nodes/midpoints, the FlatnessCertificate issued from the REAL
//! committed E1.3 grid with numbers REPRODUCING the E1.3 artifact (a
//! cross-language cross-check: Rust re-fit vs the Python-recorded values),
//! refusal on the hill (uncertifiable, band named), PrelaunchPhase
//! admission + digest sensitivity, caps.
//! Repro: cargo test -p fs-flyer --test prelaunch_battery

use fs_flyer::prelaunch::{
    CERT_MAX_SLOPE, PrelaunchPhase, TerrainGrid, issue_flatness_certificate,
};
use fs_flyer::rail::RailSpec;

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e34iii\",\"case\":\"{case}\",{payload}}}");
}

/// Minimal extractor for the committed grid JSON (regular machine-written
/// format; no serde in the workspace graph).
fn load_kdh_rows() -> Vec<Vec<f64>> {
    let text = include_str!("../../../data/wright-flyer/terrain/kill-devil-hills-17x17-v1.json");
    let after = text.split("\"rows_south_to_north\"").nth(1).unwrap();
    let block = after.split('[').skip(1).collect::<Vec<_>>().join("[");
    let mut rows = Vec::new();
    for row_text in block.split(']') {
        let vals: Vec<f64> = row_text
            .split(',')
            .filter_map(|t| t.trim().trim_start_matches('[').trim().parse::<f64>().ok())
            .collect();
        if vals.len() == 17 {
            rows.push(vals);
        }
        if rows.len() == 17 {
            break;
        }
    }
    assert_eq!(rows.len(), 17, "grid extraction must recover all 17 rows");
    rows
}

#[test]
fn bilinear_is_exact_at_nodes_and_midpoints() {
    let g = TerrainGrid::new(
        10.0,
        vec![
            vec![0.0, 2.0, 4.0],
            vec![1.0, 3.0, 5.0],
            vec![2.0, 4.0, 6.0],
        ],
    )
    .unwrap();
    assert_eq!(g.height_m(10.0, 0.0).unwrap(), 2.0, "node exact");
    assert_eq!(g.height_m(20.0, 20.0).unwrap(), 6.0, "corner exact");
    // Midpoint of the first cell: mean of its 4 nodes = (0+2+1+3)/4 = 1.5.
    assert!((g.height_m(5.0, 5.0).unwrap() - 1.5).abs() < 1e-12);
    // A planar grid reproduces the plane ANYWHERE (bilinear is exact on planes).
    assert!((g.height_m(7.5, 12.5).unwrap() - (7.5 * 0.2 + 12.5 * 0.1)).abs() < 1e-12);
    // Outside refuses.
    assert_eq!(
        g.height_m(-0.1, 0.0).unwrap_err().code,
        "terrain-query-outside-domain"
    );
    assert_eq!(
        g.height_m(20.1, 0.0).unwrap_err().code,
        "terrain-query-outside-domain"
    );
    jlog("bilinear", "\"nodes\":\"exact\",\"planes\":\"exact\"");
}

#[test]
fn certificate_reproduces_the_e13_numbers_and_refuses_the_hill() {
    let g = TerrainGrid::new(125.0, load_kdh_rows()).unwrap();
    // The E1.3 launch flat: rows 10-16, cols 3-13. The Rust re-fit must
    // reproduce the Python-recorded artifact numbers (cross-language
    // cross-check of the SAME committed data).
    let cert = issue_flatness_certificate(&g, 10, 16, 3, 13).expect("the flat must certify");
    let slope_mag = (cert.slope.0.powi(2) + cert.slope.1.powi(2)).sqrt();
    assert!(
        (slope_mag - 0.000606).abs() < 2e-5,
        "slope {slope_mag} vs E1.3 0.000606"
    );
    assert!(
        (cert.rms_residual_m - 0.801).abs() < 0.01,
        "rms {} vs E1.3 0.801",
        cert.rms_residual_m
    );
    assert!((cert.max_abs_residual_m - 3.549).abs() < 0.01);
    assert!(slope_mag < CERT_MAX_SLOPE);
    // The hill region (Big Kill Devil Hill in-frame) must REFUSE with the
    // band named — the image model may not pretend the dune is flat.
    let refusal = issue_flatness_certificate(&g, 2, 8, 5, 11).unwrap_err();
    assert_eq!(refusal.code, "flatness-uncertifiable");
    // Degenerate region refuses.
    assert_eq!(
        issue_flatness_certificate(&g, 5, 5, 0, 3).unwrap_err().code,
        "flatness-region-invalid"
    );
    jlog(
        "certificate",
        &format!(
            "\"slope\":{slope_mag},\"rms\":{},\"hill\":\"refused\"",
            cert.rms_residual_m
        ),
    );
}

#[test]
fn prelaunch_schema_admits_and_digests_sensitively() {
    let g = TerrainGrid::new(125.0, load_kdh_rows()).unwrap();
    let flatness = issue_flatness_certificate(&g, 10, 16, 3, 13).unwrap();
    let phase = PrelaunchPhase {
        rail: RailSpec {
            z_rail_m: -0.3,
            length_m: 18.29,
            hysteresis_ticks: 3,
        },
        headwind_mps: 10.73,
        controls_rad: [0.0, 0.0],
        flatness,
    };
    phase.admit().unwrap();
    let base = phase.digest();
    assert_eq!(base, phase.digest(), "digest deterministic");
    let mut windier = phase.clone();
    windier.headwind_mps = 12.07;
    assert_ne!(base, windier.digest(), "identity must see the headwind");
    // Refusals: negative headwind; rail spec propagates.
    let mut bad = phase.clone();
    bad.headwind_mps = -1.0;
    assert_eq!(bad.admit().unwrap_err().code, "prelaunch-invalid");
    let mut bad_rail = phase;
    bad_rail.rail.hysteresis_ticks = 0;
    assert_eq!(bad_rail.admit().unwrap_err().code, "rail-spec-invalid");
    jlog("schema", &format!("\"digest\":\"{}\"", &base[..16]));
}

#[test]
fn grid_caps_and_shape_gates() {
    assert_eq!(
        TerrainGrid::new(125.0, vec![vec![1.0, 2.0]])
            .unwrap_err()
            .code,
        "terrain-grid-invalid"
    );
    assert_eq!(
        TerrainGrid::new(125.0, vec![vec![1.0, 2.0], vec![3.0]])
            .unwrap_err()
            .code,
        "terrain-grid-invalid"
    );
    assert_eq!(
        TerrainGrid::new(0.0, vec![vec![1.0, 2.0], vec![3.0, 4.0]])
            .unwrap_err()
            .code,
        "terrain-grid-invalid"
    );
    assert_eq!(
        TerrainGrid::new(125.0, vec![vec![1.0, f64::NAN], vec![3.0, 4.0]])
            .unwrap_err()
            .code,
        "terrain-grid-invalid"
    );
    jlog("gates", "\"shape/spacing/finite\":\"refused\"");
}
