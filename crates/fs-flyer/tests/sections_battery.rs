//! E4.1 battery (bead wf-root-guzez.5.2): V-01 TREND holdouts against the
//! 1901 anchors (the claims the frozen dossier PERMITS — sign, ordering,
//! order-of-magnitude; the anchors are holdout: the builder never sees
//! them), lineage/independence recorded per dataset, the synthesized-
//! stall role constraint enforced, convention-gate falsifier, golden.
//! Repro: cargo test -p fs-flyer --test sections_battery

use fs_airfoil::table::SurfaceKind;
use fs_flyer::sections::{build_v1_datasets, cl_3d};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v01\",\"case\":\"{case}\",{payload}}}");
}

const RE: f64 = 6.0;
const DEG: f64 = std::f64::consts::PI / 180.0;

#[test]
fn four_datasets_build_with_lineage_and_roles() {
    let ds = build_v1_datasets().unwrap();
    assert_eq!(ds.len(), 4);
    let kinds: Vec<_> = ds.iter().map(|d| d.kind).collect();
    for want in [
        SurfaceKind::Wing,
        SurfaceKind::Canard,
        SurfaceKind::Rudder,
        SurfaceKind::Prop,
    ] {
        assert!(kinds.contains(&want), "{want:?} dataset missing");
    }
    for d in &ds {
        assert!(!d.lineage.dossier_record.is_empty() && !d.lineage.independence_group.is_empty());
        assert_eq!(
            d.lineage.ceiling, "Estimated",
            "v1 ceiling is Estimated everywhere"
        );
        d.residual.validate(1e-9).unwrap();
        jlog(
            "lineage",
            &format!(
                "\"kind\":\"{:?}\",\"record\":\"{}\",\"group\":\"{}\"",
                d.kind, d.lineage.dossier_record, d.lineage.independence_group
            ),
        );
    }
    // The synthesized-stall constraint: non-prop datasets carry the
    // prior/trend-only role text (never Wright-specific deep stall).
    let wing = ds.iter().find(|d| d.kind == SurfaceKind::Wing).unwrap();
    assert!(wing.lineage.role.contains("prior/trend role ONLY"));
}

#[test]
fn trend_holdout_anchor_12_at_5_degrees() {
    // 1901 anchor (HOLDOUT — nothing fitted): model #12, 1/20 camber, AR 6,
    // CL_wright 0.515 -> modern-equivalent ~0.659. Permitted claim: TREND —
    // correct sign and within a factor of 2 (the 1901 tunnel Re ~5e4 reads
    // low; quantitative match is FORBIDDEN by the dossier boundary).
    let ds = build_v1_datasets().unwrap();
    let wing = ds.iter().find(|d| d.kind == SurfaceKind::Wing).unwrap();
    let cl = cl_3d(wing, 5.0 * DEG, RE, 6.0).unwrap();
    let anchor = 0.659;
    assert!(cl > 0.0, "sign");
    assert!(
        cl / anchor > 0.5 && cl / anchor < 2.0,
        "order: {cl} vs anchor {anchor}"
    );
    jlog(
        "anchor-12",
        &format!("\"cl\":{cl},\"anchor_modern\":{anchor}"),
    );
}

#[test]
fn trend_holdout_anchor_7_high_alpha() {
    // Wilbur verbatim: 119% at 17.5 deg (modern-equiv ~1.52 with the 1.28
    // convention factor). Trend clause: our blended curve at 17.5 deg is
    // large (>0.8) and within a factor 2 of the anchor.
    let ds = build_v1_datasets().unwrap();
    let wing = ds.iter().find(|d| d.kind == SurfaceKind::Wing).unwrap();
    let cl = cl_3d(wing, 17.5 * DEG, RE, 6.0).unwrap();
    let anchor = 1.523;
    assert!(cl > 0.8, "high-alpha lift must remain large: {cl}");
    assert!(
        cl / anchor > 0.5 && cl / anchor < 2.0,
        "order: {cl} vs {anchor}"
    );
    jlog(
        "anchor-7",
        &format!("\"cl\":{cl},\"anchor_modern\":{anchor}"),
    );
}

#[test]
fn trend_camber_ordering_and_zero_lift_signs() {
    // The Wrights' central finding: curvature raises lift. A 1/12-camber
    // twin of the wing dataset must out-lift 1/20 at 5 deg; the symmetric
    // rudder must have zero lift at zero alpha while cambered sections
    // lift positively there.
    let ds = build_v1_datasets().unwrap();
    let wing = ds.iter().find(|d| d.kind == SurfaceKind::Wing).unwrap();
    let rudder = ds.iter().find(|d| d.kind == SurfaceKind::Rudder).unwrap();
    let mut deeper = wing.clone();
    deeper.camber_ratio = 1.0 / 12.0;
    let cl20 = cl_3d(wing, 5.0 * DEG, RE, 6.0).unwrap();
    let cl12 = cl_3d(&deeper, 5.0 * DEG, RE, 6.0).unwrap();
    assert!(cl12 > cl20, "curvature must raise lift ({cl12} vs {cl20})");
    assert_eq!(
        cl_3d(rudder, 0.0, RE, 3.0).unwrap(),
        0.0,
        "symmetric zero-lift"
    );
    assert!(
        cl_3d(wing, 0.0, RE, 6.0).unwrap() > 0.0,
        "cambered lifts at zero alpha"
    );
    // Blend continuity: no jump crossing the attached limit.
    let a = cl_3d(wing, 0.299, RE, 6.0).unwrap();
    let b = cl_3d(wing, 0.301, RE, 6.0).unwrap();
    assert!(
        (a - b).abs() < 0.05,
        "blend must be continuous ({a} vs {b})"
    );
    jlog("trends", &format!("\"cl20\":{cl20},\"cl12\":{cl12}"));
}

#[test]
fn convention_gate_falsifier() {
    // Corrupting a dataset's convention id must refuse at validation —
    // the reexpression-v1 gate runs on REAL datasets, not just fixtures.
    let ds = build_v1_datasets().unwrap();
    let mut broken = ds[0].residual.clone();
    broken.conventions.axes_id = "z-up-graphics".into();
    assert_eq!(
        broken.validate(1e-9).unwrap_err().code,
        "convention-block-mismatch"
    );
}

#[test]
fn sections_golden_digest() {
    // A 13-point wing polar sweep: exact bits (measure-then-pin).
    let ds = build_v1_datasets().unwrap();
    let wing = ds.iter().find(|d| d.kind == SurfaceKind::Wing).unwrap();
    let mut payload = Vec::new();
    for i in 0..13 {
        let alpha = (-10.0 + 5.0 * f64::from(i)) * DEG;
        let cl = cl_3d(wing, alpha, RE, 6.0).unwrap();
        payload.extend_from_slice(&cl.to_bits().to_le_bytes());
    }
    let digest = fs_blake3::hash_domain("org.frankensim.fs-flyer.v01-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "78693a7b06ca3b119ccb4f156b64f1881484d26b3197c4c85b4af452f4cb4009",
        "sections golden moved — determinism regression or an intentional \
         dataset change requiring the golden-bump protocol"
    );
}
