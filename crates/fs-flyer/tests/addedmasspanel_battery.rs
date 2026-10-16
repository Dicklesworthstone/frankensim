//! E4.2c battery (bead wf-root-guzez.5.5): the extracted cross-term
//! artifact. Identity-hashed (id stable twice, moved by any grid
//! edit); the AnalyticStrip discrepancy RECORDED and non-vacuous;
//! node exactness + the analytic derivatives verified against an
//! OFFLINE finite-difference cross-check (FD lives only in this
//! test); the runtime-FD entry refuses; out-of-domain refuses AT the
//! bounds admitting; axis caps at cap AND cap+1; artifact-id golden.
//! Repro: cargo test -p fs-flyer --test addedmasspanel_battery

use fs_flyer::addedmasspanel::{
    EXTRACTION_VERSION, ExtractionGrid, MAX_AXIS, eval_with_derivatives, extract,
    runtime_finite_difference,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-addedmasspanel\",\"case\":\"{case}\",{payload}}}");
}

fn pinned_grid() -> ExtractionGrid {
    ExtractionGrid {
        deflection_rad: vec![-0.15, 0.0, 0.15],
        height_m: vec![0.5, 1.5, 3.0],
        warp_rad: vec![-0.08, 0.0, 0.08],
    }
}

#[test]
fn extraction_is_identity_hashed_and_discrepancy_recorded() {
    let art = extract(pinned_grid()).unwrap();
    assert_eq!(art.extraction_tier, EXTRACTION_VERSION, "tier declared");
    // Identity: stable twice; ANY grid edit moves it.
    let again = extract(pinned_grid()).unwrap();
    assert_eq!(art.artifact_id, again.artifact_id, "bit-identical twice");
    let mut edited = pinned_grid();
    edited.height_m[2] = 3.1;
    assert_ne!(
        art.artifact_id,
        extract(edited).unwrap().artifact_id,
        "grid edits move the identity"
    );
    // Discrepancy vs the NOMINAL AnalyticStrip column: RECORDED and
    // non-vacuous (the deformed corners genuinely differ), per
    // component — never a totals-only norm.
    assert!(
        art.baseline_discrepancy.iter().any(|d| *d > 1e-6),
        "the deformation must move some cross component: {:?}",
        art.baseline_discrepancy
    );
    for d in &art.baseline_discrepancy {
        assert!(d.is_finite() && *d >= 0.0);
    }
    jlog(
        "identity",
        &format!(
            "\"artifact_id\":\"{}\",\"discrepancy\":{:?}",
            art.artifact_id, art.baseline_discrepancy
        ),
    );
}

#[test]
fn interpolation_is_node_exact_and_derivatives_are_analytic() {
    let art = extract(pinned_grid()).unwrap();
    let g = pinned_grid();
    // Node exactness at every grid node (per-node oracle).
    for (di, d) in g.deflection_rad.iter().enumerate() {
        for (hi, h) in g.height_m.iter().enumerate() {
            for (wi, w) in g.warp_rad.iter().enumerate() {
                let (v, _) = eval_with_derivatives(&art, *d, *h, *w).unwrap();
                for k in 0..6 {
                    let stored = art.values[di][hi][wi][k];
                    assert!(
                        (v[k] - stored).abs() <= 1e-12 * stored.abs().max(1e-12),
                        "node ({di},{hi},{wi})[{k}]: {} vs {stored}",
                        v[k]
                    );
                }
            }
        }
    }
    // Analytic derivatives vs an OFFLINE central-difference check at
    // an interior point (FD exists ONLY here, in the verifier).
    let p = (0.05, 1.0, 0.03);
    let (_, grad) = eval_with_derivatives(&art, p.0, p.1, p.2).unwrap();
    let h = 1e-6;
    for k in 0..6 {
        let fd = [
            (eval_with_derivatives(&art, p.0 + h, p.1, p.2).unwrap().0[k]
                - eval_with_derivatives(&art, p.0 - h, p.1, p.2).unwrap().0[k])
                / (2.0 * h),
            (eval_with_derivatives(&art, p.0, p.1 + h, p.2).unwrap().0[k]
                - eval_with_derivatives(&art, p.0, p.1 - h, p.2).unwrap().0[k])
                / (2.0 * h),
            (eval_with_derivatives(&art, p.0, p.1, p.2 + h).unwrap().0[k]
                - eval_with_derivatives(&art, p.0, p.1, p.2 - h).unwrap().0[k])
                / (2.0 * h),
        ];
        for a in 0..3 {
            let scale = grad[k][a].abs().max(1e-6);
            assert!(
                (grad[k][a] - fd[a]).abs() / scale < 1e-6,
                "component {k} axis {a}: analytic {} vs fd {}",
                grad[k][a],
                fd[a]
            );
        }
    }
    jlog("derivatives", "\"analytic_vs_offline_fd\":\"verified\"");
}

#[test]
fn runtime_fd_is_forbidden_and_domain_is_bounded() {
    let art = extract(pinned_grid()).unwrap();
    // The DONE-WHEN refusal: runtime finite-differencing.
    let err = runtime_finite_difference(&art, 0.0, 1.0, 0.0).unwrap_err();
    assert_eq!(err.code, "runtime-fd-forbidden");
    // Domain: AT the bounds admits; beyond refuses (no extrapolation).
    assert!(eval_with_derivatives(&art, -0.15, 0.5, -0.08).is_ok());
    assert!(eval_with_derivatives(&art, 0.15, 3.0, 0.08).is_ok());
    for (d, h, w) in [
        (0.150_000_1, 1.0, 0.0),
        (0.0, 0.499_999, 0.0),
        (0.0, 1.0, -0.080_001),
    ] {
        assert_eq!(
            eval_with_derivatives(&art, d, h, w).unwrap_err().code,
            "crossterm-out-of-domain",
            "({d},{h},{w})"
        );
    }
    jlog("refusals", &format!("\"runtime_fd\":\"{}\"", err.code));
}

#[test]
fn grid_caps_and_artifact_golden() {
    // Axis caps: 16 admits, 17 refuses; short/unordered refuse.
    let mk = |n: usize| ExtractionGrid {
        deflection_rad: (0..n).map(|i| i as f64 * 0.01).collect(),
        height_m: vec![0.5, 1.5],
        warp_rad: vec![-0.05, 0.05],
    };
    assert!(extract(mk(MAX_AXIS)).is_ok(), "AT cap");
    assert_eq!(
        extract(mk(MAX_AXIS + 1)).unwrap_err().code,
        "extraction-grid-invalid"
    );
    assert_eq!(extract(mk(1)).unwrap_err().code, "extraction-grid-invalid");
    let unordered = ExtractionGrid {
        deflection_rad: vec![0.1, 0.0],
        height_m: vec![0.5, 1.5],
        warp_rad: vec![-0.05, 0.05],
    };
    assert_eq!(
        extract(unordered).unwrap_err().code,
        "extraction-grid-invalid"
    );
    // Golden: the pinned grid's artifact id (measure-then-pin).
    let art = extract(pinned_grid()).unwrap();
    jlog(
        "golden",
        &format!("\"artifact_id\":\"{}\"", art.artifact_id),
    );
    assert_eq!(
        art.artifact_id, "40695911b35f787cb75715c8552949b451d5f2599fffc4e66504d24339ea2c81",
        "cross-term artifact moved — determinism regression or an \
         intentional extraction change requiring the golden-bump protocol"
    );
}
