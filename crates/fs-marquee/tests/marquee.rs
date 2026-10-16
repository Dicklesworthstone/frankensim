//! fs-marquee contract conformance.
//!
//! The default crate is intentionally an L6 admission/status shell. The
//! `marquee` feature exposes a smoke-tier in-process runner, but it must not
//! claim the full nightly golden lane or ledger/filesystem side effects.

#[cfg(feature = "marquee")]
use fs_marquee::study::{PlateWithHoles, StudyConfig};
use fs_marquee::{MarqueeStatus, VERSION, scope_summary, status};

fn expected_status() -> MarqueeStatus {
    if cfg!(feature = "marquee") {
        MarqueeStatus::SmokeRunnerAvailable
    } else {
        MarqueeStatus::Disabled
    }
}

#[test]
fn marquee_status_matches_feature_gate() {
    assert_eq!(status(), expected_status());
    assert!(!VERSION.is_empty());
}

#[test]
fn marquee_scope_keeps_nightly_golden_boundary_explicit() {
    let summary = scope_summary();
    assert!(summary.contains("raw SDF"));
    assert!(summary.contains("CutFEM"));
    // The smoke runner shipped; the nightly golden lane is the
    // remaining no-claim boundary.
    assert!(summary.contains("nightly golden pending"));
}

#[cfg(feature = "marquee")]
#[test]
fn default_feature_set_exposes_smoke_runner() {
    assert_eq!(status(), MarqueeStatus::SmokeRunnerAvailable);

    // These values keep the default-enabled public runner API covered without
    // making an explicit `--no-default-features` status-only build fail to
    // compile.
    let _design = PlateWithHoles {
        centers: vec![[0.5, 0.5]],
        radii: vec![0.1],
    };
    let _config = StudyConfig {
        level: 1,
        steps: 1,
        step_size: 0.0,
        area_target: 0.1,
        r_min: 0.05,
        r_max: 0.2,
    };
}

#[cfg(feature = "marquee")]
#[test]
fn marquee_runner_rejects_invalid_inputs_before_solver() {
    use fs_marquee::study::run_study;

    let design = PlateWithHoles {
        centers: Vec::new(),
        radii: Vec::new(),
    };
    let config = StudyConfig {
        level: 1,
        steps: 1,
        step_size: 1.0,
        area_target: 0.9,
        r_min: 0.05,
        r_max: 0.2,
    };

    assert!(run_study(design, &config).is_err());
}

#[cfg(feature = "marquee")]
#[test]
fn marquee_empty_design_sdf_is_total_even_when_runner_rejects_it() {
    use fs_cutfem::sdf::CutSdf;

    let design = PlateWithHoles {
        centers: Vec::new(),
        radii: Vec::new(),
    };

    let value = design.value([0.5, 0.5]);
    let gradient = design.gradient([0.5, 0.5]);

    assert!(value.is_infinite() && value.is_sign_negative());
    assert!(
        gradient
            .iter()
            .all(|component| component.abs() <= f64::EPSILON)
    );
}
