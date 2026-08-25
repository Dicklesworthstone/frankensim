//! Tests for automated mesh convergence ladder evaluation and Richardson extrapolation.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.6`

use fs_evidence::Color;
use fs_ladder::{ConvergenceEvaluator, ConvergencePlan, ConvergenceStatus, MeshRung};

#[test]
fn test_convergence_insufficient_rungs_fails_closed() {
    // 1. One rung only
    let plan1 = ConvergencePlan::new("junction_maximum", 2.0)
        .with_rung(MeshRung::new(0, "mesh_coarse", 0.04, "m", 1000, "junction_maximum", 350.0, "K"));
    let res1 = ConvergenceEvaluator::evaluate(&plan1);
    assert_eq!(res1.status, ConvergenceStatus::InsufficientRungs);
    assert!(res1.observed_order.is_none());
    assert!(res1.richardson_extrapolated_qoi.is_none());

    // 2. Two rungs only (must NEVER mint a two-point order)
    let plan2 = plan1.with_rung(MeshRung::new(1, "mesh_medium", 0.02, "m", 8000, "junction_maximum", 342.5, "K"));
    let res2 = ConvergenceEvaluator::evaluate(&plan2);
    assert_eq!(res2.status, ConvergenceStatus::InsufficientRungs);
    assert!(res2.observed_order.is_none());
    assert!(res2.richardson_extrapolated_qoi.is_none());
    assert!(matches!(res2.evidence_color, Color::Estimated { .. }));
}

#[test]
fn test_convergence_monotone_asymptotic_order_2() {
    // Synthetic order 2 convergence towards true continuum value 340.0:
    // h3 = 0.04 -> T = 340.0 + 10.0 * (0.04/0.04)^2 = 350.0
    // h2 = 0.02 -> T = 340.0 + 10.0 * (0.02/0.04)^2 = 342.5
    // h1 = 0.01 -> T = 340.0 + 10.0 * (0.01/0.04)^2 = 340.625
    let plan = ConvergencePlan::new("junction_maximum", 2.0)
        .with_rung(MeshRung::new(0, "mesh_coarse", 0.04, "m", 1000, "junction_maximum", 350.0, "K"))
        .with_rung(MeshRung::new(1, "mesh_medium", 0.02, "m", 8000, "junction_maximum", 342.5, "K"))
        .with_rung(MeshRung::new(2, "mesh_fine", 0.01, "m", 64000, "junction_maximum", 340.625, "K"));

    let res = ConvergenceEvaluator::evaluate(&plan);
    assert_eq!(res.status, ConvergenceStatus::Asymptotic);

    let observed = res.observed_order.expect("observed order fitted");
    assert!((observed - 2.0).abs() < 1e-4, "observed order should be 2.0, got {observed}");

    let extrap = res.richardson_extrapolated_qoi.expect("extrapolated value");
    assert!((extrap - 340.0).abs() < 1e-4, "extrapolated value should be 340.0, got {extrap}");

    assert!(res.discretization_error_gci.is_some());
    assert!(matches!(res.evidence_color, Color::Verified { .. }));
}

#[test]
fn test_convergence_oscillatory_refuses_extrapolation() {
    // Non-monotone sequence: 350.0 -> 335.0 -> 345.0
    let plan = ConvergencePlan::new("junction_maximum", 2.0)
        .with_rung(MeshRung::new(0, "mesh_coarse", 0.04, "m", 1000, "junction_maximum", 350.0, "K"))
        .with_rung(MeshRung::new(1, "mesh_medium", 0.02, "m", 8000, "junction_maximum", 335.0, "K"))
        .with_rung(MeshRung::new(2, "mesh_fine", 0.01, "m", 64000, "junction_maximum", 345.0, "K"));

    let res = ConvergenceEvaluator::evaluate(&plan);
    assert_eq!(res.status, ConvergenceStatus::Oscillatory);
    assert!(res.richardson_extrapolated_qoi.is_none(), "extrapolation must be refused on oscillatory ladder");
    assert!(matches!(res.evidence_color, Color::Estimated { .. }));
}

#[test]
fn test_convergence_determinism() {
    let plan1 = ConvergencePlan::new("junction_maximum", 2.0)
        .with_rung(MeshRung::new(0, "mesh_coarse", 0.04, "m", 1000, "junction_maximum", 350.0, "K"))
        .with_rung(MeshRung::new(1, "mesh_medium", 0.02, "m", 8000, "junction_maximum", 342.5, "K"))
        .with_rung(MeshRung::new(2, "mesh_fine", 0.01, "m", 64000, "junction_maximum", 340.625, "K"));

    let plan2 = plan1.clone();

    let res1 = ConvergenceEvaluator::evaluate(&plan1);
    let res2 = ConvergenceEvaluator::evaluate(&plan2);

    assert_eq!(res1.content_hash(), res2.content_hash(), "content hash bit-identical");
}
