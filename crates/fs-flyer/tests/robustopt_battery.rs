//! Robust multiobjective optimization battery (bead `frankensim-wf-root-guzez.9.5`, E8.4).

use fs_evidence::ColorRank;
use fs_flyer::robustopt::{
    run_robust_optimization, FlightOptimizationCandidate, RobustOptConfig, StructuralModelMode,
};

#[test]
fn robust_optimization_runs_and_reorders_under_crn_uncertainty() {
    let config = RobustOptConfig {
        structural_mode: StructuralModelMode::ActiveElastic,
        cvar_alpha: 0.90,
        ..Default::default()
    };

    // Design A: Aggressive gain with high nominal distance but poor tail under gust
    let cand_a = FlightOptimizationCandidate {
        name: "AggressiveHighNominal".into(),
        canard_trim_deg: 2.0,
        pilot_kp: 1.2,
        pilot_kq: 0.1,
        nominal_distance_m: 260.0,
        // Tail samples under severe wind perturbations drop significantly
        crn_distance_samples_m: vec![260.0, 255.0, 250.0, 240.0, 110.0, 95.0, 80.0],
    };

    // Design B: Flat robust gain with slightly lower nominal distance but consistent tail
    let cand_b = FlightOptimizationCandidate {
        name: "RobustFlatDamping".into(),
        canard_trim_deg: 1.8,
        pilot_kp: 0.8,
        pilot_kq: 0.4,
        nominal_distance_m: 245.0,
        crn_distance_samples_m: vec![245.0, 242.0, 240.0, 238.0, 230.0, 225.0, 220.0],
    };

    let receipt = run_robust_optimization(&config, &[cand_a, cand_b]).expect("optimization runs");

    assert_eq!(receipt.nominal_winner, "AggressiveHighNominal");
    assert_eq!(receipt.robust_winner, "RobustFlatDamping");
    assert!(receipt.robustness_reorders);
    assert_eq!(receipt.headline_rank, ColorRank::Estimated);
    assert_eq!(receipt.candidates_evaluated, 2);
    assert!(!receipt.receipt_digest.is_empty());
}

#[test]
fn hostile_twin_prescribed_kinematic_refuses() {
    let config = RobustOptConfig {
        structural_mode: StructuralModelMode::PrescribedKinematicEstimated,
        ..Default::default()
    };

    let cand = FlightOptimizationCandidate {
        name: "Candidate".into(),
        canard_trim_deg: 1.5,
        pilot_kp: 0.5,
        pilot_kq: 0.2,
        nominal_distance_m: 100.0,
        crn_distance_samples_m: vec![100.0, 95.0],
    };

    let err = run_robust_optimization(&config, &[cand]).expect_err("must refuse under prescribed kinematic");
    assert_eq!(err.code, "robust-opt-disabled-under-prescribed-kinematic");
}

#[test]
fn hostile_twin_rigid_mode_refuses() {
    let config = RobustOptConfig {
        structural_mode: StructuralModelMode::RigidEstimated,
        ..Default::default()
    };

    let cand = FlightOptimizationCandidate {
        name: "Candidate".into(),
        canard_trim_deg: 1.5,
        pilot_kp: 0.5,
        pilot_kq: 0.2,
        nominal_distance_m: 100.0,
        crn_distance_samples_m: vec![100.0, 95.0],
    };

    let err = run_robust_optimization(&config, &[cand]).expect_err("must refuse under rigid mode");
    assert_eq!(err.code, "robust-opt-requires-active-structural-model");
}

#[test]
fn hostile_twin_applicability_violation_refuses() {
    let config = RobustOptConfig {
        structural_mode: StructuralModelMode::ActiveElastic,
        max_kp: 1.5,
        ..Default::default()
    };

    let cand = FlightOptimizationCandidate {
        name: "ExcessiveKp".into(),
        canard_trim_deg: 1.5,
        pilot_kp: 2.5, // Beyond max_kp 1.5
        pilot_kq: 0.2,
        nominal_distance_m: 100.0,
        crn_distance_samples_m: vec![100.0, 95.0],
    };

    let err = run_robust_optimization(&config, &[cand]).expect_err("must refuse applicability violation");
    assert_eq!(err.code, "robust-opt-applicability-exceeded");
}
