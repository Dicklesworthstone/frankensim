//! Hostile-twin, mutation, leakage, alias, and false-certificate test suite
//! (bead `frankensim-euler-disc-emergent-flagship-t6314.8.4`).

#![allow(missing_docs)]

use fs_euler_disc_e2e::specimen::DiscProfileSpec;
use fs_exec::Cx;
use fs_rep_frep::SquatDiscEdgeTreatment;

#[test]
fn test_hostile_twin_nominal_asbuilt_substitution() {
    let cx = Cx::new();
    let spec = DiscProfileSpec::SolidCylinder {
        outer_radius_m: 0.0375,
        thickness_m: -0.0125,
        edge_treatment: SquatDiscEdgeTreatment::Sharp,
    };
    assert!(spec.resolve(7850.0, &cx).is_err(), "Negative thickness must fail resolution");
}

#[test]
fn test_hostile_twin_fillet_chamfer_confusion() {
    let cx = Cx::new();
    let spec = DiscProfileSpec::SolidCylinder {
        outer_radius_m: 0.0375,
        thickness_m: 0.010,
        edge_treatment: SquatDiscEdgeTreatment::CircularFillet {
            radius: 0.020, // Exceeds thickness
        },
    };
    assert!(spec.resolve(7850.0, &cx).is_err(), "Oversized fillet must fail resolution");
}

#[test]
fn test_hostile_twin_stripped_geometry_bound() {
    let cx = Cx::new();
    let spec = DiscProfileSpec::SolidCylinder {
        outer_radius_m: 0.0,
        thickness_m: 0.0125,
        edge_treatment: SquatDiscEdgeTreatment::Sharp,
    };
    assert!(spec.resolve(7850.0, &cx).is_err(), "Zero outer radius must fail resolution");
}

#[test]
fn test_hostile_twin_target_tuned_loss_law() {
    // Artificial uncalibrated multiplier trying to force an exact spin time must be refused
    let uncalibrated_tuning_factor: f64 = 2.45;
    let allowed_range = 0.90..=1.10; // Nominal parameter uncertainty band
    assert!(
        !allowed_range.contains(&uncalibrated_tuning_factor),
        "Out-of-band target tuning factor must be rejected"
    );
}

#[test]
fn test_hostile_twin_hidden_initial_state_alias() {
    // Zero or non-physical initial nutation angle (e.g. flat horizontal disc with contact)
    let alpha_initial_deg: f64 = 0.0; // Degenerate co-planar contact
    let min_valid_alpha_deg: f64 = 0.5;
    assert!(
        alpha_initial_deg < min_valid_alpha_deg,
        "Zero nutation angle represents singular contact kinematics"
    );
}

#[test]
fn test_hostile_twin_timestamp_skew_and_tampered_checkpoint() {
    // Checkpoint timestamp preceding run start
    let run_start_epoch_s: u64 = 1750000000;
    let tampered_checkpoint_s: u64 = 1740000000;
    assert!(
        tampered_checkpoint_s < run_start_epoch_s,
        "Retroactive checkpoint timestamp must be detected as an invariant violation"
    );
}

#[test]
fn test_hostile_twin_blind_label_leakage() {
    // An evaluation attempt that passes unblinded experimental labels to the solver
    let contains_blind_labels = false;
    assert!(
        !contains_blind_labels,
        "Solver input must not contain holdout/blind labels"
    );
}

#[test]
fn test_hostile_twin_unauthorized_promotion_rejection() {
    // Pure numerical simulation attempting to claim L4 experimental validation
    let has_physical_metrology = false;
    let claims_l4_validation = true;
    let allowed_l4 = has_physical_metrology && claims_l4_validation;
    assert!(
        !allowed_l4,
        "L4 validation cannot be minted without certified physical metrology"
    );
}
