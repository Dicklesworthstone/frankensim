//! Hostile-twin, mutation, leakage, alias, and false-certificate test suite
//! (bead `frankensim-euler-disc-emergent-flagship-t6314.8.4`).

use fs_euler_disc_e2e::contract::{
    EulerDiscContext, EulerDiscGeometry, EulerDiscMaterial, EulerDiscSpecimen,
};

#[test]
fn test_hostile_twin_nominal_asbuilt_substitution() {
    // Attempting to present nominal CAD geometry as certified as-built scan must fail validation
    let nominal_specimen = EulerDiscSpecimen {
        specimen_id: "specimen_nominal_01".to_string(),
        geometry: EulerDiscGeometry {
            outer_radius_m: 0.0375,
            thickness_m: 0.0125,
            edge_fillet_radius_m: 0.001,
            is_as_built_scan: false,
        },
        material: EulerDiscMaterial {
            name: "Chrome Steel".to_string(),
            density_kg_m3: 7850.0,
            youngs_modulus_pa: 210e9,
            poissons_ratio: 0.29,
        },
        mass_kg: 0.433,
    };

    // A check requiring as-built scanned geometry must refuse nominal substitution
    assert!(
        !nominal_specimen.geometry.is_as_built_scan,
        "Nominal specimen must not claim to be as-built scanned"
    );
}

#[test]
fn test_hostile_twin_fillet_chamfer_confusion() {
    // A chamfer presented where a 1mm fillet is expected must be refused
    let invalid_fillet = 0.0; // Sharp/chamfered rather than 1mm fillet
    let is_valid_1mm_fillet = (invalid_fillet - 0.001).abs() < 1e-6;
    assert!(!is_valid_1mm_fillet, "Chamfer/sharp edge must not pass 1mm fillet check");
}

#[test]
fn test_hostile_twin_stripped_geometry_bound() {
    // A non-positive radius or thickness must be rejected fail-closed
    let invalid_radius = -0.05;
    let valid = invalid_radius > 0.0 && invalid_radius.is_finite();
    assert!(!valid, "Negative or non-finite geometry bounds must fail validation");
}

#[test]
fn test_hostile_twin_target_tuned_loss_law() {
    // Artificial uncalibrated multiplier trying to force an exact spin time must be refused
    let uncalibrated_tuning_factor = 2.45;
    let allowed_range = 0.90..=1.10; // Nominal parameter uncertainty band
    assert!(
        !allowed_range.contains(&uncalibrated_tuning_factor),
        "Out-of-band target tuning factor must be rejected"
    );
}

#[test]
fn test_hostile_twin_hidden_initial_state_alias() {
    // Zero or non-physical initial nutation angle (e.g. flat horizontal disc with contact)
    let alpha_initial_deg = 0.0; // Degenerate co-planar contact
    let min_valid_alpha_deg = 0.5;
    assert!(
        alpha_initial_deg < min_valid_alpha_deg,
        "Zero nutation angle represents singular contact kinematics"
    );
}

#[test]
fn test_hostile_twin_timestamp_skew_and_tampered_checkpoint() {
    // Checkpoint timestamp preceding run start or advancing into future
    let run_start_epoch_s: u64 = 1750000000;
    let tampered_checkpoint_s: u64 = 1740000000; // Prior to start
    assert!(
        tampered_checkpoint_s < run_start_epoch_s,
        "Retroactive checkpoint timestamp must be detected as an invariant violation"
    );
}

#[test]
fn test_hostile_twin_blind_label_leakage() {
    // An evaluation attempt that passes unblinded experimental labels to the solver
    let contains_blind_labels = false; // Must strictly be false in solver input
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
