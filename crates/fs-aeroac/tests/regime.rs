//! l011o: 2D tonal-lock refusal, broadband admission gate, 3D operator smoke.

use fs_aeroac::regime::{
    PINNED_2D_CENTRAL_MOMENT_TONAL, SlotJet3dFollowUp, TONAL_FLATNESS_CEILING,
    TWO_D_INVERSE_CASCADE, admit_broadband_spectrum, evaluate_slot_jet_3d_operator,
    two_d_broadband_refusal,
};
use fs_aeroac::{AeroacError, SCOPE_STATEMENT};

#[test]
fn pinned_2d_central_moment_rows_are_tonal_and_in_regime() {
    assert!(PINNED_2D_CENTRAL_MOMENT_TONAL.len() >= 8);
    for row in PINNED_2D_CENTRAL_MOMENT_TONAL {
        assert!(row.ran_in_regime);
        assert!(row.reynolds > 0.0 && row.slot_half > 0.0 && row.strouhal > 0.0);
        assert!(row.flatness_upper < TONAL_FLATNESS_CEILING);
        assert!(admit_broadband_spectrum(row.flatness_upper).is_err());
    }
}

#[test]
fn two_d_refusal_names_the_inverse_cascade_and_scope_law() {
    let refusal = two_d_broadband_refusal();
    assert_eq!(refusal.mechanism, TWO_D_INVERSE_CASCADE);
    assert_eq!(refusal.scope, SCOPE_STATEMENT);
    assert!(refusal.mechanism.contains("inverse-cascade"));
    assert!(refusal.scope.contains("NOT absolute SPL"));
    assert_eq!(refusal.rows.len(), PINNED_2D_CENTRAL_MOMENT_TONAL.len());
}

#[test]
fn tonal_limit_cycle_cannot_be_cataloged_as_broadband() {
    let err = admit_broadband_spectrum(1.0e-18).expect_err("tone");
    assert!(matches!(err, AeroacError::InvalidParameter { .. }));
    admit_broadband_spectrum(0.2).expect("genuinely broadband flatness admits");
}

#[test]
fn slot_jet_3d_follow_up_operator_is_live_but_does_not_claim_broadband() {
    let spec = SlotJet3dFollowUp::minimal_central_moment();
    spec.validate().expect("spec");
    let smoke = evaluate_slot_jet_3d_operator(spec).expect("smoke");
    assert!(smoke.collision_admits);
    assert!(smoke.equilibrium_live);
    assert!(
        !smoke.broadband_demonstrated,
        "a smoke must not mint a 3D broadband claim"
    );
    let bad = SlotJet3dFollowUp { nx: 3, ..spec };
    assert!(bad.validate().is_err());
}
