//! l011o: 2D tonal-lock refusal, broadband admission gate, 3D operator smoke.

use fs_aeroac::regime::{
    PINNED_2D_CENTRAL_MOMENT_TONAL, SlotJet3dFollowUp, SpectrumClass, TONAL_FLATNESS_CEILING,
    TWO_D_INVERSE_CASCADE, admit_broadband_spectrum, classify_spectrum,
    evaluate_slot_jet_3d_operator, measure_spectral_flatness, two_d_broadband_refusal,
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

#[test]
fn spectral_flatness_classifies_tone_and_white_without_minting_a_table() {
    let mut tone = vec![1.0e-18; 64];
    tone[7] = 1.0;
    let f_tone = measure_spectral_flatness(&tone).expect("tone");
    assert!(f_tone < TONAL_FLATNESS_CEILING);
    assert!(matches!(
        classify_spectrum(&tone).expect("class"),
        SpectrumClass::Tonal { .. }
    ));

    let white = vec![1.0; 64];
    let f_white = measure_spectral_flatness(&white).expect("white");
    assert!((f_white - 1.0).abs() < 1.0e-12);
    assert!(matches!(
        classify_spectrum(&white).expect("class"),
        SpectrumClass::Broadband { .. }
    ));

    let spec = SlotJet3dFollowUp::minimal_central_moment();
    let smoke = evaluate_slot_jet_3d_operator(spec).expect("smoke");
    let still_tonal = smoke
        .incorporate_measured_spectrum(&tone, spec.broadband_flatness_floor)
        .expect("tone attach");
    assert!(!still_tonal.broadband_demonstrated);
    let white_ok = smoke
        .incorporate_measured_spectrum(&white, spec.broadband_flatness_floor)
        .expect("white attach");
    assert!(
        white_ok.broadband_demonstrated,
        "a measured white spectrum may raise the flag; a 3-D jet table is not invented"
    );
}
