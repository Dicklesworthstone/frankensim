//! Horizon trigger 11 battery (bead `frankensim-epic-addendum-xpck.5.1`):
//! boundary, mutant, refusal, and receipt tests for Proposal 11 metrology
//! partnership and reality-chart registration uncertainty gate.

use fs_govern::horizon_metrology::{
    evaluate_trigger_11, mint_trigger_11_receipt, MetrologyAgreement, MetrologyDisposition,
    Proposal11Premises, RegistrationUncertainty, ScannedPartSpecimen, Trigger11Receipt,
    Trigger11Refusal, Trigger11Verdict, MAX_UNCERTAINTY_RATIO,
};

fn valid_agreement() -> MetrologyAgreement {
    MetrologyAgreement {
        partner: "National Metrology Institute (NIST/PTB)".into(),
        agreement_ref: "AGR-2026-MET-0042".into(),
        license_terms: "CC-BY-4.0 Open-Metrology-Data".into(),
        raw_data_root_hex: "d41d8cd98f00b204e9800998ecf8427e".into(),
        is_synthetic: false,
    }
}

fn valid_specimen() -> ScannedPartSpecimen {
    ScannedPartSpecimen {
        specimen_id: "coldplate-manifold-specimen-01".into(),
        claim_class: "GD&T-ToleranceAllocation".into(),
        geometric_tolerance_m: 50e-6, // 50 microns
        coordinate_frame: "specimen_cad_origin_v1".into(),
    }
}

fn valid_registration(uncertainty_m: f64) -> RegistrationUncertainty {
    RegistrationUncertainty {
        method: "optical-tracker-fiducial-icp".into(),
        coordinate_frame: "specimen_cad_origin_v1".into(),
        is_calibrated: true,
        calibration_valid: true,
        uncertainty_95_m: uncertainty_m,
    }
}

#[test]
fn gate_activates_when_agreement_is_genuine_and_uncertainty_is_below_tolerance() {
    let premises = Proposal11Premises {
        agreement: Some(valid_agreement()),
        specimen: valid_specimen(),
        registration: valid_registration(10e-6), // 10 microns < 50 microns
        point_sensor_fallback_active: true,
    };
    assert_eq!(evaluate_trigger_11(&premises), Ok(Trigger11Verdict::Activate));

    let receipt = mint_trigger_11_receipt(Some(&premises));
    assert_eq!(receipt.disposition, MetrologyDisposition::Activate);
    assert_eq!(receipt.verdict, Trigger11Verdict::Activate);
    assert!(!receipt.point_sensor_fallback_retained);
    assert!(receipt.uncertainty_ratio < MAX_UNCERTAINTY_RATIO);
}

#[test]
fn gate_defers_when_uncertainty_exceeds_tolerance() {
    let premises = Proposal11Premises {
        agreement: Some(valid_agreement()),
        specimen: valid_specimen(), // 50 microns
        registration: valid_registration(60e-6), // 60 microns > 50 microns
        point_sensor_fallback_active: true,
    };
    assert_eq!(evaluate_trigger_11(&premises), Ok(Trigger11Verdict::FallbackPointSensors));

    let receipt = mint_trigger_11_receipt(Some(&premises));
    assert_eq!(receipt.disposition, MetrologyDisposition::Defer);
    assert_eq!(receipt.verdict, Trigger11Verdict::FallbackPointSensors);
    assert!(receipt.point_sensor_fallback_retained);
    assert!(receipt.uncertainty_ratio > MAX_UNCERTAINTY_RATIO);
}

#[test]
fn gate_refuses_synthetic_data() {
    let mut ag = valid_agreement();
    ag.is_synthetic = true; // Synthetic data trying to masquerade as empirical
    let premises = Proposal11Premises {
        agreement: Some(ag),
        specimen: valid_specimen(),
        registration: valid_registration(5e-6),
        point_sensor_fallback_active: true,
    };
    assert_eq!(
        evaluate_trigger_11(&premises),
        Err(Trigger11Refusal::SyntheticDataForbidden)
    );
}

#[test]
fn gate_refuses_missing_agreement() {
    let premises = Proposal11Premises {
        agreement: None,
        specimen: valid_specimen(),
        registration: valid_registration(5e-6),
        point_sensor_fallback_active: true,
    };
    assert_eq!(
        evaluate_trigger_11(&premises),
        Err(Trigger11Refusal::MissingAgreement)
    );
}

#[test]
fn gate_refuses_uncalibrated_instrument() {
    let mut reg = valid_registration(5e-6);
    reg.is_calibrated = false;
    let premises = Proposal11Premises {
        agreement: Some(valid_agreement()),
        specimen: valid_specimen(),
        registration: reg,
        point_sensor_fallback_active: true,
    };
    assert_eq!(
        evaluate_trigger_11(&premises),
        Err(Trigger11Refusal::CalibrationInvalid)
    );
}

#[test]
fn gate_refuses_frame_mismatch() {
    let mut reg = valid_registration(5e-6);
    reg.coordinate_frame = "optical_scanner_local_raw".into(); // Frame mismatch
    let premises = Proposal11Premises {
        agreement: Some(valid_agreement()),
        specimen: valid_specimen(),
        registration: reg,
        point_sensor_fallback_active: true,
    };
    assert!(matches!(
        evaluate_trigger_11(&premises),
        Err(Trigger11Refusal::FrameMismatch { .. })
    ));
}

#[test]
fn gate_refuses_non_positive_quantities() {
    let mut spec = valid_specimen();
    spec.geometric_tolerance_m = -1e-5;
    let premises = Proposal11Premises {
        agreement: Some(valid_agreement()),
        specimen: spec,
        registration: valid_registration(5e-6),
        point_sensor_fallback_active: true,
    };
    assert!(matches!(
        evaluate_trigger_11(&premises),
        Err(Trigger11Refusal::NonPositiveTolerance { .. })
    ));
}

#[test]
fn mint_receipt_returns_nodata_when_premises_absent() {
    let receipt: Trigger11Receipt = mint_trigger_11_receipt(None);
    assert_eq!(receipt.disposition, MetrologyDisposition::NoData);
    assert_eq!(receipt.verdict, Trigger11Verdict::FallbackPointSensors);
    assert!(receipt.point_sensor_fallback_retained);
    assert!(receipt.uncertainty_ratio.is_nan());
}
