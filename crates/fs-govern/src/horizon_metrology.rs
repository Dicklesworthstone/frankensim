//! Horizon trigger 11 (bead `frankensim-epic-addendum-xpck.5.1`): the
//! metrology-partnership and reality-chart registration-uncertainty gate for
//! Proposal 11 ("Reality as a chart", [`crate::proposals`]).
//!
//! Two conjunctive activation premises are evaluated:
//! 1. An authorized metrology partnership or retained real-data agreement exists
//!    with genuine non-synthetic provenance.
//! 2. Independently measured registration uncertainty ($u_{\text{reg}}$ at 95%
//!    confidence) on realistic scanned parts is demonstrably below the proposed
//!    geometric deviation threshold ($\delta_{\text{geom}}$) under active calibration.
//!
//! If either premise fails, the declared point-sensor assimilation fallback is
//! retained, and the disposition is [`MetrologyDisposition::Defer`] or
//! [`MetrologyDisposition::NoData`] rather than activating reality-chart authority.

use fs_blake3::hash_bytes;

/// Minimum margin by which registration uncertainty must sit below the geometric tolerance.
/// $u_{\text{reg}} < \delta_{\text{geom}}$ ($u_{\text{reg}} / \delta_{\text{geom}} < 1.0$).
pub const MAX_UNCERTAINTY_RATIO: f64 = 1.0;

/// Metrology partnership or retained-data agreement attestation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MetrologyAgreement {
    /// Partner institution or facility.
    pub partner: String,
    /// Formal agreement or contract reference.
    pub agreement_ref: String,
    /// License or data-rights term.
    pub license_terms: String,
    /// Content-addressed root of retained raw scan data.
    pub raw_data_root_hex: String,
    /// True if this is a synthetic mock (forbidden for live activation).
    pub is_synthetic: bool,
}

/// Scanned part specimen metadata and intended tolerance threshold.
#[derive(Debug, Clone, PartialEq)]
pub struct ScannedPartSpecimen {
    /// Unique specimen label.
    pub specimen_id: String,
    /// Intended claim class (e.g. "GD&T-ToleranceAllocation", "AsBuilt-Deformation").
    pub claim_class: String,
    /// Proposed geometric deviation / tolerance threshold [m] ($\delta_{\text{geom}}$).
    pub geometric_tolerance_m: f64,
    /// Reference coordinate frame name.
    pub coordinate_frame: String,
}

/// Independently measured registration uncertainty and metrology chain state.
#[derive(Debug, Clone, PartialEq)]
pub struct RegistrationUncertainty {
    /// Registration / scan method (e.g. "optical-tracker-icp", "sub-voxel-ct").
    pub method: String,
    /// Measurement coordinate frame name (must match specimen frame).
    pub coordinate_frame: String,
    /// True if calibration is current and traceable.
    pub is_calibrated: bool,
    /// Calibration validity attestation.
    pub calibration_valid: bool,
    /// 95th percentile registration uncertainty [m] ($u_{\text{reg}}$).
    pub uncertainty_95_m: f64,
}

/// Input premises for Proposal 11 activation decision.
#[derive(Debug, Clone, PartialEq)]
pub struct Proposal11Premises {
    /// Optional authorized partnership agreement.
    pub agreement: Option<MetrologyAgreement>,
    /// Scanned specimen definition.
    pub specimen: ScannedPartSpecimen,
    /// Measured registration uncertainty.
    pub registration: RegistrationUncertainty,
    /// Whether point-sensor assimilation fallback is declared and operational.
    pub point_sensor_fallback_active: bool,
}

/// Typed refusals for malformed or inadmissible premises.
#[derive(Debug, Clone, PartialEq)]
pub enum Trigger11Refusal {
    /// Agreement is absent.
    MissingAgreement,
    /// Synthetic data cannot satisfy an empirical activation gate.
    SyntheticDataForbidden,
    /// Instrument is uncalibrated or calibration is invalid.
    CalibrationInvalid,
    /// Tolerance threshold is non-positive or non-finite.
    NonPositiveTolerance { val: f64 },
    /// Registration uncertainty is non-positive or non-finite.
    NonPositiveUncertainty { val: f64 },
    /// Specimen and registration coordinate frames do not match.
    FrameMismatch { specimen: String, registration: String },
    /// Raw data root is empty or invalid.
    EmptyDataRoot,
}

/// Activation verdict for Proposal 11.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trigger11Verdict {
    /// Both premises satisfied: reality-chart authority activated.
    Activate,
    /// Premise(s) not satisfied: retain point-sensor fallback.
    FallbackPointSensors,
}

/// Overall population disposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetrologyDisposition {
    /// Ready for live reality-chart activation.
    Activate,
    /// Evaluated and deferred (premises not met; fallback active).
    Defer,
    /// No real data or partnership present yet.
    NoData,
}

/// Immutable receipt of a Trigger 11 evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct Trigger11Receipt {
    pub proposal: &'static str,
    pub disposition: MetrologyDisposition,
    pub verdict: Trigger11Verdict,
    pub uncertainty_ratio: f64,
    pub point_sensor_fallback_retained: bool,
    pub receipt_hash: String,
    pub reason: String,
}

/// Evaluate Proposal 11 activation premises.
///
/// # Errors
/// Returns [`Trigger11Refusal`] if any input parameter or invariant is violated.
pub fn evaluate_trigger_11(premises: &Proposal11Premises) -> Result<Trigger11Verdict, Trigger11Refusal> {
    let Some(agreement) = &premises.agreement else {
        return Err(Trigger11Refusal::MissingAgreement);
    };
    if agreement.is_synthetic {
        return Err(Trigger11Refusal::SyntheticDataForbidden);
    }
    if agreement.raw_data_root_hex.trim().is_empty() {
        return Err(Trigger11Refusal::EmptyDataRoot);
    }
    if !premises.registration.is_calibrated || !premises.registration.calibration_valid {
        return Err(Trigger11Refusal::CalibrationInvalid);
    }
    if !premises.specimen.geometric_tolerance_m.is_finite() || premises.specimen.geometric_tolerance_m <= 0.0 {
        return Err(Trigger11Refusal::NonPositiveTolerance { val: premises.specimen.geometric_tolerance_m });
    }
    if !premises.registration.uncertainty_95_m.is_finite() || premises.registration.uncertainty_95_m <= 0.0 {
        return Err(Trigger11Refusal::NonPositiveUncertainty { val: premises.registration.uncertainty_95_m });
    }
    if premises.specimen.coordinate_frame != premises.registration.coordinate_frame {
        return Err(Trigger11Refusal::FrameMismatch {
            specimen: premises.specimen.coordinate_frame.clone(),
            registration: premises.registration.coordinate_frame.clone(),
        });
    }

    let ratio = premises.registration.uncertainty_95_m / premises.specimen.geometric_tolerance_m;
    if ratio < MAX_UNCERTAINTY_RATIO {
        Ok(Trigger11Verdict::Activate)
    } else {
        Ok(Trigger11Verdict::FallbackPointSensors)
    }
}

/// Mint an immutable decision receipt for Proposal 11.
#[must_use]
pub fn mint_trigger_11_receipt(
    premises_opt: Option<&Proposal11Premises>,
) -> Trigger11Receipt {
    let Some(premises) = premises_opt else {
        let hash = hash_bytes(b"org.frankensim.horizon-trigger-11.nodata.v1").to_hex();
        return Trigger11Receipt {
            proposal: "11",
            disposition: MetrologyDisposition::NoData,
            verdict: Trigger11Verdict::FallbackPointSensors,
            uncertainty_ratio: f64::NAN,
            point_sensor_fallback_retained: true,
            receipt_hash: hash,
            reason: "no authorized metrology partnership or retained scan data exists in the program yet".into(),
        };
    };

    match evaluate_trigger_11(premises) {
        Ok(Trigger11Verdict::Activate) => {
            let ratio = premises.registration.uncertainty_95_m / premises.specimen.geometric_tolerance_m;
            let mut payload = Vec::new();
            payload.extend_from_slice(b"org.frankensim.horizon-trigger-11.activate.v1");
            payload.extend_from_slice(ratio.to_le_bytes().as_slice());
            if let Some(ag) = &premises.agreement {
                payload.extend_from_slice(ag.raw_data_root_hex.as_bytes());
            }
            let hash = hash_bytes(&payload).to_hex();
            Trigger11Receipt {
                proposal: "11",
                disposition: MetrologyDisposition::Activate,
                verdict: Trigger11Verdict::Activate,
                uncertainty_ratio: ratio,
                point_sensor_fallback_retained: false,
                receipt_hash: hash,
                reason: format!("registration uncertainty ({:.2e} m) is strictly below geometric tolerance ({:.2e} m; ratio {:.3}) under authorized agreement", premises.registration.uncertainty_95_m, premises.specimen.geometric_tolerance_m, ratio),
            }
        }
        Ok(Trigger11Verdict::FallbackPointSensors) => {
            let ratio = premises.registration.uncertainty_95_m / premises.specimen.geometric_tolerance_m;
            let hash = hash_bytes(b"org.frankensim.horizon-trigger-11.defer.v1").to_hex();
            Trigger11Receipt {
                proposal: "11",
                disposition: MetrologyDisposition::Defer,
                verdict: Trigger11Verdict::FallbackPointSensors,
                uncertainty_ratio: ratio,
                point_sensor_fallback_retained: true,
                receipt_hash: hash,
                reason: format!("registration uncertainty ({:.2e} m) exceeds or equals tolerance ({:.2e} m; ratio {:.3}); retaining point-sensor fallback", premises.registration.uncertainty_95_m, premises.specimen.geometric_tolerance_m, ratio),
            }
        }
        Err(refusal) => {
            let hash = hash_bytes(format!("{:?}", refusal).as_bytes()).to_hex();
            Trigger11Receipt {
                proposal: "11",
                disposition: MetrologyDisposition::Defer,
                verdict: Trigger11Verdict::FallbackPointSensors,
                uncertainty_ratio: f64::NAN,
                point_sensor_fallback_retained: true,
                receipt_hash: hash,
                reason: format!("inadmissible premises ({refusal:?}); retaining point-sensor fallback"),
            }
        }
    }
}
