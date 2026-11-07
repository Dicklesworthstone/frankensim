//! Explicit Walker mean-stress correction for the original S-N reference R.
//!
//! Equating Walker equivalent stress at the observed and reference ratios gives
//! a_ref = a * ((1 - R_ref) / (1 - R))^(1 - gamma).
//! See NASA/TM-2005-213886/REV2, NASALIFE, section 2.3, equations (23)-(25):
//! <https://ntrs.nasa.gov/citations/20140010774>.
//!
//! Gamma and its load-ratio support are supplied calibration, never inferred
//! from tensile strength or alloy names. Only the shared S-N finite-life domain
//! is used after correction. This is an elastic stress-life estimate, not an
//! elastoplastic strain-life model, static-strength check or crack-growth law.

use std::sync::Arc;

use fs_matdb::{
    ClaimSelection, MaterialAnswer, MaterialCard, PropertyClaim, PropertyKey, PropertyValue,
    QueryPoint,
};
use fs_math::det;
use fs_qty::{Dims, QuantitySpec, Temperature};

use super::source::{SnSourceError, SnSourceRequest, resolve_log_log_sn_at_stress_ratio};
use super::{FatigueError, MeanStressCorrection, SnCurve};

/// One explicitly calibrated Walker exponent and its CLOSED load-ratio domain.
/// Gamma=1 neglects mean stress; gamma=0 uses maximum stress; gamma=1/2 gives
/// the elastic Smith-Watson-Topper equivalent amplitude when R_ref=-1.
/// This first policy admits -1 <= R < 1 and gives no compressive-mean life credit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WalkerParameters {
    pub gamma: f64,
    pub minimum_load_ratio: f64,
    pub maximum_load_ratio: f64,
}

impl WalkerParameters {
    pub(crate) fn validate(self) -> Result<(), FatigueError> {
        if !self.gamma.is_finite()
            || !(0.0..=1.0).contains(&self.gamma)
            || !self.minimum_load_ratio.is_finite()
            || !self.maximum_load_ratio.is_finite()
            || self.minimum_load_ratio < -1.0
            || self.minimum_load_ratio >= self.maximum_load_ratio
            || self.maximum_load_ratio >= 1.0
        {
            return Err(FatigueError::InvalidCalibration(
                "Walker needs finite 0 <= gamma <= 1 and -1 <= minimum R < maximum R < 1",
            ));
        }
        Ok(())
    }

    pub(crate) fn validate_reference(self, reference_load_ratio: f64) -> Result<(), FatigueError> {
        self.validate()?;
        if !reference_load_ratio.is_finite()
            || !(-1.0..1.0).contains(&reference_load_ratio)
            || reference_load_ratio < self.minimum_load_ratio
            || reference_load_ratio > self.maximum_load_ratio
        {
            return Err(FatigueError::InvalidCalibration(
                "the original S-N reference R must lie inside the Walker calibration domain",
            ));
        }
        Ok(())
    }

    /// Amplitude at the supplied S-N reference ratio producing the same Walker
    /// parameter. Both cycle extrema and the resulting amplitude must remain
    /// representable. Zero amplitude is an observation with no fatigue cycle.
    ///
    /// Normalization avoids multiplying powers of dimensional stresses.
    /// Each positive-amplitude cycle must have a tensile maximum and remain in
    /// the declared R support, including when gamma=1. Only roundoff tolerance
    /// is allowed at domain boundaries; no extrapolation or clamping is used.
    pub fn equivalent_amplitude(
        self,
        amplitude_pa: f64,
        mean_pa: f64,
        reference_load_ratio: f64,
    ) -> Result<f64, FatigueError> {
        self.validate_reference(reference_load_ratio)?;
        if !amplitude_pa.is_finite() || amplitude_pa < 0.0 || !mean_pa.is_finite() {
            return Err(FatigueError::Arithmetic);
        }
        if amplitude_pa == 0.0 {
            return Ok(0.0);
        }
        let scale = amplitude_pa.max(mean_pa.abs());
        let amplitude = amplitude_pa / scale;
        let mean = mean_pa / scale;
        let maximum = mean + amplitude;
        let minimum = mean - amplitude;
        if maximum <= 0.0 {
            return Err(FatigueError::NoTensileStress {
                maximum_pa: mean_pa + amplitude_pa,
            });
        }
        if amplitude == 0.0 || maximum > f64::MAX / scale || -minimum > f64::MAX / scale {
            return Err(FatigueError::Arithmetic);
        }
        // Compare 1-R = 2a/sigma_max directly. Subtracting nearly equal
        // extrema, or an absolute tolerance in R, loses the supported range
        // when R approaches one and can hide a large amplitude correction.
        let complement = 2.0 * amplitude / maximum;
        let lower = 1.0 - self.maximum_load_ratio;
        let upper = 1.0 - self.minimum_load_ratio;
        let below = complement < lower && lower - complement > 32.0 * f64::EPSILON * lower;
        let above = complement > upper && complement - upper > 32.0 * f64::EPSILON * upper;
        if !complement.is_finite() || below || above {
            return Err(FatigueError::StressRatioOutsideWalkerDomain {
                actual_ratio: 1.0 - complement,
                minimum_ratio: self.minimum_load_ratio,
                maximum_ratio: self.maximum_load_ratio,
            });
        }
        if self.gamma == 1.0 {
            return Ok(amplitude_pa);
        }
        let reference_peak = 0.5 * (1.0 - reference_load_ratio) * maximum;
        // Preserve original knots only when the possible AMPLITUDE correction
        // is itself roundoff-sized, including at reference ratios near one.
        if (reference_peak - amplitude).abs() <= 32.0 * f64::EPSILON * amplitude {
            return Ok(amplitude_pa);
        }
        let normalized = if self.gamma == 0.0 {
            reference_peak
        } else if self.gamma == 0.5 {
            det::sqrt(amplitude) * det::sqrt(reference_peak)
        } else {
            det::exp(self.gamma * det::ln(amplitude) + (1.0 - self.gamma) * det::ln(reference_peak))
        };
        let result = scale * normalized;
        if !result.is_finite() || result <= 0.0 {
            return Err(FatigueError::Arithmetic);
        }
        Ok(result)
    }
}

/// Deliberately distinct from the crack-growth Walker gamma property.
pub fn walker_gamma_property_key() -> PropertyKey {
    PropertyKey::new("fatigue-walker-gamma", Dims::NONE)
}

/// Original gamma claim and answer attached to the SAME immutable S-N curve.
/// Card, state, full typed query and S-N observations remain in curve.source().
#[derive(Debug, Clone, PartialEq)]
pub struct WalkerSnCalibrationSource {
    parameters: WalkerParameters,
    claim: PropertyClaim,
    answer: MaterialAnswer,
    selection: ClaimSelection,
}
impl WalkerSnCalibrationSource {
    pub const fn parameters(&self) -> WalkerParameters {
        self.parameters
    }
    pub fn claim(&self) -> &PropertyClaim {
        &self.claim
    }
    pub fn answer(&self) -> &MaterialAnswer {
        &self.answer
    }
    pub const fn selection(&self) -> ClaimSelection {
        self.selection
    }
}

fn require(condition: bool, reason: &'static str) -> Result<(), SnSourceError> {
    if condition {
        Ok(())
    } else {
        Err(SnSourceError::Invalid(reason))
    }
}

/// Resolve original fixed-R S-N observations and a separately named Walker
/// gamma fit on the same card. The returned policy must travel with the curve.
///
/// Gamma may be a fixed scalar or an original temperature-only curve. Its
/// provenance and observation population must match the original S-N claim.
/// The complete typed validity domain must match too, except that gamma declares
/// its finite R interval while the S-N table retains its exact reference R.
/// This admits no gamma dependence on life, R, pressure or other coordinates;
/// a rectangular validity box cannot stand in for such a fitted surface.
///
/// Every original S-N knot and gamma receipt remains inspectable. A counter
/// refuses a different gamma or R domain when this source binding is present.
pub fn resolve_walker_sn(
    card: &MaterialCard,
    point: &QueryPoint,
    request: &SnSourceRequest,
    sn_selection: ClaimSelection,
    gamma_selection: ClaimSelection,
) -> Result<(SnCurve, MeanStressCorrection), SnSourceError> {
    let mut curve = resolve_log_log_sn_at_stress_ratio(card, point, request, sn_selection)?;
    let key = walker_gamma_property_key();
    let answer = match gamma_selection {
        ClaimSelection::Policy(policy) => card.claims().query_typed(&key, point, policy)?,
        ClaimSelection::Pinned(pin) => card.claims().query_pinned_typed(&key, point, pin)?,
    };
    card.claims().verify_receipt(&answer.receipt)?;
    let claim = card
        .claims()
        .claim(answer.receipt.selected)
        .ok_or(SnSourceError::Invalid(
            "original S-N Walker gamma claim is absent",
        ))?;
    let dimensionless = QuantitySpec::dimensional(Dims::NONE);
    require(
        claim
            .validity
            .axis_quantities()
            .get(&request.stress_ratio_axis)
            == Some(&dimensionless),
        "Walker R support must retain its explicit dimensionless coordinate",
    )?;
    let &(minimum_load_ratio, maximum_load_ratio) = claim
        .validity
        .bounds()
        .get(&request.stress_ratio_axis)
        .ok_or(SnSourceError::Invalid(
            "Walker gamma needs an explicit load-ratio support interval",
        ))?;
    require(
        match &claim.value {
            PropertyValue::Scalar { .. } => true,
            PropertyValue::Curve {
                abscissa,
                abscissa_dims,
                ..
            } => abscissa == &request.temperature_axis && *abscissa_dims == Temperature::DIMS,
            _ => false,
        },
        "Walker gamma must be a fixed scalar or a function only of the original temperature",
    )?;
    let sn = curve.source().ok_or(SnSourceError::Invalid(
        "Walker needs the original source-resolved S-N observations",
    ))?;
    let at_reference = claim.validity.clone().with_quantity(
        &request.stress_ratio_axis,
        dimensionless,
        curve.reference_stress_ratio(),
        curve.reference_stress_ratio(),
    );
    require(
        claim.provenance == sn.claim().provenance
            && claim.observations == sn.claim().observations
            && &at_reference == &sn.claim().validity,
        "Walker gamma and S-N data must share provenance, population and full typed support except R",
    )?;
    let parameters = WalkerParameters {
        gamma: answer.evidence.value.value,
        minimum_load_ratio,
        maximum_load_ratio,
    };
    parameters.validate_reference(curve.reference_stress_ratio())?;
    curve.walker_source = Some(Arc::new(WalkerSnCalibrationSource {
        parameters,
        claim: claim.clone(),
        answer,
        selection: gamma_selection,
    }));
    Ok((curve, MeanStressCorrection::Walker(parameters)))
}
