//! Fit declared excess boundary resistance to the existing passive Foster model.
//!
//! The fit acts on Re Z in excess of the caller's base R-L-C, not on an output
//! waveform. fs-phs owns coefficient identification; fs-vfit owns its persistent
//! wave-port realization. A positive fallback from the fitter is NOT sufficient
//! for acceptance: every training and separate check sample must meet the
//! caller's error allowance. The dispersive imaginary part is implied by this
//! causal model, not fitted or validated by real-resistance samples. No material
//! source authority or continuous-band error guarantee is inferred.

use crate::acoustic_realize::AcousticRealizeError;
use super::network::NetworkNode;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{RelaxationImpedanceSpec, RelaxationTerm, MAX_RELAXATION_BRANCHES};

/// An admitted finite-sample fit, not a certificate for an entire frequency band.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FittedBoundaryLoss {
    /// Base R-L-C plus the identified excess-loss terms, in original pole order.
    pub load: RelaxationImpedanceSpec,
    /// Training-band endpoints [rad/s]; no extrapolation accuracy is claimed.
    pub angular_frequency_range: [f64; 2],
    /// Maximum observed absolute Re Z discrepancy over training AND check rows.
    pub max_checked_error_pa_s_m3: f64,
}
impl FittedBoundaryLoss {
    /// Terminal consumed directly by the existing stateful ApertureNetwork.
    #[must_use]
    pub const fn termination(&self) -> NetworkNode { NetworkNode::Relaxation { load: self.load } }
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

/// Identify up to eight passive relaxation branches using the existing fitter.
/// Each row is (angular frequency [rad/s], EXCESS Re Z [Pa s/m^3]). A cavity's
/// constant resistance must not be included again in these excess samples.
/// Training frequencies are also the admitted poles; check rows must be
/// distinct from them, ordered, and inside their span. All samples are finite.
///
/// Admission at each row is |prediction-target| <= absolute + relative*target.
/// The branch count fixes memory/work before runtime construction. This models
/// boundary loss only; the connected tube sections remain lossless.
///
/// # Errors
/// Invalid table/tolerances, identification failure, or any sample outside its
/// error allowance (including a positive but inaccurate fallback from fs-phs).
pub fn fit_boundary_loss(
    base: SeriesImpedanceSpec,
    training: &[(f64, f64)],
    checks: &[(f64, f64)],
    absolute_tolerance_pa_s_m3: f64,
    relative_tolerance: f64,
) -> Result<FittedBoundaryLoss, AcousticRealizeError> {
    base.validate().map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    if !(2..=MAX_RELAXATION_BRANCHES).contains(&training.len())
        || checks.is_empty() || checks.len() > 1024
        || !absolute_tolerance_pa_s_m3.is_finite() || absolute_tolerance_pa_s_m3 < 0.0
        || !relative_tolerance.is_finite() || !(0.0..=1.0).contains(&relative_tolerance)
    {
        return Err(invalid("boundary loss needs 2..=8 training rows, 1..=1024 checks and finite explicit tolerances"));
    }
    for samples in [training, checks] {
        if samples.iter().any(|&(w, r)| !w.is_finite() || w <= 0.0 || !r.is_finite() || r < 0.0)
            || samples.windows(2).any(|p| p[0].0 >= p[1].0)
        {
            return Err(invalid("boundary loss samples require increasing finite frequencies and nonnegative excess resistance"));
        }
    }
    let range = [training[0].0, training[training.len() - 1].0];
    if checks.iter().any(|&(w, _)| w < range[0] || w > range[1] || training.iter().any(|&(t, _)| w == t)) {
        return Err(invalid("boundary loss checks must be separate in-band frequencies, not reused training rows"));
    }
    // Reject overflow/underflow in the owner's squared-frequency matrix before
    // invoking it. Do not let its positive fallback hide an invalid table scale.
    if training.iter().any(|&(w, _)| !((w * w) * 2.0).is_finite() || w * w == 0.0) {
        return Err(invalid("boundary loss frequencies exceed the fitter's representable scale"));
    }
    let fitted = fs_phs::foster_match_re(training, training.len())
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let mut terms = [RelaxationTerm::default(); MAX_RELAXATION_BRANCHES];
    for (slot, &(r, w)) in terms.iter_mut().zip(&fitted) {
        *slot = RelaxationTerm { resistance_pa_s_m3: r, rate_per_s: w };
    }
    let load = RelaxationImpedanceSpec::new(base, &terms[..training.len()])
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let mut max_error = 0.0_f64;
    for &(w, target) in training.iter().chain(checks) {
        let mut predicted = 0.0;
        for term in load.terms() {
            let scale = w.max(term.rate_per_s);
            let x = w / scale;
            let p = term.rate_per_s / scale;
            predicted += term.resistance_pa_s_m3 * (x * x / (x * x + p * p));
        }
        let error = (predicted - target).abs();
        let allowance = absolute_tolerance_pa_s_m3 + relative_tolerance * target;
        if !predicted.is_finite() || !allowance.is_finite() || error > allowance {
            return Err(invalid("passive boundary loss fit exceeds the declared training or check error allowance"));
        }
        max_error = max_error.max(error);
    }
    Ok(FittedBoundaryLoss { load, angular_frequency_range: range, max_checked_error_pa_s_m3: max_error })
}
