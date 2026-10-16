//! Modal-velocity input bundle for the existing cold vector-fitting owner.
//! These fits stay continuous: no prewarp or Tustin coefficients enter an
//! augmented exact-ZOH observer. Producer uncertainty remains separate.
use super::{
    BroadbandRadiationAuthority, MAX_BROADBAND_FREQUENCIES, MAX_BROADBAND_INPUTS,
    MAX_BROADBAND_TRANSFER_SAMPLES,
};
use crate::modal_acoustic_time::MAX_TIME_DOMAIN_ACOUSTIC_MODES;
use fs_math::{c64::C64, det};
use fs_vfit::model::RationalModel;
use fs_vfit::{FitOptions, WeightPreset, vector_fit};

/// Maximum per-mode continuous fit order; no high-order model is truncated.
pub const MAX_MODAL_TRANSFER_ORDER: usize = 128;
/// Window half-length for the explicit fractional-delay approximation axis.
pub const MAX_MODAL_DELAY_KERNEL_HALF_LEN: usize = 4096;

/// Computed receiver pressure per generalized modal velocity, exp(-i*w*t).
#[derive(Clone, Debug)]
pub struct ModalTransferSamples {
    /// Declared computed operator/run identity, retained without authority upgrade.
    pub source_id: String,
    /// Strictly increasing physical angular frequencies [rad/s].
    pub omega_rad_s: Vec<f64>,
    /// Mode-major H [Pa s/(m sqrt(kg))], including the declared pure delay.
    pub per_mode: Vec<Vec<C64>>,
    /// Independent, disjoint frequencies and their complete modal references.
    pub held_out: Vec<(f64, Vec<C64>)>,
    /// Explicit propagation delay [s], removed before continuous fitting.
    pub propagation_delay_s: f64,
    /// Producer's maximum declared relative discretization estimate. It is
    /// recorded, not subtracted from or absorbed into the realization budget.
    pub producer_relative_error_estimate: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BroadbandModalControls {
    /// Complete continuous pole count per mode; no poles are truncated later.
    pub fit_order: usize,
    /// Bounded cold vector-fit relocation budget, at most 128 updates.
    pub fit_iterations: usize,
    /// Explicit sampled-delay approximation axis, at most 4096 taps per side.
    pub delay_kernel_half_len: usize,
    /// Maximum sampled relative error including the actual delay kernel.
    pub max_held_out_rel_error: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub enum BroadbandModalError {
    TooManyModes {
        found: usize,
        cap: usize,
    },
    Cap {
        what: &'static str,
        found: usize,
        cap: usize,
    },
    Invalid(&'static str),
    HeldOut {
        max_rel: f64,
        budget: f64,
    },
    Fit(fs_vfit::VfError),
    Modal(crate::modal_acoustic_time::ModalAcousticTimeError),
    Exact(fs_phs::PhsError),
    NonFinite {
        what: &'static str,
    },
}
impl core::fmt::Display for BroadbandModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "broadband modal observer: {self:?}")
    }
}
impl std::error::Error for BroadbandModalError {}

/// Sampled operator error is an estimate, not a rigorous continuum enclosure.
#[derive(Clone, Debug)]
pub struct BroadbandModalReport {
    /// Original producer identity, not an identity inferred by the receiver.
    pub source_id: String,
    /// Always EstimateOnly: withheld samples do not enclose the continuum.
    pub authority: BroadbandRadiationAuthority,
    /// Original producer's discretization estimate, separate from fit error.
    pub producer_relative_error_estimate: f64,
    /// Actual caller declaration used in this build.
    pub controls: BroadbandModalControls,
    /// Exact observer/mechanical clock [Hz].
    pub sample_rate_hz: u32,
    /// Maximum over independent, withheld frequency/mode pairs.
    pub held_out_max_rel_error: f64,
    /// Physical delay times the clock, retaining the fractional part.
    pub delay_samples: f64,
    /// Constant/no-delay branch preserves the original quadrature operation order.
    pub constant_narrow_band_compatibility: bool,
    /// Actual modal-input batches, each at most the original 256-input cap.
    pub input_batches: usize,
}
pub(crate) struct ContinuousModalFits {
    pub models: Vec<RationalModel>,
    pub constant: Option<Vec<C64>>,
    pub report: BroadbandModalReport,
}

/// Allocation-free admission of every retained mode and sample; 4096 modal
/// inputs are processed in <=256-input batches, not truncated to 256 modes.
pub fn admit_modal_samples(
    modes: usize,
    training: usize,
    held_out: usize,
) -> Result<(), BroadbandModalError> {
    if modes > MAX_TIME_DOMAIN_ACOUSTIC_MODES {
        return Err(BroadbandModalError::TooManyModes {
            found: modes,
            cap: MAX_TIME_DOMAIN_ACOUSTIC_MODES,
        });
    }
    if modes == 0 || training == 0 || held_out == 0 {
        return Err(BroadbandModalError::Invalid(
            "nonempty modes, training and independent held-out grids required",
        ));
    }
    for (what, found) in [
        ("training frequencies", training),
        ("held-out frequencies", held_out),
    ] {
        if found > MAX_BROADBAND_FREQUENCIES {
            return Err(BroadbandModalError::Cap {
                what,
                found,
                cap: MAX_BROADBAND_FREQUENCIES,
            });
        }
    }
    let cells = modes
        .checked_mul(training + held_out)
        .ok_or(BroadbandModalError::Invalid("transfer-cell shape overflow"))?;
    if cells > MAX_BROADBAND_TRANSFER_SAMPLES {
        return Err(BroadbandModalError::Cap {
            what: "modal transfer cells",
            found: cells,
            cap: MAX_BROADBAND_TRANSFER_SAMPLES,
        });
    }
    Ok(())
}
pub(crate) fn fit_modal_velocity_transfers(
    s: &ModalTransferSamples,
    c: BroadbandModalControls,
    rate: u32,
) -> Result<ContinuousModalFits, BroadbandModalError> {
    admit_modal_samples(s.per_mode.len(), s.omega_rad_s.len(), s.held_out.len())?;
    if rate == 0
        || s.source_id.trim().is_empty()
        || !s.propagation_delay_s.is_finite()
        || s.propagation_delay_s < 0.0
        || !s.producer_relative_error_estimate.is_finite()
        || s.producer_relative_error_estimate < 0.0
        || c.fit_order == 0
        || c.fit_order > MAX_MODAL_TRANSFER_ORDER
        || c.fit_iterations == 0
        || c.fit_iterations > 128
        || c.delay_kernel_half_len == 0
        || c.delay_kernel_half_len > MAX_MODAL_DELAY_KERNEL_HALF_LEN
        || !c.max_held_out_rel_error.is_finite()
        || c.max_held_out_rel_error <= 0.0
    {
        return Err(BroadbandModalError::Invalid(
            "explicit finite source, delay, producer estimate, fit/delay controls and error budget required",
        ));
    }
    let nyquist = core::f64::consts::PI * f64::from(rate);
    let valid_frequency = |w: f64| w.is_finite() && w > 0.0 && w < nyquist;
    if s.omega_rad_s.iter().any(|w| !valid_frequency(*w))
        || s.omega_rad_s.windows(2).any(|w| w[0] >= w[1])
        || s.held_out.iter().any(|(w, row)| {
            !valid_frequency(*w)
                || row.len() != s.per_mode.len()
                || row.iter().any(|v| !v.re.is_finite() || !v.im.is_finite())
                || s.omega_rad_s.binary_search_by(|x| x.total_cmp(w)).is_ok()
        })
        || s.held_out.windows(2).any(|w| w[0].0 >= w[1].0)
        || s.per_mode.iter().any(|row| {
            row.len() != s.omega_rad_s.len()
                || row.iter().any(|v| !v.re.is_finite() || !v.im.is_finite())
        })
    {
        return Err(BroadbandModalError::Invalid(
            "finite complete modal samples on positive increasing disjoint sub-Nyquist grids required",
        ));
    }
    let low = s.omega_rad_s[0];
    let high = *s
        .omega_rad_s
        .last()
        .ok_or(BroadbandModalError::Invalid("empty grid"))?;
    if s.held_out.iter().any(|(w, _)| *w < low || *w > high) {
        return Err(BroadbandModalError::Invalid(
            "held-out frequencies must lie inside the training band",
        ));
    }
    let constant: Vec<C64> = s.per_mode.iter().map(|row| row[0]).collect();
    let legacy = s.propagation_delay_s == 0.0
        && s.per_mode.iter().enumerate().all(|(m, row)| {
            row.iter().all(|h| *h == constant[m])
                && s.held_out.iter().all(|(_, h)| h[m] == constant[m])
        });
    let report = BroadbandModalReport {
        source_id: s.source_id.clone(),
        authority: BroadbandRadiationAuthority::EstimateOnly,
        producer_relative_error_estimate: s.producer_relative_error_estimate,
        controls: c,
        sample_rate_hz: rate,
        held_out_max_rel_error: 0.0,
        delay_samples: s.propagation_delay_s * f64::from(rate),
        constant_narrow_band_compatibility: legacy,
        input_batches: s.per_mode.len().div_ceil(MAX_BROADBAND_INPUTS),
    };
    if legacy {
        return Ok(ContinuousModalFits {
            models: vec![],
            constant: Some(constant),
            report,
        });
    }
    let opts = FitOptions {
        order: c.fit_order,
        iterations: c.fit_iterations,
        weights: WeightPreset::InverseMagnitude,
        fit_e: false,
        fit_d: true,
    };
    let mut models = Vec::with_capacity(s.per_mode.len());
    for batch in s.per_mode.chunks(MAX_BROADBAND_INPUTS) {
        for row in batch {
            // H^- uses exp(-iwt); remove exp(+iwD), then conjugate for s=+iw.
            let response: Vec<_> = row
                .iter()
                .zip(&s.omega_rad_s)
                .map(|(&h, &w)| {
                    (h * C64::new(
                        det::cos(w * s.propagation_delay_s),
                        -det::sin(w * s.propagation_delay_s),
                    ))
                    .conj()
                })
                .collect();
            if response.iter().all(|h| *h == C64::ZERO) {
                models.push(RationalModel {
                    terms: vec![],
                    d: 0.0,
                    e: 0.0,
                });
                continue;
            }
            let fit =
                vector_fit(&s.omega_rad_s, &response, &opts).map_err(BroadbandModalError::Fit)?;
            if !fit.model.is_stable() || fit.model.e != 0.0 || !fit.model.d.is_finite() {
                return Err(BroadbandModalError::Invalid(
                    "observer fit must be stable, real and proper",
                ));
            }
            models.push(fit.model);
        }
    }
    Ok(ContinuousModalFits {
        models,
        constant: None,
        report,
    })
}
