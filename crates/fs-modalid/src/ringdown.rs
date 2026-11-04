//! Map the shared time-domain damped-pole fit into physical modal tables.
//!
//! A positive-imaginary Laplace pole is `p = -alpha + i*2*pi*f_d`.
//! Its natural frequency is `|p|/(2*pi)` and damping ratio is `-Re(p)/|p|`.
//! For the fit's peak `exp(-i*omega*t)` phasor `a+i*b`, the coefficient of
//! `exp(p*t)` in the real response is `(a-i*b)/2`. These observation residues
//! preserve sensor shape and phase, but depend on the initial excitation:
//! they are not force-normalized FRF residues or identified modal masses.

use crate::{ModalIdError, Mode, correct_exponential_window_pole};
use core::{f64::consts::TAU, mem::size_of};
use fs_fft::harmonic::{
    DampedToneBounds, DampedTonePole, DampedTonePoleChannelsFit, DampedToneRefinementBudget,
    DampedToneRefinementError, HarmonicFitError, refine_damped_tone_poles_channels,
};
use fs_math::c64::C64;

/// A local ringdown fit and its physical mode-table interpretation.
#[derive(Debug, Clone)]
pub struct RingdownIdentification {
    /// Modes sorted by corrected natural frequency, with sensor channels in
    /// input order. Both confidence-interval fields are NaN: none is estimated.
    pub modes: Vec<Mode>,
    /// For each sorted mode, its pole/phasor index in [`Self::fit`] and in the
    /// caller's seeds and bounds. Sorting never changes the raw fit's order.
    pub source_pole_indices: Vec<usize>,
    /// Unmodified fit of the supplied records, including any declared window.
    /// Poles, amplitudes, residuals, conditioning and corrections remain raw.
    pub fit: DampedTonePoleChannelsFit,
}

/// Typed refusal from the original fit or the physical interpretation.
#[derive(Debug, Clone, PartialEq)]
pub enum RingdownError {
    /// Fit admission, storage, conditioning or convergence refusal. Extra mode
    /// table storage is reserved first; a downstream budget error reports the
    /// remaining fit budget after that reservation.
    Fit(DampedToneRefinementError),
    /// Invalid window declaration or an unrepresentable physical modal pole.
    Mapping(ModalIdError),
}

impl core::fmt::Display for RingdownError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Fit(error) => write!(f, "ringdown fit: {error}"),
            Self::Mapping(error) => write!(f, "ringdown mode mapping: {error}"),
        }
    }
}

impl std::error::Error for RingdownError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fit(error) => Some(error),
            Self::Mapping(error) => Some(error),
        }
    }
}

impl From<DampedToneRefinementError> for RingdownError {
    fn from(error: DampedToneRefinementError) -> Self {
        Self::Fit(error)
    }
}

impl From<ModalIdError> for RingdownError {
    fn from(error: ModalIdError) -> Self {
        Self::Mapping(error)
    }
}

/// Identify physical modes from synchronized, calibrated free-decay channels.
///
/// All channels jointly refine the caller's damped-frequency/decay seeds
/// within the declared bounds. There is no channel averaging, harmonic-grid
/// rounding or independently fitted pole matching. Samples share one clock
/// and first-sample time; sensor units must be compatible for shape comparison.
///
/// `window_tau = Some(tau)` declares that the caller has already multiplied
/// the response by `exp(-t/tau)`. The existing pole-correction owner removes
/// that decay from the reported modes, retaining its passive clamp to zero
/// when the declared window explains more decay than the fit. Residues and
/// the raw fit are unchanged. No window is applied to the supplied records.
///
/// This is a bounded local fit at the supplied model order, not mode discovery,
/// an SNR gate, a stabilization ladder or a statistical uncertainty estimate.
/// Frequency and damping intervals are explicitly NaN. Unobserved, unresolved
/// or ill-conditioned poles are subject to the shared fit's typed refusals.
///
/// The fit's workspace budget also covers the extra mode, shape, index and
/// sorting storage. Caller-owned samples and allocator metadata are excluded,
/// as in [`fs_fft::harmonic::HarmonicFitBudget`].
///
/// # Errors
/// [`RingdownError::Fit`] preserves the underlying fit refusal, including an
/// insufficient storage budget. [`RingdownError::Mapping`] reports an invalid
/// window before any fit, or a nonrepresentable physical pole afterwards.
pub fn identify_ringdown(
    channels: &[&[f64]],
    timestep_s: f64,
    initial_poles: &[DampedTonePole],
    bounds: &[DampedToneBounds],
    mut budget: DampedToneRefinementBudget,
    window_tau: Option<f64>,
) -> Result<RingdownIdentification, RingdownError> {
    if let Some(tau) = window_tau {
        // Admit the declaration before fitting; this invokes the same owner
        // that will translate every fitted pole without constructing data.
        correct_exponential_window_pole(C64::new(0.0, 1.0), tau)?;
    }
    // Ordered tuples coexist with the final modes/indices while moving them
    // out. Shape buffers move, never clone. Use checked sizes before any of
    // these allocations or the much larger shared numerical fit is entered.
    let extra_bytes = channels
        .len()
        .checked_mul(size_of::<C64>())
        .and_then(|shape| {
            shape.checked_add(size_of::<(Mode, usize)>() + size_of::<Mode>() + size_of::<usize>())
        })
        .and_then(|per_mode| per_mode.checked_mul(initial_poles.len()));
    let storage_error = || {
        RingdownError::Fit(DampedToneRefinementError::Fit(
            HarmonicFitError::BudgetExceeded {
                required_bytes: extra_bytes.unwrap_or(usize::MAX),
                max_bytes: budget.fit.max_workspace_bytes,
            },
        ))
    };
    let extra_bytes = extra_bytes.ok_or_else(storage_error)?;
    budget.fit.max_workspace_bytes = budget
        .fit
        .max_workspace_bytes
        .checked_sub(extra_bytes)
        .ok_or_else(storage_error)?;
    let fit = refine_damped_tone_poles_channels(channels, timestep_s, initial_poles, bounds, budget)?;
    let mut ordered = Vec::with_capacity(fit.poles.len());
    for (source, pole) in fit.poles.iter().enumerate() {
        let measured = C64::new(-pole.decay_rate_per_s, TAU * pole.frequency_hz);
        let physical = match window_tau {
            Some(tau) => correct_exponential_window_pole(measured, tau)?,
            None => measured,
        };
        let natural_omega = physical.abs();
        if !natural_omega.is_finite() || natural_omega <= 0.0 {
            return Err(ModalIdError::Shape {
                what: "ringdown natural frequency must be positive and representable",
            }
            .into());
        }
        let mut shape = Vec::with_capacity(channels.len());
        for channel in &fit.fits {
            shape.push(channel.phasors[source].conj().scale(0.5));
        }
        ordered.push((
            Mode {
                frequency_hz: natural_omega / TAU,
                damping_ratio: -physical.re / natural_omega,
                shape,
                frequency_ci_hz: f64::NAN,
                damping_ci: f64::NAN,
            },
            source,
        ));
    }
    // An in-place unstable sort uses no allocation. Source index breaks exact
    // ties deterministically without claiming repeated-mode identifiability.
    ordered.sort_unstable_by(|a, b| {
        a.0.frequency_hz
            .total_cmp(&b.0.frequency_hz)
            .then_with(|| a.1.cmp(&b.1))
    });
    let mut modes = Vec::with_capacity(ordered.len());
    let mut source_pole_indices = Vec::with_capacity(ordered.len());
    for (mode, source) in ordered {
        modes.push(mode);
        source_pole_indices.push(source);
    }
    Ok(RingdownIdentification {
        modes,
        source_pole_indices,
        fit,
    })
}
