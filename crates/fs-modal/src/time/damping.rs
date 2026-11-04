//! Authoritative physical viscous rates and legacy ratio admission.

use super::{
    ExactZohKernel, HarmonicForceTransition, ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeError, ModalAcousticTimeModel,
    advance_exact_zoh_cached,
};
use std::sync::Arc;

pub(super) fn rate_from_ratio(mode: ModalAcousticMode) -> f64 {
    // Preserve the established product for ordinary coefficients. Regroup only
    // when doubling a finite ratio overflows before a soft frequency rescales it.
    let rate = 2.0 * mode.damping_ratio * mode.angular_frequency_rad_s;
    if rate.is_finite() {
        rate
    } else {
        2.0 * (mode.damping_ratio * mode.angular_frequency_rad_s)
    }
}

pub(super) fn diagnostic_ratio(omega: f64, rate: f64) -> f64 {
    if rate == 0.0 {
        return 0.0;
    }
    let ratio = rate / (2.0 * omega);
    if ratio.is_finite() && ratio > 0.0 {
        ratio
    } else {
        // A deliberately nonfinite legacy diagnostic makes a missed consumer
        // refuse instead of silently dropping physical damping. The dynamics
        // and all power integrals use the separately admitted finite rate.
        f64::INFINITY
    }
}

pub(super) fn validate_rates(rates: &[f64], count: usize) -> Result<(), ModalAcousticTimeError> {
    if rates.len() != count || rates.iter().any(|rate| !rate.is_finite() || *rate < 0.0) {
        return Err(ModalAcousticTimeError::InvalidInput {
            what: "modal damping needs one finite nonnegative physical rate per coordinate",
        });
    }
    Ok(())
}

impl ModalAcousticTimeModel {
    /// Admit physical coefficients `c` [1/s] in `q''+c*q'+omega^2*q=F`.
    /// For dimensional viscous force `-b*v` on mass `m`, supply `c=b/m`.
    /// The supplied modes define frequencies and pressure transfers; their
    /// damping-ratio fields are replaced with diagnostics derived from `c`.
    /// Rates replace the complete damping law; they are not added to ratios.
    ///
    /// Zero stiffness is admitted with positive drag and zero pressure transfer.
    /// No small stiffness, damping-ratio conversion, or free-body acoustic
    /// transfer is inferred. All states start at rest, and the existing state,
    /// energy, observation and transactional balance budgets apply unchanged.
    ///
    /// [`Self::damping_rates_per_s`] is the authoritative physical metadata.
    /// A positive rate whose ratio is undefined or unrepresentable is marked
    /// by positive infinity in [`Self::modes`], never by a misleading zero.
    /// Internal transitions and power integrals do not use that diagnostic.
    ///
    /// # Errors
    /// Invalid clock, mode count, frequencies, pressure transfers, rates or
    /// budgets, or natural frequencies above the admitted Nyquist guard.
    pub fn try_new_with_damping_rates(
        sample_rate_hz: u32,
        modes: Vec<ModalAcousticMode>,
        damping_rates_per_s: &[f64],
        budget: ModalAcousticTimeBudget,
    ) -> Result<Self, ModalAcousticTimeError> {
        Self::try_new_impl(
            sample_rate_hz,
            modes,
            budget,
            true,
            Some(damping_rates_per_s),
        )
    }

    /// Physical viscous coefficients `c` [1/s], one per admitted coordinate.
    /// These are authoritative for every bank, including legacy ratio inputs,
    /// zero-stiffness coordinates and very soft modes. Use them for generator
    /// assembly, force/heat evaluation, parameter updates and cache identity;
    /// reconstructing `2*zeta*omega` from diagnostic ratios can be undefined.
    #[must_use]
    pub fn damping_rates_per_s(&self) -> &[f64] {
        &self.transition_kernel.damping_rates_per_s
    }

    /// Physical viscous coefficient of one retained coordinate [1/s].
    ///
    /// # Panics
    /// An index outside the admitted mode set, like indexing [`Self::modes`].
    #[must_use]
    pub fn damping_rate_per_s(&self, mode_index: usize) -> f64 {
        self.transition_kernel.damping_rates_per_s[mode_index]
    }

    /// Replace complete physical viscous rates in the same admitted coordinates.
    /// States, frequencies, observer transfers, clock and budgets are retained;
    /// no motion or stored energy changes. Changed rates obtain a new immutable
    /// kernel, so preparation tied to old mechanics cannot be reused silently.
    /// Identical rates share the existing kernel and its prepared caches.
    ///
    /// # Errors
    /// Missing, negative or nonfinite rates. The source model is never mutated.
    pub fn with_damping_rates(&self, rates: &[f64]) -> Result<Self, ModalAcousticTimeError> {
        validate_rates(rates, self.modes.len())?;
        let mut candidate = self.clone();
        if rates != self.damping_rates_per_s() {
            for (mode, &rate) in candidate.modes.iter_mut().zip(rates) {
                mode.damping_ratio = diagnostic_ratio(mode.angular_frequency_rad_s, rate);
            }
            candidate.transition_kernel = Arc::new(ExactZohKernel::with_rates(
                &candidate.modes,
                self.sample_period_s,
                rates,
            ));
        }
        Ok(candidate)
    }

    /// Preview one coordinate at any positive finite duration using its admitted
    /// physical damping law. Full/half nominal intervals reuse prepared motion;
    /// other durations prepare the same exact transition. This pure operation
    /// applies no state/energy budget and never advances the accepted model.
    ///
    /// # Errors
    /// Invalid mode index, duration or inputs, or unrepresentable endpoint motion.
    pub fn advance_exact_zoh_duration_step(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        force: f64,
        duration_s: f64,
    ) -> Result<ModalAcousticState, ModalAcousticTimeError> {
        let Some(&mode) = self.modes.get(mode_index) else {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "modal preview index is out of range",
            });
        };
        if duration_s <= 0.0
            || ![
                duration_s,
                state.displacement_m_sqrt_kg,
                state.velocity_m_sqrt_kg_per_s,
                force,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "modal preview requires positive finite duration and finite state/force",
            });
        }
        let transition = match self.transition_kernel.for_duration(&self.modes, duration_s) {
            Some(prepared) => prepared[mode_index].transition,
            None => self
                .transition_kernel
                .transition(mode_index, mode, duration_s),
        };
        let next = advance_exact_zoh_cached(mode, state, force, transition);
        if !next.displacement_m_sqrt_kg.is_finite() || !next.velocity_m_sqrt_kg_per_s.is_finite() {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "modal preview motion must be representable",
            });
        }
        Ok(next)
    }

    pub(super) fn damping_over_interval(&self, mode_index: usize, duration_s: f64) -> f64 {
        if self.transition_kernel.explicit_rates {
            self.damping_rate_per_s(mode_index) * duration_s
        } else {
            let mode = self.modes[mode_index];
            // Retain the legacy ratio path's original floating-point product.
            2.0 * (mode.damping_ratio * (mode.angular_frequency_rad_s * duration_s))
        }
    }

    pub(super) fn prepare_harmonic_mode(
        &self,
        mode_index: usize,
        forcing_angular_frequency_rad_s: f64,
        duration_s: f64,
    ) -> Result<HarmonicForceTransition, ModalAcousticTimeError> {
        let mode = self.modes[mode_index];
        if self.transition_kernel.explicit_rates {
            HarmonicForceTransition::try_new_with_damping_rate(
                mode.angular_frequency_rad_s,
                self.damping_rate_per_s(mode_index),
                forcing_angular_frequency_rad_s,
                duration_s,
            )
        } else {
            HarmonicForceTransition::try_new(
                mode.angular_frequency_rad_s,
                mode.damping_ratio,
                forcing_angular_frequency_rad_s,
                duration_s,
            )
        }
    }
}
