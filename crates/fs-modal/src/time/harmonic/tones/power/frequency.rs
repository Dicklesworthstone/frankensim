//! A source-clock variation is an exact affine envelope on the original tone.

use super::gradient::{
    ModalTonePowerDerivative, ModalTonePowerGradientBudget, admit_layout, finite, poll,
};
use super::{ModalAcousticState, ModalAcousticTimeError, TonePower, invalid};
use crate::time::harmonic::tones::{
    MAX_HARMONIC_TONES, ModalToneBudget, PreparedModalTones, admit_cost,
};

/// Complete coherent work/heat derivatives for independently authored source
/// angular frequencies, at fixed interval-start quadratures and mechanics.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalTonePowerFrequencyGradient {
    sources: [ModalTonePowerDerivative; MAX_HARMONIC_TONES],
    count: usize,
}

impl ModalTonePowerFrequencyGradient {
    /// Original input order, including coincident frequencies. Both entries in
    /// each pair have units J*s/rad. A shared fundamental Omega_j=n_j*Omega_0
    /// requires the caller's sum of n_j times these independent derivatives.
    pub fn source_derivatives(&self) -> &[ModalTonePowerDerivative] {
        &self.sources[..self.count]
    }
}

/// Original constant-tone bank plus exact affine power companion factors.
///
/// dF_j/dOmega_j = t*(S_j*cos(Omega_j*t)-C_j*sin(Omega_j*t)). This is the
/// original affine forcing law with zero start and h*(S_j,-C_j) end variation.
/// Its full power transpose includes both the explicit source-work term and
/// the changed coherent velocity, retaining all cross-frequency contributions.
/// No source-frequency perturbation or separately differentiated solver is used.
///
/// Warm queries allocate nothing, call no transcendentals and poll cancellation
/// between factor rows and original source contractions. Duration, mechanics,
/// source initial phase and mass/basis normalization stay fixed. An affine BASE
/// would require quadratic-envelope power factors and is explicitly refused.
#[derive(Debug)]
pub struct PreparedModalTonePowerFrequencyGradient<'a> {
    bank: &'a PreparedModalTones,
    companions: Box<[TonePower]>,
}

impl PreparedModalTones {
    /// Prepare one affine power companion per original mode. The original bank,
    /// motion coefficients, accepted state and model binding remain unchanged.
    ///
    /// The existing affine bank admission is a conservative cold byte/work
    /// bound for these companions; no duplicate motion or model is constructed.
    /// It excludes the borrowed original bank. Query admission separately uses
    /// the companion's full affine dimension and the original source count.
    /// Cancellation is polled before allocation and between bounded companion
    /// preparations. Refusal publishes no partial cache or physical state.
    pub fn prepare_power_frequency_gradients(
        &self,
        preparation_budget: ModalToneBudget,
        query_budget: ModalTonePowerGradientBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedModalTonePowerFrequencyGradient<'_>, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        if self.affine {
            return Err(invalid(
                "source-frequency power gradients require constant-envelope tones",
            ));
        }
        if self.tone_count() > preparation_budget.maximum_tones {
            return Err(ModalAcousticTimeError::BudgetExceeded {
                what: "power-frequency companion tone count",
                value: self.tone_count() as f64,
                limit: preparation_budget.maximum_tones as f64,
            });
        }
        admit_layout(self.unique.len(), self.tone_count(), true, query_budget)?;
        admit_cost(
            self.mode_count(),
            |index| {
                let motion = &self.modes[index].motion[0];
                (motion.angular_frequency_rad_s, motion.transition.drag)
            },
            &self.unique,
            self.tone_count(),
            self.duration_s,
            preparation_budget,
            true,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let mut companions = Vec::with_capacity(self.mode_count());
        for mode in &self.modes {
            poll(&mut cancelled)?;
            let motion = &mode.motion[0];
            companions.push(TonePower::try_new_affine(
                motion.angular_frequency_rad_s,
                motion.transition.drag,
                &self.unique,
                self.duration_s,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(PreparedModalTonePowerFrequencyGradient {
            bank: self,
            companions: companions.into_boxed_slice(),
        })
    }
}

impl PreparedModalTonePowerFrequencyGradient<'_> {
    /// Coherent power derivatives at fixed local initial state and source phase.
    /// Exact duplicates share nominal factors, while each source retains its
    /// own original quadratures in the final clock contraction. Detuning one
    /// cancelled source therefore remains a live physical direction.
    ///
    /// A multi-interval reverse sweep also adds the motion clock derivative and
    /// the source-phase advancement chain through pullback_force_phase. Obtain
    /// state/source amplitude covectors from prepare_power_gradients. Local
    /// source-frequency partials alone do not include earlier phase history.
    pub fn gradients_mode(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        cosine: &[f64],
        sine: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalTonePowerFrequencyGradient, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        let companion = self
            .companions
            .get(mode_index)
            .ok_or_else(|| invalid("tone mode index must be admitted"))?;
        let (c, s) = self.bank.collect_forces(cosine, sine)?;
        let k = self.bank.unique.len();
        let zero = [0.0; MAX_HARMONIC_TONES];
        let gradient = companion.input_gradients(
            state,
            &c[..k],
            &s[..k],
            Some((&zero[..k], &zero[..k])),
            &mut cancelled,
        )?;
        let mut sources = [ModalTonePowerDerivative::default(); MAX_HARMONIC_TONES];
        for (source, &unique) in self.bank.source_to_unique.iter().enumerate() {
            poll(&mut cancelled)?;
            let (c_end, s_end) = (gradient[4 + 4 * unique], gradient[5 + 4 * unique]);
            // Contract original source amplitudes AFTER nominal coalescing. Do
            // not multiply h*C first: an individually large virtual envelope can
            // overflow while the admitted physical derivative is finite or zero.
            sources[source] = finite(ModalTonePowerDerivative {
                input_work: self.bank.duration_s
                    * sine[source].mul_add(c_end.input_work, -cosine[source] * s_end.input_work),
                viscous_dissipation: self.bank.duration_s
                    * sine[source].mul_add(
                        c_end.viscous_dissipation,
                        -cosine[source] * s_end.viscous_dissipation,
                    ),
            })?;
        }
        poll(&mut cancelled)?;
        Ok(ModalTonePowerFrequencyGradient {
            sources,
            count: self.bank.tone_count(),
        })
    }
}
