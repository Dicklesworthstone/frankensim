//! Project actual activation/decay events through the original evaluated response.

use std::sync::Arc;

use super::{BurnupDecayRadiation, DecayRadiationStep, invalid, visit};
use crate::burnup::stochastic::StochasticActivationStep;
use crate::burnup::{BurnupModel, BurnupRate};
use crate::decay::stochastic::MAX_STOCHASTIC_NUCLEI;
use crate::endf::EndfError;
use crate::endf::spectrum::PositiveIntegral;
use crate::{FissionError, nonnegative};

impl BurnupDecayRadiation {
    /// Expected radioactive spectral energy conditional on actual integer events.
    ///
    /// Capture remains a separate prompt-energy account. Its outward radiation
    /// is not assigned a daughter-decay spectrum. No random stream is consumed,
    /// no trajectory is evolved again, and no photon multiplicity is sampled.
    /// Neutron-producing events are refused because this stochastic report does
    /// not represent their distinct event and prompt-energy accounts.
    pub fn from_stochastic_activation_step(
        &self,
        model: &Arc<BurnupModel>,
        step: &StochasticActivationStep,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        checkpoint()?;
        if !Arc::ptr_eq(model, &self.model) {
            return Err(invalid(
                "radiation report belongs to a different stochastic activation model",
            ));
        }
        if !nonnegative(step.start_time_s)
            || !nonnegative(step.end_time_s)
            || step.end_time_s < step.start_time_s
            || step.duration_s != step.end_time_s - step.start_time_s
            || step.reaction_counts.len() != model.reactions().len()
            || (step.duration_s == 0.0 && step.events != 0)
        {
            return Err(invalid(
                "matching stochastic activation interval and original reaction counts required",
            ));
        }
        if step.events > MAX_STOCHASTIC_NUCLEI {
            return Err(FissionError::Budget("exact stochastic radiation event ceiling").into());
        }
        let mut work = 0;
        let mut events = 0_u64;
        let mut decays = 0_u64;
        let mut captures = 0_u64;
        let mut counts = vec![0_u64; model.nuclides().len()];
        let mut original = PositiveIntegral::new();
        for (reaction, &count) in model.reactions().iter().zip(&step.reaction_counts) {
            visit(&mut work, self.source.limits, &mut checkpoint)?;
            events = events
                .checked_add(count)
                .filter(|&total| total <= MAX_STOCHASTIC_NUCLEI)
                .ok_or(FissionError::Budget(
                    "exact stochastic radiation event ceiling",
                ))?;
            match reaction.rate {
                BurnupRate::Decay { .. } => {
                    // Each subtotal is bounded by the checked whole-report sum.
                    decays += count;
                    counts[reaction.parent] += count;
                    original.add_product(
                        1.0,
                        count as f64,
                        reaction.energy_per_event.radiation_j,
                        1.0,
                    )?;
                }
                BurnupRate::Capture { .. } | BurnupRate::TabulatedCapture { .. } => {
                    captures += count;
                }
                BurnupRate::NeutronProduction { .. }
                | BurnupRate::TabulatedNeutronProduction { .. } => {
                    if count != 0 {
                        return Err(invalid(
                            "stochastic activation report cannot contain neutron-producing events",
                        ));
                    }
                }
                BurnupRate::Fission { .. } | BurnupRate::TabulatedFission { .. } => {
                    if count != 0 {
                        return Err(invalid(
                            "stochastic activation report cannot contain fission events",
                        ));
                    }
                }
            }
        }
        if (events, decays, captures) != (step.events, step.decays, step.captures) {
            return Err(invalid(
                "stochastic activation totals and reaction counts disagree",
            ));
        }
        let counts: Vec<_> = counts.into_iter().map(|count| count as f64).collect();
        self.source.apply(
            step.duration_s,
            &counts,
            true,
            Some(original.value()?),
            &mut work,
            &mut checkpoint,
        )
    }
}
