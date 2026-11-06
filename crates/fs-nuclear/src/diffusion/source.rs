//! Uniform external neutron and precursor production during a spatial step.
//!
//! Source vectors are extensive rates [expected particles/s], in the same
//! cell-major order as the inventories. They enter continuously throughout
//! the interval, rather than being added at its beginning or endpoint. The
//! original reaction source owner retains the age-dependent population,
//! removal and residence responses, including populations in zero-loss nodes.
//! A source history that changes in time requires explicitly split intervals.

use fs_la::decay::reaction::{ReactionAdvance, ReactionPlan, source::ReactionSourcePlan};

use super::{
    Accounts, DiffusionBudget, DiffusionError, DiffusionLimits, DiffusionModel, DiffusionState,
    DiffusionStep, FissionError, SOURCE_NEUTRONS, SOURCE_PRECURSORS, filled, nonnegative,
    owner_call,
};

impl DiffusionState {
    /// Independent candidate with uniform external injection over the interval.
    /// Both vectors must have exactly the corresponding inventory lengths,
    /// including an empty precursor vector when no families are admitted.
    pub fn trial_with_sources(
        &self,
        duration_s: f64,
        neutron_rates_s: &[f64],
        precursor_rates_s: &[f64],
        budget: DiffusionBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, DiffusionStep), DiffusionError> {
        self.trial_from_sources(
            duration_s,
            Some(SourceRates {
                neutrons: neutron_rates_s,
                precursors: precursor_rates_s,
            }),
            budget,
            checkpoint,
        )
    }

    /// Publish sources, transport and reactions only after every account passes.
    pub fn advance_with_sources(
        &mut self,
        duration_s: f64,
        neutron_rates_s: &[f64],
        precursor_rates_s: &[f64],
        budget: DiffusionBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DiffusionStep, DiffusionError> {
        let (candidate, step) = self.trial_with_sources(
            duration_s,
            neutron_rates_s,
            precursor_rates_s,
            budget,
            checkpoint,
        )?;
        *self = candidate;
        Ok(step)
    }
}

pub(super) struct SourceRates<'a> {
    pub neutrons: &'a [f64],
    pub precursors: &'a [f64],
}

pub(super) fn injected_counts(
    source: Option<SourceRates<'_>>,
    neutron_nodes: usize,
    nodes: usize,
    duration_s: f64,
    accounts: &mut Accounts,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<Vec<f64>, DiffusionError> {
    let mut injected = filled(nodes, 0.0)?;
    if let Some(source) = source {
        if source.neutrons.len() != neutron_nodes
            || source.precursors.len() != nodes - neutron_nodes
            || !source
                .neutrons
                .iter()
                .chain(source.precursors)
                .copied()
                .all(nonnegative)
        {
            return Err(DiffusionError::Invalid(
                "one finite nonnegative external rate per neutron and precursor node required",
            ));
        }
        for (node, &rate) in source.neutrons.iter().chain(source.precursors).enumerate() {
            checkpoint()?;
            let count = accounts.multiply(rate, duration_s)?;
            injected[node] = count;
            accounts.total(
                if node < neutron_nodes {
                    SOURCE_NEUTRONS
                } else {
                    SOURCE_PRECURSORS
                },
                count,
            )?;
        }
    }
    Ok(injected)
}

/// The source-capable map can also serve a subsequent zero-source interval.
/// Injection magnitudes never invalidate preparation for the same operator/dt.
#[derive(Clone, Debug)]
pub(super) enum PreparedPlan {
    Closed(ReactionPlan),
    Source(ReactionSourcePlan),
}

impl PreparedPlan {
    pub(super) fn prepare(
        model: &DiffusionModel,
        duration_s: f64,
        limits: DiffusionLimits,
        with_sources: bool,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DiffusionError> {
        if with_sources {
            Ok(Self::Source(owner_call(checkpoint, |cancelled| {
                ReactionSourcePlan::prepare_with_population_integrals(
                    &model.removal_rates,
                    &model.edges,
                    duration_s,
                    limits.propagation,
                    cancelled,
                )
            })?))
        } else {
            Ok(Self::Closed(owner_call(checkpoint, |cancelled| {
                ReactionPlan::prepare_with_population_integrals(
                    &model.removal_rates,
                    &model.edges,
                    duration_s,
                    limits.propagation,
                    cancelled,
                )
            })?))
        }
    }

    pub(super) fn duration_s(&self) -> f64 {
        match self {
            Self::Closed(plan) => plan.duration_s(),
            Self::Source(plan) => plan.homogeneous_plan().duration_s(),
        }
    }

    pub(super) fn supports_sources(&self) -> bool {
        matches!(self, Self::Source(_))
    }

    pub(super) fn apply(
        &self,
        initial: &[f64],
        injected: &[f64],
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<ReactionAdvance, DiffusionError> {
        owner_call(checkpoint, |cancelled| match self {
            Self::Closed(plan) => plan.apply_with_cancel(initial, cancelled),
            Self::Source(plan) => plan.apply_with_cancel(initial, injected, cancelled),
        })
    }

    pub(super) fn integrals(
        &self,
        initial: &[f64],
        injected: &[f64],
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<f64>, DiffusionError> {
        owner_call(checkpoint, |cancelled| match self {
            Self::Closed(plan) => plan.integrate_populations_with_cancel(initial, cancelled),
            Self::Source(plan) => {
                plan.integrate_populations_with_cancel(initial, injected, cancelled)
            }
        })
    }
}

#[cfg(test)]
mod tests;
