//! Reconcile evaluated spectral shape with the original radioactive energy account.
//!
//! The selected spectrum must cover all outward energy of every emitting parent
//! within an explicit caller-supplied discrepancy allowance. Another outward
//! decay-heat category is never silently replaced by this spectrum. Accepted
//! spectral weights are scaled to the original parent-average radiation energy;
//! the adjustment is reported separately from numerical error and creates no
//! additional heat. This does not establish nuclear-data uncertainty or transport
//! accuracy, infer a missing spectrum, or replace the original decay/burnup owner.

use std::sync::Arc;

use crate::burnup::{BurnupModel, BurnupRate, BurnupStep};
use crate::decay::stochastic::{MAX_STOCHASTIC_NUCLEI, StochasticDecayStep};
use crate::decay::{DecayModel, DecayStep};
use crate::endf::EndfError;
use crate::endf::decay::model::emissions::{
    DecayEmissionCoverage, DecayEmissionLimits, PreparedDecayEmission,
};
use crate::endf::decay::spectra::DecayRadiation;
use crate::endf::spectrum::PositiveIntegral;
use crate::{Accumulator, FissionError, nonnegative};

mod stochastic_activation;

/// Explicit spectral-data agreement allowance, not a numerical-error tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayRadiationClosure {
    /// Absolute allowance for each parent's outward spectral energy [J/decay].
    pub absolute_j_per_decay: f64,
    /// Relative allowance against the original outward energy of that parent.
    pub relative: f64,
}

/// Integrated radiation partition; all energies are outward source joules.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayRadiationStep {
    pub duration_s: f64,
    pub radiation: DecayRadiation,
    pub group_edges_j: Arc<[f64]>,
    pub emitted_energy_j: Vec<f64>,
    /// Original radioactive radiation allocated here. Other burnup reactions
    /// retain their own independent energy accounts.
    pub original_radiation_j: f64,
    /// Sum of absolute spectral-to-aggregate adjustments over actual decay
    /// histories. This is an observable data discrepancy, not additional energy
    /// generation, an uncertainty bound or a numerical-defect allowance.
    pub reconciled_spectral_energy_j: f64,
}

/// Reusable spectral partition sharing the original evaluation and decay model.
#[derive(Clone, Debug)]
pub struct PreparedDecayRadiation {
    source: Arc<PreparedDecayEmission>,
    edges: Arc<[f64]>,
    energy: Arc<[f64]>,
    original: Arc<[f64]>,
    discrepancy: Arc<[f64]>,
    branch_sums: Arc<[f64]>,
    limits: DecayEmissionLimits,
}

impl PreparedDecayRadiation {
    pub fn new(
        source: Arc<PreparedDecayEmission>,
        closure: DecayRadiationClosure,
        limits: DecayEmissionLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        if !nonnegative(closure.absolute_j_per_decay) || !nonnegative(closure.relative) {
            return Err(invalid(
                "finite nonnegative spectral closure allowances required",
            ));
        }
        let category = category(source.radiation())?;
        let assembly = source.assembly();
        let parents = assembly.nuclides().len();
        let groups = source.group_edges_j().len() - 1;
        if groups > limits.max_groups
            || parents
                .checked_mul(groups + 3)
                .is_none_or(|n| n > limits.max_group_values)
        {
            return Err(FissionError::Budget("decay radiation coefficient storage").into());
        }
        let mut original = vec![PositiveIntegral::new(); parents];
        let mut branch_sums = vec![Accumulator::new(0.0); parents];
        let mut work = 0;
        for branch in assembly.model().branches() {
            visit(&mut work, limits, &mut checkpoint)?;
            original[branch.parent].add_product(
                1.0,
                branch.probability,
                branch.energy_per_decay.radiation_j,
                1.0,
            )?;
            branch_sums[branch.parent].add(branch.probability)?;
        }
        let original: Vec<_> = original
            .into_iter()
            .map(PositiveIntegral::value)
            .collect::<Result<_, _>>()?;
        let mut energy = Vec::with_capacity(parents * groups);
        let mut discrepancy = Vec::with_capacity(parents);
        for (parent, &target) in original.iter().enumerate() {
            visit(&mut work, limits, &mut checkpoint)?;
            let routing = assembly.routing()[parent];
            let fractions = [
                routing.light_particles.radiation,
                routing.electromagnetic.radiation,
                routing.heavy_particles.radiation,
            ];
            if target > 0.0 {
                for (other, &fraction) in fractions.iter().enumerate() {
                    let average = assembly.evaluations()[parent].average_energies()[other]
                        .value_j
                        .ok_or_else(|| invalid("original outward energy category is unknown"))?;
                    if other != category && average > 0.0 && fraction > 0.0 {
                        return Err(invalid(
                            "selected spectrum does not cover another outward energy category",
                        ));
                    }
                }
                if source.coverage()[parent] != DecayEmissionCoverage::Reported {
                    return Err(invalid(
                        "an emitting parent lacks the selected radiation spectrum",
                    ));
                }
            }
            let raw = source
                .energy_per_decay_j(parent)
                .expect("original parent order");
            let mut spectral = PositiveIntegral::new();
            for &coefficient in raw {
                visit(&mut work, limits, &mut checkpoint)?;
                if target > 0.0 {
                    spectral.add_product(1.0, coefficient, fractions[category], 1.0)?;
                }
            }
            let spectral_j = spectral.value()?;
            let difference = (spectral_j - target).abs();
            if target > 0.0 {
                if spectral_j == 0.0 {
                    return Err(invalid(
                        "positive original radiation requires positive selected spectral energy",
                    ));
                }
                let mut allowance = PositiveIntegral::new();
                allowance.add_product(1.0, closure.absolute_j_per_decay, 1.0, 1.0)?;
                allowance.add_product(1.0, closure.relative, target, 1.0)?;
                if difference > allowance.value()? {
                    return Err(invalid(
                        "selected spectral energy disagrees with original radiation beyond closure allowance",
                    ));
                }
            }
            discrepancy.push(difference);
            for &coefficient in raw {
                visit(&mut work, limits, &mut checkpoint)?;
                let mut reconciled = PositiveIntegral::new();
                let value = if target == 0.0 {
                    0.0
                } else {
                    reconciled.add_triple_product(
                        1.0,
                        coefficient,
                        fractions[category],
                        target,
                        1.0,
                    )?;
                    reconciled.ratio(spectral)?
                };
                if coefficient.is_subnormal() && value > coefficient {
                    return Err(FissionError::Arithmetic.into());
                }
                energy.push(value);
            }
        }
        checkpoint()?;
        Ok(Self {
            edges: source.group_edges_j().into(),
            source,
            energy: energy.into(),
            original: original.into(),
            discrepancy: discrepancy.into(),
            branch_sums: branch_sums
                .into_iter()
                .map(|sum| sum.value())
                .collect::<Vec<_>>()
                .into(),
            limits,
        })
    }

    pub fn source(&self) -> &Arc<PreparedDecayEmission> {
        &self.source
    }
    pub fn model(&self) -> &Arc<DecayModel> {
        self.source.assembly().model()
    }
    pub fn radiation(&self) -> DecayRadiation {
        self.source.radiation()
    }
    pub fn group_edges_j(&self) -> &[f64] {
        &self.edges
    }
    pub fn discrepancy_per_decay_j(&self) -> &[f64] {
        &self.discrepancy
    }
    pub fn original_radiation_per_decay_j(&self) -> &[f64] {
        &self.original
    }
    pub fn energy_per_decay_j(&self, parent: usize) -> Option<&[f64]> {
        let groups = self.edges.len() - 1;
        (parent < self.original.len()).then(|| &self.energy[parent * groups..(parent + 1) * groups])
    }

    /// The caller supplies the producing inventory's original model Arc. A
    /// report is an ordinary caller-supplied value, not an authenticated receipt.
    pub fn from_decay_step(
        &self,
        model: &Arc<DecayModel>,
        step: &DecayStep,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        checkpoint()?;
        if !Arc::ptr_eq(model, self.model()) {
            return Err(invalid(
                "radiation report belongs to a different decay model",
            ));
        }
        self.apply(
            step.duration_s,
            &step.decays,
            false,
            None,
            &mut 0,
            &mut checkpoint,
        )
    }

    /// Expected spectral energy conditional on an actual integer decay history.
    /// Original energy uses integer branch counts and the model's energies per
    /// event; this does not sample photons, angular correlations or decay
    /// cascades. It consumes no randomness and does not advance a trajectory.
    /// Original counts must agree exactly and stay within the stochastic owner's
    /// 2^53 event ceiling before conversion to binary64.
    pub fn from_stochastic_step(
        &self,
        model: &Arc<DecayModel>,
        step: &StochasticDecayStep,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        checkpoint()?;
        if !Arc::ptr_eq(model, self.model()) {
            return Err(invalid(
                "radiation report belongs to a different stochastic decay model",
            ));
        }
        if !nonnegative(step.start_time_s)
            || !nonnegative(step.end_time_s)
            || step.end_time_s < step.start_time_s
            || step.duration_s != step.end_time_s - step.start_time_s
            || step.decays.len() != self.original.len()
            || step.branch_decays.len() != model.branches().len()
        {
            return Err(invalid(
                "matching stochastic interval and original-order count shapes required",
            ));
        }
        if step.events > MAX_STOCHASTIC_NUCLEI {
            return Err(
                FissionError::Budget("exact stochastic radiation event-count ceiling").into(),
            );
        }
        let mut work = 0;
        let mut counted = vec![0_u64; self.original.len()];
        let mut events = 0_u64;
        let mut original = PositiveIntegral::new();
        for (branch, &count) in model.branches().iter().zip(&step.branch_decays) {
            visit(&mut work, self.limits, &mut checkpoint)?;
            if count > MAX_STOCHASTIC_NUCLEI {
                return Err(
                    FissionError::Budget("exact stochastic radiation event-count ceiling").into(),
                );
            }
            events = events
                .checked_add(count)
                .filter(|&sum| sum <= MAX_STOCHASTIC_NUCLEI)
                .ok_or(FissionError::Budget(
                    "exact stochastic radiation event-count ceiling",
                ))?;
            counted[branch.parent] += count; // bounded by the checked global sum
            original.add_product(1.0, count as f64, branch.energy_per_decay.radiation_j, 1.0)?;
        }
        if events != step.events {
            return Err(invalid("stochastic total and branch decay counts disagree"));
        }
        let mut decays = Vec::with_capacity(counted.len());
        for (&actual, &reported) in counted.iter().zip(&step.decays) {
            visit(&mut work, self.limits, &mut checkpoint)?;
            if actual != reported {
                return Err(invalid(
                    "stochastic nuclide and branch decay counts disagree",
                ));
            }
            decays.push(actual as f64);
        }
        self.apply(
            step.duration_s,
            &decays,
            true,
            Some(original.value()?),
            &mut work,
            &mut checkpoint,
        )
    }

    /// Validate the original evaluated decay channels once. Extra capture and
    /// fission channels are allowed and retain their own radiation accounts.
    pub fn bind_burnup(
        &self,
        model: Arc<BurnupModel>,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupDecayRadiation, EndfError> {
        checkpoint()?;
        let mut work = self
            .model()
            .branches()
            .len()
            .checked_add(model.reactions().len())
            .ok_or(FissionError::Budget("decay radiation binding work"))?;
        if work > self.limits.max_work {
            return Err(FissionError::Budget("decay radiation binding work").into());
        }
        self.source
            .bind_burnup(Arc::clone(&model), &mut checkpoint)?;
        let mut reactions = Vec::new();
        for (index, reaction) in model.reactions().iter().enumerate() {
            visit(&mut work, self.limits, &mut checkpoint)?;
            if matches!(reaction.rate, BurnupRate::Decay { .. }) {
                reactions.push((
                    index,
                    reaction.parent,
                    reaction.energy_per_event.radiation_j,
                ));
            }
        }
        checkpoint()?;
        Ok(BurnupDecayRadiation {
            source: self.clone(),
            model,
            reactions: reactions.into(),
        })
    }

    fn apply(
        &self,
        duration_s: f64,
        decays: &[f64],
        branch_counts: bool,
        authoritative_j: Option<f64>,
        work: &mut usize,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        if !nonnegative(duration_s)
            || decays.len() != self.original.len()
            || !decays.iter().copied().all(nonnegative)
            || (duration_s == 0.0 && decays.iter().any(|&count| count != 0.0))
        {
            return Err(invalid(
                "matching finite decay counts and duration required",
            ));
        }
        let groups = self.edges.len() - 1;
        let mut energy = vec![PositiveIntegral::new(); groups];
        let mut original = PositiveIntegral::new();
        let mut discrepancy = PositiveIntegral::new();
        for (parent, &count) in decays.iter().enumerate() {
            visit(work, self.limits, checkpoint)?;
            if count == 0.0 {
                continue;
            }
            if self.branch_sums[parent] == 0.0 {
                return Err(invalid("a stable nuclide cannot have decay events"));
            }
            // Burnup reports count individual branches. Their admitted sum is
            // near one but is not silently assumed to equal one exactly.
            let divisor = if branch_counts {
                self.branch_sums[parent]
            } else {
                1.0
            };
            if count > divisor && self.original[parent].is_subnormal() {
                return Err(FissionError::Arithmetic.into());
            }
            original.add_product(1.0, count, self.original[parent], divisor)?;
            discrepancy.add_product(1.0, count, self.discrepancy[parent], divisor)?;
            for (group, &coefficient) in energy
                .iter_mut()
                .zip(self.energy_per_decay_j(parent).unwrap())
            {
                visit(work, self.limits, checkpoint)?;
                if count > divisor && coefficient.is_subnormal() {
                    return Err(FissionError::Arithmetic.into());
                }
                group.add_product(1.0, count, coefficient, divisor)?;
            }
        }
        let result = DecayRadiationStep {
            duration_s,
            radiation: self.radiation(),
            group_edges_j: Arc::clone(&self.edges),
            emitted_energy_j: energy
                .into_iter()
                .map(PositiveIntegral::value)
                .collect::<Result<_, _>>()?,
            original_radiation_j: match authoritative_j {
                Some(value) => value,
                None => original.value()?,
            },
            reconciled_spectral_energy_j: discrepancy.value()?,
        };
        checkpoint()?;
        Ok(result)
    }
}

/// Original grouped decay partition bound to a full reaction model.
#[derive(Clone, Debug)]
pub struct BurnupDecayRadiation {
    source: PreparedDecayRadiation,
    model: Arc<BurnupModel>,
    reactions: Arc<[(usize, usize, f64)]>,
}

impl BurnupDecayRadiation {
    pub fn source(&self) -> &PreparedDecayRadiation {
        &self.source
    }
    pub fn model(&self) -> &Arc<BurnupModel> {
        &self.model
    }

    pub fn from_burnup_step(
        &self,
        model: &Arc<BurnupModel>,
        step: &BurnupStep,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        checkpoint()?;
        if !Arc::ptr_eq(model, &self.model) {
            return Err(invalid(
                "radiation report belongs to a different burnup model",
            ));
        }
        if !nonnegative(step.duration_s)
            || step.reaction_counts.len() != self.model.reactions().len()
            || !step.reaction_counts.iter().copied().all(nonnegative)
            || (step.duration_s == 0.0 && step.reaction_counts.iter().any(|&count| count != 0.0))
        {
            return Err(invalid(
                "matching finite original-order reaction counts and duration required",
            ));
        }
        let mut counts = vec![Accumulator::new(0.0); self.model.nuclides().len()];
        let mut original = PositiveIntegral::new();
        let mut work = 0;
        for &(reaction, parent, energy) in self.reactions.iter() {
            visit(&mut work, self.source.limits, &mut checkpoint)?;
            counts[parent].add(step.reaction_counts[reaction])?;
            original.add_product(1.0, step.reaction_counts[reaction], energy, 1.0)?;
        }
        let decays: Vec<_> = counts.into_iter().map(|sum| sum.value()).collect();
        self.source.apply(
            step.duration_s,
            &decays,
            true,
            Some(original.value()?),
            &mut work,
            &mut checkpoint,
        )
    }
}

fn category(radiation: DecayRadiation) -> Result<usize, EndfError> {
    match radiation {
        DecayRadiation::Gamma | DecayRadiation::XRay => Ok(1),
        DecayRadiation::BetaMinus
        | DecayRadiation::ElectronCapturePositron
        | DecayRadiation::Electron => Ok(0),
        DecayRadiation::Alpha | DecayRadiation::Neutron | DecayRadiation::Proton => Ok(2),
        _ => Err(EndfError::Unsupported(
            "unresolved outward decay radiation category",
        )),
    }
}

fn visit(
    work: &mut usize,
    limits: DecayEmissionLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    if *work >= limits.max_work {
        return Err(FissionError::Budget("decay radiation work").into());
    }
    *work += 1;
    if *work % 64 == 0 {
        checkpoint()?;
    }
    Ok(())
}

fn invalid(message: &'static str) -> EndfError {
    FissionError::InvalidInput(message).into()
}

#[cfg(test)]
pub(crate) mod tests;
