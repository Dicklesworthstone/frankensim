//! Finite fusion-neutron activation with separately owned energy accounts.
//!
//! The original spatial burnup owner advances neutron populations, target
//! depletion, activation products and their radioactive decay together. Kinetic
//! accounting consumes that SAME accepted neutron step; there is no second
//! neutron solve. Neutron kinetic transfers, explicit signed prompt Q and
//! separately authored inventory decay energy retain independent accounts.
//!
//! Every `BurnupReaction::energy_per_event` supplied here must exclude incident
//! neutron kinetic energy. A reaction with Q in the neutron energy model must
//! have zero prompt energy in burnup, as must every neutron-producing reaction,
//! including prompt fission;
//! double ownership refuses. Aggregate induced-fission precursors additionally
//! require an explicit unrepresented-product boundary. Their separately authored
//! future-neutron conversion becomes stored precursor energy at birth; prompt Q
//! remains unchanged. Positive retained fission products refuse that boundary.
//! Capture entries also exclude subsequent product
//! decay energy; the product's actual decay events release that energy later. No ENDF
//! Q, product yield, charged-particle stopping or gamma deposition is inferred.
//! Unbound neutron reactions remain explicitly fixed background interactions;
//! their energy declaration creates no tracked activation or target depletion.
//! When finite radioactive emitters are bound, their emitted neutron energy is
//! an internal transfer debited only from original decay radiation. Capture
//! radiation cannot fund it. Total kinetic injection remains diagnostic;
//! external injection and internal births are reported separately.
//!
//! The original spatial owner's interval-mean neutron flux and target population
//! approximation remains in force. Reaction-count convergence does not resolve
//! their within-step time correlation; refine physical intervals when needed.

use std::sync::Arc;

use super::{
    FusionNeutronSource, NeutronKineticError, NeutronKineticModel, NeutronKineticState,
    NeutronKineticStep, checked,
};
use crate::burnup::{BurnupInventory, BurnupRate};
use crate::diffusion::burnup::{
    DiffusionBurnupBudget, DiffusionBurnupError, DiffusionBurnupSources, DiffusionBurnupState,
    DiffusionBurnupStep,
};
use crate::diffusion::{DiffusionModel, DiffusionState};
use crate::{Accumulator, EnergyChannels, FissionError, nonnegative, sum};

mod decay;
use decay::DecayPartitions;

/// Original coupled particle/nuclide budgets and one combined energy allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FusionNeutronActivationBudget {
    /// Unchanged native diffusion, burnup and reaction-consistency limits.
    pub spatial: DiffusionBurnupBudget,
    /// Absolute sum of kinetic, nuclear, export and reaction-coupling defects
    /// [J/step]. This is a numerical allowance, not a physical uncertainty.
    pub energy_tolerance_j: f64,
}

/// One accepted neutron and finite-material history with explicit source energy.
#[derive(Clone, Debug)]
pub struct FusionNeutronActivationState {
    spatial: DiffusionBurnupState,
    // The native neutron checkpoint is cloned, not advanced independently. Its
    // shared model and arrays are exactly the accepted spatial checkpoint.
    kinetic: NeutronKineticState,
    bound_reactions: Arc<[bool]>,
    maximum_reaction_energy_j: Arc<[f64]>,
    decay_partitions: Option<Arc<DecayPartitions>>,
    nuclear_released: Vec<[Accumulator; 3]>,
    released: Vec<[Accumulator; 3]>,
}

/// Native particle/material results and disjoint energy exports for one interval.
#[derive(Clone, Debug)]
pub struct FusionNeutronActivationStep {
    /// Original finite-target/product report, including reaction-count defects.
    pub spatial: DiffusionBurnupStep,
    /// Incident/outgoing kinetic energy and signed prompt Q of exactly
    /// `spatial.neutrons`.
    /// Its neutron report is a retained copy, not another neutron solve.
    pub kinetic: NeutronKineticStep,
    /// Sum of the original per-cell burnup energy channels [J]. Excludes both
    /// incident kinetic energy and prompt Q already owned by `kinetic`.
    pub nuclear_released: EnergyChannels,
    /// Actual internal births from the original finite radioactive stock.
    pub decay_neutrons: f64,
    /// Internally born delayed families, in the same finite decay history.
    pub decay_precursors: f64,
    /// Original decay radiation transferred into neutron kinetic storage [J].
    /// Already removed from combined exports; do not add it as material heat.
    pub transferred_decay_neutron_energy_j: f64,
    /// Original decay radiation reserved in newly born delayed families [J].
    pub reserved_decay_precursor_energy_j: f64,
    /// Accepted external injection, excluding internally emitted neutrons.
    pub external_injected_neutrons: f64,
    /// External family births, excluding internally reserved decay cohorts.
    pub external_injected_precursors: f64,
    /// Kinetic energy supplied from outside this nuclear system [J].
    pub external_injected_kinetic_energy_j: f64,
    /// Stored future-neutron energy supplied by external precursors [J].
    pub external_injected_precursor_energy_j: f64,
    /// Noncancelling external/internal injection reconciliation [neutrons].
    pub injection_number_roundoff: f64,
    /// Net signed nuclear conversion: original inventory event releases plus
    /// kinetic prompt Q and explicitly funded induced-precursor storage [J].
    /// Endothermic reactions may make this negative. Subsequent represented
    /// product-decay energy enters only when its original decay events occur.
    pub nuclear_generated_energy_j: f64,
    /// Neutron-associated plus inventory exports in diffusion-cell order [J].
    pub cell_released: Vec<EnergyChannels>,
    /// Disjoint combined exports [J]. Only deposited_j is local material heat.
    pub released: EnergyChannels,
    /// Conservative energy weighting of native neutron/material count defects:
    /// each cell's defect times its largest absolute event-energy coefficient.
    pub reaction_energy_defect_j: f64,
    /// Final neutron and precursor storage - initial + all exports - external energy
    /// - generated nuclear-event energy [J].
    pub energy_balance_residual_j: f64,
    /// Noncancelling energy-account and coupling defect total [J].
    pub energy_roundoff_j: f64,
}

impl FusionNeutronActivationState {
    /// Bind an existing finite-material checkpoint to its original energy model.
    /// Only reaction coefficients bound to the original burnup reactions may
    /// differ from the original diffusion model. Group identity, effective
    /// speeds, geometry, reaction order/kinds, fixed backgrounds, links and
    /// leakage must match exactly. No input history is reset or replayed.
    ///
    /// Energy exports tracked by this adapter start at zero at binding; earlier
    /// native inventory accounts remain accessible through `spatial()`.
    /// Authored reaction energy must exclude incident neutron kinetic energy.
    pub fn new(
        model: Arc<NeutronKineticModel>,
        spatial: DiffusionBurnupState,
    ) -> Result<Self, FusionNeutronActivationError> {
        Self::prepare(model, spatial, false)
    }

    /// Admit funded aggregate induced-fission families whose entire product
    /// cohort lies outside this finite nuclide inventory.
    ///
    /// Each positive induced-family source must already carry the kinetic
    /// model's explicit future-neutron conversion declaration. Every bound
    /// fission channel producing these families must retain no positive nuclide
    /// products, including stable products; deleting only an overlapping neutron
    /// emission channel would still duplicate fragment energy or ownership.
    /// Unbound reactions declare fixed background sources with unrepresented
    /// product cohorts. Unrelated radioactive inventories remain admissible.
    ///
    /// The declaration applies to dormant and active original channels alike.
    /// Neither finite stock nor an accepted precursor/transport history is reset.
    pub fn new_with_unrepresented_fission_precursors(
        model: Arc<NeutronKineticModel>,
        spatial: DiffusionBurnupState,
    ) -> Result<Self, FusionNeutronActivationError> {
        Self::prepare(model, spatial, true)
    }

    fn prepare(
        model: Arc<NeutronKineticModel>,
        spatial: DiffusionBurnupState,
        unrepresented_fission_products: bool,
    ) -> Result<Self, FusionNeutronActivationError> {
        if !unrepresented_fission_products
            && model
                .fission_precursor_energy_j()
                .iter()
                .any(|&energy| energy > 0.0)
        {
            return Err(FusionNeutronActivationError::Invalid(
                "aggregate induced-fission families require an explicit unrepresented-product boundary",
            ));
        }
        let mut bound_reactions = vec![false; model.diffusion().reactions().len()];
        for binding in spatial.bindings() {
            let bound = bound_reactions.get_mut(binding.diffusion_reaction).ok_or(
                FusionNeutronActivationError::Invalid(
                    "activation binding outside original kinetic reaction layout",
                ),
            )?;
            *bound = true;
        }
        admit_topology(
            model.diffusion(),
            spatial.neutrons().model(),
            &bound_reactions,
        )?;
        model.admit_active_reactions(spatial.neutrons().model())?;
        let decay_partitions = DecayPartitions::new(&model, &spatial)?;
        let cells = spatial.inventories().len();
        let mut maximum_reaction_energy_j = spatial
            .inventories()
            .iter()
            .map(|inventory| {
                inventory
                    .model()
                    .reactions()
                    .iter()
                    .filter(|reaction| {
                        decay_partitions.is_some()
                            || !matches!(reaction.rate, BurnupRate::Decay { .. })
                    })
                    .map(|reaction| reaction.energy_per_event.total_j())
                    .fold(0.0_f64, f64::max)
            })
            .collect::<Vec<_>>();
        for binding in spatial.bindings() {
            let neutron = &model.diffusion().reactions()[binding.diffusion_reaction];
            let nuclear =
                &spatial.inventories()[neutron.cell].model().reactions()[binding.burnup_reaction];
            let q = model.reaction_energies()[binding.diffusion_reaction].nuclear_q_j;
            let inventory_q = nuclear.energy_per_event.total_j();
            let precursor_energy = model.fission_precursor_energy_j()[binding.diffusion_reaction];
            if precursor_energy > 0.0
                && nuclear
                    .products
                    .iter()
                    .any(|product| product.yield_per_event > 0.0)
            {
                return Err(FusionNeutronActivationError::Invalid(
                    "aggregate induced-fission families require every positive retained fission product to be unrepresented",
                ));
            }
            if (q != 0.0
                || matches!(
                    neutron.kind,
                    crate::diffusion::DiffusionReactionKind::NeutronProduction { .. }
                        | crate::diffusion::DiffusionReactionKind::Fission { .. }
                ))
                && inventory_q != 0.0
            {
                return Err(FusionNeutronActivationError::Invalid(
                    "prompt Q in neutron production or the neutron energy model requires zero burnup event energy",
                ));
            }
            maximum_reaction_energy_j[neutron.cell] = maximum_reaction_energy_j[neutron.cell]
                .max(checked(q.abs() + precursor_energy + inventory_q)?);
        }
        if decay_partitions.is_some() {
            let maximum_neutron_energy = model
                .group_energies_j()
                .iter()
                .copied()
                .fold(0.0_f64, f64::max);
            for energy in &mut maximum_reaction_energy_j {
                *energy = checked(*energy + maximum_neutron_energy)?;
            }
        }
        let kinetic = NeutronKineticState {
            model,
            neutrons: spatial.neutrons().clone(),
            released: vec![[Accumulator::new(0.0); 3]; cells],
            nuclear_reaction_energy: vec![Accumulator::new(0.0); cells],
            induced_precursor_energy: vec![Accumulator::new(0.0); cells],
            exported_precursor_energy: vec![Accumulator::new(0.0); cells],
        };
        kinetic.stored_energy_j()?;
        Ok(Self {
            spatial,
            kinetic,
            bound_reactions: bound_reactions.into(),
            maximum_reaction_energy_j: maximum_reaction_energy_j.into(),
            decay_partitions,
            nuclear_released: vec![[Accumulator::new(0.0); 3]; cells],
            released: vec![[Accumulator::new(0.0); 3]; cells],
        })
    }

    /// Original energy/source interpretation, retained as cross sections evolve.
    pub fn model(&self) -> &Arc<NeutronKineticModel> {
        self.kinetic.model()
    }

    /// Original spatial and product-inventory owner at this accepted checkpoint.
    pub const fn spatial(&self) -> &DiffusionBurnupState {
        &self.spatial
    }

    /// Neutron checkpoint from the original coupled owner.
    pub fn neutrons(&self) -> &DiffusionState {
        self.spatial.neutrons()
    }

    /// Finite targets and products, including their accepted irradiation history.
    pub fn inventories(&self) -> &[BurnupInventory] {
        self.spatial.inventories()
    }

    /// Retained neutron kinetic energy [J], without a future nuclear-Q reservoir.
    pub fn kinetic_energy_j(&self) -> Result<f64, NeutronKineticError> {
        self.kinetic.kinetic_energy_j()
    }

    /// Future-neutron energy retained in the original precursor populations [J].
    pub fn precursor_energy_j(&self) -> Result<f64, NeutronKineticError> {
        self.kinetic.precursor_energy_j()
    }

    /// Current neutron kinetic plus reserved precursor energy [J].
    pub fn stored_energy_j(&self) -> Result<f64, NeutronKineticError> {
        self.kinetic.stored_energy_j()
    }

    /// Explicit induced-fission future-neutron conversion accepted since binding [J].
    /// It is already stored in families or converted to neutrons, not extra heat.
    pub fn cumulative_induced_precursor_energy_j(&self) -> Result<f64, NeutronKineticError> {
        self.kinetic.cumulative_induced_precursor_energy_j()
    }

    /// Per-cell induced-fission future-neutron conversion since binding [J].
    pub fn cell_induced_precursor_energy_j(&self, cell: usize) -> Option<f64> {
        self.kinetic.cell_induced_precursor_energy_j(cell)
    }

    /// Stored energy carried out by precursor flow since binding [J].
    pub fn cumulative_exported_precursor_energy_j(&self) -> Result<f64, NeutronKineticError> {
        self.kinetic.cumulative_exported_precursor_energy_j()
    }

    pub fn cell_exported_precursor_energy_j(&self, cell: usize) -> Option<f64> {
        self.kinetic.cell_exported_precursor_energy_j(cell)
    }

    /// Neutron-associated exports, including prompt Q, since binding.
    pub fn cell_kinetic_released(&self, cell: usize) -> Option<EnergyChannels> {
        self.kinetic.cell_released(cell)
    }

    /// Original inventory-event exports for one cell since binding. Prompt Q
    /// owned by the neutron model is already in the neutron-associated exports.
    pub fn cell_nuclear_released(&self, cell: usize) -> Option<EnergyChannels> {
        self.nuclear_released
            .get(cell)
            .copied()
            .map(crate::channel_values)
    }

    /// Separate signed prompt Q owned by the original neutron energy model.
    pub fn cell_prompt_nuclear_energy_j(&self, cell: usize) -> Option<f64> {
        self.kinetic.cell_nuclear_reaction_energy_j(cell)
    }

    /// Combined exports for one original cell since this adapter was bound.
    pub fn cell_released(&self, cell: usize) -> Option<EnergyChannels> {
        self.released.get(cell).copied().map(crate::channel_values)
    }

    /// All combined accepted exports in original diffusion-cell order.
    pub fn cumulative_cell_released(&self) -> Vec<EnergyChannels> {
        self.released
            .iter()
            .copied()
            .map(crate::channel_values)
            .collect()
    }

    /// Evolve constant extensive neutron rates [particles/s] on the original
    /// coupled neutron/nuclide owner. There is no implicit target supply. Zero
    /// rates continue absorption and product decay after fusion switches off.
    /// Temperatures belong to the original inventory models and are frozen over
    /// this interval, including evaluation of their processed reaction tables.
    pub fn trial_with_rates(
        &self,
        duration_s: f64,
        temperatures_k: &[f64],
        neutron_rates_s: &[f64],
        budget: FusionNeutronActivationBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, FusionNeutronActivationStep), FusionNeutronActivationError> {
        self.trial_with_sources(
            duration_s,
            temperatures_k,
            neutron_rates_s,
            &vec![0.0; self.neutrons().precursor_counts().len()],
            budget,
            checkpoint,
        )
    }

    /// Evolve supplied neutron and delayed-family rates on the same finite stock.
    /// External precursors carry their explicitly admitted stored neutron energy;
    /// internally born families reserve energy from their original decay events.
    pub fn trial_with_sources(
        &self,
        duration_s: f64,
        temperatures_k: &[f64],
        neutron_rates_s: &[f64],
        precursor_rates_s: &[f64],
        budget: FusionNeutronActivationBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, FusionNeutronActivationStep), FusionNeutronActivationError> {
        checkpoint()?;
        if !nonnegative(budget.energy_tolerance_j) {
            return Err(FusionNeutronActivationError::Invalid(
                "finite nonnegative combined energy allowance required",
            ));
        }
        if neutron_rates_s.len() != self.neutrons().neutron_counts().len()
            || precursor_rates_s.len() != self.neutrons().precursor_counts().len()
        {
            return Err(FusionNeutronActivationError::Invalid(
                "one external neutron and precursor rate per original cell/group/family required",
            ));
        }
        let sources = DiffusionBurnupSources {
            neutron_rates_s: neutron_rates_s.to_vec(),
            precursor_rates_s: precursor_rates_s.to_vec(),
            nuclide_rates_s: self
                .inventories()
                .iter()
                .map(|inventory| vec![0.0; inventory.populations().len()])
                .collect(),
        };
        let (spatial, spatial_step) = self.spatial.trial_with_sources(
            duration_s,
            temperatures_k,
            &sources,
            budget.spatial,
            &mut checkpoint,
        )?;
        admit_topology(
            self.model().diffusion(),
            spatial.neutrons().model(),
            &self.bound_reactions,
        )?;
        let (kinetic, kinetic_step) = self.kinetic.account_accepted_diffusion(
            spatial.neutrons().clone(),
            spatial_step.neutrons.clone(),
            budget.energy_tolerance_j,
            &mut checkpoint,
        )?;
        let transfer = self
            .decay_partitions
            .as_ref()
            .map(|partitions| {
                partitions.account(
                    self.model(),
                    self.inventories(),
                    &spatial_step,
                    duration_s,
                    neutron_rates_s,
                    precursor_rates_s,
                    &mut checkpoint,
                )
            })
            .transpose()?;
        let injection_number_roundoff =
            transfer.as_ref().map_or(0.0, |value| value.number_roundoff);
        if checked(spatial_step.neutrons.population_roundoff + injection_number_roundoff)?
            > budget.spatial.diffusion.population_tolerance
        {
            return Err(FusionNeutronActivationError::Balance(
                "external/internal neutron injection exceeds allowance",
            ));
        }
        let external_injected_neutrons = transfer
            .as_ref()
            .map_or(spatial_step.neutrons.injected_neutrons, |value| {
                value.external_neutrons
            });
        let external_injected_kinetic_energy_j = transfer
            .as_ref()
            .map_or(kinetic_step.injected_kinetic_energy_j, |value| {
                value.external_energy_j
            });
        let external_injected_precursors = transfer
            .as_ref()
            .map_or(spatial_step.neutrons.injected_precursors, |value| {
                value.external_precursors
            });
        let external_injected_precursor_energy_j = transfer
            .as_ref()
            .map_or(kinetic_step.injected_precursor_energy_j, |value| {
                value.external_precursor_energy_j
            });
        let mut defect = Accumulator::new(kinetic_step.energy_roundoff_j);
        if let Some(transfer) = &transfer {
            defect.add(transfer.energy_roundoff_j)?;
        }
        let mut generated = Accumulator::new(kinetic_step.nuclear_reaction_energy_j);
        defect.add(
            generated
                .add(kinetic_step.induced_precursor_energy_j)?
                .abs(),
        )?;
        let mut reaction_defect = Accumulator::new(0.0);
        let mut nuclear_totals = [Accumulator::new(0.0); 3];
        let mut totals = [Accumulator::new(0.0); 3];
        let mut nuclear_released = self.nuclear_released.clone();
        let mut released = self.released.clone();
        let mut cell_released = Vec::with_capacity(spatial.inventories().len());
        for (cell, nuclear) in spatial_step.burnup.iter().enumerate() {
            checkpoint()?;
            defect.add(nuclear.energy_roundoff_j)?;
            defect.add(generated.add(nuclear.generated_energy_j)?.abs())?;
            let count_defect = spatial_step.reaction_defect_by_cell[cell];
            let event_energy = self.maximum_reaction_energy_j[cell];
            let coupling = checked(count_defect * event_energy)?;
            reaction_defect.add(coupling)?;
            defect.add(count_defect.mul_add(event_energy, -coupling).abs())?;
            let nuclear_channels = channels(nuclear.released);
            let mut exported_nuclear_channels = nuclear_channels;
            if let Some(transfer) = &transfer {
                exported_nuclear_channels[1] = transfer.remaining_radiation_j[cell];
            }
            let kinetic_channels = channels(kinetic_step.cell_released[cell]);
            let mut combined = [Accumulator::new(0.0); 3];
            for channel in 0..3 {
                defect.add(
                    combined[channel]
                        .add(exported_nuclear_channels[channel])?
                        .abs(),
                )?;
                defect.add(combined[channel].add(kinetic_channels[channel])?.abs())?;
                defect.add(
                    nuclear_released[cell][channel]
                        .add(nuclear_channels[channel])?
                        .abs(),
                )?;
                defect.add(
                    released[cell][channel]
                        .add(combined[channel].value())?
                        .abs(),
                )?;
                defect.add(
                    nuclear_totals[channel]
                        .add(nuclear_channels[channel])?
                        .abs(),
                )?;
                defect.add(totals[channel].add(combined[channel].value())?.abs())?;
            }
            crate::channel_values(nuclear_released[cell]).validate()?;
            crate::channel_values(released[cell]).validate()?;
            cell_released.push(crate::channel_values(combined));
        }
        defect.add(reaction_defect.value())?;
        let exports = crate::channel_values(totals);
        let energy_balance_residual_j = checked(sum(&[
            kinetic_step.final_kinetic_energy_j,
            -kinetic_step.initial_kinetic_energy_j,
            kinetic_step.final_precursor_energy_j,
            -kinetic_step.initial_precursor_energy_j,
            kinetic_step.exported_precursor_energy_j,
            exports.deposited_j,
            exports.radiation_j,
            exports.escaped_j,
            -external_injected_kinetic_energy_j,
            -external_injected_precursor_energy_j,
            -generated.value(),
        ]))?;
        defect.add(energy_balance_residual_j.abs())?;
        if defect.value() > budget.energy_tolerance_j {
            return Err(FusionNeutronActivationError::Balance(
                "combined kinetic, nuclear and reaction-coupling energy account exceeds allowance",
            ));
        }
        let report = FusionNeutronActivationStep {
            spatial: spatial_step,
            kinetic: kinetic_step,
            nuclear_released: crate::channel_values(nuclear_totals),
            decay_neutrons: transfer.as_ref().map_or(0.0, |value| value.decay_neutrons),
            decay_precursors: transfer
                .as_ref()
                .map_or(0.0, |value| value.decay_precursors),
            transferred_decay_neutron_energy_j: transfer
                .as_ref()
                .map_or(0.0, |value| value.transferred_energy_j),
            reserved_decay_precursor_energy_j: transfer
                .as_ref()
                .map_or(0.0, |value| value.reserved_precursor_energy_j),
            external_injected_neutrons,
            external_injected_precursors,
            external_injected_kinetic_energy_j,
            external_injected_precursor_energy_j,
            injection_number_roundoff,
            nuclear_generated_energy_j: generated.value(),
            cell_released,
            released: exports,
            reaction_energy_defect_j: reaction_defect.value(),
            energy_balance_residual_j,
            energy_roundoff_j: defect.value(),
        };
        checkpoint()?;
        Ok((
            Self {
                spatial,
                kinetic,
                bound_reactions: Arc::clone(&self.bound_reactions),
                maximum_reaction_energy_j: Arc::clone(&self.maximum_reaction_energy_j),
                decay_partitions: self.decay_partitions.clone(),
                nuclear_released,
                released,
            },
            report,
        ))
    }

    /// Evolve complete original fusion emissions over their shared interval.
    /// Their channel-mean, interval-uniform projection remains explicit, with
    /// the original emissions retained in the returned kinetic report.
    pub fn trial_with_fusion(
        &self,
        duration_s: f64,
        temperatures_k: &[f64],
        sources: &[FusionNeutronSource],
        budget: FusionNeutronActivationBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, FusionNeutronActivationStep), FusionNeutronActivationError> {
        checkpoint()?;
        let mut rates = vec![Accumulator::new(0.0); self.neutrons().neutron_counts().len()];
        let mut source_count = Accumulator::new(0.0);
        let mut source_energy = Accumulator::new(0.0);
        for source in sources {
            checkpoint()?;
            if source.model().as_ref() != self.model().as_ref()
                || source.emission().duration().value() != duration_s
            {
                return Err(FusionNeutronActivationError::Invalid(
                    "fusion source must match the original energy model and complete interval",
                ));
            }
            source_count.add(source.emission().total_particle_count())?;
            source_energy.add(source.emission().total_kinetic_energy().value())?;
            for (rate, &value) in rates.iter_mut().zip(source.neutron_rates_s()) {
                rate.add(value)?;
            }
        }
        let rates: Vec<_> = rates.into_iter().map(Accumulator::value).collect();
        let (next, mut report) =
            self.trial_with_rates(duration_s, temperatures_k, &rates, budget, &mut checkpoint)?;
        let number_residual = checked(report.external_injected_neutrons - source_count.value())?;
        let energy_residual =
            checked(report.external_injected_kinetic_energy_j - source_energy.value())?;
        report.kinetic.fusion_number_residual = number_residual;
        report.kinetic.fusion_energy_residual_j = energy_residual;
        report.kinetic.energy_roundoff_j =
            checked(report.kinetic.energy_roundoff_j + energy_residual.abs())?;
        report.energy_roundoff_j = checked(report.energy_roundoff_j + energy_residual.abs())?;
        if number_residual.abs()
            + report.spatial.neutrons.population_roundoff
            + report.injection_number_roundoff
            > budget.spatial.diffusion.population_tolerance
            || report.energy_roundoff_j > budget.energy_tolerance_j
        {
            return Err(FusionNeutronActivationError::Balance(
                "fusion-to-activation source number or energy account exceeds allowance",
            ));
        }
        report.kinetic.fusion_sources = sources.to_vec();
        checkpoint()?;
        Ok((next, report))
    }
}

/// Atomic source, topology, native propagation or combined-energy refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum FusionNeutronActivationError {
    Invalid(&'static str),
    Balance(&'static str),
    Kinetic(NeutronKineticError),
    Spatial(DiffusionBurnupError),
    Source(FissionError),
}

impl From<NeutronKineticError> for FusionNeutronActivationError {
    fn from(error: NeutronKineticError) -> Self {
        Self::Kinetic(error)
    }
}
impl From<DiffusionBurnupError> for FusionNeutronActivationError {
    fn from(error: DiffusionBurnupError) -> Self {
        Self::Spatial(error)
    }
}
impl From<FissionError> for FusionNeutronActivationError {
    fn from(error: FissionError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for FusionNeutronActivationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fusion neutron activation: {self:?}")
    }
}
impl std::error::Error for FusionNeutronActivationError {}

fn admit_topology(
    original: &DiffusionModel,
    current: &DiffusionModel,
    bound_reactions: &[bool],
) -> Result<(), FusionNeutronActivationError> {
    let a = original.data();
    let b = current.data();
    if a.data_id != b.data_id
        || a.groups != b.groups
        || a.speeds_m_s != b.speeds_m_s
        || a.cell_volumes_m3 != b.cell_volumes_m3
        || a.delayed_groups != b.delayed_groups
        || a.links != b.links
        || a.leakage != b.leakage
        || original.precursor_transport() != current.precursor_transport()
        || a.reactions.len() != b.reactions.len()
        || a.reactions.len() != bound_reactions.len()
    {
        return Err(FusionNeutronActivationError::Invalid(
            "activation must retain original group, speed, geometry and transport identities",
        ));
    }
    for ((a, b), &bound) in a.reactions.iter().zip(&b.reactions).zip(bound_reactions) {
        if a.data_id != b.data_id
            || a.cell != b.cell
            || a.incident_group != b.incident_group
            || a.kind != b.kind
            || (!bound && a.macroscopic_cross_section_m_inv != b.macroscopic_cross_section_m_inv)
        {
            return Err(FusionNeutronActivationError::Invalid(
                "only original burnup-bound capture coefficients may change",
            ));
        }
    }
    Ok(())
}

fn channels(value: EnergyChannels) -> [f64; 3] {
    [value.deposited_j, value.radiation_j, value.escaped_j]
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod precursor_tests;

#[cfg(test)]
mod induced_tests;

#[cfg(test)]
pub(crate) use induced_tests::material_induced_precursor_fixture;

#[cfg(test)]
pub(crate) use precursor_tests::{
    expected_history as precursor_expected_history, material_precursor_fixture,
};
