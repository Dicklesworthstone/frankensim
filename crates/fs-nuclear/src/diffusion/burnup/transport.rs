//! Optional joint material movement inside the original spatial Picard solve.

use super::*;
use crate::burnup::transport::{
    BurnupTransportBudget, BurnupTransportLimits, BurnupTransportState, BurnupTransportStep,
};
use crate::irradiation::GroupFlux;

/// Movement and residence from the existing joint depletion owner. Path order
/// is the original `BurnupTransportState` order; these are actual expected
/// nuclei, not neutron amplitudes. Nuclear events remain in the per-cell reports.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionBurnupTransportStep {
    /// Cell-major, then original nuclide order [nuclei s].
    pub integrated_populations_s: Vec<Vec<f64>>,
    pub transfer_counts: Vec<Vec<f64>>,
    pub outflow_counts: Vec<Vec<f64>>,
    /// Whole-cohort balances with internal movement cancelled [nuclei].
    pub cohort_balance_residuals: Vec<f64>,
}

// There is one population/history owner, never a second synchronized copy.
#[derive(Clone, Debug)]
pub(super) enum InventoryOwner {
    Independent(Vec<BurnupInventory>),
    Transported {
        state: BurnupTransportState,
        limits: BurnupTransportLimits,
    },
}

pub(super) struct InventoryTrial {
    pub owner: InventoryOwner,
    pub reports: Vec<BurnupStep>,
    pub means: Vec<Vec<f64>>,
    pub produced: Vec<Vec<f64>>,
    pub produced_by_cell: Vec<f64>,
    pub transport: Option<DiffusionBurnupTransportStep>,
}

impl InventoryTrial {
    /// Preserve the original transport report when its trajectory was solved
    /// either by its native map or inside the larger neutron/nuclide block.
    pub(super) fn from_transport(
        initial: &BurnupTransportState,
        candidate: BurnupTransportState,
        limits: BurnupTransportLimits,
        report: BurnupTransportStep,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DiffusionBurnupError> {
        let mut means = Vec::with_capacity(initial.inventories().len());
        for (inventory, residence) in initial
            .inventories()
            .iter()
            .zip(&report.integrated_populations_s)
        {
            let mut mean = Vec::with_capacity(residence.len());
            for (&population, &integral) in inventory.populations().iter().zip(residence) {
                checkpoint()?;
                // Residence retains continuous ingrowth and movement even
                // without a local removal hazard. Endpoints cannot replace it.
                mean.push(if report.duration_s == 0.0 {
                    population
                } else {
                    mean_population(integral, 1.0, report.duration_s)?
                });
            }
            means.push(mean);
        }
        let BurnupTransportStep {
            burnup,
            integrated_populations_s,
            produced_nuclei,
            produced_nuclei_by_compartment,
            transfer_counts,
            outflow_counts,
            cohort_balance_residuals,
            ..
        } = report;
        Ok(Self {
            owner: InventoryOwner::Transported {
                state: candidate,
                limits,
            },
            reports: burnup,
            means,
            produced: produced_nuclei,
            produced_by_cell: produced_nuclei_by_compartment,
            transport: Some(DiffusionBurnupTransportStep {
                integrated_populations_s,
                transfer_counts,
                outflow_counts,
                cohort_balance_residuals,
            }),
        })
    }
}

impl DiffusionBurnupState {
    /// Couple the existing finite-material transport owner to evolving spatial
    /// neutrons. Original cell volumes and ordering must match exactly. Every
    /// neutron-induced reaction retains the same complete binding requirements.
    ///
    /// The supplied limits admit the initial joint owner. Its compartment,
    /// path and coefficient ceilings remain fixed; each subsequent step uses
    /// that call's `budget.burnup.limits` for the complete joint propagation and
    /// its burnup account allowances for the whole cohort. Thus callers can
    /// tighten numerical work without rebuilding a physical history.
    ///
    /// Transport rates stay at the supplied state's accepted values. Carrier
    /// fluid dynamics, bulk mass and sensible/latent enthalpy are not inferred.
    pub fn new_with_transport(
        neutrons: DiffusionState,
        inventories: BurnupTransportState,
        bindings: Vec<ReactionBinding>,
        limits: BurnupTransportLimits,
    ) -> Result<Self, DiffusionBurnupError> {
        if inventories.volumes_m3() != neutrons.model().cell_volumes_m3() {
            return Err(DiffusionBurnupError::Invalid(
                "transport compartments must match original spatial cell volumes and order",
            ));
        }
        inventories.admit(inventories.transport(), limits, &mut || Ok(()))?;
        let state = Self {
            neutrons,
            inventories: InventoryOwner::Transported {
                state: inventories,
                limits,
            },
            bindings: bindings.into(),
            decay_neutrons: None,
        };
        state.validate_bindings()?;
        Ok(state)
    }

    /// The original accepted joint owner, including preparation and lifetime
    /// movement accounts. Independent-cell states have no implicit transport.
    #[must_use]
    pub fn transport(&self) -> Option<&BurnupTransportState> {
        match &self.inventories {
            InventoryOwner::Independent(_) => None,
            InventoryOwner::Transported { state, .. } => Some(state),
        }
    }

    /// Admit the complete joint shape under the current numerical work limits,
    /// including native thermal bindings and idle advances, without a solve.
    pub(crate) fn admit_transport(
        &self,
        budget: BurnupBudget,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), DiffusionBurnupError> {
        if let InventoryOwner::Transported { state, limits } = &self.inventories {
            state.admit(
                state.transport(),
                BurnupTransportLimits {
                    burnup: budget.limits,
                    ..*limits
                },
                checkpoint,
            )?;
        }
        Ok(())
    }
}

impl InventoryOwner {
    pub(super) fn inventories(&self) -> &[BurnupInventory] {
        match self {
            Self::Independent(inventories) => inventories,
            Self::Transported { state, .. } => state.inventories(),
        }
    }

    pub(super) fn trial(
        &self,
        duration_s: f64,
        fluxes: &[GroupFlux],
        temperatures_k: &[f64],
        sources: Option<&[Vec<f64>]>,
        budget: BurnupBudget,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<InventoryTrial, DiffusionBurnupError> {
        if let Self::Transported { state, limits } = self {
            let joint_budget = BurnupTransportBudget {
                limits: BurnupTransportLimits {
                    burnup: budget.limits,
                    ..*limits
                },
                nuclei_tolerance: budget.nuclei_tolerance,
                energy_tolerance_j: budget.energy_tolerance_j,
            };
            let (candidate, report) = if let Some(sources) = sources {
                state.trial_with_inputs(
                    duration_s,
                    fluxes,
                    temperatures_k,
                    state.transport(),
                    sources,
                    joint_budget,
                    &mut *checkpoint,
                )?
            } else {
                state.trial(
                    duration_s,
                    fluxes,
                    temperatures_k,
                    joint_budget,
                    &mut *checkpoint,
                )?
            };
            return InventoryTrial::from_transport(state, candidate, *limits, report, checkpoint);
        }

        let inventories = self.inventories();
        let mut candidates = Vec::with_capacity(inventories.len());
        let mut reports = Vec::with_capacity(inventories.len());
        let mut means = Vec::with_capacity(inventories.len());
        let mut produced = Vec::with_capacity(inventories.len());
        let mut produced_by_cell = Vec::with_capacity(inventories.len());
        for (cell, inventory) in inventories.iter().enumerate() {
            checkpoint()?;
            let flux = &fluxes[cell];
            let rates = sources
                .map(|sources| &sources[cell])
                .filter(|rates| rates.iter().any(|&rate| rate > 0.0));
            let (candidate, mut report, supplied) = if let Some(rates) = rates {
                let (candidate, production) = inventory.trial_with_production(
                    duration_s,
                    flux,
                    temperatures_k[cell],
                    rates,
                    budget,
                    &mut *checkpoint,
                )?;
                (candidate, production.burnup, production.produced_nuclei)
            } else {
                let (candidate, report) = inventory.trial(
                    duration_s,
                    flux,
                    temperatures_k[cell],
                    budget,
                    &mut *checkpoint,
                )?;
                (candidate, report, vec![0.0; inventory.populations().len()])
            };
            let produced_total = export_supplied_count(&supplied, &mut report, budget, checkpoint)?;
            let rates = inventory.model().rate_constants(
                flux,
                temperatures_k[cell],
                budget.limits,
                &mut *checkpoint,
            )?;
            let mut removal = vec![Accumulator::new(0.0); inventory.populations().len()];
            let mut hazards = vec![Accumulator::new(0.0); inventory.populations().len()];
            for ((reaction, &rate), &count) in inventory
                .model()
                .reactions()
                .iter()
                .zip(&rates)
                .zip(&report.reaction_counts)
            {
                checkpoint()?;
                removal[reaction.parent].add(count)?;
                hazards[reaction.parent].add(rate)?;
            }
            let mut mean = inventory.populations().to_vec();
            if duration_s > 0.0 {
                for parent in 0..mean.len() {
                    checkpoint()?;
                    mean[parent] = if hazards[parent].value() > 0.0 {
                        mean_population(
                            removal[parent].value(),
                            hazards[parent].value(),
                            duration_s,
                        )?
                    } else {
                        // Independent zero-hazard nodes have only direct feed
                        // relevant to this frozen neutron operator.
                        supplied_mean(inventory.populations()[parent], supplied[parent])?
                    };
                }
            }
            candidates.push(candidate);
            reports.push(report);
            means.push(mean);
            produced.push(supplied);
            produced_by_cell.push(produced_total);
        }
        Ok(InventoryTrial {
            owner: Self::Independent(candidates),
            reports,
            means,
            produced,
            produced_by_cell,
            transport: None,
        })
    }
}
