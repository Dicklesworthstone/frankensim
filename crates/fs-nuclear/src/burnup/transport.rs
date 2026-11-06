//! Finite radioactive inventories moving between well-mixed compartments.
//!
//! A directed rate k [1/s] transfers k*N nuclei/s from a finite donor. For a
//! prescribed carrier flow Q, the well-mixed coefficient is Q/V_donor. Rates
//! may differ by nuclide, allowing explicit selective extraction. They remain
//! fixed during one admitted interval; this owner does not solve carrier flow,
//! bulk mass, material composition, or transported sensible/latent enthalpy.
//!
//! Local decay/capture/fission and transport share one positive block operator.
//! Newly supplied or reaction-created nuclei can move and react again during the
//! SAME interval. Integrated residence attributes actual nuclear events and
//! directed transfers separately; movement is never called a decay or production.
//! Event energy belongs to the compartment where that event occurs. Already
//! emitted heat and lifetime reaction histories stay with their original owners.
//!
//! The original BurnupModel, BurnupInventory ledger and fs_la reaction propagators
//! are reused. The flattened (compartment, nuclide) layout retains the numerical
//! owner's 64-node bound. There is no second nuclear solver or inferred isotope
//! data. Coupled block depletion and donor-proportional transfers are described
//! at https://docs.openmc.org/en/stable/methods/depletion.html#transfer-rates .

use std::sync::Arc;

use fs_la::decay::{
    MAX_DECAY_EDGES, MAX_DECAY_NODES,
    reaction::{ReactionAdvance, ReactionEdge, ReactionPlan, source::ReactionSourcePlan},
};

use super::{
    Accumulator, BurnupBudget, BurnupError, BurnupInventory, BurnupLimits, BurnupModel, BurnupStep,
    TransportAccounts, add, checked, nonnegative_bound, product, projection_error,
};
use crate::irradiation::GroupFlux;
use crate::{FissionError, nonnegative, positive, sum};

/// Directed transfer of the same nuclide between two represented compartments.
/// Each coefficient is a finite nonnegative per-donor-nucleus rate [1/s].
/// Repeated endpoints are distinct additive paths with separate lifetime counts.
#[derive(Clone, Debug, PartialEq)]
pub struct NuclideTransfer {
    pub from_compartment: usize,
    pub to_compartment: usize,
    pub rates_per_s: Vec<f64>,
}

/// Removal into an untracked exterior. Later exterior decays release no energy
/// in a retained compartment. An outflow is not a nuclear reaction or radiation.
#[derive(Clone, Debug, PartialEq)]
pub struct NuclideOutflow {
    pub compartment: usize,
    pub rates_per_s: Vec<f64>,
}

/// Authored path order is stable across rate changes and lifetime accounts.
/// Empty paths keep nuclei in their original compartment.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NuclideTransport {
    pub transfers: Vec<NuclideTransfer>,
    pub outflows: Vec<NuclideOutflow>,
}

impl NuclideTransport {
    /// Rate values may change while these endpoints retain their original order.
    pub fn same_topology(&self, other: &Self) -> bool {
        self.transfers.len() == other.transfers.len()
            && self.outflows.len() == other.outflows.len()
            && self.transfers.iter().zip(&other.transfers).all(|(a, b)| {
                (a.from_compartment, a.to_compartment) == (b.from_compartment, b.to_compartment)
            })
            && self
                .outflows
                .iter()
                .zip(&other.outflows)
                .all(|(a, b)| a.compartment == b.compartment)
    }
}

/// Bounds cover the complete cohort, in addition to each original model's limits.
/// The propagation node/edge/work bounds apply to the flattened joint operator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BurnupTransportLimits {
    pub burnup: BurnupLimits,
    pub max_compartments: usize,
    pub max_transfers: usize,
    pub max_outflows: usize,
    /// Cell reaction rates plus every per-nuclide transfer/outflow coefficient.
    pub max_rate_coefficients: usize,
}

impl Default for BurnupTransportLimits {
    fn default() -> Self {
        Self {
            burnup: BurnupLimits::default(),
            max_compartments: MAX_DECAY_NODES,
            max_transfers: MAX_DECAY_EDGES,
            max_outflows: MAX_DECAY_EDGES,
            max_rate_coefficients: 65_536,
        }
    }
}

/// Noncancelling extensive account allowances for the entire accepted cohort.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BurnupTransportBudget {
    pub limits: BurnupTransportLimits,
    pub nuclei_tolerance: f64,
    pub energy_tolerance_j: f64,
}

impl BurnupTransportBudget {
    fn local(self) -> Result<BurnupBudget, BurnupError> {
        if !nonnegative(self.nuclei_tolerance) || !nonnegative(self.energy_tolerance_j) {
            return Err(BurnupError::Invalid(
                "finite nonnegative whole-cohort account tolerances required",
            ));
        }
        Ok(BurnupBudget {
            limits: self.limits.burnup,
            nuclei_tolerance: self.nuclei_tolerance,
            energy_tolerance_j: self.energy_tolerance_j,
        })
    }
}

#[derive(Clone, Debug)]
enum JointPlan {
    Closed(ReactionPlan),
    Source(ReactionSourcePlan),
}

impl JointPlan {
    fn prepare(
        rates: &[f64],
        edges: &[ReactionEdge],
        duration_s: f64,
        limits: BurnupTransportLimits,
        with_sources: bool,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, BurnupError> {
        if with_sources {
            Ok(Self::Source(owner_call(checkpoint, |cancelled| {
                ReactionSourcePlan::prepare_with_population_integrals(
                    rates,
                    edges,
                    duration_s,
                    limits.burnup.propagation,
                    cancelled,
                )
            })?))
        } else {
            Ok(Self::Closed(owner_call(checkpoint, |cancelled| {
                ReactionPlan::prepare_with_population_integrals(
                    rates,
                    edges,
                    duration_s,
                    limits.burnup.propagation,
                    cancelled,
                )
            })?))
        }
    }

    fn homogeneous(&self) -> &ReactionPlan {
        match self {
            Self::Closed(plan) => plan,
            Self::Source(plan) => plan.homogeneous_plan(),
        }
    }

    fn supports_sources(&self) -> bool {
        matches!(self, Self::Source(_))
    }

    fn apply(
        &self,
        initial: &[f64],
        supplied: &[f64],
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(ReactionAdvance, Vec<f64>), BurnupError> {
        let evolved = owner_call(checkpoint, |cancelled| match self {
            Self::Closed(plan) => plan.apply_with_cancel(initial, cancelled),
            Self::Source(plan) => plan.apply_with_cancel(initial, supplied, cancelled),
        })?;
        let residence = owner_call(checkpoint, |cancelled| match self {
            Self::Closed(plan) => plan.integrate_populations_with_cancel(initial, cancelled),
            Self::Source(plan) => {
                plan.integrate_populations_with_cancel(initial, supplied, cancelled)
            }
        })?;
        Ok((evolved, residence))
    }
}

/// One immutable joint map shared across cohort clones. Injection magnitudes
/// and initial populations do not invalidate it. Current flux/temperature domains
/// are still checked before reuse, including equivalent complete rate vectors.
#[derive(Clone, Debug)]
pub struct PreparedBurnupTransport {
    model: Arc<BurnupModel>,
    transport: Arc<NuclideTransport>,
    rates: Vec<Vec<f64>>,
    limits: BurnupTransportLimits,
    plan: JointPlan,
}

impl PreparedBurnupTransport {
    pub fn model(&self) -> &Arc<BurnupModel> {
        &self.model
    }
    pub fn transport(&self) -> &NuclideTransport {
        &self.transport
    }
    pub fn reaction_rate_constants(&self) -> &[Vec<f64>] {
        &self.rates
    }
    pub fn duration_s(&self) -> f64 {
        self.plan.homogeneous().duration_s()
    }
    pub fn preparation_work(&self) -> u64 {
        self.plan.homogeneous().preparation_work()
    }
    pub fn supports_sources(&self) -> bool {
        self.plan.supports_sources()
    }
}

/// Original inventories plus joint transport history, in fixed compartment order.
/// All inventories use the same exact admitted model Arc. Clones share only
/// immutable data/maps; every population, reaction and flow counter is independent.
#[derive(Clone, Debug)]
pub struct BurnupTransportState {
    model: Arc<BurnupModel>,
    inventories: Vec<BurnupInventory>,
    volumes_m3: Arc<[f64]>,
    transport: Arc<NuclideTransport>,
    transfers: Vec<Vec<Accumulator>>,
    outflows: Vec<Vec<Accumulator>>,
    prepared: Option<Arc<PreparedBurnupTransport>>,
}

impl PartialEq for BurnupTransportState {
    fn eq(&self, other: &Self) -> bool {
        self.model == other.model
            && self.inventories == other.inventories
            && self.volumes_m3 == other.volumes_m3
            && self.transport == other.transport
            && self.transfers == other.transfers
            && self.outflows == other.outflows
    }
}

/// Physical event reports and separately integrated movement/source accounts.
/// Scalar whole-cohort defect totals are reductions of the per-cell BurnupStep
/// defects. Collective balance defects are charged once to compartment zero.
#[derive(Clone, Debug, PartialEq)]
pub struct BurnupTransportStep {
    pub duration_s: f64,
    pub burnup: Vec<BurnupStep>,
    /// Integral of each compartment's original nuclide population over this
    /// interval [nuclei s]. Includes same-interval ingrowth, supply and movement,
    /// also for nuclei with zero local nuclear or transport removal rate.
    pub integrated_populations_s: Vec<Vec<f64>>,
    pub produced_nuclei: Vec<Vec<f64>>,
    pub produced_nuclei_by_compartment: Vec<f64>,
    /// Path-major then original nuclide order; repeated circulation counts.
    pub transfer_counts: Vec<Vec<f64>>,
    pub outflow_counts: Vec<Vec<f64>>,
    /// Final-initial + physical removals - nuclear products - supplied + exterior
    /// outflow, summed over compartments for each original nuclide.
    pub cohort_balance_residuals: Vec<f64>,
    pub population_roundoff_nuclei: f64,
    pub energy_roundoff_j: f64,
}

impl BurnupTransportState {
    /// Bind finite current inventories and positive fixed control volumes.
    /// Earlier physical histories are preserved, not redistributed with material.
    pub fn new(
        inventories: Vec<BurnupInventory>,
        volumes_m3: Vec<f64>,
        transport: NuclideTransport,
        limits: BurnupTransportLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, BurnupError> {
        checkpoint()?;
        let model = Arc::clone(
            inventories
                .first()
                .ok_or(BurnupError::Invalid(
                    "a nonempty compartment inventory is required",
                ))?
                .model(),
        );
        let n = model.nuclides().len();
        // Check the bounded shape before converting caller-owned volumes into an Arc.
        let cells = inventories.len();
        let nodes = cells
            .checked_mul(n)
            .ok_or(BurnupError::Budget("joint inventory nodes"))?;
        if cells > limits.max_compartments
            || nodes > limits.burnup.propagation.maximum_nodes
            || nodes > MAX_DECAY_NODES
        {
            return Err(BurnupError::Budget("joint compartment/nuclide nodes"));
        }
        if volumes_m3.len() != cells {
            return Err(BurnupError::Invalid(
                "fixed compartment volumes and transfer path ordering required",
            ));
        }
        let mut state = Self {
            model,
            inventories,
            volumes_m3: volumes_m3.into(),
            transport: Arc::new(transport),
            transfers: Vec::new(),
            outflows: Vec::new(),
            prepared: None,
        };
        state.admit(&state.transport, limits, &mut checkpoint)?;
        // Admit the numerical owner's complete flattened shape/work configuration.
        let zero = vec![0.0; state.inventories.len() * n];
        owner_call(&mut checkpoint, |cancelled| {
            ReactionPlan::prepare(&zero, &[], 0.0, limits.burnup.propagation, cancelled)
        })?;
        state.transfers = vec![vec![Accumulator::new(0.0); n]; state.transport.transfers.len()];
        state.outflows = vec![vec![Accumulator::new(0.0); n]; state.transport.outflows.len()];
        checkpoint()?;
        Ok(state)
    }

    pub fn model(&self) -> &Arc<BurnupModel> {
        &self.model
    }
    pub fn inventories(&self) -> &[BurnupInventory] {
        &self.inventories
    }
    pub fn volumes_m3(&self) -> &[f64] {
        &self.volumes_m3
    }
    pub fn transport(&self) -> &NuclideTransport {
        &self.transport
    }
    pub fn prepared(&self) -> Option<&Arc<PreparedBurnupTransport>> {
        self.prepared.as_ref()
    }
    pub fn cumulative_transfer_counts(&self) -> Vec<Vec<f64>> {
        self.transfers
            .iter()
            .map(|row| row.iter().map(|value| value.value()).collect())
            .collect()
    }
    pub fn cumulative_outflow_counts(&self) -> Vec<Vec<f64>> {
        self.outflows
            .iter()
            .map(|row| row.iter().map(|value| value.value()).collect())
            .collect()
    }

    /// Admit a proposed rate/source history without preparing or applying a map.
    /// Flow endpoints/order are fixed; rates may change independently of feed.
    pub fn validate_inputs(
        &self,
        transport: &NuclideTransport,
        production_rates_s: &[Vec<f64>],
        limits: BurnupTransportLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), BurnupError> {
        self.admit(transport, limits, &mut checkpoint)?;
        self.admit_sources(production_rates_s, &mut checkpoint)
    }

    /// Independent candidate using current transport and zero external supply.
    pub fn trial(
        &self,
        duration_s: f64,
        fluxes: &[GroupFlux],
        temperatures_k: &[f64],
        budget: BurnupTransportBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, BurnupTransportStep), BurnupError> {
        self.trial_impl(
            duration_s,
            fluxes,
            temperatures_k,
            &self.transport,
            None,
            budget,
            checkpoint,
        )
    }

    /// Evolve all reactions, movement and constant supply on the same interval.
    /// Supply needs one full row per compartment, including explicit zero rows.
    /// Callers split at changes in flows, sources, fluxes or material conditions.
    pub fn trial_with_inputs(
        &self,
        duration_s: f64,
        fluxes: &[GroupFlux],
        temperatures_k: &[f64],
        transport: &NuclideTransport,
        production_rates_s: &[Vec<f64>],
        budget: BurnupTransportBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, BurnupTransportStep), BurnupError> {
        self.trial_impl(
            duration_s,
            fluxes,
            temperatures_k,
            transport,
            Some(production_rates_s),
            budget,
            checkpoint,
        )
    }

    /// Atomic publication of a homogeneous finite-inventory transport step.
    pub fn advance(
        &mut self,
        duration_s: f64,
        fluxes: &[GroupFlux],
        temperatures_k: &[f64],
        budget: BurnupTransportBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupTransportStep, BurnupError> {
        let (candidate, step) =
            self.trial(duration_s, fluxes, temperatures_k, budget, checkpoint)?;
        *self = candidate;
        Ok(step)
    }

    /// Atomic publication including changed transport rates and external supply.
    pub fn advance_with_inputs(
        &mut self,
        duration_s: f64,
        fluxes: &[GroupFlux],
        temperatures_k: &[f64],
        transport: &NuclideTransport,
        production_rates_s: &[Vec<f64>],
        budget: BurnupTransportBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupTransportStep, BurnupError> {
        let (candidate, step) = self.trial_with_inputs(
            duration_s,
            fluxes,
            temperatures_k,
            transport,
            production_rates_s,
            budget,
            checkpoint,
        )?;
        *self = candidate;
        Ok(step)
    }

    pub(crate) fn admit(
        &self,
        transport: &NuclideTransport,
        limits: BurnupTransportLimits,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), BurnupError> {
        checkpoint()?;
        let cells = self.inventories.len();
        let n = self.model.nuclides().len();
        let nodes = cells
            .checked_mul(n)
            .ok_or(BurnupError::Budget("joint inventory nodes"))?;
        if cells == 0
            || cells > limits.max_compartments
            || nodes > limits.burnup.propagation.maximum_nodes
            || nodes > MAX_DECAY_NODES
        {
            return Err(BurnupError::Budget("joint compartment/nuclide nodes"));
        }
        self.model.admit_limits(limits.burnup)?;
        if self.volumes_m3.len() != cells || !transport.same_topology(&self.transport) {
            return Err(BurnupError::Invalid(
                "fixed compartment volumes and transfer path ordering required",
            ));
        }
        if transport.transfers.len() > limits.max_transfers
            || transport.outflows.len() > limits.max_outflows
        {
            return Err(BurnupError::Budget("nuclide transport paths"));
        }
        let flow_count = transport
            .transfers
            .len()
            .checked_add(transport.outflows.len())
            .ok_or(BurnupError::Budget("transport coefficient count"))?;
        let coefficients = cells
            .checked_mul(self.model.reactions().len())
            .and_then(|count| {
                flow_count
                    .checked_mul(n)
                    .and_then(|flows| count.checked_add(flows))
            })
            .ok_or(BurnupError::Budget("transport coefficient count"))?;
        if coefficients > limits.max_rate_coefficients {
            return Err(BurnupError::Budget("joint reaction/transport coefficients"));
        }
        let edges = cells
            .checked_mul(self.model.products)
            .and_then(|count| {
                transport
                    .transfers
                    .len()
                    .checked_mul(n)
                    .and_then(|flows| count.checked_add(flows))
            })
            .ok_or(BurnupError::Budget("joint reaction/transport edges"))?;
        if edges > limits.burnup.propagation.maximum_edges || edges > MAX_DECAY_EDGES {
            return Err(BurnupError::Budget("joint reaction/transport edges"));
        }
        for (inventory, &volume) in self.inventories.iter().zip(self.volumes_m3.iter()) {
            checkpoint()?;
            if !Arc::ptr_eq(inventory.model(), &self.model) || !positive(volume) {
                return Err(BurnupError::Invalid(
                    "all compartments require the same exact burnup model and positive volumes",
                ));
            }
        }
        for transfer in &transport.transfers {
            checkpoint()?;
            if transfer.from_compartment >= cells
                || transfer.to_compartment >= cells
                || transfer.from_compartment == transfer.to_compartment
            {
                return Err(BurnupError::Invalid("invalid directed nuclide transfer"));
            }
            admit_rate_row(&transfer.rates_per_s, n, checkpoint)?;
        }
        for outflow in &transport.outflows {
            checkpoint()?;
            if outflow.compartment >= cells {
                return Err(BurnupError::Invalid("invalid nuclide outflow"));
            }
            admit_rate_row(&outflow.rates_per_s, n, checkpoint)?;
        }
        checkpoint()?;
        Ok(())
    }

    fn admit_sources(
        &self,
        production_rates_s: &[Vec<f64>],
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), BurnupError> {
        if production_rates_s.len() != self.inventories.len() {
            return Err(BurnupError::Invalid(
                "one production-rate row per compartment required",
            ));
        }
        for row in production_rates_s {
            admit_rate_row(row, self.model.nuclides().len(), checkpoint)?;
        }
        Ok(())
    }

    fn trial_impl(
        &self,
        duration_s: f64,
        fluxes: &[GroupFlux],
        temperatures_k: &[f64],
        transport: &NuclideTransport,
        production_rates_s: Option<&[Vec<f64>]>,
        budget: BurnupTransportBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, BurnupTransportStep), BurnupError> {
        checkpoint()?;
        let local_budget = budget.local()?;
        self.admit(transport, budget.limits, &mut checkpoint)?;
        let cells = self.inventories.len();
        let n = self.model.nuclides().len();
        if !nonnegative(duration_s) || fluxes.len() != cells || temperatures_k.len() != cells {
            return Err(BurnupError::Invalid(
                "finite nonnegative duration and one flux/temperature per compartment required",
            ));
        }
        if let Some(source) = production_rates_s {
            self.admit_sources(source, &mut checkpoint)?;
        }
        let with_sources =
            production_rates_s.is_some_and(|rows| rows.iter().flatten().any(|&rate| rate > 0.0));
        let mut rates = Vec::with_capacity(cells);
        for (flux, &temperature) in fluxes.iter().zip(temperatures_k) {
            rates.push(self.model.rate_constants(
                flux,
                temperature,
                budget.limits.burnup,
                &mut checkpoint,
            )?);
        }
        let prepared = if let Some(prepared) = self.prepared.as_ref().filter(|prepared| {
            Arc::ptr_eq(&prepared.model, &self.model)
                && prepared.duration_s().to_bits() == duration_s.to_bits()
                && prepared.limits == budget.limits
                && prepared.rates == rates
                && *prepared.transport == *transport
                && (!with_sources || prepared.supports_sources())
        }) {
            Arc::clone(prepared)
        } else {
            Arc::new(self.prepare(
                duration_s,
                rates,
                transport,
                budget.limits,
                with_sources,
                &mut checkpoint,
            )?)
        };
        let mut initial = Vec::with_capacity(cells * n);
        let mut supplied = Vec::with_capacity(cells * n);
        let mut produced_nuclei = Vec::with_capacity(cells);
        let mut produced_nuclei_by_compartment = Vec::with_capacity(cells);
        let mut defects = vec![Accumulator::new(0.0); cells];
        for (cell, inventory) in self.inventories.iter().enumerate() {
            let mut row = Vec::with_capacity(n);
            let mut total = Accumulator::new(0.0);
            for (nuclide, &population) in inventory.populations().iter().enumerate() {
                checkpoint()?;
                initial.push(population);
                let rate = production_rates_s.map_or(0.0, |rows| rows[cell][nuclide]);
                let count = product(rate, duration_s)?;
                charge(&mut defects[cell], rate.mul_add(duration_s, -count).abs())?;
                charge(&mut defects[cell], add(&mut total, count)?.abs())?;
                row.push(count);
                supplied.push(count);
            }
            charge(&mut defects[cell], projection_error(total)?)?;
            produced_nuclei_by_compartment.push(checked(total.value())?);
            produced_nuclei.push(row);
        }
        let (evolved, residence) = prepared.plan.apply(&initial, &supplied, &mut checkpoint)?;
        let mut next = self.clone();
        let mut imports = vec![vec![Accumulator::new(0.0); n]; cells];
        let mut exports = vec![vec![Accumulator::new(0.0); n]; cells];
        let mut transfer_counts = Vec::with_capacity(transport.transfers.len());
        let mut outflow_counts = Vec::with_capacity(transport.outflows.len());
        for (path, transfer) in transport.transfers.iter().enumerate() {
            let from = transfer.from_compartment;
            let to = transfer.to_compartment;
            let mut counts = Vec::with_capacity(n);
            for (nuclide, &rate) in transfer.rates_per_s.iter().enumerate() {
                checkpoint()?;
                let integral = residence[from * n + nuclide];
                let count = product(rate, integral)?;
                let error = rate.mul_add(integral, -count).abs();
                charge(&mut defects[from], error)?;
                charge(&mut defects[to], error)?;
                charge(
                    &mut defects[from],
                    add(&mut exports[from][nuclide], count)?.abs(),
                )?;
                charge(
                    &mut defects[to],
                    add(&mut imports[to][nuclide], count)?.abs(),
                )?;
                charge(
                    &mut defects[from],
                    add(&mut next.transfers[path][nuclide], count)?.abs(),
                )?;
                charge(
                    &mut defects[from],
                    projection_error(next.transfers[path][nuclide])?,
                )?;
                counts.push(count);
            }
            transfer_counts.push(counts);
        }
        for (path, outflow) in transport.outflows.iter().enumerate() {
            let cell = outflow.compartment;
            let mut counts = Vec::with_capacity(n);
            for (nuclide, &rate) in outflow.rates_per_s.iter().enumerate() {
                checkpoint()?;
                let integral = residence[cell * n + nuclide];
                let count = product(rate, integral)?;
                charge(&mut defects[cell], rate.mul_add(integral, -count).abs())?;
                charge(
                    &mut defects[cell],
                    add(&mut exports[cell][nuclide], count)?.abs(),
                )?;
                charge(
                    &mut defects[cell],
                    add(&mut next.outflows[path][nuclide], count)?.abs(),
                )?;
                charge(
                    &mut defects[cell],
                    projection_error(next.outflows[path][nuclide])?,
                )?;
                counts.push(count);
            }
            outflow_counts.push(counts);
        }
        let mut burnup = Vec::with_capacity(cells);
        for cell in 0..cells {
            checkpoint()?;
            let imported = export_row(&imports[cell], &mut defects[cell])?;
            let exported = export_row(&exports[cell], &mut defects[cell])?;
            let mut removals = vec![Accumulator::new(0.0); n];
            let mut counts = Vec::with_capacity(self.model.reactions().len());
            let mut count_errors = Vec::with_capacity(self.model.reactions().len());
            for (reaction, &rate) in self.model.reactions().iter().zip(&prepared.rates[cell]) {
                checkpoint()?;
                let integral = residence[cell * n + reaction.parent];
                let count = product(rate, integral)?;
                counts.push(count);
                count_errors.push(rate.mul_add(integral, -count).abs());
                charge(
                    &mut defects[cell],
                    add(&mut removals[reaction.parent], count)?.abs(),
                )?;
            }
            let physical_removals = export_row(&removals, &mut defects[cell])?;
            for nuclide in 0..n {
                // Independent numerical removal counts include BOTH nuclear
                // events and transport departures. Check their attribution.
                let residual = checked(sum(&[
                    evolved.removals[cell * n + nuclide],
                    -removals[nuclide].hi,
                    -removals[nuclide].lo,
                    -exports[cell][nuclide].hi,
                    -exports[cell][nuclide].lo,
                ]))?;
                charge(&mut defects[cell], residual.abs())?;
            }
            let (mut inventory, step) = self.inventories[cell].account_events(
                evolved.populations[cell * n..(cell + 1) * n].to_vec(),
                physical_removals,
                duration_s,
                &[],
                Some(TransportAccounts {
                    reaction_counts: &counts,
                    reaction_roundoff: &count_errors,
                    imported: &imported,
                    exported: &exported,
                }),
                Some(&produced_nuclei[cell]),
                nonnegative_bound(defects[cell])?,
                local_budget,
                &mut checkpoint,
            )?;
            // The accepted cohort owns one map; retain no redundant local maps.
            inventory.prepared = None;
            inventory.production_prepared = None;
            next.inventories[cell] = inventory;
            burnup.push(step);
        }
        let (cohort_balance_residuals, collective_defect) = self.cohort_balance(
            &next,
            &burnup,
            &produced_nuclei,
            &outflow_counts,
            transport,
            &mut checkpoint,
        )?;
        let mut first_defect = Accumulator::new(burnup[0].population_roundoff_nuclei);
        charge(&mut first_defect, collective_defect)?;
        burnup[0].population_roundoff_nuclei = nonnegative_bound(first_defect)?;
        let population_roundoff_nuclei =
            defect_sum(burnup.iter().map(|step| step.population_roundoff_nuclei))?;
        let energy_roundoff_j = defect_sum(burnup.iter().map(|step| step.energy_roundoff_j))?;
        if population_roundoff_nuclei > budget.nuclei_tolerance {
            return Err(BurnupError::Budget(
                "joint transport/population accounting roundoff",
            ));
        }
        if energy_roundoff_j > budget.energy_tolerance_j {
            return Err(BurnupError::Budget(
                "joint reaction energy accounting roundoff",
            ));
        }
        next.transport = Arc::clone(&prepared.transport);
        next.prepared = Some(prepared);
        checkpoint()?;
        Ok((
            next,
            BurnupTransportStep {
                duration_s,
                burnup,
                integrated_populations_s: residence.chunks_exact(n).map(<[f64]>::to_vec).collect(),
                produced_nuclei,
                produced_nuclei_by_compartment,
                transfer_counts,
                outflow_counts,
                cohort_balance_residuals,
                population_roundoff_nuclei,
                energy_roundoff_j,
            },
        ))
    }

    fn prepare(
        &self,
        duration_s: f64,
        rates: Vec<Vec<f64>>,
        transport: &NuclideTransport,
        limits: BurnupTransportLimits,
        with_sources: bool,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<PreparedBurnupTransport, BurnupError> {
        let n = self.model.nuclides().len();
        let cells = self.inventories.len();
        let mut channels = Vec::with_capacity(cells);
        let mut totals = vec![Accumulator::new(0.0); cells * n];
        for (cell, rates) in rates.iter().enumerate() {
            let row = self.model.prepare_channels(rates, &mut *checkpoint)?;
            for (nuclide, &rate) in row.removal_rates.iter().enumerate() {
                totals[cell * n + nuclide] = Accumulator::new(rate);
            }
            channels.push(row);
        }
        for transfer in &transport.transfers {
            for (nuclide, &rate) in transfer.rates_per_s.iter().enumerate() {
                checkpoint()?;
                add(&mut totals[transfer.from_compartment * n + nuclide], rate)?;
            }
        }
        for outflow in &transport.outflows {
            for (nuclide, &rate) in outflow.rates_per_s.iter().enumerate() {
                checkpoint()?;
                add(&mut totals[outflow.compartment * n + nuclide], rate)?;
            }
        }
        let removals: Vec<_> = totals
            .iter()
            .map(|value| checked(value.value()))
            .collect::<Result<_, _>>()?;
        let mut edges = Vec::new();
        for (cell, row) in channels.iter().enumerate() {
            for edge in &row.edges {
                checkpoint()?;
                edges.push(ReactionEdge {
                    parent: cell * n + edge.parent,
                    daughter: cell * n + edge.daughter,
                    yield_per_removal: scale_yield(
                        edge.yield_per_removal,
                        row.removal_rates[edge.parent],
                        removals[cell * n + edge.parent],
                    )?,
                });
            }
        }
        for transfer in &transport.transfers {
            for (nuclide, &rate) in transfer.rates_per_s.iter().enumerate() {
                checkpoint()?;
                if rate == 0.0 {
                    continue;
                }
                let parent = transfer.from_compartment * n + nuclide;
                let weight = rate / removals[parent];
                if !positive(weight) {
                    return Err(BurnupError::Arithmetic(
                        "positive transfer weight is unrepresentable",
                    ));
                }
                edges.push(ReactionEdge {
                    parent,
                    daughter: transfer.to_compartment * n + nuclide,
                    yield_per_removal: weight,
                });
            }
        }
        let plan = JointPlan::prepare(
            &removals,
            &edges,
            duration_s,
            limits,
            with_sources,
            checkpoint,
        )?;
        Ok(PreparedBurnupTransport {
            model: Arc::clone(&self.model),
            transport: Arc::new(transport.clone()),
            rates,
            limits,
            plan,
        })
    }

    fn cohort_balance(
        &self,
        next: &Self,
        reports: &[BurnupStep],
        supplied: &[Vec<f64>],
        outflows: &[Vec<f64>],
        transport: &NuclideTransport,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Vec<f64>, f64), BurnupError> {
        let n = self.model.nuclides().len();
        let mut balances = vec![Accumulator::new(0.0); n];
        let mut error = Accumulator::new(0.0);
        for (cell, report) in reports.iter().enumerate() {
            for nuclide in 0..n {
                checkpoint()?;
                for value in [
                    next.inventories[cell].populations()[nuclide],
                    -self.inventories[cell].populations()[nuclide],
                    report.parent_removals[nuclide],
                    -supplied[cell][nuclide],
                ] {
                    charge(&mut error, add(&mut balances[nuclide], value)?.abs())?;
                }
            }
            for (reaction, &count) in self.model.reactions().iter().zip(&report.reaction_counts) {
                for yield_ in &reaction.products {
                    checkpoint()?;
                    let value = product(count, yield_.yield_per_event)?;
                    charge(
                        &mut error,
                        add(&mut balances[yield_.nuclide], -value)?.abs(),
                    )?;
                }
            }
        }
        for (path, counts) in transport.outflows.iter().zip(outflows) {
            checkpoint()?;
            let _ = path; // Path identity was admitted; all external counts close this cohort.
            for (balance, &count) in balances.iter_mut().zip(counts) {
                charge(&mut error, add(balance, count)?.abs())?;
            }
        }
        let mut residuals = Vec::with_capacity(n);
        for balance in balances {
            let value = checked(balance.value())?;
            charge(&mut error, value.abs())?;
            charge(&mut error, projection_error(balance)?)?;
            residuals.push(value);
        }
        Ok((residuals, nonnegative_bound(error)?))
    }
}

fn admit_rate_row(
    rates: &[f64],
    n: usize,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), BurnupError> {
    if rates.len() != n {
        return Err(BurnupError::Invalid(
            "one nonnegative rate per original nuclide required",
        ));
    }
    for &rate in rates {
        checkpoint()?;
        if !nonnegative(rate) {
            return Err(BurnupError::Invalid(
                "finite nonnegative transport/production rates required",
            ));
        }
    }
    Ok(())
}

fn charge(total: &mut Accumulator, value: f64) -> Result<(), BurnupError> {
    add(total, value)?;
    Ok(())
}

fn export_row(row: &[Accumulator], error: &mut Accumulator) -> Result<Vec<f64>, BurnupError> {
    row.iter()
        .map(|&value| {
            charge(error, projection_error(value)?)?;
            checked(value.value())
        })
        .collect()
}

fn defect_sum(values: impl Iterator<Item = f64>) -> Result<f64, BurnupError> {
    let mut total = Accumulator::new(0.0);
    for value in values {
        charge(&mut total, value)?;
    }
    nonnegative_bound(total)
}

// Choose an ordering that represents the positive normalized yield without
// premature product overflow or ratio underflow. All coefficients remain finite.
fn scale_yield(yield_: f64, physical_rate: f64, total_rate: f64) -> Result<f64, BurnupError> {
    for value in [
        yield_ * (physical_rate / total_rate),
        (yield_ / total_rate) * physical_rate,
        (yield_ * physical_rate) / total_rate,
    ] {
        if positive(value) {
            return Ok(value);
        }
    }
    Err(BurnupError::Arithmetic(
        "positive joint reaction yield is unrepresentable",
    ))
}

fn owner_call<T>(
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    operation: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, fs_la::decay::DecayError>,
) -> Result<T, BurnupError> {
    let mut failure = None;
    let result = operation(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            failure = Some(error);
            true
        }
    });
    if let Some(error) = failure {
        return Err(error.into());
    }
    result.map_err(BurnupError::from)
}

#[cfg(test)]
mod tests;
