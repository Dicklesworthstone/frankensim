//! Bounded finite-volume multigroup neutron diffusion and delayed release.
//!
//! State entries are extensive expected neutron and precursor counts. In cell
//! `c`, group flux is `phi[c,g] = v[g] N[c,g] / V[c]`. A reaction removes a
//! neutron at rate `Sigma[c,g] v[g]`; explicit fission multiplicities produce
//! prompt neutrons and precursor families. An interface carries the net current
//! `K[g] (phi[a,g] - phi[b,g])`, implemented as two conservative positive
//! transfers. `K` has units m², for example `A / (d_a/D_a + d_b/D_b)` for a
//! planar face with diffusion coefficients D [m] and centre-to-face lengths d.
//!
//! The original `fs_la::decay::reaction` owner propagates the frozen positive
//! operator and its population integrals. Thus counts and fluences include
//! particles born during the interval, even in a group with zero removal.
//! Preparation uses bounded dense maps, at most the numerical owner's 64 nodes
//! including precursors, and a state reuses one immutable interval preparation.
//! This is a small spatial diffusion model, not a scalable angular transport,
//! criticality, resonance-processing or evaluated-data-validation solver.
//! An optional prescribed flow transports delayed precursors conservatively;
//! continuous external sources retain their actual injection ages. Neither
//! option evolves bulk fluid, material heat or whole-fuel circulation.
//!
//! The conventional multigroup/precursor equations and diffusion limitations
//! are described in <https://arfc.github.io/moltres/getting_started/theory.html>.
//! G1 manufactured solutions and G3 conservation/refinement checks cover the
//! admitted discrete equations, not their fidelity for a particular reactor.

use std::sync::Arc;

use fs_la::decay::{
    DecayLimits, MAX_DECAY_EDGES, MAX_DECAY_NODES,
    reaction::{ReactionEdge, ReactionPlan},
};

use crate::irradiation::{GroupFlux, NeutronGroups};
use crate::{Accumulator, FissionError, nonnegative, positive, sum};

/// Composition-dependent neutron/depletion coupling on this same geometry.
pub mod burnup;
/// Uniform external injection through the same integrated reaction history.
pub mod source;
/// Explicit delayed-precursor material flow and outflow.
pub mod transport;

pub use transport::{PrecursorOutflow, PrecursorTransfer, PrecursorTransport};

use source::{PreparedPlan, SourceRates};

/// An explicitly identified precursor family. Each decay emits one delayed
/// neutron distributed over the supplied spectrum in its current cell.
#[derive(Clone, Debug, PartialEq)]
pub struct DelayedGroup {
    pub data_id: String,
    pub decay_constant_s_inv: f64,
    /// Nonnegative group fractions summing to one within 64 machine epsilons.
    /// Admission never changes the supplied weights.
    pub spectrum: Vec<f64>,
}

/// Reaction products are expected multiplicities, not branch probabilities.
#[derive(Clone, Debug, PartialEq)]
pub enum DiffusionReactionKind {
    Capture,
    Scatter {
        outgoing_group: usize,
    },
    Fission {
        /// Prompt neutrons per fission in each output group.
        prompt_yields: Vec<f64>,
        /// New precursors per fission in each delayed family.
        delayed_yields: Vec<f64>,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionReaction {
    pub data_id: String,
    pub cell: usize,
    pub incident_group: usize,
    pub macroscopic_cross_section_m_inv: f64,
    pub kind: DiffusionReactionKind,
}

/// A reciprocal interface. Repeated cell pairs add separate authored faces.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionLink {
    pub cell_a: usize,
    pub cell_b: usize,
    pub conductance_m2: Vec<f64>,
}

/// Explicit removal to an exterior with zero imposed incoming flux. Omitted
/// boundary faces are reflecting; no vacuum extrapolation length is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct LeakageBoundary {
    pub cell: usize,
    pub conductance_m2: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionData {
    pub data_id: String,
    pub groups: Arc<NeutronGroups>,
    pub speeds_m_s: Vec<f64>,
    pub cell_volumes_m3: Vec<f64>,
    pub delayed_groups: Vec<DelayedGroup>,
    pub reactions: Vec<DiffusionReaction>,
    pub links: Vec<DiffusionLink>,
    pub leakage: Vec<LeakageBoundary>,
}

/// Data bounds and the original propagator's numerical/work admission limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffusionLimits {
    pub max_cells: usize,
    pub max_neutron_groups: usize,
    pub max_delayed_groups: usize,
    pub max_reactions: usize,
    /// Neutron interfaces plus directed precursor transfer paths.
    pub max_links: usize,
    /// Neutron leakage faces plus precursor outflow paths.
    pub max_leakage_boundaries: usize,
    /// Authored speeds, volumes, decay constants, spectra, cross sections,
    /// fission multiplicities, conductances and precursor flow rates combined.
    pub max_coefficients: usize,
    pub propagation: DecayLimits,
}

impl Default for DiffusionLimits {
    fn default() -> Self {
        Self {
            max_cells: MAX_DECAY_NODES,
            max_neutron_groups: MAX_DECAY_NODES,
            max_delayed_groups: MAX_DECAY_NODES,
            max_reactions: MAX_DECAY_EDGES,
            max_links: MAX_DECAY_EDGES,
            max_leakage_boundaries: MAX_DECAY_EDGES,
            max_coefficients: 65_536,
            propagation: DecayLimits::default(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffusionBudget {
    pub limits: DiffusionLimits,
    /// Absolute sum of node/count attribution defects and measured account
    /// rounding [expected particles]. Not a certified propagator error bound.
    pub population_tolerance: f64,
}

/// Immutable admitted geometry, data and first-order reaction operator.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionModel {
    data: DiffusionData,
    precursor_transport: PrecursorTransport,
    neutron_nodes: usize,
    nodes: usize,
    coefficients: usize,
    flux_per_neutron: Vec<f64>,
    reaction_rates: Vec<f64>,
    link_rates: Vec<[f64; 2]>,
    leakage_rates: Vec<f64>,
    precursor_transfer_rates: Vec<f64>,
    precursor_outflow_rates: Vec<f64>,
    removal_rates: Vec<f64>,
    edges: Vec<ReactionEdge>,
}

impl DiffusionModel {
    pub fn new(data: DiffusionData, limits: DiffusionLimits) -> Result<Self, DiffusionError> {
        Self::new_with_checkpoint(data, limits, || Ok(()))
    }

    /// Admit a complete operator without publishing partially validated data.
    pub fn new_with_checkpoint(
        data: DiffusionData,
        limits: DiffusionLimits,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DiffusionError> {
        Self::new_with_transport_and_checkpoint(
            data,
            PrecursorTransport::default(),
            limits,
            checkpoint,
        )
    }

    /// Admit prescribed delayed-precursor flow alongside the original data.
    pub fn new_with_transport(
        data: DiffusionData,
        transport: PrecursorTransport,
        limits: DiffusionLimits,
    ) -> Result<Self, DiffusionError> {
        Self::new_with_transport_and_checkpoint(data, transport, limits, || Ok(()))
    }

    pub fn new_with_transport_and_checkpoint(
        data: DiffusionData,
        transport: PrecursorTransport,
        limits: DiffusionLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DiffusionError> {
        checkpoint()?;
        let cells = data.cell_volumes_m3.len();
        let groups = data.groups.len();
        let families = data.delayed_groups.len();
        shape_limits(&data, limits)?;
        transport.check_limits(&data, limits)?;
        if data.data_id.trim().is_empty() || cells == 0 || data.speeds_m_s.len() != groups {
            return Err(DiffusionError::Invalid(
                "identified cells, positive volumes and one positive speed per group required",
            ));
        }
        let neutron_nodes = count_product(cells, groups)?;
        let nodes = neutron_nodes
            .checked_add(count_product(cells, families)?)
            .ok_or(DiffusionError::Budget("diffusion nodes"))?;
        node_limit(nodes, limits)?;
        if !data.speeds_m_s.iter().copied().all(positive)
            || !data.cell_volumes_m3.iter().copied().all(positive)
        {
            return Err(DiffusionError::Invalid(
                "positive finite group speeds and cell volumes required",
            ));
        }
        let mut coefficients = groups
            .checked_add(cells)
            .ok_or(DiffusionError::Budget("diffusion coefficients"))?;
        coefficient_limit(coefficients, limits)?;
        for family in &data.delayed_groups {
            checkpoint()?;
            coefficients = add_coefficients(coefficients, family.spectrum.len(), limits)?;
            coefficients = add_coefficients(coefficients, 1, limits)?;
            if family.data_id.trim().is_empty()
                || !positive(family.decay_constant_s_inv)
                || family.spectrum.len() != groups
                || !family.spectrum.iter().copied().all(nonnegative)
                || !positive(sum(&family.spectrum))
                || (sum(&family.spectrum) - 1.0).abs() > 64.0 * f64::EPSILON
            {
                return Err(DiffusionError::Invalid(
                    "identified delayed families require positive decay and a normalized spectrum",
                ));
            }
        }
        for reaction in &data.reactions {
            checkpoint()?;
            coefficients = add_coefficients(coefficients, 1, limits)?;
            if reaction.data_id.trim().is_empty()
                || reaction.cell >= cells
                || reaction.incident_group >= groups
                || !nonnegative(reaction.macroscopic_cross_section_m_inv)
            {
                return Err(DiffusionError::Invalid("invalid diffusion reaction"));
            }
            match &reaction.kind {
                DiffusionReactionKind::Capture => {}
                DiffusionReactionKind::Scatter { outgoing_group } => {
                    if *outgoing_group >= groups {
                        return Err(DiffusionError::Invalid("invalid scattering destination"));
                    }
                }
                DiffusionReactionKind::Fission {
                    prompt_yields,
                    delayed_yields,
                } => {
                    coefficients = add_coefficients(coefficients, prompt_yields.len(), limits)?;
                    coefficients = add_coefficients(coefficients, delayed_yields.len(), limits)?;
                    if prompt_yields.len() != groups
                        || delayed_yields.len() != families
                        || !prompt_yields
                            .iter()
                            .chain(delayed_yields)
                            .copied()
                            .all(nonnegative)
                    {
                        return Err(DiffusionError::Invalid(
                            "fission requires explicit nonnegative yields per group and family",
                        ));
                    }
                }
            }
        }
        for link in &data.links {
            checkpoint()?;
            coefficients = add_coefficients(coefficients, link.conductance_m2.len(), limits)?;
            if link.cell_a >= cells
                || link.cell_b >= cells
                || link.cell_a == link.cell_b
                || link.conductance_m2.len() != groups
                || !link.conductance_m2.iter().copied().all(nonnegative)
            {
                return Err(DiffusionError::Invalid("invalid diffusion interface"));
            }
        }
        for boundary in &data.leakage {
            checkpoint()?;
            coefficients = add_coefficients(coefficients, boundary.conductance_m2.len(), limits)?;
            if boundary.cell >= cells
                || boundary.conductance_m2.len() != groups
                || !boundary.conductance_m2.iter().copied().all(nonnegative)
            {
                return Err(DiffusionError::Invalid("invalid leakage boundary"));
            }
        }

        coefficients = add_coefficients(coefficients, transport.transfers.len(), limits)?;
        coefficients = add_coefficients(coefficients, transport.outflows.len(), limits)?;
        let transport_rates = transport.rates(&data, &mut checkpoint)?;

        let mut removals = filled(nodes, Accumulator::new(0.0))?;
        let mut flux_per_neutron = reserved(neutron_nodes)?;
        for &volume in &data.cell_volumes_m3 {
            checkpoint()?;
            for &speed in &data.speeds_m_s {
                flux_per_neutron.push(quotient(speed, volume)?);
            }
        }
        let mut reaction_rates = reserved(data.reactions.len())?;
        for reaction in &data.reactions {
            checkpoint()?;
            let rate = product(
                reaction.macroscopic_cross_section_m_inv,
                data.speeds_m_s[reaction.incident_group],
            )?;
            removals[reaction.cell * groups + reaction.incident_group].add(rate)?;
            reaction_rates.push(rate);
        }
        let mut link_rates = reserved(count_product(data.links.len(), groups)?)?;
        for link in &data.links {
            checkpoint()?;
            for (group, &conductance) in link.conductance_m2.iter().enumerate() {
                let a = link.cell_a * groups + group;
                let b = link.cell_b * groups + group;
                let rates = [
                    product(conductance, flux_per_neutron[a])?,
                    product(conductance, flux_per_neutron[b])?,
                ];
                removals[a].add(rates[0])?;
                removals[b].add(rates[1])?;
                link_rates.push(rates);
            }
        }
        let mut leakage_rates = reserved(count_product(data.leakage.len(), groups)?)?;
        for boundary in &data.leakage {
            checkpoint()?;
            for (group, &conductance) in boundary.conductance_m2.iter().enumerate() {
                let node = boundary.cell * groups + group;
                let rate = product(conductance, flux_per_neutron[node])?;
                removals[node].add(rate)?;
                leakage_rates.push(rate);
            }
        }
        for cell in 0..cells {
            checkpoint()?;
            for (family, data) in data.delayed_groups.iter().enumerate() {
                removals[neutron_nodes + cell * families + family]
                    .add(data.decay_constant_s_inv)?;
            }
        }
        for (transfer, &rate) in transport.transfers.iter().zip(&transport_rates.transfers) {
            checkpoint()?;
            for family in 0..families {
                removals[neutron_nodes + transfer.from_cell * families + family].add(rate)?;
            }
        }
        for (outflow, &rate) in transport.outflows.iter().zip(&transport_rates.outflows) {
            checkpoint()?;
            for family in 0..families {
                removals[neutron_nodes + outflow.cell * families + family].add(rate)?;
            }
        }
        let removal_rates: Vec<_> = removals.into_iter().map(Accumulator::value).collect();
        let mut edges = Vec::new();
        for (reaction, &rate) in data.reactions.iter().zip(&reaction_rates) {
            checkpoint()?;
            let source = reaction.cell * groups + reaction.incident_group;
            let fraction = fraction(rate, removal_rates[source])?;
            match &reaction.kind {
                DiffusionReactionKind::Capture => {}
                DiffusionReactionKind::Scatter { outgoing_group } => push_edge(
                    &mut edges,
                    source,
                    reaction.cell * groups + outgoing_group,
                    fraction,
                    1.0,
                    limits,
                )?,
                DiffusionReactionKind::Fission {
                    prompt_yields,
                    delayed_yields,
                } => {
                    for (group, &yield_) in prompt_yields.iter().enumerate() {
                        push_edge(
                            &mut edges,
                            source,
                            reaction.cell * groups + group,
                            fraction,
                            yield_,
                            limits,
                        )?;
                    }
                    for (family, &yield_) in delayed_yields.iter().enumerate() {
                        push_edge(
                            &mut edges,
                            source,
                            neutron_nodes + reaction.cell * families + family,
                            fraction,
                            yield_,
                            limits,
                        )?;
                    }
                }
            }
        }
        for (link_index, link) in data.links.iter().enumerate() {
            checkpoint()?;
            for group in 0..groups {
                let a = link.cell_a * groups + group;
                let b = link.cell_b * groups + group;
                let rates = link_rates[link_index * groups + group];
                push_edge(
                    &mut edges,
                    a,
                    b,
                    fraction(rates[0], removal_rates[a])?,
                    1.0,
                    limits,
                )?;
                push_edge(
                    &mut edges,
                    b,
                    a,
                    fraction(rates[1], removal_rates[b])?,
                    1.0,
                    limits,
                )?;
            }
        }
        for cell in 0..cells {
            checkpoint()?;
            for (family, data) in data.delayed_groups.iter().enumerate() {
                let source = neutron_nodes + cell * families + family;
                let decay_fraction = fraction(data.decay_constant_s_inv, removal_rates[source])?;
                for (group, &weight) in data.spectrum.iter().enumerate() {
                    push_edge(
                        &mut edges,
                        source,
                        cell * groups + group,
                        decay_fraction,
                        weight,
                        limits,
                    )?;
                }
            }
        }
        for (transfer, &rate) in transport.transfers.iter().zip(&transport_rates.transfers) {
            checkpoint()?;
            for family in 0..families {
                let source = neutron_nodes + transfer.from_cell * families + family;
                let destination = neutron_nodes + transfer.to_cell * families + family;
                push_edge(
                    &mut edges,
                    source,
                    destination,
                    fraction(rate, removal_rates[source])?,
                    1.0,
                    limits,
                )?;
            }
        }
        // Original-owner admission checks every numerical limit and topology;
        // a zero interval avoids performing an unused physical propagation.
        owner_call(&mut checkpoint, |cancelled| {
            ReactionPlan::prepare(&removal_rates, &edges, 0.0, limits.propagation, cancelled)
        })?;
        checkpoint()?;
        Ok(Self {
            data,
            precursor_transport: transport,
            neutron_nodes,
            nodes,
            coefficients,
            flux_per_neutron,
            reaction_rates,
            link_rates,
            leakage_rates,
            precursor_transfer_rates: transport_rates.transfers,
            precursor_outflow_rates: transport_rates.outflows,
            removal_rates,
            edges,
        })
    }

    #[must_use]
    pub fn data(&self) -> &DiffusionData {
        &self.data
    }
    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data.data_id
    }
    #[must_use]
    pub fn groups(&self) -> &Arc<NeutronGroups> {
        &self.data.groups
    }
    #[must_use]
    pub fn speeds_m_s(&self) -> &[f64] {
        &self.data.speeds_m_s
    }
    #[must_use]
    pub fn cell_volumes_m3(&self) -> &[f64] {
        &self.data.cell_volumes_m3
    }
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.data.cell_volumes_m3.len()
    }
    #[must_use]
    pub fn delayed_groups(&self) -> &[DelayedGroup] {
        &self.data.delayed_groups
    }
    #[must_use]
    pub fn reactions(&self) -> &[DiffusionReaction] {
        &self.data.reactions
    }
    #[must_use]
    pub fn links(&self) -> &[DiffusionLink] {
        &self.data.links
    }
    #[must_use]
    pub fn leakage(&self) -> &[LeakageBoundary] {
        &self.data.leakage
    }
    #[must_use]
    pub fn reaction_rate_constants(&self) -> &[f64] {
        &self.reaction_rates
    }

    #[must_use]
    pub fn precursor_transport(&self) -> &PrecursorTransport {
        &self.precursor_transport
    }

    /// Rebuild material coefficients while retaining the admitted precursor
    /// flow. Callers then use `rebind_model` to preserve the inventory only if
    /// geometry and ordered physical identities still agree.
    pub fn rebuild_with_data(
        &self,
        data: DiffusionData,
        limits: DiffusionLimits,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DiffusionError> {
        Self::new_with_transport_and_checkpoint(
            data,
            self.precursor_transport.clone(),
            limits,
            checkpoint,
        )
    }

    fn check_limits(&self, limits: DiffusionLimits) -> Result<(), DiffusionError> {
        shape_limits(&self.data, limits)?;
        self.precursor_transport.check_limits(&self.data, limits)?;
        node_limit(self.nodes, limits)?;
        coefficient_limit(self.coefficients, limits)?;
        if self.edges.len() > limits.propagation.maximum_edges {
            return Err(DiffusionError::Budget("diffusion product edges"));
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct PreparedDiffusion {
    model: Arc<DiffusionModel>,
    limits: DiffusionLimits,
    plan: PreparedPlan,
}

/// A complete accepted spatial inventory. Clones retain independent counts
/// while sharing immutable geometry and an optional unchanged-interval map.
#[derive(Clone, Debug)]
pub struct DiffusionState {
    model: Arc<DiffusionModel>,
    neutron_counts: Vec<f64>,
    precursor_counts: Vec<f64>,
    prepared: Option<Arc<PreparedDiffusion>>,
}

impl DiffusionState {
    pub fn new(
        model: Arc<DiffusionModel>,
        neutron_counts: Vec<f64>,
        precursor_counts: Vec<f64>,
    ) -> Result<Self, DiffusionError> {
        if neutron_counts.len() != model.neutron_nodes
            || precursor_counts.len() != model.nodes - model.neutron_nodes
            || !neutron_counts
                .iter()
                .chain(&precursor_counts)
                .copied()
                .all(nonnegative)
        {
            return Err(DiffusionError::Invalid(
                "one nonnegative count per cell/group and cell/precursor required",
            ));
        }
        let mut total = Accumulator::new(0.0);
        for &count in neutron_counts.iter().chain(&precursor_counts) {
            total.add(count)?;
        }
        checked(total.value())?;
        Ok(Self {
            model,
            neutron_counts,
            precursor_counts,
            prepared: None,
        })
    }

    #[must_use]
    pub fn model(&self) -> &Arc<DiffusionModel> {
        &self.model
    }
    #[must_use]
    pub fn neutron_counts(&self) -> &[f64] {
        &self.neutron_counts
    }
    #[must_use]
    pub fn precursor_counts(&self) -> &[f64] {
        &self.precursor_counts
    }

    /// Instantaneous scalar flux from the extensive neutron inventory.
    pub fn group_flux(&self, cell: usize) -> Result<GroupFlux, DiffusionError> {
        if cell >= self.model.cell_count() {
            return Err(DiffusionError::Invalid("diffusion cell outside geometry"));
        }
        let groups = self.model.groups().len();
        let mut values = reserved(groups)?;
        for group in 0..groups {
            let node = cell * groups + group;
            values.push(product(
                self.neutron_counts[node],
                self.model.flux_per_neutron[node],
            )?);
        }
        Ok(GroupFlux::new(Arc::clone(self.model.groups()), values)?)
    }

    /// Rebind frozen material/transport coefficients without reinterpreting an
    /// accepted inventory. Cell volumes, ordered face endpoints, energy grid,
    /// speeds and complete ordered precursor identities/data remain unchanged.
    pub fn rebind_model(&self, model: Arc<DiffusionModel>) -> Result<Self, DiffusionError> {
        if self.model.groups() != model.groups()
            || self.model.speeds_m_s() != model.speeds_m_s()
            || self.model.cell_volumes_m3() != model.cell_volumes_m3()
            || self.model.delayed_groups() != model.delayed_groups()
            || !self
                .model
                .precursor_transport
                .same_topology(&model.precursor_transport)
            || self.model.links().len() != model.links().len()
            || !self
                .model
                .links()
                .iter()
                .zip(model.links())
                .all(|(a, b)| (a.cell_a, a.cell_b) == (b.cell_a, b.cell_b))
            || self.model.leakage().len() != model.leakage().len()
            || !self
                .model
                .leakage()
                .iter()
                .zip(model.leakage())
                .all(|(a, b)| a.cell == b.cell)
        {
            return Err(DiffusionError::Invalid(
                "coefficient rebinding requires unchanged geometry, grid, speeds and precursors",
            ));
        }
        let prepared = if Arc::ptr_eq(&self.model, &model) {
            self.prepared.clone()
        } else {
            None
        };
        Ok(Self {
            model,
            neutron_counts: self.neutron_counts.clone(),
            precursor_counts: self.precursor_counts.clone(),
            prepared,
        })
    }

    /// Prepare a complete independent candidate for unchanged coefficients.
    /// Neither cancellation nor numerical/accounting refusal changes `self`.
    pub fn trial(
        &self,
        duration_s: f64,
        budget: DiffusionBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, DiffusionStep), DiffusionError> {
        self.trial_from_sources(duration_s, None, budget, checkpoint)
    }

    fn trial_from_sources(
        &self,
        duration_s: f64,
        sources: Option<SourceRates<'_>>,
        budget: DiffusionBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, DiffusionStep), DiffusionError> {
        checkpoint()?;
        if !nonnegative(duration_s) || !nonnegative(budget.population_tolerance) {
            return Err(DiffusionError::Invalid(
                "finite nonnegative duration and population tolerance required",
            ));
        }
        self.model.check_limits(budget.limits)?;
        let mut accounts = Accounts::new(self.model.nodes)?;
        let mut injected = source::injected_counts(
            sources,
            self.model.neutron_nodes,
            self.model.nodes,
            duration_s,
            &mut accounts,
            &mut checkpoint,
        )?;
        let with_sources = injected.iter().any(|&count| count > 0.0);
        let prepared = if let Some(prepared) = self.prepared.as_ref().filter(|prepared| {
            Arc::ptr_eq(&self.model, &prepared.model)
                && prepared.plan.duration_s().to_bits() == duration_s.to_bits()
                && prepared.limits == budget.limits
                && (!with_sources || prepared.plan.supports_sources())
        }) {
            Arc::clone(prepared)
        } else {
            let plan = PreparedPlan::prepare(
                &self.model,
                duration_s,
                budget.limits,
                with_sources,
                &mut checkpoint,
            )?;
            Arc::new(PreparedDiffusion {
                model: Arc::clone(&self.model),
                limits: budget.limits,
                plan,
            })
        };
        let mut initial = reserved(self.model.nodes)?;
        initial.extend_from_slice(&self.neutron_counts);
        initial.extend_from_slice(&self.precursor_counts);
        let evolved = prepared.plan.apply(&initial, &injected, &mut checkpoint)?;
        let integrals = prepared
            .plan
            .integrals(&initial, &injected, &mut checkpoint)?;
        let groups = self.model.groups().len();
        let families = self.model.delayed_groups().len();
        let neutron_nodes = self.model.neutron_nodes;
        let mut reaction_counts = reserved(self.model.reactions().len())?;
        for (reaction, &rate) in self
            .model
            .reactions()
            .iter()
            .zip(&self.model.reaction_rates)
        {
            checkpoint()?;
            let source = reaction.cell * groups + reaction.incident_group;
            let count = accounts.multiply(rate, integrals[source])?;
            accounts.remove(source, count)?;
            reaction_counts.push(count);
            match &reaction.kind {
                DiffusionReactionKind::Capture => accounts.total(CAPTURES, count)?,
                DiffusionReactionKind::Scatter { outgoing_group } => {
                    accounts.total(SCATTERS, count)?;
                    accounts.produce(reaction.cell * groups + outgoing_group, count, 1.0)?;
                }
                DiffusionReactionKind::Fission {
                    prompt_yields,
                    delayed_yields,
                } => {
                    accounts.total(FISSIONS, count)?;
                    for (group, &yield_) in prompt_yields.iter().enumerate() {
                        let born =
                            accounts.produce(reaction.cell * groups + group, count, yield_)?;
                        accounts.total(PROMPT, born)?;
                    }
                    for (family, &yield_) in delayed_yields.iter().enumerate() {
                        let born = accounts.produce(
                            neutron_nodes + reaction.cell * families + family,
                            count,
                            yield_,
                        )?;
                        accounts.total(PRECURSORS, born)?;
                    }
                }
            }
        }
        let mut link_a_to_b_counts = reserved(self.model.link_rates.len())?;
        let mut link_b_to_a_counts = reserved(self.model.link_rates.len())?;
        for (link_index, link) in self.model.links().iter().enumerate() {
            checkpoint()?;
            for group in 0..groups {
                let a = link.cell_a * groups + group;
                let b = link.cell_b * groups + group;
                let rates = self.model.link_rates[link_index * groups + group];
                let forward = accounts.multiply(rates[0], integrals[a])?;
                let backward = accounts.multiply(rates[1], integrals[b])?;
                accounts.remove(a, forward)?;
                accounts.produce(b, forward, 1.0)?;
                accounts.remove(b, backward)?;
                accounts.produce(a, backward, 1.0)?;
                link_a_to_b_counts.push(forward);
                link_b_to_a_counts.push(backward);
            }
        }
        let mut leakage_counts = reserved(self.model.leakage_rates.len())?;
        for (boundary_index, boundary) in self.model.leakage().iter().enumerate() {
            checkpoint()?;
            for group in 0..groups {
                let source = boundary.cell * groups + group;
                let rate = self.model.leakage_rates[boundary_index * groups + group];
                let count = accounts.multiply(rate, integrals[source])?;
                accounts.remove(source, count)?;
                accounts.total(LEAKED, count)?;
                leakage_counts.push(count);
            }
        }
        let mut delayed_decays = reserved(self.precursor_counts.len())?;
        for cell in 0..self.model.cell_count() {
            checkpoint()?;
            for (family, data) in self.model.delayed_groups().iter().enumerate() {
                let source = neutron_nodes + cell * families + family;
                let count = accounts.multiply(data.decay_constant_s_inv, integrals[source])?;
                accounts.remove(source, count)?;
                accounts.total(DECAYS, count)?;
                delayed_decays.push(count);
                for (group, &weight) in data.spectrum.iter().enumerate() {
                    let born = accounts.produce(cell * groups + group, count, weight)?;
                    accounts.total(DELAYED, born)?;
                }
            }
        }

        let mut precursor_transfer_counts = reserved(count_product(
            self.model.precursor_transport.transfers.len(),
            families,
        )?)?;
        for (transfer, &rate) in self
            .model
            .precursor_transport
            .transfers
            .iter()
            .zip(&self.model.precursor_transfer_rates)
        {
            checkpoint()?;
            for family in 0..families {
                let source = neutron_nodes + transfer.from_cell * families + family;
                let destination = neutron_nodes + transfer.to_cell * families + family;
                let count = accounts.multiply(rate, integrals[source])?;
                accounts.remove(source, count)?;
                accounts.produce(destination, count, 1.0)?;
                precursor_transfer_counts.push(count);
            }
        }
        let mut precursor_outflow_counts = reserved(count_product(
            self.model.precursor_transport.outflows.len(),
            families,
        )?)?;
        for (outflow, &rate) in self
            .model
            .precursor_transport
            .outflows
            .iter()
            .zip(&self.model.precursor_outflow_rates)
        {
            checkpoint()?;
            for family in 0..families {
                let source = neutron_nodes + outflow.cell * families + family;
                let count = accounts.multiply(rate, integrals[source])?;
                accounts.remove(source, count)?;
                accounts.total(ESCAPED_PRECURSORS, count)?;
                precursor_outflow_counts.push(count);
            }
        }

        let mut population_balance_residuals = reserved(self.model.nodes)?;
        for node in 0..self.model.nodes {
            checkpoint()?;
            let incoming = accounts.incoming[node];
            let outgoing = accounts.outgoing[node];
            let residual = checked(sum(&[
                evolved.populations[node],
                -initial[node],
                -injected[node],
                outgoing.hi,
                outgoing.lo,
                -incoming.hi,
                -incoming.lo,
            ]))?;
            population_balance_residuals.push(residual);
            accounts.error.add(residual.abs())?;
            // Independently integrated removals must agree with the explicit
            // physical channels obtained from the population-integral map.
            accounts
                .error
                .add(checked(sum(&[evolved.removals[node], -outgoing.hi, -outgoing.lo]))?.abs())?;
        }
        let mut neutron_balance = Accumulator::new(0.0);
        let mut precursor_balance = Accumulator::new(0.0);
        for node in 0..self.model.nodes {
            let account = if node < neutron_nodes {
                &mut neutron_balance
            } else {
                &mut precursor_balance
            };
            account_add(account, evolved.populations[node], &mut accounts.error)?;
            account_add(account, -initial[node], &mut accounts.error)?;
        }
        for (kind, sign) in [
            (FISSIONS, 1.0),
            (CAPTURES, 1.0),
            (LEAKED, 1.0),
            (PROMPT, -1.0),
            (DELAYED, -1.0),
            (SOURCE_NEUTRONS, -1.0),
        ] {
            account_add(
                &mut neutron_balance,
                sign * accounts.totals[kind].hi,
                &mut accounts.error,
            )?;
            account_add(
                &mut neutron_balance,
                sign * accounts.totals[kind].lo,
                &mut accounts.error,
            )?;
        }
        for (kind, sign) in [
            (PRECURSORS, -1.0),
            (DECAYS, 1.0),
            (SOURCE_PRECURSORS, -1.0),
            (ESCAPED_PRECURSORS, 1.0),
        ] {
            account_add(
                &mut precursor_balance,
                sign * accounts.totals[kind].hi,
                &mut accounts.error,
            )?;
            account_add(
                &mut precursor_balance,
                sign * accounts.totals[kind].lo,
                &mut accounts.error,
            )?;
        }
        let neutron_balance_residual = checked(neutron_balance.value())?;
        let precursor_balance_residual = checked(precursor_balance.value())?;
        accounts.error.add(neutron_balance_residual.abs())?;
        accounts.error.add(precursor_balance_residual.abs())?;
        for &account in &accounts.totals {
            accounts
                .error
                .add(checked(sum(&[account.value(), -account.hi, -account.lo]))?.abs())?;
        }
        let population_roundoff = upward_account_bound(accounts.error)?;
        if population_roundoff > budget.population_tolerance {
            return Err(DiffusionError::Budget(
                "diffusion population and event accounting",
            ));
        }

        let mut integrated_fluence_m2 = reserved(neutron_nodes)?;
        let mut mean_fluxes = reserved(self.model.cell_count())?;
        for cell in 0..self.model.cell_count() {
            checkpoint()?;
            let mut mean = reserved(groups)?;
            for group in 0..groups {
                let node = cell * groups + group;
                let fluence = product(self.model.flux_per_neutron[node], integrals[node])?;
                integrated_fluence_m2.push(fluence);
                mean.push(if duration_s == 0.0 {
                    product(self.model.flux_per_neutron[node], initial[node])?
                } else {
                    quotient(fluence, duration_s)?
                });
            }
            mean_fluxes.push(GroupFlux::new(Arc::clone(self.model.groups()), mean)?);
        }
        let mut neutron_counts = evolved.populations;
        let precursor_counts = neutron_counts.split_off(neutron_nodes);
        let injected_precursor_counts = injected.split_off(neutron_nodes);
        let next = Self {
            model: Arc::clone(&self.model),
            neutron_counts,
            precursor_counts,
            prepared: Some(prepared),
        };
        checkpoint()?;
        Ok((
            next,
            DiffusionStep {
                duration_s,
                integrated_fluence_m2,
                mean_fluxes,
                reaction_counts,
                link_a_to_b_counts,
                link_b_to_a_counts,
                leakage_counts,
                delayed_decays,
                injected_neutron_counts: injected,
                injected_precursor_counts,
                precursor_transfer_counts,
                precursor_outflow_counts,
                injected_neutrons: accounts.totals[SOURCE_NEUTRONS].value(),
                injected_precursors: accounts.totals[SOURCE_PRECURSORS].value(),
                escaped_precursors: accounts.totals[ESCAPED_PRECURSORS].value(),
                fissions: accounts.totals[FISSIONS].value(),
                captures: accounts.totals[CAPTURES].value(),
                scatters: accounts.totals[SCATTERS].value(),
                prompt_neutrons: accounts.totals[PROMPT].value(),
                delayed_precursors_produced: accounts.totals[PRECURSORS].value(),
                precursor_decays: accounts.totals[DECAYS].value(),
                delayed_neutrons: accounts.totals[DELAYED].value(),
                leaked_neutrons: accounts.totals[LEAKED].value(),
                population_balance_residuals,
                neutron_balance_residual,
                precursor_balance_residual,
                population_roundoff,
            },
        ))
    }

    /// Commit only a fully accepted candidate.
    pub fn advance(
        &mut self,
        duration_s: f64,
        budget: DiffusionBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DiffusionStep, DiffusionError> {
        let (candidate, step) = self.trial(duration_s, budget, checkpoint)?;
        *self = candidate;
        Ok(step)
    }
}

/// All flat spatial vectors use cell-major order. Reaction counts follow the
/// authored reaction order; interface/boundary vectors use face-major groups.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionStep {
    pub duration_s: f64,
    /// Integral of scalar group flux over the interval [1/m²].
    pub integrated_fluence_m2: Vec<f64>,
    /// Fluence / duration. For zero duration, the initial instantaneous flux.
    pub mean_fluxes: Vec<GroupFlux>,
    pub reaction_counts: Vec<f64>,
    pub link_a_to_b_counts: Vec<f64>,
    pub link_b_to_a_counts: Vec<f64>,
    pub leakage_counts: Vec<f64>,
    /// Actual precursor disappearances, cell-major then family.
    pub delayed_decays: Vec<f64>,
    /// Uniform external injection totals, in the corresponding state order.
    pub injected_neutron_counts: Vec<f64>,
    pub injected_precursor_counts: Vec<f64>,
    /// Directed transfer-major then delayed-family counts.
    pub precursor_transfer_counts: Vec<f64>,
    /// Outflow-major then delayed-family counts; no later exterior decay is inferred.
    pub precursor_outflow_counts: Vec<f64>,
    pub injected_neutrons: f64,
    pub injected_precursors: f64,
    pub escaped_precursors: f64,
    pub fissions: f64,
    pub captures: f64,
    pub scatters: f64,
    pub prompt_neutrons: f64,
    pub delayed_precursors_produced: f64,
    pub precursor_decays: f64,
    pub delayed_neutrons: f64,
    pub leaked_neutrons: f64,
    /// N1 - N0 - injection + removals - products, neutrons then precursors.
    pub population_balance_residuals: Vec<f64>,
    pub neutron_balance_residual: f64,
    pub precursor_balance_residual: f64,
    /// Noncancelling absolute account defects, not a forward-error certificate.
    pub population_roundoff: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DiffusionError {
    Invalid(&'static str),
    Budget(&'static str),
    Arithmetic(&'static str),
    Fission(FissionError),
    Propagation(fs_la::decay::DecayError),
}

impl From<FissionError> for DiffusionError {
    fn from(value: FissionError) -> Self {
        Self::Fission(value)
    }
}
impl From<fs_la::decay::DecayError> for DiffusionError {
    fn from(value: fs_la::decay::DecayError) -> Self {
        Self::Propagation(value)
    }
}
impl core::fmt::Display for DiffusionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial neutron diffusion refused: {self:?}")
    }
}
impl std::error::Error for DiffusionError {}

const FISSIONS: usize = 0;
const CAPTURES: usize = 1;
const SCATTERS: usize = 2;
const PROMPT: usize = 3;
const PRECURSORS: usize = 4;
const DECAYS: usize = 5;
const DELAYED: usize = 6;
const LEAKED: usize = 7;
const SOURCE_NEUTRONS: usize = 8;
const SOURCE_PRECURSORS: usize = 9;
const ESCAPED_PRECURSORS: usize = 10;

struct Accounts {
    incoming: Vec<Accumulator>,
    outgoing: Vec<Accumulator>,
    totals: [Accumulator; 11],
    error: Accumulator,
}

impl Accounts {
    fn new(nodes: usize) -> Result<Self, DiffusionError> {
        Ok(Self {
            incoming: filled(nodes, Accumulator::new(0.0))?,
            outgoing: filled(nodes, Accumulator::new(0.0))?,
            totals: [Accumulator::new(0.0); 11],
            error: Accumulator::new(0.0),
        })
    }
    fn multiply(&mut self, a: f64, b: f64) -> Result<f64, DiffusionError> {
        let value = product(a, b)?;
        self.error.add(a.mul_add(b, -value).abs())?;
        Ok(value)
    }
    fn remove(&mut self, node: usize, count: f64) -> Result<(), DiffusionError> {
        account_add(&mut self.outgoing[node], count, &mut self.error)
    }
    fn produce(&mut self, node: usize, count: f64, yield_: f64) -> Result<f64, DiffusionError> {
        let value = self.multiply(count, yield_)?;
        account_add(&mut self.incoming[node], value, &mut self.error)?;
        Ok(value)
    }
    fn total(&mut self, kind: usize, value: f64) -> Result<(), DiffusionError> {
        account_add(&mut self.totals[kind], value, &mut self.error)
    }
}

fn account_add(
    account: &mut Accumulator,
    value: f64,
    error: &mut Accumulator,
) -> Result<(), DiffusionError> {
    error.add(account.add(value)?.abs())?;
    Ok(())
}

fn shape_limits(data: &DiffusionData, limits: DiffusionLimits) -> Result<(), DiffusionError> {
    for (actual, maximum, name) in [
        (
            data.cell_volumes_m3.len(),
            limits.max_cells,
            "diffusion cells",
        ),
        (
            data.groups.len(),
            limits.max_neutron_groups,
            "diffusion groups",
        ),
        (
            data.delayed_groups.len(),
            limits.max_delayed_groups,
            "diffusion delayed families",
        ),
        (
            data.reactions.len(),
            limits.max_reactions,
            "diffusion reactions",
        ),
        (data.links.len(), limits.max_links, "diffusion interfaces"),
        (
            data.leakage.len(),
            limits.max_leakage_boundaries,
            "diffusion boundaries",
        ),
    ] {
        if actual > maximum {
            return Err(DiffusionError::Budget(name));
        }
    }
    Ok(())
}

fn node_limit(nodes: usize, limits: DiffusionLimits) -> Result<(), DiffusionError> {
    if nodes > limits.propagation.maximum_nodes || nodes > MAX_DECAY_NODES {
        Err(DiffusionError::Budget("diffusion population nodes"))
    } else {
        Ok(())
    }
}

fn count_product(a: usize, b: usize) -> Result<usize, DiffusionError> {
    a.checked_mul(b)
        .ok_or(DiffusionError::Budget("diffusion dimensions"))
}

fn coefficient_limit(count: usize, limits: DiffusionLimits) -> Result<(), DiffusionError> {
    if count > limits.max_coefficients {
        Err(DiffusionError::Budget("diffusion coefficients"))
    } else {
        Ok(())
    }
}

fn add_coefficients(
    count: usize,
    extra: usize,
    limits: DiffusionLimits,
) -> Result<usize, DiffusionError> {
    let count = count
        .checked_add(extra)
        .ok_or(DiffusionError::Budget("diffusion coefficients"))?;
    coefficient_limit(count, limits)?;
    Ok(count)
}

fn reserved<T>(capacity: usize) -> Result<Vec<T>, DiffusionError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|_| DiffusionError::Budget("diffusion memory"))?;
    Ok(result)
}

fn filled<T: Clone>(count: usize, value: T) -> Result<Vec<T>, DiffusionError> {
    let mut result = reserved(count)?;
    result.resize(count, value);
    Ok(result)
}

fn checked(value: f64) -> Result<f64, DiffusionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DiffusionError::Arithmetic("nonfinite diffusion quantity"))
    }
}

// An exported nonnegative defect allowance must not round below its retained
// low part and then pass a budget equal to that rounded-down observation.
fn upward_account_bound(account: Accumulator) -> Result<f64, DiffusionError> {
    let rounded = checked(account.value())?;
    let remainder = checked(sum(&[account.hi, -rounded, account.lo]))?;
    checked(if remainder > 0.0 {
        rounded.next_up()
    } else {
        rounded
    })
}

fn product(a: f64, b: f64) -> Result<f64, DiffusionError> {
    let value = checked(a * b)?;
    if a > 0.0 && b > 0.0 && value == 0.0 {
        Err(DiffusionError::Arithmetic(
            "positive diffusion product underflow",
        ))
    } else {
        Ok(value)
    }
}

fn quotient(a: f64, b: f64) -> Result<f64, DiffusionError> {
    let value = checked(a / b)?;
    if a > 0.0 && value == 0.0 {
        Err(DiffusionError::Arithmetic(
            "positive diffusion ratio underflow",
        ))
    } else {
        Ok(value)
    }
}

fn fraction(rate: f64, total: f64) -> Result<f64, DiffusionError> {
    if rate == 0.0 {
        return Ok(0.0);
    }
    let value = quotient(rate, total)?;
    if !positive(value) || value > 1.0 {
        Err(DiffusionError::Arithmetic(
            "invalid competing diffusion fraction",
        ))
    } else {
        Ok(value)
    }
}

fn push_edge(
    edges: &mut Vec<ReactionEdge>,
    parent: usize,
    daughter: usize,
    fraction: f64,
    yield_: f64,
    limits: DiffusionLimits,
) -> Result<(), DiffusionError> {
    if fraction == 0.0 || yield_ == 0.0 {
        return Ok(());
    }
    if edges.len() >= limits.propagation.maximum_edges || edges.len() >= MAX_DECAY_EDGES {
        return Err(DiffusionError::Budget("diffusion product edges"));
    }
    edges
        .try_reserve(1)
        .map_err(|_| DiffusionError::Budget("diffusion edge memory"))?;
    edges.push(ReactionEdge {
        parent,
        daughter,
        yield_per_removal: product(fraction, yield_)?,
    });
    Ok(())
}

/// Preserve the caller's checkpoint error across the numerical owner's bool API.
fn owner_call<T>(
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    call: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, fs_la::decay::DecayError>,
) -> Result<T, DiffusionError> {
    let mut checkpoint_error = None;
    let result = call(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            checkpoint_error = Some(error);
            true
        }
    });
    if let Some(error) = checkpoint_error {
        return Err(error.into());
    }
    Ok(result?)
}

#[cfg(test)]
mod tests;
