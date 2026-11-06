//! Constant first-order chemistry on the original constrained reactor.
//!
//! Select [`super::IsochoricKineticsMethod::ExactFirstOrder`] explicitly. Each
//! irreversible channel has one order-one reactant and an authored constant
//! Arrhenius coefficient. Competing channels add their hazards, while the
//! original product multiplicities remain in the positive reaction generator.
//! No reverse law, temperature dependence, collider or source data is inferred.
//!
//! The shared L1 owner evolves populations and integrates their residence.
//! Channel j accumulates `k_j * integral(N_parent dt)`, including reactants born
//! during the interval. Counts are never inferred from endpoint differences or
//! from a rounded fraction of a competing channel's counts. The original
//! reactor then resolves U,V or H,p with those events and original calorics.
//! `maximum_step_s` bounds the spacing of caloric/entropy admission checks;
//! those checks do not certify every intermediate point of the trajectory.
//!
//! A prepared reactor caches one interval map by its exact duration. Applying
//! the shared checked maps uses bounded temporary arrays, unlike the existing
//! allocation-free implicit stages. Floating-point map accuracy, dense resource
//! limits and underflow refusals remain those of `fs_la::decay::reaction`.

use super::{IsochoricKineticsError, Sum, poll};
use fs_la::decay::reaction::{ReactionEdge, ReactionPlan};
use fs_la::decay::{DecayError, MAX_DECAY_EDGES, MAX_DECAY_NODES};
use fs_thermochem::kinetics::{ForwardRateV1, IdealGasKineticsV1, PressureDependenceV1};

/// Existing shared propagator node, edge, series, work and accuracy limits.
pub use fs_la::decay::DecayLimits as FirstOrderLimits;

struct Channel {
    parent: usize,
    hazard_per_s: f64,
}

pub(super) struct PreparedFirstOrder {
    rates: Vec<f64>,
    edges: Vec<ReactionEdge>,
    channels: Vec<Channel>,
    limits: FirstOrderLimits,
    cached: ReactionPlan,
}

impl PreparedFirstOrder {
    pub(super) fn new(
        mechanism: &IdealGasKineticsV1,
        limits: FirstOrderLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, IsochoricKineticsError> {
        poll(cancelled)?;
        let species = mechanism.stoichiometry().species();
        let maximum_nodes = limits.maximum_nodes.min(MAX_DECAY_NODES);
        if species.len() > maximum_nodes {
            return Err(DecayError::Budget {
                quantity: "first-order chemical species",
                required: species.len() as u64,
                maximum: maximum_nodes as u64,
            }
            .into());
        }
        let mut sums = vec![Sum::default(); species.len()];
        let mut channels = Vec::with_capacity(mechanism.reactions().len());
        let mut edge_count = 0;
        for reaction in mechanism.reactions() {
            poll(cancelled)?;
            let ForwardRateV1::Arrhenius(law) = reaction.forward_rate() else {
                return Err(IsochoricKineticsError::Invalid(
                    "exact first-order chemistry requires an Arrhenius forward law",
                ));
            };
            if reaction.is_reversible()
                || reaction.reactants().len() != 1
                || reaction.reactants()[0].1 != 1
                || law.temperature_exponent() != 0.0
                || law.activation_energy().value() != 0.0
                || !matches!(
                    reaction.pressure_dependence(),
                    PressureDependenceV1::Independent
                )
            {
                return Err(IsochoricKineticsError::Invalid(
                    "exact first-order chemistry requires irreversible constant single-reactant channels",
                ));
            }
            let parent = species
                .binary_search(&reaction.reactants()[0].0)
                .expect("original mechanism admitted its reactant identity");
            // The source rate uses c*=1 mol/m3. With one order-one reactant,
            // its numerical coefficient is the inverse-second hazard.
            let hazard_per_s = law.reference_rate().value();
            sums[parent].add(hazard_per_s);
            channels.push(Channel {
                parent,
                hazard_per_s,
            });
            if hazard_per_s > 0.0 {
                edge_count += reaction.products().len();
            }
        }
        let maximum_edges = limits.maximum_edges.min(MAX_DECAY_EDGES);
        if edge_count > maximum_edges {
            return Err(DecayError::Budget {
                quantity: "first-order chemical products",
                required: edge_count as u64,
                maximum: maximum_edges as u64,
            }
            .into());
        }
        let rates: Vec<_> = sums.iter().map(Sum::value).collect();
        if rates.iter().any(|rate| !rate.is_finite()) {
            return Err(DecayError::Unrepresentable("summed chemical removal hazard").into());
        }
        let mut edges = Vec::with_capacity(edge_count);
        for (reaction, channel) in mechanism.reactions().iter().zip(&channels) {
            poll(cancelled)?;
            if channel.hazard_per_s == 0.0 {
                continue;
            }
            let fraction = channel.hazard_per_s / rates[channel.parent];
            if fraction == 0.0 {
                return Err(
                    DecayError::Unrepresentable("competing chemical channel fraction").into(),
                );
            }
            for (product, multiplicity) in reaction.products() {
                edges.push(ReactionEdge {
                    parent: channel.parent,
                    daughter: species
                        .binary_search(product)
                        .expect("original mechanism admitted its product identity"),
                    yield_per_removal: fraction * f64::from(*multiplicity),
                });
            }
        }
        // The zero interval validates the shared resource policy at construction
        // and is also the correct cached map for a zero-duration application.
        let cached = ReactionPlan::prepare_with_population_integrals(
            &rates,
            &edges,
            0.0,
            limits,
            &mut *cancelled,
        )?;
        Ok(Self {
            rates,
            edges,
            channels,
            limits,
            cached,
        })
    }

    pub(super) fn binding_bytes(&self) -> usize {
        self.rates.capacity() * core::mem::size_of::<f64>()
            + self.edges.capacity() * core::mem::size_of::<ReactionEdge>()
            + self.channels.capacity() * core::mem::size_of::<Channel>()
    }

    /// Input and output are amounts divided by the same immutable reference
    /// volume. Thus this generator applies to both actual U,V and H,p paths.
    pub(super) fn advance(
        &mut self,
        mechanism: &IdealGasKineticsV1,
        concentrations: &[f64],
        duration_s: f64,
        extents: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, IsochoricKineticsError> {
        poll(cancelled)?;
        if extents.len() != self.channels.len() {
            return Err(IsochoricKineticsError::Invalid(
                "first-order extent workspace does not match the mechanism",
            ));
        }
        if self.cached.duration_s().to_bits() != duration_s.to_bits() {
            let candidate = ReactionPlan::prepare_with_population_integrals(
                &self.rates,
                &self.edges,
                duration_s,
                self.limits,
                &mut *cancelled,
            )?;
            self.cached = candidate;
        }
        let endpoint = self
            .cached
            .apply_with_cancel(concentrations, &mut *cancelled)?;
        let residence = self
            .cached
            .integrate_populations_with_cancel(concentrations, &mut *cancelled)?;
        let mut by_parent = [Sum::default(); MAX_DECAY_NODES];
        for (extent, channel) in extents.iter_mut().zip(&self.channels) {
            poll(cancelled)?;
            let amount_time = residence[channel.parent];
            let event_count = channel.hazard_per_s * amount_time;
            if !event_count.is_finite()
                || (channel.hazard_per_s > 0.0 && amount_time > 0.0 && event_count == 0.0)
            {
                return Err(
                    DecayError::Unrepresentable("integrated chemical channel events").into(),
                );
            }
            *extent = event_count;
            by_parent[channel.parent].add(event_count);
        }
        // Check the independent channel attribution against the same positive
        // owner's integrated removals, without subtracting endpoint populations.
        for (i, sum) in by_parent.iter().enumerate().take(self.rates.len()) {
            let count = sum.value();
            let removed = endpoint.removals[i];
            let residual = count - removed;
            let allowance = self.limits.relative_tolerance * count.max(removed);
            if !residual.is_finite() || residual.abs() > allowance {
                return Err(DecayError::Balance {
                    row: i,
                    column: self.rates.len(),
                    residual,
                    allowance,
                }
                .into());
            }
        }
        // Check original chemical stoichiometry on the represented amount
        // scale, with an explicit counter-quantization allowance. In a late
        // A->B->C chain both counters can round to one while B still has a
        // resolved positive tail; that tail must not be erased or rejected.
        // A relative tolerance times enormous cyclic event totals would hide
        // much larger errors, so only counterproduct ULPs augment this check.
        for i in 0..concentrations.len() {
            poll(cancelled)?;
            let mut discrepancy = Sum::default();
            let mut counter_roundoff = Sum::default();
            discrepancy.add(endpoint.populations[i]);
            discrepancy.add(-concentrations[i]);
            for (j, &extent) in extents.iter().enumerate() {
                let nu = mechanism
                    .stoichiometry()
                    .get(i, j)
                    .expect("original chemical stoichiometry has the admitted shape")
                    as f64;
                let term = nu * extent;
                discrepancy.add(-term);
                if term != 0.0 {
                    let magnitude = term.abs();
                    counter_roundoff.add(2.0 * (magnitude.next_up() - magnitude));
                }
            }
            let residual = discrepancy.value();
            let allowance = self.limits.relative_tolerance
                * endpoint.populations[i].max(concentrations[i])
                + counter_roundoff.value();
            if !residual.is_finite() || !allowance.is_finite() || residual.abs() > allowance {
                return Err(DecayError::Balance {
                    row: i,
                    column: concentrations.len(),
                    residual,
                    allowance,
                }
                .into());
            }
        }
        poll(cancelled)?;
        Ok(endpoint.populations)
    }
}
