//! Product-resolved irradiation, activation and radioactive decay.
//!
//! Each explicitly authored reaction removes one parent and produces its
//! declared expected yields. The numerical owner [`fs_la::decay::reaction`]
//! evolves all products together, including products born and transformed in the
//! same interval, and integrates event counts independently of endpoint losses.
//! Fission yields are multiplicities, not probabilities, and are never normalized.
//! Empty product lists explicitly export the products from this tracked network.
//!
//! A reaction's energy is released by THAT event. Fission and capture entries
//! contain prompt energy only; subsequent product decays own their own energy.
//! Do not also add an effective fission decay-heat fit or an all-inclusive fission
//! Q value for the same products. These populations retain the irradiation and
//! cooling history without a second effective energy reservoir.
//!
//! Cross sections are either supplied fixed vectors or explicit processed tables
//! evaluated at the interval temperature without extrapolation. The resulting
//! rates and product yields are fixed over that admitted interval. Different
//! group-dependent product yields can be authored as disjoint partial reaction
//! channels; their cross sections must sum to the intended total, not duplicate
//! it. A flux is absolute, not a normalized spectrum. This module does not infer
//! nuclear data, neutron multiplication, self-shielding, radiation deposition,
//! emitted-particle mass/charge accounts or composition-dependent calorics.
//! [`production`] adds explicitly prescribed nuclide supply throughout the same
//! frozen interval, retaining the subsequent reaction ages of supplied nuclei.

use std::{collections::BTreeSet, sync::Arc};

use fs_la::decay::{
    DecayLimits,
    reaction::{ReactionAdvance, ReactionEdge, ReactionPlan},
};
use fs_thermochem::decay::DecayLawV1;

use crate::cross_sections::{CrossSectionTable, NeutronReaction};
use crate::irradiation::{GroupFlux, NeutronGroups, NuclideId};
use crate::{
    Accumulator, EnergyChannels, FissionError, channel_values, nonnegative, positive, sum,
};

/// Continuous nuclide supply into the original product-resolved reaction owner.
pub mod production;

/// Conservative finite-inventory exchange between material compartments.
pub mod transport;

/// Integer capture/decay histories with explicit scheduled neutron irradiation.
pub mod stochastic;

/// A physical identity and its admitted total spontaneous-decay law.
#[derive(Clone, Debug, PartialEq)]
pub struct BurnupNuclide {
    pub id: NuclideId,
    pub law: DecayLawV1,
}

/// One retained product. Yields may exceed one and need not sum to one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BurnupProduct {
    pub nuclide: usize,
    pub yield_per_event: f64,
}

/// A disjoint reaction channel; microscopic cross sections are in square metres.
#[derive(Clone, Debug, PartialEq)]
pub enum BurnupRate {
    /// Positive branching probability of the parent's total spontaneous law.
    /// All decay branches of an unstable parent must sum to one.
    Decay {
        branch_fraction: f64,
    },
    Fission {
        cross_sections_m2: Vec<f64>,
    },
    /// Nonfission removal of the parent; elastic scattering is not capture.
    Capture {
        cross_sections_m2: Vec<f64>,
    },
    /// Nonfission parent removal with outgoing neutrons. These coefficients
    /// evolve the finite target and explicit nuclide products under the supplied
    /// flux; outgoing neutron group yields belong to the coupled transport owner.
    NeutronProduction {
        cross_sections_m2: Vec<f64>,
    },
    /// An identified processed fission table, evaluated at the actual interval T.
    /// Partial-channel tables are required when different channels have distinct yields.
    TabulatedFission {
        cross_sections: Arc<CrossSectionTable>,
    },
    /// A processed nonfission-removal table with the same parent identity and grid.
    /// Radiative capture and the explicit charged-particle absorption channels
    /// retain their original table identity. Products and immediate event energy
    /// are authored separately; no secondary neutrons or negative-Q exchange
    /// are supplied by this parent-removal law. `captures` counts all such events.
    TabulatedCapture {
        cross_sections: Arc<CrossSectionTable>,
    },
    /// An exclusive processed neutron-producing table, including explicit
    /// neutron/triton channels. Its fixed emitted multiplicity remains identified,
    /// but no neutron spectrum or nuclide product is inferred.
    /// This channel is distinct from absorption and fission in every event account.
    TabulatedNeutronProduction {
        cross_sections: Arc<CrossSectionTable>,
    },
}

/// Explicit products and immediate energy of one partial reaction channel.
#[derive(Clone, Debug, PartialEq)]
pub struct BurnupReaction {
    pub data_id: String,
    pub parent: usize,
    pub rate: BurnupRate,
    pub products: Vec<BurnupProduct>,
    /// J/event for this channel, excluding subsequent product decay energy.
    pub energy_per_event: EnergyChannels,
}

/// Caller-owned, source-identified nuclear data. No isotope table is supplied.
#[derive(Clone, Debug, PartialEq)]
pub struct BurnupData {
    pub data_id: String,
    pub groups: Arc<NeutronGroups>,
    pub temperature_min_k: f64,
    pub temperature_max_k: f64,
    pub nuclides: Vec<BurnupNuclide>,
    pub reactions: Vec<BurnupReaction>,
}

/// Explicit bounds for data work and the original dense numerical propagator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BurnupLimits {
    pub max_neutron_groups: usize,
    pub max_reactions: usize,
    /// Total authored retained products across all reaction channels.
    pub max_products: usize,
    pub propagation: DecayLimits,
}

impl Default for BurnupLimits {
    fn default() -> Self {
        Self {
            max_neutron_groups: 256,
            max_reactions: fs_la::decay::MAX_DECAY_EDGES,
            max_products: fs_la::decay::MAX_DECAY_EDGES,
            propagation: DecayLimits::default(),
        }
    }
}

/// Work limits and noncancelling extensive accounting allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BurnupBudget {
    pub limits: BurnupLimits,
    /// Absolute population, reaction-count and product-account defects [nuclei].
    pub nuclei_tolerance: f64,
    /// Absolute event-energy and cumulative-account defects [J].
    pub energy_tolerance_j: f64,
}

impl BurnupBudget {
    pub(crate) fn validate(self) -> Result<(), BurnupError> {
        if !nonnegative(self.nuclei_tolerance) || !nonnegative(self.energy_tolerance_j) {
            return Err(BurnupError::Invalid(
                "finite nonnegative burnup account tolerances required",
            ));
        }
        Ok(())
    }
}

/// Immutable admitted topology, identities, reaction yields and energy channels.
#[derive(Clone, Debug, PartialEq)]
pub struct BurnupModel {
    data: BurnupData,
    products: usize,
}

// Both homogeneous and production maps retain the same physical channel split.
pub(crate) struct ReactionChannels {
    pub(crate) removal_rates: Vec<f64>,
    pub(crate) fractions: Vec<f64>,
    pub(crate) edges: Vec<ReactionEdge>,
}

// Joint compartment propagation supplies physical reaction counts from each
// node's integrated residence. Imports/exports are movement, never reactions or
// externally produced nuclei. The same original physical ledger admits them.
#[derive(Clone, Copy)]
struct TransportAccounts<'a> {
    reaction_counts: &'a [f64],
    reaction_roundoff: &'a [f64],
    imported: &'a [f64],
    exported: &'a [f64],
}

impl BurnupModel {
    /// Admit explicit data and bounded topology. Yields are neither renormalized
    /// nor interpreted as chemical stoichiometry. A label is not data validation.
    pub fn new(data: BurnupData, limits: BurnupLimits) -> Result<Self, BurnupError> {
        if data.data_id.trim().is_empty()
            || data.nuclides.is_empty()
            || !positive(data.temperature_min_k)
            || !positive(data.temperature_max_k)
            || data.temperature_max_k < data.temperature_min_k
        {
            return Err(BurnupError::Invalid(
                "identified nonempty nuclear data and a finite temperature domain required",
            ));
        }
        let products = data.reactions.iter().try_fold(0_usize, |n, reaction| {
            n.checked_add(reaction.products.len())
                .ok_or(BurnupError::Budget("product count"))
        })?;
        let model = Self { data, products };
        model.admit_limits(limits)?;
        let mut identities = BTreeSet::new();
        for nuclide in &model.data.nuclides {
            if !identities.insert(nuclide.id) {
                return Err(BurnupError::Invalid("nuclide identities must be unique"));
            }
        }
        let mut branches = vec![Accumulator::new(0.0); model.data.nuclides.len()];
        for reaction in &model.data.reactions {
            if reaction.data_id.trim().is_empty() || reaction.parent >= model.data.nuclides.len() {
                return Err(BurnupError::Invalid(
                    "identified reaction with a represented parent required",
                ));
            }
            reaction.energy_per_event.validate()?;
            match &reaction.rate {
                BurnupRate::Decay { branch_fraction } => {
                    if !positive(*branch_fraction)
                        || *branch_fraction > 1.0
                        || model.data.nuclides[reaction.parent]
                            .law
                            .decay_constant()
                            .value()
                            == 0.0
                    {
                        return Err(BurnupError::Invalid(
                            "decay branches require an unstable parent and positive probability",
                        ));
                    }
                    add(&mut branches[reaction.parent], *branch_fraction)?;
                }
                BurnupRate::Fission { cross_sections_m2 }
                | BurnupRate::Capture { cross_sections_m2 }
                | BurnupRate::NeutronProduction { cross_sections_m2 } => {
                    if cross_sections_m2.len() != model.data.groups.len()
                        || !cross_sections_m2.iter().copied().all(nonnegative)
                    {
                        return Err(BurnupError::Invalid(
                            "one finite nonnegative microscopic cross section per group required",
                        ));
                    }
                }
                BurnupRate::TabulatedFission { cross_sections } => {
                    model.admit_table(cross_sections, reaction.parent, NeutronReaction::Fission)?;
                }
                BurnupRate::TabulatedCapture { cross_sections } => {
                    model.admit_table(cross_sections, reaction.parent, NeutronReaction::Capture)?;
                }
                BurnupRate::TabulatedNeutronProduction { cross_sections } => {
                    if cross_sections
                        .reaction()
                        .emitted_neutrons_per_event()
                        .is_none()
                    {
                        return Err(BurnupError::Invalid(
                            "neutron production requires an exclusive neutron-producing table",
                        ));
                    }
                    model.admit_table(
                        cross_sections,
                        reaction.parent,
                        cross_sections.reaction(),
                    )?;
                }
            }
            let mut products = BTreeSet::new();
            for product in &reaction.products {
                if product.nuclide >= model.data.nuclides.len()
                    || !positive(product.yield_per_event)
                    || !products.insert(product.nuclide)
                {
                    return Err(BurnupError::Invalid(
                        "products require unique represented identities and finite positive yields",
                    ));
                }
            }
        }
        for (nuclide, branches) in model.data.nuclides.iter().zip(branches) {
            if nuclide.law.decay_constant().value() > 0.0
                && (branches.value() - 1.0).abs() > 32.0 * f64::EPSILON
            {
                return Err(BurnupError::Invalid(
                    "all spontaneous decay branches must be explicit and sum to one",
                ));
            }
        }
        // The original numerical owner admits its own resource limits and shape.
        // Physical channel topology was checked above; rates are supplied later.
        ReactionPlan::prepare(
            &vec![0.0; model.data.nuclides.len()],
            &[],
            0.0,
            limits.propagation,
            || false,
        )?;
        Ok(model)
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
    pub fn nuclides(&self) -> &[BurnupNuclide] {
        &self.data.nuclides
    }
    #[must_use]
    pub fn reactions(&self) -> &[BurnupReaction] {
        &self.data.reactions
    }
    #[must_use]
    pub const fn temperature_min_k(&self) -> f64 {
        self.data.temperature_min_k
    }
    #[must_use]
    pub const fn temperature_max_k(&self) -> f64 {
        self.data.temperature_max_k
    }

    /// Every interval temperature must lie inside the model's declared domain.
    /// Construction checks that every supplied table covers this entire domain.
    pub fn validate_temperature(&self, temperature_k: f64) -> Result<(), BurnupError> {
        if !positive(temperature_k)
            || temperature_k < self.temperature_min_k()
            || temperature_k > self.temperature_max_k()
        {
            return Err(BurnupError::Invalid(
                "temperature is outside the burnup data domain",
            ));
        }
        Ok(())
    }

    fn admit_limits(&self, limits: BurnupLimits) -> Result<(), BurnupError> {
        if self.data.groups.len() > limits.max_neutron_groups {
            return Err(BurnupError::Budget("neutron groups"));
        }
        if self.data.reactions.len() > limits.max_reactions {
            return Err(BurnupError::Budget("reaction channels"));
        }
        if self.products > limits.max_products
            || self.products > limits.propagation.maximum_edges
            || self.products > fs_la::decay::MAX_DECAY_EDGES
        {
            return Err(BurnupError::Budget("retained products"));
        }
        if self.data.nuclides.len() > limits.propagation.maximum_nodes
            || self.data.nuclides.len() > fs_la::decay::MAX_DECAY_NODES
        {
            return Err(BurnupError::Budget("nuclides"));
        }
        Ok(())
    }

    fn admit_table(
        &self,
        table: &CrossSectionTable,
        parent: usize,
        reaction: NeutronReaction,
    ) -> Result<(), BurnupError> {
        let same_kind = if reaction == NeutronReaction::Capture {
            table.reaction().is_nonfission_absorption()
        } else {
            table.reaction() == reaction
        };
        if table.nuclide() != self.data.nuclides[parent].id || !same_kind {
            return Err(BurnupError::Invalid(
                "tabulated reaction must match the parent nuclide and reaction kind",
            ));
        }
        if !Arc::ptr_eq(self.groups(), table.groups()) && self.groups() != table.groups() {
            return Err(BurnupError::Invalid(
                "tabulated reaction must use the burnup model's exact neutron groups",
            ));
        }
        if table.temperature_min_k() > self.temperature_min_k()
            || table.temperature_max_k() < self.temperature_max_k()
        {
            return Err(BurnupError::Invalid(
                "tabulated reaction must cover the entire burnup temperature domain",
            ));
        }
        Ok(())
    }

    /// Per-parent rate constant for each original reaction [1/s]. These rounded
    /// rates define the frozen operator. Only explicitly supplied processed
    /// tables interpolate microscopic coefficients at the current temperature.
    pub fn rate_constants(
        &self,
        flux: &GroupFlux,
        temperature_k: f64,
        limits: BurnupLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<f64>, BurnupError> {
        checkpoint()?;
        self.admit_limits(limits)?;
        self.validate_temperature(temperature_k)?;
        if !Arc::ptr_eq(self.groups(), flux.groups()) && self.groups() != flux.groups() {
            return Err(BurnupError::Invalid(
                "burnup flux and microscopic data groups must match exactly",
            ));
        }
        let mut rates = Vec::with_capacity(self.data.reactions.len());
        for reaction in &self.data.reactions {
            checkpoint()?;
            let rate = match &reaction.rate {
                BurnupRate::Decay { branch_fraction } => product(
                    self.data.nuclides[reaction.parent]
                        .law
                        .decay_constant()
                        .value(),
                    *branch_fraction,
                )?,
                BurnupRate::Fission { cross_sections_m2 }
                | BurnupRate::Capture { cross_sections_m2 }
                | BurnupRate::NeutronProduction { cross_sections_m2 } => {
                    let mut rate = Accumulator::new(0.0);
                    for (&cross_section, &flux) in cross_sections_m2.iter().zip(flux.values_m2_s())
                    {
                        checkpoint()?;
                        add(&mut rate, product(cross_section, flux)?)?;
                    }
                    checked(rate.value())?
                }
                BurnupRate::TabulatedFission { cross_sections }
                | BurnupRate::TabulatedCapture { cross_sections }
                | BurnupRate::TabulatedNeutronProduction { cross_sections } => cross_sections
                    .rate_constant(
                        flux,
                        temperature_k,
                        limits.max_neutron_groups,
                        &mut checkpoint,
                    )?,
            };
            rates.push(rate);
        }
        checkpoint()?;
        Ok(rates)
    }

    /// Prepare one frozen-flux interval and reuse it across inventories bound to
    /// this exact immutable model. Matrix application is quadratic in nuclides.
    pub fn prepare(
        self: &Arc<Self>,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        limits: BurnupLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<PreparedBurnup, BurnupError> {
        let rates = self.rate_constants(flux, temperature_k, limits, &mut checkpoint)?;
        self.prepare_rates(duration_s, flux, temperature_k, limits, rates, checkpoint)
    }

    // One channel compiler serves both numerical map families. Rates have
    // already been admitted and collapsed at the actual flux and temperature.
    pub(crate) fn prepare_channels(
        &self,
        rates: &[f64],
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<ReactionChannels, BurnupError> {
        let mut removals = vec![Accumulator::new(0.0); self.data.nuclides.len()];
        for (reaction, &rate) in self.data.reactions.iter().zip(rates) {
            checkpoint()?;
            add(&mut removals[reaction.parent], rate)?;
        }
        let removal_rates: Vec<_> = removals.iter().map(|rate| rate.value()).collect();
        let mut fractions = Vec::with_capacity(rates.len());
        let mut edges = Vec::with_capacity(self.products);
        for (reaction, &rate) in self.data.reactions.iter().zip(rates) {
            checkpoint()?;
            let fraction = if rate == 0.0 {
                0.0
            } else {
                let ratio = rate / removal_rates[reaction.parent];
                if !positive(ratio) {
                    return Err(BurnupError::Arithmetic(
                        "a positive reaction fraction is unrepresentable",
                    ));
                }
                ratio
            };
            fractions.push(fraction);
            if fraction != 0.0 {
                for yield_ in &reaction.products {
                    checkpoint()?;
                    edges.push(ReactionEdge {
                        parent: reaction.parent,
                        daughter: yield_.nuclide,
                        yield_per_removal: product(fraction, yield_.yield_per_event)?,
                    });
                }
            }
        }
        Ok(ReactionChannels {
            removal_rates,
            fractions,
            edges,
        })
    }

    // Public preparation and automatic cache misses use the already collapsed
    // channel rates, without repeating table interpolation or neutron-group work.
    fn prepare_rates(
        self: &Arc<Self>,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        limits: BurnupLimits,
        rates: Vec<f64>,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<PreparedBurnup, BurnupError> {
        let channels = self.prepare_channels(&rates, &mut checkpoint)?;
        let mut checkpoint_error = None;
        let plan = ReactionPlan::prepare(
            &channels.removal_rates,
            &channels.edges,
            duration_s,
            limits.propagation,
            || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    checkpoint_error = Some(error);
                    true
                }
            },
        );
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        Ok(PreparedBurnup {
            model: Arc::clone(self),
            plan: plan?,
            flux: flux.clone(),
            temperature_k,
            limits,
            rates,
            fractions: channels.fractions,
        })
    }
}

/// Reusable interval maps and reaction attribution with exact model binding.
/// Flux and temperature record the original preparation point. An inventory's
/// cache may reuse these maps at other conditions only when ALL collapsed
/// reaction rates are identical, preserving the channel-by-channel attribution.
#[derive(Clone, Debug)]
pub struct PreparedBurnup {
    model: Arc<BurnupModel>,
    plan: ReactionPlan,
    flux: GroupFlux,
    temperature_k: f64,
    limits: BurnupLimits,
    rates: Vec<f64>,
    fractions: Vec<f64>,
}

impl PreparedBurnup {
    #[must_use]
    pub fn model(&self) -> &Arc<BurnupModel> {
        &self.model
    }
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.plan.duration_s()
    }
    #[must_use]
    pub fn flux(&self) -> &GroupFlux {
        &self.flux
    }
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    #[must_use]
    pub fn reaction_rate_constants(&self) -> &[f64] {
        &self.rates
    }
    #[must_use]
    pub fn propagator(&self) -> &ReactionPlan {
        &self.plan
    }
}

/// Accepted population history and compensated cumulative reaction/energy accounts.
#[derive(Clone, Debug)]
pub struct BurnupInventory {
    model: Arc<BurnupModel>,
    populations: Vec<f64>,
    reactions: Vec<Accumulator>,
    produced: Vec<Accumulator>,
    released: [Accumulator; 3],
    // One immutable map is enough for the common repeated-step case. This is
    // an optimization, not physical history, and is excluded from equality.
    prepared: Option<Arc<PreparedBurnup>>,
    production_prepared: Option<Arc<production::PreparedBurnupProduction>>,
}

impl PartialEq for BurnupInventory {
    fn eq(&self, other: &Self) -> bool {
        self.model == other.model
            && self.populations == other.populations
            && self.reactions == other.reactions
            && self.produced == other.produced
            && self.released == other.released
    }
}

impl BurnupInventory {
    /// Supply every initial expected population, including existing products.
    /// No zero-product, equilibrium or previously irradiated history is assumed.
    pub fn new(model: Arc<BurnupModel>, populations: Vec<f64>) -> Result<Self, BurnupError> {
        if populations.len() != model.nuclides().len()
            || !populations.iter().copied().all(nonnegative)
            || !sum(&populations).is_finite()
        {
            return Err(BurnupError::Invalid(
                "one finite nonnegative population per nuclide required",
            ));
        }
        let reactions = vec![Accumulator::new(0.0); model.reactions().len()];
        let produced = vec![Accumulator::new(0.0); populations.len()];
        Ok(Self {
            model,
            populations,
            reactions,
            produced,
            released: [Accumulator::new(0.0); 3],
            prepared: None,
            production_prepared: None,
        })
    }

    #[must_use]
    pub fn model(&self) -> &Arc<BurnupModel> {
        &self.model
    }
    #[must_use]
    pub fn populations(&self) -> &[f64] {
        &self.populations
    }
    /// Lifetime external production, separate from reaction-created products.
    #[must_use]
    pub fn cumulative_produced_nuclei(&self) -> Vec<f64> {
        self.produced.iter().map(|count| count.value()).collect()
    }
    #[must_use]
    pub fn cumulative_reaction_counts(&self) -> Vec<f64> {
        self.reactions.iter().map(|count| count.value()).collect()
    }
    #[must_use]
    pub fn cumulative_released(&self) -> EnergyChannels {
        channel_values(self.released)
    }

    /// Current event rates [1/s], for diagnostics. Use integrated interval counts
    /// for heat; endpoint rates times duration omit within-interval production.
    pub fn reaction_rates(
        &self,
        flux: &GroupFlux,
        temperature_k: f64,
        limits: BurnupLimits,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<f64>, BurnupError> {
        self.model
            .rate_constants(flux, temperature_k, limits, checkpoint)?
            .iter()
            .zip(self.model.reactions())
            .map(|(&rate, reaction)| product(rate, self.populations[reaction.parent]))
            .collect()
    }

    pub fn trial(
        &self,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        budget: BurnupBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, BurnupStep), BurnupError> {
        checkpoint()?;
        budget.validate()?;
        // Reevaluate the current operator: matching dt/flux is insufficient once
        // cross sections depend on temperature. Compare individual channel rates,
        // not just total removal, because product yields and event energies differ.
        let rates =
            self.model
                .rate_constants(flux, temperature_k, budget.limits, &mut checkpoint)?;
        let prepared = if let Some(prepared) = self.prepared.as_ref().filter(|prepared| {
            Arc::ptr_eq(&self.model, prepared.model())
                && prepared.duration_s().to_bits() == duration_s.to_bits()
                && prepared.limits == budget.limits
                && prepared.reaction_rate_constants() == rates.as_slice()
        }) {
            Arc::clone(prepared)
        } else {
            Arc::new(self.model.prepare_rates(
                duration_s,
                flux,
                temperature_k,
                budget.limits,
                rates,
                &mut checkpoint,
            )?)
        };
        let (mut next, report) = self.trial_prepared(&prepared, budget, checkpoint)?;
        next.prepared = Some(prepared);
        Ok((next, report))
    }

    /// Apply without matrix preparation. Numerical/data limits must match the
    /// preparation; extensive accounting tolerances may vary with inventory size.
    pub fn trial_prepared(
        &self,
        prepared: &PreparedBurnup,
        budget: BurnupBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, BurnupStep), BurnupError> {
        checkpoint()?;
        budget.validate()?;
        if !Arc::ptr_eq(&self.model, prepared.model()) || budget.limits != prepared.limits {
            return Err(BurnupError::Invalid(
                "prepared burnup requires the original model and work limits",
            ));
        }
        let mut checkpoint_error = None;
        let evolved = prepared
            .plan
            .apply_with_cancel(&self.populations, || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    checkpoint_error = Some(error);
                    true
                }
            });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let evolved = evolved?;
        let (mut next, report) = self.account_evolution(
            evolved,
            prepared.duration_s(),
            &prepared.fractions,
            None,
            0.0,
            budget,
            &mut checkpoint,
        )?;
        next.production_prepared = None;
        Ok((next, report))
    }

    // Both map families share the physical channel/product/energy ledger. A
    // supplied nucleus enters this account once, separately from reaction yields.
    pub(crate) fn account_evolution(
        &self,
        evolved: ReactionAdvance,
        duration_s: f64,
        fractions: &[f64],
        produced: Option<&[f64]>,
        production_roundoff_nuclei: f64,
        budget: BurnupBudget,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, BurnupStep), BurnupError> {
        self.account_events(
            evolved.populations,
            evolved.removals,
            duration_s,
            fractions,
            None,
            produced,
            production_roundoff_nuclei,
            budget,
            checkpoint,
        )
    }

    // All physical reaction products, energy and lifetime histories have one
    // owner. A joint transport solve changes only how event counts and external
    // movement are supplied to this ledger, not their physical interpretation.
    fn account_events(
        &self,
        populations: Vec<f64>,
        parent_removals: Vec<f64>,
        duration_s: f64,
        fractions: &[f64],
        transport: Option<TransportAccounts<'_>>,
        produced: Option<&[f64]>,
        input_roundoff_nuclei: f64,
        budget: BurnupBudget,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, BurnupStep), BurnupError> {
        let mut next = self.clone();
        let n = self.populations.len();
        let mut incoming = vec![Accumulator::new(0.0); n];
        let mut outgoing = vec![Accumulator::new(0.0); n];
        let mut counts = [Accumulator::new(0.0); 4];
        let mut generated = [Accumulator::new(0.0); 4];
        let mut channels = [Accumulator::new(0.0); 3];
        let mut nuclei_error = Accumulator::new(input_roundoff_nuclei);
        let mut energy_error = Accumulator::new(0.0);
        let mut reaction_counts = Vec::with_capacity(self.model.reactions().len());
        for (index, reaction) in self.model.reactions().iter().enumerate() {
            checkpoint()?;
            let (count, count_error) = if let Some(accounts) = transport {
                (
                    accounts.reaction_counts[index],
                    accounts.reaction_roundoff[index],
                )
            } else {
                let removal = parent_removals[reaction.parent];
                let fraction = fractions[index];
                let count = product(removal, fraction)?;
                (count, removal.mul_add(fraction, -count).abs())
            };
            add(&mut nuclei_error, count_error)?;
            let error = add(&mut outgoing[reaction.parent], count)?.abs();
            add(&mut nuclei_error, error)?;
            let error = add(&mut next.reactions[index], count)?.abs();
            add(&mut nuclei_error, error)?;
            let kind = match reaction.rate {
                BurnupRate::Fission { .. } | BurnupRate::TabulatedFission { .. } => 0,
                BurnupRate::Capture { .. } | BurnupRate::TabulatedCapture { .. } => 1,
                BurnupRate::Decay { .. } => 2,
                BurnupRate::NeutronProduction { .. }
                | BurnupRate::TabulatedNeutronProduction { .. } => 3,
            };
            let error = add(&mut counts[kind], count)?.abs();
            add(&mut nuclei_error, error)?;
            reaction_counts.push(count);
            for yield_ in &reaction.products {
                checkpoint()?;
                let population = product(count, yield_.yield_per_event)?;
                add(
                    &mut nuclei_error,
                    count.mul_add(yield_.yield_per_event, -population).abs(),
                )?;
                let error = add(&mut incoming[yield_.nuclide], population)?.abs();
                add(&mut nuclei_error, error)?;
            }
            let q = reaction.energy_per_event.total_j();
            add(&mut energy_error, product(count_error, q)?)?;
            let event_energy = product(count, q)?;
            add(&mut energy_error, count.mul_add(q, -event_energy).abs())?;
            let error = add(&mut generated[kind], event_energy)?.abs();
            add(&mut energy_error, error)?;
            for (channel, coefficient) in channels
                .iter_mut()
                .zip(channel_array(reaction.energy_per_event))
            {
                let energy = product(count, coefficient)?;
                add(&mut energy_error, count.mul_add(coefficient, -energy).abs())?;
                let error = add(channel, energy)?.abs();
                add(&mut energy_error, error)?;
            }
        }
        let mut population_balance_residuals = Vec::with_capacity(n);
        for i in 0..n {
            checkpoint()?;
            let residual = checked(sum(&[
                populations[i],
                -self.populations[i],
                outgoing[i].hi,
                outgoing[i].lo,
                -incoming[i].hi,
                -incoming[i].lo,
                -produced.map_or(0.0, |counts| counts[i]),
                -transport.map_or(0.0, |accounts| accounts.imported[i]),
                transport.map_or(0.0, |accounts| accounts.exported[i]),
            ]))?;
            population_balance_residuals.push(residual);
            add(&mut nuclei_error, residual.abs())?;
            let attribution_residual =
                checked(sum(&[parent_removals[i], -outgoing[i].hi, -outgoing[i].lo]))?;
            add(&mut nuclei_error, attribution_residual.abs())?;
            if let Some(produced) = produced {
                let error = add(&mut next.produced[i], produced[i])?.abs();
                add(&mut nuclei_error, error)?;
            }
        }
        for &value in &counts {
            add(&mut nuclei_error, projection_error(value)?)?;
        }
        for &value in channels.iter().chain(&generated) {
            add(&mut energy_error, projection_error(value)?)?;
        }
        let released = channel_values(channels);
        checked(released.total_j())?;
        for (account, energy) in next.released.iter_mut().zip(channel_array(released)) {
            let error = add(account, energy)?.abs();
            add(&mut energy_error, error)?;
        }
        checked(next.cumulative_released().total_j())?;
        let generated_energy_j = checked(sum(&generated.map(|value| value.value())))?;
        let energy_residual_j = checked(sum(&[
            generated_energy_j,
            -released.deposited_j,
            -released.radiation_j,
            -released.escaped_j,
        ]))?;
        add(&mut energy_error, energy_residual_j.abs())?;
        let population_roundoff_nuclei = nonnegative_bound(nuclei_error)?;
        let energy_roundoff_j = nonnegative_bound(energy_error)?;
        if population_roundoff_nuclei > budget.nuclei_tolerance {
            return Err(BurnupError::Budget(
                "population and reaction accounting roundoff",
            ));
        }
        if energy_roundoff_j > budget.energy_tolerance_j {
            return Err(BurnupError::Budget("event energy accounting roundoff"));
        }
        next.populations = populations;
        checkpoint()?;
        Ok((
            next,
            BurnupStep {
                duration_s,
                reaction_counts,
                parent_removals,
                fissions: counts[0].value(),
                captures: counts[1].value(),
                decays: counts[2].value(),
                neutron_producing_reactions: counts[3].value(),
                released,
                fission_energy_j: generated[0].value(),
                capture_energy_j: generated[1].value(),
                decay_energy_j: generated[2].value(),
                neutron_production_energy_j: generated[3].value(),
                generated_energy_j,
                energy_residual_j,
                energy_roundoff_j,
                population_balance_residuals,
                population_roundoff_nuclei,
            },
        ))
    }

    /// Publish only a complete accepted candidate; every refusal preserves self.
    pub fn advance(
        &mut self,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        budget: BurnupBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupStep, BurnupError> {
        let (next, report) = self.trial(duration_s, flux, temperature_k, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }

    pub fn advance_prepared(
        &mut self,
        prepared: &PreparedBurnup,
        budget: BurnupBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupStep, BurnupError> {
        let (next, report) = self.trial_prepared(prepared, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }
}

/// Integrated event and energy accounts in the original reaction/nuclide order.
/// Population residuals are nodewise production-minus-removal identities; total
/// nucleus count is not conserved when a reaction has multiple products.
#[derive(Clone, Debug, PartialEq)]
pub struct BurnupStep {
    pub duration_s: f64,
    pub reaction_counts: Vec<f64>,
    pub parent_removals: Vec<f64>,
    pub fissions: f64,
    pub captures: f64,
    pub decays: f64,
    /// Nonfission parent-removal events with outgoing neutrons, not the number
    /// of emitted neutrons. Their explicit yields belong to neutron transport.
    pub neutron_producing_reactions: f64,
    pub released: EnergyChannels,
    pub fission_energy_j: f64,
    pub capture_energy_j: f64,
    pub decay_energy_j: f64,
    /// Separately authored energy released by neutron-producing events [J].
    /// Signed prompt Q and incident/outgoing kinetic energy, when assigned to a
    /// coupled neutron-energy owner, must be excluded from this channel.
    pub neutron_production_energy_j: f64,
    /// Sum of energy released by actual events, not eventual energy of unreacted
    /// products. Product inventories, not an extra fit, own future decay heat.
    pub generated_energy_j: f64,
    pub energy_residual_j: f64,
    pub energy_roundoff_j: f64,
    pub population_balance_residuals: Vec<f64>,
    pub population_roundoff_nuclei: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BurnupError {
    Invalid(&'static str),
    Budget(&'static str),
    Arithmetic(&'static str),
    Fission(FissionError),
    Propagation(fs_la::decay::DecayError),
}

impl From<FissionError> for BurnupError {
    fn from(error: FissionError) -> Self {
        Self::Fission(error)
    }
}
impl From<fs_la::decay::DecayError> for BurnupError {
    fn from(error: fs_la::decay::DecayError) -> Self {
        Self::Propagation(error)
    }
}
impl core::fmt::Display for BurnupError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "product-resolved burnup refused: {self:?}")
    }
}
impl std::error::Error for BurnupError {}

fn channel_array(energy: EnergyChannels) -> [f64; 3] {
    [energy.deposited_j, energy.radiation_j, energy.escaped_j]
}
fn checked(value: f64) -> Result<f64, BurnupError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BurnupError::Arithmetic("nonfinite nuclear account"))
    }
}
fn product(a: f64, b: f64) -> Result<f64, BurnupError> {
    let value = checked(a * b)?;
    if a > 0.0 && b > 0.0 && value == 0.0 {
        Err(BurnupError::Arithmetic(
            "a positive nuclear quantity is unrepresentable",
        ))
    } else {
        Ok(value)
    }
}
fn add(target: &mut Accumulator, value: f64) -> Result<f64, BurnupError> {
    Ok(target.add(value)?)
}
fn projection_error(value: Accumulator) -> Result<f64, BurnupError> {
    Ok(checked(sum(&[value.value(), -value.hi, -value.lo]))?.abs())
}

// Preserve a positive low part at the exact allowance boundary instead of
// rounding a noncancelling defect sum downward into an apparently valid step.
fn nonnegative_bound(value: Accumulator) -> Result<f64, BurnupError> {
    let rounded = checked(value.value())?;
    if checked(sum(&[value.hi, value.lo, -rounded]))? > 0.0 {
        checked(rounded.next_up())
    } else {
        Ok(rounded)
    }
}

#[cfg(test)]
mod tests;
