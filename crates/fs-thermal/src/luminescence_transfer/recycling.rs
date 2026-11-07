//! Photon-number recycling through coupled singlet/triplet excitation stores.
//!
//! Spectral photons retain their actual energy. An absorbed photon produces
//! at most one explicitly selected receiver excitation, and only when it can
//! pay that state's energy. Source Stokes loss, receiver relaxation, passive
//! absorption, nonradiative decay and crossing heat remain distinct accounts.
//! The original positive reaction/source owner integrates all repeated visits
//! and uniform-pump ages together, without a bounce cutoff or fake lifetime.
//!
//! Rates, photon spectra and transport stay fixed throughout one interval.
//! Photon flight time is neglected on this kinetic clock. Geometry, source
//! validity, temperature refresh and nonlinear population processes belong to
//! the caller. Numerical admission retains the existing 64-state bound, hence
//! at most 32 active singlet/triplet cells. Passive absorbers add no states.

use super::{
    CompensatedSum, SingletTripletDecay, SingletTripletError, finite, nonnegative, positive_product,
};
use fs_la::decay::reaction::{ReactionEdge, source::ReactionSourcePlan};
use fs_la::decay::{DecayError, MAX_DECAY_EDGES, MAX_DECAY_NODES};

/// Existing bounded positive-evolution accuracy and resource controls.
pub use fs_la::decay::DecayLimits as RecyclingLimits;

/// Positive exponential-midpoint coupling to authored triplet-pair annihilation.
#[path = "recycling/annihilation.rs"]
pub mod annihilation;
pub use annihilation::{AnnihilationLimits, AnnihilationRecyclingStep};

/// One conditional emitted photon bin from one state of one active cell.
/// State 0 is singlet and state 1 triplet. Weights sum to one per emitting state;
/// transport fractions sum to one per channel, including absorption edges.
/// A state with zero radiative rate may omit its unused spectrum entirely.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EmissionChannel {
    pub source_cell: usize,
    pub source_state: usize,
    pub photon_energy_j: f64,
    pub photon_fraction: f64,
    pub escape_fraction: f64,
    pub unresolved_fraction: f64,
}

/// A transport absorption branch for one emitted photon channel. Duplicate
/// endpoints add without renormalization; self-return is supported.
/// `excitation_fractions` are conditional probabilities per absorbed photon,
/// ordered [singlet, triplet], with sum at most one. Their complement thermalizes.
/// Each positive branch requires photon energy >= its receiver-state energy.
///
/// A passive receiver has `receiver_cell: None` and both excitation fractions
/// zero. Its external thermal-cell identity stays with the caller; absorption
/// energy and heat are returned in this original edge order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonAbsorption {
    pub channel: usize,
    pub receiver_cell: Option<usize>,
    pub absorbed_fraction: f64,
    pub excitation_fractions: [f64; 2],
}

/// Immutable source/transport topology and a reusable fixed-duration propagator.
/// Counts, pump rates and accepted physical histories are never cached here.
#[derive(Clone, Debug)]
pub struct SingletTripletRecyclingPlan {
    models: Vec<SingletTripletDecay>,
    channels: Vec<EmissionChannel>,
    absorption: Vec<PhotonAbsorption>,
    duration_s: f64,
    fraction_roundoff: f64,
    limits: RecyclingLimits,
    evolution: ReactionSourcePlan,
}

/// Whole-interval candidate, including repeated internal photon traffic.
/// Cell/state arrays are ordered [singlet, triplet]; channel and absorption
/// arrays retain their original input order. Every energy field is in joules.
#[derive(Clone, Debug, PartialEq)]
pub struct SingletTripletRecyclingStep {
    pub populations: Vec<[f64; 2]>,
    pub injected_excitations: Vec<[f64; 2]>,
    pub radiative_events: Vec<[f64; 2]>,
    pub nonradiative_events: Vec<[f64; 2]>,
    pub crossing_events: Vec<[f64; 2]>,
    pub stored_energy_j: Vec<[f64; 2]>,
    pub pumped_energy_j: Vec<[f64; 2]>,
    /// State-energy release before the photon-spectrum partition. This is
    /// diagnostic internal traffic, not an additional external-energy debit.
    pub radiative_decay_energy_j: Vec<[f64; 2]>,
    pub nonradiative_heat_j: Vec<[f64; 2]>,
    /// Source-cell Stokes heat, separated from receiver absorption heat.
    pub relaxation_heat_j: Vec<[f64; 2]>,
    pub crossing_heat_j: Vec<f64>,
    pub reverse_crossing_heat_absorbed_j: Vec<f64>,
    /// Integrated incoming re-excitations and their deposited state energy.
    /// These are internal traffic, not additional externally pumped energy.
    pub recycled_excitations: Vec<[f64; 2]>,
    pub recycled_energy_j: Vec<[f64; 2]>,
    /// Spectral photon energy per emission channel, including repeated emission.
    pub emitted_energy_j: Vec<f64>,
    pub escaped_energy_j: Vec<f64>,
    /// Unresolved transport stays unresolved; it is never assigned to heat.
    pub unresolved_energy_j: Vec<f64>,
    /// Full incident photon energy and heat, respectively, per absorption edge.
    /// For active receivers their difference is deposited excitation energy.
    pub absorbed_energy_j: Vec<f64>,
    pub absorption_heat_j: Vec<f64>,
    /// Initial store + external pump + reverse-crossing bath heat - final store
    /// - all heat - escaped photons - unresolved photons. Internal emitted,
    /// absorbed and recycled energies are deliberately excluded from this sum.
    pub residual_j: f64,
}

impl SingletTripletRecyclingPlan {
    /// Prepare all crossing and radiative feedback in one positive network.
    /// The explicit fraction allowance only admits supplied partitions; no
    /// weight, activation probability or transport fraction is normalized.
    /// Physical population and external-energy balances are checked again on
    /// every candidate, without using gross recycling traffic as error scale.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        models: &[SingletTripletDecay],
        channels: &[EmissionChannel],
        absorption: &[PhotonAbsorption],
        duration_s: f64,
        fraction_roundoff: f64,
        limits: RecyclingLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SingletTripletError> {
        poll(&mut cancelled)?;
        let cells = models.len();
        let nodes = cells.checked_mul(2).ok_or(SingletTripletError::Invalid(
            "singlet/triplet recycling state count",
        ))?;
        cap(
            "singlet/triplet recycling states",
            nodes,
            limits.maximum_nodes.min(MAX_DECAY_NODES),
        )?;
        cap(
            "spectral photon channels",
            channels.len(),
            limits.maximum_edges.min(MAX_DECAY_EDGES),
        )?;
        cap(
            "spectral absorption links",
            absorption.len(),
            limits.maximum_edges.min(MAX_DECAY_EDGES),
        )?;
        if cells == 0
            || !nonnegative(duration_s)
            || !nonnegative(fraction_roundoff)
            || fraction_roundoff >= 1.0
        {
            return Err(SingletTripletError::Invalid(
                "recycling cells, duration or fraction allowance",
            ));
        }
        let mut rates = Vec::with_capacity(nodes);
        let mut spectrum_sums = sums(nodes);
        let mut transport_sums = sums(channels.len());
        let mut feedback = sums(nodes * nodes);
        let mut has_photon_loss = vec![false; nodes];
        for (cell, model) in models.iter().enumerate() {
            poll(&mut cancelled)?;
            for state in 0..2 {
                let parent = 2 * cell + state;
                rates.push(model.total_rates_s[state]);
                feedback[(2 * cell + 1 - state) * nodes + parent].add(model.fractions[state][2])?;
            }
        }
        for (index, channel) in channels.iter().enumerate() {
            poll(&mut cancelled)?;
            if channel.source_cell >= cells
                || channel.source_state >= 2
                || !channel.photon_energy_j.is_finite()
                || channel.photon_energy_j <= 0.0
                || !fraction(channel.photon_fraction)
                || !fraction(channel.escape_fraction)
                || !fraction(channel.unresolved_fraction)
            {
                return Err(SingletTripletError::Invalid(
                    "emitted spectral photon channel",
                ));
            }
            let parent = 2 * channel.source_cell + channel.source_state;
            if channel.photon_fraction > 0.0
                && channel.photon_energy_j
                    > models[channel.source_cell].energies_j[channel.source_state]
            {
                return Err(SingletTripletError::Invalid(
                    "emitted photon exceeds its source excitation energy",
                ));
            }
            spectrum_sums[parent].add(channel.photon_fraction)?;
            transport_sums[index].add(-1.0)?;
            transport_sums[index].add(channel.escape_fraction)?;
            transport_sums[index].add(channel.unresolved_fraction)?;
            has_photon_loss[parent] |= channel.photon_fraction > 0.0
                && (channel.escape_fraction > 0.0 || channel.unresolved_fraction > 0.0);
        }
        for edge in absorption {
            poll(&mut cancelled)?;
            if edge.channel >= channels.len()
                || !fraction(edge.absorbed_fraction)
                || !edge.excitation_fractions.into_iter().all(fraction)
            {
                return Err(SingletTripletError::Invalid("photon absorption branch"));
            }
            let activation = activation_sum(edge)?;
            if activation > 1.0 {
                return Err(SingletTripletError::Invalid(
                    "one absorbed photon may create at most one excitation",
                ));
            }
            let channel = channels[edge.channel];
            let source = 2 * channel.source_cell + channel.source_state;
            transport_sums[edge.channel].add(edge.absorbed_fraction)?;
            has_photon_loss[source] |=
                channel.photon_fraction > 0.0 && edge.absorbed_fraction > 0.0 && activation < 1.0;
            match edge.receiver_cell {
                None if activation != 0.0 => {
                    return Err(SingletTripletError::Invalid(
                        "passive photon receivers cannot create excitations",
                    ));
                }
                None => {}
                Some(receiver) => {
                    if receiver >= cells {
                        return Err(SingletTripletError::Invalid("active photon receiver"));
                    }
                    for state in 0..2 {
                        let branch = edge.excitation_fractions[state];
                        if branch > 0.0
                            && channel.photon_energy_j < models[receiver].energies_j[state]
                        {
                            return Err(SingletTripletError::Invalid(
                                "absorbed photon cannot pay for the selected receiver state",
                            ));
                        }
                        // An inactive receiving branch creates no reaction
                        // coefficient. Do not underflow an unused intermediate
                        // that actual emitted photon counts could still scale.
                        if branch == 0.0
                            || edge.absorbed_fraction == 0.0
                            || channel.photon_fraction == 0.0
                            || models[channel.source_cell].radiative_rates_s[channel.source_state]
                                == 0.0
                        {
                            continue;
                        }
                        let yield_per_removal = positive_product(
                            positive_product(
                                positive_product(
                                    models[channel.source_cell].fractions[channel.source_state][0],
                                    channel.photon_fraction,
                                    "spectral radiative yield",
                                )?,
                                edge.absorbed_fraction,
                                "photon absorption yield",
                            )?,
                            branch,
                            "photon re-excitation yield",
                        )?;
                        feedback[(2 * receiver + state) * nodes + source].add(yield_per_removal)?;
                    }
                }
            }
        }
        for error in &transport_sums {
            if total(error, "photon transport partition")?.abs() > fraction_roundoff {
                return Err(SingletTripletError::Invalid(
                    "spectral absorption, escape and unresolved fractions must sum to one",
                ));
            }
        }
        let mut edges = Vec::new();
        for source in 0..nodes {
            poll(&mut cancelled)?;
            let model = models[source / 2];
            let state = source % 2;
            let weight = total(&spectrum_sums[source], "photon spectrum partition")?;
            if (model.radiative_rates_s[state] > 0.0 || weight > 0.0)
                && (weight - 1.0).abs() > fraction_roundoff
            {
                return Err(SingletTripletError::Invalid(
                    "conditional photon fractions must sum to one for every emitting state",
                ));
            }
            let mut net_yield = CompensatedSum::default();
            net_yield.add(-1.0)?;
            for receiver in 0..nodes {
                let yield_per_removal = total(
                    &feedback[receiver * nodes + source],
                    "combined excitation yield",
                )?;
                net_yield.add(yield_per_removal)?;
                if yield_per_removal > 0.0 {
                    edges.push(ReactionEdge {
                        parent: source,
                        daughter: receiver,
                        yield_per_removal,
                    });
                }
            }
            let positive_loss = model.nonradiative_rates_s[state] > 0.0
                || (model.radiative_rates_s[state] > 0.0 && has_photon_loss[source]);
            if positive_loss && total(&net_yield, "net excitation yield")? >= 0.0 {
                return Err(SingletTripletError::Unrepresentable(
                    "positive excitation loss erased by photon recycling",
                ));
            }
        }
        let evolution =
            ReactionSourcePlan::prepare(&rates, &edges, duration_s, limits, &mut cancelled)?;
        Ok(Self {
            models: models.to_vec(),
            channels: channels.to_vec(),
            absorption: absorption.to_vec(),
            duration_s,
            fraction_roundoff,
            limits,
            evolution,
        })
    }

    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub fn cell_count(&self) -> usize {
        self.models.len()
    }
    pub fn models(&self) -> &[SingletTripletDecay] {
        &self.models
    }
    pub fn channels(&self) -> &[EmissionChannel] {
        &self.channels
    }
    pub fn absorption_edges(&self) -> &[PhotonAbsorption] {
        &self.absorption
    }

    /// Reprepare this exact immutable topology and coefficients for another
    /// duration through the same numerical owner and original resource limits.
    pub fn with_duration(
        &self,
        duration_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SingletTripletError> {
        Self::prepare(
            &self.models,
            &self.channels,
            &self.absorption,
            duration_s,
            self.fraction_roundoff,
            self.limits,
            cancelled,
        )
    }

    /// Advance actual excitation populations under constant production [1/s].
    /// External pump photon absorption and its initial Stokes loss belong to
    /// the caller; these supplied counts enter already relaxed emitting states.
    pub fn advance(
        &self,
        populations: &[[f64; 2]],
        pump_rates_s: &[[f64; 2]],
    ) -> Result<SingletTripletRecyclingStep, SingletTripletError> {
        self.advance_with_cancel(populations, pump_rates_s, || false)
    }

    /// Transactional candidate with bounded cancellation checkpoints.
    pub fn advance_with_cancel(
        &self,
        populations: &[[f64; 2]],
        pump_rates_s: &[[f64; 2]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SingletTripletRecyclingStep, SingletTripletError> {
        poll(&mut cancelled)?;
        let cells = self.cell_count();
        if populations.len() != cells
            || pump_rates_s.len() != cells
            || !populations
                .iter()
                .chain(pump_rates_s)
                .flatten()
                .copied()
                .all(nonnegative)
        {
            return Err(SingletTripletError::Invalid(
                "one nonnegative singlet/triplet population and pump pair per active cell required",
            ));
        }
        let initial: Vec<_> = populations.iter().flatten().copied().collect();
        let mut injected = Vec::with_capacity(2 * cells);
        for rate in pump_rates_s.iter().flatten() {
            poll(&mut cancelled)?;
            injected.push(positive_product(
                *rate,
                self.duration_s,
                "injected excitation count",
            )?);
        }
        let result = self
            .evolution
            .apply_with_cancel(&initial, &injected, &mut cancelled)?;
        let channels = self.channels.len();
        let links = self.absorption.len();
        let mut step = SingletTripletRecyclingStep {
            populations: vec![[0.0; 2]; cells],
            injected_excitations: vec![[0.0; 2]; cells],
            radiative_events: vec![[0.0; 2]; cells],
            nonradiative_events: vec![[0.0; 2]; cells],
            crossing_events: vec![[0.0; 2]; cells],
            stored_energy_j: vec![[0.0; 2]; cells],
            pumped_energy_j: vec![[0.0; 2]; cells],
            radiative_decay_energy_j: vec![[0.0; 2]; cells],
            nonradiative_heat_j: vec![[0.0; 2]; cells],
            relaxation_heat_j: vec![[0.0; 2]; cells],
            crossing_heat_j: vec![0.0; cells],
            reverse_crossing_heat_absorbed_j: vec![0.0; cells],
            recycled_excitations: vec![[0.0; 2]; cells],
            recycled_energy_j: vec![[0.0; 2]; cells],
            emitted_energy_j: vec![0.0; channels],
            escaped_energy_j: vec![0.0; channels],
            unresolved_energy_j: vec![0.0; channels],
            absorbed_energy_j: vec![0.0; links],
            absorption_heat_j: vec![0.0; links],
            residual_j: 0.0,
        };
        let mut population_losses = CompensatedSum::default();
        for (cell, model) in self.models.iter().enumerate() {
            poll(&mut cancelled)?;
            for state in 0..2 {
                let index = 2 * cell + state;
                let energy = model.energies_j[state];
                let removals = result.removals[index];
                step.populations[cell][state] = result.populations[index];
                step.injected_excitations[cell][state] = injected[index];
                step.radiative_events[cell][state] =
                    positive_product(removals, model.fractions[state][0], "radiative events")?;
                step.nonradiative_events[cell][state] =
                    positive_product(removals, model.fractions[state][1], "nonradiative events")?;
                step.crossing_events[cell][state] =
                    positive_product(removals, model.fractions[state][2], "crossing events")?;
                step.stored_energy_j[cell][state] = positive_product(
                    result.populations[index],
                    energy,
                    "stored excitation energy",
                )?;
                step.pumped_energy_j[cell][state] =
                    positive_product(injected[index], energy, "pumped excitation energy")?;
                step.radiative_decay_energy_j[cell][state] = positive_product(
                    step.radiative_events[cell][state],
                    energy,
                    "radiative state energy",
                )?;
                step.nonradiative_heat_j[cell][state] = positive_product(
                    step.nonradiative_events[cell][state],
                    energy,
                    "nonradiative heat",
                )?;
                population_losses.add(step.nonradiative_events[cell][state])?;
            }
            let gap = model.energies_j[0] - model.energies_j[1];
            step.crossing_heat_j[cell] =
                positive_product(step.crossing_events[cell][0], gap, "crossing heat")?;
            step.reverse_crossing_heat_absorbed_j[cell] = positive_product(
                step.crossing_events[cell][1],
                gap,
                "reverse-crossing bath heat",
            )?;
        }
        let mut photon_events = vec![0.0; channels];
        let mut relaxation = sums(2 * cells);
        let mut incoming = sums(2 * cells);
        for (index, channel) in self.channels.iter().enumerate() {
            poll(&mut cancelled)?;
            let cell = channel.source_cell;
            let state = channel.source_state;
            let events = positive_product(
                step.radiative_events[cell][state],
                channel.photon_fraction,
                "spectral photons",
            )?;
            photon_events[index] = events;
            step.emitted_energy_j[index] =
                positive_product(events, channel.photon_energy_j, "emitted photon energy")?;
            if channel.photon_fraction > 0.0 {
                relaxation[2 * cell + state].add(positive_product(
                    events,
                    self.models[cell].energies_j[state] - channel.photon_energy_j,
                    "source photon relaxation heat",
                )?)?;
            }
            for (fraction, output) in [
                (channel.escape_fraction, &mut step.escaped_energy_j[index]),
                (
                    channel.unresolved_fraction,
                    &mut step.unresolved_energy_j[index],
                ),
            ] {
                let lost = positive_product(events, fraction, "outgoing photon count")?;
                population_losses.add(lost)?;
                *output =
                    positive_product(lost, channel.photon_energy_j, "outgoing photon energy")?;
            }
        }
        for (index, edge) in self.absorption.iter().enumerate() {
            poll(&mut cancelled)?;
            let photon_energy = self.channels[edge.channel].photon_energy_j;
            let events = positive_product(
                photon_events[edge.channel],
                edge.absorbed_fraction,
                "absorbed photon count",
            )?;
            step.absorbed_energy_j[index] =
                positive_product(events, photon_energy, "absorbed photon energy")?;
            let lost = positive_product(
                events,
                1.0 - activation_sum(edge)?,
                "thermalized photon count",
            )?;
            population_losses.add(lost)?;
            let mut heat = CompensatedSum::default();
            heat.add(positive_product(
                lost,
                photon_energy,
                "unconverted photon heat",
            )?)?;
            if let Some(receiver) = edge.receiver_cell {
                for state in 0..2 {
                    let activated = positive_product(
                        events,
                        edge.excitation_fractions[state],
                        "recycled excitation count",
                    )?;
                    incoming[2 * receiver + state].add(activated)?;
                    if edge.excitation_fractions[state] > 0.0 {
                        heat.add(positive_product(
                            activated,
                            photon_energy - self.models[receiver].energies_j[state],
                            "receiver photon relaxation heat",
                        )?)?;
                    }
                }
            }
            step.absorption_heat_j[index] = total(&heat, "absorbed photon heat")?;
        }
        for (cell, model) in self.models.iter().enumerate() {
            poll(&mut cancelled)?;
            for state in 0..2 {
                step.relaxation_heat_j[cell][state] = total(
                    &relaxation[2 * cell + state],
                    "source spectral relaxation heat",
                )?;
                step.recycled_excitations[cell][state] =
                    total(&incoming[2 * cell + state], "recycled excitation count")?;
                step.recycled_energy_j[cell][state] = positive_product(
                    step.recycled_excitations[cell][state],
                    model.energies_j[state],
                    "recycled state energy",
                )?;
            }
        }
        admit_population_balance(
            populations,
            &step,
            total(&population_losses, "external excitation loss")?,
            self.limits.relative_tolerance,
        )?;
        step.residual_j = admit_energy_balance(
            populations,
            &self.models,
            &step,
            self.limits.relative_tolerance,
        )?;
        poll(&mut cancelled)?;
        Ok(step)
    }
}

fn admit_population_balance(
    initial: &[[f64; 2]],
    step: &SingletTripletRecyclingStep,
    losses: f64,
    tolerance: f64,
) -> Result<(), SingletTripletError> {
    let mut balance = CompensatedSum::default();
    let mut scale = CompensatedSum::default();
    for (cell, old) in initial.iter().enumerate() {
        for state in 0..2 {
            for term in [
                step.populations[cell][state],
                -old[state],
                -step.injected_excitations[cell][state],
            ] {
                balance.add(term)?;
                scale.add(term.abs())?;
            }
        }
    }
    balance.add(losses)?;
    scale.add(losses)?;
    let residual = total(&balance, "external excitation balance")?;
    let allowance = finite(
        tolerance * total(&scale, "excitation balance scale")?,
        "excitation balance allowance",
    )?;
    if residual.abs() > allowance {
        return Err(SingletTripletError::PopulationBalance {
            residual,
            allowance,
        });
    }
    Ok(())
}

fn admit_energy_balance(
    initial: &[[f64; 2]],
    models: &[SingletTripletDecay],
    step: &SingletTripletRecyclingStep,
    tolerance: f64,
) -> Result<f64, SingletTripletError> {
    let mut balance = CompensatedSum::default();
    let mut scale = CompensatedSum::default();
    let mut thermal = CompensatedSum::default();
    for (cell, model) in models.iter().enumerate() {
        for state in 0..2 {
            for term in [
                positive_product(
                    initial[cell][state],
                    model.energies_j[state],
                    "initial state energy",
                )?,
                step.pumped_energy_j[cell][state],
                -step.stored_energy_j[cell][state],
            ] {
                balance.add(term)?;
                scale.add(term.abs())?;
            }
            thermal.add(step.nonradiative_heat_j[cell][state])?;
            thermal.add(step.relaxation_heat_j[cell][state])?;
        }
        thermal.add(step.crossing_heat_j[cell])?;
        thermal.add(-step.reverse_crossing_heat_absorbed_j[cell])?;
    }
    for &heat in &step.absorption_heat_j {
        thermal.add(heat)?;
    }
    // Bath withdrawal can fund a photon that returns the same gap as receiver
    // heat. Those repeated thermal/photon cycles must not inflate the external
    // energy allowance. Retain every local heat account above, but audit their
    // compensated NET exchange with the thermal owner only once.
    let net_heat = total(&thermal, "net thermal energy exchange")?;
    balance.add(-net_heat)?;
    scale.add(net_heat.abs())?;
    for &loss in step
        .escaped_energy_j
        .iter()
        .chain(&step.unresolved_energy_j)
    {
        balance.add(-loss)?;
        scale.add(loss)?;
    }
    let residual_j = total(&balance, "external photon energy balance")?;
    let allowance_j = finite(
        tolerance * total(&scale, "external energy balance scale")?,
        "energy balance allowance",
    )?;
    if residual_j.abs() > allowance_j {
        return Err(SingletTripletError::EnergyBalance {
            residual_j,
            allowance_j,
        });
    }
    Ok(residual_j)
}

fn sums(count: usize) -> Vec<CompensatedSum> {
    (0..count).map(|_| CompensatedSum::default()).collect()
}
fn total(sum: &CompensatedSum, quantity: &'static str) -> Result<f64, SingletTripletError> {
    finite(sum.value + sum.correction, quantity)
}
fn fraction(value: f64) -> bool {
    nonnegative(value) && value <= 1.0
}
fn activation_sum(edge: &PhotonAbsorption) -> Result<f64, SingletTripletError> {
    let mut sum = CompensatedSum::default();
    for value in edge.excitation_fractions {
        sum.add(value)?;
    }
    total(&sum, "absorbed-photon excitation branching")
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SingletTripletError> {
    if cancelled() {
        Err(DecayError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn cap(quantity: &'static str, required: usize, maximum: usize) -> Result<(), SingletTripletError> {
    if required > maximum {
        Err(DecayError::Budget {
            quantity,
            required: required as u64,
            maximum: maximum as u64,
        }
        .into())
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "recycling/tests.rs"]
mod tests;
