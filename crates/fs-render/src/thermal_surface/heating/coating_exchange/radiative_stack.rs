//! Ordered radiative shields and transmitting specimens on the native timeline.
//!
//! Every body owns its original caloric store and both physical faces. One
//! common optical solve includes screening and all repeated reflections across
//! the whole assembly. Pairwise native transfers carry internal heat; only
//! body/bath channels contribute external heat. Exterior throughput is separate.

use super::*;
use crate::thermal_surface::coating_radiative_stack::{
    CoatingRadiativeStackBudget, CoatingRadiativeStackRates, CoatingRadiativeStackResponse,
    MAX_RADIATIVE_STACK_BODIES,
};
use fs_thermal::enthalpy::HeatTransfer;

/// One internally isothermal body in an ordered planar assembly. The selected
/// port of `stack` is its LEFT face. This orientation is common to every body;
/// it does not alternate with the side of a particular gap. Both surrounding
/// half spaces must be vacuum throughout the original supported source band.
#[derive(Clone, Debug)]
pub struct CoatingRadiativeBody {
    pub source_identity: ContentHash,
    pub receiver: SurfaceHeatingReceiver,
    pub stack: IsothermalCoatingStack,
}

/// Coextensive bodies ordered physically from left to right. Order determines
/// screening and is never sorted. Each body occupies two faces of `area_m2`.
/// Distinct assemblies declare disjoint patches; scalar area cannot prove that.
/// `NativeAmbient` follows the first body's boundary on the left and the last
/// body's boundary on the right, including native boundary-change events.
#[derive(Clone, Debug)]
pub struct CoatingRadiativeStack {
    pub stack_id: u64,
    pub bodies: Vec<CoatingRadiativeBody>,
    pub area_m2: f64,
    pub baths: [CoatingExchangeBath; 2],
    pub quadrature: CoatingExchangeQuadrature,
}

#[derive(Clone, Copy, Debug)]
pub struct NativeCoatingRadiativeStackBudget {
    pub max_stacks: usize,
    /// Sum of bodies over all assemblies; each also obeys the optical hard cap.
    pub max_total_bodies: usize,
    /// Distinct original optical-law rows, shared across bodies and assemblies.
    pub max_total_optical_knots: usize,
    /// Complete assembly rate resolutions within one native advance.
    pub max_source_evaluations: usize,
    /// Two original coherent face solves per body/wavelength/cosine node.
    pub max_optical_evaluations: usize,
    /// Sum of N_lambda * (N_body + 1) * (N_body + 2) / 2 over assemblies.
    pub max_channel_spectral_nodes: usize,
}

/// Accepted heat already applied to the original native caloric stores.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NativeCoatingRadiativeStackEnergy {
    pub duration_s: f64,
    /// Keys are (assembly ID, first body index, second body index), first < second.
    /// Positive energy flows from the first body to the second body.
    pub by_stack_pair_j: BTreeMap<(u64, usize, usize), f64>,
    /// Left-to-right exterior throughput; never included in cell heating.
    pub by_stack_bath_transfer_j: BTreeMap<u64, f64>,
    pub by_cell_internal_j: BTreeMap<CellId, f64>,
    pub by_cell_external_j: BTreeMap<CellId, f64>,
    pub thermal_external_j: f64,
    /// Sum of signed internal cell transfers; external heat is separate.
    pub balance_residual_j: f64,
}

#[derive(Debug)]
struct StackModel {
    stacks: Vec<CoatingRadiativeStack>,
    cells: BTreeMap<CellId, CellBinding>,
    budget: NativeCoatingRadiativeStackBudget,
}

#[derive(Clone, Debug)]
struct PreparedStack {
    states: Vec<[u64; 2]>,
    baths: [Option<f64>; 2],
    emitters: Vec<CoatingThermalEmission>,
    response: Arc<CoatingRadiativeStackResponse>,
    rates: CoatingRadiativeStackRates,
}

#[derive(Clone, Debug, Default)]
struct StackEnergySums {
    duration: CheckedSum,
    pairs: BTreeMap<(u64, usize, usize), CheckedSum>,
    baths: BTreeMap<u64, CheckedSum>,
    internal: BTreeMap<CellId, CheckedSum>,
    external: BTreeMap<CellId, CheckedSum>,
    total_external: CheckedSum,
}

/// Native thermal evolution of an ordered assembly of coated bodies.
///
/// Use `DepositionSimulation::advance_observed`, alongside ordinary conduction
/// and disjoint heat sources. Original temperature/enthalpy laws drive each
/// body's optics. Cached immutable spectral responses are shared across clones;
/// candidate rates, work and accepted energy stay private until advance succeeds.
///
/// This owner covers both faces and the common source band of every body. It
/// assumes phase-averaged, specular, propagating vacuum modes between parallel
/// coextensive bodies. It infers no contact, finite-edge view factor, coherent
/// gap interference, near-field transport or internal nonisothermal film heat.
#[derive(Clone, Debug)]
pub struct NativeCoatingRadiativeStack {
    model: Arc<StackModel>,
    identity: ContentHash,
    time_s: f64,
    prepared_at_s: f64,
    stacks: Vec<PreparedStack>,
    sources_w: BTreeMap<CellId, f64>,
    transfers: Vec<HeatTransfer>,
    conductances_w_k: BTreeMap<CellId, f64>,
    work: NativeCoatingExchangeWork,
    sums: StackEnergySums,
    energy: NativeCoatingRadiativeStackEnergy,
}

impl NativeCoatingRadiativeStack {
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        stacks: &[CoatingRadiativeStack],
        budget: NativeCoatingRadiativeStackBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, NativeCoatingExchangeError> {
        checkpoint()?;
        if stacks.is_empty()
            || stacks.len() > budget.max_stacks
            || stacks.len() > 4096
            || stacks.len() > budget.max_source_evaluations
        {
            return Err(NativeCoatingExchangeError::Budget);
        }
        // Bound the complete retained model before copying any authored laws.
        let mut bodies = 0usize;
        let mut nodes = 0usize;
        let mut solves = 0usize;
        for stack in stacks {
            checkpoint()?;
            let count = stack.bodies.len();
            if count == 0
                || count > MAX_RADIATIVE_STACK_BODIES
                || stack
                    .bodies
                    .iter()
                    .any(|body| body.stack.films().len() > MAX_COHERENT_STACK_LAYERS)
            {
                return Err(NativeCoatingExchangeError::Budget);
            }
            bodies = bodies
                .checked_add(count)
                .ok_or(NativeCoatingExchangeError::Budget)?;
            let channels = (count + 1) * (count + 2) / 2;
            nodes = nodes
                .checked_add(
                    stack
                        .quadrature
                        .wavelength_bins
                        .checked_mul(channels)
                        .ok_or(NativeCoatingExchangeError::Budget)?,
                )
                .ok_or(NativeCoatingExchangeError::Budget)?;
            solves = solves
                .checked_add(optical_solves(stack)?)
                .ok_or(NativeCoatingExchangeError::Budget)?;
            if bodies > budget.max_total_bodies
                || nodes > budget.max_channel_spectral_nodes
                || solves > budget.max_optical_evaluations
            {
                return Err(NativeCoatingExchangeError::Budget);
            }
        }
        let mut ordered = stacks.to_vec();
        ordered.sort_by_key(|stack| stack.stack_id);
        let mut cells = BTreeMap::new();
        let mut laws = BTreeSet::new();
        let mut rows = 0;
        for (index, stack) in ordered.iter().enumerate() {
            checkpoint()?;
            if index > 0 && ordered[index - 1].stack_id == stack.stack_id {
                return Err(invalid("each radiative assembly requires a unique ID"));
            }
            if !positive(stack.area_m2) || !positive(2.0 * stack.area_m2) {
                return Err(invalid(
                    "a radiative assembly needs finite positive face area",
                ));
            }
            for bath in stack.baths {
                if matches!(bath, CoatingExchangeBath::Blackbody { temperature_k } if !positive(temperature_k))
                {
                    return Err(invalid(
                        "a blackbody bath requires positive finite temperature",
                    ));
                }
            }
            let mut receivers = BTreeSet::new();
            for body in &stack.bodies {
                checkpoint()?;
                if body.source_identity == ContentHash([0; 32]) {
                    return Err(invalid(
                        "each radiative body requires a nonzero source identity",
                    ));
                }
                if !matches!(body.stack, IsothermalCoatingStack::Transmitting { .. }) {
                    return Err(CoatingExchangeError::VacuumTransmittingRequired.into());
                }
                let id = body.receiver.thermal_cell;
                if !receivers.insert(id) {
                    return Err(invalid(
                        "radiative bodies in one assembly need distinct receivers",
                    ));
                }
                if let std::collections::btree_map::Entry::Vacant(entry) = cells.entry(id) {
                    entry.insert(CellBinding::capture(thermal, id)?);
                }
                admit_law_rows(
                    &body.stack,
                    &mut laws,
                    &mut rows,
                    budget.max_total_optical_knots,
                )?;
            }
        }
        let mut owner = Self {
            model: Arc::new(StackModel {
                stacks: ordered,
                cells,
                budget,
            }),
            identity: ContentHash([0; 32]),
            time_s: thermal.time_s(),
            prepared_at_s: thermal.time_s(),
            stacks: Vec::new(),
            sources_w: BTreeMap::new(),
            transfers: Vec::new(),
            conductances_w_k: BTreeMap::new(),
            work: NativeCoatingExchangeWork::default(),
            sums: StackEnergySums::default(),
            energy: NativeCoatingRadiativeStackEnergy::default(),
        };
        owner.refresh(thermal, &mut checkpoint)?;
        let mut hash = DomainHasher::new("org.frankensim.render.native-coating-radiative-stack.v1");
        hash.update(&thermal.time_s().to_bits().to_le_bytes());
        for (&id, binding) in &owner.model.cells {
            binding.hash(id, &mut hash);
        }
        for (stack, prepared) in owner.model.stacks.iter().zip(&owner.stacks) {
            hash.update(&stack.stack_id.to_le_bytes());
            hash.update(prepared.response.identity().as_bytes());
            for (body, emitter) in stack.bodies.iter().zip(&prepared.emitters) {
                hash.update(&body.receiver.thermal_cell.0.to_le_bytes());
                hash.update(&body.receiver.material_id.0.to_le_bytes());
                for t in [
                    body.receiver.temperature_bounds_k.0,
                    body.receiver.temperature_bounds_k.1,
                ] {
                    hash.update(&t.to_bits().to_le_bytes());
                }
                hash.update(emitter.identity().as_bytes());
                for film in body.stack.films() {
                    if let CoatingOpticsModel::Law(law) = film {
                        hash.update(law.identity().as_bytes());
                    }
                }
            }
            for (bath, resolved) in stack.baths.iter().zip(prepared.baths) {
                hash.update(&[match bath {
                    CoatingExchangeBath::Vacuum => 0,
                    CoatingExchangeBath::Blackbody { .. } => 1,
                    CoatingExchangeBath::NativeAmbient => 2,
                }]);
                if let Some(t) = resolved {
                    hash.update(&t.to_bits().to_le_bytes());
                }
            }
        }
        owner.identity = hash.finalize();
        checkpoint()?;
        Ok(owner)
    }

    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn work(&self) -> NativeCoatingExchangeWork {
        self.work
    }
    pub fn sources_w(&self) -> &BTreeMap<CellId, f64> {
        &self.sources_w
    }
    pub fn heat_transfers(&self) -> &[HeatTransfer] {
        &self.transfers
    }
    pub fn conductance_bounds_w_k(&self) -> &BTreeMap<CellId, f64> {
        &self.conductances_w_k
    }
    pub fn energy(&self) -> &NativeCoatingRadiativeStackEnergy {
        &self.energy
    }
    pub fn emitter(&self, stack_id: u64, body: usize) -> Option<&CoatingThermalEmission> {
        self.stack(stack_id)?.emitters.get(body)
    }
    pub fn rates(&self, stack_id: u64) -> Option<&CoatingRadiativeStackRates> {
        self.stack(stack_id).map(|stack| &stack.rates)
    }
    fn stack(&self, stack_id: u64) -> Option<&PreparedStack> {
        self.model
            .stacks
            .binary_search_by_key(&stack_id, |stack| stack.stack_id)
            .ok()
            .map(|index| &self.stacks[index])
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, NativeCoatingExchangeError> {
        let mut temperatures = BTreeMap::new();
        for (&id, binding) in &self.model.cells {
            check()?;
            temperatures.insert(id, binding.temperature(thermal, id)?);
        }
        let selected: Vec<_> = self.model.cells.keys().copied().collect();
        let available = thermal.exposed_areas_m2(&selected, &mut *check)?;
        let mut occupied = BTreeMap::<CellId, CheckedSum>::new();
        for stack in &self.model.stacks {
            for body in &stack.bodies {
                check()?;
                let receiver = body.receiver;
                let (lo, hi) = receiver.temperature_bounds_k;
                if !positive(lo)
                    || !positive(hi)
                    || hi < lo
                    || self.model.cells[&receiver.thermal_cell].material != receiver.material_id
                {
                    return Err(invalid(
                        "radiative receiver requires its original material and positive ordered temperature support",
                    ));
                }
                let t = temperatures[&receiver.thermal_cell];
                if !positive(t) || t < lo || t > hi {
                    return Err(NativeCoatingEmissionError::Binding(
                        SurfaceHeatingError::TemperatureDomain {
                            thermal_cell: receiver.thermal_cell,
                            temperature_k: t,
                            minimum_k: lo,
                            maximum_k: hi,
                        },
                    )
                    .into());
                }
                occupied
                    .entry(receiver.thermal_cell)
                    .or_default()
                    .add(2.0 * stack.area_m2)?;
            }
        }
        for (cell, area) in occupied {
            if area.value()? > available[&cell] {
                return Err(invalid(
                    "both faces of each radiative body must fit current exposed area",
                ));
            }
        }
        Ok(temperatures)
    }

    fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), NativeCoatingExchangeError> {
        let temperatures = self.validate(thermal, check)?;
        let mut prepared = Vec::with_capacity(self.model.stacks.len());
        let mut transfers = Vec::new();
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        let mut conductances = BTreeMap::<CellId, CheckedSum>::new();
        for (index, stack) in self.model.stacks.iter().enumerate() {
            check()?;
            self.work.source_evaluations = self
                .work
                .source_evaluations
                .checked_add(1)
                .ok_or(NativeCoatingExchangeError::Budget)?;
            if self.work.source_evaluations > self.model.budget.max_source_evaluations {
                return Err(NativeCoatingExchangeError::Budget);
            }
            let states: Vec<_> = stack
                .bodies
                .iter()
                .map(|body| {
                    let id = body.receiver.thermal_cell;
                    [
                        temperatures[&id].to_bits(),
                        thermal.cells()[&id].specific_h_j_kg.to_bits(),
                    ]
                })
                .collect();
            let baths = core::array::from_fn(|side| match stack.baths[side] {
                CoatingExchangeBath::Vacuum => None,
                CoatingExchangeBath::Blackbody { temperature_k } => Some(temperature_k),
                CoatingExchangeBath::NativeAmbient => {
                    let endpoint = if side == 0 { 0 } else { stack.bodies.len() - 1 };
                    Some(
                        thermal.cells()[&stack.bodies[endpoint].receiver.thermal_cell]
                            .boundary
                            .ambient_k,
                    )
                }
            });
            let old = self.stacks.get(index);
            let mut emitters = Vec::with_capacity(stack.bodies.len());
            for (body_index, body) in stack.bodies.iter().enumerate() {
                check()?;
                emitters.push(
                    if let Some(old) =
                        old.filter(|old| old.states[body_index] == states[body_index])
                    {
                        old.emitters[body_index].clone()
                    } else {
                        resolve_stack_emitter(
                            &body.stack,
                            body.receiver,
                            body.source_identity,
                            thermal,
                            self.model.budget.max_total_optical_knots,
                            check,
                        )?
                    },
                );
            }
            let response = if let Some(old) = old.filter(|old| {
                old.response
                    .optical_identities()
                    .iter()
                    .copied()
                    .eq(emitters
                        .iter()
                        .map(CoatingThermalEmission::optical_identity))
            }) {
                Arc::clone(&old.response)
            } else {
                let evaluations = optical_solves(stack)?;
                self.work.optical_evaluations = self
                    .work
                    .optical_evaluations
                    .checked_add(evaluations)
                    .ok_or(NativeCoatingExchangeError::Budget)?;
                if self.work.optical_evaluations > self.model.budget.max_optical_evaluations {
                    return Err(NativeCoatingExchangeError::Budget);
                }
                Arc::new(optical(check, |checkpoint| {
                    CoatingRadiativeStackResponse::from_emitters(
                        &emitters,
                        stack.area_m2,
                        stack.quadrature,
                        CoatingRadiativeStackBudget {
                            max_bodies: MAX_RADIATIVE_STACK_BODIES,
                            max_channel_spectral_nodes: self
                                .model
                                .budget
                                .max_channel_spectral_nodes,
                            max_optical_evaluations: evaluations,
                        },
                        checkpoint,
                    )
                })?)
            };
            let rates =
                if let Some(old) = old.filter(|old| old.states == states && old.baths == baths) {
                    old.rates.clone()
                } else {
                    let mut body_temperatures = [0.0; MAX_RADIATIVE_STACK_BODIES];
                    for (t, body) in body_temperatures.iter_mut().zip(&stack.bodies) {
                        *t = temperatures[&body.receiver.thermal_cell];
                    }
                    optical(check, |checkpoint| {
                        response.evaluate(
                            &body_temperatures[..stack.bodies.len()],
                            baths[0],
                            baths[1],
                            checkpoint,
                        )
                    })?
                };
            for pair in rates.body_pairs() {
                let first = stack.bodies[pair.first].receiver.thermal_cell;
                let second = stack.bodies[pair.second].receiver.thermal_cell;
                transfers.push(HeatTransfer {
                    a: first,
                    b: second,
                    power_w: pair.exchange.outward_power_w,
                });
                conductances
                    .entry(first)
                    .or_default()
                    .add(pair.exchange.secant_w_k)?;
                conductances
                    .entry(second)
                    .or_default()
                    .add(pair.exchange.secant_w_k)?;
            }
            for (body, bath_rates) in stack.bodies.iter().zip(rates.body_baths()) {
                let id = body.receiver.thermal_cell;
                for rate in bath_rates {
                    sources.entry(id).or_default().add(-rate.outward_power_w)?;
                    conductances.entry(id).or_default().add(rate.secant_w_k)?;
                }
            }
            prepared.push(PreparedStack {
                states,
                baths,
                emitters,
                response,
                rates,
            });
        }
        self.sources_w = values(&sources)?;
        self.conductances_w_k = values(&conductances)?;
        self.stacks = prepared;
        self.transfers = transfers;
        self.prepared_at_s = thermal.time_s();
        check()?;
        Ok(())
    }
}

impl DepositionObserver for NativeCoatingRadiativeStack {
    type Error = NativeCoatingExchangeError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.work = NativeCoatingExchangeWork::default();
        Ok(())
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        if thermal.time_s() != self.time_s {
            return Err(invalid("radiative assembly clock mismatch"));
        }
        self.refresh(thermal, check)?;
        Ok(self.sources_w.clone())
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        check()?;
        // Native sources precede transfers on the same unmodified start state.
        // Reusing preparation avoids resolving the entire optical network twice.
        self.require_prepared(thermal)?;
        Ok(self.transfers.clone())
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.require_prepared(thermal)?;
        Ok(thermal.limit_timestep_for_additional_conductances(
            proposed_s,
            &self.conductances_w_k,
            check,
        )?)
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if thermal.time_s() != self.time_s {
            return Err(invalid("radiative assembly event clock mismatch"));
        }
        self.refresh(thermal, &mut || Ok(()))
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if start_s != self.time_s
            || self.prepared_at_s != start_s
            || !positive(duration_s)
            || thermal.time_s() <= start_s
        {
            return Err(invalid(
                "radiative energy requires its matching accepted native interval",
            ));
        }
        let mut sums = self.sums.clone();
        sums.duration.add(duration_s)?;
        for (stack, prepared) in self.model.stacks.iter().zip(&self.stacks) {
            for pair in prepared.rates.body_pairs() {
                check()?;
                let energy = finite(pair.exchange.outward_power_w * duration_s)?;
                sums.pairs
                    .entry((stack.stack_id, pair.first, pair.second))
                    .or_default()
                    .add(energy)?;
                let first = stack.bodies[pair.first].receiver.thermal_cell;
                let second = stack.bodies[pair.second].receiver.thermal_cell;
                sums.internal.entry(first).or_default().add(-energy)?;
                sums.internal.entry(second).or_default().add(energy)?;
            }
            sums.baths.entry(stack.stack_id).or_default().add(finite(
                prepared.rates.bath_transfer.outward_power_w * duration_s,
            )?)?;
        }
        for (&id, &watts) in &self.sources_w {
            check()?;
            let energy = finite(watts * duration_s)?;
            sums.external.entry(id).or_default().add(energy)?;
            sums.total_external.add(energy)?;
        }
        // Accepted-start watts own the interval. Endpoint changes prepare the
        // next interval only; transmitted exterior light never heats a body.
        self.refresh(thermal, check)?;
        let internal = values(&sums.internal)?;
        let mut residual = CheckedSum::default();
        for &joules in internal.values() {
            residual.add(joules)?;
        }
        let energy = NativeCoatingRadiativeStackEnergy {
            duration_s: sums.duration.value()?,
            by_stack_pair_j: values(&sums.pairs)?,
            by_stack_bath_transfer_j: values(&sums.baths)?,
            by_cell_internal_j: internal,
            by_cell_external_j: values(&sums.external)?,
            thermal_external_j: sums.total_external.value()?,
            balance_residual_j: residual.value()?,
        };
        check()?;
        self.sums = sums;
        self.energy = energy;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

impl NativeCoatingRadiativeStack {
    fn require_prepared(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), NativeCoatingExchangeError> {
        if thermal.time_s() != self.time_s || thermal.time_s() != self.prepared_at_s {
            Err(invalid(
                "radiative assembly requires the prepared native start state",
            ))
        } else {
            Ok(())
        }
    }
}

fn optical_solves(stack: &CoatingRadiativeStack) -> Result<usize, NativeCoatingExchangeError> {
    stack
        .quadrature
        .wavelength_bins
        .checked_mul(stack.quadrature.cosine_bins)
        .and_then(|n| n.checked_mul(stack.bodies.len()))
        .and_then(|n| n.checked_mul(2))
        .ok_or(NativeCoatingExchangeError::Budget)
}

fn values<K: Ord + Copy>(
    sums: &BTreeMap<K, CheckedSum>,
) -> Result<BTreeMap<K, f64>, OpticalHeatingError> {
    sums.iter()
        .map(|(&key, sum)| Ok((key, sum.value()?)))
        .collect()
}
