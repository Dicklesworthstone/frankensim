//! Complete finite-band radiation for two transmitting, isothermal specimens.
//!
//! The two selected ports face a phase-averaged far-field vacuum gap; their
//! opposite ports face independent exterior baths. The original coherent stack
//! owns each body's optics. One native internal transfer and four body/bath
//! channels apply the resolved heat once on the original thermal timeline.
//! Direct bath-to-bath transmission is diagnostic and heats neither body.

use super::*;
use crate::thermal_surface::coating_transmitting_pair::{
    CoatingTransmittingPairBudget, CoatingTransmittingPairRates, CoatingTransmittingPairResponse,
};
use fs_thermal::enthalpy::HeatTransfer;

/// One internally isothermal body. `stack` must be transmitting with vacuum on
/// both sides throughout its source band. Its selected observation direction
/// faces the gap; the original stack determines its opposite exterior face.
#[derive(Clone, Debug)]
pub struct CoatingTransmittingSurface {
    pub source_identity: ContentHash,
    pub receiver: SurfaceHeatingReceiver,
    pub stack: IsothermalCoatingStack,
}

/// Parallel, coextensive patches with unit view factor. Each body owns TWO
/// physical faces of `area_m2`, charged against its current native exposed area.
/// `baths[0]` is outside the first body, `baths[1]` outside the second. Each
/// `NativeAmbient` follows that exterior body's original boundary, even though
/// radiation from either bath may pass through and heat the other body.
/// Distinct pairs declare disjoint patches; scalar areas cannot prove this.
#[derive(Clone, Debug)]
pub struct CoatingTransmittingPair {
    pub pair_id: u64,
    pub surfaces: [CoatingTransmittingSurface; 2],
    pub area_m2: f64,
    pub baths: [CoatingExchangeBath; 2],
    pub quadrature: CoatingExchangeQuadrature,
}

#[derive(Clone, Copy, Debug)]
pub struct NativeCoatingTransmittingBudget {
    pub max_pairs: usize,
    /// Distinct original optical-law rows across all bodies and pairs.
    pub max_total_optical_knots: usize,
    /// Pair rate resolutions in one advance, including events/endpoint checks.
    pub max_source_evaluations: usize,
    /// Four original coherent solves per wavelength/cosine preparation node.
    pub max_optical_evaluations: usize,
    /// Six retained spectral channels per wavelength, summed over all pairs.
    pub max_channel_spectral_nodes: usize,
}

/// Accepted radiation already applied to the native caloric stores.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NativeCoatingTransmittingEnergy {
    pub duration_s: f64,
    pub by_pair_first_to_second_j: BTreeMap<u64, f64>,
    /// Pure exterior throughput from bath 0 to bath 1. NEVER cell heat.
    pub by_pair_bath_transfer_j: BTreeMap<u64, f64>,
    pub by_cell_internal_j: BTreeMap<CellId, f64>,
    pub by_cell_external_j: BTreeMap<CellId, f64>,
    /// Signed body/bath heat already included in native `external_heat_j`.
    pub thermal_external_j: f64,
    /// Sum of signed internal cell transfers; external heat is separate.
    pub balance_residual_j: f64,
}

#[derive(Debug)]
struct PairModel {
    pairs: Vec<CoatingTransmittingPair>,
    cells: BTreeMap<CellId, CellBinding>,
    budget: NativeCoatingTransmittingBudget,
}
#[derive(Clone, Debug)]
struct PreparedPair {
    states: [[u64; 2]; 2],
    baths: [Option<f64>; 2],
    emitters: [CoatingThermalEmission; 2],
    response: Arc<CoatingTransmittingPairResponse>,
    rates: CoatingTransmittingPairRates,
}
#[derive(Clone, Debug, Default)]
struct PairEnergySums {
    duration: CheckedSum,
    pairs: Vec<CheckedSum>,
    baths: Vec<CheckedSum>,
    internal: BTreeMap<CellId, CheckedSum>,
    external: BTreeMap<CellId, CheckedSum>,
    total_external: CheckedSum,
}

/// Complete two-body radiative coupling within the common original source band.
///
/// Use `DepositionSimulation::advance_observed`, composed with ordinary native
/// conduction and other disjoint sources. Bodies have independent T/h optical
/// laws and caloric stores. Immutable optical responses are shared across clones
/// and reused while both optical identities remain unchanged. Candidate state,
/// work and accepted heat remain private until the complete advance succeeds.
///
/// This owner is the sole radiation owner for both faces and the common band.
/// It assumes propagating vacuum modes and phase-averaged inter-body crossings,
/// with coherent optics inside each body. It infers no finite-edge view factor,
/// gap interference, evanescent transfer or internal nonisothermal film exchange.
#[derive(Clone, Debug)]
pub struct NativeCoatingTransmittingExchange {
    model: Arc<PairModel>,
    identity: ContentHash,
    time_s: f64,
    prepared_at_s: f64,
    pairs: Vec<PreparedPair>,
    sources_w: BTreeMap<CellId, f64>,
    transfers: Vec<HeatTransfer>,
    conductances_w_k: BTreeMap<CellId, f64>,
    work: NativeCoatingExchangeWork,
    sums: PairEnergySums,
    energy: NativeCoatingTransmittingEnergy,
}

impl NativeCoatingTransmittingExchange {
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        pairs: &[CoatingTransmittingPair],
        budget: NativeCoatingTransmittingBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, NativeCoatingExchangeError> {
        checkpoint()?;
        if pairs.is_empty()
            || pairs.len() > budget.max_pairs
            || pairs.len() > 4096
            || budget.max_source_evaluations == 0
        {
            return Err(NativeCoatingExchangeError::Budget);
        }
        let mut retained = 0usize;
        for pair in pairs {
            checkpoint()?;
            if pair
                .surfaces
                .iter()
                .any(|side| side.stack.films().len() > MAX_COHERENT_STACK_LAYERS)
            {
                return Err(NativeCoatingExchangeError::Budget);
            }
            retained = retained
                .checked_add(
                    pair.quadrature
                        .wavelength_bins
                        .checked_mul(6)
                        .ok_or(NativeCoatingExchangeError::Budget)?,
                )
                .ok_or(NativeCoatingExchangeError::Budget)?;
            if retained > budget.max_channel_spectral_nodes {
                return Err(NativeCoatingExchangeError::Budget);
            }
        }
        let mut ordered = pairs.to_vec();
        ordered.sort_by_key(|pair| pair.pair_id);
        let mut cells = BTreeMap::new();
        let mut laws = BTreeSet::new();
        let mut rows = 0;
        for (index, pair) in ordered.iter().enumerate() {
            checkpoint()?;
            if index > 0 && ordered[index - 1].pair_id == pair.pair_id {
                return Err(invalid(
                    "each transmitting coating pair requires a unique ID",
                ));
            }
            if !positive(pair.area_m2)
                || !positive(2.0 * pair.area_m2)
                || pair.surfaces[0].receiver.thermal_cell == pair.surfaces[1].receiver.thermal_cell
            {
                return Err(invalid(
                    "a transmitting pair needs finite positive face area and distinct receivers",
                ));
            }
            for bath in pair.baths {
                if matches!(bath, CoatingExchangeBath::Blackbody { temperature_k } if !positive(temperature_k))
                {
                    return Err(invalid(
                        "a blackbody bath requires a positive finite temperature",
                    ));
                }
            }
            for side in &pair.surfaces {
                if side.source_identity == ContentHash([0; 32]) {
                    return Err(invalid(
                        "each transmitting body requires a nonzero source identity",
                    ));
                }
                if !matches!(side.stack, IsothermalCoatingStack::Transmitting { .. }) {
                    return Err(CoatingExchangeError::VacuumTransmittingRequired.into());
                }
                let id = side.receiver.thermal_cell;
                if let std::collections::btree_map::Entry::Vacant(entry) = cells.entry(id) {
                    entry.insert(CellBinding::capture(thermal, id)?);
                }
                admit_law_rows(
                    &side.stack,
                    &mut laws,
                    &mut rows,
                    budget.max_total_optical_knots,
                )?;
            }
        }
        let mut owner = Self {
            model: Arc::new(PairModel {
                pairs: ordered,
                cells,
                budget,
            }),
            identity: ContentHash([0; 32]),
            time_s: thermal.time_s(),
            prepared_at_s: thermal.time_s(),
            pairs: Vec::new(),
            sources_w: BTreeMap::new(),
            transfers: Vec::new(),
            conductances_w_k: BTreeMap::new(),
            work: NativeCoatingExchangeWork::default(),
            sums: PairEnergySums {
                pairs: vec![CheckedSum::default(); pairs.len()],
                baths: vec![CheckedSum::default(); pairs.len()],
                ..PairEnergySums::default()
            },
            energy: NativeCoatingTransmittingEnergy::default(),
        };
        owner.refresh(thermal, &mut checkpoint)?;
        let mut hash =
            DomainHasher::new("org.frankensim.render.native-coating-transmitting-pair.v1");
        hash.update(&thermal.time_s().to_bits().to_le_bytes());
        for (&id, binding) in &owner.model.cells {
            binding.hash(id, &mut hash);
        }
        for (pair, prepared) in owner.model.pairs.iter().zip(&owner.pairs) {
            hash.update(&pair.pair_id.to_le_bytes());
            hash.update(prepared.response.identity().as_bytes());
            for (side, emitter) in pair.surfaces.iter().zip(&prepared.emitters) {
                hash.update(&side.receiver.thermal_cell.0.to_le_bytes());
                hash.update(&side.receiver.material_id.0.to_le_bytes());
                for t in [
                    side.receiver.temperature_bounds_k.0,
                    side.receiver.temperature_bounds_k.1,
                ] {
                    hash.update(&t.to_bits().to_le_bytes());
                }
                hash.update(emitter.identity().as_bytes());
                for film in side.stack.films() {
                    if let CoatingOpticsModel::Law(law) = film {
                        hash.update(law.identity().as_bytes());
                    }
                }
            }
            for (bath, resolved) in pair.baths.iter().zip(prepared.baths) {
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
    pub fn energy(&self) -> &NativeCoatingTransmittingEnergy {
        &self.energy
    }
    pub fn emitter(&self, pair_id: u64, side: usize) -> Option<&CoatingThermalEmission> {
        self.pair(pair_id)?.emitters.get(side)
    }
    pub fn rates(&self, pair_id: u64) -> Option<&CoatingTransmittingPairRates> {
        self.pair(pair_id).map(|pair| &pair.rates)
    }
    fn pair(&self, pair_id: u64) -> Option<&PreparedPair> {
        self.model
            .pairs
            .binary_search_by_key(&pair_id, |pair| pair.pair_id)
            .ok()
            .map(|index| &self.pairs[index])
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
        for pair in &self.model.pairs {
            for side in &pair.surfaces {
                check()?;
                let receiver = side.receiver;
                let (lo, hi) = receiver.temperature_bounds_k;
                if !positive(lo)
                    || !positive(hi)
                    || hi < lo
                    || self.model.cells[&receiver.thermal_cell].material != receiver.material_id
                {
                    return Err(invalid(
                        "transmitting receiver requires its original material and positive ordered temperature support",
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
                    .add(2.0 * pair.area_m2)?;
            }
        }
        for (cell, area) in occupied {
            if area.value()? > available[&cell] {
                return Err(invalid(
                    "both transmitting faces must fit current native exposed area",
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
        let mut pairs = Vec::with_capacity(self.model.pairs.len());
        let mut transfers = Vec::with_capacity(self.model.pairs.len());
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        let mut conductances = BTreeMap::<CellId, CheckedSum>::new();
        for (index, pair) in self.model.pairs.iter().enumerate() {
            check()?;
            self.work.source_evaluations = self
                .work
                .source_evaluations
                .checked_add(1)
                .ok_or(NativeCoatingExchangeError::Budget)?;
            if self.work.source_evaluations > self.model.budget.max_source_evaluations {
                return Err(NativeCoatingExchangeError::Budget);
            }
            let ids = pair
                .surfaces
                .each_ref()
                .map(|side| side.receiver.thermal_cell);
            let states = ids.map(|id| {
                [
                    temperatures[&id].to_bits(),
                    thermal.cells()[&id].specific_h_j_kg.to_bits(),
                ]
            });
            let baths = core::array::from_fn(|side| match pair.baths[side] {
                CoatingExchangeBath::Vacuum => None,
                CoatingExchangeBath::Blackbody { temperature_k } => Some(temperature_k),
                CoatingExchangeBath::NativeAmbient => {
                    Some(thermal.cells()[&ids[side]].boundary.ambient_k)
                }
            });
            let old = self.pairs.get(index);
            let mut resolve = |side: usize| {
                if let Some(old) = old.filter(|old| old.states[side] == states[side]) {
                    Ok(old.emitters[side].clone())
                } else {
                    let s = &pair.surfaces[side];
                    resolve_stack_emitter(
                        &s.stack,
                        s.receiver,
                        s.source_identity,
                        thermal,
                        self.model.budget.max_total_optical_knots,
                        check,
                    )
                }
            };
            let emitters = [resolve(0)?, resolve(1)?];
            let optical_ids = emitters
                .each_ref()
                .map(CoatingThermalEmission::optical_identity);
            let response = if let Some(old) =
                old.filter(|old| old.response.optical_identities() == optical_ids)
            {
                Arc::clone(&old.response)
            } else {
                let evaluations = pair
                    .quadrature
                    .wavelength_bins
                    .checked_mul(pair.quadrature.cosine_bins)
                    .and_then(|n| n.checked_mul(4))
                    .ok_or(NativeCoatingExchangeError::Budget)?;
                self.work.optical_evaluations = self
                    .work
                    .optical_evaluations
                    .checked_add(evaluations)
                    .ok_or(NativeCoatingExchangeError::Budget)?;
                if self.work.optical_evaluations > self.model.budget.max_optical_evaluations {
                    return Err(NativeCoatingExchangeError::Budget);
                }
                Arc::new(optical(check, |checkpoint| {
                    CoatingTransmittingPairResponse::from_emitters(
                        &emitters[0],
                        &emitters[1],
                        pair.area_m2,
                        pair.quadrature,
                        CoatingTransmittingPairBudget {
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
                    old.rates
                } else {
                    optical(check, |checkpoint| {
                        response.evaluate(
                            temperatures[&ids[0]],
                            temperatures[&ids[1]],
                            baths[0],
                            baths[1],
                            checkpoint,
                        )
                    })?
                };
            transfers.push(HeatTransfer {
                a: ids[0],
                b: ids[1],
                power_w: rates.pair.outward_power_w,
            });
            for (id, bath_rates) in ids.into_iter().zip([rates.first_baths, rates.second_baths]) {
                let source = sources.entry(id).or_default();
                let conductance = conductances.entry(id).or_default();
                conductance.add(rates.pair.secant_w_k)?;
                for rate in bath_rates {
                    source.add(-rate.outward_power_w)?;
                    conductance.add(rate.secant_w_k)?;
                }
            }
            pairs.push(PreparedPair {
                states,
                baths,
                emitters,
                response,
                rates,
            });
        }
        self.sources_w = values(&sources)?;
        self.conductances_w_k = values(&conductances)?;
        self.pairs = pairs;
        self.transfers = transfers;
        self.prepared_at_s = thermal.time_s();
        check()?;
        Ok(())
    }
}

impl DepositionObserver for NativeCoatingTransmittingExchange {
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
            return Err(invalid("transmitting pair clock mismatch"));
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
        // Native sources precede transfers on the same unmodified state. Reuse
        // that preparation rather than charging/solving all six channels twice.
        if thermal.time_s() != self.time_s || thermal.time_s() != self.prepared_at_s {
            return Err(invalid(
                "transmitting transfers require the prepared native start state",
            ));
        }
        Ok(self.transfers.clone())
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        if thermal.time_s() != self.time_s || thermal.time_s() != self.prepared_at_s {
            return Err(invalid(
                "transmitting timestep requires the prepared native start state",
            ));
        }
        Ok(thermal.limit_timestep_for_additional_conductances(
            proposed_s,
            &self.conductances_w_k,
            check,
        )?)
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if thermal.time_s() != self.time_s {
            return Err(invalid("transmitting pair event clock mismatch"));
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
                "transmitting energy requires its matching accepted native interval",
            ));
        }
        let mut sums = self.sums.clone();
        sums.duration.add(duration_s)?;
        for (index, transfer) in self.transfers.iter().enumerate() {
            check()?;
            let energy = finite(transfer.power_w * duration_s)?;
            sums.pairs[index].add(energy)?;
            sums.internal.entry(transfer.a).or_default().add(-energy)?;
            sums.internal.entry(transfer.b).or_default().add(energy)?;
            sums.baths[index].add(finite(
                self.pairs[index].rates.bath_transfer.outward_power_w * duration_s,
            )?)?;
        }
        for (&id, &watts) in &self.sources_w {
            let energy = finite(watts * duration_s)?;
            sums.external.entry(id).or_default().add(energy)?;
            sums.total_external.add(energy)?;
        }
        // Accepted-start watts own the interval; endpoint states only prepare
        // the next one. Exterior throughput never enters either cell account.
        self.refresh(thermal, check)?;
        let internal = values(&sums.internal)?;
        let mut residual = CheckedSum::default();
        for &joules in internal.values() {
            residual.add(joules)?;
        }
        let by_pair = |entries: &[CheckedSum]| {
            self.model
                .pairs
                .iter()
                .zip(entries)
                .map(|(pair, sum)| Ok((pair.pair_id, sum.value()?)))
                .collect::<Result<BTreeMap<_, _>, OpticalHeatingError>>()
        };
        let energy = NativeCoatingTransmittingEnergy {
            duration_s: sums.duration.value()?,
            by_pair_first_to_second_j: by_pair(&sums.pairs)?,
            by_pair_bath_transfer_j: by_pair(&sums.baths)?,
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

fn values(
    sums: &BTreeMap<CellId, CheckedSum>,
) -> Result<BTreeMap<CellId, f64>, OpticalHeatingError> {
    sums.iter()
        .map(|(&id, sum)| Ok((id, sum.value()?)))
        .collect()
}
