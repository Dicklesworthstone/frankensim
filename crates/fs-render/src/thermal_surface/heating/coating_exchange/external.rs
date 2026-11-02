//! Layer-resolved EXTERIOR radiation on the native thermal timeline.
//!
//! Each selected vacuum port supplies escaped thermal emission minus absorption
//! of its authored incoming bath. Films retain their own temperature, enthalpy,
//! material and caloric storage. This is a radiative boundary component: internal
//! radiative redistribution between different-temperature absorbing layers is
//! separate physics and is not inferred from exterior absorptance. Conductive
//! transfers and other supported owners compose through the native observer.
//!
//! Author both ports of a vacuum/vacuum transmitting specimen when both are
//! exposed. For one isothermal specimen, or one absorbing isothermal film with
//! lossless surroundings, those exterior ports supply its complete far-field
//! radiative boundary within the original optical band. Transmitted bath light
//! does not heat a lossless film. No near-field or bolometric claim is made.

use super::*;
use crate::thermal_surface::coating_external::{
    CoatingExternalBathBudget, CoatingExternalBathRates, CoatingExternalBathResponse,
};

/// Original ordered films, each with its own native receiver and optical law.
#[derive(Clone, Debug)]
pub enum ExternalCoatingStack {
    Opaque {
        films: Vec<StatefulCoatingLayer>,
        substrate_receiver: SurfaceHeatingReceiver,
        substrate: SubstrateOpticsModel,
    },
    /// The SELECTED port is vacuum. `opposite_ior` belongs to the other side.
    /// A second vacuum port requires an actual vacuum/vacuum specimen; flipping
    /// `reverse` alone cannot turn a dielectric exterior into vacuum.
    Transmitting {
        films: Vec<StatefulCoatingLayer>,
        opposite_ior: DielectricIor,
        reverse: bool,
    },
}
impl ExternalCoatingStack {
    fn films(&self) -> &[StatefulCoatingLayer] {
        match self {
            Self::Opaque { films, .. } | Self::Transmitting { films, .. } => films,
        }
    }
    fn substrate(&self) -> Option<(&SubstrateOpticsModel, SurfaceHeatingReceiver)> {
        match self {
            Self::Opaque {
                substrate,
                substrate_receiver,
                ..
            } => Some((substrate, *substrate_receiver)),
            Self::Transmitting { .. } => None,
        }
    }
    fn receivers(&self) -> impl Iterator<Item = SurfaceHeatingReceiver> + '_ {
        self.films()
            .iter()
            .map(|film| film.receiver)
            .chain(self.substrate().map(|(_, receiver)| receiver))
    }
}

/// One selected exterior patch, charged once against its physical surface.
///
/// `surface_cell` owns exposed area and, for `NativeAmbient`, bath temperature.
/// It may differ from the receiving film cells. Buried films need no native
/// exposed area; repeating a receiver sums its independently resolved channel
/// powers. Patches must be geometrically disjoint, or represent distinct ports
/// whose areas the surface geometry actually includes. Scalar area admission
/// cannot prove that declaration or detect a second owner of the same band.
#[derive(Clone, Debug)]
pub struct CoatingExternalBathSource {
    pub source_id: u64,
    pub source_identity: ContentHash,
    pub surface_cell: CellId,
    pub area_m2: f64,
    pub stack: ExternalCoatingStack,
    pub bath: CoatingExchangeBath,
    pub quadrature: CoatingExchangeQuadrature,
}

#[derive(Clone, Copy, Debug)]
pub struct NativeCoatingExternalBudget {
    /// Shared patch, unique original-law rows, rate and optical-work limits.
    pub exchange: NativeCoatingExchangeBudget,
    /// Sum of wavelength nodes times (films + opaque substrate) over all ports.
    pub max_channel_spectral_nodes: usize,
}

#[derive(Debug)]
struct ExternalModel {
    sources: Vec<CoatingExternalBathSource>,
    cells: BTreeMap<CellId, CellBinding>,
    budget: NativeCoatingExternalBudget,
}
#[derive(Clone, Debug)]
struct ExternalPatch {
    state: Vec<[u64; 2]>,
    bath: Option<f64>,
    emitter: CoatingThermalEmission,
    response: Arc<CoatingExternalBathResponse>,
    rates: CoatingExternalBathRates,
}

/// Composable exterior bath heat, with separate original native storage per layer.
///
/// Optical quadrature is shared across all channels of each port and cached by
/// resolved optical identity. Planck rates reuse that response as temperatures
/// change. Clones share immutable models/responses while candidate state and
/// energy accounts remain private to `advance_observed`.
#[derive(Clone, Debug)]
pub struct NativeCoatingExternalBathExchange {
    model: Arc<ExternalModel>,
    identity: ContentHash,
    time_s: f64,
    prepared_at_s: f64,
    patches: Vec<ExternalPatch>,
    sources_w: BTreeMap<CellId, f64>,
    conductances_w_k: BTreeMap<CellId, f64>,
    work: NativeCoatingExchangeWork,
    sums: EnergySums,
    cell_sums: BTreeMap<CellId, CheckedSum>,
    energy: NativeCoatingExchangeEnergy,
    energy_by_cell_j: BTreeMap<CellId, f64>,
}

impl NativeCoatingExternalBathExchange {
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        sources: &[CoatingExternalBathSource],
        budget: NativeCoatingExternalBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, NativeCoatingExchangeError> {
        checkpoint()?;
        if sources.is_empty()
            || sources.len() > budget.exchange.max_sources
            || sources.len() > 4096
            || budget.exchange.max_source_evaluations == 0
        {
            return Err(NativeCoatingExchangeError::Budget);
        }
        if sources
            .iter()
            .any(|source| source.stack.films().len() > MAX_COHERENT_STACK_LAYERS)
        {
            return Err(NativeCoatingExchangeError::Budget);
        }
        let mut ordered = sources.to_vec();
        ordered.sort_by_key(|source| source.source_id);
        let mut cells = BTreeMap::new();
        let mut laws = BTreeSet::new();
        let mut rows = 0usize;
        let mut retained = 0usize;
        for (i, source) in ordered.iter().enumerate() {
            checkpoint()?;
            if i > 0 && ordered[i - 1].source_id == source.source_id {
                return Err(invalid(
                    "each exterior coating port needs a unique source ID",
                ));
            }
            if !positive(source.area_m2) || source.source_identity == ContentHash([0; 32]) {
                return Err(invalid(
                    "positive finite exterior area and nonzero source identity required",
                ));
            }
            if matches!(source.bath, CoatingExchangeBath::Blackbody { temperature_k } if !positive(temperature_k))
            {
                return Err(invalid(
                    "a blackbody bath requires a positive finite temperature",
                ));
            }
            if source.stack.films().len() > MAX_COHERENT_STACK_LAYERS {
                return Err(NativeCoatingExchangeError::Budget);
            }
            let channels =
                source.stack.films().len() + usize::from(source.stack.substrate().is_some());
            retained = retained
                .checked_add(
                    channels
                        .checked_mul(source.quadrature.wavelength_bins)
                        .ok_or(NativeCoatingExchangeError::Budget)?,
                )
                .ok_or(NativeCoatingExchangeError::Budget)?;
            if retained > budget.max_channel_spectral_nodes {
                return Err(NativeCoatingExchangeError::Budget);
            }
            for id in core::iter::once(source.surface_cell).chain(
                source
                    .stack
                    .receivers()
                    .map(|receiver| receiver.thermal_cell),
            ) {
                if let std::collections::btree_map::Entry::Vacant(entry) = cells.entry(id) {
                    entry.insert(CellBinding::capture(thermal, id)?);
                }
            }
            admit_optical_law_rows(
                source.stack.films().iter().map(|film| &film.optics),
                source.stack.substrate().map(|(optics, _)| optics),
                &mut laws,
                &mut rows,
                budget.exchange.max_total_optical_knots,
            )?;
        }
        let mut owner = Self {
            model: Arc::new(ExternalModel {
                sources: ordered,
                cells,
                budget,
            }),
            identity: ContentHash([0; 32]),
            time_s: thermal.time_s(),
            prepared_at_s: thermal.time_s(),
            patches: Vec::new(),
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            work: NativeCoatingExchangeWork::default(),
            sums: EnergySums {
                patches: vec![CheckedSum::default(); sources.len()],
                ..EnergySums::default()
            },
            cell_sums: BTreeMap::new(),
            energy: NativeCoatingExchangeEnergy::default(),
            energy_by_cell_j: BTreeMap::new(),
        };
        owner.refresh(thermal, &mut checkpoint)?;
        let mut hash = DomainHasher::new("org.frankensim.render.native-coating-external-bath.v1");
        hash.update(&thermal.time_s().to_bits().to_le_bytes());
        for (&id, binding) in &owner.model.cells {
            binding.hash(id, &mut hash);
        }
        for (source, patch) in owner.model.sources.iter().zip(&owner.patches) {
            hash.update(&source.source_id.to_le_bytes());
            hash.update(&source.surface_cell.0.to_le_bytes());
            for receiver in source.stack.receivers() {
                hash.update(&receiver.thermal_cell.0.to_le_bytes());
                hash.update(&receiver.material_id.0.to_le_bytes());
                hash.update(&receiver.temperature_bounds_k.0.to_bits().to_le_bytes());
                hash.update(&receiver.temperature_bounds_k.1.to_bits().to_le_bytes());
            }
            hash.update(patch.response.identity().as_bytes());
            hash.update(patch.emitter.identity().as_bytes());
            match source.bath {
                CoatingExchangeBath::Vacuum => hash.update(&[0]),
                CoatingExchangeBath::Blackbody { temperature_k } => {
                    hash.update(&[1]);
                    hash.update(&temperature_k.to_bits().to_le_bytes());
                }
                CoatingExchangeBath::NativeAmbient => hash.update(&[2]),
            }
            for film in source.stack.films() {
                if let CoatingOpticsModel::Law(law) = &film.optics {
                    hash.update(law.identity().as_bytes());
                }
            }
            if let Some((SubstrateOpticsModel::Law(law), _)) = source.stack.substrate() {
                hash.update(law.identity().as_bytes());
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
    pub fn conductance_bounds_w_k(&self) -> &BTreeMap<CellId, f64> {
        &self.conductances_w_k
    }
    pub fn energy(&self) -> &NativeCoatingExchangeEnergy {
        &self.energy
    }
    /// Signed exterior heat already applied to each original native store.
    pub fn energy_by_cell_j(&self) -> &BTreeMap<CellId, f64> {
        &self.energy_by_cell_j
    }
    pub fn emitter(&self, source_id: u64) -> Option<&CoatingThermalEmission> {
        self.patch(source_id).map(|patch| &patch.emitter)
    }
    /// Latest rates preserve authored film order even for the reverse port.
    pub fn rates(&self, source_id: u64) -> Option<&CoatingExternalBathRates> {
        self.patch(source_id).map(|patch| &patch.rates)
    }
    fn patch(&self, source_id: u64) -> Option<&ExternalPatch> {
        self.model
            .sources
            .binary_search_by_key(&source_id, |source| source.source_id)
            .ok()
            .map(|i| &self.patches[i])
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, NativeCoatingExchangeError> {
        check()?;
        let mut temperatures = BTreeMap::new();
        for (&id, binding) in &self.model.cells {
            check()?;
            temperatures.insert(id, binding.temperature(thermal, id)?);
        }
        let mut occupied = BTreeMap::<CellId, CheckedSum>::new();
        for source in &self.model.sources {
            check()?;
            occupied
                .entry(source.surface_cell)
                .or_default()
                .add(source.area_m2)?;
            for receiver in source.stack.receivers() {
                let (lo, hi) = receiver.temperature_bounds_k;
                if !positive(lo)
                    || !positive(hi)
                    || lo > hi
                    || self.model.cells[&receiver.thermal_cell].material != receiver.material_id
                {
                    return Err(invalid(
                        "each exterior coating channel needs its original material and ordered temperature support",
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
            }
        }
        let selected: Vec<_> = occupied.keys().copied().collect();
        let available = thermal.exposed_areas_m2(&selected, &mut *check)?;
        for (cell, area) in occupied {
            if area.value()? > available[&cell] {
                return Err(invalid(
                    "exterior coating ports exceed the current physical exposed area",
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
        let mut patches = Vec::with_capacity(self.model.sources.len());
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        let mut conductances = BTreeMap::<CellId, CheckedSum>::new();
        for (i, source) in self.model.sources.iter().enumerate() {
            check()?;
            self.work.source_evaluations = self
                .work
                .source_evaluations
                .checked_add(1)
                .ok_or(NativeCoatingExchangeError::Budget)?;
            if self.work.source_evaluations > self.model.budget.exchange.max_source_evaluations {
                return Err(NativeCoatingExchangeError::Budget);
            }
            let state: Vec<_> = source
                .stack
                .receivers()
                .map(|receiver| {
                    let id = receiver.thermal_cell;
                    [
                        temperatures[&id].to_bits(),
                        thermal.cells()[&id].specific_h_j_kg.to_bits(),
                    ]
                })
                .collect();
            let old = self.patches.get(i);
            let same_state = old.is_some_and(|patch| patch.state == state);
            let emitter = if same_state {
                old.unwrap().emitter.clone()
            } else {
                resolve(
                    source,
                    thermal,
                    self.model.budget.exchange.max_total_optical_knots,
                    check,
                )?
            };
            let response = if let Some(old) =
                old.filter(|patch| patch.response.optical_identity() == emitter.optical_identity())
            {
                Arc::clone(&old.response)
            } else {
                let evaluations = source
                    .quadrature
                    .wavelength_bins
                    .checked_mul(source.quadrature.cosine_bins)
                    .ok_or(NativeCoatingExchangeError::Budget)?;
                self.work.optical_evaluations = self
                    .work
                    .optical_evaluations
                    .checked_add(evaluations)
                    .ok_or(NativeCoatingExchangeError::Budget)?;
                if self.work.optical_evaluations
                    > self.model.budget.exchange.max_optical_evaluations
                {
                    return Err(NativeCoatingExchangeError::Budget);
                }
                Arc::new(optical(check, |checkpoint| {
                    CoatingExternalBathResponse::from_emitter(
                        &emitter,
                        source.area_m2,
                        source.quadrature,
                        CoatingExternalBathBudget {
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
            let bath = match source.bath {
                CoatingExchangeBath::Vacuum => None,
                CoatingExchangeBath::Blackbody { temperature_k } => Some(temperature_k),
                CoatingExchangeBath::NativeAmbient => {
                    Some(thermal.cells()[&source.surface_cell].boundary.ambient_k)
                }
            };
            let rates = if let Some(old) = old.filter(|patch| same_state && patch.bath == bath) {
                old.rates.clone()
            } else {
                let film_temperatures: Vec<_> = source
                    .stack
                    .films()
                    .iter()
                    .map(|film| temperatures[&film.receiver.thermal_cell])
                    .collect();
                let substrate_temperature = source
                    .stack
                    .substrate()
                    .map(|(_, receiver)| temperatures[&receiver.thermal_cell]);
                optical(check, |checkpoint| {
                    response.evaluate(&film_temperatures, substrate_temperature, bath, checkpoint)
                })?
            };
            for (receiver, rate) in source.stack.receivers().zip(
                rates
                    .film_rates()
                    .iter()
                    .copied()
                    .chain(rates.substrate_rate()),
            ) {
                sources
                    .entry(receiver.thermal_cell)
                    .or_default()
                    .add(-rate.outward_power_w)?;
                conductances
                    .entry(receiver.thermal_cell)
                    .or_default()
                    .add(rate.secant_w_k)?;
            }
            patches.push(ExternalPatch {
                state,
                bath,
                emitter,
                response,
                rates,
            });
        }
        self.sources_w = sources
            .into_iter()
            .map(|(id, sum)| Ok((id, sum.value()?)))
            .collect::<Result<_, OpticalHeatingError>>()?;
        self.conductances_w_k = conductances
            .into_iter()
            .map(|(id, sum)| Ok((id, sum.value()?)))
            .collect::<Result<_, OpticalHeatingError>>()?;
        self.patches = patches;
        self.prepared_at_s = thermal.time_s();
        check()?;
        Ok(())
    }
}

impl DepositionObserver for NativeCoatingExternalBathExchange {
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
            return Err(invalid("exterior coating exchange clock mismatch"));
        }
        self.refresh(thermal, check)?;
        Ok(self.sources_w.clone())
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        if thermal.time_s() != self.prepared_at_s || thermal.time_s() != self.time_s {
            return Err(invalid(
                "exterior radiative timestep requires the prepared native start state",
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
            return Err(invalid("exterior coating event clock mismatch"));
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
                "exterior coating exchange requires the matching accepted native interval",
            ));
        }
        let mut sums = self.sums.clone();
        let mut cells = self.cell_sums.clone();
        sums.duration.add(duration_s)?;
        for (sum, patch) in sums.patches.iter_mut().zip(&self.patches) {
            let energy = finite(patch.rates.total_outward_power_w() * duration_s)?;
            sum.add(energy)?;
            sums.outward.add(energy)?;
        }
        for (&cell, &watts) in &self.sources_w {
            let energy = finite(watts * duration_s)?;
            sums.external.add(energy)?;
            cells.entry(cell).or_default().add(energy)?;
        }
        // Accepted-start rates own this interval; endpoint optics apply next.
        self.refresh(thermal, check)?;
        let outward = sums.outward.value()?;
        let external = sums.external.value()?;
        let energy = NativeCoatingExchangeEnergy {
            duration_s: sums.duration.value()?,
            outward_j: outward,
            thermal_external_j: external,
            balance_residual_j: finite(outward + external)?,
            by_source_outward_j: self
                .model
                .sources
                .iter()
                .zip(&sums.patches)
                .map(|(source, sum)| Ok((source.source_id, sum.value()?)))
                .collect::<Result<_, OpticalHeatingError>>()?,
        };
        let energy_by_cell_j = cells
            .iter()
            .map(|(&id, sum)| Ok((id, sum.value()?)))
            .collect::<Result<_, OpticalHeatingError>>()?;
        check()?;
        self.sums = sums;
        self.cell_sums = cells;
        self.energy = energy;
        self.energy_by_cell_j = energy_by_cell_j;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

fn resolve(
    source: &CoatingExternalBathSource,
    thermal: &DepositionSimulation<'_>,
    maximum_rows: usize,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<CoatingThermalEmission, NativeCoatingExchangeError> {
    let snapshot = match &source.stack {
        ExternalCoatingStack::Opaque {
            films,
            substrate_receiver,
            substrate,
        } => NativeCoatingEmission::try_opaque(
            thermal,
            films,
            *substrate_receiver,
            substrate,
            source.source_identity,
            maximum_rows,
            check,
        )?,
        ExternalCoatingStack::Transmitting {
            films,
            opposite_ior,
            reverse,
        } => NativeCoatingEmission::try_transmitting(
            thermal,
            films,
            *opposite_ior,
            *reverse,
            source.source_identity,
            maximum_rows,
            check,
        )?,
    };
    Ok(snapshot.emitter().clone())
}
