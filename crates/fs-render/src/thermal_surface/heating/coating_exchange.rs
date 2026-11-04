//! Coherent finite-band coating exchange on the original thermal timeline.
//!
//! One receiver is the complete internally isothermal film/substrate specimen.
//! Its authored mass and caloric law own that lumped storage; separate film
//! temperatures require internal radiative exchange and are not inferred here.
//! Each patch exchanges with its selected exterior vacuum port. Other ports or
//! spectral bands need separately composed owners. This is not bolometric heat.

/// Conservative mutual exchange between opaque native coating specimens.
pub mod gap;

/// Layer-resolved exchange with explicitly authored exterior radiation baths.
pub mod external;

/// Mutual radiation and two exterior baths for vacuum-facing transmitting bodies.
pub mod transmitting_pair;

/// Mutual screening and heat exchange through an ordered stack of coated bodies.
pub mod radiative_stack;

use super::{
    CoatingOpticsModel, NativeCoatingEmission, NativeCoatingEmissionError, StatefulCoatingLayer,
    SubstrateOpticsModel, SurfaceHeatingError, SurfaceHeatingReceiver,
};
use crate::conductor::MAX_COHERENT_STACK_LAYERS;
use crate::dielectric::DielectricIor;
use crate::thermal_surface::CoatingThermalEmission;
use crate::thermal_surface::coating_exchange::{
    CoatingExchangeBudget, CoatingExchangeError, CoatingExchangeQuadrature, CoatingExchangeResponse,
};
use crate::volumes::heating::{CheckedSum, OpticalHeatingError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
    MaterialId, ThermalStorageBasis,
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// Incident thermal radiation in the selected vacuum half space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CoatingExchangeBath {
    Vacuum,
    Blackbody {
        temperature_k: f64,
    },
    /// Follow the original cell boundary's ambient temperature, including events.
    NativeAmbient,
}

/// Optical films share one native caloric temperature and energy coordinate.
#[derive(Clone, Debug)]
pub enum IsothermalCoatingStack {
    Opaque {
        films: Vec<CoatingOpticsModel>,
        substrate: SubstrateOpticsModel,
    },
    Transmitting {
        films: Vec<CoatingOpticsModel>,
        opposite_ior: DielectricIor,
        reverse: bool,
    },
}
impl IsothermalCoatingStack {
    fn films(&self) -> &[CoatingOpticsModel] {
        match self {
            Self::Opaque { films, .. } | Self::Transmitting { films, .. } => films,
        }
    }
}

/// One physically disjoint exposed patch, with its own selected emitting port.
/// All patches assigned to a cell must fit its CURRENT native exposed area.
/// Ambient convection may coexist; never also attach another radiation owner
/// to the same patch and wavelength band. Scalar area admission cannot prove
/// geometric non-overlap or distinguish the two faces of a transmitting film.
#[derive(Clone, Debug)]
pub struct CoatingExchangeSource {
    pub source_id: u64,
    pub source_identity: ContentHash,
    pub receiver: SurfaceHeatingReceiver,
    pub area_m2: f64,
    pub stack: IsothermalCoatingStack,
    pub bath: CoatingExchangeBath,
    pub quadrature: CoatingExchangeQuadrature,
}

#[derive(Clone, Copy, Debug)]
pub struct NativeCoatingExchangeBudget {
    pub max_sources: usize,
    /// Distinct film/substrate law rows or retained material property queries
    /// across all patches. A shared law is counted once.
    pub max_total_optical_knots: usize,
    /// Rate evaluations across a complete advance, including endpoint checks.
    pub max_source_evaluations: usize,
    /// Coherent angle/wavelength solves across a complete advance.
    pub max_optical_evaluations: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NativeCoatingExchangeWork {
    pub source_evaluations: usize,
    pub optical_evaluations: usize,
}

/// Cumulative accepted energy, additional to no other thermal store or audit.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NativeCoatingExchangeEnergy {
    pub duration_s: f64,
    /// Positive when radiation leaves the specimen; negative for bath heating.
    pub outward_j: f64,
    /// Signed source energy already included in native external_heat_j.
    pub thermal_external_j: f64,
    pub balance_residual_j: f64,
    pub by_source_outward_j: BTreeMap<u64, f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum NativeCoatingExchangeError {
    Invalid(&'static str),
    Budget,
    Binding(NativeCoatingEmissionError),
    Exchange(CoatingExchangeError),
    Thermal(DepositionError),
    Accounting(OpticalHeatingError),
}
impl From<NativeCoatingEmissionError> for NativeCoatingExchangeError {
    fn from(value: NativeCoatingEmissionError) -> Self {
        Self::Binding(value)
    }
}
impl From<CoatingExchangeError> for NativeCoatingExchangeError {
    fn from(value: CoatingExchangeError) -> Self {
        Self::Exchange(value)
    }
}
impl From<DepositionError> for NativeCoatingExchangeError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}
impl From<OpticalHeatingError> for NativeCoatingExchangeError {
    fn from(value: OpticalHeatingError) -> Self {
        Self::Accounting(value)
    }
}
impl core::fmt::Display for NativeCoatingExchangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native coating exchange refused: {self:?}")
    }
}
impl std::error::Error for NativeCoatingExchangeError {}

#[derive(Clone, Debug)]
struct CellBinding {
    material: MaterialId,
    mass_kg: f64,
    geometry: CellGeometry,
    caloric: Option<[u8; 32]>,
    storage: ThermalStorageBasis,
}
impl CellBinding {
    fn capture(
        thermal: &DepositionSimulation<'_>,
        id: CellId,
    ) -> Result<Self, NativeCoatingExchangeError> {
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or_else(|| invalid("radiative receiver must already be active"))?;
        Ok(Self {
            material: cell.material_id,
            mass_kg: cell.mass_kg,
            geometry: cell.geometry,
            caloric: thermal.material_caloric_identity(cell.material_id)?,
            storage: thermal.material_storage_basis(cell.material_id)?,
        })
    }

    fn temperature(
        &self,
        thermal: &DepositionSimulation<'_>,
        id: CellId,
    ) -> Result<f64, NativeCoatingExchangeError> {
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or_else(|| invalid("radiative receiver is inactive"))?;
        if cell.material_id != self.material
            || cell.mass_kg != self.mass_kg
            || cell.geometry != self.geometry
            || thermal.material_caloric_identity(cell.material_id)? != self.caloric
            || thermal.material_storage_basis(cell.material_id)? != self.storage
        {
            return Err(invalid(
                "radiative receiver must retain its native material, caloric law, mass and geometry",
            ));
        }
        Ok(thermal.temperature_k(id)?)
    }

    fn hash(&self, id: CellId, hash: &mut DomainHasher) {
        hash.update(&id.0.to_le_bytes());
        hash.update(&self.material.0.to_le_bytes());
        for value in [
            self.mass_kg,
            self.geometry.volume_m3,
            self.geometry.surface_area_m2,
        ]
        .into_iter()
        .chain(self.geometry.start_m)
        .chain(self.geometry.end_m)
        {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.update(&[u8::from(
            self.storage == ThermalStorageBasis::InternalEnergy,
        )]);
        hash.update(&[u8::from(self.caloric.is_some())]);
        if let Some(caloric) = self.caloric {
            hash.update(&caloric);
        }
    }
}
#[derive(Debug)]
struct Model {
    sources: Vec<CoatingExchangeSource>,
    cells: BTreeMap<CellId, CellBinding>,
    budget: NativeCoatingExchangeBudget,
}
#[derive(Clone, Debug)]
struct PreparedPatch {
    state: [u64; 2],
    emitter: CoatingThermalEmission,
    response: Arc<CoatingExchangeResponse>,
    outward_w: f64,
}
#[derive(Clone, Debug, Default)]
struct EnergySums {
    duration: CheckedSum,
    outward: CheckedSum,
    external: CheckedSum,
    patches: Vec<CheckedSum>,
}

/// Native signed heat sources and their composed radiative step restriction.
///
/// Clones share immutable laws and optical quadrature responses; candidate
/// rates, caches and energy sums remain private to the native transactional
/// advance. Use with `DepositionSimulation::advance_observed`, or compose with
/// existing optical illumination through `fs_thermal::observers`.
#[derive(Clone, Debug)]
pub struct NativeCoatingExchange {
    model: Arc<Model>,
    identity: ContentHash,
    time_s: f64,
    prepared_at_s: f64,
    patches: Vec<PreparedPatch>,
    sources_w: BTreeMap<CellId, f64>,
    conductances_w_k: BTreeMap<CellId, f64>,
    work: NativeCoatingExchangeWork,
    sums: EnergySums,
    energy: NativeCoatingExchangeEnergy,
}

impl NativeCoatingExchange {
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        sources: &[CoatingExchangeSource],
        budget: NativeCoatingExchangeBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, NativeCoatingExchangeError> {
        checkpoint()?;
        if sources.is_empty()
            || sources.len() > budget.max_sources
            || sources.len() > 4096
            || budget.max_source_evaluations == 0
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
        for (i, source) in ordered.iter().enumerate() {
            checkpoint()?;
            if i > 0 && ordered[i - 1].source_id == source.source_id {
                return Err(invalid(
                    "each physical radiative patch needs a unique source ID",
                ));
            }
            if !positive(source.area_m2) || source.source_identity == ContentHash([0; 32]) {
                return Err(invalid(
                    "positive finite patch area and nonzero source identity required",
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
            let id = source.receiver.thermal_cell;
            cells
                .entry(id)
                .or_insert(CellBinding::capture(thermal, id)?);
            admit_law_rows(
                &source.stack,
                &mut laws,
                &mut rows,
                budget.max_total_optical_knots,
            )?;
        }
        let mut owner = Self {
            model: Arc::new(Model {
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
            energy: NativeCoatingExchangeEnergy::default(),
        };
        owner.refresh(thermal, &mut checkpoint)?;
        let mut hash = DomainHasher::new("org.frankensim.render.native-coating-exchange.v1");
        hash.update(&thermal.time_s().to_bits().to_le_bytes());
        for (&id, binding) in &owner.model.cells {
            binding.hash(id, &mut hash);
        }
        for (source, patch) in owner.model.sources.iter().zip(&owner.patches) {
            for id in [
                source.source_id,
                source.receiver.thermal_cell.0,
                source.receiver.material_id.0,
            ] {
                hash.update(&id.to_le_bytes());
            }
            for x in [
                source.receiver.temperature_bounds_k.0,
                source.receiver.temperature_bounds_k.1,
            ] {
                hash.update(&x.to_bits().to_le_bytes());
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
                if let CoatingOpticsModel::Law(law) = film {
                    hash.update(law.identity().as_bytes());
                }
            }
            if let IsothermalCoatingStack::Opaque {
                substrate: SubstrateOpticsModel::Law(law),
                ..
            } = &source.stack
            {
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
    /// Latest resolved optical/thermal state, suitable for directional observation.
    pub fn emitter(&self, source_id: u64) -> Option<&CoatingThermalEmission> {
        self.model
            .sources
            .binary_search_by_key(&source_id, |source| source.source_id)
            .ok()
            .map(|i| &self.patches[i].emitter)
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
        let selected: Vec<_> = self.model.cells.keys().copied().collect();
        let available = thermal.exposed_areas_m2(&selected, &mut *check)?;
        let mut occupied = BTreeMap::<CellId, CheckedSum>::new();
        for source in &self.model.sources {
            check()?;
            let receiver = source.receiver;
            let (lo, hi) = receiver.temperature_bounds_k;
            if !positive(lo)
                || !positive(hi)
                || lo > hi
                || self.model.cells[&receiver.thermal_cell].material != receiver.material_id
            {
                return Err(invalid(
                    "receiver requires the original material and positive ordered temperature support",
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
                .add(source.area_m2)?;
        }
        for (cell, area) in occupied {
            if area.value()? > available[&cell] {
                return Err(invalid(
                    "radiative patches exceed the current native exposed area",
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
            if self.work.source_evaluations > self.model.budget.max_source_evaluations {
                return Err(NativeCoatingExchangeError::Budget);
            }
            let cell = source.receiver.thermal_cell;
            let t = temperatures[&cell];
            let state = [
                t.to_bits(),
                thermal.cells()[&cell].specific_h_j_kg.to_bits(),
            ];
            let old = self.patches.get(i);
            let emitter = if let Some(old) = old.filter(|old| old.state == state) {
                old.emitter.clone()
            } else {
                resolve_emitter(
                    source,
                    thermal,
                    self.model.budget.max_total_optical_knots,
                    check,
                )?
            };
            let response = if let Some(old) =
                old.filter(|old| old.response.optical_identity() == emitter.optical_identity())
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
                if self.work.optical_evaluations > self.model.budget.max_optical_evaluations {
                    return Err(NativeCoatingExchangeError::Budget);
                }
                Arc::new(optical(check, |checkpoint| {
                    CoatingExchangeResponse::from_emitter(
                        &emitter,
                        source.area_m2,
                        source.quadrature,
                        CoatingExchangeBudget {
                            max_spectral_nodes: source.quadrature.wavelength_bins,
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
                    Some(thermal.cells()[&cell].boundary.ambient_k)
                }
            };
            let rates = optical(check, |checkpoint| response.evaluate(t, bath, checkpoint))?;
            sources
                .entry(cell)
                .or_default()
                .add(-rates.outward_power_w)?;
            conductances
                .entry(cell)
                .or_default()
                .add(rates.secant_w_k)?;
            patches.push(PreparedPatch {
                state,
                emitter,
                response,
                outward_w: rates.outward_power_w,
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

impl DepositionObserver for NativeCoatingExchange {
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
            return Err(invalid("coating exchange clock mismatch"));
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
                "radiative timestep requires the prepared native start state",
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
            return Err(invalid("coating exchange event clock mismatch"));
        }
        // Events can cover an emitting patch or change its native ambient bath.
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
                "coating exchange requires the matching accepted native interval",
            ));
        }
        let mut sums = self.sums.clone();
        sums.duration.add(duration_s)?;
        for (sum, patch) in sums.patches.iter_mut().zip(&self.patches) {
            let energy = finite(patch.outward_w * duration_s)?;
            sum.add(energy)?;
            sums.outward.add(energy)?;
        }
        for &watts in self.sources_w.values() {
            sums.external.add(finite(watts * duration_s)?)?;
        }
        // Charge the old accepted-start rates BEFORE refreshing endpoint optics.
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
        check()?;
        self.sums = sums;
        self.energy = energy;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

fn resolve_emitter(
    source: &CoatingExchangeSource,
    thermal: &DepositionSimulation<'_>,
    maximum_rows: usize,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<CoatingThermalEmission, NativeCoatingExchangeError> {
    resolve_stack_emitter(
        &source.stack,
        source.receiver,
        source.source_identity,
        thermal,
        maximum_rows,
        check,
    )
}

fn resolve_stack_emitter(
    stack: &IsothermalCoatingStack,
    receiver: SurfaceHeatingReceiver,
    source_identity: ContentHash,
    thermal: &DepositionSimulation<'_>,
    maximum_rows: usize,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<CoatingThermalEmission, NativeCoatingExchangeError> {
    let films: Vec<_> = stack
        .films()
        .iter()
        .map(|optics| StatefulCoatingLayer {
            receiver,
            optics: optics.clone(),
        })
        .collect();
    let snapshot = match stack {
        IsothermalCoatingStack::Opaque { substrate, .. } => NativeCoatingEmission::try_opaque(
            thermal,
            &films,
            receiver,
            substrate,
            source_identity,
            maximum_rows,
            check,
        )?,
        IsothermalCoatingStack::Transmitting {
            opposite_ior,
            reverse,
            ..
        } => NativeCoatingEmission::try_transmitting(
            thermal,
            &films,
            *opposite_ior,
            *reverse,
            source_identity,
            maximum_rows,
            check,
        )?,
    };
    Ok(snapshot.emitter().clone())
}

fn admit_law_rows(
    stack: &IsothermalCoatingStack,
    laws: &mut BTreeSet<(u8, [u8; 32])>,
    rows: &mut usize,
    maximum: usize,
) -> Result<(), NativeCoatingExchangeError> {
    admit_optical_law_rows(
        stack.films().iter(),
        match stack {
            IsothermalCoatingStack::Opaque { substrate, .. } => Some(substrate),
            IsothermalCoatingStack::Transmitting { .. } => None,
        },
        laws,
        rows,
        maximum,
    )
}

fn admit_optical_law_rows<'a>(
    films: impl IntoIterator<Item = &'a CoatingOpticsModel>,
    substrate: Option<&SubstrateOpticsModel>,
    laws: &mut BTreeSet<(u8, [u8; 32])>,
    rows: &mut usize,
    maximum: usize,
) -> Result<(), NativeCoatingExchangeError> {
    for film in films {
        if let CoatingOpticsModel::Law(law) = film {
            if laws.insert((0, *law.identity().as_bytes())) {
                *rows = rows
                    .checked_add(law.row_count())
                    .ok_or(NativeCoatingExchangeError::Budget)?;
            }
        }
    }
    if let Some(SubstrateOpticsModel::Law(law)) = substrate {
        if laws.insert((1, *law.identity().as_bytes())) {
            *rows = rows
                .checked_add(law.budget_size())
                .ok_or(NativeCoatingExchangeError::Budget)?;
        }
    }
    if *rows > maximum {
        return Err(NativeCoatingExchangeError::Budget);
    }
    Ok(())
}

// Preserve an original checkpoint error across the optical owner's own error type.
fn optical<T>(
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(
        &mut dyn FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<T, CoatingExchangeError>,
) -> Result<T, NativeCoatingExchangeError> {
    let mut refusal = None;
    let result = operation(&mut || match check() {
        Ok(()) => Ok(()),
        Err(error) => {
            refusal = Some(error);
            Err(CoatingExchangeError::Cancelled)
        }
    });
    if let Some(error) = refusal {
        Err(error.into())
    } else {
        Ok(result?)
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn invalid(reason: &'static str) -> NativeCoatingExchangeError {
    NativeCoatingExchangeError::Invalid(reason)
}
fn finite(value: f64) -> Result<f64, NativeCoatingExchangeError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable radiative heat account"))
    }
}
