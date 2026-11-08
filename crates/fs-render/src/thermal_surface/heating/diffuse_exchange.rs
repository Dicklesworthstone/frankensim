//! Spectral diffuse emission and absorption on the original thermal timeline.
//!
//! The accepted native T/h resolves one [`SurfaceDiffuseLaw`] for BOTH cooling
//! and observed light. Each physically disjoint exposed patch receives its own
//! isotropic blackbody bath or vacuum. Do not assign another radiation owner to
//! that same patch and wavelength band. Native caloric storage, conduction,
//! convection, source events and transactional acceptance remain authoritative.

use super::diffuse_feedback::SurfaceDiffuseLaw;
use super::{SurfaceHeatingError, SurfaceHeatingReceiver};
use crate::thermal_surface::SpectralDiffuseThermalSurface;
use crate::thermal_surface::diffuse_exchange::{
    DiffuseExchangeBudget, DiffuseExchangeError, DiffuseExchangeQuadrature, DiffuseExchangeRates,
    DiffuseExchangeResponse,
};
use crate::volumes::heating::{CheckedSum, OpticalHeatingError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
    MaterialId, ThermalStorageBasis,
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// Incoming radiation fills this patch's outgoing vacuum hemisphere.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DiffuseExchangeBath {
    Vacuum,
    Blackbody {
        temperature_k: f64,
    },
    /// Follow the original cell's ambient temperature, including native events.
    NativeAmbient,
}

#[derive(Clone, Debug)]
pub struct DiffuseExchangeSource {
    pub source_id: u64,
    pub source_identity: ContentHash,
    pub receiver: SurfaceHeatingReceiver,
    /// One physical emitting-side area, not volume or projected camera area.
    pub area_m2: f64,
    pub law: Arc<SurfaceDiffuseLaw>,
    pub bath: DiffuseExchangeBath,
    pub quadrature: DiffuseExchangeQuadrature,
    /// Optional absolute bound on unmodeled-wavelength exchange, in W.
    /// None explicitly permits a partial-band calculation; the omitted bound
    /// remains reported. This allowance never admits quadrature/time error.
    pub maximum_omitted_exchange_w: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct NativeDiffuseExchangeBudget {
    pub max_sources: usize,
    pub max_total_law_rows: usize,
    pub max_total_spectral_nodes: usize,
    /// Rate evaluations over the whole candidate advance, including endpoints.
    pub max_source_evaluations: usize,
    pub max_spectral_evaluations: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NativeDiffuseExchangeWork {
    pub source_evaluations: usize,
    pub spectral_evaluations: usize,
    /// Quadrature grids prepared; state-law changes reuse the original grid.
    pub response_preparations: usize,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct NativeDiffuseExchangeEnergy {
    pub duration_s: f64,
    /// Signed outward heat; negative for net bath heating.
    pub outward_j: f64,
    /// Already included once in native external_heat_j; do not add again.
    pub thermal_external_j: f64,
    pub balance_residual_j: f64,
    /// Sum of omitted-band bounds at the same frozen start states and durations
    /// used for applied heat. This is not a continuous-time error certificate.
    pub omitted_exchange_upper_bound_j: f64,
    pub by_source_outward_j: BTreeMap<u64, f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum NativeDiffuseExchangeError {
    Invalid(&'static str),
    Budget,
    Coverage {
        source_id: u64,
        upper_bound_w: f64,
        maximum_w: f64,
    },
    Binding(SurfaceHeatingError),
    Exchange(DiffuseExchangeError),
    Thermal(DepositionError),
    Accounting(OpticalHeatingError),
}
impl From<SurfaceHeatingError> for NativeDiffuseExchangeError {
    fn from(value: SurfaceHeatingError) -> Self {
        Self::Binding(value)
    }
}
impl From<DiffuseExchangeError> for NativeDiffuseExchangeError {
    fn from(value: DiffuseExchangeError) -> Self {
        Self::Exchange(value)
    }
}
impl From<DepositionError> for NativeDiffuseExchangeError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}
impl From<OpticalHeatingError> for NativeDiffuseExchangeError {
    fn from(value: OpticalHeatingError) -> Self {
        Self::Accounting(value)
    }
}
impl core::fmt::Display for NativeDiffuseExchangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native diffuse exchange refused: {self:?}")
    }
}
impl std::error::Error for NativeDiffuseExchangeError {}

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
    ) -> Result<Self, NativeDiffuseExchangeError> {
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or_else(|| invalid("diffuse radiative receiver must be active"))?;
        Ok(Self {
            material: cell.material_id,
            mass_kg: cell.mass_kg,
            geometry: cell.geometry,
            caloric: thermal.material_caloric_identity(cell.material_id)?,
            storage: thermal.material_storage_basis(cell.material_id)?,
        })
    }
    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        id: CellId,
    ) -> Result<(), NativeDiffuseExchangeError> {
        let actual = Self::capture(thermal, id)?;
        if actual.material != self.material
            || actual.mass_kg != self.mass_kg
            || actual.geometry != self.geometry
            || actual.caloric != self.caloric
            || actual.storage != self.storage
        {
            return Err(invalid(
                "diffuse receiver must retain its native material, caloric law, mass and geometry",
            ));
        }
        Ok(())
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
        hash.update(&[match self.storage {
            ThermalStorageBasis::Enthalpy => 0,
            ThermalStorageBasis::InternalEnergy => 1,
            ThermalStorageBasis::DeclaredSensibleStorage => 2,
        }]);
        hash.update(&[u8::from(self.caloric.is_some())]);
        if let Some(caloric) = self.caloric {
            hash.update(&caloric);
        }
    }
}

#[derive(Debug)]
struct Model {
    sources: Vec<DiffuseExchangeSource>,
    cells: BTreeMap<CellId, CellBinding>,
    budget: NativeDiffuseExchangeBudget,
}
#[derive(Clone, Debug)]
struct Patch {
    emitter: SpectralDiffuseThermalSurface,
    response: DiffuseExchangeResponse,
    rates: DiffuseExchangeRates,
}
#[derive(Clone, Debug, Default)]
struct EnergySums {
    duration: CheckedSum,
    outward: CheckedSum,
    external: CheckedSum,
    omitted: CheckedSum,
    patches: Vec<CheckedSum>,
}

/// Conservative source adapter for `DepositionSimulation::advance_observed`.
///
/// The native explicit owner applies accepted-start rates and alone advances
/// enthalpy, including latent storage. Frozen-emissivity conductances compose
/// with native step restrictions; they do not bound a state law's derivative.
/// Refine max_dt independently of wavelength quadrature. Candidate clones share
/// only immutable laws/node grids, so cancellation publishes no partial heat.
#[derive(Clone, Debug)]
pub struct NativeDiffuseExchange {
    model: Arc<Model>,
    identity: ContentHash,
    time_s: f64,
    prepared_at_s: f64,
    patches: Vec<Patch>,
    sources_w: BTreeMap<CellId, f64>,
    conductances_w_k: BTreeMap<CellId, f64>,
    work: NativeDiffuseExchangeWork,
    sums: EnergySums,
    energy: NativeDiffuseExchangeEnergy,
}

impl NativeDiffuseExchange {
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        sources: &[DiffuseExchangeSource],
        budget: NativeDiffuseExchangeBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, NativeDiffuseExchangeError> {
        checkpoint()?;
        if sources.is_empty()
            || sources.len() > budget.max_sources
            || sources.len() > 4096
            || budget.max_source_evaluations == 0
            || budget.max_spectral_evaluations == 0
        {
            return Err(NativeDiffuseExchangeError::Budget);
        }
        let mut ordered = sources.to_vec();
        ordered.sort_by_key(|source| source.source_id);
        let mut cells = BTreeMap::new();
        let mut laws = BTreeSet::new();
        let mut rows = 0usize;
        for (index, source) in ordered.iter().enumerate() {
            checkpoint()?;
            if index > 0 && ordered[index - 1].source_id == source.source_id {
                return Err(invalid(
                    "each physical diffuse patch needs a unique source ID",
                ));
            }
            if !positive(source.area_m2)
                || source.source_identity == ContentHash([0; 32])
                || source
                    .maximum_omitted_exchange_w
                    .is_some_and(|limit| !limit.is_finite() || limit < 0.0)
                || matches!(source.bath, DiffuseExchangeBath::Blackbody { temperature_k } if !positive(temperature_k))
            {
                return Err(invalid(
                    "positive finite area/bath, source identity and nonnegative coverage allowance required",
                ));
            }
            if laws.insert(source.law.identity()) {
                rows = rows
                    .checked_add(source.law.row_count())
                    .ok_or(NativeDiffuseExchangeError::Budget)?;
                if rows > budget.max_total_law_rows {
                    return Err(NativeDiffuseExchangeError::Budget);
                }
            }
            let id = source.receiver.thermal_cell;
            cells
                .entry(id)
                .or_insert(CellBinding::capture(thermal, id)?);
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
            work: NativeDiffuseExchangeWork::default(),
            sums: EnergySums {
                patches: vec![CheckedSum::default(); sources.len()],
                ..EnergySums::default()
            },
            energy: NativeDiffuseExchangeEnergy::default(),
        };
        owner.refresh(thermal, &mut checkpoint)?;
        let mut hash = DomainHasher::new("org.frankensim.render.native-diffuse-exchange.v1");
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
            hash.update(source.source_identity.as_bytes());
            hash.update(source.law.identity().as_bytes());
            hash.update(patch.response.identity().as_bytes());
            hash.update(patch.emitter.identity().as_bytes());
            for value in [
                source.receiver.temperature_bounds_k.0,
                source.receiver.temperature_bounds_k.1,
            ] {
                hash.update(&value.to_bits().to_le_bytes());
            }
            match source.bath {
                DiffuseExchangeBath::Vacuum => hash.update(&[0]),
                DiffuseExchangeBath::Blackbody { temperature_k } => {
                    hash.update(&[1]);
                    hash.update(&temperature_k.to_bits().to_le_bytes());
                }
                DiffuseExchangeBath::NativeAmbient => hash.update(&[2]),
            }
            hash.update(&[u8::from(source.maximum_omitted_exchange_w.is_some())]);
            if let Some(limit) = source.maximum_omitted_exchange_w {
                hash.update(&limit.to_bits().to_le_bytes());
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
    pub const fn work(&self) -> NativeDiffuseExchangeWork {
        self.work
    }
    pub fn sources_w(&self) -> &BTreeMap<CellId, f64> {
        &self.sources_w
    }
    pub fn conductance_bounds_w_k(&self) -> &BTreeMap<CellId, f64> {
        &self.conductances_w_k
    }
    pub fn energy(&self) -> &NativeDiffuseExchangeEnergy {
        &self.energy
    }
    /// Accepted native source used for cooling and ready for spectral rendering.
    pub fn emitter(&self, source_id: u64) -> Option<&SpectralDiffuseThermalSurface> {
        self.patch(source_id).map(|patch| &patch.emitter)
    }
    pub fn rates(&self, source_id: u64) -> Option<&DiffuseExchangeRates> {
        self.patch(source_id).map(|patch| &patch.rates)
    }
    fn patch(&self, source_id: u64) -> Option<&Patch> {
        self.model
            .sources
            .binary_search_by_key(&source_id, |source| source.source_id)
            .ok()
            .map(|index| &self.patches[index])
    }

    fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), NativeDiffuseExchangeError> {
        check()?;
        for (&id, binding) in &self.model.cells {
            check()?;
            binding.validate(thermal, id)?;
        }
        let selected: Vec<_> = self.model.cells.keys().copied().collect();
        let available = thermal.exposed_areas_m2(&selected, &mut *check)?;
        let mut occupied = BTreeMap::<CellId, CheckedSum>::new();
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        let mut conductances = BTreeMap::<CellId, CheckedSum>::new();
        let mut patches = Vec::with_capacity(self.model.sources.len());
        let mut retained_nodes = 0usize;
        for (index, source) in self.model.sources.iter().enumerate() {
            check()?;
            self.work.source_evaluations = self
                .work
                .source_evaluations
                .checked_add(1)
                .ok_or(NativeDiffuseExchangeError::Budget)?;
            if self.work.source_evaluations > self.model.budget.max_source_evaluations {
                return Err(NativeDiffuseExchangeError::Budget);
            }
            let receiver = source.receiver;
            let (lo, hi) = receiver.temperature_bounds_k;
            if !positive(lo)
                || !positive(hi)
                || lo > hi
                || self.model.cells[&receiver.thermal_cell].material != receiver.material_id
            {
                return Err(invalid(
                    "receiver requires its native material and positive ordered kelvin support",
                ));
            }
            let emitter =
                source
                    .law
                    .resolve_surface(thermal, receiver.thermal_cell, receiver.material_id)?;
            let t = emitter.temperature_k();
            if t < lo || t > hi {
                return Err(SurfaceHeatingError::TemperatureDomain {
                    thermal_cell: receiver.thermal_cell,
                    temperature_k: t,
                    minimum_k: lo,
                    maximum_k: hi,
                }
                .into());
            }
            occupied
                .entry(receiver.thermal_cell)
                .or_default()
                .add(source.area_m2)?;
            let response = if let Some(old) = self.patches.get(index) {
                old.response.with_emissivity(emitter.emissivity())?
            } else {
                self.work.response_preparations += 1;
                exchange(check, |checkpoint| {
                    DiffuseExchangeResponse::from_emissivity(
                        emitter.emissivity(),
                        source.area_m2,
                        source.quadrature,
                        DiffuseExchangeBudget {
                            max_spectral_nodes: self
                                .model
                                .budget
                                .max_total_spectral_nodes
                                .saturating_sub(retained_nodes),
                        },
                        checkpoint,
                    )
                })?
            };
            retained_nodes = retained_nodes
                .checked_add(response.spectral_node_count())
                .ok_or(NativeDiffuseExchangeError::Budget)?;
            self.work.spectral_evaluations = self
                .work
                .spectral_evaluations
                .checked_add(response.spectral_node_count())
                .ok_or(NativeDiffuseExchangeError::Budget)?;
            if retained_nodes > self.model.budget.max_total_spectral_nodes
                || self.work.spectral_evaluations > self.model.budget.max_spectral_evaluations
            {
                return Err(NativeDiffuseExchangeError::Budget);
            }
            let bath = match source.bath {
                DiffuseExchangeBath::Vacuum => None,
                DiffuseExchangeBath::Blackbody { temperature_k } => Some(temperature_k),
                DiffuseExchangeBath::NativeAmbient => {
                    Some(thermal.cells()[&receiver.thermal_cell].boundary.ambient_k)
                }
            };
            let rates = exchange(check, |checkpoint| response.evaluate(t, bath, checkpoint))?;
            if let Some(maximum_w) = source.maximum_omitted_exchange_w {
                if rates.omitted_exchange_upper_bound_w > maximum_w {
                    return Err(NativeDiffuseExchangeError::Coverage {
                        source_id: source.source_id,
                        upper_bound_w: rates.omitted_exchange_upper_bound_w,
                        maximum_w,
                    });
                }
            }
            sources
                .entry(receiver.thermal_cell)
                .or_default()
                .add(-rates.outward_power_w)?;
            conductances
                .entry(receiver.thermal_cell)
                .or_default()
                .add(rates.secant_w_k)?;
            patches.push(Patch {
                emitter,
                response,
                rates,
            });
        }
        for (id, area) in occupied {
            if area.value()? > available[&id] {
                return Err(invalid(
                    "diffuse patches exceed the current native exposed area",
                ));
            }
        }
        self.sources_w = sources
            .into_iter()
            .map(|(id, sum)| Ok((id, sum.value()?)))
            .collect::<Result<_, OpticalHeatingError>>()?;
        self.conductances_w_k = conductances
            .into_iter()
            .map(|(id, sum)| Ok((id, sum.value()?)))
            .collect::<Result<_, OpticalHeatingError>>()?;
        check()?;
        self.patches = patches;
        self.prepared_at_s = thermal.time_s();
        Ok(())
    }
}

impl DepositionObserver for NativeDiffuseExchange {
    type Error = NativeDiffuseExchangeError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.work = NativeDiffuseExchangeWork::default();
        Ok(())
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        if thermal.time_s() != self.time_s {
            return Err(invalid("diffuse exchange clock mismatch"));
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
        if thermal.time_s() != self.time_s || thermal.time_s() != self.prepared_at_s {
            return Err(invalid(
                "diffuse timestep requires the prepared native start state",
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
            return Err(invalid("diffuse event clock mismatch"));
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
                "diffuse exchange requires its matching accepted native interval",
            ));
        }
        let mut sums = self.sums.clone();
        sums.duration.add(duration_s)?;
        for (sum, patch) in sums.patches.iter_mut().zip(&self.patches) {
            let heat = finite(patch.rates.outward_power_w * duration_s)?;
            sum.add(heat)?;
            sums.outward.add(heat)?;
            sums.omitted.add(finite(
                patch.rates.omitted_exchange_upper_bound_w * duration_s,
            )?)?;
        }
        for &watts in self.sources_w.values() {
            sums.external.add(finite(watts * duration_s)?)?;
        }
        // Charge exactly the rates actually applied BEFORE resolving the new T/h.
        self.refresh(thermal, check)?;
        let outward_j = sums.outward.value()?;
        let thermal_external_j = sums.external.value()?;
        let energy = NativeDiffuseExchangeEnergy {
            duration_s: sums.duration.value()?,
            outward_j,
            thermal_external_j,
            balance_residual_j: finite(outward_j + thermal_external_j)?,
            omitted_exchange_upper_bound_j: sums.omitted.value()?,
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

fn exchange<T>(
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(
        &mut dyn FnMut() -> Result<(), DiffuseExchangeError>,
    ) -> Result<T, DiffuseExchangeError>,
) -> Result<T, NativeDiffuseExchangeError> {
    let mut cancellation = None;
    let result = operation(&mut || match check() {
        Ok(()) => Ok(()),
        Err(error) => {
            cancellation = Some(error);
            Err(DiffuseExchangeError::Cancelled)
        }
    });
    match cancellation {
        Some(error) => Err(error.into()),
        None => result.map_err(Into::into),
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn invalid(reason: &'static str) -> NativeDiffuseExchangeError {
    NativeDiffuseExchangeError::Invalid(reason)
}
fn finite(value: f64) -> Result<f64, NativeDiffuseExchangeError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
