//! Source-backed planar control volumes with independent cell temperatures.
//!
//! Preparation subdivides each homogeneous layer uniformly. Storage comes from
//! its original capacity source; a face conducts through the two incident half
//! cells and an explicit interlayer contact. The original thermal timeline owns
//! all evolution and energy/entropy accounts. Recompute `step_coefficients` on
//! EVERY implicit candidate; initial event conductances are only initial values.
//!
//! Spatial and temporal resolution remain independent caller choices. Only the
//! two axial faces are represented; lateral surfaces are insulated. Fixed mass,
//! geometry, composition and caloric constraints are required. No phase change,
//! moisture evolution, expansion work or automatic coupon matching is inferred.

use super::{LayerStackSourceQuery, PreparedLayerStackConduction};
use crate::thermal::{
    conduction::{LayerConductionConfig, LayerConductionResponse, TwoLayerConduction},
    transport::{
        ThermalTransportError,
        storage::{
            HeatStorageAmount, HeatStorageComponentQuery, HeatStorageThermalLaw,
            PreparedHeatStorage,
        },
    },
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, DepositionSimulation, Event,
    MaterialId, ThermalLaw, ThermalLink, TimedEvent, backward_euler::StepCoefficients,
};
use std::collections::BTreeMap;

/// Bound retained geometry and per-candidate edge work before allocation.
pub const MAX_LAYERED_THERMAL_CELLS: usize = 512;

/// Whole-layer physical inputs; capacity amount is divided among equal cells.
///
/// A mass capacity must name the supplied whole-layer mass. A volume capacity
/// must name area times whole-layer thickness. Mass remains independently
/// supplied for volumetric storage; it is never manufactured from nominal data.
/// The caller admits compatibility of the caloric and conductivity populations
/// and work constraints. A common card and context do not prove a common coupon.
#[derive(Clone, Copy, Debug)]
pub struct TransientLayerSourceQuery<'a> {
    pub conduction: LayerStackSourceQuery<'a>,
    pub capacity: HeatStorageComponentQuery<'a>,
    pub mass_kg: f64,
    pub cell_count: usize,
}

#[derive(Clone, Debug)]
struct Cell {
    layer: usize,
    width: f64,
    geometry: CellGeometry,
}

/// One source preparation per layer, shared by all its spatial control volumes.
/// Preparation allocates O(layers + cells + retained source knots); capacity
/// curves are not duplicated for each subdivision. The source support is the
/// common interval of every capacity and conductivity path.
#[derive(Clone, Debug)]
pub struct PreparedLayeredThermalBody {
    conduction: PreparedLayerStackConduction,
    storages: Vec<PreparedHeatStorage>,
    masses: Vec<f64>,
    counts: Vec<usize>,
    cells: Vec<Cell>,
    datum: f64,
    span: [f64; 2],
    total_mass: f64,
    total_volume: f64,
    identity: ContentHash,
}

impl PreparedLayeredThermalBody {
    pub fn try_new(
        queries: &[TransientLayerSourceQuery<'_>],
        contact_resistances_m2_k_w: &[f64],
        area_m2: f64,
        datum_temperature_k: f64,
    ) -> Result<Self, ThermalTransportError> {
        let count = queries.iter().try_fold(0usize, |n, q| {
            if q.cell_count == 0 {
                return Err(invalid("each thermal layer needs at least one cell"));
            }
            n.checked_add(q.cell_count)
                .filter(|&v| v <= MAX_LAYERED_THERMAL_CELLS)
                .ok_or_else(|| invalid("layered thermal body exceeds its 512-cell budget"))
        })?;
        let conduction = PreparedLayerStackConduction::try_new(
            &queries.iter().map(|q| q.conduction).collect::<Vec<_>>(),
            contact_resistances_m2_k_w,
            area_m2,
        )?;
        let mut storages = Vec::with_capacity(queries.len());
        let mut masses = Vec::with_capacity(queries.len());
        let mut counts = Vec::with_capacity(queries.len());
        let mut cells = Vec::with_capacity(count);
        let mut span = conduction.temperature_span_k();
        let mut total_mass = 0.0;
        let mut total_volume = 0.0;
        let mut mass_correction = 0.0;
        let mut volume_correction = 0.0;
        let mut x = 0.0;
        for (layer, q) in queries.iter().enumerate() {
            let volume = area_m2 * q.conduction.thickness_m;
            if !positive(q.mass_kg) || !positive(volume) {
                return Err(invalid("layer mass and volume must be positive and finite"));
            }
            if q.capacity.card.content_hash() != q.conduction.card.content_hash()
                || q.capacity.temperature_axis != q.conduction.temperature_axis
            {
                return Err(invalid(
                    "layer capacity and conductivity need the same original card and temperature axis",
                ));
            }
            for (axis, value) in q.capacity.lower.axes() {
                if axis != q.capacity.temperature_axis
                    && let Some(other) = q.conduction.lower.axes().get(axis)
                    && (value != other
                        || q.capacity.lower.axis_quantities().get(axis)
                            != q.conduction.lower.axis_quantities().get(axis))
                {
                    return Err(invalid(
                        "layer sources disagree on a shared fixed condition",
                    ));
                }
            }
            let amount = match q.capacity.amount {
                HeatStorageAmount::MassKg(mass) if same(mass, q.mass_kg) => {
                    HeatStorageAmount::MassKg(mass / q.cell_count as f64)
                }
                HeatStorageAmount::ReferenceVolumeM3(v) if same(v, volume) => {
                    HeatStorageAmount::ReferenceVolumeM3(v / q.cell_count as f64)
                }
                _ => {
                    return Err(invalid(
                        "capacity amount must equal whole-layer mass or geometric volume",
                    ));
                }
            };
            let storage = PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
                amount,
                ..q.capacity
            }])?;
            let supported = storage.temperature_span_k();
            span[0] = span[0].max(supported[0]);
            span[1] = span[1].min(supported[1]);
            let mass = q.mass_kg / q.cell_count as f64;
            let width = q.conduction.thickness_m / q.cell_count as f64;
            let cell_volume = volume / q.cell_count as f64;
            if !positive(mass)
                || !positive(width)
                || !positive(cell_volume)
                || !positive(2.0 * area_m2)
            {
                return Err(invalid("layer subdivision is not representable"));
            }
            // Validate the exact original storage/mass conversion now, including
            // finite enthalpy range and the thermal owner's slope requirement.
            HeatStorageThermalLaw::try_new(&storage, mass, datum_temperature_k)?;
            for _ in 0..q.cell_count {
                let end = x + width;
                if !end.is_finite() || end <= x {
                    return Err(invalid("cell position is not representable"));
                }
                cells.push(Cell {
                    layer,
                    width,
                    geometry: CellGeometry {
                        start_m: [x, 0.0, 0.0],
                        end_m: [end, 0.0, 0.0],
                        volume_m3: cell_volume,
                        // The planar model exposes axial faces only.
                        surface_area_m2: 2.0 * area_m2,
                    },
                });
                x = end;
            }
            accumulate(&mut total_mass, &mut mass_correction, q.mass_kg)?;
            accumulate(&mut total_volume, &mut volume_correction, volume)?;
            storages.push(storage);
            masses.push(mass);
            counts.push(q.cell_count);
        }
        if span[0] >= span[1]
            || !datum_temperature_k.is_finite()
            || datum_temperature_k < span[0]
            || datum_temperature_k > span[1]
        {
            return Err(invalid(
                "layered storage and conduction need a common interval containing the datum",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.layered-thermal-body.v1");
        hash.update(conduction.identity().as_bytes());
        hash.update(&datum_temperature_k.to_bits().to_le_bytes());
        for i in 0..storages.len() {
            hash.update(storages[i].identity().as_bytes());
            hash.update(&masses[i].to_bits().to_le_bytes());
            hash.update(&(counts[i] as u64).to_le_bytes());
        }
        Ok(Self {
            conduction,
            storages,
            masses,
            counts,
            cells,
            datum: datum_temperature_k,
            span,
            total_mass,
            total_volume,
            identity: hash.finalize(),
        })
    }

    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn conduction(&self) -> &PreparedLayerStackConduction {
        &self.conduction
    }
    /// One PER-CELL capacity assembly per layer; subdivisions share its law.
    pub fn storages(&self) -> &[PreparedHeatStorage] {
        &self.storages
    }
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
    pub fn layer_cell_counts(&self) -> &[usize] {
        &self.counts
    }
    pub const fn temperature_span_k(&self) -> [f64; 2] {
        self.span
    }
    pub const fn total_mass_kg(&self) -> f64 {
        self.total_mass
    }
    pub const fn total_volume_m3(&self) -> f64 {
        self.total_volume
    }

    /// Keep this vector alive while its original borrowed laws serve a timeline.
    pub fn thermal_laws(
        &self,
        sensible_entropy: bool,
    ) -> Result<Vec<HeatStorageThermalLaw<'_>>, ThermalTransportError> {
        self.storages
            .iter()
            .zip(&self.masses)
            .map(|(storage, &mass)| {
                let law = HeatStorageThermalLaw::try_new(storage, mass, self.datum)?;
                Ok(if sensible_entropy {
                    law.with_sensible_entropy()
                } else {
                    law
                })
            })
            .collect()
    }

    pub fn bind<'a>(
        &'a self,
        laws: &'a [HeatStorageThermalLaw<'a>],
        ids: LayeredThermalIds,
    ) -> Result<BoundLayeredThermalBody<'a>, DepositionError> {
        if laws.len() != self.storages.len()
            || laws.iter().enumerate().any(|(i, law)| {
                !std::ptr::eq(law.storage(), &self.storages[i])
                    || law.reference_mass_kg() != self.masses[i]
                    || law.datum_temperature_k() != self.datum
            })
        {
            return Err(bad(
                "binding needs the original per-layer storage laws, masses and datum",
            ));
        }
        ids.first_cell
            .checked_add((self.cell_count() - 1) as u64)
            .ok_or_else(|| bad("cell identifiers overflow"))?;
        ids.first_material
            .checked_add((laws.len() - 1) as u64)
            .ok_or_else(|| bad("material identifiers overflow"))?;
        if self.cell_count() > 1 {
            ids.first_link
                .checked_add((self.cell_count() - 2) as u64)
                .ok_or_else(|| bad("link identifiers overflow"))?;
        }
        let mut hash =
            DomainHasher::new("org.frankensim.fs-material.bound-layered-thermal-body.v1");
        hash.update(self.identity.as_bytes());
        for law in laws {
            hash.update(law.identity().as_bytes());
        }
        for id in [
            ids.first_cell,
            ids.first_material,
            ids.first_link,
            ids.body_id,
        ] {
            hash.update(&id.to_le_bytes());
        }
        Ok(BoundLayeredThermalBody {
            body: self,
            laws,
            ids,
            identity: hash.finalize(),
        })
    }
}

/// Contiguous checked IDs for embedding the body in an existing thermal network.
#[derive(Clone, Copy, Debug)]
pub struct LayeredThermalIds {
    pub first_cell: u64,
    pub first_material: u64,
    pub first_link: u64,
    pub body_id: u64,
}

/// Immutable binding to the exact original law objects used by the timeline.
pub struct BoundLayeredThermalBody<'a> {
    body: &'a PreparedLayeredThermalBody,
    laws: &'a [HeatStorageThermalLaw<'a>],
    ids: LayeredThermalIds,
    identity: ContentHash,
}

impl BoundLayeredThermalBody<'_> {
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    pub fn materials(&self) -> BTreeMap<MaterialId, &dyn ThermalLaw> {
        self.laws
            .iter()
            .enumerate()
            .map(|(i, law)| {
                (
                    MaterialId(self.ids.first_material + i as u64),
                    law as &dyn ThermalLaw,
                )
            })
            .collect()
    }
    pub fn cell_ids(&self) -> Vec<CellId> {
        (0..self.body.cell_count())
            .map(|i| self.cell_id(i))
            .collect()
    }

    /// Birth all cells before their contacts, initially insulated. The caller
    /// may append ordinary Boundary events for heat sources and end reservoirs.
    /// Interior axial faces are covered by the original owner's contact rules.
    pub fn birth_events(
        &self,
        time_s: f64,
        initial_temperatures_k: &[f64],
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<TimedEvent>, DepositionError> {
        checkpoint()?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(bad("birth time must be nonnegative and finite"));
        }
        self.admit_temperatures(initial_temperatures_k)?;
        let responses = self.responses(initial_temperatures_k, config, &mut checkpoint)?;
        let mut events = Vec::with_capacity(2 * self.body.cell_count() - 1);
        for (i, cell) in self.body.cells.iter().enumerate() {
            checkpoint()?;
            let law = &self.laws[cell.layer];
            let deposited = DepositedCell {
                id: self.cell_id(i),
                segment_id: self.body.conduction.layers[cell.layer].id,
                body_id: self.ids.body_id,
                material_id: MaterialId(self.ids.first_material + cell.layer as u64),
                mass_kg: self.body.masses[cell.layer],
                specific_h_j_kg: law
                    .specific_storage_j_kg(initial_temperatures_k[i])
                    .map_err(|_| DepositionError::MaterialDomain)?,
                geometry: cell.geometry,
                boundary: Boundary {
                    ambient_k: self.body.datum,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            };
            law.admit_cell(&deposited)?;
            events.push(TimedEvent {
                time_s,
                event: Event::Birth(deposited),
            });
        }
        for (i, response) in responses.iter().enumerate() {
            checkpoint()?;
            events.push(TimedEvent {
                time_s,
                event: Event::Contact(ThermalLink {
                    id: self.ids.first_link + i as u64,
                    a: self.cell_id(i),
                    b: self.cell_id(i + 1),
                    area_m2: self.body.conduction.area,
                    conductance_w_k: response.conductance_w_k,
                }),
            });
        }
        checkpoint()?;
        Ok(events)
    }

    /// Evaluate at the original owner's candidate state. No state or source
    /// account changes here. Replacements retain that owner's contact coverage;
    /// each edge still exchanges the same signed joules on its two endpoints.
    /// The caller owns nonlinear solution, acceptance and time-step refinement.
    /// Sparse result maps/edge diagnostics allocate O(cells); scalar edge kernels
    /// allocate nothing and all work has layer/cell/iteration bounds.
    pub fn step_coefficients(
        &self,
        simulation: &DepositionSimulation<'_>,
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<StepCoefficients<Vec<LayerConductionResponse>>, DepositionError> {
        checkpoint()?;
        let mut temperatures = Vec::with_capacity(self.body.cell_count());
        for (i, prepared) in self.body.cells.iter().enumerate() {
            checkpoint()?;
            let cell = simulation
                .cells()
                .get(&self.cell_id(i))
                .ok_or_else(|| bad("layered body cell is not active"))?;
            let material = MaterialId(self.ids.first_material + prepared.layer as u64);
            if cell.material_id != material
                || cell.segment_id != self.body.conduction.layers[prepared.layer].id
                || cell.body_id != self.ids.body_id
                || cell.geometry != prepared.geometry
                || !simulation.uses_material_law(material, &self.laws[prepared.layer])
            {
                return Err(bad(
                    "timeline cell no longer matches the prepared layer and original law",
                ));
            }
            self.laws[prepared.layer].admit_cell(cell)?;
            temperatures.push(simulation.temperature_k(cell.id)?);
        }
        self.admit_temperatures(&temperatures)?;
        for i in 0..self.body.cell_count() - 1 {
            checkpoint()?;
            let link = simulation
                .links()
                .get(&(self.ids.first_link + i as u64))
                .ok_or_else(|| bad("layered body contact is not active"))?;
            if link.a != self.cell_id(i)
                || link.b != self.cell_id(i + 1)
                || link.area_m2 != self.body.conduction.area
            {
                return Err(bad(
                    "timeline contact differs from the prepared physical ordering or area",
                ));
            }
        }
        let responses = self.responses(&temperatures, config, &mut checkpoint)?;
        let conductances_w_k = responses
            .iter()
            .enumerate()
            .map(|(i, r)| (self.ids.first_link + i as u64, r.conductance_w_k))
            .collect();
        checkpoint()?;
        Ok(StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k,
            payload: responses,
        })
    }

    fn cell_id(&self, i: usize) -> CellId {
        CellId(self.ids.first_cell + i as u64)
    }

    fn admit_temperatures(&self, temperatures: &[f64]) -> Result<(), DepositionError> {
        if temperatures.len() != self.body.cell_count()
            || temperatures
                .iter()
                .any(|t| !t.is_finite() || *t < self.body.span[0] || *t > self.body.span[1])
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }

    fn responses(
        &self,
        temperatures: &[f64],
        config: LayerConductionConfig,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<LayerConductionResponse>, DepositionError> {
        if config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(bad(
                "invalid layer conduction tolerance or iteration budget",
            ));
        }
        let mut responses = Vec::with_capacity(self.body.cell_count() - 1);
        for (i, cells) in self.body.cells.windows(2).enumerate() {
            checkpoint()?;
            let [a, b] = [&cells[0], &cells[1]];
            let contact = if a.layer == b.layer {
                0.0
            } else {
                self.body.conduction.contacts[a.layer]
            };
            let edge = TwoLayerConduction::new(
                [
                    &self.body.conduction.curves[a.layer],
                    &self.body.conduction.curves[b.layer],
                ],
                [0.5 * a.width, 0.5 * b.width],
                self.body.conduction.area,
                contact,
            )?;
            responses.push(edge.evaluate(
                [temperatures[i], temperatures[i + 1]],
                config,
                &mut *checkpoint,
            )?);
        }
        Ok(responses)
    }
}

fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn same(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}
fn accumulate(
    total: &mut f64,
    correction: &mut f64,
    value: f64,
) -> Result<(), ThermalTransportError> {
    let adjusted = value - *correction;
    let next = *total + adjusted;
    *correction = (next - *total) - adjusted;
    *total = next;
    if !positive(next) {
        return Err(invalid("layer totals are not representable"));
    }
    Ok(())
}
fn invalid(reason: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(reason)
}
fn bad(reason: &'static str) -> DepositionError {
    DepositionError::InvalidInput(reason)
}
