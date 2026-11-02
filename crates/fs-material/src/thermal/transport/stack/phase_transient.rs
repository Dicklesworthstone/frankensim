//! Distributed phase-changing layers on the original enthalpy timeline.
//!
//! Each cell retains its supplied mass and exact caloric h coordinate. The
//! original source-compiled k(T,f) supplies both sensible and latent tangents;
//! no apparent heat capacity or phase mixing rule is introduced here.
//! The native conduction observer and implicit/adaptive chain own all fluxes,
//! Newton steps, events, entropy and energy acceptance. Cellwise k and fixed
//! reference geometry are finite-volume approximations requiring refinement.

use super::transient::{LayeredThermalIds, MAX_LAYERED_THERMAL_CELLS};
use super::{LayerStackGeometry, MAX_CONDUCTION_LAYERS};
use crate::thermal::conduction::{LayerGeometry, phase::PreparedPhaseConductivity};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, DepositionSimulation, Event,
    MaterialId, Reservoir, ThermalLaw, ThermalLink, TimedEvent,
    conduction::{
        ConductionObserver, ConductiveLink, ConductiveReservoir, ConductivityLaw,
        implicit::{ImplicitConductionChain, ImplicitConductionError, ImplicitConductionLaw},
    },
    series_contact_conductance_w_k,
};
use std::collections::BTreeSet;

/// One homogeneous source-bound layer. Mass is the independently declared
/// whole-layer reference inventory; equilibrium phase density never changes it.
/// The source retains its original caloric object, claim and complete h domain.
#[derive(Clone, Copy, Debug)]
pub struct PhaseLayerSource<'a> {
    pub id: u64,
    pub thickness_m: f64,
    pub mass_kg: f64,
    pub cell_count: usize,
    pub conductivity: &'a PreparedPhaseConductivity<'a>,
}

#[derive(Clone, Debug)]
struct Cell {
    layer: usize,
    mass: f64,
    centroid: f64,
    geometry: CellGeometry,
}

/// Equal-volume phase control volumes with original source laws shared by layer.
/// The two physical faces are represented; planar lateral surfaces and cylinder
/// endcaps are insulated. Pressure/composition and geometry remain fixed. This
/// is reference-mass transport, not expansion work, advection or phase kinetics.
#[derive(Clone, Debug)]
pub struct PreparedPhaseLayeredThermalBody<'a> {
    layers: Vec<PhaseLayerSource<'a>>,
    geometry: LayerStackGeometry,
    cells: Vec<Cell>,
    links: Vec<LayerGeometry>,
    contacts: Vec<f64>,
    total_mass: f64,
    total_volume: f64,
    identity: ContentHash,
}

impl<'a> PreparedPhaseLayeredThermalBody<'a> {
    pub fn try_new(
        layers: &[PhaseLayerSource<'a>],
        contact_resistances_m2_k_w: &[f64],
        geometry: LayerStackGeometry,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if layers.is_empty()
            || layers.len() > MAX_CONDUCTION_LAYERS
            || contact_resistances_m2_k_w.len() != layers.len() - 1
            || contact_resistances_m2_k_w.iter().any(|&r| !nonnegative(r))
        {
            return Err(bad(
                "one bounded ordered phase stack and one nonnegative contact per interface required",
            ));
        }
        geometry.reference_area_m2().map_err(bad)?;
        let count = layers.iter().try_fold(0usize, |n, layer| {
            if layer.cell_count == 0 || !positive(layer.thickness_m) || !positive(layer.mass_kg) {
                return Err(bad(
                    "phase layer needs positive thickness, mass and cell count",
                ));
            }
            n.checked_add(layer.cell_count)
                .filter(|&n| n <= MAX_LAYERED_THERMAL_CELLS)
                .ok_or(DepositionError::Budget("phase layered thermal cells"))
        })?;
        let mut ids = BTreeSet::new();
        let mut cells = Vec::with_capacity(count);
        let mut position = geometry.start_m();
        let (mut total_mass, mut total_volume) = (0.0, 0.0);
        let mut hash =
            DomainHasher::new("org.frankensim.fs-material.phase-layered-thermal-body.v1");
        match geometry {
            LayerStackGeometry::Planar { area_m2 } => {
                hash.update(&[0]);
                hash.update(&area_m2.to_bits().to_le_bytes());
            }
            LayerStackGeometry::Cylindrical {
                inner_radius_m,
                length_m,
            } => {
                hash.update(&[1]);
                hash.update(&inner_radius_m.to_bits().to_le_bytes());
                hash.update(&length_m.to_bits().to_le_bytes());
            }
            LayerStackGeometry::Spherical { inner_radius_m } => {
                hash.update(&[2]);
                hash.update(&inner_radius_m.to_bits().to_le_bytes());
            }
        }
        for (index, layer) in layers.iter().enumerate() {
            checkpoint()?;
            if !ids.insert(layer.id) {
                return Err(bad("phase layer identifiers must be unique"));
            }
            let end = finite(position + layer.thickness_m)?;
            if end <= position {
                return Err(bad("phase layer width is below coordinate resolution"));
            }
            let bounds = [position, end];
            let volume = geometry.volume_m3(bounds).map_err(bad)?;
            let mass = layer.mass_kg / layer.cell_count as f64;
            let cell_volume = volume / layer.cell_count as f64;
            if !positive(mass) || !positive(cell_volume) {
                return Err(bad("phase cell inventory is not representable"));
            }
            let [lo, hi] = layer.conductivity.enthalpy_bounds_j_kg();
            if !lo.is_finite() || !hi.is_finite() || lo >= hi {
                return Err(DepositionError::MaterialDomain);
            }
            for local in 0..layer.cell_count {
                checkpoint()?;
                let next = geometry
                    .volume_fraction_position_m(
                        bounds,
                        (local + 1) as f64 / layer.cell_count as f64,
                    )
                    .map_err(bad)?;
                let centroid = geometry.centroid_m([position, next]).map_err(bad)?;
                cells.push(Cell {
                    layer: index,
                    mass,
                    centroid,
                    geometry: CellGeometry {
                        start_m: [position, 0.0, 0.0],
                        end_m: [next, 0.0, 0.0],
                        volume_m3: cell_volume,
                        surface_area_m2: geometry
                            .cell_surface_area_m2([position, next])
                            .map_err(bad)?,
                    },
                });
                position = next;
            }
            total_mass = finite(total_mass + layer.mass_kg)?;
            total_volume = finite(total_volume + volume)?;
            hash.update(&layer.id.to_le_bytes());
            hash.update(&layer.thickness_m.to_bits().to_le_bytes());
            hash.update(&layer.mass_kg.to_bits().to_le_bytes());
            hash.update(&(layer.cell_count as u64).to_le_bytes());
            hash.update(layer.conductivity.identity().as_bytes());
        }
        let mut links = Vec::with_capacity(count - 1);
        let mut contacts = Vec::with_capacity(count - 1);
        for pair in cells.windows(2) {
            checkpoint()?;
            let [a, b] = [&pair[0], &pair[1]];
            let face = a.geometry.end_m[0];
            let traces = [a.centroid, b.centroid];
            let link = match geometry {
                LayerStackGeometry::Planar { area_m2 } => {
                    LayerGeometry::planar([face - traces[0], traces[1] - face], area_m2)
                }
                LayerStackGeometry::Cylindrical { length_m, .. } => {
                    LayerGeometry::cylindrical(traces, face, length_m)
                }
                LayerStackGeometry::Spherical { .. } => LayerGeometry::spherical(traces, face),
            }?;
            links.push(link);
            contacts.push(if a.layer == b.layer {
                0.0
            } else {
                contact_resistances_m2_k_w[a.layer]
            });
        }
        for &contact in contact_resistances_m2_k_w {
            hash.update(&contact.to_bits().to_le_bytes());
        }
        checkpoint()?;
        Ok(Self {
            layers: layers.to_vec(),
            geometry,
            cells,
            links,
            contacts,
            total_mass,
            total_volume,
            identity: hash.finalize(),
        })
    }

    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn geometry(&self) -> LayerStackGeometry {
        self.geometry
    }
    pub fn layers(&self) -> &[PhaseLayerSource<'a>] {
        &self.layers
    }
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
    pub fn cell_geometries(&self) -> impl ExactSizeIterator<Item = &CellGeometry> {
        self.cells.iter().map(|cell| &cell.geometry)
    }
    pub fn cell_temperature_positions_m(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.cells.iter().map(|cell| cell.centroid)
    }
    pub fn link_geometries(&self) -> &[LayerGeometry] {
        &self.links
    }
    pub const fn total_mass_kg(&self) -> f64 {
        self.total_mass
    }
    pub const fn total_volume_m3(&self) -> f64 {
        self.total_volume
    }

    /// Bind optional reservoirs to the physical inner/outer faces. None is an
    /// insulated face; Some(0) prescribes the FACE temperature, with its bulk
    /// half-cell resistance retained. Positive values are areal film/contact
    /// resistances. One original cell supports at most one reservoir.
    pub fn bind(
        &self,
        ids: LayeredThermalIds,
        end_surface_resistances_m2_k_w: [Option<f64>; 2],
    ) -> Result<BoundPhaseLayeredThermalBody<'_, 'a>, DepositionError> {
        let n = self.cell_count();
        ids.first_cell
            .checked_add((n - 1) as u64)
            .ok_or(bad("phase cell identifiers overflow"))?;
        ids.first_material
            .checked_add((self.layers.len() - 1) as u64)
            .ok_or(bad("phase material identifiers overflow"))?;
        if n > 1 {
            ids.first_link
                .checked_add((n - 2) as u64)
                .ok_or(bad("phase link identifiers overflow"))?;
        }
        if end_surface_resistances_m2_k_w
            .iter()
            .flatten()
            .any(|&r| !nonnegative(r))
            || n == 1 && end_surface_resistances_m2_k_w.iter().all(Option::is_some)
        {
            return Err(bad(
                "nonnegative face resistances and at most one reservoir per cell required",
            ));
        }
        let mut reservoirs = [None; 2];
        for side in 0..2 {
            if let Some(resistance) = end_surface_resistances_m2_k_w[side] {
                let index = if side == 0 { 0 } else { n - 1 };
                let cell = &self.cells[index];
                let face = if side == 0 {
                    cell.geometry.start_m[0]
                } else {
                    cell.geometry.end_m[0]
                };
                let area = self.geometry.face_area_m2(face).map_err(bad)?;
                // Same original metric leaf; these resistance coordinates
                // never replace physical width, mass or volume.
                let distance = finite(
                    self.geometry
                        .resistance_length_m([face.min(cell.centroid), face.max(cell.centroid)])
                        .map_err(bad)?
                        * finite(area / self.geometry.reference_area_m2().map_err(bad)?)?,
                )?;
                if !positive(distance) {
                    return Err(bad(
                        "phase boundary half-cell resistance is not representable",
                    ));
                }
                reservoirs[side] = Some(Face {
                    index,
                    area,
                    distance,
                    resistance,
                });
            }
        }
        Ok(BoundPhaseLayeredThermalBody {
            body: self,
            ids,
            reservoirs,
        })
    }
}

#[derive(Clone, Copy)]
struct Face {
    index: usize,
    area: f64,
    distance: f64,
    resistance: f64,
}

/// Exact source objects and physical reference geometry bound to timeline IDs.
pub struct BoundPhaseLayeredThermalBody<'b, 'a> {
    body: &'b PreparedPhaseLayeredThermalBody<'a>,
    ids: LayeredThermalIds,
    reservoirs: [Option<Face>; 2],
}
impl<'a> BoundPhaseLayeredThermalBody<'_, 'a> {
    pub fn prepared_body(&self) -> &PreparedPhaseLayeredThermalBody<'a> {
        self.body
    }
    pub fn materials(&self) -> Vec<(MaterialId, &'a dyn ThermalLaw)> {
        self.body
            .layers
            .iter()
            .enumerate()
            .map(|(i, layer)| {
                (
                    MaterialId(self.ids.first_material + i as u64),
                    layer.conductivity.caloric_curve() as &dyn ThermalLaw,
                )
            })
            .collect()
    }
    pub fn cell_ids(&self) -> Vec<CellId> {
        (0..self.body.cell_count())
            .map(|i| self.cell_id(i))
            .collect()
    }

    /// Birth explicit specific enthalpies, so temperature on a latent plateau
    /// never invents its phase fraction. Reservoir temperatures are required
    /// exactly for the faces selected by bind. Cells precede their contacts.
    pub fn birth_events(
        &self,
        time_s: f64,
        specific_enthalpies_j_kg: &[f64],
        reservoir_temperatures_k: [Option<f64>; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<TimedEvent>, DepositionError> {
        checkpoint()?;
        if !nonnegative(time_s)
            || specific_enthalpies_j_kg.len() != self.body.cell_count()
            || (0..2).any(|side| {
                self.reservoirs[side].is_some() != reservoir_temperatures_k[side].is_some()
            })
            || reservoir_temperatures_k
                .iter()
                .flatten()
                .any(|&t| !positive(t))
        {
            return Err(bad(
                "phase birth needs one h per cell and exactly the declared physical face reservoirs",
            ));
        }
        let mut states = Vec::with_capacity(self.body.cell_count());
        for (cell, &h) in self.body.cells.iter().zip(specific_enthalpies_j_kg) {
            checkpoint()?;
            states.push(
                self.body.layers[cell.layer]
                    .conductivity
                    .state_at_specific_enthalpy(h)
                    .map_err(|_| DepositionError::MaterialDomain)?,
            );
        }
        let mut events = Vec::with_capacity(2 * self.body.cell_count() - 1);
        for (i, cell) in self.body.cells.iter().enumerate() {
            checkpoint()?;
            let mut reservoir = None;
            for (side, face) in self.reservoirs.iter().enumerate() {
                if let Some(face) = face.filter(|face| face.index == i) {
                    let k = states[i].conductivity_w_m_k();
                    let conductance = if k == 0.0 {
                        0.0
                    } else {
                        let r = finite(face.distance / k + face.resistance)?;
                        let g = finite(face.area / r)?;
                        if !positive(r) || !positive(g) {
                            return Err(DepositionError::Arithmetic);
                        }
                        g
                    };
                    reservoir = Some(Reservoir {
                        area_m2: face.area,
                        conductance_w_k: conductance,
                        temperature_k: reservoir_temperatures_k[side]
                            .expect("face temperature admitted"),
                    });
                }
            }
            events.push(TimedEvent {
                time_s,
                event: Event::Birth(DepositedCell {
                    id: self.cell_id(i),
                    segment_id: self.body.layers[cell.layer].id,
                    body_id: self.ids.body_id,
                    material_id: MaterialId(self.ids.first_material + cell.layer as u64),
                    mass_kg: cell.mass,
                    specific_h_j_kg: specific_enthalpies_j_kg[i],
                    geometry: cell.geometry,
                    boundary: Boundary {
                        ambient_k: states[i].temperature_k(),
                        htc_w_m2_k: 0.0,
                        reservoir,
                        source_w: 0.0,
                    },
                }),
            });
        }
        for (i, geometry) in self.body.links.iter().enumerate() {
            checkpoint()?;
            events.push(TimedEvent {
                time_s,
                event: Event::Contact(ThermalLink {
                    id: self.ids.first_link + i as u64,
                    a: self.cell_id(i),
                    b: self.cell_id(i + 1),
                    area_m2: geometry.contact_area_m2(),
                    conductance_w_k: series_contact_conductance_w_k(
                        geometry.contact_area_m2(),
                        geometry.equivalent_distances_m(),
                        [
                            states[i].conductivity_w_m_k(),
                            states[i + 1].conductivity_w_m_k(),
                        ],
                        self.body.contacts[i],
                    )?,
                }),
            });
        }
        checkpoint()?;
        Ok(events)
    }

    /// Prepare the original explicit observer; it also supplies the geometric
    /// coefficients to the native implicit chain. Other observer composition
    /// retains the original rule that each physical coefficient has one owner.
    pub fn conduction_observer(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ConductionObserver<'a>, DepositionError> {
        self.admit_timeline(thermal, &mut checkpoint)?;
        let laws: Vec<_> = self
            .body
            .layers
            .iter()
            .enumerate()
            .map(|(i, layer)| {
                (
                    MaterialId(self.ids.first_material + i as u64),
                    layer.conductivity as &'a dyn ConductivityLaw,
                )
            })
            .collect();
        let links: Vec<_> = self
            .body
            .links
            .iter()
            .enumerate()
            .map(|(i, geometry)| ConductiveLink {
                link_id: self.ids.first_link + i as u64,
                half_distances_m: geometry.equivalent_distances_m(),
                contact_resistance_m2_k_w: self.body.contacts[i],
            })
            .collect();
        let reservoirs: Vec<_> = self
            .reservoirs
            .iter()
            .flatten()
            .map(|face| ConductiveReservoir {
                cell_id: self.cell_id(face.index),
                half_distance_m: face.distance,
                surface_resistance_m2_k_w: face.resistance,
            })
            .collect();
        ConductionObserver::new(thermal, &laws, &links, &reservoirs, checkpoint)
    }

    /// Use the existing full enthalpy residual, zero latent dT/dh, analytic
    /// source dk/dh, tridiagonal Newton, exact-time events and adaptive owner.
    /// This direct entry requires the complete active timeline to be this body.
    /// Its returned chain can attach original two-sided contact laws through
    /// with_contacts; selected authored resistances must then be explicit zero.
    pub fn implicit_chain(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitConductionChain<'a>, ImplicitConductionError> {
        if thermal.cells().len() != self.body.cell_count()
            || thermal.links().len() != self.body.links.len()
        {
            return Err(bad("phase chain requires this complete ordered thermal body").into());
        }
        let observer = self.conduction_observer(thermal, &mut checkpoint)?;
        let laws: Vec<_> = self
            .body
            .layers
            .iter()
            .enumerate()
            .map(|(i, layer)| {
                (
                    MaterialId(self.ids.first_material + i as u64),
                    layer.conductivity as &'a dyn ImplicitConductionLaw,
                )
            })
            .collect();
        ImplicitConductionChain::new(thermal, &observer, &laws, checkpoint)
    }

    fn admit_timeline(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        for (i, expected) in self.body.cells.iter().enumerate() {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&self.cell_id(i))
                .ok_or(bad("phase cell is not active"))?;
            let layer = self.body.layers[expected.layer];
            let material = MaterialId(self.ids.first_material + expected.layer as u64);
            if cell.segment_id != layer.id
                || cell.body_id != self.ids.body_id
                || cell.material_id != material
                || cell.mass_kg != expected.mass
                || cell.geometry != expected.geometry
                || !thermal.uses_material_law(material, layer.conductivity.caloric_curve())
            {
                return Err(bad(
                    "phase cell differs from its original source, mass or geometry",
                ));
            }
            layer
                .conductivity
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(|_| DepositionError::MaterialDomain)?;
            let face = self
                .reservoirs
                .iter()
                .flatten()
                .find(|face| face.index == i);
            match (face, cell.boundary.reservoir) {
                (None, None) => {}
                (Some(face), Some(reservoir)) if face.area == reservoir.area_m2 => {}
                _ => {
                    return Err(bad(
                        "phase reservoir differs from its original physical face",
                    ));
                }
            }
        }
        for (i, expected) in self.body.links.iter().enumerate() {
            checkpoint()?;
            let link = thermal
                .links()
                .get(&(self.ids.first_link + i as u64))
                .ok_or(bad("phase layer contact is not active"))?;
            if link.a != self.cell_id(i)
                || link.b != self.cell_id(i + 1)
                || link.area_m2 != expected.contact_area_m2()
            {
                return Err(bad(
                    "phase contact differs from its original ordering or area",
                ));
            }
        }
        Ok(())
    }
    fn cell_id(&self, i: usize) -> CellId {
        CellId(self.ids.first_cell + i as u64)
    }
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn nonnegative(v: f64) -> bool {
    v.is_finite() && v >= 0.0
}
fn finite(v: f64) -> Result<f64, DepositionError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
fn bad(message: &'static str) -> DepositionError {
    DepositionError::InvalidInput(message)
}
