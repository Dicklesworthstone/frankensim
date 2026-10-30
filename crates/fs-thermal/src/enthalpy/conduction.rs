//! State-dependent spatial conduction on the existing enthalpy timeline.
//!
//! Each finite-volume link has two declared centre-to-face distances and an
//! optional areal contact resistance. Its conductance is
//! `A / (d_a/k_a(h_a) + R_contact + d_b/k_b(h_b))`. Each resolved surface
//! reservoir has `G = A / (d/k(h) + R_surface)`: zero surface resistance places
//! the prescribed temperature at the FACE, not at the cell centre; `1/htc`
//! represents a finite film resistance. The reservoir supplies its actual
//! temperature and covered area. Insulation is represented by zero conductivity
//! or no reservoir, never an infinite or silently clipped resistance.
//!
//! This observer reevaluates conductivity once per participating cell on every
//! actual explicit substep. The original timeline owns pairwise signed heat,
//! storage, CFL, event splitting, cancellation and acceptance. It also supplies
//! entropy when `advance_observed_with_entropy` is selected. No caloric law,
//! material mixing rule, phase geometry or second stepping kernel is inferred.
//! Reference geometry and mass remain frozen. The coefficient-freezing time
//! error and cell-resolution error require independent refinement.

use super::*;

/// Immutable, explicitly bound scalar conductivity on an enthalpy domain.
/// A physical material owner supplies pressure, phase fractions, source data and
/// any effective-medium closure. Zero is an explicit insulating state.
pub trait ConductivityLaw {
    fn conductivity_w_m_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError>;
}

/// Geometry normal to one original link. Area and endpoint IDs come from that
/// existing link, and remain bound for the lifetime of this prepared observer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductiveLink {
    pub link_id: u64,
    pub half_distances_m: [f64; 2],
    pub contact_resistance_m2_k_w: f64,
}

/// Centre-to-face conduction followed by a declared external film/contact.
/// The original cell must already have a reservoir whose area is this face.
/// Its temperature may subsequently change through an ordinary boundary event.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductiveReservoir {
    pub cell_id: CellId,
    pub half_distance_m: f64,
    /// Zero: prescribed face temperature. Positive: finite film/contact.
    pub surface_resistance_m2_k_w: f64,
}

#[derive(Clone, Copy)]
struct BoundCell {
    material_id: MaterialId,
    body_id: u64,
    segment_id: u64,
    mass_kg: f64,
    geometry: CellGeometry,
}

#[derive(Clone, Copy)]
struct BoundLink {
    geometry: ConductiveLink,
    a: CellId,
    b: CellId,
    area_m2: f64,
}

#[derive(Clone, Copy)]
struct BoundReservoir {
    geometry: ConductiveReservoir,
    area_m2: f64,
}

/// Prepared geometry and source bindings for the original sparse thermal graph.
/// This may own a subset of links/reservoirs and compose with other observers;
/// each overridden coefficient must have only one owner. It stores O(cells +
/// faces) data and performs O((cells + faces) log(cells)) work per substep.
#[derive(Clone)]
pub struct ConductionObserver<'a> {
    laws: BTreeMap<MaterialId, &'a dyn ConductivityLaw>,
    cells: BTreeMap<CellId, BoundCell>,
    links: BTreeMap<u64, BoundLink>,
    reservoirs: BTreeMap<CellId, BoundReservoir>,
    conductivity: BTreeMap<CellId, f64>,
    reservoir_values: BTreeMap<CellId, f64>,
    prepared: bool,
}

impl<'a> ConductionObserver<'a> {
    /// Bind already-active cells/faces, validating original geometry and the
    /// initial source domain. Cell/link budgets are inherited from the admitted
    /// timeline. No laws or reference densities are selected by material name.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        laws: &[(MaterialId, &'a dyn ConductivityLaw)],
        links: &[ConductiveLink],
        reservoirs: &[ConductiveReservoir],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if links.len() > thermal.links.len()
            || reservoirs.len() > thermal.cells.len()
            || laws.len() > thermal.materials.len()
        {
            return Err(DepositionError::InvalidInput(
                "spatial conduction bindings exceed active thermal owner",
            ));
        }
        let mut result = Self {
            laws: BTreeMap::new(),
            cells: BTreeMap::new(),
            links: BTreeMap::new(),
            reservoirs: BTreeMap::new(),
            conductivity: BTreeMap::new(),
            reservoir_values: BTreeMap::new(),
            prepared: false,
        };
        for (i, &(id, law)) in laws.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            if !thermal.materials.contains_key(&id) || result.laws.insert(id, law).is_some() {
                return Err(DepositionError::InvalidInput(
                    "unknown or repeated conductivity material",
                ));
            }
        }
        for (i, &geometry) in links.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            if !geometry.half_distances_m.into_iter().all(positive)
                || !nonnegative(geometry.contact_resistance_m2_k_w)
            {
                return Err(DepositionError::InvalidInput(
                    "positive half-cell distances and nonnegative contact resistance required",
                ));
            }
            let link =
                thermal
                    .links
                    .get(&geometry.link_id)
                    .ok_or(DepositionError::InvalidInput(
                        "inactive spatial conduction link",
                    ))?;
            let bound = BoundLink {
                geometry,
                a: link.a,
                b: link.b,
                area_m2: link.area_m2,
            };
            if result.links.insert(link.id, bound).is_some() {
                return Err(DepositionError::InvalidInput(
                    "repeated spatial conduction link",
                ));
            }
            result.bind_cell(thermal, link.a)?;
            result.bind_cell(thermal, link.b)?;
        }
        for (i, &geometry) in reservoirs.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            if !positive(geometry.half_distance_m)
                || !nonnegative(geometry.surface_resistance_m2_k_w)
            {
                return Err(DepositionError::InvalidInput(
                    "positive reservoir half-cell distance and nonnegative surface resistance required",
                ));
            }
            let reservoir = thermal
                .cells
                .get(&geometry.cell_id)
                .and_then(|c| c.boundary.reservoir)
                .ok_or(DepositionError::InvalidInput(
                    "inactive spatial conduction reservoir",
                ))?;
            let bound = BoundReservoir {
                geometry,
                area_m2: reservoir.area_m2,
            };
            if result.reservoirs.insert(geometry.cell_id, bound).is_some() {
                return Err(DepositionError::InvalidInput(
                    "repeated spatial conduction reservoir",
                ));
            }
            result.bind_cell(thermal, geometry.cell_id)?;
        }
        result.prepare(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(result)
    }

    fn bind_cell(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        id: CellId,
    ) -> Result<(), DepositionError> {
        let cell = thermal.cells.get(&id).ok_or(DepositionError::InvalidInput(
            "inactive spatial conduction cell",
        ))?;
        if !self.laws.contains_key(&cell.material_id) {
            return Err(DepositionError::InvalidInput(
                "missing explicitly assigned conductivity law",
            ));
        }
        self.cells.entry(id).or_insert(BoundCell {
            material_id: cell.material_id,
            body_id: cell.body_id,
            segment_id: cell.segment_id,
            mass_kg: cell.mass_kg,
            geometry: cell.geometry,
        });
        Ok(())
    }

    fn prepare(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, DepositionError> {
        self.prepared = false;
        checkpoint()?;
        for (i, (&id, binding)) in self.cells.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let cell = thermal.cells.get(&id).ok_or(DepositionError::InvalidInput(
                "bound conduction cell disappeared",
            ))?;
            if cell.material_id != binding.material_id
                || cell.body_id != binding.body_id
                || cell.segment_id != binding.segment_id
                || cell.mass_kg != binding.mass_kg
                || cell.geometry != binding.geometry
            {
                return Err(DepositionError::InvalidInput(
                    "spatial conduction cell binding changed",
                ));
            }
            let k = self.laws[&binding.material_id].conductivity_w_m_k(cell.specific_h_j_kg)?;
            if !nonnegative(k) {
                return Err(DepositionError::MaterialDomain);
            }
            self.conductivity.insert(id, k);
        }
        let mut links = BTreeMap::new();
        for (i, (&id, binding)) in self.links.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let link = thermal.links.get(&id).ok_or(DepositionError::InvalidInput(
                "bound conduction link disappeared",
            ))?;
            if link.a != binding.a || link.b != binding.b || link.area_m2 != binding.area_m2 {
                return Err(DepositionError::InvalidInput(
                    "spatial conduction link binding changed",
                ));
            }
            let geometry = binding.geometry;
            let g = series_contact_conductance_w_k(
                binding.area_m2,
                geometry.half_distances_m,
                [self.conductivity[&binding.a], self.conductivity[&binding.b]],
                geometry.contact_resistance_m2_k_w,
            )?;
            links.insert(id, g);
        }
        for (i, (&id, binding)) in self.reservoirs.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let reservoir =
                thermal.cells[&id]
                    .boundary
                    .reservoir
                    .ok_or(DepositionError::InvalidInput(
                        "bound conduction reservoir disappeared",
                    ))?;
            if reservoir.area_m2 != binding.area_m2 {
                return Err(DepositionError::InvalidInput(
                    "spatial conduction reservoir area changed",
                ));
            }
            let k = self.conductivity[&id];
            let geometry = binding.geometry;
            let g = if k == 0.0 || binding.area_m2 == 0.0 {
                0.0
            } else {
                let bulk = finite(geometry.half_distance_m / k)?;
                let resistance = finite(bulk + geometry.surface_resistance_m2_k_w)?;
                let value = finite(binding.area_m2 / resistance)?;
                if !positive(bulk) || !positive(resistance) || !positive(value) {
                    return Err(DepositionError::Arithmetic);
                }
                value
            };
            self.reservoir_values.insert(id, g);
        }
        checkpoint()?;
        self.prepared = true;
        Ok(links)
    }

    /// Conductivities used at the most recent starting-state evaluation.
    /// They are diagnostics of that substep, not final-state material queries.
    pub fn last_conductivities_w_m_k(&self) -> &BTreeMap<CellId, f64> {
        &self.conductivity
    }
}

impl DepositionObserver for ConductionObserver<'_> {
    type Error = DepositionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.prepared = false;
        Ok(())
    }
    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.prepare(thermal, checkpoint)
    }
    fn reservoir_conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        if !self.prepared {
            return Err(DepositionError::InvalidInput(
                "evaluate conduction coefficients before reservoirs",
            ));
        }
        Ok(self.reservoir_values.clone())
    }
    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.prepared = false;
        match event {
            Event::Contact(link) => {
                if self
                    .links
                    .get(&link.id)
                    .is_some_and(|b| b.a != link.a || b.b != link.b || b.area_m2 != link.area_m2)
                {
                    return Err(DepositionError::InvalidInput(
                        "event changes bound conduction face",
                    ));
                }
            }
            Event::Separation(id) if self.links.contains_key(id) => {
                return Err(DepositionError::InvalidInput(
                    "event removes bound conduction face",
                ));
            }
            Event::Boundary(id, boundary) => {
                if let Some(binding) = self.reservoirs.get(id) {
                    if boundary
                        .reservoir
                        .is_none_or(|r| r.area_m2 != binding.area_m2)
                    {
                        return Err(DepositionError::InvalidInput(
                            "event changes bound reservoir face",
                        ));
                    }
                }
            }
            _ => {}
        }
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.prepared = false;
        checkpoint()
    }
}

#[cfg(test)]
#[path = "conduction/tests.rs"]
mod tests;
