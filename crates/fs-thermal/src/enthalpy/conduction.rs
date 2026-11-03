//! State-dependent spatial conduction on the existing enthalpy timeline.
//!
//! Each finite-volume link has two declared centre-to-face distances and an
//! optional areal contact resistance. Its conductance is
//! `A / (d_a/k_a(h_a) + R_contact + d_b/k_b(h_b))`. Each resolved surface
//! reservoir has `G = A / (d/k(h) + R_surface)`: zero surface resistance places
//! the prescribed temperature at the FACE, not at the cell centre; `1/htc`
//! represents a finite film resistance. The reservoir supplies its actual
//! temperature and covered area. An explicitly selected ambient face uses the
//! same law with its original exposed area and its separate ambient bath.
//! Insulation is represented by zero conductivity or an unbound insulated
//! channel, never an infinite or silently clipped resistance.
//!
//! This observer reevaluates conductivity once per participating cell on every
//! actual explicit substep. The original timeline owns pairwise signed heat,
//! storage, CFL, event splitting, cancellation and acceptance. It also supplies
//! entropy when `advance_observed_with_entropy` is selected. No caloric law,
//! material mixing rule, phase geometry or second stepping kernel is inferred.
//! Reference geometry and mass remain frozen by default. Explicit planar phase
//! layers may instead grow or retreat at constant density and footprint: their
//! current mass determines every normal centre-to-face distance. This changes
//! thermal resistance, not the stored reference mesh or contact topology. The
//! coefficient-freezing time error and cell-resolution error require independent
//! refinement.

use super::*;

/// Structured implicit transport using the original backward-Euler residual.
#[path = "conduction/implicit.rs"]
pub mod implicit;

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

/// Centre-to-face material conduction on the original exposed ambient area.
/// The ambient bath remains the cell's authored ambient_k, including events.
/// This explicit binding owns total ambient G; authored htc is not an extra film.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductiveAmbient {
    pub cell_id: CellId,
    pub half_distance_m: f64,
    /// Zero: prescribed face temperature. Positive: finite film/contact.
    pub surface_resistance_m2_k_w: f64,
}

/// Explicit uniform planar layer whose thickness follows its conserved mass.
///
/// All bound transport channels on this cell must be normal to one of its two
/// fixed footprint faces. The initial reference volume, mass, density, total
/// surface area and authored half-distances must describe that same layer.
/// The caller owns this geometric/constitutive assertion; no orientation,
/// porosity, density evolution, sidewall area, nucleation or moving mesh is
/// inferred. The layer must retain strictly positive mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarPhaseLayer {
    pub cell_id: CellId,
    /// Fixed area of EACH normal face [m²], not their sum.
    pub footprint_area_m2: f64,
    /// Constant physical bulk density of this layer [kg/m³].
    pub density_kg_m3: f64,
}

#[derive(Clone, Copy)]
struct BoundPhaseLayer {
    layer: PlanarPhaseLayer,
    half_distance_m: f64,
    exposed_area_m2: f64,
    reservoir_area_m2: Option<f64>,
}

#[derive(Clone, Copy)]
struct BoundAmbient {
    geometry: ConductiveAmbient,
    area_m2: f64,
}

#[derive(Clone, Copy)]
struct BoundCell {
    material_id: MaterialId,
    body_id: u64,
    segment_id: u64,
    mass_kg: f64,
    geometry: CellGeometry,
    phase_layer: Option<BoundPhaseLayer>,
}

#[derive(Clone, Copy)]
struct BoundLink {
    geometry: ConductiveLink,
    a: CellId,
    b: CellId,
    area_m2: f64,
    covered_area_m2: [f64; 2],
}

#[derive(Clone, Copy)]
struct BoundReservoir {
    geometry: ConductiveReservoir,
    area_m2: f64,
}

/// Prepared geometry and source bindings for the original sparse thermal graph.
/// This may own a subset of links, reservoir and ambient faces and compose with other observers;
/// each overridden coefficient must have only one owner. It stores O(cells +
/// faces) data and performs O((cells + faces) log(cells)) work per substep.
#[derive(Clone)]
pub struct ConductionObserver<'a> {
    laws: BTreeMap<MaterialId, &'a dyn ConductivityLaw>,
    cells: BTreeMap<CellId, BoundCell>,
    links: BTreeMap<u64, BoundLink>,
    reservoirs: BTreeMap<CellId, BoundReservoir>,
    ambient_faces: BTreeMap<CellId, BoundAmbient>,
    conductivity: BTreeMap<CellId, f64>,
    reservoir_values: BTreeMap<CellId, f64>,
    ambient_values: BTreeMap<CellId, f64>,
    has_phase_layers: bool,
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
            ambient_faces: BTreeMap::new(),
            conductivity: BTreeMap::new(),
            reservoir_values: BTreeMap::new(),
            ambient_values: BTreeMap::new(),
            has_phase_layers: false,
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
                covered_area_m2: thermal.contact_areas(link),
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

    /// Let selected already-bound cells change mass at fixed planar footprint.
    ///
    /// Before each actual heat substep, the current accepted mass sets
    /// `d = (mass / density) / (2 * footprint)` in the ORIGINAL series-resistance
    /// laws. Both link sides, reservoir faces and ambient faces use the same
    /// refreshed distance and conductivity. The original heat owner continues
    /// to use that same mass for storage and CFL, so thicker frost both stores
    /// more energy and conducts less heat for a fixed conductivity.
    ///
    /// This is a caller-selected one-dimensional layer approximation. Reference
    /// geometry and face coverage stay fixed; no volume is silently reassigned
    /// to gas, and no mechanical displacement or pressure work is supplied.
    /// Unselected cells retain the strict fixed-mass behavior.
    pub fn with_planar_phase_layers(
        mut self,
        thermal: &DepositionSimulation<'_>,
        layers: &[PlanarPhaseLayer],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if layers.len() > self.cells.len() {
            return Err(DepositionError::InvalidInput(
                "planar phase bindings exceed existing conduction cells",
            ));
        }
        let ids = layers.iter().map(|layer| layer.cell_id).collect::<Vec<_>>();
        let exposed = thermal.exposed_areas_m2(&ids, &mut checkpoint)?;
        for &layer in layers {
            checkpoint()?;
            let binding = self
                .cells
                .get(&layer.cell_id)
                .ok_or(DepositionError::InvalidInput(
                    "planar phase layer needs an existing conduction cell",
                ))?;
            if binding.phase_layer.is_some()
                || !positive(layer.footprint_area_m2)
                || !positive(layer.density_kg_m3)
                || thermal.materials[&binding.material_id].storage_basis()
                    != ThermalStorageBasis::Enthalpy
            {
                return Err(DepositionError::InvalidInput(
                    "unique positive planar geometry and enthalpy storage required",
                ));
            }
            let half_distance_m = layer.half_distance(binding.mass_kg)?;
            let mass = finite(layer.density_kg_m3 * binding.geometry.volume_m3)?;
            let area = finite(2.0 * layer.footprint_area_m2)?;
            if !same_measure(mass, binding.mass_kg)
                || !same_measure(area, binding.geometry.surface_area_m2)
            {
                return Err(DepositionError::InvalidInput(
                    "planar phase reference mass, volume and two-face area disagree",
                ));
            }
            let exposed_area_m2 = exposed[&layer.cell_id];
            let reservoir_area_m2 = thermal.cells[&layer.cell_id]
                .boundary
                .reservoir
                .map(|r| r.area_m2);
            if reservoir_area_m2.is_some_and(|a| {
                a > layer.footprint_area_m2 && !same_measure(a, layer.footprint_area_m2)
            }) {
                return Err(DepositionError::InvalidInput(
                    "reservoir exceeds one planar footprint face",
                ));
            }
            self.cells
                .get_mut(&layer.cell_id)
                .expect("admitted phase cell")
                .phase_layer = Some(BoundPhaseLayer {
                layer,
                half_distance_m,
                exposed_area_m2,
                reservoir_area_m2,
            });
            self.has_phase_layers = true;
        }
        self.admit_planar_channels(thermal, &mut checkpoint)?;
        self.prepare(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }

    /// Bind selected complete exposed ambient faces to the same material owner.
    /// Coverage, reference geometry and material identities remain fixed. Bath
    /// temperatures may change through boundary events. Repeated cells refuse;
    /// calling again may add previously unbound cells. A zero conductivity is
    /// exact insulation; zero surface resistance is a prescribed physical face.
    pub fn with_ambient_faces(
        mut self,
        thermal: &DepositionSimulation<'_>,
        faces: &[ConductiveAmbient],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if faces.len() > thermal.cells.len().saturating_sub(self.ambient_faces.len()) {
            return Err(DepositionError::InvalidInput(
                "ambient conduction bindings exceed active thermal owner",
            ));
        }
        let ids = faces.iter().map(|face| face.cell_id).collect::<Vec<_>>();
        let areas = thermal.exposed_areas_m2(&ids, &mut checkpoint)?;
        for face in faces {
            checkpoint()?;
            let area_m2 = areas[&face.cell_id];
            if !positive(face.half_distance_m)
                || !nonnegative(face.surface_resistance_m2_k_w)
                || !positive(area_m2)
                || self.ambient_faces.contains_key(&face.cell_id)
            {
                return Err(DepositionError::InvalidInput(
                    "ambient conduction needs a unique exposed face and valid material/film geometry",
                ));
            }
            self.bind_cell(thermal, face.cell_id)?;
            if let Some(phase) = self.cells[&face.cell_id].phase_layer {
                if !same_measure(
                    face.half_distance_m,
                    phase
                        .layer
                        .half_distance(self.cells[&face.cell_id].mass_kg)?,
                ) || area_m2 > 2.0 * phase.layer.footprint_area_m2
                    && !same_measure(area_m2, 2.0 * phase.layer.footprint_area_m2)
                {
                    return Err(DepositionError::InvalidInput(
                        "ambient face disagrees with the planar layer reference geometry",
                    ));
                }
            }
            self.ambient_faces.insert(
                face.cell_id,
                BoundAmbient {
                    geometry: *face,
                    area_m2,
                },
            );
        }
        self.prepare(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }

    // One native coverage sweep, preserving per-side contact measures and the
    // original subtraction order, instead of a link scan for every face.
    fn admit_ambient_areas(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        if self.ambient_faces.is_empty() && !self.has_phase_layers {
            return Ok(());
        }
        let mut ids = self.ambient_faces.keys().copied().collect::<Vec<_>>();
        if self.has_phase_layers {
            ids.extend(
                self.cells
                    .iter()
                    .filter_map(|(&id, cell)| cell.phase_layer.map(|_| id)),
            );
            ids.sort_unstable();
            ids.dedup();
        }
        let areas = thermal.exposed_areas_m2(&ids, &mut *checkpoint)?;
        for (&id, face) in &self.ambient_faces {
            checkpoint()?;
            if areas[&id] != face.area_m2 {
                return Err(DepositionError::InvalidInput(
                    "bound ambient conduction coverage changed",
                ));
            }
        }
        if self.has_phase_layers {
            for (&id, cell) in &self.cells {
                if let Some(phase) = cell.phase_layer {
                    checkpoint()?;
                    if areas[&id] != phase.exposed_area_m2 {
                        return Err(DepositionError::InvalidInput(
                            "bound planar layer face coverage changed",
                        ));
                    }
                }
            }
        }
        Ok(())
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
            phase_layer: None,
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
        for (i, (&id, binding)) in self.cells.iter_mut().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let cell = thermal.cells.get(&id).ok_or(DepositionError::InvalidInput(
                "bound conduction cell disappeared",
            ))?;
            if cell.material_id != binding.material_id
                || cell.body_id != binding.body_id
                || cell.segment_id != binding.segment_id
                || binding.phase_layer.is_none() && cell.mass_kg != binding.mass_kg
                || cell.geometry != binding.geometry
            {
                return Err(DepositionError::InvalidInput(
                    "spatial conduction cell binding changed",
                ));
            }
            if let Some(phase) = &mut binding.phase_layer {
                if cell.boundary.reservoir.map(|r| r.area_m2) != phase.reservoir_area_m2 {
                    return Err(DepositionError::InvalidInput(
                        "bound planar layer reservoir coverage changed",
                    ));
                }
                phase.half_distance_m = phase.layer.half_distance(cell.mass_kg)?;
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
            if self.has_phase_layers
                && (self.cells[&link.a].phase_layer.is_some()
                    || self.cells[&link.b].phase_layer.is_some())
                && thermal.contact_areas(link) != binding.covered_area_m2
            {
                return Err(DepositionError::InvalidInput(
                    "bound planar conduction contact coverage changed",
                ));
            }
            let geometry = binding.geometry;
            let g = series_contact_conductance_w_k(
                binding.area_m2,
                [
                    self.half_distance(binding.a, geometry.half_distances_m[0]),
                    self.half_distance(binding.b, geometry.half_distances_m[1]),
                ],
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
            let g = surface_conductance_w_k(
                binding.area_m2,
                self.half_distance(id, geometry.half_distance_m),
                k,
                geometry.surface_resistance_m2_k_w,
            )?;
            self.reservoir_values.insert(id, g);
        }
        self.admit_ambient_areas(thermal, checkpoint)?;
        for (&id, binding) in &self.ambient_faces {
            checkpoint()?;
            let geometry = binding.geometry;
            let g = surface_conductance_w_k(
                binding.area_m2,
                self.half_distance(id, geometry.half_distance_m),
                self.conductivity[&id],
                geometry.surface_resistance_m2_k_w,
            )?;
            self.ambient_values.insert(id, g);
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

    /// Last starting-state half-thickness of an explicitly selected phase layer.
    /// The reference cell geometry is unchanged; this is the transport distance.
    pub fn last_planar_half_distance_m(&self, id: CellId) -> Option<f64> {
        self.cells
            .get(&id)?
            .phase_layer
            .map(|phase| phase.half_distance_m)
    }

    fn half_distance(&self, id: CellId, reference: f64) -> f64 {
        if self.has_phase_layers {
            self.last_planar_half_distance_m(id).unwrap_or(reference)
        } else {
            reference
        }
    }

    fn admit_planar_channels(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        let admit = |id, distance: f64, area: f64| {
            if let Some(binding) = self.cells.get(&id) {
                if let Some(phase) = binding.phase_layer {
                    let footprint = phase.layer.footprint_area_m2;
                    if !same_measure(distance, phase.layer.half_distance(binding.mass_kg)?)
                        || (area > footprint && !same_measure(area, footprint))
                    {
                        return Err(DepositionError::InvalidInput(
                            "conduction channel disagrees with planar layer half-thickness or footprint",
                        ));
                    }
                }
            }
            Ok(())
        };
        // One canonical pass over faces; never scan every link for each layer.
        for (&id, bound) in &self.links {
            checkpoint()?;
            let areas = thermal.contact_areas(&thermal.links[&id]);
            for (side, cell_id) in [bound.a, bound.b].into_iter().enumerate() {
                admit(
                    cell_id,
                    bound.geometry.half_distances_m[side],
                    bound.area_m2,
                )?;
                admit(cell_id, bound.geometry.half_distances_m[side], areas[side])?;
            }
        }
        for (&id, bound) in &self.reservoirs {
            checkpoint()?;
            admit(id, bound.geometry.half_distance_m, bound.area_m2)?;
        }
        for (&id, bound) in &self.ambient_faces {
            checkpoint()?;
            // One ambient channel aggregates all remaining portions of both
            // normal faces at the same bath temperature.
            admit(id, bound.geometry.half_distance_m, 0.5 * bound.area_m2)?;
        }
        Ok(())
    }
}

fn same_measure(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}

impl PlanarPhaseLayer {
    fn half_distance(self, mass: f64) -> Result<f64, DepositionError> {
        let volume = finite(mass / self.density_kg_m3)?;
        let distance = finite((0.5 * volume) / self.footprint_area_m2)?;
        if !positive(mass) || !positive(volume) || !positive(distance) {
            return Err(DepositionError::InvalidInput(
                "planar phase layer needs representable positive mass and thickness",
            ));
        }
        Ok(distance)
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
    fn ambient_conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        if !self.prepared {
            return Err(DepositionError::InvalidInput(
                "evaluate conduction coefficients before ambient faces",
            ));
        }
        Ok(self.ambient_values.clone())
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
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
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
        // Event hooks receive the already-updated candidate. Check even at the
        // final endpoint, when no later coefficient query would expose a change.
        self.admit_ambient_areas(thermal, &mut || Ok(()))?;
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

// Shared by the explicit observer and implicit endpoint. This is the same
// original finite half-cell/film resistance arithmetic for either bath channel.
fn surface_conductance_w_k(
    area: f64,
    distance: f64,
    conductivity: f64,
    surface_resistance: f64,
) -> Result<f64, DepositionError> {
    if conductivity == 0.0 || area == 0.0 {
        return Ok(0.0);
    }
    let bulk = finite(distance / conductivity)?;
    let resistance = finite(bulk + surface_resistance)?;
    let value = finite(area / resistance)?;
    if !positive(bulk) || !positive(resistance) || !positive(value) {
        return Err(DepositionError::Arithmetic);
    }
    Ok(value)
}

#[cfg(test)]
#[path = "conduction/tests.rs"]
mod tests;
