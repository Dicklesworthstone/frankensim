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
//! layers may instead grow or retreat at fixed footprint: their current mass
//! and either explicit constant density or the original caloric density law
//! determine every normal centre-to-face distance. This changes
//! thermal resistance, not the stored reference mesh or contact topology. The
//! coefficient-freezing time error and cell-resolution error require independent
//! refinement.

use super::*;

/// Structured implicit transport using the original backward-Euler residual.
#[path = "conduction/implicit.rs"]
pub mod implicit;

/// Explicit current face metrics while reference inventory stays immutable.
#[path = "conduction/geometry.rs"]
pub mod geometry;
pub use geometry::{CurrentConductionGeometry, CurrentConductiveBoundary, CurrentConductiveLink};

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

/// A fixed-footprint layer whose current volume is `mass / rho(h)`.
///
/// The registered thermal law must supply its original bulk density at the
/// same enthalpy as heat storage. Phase fractions and their mixture-volume law
/// remain with that material owner. This explicitly selected one-dimensional
/// geometry permits normal expansion/contraction; it supplies no lateral
/// strain, stress, moving mesh, gas displacement or new pressure law. The caller
/// must admit that shape approximation and the source's pressure constraint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarDensityLayer {
    pub cell_id: CellId,
    /// Fixed area of EACH normal face [m²].
    pub footprint_area_m2: f64,
}

#[derive(Clone, Copy)]
struct BoundPhaseLayer {
    layer: PlanarPhaseLayer,
    uses_material_density: bool,
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

/// A contact whose heat another owner carries: it keeps its covered area in
/// the thermal owner and has no conduction term here.
#[derive(Clone, Copy)]
struct OutsideLink {
    a: CellId,
    b: CellId,
    area_m2: f64,
}

impl OutsideLink {
    /// The contact is still the one that was named.
    fn admit(&self, thermal: &DepositionSimulation<'_>, id: u64) -> Result<(), DepositionError> {
        let link = thermal.links.get(&id).ok_or(DepositionError::InvalidInput(
            "link carried outside disappeared",
        ))?;
        if link.a != self.a || link.b != self.b || link.area_m2 != self.area_m2 {
            return Err(DepositionError::InvalidInput(
                "link carried outside changed",
            ));
        }
        Ok(())
    }
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
    /// Contacts named by [`Self::with_links_carried_outside`].
    outside: BTreeMap<u64, OutsideLink>,
    reservoirs: BTreeMap<CellId, BoundReservoir>,
    ambient_faces: BTreeMap<CellId, BoundAmbient>,
    conductivity: BTreeMap<CellId, f64>,
    reservoir_values: BTreeMap<CellId, f64>,
    ambient_values: BTreeMap<CellId, f64>,
    has_phase_layers: bool,
    has_material_density_layers: bool,
    // Complete current transport metrics; original face bindings stay intact.
    current_geometry: Option<CurrentConductionGeometry>,
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
            outside: BTreeMap::new(),
            reservoirs: BTreeMap::new(),
            ambient_faces: BTreeMap::new(),
            conductivity: BTreeMap::new(),
            reservoir_values: BTreeMap::new(),
            ambient_values: BTreeMap::new(),
            has_phase_layers: false,
            has_material_density_layers: false,
            current_geometry: None,
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
        if self.current_geometry.is_some() && !layers.is_empty() {
            return Err(DepositionError::InvalidInput(
                "mechanical face metrics and planar phase geometry need one explicit combined owner",
            ));
        }
        if !self.outside.is_empty() && !layers.is_empty() {
            return Err(DepositionError::InvalidInput(
                "planar phase geometry is not admitted beside links carried outside",
            ));
        }
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
                uses_material_density: false,
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

    /// Refresh planar resistance from the original material density as it heats,
    /// cools or changes phase. Storage and CFL still use the actual cell mass;
    /// density never recreates mass from the frozen reference volume.
    ///
    /// At binding, `rho(h) * reference_volume` must agree with the cell mass.
    /// Every later starting state uses `d = mass / (2 * rho(h) * footprint)`
    /// on links, reservoir faces and ambient faces alike. Endpoint density is
    /// checked before an advance can publish, including at latent plateaus.
    /// The existing enthalpy owner accounts heat under the material's pressure
    /// constraint; this adapter adds no second pressure-work term.
    ///
    /// Fixed-density and material-density layers may coexist on separate cells.
    /// Unknown density, unsupported states, changed reference geometry or face
    /// coverage refuse. Implicit chains/graphs admit fixed-mass layers when
    /// their original tangent law also supplies `d(1/rho)/dh`; their endpoint
    /// geometry and Newton derivative then evolve together. Temporal and
    /// spatial refinement remain necessary for either integration method.
    pub fn with_material_density_layers(
        self,
        thermal: &DepositionSimulation<'_>,
        layers: &[PlanarDensityLayer],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if layers.len() > self.cells.len() {
            return Err(DepositionError::InvalidInput(
                "material density layers exceed existing conduction cells",
            ));
        }
        let mut reference = Vec::with_capacity(layers.len());
        for layer in layers {
            checkpoint()?;
            let cell = thermal.cells.get(&layer.cell_id).ok_or(
                DepositionError::InvalidInput("inactive material density layer"),
            )?;
            reference.push(PlanarPhaseLayer {
                cell_id: layer.cell_id,
                footprint_area_m2: layer.footprint_area_m2,
                density_kg_m3: material_density(thermal, cell)?,
            });
        }
        let mut result = self.with_planar_phase_layers(thermal, &reference, &mut checkpoint)?;
        for layer in layers {
            result.cells.get_mut(&layer.cell_id)
                .expect("admitted planar cell")
                .phase_layer.as_mut().expect("admitted planar layer")
                .uses_material_density = true;
        }
        result.has_material_density_layers |= !layers.is_empty();
        result.prepare(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(result)
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
        if self.current_geometry.is_some() && !faces.is_empty() {
            return Err(DepositionError::InvalidInput(
                "bind every ambient face before supplying complete current metrics",
            ));
        }
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

    /// Name contact links whose heat another owner carries.
    ///
    /// A named link stays in the thermal owner: its two cells are still in
    /// contact, and the area it covers is still taken from their exposed
    /// faces. This binding gives it no conduction. The explicit timeline is
    /// handed a conductance of exactly zero for it in place of its authored
    /// one, and an implicit chain or graph prepared on this binding has no
    /// heat term and no tangent entry for it. What crosses the contact is
    /// then what the other owner puts into the two cells as sources: a
    /// coupling that generates heat on the interface, or a composition graph
    /// that joins the cells through their ports.
    ///
    /// A link has one owner. Refused: a link this binding already gives a
    /// conductance, an inactive link, a link named twice. An implicit contact
    /// or link law on a named link is refused by its own attachment, which
    /// finds no bound link there; the direct chain refuses a named link as a
    /// neighbouring pair without one. A contact or separation event that
    /// changes a named link is refused like one on a bound link. Later calls
    /// add links.
    ///
    /// Planar phase layers and current mechanical metrics are refused on a
    /// binding with named links, and named links on a binding with either:
    /// their coverage rules have not been extended to a contact this binding
    /// does not conduct through.
    pub fn with_links_carried_outside(
        mut self,
        thermal: &DepositionSimulation<'_>,
        link_ids: &[u64],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if !link_ids.is_empty() && (self.has_phase_layers || self.current_geometry.is_some()) {
            return Err(DepositionError::InvalidInput(
                "links carried outside are not admitted beside planar layers or current metrics",
            ));
        }
        // A repeated link ends the loop, so it runs at most once more than
        // the owner has links.
        for (i, &id) in link_ids.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let link = thermal.links.get(&id).ok_or(DepositionError::InvalidInput(
                "inactive link carried outside",
            ))?;
            if self.links.contains_key(&id) {
                return Err(DepositionError::InvalidInput(
                    "a link carried outside is also bound to a conductance",
                ));
            }
            let named = OutsideLink {
                a: link.a,
                b: link.b,
                area_m2: link.area_m2,
            };
            if self.outside.insert(id, named).is_some() {
                return Err(DepositionError::InvalidInput(
                    "repeated link carried outside",
                ));
            }
        }
        self.prepare(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }

    // Write each link carried outside into a map of link conductances as
    // exactly zero, after checking it is still the contact that was named.
    // The owner then keeps the contact's covered area and moves no heat
    // across it. This is the one place the zero is written: the explicit
    // timeline's coefficient map and the implicit endpoint's both come here,
    // so the two kernels cannot disagree about a named link.
    fn carry_outside(
        &self,
        thermal: &DepositionSimulation<'_>,
        conductances_w_k: &mut BTreeMap<u64, f64>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        for (i, (&id, named)) in self.outside.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            named.admit(thermal, id)?;
            conductances_w_k.insert(id, 0.0);
        }
        Ok(())
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
                let mut layer = phase.layer;
                if phase.uses_material_density {
                    layer.density_kg_m3 = material_density(thermal, cell)?;
                }
                phase.half_distance_m = layer.half_distance(cell.mass_kg)?;
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
            if (self.current_geometry.is_some()
                || self.has_phase_layers
                    && (self.cells[&link.a].phase_layer.is_some()
                        || self.cells[&link.b].phase_layer.is_some()))
                && thermal.contact_areas(link) != binding.covered_area_m2
            {
                return Err(DepositionError::InvalidInput(
                    "bound planar conduction contact coverage changed",
                ));
            }
            let geometry = binding.geometry;
            let metric = self.link_metrics(binding);
            let g = series_contact_conductance_w_k(
                metric.area_m2,
                [
                    self.half_distance(binding.a, metric.half_distances_m[0]),
                    self.half_distance(binding.b, metric.half_distances_m[1]),
                ],
                [self.conductivity[&binding.a], self.conductivity[&binding.b]],
                geometry.contact_resistance_m2_k_w,
            )?;
            links.insert(id, g);
        }
        self.carry_outside(thermal, &mut links, checkpoint)?;
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
            let metric = self.reservoir_metrics(binding);
            let g = surface_conductance_w_k(
                metric.area_m2,
                self.half_distance(id, metric.half_distance_m),
                k,
                geometry.surface_resistance_m2_k_w,
            )?;
            self.reservoir_values.insert(id, g);
        }
        self.admit_ambient_areas(thermal, checkpoint)?;
        for (&id, binding) in &self.ambient_faces {
            checkpoint()?;
            let geometry = binding.geometry;
            let metric = self.ambient_metrics(binding);
            let g = surface_conductance_w_k(
                metric.area_m2,
                self.half_distance(id, metric.half_distance_m),
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

fn material_density(
    thermal: &DepositionSimulation<'_>,
    cell: &DepositedCell,
) -> Result<f64, DepositionError> {
    let law = thermal.materials[&cell.material_id];
    if law.storage_basis() != ThermalStorageBasis::Enthalpy {
        return Err(DepositionError::InvalidInput(
            "material density layers require enthalpy storage",
        ));
    }
    let density = law
        .bulk_density_kg_m3(cell.specific_h_j_kg)?
        .ok_or(DepositionError::InvalidInput(
            "material density layer requires the original caloric density law",
        ))?;
    if !positive(density) {
        return Err(DepositionError::MaterialDomain);
    }
    Ok(density)
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
                    || self.outside.get(&link.id).is_some_and(|b| {
                        b.a != link.a || b.b != link.b || b.area_m2 != link.area_m2
                    })
                {
                    return Err(DepositionError::InvalidInput(
                        "event changes bound conduction face",
                    ));
                }
            }
            Event::Separation(id)
                if self.links.contains_key(id) || self.outside.contains_key(id) =>
            {
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
        if self.current_geometry.is_some() {
            // Also check a final-time topology/coverage event, even when no
            // subsequent heat substep would otherwise revisit every binding.
            self.prepare(thermal, &mut || Ok(()))?;
            self.prepared = false;
        }
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.prepared = false;
        checkpoint()?;
        if self.has_material_density_layers {
            for (&id, binding) in &self.cells {
                if let Some(phase) = binding.phase_layer {
                    if phase.uses_material_density {
                        checkpoint()?;
                        let cell = thermal.cells.get(&id).ok_or(
                            DepositionError::InvalidInput("material density layer disappeared"),
                        )?;
                        if cell.material_id != binding.material_id
                            || cell.body_id != binding.body_id
                            || cell.segment_id != binding.segment_id
                            || cell.geometry != binding.geometry
                        {
                            return Err(DepositionError::InvalidInput(
                                "spatial conduction cell binding changed",
                            ));
                        }
                        let mut layer = phase.layer;
                        layer.density_kg_m3 = material_density(thermal, cell)?;
                        layer.half_distance(cell.mass_kg)?;
                    }
                }
            }
        }
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

