//! Finite thermal storage and explicit surface geometry for a straight fin.
//!
//! This module lowers one authored fin into the original thermal cells and
//! links. It owns no temperature field, constitutive law, heat source or time
//! integrator. Cell masses and caloric identities remain explicit; neither a
//! steady fin efficiency nor the base body's mass substitutes for fin storage.
//! The fixed section is isothermal transversely. Its axial mesh and that model
//! assumption require separate physical admission by the material/film owner.

use fs_convection::fin::StraightFinGeometry;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, Event, MaterialId, ThermalLink,
    TimedEvent, conduction::ConductiveLink,
};
use std::collections::BTreeSet;

/// Explicit physical inventory of one axial control volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientFinCell {
    pub id: CellId,
    pub material_id: MaterialId,
    pub mass_kg: f64,
    pub specific_h_j_kg: f64,
    /// Original authored source watts and ambient temperature are retained.
    /// Ambient h must be zero and a separate reservoir must be absent: all
    /// selected convection belongs to the returned finite-fluid links.
    pub boundary: Boundary,
}

/// One physical fin, with its root and axis in the caller's spatial frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinPlacement {
    pub root_m: [f64; 3],
    /// Explicit unit vector; this constructor never silently normalizes it.
    pub direction: [f64; 3],
    pub body_id: u64,
    pub segment_id: u64,
}

/// Existing base/fluid endpoints and explicitly identified new fin topology.
///
/// Axial edges must run strictly from zero to the exact geometry length. Each
/// inventory entry is one cell between consecutive edges, and may retain its
/// own material and caloric chart. Base and fluid cells are supplied elsewhere;
/// their geometry must accommodate the returned contact areas. No base bare
/// surface, fluid storage, density, conductivity, or convection is inferred.
pub struct TransientFinInput<'a> {
    pub geometry: StraightFinGeometry,
    pub axial_edges_m: &'a [f64],
    pub inventory: &'a [TransientFinCell],
    pub placement: FinPlacement,
    pub base_cell: CellId,
    pub fluid_cell: CellId,
    /// Root link is ordered [base, first fin cell].
    pub root_link: u64,
    /// Interior links in axial order, each ordered [upstream, downstream].
    pub axial_links: &'a [u64],
    /// One lateral link per cell, ordered [fin cell, fluid].
    pub side_links: &'a [u64],
    /// Independent tip link [last fin cell, fluid]; None is an insulated tip.
    pub tip_link: Option<u64>,
    /// Physical base temperature trace to root face distance, strictly positive.
    pub base_half_distance_m: f64,
    /// Explicit root contact resistance, in series with both material paths.
    pub root_contact_resistance_m2_k_w: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TransientFinBudget {
    pub maximum_cells: usize,
    pub maximum_links: usize,
}

/// Geometry and original native inventory for one physical fin.
///
/// Root and internal contacts each cover the actual cross section on both
/// sides; lateral and optional tip contacts cover only their own surfaces.
/// Native cell surface area counts the authored exposed lateral perimeter
/// plus both axial faces. Omitted insulated lateral faces create no heat path.
/// Every link starts with zero conductance and requires its explicit material
/// or film owner. Bare base cooling must be a separate, disjoint surface.
#[derive(Clone, Debug)]
pub struct TransientFinTopology {
    geometry: StraightFinGeometry,
    axial_edges_m: Vec<f64>,
    base_cell: CellId,
    fluid_cell: CellId,
    cells: Vec<DepositedCell>,
    links: Vec<ThermalLink>,
    conductive_links: Vec<ConductiveLink>,
    side_links: Vec<u64>,
    tip_link: Option<u64>,
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite_positive(value: f64) -> Result<f64, DepositionError> {
    if positive(value) {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
// True when a cell surface holds both axial faces and its lateral face in
// every order the native owner can subtract or sum them (it subtracts in
// link-ID order on admission and sums per cell in its step).
fn holds_full_cover(surface: f64, face: f64, side: f64) -> bool {
    [[face, face, side], [face, side, face], [side, face, face]]
        .iter()
        .all(|order| {
            order.iter().fold(surface, |free, &covered| free - covered) >= 0.0
                && surface - order.iter().fold(0.0, |sum, &covered| sum + covered) >= 0.0
        })
}

impl TransientFinTopology {
    /// Construct bounded native geometry without changing an accepted state.
    /// Duplicate IDs, aliasing to base/fluid, invalid inventory or an unresolved
    /// spatial segment refuse before returning any topology. Source/caloric
    /// admission belongs to the original timeline and prepared source adapter.
    pub fn new(
        input: TransientFinInput<'_>,
        budget: TransientFinBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        let n = input.inventory.len();
        let link_count = n
            .checked_mul(2)
            .and_then(|count| count.checked_add(usize::from(input.tip_link.is_some())))
            .ok_or(DepositionError::Arithmetic)?;
        if n > budget.maximum_cells || link_count > budget.maximum_links {
            return Err(DepositionError::Budget("transient fin topology"));
        }
        if n == 0
            || input.axial_edges_m.len() != n.checked_add(1).ok_or(DepositionError::Arithmetic)?
            || input.axial_links.len() != n - 1
            || input.side_links.len() != n
            || input.axial_edges_m[0] != 0.0
            || input.axial_edges_m[n] != input.geometry.length_m()
            || input.base_cell == input.fluid_cell
            || !positive(input.base_half_distance_m)
            || !nonnegative(input.root_contact_resistance_m2_k_w)
        {
            return Err(DepositionError::InvalidInput(
                "a fin needs complete explicit inventory, ordered axial geometry and distinct base/fluid endpoints",
            ));
        }
        let placement = input.placement;
        let square: f64 = placement.direction.iter().map(|x| x * x).sum();
        if !placement.root_m.iter().all(|x| x.is_finite())
            || !placement.direction.iter().all(|x| x.is_finite())
            || !square.is_finite()
            || (square - 1.0).abs() > 128.0 * f64::EPSILON
        {
            return Err(DepositionError::InvalidInput(
                "finite fin root and physical unit axis required",
            ));
        }
        let mut cells_seen = BTreeSet::from([input.base_cell, input.fluid_cell]);
        let mut links_seen = BTreeSet::new();
        for id in std::iter::once(input.root_link)
            .chain(input.axial_links.iter().copied())
            .chain(input.side_links.iter().copied())
            .chain(input.tip_link)
        {
            checkpoint()?;
            if !links_seen.insert(id) {
                return Err(DepositionError::InvalidInput("duplicate fin link identity"));
            }
        }
        for (i, cell) in input.inventory.iter().enumerate() {
            checkpoint()?;
            let a = input.axial_edges_m[i];
            let b = input.axial_edges_m[i + 1];
            if !nonnegative(a)
                || !positive(b)
                || b <= a
                || !cells_seen.insert(cell.id)
                || !positive(cell.mass_kg)
                || !cell.specific_h_j_kg.is_finite()
                || !positive(cell.boundary.ambient_k)
                || cell.boundary.htc_w_m2_k != 0.0
                || cell.boundary.reservoir.is_some()
                || !cell.boundary.source_w.is_finite()
            {
                return Err(DepositionError::InvalidInput(
                    "fin cells need unique IDs, positive masses, finite stored energies and disjoint finite-fluid boundaries",
                ));
            }
        }
        let area = input.geometry.cross_section_m2();
        let perimeter = input.geometry.exposed_perimeter_m();
        let mut cells = Vec::with_capacity(n);
        let mut links = Vec::with_capacity(link_count);
        let mut conductive_links = Vec::with_capacity(n);
        let mut halves = Vec::with_capacity(n);
        for (i, inventory) in input.inventory.iter().enumerate() {
            checkpoint()?;
            let a = input.axial_edges_m[i];
            let b = input.axial_edges_m[i + 1];
            let dx = finite_positive(b - a)?;
            let half = finite_positive(0.5 * dx)?;
            let volume = finite_positive(area * dx)?;
            let side_area = finite_positive(perimeter * dx)?;
            // The native owner admits coverage with no rounding allowance, and a
            // cell can be covered on its lateral face and both axial faces. The
            // stored surface is the smallest represented value that holds them.
            let mut surface_area = finite_positive(side_area + 2.0 * area)?;
            for _ in 0..8 {
                if holds_full_cover(surface_area, area, side_area) {
                    break;
                }
                surface_area = finite_positive(surface_area.next_up())?;
            }
            if !holds_full_cover(surface_area, area, side_area) {
                return Err(DepositionError::Arithmetic);
            }
            let at = |x: f64| {
                std::array::from_fn(|axis| {
                    placement.direction[axis].mul_add(x, placement.root_m[axis])
                })
            };
            let start_m: [f64; 3] = at(a);
            let end_m: [f64; 3] = at(b);
            if !start_m.iter().chain(&end_m).all(|x| x.is_finite()) || start_m == end_m {
                return Err(DepositionError::InvalidInput(
                    "fin segment cannot be represented in the supplied spatial frame",
                ));
            }
            halves.push(half);
            cells.push(DepositedCell {
                id: inventory.id,
                segment_id: placement.segment_id,
                body_id: placement.body_id,
                material_id: inventory.material_id,
                mass_kg: inventory.mass_kg,
                specific_h_j_kg: inventory.specific_h_j_kg,
                geometry: CellGeometry {
                    start_m,
                    end_m,
                    volume_m3: volume,
                    surface_area_m2: surface_area,
                },
                boundary: inventory.boundary,
            });
            links.push(ThermalLink {
                id: input.side_links[i],
                a: inventory.id,
                b: input.fluid_cell,
                area_m2: side_area,
                conductance_w_k: 0.0,
            });
        }
        links.push(ThermalLink {
            id: input.root_link,
            a: input.base_cell,
            b: cells[0].id,
            area_m2: area,
            conductance_w_k: 0.0,
        });
        conductive_links.push(ConductiveLink {
            link_id: input.root_link,
            half_distances_m: [input.base_half_distance_m, halves[0]],
            contact_resistance_m2_k_w: input.root_contact_resistance_m2_k_w,
        });
        for (i, &id) in input.axial_links.iter().enumerate() {
            checkpoint()?;
            links.push(ThermalLink {
                id,
                a: cells[i].id,
                b: cells[i + 1].id,
                area_m2: area,
                conductance_w_k: 0.0,
            });
            conductive_links.push(ConductiveLink {
                link_id: id,
                half_distances_m: [halves[i], halves[i + 1]],
                contact_resistance_m2_k_w: 0.0,
            });
        }
        if let Some(id) = input.tip_link {
            links.push(ThermalLink {
                id,
                a: cells[n - 1].id,
                b: input.fluid_cell,
                area_m2: area,
                conductance_w_k: 0.0,
            });
        }
        checkpoint()?;
        Ok(Self {
            geometry: input.geometry,
            axial_edges_m: input.axial_edges_m.to_vec(),
            base_cell: input.base_cell,
            fluid_cell: input.fluid_cell,
            cells,
            links,
            conductive_links,
            side_links: input.side_links.to_vec(),
            tip_link: input.tip_link,
        })
    }

    pub const fn geometry(&self) -> StraightFinGeometry {
        self.geometry
    }
    pub fn axial_edges_m(&self) -> &[f64] {
        &self.axial_edges_m
    }
    pub const fn base_cell(&self) -> CellId {
        self.base_cell
    }
    pub const fn fluid_cell(&self) -> CellId {
        self.fluid_cell
    }
    /// Native fin cells in axial order; no base or fluid inventory is included.
    pub fn cells(&self) -> &[DepositedCell] {
        &self.cells
    }
    pub fn links(&self) -> &[ThermalLink] {
        &self.links
    }
    /// Root first, then internal links in increasing axial position.
    pub fn conductive_links(&self) -> &[ConductiveLink] {
        &self.conductive_links
    }
    pub fn root_link(&self) -> u64 {
        self.conductive_links[0].link_id
    }
    pub fn side_links(&self) -> &[u64] {
        &self.side_links
    }
    pub const fn tip_link(&self) -> Option<u64> {
        self.tip_link
    }

    /// Native births followed by native contacts, all at one authored time.
    /// Base and fluid must already exist by the first contact. Insert these
    /// events into the caller's ordered schedule; its original timeline admits
    /// all calorics, IDs, surface coverage, storage and event budgets.
    pub fn events_at(
        &self,
        time_s: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<TimedEvent>, DepositionError> {
        checkpoint()?;
        if !nonnegative(time_s) {
            return Err(DepositionError::InvalidInput(
                "finite nonnegative fin activation time required",
            ));
        }
        let count = self
            .cells
            .len()
            .checked_add(self.links.len())
            .ok_or(DepositionError::Arithmetic)?;
        let mut events = Vec::with_capacity(count);
        for cell in &self.cells {
            checkpoint()?;
            events.push(TimedEvent {
                time_s,
                event: Event::Birth(cell.clone()),
            });
        }
        for &link in &self.links {
            checkpoint()?;
            events.push(TimedEvent {
                time_s,
                event: Event::Contact(link),
            });
        }
        checkpoint()?;
        Ok(events)
    }
}

#[cfg(test)]
mod tests;
