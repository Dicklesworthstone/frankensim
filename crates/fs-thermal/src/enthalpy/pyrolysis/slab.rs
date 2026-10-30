//! Conservative moving geometry for a homogeneous reactive condensed slab.
//!
//! The original reactive cells conduct heat and release volatiles. At each
//! geometry interval their retained mass determines length through the supplied
//! incompressible component densities. The back face stays fixed. A limited linear
//! reconstruction of specific enthalpy in retained-mass coordinate is integrated
//! onto equal-mass cells; mass and enthalpy, never temperature, are transferred.
//! Original calorics invert the result and original series links use the new
//! lengths. Cell IDs are spatial slots after rezoning, not material particles.
//! Only the two through-thickness faces exchange external heat; lateral faces
//! are excluded/adiabatic in this one-dimensional slab model.
//!
//! The density/length relation follows the single-material limit of the FDS
//! Technical Reference Guide, "Shrinking and Swelling Materials":
//! <https://github.com/firemodels/fds/blob/a26081f3b24915552fab7ebf078185af817d9043/Manuals/FDS_Technical_Reference_Guide/Solid_Chapter.tex>.
//! Fixed-back anchoring and the conservative limited reconstruction are explicit
//! choices of this implementation, not an exact port of FDS. No material is
//! discarded or retained as a numerical seed. [`RecedingSlab::with_char`] adds
//! retained char with separate calorics, density and conductivity; its explicit
//! additive-volume law is not the general FDS shrinking/swelling algorithm.
//! Pore transport, char oxidation, mechanical work,
//! thermal expansion, gas-domain displacement and complete burnout are separate
//! models. Geometry is frozen during each implicit interval; the composition
//! remains first order and remapping can smooth a temperature profile.

use super::implicit::{ImplicitError, ImplicitSolve, ImplicitWorkspace};
use super::{
    Boundary, CellId, DepositionError, PyrolysisAudit, PyrolysisControls, PyrolysisSimulation, Sum,
    VolatileRelease, empty_thermal_audit, finite, nonnegative, positive,
};
use crate::enthalpy::series_conductance_w_k;
use std::collections::BTreeMap;

/// Bounded cell axis, matching the original dense implicit reactive owner.
pub const MAX_RECEDING_SLAB_CELLS: usize = 128;

/// Explicit homogeneous condensed properties on the source's admitted domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabGeometry<'a> {
    /// Constant cross-section normal to the original cells' x axis, m².
    pub area_m2: f64,
    /// Virgin condensed density, kg/m³; no evolving porosity or thermal expansion.
    pub density_kg_m3: f64,
    /// Constant virgin through-thickness conductivity, W/(m K).
    pub conductivity_w_m_k: f64,
    /// Density/conductivity and incompressible shrinking-model provenance.
    pub source: &'a str,
}

/// Explicit effective-conductivity closure for the two condensed components.
/// Both use their actual volume fractions, never their mass fractions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConductivityMixing {
    /// Volume-weighted arithmetic mean: components parallel to heat flow.
    Parallel,
    /// Volume-weighted harmonic mean: components in series along heat flow.
    Series,
}

/// Caller-sourced char properties for the original single-step reaction.
///
/// The retained volume is `m_virgin/rho_virgin + m_char/rho_char`. Densities
/// may be effective bulk values, but this does not evolve porosity, pore-gas
/// pressure or permeability. Expansion gives negative signed recession.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CharProperties<'a> {
    /// Constant retained-char density, kg/m³.
    pub density_kg_m3: f64,
    /// Constant char through-thickness conductivity, W/(m K).
    pub conductivity_w_m_k: f64,
    /// Declared unresolved component arrangement for heat transport.
    pub conductivity_mixing: ConductivityMixing,
    /// Density, conductivity and mixing-closure applicability/provenance.
    pub source: &'a str,
}

/// Geometry/time and numerical accounting gates, not physical error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabSolve {
    /// Maximum accepted geometry intervals in one call, 1 through 65536.
    pub maximum_geometry_steps: usize,
    /// Maximum whole-interval halvings for excessive mass/volume change, 0..32.
    pub maximum_halvings: usize,
    /// Largest admitted loss from any one cell in a frozen-geometry interval.
    /// Must lie in (0, 1); the default keeps at least 75% until rezoning.
    pub maximum_cell_mass_loss_fraction: f64,
    /// Largest absolute fractional cell-volume change in a frozen interval.
    /// Applies to shrinkage and swelling; must lie in (0, 1).
    pub maximum_cell_volume_change_fraction: f64,
    /// Smallest positive attempted geometry interval, s.
    pub minimum_interval_s: f64,
    /// Whole-call sum of absolute remap mass residuals, kg. Also bounds the
    /// numerical mass equivalent of rounded geometric cell lengths.
    pub mass_tolerance_kg: f64,
    /// Whole-call sum of absolute remap energy residuals, J.
    pub energy_tolerance_j: f64,
}
impl Default for SlabSolve {
    fn default() -> Self {
        Self {
            maximum_geometry_steps: 8192,
            maximum_halvings: 24,
            maximum_cell_mass_loss_fraction: 0.25,
            maximum_cell_volume_change_fraction: 0.25,
            minimum_interval_s: 1e-10,
            mass_tolerance_kg: 1e-12,
            energy_tolerance_j: 1e-6,
        }
    }
}

/// Original physical accounts plus separate conservative-remap diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabReport {
    /// Original thermal/reaction account, accumulated over accepted intervals.
    /// Actual remap storage changes and residuals are included once. Release
    /// integrals are unchanged; their IDs name the spatial slots at production.
    pub pyrolysis: PyrolysisAudit,
    /// Accepted frozen-geometry intervals.
    pub geometry_steps: usize,
    /// Discarded complete intervals whose local mass/volume change was too large.
    pub rejected_geometry_attempts: usize,
    /// Intervals that actually changed geometry/reconstructed the mesh.
    pub remaps: usize,
    /// Largest accepted local fractional mass loss before remapping.
    pub maximum_cell_mass_loss_fraction: f64,
    /// Largest accepted absolute component-volume change before remapping.
    pub maximum_cell_volume_change_fraction: f64,
    /// Signed remap mass change, kg; not an emitted or retained source.
    pub remap_mass_residual_kg: f64,
    /// Signed remap virgin-mass change, kg; not chemical conversion.
    pub remap_virgin_mass_residual_kg: f64,
    /// Signed remap char-mass change, kg; not chemical conversion.
    pub remap_char_mass_residual_kg: f64,
    /// Sum of absolute virgin and char remap defects across all intervals, kg.
    pub absolute_remap_component_mass_residual_kg: f64,
    /// Signed remap energy change, J; not a heat source.
    pub remap_energy_residual_j: f64,
    /// Sum of absolute remap mass changes, kg.
    pub absolute_remap_mass_residual_kg: f64,
    /// Sum of absolute remap energy changes, J.
    pub absolute_remap_energy_residual_j: f64,
    /// Sum of overlap mass times the absolute donor/receiver enthalpy difference,
    /// J. For char, component formation references are removed first. This
    /// records redistribution/smoothing, not a bound on PDE error.
    pub remap_enthalpy_redistribution_j: f64,
    /// Largest rounded-volume error times the largest component density, kg.
    pub maximum_geometry_mass_residual_kg: f64,
    /// Actual final distance between front and fixed back, m.
    pub thickness_m: f64,
    /// Actual front motion from its initially admitted position, m.
    pub recession_m: f64,
}

/// Source/domain/numerical refusal. No original state is changed on failure.
#[derive(Clone, Debug, PartialEq)]
pub enum SlabError {
    /// Original geometry, material, arithmetic or balance refusal.
    Thermal(DepositionError),
    /// Original coupled conduction/reaction solve refused.
    Implicit(ImplicitError),
    /// Bounded geometry intervals or halvings exhausted.
    GeometrySteps,
    /// Cancellation publishes no solid or geometry candidate.
    Cancelled,
}
impl core::fmt::Display for SlabError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "receding pyrolysis slab: {self:?}")
    }
}
impl std::error::Error for SlabError {}
impl From<DepositionError> for SlabError {
    fn from(e: DepositionError) -> Self {
        if e == DepositionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermal(e)
        }
    }
}
impl From<ImplicitError> for SlabError {
    fn from(e: ImplicitError) -> Self {
        if matches!(e, ImplicitError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Implicit(e)
        }
    }
}

/// Original reactive solid with one explicit, fixed-back moving geometry.
#[derive(Clone)]
pub struct RecedingSlab<'a> {
    solid: PyrolysisSimulation<'a>,
    geometry: SlabGeometry<'a>,
    char: Option<CharProperties<'a>>,
    ids: Vec<CellId>,
    link_ids: Vec<u64>,
    initial_front_m: f64,
    back_m: f64,
}

fn close(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 512.0 * f64::EPSILON * a.abs().max(b.abs())
}

impl<'a> RecedingSlab<'a> {
    /// Adopt the ORIGINAL activated x-aligned slab without changing its state.
    /// This constructor admits only zero char yield. Current density and each
    /// series link must agree with the explicit
    /// properties. Only the two outer slots may have external boundaries.
    pub fn new(
        solid: PyrolysisSimulation<'a>,
        geometry: SlabGeometry<'a>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SlabError> {
        Self::adopt(solid, geometry, None, checkpoint)
    }

    /// Adopt a slab with retained char and an explicit additive-volume closure.
    /// The original material must supply nonzero char yield and char calorics.
    /// All initial component volumes and links must match the supplied closure.
    pub fn with_char(
        solid: PyrolysisSimulation<'a>,
        geometry: SlabGeometry<'a>,
        char: CharProperties<'a>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SlabError> {
        Self::adopt(solid, geometry, Some(char), checkpoint)
    }

    fn adopt(
        solid: PyrolysisSimulation<'a>,
        geometry: SlabGeometry<'a>,
        char: Option<CharProperties<'a>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SlabError> {
        checkpoint()?;
        let n = solid.thermal.cells.len();
        if !(1..=MAX_RECEDING_SLAB_CELLS).contains(&n)
            || !positive(geometry.area_m2)
            || !positive(geometry.density_kg_m3)
            || !positive(geometry.conductivity_w_m_k)
            || geometry.source.trim().is_empty()
            || char.is_some_and(|c| {
                !positive(c.density_kg_m3)
                    || !positive(c.conductivity_w_m_k)
                    || c.source.trim().is_empty()
            })
            || !solid.thermal.covered_areas.is_empty()
            || solid.thermal.links.len() != n - 1
        {
            return Err(DepositionError::InvalidInput(
                "receding slab requires bounded cells, explicit positive properties and a plain chain",
            ).into());
        }
        let mut cells: Vec<_> = solid.thermal.cells.values().collect();
        cells.sort_by(|a, b| a.geometry.start_m[0].total_cmp(&b.geometry.start_m[0]));
        let first = cells[0];
        let material = solid.materials[&first.material_id];
        if (char.is_none() && material.char_mass_yield != 0.0)
            || (char.is_some()
                && (material.char_mass_yield <= 0.0 || material.char_caloric.is_none()))
        {
            return Err(DepositionError::InvalidInput(
                "retained char requires explicit char geometry and original char calorics",
            )
            .into());
        }
        let mut ids = Vec::with_capacity(n);
        for (i, cell) in cells.iter().enumerate() {
            checkpoint()?;
            let g = cell.geometry;
            let dx = g.end_m[0] - g.start_m[0];
            let inventory = solid.inventories[&cell.id];
            if cell.material_id != first.material_id
                || cell.body_id != first.body_id
                || cell.segment_id != first.segment_id
                || !close(
                    inventory.virgin_mass_kg + inventory.char_mass_kg,
                    cell.mass_kg,
                )
                || (char.is_none() && inventory.char_mass_kg != 0.0)
                || !positive(dx)
                || g.start_m[1..] != first.geometry.start_m[1..]
                || g.end_m[1..] != first.geometry.start_m[1..]
                || (i > 0 && cells[i - 1].geometry.end_m != g.start_m)
                || !close(g.volume_m3, geometry.area_m2 * dx)
                || !close(g.surface_area_m2, 2.0 * geometry.area_m2)
                || !close(component_volume(inventory, geometry, char)?, g.volume_m3)
                || (i > 0 && i + 1 < n && !insulated(cell.boundary))
            {
                return Err(DepositionError::InvalidInput(
                    "receding slab requires contiguous homogeneous material cells and outer-face sources",
                ).into());
            }
            ids.push(cell.id);
        }
        let mut link_ids = Vec::with_capacity(n - 1);
        for pair in cells.windows(2) {
            checkpoint()?;
            let mut links = solid.thermal.links.values().filter(|link| {
                (link.a == pair[0].id && link.b == pair[1].id)
                    || (link.b == pair[0].id && link.a == pair[1].id)
            });
            let link = links
                .next()
                .ok_or(DepositionError::InvalidInput("missing slab series link"))?;
            let dx = [pair[0], pair[1]]
                .map(|cell| 0.5 * (cell.geometry.end_m[0] - cell.geometry.start_m[0]));
            let k = [
                component_conductivity(solid.inventories[&pair[0].id], geometry, char)?,
                component_conductivity(solid.inventories[&pair[1].id], geometry, char)?,
            ];
            let expected = series_conductance_w_k(geometry.area_m2, dx, k)?;
            if links.next().is_some()
                || !close(link.area_m2, geometry.area_m2)
                || !close(link.conductance_w_k, expected)
            {
                return Err(DepositionError::InvalidInput(
                    "slab link differs from original series law",
                )
                .into());
            }
            link_ids.push(link.id);
        }
        let initial_front_m = cells[0].geometry.start_m[0];
        let back_m = cells[n - 1].geometry.end_m[0];
        checkpoint()?;
        Ok(Self {
            solid,
            geometry,
            char,
            ids,
            link_ids,
            initial_front_m,
            back_m,
        })
    }

    /// Original accepted thermal/reaction owner, with its current geometry.
    pub fn state(&self) -> &PyrolysisSimulation<'a> {
        &self.solid
    }
    /// Retained source and virgin material/geometric assumptions.
    pub const fn geometry(&self) -> SlabGeometry<'a> {
        self.geometry
    }
    /// Explicit retained-char closure, if this is a charring slab.
    pub const fn char_properties(&self) -> Option<CharProperties<'a>> {
        self.char
    }
    /// Current effective conductivity from the declared component-volume mix.
    pub fn cell_conductivity_w_m_k(&self, id: CellId) -> Result<f64, SlabError> {
        let inventory = *self
            .solid
            .inventories
            .get(&id)
            .ok_or(DepositionError::InvalidInput("unknown receding slab cell"))?;
        Ok(component_conductivity(inventory, self.geometry, self.char)?)
    }
    /// Actual accepted time, s.
    pub fn time_s(&self) -> f64 {
        self.solid.time_s()
    }
    /// Actual retained solid mass, kg.
    pub fn mass_kg(&self) -> f64 {
        self.solid.mass_kg()
    }
    /// Original physical solid enthalpy storage, J.
    pub fn energy_j(&self) -> f64 {
        self.solid.energy_j()
    }
    /// Cell identity of the moving front spatial slot.
    pub fn front_cell_id(&self) -> CellId {
        self.ids[0]
    }
    /// Cell identity of the fixed back spatial slot.
    pub fn back_cell_id(&self) -> CellId {
        self.ids[self.ids.len() - 1]
    }
    /// Actual accepted front coordinate, m.
    pub fn front_m(&self) -> f64 {
        self.solid.cells()[&self.ids[0]].geometry.start_m[0]
    }
    /// Fixed accepted back coordinate, m.
    pub const fn back_m(&self) -> f64 {
        self.back_m
    }
    /// Actual accepted thickness, m.
    pub fn thickness_m(&self) -> f64 {
        self.back_m - self.front_m()
    }
    /// Front displacement from the initially admitted state, m; negative for swelling.
    pub fn recession_m(&self) -> f64 {
        self.front_m() - self.initial_front_m
    }
    /// Actual front half-cell conductive resistance through original geometry.
    pub fn front_half_cell_conductance_w_k(&self) -> Result<f64, SlabError> {
        let cell = &self.solid.cells()[&self.front_cell_id()];
        let dx = cell.geometry.end_m[0] - cell.geometry.start_m[0];
        Ok(series_conductance_w_k(
            self.geometry.area_m2,
            [0.25 * dx; 2],
            [self.cell_conductivity_w_m_k(cell.id)?; 2],
        )?)
    }

    /// Change an outer-face boundary through original transactional admission.
    pub fn set_boundary(
        &mut self,
        id: CellId,
        boundary: Boundary,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SlabError> {
        if id != self.front_cell_id() && id != self.back_cell_id() {
            return Err(DepositionError::InvalidInput(
                "receding slab boundary must be an outer slot",
            )
            .into());
        }
        self.solid.set_boundary(id, boundary, checkpoint)?;
        Ok(())
    }

    /// Advance original implicit heat/reaction, then conservative geometry.
    ///
    /// The maximum interval is the original `maximum_dt_s`. Excessive actual
    /// local mass loss or component-volume change retries the entire interval
    /// from accepted state; releases
    /// from rejected attempts never enter the returned account. Material and
    /// numerical errors from the original implicit owner propagate unchanged.
    pub fn advance_implicit_to(
        &self,
        target_s: f64,
        controls: PyrolysisControls,
        implicit: ImplicitSolve,
        geometry_controls: SlabSolve,
        workspace: &mut ImplicitWorkspace,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, SlabReport), SlabError> {
        checkpoint()?;
        admit_controls(target_s, self.time_s(), controls, geometry_controls)?;
        let mut next = self.clone();
        let mut report = self.empty_report();
        while next.time_s() < target_s {
            checkpoint()?;
            if report.geometry_steps >= geometry_controls.maximum_geometry_steps {
                return Err(SlabError::GeometrySteps);
            }
            let remaining = target_s - next.time_s();
            let mut dt = remaining.min(controls.maximum_dt_s);
            let mut accepted = None;
            for attempt in 0..=geometry_controls.maximum_halvings {
                checkpoint()?;
                if !positive(dt) || dt < geometry_controls.minimum_interval_s {
                    return Err(SlabError::GeometrySteps);
                }
                let endpoint = if dt == remaining {
                    target_s
                } else {
                    next.time_s() + dt
                };
                if endpoint <= next.time_s() || !endpoint.is_finite() {
                    return Err(DepositionError::Arithmetic.into());
                }
                let (candidate, implicit_report) = next.solid.advance_implicit_to(
                    endpoint,
                    controls,
                    implicit,
                    workspace,
                    &mut checkpoint,
                )?;
                let maximum_loss = next
                    .ids
                    .iter()
                    .map(|id| {
                        (next.solid.cells()[id].mass_kg - candidate.cells()[id].mass_kg)
                            / next.solid.cells()[id].mass_kg
                    })
                    .fold(0.0, f64::max);
                let mut maximum_volume_change = 0.0_f64;
                for &id in &next.ids {
                    let old =
                        component_volume(next.solid.inventories[&id], next.geometry, next.char)?;
                    let new =
                        component_volume(candidate.inventories[&id], next.geometry, next.char)?;
                    maximum_volume_change = maximum_volume_change.max(((new - old) / old).abs());
                }
                if maximum_loss <= geometry_controls.maximum_cell_mass_loss_fraction
                    && maximum_volume_change
                        <= geometry_controls.maximum_cell_volume_change_fraction
                {
                    accepted = Some((
                        candidate,
                        implicit_report.pyrolysis,
                        maximum_loss,
                        maximum_volume_change,
                    ));
                    break;
                }
                report.rejected_geometry_attempts += 1;
                if attempt == geometry_controls.maximum_halvings {
                    return Err(SlabError::GeometrySteps);
                }
                dt *= 0.5;
            }
            let (candidate, audit, maximum_loss, maximum_volume_change) =
                accepted.ok_or(SlabError::GeometrySteps)?;
            let changed = next.ids.iter().any(|id| {
                next.solid.inventories[id].virgin_mass_kg
                    != candidate.inventories[id].virgin_mass_kg
                    || next.solid.inventories[id].char_mass_kg
                        != candidate.inventories[id].char_mass_kg
            });
            next.solid = candidate;
            accumulate(&mut report.pyrolysis, audit)?;
            if changed {
                let remap = next.rezone(geometry_controls, &mut checkpoint)?;
                report.remaps += 1;
                report.remap_mass_residual_kg =
                    finite(report.remap_mass_residual_kg + remap.mass_kg)?;
                report.remap_virgin_mass_residual_kg =
                    finite(report.remap_virgin_mass_residual_kg + remap.virgin_mass_kg)?;
                report.remap_char_mass_residual_kg =
                    finite(report.remap_char_mass_residual_kg + remap.char_mass_kg)?;
                report.absolute_remap_component_mass_residual_kg = finite(
                    report.absolute_remap_component_mass_residual_kg
                        + remap.virgin_mass_kg.abs()
                        + remap.char_mass_kg.abs(),
                )?;
                report.remap_energy_residual_j =
                    finite(report.remap_energy_residual_j + remap.energy_j)?;
                report.absolute_remap_mass_residual_kg =
                    finite(report.absolute_remap_mass_residual_kg + remap.mass_kg.abs())?;
                report.absolute_remap_energy_residual_j =
                    finite(report.absolute_remap_energy_residual_j + remap.energy_j.abs())?;
                report.remap_enthalpy_redistribution_j =
                    finite(report.remap_enthalpy_redistribution_j + remap.redistribution_j)?;
                report.maximum_geometry_mass_residual_kg = report
                    .maximum_geometry_mass_residual_kg
                    .max(remap.geometry_mass_kg);
                report.pyrolysis.mass_residual_kg =
                    finite(report.pyrolysis.mass_residual_kg + remap.mass_kg)?;
                report.pyrolysis.thermal.stored_energy_change_j =
                    finite(report.pyrolysis.thermal.stored_energy_change_j + remap.energy_j)?;
                report.pyrolysis.thermal.residual_j =
                    finite(report.pyrolysis.thermal.residual_j + remap.energy_j)?;
                report.pyrolysis.absolute_energy_residual_j =
                    finite(report.pyrolysis.absolute_energy_residual_j + remap.energy_j.abs())?;
                report.pyrolysis.temperature_evaluations += remap.temperature_evaluations;
            }
            report.geometry_steps += 1;
            report.maximum_cell_mass_loss_fraction =
                report.maximum_cell_mass_loss_fraction.max(maximum_loss);
            report.maximum_cell_volume_change_fraction = report
                .maximum_cell_volume_change_fraction
                .max(maximum_volume_change);
            if report.absolute_remap_mass_residual_kg > geometry_controls.mass_tolerance_kg
                || report.absolute_remap_component_mass_residual_kg
                    > geometry_controls.mass_tolerance_kg
                || report.absolute_remap_energy_residual_j > geometry_controls.energy_tolerance_j
                || report.pyrolysis.mass_residual_kg.abs()
                    + report.pyrolysis.kinetic_mass_roundoff_kg
                    > controls.mass_tolerance_kg
                || report.pyrolysis.absolute_energy_residual_j
                    + report.pyrolysis.kinetic_energy_roundoff_j
                    > controls.total_energy_tolerance_j
            {
                return Err(DepositionError::Budget(
                    "receding slab cumulative conservation discrepancy",
                )
                .into());
            }
        }
        report.thickness_m = next.thickness_m();
        report.recession_m = next.recession_m();
        checkpoint()?;
        Ok((next, report))
    }

    fn empty_report(&self) -> SlabReport {
        SlabReport {
            pyrolysis: PyrolysisAudit {
                thermal: empty_thermal_audit(),
                releases: self
                    .ids
                    .iter()
                    .map(|&id| (id, VolatileRelease::default()))
                    .collect(),
                volatile_mass_kg: 0.0,
                volatile_enthalpy_j: 0.0,
                mass_residual_kg: 0.0,
                maximum_cell_energy_residual_j: 0.0,
                absolute_energy_residual_j: 0.0,
                kinetic_mass_roundoff_kg: 0.0,
                kinetic_energy_roundoff_j: 0.0,
                temperature_evaluations: 0,
            },
            geometry_steps: 0,
            rejected_geometry_attempts: 0,
            remaps: 0,
            maximum_cell_mass_loss_fraction: 0.0,
            maximum_cell_volume_change_fraction: 0.0,
            remap_mass_residual_kg: 0.0,
            remap_virgin_mass_residual_kg: 0.0,
            remap_char_mass_residual_kg: 0.0,
            absolute_remap_component_mass_residual_kg: 0.0,
            remap_energy_residual_j: 0.0,
            absolute_remap_mass_residual_kg: 0.0,
            absolute_remap_energy_residual_j: 0.0,
            remap_enthalpy_redistribution_j: 0.0,
            maximum_geometry_mass_residual_kg: 0.0,
            thickness_m: self.thickness_m(),
            recession_m: self.recession_m(),
        }
    }

    // State is already an unpublished implicit endpoint. Rezone in place,
    // preserving the original cell/link set and every release integral.
    fn rezone(
        &mut self,
        controls: SlabSolve,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Remap, SlabError> {
        if self.char.is_some() {
            return self.rezone_char(controls, checkpoint);
        }
        let n = self.ids.len();
        let mut edges = Vec::with_capacity(n + 1);
        let mut h = Vec::with_capacity(n);
        let mut total = Sum::default();
        edges.push(0.0);
        for &id in &self.ids {
            checkpoint()?;
            let cell = &self.solid.cells()[&id];
            total.add(cell.mass_kg)?;
            if total.value() <= *edges.last().unwrap() {
                return Err(DepositionError::Arithmetic.into());
            }
            edges.push(total.value());
            h.push(cell.specific_h_j_kg);
        }
        let mass_before = total.value();
        let reference_h = h[0];
        let mut old_relative_energy = Sum::default();
        for (i, &id) in self.ids.iter().enumerate() {
            old_relative_energy.add(finite(
                self.solid.cells()[&id].mass_kg * (h[i] - reference_h),
            )?)?;
        }
        let slopes = limited_slopes(&edges, &h)?;
        let mut new_m = Vec::with_capacity(n);
        let mut new_h = Vec::with_capacity(n);
        let mut source = 0usize;
        let mut redistribution = Sum::default();
        for j in 0..n {
            checkpoint()?;
            let left = mass_before * (j as f64 / n as f64);
            let right = if j + 1 == n {
                mass_before
            } else {
                mass_before * ((j + 1) as f64 / n as f64)
            };
            if right <= left {
                return Err(DepositionError::Arithmetic.into());
            }
            while source + 1 < n && edges[source + 1] <= left {
                source += 1;
            }
            let first = source;
            let mut target_mass = Sum::default();
            let mut target_energy = Sum::default();
            let mut last = source;
            while source < n && edges[source] < right {
                let a = left.max(edges[source]);
                let b = right.min(edges[source + 1]);
                if b > a {
                    let overlap = b - a;
                    let center_offset = 0.5 * ((a - edges[source]) + (b - edges[source]))
                        - 0.5 * (edges[source + 1] - edges[source]);
                    let average =
                        finite((h[source] - reference_h) + slopes[source] * center_offset)?;
                    target_mass.add(overlap)?;
                    target_energy.add(finite(overlap * average)?)?;
                    last = source;
                }
                if edges[source + 1] >= right {
                    break;
                }
                source += 1;
            }
            let mass = target_mass.value();
            if !positive(mass) {
                return Err(DepositionError::Arithmetic.into());
            }
            let enthalpy = finite(reference_h + target_energy.value() / mass)?;
            for i in first..=last {
                let overlap = right.min(edges[i + 1]) - left.max(edges[i]);
                if overlap > 0.0 {
                    redistribution.add(finite(overlap * (h[i] - enthalpy).abs())?)?;
                }
            }
            new_m.push(mass);
            new_h.push(enthalpy);
        }
        let mut new_total = Sum::default();
        let mut new_relative_energy = Sum::default();
        for i in 0..n {
            new_total.add(new_m[i])?;
            new_relative_energy.add(finite(new_m[i] * (new_h[i] - reference_h))?)?;
        }
        let mass_residual = finite(new_total.value() - mass_before)?;
        let energy_residual = finite(
            new_relative_energy.value() - old_relative_energy.value() + mass_residual * reference_h,
        )?;
        if mass_residual.abs() > controls.mass_tolerance_kg
            || energy_residual.abs() > controls.energy_tolerance_j
        {
            return Err(DepositionError::Budget("conservative slab remap discrepancy").into());
        }
        let density_area = finite(self.geometry.density_kg_m3 * self.geometry.area_m2)?;
        let mut prefix = Sum::default();
        let mut coordinate = finite(self.back_m - new_total.value() / density_area)?;
        let mut geometry_mass = 0.0_f64;
        for (i, &id) in self.ids.iter().enumerate() {
            checkpoint()?;
            prefix.add(new_m[i])?;
            let end = if i + 1 == n {
                self.back_m
            } else {
                finite(self.back_m - (new_total.value() - prefix.value()) / density_area)?
            };
            if end <= coordinate {
                return Err(DepositionError::Arithmetic.into());
            }
            let cell = self.solid.thermal.cells.get_mut(&id).unwrap();
            let model = self.solid.materials[&cell.material_id];
            let temperature = model.virgin.temperature_k(new_h[i])?;
            model.enthalpies(temperature)?;
            geometry_mass =
                geometry_mass.max(finite(density_area * (end - coordinate) - new_m[i])?.abs());
            cell.mass_kg = new_m[i];
            cell.specific_h_j_kg = new_h[i];
            cell.geometry.start_m[0] = coordinate;
            cell.geometry.end_m[0] = end;
            cell.geometry.volume_m3 = finite(new_m[i] / self.geometry.density_kg_m3)?;
            let inventory = self.solid.inventories.get_mut(&id).unwrap();
            inventory.virgin_mass_kg = new_m[i];
            inventory.temperature_k = temperature;
            coordinate = end;
        }
        if geometry_mass > controls.mass_tolerance_kg {
            return Err(DepositionError::Budget(
                "rounded slab geometry exceeds mass-equivalent budget",
            )
            .into());
        }
        let mut conductances = BTreeMap::new();
        for (i, &link_id) in self.link_ids.iter().enumerate() {
            checkpoint()?;
            let pair = [self.ids[i], self.ids[i + 1]];
            let halves = pair.map(|id| {
                let g = self.solid.thermal.cells[&id].geometry;
                0.5 * (g.end_m[0] - g.start_m[0])
            });
            conductances.insert(
                link_id,
                series_conductance_w_k(
                    self.geometry.area_m2,
                    halves,
                    [self.geometry.conductivity_w_m_k; 2],
                )?,
            );
        }
        self.solid
            .thermal
            .set_step_conductances(&conductances, checkpoint)?;
        checkpoint()?;
        Ok(Remap {
            mass_kg: mass_residual,
            virgin_mass_kg: mass_residual,
            char_mass_kg: 0.0,
            energy_j: energy_residual,
            redistribution_j: redistribution.value(),
            geometry_mass_kg: geometry_mass,
            temperature_evaluations: 0,
        })
    }
}

fn insulated(boundary: Boundary) -> bool {
    boundary.htc_w_m2_k == 0.0 && boundary.reservoir.is_none() && boundary.source_w == 0.0
}
fn admit_controls(
    target: f64,
    initial: f64,
    thermal: PyrolysisControls,
    slab: SlabSolve,
) -> Result<(), DepositionError> {
    if !target.is_finite()
        || target < initial
        || !positive(thermal.maximum_dt_s)
        || !(4..=256).contains(&thermal.maximum_temperature_iterations)
        || !nonnegative(thermal.cell_energy_tolerance_j)
        || !nonnegative(thermal.total_energy_tolerance_j)
        || !nonnegative(thermal.mass_tolerance_kg)
        || !(1..=65_536).contains(&slab.maximum_geometry_steps)
        || slab.maximum_halvings > 32
        || !positive(slab.maximum_cell_mass_loss_fraction)
        || slab.maximum_cell_mass_loss_fraction >= 1.0
        || !positive(slab.maximum_cell_volume_change_fraction)
        || slab.maximum_cell_volume_change_fraction >= 1.0
        || !positive(slab.minimum_interval_s)
        || !nonnegative(slab.mass_tolerance_kg)
        || !nonnegative(slab.energy_tolerance_j)
    {
        Err(DepositionError::InvalidInput(
            "invalid shrinking-slab controls",
        ))
    } else {
        Ok(())
    }
}

struct Remap {
    mass_kg: f64,
    virgin_mass_kg: f64,
    char_mass_kg: f64,
    energy_j: f64,
    redistribution_j: f64,
    geometry_mass_kg: f64,
    temperature_evaluations: usize,
}

fn component_volume(
    inventory: super::PyrolysisCell,
    geometry: SlabGeometry<'_>,
    char: Option<CharProperties<'_>>,
) -> Result<f64, DepositionError> {
    let char_volume = match char {
        Some(c) => inventory.char_mass_kg / c.density_kg_m3,
        None if inventory.char_mass_kg == 0.0 => 0.0,
        None => return Err(DepositionError::InvalidInput("char volume is unspecified")),
    };
    let volume = finite(inventory.virgin_mass_kg / geometry.density_kg_m3 + char_volume)?;
    if !positive(volume) {
        return Err(DepositionError::Arithmetic);
    }
    Ok(volume)
}

fn component_conductivity(
    inventory: super::PyrolysisCell,
    geometry: SlabGeometry<'_>,
    char: Option<CharProperties<'_>>,
) -> Result<f64, DepositionError> {
    let Some(char) = char else {
        return Ok(geometry.conductivity_w_m_k);
    };
    let volume = component_volume(inventory, geometry, Some(char))?;
    let virgin_volume = inventory.virgin_mass_kg / geometry.density_kg_m3;
    let char_volume = inventory.char_mass_kg / char.density_kg_m3;
    let conductivity = finite(match char.conductivity_mixing {
        ConductivityMixing::Parallel => {
            (virgin_volume * geometry.conductivity_w_m_k + char_volume * char.conductivity_w_m_k)
                / volume
        }
        ConductivityMixing::Series => {
            volume
                / (virgin_volume / geometry.conductivity_w_m_k
                    + char_volume / char.conductivity_w_m_k)
        }
    })?;
    if !positive(conductivity) {
        return Err(DepositionError::Arithmetic);
    }
    Ok(conductivity)
}

#[path = "slab/charring.rs"]
mod charring;

// Linear reconstruction is conservative in each donor mass interval. MC
// slopes are additionally limited at BOTH donor edges, including nonuniform
// intervals. The endpoint cells use zero slope rather than inventing ghosts.
fn limited_slopes(edges: &[f64], h: &[f64]) -> Result<Vec<f64>, DepositionError> {
    let n = h.len();
    let mut slopes = vec![0.0; n];
    for i in 1..n.saturating_sub(1) {
        let widths = [
            edges[i] - edges[i - 1],
            edges[i + 1] - edges[i],
            edges[i + 2] - edges[i + 1],
        ];
        let dl = finite((h[i] - h[i - 1]) / (0.5 * (widths[0] + widths[1])))?;
        let dr = finite((h[i + 1] - h[i]) / (0.5 * (widths[1] + widths[2])))?;
        if dl * dr > 0.0 {
            let central =
                finite((h[i + 1] - h[i - 1]) / (0.5 * widths[0] + widths[1] + 0.5 * widths[2]))?;
            let edge_bound = 2.0 * (h[i] - h[i - 1]).abs().min((h[i + 1] - h[i]).abs()) / widths[1];
            let magnitude = central
                .abs()
                .min(2.0 * dl.abs())
                .min(2.0 * dr.abs())
                .min(edge_bound);
            slopes[i] = finite(magnitude.copysign(dl))?;
        }
    }
    Ok(slopes)
}

fn accumulate(total: &mut PyrolysisAudit, step: PyrolysisAudit) -> Result<(), DepositionError> {
    total.thermal.external_heat_j =
        finite(total.thermal.external_heat_j + step.thermal.external_heat_j)?;
    total.thermal.stored_energy_change_j =
        finite(total.thermal.stored_energy_change_j + step.thermal.stored_energy_change_j)?;
    total.thermal.residual_j = finite(total.thermal.residual_j + step.thermal.residual_j)?;
    total.thermal.substeps += step.thermal.substeps;
    total.thermal.events += step.thermal.events;
    total.volatile_mass_kg = finite(total.volatile_mass_kg + step.volatile_mass_kg)?;
    total.volatile_enthalpy_j = finite(total.volatile_enthalpy_j + step.volatile_enthalpy_j)?;
    total.mass_residual_kg = finite(total.mass_residual_kg + step.mass_residual_kg)?;
    total.maximum_cell_energy_residual_j = total
        .maximum_cell_energy_residual_j
        .max(step.maximum_cell_energy_residual_j);
    total.absolute_energy_residual_j =
        finite(total.absolute_energy_residual_j + step.absolute_energy_residual_j)?;
    total.kinetic_mass_roundoff_kg =
        finite(total.kinetic_mass_roundoff_kg + step.kinetic_mass_roundoff_kg)?;
    total.kinetic_energy_roundoff_j =
        finite(total.kinetic_energy_roundoff_j + step.kinetic_energy_roundoff_j)?;
    total.temperature_evaluations += step.temperature_evaluations;
    for (id, release) in step.releases {
        let out = total
            .releases
            .get_mut(&id)
            .ok_or(DepositionError::InvalidInput(
                "slab release changed source slots",
            ))?;
        out.mass_kg = finite(out.mass_kg + release.mass_kg)?;
        out.enthalpy_j = finite(out.enthalpy_j + release.enthalpy_j)?;
        out.consumed_virgin_mass_kg =
            finite(out.consumed_virgin_mass_kg + release.consumed_virgin_mass_kg)?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "slab/tests.rs"]
mod tests;
