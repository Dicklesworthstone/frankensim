//! Couple the equilibrium wall law to the original reciprocal face flux.
//! The correction replaces the fluid half-cell resistance, not the solid
//! resistance, and stays in series with any independently declared contact.

use std::collections::BTreeMap;
use super::{WallHeatConfig, WallHeatTransfer, equilibrium_wall_heat_transfer, invalid};
use crate::conjugate::{
    ChtError, domain::{FluidProperties, SolidMaterial, VoxelDomain},
    energy::{EnergyConfig, EnergySolution, FaceContactResistance, ThermalSetup, solve_energy},
    simple::FvFlow, poll,
};
use crate::d3q19::Face3;
use fs_exec::CancelGate;

/// One explicitly selected stationary, impermeable, internal fluid/solid wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallHeatFace {
    /// Fluid cell adjacent to the physical wall.
    pub fluid_cell: usize,
    /// Local face of that fluid cell, pointing toward the solid.
    pub face: Face3,
}

/// Thermal coefficients used at a wall; no temperature is frozen here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallHeatFaceReport {
    /// Selected fluid cell and wall orientation.
    pub wall: WallHeatFace,
    /// Adjacent native solid cell; it keeps its own material and storage.
    pub solid_cell: usize,
    /// Actual wall-normal half width and tangential-speed wall law.
    pub transfer: WallHeatTransfer,
    /// Native half-cell resistance before wall integration, m^2 K/W.
    pub native_fluid_resistance_m2_k_w: f64,
    /// Added unresolved-sublayer resistance only, m^2 K/W.
    pub correction_m2_k_w: f64,
    /// Independently supplied physical contact resistance, m^2 K/W.
    pub physical_contact_m2_k_w: f64,
}

/// One preparation for an accepted flow. Reprepare from the ORIGINAL setup
/// after the flow changes; applying this to its own output is refused.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedWallHeat {
    /// Ordinary conservative energy input, also usable by the transient owner
    /// while this velocity/eddy field is held fixed.
    pub thermal: ThermalSetup,
    /// One report per selected wall, preserving input order.
    pub faces: Vec<WallHeatFaceReport>,
}

/// Joint native energy solve and the wall coefficients it actually consumed.
#[derive(Clone, Debug, PartialEq)]
pub struct WallHeatSolution {
    /// Original temperature, heat-flux and global balance reporting.
    pub energy: EnergySolution,
    /// Wall transfer model, missing sublayer and physical contact separately.
    pub faces: Vec<WallHeatFaceReport>,
}

/// Prepare wall-modeled heat transfer without modifying either input.
///
/// Interior fluid conductivity uses the accepted TANGENT eddy viscosity,
/// `rho cp nu_t / Pr_t`, rather than LVEL's momentum secant. Each selected
/// wall replaces only its fluid half-cell thermal resistance with the
/// integrated equilibrium law. The original reciprocal face term then
/// couples both temperatures implicitly, preserves solid conduction and any
/// physical contact resistance, and transfers equal/opposite watts once.
/// Only tangential cell velocity drives the stationary-wall friction law.
///
/// This explicitly selects an equilibrium thermal wall model; it does not
/// infer or certify how the supplied flow was obtained. It is intended for
/// smooth equilibrium LVEL cooling, not separated or resolved-wall LES.
/// Domain-boundary walls and compact-component links are not replaced.
///
/// # Errors
/// Invalid shapes/values, duplicate walls/contacts, non-internal or non-solid
/// walls, through-flow, conflicting supplied thermal eddy fields, compact
/// links, negative required resistance correction, quadrature or cancellation.
/// An incompatible bulk eddy field cannot be repaired with negative resistance.
pub fn prepare_wall_heat(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    flow: &FvFlow,
    thermal: &ThermalSetup,
    walls: &[WallHeatFace],
    config: WallHeatConfig,
    gate: &CancelGate,
) -> Result<PreparedWallHeat, ChtError> {
    poll(gate)?;
    config.validate()?;
    fluid.validate()?;
    let count = domain.cell_count();
    if flow.field.dims() != domain.dims() || flow.velocity_m_s.len() != count
        || flow.eddy_viscosity_m2_s.len() != count
    {
        return Err(invalid("wall heat requires the accepted flow on this exact cell grid"));
    }
    if walls.is_empty() || walls.len() > count.saturating_mul(3) {
        return Err(invalid("select at least one internal wall, within the Cartesian face budget"));
    }
    if !thermal.eddy_conductivity_w_m_k.is_empty() {
        return Err(invalid("wall heat owns turbulent conductivity; use the original unprepared thermal setup"));
    }
    let scale = fluid.volumetric_heat_capacity() / config.turbulent_prandtl;
    if !scale.is_finite() || scale <= 0.0 {
        return Err(invalid("turbulent heat coefficient is not representable"));
    }
    let mut prepared = thermal.clone();
    prepared.eddy_conductivity_w_m_k = Vec::with_capacity(count);
    for cell in 0..count {
        if cell % 64 == 0 { poll(gate)?; }
        let nu = flow.eddy_viscosity_m2_s[cell];
        if !nu.is_finite() || nu < 0.0
            || !flow.velocity_m_s[cell].iter().all(|v| v.is_finite())
            || (!domain.is_fluid(cell) && nu != 0.0)
        {
            return Err(invalid(format!("cell {cell}: invalid velocity or eddy viscosity")));
        }
        let eddy = scale * nu;
        if !eddy.is_finite() || !(fluid.conductivity_w_m_k + eddy).is_finite() {
            return Err(invalid(format!("cell {cell}: thermal conductivity overflow")));
        }
        prepared.eddy_conductivity_w_m_k.push(eddy);
    }
    // Keep independently authored contacts in their original order. A generated
    // correction on the same face augments that one physical contact once.
    let mut contacts = BTreeMap::new();
    for (index, contact) in thermal.face_contacts.iter().enumerate() {
        if index % 64 == 0 { poll(gate)?; }
        let [a, b] = contact.cells;
        if a >= count || b >= count || a == b || !(0..6).any(|f| domain.neighbor(a, f) == Some(b))
            || !contact.resistance_m2_k_w.is_finite() || contact.resistance_m2_k_w < 0.0
            || contacts.insert((a.min(b), a.max(b)), index).is_some()
        {
            return Err(invalid("original thermal contacts must be finite, adjacent and unique"));
        }
    }
    let mut selected = std::collections::BTreeSet::new();
    let mut reports = Vec::with_capacity(walls.len());
    for &wall in walls {
        poll(gate)?;
        let cell = wall.fluid_cell;
        let face = wall.face as usize;
        if cell >= count || !domain.is_fluid(cell) {
            return Err(invalid("selected wall owner must be a fluid cell"));
        }
        let solid = domain.neighbor(cell, face)
            .filter(|&other| !domain.is_fluid(other))
            .ok_or_else(|| invalid("selected wall must face an internal solid, not an opening or domain face"))?;
        let key = (cell.min(solid), cell.max(solid));
        if !selected.insert(key) {
            return Err(invalid("duplicate selected wall"));
        }
        if flow.field.outward(domain, cell, face) != 0.0 {
            return Err(invalid("stationary thermal walls cannot carry normal through-flow"));
        }
        let position = domain.coords(solid);
        if thermal.compact_components.iter().any(|part|
            (0..3).all(|axis| position[axis] >= part.lo[axis] && position[axis] < part.hi[axis]))
        {
            return Err(invalid("compact-component links own this wall's thermal transport"));
        }
        let axis = face / 2;
        let velocity = flow.velocity_m_s[cell];
        let (a, b) = (velocity[(axis + 1) % 3].abs(), velocity[(axis + 2) % 3].abs());
        let magnitude = a.max(b);
        let speed = if magnitude == 0.0 { 0.0 } else {
            magnitude * ((a / magnitude).powi(2) + (b / magnitude).powi(2)).sqrt()
        };
        let distance = 0.5 * domain.widths(cell)[axis];
        let transfer = equilibrium_wall_heat_transfer(speed, distance, fluid, config, gate)?;
        let native = distance / (fluid.conductivity_w_m_k + prepared.eddy_conductivity_w_m_k[cell]);
        let correction = transfer.resistance_m2_k_w - native;
        let roundoff = 64.0 * f64::EPSILON * transfer.resistance_m2_k_w.max(native);
        if !native.is_finite() || native <= 0.0 || correction < -roundoff {
            return Err(invalid(format!("cell {cell}: bulk transport is incompatible with a nonnegative wall correction")));
        }
        // At most roundoff is removed, never an unresolved model discrepancy.
        let correction = correction.max(0.0);
        let physical = contacts.get(&key).map_or(0.0, |&index|
            thermal.face_contacts[index].resistance_m2_k_w);
        let combined = physical + correction;
        if !combined.is_finite() {
            return Err(invalid("combined contact and thermal wall resistance overflow"));
        }
        if let Some(&index) = contacts.get(&key) {
            prepared.face_contacts[index].resistance_m2_k_w = combined;
        } else if correction > 0.0 {
            prepared.face_contacts.push(FaceContactResistance {
                cells: [cell, solid], resistance_m2_k_w: correction,
            });
        }
        reports.push(WallHeatFaceReport { wall, solid_cell: solid, transfer,
            native_fluid_resistance_m2_k_w: native, correction_m2_k_w: correction,
            physical_contact_m2_k_w: physical });
    }
    poll(gate)?;
    Ok(PreparedWallHeat { thermal: prepared, faces: reports })
}

/// Solve the original joint fluid/solid energy equation with these wall laws.
/// Solid temperatures are unknowns, not fixed sinks; the original observers
/// include the same face resistance used in assembly and global accounting.
#[allow(clippy::too_many_arguments)]
pub fn solve_energy_with_wall_heat(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FvFlow,
    thermal: &ThermalSetup,
    walls: &[WallHeatFace],
    wall_config: WallHeatConfig,
    energy_config: &EnergyConfig,
    gate: &CancelGate,
) -> Result<WallHeatSolution, ChtError> {
    let prepared = prepare_wall_heat(domain, fluid, flow, thermal, walls, wall_config, gate)?;
    let energy = solve_energy(domain, fluid, solids, &flow.field,
        &prepared.thermal, energy_config, gate)?;
    poll(gate)?;
    Ok(WallHeatSolution { energy, faces: prepared.faces })
}

#[cfg(test)]
mod tests;
