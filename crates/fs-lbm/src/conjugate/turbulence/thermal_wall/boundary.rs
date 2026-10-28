//! Exterior thermal walls lowered to the existing conservative energy assembly.
//! Like `boundary_exchange`, delegated faces use native linear exchanges; imposed
//! flux uses native watts. Reports identify these boundary contributions separately.

use std::collections::{BTreeMap, BTreeSet};
use super::{WallHeatConfig, WallHeatFace, WallHeatFaceReport, WallHeatTransfer,
    equilibrium_wall_heat_transfer, invalid, prepare_wall_heat};
use crate::conjugate::{CellSink, ChtError, EnergyConfig, EnergySolution, FluidProperties,
    FvBoundary, FvFlow, SolidMaterial, ThermalFace, ThermalSetup, VoxelDomain, poll,
    solve_energy};
use fs_exec::CancelGate;

/// One impermeable exterior fluid wall with an explicit thermal condition.
/// The original thermal domain face must be Adiabatic, delegating this selected
/// face only. Other cells on that domain face are not silently reconfigured.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryWallHeat {
    /// Fluid cell and its outward domain face.
    pub wall: WallHeatFace,
    /// Temperature at the outer side of `contact_m2_k_w`, external Convective
    /// film, or HeatFlux INTO the fluid. Other ThermalFace variants refuse.
    pub condition: ThermalFace,
    /// Additional physical coating/contact resistance outside the fluid wall
    /// layer. It has no storage and is never inferred from roughness.
    pub contact_m2_k_w: f64,
}

/// Frozen coefficients actually supplied to the native heat solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryWallHeatReport {
    /// Original wall selection and physical thermal condition.
    pub boundary: BoundaryWallHeat,
    /// Physical Cartesian face area, m^2.
    pub area_m2: f64,
    /// Wall law inferred from the tangential WALL-RELATIVE velocity.
    pub transfer: WallHeatTransfer,
    /// Contact plus external film resistance, excluding the fluid layer.
    pub exterior_resistance_m2_k_w: f64,
    /// Native linear exchange coefficient, W/K; zero for prescribed flux.
    pub conductance_w_k: f64,
}

/// Heat and fluid-side surface temperature recovered from the solved cell state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryWallHeatExchange {
    /// The coefficients used by that same solve.
    pub coefficients: BoundaryWallHeatReport,
    /// Signed heat leaving the domain, W (negative for prescribed inward heat).
    pub outward_w: f64,
    /// Fluid-side wall temperature after the modeled fluid-layer drop, K.
    pub surface_temperature_k: f64,
}

/// Internal conjugate walls and exterior walls on one accepted flow image.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedWallHeatBoundaries {
    /// Native energy input. Prepare each new flow from the ORIGINAL setup;
    /// preparing this output again refuses because it already owns eddy heat.
    pub thermal: ThermalSetup,
    /// Native solid/fluid series-resistance coefficients.
    pub internal_faces: Vec<WallHeatFaceReport>,
    /// Exterior coefficients in input order; no temperatures are frozen.
    pub boundary_faces: Vec<BoundaryWallHeatReport>,
    /// Actual quadrature evaluations in this preparation. Identical exterior
    /// speed/distance queries share a per-call cache, without quantization.
    pub quadrature_evaluations: usize,
}

impl PreparedWallHeatBoundaries {
    /// Recover signed boundary heat without advancing either owner. These watts
    /// are ALREADY in native sink/source accounting, never an additional loss.
    pub fn boundary_exchanges(&self, temperature: &[f64], gate: &CancelGate)
        -> Result<Vec<BoundaryWallHeatExchange>, ChtError>
    {
        poll(gate)?;
        if temperature.len() != self.thermal.eddy_conductivity_w_m_k.len() {
            return Err(invalid("boundary temperature field does not match the prepared grid"));
        }
        let mut exchanges = Vec::with_capacity(self.boundary_faces.len());
        for &report in &self.boundary_faces {
            poll(gate)?;
            let t = temperature[report.boundary.wall.fluid_cell];
            if !t.is_finite() || t <= 0.0 {
                return Err(invalid("boundary cell temperature must be finite and positive"));
            }
            let (outward_w, surface_temperature_k) = match report.boundary.condition {
                ThermalFace::HeatFlux(inward) => (
                    -report.area_m2 * inward,
                    t + inward * report.transfer.resistance_m2_k_w,
                ),
                ThermalFace::Temperature(reservoir)
                | ThermalFace::Convective { ambient: reservoir, .. } => {
                    let resistance = report.transfer.resistance_m2_k_w
                        + report.exterior_resistance_m2_k_w;
                    let difference = t - reservoir;
                    // Convex interpolation avoids subtracting almost all of t
                    // across the fluid layer; a prescribed wall is exact.
                    let fraction = report.exterior_resistance_m2_k_w / resistance;
                    (report.conductance_w_k * difference, reservoir + fraction * difference)
                }
                _ => return Err(invalid("prepared boundary has an unsupported condition")),
            };
            if !outward_w.is_finite() || !surface_temperature_k.is_finite()
                || surface_temperature_k <= 0.0
            {
                return Err(invalid("boundary heat or surface temperature is not representable"));
            }
            exchanges.push(BoundaryWallHeatExchange {
                coefficients: report, outward_w, surface_temperature_k,
            });
        }
        poll(gate)?;
        Ok(exchanges)
    }
}

/// A joint steady energy solution with the actual wall exchanges disaggregated.
#[derive(Clone, Debug, PartialEq)]
pub struct WallHeatBoundarySolution {
    /// Original conservative solver result. Exterior temperature/film watts
    /// are included in sink_outflow_w; prescribed inward flux is in source_w.
    pub energy: EnergySolution,
    /// Internal wall coefficients, including physical contacts separately.
    pub internal_faces: Vec<WallHeatFaceReport>,
    /// Exterior heat and fluid-side wall temperatures from the solved state.
    pub boundary_faces: Vec<BoundaryWallHeatExchange>,
}

/// Prepare selected internal and exterior walls using the same tangent eddy
/// conductivity and constant-property wall law. Exterior temperature/film laws
/// use G=A/(R_fluid+R_contact+1/h), with no duplicate native half-cell term.
/// Prescribed flux is independent of resistance; the wall law determines its
/// surface temperature, not its watts. Mechanical faces must be actual Wall
/// conditions; free-slip, openings and normal wall motion refuse.
///
/// This is an equilibrium model, including when refreshed during a transient.
/// It does not extend its applicability to separation, roughness or wall LES.
/// Inputs are unchanged on every refusal and cancellation.
#[allow(clippy::too_many_arguments)]
pub fn prepare_wall_heat_boundaries(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    flow: &FvFlow,
    thermal: &ThermalSetup,
    mechanical_faces: &[FvBoundary; 6],
    internal_walls: &[WallHeatFace],
    boundary_walls: &[BoundaryWallHeat],
    config: WallHeatConfig,
    gate: &CancelGate,
) -> Result<PreparedWallHeatBoundaries, ChtError> {
    poll(gate)?;
    let count = domain.cell_count();
    if internal_walls.is_empty() && boundary_walls.is_empty() {
        return Err(invalid("select at least one internal or exterior wall"));
    }
    if boundary_walls.len() > count.saturating_mul(6) {
        return Err(invalid("exterior wall count exceeds the Cartesian face budget"));
    }
    // One shared bulk-transport owner. An empty internal selection prepares
    // only the tangent conductivity; the exterior coefficients are added below.
    let prepared = prepare_wall_heat(domain, fluid, flow, thermal, internal_walls, config, gate)?;
    let mut next = prepared.thermal;
    if !next.power_w.is_empty() && next.power_w.len() != count {
        return Err(invalid("original heat-source array has the wrong shape"));
    }
    let mut evaluations = prepared.faces.iter()
        .map(|face| face.transfer.quadrature_evaluations).sum::<usize>();
    let mut cache = BTreeMap::<(u64, u64), WallHeatTransfer>::new();
    let mut selected = BTreeSet::new();
    let mut reports = Vec::with_capacity(boundary_walls.len());
    for &boundary in boundary_walls {
        poll(gate)?;
        let cell = boundary.wall.fluid_cell;
        let face = boundary.wall.face as usize;
        let axis = face / 2;
        if cell >= count || !domain.is_fluid(cell) || domain.neighbor(cell, face).is_some() {
            return Err(invalid("exterior wall selection must be a boundary fluid cell face"));
        }
        if !selected.insert((cell, face)) || thermal.faces[face] != ThermalFace::Adiabatic {
            return Err(invalid("duplicate or already owned exterior thermal face; delegate it as Adiabatic"));
        }
        let wall_velocity = match mechanical_faces[face] {
            FvBoundary::Wall { velocity } => velocity,
            _ => return Err(invalid("thermal wall law requires a mechanical no-slip Wall, not an opening or symmetry plane")),
        };
        if !wall_velocity.iter().all(|v| v.is_finite()) || wall_velocity[axis] != 0.0
            || flow.field.outward(domain, cell, face) != 0.0
        {
            return Err(invalid("thermal wall must have finite tangential motion and zero normal through-flow"));
        }
        let contact = boundary.contact_m2_k_w;
        if !contact.is_finite() || contact < 0.0 {
            return Err(invalid("physical wall contact must be finite and nonnegative"));
        }
        let (film, reservoir, flux) = match boundary.condition {
            ThermalFace::Temperature(t) => (0.0, Some(t), None),
            ThermalFace::Convective { h, ambient } => {
                if !h.is_finite() || h <= 0.0 || !h.recip().is_finite() || h.recip() == 0.0 {
                    return Err(invalid("external film coefficient must be finite, positive and representable"));
                }
                (h.recip(), Some(ambient), None)
            }
            ThermalFace::HeatFlux(q) if q.is_finite() => (0.0, None, Some(q)),
            _ => return Err(invalid("use Temperature, Convective or finite inward HeatFlux at an exterior wall")),
        };
        if reservoir.is_some_and(|t| !t.is_finite() || t <= 0.0) {
            return Err(invalid("wall reservoir temperature must be finite and positive"));
        }
        let velocity = flow.velocity_m_s[cell];
        let a = (velocity[(axis + 1) % 3] - wall_velocity[(axis + 1) % 3]).abs();
        let b = (velocity[(axis + 2) % 3] - wall_velocity[(axis + 2) % 3]).abs();
        let scale = a.max(b);
        let speed = if scale == 0.0 { 0.0 } else {
            scale * ((a / scale).powi(2) + (b / scale).powi(2)).sqrt()
        };
        let distance = 0.5 * domain.widths(cell)[axis];
        let key = (speed.to_bits(), distance.to_bits());
        let transfer = if let Some(&cached) = cache.get(&key) { cached } else {
            let value = equilibrium_wall_heat_transfer(speed, distance, fluid, config, gate)?;
            evaluations = evaluations.checked_add(value.quadrature_evaluations)
                .ok_or_else(|| invalid("wall quadrature work count overflow"))?;
            cache.insert(key, value);
            value
        };
        let exterior = contact + film;
        let resistance = transfer.resistance_m2_k_w + exterior;
        let area = domain.face_area(cell, axis);
        let conductance = area / resistance;
        if !exterior.is_finite() || !resistance.is_finite() || resistance <= 0.0
            || !area.is_finite() || area <= 0.0 || !conductance.is_finite() || conductance <= 0.0
        {
            return Err(invalid("exterior wall transport is not representable"));
        }
        if let Some(t) = reservoir {
            if !(conductance * t).is_finite() {
                return Err(invalid("wall reservoir coefficient overflow"));
            }
            next.cell_sinks.push(CellSink {
                cell, conductance_w_k: conductance, temperature_k: t,
            });
        }
        if let Some(q) = flux {
            if next.power_w.is_empty() { next.power_w = vec![0.0; count]; }
            let watts = area * q;
            let total = next.power_w[cell] + watts;
            if !watts.is_finite() || (q != 0.0 && watts == 0.0) || !total.is_finite() {
                return Err(invalid("imposed wall heat is not representable"));
            }
            next.power_w[cell] = total;
        }
        reports.push(BoundaryWallHeatReport {
            boundary, area_m2: area, transfer,
            exterior_resistance_m2_k_w: exterior,
            conductance_w_k: if reservoir.is_some() { conductance } else { 0.0 },
        });
    }
    poll(gate)?;
    Ok(PreparedWallHeatBoundaries {
        thermal: next, internal_faces: prepared.faces, boundary_faces: reports,
        quadrature_evaluations: evaluations,
    })
}

/// Solve the original energy equations with both internal and exterior wall
/// transport, then recover the exterior surface temperatures and signed heat.
#[allow(clippy::too_many_arguments)]
pub fn solve_energy_with_wall_heat_boundaries(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FvFlow,
    thermal: &ThermalSetup,
    mechanical_faces: &[FvBoundary; 6],
    internal_walls: &[WallHeatFace],
    boundary_walls: &[BoundaryWallHeat],
    wall_config: WallHeatConfig,
    energy_config: &EnergyConfig,
    gate: &CancelGate,
) -> Result<WallHeatBoundarySolution, ChtError> {
    let prepared = prepare_wall_heat_boundaries(domain, fluid, flow, thermal,
        mechanical_faces, internal_walls, boundary_walls, wall_config, gate)?;
    let energy = solve_energy(domain, fluid, solids, &flow.field,
        &prepared.thermal, energy_config, gate)?;
    let boundary_faces = prepared.boundary_exchanges(&energy.temperature, gate)?;
    Ok(WallHeatBoundarySolution { energy, internal_faces: prepared.internal_faces, boundary_faces })
}

#[cfg(test)]
mod tests;
