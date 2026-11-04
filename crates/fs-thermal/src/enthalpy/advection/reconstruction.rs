//! Limited axial enthalpy reconstruction on the original implicit transport.
//!
//! A cell's two one-sided axial derivatives are limited by minmod. The actual
//! signed flow selects its upstream cell and the corresponding face trace.
//! Distances are measured along the authored centreline, including bends:
//! `d_left=(L_previous+L_cell)/2`, `d_right=(L_cell+L_next)/2`.
//! Each trace lies between the centre and its adjacent sample; an extremum has
//! exactly zero slope. One identical reconstructed power debits and credits the
//! original native transfer. Sources, continuity, calorics and the implicit
//! energy residual remain owned by the existing timeline.
//!
//! This is an endpoint spatial reconstruction, not an explicit time integrator
//! or a blanket TVD/convergence certificate for a coupled nonlinear solve.
//! Smooth interior affine fields are reconstructed exactly (up to roundoff);
//! absent stencils, junctions and extrema retain the original upwind trace.
//! Nonuniform-grid limiters require actual distances rather than uniform-grid
//! jump ratios; see X. Zeng, "A general approach to enhance slope limiters on
//! non-uniform rectilinear grids", <https://arxiv.org/abs/1301.0967>.

use super::{Advected, finite, invalid};
use crate::enthalpy::backward_euler::TransportStepCoefficients;
use crate::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionSimulation, MaterialId, ThermalStorageBasis,
};
use std::collections::BTreeMap;

/// Hard ceiling on authored centre stencils in one reconstruction plan.
pub const MAX_AXIAL_STENCILS: usize = 65_536;

/// Three consecutive original control volumes in declared path order.
/// `previous` and `next` describe geometry, independently of current flow sign.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AxialCellStencil {
    pub previous: CellId,
    pub cell: CellId,
    pub next: CellId,
}

#[derive(Debug, Clone, Copy)]
struct PreparedStencil {
    previous: CellId,
    next: CellId,
    left_fraction: f64,
    right_fraction: f64,
}

#[derive(Debug, Clone, Copy)]
struct BoundCell {
    geometry: CellGeometry,
    mass_kg: f64,
    material_id: MaterialId,
}

/// Immutable numerical geometry for selected interior pipe cells.
///
/// The plan pins every referenced cell's mass, geometry and material ID. Every
/// candidate must retain those bindings. Original source/caloric and mass-flow
/// admission is still performed by `Advected::endpoint_transport` at that same
/// candidate. No junction distance, flow, density or material law is inferred.
#[derive(Debug, Clone)]
pub struct LimitedEnthalpyReconstruction {
    stencils: BTreeMap<CellId, PreparedStencil>,
    cells: BTreeMap<CellId, BoundCell>,
}

impl LimitedEnthalpyReconstruction {
    /// Prepare exact contiguous centreline geometry and a bounded set of unique
    /// interior centres. Empty plans preserve the existing upwind path exactly.
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        stencils: &[AxialCellStencil],
        maximum_stencils: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if maximum_stencils > MAX_AXIAL_STENCILS || stencils.len() > maximum_stencils {
            return Err(DepositionError::Budget("axial enthalpy stencils"));
        }
        let mut result = Self {
            stencils: BTreeMap::new(),
            cells: BTreeMap::new(),
        };
        for stencil in stencils {
            checkpoint()?;
            if stencil.previous == stencil.cell
                || stencil.next == stencil.cell
                || stencil.previous == stencil.next
                || result.stencils.contains_key(&stencil.cell)
            {
                return Err(invalid(
                    "axial stencils require distinct cells and unique centres",
                ));
            }
            let cells = [stencil.previous, stencil.cell, stencil.next].map(|id| {
                thermal
                    .cells()
                    .get(&id)
                    .ok_or(invalid("inactive axial stencil cell"))
            });
            let [previous, centre, next] = cells;
            let cells = [previous?, centre?, next?];
            let geometry = cells.map(|cell| cell.geometry);
            let mut material = None;
            for cell in cells {
                if thermal.material_storage_basis(cell.material_id)?
                    != ThermalStorageBasis::Enthalpy
                    || material.is_some_and(|expected| expected != cell.material_id)
                {
                    return Err(invalid(
                        "axial reconstruction requires one original enthalpy material",
                    ));
                }
                material = Some(cell.material_id);
                result.cells.insert(
                    cell.id,
                    BoundCell {
                        geometry: cell.geometry,
                        mass_kg: cell.mass_kg,
                        material_id: cell.material_id,
                    },
                );
            }
            if geometry[0].end_m != geometry[1].start_m || geometry[1].end_m != geometry[2].start_m
            {
                return Err(invalid(
                    "axial stencil centreline sections are not contiguous",
                ));
            }
            let previous = length(geometry[0])?;
            let centre = length(geometry[1])?;
            let next = length(geometry[2])?;
            let half = 0.5 * centre;
            let left = 0.5 * previous + half;
            let right = half + 0.5 * next;
            if !positive(half) || !positive(left) || !positive(right) {
                return Err(invalid(
                    "axial reconstruction distances are unrepresentable",
                ));
            }
            result.stencils.insert(
                stencil.cell,
                PreparedStencil {
                    previous: stencil.previous,
                    next: stencil.next,
                    left_fraction: half / left,
                    right_fraction: half / right,
                },
            );
        }
        checkpoint()?;
        Ok(result)
    }

    pub fn len(&self) -> usize {
        self.stencils.len()
    }

    pub fn is_empty(&self) -> bool {
        self.stencils.is_empty()
    }

    fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        flow: &Advected<()>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        if flow.liquid_carrier.is_some() {
            return Err(invalid(
                "bulk-enthalpy reconstruction has no liquid-carrier phase closure",
            ));
        }
        for (&id, bound) in &self.cells {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(invalid("axial stencil lost an original cell"))?;
            if cell.geometry != bound.geometry
                || cell.mass_kg != bound.mass_kg
                || cell.material_id != bound.material_id
                || thermal.material_storage_basis(cell.material_id)?
                    != ThermalStorageBasis::Enthalpy
            {
                return Err(invalid(
                    "axial reconstruction changed original storage or geometry",
                ));
            }
        }
        if flow.stationary_species.is_some()
            || !flow.selective_outlets.is_empty()
            || flow.entropy.is_some()
        {
            return Err(invalid(
                "limited enthalpy reconstruction has no species or entropy-mixing closure",
            ));
        }
        let mut neighbors: BTreeMap<_, u8> = self.stencils.keys().map(|&id| (id, 0)).collect();
        for edge in &flow.flows {
            checkpoint()?;
            for (centre, other) in [(edge.a, edge.b), (edge.b, edge.a)] {
                if let Some(stencil) = self.stencils.get(&centre) {
                    let bit = if other == stencil.previous {
                        1
                    } else if other == stencil.next {
                        2
                    } else {
                        return Err(invalid(
                            "axial reconstruction cannot infer a junction slope",
                        ));
                    };
                    *neighbors
                        .get_mut(&centre)
                        .ok_or(DepositionError::Arithmetic)? |= bit;
                }
            }
        }
        if neighbors.values().any(|&mask| mask != 3) {
            return Err(invalid(
                "axial reconstruction requires both original neighbor connections",
            ));
        }
        for boundary in &flow.boundaries {
            checkpoint()?;
            if self.stencils.contains_key(&boundary.identity().1) {
                return Err(invalid(
                    "axial reconstruction cannot infer an external-stream gradient",
                ));
            }
        }
        Ok(())
    }

    fn face_increment(
        &self,
        thermal: &DepositionSimulation<'_>,
        up: CellId,
        down: CellId,
    ) -> Result<f64, DepositionError> {
        let Some(stencil) = self.stencils.get(&up) else {
            return Ok(0.0);
        };
        let centre = thermal.cells()[&up].specific_h_j_kg;
        let left = finite(centre - thermal.cells()[&stencil.previous].specific_h_j_kg)?;
        let right = finite(thermal.cells()[&stencil.next].specific_h_j_kg - centre)?;
        // Limit physical derivatives times the SAME half-width. Fractions lie
        // in [0,1], avoiding overflow from dividing a finite jump by tiny length.
        let slope_half = minmod(left * stencil.left_fraction, right * stencil.right_fraction);
        let increment = if down == stencil.next {
            slope_half
        } else if down == stencil.previous {
            -slope_half
        } else {
            return Err(invalid(
                "reconstructed face is outside its original axial stencil",
            ));
        };
        if increment == 0.0 {
            return Ok(0.0);
        }
        let face_h = finite(centre + increment)?;
        let law = thermal
            .materials
            .get(&thermal.cells()[&up].material_id)
            .ok_or(invalid("unregistered axial caloric law"))?;
        if !positive(law.temperature_k(face_h)?) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(increment)
    }
}

impl Advected<()> {
    /// Evaluate limited internal face enthalpies on the SAME implicit residual.
    ///
    /// This first runs the original endpoint query, preserving source support,
    /// local mass balance, component-datum subtraction and boundary heat. It then
    /// changes only the selected internal donor powers, equally and oppositely.
    /// Flow reversal uses the newly upstream cell's actual opposite face.
    ///
    /// No observer/CFL, entropy-mixing account, physical step, nonlinear solver
    /// or accepted-state clamp is supplied. Feed the result to the original
    /// `BackwardEulerStep::evaluate_with_transport` and its acceptance gates.
    /// Spatial and temporal refinement remain independent obligations.
    pub fn endpoint_transport_with_reconstruction(
        &self,
        thermal: &DepositionSimulation<'_>,
        reconstruction: &LimitedEnthalpyReconstruction,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<TransportStepCoefficients<()>, DepositionError> {
        checkpoint()?;
        if reconstruction.is_empty() {
            return self.endpoint_transport(thermal, checkpoint);
        }
        reconstruction.admit(thermal, self, &mut checkpoint)?;
        let mut endpoint = self.endpoint_transport(thermal, &mut checkpoint)?;
        let mut transfers = endpoint.transfers.iter_mut();
        for edge in &self.flows {
            checkpoint()?;
            let Some((up, down, rate)) = edge.directed() else {
                continue;
            };
            let transfer = transfers.next().ok_or(DepositionError::Arithmetic)?;
            if transfer.a != up || transfer.b != down {
                return Err(invalid("native upwind transfer order changed"));
            }
            let increment = reconstruction.face_increment(thermal, up, down)?;
            if increment != 0.0 {
                transfer.power_w = finite(transfer.power_w + finite(rate * increment)?)?;
            }
        }
        if transfers.next().is_some() {
            return Err(invalid("native upwind transfer count changed"));
        }
        checkpoint()?;
        Ok(endpoint)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn length(geometry: CellGeometry) -> Result<f64, DepositionError> {
    if geometry
        .start_m
        .iter()
        .chain(&geometry.end_m)
        .any(|v| !v.is_finite())
    {
        return Err(invalid(
            "axial reconstruction needs finite centreline geometry",
        ));
    }
    let mut delta = [0.0; 3];
    for axis in 0..3 {
        delta[axis] = finite(geometry.end_m[axis] - geometry.start_m[axis])?;
    }
    let scale = delta.iter().fold(0.0_f64, |a, value| a.max(value.abs()));
    if scale == 0.0 {
        return Err(invalid(
            "axial reconstruction needs positive centreline lengths",
        ));
    }
    let result = scale * fs_math::det::sqrt(delta.iter().map(|v| (v / scale).powi(2)).sum());
    if !positive(result) {
        return Err(invalid("axial reconstruction length is unrepresentable"));
    }
    Ok(result)
}

fn minmod(left: f64, right: f64) -> f64 {
    if (left > 0.0 && right > 0.0) || (left < 0.0 && right < 0.0) {
        left.abs().min(right.abs()).copysign(left)
    } else {
        0.0
    }
}

#[cfg(test)]
#[path = "reconstruction/tests.rs"]
mod tests;
