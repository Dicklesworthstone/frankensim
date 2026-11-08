//! Limited axial enthalpy reconstruction on the original thermal transport.
//!
//! A cell's two one-sided axial derivatives are limited by minmod. The actual
//! signed flow selects its upstream cell and the corresponding face trace.
//! Distances are measured along the authored centreline, including bends:
//! `d_left=(L_previous+L_cell)/2`, `d_right=(L_cell+L_next)/2`.
//! Each trace lies between the centre and its adjacent sample; an extremum has
//! exactly zero slope. One identical reconstructed power debits and credits the
//! original native transfer. Sources, continuity, calorics and the thermal
//! energy residual remain owned by the existing timeline.
//!
//! [`ReconstructedAdvected`] applies the same traces on the explicit timeline
//! with a combined conduction/advection stability bound. The endpoint query
//! supplies the original implicit residual. Neither is a blanket TVD/convergence
//! certificate for a coupled nonlinear solve.
//! Without the optional thermal entropy limiter, smooth interior affine fields
//! are reconstructed exactly (up to roundoff);
//! absent stencils, junctions and extrema retain the original upwind trace.
//! A stationary-skeleton liquid carrier reconstructs h_liquid=max(h,h_liquidus),
//! preserving liquid superheat without transporting the solid's latent fraction.
//! A dry stencil neighbor supplies no liquid gradient and retains donor upwind.
//! The bulk thermal entropy opt-in may reduce the same face correction using
//! the original caloric entropy law, while preserving its paired heat transfer.
//! Nonuniform-grid limiters require actual distances rather than uniform-grid
//! jump ratios; see X. Zeng, "A general approach to enhance slope limiters on
//! non-uniform rectilinear grids", <https://arxiv.org/abs/1301.0967>.

use super::entropy::{AdvectionEntropyRates, EndpointAdvectionEntropy, Sum};
use super::{Advected, finite, invalid};
use crate::enthalpy::backward_euler::TransportStepCoefficients;
use crate::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionSimulation, MaterialId, ThermalStorageBasis,
};
use std::collections::BTreeMap;

#[path = "reconstruction/entropy.rs"]
mod entropy;
#[path = "reconstruction/observer.rs"]
mod observer;
pub use observer::ReconstructedAdvected;

/// Original endpoint temperature and dT/dh, supplied by the native implicit
/// caloric owner. Explicit flux queries do not need a derivative callback.
pub(in crate::enthalpy) type ReconstructionCalorics<'a> =
    dyn Fn(CellId) -> Result<(f64, f64), DepositionError> + 'a;

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

struct FaceIncrement {
    value: f64,
    // Derivatives of the correction in carried-enthalpy coordinates. Ordinary
    // minmod uses donor/neighbor; the entropy bound uses donor/receiver.
    // Each column's own carrier slope converts this to stored mixture h.
    partials: Option<[(CellId, f64); 2]>,
    // Set by the original entropy query only for unequal-temperature faces.
    // Entropy-neutral latent traces retain their original flux arithmetic.
    guard_power_roundoff: bool,
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
    thermal_entropy_limit: bool,
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
            thermal_entropy_limit: false,
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

    /// Whether an original axial stencil requires this undirected connection.
    ///
    /// A live flow owner can retain these genuine geometric connections at
    /// zero mass rate, so stopping or reversing flow does not erase a required
    /// neighbor. Unrelated dormant transverse faces need not become axial
    /// junctions. This query supplies no flux and does not authorize a new
    /// geometric connection; the native flow and admission still own those.
    pub fn requires_neighbor_connection(&self, a: CellId, b: CellId) -> bool {
        [(a, b), (b, a)].into_iter().any(|(centre, other)| {
            self.stencils
                .get(&centre)
                .is_some_and(|stencil| other == stencil.previous || other == stencil.next)
        })
    }

    /// Limit each reconstructed bulk face by its original thermal entropy.
    ///
    /// The existing minmod trace is retained whenever its directed face
    /// production is nonnegative within finite-increment roundoff. Otherwise
    /// only its conservative face correction is reduced toward the original
    /// donor flux. Source entropy differences, actual temperatures and analytic
    /// caloric tangents supply the active bound; no heat capacity is averaged.
    /// Equal-temperature latent plateaus retain their original reconstruction.
    ///
    /// This opt-in requires ordinary same-material bulk transport and original
    /// caloric entropy support when a face correction is present. Stationary
    /// mixtures, selective streams and stationary-skeleton liquid carriage need
    /// separate closures. It is a local semidiscrete condition, not a guarantee
    /// about a finite time step, coupled work, or a general maximum principle.
    pub fn with_thermal_entropy_limit(mut self) -> Self {
        self.thermal_entropy_limit = true;
        self
    }

    pub(in crate::enthalpy) fn has_thermal_entropy_limit(&self) -> bool {
        self.thermal_entropy_limit
    }

    fn admit<O>(
        &self,
        thermal: &DepositionSimulation<'_>,
        flow: &Advected<O>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.admit_cells(thermal, checkpoint)?;
        self.admit_flow(flow, checkpoint)
    }

    fn admit_cells(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
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
        Ok(())
    }

    fn admit_flow<O>(
        &self,
        flow: &Advected<O>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.admit_flow_model(flow)?;
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

    fn admit_flow_model<O>(&self, flow: &Advected<O>) -> Result<(), DepositionError> {
        if flow.stationary_species.is_some()
            || !flow.selective_outlets.is_empty()
            || flow.entropy.is_some()
        {
            return Err(invalid(
                "limited enthalpy reconstruction has no species or entropy-mixing closure",
            ));
        }
        if self.thermal_entropy_limit && flow.liquid_carrier.is_some() {
            return Err(invalid(
                "thermal entropy face limiting requires ordinary bulk carriage",
            ));
        }
        Ok(())
    }

    fn face_increment<O>(
        &self,
        thermal: &DepositionSimulation<'_>,
        flow: &Advected<O>,
        up: CellId,
        down: CellId,
    ) -> Result<FaceIncrement, DepositionError> {
        let Some(stencil) = self.stencils.get(&up) else {
            return Ok(FaceIncrement {
                value: 0.0,
                partials: None,
                guard_power_roundoff: false,
            });
        };
        let centre = flow.carried_h(&thermal.cells()[&up])?;
        let previous = &thermal.cells()[&stencil.previous];
        let next = &thermal.cells()[&stencil.next];
        if let Some(carrier) = flow.liquid_carrier() {
            let solidus = carrier.plateau_enthalpy_bounds_j_kg()[0];
            if previous.specific_h_j_kg <= solidus || next.specific_h_j_kg <= solidus {
                // The native endpoint has already refused nonzero flow through
                // a dry cell. A zero-flow geometric neighbor has no moving
                // liquid state from which to infer a slope.
                return Ok(FaceIncrement {
                    value: 0.0,
                    partials: None,
                    guard_power_roundoff: false,
                });
            }
        }
        let left = finite(centre - flow.carried_h(previous)?)?;
        let right = finite(flow.carried_h(next)? - centre)?;
        // Limit physical derivatives times the SAME half-width. Fractions lie
        // in [0,1], avoiding overflow from dividing a finite jump by tiny length.
        let limited_left = left * stencil.left_fraction;
        let limited_right = right * stencil.right_fraction;
        let slope_half = minmod(limited_left, limited_right);
        let direction = if down == stencil.next {
            1.0
        } else if down == stencil.previous {
            -1.0
        } else {
            return Err(invalid(
                "reconstructed face is outside its original axial stencil",
            ));
        };
        let increment = direction * slope_half;
        if increment == 0.0 {
            return Ok(FaceIncrement {
                value: 0.0,
                partials: None,
                guard_power_roundoff: false,
            });
        }
        let face_h = finite(centre + increment)?;
        let law = thermal
            .materials
            .get(&thermal.cells()[&up].material_id)
            .ok_or(invalid("unregistered axial caloric law"))?;
        if !positive(law.temperature_k(face_h)?) {
            return Err(DepositionError::MaterialDomain);
        }
        // At equal nonzero slopes choose the left active derivative. At a
        // sign switch the zero branch above supplies a generalized derivative.
        // Minmod is piecewise linear, not classically differentiable there.
        let (neighbor, partial) = if limited_left.abs() <= limited_right.abs() {
            (stencil.previous, -direction * stencil.left_fraction)
        } else {
            (stencil.next, direction * stencil.right_fraction)
        };
        Ok(FaceIncrement {
            value: increment,
            partials: Some([(neighbor, partial), (up, -partial)]),
            guard_power_roundoff: false,
        })
    }
}

impl Advected<()> {
    /// Evaluate limited internal face enthalpies on the SAME implicit residual.
    ///
    /// This first runs the original endpoint query, preserving source support,
    /// local mass balance, component-datum subtraction and boundary heat. It then
    /// changes only the selected internal donor powers, equally and oppositely.
    /// Flow reversal uses the newly upstream cell's actual opposite face.
    /// With a liquid carrier, traces and slopes use that original liquid h;
    /// mixture h remains the storage coordinate. Dry stencil samples retain
    /// first-order upwind, and nonzero flow through a solid cell still refuses.
    ///
    /// No observer/CFL, entropy-mixing account, physical step, nonlinear solver
    /// or accepted-state clamp is supplied by this query. Feed the result to
    /// the original `BackwardEulerStep::evaluate_with_transport`, or attach
    /// this plan with `ImplicitConductionGraph::with_implicit_advection_reconstruction`
    /// to use the native nonlinear solve and its acceptance gates.
    /// Spatial and temporal refinement remain independent obligations.
    pub fn endpoint_transport_with_reconstruction(
        &self,
        thermal: &DepositionSimulation<'_>,
        reconstruction: &LimitedEnthalpyReconstruction,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<TransportStepCoefficients<()>, DepositionError> {
        self.reconstructed_endpoint(thermal, reconstruction, false, checkpoint)
            .map(|(transport, _)| transport)
    }

    /// The same conservative face flux with original advective entropy.
    /// Nonempty reconstruction requires `with_thermal_entropy_limit`. Boundary
    /// streams remain native upwind; internal face corrections change only
    /// spatial mixing. Actual transport arithmetic retains its separate balance
    /// residue, and no thermal or observer history is advanced by this query.
    pub fn endpoint_transport_with_reconstruction_and_entropy(
        &self,
        thermal: &DepositionSimulation<'_>,
        reconstruction: &LimitedEnthalpyReconstruction,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(TransportStepCoefficients<()>, AdvectionEntropyRates), DepositionError> {
        let (transport, _, entropy) = self.reconstructed_endpoint_with_calorics(
            thermal,
            reconstruction,
            false,
            true,
            None,
            checkpoint,
        )?;
        Ok((transport, entropy.ok_or(DepositionError::Arithmetic)?.rates))
    }

    /// The native implicit owner can request the active face-correction
    /// tangent in the same pass as the physical flux. Entries differentiate
    /// heat into a row with respect to a column's specific enthalpy [kg/s].
    /// The original upwind/component-reference tangent remains separate.
    pub(in crate::enthalpy) fn reconstructed_endpoint(
        &self,
        thermal: &DepositionSimulation<'_>,
        reconstruction: &LimitedEnthalpyReconstruction,
        differentiate: bool,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(TransportStepCoefficients<()>, Vec<(CellId, CellId, f64)>), DepositionError> {
        self.reconstructed_endpoint_with_calorics(
            thermal,
            reconstruction,
            differentiate,
            false,
            None,
            checkpoint,
        )
        .map(|(transport, tangent, _)| (transport, tangent))
    }

    pub(in crate::enthalpy) fn reconstructed_endpoint_with_calorics(
        &self,
        thermal: &DepositionSimulation<'_>,
        reconstruction: &LimitedEnthalpyReconstruction,
        differentiate: bool,
        with_entropy: bool,
        calorics: Option<&ReconstructionCalorics<'_>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            TransportStepCoefficients<()>,
            Vec<(CellId, CellId, f64)>,
            Option<EndpointAdvectionEntropy>,
        ),
        DepositionError,
    > {
        checkpoint()?;
        if !reconstruction.is_empty() {
            reconstruction.admit(thermal, self, &mut checkpoint)?;
            if with_entropy && !reconstruction.thermal_entropy_limit {
                return Err(invalid(
                    "reconstructed entropy requires thermal entropy face limiting",
                ));
            }
        }
        let (mut endpoint, mut entropy) = if with_entropy {
            let (transport, entropy) =
                self.endpoint_transport_entropy_state(thermal, &mut checkpoint)?;
            (transport, Some(entropy))
        } else {
            (self.endpoint_transport(thermal, &mut checkpoint)?, None)
        };
        if reconstruction.is_empty() {
            return Ok((endpoint, Vec::new(), entropy));
        }
        let mut tangent = Vec::new();
        let mut entropy_correction = Sum::default();
        if differentiate {
            let capacity = self
                .flows
                .len()
                .checked_mul(4)
                .ok_or(DepositionError::Budget("axial enthalpy tangent"))?;
            tangent
                .try_reserve(capacity)
                .map_err(|_| DepositionError::Budget("axial enthalpy tangent"))?;
        }
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
            let mut increment = reconstruction.face_increment(thermal, self, up, down)?;
            if reconstruction.thermal_entropy_limit && increment.value != 0.0 {
                increment =
                    entropy::limit_face(thermal, up, down, increment, differentiate, calorics)?;
            }
            if increment.value != 0.0 {
                let old_power = transfer.power_w;
                let correction_power = finite(rate * increment.value)?;
                transfer.power_w = if increment.guard_power_roundoff {
                    entropy::apply_correction(old_power, correction_power)?
                } else {
                    finite(old_power + correction_power)?
                };
                if increment.guard_power_roundoff && transfer.power_w == old_power {
                    // A correction below the actual flux representation's
                    // resolution retains the donor flux and donor tangent.
                    increment.partials = None;
                }
                if let Some(entropy) = &entropy {
                    let ta = entropy.temperatures[&up];
                    let tb = entropy.temperatures[&down];
                    // Classify the actually applied paired correction, including
                    // its rounding against the original component-reference flux.
                    entropy_correction.add(finite(
                        (transfer.power_w - old_power) * ((ta - tb) / ta / tb),
                    )?)?;
                }
            }
            if let Some(partials) = increment.partials.filter(|_| differentiate) {
                // One power, two opposite physical accounts. A liquid trace
                // differentiates the moving liquid's h, not the solid's latent
                // inventory: donor and neighbor can have different slopes.
                for (col, partial) in partials {
                    if partial == 0.0 {
                        continue;
                    }
                    let derivative = finite(rate * partial)?;
                    if derivative == 0.0 {
                        return Err(DepositionError::Arithmetic);
                    }
                    let slope = self.implicit_carried_enthalpy_slope(&thermal.cells()[&col])?;
                    if slope != 0.0 {
                        let value = finite(derivative * slope)?;
                        if value == 0.0 {
                            return Err(DepositionError::Arithmetic);
                        }
                        tangent.extend([(up, col, -value), (down, col, value)]);
                    }
                }
            }
        }
        if transfers.next().is_some() {
            return Err(invalid("native upwind transfer count changed"));
        }
        if let Some(entropy) = &mut entropy {
            entropy.rates.upwind_mixing_w_k =
                finite(entropy.rates.upwind_mixing_w_k + entropy_correction.total()?)?;
            entropy.reconcile_transport(&endpoint, &mut checkpoint)?;
        }
        checkpoint()?;
        Ok((endpoint, tangent, entropy))
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
