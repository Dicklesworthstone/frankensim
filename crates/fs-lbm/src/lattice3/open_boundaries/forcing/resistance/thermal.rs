//! Native D3Q7 heat transport with the actual implicit porous velocity epochs.
//!
//! Heat uses the existing collision/history, reservoir, finite-wall and tensor
//! diffusion owners. Both accepted fields and their complete receipts publish
//! together; no force-only step can detach the resistance history.

use super::{DarcyForchheimer3, ResistanceOpenEndpointExchange3, ResistanceOpenGrid3};
use crate::d3q19::ResistanceMomentumExchange3;
use crate::lattice3::open_boundaries::thermal::OpenThermalError3;
use crate::lattice3::{EquilibriumLattice3, GridError3, macro_cell};
use crate::thermal::open::{
    ThermalOpenBalance3, ThermalOpenBoundaries3, ThermalOpenError3, ThermalOpenField3,
    ThermalOpenLimits3,
};
use crate::thermal::tensor::DiffusivityField;
use fs_exec::Cx;

mod endpoint;
pub use endpoint::ResistanceOpenEndpointThermalStep3;

/// Independently measured porous momentum/work and fixed-capacity heat.
/// Mechanical drag, skeleton and wall work retain their original receipts;
/// no constitutive heat partition is inferred from an isothermal momentum step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResistanceOpenThermalStep3 {
    pub momentum: ResistanceMomentumExchange3,
    pub heat: ThermalOpenBalance3,
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> ResistanceOpenGrid3<L, Q> {
    /// Bind original D3Q7 storage to this owner's topology, accepted time and
    /// physical velocity, including each cell's retained implicit half force.
    ///
    /// Heat remains a passive scalar with common fixed reference heat capacity.
    /// The existing boundary ownership and kinetic/temperature/storage budgets
    /// apply. Moving finite walls still require solid heat advection and refuse;
    /// stationary finite walls, immersed storage and their original conduction
    /// operators remain available through the returned thermal field.
    pub fn prepare_thermal_field(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
        thermal_relaxation: f64,
        boundaries: ThermalOpenBoundaries3,
        limits: ThermalOpenLimits3,
    ) -> Result<ThermalOpenField3, OpenThermalError3<L::Error>> {
        self.flow.prepare_thermal_field(
            cx, temperature, thermal_relaxation, boundaries, limits,
        )
    }

    /// Advance porous momentum and scalar heat as one native transaction.
    ///
    /// Opening heat collision, advection and the first source half-kick use the
    /// accepted physical velocity. The closing source half-kick uses the staged
    /// physical velocity with the force resolved by the original implicit law.
    /// Reading raw population momentum or holding the opening force at that
    /// endpoint would produce a different velocity and is not used.
    ///
    /// `heat_per_step` is a signed theta increment in x-fast cell order; solid
    /// entries must be zero. Optional thermal tau has the original per-cell
    /// meaning. Reservoir and finite-wall heat are measured from actual thermal
    /// populations. The complete porous force/work/wall/periodic receipt remains
    /// independent. No drag, skeleton or moving-wall work is silently deposited.
    ///
    /// Porous laws, geometry and storage capacity remain held. The caller may
    /// supply a momentum collision field; no temperature feedback is inferred.
    /// This shares the existing partitioned source ordering, without a new
    /// temporal-order, compressible enthalpy or total-fluid-energy claim.
    pub fn step_serial_with_thermal(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
    ) -> Result<ResistanceOpenThermalStep3, OpenThermalError3<L::Error>> {
        self.step_with_thermal_poll(
            thermal, collisions, heat_per_step, thermal_relaxation, None,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Same transaction with the existing uniform or per-cell SPD diffusivity
    /// tensor. Off-diagonal thermal transport, finite-wall conduction, reference
    /// storage and endpoint source corrections retain their original meanings.
    /// Tensor entries use lattice diffusivity on the fixed Cartesian axes;
    /// the porous material axes do not rotate the independently supplied tensor.
    pub fn step_serial_with_thermal_tensor_diffusion(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        diffusion: DiffusivityField<'_, 3>,
    ) -> Result<ResistanceOpenThermalStep3, OpenThermalError3<L::Error>> {
        self.step_with_thermal_poll(
            thermal, collisions, heat_per_step, None, Some(diffusion),
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    fn step_with_thermal_poll(
        &mut self,
        thermal: &mut ThermalOpenField3,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        diffusion: Option<DiffusivityField<'_, 3>>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<ResistanceOpenThermalStep3, OpenThermalError3<L::Error>> {
        let (momentum, _, step, heat) = self.stage_with_thermal_poll(
            thermal, collisions, heat_per_step, thermal_relaxation, diffusion, None, poll,
        )?;
        poll()?;
        self.publish_step(step, momentum);
        thermal.publish(step, heat);
        Ok(ResistanceOpenThermalStep3 { momentum, heat })
    }

    fn stage_with_thermal_poll(
        &mut self,
        thermal: &mut ThermalOpenField3,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        diffusion: Option<DiffusivityField<'_, 3>>,
        next_laws: Option<(&[DarcyForchheimer3], f64)>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(ResistanceMomentumExchange3, ResistanceOpenEndpointExchange3, usize, ThermalOpenBalance3), OpenThermalError3<L::Error>> {
        poll()?;
        if thermal.dims != self.flow.inner.grid.dims
            || thermal.solid != self.flow.inner.grid.solid
            || thermal.completed_steps() != self.flow.completed_steps()
        {
            return Err(ThermalOpenError3::Invalid(
                "porous and thermal geometry or accepted time differ",
            ).into());
        }
        self.flow.inner.admit_thermal_faces(thermal.boundaries())?;
        for cell in 0..thermal.velocities.len() {
            if cell.is_multiple_of(64) {
                poll()?;
            }
            thermal.velocities[cell] = if thermal.solid[cell] {
                [0.0; 3]
            } else {
                self.flow.macroscopic(cell)?.1
            };
        }
        let candidate = if let Some(diffusion) = diffusion {
            thermal.begin_step_tensor::<OpenThermalError3<L::Error>>(
                heat_per_step, diffusion, || poll().map_err(Into::into),
            )?
        } else {
            thermal.begin_step::<OpenThermalError3<L::Error>>(
                heat_per_step, thermal_relaxation, || poll().map_err(Into::into),
            )?
        };
        // The same porous collision/stream, implicit endpoints and independent
        // work/wall/periodic gates run once, stopping before physical publication.
        let (momentum, endpoint, step) =
            self.stage_with_laws::<fs_exec::TilePool>(None, collisions, next_laws, poll)?;
        for cell in 0..thermal.velocities.len() {
            if cell.is_multiple_of(64) {
                poll()?;
            }
            thermal.velocities[cell] = if thermal.solid[cell] {
                [0.0; 3]
            } else {
                macro_cell::<L, Q>(
                    &self.flow.inner.grid.next[cell],
                    self.flow.pending_forces[cell],
                    cell,
                    self.flow.inner.grid.budget.maximum_mach,
                )?.1
            };
        }
        let heat = thermal.finish_step::<OpenThermalError3<L::Error>>(
            heat_per_step, candidate, || poll().map_err(Into::into),
        )?;
        poll()?;
        Ok((momentum, endpoint, step, heat))
    }
}

#[cfg(test)]
mod tests;
