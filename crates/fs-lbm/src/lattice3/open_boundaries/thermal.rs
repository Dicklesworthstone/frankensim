//! Joint native momentum / D3Q7 thermal publication on open grids.
//!
//! Momentum retains its original collision, force correction, reservoirs and
//! Mach admission. The thermal carrier reads those physical velocities and
//! uses the existing scalar collision, flux history and symmetric source
//! kicks. Neither field, finite thermal wall, receipt nor time counter changes
//! on any refused or cancelled step. Scratch buffers may contain a prefix.

/// Temperature-dependent buoyancy with retained native force history.
pub mod buoyancy;

use super::{
    EquilibriumLattice3, Face3, FaceBoundary3, GridError3, OpenBoundaryGrid3,
    OpenBoundaryMomentumTransport3, macro_cell,
};
use crate::thermal::open::{
    ThermalOpenBalance3, ThermalOpenBoundaries3, ThermalOpenError3, ThermalOpenFace3,
    ThermalOpenField3, ThermalOpenLimits3,
};
use crate::thermal::tensor::DiffusivityField;
use fs_exec::Cx;

#[derive(Debug, Clone)]
pub enum OpenThermalError3<E> {
    Flow(GridError3<E>),
    Thermal(ThermalOpenError3),
}
impl<E> From<GridError3<E>> for OpenThermalError3<E> {
    fn from(error: GridError3<E>) -> Self {
        Self::Flow(error)
    }
}
impl<E> From<ThermalOpenError3> for OpenThermalError3<E> {
    fn from(error: ThermalOpenError3) -> Self {
        Self::Thermal(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for OpenThermalError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native open flow/thermal step: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for OpenThermalError3<E> {}

/// Independently measured native momentum and reference-capacity heat fluxes.
#[derive(Debug, Clone)]
pub struct OpenThermalStep3 {
    pub momentum: OpenBoundaryMomentumTransport3,
    pub heat: ThermalOpenBalance3,
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    /// Bind a D3Q7 temperature field to the current native topology and time.
    /// `temperature` is x-fast and covers all grid cells; solid entries carry
    /// no thermal storage. External finite walls receive independent state
    /// for every normal fluid link, in the same x-fast surface ordering.
    ///
    /// Every native velocity/pressure face requires explicit reservoir heat
    /// data, including outlets that could reverse. Periodic and fixed-geometry
    /// wall face kinds must agree, including tangentially moving walls. Moving
    /// faces admit stateless thermal laws; finite solid storage would require
    /// solid heat advection and refuses. The scalar |u_axis| < 1/4 kinetic
    /// envelope is additional to the original native Mach budget.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_thermal_field(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
        thermal_relaxation: f64,
        boundaries: ThermalOpenBoundaries3,
        limits: ThermalOpenLimits3,
    ) -> Result<ThermalOpenField3, OpenThermalError3<L::Error>> {
        self.admit_thermal_faces(boundaries)?;
        ThermalOpenField3::try_new(
            self.grid.dims,
            &self.grid.solid,
            temperature,
            thermal_relaxation,
            boundaries,
            limits,
            self.grid.completed_steps,
            |cell| {
                self.macroscopic(cell)
                    .map(|(_, velocity)| velocity)
                    .map_err(Into::into)
            },
            || {
                cx.checkpoint()
                    .map_err(|_| OpenThermalError3::Flow(GridError3::Cancelled))
            },
        )
    }

    /// Advance native flow and heat with one final publication point.
    ///
    /// `heat_per_step` is a frozen signed temperature increment over one step;
    /// solid grid entries must be zero. Optional thermal relaxation is one
    /// admitted scalar tau per grid cell. Sources use the original equilibrium
    /// half-kicks and flux-history rebinding, first at accepted native velocity
    /// and finally at the admitted native candidate velocity.
    ///
    /// This carrier transports a passive, fixed reference-capacity scalar.
    /// The held native body force and caller's collision laws are unchanged;
    /// no hidden buoyancy model, material law or compressible enthalpy equation
    /// is introduced. The full partitioned method has no new temporal-order or
    /// nonreflecting-outlet claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_thermal(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        origin: [f64; 3],
    ) -> Result<OpenThermalStep3, OpenThermalError3<L::Error>> {
        self.step_serial_with_thermal_impl(
            thermal,
            cx,
            collisions,
            heat_per_step,
            thermal_relaxation,
            None,
            origin,
        )
    }

    /// Advance native flow and heat with a full SPD thermal diffusivity tensor.
    ///
    /// Uniform or per-cell tensors use the original D3Q7 flux-relaxation block
    /// in the fixed Cartesian grid axes. Off-diagonal entries retain oblique
    /// heat transport. Supply `alpha_lattice = K/(rho*Cp) * dt/dx^2` for the
    /// carrier's common reference volumetric heat capacity on a cubic grid;
    /// variable storage and material interface resistance are not inferred.
    /// Per-cell fields cover every grid cell, including unused solid entries.
    ///
    /// The field is borrowed for this step only. The carrier's original tau
    /// still relaxes even kinetic moments. Sources, advective-flux correction,
    /// reservoir/wall heat, finite-wall storage and atomic publication use the
    /// same owner as [`Self::step_serial_with_thermal`]. A later scalar step
    /// retains its original coefficients. SPD alone guarantees neither
    /// population positivity nor coupled stability; refinement remains required.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_thermal_tensor_diffusion(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        diffusion: DiffusivityField<'_, 3>,
        origin: [f64; 3],
    ) -> Result<OpenThermalStep3, OpenThermalError3<L::Error>> {
        self.step_serial_with_thermal_impl(
            thermal,
            cx,
            collisions,
            heat_per_step,
            None,
            Some(diffusion),
            origin,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_serial_with_thermal_impl(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        diffusion: Option<DiffusivityField<'_, 3>>,
        origin: [f64; 3],
    ) -> Result<OpenThermalStep3, OpenThermalError3<L::Error>> {
        let poll = || -> Result<(), OpenThermalError3<L::Error>> {
            cx.checkpoint()
                .map_err(|_| OpenThermalError3::Flow(GridError3::Cancelled))
        };
        poll()?;
        if thermal.dims != self.grid.dims
            || thermal.solid != self.grid.solid
            || thermal.completed_steps() != self.grid.completed_steps
        {
            return Err(ThermalOpenError3::Invalid(
                "native and thermal geometry or accepted time differ",
            )
            .into());
        }
        self.admit_thermal_faces(thermal.boundaries())?;
        for cell in 0..self.grid.f.len() {
            if cell % 64 == 0 {
                poll()?;
            }
            thermal.velocities[cell] = if self.grid.solid[cell] {
                [0.0; 3]
            } else {
                self.macroscopic(cell)?.1
            };
        }
        let thermal_candidate = if let Some(diffusion) = diffusion {
            thermal.begin_step_tensor(heat_per_step, diffusion, poll)?
        } else {
            thermal.begin_step(heat_per_step, thermal_relaxation, poll)?
        };
        let (momentum, step) = self.stage_serial_with_plan(cx, self.plan, collisions, origin)?;
        for cell in 0..self.grid.next.len() {
            if cell % 64 == 0 {
                poll()?;
            }
            thermal.velocities[cell] = if self.grid.solid[cell] {
                [0.0; 3]
            } else {
                macro_cell::<L, Q>(
                    &self.grid.next[cell],
                    self.grid.force,
                    cell,
                    self.grid.budget.maximum_mach,
                )?
                .1
            };
        }
        let heat = thermal.finish_step(heat_per_step, thermal_candidate, poll)?;
        poll()?;
        core::mem::swap(&mut self.grid.f, &mut self.grid.next);
        self.grid.completed_steps = step;
        self.publish_wall_exchange();
        thermal.publish(step, heat);
        Ok(OpenThermalStep3 { momentum, heat })
    }

    fn admit_thermal_faces(
        &self,
        thermal: ThermalOpenBoundaries3,
    ) -> Result<(), OpenThermalError3<L::Error>> {
        for face in Face3::ALL {
            if matches!(self.plan.boundaries.face(face), FaceBoundary3::Wall { velocity } if velocity != [0.0; 3])
                && matches!(
                    thermal.faces[face as usize],
                    ThermalOpenFace3::Wall(crate::thermal::boundary::ThermalWall::FiniteSolid(_))
                )
            {
                return Err(ThermalOpenError3::Invalid(
                    "moving finite-wall storage requires solid heat advection",
                )
                .into());
            }
            match (
                self.plan.boundaries.face(face),
                thermal.faces[face as usize],
            ) {
                (FaceBoundary3::Periodic, ThermalOpenFace3::Periodic)
                | (FaceBoundary3::Wall { .. }, ThermalOpenFace3::Wall(_))
                | (
                    FaceBoundary3::Velocity { .. } | FaceBoundary3::Pressure { .. },
                    ThermalOpenFace3::Reservoir { .. },
                ) => {}
                _ => {
                    return Err(ThermalOpenError3::Invalid(
                        "native and thermal face ownership differ",
                    )
                    .into());
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
