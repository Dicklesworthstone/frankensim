//! Passive heat transport on native spatial and time-dependent force epochs.
//!
//! The force owner interprets momentum at both accepted times. Its original
//! scalar companion owns convection, molecular diffusion, wall storage and
//! reservoir heat; no body-work term is silently added to sensible heat.

use super::super::forcing::ForcedOpenBoundaryGrid3;
use super::{
    BoundarySpec3, DiffusivityField, EquilibriumLattice3, FacePlan, GridError3,
    OpenThermalDiffusion3, OpenThermalError3, OpenThermalStep3,
    ThermalOpenBoundaries3, ThermalOpenError3, ThermalOpenField3, ThermalOpenLimits3, macro_cell,
};
use fs_exec::Cx;

impl<L: EquilibriumLattice3<Q>, const Q: usize> ForcedOpenBoundaryGrid3<L, Q> {
    /// Bind D3Q7 heat to the original native topology, accepted time and local
    /// physical velocities. Each velocity uses its own accepted force density,
    /// including the native half-source correction. Thermal boundary ownership,
    /// reference capacity and numerical/resource limits are the same as on
    /// [`super::OpenBoundaryGrid3::prepare_thermal_field`].
    pub fn prepare_thermal_field(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
        thermal_relaxation: f64,
        boundaries: ThermalOpenBoundaries3,
        limits: ThermalOpenLimits3,
    ) -> Result<ThermalOpenField3, OpenThermalError3<L::Error>> {
        self.prepare_thermal_field_with_poll(
            temperature,
            thermal_relaxation,
            boundaries,
            limits,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    fn prepare_thermal_field_with_poll(
        &self,
        temperature: &[f64],
        thermal_relaxation: f64,
        boundaries: ThermalOpenBoundaries3,
        limits: ThermalOpenLimits3,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<ThermalOpenField3, OpenThermalError3<L::Error>> {
        self.inner.admit_thermal_faces(boundaries)?;
        ThermalOpenField3::try_new(
            self.inner.grid.dims,
            &self.inner.grid.solid,
            temperature,
            thermal_relaxation,
            boundaries,
            limits,
            self.completed_steps(),
            |cell| self.macroscopic(cell).map(|(_, u)| u).map_err(Into::into),
            || poll().map_err(Into::into),
        )
    }

    /// Hold the accepted spatial force while advancing native flow and heat.
    /// Heat advection uses the opening physical velocity; the final prescribed
    /// source half-kick uses the native endpoint physical velocity. Finite walls,
    /// reservoirs, flux history and both time counters publish together.
    ///
    /// This is fixed-reference-capacity temperature transport, not a total-energy
    /// or variable-density caloric model. Mechanical work and buoyancy feedback
    /// are not inferred from the force field.
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
        self.step_with_thermal_poll(
            thermal,
            None,
            collisions,
            heat_per_step,
            thermal_relaxation,
            None,
            origin,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Advance with an authored spatial force at the next accepted time.
    /// Collision retains F_n; native reservoir reconstruction and final thermal
    /// source kicks use F_(n+1). The body impulse remains the native trapezoidal
    /// force integral, distinct from measured boundary momentum and thermal heat.
    ///
    /// On success `forces_at_end` receives the old accepted force buffer for
    /// reuse. On refusal or cancellation it is unchanged, and accepted flow,
    /// force, heat, finite walls, histories and counters all remain unchanged.
    /// No source-coordinate reframe or extra impulse occurs during a step.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_thermal(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        forces_at_end: &mut Vec<[f64; 3]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        origin: [f64; 3],
    ) -> Result<OpenThermalStep3, OpenThermalError3<L::Error>> {
        self.step_with_thermal_poll(
            thermal,
            Some(forces_at_end),
            collisions,
            heat_per_step,
            thermal_relaxation,
            None,
            origin,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Hold spatial force with the original full-tensor D3Q7 diffusion kernel.
    /// Tensor entries are lattice diffusivities on the fixed Cartesian axes;
    /// per-cell fields include solid entries. Storage and even-moment relaxation
    /// retain the scalar carrier's original definitions.
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
        self.step_with_thermal_poll(
            thermal,
            None,
            collisions,
            heat_per_step,
            None,
            Some(diffusion),
            origin,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Change spatial force with full-tensor thermal diffusion and one common
    /// publication point. Force-buffer exchange and all other semantics match
    /// [`Self::advance_serial_with_thermal`].
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_thermal_tensor_diffusion(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        forces_at_end: &mut Vec<[f64; 3]>,
        heat_per_step: &[f64],
        diffusion: DiffusivityField<'_, 3>,
        origin: [f64; 3],
    ) -> Result<OpenThermalStep3, OpenThermalError3<L::Error>> {
        self.step_with_thermal_poll(
            thermal,
            Some(forces_at_end),
            collisions,
            heat_per_step,
            None,
            Some(diffusion),
            origin,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Advance proposed reservoirs, force epochs and heat in one transaction.
    ///
    /// None holds the accepted spatial force without copying. Some supplies
    /// its endpoint buffer, which receives the old accepted force only after
    /// full thermal acceptance, as in advance_serial_with_thermal. Scalar and
    /// tensor diffusion preserve their existing kernels, physical velocity
    /// epochs, wall histories and separate thermal/momentum accounts. A late
    /// heat or cancellation refusal also retains the previous boundary targets.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_thermal_and_open_boundaries(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        boundaries: BoundarySpec3,
        collisions: Option<&[L::Collision]>,
        forces_at_end: Option<&mut Vec<[f64; 3]>>,
        heat_per_step: &[f64],
        diffusion: OpenThermalDiffusion3<'_>,
        origin: [f64; 3],
    ) -> Result<OpenThermalStep3, OpenThermalError3<L::Error>> {
        let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let plan = self.inner.prepare_open_boundaries(boundaries, &mut poll)?;
        let (relaxation, tensor) = diffusion.parts();
        self.step_with_thermal_plan_poll(
            thermal, plan, forces_at_end, collisions, heat_per_step, relaxation,
            tensor, origin, &mut poll,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_thermal_poll(
        &mut self,
        thermal: &mut ThermalOpenField3,
        forces_at_end: Option<&mut Vec<[f64; 3]>>,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        diffusion: Option<DiffusivityField<'_, 3>>,
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<OpenThermalStep3, OpenThermalError3<L::Error>> {
        self.step_with_thermal_plan_poll(
            thermal, self.inner.plan, forces_at_end, collisions, heat_per_step,
            thermal_relaxation, diffusion, origin, poll,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_thermal_plan_poll(
        &mut self,
        thermal: &mut ThermalOpenField3,
        plan: FacePlan,
        forces_at_end: Option<&mut Vec<[f64; 3]>>,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        diffusion: Option<DiffusivityField<'_, 3>>,
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<OpenThermalStep3, OpenThermalError3<L::Error>> {
        poll()?;
        if thermal.dims != self.inner.grid.dims
            || thermal.solid != self.inner.grid.solid
            || thermal.completed_steps() != self.completed_steps()
        {
            return Err(ThermalOpenError3::Invalid(
                "native and thermal geometry or accepted time differ",
            )
            .into());
        }
        self.inner.admit_thermal_faces_with_plan(plan, thermal.boundaries())?;
        for cell in 0..self.inner.grid.f.len() {
            if cell % 64 == 0 {
                poll()?;
            }
            thermal.velocities[cell] = if self.inner.grid.solid[cell] {
                [0.0; 3]
            } else {
                self.macroscopic(cell)?.1
            };
        }
        let candidate = if let Some(diffusion) = diffusion {
            thermal.begin_step_tensor::<OpenThermalError3<L::Error>>(
                heat_per_step,
                diffusion,
                || poll().map_err(Into::into),
            )?
        } else {
            thermal.begin_step::<OpenThermalError3<L::Error>>(
                heat_per_step,
                thermal_relaxation,
                || poll().map_err(Into::into),
            )?
        };
        let endpoint = forces_at_end
            .as_deref()
            .map_or(self.forces.as_slice(), Vec::as_slice);
        let (momentum, step) = self.inner.stage_serial_with_forces(
            plan,
            collisions,
            Some(&self.forces),
            Some(endpoint),
            origin,
            poll,
        )?;
        for (cell, &force) in endpoint.iter().enumerate() {
            if cell % 64 == 0 {
                poll()?;
            }
            thermal.velocities[cell] = if self.inner.grid.solid[cell] {
                [0.0; 3]
            } else {
                macro_cell::<L, Q>(
                    &self.inner.grid.next[cell],
                    force,
                    cell,
                    self.inner.grid.budget.maximum_mach,
                )?
                .1
            };
        }
        let heat =
            thermal.finish_step::<OpenThermalError3<L::Error>>(heat_per_step, candidate, || {
                poll().map_err(Into::into)
            })?;
        poll()?;
        if let Some(forces_at_end) = forces_at_end {
            self.publish_with_force_swap(step, forces_at_end);
        } else {
            core::mem::swap(&mut self.inner.grid.f, &mut self.inner.grid.next);
            self.inner.grid.completed_steps = step;
            self.inner.publish_wall_exchange();
        }
        self.inner.plan = plan;
        thermal.publish(step, heat);
        Ok(OpenThermalStep3 { momentum, heat })
    }
}

#[cfg(test)]
mod tests;
