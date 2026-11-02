//! Synchronous native Verlet stages around one shared endpoint waterline.
use super::{FiniteTankBodiesLoad, FiniteTankBody};
use crate::hydrostatic::{
    FiniteHydrostaticTank, HydrostaticError, HydrostaticStepControls, HydrostaticStepError, finite,
    poll,
};
use fs_material::fluid_transport::FluidTransport;
use fs_mbd::{Gravity, RigidBodyIntegrator, RigidBodyState, StepReceipt};

/// One atomically proposed batch of native steps and one finite-fluid account.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteTankBodiesStep {
    /// Complete native per-body steps in the original binding order.
    pub dynamics: Vec<StepReceipt>,
    /// One initial surface and all original buoyancy wrenches.
    pub load_before: FiniteTankBodiesLoad,
    /// One surface at ALL drifted poses and all final-half-kick wrenches.
    pub load_after: FiniteTankBodiesLoad,
    /// Sum of native mechanical energies plus one fluid potential before [J].
    pub total_energy_before_j: f64,
    /// The same total at the proposed complete endpoint [J].
    pub total_energy_after_j: f64,
    /// Signed numerical storage change; not fluid damping or heat [J].
    pub energy_defect_j: f64,
}

impl<F: FluidTransport + ?Sized> FiniteHydrostaticTank<'_, F> {
    /// Advance all bodies against the SAME liquid surface at each endpoint.
    ///
    /// The original native Verlet owner half-kicks and drifts each body once.
    /// One collective endpoint solve supplies their final half-kicks. This
    /// retains cross-body hydrostatic restoring forces and the full shared
    /// displacement-potential cross term without a separate motion integrator.
    ///
    /// The absolute energy gate applies to the whole batch: body energies plus
    /// the liquid potential exactly once. A cancellation, frame/geometry/force
    /// refusal, overflow or energy failure publishes no completed batch.
    /// Callers accept all returned states together, keep the same bindings,
    /// inventory and source between steps, and ensure nonintersecting motion.
    /// This adds no collision detection, fluid inertia, waves or hydrodynamic
    /// damping. Timestep refinement remains required for this quasistatic model.
    pub fn step_bodies(
        &self,
        bodies: &[FiniteTankBody<'_>],
        states: &[RigidBodyState],
        controls: HydrostaticStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteTankBodiesStep, HydrostaticStepError> {
        step_bodies_with(
            self.gravity(),
            bodies,
            states,
            controls,
            &mut cancelled,
            |states, cancel| self.load_bodies(bodies, states, cancel),
        )
    }
}

/// Shared native motion and energy gate for uniform or layered finite liquids.
pub(in crate::hydrostatic) fn step_bodies_with<C: FnMut() -> bool>(
    gravity: Gravity,
    bodies: &[FiniteTankBody<'_>],
    states: &[RigidBodyState],
    controls: HydrostaticStepControls,
    cancelled: &mut C,
    mut load: impl FnMut(&[RigidBodyState], &mut C) -> Result<FiniteTankBodiesLoad, HydrostaticError>,
) -> Result<FiniteTankBodiesStep, HydrostaticStepError> {
    poll(cancelled)?;
    if !controls.duration_s.is_finite()
        || controls.duration_s <= 0.0
        || !controls.maximum_energy_defect_j.is_finite()
        || controls.maximum_energy_defect_j <= 0.0
    {
        return Err(HydrostaticError::Invalid(
            "positive finite timestep and absolute energy allowance required",
        )
        .into());
    }
    let load_before = load(states, cancelled)?;
    let integrator = RigidBodyIntegrator::new(gravity);
    let mut stages = Vec::with_capacity(bodies.len());
    let mut drifted = Vec::with_capacity(bodies.len());
    for ((body, &state), load) in bodies.iter().zip(states).zip(&load_before.bodies) {
        poll(cancelled)?;
        let stage = integrator
            .begin_verlet(
                state,
                body.frame.properties(),
                load.wrench(),
                controls.duration_s,
            )
            .map_err(HydrostaticError::from)?;
        drifted.push(stage.drifted_state());
        stages.push(stage);
    }
    let load_after = load(&drifted, cancelled)?;
    let mut dynamics = Vec::with_capacity(bodies.len());
    let mut before = 0.0;
    let mut after = 0.0;
    let mut defect = 0.0;
    for (stage, load) in stages.into_iter().zip(&load_after.bodies) {
        poll(cancelled)?;
        let step = stage
            .finish(load.wrench())
            .map_err(HydrostaticError::from)?;
        let a = step.diagnostics_before;
        let b = step.diagnostics_after;
        before = finite(before + a.mechanical_energy)?;
        after = finite(after + b.mechanical_energy)?;
        defect = finite(
            defect
                + ((b.translational_kinetic_energy - a.translational_kinetic_energy)
                    + (b.rotational_kinetic_energy - a.rotational_kinetic_energy)
                    + (b.gravitational_potential_energy - a.gravitational_potential_energy)),
        )?;
        dynamics.push(step);
    }
    let fluid_before = load_before.surface.liquid_gravitational_energy_change_j;
    let fluid_after = load_after.surface.liquid_gravitational_energy_change_j;
    let total_energy_before_j = finite(before + fluid_before)?;
    let total_energy_after_j = finite(after + fluid_after)?;
    let energy_defect_j = finite(defect + (fluid_after - fluid_before))?;
    if energy_defect_j.abs() > controls.maximum_energy_defect_j {
        return Err(HydrostaticStepError::EnergyDefect {
            defect_j: energy_defect_j,
            maximum_j: controls.maximum_energy_defect_j,
        });
    }
    poll(cancelled)?;
    Ok(FiniteTankBodiesStep {
        dynamics,
        load_before,
        load_after,
        total_energy_before_j,
        total_energy_after_j,
        energy_defect_j,
    })
}
