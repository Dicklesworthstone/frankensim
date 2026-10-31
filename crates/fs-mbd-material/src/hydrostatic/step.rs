//! Re-evaluated hydrostatic loading around the existing rigid-body time step.
use super::{DisplacementMesh, HydrostaticError, HydrostaticLoad, HydrostaticMedium, finite, poll};
use crate::principal::{PrincipalBody, PrincipalFrame};
use fs_material::fluid_transport::FluidTransport;
use fs_mbd::{Gravity, RigidBodyIntegrator, RigidBodyState, StepReceipt};

/// One explicit physical timestep and an absolute numerical-energy allowance.
/// This is not a source uncertainty, stability proof or physical damping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HydrostaticStepControls {
    /// Positive physical interval [s].
    pub duration_s: f64,
    /// Maximum absolute change in kinetic + gravity + pressure potential [J].
    /// Relative normalization to a large arbitrary potential datum is not used.
    pub maximum_energy_defect_j: f64,
}

/// Failed candidates never replace the caller's accepted state.
#[derive(Clone, Debug, PartialEq)]
pub enum HydrostaticStepError {
    /// Existing material, geometry, frame, dynamics or cancellation refusal.
    Physics(HydrostaticError),
    /// Measured energy error exceeded the caller's absolute budget.
    EnergyDefect {
        /// Signed actual kinetic + gravity + pressure-potential change [J].
        defect_j: f64,
        /// Caller allowance [J]. Reduce the physical step to improve accuracy.
        maximum_j: f64,
    },
}
impl core::fmt::Display for HydrostaticStepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "hydrostatic time step: {self:?}")
    }
}
impl std::error::Error for HydrostaticStepError {}
impl From<HydrostaticError> for HydrostaticStepError {
    fn from(e: HydrostaticError) -> Self {
        Self::Physics(e)
    }
}

/// Actual rigid-body step plus re-evaluated loads at both accepted endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HydrostaticStep {
    /// Original integrator's before/after states and dynamics diagnostics.
    pub dynamics: StepReceipt,
    /// Force used for the initial half-kick, from the original state.
    pub load_before: HydrostaticLoad,
    /// Force at the returned geometry, used for the final half-kick.
    pub load_after: HydrostaticLoad,
    /// Kinetic + body gravity + fluid-displacement potential before the step [J].
    pub total_energy_before_j: f64,
    /// The same energy at the actual candidate state [J].
    pub total_energy_after_j: f64,
    /// Signed numerical energy change, NOT physical fluid dissipation [J].
    pub energy_defect_j: f64,
}

impl<F: FluidTransport + ?Sized> HydrostaticMedium<'_, F> {
    /// Advance a material-defined rigid body under gravity and buoyancy alone.
    /// Reuses `RigidBodyIntegrator::step_verlet_with`; no new momentum, attitude
    /// or contact integrator is introduced. Buoyancy is recalculated from the
    /// actual immersed geometry every step, so source mass changes the motion
    /// and the resulting motion changes both force and torque.
    ///
    /// This is a second-order symmetric endpoint-load scheme. The original
    /// momentum-preserving free-rotation map lies between two force/torque
    /// half-kicks; buoyancy at the actual drifted pose supplies the second.
    /// The hydrostatic reservoir is conservative, but its time discretization
    /// is NOT energy-exact. The
    /// explicit absolute defect gate measures K + U_gravity + U_pressure;
    /// passing a local gate is not a global stability/accuracy certificate.
    /// There is no fluid added mass, drag or wave radiation. Small-amplitude
    /// heave periods are those of this declared quasi-static-pressure model.
    ///
    /// Two bounded geometry passes surround the original rigid-body step.
    /// Cancellation, failed geometry/frame admission and an energy refusal
    /// return no candidate. The caller commits `dynamics.state_after` only on
    /// success. Keep the same body, envelope, fluid and surface while advancing:
    /// changing those data is model replacement, not a conservative event.
    pub fn step(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        controls: HydrostaticStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticStep, HydrostaticStepError> {
        self.step_frame(mesh, body.frame(), state, controls, &mut cancelled)
    }

    pub(crate) fn step_frame(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalFrame,
        state: RigidBodyState,
        controls: HydrostaticStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticStep, HydrostaticStepError> {
        step_with_load(
            self.gravity(),
            body,
            state,
            controls,
            &mut cancelled,
            |state, cancel| self.load_frame(mesh, body, state, cancel),
        )
    }
}

pub(super) fn step_with_load<C: FnMut() -> bool>(
    gravity: Gravity,
    body: &PrincipalFrame,
    state: RigidBodyState,
    controls: HydrostaticStepControls,
    cancelled: &mut C,
    mut load: impl FnMut(RigidBodyState, &mut C) -> Result<HydrostaticLoad, HydrostaticError>,
) -> Result<HydrostaticStep, HydrostaticStepError> {
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
    let load_before = load(state, cancelled)?;
    poll(cancelled)?;
    let (dynamics, load_after) = RigidBodyIntegrator::new(gravity).step_verlet_with(
        state,
        body.properties(),
        load_before.wrench(),
        controls.duration_s,
        |drifted| {
            let after = load(drifted, cancelled)?;
            Ok::<_, HydrostaticError>((after.wrench(), after))
        },
    )?;
    let before = dynamics.diagnostics_before;
    let after = dynamics.diagnostics_after;
    let total_energy_before_j =
        finite(before.mechanical_energy + load_before.pressure_potential_j())?;
    let total_energy_after_j = finite(after.mechanical_energy + load_after.pressure_potential_j())?;
    // Separate storage differences before summing. Do not first subtract
    // totals containing a potentially large common gravity reference.
    let energy_defect_j = finite(
        (after.translational_kinetic_energy - before.translational_kinetic_energy)
            + (after.rotational_kinetic_energy - before.rotational_kinetic_energy)
            + (after.gravitational_potential_energy - before.gravitational_potential_energy)
            + (load_after.pressure_potential_j() - load_before.pressure_potential_j()),
    )?;
    if energy_defect_j.abs() > controls.maximum_energy_defect_j {
        return Err(HydrostaticStepError::EnergyDefect {
            defect_j: energy_defect_j,
            maximum_j: controls.maximum_energy_defect_j,
        });
    }
    poll(cancelled)?;
    Ok(HydrostaticStep {
        dynamics,
        load_before,
        load_after,
        total_energy_before_j,
        total_energy_after_j,
        energy_defect_j,
    })
}
