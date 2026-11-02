//! Native Verlet stages for configuration-coupled endpoint forces.
use super::{
    DynamicsDiagnostics, DynamicsError, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState,
    StepReceipt, Vec3, Wrench, free_rotation_step,
};

/// An uncommitted native half-kick and full drift, with its original inputs.
///
/// A composing owner may prepare several bodies, evaluate one force field at
/// all their drifted poses, and then finish each stage. No original history is
/// modified. The composing owner must accept or refuse the complete batch.
/// Endpoint loads must depend only on configuration for the stated second-order
/// Verlet method; this stage does not supply a velocity-coupled force solve.
#[derive(Debug)]
pub struct VerletStage {
    integrator: RigidBodyIntegrator,
    properties: MassProperties,
    state_before: RigidBodyState,
    diagnostics_before: DynamicsDiagnostics,
    drifted: RigidBodyState,
    weight_world: Vec3,
    duration_seconds: f64,
    half: f64,
}

impl RigidBodyIntegrator {
    /// Prepare the original Verlet half-kick and native free translation/rotation.
    /// The initial wrench excludes gravity and uses the input body's coordinates.
    /// Neither this stage nor its endpoint view is a completed accepted step.
    pub fn begin_verlet(
        self,
        state: RigidBodyState,
        properties: MassProperties,
        wrench_before: Wrench,
        duration_seconds: f64,
    ) -> Result<VerletStage, DynamicsError> {
        if !duration_seconds.is_finite() || duration_seconds <= 0.0 {
            return Err(DynamicsError::InvalidStepDuration);
        }
        self.gravity.validate()?;
        state.validate()?;
        properties.validate()?;
        wrench_before.validate()?;
        let diagnostics_before = self.diagnostics(state, properties)?;
        let half = 0.5 * duration_seconds;
        let weight_world = self.gravity.acceleration_world().scale(properties.mass);
        let linear_mid = state
            .linear_momentum_world
            .add(wrench_before.force_world.add(weight_world).scale(half));
        let angular_kicked = state
            .angular_momentum_body
            .add(wrench_before.torque_body.scale(half));
        let (orientation_after, angular_drifted) = free_rotation_step(
            state.pose.orientation,
            angular_kicked,
            properties,
            duration_seconds,
        )?;
        let position_after = state
            .pose
            .position_world
            .add(linear_mid.scale(duration_seconds / properties.mass));
        let drifted = RigidBodyState::new(
            Pose::new(position_after, orientation_after)?,
            linear_mid,
            angular_drifted,
        )?;
        Ok(VerletStage {
            integrator: self,
            properties,
            state_before: state,
            diagnostics_before,
            drifted,
            weight_world,
            duration_seconds,
            half,
        })
    }
}

impl VerletStage {
    /// Endpoint pose with half-kicked momenta, for configuration-dependent loads.
    pub const fn drifted_state(&self) -> RigidBodyState {
        self.drifted
    }

    /// Validate the endpoint wrench and finish the original native half-kick.
    /// Force excludes gravity; torque uses this endpoint's body coordinates.
    /// Invalid arithmetic/state/force returns no completed receipt.
    pub fn finish(self, wrench_after: Wrench) -> Result<StepReceipt, DynamicsError> {
        wrench_after.validate()?;
        let state_after = RigidBodyState::new(
            self.drifted.pose,
            self.drifted.linear_momentum_world.add(
                wrench_after
                    .force_world
                    .add(self.weight_world)
                    .scale(self.half),
            ),
            self.drifted
                .angular_momentum_body
                .add(wrench_after.torque_body.scale(self.half)),
        )?;
        let diagnostics_after = self.integrator.diagnostics(state_after, self.properties)?;
        Ok(StepReceipt {
            duration_seconds: self.duration_seconds,
            state_before: self.state_before,
            state_after,
            diagnostics_before: self.diagnostics_before,
            diagnostics_after,
        })
    }
}
