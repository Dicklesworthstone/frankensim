//! Native Verlet stages for configuration-coupled endpoint forces.
use super::{
    DynamicsDiagnostics, DynamicsError, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState,
    StepReceipt, Vec3, Wrench, finite_derived, free_rotation_step,
};
use fs_math::dd::Dd;

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
    // The public stage view rounds to the original f64 state. Keep its linear
    // kick tail until position contraction and the final kick are complete;
    // two separately rounded half-kicks can otherwise lose a full impulse.
    linear_momentum_mid: [Dd; 3],
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
        let momentum_before = components(state.linear_momentum_world);
        let force_before = components(wrench_before.force_world.add(weight_world));
        let linear_momentum_mid = core::array::from_fn(|d| {
            Dd::from_f64(momentum_before[d]) + Dd::from_f64(force_before[d]) * Dd::from_f64(half)
        });
        let linear_mid = rounded(linear_momentum_mid, "verlet.linear_momentum_mid")?;
        let angular_kicked = state
            .angular_momentum_body
            .add(wrench_before.torque_body.scale(half));
        let (orientation_after, angular_drifted) = free_rotation_step(
            state.pose.orientation,
            angular_kicked,
            properties,
            duration_seconds,
        )?;
        let position_before = components(state.pose.position_world);
        let quotient = finite_derived(duration_seconds / properties.mass, "verlet.drift_scale")?;
        // Refine the finite quotient with one fused remainder. Reconstructing
        // mass*quotient separately can overflow even when dt/mass is finite
        // (for example f64::MAX / 3). Underflow retains the original scope.
        let remainder = (-quotient).mul_add(properties.mass, duration_seconds);
        let correction = finite_derived(remainder / properties.mass, "verlet.drift_scale")?;
        let drift_scale = Dd::from_pair(quotient, correction);
        let position_after = rounded(
            core::array::from_fn(|d| {
                Dd::from_f64(position_before[d]) + linear_momentum_mid[d] * drift_scale
            }),
            "verlet.position_after",
        )?;
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
            linear_momentum_mid,
            weight_world,
            duration_seconds,
            half,
        })
    }
}

impl VerletStage {
    /// Endpoint pose with rounded half-kicked momenta, for configuration-dependent
    /// loads. The private stage retains the linear kick tail for [`Self::finish`].
    pub const fn drifted_state(&self) -> RigidBodyState {
        self.drifted
    }

    /// Validate the endpoint wrench and finish the original native half-kick.
    /// Force excludes gravity; torque uses this endpoint's body coordinates.
    /// Invalid arithmetic/state/force returns no completed receipt.
    pub fn finish(self, wrench_after: Wrench) -> Result<StepReceipt, DynamicsError> {
        wrench_after.validate()?;
        let force_after = components(wrench_after.force_world.add(self.weight_world));
        let linear_momentum_after = rounded(
            core::array::from_fn(|d| {
                self.linear_momentum_mid[d] + Dd::from_f64(force_after[d]) * Dd::from_f64(self.half)
            }),
            "verlet.linear_momentum_after",
        )?;
        let state_after = RigidBodyState::new(
            self.drifted.pose,
            linear_momentum_after,
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

fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

fn rounded(values: [Dd; 3], field: &'static str) -> Result<Vec3, DynamicsError> {
    for value in values {
        finite_derived(value.hi, field)?;
        finite_derived(value.lo, field)?;
    }
    let [x, y, z] = values.map(Dd::to_f64);
    Ok(Vec3::new(x, y, z))
}
