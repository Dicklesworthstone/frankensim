//! Native Verlet stages for configuration-coupled endpoint forces.
use super::{
    DynamicsDiagnostics, DynamicsError, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState,
    StepReceipt, Vec3, Wrench, finite_derived, free_rotation_step_compensated,
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
    // The free-rotation owner uses this precision for both its attitude
    // phases and matched momentum transport. The endpoint kick must precede
    // rounding; transporting a discarded kick tail afterward is insufficient.
    angular_momentum_drifted: [Dd; 3],
    force_before_world: Vec3,
    weight_world: Vec3,
    duration_seconds: f64,
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
        let weight_world = self.gravity.acceleration_world().scale(properties.mass);
        let force_before_world = wrench_before.force_world.add(weight_world);
        let linear_momentum_mid = add_half_impulse(
            lifted(state.linear_momentum_world),
            force_before_world,
            duration_seconds,
        );
        let linear_mid = rounded(linear_momentum_mid, "verlet.linear_momentum_mid")?;
        let angular_kicked = add_half_impulse(
            lifted(state.angular_momentum_body),
            wrench_before.torque_body,
            duration_seconds,
        );
        let (orientation_after, angular_momentum_drifted) = free_rotation_step_compensated(
            state.pose.orientation,
            angular_kicked,
            properties,
            duration_seconds,
        )?;
        let position_after = drift_position(
            state.pose.position_world,
            linear_momentum_mid,
            properties.mass,
            duration_seconds,
            "verlet.drift_scale",
            "verlet.position_after",
        )?;
        let drifted = RigidBodyState::new(
            Pose::new(position_after, orientation_after)?,
            linear_mid,
            rounded(angular_momentum_drifted, "verlet.angular_momentum_drifted")?,
        )?;
        Ok(VerletStage {
            integrator: self,
            properties,
            state_before: state,
            diagnostics_before,
            drifted,
            linear_momentum_mid,
            angular_momentum_drifted,
            force_before_world,
            weight_world,
            duration_seconds,
        })
    }
}

impl VerletStage {
    /// Endpoint pose with rounded half-kicked momenta, for configuration-dependent
    /// loads. The private stage retains linear and angular kick tails for
    /// [`Self::finish`]; the angular tail also drives the native rotation phases.
    pub const fn drifted_state(&self) -> RigidBodyState {
        self.drifted
    }

    /// Validate the endpoint wrench and finish the original native half-kick.
    /// Force excludes gravity; torque uses this endpoint's body coordinates.
    /// Invalid arithmetic/state/force returns no completed receipt.
    pub fn finish(self, wrench_after: Wrench) -> Result<StepReceipt, DynamicsError> {
        wrench_after.validate()?;
        let force_after_world = wrench_after.force_world.add(self.weight_world);
        let held = if force_after_world == self.force_before_world {
            let full = add_impulse(
                lifted(self.state_before.linear_momentum_world),
                force_after_world,
                self.duration_seconds,
            );
            // A complete held impulse may be a minimum subnormal although
            // neither half is representable even in DD. Preserve that endpoint
            // without narrowing stages whose full product would overflow.
            full.iter()
                .all(|p| p.hi.is_finite() && p.lo.is_finite())
                .then_some(full)
        } else {
            None
        };
        let linear_momentum_after = rounded(
            held.unwrap_or_else(|| {
                add_half_impulse(
                    self.linear_momentum_mid,
                    force_after_world,
                    self.duration_seconds,
                )
            }),
            "verlet.linear_momentum_after",
        )?;
        let angular_momentum_after = rounded(
            add_half_impulse(
                self.angular_momentum_drifted,
                wrench_after.torque_body,
                self.duration_seconds,
            ),
            "verlet.angular_momentum_after",
        )?;
        let state_after = RigidBodyState::new(
            self.drifted.pose,
            linear_momentum_after,
            angular_momentum_after,
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

pub(super) fn lifted(value: Vec3) -> [Dd; 3] {
    components(value).map(Dd::from_f64)
}

pub(super) fn add_impulse(momentum: [Dd; 3], force: Vec3, duration: f64) -> [Dd; 3] {
    let force = components(force);
    core::array::from_fn(|d| momentum[d] + Dd::from_f64(force[d]) * Dd::from_f64(duration))
}

pub(super) fn add_half_impulse(momentum: [Dd; 3], force: Vec3, duration: f64) -> [Dd; 3] {
    let half = 0.5 * duration;
    let result = add_impulse(momentum, force, half);
    let remainder = (-2.0_f64).mul_add(half, duration);
    if remainder == 0.0 {
        return result;
    }
    // An odd subnormal duration loses half a minimum subnormal in `half`.
    // Contract that residual with the force before halving it, without ever
    // forming a possibly overflowing full force*duration product.
    let force = components(force);
    core::array::from_fn(|d| {
        result[d] + (Dd::from_f64(force[d]) * Dd::from_f64(remainder)) * Dd::from_f64(0.5)
    })
}

pub(super) fn drift_position(
    position: Vec3,
    momentum: [Dd; 3],
    mass: f64,
    duration: f64,
    scale_field: &'static str,
    position_field: &'static str,
) -> Result<Vec3, DynamicsError> {
    let position = components(position);
    let quotient = finite_derived(duration / mass, scale_field)?;
    // A fused remainder keeps finite dt/m valid even when reconstructing
    // mass*quotient would overflow. Retain the original underflow scope.
    let remainder = (-quotient).mul_add(mass, duration);
    let correction = finite_derived(remainder / mass, scale_field)?;
    let scale = Dd::from_pair(quotient, correction);
    rounded(
        core::array::from_fn(|d| Dd::from_f64(position[d]) + momentum[d] * scale),
        position_field,
    )
}

pub(super) fn rounded(values: [Dd; 3], field: &'static str) -> Result<Vec3, DynamicsError> {
    for value in values {
        finite_derived(value.hi, field)?;
        finite_derived(value.lo, field)?;
    }
    let [x, y, z] = values.map(Dd::to_f64);
    Ok(Vec3::new(x, y, z))
}
