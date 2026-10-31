//! Spatial loads on the original full-inertia rigid-body state.
//!
//! A force fixed in the world and applied to a body-fixed point has moment
//! `(R r) cross F`. Rotating its initial moment into the body once would instead
//! create a follower torque. Both endpoint kicks here use the actual attitude,
//! around the existing Hamiltonian free-rotation map.
//!
//! Preparation reduces any admitted point-force set to its resultant and three
//! first-moment vectors `C_j = sum_i r_ij F_i`. The exact moment at ANY attitude
//! is `sum_j (R e_j) cross C_j`. Thus stepping has constant storage and work,
//! without per-point allocations, traversal, or an inertia factorization.
//! The same moments give the exact work of the held forces between two poses.
//!
//! Force/moment and power conventions follow Lynch and Park, Modern Robotics,
//! section 3.4: <https://modernrobotics.northwestern.edu/nu-gm-book-resource/3-4-wrenches/>.
//! This is a smooth held-load model. Drag, contact, springs, changing fields,
//! deforming attachments and other state-dependent loads require their own
//! stage evaluations; preparing them once does not make their motion second order.

use super::{
    AdvanceOutcome, DynamicsError, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState,
    StepReceipt, UnitQuaternion, Vec3, Wrench, finite_derived,
};

/// Maximum point-force count admitted by one bounded preparation.
pub const MAX_WORLD_POINT_FORCES: usize = 16_384;

/// A world-held force applied at a material point of the rigid body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WorldPointForce {
    /// Arm from the center of mass, in the declared body coordinates, in m.
    pub arm_body: Vec3,
    /// Force fixed in the world frame, in N.
    pub force_world: Vec3,
}

/// A refusal leaves no partially prepared load or modified body state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpatialLoadError {
    /// The native dynamics owner refused a non-finite or unrepresentable value.
    Dynamics(DynamicsError),
    /// The supplied point-force count exceeds the preparation bound.
    TooManyPointForces,
    /// Cancellation observed before preparation completed.
    Cancelled,
}

impl From<DynamicsError> for SpatialLoadError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}

impl core::fmt::Display for SpatialLoadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Dynamics(error) => write!(f, "spatial load: {error}"),
            Self::TooManyPointForces => {
                write!(f, "spatial load exceeds {MAX_WORLD_POINT_FORCES} points")
            }
            Self::Cancelled => f.write_str("spatial load preparation cancelled"),
        }
    }
}

impl std::error::Error for SpatialLoadError {}

/// Constant-size, immutable loads, reusable at every orientation of one body.
///
/// Torques are free couples about the center of mass, in their named frames.
/// The body axes must be the same axes used by the original mass properties,
/// angular momentum and material-point geometry; principal axes are not inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedSpatialLoads {
    force_world: Vec3,
    torque_body: Vec3,
    torque_world: Vec3,
    point_force_moments_world: [Vec3; 3],
}

impl PreparedSpatialLoads {
    /// No external force or torque. Integrator gravity is still applied once.
    pub const ZERO: Self = Self {
        force_world: Vec3::ZERO,
        torque_body: Vec3::ZERO,
        torque_world: Vec3::ZERO,
        point_force_moments_world: [Vec3::ZERO; 3],
    };

    /// Prepare held COM force, body-held couple, world-held couple, and point
    /// forces. `held` uses the original [`Wrench`] conventions. Gravity belongs
    /// on the integrator and must not also be included here.
    ///
    /// Reductions are compensated in supplied order and use no heap storage.
    /// Cancellation is checked initially, every 64 points, and before return.
    /// All inputs and products must be finite and representable.
    pub fn new(
        held: Wrench,
        torque_world: Vec3,
        points: &[WorldPointForce],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SpatialLoadError> {
        if cancelled() {
            return Err(SpatialLoadError::Cancelled);
        }
        if points.len() > MAX_WORLD_POINT_FORCES {
            return Err(SpatialLoadError::TooManyPointForces);
        }
        held.validate()?;
        finite_vec(torque_world, "spatial_load.torque_world")?;
        let mut sums = [0.0; 12];
        let mut corrections = [0.0; 12];
        sums[..3].copy_from_slice(&components(held.force_world));
        for (index, point) in points.iter().enumerate() {
            if index % 64 == 0 && cancelled() {
                return Err(SpatialLoadError::Cancelled);
            }
            finite_vec(point.arm_body, "spatial_load.arm_body")?;
            finite_vec(point.force_world, "spatial_load.force_world")?;
            let force = components(point.force_world);
            let arm = components(point.arm_body);
            let mut values = [0.0; 12];
            values[..3].copy_from_slice(&force);
            for j in 0..3 {
                for k in 0..3 {
                    values[3 + 3 * j + k] = arm[j] * force[k];
                }
            }
            for k in 0..12 {
                let value = finite_derived(values[k], "spatial_load.first_moment")?;
                let old = sums[k];
                let next = finite_derived(old + value, "spatial_load.reduction")?;
                corrections[k] += if old.abs() >= value.abs() {
                    (old - next) + value
                } else {
                    (value - next) + old
                };
                sums[k] = next;
            }
        }
        for k in 0..12 {
            sums[k] = finite_derived(sums[k] + corrections[k], "spatial_load.reduction")?;
        }
        if cancelled() {
            return Err(SpatialLoadError::Cancelled);
        }
        Ok(Self {
            force_world: Vec3::new(sums[0], sums[1], sums[2]),
            torque_body: held.torque_body,
            torque_world,
            point_force_moments_world: core::array::from_fn(|j| {
                Vec3::new(sums[3 + 3 * j], sums[4 + 3 * j], sums[5 + 3 * j])
            }),
        })
    }

    /// Resultant external force in N, excluding integrator gravity.
    #[must_use]
    pub const fn force_world(&self) -> Vec3 {
        self.force_world
    }

    /// Instantaneous native wrench at the specified attitude.
    ///
    /// This query is useful for force/power inspection. Freezing this returned
    /// body-frame moment over a finite step changes the physical load; use
    /// [`RigidBodyIntegrator::step_spatial`] to evolve the original loads.
    pub fn wrench_at(&self, orientation: UnitQuaternion) -> Result<Wrench, DynamicsError> {
        let couple = self.free_couple_body(orientation)?;
        let mut point_torque = Vec3::ZERO;
        for (axis, moment) in basis_world(orientation)
            .into_iter()
            .zip(self.point_force_moments_world)
        {
            point_torque = point_torque.add(axis.cross(moment));
        }
        let torque_body = couple.add(orientation.rotate_world_to_body(point_torque));
        finite_vec(torque_body, "spatial_load.torque_body")?;
        Ok(Wrench {
            force_world: self.force_world,
            torque_body,
        })
    }

    /// Exact work [J] of the held COM and material-point forces between poses.
    /// Free-couple work and uniform gravitational work are excluded.
    pub fn force_work(&self, before: Pose, after: Pose) -> Result<f64, DynamicsError> {
        before.validate()?;
        after.validate()?;
        let mut work = self
            .force_world
            .dot(after.position_world.sub(before.position_world));
        let axes_before = basis_world(before.orientation);
        let axes_after = basis_world(after.orientation);
        for j in 0..3 {
            work += self.point_force_moments_world[j].dot(axes_after[j].sub(axes_before[j]));
        }
        finite_derived(work, "spatial_load.force_work")
    }

    fn free_couple_body(&self, orientation: UnitQuaternion) -> Result<Vec3, DynamicsError> {
        let result = self
            .torque_body
            .add(orientation.rotate_world_to_body(self.torque_world));
        finite_vec(result, "spatial_load.free_couple")?;
        Ok(result)
    }
}

/// Accepted native motion and independently reconstructable load balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialStepReceipt {
    /// Original before/after state and mechanical-energy/momentum diagnostics.
    pub step: StepReceipt,
    /// Applied linear impulse, including uniform gravity, in N s.
    pub linear_impulse_world: Vec3,
    /// Applied angular impulse about the fixed world origin, including gravity,
    /// from the two actual endpoint kicks, in N m s.
    pub angular_impulse_world: Vec3,
    /// Actual linear momentum change minus the applied impulse, in N s.
    pub linear_momentum_residual_world: Vec3,
    /// Actual total angular momentum change minus the applied impulse, in N m s.
    pub angular_momentum_residual_world: Vec3,
    /// Exact endpoint work of the held forces, excluding gravity, in J.
    pub force_work_j: f64,
    /// Free-couple work from the two half-kick power quadratures, in J.
    /// This is a second-order estimate, not an exact continuous work integral.
    pub couple_work_j: f64,
    /// Mechanical-energy change minus force and couple work, in J. This is a
    /// numerical defect, never physical friction, heat, or hidden damping.
    pub work_energy_residual_j: f64,
}

impl RigidBodyIntegrator {
    /// Advance the original full-inertia body under spatially held loads.
    ///
    /// The original [`Self::step_verlet_with`] owner supplies the endpoint
    /// momentum half-kicks and full-inertia free-rotation map. Torque
    /// is reevaluated at each actual endpoint orientation, including rotating
    /// material arms. A pure world-held couple gives `delta(R L) = dt tau`
    /// to roundoff. Point forces use the symmetric endpoint moment quadrature.
    /// Translation is exact under the constant resultant and uniform gravity.
    ///
    /// The composition is second order for these held loads; conservative point
    /// forces retain Hamiltonian splitting. Free couples need not be conservative.
    /// Exact forced energy conservation, contacts and constraints are not claimed.
    /// Fixed-size work is allocation-free; no nonlinear solve or new inertia
    /// factorization occurs. Refusal returns no partially advanced state.
    pub fn step_spatial(
        self,
        state: RigidBodyState,
        properties: MassProperties,
        loads: &PreparedSpatialLoads,
        duration_seconds: f64,
    ) -> Result<SpatialStepReceipt, DynamicsError> {
        let wrench_before = loads.wrench_at(state.pose.orientation)?;
        let (step, (drifted_state, wrench_after)) = self.step_verlet_with(
            state,
            properties,
            wrench_before,
            duration_seconds,
            |drifted: RigidBodyState| -> Result<_, DynamicsError> {
                let wrench = loads.wrench_at(drifted.pose.orientation)?;
                Ok((wrench, (drifted, wrench)))
            },
        )?;
        let before = step.diagnostics_before;
        let after = step.diagnostics_after;
        let state_after = step.state_after;
        let half = 0.5 * duration_seconds;
        let force = loads
            .force_world
            .add(self.gravity.acceleration_world().scale(properties.mass));
        finite_vec(force, "spatial_step.total_force")?;
        let linear_impulse = force.scale(duration_seconds);
        let position_after = state_after.pose.position_world;
        let orientation_before = state.pose.orientation;
        let orientation_after = state_after.pose.orientation;
        let torque_before = wrench_before.torque_body;
        let torque_after = wrench_after.torque_body;
        let kicked = state.angular_momentum_body.add(torque_before.scale(half));
        let drifted = drifted_state.angular_momentum_body;
        let momentum_after = state_after.angular_momentum_body;

        let couple_before = loads.free_couple_body(orientation_before)?;
        let couple_after = loads.free_couple_body(orientation_after)?;
        let omega_first_kick = properties.angular_velocity_body_checked(
            state
                .angular_momentum_body
                .scale(0.5)
                .add(kicked.scale(0.5)),
        )?;
        let omega_last_kick = properties
            .angular_velocity_body_checked(drifted.scale(0.5).add(momentum_after.scale(0.5)))?;
        let couple_work_j = finite_derived(
            half * couple_before.dot(omega_first_kick) + half * couple_after.dot(omega_last_kick),
            "spatial_step.couple_work",
        )?;
        let force_work_j = loads.force_work(state.pose, state_after.pose)?;
        let work_energy_residual_j = finite_derived(
            after.mechanical_energy - before.mechanical_energy - force_work_j - couple_work_j,
            "spatial_step.work_energy_residual",
        )?;
        let angular_impulse = orientation_before
            .rotate_body_to_world(torque_before)
            .scale(half)
            .add(
                orientation_after
                    .rotate_body_to_world(torque_after)
                    .scale(half),
            )
            .add(
                state
                    .pose
                    .position_world
                    .scale(0.5)
                    .add(position_after.scale(0.5))
                    .cross(linear_impulse),
            );
        let linear_residual = after
            .linear_momentum_world
            .sub(before.linear_momentum_world)
            .sub(linear_impulse);
        let angular_residual = after
            .angular_momentum_world
            .sub(before.angular_momentum_world)
            .sub(angular_impulse);
        finite_vec(linear_impulse, "spatial_step.linear_impulse")?;
        finite_vec(angular_impulse, "spatial_step.angular_impulse")?;
        finite_vec(linear_residual, "spatial_step.linear_residual")?;
        finite_vec(angular_residual, "spatial_step.angular_residual")?;
        Ok(SpatialStepReceipt {
            step,
            linear_impulse_world: linear_impulse,
            angular_impulse_world: angular_impulse,
            linear_momentum_residual_world: linear_residual,
            angular_momentum_residual_world: angular_residual,
            force_work_j,
            couple_work_j,
            work_energy_residual_j,
        })
    }

    /// Advance complete spatial-load steps, checking cancellation before each
    /// bounded step. The mutable state retains only successfully completed steps.
    pub fn advance_spatial(
        self,
        state: &mut RigidBodyState,
        properties: MassProperties,
        loads: &PreparedSpatialLoads,
        duration_seconds: f64,
        step_count: usize,
        mut is_cancelled: impl FnMut(usize) -> bool,
    ) -> Result<AdvanceOutcome, DynamicsError> {
        for completed_steps in 0..step_count {
            if is_cancelled(completed_steps) {
                return Ok(AdvanceOutcome::Cancelled {
                    completed_steps,
                    final_diagnostics: self.diagnostics(*state, properties)?,
                });
            }
            *state = self
                .step_spatial(*state, properties, loads, duration_seconds)?
                .step
                .state_after;
        }
        Ok(AdvanceOutcome::Completed {
            completed_steps: step_count,
            final_diagnostics: self.diagnostics(*state, properties)?,
        })
    }
}

fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

fn basis_world(orientation: UnitQuaternion) -> [Vec3; 3] {
    [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]
    .map(|axis| orientation.rotate_body_to_world(axis))
}

fn finite_vec(value: Vec3, field: &'static str) -> Result<(), DynamicsError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(DynamicsError::NonFinite(field))
    }
}
