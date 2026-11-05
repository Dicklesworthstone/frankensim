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

/// Work-conjugate Cayley stepping of the original prepared spatial loads.
pub mod midpoint;

/// Momentum-conserving internal wrench impulses between original finite bodies.
pub mod internal_impulse;

/// Finite internally actuated trajectories on two original native rigid bodies.
pub mod internal_motion;

use super::{
    AdvanceOutcome, DynamicsError, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState,
    StepReceipt, UnitQuaternion, Vec3, Wrench, finite_derived,
};
use fs_math::dd::Dd;

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
    point_force_moments_world: [CompensatedVec3; 3],
}

// Keep the product tails as well as the reduction tails until a physical
// resultant is formed. Rounding r_j F first can erase an entire finite couple:
// (2^27 + 1)(2^27 - 1) - 2^54 is -1, although both rounded products agree.
// Three of these vectors suffice at every attitude, regardless of point count.
#[derive(Clone, Copy, Debug, PartialEq)]
struct CompensatedVec3([Dd; 3]);

impl CompensatedVec3 {
    const ZERO: Self = Self([Dd::ZERO; 3]);

    fn from_vec(value: Vec3) -> Self {
        Self(components(value).map(Dd::from_f64))
    }

    fn checked(self, field: &'static str) -> Result<Self, DynamicsError> {
        for component in self.0 {
            finite_derived(component.hi, field)?;
            finite_derived(component.lo, field)?;
        }
        Ok(self)
    }

    fn rounded(self, field: &'static str) -> Result<Vec3, DynamicsError> {
        let [x, y, z] = self.checked(field)?.0.map(Dd::to_f64);
        Ok(Vec3::new(x, y, z))
    }

    fn add(self, other: Self) -> Self {
        Self(core::array::from_fn(|j| self.0[j] + other.0[j]))
    }

    fn sub(self, other: Self) -> Self {
        Self(core::array::from_fn(|j| self.0[j] - other.0[j]))
    }

    fn scale(self, value: f64) -> Self {
        Self(self.0.map(|component| component * Dd::from_f64(value)))
    }

    fn dot(self, other: Self) -> Dd {
        self.0[0] * other.0[0] + self.0[1] * other.0[1] + self.0[2] * other.0[2]
    }

    fn cross(self, other: Self) -> Self {
        Self(core::array::from_fn(|j| {
            let k = (j + 1) % 3;
            let l = (j + 2) % 3;
            self.0[k] * other.0[l] - self.0[l] * other.0[k]
        }))
    }

    fn expressed_in(self, axes: [Vec3; 3]) -> Self {
        Self(axes.map(|axis| self.dot(Self::from_vec(axis))))
    }

    fn norm_bound(self, field: &'static str) -> Result<f64, DynamicsError> {
        let [hx, hy, hz] = self.0.map(|component| component.hi);
        let [lx, ly, lz] = self.0.map(|component| component.lo);
        finite_derived(
            Vec3::new(hx, hy, hz).stable_norm(field)? + Vec3::new(lx, ly, lz).stable_norm(field)?,
            field,
        )
    }
}

impl PreparedSpatialLoads {
    /// No external force or torque. Integrator gravity is still applied once.
    pub const ZERO: Self = Self {
        force_world: Vec3::ZERO,
        torque_body: Vec3::ZERO,
        torque_world: Vec3::ZERO,
        point_force_moments_world: [CompensatedVec3::ZERO; 3],
    };

    /// Prepare held COM force, body-held couple, world-held couple, and point
    /// forces. `held` uses the original [`Wrench`] conventions. Gravity belongs
    /// on the integrator and must not also be included here.
    ///
    /// Products and reductions are compensated in supplied order and use no
    /// heap storage. First-moment tails are retained through torque and work
    /// evaluation so opposing large point loads can retain a small couple.
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
        let mut force_world = CompensatedVec3::from_vec(held.force_world);
        let mut point_force_moments_world = [CompensatedVec3::ZERO; 3];
        for (index, point) in points.iter().enumerate() {
            if index % 64 == 0 && cancelled() {
                return Err(SpatialLoadError::Cancelled);
            }
            finite_vec(point.arm_body, "spatial_load.arm_body")?;
            finite_vec(point.force_world, "spatial_load.force_world")?;
            let force = CompensatedVec3::from_vec(point.force_world);
            force_world = force_world.add(force).checked("spatial_load.reduction")?;
            for (j, arm) in components(point.arm_body).into_iter().enumerate() {
                let moment = force.scale(arm).checked("spatial_load.first_moment")?;
                point_force_moments_world[j] = point_force_moments_world[j]
                    .add(moment)
                    .checked("spatial_load.reduction")?;
            }
        }
        if cancelled() {
            return Err(SpatialLoadError::Cancelled);
        }
        Ok(Self {
            force_world: force_world.rounded("spatial_load.reduction")?,
            torque_body: held.torque_body,
            torque_world,
            point_force_moments_world,
        })
    }

    /// Resultant external force in N, excluding integrator gravity.
    #[must_use]
    pub const fn force_world(&self) -> Vec3 {
        self.force_world
    }

    /// Orientation-independent upper bound on the applied COM torque [N m].
    ///
    /// The triangle bound includes both free couples and the three prepared
    /// point-force moments. It is reusable without revisiting the input force
    /// set. Constrained steppers use it to bound angular travel before trying
    /// a new attitude; freezing the initial torque magnitude is not a bound.
    ///
    /// # Errors
    /// Refuses a finite load set whose bound cannot be represented in `f64`.
    pub fn torque_norm_bound(&self) -> Result<f64, DynamicsError> {
        let mut bound = self.torque_body.stable_norm("spatial_load.torque_bound")?;
        bound = finite_derived(
            bound + self.torque_world.stable_norm("spatial_load.torque_bound")?,
            "spatial_load.torque_bound",
        )?;
        for vector in self.point_force_moments_world {
            bound = finite_derived(
                bound + vector.norm_bound("spatial_load.torque_bound")?,
                "spatial_load.torque_bound",
            )?;
        }
        Ok(bound)
    }

    /// Instantaneous native wrench at the specified attitude.
    ///
    /// This query is useful for force/power inspection. Freezing this returned
    /// body-frame moment over a finite step changes the physical load; use
    /// [`RigidBodyIntegrator::step_spatial`] to evolve the original loads.
    pub fn wrench_at(&self, orientation: UnitQuaternion) -> Result<Wrench, DynamicsError> {
        let axes = basis_world(orientation);
        let mut torque_world = CompensatedVec3::from_vec(self.torque_world);
        for (axis, moment) in axes.into_iter().zip(self.point_force_moments_world) {
            torque_world = torque_world.add(CompensatedVec3::from_vec(axis).cross(moment));
        }
        let torque_body = torque_world
            .expressed_in(axes)
            .add(CompensatedVec3::from_vec(self.torque_body))
            .rounded("spatial_load.torque_body")?;
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
        let mut work = CompensatedVec3::from_vec(self.force_world).dot(
            CompensatedVec3::from_vec(after.position_world)
                .sub(CompensatedVec3::from_vec(before.position_world)),
        );
        let axes_before = basis_world(before.orientation);
        let axes_after = basis_world(after.orientation);
        for j in 0..3 {
            work = work
                + self.point_force_moments_world[j].dot(
                    CompensatedVec3::from_vec(axes_after[j])
                        .sub(CompensatedVec3::from_vec(axes_before[j])),
                );
        }
        finite_derived(work.lo, "spatial_load.force_work")?;
        finite_derived(work.to_f64(), "spatial_load.force_work")
    }

    fn free_couple_body(&self, orientation: UnitQuaternion) -> Result<Vec3, DynamicsError> {
        CompensatedVec3::from_vec(self.torque_world)
            .expressed_in(basis_world(orientation))
            .add(CompensatedVec3::from_vec(self.torque_body))
            .rounded("spatial_load.free_couple")
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
    /// from the selected integrator's load quadrature, in N m s.
    pub angular_impulse_world: Vec3,
    /// Actual linear momentum change minus the applied impulse, in N s.
    pub linear_momentum_residual_world: Vec3,
    /// Actual total angular momentum change minus the applied impulse, in N m s.
    pub angular_momentum_residual_world: Vec3,
    /// Exact endpoint work of the held forces, excluding gravity, in J.
    pub force_work_j: f64,
    /// Free-couple work in J: half-kick power quadrature for `step_spatial`, or
    /// the conjugate Cayley increment for `step_spatial_midpoint`. Both are
    /// second-order estimates of continuous free-couple work.
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
