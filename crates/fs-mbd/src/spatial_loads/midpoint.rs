//! Held point forces conjugate to actual finite material-point displacement.
//!
//! With `c = dt I^-1 Lbar`, `R1 = R0 cay(c)`, and `Rbar = (R0+R1)/2`,
//! `R1-R0 = Rbar [c]x`. The mean spatial moment of a held point force is
//! `(Rbar r) cross F`; pulling it back with `Rbar^-1` makes its native torque
//! work exactly `F dot (R1 r - R0 r)`. A normalized midpoint rotation would
//! break both this identity and the native spatial angular-impulse balance.

use super::{
    DynamicsError, MassProperties, PreparedSpatialLoads, RigidBodyIntegrator, RigidBodyState,
    SpatialStepReceipt, Vec3, Wrench, basis_world, components, finite_derived, finite_vec,
};
use crate::solve_rotation_3;

const BODY_AXES: [Vec3; 3] = [
    Vec3::new(1.0, 0.0, 0.0),
    Vec3::new(0.0, 1.0, 0.0),
    Vec3::new(0.0, 0.0, 1.0),
];

/// Explicit nonlinear and finite-rotation limits for one held-load step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialMidpointConfig {
    /// Maximum Newton corrections, at most 64. Zero permits an exact initial root.
    pub maximum_iterations: usize,
    /// Maximum trials per Newton correction, in `1..=32`.
    pub maximum_backtracks: usize,
    /// Absolute tolerance of the body angular-momentum midpoint equation, N m s.
    pub angular_momentum_tolerance_nms: f64,
    /// Relative tolerance of that equation, finite and in `(0, 0.01]`.
    pub relative_tolerance: f64,
    /// Maximum norm of `dt omega_mid`, finite and in `(0, 2]`. At 2 the
    /// actual Cayley rotation is pi/2; larger steps must be subdivided.
    pub maximum_cayley_increment: f64,
}

/// Refusal publishes no body state or partially accepted work account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpatialMidpointError {
    /// A native state, inertia, force, or finite arithmetic check failed.
    Dynamics(DynamicsError),
    /// An iteration, tolerance or rotation limit was invalid.
    InvalidConfig,
    /// The trial requires a larger Cayley increment than explicitly admitted.
    IncrementLimit,
    /// The bounded solve, or its replay through the original native owner, failed.
    NotConverged {
        /// Number of completed Newton corrections.
        iterations: usize,
        /// Maximum component of the angular-momentum equation residual, N m s.
        residual_nms: f64,
    },
    /// Cancellation was observed before final publication.
    Cancelled,
}

impl From<DynamicsError> for SpatialMidpointError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}

impl core::fmt::Display for SpatialMidpointError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Dynamics(error) => write!(f, "spatial midpoint: {error}"),
            Self::InvalidConfig => f.write_str("invalid spatial midpoint limits"),
            Self::IncrementLimit => f.write_str("spatial midpoint rotation limit; reduce step"),
            Self::NotConverged {
                iterations,
                residual_nms,
            } => write!(
                f,
                "spatial midpoint did not converge after {iterations} corrections ({residual_nms:e} N m s); reduce step"
            ),
            Self::Cancelled => f.write_str("spatial midpoint cancelled"),
        }
    }
}

impl std::error::Error for SpatialMidpointError {}

/// Original native motion plus its work-conjugate load and solve diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialMidpointReceipt {
    /// Native state, exact held-force work, conjugate couple work and balances.
    pub motion: SpatialStepReceipt,
    /// Wrench supplied to the original `step_cayley` owner. The body torque is
    /// the finite-step pullback, not the instantaneous torque at either endpoint.
    pub native_wrench: Wrench,
    /// `dt I^-1 (L0+L1)/2` in the original body axes, dimensionless radians.
    pub cayley_increment_body: Vec3,
    /// Number of accepted Newton corrections.
    pub iterations: usize,
    /// Number of fixed-size midpoint evaluations, including final native replay.
    pub evaluations: usize,
}

impl RigidBodyIntegrator {
    /// Advance held spatial loads with actual material-point work and full inertia.
    ///
    /// A three-variable Newton solve with an analytic Jacobian uses the original
    /// prepared first moments; cost and storage do not grow with the point count.
    /// Its accepted torque is replayed through [`Self::step_cayley`], and the
    /// original returned endpoints must satisfy the same physical load equation.
    /// Every trial is immutable, allocation-free and cancellation-polled.
    ///
    /// Held forces satisfy `delta(K+V_gravity) = force_work + couple_work`
    /// to nonlinear-solve tolerance and roundoff. Free-couple work uses the
    /// conjugate Cayley increment; it is not an exact continuous rotation integral.
    /// Phase and trajectory accuracy remain second order and need refinement.
    /// No symplecticity, global uniqueness, contact law or hidden heat is claimed.
    pub fn step_spatial_midpoint(
        self,
        state: RigidBodyState,
        properties: MassProperties,
        loads: &PreparedSpatialLoads,
        duration_seconds: f64,
        config: SpatialMidpointConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialMidpointReceipt, SpatialMidpointError> {
        poll(&mut cancelled)?;
        config.validate()?;
        if !duration_seconds.is_finite() || duration_seconds <= 0.0 {
            return Err(DynamicsError::InvalidStepDuration.into());
        }
        self.diagnostics(state, properties)?;
        let rotation = state.pose.orientation;
        let system = MidpointSystem {
            properties,
            before: state.angular_momentum_body,
            body_couple: loads.torque_body,
            world_couple_initial: rotation.rotate_world_to_body(loads.torque_world),
            moments_initial: loads
                .point_force_moments_world
                .map(|v| rotation.rotate_world_to_body(v)),
            inverse_inertia: [
                properties.angular_velocity_body_checked(Vec3::new(1.0, 0.0, 0.0))?,
                properties.angular_velocity_body_checked(Vec3::new(0.0, 1.0, 0.0))?,
                properties.angular_velocity_body_checked(Vec3::new(0.0, 0.0, 1.0))?,
            ],
            dt: duration_seconds,
            config,
        };
        let torque0 = loads.wrench_at(rotation)?.torque_body;
        let mut midpoint = system.before.add(torque0.scale(0.5 * system.dt));
        let mut evaluations = 1;
        let mut value = system.evaluate(midpoint)?;
        let mut iterations = 0;
        while !system.converged(midpoint, &value) {
            poll(&mut cancelled)?;
            if iterations == config.maximum_iterations {
                return Err(SpatialMidpointError::NotConverged {
                    iterations,
                    residual_nms: max_component(value.residual),
                });
            }
            let delta = solve_rotation_3(value.jacobian, components(value.residual.scale(-1.0)))?;
            let delta = Vec3::new(delta[0], delta[1], delta[2]);
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..config.maximum_backtracks {
                poll(&mut cancelled)?;
                let trial = midpoint.add(delta.scale(fraction));
                evaluations += 1;
                if let Ok(next) = system.evaluate(trial) {
                    if max_component(next.residual) < max_component(value.residual) {
                        accepted = Some((trial, next));
                        break;
                    }
                }
                fraction *= 0.5;
            }
            (midpoint, value) = accepted.ok_or(SpatialMidpointError::NotConverged {
                iterations,
                residual_nms: max_component(value.residual),
            })?;
            iterations += 1;
        }
        poll(&mut cancelled)?;
        let native_wrench = Wrench {
            force_world: loads.force_world,
            torque_body: value.torque,
        };
        let step = self.step_cayley(state, properties, native_wrench, duration_seconds)?;
        let after = step.state_after;
        let actual_midpoint = system
            .before
            .scale(0.5)
            .add(after.angular_momentum_body.scale(0.5));
        evaluations += 1;
        let actual = system.evaluate(actual_midpoint)?;
        if !system.converged(actual_midpoint, &actual) {
            return Err(SpatialMidpointError::NotConverged {
                iterations,
                residual_nms: max_component(actual.residual),
            });
        }
        let omega_mid = properties.angular_velocity_body_checked(actual_midpoint)?;
        let cayley_increment_body = omega_mid.scale(duration_seconds);
        let cayley_world = rotation.rotate_body_to_world(cayley_increment_body);
        let axes0 = basis_world(rotation);
        let axes1 = basis_world(after.pose.orientation);
        let mut moment_world = loads.torque_world.add(
            rotation
                .rotate_body_to_world(loads.torque_body)
                .scale(0.5)
                .add(
                    after
                        .pose
                        .orientation
                        .rotate_body_to_world(loads.torque_body)
                        .scale(0.5),
                ),
        );
        for j in 0..3 {
            moment_world = moment_world.add(
                axes0[j]
                    .scale(0.5)
                    .add(axes1[j].scale(0.5))
                    .cross(loads.point_force_moments_world[j]),
            );
        }
        let weight = self.gravity.acceleration_world().scale(properties.mass);
        let linear_impulse = loads.force_world.add(weight).scale(duration_seconds);
        let position_mid = state
            .pose
            .position_world
            .scale(0.5)
            .add(after.pose.position_world.scale(0.5));
        let angular_impulse = moment_world
            .scale(duration_seconds)
            .add(position_mid.cross(linear_impulse));
        let linear_residual = after
            .linear_momentum_world
            .sub(state.linear_momentum_world)
            .sub(linear_impulse);
        let angular_residual = step
            .diagnostics_after
            .angular_momentum_world
            .sub(step.diagnostics_before.angular_momentum_world)
            .sub(angular_impulse);
        let force_work_j = loads.force_work(state.pose, after.pose)?;
        let couple_work_j = finite_derived(
            loads.torque_body.dot(cayley_increment_body) + loads.torque_world.dot(cayley_world),
            "spatial_midpoint.couple_work",
        )?;
        // Differential quadratic energy avoids subtracting two large energies.
        let velocity_mid = state
            .center_of_mass_velocity_world(properties)?
            .scale(0.5)
            .add(after.center_of_mass_velocity_world(properties)?.scale(0.5));
        let energy_change = velocity_mid
            .dot(after.linear_momentum_world.sub(state.linear_momentum_world))
            + omega_mid.dot(after.angular_momentum_body.sub(state.angular_momentum_body))
            - weight.dot(after.pose.position_world.sub(state.pose.position_world));
        let work_energy_residual_j = finite_derived(
            energy_change - force_work_j - couple_work_j,
            "spatial_midpoint.work_energy_residual",
        )?;
        for vector in [
            linear_impulse,
            angular_impulse,
            linear_residual,
            angular_residual,
        ] {
            finite_vec(vector, "spatial_midpoint.impulse_balance")?;
        }
        poll(&mut cancelled)?;
        Ok(SpatialMidpointReceipt {
            motion: SpatialStepReceipt {
                step,
                linear_impulse_world: linear_impulse,
                angular_impulse_world: angular_impulse,
                linear_momentum_residual_world: linear_residual,
                angular_momentum_residual_world: angular_residual,
                force_work_j,
                couple_work_j,
                work_energy_residual_j,
            },
            native_wrench,
            cayley_increment_body,
            iterations,
            evaluations,
        })
    }
}

struct MidpointSystem {
    properties: MassProperties,
    before: Vec3,
    body_couple: Vec3,
    world_couple_initial: Vec3,
    moments_initial: [Vec3; 3],
    inverse_inertia: [Vec3; 3],
    dt: f64,
    config: SpatialMidpointConfig,
}

struct Evaluation {
    torque: Vec3,
    residual: Vec3,
    jacobian: [[f64; 3]; 3],
}

impl MidpointSystem {
    fn evaluate(&self, midpoint: Vec3) -> Result<Evaluation, SpatialMidpointError> {
        let omega = self.properties.angular_velocity_body_checked(midpoint)?;
        let c = omega.scale(self.dt);
        if !c.is_finite() || c.norm_squared() > self.config.maximum_cayley_increment.powi(2) {
            return Err(SpatialMidpointError::IncrementLimit);
        }
        let half_c = c.scale(0.5);
        let inverse_det = 1.0 / (1.0 + half_c.norm_squared());
        // B = (I - [c]/2)^-1, hence Rbar = R0 B. No factorization or trig.
        let b = |v: Vec3| {
            v.add(half_c.cross(v))
                .add(half_c.scale(half_c.dot(v)))
                .scale(inverse_det)
        };
        let axes = BODY_AXES.map(b);
        let mut moment = self.world_couple_initial;
        for j in 0..3 {
            moment = moment.add(axes[j].cross(self.moments_initial[j]));
        }
        let torque = self.body_couple.add(moment.sub(half_c.cross(moment)));
        let residual = midpoint
            .sub(self.before)
            .sub(midpoint.cross(omega).add(torque).scale(0.5 * self.dt));
        let mut jacobian = [[0.0; 3]; 3];
        for (j, axis) in BODY_AXES.into_iter().enumerate() {
            let domega = self.inverse_inertia[j];
            let dc = domega.scale(self.dt);
            let mut dmoment = Vec3::ZERO;
            for k in 0..3 {
                dmoment =
                    dmoment.add(b(dc.cross(axes[k]).scale(0.5)).cross(self.moments_initial[k]));
            }
            let dtorque = dmoment
                .sub(half_c.cross(dmoment))
                .sub(dc.cross(moment).scale(0.5));
            let column = axis.sub(
                axis.cross(omega)
                    .add(midpoint.cross(domega))
                    .add(dtorque)
                    .scale(0.5 * self.dt),
            );
            for (i, component) in components(column).into_iter().enumerate() {
                jacobian[i][j] = finite_derived(component, "spatial_midpoint.jacobian")?;
            }
        }
        finite_vec(torque, "spatial_midpoint.torque")?;
        finite_vec(residual, "spatial_midpoint.residual")?;
        Ok(Evaluation {
            torque,
            residual,
            jacobian,
        })
    }

    fn converged(&self, midpoint: Vec3, value: &Evaluation) -> bool {
        let scale = max_component(midpoint)
            .max(max_component(self.before))
            .max(max_component(value.torque.scale(0.5 * self.dt)));
        let tolerance =
            self.config.angular_momentum_tolerance_nms + self.config.relative_tolerance * scale;
        tolerance.is_finite() && max_component(value.residual) <= tolerance
    }
}

impl SpatialMidpointConfig {
    fn validate(self) -> Result<(), SpatialMidpointError> {
        if self.maximum_iterations > 64
            || !(1..=32).contains(&self.maximum_backtracks)
            || !self.angular_momentum_tolerance_nms.is_finite()
            || self.angular_momentum_tolerance_nms < 0.0
            || !self.relative_tolerance.is_finite()
            || self.relative_tolerance <= 0.0
            || self.relative_tolerance > 0.01
            || !self.maximum_cayley_increment.is_finite()
            || self.maximum_cayley_increment <= 0.0
            || self.maximum_cayley_increment > 2.0
        {
            return Err(SpatialMidpointError::InvalidConfig);
        }
        Ok(())
    }
}

fn max_component(v: Vec3) -> f64 {
    v.x.abs().max(v.y.abs()).max(v.z.abs())
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SpatialMidpointError> {
    if cancelled() {
        Err(SpatialMidpointError::Cancelled)
    } else {
        Ok(())
    }
}
