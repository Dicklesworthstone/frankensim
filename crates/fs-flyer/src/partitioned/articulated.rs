//! Coupled moving dry links and fluid, using the native articulated owner.
//!
//! Every midpoint has one native motion map for dry mass/bias, attached fluid
//! strips and applied forces. The full body/joint mass is solved together;
//! moving dry COMs, link orientation, serial-joint coupling and base recoil are
//! retained. This is a bounded Lie midpoint lane, not the dry owner's separate
//! momentum-coordinate integrator. General nonlinear energy and momentum
//! balance converge with timestep refinement; exact finite-step preservation,
//! impacts, joint stops and closed-loop constraints are not claimed.

use super::{MAX_CERT_STATES, MemoryState, refuse};
use crate::Refusal;
use crate::addedmass::joint::{admit_fluid, native_index, solve_articulated_with_diagonal};
use crate::addedmass::{AeroGeneralizedLoads, MAX_CONTROL_COORDS};
use fs_ga::{Se3, Twist, Vec3, Wrench};
use fs_mbd::articulated::generalized::{
    GeneralizedError, GeneralizedKinematics, free_floating_terms,
};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{ArticulatedModel, BaseState, forward_kinematics};

/// Link cap for the compact aircraft coupling lane, including fixed links.
pub const MAX_ARTICULATED_LINKS: usize = 256;

/// Complete native free-base state and the existing linear aero-memory states.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedState {
    /// Native base pose/body twist and compact joint positions/rates.
    pub motion: FreeFloatingMotionState,
    /// Linear states whose declared inputs are held over the tick.
    pub memory: Vec<MemoryState>,
}

/// Passive linear impedance of one scalar native joint, without extra inertia.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointImpedance {
    /// Nonnegative effort/rate coefficient: N m s/rad or N s/m.
    pub damping: f64,
    /// Nonnegative effort/displacement coefficient: N m/rad or N/m.
    pub stiffness: f64,
    /// Unloaded joint coordinate, in radians or metres.
    pub rest_position: f64,
}

/// Physical residual budgets for the mixed rotational/translational unknowns.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArticulatedStepConfig<'a> {
    /// Positive solve budget, no greater than `coupled::MAX_COUPLED_ITERATIONS`.
    pub max_iterations: usize,
    /// Absolute base linear-velocity residual [m/s].
    pub linear_tolerance_mps: f64,
    /// Absolute base angular-velocity residual [rad/s].
    pub angular_tolerance_rad_s: f64,
    /// One positive rate tolerance per joint, in its native rad/s or m/s units.
    pub joint_velocity_tolerances: &'a [f64],
}

/// A pure midpoint load evaluation on the shared native motion map.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedLoads {
    /// Added mass and complete fluid bias at this same stage. Its applied-load
    /// fields must be zero; the fields below own all noninertial loads.
    pub fluid: AeroGeneralizedLoads,
    /// One external wrench per link, in world axes, with torque about THAT
    /// LINK'S ORIGIN. COM-applied forces must include their origin moment.
    pub world_wrenches: Vec<Wrench>,
    /// Explicit uniform gravitational acceleration applied to every dry COM.
    /// Fluid hydrostatic/buoyancy forces, if desired, belong in world_wrenches.
    pub gravity_world: Vec3,
    /// Internal joint efforts excluding the declared linear spring/damper.
    pub joint_effort: Vec<f64>,
}

/// Accepted coupled tick; no borrowed input state is mutated on refusal.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedStep {
    /// Full next state, including exact full-tick memory transitions.
    pub state: ArticulatedState,
    /// Number of simultaneous mass solves used.
    pub iterations: usize,
    /// Largest velocity residual divided by its declared physical budget.
    pub normalized_residual: f64,
}

fn invalid(message: impl ToString) -> Refusal {
    refuse(
        "articulated-step-invalid",
        message.to_string(),
        "inspect the admitted model, loads, state and physical step budgets",
    )
}

fn cancelled() -> Refusal {
    refuse(
        "articulated-step-cancelled",
        "coupled step cancelled before publication".into(),
        "resume from the unchanged input state",
    )
}

fn native_error(error: GeneralizedError) -> Refusal {
    if error == GeneralizedError::Cancelled {
        cancelled()
    } else {
        invalid(error)
    }
}

fn twist(u: &[f64]) -> Twist {
    Twist::new(Vec3::new(u[0], u[1], u[2]), Vec3::new(u[3], u[4], u[5]))
}

fn pose_at(pose: Se3, velocity: Twist, dt: f64) -> Result<Se3, Refusal> {
    pose.compose(Se3::exp(velocity.scale(dt)).map_err(invalid)?)
        .map_err(invalid)
}

fn admit_state(model: &ArticulatedModel, state: &FreeFloatingMotionState) -> Result<(), Refusal> {
    model.free_floating_complexity().map_err(invalid)?;
    forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )
    .map_err(invalid)?;
    Ok(())
}

/// Advance an articulated body and its attached fluid through one joint solve.
///
/// The callback receives the actual midpoint motion and the SAME native link
/// kinematics used for dry inertia. Use `addedmass::articulated` to assemble
/// serial moving strips from this map. The callback must be a pure evaluation;
/// nonlinear iteration may repeat the same time. Memory inputs are held, with
/// exact half-tick values at every evaluation and a full transition on acceptance.
///
/// Linear impedance is eliminated as `h*c + h*h*k` on each joint diagonal,
/// where `h=dt/2`; no dry joint inertia or body recoil is approximated by a
/// diagonal. Native joint position/rate and total internal-effort limits apply.
/// Remaining nonlinear terms require a converged fixed point within the budget.
///
/// # Errors
/// Oversized/non-free models, invalid states/loads/units/budgets, mismatched
/// blocks, duplicated applied fluid loads, unresolved mass, joint-limit
/// violations, iteration exhaustion or cancellation. Inputs remain unchanged.
#[allow(clippy::too_many_arguments)]
pub fn step<F>(
    model: &ArticulatedModel,
    impedance: &[JointImpedance],
    state: &ArticulatedState,
    memory_inputs: &[f64],
    t_s: f64,
    dt_s: f64,
    config: ArticulatedStepConfig<'_>,
    evaluate: F,
    mut is_cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedStep, Refusal>
where
    F: Fn(f64, &ArticulatedState, &GeneralizedKinematics) -> Result<ArticulatedLoads, Refusal>,
{
    if is_cancelled() {
        return Err(cancelled());
    }
    let nc = model.dof_count();
    if nc > MAX_CONTROL_COORDS
        || model.link_count() > MAX_ARTICULATED_LINKS
        || impedance.len() != nc
        || config.joint_velocity_tolerances.len() != nc
        || config.max_iterations == 0
        || config.max_iterations > super::coupled::MAX_COUPLED_ITERATIONS
        || !t_s.is_finite()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !(t_s + dt_s).is_finite()
        || t_s + dt_s <= t_s
        || !config
            .joint_velocity_tolerances
            .iter()
            .chain([
                &config.linear_tolerance_mps,
                &config.angular_tolerance_rad_s,
            ])
            .all(|x| x.is_finite() && *x > 0.0)
        || !impedance.iter().all(|p| {
            p.damping.is_finite()
                && p.damping >= 0.0
                && p.stiffness.is_finite()
                && p.stiffness >= 0.0
                && p.rest_position.is_finite()
        })
    {
        return Err(invalid(
            "invalid coupled model size, impedance, time or residual budget",
        ));
    }
    if memory_inputs.len() != state.memory.len()
        || state.memory.len() > MAX_CERT_STATES
        || !memory_inputs.iter().all(|x| x.is_finite())
        || !state
            .memory
            .iter()
            .all(|m| m.tau_s.is_finite() && m.tau_s > 0.0 && m.value.is_finite())
    {
        return Err(invalid(
            "one finite held input per finite positive-time-scale memory state is required",
        ));
    }
    admit_state(model, &state.motion)?;
    let n = 6 + nc;
    let half = 0.5 * dt_s;
    let initial_pose = state.motion.base.world_from_base;
    let mut initial = [0.0; 6 + MAX_CONTROL_COORDS];
    initial[..6].copy_from_slice(&state.motion.base.twist_body.to_array());
    initial[6..n].copy_from_slice(&state.motion.joint_velocity);
    let control_diagonal: Vec<_> = impedance
        .iter()
        .map(|p| half * p.damping + half * half * p.stiffness)
        .collect();
    if !control_diagonal.iter().all(|x| x.is_finite()) {
        return Err(invalid("unrepresentable implicit control diagonal"));
    }
    let memory_at = |duration| {
        state
            .memory
            .iter()
            .zip(memory_inputs)
            .map(|(m, u)| m.advanced(*u, duration))
            .collect::<Vec<_>>()
    };
    let mut midpoint = state.clone();
    midpoint.memory = memory_at(half);
    if !midpoint.memory.iter().all(|m| m.value.is_finite()) {
        return Err(invalid("unrepresentable midpoint memory state"));
    }
    let mut delta = [0.0; 6 + MAX_CONTROL_COORDS];
    let mut residual = f64::INFINITY;
    for iteration in 1..=config.max_iterations {
        let mut um = initial;
        for i in 0..n {
            um[i] += 0.5 * delta[i];
        }
        midpoint.motion.base.twist_body = twist(&um);
        midpoint.motion.base.world_from_base = pose_at(initial_pose, twist(&um), half)?;
        for k in 0..nc {
            midpoint.motion.joint_velocity[k] = um[6 + k];
            midpoint.motion.joint_position[k] = state.motion.joint_position[k] + half * um[6 + k];
        }
        let dry = free_floating_terms(
            model,
            midpoint.motion.base,
            &midpoint.motion.joint_position,
            &midpoint.motion.joint_velocity,
            &mut is_cancelled,
        )
        .map_err(native_error)?;
        let loads = evaluate(t_s + half, &midpoint, &dry.kinematics)?;
        if is_cancelled() {
            return Err(cancelled());
        }
        admit_fluid(&loads.fluid)?;
        if loads.fluid.q_control_nonaccel.len() != nc
            || loads.joint_effort.len() != nc
            || !loads.joint_effort.iter().all(|x| x.is_finite())
        {
            return Err(invalid(
                "applied loads do not match the native active joint set",
            ));
        }
        if loads
            .fluid
            .q_rigid_nonaccel
            .iter()
            .chain(&loads.fluid.q_control_nonaccel)
            .any(|x| *x != 0.0)
        {
            return Err(refuse(
                "coupled-load-ownership",
                "fluid assembly contains duplicated applied loads".into(),
                "supply applied loads only through world_wrenches, gravity_world and joint_effort",
            ));
        }
        let mut applied = dry
            .kinematics
            .project_world_wrenches(&loads.world_wrenches)
            .map_err(native_error)?;
        let gravity = dry
            .kinematics
            .project_uniform_gravity(loads.gravity_world)
            .map_err(native_error)?;
        for (force, gravity) in applied.iter_mut().zip(gravity) {
            *force += gravity;
        }
        let internal_effort: Vec<_> = impedance
            .iter()
            .enumerate()
            .map(|(k, p)| {
                loads.joint_effort[k]
                    - p.damping * um[6 + k]
                    - p.stiffness * (midpoint.motion.joint_position[k] - p.rest_position)
            })
            .collect();
        model
            .validate_applied_effort(&internal_effort)
            .map_err(invalid)?;
        let mut fluid = loads.fluid;
        for i in 0..6 {
            fluid.q_rigid_nonaccel[i] = applied[native_index(i)];
        }
        for k in 0..nc {
            let p = impedance[k];
            fluid.q_control_nonaccel[k] = applied[6 + k] + loads.joint_effort[k]
                - p.damping * initial[6 + k]
                - p.stiffness
                    * (state.motion.joint_position[k] - p.rest_position + half * initial[6 + k]);
        }
        let (base_accel, joint_accel) =
            solve_articulated_with_diagonal(&dry, &fluid, &control_diagonal)?;
        let mut next = [0.0; 6 + MAX_CONTROL_COORDS];
        for i in 0..6 {
            next[native_index(i)] = dt_s * base_accel[i];
        }
        for k in 0..nc {
            next[6 + k] = dt_s * joint_accel[k];
        }
        if !next[..n].iter().all(|x| x.is_finite()) {
            return Err(invalid("unrepresentable coupled velocity increment"));
        }
        residual = 0.0_f64;
        for i in 0..n {
            let budget = if i < 3 {
                config.angular_tolerance_rad_s
            } else if i < 6 {
                config.linear_tolerance_mps
            } else {
                config.joint_velocity_tolerances[i - 6]
            };
            residual = residual.max((next[i] - delta[i]).abs() / budget);
        }
        if residual <= 1.0 {
            // Publish the candidate whose midpoint was actually evaluated.
            let mut accepted = midpoint;
            accepted.motion.base.world_from_base = pose_at(initial_pose, twist(&um), dt_s)?;
            for i in 0..6 {
                initial[i] += delta[i];
            }
            accepted.motion.base.twist_body = twist(&initial);
            for k in 0..nc {
                accepted.motion.joint_position[k] =
                    state.motion.joint_position[k] + dt_s * um[6 + k];
                accepted.motion.joint_velocity[k] = initial[6 + k] + delta[6 + k];
            }
            accepted.memory = memory_at(dt_s);
            if !accepted.memory.iter().all(|m| m.value.is_finite()) {
                return Err(invalid("unrepresentable accepted memory state"));
            }
            admit_state(model, &accepted.motion)?;
            if is_cancelled() {
                return Err(cancelled());
            }
            return Ok(ArticulatedStep {
                state: accepted,
                iterations: iteration,
                normalized_residual: residual,
            });
        }
        delta = next;
    }
    Err(refuse(
        "articulated-step-not-converged",
        format!("midpoint velocity residual/budget {residual:e}"),
        "reduce the timestep or increase the bounded midpoint iteration budget",
    ))
}
