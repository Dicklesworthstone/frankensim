//! Coupled moving dry links and fluid, using the native articulated owner.
//!
//! Every midpoint has one native motion map for dry mass/bias, attached fluid
//! strips and applied forces. The full body/joint mass is solved together;
//! moving dry COMs, link orientation, serial-joint coupling and base recoil are
//! retained. This is a bounded Lie midpoint lane, not the dry owner's separate
//! momentum-coordinate integrator. General nonlinear energy and momentum
//! balance converge with timestep refinement. [`step`] owns smooth motion;
//! [`impact`] supplies simultaneous point-contact resets and [`stops`] supplies
//! mechanical stop events and unilateral support through the same mass owner.
//! Exact finite-step conservation and closed-loop constraints are not claimed.

/// Instantaneous simultaneous normal impacts with full dry and fluid inertia.
#[path = "articulated/impact.rs"]
pub mod impact;

use super::{MAX_CERT_STATES, MemoryState, refuse};
use crate::Refusal;
use crate::addedmass::joint::{admit_fluid, native_index};
use crate::addedmass::{AeroGeneralizedLoads, MAX_CONTROL_COORDS, PreparedArticulatedMass};
use fs_ga::{Se3, Twist, Vec3, Wrench};
use fs_mbd::articulated::generalized::{
    FreeFloatingTerms, GeneralizedError, GeneralizedKinematics, free_floating_terms,
};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{ArticulatedError, ArticulatedModel, BaseState, forward_kinematics};

#[path = "articulated/stops.rs"]
pub mod stops;

/// Sustained Coulomb contact through the original coupled Lie midpoint.
#[path = "articulated/contact.rs"]
pub mod contact;

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
    match error {
        GeneralizedError::Cancelled => cancelled(),
        GeneralizedError::Dynamics(error) => dynamics_error(error),
        other => invalid(other),
    }
}

fn dynamics_error(error: ArticulatedError) -> Refusal {
    if matches!(error, ArticulatedError::JointLimitViolation { .. }) {
        refuse(
            "articulated-joint-position-limit",
            error.to_string(),
            "shorten the trial interval or resolve the mechanical stop before the native limit",
        )
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
    .map_err(dynamics_error)?;
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
    is_cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedStep, Refusal>
where
    F: Fn(f64, &ArticulatedState, &GeneralizedKinematics) -> Result<ArticulatedLoads, Refusal>,
{
    step_with_velocity_response(
        model,
        impedance,
        state,
        memory_inputs,
        t_s,
        dt_s,
        config,
        evaluate,
        |_, _, _| Ok(true),
        is_cancelled,
    )
}

/// One synchronized midpoint view for a simultaneous multi-tree load evaluation.
/// Every state and motion map belongs to the same outer iterate and stage time.
#[derive(Clone, Copy)]
pub struct ArticulatedStage<'a> {
    pub state: &'a ArticulatedState,
    pub kinematics: &'a GeneralizedKinematics,
}

// Contact consumers change only the simultaneous velocity response. Every
// consumer sees the SAME admitted midpoint geometry and complete mass factor.
#[derive(Clone, Copy)]
struct VelocityResponseContext<'a> {
    mass: &'a PreparedArticulatedMass,
    kinematics: &'a GeneralizedKinematics,
    midpoint: &'a ArticulatedState,
}

#[derive(Clone, Copy)]
struct StepSystem<'a> {
    model: &'a ArticulatedModel,
    impedance: &'a [JointImpedance],
    state: &'a ArticulatedState,
    memory_inputs: &'a [f64],
    config: ArticulatedStepConfig<'a>,
}

// The same per-tree stage is used by smooth motion, mechanical stops, and
// single/multiple-tree sustained contact. No consumer owns a second integrator.
struct StepWorkspace<'a> {
    system: StepSystem<'a>,
    initial: [f64; 6 + MAX_CONTROL_COORDS],
    delta: [f64; 6 + MAX_CONTROL_COORDS],
    midpoint_velocity: [f64; 6 + MAX_CONTROL_COORDS],
    midpoint: ArticulatedState,
    control_diagonal: Vec<f64>,
    residual: f64,
}

impl<'a> StepWorkspace<'a> {
    fn new(system: StepSystem<'a>, t_s: f64, dt_s: f64) -> Result<Self, Refusal> {
        let StepSystem {
            model,
            impedance,
            state,
            memory_inputs,
            config,
        } = system;
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
        let half = 0.5 * dt_s;
        let mut initial = [0.0; 6 + MAX_CONTROL_COORDS];
        initial[..6].copy_from_slice(&state.motion.base.twist_body.to_array());
        initial[6..6 + nc].copy_from_slice(&state.motion.joint_velocity);
        let control_diagonal: Vec<_> = impedance
            .iter()
            .map(|p| half * p.damping + half * half * p.stiffness)
            .collect();
        if !control_diagonal.iter().all(|x| x.is_finite()) {
            return Err(invalid("unrepresentable implicit control diagonal"));
        }
        let mut midpoint = state.clone();
        midpoint.memory = state
            .memory
            .iter()
            .zip(memory_inputs)
            .map(|(m, u)| m.advanced(*u, half))
            .collect();
        if !midpoint.memory.iter().all(|m| m.value.is_finite()) {
            return Err(invalid("unrepresentable midpoint memory state"));
        }
        Ok(Self {
            system,
            initial,
            delta: [0.0; 6 + MAX_CONTROL_COORDS],
            midpoint_velocity: initial,
            midpoint,
            control_diagonal,
            residual: f64::INFINITY,
        })
    }

    fn prepare(
        &mut self,
        half: f64,
        is_cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FreeFloatingTerms, Refusal> {
        let nc = self.system.model.dof_count();
        let mut um = self.initial;
        for (i, value) in um.iter_mut().enumerate().take(6 + nc) {
            *value += 0.5 * self.delta[i];
        }
        self.midpoint_velocity = um;
        self.midpoint.motion.base.twist_body = twist(&um);
        self.midpoint.motion.base.world_from_base = pose_at(
            self.system.state.motion.base.world_from_base,
            twist(&um),
            half,
        )?;
        for k in 0..nc {
            self.midpoint.motion.joint_velocity[k] = um[6 + k];
            self.midpoint.motion.joint_position[k] =
                self.system.state.motion.joint_position[k] + half * um[6 + k];
        }
        free_floating_terms(
            self.system.model,
            self.midpoint.motion.base,
            &self.midpoint.motion.joint_position,
            &self.midpoint.motion.joint_velocity,
            is_cancelled,
        )
        .map_err(native_error)
    }

    fn free_response(
        &self,
        dry: &FreeFloatingTerms,
        loads: ArticulatedLoads,
        dt_s: f64,
    ) -> Result<(PreparedArticulatedMass, Vec<f64>), Refusal> {
        let StepSystem {
            model,
            impedance,
            state,
            ..
        } = self.system;
        let nc = model.dof_count();
        let half = 0.5 * dt_s;
        let initial = &self.initial;
        let um = &self.midpoint_velocity;
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
                    - p.stiffness * (self.midpoint.motion.joint_position[k] - p.rest_position)
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
        let mass =
            PreparedArticulatedMass::with_control_diagonal(dry, &fluid, &self.control_diagonal)?;
        let rhs: Vec<_> = (0..6 + nc)
            .map(|i| {
                let applied = if i < 6 {
                    fluid.q_rigid_nonaccel[i]
                } else {
                    fluid.q_control_nonaccel[i - 6]
                };
                applied + fluid.q_added_bias[i] - dry.bias[native_index(i)]
            })
            .collect();
        let acceleration = mass.solve(&rhs)?;
        let increment = acceleration.iter().map(|a| dt_s * a).collect();
        Ok((mass, increment))
    }

    fn next(&mut self, increment: &[f64]) -> Result<[f64; 6 + MAX_CONTROL_COORDS], Refusal> {
        let nc = self.system.model.dof_count();
        let n = 6 + nc;
        let config = self.system.config;
        let mut next = [0.0; 6 + MAX_CONTROL_COORDS];
        for i in 0..6 {
            next[native_index(i)] = increment[i];
        }
        next[6..n].copy_from_slice(&increment[6..n]);
        if !next[..n].iter().all(|x| x.is_finite()) {
            return Err(invalid("unrepresentable coupled velocity increment"));
        }
        self.residual = 0.0_f64;
        for i in 0..n {
            let budget = if i < 3 {
                config.angular_tolerance_rad_s
            } else if i < 6 {
                config.linear_tolerance_mps
            } else {
                config.joint_velocity_tolerances[i - 6]
            };
            self.residual = self.residual.max((next[i] - self.delta[i]).abs() / budget);
        }
        Ok(next)
    }

    fn accept(self, dt_s: f64, iteration: usize) -> Result<ArticulatedStep, Refusal> {
        let nc = self.system.model.dof_count();
        let mut accepted = self.midpoint;
        let mut initial = self.initial;
        let um = self.midpoint_velocity;
        // Publish the candidate whose midpoint was actually evaluated.
        accepted.motion.base.world_from_base = pose_at(
            self.system.state.motion.base.world_from_base,
            twist(&um),
            dt_s,
        )?;
        for i in 0..6 {
            initial[i] += self.delta[i];
        }
        accepted.motion.base.twist_body = twist(&initial);
        for k in 0..nc {
            accepted.motion.joint_position[k] =
                self.system.state.motion.joint_position[k] + dt_s * um[6 + k];
            accepted.motion.joint_velocity[k] = initial[6 + k] + self.delta[6 + k];
        }
        accepted.memory = self
            .system
            .state
            .memory
            .iter()
            .zip(self.system.memory_inputs)
            .map(|(m, u)| m.advanced(*u, dt_s))
            .collect();
        if !accepted.memory.iter().all(|m| m.value.is_finite()) {
            return Err(invalid("unrepresentable accepted memory state"));
        }
        admit_state(self.system.model, &accepted.motion)?;
        Ok(ArticulatedStep {
            state: accepted,
            iterations: iteration,
            normalized_residual: self.residual,
        })
    }
}

fn step_systems_with_velocity_response<F, R>(
    systems: &[StepSystem<'_>],
    t_s: f64,
    dt_s: f64,
    evaluate: F,
    mut response: R,
    mut is_cancelled: impl FnMut() -> bool,
) -> Result<Vec<ArticulatedStep>, Refusal>
where
    F: Fn(f64, &[ArticulatedStage<'_>]) -> Result<Vec<ArticulatedLoads>, Refusal>,
    R: FnMut(
        &[VelocityResponseContext<'_>],
        &mut [Vec<f64>],
        &mut dyn FnMut() -> bool,
    ) -> Result<bool, Refusal>,
{
    if is_cancelled() {
        return Err(cancelled());
    }
    if systems.is_empty()
        || systems.len() > 128
        || systems
            .iter()
            .try_fold(0usize, |n, s| n.checked_add(s.model.link_count()))
            .is_none_or(|n| n > 4096)
    {
        return Err(invalid(
            "one to 128 trees and at most 4096 total links are required",
        ));
    }
    let mut work: Vec<_> = systems
        .iter()
        .copied()
        .map(|system| StepWorkspace::new(system, t_s, dt_s))
        .collect::<Result<_, _>>()?;
    let maximum_iterations = systems
        .iter()
        .map(|s| s.config.max_iterations)
        .min()
        .unwrap_or(1);
    let half = 0.5 * dt_s;
    let mut residual = f64::INFINITY;
    for iteration in 1..=maximum_iterations {
        // Prepare every candidate before evaluating any cross-system loads.
        let dry: Vec<_> = work
            .iter_mut()
            .map(|w| w.prepare(half, &mut is_cancelled))
            .collect::<Result<_, _>>()?;
        let views: Vec<_> = work
            .iter()
            .zip(&dry)
            .map(|(w, d)| ArticulatedStage {
                state: &w.midpoint,
                kinematics: &d.kinematics,
            })
            .collect();
        let loads = evaluate(t_s + half, &views)?;
        if is_cancelled() {
            return Err(cancelled());
        }
        if loads.len() != systems.len() {
            return Err(invalid(
                "one synchronized load evaluation per participating tree is required",
            ));
        }
        let evaluated: Vec<_> = work
            .iter()
            .zip(&dry)
            .zip(loads)
            .map(|((w, d), loads)| w.free_response(d, loads, dt_s))
            .collect::<Result<_, _>>()?;
        let (masses, mut increments): (Vec<_>, Vec<_>) = evaluated.into_iter().unzip();
        let contexts: Vec<_> = work
            .iter()
            .zip(&dry)
            .zip(&masses)
            .map(|((w, d), mass)| VelocityResponseContext {
                mass,
                kinematics: &d.kinematics,
                midpoint: &w.midpoint,
            })
            .collect();
        let ready = response(&contexts, &mut increments, &mut is_cancelled)?;
        let next: Vec<_> = work
            .iter_mut()
            .zip(&increments)
            .map(|(w, increment)| w.next(increment))
            .collect::<Result<_, _>>()?;
        residual = work.iter().map(|w| w.residual).fold(0.0_f64, f64::max);
        if residual <= 1.0 && ready {
            let steps = work
                .into_iter()
                .map(|w| w.accept(dt_s, iteration))
                .collect::<Result<_, _>>()?;
            if is_cancelled() {
                return Err(cancelled());
            }
            return Ok(steps);
        }
        for (w, delta) in work.iter_mut().zip(next) {
            w.delta = delta;
        }
    }
    Err(refuse(
        "articulated-step-not-converged",
        format!("midpoint velocity residual/budget {residual:e}"),
        "reduce the timestep or increase the bounded midpoint iteration budget",
    ))
}

#[allow(clippy::too_many_arguments)]
fn step_with_velocity_response<F, R>(
    model: &ArticulatedModel,
    impedance: &[JointImpedance],
    state: &ArticulatedState,
    memory_inputs: &[f64],
    t_s: f64,
    dt_s: f64,
    config: ArticulatedStepConfig<'_>,
    evaluate: F,
    mut response: R,
    is_cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedStep, Refusal>
where
    F: Fn(f64, &ArticulatedState, &GeneralizedKinematics) -> Result<ArticulatedLoads, Refusal>,
    R: FnMut(
        VelocityResponseContext<'_>,
        &mut [f64],
        &mut dyn FnMut() -> bool,
    ) -> Result<bool, Refusal>,
{
    let systems = [StepSystem {
        model,
        impedance,
        state,
        memory_inputs,
        config,
    }];
    step_systems_with_velocity_response(
        &systems,
        t_s,
        dt_s,
        |time, stages| Ok(vec![evaluate(time, stages[0].state, stages[0].kinematics)?]),
        |contexts, increments, cancel| response(contexts[0], &mut increments[0], cancel),
        is_cancelled,
    )?
    .into_iter()
    .next()
    .ok_or_else(|| invalid("missing accepted tree"))
}
