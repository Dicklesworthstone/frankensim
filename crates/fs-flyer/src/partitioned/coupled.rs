//! Simultaneous body/hinge midpoint motion through the original joint mass solve.
//!
//! The body velocities and hinge rate are one seven-coordinate unknown. Fluid
//! geometry, momentum bias and applied loads are reevaluated at its midpoint;
//! cross-block accelerations are never lagged into a separate rigid-body step.
//! Linear hinge stiffness/damping enter the joint matrix analytically, so their
//! stiffness alone does not limit the fixed-point iteration. Nonlinear fluid,
//! gyroscopic and load coupling still requires convergence within the budget.
//!
//! This advances the existing constant dry-body + scalar dry-hinge model. A
//! mechanical linkage whose motion changes the dry mass matrix belongs to the
//! articulated owner; it cannot be represented by an extra diagonal inertia.
//! Lie reconstruction stays on SO(3). General nonlinear energy/world-momentum
//! conservation is a time-refinement property, not a finite-step exact claim.

use super::{HingeParams, MemoryState, PartitionedState, refuse};
use crate::Refusal;
use crate::addedmass::{AeroGeneralizedLoads, solve_joint};
use crate::spine::{Loads, RigidBody, admit_and_canonicalize_state};
use fs_ga::{GaError, Quat, So3, So3Tangent, Vec3};

/// Maximum midpoint iterations admitted for one tick.
pub const MAX_COUPLED_ITERATIONS: usize = 64;

/// Absolute velocity residual budgets in their physical units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledStepConfig {
    /// Positive iteration budget, at most [`MAX_COUPLED_ITERATIONS`].
    pub max_iterations: usize,
    /// Body linear-velocity residual [m/s].
    pub linear_tolerance_mps: f64,
    /// Body angular-velocity residual [rad/s].
    pub angular_tolerance_rad_s: f64,
    /// Hinge-rate residual [rad/s].
    pub hinge_tolerance_rad_s: f64,
}

/// One midpoint evaluation with explicit load ownership.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledLoads {
    /// Fluid mass blocks and their complete momentum/configuration bias.
    /// Its nonacceleration-load fields must be zero: applied loads enter below.
    pub fluid: AeroGeneralizedLoads,
    /// Applied world force and body moment, excluding all inertial bias.
    pub rigid: Loads,
    /// Applied generalized hinge torque [N m], excluding the spring/damper.
    pub hinge_torque_nm: f64,
}

/// An accepted, fully coupled tick; refusal leaves the borrowed state untouched.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledStep {
    /// Complete next state, including the original exact memory transitions.
    pub state: PartitionedState,
    /// Number of joint solves used.
    pub iterations: usize,
    /// Largest velocity residual divided by its corresponding declared budget.
    pub normalized_residual: f64,
}

fn geometry(error: GaError) -> Refusal {
    refuse(
        "coupled-state-invalid",
        error.to_string(),
        "reduce the timestep or inspect the state",
    )
}

fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn orientation(rotation: So3, omega: Vec3, dt: f64) -> Result<So3, Refusal> {
    rotation
        .compose(So3::exp(So3Tangent::new(omega.scale(dt))).map_err(geometry)?)
        .map_err(geometry)
}

fn quaternion(rotation: So3) -> [f64; 4] {
    let q = rotation.as_quat();
    [q.w, q.x, q.y, q.z]
}

fn admit_fluid(loads: &CoupledLoads) -> Result<(), Refusal> {
    let f = &loads.fluid;
    if f.q_control_nonaccel.len() != 1
        || f.m_added_rc.len() != 1
        || f.m_added_cc.len() != 1
        || f.m_added_cc[0].len() != 1
        || f.q_added_bias.len() != 7
    {
        return Err(refuse(
            "coupled-load-shape",
            "one active hinge is required".into(),
            "assemble one control coordinate, including a seven-entry bias",
        ));
    }
    if f.q_rigid_nonaccel
        .iter()
        .chain(&f.q_control_nonaccel)
        .any(|x| *x != 0.0)
    {
        return Err(refuse(
            "coupled-load-ownership",
            "fluid assembly contains applied loads".into(),
            "supply applied loads through rigid and hinge_torque_nm only",
        ));
    }
    if !f
        .m_added_rr
        .iter()
        .flatten()
        .chain(f.m_added_rc.iter().flatten())
        .chain(f.m_added_cc.iter().flatten())
        .chain(&f.q_added_bias)
        .chain(&loads.rigid.force_n)
        .chain(&loads.rigid.moment_nm)
        .chain(core::iter::once(&loads.hinge_torque_nm))
        .all(|x| x.is_finite())
    {
        return Err(refuse(
            "coupled-load-invalid",
            "nonfinite mass, bias or load".into(),
            "use finite admitted inertia and loads",
        ));
    }
    for i in 0..6 {
        for j in 0..i {
            let a = f.m_added_rr[i][j];
            let b = f.m_added_rr[j][i];
            if (a - b).abs() > 64.0 * f64::EPSILON * a.abs().max(b.abs()) {
                return Err(refuse(
                    "coupled-load-invalid",
                    "nonsymmetric fluid mass".into(),
                    "use the original symmetric added-mass assembly",
                ));
            }
        }
    }
    Ok(())
}

/// Advance the existing partitioned state with a simultaneous body/hinge solve.
///
/// `hinge.inertia` is DRY structural inertia here; added hinge inertia and
/// cross blocks come only from the assembly callback. `memory_inputs` are held
/// over the tick. The callback observes their exact half-tick state, the actual
/// midpoint pose/hinge state, and `nu = [v_body, omega_body]`. It must be a pure
/// evaluation: nonlinear iterations may call it more than once at the same time.
/// [`crate::addedmass::hinged::assemble_hinged_strips`] supplies finite-deflection
/// fluid inertia and its analytic bias directly at this boundary.
///
/// Linear spring/damping are eliminated into `I + c*dt/2 + k*dt²/4` within
/// the SAME joint solve. The remaining gyroscopic/configuration/load fixed
/// point refuses when it cannot meet the velocity budgets. Only converged
/// candidates are returned; no partial memory or hinge advancement is published.
///
/// # Errors
/// Invalid states, coefficients, budgets, memory/assembly dimensions or load
/// ownership; non-PD effective mass; nonlinear iteration exhaustion.
#[allow(clippy::too_many_arguments)]
pub fn step<F>(
    body: &RigidBody,
    hinge: &HingeParams,
    state: &PartitionedState,
    memory_inputs: &[f64],
    t_s: f64,
    dt_s: f64,
    config: CoupledStepConfig,
    evaluate: F,
) -> Result<CoupledStep, Refusal>
where
    F: Fn(f64, &PartitionedState, &[f64; 6]) -> Result<CoupledLoads, Refusal>,
{
    let rigid = admit_and_canonicalize_state(body, &state.rigid, dt_s)?;
    if !t_s.is_finite()
        || !(t_s + dt_s).is_finite()
        || t_s + dt_s <= t_s
        || !state.hinge_q.is_finite()
        || !state.hinge_v.is_finite()
        || !hinge.inertia.is_finite()
        || hinge.inertia <= 0.0
        || !hinge.damping.is_finite()
        || hinge.damping < 0.0
        || !hinge.stiffness.is_finite()
        || hinge.stiffness < 0.0
        || config.max_iterations == 0
        || config.max_iterations > MAX_COUPLED_ITERATIONS
        || ![
            config.linear_tolerance_mps,
            config.angular_tolerance_rad_s,
            config.hinge_tolerance_rad_s,
        ]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
    {
        return Err(refuse(
            "coupled-input-invalid",
            "invalid state, time, coefficients or budget".into(),
            "use finite state, positive inertia/tolerances and an admitted iteration budget",
        ));
    }
    if memory_inputs.len() != state.memory.len()
        || state.memory.len() > super::MAX_CERT_STATES
        || !memory_inputs.iter().all(|x| x.is_finite())
        || !state
            .memory
            .iter()
            .all(|m| m.tau_s.is_finite() && m.tau_s > 0.0 && m.value.is_finite())
    {
        return Err(refuse(
            "memory-input-mismatch",
            "invalid memory states or inputs".into(),
            "supply one finite held input per finite positive-time-scale memory state",
        ));
    }
    let q = rigid.quat;
    let rotation = So3::try_from_quat(Quat {
        w: q[0],
        x: q[1],
        y: q[2],
        z: q[3],
    })
    .map_err(geometry)?;
    let v0 = array(
        rotation
            .inverse()
            .rotate(vector(rigid.vel_mps))
            .map_err(geometry)?,
    );
    let mut initial = [0.0; 7];
    initial[..3].copy_from_slice(&v0);
    initial[3..6].copy_from_slice(&rigid.omega_body);
    initial[6] = state.hinge_v;
    let half = 0.5 * dt_s;
    let effective_hinge = hinge.inertia + half * hinge.damping + half * half * hinge.stiffness;
    let memory_at = |dt| -> Vec<MemoryState> {
        state
            .memory
            .iter()
            .zip(memory_inputs)
            .map(|(m, u)| m.advanced(*u, dt))
            .collect()
    };
    let mut midpoint = PartitionedState {
        rigid,
        hinge_q: state.hinge_q,
        hinge_v: state.hinge_v,
        memory: memory_at(half),
    };
    let mut delta = [0.0; 7];
    let mut residual = f64::INFINITY;
    for iteration in 1..=config.max_iterations {
        let mut nu = [0.0; 6];
        for i in 0..6 {
            nu[i] = initial[i] + 0.5 * delta[i];
        }
        midpoint.hinge_v = initial[6] + 0.5 * delta[6];
        midpoint.hinge_q = state.hinge_q + half * midpoint.hinge_v;
        let omega = Vec3::new(nu[3], nu[4], nu[5]);
        let vbody = Vec3::new(nu[0], nu[1], nu[2]);
        let rmid = orientation(rotation, omega, half)?;
        let vworld = rmid.rotate(vbody).map_err(geometry)?;
        midpoint.rigid.quat = quaternion(rmid);
        midpoint.rigid.vel_mps = array(vworld);
        midpoint.rigid.omega_body = array(omega);
        midpoint.rigid.pos_m = array(vector(rigid.pos_m) + vworld.scale(half));
        admit_and_canonicalize_state(body, &midpoint.rigid, dt_s)?;
        if !midpoint.hinge_q.is_finite() || !midpoint.hinge_v.is_finite() {
            return Err(refuse(
                "coupled-state-invalid",
                "nonfinite midpoint hinge".into(),
                "reduce the step",
            ));
        }
        let stage = evaluate(t_s + half, &midpoint, &nu)?;
        admit_fluid(&stage)?;
        let mut fluid = stage.fluid;
        let force = rmid
            .inverse()
            .rotate(vector(stage.rigid.force_n))
            .map_err(geometry)?
            - omega.cross(vbody.scale(body.mass_kg));
        let spin = body.inertia_tensor().matrix().apply(omega);
        let torque = vector(stage.rigid.moment_nm) - omega.cross(spin);
        fluid.q_rigid_nonaccel[..3].copy_from_slice(&array(force));
        fluid.q_rigid_nonaccel[3..].copy_from_slice(&array(torque));
        fluid.q_control_nonaccel[0] = stage.hinge_torque_nm
            - hinge.damping * state.hinge_v
            - hinge.stiffness * (state.hinge_q + half * state.hinge_v);
        let (acceleration, control) = solve_joint(body, &[effective_hinge], &fluid)?;
        let mut next = [0.0; 7];
        for i in 0..6 {
            next[i] = dt_s * acceleration[i];
        }
        next[6] = dt_s * control[0];
        if next.iter().any(|x| !x.is_finite()) {
            return Err(refuse(
                "coupled-state-invalid",
                "nonfinite coupled increment".into(),
                "reduce the step",
            ));
        }
        residual = 0.0_f64;
        for i in 0..7 {
            let budget = if i < 3 {
                config.linear_tolerance_mps
            } else if i < 6 {
                config.angular_tolerance_rad_s
            } else {
                config.hinge_tolerance_rad_s
            };
            residual = residual.max((next[i] - delta[i]).abs() / budget);
        }
        if residual <= 1.0 {
            let rnext = orientation(rotation, omega, dt_s)?;
            let mut result = midpoint;
            result.rigid.pos_m = array(vector(rigid.pos_m) + vworld.scale(dt_s));
            result.rigid.vel_mps = array(
                rnext
                    .rotate(Vec3::new(
                        initial[0] + delta[0],
                        initial[1] + delta[1],
                        initial[2] + delta[2],
                    ))
                    .map_err(geometry)?,
            );
            result.rigid.omega_body = [
                initial[3] + delta[3],
                initial[4] + delta[4],
                initial[5] + delta[5],
            ];
            result.rigid.quat = quaternion(rnext);
            result.hinge_q = state.hinge_q + dt_s * result.hinge_v;
            result.hinge_v = initial[6] + delta[6];
            result.memory = memory_at(dt_s);
            result.rigid = admit_and_canonicalize_state(body, &result.rigid, dt_s)?;
            if !result.hinge_q.is_finite()
                || !result.hinge_v.is_finite()
                || result.memory.iter().any(|m| !m.value.is_finite())
            {
                return Err(refuse(
                    "coupled-state-invalid",
                    "nonfinite accepted state".into(),
                    "reduce the step",
                ));
            }
            return Ok(CoupledStep {
                state: result,
                iterations: iteration,
                normalized_residual: residual,
            });
        }
        delta = next;
    }
    Err(refuse(
        "coupled-step-not-converged",
        format!("midpoint velocity residual/budget {residual:e}"),
        "reduce the timestep or increase the bounded midpoint iteration budget",
    ))
}
