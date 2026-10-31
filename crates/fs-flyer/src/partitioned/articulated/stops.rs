//! Mechanical travel stops on the original coupled articulated midpoint lane.
//!
//! Instantaneous impacts share the original simultaneous normal-impact law
//! with point contacts, using the physical dry-plus-fluid mass. Resting stops
//! supply unilateral midpoint reactions through that same step's effective
//! mass, including its implicit spring/damper diagonal. Neither phase replaces
//! the mechanism by independent scalar inertias.
//!
//! The existing ballistic owner checks the midpoint step's quadratic joint
//! paths, including return after departure. Trials are reintegrated from the
//! unchanged state and shortened to an admissible-side contact. This is a
//! bounded numerical event lane, not certified CCD of the exact nonlinear
//! trajectory. The caller resolves smooth motion with `maximum_step_s`.

use super::impact::{GeneralizedNormalRow, resolve_normal_rows};
use super::{
    ArticulatedLoads, ArticulatedState, ArticulatedStepConfig, JointImpedance, admit_state,
    cancelled, native_error, step_with_velocity_response,
};
use crate::addedmass::joint::{admit_fluid, native_index};
use crate::addedmass::{MAX_CONTROL_COORDS, PreparedArticulatedMass};
use crate::{Refusal, partitioned::refuse};
use fs_contact::ballistic::{BallisticContactInput, first_ballistic_contact};
use fs_contact::multicontact::MultiImpactConfig;
use fs_ga::{Twist, Vec3};
use fs_mbd::articulated::ArticulatedModel;
use fs_mbd::articulated::generalized::{GeneralizedKinematics, free_floating_terms};
use fs_solver::orthant::active_set::{SpdOrthantError, solve_spd_orthant};

/// One joint's mechanical lower and upper stops, in native radians or metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointTravel {
    /// Compact scalar-joint index; each joint appears at most once.
    pub joint: usize,
    pub lower: f64,
    pub upper: f64,
    /// Positive contact resolution in this joint's coordinate units.
    /// No negative configuration gap is admitted or projected away.
    pub activation_tolerance: f64,
    /// Passive Newton restitution in [0,1], used only on closing impacts.
    pub restitution: f64,
}

/// Side of a mechanical travel interval; its normal points into the interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StopSide {
    Lower,
    Upper,
}

impl StopSide {
    fn sign(self) -> f64 {
        if self == Self::Lower { 1.0 } else { -1.0 }
    }
    fn index(self) -> usize {
        if self == Self::Lower { 0 } else { 1 }
    }
}

/// Bounds for one atomically accepted event window.
#[derive(Clone, Copy, Debug)]
pub struct StopStepConfig<'a> {
    pub smooth: ArticulatedStepConfig<'a>,
    /// Positive maximum smooth interval; resolve nonlinear motion and forcing.
    pub maximum_step_s: f64,
    /// Positive minimum attempted interval; unresolved smaller motion refuses.
    pub minimum_interval_s: f64,
    /// Counts every smooth attempt, including event/domain retries; at most4096.
    pub maximum_attempts: usize,
    /// At most64 simultaneous impact events; zero permits only smooth/support motion.
    pub maximum_events: usize,
    /// Original simultaneous impact and positive-definite orthant controls.
    pub contact: MultiImpactConfig,
}

/// A simultaneous stop impact at the reported numerical event time.
#[derive(Clone, Debug, PartialEq)]
pub struct StopEvent {
    pub time_s: f64,
    /// `(index in the supplied travel list, side)`, matching impulse order.
    pub boundaries: Vec<(usize, StopSide)>,
    /// Actual nonnegative activation gaps; positions are unchanged by impact.
    pub gaps: Vec<f64>,
    pub normal_impulses: Vec<f64>,
    pub kinetic_energy_before_j: f64,
    pub kinetic_energy_after_j: f64,
    pub dissipated_energy_j: f64,
    pub energy_balance_residual_j: f64,
}

/// Complete accepted window, including impact loss and resting-stop reactions.
#[derive(Clone, Debug, PartialEq)]
pub struct StopStep {
    pub state: ArticulatedState,
    pub events: Vec<StopEvent>,
    /// Integrated compressive reaction impulses `[lower,upper]` per travel row.
    /// Units are N m s for revolute joints, N s for prismatic joints.
    pub reaction_impulses: Vec<[f64; 2]>,
    /// Reactions paired with the actual accepted midpoint joint velocities.
    pub reaction_work_j: f64,
    pub impact_dissipation_j: f64,
    pub accepted_intervals: usize,
    pub attempted_intervals: usize,
}

fn invalid(message: impl ToString) -> Refusal {
    refuse(
        "joint-stop-invalid",
        message.to_string(),
        "inspect travel geometry and the explicit event, contact and smooth-step budgets",
    )
}
fn unresolved(message: impl ToString) -> Refusal {
    refuse(
        "joint-stop-unresolved",
        message.to_string(),
        "resolve the contact with smaller smooth steps or an appropriate activation/time budget",
    )
}
fn gap(stop: JointTravel, side: StopSide, q: f64) -> f64 {
    if side == StopSide::Lower {
        q - stop.lower
    } else {
        stop.upper - q
    }
}
fn velocity(state: &ArticulatedState) -> Vec<f64> {
    let native = state.motion.base.twist_body.to_array();
    let mut u = Vec::with_capacity(6 + state.motion.joint_velocity.len());
    for i in 0..6 {
        u.push(native[native_index(i)]);
    }
    u.extend_from_slice(&state.motion.joint_velocity);
    u
}
fn set_velocity(state: &mut ArticulatedState, u: &[f64]) {
    state.motion.base.twist_body =
        Twist::new(Vec3::new(u[3], u[4], u[5]), Vec3::new(u[0], u[1], u[2]));
    state.motion.joint_velocity.copy_from_slice(&u[6..]);
}
fn nearby(
    stops: &[JointTravel],
    state: &ArticulatedState,
) -> Result<Vec<(usize, StopSide)>, Refusal> {
    let mut rows = Vec::new();
    for (i, stop) in stops.iter().enumerate() {
        for side in [StopSide::Lower, StopSide::Upper] {
            let g = gap(*stop, side, state.motion.joint_position[stop.joint]);
            if !g.is_finite() || g < 0.0 {
                return Err(unresolved(
                    "a joint configuration is outside its mechanical travel interval",
                ));
            }
            if g <= stop.activation_tolerance {
                rows.push((i, side));
            }
        }
    }
    Ok(rows)
}

// Solve midpoint-velocity complementarity. The unknowns are physical reaction
// impulses, so H = J M_eff^-1 J^T / 2 and b = J(u0 + delta_free/2).
fn support(
    mass: &PreparedArticulatedMass,
    increment: &mut [f64],
    initial: &[f64],
    stops: &[JointTravel],
    active: &[(usize, StopSide)],
    config: MultiImpactConfig,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<Vec<f64>, Refusal> {
    if active.is_empty() {
        return Ok(Vec::new());
    }
    let n = increment.len();
    let count = active.len();
    let mut responses = Vec::with_capacity(count);
    let mut linear = Vec::with_capacity(count);
    let mut arithmetic_scale = Vec::with_capacity(count);
    for &(i, side) in active {
        let coordinate = 6 + stops[i].joint;
        let mut row = vec![0.0; n];
        row[coordinate] = side.sign();
        responses.push(mass.solve(&row)?);
        linear.push(side.sign() * (initial[coordinate] + 0.5 * increment[coordinate]));
        arithmetic_scale.push(1.0 + initial[coordinate].abs() + increment[coordinate].abs());
    }
    let mut hessian = vec![0.0; count * count];
    for (i, &(stop, side)) in active.iter().enumerate() {
        for j in 0..count {
            hessian[i * count + j] = 0.5 * side.sign() * responses[j][6 + stops[stop].joint];
        }
    }
    let solved = solve_spd_orthant(&hessian, &linear, config.solver, cancelled).map_err(|e| {
        if matches!(e, SpdOrthantError::Cancelled) {
            super::cancelled()
        } else {
            unresolved(e)
        }
    })?;
    for (response, impulse) in responses.iter().zip(&solved.point) {
        for (value, change) in increment.iter_mut().zip(response) {
            *value += impulse * change;
        }
        for (row, &(i, _)) in active.iter().enumerate() {
            arithmetic_scale[row] += (impulse * response[6 + stops[i].joint]).abs();
        }
    }
    for (row, &(i, side)) in active.iter().enumerate() {
        if solved.point[row] > 0.0 {
            let coordinate = 6 + stops[i].joint;
            let normal_midpoint = side.sign() * (initial[coordinate] + 0.5 * increment[coordinate]);
            if !normal_midpoint.is_finite()
                || !arithmetic_scale[row].is_finite()
                || normal_midpoint.abs() > 128.0 * f64::EPSILON * arithmetic_scale[row]
            {
                return Err(unresolved(
                    "resting reaction did not resolve its midpoint constraint",
                ));
            }
            // Enforce the solved zero-rate coordinate only within arithmetic
            // roundoff, independently of user tolerances. No position is projected and every other
            // coordinate retains the same full-mass reaction response.
            increment[coordinate] = -2.0 * initial[coordinate];
        } else if side.sign() * (initial[6 + stops[i].joint] + 0.5 * increment[6 + stops[i].joint])
            < 0.0
        {
            return Err(unresolved("the reaction solve left a closing resting stop"));
        }
    }
    Ok(solved.point)
}

/// Advance through simultaneous impacts, unilateral holding, release and return.
///
/// The load callback is the original pure articulated callback. Instantaneous
/// evaluations provide fluid inertia at the actual event; applied forces have
/// no impulsive contribution. Held aero-memory inputs advance only with accepted
/// smooth time. Mechanical stops do not replace or widen native joint limits.
/// Only a specifically typed native position-limit trial may be shortened;
/// other model, load, cancellation and solver refusals propagate unchanged.
///
/// Event geometry follows each smooth midpoint candidate's complete quadratic
/// joint path, reintegrated after shortening. Activations may occur early only
/// within the declared positive gap band. This is not an exact nonlinear path
/// certificate. Persistent positive reactions constrain their midpoint rate;
/// noncompressive reactions release automatically. Events at the final time
/// are included. Exhausted budgets, unresolved chatter or negative gaps refuse
/// the whole window, leaving every borrowed input untouched.
#[allow(clippy::too_many_arguments)]
pub fn step_with_stops<F>(
    model: &ArticulatedModel,
    impedance: &[JointImpedance],
    state: &ArticulatedState,
    memory_inputs: &[f64],
    t_s: f64,
    dt_s: f64,
    stops: &[JointTravel],
    config: StopStepConfig<'_>,
    evaluate: F,
    mut is_cancelled: impl FnMut() -> bool,
) -> Result<StopStep, Refusal>
where
    F: Fn(f64, &ArticulatedState, &GeneralizedKinematics) -> Result<ArticulatedLoads, Refusal>,
{
    if is_cancelled() {
        return Err(cancelled());
    }
    let nc = model.dof_count();
    if nc > MAX_CONTROL_COORDS
        || stops.len() > nc
        || state.motion.joint_position.len() != nc
        || state.motion.joint_velocity.len() != nc
        || !t_s.is_finite()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !(t_s + dt_s).is_finite()
        || t_s + dt_s <= t_s
        || !config.maximum_step_s.is_finite()
        || config.maximum_step_s <= 0.0
        || !config.minimum_interval_s.is_finite()
        || config.minimum_interval_s <= 0.0
        || config.minimum_interval_s > config.maximum_step_s
        || config.maximum_attempts == 0
        || config.maximum_attempts > 4096
        || config.maximum_events > 64
        || !stops.iter().enumerate().all(|(i, s)| {
            s.joint < nc
                && s.lower.is_finite()
                && s.upper.is_finite()
                && s.lower < s.upper
                && s.activation_tolerance.is_finite()
                && s.activation_tolerance > 0.0
                && s.activation_tolerance < 0.5 * (s.upper - s.lower)
                && s.restitution.is_finite()
                && (0.0..=1.0).contains(&s.restitution)
                && stops[..i].iter().all(|p| p.joint != s.joint)
        })
    {
        return Err(invalid(
            "inadmissible travel intervals, state dimensions or event budgets",
        ));
    }
    admit_state(model, &state.motion)?;
    nearby(stops, state)?;
    if !config.contact.energy_tolerance.is_finite() || config.contact.energy_tolerance < 0.0 {
        return Err(invalid("invalid contact energy budget"));
    }
    solve_spd_orthant(&[], &[], config.contact.solver, &mut is_cancelled).map_err(unresolved)?;
    let mut out = StopStep {
        state: state.clone(),
        events: Vec::new(),
        reaction_impulses: vec![[0.0; 2]; stops.len()],
        reaction_work_j: 0.0,
        impact_dissipation_j: 0.0,
        accepted_intervals: 0,
        attempted_intervals: 0,
    };
    let end = t_s + dt_s;
    let mut time = t_s;
    loop {
        if is_cancelled() {
            return Err(cancelled());
        }
        let touching = nearby(stops, &out.state)?;
        if touching
            .iter()
            .any(|&(i, side)| side.sign() * out.state.motion.joint_velocity[stops[i].joint] < 0.0)
        {
            if out.events.len() == config.maximum_events {
                return Err(unresolved("simultaneous stop-event budget exhausted"));
            }
            let dry = free_floating_terms(
                model,
                out.state.motion.base,
                &out.state.motion.joint_position,
                &out.state.motion.joint_velocity,
                &mut is_cancelled,
            )
            .map_err(native_error)?;
            let loads = evaluate(time, &out.state, &dry.kinematics)?;
            admit_fluid(&loads.fluid)?;
            if loads
                .fluid
                .q_rigid_nonaccel
                .iter()
                .chain(&loads.fluid.q_control_nonaccel)
                .any(|x| *x != 0.0)
            {
                return Err(refuse(
                    "coupled-load-ownership",
                    "fluid assembly contains applied loads at a stop".into(),
                    "supply noninertial loads through the original articulated load fields",
                ));
            }
            let mass = PreparedArticulatedMass::new(&dry, &loads.fluid)?;
            let before = velocity(&out.state);
            let rows: Vec<_> = touching
                .iter()
                .map(|&(i, side)| {
                    let stop = stops[i];
                    let mut row = vec![0.0; 6 + nc];
                    row[6 + stop.joint] = side.sign();
                    GeneralizedNormalRow {
                        jacobian: row,
                        gap: gap(stop, side, out.state.motion.joint_position[stop.joint]),
                        activation_tolerance: stop.activation_tolerance,
                        support_speed: 0.0,
                        restitution: stop.restitution,
                    }
                })
                .collect();
            let reset =
                resolve_normal_rows(&mass, &before, &rows, config.contact, &mut is_cancelled)?;
            let mut after = reset.velocity_after;
            let responses: Vec<_> = rows
                .iter()
                .map(|r| mass.solve(&r.jacobian))
                .collect::<Result<_, _>>()?;
            for (r, &(i, _)) in touching.iter().enumerate() {
                if reset.law.restitution_targets[r] == 0.0 {
                    let coordinate = 6 + stops[i].joint;
                    let scale = 1.0
                        + before[coordinate].abs()
                        + responses
                            .iter()
                            .zip(&reset.law.normal_impulses)
                            .map(|(v, p)| (v[coordinate] * p).abs())
                            .sum::<f64>();
                    if !scale.is_finite() {
                        return Err(unresolved("unrepresentable stop-reset arithmetic scale"));
                    }
                    if after[coordinate].abs() <= 128.0 * f64::EPSILON * scale {
                        after[coordinate] = 0.0;
                    } else if after[coordinate] * touching[r].1.sign() < 0.0 {
                        return Err(unresolved("impact left an unresolved closing stop"));
                    }
                }
            }
            // Re-audit the actual roundoff-canonicalized reset, including all
            // six base momenta. A loose contact budget cannot permit cleanup.
            let change: Vec<_> = after.iter().zip(&before).map(|(a, b)| a - b).collect();
            let impulse = mass.momentum(&change)?;
            let mut expected = vec![0.0; 6 + nc];
            for (row, p) in rows.iter().zip(&reset.law.normal_impulses) {
                for (x, j) in expected.iter_mut().zip(&row.jacobian) {
                    *x += p * j;
                }
            }
            let p0 = mass.momentum(&before)?;
            let scale = 1.0
                + p0.iter().map(|x| x.abs()).sum::<f64>()
                + impulse.iter().map(|x| x.abs()).sum::<f64>()
                + expected.iter().map(|x| x.abs()).sum::<f64>();
            if !scale.is_finite()
                || impulse
                    .iter()
                    .zip(&expected)
                    .any(|(a, b)| (a - b).abs() > 512.0 * f64::EPSILON * scale)
            {
                return Err(unresolved(
                    "stop reset did not retain full coupled momentum",
                ));
            }
            let kinetic_after = mass.kinetic_energy(&after)?;
            let kinetic_change: f64 = after
                .iter()
                .zip(&before)
                .zip(&impulse)
                .map(|((a, b), p)| a.midpoint(*b) * p)
                .sum();
            let balance = kinetic_change + reset.law.dissipated_energy;
            let energy_scale = 1.0
                + reset.kinetic_before.abs()
                + kinetic_after.abs()
                + reset.law.dissipated_energy;
            if !kinetic_change.is_finite()
                || !balance.is_finite()
                || !energy_scale.is_finite()
                || balance.abs()
                    > config.contact.energy_tolerance + 512.0 * f64::EPSILON * energy_scale
            {
                return Err(unresolved("stop reset energy balance did not resolve"));
            }
            set_velocity(&mut out.state, &after);
            admit_state(model, &out.state.motion)?;
            out.impact_dissipation_j += reset.law.dissipated_energy;
            if !out.impact_dissipation_j.is_finite() {
                return Err(unresolved("unrepresentable accumulated impact loss"));
            }
            out.events.push(StopEvent {
                time_s: time,
                boundaries: touching,
                gaps: rows.iter().map(|r| r.gap).collect(),
                normal_impulses: reset.law.normal_impulses,
                kinetic_energy_before_j: reset.kinetic_before,
                kinetic_energy_after_j: kinetic_after,
                dissipated_energy_j: reset.law.dissipated_energy,
                energy_balance_residual_j: balance,
            });
        }
        if time == end {
            return Ok(out);
        }
        let initial = velocity(&out.state);
        let resting: Vec<_> = nearby(stops, &out.state)?
            .into_iter()
            .filter(|&(i, _)| out.state.motion.joint_velocity[stops[i].joint] == 0.0)
            .collect();
        let mut duration = (end - time).min(config.maximum_step_s);
        // An admissible-side linear predictor avoids probing far beyond a
        // native position boundary just to discover an already approaching stop.
        for (i, stop) in stops.iter().enumerate() {
            if resting.iter().any(|&(j, _)| i == j) {
                continue;
            }
            for side in [StopSide::Lower, StopSide::Upper] {
                let query = BallisticContactInput {
                    gap: gap(*stop, side, out.state.motion.joint_position[stop.joint]),
                    relative_normal_speed: side.sign()
                        * out.state.motion.joint_velocity[stop.joint],
                    relative_normal_acceleration: 0.0,
                    duration_s: duration,
                };
                if let Some(hit) = first_ballistic_contact(query).map_err(unresolved)? {
                    if hit.time_s > 0.0 {
                        duration = duration.min(hit.time_s);
                    }
                }
            }
        }
        loop {
            if is_cancelled() {
                return Err(cancelled());
            }
            if out.attempted_intervals == config.maximum_attempts {
                return Err(unresolved("smooth/event attempt budget exhausted"));
            }
            if (duration < config.minimum_interval_s && duration < end - time)
                || time + duration <= time
            {
                return Err(unresolved(
                    "stop interval cannot make admitted time progress",
                ));
            }
            out.attempted_intervals += 1;
            let mut reaction = Vec::new();
            let mut previous_response = None;
            let trial = step_with_velocity_response(
                model,
                impedance,
                &out.state,
                memory_inputs,
                time,
                duration,
                config.smooth,
                &evaluate,
                |mass, increment, cancel| {
                    let next = support(
                        mass,
                        increment,
                        &initial,
                        stops,
                        &resting,
                        config.contact,
                        cancel,
                    )?;
                    if resting.is_empty() {
                        return Ok(true);
                    }
                    let ready = previous_response.is_some();
                    reaction = previous_response.replace(next).unwrap_or_default();
                    Ok(ready)
                },
                &mut is_cancelled,
            );
            let trial = match trial {
                Ok(value) => value,
                Err(error) if error.code == "articulated-joint-position-limit" => {
                    duration *= 0.5;
                    continue;
                }
                Err(error) => return Err(error),
            };
            let mut first_hit = duration;
            let mut geometry_valid = true;
            for (i, stop) in stops.iter().enumerate() {
                let q0 = out.state.motion.joint_position[stop.joint];
                let q1 = trial.state.motion.joint_position[stop.joint];
                let held = resting.iter().enumerate().any(|(r, &(j, _))| {
                    j == i
                        && reaction[r] > 0.0
                        && q1 == q0
                        && trial.state.motion.joint_velocity[stop.joint] == 0.0
                });
                for side in [StopSide::Lower, StopSide::Upper] {
                    let endpoint = gap(*stop, side, q1);
                    geometry_valid &= endpoint.is_finite() && endpoint >= 0.0;
                    if held {
                        continue;
                    }
                    let speed = side.sign() * out.state.motion.joint_velocity[stop.joint];
                    let acceleration = side.sign()
                        * (trial.state.motion.joint_velocity[stop.joint]
                            - out.state.motion.joint_velocity[stop.joint])
                        / duration;
                    let query = BallisticContactInput {
                        gap: gap(*stop, side, q0),
                        relative_normal_speed: speed,
                        relative_normal_acceleration: acceleration,
                        duration_s: duration,
                    };
                    if let Some(hit) = first_ballistic_contact(query).map_err(unresolved)? {
                        if hit.relative_normal_speed < 0.0 {
                            first_hit = first_hit.min(hit.time_s);
                        }
                    }
                }
            }
            if !geometry_valid || first_hit < duration {
                duration = if first_hit > 0.0 && first_hit < duration {
                    first_hit.next_down()
                } else {
                    0.5 * duration
                };
                continue;
            }
            for (r, &(i, side)) in resting.iter().enumerate() {
                let midpoint = out.state.motion.joint_velocity[stops[i].joint]
                    .midpoint(trial.state.motion.joint_velocity[stops[i].joint]);
                out.reaction_impulses[i][side.index()] += reaction[r];
                out.reaction_work_j += reaction[r] * side.sign() * midpoint;
            }
            if !out.reaction_work_j.is_finite()
                || !out
                    .reaction_impulses
                    .iter()
                    .flatten()
                    .all(|x| x.is_finite())
            {
                return Err(unresolved(
                    "unrepresentable accumulated stop-reaction receipt",
                ));
            }
            time = (time + duration).min(end);
            out.state = trial.state;
            out.accepted_intervals += 1;
            break;
        }
    }
}
