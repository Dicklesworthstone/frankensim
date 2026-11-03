//! Smooth mixed-boundary motion joined by caller-located Coulomb impacts.
//!
//! Original free and prescribed ABA motion owners advance each retained tree.
//! One native simultaneous mixed-boundary event changes velocities at fixed
//! configuration. Support paths, collision times and contact geometry belong
//! to the caller; no artificial support inertia or path interpolation is used.

/// Automatic encounters of spherical free and prescribed articulated carriers.
pub mod spheres;

use super::{
    ArticulatedTrajectoryConfig, ArticulatedTrajectoryError, ArticulatedTrajectoryEvent,
    ArticulatedTrajectoryMomentum, Impulses, event_schedule, finite_vector, poll, smooth_error,
};
use crate::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath, step_prescribed_base,
};
use crate::articulated::motion::{
    FreeFloatingMotionLoads, FreeFloatingMotionState, MAX_FREE_FLOATING_STEP_LINKS,
    step_free_floating,
};
use crate::articulated::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, forward_kinematics,
    free_floating_momentum, momentum_from_kinematics,
};
use crate::articulated_contacts::ArticulatedContactError;
use crate::articulated_contacts::friction::ArticulatedFrictionalContactResult;
use crate::articulated_contacts::mixed::{
    ArticulatedContactBoundary, MixedArticulatedContactSystem, MixedArticulatedImpulseResponse,
    resolve_mixed_articulated_frictional_contacts,
};
use fs_contact::frictional_multicontact::FrictionalContactResult;
use fs_ga::{Se3, Vec3};
use fs_math::dd::Dd;

/// Actual native state of one free or supported original articulated tree.
#[derive(Debug, Clone, PartialEq)]
pub enum MixedArticulatedTrajectoryState {
    /// Free base and joints, with physical base recoil.
    Free(FreeFloatingMotionState),
    /// Prescribed support pose/twist/acceleration, clock and movable joints.
    Prescribed(PrescribedBaseMotionState),
}

impl MixedArticulatedTrajectoryState {
    fn contact_system<'a>(
        &'a self,
        model: &'a ArticulatedModel,
    ) -> MixedArticulatedContactSystem<'a> {
        let (boundary, joint_position, joint_velocity) = match self {
            Self::Free(state) => (
                ArticulatedContactBoundary::Free(state.base),
                &state.joint_position,
                &state.joint_velocity,
            ),
            Self::Prescribed(state) => (
                ArticulatedContactBoundary::Prescribed(state.base),
                &state.joint_position,
                &state.joint_velocity,
            ),
        };
        MixedArticulatedContactSystem {
            model,
            boundary,
            joint_position,
            joint_velocity,
        }
    }
}

/// Immutable initial native state and held physical loads for one tree.
#[derive(Debug, Clone, Copy)]
pub struct MixedArticulatedTrajectorySystem<'a> {
    /// Original model; a prescribed tree may retain a movable scalar root.
    pub model: &'a ArticulatedModel,
    /// Complete state at the common window clock.
    pub state: &'a MixedArticulatedTrajectoryState,
    /// Held relative actuator efforts, link-frame wrenches and world gravity.
    pub loads: FreeFloatingMotionLoads<'a>,
}

/// One complete accepted simultaneous native event.
#[derive(Debug, Clone)]
pub struct MixedArticulatedTrajectoryEventResult {
    /// Caller-located elapsed time from the window start, in seconds.
    pub time: f64,
    /// Common native clock at the event, in seconds.
    pub time_s: f64,
    /// Original impulses, support reactions, determinacy and endpoint audits.
    pub impact: ArticulatedFrictionalContactResult<
        FrictionalContactResult,
        MixedArticulatedImpulseResponse,
    >,
}

/// Complete window, with reusable states and independently evaluated balances.
#[derive(Debug, Clone)]
pub struct MixedArticulatedTrajectoryResult {
    /// Common initial clock, in seconds.
    pub time_before_s: f64,
    /// Common accepted clock, accumulated by the original native steps.
    pub time_after_s: f64,
    /// Complete endpoints in original system order, including support clocks.
    pub states_after: Vec<MixedArticulatedTrajectoryState>,
    /// Accepted simultaneous events in chronological order.
    pub events: Vec<MixedArticulatedTrajectoryEventResult>,
    /// Number of common smooth intervals actually advanced.
    pub smooth_intervals: usize,
    /// Fixed world origin for all aggregate angular quantities.
    pub angular_origin_world: Vec3,
    /// Physical totals from original FK and link inertias at the initial states.
    pub momentum_before: ArticulatedTrajectoryMomentum,
    /// Independently recomputed totals at the actual accepted endpoints.
    pub momentum_after: ArticulatedTrajectoryMomentum,
    /// Actual native stage/event linear increments, reduced before rounding.
    /// This retains small transfers beside large counter-moving systems.
    pub linear_momentum_change_world: Vec3,
    /// Corresponding angular increments about the fixed origin, in N m s.
    pub angular_momentum_change_about_origin_world: Vec3,
    /// Sum of actual native kinetic increments, in joules.
    pub kinetic_energy_change: f64,
    /// Applied body/gravity, plane and prescribed-base impulses, in N s.
    /// Dynamic contact pairs cancel; support reactions act ON retained trees.
    pub external_linear_impulse_world: Vec3,
    /// Corresponding angular impulse about the declared origin, in N m s.
    pub external_angular_impulse_about_origin_world: Vec3,
    /// Smooth relative actuator effort/rate work, in joules.
    pub actuator_work: f64,
    /// Smooth held body-wrench and gravity work, excluding base support, J.
    pub smooth_external_work: f64,
    /// Signed native smooth base-support work, in joules.
    pub smooth_support_work: f64,
    /// Sum of determined total plane/base work at instantaneous events, J.
    /// Individual sources can remain representative under redundant contact.
    pub contact_boundary_work: f64,
    /// Total smooth external/support and instantaneous boundary work, J.
    pub external_work: f64,
    /// Admitted mechanical impact loss, without an inferred heat partition, J.
    pub contact_dissipated_energy: f64,
    /// Sum of native smooth numerical kinetic-work defects, never heat, J.
    pub smooth_kinetic_work_residual: f64,
    /// Sum of smooth actual-minus-integrated momentum increments, in N s.
    /// Prescribed-path midpoint integration need not conserve this exactly.
    pub smooth_linear_momentum_residual: Vec3,
    /// Corresponding smooth angular numerical defect about the origin, N m s.
    pub smooth_angular_momentum_residual: Vec3,
    /// Actual complete-window kinetic change minus work plus impact loss, J.
    pub kinetic_work_residual: f64,
    /// Actual complete-window linear change minus external impulse, N s.
    pub linear_momentum_residual: Vec3,
    /// Actual angular change minus external angular impulse, N m s.
    pub angular_momentum_residual: Vec3,
}

/// Advance free and prescribed-base trees through explicit simultaneous impacts.
///
/// `start_time_s` is the common physical clock. Every prescribed input must
/// already carry that exact clock. Events retain elapsed times from the existing
/// trajectory API. Each smooth interval advances the common clock once by its
/// actual duration, without resetting native time after an impact. Ordinary
/// floating-point accumulation may differ from one addition of the full window.
///
/// `support_path(system_index, accepted_state, duration)` is called once per
/// prescribed tree per smooth interval, in original system order. It supplies
/// native midpoint/endpoint samples for that interval; the initial support is
/// the accepted state and cannot be replaced. After an impact the callback sees
/// the actual updated joint rates and unchanged support. Path samples must be
/// smooth and derivative-consistent under the native spatial-acceleration
/// convention. This composition does not certify that input obligation.
///
/// The original bounded event schedule, O(links) smooth stages and small native
/// Coulomb system are reused. No collision detection, sustained contact, pose
/// projection or fictitious support mass is introduced. All input states are
/// immutable, and a late path/event/limit refusal or cancellation publishes no
/// partial window. Cancellation is polled around callbacks and native stages;
/// callback-owned side effects cannot be rolled back. Use pure path callbacks
/// when deterministic retry is required.
///
/// Smooth force quadrature defects and actual global balances remain explicit;
/// they are neither physical dissipation nor corrected momenta. Angular impulse
/// translation retains separate source/origin terms. Signed work and native
/// actual increments are reduced before rounding; global residuals do not
/// subtract large rounded whole-window totals. Both smooth owners and native
/// impulses supply compensated actual increments. Extreme native source
/// magnitudes still have finite precision.
pub fn advance_mixed_articulated_events(
    systems: &[MixedArticulatedTrajectorySystem<'_>],
    events: &[ArticulatedTrajectoryEvent<'_>],
    start_time_s: f64,
    config: ArticulatedTrajectoryConfig,
    mut support_path: impl FnMut(
        usize,
        &PrescribedBaseMotionState,
        f64,
    ) -> Result<PrescribedBasePath, PrescribedBaseMotionError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixedArticulatedTrajectoryResult, ArticulatedTrajectoryError> {
    poll(&mut cancelled)?;
    let counts = event_schedule(systems.len(), events, config)?;
    admit_states(systems, start_time_s)?;
    let mut states: Vec<_> = systems.iter().map(|system| system.state.clone()).collect();
    let momentum_before = momentum(
        systems,
        &states,
        config.angular_origin_world,
        0.0,
        &mut cancelled,
    )?;
    let mut impulses = Impulses::default();
    let mut momentum_change = Impulses::default();
    let mut smooth_momentum_defect = Impulses::default();
    let mut kinetic_change = Dd::ZERO;
    let mut actuator_work = Dd::ZERO;
    let mut smooth_external_work = Dd::ZERO;
    let mut smooth_support_work = Dd::ZERO;
    let mut contact_boundary_work = Dd::ZERO;
    let mut loss = Dd::ZERO;
    let mut smooth_defect = Dd::ZERO;
    let mut accepted_events = Vec::with_capacity(events.len());
    let mut clock = start_time_s;
    let mut start = 0.0;
    for (interval_index, &count) in counts.iter().enumerate() {
        let end = events
            .get(interval_index)
            .map_or(config.duration, |event| event.time);
        let mut previous_time = start;
        for index in 0..count {
            poll(&mut cancelled)?;
            let next_time = if index + 1 == count {
                end
            } else {
                start + (end - start) * ((index + 1) as f64 / count as f64)
            };
            let duration = next_time - previous_time;
            let midpoint_clock = clock + 0.5 * duration;
            let next_clock = clock + duration;
            if duration <= 0.0
                || duration > config.maximum_smooth_step * (1.0 + 8.0 * f64::EPSILON)
                || !next_clock.is_finite()
                || !(midpoint_clock > clock && next_clock > midpoint_clock)
            {
                return Err(ArticulatedTrajectoryError::Invalid(
                    "smooth interval or common clock is not representable",
                ));
            }
            for (system_index, (system, state)) in systems.iter().zip(&mut states).enumerate() {
                poll(&mut cancelled)?;
                match state {
                    MixedArticulatedTrajectoryState::Free(current) => {
                        let step = step_free_floating(
                            system.model,
                            current,
                            system.loads,
                            duration,
                            &mut cancelled,
                        )
                        .map_err(|error| smooth_error(system_index, next_time, error))?;
                        impulses.add(
                            step.external_linear_impulse_world,
                            step.external_angular_impulse_world,
                            config.angular_origin_world,
                        );
                        momentum_change.add(
                            step.linear_momentum_change_world,
                            step.angular_momentum_change_world,
                            config.angular_origin_world,
                        );
                        smooth_momentum_defect.add(
                            step.linear_momentum_residual_world,
                            step.angular_momentum_residual_world,
                            config.angular_origin_world,
                        );
                        kinetic_change = kinetic_change + Dd::from_f64(step.kinetic_energy_change);
                        actuator_work = actuator_work + Dd::from_f64(step.actuator_work);
                        smooth_external_work =
                            smooth_external_work + Dd::from_f64(step.external_work);
                        smooth_defect = smooth_defect + Dd::from_f64(step.kinetic_work_residual);
                        *current = step.state_after;
                    }
                    MixedArticulatedTrajectoryState::Prescribed(current) => {
                        let path = support_path(system_index, current, duration)
                            .map_err(|error| prescribed_error(system_index, next_time, error))?;
                        poll(&mut cancelled)?;
                        let step = step_prescribed_base(
                            system.model,
                            current,
                            path,
                            system.loads,
                            duration,
                            &mut cancelled,
                        )
                        .map_err(|error| prescribed_error(system_index, next_time, error))?;
                        impulses.add(
                            step.external_linear_impulse_world,
                            step.external_angular_impulse_world,
                            config.angular_origin_world,
                        );
                        impulses.add(
                            step.support_linear_impulse_world,
                            step.support_angular_impulse_world,
                            config.angular_origin_world,
                        );
                        smooth_momentum_defect.add(
                            step.linear_momentum_residual_world,
                            step.angular_momentum_residual_world,
                            config.angular_origin_world,
                        );
                        momentum_change.add(
                            step.linear_momentum_change_world,
                            step.angular_momentum_change_world,
                            config.angular_origin_world,
                        );
                        kinetic_change = kinetic_change + Dd::from_f64(step.kinetic_energy_change);
                        actuator_work = actuator_work + Dd::from_f64(step.actuator_work);
                        smooth_external_work =
                            smooth_external_work + Dd::from_f64(step.external_work);
                        smooth_support_work = smooth_support_work + Dd::from_f64(step.support_work);
                        smooth_defect = smooth_defect + Dd::from_f64(step.kinetic_work_residual);
                        *current = step.state_after;
                    }
                }
            }
            clock = next_clock;
            previous_time = next_time;
        }
        start = end;
        if let Some(event) = events.get(interval_index) {
            poll(&mut cancelled)?;
            let inputs: Vec<_> = systems
                .iter()
                .zip(&states)
                .map(|(system, state)| state.contact_system(system.model))
                .collect();
            let impact = resolve_mixed_articulated_frictional_contacts(
                &inputs,
                event.contacts,
                config.contact,
                &mut cancelled,
            )
            .map_err(|error| {
                if matches!(error, ArticulatedContactError::Cancelled) {
                    ArticulatedTrajectoryError::Cancelled
                } else {
                    ArticulatedTrajectoryError::Contact {
                        event_index: interval_index,
                        time: event.time,
                        error,
                    }
                }
            })?;
            for (state, response) in states.iter_mut().zip(&impact.responses) {
                match (state, response) {
                    (
                        MixedArticulatedTrajectoryState::Free(current),
                        MixedArticulatedImpulseResponse::Free(response),
                    ) => {
                        impulses.add(
                            response.external_linear_impulse_world,
                            response.external_angular_impulse_world,
                            config.angular_origin_world,
                        );
                        momentum_change.add(
                            response.linear_momentum_change_world,
                            response.angular_momentum_change_world,
                            config.angular_origin_world,
                        );
                        kinetic_change =
                            kinetic_change + Dd::from_f64(response.kinetic_energy_change);
                        current.base = response.base_after;
                        current
                            .joint_velocity
                            .clone_from(&response.joint_velocity_after);
                    }
                    (
                        MixedArticulatedTrajectoryState::Prescribed(current),
                        MixedArticulatedImpulseResponse::Prescribed(response),
                    ) => {
                        impulses.add(
                            response.external_impulse_world.force,
                            response.external_impulse_world.torque,
                            config.angular_origin_world,
                        );
                        impulses.add(
                            response.support_impulse_world.force,
                            response.support_impulse_world.torque,
                            config.angular_origin_world,
                        );
                        momentum_change.add(
                            response.linear_momentum_change_world,
                            response.angular_momentum_change_world,
                            config.angular_origin_world,
                        );
                        kinetic_change =
                            kinetic_change + Dd::from_f64(response.kinetic_energy_change);
                        current.base = response.base_after;
                        current
                            .joint_velocity
                            .clone_from(&response.joint_velocity_after);
                    }
                    _ => {
                        return Err(ArticulatedTrajectoryError::Invalid(
                            "native event changed boundary kind",
                        ));
                    }
                }
            }
            contact_boundary_work = contact_boundary_work + Dd::from_f64(impact.boundary_work);
            loss = loss + Dd::from_f64(impact.contact.dissipated_energy);
            accepted_events.push(MixedArticulatedTrajectoryEventResult {
                time: event.time,
                time_s: clock,
                impact,
            });
        }
    }
    let momentum_after = momentum(
        systems,
        &states,
        config.angular_origin_world,
        config.duration,
        &mut cancelled,
    )?;
    let external_work = smooth_external_work + smooth_support_work + contact_boundary_work;
    let kinetic_work_residual = (kinetic_change - actuator_work - external_work + loss).to_f64();
    let (linear_momentum_residual, angular_momentum_residual) =
        momentum_change.difference(&impulses);
    let [
        kinetic_energy_change,
        actuator_work,
        smooth_external_work,
        smooth_support_work,
        contact_boundary_work,
        external_work,
        loss,
        smooth_defect,
    ] = [
        kinetic_change,
        actuator_work,
        smooth_external_work,
        smooth_support_work,
        contact_boundary_work,
        external_work,
        loss,
        smooth_defect,
    ]
    .map(Dd::to_f64);
    if ![
        kinetic_energy_change,
        actuator_work,
        smooth_external_work,
        smooth_support_work,
        contact_boundary_work,
        external_work,
        loss,
        smooth_defect,
        kinetic_work_residual,
    ]
    .iter()
    .all(|value| value.is_finite())
        || ![
            impulses.linear,
            impulses.angular,
            smooth_momentum_defect.linear,
            smooth_momentum_defect.angular,
            momentum_change.linear,
            momentum_change.angular,
            linear_momentum_residual,
            angular_momentum_residual,
        ]
        .into_iter()
        .all(finite_vector)
    {
        return Err(ArticulatedTrajectoryError::Invalid(
            "nonfinite mixed trajectory balance",
        ));
    }
    poll(&mut cancelled)?;
    Ok(MixedArticulatedTrajectoryResult {
        time_before_s: start_time_s,
        time_after_s: clock,
        states_after: states,
        events: accepted_events,
        smooth_intervals: counts.iter().sum(),
        angular_origin_world: config.angular_origin_world,
        momentum_before,
        momentum_after,
        linear_momentum_change_world: momentum_change.linear,
        angular_momentum_change_about_origin_world: momentum_change.angular,
        kinetic_energy_change,
        external_linear_impulse_world: impulses.linear,
        external_angular_impulse_about_origin_world: impulses.angular,
        actuator_work,
        smooth_external_work,
        smooth_support_work,
        contact_boundary_work,
        external_work,
        contact_dissipated_energy: loss,
        smooth_kinetic_work_residual: smooth_defect,
        smooth_linear_momentum_residual: smooth_momentum_defect.linear,
        smooth_angular_momentum_residual: smooth_momentum_defect.angular,
        kinetic_work_residual,
        linear_momentum_residual,
        angular_momentum_residual,
    })
}

fn admit_states(
    systems: &[MixedArticulatedTrajectorySystem<'_>],
    clock: f64,
) -> Result<(), ArticulatedTrajectoryError> {
    if !clock.is_finite() {
        return Err(ArticulatedTrajectoryError::Invalid(
            "nonfinite common clock",
        ));
    }
    let mut links = 0usize;
    for system in systems {
        let input = system.state.contact_system(system.model);
        if input.joint_position.len() != system.model.dof_count()
            || input.joint_velocity.len() != system.model.dof_count()
        {
            return Err(ArticulatedTrajectoryError::Invalid(
                "joint state lengths must match the original model DOF count",
            ));
        }
        if let MixedArticulatedTrajectoryState::Prescribed(state) = system.state {
            if state.time_s != clock {
                return Err(ArticulatedTrajectoryError::Invalid(
                    "prescribed state clock must equal the common window clock",
                ));
            }
        }
        links = links.checked_add(system.model.link_count()).ok_or(
            ArticulatedTrajectoryError::Invalid("trajectory link count overflow"),
        )?;
        if links > MAX_FREE_FLOATING_STEP_LINKS {
            return Err(ArticulatedTrajectoryError::Invalid(
                "trajectory exceeds total link cap",
            ));
        }
    }
    Ok(())
}

fn momentum(
    systems: &[MixedArticulatedTrajectorySystem<'_>],
    states: &[MixedArticulatedTrajectoryState],
    origin: Vec3,
    time: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ArticulatedTrajectoryMomentum, ArticulatedTrajectoryError> {
    let mut spatial = Impulses::default();
    let mut kinetic = Dd::ZERO;
    for (index, (system, state)) in systems.iter().zip(states).enumerate() {
        poll(cancelled)?;
        let value = native_momentum(system.model, state, origin).map_err(|error| match state {
            MixedArticulatedTrajectoryState::Free(_) => smooth_error(index, time, error.into()),
            MixedArticulatedTrajectoryState::Prescribed(_) => {
                prescribed_error(index, time, error.into())
            }
        })?;
        spatial.add(
            value.linear_momentum_world,
            value.angular_momentum_world,
            Vec3::default(),
        );
        kinetic = kinetic + Dd::from_f64(value.kinetic_energy);
    }
    let total = ArticulatedTrajectoryMomentum {
        linear_world: spatial.linear,
        angular_about_origin_world: spatial.angular,
        kinetic_energy: kinetic.to_f64(),
    };
    if !finite_vector(total.linear_world)
        || !finite_vector(total.angular_about_origin_world)
        || !total.kinetic_energy.is_finite()
    {
        return Err(ArticulatedTrajectoryError::Invalid(
            "nonfinite mixed total momentum",
        ));
    }
    Ok(total)
}

fn native_momentum(
    model: &ArticulatedModel,
    state: &MixedArticulatedTrajectoryState,
    origin: Vec3,
) -> Result<ArticulatedMomentum, ArticulatedError> {
    match state {
        MixedArticulatedTrajectoryState::Free(state) => {
            let mut base = state.base;
            base.world_from_base = Se3::from_parts(
                base.world_from_base.rotation(),
                base.world_from_base.translation() - origin,
            )?;
            free_floating_momentum(model, base, &state.joint_position, &state.joint_velocity)
        }
        MixedArticulatedTrajectoryState::Prescribed(state) => {
            let mut base = state.base;
            base.world_from_base = Se3::from_parts(
                base.world_from_base.rotation(),
                base.world_from_base.translation() - origin,
            )?;
            let kinematics =
                forward_kinematics(model, base, &state.joint_position, &state.joint_velocity)?;
            momentum_from_kinematics(model, &kinematics)
        }
    }
}

fn prescribed_error(
    system_index: usize,
    time: f64,
    error: PrescribedBaseMotionError,
) -> ArticulatedTrajectoryError {
    if matches!(error, PrescribedBaseMotionError::Cancelled) {
        ArticulatedTrajectoryError::Cancelled
    } else {
        ArticulatedTrajectoryError::PrescribedSmooth {
            system_index,
            time,
            error,
        }
    }
}
