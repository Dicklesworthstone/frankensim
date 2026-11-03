//! Smooth free-flight segments joined by caller-located simultaneous impacts.
//!
//! This owner advances the original free bases and joints, applies native
//! Coulomb impulses at the supplied times, and continues with the accepted
//! velocities. Event localization, feature coverage, and persistent contact
//! remain the geometry/contact owner's responsibility; a supplied event list
//! is not a no-missed-collision certificate.

/// Automatic ballistic plane encounters for spherical carriers with balanced rotors.
pub mod sphere_planes;

/// Caller-located impacts between free trees and smoothly prescribed mechanisms.
pub mod mixed;

use super::friction::{
    ArticulatedFrictionalContact, ArticulatedFrictionalContactResult,
    resolve_articulated_frictional_contacts,
};
use super::{ArticulatedContactError, ArticulatedContactSystem};
use crate::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    MAX_FREE_FLOATING_STEP_LINKS, step_free_floating,
};
use crate::articulated::{ArticulatedModel, ArticulatedMomentum, free_floating_momentum};
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_ga::{Se3, Vec3};
use fs_math::dd::Dd;

/// Immutable initial state and loads of one independently free articulated tree.
#[derive(Debug, Clone, Copy)]
pub struct ArticulatedTrajectorySystem<'a> {
    /// Validated original tree model.
    pub model: &'a ArticulatedModel,
    /// Base pose/twist and joint coordinates/rates at the window start.
    pub state: &'a FreeFloatingMotionState,
    /// Held internal efforts, link-frame wrenches, and uniform world gravity.
    pub loads: FreeFloatingMotionLoads<'a>,
}

/// One simultaneous contact set at a caller-located time from the step start.
#[derive(Debug, Clone, Copy)]
pub struct ArticulatedTrajectoryEvent<'a> {
    /// Seconds in [0, duration]. Times must be strictly increasing; contacts at
    /// the same time belong to one simultaneous solve, never sequential events.
    pub time: f64,
    /// Material points, normals, and prescribed-plane positions at THIS event.
    /// A moving support's event position is supplied explicitly, not advanced
    /// again by this owner.
    pub contacts: &'a [ArticulatedFrictionalContact],
}

/// Fixed bounds for one atomic trajectory window.
#[derive(Debug, Clone, Copy)]
pub struct ArticulatedTrajectoryConfig {
    /// Finite positive window duration, in seconds.
    pub duration: f64,
    /// Each event-free interval is split uniformly into steps no larger than
    /// this value. It is an accuracy control, not an error estimate.
    pub maximum_smooth_step: f64,
    /// Total common smooth intervals across the window, in 1..=4096.
    /// Each interval advances every supplied system once.
    pub maximum_smooth_intervals: usize,
    /// Fixed world-coordinate origin for all aggregate angular quantities.
    pub angular_origin_world: Vec3,
    /// Original simultaneous Coulomb solver controls, unchanged at every event.
    pub contact: FrictionalContactConfig,
}

/// Physical aggregate momentum evaluated from the actual original-tree states.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArticulatedTrajectoryMomentum {
    /// Total world linear momentum, in kg m/s.
    pub linear_world: Vec3,
    /// Spin plus orbital momentum about the fixed declared origin, in kg m²/s.
    pub angular_about_origin_world: Vec3,
    /// Total link kinetic energy, in joules.
    pub kinetic_energy: f64,
}

/// Full native evidence for one accepted event, in time order.
#[derive(Debug, Clone)]
pub struct ArticulatedTrajectoryEventResult {
    /// Accepted event time from the window start, in seconds.
    pub time: f64,
    /// Original simultaneous event result, including native post-impact states.
    pub impact: ArticulatedFrictionalContactResult,
}

/// Only a complete, admitted window publishes updated physical states.
#[derive(Debug, Clone)]
pub struct ArticulatedTrajectoryResult {
    /// Reusable endpoint states in input-system order.
    pub states_after: Vec<FreeFloatingMotionState>,
    /// Accepted native contact results in time order.
    pub events: Vec<ArticulatedTrajectoryEventResult>,
    /// Number of common smooth intervals actually advanced.
    pub smooth_intervals: usize,
    /// Fixed world-coordinate origin used by aggregate angular quantities.
    pub angular_origin_world: Vec3,
    /// Independently evaluated initial physical totals.
    pub momentum_before: ArticulatedTrajectoryMomentum,
    /// Independently evaluated accepted endpoint physical totals.
    pub momentum_after: ArticulatedTrajectoryMomentum,
    /// Held-load/gravity impulse plus the net prescribed-support impulse.
    /// Opposite impulses on dynamic counterparts cancel in the global sum.
    pub external_linear_impulse_world: Vec3,
    /// Total external angular impulse about the declared origin, in N m s.
    pub external_angular_impulse_about_origin_world: Vec3,
    /// Smooth internal actuator effort/rate work, in joules.
    pub actuator_work: f64,
    /// Smooth external work plus the native signed prescribed-support work.
    pub external_work: f64,
    /// Native total impact loss; no heat/deformation partition is inferred.
    pub contact_dissipated_energy: f64,
    /// Sum of the smooth integrator's measured numerical work defects.
    pub smooth_kinetic_work_residual: f64,
    /// Actual native stage/event kinetic increments minus work plus impact loss.
    /// Large rounded endpoint energies are not subtracted. Smooth numerical
    /// integration error remains visible in this quantity.
    pub kinetic_work_residual: f64,
    /// Actual native stage/event linear increments minus external impulse, N s.
    pub linear_momentum_residual: Vec3,
    /// Actual stage/event angular increments minus external impulse, N m s.
    pub angular_momentum_residual: Vec3,
}

/// Complete-window refusal with no partially advanced state publication.
#[derive(Debug)]
pub enum ArticulatedTrajectoryError {
    /// Invalid size, event ordering, finite-value, or representability bound.
    Invalid(&'static str),
    /// The declared fixed smooth-interval budget cannot cover this window.
    IntervalBudget {
        /// Lower bound on the required interval count.
        required_at_least: usize,
        /// Caller-declared maximum interval count.
        maximum: usize,
    },
    /// Native physical-state or smooth-step admission refused.
    Smooth {
        /// Original system index.
        system_index: usize,
        /// Attempted stage endpoint time from the window start, in seconds.
        time: f64,
        /// Unmodified original dynamics error.
        error: FreeFloatingMotionError,
    },
    /// Original prescribed-support path or smooth-step admission refused.
    PrescribedSmooth {
        /// Original system index.
        system_index: usize,
        /// Attempted stage endpoint time from the window start, in seconds.
        time: f64,
        /// Unmodified support-path or native dynamics error.
        error: crate::articulated::motion::prescribed::PrescribedBaseMotionError,
    },
    /// A native simultaneous contact event refused.
    Contact {
        /// Index in the strictly time-ordered event list.
        event_index: usize,
        /// Supplied event time from the window start, in seconds.
        time: f64,
        /// Unmodified original contact error.
        error: ArticulatedContactError,
    },
    /// Cancellation observed before the complete window could be published.
    Cancelled,
}
impl core::fmt::Display for ArticulatedTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "articulated contact trajectory: {self:?}")
    }
}
impl std::error::Error for ArticulatedTrajectoryError {}

/// Advance free articulated systems through explicit simultaneous contact events.
///
/// Native second-order smooth steps apply each held load once over its interval;
/// native impact steps change velocities at fixed poses. Internal joint forces
/// and dynamic contact pairs retain physical base recoil. An event at either
/// endpoint is allowed. Equal-time contact sets must be combined by the caller.
/// Zero events and a single smooth interval reproduce the native smooth state.
///
/// At most 128 trees, 4096 total links, 64 events, and 4096 common smooth
/// intervals are admitted. Contact work retains the native 21-contact cap and
/// contact-space solve; no generalized dense mass matrix is assembled. Returned
/// event receipts retain native determinacy and passivity evidence. All scratch
/// is private, so a later event failure or cancellation publishes no partial
/// trajectory. Cancellation is polled between systems/intervals/events and by
/// each native numerical owner.
///
/// This is an event-driven composition, not detection, constraint projection,
/// continuous contact, or a claim that the supplied event times are exact.
/// A missed contact or inadequate smooth resolution is not repaired here.
/// Angular impulses are translated from native world-origin receipts; ordinary
/// floating-point cancellation can limit extremely distant-origin balances.
/// Signed source/work sums and actual native stage/event increments are reduced
/// before rounding. Endpoint physical totals remain independent diagnostics;
/// their subtraction does not determine the complete-window residuals.
pub fn advance_articulated_events(
    systems: &[ArticulatedTrajectorySystem<'_>],
    events: &[ArticulatedTrajectoryEvent<'_>],
    config: ArticulatedTrajectoryConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedTrajectoryResult, ArticulatedTrajectoryError> {
    poll(&mut cancelled)?;
    let counts = admission(systems, events, config)?;
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
    let mut kinetic_change = Dd::ZERO;
    let mut actuator_work = Dd::ZERO;
    let mut external_work = Dd::ZERO;
    let mut loss = Dd::ZERO;
    let mut smooth_defect = Dd::ZERO;
    let mut accepted_events = Vec::with_capacity(events.len());
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
            if duration <= 0.0 || duration > config.maximum_smooth_step * (1.0 + 8.0 * f64::EPSILON)
            {
                return Err(ArticulatedTrajectoryError::Invalid(
                    "smooth interval is not representable within its step bound",
                ));
            }
            for (system_index, system) in systems.iter().enumerate() {
                let step = step_free_floating(
                    system.model,
                    &states[system_index],
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
                kinetic_change = kinetic_change + Dd::from_f64(step.kinetic_energy_change);
                actuator_work = actuator_work + Dd::from_f64(step.actuator_work);
                external_work = external_work + Dd::from_f64(step.external_work);
                smooth_defect = smooth_defect + Dd::from_f64(step.kinetic_work_residual);
                states[system_index] = step.state_after;
            }
            previous_time = next_time;
        }
        start = end;
        if let Some(event) = events.get(interval_index) {
            poll(&mut cancelled)?;
            let inputs: Vec<_> = systems
                .iter()
                .zip(&states)
                .map(|(system, state)| ArticulatedContactSystem {
                    model: system.model,
                    base: state.base,
                    joint_position: &state.joint_position,
                    joint_velocity: &state.joint_velocity,
                })
                .collect();
            let impact = resolve_articulated_frictional_contacts(
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
                kinetic_change = kinetic_change + Dd::from_f64(response.kinetic_energy_change);
                state.base = response.base_after;
                state
                    .joint_velocity
                    .clone_from(&response.joint_velocity_after);
            }
            external_work = external_work + Dd::from_f64(impact.plane_work);
            loss = loss + Dd::from_f64(impact.contact.dissipated_energy);
            accepted_events.push(ArticulatedTrajectoryEventResult {
                time: event.time,
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
    let kinetic_work_residual = (kinetic_change - actuator_work - external_work + loss).to_f64();
    let (linear_momentum_residual, angular_momentum_residual) =
        momentum_change.difference(&impulses);
    let [actuator_work, external_work, loss, smooth_defect] =
        [actuator_work, external_work, loss, smooth_defect].map(Dd::to_f64);
    if ![
        actuator_work,
        external_work,
        loss,
        smooth_defect,
        kinetic_work_residual,
    ]
    .iter()
    .all(|value| value.is_finite())
        || !finite_vector(impulses.linear)
        || !finite_vector(impulses.angular)
        || !finite_vector(momentum_change.linear)
        || !finite_vector(momentum_change.angular)
        || !finite_vector(linear_momentum_residual)
        || !finite_vector(angular_momentum_residual)
    {
        return Err(ArticulatedTrajectoryError::Invalid(
            "nonfinite trajectory balance",
        ));
    }
    poll(&mut cancelled)?;
    Ok(ArticulatedTrajectoryResult {
        states_after: states,
        events: accepted_events,
        smooth_intervals: counts.iter().sum(),
        angular_origin_world: config.angular_origin_world,
        momentum_before,
        momentum_after,
        external_linear_impulse_world: impulses.linear,
        external_angular_impulse_about_origin_world: impulses.angular,
        actuator_work,
        external_work,
        contact_dissipated_energy: loss,
        smooth_kinetic_work_residual: smooth_defect,
        kinetic_work_residual,
        linear_momentum_residual,
        angular_momentum_residual,
    })
}

fn admission(
    systems: &[ArticulatedTrajectorySystem<'_>],
    events: &[ArticulatedTrajectoryEvent<'_>],
    config: ArticulatedTrajectoryConfig,
) -> Result<Vec<usize>, ArticulatedTrajectoryError> {
    let counts = event_schedule(systems.len(), events, config)?;
    let mut links = 0usize;
    for system in systems {
        if system.state.joint_position.len() != system.model.dof_count()
            || system.state.joint_velocity.len() != system.model.dof_count()
        {
            return Err(ArticulatedTrajectoryError::Invalid(
                "joint state lengths must match the original model DOF count",
            ));
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
    Ok(counts)
}

fn event_schedule(
    system_count: usize,
    events: &[ArticulatedTrajectoryEvent<'_>],
    config: ArticulatedTrajectoryConfig,
) -> Result<Vec<usize>, ArticulatedTrajectoryError> {
    if system_count == 0
        || system_count > 128
        || events.len() > 64
        || !(1..=4096).contains(&config.maximum_smooth_intervals)
        || !config.duration.is_finite()
        || config.duration <= 0.0
        || !config.maximum_smooth_step.is_finite()
        || config.maximum_smooth_step <= 0.0
        || !finite_vector(config.angular_origin_world)
    {
        return Err(ArticulatedTrajectoryError::Invalid(
            "invalid trajectory bounds",
        ));
    }
    let mut counts = Vec::with_capacity(events.len() + 1);
    let mut previous = 0.0;
    let mut total = 0usize;
    for index in 0..=events.len() {
        let time = if let Some(event) = events.get(index) {
            if !event.time.is_finite()
                || event.time < 0.0
                || event.time > config.duration
                || (index != 0 && event.time <= previous)
                || event.contacts.is_empty()
            {
                return Err(ArticulatedTrajectoryError::Invalid(
                    "invalid or non-increasing contact event",
                ));
            }
            event.time
        } else {
            config.duration
        };
        let needed = if time > previous {
            ((time - previous) / config.maximum_smooth_step)
                .ceil()
                .max(1.0)
        } else {
            0.0
        };
        if !needed.is_finite() || needed > 4096.0 {
            return Err(ArticulatedTrajectoryError::IntervalBudget {
                required_at_least: 4097,
                maximum: config.maximum_smooth_intervals,
            });
        }
        let count = needed as usize;
        total += count;
        if total > config.maximum_smooth_intervals {
            return Err(ArticulatedTrajectoryError::IntervalBudget {
                required_at_least: total,
                maximum: config.maximum_smooth_intervals,
            });
        }
        counts.push(count);
        previous = time;
    }
    Ok(counts)
}

fn momentum(
    systems: &[ArticulatedTrajectorySystem<'_>],
    states: &[FreeFloatingMotionState],
    origin: Vec3,
    time: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ArticulatedTrajectoryMomentum, ArticulatedTrajectoryError> {
    let mut spatial = Impulses::default();
    let mut kinetic = Dd::ZERO;
    for (index, (system, state)) in systems.iter().zip(states).enumerate() {
        poll(cancelled)?;
        // Move the coordinate origin before evaluating native link momenta;
        // avoid subtracting two large whole-system orbital angular momenta.
        let mut base = state.base;
        base.world_from_base = Se3::from_parts(
            base.world_from_base.rotation(),
            base.world_from_base.translation() - origin,
        )
        .map_err(|error| smooth_error(index, time, error.into()))?;
        let value: ArticulatedMomentum = free_floating_momentum(
            system.model,
            base,
            &state.joint_position,
            &state.joint_velocity,
        )
        .map_err(|error| smooth_error(index, time, error.into()))?;
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
            "nonfinite total momentum",
        ));
    }
    Ok(total)
}

#[derive(Default)]
struct Impulses {
    linear: Vec3,
    angular: Vec3,
    linear_sum: [Dd; 3],
    angular_sum: [Dd; 3],
}
impl Impulses {
    fn add(&mut self, linear: Vec3, angular: Vec3, origin: Vec3) {
        let p = [linear.x, linear.y, linear.z].map(Dd::from_f64);
        let l = [angular.x, angular.y, angular.z].map(Dd::from_f64);
        let x = [origin.x, origin.y, origin.z].map(Dd::from_f64);
        for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
            self.linear_sum[i] = self.linear_sum[i] + p[i];
            // Retain the original angular source and each origin-shift term;
            // their cancellation must not erase an independent small couple.
            self.angular_sum[i] = self.angular_sum[i] + l[i] - x[j] * p[k] + x[k] * p[j];
        }
        self.linear = Self::vector(self.linear_sum);
        self.angular = Self::vector(self.angular_sum);
    }

    fn difference(&self, other: &Self) -> (Vec3, Vec3) {
        (
            Self::vector(core::array::from_fn(|i| {
                self.linear_sum[i] - other.linear_sum[i]
            })),
            Self::vector(core::array::from_fn(|i| {
                self.angular_sum[i] - other.angular_sum[i]
            })),
        )
    }

    fn vector(value: [Dd; 3]) -> Vec3 {
        Vec3::new(value[0].to_f64(), value[1].to_f64(), value[2].to_f64())
    }
}
fn finite_vector(value: Vec3) -> bool {
    value.x.is_finite() && value.y.is_finite() && value.z.is_finite()
}
fn smooth_error(
    system_index: usize,
    time: f64,
    error: FreeFloatingMotionError,
) -> ArticulatedTrajectoryError {
    if matches!(error, FreeFloatingMotionError::Cancelled) {
        ArticulatedTrajectoryError::Cancelled
    } else {
        ArticulatedTrajectoryError::Smooth {
            system_index,
            time,
            error,
        }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ArticulatedTrajectoryError> {
    if cancelled() {
        Err(ArticulatedTrajectoryError::Cancelled)
    } else {
        Ok(())
    }
}
