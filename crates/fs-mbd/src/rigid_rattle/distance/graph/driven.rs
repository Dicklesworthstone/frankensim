//! Prescribed positive-length actuators on the original simultaneous graph.

use super::{
    ActionReactionImpulseReceipt, Dd, DistanceGraphControls, DistanceGraphError, DistanceGraphStep,
    DistanceLink, DistanceStepError, MAX_DISTANCE_ROWS, MassProperties, RigidBodyIntegrator,
    RigidBodyState, Wrench, advance_between, limit, planned_bytes, poll, positive, reserved,
};
use crate::{PointKinematics, Vec3};

/// Held range acceleration of one bilateral actuator for the accepted interval.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct DistanceDrive {
    /// Initial signed rate of increasing endpoint separation [m/s].
    pub range_rate_mps: f64,
    /// Constant signed range acceleration [m/s²].
    pub range_acceleration_mps2: f64,
}

/// Continuation and work of one actuator, in the original link order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DrivenDistanceLinkStep {
    /// Original endpoints with the represented final commanded length [m].
    pub link_after: DistanceLink,
    /// Final range rate and the unchanged held acceleration for continuation.
    pub drive_after: DistanceDrive,
    /// Represented endpoint length minus the original polynomial value [m].
    pub length_representation_residual_m: f64,
    /// Represented endpoint range rate minus the original polynomial rate [m/s].
    pub rate_representation_residual_mps: f64,
    /// Signed work delivered by this actuator [J], `j_start*r_start + j_end*r_end`.
    pub actuator_work_j: f64,
    /// Actual ordered pair work minus actuator work [J], evaluated directly
    /// from relative point velocities. Cross terms follow native link order.
    pub constraint_work_defect_j: f64,
}

/// Atomic prescribed-length evolution, including finite-body recoil and work.
#[derive(Debug, Clone, PartialEq)]
pub struct DrivenDistanceGraphStep {
    /// Original simultaneous native step, with residuals against the commanded
    /// endpoint polynomial. Its payload bound describes the native solve alone.
    pub step: DistanceGraphStep,
    /// Updated positive lengths, rates, representation errors and per-link work.
    pub links: Vec<DrivenDistanceLinkStep>,
    /// Sum of signed actuator work, accumulated before rounding [J].
    pub actuator_work_j: f64,
    /// Sum of the independently paired numerical constraint-work defect [J].
    /// This is not dissipation or heat and is never projected out of the bodies.
    pub constraint_work_defect_j: f64,
    /// Conservative owned-payload bound including the continuation receipts.
    pub planned_live_bytes: usize,
}

/// Drive simultaneous internal rods with `length(t)=length(0)+rate*t+a*t²/2`.
///
/// Both bodies retain their actual inertia and recoil. Initial/final central
/// impulses are solved together by the original graph owner; they conserve the
/// same total world momenta as static rods, while the prescribed actuators may
/// supply or recover energy. Each event's scalar impulse is paired with the
/// range rate at that event. Work defect is evaluated in relative velocities,
/// without subtracting large laboratory-frame kinetic energies.
///
/// The prescribed quadratic must stay above the position tolerance throughout
/// the interval. Endpoint length/rate representation errors must fit their
/// physical tolerances; the solver retains those errors against the original
/// polynomial. Nonzero displacement that rounds to unchanged length refuses.
/// Use `links[*].link_after` and `drive_after` for subsequent accepted intervals.
///
/// All native graph resource, regularity and local angular/impulse bounds still
/// apply. The static graph API is unchanged. Driven rate admission and work use
/// the separate native COM and spin components even for zero motion commands;
/// a boost cannot hide spin by rounding the convenience point velocity first.
/// This is an ideal bilateral actuator with an external energy supply; no
/// motor/storage/efficiency model, unilateral cable, impact, exact energy law,
/// global branch enclosure or arbitrary time callback is implied.
#[allow(clippy::too_many_arguments)]
pub fn advance_driven_distance_graph(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[DistanceLink],
    drives: &[DistanceDrive],
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    duration_s: f64,
    controls: DistanceGraphControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DrivenDistanceGraphStep, DistanceGraphError> {
    poll(&mut cancelled)?;
    if links.len() != drives.len()
        || !positive(duration_s)
        || !positive(controls.position_tolerance_m)
        || !positive(controls.velocity_tolerance_mps)
    {
        return Err(DistanceStepError::Invalid("prescribed graph shape or controls").into());
    }
    limit(
        "bodies",
        states.len() as u128,
        controls.maximum_bodies as u128,
    )?;
    limit(
        "constraints",
        links.len() as u128,
        controls.maximum_constraints.min(MAX_DISTANCE_ROWS) as u128,
    )?;
    let extra = links.len() as u128 * core::mem::size_of::<DrivenDistanceLinkStep>() as u128;
    let bytes = planned_bytes(states.len(), links.len()) + extra;
    limit("live bytes", bytes, controls.maximum_live_bytes as u128)?;
    // The existing hard row bound also bounds this allocation-free command bank.
    let mut targets = [Target::default(); MAX_DISTANCE_ROWS];
    for (i, (link, drive)) in links.iter().zip(drives).enumerate() {
        poll(&mut cancelled)?;
        targets[i] = Target::prepare(i, *link, *drive, duration_s, controls)?;
    }
    let mut native_controls = controls;
    native_controls.maximum_live_bytes -= extra as usize;
    let (step, work) = advance_between(
        states,
        properties,
        links,
        integrator,
        wrenches,
        duration_s,
        native_controls,
        Some(&targets[..links.len()]),
        &mut cancelled,
    )?;
    let mut continuation = reserved(links.len())?;
    let mut source = Dd::ZERO;
    let mut defect = Dd::ZERO;
    for (i, (link, drive)) in links.iter().zip(drives).enumerate() {
        poll(&mut cancelled)?;
        let target = targets[i];
        let account = work.links[i];
        source = source + account.source;
        defect = defect + account.defect;
        continuation.push(DrivenDistanceLinkStep {
            link_after: DistanceLink {
                length_m: target.length,
                ..*link
            },
            drive_after: DistanceDrive {
                range_rate_mps: target.final_rate.to_f64(),
                ..*drive
            },
            length_representation_residual_m: target.length_residual.to_f64(),
            rate_representation_residual_mps: (Dd::from_f64(target.final_rate.to_f64())
                - target.final_rate)
                .to_f64(),
            actuator_work_j: account.source.to_f64(),
            constraint_work_defect_j: account.defect.to_f64(),
        });
    }
    if !finite(source) || !finite(defect) {
        return Err(DistanceStepError::Invalid("prescribed graph work sum").into());
    }
    poll(&mut cancelled)?;
    Ok(DrivenDistanceGraphStep {
        step,
        links: continuation,
        actuator_work_j: source.to_f64(),
        constraint_work_defect_j: defect.to_f64(),
        planned_live_bytes: bytes as usize,
    })
}

#[derive(Clone, Copy, Default)]
pub(super) struct Target {
    pub length: f64,
    pub length_residual: Dd,
    pub initial_rate: f64,
    pub final_rate: Dd,
}
impl Target {
    fn prepare(
        index: usize,
        link: DistanceLink,
        drive: DistanceDrive,
        duration: f64,
        controls: DistanceGraphControls,
    ) -> Result<Self, DistanceGraphError> {
        let error = |reason| DistanceGraphError::PrescribedMotion {
            constraint: index,
            reason,
        };
        if !positive(link.length_m)
            || !drive.range_rate_mps.is_finite()
            || !drive.range_acceleration_mps2.is_finite()
        {
            return Err(error("finite positive length and finite motion required"));
        }
        let h = Dd::from_f64(duration);
        let start = Dd::from_f64(link.length_m);
        let rate = Dd::from_f64(drive.range_rate_mps);
        let acceleration = Dd::from_f64(drive.range_acceleration_mps2);
        let delta_rate = acceleration * h;
        let final_rate = rate + delta_rate;
        let delta_length = h * (rate + Dd::from_f64(0.5) * delta_rate);
        let final_length = start + delta_length;
        let length = final_length.to_f64();
        let length_residual = Dd::from_f64(length) - final_length;
        let rate_residual = Dd::from_f64(final_rate.to_f64()) - final_rate;
        if !finite(final_length)
            || !finite(final_rate)
            || !finite(length_residual)
            || !finite(rate_residual)
            || length_residual.abs().to_f64() > controls.position_tolerance_m
            || rate_residual.abs().to_f64() > controls.velocity_tolerance_mps
            || (delta_length != Dd::ZERO && length == link.length_m)
        {
            return Err(error("endpoint prescription is unrepresentable"));
        }
        let tolerance = Dd::from_f64(controls.position_tolerance_m);
        if !tolerance.lt(start)
            || !tolerance.lt(final_length)
            || length <= controls.position_tolerance_m
        {
            return Err(error("prescribed length reaches the position tolerance"));
        }
        // Endpoint positivity alone misses a contracting-then-extending path.
        if drive.range_acceleration_mps2 > 0.0
            && drive.range_rate_mps < 0.0
            && Dd::ZERO.lt(final_rate)
        {
            let turn = -rate / acceleration;
            let minimum = start + turn * (rate + Dd::from_f64(0.5) * acceleration * turn);
            if !finite(minimum) || !tolerance.lt(minimum) {
                return Err(error("prescribed length collapses inside the interval"));
            }
        }
        Ok(Self {
            length,
            length_residual,
            initial_rate: drive.range_rate_mps,
            final_rate,
        })
    }
}

#[derive(Clone, Copy, Default)]
pub(super) struct LinkWork {
    pub source: Dd,
    pub defect: Dd,
}
pub(super) struct Work {
    pub links: [LinkWork; MAX_DISTANCE_ROWS],
}
impl Default for Work {
    fn default() -> Self {
        Self {
            links: [LinkWork::default(); MAX_DISTANCE_ROWS],
        }
    }
}
impl Work {
    pub fn add(
        &mut self,
        events: &[ActionReactionImpulseReceipt],
        impulses: &[f64],
        targets: &[Target],
        endpoint: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DistanceGraphError> {
        for (i, ((event, impulse), target)) in events.iter().zip(impulses).zip(targets).enumerate()
        {
            poll(cancelled)?;
            let rate = if endpoint {
                target.final_rate
            } else {
                Dd::from_f64(target.initial_rate)
            };
            let source = Dd::from_f64(*impulse) * rate;
            let before = relative_velocity(
                event.first.point_kinematics_before,
                event.second.point_kinematics_before,
            );
            let after = relative_velocity(
                event.first.point_kinematics_after,
                event.second.point_kinematics_after,
            );
            let j = event.second.impulse_world;
            let mut relative_work = Dd::ZERO;
            for ((before, after), j) in before.into_iter().zip(after).zip([j.x, j.y, j.z]) {
                relative_work = relative_work
                    + (before * Dd::from_f64(0.5) + after * Dd::from_f64(0.5)) * Dd::from_f64(j);
            }
            let work = &mut self.links[i];
            work.source = work.source + source;
            work.defect = work.defect + (relative_work - source);
            if !finite(work.source) || !finite(work.defect) {
                return Err(DistanceStepError::Invalid("prescribed graph impulse work").into());
            }
        }
        Ok(())
    }
}

// Cancel common translation before adding local rotational motion. Subtracting
// the convenience point velocities would already have discarded small spin
// beside a large common boost, before any double-double pairing could save it.
fn relative_velocity(first: PointKinematics, second: PointKinematics) -> [Dd; 3] {
    let components = |v: Vec3| [v.x, v.y, v.z].map(Dd::from_f64);
    let a = components(first.center_of_mass_velocity_world);
    let b = components(second.center_of_mass_velocity_world);
    let wa = components(first.angular_velocity_world);
    let wb = components(second.angular_velocity_world);
    let ra = components(first.arm_world);
    let rb = components(second.arm_world);
    core::array::from_fn(|k| {
        let i = (k + 1) % 3;
        let j = (k + 2) % 3;
        (b[k] - a[k]) + (wb[i] * rb[j] - wb[j] * rb[i]) - (wa[i] * ra[j] - wa[j] * ra[i])
    })
}
pub(super) fn range_rate(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    points: [Vec3; 2],
    normal: Vec3,
) -> Result<Dd, DistanceStepError> {
    let velocity = relative_velocity(
        states[0].point_kinematics(properties[0], points[0])?,
        states[1].point_kinematics(properties[1], points[1])?,
    );
    let result = velocity
        .into_iter()
        .zip([normal.x, normal.y, normal.z])
        .fold(Dd::ZERO, |rate, (v, n)| rate + v * Dd::from_f64(n));
    if !finite(result) {
        return Err(DistanceStepError::Invalid("prescribed native range rate"));
    }
    Ok(result)
}
fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
