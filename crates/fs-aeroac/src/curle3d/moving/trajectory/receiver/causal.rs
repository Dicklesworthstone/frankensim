//! Exclude whole source-span ranges before expensive pair-hull admission.
//!
//! For strictly subsonic paths, F(tau,t) = tau + |x(t)-y(tau)|/c - t
//! increases with tau and decreases with t. Thus a source knot whose residual
//! is below -tolerance at the first receiver time excludes EVERY earlier
//! emission; a knot above +tolerance at the last time excludes EVERY later
//! emission. Outward intervals, not approximate solved roots, decide exclusion.
//! An inconclusive sign retains extra spans. No physical path is extrapolated.

use super::super::path::KinematicSegment;
use super::super::TimedSourceKinematics;
use super::{PreparedMovingReceiverSource, distance_hull, poll};
use crate::AeroacError;
use fs_ivl::Interval;

pub(super) fn source_spans(
    owner: &PreparedMovingReceiverSource<'_>,
    receiver: &KinematicSegment,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<core::ops::Range<usize>, AeroacError> {
    let r = &receiver.polynomial;
    let first_time = owner.config.receiver_time_interval_s[0].max(r.first.source_time_s);
    let last_time = owner.config.receiver_time_interval_s[1].min(r.last.source_time_s);
    let transform = |time| {
        let position = position_at(receiver, time);
        owner.image.map_or(position, |image| image.point_bounds(position))
    };
    let first_position = transform(first_time);
    let last_position = transform(last_time);
    let c = Interval::point(owner.source.medium().sound_speed_m_s());
    let tolerance = owner.config.maximum_time_residual_s;
    let spans = &owner.source.segments;
    source_span_range(
        spans.len(),
        tolerance,
        |mid| knot_residual(spans[mid].polynomial.last, first_position, first_time, c),
        |mid| knot_residual(spans[mid].polynomial.first, last_position, last_time, c),
        cancelled,
    )
}

// The endpoint residual providers preserve each owner's original receiver
// geometry. Only decisive outward signs remove a prefix or suffix; uncertain
// signs need not themselves form a monotone predicate.
pub(in crate::curle3d::moving::trajectory) fn source_span_range(
    count: usize,
    tolerance: f64,
    mut last_knot_at_first_receiver: impl FnMut(usize) -> Interval,
    mut first_knot_at_last_receiver: impl FnMut(usize) -> Interval,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<core::ops::Range<usize>, AeroacError> {
    let mut lo = 0;
    let mut hi = count;
    while lo < hi {
        poll(cancelled)?;
        let mid = lo + (hi - lo) / 2;
        if last_knot_at_first_receiver(mid).hi() < -tolerance {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let first = lo;
    hi = count;
    while lo < hi {
        poll(cancelled)?;
        let mid = lo + (hi - lo) / 2;
        if first_knot_at_last_receiver(mid).lo() > tolerance {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    Ok(first..lo)
}

pub(in crate::curle3d::moving::trajectory) fn knot_residual(
    source: TimedSourceKinematics,
    receiver_position: [Interval; 3],
    receiver_time: f64,
    sound_speed: Interval,
) -> Interval {
    let displacement = core::array::from_fn(|axis| {
        receiver_position[axis] - Interval::point(source.kinematics.position_m[axis])
    });
    Interval::point(source.source_time_s) - Interval::point(receiver_time)
        + distance_hull(displacement) / sound_speed
}

// Intersected receiver windows can begin/end INSIDE an original span. Evaluate
// that same Bernstein polynomial in interval arithmetic; endpoint samples and
// a rounded Newton solution are not continuous bounds on its interior.
fn position_at(segment: &KinematicSegment, time: f64) -> [Interval; 3] {
    let p = &segment.polynomial;
    if time == p.first.source_time_s {
        return p.first.kinematics.position_m.map(Interval::point);
    }
    if time == p.last.source_time_s {
        return p.last.kinematics.position_m.map(Interval::point);
    }
    let unit = Interval::new(0.0, 1.0);
    let duration = Interval::point(p.last.source_time_s) - Interval::point(p.first.source_time_s);
    let fraction = (Interval::point(time) - Interval::point(p.first.source_time_s)) / duration;
    let u = fraction.intersect(unit).unwrap_or(unit);
    let complement = Interval::point(1.0) - u;
    let mut controls = segment.position_bounds;
    for width in (1..6).rev() {
        for index in 0..width {
            for axis in 0..3 {
                controls[index][axis] = complement * controls[index][axis]
                    + u * controls[index + 1][axis];
            }
        }
    }
    core::array::from_fn(|axis| {
        Interval::point(p.first.kinematics.position_m[axis]) + controls[0][axis]
    })
}
