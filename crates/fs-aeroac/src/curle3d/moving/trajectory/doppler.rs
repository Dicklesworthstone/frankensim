//! Directional rates on the original fixed-microphone admission leaves.
//!
//! Preparation caches one upper bound per ORIGINAL source span. Finite-window
//! queries reuse the monotone outward causal search, so unrelated future spans
//! do not impose their Doppler expansion. A partially overlapping span retains
//! its complete continuous bound; no endpoint-only velocity shortcut is used.

use super::receiver::causal::{knot_residual, source_span_range};
use super::receiver::{distance_hull, doppler::directional_maximum_rate};
use super::{BoundNode, MovingCompactObserver, PreparedSourceTrajectory, invalid};
use crate::AeroacError;
use fs_ivl::Interval;

pub(super) fn independent_maximum_rate(mach: f64) -> f64 {
    if mach == 0.0 {
        1.0
    } else {
        let one = Interval::point(1.0);
        (one / (one - Interval::point(mach))).hi()
    }
}

pub(super) fn maximum_rate(
    observer: &MovingCompactObserver,
    origin: [f64; 3],
    node: &BoundNode,
    maximum_mach: f64,
) -> f64 {
    if maximum_mach == 0.0 {
        return 1.0;
    }
    let c = Interval::point(observer.medium().sound_speed_m_s());
    let displacement = core::array::from_fn(|axis| {
        let mut hull = node.position[0][axis];
        for control in &node.position[1..] {
            hull = hull.hull(control[axis]);
        }
        Interval::point(observer.observer_position_m()[axis])
            - Interval::point(origin[axis]) - hull
    });
    let velocity = core::array::from_fn(|axis| {
        let mut hull = node.velocity[0][axis];
        for control in &node.velocity[1..] {
            hull = hull.hull(control[axis]);
        }
        hull / c
    });
    directional_maximum_rate(
        displacement,
        distance_hull(displacement),
        velocity,
        [Interval::point(0.0); 3],
        maximum_mach,
        0.0,
        independent_maximum_rate(maximum_mach),
    )
}

pub(super) fn admit_window(
    window: [f64; 2],
    coverage: [f64; 2],
    tolerance: f64,
) -> Result<(), AeroacError> {
    if window.iter().any(|time| !time.is_finite())
        || window[0] > window[1]
        || window[0] < coverage[0]
        || window[1] > coverage[1]
        || !tolerance.is_finite()
        || tolerance < 0.0
    {
        return Err(invalid(
            "Doppler window exceeds finite receiver coverage or residual domain",
        ));
    }
    Ok(())
}

pub(super) fn source_spans(
    trajectory: &PreparedSourceTrajectory,
    observer: &MovingCompactObserver,
    window: [f64; 2],
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<core::ops::Range<usize>, AeroacError> {
    let spans = &trajectory.segments;
    let position = observer.observer_position_m().map(Interval::point);
    let c = Interval::point(observer.medium().sound_speed_m_s());
    source_span_range(
        spans.len(),
        tolerance,
        |mid| knot_residual(spans[mid].polynomial.last, position, window[0], c),
        |mid| knot_residual(spans[mid].polynomial.first, position, window[1], c),
        cancelled,
    )
}
