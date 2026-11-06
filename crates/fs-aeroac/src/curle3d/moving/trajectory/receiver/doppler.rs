//! Continuous Doppler admission using both accepted velocity histories.
//!
//! At a causal event, a = v_source(tau)/c, b = v_receiver(t)/c and |n| = 1:
//!   rate = (1 - n.b)/(1 - n.a)
//!        = 1 + n.(a-b)/(1 - n.a)
//!       <= 1 + sup|a-b|/(1 - M_source).
//! Both velocity boxes enclose their COMPLETE native interpolants, including
//! all accepted Bernstein subdivisions and exact knots. They are independent
//! event-time ranges: equal simultaneous velocities or equal endpoint samples
//! do NOT imply equal velocities at emission and reception when accelerating.
//!
//! Intersect this upper bound with (1+M_receiver)/(1-M_source). Both bounds
//! are outward-rounded and valid in every direction; either can be tighter.
//! This removes the false opposing-motion expansion for common translation,
//! including fixed laboratory endpoints in uniform wind. It never relaxes
//! source sampling, receiver Nyquist, clearance, root or amplitude-modulation
//! policies. Original pressure coefficients and retarded roots are unchanged.

use super::{PreparedKinematicTrajectory, ReceiverImage};
use fs_ivl::Interval;

pub(super) fn maximum_rate(
    source: &PreparedKinematicTrajectory,
    receiver: &PreparedKinematicTrajectory,
    image: Option<ReceiverImage<'_>>,
) -> f64 {
    if source.maximum_mach() == 0.0 && receiver.maximum_mach() == 0.0 {
        return 1.0;
    }
    let one = Interval::point(1.0);
    let denominator = one - Interval::point(source.maximum_mach());
    let receiver_mach = image.map_or(receiver.maximum_mach(), |image| {
        image.maximum_mach(receiver.maximum_mach())
    });
    if !receiver_mach.is_finite() || receiver_mach >= 1.0 {
        return f64::NAN;
    }
    let independent = ((one + Interval::point(receiver_mach)) / denominator).hi();
    let receiver_velocity = image.map_or(receiver.velocity_mach_bounds, |image| {
        image.vector_bounds(receiver.velocity_mach_bounds)
    });
    // Dimensionless components avoid forming large physical velocity
    // differences or squaring sound speeds. Each path is already subsonic
    // in the SAME medium/frame, checked by the original pair constructor.
    let mut square = Interval::point(0.0);
    for axis in 0..3 {
        let difference = (source.velocity_mach_bounds[axis] - receiver_velocity[axis]).abs();
        square = square + difference * difference;
    }
    let relative = (one + square.sqrt() / denominator).hi();
    if relative.is_finite() && relative >= 1.0 {
        independent.min(relative)
    } else {
        // A loose or unrepresentable auxiliary bound cannot make a previously
        // admissible trajectory fail or waive the original speed-only bound.
        independent
    }
}

// A directional refinement on ONE causally possible clearance leaf. The
// velocity hulls enclose the complete original source/receiver spans, not
// simultaneous samples. Geometry and velocity remain independent enclosures.
pub(super) fn directional_maximum_rate(
    displacement: [Interval; 3],
    distance: Interval,
    source_velocity: [Interval; 3],
    receiver_velocity: [Interval; 3],
    source_mach: f64,
    receiver_mach: f64,
    fallback: f64,
) -> f64 {
    if distance.lo() <= 0.0 || !distance.hi().is_finite() {
        return fallback;
    }
    let Some(direction) = direction_bounds(displacement, distance) else {
        return fallback;
    };
    let one = Interval::point(1.0);
    let source_speed = Interval::new(-source_mach, source_mach);
    let receiver_speed = Interval::new(-receiver_mach, receiver_mach);
    let mut source_radial = Interval::point(0.0);
    let mut receiver_radial = Interval::point(0.0);
    let mut relative_radial = Interval::point(0.0);
    for axis in 0..3 {
        source_radial = source_radial + direction[axis] * source_velocity[axis];
        receiver_radial = receiver_radial + direction[axis] * receiver_velocity[axis];
        relative_radial = relative_radial
            + direction[axis] * (source_velocity[axis] - receiver_velocity[axis]);
    }
    let Some(source_radial) = source_radial.intersect(source_speed) else {
        return fallback;
    };
    let Some(receiver_radial) = receiver_radial.intersect(receiver_speed) else {
        return fallback;
    };
    let denominator = one - source_radial;
    if denominator.lo() <= 0.0 {
        return fallback;
    }
    // The signed relative form preserves cancellation for common motion.
    // Unlike the direction-free norm bound it also preserves recession.
    let direct = ((one - receiver_radial) / denominator).hi();
    let relative = (one + relative_radial / denominator).hi();
    let mut upper = fallback;
    for candidate in [direct, relative] {
        if candidate.is_finite() && candidate > 0.0 {
            upper = upper.min(candidate);
        }
    }
    upper
}

fn direction_bounds(displacement: [Interval; 3], distance: Interval) -> Option<[Interval; 3]> {
    let unit = Interval::new(-1.0, 1.0);
    let one = Interval::point(1.0);
    let mut direction = [unit; 3];
    for axis in 0..3 {
        direction[axis] = (displacement[axis] / distance).intersect(unit)?;
        if displacement[axis].contains_zero() {
            continue;
        }
        // n_i = sign(d_i)/sqrt(1 + sum_{j != i}(d_j/d_i)^2).
        // This retains the correlation between axial distance and its norm:
        // an axial receding path has n_i = +/-1 even over a long time span.
        let mut square = one;
        for other in 0..3 {
            if other != axis {
                let ratio = (displacement[other] / displacement[axis]).abs();
                square = square + ratio * ratio;
            }
        }
        let magnitude = one / square.sqrt();
        let signed = if displacement[axis].lo() > 0.0 {
            magnitude
        } else {
            -magnitude
        };
        direction[axis] = direction[axis].intersect(signed)?;
    }
    Some(direction)
}
