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

use super::PreparedKinematicTrajectory;
use fs_ivl::Interval;

pub(super) fn maximum_rate(
    source: &PreparedKinematicTrajectory,
    receiver: &PreparedKinematicTrajectory,
) -> f64 {
    if source.maximum_mach() == 0.0 && receiver.maximum_mach() == 0.0 {
        return 1.0;
    }
    let one = Interval::point(1.0);
    let denominator = one - Interval::point(source.maximum_mach());
    let independent = ((one + Interval::point(receiver.maximum_mach())) / denominator).hi();
    // Dimensionless components avoid forming large physical velocity
    // differences or squaring sound speeds. Each path is already subsonic
    // in the SAME medium/frame, checked by the original pair constructor.
    let mut square = Interval::point(0.0);
    for axis in 0..3 {
        let difference = (source.velocity_mach_bounds[axis]
            - receiver.velocity_mach_bounds[axis])
            .abs();
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
