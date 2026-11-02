//! Entire convolution of a linear rotating-quadrature envelope.
use super::{DampedTransition, HarmonicForceTransition, double_force, invalid, small_force};
use crate::time::ModalAcousticTimeError;
use fs_math::det;

// Gains for (t/h)*(Dc*cos(Omega*t)+Ds*sin(Omega*t)), in the same
// [q/h^2,q/h^2,v/h,v/h] convention as the constant harmonic owner.
pub(super) fn delta_gains(
    base: &HarmonicForceTransition,
    frequency: f64,
) -> Result<[f64; 4], ModalAcousticTimeError> {
    let mut angle = base.transition.angle;
    let mut drag = base.transition.drag;
    let mut phase = frequency * base.duration_s;
    let mut doublings = 0;
    while angle.max(drag).max(phase) > 0.25 {
        angle *= 0.5;
        drag *= 0.5;
        phase *= 0.5;
        doublings += 1;
    }
    let mut flow = DampedTransition::small(angle, drag);
    let mut force = small_force(angle, drag, phase, false);
    let mut delta = small_delta(angle, drag, phase);
    for _ in 0..doublings {
        delta = double_delta(delta, force, flow, phase);
        force = double_force(force, flow, phase, false);
        flow.double();
        phase *= 2.0;
    }
    if delta.into_iter().any(|value| !value.is_finite()) {
        return Err(invalid("affine tone convolution must be representable"));
    }
    Ok(delta)
}

pub(super) fn small_delta(angle: f64, drag: f64, phase: f64) -> [f64; 4] {
    let mut gains = [0.0; 4];
    for column in 0..2 {
        // Dimensionless entire action on [q,h*v,h^2*C,h^2*S,h^2*Dc,h^2*Ds].
        // C'=phase*S+Dc and S'=-phase*C+Ds; the deltas also rotate.
        let mut term = [0.0; 6];
        term[4 + column] = 1.0;
        for order in 1..=24 {
            let inverse = 1.0 / f64::from(order);
            term = [
                term[1] * inverse,
                (-angle * angle * term[0] - drag * term[1] + term[2]) * inverse,
                (phase * term[3] + term[4]) * inverse,
                (-phase * term[2] + term[5]) * inverse,
                phase * term[5] * inverse,
                -phase * term[4] * inverse,
            ];
            gains[column] += term[0];
            gains[2 + column] += term[1];
        }
    }
    gains
}

pub(super) fn double_delta(
    delta: [f64; 4],
    force: [f64; 4],
    flow: DampedTransition,
    phase: f64,
) -> [f64; 4] {
    let (sine, cosine) = det::sin_cos(phase);
    let [qc, qs, vc, vs] = delta;
    let restoring = -(flow.angle * flow.b) * flow.angle;
    // On each half the envelope increment is Delta/2. At the second
    // half's start it also has the constant amplitude R*Delta/2.
    let (cq, sq) = (qc + force[0], qs + force[1]);
    let (cv, sv) = (vc + force[2], vs + force[3]);
    [
        0.125 * (flow.position_retention * qc + flow.b * vc + cq * cosine - sq * sine),
        0.125 * (flow.position_retention * qs + flow.b * vs + cq * sine + sq * cosine),
        0.25 * (restoring * qc + flow.velocity_retention * vc + cv * cosine - sv * sine),
        0.25 * (restoring * qs + flow.velocity_retention * vs + cv * sine + sv * cosine),
    ]
}
