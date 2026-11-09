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
    delta_gains_for_phase(base, frequency * base.duration_s)
}

// Sensitivities retain the original unwrapped source phase, without recovering
// a frequency by division only to multiply it by the same duration again.
pub(super) fn delta_gains_for_phase(
    base: &HarmonicForceTransition,
    mut phase: f64,
) -> Result<[f64; 4], ModalAcousticTimeError> {
    let mut angle = base.transition.angle;
    let mut drag = base.transition.drag;
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

/// Complementary endpoint moments used by affine-envelope pullbacks.
/// `u=t/h`: start integrates `1-u`, clock_start integrates `u*(1-u)` and
/// clock_end integrates `u*u`, all multiplying the rotating quadratures.
/// Keeping the complementary moments separately avoids subtracting nearly
/// equal constant/ramp gains after a long physical decay.
#[derive(Clone, Copy, Debug)]
pub(super) struct EndpointMoments {
    pub(super) start: [f64; 4],
    pub(super) clock_start: [f64; 4],
    pub(super) clock_end: [f64; 4],
}

impl EndpointMoments {
    pub(super) fn small(
        angle: f64,
        drag: f64,
        phase: f64,
        force: [f64; 4],
        delta: [f64; 4],
    ) -> Self {
        let mut quadratic = [0.0; 4];
        for column in 0..2 {
            // The original rotating linear-envelope action with one more
            // polynomial level: D'=phase*J*D+2*Q, Q'=phase*J*Q.
            let mut term = [0.0; 8];
            term[6 + column] = 1.0;
            for order in 1..=24 {
                let inverse = 1.0 / f64::from(order);
                term = [
                    term[1] * inverse,
                    (-angle * angle * term[0] - drag * term[1] + term[2]) * inverse,
                    (phase * term[3] + term[4]) * inverse,
                    (-phase * term[2] + term[5]) * inverse,
                    (phase * term[5] + 2.0 * term[6]) * inverse,
                    (-phase * term[4] + 2.0 * term[7]) * inverse,
                    phase * term[7] * inverse,
                    -phase * term[6] * inverse,
                ];
                quadratic[column] += term[0];
                quadratic[2 + column] += term[1];
            }
        }
        // These differences are well separated on the resolved small step;
        // subsequent doublings never form a complementary moment by subtraction.
        Self {
            start: std::array::from_fn(|i| force[i] - delta[i]),
            clock_start: std::array::from_fn(|i| delta[i] - quadratic[i]),
            clock_end: quadratic,
        }
    }

    pub(super) fn double(
        self,
        force: [f64; 4],
        delta: [f64; 4],
        flow: DampedTransition,
        phase: f64,
    ) -> Self {
        let (sine, cosine) = det::sin_cos(phase);
        let restoring = -(flow.angle * flow.b) * flow.angle;
        let combine = |first: [f64; 4], second: [f64; 4], scale: f64| {
            [
                scale
                    * (flow.position_retention * first[0] + flow.b * first[2] + cosine * second[0]
                        - sine * second[1]),
                scale
                    * (flow.position_retention * first[1]
                        + flow.b * first[3]
                        + sine * second[0]
                        + cosine * second[1]),
                2.0 * scale
                    * (restoring * first[0]
                        + flow.velocity_retention * first[2]
                        + cosine * second[2]
                        - sine * second[3]),
                2.0 * scale
                    * (restoring * first[1]
                        + flow.velocity_retention * first[3]
                        + sine * second[2]
                        + cosine * second[3]),
            ]
        };
        Self {
            // 1-u is (1+(1-v))/2 on the first half and (1-v)/2
            // on the second, whose carrier has advanced by phase.
            start: combine(
                std::array::from_fn(|i| force[i] + self.start[i]),
                self.start,
                0.125,
            ),
            // u*(1-u) is (v+v*(1-v))/4, then ((1-v)+v*(1-v))/4.
            clock_start: combine(
                std::array::from_fn(|i| delta[i] + self.clock_start[i]),
                std::array::from_fn(|i| self.start[i] + self.clock_start[i]),
                0.0625,
            ),
            // u^2 is v^2/4, then (1+2*v+v^2)/4.
            clock_end: combine(
                self.clock_end,
                std::array::from_fn(|i| force[i] + 2.0 * delta[i] + self.clock_end[i]),
                0.0625,
            ),
        }
    }
}
