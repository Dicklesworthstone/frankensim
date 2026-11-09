//! Tangents of the original rotating affine-envelope convolution.

use super::{
    DampedTransition, HarmonicForceTransition, ModalAcousticTimeError, MomentDerivative, Parameter,
    PreparedHarmonicForceSensitivity, double_derivative, double_force, invalid, small_derivative,
    small_force, small_ramp_moments,
};
use crate::time::harmonic::affine::{EndpointMoments, double_delta, small_delta};
use fs_math::det;

/// Extra cold coefficients only for affine banks. The original scalar cache
/// still owns the homogeneous map, constant drive and near-DC imbalance.
#[derive(Clone, Copy, Debug)]
pub(in crate::time::harmonic) struct AffineForceSensitivity {
    stiffness: [f64; 4],
    damping: [f64; 4],
    moments: EndpointMoments,
}

impl AffineForceSensitivity {
    pub(in crate::time::harmonic) fn prepare(
        motion: &HarmonicForceTransition,
    ) -> Result<Self, ModalAcousticTimeError> {
        let (mut angle, mut drag, mut phase) = (
            motion.transition.angle,
            motion.transition.drag,
            motion.forcing_phase_rad,
        );
        let mut doublings = 0;
        while angle.max(drag).max(phase) > 0.25 {
            angle *= 0.5;
            drag *= 0.5;
            phase *= 0.5;
            doublings += 1;
        }
        let mut flow = DampedTransition::small(angle, drag);
        let mut ramp = small_ramp_moments(angle, drag).0;
        let mut stiffness = MomentDerivative::small(flow, Parameter::Stiffness);
        let mut damping = MomentDerivative::small(flow, Parameter::Damping);
        // The affine doubling needs the complete constant-force columns, even
        // when the original nominal motion uses a separate cos-1 remainder.
        let mut force = small_force(angle, drag, phase, false);
        let mut force_k = small_derivative(angle, drag, phase, false, Parameter::Stiffness);
        let mut force_c = small_derivative(angle, drag, phase, false, Parameter::Damping);
        let mut delta = small_delta(angle, drag, phase);
        let mut delta_k = small_delta_derivative(angle, drag, phase, Parameter::Stiffness);
        let mut delta_c = small_delta_derivative(angle, drag, phase, Parameter::Damping);
        let mut moments = EndpointMoments::small(angle, drag, phase, force, delta);
        for _ in 0..doublings {
            delta_k = double_delta_derivative(
                delta_k,
                delta,
                force_k,
                flow,
                stiffness,
                phase,
                Parameter::Stiffness,
            );
            delta_c = double_delta_derivative(
                delta_c,
                delta,
                force_c,
                flow,
                damping,
                phase,
                Parameter::Damping,
            );
            force_k = double_derivative(
                force_k,
                force,
                flow,
                stiffness,
                phase,
                false,
                Parameter::Stiffness,
            );
            force_c = double_derivative(
                force_c,
                force,
                flow,
                damping,
                phase,
                false,
                Parameter::Damping,
            );
            moments = moments.double(force, delta, flow, phase);
            delta = double_delta(delta, force, flow, phase);
            force = double_force(force, flow, phase, false);
            let before = flow;
            flow.double();
            stiffness = stiffness.double(before, flow, ramp, Parameter::Stiffness);
            damping = damping.double(before, flow, ramp, Parameter::Damping);
            ramp = ((1.0 + before.position_retention) * ramp + (1.0 + before.b) * before.g) * 0.125;
            phase *= 2.0;
        }
        if delta_k
            .into_iter()
            .chain(delta_c)
            .chain(moments.start)
            .chain(moments.clock_start)
            .chain(moments.clock_end)
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "affine harmonic sensitivities must be representable",
            ));
        }
        Ok(Self {
            stiffness: delta_k,
            damping: delta_c,
            moments,
        })
    }

    // Add this delta-only contribution to the original constant source before
    // the shared mechanical contraction. Do not advance another initial state.
    pub(in crate::time::harmonic) fn parameter_forcing(
        &self,
        duration: f64,
        delta_cosine: f64,
        delta_sine: f64,
    ) -> [f64; 5] {
        let (c, s) = (duration * delta_cosine, duration * delta_sine);
        let [kqc, kqs, kvc, kvs] = self.stiffness;
        let [cqc, cqs, cvc, cvs] = self.damping;
        [
            0.0,
            kqc * c + kqs * s,
            kvc * c + kvs * s,
            cqc * c + cqs * s,
            cvc * c + cvs * s,
        ]
    }

    pub(in crate::time::harmonic) fn clock_partials(
        &self,
        scalar: &PreparedHarmonicForceSensitivity,
        cosine_start: f64,
        sine_start: f64,
        cosine_end: f64,
        sine_end: f64,
    ) -> Result<[f64; 4], ModalAcousticTimeError> {
        let h = scalar.motion.duration_s;
        let (c0, s0, c1, s1) = (
            h * cosine_start,
            h * sine_start,
            h * cosine_end,
            h * sine_end,
        );
        if ![c0, s0, c1, s1].into_iter().all(f64::is_finite) {
            return Err(invalid("affine source impulses must be representable"));
        }
        let contract = |start: [f64; 4], end: [f64; 4], row: usize| {
            start[row].mul_add(s0, -start[row + 1] * c0) + end[row].mul_add(s1, -end[row + 1] * c1)
        };
        // dF/dOmega = h*u*((1-u)*[S0,-C0]+u*[S1,-C1]) in the
        // rotating carrier basis. No individual endpoint difference is needed.
        let result = [
            h * (h * contract(self.moments.clock_start, self.moments.clock_end, 0)),
            h * contract(self.moments.clock_start, self.moments.clock_end, 2),
            h * contract(self.moments.start, scalar.force_envelope, 0),
            contract(self.moments.start, scalar.force_envelope, 2),
        ];
        if !result.into_iter().all(f64::is_finite) {
            return Err(invalid(
                "affine source-clock partials must be representable",
            ));
        }
        Ok(result)
    }

    pub(in crate::time::harmonic) fn input_cotangents(
        &self,
        scalar: &PreparedHarmonicForceSensitivity,
        q_bar: f64,
        v_bar: f64,
    ) -> Result<[f64; 4], ModalAcousticTimeError> {
        let h = scalar.motion.duration_s;
        let qb = h * q_bar;
        let contract =
            |gain: [f64; 4], column: usize| h * gain[column].mul_add(qb, gain[2 + column] * v_bar);
        let result = [
            contract(self.moments.start, 0),
            contract(self.moments.start, 1),
            contract(scalar.force_envelope, 0),
            contract(scalar.force_envelope, 1),
        ];
        if !result.into_iter().all(f64::is_finite) {
            return Err(invalid(
                "affine endpoint-force cotangents must be representable",
            ));
        }
        Ok(result)
    }
}

fn small_delta_derivative(angle: f64, drag: f64, phase: f64, parameter: Parameter) -> [f64; 4] {
    let (da, dd, _) = parameter.seeds();
    let mut gains = [0.0; 4];
    for column in 0..2 {
        let mut term = [0.0; 6];
        term[4 + column] = 1.0;
        let mut derivative = [0.0; 2];
        for order in 1..=24 {
            let inverse = 1.0 / f64::from(order);
            derivative = [
                derivative[1] * inverse,
                (-da * term[0]
                    - angle * angle * derivative[0]
                    - dd * term[1]
                    - drag * derivative[1])
                    * inverse,
            ];
            term = [
                term[1] * inverse,
                (-angle * angle * term[0] - drag * term[1] + term[2]) * inverse,
                (phase * term[3] + term[4]) * inverse,
                (-phase * term[2] + term[5]) * inverse,
                phase * term[5] * inverse,
                -phase * term[4] * inverse,
            ];
            gains[column] += derivative[0];
            gains[2 + column] += derivative[1];
        }
    }
    gains
}

fn double_delta_derivative(
    derivative: [f64; 4],
    delta: [f64; 4],
    force_derivative: [f64; 4],
    flow: DampedTransition,
    tangent: MomentDerivative,
    phase: f64,
    parameter: Parameter,
) -> [f64; 4] {
    let (da, _, scale) = parameter.seeds();
    let (sine, cosine) = det::sin_cos(phase);
    let [qc, qs, vc, vs] = delta;
    let [dqc, dqs, dvc, dvs] = derivative;
    let [fqc, fqs, fvc, fvs] = force_derivative;
    let restoring = -(flow.angle * flow.b) * flow.angle;
    let dr = -da * flow.b - flow.angle * (flow.angle * tangent.b);
    [
        0.125
            * scale
            * (tangent.position * qc
                + flow.position_retention * dqc
                + tangent.b * vc
                + flow.b * dvc
                + (dqc + fqc) * cosine
                - (dqs + fqs) * sine),
        0.125
            * scale
            * (tangent.position * qs
                + flow.position_retention * dqs
                + tangent.b * vs
                + flow.b * dvs
                + (dqc + fqc) * sine
                + (dqs + fqs) * cosine),
        0.25 * scale
            * (dr * qc
                + restoring * dqc
                + tangent.velocity * vc
                + flow.velocity_retention * dvc
                + (dvc + fvc) * cosine
                - (dvs + fvs) * sine),
        0.25 * scale
            * (dr * qs
                + restoring * dqs
                + tangent.velocity * vs
                + flow.velocity_retention * dvs
                + (dvc + fvc) * sine
                + (dvs + fvs) * cosine),
    ]
}
