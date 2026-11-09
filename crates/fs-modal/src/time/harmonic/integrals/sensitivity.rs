//! Tangents of the original resolved work/heat integrals and their doubling.
//!
//! The conditioning coordinate kappa stays fixed during differentiation. It
//! cancels from the physical integral; differentiating its max() would add an
//! artificial branch to a smooth physical observable.

use super::{
    DampedTransition, HarmonicForceTransition, HarmonicIntegrals, Matrix, ModalAcousticState,
    ModalAcousticTimeError, NODES, augmented_flow, double_force, invalid, multiply, small_force,
};
use crate::time::harmonic::HarmonicForcePullback;

/// Derivatives of each independent physical power integral with respect to the
/// original initial state, force quadratures, k, c, drive frequency and phase.
/// Units are joules per unit of the corresponding input. A zero damping rate
/// still has a nonzero right derivative of heat when the body is moving.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicPowerGradient {
    pub input_work: HarmonicForcePullback,
    pub viscous_dissipation: HarmonicForcePullback,
}

// Three derivatives are with respect to full-interval k*h^2, c*h and Omega*h.
// Both observable forms use the primal's balanced four-dimensional chart.
pub(in crate::time::harmonic) struct HarmonicIntegralSensitivity {
    work: [Matrix; 3],
    heat: [Matrix; 3],
}

impl HarmonicIntegralSensitivity {
    pub(in crate::time::harmonic) fn prepare(
        motion: &HarmonicForceTransition,
        integrals: &HarmonicIntegrals,
    ) -> Result<Self, ModalAcousticTimeError> {
        let mut angle = motion.transition.angle;
        let mut drag = motion.transition.drag;
        let mut phase = motion.forcing_phase_rad;
        let mut scaled_time = integrals.scale;
        let mut fraction = 1.0;
        let mut doublings = 0;
        while drag.max(scaled_time).max(phase) > 0.125 {
            angle *= 0.5;
            drag *= 0.5;
            phase *= 0.5;
            scaled_time *= 0.5;
            fraction *= 0.5;
            doublings += 1;
        }
        let frequency_ratio = motion.transition.angle / integrals.scale;
        let mut work = [[0.0; 4]; 4];
        let mut heat = [[0.0; 4]; 4];
        let mut result = Self {
            work: [[[0.0; 4]; 4]; 3],
            heat: [[[0.0; 4]; 4]; 3],
        };
        // Use the primal's short Gauss rule and continuous transition. Only its
        // analytic tangent is new; no perturbed solve or time sampling is used.
        for (node, weight) in NODES {
            let flow = DampedTransition::small(angle * node, drag * node);
            let force = small_force(angle * node, drag * node, phase * node, false);
            let e = augmented_flow(
                flow,
                force,
                phase * node,
                scaled_time * node,
                frequency_ratio,
            );
            let de = small_action_derivatives(
                angle * node,
                drag * node,
                phase * node,
                scaled_time * node,
                fraction * node,
                integrals.scale,
            );
            for i in 0..4 {
                for j in i..4 {
                    work[i][j] +=
                        0.5 * scaled_time * weight * (e[2][i] * e[0][j] + e[0][i] * e[2][j]);
                    heat[i][j] += drag * weight * e[0][i] * e[0][j];
                    work[j][i] = work[i][j];
                    heat[j][i] = heat[i][j];
                    for p in 0..3 {
                        result.work[p][i][j] += 0.5
                            * scaled_time
                            * weight
                            * (de[p][2][i] * e[0][j]
                                + e[2][i] * de[p][0][j]
                                + de[p][0][i] * e[2][j]
                                + e[0][i] * de[p][2][j]);
                        result.heat[p][i][j] +=
                            drag * weight * (de[p][0][i] * e[0][j] + e[0][i] * de[p][0][j]);
                        if p == 1 {
                            // Differentiate c*v^2 itself, including c=0.
                            result.heat[p][i][j] += fraction * weight * e[0][i] * e[0][j];
                        }
                        result.work[p][j][i] = result.work[p][i][j];
                        result.heat[p][j][i] = result.heat[p][i][j];
                    }
                }
            }
        }
        let mut flow = DampedTransition::small(angle, drag);
        let mut force = small_force(angle, drag, phase, false);
        let mut de =
            small_action_derivatives(angle, drag, phase, scaled_time, fraction, integrals.scale);
        for _ in 0..doublings {
            let e = augmented_flow(flow, force, phase, scaled_time, frequency_ratio);
            for p in 0..3 {
                // D(Q + E^T Q E), with all three chain-rule contributions.
                result.work[p] = double_form_direction(work, result.work[p], e, de[p]);
                result.heat[p] = double_form_direction(heat, result.heat[p], e, de[p]);
                let left = multiply(de[p], e);
                let right = multiply(e, de[p]);
                for i in 0..4 {
                    for j in 0..4 {
                        de[p][i][j] = left[i][j] + right[i][j];
                    }
                }
            }
            work = double_form(work, e);
            heat = double_form(heat, e);
            force = double_force(force, flow, phase, false);
            flow.double();
            phase *= 2.0;
            scaled_time *= 2.0;
        }
        if result
            .work
            .iter()
            .chain(&result.heat)
            .flatten()
            .flatten()
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "harmonic power sensitivities must be representable",
            ));
        }
        Ok(result)
    }

    pub(in crate::time::harmonic) fn gradients(
        &self,
        base: &HarmonicIntegrals,
        state: ModalAcousticState,
        cosine: f64,
        sine: f64,
    ) -> Result<HarmonicPowerGradient, ModalAcousticTimeError> {
        let x = base.coordinates(state, cosine, sine)?;
        let mut work_gradient = [0.0; 4];
        let mut heat_gradient = [0.0; 4];
        let heat_action = base.heat_root.map(|row| dot(row, x));
        for i in 0..4 {
            work_gradient[i] = 2.0 * dot(base.work[i], x);
            heat_gradient[i] =
                2.0 * dot(core::array::from_fn(|j| base.heat_root[j][i]), heat_action);
        }
        Ok(HarmonicPowerGradient {
            input_work: physical_gradient(base, state, cosine, sine, x, work_gradient, &self.work)?,
            viscous_dissipation: physical_gradient(
                base,
                state,
                cosine,
                sine,
                x,
                heat_gradient,
                &self.heat,
            )?,
        })
    }
}

fn physical_gradient(
    base: &HarmonicIntegrals,
    state: ModalAcousticState,
    cosine: f64,
    sine: f64,
    x: [f64; 4],
    g: [f64; 4],
    directions: &[Matrix; 3],
) -> Result<HarmonicForcePullback, ModalAcousticTimeError> {
    let h = base.duration;
    let h_over_scale = h / base.scale;
    let cosine_force = (g[1] + g[2]) * h_over_scale;
    let sine_force = g[3] * h_over_scale;
    let result = HarmonicForcePullback {
        initial_displacement: -((g[1] * base.omega) * (base.omega * h)) / base.scale,
        initial_velocity: g[0],
        cosine_force,
        sine_force,
        stiffness: (quadratic(directions[0], x) * h) * h
            - (g[1] * h_over_scale) * state.displacement_m_sqrt_kg,
        damping_rate: quadratic(directions[1], x) * h,
        forcing_angular_frequency: quadratic(directions[2], x) * h,
        initial_phase: sine * cosine_force - cosine * sine_force,
    };
    if ![
        result.initial_displacement,
        result.initial_velocity,
        result.cosine_force,
        result.sine_force,
        result.stiffness,
        result.damping_rate,
        result.forcing_angular_frequency,
        result.initial_phase,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(invalid(
            "physical harmonic power derivatives must be finite",
        ));
    }
    Ok(result)
}

// Frechet action of the same short augmented flow. The original entire
// exponential's differentiated recurrence has no pole/detuning denominators.
fn small_action_derivatives(
    angle: f64,
    drag: f64,
    phase: f64,
    scaled_time: f64,
    fraction: f64,
    scale: f64,
) -> [Matrix; 3] {
    let a = [
        [-drag, scaled_time, 0.0, 0.0],
        [-(angle / scaled_time) * angle, 0.0, 0.0, phase],
        [0.0, 0.0, 0.0, phase],
        [0.0, 0.0, -phase, 0.0],
    ];
    let mut da = [[[0.0; 4]; 4]; 3];
    da[0][1][0] = -fraction / scale;
    da[1][0][0] = -fraction;
    da[2][1][3] = fraction;
    da[2][2][3] = fraction;
    da[2][3][2] = -fraction;
    let mut term = [[0.0; 4]; 4];
    for (i, row) in term.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let mut dt = [[[0.0; 4]; 4]; 3];
    let mut result = dt;
    for order in 1..=24 {
        let inverse = 1.0 / f64::from(order);
        for p in 0..3 {
            let old = multiply(a, dt[p]);
            let source = multiply(da[p], term);
            for i in 0..4 {
                for j in 0..4 {
                    dt[p][i][j] = (old[i][j] + source[i][j]) * inverse;
                    result[p][i][j] += dt[p][i][j];
                }
            }
        }
        term = multiply(a, term);
        for row in &mut term {
            for v in row {
                *v *= inverse;
            }
        }
    }
    result
}

fn double_form(q: Matrix, e: Matrix) -> Matrix {
    let right = multiply(q, e);
    let mut result = q;
    for i in 0..4 {
        for j in i..4 {
            result[i][j] += dot(
                core::array::from_fn(|k| e[k][i]),
                core::array::from_fn(|k| right[k][j]),
            );
            result[j][i] = result[i][j];
        }
    }
    result
}

fn double_form_direction(q: Matrix, dq: Matrix, e: Matrix, de: Matrix) -> Matrix {
    let right = multiply(q, e);
    let dq_e = multiply(dq, e);
    let q_de = multiply(q, de);
    let mut result = dq;
    for i in 0..4 {
        for j in i..4 {
            result[i][j] += dot(
                core::array::from_fn(|k| de[k][i]),
                core::array::from_fn(|k| right[k][j]),
            );
            result[i][j] += dot(
                core::array::from_fn(|k| e[k][i]),
                core::array::from_fn(|k| dq_e[k][j] + q_de[k][j]),
            );
            result[j][i] = result[i][j];
        }
    }
    result
}

fn quadratic(q: Matrix, x: [f64; 4]) -> f64 {
    dot(x, q.map(|row| dot(row, x)))
}

fn dot(a: [f64; 4], b: [f64; 4]) -> f64 {
    let (mut sum, mut correction) = (0.0, 0.0);
    for i in 0..4 {
        let value = a[i] * b[i] - correction;
        let next = sum + value;
        correction = (next - sum) - value;
        sum = next;
    }
    sum
}
