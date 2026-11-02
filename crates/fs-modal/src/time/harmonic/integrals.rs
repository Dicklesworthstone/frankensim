//! Independent finite-interval force work and positive viscous heat.

use super::super::{ModalAcousticState, ModalAcousticTimeError};
use super::{DampedTransition, HarmonicForceTransition, double_force, invalid, small_force};
use fs_math::det;

type Matrix = [[f64; 4]; 4];

#[derive(Clone, Debug)]
pub(super) struct HarmonicIntegrals {
    duration: f64,
    omega: f64,
    scale: f64,
    heat_root: Matrix,
    work: Matrix,
    roundoff: f64,
}

impl HarmonicIntegrals {
    pub fn try_new(
        step: &HarmonicForceTransition,
        forcing_frequency: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        let mut angle = step.transition.angle;
        let mut drag = step.transition.drag;
        let mut phase = forcing_frequency * step.duration_s;
        let scale = angle.max(phase).max(1.0);
        let frequency_ratio = angle / scale;
        let mut scaled_time = scale;
        let mut doublings = 0;
        // Coordinates are [v, (F-omega^2*q)/kappa, F/kappa, Fs/kappa],
        // kappa=max(omega,Omega,1/h). All four have velocity units. On this
        // short interval the dimensionless generator's 1/inf norms <= 1/4.
        while drag.max(scaled_time).max(phase) > 0.125 {
            angle *= 0.5;
            drag *= 0.5;
            phase *= 0.5;
            scaled_time *= 0.5;
            doublings += 1;
        }
        let mut result = Self {
            duration: step.duration_s,
            omega: step.angular_frequency_rad_s,
            scale,
            heat_root: [[0.0; 4]; 4],
            work: [[0.0; 4]; 4],
            roundoff: (1024.0 + 256.0 * f64::from(doublings)) * f64::EPSILON,
        };
        // Gauss-Legendre on the resolved short interval only. Its degree-15
        // rule integrates products of these analytic trajectories to below
        // f64 rounding at the declared generator norm. This is not interval
        // certification: the returned allowances cover floating arithmetic.
        const NODES: [(f64, f64); 8] = [
            (0.019855071751231884, 0.05061426814518813),
            (0.10166676129318664, 0.11119051722668724),
            (0.2372337950418355, 0.15685332293894366),
            (0.4082826787521751, 0.181341891689181),
            (0.5917173212478248, 0.181341891689181),
            (0.7627662049581645, 0.15685332293894366),
            (0.8983332387068134, 0.11119051722668724),
            (0.9801449282487681, 0.05061426814518813),
        ];
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
            if drag != 0.0 {
                let factor = det::sqrt(drag * weight);
                result.append_heat_row(e[0].map(|x| factor * x));
            }
            for i in 0..4 {
                for j in i..4 {
                    let entry =
                        0.5 * scaled_time * weight * (e[2][i] * e[0][j] + e[0][i] * e[2][j]);
                    result.work[i][j] += entry;
                    result.work[j][i] = result.work[i][j];
                }
            }
        }
        let mut flow = DampedTransition::small(angle, drag);
        let mut force = small_force(angle, drag, phase, false);
        for _ in 0..doublings {
            let e = augmented_flow(flow, force, phase, scaled_time, frequency_ratio);
            // Integral(2h) = Integral(h) + E(h)^T Integral(h) E(h).
            // Heat is QR of [R;R E], so positivity needs no clipping or
            // subtraction of mechanical energies. Work is a separate signed
            // force/velocity integral, including source quadrature cross terms.
            let rows = multiply(result.heat_root, e);
            for row in rows {
                result.append_heat_row(row);
            }
            let right = multiply(result.work, e);
            for i in 0..4 {
                for j in i..4 {
                    let mut entry = 0.0;
                    for k in 0..4 {
                        entry += e[k][i] * right[k][j];
                    }
                    result.work[i][j] += entry;
                    result.work[j][i] = result.work[i][j];
                }
            }
            force = double_force(force, flow, phase, false);
            flow.double();
            phase *= 2.0;
            scaled_time *= 2.0;
        }
        if result
            .heat_root
            .iter()
            .chain(&result.work)
            .flatten()
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "harmonic power integral coefficients must be representable",
            ));
        }
        Ok(result)
    }

    // (work, work allowance, heat, heat allowance), all in joules.
    pub fn evaluate(
        &self,
        state: ModalAcousticState,
        cosine_force: f64,
        sine_force: f64,
    ) -> Result<(f64, f64, f64, f64), ModalAcousticTimeError> {
        let fc = self.duration * cosine_force;
        let fs = self.duration * sine_force;
        let restoring = (self.omega * state.displacement_m_sqrt_kg) * (self.omega * self.duration);
        let x = [
            state.velocity_m_sqrt_kg_per_s,
            (fc - restoring) / self.scale,
            fc / self.scale,
            fs / self.scale,
        ];
        if !x.into_iter().all(f64::is_finite) {
            return Err(invalid("harmonic power coordinates must be finite"));
        }
        let mut heat = 0.0;
        let mut heat_scale = 0.0;
        for row in &self.heat_root {
            let mut value = 0.0;
            let mut magnitude = 0.0;
            for i in 0..4 {
                let term = row[i] * x[i];
                value += term;
                magnitude += term.abs();
            }
            heat += value * value;
            heat_scale += magnitude * magnitude;
        }
        let mut work = 0.0;
        let mut work_scale = 0.0;
        for i in 0..4 {
            for j in i..4 {
                let term = self.work[i][j] * x[i] * x[j] * if i == j { 1.0 } else { 2.0 };
                work += term;
                work_scale += term.abs();
            }
        }
        let work_allowance = self.roundoff * work_scale;
        let heat_allowance = self.roundoff * heat_scale;
        if ![work, work_allowance, heat, heat_allowance]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid("harmonic work and heat must be representable"));
        }
        Ok((work, work_allowance, heat, heat_allowance))
    }

    fn append_heat_row(&mut self, mut row: [f64; 4]) {
        for i in 0..4 {
            let old = self.heat_root[i][i];
            let scale = old.abs().max(row[i].abs());
            if scale == 0.0 {
                continue;
            }
            let norm = scale * det::sqrt((old / scale).powi(2) + (row[i] / scale).powi(2));
            let (c, s) = (old / norm, row[i] / norm);
            self.heat_root[i][i] = norm;
            for (j, value) in row.iter_mut().enumerate().skip(i + 1) {
                let old = self.heat_root[i][j];
                self.heat_root[i][j] = c * old + s * *value;
                *value = -s * old + c * *value;
            }
        }
    }
}

fn augmented_flow(
    flow: DampedTransition,
    force: [f64; 4],
    phase: f64,
    scaled_time: f64,
    frequency_ratio: f64,
) -> Matrix {
    let [qc, qs, vc, vs] = force;
    let (s, c) = det::sin_cos(phase);
    [
        [
            flow.velocity_retention,
            scaled_time * flow.b,
            -scaled_time * (phase * qs),
            scaled_time * (phase * qc),
        ],
        [
            -frequency_ratio * (flow.angle * flow.b),
            flow.position_retention,
            -phase * (vs + flow.drag * qs),
            phase * (vc + flow.drag * qc),
        ],
        [0.0, 0.0, c, s],
        [0.0, 0.0, -s, c],
    ]
}

fn multiply(a: Matrix, b: Matrix) -> Matrix {
    let mut product = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                product[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    product
}
