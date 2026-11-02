//! Cross-frequency power forms, using the scalar entire flow and heat QR owner.

use super::super::integrals::{NODES, augmented_flow, heat_rotation};
use super::super::{DampedTransition, double_force, invalid, small_force};
use super::{MAX_COORDINATES, Sum};
use crate::time::{ModalAcousticState, ModalAcousticTimeError};
use fs_math::det;

#[derive(Clone, Debug)]
pub(super) struct TonePower {
    h: f64,
    omega: f64,
    scale: f64,
    heat_scale: f64,
    work_scale: f64,
    heat_root: Box<[f64]>,
    force_rows: Box<[f64]>,
    roundoff: f64,
}

// (resolved angle, drag, scaled time, original coordinate scale, doublings).
// x=[v,(sum Fc-omega^2*q)/kappa,Fc_0/kappa,Fs_0/kappa,...],
// kappa=max(omega,max Omega,1/h). The scaling includes tone count so BOTH
// generator norms remain resolved even though its residual row sums tones.
pub(super) fn scaling(
    omega: f64,
    zeta: f64,
    frequencies: &[f64],
    h: f64,
) -> Result<(f64, f64, f64, f64, usize), ModalAcousticTimeError> {
    let mut angle = omega * h;
    let mut drag = 2.0 * (zeta * angle);
    let mut scale = angle.max(1.0);
    for &frequency in frequencies {
        let phase = frequency * h;
        if !phase.is_finite() {
            return Err(invalid("tone phase must be representable"));
        }
        scale = scale.max(phase);
    }
    if !angle.is_finite() || !drag.is_finite() {
        return Err(invalid(
            "tone modal angle and damping must be representable",
        ));
    }
    let mut time = scale;
    let mut doublings = 0;
    let threshold = 0.125 / (frequencies.len() + 2) as f64;
    while drag.max(time) > threshold {
        angle *= 0.5;
        drag *= 0.5;
        time *= 0.5;
        doublings += 1;
    }
    Ok((angle, drag, time, scale, doublings))
}

impl TonePower {
    pub(super) fn try_new(
        omega: f64,
        zeta: f64,
        frequencies: &[f64],
        h: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        let (angle, drag, mut time, scale, doublings) = scaling(omega, zeta, frequencies, h)?;
        let ratio = (omega * h) / scale;
        let n = 2 + 2 * frequencies.len();
        let mut phases = frequencies
            .iter()
            .map(|frequency| frequency * h)
            .collect::<Vec<_>>();
        for phase in &mut phases {
            for _ in 0..doublings {
                *phase *= 0.5;
            }
        }
        let mut fraction: f64 = 1.0;
        for _ in 0..doublings {
            fraction *= 0.5;
        }
        let mut result = Self {
            h,
            omega,
            scale,
            heat_scale: det::sqrt(2.0 * (zeta * (omega * h))),
            work_scale: det::sqrt(scale),
            heat_root: vec![0.0; n * n].into_boxed_slice(),
            force_rows: vec![0.0; n * n].into_boxed_slice(),
            roundoff: (1024.0 + 256.0 * doublings as f64) * (n as f64 / 4.0) * f64::EPSILON,
        };
        let mut e = vec![0.0; n * n];
        let mut product = vec![0.0; n * n];
        let mut force_product = vec![0.0; n * n];
        let mut forces = vec![[0.0; 4]; frequencies.len()];
        let mut row = [0.0; MAX_COORDINATES];
        let mut phase_row = [0.0; MAX_COORDINATES];
        // The shared degree-15 Gaussian rule is used only after resolving the
        // full augmented generator. Long beating/oscillatory windows compose
        // these local integrals; they are never sampled with only eight nodes.
        for (node, weight) in NODES {
            let flow = DampedTransition::small(angle * node, drag * node);
            for ((force, &phase), pair) in forces
                .iter_mut()
                .zip(&phases)
                .zip(phase_row[2..n].chunks_exact_mut(2))
            {
                *force = small_force(angle * node, drag * node, phase * node, false);
                let (s, c) = det::sin_cos(phase * node);
                pair[0] = c;
                pair[1] = s;
            }
            fill_flow(&mut e, flow, &forces, &phases, node, time * node, ratio);
            let factor = det::sqrt(fraction) * det::sqrt(weight);
            for j in 0..n {
                row[j] = factor * e[j];
                phase_row[j] *= factor;
            }
            result.append_pair(&mut row[..n], &mut phase_row[..n]);
        }
        let mut flow = DampedTransition::small(angle, drag);
        for (force, &phase) in forces.iter_mut().zip(&phases) {
            *force = small_force(angle, drag, phase, false);
        }
        for _ in 0..doublings {
            fill_flow(&mut e, flow, &forces, &phases, 1.0, time, ratio);
            // E has two dense rows plus independent 2x2 source rotations.
            // Exploit that structure for O(n^2) right multiplication rather
            // than treating the augmented flow as a dense generic matrix.
            right_multiply(&result.heat_root, &e, &mut product, n);
            right_multiply(&result.force_rows, &e, &mut force_product, n);
            for (velocity, force) in product
                .chunks_exact_mut(n)
                .zip(force_product.chunks_exact_mut(n))
            {
                result.append_pair(velocity, force);
            }
            for (force, phase) in forces.iter_mut().zip(&mut phases) {
                *force = double_force(*force, flow, *phase, false);
                *phase *= 2.0;
            }
            flow.double();
            time *= 2.0;
        }
        if result
            .heat_root
            .iter()
            .chain(&result.force_rows)
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "simultaneous-tone power coefficients must be representable",
            ));
        }
        Ok(result)
    }

    // Apply the SAME orthogonal rotations to sampled velocity and force rows.
    // R^T R is positive velocity power; R^T T is the independent F*v form.
    // Keep the factors, not the expanded signed quadratic: close opposing
    // tones cancel before products, preserving their small physical beat work.
    // A discarded zero-velocity residual has zero cross power and stays zero
    // under right composition with E, so no force-only factor is needed.
    fn append_pair(&mut self, velocity: &mut [f64], force: &mut [f64]) {
        let n = velocity.len();
        for i in 0..n {
            let Some((norm, c, s)) = heat_rotation(self.heat_root[i * n + i], velocity[i]) else {
                continue;
            };
            self.heat_root[i * n + i] = norm;
            for j in i + 1..n {
                let old = self.heat_root[i * n + j];
                self.heat_root[i * n + j] = c * old + s * velocity[j];
                velocity[j] = -s * old + c * velocity[j];
            }
            for j in 0..n {
                let old = self.force_rows[i * n + j];
                self.force_rows[i * n + j] = c * old + s * force[j];
                force[j] = -s * old + c * force[j];
            }
        }
    }

    // Complete signed work and positive heat, never a sum of per-tone heats.
    pub(super) fn evaluate(
        &self,
        state: ModalAcousticState,
        cosine: &[f64],
        sine: &[f64],
    ) -> Result<(f64, f64, f64, f64), ModalAcousticTimeError> {
        let n = 2 + 2 * cosine.len();
        let mut x = [0.0; MAX_COORDINATES];
        let mut total = Sum::default();
        x[0] = state.velocity_m_sqrt_kg_per_s;
        for (tone, (&c, &s)) in cosine.iter().zip(sine).enumerate() {
            let c = self.h * c;
            let s = self.h * s;
            total.add(c);
            x[2 + 2 * tone] = c / self.scale;
            x[3 + 2 * tone] = s / self.scale;
        }
        let restoring = (self.omega * state.displacement_m_sqrt_kg) * (self.omega * self.h);
        x[1] = (total.value() - restoring) / self.scale;
        let x = &x[..n];
        if x.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "simultaneous-tone power coordinates must be finite",
            ));
        }
        let mut heat = Sum::default();
        let mut heat_scale = 0.0;
        let mut work = Sum::default();
        let mut work_scale = 0.0;
        for (velocity, force) in self
            .heat_root
            .chunks_exact(n)
            .zip(self.force_rows.chunks_exact(n))
        {
            let mut v = Sum::default();
            let mut f = Sum::default();
            let (mut v_magnitude, mut f_magnitude) = (0.0, 0.0);
            for ((&vr, &fr), &coordinate) in velocity.iter().zip(force).zip(x) {
                v.add(vr * coordinate);
                f.add(fr * coordinate);
                v_magnitude += (vr * coordinate).abs();
                f_magnitude += (fr * coordinate).abs();
            }
            let viscous = self.heat_scale * v.value();
            heat.add(viscous * viscous);
            let heat_magnitude = self.heat_scale * v_magnitude;
            heat_scale += heat_magnitude * heat_magnitude;
            work.add((self.work_scale * v.value()) * (self.work_scale * f.value()));
            work_scale += (self.work_scale * v_magnitude) * (self.work_scale * f_magnitude);
        }
        let values = (
            work.value(),
            self.roundoff * work_scale,
            heat.value(),
            self.roundoff * heat_scale,
        );
        if ![values.0, values.1, values.2, values.3]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid(
                "simultaneous-tone work and heat must be representable",
            ));
        }
        Ok(values)
    }
}

fn fill_flow(
    e: &mut [f64],
    flow: DampedTransition,
    forces: &[[f64; 4]],
    phases: &[f64],
    phase_fraction: f64,
    time: f64,
    ratio: f64,
) {
    let n = 2 + 2 * phases.len();
    e.fill(0.0);
    for (tone, (&force, &phase)) in forces.iter().zip(phases).enumerate() {
        let block = augmented_flow(flow, force, phase * phase_fraction, time, ratio);
        let c = 2 + 2 * tone;
        for i in 0..2 {
            e[i * n] = block[i][0];
            e[i * n + 1] = block[i][1];
            e[i * n + c] = block[i][2];
            e[i * n + c + 1] = block[i][3];
        }
        e[c * n + c] = block[2][2];
        e[c * n + c + 1] = block[2][3];
        e[(c + 1) * n + c] = block[3][2];
        e[(c + 1) * n + c + 1] = block[3][3];
    }
}

fn right_multiply(a: &[f64], e: &[f64], product: &mut [f64], n: usize) {
    for i in 0..n {
        for j in 0..n {
            let mut value = a[i * n] * e[j] + a[i * n + 1] * e[n + j];
            if j >= 2 {
                let pair = 2 + 2 * ((j - 2) / 2);
                value +=
                    a[i * n + pair] * e[pair * n + j] + a[i * n + pair + 1] * e[(pair + 1) * n + j];
            }
            product[i * n + j] = value;
        }
    }
}
