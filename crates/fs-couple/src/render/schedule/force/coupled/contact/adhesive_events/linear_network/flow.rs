//! Matrix-free entire flow and independent physical resistance integrals.
use super::{
    Error, LinearNetworkEventFlow, LinearNetworkEventStep, LinearNetworkPhase, ORDER, apply, dot,
    finite, poll,
};
use fs_exec::CancelGate;
use fs_math::det;

impl LinearNetworkEventFlow {
    pub(super) fn prepare_cell(
        &mut self,
        external: &[f64],
        phase: LinearNetworkPhase,
        kinetic: f64,
        remaining: f64,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        let n = self.q.len();
        let dynamic_scale = det::sqrt(self.stiffness_norm).max(self.damping_norm);
        let scale = if dynamic_scale == 0.0 {
            1.0 / remaining
        } else {
            dynamic_scale
        };
        let projection = if phase == LinearNetworkPhase::Sticking {
            self.projection_norm
        } else {
            1.0
        };
        let rate =
            finite(scale.max(projection * (self.stiffness_norm / scale + self.damping_norm)))?;
        let h = if dynamic_scale == 0.0 {
            remaining
        } else {
            remaining.min(0.25 / rate)
        };
        if !h.is_finite() || h <= 0.0 {
            return Err(Error::Budget("resolved trajectory cell"));
        }
        self.cell_scale = scale;
        self.cell_rate = rate;
        self.cell_h = h;
        self.q_coeff[..n].copy_from_slice(&self.q);
        self.v_coeff[..n].copy_from_slice(&self.v);
        self.velocity_coeff.fill(0.0);
        self.reaction_coeff.fill(0.0);
        let mut forcing_norm: f64 = 0.0;
        for j in 0..=ORDER {
            poll(gate)?;
            let q = &self.q_coeff[j * n..(j + 1) * n];
            let v = &self.v_coeff[j * n..(j + 1) * n];
            apply(
                &self.stiffness,
                &self.damping,
                &self.links,
                q,
                v,
                if j == 0 { Some(external) } else { None },
                &mut self.accel,
            )?;
            self.operator_evaluations += 1;
            if self.operator_evaluations > self.config.maximum_operator_evaluations {
                return Err(Error::Budget("operator evaluations"));
            }
            self.velocity_coeff[j] = dot(&self.tangent, v)?;
            self.reaction_coeff[j] = finite(dot(&self.tangent, &self.accel)? / self.norm_b2)?;
            let reaction = match phase {
                LinearNetworkPhase::Sticking => self.reaction_coeff[j],
                LinearNetworkPhase::Sliding { direction } if j == 0 => direction as f64 * kinetic,
                _ => 0.0,
            };
            for i in 0..n {
                self.accel[i] = finite(self.accel[i] - self.tangent[i] * reaction)?;
            }
            if j == 0 {
                forcing_norm = self.accel.iter().map(|a| a.abs()).fold(0.0, f64::max);
            }
            if j < ORDER {
                let factor = h / (j + 1) as f64;
                for i in 0..n {
                    self.q_coeff[(j + 1) * n + i] = finite(factor * self.v_coeff[j * n + i])?;
                    self.v_coeff[(j + 1) * n + i] = finite(factor * self.accel[i])?;
                }
            }
        }
        // An exactly stationary constrained state has constant reaction and no
        // unresolved clock scale. Keep its complete held interval in one cell
        // even when the original network has very high elastic frequencies.
        let stationary = self.q_coeff[n..]
            .iter()
            .chain(&self.v_coeff[n..])
            .all(|x| *x == 0.0);
        if stationary {
            self.cell_h = remaining;
        }
        let h = self.cell_h;
        let norm = self
            .q
            .iter()
            .map(|q| scale * q.abs())
            .chain(self.v.iter().map(|v| v.abs()))
            .fold(0.0, f64::max);
        let mut tail = 1.0;
        if !stationary {
            for j in 1..=ORDER {
                tail *= h * rate / j as f64;
            }
        }
        // A nilpotent free-body operator terminates after q's quadratic term.
        if dynamic_scale == 0.0 || stationary {
            tail = 0.0;
        }
        let arithmetic = 32.0 * f64::EPSILON * (n + self.links.len() + ORDER) as f64;
        self.cell_state_error = finite((3.0 * tail + arithmetic) * (norm + h * forcing_norm))?;
        Ok(())
    }

    fn evaluate(&mut self, z: f64) -> Result<(), Error> {
        let n = self.q.len();
        for i in 0..n {
            let mut q = 0.0;
            let mut v = 0.0;
            for j in (0..=ORDER).rev() {
                q = q * z + self.q_coeff[j * n + i];
                v = v * z + self.v_coeff[j * n + i];
            }
            self.node_q[i] = finite(q)?;
            self.node_v[i] = finite(v)?;
        }
        Ok(())
    }

    pub(super) fn integrate_cell(
        &mut self,
        z: f64,
        time: f64,
        phase: LinearNetworkPhase,
        kinetic: f64,
        external: &[f64],
        step: &mut LinearNetworkEventStep,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        let n = self.q.len();
        let duration = finite(self.cell_h * z)?;
        // Bound every original state and connection force over the accepted
        // polynomial segment by its absolute coefficient sum.
        for i in 0..n {
            let mut qb = 0.0;
            let mut vb = 0.0;
            for j in (0..=ORDER).rev() {
                qb = qb * z + self.q_coeff[j * n + i].abs();
                vb = vb * z + self.v_coeff[j * n + i].abs();
            }
            if qb > self.config.maximum_abs_modal_coordinate
                || vb > self.config.maximum_abs_modal_velocity
            {
                return Err(Error::Budget("complete cell modal trajectory"));
            }
        }
        for (i, link) in self.links.iter().enumerate() {
            let mut bound = 0.0;
            for j in (0..=ORDER).rev() {
                let q = dot(&link.column, &self.q_coeff[j * n..(j + 1) * n])?
                    - if j == 0 { link.rest_extension_m } else { 0.0 };
                let v = dot(&link.column, &self.v_coeff[j * n..(j + 1) * n])?;
                bound = finite(
                    bound * z + (link.stiffness_n_per_m * q + link.damping_n_s_per_m * v).abs(),
                )?;
            }
            let column_norm = link.column.iter().map(|b| b.abs()).sum::<f64>();
            let uncertainty = column_norm
                * self.cell_state_error
                * (link.stiffness_n_per_m / self.cell_scale + link.damping_n_s_per_m)
                + 64.0 * f64::EPSILON * bound;
            self.link_force_bounds[i] = self.link_force_bounds[i].max(finite(bound + uncertainty)?);
        }
        const GAUSS: [(f64, f64); 8] = [
            (0.019855071751231884, 0.05061426814518813),
            (0.10166676129318664, 0.11119051722668724),
            (0.2372337950418355, 0.15685332293894366),
            (0.4082826787521751, 0.181341891689181),
            (0.5917173212478248, 0.181341891689181),
            (0.7627662049581645, 0.15685332293894366),
            (0.8983332387068134, 0.11119051722668724),
            (0.9801449282487681, 0.05061426814518813),
        ];
        let mut cell_modal_loss = 0.0;
        let mut cell_link_loss = 0.0;
        for (node, weight) in GAUSS {
            poll(gate)?;
            self.evaluate(z * node)?;
            let measure = duration * weight;
            cell_modal_loss = finite(
                cell_modal_loss
                    + measure
                        * self
                            .damping
                            .iter()
                            .zip(&self.node_v)
                            .map(|(c, v)| c * v * v)
                            .sum::<f64>(),
            )?;
            for (i, link) in self.links.iter().enumerate() {
                let speed = dot(&link.column, &self.node_v)?;
                let loss = finite(measure * link.damping_n_s_per_m * speed * speed)?;
                self.link_losses[i] = finite(self.link_losses[i] + loss)?;
                cell_link_loss = finite(cell_link_loss + loss)?;
            }
        }
        step.modal_damping_dissipation_j =
            finite(step.modal_damping_dissipation_j + cell_modal_loss)?;
        let mut integral = 0.0;
        let mut impulse = 0.0;
        let mut moment = 0.0;
        for j in (0..=ORDER).rev() {
            integral = integral * z + self.velocity_coeff[j] / (j + 1) as f64;
            impulse = impulse * z + self.reaction_coeff[j] / (j + 1) as f64;
            moment = moment * z + self.reaction_coeff[j] / (j + 2) as f64;
        }
        let displacement = finite(duration * integral)?;
        let (j, m) = match phase {
            LinearNetworkPhase::Sticking => (
                finite(duration * impulse)?,
                finite(self.cell_h * self.cell_h * z * z * moment)?,
            ),
            LinearNetworkPhase::Sliding { direction } => {
                let slip = finite(direction as f64 * displacement)?;
                // A chosen polynomial segment has no interior stop. Negative
                // signed travel means that event isolation is unresolved.
                if slip < -self.norm_b1 * self.cell_state_error / self.cell_scale {
                    return Err(Error::EventUnresolved);
                }
                let physical_slip = slip.max(0.0);
                step.accumulated_slip_m = finite(step.accumulated_slip_m + physical_slip)?;
                step.friction_dissipation_j =
                    finite(step.friction_dissipation_j + kinetic * physical_slip)?;
                (
                    finite(direction as f64 * kinetic * duration)?,
                    finite(direction as f64 * kinetic * 0.5 * duration * duration)?,
                )
            }
            LinearNetworkPhase::Fresh => return Err(Error::PhaseMismatch),
        };
        step.contact_impulse_n_s = finite(step.contact_impulse_n_s + j)?;
        step.contact_time_moment_n_s2 = finite(step.contact_time_moment_n_s2 + time * j + m)?;
        self.evaluate(z)?;
        let qmax = self
            .q
            .iter()
            .chain(&self.node_q)
            .map(|q| q.abs())
            .fold(0.0, f64::max);
        let vmax = self
            .v
            .iter()
            .chain(&self.node_v)
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        let mut rest_bound = 0.0;
        for link in &self.links {
            rest_bound = finite(
                rest_bound
                    + link.stiffness_n_per_m
                        * link.rest_extension_m.abs()
                        * link.column.iter().map(|b| b.abs()).fold(0.0, f64::max),
            )?;
        }
        let mut qerr = self.cell_state_error / self.cell_scale;
        let mut verr = self.cell_state_error;
        if phase == LinearNetworkPhase::Sticking {
            let drift = finite(
                self.tangent
                    .iter()
                    .zip(self.node_q.iter().zip(&self.q))
                    .map(|(b, (next, old))| b * (next - old))
                    .sum(),
            )?;
            let correction = finite(drift / self.norm_b2)?;
            let velocity = dot(&self.tangent, &self.node_v)? / self.norm_b2;
            for i in 0..n {
                let dq = self.tangent[i] * correction;
                let dv = self.tangent[i] * velocity;
                self.node_q[i] = finite(self.node_q[i] - dq)?;
                self.node_v[i] = finite(self.node_v[i] - dv)?;
                qerr = qerr.max(dq.abs());
                verr = verr.max(dv.abs());
            }
        }
        let energy_error = n as f64
            * (qerr * (self.stiffness_norm * qmax + rest_bound)
                + verr * vmax
                + 0.5 * (self.stiffness_norm * qerr * qerr + verr * verr));
        let work_error =
            external.iter().map(|f| f.abs()).sum::<f64>() * qerr + kinetic * self.norm_b1 * qerr;
        let resistance_error =
            duration * self.damping_norm * n as f64 * (2.0 * vmax * verr + verr * verr)
                + 256.0 * f64::EPSILON * (cell_modal_loss + cell_link_loss)
                + 1e-16
                    * duration
                    * self.damping_norm
                    * n as f64
                    * (vmax + self.cell_h * self.cell_rate * vmax).powi(2);
        step.numerical_work_uncertainty_j = finite(
            step.numerical_work_uncertainty_j + energy_error + work_error + resistance_error,
        )?;
        self.q.copy_from_slice(&self.node_q);
        self.v.copy_from_slice(&self.node_v);
        self.check_state()?;
        Ok(())
    }
}
