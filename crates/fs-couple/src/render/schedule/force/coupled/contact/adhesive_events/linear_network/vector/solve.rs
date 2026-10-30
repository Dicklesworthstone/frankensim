//! Original full operator condensation, disk residual and discrete work account.
use super::{
    Error, ROWS, VectorFrictionRegime, VectorLinearNetworkFlow, VectorLinearNetworkStep, dot,
    finite, norm2, poll,
};
use fs_exec::CancelGate;

impl VectorLinearNetworkFlow {
    pub(super) fn prepare(&mut self, h: f64, gate: &CancelGate) -> Result<(), Error> {
        if self.cached_h == Some(h) {
            return Ok(());
        }
        self.cached_h = None;
        let n = self.q.len();
        let m = self.columns.len();
        self.matrix.fill(0.0);
        for i in 0..n {
            self.matrix[i * n + i] = finite(1.0 + h * self.damping[i] + h * h * self.stiffness[i])?;
        }
        for l in &self.links {
            poll(gate)?;
            let coefficient = finite(h * l.damping_n_s_per_m + h * h * l.stiffness_n_per_m)?;
            for i in 0..n {
                for j in 0..n {
                    self.matrix[i * n + j] =
                        finite(self.matrix[i * n + j] + coefficient * l.column[i] * l.column[j])?;
                }
            }
        }
        poll(gate)?;
        self.mechanical.factor(&self.matrix)?;
        poll(gate)?;
        for (j, b) in self.columns.iter().enumerate() {
            self.mechanical
                .solve_factored_into(b, &mut self.response_columns[j * n..(j + 1) * n])?;
            poll(gate)?;
        }
        for i in 0..m {
            for j in 0..m {
                self.mobility[i * m + j] =
                    finite(h * dot(&self.columns[i], &self.response_columns[j * n..(j + 1) * n])?)?;
            }
        }
        // Symmetry belongs to the original reciprocal mechanical operator.
        // The restored uncondensed equation is independently audited each cell.
        for i in 0..m {
            for j in 0..i {
                let mean = finite(0.5 * (self.mobility[i * m + j] + self.mobility[j * m + i]))?;
                self.mobility[i * m + j] = mean;
                self.mobility[j * m + i] = mean;
            }
        }
        self.lipschitz = 0.0;
        for i in 0..m {
            let sum = finite(
                self.mobility[i * m..(i + 1) * m]
                    .iter()
                    .map(|v| v.abs())
                    .sum(),
            )?;
            self.lipschitz = self.lipschitz.max(sum);
        }
        if self.lipschitz <= 0.0 {
            return Err(Error::Invalid("positive tangent mobility"));
        }
        self.cached_h = Some(h);
        Ok(())
    }

    pub(super) fn cell(
        &mut self,
        external: &[f64],
        h: f64,
        capacities: &[f64],
        time: f64,
        step: &mut VectorLinearNetworkStep,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        let n = self.q.len();
        let m = self.columns.len();
        for i in 0..n {
            self.rhs[i] = finite(self.v[i] + h * (external[i] - self.stiffness[i] * self.q[i]))?;
        }
        for l in &self.links {
            let force =
                finite(l.stiffness_n_per_m * (dot(&l.column, &self.q)? - l.rest_extension_m))?;
            for i in 0..n {
                self.rhs[i] = finite(self.rhs[i] - h * l.column[i] * force)?;
            }
        }
        self.mechanical
            .solve_factored_into(&self.rhs, &mut self.free_v)?;
        for i in 0..m {
            self.free_slip[i] = dot(&self.columns[i], &self.free_v)?;
        }
        self.solve_disks(capacities, step, gate)?;
        for i in 0..n {
            let mut v = self.free_v[i];
            for j in 0..m {
                v -= h * self.response_columns[j * n + i] * self.force[j];
            }
            self.next_v[i] = finite(v)?;
        }
        // Audit the actual complete original equation, not only its condensed
        // contact residual. Its work defect enters numerical uncertainty.
        let mut equation_work_error = 0.0;
        let mut arithmetic_scale = 0.0;
        for i in 0..n {
            let mut force = 0.0;
            for j in 0..m {
                force += self.columns[j][i] * self.force[j];
            }
            let av = dot(&self.matrix[i * n..(i + 1) * n], &self.next_v)?;
            let defect = finite(av - self.rhs[i] + h * force)?;
            equation_work_error += defect.abs() * self.next_v[i].abs();
            arithmetic_scale += av.abs() * self.next_v[i].abs()
                + self.rhs[i].abs() * self.next_v[i].abs()
                + h * force.abs() * self.next_v[i].abs();
        }
        let mut modal_loss = 0.0;
        let mut numerical_loss = 0.0;
        for i in 0..n {
            let dq = h * self.next_v[i];
            let dv = self.next_v[i] - self.v[i];
            modal_loss += h * self.damping[i] * self.next_v[i] * self.next_v[i];
            numerical_loss += 0.5 * (dv * dv + self.stiffness[i] * dq * dq);
        }
        for (i, l) in self.links.iter().enumerate() {
            let speed = dot(&l.column, &self.next_v)?;
            let loss = finite(h * l.damping_n_s_per_m * speed * speed)?;
            let numerical = finite(0.5 * l.stiffness_n_per_m * (h * speed) * (h * speed))?;
            self.link_losses[i] = finite(self.link_losses[i] + loss)?;
            self.link_numerical_losses[i] = finite(self.link_numerical_losses[i] + numerical)?;
            numerical_loss += numerical;
        }
        let mut law_work_error = 0.0;
        for j in 0..self.contact_count() {
            let u = [
                dot(&self.columns[2 * j], &self.next_v)?,
                dot(&self.columns[2 * j + 1], &self.next_v)?,
            ];
            let t = [self.force[2 * j], self.force[2 * j + 1]];
            let speed = norm2(u);
            let traction = norm2(t);
            let ftol = self.config.force_absolute_tolerance_n
                + self.config.force_relative_tolerance * capacities[j];
            if traction > capacities[j] + ftol {
                return Err(Error::NotConverged);
            }
            let work = finite(h * (t[0] * u[0] + t[1] * u[1]))?;
            let friction = finite(h * capacities[j] * speed)?;
            law_work_error += (friction - work).abs();
            // A solved interior force is a sticking disk multiplier. Its tiny
            // velocity residual is measured, never erased by a velocity clip.
            step.regimes[j] = if capacities[j] == 0.0 {
                VectorFrictionRegime::Inactive
            } else if speed <= self.config.velocity_absolute_tolerance_m_s
                && traction <= capacities[j] + ftol
            {
                VectorFrictionRegime::Sticking
            } else {
                VectorFrictionRegime::Sliding
            };
            step.contact_friction_dissipation_j[j] =
                finite(step.contact_friction_dissipation_j[j] + friction)?;
            step.accumulated_slip_m[j] = finite(step.accumulated_slip_m[j] + h * speed)?;
            step.maximum_contact_speed_m_s[j] = step.maximum_contact_speed_m_s[j].max(speed);
            step.maximum_contact_force_n[j] = step.maximum_contact_force_n[j].max(traction);
            for a in 0..2 {
                let impulse = h * t[a];
                step.contact_impulse_n_s[j][a] = finite(step.contact_impulse_n_s[j][a] + impulse)?;
                step.contact_time_moment_n_s2[j][a] =
                    finite(step.contact_time_moment_n_s2[j][a] + (time + 0.5 * h) * impulse)?;
            }
        }
        step.modal_damping_dissipation_j = finite(step.modal_damping_dissipation_j + modal_loss)?;
        step.numerical_dissipation_j = finite(step.numerical_dissipation_j + numerical_loss)?;
        step.numerical_work_uncertainty_j = finite(
            step.numerical_work_uncertainty_j
                + equation_work_error
                + law_work_error
                + 32.0
                    * f64::EPSILON
                    * (n + self.links.len() + m) as f64
                    * (arithmetic_scale + modal_loss + numerical_loss + law_work_error),
        )?;
        for i in 0..n {
            self.q[i] = finite(self.q[i] + h * self.next_v[i])?;
            self.v[i] = self.next_v[i];
        }
        poll(gate)?;
        Ok(())
    }

    fn solve_disks(
        &mut self,
        capacities: &[f64],
        step: &mut VectorLinearNetworkStep,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        let m = self.columns.len();
        loop {
            poll(gate)?;
            let norm = residual(
                &self.mobility,
                self.lipschitz,
                &self.free_slip,
                &self.force,
                capacities,
                &mut self.residual,
                Some(&mut self.jacobian),
            )?;
            let mut converged = true;
            for (j, capacity) in capacities.iter().enumerate() {
                let r = norm2([self.residual[2 * j], self.residual[2 * j + 1]]);
                let ftol = self.config.force_absolute_tolerance_n
                    + self.config.force_relative_tolerance * capacity;
                converged &=
                    r <= ftol && self.lipschitz * r <= self.config.velocity_absolute_tolerance_m_s;
            }
            if converged {
                return Ok(());
            }
            if step.nonlinear_iterations >= self.config.maximum_nonlinear_iterations {
                return Err(Error::Budget("joint disk nonlinear iterations"));
            }
            step.nonlinear_iterations += 1;
            for i in 0..m {
                self.trial_residual[i] = -self.residual[i];
            }
            self.contact_solve.solve_into(
                &self.jacobian[..m * m],
                &self.trial_residual[..m],
                &mut self.direction[..m],
            )?;
            let mut accepted = false;
            let mut alpha = 1.0;
            for _ in 0..self.config.maximum_backtracks {
                poll(gate)?;
                for i in 0..m {
                    self.trial_force[i] = finite(self.force[i] + alpha * self.direction[i])?;
                }
                let candidate = residual(
                    &self.mobility,
                    self.lipschitz,
                    &self.free_slip,
                    &self.trial_force,
                    capacities,
                    &mut self.trial_residual,
                    None,
                )?;
                if candidate < norm * (1.0 - 1e-4 * alpha) {
                    self.force[..m].copy_from_slice(&self.trial_force[..m]);
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            if !accepted {
                // The gradient-projection map of an SPD quadratic with this
                // row-sum step is nonexpansive. This fallback is a genuine
                // disk iteration, not a regularized friction law.
                for i in 0..m {
                    self.trial_force[i] = self.force[i] - self.residual[i];
                }
                let candidate = residual(
                    &self.mobility,
                    self.lipschitz,
                    &self.free_slip,
                    &self.trial_force,
                    capacities,
                    &mut self.trial_residual,
                    None,
                )?;
                if candidate >= norm {
                    return Err(Error::NotConverged);
                }
                self.force[..m].copy_from_slice(&self.trial_force[..m]);
            }
        }
    }
}

/// R(t)=t-proj_disk(t+(w-Ht)/L), with its analytic block derivative.
fn residual(
    h: &[f64; ROWS * ROWS],
    l: f64,
    w: &[f64; ROWS],
    t: &[f64; ROWS],
    capacities: &[f64],
    out: &mut [f64; ROWS],
    mut jacobian: Option<&mut [f64; ROWS * ROWS]>,
) -> Result<f64, Error> {
    let m = 2 * capacities.len();
    let mut norm = 0.0;
    for (j, capacity) in capacities.iter().enumerate() {
        let rows = [2 * j, 2 * j + 1];
        let mut z = [0.0; 2];
        for a in 0..2 {
            let row = rows[a];
            let mut u = w[row];
            for k in 0..m {
                u -= h[row * m + k] * t[k];
            }
            z[a] = finite(t[row] + u / l)?;
        }
        let size = finite(norm2(z))?;
        let (p, derivative) = if *capacity == 0.0 {
            ([0.0; 2], [[0.0; 2]; 2])
        } else if size <= *capacity {
            (z, [[1.0, 0.0], [0.0, 1.0]])
        } else {
            let d = [z[0] / size, z[1] / size];
            let scale = capacity / size;
            (
                [capacity * d[0], capacity * d[1]],
                [
                    [scale * (1.0 - d[0] * d[0]), -scale * d[0] * d[1]],
                    [-scale * d[1] * d[0], scale * (1.0 - d[1] * d[1])],
                ],
            )
        };
        for a in 0..2 {
            let row = rows[a];
            out[row] = finite(t[row] - p[a])?;
            norm += out[row] * out[row];
            if let Some(jac) = jacobian.as_deref_mut() {
                for k in 0..m {
                    let mut value = if row == k { 1.0 } else { 0.0 };
                    for b in 0..2 {
                        value -= derivative[a][b]
                            * ((if rows[b] == k { 1.0 } else { 0.0 }) - h[rows[b] * m + k] / l);
                    }
                    jac[row * m + k] = finite(value)?;
                }
            }
        }
    }
    finite(norm)
}
