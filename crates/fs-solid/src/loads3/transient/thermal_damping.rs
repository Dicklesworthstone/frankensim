//! Passive Rayleigh damping while the source-resolved stiffness evolves.
//!
//! q1-q0 = h*vbar,
//! M*(v1-v0)/h + Cbar*vbar + (K0*q0+K1*q1)/2 = fbar,
//! Cbar = eta*M + zeta*(K0+K1)/2.
//!
//! Thus D = h*vbar'Cbar*vbar >= 0 for the admitted elastic operators.
//! This is NOT (C0*v0+C1*v1)/2, which can supply artificial work.
//! The original signed stiffness/eigenstrain accounts remain unchanged.
//! Endpoint acceleration is recovered from the actual endpoint equation, not
//! from the midpoint momentum increment. All solves use the existing PCG owner.

use super::*;

pub(super) struct Endpoint {
    pub(super) displacement: Vec<f64>,
    pub(super) velocity: Vec<f64>,
    pub(super) acceleration: Vec<f64>,
    pub(super) iterations: usize,
    pub(super) force_residual_n: f64,
}

/// M + h*Cbar/2 + h^2*K1/4, without an additional assembled matrix.
struct MidpointOperator<'a, 'cx> {
    mass: &'a Csr,
    old: &'a Csr,
    next: &'a Csr,
    cm: f64,
    c0: f64,
    c1: f64,
    cx: &'a Cx<'cx>,
}
impl LinearOp for MidpointOperator<'_, '_> {
    fn n(&self) -> usize { self.mass.nrows() }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let mut work = 0usize;
        for (row, out) in y.iter_mut().enumerate() {
            let mut value = 0.0;
            for (matrix, coefficient) in [
                (self.mass, self.cm), (self.old, self.c0), (self.next, self.c1),
            ] {
                if coefficient == 0.0 { continue; }
                let (columns, values) = matrix.row(row);
                let mut applied = 0.0;
                for (&column, &entry) in columns.iter().zip(values) {
                    if work % 256 == 0 && self.cx.checkpoint().is_err() {
                        // Poison the output; the checked solve polls before use.
                        *out = f64::NAN;
                        return;
                    }
                    applied = entry.mul_add(x[column], applied);
                    work += 1;
                }
                value = coefficient.mul_add(applied, value);
            }
            *out = value;
        }
    }
}

impl TetTransientRun3<'_, '_, '_, '_> {
    fn old_stiffness_action(&self, old: &Csr, x: &[f64], cx: &Cx<'_>)
        -> Result<Vec<f64>, TransientError3> {
        let mut result = vec![0.0; x.len()];
        Operator { m: &self.operators.mass, k: old, cm: 0.0, ck: 1.0, cx }
            .apply(x, &mut result);
        check(cx)?;
        if result.iter().any(|v| !v.is_finite()) {
            return Err(TransientError3::Invalid("nonfinite previous stiffness action"));
        }
        Ok(result)
    }

    pub(super) fn damped_material_endpoint(
        &self, old: &Csr, force: &[f64], thermal: &[f64], h: f64, b: f64, cx: &Cx<'_>,
    ) -> Result<Endpoint, TransientError3> {
        check(cx)?;
        let eta = self.config.mass_damping_per_s;
        let zeta = self.config.stiffness_damping_s;
        let q0 = &self.state.displacement;
        let v0 = &self.state.velocity;
        let old_q = self.old_stiffness_action(old, q0, cx)?;
        let next_q = self.apply(0.0, 1.0, q0, cx)?;
        let old_v = self.old_stiffness_action(old, v0, cx)?;
        let next_v = self.apply(0.0, 1.0, v0, cx)?;
        let mass_v = self.apply(1.0, 0.0, v0, cx)?;
        let mut rhs = Vec::with_capacity(v0.len());
        for i in 0..v0.len() {
            if i % 256 == 0 { check(cx)?; }
            rhs.push(number(0.5*force[i] + 0.5*self.force[i]
                + 0.5*thermal[i] + 0.5*self.thermal[i]
                - eta*mass_v[i] - zeta*(0.5*old_v[i] + 0.5*next_v[i])
                - 0.5*old_q[i] - 0.5*next_q[i] - 0.5*h*next_v[i])?);
        }
        let op = MidpointOperator {
            mass: &self.operators.mass, old, next: &self.operators.stiffness,
            cm: number(1.0 + 0.5*h*eta)?, c0: number(0.25*h*zeta)?,
            c1: number(b + 0.25*h*zeta)?, cx,
        };
        let (rate, first_iterations, first_residual) = solve_with_diagonal(
            &op, &rhs, self.config, cx,
            |i| op.cm*op.mass.get(i,i) + op.c0*op.old.get(i,i) + op.c1*op.next.get(i,i),
        )?;
        let mut q = Vec::with_capacity(v0.len());
        let mut v = Vec::with_capacity(v0.len());
        for i in 0..v0.len() {
            if i % 256 == 0 { check(cx)?; }
            let velocity = number(v0[i] + h*rate[i])?;
            q.push(number(q0[i] + h*(0.5*v0[i] + 0.5*velocity))?);
            v.push(velocity);
        }
        // The next step, reaction recovery and finite load events need physical
        // a1, not the midpoint rate. Recover it at K1 and C1 with the same M.
        let ku = self.apply(0.0, 1.0, &q, cx)?;
        let cv = self.apply(eta, zeta, &v, cx)?;
        for i in 0..rhs.len() {
            if i % 256 == 0 { check(cx)?; }
            rhs[i] = number(force[i] + thermal[i] - ku[i] - cv[i])?;
        }
        let mut remaining = self.config;
        remaining.solve.maximum_iterations = remaining.solve.maximum_iterations
            .checked_sub(first_iterations).ok_or(TransientError3::Budget)?;
        let (acceleration, second_iterations, second_residual) =
            solve(&self.operator(1.0, 0.0, cx), &rhs, remaining)?;
        check(cx)?;
        Ok(Endpoint {
            displacement: q, velocity: v, acceleration,
            iterations: first_iterations.checked_add(second_iterations)
                .ok_or(TransientError3::Budget)?,
            force_residual_n: first_residual.max(second_residual),
        })
    }

    pub(super) fn material_damping_loss(
        &self, old: &Csr, average_v: &[f64], h: f64, cx: &Cx<'_>,
    ) -> Result<f64, TransientError3> {
        let old_action = self.old_stiffness_action(old, average_v, cx)?;
        let next_action = self.apply(0.0, 1.0, average_v, cx)?;
        let mass_action = self.apply(1.0, 0.0, average_v, cx)?;
        let mut power = 0.0;
        for i in 0..average_v.len() {
            if i % 256 == 0 { check(cx)?; }
            let damping = number(self.config.mass_damping_per_s*mass_action[i]
                + self.config.stiffness_damping_s*(0.5*old_action[i] + 0.5*next_action[i]))?;
            power = number(average_v[i].mul_add(damping, power))?;
        }
        number(h*power)
    }
}
