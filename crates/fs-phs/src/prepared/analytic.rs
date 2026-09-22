//! Analytic derivative of the SAME Gonzalez discrete gradient, not midpoint
//! linearization. All scratch belongs to StepWorkspace; physical storage and
//! candidate outputs are untouched. The Jacobian is generally nonsymmetric.
use super::{StepWorkspace, PreparedStepError, dimensions, norm};
use crate::PortHamiltonian;
use super::dissipation::Dissipation;

impl StepWorkspace {
    #[allow(clippy::type_complexity)]
    pub(super) fn analytic_jacobian_into<F>(
        &mut self, sys: &PortHamiltonian, x0: &[f64], dt: f64,
        hessian: &dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
        dissipation: Option<Dissipation<'_>>, poll: &mut F,
    ) -> Result<(), PreparedStepError>
    where F: FnMut() -> Result<(), PreparedStepError> {
        let n = self.n;
        for i in 0..n {
            self.midpoint[i] = f64::midpoint(x0[i], self.x[i]);
            self.delta[i] = self.x[i] - x0[i];
        }
        self.effort.fill(0.0);
        sys.storage.gradient(&self.midpoint, &mut self.effort);
        norm(&self.effort)?;
        let length2 = self.delta.iter().map(|d| d*d).sum::<f64>();
        let scale = x0.iter().chain(&self.x)
            .fold(f64::MIN_POSITIVE, |s, &x| s.max(x.abs()));
        let floor = 1.0e-14 * scale;
        // Exactly the same cancellation guard as the canonical DG kernel.
        let corrected = length2 > floor*floor;
        let alpha = if corrected {
            self.minus.fill(0.0);
            sys.storage.gradient(&self.x, &mut self.minus);
            norm(&self.minus)?;
            let mid_dot = self.effort.iter().zip(&self.delta).map(|(g,d)| g*d).sum::<f64>();
            (sys.hamiltonian(&self.x) - sys.hamiltonian(x0) - mid_dot) / length2
        } else { 0.0 };
        norm(&[length2, alpha])?;
        if dissipation.is_some() {
            // A nonlinear port is conjugate to the COMPLETE Gonzalez effort,
            // not just grad H at the midpoint. Keep the latter for d(alpha).
            for i in 0..n { self.nonlinear_effort[i] = self.effort[i] + alpha*self.delta[i]; }
            norm(&self.nonlinear_effort)?;
        }
        self.trial.fill(0.0);
        for col in 0..n {
            poll()?;
            // One analytic Hessian action instead of two whole residual probes.
            self.trial[col] = 1.0;
            self.plus.fill(0.0);
            if !hessian(&self.midpoint, &self.trial, &mut self.plus) {
                return Err(dimensions("storage lacks an analytic Hessian-vector action").into());
            }
            self.trial[col] = 0.0;
            norm(&self.plus)?;
            let derivative = if corrected {
                let h_dot = self.plus.iter().zip(&self.delta).map(|(h,d)| h*d).sum::<f64>();
                (self.minus[col] - self.effort[col] - 0.5*h_dot
                    - 2.0*alpha*self.delta[col]) / length2
            } else { 0.0 };
            for row in 0..n {
                self.plus[row] = 0.5*self.plus[row] + self.delta[row]*derivative
                    + if row == col { alpha } else { 0.0 };
            }
            norm(&self.plus)?;
            if let Some(port) = dissipation {
                self.trial[col] = 0.5; // d(midpoint)/d(x1_col)
                let result = port.directional(&self.midpoint, &self.nonlinear_effort,
                    &self.trial, &self.plus, &mut self.nonlinear_loss);
                self.trial[col] = 0.0;
                result?;
            }
            for row in 0..n {
                let mut flow = 0.0;
                for &k in self.flow.row(row) {
                    flow += (sys.j[row*n+k] - sys.r[row*n+k]) * self.plus[k];
                }
                if dissipation.is_some() { flow -= self.nonlinear_loss[row]; }
                self.jacobian[row*n+col] = (if row == col { 1.0 } else { 0.0 }) - dt*flow;
            }
        }
        norm(&self.jacobian)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
