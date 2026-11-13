//! Original-coordinate contractions; no fit, inverse, integrator or storage law.

use fs_math::dd::Dd;

/// Immutable copies of the source realization and the physical connection.
/// Construction/shape/resource admission belongs to the enclosing adapter.
pub(super) struct Boundary {
    pub host: usize,
    pub states: usize,
    pub ports: usize,
    pub columns: Vec<f64>,
    pub a: Vec<f64>,
    pub b: Vec<f64>,
    pub c: Vec<f64>,
    pub d: Vec<f64>,
}

fn product(sum: Dd, a: f64, b: f64) -> Dd {
    if a == 0.0 || b == 0.0 {
        sum
    } else {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    }
}

impl Boundary {
    /// Constant source-to-host tangent row -S C, prepared once. Original
    /// source A is copied directly into the joint Jacobian, never rediscovered
    /// by repeated basis evaluations of a dense q-by-q state rate.
    pub fn state_reaction_row(&self, row: usize, out: &mut [f64]) -> bool {
        if row >= self.host || out.len() != self.states {
            return false;
        }
        for (column, value) in out.iter_mut().enumerate() {
            let mut sum = Dd::ZERO;
            for port in 0..self.ports {
                sum = product(
                    sum,
                    -self.columns[row * self.ports + port],
                    self.c[port * self.states + column],
                );
            }
            *value = sum.to_f64();
        }
        out.iter().all(|v| v.is_finite())
    }

    /// y = S^T grad H in the original source port order.
    pub fn flow_into(&self, effort: &[f64], flow: &mut [f64]) -> bool {
        if effort.len() != self.host || flow.len() != self.ports {
            return false;
        }
        for (port, value) in flow.iter_mut().enumerate() {
            let mut sum = Dd::ZERO;
            for (row, &e) in effort.iter().enumerate() {
                sum = product(sum, self.columns[row * self.ports + port], e);
            }
            *value = sum.to_f64();
        }
        flow.iter().all(|v| v.is_finite())
    }

    /// Opposing boundary effort C z + D y. Singular/zero D needs no special case.
    pub fn effort_into(&self, state: &[f64], flow: &[f64], out: &mut [f64]) -> bool {
        if state.len() != self.states || flow.len() != self.ports || out.len() != self.ports {
            return false;
        }
        for (row, value) in out.iter_mut().enumerate() {
            let mut sum = Dd::ZERO;
            for (column, &z) in state.iter().enumerate() {
                sum = product(sum, self.c[row * self.states + column], z);
            }
            for (column, &y) in flow.iter().enumerate() {
                sum = product(sum, self.d[row * self.ports + column], y);
            }
            *value = sum.to_f64();
        }
        out.iter().all(|v| v.is_finite())
    }

    /// D dy for an original host-state perturbation; avoids visiting source
    /// coordinates that are identically zero in this Jacobian block.
    pub fn direct_effort_into(&self, flow: &[f64], out: &mut [f64]) -> bool {
        if flow.len() != self.ports || out.len() != self.ports {
            return false;
        }
        for (row, value) in out.iter_mut().enumerate() {
            let mut sum = Dd::ZERO;
            for (column, &y) in flow.iter().enumerate() {
                sum = product(sum, self.d[row * self.ports + column], y);
            }
            *value = sum.to_f64();
        }
        out.iter().all(|v| v.is_finite())
    }

    /// B dy in the original source coordinates, without a zero-state A sweep.
    pub fn driven_rate_into(&self, flow: &[f64], out: &mut [f64]) -> bool {
        if flow.len() != self.ports || out.len() != self.states {
            return false;
        }
        for (row, value) in out.iter_mut().enumerate() {
            let mut sum = Dd::ZERO;
            for (column, &y) in flow.iter().enumerate() {
                sum = product(sum, self.b[row * self.ports + column], y);
            }
            *value = sum.to_f64();
        }
        out.iter().all(|v| v.is_finite())
    }

    /// Every original pole coordinate evolves by A z + B y.
    pub fn rate_into(&self, state: &[f64], flow: &[f64], out: &mut [f64]) -> bool {
        if state.len() != self.states || flow.len() != self.ports || out.len() != self.states {
            return false;
        }
        for (row, value) in out.iter_mut().enumerate() {
            let mut sum = Dd::ZERO;
            for (column, &z) in state.iter().enumerate() {
                sum = product(sum, self.a[row * self.states + column], z);
            }
            for (column, &y) in flow.iter().enumerate() {
                sum = product(sum, self.b[row * self.ports + column], y);
            }
            *value = sum.to_f64();
        }
        out.iter().all(|v| v.is_finite())
    }

    /// Host-side boundary work: retain both words before the effort contraction.
    /// This is intentionally independent of the rounded flow seen by the load.
    pub fn host_power(&self, host_effort: &[f64], boundary_effort: &[f64]) -> f64 {
        let mut sum = Dd::ZERO;
        for (row, &h) in host_effort.iter().enumerate() {
            let mut reaction = Dd::ZERO;
            for (port, &e) in boundary_effort.iter().enumerate() {
                reaction = product(reaction, self.columns[row * self.ports + port], e);
            }
            sum = sum + reaction * Dd::from_f64(h);
        }
        sum.to_f64()
    }
}

pub(super) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .fold(Dd::ZERO, |sum, (&a, &b)| product(sum, a, b))
        .to_f64()
}

#[cfg(test)]
#[path = "operator/tests.rs"]
mod tests;
