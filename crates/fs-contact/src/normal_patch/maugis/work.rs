//! Branch work uses P(t) d(delta)/dt; every integrand evaluation is closed form.
//! Adaptive Simpson reports an error estimate, never an interval certificate.
use super::{MaugisError, SourceMaugisSphereContact, kernel};

/// Explicit absolute/relative work accuracy and deterministic evaluation budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisWorkBudget {
    /// Positive absolute tolerance, J.
    pub absolute_tolerance_j: f64,
    /// Nonnegative relative tolerance, smaller than one.
    pub relative_tolerance: f64,
    /// Maximum dimensionless integrand evaluations (at least five).
    pub maximum_evaluations: usize,
}

/// Reversible work along the contacting branch, excluding jumps to other branches.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisWorkEstimate {
    /// Signed integral F d(delta), J.
    pub work_j: f64,
    /// Sum of adaptive Simpson local error estimates, J; not a certified bound.
    pub estimated_absolute_error_j: f64,
    /// Integrand evaluations actually used.
    pub evaluations: usize,
}
#[derive(Clone, Copy)]
struct Panel {
    a: f64,
    b: f64,
    fa: f64,
    fm: f64,
    fb: f64,
    s: f64,
    tol: f64,
    depth: u8,
}
struct Integrator {
    lambda: f64,
    count: usize,
    limit: usize,
}
impl Integrator {
    fn f(&mut self, t: f64) -> Result<f64, MaugisError> {
        if self.count >= self.limit {
            return Err(MaugisError::WorkBudgetExceeded);
        }
        self.count += 1;
        let s = kernel::state(self.lambda, t)?;
        let f = s.p * s.dd;
        if !f.is_finite() {
            return Err(MaugisError::Unrepresentable {
                field: "branch work integrand",
            });
        }
        Ok(f)
    }
    fn integrate(&mut self, p: Panel) -> Result<(f64, f64), MaugisError> {
        let m = 0.5 * p.a + 0.5 * p.b;
        let l = 0.5 * p.a + 0.5 * m;
        let r = 0.5 * m + 0.5 * p.b;
        if l == p.a || l == m || r == m || r == p.b || p.depth >= 48 {
            return Err(MaugisError::WorkBudgetExceeded);
        }
        let fl = self.f(l)?;
        let fr = self.f(r)?;
        let left = (m - p.a) * (p.fa + 4.0 * fl + p.fm) / 6.0;
        let right = (p.b - m) * (p.fm + 4.0 * fr + p.fb) / 6.0;
        let delta = left + right - p.s;
        let error = delta.abs() / 15.0;
        if error <= p.tol {
            return Ok((left + right + delta / 15.0, error));
        }
        let (vl, el) = self.integrate(Panel {
            a: p.a,
            b: m,
            fa: p.fa,
            fm: fl,
            fb: p.fm,
            s: left,
            tol: p.tol * 0.5,
            depth: p.depth + 1,
        })?;
        let (vr, er) = self.integrate(Panel {
            a: m,
            b: p.b,
            fa: p.fm,
            fm: fr,
            fb: p.fb,
            s: right,
            tol: p.tol * 0.5,
            depth: p.depth + 1,
        })?;
        Ok((vl + vr, el + er))
    }
}
impl SourceMaugisSphereContact {
    /// Integrate the original force over a smooth stable contact branch. This
    /// supplies work differences, not a detached-reference absolute potential.
    pub fn reversible_work_between(
        &self,
        from_m: f64,
        to_m: f64,
        budget: MaugisWorkBudget,
    ) -> Result<MaugisWorkEstimate, MaugisError> {
        if !budget.absolute_tolerance_j.is_finite()
            || budget.absolute_tolerance_j <= 0.0
            || !budget.relative_tolerance.is_finite()
            || budget.relative_tolerance < 0.0
            || budget.relative_tolerance >= 1.0
            || budget.maximum_evaluations < 5
        {
            return Err(MaugisError::InvalidInput {
                field: "branch work budget",
            });
        }
        let from = self.state_at_indentation(from_m)?;
        let to = self.state_at_indentation(to_m)?;
        if from_m == to_m {
            return Ok(MaugisWorkEstimate {
                work_j: 0.0,
                estimated_absolute_error_j: 0.0,
                evaluations: 0,
            });
        }
        if from.t == to.t {
            return Err(MaugisError::Unrepresentable {
                field: "distinct work endpoints collapse",
            });
        }
        let scale = self.force_scale_n * self.indentation_scale_m;
        if !scale.is_finite() || scale <= 0.0 {
            return Err(MaugisError::Unrepresentable {
                field: "branch work scale",
            });
        }
        let sign = if to.t > from.t { 1.0 } else { -1.0 };
        let a = from.t.min(to.t);
        let b = from.t.max(to.t);
        let mut q = Integrator {
            lambda: self.lambda,
            count: 0,
            limit: budget.maximum_evaluations,
        };
        let fa = q.f(a)?;
        let fm = q.f(0.5 * a + 0.5 * b)?;
        let fb = q.f(b)?;
        let s = (b - a) * (fa + 4.0 * fm + fb) / 6.0;
        let tol = budget.absolute_tolerance_j / scale + budget.relative_tolerance * s.abs();
        if !tol.is_finite() || tol <= 0.0 {
            return Err(MaugisError::Unrepresentable {
                field: "scaled branch work tolerance",
            });
        }
        let (value, error) = q.integrate(Panel {
            a,
            b,
            fa,
            fm,
            fb,
            s,
            tol,
            depth: 0,
        })?;
        let work_j = sign * value * scale;
        let estimated_absolute_error_j = error * scale;
        if !work_j.is_finite() || !estimated_absolute_error_j.is_finite() {
            return Err(MaugisError::Unrepresentable {
                field: "branch work result",
            });
        }
        Ok(MaugisWorkEstimate {
            work_j,
            estimated_absolute_error_j,
            evaluations: q.count,
        })
    }
}
