//! The gradient gate primitive (plan §8.7: "a solver without a passing
//! gradient check cannot merge"). [`gradcheck`] compares the forward-dual
//! gradient against central finite differences with scale-aware
//! tolerances and returns a structured verdict — CI (the merge-gate bead)
//! calls this and fails the lane on `pass == false`.
//!
//! What it catches: exactly the derivative-killing bug class — value
//! round-trips through `Real::value()`/`from_f64` (which silently zero
//! the dual channel), hand-written vjps that drift from the primal, and
//! branch conventions that disagree across paths. Tested on a seeded
//! specimen of that class.

use crate::dual::Dual64;

/// Verdict + evidence from a gradient check.
#[derive(Debug, Clone)]
pub struct GradCheckReport {
    /// The dual (analytic) gradient.
    pub grad: Vec<f64>,
    /// Central finite-difference gradient.
    pub fd: Vec<f64>,
    /// max over components of |grad − fd| / max(|grad|, |fd|, 1).
    pub max_rel_err: f64,
    /// Component achieving `max_rel_err`.
    pub worst_index: usize,
    /// Tolerance used.
    pub tol: f64,
    /// The gate verdict.
    pub pass: bool,
}

impl core::fmt::Display for GradCheckReport {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "{{\"check\":\"gradcheck\",\"pass\":{},\"max_rel_err\":{:.3e},\"worst_index\":{},\"tol\":{:.1e}}}",
            self.pass, self.max_rel_err, self.worst_index, self.tol
        )
    }
}

/// Check the dual gradient of `f` at `x` against central differences.
///
/// `f` must be generic over [`Real`] in the sense that a SINGLE closure
/// evaluates both the dual gradient (N-lane seeding) and the f64 FD
/// samples (constant duals; the primal channel is bit-exact, fs-ad's
/// proven property). `tol` is the scale-aware relative tolerance —
/// 5e-6 is a sound default for central differences at h = cbrt(eps).
#[must_use]
pub fn gradcheck<const N: usize, F>(f: F, x: [f64; N], tol: f64) -> GradCheckReport
where
    F: Fn([Dual64<N>; N]) -> Dual64<N>,
{
    // Dual gradient in one N-lane evaluation.
    let mut vars = [Dual64::<N>::constant(0.0); N];
    for (i, v) in vars.iter_mut().enumerate() {
        *v = Dual64::variable(x[i], i);
    }
    let out = f(vars);
    let grad: Vec<f64> = out.eps.to_vec();
    // Central differences at the classic h ~ cbrt(eps)·scale.
    let eval = |xs: [f64; N]| -> f64 {
        let mut c = [Dual64::<N>::constant(0.0); N];
        for (i, v) in c.iter_mut().enumerate() {
            *v = Dual64::constant(xs[i]);
        }
        f(c).re
    };
    let mut fd = vec![0.0f64; N];
    for k in 0..N {
        let h = 6e-6 * x[k].abs().max(1.0);
        let mut xp = x;
        let mut xm = x;
        xp[k] += h;
        xm[k] -= h;
        fd[k] = (eval(xp) - eval(xm)) / (2.0 * h);
    }
    let (mut max_rel_err, mut worst_index) = (0.0f64, 0usize);
    for k in 0..N {
        let scale = grad[k].abs().max(fd[k].abs()).max(1.0);
        let rel = (grad[k] - fd[k]).abs() / scale;
        if rel > max_rel_err {
            max_rel_err = rel;
            worst_index = k;
        }
    }
    GradCheckReport { grad, fd, max_rel_err, worst_index, tol, pass: max_rel_err <= tol }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Real;

    #[test]
    fn correct_gradient_passes() {
        let rep = gradcheck(
            |[x, y]: [Dual64<2>; 2]| (x * y).sin() + x.exp() * y.tanh(),
            [0.7, -0.4],
            5e-6,
        );
        assert!(rep.pass, "correct derivative must pass: {rep}");
        assert!(rep.max_rel_err < 1e-8, "central FD should be much tighter: {rep}");
    }

    #[test]
    fn derivative_killing_roundtrip_is_caught() {
        // The classic bug: routing a value through value()/from_f64 zeroes
        // the dual channel — the primal is unchanged (FD sees the true
        // slope) but the dual gradient is wrong. The gate MUST fail this.
        let rep = gradcheck(
            |[x]: [Dual64<1>; 1]| {
                let leaked = Dual64::<1>::from_f64(x.value()); // dual lost!
                leaked * x // dual grad = x, true grad = 2x
            },
            [1.3],
            5e-6,
        );
        assert!(!rep.pass, "the round-trip bug must be caught: {rep}");
        assert!(rep.max_rel_err > 0.3, "the error is O(1), not noise: {rep}");
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"gradcheck-gate\",\"verdict\":\"pass\",\"detail\":\"kills value() round-trip bug: {rep}\"}}"
        );
    }

    #[test]
    fn report_serializes_as_json_line() {
        let rep = gradcheck(|[x]: [Dual64<1>; 1]| x * x, [2.0], 5e-6);
        let line = format!("{rep}");
        assert!(line.starts_with('{') && line.ends_with('}'), "JSON-ish line: {line}");
        assert!(line.contains("\"pass\":true"));
    }
}
