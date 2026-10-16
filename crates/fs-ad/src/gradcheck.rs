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
    /// The dual (analytic) gradient, or its validated finite prefix on refusal.
    pub grad: Vec<f64>,
    /// Central finite-difference gradient, or its validated finite prefix on
    /// refusal.
    pub fd: Vec<f64>,
    /// max over components of |grad − fd| / max(|grad|, |fd|, 1), or
    /// positive infinity when the gate refuses invalid evidence.
    pub max_rel_err: f64,
    /// Component achieving `max_rel_err`, or the first component whose
    /// evidence refused. [`usize::MAX`] denotes a global refusal before any
    /// component applied.
    pub worst_index: usize,
    /// Tolerance used.
    pub tol: f64,
    /// The gate verdict.
    pub pass: bool,
}

fn json_f64(value: f64, precision: usize) -> String {
    if value.is_finite() {
        format!("{value:.precision$e}")
    } else {
        "null".to_owned()
    }
}

fn json_index(value: usize) -> String {
    if value == usize::MAX {
        "null".to_owned()
    } else {
        value.to_string()
    }
}

fn fail_closed(grad: Vec<f64>, fd: Vec<f64>, worst_index: usize, tol: f64) -> GradCheckReport {
    GradCheckReport {
        grad,
        fd,
        max_rel_err: f64::INFINITY,
        worst_index,
        tol,
        pass: false,
    }
}

impl core::fmt::Display for GradCheckReport {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "{{\"check\":\"gradcheck\",\"pass\":{},\"max_rel_err\":{},\"worst_index\":{},\"tol\":{}}}",
            self.pass,
            json_f64(self.max_rel_err, 3),
            json_index(self.worst_index),
            json_f64(self.tol, 1)
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
/// Empty checks, negative/non-finite tolerances, non-finite points,
/// callback evidence, perturbations, or arithmetic intermediates fail closed.
#[must_use]
pub fn gradcheck<const N: usize, F>(f: F, x: [f64; N], tol: f64) -> GradCheckReport
where
    F: Fn([Dual64<N>; N]) -> Dual64<N>,
{
    if N == 0 || !tol.is_finite() || tol < 0.0 {
        return fail_closed(Vec::new(), Vec::new(), usize::MAX, tol);
    }
    if let Some(index) = x.iter().position(|value| !value.is_finite()) {
        return fail_closed(Vec::new(), Vec::new(), index, tol);
    }

    // Dual gradient in one N-lane evaluation.
    let mut vars = [Dual64::<N>::constant(0.0); N];
    for (i, v) in vars.iter_mut().enumerate() {
        *v = Dual64::variable(x[i], i);
    }
    let out = f(vars);
    if !out.re.is_finite() {
        return fail_closed(Vec::new(), Vec::new(), usize::MAX, tol);
    }
    let mut grad = Vec::with_capacity(N);
    for (index, value) in out.eps.into_iter().enumerate() {
        if !value.is_finite() {
            return fail_closed(grad, Vec::new(), index, tol);
        }
        grad.push(value);
    }
    // Central differences at the classic h ~ cbrt(eps)·scale.
    let eval = |xs: [f64; N]| -> Dual64<N> {
        let mut c = [Dual64::<N>::constant(0.0); N];
        for (i, v) in c.iter_mut().enumerate() {
            *v = Dual64::constant(xs[i]);
        }
        f(c)
    };
    let mut fd = Vec::with_capacity(N);
    for k in 0..N {
        let coordinate_scale = x[k].abs().max(1.0);
        let h = 6e-6 * coordinate_scale;
        let denominator = 2.0 * h;
        if !coordinate_scale.is_finite()
            || !h.is_finite()
            || h <= 0.0
            || !denominator.is_finite()
            || denominator <= 0.0
        {
            return fail_closed(grad, fd, k, tol);
        }
        let mut xp = x;
        let mut xm = x;
        xp[k] += h;
        xm[k] -= h;
        if !xp[k].is_finite()
            || !xm[k].is_finite()
            || xp[k].to_bits() == x[k].to_bits()
            || xm[k].to_bits() == x[k].to_bits()
            || xp[k].to_bits() == xm[k].to_bits()
        {
            return fail_closed(grad, fd, k, tol);
        }
        let plus = eval(xp);
        if !plus.re.is_finite() {
            return fail_closed(grad, fd, k, tol);
        }
        let minus = eval(xm);
        if !minus.re.is_finite() {
            return fail_closed(grad, fd, k, tol);
        }
        let numerator = plus.re - minus.re;
        if !numerator.is_finite() {
            return fail_closed(grad, fd, k, tol);
        }
        let derivative = numerator / denominator;
        if !derivative.is_finite() {
            return fail_closed(grad, fd, k, tol);
        }
        fd.push(derivative);
    }
    let (mut max_rel_err, mut worst_index) = (0.0f64, 0usize);
    for k in 0..N {
        let scale = grad[k].abs().max(fd[k].abs()).max(1.0);
        let error = (grad[k] - fd[k]).abs();
        if !scale.is_finite() || !error.is_finite() {
            return fail_closed(grad, fd, k, tol);
        }
        let rel = error / scale;
        if !rel.is_finite() {
            return fail_closed(grad, fd, k, tol);
        }
        if rel > max_rel_err {
            max_rel_err = rel;
            worst_index = k;
        }
    }
    let complete = grad.len() == N && fd.len() == N;
    GradCheckReport {
        grad,
        fd,
        max_rel_err,
        worst_index,
        tol,
        pass: complete && max_rel_err.is_finite() && max_rel_err <= tol,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Real;

    fn assert_refused(report: &GradCheckReport) {
        assert!(!report.pass, "invalid evidence must fail closed: {report}");
        assert!(report.max_rel_err.is_infinite());
        assert!(report.max_rel_err.is_sign_positive());
        assert!(report.grad.iter().all(|value| value.is_finite()));
        assert!(report.fd.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn correct_gradient_passes() {
        let rep = gradcheck(
            |[x, y]: [Dual64<2>; 2]| (x * y).sin() + x.exp() * y.tanh(),
            [0.7, -0.4],
            5e-6,
        );
        assert!(rep.pass, "correct derivative must pass: {rep}");
        assert!(
            rep.max_rel_err < 1e-8,
            "central FD should be much tighter: {rep}"
        );
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
    fn non_finite_inputs_and_callback_evidence_fail_closed() {
        for non_finite in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_refused(&gradcheck(
                |_: [Dual64<1>; 1]| Dual64::<1>::constant(non_finite),
                [0.0],
                5e-6,
            ));
            assert_refused(&gradcheck(
                |_: [Dual64<1>; 1]| Dual64 {
                    re: 0.0,
                    eps: [non_finite],
                },
                [0.0],
                5e-6,
            ));
            assert_refused(&gradcheck(
                |[x]: [Dual64<1>; 1]| {
                    if x.eps[0] == 0.0 {
                        Dual64::<1>::constant(non_finite)
                    } else {
                        x
                    }
                },
                [0.0],
                5e-6,
            ));
            assert_refused(&gradcheck(
                |[x]: [Dual64<1>; 1]| {
                    if x.eps[0] == 0.0 && x.re.is_sign_negative() {
                        Dual64::<1>::constant(non_finite)
                    } else {
                        x
                    }
                },
                [0.0],
                5e-6,
            ));
            assert_refused(&gradcheck(|[x]: [Dual64<1>; 1]| x * x, [non_finite], 5e-6));
            assert_refused(&gradcheck(|[x]: [Dual64<1>; 1]| x * x, [0.0], non_finite));
        }
        assert_refused(&gradcheck(|[x]: [Dual64<1>; 1]| x * x, [0.0], -5e-6));

        let exact = gradcheck(
            |[x]: [Dual64<1>; 1]| {
                if x.eps[0] == 0.0 {
                    Dual64 {
                        re: x.re,
                        eps: [f64::NAN],
                    }
                } else {
                    x
                }
            },
            [0.0],
            0.0,
        );
        assert!(exact.pass, "zero tolerance retains exact-check behavior");
        assert_eq!(exact.max_rel_err.to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn vacuous_and_non_finite_intermediate_checks_fail_closed() {
        let calls = core::cell::Cell::new(0usize);
        assert_refused(&gradcheck::<0, _>(
            |_: [Dual64<0>; 0]| {
                calls.set(calls.get() + 1);
                Dual64::<0>::constant(0.0)
            },
            [],
            5e-6,
        ));
        assert_refused(&gradcheck(
            |[x]: [Dual64<1>; 1]| {
                calls.set(calls.get() + 1);
                x
            },
            [f64::NAN],
            5e-6,
        ));
        assert_refused(&gradcheck(
            |[x]: [Dual64<1>; 1]| {
                calls.set(calls.get() + 1);
                x
            },
            [0.0],
            f64::INFINITY,
        ));
        assert_eq!(calls.get(), 0, "preflight refusals must not call f");

        assert_refused(&gradcheck(|[x]: [Dual64<1>; 1]| x, [f64::MAX], 5e-6));
        assert_refused(&gradcheck(
            |[x]: [Dual64<1>; 1]| {
                Dual64::<1>::constant(if x.re.is_sign_negative() {
                    -f64::MAX
                } else {
                    f64::MAX
                })
            },
            [0.0],
            5e-6,
        ));
        assert_refused(&gradcheck(
            |[x]: [Dual64<1>; 1]| {
                Dual64::<1>::constant(if x.re.is_sign_negative() {
                    -5.0e303
                } else {
                    5.0e303
                })
            },
            [0.0],
            5e-6,
        ));
        assert_refused(&gradcheck(
            |[x]: [Dual64<1>; 1]| {
                if x.eps[0] != 0.0 {
                    Dual64 {
                        re: 0.0,
                        eps: [1.0e308],
                    }
                } else {
                    Dual64::<1>::constant(if x.re.is_sign_negative() {
                        6.0e302
                    } else {
                        -6.0e302
                    })
                }
            },
            [0.0],
            5e-6,
        ));
    }

    #[test]
    fn report_serializes_as_json_line() {
        let rep = gradcheck(|[x]: [Dual64<1>; 1]| x * x, [2.0], 5e-6);
        let line = format!("{rep}");
        assert!(
            line.starts_with('{') && line.ends_with('}'),
            "JSON-ish line: {line}"
        );
        assert!(line.contains("\"pass\":true"));
        assert_eq!(
            line,
            format!(
                "{{\"check\":\"gradcheck\",\"pass\":true,\"max_rel_err\":{:.3e},\"worst_index\":{},\"tol\":{:.1e}}}",
                rep.max_rel_err, rep.worst_index, rep.tol
            )
        );

        let refused = gradcheck(|[x]: [Dual64<1>; 1]| x * x, [0.0], f64::NAN);
        assert_eq!(
            format!("{refused}"),
            "{\"check\":\"gradcheck\",\"pass\":false,\"max_rel_err\":null,\"worst_index\":null,\"tol\":null}"
        );
    }
}
