//! Strong-Wolfe line search (bracket + zoom, Nocedal–Wright form)
//! with DETERMINISTIC control flow: fixed expansion factor, bisection
//! zoom, hard iteration caps — the same seed always walks the same
//! trajectory. Evaluation counts are returned so budgets stay honest.

/// Line-search outcome.
#[derive(Debug, Clone)]
pub struct WolfeOutcome {
    /// Accepted step length (0 on failure).
    pub alpha: f64,
    /// φ(α) at the accepted step.
    pub f_new: f64,
    /// Function+gradient evaluations spent.
    pub evals: usize,
    /// Whether both Wolfe conditions were certified.
    pub success: bool,
}

/// Strong Wolfe on φ(α) = f(x + α·d): sufficient decrease
/// φ(α) ≤ φ(0) + c₁·α·φ′(0) and curvature |φ′(α)| ≤ c₂·|φ′(0)|.
/// `phi` returns (value, derivative) at a step length.
pub fn strong_wolfe(
    phi: &mut dyn FnMut(f64) -> (f64, f64),
    f0: f64,
    dphi0: f64,
    alpha_init: f64,
    c1: f64,
    c2: f64,
) -> WolfeOutcome {
    strong_wolfe_with_budget(phi, f0, dphi0, alpha_init, c1, c2, usize::MAX)
}

/// Strong Wolfe with a hard upper bound on callback evaluations.
///
/// Unlike a post-hoc accounting check, this function checks the bound before
/// every call to `phi`, including calls made by zoom. A zero budget is a
/// well-formed unsuccessful search that performs no callback. Callers that
/// supplied a finite budget can distinguish exhaustion from an ordinary line
/// search failure by checking `!outcome.success && outcome.evals == max_evals`.
pub(crate) fn strong_wolfe_with_budget(
    phi: &mut dyn FnMut(f64) -> (f64, f64),
    f0: f64,
    dphi0: f64,
    alpha_init: f64,
    c1: f64,
    c2: f64,
    max_evals: usize,
) -> WolfeOutcome {
    assert!(
        dphi0.is_finite() && dphi0 < 0.0,
        "line search needs a descent direction (phi'(0) = {dphi0})"
    );
    assert!(f0.is_finite(), "line-search origin value must be finite");
    assert!(
        alpha_init.is_finite() && alpha_init > 0.0,
        "initial line-search step must be finite and positive"
    );
    assert!(
        c1.is_finite() && c2.is_finite() && 0.0 < c1 && c1 < c2 && c2 < 1.0,
        "strong-Wolfe constants must satisfy 0 < c1 < c2 < 1"
    );
    let mut evals = 0usize;
    let mut alpha_prev = 0.0f64;
    let mut f_prev = f0;
    let mut alpha = alpha_init;
    let max_expand = 20usize;
    for i in 0..max_expand {
        if evals == max_evals {
            return failed(f0, evals);
        }
        let (f_a, d_a) = phi(alpha);
        evals += 1;
        assert!(
            f_a.is_finite() && d_a.is_finite(),
            "line-search callback must return finite value and derivative"
        );
        if f_a > c1.mul_add(alpha * dphi0, f0) || (i > 0 && f_a >= f_prev) {
            return zoom(
                phi, f0, dphi0, alpha_prev, f_prev, alpha, c1, c2, evals, max_evals,
            );
        }
        if d_a.abs() <= c2 * dphi0.abs() {
            return WolfeOutcome {
                alpha,
                f_new: f_a,
                evals,
                success: true,
            };
        }
        if d_a >= 0.0 {
            return zoom(
                phi, f0, dphi0, alpha, f_a, alpha_prev, c1, c2, evals, max_evals,
            );
        }
        alpha_prev = alpha;
        f_prev = f_a;
        alpha *= 2.0;
        if !alpha.is_finite() {
            return failed(f0, evals);
        }
    }
    failed(f0, evals)
}

fn failed(f0: f64, evals: usize) -> WolfeOutcome {
    WolfeOutcome {
        alpha: 0.0,
        f_new: f0,
        evals,
        success: false,
    }
}

#[allow(clippy::too_many_arguments)]
fn zoom(
    phi: &mut dyn FnMut(f64) -> (f64, f64),
    f0: f64,
    dphi0: f64,
    mut lo: f64,
    mut f_lo: f64,
    mut hi: f64,
    c1: f64,
    c2: f64,
    mut evals: usize,
    max_evals: usize,
) -> WolfeOutcome {
    for _ in 0..40 {
        if evals == max_evals {
            return failed(f0, evals);
        }
        let alpha = f64::midpoint(lo, hi);
        let (f_a, d_a) = phi(alpha);
        evals += 1;
        assert!(
            f_a.is_finite() && d_a.is_finite(),
            "line-search callback must return finite value and derivative"
        );
        if f_a > c1.mul_add(alpha * dphi0, f0) || f_a >= f_lo {
            hi = alpha;
        } else {
            if d_a.abs() <= c2 * dphi0.abs() {
                return WolfeOutcome {
                    alpha,
                    f_new: f_a,
                    evals,
                    success: true,
                };
            }
            if d_a * (hi - lo) >= 0.0 {
                hi = lo;
            }
            lo = alpha;
            f_lo = f_a;
        }
        if (hi - lo).abs() < 1e-16 {
            break;
        }
    }
    failed(f0, evals)
}

#[cfg(test)]
mod tests {
    use super::strong_wolfe_with_budget;

    #[test]
    fn hard_budget_zero_never_calls_curve() {
        let mut calls = 0usize;
        let mut phi = |_alpha: f64| {
            calls += 1;
            (0.0, -1.0)
        };
        let outcome = strong_wolfe_with_budget(&mut phi, 1.0, -1.0, 1.0, 1e-4, 0.9, 0);
        assert!(!outcome.success);
        assert_eq!(outcome.evals, 0);
        assert_eq!(calls, 0);
    }

    #[test]
    fn hard_budget_covers_expansion_and_zoom() {
        let mut calls = 0usize;
        let mut phi = |_alpha: f64| {
            calls += 1;
            // No sufficient decrease, forcing zoom until the hard cap.
            (1.0, -1.0)
        };
        let outcome = strong_wolfe_with_budget(&mut phi, 1.0, -1.0, 1.0, 1e-4, 0.9, 3);
        assert!(!outcome.success);
        assert_eq!(outcome.evals, 3);
        assert_eq!(calls, 3);
    }

    #[test]
    fn expansion_overflow_fails_before_exposing_nonfinite_alpha() {
        let mut calls = 0usize;
        let mut phi = |alpha: f64| {
            assert!(alpha.is_finite());
            calls += 1;
            (0.0, -1.0)
        };
        let outcome = strong_wolfe_with_budget(&mut phi, 1.0, -1.0, f64::MAX, 1e-4, 0.9, 10);
        assert!(!outcome.success);
        assert_eq!(outcome.evals, 1);
        assert_eq!(calls, 1);
    }
}
