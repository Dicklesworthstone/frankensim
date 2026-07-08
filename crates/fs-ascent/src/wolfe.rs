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
    assert!(
        dphi0 < 0.0,
        "line search needs a descent direction (phi'(0) = {dphi0})"
    );
    let mut evals = 0usize;
    let mut alpha_prev = 0.0f64;
    let mut f_prev = f0;
    let mut alpha = alpha_init;
    let max_expand = 20usize;
    for i in 0..max_expand {
        let (f_a, d_a) = phi(alpha);
        evals += 1;
        if f_a > c1.mul_add(alpha * dphi0, f0) || (i > 0 && f_a >= f_prev) {
            return zoom(phi, f0, dphi0, alpha_prev, f_prev, alpha, c1, c2, evals);
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
            return zoom(phi, f0, dphi0, alpha, f_a, alpha_prev, c1, c2, evals);
        }
        alpha_prev = alpha;
        f_prev = f_a;
        alpha *= 2.0;
    }
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
) -> WolfeOutcome {
    for _ in 0..40 {
        let alpha = f64::midpoint(lo, hi);
        let (f_a, d_a) = phi(alpha);
        evals += 1;
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
    WolfeOutcome {
        alpha: 0.0,
        f_new: f0,
        evals,
        success: false,
    }
}
