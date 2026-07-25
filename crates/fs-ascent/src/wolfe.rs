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
            admissible_probe(f_a, d_a),
            "line-search callback must return a finite or +inf value and a finite derivative"
        );
        // `f_a == +inf` lands here and is necessarily greater than the finite
        // sufficient-decrease threshold, so an infeasible probe is bracketed
        // into [alpha_prev, alpha] rather than accepted. `f_prev` therefore
        // never becomes infinite: this branch returns before the `f_prev = f_a`
        // update below, which keeps zoom's `lo` endpoint finite.
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

/// Whether a probe `(φ(α), φ′(α))` is an admissible line-search observation.
///
/// `φ(α) = +∞` IS admissible and means "this trial step is outside the
/// objective's domain, so it is too long". That is the documented barrier
/// convention: [`crate::interior::interior_point`] returns `+∞` with a finite
/// all-zero gradient for an infeasible sample precisely so the search shrinks
/// back into the interior. Both probe sites below then reject `+∞` on the
/// sufficient-decrease test and bisect toward the feasible endpoint, which is
/// the behaviour that comment already promises.
///
/// `NaN` and `-∞` are NOT admissible and still panic. Neither carries
/// "too long" information: `NaN` poisons every ordering comparison the bracket
/// and zoom logic depend on, and `-∞` asserts an objective unbounded below,
/// which is a modelling error rather than a step-length one. A non-finite
/// derivative is likewise refused — a broken callback must fail loudly instead
/// of steering the search.
fn admissible_probe(f_a: f64, d_a: f64) -> bool {
    (f_a.is_finite() || f_a == f64::INFINITY) && d_a.is_finite()
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
            admissible_probe(f_a, d_a),
            "line-search callback must return a finite or +inf value and a finite derivative"
        );
        // An infeasible (`+inf`) midpoint takes this branch, moving `hi` down to
        // it and bisecting back toward the feasible `lo`. `f_lo` is only ever
        // assigned in the `else` arm, so the interval keeps one finite endpoint.
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

    /// A log-barrier-shaped curve: a smooth quadratic on the feasible side and
    /// `+inf` once the trial step leaves the domain at `alpha >= 0.5`. This is
    /// the exact shape `interior::interior_point` presents (bead frankensim-35yeu).
    ///
    /// `phi(a) = 0.5*(a - 0.25)^2`, so `phi(0) = 0.03125` and `phi'(0) = -0.25`.
    fn barrier_curve(alpha: f64) -> (f64, f64) {
        if alpha >= 0.5 {
            return (f64::INFINITY, 0.0);
        }
        let shifted = alpha - 0.25;
        (0.5 * shifted * shifted, shifted)
    }

    #[test]
    fn infeasible_probe_backtracks_into_the_interior_instead_of_panicking() {
        let mut calls = Vec::new();
        let mut phi = |alpha: f64| {
            calls.push(alpha);
            barrier_curve(alpha)
        };
        // The initial step is deliberately infeasible.
        let outcome =
            strong_wolfe_with_budget(&mut phi, 0.03125, -0.25, 1.0, 1e-4, 0.9, usize::MAX);

        // Not merely "did not panic": the search must certify a step strictly
        // inside the feasible region, at the curve's interior minimiser.
        assert!(
            outcome.success,
            "infeasible bracket must still certify a step"
        );
        assert!(
            outcome.alpha > 0.0 && outcome.alpha < 0.5,
            "accepted step must be strictly feasible, got {}",
            outcome.alpha
        );
        assert_eq!(outcome.alpha, 0.25);
        assert_eq!(outcome.f_new, 0.0);
        // alpha=1.0 (+inf) brackets, zoom bisects to 0.5 (+inf), then 0.25.
        assert_eq!(calls, vec![1.0, 0.5, 0.25]);
        assert_eq!(outcome.evals, 3);
    }

    #[test]
    #[should_panic(expected = "finite or +inf value and a finite derivative")]
    fn nan_probe_value_still_refuses() {
        let mut phi = |_alpha: f64| (f64::NAN, -1.0);
        let _ = strong_wolfe_with_budget(&mut phi, 1.0, -1.0, 1.0, 1e-4, 0.9, usize::MAX);
    }

    #[test]
    #[should_panic(expected = "finite or +inf value and a finite derivative")]
    fn negative_infinite_probe_value_still_refuses() {
        let mut phi = |_alpha: f64| (f64::NEG_INFINITY, -1.0);
        let _ = strong_wolfe_with_budget(&mut phi, 1.0, -1.0, 1.0, 1e-4, 0.9, usize::MAX);
    }

    #[test]
    #[should_panic(expected = "finite or +inf value and a finite derivative")]
    fn nonfinite_probe_derivative_still_refuses() {
        let mut phi = |_alpha: f64| (0.0, f64::NAN);
        let _ = strong_wolfe_with_budget(&mut phi, 1.0, -1.0, 1.0, 1e-4, 0.9, usize::MAX);
    }

    #[test]
    #[should_panic(expected = "finite or +inf value and a finite derivative")]
    fn nan_probe_inside_zoom_still_refuses() {
        // Feasible-looking bracket, then a poisoned midpoint once zoom starts.
        let mut probes = 0usize;
        let mut phi = |_alpha: f64| {
            probes += 1;
            if probes == 1 {
                (5.0, -1.0)
            } else {
                (f64::NAN, -1.0)
            }
        };
        let _ = strong_wolfe_with_budget(&mut phi, 1.0, -1.0, 1.0, 1e-4, 0.9, usize::MAX);
    }

    #[test]
    fn expansion_overflow_fails_before_exposing_nonfinite_alpha() {
        let mut calls = 0usize;
        // The probe value must be low enough to PASS sufficient decrease at
        // `alpha = f64::MAX`, or the bracket zooms instead of expanding and this
        // test never reaches the overflow path it is named for. The threshold
        // there is `c1 * (f64::MAX * -1.0) + 1.0`, about -1.8e304, so a modest
        // value like `0.0` fails the test's own premise: it exceeds the
        // threshold, triggers zoom, and burns the whole eval budget on
        // bisection. `-f64::MAX` clears it, while `d_a = -1.0` keeps curvature
        // unsatisfied and the direction descending, so the loop expands.
        let mut phi = |alpha: f64| {
            assert!(alpha.is_finite(), "phi must never see a non-finite step");
            calls += 1;
            (-f64::MAX, -1.0)
        };
        let outcome = strong_wolfe_with_budget(&mut phi, 1.0, -1.0, f64::MAX, 1e-4, 0.9, 10);
        assert!(!outcome.success);
        // One probe, then `alpha *= 2.0` overflows to +inf and the search bails
        // out rather than handing a non-finite step to the callback.
        assert_eq!(outcome.evals, 1);
        assert_eq!(calls, 1);
    }
}
