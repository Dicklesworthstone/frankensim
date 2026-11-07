//! Residual-based globalization of the native conductivity Picard update.
//!
//! The caller evaluates the original physical operator on the segment from
//! its current iterate to the full proposal. This helper selects only a step;
//! it changes no material law, tolerance, history or linear-solve budget.

/// At most forty original-operator evaluations, with no linear solves.
/// A small step never constitutes convergence: accept only the requested
/// residual tolerance or sufficient decrease from the current residual.
pub(super) fn backtrack<T, E>(
    current_residual: f64,
    tolerance: f64,
    mut evaluate: impl FnMut(f64) -> Result<(T, f64), E>,
) -> Result<Option<(T, f64)>, E> {
    let mut fraction = 1.0;
    for _ in 0..40 {
        let (candidate, residual) = evaluate(fraction)?;
        if residual.is_finite()
            && (residual <= tolerance
                || residual < (1.0 - 1e-4 * fraction) * current_residual)
        {
            return Ok(Some((candidate, residual)));
        }
        fraction *= 0.5;
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::backtrack;

    // Independent exact Fourier integral: k=1 to 350 K, a linear rise to
    // k=100 at 351 K, then k=100. The boundary reservoir is fixed at 300 K.
    fn potential(t: f64) -> f64 {
        let ramp = (t - 350.0).clamp(0.0, 1.0);
        t.min(350.0) - 300.0 + ramp + 49.5 * ramp * ramp
            + 100.0 * (t - 351.0).max(0.0)
    }

    #[test]
    fn sharp_conductivity_converges_using_actual_fourier_residual() {
        let power = 2.0 * potential(350.5);
        let merit = |t| {
            let heat = 2.0 * potential(t);
            (heat - power).abs() / (heat.abs() + power)
        };
        let proposal = |t| 300.0 + power * (t - 300.0) / (2.0 * potential(t));
        let mut t = 300.0 + power / 2.0;
        let initial = t;
        // The undamped production Picard equation has an exact two-cycle.
        assert!((proposal(proposal(t)) - t).abs() < 1e-12);
        let mut residual = merit(t);
        let mut damped = false;
        for _ in 0..200 {
            if residual <= 1e-12 { break; }
            let full = proposal(t);
            let (next, next_residual) = backtrack(residual, 1e-12, |fraction| {
                damped |= fraction < 1.0;
                let candidate = fraction.mul_add(full - t, t);
                Ok::<_, ()>((candidate, merit(candidate)))
            }).unwrap().expect("a descending physical step");
            assert!(next_residual < residual);
            (t, residual) = (next, next_residual);
        }
        assert!(damped, "the undamped cycle starts at {initial}");
        assert!(residual <= 1e-12, "T={t}, residual={residual}");
        assert!((t - 350.5).abs() < 1e-10);
        assert!((2.0 * potential(t) - power).abs() < 1e-9);
    }

    #[test]
    fn no_descent_is_bounded_and_evaluation_errors_propagate() {
        let mut calls = 0;
        let refused = backtrack(0.5, 1e-12, |_| {
            calls += 1;
            Ok::<_, ()>(((), 0.5))
        }).unwrap();
        assert!(refused.is_none());
        assert_eq!(calls, 40);
        let error = backtrack::<(), _>(0.5, 1e-12, |_| Err("cancelled"));
        assert_eq!(error, Err("cancelled"));
    }
}
