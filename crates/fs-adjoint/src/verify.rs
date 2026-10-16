//! The gradient-verification gate: every adjoint gradient is checked
//! against central finite differences along random directions before
//! it is trusted (AGENTS.md: "gradients must be verified against
//! independent checks"). ci-gauntlet wires this so a solver without a
//! passing gradient check cannot merge; the helper returns a verdict
//! object with the worst direction, not just a boolean.

/// Outcome of a gradient verification.
#[derive(Debug, Clone)]
pub struct GradientVerdict {
    /// Worst relative error across the probed directions.
    pub max_rel_err: f64,
    /// Directional derivatives (analytic, finite-difference) pairs.
    pub pairs: Vec<(f64, f64)>,
    /// How many probed directions carried SIGNAL — directions where the
    /// larger of `|analytic|` and `|fd|` cleared the `1e-12·‖d‖∞` floor.
    ///
    /// A direction whose analytic and finite-difference derivatives are
    /// BOTH at that floor (both exactly `0.0` is the common case) has
    /// zero discriminating power: the comparison `|analytic − fd| /
    /// scale` is `0/scale = 0` — a perfect score — whether the gradient
    /// is right or silently zero. Such a probe is not weak evidence, it
    /// is no evidence, so it cannot contribute to a passing verdict; a
    /// verification with no informative direction fails closed with the
    /// `+∞` sentinel and `informative_directions == 0`.
    pub informative_directions: usize,
    /// Verdict at the supplied tolerance.
    pub pass: bool,
}

fn fail_closed(pairs: Vec<(f64, f64)>) -> GradientVerdict {
    GradientVerdict {
        max_rel_err: f64::INFINITY,
        pairs,
        informative_directions: 0,
        pass: false,
    }
}

/// Verify a claimed gradient of `j` at `p` against central FD along
/// the supplied directions (callers pass deterministic keyed-stream
/// directions). `eps` is the literal scalar step in `p ± eps * d`;
/// the perturbation is not rescaled by the point or direction norm.
/// The relative-error floor is `1e-12 * ‖d‖∞`, making that floor
/// homogeneous under paired rescalings of `d` and the inverse `eps`.
///
/// The same floor decides whether a direction is INFORMATIVE. A probe
/// whose analytic and finite-difference directional derivatives are
/// both at or below it (the silent-zero seam: a broken forward path
/// yields a zero cotangent AND two bit-identical objective values)
/// scores a relative error of `0.0` no matter what the gradient is, so
/// it is refused as vacuous evidence rather than counted as a pass —
/// see [`GradientVerdict::informative_directions`].
pub fn verify_gradient(
    j: &dyn Fn(&[f64]) -> f64,
    p: &[f64],
    gradient: &[f64],
    directions: &[Vec<f64>],
    eps: f64,
    tol: f64,
) -> GradientVerdict {
    assert_eq!(p.len(), gradient.len(), "gradient length mismatch");
    for d in directions {
        assert_eq!(d.len(), p.len(), "direction length mismatch");
    }
    if !eps.is_finite()
        || eps <= 0.0
        || !tol.is_finite()
        || tol <= 0.0
        || !p.iter().all(|value| value.is_finite())
        || !gradient.iter().all(|value| value.is_finite())
        || directions.is_empty()
    {
        return fail_closed(Vec::new());
    }

    let mut pairs = Vec::with_capacity(directions.len());
    let mut worst = 0.0f64;
    let mut informative_directions = 0usize;
    for d in directions {
        if !d.iter().all(|value| value.is_finite()) || d.iter().all(|value| *value == 0.0) {
            return fail_closed(pairs);
        }
        let direction_magnitude = d.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
        if !direction_magnitude.is_finite() || direction_magnitude == 0.0 {
            return fail_closed(pairs);
        }
        let mut analytic = 0.0;
        for (&gradient_value, &direction_value) in gradient.iter().zip(d) {
            let term = gradient_value * direction_value;
            if !term.is_finite() {
                return fail_closed(pairs);
            }
            analytic += term;
            if !analytic.is_finite() {
                return fail_closed(pairs);
            }
        }
        let mut plus = p.to_vec();
        let mut minus = p.to_vec();
        for i in 0..p.len() {
            let step = eps * d[i];
            if !step.is_finite() {
                return fail_closed(pairs);
            }
            let plus_value = p[i] + step;
            let minus_value = p[i] - step;
            if !plus_value.is_finite() || !minus_value.is_finite() {
                return fail_closed(pairs);
            }
            if d[i] != 0.0
                && (plus_value.to_bits() == p[i].to_bits()
                    || minus_value.to_bits() == p[i].to_bits())
            {
                return fail_closed(pairs);
            }
            plus[i] = plus_value;
            minus[i] = minus_value;
        }
        let j_plus = j(&plus);
        let j_minus = j(&minus);
        if !j_plus.is_finite() || !j_minus.is_finite() {
            return fail_closed(pairs);
        }
        let fd_numerator = j_plus - j_minus;
        let fd_denominator = 2.0 * eps;
        if !fd_numerator.is_finite() || !fd_denominator.is_finite() {
            return fail_closed(pairs);
        }
        let fd = fd_numerator / fd_denominator;
        if !fd.is_finite() {
            return fail_closed(pairs);
        }
        let scale_floor = 1e-12 * direction_magnitude;
        let signal = analytic.abs().max(fd.abs());
        if signal > scale_floor {
            informative_directions += 1;
        }
        let scale = signal.max(scale_floor);
        let error_numerator = (analytic - fd).abs();
        if !scale.is_finite() || !error_numerator.is_finite() {
            return fail_closed(pairs);
        }
        let relative_error = error_numerator / scale;
        if !relative_error.is_finite() {
            return fail_closed(pairs);
        }
        worst = worst.max(relative_error);
        pairs.push((analytic, fd));
    }
    if informative_directions == 0 {
        // Every probe was bit-insensitive: the objective did not move and
        // the claimed gradient claims no movement either. The arithmetic
        // would report `max_rel_err = 0.0` — the PERFECT score — from an
        // experiment that cannot tell a correct gradient from a silently
        // zero one. Refuse with the same +∞ sentinel used for corrupted
        // evidence: a gate that cannot fail is not a gate.
        return fail_closed(pairs);
    }
    let pass = !pairs.is_empty()
        && pairs.len() == directions.len()
        && informative_directions > 0
        && worst.is_finite()
        && worst < tol;
    GradientVerdict {
        max_rel_err: worst,
        pairs,
        informative_directions,
        pass,
    }
}
