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
    /// Verdict at the supplied tolerance.
    pub pass: bool,
}

/// Verify a claimed gradient of `j` at `p` against central FD along
/// the supplied directions (callers pass deterministic keyed-stream
/// directions). `eps` is the FD step (scaled per direction by
/// ‖p‖/‖d‖ internally).
pub fn verify_gradient(
    j: &dyn Fn(&[f64]) -> f64,
    p: &[f64],
    gradient: &[f64],
    directions: &[Vec<f64>],
    eps: f64,
    tol: f64,
) -> GradientVerdict {
    assert_eq!(p.len(), gradient.len(), "gradient length mismatch");
    let mut pairs = Vec::with_capacity(directions.len());
    let mut worst = 0.0f64;
    for d in directions {
        assert_eq!(d.len(), p.len(), "direction length mismatch");
        let analytic: f64 = gradient.iter().zip(d).map(|(g, di)| g * di).sum();
        let mut plus = p.to_vec();
        let mut minus = p.to_vec();
        for i in 0..p.len() {
            plus[i] += eps * d[i];
            minus[i] -= eps * d[i];
        }
        let fd = (j(&plus) - j(&minus)) / (2.0 * eps);
        let scale = analytic.abs().max(fd.abs()).max(1e-12);
        worst = worst.max((analytic - fd).abs() / scale);
        pairs.push((analytic, fd));
    }
    GradientVerdict {
        max_rel_err: worst,
        pairs,
        pass: worst < tol,
    }
}
