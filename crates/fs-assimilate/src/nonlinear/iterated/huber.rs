//! Huber loss for explicitly independent, standardized scalar residuals.
//!
//! Huber, "Robust Estimation of a Location Parameter" (1964),
//! DOI: 10.1214/aoms/1177703732, defines the quadratic/linear loss. In the
//! convention used here, `rho(r) = r²/2` for `|r| <= delta` and otherwise
//! `delta * (|r| - delta/2)`. Differentiation gives the bounded score
//! `psi(r) = clamp(r, -delta, delta)`. An IRLS proposal uses
//! `w(r) = psi(r)/r = min(1, delta/|r|)`, with `w(0) = 1`.
//!
//! The driver must evaluate the actual loss and score for descent and
//! stationarity. The weighted quadratic is only a local proposal surrogate;
//! its inverse curvature is not a calibrated posterior covariance.

use crate::AssimError;

use super::{finite, invalid};

/// A caller-declared Huber transition in standardized residual units.
///
/// Residuals up to this magnitude retain the quadratic loss. Beyond it, the
/// loss grows linearly and the score is bounded by the transition. This does
/// not estimate the noise scale or establish a contamination probability.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HuberLoss {
    transition: f64,
}

impl HuberLoss {
    /// Declare a strictly positive, finite transition without choosing a
    /// statistical tuning rule on the caller's behalf.
    pub fn new(transition: f64) -> Result<Self, AssimError> {
        if !transition.is_finite() || transition <= 0.0 {
            return Err(invalid("Huber transition"));
        }
        Ok(Self { transition })
    }

    /// Transition from quadratic to linear loss in standardized residual units.
    #[must_use]
    pub const fn transition(self) -> f64 {
        self.transition
    }

    /// Actual loss, bounded score, and square root of the IRLS weight.
    pub(super) fn terms(self, residual: f64) -> Result<(f64, f64, f64), AssimError> {
        let residual = finite(residual, "Huber standardized residual")?;
        let magnitude = residual.abs();
        let delta = self.transition;
        if magnitude <= delta {
            // Halve before multiplying: r² may overflow while r²/2 is finite.
            let objective = finite((0.5 * magnitude) * magnitude, "Huber objective")?;
            Ok((objective, residual, 1.0))
        } else {
            // Subtract before multiplying: delta*|r| may overflow while this
            // complete tail loss remains finite.
            let objective = finite(delta * (magnitude - 0.5 * delta), "Huber objective")?;
            let influence = delta.copysign(residual);
            // Never form delta/|r| first: it can underflow although its square
            // root, and the weighted physical Jacobian, remain representable.
            let sqrt_weight = finite(delta.sqrt() / magnitude.sqrt(), "Huber root weight")?;
            Ok((objective, influence, sqrt_weight))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::HuberLoss;

    #[test]
    fn piecewise_loss_and_score_join_at_both_transitions() {
        let loss = HuberLoss::new(1.0).unwrap();
        for (magnitude, expected_loss, expected_score, expected_root) in [
            (0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64),
            (0.5, 0.125, 0.5, 1.0),
            (1.0, 0.5, 1.0, 1.0),
            (4.0, 3.5, 1.0, 0.5),
        ] {
            for sign in [-1.0, 1.0] {
                let (objective, influence, root) = loss.terms(sign * magnitude).unwrap();
                assert_eq!(objective, expected_loss);
                assert_eq!(influence, sign * expected_score);
                assert_eq!(root, expected_root);
            }
        }
        assert_eq!(loss.terms(-0.0).unwrap().1.to_bits(), 0);
    }

    #[test]
    fn extreme_terms_keep_finite_objectives_and_nonzero_root_weights() {
        // The weight is 2^-1600 and rounds to zero; its root is exactly 2^-800.
        let delta = 2.0_f64.powi(-1000);
        let magnitude = 2.0_f64.powi(600);
        assert_eq!(delta / magnitude, 0.0);
        for sign in [-1.0, 1.0] {
            let (objective, influence, root) = HuberLoss::new(delta)
                .unwrap()
                .terms(sign * magnitude)
                .unwrap();
            assert_eq!(objective, 2.0_f64.powi(-400));
            assert_eq!(influence, sign * delta);
            assert_eq!(root, 2.0_f64.powi(-800));
        }

        // Both unsafely ordered products overflow, but the complete losses fit.
        let delta = 2.0_f64.powi(512);
        let loss = HuberLoss::new(delta).unwrap();
        assert!((delta * delta).is_infinite());
        assert_eq!(loss.terms(delta).unwrap().0, 2.0_f64.powi(1023));
        assert_eq!(
            loss.terms(1.25 * delta).unwrap().0,
            1.5 * 2.0_f64.powi(1023)
        );

        let smallest = HuberLoss::new(f64::from_bits(1)).unwrap();
        let (objective, influence, root) = smallest.terms(f64::MAX).unwrap();
        assert!(objective.is_finite() && objective > 0.0);
        assert_eq!(influence, f64::from_bits(1));
        assert!(root.is_finite() && root > 0.0);
    }

    #[test]
    fn invalid_transitions_and_nonfinite_terms_refuse() {
        for transition in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(HuberLoss::new(transition).is_err());
        }
        let loss = HuberLoss::new(1.5).unwrap();
        assert_eq!(loss.transition(), 1.5);
        for residual in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(loss.terms(residual).is_err());
        }
        assert!(HuberLoss::new(f64::MAX).unwrap().terms(f64::MAX).is_err());
    }
}
