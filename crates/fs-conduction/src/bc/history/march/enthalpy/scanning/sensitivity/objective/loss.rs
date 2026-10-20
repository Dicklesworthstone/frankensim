//! Fixed residual losses for physical scan measurements.
//!
//! With r=(prediction-target)/scale, Huber is r²/2 for |r|<=delta and
//! delta*(|r|-delta/2) otherwise. Its derivative is r or signed delta.
//! This is the Huber (1964) / scipy.special.huber convention, also used by
//! fs-ascent's transient sensor adapter. There is deliberately no dependency
//! from the conduction layer back to that higher-level calibration adapter.
//!
//! Only the scalar fitting objective changes. The original measured field,
//! spatial transpose, phase chart, physical solver and hard constraints do not.
//! Raw residuals remain visible, including large ones. A robust fit is not an
//! outlier classification, uncertainty estimate or physical-validity certificate.

use super::{ScanObservable, ScanSensitivityError, ScanStep, finite, invalid};
use super::spatial::{ScanObservationSpec, sealed};
use crate::ConductionMesh;
use fs_exec::Cx;

/// Dimensionless loss after dividing a physical measurement residual by scale.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ScanResidualLoss {
    /// The existing least-squares convention, r²/2.
    #[default]
    Quadratic,
    /// Explicit positive threshold in standardized-residual units, not kelvin,
    /// J/kg, percent, a source-power limit or an inferred noise distribution.
    Huber { threshold: f64 },
}

impl ScanResidualLoss {
    /// Refuse a nonpositive or nonfinite Huber threshold, even at zero residual.
    pub fn validate(self) -> Result<(), ScanSensitivityError> {
        if let Self::Huber { threshold } = self {
            if !threshold.is_finite() || threshold <= 0.0 {
                return Err(invalid("scan Huber threshold must be finite and strictly positive"));
            }
        }
        Ok(())
    }

    /// Return (rho(r), d rho/dr), without altering or classifying r.
    /// Divide the derivative by the original physical residual scale exactly
    /// once before scattering it through the observation's spatial transpose.
    /// Huber is C1 at its threshold; no second derivative is claimed there.
    ///
    /// # Errors
    /// Invalid policy, nonfinite residual, or unrepresentable loss/derivative.
    /// The linear branch never squares a large residual as an intermediate.
    pub fn value_derivative(self, residual: f64) -> Result<(f64, f64), ScanSensitivityError> {
        self.validate()?;
        finite(residual)?;
        let (value, derivative) = match self {
            Self::Huber { threshold } if residual.abs() > threshold =>
                (threshold * (residual.abs() - 0.5 * threshold), threshold.copysign(residual)),
            _ => ((0.5 * residual) * residual, residual),
        };
        Ok((finite(value)?, finite(derivative)?))
    }
}

/// Attach an immutable fitting loss to either supported observation layout.
/// Wrapping preserves its physical units, target, scale, mesh identity and
/// exact spatial row. The outer policy replaces an inner policy; it does not
/// compose losses or reinterpret a spatial mean as separate node residuals.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanObservationWithLoss<O> {
    observation: O,
    loss: ScanResidualLoss,
}

impl<O: ScanObservationSpec> ScanObservationWithLoss<O> {
    /// Admit the fixed loss. Original endpoint/readout/target/scale admission
    /// remains with the shared objective, before any physical solve.
    pub fn new(observation: O, loss: ScanResidualLoss) -> Result<Self, ScanSensitivityError> {
        loss.validate()?;
        Ok(Self { observation, loss })
    }

    /// The original immutable physical measurement declaration.
    pub fn observation(&self) -> &O { &self.observation }

    /// The fixed policy; changing it requires constructing another observation.
    pub const fn residual_loss(&self) -> ScanResidualLoss { self.loss }
}

impl<O: ScanObservationSpec> sealed::Sealed for ScanObservationWithLoss<O> {}

impl<O: ScanObservationSpec> ScanObservationSpec for ScanObservationWithLoss<O> {
    fn endpoint(&self) -> usize { self.observation.endpoint() }
    fn observable(&self) -> ScanObservable { self.observation.observable() }
    fn target(&self) -> f64 { self.observation.target() }
    fn scale(&self) -> f64 { self.observation.scale() }
    fn loss(&self) -> ScanResidualLoss { self.loss }
    fn validate_readout(&self, mesh: &ConductionMesh) -> Result<(), ScanSensitivityError> {
        self.observation.validate_readout(mesh)
    }
    fn prediction(&self, cx: &Cx<'_>, report: &ScanStep) -> Result<f64, ScanSensitivityError> {
        self.observation.prediction(cx, report)
    }
    fn add_seed(&self, cx: &Cx<'_>, seed: &mut [f64], value: f64) -> Result<(), ScanSensitivityError> {
        self.observation.add_seed(cx, seed, value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_huber_retains_quadratic_core_linear_tails_and_continuous_gradient() {
        let loss = ScanResidualLoss::Huber { threshold: 2.0 };
        for (r, value, derivative) in [(-4.0, 6.0, -2.0), (-2.0, 2.0, -2.0),
            (-1.0, 0.5, -1.0), (0.0, 0.0, 0.0), (1.0, 0.5, 1.0),
            (2.0, 2.0, 2.0), (4.0, 6.0, 2.0)]
        {
            assert_eq!(loss.value_derivative(r).unwrap(), (value, derivative));
            let epsilon = 1e-5;
            let fd = (loss.value_derivative(r + epsilon).unwrap().0
                - loss.value_derivative(r - epsilon).unwrap().0) / (2.0 * epsilon);
            assert!((fd - derivative).abs() < 3e-6);
        }
        assert!(loss.value_derivative(1e300).unwrap().0.is_finite());
        assert!(ScanResidualLoss::Quadratic.value_derivative(1e300).is_err());
    }

    #[test]
    fn g0_invalid_policy_or_data_is_not_hidden_by_robust_scoring() {
        for threshold in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(ScanResidualLoss::Huber { threshold }.value_derivative(0.0).is_err());
        }
        for loss in [ScanResidualLoss::Quadratic, ScanResidualLoss::Huber { threshold: 1.0 }] {
            for r in [f64::NAN, f64::NEG_INFINITY, f64::INFINITY] {
                assert!(loss.value_derivative(r).is_err());
            }
        }
        assert!(ScanResidualLoss::Huber { threshold: f64::MAX }
            .value_derivative(f64::MAX).is_err());
    }

    #[test]
    fn g3_wrapper_keeps_physical_measurement_and_original_quadratic_default() {
        let row = super::super::ScanObservation { endpoint: 2, vertex: 3,
            observable: ScanObservable::SpecificEnthalpy, target: 157.0, scale: 10.0 };
        let wrapped = ScanObservationWithLoss::new(row, ScanResidualLoss::Huber { threshold: 1.5 }).unwrap();
        assert_eq!(wrapped.observation(), &row);
        assert_eq!(wrapped.endpoint(), row.endpoint);
        assert_eq!(wrapped.observable(), row.observable);
        assert_eq!(wrapped.target(), row.target);
        assert_eq!(wrapped.scale(), row.scale);
        assert_eq!(row.loss(), ScanResidualLoss::Quadratic);
        assert_eq!(wrapped.loss().value_derivative(-3.0).unwrap(), (3.375, -1.5));
        for r in [-3.0, 0.0, 0.25, 10.0] {
            assert_eq!(row.loss().value_derivative(r).unwrap(), ((0.5*r)*r, r));
        }
        assert!(ScanObservationWithLoss::new(row, ScanResidualLoss::Huber { threshold: 0.0 }).is_err());
    }
}
