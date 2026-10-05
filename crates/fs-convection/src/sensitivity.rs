//! Local Reynolds derivatives of admitted correlation evaluations.
//! The formula, group point and domain stay with their physical owner.

use super::{CorrelationError, CorrelationId, NusseltEvaluation, det};

impl NusseltEvaluation {
    /// `d ln(Nu) / d ln(Re)` at fixed Pr, length ratio, aspect and heat-flow
    /// direction. At fixed geometry/fluid properties this is also dln(h)/dln(Re).
    ///
    /// Supports the smooth duct, flat-plate and cylinder formulas. Returns
    /// None for the piecewise developing table, natural convection, or a
    /// Reynolds-dependent validity-axis endpoint. A caller must not interpret
    /// None as zero. Analytic laminar constants return Some(0); this does not
    /// make a convection-only derivative a total fan-speed derivative.
    /// Capacity and upstream air references also change with fan speed.
    ///
    /// This is a local Estimated formula derivative, not an interval or model
    /// validation certificate. No perturbed evaluator or finite-difference
    /// stencil supplies the value.
    ///
    /// # Errors
    /// Invalid groups or nonfinite derivative arithmetic.
    pub fn reynolds_elasticity(&self) -> Result<Option<f64>, CorrelationError> {
        if !self.evidence().model.in_domain { return Ok(None); }
        // These groups change proportionally to Re with the other inputs fixed.
        // Pr/aspect/L_over_Dh boundaries do not by themselves block this partial.
        for axis in ["Re", "Pe", "Gz"] {
            if let Some(&(low, high)) = self.card().model.validity.bounds().get(axis) {
                let value = self.groups().get(axis).copied().ok_or(CorrelationError::InvalidGroup {
                    axis, value_bits: f64::NAN.to_bits(),
                })?;
                if value <= low || value >= high { return Ok(None); }
            }
        }
        let group = |axis: &'static str| self.groups().get(axis).copied()
            .filter(|v| v.is_finite() && *v > 0.0)
            .ok_or(CorrelationError::InvalidGroup { axis, value_bits: f64::NAN.to_bits() });
        let nu = self.evidence().value;
        let value = match self.card().id {
            CorrelationId::CircularDuctLaminarCwt | CorrelationId::CircularDuctLaminarChf
            | CorrelationId::RectangularDuctLaminarCwt | CorrelationId::RectangularDuctLaminarChf => 0.0,
            CorrelationId::CircularDuctHausen => {
                let gz = group("Gz")?;
                let b = 0.04*det::pow(gz, 2.0/3.0);
                (0.0668*gz/(1.0+b))/nu*(1.0+b/3.0)/(1.0+b)
            }
            CorrelationId::DittusBoelter => 0.8,
            CorrelationId::Gnielinski => {
                let re = group("Re")?; let pr = group("Pr")?;
                let z = 0.79*det::ln(re)-1.64;
                let c = 12.7/(det::sqrt(8.0)*z)*(det::pow(pr, 2.0/3.0)-1.0);
                re/(re-1000.0) - 1.58/z*(1.0-0.5*c/(1.0+c))
            }
            CorrelationId::FlatPlateLaminarAverage => 0.5,
            CorrelationId::FlatPlateTurbulentAverage => {
                let term = 0.037*det::pow(group("Re")?, 0.8);
                0.8*term/(term-871.0)
            }
            CorrelationId::ChurchillBernsteinCylinder => {
                let b = det::pow(group("Re")?/282000.0, 5.0/8.0);
                (nu-0.3)/nu*0.5*(1.0+b/(1.0+b))
            }
            _ => return Ok(None),
        };
        if !value.is_finite() { return Err(CorrelationError::NonFiniteResult {
            stage: "Reynolds elasticity", value_bits: value.to_bits(),
        }); }
        Ok(Some(value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CorrelationInputs, ThermalDirection, evaluate};

    #[test]
    fn reynolds_elasticity_matches_the_actual_card_evaluator() {
        let step = 1e-5_f64;
        for (card,re,pr) in [
            (CorrelationId::CircularDuctLaminarCwt,1000.0,0.72),
            (CorrelationId::CircularDuctLaminarChf,1000.0,0.72),
            (CorrelationId::RectangularDuctLaminarCwt,1000.0,0.72),
            (CorrelationId::RectangularDuctLaminarChf,1000.0,0.72),
            (CorrelationId::CircularDuctHausen,100.0,0.72),
            (CorrelationId::CircularDuctHausen,1000.0,0.72),
            (CorrelationId::DittusBoelter,50_000.0,0.72),
            (CorrelationId::Gnielinski,20_000.0,0.72),
            (CorrelationId::Gnielinski,20_000.0,7.0),
            (CorrelationId::FlatPlateLaminarAverage,100_000.0,0.72),
            (CorrelationId::FlatPlateTurbulentAverage,1_000_000.0,0.72),
            (CorrelationId::ChurchillBernsteinCylinder,10_000.0,0.72),
        ] {
            for direction in [ThermalDirection::HeatingFluid,ThermalDirection::CoolingFluid] {
                let run = |r| evaluate(card,CorrelationInputs::forced(r,pr)
                    .with_length_ratio(1000.0).with_aspect_ratio(0.5).with_direction(direction)).unwrap();
                let actual = run(re).reynolds_elasticity().unwrap().unwrap();
                let expected = (run(re*step.exp()).evidence().value.ln()
                    - run(re*(-step).exp()).evidence().value.ln())/(2.0*step);
                assert!((actual-expected).abs() < 1e-8, "{card:?}: {actual} != {expected}");
            }
        }
    }
    #[test]
    fn reynolds_elasticity_does_not_invent_table_or_endpoint_derivatives() {
        let run = |r| evaluate(CorrelationId::Gnielinski,CorrelationInputs::forced(r,0.72)
            .with_length_ratio(1000.0)).unwrap();
        let nominal = run(20_000.0);
        let (low,high) = nominal.card().model.validity.bounds()["Re"];
        for endpoint in [low,high] {
            assert_eq!(run(endpoint).reynolds_elasticity().unwrap(),None);
        }
        let table = evaluate(CorrelationId::RectangularDuctLaminarCwtDevelopingPr072,
            CorrelationInputs::forced(1000.0,0.72).with_length_ratio(24.0)
                .with_aspect_ratio(0.5)).unwrap();
        assert_eq!(table.reynolds_elasticity().unwrap(),None);
        let natural = evaluate(CorrelationId::ChurchillChuVerticalPlate,
            CorrelationInputs::natural(1e6,0.72)).unwrap();
        assert_eq!(natural.reynolds_elasticity().unwrap(),None);
    }
}
