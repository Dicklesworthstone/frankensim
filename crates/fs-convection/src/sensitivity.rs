//! Local input derivatives of admitted correlation evaluations.
//! The formula, group point and domain stay with their physical owner.

use super::{
    CorrelationError, CorrelationId, CorrelationInputs, NusseltEvaluation,
    ProvenanceHash, ThermalDirection, det, evaluate,
};

const INPUT_AXES: [&str; 5] = ["Re", "Pr", "L_over_Dh", "aspect_ratio", "Ra"];

impl CorrelationInputs {
    /// Evaluate a card and attach its analytic input partials to the returned
    /// evidence's `sensitivity.d_qoi` map.
    ///
    /// Entries are **absolute** derivatives `dNu/dx` for the supplied
    /// independent inputs `Re`, `Pr`, `L_over_Dh`, `aspect_ratio`, and `Ra`.
    /// Pe and Gz are dependent groups: their chain rules are included, never
    /// exposed as extra independent coordinates. The heating/cooling branch
    /// stays fixed, including at Pr = 1 where its two values coincide but its
    /// derivatives do not.
    ///
    /// A missing entry means unavailable, NOT zero. Explicit zeros describe
    /// admitted, locally constant directions. Partials that move a validity
    /// axis at its endpoint are omitted; a fixed, unrelated endpoint does not
    /// invalidate the other partials. The developing rectangular table does
    /// not advertise smooth input derivatives, including inside its bridge.
    ///
    /// These are local formula estimates, not error enclosures or physical
    /// validation. The nominal card, source region, and numerical status are
    /// unchanged. Use `evaluate` when only the primal value is required.
    ///
    /// # Errors
    /// The ordinary evaluator's admission errors, or nonfinite derivative
    /// arithmetic. No finite-difference evaluations supply a derivative.
    pub fn evaluate_with_sensitivities(
        self,
        id: CorrelationId,
    ) -> Result<NusseltEvaluation, CorrelationError> {
        let mut evaluation = evaluate(id, self)?;
        for axis in INPUT_AXES {
            let Some(&input) = evaluation.groups().get(axis) else {
                continue;
            };
            if let Some(elasticity) = input_elasticity(&evaluation, self.direction, axis)? {
                let derivative = finite_derivative(
                    evaluation.evidence().value * elasticity / input,
                    "Nusselt input derivative",
                )?;
                evaluation
                    .evidence
                    .sensitivity
                    .d_qoi
                    .insert(axis.to_string(), derivative);
            }
        }
        evaluation.evidence.provenance = ProvenanceHash::chain(
            "convection-input-sensitivities-v1",
            &[evaluation.evidence.provenance],
        );
        Ok(evaluation)
    }
}

fn finite_derivative(value: f64, stage: &'static str) -> Result<f64, CorrelationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CorrelationError::NonFiniteResult {
            stage,
            value_bits: value.to_bits(),
        })
    }
}

fn input_is_interior(evaluation: &NusseltEvaluation, axis: &'static str) -> bool {
    let moving_axes: &[&str] = match axis {
        "Re" => &["Re", "Pe", "Gz"],
        "Pr" => &["Pr", "Pe", "Gz"],
        "L_over_Dh" => &["L_over_Dh", "Gz"],
        "aspect_ratio" => &["aspect_ratio"],
        "Ra" => &["Ra"],
        _ => return false,
    };
    evaluation.evidence().model.in_domain
        && moving_axes.iter().all(|&moving| {
            evaluation.card().model.validity.bounds().get(moving).is_none_or(|&(low, high)| {
                evaluation.groups().get(moving).is_some_and(|&value| value > low && value < high)
            })
        })
}

fn input_elasticity(
    evaluation: &NusseltEvaluation,
    direction: ThermalDirection,
    axis: &'static str,
) -> Result<Option<f64>, CorrelationError> {
    if !input_is_interior(evaluation, axis) {
        return Ok(None);
    }
    // Retain the public Reynolds helper's support and refusal semantics.
    if axis == "Re" {
        return evaluation.reynolds_elasticity();
    }
    let group = |name: &str| evaluation.groups()[name];
    let nu = evaluation.evidence().value;
    let value = match evaluation.card().id {
        CorrelationId::CircularDuctLaminarCwt | CorrelationId::CircularDuctLaminarChf => 0.0,
        CorrelationId::CircularDuctHausen => match axis {
            "Pr" | "L_over_Dh" => {
                let gz = group("Gz");
                let b = 0.04 * det::pow(gz, 2.0 / 3.0);
                let elasticity = (0.0668 * gz / (1.0 + b)) / nu
                    * (1.0 + b / 3.0) / (1.0 + b);
                if axis == "Pr" { elasticity } else { -elasticity }
            }
            _ => 0.0,
        },
        CorrelationId::RectangularDuctLaminarCwt => {
            if axis == "aspect_ratio" {
                let a = group(axis);
                let slope = (-2.740_f64)
                    .mul_add(a, 10.808)
                    .mul_add(a, -15.357)
                    .mul_add(a, 9.940)
                    .mul_add(a, -2.610);
                7.541 * a * slope / nu
            } else {
                0.0
            }
        }
        CorrelationId::RectangularDuctLaminarChf => {
            if axis == "aspect_ratio" {
                let a = group(axis);
                let slope = (-0.9305_f64)
                    .mul_add(a, 4.2312)
                    .mul_add(a, -7.4295)
                    .mul_add(a, 6.1706)
                    .mul_add(a, -2.0421);
                8.235 * a * slope / nu
            } else {
                0.0
            }
        }
        CorrelationId::RectangularDuctLaminarCwtDevelopingPr072 => return Ok(None),
        CorrelationId::DittusBoelter => {
            if axis == "Pr" {
                match direction {
                    ThermalDirection::HeatingFluid => 0.4,
                    ThermalDirection::CoolingFluid => 0.3,
                }
            } else {
                0.0
            }
        }
        CorrelationId::Gnielinski => {
            if axis == "Pr" {
                let z = 0.79_f64.mul_add(det::ln(group("Re")), -1.64);
                let friction = 1.0 / det::pow(z, 2.0);
                let correction = 12.7 * det::sqrt(friction / 8.0);
                let pr_power = det::pow(group("Pr"), 2.0 / 3.0);
                1.0 - (2.0 / 3.0) * correction * pr_power
                    / (1.0 + correction * (pr_power - 1.0))
            } else {
                0.0
            }
        }
        CorrelationId::FlatPlateLaminarAverage | CorrelationId::FlatPlateTurbulentAverage => {
            if axis == "Pr" { 1.0 / 3.0 } else { 0.0 }
        }
        CorrelationId::ChurchillBernsteinCylinder => {
            if axis == "Pr" {
                let b = det::pow(0.4 / group("Pr"), 2.0 / 3.0);
                (nu - 0.3) / nu * (1.0 / 3.0 + b / (6.0 * (1.0 + b)))
            } else {
                0.0
            }
        }
        CorrelationId::ChurchillChuVerticalPlate => {
            let b = det::pow(0.492 / group("Pr"), 9.0 / 16.0);
            let term = 0.387 * det::pow(group("Ra"), 1.0 / 6.0)
                / det::pow(1.0 + b, 8.0 / 27.0);
            let ra_elasticity = term / (3.0 * (0.825 + term));
            match axis {
                "Ra" => ra_elasticity,
                "Pr" => ra_elasticity * b / (1.0 + b),
                _ => 0.0,
            }
        }
    };
    finite_derivative(value, "Nusselt input elasticity").map(Some)
}

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

    fn scaled_input(mut inputs: CorrelationInputs, axis: &str, factor: f64) -> CorrelationInputs {
        let value = match axis {
            "Re" => &mut inputs.reynolds,
            "Pr" => &mut inputs.prandtl,
            "L_over_Dh" => &mut inputs.length_over_hydraulic_diameter,
            "aspect_ratio" => &mut inputs.aspect_ratio,
            "Ra" => &mut inputs.rayleigh,
            _ => panic!("not an independent input: {axis}"),
        };
        *value.as_mut().expect("supplied input") *= factor;
        inputs
    }

    #[test]
    fn all_smooth_input_partials_match_the_actual_evaluator() {
        let step = 1e-5_f64;
        for (id, re) in [
            (CorrelationId::CircularDuctLaminarCwt, 1_000.0),
            (CorrelationId::CircularDuctLaminarChf, 1_000.0),
            (CorrelationId::CircularDuctHausen, 1_000.0),
            (CorrelationId::RectangularDuctLaminarCwt, 1_000.0),
            (CorrelationId::RectangularDuctLaminarChf, 1_000.0),
            (CorrelationId::DittusBoelter, 50_000.0),
            (CorrelationId::Gnielinski, 20_000.0),
            (CorrelationId::FlatPlateLaminarAverage, 100_000.0),
            (CorrelationId::FlatPlateTurbulentAverage, 1_000_000.0),
            (CorrelationId::ChurchillBernsteinCylinder, 10_000.0),
            (CorrelationId::ChurchillChuVerticalPlate, 0.0),
        ] {
            for pr in [0.72, 1.0, 7.0] {
                for direction in [ThermalDirection::HeatingFluid, ThermalDirection::CoolingFluid] {
                    let inputs = if id == CorrelationId::ChurchillChuVerticalPlate {
                        CorrelationInputs::natural(1e6, pr)
                    } else {
                        CorrelationInputs::forced(re, pr)
                            .with_length_ratio(200.0)
                            .with_aspect_ratio(0.5)
                    }.with_direction(direction);
                    let evaluated = inputs.evaluate_with_sensitivities(id).unwrap();
                    let primal = evaluate(id, inputs).unwrap();
                    assert_eq!(evaluated.evidence().value.to_bits(), primal.evidence().value.to_bits());
                    assert_eq!(evaluated.evidence().model, primal.evidence().model);
                    assert_eq!(evaluated.evidence().numerical, primal.evidence().numerical);
                    assert_eq!(evaluated.source_region(), primal.source_region());
                    assert_ne!(evaluated.evidence().provenance, primal.evidence().provenance);
                    let partials = &evaluated.evidence().sensitivity.d_qoi;
                    for axis in INPUT_AXES {
                        let Some(&input) = evaluated.groups().get(axis) else { continue; };
                        let actual = partials[axis] * input / evaluated.evidence().value;
                        let run = |factor| evaluate(id, scaled_input(inputs, axis, factor)).unwrap();
                        let expected = (run(step.exp()).evidence().value.ln()
                            - run((-step).exp()).evidence().value.ln()) / (2.0 * step);
                        assert!((actual - expected).abs() < 1e-8,
                            "{id:?} {axis} Pr={pr}: {actual} != {expected}");
                    }
                    assert!(!partials.contains_key("Pe"));
                    assert!(!partials.contains_key("Gz"));
                }
            }
        }
    }

    #[test]
    fn prandtl_one_keeps_distinct_heating_and_cooling_slopes() {
        let inputs = CorrelationInputs::forced(50_000.0, 1.0).with_length_ratio(200.0);
        let heat = inputs.evaluate_with_sensitivities(CorrelationId::DittusBoelter).unwrap();
        let cool = inputs.with_direction(ThermalDirection::CoolingFluid)
            .evaluate_with_sensitivities(CorrelationId::DittusBoelter).unwrap();
        assert_eq!(heat.evidence().value.to_bits(), cool.evidence().value.to_bits());
        assert_eq!(heat.evidence().sensitivity.d_qoi["Pr"], 0.4 * heat.evidence().value);
        assert_eq!(cool.evidence().sensitivity.d_qoi["Pr"], 0.3 * cool.evidence().value);
    }

    #[test]
    fn only_partials_that_move_a_domain_endpoint_are_omitted() {
        let inputs = CorrelationInputs::forced(20_000.0, 0.5).with_length_ratio(10.0);
        let value = inputs.evaluate_with_sensitivities(CorrelationId::Gnielinski).unwrap();
        let partials = &value.evidence().sensitivity.d_qoi;
        assert!(partials.contains_key("Re"));
        assert!(!partials.contains_key("Pr"));
        assert!(!partials.contains_key("L_over_Dh"));

        let value = CorrelationInputs::natural(0.1, 0.72)
            .evaluate_with_sensitivities(CorrelationId::ChurchillChuVerticalPlate).unwrap();
        assert!(value.evidence().sensitivity.d_qoi.contains_key("Pr"));
        assert!(!value.evidence().sensitivity.d_qoi.contains_key("Ra"));

        let value = CorrelationInputs::forced(1_000.0, 0.72).with_length_ratio(200.0)
            .with_aspect_ratio(1.0)
            .evaluate_with_sensitivities(CorrelationId::RectangularDuctLaminarCwt).unwrap();
        assert_eq!(value.evidence().sensitivity.d_qoi["Re"], 0.0);
        assert!(!value.evidence().sensitivity.d_qoi.contains_key("aspect_ratio"));
    }

    #[test]
    fn dependent_group_endpoints_block_all_inputs_that_move_them() {
        let mut value = evaluate(CorrelationId::CircularDuctHausen,
            CorrelationInputs::forced(1_000.0, 0.72).with_length_ratio(200.0)).unwrap();
        // Narrow the card to place Gz at its endpoint while each primitive
        // input remains interior; this isolates the derived-axis chain rule.
        value.card.model.validity = value.card.model.validity.clone()
            .with("Gz", 0.01, value.groups()["Gz"]);
        for axis in ["Re", "Pr", "L_over_Dh"] {
            assert_eq!(input_elasticity(&value, ThermalDirection::HeatingFluid, axis).unwrap(), None);
        }
    }

    #[test]
    fn table_and_bridge_do_not_acquire_fabricated_input_derivatives() {
        for length_ratio in [24.0, 200.0] {
            let value = CorrelationInputs::forced(1_000.0, 0.72)
                .with_length_ratio(length_ratio).with_aspect_ratio(0.5)
                .evaluate_with_sensitivities(CorrelationId::RectangularDuctLaminarCwtDevelopingPr072)
                .unwrap();
            assert!(value.evidence().sensitivity.d_qoi.is_empty());
        }
        assert!(CorrelationInputs::forced(f64::NAN, 0.72).with_length_ratio(200.0)
            .evaluate_with_sensitivities(CorrelationId::Gnielinski).is_err());
        assert!(CorrelationInputs::forced(1_000.0, 0.72).with_length_ratio(200.0)
            .evaluate_with_sensitivities(CorrelationId::Gnielinski).is_err());
    }

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
