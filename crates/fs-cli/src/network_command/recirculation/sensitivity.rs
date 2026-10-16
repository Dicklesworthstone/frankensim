//! Return-fraction controls recovered from the existing total fresh-air adjoint.
//!
//! For supply i, x_i = (1 - sum_j r_ij) F_i + sum_j r_ij T_j.
//! The total fresh-temperature gradient is (1 - sum_j r_ij) lambda_i.
//! Hence dJ/dr_ij = lambda_i (T_j - F_i), including the same implicit
//! air/solid/radiation feedback. This identity is only used for the existing
//! temperature objectives, which have no explicit fresh-temperature term.
use super::*;

impl Policy {
    pub(super) fn gradients(&self, request: &Request, evaluated: &Evaluation)
        -> Result<Option<Vec<f64>>>
    {
        let Some(gradient) = &evaluated.gradient else { return Ok(None); };
        let march = &evaluated.coupled.transport;
        if gradient.inlets.len() != march.node_temperatures_k.len() {
            return Err(producer("return-fraction adjoint has the wrong inlet dimension"));
        }
        let mut totals = BTreeMap::<usize, f64>::new();
        for link in &self.links {
            let total = totals.entry(link.supply_node).or_default();
            *total = checked(*total + link.fraction)?;
        }
        let mut values = Vec::with_capacity(self.links.len());
        for link in &self.links {
            let fresh = request.inlets.iter().find(|inlet| inlet.node == link.supply_node)
                .ok_or_else(|| bad("return-fraction adjoint lacks fresh supply"))?.temperature.value();
            let returned = march.node_temperatures_k.get(link.return_node).copied().flatten()
                .ok_or_else(|| bad("return-fraction adjoint lacks transported exhaust"))?;
            let makeup = checked(1.0 - totals[&link.supply_node])?;
            if makeup <= 0.0 { return Err(bad("return-fraction adjoint needs positive fresh makeup")); }
            let inlet = *gradient.inlets.get(link.supply_node)
                .ok_or_else(|| bad("return-fraction adjoint lacks supply slot"))?;
            let multiplier = checked(inlet / makeup)?;
            values.push(checked(multiplier * checked(returned - fresh)?)?);
        }
        Ok(Some(values))
    }

    pub(super) fn sensitivity_json(&self, request: &Request, evaluated: &Evaluation) -> Result<String> {
        let Some(values) = self.gradients(request, evaluated)? else { return Ok("null".into()); };
        let rows = self.links.iter().zip(values).map(|(link, value)| Ok(format!(
            "{{\"supply_node\":{},\"return_node\":{},\"fraction\":{},\"dobjective_dfraction_k\":{}}}",
            link.supply_node, link.return_node, num(link.fraction)?, num(value)?)))
            .collect::<Result<Vec<_>>>()?.join(",");
        Ok(format!("{{\"method\":\"implicit-return-mixer-adjoint\",\"links\":[{rows}],\"scope\":\"selected steady temperature objective; full feedback retained from its total fresh-inlet adjoint; fraction controls displace fresh makeup at fixed hydraulic flows and all other model inputs; derivative per unit fraction, not per percent; zero-fraction values are right-hand model derivatives, usable only when makeup and exhaust-draw constraints admit the perturbation; peak objectives retain the producer's active-vertex limitation; no finite-change, global monotonicity or hardware-compliance guarantee\"}}"))
    }
}

fn checked(value: f64) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(producer("nonfinite return-fraction sensitivity")) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::super::tests::{with_cx, close};

    fn request() -> Request {
        Request::parse(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../examples/cooling-network/recirculated-slab.json"))).unwrap()
    }
    fn evaluate(request: &Request, gradient: bool) -> Evaluation {
        with_cx(|cx| {
            let flow = request.flow(cx).unwrap();
            let h = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
            match &request.radiation {
                Some(policy) => policy.evaluate(request,cx,&flow,&h,gradient).unwrap().value,
                None => request.evaluate(cx, &flow, &h, gradient).unwrap(),
            }
        })
    }

    #[test]
    fn fraction_adjoints_match_complete_perturbed_fem_with_both_signs() {
        let mut r = request();
        let value = evaluate(&r, true);
        let derivatives = r.recirculation.as_ref().unwrap().gradients(&r, &value).unwrap().unwrap();
        assert!(derivatives[0] < 0.0 && derivatives[1] > 0.0);
        for (i, expected) in derivatives.into_iter().enumerate() {
            let original = r.recirculation.as_ref().unwrap().links[i].fraction;
            let delta = 1e-5;
            r.recirculation.as_mut().unwrap().links[i].fraction = original + delta;
            let plus = evaluate(&r, false).objective;
            r.recirculation.as_mut().unwrap().links[i].fraction = original - delta;
            let minus = evaluate(&r, false).objective;
            r.recirculation.as_mut().unwrap().links[i].fraction = original;
            close(expected, (plus-minus)/(2.0*delta), 1e-4);
        }
        let output = execute(&r, &CancelGate::new_clock_free()).unwrap();
        assert_eq!(J::parse(&output).unwrap().path(&["recirculation_sensitivity", "links"])
            .unwrap().as_array().unwrap().len(), 2);
    }

    #[test]
    fn zero_returns_have_one_sided_sensitivity_without_inventing_a_primal_return() {
        let mut r = request();
        for link in &mut r.recirculation.as_mut().unwrap().links { link.fraction = 0.0; }
        let value = evaluate(&r, true);
        assert!(value.coupled.transport.recirculation.is_none());
        let expected = r.recirculation.as_ref().unwrap().gradients(&r, &value).unwrap().unwrap();
        for (i, derivative) in expected.into_iter().enumerate() {
            r.recirculation.as_mut().unwrap().links[i].fraction = 1e-6;
            let perturbed = evaluate(&r, false);
            close(derivative, (perturbed.objective-value.objective)/1e-6, 2e-3);
            r.recirculation.as_mut().unwrap().links[i].fraction = 0.0;
        }
    }
    #[test]
    fn fraction_derivative_keeps_nonlinear_radiation_feedback() {
        let mut r = request();
        let policy = J::parse(r#"{"max_iterations":200,"temperature_tolerance_k":1e-10,
            "relaxation":1,"surfaces":[{"surface":"first-face","emissivity":0.8,
            "ambient_temperature_k":300,"source":"declared regression surface"}]}"#).unwrap();
        r.radiation = Some(radiation::Policy::parse(&policy,&J::Null,&r.surfaces).unwrap());
        let value = evaluate(&r,true);
        let expected = r.recirculation.as_ref().unwrap().gradients(&r,&value).unwrap().unwrap()[0];
        let delta = 1e-5;
        r.recirculation.as_mut().unwrap().links[0].fraction = 0.6+delta;
        let plus = evaluate(&r,false).objective;
        r.recirculation.as_mut().unwrap().links[0].fraction = 0.6-delta;
        let minus = evaluate(&r,false).objective;
        close(expected,(plus-minus)/(2.0*delta),2e-4);
    }

}
