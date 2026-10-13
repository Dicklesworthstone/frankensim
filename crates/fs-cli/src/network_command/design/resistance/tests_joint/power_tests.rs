use super::*;

#[test]
fn power_cost_uses_both_product_rule_terms_without_changing_thermal_gradients() {
    with_cx(|cx| {
        let mut request = model(true);
        let original = request.graph.clone();
        let text = policy(303.0, true, false);
        let plain = Specification::parse(&text, &request).unwrap();
        let text = text.replace("\"max_evaluations\":240", "\"air_power_cost\":{\"weight_per_w\":10},\"max_evaluations\":240");
        let design = Specification::parse(&text, &request).unwrap();
        let point = [0.1, -0.2, 0.15];
        let (unweighted, _, _) = evaluate(&mut request, &original, &plain, &point, cx).unwrap();
        let (sample, flow, _) = evaluate(&mut request, &original, &design, &point, cx).unwrap();
        assert_eq!(sample.constraint, unweighted.constraint);
        assert_eq!(sample.constraint_gradient, unweighted.constraint_gradient);
        let air_power = power::value(&request, &flow).unwrap();
        close(sample.objective - unweighted.objective, 10.0 * air_power, 1e-12);
        close(sample.gradient[2] - unweighted.gradient[2], 30.0 * air_power, 1e-10);
        for i in 0..3 {
            let h = 1e-4;
            let mut plus = point; plus[i] += h;
            let mut minus = point; minus[i] -= h;
            let (a, _, _) = evaluate(&mut request, &original, &design, &plus, cx).unwrap();
            let (b, _, _) = evaluate(&mut request, &original, &design, &minus, cx).unwrap();
            close(sample.gradient[i], (a.objective-b.objective)/(2.0*h), 2e-5);
        }
    });
}

#[test]
fn air_power_preference_returns_a_lower_power_replayed_feasible_design() {
    with_cx(|cx| {
        let request = model(false);
        let original = request.flow(cx).unwrap();
        let original_power = power::value(&request, &original).unwrap();
        let text = policy(1000.0, false, true).replace("\"max_evaluations\":240",
            "\"air_power_cost\":{\"weight_per_w\":10},\"max_evaluations\":240");
        let design = Specification::parse(&text, &request).unwrap();
        let output = J::parse(&solve(request, design, cx).unwrap()).unwrap();
        let report = output.get("resistance_design").unwrap();
        let cost = report.get("air_power_cost").unwrap();
        let selected = cost.f64_field("selected_air_power_w").unwrap();
        assert!(selected < original_power);
        close(selected, output.path(&["fan", "air_power_w"]).unwrap().as_f64().unwrap(), 1e-13);
        assert_eq!(cost.get("electrical_power_w"), Some(&J::Null));
        assert!(report.f64_field("selected_cost").unwrap() < 10.0 * original_power);
        assert!(report.f64_field("selected_temperature_k").unwrap() <= 1000.0);
    });
}

#[test]
fn power_cost_requires_a_fan_and_positive_explicit_weight() {
    let request = model(false);
    for text in [r#"{"weight_per_w":0}"#, r#"{"weight_per_w":-1}"#,
        r#"{"weight_per_w":1,"efficiency":0.5}"#] {
        assert!(power::parse(Some(&J::parse(text).unwrap()), &request).is_err());
    }
    let mut request = model(false); request.fan = None;
    assert!(power::parse(Some(&J::parse(r#"{"weight_per_w":1}"#).unwrap()), &request).is_err());
    assert_eq!(power::parse(None, &request).unwrap(), 0.0);
}
