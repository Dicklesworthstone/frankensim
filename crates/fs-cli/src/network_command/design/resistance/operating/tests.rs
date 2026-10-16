use super::*;
use fs_exec::{Budget, ExecMode, StreamKey};

const MODEL: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/resistance-model.json"));
fn model() -> Request { Request::parse(MODEL).unwrap() }
fn policy(limit: f64, evaluations: usize) -> String {
    format!(r#"{{"schema":"frankensim.cooling-resistance-design.v1",
        "temperature_limit_k":{limit},"temperature_scale_k":0.1,
        "controls":[{{"branch":"first-stream","minimum_multiplier":0.25,"maximum_multiplier":4,"weight":1}},
                    {{"branch":"fresh-bypass","minimum_multiplier":0.25,"maximum_multiplier":4,"weight":0.3}}],
        "max_evaluations":{evaluations},"max_iterations":80,"max_multiplier_updates":40,
        "max_backtracks":24,"kkt_tolerance":0.000001,
        "operating_envelope":{{"smoothing_k":0.02,"power_relative_tolerance":1e-12,"cases":[
            {{"name":"hot","source_multiplier":1.3,"inlet_temperatures":[{{"node":0,"temperature_k":301}}]}},
            {{"name":"light","source_multiplier":0.8}}]}}}}"#)
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 41, kernel_id: 718, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic)))
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance * b.abs().max(1.0), "{a:.16e} != {b:.16e}");
}

#[test]
fn envelope_is_unnormalized_and_differentiable_at_ties() {
    let (value, weights) = envelope(&[300.0; 3], 0.02).unwrap();
    close(value, 300.0 + 0.02 * 3.0_f64.ln(), 1e-14);
    for weight in weights { close(weight, 1.0/3.0, 1e-14); }
    for tau in [1e-10, 0.02, 10.0] {
        let temperatures = [300.0, 300.1, 299.0];
        let (value, weights) = envelope(&temperatures, tau).unwrap();
        assert!(value >= 300.1);
        assert!(value <= 300.1 + tau * 3.0_f64.ln() + 1e-12);
        close(weights.iter().sum(), 1.0, 1e-14);
    }
    let (value, _) = envelope(&[300.0, 301.0], 1e-300).unwrap();
    assert_eq!(value, 301.0);
    assert!(envelope(&[], 0.1).is_err());
    assert!(envelope(&[300.0], 0.0).is_err());
    assert!(envelope(&[f64::NAN], 0.1).is_err());
}

#[test]
fn case_sources_use_original_power_footprints_and_restore_nominal() {
    let mut request = model();
    let original = request.graph.clone();
    let design = Specification::parse(&policy(1000.0, 6), &request).unwrap();
    let cases = design.operating_envelope.as_ref().unwrap();
    assert_eq!(cases.len(), 3);
    for (index, watts) in [(0, 1.0), (1, 1.3), (2, 0.8)] {
        cases.cases[index].loads.install(&mut request);
        let audit = request.solid_data.power.as_ref().unwrap();
        close(audit.delivered_total_w(), watts, 1e-12);
        assert_eq!(request.graph.branches(), original.branches());
        assert_eq!(request.solid_data.component_map.as_ref().unwrap().components()[0].vertices(), &[4]);
    }
    with_cx(|cx| {
        let a = evaluate_candidate(&mut request, &original, &design, &[0.1, -0.2], cx, false).unwrap();
        let b = evaluate_candidate(&mut request, &original, &design, &[0.1, -0.2], cx, false).unwrap();
        assert_eq!(a.sample, b.sample);
        assert_eq!(a.signatures, b.signatures);
        assert_eq!(request.inlets[0].temperature.value(), 300.0);
        close(request.solid_data.power.as_ref().unwrap().delivered_total_w(), 1.0, 1e-12);
    });
}

#[test]
fn production_case_gradient_matches_complete_perturbed_sweeps() {
    with_cx(|cx| {
        let mut request = model();
        let original = request.graph.clone();
        let design = Specification::parse(&policy(304.0, 240), &request).unwrap();
        let point = [0.1, -0.2];
        let result = evaluate_candidate(&mut request, &original, &design, &point, cx, false).unwrap();
        for index in 0..2 {
            let h = 1e-3;
            let mut plus = point; plus[index] += h;
            let mut minus = point; minus[index] -= h;
            let a = evaluate_candidate(&mut request, &original, &design, &plus, cx, false).unwrap();
            let b = evaluate_candidate(&mut request, &original, &design, &minus, cx, false).unwrap();
            close(result.sample.constraint_gradient[index], (a.sample.constraint-b.sample.constraint)/(2.0*h), 3e-4);
        }
        let nominal = result.signatures[0].constraint_temperature;
        let hot = result.signatures[1].constraint_temperature;
        assert!(hot > nominal);
        request.graph = original.clone();
        let failing = Specification::parse(&policy(0.5*(nominal+hot), 6), &request).unwrap();
        let result = evaluate_candidate(&mut request, &original, &failing, &point, cx, false).unwrap();
        assert!(result.physical.objective < failing.limit);
        assert!(!result.passes);
        assert!(result.sample.constraint > 0.0);
    });
}

#[test]
fn every_replayed_case_is_reported_and_charged_without_multiplying_cost() {
    with_cx(|cx| {
        let request = model();
        let text = policy(1000.0, 6);
        let design = Specification::parse(&text, &request).unwrap();
        assert_eq!(design.options.max_evaluations, 1);
        let output = J::parse(&solve(request, design, cx).unwrap()).unwrap();
        let design = output.get("resistance_design").unwrap();
        assert_eq!(design.f64_field("complete_evaluations"), Some(6.0));
        assert_eq!(design.f64_field("replay_evaluations"), Some(3.0));
        assert_eq!(design.f64_field("selected_cost"), Some(0.0));
        assert_eq!(design.str_field("status"), Some("locally-converged"));
        let report = design.get("operating_envelope").unwrap();
        let cases = report.get("cases").unwrap().as_array().unwrap();
        assert_eq!(cases.len(), 3);
        for (case, watts) in cases.iter().zip([1.0, 1.3, 0.8]) {
            let physical = case.get("result").unwrap();
            close(physical.f64_field("source_w").unwrap(), watts, 1e-8);
            assert!(physical.get("solid_temperatures_k").unwrap().as_array().is_some());
        }
        assert_eq!(cases[0].get("result"), Some(&output_without_design(&output)));
        close(output.f64_field("source_w").unwrap(), 1.0, 1e-8);

        let mut request = model();
        let original = request.graph.clone();
        let design = Specification::parse(&text.replace("\"operating_envelope\":", "\"air_power_cost\":{\"weight_per_w\":10},\"operating_envelope\":"), &request).unwrap();
        let point = [0.1, -0.2];
        let one = evaluate(&mut request, &original, &design, &point, cx).unwrap().0;
        let all = evaluate_candidate(&mut request, &original, &design, &point, cx, false).unwrap();
        assert_eq!(one.objective, all.sample.objective);
        assert_eq!(one.gradient, all.sample.gradient);
    });
}
fn output_without_design(output: &J) -> J {
    let J::Object(fields) = output else { panic!("expected object") };
    J::Object(fields.iter().filter(|(key,_)| key != "resistance_design").cloned().collect())
}

#[test]
fn case_policy_and_total_evaluation_budget_refuse_ambiguous_inputs() {
    let request = model();
    let text = policy(304.0, 6);
    for (from, to) in [
        ("\"max_evaluations\":6", "\"max_evaluations\":5"),
        ("\"smoothing_k\":0.02", "\"smoothing_k\":0"),
        ("\"power_relative_tolerance\":1e-12", "\"power_relative_tolerance\":1"),
        ("\"name\":\"hot\"", "\"name\":\"nominal\""),
        ("\"name\":\"light\"", "\"name\":\"hot\""),
        ("\"node\":0", "\"node\":2"),
        ("\"temperature_k\":301", "\"temperature_k\":0"),
        ("\"source_multiplier\":1.3", "\"source_multiplier\":-1"),
    ] {
        assert!(text.contains(from));
        assert!(Specification::parse(&text.replace(from, to), &request).is_err(), "{to}");
    }
    with_cx(|cx| {
        let design = Specification::parse(&policy(1.0, 6), &request).unwrap();
        assert_eq!(solve(request, design, cx).unwrap_err().code, "cooling-network-design-no-feasible-sample");
    });
}

#[test]
fn cancelled_case_sweep_returns_no_candidate() {
    let mut request = model();
    let original = request.graph.clone();
    let design = Specification::parse(&policy(1000.0, 6), &request).unwrap();
    let gate = CancelGate::new_clock_free(); gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 41, kernel_id: 718, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        assert!(evaluate_candidate(&mut request, &original, &design, &[0.0; 2], &cx, true).is_err());
    });
}

#[test]
fn component_cases_move_heat_between_footprints_without_changing_geometry() {
    let text = MODEL.replace("\"watts\": 1, \"vertices\": [4]", "\"watts\": 0.5, \"vertices\": [4]")
        .replace("\"vertices\": [4] }", "\"vertices\": [4] }, {\"name\":\"memory\",\"watts\":0.5,\"vertices\":[7]}");
    assert!(text.contains("memory"));
    let request = Request::parse(&text).unwrap();
    let value = J::parse(r#"{"smoothing_k":0.02,"power_relative_tolerance":1e-12,"cases":[
        {"name":"chip-only","component_powers_w":{"chip":1,"memory":0}},
        {"name":"memory-only","component_powers_w":{"chip":0,"memory":1}}]}"#).unwrap();
    let cases = Cases::parse(&value, &request, 1000.0, None).unwrap();
    for (index, name) in [(1, "chip"), (2, "memory")] {
        let loads = &cases.cases[index].loads;
        close(loads.audit.as_ref().unwrap().delivered_total_w(), 1.0, 1e-12);
        let active = loads.map.as_ref().unwrap().components().iter().find(|row| row.watts() == 1.0).unwrap();
        assert_eq!(active.name(), name);
        assert_eq!(active.vertices(), if name == "chip" { &[4][..] } else { &[7][..] });
        assert_eq!(active.uncertainty(), PowerUncertainty::Unstated);
    }
    assert_ne!(cases.cases[1].loads.nodal, cases.cases[2].loads.nodal);
    close(request.solid_data.component_map.as_ref().unwrap().components()[0].watts(), 0.5, 0.0);
    let partial = J::parse(r#"{"chip":0.25}"#).unwrap();
    let loads = Loads::capture(&request).with_component_powers(&partial, &request.mesh, 1e-12).unwrap();
    close(loads.audit.as_ref().unwrap().delivered_total_w(), 0.75, 1e-12);
    for bad in [r#"{}"#, r#"{"missing":1}"#, r#"{"chip":-1}"#] {
        assert!(Loads::capture(&request).with_component_powers(&J::parse(bad).unwrap(), &request.mesh, 1e-12).is_err());
    }
    let mixed = J::parse(r#"{"smoothing_k":0.02,"power_relative_tolerance":1e-12,"cases":[
        {"name":"ambiguous","source_multiplier":1,"component_powers_w":{"chip":1}}]}"#).unwrap();
    assert!(Cases::parse(&mixed, &request, 1000.0, None).is_err());
}

#[test]
fn failed_case_sweep_restores_nominal_sources_and_supplies() {
    with_cx(|cx| {
        let mut request = model();
        let original = request.graph.clone();
        let mut design = Specification::parse(&policy(1000.0, 6), &request).unwrap();
        // Corrupt a private test-only case to force a real producer refusal
        // after nominal succeeded. The public case parser cannot create this.
        design.operating_envelope.as_mut().unwrap().cases[1].loads.inlets.clear();
        assert!(evaluate_candidate(&mut request, &original, &design, &[0.0;2], cx, true).is_err());
        assert_eq!(request.inlets[0].temperature.value(), 300.0);
        close(request.solid_data.power.as_ref().unwrap().delivered_total_w(), 1.0, 1e-12);
    });
}

fn requirement_policy(limit: f64, local_limit: f64, evaluations: usize) -> String {
    format!(r#"{{"schema":"frankensim.cooling-resistance-design.v1",
        "temperature_limit_k":{limit},"temperature_scale_k":0.1,
        "controls":[{{"branch":"first-stream","minimum_multiplier":0.25,"maximum_multiplier":4,"weight":1}},
                    {{"branch":"fresh-bypass","minimum_multiplier":0.25,"maximum_multiplier":4,"weight":0.3}}],
        "max_evaluations":{evaluations},"max_iterations":80,"max_multiplier_updates":40,
        "max_backtracks":24,"kkt_tolerance":0.000001,
        "operating_envelope":{{"smoothing_k":0.02,"power_relative_tolerance":1e-12,"cases":[
          {{"name":"chip","objective":{{"max_vertices":[4,7]}},"maximum_smoothing_k":0.03,
            "source_multiplier":1.2,"temperature_limit_k":{local_limit}}},
          {{"name":"sink","objective":{{"mean_wall_region":"last-face"}},
            "source_multiplier":0.7,"temperature_limit_k":{local_limit}}}]}}}}"#)
}

#[test]
fn mixed_requirements_use_their_own_seeds_limits_and_replay_reports() {
    with_cx(|cx| {
        let mut request = model();
        let original = request.graph.clone();
        let design = Specification::parse(&requirement_policy(1000.0, 900.0, 6), &request).unwrap();
        let candidate = evaluate_candidate(&mut request, &original, &design, &[0.0;2], cx, true).unwrap();
        assert!(candidate.passes);
        assert!(request.objective.is_mean());
        assert_eq!(request.objective.region(), Some("first-face"));
        let report = J::parse(candidate.report.as_ref().unwrap()).unwrap();
        let cases = report.get("cases").unwrap().as_array().unwrap();
        assert_eq!(cases[0].f64_field("temperature_limit_k"), Some(1000.0));
        assert_eq!(cases[1].f64_field("temperature_limit_k"), Some(900.0));
        let chip = cases[1].path(&["result","objective"]).unwrap();
        assert_eq!(chip.str_field("kind"), Some("smooth_upper_envelope"));
        assert_eq!(chip.str_field("maximum_selection"), Some("max_vertices"));
        assert_eq!(chip.f64_field("maximum_smoothing_k"), Some(0.03));
        let sink = cases[2].path(&["result","objective"]).unwrap();
        assert_eq!(sink.str_field("kind"), Some("mean_wall_region"));
        assert_eq!(sink.str_field("region"), Some("last-face"));
        close(candidate.signatures[1].sample.constraint,
            (candidate.signatures[1].constraint_temperature-900.0)/0.1, 1e-14);
        let output = J::parse(&solve(request, design, cx).unwrap()).unwrap();
        assert!(output.path(&["resistance_design","operating_envelope","cases"]).unwrap().as_array().unwrap()
            .iter().all(|row| row.path(&["result","solid_temperatures_k"]).is_some()));
    });
}

#[test]
fn cooler_component_with_a_tighter_limit_cannot_hide_behind_absolute_temperature() {
    with_cx(|cx| {
        let mut request = model();
        let original = request.graph.clone();
        let mut design = Specification::parse(&requirement_policy(1000.0, 900.0, 240), &request).unwrap();
        let first = evaluate_candidate(&mut request, &original, &design, &[0.0;2], cx, false).unwrap();
        assert!(first.signatures[2].constraint_temperature < first.signatures[1].constraint_temperature);
        // All other requirements pass comfortably, while the cooler sink has
        // a positive margin. The hottest absolute temperature is irrelevant.
        let cases = &mut design.operating_envelope.as_mut().unwrap().cases;
        cases[1].limit = first.signatures[1].constraint_temperature + 1.0;
        cases[2].limit = first.signatures[2].constraint_temperature - 0.01;
        let candidate = evaluate_candidate(&mut request, &original, &design, &[0.0;2], cx, true).unwrap();
        assert!(!candidate.passes);
        assert!(candidate.sample.constraint >= 0.099999);
        let report = J::parse(candidate.report.as_ref().unwrap()).unwrap();
        assert!(report.f64_field("constraint_excess_k").unwrap() > 0.0);
        let weights = report.get("case_weights").unwrap().as_array().unwrap();
        assert!(weights[2].as_f64().unwrap() > 0.999);
        // Check every physical control, with the case limits held fixed.
        for index in 0..2 {
            let mut plus=[0.0;2]; plus[index]=1e-3;
            let mut minus=[0.0;2]; minus[index]=-1e-3;
            let a = evaluate_candidate(&mut request,&original,&design,&plus,cx,false).unwrap();
            let b = evaluate_candidate(&mut request,&original,&design,&minus,cx,false).unwrap();
            close(candidate.sample.constraint_gradient[index], (a.sample.constraint-b.sample.constraint)/0.002, 3e-4);
        }
    });
}

#[test]
fn mixed_objective_admission_does_not_invent_smoothing_or_disable_gradients() {
    let request = model();
    let text = requirement_policy(1000.0, 900.0, 6);
    for (from, to) in [
        (",\"maximum_smoothing_k\":0.03", ""),
        ("\"temperature_limit_k\":900", "\"temperature_limit_k\":0"),
        ("\"max_vertices\":[4,7]", "\"max_vertices\":[]"),
        ("\"max_vertices\":[4,7]", "\"max_vertices\":[4,7],\"gradient\":false"),
        ("\"mean_wall_region\":\"last-face\"", "\"mean_wall_region\":\"unknown\""),
        ("\"name\":\"sink\"", "\"name\":\"sink\",\"maximum_smoothing_k\":0.1"),
    ] {
        assert!(text.contains(from));
        assert!(Specification::parse(&text.replace(from,to), &request).is_err(), "{to}");
    }
    // A mean case under a peak nominal must not inherit peak semantics. A
    // peak case without its own tau deliberately inherits the declared tau.
    let mut request = model();
    request.objective = objective::Objective::parse(&J::parse(r#"{"max_solid_temperature":true}"#).unwrap(),
        &request.surfaces,&request.mesh).unwrap();
    let value = J::parse(r#"{"smoothing_k":0.02,"power_relative_tolerance":1e-12,"cases":[
        {"name":"mean","objective":{"mean_wall_region":"first-face"}},
        {"name":"peak","objective":{"max_vertices":[4]}}]}"#).unwrap();
    let cases = Cases::parse(&value,&request,1000.0,Some(0.04)).unwrap();
    assert!(cases.cases[1].loads.objective.is_mean());
    match &cases.cases[2].loads.objective {
        objective::Objective::Maximum { smoothing_k, .. } => assert_eq!(*smoothing_k,Some(0.04)),
        _ => panic!("peak case lost its selection"),
    }
}

#[test]
fn optimizer_replays_a_nontrivial_design_against_each_component_requirement() {
    with_cx(|cx| {
        let mut request = model();
        let original = request.graph.clone();
        let mut design = Specification::parse(&requirement_policy(1000.0,900.0,720),&request).unwrap();
        let initial=evaluate_candidate(&mut request,&original,&design,&[0.0;2],cx,false).unwrap();
        let improved=evaluate_candidate(&mut request,&original,&design,&[-0.7,0.7],cx,false).unwrap();
        assert!(improved.signatures[1].constraint_temperature < initial.signatures[1].constraint_temperature);
        let limit=0.5*(improved.signatures[1].constraint_temperature+initial.signatures[1].constraint_temperature);
        design.operating_envelope.as_mut().unwrap().cases[1].limit=limit;
        request.graph=original;
        let result=J::parse(&solve(request,design,cx).unwrap()).unwrap();
        let design=result.get("resistance_design").unwrap();
        assert!(design.f64_field("selected_cost").unwrap()>0.0);
        assert!(design.f64_field("complete_evaluations").unwrap()<=720.0);
        for case in design.path(&["operating_envelope","cases"]).unwrap().as_array().unwrap() {
            let limit=case.f64_field("temperature_limit_k").unwrap();
            let objective=case.path(&["result","objective"]).unwrap();
            assert!(objective.f64_field("value_k").unwrap()<=limit);
            if let Some(actual)=objective.f64_field("actual_maximum_k") { assert!(actual<=limit); }
        }
    });
}
