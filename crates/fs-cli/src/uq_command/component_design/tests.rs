use super::*;
// Compacted so fixture edits are independent of the example's formatting.
static BASE: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| crate::uq_command::json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/adjoint-component-contact-pulse.json"))));
const SPEC: &str = r#"{"schema":"frankensim.cooling-component-design.v1","units":"SI",
    "temperature_limit_k":304,"power_tolerance_w":0.0001,"temperature_tolerance_k":0.0001,
    "max_evaluations":128,"max_total_steps":100000,"wall_seconds":60,
    "priority":[{"component":"chip","interval":0,"min_power_w":0,"max_power_w":8},
                {"component":"memory","interval":0,"min_power_w":0,"max_power_w":8}]}"#;
fn plan() -> Plan { Plan::parse(&J::parse(&BASE).unwrap(),&J::parse(SPEC).unwrap()).unwrap() }
fn deadline() -> Instant { Instant::now()+Duration::from_secs(60) }
fn fake(plan: &Plan, values: &[f64]) -> Evaluation {
    Evaluation {document:J::Null,peak:300.0+2.0*values[0]+0.25*values[1],peak_time:20.0,
        steps:plan.planned_steps,solves:7,slopes:vec![Some(2.0),Some(0.25)],constraints:None}
}
fn interval<'a>(base: &'a J, index: usize) -> &'a J {
    &base.path(&["transient","intervals"]).unwrap().as_array().unwrap()[index]
}

#[test]
fn controls_keep_footprints_other_component_watts_and_absolute_overrides() {
    let base=J::parse(&BASE).unwrap();
    let spec=J::parse(&SPEC.replace("\"component\":\"memory\",\"interval\":0",
        "\"component\":\"standby\",\"interval\":1")).unwrap();
    let plan=Plan::parse(&base,&spec).unwrap();
    let request=plan.request(&[4.0,2.0]).unwrap();
    for key in ["solid","hydraulics","air","radiation","budgets","objective"] {
        assert_eq!(request.get(key),base.get(key));
    }
    assert_eq!(interval(&request,0).path(&["component_powers_w","memory"]).unwrap().as_f64(),Some(3.75));
    assert_eq!(interval(&request,0).path(&["component_powers_w","chip"]).unwrap().as_f64(),Some(4.0));
    assert_eq!(interval(&request,1).path(&["component_powers_w","chip"]).unwrap().as_f64(),Some(0.0));
    assert_eq!(interval(&request,1).path(&["component_powers_w","standby"]).unwrap().as_f64(),Some(2.0));
    assert_eq!(plan.planned_steps,14);
    assert_eq!(base,J::parse(&BASE).unwrap());
    assert!(plan.request(&[9.0,2.0]).is_err());
    assert!(plan.request(&[f64::NAN,2.0]).is_err());
}

#[test]
fn incompatible_or_ambiguous_policies_refuse_before_any_candidate() {
    let base=J::parse(&BASE).unwrap();
    for text in [SPEC.replace("\"memory\"","\"chip\""),SPEC.replace("\"memory\"","\"unknown\""),
        SPEC.replace("\"interval\":0","\"interval\":99"),SPEC.replace("\"min_power_w\":0","\"min_power_w\":-1"),
        SPEC.replace("\"max_evaluations\":128","\"max_evaluations\":0")] {
        assert!(Plan::parse(&base,&J::parse(&text).unwrap()).is_err());
    }
    let spec=J::parse(SPEC).unwrap();
    for text in [BASE.replace("\"qoi\":\"sampled-peak\"","\"qoi\":\"final\""),
        BASE.replace("\"component_power\":true","\"component_power\":false"),
        BASE.replace("\"max_step_s\":2","\"adaptive\":{},\"max_step_s\":2"),
        BASE.replace("\"max_step_s\":2","\"time_convergence\":{},\"max_step_s\":2")] {
        assert!(Plan::parse(&J::parse(&text).unwrap(),&spec).is_err());
    }
}

#[test]
fn explicit_interval_counts_admit_the_same_nested_grid_as_the_producer() {
    let text=BASE.replace("\"max_step_s\":2","\"max_step_s\":1")
        .replace("\"duration_s\":6","\"duration_s\":2.3,\"steps\":6");
    let p=Plan::parse(&J::parse(&text).unwrap(),&J::parse(SPEC).unwrap()).unwrap();
    assert_eq!(p.planned_steps,28);
    assert!(Plan::parse(&J::parse(&text.replace("\"steps\":6","\"steps\":2")).unwrap(),&J::parse(SPEC).unwrap()).is_err());
}

#[test]
fn priority_order_is_physical_policy_and_every_returned_vector_was_evaluated() {
    let p=plan();
    let result=search::allocate(&p,deadline(),|values|Ok(fake(&p,values))).unwrap();
    assert!(result.reason.is_none());assert_eq!(result.completed,2);
    assert!(result.values[0]>1.9999 && result.values[0]<=2.0);
    assert!(result.values[1]<0.001);
    assert!(result.passing.peak<=p.limit);
    assert!(result.history.iter().any(|row|row.get("power_w")==Some(&number_array(&result.values).unwrap())));
    let reverse=search::allocate(&p,deadline(),|v| {
        let mut e=fake(&p,&[v[1],v[0]]);e.slopes=vec![Some(0.25),Some(2.0)];Ok(e)
    }).unwrap();
    assert_eq!(reverse.values[0],8.0);
    assert!((reverse.values[1]-1.0).abs()<0.0001);
}

#[test]
fn cumulative_limits_preserve_only_the_last_complete_passing_allocation() {
    let mut p=plan();p.max_evaluations=1;
    let result=search::allocate(&p,deadline(),|v|Ok(fake(&p,v))).unwrap();
    assert!(result.reason.is_some());assert_eq!(result.completed,0);assert_eq!(result.attempted,1);
    assert_eq!(result.values,vec![0.0,0.0]);assert_eq!(result.history.len(),1);
    p.max_evaluations=128;p.max_total_steps=p.planned_steps;
    let result=search::allocate(&p,deadline(),|v|Ok(fake(&p,v))).unwrap();
    assert_eq!(result.steps,p.planned_steps);assert_eq!(result.attempted,1);assert!(result.reason.is_some());
    p.max_total_steps=100000;
    let mut count=0;
    let result=search::allocate(&p,deadline(),|v| {count+=1;if count==2 {Err(budget("interrupt"))} else {Ok(fake(&p,v))}}).unwrap();
    assert_eq!(result.values,vec![0.0,0.0]);assert_eq!(result.attempted,2);assert_eq!(result.history.len(),1);
    let mut count=0;
    assert!(search::allocate(&p,deadline(),|v| {count+=1;if count==2 {Err(model_failure("bad physics"))} else {Ok(fake(&p,v))}}).is_err());
}

#[test]
fn gradient_proposals_handle_zero_watts_and_reject_wrong_signs_without_predicting_feasibility() {
    let p=plan();
    let result=search::allocate(&p,deadline(),|v| {
        let mut e=fake(&p,v);e.slopes=vec![Some(-100.0),None];Ok(e)
    }).unwrap();
    assert!(result.reason.is_none());assert_eq!(result.newton_trials,0);
    assert!(result.passing.peak<=p.limit);
    assert!(search::allocate(&p,Instant::now(),|v|Ok(fake(&p,v))).is_err());
}

#[test]
fn result_binding_rejects_wrong_peak_watts_cycles_or_ambiguous_gradients() {
    let p=plan();
    let data=J::parse(r#"{"schema":"frankensim.cooling-network.result.v1",
        "objective":{"value_k":301},"transient":{"sampled_peak_objective_k":310},
        "repeated_cycles":{"status":"fixed-count-complete","cycles_completed":2,
          "total_accepted_steps":14,"total_solid_solves":7,"sampled_peak_objective_k":310,
          "sampled_peak_time_s":20,"adjoint":{"method":"discrete-backward-euler-coupled-adjoint",
            "qoi":"sampled-peak","value_k":310,"time_s":20,"cycles":2,
            "component_power_sensitivities":{"intervals":[
              {"interval":0,"rows":[{"component":"chip","applied_power_w":0,"dtemperature_dpower_w_k_per_w":1},
                {"component":"memory","applied_power_w":0,"dtemperature_dpower_w_k_per_w":0.25}]},
              {"interval":1,"rows":[]}]}}}}"#).unwrap();
    assert_eq!(inspect(&p,&[0.0,0.0],data.clone()).unwrap().slopes,vec![Some(1.0),Some(0.25)]);
    for (key,value) in [("value_k",311.0),("time_s",21.0),("cycles",1.0)] {
        let mut changed=data.clone();
        let adjoint=input::member_mut(input::member_mut(&mut changed,"repeated_cycles").unwrap(),"adjoint").unwrap();
        input::put(adjoint,key,number_value(value).unwrap()).unwrap();
        assert!(inspect(&p,&[0.0,0.0],changed).is_err());
    }
    assert!(inspect(&p,&[1.0,0.0],data.clone()).is_err());
    let encoded=serialize(&data).unwrap();
    let ambiguous=J::parse(&encoded.replace("\"component\":\"memory\"","\"component\":\"chip\"")).unwrap();
    assert!(inspect(&p,&[0.0,0.0],ambiguous).is_err());
}

#[test]
fn expired_result_serialization_preserves_the_allocation_but_not_success_status() {
    let doc=J::parse(r#"{"status":"priority-allocation-complete","selected":[3,0],"reason":null}"#).unwrap();
    let (code,text)=publish(doc,true,Instant::now()).unwrap();
    assert_eq!(code,exit::BUDGET);
    let doc=J::parse(&text).unwrap();
    assert_eq!(doc.str_field("status"),Some("budget-exhausted"));
    assert_eq!(doc.get("selected").unwrap().as_array().unwrap().len(),2);
}

#[path = "multi_limit_tests.rs"]
mod multi_limits;
