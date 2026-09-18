//! Exercise the actual allocator and real cooling subprocesses, not fake peaks.
use super::*;
const MULTI: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/allocate-multiple-temperature-limits.json"));
fn phase(result: &J) -> &J { result.get("repeated_cycles").unwrap_or_else(|| result.get("transient").unwrap()) }
fn limit_rows(result: &J) -> &[J] {
    result.path(&["thermal_constraints","rows"]).unwrap().as_array().unwrap()
}
fn limit<'a>(result: &'a J, name: &str) -> &'a J {
    limit_rows(result).iter().find(|r| r.str_field("name")==Some(name)).unwrap()
}
fn verify_limits(dir: &Scratch, result: &J) {
    verify_replay(dir,result);
    assert_eq!(result.path(&["thermal_constraints","all_passed"]),Some(&J::Bool(true)));
    for row in limit_rows(result) {
        assert!(n(row,"sampled_peak_k")<=n(row,"temperature_limit_k"));
        let mut request=result.get("resolved_request").unwrap().clone();
        put(&mut request,"objective",row.get("objective").unwrap().clone());
        put(member(&mut request,"transient"),"temperature_limit_k",row.get("temperature_limit_k").unwrap().clone());
        let replay=dir.cooling(&request);
        assert_eq!(replay.get("solid_temperatures_k"),result.path(&["cooling_result","solid_temperatures_k"]));
        near(n(phase(&replay),"sampled_peak_objective_k"),n(row,"sampled_peak_k"),1e-10);
        near(n(phase(&replay),"sampled_peak_time_s"),n(row,"sampled_peak_time_s"),0.0);
        if let Some(adjoint)=phase(&replay).get("adjoint").filter(|v| v.as_object().is_some()) {
            let intervals=adjoint.path(&["component_power_sensitivities","intervals"]).unwrap().as_array().unwrap();
            let derivatives=row.get("dpeak_dcontrolled_power_w_k_per_w").unwrap().as_array().unwrap();
            for (axis,derivative) in result.get("selected").unwrap().as_array().unwrap().iter().zip(derivatives) {
                let interval=intervals.iter().find(|r| r.get("interval")==axis.get("interval")).unwrap();
                let component=interval.get("rows").unwrap().as_array().unwrap().iter()
                    .find(|r| r.str_field("component")==axis.str_field("component")).unwrap();
                near(derivative.as_f64().unwrap(),n(component,"dtemperature_dpower_w_k_per_w"),1e-10);
            }
        }
    }
}

#[test]
fn a_cooler_memory_limit_restricts_an_otherwise_passing_global_peak() {
    let dir=Scratch::new();let base=J::parse(BASE).unwrap();let spec=J::parse(MULTI).unwrap();
    let result=success(&dir.allocate(&base,&spec));verify_limits(&dir,&result);
    assert_eq!(result.path(&["thermal_constraints","active_constraint"]).unwrap().as_str(),Some("memory-limit"));
    near(selected(&result,"chip"),10.1236293473,0.0011);
    assert_eq!(selected(&result,"standby"),3.0);assert_eq!(selected(&result,"memory"),6.0);
    assert!(n(&result,"selected_sampled_peak_k")<302.0);
    assert!(n(limit(&result,"memory-limit"),"sampled_peak_k")<n(&result,"selected_sampled_peak_k"));
    assert_eq!(n(&result,"trajectories_per_candidate"),3.0);
    assert_eq!(n(&result,"trajectory_evaluations_completed"),3.0*n(&result,"evaluations_completed"));
    assert_eq!(n(&result,"total_completed_trajectory_steps"),42.0*n(&result,"evaluations_completed"));
    let mut only_global=spec.clone();remove(&mut only_global,"thermal_constraints");
    let unchecked=success(&dir.allocate(&base,&only_global));assert_eq!(selected(&unchecked,"chip"),24.0);
    let mut memory=unchecked.get("resolved_request").unwrap().clone();
    put(&mut memory,"objective",limit(&result,"memory-limit").get("objective").unwrap().clone());
    let memory=dir.cooling(&memory);
    assert!(n(phase(&memory),"sampled_peak_objective_k")>301.4);
    let mut reordered=spec;rows(member(&mut reordered,"thermal_constraints")).reverse();
    let reordered=success(&dir.allocate(&base,&reordered));
    assert_eq!(result.get("selected"),reordered.get("selected"));
}

#[test]
fn derivative_free_and_changed_active_constraints_still_require_every_limit() {
    let dir=Scratch::new();let mut base=J::parse(BASE).unwrap();let mut spec=J::parse(MULTI).unwrap();
    remove(member(&mut base,"transient"),"adjoint");
    let no_adjoint=success(&dir.allocate(&base,&spec));verify_limits(&dir,&no_adjoint);
    assert_eq!(n(&no_adjoint,"newton_trials"),0.0);
    assert_eq!(no_adjoint.path(&["cooling_result","repeated_cycles","adjoint"]),Some(&J::Null));
    for row in limit_rows(&no_adjoint) {
        assert!(row.get("dpeak_dcontrolled_power_w_k_per_w").unwrap().as_array().unwrap().iter().all(|v|v==&J::Null));
    }
    near(selected(&no_adjoint,"chip"),10.1236293473,0.0011);
    put(&mut rows(member(&mut spec,"thermal_constraints"))[0],"temperature_limit_k",number(301.5));
    let base=J::parse(BASE).unwrap();let changed=success(&dir.allocate(&base,&spec));verify_limits(&dir,&changed);
    assert_eq!(changed.path(&["thermal_constraints","active_constraint"]).unwrap().as_str(),Some("chip-limit"));
    near(selected(&changed,"chip"),8.259831758,0.0011);
    assert!(selected(&changed,"chip")<selected(&no_adjoint,"chip")-1.0);
}

#[test]
fn component_limits_keep_their_own_peak_times_and_include_initial_wall_maxima() {
    let dir=Scratch::new();let mut base=J::parse(BASE).unwrap();
    let power=member(member(&mut base,"solid"),"component_power");
    for row in rows(member(power,"components")) {
        if row.str_field("name")==Some("memory") {put(row,"watts",number(0.0));}
    }
    put(power,"total_w",number(12.0));
    let mut spec=J::parse(MULTI).unwrap();
    put(&mut spec,"priority",J::parse(r#"[{"component":"standby","interval":0,"min_power_w":0,"max_power_w":1}]"#).unwrap());
    put(&mut spec,"thermal_constraints",J::parse(r#"[
      {"name":"memory","component":"memory","temperature_limit_k":304},
      {"name":"wall","objective":{"mean_wall_region":"last-face"},"temperature_limit_k":301}]"#).unwrap());
    let result=success(&dir.allocate(&base,&spec));verify_limits(&dir,&result);
    assert_eq!(n(limit(&result,"primary"),"sampled_peak_time_s"),20.0);
    assert_eq!(n(limit(&result,"memory"),"sampled_peak_time_s"),28.0);
    assert_eq!(n(limit(&result,"wall"),"sampled_peak_time_s"),0.0);
    assert_eq!(n(limit(&result,"wall"),"sampled_peak_k"),300.0);
    put(&mut rows(member(&mut spec,"thermal_constraints"))[1],"temperature_limit_k",number(299.0));
    let refused=dir.allocate(&base,&spec);assert!(!refused.status.success());assert!(refused.stdout.is_empty());
    assert!(String::from_utf8_lossy(&refused.stderr).contains("wall"));
}

#[test]
fn all_constraints_share_step_limits_and_bad_selectors_cannot_be_ignored() {
    let dir=Scratch::new();let mut base=J::parse(BASE).unwrap();remove(member(&mut base,"transient"),"adjoint");
    let mut spec=J::parse(MULTI).unwrap();put(&mut spec,"max_total_steps",number(14.0));
    let too_small=dir.allocate(&base,&spec);assert_eq!(too_small.status.code(),Some(6));assert!(too_small.stdout.is_empty());
    put(&mut spec,"max_total_steps",number(42.0));
    let partial=dir.allocate(&base,&spec);assert_eq!(partial.status.code(),Some(6));
    let result=parsed(&partial);verify_limits(&dir,&result);
    assert_eq!(result.str_field("status"),Some("budget-exhausted"));
    assert_eq!(n(&result,"evaluations_completed"),1.0);assert_eq!(n(&result,"trajectory_evaluations_completed"),3.0);
    assert_eq!(n(&result,"total_completed_trajectory_steps"),42.0);assert_eq!(selected(&result,"chip"),0.0);
    for constraints in [r#"[]"#,r#"[{"name":"primary","component":"chip","temperature_limit_k":302}]"#,
        r#"[{"name":"bad","component":"missing","temperature_limit_k":302}]"#,
        r#"[{"name":"bad","objective":{"max_vertices":[99999]},"temperature_limit_k":302}]"#] {
        let mut invalid=J::parse(MULTI).unwrap();put(&mut invalid,"thermal_constraints",J::parse(constraints).unwrap());
        let failed=dir.allocate(&base,&invalid);assert!(!failed.status.success());assert!(failed.stdout.is_empty());
    }
}

#[test]
fn enclosure_limits_replay_reflection_and_use_each_observations_own_adjoint() {
    let dir=Scratch::new();let base=J::parse(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/adjoint-enclosure-pulse.json"))).unwrap();
    let mut spec=J::parse(MULTI).unwrap();put(&mut spec,"temperature_limit_k",number(320.0));
    put(&mut spec,"priority",J::parse(r#"[{"component":"heater","interval":0,"min_power_w":1,"max_power_w":2}]"#).unwrap());
    put(&mut spec,"thermal_constraints",J::parse(r#"[{"name":"receiver","objective":{"max_wall_region":"receiver"},"temperature_limit_k":320}]"#).unwrap());
    let result=success(&dir.allocate(&base,&spec));verify_limits(&dir,&result);
    assert_eq!(selected(&result,"heater"),2.0);assert_eq!(n(&result,"trajectory_evaluations_completed"),4.0);
    assert_eq!(result.path(&["cooling_result","radiation","model"]).unwrap().as_str(),Some("closed-gray-diffuse-enclosure"));
    near(n(result.path(&["cooling_result","radiation"]).unwrap(),"radiative_out_w"),0.0,1e-7);
}
