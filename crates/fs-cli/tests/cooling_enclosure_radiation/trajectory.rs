//! Complete-trajectory consumers of the shared enclosure response transpose.
use super::*;
const PULSE:&str=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/enclosure-radiation-pulse.json"));
fn pulse(qoi:&str)->J {
    let mut input=J::parse(PULSE).unwrap();
    put(member(&mut input,"transient"),"adjoint",J::parse(&format!(
        r#"{{"qoi":"{qoi}","max_checkpoint_bytes":1048576,"component_power":true}}"#)).unwrap());
    input
}
fn phase(result:&J)->&J {result.get("repeated_cycles").unwrap_or_else(||result.get("transient").unwrap())}
fn adjoint(result:&J)->&J {phase(result).get("adjoint").unwrap()}
fn objective(result:&J,qoi:&str)->f64 {
    if qoi=="final" {value(result)} else {n(phase(result),"sampled_peak_objective_k")}
}
fn trajectory_primal(input:&J)->J {
    let mut input=primal(input);remove(member(&mut input,"transient"),"adjoint");input
}
fn interval(input:&mut J,index:usize)->&mut J {&mut rows(member(member(input,"transient"),"intervals"))[index]}
fn interval_gradient(result:&J,index:usize,key:&str)->f64 {
    n(&adjoint(result).get("intervals").unwrap().as_array().unwrap()[index],key)
}
fn finish_gradient(result:&J,name:&str)->f64 {
    let row=adjoint(result).path(&["radiation","surfaces"]).unwrap().as_array().unwrap().iter()
        .find(|r|r.str_field("surface")==Some(name)).unwrap();
    assert!(row.get("dtemperature_dambient_temperature").is_none());
    n(row,"dtemperature_dlog_emissivity_k")
}
fn perturb(input:&mut J,control:usize,amount:f64) {
    match control {
        0|1=>{
            let row=&mut rows(member(enclosure(input),"surfaces"))[control];
            let epsilon=n(row,"emissivity");put(row,"emissivity",num(epsilon*amount.exp()));
        }
        2=>{let row=interval(input,0);let p=n(row,"power_scale");put(row,"power_scale",num(p*(1.0+amount)));}
        3|4=>{let row=interval(input,control-3);let s=n(row,"fan_speed_ratio");put(row,"fan_speed_ratio",num(s*amount.exp()));}
        5=>{
            for v in rows(member(member(input,"transient"),"element_heat_capacities_j_m3_k")) {
                *v=num(v.as_f64().unwrap()*(1.0+amount));
            }
        }
        6=>{let s=member(input,"transient");let t=n(s,"initial_temperature_k");put(s,"initial_temperature_k",num(t+amount));}
        _=>{let fan=member(member(input,"hydraulics"),"fan");let t=n(fan,"temperature_k");put(fan,"temperature_k",num(t+amount));}
    }
}

#[test]
fn repeated_enclosure_adjoint_matches_complete_trajectory_perturbations_without_changing_forward_history() {
    for qoi in ["final","sampled-peak"] {
        let input=pulse(qoi);let result=run(&input);let plain=run(&trajectory_primal(&input));
        assert_eq!(result.get("solid_temperatures_k"),plain.get("solid_temperatures_k"));
        assert_eq!(result.get("radiation"),plain.get("radiation"));
        assert_eq!(phase(&result).get("cycles"),phase(&plain).get("cycles"));
        for key in ["stored_energy_change_j","input_energy_j","air_energy_gain_j","radiative_energy_loss_j","forward_solid_solves"] {
            assert_eq!(phase(&result).get(key),phase(&plain).get(key));
        }
        assert_eq!(result.path(&["transient","adjoint"]),Some(&J::Null));
        assert!(n(adjoint(&result),"reconstruction_solid_solves")>0.0);
        let derivatives=[finish_gradient(&result,"emitter"),finish_gradient(&result,"receiver"),
            interval_gradient(&result,0,"dtemperature_dpower_multiplier_k"),
            interval_gradient(&result,0,"dtemperature_dlog_fan_speed_ratio_k"),
            interval_gradient(&result,1,"dtemperature_dlog_fan_speed_ratio_k"),
            n(adjoint(&result),"dtemperature_dcapacity_multiplier_k"),
            n(adjoint(&result),"dtemperature_duniform_initial_k"),
            adjoint(&result).get("dtemperature_dinlet_temperatures").unwrap().as_array().unwrap()[0].as_f64().unwrap()];
        for (control,&gradient) in derivatives.iter().enumerate() {
            let h=1e-4;let mut plus=trajectory_primal(&input);let mut minus=plus.clone();
            perturb(&mut plus,control,h);perturb(&mut minus,control,-h);
            near(gradient,(objective(&run(&plus),qoi)-objective(&run(&minus),qoi))/(2.0*h),5e-5);
        }
        let component=&adjoint(&result).path(&["component_power_sensitivities","intervals"]).unwrap()
            .as_array().unwrap()[0].get("rows").unwrap().as_array().unwrap()[0];
        near(5.0*n(component,"dtemperature_dpower_w_k_per_w"),derivatives[2],1e-9);
    }
}

#[test]
fn repeated_enclosure_controls_equal_unrolled_history_and_later_occurrences_do_not_reach_backwards() {
    let input=pulse("sampled-peak");let result=run(&input);
    let mut unrolled=input.clone();let s=member(&mut unrolled,"transient");remove(s,"repeat");
    let phases=rows(member(s,"intervals"));phases.extend(phases.clone());
    let unrolled=run(&unrolled);
    assert_eq!(result.get("solid_temperatures_k"),unrolled.get("solid_temperatures_k"));
    near(objective(&result,"sampled-peak"),objective(&unrolled,"sampled-peak"),1e-10);
    for name in ["emitter","receiver"] {near(finish_gradient(&result,name),finish_gradient(&unrolled,name),1e-8);}
    for i in 0..2 {for key in ["dtemperature_dpower_multiplier_k","dtemperature_dlog_fan_speed_ratio_k"] {
        near(interval_gradient(&result,i,key),interval_gradient(&unrolled,i,key)+interval_gradient(&unrolled,i+2,key),1e-8);
    }}
    assert_eq!(interval_gradient(&unrolled,3,"dtemperature_dlog_fan_speed_ratio_k"),0.0);
    assert!(interval_gradient(&result,1,"dtemperature_dlog_fan_speed_ratio_k").abs()>1e-5);
    let mut single=input;remove(member(&mut single,"transient"),"repeat");
    let single=run(&single);
    assert_eq!(interval_gradient(&single,1,"dtemperature_dlog_fan_speed_ratio_k"),0.0);
}

#[test]
fn initial_enclosure_peak_has_zero_finish_controls_and_memory_or_reverse_exhaustion_refuses() {
    let mut input=pulse("sampled-peak");
    put(&mut input,"objective",J::parse(r#"{"mean_wall_region":"emitter","gradient":false}"#).unwrap());
    put(member(&mut input,"transient"),"initial_temperature_k",num(350.0));
    for row in rows(member(member(&mut input,"transient"),"intervals")) {put(row,"power_scale",num(0.0));}
    let result=run(&input);
    assert_eq!(n(adjoint(&result),"state_index"),0.0);
    assert_eq!(n(adjoint(&result),"reconstructed_solid_endpoints"),0.0);
    for name in ["emitter","receiver"] {assert_eq!(finish_gradient(&result,name),0.0);}
    for mode in 0..2 {
        let mut rejected=pulse("final");
        if mode==0 {put(member(member(&mut rejected,"transient"),"adjoint"),"max_checkpoint_bytes",num(1.0));}
        else {put(member(&mut rejected,"budgets"),"derivative_iterations",num(1.0));}
        let result=output(&rejected);assert_eq!(result.status.code(),Some(6));assert!(result.stdout.is_empty());
    }
}

#[test]
fn enclosure_workload_sizing_retains_the_selected_total_trajectory_adjoint() {
    let mut input=pulse("sampled-peak");let nominal=run(&trajectory_primal(&input));
    let s=member(&mut input,"transient");put(s,"temperature_limit_k",num(objective(&nominal,"sampled-peak")));
    put(s,"power_design",J::parse(r#"{"min_power_multiplier":0.8,"max_power_multiplier":1.2,"power_multiplier_tolerance":0.5,"temperature_tolerance_k":10,"max_evaluations":3}"#).unwrap());
    let result=run(&input);let decision=result.get("transient_power_design").unwrap();
    assert_eq!(decision.str_field("search_method"),Some("safeguarded-adjoint-newton-bisection"));
    assert_eq!(decision.str_field("status"),Some("target-bracketed"));
    let scale=n(decision,"selected_power_multiplier");near(scale,0.8,1e-15);
    remove(member(&mut input,"transient"),"power_design");
    for row in rows(member(member(&mut input,"transient"),"intervals")) {
        let p=n(row,"power_scale");put(row,"power_scale",num(scale*p));
    }
    let replay=run(&input);
    assert_eq!(result.get("solid_temperatures_k"),replay.get("solid_temperatures_k"));
    assert_eq!(adjoint(&result),adjoint(&replay));
    assert!(objective(&result,"sampled-peak")<=objective(&nominal,"sampled-peak"));
}
