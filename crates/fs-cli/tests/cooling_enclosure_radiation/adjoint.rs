//! Actual command consumers: complete perturbed solves, search and adaptation.
use super::*;
fn value(result:&J)->f64 {n(result.get("objective").unwrap(),"value_k")}
fn differentiated()->J {
    let mut input=J::parse(BASE).unwrap();
    put(&mut input,"objective",J::parse(r#"{"mean_wall_region":"emitter","gradient":true}"#).unwrap());
    input
}
fn primal(input:&J)->J {
    let mut result=input.clone();put(member(&mut result,"objective"),"gradient",J::Bool(false));result
}
fn finish_gradient(result:&J,name:&str)->f64 {
    let row=result.path(&["radiation","adjoint","surfaces"]).unwrap().as_array().unwrap().iter()
        .find(|r|r.str_field("surface")==Some(name)).unwrap();
    n(row,"dobjective_dlog_emissivity_k")
}
fn surface<'a>(input:&'a mut J,name:&str)->&'a mut J {
    rows(member(member(input,"solid"),"surfaces")).iter_mut()
        .find(|s|s.str_field("name")==Some(name)).unwrap()
}
fn nonlinear(input:&mut J) {
    for row in rows(member(member(input,"solid"),"materials")) {
        remove(row,"conductivity_w_m_k");
        put(row,"conductivity_curve",J::parse(r#"{"temperature_k":[250,450],"conductivity_w_m_k":[0.5,4.5]}"#).unwrap());
    }
}

#[test]
fn enclosure_finish_inlet_coefficient_and_speed_gradients_use_complete_physics() {
    for nonlinear_material in [false,true] {
        let mut input=differentiated();if nonlinear_material {nonlinear(&mut input);}
        let result=run(&input);let plain=run(&primal(&input));
        assert_eq!(result.get("solid_temperatures_k"),plain.get("solid_temperatures_k"));
        let rad=result.get("radiation").unwrap();let old=plain.get("radiation").unwrap();
        assert_eq!(rad.get("surfaces"),old.get("surfaces"));
        assert_eq!(rad.get("forward_solid_solves"),old.get("forward_solid_solves"));
        assert_eq!(n(rad,"total_solid_solves"),n(old,"total_solid_solves")+1.0);
        let h=1e-4_f64;
        for (name,e) in [("emitter",0.8),("receiver",0.6)] {
            let mut plus=primal(&input);let mut minus=primal(&input);
            for (trial,factor) in [(&mut plus,h.exp()),(&mut minus,(-h).exp())] {
                let row=rows(member(enclosure(trial),"surfaces")).iter_mut()
                    .find(|r|r.str_field("surface")==Some(name)).unwrap();
                put(row,"emissivity",num(e*factor));
            }
            near(finish_gradient(&result,name),(value(&run(&plus))-value(&run(&minus)))/(2.0*h),3e-5);
            let mut plus=primal(&input);let mut minus=primal(&input);
            let coefficient=n(surface(&mut plus,name),"htc_w_m2_k");
            put(surface(&mut plus,name),"htc_w_m2_k",num(coefficient*h.exp()));
            put(surface(&mut minus,name),"htc_w_m2_k",num(coefficient*(-h).exp()));
            let wall=result.get("walls").unwrap().as_array().unwrap().iter()
                .find(|r|r.str_field("region")==Some(name)).unwrap();
            near(n(wall,"dobjective_dlog_htc"),(value(&run(&plus))-value(&run(&minus)))/(2.0*h),3e-5);
        }
        let mut plus=primal(&input);let mut minus=primal(&input);
        put(member(member(&mut plus,"hydraulics"),"fan"),"temperature_k",num(300.001));
        put(member(member(&mut minus,"hydraulics"),"fan"),"temperature_k",num(299.999));
        let inlet=result.get("dobjective_dinlet_k").unwrap().as_array().unwrap()[0].as_f64().unwrap();
        near(inlet,(value(&run(&plus))-value(&run(&minus)))/0.002,3e-5);
        let mut plus=primal(&input);let mut minus=primal(&input);
        put(member(member(&mut plus,"hydraulics"),"fan"),"speed_ratio",num(h.exp()));
        put(member(member(&mut minus,"hydraulics"),"fan"),"speed_ratio",num((-h).exp()));
        near(n(result.get("fan_speed_sensitivity").unwrap(),"dobjective_dlog_speed_ratio_k"),
            (value(&run(&plus))-value(&run(&minus)))/(2.0*h),3e-5);
        near(n(rad,"radiative_out_w"),0.0,1e-7);
    }
}

#[test]
fn black_finish_derivatives_are_one_sided_and_matrix_axis_order_replays() {
    let mut input=differentiated();finish(&mut input,1.0,1.0);
    let original=output(&input);let result=success(&original);
    let mut minus=primal(&input);let h=1e-5_f64;finish(&mut minus,(-h).exp(),1.0);
    near(finish_gradient(&result,"emitter"),(value(&result)-value(&run(&minus)))/h,1e-3);
    let mut reordered=input;
    rows(member(enclosure(&mut reordered),"surfaces")).reverse();
    let matrix=rows(member(enclosure(&mut reordered),"view_factors"));
    matrix.reverse();for row in matrix {rows(row).reverse();}
    let replay=output(&reordered);success(&replay);assert_eq!(original.stdout,replay.stdout);
}

#[test]
fn nonmatching_contact_receives_the_total_reflected_radiation_adjoint() {
    let mut input=J::parse(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/nonmatching-contact-hotspot.json"))).unwrap();
    let mut radiation=J::parse(BASE).unwrap().get("radiation").unwrap().clone();
    let patches=rows(member(member(&mut radiation,"enclosure"),"surfaces"));
    put(&mut patches[0],"surface",J::Str("first-face".into()));
    put(&mut patches[1],"surface",J::Str("last-face".into()));
    put(&mut input,"radiation",radiation);nonlinear(&mut input);
    let result=run(&input);
    let row=&result.path(&["contact_sensitivities","rows"]).unwrap().as_array().unwrap()[0];
    let h=1e-4_f64;let mut plus=primal(&input);let mut minus=primal(&input);
    for (trial,factor) in [(&mut plus,h.exp()),(&mut minus,(-h).exp())] {
        let row=&mut rows(member(member(trial,"solid"),"contacts"))[0];
        put(row,"resistance_m2_k_w",num(0.01*factor));
    }
    near(n(row,"dobjective_dlog_resistance_k"),(value(&run(&plus))-value(&run(&minus)))/(2.0*h),3e-5);
    assert_eq!(result.get("solid_temperatures_k"),run(&primal(&input)).get("solid_temperatures_k"));
}

#[test]
fn adjoint_guided_enclosure_fan_sizing_replays_the_selected_complete_candidate() {
    let mut input=differentiated();
    let mut low=primal(&input);let mut high=low.clone();
    put(member(member(&mut low,"hydraulics"),"fan"),"speed_ratio",num(0.6));
    put(member(member(&mut high,"hydraulics"),"fan"),"speed_ratio",num(1.8));
    let lo=value(&run(&low));let hi=value(&run(&high));assert!(lo>hi);
    let target=0.5*(lo+hi);
    put(&mut input,"fan_speed_design",J::parse(&format!(r#"{{"min_speed_ratio":0.6,"max_speed_ratio":1.8,"temperature_limit_k":{target},"speed_ratio_tolerance":0.001,"temperature_tolerance_k":0.001,"max_evaluations":64}}"#)).unwrap());
    let result=run(&input);let design=result.get("fan_speed_design").unwrap();
    assert!(value(&result)<=target);assert!(n(design,"newton_trials")>0.0);
    assert!(n(design,"speed_bracket_width")<=0.001);
    let selected=n(design,"selected_speed_ratio");
    remove(&mut input,"fan_speed_design");
    put(member(member(&mut input,"hydraulics"),"fan"),"speed_ratio",num(selected));
    let replay=run(&input);
    assert_eq!(result.get("solid_temperatures_k"),replay.get("solid_temperatures_k"));
    assert_eq!(result.get("fan_speed_sensitivity"),replay.get("fan_speed_sensitivity"));
    assert_eq!(result.path(&["radiation","adjoint"]),replay.path(&["radiation","adjoint"]));
}

#[test]
fn enclosure_goal_marking_requires_global_confirmation_and_preserves_final_physics() {
    let mut input=J::parse(BASE).unwrap();
    put(&mut input,"mesh_convergence",J::parse(r#"{"strategy":"goal-recovery","marking_fraction":0.5,"max_refinements":3,"consecutive_passes":2,"temperature_tolerance_k":2,"max_vertices":10000,"max_tetrahedra":10000}"#).unwrap());
    let result=run(&input);let mesh=result.get("mesh_convergence").unwrap();
    assert_eq!(mesh.get("global_confirmation"),Some(&J::Bool(true)));
    assert!(n(mesh,"total_adjoint_sweeps")>0.0);
    assert_eq!(result.path(&["radiation","adjoint"]),Some(&J::Null));
    assert_eq!(result.get("contact_sensitivities"),Some(&J::Null));
    let resolved=mesh.get("resolved_request").unwrap();
    assert_eq!(resolved.get("radiation"),input.get("radiation"));
    let replay=run(resolved);
    assert_eq!(result.get("solid_temperatures_k"),replay.get("solid_temperatures_k"));
    assert_eq!(result.path(&["radiation","surfaces"]),replay.path(&["radiation","surfaces"]));
    let history=mesh.get("history").unwrap().as_array().unwrap();
    assert_eq!(history.last().unwrap().str_field("arrived_by"),Some("uniform"));
    let mut insufficient=input;put(member(&mut insufficient,"mesh_convergence"),"max_refinements",num(2.0));
    let attempted=output(&insufficient);
    if attempted.status.success() {
        let doc=success(&attempted);let h=doc.path(&["mesh_convergence","history"]).unwrap().as_array().unwrap();
        assert_eq!(h.last().unwrap().str_field("arrived_by"),Some("uniform"));
    } else {assert!(attempted.stdout.is_empty());}
}

#[test]
fn derivative_budgets_and_adaptive_trajectory_derivatives_publish_no_partial_gradient() {
    let mut input=differentiated();put(member(&mut input,"budgets"),"derivative_iterations",num(1.0));
    let result=output(&input);assert_eq!(result.status.code(),Some(6));assert!(result.stdout.is_empty());
    assert!(String::from_utf8_lossy(&result.stderr).contains("adjoint"));
    run(&primal(&input));
    let mut transient=primal(&input);
    put(&mut transient,"transient",J::parse(r#"{"initial_temperature_k":300,"volumetric_heat_capacity_j_m3_k":1000000,"max_step_s":1,"max_steps":2,"intervals":[{"duration_s":1,"power_scale":1,"fan_speed_ratio":1}],"adjoint":{"qoi":"final","max_checkpoint_bytes":1048576},"adaptive":{"absolute_tolerance_k":0.01,"relative_tolerance":0,"minimum_trial_step_s":0.001,"max_trials":10}}"#).unwrap());
    let result=output(&transient);assert!(!result.status.success());assert!(result.stdout.is_empty());
    assert!(String::from_utf8_lossy(&result.stderr).contains("fixed timesteps"));
}

#[path="trajectory.rs"]
mod trajectory;
