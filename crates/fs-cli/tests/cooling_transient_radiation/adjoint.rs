//! Whole-command derivatives, not differences of an isolated boundary law.
//! Every finite-difference side reruns the complete cooling trajectory.
use super::*;

fn input(qoi: &str, cycles: usize) -> J {
    let mut root = short();
    let schedule = member(&mut root,"transient");
    put(schedule,"adjoint",J::parse(&format!(
        r#"{{"qoi":"{qoi}","max_checkpoint_bytes":1048576}}"#)).unwrap());
    if cycles > 1 {
        put(schedule,"repeat",J::parse(&format!(
            r#"{{"cycles":{cycles},"max_total_steps":100}}"#)).unwrap());
    }
    root
}
fn trajectory(result: &J) -> &J {
    result.get("repeated_cycles").unwrap_or_else(|| result.get("transient").unwrap())
}
fn gradient(result: &J) -> &J { trajectory(result).get("adjoint").unwrap() }
fn patch<'a>(result: &'a J, name: &str) -> &'a J {
    gradient(result).path(&["radiation","surfaces"]).unwrap().as_array().unwrap()
        .iter().find(|row| row.str_field("surface") == Some(name)).unwrap()
}
fn patch_input<'a>(root: &'a mut J, name: &str) -> &'a mut J {
    let J::Array(rows) = member(member(root,"radiation"),"surfaces") else { panic!() };
    rows.iter_mut().find(|row| row.str_field("surface") == Some(name)).unwrap()
}
fn interval(root: &mut J, index: usize) -> &mut J {
    let J::Array(rows) = member(member(root,"transient"),"intervals") else { panic!() };
    &mut rows[index]
}
fn observable(result: &J, peak: bool) -> f64 {
    if peak { n(trajectory(result),"sampled_peak_objective_k") }
    else { n(result.get("objective").unwrap(),"value_k") }
}
fn difference(root: &J, peak: bool, width: f64, mut perturb: impl FnMut(&mut J,f64)) -> f64 {
    let mut values = [0.0;2];
    for (slot, delta) in values.iter_mut().zip([-width,width]) {
        let mut trial = root.clone();
        remove(member(&mut trial,"transient"),"adjoint");
        perturb(&mut trial,delta);
        *slot = observable(&run(&trial),peak);
    }
    (values[1]-values[0])/(2.0*width)
}
fn derivative_close(actual: f64, reference: f64) {
    near(actual,reference,3e-5+2e-4*reference.abs());
}

#[test]
fn complete_nonlinear_contact_trajectories_match_radiation_and_storage_derivatives() {
    for peak in [false,true] {
        let root = input(if peak {"sampled-peak"} else {"final"},1);
        let result = run(&root);
        let g = gradient(&result);
        near(n(g,"value_k"),observable(&result,peak),1e-10);
        for name in ["first-face","last-face"] {
            let epsilon = n(patch_input(&mut root.clone(),name),"emissivity");
            let ambient = n(patch_input(&mut root.clone(),name),"ambient_temperature_k");
            let numeric = difference(&root,peak,0.002,|trial,d|
                put(patch_input(trial,name),"emissivity",num(epsilon*d.exp())));
            derivative_close(n(patch(&result,name),"dtemperature_dlog_emissivity_k"),numeric);
            derivative_close(n(patch(&result,name),"dtemperature_demissivity_k")*epsilon,numeric);
            let numeric = difference(&root,peak,0.01,|trial,d|
                put(patch_input(trial,name),"ambient_temperature_k",num(ambient+d)));
            derivative_close(n(patch(&result,name),"dtemperature_dambient_temperature"),numeric);
        }
        let numeric = difference(&root,peak,0.002,|trial,d|
            put(member(trial,"transient"),"initial_temperature_k",num(300.0+d)));
        derivative_close(n(g,"dtemperature_duniform_initial_k"),numeric);
        let numeric = difference(&root,peak,0.002,|trial,d| {
            let J::Array(values) = member(member(trial,"transient"),"element_heat_capacities_j_m3_k") else { panic!() };
            for value in values { *value = num(value.as_f64().unwrap()*(1.0+d)); }
        });
        derivative_close(n(g,"dtemperature_dcapacity_multiplier_k"),numeric);
        let rows = g.get("intervals").unwrap().as_array().unwrap();
        let numeric = difference(&root,peak,0.002,|trial,d|
            put(interval(trial,0),"power_scale",num(1.0+d)));
        derivative_close(n(&rows[0],"dtemperature_dpower_multiplier_k"),numeric);
        let numeric = difference(&root,peak,0.002,|trial,d|
            put(interval(trial,0),"fan_speed_ratio",num(d.exp())));
        derivative_close(n(&rows[0],"dtemperature_dlog_fan_speed_ratio_k"),numeric);
    }
}

#[test]
fn requesting_the_adjoint_keeps_forward_fields_and_heat_history_unchanged() {
    let root = input("sampled-peak",1);
    let mut plain = root.clone(); remove(member(&mut plain,"transient"),"adjoint");
    let forward = run(&plain); let result = run(&root); let g = gradient(&result);
    assert_eq!(result.get("solid_temperatures_k"),forward.get("solid_temperatures_k"));
    assert_eq!(result.path(&["transient","history"]),forward.path(&["transient","history"]));
    assert_eq!(result.get("radiation"),forward.get("radiation"));
    let t = trajectory(&result);
    near(n(t,"total_solid_solves"),n(t,"forward_solid_solves")+n(g,"reconstruction_solid_solves"),0.0);
    assert!(n(g,"reconstruction_solid_solves")>n(g,"reconstructed_solid_endpoints"));
    assert_eq!(n(t,"forward_solid_solves"),n(trajectory(&forward),"forward_solid_solves"));
    assert!(n(g,"time_s")<=4.0,"fixture peak must precede cooldown control");
    assert_eq!(n(&g.get("intervals").unwrap().as_array().unwrap()[1],"dtemperature_dlog_fan_speed_ratio_k"),0.0);
    energy(&result,"transient");
}

#[test]
fn repeated_and_unrolled_radiative_adjoints_share_the_same_chronological_history() {
    let repeated = input("final",2);
    let mut unrolled = repeated.clone(); remove(member(&mut unrolled,"transient"),"repeat");
    let J::Array(rows) = member(member(&mut unrolled,"transient"),"intervals") else { panic!() };
    rows.extend(rows.clone());
    let a = run(&repeated); let b = run(&unrolled);
    assert_eq!(a.get("solid_temperatures_k"),b.get("solid_temperatures_k"));
    assert_eq!(a.path(&["transient","adjoint"]),Some(&J::Null));
    let ga = gradient(&a); let gb = gradient(&b);
    for field in ["value_k","dtemperature_duniform_initial_k","dtemperature_dcapacity_multiplier_k"] {
        derivative_close(n(ga,field),n(gb,field));
    }
    for name in ["first-face","last-face"] {
        for field in ["dtemperature_dlog_emissivity_k","dtemperature_dambient_temperature"] {
            derivative_close(n(patch(&a,name),field),n(patch(&b,name),field));
        }
    }
    let pa = ga.get("intervals").unwrap().as_array().unwrap();
    let pb = gb.get("intervals").unwrap().as_array().unwrap();
    for i in 0..2 {
        for field in ["dtemperature_dpower_multiplier_k","dtemperature_dlog_fan_speed_ratio_k"] {
            derivative_close(n(&pa[i],field),n(&pb[i],field)+n(&pb[i+2],field));
        }
    }
    let numeric = difference(&repeated,false,0.01,|trial,d|
        put(patch_input(trial,"first-face"),"ambient_temperature_k",num(280.0+d)));
    derivative_close(n(patch(&a,"first-face"),"dtemperature_dambient_temperature"),numeric);
    energy(&a,"repeated_cycles");
}

#[test]
fn an_initial_maximum_has_no_radiation_or_future_workload_sensitivity() {
    let mut root = input("sampled-peak",2);
    // A consistent capacity matrix does not assert an interior-node maximum
    // principle. The directly cooled patch mean is verified to decrease here.
    put(&mut root,"objective",J::parse(r#"{"mean_wall_region":"first-face","gradient":false}"#).unwrap());
    put(member(&mut root,"transient"),"initial_temperature_k",num(350.0));
    for i in 0..2 { put(interval(&mut root,i),"power_scale",num(0.0)); }
    let result = run(&root); let g = gradient(&result);
    assert_eq!(n(g,"state_index"),0.0);
    assert_eq!(n(g,"reconstruction_solid_solves"),0.0);
    near(n(g,"dtemperature_duniform_initial_k"),1.0,1e-14);
    for name in ["first-face","last-face"] {
        assert_eq!(n(patch(&result,name),"dtemperature_dlog_emissivity_k"),0.0);
        assert_eq!(n(patch(&result,name),"dtemperature_dambient_temperature"),0.0);
    }
}

#[test]
fn hot_surroundings_and_reordered_mean_objectives_use_the_actual_surface() {
    let mut root = input("final",1);
    put(&mut root,"objective",J::parse(r#"{"mean_wall_region":"last-face","gradient":false}"#).unwrap());
    put(patch_input(&mut root,"last-face"),"ambient_temperature_k",num(350.0));
    let result = run(&root);
    let numeric = difference(&root,false,0.01,|trial,d|
        put(patch_input(trial,"last-face"),"ambient_temperature_k",num(350.0+d)));
    derivative_close(n(patch(&result,"last-face"),"dtemperature_dambient_temperature"),numeric);
    let J::Array(surfaces) = member(member(&mut root,"solid"),"surfaces") else { panic!() };
    surfaces.reverse();
    let reordered = run(&root);
    near(n(gradient(&result),"value_k"),n(gradient(&reordered),"value_k"),1e-10);
    derivative_close(n(patch(&result,"last-face"),"dtemperature_dambient_temperature"),
        n(patch(&reordered,"last-face"),"dtemperature_dambient_temperature"));
}

#[test]
fn reverse_budgets_and_unsupported_time_policies_do_not_publish_partial_gradients() {
    let root = input("final",1);
    let mut memory = root.clone();
    put(member(member(&mut memory,"transient"),"adjoint"),"max_checkpoint_bytes",num(1.0));
    let result = output(&memory); assert_eq!(result.status.code(),Some(6)); assert!(result.stdout.is_empty());
    let mut derivative = root.clone();
    put(member(&mut derivative,"budgets"),"derivative_iterations",num(1.0));
    let result = output(&derivative); assert!(!result.status.success()); assert!(result.stdout.is_empty());
    let mut adaptive = root.clone();
    put(member(&mut adaptive,"transient"),"adaptive",J::parse(r#"{"absolute_tolerance_k":0.01,"relative_tolerance":0,"minimum_trial_step_s":0.001,"max_trials":100}"#).unwrap());
    let result = output(&adaptive); assert!(!result.status.success()); assert!(result.stdout.is_empty());
    assert!(String::from_utf8_lossy(&result.stderr).contains("fixed timesteps"));
}

#[test]
fn gradient_guided_workload_sizing_replays_the_selected_radiating_candidate() {
    let mut root = input("sampled-peak",1);
    let nominal = run(&root);
    let limit = 300.0+0.7*(n(trajectory(&nominal),"sampled_peak_objective_k")-300.0);
    put(member(&mut root,"transient"),"temperature_limit_k",num(limit));
    put(member(&mut root,"transient"),"power_design",J::parse(r#"{"min_power_multiplier":0,"max_power_multiplier":1,"power_multiplier_tolerance":0.002,"temperature_tolerance_k":0.01,"max_evaluations":40}"#).unwrap());
    let designed = run(&root);
    let report = designed.get("transient_power_design").unwrap();
    assert_eq!(report.str_field("search_method"),Some("safeguarded-adjoint-newton-bisection"));
    let multiplier = n(report,"selected_power_multiplier");
    assert!(multiplier>0.0 && multiplier<1.0);
    assert!(n(report,"sampled_peak_objective_k")<=limit);
    remove(member(&mut root,"transient"),"power_design");
    put(interval(&mut root,0),"power_scale",num(multiplier));
    let replay = run(&root);
    assert_eq!(designed.get("solid_temperatures_k"),replay.get("solid_temperatures_k"));
    assert_eq!(gradient(&designed),gradient(&replay));
    energy(&designed,"transient");
}
