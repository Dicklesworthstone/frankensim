//! Complete-command checks: every difference reruns the nonlinear radiation,
//! air, and solid equations. No frozen reference or effective-h surrogate.
use super::*;

const CORRELATED: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/fan-correlated-hotspot.json"));

fn derivative_input(base: &str) -> J {
    let mut root=J::parse(base).unwrap();
    put(member(&mut root,"objective"),"gradient",J::Bool(true));
    put(member(&mut root,"tolerances"),"temperature_k",number(1e-10));
    put(member(&mut root,"radiation"),"temperature_tolerance_k",number(1e-11));
    root
}
fn patch_mut<'a>(root:&'a mut J,name:&str)->&'a mut J {
    let J::Array(rows)=member(member(root,"radiation"),"surfaces") else {panic!()};
    rows.iter_mut().find(|r|r.str_field("surface")==Some(name)).unwrap()
}
fn patch_gradient<'a>(result:&'a J,name:&str)->&'a J {
    result.path(&["radiation","adjoint","surfaces"]).unwrap().as_array().unwrap().iter()
        .find(|r|r.str_field("surface")==Some(name)).unwrap()
}
fn objective_value(input:&J)->f64 {
    let mut input=input.clone();put(member(&mut input,"objective"),"gradient",J::Bool(false));
    n(run(&input).get("objective").unwrap(),"value_k")
}
fn central(input:&J,step:f64,change:impl Fn(&mut J,f64))->f64 {
    let mut plus=input.clone();let mut minus=input.clone();
    change(&mut plus,step);change(&mut minus,-step);
    (objective_value(&plus)-objective_value(&minus))/(2.0*step)
}
fn checked_gradient(actual:f64,expected:f64) {
    near(actual,expected,2e-5*(1.0+expected.abs()));
}
fn wall<'a>(result:&'a J,name:&str)->&'a J {
    result.get("walls").unwrap().as_array().unwrap().iter()
        .find(|r|r.str_field("region")==Some(name)).unwrap()
}

#[test]
fn emissivity_and_reservoir_gradients_close_the_complete_nonlinear_problem() {
    let input=derivative_input(CONTACT);let result=run(&input);
    for name in ["first-face","last-face"] {
        let gradient=patch_gradient(&result,name);
        let epsilon=central(&input,1e-4,|root,delta| {
            let patch=patch_mut(root,name);let value=n(patch,"emissivity");
            put(patch,"emissivity",number(value*delta.exp()));
        });
        checked_gradient(n(gradient,"dobjective_dlog_emissivity_k"),epsilon);
        let ambient=central(&input,0.002,|root,delta| {
            let patch=patch_mut(root,name);let value=n(patch,"ambient_temperature_k");
            put(patch,"ambient_temperature_k",number(value+delta));
        });
        checked_gradient(n(gradient,"dobjective_dambient_temperature"),ambient);
        assert!(ambient>0.0,"warmer surroundings must influence this fixture");
    }
    let adjoint=result.path(&["radiation","adjoint"]).unwrap();
    assert!(n(adjoint,"equation_residual")<=n(adjoint,"equation_threshold"));
    assert_eq!(adjoint.f64_field("reconstruction_solid_solves"),Some(1.0));
}

#[test]
fn thermal_and_contact_controls_use_the_total_radiation_adjoint() {
    let input=derivative_input(CONTACT);let result=run(&input);
    let inlet=central(&input,0.002,|root,delta| {
        let fan=member(member(root,"hydraulics"),"fan");let value=n(fan,"temperature_k");
        put(fan,"temperature_k",number(value+delta));
    });
    checked_gradient(result.get("dobjective_dinlet_k").unwrap().as_array().unwrap()[0].as_f64().unwrap(),inlet);
    for name in ["first-face","last-face"] {
        let expected=central(&input,1e-4,|root,delta| {
            let J::Array(rows)=member(member(root,"solid"),"surfaces") else {panic!()};
            let row=rows.iter_mut().find(|s|s.str_field("name")==Some(name)).unwrap();
            let value=n(row,"htc_w_m2_k");put(row,"htc_w_m2_k",number(value*delta.exp()));
        });
        checked_gradient(n(wall(&result,name),"dobjective_dlog_htc"),expected);
    }
    let contact=central(&input,1e-4,|root,delta| {
        let J::Array(rows)=member(member(root,"solid"),"contacts") else {panic!()};
        let value=n(&rows[0],"resistance_m2_k_w");
        put(&mut rows[0],"resistance_m2_k_w",number(value*delta.exp()));
    });
    let row=&result.path(&["contact_sensitivities","rows"]).unwrap().as_array().unwrap()[0];
    checked_gradient(n(row,"dobjective_dlog_resistance_k"),contact);
}

#[test]
fn fan_gradient_includes_radiation_air_capacity_and_smooth_convection() {
    let mut correlated=J::parse(CORRELATED).unwrap();
    let patches=J::parse(CONTACT).unwrap().get("radiation").unwrap().clone();
    put(&mut correlated,"radiation",patches);
    for input in [derivative_input(CONTACT),derivative_input(&text(&correlated))] {
        let result=run(&input);
        let expected=central(&input,1e-4,|root,delta| {
            let fan=member(member(root,"hydraulics"),"fan");let speed=n(fan,"speed_ratio");
            put(fan,"speed_ratio",number(speed*delta.exp()));
        });
        let report=result.get("fan_speed_sensitivity").unwrap();
        assert_eq!(report.str_field("status"),Some("available"));
        checked_gradient(n(report,"dobjective_dlog_speed_ratio_k"),expected);
        near(n(report,"dobjective_dlog_speed_ratio_k"),
            n(report,"capacity_contribution_k")+n(report,"convection_contribution_k"),1e-12);
    }
}

#[test]
fn requesting_adjoint_leaves_forward_temperatures_and_energy_unchanged() {
    let input=derivative_input(CONTACT);let with=run(&input);
    let mut plain=input.clone();put(member(&mut plain,"objective"),"gradient",J::Bool(false));
    let without=run(&plain);
    for key in ["solid_temperatures_k","objective","source_w","robin_out_w","coupling_iterations","branches"] {
        assert_eq!(with.get(key),without.get(key),"forward field {key} changed");
    }
    let a=with.get("radiation").unwrap();let b=without.get("radiation").unwrap();
    assert_eq!(a.get("surfaces"),b.get("surfaces"));
    assert_eq!(b.get("adjoint"),Some(&J::Null));
    assert_eq!(b.f64_field("reconstruction_solid_solves"),Some(0.0));
    near(n(a,"solid_solves"),n(b,"solid_solves")+1.0,0.0);
    assert_eq!(a.get("forward_solid_solves"),b.get("forward_solid_solves"));
}

#[test]
fn mean_objective_and_patch_order_retain_the_same_total_derivative() {
    let mut input=derivative_input(CONTACT);
    put(&mut input,"objective",J::parse(r#"{"mean_wall_region":"last-face","gradient":true}"#).unwrap());
    let before=output(&input);assert!(before.status.success(),"{}",String::from_utf8_lossy(&before.stderr));
    let result=J::parse(std::str::from_utf8(&before.stdout).unwrap()).unwrap();
    let expected=central(&input,0.002,|root,delta| {
        let patch=patch_mut(root,"first-face");let value=n(patch,"ambient_temperature_k");
        put(patch,"ambient_temperature_k",number(value+delta));
    });
    checked_gradient(n(patch_gradient(&result,"first-face"),"dobjective_dambient_temperature"),expected);
    let J::Array(rows)=member(member(&mut input,"radiation"),"surfaces") else {panic!()};rows.reverse();
    let after=output(&input);assert!(after.status.success());assert_eq!(before.stdout,after.stdout);
}

#[test]
fn disappearing_radiation_recovers_the_existing_coupled_thermal_gradients() {
    let mut input=derivative_input(CONTACT);
    for name in ["first-face","last-face"] {put(patch_mut(&mut input,name),"emissivity",number(1e-10));}
    let radiative=run(&input);remove(&mut input,"radiation");let ordinary=run(&input);
    for name in ["first-face","last-face"] {
        checked_gradient(n(wall(&radiative,name),"dobjective_dlog_htc"),n(wall(&ordinary,name),"dobjective_dlog_htc"));
    }
    checked_gradient(n(radiative.get("fan_speed_sensitivity").unwrap(),"dobjective_dlog_speed_ratio_k"),
        n(ordinary.get("fan_speed_sensitivity").unwrap(),"dobjective_dlog_speed_ratio_k"));
}

#[test]
fn exhausted_derivatives_never_fall_back_to_frozen_radiation() {
    let mut input=derivative_input(CONTACT);
    put(member(&mut input,"budgets"),"derivative_iterations",number(1.0));
    let failed=output(&input);assert_eq!(failed.status.code(),Some(6));assert!(failed.stdout.is_empty());
    // Primal-only work does not consume the missing derivative budget.
    put(member(&mut input,"objective"),"gradient",J::Bool(false));
    assert!(output(&input).status.success());
}
