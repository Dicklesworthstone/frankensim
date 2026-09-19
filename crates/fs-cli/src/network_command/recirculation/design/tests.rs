use super::*;
use super::super::super::tests::{with_cx, close};

fn input(target: usize, limit: f64, gradient: bool, evaluations: usize) -> J {
    let mut root = J::parse(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/recirculated-slab.json"))).unwrap();
    *member_mut(member_mut(&mut root,"objective").unwrap(),"gradient").unwrap() = J::Bool(gradient);
    let design = J::parse(&format!(r#"{{"supply_node":{target},"return_node":3,
        "min_fraction":0,"max_fraction":0.9,"temperature_limit_k":{limit},
        "fraction_tolerance":1e-5,"temperature_tolerance_k":1e-5,
        "max_evaluations":{evaluations}}}"#)).unwrap();
    members(member_mut(&mut root,"recirculation").unwrap()).unwrap().push(("design".into(),design));
    root
}
fn request(target: usize, limit: f64, gradient: bool, evaluations: usize) -> Request {
    Request::parse(&encode(&input(target,limit,gradient,evaluations)).unwrap()).unwrap()
}

#[test]
fn search_keeps_the_actual_passing_field_on_either_side_and_replays() {
    for (target, limit, gradient) in [(0,318.0,true),(1,320.0,true),(0,318.0,false)] {
        let r = request(target,limit,gradient,80);
        let output = execute(&r,&CancelGate::new_clock_free()).unwrap();
        let result = J::parse(&output).unwrap();
        let design = result.get("recirculation_design").unwrap();
        assert_eq!(design.str_field("status"),Some("target-bracketed"));
        let temperature = result.path(&["objective","value_k"]).unwrap().as_f64().unwrap();
        assert!(temperature <= limit && limit-temperature <= 1e-5);
        assert!(design.f64_field("fraction_bracket_width").unwrap() <= 1e-5);
        assert!(design.path(&["failed_endpoint","temperature_k"]).unwrap().as_f64().unwrap() > limit);
        let resolved = design.get("resolved_request").unwrap();
        assert!(resolved.path(&["recirculation","design"]).is_none());
        let replay = Request::parse(&encode(resolved).unwrap()).unwrap();
        let replay = J::parse(&execute(&replay,&CancelGate::new_clock_free()).unwrap()).unwrap();
        close(replay.path(&["objective","value_k"]).unwrap().as_f64().unwrap(),temperature,1e-7);
        let selected = design.f64_field("selected_fraction").unwrap();
        let failed = design.path(&["failed_endpoint","fraction"]).unwrap().as_f64().unwrap();
        assert_eq!(selected > failed, target == 0);
        if gradient {
            assert!(result.path(&["recirculation_sensitivity","links"]).is_some());
        } else { assert_eq!(result.get("recirculation_sensitivity"),Some(&J::Null)); }
    }
}

#[test]
fn both_passing_bounds_are_distinguished_from_an_unbracketed_or_exhausted_search() {
    let r = request(0,400.0,true,2);
    let result = J::parse(&execute(&r,&CancelGate::new_clock_free()).unwrap()).unwrap();
    let study = result.get("recirculation_design").unwrap();
    assert_eq!(study.str_field("status"),Some("both-bounds-feasible"));
    close(study.f64_field("selected_fraction").unwrap(),0.9,0.0);
    assert_eq!(study.get("failed_endpoint"),Some(&J::Null));
    let refused = execute(&request(0,200.0,false,20),&CancelGate::new_clock_free()).unwrap_err();
    assert_eq!(refused.code,"cooling-network-design-bracket");
    let exhausted = execute(&request(0,318.0,true,2),&CancelGate::new_clock_free()).unwrap_err();
    assert_eq!(exhausted.code,"cooling-network-design-budget");
}

#[test]
fn unresolved_designs_cannot_be_silently_consumed_as_fixed_fraction_runs() {
    let r = request(0,318.0,true,80);
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let h = r.surfaces.iter().map(|s|(s.name.clone(),s.h)).collect();
        assert!(r.transport(cx,&flow,&h).is_err());
    });
    let gate=CancelGate::new_clock_free(); gate.request();
    assert!(execute(&r,&gate).is_err());
    for (field,value) in [("max_fraction",J::Number{value:1.0,raw:"1".into()}),
        ("supply_node",J::Number{value:2.0,raw:"2".into()}),
        ("max_evaluations",J::Number{value:1.0,raw:"1".into()})] {
        let mut root=input(0,318.0,true,80);
        let design=member_mut(member_mut(&mut root,"recirculation").unwrap(),"design").unwrap();
        *member_mut(design,field).unwrap()=value;
        assert!(Request::parse(&encode(&root).unwrap()).is_err());
    }
}
