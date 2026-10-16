use super::*;
use super::super::super::tests::{close, with_cx};

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/transient-component-workloads.json"));
fn request() -> Request { Request::parse(FIXTURE).expect("named component schedule") }
fn named(cpu: f64, gpu: f64) -> Workload {
    Workload::Components(BTreeMap::from([("cpu".into(), cpu), ("gpu".into(), gpu)]))
}
fn first_interval(r: &mut Request) {
    let schedule = r.transient.as_mut().unwrap();
    schedule.intervals.truncate(1);
    schedule.total_steps = schedule.intervals[0].steps;
}
fn run(r: &Request) -> J {
    J::parse(&execute(r, &CancelGate::new_clock_free()).expect("actual transient solve")).unwrap()
}

#[test]
fn independent_workloads_move_the_hotspot_and_preserve_window_energy() {
    let r = request();
    let doc = run(&r);
    let result = doc.get("transient").unwrap();
    let history = result.get("history").unwrap().as_array().unwrap();
    assert_eq!(history[15].f64_field("active_vertex"), Some(4.0));
    assert_eq!(history[30].f64_field("active_vertex"), Some(13.0));
    close(result.f64_field("sampled_peak_objective_k").unwrap(), 310.4916199387785, 5e-5);
    close(result.f64_field("sampled_peak_time_s").unwrap(), 60.0, 1e-12);
    close(doc.path(&["objective", "value_k"]).unwrap().as_f64().unwrap(), 305.7727254325046, 5e-5);
    close(result.f64_field("input_energy_j").unwrap(), 1200.0, 1e-6);
    close(result.f64_field("stored_energy_change_j").unwrap(), 1027.9643694287352, 5e-5);
    close(result.f64_field("air_energy_gain_j").unwrap(), 172.03563082139772, 5e-5);
    assert!(result.f64_field("energy_residual_j").unwrap().abs() <= r.limits.heat * 150.0);
    assert_eq!(history[16].get("power_scale"), Some(&J::Null));
    assert_eq!(history[16].path(&["component_powers_w", "gpu"]).unwrap().as_f64(), Some(20.0));
    assert_eq!(history[16].path(&["component_powers_w", "cpu"]).unwrap().as_f64(), Some(0.0));
    assert!(history[1..31].iter().all(|s| (s.f64_field("source_w").unwrap() - 20.0).abs() < 1e-8));
    assert!(history[31..].iter().all(|s| s.f64_field("source_w").unwrap() == 0.0));
}

#[test]
fn identical_total_power_does_not_erase_component_location() {
    let mut cpu = request();
    first_interval(&mut cpu);
    let mut gpu = request();
    first_interval(&mut gpu);
    gpu.transient.as_mut().unwrap().intervals[0].workload = named(0.0, 20.0);
    let a = run(&cpu);
    let b = run(&gpu);
    close(a.path(&["objective", "value_k"]).unwrap().as_f64().unwrap(), 304.11663335714235, 5e-5);
    close(b.path(&["objective", "value_k"]).unwrap().as_f64().unwrap(), 309.86828057433263, 5e-5);
    close(a.path(&["transient", "input_energy_j"]).unwrap().as_f64().unwrap(),
        b.path(&["transient", "input_energy_j"]).unwrap().as_f64().unwrap(), 1e-6);
    // gpu starts at zero nominal watts; rebuilding from its footprint must
    // still turn it on, without dividing by the baseline component power.
    assert_eq!(gpu.solid_data.component_map.as_ref().unwrap().components()[1].watts(), 0.0);
}

#[test]
fn named_source_matches_legacy_scaling_and_reports_unambiguous_controls() {
    let r = request();
    with_cx(|cx| {
        let legacy = Workload::Scale(0.5).prepare(&r, cx).unwrap();
        let independent = named(10.0, 0.0).prepare(&r, cx).unwrap();
        for vertex in 0..r.mesh.vertex_count() {
            close(legacy.source.at(vertex), independent.source.at(vertex), 1e-10);
        }
        close(legacy.expected_power_w.unwrap(), independent.expected_power_w.unwrap(), 1e-12);
        let zero = named(0.0, 0.0).prepare(&r, cx).unwrap();
        assert_eq!(zero.expected_power_w, Some(0.0));
        for vertex in 0..r.mesh.vertex_count() { assert_eq!(zero.source.at(vertex), 0.0); }
        for control in [Workload::Scale(0.5), named(10.0, 0.0)] {
            let doc = J::parse(&format!("{{{}}}", control.render().unwrap())).unwrap();
            assert_ne!(doc.get("power_scale") == Some(&J::Null),
                doc.get("component_powers_w") == Some(&J::Null));
        }
    });
}

#[test]
fn missing_unknown_negative_and_ambiguous_workloads_refuse_before_physics() {
    for text in [r#"{}"#, r#"{"component_powers_w":{}}"#,
        r#"{"component_powers_w":{"cpu":-1,"gpu":0}}"#,
        r#"{"power_scale":1,"component_powers_w":{"cpu":1}}"#,
        r#"{"component_powers_w":[1,2]}"#] {
        assert!(Workload::parse(&J::parse(text).unwrap()).is_err(), "accepted {text}");
    }
    let mut r = request();
    r.limits.coupling = 1;
    r.transient.as_mut().unwrap().intervals[1].workload =
        Workload::Components(BTreeMap::from([("cpu".into(), 0.0), ("unknown".into(), 20.0)]));
    // The later bad workload must be diagnosed before the first step could
    // fail its deliberately insufficient coupling budget.
    assert_eq!(execute(&r, &CancelGate::new_clock_free()).unwrap_err().code, "cooling-network-input");
    with_cx(|cx| {
        assert!(Workload::Components(BTreeMap::from([("cpu".into(), 1.0)])).validate(&r, cx).is_err());
        r.solid_data.component_map = None;
        assert!(named(1.0, 2.0).prepare(&r, cx).is_err());
    });
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(execute(&request(), &gate).is_err());
}
