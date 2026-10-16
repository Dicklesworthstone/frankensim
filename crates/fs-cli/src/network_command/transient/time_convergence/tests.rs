use super::*;

fn with_cx<T>(run: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(&gate,arena,StreamKey { seed:1,kernel_id:717,tile:0,iteration:0 },
            Budget::INFINITE,ExecMode::Deterministic))
    })
}
fn trace(rows: &[(f64,[f64;2])]) -> Trace {
    let mut result = Trace::new(2,rows.len()-1).unwrap();
    for &(time,field) in rows { result.record(time,&field).unwrap(); }
    result
}

#[test]
fn common_endpoint_comparison_checks_every_node_not_only_final_or_maximum() {
    let coarse = trace(&[(0.0,[300.0,300.0]),(1.0,[302.0,301.0]),(2.0,[303.0,300.0])]);
    let fine = trace(&[(0.0,[300.0,300.0]),(0.5,[301.0,300.0]),
        (1.0,[302.0,300.5]),(1.5,[302.5,300.0]),(2.0,[303.0,300.0])]);
    // Final fields and the all-time maximum agree, but an earlier nonmaximum
    // node does not. A final-only or scalar-peak test would wrongly pass.
    with_cx(|cx| assert_eq!(fine.compare(&coarse,cx).unwrap(),0.5));
}

#[test]
fn shifted_common_times_and_incomplete_fields_are_not_interpolated_or_ignored() {
    let coarse = trace(&[(0.0,[300.0;2]),(1.0,[301.0;2])]);
    let fine = trace(&[(0.0,[300.0;2]),(0.5,[300.5;2]),
        (f64::from_bits(1.0_f64.to_bits()+1),[301.0;2])]);
    with_cx(|cx| assert!(fine.compare(&coarse,cx).is_err()));
    let mut incomplete = Trace::new(2,2).unwrap();
    incomplete.record(0.0,&[300.0;2]).unwrap();
    assert!(incomplete.complete().is_err());
    assert!(incomplete.record(0.0,&[300.0;2]).is_err());
    assert!(incomplete.record(1.0,&[f64::NAN,300.0]).is_err());
    assert!(Trace::values(usize::MAX,2).is_err());
    assert!(Trace::bytes(usize::MAX).is_err());
}

#[test]
fn interval_doubling_keeps_exact_common_times_for_nondyadic_durations() {
    // Halving a global max_step and repeating ceil(duration/max_step) does
    // NOT necessarily double an interval's subdivision count. Keep counts.
    assert_eq!((2.3_f64/1.0).ceil(),3.0);
    assert_eq!((2.3_f64/0.5).ceil(),5.0);
    let mut start = 0.0;
    for (duration,steps) in [(2.3,3),(3.1,4),(0.37,2)] {
        let end = start+duration;
        for i in 1..=steps {
            let old = if i == steps {end} else {start+duration*(i as f64/steps as f64)};
            let new = if i == steps {end} else {start+duration*((2*i) as f64/(2*steps) as f64)};
            assert_eq!(old.to_bits(),new.to_bits());
        }
        start = end;
    }
}

#[test]
fn time_studies_require_a_fixed_horizon_and_two_successive_comparisons() {
    let policy = J::parse(r#"{"max_refinements":4,"consecutive_passes":2,"temperature_tolerance_k":0.01,"max_total_steps":1000,"max_trace_bytes":1048576}"#).unwrap();
    let empty = J::parse("{}").unwrap();
    assert!(Config::parse(&policy,&empty,1000).is_ok());
    for specification in [r#"{"adaptive":{}}"#,r#"{"adjoint":{}}"#,
        r#"{"power_design":{}}"#,r#"{"fan_speed_design":{}}"#,
        r#"{"repeat":{"until_periodic":{}}}"#,r#"{"repeat":{"cycles":2,"fan_controller":{}}}"#] {
        assert!(Config::parse(&policy,&J::parse(specification).unwrap(),1000).is_err());
    }
    let weak = J::parse(r#"{"max_refinements":4,"consecutive_passes":1,"temperature_tolerance_k":0.01,"max_total_steps":1000,"max_trace_bytes":1048576}"#).unwrap();
    assert!(Config::parse(&weak,&empty,1000).is_err());
}
