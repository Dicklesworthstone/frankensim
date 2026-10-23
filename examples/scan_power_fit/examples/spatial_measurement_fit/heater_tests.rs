use super::*;

const MODEL: &str = include_str!("../specimens/rectangular.scan");
const PROGRAM: &str = include_str!("../specimens/multiple.heaters");
const DATA: &str = include_str!("../specimens/multiple_measurements.csv");

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate=CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&gate,&Cx::new(&gate,arena,StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic));
    });
}
fn selected(result: &Fit) -> Vec<f64> {
    result.schedules.iter().flat_map(|s|s.powers_w.iter().zip(&s.fitted)
        .filter_map(|(&p,&f)|f.then_some(p))).collect()
}

#[test]
fn multi_heater_inverse_fit_keeps_known_heat_and_fits_three_independent_stages() {
    with_cx(|_,cx| {
        let model=ScanSpecimen::from_text(cx,MODEL,SPECIMEN_BUDGET).unwrap();
        for substeps in [1,4] {
            let program=ScanHeatingProgram::from_text(cx,&model.mesh,model.end_time_s,PROGRAM,HEATER_BUDGET).unwrap();
            let fit=fit_on_program(cx,&model,&program,DATA,ScanResidualLoss::Quadratic,substeps).unwrap();
            assert!(fit.converged && fit.violation<=1e-8 && fit.objective<1e-10);
            for (actual,expected) in selected(&fit).into_iter().zip([8.,16.,24.]) {
                assert!((actual-expected).abs()<1e-4,"{actual} != {expected}");
            }
            assert_eq!(fit.schedules.len(),3);
            assert_eq!(fit.schedules[2].powers_w,vec![2.,0.]);
            assert_eq!(fit.schedules[0].times_s,vec![0.,1.,2.]);
            assert_eq!(fit.schedules[1].times_s,vec![0.,2.,3.,4.]);
            for schedule in &fit.schedules { assert_eq!(schedule.powers_w.last().unwrap().to_bits(),0.0_f64.to_bits()); }
            assert!(fit.evaluations<=128);
        }
        // The very same observations can also be fitted under the WRONG assumed
        // background heat, but the recovered physical powers then change.
        let wrong=PROGRAM.replace("0.5 0.5 2 fixed","0.5 0.5 0 fixed");
        let program=ScanHeatingProgram::from_text(cx,&model.mesh,model.end_time_s,&wrong,HEATER_BUDGET).unwrap();
        let fit=fit_on_program(cx,&model,&program,DATA,ScanResidualLoss::Quadratic,4).unwrap();
        assert!(fit.converged && fit.objective<1e-10);
        for (actual,expected) in selected(&fit).into_iter().zip([10.,18.,28.]) {
            assert!((actual-expected).abs()<1e-4);
        }
    });
}

#[test]
fn all_heater_events_enter_the_grid_before_exposure_and_sensor_weights() {
    with_cx(|gate,cx| {
        let model=ScanSpecimen::from_text(cx,MODEL,SPECIMEN_BUDGET).unwrap();
        let events=[0.25,0.75,1.,1.125,2.,3.,4.];
        let text=format!("{HEADER}\n1..4,specific_enthalpy_j_kg,160,1,volume,0;1;2;3;4;5\nlag:0.25:2:0.5:150,specific_enthalpy_j_kg,154,1,vertex,0\n");
        let data=parse_for_events(cx,&model.mesh,&text,2,4.,&events).unwrap();
        for event in events { assert!(data.endpoints.iter().any(|t|t.to_bits()==event.to_bits())); }
        assert_eq!(data.endpoints.len(),14);
        assert_eq!(data.endpoints.last(),Some(&4.));
        assert_eq!(data.rows[0].endpoint,13);
        let expected:usize=data.rows.iter().zip(&data.starts).map(|(r,s)|r.endpoint-s.unwrap()+1).sum();
        assert_eq!(data.contributions,expected);
        assert_eq!(data.sensor_lags[1],Some(SensorLag {time_constant_s:0.5,initial_reading:150.}));
        for invalid in [vec![0.,1.,4.],vec![1.,1.,4.],vec![2.,1.],vec![5.],vec![f64::NAN]] {
            assert!(parse_for_events(cx,&model.mesh,DATA,1,4.,&invalid).is_err());
        }
        gate.request();
        assert!(parse_for_events(cx,&model.mesh,DATA,1,4.,&[1.,2.,4.]).is_err());
    });
}

#[test]
fn program_options_never_silently_override_starting_powers_or_legacy_behavior() {
    let parse=|args:&[&str]|options(&args.iter().map(|s|s.to_string()).collect::<Vec<_>>());
    let old=parse(&["measurements.csv","0","0"]).unwrap();
    assert!(old.heaters.is_none());
    assert_eq!(old.initial,[0.,0.]);
    let input=parse(&["measurements.csv","--heaters","input.heaters","--specimen","model.scan",
        "--substeps","4","--huber","1.5"]).unwrap();
    assert_eq!(input.heaters.as_deref(),Some("input.heaters"));
    assert_eq!(input.specimen.as_deref(),Some("model.scan"));
    assert_eq!(input.substeps,4);
    assert_eq!(input.loss,ScanResidualLoss::Huber {threshold:1.5});
    for args in [vec!["data.csv","--heaters"],vec!["data.csv","--heaters",""],
        vec!["data.csv","0","0","--heaters","input.heaters"],
        vec!["data.csv","--heaters","a","--heaters","b"],
        vec!["data.csv","--heaters","--substeps","4"]] {
        assert!(parse(&args).is_err());
    }
}
