use super::*;

const MODEL: &str = include_str!("../specimens/rectangular.scan");
const PROGRAM: &str = include_str!("../specimens/multiple.heaters");
const DATA: &str = include_str!("../specimens/multiple_measurements.csv");
const REQUIREMENTS: &str = include_str!("../specimens/multiple.requirements");

fn with_program(f: impl FnOnce(&CancelGate, &Cx<'_>, &ScanSpecimen, &ScanHeatingProgram)) {
    let gate=CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic);
        let model=ScanSpecimen::from_text(&cx,MODEL,SPECIMEN_BUDGET).unwrap();
        let program=ScanHeatingProgram::from_text(&cx,&model.mesh,model.end_time_s,PROGRAM,HEATER_BUDGET).unwrap();
        f(&gate,&cx,&model,&program);
    });
}
fn selected(result: &Fit) -> Vec<f64> {
    result.schedules.iter().flat_map(|s|s.powers_w.iter().zip(&s.fitted)
        .filter_map(|(&p,&f)|f.then_some(p))).collect()
}

#[test]
fn energy_cap_includes_fixed_heaters_and_returns_a_feasible_measurement_mismatch() {
    with_program(|_,cx,model,program| {
        for substeps in [1,4] {
            let result=fit_with_requirements(cx,model,program,DATA,ScanResidualLoss::Quadratic,
                substeps,Some(REQUIREMENTS)).unwrap();
            assert!(result.converged && result.violation<=1e-8);
            for (actual,expected) in selected(&result).into_iter().zip([8.,16.,16.]) {
                assert!((actual-expected).abs()<2e-3,"{actual} != {expected}");
            }
            assert!((result.objective-4.).abs()<1e-5);
            assert!((result.constraints.absorbed_energy_j[0]-48.).abs()<1e-7);
            assert!((result.measurements[2].3-162.).abs()<1e-6);
            assert!((result.measurements[3].3-0.62).abs()<1e-8);
            assert_eq!(result.schedules[2].powers_w,vec![2.,0.]);
            assert_eq!(result.constraints.bound_values.len(),1);
            assert_eq!(result.constraints.state_values.len(),1);
            assert_eq!(result.constraints.peak_values.len(),1);
            assert_eq!(result.constraints.power_slew_values_w_s.len(),3);
            assert_eq!(result.constraints.path_speed_values.len(),3);
            for schedule in &result.schedules { assert_eq!(schedule.powers_w.last().unwrap().to_bits(),0.0_f64.to_bits()); }
            assert!(result.evaluations<=128);
        }
    });
}

#[test]
fn actual_state_and_fixed_shutdown_rate_limits_change_the_fitted_schedule() {
    with_program(|_,cx,model,program| {
        let limits="FS_SCAN_REQUIREMENTS_V1\nstate 4 0 liquid_mass_fraction * 0.6 1\n";
        let result=fit_with_requirements(cx,model,program,DATA,ScanResidualLoss::Quadratic,2,Some(limits)).unwrap();
        assert!(result.converged && result.violation<=1e-8);
        assert!(result.constraints.state_values[0]<=0.6+1e-8);
        assert!((result.objective-16.).abs()<1e-4);
        for (actual,expected) in selected(&result).into_iter().zip([8.,16.,8.]) {
            assert!((actual-expected).abs()<3e-3);
        }
        let limits="FS_SCAN_REQUIREMENTS_V1\nslew 1 50 10 1\n";
        let result=fit_with_requirements(cx,model,program,DATA,ScanResidualLoss::Quadratic,4,Some(limits)).unwrap();
        assert!(result.converged && result.violation<=1e-8 && result.objective>4.0);
        for (actual,expected) in selected(&result).into_iter().zip([8.,76./3.,10.]) {
            assert!((actual-expected).abs()<3e-3,"{actual} != {expected}");
        }
        assert!(result.constraints.power_slew_values_w_s.last().unwrap()>=&(-10.-1e-8));
        assert_eq!(result.schedules[1].powers_w[3],0.);
    });
}

#[test]
fn constraint_clock_is_bound_before_spatial_exposure_and_sensor_coefficients() {
    with_program(|_,cx,model,program| {
        let limits="FS_SCAN_REQUIREMENTS_V1\nstate 1.5 0 specific_enthalpy_j_kg * 180 1\npeak 1.5 2.5 liquid_mass_fraction 0.9 1 0.001 all\n";
        let text=format!("{DATA}1..2,specific_enthalpy_j_kg,154.75,1,volume,0;1;2;3;4;5\nlag:1:2:0.5:152.5,specific_enthalpy_j_kg,{:.17},1,volume,0;1;2;3;4;5\n",
            157.-2.25*(1.-(-2.0_f64).exp()));
        let result=fit_with_requirements(cx,model,program,&text,ScanResidualLoss::Quadratic,2,Some(limits)).unwrap();
        assert!(result.converged && result.violation<=1e-8 && result.objective<1e-10);
        assert!((result.constraints.state_values[0]-154.75).abs()<1e-5);
        assert_eq!(result.exposure_starts_s[4],Some(1.));
        assert_eq!(result.exposure_starts_s[5],Some(1.));
        assert!(result.sensor_lags[5].is_some());
        assert!(result.constraints.peak_smoothing_gaps[0]>0.);
        assert!(result.constraints.peak_values[0]>=result.constraints.sampled_peak_values[0]);
    });
}

#[test]
fn requirement_options_and_unsupported_inputs_never_fall_back_to_unconstrained_fitting() {
    let parse=|args:&[&str]|options(&args.iter().map(|s|s.to_string()).collect::<Vec<_>>());
    assert!(parse(&["measurements.csv"]).unwrap().requirements.is_none());
    let chosen=parse(&["measurements.csv","--requirements","rules.txt","--heaters","heaters.txt"]).unwrap();
    assert_eq!(chosen.requirements.as_deref(),Some("rules.txt"));
    for args in [vec!["data.csv","--requirements"],vec!["data.csv","--requirements",""],
        vec!["data.csv","--requirements","--huber","1"],
        vec!["data.csv","--requirements","a","--requirements","b"]] {
        assert!(parse(&args).is_err());
    }
    with_program(|gate,cx,model,program| {
        for text in ["FS_SCAN_REQUIREMENTS_V1", "FS_SCAN_REQUIREMENTS_V1\nslew 0 1 1 1\n",
            "FS_SCAN_REQUIREMENTS_V1\npower 2 0 * 1 1\n", "FS_SCAN_REQUIREMENTS_V1\nstate 5 0 temperature_k * 311 1\n"] {
            assert!(fit_with_requirements(cx,model,program,DATA,ScanResidualLoss::Quadratic,1,Some(text)).is_err());
        }
        gate.request();
        assert!(fit_with_requirements(cx,model,program,DATA,ScanResidualLoss::Quadratic,1,Some(REQUIREMENTS)).is_err());
    });
}
