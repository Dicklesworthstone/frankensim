use super::*;
use fs_dfo::differential_evolution::DifferentialEvolutionStop;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate=CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&gate,&Cx::new(&gate,arena,StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic));
    });
}
fn method(attempts: usize, target: f64) -> Method {
    Method::Population(PopulationOptions {maximum_w:40.,attempts,seed:17,objective_target:target,workers:None})
}
fn melting_data() -> String {
    format!("{HEADER}\n1,liquid_mass_fraction,0.1,0.1,volume,0;1;2;3;4;5\n2,liquid_mass_fraction,0.3,0.1,volume,0;1;2;3;4;5\n")
}
fn model(cx: &Cx<'_>, h: f64) -> ScanSpecimen {
    let mut model=ScanSpecimen::synthetic_reference(cx).unwrap();
    model.initial_specific_enthalpy_j_kg=h;
    // All eight vertices see equal kernel values: exact uniform reference,
    // not an assumption that a general moving footprint heats uniformly.
    model.first_center_m=[0.5;3]; model.second_center_m=[0.5;3];
    model
}

#[test]
fn population_fits_phase_corner_and_cold_start_without_requesting_an_adjoint() {
    with_cx(|_,cx| {
        for (initial_h,expected) in [(100.,[10.,20.]),(80.,[30.,20.])] {
            let model=if initial_h==80. {
                ScanSpecimen::from_text(cx,include_str!("../specimens/cold.scan"),SPECIMEN_BUDGET).unwrap()
            } else {model(cx,initial_h)};
            let program=ScanHeatingProgram::from_specimen(cx,&model,[0.,0.]).unwrap();
            let data=melting_data();
            if initial_h==100. {
                assert!(fit_with_method(cx,&model,&program,&data,ScanResidualLoss::Quadratic,
                    1,None,Method::Sqp).is_err());
            }
            let fit=fit_with_method(cx,&model,&program,&data,ScanResidualLoss::Quadratic,
                1,None,method(2400,1e-8)).unwrap();
            assert!(fit.accepted(),"objective={} violation={}",fit.objective,fit.violation);
            assert!(!fit.converged,"population does not claim SQP stationarity");
            for (&power,expected) in fit.powers.iter().zip(expected) {
                assert!((power-expected).abs()<0.01,"{power} != {expected}");
            }
            assert_eq!(fit.powers[2].to_bits(),0.0_f64.to_bits());
            let status=fit.population.as_ref().unwrap();
            assert_eq!(status.stop,DifferentialEvolutionStop::TargetReached);
            assert_eq!(fit.evaluations,status.attempts+1);
            assert!(fit.evaluations<=2401);
        }
    });
}

#[test]
fn population_keeps_hard_energy_and_reports_unrecovered_targets_after_budget_exhaustion() {
    with_cx(|_,cx| {
        let model=model(cx,100.);
        let program=ScanHeatingProgram::from_specimen(cx,&model,[0.,0.]).unwrap();
        let fit=fit_with_method(cx,&model,&program,&melting_data(),ScanResidualLoss::Quadratic,1,
            Some("FS_SCAN_REQUIREMENTS_V1\nenergy 0 2 20 1\n"),method(1200,1e-8)).unwrap();
        assert!(!fit.accepted());
        assert!(fit.violation<=search::FEASIBILITY_TOLERANCE);
        assert!(fit.constraints.absorbed_energy_j[0]<=20.+1e-8);
        assert!(fit.objective>0.49);
        assert!(fit.measurements[1].3<=0.2+1e-8);
        let status=fit.population.as_ref().unwrap();
        assert_eq!(status.stop,DifferentialEvolutionStop::EvaluationLimit);
        assert_eq!(status.attempts,1200);
        assert_eq!(fit.evaluations,1201);
    });
}

#[test]
fn population_replays_every_requirement_and_retains_fixed_background_and_shutdowns() {
    with_cx(|_,cx| {
        let model=ScanSpecimen::from_text(cx,include_str!("../specimens/rectangular.scan"),SPECIMEN_BUDGET).unwrap();
        let program=ScanHeatingProgram::from_text(cx,&model.mesh,model.end_time_s,
            include_str!("../specimens/multiple.heaters"),HEATER_BUDGET).unwrap();
        let fit=fit_with_method(cx,&model,&program,include_str!("../specimens/multiple_measurements.csv"),
            ScanResidualLoss::Huber {threshold:1.},2,Some(include_str!("../specimens/multiple.requirements")),
            method(80,0.)).unwrap();
        assert!(!fit.accepted());
        assert!(fit.violation<=search::FEASIBILITY_TOLERANCE);
        assert_eq!(fit.schedules[2].powers_w,vec![2.,0.]);
        assert_eq!(fit.constraints.bound_values.len(),1);
        assert_eq!(fit.constraints.state_values.len(),1);
        assert_eq!(fit.constraints.peak_values.len(),1);
        assert_eq!(fit.constraints.power_slew_values_w_s.len(),3);
        assert_eq!(fit.constraints.path_speed_values.len(),3);
        assert!(fit.constraints.absorbed_energy_j[0]<=48.+1e-8);
        assert_eq!(fit.evaluations,81);
        for schedule in &fit.schedules { assert_eq!(schedule.powers_w.last().unwrap().to_bits(),0.0_f64.to_bits()); }
    });
}

#[test]
fn population_options_boxes_physical_failures_and_cancellation_do_not_fall_back() {
    let parse=|args:&[&str]|options(&args.iter().map(|s|s.to_string()).collect::<Vec<_>>());
    assert_eq!(parse(&["data.csv"]).unwrap().method,Method::Sqp);
    assert_eq!(parse(&["data.csv","--population","40","200","--seed","9","--objective-target","0.01"]).unwrap().method,
        Method::Population(PopulationOptions {maximum_w:40.,attempts:200,seed:9,objective_target:0.01,workers:None}));
    for args in [vec!["data.csv","--population","40"],vec!["data.csv","--population","0","10"],
        vec!["data.csv","--population","NaN","10"],vec!["data.csv","--population","40","0"],
        vec!["data.csv","--population","40","100001"],vec!["data.csv","--seed","17"],
        vec!["data.csv","--objective-target","1"],
        vec!["data.csv","--population","40","10","--objective-target","-1"],
        vec!["data.csv","--population","40","10","--population","40","20"]] {
        assert!(parse(&args).is_err());
    }
    with_cx(|gate,cx| {
        let model=model(cx,100.);
        let outside=ScanHeatingProgram::from_specimen(cx,&model,[41.,0.]).unwrap();
        assert!(fit_with_method(cx,&model,&outside,&melting_data(),ScanResidualLoss::Quadratic,
            1,None,method(20,1e-8)).is_err());
        // Reference scales change optimizer coordinates, never the watt box.
        let text="FS_SCAN_HEATERS_V1\nheater 0 hold 2 2 2\nknot 0 0 0.5 0.5 0.5 0 fit 2\nknot 0 1 0.5 0.5 0.5 0 fit 0.5\nknot 0 2 0.5 0.5 0.5 0 fixed\n";
        let program=ScanHeatingProgram::from_text(cx,&model.mesh,model.end_time_s,text,HEATER_BUDGET).unwrap();
        let Method::Population(settings)=method(20,1e-8) else {unreachable!()};
        assert_eq!(settings.bounds(cx,&program).unwrap(),vec![[0.,20.],[0.,80.]]);
        let mut invalid=model;
        invalid.initial_specific_enthalpy_j_kg=1000.;
        assert!(fit_with_method(cx,&invalid,&program,&melting_data(),ScanResidualLoss::Quadratic,
            1,None,method(20,1e-8)).is_err());
        gate.request();
        assert!(settings.bounds(cx,&program).is_err());
    });
}

#[test]
fn parallel_csv_evaluations_preserve_all_replayed_rows_and_do_not_multiply_the_budget() {
    with_cx(|_,cx| {
        let model=ScanSpecimen::from_text(cx,include_str!("../specimens/rectangular.scan"),SPECIMEN_BUDGET).unwrap();
        let program=ScanHeatingProgram::from_text(cx,&model.mesh,model.end_time_s,
            include_str!("../specimens/multiple.heaters"),HEATER_BUDGET).unwrap();
        let data=format!("{}1..2,specific_enthalpy_j_kg,154.75,1,volume,0;1;2;3;4;5\nlag:1:2:0.5:152.5,specific_enthalpy_j_kg,154.9,1,volume,0;1;2;3;4;5\n",
            include_str!("../specimens/multiple_measurements.csv"));
        let requirements=Some(include_str!("../specimens/multiple.requirements"));
        let serial=fit_with_method(cx,&model,&program,&data,ScanResidualLoss::Huber {threshold:1.},2,
            requirements,method(80,0.)).unwrap();
        let bits=|values:&[f64]|values.iter().map(|v|v.to_bits()).collect::<Vec<_>>();
        for workers in [1,2,4] {
            let Method::Population(mut settings)=method(80,0.) else {unreachable!()};
            settings.workers=Some(workers);
            let parallel=fit_with_method(cx,&model,&program,&data,ScanResidualLoss::Huber {threshold:1.},2,
                requirements,Method::Population(settings)).unwrap();
            assert_eq!(parallel.objective.to_bits(),serial.objective.to_bits());
            assert_eq!(bits(&parallel.constraints.residuals),bits(&serial.constraints.residuals));
            assert_eq!(parallel.schedules,serial.schedules);
            for (a,b) in parallel.schedules.iter().zip(&serial.schedules) {
                assert_eq!(bits(&a.powers_w),bits(&b.powers_w));
            }
            for (a,b) in parallel.measurements.iter().zip(&serial.measurements) {
                assert_eq!(a.3.to_bits(),b.3.to_bits());
                assert_eq!(a.4.to_bits(),b.4.to_bits());
            }
            assert_eq!(parallel.evaluations,81);
            assert_eq!(parallel.population.as_ref().unwrap().attempts,80);
            assert_eq!(parallel.population.as_ref().unwrap().settings.workers,Some(workers));
            assert!(!parallel.accepted());
        }
    });
}

#[test]
fn worker_options_require_population_and_unrepresentable_boxes_or_failed_physics_refuse() {
    let parse=|args:&[&str]|options(&args.iter().map(|s|s.to_string()).collect::<Vec<_>>());
    let Method::Population(settings)=parse(&["data.csv","--population","40","80","--workers","4"]).unwrap().method
        else {panic!("population option lost")};
    assert_eq!(settings.workers,Some(4));
    for args in [vec!["data.csv","--workers","4"],vec!["data.csv","--population","40","80","--workers"],
        vec!["data.csv","--population","40","80","--workers","0"],
        vec!["data.csv","--population","40","80","--workers","129"],
        vec!["data.csv","--population","40","80","--workers","2","--workers","4"]] {
        assert!(parse(&args).is_err());
    }
    with_cx(|gate,cx| {
        let mut model=model(cx,100.);
        let text="FS_SCAN_HEATERS_V1\nheater 0 hold 2 2 2\nknot 0 0 0.5 0.5 0.5 0 fit 1e-320\nknot 0 2 0.5 0.5 0.5 0 fixed\n";
        let tiny=ScanHeatingProgram::from_text(cx,&model.mesh,model.end_time_s,text,HEATER_BUDGET).unwrap();
        assert!(settings.bounds(cx,&tiny).is_err());
        let program=ScanHeatingProgram::from_specimen(cx,&model,[0.,0.]).unwrap();
        model.initial_specific_enthalpy_j_kg=1000.;
        assert!(fit_with_method(cx,&model,&program,&melting_data(),ScanResidualLoss::Quadratic,1,None,
            Method::Population(settings)).is_err());
        gate.request();
        assert!(fit_with_method(cx,&model,&program,&melting_data(),ScanResidualLoss::Quadratic,1,None,
            Method::Population(settings)).is_err());
    });
}
