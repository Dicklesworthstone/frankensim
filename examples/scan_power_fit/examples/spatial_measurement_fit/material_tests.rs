use super::*;

const SPECIMEN: &str = include_str!("../specimens/two_materials.scan");
const MATERIALS: &str = include_str!("../specimens/two_materials.materials");
const HEATERS: &str = include_str!("../specimens/two_materials.heaters");
const DATA: &str = include_str!("../specimens/two_materials_measurements.csv");
fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&gate, &Cx::new(&gate, arena,
            StreamKey { seed:76, kernel_id:51, tile:0, iteration:0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}
fn method(attempts: usize, workers: Option<usize>) -> Method {
    Method::Population(PopulationOptions {
        maximum_w:12., attempts, seed:17, objective_target:0., workers,
    })
}

#[test]
fn distinct_material_masses_and_phase_charts_recover_all_four_heater_powers() {
    with_cx(|_, cx| {
        let model = ScanSpecimen::from_text(cx, SPECIMEN, SPECIMEN_BUDGET).unwrap();
        let program = ScanHeatingProgram::from_text(cx, &model.mesh, model.end_time_s,
            HEATERS, HEATER_BUDGET).unwrap();
        let fit = fit_with_materials(cx, &model, &program, DATA,
            ScanResidualLoss::Huber { threshold: 4. }, 2, None, Method::Sqp, Some(MATERIALS)).unwrap();
        assert!(fit.accepted() && fit.objective < 1e-10, "{}", fit.objective);
        assert_eq!(fit.material_count, 2);
        assert!((fit.reference_mass_kg-3.).abs() < 1e-14);
        for (row, expected) in fit.schedules.iter().zip([[2.,5.,0.], [6.,10.,0.]]) {
            for (&power, target) in row.powers_w.iter().zip(expected) {
                assert!((power-target).abs() < 1e-4, "{power} != {target}");
            }
            assert_eq!(row.powers_w.last().unwrap().to_bits(), 0.0_f64.to_bits());
        }
        let lighter = MATERIALS.replace("material 20 4 12 200", "material 20 4 6 200");
        let fit = fit_with_materials(cx, &model, &program, DATA,
            ScanResidualLoss::Quadratic, 1, None, Method::Sqp, Some(&lighter)).unwrap();
        assert!(fit.accepted() && fit.objective < 1e-10);
        assert!((fit.reference_mass_kg-2.).abs() < 1e-14);
        assert!((fit.schedules[1].powers_w[0]-3.).abs() < 1e-4);
        assert!((fit.schedules[1].powers_w[1]-5.).abs() < 1e-4);
    });
}

#[test]
fn heterogeneous_conductivity_and_nodal_initial_history_reach_the_actual_solver() {
    with_cx(|_, cx| {
        let model = ScanSpecimen::from_text(cx, SPECIMEN, SPECIMEN_BUDGET).unwrap();
        let program = ScanHeatingProgram::from_text(cx, &model.mesh, model.end_time_s,
            HEATERS, HEATER_BUDGET).unwrap();
        let data = format!("{HEADER}\n0.01,specific_enthalpy_j_kg,150,1,vertex,0\n0.01,specific_enthalpy_j_kg,50,1,vertex,4\n");
        for k in [4.,40.] {
            let materials = format!("{}\ninitial 4 50\n", MATERIALS.replace(
                "material 20 4 12 200", &format!("material 20 {k} 12 200")));
            // One search attempt evaluates the exact supplied zero-power start;
            // a separate replay consumes the reserved second physical evaluation.
            let fit = fit_with_materials(cx, &model, &program, &data,
                ScanResidualLoss::Quadratic, 1, None, method(1,None), Some(&materials)).unwrap();
            let expected = (50. + 10.*k*0.01)/(1. + 0.1*k*0.01);
            assert!((fit.measurements[1].3-expected).abs() < 1e-7,
                "k={k}, prediction={}, expected={expected}",fit.measurements[1].3);
            assert!((fit.measurements[0].3-150.).abs() < 1e-10);
            assert_eq!(fit.evaluations,2);
        }
    });
}

#[test]
fn parallel_population_keeps_original_materials_requirements_and_shared_accounting() {
    with_cx(|_, cx| {
        let model = ScanSpecimen::from_text(cx, SPECIMEN, SPECIMEN_BUDGET).unwrap();
        let program = ScanHeatingProgram::from_text(cx, &model.mesh, model.end_time_s,
            HEATERS, HEATER_BUDGET).unwrap();
        let requirements = Some("FS_SCAN_REQUIREMENTS_V1\nenergy 0 2 100 1\nstate 2 4 temperature_k * 331 1\n");
        let serial = fit_with_materials(cx,&model,&program,DATA,ScanResidualLoss::Quadratic,
            2,requirements,method(48,None),Some(MATERIALS)).unwrap();
        for workers in [1,2] {
            let parallel = fit_with_materials(cx,&model,&program,DATA,ScanResidualLoss::Quadratic,
                2,requirements,method(48,Some(workers)),Some(MATERIALS)).unwrap();
            assert_eq!(parallel.schedules,serial.schedules);
            assert_eq!(parallel.objective.to_bits(),serial.objective.to_bits());
            assert_eq!(parallel.constraints.residuals.iter().map(|v|v.to_bits()).collect::<Vec<_>>(),
                serial.constraints.residuals.iter().map(|v|v.to_bits()).collect::<Vec<_>>());
            for (a,b) in parallel.measurements.iter().zip(&serial.measurements) {
                assert_eq!(a.3.to_bits(),b.3.to_bits());
            }
            assert_eq!(parallel.evaluations,49);
            assert_eq!(parallel.material_count,2);
            assert!((parallel.constraints.state_values[0]-330.).abs()<1e-10);
        }
    });
}

#[test]
fn material_flags_invalid_maps_and_cancellation_never_use_homogeneous_fallback() {
    let parse = |args:&[&str]|options(&args.iter().map(|s|s.to_string()).collect::<Vec<_>>());
    assert!(parse(&["data.csv"]).unwrap().materials.is_none());
    assert_eq!(parse(&["data.csv","--materials","body.materials"]).unwrap().materials.as_deref(),Some("body.materials"));
    for args in [vec!["data.csv","--materials"],vec!["data.csv","--materials",""],
        vec!["data.csv","--materials","a","--materials","b"]] { assert!(parse(&args).is_err()); }
    with_cx(|gate,cx| {
        let model=ScanSpecimen::from_text(cx,SPECIMEN,SPECIMEN_BUDGET).unwrap();
        let program=ScanHeatingProgram::from_text(cx,&model.mesh,model.end_time_s,HEATERS,HEATER_BUDGET).unwrap();
        let incomplete=MATERIALS.replace("element 1 20","");
        for method in [Method::Sqp,method(1,None)] {
            assert!(fit_with_materials(cx,&model,&program,DATA,ScanResidualLoss::Quadratic,
                1,None,method,Some(&incomplete)).is_err());
        }
        gate.request();
        assert!(fit_with_materials(cx,&model,&program,DATA,ScanResidualLoss::Quadratic,
            1,None,method(1,None),Some(MATERIALS)).is_err());
    });
}
