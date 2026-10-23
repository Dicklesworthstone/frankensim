use super::*;
const MODEL: &str = include_str!("../specimens/rectangular.scan");
const DATA: &str = include_str!("../specimens/rectangular_measurements.csv");
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate=CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate,arena,StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic));
    });
}
fn assert_fit(result: &Fit, expected: [f64; 2]) {
    assert!(result.converged && result.violation<=1e-8 && result.objective<1e-10);
    for (&power,target) in result.powers.iter().zip(expected) {
        assert!((power-target).abs()<2e-4,"{power} != {target}");
    }
    assert_eq!(result.powers[2].to_bits(),0.0_f64.to_bits());
}

#[test]
fn supplied_mass_clock_exposures_and_sensor_lag_share_one_physical_fit() {
    with_cx(|cx| {
        let model=ScanSpecimen::from_text(cx,MODEL,SPECIMEN_BUDGET).unwrap();
        for substeps in [1,4] {
            let result=fit_on_specimen(cx,&model,DATA,[0.,0.],ScanResidualLoss::Quadratic,substeps).unwrap();
            assert_fit(&result,[8.,20.]);
            assert_eq!(result.source_times_s,[0.,1.5,4.]);
            assert_eq!((result.mesh_vertices,result.mesh_elements),(8,6));
            assert!((result.reference_mass_kg-4.).abs()<1e-14);
            assert_eq!(result.measurements.len(),5);
            assert_eq!(result.exposure_starts_s,vec![None,None,Some(1.5),Some(1.5),None]);
        }
        let heavier=ScanSpecimen::from_text(cx,
            &MODEL.replace("reference_density_kg_m3 2","reference_density_kg_m3 4"),SPECIMEN_BUDGET).unwrap();
        let result=fit_on_specimen(cx,&heavier,DATA,[0.,0.],ScanResidualLoss::Quadratic,1).unwrap();
        assert_fit(&result,[16.,40.]);
        assert!((result.reference_mass_kg-8.).abs()<1e-14);
    });
}

#[test]
fn declared_exterior_heat_loss_changes_required_power_not_the_measurements() {
    with_cx(|cx| {
        // All nodes stay on the 310 K plateau. Surface area is 10 m^2, so
        // h=0.1 W/(m^2 K), ambient=300 K removes 10 W throughout both intervals.
        let model=ScanSpecimen::from_text(cx,&MODEL.replace("exterior adiabatic",
            "exterior robin 0.1 300"),SPECIMEN_BUDGET).unwrap();
        let result=fit_on_specimen(cx,&model,DATA,[0.,0.],ScanResidualLoss::Quadratic,4).unwrap();
        assert_fit(&result,[18.,30.]);
        let adiabatic=ScanSpecimen::from_text(cx,MODEL,SPECIMEN_BUDGET).unwrap();
        let wrong=fit_on_specimen(cx,&adiabatic,DATA,[0.,0.],ScanResidualLoss::Quadratic,4).unwrap();
        assert_fit(&wrong,[8.,20.]);
        assert!(result.powers.iter().zip(&wrong.powers).take(2).all(|(a,b)|(a-b-10.).abs()<2e-4));
    });
}

#[test]
fn fitting_uses_supplied_vertex_and_element_counts_not_the_reference_cube() {
    with_cx(|cx| {
        let prefix=MODEL.split("vertex_m").next().unwrap();
        let text=format!("{prefix}vertex_m 0 0 0\nvertex_m 1 0 0\nvertex_m 0 1 0\nvertex_m 0 0 1\nvertex_m 0 0 -1\ntet 0 1 2 3\ntet 0 2 1 4\n");
        let model=ScanSpecimen::from_text(cx,&text,SPECIMEN_BUDGET).unwrap();
        let data=format!("{HEADER}\n1.5,specific_enthalpy_j_kg,151.125,1,volume,0;1\n4,specific_enthalpy_j_kg,153.9375,1,volume,0;1\n");
        let result=fit_on_specimen(cx,&model,&data,[0.,0.],ScanResidualLoss::Quadratic,1).unwrap();
        assert_fit(&result,[0.5,0.75]);
        assert_eq!((result.mesh_vertices,result.mesh_elements),(5,2));
        assert!((result.reference_mass_kg-2./3.).abs()<1e-14);
    });
}

#[test]
fn specimen_options_and_measurement_clock_are_explicit_and_bounded() {
    let args=["data.csv","0","0","--specimen","part.scan","--substeps","4","--huber","1"].map(str::to_string);
    let input=options(&args).unwrap();
    assert_eq!(input.specimen.as_deref(),Some("part.scan"));
    assert_eq!(input.initial,[0.,0.]);
    assert_eq!(input.substeps,4);
    assert!(options(&["data.csv".into()]).unwrap().specimen.is_none());
    for args in [vec!["data.csv","--specimen"],vec!["data.csv","--specimen",""],
        vec!["data.csv","--specimen","--huber","1"],
        vec!["data.csv","--specimen","a","--specimen","b"]] {
        assert!(options(&args.into_iter().map(str::to_string).collect::<Vec<_>>()).is_err());
    }
    with_cx(|cx| {
        let model=ScanSpecimen::from_text(cx,MODEL,SPECIMEN_BUDGET).unwrap();
        let data=parse_for_clock(cx,&model.mesh,DATA,4,[1.5,4.]).unwrap();
        assert_eq!(data.endpoints,vec![0.375,0.75,1.125,1.5,2.125,2.75,3.375,4.]);
        for line in ["4.1,specific_enthalpy_j_kg,165.5,1,vertex,0",
            "4,specific_enthalpy_j_kg,165.5,1,vertex,8"] {
            assert!(parse_for_clock(cx,&model.mesh,&format!("{HEADER}\n{line}\n"),1,[1.5,4.]).is_err());
        }
    });
}
