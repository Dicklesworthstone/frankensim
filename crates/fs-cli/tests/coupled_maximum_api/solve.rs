//! Public AirPath-to-corrected-field-to-fresh-air-state path.
use super::*;
use fs_airflow::conjugate::goal::maximum::solve_linear_maximum;
use fs_conduction::adjoint::{LinearGoalSolveConfig, LinearGoalStop};

fn control() -> LinearGoalSolveConfig {
    LinearGoalSolveConfig { absolute_tolerance: 1e-6, max_primal_iterations: 160,
        check_every: 8, max_defect_corrections: 2 }
}

#[test]
fn corrected_air_bundle_contains_no_stale_initial_references_or_heat_rates() {
    let fixture=Fixture::new();
    with_cx(|_,cx| {
        let (linear,solid,feedback)=configs();
        let initial=vec![300.0;fixture.mesh.vertex_count()];
        let vertices:Vec<_>=(0..initial.len()).collect();
        let paths=[fixture.path(330.0)];
        let stale=paths[0].march(&[300.0]).unwrap();
        assert!(stale.total_heat_rate_w.abs()>1.0);
        let result=solve_linear_maximum(cx,fixture.problem(),None,&paths,linear,
            &initial,&vertices,solid,feedback,control()).unwrap();
        assert_eq!(result.solid.stop,LinearGoalStop::GoalTolerance,"{result:?}");
        assert!(result.solid.primal_iterations>0);
        assert!((result.solid.analysis.nominal_k()-330.0).abs()<1e-6);
        assert!((result.wall_temperatures_k[0]-330.0).abs()<1e-6);
        assert!(result.air[0].total_heat_rate_w.abs()<1e-5);
        assert!((result.air[0].outlet_temperature_k-330.0).abs()<1e-6);
        assert_eq!(result.air[0],paths[0].march(&result.wall_temperatures_k).unwrap());
        assert_ne!(result.air[0],stale);
        assert_eq!(result,solve_linear_maximum(cx,fixture.problem(),None,&paths,linear,
            &initial,&vertices,solid,feedback,control()).unwrap());
    });
}

#[test]
fn unfinished_air_solve_keeps_non_success_and_rebuilds_its_returned_best_field() {
    let fixture=Fixture::new();
    with_cx(|_,cx| {
        let (linear,solid,feedback)=configs();
        let initial=vec![300.0;fixture.mesh.vertex_count()];
        let vertices:Vec<_>=(0..initial.len()).collect();
        let paths=[fixture.path(330.0)];
        for cap in [0,1] {
            let mut config=control();config.max_primal_iterations=cap;config.absolute_tolerance=1e-28;
            let result=solve_linear_maximum(cx,fixture.problem(),None,&paths,linear,
                &initial,&vertices,solid,feedback,config).unwrap();
            assert_ne!(result.solid.stop,LinearGoalStop::GoalTolerance);
            assert!(result.solid.primal_iterations<=cap);
            assert_eq!(result.air[0],paths[0].march(&result.wall_temperatures_k).unwrap());
            let prepared=prepare_linear_maximum(cx,fixture.problem(),None,&paths,linear,
                &initial,solid,feedback).unwrap();
            assert_eq!(result.wall_temperatures_k,prepared.wall_mean_temperatures(cx,&result.solid.temperature).unwrap());
            assert_eq!(result.solid.analysis,prepared.analyze_maximum(cx,&result.solid.temperature,&vertices).unwrap());
            if cap==0 { assert_eq!(result.solid.temperature,initial); }
        }
    });
}

#[test]
fn separate_inlets_and_ordered_segments_are_recomputed_from_the_corrected_solid() {
    // Split the existing cube boundary without changing its geometry or material.
    let mut fixture=Fixture::new();
    fixture.boundary=ThermalBoundaryBuilder::new(&fixture.mesh)
        .region("upstream",|f|f.centroid[0]<0.5,ThermalBc::robin(2.0,300.0).unwrap()).unwrap()
        .region("downstream",|f|f.centroid[0]>=0.5,ThermalBc::robin(2.0,300.0).unwrap()).unwrap().finish().unwrap();
    let areas=[true,false].map(|left|fixture.mesh.boundary().iter()
        .filter(|f|(f.centroid[0]<0.5)==left).map(|f|f.area).sum());
    with_cx(|_,cx| {
        let (linear,solid,feedback)=configs();
        let initial=vec![300.0;fixture.mesh.vertex_count()];
        let vertices:Vec<_>=(0..initial.len()).collect();
        let paths=[path(325.0,24.0,&[("upstream",areas[0],2.0)]),
            path(335.0,30.0,&[("downstream",areas[1],2.0)])];
        let result=solve_linear_maximum(cx,fixture.problem(),None,&paths,linear,
            &initial,&vertices,solid,feedback,control()).unwrap();
        assert_eq!(result.solid.stop,LinearGoalStop::GoalTolerance,"{result:?}");
        assert_eq!(result.air.len(),2);
        assert_eq!(result.air[0].segments[0].inlet_temperature_k,325.0);
        assert_eq!(result.air[1].segments[0].inlet_temperature_k,335.0);
        for i in 0..2 { assert_eq!(result.air[i],paths[i].march(&result.wall_temperatures_k[i..=i]).unwrap()); }
        // With no source, two independent reservoirs exchange heat through the
        // solid. This is a measured numerical cross-check, not an energy certificate.
        assert!((result.air[0].total_heat_rate_w+result.air[1].total_heat_rate_w).abs()<1e-5);
        let ordered=[path(330.0,24.0,&[("upstream",areas[0],2.0),("downstream",areas[1],2.0)])];
        let result=solve_linear_maximum(cx,fixture.problem(),None,&ordered,linear,
            &initial,&vertices,solid,feedback,control()).unwrap();
        assert_eq!(result.solid.stop,LinearGoalStop::GoalTolerance);
        assert_eq!(result.air[0],ordered[0].march(&result.wall_temperatures_k).unwrap());
        assert_eq!(result.air[0].segments[1].inlet_temperature_k,result.air[0].segments[0].outlet_temperature_k);
    });
}
