use super::*;
use fs_ascent::scan_design::constraints::{ScanEnergyCap, ScanControlBound};
use fs_ascent::scan_design::constraints::value::{ValueConstrainedScanStudy, ScanValueConstraintBudget};
use fs_ascent::{SqpStop, SqpError};

fn physical_limits() -> (Vec<ScanEnergyCap>, Vec<ScanStateBound>, Vec<ScanPeakBound>) {
    (vec![ScanEnergyCap { start_time_s:0., end_time_s:2., maximum_j:5., scale_j:1. }],
     vec![ScanStateBound { response:ScanResponse { endpoint:0, vertex:0,
         observable:ScanObservable::Temperature }, lower:None, upper:Some(311.), scale:1. }],
     vec![ScanPeakBound { response:ScanPeakResponse { first_endpoint:0, last_endpoint:1,
         vertices:(0..8).collect(), observable:ScanObservable::LiquidMassFraction,
         smoothing_scale:0.001 }, maximum:0.9, scale:1. }])
}
fn response_budget() -> ScanResponseBudget {
    ScanResponseBudget { max_responses:2, max_gradient_entries:64, max_transpose_solves:5 }
}
fn value_budget() -> ScanValueConstraintBudget {
    ScanValueConstraintBudget { max_rows:8, max_energy_weights:2,
        responses:ScanValueResponseBudget { max_responses:2, max_peak_samples:16 } }
}

#[test]
fn g3_handoff_keeps_original_rows_and_appends_exact_optimizer_box_derivatives() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(),&Environment::earth_lab(),vec![],0.).unwrap();
        let make=|| ScanPowerStudy::new_joint_with_power_coordinates(cx,sources([1.,1.,0.]),
            controls(),vec![],ScanPowerCoordinates::NonnegativeLinear,
            objective(observations(cx,&m),false),factory(&boundary,&k,&c),20).unwrap();
        let (energy,states,peaks)=physical_limits();
        let direct=ConstrainedScanStudy::new_with_peak_bounds(cx,make(),vec![],energy.clone(),
            states.clone(),peaks.clone(),response_budget(),16,16).unwrap();
        let mut owner=ValueConstrainedScanStudy::new(cx,make(),vec![],energy,states,peaks,value_budget()).unwrap();
        let point=owner.initial_coordinates().to_vec();
        let _=owner.evaluate_at(cx,&point).unwrap();
        let mut session=owner.begin_sqp_at(cx,&point,&[[0.,5.];2],response_budget(),16).unwrap();
        assert_eq!(session.physical_evaluations(),2);
        assert_eq!(session.sample().f,direct.sample().f);
        assert_eq!(session.sample().gradient,direct.sample().gradient);
        assert_eq!(&session.sample().ci[..5],direct.sample().ci.as_slice());
        assert_eq!(&session.sample().ji[..10],direct.sample().ji.as_slice());
        assert_eq!(&session.sample().ji[10..],&[-1.,0.,1.,0.,0.,-1.,0.,1.]);
        let replay=session.replay(cx).unwrap();
        assert_eq!(replay.constraints.residuals,direct.sample().ci);
        assert_eq!(replay.box_residuals,vec![-point[0],point[0]-5.,-point[1],point[1]-5.]);
        assert_eq!(replay.sources[0].absorbed_power_values(3).unwrap()[2],0.);
        drop(session);
        assert_eq!(owner.evaluations(),3);
    });
}

#[test]
fn g1_sqp_polish_reaches_energy_limited_fit_without_replacing_the_physics() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(),&Environment::earth_lab(),vec![],0.).unwrap();
        let readout=ScanReadout::volume_mean(cx,&m,&[0,1,2,3,4,5],6).unwrap();
        let rows:Vec<_>=[152.,157.].into_iter().enumerate().map(|(endpoint,target)|
            ScanSpatialObservation { endpoint,readout:readout.clone(),
                observable:ScanObservable::SpecificEnthalpy,target,scale:10. }).collect();
        let study=ScanPowerStudy::new_joint_with_power_coordinates(cx,sources([1.,1.,0.]),
            vec![ScanPowerControl {source:0,value:0,reference_w:1.},
                 ScanPowerControl {source:0,value:1,reference_w:1.}],vec![],
            ScanPowerCoordinates::NonnegativeLinear,objective(rows,false),factory(&boundary,&k,&c),100).unwrap();
        let (energy,states,peaks)=physical_limits();
        let mut owner=ValueConstrainedScanStudy::new(cx,study,vec![],energy,states,peaks,value_budget()).unwrap();
        let mut session=owner.begin_sqp_at(cx,&[1.,1.],&[[0.,50.];2],response_budget(),16).unwrap();
        let report=session.advance(cx,1e-8,100,80).unwrap();
        assert_eq!(report.stop,SqpStop::Converged);
        let replay=session.replay(cx).unwrap();
        assert!(replay.constraints.is_feasible(1e-8).unwrap());
        assert!((replay.physical.powers_w[0]-2.).abs()<1e-4);
        assert!((replay.physical.powers_w[1]-3.).abs()<1e-4);
        assert!((replay.physical.physical.value-0.02).abs()<1e-8);
        assert!(replay.box_residuals.iter().all(|&g|g<=0.));
    });
}

#[test]
fn g0_sqp_budget_reserves_replay_and_cancellation_retains_the_owner() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|gate,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(),&Environment::earth_lab(),vec![],0.).unwrap();
        let study=ScanPowerStudy::new_joint_with_power_coordinates(cx,sources([1.,1.,0.]),controls(),
            vec![],ScanPowerCoordinates::NonnegativeLinear,objective(observations(cx,&m),false),
            factory(&boundary,&k,&c),2).unwrap();
        let mut owner=ValueConstrainedScanStudy::new(cx,study,
            vec![ScanControlBound {coordinate:0,lower:None,upper:Some(5.),scale:1.}],
            vec![],vec![],vec![],value_budget()).unwrap();
        let point=owner.initial_coordinates().to_vec();
        assert!(owner.begin_sqp_at(cx,&point,&[[0.,5.];2],response_budget(),2).is_err());
        assert_eq!(owner.evaluations(),0);
        {
            let mut session=owner.begin_sqp_at(cx,&point,&[[0.,5.];2],response_budget(),16).unwrap();
            assert_eq!(session.advance(cx,1e-8,100,100).unwrap().stop,SqpStop::EvaluationLimit);
            assert_eq!(session.physical_evaluations(),1);
            session.replay(cx).unwrap();
            assert_eq!(session.physical_evaluations(),2);
        }
        owner.extend_evaluation_budget(5).unwrap();
        owner.evaluate_at(cx,&point).unwrap();
        gate.request();
        assert!(owner.begin_sqp_at(cx,&point,&[],response_budget(),16).is_err());
        assert_eq!(owner.evaluations(),3);
    });
}

#[test]
fn g1_phase_corner_refusal_leaves_value_only_search_available() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(),&Environment::earth_lab(),vec![],0.).unwrap();
        let create=|cx:&Cx<'_>,sources| -> std::result::Result<_,ScanSensitivityError> {
            let run=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {
                boundary:&boundary,material:&k,element_materials:None,interfaces:None,
                source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
                EnthalpyHistoryMaterials {table:&[ReferenceEnthalpyMaterial {curve:&c,reference_density_kg_m3:1.}],
                    of_element:&[0;6],phase_conductivity:None},&[100.;8],
                EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:8,max_elements:6},max_steps:2})?;
            Ok(ScanningEnthalpyRun::new(cx,run,sources,1,None)?)
        };
        let spec=objective(vec![ScanObservation {endpoint:1,vertex:0,
            observable:ScanObservable::LiquidMassFraction,target:0.5,scale:0.1}],false);
        let study=ScanPowerStudy::new_joint_with_power_coordinates(cx,sources([0.;3]),controls(),vec![],
            ScanPowerCoordinates::NonnegativeLinear,spec,create,4).unwrap();
        let mut owner=ValueConstrainedScanStudy::new(cx,study,vec![],vec![],vec![],vec![],value_budget()).unwrap();
        assert!(matches!(owner.begin_sqp_at(cx,&[0.,0.],&[[0.,5.];2],response_budget(),16),
            Err(SqpError::Evaluation(_))));
        assert_eq!(owner.evaluations(),1);
        let (value,rows)=owner.evaluate_at(cx,&[0.,0.]).unwrap();
        assert_eq!(value.physical.value,12.5);
        assert!(rows.is_feasible(0.).unwrap());
        assert_eq!(owner.evaluations(),2);
    });
}
