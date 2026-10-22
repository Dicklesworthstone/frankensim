//! Finite-exposure measurements through the original thermal study and adjoints.
use super::*;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::temporal::ScanTemporalObservation;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::loss::{ScanObservationWithLoss, ScanResidualLoss};
use fs_ascent::scan_design::constraints::value::{ValueConstrainedScanStudy, ScanValueConstraintBudget};
use fs_ascent::scan_design::value::ScanValueResponseBudget;

#[test]
fn g3_window_means_share_full_short_and_value_paths_with_complete_control_gradients() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        for loss in [ScanResidualLoss::Quadratic, ScanResidualLoss::Huber {threshold:0.1}] {
            let mut comparison=None;
            for short in [false,true] {
                let rows:Vec<_>=observations(cx,&m).into_iter().map(|row| {
                    let window=if row.endpoint==0 {ScanTemporalObservation::instantaneous(row)}
                        else {ScanTemporalObservation::trapezoidal_mean(cx,row,&[1.,2.],0,2).unwrap()};
                    ScanObservationWithLoss::new(window,loss).unwrap()
                }).collect();
                let mut study=ScanPowerStudy::new_joint(cx,sources([2.,3.,0.]),controls(),vec![
                    ScanGeometryControl::CenterValue {source:0,axis:0,value:0,reference_m:0.,scale_m:0.5},
                    ScanGeometryControl::HalfWidth {source:0,axis:0,reference_m:1.}],
                    objective(rows,short),factory(&boundary,&k,&c),32).unwrap();
                if short {study.set_checkpointing(cx,checkpoints()).unwrap();}
                let point=study.initial_coordinates().to_vec();
                let adjoint=study.evaluate(cx,&point).unwrap();
                let primal=study.evaluate_value(cx,&point).unwrap();
                assert_eq!(adjoint.physical.value.to_bits(),primal.physical.value.to_bits());
                assert_eq!(adjoint.physical.predictions,primal.physical.predictions);
                for j in 0..point.len() {
                    let (mut plus,mut minus)=(point.clone(),point.clone());
                    plus[j]+=1e-5;minus[j]-=1e-5;
                    let fd=(study.evaluate_value(cx,&plus).unwrap().physical.value
                        -study.evaluate_value(cx,&minus).unwrap().physical.value)/2e-5;
                    close(adjoint.gradient[j],fd);
                }
                let actual=(adjoint.physical.value,adjoint.physical.predictions,
                    adjoint.gradient,adjoint.physical.initial_specific_enthalpy);
                if let Some(expected)=&comparison {assert_eq!(&actual,expected);}
                else {comparison=Some(actual);}
            }
        }
    });
}

#[test]
fn g1_time_mean_fit_recovers_powers_instead_of_treating_exposure_as_a_final_sample() {
    let m=mesh();let k=ConductivityModel::isotropic_declared(2.).unwrap();let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let readout=ScanReadout::volume_mean(cx,&m,&[0,1,2,3,4,5],6).unwrap();
        for averaged in [true,false] {
            let rows:Vec<_>=[152.,154.5].into_iter().enumerate().map(|(endpoint,target)| {
                let row=ScanSpatialObservation {endpoint,readout:readout.clone(),
                    observable:ScanObservable::SpecificEnthalpy,target,scale:10.};
                if averaged && endpoint==1 {ScanTemporalObservation::trapezoidal_mean(cx,row,&[1.,2.],0,2).unwrap()}
                else {ScanTemporalObservation::instantaneous(row)}
            }).collect();
            let mut study=ScanPowerStudy::new(cx,sources([1.,1.,0.]),controls(),objective(rows,true),
                factory(&boundary,&k,&c),160).unwrap();
            study.set_checkpointing(cx,checkpoints()).unwrap();
            let mut optimizer=study.begin_optimizer(cx,6).unwrap();
            study.advance_optimizer(cx,&mut optimizer,&StopRule::GradNorm(1e-9),100).unwrap();
            let fitted=study.evaluate_value(cx,&optimizer.x).unwrap();
            assert!(fitted.physical.value<1e-12);
            close(fitted.powers_w[0],2.);
            close(fitted.powers_w[1],if averaged {5.} else {2.5});
            assert_eq!(study.sources_at(cx,&optimizer.x).unwrap()[0].absorbed_power_values(3).unwrap()[2],0.);
        }
    });
}

#[test]
fn g3_temporal_objectives_do_not_average_away_hard_point_and_peak_requirements() {
    let m=mesh();let k=ConductivityModel::isotropic_declared(2.).unwrap();let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let row=ScanSpatialObservation {endpoint:1,
            readout:ScanReadout::volume_mean(cx,&m,&[0,1,2,3,4,5],6).unwrap(),
            observable:ScanObservable::SpecificEnthalpy,target:154.5,scale:10.};
        let mut rows=Vec::new();
        for averaged in [false,true] {
            let observation=if averaged {ScanTemporalObservation::trapezoidal_mean(cx,row.clone(),&[1.,2.],0,2).unwrap()}
                else {ScanTemporalObservation::instantaneous(row.clone())};
            let study=ScanPowerStudy::new(cx,sources([2.,5.,0.]),controls(),objective(vec![observation],false),
                factory(&boundary,&k,&c),4).unwrap();
            let mut study=ValueConstrainedScanStudy::new(cx,study,vec![],vec![],
                vec![ScanStateBound {response:ScanResponse {endpoint:1,vertex:0,observable:ScanObservable::SpecificEnthalpy},
                    lower:None,upper:Some(155.),scale:1.}],
                vec![ScanPeakBound {response:ScanPeakResponse {first_endpoint:0,last_endpoint:1,vertices:(0..8).collect(),
                    observable:ScanObservable::SpecificEnthalpy,smoothing_scale:0.01},maximum:155.,scale:1.}],
                ScanValueConstraintBudget {max_rows:2,max_energy_weights:0,
                    responses:ScanValueResponseBudget {max_responses:2,max_peak_samples:16}}).unwrap();
            let point=study.initial_coordinates().to_vec();
            let (physical,constraints)=study.evaluate_at(cx,&point).unwrap();
            if averaged {assert!(physical.physical.value<1e-20);} else {assert!(physical.physical.value>0.01);}
            assert!(!constraints.is_feasible(0.).unwrap());
            rows.push(constraints);
        }
        assert_eq!(rows[0],rows[1]);
    });
}

#[test]
fn g0_time_sample_budget_and_grid_mismatch_refuse_without_bypassing_cancellation() {
    let m=mesh();let k=ConductivityModel::isotropic_declared(2.).unwrap();let c=curve();
    with_cx(|gate,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let row=ScanObservation {endpoint:1,vertex:0,observable:ScanObservable::SpecificEnthalpy,target:154.5,scale:10.};
        for wrong_grid in [false,true] {
            let observation=ScanTemporalObservation::trapezoidal_mean(cx,row,&[1.,2.],0,2).unwrap();
            let mut spec=objective(vec![observation],false);
            if wrong_grid {spec.endpoints_s[0]=0.5;} else {spec.budget.max_observations=1;}
            let mut study=ScanPowerStudy::new(cx,sources([2.,5.,0.]),controls(),spec,factory(&boundary,&k,&c),2).unwrap();
            let point=study.initial_coordinates().to_vec();
            assert!(study.evaluate_value(cx,&point).is_err());
            assert_eq!(study.evaluations(),1);
        }
        let observation=ScanTemporalObservation::trapezoidal_mean(cx,row,&[1.,2.],0,2).unwrap();
        let mut study=ScanPowerStudy::new(cx,sources([2.,5.,0.]),controls(),objective(vec![observation],false),
            factory(&boundary,&k,&c),2).unwrap();
        let point=study.initial_coordinates().to_vec();gate.request();
        assert!(study.evaluate_value(cx,&point).is_err());
        assert_eq!(study.evaluations(),0);
    });
}

#[test]
fn g3_sensor_lag_retains_full_checkpoint_and_value_control_gradients() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        for loss in [ScanResidualLoss::Quadratic,ScanResidualLoss::Huber {threshold:0.1}] {
            let mut comparison=None;
            for short in [false,true] {
                let rows:Vec<_>=observations(cx,&m).into_iter().map(|row| {
                    let response=if row.endpoint==0 {ScanTemporalObservation::instantaneous(row)}
                        else {
                            let initial=match row.observable {
                                ScanObservable::Temperature=>308.,
                                ScanObservable::SpecificEnthalpy=>150.,
                                ScanObservable::LiquidMassFraction=>0.5,
                            };
                            ScanTemporalObservation::first_order_lag(cx,row,&[1.,2.],0,0.5,initial,2).unwrap()
                        };
                    ScanObservationWithLoss::new(response,loss).unwrap()
                }).collect();
                let mut study=ScanPowerStudy::new_joint(cx,sources([2.,3.,0.]),controls(),vec![
                    ScanGeometryControl::CenterValue {source:0,axis:0,value:0,reference_m:0.,scale_m:0.5},
                    ScanGeometryControl::HalfWidth {source:0,axis:0,reference_m:1.}],
                    objective(rows,short),factory(&boundary,&k,&c),32).unwrap();
                if short {study.set_checkpointing(cx,checkpoints()).unwrap();}
                let point=study.initial_coordinates().to_vec();
                let adjoint=study.evaluate(cx,&point).unwrap();
                let primal=study.evaluate_value(cx,&point).unwrap();
                assert_eq!(adjoint.physical.value.to_bits(),primal.physical.value.to_bits());
                assert_eq!(adjoint.physical.predictions,primal.physical.predictions);
                for j in 0..point.len() {
                    let (mut plus,mut minus)=(point.clone(),point.clone());
                    plus[j]+=1e-5; minus[j]-=1e-5;
                    let fd=(study.evaluate_value(cx,&plus).unwrap().physical.value
                        -study.evaluate_value(cx,&minus).unwrap().physical.value)/2e-5;
                    close(adjoint.gradient[j],fd);
                }
                let actual=(adjoint.physical.value,adjoint.physical.predictions,
                    adjoint.gradient,adjoint.physical.initial_specific_enthalpy);
                if let Some(expected)=&comparison {assert_eq!(&actual,expected);}
                else {comparison=Some(actual);}
            }
        }
    });
}

#[test]
fn g1_sensor_lag_recovers_heating_instead_of_fitting_delayed_readings_as_material_states() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let readout=ScanReadout::volume_mean(cx,&m,&[0,1,2,3,4,5],6).unwrap();
        let delayed=152.+5.*(1.-0.5*(1.-(-2.0_f64).exp()));
        for modeled in [true,false] {
            let rows:Vec<_>=[152.,delayed].into_iter().enumerate().map(|(endpoint,target)| {
                let row=ScanSpatialObservation {endpoint,readout:readout.clone(),
                    observable:ScanObservable::SpecificEnthalpy,target,scale:10.};
                if modeled && endpoint==1 {
                    ScanTemporalObservation::first_order_lag(cx,row,&[1.,2.],0,0.5,152.,2).unwrap()
                } else {ScanTemporalObservation::instantaneous(row)}
            }).collect();
            let mut study=ScanPowerStudy::new(cx,sources([1.,1.,0.]),controls(),objective(rows,true),
                factory(&boundary,&k,&c),160).unwrap();
            study.set_checkpointing(cx,checkpoints()).unwrap();
            let mut optimizer=study.begin_optimizer(cx,6).unwrap();
            study.advance_optimizer(cx,&mut optimizer,&StopRule::GradNorm(1e-9),100).unwrap();
            let fitted=study.evaluate_value(cx,&optimizer.x).unwrap();
            assert!(fitted.physical.value<1e-12);
            close(fitted.powers_w[0],2.);
            close(fitted.powers_w[1],if modeled {5.} else {delayed-152.});
            assert_eq!(study.sources_at(cx,&optimizer.x).unwrap()[0].absorbed_power_values(3).unwrap()[2],0.);
        }
    });
}

#[test]
fn g3_sensor_readings_do_not_replace_true_thermal_limits() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let row=ScanSpatialObservation {endpoint:1,
            readout:ScanReadout::volume_mean(cx,&m,&[0,1,2,3,4,5],6).unwrap(),
            observable:ScanObservable::SpecificEnthalpy,target:154.,scale:10.};
        let mut limits=Vec::new();
        for modeled in [false,true] {
            let observation=if modeled {
                ScanTemporalObservation::first_order_lag(cx,row.clone(),&[1.,2.],0,0.5,152.,2).unwrap()
            } else {ScanTemporalObservation::instantaneous(row.clone())};
            let study=ScanPowerStudy::new(cx,sources([2.,5.,0.]),controls(),objective(vec![observation],false),
                factory(&boundary,&k,&c),2).unwrap();
            let mut study=ValueConstrainedScanStudy::new(cx,study,vec![],vec![],vec![],
                vec![ScanPeakBound {response:ScanPeakResponse {first_endpoint:0,last_endpoint:1,vertices:(0..8).collect(),
                    observable:ScanObservable::SpecificEnthalpy,smoothing_scale:0.01},maximum:155.,scale:1.}],
                ScanValueConstraintBudget {max_rows:1,max_energy_weights:0,
                    responses:ScanValueResponseBudget {max_responses:1,max_peak_samples:16}}).unwrap();
            let point=study.initial_coordinates().to_vec();
            let (_,constraint)=study.evaluate_at(cx,&point).unwrap();
            assert!(!constraint.is_feasible(0.).unwrap());
            limits.push(constraint);
        }
        assert_eq!(limits[0],limits[1]);
    });
}
