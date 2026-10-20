//! Reuse the original spatial specimen to test robust fitting, not a new solver.
use super::*;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::loss::{
    ScanObservationWithLoss, ScanResidualLoss};

#[test]
fn g3_robust_spatial_gradient_preserves_raw_data_and_matches_full_checkpointed_paths() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let mut compared=Vec::new();
        for short in [false,true] {
            let rows=observations(cx,&m);
            let wrapped=rows.iter().cloned().enumerate().map(|(i,row)|ScanObservationWithLoss::new(row,
                if i==2 {ScanResidualLoss::Quadratic} else {ScanResidualLoss::Huber {threshold:0.1}}).unwrap()).collect();
            let geometry=vec![ScanGeometryControl::CenterValue {source:0,axis:0,value:0,reference_m:0.,scale_m:0.5},
                ScanGeometryControl::HalfWidth {source:0,axis:0,reference_m:1.}];
            let mut study=ScanPowerStudy::new_joint(cx,sources([2.,3.,0.]),controls(),geometry,
                objective(wrapped,short),factory(&boundary,&k,&c),20).unwrap();
            if short {study.set_checkpointing(cx,checkpoints()).unwrap();}
            let point=study.initial_coordinates().to_vec(); let result=study.evaluate(cx,&point).unwrap();
            let manual:f64=result.physical.standardized_residuals.iter().enumerate().map(|(i,&r)|
                (if i==2 {ScanResidualLoss::Quadratic} else {ScanResidualLoss::Huber {threshold:0.1}})
                    .value_derivative(r).unwrap().0).sum();
            assert_eq!(result.physical.value,manual);
            for (row,(&prediction,&r)) in rows.iter().zip(result.physical.predictions.iter()
                .zip(&result.physical.standardized_residuals)) {assert_eq!(r,(prediction-row.target)/row.scale);}
            assert!(result.physical.standardized_residuals.iter().any(|r|r.abs()>0.1));
            for j in 0..point.len() {
                let (mut plus,mut minus)=(point.clone(),point.clone());plus[j]+=1e-5;minus[j]-=1e-5;
                let fd=(study.evaluate(cx,&plus).unwrap().physical.value-study.evaluate(cx,&minus).unwrap().physical.value)/2e-5;
                close(result.gradient[j],fd);
            }
            compared.push((result.physical.value,result.physical.predictions,result.gradient,
                result.physical.initial_specific_enthalpy));
        }
        assert_eq!(compared[0],compared[1]);
    });
}

#[test]
fn g3_quadratic_wrapper_is_exact_and_huber_does_not_relax_point_or_peak_limits() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let rows=observations(cx,&m);
        let mut plain=ScanPowerStudy::new(cx,sources([2.,3.,0.]),controls(),objective(rows.clone(),false),
            factory(&boundary,&k,&c),4).unwrap();
        let point=plain.initial_coordinates().to_vec();let original=plain.evaluate(cx,&point).unwrap();
        let mut limits=Vec::new();
        for loss in [ScanResidualLoss::Quadratic,ScanResidualLoss::Huber {threshold:0.1}] {
            let wrapped=rows.iter().cloned().map(|r|ScanObservationWithLoss::new(r,loss).unwrap()).collect();
            let study=ScanPowerStudy::new(cx,sources([2.,3.,0.]),controls(),objective(wrapped,false),
                factory(&boundary,&k,&c),4).unwrap();
            let mut study=ConstrainedScanStudy::new_with_peak_bounds(cx,study,vec![],vec![],vec![ScanStateBound {
                response:ScanResponse {endpoint:0,vertex:0,observable:ScanObservable::SpecificEnthalpy},
                lower:None,upper:Some(151.),scale:1.}],vec![ScanPeakBound {
                response:ScanPeakResponse {first_endpoint:0,last_endpoint:1,vertices:(0..8).collect(),
                    observable:ScanObservable::LiquidMassFraction,smoothing_scale:0.01},
                maximum:0.51,scale:0.1}],ScanResponseBudget {max_responses:2,max_gradient_entries:64,
                    max_transpose_solves:5},16,8).unwrap();
            let (physical,requirements)=study.replay_with_constraints(cx).unwrap();
            assert_eq!(original.physical.predictions,physical.physical.predictions);
            assert_eq!(original.physical.standardized_residuals,physical.physical.standardized_residuals);
            if loss==ScanResidualLoss::Quadratic {
                assert_eq!(original.physical.value,physical.physical.value);
                assert_eq!(original.gradient,physical.gradient);
            } else {assert!(physical.physical.value<original.physical.value);}
            assert!(requirements.max_scaled_violation()>0.);
            limits.push((requirements.state_values,requirements.residuals,requirements.jacobian));
        }
        assert_eq!(limits[0],limits[1]);
    });
}

#[test]
fn g0_robust_wrapper_retains_foreign_mesh_and_cancellation_refusals() {
    let m=mesh();let other=mesh();let k=ConductivityModel::isotropic_declared(2.).unwrap();let c=curve();
    with_cx(|gate,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let wrapped=observations(cx,&other).into_iter().map(|r|ScanObservationWithLoss::new(r,
            ScanResidualLoss::Huber {threshold:1.}).unwrap()).collect();
        let mut study=ScanPowerStudy::new(cx,sources([2.,3.,0.]),controls(),objective(wrapped,false),
            factory(&boundary,&k,&c),4).unwrap();
        let point=study.initial_coordinates().to_vec();assert!(study.evaluate(cx,&point).is_err());
        let spent=study.evaluations();
        gate.request();
        assert!(study.evaluate(cx,&point).is_err());
        assert_eq!(study.evaluations(),spent);
    });
}
