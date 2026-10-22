use super::*;
use fs_ascent::scan_design::constraints::{ScanActuatorBounds, ScanPathSpeedBound,
    ScanPathSpeedBudget, ScanPowerSlewBound, ScanEnergyCap};
use fs_ascent::scan_design::constraints::value::{ValueConstrainedScanStudy, ScanValueConstraintBudget};
use fs_dfo::differential_evolution::DifferentialEvolutionConfig;
use scan_power_fit::population::PopulationScanStudy;

fn path_source(middle: f64, y_middle: f64, asynchronous: bool, powers: [f64; 3]) -> MovingHeatSource {
    let length = Dims([1,0,0,0,0,0]);
    let y = if asynchronous { TimeSignal::Table { times:vec![0.5,1.5], values:vec![0.2,0.6],
        dims:length, interp:Interp::Linear } } else { TimeSignal::Table {
            times:vec![0.,1.,2.], values:vec![0.5,y_middle,0.5], dims:length, interp:Interp::Linear } };
    MovingHeatSource::new([
        TimeSignal::Table { times:vec![0.,1.,2.], values:vec![0.25,middle,0.25], dims:length, interp:Interp::Linear },
        y, TimeSignal::Constant(QtyAny::new(0.5,length))], [1.4,1.2,1.1],
        TimeSignal::Table { times:vec![0.,1.,2.], values:powers.to_vec(), dims:COMPONENT_POWER_DIMS,
            interp:Interp::Linear }, 0., MovingHeatBudget { max_vertices:8, max_elements:6,
                min_active_vertices:1, power_relative_tolerance:1e-11 }).unwrap()
}
fn immutable_factory<'a>(boundary: &'a ThermalBoundaryHistory<'a>, material: &'a ConductivityModel,
    curve: &'a EquilibriumEnthalpyPhaseCurve)
    -> impl Fn(&Cx<'_>, Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'a,'a,'a>, ScanSensitivityError> + Sync + 'a
{
    move |cx, sources| {
        let run = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem { boundary, material,
            element_materials:None, interfaces:None, source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
            EnthalpyHistoryMaterials { table:&[ReferenceEnthalpyMaterial { curve, reference_density_kg_m3:1. }],
                of_element:&[0;6], phase_conductivity:None }, &[150.;8],
            EnthalpyHistoryBudget { spatial:EnthalpyBudget { max_vertices:8, max_elements:6 }, max_steps:2 })?;
        Ok(ScanningEnthalpyRun::new(cx, run, sources, 1, None)?)
    }
}
fn actuators(maximum: f64, with_power: bool) -> ScanActuatorBounds {
    ScanActuatorBounds {
        power_slew: if with_power { vec![ScanPowerSlewBound { source:0,
            maximum_rise_w_s:5., maximum_fall_w_s:5., scale_w_s:2. }] } else {vec![]},
        max_power_segments:2, path_speed:vec![ScanPathSpeedBound { source:0, maximum_m_s:maximum, scale_m_s:0.5 }],
        path_budget:ScanPathSpeedBudget { max_knots_per_source:6, max_segments:4 },
    }
}
fn value_budget() -> ScanValueConstraintBudget {
    ScanValueConstraintBudget { max_rows:32, max_energy_weights:2,
        responses:ScanValueResponseBudget { max_responses:0, max_peak_samples:0 } }
}
fn response_budget() -> ScanResponseBudget {
    ScanResponseBudget { max_responses:0, max_gradient_entries:0, max_transpose_solves:2 }
}

#[test]
fn g1_fixed_off_heater_cannot_hide_diagonal_speed_in_separate_axis_limits() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let study=ScanPowerStudy::new_joint(cx, vec![path_source(0.65,0.5,true,[0.;3])],vec![],
            vec![ScanGeometryControl::HalfWidth {source:0,axis:0,reference_m:1.}],
            objective(vec![ScanObservation { endpoint:0,vertex:0,observable:ScanObservable::SpecificEnthalpy,
                target:150.,scale:1. }],false),immutable_factory(&boundary,&k,&c),4).unwrap();
        let study=ValueConstrainedScanStudy::new_with_actuator_bounds(cx,study,vec![],vec![],vec![],vec![],
            actuators(0.5,false),value_budget()).unwrap();
        let result=study.constraints_at(cx,study.initial_coordinates()).unwrap();
        assert_eq!(result.path_speed_values.len(),4);
        assert!(!result.is_feasible(0.).unwrap());
        let diagonal=&result.path_speed_values[1];
        assert!(diagonal.velocity_m_s.iter().all(|v|v.abs()<0.5));
        close(diagonal.speed_m_s,0.4*2.0_f64.sqrt());
        assert!(diagonal.speed_m_s>0.5);
        let changed_width=study.constraints_at(cx,&[2.0_f64.ln()]).unwrap();
        assert_eq!(changed_width,result);
        assert_eq!(study.evaluations(),0);
    });
}

#[test]
fn g3_speed_and_power_rows_keep_jacobians_batches_and_sqp_handoff() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        for policy in [ScanPowerCoordinates::PositiveLog,ScanPowerCoordinates::NonnegativeLinear] {
            let make=|| ScanPowerStudy::new_joint_with_power_coordinates(cx,
                vec![path_source(0.6,0.7,false,[2.,3.,0.])],controls(),vec![
                    ScanGeometryControl::CenterValue {source:0,axis:0,value:1,reference_m:0.,scale_m:0.25},
                    ScanGeometryControl::CenterValue {source:0,axis:1,value:1,reference_m:0.,scale_m:0.4},
                    ScanGeometryControl::HalfWidth {source:0,axis:0,reference_m:1.}],policy,
                objective(observations(cx,&m),false),immutable_factory(&boundary,&k,&c),64).unwrap();
            let energy=vec![ScanEnergyCap {start_time_s:0.,end_time_s:2.,maximum_j:20.,scale_j:5.}];
            let mut value=ValueConstrainedScanStudy::new_with_actuator_bounds(cx,make(),vec![],energy.clone(),
                vec![],vec![],actuators(0.8,true),value_budget()).unwrap();
            let mut sqp=ConstrainedScanStudy::new_with_actuator_bounds(cx,make(),vec![],energy,
                vec![],vec![],response_budget(),0,actuators(0.8,true),64).unwrap();
            let point=value.initial_coordinates().to_vec();
            let (physical,rows)=value.evaluate_at(cx,&point).unwrap();
            let (_,tangent)=sqp.evaluate_at(cx,&point).unwrap();
            assert_eq!(tangent.residuals,rows.residuals);
            assert_eq!(tangent.path_speed_values,rows.path_speed_values);
            assert_eq!(rows.power_slew_values_w_s,vec![1.,-3.]);
            assert_eq!(rows.path_speed_values.len(),2);
            for j in 0..point.len() {
                let mut plus=point.clone(); let mut minus=point.clone();
                plus[j]+=1e-5; minus[j]-=1e-5;
                let p=value.constraints_at(cx,&plus).unwrap(); let m=value.constraints_at(cx,&minus).unwrap();
                for i in 0..rows.residuals.len() {
                    close(tangent.jacobian[i*point.len()+j],(p.residuals[i]-m.residuals[i])/2e-5);
                }
            }
            {
                let batch=value.value_batch(cx,1).unwrap();
                let (again,constraints)=batch.evaluate_at(cx,&point).unwrap();
                assert_eq!(constraints,rows);
                assert_eq!(again.physical.predictions,physical.physical.predictions);
            }
            let mut polish=value.begin_sqp_at(cx,&point,&[],response_budget(),64).unwrap();
            assert_eq!(polish.original_constraint_count(),rows.residuals.len());
            assert_eq!(polish.sample().ci,rows.residuals);
            let replay=polish.replay(cx).unwrap();
            assert_eq!(replay.constraints,rows);
        }
    });
}

#[test]
fn g1_population_moves_the_heater_but_stops_at_the_physical_speed_boundary() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let mut spec=objective(vec![ScanObservation {endpoint:0,vertex:7,
            observable:ScanObservable::SpecificEnthalpy,target:150.,scale:1.}],false);
        let mut truth=immutable_factory(&boundary,&k,&c)(cx,vec![path_source(0.75,0.5,false,[2.,2.,0.])]).unwrap();
        spec.observations[0].target=truth.advance_event_step(cx,1.,1.,spec.step).unwrap()
            .report.solution().specific_enthalpy_j_kg[7];
        let study=ScanPowerStudy::new_joint(cx,vec![path_source(0.25,0.5,false,[2.,2.,0.])],vec![],
            vec![ScanGeometryControl::CenterValue {source:0,axis:0,value:1,reference_m:0.,scale_m:1.}],
            spec,immutable_factory(&boundary,&k,&c),385).unwrap();
        let study=ValueConstrainedScanStudy::new_with_actuator_bounds(cx,study,vec![],vec![],vec![],vec![],
            actuators(0.25,false),value_budget()).unwrap();
        let mut driver=PopulationScanStudy::new(cx,study,&[[0.25,0.85]],DifferentialEvolutionConfig {
            population_size:8,differential_weight:0.8,crossover_probability:1.,seed:17,
            max_evaluations:384,max_population_entries:8,target:Some(1e-12)},1e-10).unwrap();
        driver.advance(cx,384).unwrap();
        let replay=driver.replay_best(cx).unwrap();
        assert!(replay.score.is_feasible());
        assert!(!replay.target_recovered());
        assert!((replay.coordinates[0]-0.5).abs()<0.002);
        assert!(replay.constraints.path_speed_values.iter().all(|v|v.speed_m_s<=0.250001));
        assert_eq!(replay.sources[0].absorbed_power_values(3).unwrap()[2],0.);
    });
}

#[test]
fn g0_g4_path_caps_refuse_invalid_declarations_and_preserve_cancellation() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|gate,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let make=|| ScanPowerStudy::new_joint(cx,vec![path_source(0.5,0.5,false,[0.;3])],vec![],
            vec![ScanGeometryControl::CenterValue {source:0,axis:0,value:1,reference_m:0.,scale_m:1.}],
            objective(vec![ScanObservation {endpoint:0,vertex:0,observable:ScanObservable::SpecificEnthalpy,
                target:150.,scale:1.}],false),immutable_factory(&boundary,&k,&c),8).unwrap();
        let admit=|a,b| ValueConstrainedScanStudy::new_with_actuator_bounds(cx,make(),vec![],vec![],vec![],vec![],a,b);
        for bad in [-1.,f64::NAN,f64::INFINITY] { assert!(admit(actuators(bad,false),value_budget()).is_err()); }
        let mut duplicate=actuators(0.5,false); duplicate.path_speed.push(duplicate.path_speed[0]);
        assert!(admit(duplicate,value_budget()).is_err());
        let mut budget=actuators(0.5,false); budget.path_budget.max_segments=1;
        assert!(admit(budget,value_budget()).is_err());
        assert!(admit(actuators(0.5,false),ScanValueConstraintBudget {max_rows:1,..value_budget()}).is_err());
        let mut valid=admit(actuators(0.5,false),value_budget()).unwrap();
        gate.request();
        assert!(valid.evaluate_at(cx,&[0.5]).is_err());
        assert_eq!(valid.evaluations(),0);
    });
}
