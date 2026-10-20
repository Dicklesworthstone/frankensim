//! Same-model spatial recovery fixture, not independently calibrated physics.
use super::*;
use fs_ascent::scan_design::ScanGeometryControl;
use fs_conduction::bc::history::march::moving::sensitivity::geometry::MovingHeatGeometryControl as G;

const DEFAULT_TARGET: [f64;5] = [0.2,0.8,1.15,2.,3.];
const INITIAL: [f64;5] = [0.3,0.65,1.4,1.7,2.6];
const ENDPOINTS: [f64;4] = [0.2,0.4,0.6,0.8];

#[derive(Clone,Copy)]
enum Mode { Joint, GeometryOnly, PowerOnly }
struct GeometryFit {
    report: LbfgsReport,
    physical: [f64;5],
    loss: f64,
    evaluations: usize,
    max_gradient_error: f64,
    final_temperatures: Vec<f64>,
}
fn heater(p: [f64;5]) -> Result<MovingHeatSource> {
    if p.iter().any(|v|!v.is_finite()) || p[2] <= 0. || p[3] <= 0. || p[4] <= 0. {
        return Err("geometry target requires finite centers and positive width/powers".into());
    }
    let length=Dims([1,0,0,0,0,0]);
    Ok(MovingHeatSource::new([
        TimeSignal::Ramp {t_start:0.,t_end:0.8,from:QtyAny::new(p[0],length),to:QtyAny::new(p[1],length)},
        TimeSignal::Constant(QtyAny::new(0.45,length)),TimeSignal::Constant(QtyAny::new(0.5,length)),
    ],[p[2],1.6,1.4],TimeSignal::Table {times:vec![0.,0.4,0.8],values:vec![p[3],p[4],0.],
        dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},0.,MovingHeatBudget {
        max_vertices:27,max_elements:48,min_active_vertices:4,power_relative_tolerance:1e-11})?)
}
fn fit_geometry(cx:&Cx<'_>, target:[f64;5], initial:[f64;5], mode:Mode, audit:bool) -> Result<GeometryFit> {
    let (complex,positions)=fs_conduction::fixtures::box_grid([2;3],[1.;3]);
    let mesh=ConductionMesh::new(complex,positions)?;
    let material=ConductivityModel::isotropic_declared(2.)?;
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([109;32]),
        [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
            .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
                liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect())?;
    let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface",ThermalBc::adiabatic())?,&Environment::earth_lab(),vec![],0.)?;
    let factory=|cx:&Cx<'_>,sources:Vec<MovingHeatSource>|
        ->std::result::Result<ScanningEnthalpyRun<'_, '_, '_>,ScanSensitivityError> {
        let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:&boundary,material:&material,
            element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
            EnthalpyHistoryMaterials {table:&[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}],
                of_element:&[0;48],phase_conductivity:None},&[150.;27],
            EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:27,max_elements:48},max_steps:4})?;
        Ok(ScanningEnthalpyRun::new(cx,base,sources,1,None)?)
    };
    let step=EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,
        relative_tolerance:1e-12,linear_restart:27,max_linear_cycles:16,forcing_minimum:1e-12,
        forcing_maximum:1e-3,..Default::default()},max_newton_iterations:64,energy_tolerance_j:1e-8};
    // Generate synthetic targets once. Candidate objectives below run the
    // original forward/adjoint owners; they receive no target parameters.
    let mut truth=factory(cx,vec![heater(target)?])?;
    let mut observations=Vec::new();
    for (endpoint,&time) in ENDPOINTS.iter().enumerate() {
        let solved=truth.advance_event_step(cx,time,time-truth.conduction().time_s(),step)?;
        if solved.end_time_s!=time {return Err("target grid missed a source event".into());}
        for (vertex,&h) in solved.report.solution().specific_enthalpy_j_kg.iter().enumerate() {
            observations.push(ScanObservation {endpoint,vertex,observable:ScanObservable::SpecificEnthalpy,target:h,scale:5.});
        }
    }
    let objective=ScanPowerObjective {endpoints_s:ENDPOINTS.to_vec(),observations,
        sampling:ScanObjectiveSampling::EventLeft,step,
        transpose:LinearConfig {tolerance:1e-10,max_iterations:1000,restart:27},
        budget:ScanObjectiveBudget {max_steps:4,max_observations:108,max_retained_vertices:108,
            max_retained_source_vertices:108,max_schedule_controls:32,max_retained_radiation_entries:0},
        source_reference_time_s:0.};
    let powers=if matches!(mode,Mode::GeometryOnly) {vec![]} else {
        vec![ScanPowerControl {source:0,value:0,reference_w:1.},ScanPowerControl {source:0,value:1,reference_w:1.}]
    };
    let geometry=if matches!(mode,Mode::PowerOnly) {vec![]} else {vec![
        ScanGeometryControl::CenterValue {source:0,axis:0,value:0,reference_m:initial[0],scale_m:0.25},
        ScanGeometryControl::CenterValue {source:0,axis:0,value:1,reference_m:initial[1],scale_m:0.25},
        ScanGeometryControl::HalfWidth {source:0,axis:0,reference_m:initial[2]},
    ]};
    let mut study=ScanPowerStudy::new_joint(cx,vec![heater(initial)?],powers,geometry,objective,factory,255)?;
    let point=study.initial_coordinates().to_vec();
    let mut max_gradient_error=0.0_f64;
    if audit {
        let evaluated=study.evaluate(cx,&point)?;
        for index in 0..point.len() {
            let (mut plus,mut minus)=(point.clone(),point.clone());plus[index]+=1e-5;minus[index]-=1e-5;
            let fd=(study.evaluate(cx,&plus)?.physical.value-study.evaluate(cx,&minus)?.physical.value)/2e-5;
            let analytic=evaluated.gradient[index];
            max_gradient_error=max_gradient_error.max((analytic-fd).abs()/(1.+analytic.abs().max(fd.abs())));
        }
        let accepted_evaluations=study.evaluations();
        assert!(study.sources_at(cx,&[]).is_err());
        assert_eq!(study.evaluations(),accepted_evaluations);
    }
    let mut state=study.begin_optimizer(cx,8)?;
    let report=study.advance_optimizer(cx,&mut state,&StopRule::GradNorm(1e-9),150)?;
    let count=study.evaluations();let exported=study.sources_at(cx,&state.x)?;
    assert_eq!(study.evaluations(),count);
    let source=&exported[0];let p=source.absorbed_power_values(3)?;
    assert_eq!(p[2],0.);
    let physical=[source.geometry_value(G::CenterValue {axis:0,value:0})?,
        source.geometry_value(G::CenterValue {axis:0,value:1})?,source.geometry_value(G::HalfWidth {axis:0})?,p[0],p[1]];
    assert_eq!(source.geometry_value(G::CenterValue {axis:1,value:0})?,0.45);
    assert_eq!(source.geometry_value(G::HalfWidth {axis:1})?,1.6);
    study.extend_evaluation_budget(256)?;
    let replay=study.evaluate(cx,&state.x)?;
    if (replay.physical.value-state.f).abs()>1e-10 {return Err("joint fit did not replay its accepted objective".into());}
    Ok(GeometryFit {report,physical,loss:replay.physical.value,evaluations:study.evaluations(),
        max_gradient_error,final_temperatures:replay.physical.final_step.solution().temperature.clone()})
}

pub(super) fn run_cli(args:&[String]) -> Result<()> {
    if ![0,5].contains(&args.len()) {return Err("--geometry [target_start_m target_end_m target_width_x_m target_first_W target_second_W]".into());}
    let mut target=DEFAULT_TARGET;
    for (value,arg) in target.iter_mut().zip(args) {*value=arg.parse()?;}
    let fitted=with_context(|cx|fit_geometry(cx,target,INITIAL,Mode::Joint,false))?;
    println!("# synthetic same-model geometry recovery; no experimental calibration");
    println!("# status={:?} objective={:.12e} gradient_norm={:.12e} evaluations={}",
        fitted.report.reason,fitted.loss,fitted.report.grad_norm,fitted.evaluations);
    println!("start_x_m,end_x_m,half_width_x_m,first_power_w,second_power_w");
    println!("{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",fitted.physical[0],fitted.physical[1],
        fitted.physical[2],fitted.physical[3],fitted.physical[4]);
    if fitted.report.reason!=StopReason::GradNorm {return Err("joint fit stopped before its gradient threshold".into());}
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g3_joint_path_width_and_power_gradients_match_the_full_objective_and_recover_spatial_targets() {
        let fitted=with_context(|cx|fit_geometry(cx,DEFAULT_TARGET,INITIAL,Mode::Joint,true)).unwrap();
        assert!(fitted.max_gradient_error<1e-6);
        assert_eq!(fitted.report.reason,StopReason::GradNorm);
        assert!(fitted.loss<1e-12);
        for (actual,target) in fitted.physical.iter().zip(DEFAULT_TARGET) {assert!((actual-target).abs()<1e-4);}
        assert!(fitted.final_temperatures.iter().all(|&t|t==310.));
        let powers=with_context(|cx|fit_geometry(cx,DEFAULT_TARGET,INITIAL,Mode::PowerOnly,false)).unwrap();
        assert!(powers.loss>1e-4,"power-only fit must not reproduce the spatial target");
    }
    #[test]
    fn g1_geometry_only_fit_retains_every_authored_power_value() {
        let mut initial=INITIAL;initial[3]=DEFAULT_TARGET[3];initial[4]=DEFAULT_TARGET[4];
        let fitted=with_context(|cx|fit_geometry(cx,DEFAULT_TARGET,initial,Mode::GeometryOnly,false)).unwrap();
        assert!(fitted.loss<1e-12);assert_eq!(&fitted.physical[3..],&DEFAULT_TARGET[3..]);
        for (actual,target) in fitted.physical.iter().zip(DEFAULT_TARGET) {assert!((actual-target).abs()<1e-4);}
    }
    #[test]
    fn g0_bad_geometry_cli_arguments_do_not_start_an_inverse_run() {
        assert!(run_cli(&["0.2".into()]).is_err());
        assert!(heater([0.2,0.8,0.,2.,3.]).is_err());
        assert!(heater([f64::NAN,0.8,1.,2.,3.]).is_err());
    }
}
