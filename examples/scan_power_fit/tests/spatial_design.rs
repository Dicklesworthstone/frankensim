//! Spatial observations through the original L-BFGS/SQP scan study.
use fs_ascent::scan_design::{ScanGeometryControl, ScanPowerControl, ScanPowerCoordinates,
    ScanPowerObjective, ScanPowerStudy};
use fs_ascent::scan_design::constraints::{ConstrainedScanStudy, ScanPeakBound, ScanStateBound,
    ScanResponse, ScanResponseBudget};
use fs_ascent::{StopReason, StopRule};
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun, EnthalpyHistoryMaterials,
    EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::{ScanSensitivityError,
    objective::{ScanObservable, ScanObservation, ScanObjectiveSampling, ScanObjectiveBudget}};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::spatial::{
    ScanReadout, ScanSpatialObservation, evaluate_scan_observations};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::responses::peaks::ScanPeakResponse;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::checkpoint::ScanCheckpointBudget;
use fs_conduction::bc::history::march::moving::{MovingHeatSource, MovingHeatBudget};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionMesh, ConductivityModel, ScalarField, ThermalBc,
    ThermalBoundaryBuilder, LinearConfig, COMPONENT_POWER_DIMS};
use fs_exec::{Cx, Budget, CancelGate, ExecMode, StreamKey};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_qty::{Dims, QtyAny};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{TimeSignal, Interp};
use fs_solver::NewtonKrylovConfig;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&gate, &Cx::new(&gate, arena,
        StreamKey { seed: 76, kernel_id: 51, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)));
}
fn mesh() -> ConductionMesh {
    let (complex, positions) = fs_conduction::fixtures::box_grid([1,1,1], [1.;3]);
    ConductionMesh::new(complex, positions).unwrap()
}
fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([76;32]),
        [(-200.,280.,0.), (100.,310.,0.), (200.,310.,1.), (900.,380.,1.)].into_iter()
            .map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h, temperature_k:t,
                liquid_mass_fraction:f, bulk_density_kg_m3:1. }).collect()).unwrap()
}
fn sources(p: [f64;3]) -> Vec<MovingHeatSource> {
    let length = Dims([1,0,0,0,0,0]);
    vec![MovingHeatSource::new([
        TimeSignal::Table { times:vec![0.,1.,2.], values:vec![0.25,0.75,0.75], dims:length, interp:Interp::Hold },
        TimeSignal::Constant(QtyAny::new(0.5,length)),
        TimeSignal::Constant(QtyAny::new(0.5,length))], [1.4,1.2,1.1],
        TimeSignal::Table { times:vec![0.,1.,2.], values:p.to_vec(), dims:COMPONENT_POWER_DIMS,
            interp:Interp::Hold }, 0., MovingHeatBudget { max_vertices:8, max_elements:6,
            min_active_vertices:1, power_relative_tolerance:1e-11 }).unwrap()]
}
type Factory<'a> = Box<dyn FnMut(&Cx<'_>, Vec<MovingHeatSource>)
    -> Result<ScanningEnthalpyRun<'a,'a,'a>, ScanSensitivityError> + 'a>;
fn factory<'a>(boundary: &'a ThermalBoundaryHistory<'a>, material: &'a ConductivityModel,
    curve: &'a EquilibriumEnthalpyPhaseCurve) -> Factory<'a>
{
    Box::new(move |cx, sources| {
        let run = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {
            boundary, material, element_materials:None, interfaces:None,
            source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
            EnthalpyHistoryMaterials { table:&[ReferenceEnthalpyMaterial {
                curve, reference_density_kg_m3:1. }], of_element:&[0;6], phase_conductivity:None },
            &[150.;8], EnthalpyHistoryBudget { spatial:EnthalpyBudget { max_vertices:8, max_elements:6 },
                max_steps:2 })?;
        Ok(ScanningEnthalpyRun::new(cx, run, sources, 1, None)?)
    })
}
fn controls() -> Vec<ScanPowerControl> {
    vec![ScanPowerControl { source:0, value:0, reference_w:2. },
        ScanPowerControl { source:0, value:1, reference_w:3. }]
}
fn objective<O>(observations: Vec<O>, short: bool) -> ScanPowerObjective<O> {
    let retained = if short {1} else {2};
    ScanPowerObjective { endpoints_s:vec![1.,2.], observations,
        sampling:ScanObjectiveSampling::EventLeft,
        step:EnthalpyStepConfig { newton:NewtonKrylovConfig { absolute_tolerance:1e-11,
            relative_tolerance:1e-12, linear_restart:8, max_linear_cycles:16,
            forcing_minimum:1e-12, forcing_maximum:1e-3, ..Default::default() },
            max_newton_iterations:64, energy_tolerance_j:1e-8 },
        transpose:LinearConfig { tolerance:1e-10, max_iterations:1000, restart:8 },
        budget:ScanObjectiveBudget { max_steps:2, max_observations:16,
            max_retained_vertices:8*retained, max_retained_source_vertices:8*retained,
            max_schedule_controls:16, max_retained_radiation_entries:0 },
        source_reference_time_s:0. }
}
fn checkpoints() -> ScanCheckpointBudget {
    ScanCheckpointBudget { segment_steps:1, max_checkpoint_vertices:24, max_forward_solves:4 }
}
fn observations<'a>(cx: &Cx<'_>, m: &'a ConductionMesh) -> Vec<ScanSpatialObservation<'a>> {
    vec![ScanSpatialObservation { endpoint:0,
        readout:ScanReadout::volume_mean(cx,m,&[0,1],2).unwrap(),
        observable:ScanObservable::SpecificEnthalpy, target:153., scale:10. },
        ScanSpatialObservation { endpoint:1,
            readout:ScanReadout::weighted_mean(cx,m,&[(0,1.),(3,2.),(7,1.)],3).unwrap(),
            observable:ScanObservable::LiquidMassFraction, target:0.57, scale:0.1 },
        ScanSpatialObservation { endpoint:1,
            readout:ScanReadout::in_tetrahedron(cx,m,0,[1.,2.,3.,4.]).unwrap(),
            observable:ScanObservable::Temperature, target:309., scale:2. }]
}
fn close(a: f64, b: f64) { assert!((a-b).abs() <= 2e-6*(1.+a.abs().max(b.abs())), "{a} != {b}"); }

#[test]
fn g3_spatial_study_matches_original_objective_and_complete_coordinate_gradients() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        for short in [false,true] {
            let rows=observations(cx,&m); let spec=objective(rows.clone(),short);
            let mut study=ScanPowerStudy::new_joint(cx,sources([2.,3.,0.]),controls(),vec![
                ScanGeometryControl::CenterValue { source:0,axis:0,value:0,reference_m:0.,scale_m:0.5 },
                ScanGeometryControl::HalfWidth { source:0,axis:0,reference_m:1. }],
                spec,factory(&boundary,&k,&c),20).unwrap();
            if short { study.set_checkpointing(cx,checkpoints()).unwrap(); }
            let point=study.initial_coordinates().to_vec(); let actual=study.evaluate(cx,&point).unwrap();
            let direct=objective(rows.clone(),short);
            let controlled=study.sources_at(cx,&point).unwrap();
            let physical=evaluate_scan_observations(cx,factory(&boundary,&k,&c)(cx,controlled).unwrap(),
                &direct.endpoints_s,&rows,direct.sampling,direct.step,direct.transpose,direct.budget,
                &[],&[],ScanResponseBudget {max_responses:0,max_gradient_entries:0,max_transpose_solves:2},
                0,short.then_some(checkpoints())).unwrap().0;
            assert_eq!(actual.physical.value,physical.value);
            assert_eq!(actual.physical.predictions,physical.predictions);
            assert_eq!(actual.physical.sources,physical.sources);
            assert_eq!(actual.physical.initial_specific_enthalpy,physical.initial_specific_enthalpy);
            assert_eq!(actual.checkpoint_work.is_some(),short);
            for j in 0..point.len() {
                let (mut plus,mut minus)=(point.clone(),point.clone());
                plus[j]+=1e-5; minus[j]-=1e-5;
                let fd=(study.evaluate(cx,&plus).unwrap().physical.value
                    -study.evaluate(cx,&minus).unwrap().physical.value)/2e-5;
                close(actual.gradient[j],fd);
            }
            let spent=study.evaluations(); let exported=study.sources_at(cx,&point).unwrap();
            assert_eq!(exported[0].absorbed_power_values(3).unwrap()[2],0.);
            assert_eq!(study.evaluations(),spent);
        }
    });
}

#[test]
fn g1_spatial_lbfgs_recovers_heating_without_replacing_one_sensor_by_node_targets() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let mut rows=observations(cx,&m);
        let spec=objective(rows.clone(),false);
        let mut truth=factory(&boundary,&k,&c)(cx,sources([2.,5.,0.])).unwrap();
        for (index,time) in [1.,2.].into_iter().enumerate() {
            let frame=truth.advance_event_step(cx,time,time-truth.conduction().time_s(),spec.step).unwrap();
            for row in rows.iter_mut().filter(|r|r.endpoint==index) {
                let solution=frame.report.solution();
                let field=match row.observable { ScanObservable::Temperature=>&solution.temperature,
                    ScanObservable::SpecificEnthalpy=>&solution.specific_enthalpy_j_kg,
                    ScanObservable::LiquidMassFraction=>&solution.liquid_mass_fraction };
                row.target=row.readout.evaluate(cx,field).unwrap();
            }
        }
        for short in [false,true] {
            let mut study=ScanPowerStudy::new(cx,sources([1.,1.,0.]),controls(),objective(rows.clone(),short),
                factory(&boundary,&k,&c),1).unwrap();
            if short { study.set_checkpointing(cx,checkpoints()).unwrap(); }
            let mut state=study.begin_optimizer(cx,6).unwrap(); let initial=state.x.clone();
            assert_eq!(study.advance_optimizer(cx,&mut state,&StopRule::GradNorm(1e-9),80).unwrap().reason,StopReason::Budget);
            assert_eq!(state.x,initial);
            study.extend_evaluation_budget(160).unwrap();
            study.advance_optimizer(cx,&mut state,&StopRule::GradNorm(1e-9),100).unwrap();
            let fitted=study.evaluate(cx,&state.x).unwrap();
            assert!(fitted.physical.value<1e-12);
            close(fitted.powers_w[0],2.); close(fitted.powers_w[1],5.);
            assert!(fitted.physical.final_step.solution().specific_enthalpy_j_kg
                .windows(2).any(|x| (x[0]-x[1]).abs()>1e-3));
        }
    });
}

#[test]
fn g3_spatial_sqp_keeps_point_peak_limits_zero_power_gradients_and_replay() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let mut outputs=Vec::new();
        for short in [false,true] {
            let mut study=ScanPowerStudy::new_joint_with_power_coordinates(cx,sources([0.,0.,0.]),controls(),vec![],
                ScanPowerCoordinates::NonnegativeLinear,objective(observations(cx,&m),short),
                factory(&boundary,&k,&c),32).unwrap();
            if short {study.set_checkpointing(cx,checkpoints()).unwrap();}
            let mut study=ConstrainedScanStudy::new_with_peak_bounds(cx,study,vec![],vec![],vec![ScanStateBound {
                response:ScanResponse { endpoint:0,vertex:0,observable:ScanObservable::SpecificEnthalpy },
                lower:None,upper:Some(158.),scale:10. }],vec![ScanPeakBound {
                response:ScanPeakResponse {first_endpoint:0,last_endpoint:1,vertices:(0..8).collect(),
                    observable:ScanObservable::LiquidMassFraction,smoothing_scale:0.01},
                maximum:0.8,scale:0.1}],ScanResponseBudget {max_responses:2,max_gradient_entries:64,
                    max_transpose_solves:5},16,8).unwrap();
            let accepted=study.point().to_vec();
            let (actual,limits)=study.evaluate_at(cx,&[0.,0.]).unwrap();
            assert!(actual.gradient.iter().any(|g|g.abs()>1e-5));
            assert!(study.constraints_at(cx,&[0.,0.]).is_err());
            for j in 0..2 {
                let mut plus=vec![0.,0.]; plus[j]=1e-6;
                let (p,pc)=study.evaluate_at(cx,&plus).unwrap();
                close(actual.gradient[j],(p.physical.value-actual.physical.value)/1e-6);
                for r in 0..limits.residuals.len() {
                    close(limits.jacobian[r*2+j],(pc.residuals[r]-limits.residuals[r])/1e-6);
                }
            }
            let spent=study.evaluations();
            assert!(study.evaluate_at(cx,&[-1.,0.]).is_err());
            assert_eq!(study.point(),accepted.as_slice());
            assert_eq!(study.evaluations(),spent);
            outputs.push((actual.gradient,limits.residuals,limits.jacobian));
        }
        assert_eq!(outputs[0],outputs[1]);
    });
}

#[test]
fn g0_foreign_sensor_mesh_and_cancellation_never_publish_a_design() {
    let m=mesh(); let foreign=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|gate,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let mut study=ScanPowerStudy::new(cx,sources([1.,1.,0.]),controls(),objective(observations(cx,&foreign),false),
            factory(&boundary,&k,&c),3).unwrap();
        let initial=study.initial_coordinates().to_vec();
        assert!(study.evaluate(cx,&initial).is_err());
        assert_eq!(study.evaluations(),1);
        assert_eq!(study.initial_coordinates(),initial.as_slice());
        gate.request();
        assert!(study.evaluate(cx,&initial).is_err());
        assert_eq!(study.evaluations(),1);
    });
}

#[test]
fn g3_spatial_vertex_reduction_preserves_the_existing_default_observation_api() {
    let m=mesh(); let k=ConductivityModel::isotropic_declared(2.).unwrap(); let c=curve();
    with_cx(|_,cx| {
        let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&m)
            .remainder("surface",ThermalBc::adiabatic()).unwrap(),&Environment::earth_lab(),vec![],0.).unwrap();
        let rows=vec![ScanObservation {endpoint:1,vertex:3,observable:ScanObservable::SpecificEnthalpy,target:157.,scale:10.}];
        let native_spec: ScanPowerObjective=objective(rows,false);
        let spatial_spec=objective(vec![ScanSpatialObservation {endpoint:1,
            readout:ScanReadout::vertex(cx,&m,3).unwrap(),observable:ScanObservable::SpecificEnthalpy,target:157.,scale:10.}],false);
        let mut vertex=ScanPowerStudy::new(cx,sources([2.,3.,0.]),controls(),native_spec,factory(&boundary,&k,&c),1).unwrap();
        let mut spatial=ScanPowerStudy::new(cx,sources([2.,3.,0.]),controls(),spatial_spec,factory(&boundary,&k,&c),1).unwrap();
        let point=vertex.initial_coordinates().to_vec();
        let a=vertex.evaluate(cx,&point).unwrap(); let b=spatial.evaluate(cx,&point).unwrap();
        assert_eq!(a.physical.value,b.physical.value); assert_eq!(a.gradient,b.gradient);
        assert_eq!(a.physical.predictions,b.physical.predictions);
    });
}

#[path = "spatial_design/robust.rs"]
mod robust;
