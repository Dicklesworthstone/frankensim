//! Checkpointed gradients through the original thermal and radiation solvers.
mod support;
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory,ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun,EnthalpyHistoryMaterials,EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::{ScanningEnthalpyRun,ScanRadiation};
use fs_conduction::bc::history::march::moving::{MovingHeatBudget,MovingHeatSource};
use fs_conduction::bc::history::march::radiation::RadiationRegionHistory;
use fs_conduction::transient::enthalpy::{EnthalpyBudget,EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionMesh,ConductivityModel,ConductivityTable,ScalarField,ThermalBc,
    ThermalBoundaryBuilder,LinearConfig,AmbientRadiationConfig,SurfaceEmissivity,
    COMPONENT_POWER_DIMS,TEMPERATURE_DIMS,EMISSIVITY_DIMS,SURFACE_EMISSIVITY_PROPERTY};
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve,EnthalpyPhaseKnot};
use fs_matdb::{ClaimSet,PropertyClaim,PropertyKey,PropertyValue,UncertaintyModel,InterpolationPolicy,
    Provenance,MaterialCard,MaterialStateId,SelectionPolicy};
use fs_qty::{Dims,QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{TimeSignal,Interp};
use fs_solver::NewtonKrylovConfig;
use support::{with_cx,with_cancelled_cx};
fn constant(v:f64,d:Dims)->TimeSignal{TimeSignal::Constant(QtyAny::new(v,d))}
fn mesh()->ConductionMesh{ConductionMesh::new(TetComplex::from_tets(4,vec![[0,1,2,3]]),
    vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]).unwrap()}
fn chart()->EquilibriumEnthalpyPhaseCurve{
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([74;32]),
        [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
            .map(|(h,t,f)|EnthalpyPhaseKnot{specific_enthalpy_j_kg:h,temperature_k:t,
                liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect()).unwrap()
}
fn config()->EnthalpyStepConfig{EnthalpyStepConfig{newton:NewtonKrylovConfig{
    absolute_tolerance:1e-12,relative_tolerance:1e-12,linear_restart:4,max_linear_cycles:16,
    forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
    max_newton_iterations:64,energy_tolerance_j:1e-9}}
fn linear()->LinearConfig{LinearConfig{tolerance:1e-11,max_iterations:1000,restart:4}}
fn heater(p:[f64;7])->MovingHeatSource{MovingHeatSource::new(
    [p[0],p[1],p[2]].map(|v|constant(v,Dims([1,0,0,0,0,0]))),[p[3],p[4],p[5]],
    constant(p[6],COMPONENT_POWER_DIMS),0.,MovingHeatBudget{max_vertices:4,max_elements:1,
        min_active_vertices:2,power_relative_tolerance:1e-12}).unwrap()}
fn history<'m>(cx:&Cx<'_>,m:&'m ConductionMesh)->ThermalBoundaryHistory<'m>{
    ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(m).remainder("surface",ThermalBc::adiabatic()).unwrap(),
        &Environment::earth_lab(),vec![],0.).unwrap()
}
fn radiation()->ScanRadiation{
    let mut claims=ClaimSet::new();claims.insert_claim(PropertyClaim{
        key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar{value:0.2,dims:EMISSIVITY_DIMS},
        validity:ValidityDomain::unconstrained().with("T",280.,380.),uncertainty:UncertaintyModel::Unstated,
        interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
        provenance:Provenance{source:"synthetic scan gradient fixture".into(),license:"test-only".into(),artifact:None}}).unwrap();
    let card=MaterialCard::assemble(MaterialStateId{chemistry:"gray fixture".into(),phase:"solid".into(),
        process:"synthetic".into(),revision:0},claims,vec![]).unwrap();
    ScanRadiation{patches:vec![RadiationRegionHistory{region:"surface".into(),
        emissivity:SurfaceEmissivity::from_card("surface",&card,300.,SelectionPolicy::SingleClaimOnly).unwrap(),
        ambient_temperature:constant(293.,TEMPERATURE_DIMS)}],config:AmbientRadiationConfig::default()}
}
fn run<'h,'m,'c>(cx:&Cx<'_>,history:&'h ThermalBoundaryHistory<'m>,material:&'h ConductivityModel,
    curve:&'c EquilibriumEnthalpyPhaseCurve,old:&[f64],sources:Vec<MovingHeatSource>,rad:bool)
    ->ScanningEnthalpyRun<'h,'m,'c>{
    let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem{boundary:history,material,
        element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
        EnthalpyHistoryMaterials{table:&[ReferenceEnthalpyMaterial{curve,reference_density_kg_m3:1.}],
            of_element:&[0],phase_conductivity:None},old,
        EnthalpyHistoryBudget{spatial:EnthalpyBudget{max_vertices:4,max_elements:1},max_steps:8}).unwrap();
    ScanningEnthalpyRun::new(cx,base,sources,4,rad.then(radiation)).unwrap()
}

// Checkpoint policy must not change the original thermal objective or adjoint.
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
    evaluate_scan_objective, ScanObjectiveBudget, ScanObjectiveEvaluation,
    ScanObjectiveSampling, ScanObservable, ScanObservation,
    checkpoint::{evaluate_scan_objective_checkpointed, ScanCheckpointBudget},
};

fn budget(steps: usize, retained: usize) -> ScanObjectiveBudget {
    ScanObjectiveBudget { max_steps: steps, max_observations: 16,
        max_retained_vertices: 4*retained, max_retained_source_vertices: 8*retained,
        max_schedule_controls: 32, max_retained_radiation_entries: 1024*retained }
}
fn checkpoint_budget(steps: usize, segment: usize) -> ScanCheckpointBudget {
    ScanCheckpointBudget { segment_steps: segment,
        max_checkpoint_vertices: (steps.div_ceil(segment)+1)*4, max_forward_solves: 2*steps }
}
fn observations(last: usize) -> Vec<ScanObservation> {
    vec![ScanObservation {endpoint: 0, vertex: 0, observable: ScanObservable::Temperature, target: 303., scale: 2.},
        ScanObservation {endpoint: last/2, vertex: 1, observable: ScanObservable::SpecificEnthalpy, target: 40., scale: 10.},
        ScanObservation {endpoint: last, vertex: 2, observable: ScanObservable::Temperature, target: 307., scale: 3.}]
}
fn same(a: &ScanObjectiveEvaluation, b: &ScanObjectiveEvaluation) {
    assert_eq!(a.value.to_bits(), b.value.to_bits());
    assert_eq!(a.predictions, b.predictions);
    assert_eq!(a.standardized_residuals, b.standardized_residuals);
    assert_eq!(a.sources, b.sources);
    assert_eq!(a.initial_specific_enthalpy, b.initial_specific_enthalpy);
    assert_eq!(a.transpose_iterations, b.transpose_iterations);
    assert_eq!(a.max_transpose_relative_residual, b.max_transpose_relative_residual);
    assert_eq!(a.final_step.solution().specific_enthalpy_j_kg, b.final_step.solution().specific_enthalpy_j_kg);
    assert_eq!(a.final_step.solution().temperature, b.final_step.solution().temperature);
    assert_eq!(a.final_step.solution().liquid_mass_fraction, b.final_step.solution().liquid_mass_fraction);
}

#[test]
fn g5_segments_preserve_full_tape_gradients_across_phase_branches_and_radiation() {
    let mesh = mesh(); let curve = chart();
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(280.,1.),(380.,3.)]).unwrap());
    with_cx(|cx| {
        let hist = history(cx, &mesh);
        let sources = vec![heater([0.2,0.25,0.3,1.4,1.3,1.2,3.]), heater([0.6,0.1,0.2,1.5,1.4,1.3,1.5])];
        let times = [0.001,0.002,0.004,0.006,0.009]; let obs = observations(4);
        for old in [20.,150.,500.] {
            for rad in [false,true] {
                let make = || run(cx,&hist,&material,&curve,&[old;4],sources.clone(),rad);
                let full = evaluate_scan_objective(cx,make(),&times,&obs,
                    ScanObjectiveSampling::RightEndpoint,config(),linear(),budget(5,5)).unwrap();
                for segment in [1,2,5] {
                    let short = evaluate_scan_objective_checkpointed(cx,make(),&times,&obs,
                        ScanObjectiveSampling::RightEndpoint,config(),linear(),budget(5,segment),
                        checkpoint_budget(5,segment)).unwrap();
                    same(&full,&short.objective);
                    assert_eq!(short.work.original_forward_solves,5);
                    assert_eq!(short.work.recomputed_forward_solves,5);
                    assert_eq!(short.work.peak_retained_linearizations,segment);
                    assert_eq!(short.work.checkpoint_states,5usize.div_ceil(segment)+1);
                }
            }
        }
    });
}

#[test]
fn g1_short_tape_admission_runs_a_trajectory_the_full_tape_refuses() {
    let mesh = mesh(); let curve = chart(); let material = ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let hist = history(cx,&mesh);
        let times: Vec<_> = (1..=32).map(|i| f64::from(i)*0.00025).collect();
        let obs = observations(31);
        let make = || {
            let mut run = run(cx,&hist,&material,&curve,&[150.;4],vec![heater([0.2,0.25,0.3,1.4,1.3,1.2,3.])],false);
            run.extend_step_budget(32).unwrap(); run
        };
        assert!(evaluate_scan_objective(cx,make(),&times,&obs,ScanObjectiveSampling::RightEndpoint,
            config(),linear(),budget(32,4)).is_err());
        let short = evaluate_scan_objective_checkpointed(cx,make(),&times,&obs,
            ScanObjectiveSampling::RightEndpoint,config(),linear(),budget(32,4),checkpoint_budget(32,4)).unwrap();
        assert_eq!(short.work.peak_retained_linearizations,4);
        assert_eq!(short.work.checkpoint_states,9);
        assert_eq!(short.objective.endpoint_count,32);
        assert!(short.objective.sources[0].absorbed_power_values_w[0] != 0.);
    });
}

#[test]
fn g5_event_sides_and_an_already_advanced_initial_state_survive_replay() {
    let mesh = mesh(); let curve = chart(); let material = ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let hist = history(cx,&mesh); let p = [0.2,0.25,0.3,1.4,1.3,1.2,3.];
        let held = MovingHeatSource::new([p[0],p[1],p[2]].map(|v|constant(v,Dims([1,0,0,0,0,0]))),
            [p[3],p[4],p[5]],TimeSignal::Table {times:vec![0.,0.004,0.008],values:vec![3.,0.,1.],
                dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},0.,heater(p).budget()).unwrap();
        let times = [0.002,0.004,0.006,0.008]; let obs = observations(3);
        let make = || {
            let mut r = run(cx,&hist,&material,&curve,&[150.;4],vec![held.clone()],false);
            r.advance_to(cx,0.001,config()).unwrap(); r
        };
        for side in [ScanObjectiveSampling::RightEndpoint,ScanObjectiveSampling::EventLeft] {
            let full = evaluate_scan_objective(cx,make(),&times,&obs,side,config(),linear(),budget(4,4)).unwrap();
            let short = evaluate_scan_objective_checkpointed(cx,make(),&times,&obs,side,
                config(),linear(),budget(4,2),checkpoint_budget(4,2)).unwrap();
            same(&full,&short.objective);
        }
        assert!(evaluate_scan_objective_checkpointed(cx,make(),&[0.003,0.008],&observations(1),
            ScanObjectiveSampling::EventLeft,config(),linear(),budget(2,1),checkpoint_budget(2,1)).is_err());
    });
}

#[test]
fn g0_checkpoint_work_memory_and_cancellation_refusals_return_no_objective() {
    let mesh = mesh(); let curve = chart(); let material = ConductivityModel::isotropic_declared(2.).unwrap();
    let hist = with_cx(|cx|history(cx,&mesh));
    let make = |cx:&Cx<'_>|run(cx,&hist,&material,&curve,&[150.;4],vec![heater([0.2,0.25,0.3,1.4,1.3,1.2,3.])],false);
    with_cx(|cx| {
        let good = checkpoint_budget(2,1);
        for bad in [ScanCheckpointBudget {segment_steps:0,..good},
            ScanCheckpointBudget {max_checkpoint_vertices:11,..good},
            ScanCheckpointBudget {max_forward_solves:3,..good}]
        {
            assert!(evaluate_scan_objective_checkpointed(cx,make(cx),&[0.001,0.002],&observations(1),
                ScanObjectiveSampling::RightEndpoint,config(),linear(),budget(2,1),bad).is_err());
        }
    });
    let ready = with_cx(make);
    with_cancelled_cx(|cx|assert!(evaluate_scan_objective_checkpointed(cx,ready,&[0.001,0.002],&observations(1),
        ScanObjectiveSampling::RightEndpoint,config(),linear(),budget(2,1),checkpoint_budget(2,1)).is_err()));
}

#[test]
fn g5_point_and_peak_constraints_share_replay_without_changing_their_gradients() {
    use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
        checkpoint::evaluate_scan_objective_checkpointed_with_peaks,
        responses::{ScanResponse,ScanResponseBudget,peaks::{ScanPeakResponse,evaluate_scan_objective_with_peaks}},
    };
    let mesh = mesh(); let curve = chart();
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(280.,1.),(380.,3.)]).unwrap());
    with_cx(|cx| {
        let hist=history(cx,&mesh);let times=[0.001,0.002,0.004,0.006,0.009];let obs=observations(4);
        let rows=[ScanResponse {endpoint:1,vertex:0,observable:ScanObservable::Temperature},
            ScanResponse {endpoint:3,vertex:1,observable:ScanObservable::SpecificEnthalpy}];
        let peaks=[ScanPeakResponse {first_endpoint:0,last_endpoint:4,vertices:vec![0,1,2,3],
            observable:ScanObservable::Temperature,smoothing_scale:0.2},
            ScanPeakResponse {first_endpoint:1,last_endpoint:3,vertices:vec![0,2],
                observable:ScanObservable::SpecificEnthalpy,smoothing_scale:0.5}];
        let extra=ScanResponseBudget {max_responses:4,max_gradient_entries:128,max_transpose_solves:20};
        for initial in [20.,150.] {
            for rad in [false,true] {
                let make=||run(cx,&hist,&material,&curve,&[initial;4],
                    vec![heater([0.2,0.25,0.3,1.4,1.3,1.2,3.])],rad);
                let (objective,full_rows,full_peaks)=evaluate_scan_objective_with_peaks(cx,make(),&times,&obs,
                    ScanObjectiveSampling::RightEndpoint,config(),linear(),budget(5,5),&rows,&peaks,extra,26).unwrap();
                for segment in [1,2,5] {
                    let replay=evaluate_scan_objective_checkpointed_with_peaks(cx,make(),&times,&obs,
                        ScanObjectiveSampling::RightEndpoint,config(),linear(),budget(5,segment),&rows,&peaks,
                        extra,26,checkpoint_budget(5,segment)).unwrap();
                    same(&objective,&replay.objective);
                    assert_eq!(replay.work.recomputed_forward_solves,5);
                    for (a,b) in full_rows.iter().zip(&replay.responses) {
                        assert_eq!(a.value.to_bits(),b.value.to_bits());
                        assert_eq!(a.sources,b.sources);
                        assert_eq!(a.initial_specific_enthalpy,b.initial_specific_enthalpy);
                        assert_eq!(a.transpose_iterations,b.transpose_iterations);
                    }
                    for (a,b) in full_peaks.iter().zip(&replay.peaks) {
                        assert_eq!(a.value.to_bits(),b.value.to_bits());
                        assert_eq!(a.sampled_maximum,b.sampled_maximum);
                        assert_eq!(a.smoothing_gap,b.smoothing_gap);
                        assert_eq!(a.maximum_smoothing_gap,b.maximum_smoothing_gap);
                        assert_eq!(a.sources,b.sources);
                        assert_eq!(a.initial_specific_enthalpy,b.initial_specific_enthalpy);
                        assert_eq!(a.transpose_iterations,b.transpose_iterations);
                    }
                }
                assert!(evaluate_scan_objective_checkpointed_with_peaks(cx,make(),&times,&obs,
                    ScanObjectiveSampling::RightEndpoint,config(),linear(),budget(5,2),&rows,&peaks,
                    extra,25,checkpoint_budget(5,2)).is_err());
            }
        }
    });
}
