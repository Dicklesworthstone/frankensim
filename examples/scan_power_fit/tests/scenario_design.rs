//! Finite-case objectives and constraints use the original thermal producer.
use fs_ascent::scan_design::{ScanPowerControl, ScanPowerCoordinates, ScanPowerObjective, ScanPowerStudy};
use fs_ascent::scan_design::constraints::{ScanResponse, ScanStateBound};
use fs_ascent::scan_design::constraints::value::{ValueConstrainedScanStudy, ScanValueConstraintBudget};
use fs_ascent::scan_design::constraints::value::scenarios::*;
use fs_ascent::scan_design::value::ScanValueResponseBudget;
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun, EnthalpyHistoryMaterials, EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::{ScanSensitivityError,
    objective::{ScanObservable, ScanObservation, ScanObjectiveSampling, ScanObjectiveBudget}};
use fs_conduction::bc::history::march::moving::{MovingHeatSource, MovingHeatBudget};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionError, ConductionMesh, ConductivityModel, ScalarField,
    ThermalBc, ThermalBoundaryBuilder, LinearConfig, COMPONENT_POWER_DIMS};
use fs_exec::{Cx, Budget, CancelGate, ExecMode, StreamKey};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{scenario::Environment, signal::{TimeSignal, Interp}};
use fs_solver::NewtonKrylovConfig;

type Factory<'a> = Box<dyn Fn(&Cx<'_>, Vec<MovingHeatSource>)
    -> Result<ScanningEnthalpyRun<'a, 'a, 'a>, ScanSensitivityError> + Sync + 'a>;
fn with_cx<T>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> T) -> T {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena|
        f(&gate, &Cx::new(&gate, arena, StreamKey { seed: 92, kernel_id: 51, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic)))
}
fn fixture() -> (ConductionMesh, ConductivityModel, EquilibriumEnthalpyPhaseCurve) {
    let (complex, points) = fs_conduction::fixtures::box_grid([1,1,1], [1.;3]);
    (ConductionMesh::new(complex, points).unwrap(), ConductivityModel::isotropic_declared(2.).unwrap(),
        EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([92;32]),
            [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
                .map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h, temperature_k:t,
                    liquid_mass_fraction:f, bulk_density_kg_m3:1. }).collect()).unwrap())
}
fn boundary<'a>(cx: &Cx<'_>, mesh: &'a ConductionMesh) -> ThermalBoundaryHistory<'a> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .remainder("surface", ThermalBc::adiabatic()).unwrap(), &Environment::earth_lab(), vec![], 0.).unwrap()
}
fn case<'a>(cx: &Cx<'_>, b: &'a ThermalBoundaryHistory<'a>, k: &'a ConductivityModel,
    c: &'a EquilibriumEnthalpyPhaseCurve, name: &str, h: f64, weight: f64,
    reference: f64, allowance: usize, fail: bool) -> ScanScenario<Factory<'a>>
{
    let source = MovingHeatSource::new(std::array::from_fn(|_| TimeSignal::Constant(
        QtyAny::new(0.5, Dims([1,0,0,0,0,0])))), [2.;3],
        TimeSignal::Table { times:vec![0.,1.,2.], values:vec![0.;3], dims:COMPONENT_POWER_DIMS,
            interp:Interp::Hold }, 0., MovingHeatBudget { max_vertices:8, max_elements:6,
                min_active_vertices:1, power_relative_tolerance:1e-11 }).unwrap();
    let factory: Factory<'a> = Box::new(move |cx, sources| {
        if fail { return Err(ConductionError::Config { parameter:"case test", what:"declared failure".into() }.into()); }
        let run = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem { boundary:b, material:k,
            element_materials:None, interfaces:None, source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
            EnthalpyHistoryMaterials { table:&[ReferenceEnthalpyMaterial { curve:c, reference_density_kg_m3:1. }],
                of_element:&[0;6], phase_conductivity:None }, &[h;8],
            EnthalpyHistoryBudget { spatial:EnthalpyBudget {max_vertices:8,max_elements:6}, max_steps:2 })?;
        Ok(ScanningEnthalpyRun::new(cx, run, sources, 1, None)?)
    });
    let study = ScanPowerStudy::new_joint_with_power_coordinates(cx, vec![source],
        vec![ScanPowerControl {source:0,value:0,reference_w:reference},
             ScanPowerControl {source:0,value:1,reference_w:reference}], vec![],
        ScanPowerCoordinates::NonnegativeLinear, ScanPowerObjective {
            endpoints_s:vec![1.,2.], observations:[152.,157.].into_iter().enumerate().map(|(endpoint,target)|
                ScanObservation {endpoint,vertex:0,observable:ScanObservable::SpecificEnthalpy,target,scale:10.}).collect(),
            sampling:ScanObjectiveSampling::EventLeft,
            step:EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,
                linear_restart:8,max_linear_cycles:16,forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
                max_newton_iterations:64,energy_tolerance_j:1e-8},
            transpose:LinearConfig {tolerance:1e-10,max_iterations:1000,restart:8},
            budget:ScanObjectiveBudget {max_steps:2,max_observations:2,max_retained_vertices:0,
                max_retained_source_vertices:0,max_schedule_controls:16,max_retained_radiation_entries:0},
            source_reference_time_s:0. }, factory, allowance).unwrap();
    ScanScenario {name:name.into(),weight,study:ValueConstrainedScanStudy::new(cx,study,vec![],vec![],
        vec![ScanStateBound { response:ScanResponse {endpoint:1,vertex:0,observable:ScanObservable::SpecificEnthalpy},
            lower:None,upper:Some(155.),scale:10.}],vec![],
        ScanValueConstraintBudget {max_rows:3,max_energy_weights:0,
            responses:ScanValueResponseBudget {max_responses:1,max_peak_samples:0}}).unwrap()}
}
fn budget() -> ScanScenarioBudget {
    ScanScenarioBudget {max_cases:2,max_constraint_rows:6,max_retained_vertices:16,max_name_bytes:64}
}

#[test]
fn g3_case_aggregation_preserves_each_original_prediction_and_violation() {
    let (mesh,k,c)=fixture();
    with_cx(|_,cx| {
        let b=boundary(cx,&mesh);
        for policy in [ScanScenarioObjective::WeightedMean,ScanScenarioObjective::WorstCase] {
            let mut a=case(cx,&b,&k,&c,"warm",150.,1.,1.,4,false);
            let mut d=case(cx,&b,&k,&c,"cold",140.,3.,1.,4,false);
            let expected_a=a.study.evaluate_at(cx,&[2.,5.]).unwrap();
            let expected_d=d.study.evaluate_at(cx,&[2.,5.]).unwrap();
            let mut cases=ScanScenarioStudy::new(cx,vec![a,d],policy,budget()).unwrap();
            let result=cases.evaluate_at(cx,&[2.,5.]).unwrap();
            assert_eq!(result.cases[0].physical.physical.predictions,expected_a.0.physical.predictions);
            assert_eq!(result.cases[1].physical.physical.predictions,expected_d.0.physical.predictions);
            assert_eq!(result.cases[0].constraints,expected_a.1);
            assert_eq!(result.cases[1].constraints,expected_d.1);
            assert_eq!(result.constraint_offsets,vec![0,3,6]);
            assert_eq!(&result.residuals[..3],result.cases[0].constraints.residuals.as_slice());
            assert!(result.residuals[2]>0.19 && result.residuals[5] < -0.79);
            let expected=if policy==ScanScenarioObjective::WeightedMean {0.75} else {1.0};
            assert!((result.objective-expected).abs()<1e-9);
            assert_eq!(result.worst_objective_case,1);
            assert_eq!(cases.evaluations(),vec![2,2]);
            assert_eq!(cases.sources_at(cx,&[2.,5.]).unwrap()[0].absorbed_power_values(3).unwrap(),vec![2.,5.,0.]);
        }
    });
}

#[test]
fn g0_common_dimensions_do_not_admit_different_source_maps_or_bad_case_sets() {
    let (mesh,k,c)=fixture();
    with_cx(|_,cx| {
        let b=boundary(cx,&mesh);
        for (second_name,weight,reference) in [("warm",1.,1.),("cold",0.,1.),("cold",1.,2.)] {
            let cases=vec![case(cx,&b,&k,&c,"warm",150.,1.,1.,1,false),
                case(cx,&b,&k,&c,second_name,140.,weight,reference,1,false)];
            assert!(ScanScenarioStudy::new(cx,cases,ScanScenarioObjective::WeightedMean,budget()).is_err());
        }
        let cases=vec![case(cx,&b,&k,&c,"warm",150.,1.,1.,1,false),case(cx,&b,&k,&c,"cold",140.,1.,1.,1,false)];
        assert!(ScanScenarioStudy::new(cx,cases,ScanScenarioObjective::WorstCase,
            ScanScenarioBudget {max_constraint_rows:5,..budget()}).is_err());
    });
}

#[test]
fn g4_partial_physics_failure_and_known_budget_exhaustion_preserve_work() {
    let (mesh,k,c)=fixture();
    with_cx(|gate,cx| {
        let b=boundary(cx,&mesh);
        let mut cases=ScanScenarioStudy::new(cx,vec![case(cx,&b,&k,&c,"good",150.,1.,1.,1,false),
            case(cx,&b,&k,&c,"fails",140.,1.,1.,1,true)],ScanScenarioObjective::WeightedMean,budget()).unwrap();
        assert!(matches!(cases.evaluate_at(cx,&[2.,5.]),Err(ScanScenarioError::Physical {case:1,..})));
        assert_eq!(cases.evaluations(),vec![1,1]);
        assert!(cases.evaluate_at(cx,&[2.,5.]).is_err());
        assert_eq!(cases.evaluations(),vec![1,1]);
        assert!(cases.extend_evaluation_budgets(&[2,0]).is_err());
        assert_eq!(cases.evaluation_limits(),vec![1,1]);
        cases.extend_evaluation_budgets(&[2,2]).unwrap();
        assert_eq!(cases.remaining_evaluations(),1);
        gate.request();
        assert!(matches!(cases.evaluate_at(cx,&[2.,5.]),Err(ScanScenarioError::Cancelled)));
        assert_eq!(cases.evaluations(),vec![1,1]);
    });
}

#[test]
#[cfg(not(target_arch="wasm32"))]
fn g5_shared_case_batches_match_serial_reports_and_settle_each_budget() {
    let (mesh,k,c)=fixture();
    with_cx(|_,cx| {
        let b=boundary(cx,&mesh);
        let build=|| ScanScenarioStudy::new(cx,vec![case(cx,&b,&k,&c,"warm",150.,1.,1.,4,false),
            case(cx,&b,&k,&c,"cold",140.,1.,1.,4,false)],ScanScenarioObjective::WorstCase,budget()).unwrap();
        let mut serial=build();
        let expected=serial.evaluate_at(cx,&[2.,5.]).unwrap();
        let mut parallel=build();
        {
            let batch=parallel.value_batch(cx,3).unwrap();
            std::thread::scope(|scope| {
                let shared=&batch;
                let jobs:Vec<_>=(0..2).map(|_|scope.spawn(move ||with_cx(|_,worker| {
                    let value=shared.evaluate_at(worker,&[2.,5.]).unwrap();
                    (value.objective,value.residuals)
                }))).collect();
                for job in jobs {
                    let (value,rows)=job.join().unwrap();
                    assert_eq!(value.to_bits(),expected.objective.to_bits());
                    assert_eq!(rows,expected.residuals);
                }
            });
        }
        assert_eq!(parallel.evaluations(),vec![2,2]);
        assert_eq!(parallel.remaining_evaluations(),2);
        let replay=parallel.evaluate_at(cx,&[2.,5.]).unwrap();
        assert_eq!(replay.objective.to_bits(),expected.objective.to_bits());
    });
}
