//! One heating schedule across explicit starting-enthalpy cases.
//! These are synthetic uniform P1 specimens, not a sampled probability law.
use scan_power_fit::scenarios::{ScenarioPopulationStudy, ScenarioPopulationReplay};
use fs_ascent::scan_design::{ScanPowerControl, ScanPowerCoordinates, ScanPowerObjective, ScanPowerStudy};
use fs_ascent::scan_design::constraints::{ScanEnergyCap, ScanStateBound, ScanPeakBound, ScanResponse, ScanPeakResponse};
use fs_ascent::scan_design::constraints::value::{ValueConstrainedScanStudy, ScanValueConstraintBudget};
use fs_ascent::scan_design::constraints::value::scenarios::{ScanScenario, ScanScenarioStudy,
    ScanScenarioBudget, ScanScenarioObjective};
use fs_ascent::scan_design::value::ScanValueResponseBudget;
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun, EnthalpyHistoryMaterials, EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::{ScanSensitivityError,
    objective::{ScanObservable, ScanObjectiveSampling, ScanObjectiveBudget}};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::spatial::{ScanReadout, ScanSpatialObservation};
use fs_conduction::bc::history::march::moving::{MovingHeatSource, MovingHeatBudget};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionMesh, ConductivityModel, ScalarField, ThermalBc,
    ThermalBoundaryBuilder, LinearConfig, COMPONENT_POWER_DIMS};
use fs_dfo::differential_evolution::{DifferentialEvolutionConfig, DifferentialEvolutionStop};
use fs_exec::{Cx, Budget, CancelGate, ExecMode, StreamKey, TilePool};
use fs_alloc::OperationMemoryLease;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{scenario::Environment, signal::{TimeSignal, Interp}};
use fs_solver::NewtonKrylovConfig;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const TARGET:f64=1e-8;
const TOLERANCE:f64=1e-9;
const USAGE:&str="scenario_scan_fit [--cases h0,h1,...] [--weights w0,w1,...] [--worst-case] [--workers 1..128] [energy_cap_J [peak_fraction [search_attempts]]]";
#[derive(Clone)]
struct Options {enthalpies:Vec<f64>,weights:Vec<f64>,policy:ScanScenarioObjective,
    workers:Option<usize>,energy_j:f64,peak_fraction:f64,attempts:usize}
fn list(s:&str)->Result<Vec<f64>> {s.split(',').map(|x|x.parse().map_err(Into::into)).collect()}
fn options(args:&[String])->Result<Options> {
    let (mut hs,mut weights,mut workers)=(None,None,None);
    let mut worst=false;
    let mut values=Vec::new();
    let mut i=0;
    while i<args.len() {
        match args[i].as_str() {
            "--cases" if hs.is_none()=>{i+=1;hs=Some(list(args.get(i).ok_or(USAGE)?)?);},
            "--weights" if weights.is_none()=>{i+=1;weights=Some(list(args.get(i).ok_or(USAGE)?)?);},
            "--workers" if workers.is_none()=>{i+=1;workers=Some(args.get(i).ok_or(USAGE)?.parse::<usize>()?);},
            "--worst-case" if !worst=>worst=true,
            value if !value.starts_with("--")=>values.push(value),
            _=>return Err(USAGE.into()),
        }
        i+=1;
    }
    let enthalpies=hs.unwrap_or_else(||vec![80.,90.]);
    let weights=weights.unwrap_or_else(||vec![1.;enthalpies.len()]);
    if values.len()>3 || enthalpies.is_empty() || enthalpies.len()>16 || weights.len()!=enthalpies.len()
        || enthalpies.iter().any(|v|!v.is_finite() || !(-100.0..=200.0).contains(v))
        || weights.iter().any(|v|!v.is_finite() || *v<=0.0)
        || workers.is_some_and(|n|!(1..=128).contains(&n)) {return Err(USAGE.into());}
    let out=Options {enthalpies,weights,workers,policy:if worst {ScanScenarioObjective::WorstCase} else {ScanScenarioObjective::WeightedMean},
        energy_j:values.first().map(|s|s.parse()).transpose()?.unwrap_or(60.),
        peak_fraction:values.get(1).map(|s|s.parse()).transpose()?.unwrap_or(0.5),
        attempts:values.get(2).map(|s|s.parse()).transpose()?.unwrap_or(2400)};
    if !out.energy_j.is_finite() || !(0.0..=100.0).contains(&out.energy_j)
        || !out.peak_fraction.is_finite() || !(0.0..=1.0).contains(&out.peak_fraction)
        || !(32..=10_000).contains(&out.attempts) {return Err(USAGE.into());}
    Ok(out)
}

fn fit(cx:&Cx<'_>,input:&Options,partitioned:bool)->Result<ScenarioPopulationReplay> {
    let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
    let mesh=ConductionMesh::new(complex,positions)?;
    let material=ConductivityModel::isotropic_declared(2.)?;
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([92;32]),
        [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
            .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
                liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect())?;
    let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface",ThermalBc::adiabatic())?,&Environment::earth_lab(),vec![],0.)?;
    let readout=ScanReadout::volume_mean(cx,&mesh,&[0,1,2,3,4,5],6)?;
    let observations:Vec<_>=[0.1,0.3].into_iter().enumerate().map(|(endpoint,target)|ScanSpatialObservation {
        endpoint,readout:readout.clone(),observable:ScanObservable::LiquidMassFraction,target,scale:0.1}).collect();
    let source=MovingHeatSource::new(std::array::from_fn(|_|TimeSignal::Constant(QtyAny::new(0.5,Dims([1,0,0,0,0,0])))),
        [2.;3],TimeSignal::Table {times:vec![0.,1.,2.],values:vec![0.;3],dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},
        0.,MovingHeatBudget {max_vertices:8,max_elements:6,min_active_vertices:1,power_relative_tolerance:1e-11})?;
    let initial_cap=if partitioned {7} else {input.attempts};
    let mut cases=Vec::new();
    for (index,(&h,&weight)) in input.enthalpies.iter().zip(&input.weights).enumerate() {
        let (b,k,c)=(&boundary,&material,&curve);
        let factory=move |cx:&Cx<'_>,sources|->std::result::Result<_,ScanSensitivityError> {
            let run=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:b,material:k,
                element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
                EnthalpyHistoryMaterials {table:&[ReferenceEnthalpyMaterial {curve:c,reference_density_kg_m3:1.}],
                    of_element:&[0;6],phase_conductivity:None},&[h;8],
                EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:8,max_elements:6},max_steps:2})?;
            Ok(ScanningEnthalpyRun::new(cx,run,sources,1,None)?)
        };
        let study=ScanPowerStudy::new_joint_with_power_coordinates(cx,vec![source.clone()],
            vec![ScanPowerControl {source:0,value:0,reference_w:1.},ScanPowerControl {source:0,value:1,reference_w:1.}],
            vec![],ScanPowerCoordinates::NonnegativeLinear,ScanPowerObjective {
                endpoints_s:vec![1.,2.],observations:observations.clone(),sampling:ScanObjectiveSampling::EventLeft,
                step:EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,
                    linear_restart:8,max_linear_cycles:16,forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
                    max_newton_iterations:64,energy_tolerance_j:1e-8},
                transpose:LinearConfig {tolerance:1e-10,max_iterations:1000,restart:8},
                budget:ScanObjectiveBudget {max_steps:2,max_observations:2,max_retained_vertices:0,
                    max_retained_source_vertices:0,max_schedule_controls:16,max_retained_radiation_entries:0},
                source_reference_time_s:0.},factory,initial_cap+1)?;
        let study=ValueConstrainedScanStudy::new(cx,study,vec![],
            vec![ScanEnergyCap {start_time_s:0.,end_time_s:2.,maximum_j:input.energy_j,scale_j:1.}],
            vec![ScanStateBound {response:ScanResponse {endpoint:0,vertex:0,observable:ScanObservable::Temperature},
                lower:None,upper:Some(311.),scale:1.}],
            vec![ScanPeakBound {response:ScanPeakResponse {first_endpoint:0,last_endpoint:1,vertices:(0..8).collect(),
                observable:ScanObservable::LiquidMassFraction,smoothing_scale:0.001},maximum:input.peak_fraction,scale:1.}],
            ScanValueConstraintBudget {max_rows:5,max_energy_weights:2,
                responses:ScanValueResponseBudget {max_responses:2,max_peak_samples:16}})?;
        cases.push(ScanScenario {name:format!("case-{index}"),weight,study});
    }
    let study=ScanScenarioStudy::new(cx,cases,input.policy,ScanScenarioBudget {
        max_cases:16,max_constraint_rows:80,max_retained_vertices:128,max_name_bytes:256})?;
    let mut driver=ScenarioPopulationStudy::new(cx,study,&[[0.,50.];2],DifferentialEvolutionConfig {
        population_size:24,differential_weight:0.8,crossover_probability:0.9,seed:17,
        max_evaluations:initial_cap,max_population_entries:48,target:Some(TARGET)},TOLERANCE)?;
    let pool=input.workers.map(|workers|TilePool::for_host(workers,17));
    let lease=OperationMemoryLease::bounded(64<<20);
    if partitioned {
        match &pool {Some(pool)=>{driver.advance_parallel(cx,pool,&lease,7,24)?;},None=>{driver.advance(cx,7)?;}}
        driver.extend_evaluation_budgets(input.attempts,&vec![input.attempts+1;input.enthalpies.len()])?;
    }
    loop {
        let additional=if partitioned {13} else {input.attempts};
        let report=match &pool {Some(pool)=>driver.advance_parallel(cx,pool,&lease,additional,24)?,
            None=>driver.advance(cx,additional)?};
        if report.remaining_physical_search_candidates==0 || report.search.stop!=DifferentialEvolutionStop::Paused {break;}
    }
    let result=driver.replay_best(cx)?;
    if result.physical_evaluations.iter().any(|&n|n!=result.search_evaluations+1) {
        return Err("physical cases did not retain their search work and individual replay".into());
    }
    if result.sources[0].absorbed_power_values(3)?[2].to_bits()!=0.0_f64.to_bits() {
        return Err("common source export changed the fixed shutdown".into());
    }
    Ok(result)
}
fn with_cx<T>(f:impl FnOnce(&CancelGate,&Cx<'_>)->Result<T>)->Result<T> {
    let gate=CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena|
        f(&gate,&Cx::new(&gate,arena,StreamKey {seed:92,kernel_id:51,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn main()->Result<()> {
    let input=options(&std::env::args().skip(1).collect::<Vec<_>>())?;
    let result=with_cx(|_,cx|fit(cx,&input,true))?;
    println!("# finite synthetic cases; no probability, continuum or optimality certificate");
    println!("# policy={:?}; feasible={}; aggregate_target_met={}; every_case_target_met={}; objective={:.17e}",
        input.policy,result.score.is_feasible(),result.aggregate_target_met(),result.every_case_target_met(),result.evaluation.objective);
    println!("# source_watts={:?}; search_attempts={}; physical_evaluations_by_case={:?}",
        result.coordinates,result.search_evaluations,result.physical_evaluations);
    println!("case,initial_h_j_kg,weight,loss,first_fraction,second_fraction,scaled_violation,sampled_peak,envelope");
    for (i,case) in result.evaluation.cases.iter().enumerate() {
        println!("{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            case.name,input.enthalpies[i],case.weight,case.physical.physical.value,
            case.physical.physical.predictions[0],case.physical.physical.predictions[1],
            case.constraints.max_scaled_violation(),case.constraints.sampled_peak_values[0],case.constraints.peak_values[0]);
    }
    if !result.score.is_feasible() {return Err("no feasible sampled design found; not an infeasibility proof".into());}
    if !result.every_case_target_met() {return Err("feasible compromise still misses individual case targets; not exact recovery".into());}
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g1_shared_design_is_a_feasible_compromise_not_nominal_target_recovery() {
        let input=options(&[]).unwrap();
        let r=with_cx(|_,cx|fit(cx,&input,false)).unwrap();
        assert!(r.score.is_feasible());
        assert!(!r.every_case_target_met());
        assert!((r.coordinates[0]-25.).abs()<0.1 && (r.coordinates[1]-20.).abs()<0.1);
        assert!((r.evaluation.objective-0.25).abs()<1e-4);
        assert!(r.evaluation.cases.iter().all(|c|c.constraints.is_feasible(TOLERANCE).unwrap()));
    }
    #[test]
    fn g1_worst_case_selection_and_each_hot_spot_limit_remain_physical() {
        let mut input=options(&[]).unwrap();
        input.policy=ScanScenarioObjective::WorstCase;
        let r=with_cx(|_,cx|fit(cx,&input,false)).unwrap();
        assert!(r.score.is_feasible() && !r.every_case_target_met());
        assert!((r.evaluation.objective-0.25).abs()<0.01);
        input.peak_fraction=0.32;
        let limited=with_cx(|_,cx|fit(cx,&input,false)).unwrap();
        assert!(limited.score.is_feasible());
        assert!(limited.evaluation.cases.iter().all(|c|c.constraints.peak_values[0]<=0.32+TOLERANCE));
        assert!(limited.evaluation.objective>0.34);
    }
    #[test]
    #[cfg(not(target_arch="wasm32"))]
    fn g5_fixed_work_case_search_replays_across_workers_and_budget_extensions() {
        let mut input=options(&[]).unwrap(); input.attempts=80;
        let serial=with_cx(|_,cx|fit(cx,&input,false)).unwrap();
        for workers in [None,Some(1),Some(4)] {
            input.workers=workers;
            let other=with_cx(|_,cx|fit(cx,&input,true)).unwrap();
            assert_eq!(other.coordinates,serial.coordinates);
            assert_eq!(other.score,serial.score);
            assert_eq!(other.evaluation.residuals,serial.evaluation.residuals);
            assert_eq!(other.physical_evaluations,vec![81,81]);
        }
    }
    #[test]
    fn g0_case_lists_weights_and_cancellation_are_explicit() {
        for args in [vec!["--cases",""],vec!["--cases","NaN"],vec!["--weights","1"],
            vec!["--weights","1,0"],vec!["--workers","0"],vec!["--worst-case","--worst-case"]] {
            assert!(options(&args.into_iter().map(str::to_string).collect::<Vec<_>>()).is_err());
        }
        let input=options(&["--cases".into(),"80".into(),"--weights".into(),"2".into()]).unwrap();
        assert_eq!(input.enthalpies,vec![80.]);
        with_cx(|gate,cx|{gate.request();assert!(fit(cx,&input,false).is_err());Ok(())}).unwrap();
    }
}
