//! Derivative-free fitting from a physical phase-chart corner and zero heating.
//! Synthetic fixed unit cube, density 1 kg/m^3 and latent heat 100 J/kg.
//! Uses the original P1 enthalpy solver and spatial observation owner. No
//! lumped surrogate, finite-difference gradient or failed-solve penalty is used.
use fs_ascent::scan_design::{ScanPowerControl, ScanPowerCoordinates, ScanPowerObjective, ScanPowerStudy};
use fs_ascent::scan_design::value::ScanPowerValueEvaluation;
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
use fs_dfo::neldermead::state::{NelderMeadConfig, NelderMeadReport, NelderMeadState, NelderMeadStop};
use fs_exec::{Cx, Budget, CancelGate, ExecMode, StreamKey};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_qty::{Dims, QtyAny};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{TimeSignal, Interp};
use fs_solver::NewtonKrylovConfig;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const TARGET: f64 = 1e-12;
const USAGE: &str = "simplex_scan_fit [first_fraction second_fraction [maximum_W total_evaluations]]";

#[derive(Clone, Copy)]
struct Options { targets:[f64;2], maximum_w:f64, evaluations:usize }
fn options(args: &[String]) -> Result<Options> {
    if !matches!(args.len(),0|2|4) { return Err(USAGE.into()); }
    let targets: [f64;2]=if args.is_empty() {[0.1,0.3]} else {[args[0].parse()?,args[1].parse()?]};
    let (maximum_w,evaluations): (f64,usize)=if args.len()==4 {(args[2].parse()?,args[3].parse()?)} else {(50.0,400)};
    if targets.iter().any(|v| !v.is_finite() || !(0.0..=1.0).contains(v))
        || !maximum_w.is_finite() || maximum_w <= 0.0 || maximum_w > 100.0
        || !(10..=10_000).contains(&evaluations)
    { return Err("fractions must be in [0,1], maximum watts in (0,100], total evaluations in [10,10000]".into()); }
    Ok(Options { targets, maximum_w, evaluations })
}

fn source() -> Result<MovingHeatSource> {
    let length=Dims([1,0,0,0,0,0]);
    Ok(MovingHeatSource::new([
        TimeSignal::Table { times:vec![0.,1.,2.], values:vec![0.25,0.75,0.75], dims:length, interp:Interp::Hold },
        TimeSignal::Constant(QtyAny::new(0.5,length)),
        TimeSignal::Constant(QtyAny::new(0.5,length))], [1.4,1.2,1.1],
        TimeSignal::Table { times:vec![0.,1.,2.], values:vec![0.;3], dims:COMPONENT_POWER_DIMS,
            interp:Interp::Hold }, 0., MovingHeatBudget { max_vertices:8, max_elements:6,
            min_active_vertices:1, power_relative_tolerance:1e-11 })?)
}

struct Outcome {
    search: NelderMeadReport,
    replay: ScanPowerValueEvaluation,
    spent: usize,
    sources: Vec<MovingHeatSource>,
}
impl Outcome {
    fn target_met(&self) -> bool { self.replay.physical.value <= TARGET }
}

fn fit(cx: &Cx<'_>, input: Options, pause: bool) -> Result<Outcome> {
    let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
    let mesh=ConductionMesh::new(complex,positions)?;
    let material=ConductivityModel::isotropic_declared(2.)?;
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([92;32]),
        [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
            .map(|(h,t,f)|EnthalpyPhaseKnot { specific_enthalpy_j_kg:h, temperature_k:t,
                liquid_mass_fraction:f, bulk_density_kg_m3:1. }).collect())?;
    let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface",ThermalBc::adiabatic())?,&Environment::earth_lab(),vec![],0.)?;
    let readout=ScanReadout::volume_mean(cx,&mesh,&[0,1,2,3,4,5],6)?;
    let observations: Vec<ScanSpatialObservation<'_>>=input.targets.into_iter().enumerate().map(|(endpoint,target)|
        ScanSpatialObservation { endpoint, readout:readout.clone(),
            observable:ScanObservable::LiquidMassFraction, target, scale:0.1 }).collect();
    let factory=|cx:&Cx<'_>,sources| -> std::result::Result<_,ScanSensitivityError> {
        let run=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {
            boundary:&boundary, material:&material, element_materials:None, interfaces:None,
            source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
            EnthalpyHistoryMaterials { table:&[ReferenceEnthalpyMaterial {
                curve:&curve, reference_density_kg_m3:1. }], of_element:&[0;6], phase_conductivity:None },
            // Exactly at the solidus: a valid state, not a two-sided adjoint point.
            &[100.;8],EnthalpyHistoryBudget {
                spatial:EnthalpyBudget { max_vertices:8, max_elements:6 },max_steps:2 })?;
        Ok(ScanningEnthalpyRun::new(cx,run,sources,1,None)?)
    };
    let initial_cap=if pause {8} else {input.evaluations-1};
    let mut study=ScanPowerStudy::new_joint_with_power_coordinates(cx,vec![source()?],
        vec![ScanPowerControl {source:0,value:0,reference_w:1.},
             ScanPowerControl {source:0,value:1,reference_w:1.}],vec![],
        ScanPowerCoordinates::NonnegativeLinear,ScanPowerObjective {
            endpoints_s:vec![1.,2.],observations,sampling:ScanObjectiveSampling::EventLeft,
            step:EnthalpyStepConfig { newton:NewtonKrylovConfig { absolute_tolerance:1e-11,
                relative_tolerance:1e-12,linear_restart:8,max_linear_cycles:16,
                forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default() },
                max_newton_iterations:64,energy_tolerance_j:1e-8 },
            // Deliberately zero retained-tangent allowance: no hidden derivative path.
            transpose:LinearConfig {tolerance:1e-10,max_iterations:1000,restart:8},
            budget:ScanObjectiveBudget {max_steps:2,max_observations:2,
                max_retained_vertices:0,max_retained_source_vertices:0,
                max_schedule_controls:16,max_retained_radiation_entries:0},source_reference_time_s:0.,
        },factory,initial_cap+1)?;
    let mut solver=NelderMeadState::new(study.initial_coordinates(),
        &[input.maximum_w.min(1.);2],NelderMeadConfig {
            max_evaluations:initial_cap,max_simplex_entries:6,target:Some(TARGET),
            x_tolerance:1e-8,f_tolerance:1e-14,bounds:Some(vec![[0.,input.maximum_w];2]),
        })?;
    let mut search=solver.advance(cx,usize::MAX,&mut |point|
        study.evaluate_value(cx,point).map(|value|value.physical.value))?;
    if pause && search.stop==NelderMeadStop::EvaluationLimit {
        // Keep both cumulative counters. One physical evaluation stays reserved
        // for independent replay of the exported best candidate.
        study.extend_evaluation_budget(input.evaluations)?;
        solver.extend_evaluation_budget(input.evaluations-1)?;
        search=solver.advance(cx,usize::MAX,&mut |point|
            study.evaluate_value(cx,point).map(|value|value.physical.value))?;
    }
    let best=search.best.as_ref().ok_or("search returned no evaluated design")?;
    let before_export=study.evaluations();
    let sources=study.sources_at(cx,&best.coordinates)?;
    if study.evaluations()!=before_export || sources[0].absorbed_power_values(3)?[2]!=0. {
        return Err("source export changed work accounting or fixed shutdown".into());
    }
    let replay=study.evaluate_value(cx,&best.coordinates)?;
    if replay.physical.value.to_bits()!=best.value.to_bits() {
        return Err("best-candidate replay changed the physical objective".into());
    }
    if study.evaluations()!=search.evaluations+1 || study.evaluations()>input.evaluations {
        return Err("callback and physical-evaluation budgets disagree".into());
    }
    Ok(Outcome { search,replay,spent:study.evaluations(),sources })
}

fn with_cx<T>(f: impl FnOnce(&CancelGate,&Cx<'_>) -> Result<T>) -> Result<T> {
    let gate=CancelGate::new();
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena|f(&gate,&Cx::new(&gate,arena,
        StreamKey {seed:92,kernel_id:51,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn main() -> Result<()> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    let input=options(&args)?;
    let outcome=with_cx(|_,cx|fit(cx,input,true))?;
    println!("# synthetic specimen; derivative-free local search, not a global optimum");
    println!("# stop={:?}; objective={:.17e}; callback_attempts={}; total_physical_evaluations={}",
        outcome.search.stop,outcome.replay.physical.value,outcome.search.evaluations,outcome.spent);
    println!("start_s,end_s,absorbed_power_w,mean_fraction_at_end,target_fraction_at_end");
    for i in 0..2 {
        println!("{},{},{:.17e},{:.17e},{:.17e}",i,i+1,outcome.replay.powers_w[i],
            outcome.replay.physical.predictions[i],input.targets[i]);
    }
    println!("2,2,{:.17e},,",outcome.sources[0].absorbed_power_values(3)?[2]);
    if !outcome.target_met() {
        return Err("search stopped without meeting the declared targets; reported best is not target recovery".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g1_simplex_fits_melting_from_zero_at_the_original_nondifferentiable_solidus() {
        let input=options(&[]).unwrap();
        let result=with_cx(|_,cx|fit(cx,input,true)).unwrap();
        assert!(result.target_met());
        assert_eq!(result.search.stop,NelderMeadStop::TargetReached);
        for (actual,expected) in result.replay.powers_w.iter().zip([10.,20.]) {
            assert!((*actual-expected).abs()<1e-3,"{actual} != {expected}");
        }
        assert!(result.replay.physical.final_step.solution().temperature.iter().all(|t|(*t-310.).abs()<1e-8));
        assert_eq!(result.spent,result.search.evaluations+1);
        assert!(result.spent<=400);
    }
    #[test]
    fn g3_resumed_physical_search_matches_uninterrupted_and_never_hides_missed_targets() {
        let input=options(&[]).unwrap();
        with_cx(|_,cx| {
            let paused=fit(cx,input,true)?;
            let direct=fit(cx,input,false)?;
            assert_eq!(paused.search,direct.search);
            assert_eq!(paused.replay.physical.predictions,direct.replay.physical.predictions);
            assert_eq!(paused.sources,direct.sources);
            let limited=fit(cx,Options {maximum_w:1.,..input},true)?;
            assert!(!limited.target_met());
            assert!(limited.replay.powers_w.iter().all(|&p|(0.0..=1.0).contains(&p)));
            Ok(())
        }).unwrap();
    }
    #[test]
    fn g0_bad_command_values_and_cancelled_physics_are_refused() {
        for args in [vec!["0.1"],vec!["20","30"],vec!["NaN","0.2"],
            vec!["0.1","0.3","0","400"],vec!["0.1","0.3","50","9"]] {
            assert!(options(&args.into_iter().map(str::to_owned).collect::<Vec<_>>()).is_err());
        }
        with_cx(|gate,cx| { gate.request();assert!(fit(cx,options(&[])?,true).is_err());Ok(()) }).unwrap();
    }
}
