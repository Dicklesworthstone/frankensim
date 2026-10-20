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
use fs_dfo::differential_evolution::{DifferentialEvolutionConfig, DifferentialEvolutionReport,
    DifferentialEvolutionState, DifferentialEvolutionStop};
use fs_exec::{Cx, Budget, CancelGate, ExecMode, StreamKey};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_qty::{Dims, QtyAny};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{TimeSignal, Interp};
use fs_solver::NewtonKrylovConfig;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const TARGET: f64 = 1e-12;
const USAGE: &str = "simplex_scan_fit [--population] [--initial-enthalpy J/kg] [--seed u64] [first_fraction second_fraction [maximum_W total_evaluations]]";

#[derive(Clone, Copy)]
struct Options {
    targets: [f64;2], maximum_w: f64, evaluations: usize,
    population: bool, initial_h: f64, seed: u64,
}
fn options(args: &[String]) -> Result<Options> {
    let mut values=Vec::new();
    let (mut population,mut initial_h,mut seed)=(false,100.0_f64,17_u64);
    let mut seen=[false;3];
    let mut i=0;
    while i<args.len() {
        match args[i].as_str() {
            "--population" if !seen[0] => { population=true;seen[0]=true; },
            "--initial-enthalpy" if !seen[1] => {
                i+=1;initial_h=args.get(i).ok_or(USAGE)?.parse()?;seen[1]=true;
            },
            "--seed" if !seen[2] => {
                i+=1;seed=args.get(i).ok_or(USAGE)?.parse()?;seen[2]=true;
            },
            value if !value.starts_with("--") => values.push(value),
            _ => return Err(USAGE.into()),
        }
        i+=1;
    }
    if !matches!(values.len(),0|2|4) || (seen[2] && !population) { return Err(USAGE.into()); }
    let targets: [f64;2]=if values.is_empty() {[0.1,0.3]} else {[values[0].parse()?,values[1].parse()?]};
    let (maximum_w,evaluations): (f64,usize)=if values.len()==4 {
        (values[2].parse()?,values[3].parse()?)
    } else { (50.0,if population {3000} else {400}) };
    if targets.iter().any(|v| !v.is_finite() || !(0.0..=1.0).contains(v))
        || !maximum_w.is_finite() || maximum_w <= 0.0 || maximum_w > 100.0
        || !(10..=10_000).contains(&evaluations)
        || !initial_h.is_finite() || !(-100.0..=200.0).contains(&initial_h)
    { return Err("fractions in [0,1], watts in (0,100], evaluations in [10,10000], initial enthalpy in [-100,200] J/kg required".into()); }
    Ok(Options { targets,maximum_w,evaluations,population,initial_h,seed })
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

#[derive(Debug, Clone, PartialEq)]
enum SearchReport { Simplex(NelderMeadReport), Population(DifferentialEvolutionReport) }
impl SearchReport {
    fn best(&self) -> Result<(&[f64],f64)> {
        match self {
            Self::Simplex(r) => r.best.as_ref().map(|p|(p.coordinates.as_slice(),p.value)),
            Self::Population(r) => r.best.as_ref().map(|p|(p.coordinates.as_slice(),p.value)),
        }.ok_or_else(||"search returned no evaluated design".into())
    }
    fn evaluations(&self) -> usize {
        match self { Self::Simplex(r)=>r.evaluations,Self::Population(r)=>r.evaluations }
    }
    fn method(&self) -> &'static str {
        match self { Self::Simplex(_)=>"Nelder-Mead",Self::Population(_)=>"DE/rand/1/bin" }
    }
    fn stop(&self) -> &'static str {
        match self {
            Self::Simplex(r) => match r.stop {
                NelderMeadStop::TargetReached=>"TargetReached",NelderMeadStop::SimplexSmall=>"SimplexSmall",
                NelderMeadStop::EvaluationLimit=>"EvaluationLimit",NelderMeadStop::Paused=>"Paused",
            },
            Self::Population(r) => match r.stop {
                DifferentialEvolutionStop::TargetReached=>"TargetReached",
                DifferentialEvolutionStop::EvaluationLimit=>"EvaluationLimit",
                DifferentialEvolutionStop::Paused=>"Paused",
            },
        }
    }
}
struct Outcome {
    search: SearchReport,
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
            // The default is the solidus; an explicit colder start introduces
            // a flat liquid-fraction response before sensible heating reaches it.
            &[input.initial_h;8],EnthalpyHistoryBudget {
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
    let search=if input.population {
        let mut solver=DifferentialEvolutionState::new(cx,&[[0.,input.maximum_w];2],
            Some(study.initial_coordinates()),DifferentialEvolutionConfig {
                population_size:24,differential_weight:0.8,crossover_probability:0.9,seed:input.seed,
                max_evaluations:initial_cap,max_population_entries:48,target:Some(TARGET),
            })?;
        let mut report=solver.advance(cx,usize::MAX,&mut |point|
            study.evaluate_value(cx,point).map(|value|value.physical.value))?;
        if pause && report.stop==DifferentialEvolutionStop::EvaluationLimit {
            study.extend_evaluation_budget(input.evaluations)?;
            solver.extend_evaluation_budget(input.evaluations-1)?;
            report=solver.advance(cx,usize::MAX,&mut |point|
                study.evaluate_value(cx,point).map(|value|value.physical.value))?;
        }
        SearchReport::Population(report)
    } else {
        let mut solver=NelderMeadState::new(study.initial_coordinates(),
            &[input.maximum_w.min(1.);2],NelderMeadConfig {
                max_evaluations:initial_cap,max_simplex_entries:6,target:Some(TARGET),
                x_tolerance:1e-8,f_tolerance:1e-14,bounds:Some(vec![[0.,input.maximum_w];2]),
            })?;
        let mut report=solver.advance(cx,usize::MAX,&mut |point|
            study.evaluate_value(cx,point).map(|value|value.physical.value))?;
        if pause && report.stop==NelderMeadStop::EvaluationLimit {
            study.extend_evaluation_budget(input.evaluations)?;
            solver.extend_evaluation_budget(input.evaluations-1)?;
            report=solver.advance(cx,usize::MAX,&mut |point|
                study.evaluate_value(cx,point).map(|value|value.physical.value))?;
        }
        SearchReport::Simplex(report)
    };
    // Both search methods reserve the same final physical replay and keep
    // original callbacks charged when their lifetime allowance is extended.
    let (best_coordinates,best_value)=search.best()?;
    let before_export=study.evaluations();
    let sources=study.sources_at(cx,best_coordinates)?;
    if study.evaluations()!=before_export || sources[0].absorbed_power_values(3)?[2]!=0. {
        return Err("source export changed work accounting or fixed shutdown".into());
    }
    let replay=study.evaluate_value(cx,best_coordinates)?;
    if replay.physical.value.to_bits()!=best_value.to_bits() {
        return Err("best-candidate replay changed the physical objective".into());
    }
    if study.evaluations()!=search.evaluations()+1 || study.evaluations()>input.evaluations {
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
    println!("# synthetic specimen; search={}; no global-optimality guarantee; initial_h_j_kg={}; seed={}",
        outcome.search.method(),input.initial_h,input.seed);
    println!("# stop={}; objective={:.17e}; callback_attempts={}; total_physical_evaluations={}",
        outcome.search.stop(),outcome.replay.physical.value,outcome.search.evaluations(),outcome.spent);
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
        assert!(matches!(&result.search,SearchReport::Simplex(r) if r.stop==NelderMeadStop::TargetReached));
        for (actual,expected) in result.replay.powers_w.iter().zip([10.,20.]) {
            assert!((*actual-expected).abs()<1e-3,"{actual} != {expected}");
        }
        assert!(result.replay.physical.final_step.solution().temperature.iter().all(|t|(*t-310.).abs()<1e-8));
        assert_eq!(result.spent,result.search.evaluations()+1);
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

    #[test]
    fn g1_population_heats_through_a_flat_solid_regime_and_resumes_identically() {
        let input=options(&["--population","--initial-enthalpy","80"]
            .map(str::to_owned)).unwrap();
        with_cx(|_,cx| {
            let paused=fit(cx,input,true)?;
            let direct=fit(cx,input,false)?;
            assert!(paused.target_met());
            assert!(matches!(&paused.search,SearchReport::Population(r)
                if r.stop==DifferentialEvolutionStop::TargetReached));
            assert_eq!(paused.search,direct.search);
            assert_eq!(paused.replay.physical.predictions,direct.replay.physical.predictions);
            for (actual,expected) in paused.replay.powers_w.iter().zip([30.,20.]) {
                assert!((*actual-expected).abs()<2e-3,"{actual} != {expected}");
            }
            assert_eq!(paused.spent,paused.search.evaluations()+1);
            assert!(paused.spent<=input.evaluations);
            let local=fit(cx,Options {population:false,evaluations:400,..input},true)?;
            assert!(!local.target_met());
            assert!(local.replay.physical.value>1.);
            Ok(())
        }).unwrap();
    }
    #[test]
    fn g0_population_reports_unreachable_targets_and_rejects_ambiguous_options() {
        for args in [vec!["--population","--population"],vec!["--seed","1"],
            vec!["--population","--seed","-1"],vec!["--initial-enthalpy"],
            vec!["--initial-enthalpy","NaN"],vec!["--initial-enthalpy","900"]] {
            assert!(options(&args.into_iter().map(str::to_owned).collect::<Vec<_>>()).is_err());
        }
        let input=options(&["--population","--initial-enthalpy","80","0.1","0.3","1","65"]
            .map(str::to_owned)).unwrap();
        let result=with_cx(|_,cx|fit(cx,input,true)).unwrap();
        assert!(!result.target_met());
        assert!(matches!(&result.search,SearchReport::Population(r)
            if r.stop==DifferentialEvolutionStop::EvaluationLimit));
        assert_eq!(result.spent,65);
        assert!(result.replay.powers_w.iter().all(|&p|(0.0..=1.0).contains(&p)));
        with_cx(|gate,cx| {gate.request();assert!(fit(cx,input,true).is_err());Ok(())}).unwrap();
    }
}
