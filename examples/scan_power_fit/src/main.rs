//! Recover two absorbed-power levels using the native thermal history adjoint.
//! All geometry and materials are synthetic, not experimental calibration.
mod geometry;
mod limits;
mod state_limits;
mod peak_limits;
mod zero_power;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_ascent::{LbfgsReport, StopReason, StopRule};
use fs_ascent::scan_design::{ScanPowerControl, ScanPowerObjective, ScanPowerStudy};
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun, EnthalpyHistoryMaterials, EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::{ScanSensitivityError,
    objective::{ScanObservable, ScanObservation, ScanObjectiveSampling, ScanObjectiveBudget}};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::checkpoint::{
    ScanCheckpointBudget, ScanCheckpointWork};
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

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const USAGE: &str = "scan-power-fit [first_target_J_per_kg second_target_J_per_kg [initial_first_W initial_second_W]]\nDefaults: 152 157 1 1. Targets must increase strictly inside (150,200) J/kg.\nA synthetic unit-volume, unit-density body starts at 150 J/kg and 310 K.\nOutput is a three-row held power schedule with # diagnostic comments.\nJoint spatial recovery: --geometry [target_start_m target_end_m target_width_x_m target_first_W target_second_W]";

#[derive(Debug, Clone, Copy)]
struct FitInput { targets: [f64; 2], initial: [f64; 2] }
impl Default for FitInput {
    fn default() -> Self { Self { targets: [152., 157.], initial: [1., 1.] } }
}
impl FitInput {
    fn parse(args: &[String]) -> Result<Self> {
        if ![0, 2, 4].contains(&args.len()) { return Err(USAGE.into()); }
        let mut input = Self::default();
        for (i, arg) in args.iter().enumerate() {
            let value: f64 = arg.parse()?;
            if i < 2 { input.targets[i] = value; } else { input.initial[i-2] = value; }
        }
        input.validate()?;
        Ok(input)
    }
    fn validate(self) -> Result<()> {
        if self.targets.iter().chain(&self.initial).any(|v| !v.is_finite())
            || self.targets[0] <= 150. || self.targets[1] <= self.targets[0] || self.targets[1] >= 200.
            || self.initial.iter().any(|&v| v <= 0.) || self.initial[0]+self.initial[1] >= 50.
        { return Err("require 150 < first target < second target < 200 J/kg and positive initial powers summing to less than 50 W".into()); }
        Ok(())
    }
    fn validate_observable(self, observable: ScanObservable) -> Result<()> {
        match observable {
            ScanObservable::SpecificEnthalpy => self.validate(),
            ScanObservable::LiquidMassFraction => {
                if self.targets.iter().chain(&self.initial).any(|v| !v.is_finite())
                    || self.targets[0] <= 0.5 || self.targets[1] <= self.targets[0] || self.targets[1] >= 1.
                    || self.initial.iter().any(|&v| v <= 0.) || self.initial[0]+self.initial[1] >= 50.
                { return Err("require 0.5 < first liquid mass fraction < second < 1, and positive initial watts summing to less than 50".into()); }
                Ok(())
            }
            ScanObservable::Temperature => Err("this isothermal fixture fits enthalpy or liquid fraction, not temperature".into()),
        }
    }
}

#[derive(Debug)]
struct FitResult {
    report: LbfgsReport,
    powers: Vec<f64>,
    predictions: Vec<f64>,
    final_temperatures: Vec<f64>,
    final_fractions: Vec<f64>,
    evaluations: usize,
    loss: f64,
    checkpoint_work: Option<ScanCheckpointWork>,
}

fn fit(cx: &Cx<'_>, input: FitInput, replace_sources: bool) -> Result<FitResult> {
    fit_with_grid(cx, input, replace_sources, 2, None, ScanObservable::SpecificEnthalpy)
}

fn fit_with_grid(cx: &Cx<'_>, input: FitInput, replace_sources: bool,
    endpoints: usize, segment: Option<usize>, observable: ScanObservable) -> Result<FitResult> {
    input.validate_observable(observable)?;
    if endpoints < 2 || endpoints > 4096 || endpoints % 2 != 0
        || segment.is_some_and(|n| n == 0 || n > endpoints)
    { return Err("require an even endpoint count in 2..=4096 and a segment length in 1..=endpoints".into()); }
    let retained = segment.unwrap_or(endpoints);
    let (complex, positions) = fs_conduction::fixtures::box_grid([1,1,1], [1.;3]);
    let mesh = ConductionMesh::new(complex, positions)?;
    let material = ConductivityModel::isotropic_declared(2.)?;
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([96;32]),
        [(-200.,280.,0.), (100.,310.,0.), (200.,310.,1.), (900.,380.,1.)].into_iter()
            .map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h, temperature_k:t,
                liquid_mass_fraction:f, bulk_density_kg_m3:1. }).collect())?;
    let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface", ThermalBc::adiabatic())?, &Environment::earth_lab(), vec![], 0.)?;
    let source = MovingHeatSource::new(std::array::from_fn(|_| TimeSignal::Constant(
        QtyAny::new(0.5, Dims([1,0,0,0,0,0])))), [2.;3], TimeSignal::Table {
            times:vec![0.,1.,2.], values:vec![input.initial[0],input.initial[1],0.],
            dims:COMPONENT_POWER_DIMS, interp:Interp::Hold }, 0.,
        MovingHeatBudget { max_vertices:8, max_elements:6, min_active_vertices:1, power_relative_tolerance:1e-11 })?;
    let factory = |cx: &Cx<'_>, mut sources: Vec<MovingHeatSource>|
        -> std::result::Result<ScanningEnthalpyRun<'_, '_, '_>, ScanSensitivityError>
    {
        if replace_sources { sources[0] = sources[0].with_absorbed_power_values(&[9.,9.,0.], 0., 3)?; }
        let base = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem { boundary:&boundary, material:&material,
            element_materials:None, interfaces:None, source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
            EnthalpyHistoryMaterials { table:&[ReferenceEnthalpyMaterial { curve:&curve, reference_density_kg_m3:1. }],
                of_element:&[0;6], phase_conductivity:None }, &[150.;8],
            EnthalpyHistoryBudget { spatial:EnthalpyBudget { max_vertices:8, max_elements:6 }, max_steps:endpoints })?;
        Ok(ScanningEnthalpyRun::new(cx, base, sources, 1, None)?)
    };
    let objective = ScanPowerObjective { endpoints_s:(1..=endpoints).map(|i| 2.0*i as f64/endpoints as f64).collect(), observations:vec![
        ScanObservation { endpoint:endpoints/2-1, vertex:0, observable,
            target:input.targets[0], scale:if observable==ScanObservable::LiquidMassFraction {0.1} else {10.} },
        ScanObservation { endpoint:endpoints-1, vertex:0, observable,
            target:input.targets[1], scale:if observable==ScanObservable::LiquidMassFraction {0.1} else {10.} }],
        sampling:ScanObjectiveSampling::EventLeft,
        step:EnthalpyStepConfig { newton:NewtonKrylovConfig { absolute_tolerance:1e-11, relative_tolerance:1e-12,
            linear_restart:8, max_linear_cycles:16, forcing_minimum:1e-12, forcing_maximum:1e-3, ..Default::default() },
            max_newton_iterations:64, energy_tolerance_j:1e-8 },
        transpose:LinearConfig { tolerance:1e-10, max_iterations:1000, restart:8 },
        budget:ScanObjectiveBudget { max_steps:endpoints, max_observations:2, max_retained_vertices:8*retained,
            max_retained_source_vertices:8*retained, max_schedule_controls:16, max_retained_radiation_entries:0 },
        source_reference_time_s:0. };
    let controls = vec![ScanPowerControl { source:0, value:0, reference_w:1. },
        ScanPowerControl { source:0, value:1, reference_w:1. }];
    let mut study = ScanPowerStudy::new(cx, vec![source], controls, objective, factory, 1)?;
    if let Some(segment_steps) = segment {
        study.set_checkpointing(cx, ScanCheckpointBudget { segment_steps,
            max_checkpoint_vertices:8*(endpoints.div_ceil(segment_steps)+1), max_forward_solves:2*endpoints })?;
    }
    let mut state = study.begin_optimizer(cx, 8)?;
    let original = state.x.clone();
    let paused = study.advance_optimizer(cx, &mut state, &StopRule::GradNorm(1e-9), 80)?;
    if paused.reason != StopReason::Budget || state.x != original {
        return Err("exhausted study changed the accepted optimizer checkpoint".into());
    }
    // One of the 128 total evaluations is reserved explicitly for final replay.
    study.extend_evaluation_budget(127)?;
    let report = study.advance_optimizer(cx, &mut state, &StopRule::GradNorm(1e-9), 80)?;
    let spent = study.evaluations();
    let exported = study.sources_at(cx, &state.x)?;
    let powers = exported[0].absorbed_power_values(3)?;
    if powers[2] != 0. || study.evaluations() != spent {
        return Err("source export changed the off state or consumed an evaluation".into());
    }
    study.extend_evaluation_budget(128)?;
    let verified = study.evaluate(cx, &state.x)?;
    if (verified.physical.value - state.f).abs() > 1e-10 * (1. + state.f.abs()) {
        return Err("replayed fitted schedule differs from the accepted objective".into());
    }
    let final_temperatures = verified.physical.final_step.solution().temperature.clone();
    let final_fractions = verified.physical.final_step.solution().liquid_mass_fraction.clone();
    Ok(FitResult { report, powers, predictions:verified.physical.predictions,
        final_temperatures, final_fractions, evaluations:study.evaluations(), loss:verified.physical.value,
        checkpoint_work:verified.checkpoint_work })
}

fn with_context<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed:96, kernel_id:990, tile:0, iteration:0 }, Budget::INFINITE, ExecMode::Deterministic)))
}
fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().is_some_and(|arg|arg == "--geometry") { return geometry::run_cli(&args[1..]); }
    if args.first().is_some_and(|arg|arg == "--limits") { return limits::run_cli(&args[1..]); }
    if args.first().is_some_and(|arg|arg == "--state-limits") { return state_limits::run_cli(&args[1..]); }
    if args.first().is_some_and(|arg|arg == "--peak-limits") { return peak_limits::run_cli(&args[1..]); }
    if args.first().is_some_and(|arg|arg == "--heater-off") { return zero_power::run_cli(&args[1..]); }
    if args.first().is_some_and(|arg|arg == "--melt-fraction") {
        let targets=match &args[1..] {
            [] => [0.6,0.8], [a,b] => [a.parse()?,b.parse()?],
            _ => return Err("--melt-fraction [first_fraction second_fraction], using fractions not percentages".into()),
        };
        let input=FitInput {targets,initial:[5.,5.]};
        let result=with_context(|cx|fit_with_grid(cx,input,false,8,Some(2),ScanObservable::LiquidMassFraction))?;
        return print_result(input,result,ScanObservable::LiquidMassFraction);
    }
    if args.first().is_some_and(|arg|arg == "--checkpointed") {
        let (endpoints, segment) = match &args[1..] {
            [] => (64,8), [n,s] => (n.parse()?,s.parse()?),
            _ => return Err("--checkpointed [even_endpoint_count segment_steps]".into()),
        };
        let input = FitInput::default();
        let result = with_context(|cx| fit_with_grid(cx,input,false,endpoints,Some(segment),ScanObservable::SpecificEnthalpy))?;
        return print_result(input,result,ScanObservable::SpecificEnthalpy);
    }
    if args.len() == 1 && matches!(args[0].as_str(), "--help" | "-h") {
        println!("{USAGE}\nConstrained fitting: --limits [maximum_power_W maximum_absorbed_energy_J]");
        println!("Sampled hot-spot fitting: --peak-limits [maximum_enthalpy_J_per_kg smoothing_J_per_kg]");
        println!("Short-tape fitting: --checkpointed [even_endpoint_count segment_steps] (defaults: 64 8)");
        println!("Melt-state fitting: --melt-fraction [first_fraction second_fraction] (defaults: 0.6 0.8)");
        println!("Heater activation/off optimization: --heater-off [first_target_J_per_kg second_target_J_per_kg [initial_first_W initial_second_W]]");
        return Ok(());
    }
    let input = FitInput::parse(&args)?;
    let result = with_context(|cx| fit(cx, input, false))?;
    print_result(input,result,ScanObservable::SpecificEnthalpy)
}

fn print_result(input: FitInput, result: FitResult, observable: ScanObservable) -> Result<()> {
    println!("# status={:?} objective={:.12e} gradient_norm={:.12e} evaluations={}",
        result.report.reason, result.loss, result.report.grad_norm, result.evaluations);
    if observable==ScanObservable::LiquidMassFraction {
        println!("# target_liquid_mass_fraction={:?} fitted_liquid_mass_fraction={:?}", input.targets, result.predictions);
        println!("# final_nodal_liquid_mass_fractions={:?}", result.final_fractions);
    } else {
        println!("# target_enthalpy_J_per_kg={:?} fitted_enthalpy_J_per_kg={:?}", input.targets, result.predictions);
    }
    println!("# final_nodal_temperatures_K={:?}", result.final_temperatures);
    if let Some(work) = result.checkpoint_work {
        println!("# final_replay: original_endpoints={} recomputed_endpoints={} checkpoint_states={} peak_live_tangents={}",
            work.original_forward_solves,work.recomputed_forward_solves,work.checkpoint_states,work.peak_retained_linearizations);
    }
    println!("time_s,absorbed_power_w\n0,{:.12e}\n1,{:.12e}\n2,0", result.powers[0], result.powers[1]);
    if result.report.reason != StopReason::GradNorm {
        return Err("fit stopped before its requested gradient threshold; the printed schedule is not a convergence claim".into());
    }
    if observable==ScanObservable::LiquidMassFraction && result.loss>1e-10 {
        return Err("stationary fraction fit still misses the requested targets; saturation is not target recovery".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g5_checkpointed_optimizer_retains_the_fit_and_reports_recomputed_work() {
        let input=FitInput::default();
        let full=with_context(|cx|fit_with_grid(cx,input,false,8,None,ScanObservable::SpecificEnthalpy)).unwrap();
        let short=with_context(|cx|fit_with_grid(cx,input,false,8,Some(2),ScanObservable::SpecificEnthalpy)).unwrap();
        assert_eq!(full.powers,short.powers);
        assert_eq!(full.loss,short.loss);
        assert_eq!(full.evaluations,short.evaluations);
        assert_eq!(short.report.reason,StopReason::GradNorm);
        let work=short.checkpoint_work.unwrap();
        assert_eq!(work.original_forward_solves,8);
        assert_eq!(work.recomputed_forward_solves,8);
        assert_eq!(work.checkpoint_states,5);
        assert_eq!(work.peak_retained_linearizations,2);
        assert!(full.checkpoint_work.is_none());
        for (n,s) in [(0,1),(3,1),(2,0),(2,3),(4098,2)] {
            assert!(with_context(|cx|fit_with_grid(cx,input,false,n,Some(s),ScanObservable::SpecificEnthalpy)).is_err());
        }
    }
    #[test]
    fn g1_fraction_targets_recover_heating_levels_with_full_and_short_tapes() {
        for (targets,expected) in [([0.6,0.8],[10.,20.]),([0.55,0.7],[5.,15.])] {
            let input=FitInput {targets,initial:[5.,5.]};
            let full=with_context(|cx|fit_with_grid(cx,input,false,4,None,ScanObservable::LiquidMassFraction)).unwrap();
            let short=with_context(|cx|fit_with_grid(cx,input,false,4,Some(1),ScanObservable::LiquidMassFraction)).unwrap();
            assert_eq!(full.powers,short.powers); assert_eq!(full.loss,short.loss);
            assert_eq!(short.report.reason,StopReason::GradNorm);
            assert!(short.loss<1e-12);
            for i in 0..2 { assert!((short.powers[i]-expected[i]).abs()<1e-5); }
            assert_eq!(short.powers[2],0.);
            assert!(short.final_temperatures.iter().all(|&t|t==310.));
            assert!(short.final_fractions.iter().all(|&f|(f-targets[1]).abs()<1e-7));
        }
    }
    #[test]
    fn g0_fraction_inputs_reject_percentages_and_stationarity_is_not_target_recovery() {
        for targets in [[60.,80.],[0.4,0.8],[0.8,0.6],[0.6,1.],[f64::NAN,0.8]] {
            let input=FitInput {targets,initial:[5.,5.]};
            assert!(input.validate_observable(ScanObservable::LiquidMassFraction).is_err());
        }
        let input=FitInput {targets:[0.6,0.8],initial:[5.,5.]};
        let mut result=with_context(|cx|fit_with_grid(cx,input,false,2,None,ScanObservable::LiquidMassFraction)).unwrap();
        result.loss=1.;
        assert!(print_result(input,result,ScanObservable::LiquidMassFraction).is_err());
    }
    #[test]
    fn g1_original_adjoint_and_lbfgs_recover_both_levels_after_budget_resume() {
        let fit = with_context(|cx| fit(cx, FitInput::default(), false)).unwrap();
        assert!((fit.powers[0]-2.).abs()<1e-5 && (fit.powers[1]-5.).abs()<1e-5);
        assert_eq!(fit.powers[2], 0.);
        assert!(fit.loss<1e-12);
        assert!(fit.evaluations<=128);
        assert_eq!(fit.report.reason, StopReason::GradNorm);
        assert!(fit.final_temperatures.iter().all(|&t| t == 310.));
    }
    #[test]
    fn g1_changed_targets_and_initial_powers_reach_a_different_physical_schedule() {
        let input = FitInput { targets:[153.,159.], initial:[2.,2.] };
        let fit = with_context(|cx| fit(cx, input, false)).unwrap();
        assert!((fit.powers[0]-3.).abs()<1e-5 && (fit.powers[1]-6.).abs()<1e-5);
        assert!(fit.loss<1e-12);
        assert_eq!(fit.report.reason, StopReason::GradNorm);
        assert_eq!(fit.powers[2], 0.);
    }
    #[test]
    fn g0_factory_cannot_ignore_or_replace_the_optimizer_sources() {
        let error = with_context(|cx| fit(cx, FitInput::default(), true)).unwrap_err();
        assert!(error.to_string().contains("replaced the supplied controlled sources"));
    }
    #[test]
    fn g0_invalid_cli_data_refuse_before_any_physical_trial() {
        for args in [vec!["152"], vec!["nan","157"], vec!["152","151"], vec!["150","157"],
            vec!["152","200"], vec!["152","157","0","1"], vec!["152","157","25","25"]] {
            assert!(FitInput::parse(&args.into_iter().map(String::from).collect::<Vec<_>>()).is_err());
        }
        let parsed = FitInput::parse(&["153","159","2","2"].map(String::from)).unwrap();
        assert_eq!(parsed.targets, [153.,159.]);
        assert_eq!(parsed.initial, [2.,2.]);
    }
}
