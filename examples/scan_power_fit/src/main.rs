//! Recover two absorbed-power levels using the native thermal history adjoint.
//! All geometry and materials are synthetic, not experimental calibration.
mod geometry;
mod limits;
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
}

#[derive(Debug)]
struct FitResult {
    report: LbfgsReport,
    powers: Vec<f64>,
    predictions: Vec<f64>,
    final_temperatures: Vec<f64>,
    evaluations: usize,
    loss: f64,
}

fn fit(cx: &Cx<'_>, input: FitInput, replace_sources: bool) -> Result<FitResult> {
    input.validate()?;
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
            EnthalpyHistoryBudget { spatial:EnthalpyBudget { max_vertices:8, max_elements:6 }, max_steps:2 })?;
        Ok(ScanningEnthalpyRun::new(cx, base, sources, 1, None)?)
    };
    let objective = ScanPowerObjective { endpoints_s:vec![1.,2.], observations:vec![
        ScanObservation { endpoint:0, vertex:0, observable:ScanObservable::SpecificEnthalpy,
            target:input.targets[0], scale:10. },
        ScanObservation { endpoint:1, vertex:0, observable:ScanObservable::SpecificEnthalpy,
            target:input.targets[1], scale:10. }], sampling:ScanObjectiveSampling::EventLeft,
        step:EnthalpyStepConfig { newton:NewtonKrylovConfig { absolute_tolerance:1e-11, relative_tolerance:1e-12,
            linear_restart:8, max_linear_cycles:16, forcing_minimum:1e-12, forcing_maximum:1e-3, ..Default::default() },
            max_newton_iterations:64, energy_tolerance_j:1e-8 },
        transpose:LinearConfig { tolerance:1e-10, max_iterations:1000, restart:8 },
        budget:ScanObjectiveBudget { max_steps:2, max_observations:2, max_retained_vertices:16,
            max_retained_source_vertices:16, max_schedule_controls:16, max_retained_radiation_entries:0 },
        source_reference_time_s:0. };
    let controls = vec![ScanPowerControl { source:0, value:0, reference_w:1. },
        ScanPowerControl { source:0, value:1, reference_w:1. }];
    let mut study = ScanPowerStudy::new(cx, vec![source], controls, objective, factory, 1)?;
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
    Ok(FitResult { report, powers, predictions:verified.physical.predictions,
        final_temperatures, evaluations:study.evaluations(), loss:verified.physical.value })
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
    if args.len() == 1 && matches!(args[0].as_str(), "--help" | "-h") {
        println!("{USAGE}\nConstrained fitting: --limits [maximum_power_W maximum_absorbed_energy_J]");
        return Ok(());
    }
    let input = FitInput::parse(&args)?;
    let result = with_context(|cx| fit(cx, input, false))?;
    println!("# status={:?} objective={:.12e} gradient_norm={:.12e} evaluations={}",
        result.report.reason, result.loss, result.report.grad_norm, result.evaluations);
    println!("# target_enthalpy_J_per_kg={:?} fitted_enthalpy_J_per_kg={:?}", input.targets, result.predictions);
    println!("# final_nodal_temperatures_K={:?}", result.final_temperatures);
    println!("time_s,absorbed_power_w\n0,{:.12e}\n1,{:.12e}\n2,0", result.powers[0], result.powers[1]);
    if result.report.reason != StopReason::GradNorm {
        return Err("fit stopped before its requested gradient threshold; the printed schedule is not a convergence claim".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
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
