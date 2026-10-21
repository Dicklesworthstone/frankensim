//! Physical population fitting with energy, temperature and melt-fraction limits.
//! Synthetic uniform unit cube; all trials use the original P1 enthalpy solver.
use scan_power_fit::population::{PopulationScanReplay, PopulationScanStudy};
use fs_ascent::scan_design::{ScanPowerControl, ScanPowerCoordinates, ScanPowerObjective, ScanPowerStudy};
use fs_ascent::scan_design::constraints::{ScanEnergyCap, ScanStateBound, ScanPeakBound,
    ScanResponse, ScanPeakResponse};
use fs_ascent::scan_design::constraints::value::{ValueConstrainedScanStudy, ScanValueConstraintBudget};
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
use fs_exec::{Cx, Budget, CancelGate, ExecMode, StreamKey};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_qty::{Dims, QtyAny};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{TimeSignal, Interp};
use fs_solver::NewtonKrylovConfig;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const TARGET: f64 = 1e-8;
const TOLERANCE: f64 = 1e-9;
#[derive(Clone, Copy)]
struct Options { energy_j: f64, peak_fraction: f64, attempts: usize }
fn options(args: &[String]) -> Result<Options> {
    if args.len() > 3 { return Err("constrained_scan_fit [energy_cap_J [peak_fraction [search_attempts]]]".into()); }
    let energy_j = args.first().map(|s| s.parse()).transpose()?.unwrap_or(60.0);
    let peak_fraction = args.get(1).map(|s| s.parse()).transpose()?.unwrap_or(0.5);
    let attempts = args.get(2).map(|s| s.parse()).transpose()?.unwrap_or(2400);
    let out = Options { energy_j, peak_fraction, attempts };
    if !out.energy_j.is_finite() || out.energy_j < 0.0 || out.energy_j > 100.0
        || !out.peak_fraction.is_finite() || !(0.0..=1.0).contains(&out.peak_fraction)
        || !(32..=10_000).contains(&out.attempts)
    { return Err("energy in [0,100] J, dimensionless fraction in [0,1], and 32..10000 attempts required".into()); }
    Ok(out)
}

fn fit(cx: &Cx<'_>, input: Options, partitioned: bool) -> Result<PopulationScanReplay> {
    let (complex, positions) = fs_conduction::fixtures::box_grid([1,1,1], [1.;3]);
    let mesh = ConductionMesh::new(complex, positions)?;
    let material = ConductivityModel::isotropic_declared(2.)?;
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([92;32]),
        [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
            .map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h, temperature_k:t,
                liquid_mass_fraction:f, bulk_density_kg_m3:1. }).collect())?;
    let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface", ThermalBc::adiabatic())?, &Environment::earth_lab(), vec![], 0.)?;
    let readout = ScanReadout::volume_mean(cx, &mesh, &[0,1,2,3,4,5], 6)?;
    let observations: Vec<_> = [0.1,0.3].into_iter().enumerate().map(|(endpoint,target)|
        ScanSpatialObservation { endpoint, readout:readout.clone(),
            observable:ScanObservable::LiquidMassFraction, target, scale:0.1 }).collect();
    let length = Dims([1,0,0,0,0,0]);
    // All eight vertices have equal footprint weights: the analytic mean
    // balance is available as an oracle, not as a replacement solver.
    let source = MovingHeatSource::new(
        std::array::from_fn(|_| TimeSignal::Constant(QtyAny::new(0.5, length))), [2.;3],
        TimeSignal::Table { times:vec![0.,1.,2.], values:vec![0.;3],
            dims:COMPONENT_POWER_DIMS, interp:Interp::Hold }, 0.,
        MovingHeatBudget { max_vertices:8, max_elements:6, min_active_vertices:1,
            power_relative_tolerance:1e-11 })?;
    let factory = |cx:&Cx<'_>, sources| -> std::result::Result<_, ScanSensitivityError> {
        let run = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {
            boundary:&boundary, material:&material, element_materials:None, interfaces:None,
            source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
            EnthalpyHistoryMaterials { table:&[ReferenceEnthalpyMaterial {
                curve:&curve, reference_density_kg_m3:1. }], of_element:&[0;6], phase_conductivity:None },
            &[80.;8], EnthalpyHistoryBudget {
                spatial:EnthalpyBudget { max_vertices:8, max_elements:6 }, max_steps:2 })?;
        Ok(ScanningEnthalpyRun::new(cx, run, sources, 1, None)?)
    };
    let initial_cap = if partitioned { 7 } else { input.attempts };
    let study = ScanPowerStudy::new_joint_with_power_coordinates(cx, vec![source],
        vec![ScanPowerControl { source:0, value:0, reference_w:1. },
            ScanPowerControl { source:0, value:1, reference_w:1. }], vec![],
        ScanPowerCoordinates::NonnegativeLinear, ScanPowerObjective {
            endpoints_s:vec![1.,2.], observations, sampling:ScanObjectiveSampling::EventLeft,
            step:EnthalpyStepConfig { newton:NewtonKrylovConfig { absolute_tolerance:1e-11,
                relative_tolerance:1e-12, linear_restart:8, max_linear_cycles:16,
                forcing_minimum:1e-12, forcing_maximum:1e-3, ..Default::default() },
                max_newton_iterations:64, energy_tolerance_j:1e-8 },
            transpose:LinearConfig { tolerance:1e-10, max_iterations:1000, restart:8 },
            budget:ScanObjectiveBudget { max_steps:2, max_observations:2,
                max_retained_vertices:0, max_retained_source_vertices:0,
                max_schedule_controls:16, max_retained_radiation_entries:0 },
            source_reference_time_s:0.,
        }, factory, initial_cap + 1)?;
    let study = ValueConstrainedScanStudy::new(cx, study, vec![],
        vec![ScanEnergyCap { start_time_s:0., end_time_s:2., maximum_j:input.energy_j, scale_j:1. }],
        vec![ScanStateBound { response:ScanResponse { endpoint:0, vertex:0,
            observable:ScanObservable::Temperature }, lower:None, upper:Some(311.), scale:1. }],
        vec![ScanPeakBound { response:ScanPeakResponse { first_endpoint:0, last_endpoint:1,
            vertices:(0..8).collect(), observable:ScanObservable::LiquidMassFraction,
            smoothing_scale:0.001 }, maximum:input.peak_fraction, scale:1. }],
        ScanValueConstraintBudget { max_rows:5, max_energy_weights:2,
            responses:ScanValueResponseBudget { max_responses:2, max_peak_samples:16 } })?;
    let mut driver = PopulationScanStudy::new(cx, study, &[[0.,50.];2], DifferentialEvolutionConfig {
        population_size:24, differential_weight:0.8, crossover_probability:0.9, seed:17,
        max_evaluations:initial_cap, max_population_entries:48, target:Some(TARGET),
    }, TOLERANCE)?;
    assert_eq!(driver.physical_evaluations(), 0);
    if partitioned {
        let paused = driver.advance(cx, 7)?;
        assert_eq!(paused.search.evaluations, 7);
        driver.extend_evaluation_budgets(input.attempts, input.attempts + 1)?;
    }
    loop {
        let report = driver.advance(cx, if partitioned { 13 } else { input.attempts })?;
        if report.physical_budget_exhausted() || report.search.stop != DifferentialEvolutionStop::Paused { break; }
    }
    let replay = driver.replay_best(cx)?;
    if replay.physical_evaluations != replay.search_evaluations + 1 {
        return Err("physical work differs from search attempts plus one replay".into());
    }
    if replay.sources[0].absorbed_power_values(3)?[2].to_bits() != 0.0_f64.to_bits() {
        return Err("source export changed the fixed shutdown".into());
    }
    Ok(replay)
}

fn with_cx<T>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> Result<T>) -> Result<T> {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&gate, &Cx::new(&gate, arena,
        StreamKey { seed:92, kernel_id:51, tile:0, iteration:0 }, Budget::INFINITE, ExecMode::Deterministic)))
}
fn main() -> Result<()> {
    let input = options(&std::env::args().skip(1).collect::<Vec<_>>())?;
    let result = with_cx(|_,cx| fit(cx, input, true))?;
    println!("# synthetic P1 specimen; feasibility-first DE; no optimality certificate");
    println!("# feasible={}; target_recovered={}; objective={:.17e}; max_scaled_violation={:.17e}",
        result.score.is_feasible(), result.target_recovered(), result.score.objective(),
        result.score.max_scaled_violation());
    println!("# absorbed_energy_j={:.17e}; sampled_peak_fraction={:.17e}; peak_envelope={:.17e}; smoothing_gap={:.17e}",
        result.constraints.absorbed_energy_j[0], result.constraints.sampled_peak_values[0],
        result.constraints.peak_values[0], result.constraints.peak_smoothing_gaps[0]);
    println!("# search_attempts={}; physical_evaluations={}", result.search_evaluations, result.physical_evaluations);
    println!("start_s,end_s,power_w,mean_liquid_mass_fraction,target");
    for (i,target) in [0.1,0.3].into_iter().enumerate() {
        println!("{},{},{:.17e},{:.17e},{}", i, i+1, result.physical.powers_w[i],
            result.physical.physical.predictions[i], target);
    }
    println!("2,2,0,,");
    if !result.score.is_feasible() { return Err("search found no feasible design; this is not an infeasibility proof".into()); }
    if !result.target_recovered() { return Err("feasible best-so-far still misses fitting targets; no target recovery or optimum is claimed".into()); }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g1_cold_start_reaches_targets_without_violating_physical_limits() {
        let result = with_cx(|_,cx| fit(cx, options(&[])?, true)).unwrap();
        assert!(result.target_recovered());
        assert!(result.constraints.is_feasible(TOLERANCE).unwrap());
        for (&power, expected) in result.physical.powers_w.iter().zip([30.,20.]) {
            assert!((power - expected).abs() < 0.02);
        }
    }
    #[test]
    fn g1_tighter_energy_cap_returns_a_feasible_mismatch_not_a_false_recovery() {
        let input = Options { energy_j:45., peak_fraction:0.5, attempts:2400 };
        let result = with_cx(|_,cx| fit(cx, input, false)).unwrap();
        assert!(result.score.is_feasible());
        assert!(!result.target_recovered());
        assert!(result.constraints.absorbed_energy_j[0] <= 45. + TOLERANCE);
        assert!(result.score.objective() > 0.12 && result.score.objective() < 0.14);
    }
    #[test]
    fn g5_pause_extension_replays_the_same_population_and_physical_result() {
        let input = Options { energy_j:45., peak_fraction:0.5, attempts:80 };
        let a = with_cx(|_,cx| fit(cx, input, false)).unwrap();
        let b = with_cx(|_,cx| fit(cx, input, true)).unwrap();
        assert_eq!(a.coordinates, b.coordinates);
        assert_eq!(a.score, b.score);
        assert_eq!(a.constraints, b.constraints);
        assert_eq!(a.physical.physical.predictions, b.physical.physical.predictions);
        assert_eq!(a.physical_evaluations, b.physical_evaluations);
    }
    #[test]
    fn g0_invalid_limits_and_cancelled_trials_never_report_recovery() {
        for args in [vec!["-1"], vec!["NaN"], vec!["60","101"], vec!["60","0.5","1"]] {
            assert!(options(&args.into_iter().map(str::to_string).collect::<Vec<_>>()).is_err());
        }
        with_cx(|gate,cx| {
            gate.request();
            assert!(fit(cx, options(&[])?, false).is_err());
            Ok(())
        }).unwrap();
        let input = Options { energy_j:60., peak_fraction:0., attempts:32 };
        let result = with_cx(|_,cx| fit(cx, input, false)).unwrap();
        assert!(!result.score.is_feasible());
        assert!(!result.target_recovered());
        assert!(result.constraints.peak_values[0] > 0.);
    }
}
