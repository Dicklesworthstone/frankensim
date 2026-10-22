//! Physical population fitting with energy, temperature and melt-fraction limits.
//! Synthetic uniform unit cube; all trials use the original P1 enthalpy solver.
use scan_power_fit::population::{PopulationScanReplay, PopulationScanStudy};
use fs_ascent::scan_design::{ScanPowerControl, ScanPowerCoordinates, ScanPowerObjective, ScanPowerStudy};
use fs_ascent::scan_design::constraints::{ScanEnergyCap, ScanStateBound, ScanPeakBound,
    ScanResponse, ScanPeakResponse, ScanPowerSlewBound};
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
use fs_exec::{Cx, Budget, CancelGate, ExecMode, StreamKey, TilePool};
use fs_alloc::OperationMemoryLease;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_qty::{Dims, QtyAny};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{TimeSignal, Interp};
use fs_solver::NewtonKrylovConfig;
use scan_power_fit::population::restart::DifferentialEvolutionRestoreLimits;
#[path = "constrained_scan_fit/recovery.rs"]
mod recovery;
use recovery::{Recovery, MAX_BYTES};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const TARGET: f64 = 1e-8;
const TOLERANCE: f64 = 1e-9;
const USAGE: &str = "constrained_scan_fit [--workers 1..128] [--power-slew RISE_W_PER_S FALL_W_PER_S [--ramp-steps 8..128]] [energy_cap_J [peak_fraction [search_attempts]]] [--checkpoint FILE [--pause-after N]] [--resume FILE]";
#[derive(Clone, Copy)]
struct Options {
    energy_j: f64, peak_fraction: f64, attempts: usize, workers: Option<usize>,
    power_slew: Option<[f64; 2]>, ramp_steps: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self { energy_j:60.0, peak_fraction:0.5, attempts:2400, workers:None,
            power_slew:None, ramp_steps:16 }
    }
}
impl Options {
    fn endpoints(self) -> Vec<f64> {
        if self.power_slew.is_some() {
            (1..=3*self.ramp_steps).map(|i| i as f64 / self.ramp_steps as f64).collect()
        } else { vec![1.0, 2.0] }
    }
    fn observation_endpoints(self) -> [usize; 2] {
        if self.power_slew.is_some() { [2*self.ramp_steps-1, 3*self.ramp_steps-1] }
        else { [0, 1] }
    }
}
fn options(args: &[String]) -> Result<Options> {
    let mut values = Vec::new();
    let mut workers = None;
    let mut power_slew = None;
    let mut ramp_steps = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--workers" if workers.is_none() => {
                i += 1;
                workers = Some(args.get(i).ok_or(USAGE)?.parse::<usize>()?);
            }
            "--power-slew" if power_slew.is_none() => {
                let rise = args.get(i+1).ok_or(USAGE)?.parse::<f64>()?;
                let fall = args.get(i+2).ok_or(USAGE)?.parse::<f64>()?;
                if [rise,fall].iter().any(|v| !v.is_finite() || *v < 0.0) {
                    return Err("power slew limits must be finite nonnegative W/s".into());
                }
                power_slew = Some([rise,fall]);
                i += 2;
            }
            "--ramp-steps" if ramp_steps.is_none() => {
                i += 1;
                ramp_steps = Some(args.get(i).ok_or(USAGE)?.parse::<usize>()?);
            }
            value if !value.starts_with("--") => values.push(value),
            _ => return Err(USAGE.into()),
        }
        i += 1;
    }
    if values.len() > 3 || workers.is_some_and(|n| !(1..=128).contains(&n))
        || ramp_steps.is_some_and(|n| !(8..=128).contains(&n))
        || (ramp_steps.is_some() && power_slew.is_none())
    { return Err(USAGE.into()); }
    let energy_j = values.first().map(|s| s.parse()).transpose()?.unwrap_or(60.0);
    let peak_fraction = values.get(1).map(|s| s.parse()).transpose()?.unwrap_or(0.5);
    let attempts = values.get(2).map(|s| s.parse()).transpose()?.unwrap_or(2400);
    let out = Options { energy_j, peak_fraction, attempts, workers, power_slew,
        ramp_steps:ramp_steps.unwrap_or(16) };
    if !out.energy_j.is_finite() || out.energy_j < 0.0 || out.energy_j > 100.0
        || !out.peak_fraction.is_finite() || !(0.0..=1.0).contains(&out.peak_fraction)
        || !(32..=10_000).contains(&out.attempts)
    { return Err("energy in [0,100] J, dimensionless fraction in [0,1], and 32..10000 attempts required".into()); }
    Ok(out)
}

#[cfg(test)]
fn fit(cx: &Cx<'_>, input: Options, partitioned: bool) -> Result<PopulationScanReplay> {
    match fit_with_recovery(cx, input, partitioned, &Recovery::default())? {
        FitResult::Finished(result) => Ok(result),
        FitResult::Paused { .. } => Err("unexpected checkpoint pause".into()),
    }
}
enum FitResult {
    Finished(PopulationScanReplay),
    Paused { search_attempts: usize, physical_evaluations: usize },
}
fn fit_with_recovery(cx: &Cx<'_>, input: Options, partitioned: bool, recovery: &Recovery)
    -> Result<FitResult>
{
    let model = recovery.model(cx, input)?;
    let saved = recovery.load()?;
    let partitioned = partitioned && !recovery.active();
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
    let endpoints = input.endpoints();
    let steps = endpoints.len();
    let horizon_s = *endpoints.last().ok_or("missing thermal endpoints")?;
    let observation_endpoints = input.observation_endpoints();
    let observations: Vec<_> = observation_endpoints.into_iter().zip([0.1,0.3]).map(|(endpoint,target)|
        ScanSpatialObservation { endpoint, readout:readout.clone(),
            observable:ScanObservable::LiquidMassFraction, target, scale:0.1 }).collect();
    let length = Dims([1,0,0,0,0,0]);
    // All eight vertices have equal footprint weights: the analytic mean
    // balance is available as an oracle, not as a replacement solver.
    // Rate-limited mode starts AND ends at fixed zero watts. Its final second
    // is physically simulated; shutdown heat cannot disappear from state limits.
    let (times, values, interp, control_values) = if input.power_slew.is_some() {
        (vec![0.,1.,2.,3.], vec![0.;4], Interp::Linear, [1,2])
    } else { (vec![0.,1.,2.], vec![0.;3], Interp::Hold, [0,1]) };
    let source = MovingHeatSource::new(
        std::array::from_fn(|_| TimeSignal::Constant(QtyAny::new(0.5, length))), [2.;3],
        TimeSignal::Table { times, values, dims:COMPONENT_POWER_DIMS, interp }, 0.,
        MovingHeatBudget { max_vertices:8, max_elements:6, min_active_vertices:1,
            power_relative_tolerance:1e-11 })?;
    let factory = |cx:&Cx<'_>, sources| -> std::result::Result<_, ScanSensitivityError> {
        let run = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {
            boundary:&boundary, material:&material, element_materials:None, interfaces:None,
            source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
            EnthalpyHistoryMaterials { table:&[ReferenceEnthalpyMaterial {
                curve:&curve, reference_density_kg_m3:1. }], of_element:&[0;6], phase_conductivity:None },
            &[80.;8], EnthalpyHistoryBudget {
                spatial:EnthalpyBudget { max_vertices:8, max_elements:6 }, max_steps:steps })?;
        Ok(ScanningEnthalpyRun::new(cx, run, sources, 1, None)?)
    };
    let initial_cap = if partitioned { 7 } else { input.attempts };
    let study = ScanPowerStudy::new_joint_with_power_coordinates(cx, vec![source],
        control_values.into_iter().map(|value| ScanPowerControl {
            source:0, value, reference_w:1. }).collect(), vec![],
        ScanPowerCoordinates::NonnegativeLinear, ScanPowerObjective {
            endpoints_s:endpoints, observations, sampling:ScanObjectiveSampling::EventLeft,
            step:EnthalpyStepConfig { newton:NewtonKrylovConfig { absolute_tolerance:1e-11,
                relative_tolerance:1e-12, linear_restart:8, max_linear_cycles:16,
                forcing_minimum:1e-12, forcing_maximum:1e-3, ..Default::default() },
                max_newton_iterations:64, energy_tolerance_j:1e-8 },
            transpose:LinearConfig { tolerance:1e-10, max_iterations:1000, restart:8 },
            budget:ScanObjectiveBudget { max_steps:steps, max_observations:2,
                max_retained_vertices:0, max_retained_source_vertices:0,
                max_schedule_controls:16, max_retained_radiation_entries:0 },
            source_reference_time_s:0.,
        }, factory, initial_cap + 1)?;
    let slew_bounds = input.power_slew.map(|[rise,fall]| ScanPowerSlewBound {
        source:0, maximum_rise_w_s:rise, maximum_fall_w_s:fall, scale_w_s:1.0,
    }).into_iter().collect();
    let study = ValueConstrainedScanStudy::new_with_power_slew_bounds(cx, study, vec![],
        vec![ScanEnergyCap { start_time_s:0., end_time_s:horizon_s, maximum_j:input.energy_j, scale_j:1. }],
        vec![ScanStateBound { response:ScanResponse { endpoint:observation_endpoints[0], vertex:0,
            observable:ScanObservable::Temperature }, lower:None, upper:Some(311.), scale:1. }],
        vec![ScanPeakBound { response:ScanPeakResponse { first_endpoint:0, last_endpoint:steps-1,
            vertices:(0..8).collect(), observable:ScanObservable::LiquidMassFraction,
            smoothing_scale:0.001 }, maximum:input.peak_fraction, scale:1. }],
        slew_bounds, if input.power_slew.is_some() {3} else {0},
        ScanValueConstraintBudget { max_rows:if input.power_slew.is_some() {11} else {5}, max_energy_weights:2,
            responses:ScanValueResponseBudget { max_responses:2, max_peak_samples:8*steps } })?;
    let config = DifferentialEvolutionConfig {
        population_size:24, differential_weight:0.8, crossover_probability:0.9, seed:17,
        max_evaluations:initial_cap, max_population_entries:48, target:Some(TARGET),
    };
    let mut driver = match saved {
        Some(bytes) => PopulationScanStudy::from_checkpoint_bytes(cx, study, &bytes, model,
            &[[0.,50.];2], &config, TOLERANCE, DifferentialEvolutionRestoreLimits {
                max_bytes:MAX_BYTES, max_dimension:2, max_population_entries:48 })?,
        None => PopulationScanStudy::new(cx, study, &[[0.,50.];2], config, TOLERANCE)?,
    };
    let start_attempts = driver.search_evaluations();
    let pool = input.workers.map(|threads| TilePool::for_host(threads, 17));
    // Executor bookkeeping is leased separately from the original physical
    // solver's own spatial/work budgets; this is not a total RSS guarantee.
    let lease = OperationMemoryLease::bounded(64 << 20);
    if partitioned {
        let paused = match &pool {
            Some(pool) => driver.advance_parallel(cx, pool, &lease, 7, 24)?,
            None => driver.advance(cx, 7)?,
        };
        assert_eq!(paused.search.evaluations, 7);
        driver.extend_evaluation_budgets(input.attempts, input.attempts + 1)?;
    }
    loop {
        let additional = if recovery.checkpoint.is_some() { 24 }
            else if partitioned { 13 } else { input.attempts };
        let additional = additional.min(recovery.pause_after.map_or(usize::MAX,
            |n| n.saturating_sub(driver.search_evaluations() - start_attempts)));
        let attempted: Result<_> = match &pool {
            Some(pool) => driver.advance_parallel(cx, pool, &lease, additional, 24).map_err(Into::into),
            None => driver.advance(cx, additional).map_err(Into::into),
        };
        let report = match attempted {
            Ok(report) => report,
            Err(error) => {
                if recovery.checkpoint.is_some() {
                    // Workers are drained; a fresh context may capture completed
                    // siblings even when the simulation context is cancelled.
                    if let Err(save_error) = with_cx(|_, fresh| {
                        recovery.save(&driver.checkpoint_bytes(fresh, model, MAX_BYTES)?)
                    }) { eprintln!("checkpoint after failure: {save_error}"); }
                }
                return Err(error);
            }
        };
        if recovery.checkpoint.is_some() {
            recovery.save(&driver.checkpoint_bytes(cx, model, MAX_BYTES)?)?;
        }
        if report.physical_budget_exhausted() || report.search.stop != DifferentialEvolutionStop::Paused { break; }
        if recovery.pause_after.is_some_and(|n| driver.search_evaluations() - start_attempts >= n) {
            return Ok(FitResult::Paused { search_attempts:driver.search_evaluations(),
                physical_evaluations:driver.physical_evaluations() });
        }
    }
    let replay = match driver.replay_best(cx) {
        Ok(replay) => replay,
        Err(error) => {
            if recovery.checkpoint.is_some() {
                if let Err(save_error) = with_cx(|_, fresh| {
                    recovery.save(&driver.checkpoint_bytes(fresh, model, MAX_BYTES)?)
                }) { eprintln!("checkpoint after failed replay: {save_error}"); }
            }
            return Err(error.into());
        }
    };
    if recovery.checkpoint.is_some() {
        recovery.save(&driver.checkpoint_bytes(cx, model, MAX_BYTES)?)?;
    }
    // Cancellation may occur after a callback starts but before its physical
    // reservation. Restored accounts must stay separate, not be equated.
    if replay.search_evaluations > input.attempts
        || replay.physical_evaluations > input.attempts + 1
    {
        return Err("restored search or physical work exceeds its lifetime allowance".into());
    }
    let exported = replay.sources[0].absorbed_power_values(4)?;
    if exported.last().copied().map(f64::to_bits) != Some(0.0_f64.to_bits())
        || (input.power_slew.is_some() && exported[0].to_bits() != 0.0_f64.to_bits())
    {
        return Err("source export changed the fixed shutdown".into());
    }
    Ok(FitResult::Finished(replay))
}

fn with_cx<T>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> Result<T>) -> Result<T> {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&gate, &Cx::new(&gate, arena,
        StreamKey { seed:92, kernel_id:51, tile:0, iteration:0 }, Budget::INFINITE, ExecMode::Deterministic)))
}
fn main() -> Result<()> {
    let (args, recovery) = Recovery::parse(std::env::args().skip(1).collect())?;
    let input = options(&args)?;
    let result = match with_cx(|_,cx| fit_with_recovery(cx, input, true, &recovery))? {
        FitResult::Finished(result) => result,
        FitResult::Paused { search_attempts, physical_evaluations } => {
            println!("# checkpoint_saved; paused=true; search_attempts={search_attempts}; physical_evaluations={physical_evaluations}; no recovery claim");
            return Ok(());
        }
    };
    println!("# synthetic P1 specimen; feasibility-first DE; no optimality certificate");
    println!("# evaluation_mode={}; requested_workers={}; maximum_wave=24",
        if input.workers.is_some() { "parallel" } else { "serial" }, input.workers.unwrap_or(1));
    println!("# feasible={}; target_recovered={}; objective={:.17e}; max_scaled_violation={:.17e}",
        result.score.is_feasible(), result.target_recovered(), result.score.objective(),
        result.score.max_scaled_violation());
    println!("# absorbed_energy_j={:.17e}; sampled_peak_fraction={:.17e}; peak_envelope={:.17e}; smoothing_gap={:.17e}",
        result.constraints.absorbed_energy_j[0], result.constraints.sampled_peak_values[0],
        result.constraints.peak_values[0], result.constraints.peak_smoothing_gaps[0]);
    println!("# search_attempts={}; physical_evaluations={}", result.search_evaluations, result.physical_evaluations);
    if let Some([rise,fall]) = input.power_slew {
        println!("# power_interpolation=linear; rise_limit_w_s={rise}; fall_limit_w_s={fall}; steps_per_second={}", input.ramp_steps);
        println!("# signed_segment_rates_w_s={:?}", result.constraints.power_slew_values_w_s);
        println!("time_s,power_w");
        for (time,power) in result.sources[0].absorbed_power_values(4)?.into_iter().enumerate() {
            println!("{time},{power:.17e}");
        }
        println!("# thermal_measurements");
        println!("time_s,mean_liquid_mass_fraction,target");
        for (i,target) in [0.1,0.3].into_iter().enumerate() {
            println!("{},{:.17e},{}", i+2, result.physical.physical.predictions[i], target);
        }
    } else {
        println!("start_s,end_s,power_w,mean_liquid_mass_fraction,target");
        for (i,target) in [0.1,0.3].into_iter().enumerate() {
            println!("{},{},{:.17e},{:.17e},{}", i, i+1, result.physical.powers_w[i],
                result.physical.physical.predictions[i], target);
        }
        println!("2,2,0,,");
    }
    if !result.score.is_feasible() { return Err("search found no feasible design; this is not an infeasibility proof".into()); }
    if !result.target_recovered() { return Err("feasible best-so-far still misses fitting targets; no target recovery or optimum is claimed".into()); }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    struct SavedRun(std::path::PathBuf);
    impl SavedRun {
        fn new() -> Self {
            static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let path = std::env::temp_dir().join(format!("fs-constrained-restart-{}-{}-{}",
                std::process::id(), std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos(),
                NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
        fn path(&self) -> std::path::PathBuf { self.0.join("search.checkpoint") }
    }
    impl Drop for SavedRun {
        fn drop(&mut self) { let _ = std::fs::remove_dir_all(&self.0); }
    }

    #[test]
    fn g5_saved_thermal_search_reconstructs_model_and_both_work_accounts() {
        let input = Options { energy_j:45., attempts:80, ..Options::default() };
        let expected = with_cx(|_,cx| fit(cx, input, false)).unwrap();
        let saved = SavedRun::new();
        let paused = with_cx(|_,cx| fit_with_recovery(cx, input, false, &Recovery {
            checkpoint:Some(saved.path()), pause_after:Some(7), ..Default::default()
        })).unwrap();
        assert!(matches!(paused, FitResult::Paused { search_attempts:7, physical_evaluations:7 }));
        let bytes = std::fs::read(saved.path()).unwrap();
        assert_eq!(u64::from_le_bytes(bytes[40..48].try_into().unwrap()), 7);
        let restored = with_cx(|_,cx| fit_with_recovery(cx, input, false, &Recovery {
            resume:Some(saved.path()), ..Default::default()
        })).unwrap();
        let FitResult::Finished(actual) = restored else { panic!("unexpected restored pause") };
        assert_eq!(actual.coordinates, expected.coordinates);
        assert_eq!(actual.score, expected.score);
        assert_eq!(actual.constraints, expected.constraints);
        assert_eq!(actual.physical.physical.predictions, expected.physical.physical.predictions);
        assert_eq!(actual.search_evaluations, 80);
        assert_eq!(actual.physical_evaluations, 81);
    }

    #[test]
    #[cfg(not(target_arch = "wasm32"))]
    fn g5_checkpoint_resume_may_change_worker_count_without_changing_fixed_work_physics() {
        let input = Options { energy_j:45., attempts:80, workers:Some(4), ..Options::default() };
        let expected = with_cx(|_,cx| fit(cx, input, false)).unwrap();
        let saved = SavedRun::new();
        let paused = with_cx(|_,cx| fit_with_recovery(cx, input, false, &Recovery {
            checkpoint:Some(saved.path()), pause_after:Some(29), ..Default::default()
        })).unwrap();
        assert!(matches!(paused, FitResult::Paused { search_attempts:29, physical_evaluations:29 }));
        let restored = with_cx(|_,cx| fit_with_recovery(cx, Options { workers:Some(1), ..input }, false,
            &Recovery { resume:Some(saved.path()), ..Default::default() })).unwrap();
        let FitResult::Finished(actual) = restored else { panic!("unexpected restored pause") };
        assert_eq!(actual.coordinates, expected.coordinates);
        assert_eq!(actual.score, expected.score);
        assert_eq!(actual.constraints, expected.constraints);
        assert_eq!(actual.search_evaluations, expected.search_evaluations);
        assert_eq!(actual.physical_evaluations, expected.physical_evaluations);
    }

    #[test]
    fn g0_recovery_rejects_changed_physics_changed_work_caps_and_corrupt_files() {
        let input = Options { energy_j:45., attempts:80, ..Options::default() };
        let saved = SavedRun::new();
        with_cx(|_,cx| fit_with_recovery(cx, input, false, &Recovery {
            checkpoint:Some(saved.path()), pause_after:Some(7), ..Default::default()
        })).unwrap();
        let load = Recovery { resume:Some(saved.path()), ..Default::default() };
        for changed in [Options { energy_j:46., ..input }, Options { peak_fraction:0.6, ..input },
            Options { attempts:81, ..input }] {
            assert!(with_cx(|_,cx| fit_with_recovery(cx, changed, false, &load)).is_err());
        }
        let mut bytes = std::fs::read(saved.path()).unwrap();
        bytes[40] ^= 1;
        std::fs::write(saved.path(), bytes).unwrap();
        assert!(with_cx(|_,cx| fit_with_recovery(cx, input, false, &load)).is_err());
        for args in [vec!["--pause-after", "0"], vec!["--pause-after", "7"],
            vec!["--checkpoint"], vec!["--resume"], vec!["--resume", "a", "--resume", "b"]] {
            assert!(Recovery::parse(args.into_iter().map(str::to_string).collect()).is_err());
        }
    }

    #[test]
    fn g4_completed_checkpoint_does_not_refund_the_final_physical_replay() {
        let input = Options { energy_j:45., attempts:32, ..Options::default() };
        let saved = SavedRun::new();
        let done = with_cx(|_,cx| fit_with_recovery(cx, input, false, &Recovery {
            checkpoint:Some(saved.path()), ..Default::default()
        })).unwrap();
        let FitResult::Finished(result) = done else { panic!("unexpected pause") };
        assert_eq!(result.physical_evaluations, 33);
        let bytes = std::fs::read(saved.path()).unwrap();
        assert_eq!(u64::from_le_bytes(bytes[40..48].try_into().unwrap()), 33);
        // Loading may succeed, but no new physical replay is affordable. The
        // loader must not silently turn the terminal counter back into zero.
        assert!(with_cx(|_,cx| fit_with_recovery(cx, input, false, &Recovery {
            resume:Some(saved.path()), ..Default::default()
        })).is_err());
    }
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
        let input = Options { energy_j:45., ..Options::default() };
        let result = with_cx(|_,cx| fit(cx, input, false)).unwrap();
        assert!(result.score.is_feasible());
        assert!(!result.target_recovered());
        assert!(result.constraints.absorbed_energy_j[0] <= 45. + TOLERANCE);
        assert!(result.score.objective() > 0.12 && result.score.objective() < 0.14);
    }
    #[test]
    fn g5_pause_extension_replays_the_same_population_and_physical_result() {
        let input = Options { energy_j:45., attempts:80, ..Options::default() };
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
        let input = Options { peak_fraction:0., attempts:32, ..Options::default() };
        let result = with_cx(|_,cx| fit(cx, input, false)).unwrap();
        assert!(!result.score.is_feasible());
        assert!(!result.target_recovered());
        assert!(result.constraints.peak_values[0] > 0.);
    }

    #[test]
    #[cfg(not(target_arch = "wasm32"))]
    fn g5_worker_counts_and_pause_partitions_preserve_constrained_physics_at_fixed_work() {
        // The target cannot be reached under this energy cap, so every mode
        // spends the same prefix; between-wave target stopping is not compared.
        let input = Options { energy_j:45., attempts:80, ..Options::default() };
        let serial = with_cx(|_,cx| fit(cx, input, false)).unwrap();
        for workers in [1,4] {
            for partitioned in [false,true] {
                let parallel = with_cx(|_,cx| fit(cx,
                    Options { workers:Some(workers), ..input }, partitioned)).unwrap();
                assert_eq!(parallel.coordinates, serial.coordinates);
                assert_eq!(parallel.score, serial.score);
                assert_eq!(parallel.constraints, serial.constraints);
                assert_eq!(parallel.physical.physical.predictions, serial.physical.physical.predictions);
                assert_eq!(parallel.search_evaluations, 80);
                assert_eq!(parallel.physical_evaluations, 81);
            }
        }
    }

    #[test]
    fn g0_worker_options_are_explicit_bounded_and_do_not_change_serial_defaults() {
        assert_eq!(options(&[]).unwrap().workers, None);
        let args = ["45", "--workers", "4", "0.5", "80"].map(str::to_string);
        let parsed = options(&args).unwrap();
        assert_eq!(parsed.workers, Some(4));
        assert_eq!(parsed.energy_j, 45.);
        assert_eq!(parsed.attempts, 80);
        for args in [vec!["--workers"], vec!["--workers", "0"], vec!["--workers", "129"],
            vec!["--workers", "2", "--workers", "4"], vec!["--worker", "4"]] {
            assert!(options(&args.into_iter().map(str::to_string).collect::<Vec<_>>()).is_err());
        }
    }

    #[test]
    fn g0_ramped_mode_requires_explicit_finite_rates_and_an_event_aligned_grid() {
        let original = options(&[]).unwrap();
        assert!(original.power_slew.is_none());
        assert_eq!(original.endpoints(), vec![1.0,2.0]);
        for n in [8,16,31,128] {
            let args = ["--power-slew".into(), "50".into(), "30".into(),
                "--ramp-steps".into(), n.to_string()];
            let input = options(&args).unwrap();
            assert_eq!(input.power_slew, Some([50.0,30.0]));
            let times = input.endpoints();
            assert_eq!(times.len(), 3*n);
            assert!(times.windows(2).all(|t| t[0] < t[1]));
            assert_eq!([times[n-1],times[2*n-1],times[3*n-1]], [1.0,2.0,3.0]);
            assert_eq!(input.observation_endpoints(), [2*n-1,3*n-1]);
        }
        for args in [vec!["--power-slew"], vec!["--power-slew","50"],
            vec!["--power-slew","NaN","30"], vec!["--power-slew","50","-1"],
            vec!["--power-slew","50","inf"], vec!["--ramp-steps","16"],
            vec!["--power-slew","50","30","--power-slew","50","30"],
            vec!["--power-slew","50","30","--ramp-steps","7"],
            vec!["--power-slew","50","30","--ramp-steps","129"]] {
            assert!(options(&args.into_iter().map(str::to_string).collect::<Vec<_>>()).is_err());
        }
    }

    #[test]
    fn g1_ramped_physical_replay_includes_shutdown_heat_and_every_rate_row() {
        for n in [8,16] {
            let input = Options { attempts:80, power_slew:Some([50.0,30.0]),
                ramp_steps:n, ..Options::default() };
            let result = with_cx(|_,cx| fit(cx, input, true)).unwrap();
            let p = &result.physical.powers_w;
            // Exact uniform-specimen balance for this right-endpoint grid.
            // The final second adds (N-1)/(2N) * P2, not zero joules.
            let a = (n+1) as f64 / (2*n) as f64;
            let h2 = 80.0 + p[0] + a*p[1];
            let h3 = 80.0 + p[0] + p[1];
            for (&actual, h) in result.physical.physical.predictions.iter().zip([h2,h3]) {
                assert!((actual - ((h-100.0)/100.0).clamp(0.0,1.0)).abs() < 1e-7);
            }
            assert_eq!(result.physical.physical.endpoint_count, 3*n);
            for &h in &result.physical.physical.final_step.solution().specific_enthalpy_j_kg {
                assert!((h-h3).abs() < 1e-6);
            }
            assert!((result.constraints.absorbed_energy_j[0] - p[0] - p[1]).abs() < 1e-10);
            assert_eq!(result.constraints.power_slew_values_w_s, vec![p[0],p[1]-p[0],-p[1]]);
            assert_eq!(result.constraints.residuals.len(), 11);
            assert!(result.constraints.is_feasible(TOLERANCE).unwrap());
            let exported = result.sources[0].absorbed_power_values(4).unwrap();
            assert_eq!(exported[0].to_bits(), 0.0_f64.to_bits());
            assert_eq!(exported[3].to_bits(), 0.0_f64.to_bits());
        }
    }

    #[test]
    fn g5_ramped_restart_retains_rates_grid_and_shutdown_through_worker_changes() {
        let input = Options { attempts:80, power_slew:Some([50.0,30.0]),
            ramp_steps:8, ..Options::default() };
        let expected = with_cx(|_,cx| fit(cx, input, false)).unwrap();
        let saved = SavedRun::new();
        let paused = with_cx(|_,cx| fit_with_recovery(cx, input, false, &Recovery {
            checkpoint:Some(saved.path()), pause_after:Some(7), ..Default::default()
        })).unwrap();
        assert!(matches!(paused, FitResult::Paused { search_attempts:7, .. }));
        let load = Recovery { resume:Some(saved.path()), ..Default::default() };
        for changed in [Options { power_slew:Some([49.0,30.0]), ..input },
            Options { power_slew:Some([50.0,31.0]), ..input },
            Options { ramp_steps:16, ..input }, Options { power_slew:None, ..input }] {
            assert!(with_cx(|_,cx| fit_with_recovery(cx, changed, false, &load)).is_err());
        }
        let workers = if cfg!(target_arch = "wasm32") {None} else {Some(2)};
        let resumed = with_cx(|_,cx| fit_with_recovery(cx,
            Options { workers, ..input }, false, &load)).unwrap();
        let FitResult::Finished(actual) = resumed else { panic!("unexpected resumed pause") };
        assert_eq!(actual.coordinates, expected.coordinates);
        assert_eq!(actual.score, expected.score);
        assert_eq!(actual.constraints, expected.constraints);
        assert_eq!(actual.physical.physical.predictions, expected.physical.physical.predictions);
        assert_eq!(actual.physical_evaluations, expected.physical_evaluations);
    }
}
