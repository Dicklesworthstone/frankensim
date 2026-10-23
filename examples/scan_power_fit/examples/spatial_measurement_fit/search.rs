//! Select the existing constrained engine without changing the physical model.
//! Population search deliberately requests values only: a constitutive corner
//! is not a failed physical solve merely because its classical derivative fails.
use super::*;
use scan_power_fit::population::PopulationScanStudy;
use scan_power_fit::requirements::BoundScanRequirements;
use fs_ascent::scan_design::constraints::value::{ScanConstraintValues,
    ScanValueConstraintBudget, ValueConstrainedScanStudy};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
    spatial::ScanObservationSpec, value::responses::ScanValueResponseBudget};
use fs_dfo::differential_evolution::{DifferentialEvolutionConfig, DifferentialEvolutionStop};
use fs_alloc::OperationMemoryLease;
use fs_exec::TilePool;

pub(super) const FEASIBILITY_TOLERANCE: f64 = 1e-8;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct PopulationOptions {
    /// Physical watts, applied to each selected power, not optimizer coordinates.
    pub maximum_w: f64,
    pub attempts: usize,
    pub seed: u64,
    /// Threshold on the declared sum of residual losses, not a KKT tolerance.
    pub objective_target: f64,
    /// None is serial; Some selects bounded waves on the original executor.
    pub workers: Option<usize>,
}
impl PopulationOptions {
    pub fn validate(self) -> Result<()> {
        if !self.maximum_w.is_finite() || self.maximum_w <= 0.0
            || !(1..=100_000).contains(&self.attempts)
            || !self.objective_target.is_finite() || self.objective_target < 0.0
            || self.workers.is_some_and(|n| !(1..=128).contains(&n)) {
            return Err("population search needs positive finite maximum watts, 1..100000 attempts, a finite nonnegative objective target and 1..128 workers when supplied".into());
        }
        Ok(())
    }

    // Retain each original power scale. Clipping an out-of-box initial design
    // would silently change the supplied heater program and is not permitted.
    pub fn bounds(self, cx: &Cx<'_>, program: &ScanHeatingProgram) -> Result<Vec<[f64; 2]>> {
        self.validate()?;
        let mut bounds = Vec::with_capacity(program.controls().len());
        for control in program.controls() {
            cx.checkpoint().map_err(|_| "population box admission cancelled")?;
            let values = program.sources()[control.source]
                .absorbed_power_values(program.authored_control_count())?;
            let initial = values[control.value];
            let mut upper = self.maximum_w / control.reference_w;
            if !upper.is_finite() || upper <= 0.0 {
                return Err("population watt box is not representable in the original power coordinates".into());
            }
            if upper * control.reference_w > self.maximum_w { upper = upper.next_down(); }
            if !upper.is_finite() || upper <= 0.0 || initial > self.maximum_w {
                return Err("population box must contain every selected initial power and have representable scaled bounds".into());
            }
            if initial / control.reference_w > upper {
                return Err("selected initial watts cannot be represented inside the declared population box".into());
            }
            bounds.push([0.0, upper]);
        }
        Ok(bounds)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Method { Sqp, Population(PopulationOptions) }
impl Method {
    pub fn evaluation_limit(self) -> usize {
        match self { Self::Sqp => 127, Self::Population(p) => p.attempts + 1 }
    }
}

/// Consume only engine options; the existing parser owns all physical inputs.
/// Unknown arguments stay untouched and are rejected by that parser as before.
pub(super) fn select(args: &[String]) -> Result<(Vec<String>, Method)> {
    let mut remaining=Vec::new();
    let mut population=None; let mut seed=None; let mut target=None; let mut workers=None;
    let mut i=0;
    while i<args.len() {
        match args[i].as_str() {
            "--population" => {
                if population.is_some() { return Err("duplicate --population".into()); }
                i+=1; let maximum_w=finite(args.get(i).ok_or("--population needs MAX_W ATTEMPTS")?)?;
                i+=1; let attempts=args.get(i).ok_or("--population needs MAX_W ATTEMPTS")?.parse::<usize>()?;
                population=Some((maximum_w,attempts));
            }
            "--seed" => {
                if seed.is_some() { return Err("duplicate --seed".into()); }
                i+=1; seed=Some(args.get(i).ok_or("--seed needs U64")?.parse::<u64>()?);
            }
            "--objective-target" => {
                if target.is_some() { return Err("duplicate --objective-target".into()); }
                i+=1; target=Some(finite(args.get(i).ok_or("--objective-target needs LOSS")?)?);
            }
            "--workers" => {
                if workers.is_some() { return Err("duplicate --workers".into()); }
                i+=1; workers=Some(args.get(i).ok_or("--workers needs a count in 1..128")?.parse::<usize>()?);
            }
            _ => remaining.push(args[i].clone()),
        }
        i+=1;
    }
    let method=match population {
        Some((maximum_w,attempts)) => {
            let settings=PopulationOptions {maximum_w,attempts,seed:seed.unwrap_or(17),
                objective_target:target.unwrap_or(1e-8),workers};
            settings.validate()?;
            Method::Population(settings)
        }
        None => {
            if seed.is_some() || target.is_some() || workers.is_some() {
                return Err("--seed, --objective-target and --workers require --population".into());
            }
            Method::Sqp
        }
    };
    Ok((remaining,method))
}

#[derive(Debug)]
pub(super) struct PopulationStatus {
    pub stop: DifferentialEvolutionStop,
    pub attempts: usize,
    pub objective_target: f64,
    pub target_recovered: bool,
    pub settings: PopulationOptions,
}

pub(super) fn print_status(status: &PopulationStatus, objective: f64, violation: f64, evaluations: usize) {
    println!("# optimizer=population stop={:?} search_attempts={} physical_evaluations={} objective={objective:.17e} objective_target={:.17e} target_recovered={} max_scaled_violation={violation:.17e}; no_stationarity_claim=true",
        status.stop,status.attempts,evaluations,status.objective_target,status.target_recovered);
    println!("# population_maximum_w={:.17e} seed={} evaluation_mode={} requested_workers={} maximum_parallel_wave=24",
        status.settings.maximum_w,status.settings.seed,
        if status.settings.workers.is_some() {"parallel"} else {"serial"},status.settings.workers.unwrap_or(1));
}

pub(super) struct Outcome {
    pub sources: Vec<MovingHeatSource>,
    pub value: f64,
    pub predictions: Vec<f64>,
    pub residuals: Vec<f64>,
    pub constraints: ScanConstraintValues,
    pub evaluations: usize,
    /// Only SQP supplies a constrained-stationarity test. Population success
    /// means feasible target attainment, never derivative-based convergence.
    pub sqp_converged: bool,
    pub population: Option<PopulationStatus>,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve<'h, 'm, 'c, F, O: ScanObservationSpec + Sync>(cx: &Cx<'_>,
    program: &ScanHeatingProgram, mut study: ScanPowerStudy<F, O>,
    required: BoundScanRequirements, method: Method, checkpoints: ScanCheckpointBudget)
    -> Result<Outcome>
where F: Fn(&Cx<'_>, Vec<MovingHeatSource>)
    -> std::result::Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError> + Sync
{
    match method {
        Method::Sqp => {
            study.set_checkpointing(cx, checkpoints)?;
            let mut study = ConstrainedScanStudy::new_with_actuator_bounds(cx, study,
                required.powers, required.energy, required.states, required.peaks,
                required.response_budget, required.max_peak_samples, required.actuators,
                required.kkt_dimension.max(8))?;
            let report = study.advance_optimizer(cx, FEASIBILITY_TOLERANCE, 100)?;
            let sources = study.sources(cx)?;
            study.extend_evaluation_budget(128)?;
            let (replay, limits) = study.replay_with_constraints(cx)?;
            if (report.solution.f - replay.physical.value).abs() > 1e-10 * (1.0 + report.solution.f.abs()) {
                return Err("replayed sources changed the accepted fitting objective".into());
            }
            if limits.residuals.len() != study.sample().ci.len()
                || limits.residuals.iter().zip(&study.sample().ci).any(|(a,b)| a.to_bits() != b.to_bits()) {
                return Err("physical replay changed the accepted constraint rows".into());
            }
            Ok(Outcome { sources, value: replay.physical.value,
                predictions: replay.physical.predictions, residuals: replay.physical.standardized_residuals,
                constraints: values(limits), evaluations: study.evaluations(),
                sqp_converged: report.solution.converged, population: None })
        }
        Method::Population(options) => {
            let bounds = options.bounds(cx, program)?;
            let dimension = program.controls().len();
            let budget = ScanValueConstraintBudget {
                max_rows: required.kkt_dimension.checked_sub(dimension).ok_or("invalid constraint row allowance")?,
                max_energy_weights: required.energy.len().checked_mul(dimension).ok_or("energy allowance overflow")?,
                responses: ScanValueResponseBudget { max_responses: required.response_budget.max_responses,
                    max_peak_samples: required.max_peak_samples },
            };
            let study = ValueConstrainedScanStudy::new_with_actuator_bounds(cx, study,
                required.powers, required.energy, required.states, required.peaks, required.actuators, budget)?;
            let population_size = (4 * dimension).max(24);
            let mut study = PopulationScanStudy::new(cx, study, &bounds, DifferentialEvolutionConfig {
                population_size, differential_weight: 0.8, crossover_probability: 0.9,
                seed: options.seed, max_evaluations: options.attempts,
                max_population_entries: population_size * dimension, target: Some(options.objective_target),
            }, FEASIBILITY_TOLERANCE)?;
            let report = match options.workers {
                None => study.advance(cx, options.attempts)?,
                Some(workers) => {
                    let pool = TilePool::for_host(workers, options.seed);
                    // Executor metadata only: each original solver retains its
                    // own workspace/heap limits, not this lease as a total cap.
                    let lease = OperationMemoryLease::bounded(64 << 20);
                    study.advance_parallel(cx, &pool, &lease, options.attempts, 24)?
                }
            };
            // The owner has reserved one actual solve, separately from the DE
            // callback budget. Every declared physical constraint is replayed.
            let replay = study.replay_best(cx)?;
            let status = PopulationStatus { stop: report.search.stop,
                attempts: replay.search_evaluations, objective_target: options.objective_target,
                target_recovered: replay.target_recovered(), settings: options };
            Ok(Outcome { sources: replay.sources, value: replay.physical.physical.value,
                predictions: replay.physical.physical.predictions,
                residuals: replay.physical.physical.standardized_residuals,
                constraints: replay.constraints, evaluations: replay.physical_evaluations,
                sqp_converged: false, population: Some(status) })
        }
    }
}

// Discard derivatives only when presenting an SQP replay. Value-only results
// never acquire a fabricated zero Jacobian or fake transpose-work report.
fn values(v: ScanConstraintEvaluation) -> ScanConstraintValues {
    ScanConstraintValues { nonnegative_power_values_w: v.nonnegative_power_values_w,
        bound_values: v.bound_values, absorbed_energy_j: v.absorbed_energy_j,
        power_slew_values_w_s: v.power_slew_values_w_s, path_speed_values: v.path_speed_values,
        state_values: v.state_values, peak_values: v.peak_values,
        sampled_peak_values: v.sampled_peak_values, peak_smoothing_gaps: v.peak_smoothing_gaps,
        residuals: v.residuals }
}
