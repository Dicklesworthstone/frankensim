//! Step-doubling error control on the existing closed electrothermal step.
//!
//! Compare one full step with two half steps in actual cell joules, not Kelvin.
//! Accept the two conservative half steps without Richardson extrapolation.
//! Coarse trajectories, rejected halves and failed Newton trials spend compute
//! but never contribute heat, source work, or an accepted physical history.

use super::*;

/// Interval-wide controls. No tolerance is inferred from a material name.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveJouleConfig {
    /// Existing nonlinear policy. Absolute Newton and physical energy allowances
    /// are apportioned by each attempted step's fraction of the full interval.
    pub step: ImplicitJouleConfig,
    /// Maximum span of a full-step/two-half-step comparison [s].
    pub maximum_step_s: f64,
    /// Minimum comparison span allowed after rejection [s]. An initial/final
    /// interval shorter than this may be tried once, but cannot be subdivided.
    pub minimum_step_s: f64,
    /// Maximum calls to the implicit owner, including coarse/failed/rejected work.
    pub max_step_attempts: usize,
    /// Per-cell budget for the SUM of accepted step-doubling indicators [J].
    /// Each comparison receives dt/interval of this allowance. This is a local
    /// error estimate, not a certified global-error or spatial-error bound.
    pub cell_error_budget_j: f64,
}

/// Why an entire proposed interval was refused. No partial state accompanies it.
#[derive(Clone, Debug, PartialEq)]
pub enum AdaptiveJouleError {
    /// Invalid physical-time, accuracy, or work controls.
    Invalid(&'static str),
    /// Original material/current/thermal/Newton refusal, never replaced by data.
    Step(ImplicitJouleError),
    /// Nonlinear work or error estimation needs a smaller unadmitted timestep.
    StepFloor {
        /// Start of the rejected comparison [s].
        time_s: f64,
        /// Rejected span [s].
        span_s: f64,
        /// Original nonconvergence, or None for a temporal-error rejection.
        nonlinear: Option<ImplicitJouleError>,
    },
    /// The interval exhausted its declared number of implicit attempts.
    AttemptBudget,
    /// Accepted substep accounting did not close over the entire interval.
    Balance(&'static str),
}

/// Work survives failure for diagnosis, but tentative physical history does not.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveJouleFailure {
    /// Original refusal or independent interval gate.
    pub cause: AdaptiveJouleError,
    /// Calls already spent, including discarded and failed steps.
    pub step_attempts: usize,
    /// Current solves already spent, including numerical Jacobian probes.
    pub current_solves: usize,
}
impl core::fmt::Display for AdaptiveJouleFailure {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive electrothermal: {:?} ({} steps, {} current solves)",
            self.cause, self.step_attempts, self.current_solves)
    }
}
impl std::error::Error for AdaptiveJouleFailure {}

/// Accepted fine trajectory's energy accounts; no rejected work is physical heat.
#[derive(Debug)]
pub struct AdaptiveJouleReport<F, I, J> {
    /// Last accepted half step, including its fresh endpoint material/current
    /// evaluations. Its energy fields describe ONLY that half step.
    pub last_step: ImplicitJouleReport<F, I, J>,
    /// Chronological accepted half-step endpoints, bounded by max_step_attempts.
    pub accepted_times_s: Vec<f64>,
    /// Coarse, fine and failed calls to the implicit owner.
    pub step_attempts: usize,
    /// All current work, not just the final accepted trajectory.
    pub current_solves: usize,
    /// Comparisons discarded for temporal error or nonlinear nonconvergence.
    pub rejected_intervals: usize,
    /// Newton iterations on the accepted fine trajectory, excluding trial work.
    pub accepted_newton_iterations: usize,
    /// Largest physical cell residual among the accepted half steps [J].
    pub maximum_cell_residual_j: f64,
    /// Sum of accepted indicators per cell, in stable thermal CellId order [J].
    pub estimated_cell_errors_j: Vec<(CellId, f64)>,
    /// Total accepted boundary/additional/Joule heat, counted once [J].
    pub external_heat_j: f64,
    /// Total accepted electrical source work [J].
    pub electrode_energy_j: f64,
    /// Total accepted irreversible electrical heat [J].
    pub joule_energy_j: f64,
    /// Accepted Joule heat actually handed to receiving cells [J].
    pub deposited_joule_j: f64,
    /// Actual final minus initial cell enthalpy storage [J].
    pub stored_energy_change_j: f64,
    /// Actual storage change minus accepted external heat [J].
    pub energy_residual_j: f64,
}

/// Advance a fixed-topology, continuous-drive interval with material-sensitive
/// temporal error control. The original fixed-step API remains unchanged.
///
/// All four material/source callbacks have the same semantics as
/// `step_joule_enthalpy_implicit`. `initial_guess(state, dt)` supplies a PURE
/// joule-increment guess in sorted thermal CellId order for EACH attempted step;
/// it must not debit any source. Trial and accepted states use identical laws.
/// The initial state's next scheduled event may equal target, but not precede
/// it. Split explicitly at drive changes and events; target-time events remain
/// pending. Topology, reference mass, material charts and geometry are unchanged.
///
/// Only Newton nonconvergence and excessive step-doubling error trigger halving.
/// Missing/unsupported material, failed physical acceptance, cancellation and
/// work exhaustion refuse immediately. The maximum current-solve budget in
/// `input` covers the entire call, including all rejected candidates. Every
/// return is transactional: late failure exposes neither an earlier accepted
/// prefix nor any emitted observation. Callbacks must not mutate external state.
///
/// The error indicator m*abs(h_two_halves-h_full) includes latent storage even
/// when temperatures are identical. Accepted half-step results are NOT
/// extrapolated, preserving the existing local energy gates and source records.
/// Error-estimate refinement and spatial/physical validation remain separate.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_joule_enthalpy_adaptive<'t, F, I, J>(
    thermal: &DepositionSimulation<'t>, target_s: f64,
    input: JouleEnthalpyInput<'_>, link_ids: &[u64], gate: &CancelGate,
    config: AdaptiveJouleConfig,
    initial_guess: impl Fn(&DepositionSimulation<'t>, f64) -> Result<Vec<f64>, JouleEnthalpyError>,
    resolve: impl Fn(usize, &DepositedCell, f64) -> Result<EnthalpyTransport<F>, JouleEnthalpyError>,
    resolve_link: impl Fn(&ContactPatch, [&EnthalpyTransportSample<F>; 2]) -> Result<(f64,I), JouleEnthalpyError>,
    resolve_interface: impl Fn(usize, &ElectricalContact, [&EnthalpyTransportSample<F>; 2])
        -> Result<ElectricalInterfaceResponse<J>, JouleEnthalpyError>,
    external_sources: impl Fn(&DepositionSimulation<'t>) -> Result<BTreeMap<CellId,f64>, JouleEnthalpyError>,
) -> Result<(DepositionSimulation<'t>, AdaptiveJouleReport<F,I,J>), AdaptiveJouleFailure> {
    let attempts = Cell::new(0usize);
    let work = Cell::new(0usize);
    let fail = |cause| AdaptiveJouleFailure { cause, step_attempts:attempts.get(), current_solves:work.get() };
    let physical = |e: JouleEnthalpyError| fail(AdaptiveJouleError::Step(ImplicitJouleError::Physical(e)));
    let thermal_failure = |e: DepositionError| physical(e.into());
    checkpoint(gate).map_err(&thermal_failure)?;
    let duration = target_s-thermal.time_s();
    if !duration.is_finite() || duration <= 0.0
        || !config.maximum_step_s.is_finite() || config.maximum_step_s <= 0.0
        || !config.minimum_step_s.is_finite() || config.minimum_step_s <= 0.0
        || config.minimum_step_s > config.maximum_step_s || config.max_step_attempts < 3
        || config.step.max_unknowns == 0 || thermal.cells().len() > config.step.max_unknowns
        || !config.cell_error_budget_j.is_finite() || config.cell_error_budget_j <= 0.0 {
        return Err(fail(AdaptiveJouleError::Invalid("invalid interval, step bounds, attempt count or cell error budget")));
    }
    config.step.newton.validate().map_err(|e|fail(AdaptiveJouleError::Step(e.into())))?;
    // Reject pending events and inconsistent maps BEFORE any numerical trial.
    BackwardEulerStep::new(thermal,target_s,||checkpoint(gate)).map_err(|e|physical(e.into()))?;
    admit_map(thermal,input,gate).map_err(&physical)?;
    admit_interfaces(input,gate).map_err(&physical)?;
    let attempt = |state: &DepositionSimulation<'t>, end: f64| {
        checkpoint(gate).map_err(&thermal_failure)?;
        if attempts.get() >= config.max_step_attempts { return Err(fail(AdaptiveJouleError::AttemptBudget)); }
        let remaining = input.max_current_solves-work.get();
        if remaining == 0 { return Err(physical(JouleEnthalpyError::CurrentBudget)); }
        attempts.set(attempts.get()+1);
        let dt = end-state.time_s();
        let fraction = dt/duration;
        let mut step_config = config.step;
        step_config.cell_energy_tolerance_j *= fraction;
        step_config.total_energy_tolerance_j *= fraction;
        step_config.newton.absolute_tolerance *= fraction;
        let guess = initial_guess(state,dt).map_err(&physical)?;
        let spent = Cell::new(0);
        let result = step_counted(state,end,JouleEnthalpyInput {max_current_solves:remaining,..input},
            link_ids,gate,step_config,&guess,&resolve,&resolve_link,&resolve_interface,&external_sources,&spent);
        work.set(work.get()+spent.get()); // the inner owner cannot exceed remaining
        result.map_err(|e|fail(AdaptiveJouleError::Step(e)))
    };
    let mut state = thermal.clone();
    let mut span = config.maximum_step_s.min(duration);
    let mut accepted_times_s = Vec::new();
    let mut errors: Vec<_> = thermal.cells().keys().map(|&id|(id,0.0_f64)).collect();
    let mut last = None;
    let mut rejected = 0usize;
    let mut accepted_newton_iterations = 0usize;
    let mut maximum_cell_residual_j = 0.0_f64;
    let (mut external, mut electrode, mut joule, mut deposited) = (0.0,0.0,0.0,0.0);
    while state.time_s() < target_s {
        checkpoint(gate).map_err(&thermal_failure)?;
        let start = state.time_s();
        let remaining = target_s-start;
        let end = if span >= remaining {target_s} else {start+span};
        let dt = end-start;
        let mid = start+0.5*dt;
        if !(start < mid && mid < end && end <= target_s) {
            return Err(fail(AdaptiveJouleError::Invalid("time interval cannot represent two ordered half steps")));
        }
        let comparison: Result<_, AdaptiveJouleFailure> = (|| {
            let (coarse,_) = attempt(&state,end)?;
            let (half,r1) = attempt(&state,mid)?;
            let (fine,r2) = attempt(&half,end)?;
            Ok((coarse,fine,r1,r2))
        })();
        let mut nonlinear = None;
        match comparison {
            Ok((coarse,fine,r1,r2)) => {
                let allowance = config.cell_error_budget_j*(dt/duration);
                if !allowance.is_finite() || allowance <= 0.0 {
                    return Err(fail(AdaptiveJouleError::Invalid("unrepresentable cell error allowance")));
                }
                let mut indicators = Vec::with_capacity(errors.len());
                let mut ratio = 0.0_f64;
                for (i,(id,_)) in errors.iter().enumerate() {
                    if i%256 == 0 {checkpoint(gate).map_err(&thermal_failure)?;}
                    let a = &coarse.cells()[id]; let b = &fine.cells()[id];
                    let difference = finite(a.mass_kg*(b.specific_h_j_kg-a.specific_h_j_kg).abs()).map_err(&physical)?;
                    ratio = ratio.max(difference/allowance);
                    indicators.push(difference);
                }
                if ratio <= 1.0 {
                    for ((_,sum),value) in errors.iter_mut().zip(indicators) { *sum=finite(*sum+value).map_err(&physical)?; }
                    external=finite(external+r1.thermal.external_heat_j+r2.thermal.external_heat_j).map_err(&physical)?;
                    electrode=finite(electrode+r1.electrode_energy_j+r2.electrode_energy_j).map_err(&physical)?;
                    joule=finite(joule+r1.joule_energy_j+r2.joule_energy_j).map_err(&physical)?;
                    deposited=finite(deposited+r1.deposited_joule_j+r2.deposited_joule_j).map_err(&physical)?;
                    accepted_newton_iterations = accepted_newton_iterations
                        .checked_add(r1.nonlinear.iterations).and_then(|n|n.checked_add(r2.nonlinear.iterations))
                        .ok_or_else(||fail(AdaptiveJouleError::Balance("accepted iteration counter overflow")))?;
                    maximum_cell_residual_j = maximum_cell_residual_j
                        .max(r1.maximum_cell_residual_j).max(r2.maximum_cell_residual_j);
                    accepted_times_s.extend([mid,end]);
                    state=fine; last=Some(r2);
                    span=if ratio<=0.125 {dt.min(config.maximum_step_s/2.0)*2.0} else {dt};
                    continue;
                }
            }
            Err(AdaptiveJouleFailure {cause:AdaptiveJouleError::Step(e@ImplicitJouleError::NotConverged(_)),..}) => {
                nonlinear=Some(e);
            }
            Err(e) => return Err(e),
        }
        rejected+=1;
        span=0.5*dt;
        if span < config.minimum_step_s {
            return Err(fail(AdaptiveJouleError::StepFloor {time_s:start,span_s:dt,nonlinear}));
        }
    }
    let mut stored = 0.0;
    for (i,(id,initial)) in thermal.cells().iter().enumerate() {
        if i%256 == 0 {checkpoint(gate).map_err(&thermal_failure)?;}
        let delta = finite(initial.mass_kg*(state.cells()[id].specific_h_j_kg-initial.specific_h_j_kg)).map_err(&physical)?;
        stored=finite(stored+delta).map_err(&physical)?;
    }
    let residual=finite(stored-external).map_err(&physical)?;
    if residual.abs()>config.step.total_energy_tolerance_j
        || (electrode-joule).abs()>input.energy_tolerance_j || (joule-deposited).abs()>input.energy_tolerance_j {
        return Err(fail(AdaptiveJouleError::Balance("whole-interval storage/source/transfer balance")));
    }
    // Allow only floating summation roundoff in the partitioned indicator budget.
    if errors.iter().any(|(_,e)|*e>config.cell_error_budget_j*(1.0+64.0*f64::EPSILON*accepted_times_s.len() as f64)) {
        return Err(fail(AdaptiveJouleError::Balance("whole-interval temporal indicator budget")));
    }
    checkpoint(gate).map_err(&thermal_failure)?;
    Ok((state,AdaptiveJouleReport {
        last_step:last.ok_or_else(||fail(AdaptiveJouleError::Invalid("empty accepted trajectory")))?,
        accepted_times_s,step_attempts:attempts.get(),current_solves:work.get(),
        rejected_intervals:rejected,estimated_cell_errors_j:errors,
        accepted_newton_iterations,maximum_cell_residual_j,
        external_heat_j:external,electrode_energy_j:electrode,joule_energy_j:joule,
        deposited_joule_j:deposited,stored_energy_change_j:stored,energy_residual_j:residual,
    }))
}
