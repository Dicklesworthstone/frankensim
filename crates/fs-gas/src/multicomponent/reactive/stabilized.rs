//! RKL2 composition through the original stationary reacting-flow owner.

use super::*;
use crate::multicomponent::dissipation::{
    StabilizedDissipationControl, StabilizedDissipationReport,
};

/// Accepted molecular halves of one coupled RKL2/chemistry/advection step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitStabilizedDissipationReport {
    /// Original RKL2 work and physical balances before chemistry/advection.
    pub first: StabilizedDissipationReport,
    /// Original RKL2 work and physical balances after chemistry/advection.
    pub second: StabilizedDissipationReport,
}

/// One accepted step through the original coupled composition and ledgers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedReactiveStepReport {
    /// Original chemistry, advection, boundary and source accounting.
    pub coupled: ReactiveStepReport,
    /// Accepted RKL2 molecular-half diagnostics.
    pub dissipation: SplitStabilizedDissipationReport,
    /// All actual RHS calls, including discarded coupled attempts and stages.
    pub operator_evaluations: usize,
}

/// Complete-duration work bounds for stationary reacting molecular evolution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedReactiveIntervalControl {
    /// Original acoustic, reconstruction and whole-step retry controls.
    pub step: ReactiveStepControl,
    /// Original RKL2 stage and macrostep controls. Its RHS budget is shared by
    /// the entire interval, both molecular halves and every rejected trial.
    pub dissipation: StabilizedDissipationControl,
    /// Total coupled attempts, including restarts, in [1,1,000,000]. This is
    /// distinct from molecular stages and the original per-cell chemistry work.
    pub maximum_coupled_attempts: usize,
}
impl Default for StabilizedReactiveIntervalControl {
    fn default() -> Self {
        Self {
            step: ReactiveStepControl::default(),
            dissipation: StabilizedDissipationControl::default(),
            maximum_coupled_attempts: 10_000,
        }
    }
}
impl StabilizedReactiveIntervalControl {
    /// Admit all time-integration and aggregate work controls before evolution.
    pub fn validate(self) -> Result<(), ReactiveError> {
        self.step.validate()?;
        self.dissipation.validate()?;
        if !(1..=1_000_000).contains(&self.maximum_coupled_attempts) {
            return Err(ReactiveError::Invalid("reacting interval attempt budget"));
        }
        Ok(())
    }
}

/// Accepted complete-window work and endpoint balances, without a step trace.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedReactiveIntervalReport {
    /// Compensated sum of actual accepted physical step durations.
    pub elapsed: Time,
    /// Original physical clock.
    pub initial_time: Time,
    /// Exact represented target reached by an actually integrated last step.
    pub final_time: Time,
    /// Completed coupled steps.
    pub coupled_steps: usize,
    /// Discarded complete coupled attempts before those accepted steps.
    pub coupling_retries: usize,
    /// Every actual molecular RHS call, including rejected/abandoned work.
    pub operator_evaluations: usize,
    /// Largest RKL2 stage count over accepted molecular halves.
    pub maximum_stages_used: usize,
    /// Largest accepted RKL2 parabolic number divided by its stage extension.
    pub maximum_scaled_parabolic_number: f64,
    /// Accepted chemistry work and source changes, summed over the window.
    pub chemistry: ChemistryReport,
    /// Time/area-integrated conservative advection through both boundaries.
    pub advection_boundary_change: Totals,
    /// Original molecular species/enthalpy/stress boundary transfers.
    pub molecular_boundary_change: Totals,
    /// Combined advection and molecular transfer into the left end.
    pub left_boundary_change: Totals,
    /// Combined advection and molecular transfer into the right end.
    pub right_boundary_change: Totals,
    /// Original accepted area-wall pressure impulses.
    pub area_momentum_impulse: f64,
    /// Original field integrals.
    pub before: Totals,
    /// Final field integrals.
    pub after: Totals,
    /// Endpoint difference minus original boundary and chemical source changes.
    pub conservation_residual: Totals,
    /// Total mass difference minus boundary transport [kg].
    pub mass_residual_kg: f64,
    /// Physical total-energy balance including original species references [J].
    pub absolute_energy_residual_j: f64,
}

pub(super) enum MolecularMethod<'a> {
    Explicit(DissipationControl),
    Stabilized {
        controls: StabilizedDissipationControl,
        evaluations: &'a mut usize,
    },
}
pub(super) struct MolecularHalf {
    physical: DissipationReport,
    stabilized: Option<StabilizedDissipationReport>,
}
impl MolecularMethod<'_> {
    pub(super) fn advance<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        flow: &mut MixtureEuler1d,
        duration: Time,
        reservoirs: [Option<Cell<'_>>; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularHalf, ReactiveError> {
        match self {
            Self::Explicit(controls) => Ok(MolecularHalf {
                physical: operator.advance_with_reservoirs_for(
                    flow, duration, *controls, reservoirs, cancelled,
                )?,
                stabilized: None,
            }),
            Self::Stabilized {
                controls,
                evaluations,
            } => {
                let report = operator.advance_stabilized_with_reservoirs_budget(
                    flow,
                    duration,
                    *controls,
                    reservoirs,
                    evaluations,
                    cancelled,
                )?;
                Ok(MolecularHalf {
                    physical: report.transport,
                    stabilized: Some(report),
                })
            }
        }
    }
}
pub(super) fn split_halves(
    first: Option<MolecularHalf>,
    second: Option<MolecularHalf>,
) -> (
    Option<SplitDissipationReport>,
    Option<SplitStabilizedDissipationReport>,
) {
    match first.zip(second) {
        Some((first, second)) => (
            Some(SplitDissipationReport {
                first: first.physical,
                second: second.physical,
            }),
            first
                .stabilized
                .zip(second.stabilized)
                .map(|(first, second)| SplitStabilizedDissipationReport { first, second }),
        ),
        None => (None, None),
    }
}

pub(super) struct CouplingBudget {
    maximum: usize,
    attempts: usize,
    completed: usize,
    start: f64,
}
impl CouplingBudget {
    pub(super) fn reserve(&mut self, time: f64) -> Result<(), ReactiveError> {
        if self.attempts >= self.maximum {
            return Err(ReactiveError::IntervalWorkLimit {
                elapsed_time_s: time - self.start,
                completed_steps: self.completed,
                attempts: self.attempts,
            });
        }
        self.attempts += 1;
        Ok(())
    }
}

impl ReactingEuler1d {
    /// Opt into RKL2 for both molecular halves of the original D-C-A-C-D step.
    ///
    /// Acoustic/heat-release restarts replay the complete original composition.
    /// The molecular RHS budget counts both halves and every discarded retry;
    /// exhaustion is terminal, so shortening never creates a fresh work budget.
    /// The original chemistry method, source interval and balance owner remain
    /// in control. Symmetric splitting alone is no global accuracy or entropy
    /// certificate; RKL2 retains its documented nonlinear-stability limitations.
    pub fn step_with_stabilized_dissipation<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveStepReport, ReactiveError> {
        poll(&mut cancelled)?;
        dissipation_controls.validate()?;
        self.stabilized_step(
            operator,
            maximum_duration,
            controls,
            dissipation_controls,
            [None; 2],
            &mut 0,
            None,
            &mut cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn stabilized_step<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        reservoirs: [Option<Cell<'_>>; 2],
        evaluations: &mut usize,
        budget: Option<&mut CouplingBudget>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveStepReport, ReactiveError> {
        let before = *evaluations;
        let (coupled, dissipation) = self.step_impl(
            maximum_duration,
            controls,
            Some(operator),
            MolecularMethod::Stabilized {
                controls: dissipation_controls,
                evaluations,
            },
            reservoirs,
            budget,
            cancelled,
        )?;
        Ok(StabilizedReactiveStepReport {
            coupled,
            dissipation: dissipation.expect("accepted stabilized composition has both halves"),
            operator_evaluations: *evaluations - before,
        })
    }

    /// Complete a nonnegative duration through the original coupled step owner.
    ///
    /// Reject a duration whose addition to the current clock changes its span
    /// by more than 1e-10 relative to that duration. This prevents a coarse
    /// absolute clock from silently changing the requested physical evolution.
    /// Use [`Self::advance_stabilized_to`] to explicitly select a target time.
    pub fn advance_stabilized_for<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        duration: Time,
        controls: StabilizedReactiveIntervalControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveIntervalReport, ReactiveError> {
        self.advance_stabilized_for_impl(operator, duration, controls, [None; 2], cancelled)
    }

    pub(super) fn advance_stabilized_for_impl<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        duration: Time,
        controls: StabilizedReactiveIntervalControl,
        reservoirs: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveIntervalReport, ReactiveError> {
        poll(&mut cancelled)?;
        let start = self.accepted.time;
        let target = start + duration.value();
        if !duration.value().is_finite()
            || duration.value() < 0.0
            || !target.is_finite()
            || !same_duration(target - start, duration.value())
        {
            return Err(ReactiveError::Invalid(
                "reacting duration is not representable at this clock",
            ));
        }
        self.advance_stabilized_to_impl(
            operator,
            Time::new(target),
            controls,
            reservoirs,
            cancelled,
        )
    }

    /// Reach an explicit represented target using molecular/chemistry/advection.
    ///
    /// Each terminal remainder is actually solved. No clock snap, epsilon-sized
    /// skipped interval or interpolation is used. Actual step durations and all
    /// source/boundary transfers accumulate in fixed order with compensation.
    /// The entire window restores the original field and clock on any late
    /// cancellation, source refusal, stagnant clock or exhausted work budget.
    ///
    /// A single conserved-state backup is lazily allocated and reused. Geometry
    /// is stationary; this method does not invoke moving-mesh/piston owners.
    /// Caller-selected macrostep and chemistry tolerances govern accuracy;
    /// this composition provides no global temporal error estimator.
    pub fn advance_stabilized_to<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        target: Time,
        controls: StabilizedReactiveIntervalControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveIntervalReport, ReactiveError> {
        self.advance_stabilized_to_impl(operator, target, controls, [None; 2], cancelled)
    }

    pub(super) fn advance_stabilized_to_impl<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        target: Time,
        controls: StabilizedReactiveIntervalControl,
        reservoirs: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveIntervalReport, ReactiveError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        let start = self.accepted.time;
        if !target.value().is_finite() || target.value() < start {
            return Err(ReactiveError::Invalid(
                "finite nondecreasing reacting target required",
            ));
        }
        reservoir::admit(&self.accepted, reservoirs, &mut cancelled)?;
        if self.interval_initial.len() != self.accepted.states.len() {
            self.interval_initial = allocated(self.accepted.states.len(), 0.0)?;
        }
        self.interval_initial.copy_from_slice(&self.accepted.states);
        let before = self.accepted.totals();
        let n = self.accepted.gas.species_count();
        let mut budget = CouplingBudget {
            maximum: controls.maximum_coupled_attempts,
            attempts: 0,
            completed: 0,
            start,
        };
        let mut evaluations = 0;
        let mut ledger = IntervalLedger::new();
        let outcome = (|| {
            if target.value() == start {
                // A no-op still admits the exact original operator and source.
                operator.advance_stabilized_with_reservoirs_budget(
                    &mut self.accepted,
                    Time::new(0.0),
                    controls.dissipation,
                    reservoirs,
                    &mut evaluations,
                    &mut cancelled,
                )?;
            }
            while self.accepted.time < target.value() {
                poll(&mut cancelled)?;
                let previous_time = self.accepted.time;
                let remaining = target.value() - previous_time;
                let step = self.stabilized_step(
                    operator,
                    Time::new(remaining),
                    controls.step,
                    controls.dissipation,
                    reservoirs,
                    &mut evaluations,
                    Some(&mut budget),
                    &mut cancelled,
                )?;
                if !positive(step.coupled.elapsed.value())
                    || step.coupled.elapsed.value() > remaining
                    || self.accepted.time <= previous_time
                    || self.accepted.time > target.value()
                {
                    return Err(ReactiveError::Invalid(
                        "reacting interval clock failed to progress within target",
                    ));
                }
                budget.completed += 1;
                ledger.add(step, n)?;
            }
            if !same_duration(ledger.elapsed.value, target.value() - start) {
                return Err(ReactiveError::Invalid(
                    "reacting clock span differs from integrated durations",
                ));
            }
            let after = self.accepted.totals();
            let complete = ledger.finish(
                &self.accepted.gas,
                before,
                after,
                start,
                target,
                &budget,
                evaluations,
            )?;
            poll(&mut cancelled)?;
            Ok::<_, ReactiveError>(complete)
        })();
        if outcome.is_err() {
            self.accepted.states.copy_from_slice(&self.interval_initial);
            self.accepted.time = start;
        }
        outcome
    }
}

// Admission of represented clock arithmetic, not a temporal error estimator.
const MAXIMUM_CLOCK_ROUNDOFF_FRACTION: f64 = 1e-10;

fn same_duration(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= MAXIMUM_CLOCK_ROUNDOFF_FRACTION * b.abs()
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let y = value - self.correction;
        let next = self.value + y;
        self.correction = (next - self.value) - y;
        self.value = next;
    }
}
#[derive(Clone, Copy)]
struct Ledger([Sum; MAX_COMPONENTS]);
impl Ledger {
    fn add(&mut self, total: Totals, n: usize) {
        for k in 0..n + 2 {
            self.0[k].add(total.get(k));
        }
    }
    fn total(&self, n: usize) -> Result<Totals, ReactiveError> {
        let mut values = [0.0; MAX_COMPONENTS];
        for (k, value) in values[..n + 2].iter_mut().enumerate() {
            *value = self.0[k].value;
            if !value.is_finite() {
                return Err(ReactiveError::Invalid(
                    "non-finite reacting interval ledger",
                ));
            }
        }
        Ok(Totals::from_values(&values[..n + 2], n))
    }
}
struct IntervalLedger {
    elapsed: Sum,
    // Advection, molecular, combined left, combined right.
    boundaries: [Ledger; 4],
    chemistry_sources: Ledger,
    entropy: Sum,
    area_impulse: Sum,
    chemistry: ChemistryReport,
    maximum_stages: usize,
    maximum_scaled_parabolic: f64,
}
impl IntervalLedger {
    fn new() -> Self {
        Self {
            elapsed: Sum::default(),
            boundaries: [Ledger([Sum::default(); MAX_COMPONENTS]); 4],
            chemistry_sources: Ledger([Sum::default(); MAX_COMPONENTS]),
            entropy: Sum::default(),
            area_impulse: Sum::default(),
            chemistry: ChemistryReport::default(),
            maximum_stages: 0,
            maximum_scaled_parabolic: 0.0,
        }
    }
    fn add(&mut self, step: StabilizedReactiveStepReport, n: usize) -> Result<(), ReactiveError> {
        let report = step.coupled;
        self.elapsed.add(report.elapsed.value());
        self.boundaries[0].add(report.transport.boundary_change, n);
        self.boundaries[2].add(report.transport.left_boundary_change, n);
        self.boundaries[3].add(report.transport.right_boundary_change, n);
        self.area_impulse
            .add(report.transport.area_momentum_impulse);
        for half in [step.dissipation.first, step.dissipation.second] {
            self.boundaries[1].add(half.transport.boundary_change, n);
            self.boundaries[2].add(half.transport.left_boundary_change, n);
            self.boundaries[3].add(half.transport.right_boundary_change, n);
            self.maximum_stages = self.maximum_stages.max(half.maximum_stages_used);
            self.maximum_scaled_parabolic = self
                .maximum_scaled_parabolic
                .max(half.maximum_scaled_parabolic_number);
        }
        let c = report.chemistry;
        for (total, increment) in [
            (&mut self.chemistry.cell_solves, c.cell_solves),
            (&mut self.chemistry.accepted_steps, c.accepted_steps),
            (&mut self.chemistry.rejected_steps, c.rejected_steps),
            (
                &mut self.chemistry.nonlinear_iterations,
                c.nonlinear_iterations,
            ),
        ] {
            *total = total.checked_add(increment).ok_or(ReactiveError::Invalid(
                "reacting interval chemistry work overflow",
            ))?;
        }
        self.chemistry.maximum_energy_residual_j = self
            .chemistry
            .maximum_energy_residual_j
            .max(c.maximum_energy_residual_j);
        self.chemistry.maximum_reconstruction_residual_j = self
            .chemistry
            .maximum_reconstruction_residual_j
            .max(c.maximum_reconstruction_residual_j);
        self.entropy.add(c.entropy_change_j_k);
        for k in 0..n {
            self.chemistry_sources.0[k].add(c.species_change_kg[k]);
        }
        self.chemistry_sources.0[n + 1].add(c.referenced_energy_change_j);
        Ok(())
    }
    #[allow(clippy::too_many_arguments)]
    fn finish(
        mut self,
        gas: &super::super::MixtureGas,
        before: Totals,
        after: Totals,
        start: f64,
        target: Time,
        budget: &CouplingBudget,
        evaluations: usize,
    ) -> Result<StabilizedReactiveIntervalReport, ReactiveError> {
        let n = gas.species_count();
        let advection = self.boundaries[0].total(n)?;
        let molecular = self.boundaries[1].total(n)?;
        let left = self.boundaries[2].total(n)?;
        let right = self.boundaries[3].total(n)?;
        let sources = self.chemistry_sources.total(n)?;
        self.chemistry.species_change_kg = sources.species_masses;
        self.chemistry.referenced_energy_change_j = sources.energy;
        self.chemistry.entropy_change_j_k = self.entropy.value;
        let (residual, mass_residual, energy_residual) = coupled_balance(
            gas,
            before,
            after,
            advection,
            molecular,
            self.area_impulse.value,
            &self.chemistry,
        )?;
        Ok(StabilizedReactiveIntervalReport {
            elapsed: Time::new(self.elapsed.value),
            initial_time: Time::new(start),
            final_time: target,
            coupled_steps: budget.completed,
            coupling_retries: budget.attempts - budget.completed,
            operator_evaluations: evaluations,
            maximum_stages_used: self.maximum_stages,
            maximum_scaled_parabolic_number: self.maximum_scaled_parabolic,
            chemistry: self.chemistry,
            advection_boundary_change: advection,
            molecular_boundary_change: molecular,
            left_boundary_change: left,
            right_boundary_change: right,
            area_momentum_impulse: self.area_impulse.value,
            before,
            after,
            conservation_residual: residual,
            mass_residual_kg: mass_residual,
            absolute_energy_residual_j: energy_residual,
        })
    }
}
