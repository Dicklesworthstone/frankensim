//! Same-area independent-load descent using the existing multi-load kernel.
//!
//! Projection is applied BEFORE solving each candidate under every load. Neither
//! an overfilled starting strip nor an unevaluated proposal is an improvement
//! baseline. Worst-case selection is recomputed on the projected field, rather
//! than trusting the previously active load. Weights are not probabilities.

use super::*;
use crate::{RobustSampledStressEvaluation, SampledStressLimit};
use crate::volume::{
    VolumeProjectionReport, VolumeProjectionSettings, VolumeProjectionStage,
    project_material_volume, project_material_volume_controlled,
};

mod stress;
mod checkpoint;
mod restoration;
mod resolution;

/// Bounded candidate search and actual case-solve allowance.
#[derive(Debug, Clone, Copy)]
pub struct MultiLoadProjectedSettings {
    /// Maximum projected candidates per accepted update, in `1..=64`.
    pub max_candidates: usize,
    /// Interface-travel contraction after each rejected candidate, in `(0,1)`.
    pub contraction: f64,
    /// Strict relative decrease of the selected aggregate, in `[0,1)`.
    pub min_relative_improvement: f64,
    /// Total case solves that may start, including baseline, failures and retries.
    pub max_solves: usize,
}

impl Default for MultiLoadProjectedSettings {
    fn default() -> Self {
        Self {
            max_candidates: 6,
            contraction: 0.5,
            min_relative_improvement: 1e-8,
            max_solves: 4096,
        }
    }
}

/// Cooperative boundaries. `advance_one_polling` also polls inside candidate
/// CG solves. Assembly, smoothing and area evaluations remain indivisible;
/// neither entry point gives a hard wall-time guarantee.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiLoadProjectedStage {
    /// Before sensitivity construction or any candidate work for an update.
    Direction,
    /// Before evolving a candidate from the retained geometry.
    Prepare(usize),
    /// One existing numerical-area projection boundary.
    Projection(usize, VolumeProjectionStage),
    /// Before (`complete=false`) or after one independent case solve.
    CaseSolve { candidate: usize, case: usize, complete: bool },
    /// Candidate CG progress, including true-residual correction passes. Only
    /// emitted by `advance_one_polling`; iterations are cumulative per case.
    CaseIterations { candidate: usize, case: usize, iterations: usize },
    /// Before a cell's stress probes on an already solved independent load case.
    StressCell { candidate: usize, case: usize, cell: usize },
    /// Complete projected-field family and acceptance gates passed, before mutation.
    Publish(usize),
}

/// Numerical metrics bound to one exact, fully evaluated nodal field.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiLoadProjectedState {
    /// Selected weighted-sum or worst-weighted objective.
    pub objective: f64,
    /// Unweighted compliance of every input case, in declaration order.
    pub case_compliances: Vec<f64>,
    /// First maximizing case in worst-weighted mode; `None` for weighted sum.
    pub active_case: Option<usize>,
    /// Numerical cut-quadrature area, not a continuum enclosure.
    pub volume: f64,
    /// Fingerprint of the exact retained level-set bits.
    pub snapshot: u64,
}

impl MultiLoadProjectedState {
    fn of(state: &MultiState) -> Self {
        Self {
            objective: state.objective,
            case_compliances: state.compliances.clone(),
            active_case: state.active,
            volume: state.volume,
            snapshot: fnv(&state.phi),
        }
    }
}

/// One candidate's actual projection and complete-family evaluation, if any.
#[derive(Debug, Clone)]
pub struct MultiLoadProjectedAttempt {
    /// Zero-based candidate index within the current update.
    pub index: usize,
    /// Interface travel given to the existing evolution kernel.
    pub move_cells: f64,
    /// Numerical projection evidence, present only after successful projection.
    pub projection: Option<VolumeProjectionReport>,
    /// Never present for a partly solved scenario family.
    pub state: Option<MultiLoadProjectedState>,
    /// Complete-family sampled stress; absent when disabled or sampling refuses.
    pub stress: Option<RobustSampledStressEvaluation>,
    /// Numerical refusal or failed acceptance gate, not a convergence claim.
    pub refusal: Option<String>,
}

/// One accepted update; proposal diagnostics do not certify the projected field.
#[derive(Debug, Clone)]
pub struct MultiLoadProjectedStep {
    /// Global accepted-update ordinal.
    pub iteration: usize,
    /// This update reduces stress violation, not necessarily compliance.
    pub restoration: bool,
    /// Fully evaluated area-feasible previous state.
    pub previous: MultiLoadProjectedState,
    /// Fully evaluated same-area accepted state. During restoration its stress
    /// may still exceed the declared limit; inspect `restoration` and `stress`.
    pub state: MultiLoadProjectedState,
    /// Complete-family stress of the accepted state, when a limit is installed.
    pub stress: Option<RobustSampledStressEvaluation>,
    /// Projection applied before the accepted state's independent load solves.
    pub projection: VolumeProjectionReport,
    /// Unprojected evolution audit, deliberately separate from final metrics.
    pub proposal_audit: RedistanceAudit,
    /// Hole insertions attempted by the unprojected proposal.
    pub proposal_events: Vec<NucleationEvent>,
    /// Load-pad assignments made before projection restored fixed nodes.
    pub proposal_load_pad_nodes: usize,
    /// All attempted candidates through the accepted one.
    pub attempts: Vec<MultiLoadProjectedAttempt>,
}

/// Terminal labels express bounded work, not convergence or optimality.
#[derive(Debug, Clone)]
pub enum MultiLoadProjectedProgress {
    /// One accepted same-area state replaced the previous one.
    Accepted(Box<MultiLoadProjectedStep>),
    /// The requested accepted-update count has been reached.
    IterationLimit,
    /// The candidate family found no admissible strict aggregate decrease.
    NoDescent(Vec<MultiLoadProjectedAttempt>),
    /// Too few solves remain for another COMPLETE scenario family.
    SolveBudget(Vec<MultiLoadProjectedAttempt>),
}

/// Stateful projected descent. Geometry and cached independent PDE solutions
/// change atomically; spent work is intentionally NOT refunded on interruption.
/// Repeated `advance_one` calls retain the global nucleation ordinal and AL
/// search multiplier. No new elasticity, sensitivity or advection kernel exists.
pub struct MultiLoadProjectedOptimizer {
    kernel: Kernel,
    current: MultiState,
    baseline: MultiLoadProjectedState,
    baseline_geometry: GridSdf,
    fixed: Vec<(usize, f64)>,
    projection: VolumeProjectionSettings,
    controls: MultiLoadProjectedSettings,
    next_iteration: usize,
    ell: f64,
    solves_started: usize,
    stress_limit: Option<SampledStressLimit>,
    baseline_stress: Option<RobustSampledStressEvaluation>,
    current_stress: Option<RobustSampledStressEvaluation>,
    restoration_reduction: Option<f64>,
    restoration_updates: usize,
}

impl MultiLoadProjectedOptimizer {
    /// Restore numerical-area feasibility and solve every declared load before
    /// establishing an objective baseline. Inputs are consumed, never published
    /// as a successful optimizer when projection or any baseline case refuses.
    ///
    /// # Errors
    /// Refuses invalid controls, a baseline exceeding the solve allowance,
    /// unattainable fixed-node area, unsupported loads or a canonical PDE error.
    /// Baseline construction is synchronous and is not cancellable internally.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mut geometry: GridSdf,
        load_cases: &[RobustLoadCase],
        settings: OptimizeSettings,
        aggregate: RobustAggregate,
        fixed: Vec<(usize, f64)>,
        projection: VolumeProjectionSettings,
        controls: MultiLoadProjectedSettings,
    ) -> Result<Self, CutFemError> {
        if !(1..=8).contains(&settings.level)
            || !(1..=10_000).contains(&settings.iterations)
            || !(1..=64).contains(&load_cases.len())
            || !(1..=64).contains(&controls.max_candidates)
            || projection.target.to_bits() != settings.volfrac.to_bits()
            || !(settings.move_cells.is_finite() && settings.move_cells > 0.0)
            || !(controls.contraction.is_finite() && controls.contraction > 0.0
                && controls.contraction < 1.0)
            || !(controls.min_relative_improvement.is_finite()
                && (0.0..1.0).contains(&controls.min_relative_improvement))
            || controls.max_solves < load_cases.len()
        {
            return Err(invalid("projected multi-load descent needs level 1..=8, updates 1..=10000, cases/candidates 1..=64, equal area targets, positive travel, contraction (0,1), decrease [0,1), and a complete baseline solve allowance"));
        }
        // Check every evolution/material input before projection work.
        validate(&geometry, load_cases, settings)?;
        material(settings)?;
        project_material_volume(&mut geometry, settings.level, &fixed, projection)?;
        let kernel = Kernel::new(&geometry, load_cases, settings, aggregate)?;
        let current = kernel.evaluate(geometry)?;
        if (current.volume - projection.target).abs() > projection.tolerance {
            return Err(invalid("projected multi-load baseline failed the independent area gate"));
        }
        let baseline = MultiLoadProjectedState::of(&current);
        let baseline_geometry = current.phi.clone();
        Ok(Self {
            kernel, current, baseline, baseline_geometry, fixed, projection, controls,
            next_iteration: 0, ell: settings.ell0, solves_started: load_cases.len(),
            stress_limit: None, baseline_stress: None, current_stress: None,
            restoration_reduction: None, restoration_updates: 0,
        })
    }

    /// Exact accepted geometry, not the last attempted proposal.
    #[must_use]
    pub fn geometry(&self) -> &GridSdf { &self.current.phi }

    /// Same-material initial baseline, unchanged by subsequent updates.
    #[must_use]
    pub fn baseline(&self) -> &MultiLoadProjectedState { &self.baseline }

    /// Metrics of the exact accepted geometry under EVERY declared case.
    #[must_use]
    pub fn current(&self) -> MultiLoadProjectedState { MultiLoadProjectedState::of(&self.current) }

    /// Global ordinal of the next accepted update.
    #[must_use]
    pub const fn next_iteration(&self) -> usize { self.next_iteration }

    /// Actual case-solve attempts, including refused/cancelled candidates.
    #[must_use]
    pub const fn solves_started(&self) -> usize { self.solves_started }

    /// Immutable total case-solve allowance.
    #[must_use]
    pub const fn max_solves(&self) -> usize { self.controls.max_solves }

    /// Attempt one same-area aggregate-decreasing update.
    ///
    /// # Errors
    /// Propagates sensitivity or internal state refusal without changing geometry.
    pub fn advance_one(&mut self) -> Result<MultiLoadProjectedProgress, CutFemError> {
        match self.advance_one_controlled(|_| ControlFlow::<Infallible>::Continue(()))? {
            ControlFlow::Continue(progress) => Ok(progress),
            ControlFlow::Break(never) => match never {},
        }
    }

    /// Interruption never changes accepted geometry, ordinal or multiplier.
    /// Already started case solves remain charged. A complete family is reserved
    /// before an attempt; a late refused or interrupted case cannot publish an
    /// aggregate from the earlier subset, including zero-weight scenarios.
    ///
    /// # Errors
    /// Propagates sensitivity or internal state refusal; individual candidate
    /// refusals are retained in attempt records and permit bounded backtracking.
    pub fn advance_one_controlled<B>(
        &mut self,
        control: impl FnMut(MultiLoadProjectedStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, MultiLoadProjectedProgress>, CutFemError> {
        self.advance_one_scheduled(None, control)
    }

    /// Advance using the original CutFEM solver with at most `poll_iters`
    /// additional CG iterations between callbacks, including correction passes.
    /// Every previous projection, stress, case and publication boundary remains.
    /// Interrupted candidate fields never replace the accepted solution family;
    /// their attempted case solves remain charged exactly once per start.
    ///
    /// This is a per-call scheduling choice, not a new physical/study policy.
    /// It is not serialized; a checkpoint can resume using either entry point.
    /// Baseline construction and checkpoint recovery are still synchronous.
    /// Assembly, reductions, sparse applies, smoothing and area evaluations are
    /// not preemptible. No hard deadline or mid-CG checkpoint is promised.
    ///
    /// # Errors
    /// Refuses zero polling before work; otherwise follows the same numerical
    /// refusal and complete-family publication rules as `advance_one_controlled`.
    pub fn advance_one_polling<B>(
        &mut self,
        poll_iters: usize,
        control: impl FnMut(MultiLoadProjectedStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, MultiLoadProjectedProgress>, CutFemError> {
        if poll_iters == 0 {
            return Err(invalid("multi-load CG polling interval must be positive"));
        }
        self.advance_one_scheduled(Some(poll_iters), control)
    }

    fn advance_one_scheduled<B>(
        &mut self,
        poll_iters: Option<usize>,
        control: impl FnMut(MultiLoadProjectedStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, MultiLoadProjectedProgress>, CutFemError> {
        self.advance_one_admitted_scheduled(poll_iters, 1,
            |_, _, _, _, _| Ok(ControlFlow::Continue(None)), control)
    }

    // One proposal/projection/equilibrium/publication owner for ordinary and
    // mesh-checked descent. Additional solves use the SAME lifetime allowance.
    fn advance_one_admitted_scheduled<B, C>(
        &mut self,
        poll_iters: Option<usize>,
        levels_per_candidate: usize,
        mut admit: impl FnMut(usize, &MultiState, Option<&RobustSampledStressEvaluation>,
            &mut usize, &mut C) -> Result<ControlFlow<B, Option<String>>, CutFemError>,
        mut control: C,
    ) -> Result<ControlFlow<B, MultiLoadProjectedProgress>, CutFemError>
    where C: FnMut(MultiLoadProjectedStage) -> ControlFlow<B>,
    {
        if self.next_iteration == self.kernel.settings.iterations {
            return Ok(ControlFlow::Continue(MultiLoadProjectedProgress::IterationLimit));
        }
        let count = self.kernel.load_cases.len().checked_mul(levels_per_candidate)
            .filter(|&count| count > 0)
            .ok_or_else(|| invalid("invalid complete mesh/load-family solve count"))?;
        if self.controls.max_solves - self.solves_started < count {
            return Ok(ControlFlow::Continue(MultiLoadProjectedProgress::SolveBudget(Vec::new())));
        }
        if let ControlFlow::Break(reason) = control(MultiLoadProjectedStage::Direction) {
            return Ok(ControlFlow::Break(reason));
        }
        let restoring = self.is_restoring_stress();
        let direction = if restoring {
            let case = self.current_stress.as_ref()
                .ok_or_else(|| invalid("restoration requires complete current stress"))?.worst_stress_case;
            self.kernel.direction_for_case(&self.current, self.next_iteration, Some(case))?
        } else {
            self.kernel.direction(&self.current, self.next_iteration)?
        };
        let mut attempts = Vec::with_capacity(self.controls.max_candidates);
        let mut move_cells = self.kernel.settings.move_cells;
        for index in 0..self.controls.max_candidates {
            if self.controls.max_solves - self.solves_started < count {
                return Ok(ControlFlow::Continue(MultiLoadProjectedProgress::SolveBudget(attempts)));
            }
            if let ControlFlow::Break(reason) = control(MultiLoadProjectedStage::Prepare(index)) {
                return Ok(ControlFlow::Break(reason));
            }
            let mut attempt = MultiLoadProjectedAttempt {
                index, move_cells, projection: None, state: None, stress: None, refusal: None,
            };
            let trial = self.kernel.propose_with_move(&self.current, &direction, self.ell, move_cells);
            move_cells *= self.controls.contraction;
            let Trial { mut phi, audit, events, load_pad_nodes } = match trial {
                Ok(trial) => trial,
                Err(error) => {
                    attempt.refusal = Some(format!("proposal: {error}"));
                    attempts.push(attempt);
                    continue;
                }
            };
            let projection = match project_material_volume_controlled(
                &mut phi, self.kernel.settings.level, &self.fixed, self.projection,
                |stage| control(MultiLoadProjectedStage::Projection(index, stage)),
            ) {
                Ok(ControlFlow::Continue(report)) => report,
                Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                Err(error) => {
                    attempt.refusal = Some(format!("projection: {error}"));
                    attempts.push(attempt);
                    continue;
                }
            };
            attempt.projection = Some(projection);
            let spent = &mut self.solves_started;
            let candidate = self.kernel.evaluate_scheduled(phi, poll_iters, |case, progress| {
                let stage = match progress {
                    CaseProgress::Start => MultiLoadProjectedStage::CaseSolve {
                        candidate: index, case, complete: false,
                    },
                    CaseProgress::Complete => MultiLoadProjectedStage::CaseSolve {
                        candidate: index, case, complete: true,
                    },
                    CaseProgress::Iterations(iterations) => MultiLoadProjectedStage::CaseIterations {
                        candidate: index, case, iterations,
                    },
                };
                if let ControlFlow::Break(reason) = control(stage) {
                    return ControlFlow::Break(reason);
                }
                // The complete family was admitted before any proposal work.
                // Charge BEFORE the solve, so numerical refusal spends budget.
                if matches!(progress, CaseProgress::Start) { *spent += 1; }
                ControlFlow::Continue(())
            });
            let candidate = match candidate {
                Ok(ControlFlow::Continue(candidate)) => candidate,
                Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                Err(error) => {
                    attempt.refusal = Some(format!("projected case family: {error}"));
                    attempts.push(attempt);
                    continue;
                }
            };
            let state = MultiLoadProjectedState::of(&candidate);
            attempt.state = Some(state.clone());
            let stress = if self.stress_limit.is_some() {
                match self.sample_stress_controlled(&candidate, |case, cell| {
                    control(MultiLoadProjectedStage::StressCell { candidate: index, case, cell })
                }) {
                    Ok(ControlFlow::Continue(evaluation)) => Some(evaluation),
                    Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                    Err(error) => {
                        attempt.refusal = Some(format!("projected stress family: {error}"));
                        attempts.push(attempt);
                        continue;
                    }
                }
            } else { None };
            attempt.stress = stress.clone();
            if (state.volume - self.projection.target).abs() > self.projection.tolerance {
                attempt.refusal = Some("projected area gate failed".into());
            } else if let Err(error) = self.require_candidate(&state, stress.as_ref()) {
                attempt.refusal = Some(error.to_string());
            } else {
                match admit(index, &candidate, stress.as_ref(), &mut self.solves_started, &mut control)? {
                    ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
                    ControlFlow::Continue(Some(reason)) => {
                        attempt.refusal = Some(reason);
                        attempts.push(attempt);
                        continue;
                    }
                    ControlFlow::Continue(None) => {}
                }
                let next_ell = self.ell + self.kernel.settings.mu_al
                    * direction.mean_energy.abs().max(1e-30)
                    * (candidate.volume - self.kernel.settings.volfrac) / self.kernel.settings.volfrac;
                if !next_ell.is_finite() {
                    return Err(invalid("projected multi-load search multiplier overflowed"));
                }
                if let ControlFlow::Break(reason) = control(MultiLoadProjectedStage::Publish(index)) {
                    return Ok(ControlFlow::Break(reason));
                }
                attempts.push(attempt);
                let step = MultiLoadProjectedStep {
                    iteration: self.next_iteration,
                    restoration: restoring,
                    previous: self.current(), state, stress: stress.clone(), projection,
                    proposal_audit: audit, proposal_events: events,
                    proposal_load_pad_nodes: load_pad_nodes, attempts,
                };
                self.current = candidate;
                self.current_stress = stress;
                self.ell = next_ell.max(0.0);
                self.next_iteration += 1;
                if restoring { self.restoration_updates += 1; }
                return Ok(ControlFlow::Continue(MultiLoadProjectedProgress::Accepted(Box::new(step))));
            }
            attempts.push(attempt);
        }
        Ok(ControlFlow::Continue(MultiLoadProjectedProgress::NoDescent(attempts)))
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod polling_tests;
