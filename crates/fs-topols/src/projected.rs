//! Feasible-baseline, hard-area-constrained elasticity descent on a level set.
//!
//! The existing checkpoint engine supplies each proposal (CutFEM elasticity,
//! Sobolev velocity, WENO evolution and its global nucleation schedule). This
//! owner projects the COMPLETE proposed field, re-solves that different field,
//! and accepts only a measured compliance decrease at the same numerical area.
//! Proposal rows/audits remain explicitly distinct from accepted-state evidence.
//! The proposal engine's AL multiplier is search state, NOT a KKT multiplier for
//! the projected geometry. No Armijo theorem, optimum, physical validation or
//! continuum-volume certificate is claimed.

use crate::checkpoint::{CheckpointStage, OptimizeCheckpoint};
use crate::volume::{
    VolumeProjectionReport, VolumeProjectionSettings, VolumeProjectionStage,
    project_material_volume, project_material_volume_controlled,
};
use crate::{
    Cantilever, EvaluatedFinalState, GridSdf, OptimizeReport, OptimizeSettings,
    evaluate_compliance_design,
};
use fs_cutfem::CutFemError;
use std::convert::Infallible;
use std::ops::ControlFlow;

/// Per-update search and cooperative-solver work controls.
#[derive(Debug, Clone, Copy)]
pub struct ProjectedSettings {
    /// Complete evolution/projection/evaluation attempts per accepted update.
    pub max_candidates: usize,
    /// Contraction of the proposal's interface travel after a rejected attempt.
    pub contraction: f64,
    /// Strict relative decrease required against the CURRENT feasible design.
    pub min_relative_improvement: f64,
    /// CG iterations between cancellation callbacks in the proposal solves.
    pub poll_iters: usize,
}

impl Default for ProjectedSettings {
    fn default() -> Self {
        Self { max_candidates: 6, contraction: 0.5, min_relative_improvement: 1e-8, poll_iters: 32 }
    }
}

/// Cooperative stages. The final re-evaluation is checked at its boundaries,
/// not internally preempted. One quadrature evaluation is also indivisible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectedStage {
    /// Before staging an attempt; candidate indices are zero-based.
    Prepare(usize),
    /// An existing checkpoint-engine boundary or CG poll.
    Proposal(usize, CheckpointStage),
    /// An area-projection boundary.
    Projection(usize, VolumeProjectionStage),
    /// Before independently solving the projected geometry.
    Evaluate(usize),
    /// After the independent solve, including candidates that will be rejected.
    Evaluated(usize),
    /// All gates passed, but durable accepted state is still unchanged.
    Publish(usize),
}

/// Retained outcome of a bounded candidate attempt.
#[derive(Debug, Clone)]
pub struct ProjectedAttempt {
    /// Candidate index within the current update.
    pub index: usize,
    /// Interface travel supplied to the proposal engine.
    pub move_cells: f64,
    /// Projection evidence when projection succeeded.
    pub projection: Option<VolumeProjectionReport>,
    /// Independent evaluation of the projected field, not the proposal field.
    pub state: Option<EvaluatedFinalState>,
    /// Numerical refusal or an explicit failed acceptance gate.
    pub refusal: Option<String>,
}

/// One accepted same-area update with separate proposal and publication data.
#[derive(Debug, Clone)]
pub struct ProjectedStep {
    /// Global zero-based ordinal of this accepted update.
    pub iteration: usize,
    /// Independently evaluated feasible state before the update.
    pub previous: EvaluatedFinalState,
    /// Independently evaluated exact geometry now held by the optimizer.
    pub state: EvaluatedFinalState,
    /// Area projection used by the accepted candidate.
    pub projection: VolumeProjectionReport,
    /// Unprojected proposal evidence; never the final geometry's certificate.
    pub proposal: OptimizeReport,
    /// Attempts through the first accepted candidate, in deterministic order.
    pub attempts: Vec<ProjectedAttempt>,
}

/// An iteration limit is not convergence; a search refusal is not an optimum.
#[derive(Debug, Clone)]
pub enum ProjectedProgress {
    /// A fully projected and independently evaluated candidate was committed.
    Accepted(Box<ProjectedStep>),
    /// The declared accepted-update count has already been reached.
    IterationLimit,
    /// All attempts failed or lacked sufficient decrease; state is unchanged.
    NoDescent(Vec<ProjectedAttempt>),
}

/// Stateful same-material descent, retaining only independently evaluated fields.
#[derive(Debug, Clone)]
pub struct ProjectedOptimizer {
    checkpoint: OptimizeCheckpoint,
    fixed: Vec<(usize, f64)>,
    projection: VolumeProjectionSettings,
    controls: ProjectedSettings,
    baseline: EvaluatedFinalState,
    current: EvaluatedFinalState,
}

fn refused(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn validate_controls(
    settings: OptimizeSettings,
    projection: VolumeProjectionSettings,
    controls: ProjectedSettings,
) -> Result<(), CutFemError> {
    if projection.target.to_bits() != settings.volfrac.to_bits() {
        return Err(refused("projected descent and volume projection must declare the same material target"));
    }
    if !(1..=8).contains(&settings.level)
        || settings.iterations == 0 || settings.iterations > 10_000
        || !(settings.move_cells.is_finite() && settings.move_cells > 0.0)
        || !(1..=64).contains(&controls.max_candidates)
        || !(controls.contraction.is_finite() && controls.contraction > 0.0 && controls.contraction < 1.0)
        || !(controls.min_relative_improvement.is_finite()
            && (0.0..1.0).contains(&controls.min_relative_improvement))
        || controls.poll_iters == 0
    {
        return Err(refused("projected descent requires level in 1..=8, 1..=10000 updates, positive move/poll controls, 1..=64 candidates, contraction in (0,1), and relative decrease in [0,1)"));
    }
    Ok(())
}

impl ProjectedOptimizer {
    /// Establish a FEASIBLE baseline before counting any objective improvement.
    ///
    /// Imposes the prescribed fixed nodes, projects area and independently solves
    /// the resulting geometry. Projection is feasibility restoration, not an
    /// optimization win against the possibly overfilled supplied design.
    ///
    /// # Errors
    /// Refuses invalid controls, unattainable projection or a failed baseline PDE.
    pub fn new(
        mut geometry: GridSdf,
        fixture: Cantilever,
        settings: OptimizeSettings,
        fixed: Vec<(usize, f64)>,
        projection: VolumeProjectionSettings,
        controls: ProjectedSettings,
    ) -> Result<Self, CutFemError> {
        validate_controls(settings, projection, controls)?;
        if geometry.n() != (1usize << settings.level) {
            return Err(refused("projected descent requires a level-matched input lattice"));
        }
        // The existing admission path validates ALL proposal controls before
        // running either area quadrature or the baseline PDE.
        let _ = OptimizeCheckpoint::new(geometry.clone(), fixture, settings)?;
        project_material_volume(&mut geometry, settings.level, &fixed, projection)?;
        let checkpoint = OptimizeCheckpoint::new(geometry, fixture, settings)?;
        Self::from_checkpoint(checkpoint, fixed, projection, controls)
    }

    /// Continue an already feasible retained state WITHOUT re-projecting it.
    ///
    /// Re-solves the retained geometry, checks fixed-node bits and numerical
    /// area, then resumes its global ordinal and proposal multiplier exactly.
    /// `baseline()` is the start of THIS segment, not a fabricated historical
    /// initial design. The caller retains prior segment evidence separately.
    ///
    /// # Errors
    /// Refuses mismatched controls, fixed values, area or canonical PDE failure.
    pub fn from_checkpoint(
        checkpoint: OptimizeCheckpoint,
        fixed: Vec<(usize, f64)>,
        projection: VolumeProjectionSettings,
        controls: ProjectedSettings,
    ) -> Result<Self, CutFemError> {
        validate_controls(checkpoint.settings(), projection, controls)?;
        // Validate projection admission on a disposable copy. A resumed state
        // must already be feasible; no changed field may be quietly substituted.
        let mut checked = checkpoint.geometry().clone();
        project_material_volume(&mut checked, checkpoint.settings().level, &fixed, projection)?;
        if checked.nodes().iter().zip(checkpoint.geometry().nodes())
            .any(|(left, right)| left.to_bits() != right.to_bits())
        {
            return Err(refused("resumed projected state does not already satisfy fixed nodes and material area"));
        }
        let current = evaluate_compliance_design(
            checkpoint.geometry(), checkpoint.fixture(), checkpoint.settings(),
        )?;
        if (current.volume - projection.target).abs() > projection.tolerance {
            return Err(refused("canonical baseline area disagrees with the projected material constraint"));
        }
        Ok(Self { checkpoint, fixed, projection, controls, baseline: current, current })
    }

    /// Exact feasible accepted geometry, global ordinal and proposal multiplier.
    #[must_use]
    pub fn checkpoint(&self) -> &OptimizeCheckpoint { &self.checkpoint }

    /// Independently solved same-material baseline for this segment.
    #[must_use]
    pub const fn baseline(&self) -> EvaluatedFinalState { self.baseline }

    /// Independently solved state of the exact retained geometry.
    #[must_use]
    pub const fn current(&self) -> EvaluatedFinalState { self.current }

    /// Attempt one accepted update with the declared bounded candidate family.
    ///
    /// # Errors
    /// Propagates an internal checkpoint admission failure. Numerical candidate
    /// refusals are retained in `NoDescent`, not promoted to convergence.
    pub fn advance_one(&mut self) -> Result<ProjectedProgress, CutFemError> {
        match self.advance_one_controlled(|_| ControlFlow::<Infallible>::Continue(()))? {
            ControlFlow::Continue(progress) => Ok(progress),
            ControlFlow::Break(never) => match never {},
        }
    }

    /// Controlled update; interruption at ANY boundary leaves accepted state
    /// unchanged, including interruption after a complete projected-field solve.
    /// Earlier accepted updates remain available through `checkpoint()`.
    ///
    /// # Errors
    /// Propagates an internal checkpoint admission failure.
    pub fn advance_one_controlled<B>(
        &mut self,
        control: impl FnMut(ProjectedStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, ProjectedProgress>, CutFemError> {
        self.advance_one_admitted(|_, _, _| Ok(None), control)
    }

    // The stress-constrained consumer shares the exact proposal/projection
    // search, rather than accepting an unconstrained step and checking it later.
    // An admission refusal consumes this candidate and contracts the next move;
    // the last accepted checkpoint and AL multiplier stay untouched.
    pub(crate) fn advance_one_admitted<B>(
        &mut self,
        mut admit: impl FnMut(usize, &GridSdf, EvaluatedFinalState) -> Result<Option<String>, CutFemError>,
        mut control: impl FnMut(ProjectedStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, ProjectedProgress>, CutFemError> {
        if self.checkpoint.is_complete() {
            return Ok(ControlFlow::Continue(ProjectedProgress::IterationLimit));
        }
        let settings = self.checkpoint.settings();
        let fixture = self.checkpoint.fixture();
        let ordinal = self.checkpoint.next_iteration();
        let mut attempts = Vec::with_capacity(self.controls.max_candidates);
        let mut move_cells = settings.move_cells;
        let limit = self.current.compliance * (1.0 - self.controls.min_relative_improvement);
        for index in 0..self.controls.max_candidates {
            if let ControlFlow::Break(reason) = control(ProjectedStage::Prepare(index)) {
                return Ok(ControlFlow::Break(reason));
            }
            let mut attempt = ProjectedAttempt {
                index, move_cells, projection: None, state: None, refusal: None,
            };
            let trial_settings = OptimizeSettings { move_cells, ..settings };
            let mut trial = OptimizeCheckpoint::restore(
                self.checkpoint.geometry().clone(), fixture, trial_settings,
                ordinal, self.checkpoint.ell(),
            )?;
            let proposal = match trial.advance_one_controlled(self.controls.poll_iters, |stage| {
                control(ProjectedStage::Proposal(index, stage))
            }) {
                Ok(ControlFlow::Continue(Some(report))) => report,
                Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                Ok(ControlFlow::Continue(None)) => return Err(refused("proposal unexpectedly reached the iteration limit")),
                Err(error) => {
                    attempt.refusal = Some(format!("proposal: {error}"));
                    attempts.push(attempt);
                    move_cells *= self.controls.contraction;
                    continue;
                }
            };
            let proposal_ell = trial.ell();
            let mut geometry = trial.into_geometry();
            let projection = match project_material_volume_controlled(
                &mut geometry, settings.level, &self.fixed, self.projection,
                |stage| control(ProjectedStage::Projection(index, stage)),
            ) {
                Ok(ControlFlow::Continue(report)) => report,
                Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                Err(error) => {
                    attempt.refusal = Some(format!("projection: {error}"));
                    attempts.push(attempt);
                    move_cells *= self.controls.contraction;
                    continue;
                }
            };
            attempt.projection = Some(projection);
            if let ControlFlow::Break(reason) = control(ProjectedStage::Evaluate(index)) {
                return Ok(ControlFlow::Break(reason));
            }
            let state = match evaluate_compliance_design(&geometry, fixture, settings) {
                Ok(state) => state,
                Err(error) => {
                    attempt.refusal = Some(format!("projected-state solve: {error}"));
                    attempts.push(attempt);
                    move_cells *= self.controls.contraction;
                    continue;
                }
            };
            if let ControlFlow::Break(reason) = control(ProjectedStage::Evaluated(index)) {
                return Ok(ControlFlow::Break(reason));
            }
            attempt.state = Some(state);
            if (state.volume - self.projection.target).abs() > self.projection.tolerance {
                attempt.refusal = Some("independent area gate failed".to_string());
            } else if !(state.compliance < limit) {
                attempt.refusal = Some("insufficient same-material compliance decrease".to_string());
            } else {
                attempt.refusal = match admit(index, &geometry, state) {
                    Ok(reason) => reason,
                    Err(error) => Some(format!("candidate constraint: {error}")),
                };
                if attempt.refusal.is_none() {
                    let accepted = OptimizeCheckpoint::restore(
                        geometry, fixture, settings, ordinal + 1, proposal_ell,
                    )?;
                    if let ControlFlow::Break(reason) = control(ProjectedStage::Publish(index)) {
                        return Ok(ControlFlow::Break(reason));
                    }
                    attempts.push(attempt);
                    let step = ProjectedStep {
                        iteration: ordinal, previous: self.current, state, projection, proposal, attempts,
                    };
                    self.checkpoint = accepted;
                    self.current = state;
                    return Ok(ControlFlow::Continue(ProjectedProgress::Accepted(Box::new(step))));
                }
            }
            attempts.push(attempt);
            move_cells *= self.controls.contraction;
        }
        Ok(ControlFlow::Continue(ProjectedProgress::NoDescent(attempts)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup() -> (GridSdf, OptimizeSettings, Vec<(usize, f64)>, VolumeProjectionSettings) {
        let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
        let settings = OptimizeSettings {
            level: 3, iterations: 2, volfrac: 0.6, move_cells: 0.1,
            nucleation_period: 0, ..OptimizeSettings::default()
        };
        let fixed = phi.nodes().iter().copied().enumerate()
            .filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect();
        let projection = VolumeProjectionSettings {
            target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
        };
        (phi, settings, fixed, projection)
    }

    fn optimizer() -> ProjectedOptimizer {
        let (phi, settings, fixed, projection) = setup();
        ProjectedOptimizer::new(phi, Cantilever { load: 1.0, band: 0.125 }, settings,
            fixed, projection, ProjectedSettings { max_candidates: 8, ..ProjectedSettings::default() })
            .expect("feasible baseline")
    }

    #[test]
    fn g0_baseline_is_feasible_and_matches_an_independent_endpoint_solve() {
        let optimizer = optimizer();
        assert!((optimizer.baseline().volume - 0.6).abs() <= 1e-4);
        assert_eq!(optimizer.checkpoint().next_iteration(), 0);
        let replay = evaluate_compliance_design(optimizer.checkpoint().geometry(),
            optimizer.checkpoint().fixture(), optimizer.checkpoint().settings()).unwrap();
        assert_eq!(replay.snapshot, optimizer.current().snapshot);
        assert_eq!(replay.compliance.to_bits(), optimizer.current().compliance.to_bits());
    }

    #[test]
    fn g4_cancelled_proposal_and_restored_state_retain_exact_feasible_geometry() {
        let mut optimizer = optimizer();
        let before = optimizer.checkpoint().geometry().nodes().to_vec();
        let result = optimizer.advance_one_controlled(|stage| {
            if matches!(stage, ProjectedStage::Proposal(_, CheckpointStage::Prepare)) {
                ControlFlow::Break("cancel")
            } else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break("cancel")));
        assert_eq!(optimizer.checkpoint().geometry().nodes(), before.as_slice());
        assert_eq!(optimizer.checkpoint().next_iteration(), 0);
        let (_, _, fixed, projection) = setup();
        let restored = ProjectedOptimizer::from_checkpoint(optimizer.checkpoint().clone(),
            fixed, projection, ProjectedSettings::default()).unwrap();
        assert_eq!(restored.current().snapshot, optimizer.current().snapshot);
        assert_eq!(restored.current().compliance.to_bits(), optimizer.current().compliance.to_bits());
    }

    #[test]
    fn g3_an_accepted_step_reduces_compliance_at_the_same_material_budget() {
        let mut optimizer = optimizer();
        let before = optimizer.current();
        let ProjectedProgress::Accepted(step) = optimizer.advance_one().expect("bounded update") else {
            panic!("the uniform-beam fixture must yield an informative accepted update");
        };
        assert!(step.state.compliance < before.compliance);
        assert!((step.state.volume - before.volume).abs() <= 2e-4);
        assert_eq!(step.state.snapshot, optimizer.current().snapshot);
        assert_eq!(optimizer.checkpoint().next_iteration(), 1);
        let replay = evaluate_compliance_design(optimizer.checkpoint().geometry(),
            optimizer.checkpoint().fixture(), optimizer.checkpoint().settings()).unwrap();
        assert_eq!(replay.snapshot, step.state.snapshot);
        assert_eq!(replay.compliance.to_bits(), step.state.compliance.to_bits());
        let mut continued = optimizer.clone();
        let (_, _, fixed, projection) = setup();
        let mut resumed = ProjectedOptimizer::from_checkpoint(optimizer.checkpoint().clone(),
            fixed, projection, optimizer.controls).unwrap();
        let a = continued.advance_one().unwrap();
        let b = resumed.advance_one().unwrap();
        assert_eq!(std::mem::discriminant(&a), std::mem::discriminant(&b));
        assert_eq!(continued.current().snapshot, resumed.current().snapshot);
        assert_eq!(continued.current().compliance.to_bits(), resumed.current().compliance.to_bits());
    }

    #[test]
    fn candidate_constraint_refusal_retries_without_publishing() {
        let mut optimizer = optimizer();
        let before = optimizer.clone();
        let mut screened = 0;
        let result = optimizer.advance_one_admitted(|_, _, state| {
            assert!(state.compliance < before.current().compliance);
            screened += 1;
            Ok(Some("constraint falsifier".to_string()))
        }, |_| ControlFlow::<Infallible>::Continue(())).unwrap();
        let ControlFlow::Continue(ProjectedProgress::NoDescent(attempts)) = result else {
            panic!("a refused candidate cannot be published");
        };
        assert!(screened > 0);
        assert_eq!(attempts.len(), optimizer.controls.max_candidates);
        assert!(attempts.iter().any(|a| a.refusal.as_deref() == Some("constraint falsifier")));
        assert!(attempts.windows(2).all(|a| a[1].move_cells < a[0].move_cells));
        assert_eq!(optimizer.checkpoint().geometry().nodes(), before.checkpoint().geometry().nodes());
        assert_eq!(optimizer.checkpoint().ell().to_bits(), before.checkpoint().ell().to_bits());
        assert_eq!(optimizer.checkpoint().next_iteration(), before.checkpoint().next_iteration());
        assert_eq!(optimizer.current().compliance.to_bits(), before.current().compliance.to_bits());
    }
}
