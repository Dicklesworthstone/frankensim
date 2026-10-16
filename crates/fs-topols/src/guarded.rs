//! Bounded whole-trajectory acceptance for the level-set elasticity optimizer.
//!
//! Each candidate starts from the identical admitted geometry, runs the existing
//! optimizer with a contracted interface-move limit, and is independently
//! re-solved by [`crate::optimize_compliance_evaluated`]. Selection therefore
//! uses the actual returned geometry's PDE objective and cut-quadrature area.
//! This is a deterministic bounded safeguard, not an Armijo proof, KKT test, or
//! global-optimality certificate.

use crate::{
    Cantilever, EvaluatedFinalState, EvaluatedOptimizeReport, GridSdf, OptimizeReport,
    OptimizeSettings, evaluate_compliance_design, optimize_compliance_evaluated,
};
use fs_cutfem::CutFemError;

/// Bounded candidate-family controls.
#[derive(Debug, Clone, Copy)]
pub struct GuardedSettings {
    /// Maximum complete optimization trajectories to evaluate.
    pub max_candidates: usize,
    /// Geometric contraction applied to `move_cells` after each candidate.
    pub contraction: f64,
    /// Numerical allowance above the requested material area fraction.
    pub volume_tolerance: f64,
    /// When the initial design is already volume-feasible, require at least this
    /// relative compliance reduction before a candidate can be published.
    pub min_relative_improvement: f64,
}

impl Default for GuardedSettings {
    fn default() -> Self {
        Self {
            max_candidates: 6,
            contraction: 0.5,
            volume_tolerance: 1e-3,
            min_relative_improvement: 0.0,
        }
    }
}

/// Outcome of one bounded candidate trajectory.
#[derive(Debug, Clone)]
pub struct GuardedCandidate {
    /// Zero-based deterministic candidate index.
    pub index: usize,
    /// Interface move limit supplied to this complete trajectory.
    pub move_cells: f64,
    /// Independently re-solved final design, when evolution and publication succeeded.
    pub final_state: Option<EvaluatedFinalState>,
    /// Whether the independently measured area met the declared upper limit.
    pub volume_feasible: bool,
    /// Whether the independently measured compliance met the improvement gate.
    pub improvement_gate: bool,
    /// Typed refusal rendered for retained diagnostics, if this candidate refused.
    pub refusal: Option<String>,
}

/// Why the bounded search stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GuardedStop {
    /// A feasible candidate satisfying the improvement rule was published.
    Accepted,
    /// Every successfully evaluated candidate violated the area limit.
    NoFeasibleCandidate,
    /// Feasible candidates existed but none met the requested compliance reduction.
    NoImprovingCandidate,
    /// Every candidate refused before an authoritative final state existed.
    AllCandidatesRefused,
}

/// Deterministic bounded-search evidence.
#[derive(Debug, Clone)]
pub struct GuardedOptimizeReport {
    /// Independently evaluated starting design.
    pub baseline: EvaluatedFinalState,
    /// Candidate attempts in deterministic move-size order.
    pub candidates: Vec<GuardedCandidate>,
    /// Search termination reason. Only `Accepted` mutates caller geometry.
    pub stop: GuardedStop,
    /// Accepted independently evaluated state, if any.
    pub accepted: Option<EvaluatedFinalState>,
    /// Accepted legacy trajectory evidence, if any. Its final row is not used for selection.
    pub trajectory: Option<OptimizeReport>,
    /// Move limit used by the accepted trajectory.
    pub accepted_move_cells: Option<f64>,
}

fn invalid_input(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn validate_controls(settings: OptimizeSettings, guarded: GuardedSettings) -> Result<(), CutFemError> {
    if settings.iterations == 0 {
        return Err(invalid_input("guarded optimization requires at least one evolution iteration"));
    }
    if !(settings.move_cells.is_finite() && settings.move_cells > 0.0) {
        return Err(invalid_input("guarded optimization requires finite positive move_cells"));
    }
    if !(settings.volfrac.is_finite() && settings.volfrac > 0.0 && settings.volfrac <= 1.0) {
        return Err(invalid_input("guarded optimization requires volfrac in (0, 1]"));
    }
    if guarded.max_candidates == 0 || guarded.max_candidates > 64 {
        return Err(invalid_input("guarded max_candidates must lie in [1, 64]"));
    }
    if !(guarded.contraction.is_finite() && guarded.contraction > 0.0 && guarded.contraction < 1.0) {
        return Err(invalid_input("guarded contraction must lie strictly between zero and one"));
    }
    if !(guarded.volume_tolerance.is_finite()
        && guarded.volume_tolerance >= 0.0
        && guarded.volume_tolerance <= 1.0)
    {
        return Err(invalid_input("guarded volume_tolerance must lie in [0, 1]"));
    }
    if !(guarded.min_relative_improvement.is_finite()
        && guarded.min_relative_improvement >= 0.0
        && guarded.min_relative_improvement < 1.0)
    {
        return Err(invalid_input("guarded min_relative_improvement must lie in [0, 1)"));
    }
    Ok(())
}

/// Run a bounded family of complete trajectories and publish the best candidate
/// selected from independently re-solved final states.
///
/// All candidates start from the same exact input geometry. Candidate `k+1` uses
/// `move_cells[k+1] = contraction * move_cells[k]`; every other optimization
/// setting, including nucleation schedule, remains unchanged. The caller's
/// geometry is unchanged unless a final candidate is accepted.
///
/// A candidate is area-feasible when `volume <= volfrac + volume_tolerance`.
/// If the baseline is already feasible, candidates must additionally reduce
/// compliance by `min_relative_improvement`. Among qualifying candidates, the
/// smallest independently evaluated compliance wins; ties retain the earlier,
/// larger-step candidate deterministically.
///
/// # Errors
/// Returns typed input/admission errors before candidate work. Individual
/// candidate evolution/solve refusals are retained in the report and do not
/// abort the bounded family.
pub fn optimize_compliance_guarded(
    phi: &mut GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
    guarded: GuardedSettings,
) -> Result<GuardedOptimizeReport, CutFemError> {
    validate_controls(settings, guarded)?;
    let baseline = evaluate_compliance_design(phi, fixture, settings)?;
    let limit = settings.volfrac + guarded.volume_tolerance;
    let baseline_feasible = baseline.volume <= limit;
    let improvement_limit = baseline.compliance * (1.0 - guarded.min_relative_improvement);
    let origin = phi.clone();
    let mut candidates = Vec::with_capacity(guarded.max_candidates);
    let mut best: Option<(GridSdf, EvaluatedOptimizeReport, f64)> = None;
    let mut move_cells = settings.move_cells;
    let mut any_final = false;
    let mut any_feasible = false;

    for index in 0..guarded.max_candidates {
        let mut trial_settings = settings;
        trial_settings.move_cells = move_cells;
        let mut trial = origin.clone();
        match optimize_compliance_evaluated(&mut trial, fixture, trial_settings) {
            Ok(report) => {
                any_final = true;
                let state = report.final_state;
                let volume_feasible = state.volume <= limit;
                any_feasible |= volume_feasible;
                let improvement_gate = !baseline_feasible || state.compliance <= improvement_limit;
                if volume_feasible && improvement_gate {
                    let replace = best.as_ref().is_none_or(|(_, current, _)| {
                        state.compliance < current.final_state.compliance
                    });
                    if replace {
                        best = Some((trial, report.clone(), move_cells));
                    }
                }
                candidates.push(GuardedCandidate {
                    index,
                    move_cells,
                    final_state: Some(state),
                    volume_feasible,
                    improvement_gate,
                    refusal: None,
                });
            }
            Err(error) => candidates.push(GuardedCandidate {
                index,
                move_cells,
                final_state: None,
                volume_feasible: false,
                improvement_gate: false,
                refusal: Some(format!("{error:?}")),
            }),
        }
        move_cells *= guarded.contraction;
    }

    if let Some((accepted_phi, report, accepted_move_cells)) = best {
        let accepted = report.final_state;
        *phi = accepted_phi;
        return Ok(GuardedOptimizeReport {
            baseline,
            candidates,
            stop: GuardedStop::Accepted,
            accepted: Some(accepted),
            trajectory: Some(report.trajectory),
            accepted_move_cells: Some(accepted_move_cells),
        });
    }

    let stop = if !any_final {
        GuardedStop::AllCandidatesRefused
    } else if !any_feasible {
        GuardedStop::NoFeasibleCandidate
    } else {
        GuardedStop::NoImprovingCandidate
    };
    Ok(GuardedOptimizeReport {
        baseline,
        candidates,
        stop,
        accepted: None,
        trajectory: None,
        accepted_move_cells: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn beam(level: u32) -> GridSdf {
        GridSdf::from_fn(1usize << level, &|_, y| (y - 0.5).abs() - 0.35)
    }

    #[test]
    fn invalid_guard_controls_refuse_without_mutation() {
        let mut phi = beam(3);
        let before = phi.nodes().to_vec();
        let settings = OptimizeSettings { level: 3, iterations: 1, ..OptimizeSettings::default() };
        let controls = GuardedSettings { max_candidates: 0, ..GuardedSettings::default() };
        assert!(optimize_compliance_guarded(
            &mut phi,
            Cantilever { load: 1.0, band: 0.125 },
            settings,
            controls,
        ).is_err());
        assert_eq!(phi.nodes(), before.as_slice());
    }

    #[test]
    fn zero_iteration_request_is_not_laundered_into_a_guarded_success() {
        let mut phi = beam(3);
        let before = phi.nodes().to_vec();
        let settings = OptimizeSettings { level: 3, iterations: 0, ..OptimizeSettings::default() };
        assert!(optimize_compliance_guarded(
            &mut phi,
            Cantilever { load: 1.0, band: 0.125 },
            settings,
            GuardedSettings::default(),
        ).is_err());
        assert_eq!(phi.nodes(), before.as_slice());
    }
}
