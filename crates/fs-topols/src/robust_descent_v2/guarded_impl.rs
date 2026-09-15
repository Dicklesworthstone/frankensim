use super::*;

fn validate_guarded(settings: OptimizeSettings, guarded: GuardedSettings) -> Result<(), CutFemError> {
    if settings.iterations == 0 {
        return Err(invalid("multi-load guarded optimization requires at least one evolution iteration"));
    }
    if guarded.max_candidates == 0 || guarded.max_candidates > 64 {
        return Err(invalid("multi-load guarded max_candidates must lie in [1, 64]"));
    }
    if !(guarded.contraction.is_finite() && guarded.contraction > 0.0 && guarded.contraction < 1.0) {
        return Err(invalid("multi-load guarded contraction must lie strictly between zero and one"));
    }
    if !(guarded.volume_tolerance.is_finite()
        && guarded.volume_tolerance >= 0.0
        && guarded.volume_tolerance <= 1.0)
    {
        return Err(invalid("multi-load guarded volume_tolerance must lie in [0, 1]"));
    }
    if !(guarded.min_relative_improvement.is_finite()
        && guarded.min_relative_improvement >= 0.0
        && guarded.min_relative_improvement < 1.0)
    {
        return Err(invalid("multi-load guarded min_relative_improvement must lie in [0, 1)"));
    }
    Ok(())
}

/// Bounded transactional candidate search whose candidate generation already
/// uses simultaneous independent-load descent. Every completed candidate is
/// independently replayed once more under all declared scenarios before it can
/// be published.
///
/// # Errors
/// Refuses malformed controls or baseline evaluation. Individual candidate
/// evolution/replay refusals are retained and never mutate caller geometry.
pub fn optimize_compliance_multi_load_guarded(
    phi: &mut GridSdf,
    load_cases: &[RobustLoadCase],
    settings: OptimizeSettings,
    guarded: GuardedSettings,
    aggregate: RobustAggregate,
) -> Result<RobustOptimizeReport, CutFemError> {
    validate_guarded(settings, guarded)?;
    let baseline = evaluate_robust_design(phi, load_cases, settings, aggregate)?;
    let limit = settings.volfrac + guarded.volume_tolerance;
    let baseline_feasible = baseline.volume <= limit;
    let improvement_limit = baseline.objective * (1.0 - guarded.min_relative_improvement);
    let origin = phi.clone();
    let mut candidates = Vec::with_capacity(guarded.max_candidates);
    let mut best: Option<(GridSdf, RobustEvaluation, OptimizeReport, f64)> = None;
    let mut move_cells = settings.move_cells;
    let mut any_final = false;
    let mut any_feasible = false;

    for index in 0..guarded.max_candidates {
        let mut trial_settings = settings;
        trial_settings.move_cells = move_cells;
        let mut trial = origin.clone();
        match optimize_compliance_multi_load(&mut trial, load_cases, trial_settings, aggregate) {
            Ok(descent) => match evaluate_robust_design(&trial, load_cases, settings, aggregate) {
                Ok(evaluation) => {
                    any_final = true;
                    let volume_feasible = evaluation.volume <= limit;
                    any_feasible |= volume_feasible;
                    let improvement_gate = !baseline_feasible || evaluation.objective <= improvement_limit;
                    if volume_feasible && improvement_gate {
                        let replace = best.as_ref().is_none_or(|(_, current, _, _)| {
                            evaluation.objective < current.objective
                        });
                        if replace {
                            best = Some((trial, evaluation.clone(), descent.trajectory.clone(), move_cells));
                        }
                    }
                    candidates.push(RobustCandidate {
                        index,
                        move_cells,
                        evaluation: Some(evaluation),
                        volume_feasible,
                        improvement_gate,
                        refusal: None,
                    });
                }
                Err(error) => candidates.push(RobustCandidate {
                    index,
                    move_cells,
                    evaluation: None,
                    volume_feasible: false,
                    improvement_gate: false,
                    refusal: Some(format!("{error:?}")),
                }),
            },
            Err(error) => candidates.push(RobustCandidate {
                index,
                move_cells,
                evaluation: None,
                volume_feasible: false,
                improvement_gate: false,
                refusal: Some(format!("{error:?}")),
            }),
        }
        move_cells *= guarded.contraction;
    }

    if let Some((accepted_phi, accepted, trajectory, accepted_move_cells)) = best {
        *phi = accepted_phi;
        return Ok(RobustOptimizeReport {
            baseline,
            candidates,
            stop: RobustStop::Accepted,
            accepted: Some(accepted),
            trajectory: Some(trajectory),
            accepted_move_cells: Some(accepted_move_cells),
        });
    }

    let stop = if !any_final {
        RobustStop::AllCandidatesRefused
    } else if !any_feasible {
        RobustStop::NoFeasibleCandidate
    } else {
        RobustStop::NoImprovingCandidate
    };
    Ok(RobustOptimizeReport {
        baseline,
        candidates,
        stop,
        accepted: None,
        trajectory: None,
        accepted_move_cells: None,
    })
}
