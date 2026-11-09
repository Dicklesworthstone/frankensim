//! Durable source, candidate and promotion boundaries for one native refit.
//! The invocation cap counts accepted updates across both optimizer owners.
use super::*;
use std::io::Write as _;

fn retained(old: &Loaded) -> Result<Outcome> {
    Ok(Outcome {
        pointer: format!("study-{}", old.hash.to_hex()),
        receipt: old.bytes.clone(),
        status: resume::status(old.value.str_field("status"))?,
    })
}

fn receipt_hash(out: &Outcome) -> Result<ContentHash> {
    out.pointer
        .strip_prefix("study-")
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| {
            fail(
                "cli-study-sdf3-refinement",
                "missing retained result identity",
            )
        })
}

#[allow(clippy::too_many_arguments)]
fn persist(
    spec: &Spec,
    ledger: &Ledger,
    study: &CutDensityStudy3<AdaptiveSolveSpace3>,
    state: &State,
    producer: ContentHash,
    previous: &mut Option<ContentHash>,
    limit: usize,
) -> Result<Outcome> {
    let out = output::persist(spec, ledger, study, state, producer, *previous, limit)?;
    *previous = Some(receipt_hash(&out)?);
    let phase = state
        .refinement
        .as_ref()
        .map_or("source", |c| c.phase.as_str());
    let _ = writeln!(
        std::io::stderr().lock(),
        "{{\"schema\":\"frankensim.cli.sdf3-affine-load-space-progress.v1\",\"run_id\":{},\"iterations_completed\":{},\"grid_phase\":{phase:?}}}",
        quoted(&out.pointer),
        state.optimizer.iterations
    );
    Ok(out)
}

fn allowance_finished(state: &State, spec: &Spec) -> bool {
    state.error.is_none()
        && match state.report.as_ref().map(|report| report.stop) {
            Some(ProjectedAlStop::Converged) => true,
            Some(ProjectedAlStop::IterationLimit) => state.optimizer.iterations == spec.updates,
            _ => false,
        }
}

fn recovery_allowance_finished(recovery: &resume::Recovery, spec: &Spec) -> bool {
    match recovery.optimizer_stop {
        Some(ProjectedAlStop::Converged) => true,
        Some(ProjectedAlStop::IterationLimit) => {
            recovery.checkpoint.optimizer.work.iterations == spec.updates
        }
        _ => false,
    }
}

fn source_wants_refit(run: &Computation, spec: &Spec) -> bool {
    let Some(config) = spec.affine_load_space else {
        return false;
    };
    allowance_finished(&run.state, spec)
        && config.feasible(run.state.selected())
        && spec.affine_grid_assessment.is_some_and(|policy| {
            run.state.assessment.as_ref().is_some_and(|assessment| {
                assessment.status == "checked" && !assessment.passed(policy, config.goal)
            })
        })
}

fn old_source_wants_refit(ledger: &Ledger, old: &Loaded) -> Result<bool> {
    let bytes = linked(ledger, &old.value, "report_json", "study-report-json")?;
    let report = JsonValue::parse(
        std::str::from_utf8(&bytes)
            .map_err(|e| fail("cli-study-sdf3-refinement", e.to_string()))?,
    )
    .map_err(|e| fail("cli-study-sdf3-refinement", e.to_string()))?;
    let assessment = report.get("grid_assessment").ok_or_else(|| {
        fail(
            "cli-study-sdf3-refinement",
            "missing retained grid assessment",
        )
    })?;
    match assessment.get("passed") {
        Some(JsonValue::Bool(passed)) => {
            Ok(!passed && report.get("selected_feasible") == Some(&JsonValue::Bool(true)))
        }
        _ => Err(fail(
            "cli-study-sdf3-refinement",
            "invalid retained assessment gate",
        )),
    }
}

fn failure_status(error: &Failure) -> &'static str {
    match error.exit {
        exit::CANCELLED => "cancelled",
        exit::BUDGET => "budget-exhausted",
        _ => "numerical-failure",
    }
}

fn reserve_restoration(recovery: &resume::Recovery, spec: &Spec) -> Result<()> {
    recovery.spent.require_remaining(spec)?;
    let maximum = spec
        .affine_load_space
        .ok_or_else(|| fail("cli-study-sdf3-refinement", "missing affine optimizer"))?
        .options
        .optimizer
        .max_evaluations;
    if recovery.checkpoint.restoration_cost()
        > maximum.saturating_sub(recovery.checkpoint.optimizer.work.evaluations)
    {
        return Err(Failure {
            code: "cli-study-sdf3-resume-budget",
            message: "the original working-grid evaluation allowance cannot rebuild the retained endpoints".into(),
            exit: exit::BUDGET,
        });
    }
    Ok(())
}

/// A successful candidate allowance is still provisional until the same
/// native owner freshly re-solves its complete selected physical family.
#[allow(clippy::too_many_arguments)]
fn finish_candidate(
    spec: &Spec,
    ledger: &Ledger,
    gate: &CancelGate,
    run: &mut Computation,
    producer: ContentHash,
    previous: &mut Option<ContentHash>,
    limit: usize,
) -> Result<Outcome> {
    let transition = Instant::now();
    let verified = allowance_finished(&run.state, spec)
        && spec.affine_load_space.is_some_and(|config| {
            run.state.refinement_recheck.as_ref().is_some_and(|result| {
                result.as_ref().is_ok_and(|point| {
                    point.rho == run.state.selected().rho && config.feasible(point)
                })
            })
        });
    if !verified {
        return persist(
            spec, ledger, &run.study, &run.state, producer, previous, limit,
        );
    }
    let context = run.state.refinement.as_mut().ok_or_else(|| {
        fail(
            "cli-study-sdf3-refinement",
            "candidate lost its source context",
        )
    })?;
    if context.phase != Phase::Candidate {
        return Err(fail(
            "cli-study-sdf3-refinement",
            "only a rechecked candidate can be promoted",
        ));
    }
    context.phase = Phase::Promoted;
    run.state.refinement_failure = None;
    run.state.assessment = None;
    run.state.status = "checkpointed";
    // Promotion and its exact fine-grid optimizer checkpoint survive a stop
    // during the subsequent, independently charged enriched-grid experiment.
    persist(
        spec, ledger, &run.study, &run.state, producer, previous, limit,
    )?;
    run.state.spent = run.state.spent.add(
        SolveWork::default(),
        QuadratureWork3::default(),
        transition.elapsed().as_secs_f64(),
    )?;
    if let Err(error) = assess(spec, gate, run) {
        run.state.status = failure_status(&error);
        run.state.refinement_failure = Some(error.message.clone());
        persist(
            spec, ledger, &run.study, &run.state, producer, previous, limit,
        )?;
        return Err(error);
    }
    persist(
        spec, ledger, &run.study, &run.state, producer, previous, limit,
    )
}

pub(in super::super) fn drive(
    original: &Spec,
    ledger: &Ledger,
    cap: Option<usize>,
    gate: &CancelGate,
    old: Option<&Loaded>,
) -> Result<Outcome> {
    let admission = Instant::now();
    let policy = original.affine_grid_refinement.ok_or_else(|| {
        fail(
            "cli-study-sdf3-refinement",
            "missing authored grid-refinement policy",
        )
    })?;
    policy.validate(original)?;
    let producer = checkpoint::producer_identity(gate)?;
    let recovery = old
        .map(|old| resume::recover(original, ledger, old, producer))
        .transpose()?;
    let effective = effective_spec(
        original,
        recovery.as_ref().and_then(|r| r.refinement.as_ref()),
    )?;
    let maximum = original
        .updates
        .checked_add(policy.max_updates)
        .ok_or_else(|| {
            fail(
                "cli-study-sdf3-refinement",
                "accepted-update allowance overflow",
            )
        })?;
    let limit = cap.unwrap_or(maximum).min(maximum);
    if let (Some(old), Some(recovery)) = (old, &recovery) {
        let phase = recovery.refinement.as_ref().map(|c| c.phase);
        let finished = recovery.optimizer_finished(&effective);
        let unchanged = match phase {
            Some(Phase::Promoted) => finished && recovery.assessment_checked,
            Some(Phase::Preparing) => limit == 0,
            Some(Phase::Candidate) => matches!(
                recovery.optimizer_stop,
                Some(ProjectedAlStop::Stalled | ProjectedAlStop::PenaltyLimit)
            ),
            None => {
                finished
                    && recovery.assessment_checked
                    && (!recovery_allowance_finished(recovery, &effective)
                        || !old_source_wants_refit(ledger, old)?)
            }
        };
        if unchanged {
            return retained(old);
        }
        reserve_restoration(recovery, &effective)?;
    }

    let mut previous = old.map(|old| old.hash);
    let result = (|| -> Result<Outcome> {
        let before = recovery
            .as_ref()
            .map_or(0, |r| r.checkpoint.optimizer.work.iterations);
        let mut run = compute_on_grid(
            &effective,
            gate,
            limit,
            recovery.as_ref(),
            admission.elapsed().as_secs_f64(),
            None,
            |study, state| {
                persist(
                    &effective,
                    ledger,
                    study,
                    state,
                    producer,
                    &mut previous,
                    limit,
                )?;
                Ok(())
            },
        )?;
        let transition = Instant::now();
        let used = run
            .state
            .optimizer
            .iterations
            .checked_sub(before)
            .ok_or_else(|| {
                fail(
                    "cli-study-sdf3-refinement",
                    "optimizer lost its accepted prefix",
                )
            })?;
        let remaining = limit.checked_sub(used).ok_or_else(|| {
            fail(
                "cli-study-sdf3-refinement",
                "invocation exceeded its accepted-update cap",
            )
        })?;
        match run.state.refinement.as_ref().map(|context| context.phase) {
            Some(Phase::Promoted) => {
                return persist(
                    &effective,
                    ledger,
                    &run.study,
                    &run.state,
                    producer,
                    &mut previous,
                    limit,
                );
            }
            Some(Phase::Candidate) => {
                return finish_candidate(
                    &effective,
                    ledger,
                    gate,
                    &mut run,
                    producer,
                    &mut previous,
                    limit,
                );
            }
            Some(Phase::Preparing) => {}
            None => {
                let source = persist(
                    &effective,
                    ledger,
                    &run.study,
                    &run.state,
                    producer,
                    &mut previous,
                    limit,
                )?;
                if !source_wants_refit(&run, &effective) {
                    return Ok(source);
                }
                run.state.refinement = Some(Context::new(&run, receipt_hash(&source)?)?);
            }
        }
        // Keep the source-grid AL checkpoint and exact marks before any fine
        // geometry/body-map work. A zero remaining step cap stops at this intent.
        run.state.status = "checkpointed";
        run.state.refinement_failure = None;
        let intent = persist(
            &effective,
            ledger,
            &run.study,
            &run.state,
            producer,
            &mut previous,
            limit,
        )?;
        if remaining == 0 {
            return Ok(intent);
        }
        let context = run
            .state
            .refinement
            .as_ref()
            .ok_or_else(|| {
                fail(
                    "cli-study-sdf3-refinement",
                    "missing prepared source context",
                )
            })?
            .clone();
        run.state.spent = run.state.spent.add(
            SolveWork::default(),
            QuadratureWork3::default(),
            transition.elapsed().as_secs_f64(),
        )?;
        let preparation = prepare(original, &run, context, gate);
        let prepared = match preparation.result {
            Ok(prepared) => prepared,
            Err(error) => {
                run.state.spent = preparation.spent;
                run.state.status = failure_status(&error);
                run.state.refinement_failure = Some(error.message.clone());
                persist(
                    &effective,
                    ledger,
                    &run.study,
                    &run.state,
                    producer,
                    &mut previous,
                    limit,
                )?;
                return Err(error);
            }
        };
        let candidate_admission = Instant::now();
        // A failed initial fine-grid audit has no admitted candidate state.
        // Its preceding COMPLETE seed work can still be durably charged to
        // the retained source without inventing a partial optimizer receipt.
        run.state.spent = prepared.spent;
        persist(
            &effective,
            ledger,
            &run.study,
            &run.state,
            producer,
            &mut previous,
            limit,
        )?;
        let candidate_spec = effective_spec(original, Some(&prepared.context))?;
        // The source has its own immutable artifact and the seed owns all fine
        // physics. Release the coarse numerical owner before candidate solves.
        drop(run);
        let mut candidate = compute_on_grid(
            &candidate_spec,
            gate,
            remaining,
            None,
            candidate_admission.elapsed().as_secs_f64(),
            Some(prepared),
            |study, state| {
                persist(
                    &candidate_spec,
                    ledger,
                    study,
                    state,
                    producer,
                    &mut previous,
                    limit,
                )?;
                Ok(())
            },
        )?;
        if candidate.state.optimizer.iterations > remaining {
            return Err(fail(
                "cli-study-sdf3-refinement",
                "candidate exceeded its invocation allowance",
            ));
        }
        finish_candidate(
            &candidate_spec,
            ledger,
            gate,
            &mut candidate,
            producer,
            &mut previous,
            limit,
        )
    })();
    result.map_err(|mut error| {
        if let Some(hash) = previous {
            error
                .message
                .push_str(&format!("; last durable result: study-{}", hash.to_hex()));
        }
        error
    })
}
