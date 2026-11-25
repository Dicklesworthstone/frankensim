//! Reusable durable multi-stage driver factored out of `crate::solve`
//! (bead `frankensim-rc-1109-root-3gc12.11.7`, Stage 1).
//!
//! Both the thermal cooling pipeline (`crate::solve`) and subsequent assembly
//! pipelines implement [`StageDescriptor`], [`StagePipeline`], and
//! [`ResumeAttestor`] so run-identity framing, session capability/budget
//! governance, cancellation gates, atomic per-stage ledger commits with
//! sealed checkpoints and edge seals, durable refusal recording, terminal
//! receipts, and controlled visible-op resume discovery live in one place.

#![allow(clippy::result_large_err)]

use std::collections::BTreeMap;
use std::fmt::Write as _;

use fs_blake3::{ContentHash, hash_domain};
use fs_ledger::{
    EdgeRole, ExecMode, FiveExplicits, Ledger, LedgerError, MAIN_BRANCH, OpOutcome, OpRow,
};
use fs_session::{CapabilityToken, Charge, Enforcement, Governor, SessionError, SessionId};

use crate::import::json_string;
use crate::solve::{
    CandidateOpReadError, CompletedStage, EDGE_SCAN_CAP, EvidenceReadError, EvidenceWork,
    InvocationWorkExceeded, JsonCursor, MAX_SOLVE_VISIBLE_OP_IDS, MAX_SOLVE_VISIBLE_OP_PAGES,
    SolveDriverState, SolveOutcome, SolveRefusal, SolveRunId, SolveRunStatus, StageSummary,
    VisibleOpScanError, cancelled_fresh_refusal_named, cancelled_resume_refusal_named,
    dump_trace_spans, invocation_work_refusal_named, parse_hash_string,
    read_artifact_info_controlled, read_candidate_op, read_op_edges_controlled,
    read_visible_op_page, resume_identity, resume_ledger,
};

#[cfg(test)]
mod tests;

/// Warn when consumed wall budget crosses these fractions of the grant.
pub(crate) const BUDGET_WARN_FRACTIONS: [f64; 2] = [0.5, 0.9];

/// Side artifact produced by a stage and retained as an `Out` edge beside the
/// stage receipt, with its precomputed content hash checked on persistence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct RetainedSideArtifact {
    pub(crate) kind: &'static str,
    pub(crate) artifact: ContentHash,
    pub(crate) bytes: Vec<u8>,
}

/// Ordered stage vocabulary for a durable CLI pipeline.
pub(crate) trait StageDescriptor: Copy + Eq + 'static {
    /// All stages in canonical execution order.
    fn all() -> &'static [Self];
    /// Zero-based index in [`Self::all`].
    fn ordinal(self) -> u32;
    /// Stable kebab-case stage name used in IR, receipts, progress, and diagnostics.
    fn name(self) -> &'static str;
    /// Capability verb the session token must grant for this stage.
    fn verb(self) -> &'static str;
    /// Owning bead id when a stage is reserved without an authoritative producer.
    fn gap_dependency(self) -> Option<&'static str>;
    /// Lookup by zero-based ordinal.
    fn from_ordinal(ordinal: u32) -> Option<Self> {
        Self::all().get(ordinal as usize).copied()
    }
}

/// Static wire schemas, artifact kinds, and diagnostics vocabulary for one pipeline.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PipelineSchema {
    pub(crate) driver_version: u32,
    pub(crate) stage_ir_schema: &'static str,
    pub(crate) run_receipt_schema: &'static str,
    pub(crate) progress_schema: &'static str,
    pub(crate) refusal_schema: &'static str,
    pub(crate) source_artifact_kind: &'static str,
    pub(crate) stage_state_kind: &'static str,
    pub(crate) stage_receipt_kind: &'static str,
    pub(crate) run_receipt_kind: &'static str,
    pub(crate) refusal_artifact_kind: &'static str,
    pub(crate) cores: u64,
    pub(crate) terminal_verb: &'static str,
    pub(crate) ledger_scope_prefix: &'static str,
    pub(crate) session_open_prefix: &'static str,
    pub(crate) meter_prefix: &'static str,
    pub(crate) resume_command_prefix: &'static str,
    pub(crate) budget_resume_hint: &'static str,
    pub(crate) incompatible_run_hint: &'static str,
    pub(crate) unknown_run_hint: &'static str,
    pub(crate) run_noun: &'static str,
    pub(crate) terminal_fallback_stage: &'static str,
}

/// Read-only view of the active session's rendered Five Explicits handed into
/// each stage body during execution.
pub(crate) struct StageSessionView<'a> {
    pub(crate) seed: &'a [u8; 8],
    pub(crate) versions_json: &'a str,
    pub(crate) budget_json: &'a str,
    pub(crate) capability_json: &'a str,
}

/// Output produced by one stage execution before atomic ledger commit.
pub(crate) struct StageExecution<C> {
    pub(crate) receipt_json: String,
    pub(crate) side_artifacts: Vec<RetainedSideArtifact>,
    pub(crate) on_commit: C,
}

/// Pipeline-specific hooks plugged into [`DurableStageSession::drive`].
pub(crate) trait StagePipeline {
    type Stage: StageDescriptor;
    type Context;
    type StageCommit;

    fn schema(&self) -> &'static PipelineSchema;
    fn canonical_source(&self) -> &[u8];
    fn seal_state(&self, run: SolveRunId, state: &SolveDriverState) -> Vec<u8>;
    fn declaration_gap(&self, stage: Self::Stage) -> Option<(&'static str, String, String)>;
    fn is_cancelled(&self) -> bool;
    fn pre_publish_checkpoint(&self, stage: Self::Stage) -> Result<(), ()>;
    fn execute_stage(
        &mut self,
        stage: Self::Stage,
        session: &StageSessionView<'_>,
        state: &SolveDriverState,
        context: &mut Self::Context,
    ) -> Result<StageExecution<Self::StageCommit>, SolveRefusal>;
    fn format_stage_progress(
        &self,
        _run_hex: &str,
        _stage: Self::Stage,
        _commit: &Self::StageCommit,
        _elapsed_s: f64,
    ) -> Result<Option<String>, SolveRefusal> {
        Ok(None)
    }
    fn commit_stage(&mut self, context: &mut Self::Context, commit: Self::StageCommit);
    fn link_initial_extra_inputs(&self, _ledger: &Ledger, _op: i64) -> Result<(), LedgerError> {
        Ok(())
    }
    fn link_stage_in_edges(
        &self,
        _ledger: &Ledger,
        _op: i64,
        _stage: Self::Stage,
        _state_before: &SolveDriverState,
        _context: &Self::Context,
    ) -> Result<(), LedgerError> {
        Ok(())
    }
}

/// Pipeline-specific candidate verification hooks plugged into
/// [`discover_and_attest_resume`].
pub(crate) trait ResumeAttestor {
    type Stage: StageDescriptor;
    type Verified;

    fn schema(&self) -> &'static PipelineSchema;
    fn decode_state(
        &self,
        ledger: &Ledger,
        run: SolveRunId,
        artifact: ContentHash,
        work: EvidenceWork<'_>,
    ) -> Result<SolveDriverState, SolveRefusal>;
    fn validate_candidate(
        &mut self,
        ledger: &Ledger,
        run: SolveRunId,
        state: SolveDriverState,
        state_artifact: ContentHash,
        discovery_op: i64,
        work: EvidenceWork<'_>,
    ) -> Result<Self::Verified, SolveRefusal>;
    fn completed_len(verified: &Self::Verified) -> usize;
}

/// Rendered Five Explicits carried by one durable session.
pub(crate) struct SessionExplicits {
    pub(crate) seed: [u8; 8],
    pub(crate) versions_json: String,
    pub(crate) budget_json: String,
    pub(crate) capability_json: String,
}

/// Reusable governed durable execution session shared across CLI pipelines.
pub(crate) struct DurableStageSession<'a> {
    ledger: &'a Ledger,
    clock: &'a mut dyn FnMut() -> f64,
    run: SolveRunId,
    project_hash: ContentHash,
    schema: &'static PipelineSchema,
    governor: Governor,
    session: SessionId,
    token: CapabilityToken,
    versions_json: String,
    budget_json: String,
    capability_json: String,
    seed: [u8; 8],
    progress: Vec<String>,
}

impl<'a> DurableStageSession<'a> {
    pub(crate) fn open(
        ledger: &'a Ledger,
        clock: &'a mut dyn FnMut() -> f64,
        run: SolveRunId,
        project_hash: ContentHash,
        schema: &'static PipelineSchema,
        token: CapabilityToken,
        explicits: SessionExplicits,
    ) -> Result<Self, SolveRefusal> {
        if ledger.in_transaction() {
            return Err(SolveRefusal::plain(
                "cli-solve-ledger-transaction",
                "the ledger connection is already inside a caller-owned transaction",
                "commit or roll back before solve so stage groups stay atomic",
            ));
        }
        let governor = Governor::new();
        let session = token.session;
        let open_id = governor
            .session_open_id(
                session,
                &format!("{}-{}", schema.session_open_prefix, run.to_hex()),
            )
            .map_err(|error| session_refusal(run, &error))?;
        governor
            .open_session_declared(open_id, token.clone())
            .map_err(|error| session_refusal(run, &error))?;
        Ok(Self {
            ledger,
            clock,
            run,
            project_hash,
            schema,
            governor,
            session,
            token,
            versions_json: explicits.versions_json,
            budget_json: explicits.budget_json,
            capability_json: explicits.capability_json,
            seed: explicits.seed,
            progress: Vec::new(),
        })
    }

    /// Re-charge prior consumption on resume so the run continues under its
    /// original budget grant instead of resetting it.
    pub(crate) fn restore_budget(
        &self,
        state: &SolveDriverState,
        state_artifact: ContentHash,
    ) -> Result<(), SolveRefusal> {
        if state.consumed_core_s > 0.0 || state.consumed_wall_s > 0.0 {
            let enforcement = self.charge(
                "resume-restore",
                Charge {
                    core_s: state.consumed_core_s,
                    mem_peak_bytes: 0,
                    wall_s: state.consumed_wall_s,
                },
            )?;
            if !matches!(enforcement, Enforcement::Ok) {
                return Err(SolveRefusal {
                    code: "cli-solve-resume-budget",
                    stage: None,
                    what: format!(
                        "the recorded consumption already exhausts the project budget \
                         (state artifact {})",
                        state_artifact.to_hex()
                    ),
                    fix: "raise the project budgets; the changed project starts a fresh run whose \
                          completed artifacts deduplicate by content"
                        .to_string(),
                    dependency: None,
                    run: Some(self.run.to_hex()),
                    recorded_op: None,
                });
            }
        }
        Ok(())
    }

    pub(crate) fn drain_progress(&mut self, sink: &mut Vec<String>) {
        sink.append(&mut self.progress);
    }

    #[allow(clippy::too_many_lines)]
    pub(crate) fn drive<P: StagePipeline>(
        &mut self,
        pipeline: &mut P,
        mut state: SolveDriverState,
        mut context: P::Context,
        mut predecessor_state: Option<ContentHash>,
    ) -> Result<SolveOutcome, SolveRefusal> {
        debug_assert_eq!(pipeline.schema().driver_version, self.schema.driver_version);
        let prior_stages = u32::try_from(state.completed.len()).expect("stage count fits u32");
        let mut summaries = Vec::new();
        let mut warned = [false; BUDGET_WARN_FRACTIONS.len()];
        for (index, &fraction) in BUDGET_WARN_FRACTIONS.iter().enumerate() {
            warned[index] = state.consumed_wall_s >= self.token.wall_s * fraction;
        }
        let run_hex = self.run.to_hex();
        for &stage in P::Stage::all().iter().skip(state.completed.len()) {
            if pipeline.is_cancelled() {
                return Err(self.cancelled_refusal::<P::Stage>(&state));
            }
            if let Some(dependency) = stage.gap_dependency() {
                let refusal = SolveRefusal {
                    code: "cli-solve-stage-gap",
                    stage: Some(stage.name()),
                    what: format!(
                        "stage `{}` is reserved but cannot execute until `{dependency}` \
                         supplies its authoritative producer",
                        stage.name()
                    ),
                    fix: format!(
                        "complete and verify `{dependency}`; do not substitute a skeleton \
                         stage or placeholder artifact"
                    ),
                    dependency: Some(dependency),
                    run: Some(run_hex.clone()),
                    recorded_op: None,
                };
                return Err(self.record_refusal(pipeline, &state, stage, refusal));
            }
            if let Some((code, what, fix)) = pipeline.declaration_gap(stage) {
                let refusal = SolveRefusal {
                    code,
                    stage: Some(stage.name()),
                    what,
                    fix,
                    dependency: None,
                    run: Some(run_hex.clone()),
                    recorded_op: None,
                };
                return Err(self.record_refusal(pipeline, &state, stage, refusal));
            }
            if !self.token.grants_op(stage.verb()) {
                let refusal = staged_refusal(
                    "cli-solve-capability",
                    stage.name(),
                    format!("the session token does not grant `{}`", stage.verb()),
                    "derive the token from project budgets through the solve driver",
                );
                return Err(self.record_refusal(pipeline, &state, stage, refusal));
            }
            let started = (self.clock)();
            if !started.is_finite() {
                let refusal = staged_refusal(
                    "cli-solve-budget",
                    stage.name(),
                    "the solve clock returned a non-finite stage start",
                    "supply a finite monotonic clock; no checkpoint was sealed",
                );
                return Err(self.record_refusal(pipeline, &state, stage, refusal));
            }
            let view = StageSessionView {
                seed: &self.seed,
                versions_json: &self.versions_json,
                budget_json: &self.budget_json,
                capability_json: &self.capability_json,
            };
            let execution = match pipeline.execute_stage(stage, &view, &state, &mut context) {
                Ok(produced) => produced,
                Err(refusal) if refusal.code == "cli-solve-cancelled" => {
                    return Err(self.cancelled_refusal::<P::Stage>(&state));
                }
                Err(refusal) if refusal.code == "cli-solve-work-envelope" => return Err(refusal),
                Err(refusal) => return Err(self.record_refusal(pipeline, &state, stage, refusal)),
            };
            if pipeline.pre_publish_checkpoint(stage).is_err() {
                return Err(self.cancelled_refusal::<P::Stage>(&state));
            }
            let finished = (self.clock)();
            let elapsed = finished - started;
            if !finished.is_finite() || !elapsed.is_finite() || elapsed < 0.0 {
                let refusal = staged_refusal(
                    "cli-solve-budget",
                    stage.name(),
                    format!(
                        "the solve clock did not produce a finite monotonic interval \
                         (start={started}, finish={finished})"
                    ),
                    "supply a finite monotonic clock; no checkpoint was sealed",
                );
                return Err(self.record_refusal(pipeline, &state, stage, refusal));
            }
            let custom_progress = match pipeline.format_stage_progress(
                &run_hex,
                stage,
                &execution.on_commit,
                elapsed,
            ) {
                Ok(line) => line,
                Err(refusal) => return Err(self.record_refusal(pipeline, &state, stage, refusal)),
            };
            #[allow(clippy::cast_precision_loss)]
            let charge = Charge {
                core_s: elapsed * self.schema.cores as f64,
                mem_peak_bytes: 0,
                wall_s: elapsed,
            };
            let consumed_core_s = state.consumed_core_s + charge.core_s;
            let consumed_wall_s = state.consumed_wall_s + charge.wall_s;
            if !consumed_core_s.is_finite()
                || consumed_core_s < 0.0
                || !consumed_wall_s.is_finite()
                || consumed_wall_s < 0.0
            {
                let refusal = staged_refusal(
                    "cli-solve-budget",
                    stage.name(),
                    "the accumulated solve consumption is non-finite or negative",
                    "repair the clock/budget source; no checkpoint was sealed",
                );
                return Err(self.record_refusal(pipeline, &state, stage, refusal));
            }
            state.consumed_core_s = consumed_core_s;
            state.consumed_wall_s = consumed_wall_s;
            let (op_id, receipt_hash, state_hash) = self
                .persist_stage(
                    pipeline,
                    &state,
                    stage,
                    &execution.receipt_json,
                    &execution.side_artifacts,
                    &context,
                    predecessor_state,
                )
                .map_err(|error| ledger_refusal(self.run, stage.name(), &error))?;
            pipeline.commit_stage(&mut context, execution.on_commit);
            state.completed.push(CompletedStage {
                ordinal: stage.ordinal(),
                op_id,
                receipt: receipt_hash,
            });
            predecessor_state = Some(state_hash);
            summaries.push(StageSummary {
                stage: stage.name(),
                op_id,
                receipt: receipt_hash.to_hex(),
                wall_s: elapsed,
            });
            self.progress.push(custom_progress.unwrap_or_else(|| {
                format_progress_line(
                    self.schema.progress_schema,
                    &run_hex,
                    stage.name(),
                    stage.ordinal(),
                    "ok",
                    elapsed,
                )
            }));
            let enforcement = self.charge(stage.verb(), charge)?;
            match enforcement {
                Enforcement::Ok => {
                    for (index, &fraction) in BUDGET_WARN_FRACTIONS.iter().enumerate() {
                        if !warned[index] && state.consumed_wall_s >= self.token.wall_s * fraction {
                            warned[index] = true;
                            self.progress.push(format_budget_warning_line(
                                self.schema.progress_schema,
                                &run_hex,
                                fraction,
                                state.consumed_wall_s,
                                self.token.wall_s,
                            ));
                        }
                    }
                }
                Enforcement::Throttled {
                    resource,
                    used,
                    granted,
                }
                | Enforcement::Paused {
                    resource,
                    used,
                    granted,
                    ..
                } => {
                    let status = SolveRunStatus::BudgetExceeded {
                        resource,
                        used,
                        granted,
                    };
                    let receipt = self
                        .persist_terminal::<P::Stage>(&state, &status)
                        .map_err(|error| ledger_refusal(self.run, stage.name(), &error))?;
                    return Ok(SolveOutcome {
                        run: run_hex,
                        status,
                        stages: summaries,
                        prior_stages,
                        run_receipt: Some(receipt.to_hex()),
                    });
                }
            }
        }
        let status = SolveRunStatus::Completed;
        let receipt = self
            .persist_terminal::<P::Stage>(&state, &status)
            .map_err(|error| {
                ledger_refusal(self.run, self.schema.terminal_fallback_stage, &error)
            })?;
        Ok(SolveOutcome {
            run: run_hex,
            status,
            stages: summaries,
            prior_stages,
            run_receipt: Some(receipt.to_hex()),
        })
    }

    fn charge(&self, meter_key: &str, charge: Charge) -> Result<Enforcement, SolveRefusal> {
        let report_id = self
            .governor
            .meter_report_id(
                self.session,
                &format!(
                    "{}-{}-{meter_key}",
                    self.schema.meter_prefix,
                    self.run.to_hex()
                ),
            )
            .map_err(|error| session_refusal(self.run, &error))?;
        let receipt = self
            .governor
            .charge(report_id, charge)
            .map_err(|error| session_refusal(self.run, &error))?;
        Ok(receipt.enforcement().clone())
    }

    fn cancelled_refusal<S: StageDescriptor>(&self, state: &SolveDriverState) -> SolveRefusal {
        if state.completed.is_empty() {
            return cancelled_fresh_refusal_named(self.run, S::from_ordinal(0).map(S::name));
        }
        SolveRefusal {
            code: "cli-solve-cancelled",
            stage: S::from_ordinal(
                u32::try_from(state.completed.len()).expect("stage count fits u32"),
            )
            .map(S::name),
            what: format!(
                "cancellation observed after {} completed stage(s); the completed prefix \
                 is durable",
                state.completed.len()
            ),
            fix: format!(
                "resume with `{} {} <ledger>`",
                self.schema.resume_command_prefix,
                self.run.to_hex()
            ),
            dependency: None,
            run: Some(self.run.to_hex()),
            recorded_op: None,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn persist_stage<P: StagePipeline>(
        &self,
        pipeline: &P,
        state_before: &SolveDriverState,
        stage: P::Stage,
        receipt_json: &str,
        side_artifacts: &[RetainedSideArtifact],
        context: &P::Context,
        predecessor_state: Option<ContentHash>,
    ) -> Result<(i64, ContentHash, ContentHash), LedgerError> {
        let ordinal = stage.ordinal();
        if state_before.completed.len() != ordinal as usize {
            return Err(LedgerError::Invalid {
                field: "solve_stage_prefix".to_string(),
                problem: format!(
                    "stage `{}` ordinal {ordinal} received a {}-stage prefix",
                    stage.name(),
                    state_before.completed.len()
                ),
            });
        }
        if predecessor_state.is_some() != (ordinal > 0) {
            return Err(LedgerError::Invalid {
                field: "solve_stage_predecessor".to_string(),
                problem: format!(
                    "stage `{}` ordinal {ordinal} has the wrong predecessor-state presence",
                    stage.name()
                ),
            });
        }
        let ir = format_stage_ir(self.schema, stage, self.run, self.project_hash);
        self.ledger.begin()?;
        let result = (|| -> Result<(i64, ContentHash, ContentHash), LedgerError> {
            let op = self.ledger.begin_op(
                Some(self.run.as_bytes()),
                &ir,
                &FiveExplicits {
                    seed: &self.seed,
                    versions: &self.versions_json,
                    budget: &self.budget_json,
                    capability: &self.capability_json,
                },
                i64::from(ordinal) * 2,
            )?;
            if let Some(predecessor) = predecessor_state {
                self.ledger.link(op, &predecessor, EdgeRole::In)?;
            }
            if ordinal == 0 {
                let source = self.ledger.put_artifact(
                    self.schema.source_artifact_kind,
                    pipeline.canonical_source(),
                    None,
                )?;
                self.ledger.link(op, &source.hash, EdgeRole::In)?;
                pipeline.link_initial_extra_inputs(self.ledger, op)?;
            }
            pipeline.link_stage_in_edges(self.ledger, op, stage, state_before, context)?;
            for usage in side_artifacts {
                let retained = self.ledger.put_artifact(usage.kind, &usage.bytes, None)?;
                if retained.hash != usage.artifact {
                    return Err(LedgerError::Invalid {
                        field: "solve_stage_side_artifact".to_string(),
                        problem: format!(
                            "stage `{}` precomputed artifact {} but ledger retained {}",
                            stage.name(),
                            usage.artifact.to_hex(),
                            retained.hash.to_hex()
                        ),
                    });
                }
                self.ledger.link(op, &retained.hash, EdgeRole::Out)?;
            }
            let receipt = self.ledger.put_artifact(
                self.schema.stage_receipt_kind,
                receipt_json.as_bytes(),
                None,
            )?;
            self.ledger.link(op, &receipt.hash, EdgeRole::Out)?;
            let mut sealed_state = state_before.clone();
            sealed_state.completed.push(CompletedStage {
                ordinal,
                op_id: op,
                receipt: receipt.hash,
            });
            let envelope = pipeline.seal_state(self.run, &sealed_state);
            let state_artifact =
                self.ledger
                    .put_artifact(self.schema.stage_state_kind, &envelope, None)?;
            self.ledger.link(op, &state_artifact.hash, EdgeRole::Out)?;
            let edges = self.ledger.op_artifact_edges_bounded(op, EDGE_SCAN_CAP)?;
            if edges.truncated {
                return Err(LedgerError::Invalid {
                    field: "solve_stage_edges".to_string(),
                    problem: format!(
                        "stage `{}` exceeds the {EDGE_SCAN_CAP}-edge verification cap",
                        stage.name()
                    ),
                });
            }
            self.ledger.seal_op_artifact_edges(op, edges.edges.len())?;
            self.ledger
                .finish_op(op, OpOutcome::Ok, None, i64::from(ordinal) * 2 + 1)?;
            Ok((op, receipt.hash, state_artifact.hash))
        })();
        finish_ledger_transaction(self.ledger, result)
    }

    fn persist_terminal<S: StageDescriptor>(
        &self,
        state: &SolveDriverState,
        status: &SolveRunStatus,
    ) -> Result<ContentHash, LedgerError> {
        let (status_name, detail) = match status {
            SolveRunStatus::Completed => ("completed", String::new()),
            SolveRunStatus::BudgetExceeded {
                resource,
                used,
                granted,
            } => (
                "budget-exceeded",
                format!(
                    ",\"resource\":{},\"used\":{used},\"granted\":{granted},\"resume\":{}",
                    json_string(resource),
                    json_string(self.schema.budget_resume_hint),
                ),
            ),
            SolveRunStatus::Cancelled => ("cancelled", String::new()),
        };
        let receipt_json = format_run_receipt_json::<S>(
            self.schema.run_receipt_schema,
            self.run,
            self.project_hash,
            status_name,
            &detail,
            state,
        );
        let ir = format_terminal_ir::<S>(self.schema, self.run, self.project_hash);
        self.ledger.begin()?;
        let result = (|| -> Result<ContentHash, LedgerError> {
            let op = self.ledger.begin_op(
                Some(self.run.as_bytes()),
                &ir,
                &FiveExplicits {
                    seed: &self.seed,
                    versions: &self.versions_json,
                    budget: &self.budget_json,
                    capability: &self.capability_json,
                },
                100,
            )?;
            let receipt = self.ledger.put_artifact(
                self.schema.run_receipt_kind,
                receipt_json.as_bytes(),
                None,
            )?;
            self.ledger.link(op, &receipt.hash, EdgeRole::Out)?;
            self.ledger.finish_op(op, OpOutcome::Ok, None, 101)?;
            Ok(receipt.hash)
        })();
        finish_ledger_transaction(self.ledger, result)
    }

    fn record_refusal<P: StagePipeline>(
        &self,
        pipeline: &P,
        state: &SolveDriverState,
        stage: P::Stage,
        mut refusal: SolveRefusal,
    ) -> SolveRefusal {
        refusal.run = Some(self.run.to_hex());
        let diagnostic = format!(
            "{{\"schema\":{},\"run\":{},\"stage\":{},\"code\":{},\"what\":{},\"fix\":{}{}}}",
            json_string(self.schema.refusal_schema),
            json_string(&self.run.to_hex()),
            json_string(stage.name()),
            json_string(refusal.code),
            json_string(&refusal.what),
            json_string(&refusal.fix),
            refusal.dependency.map_or(String::new(), |dependency| {
                format!(",\"dependency\":{}", json_string(dependency))
            }),
        );
        let ir = format_stage_ir(self.schema, stage, self.run, self.project_hash);
        self.ledger.begin().ok();
        let recorded = (|| -> Result<i64, LedgerError> {
            let op = self.ledger.begin_op(
                Some(self.run.as_bytes()),
                &ir,
                &FiveExplicits {
                    seed: &self.seed,
                    versions: &self.versions_json,
                    budget: &self.budget_json,
                    capability: &self.capability_json,
                },
                i64::from(stage.ordinal()) * 2,
            )?;
            let artifact = self.ledger.put_artifact(
                self.schema.refusal_artifact_kind,
                diagnostic.as_bytes(),
                None,
            )?;
            self.ledger.link(op, &artifact.hash, EdgeRole::Out)?;
            if stage.ordinal() == 0 && state.completed.is_empty() {
                let source = self.ledger.put_artifact(
                    self.schema.source_artifact_kind,
                    pipeline.canonical_source(),
                    None,
                )?;
                self.ledger.link(op, &source.hash, EdgeRole::In)?;
            }
            self.ledger.finish_op(
                op,
                OpOutcome::Error,
                Some(&diagnostic),
                i64::from(stage.ordinal()) * 2 + 1,
            )?;
            Ok(op)
        })();
        match finish_ledger_transaction(self.ledger, recorded) {
            Ok(op) => refusal.recorded_op = Some(op),
            Err(ledger_error) => {
                refusal.what = format!(
                    "{}; durable refusal recording also failed: {ledger_error}",
                    refusal.what
                );
            }
        }
        refusal
    }
}

/// Discover and independently attest the longest valid durable checkpoint for
/// `run`, refusing competing checkpoints of equal prefix length or runs
/// retained under an incompatible driver version.
#[allow(clippy::too_many_lines)]
pub(crate) fn discover_and_attest_resume<A: ResumeAttestor>(
    ledger: &Ledger,
    run: SolveRunId,
    work: EvidenceWork<'_>,
    attestor: &mut A,
) -> Result<A::Verified, SolveRefusal> {
    let schema = attestor.schema();
    let first = A::Stage::from_ordinal(0).map(StageDescriptor::name);
    let cancel_err = || cancelled_resume_refusal_named(run, first, schema);
    let work_err = |e| invocation_work_refusal_named(Some(run), first, e);
    let mut cursor = None;
    let mut page_index = 0usize;
    let mut candidate_index = 0usize;
    let mut visible_ids = 0usize;
    let mut best: Option<A::Verified> = None;
    let mut best_is_ambiguous = false;
    let mut candidates: Vec<(SolveDriverState, ContentHash, i64)> = Vec::new();
    let trace_export = std::env::var_os("FS_CLI_TRACE_EXPORT").is_some();
    let t_scan = std::time::Instant::now();
    let run_bytes: &[u8] = run.as_bytes().as_slice();
    let mut seen_for_run = 0usize;
    let mut rejections: BTreeMap<String, usize> = BTreeMap::new();
    loop {
        if page_index >= MAX_SOLVE_VISIBLE_OP_PAGES {
            return Err(work_err(InvocationWorkExceeded::VisiblePages {
                limit: MAX_SOLVE_VISIBLE_OP_PAGES,
            }));
        }
        let page =
            read_visible_op_page(ledger, cursor, work, page_index).map_err(
                |error| match error {
                    VisibleOpScanError::Cancelled => cancel_err(),
                    VisibleOpScanError::WorkEnvelope(error) => work_err(error),
                    VisibleOpScanError::Ledger(error) => {
                        resume_ledger("scanning the ledger for the run failed", error)
                    }
                },
            )?;
        page_index = page_index.saturating_add(1);
        if page.truncated != page.next_cursor.is_some() {
            return Err(resume_ledger(
                "scanning the ledger for the run failed",
                LedgerError::Invalid {
                    field: "visible_op_page".to_string(),
                    problem: "controlled page truncation and continuation disagree".to_string(),
                },
            ));
        }
        for id in page.op_ids {
            if visible_ids >= MAX_SOLVE_VISIBLE_OP_IDS {
                return Err(work_err(InvocationWorkExceeded::VisibleIds {
                    limit: MAX_SOLVE_VISIBLE_OP_IDS,
                }));
            }
            visible_ids = visible_ids.saturating_add(1);
            let this_candidate = candidate_index;
            candidate_index = candidate_index.saturating_add(1);
            let Some(candidate) = read_candidate_op(ledger, id, work, this_candidate).map_err(
                |error| match error {
                    CandidateOpReadError::Cancelled => cancel_err(),
                    CandidateOpReadError::WorkEnvelope(error) => work_err(error),
                    CandidateOpReadError::Ledger(error) => {
                        resume_ledger("scanning the ledger for the run failed", error)
                    }
                },
            )?
            else {
                continue;
            };
            let this_run = candidate.row.session.as_deref() == Some(run_bytes);
            if this_run {
                seen_for_run = seen_for_run.saturating_add(1);
            }
            let rejection = if candidate.branch != MAIN_BRANCH {
                Some("op is not on the main branch".to_string())
            } else if candidate.exec_mode != ExecMode::Deterministic {
                Some("op did not execute deterministically".to_string())
            } else {
                stage_discovery_rejection_for::<A::Stage>(&candidate.row, run, schema).err()
            };
            if let Some(reason) = rejection {
                if this_run {
                    *rejections.entry(reason).or_insert(0) += 1;
                }
                continue;
            }
            let mut saw_state = false;
            let edges =
                read_op_edges_controlled(ledger, id, work, this_candidate).map_err(|error| {
                    match error {
                        EvidenceReadError::Cancelled => cancel_err(),
                        EvidenceReadError::WorkEnvelope(error) => work_err(error),
                        EvidenceReadError::Ledger(error) => {
                            resume_ledger("reading solve operation edges failed", error)
                        }
                    }
                })?;
            if edges.truncated {
                return Err(resume_identity(format!(
                    "operation {id} exceeds the complete {EDGE_SCAN_CAP}-edge resume scan"
                )));
            }
            for (descriptor_index, edge) in edges.edges.iter().enumerate() {
                if edge.role != EdgeRole::Out {
                    continue;
                }
                let Some(info) =
                    read_artifact_info_controlled(ledger, &edge.artifact, work, descriptor_index)
                        .map_err(|error| match error {
                        EvidenceReadError::Cancelled => cancel_err(),
                        EvidenceReadError::WorkEnvelope(error) => work_err(error),
                        EvidenceReadError::Ledger(error) => resume_ledger(
                            "reading a retained driver-state descriptor failed",
                            error,
                        ),
                    })?
                else {
                    continue;
                };
                if info.kind != schema.stage_state_kind {
                    continue;
                }
                saw_state = true;
                let state = attestor.decode_state(ledger, run, edge.artifact, work)?;
                validate_state_shape_for::<A::Stage>(&state, run, schema.cores)?;
                candidates.push((state, edge.artifact, id));
            }
            if !saw_state {
                *rejections
                    .entry("op retains no driver-state artifact".to_string())
                    .or_insert(0) += 1;
            }
        }
        match page.next_cursor {
            Some(next) => {
                if visible_ids >= MAX_SOLVE_VISIBLE_OP_IDS {
                    return Err(work_err(InvocationWorkExceeded::VisibleIds {
                        limit: MAX_SOLVE_VISIBLE_OP_IDS,
                    }));
                }
                cursor = Some(next);
            }
            None => break,
        }
    }
    if trace_export {
        eprintln!(
            "TRACE export: op scan {:.3}s ({visible_ids} ops, {} candidates)",
            t_scan.elapsed().as_secs_f64(),
            candidates.len()
        );
    }
    candidates.sort_by(|left, right| {
        right
            .0
            .completed
            .len()
            .cmp(&left.0.completed.len())
            .then_with(|| left.2.cmp(&right.2))
    });
    for (state, state_artifact, discovery_op) in candidates {
        if let Some(existing) = best.as_ref()
            && A::completed_len(existing) > state.completed.len()
        {
            break;
        }
        let t_validate = std::time::Instant::now();
        let prefix_len = state.completed.len();
        let verified =
            attestor.validate_candidate(ledger, run, state, state_artifact, discovery_op, work)?;
        if trace_export {
            eprintln!(
                "TRACE export: validate prefix {prefix_len} {:.3}s",
                t_validate.elapsed().as_secs_f64()
            );
        }
        match best.as_ref() {
            Some(existing) if A::completed_len(existing) == A::completed_len(&verified) => {
                best_is_ambiguous = true;
            }
            Some(_) => {}
            None => {
                best = Some(verified);
                best_is_ambiguous = false;
            }
        }
    }
    let best = best.ok_or_else(|| {
        if seen_for_run > 0 {
            let mut reasons: Vec<String> = rejections
                .iter()
                .map(|(reason, count)| format!("{count} x {reason}"))
                .collect();
            if reasons.is_empty() {
                reasons.push("no candidate op carried an admissible driver state".to_string());
            }
            return SolveRefusal::plain(
                "cli-solve-run-incompatible",
                format!(
                    "{} `{}` has {seen_for_run} retained operation(s) in this ledger but \
                     none is admissible under driver version {}: {}",
                    schema.run_noun,
                    run.to_hex(),
                    schema.driver_version,
                    reasons.join("; ")
                ),
                schema.incompatible_run_hint,
            );
        }
        SolveRefusal::plain(
            "cli-solve-unknown-run",
            format!(
                "no {} `{}` exists in this ledger",
                schema.run_noun,
                run.to_hex()
            ),
            schema.unknown_run_hint,
        )
    })?;
    if best_is_ambiguous {
        return Err(resume_identity(format!(
            "the ledger carries competing independently valid {}-stage checkpoints for run `{}`",
            A::completed_len(&best),
            run.to_hex()
        )));
    }
    if trace_export {
        dump_trace_spans();
    }
    Ok(best)
}

pub(crate) fn stage_discovery_rejection_for<S: StageDescriptor>(
    row: &OpRow,
    run: SolveRunId,
    schema: &PipelineSchema,
) -> Result<S, String> {
    if row.id <= 0 {
        return Err("op id is not positive".to_string());
    }
    if row.session.as_deref() != Some(run.as_bytes().as_slice()) {
        return Err("op belongs to another run".to_string());
    }
    if row.outcome.as_deref() != Some("ok") {
        return Err(format!(
            "op outcome is {}, not ok",
            row.outcome.as_deref().unwrap_or("absent")
        ));
    }
    if row.diag.is_some() {
        return Err("op carries a diagnostic".to_string());
    }
    if row.ir.len() > 512 {
        return Err("op IR exceeds the discovery size cap".to_string());
    }
    let stage = parse_stage_discovery_ir_for::<S>(&row.ir, run, schema)
        .map_err(|error| format!("discovery IR rejected ({error})"))?;
    if stage.gap_dependency().is_some() {
        return Err(format!(
            "stage `{}` has no producer under this driver",
            stage.name()
        ));
    }
    if row.t_start != i64::from(stage.ordinal()) * 2
        || row.t_end != Some(i64::from(stage.ordinal()) * 2 + 1)
    {
        return Err(format!(
            "stage `{}` timing slots do not match its ordinal",
            stage.name()
        ));
    }
    Ok(stage)
}

pub(crate) fn parse_stage_discovery_ir_for<S: StageDescriptor>(
    ir: &str,
    run: SolveRunId,
    schema: &PipelineSchema,
) -> Result<S, String> {
    let mut cursor = JsonCursor::new(ir);
    cursor.expect("{\"schema\":")?;
    if cursor.parse_string()? != schema.stage_ir_schema {
        return Err("solve-stage discovery schema does not match the driver".to_string());
    }
    cursor.expect(",\"stage\":")?;
    let stage_name = cursor.parse_string()?;
    let stage = S::all()
        .iter()
        .copied()
        .find(|stage| stage.name() == stage_name)
        .ok_or_else(|| "solve-stage discovery names an unknown stage".to_string())?;
    cursor.expect(",\"ordinal\":")?;
    if cursor.parse_u64()? != u64::from(stage.ordinal()) {
        return Err("solve-stage discovery ordinal does not match its stage".to_string());
    }
    cursor.expect(",\"run\":")?;
    if cursor.parse_string()? != run.to_hex() {
        return Err("solve-stage discovery run does not match the requested run".to_string());
    }
    cursor.expect(",\"project\":")?;
    let _ = parse_hash_string(&mut cursor, "project")?;
    cursor.expect(",\"driver_version\":")?;
    if cursor.parse_u64()? != u64::from(schema.driver_version) {
        return Err("solve-stage discovery driver version is unsupported".to_string());
    }
    cursor.expect("}")?;
    cursor.finish()?;
    Ok(stage)
}

pub(crate) fn validate_state_shape_for<S: StageDescriptor>(
    state: &SolveDriverState,
    run: SolveRunId,
    cores: u64,
) -> Result<(), SolveRefusal> {
    if state.run != *run.as_bytes() {
        return Err(resume_identity(
            "the retained driver state carries a different run identity",
        ));
    }
    #[allow(clippy::float_cmp)]
    let core_wall_mismatch = cores == 1 && state.consumed_core_s != state.consumed_wall_s;
    if !state.consumed_core_s.is_finite()
        || state.consumed_core_s < 0.0
        || !state.consumed_wall_s.is_finite()
        || state.consumed_wall_s < 0.0
        || core_wall_mismatch
    {
        return Err(resume_identity(
            "the retained driver state carries non-finite or negative consumption",
        ));
    }
    if state.completed.is_empty() || state.completed.len() > S::all().len() {
        return Err(resume_identity(format!(
            "the retained driver state has an impossible {}-stage prefix",
            state.completed.len()
        )));
    }
    let mut previous_op = 0i64;
    for (index, completed) in state.completed.iter().enumerate() {
        if completed.ordinal as usize != index {
            return Err(resume_identity(format!(
                "completed stage index {index} carries ordinal {}",
                completed.ordinal
            )));
        }
        if completed.op_id <= 0 || completed.op_id <= previous_op {
            return Err(resume_identity(format!(
                "completed stage index {index} carries non-positive or non-increasing operation id {}",
                completed.op_id
            )));
        }
        previous_op = completed.op_id;
    }
    Ok(())
}

pub(crate) fn validate_checkpoint_prefix_for(
    checkpoint: &SolveDriverState,
    final_state: &SolveDriverState,
    index: usize,
    predecessor: Option<&SolveDriverState>,
    cores: u64,
) -> Result<(), SolveRefusal> {
    #[allow(clippy::float_cmp)]
    let core_wall_mismatch = cores == 1 && checkpoint.consumed_core_s != checkpoint.consumed_wall_s;
    if checkpoint.run != final_state.run
        || checkpoint.project != final_state.project
        || checkpoint.completed.as_slice() != &final_state.completed[..=index]
        || !checkpoint.consumed_core_s.is_finite()
        || checkpoint.consumed_core_s < 0.0
        || !checkpoint.consumed_wall_s.is_finite()
        || checkpoint.consumed_wall_s < 0.0
        || core_wall_mismatch
    {
        return Err(resume_identity(format!(
            "stage {index} checkpoint does not encode the exact trusted prefix"
        )));
    }
    if let Some(predecessor) = predecessor
        && (checkpoint.consumed_core_s < predecessor.consumed_core_s
            || checkpoint.consumed_wall_s < predecessor.consumed_wall_s)
    {
        return Err(resume_identity(format!(
            "stage {index} checkpoint consumption regresses from its predecessor"
        )));
    }
    Ok(())
}

pub(crate) fn derive_run_id_bytes(
    domain: &str,
    input_hash: ContentHash,
    constellation: &str,
    workspace: &str,
    seed: u64,
    cards_root: ContentHash,
    driver_version: u32,
) -> ContentHash {
    let mut preimage = Vec::with_capacity(160);
    preimage.extend_from_slice(input_hash.as_bytes());
    push_framed(&mut preimage, constellation.as_bytes());
    push_framed(&mut preimage, workspace.as_bytes());
    preimage.extend_from_slice(&seed.to_le_bytes());
    push_framed(&mut preimage, cards_root.as_bytes());
    preimage.extend_from_slice(&driver_version.to_le_bytes());
    hash_domain(domain, &preimage)
}

pub(crate) fn derive_pipeline_capability_token<S: StageDescriptor>(
    run: SolveRunId,
    wall_s: f64,
    mem_bytes: u64,
    schema: &PipelineSchema,
) -> Result<CapabilityToken, SolveRefusal> {
    if !wall_s.is_finite() || wall_s <= 0.0 {
        return Err(SolveRefusal::plain(
            "cli-solve-budget",
            format!("project solve-time budget `{wall_s}` is not a positive finite second count"),
            "declare a positive finite `:solve-time` budget",
        ));
    }
    let mut ops: Vec<String> = S::all()
        .iter()
        .map(|stage| stage.verb().to_string())
        .collect();
    ops.push(schema.terminal_verb.to_string());
    #[allow(clippy::cast_precision_loss)]
    let core_s = wall_s * schema.cores as f64;
    Ok(CapabilityToken {
        session: SessionId(run.session_u64()),
        ops,
        core_s,
        mem_bytes,
        wall_s,
        cores: schema.cores,
        ledger_scope: format!("{}-{}", schema.ledger_scope_prefix, &run.to_hex()[..16]),
    })
}

pub(crate) fn format_stage_ir<S: StageDescriptor>(
    schema: &PipelineSchema,
    stage: S,
    run: SolveRunId,
    input_hash: ContentHash,
) -> String {
    format!(
        "{{\"schema\":{},\"stage\":{},\"ordinal\":{},\"run\":{},\"project\":{},\"driver_version\":{}}}",
        json_string(schema.stage_ir_schema),
        json_string(stage.name()),
        stage.ordinal(),
        json_string(&run.to_hex()),
        json_string(&input_hash.to_hex()),
        schema.driver_version,
    )
}

pub(crate) fn format_terminal_ir<S: StageDescriptor>(
    schema: &PipelineSchema,
    run: SolveRunId,
    input_hash: ContentHash,
) -> String {
    format!(
        "{{\"schema\":{},\"stage\":\"terminal\",\"ordinal\":{},\"run\":{},\"project\":{},\"driver_version\":{}}}",
        json_string(schema.stage_ir_schema),
        S::all().len(),
        json_string(&run.to_hex()),
        json_string(&input_hash.to_hex()),
        schema.driver_version,
    )
}

pub(crate) fn format_run_receipt_json<S: StageDescriptor>(
    run_receipt_schema: &str,
    run: SolveRunId,
    input_hash: ContentHash,
    status_name: &str,
    detail: &str,
    state: &SolveDriverState,
) -> String {
    let mut stages_json = String::new();
    for (index, stage) in state.completed.iter().enumerate() {
        if index > 0 {
            stages_json.push(',');
        }
        let name = S::from_ordinal(stage.ordinal).map_or("unknown", StageDescriptor::name);
        let _ = write!(
            stages_json,
            "{{\"stage\":{},\"op\":{},\"receipt\":{}}}",
            json_string(name),
            stage.op_id,
            json_string(&stage.receipt.to_hex()),
        );
    }
    format!(
        "{{\"schema\":{},\"run\":{},\"project_hash\":{},\"status\":{}{detail},\"stages\":[{stages_json}],\"consumed_wall_s\":{},\"consumed_core_s\":{},\"no_claim\":\"stage receipts carry their own authority; this record is run bookkeeping\"}}",
        json_string(run_receipt_schema),
        json_string(&run.to_hex()),
        json_string(&input_hash.to_hex()),
        json_string(status_name),
        state.consumed_wall_s,
        state.consumed_core_s,
    )
}

pub(crate) fn format_progress_line(
    progress_schema: &str,
    run: &str,
    stage: &str,
    ordinal: u32,
    status: &str,
    wall_s: f64,
) -> String {
    format!(
        "{{\"schema\":{},\"run\":{},\"stage\":{},\"ordinal\":{ordinal},\"status\":{},\"wall_s\":{wall_s}}}",
        json_string(progress_schema),
        json_string(run),
        json_string(stage),
        json_string(status),
    )
}

pub(crate) fn format_budget_warning_line(
    progress_schema: &str,
    run: &str,
    fraction: f64,
    used: f64,
    granted: f64,
) -> String {
    format!(
        "{{\"schema\":{},\"run\":{},\"stage\":\"budget\",\"status\":\"warning\",\"fraction\":{fraction},\"used_wall_s\":{used},\"granted_wall_s\":{granted}}}",
        json_string(progress_schema),
        json_string(run),
    )
}

pub(crate) fn finish_ledger_transaction<T>(
    ledger: &Ledger,
    result: Result<T, LedgerError>,
) -> Result<T, LedgerError> {
    match result {
        Ok(value) => match ledger.commit() {
            Ok(()) => Ok(value),
            Err(error) => {
                let _ = ledger.rollback();
                Err(error)
            }
        },
        Err(error) => {
            let _ = ledger.rollback();
            Err(error)
        }
    }
}

pub(crate) fn session_refusal(run: SolveRunId, error: &SessionError) -> SolveRefusal {
    SolveRefusal {
        code: "cli-solve-session",
        stage: None,
        what: format!("the session governor refused: {error:?}"),
        fix: "the capability token derivation and governor wiring disagree; this is a \
              driver defect"
            .to_string(),
        dependency: None,
        run: Some(run.to_hex()),
        recorded_op: None,
    }
}

pub(crate) fn ledger_refusal(
    run: SolveRunId,
    stage_name: &'static str,
    error: &LedgerError,
) -> SolveRefusal {
    SolveRefusal {
        code: "cli-solve-ledger",
        stage: Some(stage_name),
        what: format!("ledger transaction failed: {error}"),
        fix: "repair the ledger or contention failure, verify ledger lint/integrity, \
              and retry"
            .to_string(),
        dependency: None,
        run: Some(run.to_hex()),
        recorded_op: None,
    }
}

pub(crate) fn staged_refusal(
    code: &'static str,
    stage_name: &'static str,
    what: impl Into<String>,
    fix: impl Into<String>,
) -> SolveRefusal {
    SolveRefusal {
        code,
        stage: Some(stage_name),
        what: what.into(),
        fix: fix.into(),
        dependency: None,
        run: None,
        recorded_op: None,
    }
}

pub(crate) fn envelope_provenance(run: SolveRunId) -> u64 {
    run.session_u64()
}

fn push_framed(preimage: &mut Vec<u8>, bytes: &[u8]) {
    preimage.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    preimage.extend_from_slice(bytes);
}
