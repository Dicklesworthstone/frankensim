//! Six-stage durable assembly pipeline (`admit -> lower -> plan -> run ->
//! probe -> report`) on `crate::stage_driver` plus `explain_assembly_run`
//! (bead `frankensim-rc-1109-root-3gc12.11.7`, Stage 2).

#![allow(clippy::result_large_err)]

pub(crate) mod cli;
mod codec;
mod report;

use std::collections::BTreeMap;
use std::f64::consts::TAU;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::identity::ContentId;
use fs_blake3::{ContentHash, hash_bytes};
use fs_exec::solver::{
    LegacySnapshotExpectationV1, LegacySnapshotLimitsV1, LegacySnapshotV1Adapter,
    LegacySolverStateV1, codec as snap_codec,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ledger::{EdgeRole, Ledger, LedgerError};
use fs_scenario::assembly::draft::{PhysicsDomain, StepPolicy};
use fs_system::elements::Motion;
use fs_system::{
    AnalysisError, AnalysisOutcome, AnalysisRequest, LinearizeConfig, LoweredSystem,
    LoweringRegistry, ModalConfig, NodeId, OutputChannel, StateMap, StationaryConfig, System,
};

use self::codec::{AdmittedAssembly, admit_assembly_source};
use self::report::{
    AssemblyReportInputs, EntitySample, EvaluatedProbe, PROBE_SERIES_KIND, RunTrace,
    build_assembly_report, evaluate_probes,
};
use crate::cards::{CardPackKind, CardPackSet, RawCardPack};
use crate::import::json_string;
use crate::json_read::JsonValue;
use crate::solve::{
    EvidenceWork, InvocationWorkLedger, SolveDriverState, SolveOutcome, SolveRefusal, SolveRunId,
    resume_identity,
};
use crate::stage_driver::{
    DurableStageSession, PipelineSchema, ResumeAttestor, RetainedSideArtifact, SessionExplicits,
    StageDescriptor, StageExecution, StagePipeline, StageSessionView,
    derive_pipeline_capability_token, derive_run_id_bytes, discover_and_attest_resume,
    envelope_provenance, staged_refusal, validate_checkpoint_prefix_for, validate_state_shape_for,
};

/// Wire version of the assembly stage driver.
pub const ASSEMBLY_DRIVER_VERSION: u32 = 1;
/// BLAKE3 domain separating assembly run identities from thermal solve runs.
pub const ASSEMBLY_RUN_IDENTITY_DOMAIN: &str = "org.frankensim.fs-cli.assembly-run-identity.v1";

const ASSEMBLY_STATE_TYPE_ID: u64 = 0x6673_636c_6961_736d; // b"fscliasm"
const ASSEMBLY_TRACE_KIND: &str = "assembly-run-trace";
const ASSEMBLY_ANALYSIS_META_KIND: &str = "assembly-analysis-mode";

const ASSEMBLY_SCHEMA: PipelineSchema = PipelineSchema {
    driver_version: ASSEMBLY_DRIVER_VERSION,
    stage_ir_schema: "frankensim.cli.assembly-stage.v1",
    run_receipt_schema: "frankensim.cli.assembly-run-receipt.v1",
    progress_schema: "frankensim.cli.assembly-progress.v1",
    refusal_schema: "frankensim.cli.assembly-refusal.v1",
    source_artifact_kind: "assembly-source",
    stage_state_kind: "assembly-stage-state",
    stage_receipt_kind: "assembly-stage-receipt",
    run_receipt_kind: "assembly-run-receipt",
    refusal_artifact_kind: "assembly-refusal",
    cores: 1,
    terminal_verb: "assembly.terminal",
    ledger_scope_prefix: "assembly",
    session_open_prefix: "assembly-open",
    meter_prefix: "assembly-meter",
    resume_command_prefix: "frankensim run --resume",
    budget_resume_hint: "raise assembly budget and re-run",
    incompatible_run_hint: "re-run with the matching assembly driver version",
    unknown_run_hint: "verify the run id and ledger path",
    run_noun: "assembly run",
    terminal_fallback_stage: "report",
};

/// Ordered six-stage assembly workflow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyStage {
    /// Parse, Five Explicits, card admission, and `DraftAssembly::admit`.
    Admit,
    /// Lower via `LoweringRegistry::with_kernel_defaults` with receipt verification.
    Lower,
    /// Seal static `ExecutionPlan` and `CausalizationAudit`.
    Plan,
    /// Step or analyze the unchanged lowered `System` and retain ledger/trace.
    Run,
    /// Evaluate declared probes honestly from retained kernel outputs.
    Probe,
    /// Build HTML/JSON report and checker-verified format-9 evidence package.
    Report,
}

impl StageDescriptor for AssemblyStage {
    fn all() -> &'static [Self] {
        &[
            Self::Admit,
            Self::Lower,
            Self::Plan,
            Self::Run,
            Self::Probe,
            Self::Report,
        ]
    }
    fn ordinal(self) -> u32 {
        self as u32
    }
    fn name(self) -> &'static str {
        ["admit", "lower", "plan", "run", "probe", "report"][self as usize]
    }
    fn verb(self) -> &'static str {
        [
            "assembly.admit",
            "assembly.lower",
            "assembly.plan",
            "assembly.run",
            "assembly.probe",
            "assembly.report",
        ][self as usize]
    }
    fn gap_dependency(self) -> Option<&'static str> {
        None
    }
}

/// Selectable analysis mode for `--analysis transient|stationary|harmonic|modal`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyAnalysisMode {
    /// Coupled transient stepping via `System::step`.
    Transient,
    /// Stationary equilibrium via `System::run_analysis(Stationary)`.
    Stationary,
    /// Linearised frequency response via `System::run_analysis(Harmonic)`.
    Harmonic,
    /// Linearised eigenmodes via `System::run_analysis(Modal)`.
    Modal,
}

impl AssemblyAnalysisMode {
    /// Canonical lowercase mode identifier.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Transient => "transient",
            Self::Stationary => "stationary",
            Self::Harmonic => "harmonic",
            Self::Modal => "modal",
        }
    }
    pub(crate) fn parse(s: &str) -> Option<Self> {
        match s {
            "transient" => Some(Self::Transient),
            "stationary" => Some(Self::Stationary),
            "harmonic" => Some(Self::Harmonic),
            "modal" => Some(Self::Modal),
            _ => None,
        }
    }
}

/// Deterministic cancellation hook for stage-boundary kill/resume testing.
#[derive(Debug, Clone, Copy, Default)]
pub struct AssemblyCancelPlan {
    /// Trigger cancellation immediately after this stage commits durably.
    pub cancel_after: Option<AssemblyStage>,
    /// Trigger cancellation right before this stage publishes to the ledger.
    pub cancel_pre_publish: Option<AssemblyStage>,
}

/// Read-only projection of a retained assembly run for `frankensim explain <run>`.
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyExplanation {
    /// Run identifier hex.
    pub run: String,
    /// Number of completed durable stages (`1..=6`).
    pub stages_completed: usize,
    /// Selected analysis mode label.
    pub analysis: String,
    /// Lowering receipt content hash hex, when `lower` completed.
    pub lowering_receipt_hash: Option<String>,
    /// Canonical lowering receipt JSON, when `lower` completed.
    pub lowering_receipt_json: Option<String>,
    /// Static execution plan hash hex, when `plan` completed.
    pub plan_hash: Option<String>,
    /// Canonical execution plan JSON, when `plan` completed.
    pub plan_json: Option<String>,
    /// Canonical causalization audit JSON, when `plan` completed.
    pub causalization_audit_json: Option<String>,
    /// Window energy ledger JSON, when transient `run` completed.
    pub window_ledger_json: Option<String>,
    /// Completed stages `(name, op_id, receipt_hash)`.
    pub stages: Vec<(String, i64, String)>,
}

struct AssemblyStateEnvelope(SolveDriverState);
impl LegacySolverStateV1 for AssemblyStateEnvelope {
    const SCHEMA_VERSION_V1: u32 = 1;
    const TYPE_ID_V1: u64 = ASSEMBLY_STATE_TYPE_ID;
    fn encode_v1(&self, enc: &mut snap_codec::Enc) {
        self.0.encode_v1(enc);
    }
    fn decode_v1(dec: &mut snap_codec::Dec<'_>) -> Result<Self, snap_codec::CodecError> {
        SolveDriverState::decode_v1(dec).map(Self)
    }
}

struct LoweredStateMap {
    coords: Vec<(NodeId, bool, String)>, // (node, is_motion_velocity, label)
}
impl LoweredStateMap {
    fn from_lowered(lowered: &LoweredSystem) -> Self {
        let mut coords = Vec::new();
        for (name, &node) in &lowered.node_map {
            if lowered.system.state::<Motion>(node).is_some() {
                // `Motion::position_m` is an output integral absent from port laws;
                // mapping only `velocity_m_s` keeps the state Jacobian nonsingular
                // (matching `Mech2Map` / `ThermoMechMap` in `fs-system`).
                coords.push((node, true, format!("{name}.velocity_m_s")));
            } else if lowered.system.state::<f64>(node).is_some() {
                coords.push((node, false, format!("{name}.state")));
            }
        }
        Self { coords }
    }
}
impl StateMap for LoweredStateMap {
    fn dim(&self) -> usize {
        self.coords.len()
    }
    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        if out.len() != self.coords.len() {
            return Err(AnalysisError::Invalid {
                what: "state length mismatch",
            });
        }
        for (slot, (node, is_motion, _)) in out.iter_mut().zip(&self.coords) {
            *slot = if *is_motion {
                system
                    .state::<Motion>(*node)
                    .ok_or(AnalysisError::Invalid {
                        what: "missing Motion",
                    })?
                    .velocity_m_s
            } else {
                *system.state::<f64>(*node).ok_or(AnalysisError::Invalid {
                    what: "missing f64",
                })?
            };
        }
        Ok(())
    }
    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError> {
        if state.len() != self.coords.len() {
            return Err(AnalysisError::Invalid {
                what: "state length mismatch",
            });
        }
        for ((node, is_motion, _), &val) in self.coords.iter().zip(state) {
            if *is_motion {
                let mut motion = *system
                    .state::<Motion>(*node)
                    .ok_or(AnalysisError::Invalid {
                        what: "missing Motion",
                    })?;
                motion.velocity_m_s = val;
                system.set_accepted_state(*node, motion)?;
            } else {
                system.set_accepted_state(*node, val)?;
            }
        }
        Ok(())
    }
    fn labels(&self) -> Vec<String> {
        self.coords.iter().map(|c| c.2.clone()).collect()
    }
}

#[derive(Default)]
struct AssemblyCtx {
    lowered: Option<LoweredSystem>,
    lowering_hash: Option<ContentHash>,
    plan_hash: Option<String>,
    plan_artifact: Option<ContentHash>,
    ledger_artifact: Option<ContentHash>,
    energy_residual_j: Option<f64>,
    trace: Option<RunTrace>,
    probe_series_artifact: Option<ContentHash>,
    probes: Vec<EvaluatedProbe>,
    pending_metrics: Vec<(i64, String, f64)>,
}

enum StageCommit {
    Admit,
    Lower {
        lowered: Box<LoweredSystem>,
        receipt_hash: ContentHash,
    },
    Plan {
        plan_hash: String,
        plan_artifact: ContentHash,
    },
    Run {
        trace: RunTrace,
        ledger_artifact: Option<ContentHash>,
        energy_residual_j: Option<f64>,
    },
    Probe {
        series_artifact: ContentHash,
        probes: Vec<EvaluatedProbe>,
    },
    Report,
}

struct AssemblyPipeline<'a> {
    admitted: &'a AdmittedAssembly,
    packs: &'a CardPackSet,
    analysis: AssemblyAnalysisMode,
    run: SolveRunId,
    gate: CancelGate,
    cancel_plan: AssemblyCancelPlan,
}

impl AssemblyPipeline<'_> {
    fn run_transient(
        &self,
        lowered: &mut LoweredSystem,
    ) -> Result<(RunTrace, fs_system::LedgerExport), SolveRefusal> {
        let tp = &self.admitted.draft.input().time_plan;
        let dt_s = match tp.step_policy {
            StepPolicy::Fixed { dt } => dt.value,
            StepPolicy::Adaptive { dt_initial, .. } => dt_initial.value,
        };
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let steps = ((tp.horizon.value / dt_s).round().max(1.0) as usize).min(200_000);
        let mut times_s = Vec::with_capacity(steps + 1);
        let mut entities: BTreeMap<String, Vec<EntitySample>> = self
            .admitted
            .draft
            .input()
            .parts
            .iter()
            .map(|p| (p.name.clone(), Vec::with_capacity(steps + 1)))
            .collect();
        let sample = |sys: &System, ents: &mut BTreeMap<String, Vec<EntitySample>>| {
            for part in &self.admitted.draft.input().parts {
                let mut s = EntitySample {
                    position_m: None,
                    velocity_m_s: None,
                    temperature_k: None,
                };
                for &node in lowered.entity_nodes.get(&part.entity).into_iter().flatten() {
                    if let Some(m) = sys.state::<Motion>(node) {
                        s.position_m = Some(m.position_m);
                        s.velocity_m_s = Some(m.velocity_m_s);
                    } else if let Some(&v) = sys.state::<f64>(node) {
                        if part.physics.contains(&PhysicsDomain::Thermal) {
                            s.temperature_k = Some(v);
                        } else if s.position_m.is_none() {
                            s.position_m = Some(v);
                        }
                    }
                }
                if let Some(vec) = ents.get_mut(&part.name) {
                    vec.push(s);
                }
            }
        };
        times_s.push(lowered.system.time_s());
        sample(&lowered.system, &mut entities);
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(
                &self.gate,
                arena,
                StreamKey {
                    seed: self.admitted.seed,
                    kernel_id: 117,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            for _ in 0..steps {
                if self.gate.is_requested() {
                    return Err(SolveRefusal::plain(
                        "cli-solve-cancelled",
                        "cancelled during stepping",
                        "resume",
                    ));
                }
                lowered.system.step(dt_s, &cx).map_err(|e| {
                    staged_refusal(
                        "cli-assembly-step",
                        "run",
                        format!("kernel step failed: {e}"),
                        "check step size",
                    )
                })?;
                times_s.push(lowered.system.time_s());
                sample(&lowered.system, &mut entities);
            }
            Ok(())
        })?;
        Ok((
            RunTrace {
                analysis: "transient".into(),
                times_s,
                entities,
                modal_modes: Vec::new(),
                harmonic_points: Vec::new(),
                stationary_norm: None,
                stationary_state: Vec::new(),
            },
            lowered.system.ledger().export(),
        ))
    }

    fn run_non_transient(&self, lowered: &mut LoweredSystem) -> Result<RunTrace, SolveRefusal> {
        let smap = LoweredStateMap::from_lowered(lowered);
        let outcome = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(
                &self.gate,
                arena,
                StreamKey {
                    seed: self.admitted.seed,
                    kernel_id: 117,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            match self.analysis {
                AssemblyAnalysisMode::Stationary => lowered.system.run_analysis(
                    AnalysisRequest::Stationary {
                        state_map: &smap,
                        config: StationaryConfig::default(),
                    },
                    &cx,
                ),
                AssemblyAnalysisMode::Modal => lowered.system.run_analysis(
                    AnalysisRequest::Modal {
                        state_map: &smap,
                        config: ModalConfig::default(),
                    },
                    &cx,
                ),
                AssemblyAnalysisMode::Harmonic => {
                    let inputs = lowered.system.boundary_channels();
                    if inputs.is_empty() {
                        return Err(AnalysisError::Invalid {
                            what: "no boundary channels available for harmonic response",
                        });
                    }
                    let outputs: Vec<OutputChannel> =
                        (0..smap.dim()).map(OutputChannel::State).collect();
                    let hz = if self.admitted.harmonic_frequencies_hz.is_empty() {
                        vec![1.0, 5.0, 10.0]
                    } else {
                        self.admitted.harmonic_frequencies_hz.clone()
                    };
                    let omegas: Vec<f64> = hz.iter().map(|f| f * TAU).collect();
                    lowered.system.run_analysis(
                        AnalysisRequest::Harmonic {
                            state_map: &smap,
                            inputs: &inputs,
                            outputs: &outputs,
                            frequencies_rad_s: &omegas,
                            config: LinearizeConfig::default(),
                        },
                        &cx,
                    )
                }
                AssemblyAnalysisMode::Transient => unreachable!(),
            }
        });
        let outcome = outcome.map_err(|err| match err {
            AnalysisError::MissingTier { analysis, node, required } => staged_refusal(
                "cli-assembly-analysis-capability",
                "run",
                format!("analysis `{analysis}` refused: component `{node}` lacks required capability tier `{required}`"),
                format!("implement `{required}` on `{node}` or choose a supported analysis"),
            ),
            other => staged_refusal("cli-assembly-analysis", "run", format!("analysis `{}` failed: {other}", self.analysis.as_str()), "check analysis setup"),
        })?;
        let mut trace = RunTrace {
            analysis: self.analysis.as_str().into(),
            times_s: Vec::new(),
            entities: BTreeMap::new(),
            modal_modes: Vec::new(),
            harmonic_points: Vec::new(),
            stationary_norm: None,
            stationary_state: Vec::new(),
        };
        match outcome {
            AnalysisOutcome::Modal(m) => {
                trace.modal_modes = m
                    .modes
                    .iter()
                    .take(self.admitted.modal_max_modes)
                    .map(|x| (x.natural_frequency_rad_s, x.frequency_hz, x.damping_ratio))
                    .collect();
            }
            AnalysisOutcome::Harmonic(h) => {
                for (fi, &w) in h.frequencies_rad_s.iter().enumerate() {
                    let c = h.continuous(fi, 0, 0);
                    trace
                        .harmonic_points
                        .push((w / TAU, c.norm_sq().sqrt(), c.im.atan2(c.re)));
                }
            }
            AnalysisOutcome::Stationary(s) => {
                trace.stationary_norm = Some(s.residual_norm);
                trace.stationary_state = smap.labels().into_iter().zip(s.state).collect();
            }
            _ => {}
        }
        Ok(trace)
    }
}

impl StagePipeline for AssemblyPipeline<'_> {
    type Stage = AssemblyStage;
    type Context = AssemblyCtx;
    type StageCommit = StageCommit;

    fn schema(&self) -> &'static PipelineSchema {
        &ASSEMBLY_SCHEMA
    }
    fn canonical_source(&self) -> &[u8] {
        &self.admitted.canonical_source
    }
    fn seal_state(&self, run: SolveRunId, state: &SolveDriverState) -> Vec<u8> {
        LegacySnapshotV1Adapter::<AssemblyStateEnvelope>::seal(
            &AssemblyStateEnvelope(state.clone()),
            envelope_provenance(run),
        )
    }
    fn declaration_gap(&self, _: Self::Stage) -> Option<(&'static str, String, String)> {
        None
    }
    fn is_cancelled(&self) -> bool {
        self.gate.is_requested()
    }
    fn pre_publish_checkpoint(&self, stage: Self::Stage) -> Result<(), ()> {
        if self.cancel_plan.cancel_pre_publish == Some(stage) {
            self.gate.request();
        }
        (!self.gate.is_requested()).then_some(()).ok_or(())
    }

    #[allow(clippy::too_many_lines)]
    fn execute_stage(
        &mut self,
        stage: Self::Stage,
        _: &StageSessionView<'_>,
        state: &SolveDriverState,
        ctx: &mut Self::Context,
    ) -> Result<StageExecution<Self::StageCommit>, SolveRefusal> {
        ctx.pending_metrics.clear();
        match stage {
            AssemblyStage::Admit => {
                let bytes = self.analysis.as_str().as_bytes().to_vec();
                let artifact = hash_bytes(&bytes);
                let inp = self.admitted.draft.input();
                let receipt_json = format!(
                    "{{\"schema\":\"frankensim.cli.assembly-admit.v1\",\"assembly\":{},\"analysis\":{},\"source_hash\":{},\"parts\":{},\"joints\":{},\"interfaces\":{},\"sources\":{},\"probes\":{},\"card_packs\":{}}}",
                    json_string(&self.admitted.name),
                    json_string(self.analysis.as_str()),
                    json_string(&self.admitted.source_hash.to_hex()),
                    inp.parts.len(),
                    inp.joints.len(),
                    inp.interfaces.len(),
                    inp.sources.len(),
                    inp.probes.len(),
                    self.packs.len(),
                );
                Ok(StageExecution {
                    receipt_json,
                    side_artifacts: vec![RetainedSideArtifact {
                        kind: ASSEMBLY_ANALYSIS_META_KIND,
                        artifact,
                        bytes,
                    }],
                    on_commit: StageCommit::Admit,
                })
            }
            AssemblyStage::Lower => {
                let lowered = LoweringRegistry::with_kernel_defaults()
                    .lower(&self.admitted.draft, &self.admitted.card_store)
                    .map_err(|e| {
                        staged_refusal(
                            "cli-assembly-lower",
                            "lower",
                            format!("lowering failed: {e}"),
                            "check assembly cards and domains",
                        )
                    })?;
                lowered.receipt.verify_completeness().map_err(|e| {
                    staged_refusal(
                        "cli-assembly-lower",
                        "lower",
                        format!("lowering receipt incomplete: {e}"),
                        "ensure provenance completeness",
                    )
                })?;
                let bytes = lowered.receipt.canonical_json().into_bytes();
                let artifact = hash_bytes(&bytes);
                let receipt_hash = lowered.receipt.receipt_hash;
                #[allow(clippy::cast_precision_loss)]
                ctx.pending_metrics.extend([
                    (
                        2,
                        "lowering.components".into(),
                        lowered.receipt.components.len() as f64,
                    ),
                    (
                        2,
                        "lowering.interfaces".into(),
                        lowered.receipt.interfaces.len() as f64,
                    ),
                    (
                        2,
                        "lowering.joints".into(),
                        lowered.receipt.joints.len() as f64,
                    ),
                    (
                        2,
                        "lowering.property_receipts".into(),
                        lowered.receipt.total_property_receipts() as f64,
                    ),
                ]);
                let receipt_json = format!(
                    "{{\"schema\":\"frankensim.cli.assembly-lower.v1\",\"lowering_receipt_hash\":{},\"lowering_artifact\":{},\"components\":{},\"interfaces\":{},\"joints\":{},\"property_receipts\":{}}}",
                    json_string(&receipt_hash.to_hex()),
                    json_string(&artifact.to_hex()),
                    lowered.receipt.components.len(),
                    lowered.receipt.interfaces.len(),
                    lowered.receipt.joints.len(),
                    lowered.receipt.total_property_receipts(),
                );
                Ok(StageExecution {
                    receipt_json,
                    side_artifacts: vec![RetainedSideArtifact {
                        kind: Ledger::LOWERING_RECEIPT_ARTIFACT_KIND,
                        artifact,
                        bytes,
                    }],
                    on_commit: StageCommit::Lower {
                        lowered: Box::new(lowered),
                        receipt_hash,
                    },
                })
            }
            AssemblyStage::Plan => {
                let lowered = ctx
                    .lowered
                    .as_ref()
                    .ok_or_else(|| resume_identity("missing lowered system"))?;
                let plan_hash = lowered.plan.plan_hash_hex();
                let (pb, ab) = (
                    lowered.plan.to_canonical_json().into_bytes(),
                    lowered.causalization_audit.to_canonical_json().into_bytes(),
                );
                let (pa, aa) = (hash_bytes(&pb), hash_bytes(&ab));
                #[allow(clippy::cast_precision_loss)]
                ctx.pending_metrics.extend([
                    (4, "plan.blocks".into(), lowered.plan.block_count as f64),
                    (4, "plan.loops".into(), lowered.plan.loop_count as f64),
                    (
                        4,
                        "plan.kept_unknowns".into(),
                        lowered.plan.total_kept_unknowns() as f64,
                    ),
                    (
                        4,
                        "plan.workspace_bytes".into(),
                        lowered.plan.total_workspace_bytes as f64,
                    ),
                ]);
                let receipt_json = format!(
                    "{{\"schema\":\"frankensim.cli.assembly-plan.v1\",\"plan_hash\":{},\"plan_artifact\":{},\"causalization_artifact\":{},\"blocks\":{},\"loops\":{},\"nodes\":{},\"connections\":{}}}",
                    json_string(&plan_hash),
                    json_string(&pa.to_hex()),
                    json_string(&aa.to_hex()),
                    lowered.plan.block_count,
                    lowered.plan.loop_count,
                    lowered.plan.node_count,
                    lowered.plan.connection_count,
                );
                Ok(StageExecution {
                    receipt_json,
                    side_artifacts: vec![
                        RetainedSideArtifact {
                            kind: Ledger::EXECUTION_PLAN_ARTIFACT_KIND,
                            artifact: pa,
                            bytes: pb,
                        },
                        RetainedSideArtifact {
                            kind: Ledger::CAUSALIZATION_AUDIT_ARTIFACT_KIND,
                            artifact: aa,
                            bytes: ab,
                        },
                    ],
                    on_commit: StageCommit::Plan {
                        plan_hash,
                        plan_artifact: pa,
                    },
                })
            }
            AssemblyStage::Run => {
                let lowered = ctx
                    .lowered
                    .as_mut()
                    .ok_or_else(|| resume_identity("missing lowered system"))?;
                let (trace, sides, ledger_artifact, energy_residual_j) =
                    if self.analysis == AssemblyAnalysisMode::Transient {
                        let (trace, export) = self.run_transient(lowered)?;
                        let defect = lowered.system.ledger().energy_residual_j();
                        let wb = export.canonical_json.into_bytes();
                        let wh = hash_bytes(&wb);
                        ctx.pending_metrics
                            .extend(export.metrics.into_iter().map(|(k, v)| (6, k, v)));
                        let tb = trace.to_json().into_bytes();
                        let th = hash_bytes(&tb);
                        (
                            trace,
                            vec![
                                RetainedSideArtifact {
                                    kind: Ledger::WINDOW_LEDGER_ARTIFACT_KIND,
                                    artifact: wh,
                                    bytes: wb,
                                },
                                RetainedSideArtifact {
                                    kind: ASSEMBLY_TRACE_KIND,
                                    artifact: th,
                                    bytes: tb,
                                },
                            ],
                            Some(wh),
                            Some(defect),
                        )
                    } else {
                        let trace = self.run_non_transient(lowered)?;
                        let tb = trace.to_json().into_bytes();
                        let th = hash_bytes(&tb);
                        (
                            trace,
                            vec![RetainedSideArtifact {
                                kind: ASSEMBLY_TRACE_KIND,
                                artifact: th,
                                bytes: tb,
                            }],
                            None,
                            None,
                        )
                    };
                let receipt_json = format!(
                    "{{\"schema\":\"frankensim.cli.assembly-run.v1\",\"analysis\":{},\"trace_artifact\":{},\"window_ledger\":{},\"steps\":{}}}",
                    json_string(self.analysis.as_str()),
                    json_string(&sides.last().unwrap().artifact.to_hex()),
                    ledger_artifact.map_or_else(|| "null".into(), |h| json_string(&h.to_hex())),
                    trace.times_s.len().saturating_sub(1),
                );
                Ok(StageExecution {
                    receipt_json,
                    side_artifacts: sides,
                    on_commit: StageCommit::Run {
                        trace,
                        ledger_artifact,
                        energy_residual_j,
                    },
                })
            }
            AssemblyStage::Probe => {
                let trace = ctx
                    .trace
                    .as_ref()
                    .ok_or_else(|| resume_identity("missing run trace"))?;
                let (receipt_json, side, probes) = evaluate_probes(
                    &self.admitted.draft.input().probes,
                    &self.admitted.entity_names,
                    trace,
                )?;
                let series_artifact = side.artifact;
                Ok(StageExecution {
                    receipt_json,
                    side_artifacts: vec![side],
                    on_commit: StageCommit::Probe {
                        series_artifact,
                        probes,
                    },
                })
            }
            AssemblyStage::Report => {
                let stage_receipts: Vec<(&'static str, u32, ContentHash, String)> = state
                    .completed
                    .iter()
                    .map(|d| {
                        let st =
                            AssemblyStage::from_ordinal(d.ordinal).unwrap_or(AssemblyStage::Admit);
                        (
                            st.name(),
                            d.ordinal,
                            d.receipt,
                            format!("stage `{}` receipt", st.name()),
                        )
                    })
                    .collect();
                let (receipt_json, side_artifacts) =
                    build_assembly_report(&AssemblyReportInputs {
                        run: self.run,
                        assembly_name: &self.admitted.name,
                        source_hash: self.admitted.source_hash,
                        constellation: &self.admitted.constellation,
                        workspace: &self.admitted.workspace,
                        seed: self.admitted.seed,
                        wall_budget_s: self.admitted.wall_budget_s,
                        mem_budget_bytes: self.admitted.mem_budget_bytes,
                        analysis: self.analysis.as_str(),
                        lowering_hash: ctx
                            .lowering_hash
                            .ok_or_else(|| resume_identity("missing lowering hash"))?,
                        plan_hash: ctx
                            .plan_hash
                            .as_deref()
                            .ok_or_else(|| resume_identity("missing plan hash"))?,
                        plan_artifact: ctx
                            .plan_artifact
                            .ok_or_else(|| resume_identity("missing plan artifact"))?,
                        ledger_artifact: ctx.ledger_artifact,
                        energy_residual_j: ctx.energy_residual_j,
                        probe_series_artifact: ctx
                            .probe_series_artifact
                            .ok_or_else(|| resume_identity("missing probe series"))?,
                        probes: &ctx.probes,
                        stage_receipts: &stage_receipts,
                    })?;
                Ok(StageExecution {
                    receipt_json,
                    side_artifacts,
                    on_commit: StageCommit::Report,
                })
            }
        }
    }

    fn commit_stage(&mut self, ctx: &mut Self::Context, commit: Self::StageCommit) {
        let stage = match commit {
            StageCommit::Admit => AssemblyStage::Admit,
            StageCommit::Lower {
                lowered,
                receipt_hash,
            } => {
                ctx.lowered = Some(*lowered);
                ctx.lowering_hash = Some(receipt_hash);
                AssemblyStage::Lower
            }
            StageCommit::Plan {
                plan_hash,
                plan_artifact,
            } => {
                ctx.plan_hash = Some(plan_hash);
                ctx.plan_artifact = Some(plan_artifact);
                AssemblyStage::Plan
            }
            StageCommit::Run {
                trace,
                ledger_artifact,
                energy_residual_j,
            } => {
                ctx.trace = Some(trace);
                ctx.ledger_artifact = ledger_artifact;
                ctx.energy_residual_j = energy_residual_j;
                AssemblyStage::Run
            }
            StageCommit::Probe {
                series_artifact,
                probes,
            } => {
                ctx.probe_series_artifact = Some(series_artifact);
                ctx.probes = probes;
                AssemblyStage::Probe
            }
            StageCommit::Report => AssemblyStage::Report,
        };
        if self.cancel_plan.cancel_after == Some(stage) {
            self.gate.request();
        }
    }

    fn link_initial_extra_inputs(&self, ledger: &Ledger, op: i64) -> Result<(), LedgerError> {
        for pack in self.packs.iter() {
            let art = ledger.put_artifact(pack.kind().artifact_kind(), pack.bytes(), None)?;
            ledger.link(op, &art.hash, EdgeRole::In)?;
        }
        Ok(())
    }

    fn link_stage_in_edges(
        &self,
        ledger: &Ledger,
        op: i64,
        stage: Self::Stage,
        before: &SolveDriverState,
        ctx: &Self::Context,
    ) -> Result<(), LedgerError> {
        if let Some(prev) = before.completed.last() {
            ledger.link(op, &prev.receipt, EdgeRole::In)?;
        }
        if stage == AssemblyStage::Report {
            for done in &before.completed {
                let _ = ledger.link(op, &done.receipt, EdgeRole::In);
            }
        }
        for (t, name, val) in &ctx.pending_metrics {
            ledger.record_metric(op, *t, name, *val)?;
        }
        Ok(())
    }
}

fn find_edge(
    ledger: &Ledger,
    op: i64,
    role: EdgeRole,
    kind: &str,
) -> Result<Option<ContentHash>, SolveRefusal> {
    let edges = ledger
        .op_artifact_edges_bounded(op, 64)
        .map_err(|e| resume_identity(format!("edges error: {e}")))?;
    for edge in edges.edges {
        if edge.role == role
            && ledger
                .artifact_info(&edge.artifact)
                .map_err(|e| resume_identity(format!("info error: {e}")))?
                .is_some_and(|i| i.kind == kind)
        {
            return Ok(Some(edge.artifact));
        }
    }
    Ok(None)
}

fn read_utf8(ledger: &Ledger, hash: &ContentHash) -> Result<String, SolveRefusal> {
    let bytes = ledger
        .get_artifact(hash)
        .map_err(|e| resume_identity(format!("read error: {e}")))?
        .ok_or_else(|| resume_identity("missing artifact"))?;
    String::from_utf8(bytes).map_err(|_| resume_identity("non-UTF-8 artifact"))
}

struct VerifiedAssemblyResume {
    state: SolveDriverState,
    state_artifact: ContentHash,
    admitted: AdmittedAssembly,
    packs: CardPackSet,
    analysis: AssemblyAnalysisMode,
    ctx: AssemblyCtx,
}

struct AssemblyResumeAttestor;
impl ResumeAttestor for AssemblyResumeAttestor {
    type Stage = AssemblyStage;
    type Verified = VerifiedAssemblyResume;

    fn schema(&self) -> &'static PipelineSchema {
        &ASSEMBLY_SCHEMA
    }
    fn decode_state(
        &self,
        ledger: &Ledger,
        run: SolveRunId,
        artifact: ContentHash,
        _: EvidenceWork<'_>,
    ) -> Result<SolveDriverState, SolveRefusal> {
        let bytes = ledger
            .get_artifact(&artifact)
            .map_err(|e| resume_identity(format!("checkpoint read: {e}")))?
            .ok_or_else(|| resume_identity("missing checkpoint"))?;
        let cid = ContentId::parse_slice(artifact.as_bytes())
            .ok_or_else(|| resume_identity("bad checkpoint hash"))?;
        let exp = LegacySnapshotExpectationV1::new(
            cid,
            ASSEMBLY_STATE_TYPE_ID,
            1,
            envelope_provenance(run),
        );
        LegacySnapshotV1Adapter::<AssemblyStateEnvelope>::open_expected(
            &bytes,
            exp,
            LegacySnapshotLimitsV1::new(1024 * 1024, 64 * 1024),
            || false,
        )
        .map(|o| o.into_parts().0.0)
        .map_err(|e| resume_identity(format!("invalid assembly checkpoint: {e:?}")))
    }

    #[allow(clippy::too_many_lines)]
    fn validate_candidate(
        &mut self,
        ledger: &Ledger,
        run: SolveRunId,
        state: SolveDriverState,
        state_artifact: ContentHash,
        discovery_op: i64,
        work: EvidenceWork<'_>,
    ) -> Result<Self::Verified, SolveRefusal> {
        validate_state_shape_for::<AssemblyStage>(&state, run, ASSEMBLY_SCHEMA.cores)?;
        let op0 = state.completed[0].op_id;
        let src_h = find_edge(
            ledger,
            op0,
            EdgeRole::In,
            ASSEMBLY_SCHEMA.source_artifact_kind,
        )?
        .ok_or_else(|| resume_identity("stage 0 missing source"))?;
        let src_bytes = ledger
            .get_artifact(&src_h)
            .map_err(|e| resume_identity(format!("{e}")))?
            .ok_or_else(|| resume_identity("missing source bytes"))?;
        let meta_h = find_edge(ledger, op0, EdgeRole::Out, ASSEMBLY_ANALYSIS_META_KIND)?
            .ok_or_else(|| resume_identity("missing analysis mode"))?;
        let analysis = AssemblyAnalysisMode::parse(read_utf8(ledger, &meta_h)?.trim())
            .ok_or_else(|| resume_identity("unknown analysis mode"))?;
        let mut raw_packs = Vec::new();
        for edge in ledger
            .op_artifact_edges_bounded(op0, 64)
            .map_err(|e| resume_identity(format!("{e}")))?
            .edges
        {
            if edge.role == EdgeRole::In
                && let Some(info) = ledger
                    .artifact_info(&edge.artifact)
                    .map_err(|e| resume_identity(format!("{e}")))?
                && let Some(kind) = CardPackKind::from_artifact_kind(&info.kind)
            {
                let bytes = ledger
                    .get_artifact(&edge.artifact)
                    .map_err(|e| resume_identity(format!("{e}")))?
                    .ok_or_else(|| resume_identity("missing pack bytes"))?;
                raw_packs.push(RawCardPack {
                    kind,
                    source: edge.artifact.to_hex(),
                    bytes,
                    expect: None,
                });
            }
        }
        let packs = CardPackSet::admit(raw_packs)
            .map_err(|e| resume_identity(format!("pack set invalid: {}", e.what)))?;
        let admitted = admit_assembly_source(&src_bytes, &packs)?;

        let (mut ctx, mut prev) = (AssemblyCtx::default(), None);
        for (idx, done) in state.completed.iter().enumerate() {
            let n = ledger
                .op_artifact_edges_bounded(done.op_id, 64)
                .map_err(|e| resume_identity(format!("{e}")))?
                .edges
                .len();
            if ledger
                .op_artifact_edge_seal(done.op_id)
                .map_err(|e| resume_identity(format!("{e}")))?
                != Some(n)
            {
                return Err(resume_identity(format!("stage {idx} edge seal mismatch")));
            }
            let cp_h = find_edge(
                ledger,
                done.op_id,
                EdgeRole::Out,
                ASSEMBLY_SCHEMA.stage_state_kind,
            )?
            .ok_or_else(|| resume_identity("missing checkpoint"))?;
            let cp = self.decode_state(ledger, run, cp_h, work)?;
            validate_checkpoint_prefix_for(&cp, &state, idx, prev.as_ref(), ASSEMBLY_SCHEMA.cores)?;
            if idx + 1 == state.completed.len()
                && (discovery_op != done.op_id || cp_h != state_artifact)
            {
                return Err(resume_identity("final checkpoint mismatch"));
            }
            let stage = AssemblyStage::from_ordinal(done.ordinal)
                .ok_or_else(|| resume_identity("bad ordinal"))?;
            let parsed = JsonValue::parse(&read_utf8(ledger, &done.receipt)?)
                .map_err(|e| resume_identity(format!("bad receipt json: {e}")))?;
            match stage {
                AssemblyStage::Lower => {
                    let lh = ContentHash::from_hex(
                        parsed
                            .str_field("lowering_receipt_hash")
                            .ok_or_else(|| resume_identity("missing lowering hash"))?,
                    )
                    .ok_or_else(|| resume_identity("bad lowering hash"))?;
                    ctx.lowering_hash = Some(lh);
                    if state.completed.len() < 4 {
                        let lowered = LoweringRegistry::with_kernel_defaults()
                            .lower(&admitted.draft, &admitted.card_store)
                            .map_err(|e| resume_identity(format!("re-lowering failed: {e}")))?;
                        if lowered.receipt.receipt_hash != lh {
                            return Err(resume_identity(
                                "resume re-lowering receipt hash diverged",
                            ));
                        }
                        ctx.lowered = Some(lowered);
                    }
                }
                AssemblyStage::Plan => {
                    ctx.plan_hash = Some(
                        parsed
                            .str_field("plan_hash")
                            .ok_or_else(|| resume_identity("missing plan_hash"))?
                            .to_string(),
                    );
                    ctx.plan_artifact = Some(
                        find_edge(
                            ledger,
                            done.op_id,
                            EdgeRole::Out,
                            Ledger::EXECUTION_PLAN_ARTIFACT_KIND,
                        )?
                        .ok_or_else(|| resume_identity("missing plan artifact"))?,
                    );
                }
                AssemblyStage::Run => {
                    let ta = find_edge(ledger, done.op_id, EdgeRole::Out, ASSEMBLY_TRACE_KIND)?
                        .ok_or_else(|| resume_identity("missing trace artifact"))?;
                    ctx.trace = Some(RunTrace::from_json(&read_utf8(ledger, &ta)?)?);
                    ctx.ledger_artifact = find_edge(
                        ledger,
                        done.op_id,
                        EdgeRole::Out,
                        Ledger::WINDOW_LEDGER_ARTIFACT_KIND,
                    )?;
                    if let Some(la) = ctx.ledger_artifact {
                        ctx.energy_residual_j = JsonValue::parse(&read_utf8(ledger, &la)?)
                            .ok()
                            .and_then(|v| v.get("totals")?.f64_field("energy_residual_j"));
                    }
                }
                AssemblyStage::Probe => {
                    let pa = find_edge(ledger, done.op_id, EdgeRole::Out, PROBE_SERIES_KIND)?
                        .ok_or_else(|| resume_identity("missing probe series"))?;
                    ctx.probe_series_artifact = Some(pa);
                    let pj = JsonValue::parse(&read_utf8(ledger, &pa)?)
                        .map_err(|e| resume_identity(format!("bad probe series: {e}")))?;
                    for item in pj
                        .get("probes")
                        .and_then(JsonValue::as_array)
                        .into_iter()
                        .flatten()
                    {
                        let values = item
                            .get("values")
                            .and_then(JsonValue::as_array)
                            .map(|a| a.iter().filter_map(JsonValue::as_f64).collect())
                            .unwrap_or_default();
                        ctx.probes.push(EvaluatedProbe {
                            name: item.str_field("name").unwrap_or("probe").to_string(),
                            kind: item.str_field("kind").unwrap_or("retained").to_string(),
                            unit: item.str_field("unit").unwrap_or("SI").to_string(),
                            final_value: item.f64_field("final").unwrap_or(0.0),
                            peak_value: item.f64_field("peak").unwrap_or(0.0),
                            values,
                        });
                    }
                }
                AssemblyStage::Admit | AssemblyStage::Report => {}
            }
            prev = Some(cp);
        }
        Ok(VerifiedAssemblyResume {
            state,
            state_artifact,
            admitted,
            packs,
            analysis,
            ctx,
        })
    }
    fn completed_len(v: &Self::Verified) -> usize {
        v.state.completed.len()
    }
}

fn assembly_run_id(
    admitted: &AdmittedAssembly,
    packs: &CardPackSet,
    analysis: AssemblyAnalysisMode,
) -> SolveRunId {
    let mut pre = admitted.source_hash.as_bytes().to_vec();
    pre.extend_from_slice(analysis.as_str().as_bytes());
    let h = derive_run_id_bytes(
        ASSEMBLY_RUN_IDENTITY_DOMAIN,
        hash_bytes(&pre),
        &admitted.constellation,
        &admitted.workspace,
        admitted.seed,
        packs.root(),
        ASSEMBLY_DRIVER_VERSION,
    );
    SolveRunId::parse_hex(&h.to_hex()).expect("content hash is valid 64-hex run id")
}

fn open_session<'a>(
    ledger: &'a Ledger,
    clock: &'a mut dyn FnMut() -> f64,
    admitted: &AdmittedAssembly,
    packs: &CardPackSet,
    analysis: AssemblyAnalysisMode,
    run: SolveRunId,
) -> Result<DurableStageSession<'a>, SolveRefusal> {
    let token = derive_pipeline_capability_token::<AssemblyStage>(
        run,
        admitted.wall_budget_s,
        admitted.mem_budget_bytes,
        &ASSEMBLY_SCHEMA,
    )?;
    let explicits = SessionExplicits {
        seed: admitted.seed.to_le_bytes(),
        versions_json: format!(
            "{{\"constellation\":{},\"workspace\":{},\"driver_version\":{},\"analysis\":{}}}",
            json_string(&admitted.constellation),
            json_string(&admitted.workspace),
            ASSEMBLY_DRIVER_VERSION,
            json_string(analysis.as_str())
        ),
        budget_json: format!(
            "{{\"wall_s\":{},\"mem_bytes\":{}}}",
            admitted.wall_budget_s, admitted.mem_budget_bytes
        ),
        capability_json: format!(
            "{{\"cores\":{},\"card_pack_root\":{},\"ops_count\":{}}}",
            admitted.cores,
            json_string(&packs.root().to_hex()),
            admitted.granted_ops.len()
        ),
    };
    DurableStageSession::open(
        ledger,
        clock,
        run,
        admitted.source_hash,
        &ASSEMBLY_SCHEMA,
        token,
        explicits,
    )
}

/// Run an assembly source (`.fsim` s-expr or `fs-ir` JSON) through all six durable stages.
pub fn run_assembly(
    ledger: &Ledger,
    source: &[u8],
    packs: &CardPackSet,
    analysis: AssemblyAnalysisMode,
    clock: &mut dyn FnMut() -> f64,
    cancel_plan: AssemblyCancelPlan,
) -> Result<(SolveOutcome, Vec<String>), SolveRefusal> {
    let admitted = admit_assembly_source(source, packs)?;
    let run = assembly_run_id(&admitted, packs, analysis);
    let mut sess = open_session(ledger, clock, &admitted, packs, analysis, run)?;
    let state = SolveDriverState {
        run: *run.as_bytes(),
        project: *admitted.source_hash.as_bytes(),
        consumed_core_s: 0.0,
        consumed_wall_s: 0.0,
        completed: Vec::new(),
    };
    let mut pipe = AssemblyPipeline {
        admitted: &admitted,
        packs,
        analysis,
        run,
        gate: CancelGate::new_clock_free(),
        cancel_plan,
    };
    let res = sess.drive(&mut pipe, state, AssemblyCtx::default(), None);
    let mut progress = Vec::new();
    sess.drain_progress(&mut progress);
    res.map(|outcome| (outcome, progress))
}

/// Resume an interrupted assembly run from its longest attested checkpoint.
pub fn resume_assembly(
    ledger: &Ledger,
    run: SolveRunId,
    clock: &mut dyn FnMut() -> f64,
    cancel_plan: AssemblyCancelPlan,
) -> Result<(SolveOutcome, Vec<String>), SolveRefusal> {
    let (gate, inv) = (
        CancelGate::new_clock_free(),
        InvocationWorkLedger::default(),
    );
    let v = discover_and_attest_resume(
        ledger,
        run,
        EvidenceWork::new(&gate, None, &inv),
        &mut AssemblyResumeAttestor,
    )?;
    let mut sess = open_session(ledger, clock, &v.admitted, &v.packs, v.analysis, run)?;
    sess.restore_budget(&v.state, v.state_artifact)?;
    let mut pipe = AssemblyPipeline {
        admitted: &v.admitted,
        packs: &v.packs,
        analysis: v.analysis,
        run,
        gate: CancelGate::new_clock_free(),
        cancel_plan,
    };
    let res = sess.drive(&mut pipe, v.state, v.ctx, Some(v.state_artifact));
    let mut progress = Vec::new();
    sess.drain_progress(&mut progress);
    res.map(|outcome| (outcome, progress))
}

/// Project retained lowering receipt, execution plan, causalization audit, and
/// window energy ledger for `run` without replaying physics.
pub fn explain_assembly_run(
    ledger: &Ledger,
    run: SolveRunId,
) -> Result<AssemblyExplanation, SolveRefusal> {
    let (gate, inv) = (
        CancelGate::new_clock_free(),
        InvocationWorkLedger::default(),
    );
    let v = discover_and_attest_resume(
        ledger,
        run,
        EvidenceWork::new(&gate, None, &inv),
        &mut AssemblyResumeAttestor,
    )?;
    let (
        mut lowering_receipt_json,
        mut plan_json,
        mut causalization_audit_json,
        mut window_ledger_json,
    ) = (None, None, None, None);
    let mut stages = Vec::with_capacity(v.state.completed.len());
    for done in &v.state.completed {
        let stage = AssemblyStage::from_ordinal(done.ordinal)
            .ok_or_else(|| resume_identity("bad ordinal"))?;
        stages.push((stage.name().to_string(), done.op_id, done.receipt.to_hex()));
        match stage {
            AssemblyStage::Lower => {
                if let Some(h) = find_edge(
                    ledger,
                    done.op_id,
                    EdgeRole::Out,
                    Ledger::LOWERING_RECEIPT_ARTIFACT_KIND,
                )? {
                    lowering_receipt_json = Some(read_utf8(ledger, &h)?);
                }
            }
            AssemblyStage::Plan => {
                if let Some(h) = find_edge(
                    ledger,
                    done.op_id,
                    EdgeRole::Out,
                    Ledger::EXECUTION_PLAN_ARTIFACT_KIND,
                )? {
                    plan_json = Some(read_utf8(ledger, &h)?);
                }
                if let Some(h) = find_edge(
                    ledger,
                    done.op_id,
                    EdgeRole::Out,
                    Ledger::CAUSALIZATION_AUDIT_ARTIFACT_KIND,
                )? {
                    causalization_audit_json = Some(read_utf8(ledger, &h)?);
                }
            }
            AssemblyStage::Run => {
                if let Some(h) = find_edge(
                    ledger,
                    done.op_id,
                    EdgeRole::Out,
                    Ledger::WINDOW_LEDGER_ARTIFACT_KIND,
                )? {
                    window_ledger_json = Some(read_utf8(ledger, &h)?);
                }
            }
            _ => {}
        }
    }
    Ok(AssemblyExplanation {
        run: run.to_hex(),
        stages_completed: v.state.completed.len(),
        analysis: v.analysis.as_str().to_string(),
        lowering_receipt_hash: v.ctx.lowering_hash.map(|h| h.to_hex()),
        lowering_receipt_json,
        plan_hash: v.ctx.plan_hash,
        plan_json,
        causalization_audit_json,
        window_ledger_json,
        stages,
    })
}
