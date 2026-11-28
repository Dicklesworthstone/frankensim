//! Six-stage durable assembly pipeline (`admit -> lower -> plan -> run ->
//! probe -> report`) on `crate::stage_driver` plus `explain_assembly_run`
//! (bead `frankensim-rc-1109-root-3gc12.11.7`, Stage 2).

#![allow(clippy::result_large_err)]

mod budget;
pub(crate) mod cli;
mod report;

use budget::{
    BudgetMode, BudgetRun, EvaluatedQoiBudget, QOI_BUDGET_KIND, attach_wall_metrics,
    evaluate_qoi_budgets,
};

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
use fs_gas::component::ClosedGasParcel;
use fs_ledger::{EdgeRole, Ledger, LedgerError};
use fs_ports::component::PortSpec;
use fs_qty::Dims;
use fs_scenario::assembly::draft::{FidelityRung, PhysicsDomain, ProbeKind, StepPolicy};
use fs_system::acoustic::AcousticMicrophone;
use fs_system::elements::Motion;
use fs_system::lowering::Observed;
use fs_system::{
    AnalysisError, AnalysisOutcome, AnalysisRequest, LinearizeConfig, LoweredSystem,
    LoweringRegistry, ModalConfig, NodeId, OutputChannel, StateMap, StationaryConfig, System,
};

pub use fs_fsim::AssemblyAnalysisMode;
use fs_fsim::{AdmittedAssembly, PackCatalog, PackEntry};

fn admit_assembly_source(
    source: &[u8],
    packs: &CardPackSet,
) -> Result<AdmittedAssembly, SolveRefusal> {
    let to_entry = |p: &crate::cards::AdmittedCardPack| PackEntry {
        card: p.card(),
        identity: p.identity().to_owned(),
        bytes: p.bytes().to_vec(),
    };
    let catalog = PackCatalog::from_parts(
        packs.materials().iter().map(to_entry).collect(),
        packs.interfaces().iter().map(to_entry).collect(),
    );
    fs_fsim::admit_assembly_source(source, &catalog)
        .map_err(|r| staged_refusal(r.code, "admit", r.what, r.fix))
}
use self::report::{
    AssemblyReportInputs, EntitySample, EvaluatedProbe, PROBE_SERIES_KIND, RunTrace,
    build_assembly_report, dims_to_unit, evaluate_probes, plan_step_s,
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
    envelope_provenance, resolve_stage_consumption, staged_refusal, validate_checkpoint_prefix_for,
    validate_state_shape_for,
};

/// Wire version of the assembly stage driver.
pub const ASSEMBLY_DRIVER_VERSION: u32 = 2;
/// BLAKE3 domain separating assembly run identities from thermal solve runs.
pub const ASSEMBLY_RUN_IDENTITY_DOMAIN: &str = "org.frankensim.fs-cli.assembly-run-identity.v1";

const ASSEMBLY_STATE_TYPE_ID: u64 = 0x6673_636c_6961_736d; // b"fscliasm"
const ASSEMBLY_TRACE_KIND: &str = "assembly-run-trace";
const ASSEMBLY_ANALYSIS_META_KIND: &str = "assembly-analysis-mode";
const MAX_TRANSIENT_STEPS: usize = 200_000;

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

/// The registry every assembly lowers through: the kernel defaults, then the
/// owner groups (viscous networks, the thermal timeline, loss routes), which
/// lower only drafts that name their card families or caloric laws.
fn lowering_registry() -> LoweringRegistry {
    let mut registry = LoweringRegistry::with_kernel_defaults();
    fs_system::lowering::owners::register_owner_defaults(&mut registry);
    registry
}

/// Ordered six-stage assembly workflow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyStage {
    /// Parse, Five Explicits, card admission, and `DraftAssembly::admit`.
    Admit,
    /// Lower via `lowering_registry` with receipt verification.
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
    /// Canonical evaluated probe series JSON, when `probe` completed.
    pub probe_series_json: Option<String>,
    /// Canonical eight-term QoI uncertainty budget & verdict JSON, when requirements were declared.
    pub qoi_budget_json: Option<String>,
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

const GAS_VOL_DIMS: Dims = Dims([3, 0, 0, 0, 0, 0]);
const GAS_TEMP_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);

enum StateEntry {
    Motion(NodeId, String, Option<Dims>),
    Scalar(NodeId, String, Option<Dims>),
    Gas(NodeId, String, String),
}

struct LoweredStateMap {
    entries: Vec<StateEntry>,
}
impl LoweredStateMap {
    fn from_lowered(lowered: &LoweredSystem) -> Self {
        let mut entries = Vec::new();
        for (name, &node) in &lowered.node_map {
            let dims = lowered.system.node_state_dimensions(node);
            if lowered.system.state::<Motion>(node).is_some() {
                // `Motion::position_m` is an output integral absent from port laws;
                // mapping only `velocity_m_s` keeps the state Jacobian nonsingular
                // (matching `Mech2Map` / `ThermoMechMap` in `fs-system`).
                entries.push(StateEntry::Motion(
                    node,
                    format!("{name}.velocity_m_s"),
                    dims,
                ));
            } else if lowered.system.state::<f64>(node).is_some() {
                entries.push(StateEntry::Scalar(node, format!("{name}.state"), dims));
            } else if lowered.system.state::<ClosedGasParcel>(node).is_some() {
                entries.push(StateEntry::Gas(
                    node,
                    format!("{name}.volume_m3"),
                    format!("{name}.temperature_k"),
                ));
            }
        }
        Self { entries }
    }
    fn slot_meta(&self) -> Vec<(String, Option<Dims>)> {
        let mut out = Vec::with_capacity(self.dim());
        for e in &self.entries {
            match e {
                StateEntry::Motion(_, l, d) | StateEntry::Scalar(_, l, d) => {
                    out.push((l.clone(), *d));
                }
                StateEntry::Gas(_, lv, lt) => {
                    out.push((lv.clone(), Some(GAS_VOL_DIMS)));
                    out.push((lt.clone(), Some(GAS_TEMP_DIMS)));
                }
            }
        }
        out
    }
    fn dims(&self, i: usize) -> Option<Dims> {
        self.slot_meta().get(i).and_then(|c| c.1)
    }
}
impl StateMap for LoweredStateMap {
    fn dim(&self) -> usize {
        self.entries
            .iter()
            .map(|e| usize::from(matches!(e, StateEntry::Gas(..))) + 1)
            .sum()
    }
    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        let inv = |what| AnalysisError::Invalid { what };
        if out.len() != self.dim() {
            return Err(inv("state length mismatch"));
        }
        let mut i = 0;
        for e in &self.entries {
            match e {
                StateEntry::Motion(n, ..) => {
                    let m = system.state::<Motion>(*n).ok_or(inv("missing Motion"))?;
                    out[i] = m.velocity_m_s;
                    i += 1;
                }
                StateEntry::Scalar(n, ..) => {
                    out[i] = *system.state::<f64>(*n).ok_or(inv("missing f64"))?;
                    i += 1;
                }
                StateEntry::Gas(n, ..) => {
                    let g = system
                        .state::<ClosedGasParcel>(*n)
                        .ok_or(inv("missing gas"))?;
                    out[i] = g.volume_m3();
                    out[i + 1] = g.state().state().temperature;
                    i += 2;
                }
            }
        }
        Ok(())
    }
    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError> {
        let inv = |what| AnalysisError::Invalid { what };
        if state.len() != self.dim() {
            return Err(inv("state length mismatch"));
        }
        let mut i = 0;
        for e in &self.entries {
            match e {
                StateEntry::Motion(n, ..) => {
                    let mut m = *system.state::<Motion>(*n).ok_or(inv("missing Motion"))?;
                    m.velocity_m_s = state[i];
                    system.set_accepted_state(*n, m)?;
                    i += 1;
                }
                StateEntry::Scalar(n, ..) => {
                    system.set_accepted_state(*n, state[i])?;
                    i += 1;
                }
                StateEntry::Gas(n, ..) => {
                    let g = system
                        .state::<ClosedGasParcel>(*n)
                        .ok_or(inv("missing gas"))?;
                    let next = g
                        .at_volume_temperature(state[i], state[i + 1])
                        .map_err(|_| inv("gas"))?;
                    system.set_accepted_state(*n, next)?;
                    i += 2;
                }
            }
        }
        Ok(())
    }
    fn labels(&self) -> Vec<String> {
        self.slot_meta().into_iter().map(|(l, _)| l).collect()
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
    largest_energy_defect_j: Option<f64>,
    trace: Option<RunTrace>,
    probe_series_artifact: Option<ContentHash>,
    probes: Vec<EvaluatedProbe>,
    qoi_budget: Option<(ContentHash, EvaluatedQoiBudget)>,
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
        largest_energy_defect_j: Option<f64>,
    },
    Probe {
        series_artifact: ContentHash,
        probes: Vec<EvaluatedProbe>,
        qoi_budget: Option<(ContentHash, EvaluatedQoiBudget)>,
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
        self.run_transient_with(lowered, None)
    }

    #[allow(clippy::too_many_lines)]
    fn run_transient_with(
        &self,
        lowered: &mut LoweredSystem,
        dt_override: Option<f64>,
    ) -> Result<(RunTrace, fs_system::LedgerExport), SolveRefusal> {
        let tp = &self.admitted.draft.input().time_plan;
        let dt_s = match dt_override {
            Some(dt) => dt,
            None => match tp.step_policy {
                StepPolicy::Fixed { dt } => dt.value,
                StepPolicy::Adaptive { .. } => {
                    return Err(staged_refusal(
                        "cli-assembly-step-policy",
                        "run",
                        "the assembly runner does not yet bind an adaptive state-error metric",
                        "use the fixed step policy",
                    ));
                }
            },
        };
        let end_s = lowered.system.time_s() + tp.horizon.value;
        if !(end_s.is_finite() && end_s > lowered.system.time_s()) {
            return Err(staged_refusal(
                "cli-assembly-time-range",
                "run",
                "the requested horizon has no finite advancing endpoint",
                "choose a representable simulation duration",
            ));
        }
        let step_limit = || {
            staged_refusal(
                "cli-assembly-step-limit",
                "run",
                format!(
                    "the full {} s horizon requires more than {MAX_TRANSIENT_STEPS} accepted steps",
                    tp.horizon.value,
                ),
                "increase the fixed step size or reduce the horizon",
            )
        };
        // A count estimated from the requested step is only a capacity hint:
        // scheduled events and owners may shorten any accepted kernel step.
        let estimated_steps = (tp.horizon.value / dt_s).ceil();
        // Leave one step of slack for division rounding at the limit. The
        // accepted-step guard below is authoritative.
        if estimated_steps > (MAX_TRANSIENT_STEPS + 1) as f64 {
            return Err(step_limit());
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let steps = (estimated_steps as usize).min(MAX_TRANSIENT_STEPS);
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
                let mut s = EntitySample::default();
                // A part in a slot of an owner node (a body in a network, a cell in a timeline).
                let observer = lowered.observers.get(&part.entity);
                match observer.and_then(|read| read(sys)) {
                    Some(Observed::Body {
                        velocity_m_s: [vx, vy, vz],
                        angular_velocity_rad_s: [wx, wy, wz],
                    }) => s.motion = Some([vx, vy, vz, wx, wy, wz]),
                    Some(Observed::Cell(t)) => s.temperature_k = Some(t),
                    Some(Observed::HeatRate(q)) => s.heat_rate_w = Some(q),
                    None => {}
                }
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
                    } else if let Some(g) = sys.state::<ClosedGasParcel>(node) {
                        s.temperature_k = Some(g.state().state().temperature);
                    }
                }
                if let Some(vec) = ents.get_mut(&part.name) {
                    vec.push(s);
                }
            }
        };
        let mut mics: Vec<(String, AcousticMicrophone, Vec<f64>)> = Vec::new();
        for p in &self.admitted.draft.input().probes {
            if matches!(p.kind, ProbeKind::Microphone { .. }) {
                let mic = AcousticMicrophone::from_draft_and_lowered(
                    p,
                    &self.admitted.draft,
                    lowered,
                    dt_s,
                )
                .map_err(|why| {
                    staged_refusal(
                        "cli-assembly-probe-unsupported",
                        "run",
                        format!("microphone `{}` refused: {why}", p.name),
                        "check microphone target, position, decimation ratio, and Nyquist limits",
                    )
                })?;
                let cap = steps / mic.decimation_ratio() + 1;
                mics.push((p.name.clone(), mic, Vec::with_capacity(cap)));
            }
        }
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
            while lowered.system.time_s() < end_s {
                if self.gate.is_requested() {
                    return Err(SolveRefusal::plain(
                        "cli-solve-cancelled",
                        "cancelled during stepping",
                        "resume",
                    ));
                }
                if times_s.len() - 1 == MAX_TRANSIENT_STEPS {
                    return Err(step_limit());
                }
                let before_s = lowered.system.time_s();
                lowered.system.step_toward(end_s, dt_s, &cx).map_err(|e| {
                    staged_refusal(
                        "cli-assembly-step",
                        "run",
                        format!("kernel step failed: {e}"),
                        "check step size",
                    )
                })?;
                let after_s = lowered.system.time_s();
                if !(after_s.is_finite() && after_s > before_s) {
                    return Err(staged_refusal(
                        "cli-assembly-time-progress",
                        "run",
                        "the accepted kernel step made no finite clock progress",
                        "increase the step size or rescale the simulation time",
                    ));
                }
                if !mics.is_empty() {
                    let step_dt = after_s - before_s;
                    if (step_dt - dt_s).abs() <= 1.0e-6 * dt_s {
                        for (name, mic, series) in &mut mics {
                            if let Some(pa) = mic.advance(&lowered.system).map_err(|why| {
                                staged_refusal(
                                    "cli-assembly-step",
                                    "run",
                                    format!("microphone `{name}` failed: {why}"),
                                    "check modal stability and microphone bounds",
                                )
                            })? {
                                series.push(pa);
                            }
                        }
                    } else if after_s < end_s - 1.0e-6 * dt_s {
                        return Err(staged_refusal(
                            "cli-assembly-probe-unsupported",
                            "run",
                            "microphone requires uniform fixed steps without mid-run shortening",
                            "use a uniform fixed step plan for microphone probes",
                        ));
                    }
                }
                times_s.push(after_s);
                sample(&lowered.system, &mut entities);
            }
            Ok(())
        })?;
        let microphones = mics.into_iter().map(|(n, _, s)| (n, s)).collect();
        Ok((
            RunTrace {
                analysis: "transient".into(),
                times_s,
                entities,
                microphones,
                modal_modes: Vec::new(),
                harmonic_points: Vec::new(),
                harmonic_unit: "SI".into(),
                stationary_norm: None,
                stationary_state: Vec::new(),
            },
            lowered.system.ledger().export(),
        ))
    }

    #[allow(clippy::too_many_lines)]
    fn run_non_transient_with(
        &self,
        lowered: &mut LoweredSystem,
        dt_scale: f64,
        tight: bool,
    ) -> Result<RunTrace, SolveRefusal> {
        let smap = LoweredStateMap::from_lowered(lowered);
        let mut harmonic_unit = "SI".to_string();
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
                AssemblyAnalysisMode::Stationary => {
                    let mut cfg = StationaryConfig {
                        commit: true,
                        ..StationaryConfig::default()
                    };
                    cfg.probe_dt_s *= dt_scale;
                    if tight {
                        cfg.absolute_tolerance *= 1.0e-2;
                        cfg.monolithic.absolute_tolerance *= 1.0e-1;
                    }
                    lowered.system.run_analysis(
                        AnalysisRequest::Stationary {
                            state_map: &smap,
                            config: cfg,
                        },
                        &cx,
                    )
                }
                AssemblyAnalysisMode::Modal => {
                    let mut cfg = ModalConfig::default();
                    cfg.linearize.dt_s *= dt_scale;
                    if tight {
                        cfg.linearize.fd_eps *= 0.1;
                    }
                    lowered.system.run_analysis(
                        AnalysisRequest::Modal {
                            state_map: &smap,
                            config: cfg,
                        },
                        &cx,
                    )
                }
                AssemblyAnalysisMode::Harmonic => {
                    let inputs = lowered.system.boundary_channels();
                    if inputs.is_empty() {
                        return Err(AnalysisError::Invalid {
                            what: "no boundary channels available for harmonic response",
                        });
                    }
                    let in_d = lowered
                        .system
                        .port_spec(inputs[0].port)
                        .and_then(PortSpec::input_dimensions);
                    if let Some(u) = smap
                        .dims(0)
                        .zip(in_d)
                        .and_then(|(o, i)| o.checked_minus(i))
                        .and_then(dims_to_unit)
                    {
                        harmonic_unit = u.into();
                    }
                    let outputs: Vec<OutputChannel> =
                        (0..smap.dim()).map(OutputChannel::State).collect();
                    let hz = if self.admitted.harmonic_frequencies_hz.is_empty() {
                        vec![1.0, 5.0, 10.0]
                    } else {
                        self.admitted.harmonic_frequencies_hz.clone()
                    };
                    let omegas: Vec<f64> = hz.iter().map(|f| f * TAU).collect();
                    let mut lcfg = LinearizeConfig::default();
                    lcfg.dt_s *= dt_scale;
                    if tight {
                        lcfg.fd_eps *= 0.1;
                    }
                    lowered.system.run_analysis(
                        AnalysisRequest::Harmonic {
                            state_map: &smap,
                            inputs: &inputs,
                            outputs: &outputs,
                            frequencies_rad_s: &omegas,
                            config: lcfg,
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
            microphones: BTreeMap::new(),
            modal_modes: Vec::new(),
            harmonic_points: Vec::new(),
            harmonic_unit,
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
                trace.stationary_state = smap
                    .slot_meta()
                    .into_iter()
                    .zip(s.state)
                    .map(|((label, dims), v)| {
                        (label, v, dims.and_then(dims_to_unit).unwrap_or("SI").into())
                    })
                    .collect();
                for part in &self.admitted.draft.input().parts {
                    if let Some(Observed::HeatRate(q)) = lowered
                        .observers
                        .get(&part.entity)
                        .and_then(|read| read(&lowered.system))
                    {
                        let s = EntitySample {
                            heat_rate_w: Some(q),
                            ..Default::default()
                        };
                        trace.entities.insert(part.name.clone(), vec![s]);
                    }
                }
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
                if self.analysis != AssemblyAnalysisMode::Stationary
                    && let Some(p) = inp.parts.iter().find(|p| {
                        p.fidelity
                            .rungs
                            .iter()
                            .any(|r| matches!(r, FidelityRung::OneDField { .. }))
                    })
                {
                    let m = self.analysis.as_str();
                    return Err(staged_refusal(
                        "cli-assembly-fidelity-analysis",
                        "admit",
                        format!(
                            "part `{}` (one-d-field 1) requires --analysis stationary (got `{m}`)",
                            p.name
                        ),
                        "re-run with --analysis stationary",
                    ));
                }
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
                let lowered = lowering_registry()
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
                let (trace, sides, ledger_artifact, energy_residual_j, largest_energy_defect_j) =
                    if self.analysis == AssemblyAnalysisMode::Transient {
                        let (trace, export) = self.run_transient(lowered)?;
                        let defect = lowered.system.ledger().energy_residual_j();
                        let largest = lowered.system.ledger().largest_absolute_energy_defect_j();
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
                            Some(largest),
                        )
                    } else {
                        let trace = self.run_non_transient_with(lowered, 1.0, false)?;
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
                        largest_energy_defect_j,
                    },
                })
            }
            AssemblyStage::Probe => {
                let trace = ctx
                    .trace
                    .as_ref()
                    .ok_or_else(|| resume_identity("missing run trace"))?;
                let input = self.admitted.draft.input();
                let (receipt_json, side, probes) = evaluate_probes(
                    &input.probes,
                    &self.admitted.entity_names,
                    trace,
                    plan_step_s(&input.time_plan),
                )?;
                let series_artifact = side.artifact;
                let mut side_artifacts = vec![side];
                let active_reqs: Vec<_> = self
                    .admitted
                    .requirements
                    .iter()
                    .filter(|r| match r.analysis {
                        Some(m) => m == self.analysis,
                        None => {
                            self.analysis == AssemblyAnalysisMode::Transient
                                || probes.iter().any(|p| p.name == r.probe)
                        }
                    })
                    .collect();
                let qoi_budget = if active_reqs.is_empty() {
                    None
                } else {
                    use AssemblyAnalysisMode as Am;
                    use BudgetRun as Br;
                    let defect = ctx
                        .energy_residual_j
                        .map(f64::abs)
                        .into_iter()
                        .chain(ctx.largest_energy_defect_j.map(f64::abs))
                        .reduce(f64::max);
                    let mode = match self.analysis {
                        Am::Transient => BudgetMode::Transient {
                            energy_defect_j: defect,
                        },
                        Am::Stationary => BudgetMode::Stationary {
                            residual_norm: trace.stationary_norm,
                        },
                        Am::Modal => BudgetMode::Modal,
                        Am::Harmonic => BudgetMode::Harmonic,
                    };
                    let (qb, qs, bm) = evaluate_qoi_budgets(
                        self.admitted,
                        &active_reqs,
                        &probes,
                        ctx.lowered.as_ref().map(|l| &l.receipt),
                        mode,
                        state.consumed_wall_s,
                        |ls, run| match (self.analysis, run) {
                            (Am::Transient, Br::HalfStep(d)) => {
                                self.run_transient_with(ls, d).map(|(t, _)| t)
                            }
                            (Am::Transient, _) => self.run_transient_with(ls, None).map(|(t, _)| t),
                            (_, Br::HalfStep(_)) => self.run_non_transient_with(ls, 0.5, false),
                            (_, Br::TightSolver) => self.run_non_transient_with(ls, 1.0, true),
                            (_, Br::Nominal) => self.run_non_transient_with(ls, 1.0, false),
                        },
                    )?;
                    let qa = qs.artifact;
                    side_artifacts.push(qs);
                    ctx.pending_metrics.extend(bm);
                    Some((qa, qb))
                };
                Ok(StageExecution {
                    receipt_json,
                    side_artifacts,
                    on_commit: StageCommit::Probe {
                        series_artifact,
                        probes,
                        qoi_budget,
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
                        qoi_budget: ctx.qoi_budget.as_ref().map(|(h, b)| (*h, b)),
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
                largest_energy_defect_j,
            } => {
                ctx.trace = Some(trace);
                ctx.ledger_artifact = ledger_artifact;
                ctx.energy_residual_j = energy_residual_j;
                ctx.largest_energy_defect_j = largest_energy_defect_j;
                AssemblyStage::Run
            }
            StageCommit::Probe {
                series_artifact,
                probes,
                qoi_budget,
            } => {
                ctx.probe_series_artifact = Some(series_artifact);
                ctx.probes = probes;
                ctx.qoi_budget = qoi_budget;
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
        mut state: SolveDriverState,
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

        let (mut ctx, mut prev, mut resolved_consumption) = (AssemblyCtx::default(), None, None);
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
            let stage_consumption = resolve_stage_consumption(
                ledger,
                done.op_id,
                done.ordinal,
                &cp,
                resolved_consumption,
                ASSEMBLY_SCHEMA.cores,
            )?;
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
                    if state.completed.len() < 4
                        || (state.completed.len() == 4
                            && analysis == AssemblyAnalysisMode::Transient
                            && !admitted.requirements.is_empty())
                    {
                        let lowered = lowering_registry()
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
                    if let Some(la) = ctx.ledger_artifact
                        && let Ok(v) = JsonValue::parse(&read_utf8(ledger, &la)?)
                        && let Some(t) = v.get("totals")
                    {
                        ctx.energy_residual_j = t.f64_field("energy_residual_j");
                        ctx.largest_energy_defect_j =
                            t.f64_field("largest_absolute_energy_defect_j");
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
                    if let Some(qa) = find_edge(ledger, done.op_id, EdgeRole::Out, QOI_BUDGET_KIND)?
                    {
                        ctx.qoi_budget =
                            Some((qa, EvaluatedQoiBudget::from_json(&read_utf8(ledger, &qa)?)?));
                    }
                }
                AssemblyStage::Admit | AssemblyStage::Report => {}
            }
            prev = Some(cp);
            resolved_consumption = Some(stage_consumption);
        }
        if let Some((consumed_core_s, consumed_wall_s)) = resolved_consumption {
            state.consumed_core_s = consumed_core_s;
            state.consumed_wall_s = consumed_wall_s;
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
        mut probe_series_json,
        mut qoi_budget_json,
    ) = (None, None, None, None, None, None);
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
            AssemblyStage::Probe => {
                if let Some(h) = find_edge(ledger, done.op_id, EdgeRole::Out, PROBE_SERIES_KIND)? {
                    probe_series_json = Some(read_utf8(ledger, &h)?);
                }
                if let Some(h) = find_edge(ledger, done.op_id, EdgeRole::Out, QOI_BUDGET_KIND)? {
                    let raw = read_utf8(ledger, &h)?;
                    let m = ledger.metrics_for_op(done.op_id).unwrap_or_default();
                    qoi_budget_json = Some(attach_wall_metrics(&raw, &m));
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
        probe_series_json,
        qoi_budget_json,
        stages,
    })
}
