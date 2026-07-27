//! Deterministic solve orchestration for the `solve` product verb
//! (bead frankensim-extreal-program-f85xj.6.5, slice 1).
//!
//! The driver composes surfaces that already exist separately — fs-session
//! budgets and metering, fs-exec cancellation and the versioned legacy
//! snapshot envelope, fs-ledger operations and content-addressed artifacts —
//! into one product path. Every stage is a ledgered operation carrying the
//! Five Explicits; the driver state after each completed stage is sealed
//! through the fs-exec legacy v1 snapshot envelope and retained as a ledger
//! artifact, so a killed or cancelled run resumes from durable evidence
//! alone. Stages whose authoritative producers are still open refuse with
//! `cli-solve-stage-gap` naming the owning bead; a skeleton run is never
//! presented as an integrated solve.
//!
//! Slice-1 boundary (stated, not implied): `import-verify` and `assign`
//! execute against retained import evidence; `material-resolve`
//! (frankensim-hp7tb), `flow-network` (frankensim-frn2i), `conduction`
//! (frankensim-s93ej), and `qoi` (frankensim-s2l9v) are typed gaps.

// Refusals are cold-path values carrying complete diagnostics; the crate's
// refusal idiom (`GeometryImportRefusal`) is by-value for the same reason.
#![allow(clippy::result_large_err)]

use std::fmt::Write as _;

use fs_blake3::hash_domain;
use fs_blake3::identity::ContentId;
use fs_exec::CancelGate;
use fs_exec::solver::{
    LegacySnapshotExpectationV1, LegacySnapshotLimitsV1, LegacySnapshotV1Adapter,
    LegacySolverStateV1, codec,
};
use fs_ledger::{ContentHash, EdgeRole, FiveExplicits, Ledger, LedgerError, OpOutcome, hash_bytes};
use fs_project::{DecodedProject, ProjectSpec, geometry_source_identity};
use fs_session::{CapabilityToken, Charge, Enforcement, Governor, SessionError, SessionId};

use crate::import::{explicits, json_string};

/// Domain separating solve-run identity derivation from every other hash.
pub const SOLVE_RUN_IDENTITY_DOMAIN: &str = "org.frankensim.fs-cli.solve-run.v1";
/// Driver semantics version bound into run identity and driver state.
pub const SOLVE_DRIVER_VERSION: u32 = 1;

const SOLVE_STAGE_SCHEMA: &str = "frankensim.cli.solve-stage.v1";
const SOLVE_RUN_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-run-receipt.v1";
const IMPORT_VERIFY_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-import-verify-receipt.v1";
const ASSIGN_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-assignment-binding.v1";

const PROJECT_SOURCE_KIND: &str = "solve-project-source";
const STAGE_STATE_KIND: &str = "solve-stage-state";
const STAGE_RECEIPT_KIND: &str = "solve-stage-receipt";
const RUN_RECEIPT_KIND: &str = "solve-run-receipt";
const IMPORT_SUMMARY_KIND: &str = "geometry-import-run-receipt";

/// Historical type identity of the driver state inside the legacy v1
/// snapshot envelope (`b"fsclisol"` as big-endian bytes).
const DRIVER_STATE_TYPE_ID_V1: u64 = 0x6673_636c_6973_6f6c;
const DRIVER_STATE_SCHEMA_VERSION_V1: u32 = 1;

/// Whole-envelope cap for a retained driver-state snapshot.
const MAX_STATE_ENVELOPE_BYTES: u64 = 4 * 1024 * 1024;
const STATE_HASH_POLL_BYTES: u32 = 64 * 1024;

/// Read cap for stage receipts and import summaries consumed on resume.
const MAX_RECEIPT_READ_BYTES: u64 = 4 * 1024 * 1024;
/// Read cap when re-hashing retained import artifacts (matches the import
/// per-source cap).
const MAX_EVIDENCE_READ_BYTES: u64 = 64 * 1024 * 1024;
/// Edge scan cap per operation while locating retained evidence.
const EDGE_SCAN_CAP: usize = 512;

/// Warn when consumed wall budget crosses these fractions of the grant.
const BUDGET_WARN_FRACTIONS: [f64; 2] = [0.5, 0.9];

/// Cores leased by the slice-1 driver. The pipeline is single-threaded and
/// deterministic; the core-second grant is derived as `wall * CORES`.
const SOLVE_CORES: u64 = 1;

/// Canonical stage order. The ordinal is the array index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveStage {
    /// Verify retained import evidence against the run's pinned project.
    ImportVerify,
    /// Bind verified assignment evidence into the run's lineage.
    Assign,
    /// Resolve material/interface cards (gap: frankensim-hp7tb).
    MaterialResolve,
    /// Solve the enclosure flow network (gap: frankensim-frn2i).
    FlowNetwork,
    /// Conduction/conjugate coupling solve (gap: frankensim-s93ej).
    Conduction,
    /// QoI extraction against requirements (gap: frankensim-s2l9v).
    Qoi,
}

impl SolveStage {
    /// All stages in execution order.
    pub const ALL: [SolveStage; 6] = [
        SolveStage::ImportVerify,
        SolveStage::Assign,
        SolveStage::MaterialResolve,
        SolveStage::FlowNetwork,
        SolveStage::Conduction,
        SolveStage::Qoi,
    ];

    /// Stable kebab-case stage name used in IR, receipts, and diagnostics.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            SolveStage::ImportVerify => "import-verify",
            SolveStage::Assign => "assign",
            SolveStage::MaterialResolve => "material-resolve",
            SolveStage::FlowNetwork => "flow-network",
            SolveStage::Conduction => "conduction",
            SolveStage::Qoi => "qoi",
        }
    }

    /// Capability verb the session token must grant for this stage.
    #[must_use]
    pub const fn verb(self) -> &'static str {
        match self {
            SolveStage::ImportVerify => "solve.import-verify",
            SolveStage::Assign => "solve.assign",
            SolveStage::MaterialResolve => "solve.material-resolve",
            SolveStage::FlowNetwork => "solve.flow-network",
            SolveStage::Conduction => "solve.conduction",
            SolveStage::Qoi => "solve.qoi",
        }
    }

    /// The Beads id that owns this stage's authoritative producer when the
    /// stage is still a typed gap, or `None` when the stage executes.
    #[must_use]
    pub const fn gap_dependency(self) -> Option<&'static str> {
        match self {
            SolveStage::ImportVerify | SolveStage::Assign => None,
            SolveStage::MaterialResolve => Some("frankensim-hp7tb"),
            SolveStage::FlowNetwork => Some("frankensim-frn2i"),
            SolveStage::Conduction => Some("frankensim-s93ej"),
            SolveStage::Qoi => Some("frankensim-s2l9v"),
        }
    }

    fn ordinal(self) -> u32 {
        match self {
            SolveStage::ImportVerify => 0,
            SolveStage::Assign => 1,
            SolveStage::MaterialResolve => 2,
            SolveStage::FlowNetwork => 3,
            SolveStage::Conduction => 4,
            SolveStage::Qoi => 5,
        }
    }

    fn from_ordinal(ordinal: u32) -> Option<SolveStage> {
        SolveStage::ALL.get(ordinal as usize).copied()
    }
}

/// Content-derived run identity: the hash of the exact inputs that determine
/// the run's answer plus the driver semantics version. Budgets travel inside
/// the project hash, so raising a budget starts a new run whose completed
/// artifacts still deduplicate by content.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SolveRunId(ContentHash);

impl SolveRunId {
    /// Derive the run identity from the decoded project.
    #[must_use]
    pub fn derive(project: &DecodedProject) -> SolveRunId {
        let project_hash = project.hash();
        let spec = &project.spec;
        let mut preimage = Vec::with_capacity(128);
        preimage.extend_from_slice(project_hash.as_bytes());
        let (constellation, workspace) = spec.versions.as_ref().map_or(("", ""), |v| {
            (v.constellation.as_str(), v.workspace.as_str())
        });
        push_framed(&mut preimage, constellation.as_bytes());
        push_framed(&mut preimage, workspace.as_bytes());
        let seed = spec.seeds.as_ref().map_or(0, |s| s.root);
        preimage.extend_from_slice(&seed.to_le_bytes());
        preimage.extend_from_slice(&SOLVE_DRIVER_VERSION.to_le_bytes());
        SolveRunId(hash_domain(SOLVE_RUN_IDENTITY_DOMAIN, &preimage))
    }

    /// Parse the 64-hex user-facing run id.
    #[must_use]
    pub fn parse_hex(value: &str) -> Option<SolveRunId> {
        ContentHash::from_hex(value).map(SolveRunId)
    }

    /// Lowercase hexadecimal rendering (the user-facing run id).
    #[must_use]
    pub fn to_hex(self) -> String {
        self.0.to_hex()
    }

    /// Exact identity bytes; also the `session` value on every run op.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        self.0.as_bytes()
    }

    /// Session identity for the governor: first eight identity bytes.
    fn session_u64(self) -> u64 {
        let bytes = self.0.as_bytes();
        u64::from_le_bytes([
            bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
        ])
    }
}

/// Derive the session capability token from project budgets. The project
/// declares wall-clock seconds, peak memory bytes, and relative accuracy;
/// core-seconds and the core ceiling are driver-derived because the slice-1
/// pipeline is single-threaded.
pub(crate) fn derive_capability_token(
    spec: &ProjectSpec,
    run: SolveRunId,
) -> Result<CapabilityToken, SolveRefusal> {
    let budgets = spec.budgets.as_ref().ok_or_else(|| {
        SolveRefusal::plain(
            "project-budgets-missing",
            "project budgets are unavailable at solve",
            "run strict project validation before solve",
        )
    })?;
    let wall_s = budgets.solve_time.value;
    if !wall_s.is_finite() || wall_s <= 0.0 {
        return Err(SolveRefusal::plain(
            "cli-solve-budget",
            format!("project solve-time budget `{wall_s}` is not a positive finite second count"),
            "declare a positive finite `:solve-time` budget",
        ));
    }
    let mut ops: Vec<String> = SolveStage::ALL
        .iter()
        .map(|stage| stage.verb().to_string())
        .collect();
    ops.push("solve.terminal".to_string());
    #[allow(clippy::cast_precision_loss)]
    let core_s = wall_s * SOLVE_CORES as f64;
    let token = CapabilityToken {
        session: SessionId(run.session_u64()),
        ops,
        core_s,
        mem_bytes: budgets.memory_bytes,
        wall_s,
        cores: SOLVE_CORES,
        ledger_scope: format!("solve-{}", &run.to_hex()[..16]),
    };
    Ok(token)
}

/// One completed stage retained inside the driver state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CompletedStage {
    /// Stage ordinal in [`SolveStage::ALL`] order.
    pub ordinal: u32,
    /// Ledger operation id of the stage.
    pub op_id: i64,
    /// Content hash of the stage receipt artifact.
    pub receipt: ContentHash,
}

/// Driver state sealed into the legacy v1 snapshot envelope after each
/// completed stage. Contains no wall-clock timestamps: consumption totals
/// are budget accounting, not identity.
#[derive(Debug, Clone, PartialEq)]
pub struct SolveDriverState {
    /// Run identity bytes.
    pub run: [u8; 32],
    /// Pinned project canonical hash.
    pub project: [u8; 32],
    /// Core-seconds charged so far.
    pub consumed_core_s: f64,
    /// Wall seconds charged so far.
    pub consumed_wall_s: f64,
    /// Completed stages in order.
    pub completed: Vec<CompletedStage>,
}

impl LegacySolverStateV1 for SolveDriverState {
    const SCHEMA_VERSION_V1: u32 = DRIVER_STATE_SCHEMA_VERSION_V1;
    const TYPE_ID_V1: u64 = DRIVER_STATE_TYPE_ID_V1;

    fn encode_v1(&self, enc: &mut codec::Enc) {
        for chunk in self.run.as_chunks::<8>().0 {
            enc.put_u64(u64::from_le_bytes(*chunk));
        }
        for chunk in self.project.as_chunks::<8>().0 {
            enc.put_u64(u64::from_le_bytes(*chunk));
        }
        enc.put_f64(self.consumed_core_s);
        enc.put_f64(self.consumed_wall_s);
        enc.put_u32(u32::try_from(self.completed.len()).expect("stage count fits u32"));
        for stage in &self.completed {
            enc.put_u32(stage.ordinal);
            enc.put_u64(stage.op_id.cast_unsigned());
            for chunk in stage.receipt.as_bytes().as_chunks::<8>().0 {
                enc.put_u64(u64::from_le_bytes(*chunk));
            }
        }
    }

    fn decode_v1(dec: &mut codec::Dec<'_>) -> Result<Self, codec::CodecError> {
        let mut run = [0u8; 32];
        for chunk in run.as_chunks_mut::<8>().0 {
            *chunk = dec.get_u64()?.to_le_bytes();
        }
        let mut project = [0u8; 32];
        for chunk in project.as_chunks_mut::<8>().0 {
            *chunk = dec.get_u64()?.to_le_bytes();
        }
        let consumed_core_s = dec.get_f64()?;
        let consumed_wall_s = dec.get_f64()?;
        let count = dec.get_u32()?;
        let mut completed = Vec::new();
        for _ in 0..count {
            let ordinal = dec.get_u32()?;
            let op_id = dec.get_u64()?.cast_signed();
            let mut receipt = [0u8; 32];
            for chunk in receipt.as_chunks_mut::<8>().0 {
                *chunk = dec.get_u64()?.to_le_bytes();
            }
            completed.push(CompletedStage {
                ordinal,
                op_id,
                receipt: ContentHash::from_slice(&receipt).expect("32-byte receipt hash"),
            });
        }
        Ok(SolveDriverState {
            run,
            project,
            consumed_core_s,
            consumed_wall_s,
            completed,
        })
    }
}

/// Structured solve refusal mirroring the import refusal shape.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SolveRefusal {
    /// Stable machine code.
    pub code: &'static str,
    /// Stage under execution when refusing, if any.
    pub stage: Option<&'static str>,
    /// What was refused.
    pub what: String,
    /// Actionable fix.
    pub fix: String,
    /// Owning bead when the refusal is an unbuilt-producer gap.
    pub dependency: Option<&'static str>,
    /// Run identity when derivation succeeded before the refusal.
    pub run: Option<String>,
    /// Ledger op retaining the refusal, when recording succeeded.
    pub recorded_op: Option<i64>,
}

impl SolveRefusal {
    fn plain(code: &'static str, what: impl Into<String>, fix: impl Into<String>) -> SolveRefusal {
        SolveRefusal {
            code,
            stage: None,
            what: what.into(),
            fix: fix.into(),
            dependency: None,
            run: None,
            recorded_op: None,
        }
    }

    fn staged(
        code: &'static str,
        stage: SolveStage,
        what: impl Into<String>,
        fix: impl Into<String>,
    ) -> SolveRefusal {
        SolveRefusal {
            code,
            stage: Some(stage.name()),
            what: what.into(),
            fix: fix.into(),
            dependency: None,
            run: None,
            recorded_op: None,
        }
    }
}

impl core::fmt::Display for SolveRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}: {}", self.code, self.what)
    }
}

impl std::error::Error for SolveRefusal {}

/// Why a run stopped without refusing.
#[derive(Debug, Clone, PartialEq)]
pub enum SolveRunStatus {
    /// Every stage completed.
    Completed,
    /// Budget enforcement stopped the run after a completed stage; the
    /// completed prefix is durable and the run record names the resume path.
    BudgetExceeded {
        /// Exhausted resource name from the governor.
        resource: &'static str,
        /// Consumed amount in the resource's unit.
        used: f64,
        /// Granted amount in the resource's unit.
        granted: f64,
    },
    /// Cancellation was observed between stages; the completed prefix is
    /// durable and resumable.
    Cancelled,
}

/// One completed stage summary for rendering.
#[derive(Debug, Clone, PartialEq)]
pub struct StageSummary {
    /// Stage name.
    pub stage: &'static str,
    /// Ledger op id.
    pub op_id: i64,
    /// Stage receipt artifact hash (hex).
    pub receipt: String,
    /// Measured wall seconds for the stage (reporting only, not identity).
    pub wall_s: f64,
}

/// Outcome of a solve or resume invocation that did not refuse.
#[derive(Debug, Clone, PartialEq)]
pub struct SolveOutcome {
    /// User-facing run id (64 hex).
    pub run: String,
    /// Terminal status.
    pub status: SolveRunStatus,
    /// Stages completed by THIS invocation, in order.
    pub stages: Vec<StageSummary>,
    /// Stages completed before this invocation (resume prefix length).
    pub prior_stages: u32,
    /// Terminal run-receipt artifact hash (hex), when one was written.
    pub run_receipt: Option<String>,
}

/// In-memory context threaded between stages; reconstructible from retained
/// stage receipts on resume.
#[derive(Debug, Default)]
struct StageContext {
    /// Verified imports: (role, source identity, promoted mesh, report).
    verified_imports: Vec<VerifiedImport>,
}

#[derive(Debug, Clone)]
struct VerifiedImport {
    role: String,
    source_identity: String,
    raw_source: ContentHash,
    promotion_receipt: ContentHash,
    promoted_mesh: ContentHash,
    assignment_report: ContentHash,
}

/// Everything a running or resumed solve needs in one place.
struct SolveEngine<'a> {
    ledger: &'a Ledger,
    gate: &'a CancelGate,
    clock: &'a mut dyn FnMut() -> f64,
    spec: &'a ProjectSpec,
    canonical_source: &'a str,
    project_hash: ContentHash,
    run: SolveRunId,
    governor: Governor,
    session: SessionId,
    token: CapabilityToken,
    versions_json: String,
    budget_json: String,
    capability_json: String,
    seed: [u8; 8],
    progress: Vec<String>,
}

/// Execute a fresh solve run against a validated project.
///
/// `canonical_source` must be the project's canonical s-expression render
/// (retained as the run's input artifact so resume can re-verify identity).
/// `clock` returns monotonic seconds and is only used for budget accounting
/// and reporting, never identity.
///
/// # Errors
/// [`SolveRefusal`] with a stable code; when the refusal happened inside a
/// stage its evidence is retained in the ledger first.
pub fn run_solve(
    ledger: &Ledger,
    gate: &CancelGate,
    clock: &mut dyn FnMut() -> f64,
    project: &DecodedProject,
    progress_sink: &mut Vec<String>,
) -> Result<SolveOutcome, SolveRefusal> {
    let findings = project.findings();
    if !findings.is_empty() {
        return Err(SolveRefusal::plain(
            "cli-solve-project-invalid",
            format!("the project has {} validation findings", findings.len()),
            "run `frankensim validate` and repair every finding before solve",
        ));
    }
    let run = SolveRunId::derive(project);
    let mut engine = SolveEngine::open(ledger, gate, clock, project, run)?;
    let state = SolveDriverState {
        run: *run.as_bytes(),
        project: *engine.project_hash.as_bytes(),
        consumed_core_s: 0.0,
        consumed_wall_s: 0.0,
        completed: Vec::new(),
    };
    let outcome = engine.drive(state, StageContext::default());
    progress_sink.append(&mut engine.progress);
    outcome
}

/// Resume a previously started run from its durable state.
///
/// # Errors
/// [`SolveRefusal`]; identity mismatches between the requested run id, the
/// retained project source, and the retained driver state all refuse.
pub fn resume_solve(
    ledger: &Ledger,
    gate: &CancelGate,
    clock: &mut dyn FnMut() -> f64,
    run_id_hex: &str,
    progress_sink: &mut Vec<String>,
) -> Result<SolveOutcome, SolveRefusal> {
    let run = SolveRunId::parse_hex(run_id_hex).ok_or_else(|| {
        SolveRefusal::plain(
            "cli-solve-run-id",
            format!("`{run_id_hex}` is not a 64-hex run id"),
            "pass the run id printed by `frankensim solve`",
        )
    })?;
    let (state, state_artifact) = load_latest_state(ledger, run)?;
    let canonical_source = load_retained_project(ledger, &state)?;
    let project = fs_project::parse_sexpr(&canonical_source).map_err(|error| {
        SolveRefusal::plain(
            "cli-solve-resume-identity",
            format!(
                "the retained project source no longer parses strictly: {} ({})",
                error.code, error.detail
            ),
            "the ledger is inconsistent; verify artifact integrity",
        )
    })?;
    if !project.findings().is_empty() {
        return Err(SolveRefusal::plain(
            "cli-solve-resume-identity",
            "the retained project source no longer validates cleanly",
            "the ledger is inconsistent; verify artifact integrity",
        ));
    }
    let rederived = SolveRunId::derive(&project);
    if rederived != run {
        return Err(SolveRefusal::plain(
            "cli-solve-resume-identity",
            format!(
                "retained project re-derives run `{}` but resume requested `{}`",
                rederived.to_hex(),
                run.to_hex()
            ),
            "resuming across a changed project is refused; start a fresh solve",
        ));
    }
    if state.completed.len() >= SolveStage::ALL.len() {
        return Err(SolveRefusal::plain(
            "cli-solve-resume-complete",
            format!("run `{}` already completed every stage", run.to_hex()),
            "use `frankensim report <run-id>` once reporting ships (f85xj.6.9)",
        ));
    }
    let mut engine = SolveEngine::open(ledger, gate, clock, &project, run)?;
    // Restore budget continuity: re-charge the recorded consumption so the
    // resumed run continues under the same grant instead of resetting it.
    if state.consumed_core_s > 0.0 || state.consumed_wall_s > 0.0 {
        let enforcement = engine.charge(
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
                run: Some(run.to_hex()),
                recorded_op: None,
            });
        }
    }
    let context = StageContext::rebuild(ledger, &state)?;
    let outcome = engine.drive(state, context);
    progress_sink.append(&mut engine.progress);
    outcome
}

impl<'a> SolveEngine<'a> {
    fn open(
        ledger: &'a Ledger,
        gate: &'a CancelGate,
        clock: &'a mut dyn FnMut() -> f64,
        project: &'a DecodedProject,
        run: SolveRunId,
    ) -> Result<SolveEngine<'a>, SolveRefusal> {
        if ledger.in_transaction() {
            return Err(SolveRefusal::plain(
                "cli-solve-ledger-transaction",
                "the ledger connection is already inside a caller-owned transaction",
                "commit or roll back before solve so stage groups stay atomic",
            ));
        }
        let (versions_json, budget_json, capability_json, seed) = explicits(&project.spec)
            .map_err(|refusal| SolveRefusal {
                code: "cli-solve-project-invalid",
                stage: None,
                what: refusal.what,
                fix: refusal.fix,
                dependency: None,
                run: Some(run.to_hex()),
                recorded_op: None,
            })?;
        let token = derive_capability_token(&project.spec, run)?;
        let governor = Governor::new();
        let session = token.session;
        let open_id = governor
            .session_open_id(session, &format!("solve-open-{}", run.to_hex()))
            .map_err(|error| session_refusal(run, &error))?;
        governor
            .open_session_declared(open_id, token.clone())
            .map_err(|error| session_refusal(run, &error))?;
        Ok(SolveEngine {
            ledger,
            gate,
            clock,
            spec: &project.spec,
            canonical_source: &project.canonical,
            project_hash: project.hash(),
            run,
            governor,
            session,
            token,
            versions_json,
            budget_json,
            capability_json,
            seed,
            progress: Vec::new(),
        })
    }

    // The stage loop reads as one narrative: gate poll, gap check,
    // capability check, body, measure, persist, charge, enforce. Splitting
    // it would scatter the invariant ordering the CONTRACT documents.
    #[allow(clippy::too_many_lines)]
    fn drive(
        &mut self,
        mut state: SolveDriverState,
        mut context: StageContext,
    ) -> Result<SolveOutcome, SolveRefusal> {
        let prior_stages = u32::try_from(state.completed.len()).expect("stage count fits u32");
        let mut summaries = Vec::new();
        // A resumed run re-warns only fractions it has not yet crossed.
        let mut warned = [false; BUDGET_WARN_FRACTIONS.len()];
        for (index, &fraction) in BUDGET_WARN_FRACTIONS.iter().enumerate() {
            warned[index] = state.consumed_wall_s >= self.token.wall_s * fraction;
        }
        for stage in SolveStage::ALL.into_iter().skip(state.completed.len()) {
            if self.gate.is_requested() {
                return Err(self.cancelled_refusal(&state));
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
                    run: Some(self.run.to_hex()),
                    recorded_op: None,
                };
                return Err(self.record_refusal(&state, stage, refusal));
            }
            if !self.token.grants_op(stage.verb()) {
                let refusal = SolveRefusal::staged(
                    "cli-solve-capability",
                    stage,
                    format!("the session token does not grant `{}`", stage.verb()),
                    "derive the token from project budgets through the solve driver",
                );
                return Err(self.record_refusal(&state, stage, refusal));
            }
            let started = (self.clock)();
            let body = match stage {
                SolveStage::ImportVerify => self.stage_import_verify(&mut context),
                SolveStage::Assign => self.stage_assign(&context),
                _ => unreachable!("gap stages returned above"),
            };
            let receipt_json = match body {
                Ok(receipt) => receipt,
                Err(refusal) => return Err(self.record_refusal(&state, stage, refusal)),
            };
            let elapsed = ((self.clock)() - started).max(0.0);
            #[allow(clippy::cast_precision_loss)]
            let charge = Charge {
                core_s: elapsed * SOLVE_CORES as f64,
                mem_peak_bytes: 0,
                wall_s: elapsed,
            };
            state.consumed_core_s += charge.core_s;
            state.consumed_wall_s += charge.wall_s;
            let (op_id, receipt_hash) = self
                .persist_stage(&state, stage, &receipt_json, &context)
                .map_err(|error| self.ledger_refusal(stage, &error))?;
            state.completed.push(CompletedStage {
                ordinal: stage.ordinal(),
                op_id,
                receipt: receipt_hash,
            });
            // Re-seal the state WITH the new stage and persist it in a second
            // small transaction bound to the same op via a follow-up link is
            // impossible (lineage seals at finish); instead the state sealed
            // inside persist_stage already includes this stage. Nothing to do
            // here; the push above mirrors what persist_stage encoded.
            summaries.push(StageSummary {
                stage: stage.name(),
                op_id,
                receipt: receipt_hash.to_hex(),
                wall_s: elapsed,
            });
            self.progress.push(progress_line(
                &self.run.to_hex(),
                stage.name(),
                stage.ordinal(),
                "ok",
                elapsed,
            ));
            let enforcement = self.charge(stage.verb(), charge)?;
            match enforcement {
                Enforcement::Ok => {
                    for (index, &fraction) in BUDGET_WARN_FRACTIONS.iter().enumerate() {
                        if !warned[index] && state.consumed_wall_s >= self.token.wall_s * fraction {
                            warned[index] = true;
                            self.progress.push(budget_warning_line(
                                &self.run.to_hex(),
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
                        .persist_terminal(&state, &status)
                        .map_err(|error| self.ledger_refusal(stage, &error))?;
                    return Ok(SolveOutcome {
                        run: self.run.to_hex(),
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
            .persist_terminal(&state, &status)
            .map_err(|error| self.ledger_refusal(SolveStage::Qoi, &error))?;
        Ok(SolveOutcome {
            run: self.run.to_hex(),
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
                &format!("solve-meter-{}-{meter_key}", self.run.to_hex()),
            )
            .map_err(|error| session_refusal(self.run, &error))?;
        let receipt = self
            .governor
            .charge(report_id, charge)
            .map_err(|error| session_refusal(self.run, &error))?;
        Ok(receipt.enforcement().clone())
    }

    fn cancelled_refusal(&self, state: &SolveDriverState) -> SolveRefusal {
        SolveRefusal {
            code: "cli-solve-cancelled",
            stage: SolveStage::from_ordinal(
                u32::try_from(state.completed.len()).expect("stage count fits u32"),
            )
            .map(SolveStage::name),
            what: format!(
                "cancellation observed after {} completed stage(s); the completed prefix \
                 is durable",
                state.completed.len()
            ),
            fix: format!(
                "resume with `frankensim solve --resume {} <ledger>`",
                self.run.to_hex()
            ),
            dependency: None,
            run: Some(self.run.to_hex()),
            recorded_op: None,
        }
    }

    /// Verify retained import evidence against the pinned project.
    fn stage_import_verify(&mut self, context: &mut StageContext) -> Result<String, SolveRefusal> {
        let stage = SolveStage::ImportVerify;
        let geometry = self.spec.geometry.as_deref().unwrap_or(&[]);
        if geometry.is_empty() {
            return Err(SolveRefusal::staged(
                "cli-solve-import-evidence",
                stage,
                "the project declares no geometry to verify",
                "declare geometry rows and import them before solve",
            ));
        }
        let expected: Vec<String> = geometry.iter().map(geometry_source_identity).collect();
        let (import_op, summary) = find_import_summary(self.ledger, self.project_hash)
            .map_err(|error| self.ledger_refusal(stage, &error))?
            .ok_or_else(|| {
                SolveRefusal::staged(
                    "cli-solve-import-evidence",
                    stage,
                    format!(
                        "no completed geometry import for project `{}` exists in the ledger",
                        self.project_hash.to_hex()
                    ),
                    "run `frankensim import` for this exact project first",
                )
            })?;
        let entries = parse_import_summary(&summary).map_err(|what| {
            SolveRefusal::staged(
                "cli-solve-import-evidence",
                stage,
                what,
                "the retained import summary is not in the expected schema; re-import",
            )
        })?;
        let mut found: Vec<&str> = entries
            .iter()
            .map(|entry| entry.source_identity.as_str())
            .collect();
        found.sort_unstable();
        let mut wanted: Vec<&str> = expected.iter().map(String::as_str).collect();
        wanted.sort_unstable();
        if found != wanted {
            return Err(SolveRefusal::staged(
                "cli-solve-import-evidence",
                stage,
                format!(
                    "retained import covers sources [{}] but the project declares [{}]",
                    found.join(", "),
                    wanted.join(", ")
                ),
                "re-run `frankensim import` so every declared geometry row is retained",
            ));
        }
        // Re-hash every retained artifact so the run's authority is content
        // verification, not row presence.
        for entry in &entries {
            for (label, hash) in [
                ("raw source", entry.raw_source),
                ("promotion receipt", entry.promotion_receipt),
                ("promoted mesh", entry.promoted_mesh),
                ("assignment report", entry.assignment_report),
            ] {
                verify_artifact(self.ledger, hash, label).map_err(|what| {
                    SolveRefusal::staged(
                        "cli-solve-import-evidence",
                        stage,
                        what,
                        "the ledger no longer carries the imported evidence intact; re-import",
                    )
                })?;
            }
        }
        context.verified_imports = entries;
        let mut receipt = format!(
            "{{\"schema\":{},\"run\":{},\"project_hash\":{},\"import_op\":{import_op},\"verified\":[",
            json_string(IMPORT_VERIFY_RECEIPT_SCHEMA),
            json_string(&self.run.to_hex()),
            json_string(&self.project_hash.to_hex()),
        );
        for (index, entry) in context.verified_imports.iter().enumerate() {
            if index > 0 {
                receipt.push(',');
            }
            let _ = write!(
                receipt,
                "{{\"role\":{},\"source_identity\":{},\"raw_source\":{},\"promotion_receipt\":{},\"promoted_mesh\":{},\"assignment_report\":{}}}",
                json_string(&entry.role),
                json_string(&entry.source_identity),
                json_string(&entry.raw_source.to_hex()),
                json_string(&entry.promotion_receipt.to_hex()),
                json_string(&entry.promoted_mesh.to_hex()),
                json_string(&entry.assignment_report.to_hex()),
            );
        }
        receipt.push_str(
            "],\"authority\":\"re-hashed retained import evidence\",\"no_claim\":\"does not \
             prove the imported geometry is watertight, meshable, or physically meaningful\"}",
        );
        Ok(receipt)
    }

    /// Bind verified assignment evidence to the run's declared targets.
    fn stage_assign(&mut self, context: &StageContext) -> Result<String, SolveRefusal> {
        let stage = SolveStage::Assign;
        let assignments = self.spec.assignments.as_deref().unwrap_or(&[]);
        if assignments.is_empty() {
            return Err(SolveRefusal::staged(
                "cli-solve-assignment",
                stage,
                "the project declares no geometry assignments",
                "declare assignments before solve",
            ));
        }
        let mut rows = String::new();
        for (index, assignment) in assignments.iter().enumerate() {
            let imported = context
                .verified_imports
                .iter()
                .find(|entry| entry.role == assignment.artifact)
                .ok_or_else(|| {
                    SolveRefusal::staged(
                        "cli-solve-assignment",
                        stage,
                        format!(
                            "assignment target `{}` references geometry role `{}` with no \
                         verified import",
                            assignment.target, assignment.artifact
                        ),
                        "verify the import stage covers every assigned geometry role",
                    )
                })?;
            if index > 0 {
                rows.push(',');
            }
            let _ = write!(
                rows,
                "{{\"target\":{},\"artifact_role\":{},\"length_unit\":{},\"promoted_mesh\":{},\"assignment_report\":{}}}",
                json_string(&assignment.target),
                json_string(&assignment.artifact),
                json_string(&assignment.length_unit),
                json_string(&imported.promoted_mesh.to_hex()),
                json_string(&imported.assignment_report.to_hex()),
            );
        }
        Ok(format!(
            "{{\"schema\":{},\"run\":{},\"bindings\":[{rows}],\"authority\":\"declared targets \
             bound to verified import evidence\",\"no_claim\":\"selector re-resolution against \
             the mesh is the import's retained report; this stage does not re-run it\"}}",
            json_string(ASSIGN_RECEIPT_SCHEMA),
            json_string(&self.run.to_hex()),
        ))
    }

    /// Persist one completed stage as a ledgered op: stage receipt, sealed
    /// driver state (including this stage), and lineage links, atomically.
    fn persist_stage(
        &mut self,
        state_before: &SolveDriverState,
        stage: SolveStage,
        receipt_json: &str,
        context: &StageContext,
    ) -> Result<(i64, ContentHash), LedgerError> {
        let ordinal = stage.ordinal();
        let ir = format!(
            "{{\"schema\":{},\"stage\":{},\"ordinal\":{ordinal},\"run\":{},\"project\":{},\"driver_version\":{}}}",
            json_string(SOLVE_STAGE_SCHEMA),
            json_string(stage.name()),
            json_string(&self.run.to_hex()),
            json_string(&self.project_hash.to_hex()),
            SOLVE_DRIVER_VERSION,
        );
        self.ledger.begin()?;
        let result = (|| -> Result<(i64, ContentHash), LedgerError> {
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
            if ordinal == 0 {
                let source = self.ledger.put_artifact(
                    PROJECT_SOURCE_KIND,
                    self.canonical_source.as_bytes(),
                    None,
                )?;
                self.ledger.link(op, &source.hash, EdgeRole::In)?;
            }
            if stage == SolveStage::ImportVerify {
                for entry in &context.verified_imports {
                    self.ledger.link(op, &entry.promoted_mesh, EdgeRole::In)?;
                    self.ledger
                        .link(op, &entry.assignment_report, EdgeRole::In)?;
                }
            }
            let receipt =
                self.ledger
                    .put_artifact(STAGE_RECEIPT_KIND, receipt_json.as_bytes(), None)?;
            self.ledger.link(op, &receipt.hash, EdgeRole::Out)?;
            let mut sealed_state = state_before.clone();
            sealed_state.completed.push(CompletedStage {
                ordinal,
                op_id: op,
                receipt: receipt.hash,
            });
            let envelope = LegacySnapshotV1Adapter::<SolveDriverState>::seal(
                &sealed_state,
                envelope_provenance(self.run),
            );
            let state_artifact = self
                .ledger
                .put_artifact(STAGE_STATE_KIND, &envelope, None)?;
            self.ledger.link(op, &state_artifact.hash, EdgeRole::Out)?;
            self.ledger
                .finish_op(op, OpOutcome::Ok, None, i64::from(ordinal) * 2 + 1)?;
            Ok((op, receipt.hash))
        })();
        finish_solve_transaction(self.ledger, result)
    }

    /// Persist the run terminal receipt as its own ledgered op.
    fn persist_terminal(
        &mut self,
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
                    json_string(
                        "raising budgets changes the project identity and starts a fresh \
                         run `frankensim solve <project> <ledger>`; completed artifacts \
                         deduplicate by content"
                    ),
                ),
            ),
            SolveRunStatus::Cancelled => ("cancelled", String::new()),
        };
        let mut stages_json = String::new();
        for (index, stage) in state.completed.iter().enumerate() {
            if index > 0 {
                stages_json.push(',');
            }
            let name = SolveStage::from_ordinal(stage.ordinal).map_or("unknown", SolveStage::name);
            let _ = write!(
                stages_json,
                "{{\"stage\":{},\"op\":{},\"receipt\":{}}}",
                json_string(name),
                stage.op_id,
                json_string(&stage.receipt.to_hex()),
            );
        }
        let receipt_json = format!(
            "{{\"schema\":{},\"run\":{},\"project_hash\":{},\"status\":{}{detail},\"stages\":[{stages_json}],\"consumed_wall_s\":{},\"consumed_core_s\":{},\"no_claim\":\"stage receipts carry their own authority; this record is run bookkeeping\"}}",
            json_string(SOLVE_RUN_RECEIPT_SCHEMA),
            json_string(&self.run.to_hex()),
            json_string(&self.project_hash.to_hex()),
            json_string(status_name),
            state.consumed_wall_s,
            state.consumed_core_s,
        );
        let ir = format!(
            "{{\"schema\":{},\"stage\":\"terminal\",\"ordinal\":{},\"run\":{},\"project\":{},\"driver_version\":{}}}",
            json_string(SOLVE_STAGE_SCHEMA),
            SolveStage::ALL.len(),
            json_string(&self.run.to_hex()),
            json_string(&self.project_hash.to_hex()),
            SOLVE_DRIVER_VERSION,
        );
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
            let receipt =
                self.ledger
                    .put_artifact(RUN_RECEIPT_KIND, receipt_json.as_bytes(), None)?;
            self.ledger.link(op, &receipt.hash, EdgeRole::Out)?;
            self.ledger.finish_op(op, OpOutcome::Ok, None, 101)?;
            Ok(receipt.hash)
        })();
        finish_solve_transaction(self.ledger, result)
    }

    /// Retain a stage refusal as a ledgered error op, mirroring the import
    /// verb's durable-refusal doctrine.
    fn record_refusal(
        &mut self,
        state: &SolveDriverState,
        stage: SolveStage,
        mut refusal: SolveRefusal,
    ) -> SolveRefusal {
        refusal.run = Some(self.run.to_hex());
        let diagnostic = format!(
            "{{\"schema\":\"frankensim.cli.solve-refusal.v1\",\"run\":{},\"stage\":{},\"code\":{},\"what\":{},\"fix\":{}{}}}",
            json_string(&self.run.to_hex()),
            json_string(stage.name()),
            json_string(refusal.code),
            json_string(&refusal.what),
            json_string(&refusal.fix),
            refusal.dependency.map_or(String::new(), |dependency| {
                format!(",\"dependency\":{}", json_string(dependency))
            }),
        );
        let ir = format!(
            "{{\"schema\":{},\"stage\":{},\"ordinal\":{},\"run\":{},\"project\":{},\"driver_version\":{}}}",
            json_string(SOLVE_STAGE_SCHEMA),
            json_string(stage.name()),
            stage.ordinal(),
            json_string(&self.run.to_hex()),
            json_string(&self.project_hash.to_hex()),
            SOLVE_DRIVER_VERSION,
        );
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
            let artifact =
                self.ledger
                    .put_artifact("solve-refusal", diagnostic.as_bytes(), None)?;
            self.ledger.link(op, &artifact.hash, EdgeRole::Out)?;
            if stage.ordinal() == 0 && state.completed.is_empty() {
                let source = self.ledger.put_artifact(
                    PROJECT_SOURCE_KIND,
                    self.canonical_source.as_bytes(),
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
        match finish_solve_transaction(self.ledger, recorded) {
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

    fn ledger_refusal(&self, stage: SolveStage, error: &LedgerError) -> SolveRefusal {
        SolveRefusal {
            code: "cli-solve-ledger",
            stage: Some(stage.name()),
            what: format!("ledger transaction failed: {error}"),
            fix: "repair the ledger or contention failure, verify ledger lint/integrity, \
                  and retry"
                .to_string(),
            dependency: None,
            run: Some(self.run.to_hex()),
            recorded_op: None,
        }
    }
}

impl StageContext {
    /// Rebuild the in-memory context from retained stage receipts.
    fn rebuild(ledger: &Ledger, state: &SolveDriverState) -> Result<StageContext, SolveRefusal> {
        let mut context = StageContext::default();
        let Some(first) = state
            .completed
            .iter()
            .find(|stage| stage.ordinal == SolveStage::ImportVerify.ordinal())
        else {
            return Ok(context);
        };
        let bytes = ledger
            .get_artifact_bounded(&first.receipt, MAX_RECEIPT_READ_BYTES)
            .map_err(|error| {
                SolveRefusal::plain(
                    "cli-solve-ledger",
                    format!("reading the retained import-verify receipt failed: {error}"),
                    "verify ledger integrity and retry",
                )
            })?
            .ok_or_else(|| {
                SolveRefusal::plain(
                    "cli-solve-resume-identity",
                    "the driver state names an import-verify receipt the ledger does not carry",
                    "verify ledger integrity; the run cannot resume without its evidence",
                )
            })?;
        let text = String::from_utf8(bytes).map_err(|_| {
            SolveRefusal::plain(
                "cli-solve-resume-identity",
                "the retained import-verify receipt is not UTF-8",
                "verify ledger integrity",
            )
        })?;
        context.verified_imports = parse_import_summary(&text).map_err(|what| {
            SolveRefusal::plain(
                "cli-solve-resume-identity",
                format!("the retained import-verify receipt does not parse: {what}"),
                "verify ledger integrity",
            )
        })?;
        Ok(context)
    }
}

fn session_refusal(run: SolveRunId, error: &SessionError) -> SolveRefusal {
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

fn envelope_provenance(run: SolveRunId) -> u64 {
    run.session_u64()
}

fn push_framed(preimage: &mut Vec<u8>, bytes: &[u8]) {
    preimage.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    preimage.extend_from_slice(bytes);
}

fn progress_line(run: &str, stage: &str, ordinal: u32, status: &str, wall_s: f64) -> String {
    format!(
        "{{\"schema\":\"frankensim.cli.solve-progress.v1\",\"run\":{},\"stage\":{},\"ordinal\":{ordinal},\"status\":{},\"wall_s\":{wall_s}}}",
        json_string(run),
        json_string(stage),
        json_string(status),
    )
}

fn budget_warning_line(run: &str, fraction: f64, used: f64, granted: f64) -> String {
    format!(
        "{{\"schema\":\"frankensim.cli.solve-progress.v1\",\"run\":{},\"stage\":\"budget\",\"status\":\"warning\",\"fraction\":{fraction},\"used_wall_s\":{used},\"granted_wall_s\":{granted}}}",
        json_string(run),
    )
}

fn finish_solve_transaction<T>(
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

/// Locate the latest completed import op for this project and return its
/// summary receipt text.
fn find_import_summary(
    ledger: &Ledger,
    project_hash: ContentHash,
) -> Result<Option<(i64, String)>, LedgerError> {
    let needle = format!("\"project_hash\":{}", json_string(&project_hash.to_hex()));
    let mut ids = ledger.visible_op_ids(fs_ledger::MAIN_BRANCH, None)?;
    ids.sort_unstable_by(|a, b| b.cmp(a));
    for id in ids {
        let Some(row) = ledger.op(id)? else { continue };
        if row.outcome.as_deref() != Some("ok") {
            continue;
        }
        let edges = ledger.op_artifact_edges_bounded(id, EDGE_SCAN_CAP)?;
        for edge in &edges.edges {
            if edge.role != EdgeRole::Out {
                continue;
            }
            let Some(info) = ledger.artifact_info(&edge.artifact)? else {
                continue;
            };
            if info.kind != IMPORT_SUMMARY_KIND {
                continue;
            }
            let Some(bytes) =
                ledger.get_artifact_bounded(&edge.artifact, MAX_RECEIPT_READ_BYTES)?
            else {
                continue;
            };
            let Ok(text) = String::from_utf8(bytes) else {
                continue;
            };
            if text.contains(&needle) {
                return Ok(Some((id, text)));
            }
        }
    }
    Ok(None)
}

/// Load the latest sealed driver state for a run, walking the run's own ops.
fn load_latest_state(
    ledger: &Ledger,
    run: SolveRunId,
) -> Result<(SolveDriverState, ContentHash), SolveRefusal> {
    let ledger_error = |error: LedgerError| {
        SolveRefusal::plain(
            "cli-solve-ledger",
            format!("scanning the ledger for the run failed: {error}"),
            "verify ledger integrity and retry",
        )
    };
    let mut ids = ledger
        .visible_op_ids(fs_ledger::MAIN_BRANCH, None)
        .map_err(ledger_error)?;
    ids.sort_unstable_by(|a, b| b.cmp(a));
    let mut best: Option<(SolveDriverState, ContentHash)> = None;
    for id in ids {
        let Some(row) = ledger.op(id).map_err(ledger_error)? else {
            continue;
        };
        if row.session.as_deref() != Some(run.as_bytes().as_slice()) {
            continue;
        }
        if row.outcome.as_deref() != Some("ok") {
            continue;
        }
        let edges = ledger
            .op_artifact_edges_bounded(id, EDGE_SCAN_CAP)
            .map_err(ledger_error)?;
        for edge in &edges.edges {
            if edge.role != EdgeRole::Out {
                continue;
            }
            let Some(info) = ledger.artifact_info(&edge.artifact).map_err(ledger_error)? else {
                continue;
            };
            if info.kind != STAGE_STATE_KIND {
                continue;
            }
            let Some(bytes) = ledger
                .get_artifact_bounded(&edge.artifact, MAX_STATE_ENVELOPE_BYTES)
                .map_err(ledger_error)?
            else {
                continue;
            };
            let expectation = LegacySnapshotExpectationV1::new(
                ContentId::of_bytes(&bytes),
                DRIVER_STATE_TYPE_ID_V1,
                DRIVER_STATE_SCHEMA_VERSION_V1,
                envelope_provenance(run),
            );
            let limits =
                LegacySnapshotLimitsV1::new(MAX_STATE_ENVELOPE_BYTES, STATE_HASH_POLL_BYTES);
            let opened = LegacySnapshotV1Adapter::<SolveDriverState>::open_expected(
                &bytes,
                expectation,
                limits,
                fs_blake3::identity::NeverCancel,
            )
            .map_err(|error| {
                SolveRefusal::plain(
                    "cli-solve-resume-identity",
                    format!("the retained driver state failed envelope admission: {error:?}"),
                    "verify ledger integrity; the run cannot resume without a valid state",
                )
            })?;
            let state = opened.state().clone();
            if state.run != *run.as_bytes() {
                return Err(SolveRefusal::plain(
                    "cli-solve-resume-identity",
                    "the retained driver state carries a different run identity",
                    "verify ledger integrity",
                ));
            }
            let stage_count = state.completed.len();
            let replace = best
                .as_ref()
                .is_none_or(|(existing, _)| existing.completed.len() < stage_count);
            if replace {
                best = Some((state, edge.artifact));
            }
        }
    }
    best.ok_or_else(|| {
        SolveRefusal::plain(
            "cli-solve-unknown-run",
            format!("no solve run `{}` exists in this ledger", run.to_hex()),
            "pass the run id printed by `frankensim solve` and the same ledger path",
        )
    })
}

/// Load the retained canonical project source for a run. The run identity
/// itself is re-verified by the caller re-deriving it from the returned
/// source; this loader only follows the state's first-stage lineage.
fn load_retained_project(
    ledger: &Ledger,
    state: &SolveDriverState,
) -> Result<String, SolveRefusal> {
    let ledger_error = |error: LedgerError| {
        SolveRefusal::plain(
            "cli-solve-ledger",
            format!("reading the retained project failed: {error}"),
            "verify ledger integrity and retry",
        )
    };
    let first = state.completed.first().ok_or_else(|| {
        SolveRefusal::plain(
            "cli-solve-resume-identity",
            "the driver state records no completed stage",
            "verify ledger integrity",
        )
    })?;
    let edges = ledger
        .op_artifact_edges_bounded(first.op_id, EDGE_SCAN_CAP)
        .map_err(ledger_error)?;
    for edge in &edges.edges {
        if edge.role != EdgeRole::In {
            continue;
        }
        let Some(info) = ledger.artifact_info(&edge.artifact).map_err(ledger_error)? else {
            continue;
        };
        if info.kind != PROJECT_SOURCE_KIND {
            continue;
        }
        let Some(bytes) = ledger
            .get_artifact_bounded(&edge.artifact, crate::MAX_PROJECT_BYTES)
            .map_err(ledger_error)?
        else {
            continue;
        };
        return String::from_utf8(bytes).map_err(|_| {
            SolveRefusal::plain(
                "cli-solve-resume-identity",
                "the retained project source is not UTF-8",
                "verify ledger integrity",
            )
        });
    }
    Err(SolveRefusal::plain(
        "cli-solve-resume-identity",
        "the run's first stage op does not link a retained project source",
        "verify ledger integrity; the run cannot resume without its pinned project",
    ))
}

/// Verify an artifact exists and its bytes re-hash to the retained identity.
fn verify_artifact(ledger: &Ledger, hash: ContentHash, label: &str) -> Result<(), String> {
    let bytes = ledger
        .get_artifact_bounded(&hash, MAX_EVIDENCE_READ_BYTES)
        .map_err(|error| format!("reading the retained {label} failed: {error}"))?
        .ok_or_else(|| {
            format!(
                "the retained {label} `{}` is missing from the ledger",
                hash.to_hex()
            )
        })?;
    if hash_bytes(&bytes) != hash {
        return Err(format!(
            "the retained {label} `{}` no longer hashes to its identity",
            hash.to_hex()
        ));
    }
    Ok(())
}

/// Extract the per-source entries from an import summary or import-verify
/// receipt. Both schemas carry the same field names in the same writer
/// order; the reader is strict about shape and refuses escaped content.
fn parse_import_summary(text: &str) -> Result<Vec<VerifiedImport>, String> {
    let (_, tail) = text
        .split_once("\"artifacts\":[")
        .or_else(|| text.split_once("\"verified\":["))
        .ok_or_else(|| "no artifacts/verified array".to_string())?;
    let (body, _) = tail
        .split_once(']')
        .ok_or_else(|| "unterminated artifact array".to_string())?;
    if body.trim().is_empty() {
        return Err("empty artifact array".to_string());
    }
    let mut entries = Vec::new();
    for object in body.split("},{") {
        let object = object.trim_start_matches('{').trim_end_matches('}');
        let field = |key: &str| -> Result<String, String> {
            let marker = format!("\"{key}\":\"");
            let (_, rest) = object
                .split_once(&marker)
                .ok_or_else(|| format!("missing field `{key}`"))?;
            let (value, _) = rest
                .split_once('"')
                .ok_or_else(|| format!("unterminated field `{key}`"))?;
            if value.contains('\\') {
                return Err(format!("field `{key}` contains escapes; refusing to parse"));
            }
            Ok(value.to_string())
        };
        let hash_field = |key: &str| -> Result<ContentHash, String> {
            let value = field(key)?;
            ContentHash::from_hex(&value)
                .ok_or_else(|| format!("field `{key}` is not a 64-hex content hash"))
        };
        entries.push(VerifiedImport {
            role: field("role")?,
            source_identity: field("source_identity")?,
            raw_source: hash_field("raw_source")?,
            promotion_receipt: hash_field("promotion_receipt")?,
            promoted_mesh: hash_field("promoted_mesh")?,
            assignment_report: hash_field("assignment_report")?,
        });
    }
    Ok(entries)
}
