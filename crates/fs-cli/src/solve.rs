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
use fs_io::{MESH_ASSIGNMENT_SEMANTICS_VERSION, MeshSelector, NamedFaceGroup};
use fs_ledger::{
    ContentHash, EdgeRole, ExecMode, FiveExplicits, Ledger, LedgerError, MAIN_BRANCH,
    OpArtifactEdge, OpOutcome, OpRow, hash_bytes,
};
use fs_project::{DecodedProject, GeometryArtifact, ProjectSpec, geometry_source_identity};
use fs_session::{CapabilityToken, Charge, Enforcement, Governor, SessionError, SessionId};

use crate::import::{explicits, json_string};

/// Domain separating solve-run identity derivation from every other hash.
pub const SOLVE_RUN_IDENTITY_DOMAIN: &str = "org.frankensim.fs-cli.solve-run.v1";
/// Driver semantics version bound into run identity and driver state.
pub const SOLVE_DRIVER_VERSION: u32 = 2;

const SOLVE_STAGE_SCHEMA: &str = "frankensim.cli.solve-stage.v1";
const SOLVE_RUN_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-run-receipt.v1";
const IMPORT_VERIFY_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-import-verify-receipt.v1";
const ASSIGN_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-assignment-binding.v1";
const IMPORT_IR_SCHEMA: &str = "frankensim.cli.geometry-import.v1";
const IMPORT_SUMMARY_SCHEMA: &str = "frankensim.cli.geometry-import-receipt.v1";

const PROJECT_SOURCE_KIND: &str = "solve-project-source";
const STAGE_STATE_KIND: &str = "solve-stage-state";
const STAGE_RECEIPT_KIND: &str = "solve-stage-receipt";
const RUN_RECEIPT_KIND: &str = "solve-run-receipt";
const IMPORT_SUMMARY_KIND: &str = "geometry-import-run-receipt";
const IMPORT_RAW_KIND: &str = "geometry-source";
const IMPORT_PROMOTION_KIND: &str = "geometry-import-receipt";
const IMPORT_MESH_KIND: &str = "geometry-mesh-ply";
const IMPORT_ASSIGNMENT_KIND: &str = "geometry-assignment-report";

const IMPORT_SOURCE_LABEL_AUTHORITY: &str = "caller-reported";
const IMPORT_SUMMARY_AUTHORITY: &str = "retained-import-and-assignment-evidence";
const IMPORT_SUMMARY_NO_CLAIM: &str = "the ledger binds exact raw bytes, fs-io receipts, one \
    promoted finite tessellation, and the project assignment report; the project row's legacy \
    FNV hook and a caller-supplied source label do not authenticate custody, physical/CAD \
    sameness, continuum coverage, units, or topology beyond the retained lower-layer claims";
const IMPORT_VERIFY_AUTHORITY: &str = "re-hashed retained import evidence";
const IMPORT_VERIFY_NO_CLAIM: &str =
    "does not prove the imported geometry is watertight, meshable, or physically meaningful";
const ASSIGNMENT_REPORT_AUTHORITY: &str = "finite-tessellation-selection";
const ASSIGNMENT_REPORT_NO_CLAIM: &str = "selectors classify the supplied finite tessellation \
    only; caller-supplied source and subject identities are retained but not authenticated; no \
    between-facet, continuum-topology, self-intersection, CAD-semantic, or \
    physical-region-sameness claim is made";

/// Historical type identity of the driver state inside the legacy v1
/// snapshot envelope (`b"fsclisol"` as big-endian bytes).
const DRIVER_STATE_TYPE_ID_V1: u64 = 0x6673_636c_6973_6f6c;
const DRIVER_STATE_SCHEMA_VERSION_V1: u32 = 1;

/// Whole-envelope cap for a retained driver-state snapshot.
const MAX_STATE_ENVELOPE_BYTES: u64 = 4 * 1024 * 1024;
const STATE_HASH_POLL_BYTES: u32 = 64 * 1024;

/// Read cap for stage receipts and import summaries consumed on resume.
const MAX_RECEIPT_READ_BYTES: u64 = 4 * 1024 * 1024;
/// Materialization cap for retained mesh and assignment-report evidence that
/// solve must parse. Raw sources and opaque promotion receipts are verified
/// incrementally and are not constrained by this parse envelope.
const MAX_PARSED_EVIDENCE_BYTES: u64 = 64 * 1024 * 1024;
/// Maximum opaque lower-layer receipt bytes solve will stream for one source.
const MAX_OPAQUE_IMPORT_RECEIPT_BYTES: u64 = 4 * 1024 * 1024;
/// Absolute project-wide solve-verification admitted-input/work envelope. The
/// effective cap is the smaller of this value and the project's declared
/// memory-budget value. This byte-total preflight is not a peak-allocation
/// proof: parsed PLY validation can temporarily retain input, decoded soup, and
/// canonical re-encoding at the same time.
const MAX_TOTAL_SOLVE_EVIDENCE_BYTES: u64 = 512 * 1024 * 1024;
/// Edge scan cap per operation while locating retained evidence.
const EDGE_SCAN_CAP: usize = 1024;
/// Largest geometry set the solve evidence contract can attest completely.
/// One import operation has at most `4 * sources + 1` typed edges, so 255
/// sources fit under the ledger's 1024-edge bounded scan.
const SOLVE_MAX_IMPORT_SOURCES: usize = 255;

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
///
/// This public shape is a codec and inspection surface only. Constructing a
/// value does not make the named operation or receipt eligible for resume;
/// [`resume_solve`] independently re-attests the complete retained lineage.
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
///
/// Codec validity is not resume authority. Callers may construct and seal
/// this public state, but [`resume_solve`] treats every decoded field as
/// untrusted until the ordered operation, receipt, checkpoint, project, and
/// predecessor lineage all match records emitted by this driver.
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
    /// Exact import op selected by import verification.
    import_op: Option<i64>,
    /// Exact versioned import summary consumed by import verification.
    import_summary: Option<ContentHash>,
    /// Verified imports: (role, source identity, promoted mesh, report).
    verified_imports: Vec<VerifiedImport>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct VerifiedImport {
    role: String,
    source_identity: String,
    raw_source: ContentHash,
    promotion_receipt: ContentHash,
    promoted_mesh: ContentHash,
    assignment_report: ContentHash,
}

#[derive(Debug)]
struct ImportSummary {
    op_id: i64,
    artifact: ContentHash,
    entries: Vec<VerifiedImport>,
}

#[derive(Debug, Clone)]
struct ImportIrAttestation {
    limits: ImportIrLimits,
    sources: Vec<ImportIrSource>,
}

#[derive(Debug, Clone)]
struct ImportIrSource {
    length_unit: String,
    named_groups: Vec<NamedFaceGroup>,
    promotion_policy: ImportPromotionPolicy,
}

#[derive(Debug, Clone, Copy)]
enum ImportPromotionPolicy {
    Mesh { max_hole_edges: usize },
    FacetedStep { root_id: u64, target_h_bits: u64 },
}

#[derive(Debug)]
struct VerifiedResume {
    state: SolveDriverState,
    state_artifact: ContentHash,
    project: DecodedProject,
    context: StageContext,
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
    let outcome = engine.drive(state, StageContext::default(), None);
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
    let verified = load_latest_state(ledger, run)?;
    let VerifiedResume {
        state,
        state_artifact,
        project,
        context,
    } = verified;
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
    let outcome = engine.drive(state, context, Some(state_artifact));
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
        mut predecessor_state: Option<ContentHash>,
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
            if !started.is_finite() {
                let refusal = SolveRefusal::staged(
                    "cli-solve-budget",
                    stage,
                    "the solve clock returned a non-finite stage start",
                    "supply a finite monotonic clock; no checkpoint was sealed",
                );
                return Err(self.record_refusal(&state, stage, refusal));
            }
            let body = match stage {
                SolveStage::ImportVerify => self.stage_import_verify(&mut context),
                SolveStage::Assign => self.stage_assign(&context),
                _ => unreachable!("gap stages returned above"),
            };
            let receipt_json = match body {
                Ok(receipt) => receipt,
                Err(refusal) => return Err(self.record_refusal(&state, stage, refusal)),
            };
            let finished = (self.clock)();
            let elapsed = finished - started;
            if !finished.is_finite() || !elapsed.is_finite() || elapsed < 0.0 {
                let refusal = SolveRefusal::staged(
                    "cli-solve-budget",
                    stage,
                    format!(
                        "the solve clock did not produce a finite monotonic interval \
                         (start={started}, finish={finished})"
                    ),
                    "supply a finite monotonic clock; no checkpoint was sealed",
                );
                return Err(self.record_refusal(&state, stage, refusal));
            }
            #[allow(clippy::cast_precision_loss)]
            let charge = Charge {
                core_s: elapsed * SOLVE_CORES as f64,
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
                let refusal = SolveRefusal::staged(
                    "cli-solve-budget",
                    stage,
                    "the accumulated solve consumption is non-finite or negative",
                    "repair the clock/budget source; no checkpoint was sealed",
                );
                return Err(self.record_refusal(&state, stage, refusal));
            }
            state.consumed_core_s = consumed_core_s;
            state.consumed_wall_s = consumed_wall_s;
            let (op_id, receipt_hash, state_hash) = self
                .persist_stage(&state, stage, &receipt_json, &context, predecessor_state)
                .map_err(|error| self.ledger_refusal(stage, &error))?;
            state.completed.push(CompletedStage {
                ordinal: stage.ordinal(),
                op_id,
                receipt: receipt_hash,
            });
            predecessor_state = Some(state_hash);
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
        let summary = match find_import_summary(self.ledger, self.spec, self.project_hash) {
            Ok(Some(summary)) => summary,
            Ok(None) => {
                return Err(SolveRefusal::staged(
                    "cli-solve-import-evidence",
                    stage,
                    format!(
                        "no completed geometry import for project `{}` exists in the ledger",
                        self.project_hash.to_hex()
                    ),
                    "run `frankensim import` for this exact project first",
                ));
            }
            Err(ImportSummaryError::Unsupported(what)) => {
                return Err(SolveRefusal::staged(
                    "cli-solve-import-envelope",
                    stage,
                    what,
                    "reduce or split the retained geometry evidence to the documented solve envelope",
                ));
            }
            Err(ImportSummaryError::Ledger(error)) => {
                return Err(self.ledger_refusal(stage, &error));
            }
            Err(ImportSummaryError::Invalid(what)) => {
                return Err(SolveRefusal::staged(
                    "cli-solve-import-evidence",
                    stage,
                    format!("the retained import candidate is internally inconsistent: {what}"),
                    "re-run `frankensim import` for this exact project",
                ));
            }
        };
        let ImportSummary {
            op_id: import_op,
            artifact: import_summary,
            entries,
        } = summary;
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
        // Candidate validation above streamed or materialized every retained
        // artifact exactly once under the solve evidence envelope.
        context.import_op = Some(import_op);
        context.import_summary = Some(import_summary);
        context.verified_imports = entries;
        Ok(import_verify_receipt(
            self.run,
            self.project_hash,
            import_op,
            &context.verified_imports,
        ))
    }

    /// Bind verified assignment evidence to the run's declared targets.
    fn stage_assign(&mut self, context: &StageContext) -> Result<String, SolveRefusal> {
        assignment_receipt(self.spec, context, self.run)
    }

    /// Persist one completed stage as a ledgered op: stage receipt, sealed
    /// driver state (including this stage), and lineage links, atomically.
    fn persist_stage(
        &mut self,
        state_before: &SolveDriverState,
        stage: SolveStage,
        receipt_json: &str,
        context: &StageContext,
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
        let ir = solve_stage_ir(stage, self.run, self.project_hash);
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
                    PROJECT_SOURCE_KIND,
                    self.canonical_source.as_bytes(),
                    None,
                )?;
                self.ledger.link(op, &source.hash, EdgeRole::In)?;
            }
            if stage == SolveStage::ImportVerify {
                let import_summary =
                    context.import_summary.ok_or_else(|| LedgerError::Invalid {
                        field: "solve_import_summary".to_string(),
                        problem: "import verification did not retain its exact summary input"
                            .to_string(),
                    })?;
                self.ledger.link(op, &import_summary, EdgeRole::In)?;
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
        let ir = solve_stage_ir(stage, self.run, self.project_hash);
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

fn solve_stage_ir(stage: SolveStage, run: SolveRunId, project_hash: ContentHash) -> String {
    format!(
        "{{\"schema\":{},\"stage\":{},\"ordinal\":{},\"run\":{},\"project\":{},\"driver_version\":{}}}",
        json_string(SOLVE_STAGE_SCHEMA),
        json_string(stage.name()),
        stage.ordinal(),
        json_string(&run.to_hex()),
        json_string(&project_hash.to_hex()),
        SOLVE_DRIVER_VERSION,
    )
}

fn import_verify_receipt(
    run: SolveRunId,
    project_hash: ContentHash,
    import_op: i64,
    entries: &[VerifiedImport],
) -> String {
    let mut receipt = format!(
        "{{\"schema\":{},\"run\":{},\"project_hash\":{},\"import_op\":{import_op},\"verified\":[",
        json_string(IMPORT_VERIFY_RECEIPT_SCHEMA),
        json_string(&run.to_hex()),
        json_string(&project_hash.to_hex()),
    );
    for (index, entry) in entries.iter().enumerate() {
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
    let _ = write!(
        receipt,
        "],\"authority\":{},\"no_claim\":{}}}",
        json_string(IMPORT_VERIFY_AUTHORITY),
        json_string(IMPORT_VERIFY_NO_CLAIM),
    );
    receipt
}

fn assignment_receipt(
    spec: &ProjectSpec,
    context: &StageContext,
    run: SolveRunId,
) -> Result<String, SolveRefusal> {
    let stage = SolveStage::Assign;
    let assignments = spec.assignments.as_deref().unwrap_or(&[]);
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
        json_string(&run.to_hex()),
    ))
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

#[derive(Debug)]
enum ImportSummaryError {
    Ledger(LedgerError),
    Invalid(String),
    Unsupported(String),
}

impl From<LedgerError> for ImportSummaryError {
    fn from(error: LedgerError) -> Self {
        Self::Ledger(error)
    }
}

/// Locate the latest completed canonical import op for this exact project.
///
/// Invalid decoys are ignored; database/read failures and truncated edge
/// scans are not. A candidate becomes usable only after its versioned IR,
/// Five Explicits, complete typed edge set, and whole summary payload match
/// the import driver's writer contract.
fn find_import_summary(
    ledger: &Ledger,
    spec: &ProjectSpec,
    project_hash: ContentHash,
) -> Result<Option<ImportSummary>, ImportSummaryError> {
    let mut ids = ledger.visible_op_ids(MAIN_BRANCH, None)?;
    ids.sort_unstable_by(|a, b| b.cmp(a));
    'candidate: for id in ids {
        let Some(row) = ledger.op(id)? else {
            continue;
        };
        let attestation = match attest_import_row(ledger, spec, id, &row) {
            Ok(attestation) => attestation,
            Err(ImportSummaryError::Invalid(_)) => continue,
            Err(error @ ImportSummaryError::Unsupported(_))
            | Err(error @ ImportSummaryError::Ledger(_)) => return Err(error),
        };
        let edges = ledger.op_artifact_edges_bounded(id, EDGE_SCAN_CAP)?;
        if edges.truncated {
            return Err(ImportSummaryError::Ledger(LedgerError::Invalid {
                field: "geometry_import_edges".to_string(),
                problem: format!(
                    "operation {id} exceeds the complete {EDGE_SCAN_CAP}-edge import scan"
                ),
            }));
        }
        let expected_edge_count = attestation
            .sources
            .len()
            .checked_mul(4)
            .and_then(|count| count.checked_add(1))
            .expect("solve source cap makes import edge count representable");
        if edges.edges.len() != expected_edge_count {
            continue;
        }
        let mut raw_inputs = 0usize;
        let mut promotion_outputs = 0usize;
        let mut mesh_outputs = 0usize;
        let mut assignment_outputs = 0usize;
        let mut summary_artifact = None;
        for edge in &edges.edges {
            let Some(info) = ledger.artifact_info(&edge.artifact)? else {
                continue 'candidate;
            };
            match (edge.role, info.kind.as_str()) {
                (EdgeRole::In, IMPORT_RAW_KIND) => raw_inputs += 1,
                (EdgeRole::Out, IMPORT_PROMOTION_KIND) => promotion_outputs += 1,
                (EdgeRole::Out, IMPORT_MESH_KIND) => mesh_outputs += 1,
                (EdgeRole::Out, IMPORT_ASSIGNMENT_KIND) => assignment_outputs += 1,
                (EdgeRole::Out, IMPORT_SUMMARY_KIND) if summary_artifact.is_none() => {
                    summary_artifact = Some(edge.artifact);
                }
                _ => continue 'candidate,
            }
        }
        let source_count = attestation.sources.len();
        if raw_inputs != source_count
            || promotion_outputs != source_count
            || mesh_outputs != source_count
            || assignment_outputs != source_count
        {
            continue;
        }
        let Some(summary_artifact) = summary_artifact else {
            continue;
        };
        match validate_import_evidence(
            ledger,
            spec,
            project_hash,
            id,
            summary_artifact,
            &attestation,
        ) {
            Ok(summary) => return Ok(Some(summary)),
            Err(ImportSummaryError::Invalid(_)) => {}
            Err(error @ ImportSummaryError::Unsupported(_))
            | Err(error @ ImportSummaryError::Ledger(_)) => return Err(error),
        }
    }
    Ok(None)
}

/// Load and independently attest the latest sealed driver state for a run.
///
/// The legacy envelope expectation used while decoding is only a bounded
/// codec check: its byte identity is derived from the candidate bytes and
/// grants no authority. Resume eligibility comes exclusively from
/// [`validate_resume_candidate`], which re-attests the complete canonical
/// operation and checkpoint chain before a governor is opened.
fn load_latest_state(ledger: &Ledger, run: SolveRunId) -> Result<VerifiedResume, SolveRefusal> {
    let mut ids = ledger
        .visible_op_ids(MAIN_BRANCH, None)
        .map_err(|error| resume_ledger("scanning the ledger for the run failed", error))?;
    ids.sort_unstable_by(|a, b| b.cmp(a));
    let mut best: Option<VerifiedResume> = None;
    let mut best_is_ambiguous = false;
    for id in ids {
        let Some(row) = ledger
            .op(id)
            .map_err(|error| resume_ledger("scanning the ledger for the run failed", error))?
        else {
            continue;
        };
        if !is_supported_stage_discovery_row(&row, run) {
            continue;
        }
        let edges = resume_edges(ledger, id)?;
        for edge in &edges {
            if edge.role != EdgeRole::Out {
                continue;
            }
            let Some(info) = ledger.artifact_info(&edge.artifact).map_err(|error| {
                resume_ledger("reading a retained driver-state descriptor failed", error)
            })?
            else {
                continue;
            };
            if info.kind != STAGE_STATE_KIND {
                continue;
            }
            let state = decode_driver_state(ledger, run, edge.artifact)?;
            validate_state_shape(&state, run)?;
            if best
                .as_ref()
                .is_some_and(|existing| existing.state.completed.len() > state.completed.len())
            {
                continue;
            }
            let verified = validate_resume_candidate(ledger, run, state, edge.artifact, id)?;
            match best.as_ref() {
                Some(existing)
                    if existing.state.completed.len() == verified.state.completed.len() =>
                {
                    best_is_ambiguous = true;
                }
                Some(existing)
                    if existing.state.completed.len() > verified.state.completed.len() => {}
                _ => {
                    best = Some(verified);
                    best_is_ambiguous = false;
                }
            }
        }
    }
    let best = best.ok_or_else(|| {
        SolveRefusal::plain(
            "cli-solve-unknown-run",
            format!("no solve run `{}` exists in this ledger", run.to_hex()),
            "pass the run id printed by `frankensim solve` and the same ledger path",
        )
    })?;
    if best_is_ambiguous {
        return Err(resume_identity(format!(
            "the ledger carries competing independently valid {}-stage checkpoints for run `{}`",
            best.state.completed.len(),
            run.to_hex()
        )));
    }
    Ok(best)
}

fn is_supported_stage_discovery_row(row: &OpRow, run: SolveRunId) -> bool {
    if row.id <= 0
        || row.session.as_deref() != Some(run.as_bytes().as_slice())
        || row.outcome.as_deref() != Some("ok")
        || row.diag.is_some()
        || row.ir.len() > 512
    {
        return false;
    }
    let Ok(stage) = parse_stage_discovery_ir(&row.ir, run) else {
        return false;
    };
    stage.gap_dependency().is_none()
        && row.t_start == i64::from(stage.ordinal()) * 2
        && row.t_end == Some(i64::from(stage.ordinal()) * 2 + 1)
}

fn parse_stage_discovery_ir(ir: &str, run: SolveRunId) -> Result<SolveStage, String> {
    let mut cursor = JsonCursor::new(ir);
    cursor.expect("{\"schema\":")?;
    if cursor.parse_string()? != SOLVE_STAGE_SCHEMA {
        return Err("solve-stage discovery schema does not match the driver".to_string());
    }
    cursor.expect(",\"stage\":")?;
    let stage_name = cursor.parse_string()?;
    let stage = SolveStage::ALL
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
    if cursor.parse_u64()? != u64::from(SOLVE_DRIVER_VERSION) {
        return Err("solve-stage discovery driver version is unsupported".to_string());
    }
    cursor.expect("}")?;
    cursor.finish()?;
    Ok(stage)
}

#[allow(clippy::too_many_lines)]
fn validate_resume_candidate(
    ledger: &Ledger,
    run: SolveRunId,
    state: SolveDriverState,
    state_artifact: ContentHash,
    discovery_op: i64,
) -> Result<VerifiedResume, SolveRefusal> {
    validate_state_shape(&state, run)?;
    let first = state
        .completed
        .first()
        .expect("state shape requires a completed prefix");
    let first_edges = resume_edges(ledger, first.op_id)?;
    require_stage_edge_seal(ledger, first.op_id, first_edges.len())?;
    let (project, project_source) = load_retained_project(ledger, &first_edges)?;
    let project_hash = project.hash();
    if state.project != *project_hash.as_bytes() {
        return Err(resume_identity(
            "the retained driver state carries a different project identity",
        ));
    }
    let rederived = SolveRunId::derive(&project);
    if rederived != run {
        return Err(resume_identity(format!(
            "retained project re-derives run `{}` but resume requested `{}`",
            rederived.to_hex(),
            run.to_hex()
        )));
    }
    let (versions, budget, capability, seed) = explicits(&project.spec).map_err(|error| {
        resume_identity(format!(
            "the retained project cannot reproduce the stage Five Explicits: {}",
            error.what
        ))
    })?;

    let mut context = StageContext::default();
    let mut predecessor_state = None;
    let mut predecessor_checkpoint: Option<SolveDriverState> = None;
    for (index, completed) in state.completed.iter().enumerate() {
        let stage = SolveStage::ALL[index];
        if stage.gap_dependency().is_some() {
            return Err(resume_identity(format!(
                "driver version {SOLVE_DRIVER_VERSION} cannot have completed unavailable stage `{}`",
                stage.name()
            )));
        }
        let row = ledger
            .op(completed.op_id)
            .map_err(|error| resume_ledger("reading a completed solve operation failed", error))?
            .ok_or_else(|| {
                resume_identity(format!(
                    "completed stage {index} names missing operation {}",
                    completed.op_id
                ))
            })?;
        validate_stage_row(
            ledger,
            &row,
            stage,
            run,
            project_hash,
            &seed,
            &versions,
            &budget,
            &capability,
        )?;
        let edges = if index == 0 {
            first_edges.clone()
        } else {
            resume_edges(ledger, completed.op_id)?
        };
        require_stage_edge_seal(ledger, completed.op_id, edges.len())?;
        require_artifact_kind_resume(
            ledger,
            completed.receipt,
            STAGE_RECEIPT_KIND,
            "stage receipt",
        )?;
        if !has_edge(&edges, EdgeRole::Out, completed.receipt) {
            return Err(resume_identity(format!(
                "stage `{}` receipt {} is not an Out edge of operation {}",
                stage.name(),
                completed.receipt.to_hex(),
                completed.op_id
            )));
        }

        let checkpoint_outputs =
            artifacts_with_kind_resume(ledger, &edges, EdgeRole::Out, STAGE_STATE_KIND)?;
        if checkpoint_outputs.len() != 1 {
            return Err(resume_identity(format!(
                "stage `{}` operation {} has {} checkpoint outputs; exactly one is required",
                stage.name(),
                completed.op_id,
                checkpoint_outputs.len()
            )));
        }
        let checkpoint_hash = checkpoint_outputs[0];
        let checkpoint = decode_driver_state(ledger, run, checkpoint_hash)?;
        validate_checkpoint_prefix(&checkpoint, &state, index, predecessor_checkpoint.as_ref())?;

        if let Some(predecessor) = predecessor_state {
            if !has_edge(&edges, EdgeRole::In, predecessor) {
                return Err(resume_identity(format!(
                    "stage `{}` operation {} does not consume predecessor checkpoint {}",
                    stage.name(),
                    completed.op_id,
                    predecessor.to_hex()
                )));
            }
        }

        let receipt_text = read_text_resume(
            ledger,
            completed.receipt,
            MAX_RECEIPT_READ_BYTES,
            "stage receipt",
        )?;
        let mut expected_edges = vec![
            (EdgeRole::Out, completed.receipt),
            (EdgeRole::Out, checkpoint_hash),
        ];
        match stage {
            SolveStage::ImportVerify => {
                expected_edges.push((EdgeRole::In, project_source));
                let summary_inputs =
                    artifacts_with_kind_resume(ledger, &edges, EdgeRole::In, IMPORT_SUMMARY_KIND)?;
                if summary_inputs.len() != 1 {
                    return Err(resume_identity(format!(
                        "import-verify operation {} has {} geometry-import summary inputs; exactly one is required",
                        completed.op_id,
                        summary_inputs.len()
                    )));
                }
                let summary_hash = summary_inputs[0];
                let (import_op, receipt_entries) =
                    parse_import_verify_receipt(&receipt_text, run, project_hash).map_err(
                        |what| resume_identity(format!("invalid import-verify receipt: {what}")),
                    )?;
                let summary = validate_import_candidate(
                    ledger,
                    &project.spec,
                    project_hash,
                    import_op,
                    summary_hash,
                )
                .map_err(|error| match error {
                    ImportSummaryError::Ledger(error) => {
                        resume_ledger("re-attesting retained import evidence failed", error)
                    }
                    ImportSummaryError::Invalid(what) => resume_identity(format!(
                        "the retained import lineage is not canonical: {what}"
                    )),
                    ImportSummaryError::Unsupported(what) => resume_import_envelope(
                        run,
                        format!(
                            "the retained import lineage exceeds the solve evidence envelope: {what}"
                        ),
                    ),
                })?;
                if receipt_entries != summary.entries {
                    return Err(resume_identity(
                        "the import-verify receipt entries differ from its exact import summary",
                    ));
                }
                let expected_receipt =
                    import_verify_receipt(run, project_hash, summary.op_id, &summary.entries);
                if receipt_text != expected_receipt {
                    return Err(resume_identity(
                        "the retained import-verify receipt is not the canonical driver receipt",
                    ));
                }
                for entry in &summary.entries {
                    expected_edges.push((EdgeRole::In, entry.promoted_mesh));
                    expected_edges.push((EdgeRole::In, entry.assignment_report));
                }
                expected_edges.push((EdgeRole::In, summary_hash));
                context.import_op = Some(summary.op_id);
                context.import_summary = Some(summary.artifact);
                context.verified_imports = summary.entries;
            }
            SolveStage::Assign => {
                let expected_receipt =
                    assignment_receipt(&project.spec, &context, run).map_err(|error| {
                        resume_identity(format!(
                            "the retained assignment context cannot be reconstructed: {}",
                            error.what
                        ))
                    })?;
                if receipt_text != expected_receipt {
                    return Err(resume_identity(
                        "the retained assignment receipt is not the canonical driver receipt",
                    ));
                }
                expected_edges.push((
                    EdgeRole::In,
                    predecessor_state.expect("assign follows import-verify"),
                ));
            }
            _ => unreachable!("completed unavailable stages were refused above"),
        }
        require_exact_edges(stage, completed.op_id, &edges, &expected_edges)?;

        if index + 1 == state.completed.len() {
            if discovery_op != completed.op_id {
                return Err(resume_identity(format!(
                    "candidate checkpoint was discovered on operation {discovery_op}, not its last completed operation {}",
                    completed.op_id
                )));
            }
            if checkpoint_hash != state_artifact {
                return Err(resume_identity(format!(
                    "saved checkpoint {} is not the checkpoint Out edge of the last completed operation {}",
                    state_artifact.to_hex(),
                    completed.op_id
                )));
            }
            if checkpoint != state {
                return Err(resume_identity(
                    "the last completed operation's checkpoint does not equal the candidate state",
                ));
            }
        }
        predecessor_state = Some(checkpoint_hash);
        predecessor_checkpoint = Some(checkpoint);
    }

    Ok(VerifiedResume {
        state,
        state_artifact,
        project,
        context,
    })
}

fn validate_state_shape(state: &SolveDriverState, run: SolveRunId) -> Result<(), SolveRefusal> {
    if state.run != *run.as_bytes() {
        return Err(resume_identity(
            "the retained driver state carries a different run identity",
        ));
    }
    if !state.consumed_core_s.is_finite()
        || state.consumed_core_s < 0.0
        || !state.consumed_wall_s.is_finite()
        || state.consumed_wall_s < 0.0
        || state.consumed_core_s != state.consumed_wall_s
    {
        return Err(resume_identity(
            "the retained driver state carries non-finite or negative consumption",
        ));
    }
    if state.completed.is_empty() || state.completed.len() > SolveStage::ALL.len() {
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

#[allow(clippy::too_many_arguments)]
fn validate_stage_row(
    ledger: &Ledger,
    row: &OpRow,
    stage: SolveStage,
    run: SolveRunId,
    project_hash: ContentHash,
    seed: &[u8],
    versions: &str,
    budget: &str,
    capability: &str,
) -> Result<(), SolveRefusal> {
    let ordinal = stage.ordinal();
    if row.id <= 0
        || row.session.as_deref() != Some(run.as_bytes().as_slice())
        || row.outcome.as_deref() != Some("ok")
        || row.diag.is_some()
        || row.t_start != i64::from(ordinal) * 2
        || row.t_end != Some(i64::from(ordinal) * 2 + 1)
        || row.ir != solve_stage_ir(stage, run, project_hash)
        || row.seed != seed
        || row.versions != versions
        || row.budget != budget
        || row.capability != capability
    {
        return Err(resume_identity(format!(
            "operation {} does not match canonical stage `{}` semantics",
            row.id,
            stage.name()
        )));
    }
    let execution = ledger
        .op_execution_context(row.id)
        .map_err(|error| resume_ledger("reading solve operation execution context failed", error))?
        .ok_or_else(|| resume_identity(format!("operation {} has no execution context", row.id)))?;
    if execution.branch != MAIN_BRANCH || execution.exec_mode != ExecMode::Deterministic {
        return Err(resume_identity(format!(
            "operation {} is not a deterministic main-branch stage operation",
            row.id
        )));
    }
    ledger
        .op_content_identity(row.id)
        .map_err(|error| {
            resume_ledger("validating solve operation content identity failed", error)
        })?
        .ok_or_else(|| {
            resume_identity(format!(
                "operation {} has no typed content-identity sidecar",
                row.id
            ))
        })?;
    Ok(())
}

fn validate_checkpoint_prefix(
    checkpoint: &SolveDriverState,
    final_state: &SolveDriverState,
    index: usize,
    predecessor: Option<&SolveDriverState>,
) -> Result<(), SolveRefusal> {
    if checkpoint.run != final_state.run
        || checkpoint.project != final_state.project
        || checkpoint.completed.as_slice() != &final_state.completed[..=index]
        || !checkpoint.consumed_core_s.is_finite()
        || checkpoint.consumed_core_s < 0.0
        || !checkpoint.consumed_wall_s.is_finite()
        || checkpoint.consumed_wall_s < 0.0
        || checkpoint.consumed_core_s != checkpoint.consumed_wall_s
    {
        return Err(resume_identity(format!(
            "stage {index} checkpoint does not encode the exact trusted prefix"
        )));
    }
    if let Some(predecessor) = predecessor {
        if checkpoint.consumed_core_s < predecessor.consumed_core_s
            || checkpoint.consumed_wall_s < predecessor.consumed_wall_s
        {
            return Err(resume_identity(format!(
                "stage {index} checkpoint consumption regresses from its predecessor"
            )));
        }
    }
    Ok(())
}

fn resume_edges(ledger: &Ledger, op: i64) -> Result<Vec<OpArtifactEdge>, SolveRefusal> {
    let edges = ledger
        .op_artifact_edges_bounded(op, EDGE_SCAN_CAP)
        .map_err(|error| resume_ledger("reading solve operation edges failed", error))?;
    if edges.truncated {
        return Err(resume_identity(format!(
            "operation {op} exceeds the complete {EDGE_SCAN_CAP}-edge resume scan"
        )));
    }
    Ok(edges.edges)
}

fn require_stage_edge_seal(
    ledger: &Ledger,
    op: i64,
    edge_count: usize,
) -> Result<(), SolveRefusal> {
    let seal = ledger
        .op_artifact_edge_seal(op)
        .map_err(|error| resume_ledger("validating solve operation edge seal failed", error))?;
    if seal != Some(edge_count) {
        return Err(resume_identity(format!(
            "operation {op} lacks the driver's exact {edge_count}-edge lineage seal"
        )));
    }
    Ok(())
}

fn load_retained_project(
    ledger: &Ledger,
    edges: &[OpArtifactEdge],
) -> Result<(DecodedProject, ContentHash), SolveRefusal> {
    let mut sources = Vec::new();
    for edge in edges {
        if edge.role != EdgeRole::In {
            continue;
        }
        let info = ledger
            .artifact_info(&edge.artifact)
            .map_err(|error| resume_ledger("reading the retained project failed", error))?
            .ok_or_else(|| {
                resume_identity(format!(
                    "the first stage links missing artifact {}",
                    edge.artifact.to_hex()
                ))
            })?;
        if info.kind == PROJECT_SOURCE_KIND {
            sources.push(edge.artifact);
        }
    }
    if sources.len() != 1 {
        return Err(resume_identity(format!(
            "the first stage operation has {} retained project inputs; exactly one is required",
            sources.len()
        )));
    }
    let source_hash = sources[0];
    let bytes = ledger
        .get_artifact_bounded(&source_hash, crate::MAX_PROJECT_BYTES)
        .map_err(|error| resume_ledger("reading the retained project failed", error))?
        .ok_or_else(|| resume_identity("the retained project source is missing"))?;
    let source = String::from_utf8(bytes)
        .map_err(|_| resume_identity("the retained project source is not UTF-8"))?;
    let project = fs_project::parse_sexpr(&source).map_err(|error| {
        resume_identity(format!(
            "the retained project source no longer parses strictly: {} ({})",
            error.code, error.detail
        ))
    })?;
    if !project.findings().is_empty() || project.canonical != source {
        return Err(resume_identity(
            "the retained project source is not the exact canonical project",
        ));
    }
    Ok((project, source_hash))
}

fn resume_identity(what: impl Into<String>) -> SolveRefusal {
    SolveRefusal::plain(
        "cli-solve-resume-identity",
        what,
        "verify ledger integrity; codec validity alone cannot authorize resume",
    )
}

fn resume_import_envelope(run: SolveRunId, what: impl Into<String>) -> SolveRefusal {
    SolveRefusal {
        code: "cli-solve-import-envelope",
        stage: Some(SolveStage::ImportVerify.name()),
        what: what.into(),
        fix: "reduce or split the retained geometry evidence to the documented solve envelope"
            .to_string(),
        dependency: None,
        run: Some(run.to_hex()),
        recorded_op: None,
    }
}

fn resume_ledger(context: &str, error: LedgerError) -> SolveRefusal {
    SolveRefusal::plain(
        "cli-solve-ledger",
        format!("{context}: {error}"),
        "verify ledger integrity and retry",
    )
}

fn decode_driver_state(
    ledger: &Ledger,
    run: SolveRunId,
    artifact: ContentHash,
) -> Result<SolveDriverState, SolveRefusal> {
    let bytes = ledger
        .get_artifact_bounded(&artifact, MAX_STATE_ENVELOPE_BYTES)
        .map_err(|error| resume_ledger("reading a retained driver checkpoint failed", error))?
        .ok_or_else(|| {
            resume_identity(format!(
                "retained driver checkpoint {} is missing",
                artifact.to_hex()
            ))
        })?;
    if hash_bytes(&bytes) != artifact {
        return Err(resume_identity(format!(
            "retained driver checkpoint {} does not hash to its artifact identity",
            artifact.to_hex()
        )));
    }
    // This self-derived expectation is deliberately only a bounded codec and
    // corruption check. `validate_resume_candidate` supplies the independent
    // semantic/lineage admission that makes the decoded value usable.
    let expectation = LegacySnapshotExpectationV1::new(
        ContentId::of_bytes(&bytes),
        DRIVER_STATE_TYPE_ID_V1,
        DRIVER_STATE_SCHEMA_VERSION_V1,
        envelope_provenance(run),
    );
    let limits = LegacySnapshotLimitsV1::new(MAX_STATE_ENVELOPE_BYTES, STATE_HASH_POLL_BYTES);
    let opened = LegacySnapshotV1Adapter::<SolveDriverState>::open_expected(
        &bytes,
        expectation,
        limits,
        fs_blake3::identity::NeverCancel,
    )
    .map_err(|error| {
        resume_identity(format!(
            "the retained driver state failed bounded envelope admission: {error:?}"
        ))
    })?;
    Ok(opened.state().clone())
}

fn require_artifact_kind_resume(
    ledger: &Ledger,
    artifact: ContentHash,
    expected_kind: &str,
    label: &str,
) -> Result<(), SolveRefusal> {
    let info = ledger
        .artifact_info(&artifact)
        .map_err(|error| resume_ledger("reading an artifact descriptor failed", error))?
        .ok_or_else(|| {
            resume_identity(format!(
                "the retained {label} {} is missing",
                artifact.to_hex()
            ))
        })?;
    if info.kind != expected_kind {
        return Err(resume_identity(format!(
            "the retained {label} {} has kind `{}`, not `{expected_kind}`",
            artifact.to_hex(),
            info.kind
        )));
    }
    Ok(())
}

fn artifacts_with_kind_resume(
    ledger: &Ledger,
    edges: &[OpArtifactEdge],
    role: EdgeRole,
    kind: &str,
) -> Result<Vec<ContentHash>, SolveRefusal> {
    let mut matches = Vec::new();
    for edge in edges {
        if edge.role != role {
            continue;
        }
        let info = ledger
            .artifact_info(&edge.artifact)
            .map_err(|error| resume_ledger("reading an artifact descriptor failed", error))?
            .ok_or_else(|| {
                resume_identity(format!(
                    "operation lineage names missing artifact {}",
                    edge.artifact.to_hex()
                ))
            })?;
        if info.kind == kind {
            matches.push(edge.artifact);
        }
    }
    Ok(matches)
}

fn read_text_resume(
    ledger: &Ledger,
    artifact: ContentHash,
    cap: u64,
    label: &str,
) -> Result<String, SolveRefusal> {
    let bytes = ledger
        .get_artifact_bounded(&artifact, cap)
        .map_err(|error| resume_ledger(&format!("reading the retained {label} failed"), error))?
        .ok_or_else(|| {
            resume_identity(format!(
                "the retained {label} {} is missing",
                artifact.to_hex()
            ))
        })?;
    String::from_utf8(bytes).map_err(|_| {
        resume_identity(format!(
            "the retained {label} {} is not UTF-8",
            artifact.to_hex()
        ))
    })
}

fn has_edge(edges: &[OpArtifactEdge], role: EdgeRole, artifact: ContentHash) -> bool {
    edges
        .iter()
        .any(|edge| edge.role == role && edge.artifact == artifact)
}

fn edge_sets_match(edges: &[OpArtifactEdge], expected: &[(EdgeRole, ContentHash)]) -> bool {
    edges.len() == expected.len()
        && expected
            .iter()
            .all(|(role, artifact)| has_edge(edges, *role, *artifact))
        && edges.iter().all(|edge| {
            expected
                .iter()
                .any(|(role, artifact)| edge.role == *role && edge.artifact == *artifact)
        })
}

fn require_exact_edges(
    stage: SolveStage,
    op: i64,
    edges: &[OpArtifactEdge],
    expected: &[(EdgeRole, ContentHash)],
) -> Result<(), SolveRefusal> {
    if !edge_sets_match(edges, expected) {
        return Err(resume_identity(format!(
            "stage `{}` operation {op} has {} artifact edges, not the driver's exact {}-edge set",
            stage.name(),
            edges.len(),
            expected.len()
        )));
    }
    Ok(())
}

fn validate_import_candidate(
    ledger: &Ledger,
    spec: &ProjectSpec,
    project_hash: ContentHash,
    op: i64,
    summary_artifact: ContentHash,
) -> Result<ImportSummary, ImportSummaryError> {
    let row = ledger.op(op)?.ok_or_else(|| {
        ImportSummaryError::Invalid(format!("import summary names missing operation {op}"))
    })?;
    let attestation = attest_import_row(ledger, spec, op, &row)?;
    validate_import_evidence(
        ledger,
        spec,
        project_hash,
        op,
        summary_artifact,
        &attestation,
    )
}

fn attest_import_row(
    ledger: &Ledger,
    spec: &ProjectSpec,
    op: i64,
    row: &OpRow,
) -> Result<ImportIrAttestation, ImportSummaryError> {
    if op <= 0
        || row.id != op
        || row.session.is_some()
        || row.outcome.as_deref() != Some("ok")
        || row.diag.is_some()
        || row.t_start != 0
        || row.t_end != Some(1)
    {
        return Err(ImportSummaryError::Invalid(format!(
            "operation {op} is not a canonical completed import operation"
        )));
    }
    let attestation = validate_import_ir(&row.ir, spec).map_err(ImportSummaryError::Invalid)?;
    if attestation.sources.len() > SOLVE_MAX_IMPORT_SOURCES {
        return Err(ImportSummaryError::Unsupported(format!(
            "import operation {op} carries {} sources above the solve evidence cap {SOLVE_MAX_IMPORT_SOURCES}",
            attestation.sources.len()
        )));
    }
    let (versions, budget, capability, seed) = explicits(spec).map_err(|error| {
        ImportSummaryError::Invalid(format!(
            "project cannot reproduce import Five Explicits: {}",
            error.what
        ))
    })?;
    if row.seed.as_slice() != seed
        || row.versions != versions
        || row.budget != budget
        || row.capability != capability
    {
        return Err(ImportSummaryError::Invalid(format!(
            "import operation {op} does not carry the project's exact Five Explicits"
        )));
    }
    let execution = ledger.op_execution_context(op)?.ok_or_else(|| {
        ImportSummaryError::Invalid(format!("import operation {op} has no execution context"))
    })?;
    if execution.branch != MAIN_BRANCH || execution.exec_mode != ExecMode::Deterministic {
        return Err(ImportSummaryError::Invalid(format!(
            "import operation {op} is not deterministic on the main branch"
        )));
    }
    ledger.op_content_identity(op)?.ok_or_else(|| {
        ImportSummaryError::Invalid(format!(
            "import operation {op} has no typed content-identity sidecar"
        ))
    })?;
    Ok(attestation)
}

fn validate_import_evidence(
    ledger: &Ledger,
    spec: &ProjectSpec,
    project_hash: ContentHash,
    op: i64,
    summary_artifact: ContentHash,
    attestation: &ImportIrAttestation,
) -> Result<ImportSummary, ImportSummaryError> {
    let bounded = ledger.op_artifact_edges_bounded(op, EDGE_SCAN_CAP)?;
    if bounded.truncated {
        return Err(ImportSummaryError::Invalid(format!(
            "import operation {op} exceeds the complete {EDGE_SCAN_CAP}-edge scan"
        )));
    }
    if !has_edge(&bounded.edges, EdgeRole::Out, summary_artifact) {
        return Err(ImportSummaryError::Invalid(format!(
            "summary {} is not an Out edge of import operation {op}",
            summary_artifact.to_hex()
        )));
    }
    require_import_artifact_kind(
        ledger,
        summary_artifact,
        IMPORT_SUMMARY_KIND,
        "import summary",
    )?;
    let summary_len = import_artifact_len(ledger, summary_artifact, "import summary")?;
    if summary_len > MAX_RECEIPT_READ_BYTES {
        return Err(ImportSummaryError::Unsupported(format!(
            "import summary is {summary_len} bytes above the solve receipt envelope {MAX_RECEIPT_READ_BYTES}"
        )));
    }
    let bytes = ledger
        .get_artifact_bounded(&summary_artifact, MAX_RECEIPT_READ_BYTES)?
        .ok_or_else(|| {
            ImportSummaryError::Invalid(format!(
                "import summary {} is missing",
                summary_artifact.to_hex()
            ))
        })?;
    let text = String::from_utf8(bytes).map_err(|_| {
        ImportSummaryError::Invalid(format!(
            "import summary {} is not UTF-8",
            summary_artifact.to_hex()
        ))
    })?;
    let entries = parse_geometry_import_summary(&text, spec, project_hash)
        .map_err(ImportSummaryError::Invalid)?;
    let mut expected_edges = vec![(EdgeRole::Out, summary_artifact)];
    for entry in &entries {
        require_import_artifact_kind(ledger, entry.raw_source, IMPORT_RAW_KIND, "raw source")?;
        require_import_artifact_kind(
            ledger,
            entry.promotion_receipt,
            IMPORT_PROMOTION_KIND,
            "promotion receipt",
        )?;
        require_import_artifact_kind(
            ledger,
            entry.promoted_mesh,
            IMPORT_MESH_KIND,
            "promoted mesh",
        )?;
        require_import_artifact_kind(
            ledger,
            entry.assignment_report,
            IMPORT_ASSIGNMENT_KIND,
            "assignment report",
        )?;
        expected_edges.push((EdgeRole::In, entry.raw_source));
        expected_edges.push((EdgeRole::Out, entry.promotion_receipt));
        expected_edges.push((EdgeRole::Out, entry.promoted_mesh));
        expected_edges.push((EdgeRole::Out, entry.assignment_report));
    }
    if !edge_sets_match(&bounded.edges, &expected_edges) {
        return Err(ImportSummaryError::Invalid(format!(
            "import operation {op} has {} artifact edges, not the exact typed {}-edge set",
            bounded.edges.len(),
            expected_edges.len()
        )));
    }
    validate_import_admission_evidence(ledger, spec, summary_artifact, &entries, attestation)?;
    Ok(ImportSummary {
        op_id: op,
        artifact: summary_artifact,
        entries,
    })
}

fn require_import_artifact_kind(
    ledger: &Ledger,
    artifact: ContentHash,
    expected_kind: &str,
    label: &str,
) -> Result<(), ImportSummaryError> {
    let info = ledger.artifact_info(&artifact)?.ok_or_else(|| {
        ImportSummaryError::Invalid(format!("{label} artifact {} is missing", artifact.to_hex()))
    })?;
    if info.kind != expected_kind {
        return Err(ImportSummaryError::Invalid(format!(
            "{label} artifact {} has kind `{}`, not `{expected_kind}`",
            artifact.to_hex(),
            info.kind
        )));
    }
    Ok(())
}

#[allow(clippy::too_many_lines)]
fn validate_import_admission_evidence(
    ledger: &Ledger,
    spec: &ProjectSpec,
    summary_artifact: ContentHash,
    entries: &[VerifiedImport],
    attestation: &ImportIrAttestation,
) -> Result<(), ImportSummaryError> {
    let geometry = spec.geometry.as_deref().unwrap_or(&[]);
    if entries.len() != geometry.len() || attestation.sources.len() != geometry.len() {
        return Err(ImportSummaryError::Invalid(
            "import IR, summary, and project geometry cardinalities differ".to_string(),
        ));
    }
    let assignments = spec.assignments.as_deref().unwrap_or(&[]);
    if assignments.len() > attestation.limits.max_requests {
        return Err(ImportSummaryError::Invalid(format!(
            "project declares {} assignment requests above import IR max_requests {}",
            assignments.len(),
            attestation.limits.max_requests
        )));
    }
    let mut entity_violations = Vec::new();
    let entity_ids = spec.resolve_entities(&mut entity_violations);
    if !entity_violations.is_empty() {
        return Err(ImportSummaryError::Invalid(
            "project entity identities cannot be re-derived for assignment evidence".to_string(),
        ));
    }

    let source_cap = u64::try_from(attestation.limits.max_source_bytes).map_err(|_| {
        ImportSummaryError::Invalid(
            "import IR max_source_bytes is outside the ledger byte range".to_string(),
        )
    })?;
    let total_cap = u64::try_from(attestation.limits.max_total_source_bytes).map_err(|_| {
        ImportSummaryError::Invalid(
            "import IR max_total_source_bytes is outside the ledger byte range".to_string(),
        )
    })?;
    let project_memory_cap = spec
        .budgets
        .as_ref()
        .map(|budget| budget.memory_bytes)
        .ok_or_else(|| {
            ImportSummaryError::Invalid(
                "project has no memory budget for the solve evidence envelope".to_string(),
            )
        })?;
    let solve_total_cap = project_memory_cap.min(MAX_TOTAL_SOLVE_EVIDENCE_BYTES);
    let raw_stream_cap = source_cap.min(solve_total_cap);
    let mut preflight_raw_bytes = 0u64;
    let mut preflight_evidence_bytes =
        import_artifact_len(ledger, summary_artifact, "import summary")?;
    for (source_index, entry) in entries.iter().enumerate() {
        let raw_len = import_artifact_len(ledger, entry.raw_source, "raw source")?;
        if raw_len > source_cap {
            return Err(ImportSummaryError::Invalid(format!(
                "source {source_index} raw artifact is {raw_len} bytes above max_source_bytes {}",
                attestation.limits.max_source_bytes
            )));
        }
        preflight_raw_bytes = preflight_raw_bytes.checked_add(raw_len).ok_or_else(|| {
            ImportSummaryError::Invalid(
                "aggregate retained raw-source byte count overflowed u64".to_string(),
            )
        })?;
        let promotion_len =
            import_artifact_len(ledger, entry.promotion_receipt, "promotion receipt")?;
        if promotion_len > MAX_OPAQUE_IMPORT_RECEIPT_BYTES {
            return Err(ImportSummaryError::Unsupported(format!(
                "source {source_index} promotion receipt is {promotion_len} bytes above the solve opaque-receipt envelope {MAX_OPAQUE_IMPORT_RECEIPT_BYTES}"
            )));
        }
        let mesh_len = import_artifact_len(ledger, entry.promoted_mesh, "promoted mesh")?;
        if mesh_len > MAX_PARSED_EVIDENCE_BYTES {
            return Err(ImportSummaryError::Unsupported(format!(
                "source {source_index} promoted mesh is {mesh_len} bytes above the solve parse envelope {MAX_PARSED_EVIDENCE_BYTES}"
            )));
        }
        let report_len = import_artifact_len(ledger, entry.assignment_report, "assignment report")?;
        if report_len > MAX_PARSED_EVIDENCE_BYTES {
            return Err(ImportSummaryError::Unsupported(format!(
                "source {source_index} assignment report is {report_len} bytes above the solve parse envelope {MAX_PARSED_EVIDENCE_BYTES}"
            )));
        }
        for length in [raw_len, promotion_len, mesh_len, report_len] {
            preflight_evidence_bytes =
                preflight_evidence_bytes
                    .checked_add(length)
                    .ok_or_else(|| {
                        ImportSummaryError::Unsupported(
                            "project-wide solve evidence byte count overflowed u64".to_string(),
                        )
                    })?;
        }
    }
    if preflight_raw_bytes > total_cap {
        return Err(ImportSummaryError::Invalid(format!(
            "aggregate retained raw sources reached {preflight_raw_bytes} bytes above max_total_source_bytes {}",
            attestation.limits.max_total_source_bytes
        )));
    }
    if preflight_evidence_bytes > solve_total_cap {
        return Err(ImportSummaryError::Unsupported(format!(
            "retained import evidence totals {preflight_evidence_bytes} bytes above the effective solve envelope {solve_total_cap} (project memory {project_memory_cap}, hard cap {MAX_TOTAL_SOLVE_EVIDENCE_BYTES})"
        )));
    }

    let mut total_raw_bytes = 0u64;
    let mut total_selected_faces = 0usize;

    for (source_index, ((artifact, entry), source)) in geometry
        .iter()
        .zip(entries)
        .zip(&attestation.sources)
        .enumerate()
    {
        if entry.source_identity.len() > attestation.limits.max_label_bytes {
            return Err(ImportSummaryError::Invalid(format!(
                "source {source_index} identity exceeds max_label_bytes {}",
                attestation.limits.max_label_bytes
            )));
        }
        // These policy values are retained so they are never silently
        // discarded. Lower-layer receipts may retain the corresponding
        // repair/root/spacing values, but solve deliberately does not parse,
        // cross-check, or replay those promotion semantics. The contract
        // therefore makes no semantic claim for these values beyond exact
        // writer grammar and basic scalar admission.
        match source.promotion_policy {
            ImportPromotionPolicy::Mesh { max_hole_edges } => {
                let _unreplayed_mesh_repair_cap = max_hole_edges;
            }
            ImportPromotionPolicy::FacetedStep {
                root_id,
                target_h_bits,
            } => {
                let _unreplayed_step_policy = (root_id, target_h_bits);
            }
        }

        let raw_len = ledger
            .read_artifact_chunks_bounded(&entry.raw_source, raw_stream_cap, &mut |_| {})
            .map_err(|error| match error {
                LedgerError::ArtifactReadLimit { observed, .. } if observed > source_cap => {
                    ImportSummaryError::Invalid(format!(
                        "source {source_index} raw artifact exceeds max_source_bytes {}",
                        attestation.limits.max_source_bytes
                    ))
                }
                LedgerError::ArtifactReadLimit { observed, .. } => {
                    ImportSummaryError::Unsupported(format!(
                        "source {source_index} raw artifact is {observed} bytes above the effective solve stream envelope {raw_stream_cap}"
                    ))
                }
                error => ImportSummaryError::Ledger(error),
            })?
            .ok_or_else(|| {
                ImportSummaryError::Invalid(format!(
                    "source {source_index} raw artifact {} is missing",
                    entry.raw_source.to_hex()
                ))
            })?;
        total_raw_bytes = total_raw_bytes.checked_add(raw_len).ok_or_else(|| {
            ImportSummaryError::Invalid(
                "aggregate retained raw-source byte count overflowed u64".to_string(),
            )
        })?;
        if total_raw_bytes > total_cap {
            return Err(ImportSummaryError::Invalid(format!(
                "aggregate retained raw sources reached {total_raw_bytes} bytes above max_total_source_bytes {}",
                attestation.limits.max_total_source_bytes
            )));
        }
        ledger
            .read_artifact_chunks_bounded(
                &entry.promotion_receipt,
                MAX_OPAQUE_IMPORT_RECEIPT_BYTES,
                &mut |_| {},
            )?
            .ok_or_else(|| {
                ImportSummaryError::Invalid(format!(
                    "source {source_index} promotion receipt {} is missing",
                    entry.promotion_receipt.to_hex()
                ))
            })?;

        let mesh_bytes = read_parsed_import_artifact(
            ledger,
            entry.promoted_mesh,
            "promoted mesh",
            source_index,
        )?;
        preflight_canonical_ply(
            &mesh_bytes,
            attestation.limits.max_mesh_vertices,
            attestation.limits.max_mesh_faces,
        )
        .map_err(|what| {
            ImportSummaryError::Invalid(format!(
                "source {source_index} promoted PLY is not canonical: {what}"
            ))
        })?;
        let soup = fs_io::ply::read_ply(&mesh_bytes).map_err(|error| {
            ImportSummaryError::Invalid(format!(
                "source {source_index} promoted PLY does not parse: {error}"
            ))
        })?;
        if fs_io::ply::write_ply(&soup).as_bytes() != mesh_bytes {
            return Err(ImportSummaryError::Invalid(format!(
                "source {source_index} promoted PLY is not the exact canonical writer output"
            )));
        }
        if soup.positions.len() > attestation.limits.max_mesh_vertices
            || soup.triangles.len() > attestation.limits.max_mesh_faces
        {
            return Err(ImportSummaryError::Invalid(format!(
                "source {source_index} promoted mesh has {} vertices and {} faces above import IR caps {} and {}",
                soup.positions.len(),
                soup.triangles.len(),
                attestation.limits.max_mesh_vertices,
                attestation.limits.max_mesh_faces
            )));
        }
        for group in &source.named_groups {
            if group.faces.iter().any(|face| {
                usize::try_from(*face).map_or(true, |face| face >= soup.triangles.len())
            }) {
                return Err(ImportSummaryError::Invalid(format!(
                    "source {source_index} named group `{}` references a face outside the promoted mesh",
                    group.name
                )));
            }
        }

        let rows: Vec<_> = assignments
            .iter()
            .filter(|assignment| assignment.artifact == artifact.role)
            .collect();
        let mut expected_subjects = Vec::new();
        expected_subjects
            .try_reserve_exact(rows.len())
            .map_err(|_| {
                ImportSummaryError::Invalid(
                    "assignment subject-attestation allocation refused".to_string(),
                )
            })?;
        let mut geometric_requests = 0u64;
        for row in &rows {
            let entity = entity_ids.get(&row.target).ok_or_else(|| {
                ImportSummaryError::Invalid(format!(
                    "assignment target `{}` has no re-derived entity identity",
                    row.target
                ))
            })?;
            let subject = entity.token();
            if subject.len() > attestation.limits.max_label_bytes {
                return Err(ImportSummaryError::Invalid(format!(
                    "assignment subject `{subject}` exceeds max_label_bytes {}",
                    attestation.limits.max_label_bytes
                )));
            }
            if let MeshSelector::NamedGroup { name } = &row.selector {
                if name.len() > attestation.limits.max_label_bytes {
                    return Err(ImportSummaryError::Invalid(format!(
                        "named-group selector `{name}` exceeds max_label_bytes {}",
                        attestation.limits.max_label_bytes
                    )));
                }
            }
            if let MeshSelector::ExplicitFaceSet { faces, .. } = &row.selector {
                if faces.len() > attestation.limits.max_selected_faces {
                    return Err(ImportSummaryError::Invalid(format!(
                        "assignment `{}` explicit face set exceeds max_selected_faces {}",
                        row.target, attestation.limits.max_selected_faces
                    )));
                }
            }
            if matches!(
                &row.selector,
                MeshSelector::HalfSpace { .. }
                    | MeshSelector::Box { .. }
                    | MeshSelector::Cylinder { .. }
                    | MeshSelector::NearestDatum { .. }
            ) {
                geometric_requests = geometric_requests.checked_add(1).ok_or_else(|| {
                    ImportSummaryError::Invalid(
                        "geometric assignment request count overflowed u64".to_string(),
                    )
                })?;
            }
            expected_subjects.push((subject, row.allow_overlap));
        }
        let face_count = u64::try_from(soup.triangles.len()).map_err(|_| {
            ImportSummaryError::Invalid(
                "promoted mesh face count is outside the predicate-work range".to_string(),
            )
        })?;
        let predicate_tests = geometric_requests.checked_mul(face_count).ok_or_else(|| {
            ImportSummaryError::Invalid(
                "assignment predicate-work count overflowed u64".to_string(),
            )
        })?;
        if predicate_tests > attestation.limits.max_predicate_tests {
            return Err(ImportSummaryError::Invalid(format!(
                "source {source_index} requires {predicate_tests} predicate tests above max_predicate_tests {}",
                attestation.limits.max_predicate_tests
            )));
        }

        let report_bytes = read_parsed_import_artifact(
            ledger,
            entry.assignment_report,
            "assignment report",
            source_index,
        )?;
        let report_text = core::str::from_utf8(&report_bytes).map_err(|_| {
            ImportSummaryError::Invalid(format!(
                "source {source_index} assignment report is not UTF-8"
            ))
        })?;
        let selected = parse_assignment_report_counts(
            report_text,
            &entry.source_identity,
            &source.length_unit,
            &expected_subjects,
            soup.triangles.len(),
        )
        .map_err(|what| {
            ImportSummaryError::Invalid(format!(
                "source {source_index} assignment report is not canonical: {what}"
            ))
        })?;
        total_selected_faces = total_selected_faces.checked_add(selected).ok_or_else(|| {
            ImportSummaryError::Invalid(
                "project-wide selected-face count overflowed usize".to_string(),
            )
        })?;
        if total_selected_faces > attestation.limits.max_selected_faces {
            return Err(ImportSummaryError::Invalid(format!(
                "project-wide assignment reports select {total_selected_faces} faces above max_selected_faces {}",
                attestation.limits.max_selected_faces
            )));
        }
    }
    Ok(())
}

fn import_artifact_len(
    ledger: &Ledger,
    artifact: ContentHash,
    label: &str,
) -> Result<u64, ImportSummaryError> {
    ledger
        .artifact_info(&artifact)?
        .map(|info| info.len)
        .ok_or_else(|| {
            ImportSummaryError::Invalid(format!(
                "{label} artifact {} is missing",
                artifact.to_hex()
            ))
        })
}

fn read_parsed_import_artifact(
    ledger: &Ledger,
    artifact: ContentHash,
    label: &str,
    source_index: usize,
) -> Result<Vec<u8>, ImportSummaryError> {
    let info = ledger.artifact_info(&artifact)?.ok_or_else(|| {
        ImportSummaryError::Invalid(format!(
            "source {source_index} {label} artifact {} is missing",
            artifact.to_hex()
        ))
    })?;
    if info.len > MAX_PARSED_EVIDENCE_BYTES {
        return Err(ImportSummaryError::Unsupported(format!(
            "source {source_index} {label} is {} bytes above the solve parse envelope {MAX_PARSED_EVIDENCE_BYTES}",
            info.len
        )));
    }
    ledger
        .get_artifact_bounded(&artifact, MAX_PARSED_EVIDENCE_BYTES)?
        .ok_or_else(|| {
            ImportSummaryError::Invalid(format!(
                "source {source_index} {label} artifact {} disappeared",
                artifact.to_hex()
            ))
        })
}

fn preflight_canonical_ply(
    bytes: &[u8],
    max_vertices: usize,
    max_faces: usize,
) -> Result<(), String> {
    const END_HEADER: &[u8] = b"end_header\n";
    let header_end = bytes
        .windows(END_HEADER.len())
        .position(|window| window == END_HEADER)
        .and_then(|offset| offset.checked_add(END_HEADER.len()))
        .ok_or_else(|| "the exact `end_header\\n` terminator is missing".to_string())?;
    let header = core::str::from_utf8(&bytes[..header_end])
        .map_err(|_| "the header is not UTF-8".to_string())?;
    let header = header
        .strip_suffix('\n')
        .ok_or_else(|| "the header is not newline-terminated".to_string())?;
    let mut lines = header.split('\n');
    if lines.next() != Some("ply") || lines.next() != Some("format ascii 1.0") {
        return Err("the magic or ASCII format line differs from the fs-io writer".to_string());
    }
    let vertices = parse_canonical_ply_count(
        lines
            .next()
            .ok_or_else(|| "the vertex element line is missing".to_string())?,
        "element vertex ",
        "vertex",
    )?;
    if lines.next() != Some("property double x")
        || lines.next() != Some("property double y")
        || lines.next() != Some("property double z")
    {
        return Err("the vertex property lines differ from the fs-io writer".to_string());
    }
    let faces = parse_canonical_ply_count(
        lines
            .next()
            .ok_or_else(|| "the face element line is missing".to_string())?,
        "element face ",
        "face",
    )?;
    if lines.next() != Some("property list uchar uint vertex_indices")
        || lines.next() != Some("end_header")
        || lines.next().is_some()
    {
        return Err("the face property or header shape differs from the fs-io writer".to_string());
    }
    if vertices > max_vertices || faces > max_faces {
        return Err(format!(
            "declared counts {vertices} vertices and {faces} faces exceed frozen caps {max_vertices} and {max_faces}"
        ));
    }

    let expected_records = vertices
        .checked_add(faces)
        .ok_or_else(|| "the declared body-record count overflows usize".to_string())?;
    let body = &bytes[header_end..];
    if !body.is_empty() && body.last() != Some(&b'\n') {
        return Err("the body is not newline-terminated".to_string());
    }
    let body_records = body.iter().filter(|byte| **byte == b'\n').count();
    if body_records != expected_records {
        return Err(format!(
            "the body has {body_records} newline-terminated records but the header declares {expected_records}"
        ));
    }
    Ok(())
}

fn parse_canonical_ply_count(line: &str, prefix: &str, label: &str) -> Result<usize, String> {
    let spelling = line
        .strip_prefix(prefix)
        .ok_or_else(|| format!("the {label} element line differs from the fs-io writer"))?;
    if spelling.is_empty()
        || !spelling.bytes().all(|byte| byte.is_ascii_digit())
        || (spelling.len() > 1 && spelling.starts_with('0'))
    {
        return Err(format!(
            "the {label} count is not a canonical unsigned decimal"
        ));
    }
    let count = spelling
        .parse::<usize>()
        .map_err(|_| format!("the {label} count is outside usize"))?;
    if count.to_string() != spelling {
        return Err(format!(
            "the {label} count is not the writer's canonical decimal spelling"
        ));
    }
    Ok(count)
}

fn parse_assignment_report_counts(
    text: &str,
    expected_source: &str,
    expected_unit: &str,
    expected_assignments: &[(String, bool)],
    mesh_faces: usize,
) -> Result<usize, String> {
    let mut cursor = JsonCursor::new(text);
    cursor.expect("{\"kind\":")?;
    if cursor.parse_string()? != "mesh-assignment-receipt" {
        return Err("kind is not mesh-assignment-receipt".to_string());
    }
    cursor.expect(",\"version\":")?;
    if cursor.parse_string()? != MESH_ASSIGNMENT_SEMANTICS_VERSION {
        return Err("version does not match the fs-io writer".to_string());
    }
    cursor.expect(",\"source_identity\":")?;
    if cursor.parse_string()? != expected_source {
        return Err("source identity differs from the import summary".to_string());
    }
    cursor.expect(",\"length_unit\":")?;
    if cursor.parse_string()? != expected_unit {
        return Err("length unit differs from the import IR".to_string());
    }
    for field in [
        "source_mesh_fingerprint",
        "named_groups_fingerprint",
        "requests_fingerprint",
        "assignments_fingerprint",
    ] {
        cursor.expect(&format!(",\"{field}\":"))?;
        let fingerprint = cursor.parse_string()?;
        require_lower_hex_16(&fingerprint, field)?;
    }
    cursor.expect(",\"assignments\":[")?;
    let mut selected = 0usize;
    for (index, (expected_subject, expected_overlap)) in expected_assignments.iter().enumerate() {
        if index > 0 {
            cursor.expect(",")?;
        }
        cursor.expect("{\"subject\":")?;
        if cursor.parse_string()? != *expected_subject {
            return Err(format!("assignment row {index} has the wrong subject"));
        }
        cursor.expect(",\"selector_fingerprint\":")?;
        let selector_fingerprint = cursor.parse_string()?;
        require_lower_hex_16(&selector_fingerprint, "selector_fingerprint")?;
        cursor.expect(",\"face_count\":")?;
        let face_count = cursor.parse_usize()?;
        if face_count == 0 || face_count > mesh_faces {
            return Err(format!(
                "assignment row {index} face_count {face_count} is outside 1..={mesh_faces}"
            ));
        }
        selected = selected
            .checked_add(face_count)
            .ok_or_else(|| "assignment face_count sum overflowed usize".to_string())?;
        cursor.expect(",\"surface_area\":")?;
        let surface_area = cursor.parse_canonical_finite_f64()?;
        if surface_area <= 0.0 {
            return Err(format!(
                "assignment row {index} surface area is not positive"
            ));
        }
        cursor.expect(",\"enclosed_volume\":")?;
        if cursor.peek() == Some(b'n') {
            cursor.expect("null")?;
        } else {
            let _ = cursor.parse_canonical_finite_f64()?;
        }
        cursor.expect(",\"bounds_min\":[")?;
        parse_finite_vector3(&mut cursor)?;
        cursor.expect("],\"bounds_max\":[")?;
        parse_finite_vector3(&mut cursor)?;
        cursor.expect("],\"allow_overlap\":")?;
        cursor.expect(if *expected_overlap { "true" } else { "false" })?;
        cursor.expect("}")?;
    }
    cursor.expect("]")?;
    cursor.expect(",\"authority\":")?;
    if cursor.parse_string()? != ASSIGNMENT_REPORT_AUTHORITY {
        return Err("assignment authority differs from the writer".to_string());
    }
    cursor.expect(",\"no_claim\":")?;
    if cursor.parse_string()? != ASSIGNMENT_REPORT_NO_CLAIM {
        return Err("assignment no-claim differs from the writer".to_string());
    }
    cursor.expect("}")?;
    cursor.finish()?;
    Ok(selected)
}

fn require_lower_hex_16(value: &str, label: &str) -> Result<(), String> {
    if value.len() == 16
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        Ok(())
    } else {
        Err(format!("{label} is not 16 lowercase hex digits"))
    }
}

fn parse_finite_vector3(cursor: &mut JsonCursor<'_>) -> Result<(), String> {
    let _ = cursor.parse_canonical_finite_f64()?;
    cursor.expect(",")?;
    let _ = cursor.parse_canonical_finite_f64()?;
    cursor.expect(",")?;
    let _ = cursor.parse_canonical_finite_f64()?;
    Ok(())
}

fn validate_import_ir(ir: &str, spec: &ProjectSpec) -> Result<ImportIrAttestation, String> {
    let mut cursor = JsonCursor::new(ir);
    cursor.expect("{\"schema\":")?;
    let schema = cursor.parse_string()?;
    if schema != IMPORT_IR_SCHEMA {
        return Err(format!(
            "import IR schema is `{schema}`, not `{IMPORT_IR_SCHEMA}`"
        ));
    }
    cursor.expect(",\"project\":")?;
    let project_json = cursor.take_value()?;
    let canonical_project = fs_project::print_json(spec)
        .map_err(|error| format!("canonical project JSON failed: {error:?}"))?;
    if project_json != canonical_project {
        return Err("import IR does not embed the exact canonical project JSON".to_string());
    }
    let limits = parse_import_ir_limits(&mut cursor)?;
    cursor.expect(",\"sources\":[")?;
    let geometry = spec.geometry.as_deref().unwrap_or(&[]);
    if geometry.len() > limits.max_sources {
        return Err(format!(
            "import IR declares {} project sources above its max_sources {}",
            geometry.len(),
            limits.max_sources
        ));
    }
    let mut sources = Vec::new();
    sources
        .try_reserve_exact(geometry.len())
        .map_err(|_| "import IR source-attestation allocation refused".to_string())?;
    for (index, artifact) in geometry.iter().enumerate() {
        if index > 0 {
            cursor.expect(",")?;
        }
        cursor.expect("{\"source_identity\":")?;
        let source_identity = cursor.parse_string()?;
        let expected_identity = geometry_source_identity(artifact);
        if source_identity != expected_identity {
            return Err(format!(
                "import IR source {index} identity `{source_identity}` does not match `{expected_identity}`"
            ));
        }
        cursor.expect(",\"policy\":")?;
        sources.push(parse_import_ir_policy(
            &mut cursor,
            spec,
            artifact,
            index,
            limits,
        )?);
        cursor.expect("}")?;
    }
    cursor.expect("]}")?;
    cursor.finish()?;
    Ok(ImportIrAttestation { limits, sources })
}

#[derive(Debug, Clone, Copy)]
struct ImportIrLimits {
    max_sources: usize,
    max_source_bytes: usize,
    max_total_source_bytes: usize,
    max_mesh_vertices: usize,
    max_mesh_faces: usize,
    max_requests: usize,
    max_named_groups: usize,
    max_group_faces: usize,
    max_selected_faces: usize,
    max_predicate_tests: u64,
    max_label_bytes: usize,
}

fn parse_import_ir_limits(cursor: &mut JsonCursor<'_>) -> Result<ImportIrLimits, String> {
    cursor.expect(",\"limits\":{\"max_sources\":")?;
    let max_sources = cursor.parse_usize()?;
    cursor.expect(",\"max_source_bytes\":")?;
    let max_source_bytes = cursor.parse_usize()?;
    cursor.expect(",\"max_total_source_bytes\":")?;
    let max_total_source_bytes = cursor.parse_usize()?;
    if max_sources == 0 || max_source_bytes == 0 || max_total_source_bytes == 0 {
        return Err("import IR source-count and byte limits are not positive".to_string());
    }
    cursor.expect(",\"assignment\":{\"max_mesh_vertices\":")?;
    let max_mesh_vertices = cursor.parse_usize()?;
    cursor.expect(",\"max_mesh_faces\":")?;
    let max_mesh_faces = cursor.parse_usize()?;
    cursor.expect(",\"max_requests\":")?;
    let max_requests = cursor.parse_usize()?;
    cursor.expect(",\"max_named_groups\":")?;
    let max_named_groups = cursor.parse_usize()?;
    cursor.expect(",\"max_group_faces\":")?;
    let max_group_faces = cursor.parse_usize()?;
    cursor.expect(",\"max_selected_faces\":")?;
    let max_selected_faces = cursor.parse_usize()?;
    cursor.expect(",\"max_predicate_tests\":")?;
    let max_predicate_tests = cursor.parse_u64()?;
    cursor.expect(",\"max_label_bytes\":")?;
    let max_label_bytes = cursor.parse_usize()?;
    cursor.expect("}}")?;
    if [
        max_mesh_vertices,
        max_mesh_faces,
        max_requests,
        max_named_groups,
        max_group_faces,
        max_selected_faces,
        max_label_bytes,
    ]
    .contains(&0)
        || max_predicate_tests == 0
    {
        return Err("import IR assignment limits are not all positive".to_string());
    }
    Ok(ImportIrLimits {
        max_sources,
        max_source_bytes,
        max_total_source_bytes,
        max_mesh_vertices,
        max_mesh_faces,
        max_requests,
        max_named_groups,
        max_group_faces,
        max_selected_faces,
        max_predicate_tests,
        max_label_bytes,
    })
}

fn parse_import_ir_policy(
    cursor: &mut JsonCursor<'_>,
    spec: &ProjectSpec,
    artifact: &GeometryArtifact,
    source_index: usize,
    limits: ImportIrLimits,
) -> Result<ImportIrSource, String> {
    cursor.expect("{\"kind\":")?;
    let kind = cursor.parse_string()?;
    let source = match kind.as_str() {
        "mesh" => {
            if !matches!(artifact.format.as_str(), "stl" | "obj" | "ply") {
                return Err(format!(
                    "import IR source {source_index} uses mesh policy for project format `{}`",
                    artifact.format
                ));
            }
            cursor.expect(",\"length_unit\":")?;
            let length_unit = cursor.parse_string()?;
            validate_import_ir_unit(
                spec,
                artifact,
                source_index,
                &length_unit,
                limits.max_label_bytes,
            )?;
            cursor.expect(",\"max_hole_edges\":")?;
            let max_hole_edges = cursor.parse_usize()?;
            cursor.expect(",\"named_groups\":")?;
            let named_groups = parse_import_ir_named_groups(cursor, source_index, limits)?;
            ImportIrSource {
                length_unit,
                named_groups,
                promotion_policy: ImportPromotionPolicy::Mesh { max_hole_edges },
            }
        }
        "faceted-step" => {
            if artifact.format != "step" {
                return Err(format!(
                    "import IR source {source_index} uses faceted-step policy for project format `{}`",
                    artifact.format
                ));
            }
            cursor.expect(",\"root_id\":")?;
            let root_id = cursor.parse_u64()?;
            if root_id == 0 {
                return Err(format!(
                    "import IR source {source_index} has a zero faceted-step root"
                ));
            }
            cursor.expect(",\"length_unit\":")?;
            let length_unit = cursor.parse_string()?;
            validate_import_ir_unit(
                spec,
                artifact,
                source_index,
                &length_unit,
                limits.max_label_bytes,
            )?;
            cursor.expect(",\"target_h_bits\":")?;
            let target_h_bits = cursor.parse_string()?;
            if target_h_bits.len() != 16
                || !target_h_bits
                    .bytes()
                    .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
            {
                return Err(format!(
                    "import IR source {source_index} target_h_bits is not 16 lowercase hex digits"
                ));
            }
            let bits = u64::from_str_radix(&target_h_bits, 16).map_err(|_| {
                format!("import IR source {source_index} target_h_bits is outside u64")
            })?;
            let target_h = f64::from_bits(bits);
            if !target_h.is_finite() || target_h <= 0.0 {
                return Err(format!(
                    "import IR source {source_index} target spacing is not positive and finite"
                ));
            }
            cursor.expect(",\"named_groups\":")?;
            let named_groups = parse_import_ir_named_groups(cursor, source_index, limits)?;
            ImportIrSource {
                length_unit,
                named_groups,
                promotion_policy: ImportPromotionPolicy::FacetedStep {
                    root_id,
                    target_h_bits: bits,
                },
            }
        }
        _ => {
            return Err(format!(
                "import IR source {source_index} policy kind `{kind}` is not a versioned writer kind"
            ));
        }
    };
    cursor.expect("}")?;
    Ok(source)
}

fn validate_import_ir_unit(
    spec: &ProjectSpec,
    artifact: &GeometryArtifact,
    source_index: usize,
    unit: &str,
    max_label_bytes: usize,
) -> Result<(), String> {
    if unit.is_empty()
        || unit.trim() != unit
        || unit.len() > max_label_bytes
        || unit.chars().any(char::is_control)
    {
        return Err(format!(
            "import IR source {source_index} carries an invalid length-unit spelling"
        ));
    }
    for assignment in spec
        .assignments
        .as_deref()
        .unwrap_or(&[])
        .iter()
        .filter(|assignment| assignment.artifact == artifact.role)
    {
        if assignment.length_unit != unit {
            return Err(format!(
                "import IR source {source_index} unit `{unit}` differs from assignment `{}` unit `{}`",
                assignment.target, assignment.length_unit
            ));
        }
    }
    Ok(())
}

fn parse_import_ir_named_groups(
    cursor: &mut JsonCursor<'_>,
    source_index: usize,
    limits: ImportIrLimits,
) -> Result<Vec<NamedFaceGroup>, String> {
    cursor.expect("[")?;
    let mut groups = 0usize;
    let mut face_references = 0usize;
    let mut named_groups = Vec::new();
    if cursor.peek() != Some(b']') {
        loop {
            groups = groups
                .checked_add(1)
                .ok_or_else(|| "import IR named-group count overflowed usize".to_string())?;
            if groups > limits.max_named_groups {
                return Err(format!(
                    "import IR source {source_index} exceeds max_named_groups {}",
                    limits.max_named_groups
                ));
            }
            cursor.expect("{\"name\":")?;
            let name = cursor.parse_string()?;
            if name.is_empty()
                || name.trim() != name
                || name.len() > limits.max_label_bytes
                || name.chars().any(char::is_control)
            {
                return Err(format!(
                    "import IR source {source_index} named-group label is not nonempty, trim-canonical, control-free, and bounded"
                ));
            }
            named_groups
                .try_reserve(1)
                .map_err(|_| "import IR named-group allocation refused".to_string())?;
            cursor.expect(",\"faces\":[")?;
            if cursor.peek() == Some(b']') {
                return Err(format!(
                    "import IR source {source_index} named group has no faces"
                ));
            }
            let mut group_faces = Vec::new();
            loop {
                let face = cursor.parse_u64()?;
                let face = u32::try_from(face).map_err(|_| {
                    format!("import IR source {source_index} named-group face is outside u32")
                })?;
                face_references = face_references.checked_add(1).ok_or_else(|| {
                    "import IR named-group face count overflowed usize".to_string()
                })?;
                if face_references > limits.max_group_faces {
                    return Err(format!(
                        "import IR source {source_index} exceeds max_group_faces {}",
                        limits.max_group_faces
                    ));
                }
                group_faces
                    .try_reserve(1)
                    .map_err(|_| "import IR named-group face allocation refused".to_string())?;
                group_faces.push(face);
                match cursor.peek() {
                    Some(b',') => cursor.expect(",")?,
                    Some(b']') => break,
                    _ => {
                        return Err(cursor.problem("expected `,` or `]` after named-group face"));
                    }
                }
            }
            group_faces.sort_unstable();
            if let Some(duplicate) = group_faces
                .windows(2)
                .find(|window| window[0] == window[1])
                .map(|window| window[0])
            {
                return Err(format!(
                    "import IR source {source_index} named group repeats face {duplicate}"
                ));
            }
            named_groups.push(NamedFaceGroup {
                name,
                faces: group_faces,
            });
            cursor.expect("]}")?;
            match cursor.peek() {
                Some(b',') => cursor.expect(",")?,
                Some(b']') => break,
                _ => {
                    return Err(cursor.problem("expected `,` or `]` after import named group"));
                }
            }
        }
    }
    cursor.expect("]")?;
    let mut group_order: Vec<usize> = (0..named_groups.len()).collect();
    group_order
        .sort_unstable_by(|left, right| named_groups[*left].name.cmp(&named_groups[*right].name));
    if let Some(duplicate) = group_order
        .windows(2)
        .find(|window| named_groups[window[0]].name == named_groups[window[1]].name)
    {
        return Err(format!(
            "import IR source {source_index} repeats named group `{}`",
            named_groups[duplicate[0]].name
        ));
    }
    Ok(named_groups)
}

fn parse_geometry_import_summary(
    text: &str,
    spec: &ProjectSpec,
    project_hash: ContentHash,
) -> Result<Vec<VerifiedImport>, String> {
    let mut cursor = JsonCursor::new(text);
    cursor.expect("{\"schema\":")?;
    let schema = cursor.parse_string()?;
    if schema != IMPORT_SUMMARY_SCHEMA {
        return Err(format!(
            "summary schema is `{schema}`, not `{IMPORT_SUMMARY_SCHEMA}`"
        ));
    }
    cursor.expect(",\"project_hash\":")?;
    let retained_project = cursor.parse_string()?;
    if retained_project != project_hash.to_hex() {
        return Err(format!(
            "summary project hash `{retained_project}` does not match `{}`",
            project_hash.to_hex()
        ));
    }
    cursor.expect(",\"artifacts\":[")?;
    let mut entries = Vec::new();
    if cursor.peek() != Some(b']') {
        loop {
            entries.push(parse_geometry_import_entry(&mut cursor, project_hash)?);
            match cursor.peek() {
                Some(b',') => cursor.expect(",")?,
                Some(b']') => break,
                _ => return Err(cursor.problem("expected `,` or `]` after import artifact")),
            }
        }
    }
    cursor.expect("]")?;
    cursor.expect(",\"assignment_table\":")?;
    let _assignment_table = cursor.parse_string()?;
    cursor.expect(",\"authority\":")?;
    let authority = cursor.parse_string()?;
    if authority != IMPORT_SUMMARY_AUTHORITY {
        return Err(format!(
            "summary authority is `{authority}`, not `{IMPORT_SUMMARY_AUTHORITY}`"
        ));
    }
    cursor.expect(",\"no_claim\":")?;
    let no_claim = cursor.parse_string()?;
    if no_claim != IMPORT_SUMMARY_NO_CLAIM {
        return Err("summary no-claim boundary does not match the versioned writer".to_string());
    }
    cursor.expect("}")?;
    cursor.finish()?;

    let geometry = spec.geometry.as_deref().unwrap_or(&[]);
    if entries.len() != geometry.len() {
        return Err(format!(
            "summary carries {} artifacts but the project declares {}",
            entries.len(),
            geometry.len()
        ));
    }
    for (index, (entry, artifact)) in entries.iter().zip(geometry).enumerate() {
        let expected_identity = geometry_source_identity(artifact);
        if entry.role != artifact.role || entry.source_identity != expected_identity {
            return Err(format!(
                "summary artifact {index} does not match project role `{}` and source identity `{expected_identity}`",
                artifact.role
            ));
        }
    }
    Ok(entries)
}

fn parse_geometry_import_entry(
    cursor: &mut JsonCursor<'_>,
    project_hash: ContentHash,
) -> Result<VerifiedImport, String> {
    cursor.expect("{\"role\":")?;
    let role = cursor.parse_string()?;
    cursor.expect(",\"source_label\":")?;
    let source_label = cursor.parse_string()?;
    if source_label.is_empty()
        || source_label.len() > 4096
        || source_label.chars().any(char::is_control)
    {
        return Err("summary source label violates the import writer bound".to_string());
    }
    cursor.expect(",\"source_label_authority\":")?;
    let source_label_authority = cursor.parse_string()?;
    if source_label_authority != IMPORT_SOURCE_LABEL_AUTHORITY {
        return Err("summary source-label authority is not caller-reported".to_string());
    }
    cursor.expect(",\"source_identity\":")?;
    let source_identity = cursor.parse_string()?;
    cursor.expect(",\"raw_source\":")?;
    let raw_source = parse_hash_string(cursor, "raw_source")?;
    cursor.expect(",\"promotion_receipt\":")?;
    let promotion_receipt = parse_hash_string(cursor, "promotion_receipt")?;
    cursor.expect(",\"promoted_mesh\":")?;
    let promoted_mesh = parse_hash_string(cursor, "promoted_mesh")?;
    cursor.expect(",\"assignment_report\":")?;
    let assignment_report = parse_hash_string(cursor, "assignment_report")?;
    cursor.expect(",\"import_record\":")?;
    let import_record = cursor.parse_string()?;
    let expected_record = format!("{}:{source_identity}", project_hash.to_hex());
    if import_record != expected_record {
        return Err(format!(
            "summary import record `{import_record}` does not match `{expected_record}`"
        ));
    }
    cursor.expect("}")?;
    Ok(VerifiedImport {
        role,
        source_identity,
        raw_source,
        promotion_receipt,
        promoted_mesh,
        assignment_report,
    })
}

fn parse_import_verify_receipt(
    text: &str,
    run: SolveRunId,
    project_hash: ContentHash,
) -> Result<(i64, Vec<VerifiedImport>), String> {
    let mut cursor = JsonCursor::new(text);
    cursor.expect("{\"schema\":")?;
    let schema = cursor.parse_string()?;
    if schema != IMPORT_VERIFY_RECEIPT_SCHEMA {
        return Err(format!(
            "receipt schema is `{schema}`, not `{IMPORT_VERIFY_RECEIPT_SCHEMA}`"
        ));
    }
    cursor.expect(",\"run\":")?;
    let retained_run = cursor.parse_string()?;
    if retained_run != run.to_hex() {
        return Err(format!(
            "receipt run `{retained_run}` does not match `{}`",
            run.to_hex()
        ));
    }
    cursor.expect(",\"project_hash\":")?;
    let retained_project = cursor.parse_string()?;
    if retained_project != project_hash.to_hex() {
        return Err(format!(
            "receipt project `{retained_project}` does not match `{}`",
            project_hash.to_hex()
        ));
    }
    cursor.expect(",\"import_op\":")?;
    let import_op = cursor.parse_i64()?;
    if import_op <= 0 {
        return Err("receipt import operation id is not positive".to_string());
    }
    cursor.expect(",\"verified\":[")?;
    let mut entries = Vec::new();
    if cursor.peek() != Some(b']') {
        loop {
            entries.push(parse_verified_import_entry(&mut cursor)?);
            match cursor.peek() {
                Some(b',') => cursor.expect(",")?,
                Some(b']') => break,
                _ => return Err(cursor.problem("expected `,` or `]` after verified import")),
            }
        }
    }
    if entries.is_empty() {
        return Err("receipt verified-import array is empty".to_string());
    }
    cursor.expect("]")?;
    cursor.expect(",\"authority\":")?;
    let authority = cursor.parse_string()?;
    if authority != IMPORT_VERIFY_AUTHORITY {
        return Err("receipt authority does not match the versioned writer".to_string());
    }
    cursor.expect(",\"no_claim\":")?;
    let no_claim = cursor.parse_string()?;
    if no_claim != IMPORT_VERIFY_NO_CLAIM {
        return Err("receipt no-claim boundary does not match the versioned writer".to_string());
    }
    cursor.expect("}")?;
    cursor.finish()?;
    Ok((import_op, entries))
}

fn parse_verified_import_entry(cursor: &mut JsonCursor<'_>) -> Result<VerifiedImport, String> {
    cursor.expect("{\"role\":")?;
    let role = cursor.parse_string()?;
    cursor.expect(",\"source_identity\":")?;
    let source_identity = cursor.parse_string()?;
    cursor.expect(",\"raw_source\":")?;
    let raw_source = parse_hash_string(cursor, "raw_source")?;
    cursor.expect(",\"promotion_receipt\":")?;
    let promotion_receipt = parse_hash_string(cursor, "promotion_receipt")?;
    cursor.expect(",\"promoted_mesh\":")?;
    let promoted_mesh = parse_hash_string(cursor, "promoted_mesh")?;
    cursor.expect(",\"assignment_report\":")?;
    let assignment_report = parse_hash_string(cursor, "assignment_report")?;
    cursor.expect("}")?;
    Ok(VerifiedImport {
        role,
        source_identity,
        raw_source,
        promotion_receipt,
        promoted_mesh,
        assignment_report,
    })
}

fn parse_hash_string(cursor: &mut JsonCursor<'_>, field: &str) -> Result<ContentHash, String> {
    let value = cursor.parse_string()?;
    let hash = ContentHash::from_hex(&value)
        .ok_or_else(|| format!("field `{field}` is not a 64-hex content hash"))?;
    if value != hash.to_hex() {
        return Err(format!(
            "field `{field}` is not the canonical lowercase hash spelling"
        ));
    }
    Ok(hash)
}

struct JsonCursor<'a> {
    input: &'a str,
    pos: usize,
}

impl<'a> JsonCursor<'a> {
    fn new(input: &'a str) -> Self {
        Self { input, pos: 0 }
    }

    fn peek(&self) -> Option<u8> {
        self.input.as_bytes().get(self.pos).copied()
    }

    fn problem(&self, what: &str) -> String {
        format!("{what} at byte {}", self.pos)
    }

    fn expect(&mut self, expected: &str) -> Result<(), String> {
        if self.input[self.pos..].starts_with(expected) {
            self.pos += expected.len();
            Ok(())
        } else {
            Err(self.problem(&format!("expected `{expected}`")))
        }
    }

    fn finish(&self) -> Result<(), String> {
        if self.pos == self.input.len() {
            Ok(())
        } else {
            Err(self.problem("trailing bytes after complete JSON value"))
        }
    }

    fn parse_string(&mut self) -> Result<String, String> {
        self.expect("\"")?;
        let mut value = String::new();
        loop {
            let character = self.input[self.pos..]
                .chars()
                .next()
                .ok_or_else(|| self.problem("unterminated JSON string"))?;
            self.pos += character.len_utf8();
            match character {
                '"' => return Ok(value),
                '\\' => {
                    let escape = self.input[self.pos..]
                        .chars()
                        .next()
                        .ok_or_else(|| self.problem("unterminated JSON escape"))?;
                    self.pos += escape.len_utf8();
                    match escape {
                        '"' => value.push('"'),
                        '\\' => value.push('\\'),
                        'n' => value.push('\n'),
                        'r' => value.push('\r'),
                        't' => value.push('\t'),
                        'u' => {
                            let scalar = self.parse_hex_quad()?;
                            if scalar > 0x1f || matches!(scalar, 0x09 | 0x0a | 0x0d) {
                                return Err(self.problem(
                                    "non-canonical Unicode escape; writer emits raw scalars or named whitespace escapes",
                                ));
                            }
                            value.push(
                                char::from_u32(u32::from(scalar))
                                    .expect("control code is a valid Unicode scalar"),
                            );
                        }
                        _ => return Err(self.problem("non-canonical JSON escape")),
                    }
                }
                character if character <= '\u{1f}' => {
                    return Err(self.problem("unescaped control character in JSON string"));
                }
                character => value.push(character),
            }
        }
    }

    fn parse_hex_quad(&mut self) -> Result<u16, String> {
        let end = self
            .pos
            .checked_add(4)
            .ok_or_else(|| self.problem("Unicode escape offset overflow"))?;
        let digits = self
            .input
            .get(self.pos..end)
            .ok_or_else(|| self.problem("short Unicode escape"))?;
        if !digits
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return Err(self.problem("non-lowercase-hex Unicode escape"));
        }
        self.pos = end;
        u16::from_str_radix(digits, 16).map_err(|_| self.problem("invalid Unicode escape"))
    }

    fn parse_i64(&mut self) -> Result<i64, String> {
        let start = self.pos;
        if self.peek() == Some(b'-') {
            self.pos += 1;
        }
        match self.peek() {
            Some(b'0') => {
                self.pos += 1;
                if self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    return Err(self.problem("leading zero in JSON integer"));
                }
            }
            Some(b'1'..=b'9') => {
                self.pos += 1;
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.pos += 1;
                }
            }
            _ => return Err(self.problem("expected JSON integer")),
        }
        self.input[start..self.pos]
            .parse::<i64>()
            .map_err(|_| self.problem("JSON integer is outside i64"))
    }

    fn parse_u64(&mut self) -> Result<u64, String> {
        let start = self.pos;
        match self.peek() {
            Some(b'0') => {
                self.pos += 1;
                if self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    return Err(self.problem("leading zero in JSON unsigned integer"));
                }
            }
            Some(b'1'..=b'9') => {
                self.pos += 1;
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.pos += 1;
                }
            }
            _ => return Err(self.problem("expected JSON unsigned integer")),
        }
        self.input[start..self.pos]
            .parse::<u64>()
            .map_err(|_| self.problem("JSON unsigned integer is outside u64"))
    }

    fn parse_usize(&mut self) -> Result<usize, String> {
        let value = self.parse_u64()?;
        usize::try_from(value).map_err(|_| self.problem("JSON unsigned integer is outside usize"))
    }

    fn parse_canonical_finite_f64(&mut self) -> Result<f64, String> {
        let start = self.pos;
        self.skip_number()?;
        let spelling = &self.input[start..self.pos];
        let value = spelling
            .parse::<f64>()
            .map_err(|_| self.problem("JSON number is outside f64"))?;
        if !value.is_finite() {
            return Err(self.problem("JSON number is not finite"));
        }
        if value.to_string() != spelling {
            return Err(self.problem("JSON number is not the canonical Rust writer spelling"));
        }
        Ok(value)
    }

    fn take_value(&mut self) -> Result<&'a str, String> {
        let start = self.pos;
        self.skip_value(0)?;
        Ok(&self.input[start..self.pos])
    }

    fn skip_value(&mut self, depth: usize) -> Result<(), String> {
        if depth > 64 {
            return Err(self.problem("JSON nesting exceeds 64"));
        }
        match self.peek() {
            Some(b'"') => {
                let _ = self.parse_string()?;
                Ok(())
            }
            Some(b'{') => {
                self.pos += 1;
                if self.peek() == Some(b'}') {
                    self.pos += 1;
                    return Ok(());
                }
                loop {
                    let _ = self.parse_string()?;
                    self.expect(":")?;
                    self.skip_value(depth + 1)?;
                    match self.peek() {
                        Some(b',') => self.pos += 1,
                        Some(b'}') => {
                            self.pos += 1;
                            return Ok(());
                        }
                        _ => return Err(self.problem("expected `,` or `}` in JSON object")),
                    }
                }
            }
            Some(b'[') => {
                self.pos += 1;
                if self.peek() == Some(b']') {
                    self.pos += 1;
                    return Ok(());
                }
                loop {
                    self.skip_value(depth + 1)?;
                    match self.peek() {
                        Some(b',') => self.pos += 1,
                        Some(b']') => {
                            self.pos += 1;
                            return Ok(());
                        }
                        _ => return Err(self.problem("expected `,` or `]` in JSON array")),
                    }
                }
            }
            Some(b'-' | b'0'..=b'9') => self.skip_number(),
            Some(b't') => self.expect("true"),
            Some(b'f') => self.expect("false"),
            Some(b'n') => self.expect("null"),
            _ => Err(self.problem("expected JSON value")),
        }
    }

    fn skip_number(&mut self) -> Result<(), String> {
        if self.peek() == Some(b'-') {
            self.pos += 1;
        }
        match self.peek() {
            Some(b'0') => {
                self.pos += 1;
                if self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    return Err(self.problem("leading zero in JSON number"));
                }
            }
            Some(b'1'..=b'9') => {
                self.pos += 1;
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.pos += 1;
                }
            }
            _ => return Err(self.problem("invalid JSON number")),
        }
        if self.peek() == Some(b'.') {
            self.pos += 1;
            let fraction_start = self.pos;
            while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                self.pos += 1;
            }
            if self.pos == fraction_start {
                return Err(self.problem("empty JSON number fraction"));
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.pos += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.pos += 1;
            }
            let exponent_start = self.pos;
            while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                self.pos += 1;
            }
            if self.pos == exponent_start {
                return Err(self.problem("empty JSON number exponent"));
            }
        }
        Ok(())
    }
}
