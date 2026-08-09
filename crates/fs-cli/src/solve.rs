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
//! Executing boundary (stated, not implied): `import-verify` and `assign`
//! execute against retained import evidence, and `material-resolve` executes
//! against caller-supplied normalized card packs (bead frankensim-hp7tb);
//! `flow-network` (frankensim-frn2i), `conduction` (frankensim-s93ej), and
//! `qoi` (frankensim-s2l9v) remain typed gaps.
//!
//! Card packs are invocation inputs, so their canonical set root is bound
//! into the run identity: a different pack set is a different run, never the
//! same run with a different answer. The exact pack bytes are retained
//! against the run's first operation alongside the project source, so resume
//! recovers them from the ledger no matter which stage was interrupted, and
//! re-deriving the run identity from the recovered set is what makes that
//! recovery an attestation rather than a restatement.

// Refusals are cold-path values carrying complete diagnostics; the crate's
// refusal idiom (`GeometryImportRefusal`) is by-value for the same reason.
#![allow(clippy::result_large_err)]

use std::cell::Cell;
use std::collections::HashSet;
use std::collections::hash_map::DefaultHasher;
use std::fmt::Write as _;
use std::hash::BuildHasherDefault;
use std::ops::ControlFlow;
use std::rc::Rc;

use fs_blake3::identity::ContentId;
use fs_blake3::{hash_bytes, hash_domain};
use fs_exec::CancelGate;
use fs_exec::solver::{
    LegacySnapshotExpectationV1, LegacySnapshotLimitsV1, LegacySnapshotV1Adapter,
    LegacySnapshotV1Error, LegacySolverStateV1, codec,
};
use fs_io::{MESH_ASSIGNMENT_SEMANTICS_VERSION, MeshSelector, NamedFaceGroup};
use fs_ledger::{
    ArtifactInfo, BoundedOpArtifactEdges, CONTROLLED_ARTIFACT_TILE_LEN, ContentHash,
    ControlledOpRead, ControlledVisibleOpPage, EdgeRole, ExecMode, FiveExplicits, Ledger,
    LedgerError, MAIN_BRANCH, MAX_VISIBLE_OP_PAGE_ROWS, OpArtifactEdge, OpOutcome, OpRow,
    OpVariableField, PrehashedOpContent, VisibleOpCursor, VisibleOpPage,
};
use fs_project::{
    BindingRequirements, DecodedProject, GeometryArtifact, ProjectSpec, geometry_source_identity,
    resolve_bindings,
};
use fs_session::{CapabilityToken, Charge, Enforcement, Governor, SessionError, SessionId};

use crate::cards::{CardPackKind, CardPackSet, CardPackSetBuilder, RawCardPack};
use crate::import::{explicits, json_string};

/// Domain separating solve-run identity derivation from every other hash.
pub const SOLVE_RUN_IDENTITY_DOMAIN: &str = "org.frankensim.fs-cli.solve-run.v1";
/// Driver semantics version bound into run identity and driver state.
///
/// Bumped to 3 when `material-resolve` stopped being a typed gap: the stage
/// set a run identity stands for changed, and the run preimage gained the
/// admitted card-pack-set root. A v2 checkpoint therefore cannot be mistaken
/// for a v3 run. Bumped to 4 when `flow-network` stopped being a typed gap
/// (frankensim-frn2i.2): v3 checkpoints cannot resume into a v4 stage set.
pub const SOLVE_DRIVER_VERSION: u32 = 4;

const SOLVE_STAGE_SCHEMA: &str = "frankensim.cli.solve-stage.v1";
const SOLVE_RUN_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-run-receipt.v1";
const IMPORT_VERIFY_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-import-verify-receipt.v1";
const ASSIGN_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-assignment-binding.v1";
const MATERIAL_RESOLVE_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-material-resolve-receipt.v1";
const IMPORT_IR_SCHEMA: &str = "frankensim.cli.geometry-import.v1";
const IMPORT_SUMMARY_SCHEMA: &str = "frankensim.cli.geometry-import-receipt.v1";

const PROJECT_SOURCE_KIND: &str = "solve-project-source";
const STAGE_STATE_KIND: &str = "solve-stage-state";
const STAGE_RECEIPT_KIND: &str = "solve-stage-receipt";
const RUN_RECEIPT_KIND: &str = "solve-run-receipt";
const MATERIAL_USAGE_KIND: &str = "solve-material-usage-receipt";
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
const MATERIAL_RESOLVE_AUTHORITY: &str = "declared-binding-resolution-against-admitted-card-packs";
const FLOW_NETWORK_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-flow-network-receipt.v1";
/// Specific gas constant of dry air, J/(kg·K); used only for the declared
/// envelope-derived density estimate (see FLOW_NETWORK_NO_CLAIM).
const AIR_SPECIFIC_GAS_CONSTANT: f64 = 287.05;
const FLOW_NETWORK_AUTHORITY: &str =
    "lossless-project-lowering-plus-interval-certified-operating-point";
const FLOW_NETWORK_NO_CLAIM: &str = "the stage proves the declared fan system lowered losslessly \
    and the enclosure network produced an interval-certified nominal operating point under the \
    declared orifice/leakage models; it does not authenticate manufacturer curve data, system \
    effects, compressibility, installation effects, or any experimental validation, and the \
    envelope-derived air density is a declared ideal-gas estimate, not a measurement";
const MATERIAL_RESOLVE_NO_CLAIM: &str = "the stage proves that every declared region and \
    interface resolves to an admitted card whose selected claim covers the declared temperature \
    range, and retains that claim's replayable usage receipt; it does not authenticate the pack \
    producer, validate the claim against any external corpus, narrow or replace a claim's stated \
    uncertainty, or upgrade an Unstated uncertainty into a bound";
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
/// Maximum solve-owned byte work between evidence cancellation observations.
const EVIDENCE_POLL_BYTES: usize = CONTROLLED_ARTIFACT_TILE_LEN;
/// Longest canonical signed/unsigned 64-bit JSON integer spelling.
const MAX_JSON_INTEGER_BYTES: usize = 20;
/// Longest finite Rust `f64::to_string()` spelling. The negative smallest
/// subnormal needs 327 bytes because finite `Display` uses fixed decimal text
/// at that magnitude.
const MAX_CANONICAL_F64_BYTES: usize = 327;
/// Longest canonical decimal spelling of a `u32`.
const MAX_CANONICAL_U32_BYTES: usize = 10;
/// Separator-inclusive bound for one canonical writer vertex line.
const MAX_CANONICAL_PLY_VERTEX_LINE_BYTES: usize = MAX_CANONICAL_F64_BYTES * 3 + 2;
/// Separator-inclusive bound for one canonical writer triangle line.
const MAX_CANONICAL_PLY_FACE_LINE_BYTES: usize = MAX_CANONICAL_U32_BYTES * 3 + 4;
/// Solve-local ceiling for retained names and identities that are compared or
/// hashed outside the JSON cursor. This retains the lower-layer 4-KiB default
/// and keeps every individual scalar operation below one evidence tile.
const MAX_SOLVE_EVIDENCE_LABEL_BYTES: usize = 4096;
type SolveEvidenceSet<T> = HashSet<T, BuildHasherDefault<DefaultHasher>>;

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
/// proof: evidence materialization, incremental UTF-8 copies, and parser-owned
/// allocations can coexist within one candidate verification.
const MAX_TOTAL_SOLVE_EVIDENCE_BYTES: u64 = 512 * 1024 * 1024;
/// Invocation-wide cumulative work envelope for discovery and retained
/// evidence re-attestation. This deterministic conservative byte-equivalent
/// measure charges controlled input bytes once, accepted copy/comparison and
/// derived-output bytes, plus fixed proxies for ids/items. It is not exact CPU
/// accounting. The separate visible-id ceiling prevents a long history of tiny
/// rows from avoiding this byte envelope.
const MAX_SOLVE_INVOCATION_WORK_BYTES: u64 = 1024 * 1024 * 1024;
/// Maximum visible operation IDs examined by one solve invocation.
#[doc(hidden)]
pub const MAX_SOLVE_VISIBLE_OP_IDS: usize = 8192;
/// Maximum fixed-size descending pages consumed before an explicit refusal.
const MAX_SOLVE_VISIBLE_OP_PAGES: usize =
    MAX_SOLVE_VISIBLE_OP_IDS.div_ceil(MAX_VISIBLE_OP_PAGE_ROWS);
/// Byte-equivalent work charged for one fixed-size operation id.
const OP_ID_WORK_BYTES: u64 = 8;
/// Fixed byte-equivalent work charged for one typed operation sidecar.
const OP_CONTENT_IDENTITY_WORK_BYTES: u64 = 256;
/// Byte-equivalent work charged for one role-qualified content-hash edge.
const EDGE_ITEM_WORK_BYTES: u64 = 33;
/// Byte-equivalent work charged for one bounded project/assignment item.
const DERIVATION_ITEM_WORK_BYTES: u64 = 64;
/// Edge scan cap per operation while locating retained evidence.
const EDGE_SCAN_CAP: usize = 1024;
/// Largest distinct claim-usage receipt set one `material-resolve` operation
/// retains. Together with the card-pack ceiling this keeps the stage's typed
/// edge set inside the ledger's bounded scan.
const MAX_MATERIAL_USAGE_RECEIPTS: usize = 512;
/// Read cap for one retained card pack recovered during resume. This mirrors
/// the admission ceiling so a pack that was admissible on the fresh path
/// cannot become unreadable on the resume path.
const MAX_CARD_PACK_READ_BYTES: u64 = crate::cards::MAX_CARD_PACK_BYTES;
/// Largest geometry set the solve evidence contract can attest completely.
/// One import operation has at most `4 * sources + 1` typed edges, so 255
/// sources fit under the ledger's 1024-edge bounded scan.
const SOLVE_MAX_IMPORT_SOURCES: usize = 255;

/// Stable solve evidence-verification phases exposed only for deterministic
/// conformance-test cancellation plans.
///
/// Production callers request cancellation through [`CancelGate`] directly.
/// This enum does not grant a callback or any publication authority.
#[doc(hidden)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveEvidencePhase {
    /// Read one descending, frozen-high-water page of visible operation ids.
    VisibleOpPage,
    /// Retrieve one bounded candidate operation row.
    CandidateOpRowRead,
    /// Convert one guarded candidate text field after controlled delivery.
    CandidateOpTextConversion,
    /// Recompute one candidate operation's typed content identity.
    OperationContentIdentity,
    /// Render the project-derived Five Explicits.
    FiveExplicitsRender,
    /// Compare retained and re-rendered Five Explicits in bounded tiles.
    FiveExplicitsCompare,
    /// Render the canonical project JSON once for import-row comparison.
    CanonicalProjectRender,
    /// Validate one decoded project through its lower-layer validator.
    ProjectValidation,
    /// Derive the canonical project and solve-run identities.
    ProjectIdentityDerive,
    /// Re-derive the project's stable entity identities.
    EntityResolution,
    /// Derive assignment identities, rows, or the canonical assignment receipt.
    AssignmentDerivation,
    /// Resolve declared bindings against the admitted card library, or build
    /// the canonical material-resolve receipt.
    MaterialBindingResolution,
    /// Lower the declared cooling fan system and solve the enclosure's
    /// network operating point, or build the canonical flow-network receipt.
    FlowNetworkSolve,
    /// Materialize one retained card pack during resume.
    ///
    /// The optional plan index is the pack's canonical position in the set.
    ResumeCardPackRead,
    /// Decode and re-admit one recovered card pack during resume.
    ResumeCardPackDecode,
    /// Construct one canonical stage receipt from already verified evidence.
    ReceiptDerivation,
    /// Retrieve one bounded operation edge page.
    EdgePageRead,
    /// Retrieve one bounded artifact descriptor used by discovery or re-attestation.
    ArtifactDescriptorRead,
    /// Retrieve one fixed operation-edge lineage seal.
    EdgeSealRead,
    /// Compare one bounded retained edge set with its exact expected set.
    EdgeSetCompare,
    /// Parse the retained geometry-import operation IR.
    ImportIrParse,
    /// Compare the import IR's embedded project with a canonical rendering.
    ImportIrCanonicalCompare,
    /// Check import-IR named-group face and name uniqueness without sorting.
    ImportIrDuplicateCheck,
    /// Materialize the retained geometry-import summary.
    ImportSummaryRead,
    /// Incrementally validate and copy the retained import summary as UTF-8.
    ImportSummaryUtf8,
    /// Parse the retained geometry-import summary.
    ImportSummaryParse,
    /// Verify one retained raw source without materializing it.
    RawSourceRead,
    /// Verify one opaque lower-layer promotion receipt.
    PromotionReceiptRead,
    /// Materialize one promoted canonical PLY.
    PromotedMeshRead,
    /// Check one promoted PLY's bounded canonical header and body shape.
    PromotedMeshPreflight,
    /// Decode and range-check one canonical promoted PLY payload.
    PromotedMeshDecode,
    /// Re-scan one promoted PLY for exact fs-io writer token spellings.
    PromotedMeshEncodeCompare,
    /// Range-check retained named-group face references in bounded tiles.
    NamedGroupFaceRange,
    /// Materialize one assignment report.
    AssignmentReportRead,
    /// Incrementally validate and copy one assignment report as UTF-8.
    AssignmentReportUtf8,
    /// Parse and re-attest one assignment report.
    AssignmentReportParse,
    /// Materialize the retained canonical project during resume.
    ResumeProjectRead,
    /// Incrementally validate and copy the retained canonical project.
    ResumeProjectUtf8,
    /// Strictly parse the retained project in one bounded opaque call.
    ResumeProjectParse,
    /// Compare the parser's canonical project with retained bytes in tiles.
    ResumeProjectCanonicalCompare,
    /// Materialize a retained driver checkpoint during resume.
    ResumeStateRead,
    /// Decode and inspect a retained driver checkpoint envelope.
    ResumeStateDecode,
    /// Materialize a retained stage receipt during resume.
    ResumeStageReceiptRead,
    /// Incrementally validate and copy a retained stage receipt as UTF-8.
    ResumeStageReceiptUtf8,
    /// Parse a retained import-verification receipt during resume.
    ResumeStageReceiptParse,
    /// Compare one reconstructed canonical stage receipt with retained bytes.
    ///
    /// The optional plan index is the completed-stage index for this phase.
    ResumeStageReceiptCanonicalCompare,
    /// Final gate observation after a successful stage body and before any
    /// clock, charge, or ledger publication for that stage.
    PrePublication,
}

/// Deterministic conformance-test plan that requests the supplied gate at one
/// typed evidence-work checkpoint.
///
/// `after_units` is a phase-relative lower bound. Byte phases observe it at
/// fixed boundaries no wider than [`EVIDENCE_POLL_BYTES`]; item and opaque-call
/// phases observe it before and after each owned item/call. A plan can only
/// request the gate once and cannot run caller code.
#[doc(hidden)]
#[derive(Debug)]
pub struct SolveCancellationPlan {
    phase: SolveEvidencePhase,
    /// Optional phase-defined item index: geometry-source index for import
    /// evidence, or completed-stage index for canonical stage receipts.
    source_index: Option<usize>,
    after_units: u64,
    fired: Cell<bool>,
}

impl SolveCancellationPlan {
    /// Construct a one-shot deterministic cancellation request.
    #[doc(hidden)]
    #[must_use]
    pub const fn new(
        phase: SolveEvidencePhase,
        source_index: Option<usize>,
        after_units: u64,
    ) -> Self {
        Self {
            phase,
            source_index,
            after_units,
            fired: Cell::new(false),
        }
    }

    /// Whether this plan reached its requested checkpoint.
    #[doc(hidden)]
    #[must_use]
    pub fn fired(&self) -> bool {
        self.fired.get()
    }

    fn observe(
        &self,
        gate: &CancelGate,
        phase: SolveEvidencePhase,
        source_index: Option<usize>,
        units: u64,
    ) {
        if !self.fired.get()
            && self.phase == phase
            && self.source_index == source_index
            && units >= self.after_units
        {
            self.fired.set(true);
            gate.request();
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct EvidenceCancelled;

#[derive(Debug, Clone, Copy)]
enum InvocationWorkExceeded {
    CumulativeBytes { attempted: u64 },
    VisibleIds { limit: usize },
    VisiblePages { limit: usize },
}

#[derive(Debug, Default)]
struct InvocationWorkLedger {
    used: Cell<u64>,
}

impl InvocationWorkLedger {
    fn charge(&self, bytes: u64) -> Result<u64, InvocationWorkExceeded> {
        let attempted =
            self.used
                .get()
                .checked_add(bytes)
                .ok_or(InvocationWorkExceeded::CumulativeBytes {
                    attempted: u64::MAX,
                })?;
        if attempted > MAX_SOLVE_INVOCATION_WORK_BYTES {
            return Err(InvocationWorkExceeded::CumulativeBytes { attempted });
        }
        self.used.set(attempted);
        Ok(attempted)
    }
}

#[derive(Clone, Copy)]
struct EvidenceWork<'a> {
    gate: &'a CancelGate,
    plan: Option<&'a SolveCancellationPlan>,
    invocation: Option<&'a InvocationWorkLedger>,
}

impl<'a> EvidenceWork<'a> {
    const fn new(
        gate: &'a CancelGate,
        plan: Option<&'a SolveCancellationPlan>,
        invocation: &'a InvocationWorkLedger,
    ) -> Self {
        Self {
            gate,
            plan,
            invocation: Some(invocation),
        }
    }

    #[cfg(test)]
    const fn unmetered(gate: &'a CancelGate, plan: Option<&'a SolveCancellationPlan>) -> Self {
        Self {
            gate,
            plan,
            invocation: None,
        }
    }

    fn checkpoint(
        self,
        phase: SolveEvidencePhase,
        source_index: Option<usize>,
        units: u64,
    ) -> Result<(), EvidenceCancelled> {
        if let Some(plan) = self.plan {
            plan.observe(self.gate, phase, source_index, units);
        }
        if self.gate.is_requested() {
            Err(EvidenceCancelled)
        } else {
            Ok(())
        }
    }

    fn is_requested(self) -> bool {
        self.gate.is_requested()
    }

    fn charge(self, bytes: u64) -> Result<u64, InvocationWorkExceeded> {
        self.invocation
            .map_or(Ok(bytes), |meter| meter.charge(bytes))
    }
}

#[derive(Debug)]
enum EvidenceReadError {
    Cancelled,
    WorkEnvelope(InvocationWorkExceeded),
    Ledger(LedgerError),
}

#[derive(Debug)]
enum EvidenceUtf8Error {
    Cancelled,
    WorkEnvelope(InvocationWorkExceeded),
    Invalid(String),
}

#[derive(Debug)]
enum EvidenceCompareError {
    Cancelled,
    WorkEnvelope(InvocationWorkExceeded),
}

#[derive(Debug)]
enum ProjectRenderError {
    Cancelled,
    WorkEnvelope(InvocationWorkExceeded),
    Invalid(String),
}

#[derive(Debug)]
enum CandidateOpReadError {
    Cancelled,
    WorkEnvelope(InvocationWorkExceeded),
    Ledger(LedgerError),
}

#[derive(Debug)]
enum VisibleOpScanError {
    Cancelled,
    WorkEnvelope(InvocationWorkExceeded),
    Ledger(LedgerError),
}

#[derive(Default)]
struct CandidateOpFields {
    session: Vec<u8>,
    ir: Vec<u8>,
    seed: Vec<u8>,
    versions: Vec<u8>,
    budget: Vec<u8>,
    capability: Vec<u8>,
    diagnostic: Vec<u8>,
}

impl CandidateOpFields {
    fn field_mut(&mut self, field: OpVariableField) -> &mut Vec<u8> {
        match field {
            OpVariableField::Session => &mut self.session,
            OpVariableField::Ir => &mut self.ir,
            OpVariableField::Seed => &mut self.seed,
            OpVariableField::Versions => &mut self.versions,
            OpVariableField::Budget => &mut self.budget,
            OpVariableField::Capability => &mut self.capability,
            OpVariableField::Diagnostic => &mut self.diagnostic,
        }
    }
}

struct ControlledCandidateOp {
    row: OpRow,
    branch: i64,
    exec_mode: ExecMode,
    prehashed_content: PrehashedOpContent,
}

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
    /// Resolve declared material/interface bindings against the admitted
    /// card packs (built at frankensim-hp7tb).
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
            SolveStage::ImportVerify | SolveStage::Assign | SolveStage::MaterialResolve => None,
            SolveStage::FlowNetwork => None,
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
    /// Derive the run identity from the decoded project and the admitted
    /// card-pack set.
    ///
    /// The pack set is part of the run's answer, not an incidental
    /// invocation detail: two runs that bind different cards must not share
    /// one identity, because the driver would then see two competing
    /// checkpoints for what it believes is one run.
    #[must_use]
    pub fn derive(project: &DecodedProject, cards: &CardPackSet) -> SolveRunId {
        let project_hash = project.hash();
        Self::derive_with_project_hash(project, project_hash, cards.root())
    }

    fn derive_with_project_hash(
        project: &DecodedProject,
        project_hash: ContentHash,
        cards_root: ContentHash,
    ) -> SolveRunId {
        let spec = &project.spec;
        let mut preimage = Vec::with_capacity(160);
        preimage.extend_from_slice(project_hash.as_bytes());
        let (constellation, workspace) = spec.versions.as_ref().map_or(("", ""), |v| {
            (v.constellation.as_str(), v.workspace.as_str())
        });
        push_framed(&mut preimage, constellation.as_bytes());
        push_framed(&mut preimage, workspace.as_bytes());
        let seed = spec.seeds.as_ref().map_or(0, |s| s.root);
        preimage.extend_from_slice(&seed.to_le_bytes());
        push_framed(&mut preimage, cards_root.as_bytes());
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

#[derive(Clone, Copy)]
struct RenderedExplicitsRef<'a> {
    seed: &'a [u8],
    versions: &'a str,
    budget: &'a str,
    capability: &'a str,
    canonical_project_json: &'a str,
}

#[derive(Clone, Copy)]
struct ImportCandidateLocator {
    op: i64,
    summary_artifact: ContentHash,
}

#[derive(Clone, Copy)]
struct StageRowExpectation<'a> {
    stage: SolveStage,
    run: SolveRunId,
    project_hash: ContentHash,
    explicits: RenderedExplicitsRef<'a>,
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
    attestation: Rc<ResumeImportCache>,
    context: StageContext,
    /// Card packs recovered from the run's retained inputs and re-admitted
    /// through the same canonicalization the fresh path uses. The recovered
    /// set root must reproduce the requested run identity.
    cards: CardPackSet,
}

#[derive(Debug)]
struct ResumeImportCache {
    source_hash: ContentHash,
    project: DecodedProject,
    project_hash: ContentHash,
    versions: String,
    budget: String,
    capability: String,
    seed: [u8; 8],
    canonical_project_json: String,
}

impl ResumeImportCache {
    fn expectations(&self) -> RenderedExplicitsRef<'_> {
        RenderedExplicitsRef {
            seed: &self.seed,
            versions: &self.versions,
            budget: &self.budget,
            capability: &self.capability,
            canonical_project_json: &self.canonical_project_json,
        }
    }
}

/// Everything a running or resumed solve needs in one place.
struct SolveEngine<'a> {
    ledger: &'a Ledger,
    work: EvidenceWork<'a>,
    clock: &'a mut dyn FnMut() -> f64,
    spec: &'a ProjectSpec,
    canonical_source: &'a str,
    /// The admitted card packs backing this run. On a fresh run these are the
    /// caller's inputs; on resume they are recovered from the ledger and
    /// re-attested against the run identity before the engine opens.
    cards: &'a CardPackSet,
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
    cards: &CardPackSet,
    progress_sink: &mut Vec<String>,
) -> Result<SolveOutcome, SolveRefusal> {
    run_solve_inner(ledger, gate, clock, project, cards, progress_sink, None)
}

/// Deterministic conformance-test entry point for requesting cancellation at
/// one typed evidence-work checkpoint.
///
/// The plan can only request `gate`; it cannot execute caller code or publish
/// ledger state.
#[doc(hidden)]
pub fn run_solve_with_cancellation_plan(
    ledger: &Ledger,
    gate: &CancelGate,
    clock: &mut dyn FnMut() -> f64,
    project: &DecodedProject,
    cards: &CardPackSet,
    progress_sink: &mut Vec<String>,
    plan: &SolveCancellationPlan,
) -> Result<SolveOutcome, SolveRefusal> {
    run_solve_inner(
        ledger,
        gate,
        clock,
        project,
        cards,
        progress_sink,
        Some(plan),
    )
}

fn run_solve_inner<'a>(
    ledger: &'a Ledger,
    gate: &'a CancelGate,
    clock: &'a mut dyn FnMut() -> f64,
    project: &'a DecodedProject,
    cards: &'a CardPackSet,
    progress_sink: &mut Vec<String>,
    plan: Option<&'a SolveCancellationPlan>,
) -> Result<SolveOutcome, SolveRefusal> {
    let invocation = InvocationWorkLedger::default();
    let work = EvidenceWork::new(gate, plan, &invocation);
    work.checkpoint(SolveEvidencePhase::ProjectValidation, None, 0)
        .map_err(|_| cancelled_before_run_refusal())?;
    let findings = project.findings();
    work.charge(u64::try_from(project.canonical.len()).map_err(|_| {
        invocation_work_refusal(
            None,
            None,
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        )
    })?)
    .map_err(|error| invocation_work_refusal(None, None, error))?;
    work.checkpoint(SolveEvidencePhase::ProjectValidation, None, 1)
        .map_err(|_| cancelled_before_run_refusal())?;
    if !findings.is_empty() {
        return Err(SolveRefusal::plain(
            "cli-solve-project-invalid",
            format!("the project has {} validation findings", findings.len()),
            "run `frankensim validate` and repair every finding before solve",
        ));
    }
    work.checkpoint(SolveEvidencePhase::ProjectIdentityDerive, None, 0)
        .map_err(|_| cancelled_before_run_refusal())?;
    let project_hash = project.hash();
    let run = SolveRunId::derive_with_project_hash(project, project_hash, cards.root());
    work.charge(u64::try_from(project.canonical.len()).map_err(|_| {
        invocation_work_refusal(
            Some(run),
            None,
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        )
    })?)
    .map_err(|error| invocation_work_refusal(Some(run), None, error))?;
    work.checkpoint(SolveEvidencePhase::ProjectIdentityDerive, None, 1)
        .map_err(|_| cancelled_fresh_refusal(run, None))?;
    let mut engine =
        SolveEngine::open(ledger, work, clock, project, cards, project_hash, run, None)?;
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
    resume_solve_inner(ledger, gate, clock, run_id_hex, progress_sink, None)
}

/// Deterministic conformance-test entry point for requesting cancellation
/// while resume re-attests retained evidence.
#[doc(hidden)]
pub fn resume_solve_with_cancellation_plan(
    ledger: &Ledger,
    gate: &CancelGate,
    clock: &mut dyn FnMut() -> f64,
    run_id_hex: &str,
    progress_sink: &mut Vec<String>,
    plan: &SolveCancellationPlan,
) -> Result<SolveOutcome, SolveRefusal> {
    resume_solve_inner(ledger, gate, clock, run_id_hex, progress_sink, Some(plan))
}

fn resume_solve_inner<'a>(
    ledger: &'a Ledger,
    gate: &'a CancelGate,
    clock: &'a mut dyn FnMut() -> f64,
    run_id_hex: &str,
    progress_sink: &mut Vec<String>,
    plan: Option<&'a SolveCancellationPlan>,
) -> Result<SolveOutcome, SolveRefusal> {
    let run = SolveRunId::parse_hex(run_id_hex).ok_or_else(|| {
        SolveRefusal::plain(
            "cli-solve-run-id",
            format!("`{run_id_hex}` is not a 64-hex run id"),
            "pass the run id printed by `frankensim solve`",
        )
    })?;
    let invocation = InvocationWorkLedger::default();
    let work = EvidenceWork::new(gate, plan, &invocation);
    if work.is_requested() {
        return Err(cancelled_resume_refusal(run));
    }
    let verified = load_latest_state(ledger, run, work)?;
    let VerifiedResume {
        state,
        state_artifact,
        attestation,
        context,
        cards,
    } = verified;
    let attestation = Rc::try_unwrap(attestation).unwrap_or_else(|_| {
        unreachable!("latest-state selection releases every non-selected attestation reference")
    });
    let ResumeImportCache {
        project,
        project_hash,
        versions,
        budget,
        capability,
        seed,
        ..
    } = attestation;
    let rendered_explicits = (versions, budget, capability, seed);
    if state.completed.len() >= SolveStage::ALL.len() {
        return Err(SolveRefusal::plain(
            "cli-solve-resume-complete",
            format!("run `{}` already completed every stage", run.to_hex()),
            "use `frankensim report <run-id>` once reporting ships (f85xj.6.9)",
        ));
    }
    let mut engine = SolveEngine::open(
        ledger,
        work,
        clock,
        &project,
        &cards,
        project_hash,
        run,
        Some(rendered_explicits),
    )?;
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
        work: EvidenceWork<'a>,
        clock: &'a mut dyn FnMut() -> f64,
        project: &'a DecodedProject,
        cards: &'a CardPackSet,
        project_hash: ContentHash,
        run: SolveRunId,
        rendered_explicits: Option<(String, String, String, [u8; 8])>,
    ) -> Result<SolveEngine<'a>, SolveRefusal> {
        if ledger.in_transaction() {
            return Err(SolveRefusal::plain(
                "cli-solve-ledger-transaction",
                "the ledger connection is already inside a caller-owned transaction",
                "commit or roll back before solve so stage groups stay atomic",
            ));
        }
        let (versions_json, budget_json, capability_json, seed) =
            if let Some(rendered_explicits) = rendered_explicits {
                rendered_explicits
            } else {
                work.checkpoint(SolveEvidencePhase::FiveExplicitsRender, None, 0)
                    .map_err(|_| cancelled_fresh_refusal(run, None))?;
                let rendered_explicits = explicits(&project.spec);
                work.checkpoint(SolveEvidencePhase::FiveExplicitsRender, None, 1)
                    .map_err(|_| cancelled_fresh_refusal(run, None))?;
                let rendered_explicits = rendered_explicits.map_err(|refusal| SolveRefusal {
                    code: "cli-solve-project-invalid",
                    stage: None,
                    what: refusal.what,
                    fix: refusal.fix,
                    dependency: None,
                    run: Some(run.to_hex()),
                    recorded_op: None,
                })?;
                let explicit_bytes = rendered_explicits
                    .0
                    .len()
                    .checked_add(rendered_explicits.1.len())
                    .and_then(|bytes| bytes.checked_add(rendered_explicits.2.len()))
                    .and_then(|bytes| bytes.checked_add(rendered_explicits.3.len()))
                    .and_then(|bytes| u64::try_from(bytes).ok())
                    .ok_or_else(|| {
                        invocation_work_refusal(
                            Some(run),
                            None,
                            InvocationWorkExceeded::CumulativeBytes {
                                attempted: u64::MAX,
                            },
                        )
                    })?;
                work.charge(explicit_bytes)
                    .map_err(|error| invocation_work_refusal(Some(run), None, error))?;
                rendered_explicits
            };
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
            work,
            clock,
            spec: &project.spec,
            canonical_source: &project.canonical,
            cards,
            project_hash,
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
            if self.work.is_requested() {
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
            // Stages that retain no side evidence yield an empty usage list;
            // only `material-resolve` retains replayable claim-usage receipts.
            let body = match stage {
                SolveStage::ImportVerify => self
                    .stage_import_verify(&mut context)
                    .map(|receipt| (receipt, Vec::new())),
                SolveStage::Assign => self
                    .stage_assign(&context)
                    .map(|receipt| (receipt, Vec::new())),
                SolveStage::MaterialResolve => self.stage_material_resolve(),
                SolveStage::FlowNetwork => self
                    .stage_flow_network()
                    .map(|receipt| (receipt, Vec::new())),
                _ => unreachable!("gap stages returned above"),
            };
            let (receipt_json, usages) = match body {
                Ok(produced) => produced,
                Err(refusal) if refusal.code == "cli-solve-cancelled" => {
                    return Err(self.cancelled_refusal(&state));
                }
                Err(refusal) if refusal.code == "cli-solve-work-envelope" => return Err(refusal),
                Err(refusal) => return Err(self.record_refusal(&state, stage, refusal)),
            };
            if self
                .work
                .checkpoint(
                    SolveEvidencePhase::PrePublication,
                    None,
                    u64::from(stage.ordinal()),
                )
                .is_err()
            {
                return Err(self.cancelled_refusal(&state));
            }
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
                .persist_stage(
                    &state,
                    stage,
                    &receipt_json,
                    &usages,
                    &context,
                    predecessor_state,
                )
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
        if state.completed.is_empty() {
            return cancelled_fresh_refusal(self.run, SolveStage::from_ordinal(0));
        }
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
        let canonical_project_json =
            render_canonical_project_json(self.spec, self.work).map_err(|error| match error {
                ProjectRenderError::Cancelled => cancelled_fresh_refusal(self.run, Some(stage)),
                ProjectRenderError::WorkEnvelope(error) => {
                    invocation_work_refusal(Some(self.run), Some(stage), error)
                }
                ProjectRenderError::Invalid(problem) => SolveRefusal::staged(
                    "cli-solve-project-invalid",
                    stage,
                    problem,
                    "run `frankensim validate` and repair the canonical project renderer",
                ),
            })?;
        let expectations = RenderedExplicitsRef {
            seed: &self.seed,
            versions: &self.versions_json,
            budget: &self.budget_json,
            capability: &self.capability_json,
            canonical_project_json: &canonical_project_json,
        };
        let summary = match find_import_summary(
            self.ledger,
            self.spec,
            self.project_hash,
            self.work,
            expectations,
        ) {
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
            Err(ImportSummaryError::WorkEnvelope(error)) => {
                return Err(invocation_work_refusal(Some(self.run), Some(stage), error));
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
            Err(ImportSummaryError::Cancelled) => {
                return Err(cancelled_fresh_refusal(self.run, Some(stage)));
            }
        };
        let ImportSummary {
            op_id: import_op,
            artifact: import_summary,
            entries,
        } = summary;
        // Summary parsing has already checked every entry, in order, against
        // the same project geometry role and controlled source identity.
        // Candidate validation above streamed or materialized every retained
        // artifact exactly once under the solve evidence envelope.
        context.import_op = Some(import_op);
        context.import_summary = Some(import_summary);
        context.verified_imports = entries;
        import_verify_receipt(
            self.run,
            self.project_hash,
            import_op,
            &context.verified_imports,
            self.work,
        )
        .map_err(|error| match error {
            EvidenceCompareError::Cancelled => cancelled_fresh_refusal(self.run, Some(stage)),
            EvidenceCompareError::WorkEnvelope(error) => {
                invocation_work_refusal(Some(self.run), Some(stage), error)
            }
        })
    }

    /// Bind verified assignment evidence to the run's declared targets.
    fn stage_assign(&mut self, context: &StageContext) -> Result<String, SolveRefusal> {
        assignment_receipt(self.spec, context, self.run, self.work, false)
    }

    /// Resolve every declared material and interface binding against the
    /// admitted card packs.
    fn stage_material_resolve(&mut self) -> Result<(String, Vec<RetainedUsage>), SolveRefusal> {
        material_resolve_receipt(self.spec, self.cards, self.run, self.work, false)
    }

    fn stage_flow_network(&mut self) -> Result<String, SolveRefusal> {
        flow_network_receipt(self.spec, self.run, self.work, false)
    }

    /// Persist one completed stage as a ledgered op: stage receipt, sealed
    /// driver state (including this stage), and lineage links, atomically.
    fn persist_stage(
        &mut self,
        state_before: &SolveDriverState,
        stage: SolveStage,
        receipt_json: &str,
        usages: &[RetainedUsage],
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
                // Card packs are run inputs, so they are retained against the
                // run's first operation exactly as the project source is.
                // Retaining them here rather than at `material-resolve` is
                // what lets a run cancelled during an earlier stage resume at
                // all: the recovered set is re-attested against the run
                // identity, which binds its canonical root.
                for pack in self.cards.iter() {
                    let retained = self.ledger.put_artifact(
                        pack.kind().artifact_kind(),
                        pack.bytes(),
                        None,
                    )?;
                    self.ledger.link(op, &retained.hash, EdgeRole::In)?;
                }
            }
            if stage == SolveStage::MaterialResolve {
                for pack in self.cards.iter() {
                    self.ledger.link(op, &pack.artifact(), EdgeRole::In)?;
                }
                for usage in usages {
                    let retained =
                        self.ledger
                            .put_artifact(MATERIAL_USAGE_KIND, &usage.bytes, None)?;
                    self.ledger.link(op, &retained.hash, EdgeRole::Out)?;
                }
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
    work: EvidenceWork<'_>,
) -> Result<String, EvidenceCompareError> {
    let phase = SolveEvidencePhase::ReceiptDerivation;
    work.checkpoint(phase, None, 0)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    let mut receipt = format!(
        "{{\"schema\":{},\"run\":{},\"project_hash\":{},\"import_op\":{import_op},\"verified\":[",
        json_string(IMPORT_VERIFY_RECEIPT_SCHEMA),
        json_string(&run.to_hex()),
        json_string(&project_hash.to_hex()),
    );
    work.charge(u64::try_from(receipt.len()).map_err(|_| {
        EvidenceCompareError::WorkEnvelope(InvocationWorkExceeded::CumulativeBytes {
            attempted: u64::MAX,
        })
    })?)
    .map_err(EvidenceCompareError::WorkEnvelope)?;
    for (index, entry) in entries.iter().enumerate() {
        let units = u64::try_from(index).unwrap_or(u64::MAX);
        work.checkpoint(phase, None, units)
            .map_err(|_| EvidenceCompareError::Cancelled)?;
        let before = receipt.len();
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
        let appended =
            receipt
                .len()
                .checked_sub(before)
                .ok_or(EvidenceCompareError::WorkEnvelope(
                    InvocationWorkExceeded::CumulativeBytes {
                        attempted: u64::MAX,
                    },
                ))?;
        let charge = u64::try_from(appended)
            .ok()
            .and_then(|bytes| bytes.checked_add(DERIVATION_ITEM_WORK_BYTES))
            .ok_or(EvidenceCompareError::WorkEnvelope(
                InvocationWorkExceeded::CumulativeBytes {
                    attempted: u64::MAX,
                },
            ))?;
        work.charge(charge)
            .map_err(EvidenceCompareError::WorkEnvelope)?;
        work.checkpoint(phase, None, units.saturating_add(1))
            .map_err(|_| EvidenceCompareError::Cancelled)?;
    }
    let before = receipt.len();
    let _ = write!(
        receipt,
        "],\"authority\":{},\"no_claim\":{}}}",
        json_string(IMPORT_VERIFY_AUTHORITY),
        json_string(IMPORT_VERIFY_NO_CLAIM),
    );
    let suffix = receipt
        .len()
        .checked_sub(before)
        .ok_or(EvidenceCompareError::WorkEnvelope(
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        ))?;
    work.charge(u64::try_from(suffix).map_err(|_| {
        EvidenceCompareError::WorkEnvelope(InvocationWorkExceeded::CumulativeBytes {
            attempted: u64::MAX,
        })
    })?)
    .map_err(EvidenceCompareError::WorkEnvelope)?;
    work.checkpoint(phase, None, u64::MAX)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    Ok(receipt)
}

/// One replayable claim-usage receipt retained as its own ledger artifact.
///
/// `id` is the `fs-matdb` receipt identity (what `ClaimSet::verify_receipt`
/// replays against); `artifact` is the ledger content address of the exact
/// retained bytes. They are different hashes over different preimages and the
/// stage receipt records both, so neither can stand in for the other.
#[derive(Debug, Clone, PartialEq, Eq)]
struct RetainedUsage {
    id: String,
    artifact: ContentHash,
    bytes: Vec<u8>,
}

/// Canonical finite `f64` spelling, or `None` for a non-finite value.
///
/// The retained-receipt reader requires `value.to_string() == spelling`, so a
/// non-finite value has no admissible canonical form and must refuse rather
/// than reach a receipt.
fn canonical_f64(value: f64) -> Option<String> {
    value.is_finite().then(|| value.to_string())
}

fn material_nonfinite(stage: SolveStage, what: impl Into<String>) -> SolveRefusal {
    SolveRefusal::staged(
        "cli-solve-material-nonfinite",
        stage,
        what,
        "repair the card claim so every resolved endpoint is a finite quantity",
    )
}

fn dims_json(dims: [i8; 6]) -> String {
    let mut out = String::from("[");
    for (index, exponent) in dims.into_iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "{exponent}");
    }
    out.push(']');
    out
}

fn required_properties_json(properties: &[fs_project::RequiredProperty]) -> String {
    let mut out = String::from("[");
    for (index, required) in properties.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"property\":{},\"dims\":{}}}",
            json_string(&required.property),
            dims_json(required.dims.0)
        );
    }
    out.push(']');
    out
}

fn uncertainty_json(
    stage: SolveStage,
    context: &str,
    uncertainty: &fs_matdb::UncertaintyModel,
) -> Result<String, SolveRefusal> {
    Ok(match *uncertainty {
        fs_matdb::UncertaintyModel::Unstated => "{\"kind\":\"unstated\"}".to_string(),
        fs_matdb::UncertaintyModel::HalfWidth {
            half_width,
            confidence,
        } => {
            let (half_width, confidence) = finite_pair(stage, context, half_width, confidence)?;
            format!(
                "{{\"kind\":\"half-width\",\"half_width\":{half_width},\"confidence\":{confidence}}}"
            )
        }
        fs_matdb::UncertaintyModel::RelativeHalfWidth {
            fraction,
            confidence,
        } => {
            let (fraction, confidence) = finite_pair(stage, context, fraction, confidence)?;
            format!(
                "{{\"kind\":\"relative-half-width\",\"fraction\":{fraction},\
                 \"confidence\":{confidence}}}"
            )
        }
    })
}

fn finite_pair(
    stage: SolveStage,
    context: &str,
    first: f64,
    second: f64,
) -> Result<(String, String), SolveRefusal> {
    let first = canonical_f64(first).ok_or_else(|| {
        material_nonfinite(stage, format!("{context} carries a non-finite value"))
    })?;
    let second = canonical_f64(second).ok_or_else(|| {
        material_nonfinite(stage, format!("{context} carries a non-finite value"))
    })?;
    Ok((first, second))
}

/// Resolve the project's declared bindings against the admitted card packs
/// and render the canonical `material-resolve` stage receipt.
///
/// Deliberately a free function so resume can rebuild the byte-identical
/// receipt from re-attested evidence, exactly as [`assignment_receipt`] does.
///
/// The receipt names packs and cards by content root only. Caller paths never
/// appear, so the same content admitted from a different path reproduces the
/// same receipt — which is what lets resume, whose only inputs are the ledger
/// and the run id, reproduce it at all.
fn material_resolve_receipt(
    spec: &ProjectSpec,
    cards: &CardPackSet,
    run: SolveRunId,
    work: EvidenceWork<'_>,
    resume: bool,
) -> Result<(String, Vec<RetainedUsage>), SolveRefusal> {
    let stage = SolveStage::MaterialResolve;
    let cancelled = || {
        if resume {
            cancelled_resume_refusal(run)
        } else {
            cancelled_fresh_refusal(run, Some(stage))
        }
    };
    let phase = SolveEvidencePhase::MaterialBindingResolution;
    let mut units = 0u64;
    work.checkpoint(phase, None, units)
        .map_err(|_| cancelled())?;

    // Building the library never chooses a claim: every admitted card enters
    // with its complete claim set, and selection happens inside
    // `resolve_bindings`, which leaves a replayable receipt for what it used.
    let library = cards.library();
    for _ in 0..cards.len() {
        work.charge(DERIVATION_ITEM_WORK_BYTES)
            .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
        units = units.saturating_add(1);
        work.checkpoint(phase, None, units)
            .map_err(|_| cancelled())?;
    }

    let requirements = BindingRequirements::thermal_steady_v1();
    let resolution = resolve_bindings(spec, &library, &requirements);
    work.charge(DERIVATION_ITEM_WORK_BYTES)
        .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
    units = units.saturating_add(1);
    work.checkpoint(phase, None, units)
        .map_err(|_| cancelled())?;

    if !resolution.admissible() {
        // The binding layer owns this vocabulary. Propagating its own stable
        // code keeps one diagnostic identity across the library boundary
        // instead of flattening every cause into one CLI code.
        let first = &resolution.violations[0];
        let mut what = format!(
            "{} declared material/interface binding violation(s) against the admitted card \
             packs; first: {}",
            resolution.violations.len(),
            first.what
        );
        for violation in resolution.violations.iter().skip(1).take(7) {
            let _ = write!(what, "; also [{}] {}", violation.code, violation.what);
        }
        if resolution.violations.len() > 8 {
            let _ = write!(
                what,
                "; and {} further violation(s)",
                resolution.violations.len() - 8
            );
        }
        return Err(SolveRefusal::staged(
            first.code,
            stage,
            what,
            first.fix.clone(),
        ));
    }

    let mut usages: Vec<RetainedUsage> = Vec::new();
    let mut bindings = String::new();
    for (index, binding) in resolution.bindings.iter().enumerate() {
        if index > 0 {
            bindings.push(',');
        }
        let (target_kind, target_name) = match &binding.target {
            fs_project::BindingTarget::Region(name) => ("region", name),
            fs_project::BindingTarget::Interface(name) => ("interface", name),
        };
        let context = format!("{target_kind} `{target_name}`");
        let (range_lo, range_hi) = finite_pair(
            stage,
            &format!("{context} admitted range"),
            binding.range_lo,
            binding.range_hi,
        )?;
        let mut properties = String::new();
        for (property_index, property) in binding.properties.iter().enumerate() {
            if property_index > 0 {
                properties.push(',');
            }
            let property_context = format!("{context} property `{}`", property.property);
            let (value_lo, value_hi) = finite_pair(
                stage,
                &property_context,
                property.value_lo,
                property.value_hi,
            )?;
            let uncertainty = uncertainty_json(stage, &property_context, &property.uncertainty)?;
            for retained in [&property.receipt_lo, &property.receipt_hi] {
                let artifact = hash_bytes(&retained.bytes);
                if !usages.iter().any(|usage| usage.artifact == artifact) {
                    usages.push(RetainedUsage {
                        id: retained.receipt_hash.clone(),
                        artifact,
                        bytes: retained.bytes.clone(),
                    });
                }
            }
            let _ = write!(
                properties,
                "{{\"property\":{},\"value_lo\":{value_lo},\"value_hi\":{value_hi},\
                 \"dims\":{},\"uncertainty\":{uncertainty},\"unstated_uncertainty\":{},\
                 \"selected_claim\":{},\"regime_card\":{{\"name\":{},\"version\":{}}},\
                 \"provenance_source\":{},\"provenance_license\":{},\
                 \"usage_receipt_lo\":{{\"id\":{},\"artifact\":{}}},\
                 \"usage_receipt_hi\":{{\"id\":{},\"artifact\":{}}}}}",
                json_string(&property.property),
                dims_json(property.dims.0),
                property.unstated_uncertainty,
                json_string(&property.selected_claim),
                json_string(&property.regime_card.name),
                json_string(&property.regime_card.version),
                json_string(&property.provenance_source),
                json_string(&property.provenance_license),
                json_string(&property.receipt_lo.receipt_hash),
                json_string(&hash_bytes(&property.receipt_lo.bytes).to_hex()),
                json_string(&property.receipt_hi.receipt_hash),
                json_string(&hash_bytes(&property.receipt_hi.bytes).to_hex()),
            );
            work.charge(DERIVATION_ITEM_WORK_BYTES)
                .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
            units = units.saturating_add(1);
            work.checkpoint(phase, None, units)
                .map_err(|_| cancelled())?;
        }
        let _ = write!(
            bindings,
            "{{\"target_kind\":{},\"target\":{},\"card\":{},\"card_identity\":{},\
             \"declared_source\":{},\"declared_interface_state\":{},\
             \"range_lo\":{range_lo},\"range_hi\":{range_hi},\"pinned_claim\":{},\
             \"properties\":[{properties}]}}",
            json_string(target_kind),
            json_string(target_name),
            json_string(&binding.card),
            json_string(&binding.card_identity),
            json_string(&binding.declared_source),
            binding
                .declared_interface_state
                .as_ref()
                .map_or_else(|| "null".to_string(), |state| json_string(state)),
            binding
                .pinned_claim
                .as_ref()
                .map_or_else(|| "null".to_string(), |claim| json_string(claim)),
        );
        work.charge(DERIVATION_ITEM_WORK_BYTES)
            .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
        units = units.saturating_add(1);
        work.checkpoint(phase, None, units)
            .map_err(|_| cancelled())?;
    }

    if usages.len() > MAX_MATERIAL_USAGE_RECEIPTS {
        return Err(SolveRefusal::staged(
            "cli-solve-material-receipt-envelope",
            stage,
            format!(
                "resolving this project retains {} distinct claim-usage receipts, above the \
                 {MAX_MATERIAL_USAGE_RECEIPTS}-receipt stage envelope",
                usages.len()
            ),
            "reduce the number of declared bindings, or raise the documented stage envelope \
             deliberately",
        ));
    }

    let mut packs = String::new();
    for (index, pack) in cards.iter().enumerate() {
        if index > 0 {
            packs.push(',');
        }
        let _ = write!(
            packs,
            "{{\"kind\":{},\"root\":{},\"card\":{},\"identity\":{}}}",
            json_string(pack.kind().label()),
            json_string(&pack.root().to_hex()),
            json_string(&pack.card().to_hex()),
            json_string(pack.identity()),
        );
    }

    let mut advisories = String::new();
    for (index, advisory) in resolution.advisories.iter().enumerate() {
        if index > 0 {
            advisories.push(',');
        }
        let _ = write!(
            advisories,
            "{{\"code\":{},\"what\":{},\"note\":{}}}",
            json_string(advisory.code),
            json_string(&advisory.what),
            json_string(&advisory.note),
        );
    }

    let receipt = format!(
        "{{\"schema\":{},\"run\":{},\"pack_set_root\":{},\"packs\":[{packs}],\
         \"requirements\":{{\"temperature_axis\":{},\"material_properties\":{},\
         \"interface_properties\":{}}},\"bindings\":[{bindings}],\
         \"advisories\":[{advisories}],\"authority\":{},\"no_claim\":{}}}",
        json_string(MATERIAL_RESOLVE_RECEIPT_SCHEMA),
        json_string(&run.to_hex()),
        json_string(&cards.root().to_hex()),
        json_string(&requirements.temperature_axis),
        required_properties_json(&requirements.material_properties),
        required_properties_json(&requirements.interface_properties),
        json_string(MATERIAL_RESOLVE_AUTHORITY),
        json_string(MATERIAL_RESOLVE_NO_CLAIM),
    );
    work.charge(u64::try_from(receipt.len()).map_err(|_| {
        invocation_work_refusal(
            Some(run),
            Some(stage),
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        )
    })?)
    .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
    work.checkpoint(phase, None, u64::MAX)
        .map_err(|_| cancelled())?;
    Ok((receipt, usages))
}

fn flow_network_receipt(
    spec: &ProjectSpec,
    run: SolveRunId,
    work: EvidenceWork<'_>,
    resume: bool,
) -> Result<String, SolveRefusal> {
    let stage = SolveStage::FlowNetwork;
    let cancelled = || {
        if resume {
            cancelled_resume_refusal(run)
        } else {
            cancelled_fresh_refusal(run, Some(stage))
        }
    };
    let phase = SolveEvidencePhase::FlowNetworkSolve;
    let mut units = 0u64;
    work.checkpoint(phase, None, units)
        .map_err(|_| cancelled())?;

    let cooling = spec.cooling.as_ref().ok_or_else(|| {
        SolveRefusal::staged(
            "cli-solve-flow-network-no-cooling",
            stage,
            "the project declares no cooling section".to_string(),
            "declare the cooling section; the flow-network stage consumes declared fans, vents, and leakage",
        )
    })?;
    let fan_system = cooling.fan_system.as_ref().ok_or_else(|| {
        SolveRefusal::staged(
            "cli-solve-flow-network-no-fan-system",
            stage,
            "the cooling section declares no fan system".to_string(),
            "declare `(fan-system ...)` under cooling (schema v2); the stage never infers banks, speeds, or topology",
        )
    })?;
    let leakage = cooling.airflow_leakage.as_ref().ok_or_else(|| {
        SolveRefusal::staged(
            "cli-solve-flow-network-no-leakage",
            stage,
            "the cooling section declares no airflow leakage area".to_string(),
            "declare `(airflow-leakage :area ...)`; the stage refuses to invent the mandatory leakage branch",
        )
    })?;
    let envelope = spec.envelope.as_ref().ok_or_else(|| {
        SolveRefusal::staged(
            "cli-solve-flow-network-no-envelope",
            stage,
            "the project declares no operating envelope".to_string(),
            "declare the envelope; the air density estimate needs ambient pressure and temperature",
        )
    })?;
    units = units.saturating_add(1);
    work.checkpoint(phase, None, units)
        .map_err(|_| cancelled())?;

    // Lossless lowering of the versioned fan-system declaration (bead
    // frn2i.1): every bank is consumed exactly once and the composite
    // carries member-bound provenance.
    let lowered = fs_project::fansystem::lower_fan_system(fan_system).map_err(|error| {
        SolveRefusal::staged(
            "cli-solve-flow-network-lowering",
            stage,
            format!("fan-system lowering refused: {}", error.detail),
            error.hint,
        )
    })?;
    work.charge(DERIVATION_ITEM_WORK_BYTES)
        .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
    units = units.saturating_add(1);
    work.checkpoint(phase, None, units)
        .map_err(|_| cancelled())?;

    // Envelope-derived ideal-gas density at the ambient midpoint: a typed,
    // receipted estimate, never a measurement.
    let ambient_mid = (envelope.ambient_lo.value + envelope.ambient_hi.value) / 2.0;
    let air_density = envelope.pressure.value / (AIR_SPECIFIC_GAS_CONSTANT * ambient_mid);
    if !(air_density.is_finite() && air_density > 0.0) {
        return Err(SolveRefusal::staged(
            "cli-solve-flow-network-density",
            stage,
            format!("the envelope yields a non-positive or non-finite air density ({air_density})"),
            "check ambient temperature and pressure bounds",
        ));
    }
    let density = fs_qty::Density::new(air_density);

    // Vents lower to sharp-edged-orifice loss elements in parallel; the
    // declared leakage area is the mandatory leakage branch.
    let mut branches: Vec<fs_airflow::LossNetwork> = Vec::new();
    for vent in &cooling.vents {
        if vent.area.dims != fs_project::spec::dims::AREA {
            return Err(SolveRefusal::staged(
                "cli-solve-flow-network-vent-units",
                stage,
                format!(
                    "vent `{}` area carries dims {}",
                    vent.region,
                    vent.area.dims.unit_string()
                ),
                "vent areas must carry m^2 dimensions",
            ));
        }
        let element = fs_airflow::sharp_edged_orifice_loss(
            format!("vent:{}", vent.region),
            fs_qty::Area::new(vent.area.value),
            density,
        )
        .map_err(|error| {
            SolveRefusal::staged(
                "cli-solve-flow-network-vent",
                stage,
                format!("vent `{}` refused: {error}", vent.region),
                "declare a positive finite vent area",
            )
        })?;
        work.charge(DERIVATION_ITEM_WORK_BYTES)
            .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
        units = units.saturating_add(1);
        work.checkpoint(phase, None, units)
            .map_err(|_| cancelled())?;
        branches.push(fs_airflow::LossNetwork::Element(element));
    }
    if branches.is_empty() {
        return Err(SolveRefusal::staged(
            "cli-solve-flow-network-no-vents",
            stage,
            "the cooling section declares no vents".to_string(),
            "declare at least one vent area; the stage never invents flow paths",
        ));
    }
    let primary = if branches.len() == 1 {
        branches.pop().expect("one branch")
    } else {
        fs_airflow::LossNetwork::parallel(branches).map_err(|error| {
            SolveRefusal::staged(
                "cli-solve-flow-network-network",
                stage,
                format!("vent network composition refused: {error}"),
                "report this; validated vents cannot fail parallel composition",
            )
        })?
    };
    let leakage_element = if leakage.area.dims != fs_project::spec::dims::AREA {
        return Err(SolveRefusal::staged(
            "cli-solve-flow-network-leakage-units",
            stage,
            format!("leakage area carries dims {}", leakage.area.dims.unit_string()),
            "the leakage area must carry m^2 dimensions",
        ));
    } else {
        fs_airflow::sharp_edged_orifice_loss(
            "leakage",
            fs_qty::Area::new(leakage.area.value),
            density,
        )
        .map_err(|error| {
            SolveRefusal::staged(
                "cli-solve-flow-network-leakage",
                stage,
                format!("leakage lowering refused: {error}"),
                "declare a positive finite leakage area",
            )
        })?
    };
    let network = fs_airflow::EnclosureNetwork::new(
        primary,
        fs_airflow::LeakageElement::new(leakage_element),
    );
    units = units.saturating_add(1);
    work.checkpoint(phase, None, units)
        .map_err(|_| cancelled())?;

    let operating =
        fs_airflow::solve_operating_point(&lowered.system_bank, &network).map_err(|error| {
            SolveRefusal::staged(
                "cli-solve-flow-network-solve",
                stage,
                format!("the interval-certified operating-point solve refused: {error}"),
                "check fan domains against the network's resistance, or revisit the declaration",
            )
        })?;
    work.charge(DERIVATION_ITEM_WORK_BYTES)
        .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
    units = units.saturating_add(1);
    work.checkpoint(phase, None, units)
        .map_err(|_| cancelled())?;

    let flow = operating.nominal_root.flow;
    let pressure_lo = operating.pressure.numerical.lo;
    let pressure_hi = operating.pressure.numerical.hi;
    let receipt = format!(
        "{{\"schema\":{},\"run\":{},\"stage\":{},\"declaration\":{},\
         \"composite\":{},\"density_estimate\":{},\"vent_count\":{},\
         \"operating_point\":{{\"flow_lo\":{},\"flow_hi\":{},\"flow_mid\":{},\
         \"pressure_lo\":{},\"pressure_hi\":{}}},\"leakage_fraction\":{},\
         \"authority\":{},\"no_claim\":{}}}",
        json_string(FLOW_NETWORK_RECEIPT_SCHEMA),
        json_string(&run.to_hex()),
        json_string(stage.name()),
        json_string(&lowered.declaration_identity),
        json_string(&lowered.system_bank.curve().source().identifier),
        json_string(&format!("ideal-gas:{air_density}")),
        cooling.vents.len(),
        json_string(&flow.lo().to_string()),
        json_string(&flow.hi().to_string()),
        json_string(&flow.midpoint().to_string()),
        json_string(&pressure_lo.to_string()),
        json_string(&pressure_hi.to_string()),
        json_string(&operating.leakage_fraction.to_string()),
        json_string(FLOW_NETWORK_AUTHORITY),
        json_string(FLOW_NETWORK_NO_CLAIM),
    );
    work.charge(u64::try_from(receipt.len()).map_err(|_| {
        invocation_work_refusal(
            Some(run),
            Some(stage),
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        )
    })?)
    .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
    work.checkpoint(phase, None, u64::MAX)
        .map_err(|_| cancelled())?;
    Ok(receipt)
}

fn assignment_receipt(
    spec: &ProjectSpec,
    context: &StageContext,
    run: SolveRunId,
    work: EvidenceWork<'_>,
    resume: bool,
) -> Result<String, SolveRefusal> {
    let stage = SolveStage::Assign;
    let cancelled = || {
        if resume {
            cancelled_resume_refusal(run)
        } else {
            cancelled_fresh_refusal(run, Some(stage))
        }
    };
    work.checkpoint(SolveEvidencePhase::AssignmentDerivation, None, 0)
        .map_err(|_| cancelled())?;
    let assignments = spec.assignments.as_deref().unwrap_or(&[]);
    if assignments.is_empty() {
        work.checkpoint(SolveEvidencePhase::AssignmentDerivation, None, 1)
            .map_err(|_| cancelled())?;
        return Err(SolveRefusal::staged(
            "cli-solve-assignment",
            stage,
            "the project declares no geometry assignments",
            "declare assignments before solve",
        ));
    }
    let mut rows = String::new();
    let mut assignment_units = 0u64;
    for (index, assignment) in assignments.iter().enumerate() {
        work.checkpoint(
            SolveEvidencePhase::AssignmentDerivation,
            None,
            assignment_units,
        )
        .map_err(|_| cancelled())?;
        let mut imported = None;
        for entry in &context.verified_imports {
            let matches = evidence_bytes_equal(
                entry.role.as_bytes(),
                assignment.artifact.as_bytes(),
                work,
                SolveEvidencePhase::AssignmentDerivation,
                None,
            )
            .map_err(|error| match error {
                EvidenceCompareError::Cancelled => cancelled(),
                EvidenceCompareError::WorkEnvelope(error) => {
                    invocation_work_refusal(Some(run), Some(stage), error)
                }
            })?;
            work.charge(DERIVATION_ITEM_WORK_BYTES)
                .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
            assignment_units = assignment_units.saturating_add(1);
            work.checkpoint(
                SolveEvidencePhase::AssignmentDerivation,
                None,
                assignment_units,
            )
            .map_err(|_| cancelled())?;
            if matches {
                imported = Some(entry);
                break;
            }
        }
        work.charge(DERIVATION_ITEM_WORK_BYTES)
            .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
        assignment_units = assignment_units.saturating_add(1);
        work.checkpoint(
            SolveEvidencePhase::AssignmentDerivation,
            None,
            assignment_units,
        )
        .map_err(|_| cancelled())?;
        let imported = imported.ok_or_else(|| {
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
    let receipt = format!(
        "{{\"schema\":{},\"run\":{},\"bindings\":[{rows}],\"authority\":\"declared targets \
         bound to verified import evidence\",\"no_claim\":\"selector re-resolution against \
         the mesh is the import's retained report; this stage does not re-run it\"}}",
        json_string(ASSIGN_RECEIPT_SCHEMA),
        json_string(&run.to_hex()),
    );
    work.checkpoint(SolveEvidencePhase::AssignmentDerivation, None, u64::MAX)
        .map_err(|_| cancelled())?;
    work.charge(u64::try_from(receipt.len()).map_err(|_| {
        invocation_work_refusal(
            Some(run),
            Some(stage),
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        )
    })?)
    .map_err(|error| invocation_work_refusal(Some(run), Some(stage), error))?;
    Ok(receipt)
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
    Cancelled,
    WorkEnvelope(InvocationWorkExceeded),
    Ledger(LedgerError),
    Invalid(String),
    Unsupported(String),
}

impl From<LedgerError> for ImportSummaryError {
    fn from(error: LedgerError) -> Self {
        Self::Ledger(error)
    }
}

impl From<EvidenceReadError> for ImportSummaryError {
    fn from(error: EvidenceReadError) -> Self {
        match error {
            EvidenceReadError::Cancelled => Self::Cancelled,
            EvidenceReadError::WorkEnvelope(error) => Self::WorkEnvelope(error),
            EvidenceReadError::Ledger(error) => Self::Ledger(error),
        }
    }
}

impl From<CandidateOpReadError> for ImportSummaryError {
    fn from(error: CandidateOpReadError) -> Self {
        match error {
            CandidateOpReadError::Cancelled => Self::Cancelled,
            CandidateOpReadError::WorkEnvelope(error) => Self::WorkEnvelope(error),
            CandidateOpReadError::Ledger(error) => Self::Ledger(error),
        }
    }
}

impl From<VisibleOpScanError> for ImportSummaryError {
    fn from(error: VisibleOpScanError) -> Self {
        match error {
            VisibleOpScanError::Cancelled => Self::Cancelled,
            VisibleOpScanError::WorkEnvelope(error) => Self::WorkEnvelope(error),
            VisibleOpScanError::Ledger(error) => Self::Ledger(error),
        }
    }
}

fn materialize_evidence_artifact(
    ledger: &Ledger,
    work: EvidenceWork<'_>,
    artifact: ContentHash,
    cap: u64,
    phase: SolveEvidencePhase,
    source_index: Option<usize>,
) -> Result<Option<Vec<u8>>, EvidenceReadError> {
    work.checkpoint(phase, source_index, 0)
        .map_err(|_| EvidenceReadError::Cancelled)?;
    let mut bytes = Vec::new();
    let mut processed = 0u64;
    let controlled = ledger.read_artifact_chunks_bounded_controlled(&artifact, cap, &mut |tile| {
        if work.checkpoint(phase, source_index, processed).is_err() {
            return ControlFlow::Break(EvidenceReadError::Cancelled);
        }
        let Ok(tile_len) = u64::try_from(tile.len()) else {
            return ControlFlow::Break(EvidenceReadError::Ledger(LedgerError::Invalid {
                field: "solve evidence materialization".to_string(),
                problem: "tile length is outside the ledger byte range".to_string(),
            }));
        };
        if let Err(error) = work.charge(tile_len) {
            return ControlFlow::Break(EvidenceReadError::WorkEnvelope(error));
        }
        if bytes.try_reserve(tile.len()).is_err() {
            return ControlFlow::Break(EvidenceReadError::Ledger(LedgerError::Invalid {
                field: "solve evidence materialization".to_string(),
                problem: "allocation refused under the admitted evidence envelope".to_string(),
            }));
        }
        bytes.extend_from_slice(tile);
        let Some(next) = processed.checked_add(tile_len) else {
            return ControlFlow::Break(EvidenceReadError::Ledger(LedgerError::Invalid {
                field: "solve evidence materialization".to_string(),
                problem: "materialized byte count overflowed u64".to_string(),
            }));
        };
        processed = next;
        if work.checkpoint(phase, source_index, processed).is_err() {
            return ControlFlow::Break(EvidenceReadError::Cancelled);
        }
        ControlFlow::Continue(())
    });
    work.checkpoint(phase, source_index, processed)
        .map_err(|_| EvidenceReadError::Cancelled)?;
    let controlled = controlled.map_err(EvidenceReadError::Ledger)?;
    match controlled {
        None => Ok(None),
        Some(ControlFlow::Break(error)) => Err(error),
        Some(ControlFlow::Continue(streamed)) if streamed == processed => Ok(Some(bytes)),
        Some(ControlFlow::Continue(streamed)) => {
            Err(EvidenceReadError::Ledger(LedgerError::Invalid {
                field: "solve evidence materialization".to_string(),
                problem: format!(
                    "controlled read completed {streamed} bytes but delivered {processed}"
                ),
            }))
        }
    }
}

fn verify_evidence_artifact(
    ledger: &Ledger,
    work: EvidenceWork<'_>,
    artifact: ContentHash,
    cap: u64,
    phase: SolveEvidencePhase,
    source_index: Option<usize>,
) -> Result<Option<u64>, EvidenceReadError> {
    work.checkpoint(phase, source_index, 0)
        .map_err(|_| EvidenceReadError::Cancelled)?;
    let mut processed = 0u64;
    let controlled = ledger.read_artifact_chunks_bounded_controlled(&artifact, cap, &mut |tile| {
        if work.checkpoint(phase, source_index, processed).is_err() {
            return ControlFlow::Break(EvidenceReadError::Cancelled);
        }
        let Ok(tile_len) = u64::try_from(tile.len()) else {
            return ControlFlow::Break(EvidenceReadError::Ledger(LedgerError::Invalid {
                field: "solve evidence verification".to_string(),
                problem: "tile length is outside the ledger byte range".to_string(),
            }));
        };
        if let Err(error) = work.charge(tile_len) {
            return ControlFlow::Break(EvidenceReadError::WorkEnvelope(error));
        }
        let Some(next) = processed.checked_add(tile_len) else {
            return ControlFlow::Break(EvidenceReadError::Ledger(LedgerError::Invalid {
                field: "solve evidence verification".to_string(),
                problem: "verified byte count overflowed u64".to_string(),
            }));
        };
        processed = next;
        if work.checkpoint(phase, source_index, processed).is_err() {
            return ControlFlow::Break(EvidenceReadError::Cancelled);
        }
        ControlFlow::Continue(())
    });
    work.checkpoint(phase, source_index, processed)
        .map_err(|_| EvidenceReadError::Cancelled)?;
    let controlled = controlled.map_err(EvidenceReadError::Ledger)?;
    match controlled {
        None => Ok(None),
        Some(ControlFlow::Break(error)) => Err(error),
        Some(ControlFlow::Continue(streamed)) if streamed == processed => Ok(Some(streamed)),
        Some(ControlFlow::Continue(streamed)) => {
            Err(EvidenceReadError::Ledger(LedgerError::Invalid {
                field: "solve evidence verification".to_string(),
                problem: format!(
                    "controlled read completed {streamed} bytes but delivered {processed}"
                ),
            }))
        }
    }
}

fn evidence_utf8_string(
    bytes: &[u8],
    work: EvidenceWork<'_>,
    phase: SolveEvidencePhase,
    source_index: Option<usize>,
    label: &str,
) -> Result<String, EvidenceUtf8Error> {
    work.checkpoint(phase, source_index, 0)
        .map_err(|_| EvidenceUtf8Error::Cancelled)?;
    let mut output = String::new();
    let mut cursor = 0usize;
    while cursor < bytes.len() {
        let proposed_end = cursor
            .checked_add(EVIDENCE_POLL_BYTES)
            .unwrap_or(bytes.len())
            .min(bytes.len());
        let mut end = proposed_end;
        if end < bytes.len() {
            let mut backed_up = 0usize;
            while end > cursor && bytes[end] & 0b1100_0000 == 0b1000_0000 {
                end -= 1;
                backed_up += 1;
                if backed_up > 3 {
                    end = proposed_end;
                    break;
                }
            }
        }
        if end == cursor {
            end = proposed_end;
        }
        let tile = &bytes[cursor..end];
        let text = core::str::from_utf8(tile);
        let inspected = u64::try_from(end)
            .map_err(|_| EvidenceUtf8Error::Invalid(format!("{label} length is outside u64")))?;
        work.checkpoint(phase, source_index, inspected)
            .map_err(|_| EvidenceUtf8Error::Cancelled)?;
        let text = text.map_err(|error| {
            EvidenceUtf8Error::Invalid(format!(
                "{label} is not UTF-8 at byte {}",
                cursor.saturating_add(error.valid_up_to())
            ))
        })?;
        let reserve = output.try_reserve(text.len());
        work.checkpoint(phase, source_index, inspected)
            .map_err(|_| EvidenceUtf8Error::Cancelled)?;
        reserve.map_err(|_| {
            EvidenceUtf8Error::Invalid(format!("{label} UTF-8 output allocation was refused"))
        })?;
        let tile_len = u64::try_from(tile.len()).map_err(|_| {
            EvidenceUtf8Error::WorkEnvelope(InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            })
        })?;
        work.charge(tile_len)
            .map_err(EvidenceUtf8Error::WorkEnvelope)?;
        output.push_str(text);
        cursor = end;
        let processed = u64::try_from(cursor)
            .map_err(|_| EvidenceUtf8Error::Invalid(format!("{label} length is outside u64")))?;
        work.checkpoint(phase, source_index, processed)
            .map_err(|_| EvidenceUtf8Error::Cancelled)?;
    }
    let completed = u64::try_from(cursor)
        .map_err(|_| EvidenceUtf8Error::Invalid(format!("{label} length is outside u64")))?;
    work.checkpoint(phase, source_index, completed)
        .map_err(|_| EvidenceUtf8Error::Cancelled)?;
    Ok(output)
}

fn evidence_bytes_equal(
    left: &[u8],
    right: &[u8],
    work: EvidenceWork<'_>,
    phase: SolveEvidencePhase,
    source_index: Option<usize>,
) -> Result<bool, EvidenceCompareError> {
    work.checkpoint(phase, source_index, 0)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    if left.len() != right.len() {
        work.checkpoint(phase, source_index, 0)
            .map_err(|_| EvidenceCompareError::Cancelled)?;
        return Ok(false);
    }
    let mut compared = 0u64;
    for (left_tile, right_tile) in left
        .chunks(EVIDENCE_POLL_BYTES)
        .zip(right.chunks(EVIDENCE_POLL_BYTES))
    {
        work.checkpoint(phase, source_index, compared)
            .map_err(|_| EvidenceCompareError::Cancelled)?;
        let matches = left_tile == right_tile;
        let tile_len = u64::try_from(left_tile.len()).map_err(|_| {
            EvidenceCompareError::WorkEnvelope(InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            })
        })?;
        compared = compared
            .checked_add(tile_len)
            .ok_or(EvidenceCompareError::WorkEnvelope(
                InvocationWorkExceeded::CumulativeBytes {
                    attempted: u64::MAX,
                },
            ))?;
        work.checkpoint(phase, source_index, compared)
            .map_err(|_| EvidenceCompareError::Cancelled)?;
        work.charge(tile_len)
            .map_err(EvidenceCompareError::WorkEnvelope)?;
        if !matches {
            return Ok(false);
        }
    }
    Ok(true)
}

fn render_canonical_project_json(
    spec: &ProjectSpec,
    work: EvidenceWork<'_>,
) -> Result<String, ProjectRenderError> {
    let phase = SolveEvidencePhase::CanonicalProjectRender;
    work.checkpoint(phase, None, 0)
        .map_err(|_| ProjectRenderError::Cancelled)?;
    let rendered = fs_project::print_json(spec);
    work.checkpoint(phase, None, 1)
        .map_err(|_| ProjectRenderError::Cancelled)?;
    let rendered = rendered.map_err(|error| {
        ProjectRenderError::Invalid(format!("canonical project JSON failed: {error:?}"))
    })?;
    let bytes = u64::try_from(rendered.len()).map_err(|_| {
        ProjectRenderError::WorkEnvelope(InvocationWorkExceeded::CumulativeBytes {
            attempted: u64::MAX,
        })
    })?;
    work.charge(bytes)
        .map_err(ProjectRenderError::WorkEnvelope)?;
    Ok(rendered)
}

fn five_explicits_match(
    row: &OpRow,
    expectations: RenderedExplicitsRef<'_>,
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<bool, EvidenceCompareError> {
    let phase = SolveEvidencePhase::FiveExplicitsCompare;
    let plan_index = Some(candidate_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    let fields: [(&[u8], &[u8]); 4] = [
        (row.seed.as_slice(), expectations.seed),
        (row.versions.as_bytes(), expectations.versions.as_bytes()),
        (row.budget.as_bytes(), expectations.budget.as_bytes()),
        (
            row.capability.as_bytes(),
            expectations.capability.as_bytes(),
        ),
    ];
    let mut compared = 0u64;
    for (retained, expected) in fields {
        if retained.len() != expected.len() {
            let next = compared.saturating_add(1);
            work.checkpoint(phase, plan_index, next)
                .map_err(|_| EvidenceCompareError::Cancelled)?;
            work.charge(1).map_err(EvidenceCompareError::WorkEnvelope)?;
            return Ok(false);
        }
        for (retained_tile, expected_tile) in retained
            .chunks(EVIDENCE_POLL_BYTES)
            .zip(expected.chunks(EVIDENCE_POLL_BYTES))
        {
            work.checkpoint(phase, plan_index, compared)
                .map_err(|_| EvidenceCompareError::Cancelled)?;
            let matches = retained_tile == expected_tile;
            let tile_len = u64::try_from(retained_tile.len()).unwrap_or(u64::MAX);
            compared = compared.saturating_add(tile_len);
            work.checkpoint(phase, plan_index, compared)
                .map_err(|_| EvidenceCompareError::Cancelled)?;
            work.charge(tile_len)
                .map_err(EvidenceCompareError::WorkEnvelope)?;
            if !matches {
                return Ok(false);
            }
        }
    }
    work.checkpoint(phase, plan_index, u64::MAX)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    Ok(true)
}

fn geometry_source_identity_controlled(
    artifact: &GeometryArtifact,
    work: EvidenceWork<'_>,
    source_index: usize,
) -> Result<String, EvidenceCompareError> {
    let phase = SolveEvidencePhase::ProjectIdentityDerive;
    let plan_index = Some(source_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    let identity = geometry_source_identity(artifact);
    work.checkpoint(phase, plan_index, 1)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    let charge = u64::try_from(identity.len())
        .ok()
        .and_then(|bytes| bytes.checked_add(DERIVATION_ITEM_WORK_BYTES))
        .ok_or(EvidenceCompareError::WorkEnvelope(
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        ))?;
    work.charge(charge)
        .map_err(EvidenceCompareError::WorkEnvelope)?;
    Ok(identity)
}

fn controlled_candidate_text(
    op: i64,
    field: &'static str,
    bytes: Vec<u8>,
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<String, CandidateOpReadError> {
    let phase = SolveEvidencePhase::CandidateOpTextConversion;
    let plan_index = Some(candidate_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| CandidateOpReadError::Cancelled)?;
    let converted = String::from_utf8(bytes);
    work.checkpoint(phase, plan_index, 1)
        .map_err(|_| CandidateOpReadError::Cancelled)?;
    converted.map_err(|_| {
        CandidateOpReadError::Ledger(LedgerError::OpCorrupt {
            op,
            detail: format!(
                "controlled {field} delivery was not UTF-8 despite the guarded text-field read"
            ),
        })
    })
}

fn read_visible_op_page(
    ledger: &Ledger,
    cursor: Option<VisibleOpCursor>,
    work: EvidenceWork<'_>,
    page_index: usize,
) -> Result<VisibleOpPage, VisibleOpScanError> {
    let phase = SolveEvidencePhase::VisibleOpPage;
    let plan_index = Some(page_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| VisibleOpScanError::Cancelled)?;
    let mut examined_rows = 0u64;
    let controlled = ledger.visible_op_ids_page_controlled(
        MAIN_BRANCH,
        cursor,
        MAX_VISIBLE_OP_PAGE_ROWS,
        |progress| {
            if work
                .checkpoint(phase, plan_index, progress.examined_rows)
                .is_err()
            {
                return ControlFlow::Break(VisibleOpScanError::Cancelled);
            }
            let Some(delta) = progress.examined_rows.checked_sub(examined_rows) else {
                return ControlFlow::Break(VisibleOpScanError::Ledger(LedgerError::Invalid {
                    field: "visible_op_page.progress".to_string(),
                    problem: "controlled visible-op progress regressed".to_string(),
                }));
            };
            let Some(charge) = delta.checked_mul(OP_ID_WORK_BYTES) else {
                return ControlFlow::Break(VisibleOpScanError::WorkEnvelope(
                    InvocationWorkExceeded::CumulativeBytes {
                        attempted: u64::MAX,
                    },
                ));
            };
            if let Err(error) = work.charge(charge) {
                return ControlFlow::Break(VisibleOpScanError::WorkEnvelope(error));
            }
            examined_rows = progress.examined_rows;
            if work
                .checkpoint(phase, plan_index, progress.examined_rows)
                .is_err()
            {
                ControlFlow::Break(VisibleOpScanError::Cancelled)
            } else {
                ControlFlow::Continue(())
            }
        },
    );
    work.checkpoint(phase, plan_index, u64::MAX)
        .map_err(|_| VisibleOpScanError::Cancelled)?;
    let controlled = controlled.map_err(VisibleOpScanError::Ledger)?;
    match controlled {
        ControlledVisibleOpPage::Break { reason, .. } => Err(reason),
        ControlledVisibleOpPage::Complete(page) => Ok(page),
    }
}

fn read_op_edges_controlled(
    ledger: &Ledger,
    op: i64,
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<BoundedOpArtifactEdges, EvidenceReadError> {
    let phase = SolveEvidencePhase::EdgePageRead;
    let plan_index = Some(candidate_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| EvidenceReadError::Cancelled)?;
    let result = ledger.op_artifact_edges_bounded(op, EDGE_SCAN_CAP);
    work.checkpoint(phase, plan_index, 1)
        .map_err(|_| EvidenceReadError::Cancelled)?;
    let page = result.map_err(EvidenceReadError::Ledger)?;
    let observed_items = page
        .edges
        .len()
        .checked_add(usize::from(page.truncated))
        .and_then(|count| u64::try_from(count).ok())
        .and_then(|count| count.checked_mul(EDGE_ITEM_WORK_BYTES))
        .ok_or(EvidenceReadError::WorkEnvelope(
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        ))?;
    work.charge(observed_items)
        .map_err(EvidenceReadError::WorkEnvelope)?;
    Ok(page)
}

fn read_artifact_info_controlled(
    ledger: &Ledger,
    artifact: &ContentHash,
    work: EvidenceWork<'_>,
    descriptor_index: usize,
) -> Result<Option<ArtifactInfo>, EvidenceReadError> {
    let phase = SolveEvidencePhase::ArtifactDescriptorRead;
    let plan_index = Some(descriptor_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| EvidenceReadError::Cancelled)?;
    let result = ledger.artifact_info(artifact);
    work.checkpoint(phase, plan_index, 1)
        .map_err(|_| EvidenceReadError::Cancelled)?;
    let info = result.map_err(EvidenceReadError::Ledger)?;
    let variable_bytes = match info.as_ref() {
        Some(info) => info
            .kind
            .len()
            .checked_add(info.meta.as_ref().map_or(0, String::len))
            .ok_or(EvidenceReadError::WorkEnvelope(
                InvocationWorkExceeded::CumulativeBytes {
                    attempted: u64::MAX,
                },
            ))?,
        None => 0,
    };
    let charge = u64::try_from(variable_bytes)
        .ok()
        .and_then(|bytes| bytes.checked_add(DERIVATION_ITEM_WORK_BYTES))
        .ok_or(EvidenceReadError::WorkEnvelope(
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        ))?;
    work.charge(charge)
        .map_err(EvidenceReadError::WorkEnvelope)?;
    Ok(info)
}

fn read_candidate_op(
    ledger: &Ledger,
    op: i64,
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<Option<ControlledCandidateOp>, CandidateOpReadError> {
    let phase = SolveEvidencePhase::CandidateOpRowRead;
    let plan_index = Some(candidate_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| CandidateOpReadError::Cancelled)?;
    let mut fields = CandidateOpFields::default();
    let mut delivered = 0u64;
    let controlled = ledger.read_op_fields_controlled(op, |field, offset, tile| {
        let callback_units = delivered.saturating_add(1);
        if work.checkpoint(phase, plan_index, callback_units).is_err() {
            return ControlFlow::Break(CandidateOpReadError::Cancelled);
        }
        let output = fields.field_mut(field);
        let expected_offset = u64::try_from(output.len()).unwrap_or(u64::MAX);
        if offset != expected_offset {
            return ControlFlow::Break(CandidateOpReadError::Ledger(LedgerError::OpCorrupt {
                op,
                detail: format!(
                    "controlled {} field resumed at byte {offset}, expected {expected_offset}",
                    field.as_str()
                ),
            }));
        }
        let tile_len = u64::try_from(tile.len()).unwrap_or(u64::MAX);
        if let Err(error) = work.charge(tile_len) {
            return ControlFlow::Break(CandidateOpReadError::WorkEnvelope(error));
        }
        if output.try_reserve(tile.len()).is_err() {
            return ControlFlow::Break(CandidateOpReadError::Ledger(LedgerError::Invalid {
                field: format!("op.{field}", field = field.as_str()),
                problem: format!(
                    "allocation for controlled operation {op} field {} was refused",
                    field.as_str()
                ),
            }));
        }
        output.extend_from_slice(tile);
        delivered = delivered.saturating_add(tile_len);
        if work
            .checkpoint(phase, plan_index, delivered.saturating_add(1))
            .is_err()
        {
            ControlFlow::Break(CandidateOpReadError::Cancelled)
        } else {
            ControlFlow::Continue(())
        }
    });
    work.checkpoint(phase, plan_index, u64::MAX)
        .map_err(|_| CandidateOpReadError::Cancelled)?;
    let controlled = controlled.map_err(CandidateOpReadError::Ledger)?;
    let complete = match controlled {
        ControlledOpRead::NotFound => return Ok(None),
        ControlledOpRead::Break(error) => return Err(error),
        ControlledOpRead::Complete(complete) => complete,
    };
    let metadata = complete.metadata;
    let exact_len = |field: &'static str,
                     expected: Option<u64>,
                     observed: usize|
     -> Result<(), CandidateOpReadError> {
        let observed = u64::try_from(observed).unwrap_or(u64::MAX);
        if expected.map_or(observed == 0, |expected| expected == observed) {
            Ok(())
        } else {
            Err(CandidateOpReadError::Ledger(LedgerError::OpCorrupt {
                op,
                detail: format!(
                    "controlled {field} delivery length {observed} disagrees with metadata \
                     length {expected:?}"
                ),
            }))
        }
    };
    exact_len("session", metadata.session_len, fields.session.len())?;
    exact_len("IR", Some(metadata.ir_len), fields.ir.len())?;
    exact_len("seed", Some(metadata.seed_len), fields.seed.len())?;
    exact_len(
        "versions",
        Some(metadata.versions_len),
        fields.versions.len(),
    )?;
    exact_len("budget", Some(metadata.budget_len), fields.budget.len())?;
    exact_len(
        "capability",
        Some(metadata.capability_len),
        fields.capability.len(),
    )?;
    exact_len(
        "diagnostic",
        metadata.diagnostic_len,
        fields.diagnostic.len(),
    )?;
    let session = metadata.session_len.map(|_| fields.session);
    let diagnostic = if metadata.diagnostic_len.is_some() {
        Some(controlled_candidate_text(
            op,
            "diagnostic",
            fields.diagnostic,
            work,
            candidate_index,
        )?)
    } else {
        None
    };
    let row = OpRow {
        id: metadata.id,
        session,
        ir: controlled_candidate_text(op, "IR", fields.ir, work, candidate_index)?,
        seed: fields.seed,
        versions: controlled_candidate_text(
            op,
            "versions",
            fields.versions,
            work,
            candidate_index,
        )?,
        budget: controlled_candidate_text(op, "budget", fields.budget, work, candidate_index)?,
        capability: controlled_candidate_text(
            op,
            "capability",
            fields.capability,
            work,
            candidate_index,
        )?,
        t_start: metadata.t_start,
        t_end: metadata.t_end,
        outcome: metadata.outcome.map(|outcome| outcome.as_str().to_string()),
        diag: diagnostic,
    };
    Ok(Some(ControlledCandidateOp {
        row,
        branch: metadata.branch,
        exec_mode: metadata.exec_mode,
        prehashed_content: complete.prehashed_content,
    }))
}

fn verify_candidate_op_identity(
    ledger: &Ledger,
    candidate: &ControlledCandidateOp,
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<(), CandidateOpReadError> {
    let phase = SolveEvidencePhase::OperationContentIdentity;
    let plan_index = Some(candidate_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| CandidateOpReadError::Cancelled)?;
    let result = ledger.verify_op_content_identity_prehashed(&candidate.prehashed_content);
    work.checkpoint(phase, plan_index, u64::MAX)
        .map_err(|_| CandidateOpReadError::Cancelled)?;
    work.charge(OP_CONTENT_IDENTITY_WORK_BYTES)
        .map_err(CandidateOpReadError::WorkEnvelope)?;
    result.map_err(CandidateOpReadError::Ledger)?;
    Ok(())
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
    work: EvidenceWork<'_>,
    expectations: RenderedExplicitsRef<'_>,
) -> Result<Option<ImportSummary>, ImportSummaryError> {
    let mut cursor = None;
    let mut page_index = 0usize;
    let mut candidate_index = 0usize;
    let mut visible_ids = 0usize;
    loop {
        if page_index >= MAX_SOLVE_VISIBLE_OP_PAGES {
            return Err(ImportSummaryError::WorkEnvelope(
                InvocationWorkExceeded::VisiblePages {
                    limit: MAX_SOLVE_VISIBLE_OP_PAGES,
                },
            ));
        }
        let page = read_visible_op_page(ledger, cursor, work, page_index)?;
        page_index = page_index.saturating_add(1);
        if page.truncated != page.next_cursor.is_some() {
            return Err(ImportSummaryError::Ledger(LedgerError::Invalid {
                field: "visible_op_page".to_string(),
                problem: "controlled page truncation and continuation disagree".to_string(),
            }));
        }
        'candidate: for id in page.op_ids {
            if visible_ids >= MAX_SOLVE_VISIBLE_OP_IDS {
                return Err(ImportSummaryError::WorkEnvelope(
                    InvocationWorkExceeded::VisibleIds {
                        limit: MAX_SOLVE_VISIBLE_OP_IDS,
                    },
                ));
            }
            visible_ids = visible_ids.saturating_add(1);
            let this_candidate = candidate_index;
            candidate_index = candidate_index.saturating_add(1);
            let Some(candidate) = read_candidate_op(ledger, id, work, this_candidate)? else {
                continue;
            };
            let attestation = match attest_import_row(
                ledger,
                spec,
                id,
                &candidate,
                work,
                this_candidate,
                expectations,
            ) {
                Ok(attestation) => attestation,
                Err(ImportSummaryError::Invalid(_)) => continue,
                Err(error @ ImportSummaryError::Cancelled)
                | Err(error @ ImportSummaryError::WorkEnvelope(_))
                | Err(error @ ImportSummaryError::Unsupported(_))
                | Err(error @ ImportSummaryError::Ledger(_)) => return Err(error),
            };
            let edges = read_op_edges_controlled(ledger, id, work, this_candidate)?;
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
            for (descriptor_index, edge) in edges.edges.iter().enumerate() {
                let Some(info) =
                    read_artifact_info_controlled(ledger, &edge.artifact, work, descriptor_index)?
                else {
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
                &edges.edges,
                work,
                this_candidate,
            ) {
                Ok(summary) => return Ok(Some(summary)),
                Err(ImportSummaryError::Invalid(_)) => {}
                Err(error @ ImportSummaryError::Cancelled)
                | Err(error @ ImportSummaryError::WorkEnvelope(_))
                | Err(error @ ImportSummaryError::Unsupported(_))
                | Err(error @ ImportSummaryError::Ledger(_)) => return Err(error),
            }
        }
        match page.next_cursor {
            Some(next) => {
                if visible_ids >= MAX_SOLVE_VISIBLE_OP_IDS {
                    return Err(ImportSummaryError::WorkEnvelope(
                        InvocationWorkExceeded::VisibleIds {
                            limit: MAX_SOLVE_VISIBLE_OP_IDS,
                        },
                    ));
                }
                cursor = Some(next);
            }
            None => return Ok(None),
        }
    }
}

/// Load and independently attest the latest sealed driver state for a run.
///
/// The legacy envelope expectation used while decoding is only a bounded
/// codec check: its byte identity is derived from the candidate bytes and
/// grants no authority. Resume eligibility comes exclusively from
/// [`validate_resume_candidate`], which re-attests the complete canonical
/// operation and checkpoint chain before a governor is opened.
fn load_latest_state(
    ledger: &Ledger,
    run: SolveRunId,
    work: EvidenceWork<'_>,
) -> Result<VerifiedResume, SolveRefusal> {
    let mut cursor = None;
    let mut page_index = 0usize;
    let mut candidate_index = 0usize;
    let mut visible_ids = 0usize;
    let mut best: Option<VerifiedResume> = None;
    let mut best_is_ambiguous = false;
    let mut import_cache = None;
    loop {
        if page_index >= MAX_SOLVE_VISIBLE_OP_PAGES {
            return Err(invocation_work_refusal(
                Some(run),
                Some(SolveStage::ImportVerify),
                InvocationWorkExceeded::VisiblePages {
                    limit: MAX_SOLVE_VISIBLE_OP_PAGES,
                },
            ));
        }
        let page =
            read_visible_op_page(ledger, cursor, work, page_index).map_err(
                |error| match error {
                    VisibleOpScanError::Cancelled => cancelled_resume_refusal(run),
                    VisibleOpScanError::WorkEnvelope(error) => {
                        invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
                    }
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
                return Err(invocation_work_refusal(
                    Some(run),
                    Some(SolveStage::ImportVerify),
                    InvocationWorkExceeded::VisibleIds {
                        limit: MAX_SOLVE_VISIBLE_OP_IDS,
                    },
                ));
            }
            visible_ids = visible_ids.saturating_add(1);
            let this_candidate = candidate_index;
            candidate_index = candidate_index.saturating_add(1);
            let Some(candidate) = read_candidate_op(ledger, id, work, this_candidate).map_err(
                |error| match error {
                    CandidateOpReadError::Cancelled => cancelled_resume_refusal(run),
                    CandidateOpReadError::WorkEnvelope(error) => {
                        invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
                    }
                    CandidateOpReadError::Ledger(error) => {
                        resume_ledger("scanning the ledger for the run failed", error)
                    }
                },
            )?
            else {
                continue;
            };
            if candidate.branch != MAIN_BRANCH
                || candidate.exec_mode != ExecMode::Deterministic
                || !is_supported_stage_discovery_row(&candidate.row, run)
            {
                continue;
            }
            let edges = resume_edges(ledger, run, id, work, this_candidate)?;
            for (descriptor_index, edge) in edges.iter().enumerate() {
                if edge.role != EdgeRole::Out {
                    continue;
                }
                let Some(info) =
                    read_artifact_info_controlled(ledger, &edge.artifact, work, descriptor_index)
                        .map_err(|error| match error {
                        EvidenceReadError::Cancelled => cancelled_resume_refusal(run),
                        EvidenceReadError::WorkEnvelope(error) => invocation_work_refusal(
                            Some(run),
                            Some(SolveStage::ImportVerify),
                            error,
                        ),
                        EvidenceReadError::Ledger(error) => resume_ledger(
                            "reading a retained driver-state descriptor failed",
                            error,
                        ),
                    })?
                else {
                    continue;
                };
                if info.kind != STAGE_STATE_KIND {
                    continue;
                }
                let state = decode_driver_state(ledger, run, edge.artifact, work)?;
                validate_state_shape(&state, run)?;
                if best
                    .as_ref()
                    .is_some_and(|existing| existing.state.completed.len() > state.completed.len())
                {
                    continue;
                }
                let verified = validate_resume_candidate(
                    ledger,
                    run,
                    state,
                    edge.artifact,
                    id,
                    work,
                    &mut import_cache,
                )?;
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
        match page.next_cursor {
            Some(next) => {
                if visible_ids >= MAX_SOLVE_VISIBLE_OP_IDS {
                    return Err(invocation_work_refusal(
                        Some(run),
                        Some(SolveStage::ImportVerify),
                        InvocationWorkExceeded::VisibleIds {
                            limit: MAX_SOLVE_VISIBLE_OP_IDS,
                        },
                    ));
                }
                cursor = Some(next);
            }
            None => break,
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
    drop(import_cache);
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
    work: EvidenceWork<'_>,
    import_cache: &mut Option<Rc<ResumeImportCache>>,
) -> Result<VerifiedResume, SolveRefusal> {
    validate_state_shape(&state, run)?;
    let first = state
        .completed
        .first()
        .expect("state shape requires a completed prefix");
    let first_edges = resume_edges(ledger, run, first.op_id, work, 0)?;
    require_stage_edge_seal(ledger, run, first.op_id, first_edges.len(), work, 0)?;
    let (project_source, retained_source) =
        read_retained_project_source(ledger, run, &first_edges, work)?;
    let attestation = if let Some(cached) = import_cache
        .as_ref()
        .filter(|cached| cached.source_hash == project_source)
    {
        let canonical_matches = evidence_bytes_equal(
            cached.project.canonical.as_bytes(),
            retained_source.as_bytes(),
            work,
            SolveEvidencePhase::ResumeProjectCanonicalCompare,
            None,
        )
        .map_err(|error| match error {
            EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
            EvidenceCompareError::WorkEnvelope(error) => {
                invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
            }
        })?;
        if !canonical_matches {
            return Err(resume_identity(
                "the retained project source is not the cached canonical project",
            ));
        }
        Rc::clone(cached)
    } else {
        let attested = Rc::new(attest_retained_project(
            run,
            project_source,
            &retained_source,
            work,
        )?);
        if import_cache.is_none() {
            *import_cache = Some(Rc::clone(&attested));
        }
        attested
    };
    let project = &attestation.project;
    let project_hash = attestation.project_hash;
    let expectations = attestation.expectations();
    if state.project != *project_hash.as_bytes() {
        return Err(resume_identity(
            "the retained driver state carries a different project identity",
        ));
    }

    let mut context = StageContext::default();
    let mut recovered_cards: Option<CardPackSet> = None;
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
        let candidate = read_candidate_op(ledger, completed.op_id, work, index)
            .map_err(|error| match error {
                CandidateOpReadError::Cancelled => cancelled_resume_refusal(run),
                CandidateOpReadError::WorkEnvelope(error) => {
                    invocation_work_refusal(Some(run), Some(stage), error)
                }
                CandidateOpReadError::Ledger(error) => {
                    resume_ledger("reading a completed solve operation failed", error)
                }
            })?
            .ok_or_else(|| {
                resume_identity(format!(
                    "completed stage {index} names missing operation {}",
                    completed.op_id
                ))
            })?;
        validate_stage_row(
            ledger,
            &candidate,
            StageRowExpectation {
                stage,
                run,
                project_hash,
                explicits: expectations,
            },
            work,
            index,
        )?;
        let edges = if index == 0 {
            first_edges.clone()
        } else {
            resume_edges(ledger, run, completed.op_id, work, index)?
        };
        if index > 0 {
            require_stage_edge_seal(ledger, run, completed.op_id, edges.len(), work, index)?;
        }
        require_artifact_kind_resume(
            ledger,
            run,
            completed.receipt,
            STAGE_RECEIPT_KIND,
            "stage receipt",
            work,
            index,
        )?;
        let has_receipt =
            has_edge_controlled(&edges, EdgeRole::Out, completed.receipt, work, index).map_err(
                |error| match error {
                    EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
                    EvidenceCompareError::WorkEnvelope(error) => {
                        invocation_work_refusal(Some(run), Some(stage), error)
                    }
                },
            )?;
        if !has_receipt {
            return Err(resume_identity(format!(
                "stage `{}` receipt {} is not an Out edge of operation {}",
                stage.name(),
                completed.receipt.to_hex(),
                completed.op_id
            )));
        }

        let checkpoint_outputs =
            artifacts_with_kind_resume(ledger, run, &edges, EdgeRole::Out, STAGE_STATE_KIND, work)?;
        if checkpoint_outputs.len() != 1 {
            return Err(resume_identity(format!(
                "stage `{}` operation {} has {} checkpoint outputs; exactly one is required",
                stage.name(),
                completed.op_id,
                checkpoint_outputs.len()
            )));
        }
        let checkpoint_hash = checkpoint_outputs[0];
        let checkpoint = decode_driver_state(ledger, run, checkpoint_hash, work)?;
        validate_checkpoint_prefix(&checkpoint, &state, index, predecessor_checkpoint.as_ref())?;

        if let Some(predecessor) = predecessor_state {
            let has_predecessor =
                has_edge_controlled(&edges, EdgeRole::In, predecessor, work, index).map_err(
                    |error| match error {
                        EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
                        EvidenceCompareError::WorkEnvelope(error) => {
                            invocation_work_refusal(Some(run), Some(stage), error)
                        }
                    },
                )?;
            if !has_predecessor {
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
            run,
            completed.receipt,
            MAX_RECEIPT_READ_BYTES,
            "stage receipt",
            work,
            SolveEvidencePhase::ResumeStageReceiptRead,
            SolveEvidencePhase::ResumeStageReceiptUtf8,
        )?;
        let mut expected_edges = vec![
            (EdgeRole::Out, completed.receipt),
            (EdgeRole::Out, checkpoint_hash),
        ];
        match stage {
            SolveStage::ImportVerify => {
                expected_edges.push((EdgeRole::In, project_source));
                // Card packs are retained against this first operation, so
                // recovery happens here regardless of which later stage was
                // interrupted.
                let (packs, pack_artifacts) = recover_card_packs_resume(ledger, run, &edges, work)?;
                // The complete run preimage is only available now. Deriving
                // it here is what turns pack recovery into an attestation:
                // any added, dropped, or substituted pack moves the identity
                // and cannot present itself as this run.
                let rederived =
                    SolveRunId::derive_with_project_hash(project, project_hash, packs.root());
                if rederived != run {
                    return Err(resume_identity(format!(
                        "the retained project and card packs re-derive run `{}` but resume \
                         requested `{}`",
                        rederived.to_hex(),
                        run.to_hex()
                    )));
                }
                for artifact in pack_artifacts {
                    expected_edges.push((EdgeRole::In, artifact));
                }
                recovered_cards = Some(packs);
                let summary_inputs = artifacts_with_kind_resume(
                    ledger,
                    run,
                    &edges,
                    EdgeRole::In,
                    IMPORT_SUMMARY_KIND,
                    work,
                )?;
                if summary_inputs.len() != 1 {
                    return Err(resume_identity(format!(
                        "import-verify operation {} has {} geometry-import summary inputs; exactly one is required",
                        completed.op_id,
                        summary_inputs.len()
                    )));
                }
                let summary_hash = summary_inputs[0];
                let (import_op, _receipt_entries) = parse_import_verify_receipt(
                    &receipt_text,
                    run,
                    project_hash,
                    work,
                )
                .map_err(|error| match error {
                    ImportSummaryError::Cancelled => cancelled_resume_refusal(run),
                    ImportSummaryError::WorkEnvelope(error) => {
                        invocation_work_refusal(Some(run), Some(stage), error)
                    }
                    ImportSummaryError::Invalid(what) => {
                        resume_identity(format!("invalid import-verify receipt: {what}"))
                    }
                    ImportSummaryError::Ledger(error) => {
                        resume_ledger("parsing retained import-verify receipt failed", error)
                    }
                    ImportSummaryError::Unsupported(what) => resume_import_envelope(run, what),
                })?;
                let summary = validate_import_candidate(
                    ledger,
                    &project.spec,
                    project_hash,
                    ImportCandidateLocator {
                        op: import_op,
                        summary_artifact: summary_hash,
                    },
                    work,
                    index,
                    expectations,
                )
                .map_err(|error| match error {
                    ImportSummaryError::Cancelled => cancelled_resume_refusal(run),
                    ImportSummaryError::WorkEnvelope(error) => {
                        invocation_work_refusal(Some(run), Some(stage), error)
                    }
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
                let receipt_stage_index = Some(index);
                work.checkpoint(
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                    0,
                )
                .map_err(|_| cancelled_resume_refusal(run))?;
                let expected_receipt =
                    import_verify_receipt(run, project_hash, summary.op_id, &summary.entries, work)
                        .map_err(|error| match error {
                            EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
                            EvidenceCompareError::WorkEnvelope(error) => {
                                invocation_work_refusal(Some(run), Some(stage), error)
                            }
                        })?;
                work.checkpoint(
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                    1,
                )
                .map_err(|_| cancelled_resume_refusal(run))?;
                let receipt_matches = evidence_bytes_equal(
                    receipt_text.as_bytes(),
                    expected_receipt.as_bytes(),
                    work,
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                )
                .map_err(|error| match error {
                    EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
                    EvidenceCompareError::WorkEnvelope(error) => {
                        invocation_work_refusal(Some(run), Some(stage), error)
                    }
                })?;
                if !receipt_matches {
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
                let receipt_stage_index = Some(index);
                work.checkpoint(
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                    0,
                )
                .map_err(|_| cancelled_resume_refusal(run))?;
                let expected_receipt = assignment_receipt(&project.spec, &context, run, work, true);
                work.checkpoint(
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                    1,
                )
                .map_err(|_| cancelled_resume_refusal(run))?;
                let expected_receipt = match expected_receipt {
                    Ok(receipt) => receipt,
                    Err(error)
                        if matches!(
                            error.code,
                            "cli-solve-cancelled" | "cli-solve-work-envelope"
                        ) =>
                    {
                        return Err(error);
                    }
                    Err(error) => {
                        return Err(resume_identity(format!(
                            "the retained assignment context cannot be reconstructed: {}",
                            error.what
                        )));
                    }
                };
                let receipt_matches = evidence_bytes_equal(
                    receipt_text.as_bytes(),
                    expected_receipt.as_bytes(),
                    work,
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                )
                .map_err(|error| match error {
                    EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
                    EvidenceCompareError::WorkEnvelope(error) => {
                        invocation_work_refusal(Some(run), Some(stage), error)
                    }
                })?;
                if !receipt_matches {
                    return Err(resume_identity(
                        "the retained assignment receipt is not the canonical driver receipt",
                    ));
                }
                let predecessor = predecessor_state.ok_or_else(|| {
                    resume_identity(
                        "the retained assign stage has no verified import-stage predecessor",
                    )
                })?;
                expected_edges.push((EdgeRole::In, predecessor));
            }
            SolveStage::MaterialResolve => {
                let cards = recovered_cards.as_ref().ok_or_else(|| {
                    resume_identity(
                        "the retained material-resolve stage has no recovered card-pack set",
                    )
                })?;
                let receipt_stage_index = Some(index);
                work.checkpoint(
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                    0,
                )
                .map_err(|_| cancelled_resume_refusal(run))?;
                let rebuilt = material_resolve_receipt(&project.spec, cards, run, work, true);
                work.checkpoint(
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                    1,
                )
                .map_err(|_| cancelled_resume_refusal(run))?;
                let (expected_receipt, usages) = match rebuilt {
                    Ok(rebuilt) => rebuilt,
                    Err(error)
                        if matches!(
                            error.code,
                            "cli-solve-cancelled" | "cli-solve-work-envelope"
                        ) =>
                    {
                        return Err(error);
                    }
                    Err(error) => {
                        return Err(resume_identity(format!(
                            "the retained material bindings no longer resolve against the \
                             recovered card packs: {}",
                            error.what
                        )));
                    }
                };
                let receipt_matches = evidence_bytes_equal(
                    receipt_text.as_bytes(),
                    expected_receipt.as_bytes(),
                    work,
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                )
                .map_err(|error| match error {
                    EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
                    EvidenceCompareError::WorkEnvelope(error) => {
                        invocation_work_refusal(Some(run), Some(stage), error)
                    }
                })?;
                if !receipt_matches {
                    return Err(resume_identity(
                        "the retained material-resolve receipt is not the canonical driver receipt",
                    ));
                }
                for pack in cards.iter() {
                    expected_edges.push((EdgeRole::In, pack.artifact()));
                }
                for usage in &usages {
                    expected_edges.push((EdgeRole::Out, usage.artifact));
                }
                let predecessor = predecessor_state.ok_or_else(|| {
                    resume_identity(
                        "the retained material-resolve stage has no verified assign-stage \
                         predecessor",
                    )
                })?;
                expected_edges.push((EdgeRole::In, predecessor));
            }
            SolveStage::FlowNetwork => {
                let receipt_stage_index = Some(index);
                work.checkpoint(
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                    0,
                )
                .map_err(|_| cancelled_resume_refusal(run))?;
                let rebuilt = flow_network_receipt(&project.spec, run, work, true);
                work.checkpoint(
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                    1,
                )
                .map_err(|_| cancelled_resume_refusal(run))?;
                let expected_receipt = match rebuilt {
                    Ok(rebuilt) => rebuilt,
                    Err(error)
                        if matches!(
                            error.code,
                            "cli-solve-cancelled" | "cli-solve-work-envelope"
                        ) =>
                    {
                        return Err(error);
                    }
                    Err(error) => {
                        return Err(resume_identity(format!(
                            "the retained cooling declaration no longer lowers or solves: {}",
                            error.what
                        )));
                    }
                };
                let receipt_matches = evidence_bytes_equal(
                    receipt_text.as_bytes(),
                    expected_receipt.as_bytes(),
                    work,
                    SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
                    receipt_stage_index,
                )
                .map_err(|error| match error {
                    EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
                    EvidenceCompareError::WorkEnvelope(error) => {
                        invocation_work_refusal(Some(run), Some(stage), error)
                    }
                })?;
                if !receipt_matches {
                    return Err(resume_identity(
                        "the retained flow-network receipt is not the canonical driver receipt",
                    ));
                }
                let predecessor = predecessor_state.ok_or_else(|| {
                    resume_identity(
                        "the retained flow-network stage has no verified material-resolve-stage \
                         predecessor",
                    )
                })?;
                expected_edges.push((EdgeRole::In, predecessor));
            }
            _ => unreachable!("completed unavailable stages were refused above"),
        }
        require_exact_edges(
            run,
            stage,
            completed.op_id,
            &edges,
            &expected_edges,
            work,
            index,
        )?;

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

    let cards = recovered_cards.ok_or_else(|| {
        resume_identity(
            "the candidate has no re-attested card-pack set; every run retains its inputs against \
             its first operation",
        )
    })?;
    Ok(VerifiedResume {
        state,
        state_artifact,
        attestation,
        context,
        cards,
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

fn validate_stage_row(
    ledger: &Ledger,
    candidate: &ControlledCandidateOp,
    expected: StageRowExpectation<'_>,
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<(), SolveRefusal> {
    let StageRowExpectation {
        stage,
        run,
        project_hash,
        explicits,
    } = expected;
    let row = &candidate.row;
    let ordinal = stage.ordinal();
    if row.id <= 0
        || row.session.as_deref() != Some(run.as_bytes().as_slice())
        || row.outcome.as_deref() != Some("ok")
        || row.diag.is_some()
        || row.t_start != i64::from(ordinal) * 2
        || row.t_end != Some(i64::from(ordinal) * 2 + 1)
        || row.ir != solve_stage_ir(stage, run, project_hash)
    {
        return Err(resume_identity(format!(
            "operation {} does not match canonical stage `{}` semantics",
            row.id,
            stage.name()
        )));
    }
    let explicits_match = five_explicits_match(row, explicits, work, candidate_index).map_err(
        |error| match error {
            EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
            EvidenceCompareError::WorkEnvelope(error) => {
                invocation_work_refusal(Some(run), Some(stage), error)
            }
        },
    )?;
    if !explicits_match {
        return Err(resume_identity(format!(
            "operation {} does not match canonical stage `{}` Five Explicits",
            row.id,
            stage.name()
        )));
    }
    if candidate.branch != MAIN_BRANCH || candidate.exec_mode != ExecMode::Deterministic {
        return Err(resume_identity(format!(
            "operation {} is not a deterministic main-branch stage operation",
            row.id
        )));
    }
    verify_candidate_op_identity(ledger, candidate, work, candidate_index).map_err(|error| {
        match error {
            CandidateOpReadError::Cancelled => cancelled_resume_refusal(run),
            CandidateOpReadError::WorkEnvelope(error) => {
                invocation_work_refusal(Some(run), Some(stage), error)
            }
            CandidateOpReadError::Ledger(error) => {
                resume_ledger("validating solve operation content identity failed", error)
            }
        }
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

fn resume_edges(
    ledger: &Ledger,
    run: SolveRunId,
    op: i64,
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<Vec<OpArtifactEdge>, SolveRefusal> {
    let edges =
        read_op_edges_controlled(ledger, op, work, candidate_index).map_err(
            |error| match error {
                EvidenceReadError::Cancelled => cancelled_resume_refusal(run),
                EvidenceReadError::WorkEnvelope(error) => {
                    invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
                }
                EvidenceReadError::Ledger(error) => {
                    resume_ledger("reading solve operation edges failed", error)
                }
            },
        )?;
    if edges.truncated {
        return Err(resume_identity(format!(
            "operation {op} exceeds the complete {EDGE_SCAN_CAP}-edge resume scan"
        )));
    }
    Ok(edges.edges)
}

fn require_stage_edge_seal(
    ledger: &Ledger,
    run: SolveRunId,
    op: i64,
    edge_count: usize,
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<(), SolveRefusal> {
    let phase = SolveEvidencePhase::EdgeSealRead;
    let plan_index = Some(candidate_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| cancelled_resume_refusal(run))?;
    let seal = ledger.op_artifact_edge_seal(op);
    work.checkpoint(phase, plan_index, 1)
        .map_err(|_| cancelled_resume_refusal(run))?;
    work.charge(DERIVATION_ITEM_WORK_BYTES).map_err(|error| {
        invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
    })?;
    let seal =
        seal.map_err(|error| resume_ledger("validating solve operation edge seal failed", error))?;
    if seal != Some(edge_count) {
        return Err(resume_identity(format!(
            "operation {op} lacks the driver's exact {edge_count}-edge lineage seal"
        )));
    }
    Ok(())
}

fn read_retained_project_source(
    ledger: &Ledger,
    run: SolveRunId,
    edges: &[OpArtifactEdge],
    work: EvidenceWork<'_>,
) -> Result<(ContentHash, String), SolveRefusal> {
    let mut sources = Vec::new();
    for (descriptor_index, edge) in edges.iter().enumerate() {
        if edge.role != EdgeRole::In {
            continue;
        }
        let info = read_artifact_info_controlled(ledger, &edge.artifact, work, descriptor_index)
            .map_err(|error| match error {
                EvidenceReadError::Cancelled => cancelled_resume_refusal(run),
                EvidenceReadError::WorkEnvelope(error) => {
                    invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
                }
                EvidenceReadError::Ledger(error) => {
                    resume_ledger("reading the retained project failed", error)
                }
            })?
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
    let bytes = materialize_evidence_artifact(
        ledger,
        work,
        source_hash,
        crate::MAX_PROJECT_BYTES,
        SolveEvidencePhase::ResumeProjectRead,
        None,
    )
    .map_err(|error| match error {
        EvidenceReadError::Cancelled => cancelled_resume_refusal(run),
        EvidenceReadError::WorkEnvelope(error) => {
            invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
        }
        EvidenceReadError::Ledger(error) => {
            resume_ledger("reading the retained project failed", error)
        }
    })?
    .ok_or_else(|| resume_identity("the retained project source is missing"))?;
    let source = evidence_utf8_string(
        &bytes,
        work,
        SolveEvidencePhase::ResumeProjectUtf8,
        None,
        "retained project source",
    )
    .map_err(|error| match error {
        EvidenceUtf8Error::Cancelled => cancelled_resume_refusal(run),
        EvidenceUtf8Error::WorkEnvelope(error) => {
            invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
        }
        EvidenceUtf8Error::Invalid(problem) => resume_identity(problem),
    })?;
    Ok((source_hash, source))
}

fn attest_retained_project(
    run: SolveRunId,
    source_hash: ContentHash,
    source: &str,
    work: EvidenceWork<'_>,
) -> Result<ResumeImportCache, SolveRefusal> {
    work.checkpoint(SolveEvidencePhase::ResumeProjectParse, None, 0)
        .map_err(|_| cancelled_resume_refusal(run))?;
    let project = fs_project::parse_sexpr_migrating(source).map(|migrated| migrated.decoded);
    work.checkpoint(SolveEvidencePhase::ResumeProjectParse, None, 1)
        .map_err(|_| cancelled_resume_refusal(run))?;
    let project = project.map_err(|error| {
        resume_identity(format!(
            "the retained project source no longer parses strictly: {} ({})",
            error.code, error.detail
        ))
    })?;
    work.checkpoint(SolveEvidencePhase::ProjectValidation, None, 0)
        .map_err(|_| cancelled_resume_refusal(run))?;
    let findings = project.findings();
    work.checkpoint(SolveEvidencePhase::ProjectValidation, None, 1)
        .map_err(|_| cancelled_resume_refusal(run))?;
    work.charge(u64::try_from(project.canonical.len()).map_err(|_| {
        invocation_work_refusal(
            Some(run),
            Some(SolveStage::ImportVerify),
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        )
    })?)
    .map_err(|error| invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error))?;
    if !findings.is_empty() {
        return Err(resume_identity(
            "the retained project source has validation findings",
        ));
    }
    let canonical_matches = evidence_bytes_equal(
        project.canonical.as_bytes(),
        source.as_bytes(),
        work,
        SolveEvidencePhase::ResumeProjectCanonicalCompare,
        None,
    )
    .map_err(|error| match error {
        EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
        EvidenceCompareError::WorkEnvelope(error) => {
            invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
        }
    })?;
    if !canonical_matches {
        return Err(resume_identity(
            "the retained project source is not the exact canonical project",
        ));
    }
    work.checkpoint(SolveEvidencePhase::ProjectIdentityDerive, None, 0)
        .map_err(|_| cancelled_resume_refusal(run))?;
    let project_hash = project.hash();
    work.checkpoint(SolveEvidencePhase::ProjectIdentityDerive, None, 1)
        .map_err(|_| cancelled_resume_refusal(run))?;
    let identity_bytes = u64::try_from(project.canonical.len()).map_err(|_| {
        invocation_work_refusal(
            Some(run),
            Some(SolveStage::ImportVerify),
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        )
    })?;
    work.charge(identity_bytes).map_err(|error| {
        invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
    })?;
    // Run identity now binds the admitted card-pack set as well as the
    // project, so it cannot be re-derived here: the packs are recovered from
    // the run's first operation, which the caller re-attests. The complete
    // re-derivation lives in that stage's arm and always executes, because a
    // checkpoint exists only after the first stage completes.
    work.checkpoint(SolveEvidencePhase::FiveExplicitsRender, None, 0)
        .map_err(|_| cancelled_resume_refusal(run))?;
    let rendered_explicits = explicits(&project.spec);
    work.checkpoint(SolveEvidencePhase::FiveExplicitsRender, None, 1)
        .map_err(|_| cancelled_resume_refusal(run))?;
    let (versions, budget, capability, seed) = rendered_explicits.map_err(|error| {
        resume_identity(format!(
            "the retained project cannot reproduce the stage Five Explicits: {}",
            error.what
        ))
    })?;
    let explicit_bytes = versions
        .len()
        .checked_add(budget.len())
        .and_then(|bytes| bytes.checked_add(capability.len()))
        .and_then(|bytes| bytes.checked_add(seed.len()))
        .and_then(|bytes| u64::try_from(bytes).ok())
        .ok_or_else(|| {
            invocation_work_refusal(
                Some(run),
                Some(SolveStage::ImportVerify),
                InvocationWorkExceeded::CumulativeBytes {
                    attempted: u64::MAX,
                },
            )
        })?;
    work.charge(explicit_bytes).map_err(|error| {
        invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
    })?;
    let canonical_project_json =
        render_canonical_project_json(&project.spec, work).map_err(|error| match error {
            ProjectRenderError::Cancelled => cancelled_resume_refusal(run),
            ProjectRenderError::WorkEnvelope(error) => {
                invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
            }
            ProjectRenderError::Invalid(problem) => resume_identity(problem),
        })?;
    Ok(ResumeImportCache {
        source_hash,
        project,
        project_hash,
        versions,
        budget,
        capability,
        seed,
        canonical_project_json,
    })
}

fn resume_identity(what: impl Into<String>) -> SolveRefusal {
    SolveRefusal::plain(
        "cli-solve-resume-identity",
        what,
        "verify ledger integrity; codec validity alone cannot authorize resume",
    )
}

fn cancelled_before_run_refusal() -> SolveRefusal {
    SolveRefusal::plain(
        "cli-solve-cancelled",
        "cancellation observed while validating the project before solve; no solve publication \
         was made",
        "retry the fresh solve with a new cancellation gate",
    )
}

fn cancelled_fresh_refusal(run: SolveRunId, stage: Option<SolveStage>) -> SolveRefusal {
    SolveRefusal {
        code: "cli-solve-cancelled",
        stage: stage.map(SolveStage::name),
        what: "cancellation observed during bounded solve derivation; no new solve publication \
               was made"
            .to_string(),
        fix: "retry the fresh solve with a new cancellation gate".to_string(),
        dependency: None,
        run: Some(run.to_hex()),
        recorded_op: None,
    }
}

fn invocation_work_refusal(
    run: Option<SolveRunId>,
    stage: Option<SolveStage>,
    exceeded: InvocationWorkExceeded,
) -> SolveRefusal {
    let what = match exceeded {
        InvocationWorkExceeded::CumulativeBytes { attempted } => format!(
            "solve discovery and re-attestation attempted {attempted} cumulative work bytes above \
             the invocation envelope {MAX_SOLVE_INVOCATION_WORK_BYTES}"
        ),
        InvocationWorkExceeded::VisibleIds { limit } => format!(
            "solve discovery reached the invocation ceiling of {limit} visible operation ids \
             while additional frozen-history rows remained"
        ),
        InvocationWorkExceeded::VisiblePages { limit } => format!(
            "solve discovery reached the invocation ceiling of {limit} visible-operation pages \
             while additional frozen-history rows remained"
        ),
    };
    SolveRefusal {
        code: "cli-solve-work-envelope",
        stage: stage.map(SolveStage::name),
        what,
        fix: "compact or split unrelated ledger history, or reduce retained candidate evidence, \
              then retry; the driver will not report a false not-found result"
            .to_string(),
        dependency: None,
        run: run.map(|run| run.to_hex()),
        recorded_op: None,
    }
}

fn cancelled_resume_refusal(run: SolveRunId) -> SolveRefusal {
    SolveRefusal {
        code: "cli-solve-cancelled",
        stage: Some(SolveStage::ImportVerify.name()),
        what: "cancellation observed while re-attesting retained solve evidence; no new solve \
               publication was made"
            .to_string(),
        fix: format!(
            "retry `frankensim solve --resume {} <ledger>` with a fresh cancellation gate",
            run.to_hex()
        ),
        dependency: None,
        run: Some(run.to_hex()),
        recorded_op: None,
    }
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
    work: EvidenceWork<'_>,
) -> Result<SolveDriverState, SolveRefusal> {
    let bytes = materialize_evidence_artifact(
        ledger,
        work,
        artifact,
        MAX_STATE_ENVELOPE_BYTES,
        SolveEvidencePhase::ResumeStateRead,
        None,
    )
    .map_err(|error| match error {
        EvidenceReadError::Cancelled => cancelled_resume_refusal(run),
        EvidenceReadError::WorkEnvelope(error) => {
            invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
        }
        EvidenceReadError::Ledger(error) => {
            resume_ledger("reading a retained driver checkpoint failed", error)
        }
    })?
    .ok_or_else(|| {
        resume_identity(format!(
            "retained driver checkpoint {} is missing",
            artifact.to_hex()
        ))
    })?;
    // The completing controlled ledger read above already validates and
    // supplies the exact content identity. Parsing that fixed 32-byte identity
    // avoids a second whole-envelope hash pass; a callback stop returns before
    // this decode path.
    // This self-derived expectation is deliberately only a bounded codec and
    // corruption check. `validate_resume_candidate` supplies the independent
    // semantic/lineage admission that makes the decoded value usable.
    let expectation = LegacySnapshotExpectationV1::new(
        ContentId::parse_slice(artifact.as_bytes()).expect("ledger content hash is 32 bytes"),
        DRIVER_STATE_TYPE_ID_V1,
        DRIVER_STATE_SCHEMA_VERSION_V1,
        envelope_provenance(run),
    );
    let limits = LegacySnapshotLimitsV1::new(MAX_STATE_ENVELOPE_BYTES, STATE_HASH_POLL_BYTES);
    let mut decode_polls = 0u64;
    work.checkpoint(SolveEvidencePhase::ResumeStateDecode, None, 0)
        .map_err(|_| cancelled_resume_refusal(run))?;
    let opened = LegacySnapshotV1Adapter::<SolveDriverState>::open_expected(
        &bytes,
        expectation,
        limits,
        || {
            let units = decode_polls.saturating_mul(u64::from(STATE_HASH_POLL_BYTES));
            decode_polls = decode_polls.saturating_add(1);
            work.checkpoint(SolveEvidencePhase::ResumeStateDecode, None, units)
                .is_err()
        },
    );
    let decoded_units = u64::try_from(bytes.len()).map_err(|_| {
        resume_identity("the retained driver state length is outside the decode checkpoint range")
    })?;
    work.checkpoint(SolveEvidencePhase::ResumeStateDecode, None, decoded_units)
        .map_err(|_| cancelled_resume_refusal(run))?;
    let opened = opened.map_err(|error| {
        if matches!(error, LegacySnapshotV1Error::Cancelled { .. }) {
            cancelled_resume_refusal(run)
        } else {
            resume_identity(format!(
                "the retained driver state failed bounded envelope admission: {error:?}"
            ))
        }
    })?;
    Ok(opened.into_parts().0)
}

fn require_artifact_kind_resume(
    ledger: &Ledger,
    run: SolveRunId,
    artifact: ContentHash,
    expected_kind: &str,
    label: &str,
    work: EvidenceWork<'_>,
    descriptor_index: usize,
) -> Result<(), SolveRefusal> {
    let info = read_artifact_info_controlled(ledger, &artifact, work, descriptor_index)
        .map_err(|error| match error {
            EvidenceReadError::Cancelled => cancelled_resume_refusal(run),
            EvidenceReadError::WorkEnvelope(error) => {
                invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
            }
            EvidenceReadError::Ledger(error) => {
                resume_ledger("reading an artifact descriptor failed", error)
            }
        })?
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

/// Recover the run's admitted card packs from the retained inputs of its
/// first operation and re-admit them through the fresh path's canonicalizer.
///
/// This is recovery, not authority. What makes it an attestation is the
/// caller's follow-up check that the recovered set root reproduces the
/// requested run identity: adding, dropping, or substituting a pack moves
/// that identity, so a tampered input set cannot present itself as this run.
fn recover_card_packs_resume(
    ledger: &Ledger,
    run: SolveRunId,
    edges: &[OpArtifactEdge],
    work: EvidenceWork<'_>,
) -> Result<(CardPackSet, Vec<ContentHash>), SolveRefusal> {
    let mut builder = CardPackSetBuilder::new();
    let mut artifacts = Vec::new();
    for kind in [CardPackKind::Material, CardPackKind::Interface] {
        let retained = artifacts_with_kind_resume(
            ledger,
            run,
            edges,
            EdgeRole::In,
            kind.artifact_kind(),
            work,
        )?;
        for (index, artifact) in retained.into_iter().enumerate() {
            let bytes = materialize_evidence_artifact(
                ledger,
                work,
                artifact,
                MAX_CARD_PACK_READ_BYTES,
                SolveEvidencePhase::ResumeCardPackRead,
                Some(index),
            )
            .map_err(|error| match error {
                EvidenceReadError::Cancelled => cancelled_resume_refusal(run),
                EvidenceReadError::WorkEnvelope(error) => {
                    invocation_work_refusal(Some(run), Some(SolveStage::MaterialResolve), error)
                }
                EvidenceReadError::Ledger(error) => {
                    resume_ledger("reading a retained card pack failed", error)
                }
            })?
            .ok_or_else(|| {
                resume_identity(format!(
                    "the run's retained {} pack {} exceeds the {MAX_CARD_PACK_READ_BYTES}-byte \
                     read envelope or is missing",
                    kind.label(),
                    artifact.to_hex()
                ))
            })?;
            work.checkpoint(SolveEvidencePhase::ResumeCardPackDecode, Some(index), 0)
                .map_err(|_| cancelled_resume_refusal(run))?;
            builder
                .push(RawCardPack {
                    kind,
                    source: format!("ledger artifact {}", artifact.to_hex()),
                    bytes,
                    expect: None,
                })
                .map_err(|refusal| {
                    resume_identity(format!(
                        "a retained card pack no longer re-admits: {}",
                        refusal.what
                    ))
                })?;
            work.checkpoint(SolveEvidencePhase::ResumeCardPackDecode, Some(index), 1)
                .map_err(|_| cancelled_resume_refusal(run))?;
            artifacts.push(artifact);
        }
    }
    let set = builder.finish().map_err(|refusal| {
        resume_identity(format!(
            "the retained card-pack set no longer canonicalizes: {}",
            refusal.what
        ))
    })?;
    if set.len() != artifacts.len() {
        return Err(resume_identity(format!(
            "the run retains {} card-pack artifacts but only {} distinct packs; a retained \
             input set cannot contain duplicates",
            artifacts.len(),
            set.len()
        )));
    }
    Ok((set, artifacts))
}

fn artifacts_with_kind_resume(
    ledger: &Ledger,
    run: SolveRunId,
    edges: &[OpArtifactEdge],
    role: EdgeRole,
    kind: &str,
    work: EvidenceWork<'_>,
) -> Result<Vec<ContentHash>, SolveRefusal> {
    let mut matches = Vec::new();
    for (descriptor_index, edge) in edges.iter().enumerate() {
        if edge.role != role {
            continue;
        }
        let info = read_artifact_info_controlled(ledger, &edge.artifact, work, descriptor_index)
            .map_err(|error| match error {
                EvidenceReadError::Cancelled => cancelled_resume_refusal(run),
                EvidenceReadError::WorkEnvelope(error) => {
                    invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
                }
                EvidenceReadError::Ledger(error) => {
                    resume_ledger("reading an artifact descriptor failed", error)
                }
            })?
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
    run: SolveRunId,
    artifact: ContentHash,
    cap: u64,
    label: &str,
    work: EvidenceWork<'_>,
    read_phase: SolveEvidencePhase,
    text_phase: SolveEvidencePhase,
) -> Result<String, SolveRefusal> {
    let bytes = materialize_evidence_artifact(ledger, work, artifact, cap, read_phase, None)
        .map_err(|error| match error {
            EvidenceReadError::Cancelled => cancelled_resume_refusal(run),
            EvidenceReadError::WorkEnvelope(error) => {
                invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
            }
            EvidenceReadError::Ledger(error) => {
                resume_ledger(&format!("reading the retained {label} failed"), error)
            }
        })?
        .ok_or_else(|| {
            resume_identity(format!(
                "the retained {label} {} is missing",
                artifact.to_hex()
            ))
        })?;
    evidence_utf8_string(&bytes, work, text_phase, None, label).map_err(|error| match error {
        EvidenceUtf8Error::Cancelled => cancelled_resume_refusal(run),
        EvidenceUtf8Error::WorkEnvelope(error) => {
            invocation_work_refusal(Some(run), Some(SolveStage::ImportVerify), error)
        }
        EvidenceUtf8Error::Invalid(problem) => resume_identity(format!(
            "the retained {label} {} is invalid: {problem}",
            artifact.to_hex()
        )),
    })
}

fn has_edge_controlled(
    edges: &[OpArtifactEdge],
    role: EdgeRole,
    artifact: ContentHash,
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<bool, EvidenceCompareError> {
    let phase = SolveEvidencePhase::EdgeSetCompare;
    let plan_index = Some(candidate_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    for (index, edge) in edges.iter().enumerate() {
        let units = u64::try_from(index).unwrap_or(u64::MAX).saturating_add(1);
        work.checkpoint(phase, plan_index, units)
            .map_err(|_| EvidenceCompareError::Cancelled)?;
        let matches = edge.role == role && edge.artifact == artifact;
        work.checkpoint(phase, plan_index, units.saturating_add(1))
            .map_err(|_| EvidenceCompareError::Cancelled)?;
        work.charge(EDGE_ITEM_WORK_BYTES)
            .map_err(EvidenceCompareError::WorkEnvelope)?;
        if matches {
            work.checkpoint(phase, plan_index, u64::MAX)
                .map_err(|_| EvidenceCompareError::Cancelled)?;
            return Ok(true);
        }
    }
    work.checkpoint(phase, plan_index, u64::MAX)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    Ok(false)
}

fn edge_sets_match_controlled(
    edges: &[OpArtifactEdge],
    expected: &[(EdgeRole, ContentHash)],
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<bool, EvidenceCompareError> {
    let phase = SolveEvidencePhase::EdgeSetCompare;
    let plan_index = Some(candidate_index);
    work.checkpoint(phase, plan_index, 0)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    if edges.len() != expected.len() {
        work.checkpoint(phase, plan_index, 1)
            .map_err(|_| EvidenceCompareError::Cancelled)?;
        work.charge(1).map_err(EvidenceCompareError::WorkEnvelope)?;
        return Ok(false);
    }
    let mut comparisons = 0u64;
    for (role, artifact) in expected {
        let mut found = false;
        for edge in edges {
            work.checkpoint(phase, plan_index, comparisons.saturating_add(1))
                .map_err(|_| EvidenceCompareError::Cancelled)?;
            let matches = edge.role == *role && edge.artifact == *artifact;
            comparisons = comparisons.saturating_add(1);
            work.checkpoint(phase, plan_index, comparisons.saturating_add(1))
                .map_err(|_| EvidenceCompareError::Cancelled)?;
            work.charge(EDGE_ITEM_WORK_BYTES)
                .map_err(EvidenceCompareError::WorkEnvelope)?;
            if matches {
                found = true;
                break;
            }
        }
        if !found {
            return Ok(false);
        }
    }
    for edge in edges {
        let mut found = false;
        for (role, artifact) in expected {
            work.checkpoint(phase, plan_index, comparisons.saturating_add(1))
                .map_err(|_| EvidenceCompareError::Cancelled)?;
            let matches = edge.role == *role && edge.artifact == *artifact;
            comparisons = comparisons.saturating_add(1);
            work.checkpoint(phase, plan_index, comparisons.saturating_add(1))
                .map_err(|_| EvidenceCompareError::Cancelled)?;
            work.charge(EDGE_ITEM_WORK_BYTES)
                .map_err(EvidenceCompareError::WorkEnvelope)?;
            if matches {
                found = true;
                break;
            }
        }
        if !found {
            return Ok(false);
        }
    }
    work.checkpoint(phase, plan_index, u64::MAX)
        .map_err(|_| EvidenceCompareError::Cancelled)?;
    Ok(true)
}

fn require_exact_edges(
    run: SolveRunId,
    stage: SolveStage,
    op: i64,
    edges: &[OpArtifactEdge],
    expected: &[(EdgeRole, ContentHash)],
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<(), SolveRefusal> {
    let matches =
        edge_sets_match_controlled(edges, expected, work, candidate_index).map_err(|error| {
            match error {
                EvidenceCompareError::Cancelled => cancelled_resume_refusal(run),
                EvidenceCompareError::WorkEnvelope(error) => {
                    invocation_work_refusal(Some(run), Some(stage), error)
                }
            }
        })?;
    if !matches {
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
    locator: ImportCandidateLocator,
    work: EvidenceWork<'_>,
    candidate_index: usize,
    expectations: RenderedExplicitsRef<'_>,
) -> Result<ImportSummary, ImportSummaryError> {
    let ImportCandidateLocator {
        op,
        summary_artifact,
    } = locator;
    let candidate = read_candidate_op(ledger, op, work, candidate_index)?.ok_or_else(|| {
        ImportSummaryError::Invalid(format!("import summary names missing operation {op}"))
    })?;
    let attestation = attest_import_row(
        ledger,
        spec,
        op,
        &candidate,
        work,
        candidate_index,
        expectations,
    )?;
    let edges = read_op_edges_controlled(ledger, op, work, candidate_index)?;
    if edges.truncated {
        return Err(ImportSummaryError::Invalid(format!(
            "import operation {op} exceeds the complete {EDGE_SCAN_CAP}-edge scan"
        )));
    }
    validate_import_evidence(
        ledger,
        spec,
        project_hash,
        op,
        summary_artifact,
        &attestation,
        &edges.edges,
        work,
        candidate_index,
    )
}

fn attest_import_row(
    ledger: &Ledger,
    spec: &ProjectSpec,
    op: i64,
    candidate: &ControlledCandidateOp,
    work: EvidenceWork<'_>,
    candidate_index: usize,
    expectations: RenderedExplicitsRef<'_>,
) -> Result<ImportIrAttestation, ImportSummaryError> {
    let row = &candidate.row;
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
    let attestation = validate_import_ir(&row.ir, spec, expectations.canonical_project_json, work)?;
    if attestation.sources.len() > SOLVE_MAX_IMPORT_SOURCES {
        return Err(ImportSummaryError::Unsupported(format!(
            "import operation {op} carries {} sources above the solve evidence cap {SOLVE_MAX_IMPORT_SOURCES}",
            attestation.sources.len()
        )));
    }
    let explicits_match = five_explicits_match(row, expectations, work, candidate_index).map_err(
        |error| match error {
            EvidenceCompareError::Cancelled => ImportSummaryError::Cancelled,
            EvidenceCompareError::WorkEnvelope(error) => ImportSummaryError::WorkEnvelope(error),
        },
    )?;
    if !explicits_match {
        return Err(ImportSummaryError::Invalid(format!(
            "import operation {op} does not carry the project's exact Five Explicits"
        )));
    }
    if candidate.branch != MAIN_BRANCH || candidate.exec_mode != ExecMode::Deterministic {
        return Err(ImportSummaryError::Invalid(format!(
            "import operation {op} is not deterministic on the main branch"
        )));
    }
    verify_candidate_op_identity(ledger, candidate, work, candidate_index)?;
    Ok(attestation)
}

fn validate_import_evidence(
    ledger: &Ledger,
    spec: &ProjectSpec,
    project_hash: ContentHash,
    op: i64,
    summary_artifact: ContentHash,
    attestation: &ImportIrAttestation,
    edges: &[OpArtifactEdge],
    work: EvidenceWork<'_>,
    candidate_index: usize,
) -> Result<ImportSummary, ImportSummaryError> {
    let has_summary = has_edge_controlled(
        edges,
        EdgeRole::Out,
        summary_artifact,
        work,
        candidate_index,
    )
    .map_err(|error| match error {
        EvidenceCompareError::Cancelled => ImportSummaryError::Cancelled,
        EvidenceCompareError::WorkEnvelope(error) => ImportSummaryError::WorkEnvelope(error),
    })?;
    if !has_summary {
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
        work,
        0,
    )?;
    let summary_len = import_artifact_len(ledger, summary_artifact, "import summary", work, 0)?;
    if summary_len > MAX_RECEIPT_READ_BYTES {
        return Err(ImportSummaryError::Unsupported(format!(
            "import summary is {summary_len} bytes above the solve receipt envelope {MAX_RECEIPT_READ_BYTES}"
        )));
    }
    let bytes = materialize_evidence_artifact(
        ledger,
        work,
        summary_artifact,
        MAX_RECEIPT_READ_BYTES,
        SolveEvidencePhase::ImportSummaryRead,
        None,
    )?
    .ok_or_else(|| {
        ImportSummaryError::Invalid(format!(
            "import summary {} is missing",
            summary_artifact.to_hex()
        ))
    })?;
    let text = evidence_utf8_string(
        &bytes,
        work,
        SolveEvidencePhase::ImportSummaryUtf8,
        None,
        "import summary",
    )
    .map_err(|error| match error {
        EvidenceUtf8Error::Cancelled => ImportSummaryError::Cancelled,
        EvidenceUtf8Error::WorkEnvelope(error) => ImportSummaryError::WorkEnvelope(error),
        EvidenceUtf8Error::Invalid(problem) => ImportSummaryError::Invalid(format!(
            "import summary {} is invalid: {problem}",
            summary_artifact.to_hex()
        )),
    })?;
    let entries = parse_geometry_import_summary(&text, spec, project_hash, work)?;
    let mut expected_edges = vec![(EdgeRole::Out, summary_artifact)];
    for (source_index, entry) in entries.iter().enumerate() {
        let descriptor_base = source_index.saturating_mul(4).saturating_add(1);
        require_import_artifact_kind(
            ledger,
            entry.raw_source,
            IMPORT_RAW_KIND,
            "raw source",
            work,
            descriptor_base,
        )?;
        require_import_artifact_kind(
            ledger,
            entry.promotion_receipt,
            IMPORT_PROMOTION_KIND,
            "promotion receipt",
            work,
            descriptor_base.saturating_add(1),
        )?;
        require_import_artifact_kind(
            ledger,
            entry.promoted_mesh,
            IMPORT_MESH_KIND,
            "promoted mesh",
            work,
            descriptor_base.saturating_add(2),
        )?;
        require_import_artifact_kind(
            ledger,
            entry.assignment_report,
            IMPORT_ASSIGNMENT_KIND,
            "assignment report",
            work,
            descriptor_base.saturating_add(3),
        )?;
        expected_edges.push((EdgeRole::In, entry.raw_source));
        expected_edges.push((EdgeRole::Out, entry.promotion_receipt));
        expected_edges.push((EdgeRole::Out, entry.promoted_mesh));
        expected_edges.push((EdgeRole::Out, entry.assignment_report));
    }
    let edges_match = edge_sets_match_controlled(edges, &expected_edges, work, candidate_index)
        .map_err(|error| match error {
            EvidenceCompareError::Cancelled => ImportSummaryError::Cancelled,
            EvidenceCompareError::WorkEnvelope(error) => ImportSummaryError::WorkEnvelope(error),
        })?;
    if !edges_match {
        return Err(ImportSummaryError::Invalid(format!(
            "import operation {op} has {} artifact edges, not the exact typed {}-edge set",
            edges.len(),
            expected_edges.len()
        )));
    }
    validate_import_admission_evidence(
        ledger,
        spec,
        summary_artifact,
        &entries,
        attestation,
        work,
    )?;
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
    work: EvidenceWork<'_>,
    descriptor_index: usize,
) -> Result<(), ImportSummaryError> {
    let info = read_artifact_info_controlled(ledger, &artifact, work, descriptor_index)?
        .ok_or_else(|| {
            ImportSummaryError::Invalid(format!(
                "{label} artifact {} is missing",
                artifact.to_hex()
            ))
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
    work: EvidenceWork<'_>,
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
    work.checkpoint(SolveEvidencePhase::EntityResolution, None, 0)
        .map_err(|_| ImportSummaryError::Cancelled)?;
    let entity_ids = spec.resolve_entities(&mut entity_violations);
    work.checkpoint(SolveEvidencePhase::EntityResolution, None, 1)
        .map_err(|_| ImportSummaryError::Cancelled)?;
    let entity_items = entity_ids
        .len()
        .checked_add(entity_violations.len())
        .and_then(|count| u64::try_from(count).ok())
        .and_then(|count| count.checked_mul(DERIVATION_ITEM_WORK_BYTES))
        .ok_or(ImportSummaryError::WorkEnvelope(
            InvocationWorkExceeded::CumulativeBytes {
                attempted: u64::MAX,
            },
        ))?;
    work.charge(entity_items)
        .map_err(ImportSummaryError::WorkEnvelope)?;
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
        import_artifact_len(ledger, summary_artifact, "import summary", work, 0)?;
    for (source_index, entry) in entries.iter().enumerate() {
        let descriptor_base = source_index.saturating_mul(4).saturating_add(1);
        let raw_len = import_artifact_len(
            ledger,
            entry.raw_source,
            "raw source",
            work,
            descriptor_base,
        )?;
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
        let promotion_len = import_artifact_len(
            ledger,
            entry.promotion_receipt,
            "promotion receipt",
            work,
            descriptor_base.saturating_add(1),
        )?;
        if promotion_len > MAX_OPAQUE_IMPORT_RECEIPT_BYTES {
            return Err(ImportSummaryError::Unsupported(format!(
                "source {source_index} promotion receipt is {promotion_len} bytes above the solve opaque-receipt envelope {MAX_OPAQUE_IMPORT_RECEIPT_BYTES}"
            )));
        }
        let mesh_len = import_artifact_len(
            ledger,
            entry.promoted_mesh,
            "promoted mesh",
            work,
            descriptor_base.saturating_add(2),
        )?;
        if mesh_len > MAX_PARSED_EVIDENCE_BYTES {
            return Err(ImportSummaryError::Unsupported(format!(
                "source {source_index} promoted mesh is {mesh_len} bytes above the solve parse envelope {MAX_PARSED_EVIDENCE_BYTES}"
            )));
        }
        let report_len = import_artifact_len(
            ledger,
            entry.assignment_report,
            "assignment report",
            work,
            descriptor_base.saturating_add(3),
        )?;
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
        let effective_label_max = attestation
            .limits
            .max_label_bytes
            .min(MAX_SOLVE_EVIDENCE_LABEL_BYTES);
        if entry.source_identity.len() > effective_label_max {
            return Err(ImportSummaryError::Invalid(format!(
                "source {source_index} identity exceeds the effective solve label ceiling {effective_label_max}"
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

        let raw_len = verify_evidence_artifact(
            ledger,
            work,
            entry.raw_source,
            raw_stream_cap,
            SolveEvidencePhase::RawSourceRead,
            Some(source_index),
        )
            .map_err(|error| match error {
                EvidenceReadError::Cancelled => ImportSummaryError::Cancelled,
                EvidenceReadError::WorkEnvelope(error) => {
                    ImportSummaryError::WorkEnvelope(error)
                }
                EvidenceReadError::Ledger(LedgerError::ArtifactReadLimit {
                    observed,
                    ..
                }) if observed > source_cap => {
                    ImportSummaryError::Invalid(format!(
                        "source {source_index} raw artifact exceeds max_source_bytes {}",
                        attestation.limits.max_source_bytes
                    ))
                }
                EvidenceReadError::Ledger(LedgerError::ArtifactReadLimit {
                    observed,
                    ..
                }) => {
                    ImportSummaryError::Unsupported(format!(
                        "source {source_index} raw artifact is {observed} bytes above the effective solve stream envelope {raw_stream_cap}"
                    ))
                }
                EvidenceReadError::Ledger(error) => ImportSummaryError::Ledger(error),
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
        verify_evidence_artifact(
            ledger,
            work,
            entry.promotion_receipt,
            MAX_OPAQUE_IMPORT_RECEIPT_BYTES,
            SolveEvidencePhase::PromotionReceiptRead,
            Some(source_index),
        )?
        .ok_or_else(|| {
            ImportSummaryError::Invalid(format!(
                "source {source_index} promotion receipt {} is missing",
                entry.promotion_receipt.to_hex()
            ))
        })?;

        let mesh_bytes = read_parsed_import_artifact(
            ledger,
            work,
            entry.promoted_mesh,
            "promoted mesh",
            source_index,
            SolveEvidencePhase::PromotedMeshRead,
        )?;
        let ply_shape = match preflight_canonical_ply(
            &mesh_bytes,
            attestation.limits.max_mesh_vertices,
            attestation.limits.max_mesh_faces,
            work,
            source_index,
        ) {
            Ok(shape) => shape,
            Err(ImportSummaryError::Cancelled) => return Err(ImportSummaryError::Cancelled),
            Err(ImportSummaryError::Invalid(what)) => {
                return Err(ImportSummaryError::Invalid(format!(
                    "source {source_index} promoted PLY is not canonical: {what}"
                )));
            }
            Err(error) => return Err(error),
        };
        validate_canonical_ply_payload(&mesh_bytes, ply_shape, work, source_index).map_err(
            |error| match error {
                ImportSummaryError::Cancelled => ImportSummaryError::Cancelled,
                ImportSummaryError::Invalid(what) => ImportSummaryError::Invalid(format!(
                    "source {source_index} promoted PLY does not decode canonically: {what}"
                )),
                error => error,
            },
        )?;
        verify_canonical_ply_writer_spelling(&mesh_bytes, ply_shape, work, source_index).map_err(
            |error| match error {
                ImportSummaryError::Cancelled => ImportSummaryError::Cancelled,
                ImportSummaryError::Invalid(what) => ImportSummaryError::Invalid(format!(
                    "source {source_index} promoted PLY is not exact writer output: {what}"
                )),
                error => error,
            },
        )?;
        validate_named_group_face_ranges(
            &source.named_groups,
            ply_shape.faces,
            work,
            source_index,
        )?;

        let assignment_phase = SolveEvidencePhase::AssignmentDerivation;
        let assignment_plan_index = Some(source_index);
        work.checkpoint(assignment_phase, assignment_plan_index, 0)
            .map_err(|_| ImportSummaryError::Cancelled)?;
        let mut assignment_units = 0u64;
        let mut rows = Vec::new();
        let reserve = rows.try_reserve_exact(assignments.len());
        work.checkpoint(assignment_phase, assignment_plan_index, 1)
            .map_err(|_| ImportSummaryError::Cancelled)?;
        reserve.map_err(|_| {
            ImportSummaryError::Invalid("assignment-row filter allocation was refused".to_string())
        })?;
        for assignment in assignments {
            work.checkpoint(
                assignment_phase,
                assignment_plan_index,
                assignment_units.saturating_add(2),
            )
            .map_err(|_| ImportSummaryError::Cancelled)?;
            let selected = assignment.artifact == artifact.role;
            assignment_units = assignment_units.saturating_add(1);
            work.checkpoint(
                assignment_phase,
                assignment_plan_index,
                assignment_units.saturating_add(2),
            )
            .map_err(|_| ImportSummaryError::Cancelled)?;
            work.charge(DERIVATION_ITEM_WORK_BYTES)
                .map_err(ImportSummaryError::WorkEnvelope)?;
            if selected {
                rows.push(assignment);
            }
        }
        let mut expected_subjects = Vec::new();
        let reserve = expected_subjects.try_reserve_exact(rows.len());
        work.checkpoint(
            assignment_phase,
            assignment_plan_index,
            assignment_units.saturating_add(2),
        )
        .map_err(|_| ImportSummaryError::Cancelled)?;
        reserve.map_err(|_| {
            ImportSummaryError::Invalid(
                "assignment subject-attestation allocation refused".to_string(),
            )
        })?;
        let mut geometric_requests = 0u64;
        for row in &rows {
            work.checkpoint(
                assignment_phase,
                assignment_plan_index,
                assignment_units.saturating_add(2),
            )
            .map_err(|_| ImportSummaryError::Cancelled)?;
            let derived = (|| -> Result<(String, bool), ImportSummaryError> {
                let entity = entity_ids.get(&row.target).ok_or_else(|| {
                    ImportSummaryError::Invalid(format!(
                        "assignment target `{}` has no re-derived entity identity",
                        row.target
                    ))
                })?;
                let subject = entity.token();
                if subject.len() > effective_label_max {
                    return Err(ImportSummaryError::Invalid(format!(
                        "assignment subject `{subject}` exceeds the effective solve label ceiling {effective_label_max}"
                    )));
                }
                if let MeshSelector::NamedGroup { name } = &row.selector {
                    if name.len() > effective_label_max {
                        return Err(ImportSummaryError::Invalid(format!(
                            "named-group selector `{name}` exceeds the effective solve label ceiling {effective_label_max}"
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
                let geometric = matches!(
                    &row.selector,
                    MeshSelector::HalfSpace { .. }
                        | MeshSelector::Box { .. }
                        | MeshSelector::Cylinder { .. }
                        | MeshSelector::NearestDatum { .. }
                );
                Ok((subject, geometric))
            })();
            assignment_units = assignment_units.saturating_add(1);
            work.checkpoint(
                assignment_phase,
                assignment_plan_index,
                assignment_units.saturating_add(2),
            )
            .map_err(|_| ImportSummaryError::Cancelled)?;
            let extra_items = match &row.selector {
                MeshSelector::ExplicitFaceSet { faces, .. } => faces.len(),
                _ => 0,
            };
            let charge = u64::try_from(extra_items)
                .ok()
                .and_then(|items| items.checked_mul(4))
                .and_then(|bytes| bytes.checked_add(DERIVATION_ITEM_WORK_BYTES))
                .ok_or(ImportSummaryError::WorkEnvelope(
                    InvocationWorkExceeded::CumulativeBytes {
                        attempted: u64::MAX,
                    },
                ))?;
            work.charge(charge)
                .map_err(ImportSummaryError::WorkEnvelope)?;
            let (subject, geometric) = derived?;
            if geometric {
                geometric_requests = geometric_requests.checked_add(1).ok_or_else(|| {
                    ImportSummaryError::Invalid(
                        "geometric assignment request count overflowed u64".to_string(),
                    )
                })?;
            }
            expected_subjects.push((subject, row.allow_overlap));
        }
        work.checkpoint(assignment_phase, assignment_plan_index, u64::MAX)
            .map_err(|_| ImportSummaryError::Cancelled)?;
        let face_count = u64::try_from(ply_shape.faces).map_err(|_| {
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
            work,
            entry.assignment_report,
            "assignment report",
            source_index,
            SolveEvidencePhase::AssignmentReportRead,
        )?;
        let report_text = evidence_utf8_string(
            &report_bytes,
            work,
            SolveEvidencePhase::AssignmentReportUtf8,
            Some(source_index),
            "assignment report",
        )
        .map_err(|error| match error {
            EvidenceUtf8Error::Cancelled => ImportSummaryError::Cancelled,
            EvidenceUtf8Error::WorkEnvelope(error) => ImportSummaryError::WorkEnvelope(error),
            EvidenceUtf8Error::Invalid(problem) => ImportSummaryError::Invalid(format!(
                "source {source_index} assignment report is invalid: {problem}"
            )),
        })?;
        let selected = parse_assignment_report_counts(
            &report_text,
            &entry.source_identity,
            &source.length_unit,
            &expected_subjects,
            ply_shape.faces,
            work,
            source_index,
        )
        .map_err(|error| match error {
            ImportSummaryError::Cancelled => ImportSummaryError::Cancelled,
            ImportSummaryError::Invalid(what) => ImportSummaryError::Invalid(format!(
                "source {source_index} assignment report is not canonical: {what}"
            )),
            error => error,
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
    work: EvidenceWork<'_>,
    descriptor_index: usize,
) -> Result<u64, ImportSummaryError> {
    read_artifact_info_controlled(ledger, &artifact, work, descriptor_index)?
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
    work: EvidenceWork<'_>,
    artifact: ContentHash,
    label: &str,
    source_index: usize,
    phase: SolveEvidencePhase,
) -> Result<Vec<u8>, ImportSummaryError> {
    let info =
        read_artifact_info_controlled(ledger, &artifact, work, source_index)?.ok_or_else(|| {
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
    materialize_evidence_artifact(
        ledger,
        work,
        artifact,
        MAX_PARSED_EVIDENCE_BYTES,
        phase,
        Some(source_index),
    )?
    .ok_or_else(|| {
        ImportSummaryError::Invalid(format!(
            "source {source_index} {label} artifact {} disappeared",
            artifact.to_hex()
        ))
    })
}

#[derive(Debug, Clone, Copy)]
struct CanonicalPlyShape {
    body_start: usize,
    vertices: usize,
    faces: usize,
}

fn preflight_canonical_ply(
    bytes: &[u8],
    max_vertices: usize,
    max_faces: usize,
    work: EvidenceWork<'_>,
    source_index: usize,
) -> Result<CanonicalPlyShape, ImportSummaryError> {
    work.checkpoint(
        SolveEvidencePhase::PromotedMeshPreflight,
        Some(source_index),
        0,
    )
    .map_err(|_| ImportSummaryError::Cancelled)?;
    let result = preflight_canonical_ply_inner(bytes, max_vertices, max_faces, work, source_index);
    let completion_units = u64::try_from(bytes.len()).map_err(|_| {
        ImportSummaryError::Invalid("PLY artifact length is outside u64".to_string())
    })?;
    work.checkpoint(
        SolveEvidencePhase::PromotedMeshPreflight,
        Some(source_index),
        completion_units,
    )
    .map_err(|_| ImportSummaryError::Cancelled)?;
    result
}

fn preflight_canonical_ply_inner(
    bytes: &[u8],
    max_vertices: usize,
    max_faces: usize,
    work: EvidenceWork<'_>,
    source_index: usize,
) -> Result<CanonicalPlyShape, ImportSummaryError> {
    const END_HEADER: &[u8] = b"end_header\n";
    const MAX_CANONICAL_HEADER_BYTES: usize = 256;
    let header_prefix = &bytes[..bytes.len().min(MAX_CANONICAL_HEADER_BYTES)];
    let header_end = header_prefix
        .windows(END_HEADER.len())
        .position(|window| window == END_HEADER)
        .and_then(|offset| offset.checked_add(END_HEADER.len()))
        .ok_or_else(|| {
            ImportSummaryError::Invalid(
                "the exact bounded `end_header\\n` terminator is missing".to_string(),
            )
        })?;
    let header = core::str::from_utf8(&bytes[..header_end])
        .map_err(|_| ImportSummaryError::Invalid("the header is not UTF-8".to_string()))?;
    let header = header.strip_suffix('\n').ok_or_else(|| {
        ImportSummaryError::Invalid("the header is not newline-terminated".to_string())
    })?;
    let mut lines = header.split('\n');
    if lines.next() != Some("ply") || lines.next() != Some("format ascii 1.0") {
        return Err(ImportSummaryError::Invalid(
            "the magic or ASCII format line differs from the fs-io writer".to_string(),
        ));
    }
    let vertices = parse_canonical_ply_count(
        lines.next().ok_or_else(|| {
            ImportSummaryError::Invalid("the vertex element line is missing".to_string())
        })?,
        "element vertex ",
        "vertex",
    )
    .map_err(ImportSummaryError::Invalid)?;
    if lines.next() != Some("property double x")
        || lines.next() != Some("property double y")
        || lines.next() != Some("property double z")
    {
        return Err(ImportSummaryError::Invalid(
            "the vertex property lines differ from the fs-io writer".to_string(),
        ));
    }
    let faces = parse_canonical_ply_count(
        lines.next().ok_or_else(|| {
            ImportSummaryError::Invalid("the face element line is missing".to_string())
        })?,
        "element face ",
        "face",
    )
    .map_err(ImportSummaryError::Invalid)?;
    if lines.next() != Some("property list uchar uint vertex_indices")
        || lines.next() != Some("end_header")
        || lines.next().is_some()
    {
        return Err(ImportSummaryError::Invalid(
            "the face property or header shape differs from the fs-io writer".to_string(),
        ));
    }
    if vertices > max_vertices || faces > max_faces {
        return Err(ImportSummaryError::Invalid(format!(
            "declared counts {vertices} vertices and {faces} faces exceed frozen caps {max_vertices} and {max_faces}"
        )));
    }
    if vertices == 0 || faces == 0 {
        return Err(ImportSummaryError::Invalid(
            "the canonical promoted mesh must contain at least one vertex and one face".to_string(),
        ));
    }

    let expected_records = vertices.checked_add(faces).ok_or_else(|| {
        ImportSummaryError::Invalid("the declared body-record count overflows usize".to_string())
    })?;
    let body = &bytes[header_end..];
    if !body.is_empty() && body.last() != Some(&b'\n') {
        return Err(ImportSummaryError::Invalid(
            "the body is not newline-terminated".to_string(),
        ));
    }
    let mut body_records = 0usize;
    let mut inspected = 0u64;
    for tile in body.chunks(EVIDENCE_POLL_BYTES) {
        work.checkpoint(
            SolveEvidencePhase::PromotedMeshPreflight,
            Some(source_index),
            inspected,
        )
        .map_err(|_| ImportSummaryError::Cancelled)?;
        body_records = body_records
            .checked_add(tile.iter().filter(|byte| **byte == b'\n').count())
            .ok_or_else(|| {
                ImportSummaryError::Invalid("PLY body-record count overflowed usize".to_string())
            })?;
        inspected = inspected
            .checked_add(u64::try_from(tile.len()).map_err(|_| {
                ImportSummaryError::Invalid("PLY body length is outside u64".to_string())
            })?)
            .ok_or_else(|| {
                ImportSummaryError::Invalid("PLY body byte count overflowed u64".to_string())
            })?;
        work.checkpoint(
            SolveEvidencePhase::PromotedMeshPreflight,
            Some(source_index),
            inspected,
        )
        .map_err(|_| ImportSummaryError::Cancelled)?;
    }
    if body_records != expected_records {
        return Err(ImportSummaryError::Invalid(format!(
            "the body has {body_records} newline-terminated records but the header declares {expected_records}"
        )));
    }
    Ok(CanonicalPlyShape {
        body_start: header_end,
        vertices,
        faces,
    })
}

fn parse_canonical_ply_count(line: &str, prefix: &str, label: &str) -> Result<usize, String> {
    let spelling = line
        .strip_prefix(prefix)
        .ok_or_else(|| format!("the {label} element line differs from the fs-io writer"))?;
    if spelling.is_empty()
        || spelling.len() > MAX_JSON_INTEGER_BYTES
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

fn validate_canonical_ply_payload(
    bytes: &[u8],
    shape: CanonicalPlyShape,
    work: EvidenceWork<'_>,
    source_index: usize,
) -> Result<(), ImportSummaryError> {
    work.checkpoint(
        SolveEvidencePhase::PromotedMeshDecode,
        Some(source_index),
        0,
    )
    .map_err(|_| ImportSummaryError::Cancelled)?;
    let result = walk_canonical_ply_payload(
        bytes,
        shape,
        work,
        source_index,
        SolveEvidencePhase::PromotedMeshDecode,
        false,
    );
    let completion_units = u64::try_from(bytes.len()).map_err(|_| {
        ImportSummaryError::Invalid("PLY artifact length is outside u64".to_string())
    })?;
    work.checkpoint(
        SolveEvidencePhase::PromotedMeshDecode,
        Some(source_index),
        completion_units,
    )
    .map_err(|_| ImportSummaryError::Cancelled)?;
    result
}

fn verify_canonical_ply_writer_spelling(
    bytes: &[u8],
    shape: CanonicalPlyShape,
    work: EvidenceWork<'_>,
    source_index: usize,
) -> Result<(), ImportSummaryError> {
    work.checkpoint(
        SolveEvidencePhase::PromotedMeshEncodeCompare,
        Some(source_index),
        0,
    )
    .map_err(|_| ImportSummaryError::Cancelled)?;
    let result = walk_canonical_ply_payload(
        bytes,
        shape,
        work,
        source_index,
        SolveEvidencePhase::PromotedMeshEncodeCompare,
        true,
    );
    let completion_units = u64::try_from(bytes.len()).map_err(|_| {
        ImportSummaryError::Invalid("PLY artifact length is outside u64".to_string())
    })?;
    work.checkpoint(
        SolveEvidencePhase::PromotedMeshEncodeCompare,
        Some(source_index),
        completion_units,
    )
    .map_err(|_| ImportSummaryError::Cancelled)?;
    result
}

fn walk_canonical_ply_payload(
    bytes: &[u8],
    shape: CanonicalPlyShape,
    work: EvidenceWork<'_>,
    source_index: usize,
    phase: SolveEvidencePhase,
    require_writer_spelling: bool,
) -> Result<(), ImportSummaryError> {
    work.checkpoint(phase, Some(source_index), 0)
        .map_err(|_| ImportSummaryError::Cancelled)?;
    let body = bytes.get(shape.body_start..).ok_or_else(|| {
        ImportSummaryError::Invalid(
            "the canonical PLY body offset is outside the retained artifact".to_string(),
        )
    })?;
    let mut cursor = 0usize;
    for record in 0..shape.vertices {
        let line = next_bounded_ply_line(
            body,
            &mut cursor,
            MAX_CANONICAL_PLY_VERTEX_LINE_BYTES,
            work,
            source_index,
            phase,
            "vertex",
            record,
        );
        let current = u64::try_from(cursor).map_err(|_| {
            ImportSummaryError::Invalid("PLY body cursor is outside u64".to_string())
        })?;
        work.checkpoint(phase, Some(source_index), current)
            .map_err(|_| ImportSummaryError::Cancelled)?;
        let line = line?;
        let parsed = parse_canonical_ply_vertex_line(line, require_writer_spelling);
        work.checkpoint(phase, Some(source_index), current)
            .map_err(|_| ImportSummaryError::Cancelled)?;
        parsed.map_err(|problem| {
            ImportSummaryError::Invalid(format!("vertex record {record}: {problem}"))
        })?;
    }
    for record in 0..shape.faces {
        let line = next_bounded_ply_line(
            body,
            &mut cursor,
            MAX_CANONICAL_PLY_FACE_LINE_BYTES,
            work,
            source_index,
            phase,
            "face",
            record,
        );
        let current = u64::try_from(cursor).map_err(|_| {
            ImportSummaryError::Invalid("PLY body cursor is outside u64".to_string())
        })?;
        work.checkpoint(phase, Some(source_index), current)
            .map_err(|_| ImportSummaryError::Cancelled)?;
        let line = line?;
        let parsed = parse_canonical_ply_face_line(line, shape.vertices, require_writer_spelling);
        work.checkpoint(phase, Some(source_index), current)
            .map_err(|_| ImportSummaryError::Cancelled)?;
        parsed.map_err(|problem| {
            ImportSummaryError::Invalid(format!("face record {record}: {problem}"))
        })?;
    }
    let completion = u64::try_from(cursor)
        .map_err(|_| ImportSummaryError::Invalid("PLY body cursor is outside u64".to_string()))?;
    work.checkpoint(phase, Some(source_index), completion)
        .map_err(|_| ImportSummaryError::Cancelled)?;
    if cursor != body.len() {
        return Err(ImportSummaryError::Invalid(format!(
            "canonical PLY parser consumed {cursor} of {} body bytes",
            body.len()
        )));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn next_bounded_ply_line<'a>(
    body: &'a [u8],
    cursor: &mut usize,
    max_line_bytes: usize,
    work: EvidenceWork<'_>,
    source_index: usize,
    phase: SolveEvidencePhase,
    record_kind: &str,
    record: usize,
) -> Result<&'a str, ImportSummaryError> {
    debug_assert!(max_line_bytes < EVIDENCE_POLL_BYTES);
    let before = u64::try_from(*cursor).map_err(|_| {
        ImportSummaryError::Invalid("PLY body cursor is outside the cancellation range".to_string())
    })?;
    work.checkpoint(phase, Some(source_index), before)
        .map_err(|_| ImportSummaryError::Cancelled)?;
    let remaining = body.get(*cursor..).ok_or_else(|| {
        ImportSummaryError::Invalid("PLY body cursor advanced past the artifact".to_string())
    })?;
    let search_bytes = max_line_bytes.checked_add(1).ok_or_else(|| {
        ImportSummaryError::Invalid("PLY line-search bound overflowed usize".to_string())
    })?;
    let inspected = &remaining[..remaining.len().min(search_bytes)];
    let newline = inspected
        .iter()
        .position(|byte| *byte == b'\n')
        .ok_or_else(|| {
            if remaining.len() > max_line_bytes {
                ImportSummaryError::Invalid(format!(
                    "{record_kind} record {record} exceeds the {max_line_bytes}-byte canonical line bound"
                ))
            } else {
                ImportSummaryError::Invalid(format!(
                    "{record_kind} record {record} is not newline-terminated"
                ))
            }
        })?;
    if newline > max_line_bytes {
        return Err(ImportSummaryError::Invalid(format!(
            "{record_kind} record {record} exceeds the {max_line_bytes}-byte canonical line bound"
        )));
    }
    let line = core::str::from_utf8(&remaining[..newline]).map_err(|_| {
        ImportSummaryError::Invalid(format!(
            "{record_kind} record {record} is not ASCII-compatible UTF-8"
        ))
    })?;
    *cursor = cursor
        .checked_add(newline)
        .and_then(|offset| offset.checked_add(1))
        .ok_or_else(|| {
            ImportSummaryError::Invalid("PLY body cursor overflowed usize".to_string())
        })?;
    let after = u64::try_from(*cursor).map_err(|_| {
        ImportSummaryError::Invalid("PLY body cursor is outside the cancellation range".to_string())
    })?;
    work.checkpoint(phase, Some(source_index), after)
        .map_err(|_| ImportSummaryError::Cancelled)?;
    Ok(line)
}

fn parse_canonical_ply_vertex_line(
    line: &str,
    require_writer_spelling: bool,
) -> Result<(), String> {
    let mut fields = line.split(' ');
    let x = fields
        .next()
        .ok_or_else(|| "x coordinate is missing".to_string())?;
    let y = fields
        .next()
        .ok_or_else(|| "y coordinate is missing".to_string())?;
    let z = fields
        .next()
        .ok_or_else(|| "z coordinate is missing".to_string())?;
    if fields.next().is_some() {
        return Err("vertex line does not contain exactly three single-space fields".to_string());
    }
    for (axis, spelling) in [("x", x), ("y", y), ("z", z)] {
        if spelling.is_empty() || spelling.len() > MAX_CANONICAL_F64_BYTES {
            return Err(format!(
                "{axis} coordinate exceeds the {MAX_CANONICAL_F64_BYTES}-byte finite-f64 token bound"
            ));
        }
        let value = spelling
            .parse::<f64>()
            .map_err(|_| format!("{axis} coordinate is not an f64"))?;
        if !value.is_finite() {
            return Err(format!("{axis} coordinate is not finite"));
        }
        if require_writer_spelling && value.to_string() != spelling {
            return Err(format!(
                "{axis} coordinate is not the writer's canonical finite-f64 spelling"
            ));
        }
    }
    Ok(())
}

fn parse_canonical_ply_face_line(
    line: &str,
    vertices: usize,
    require_writer_spelling: bool,
) -> Result<(), String> {
    let mut fields = line.split(' ');
    if fields.next() != Some("3") {
        return Err("face line does not begin with the canonical triangle arity `3`".to_string());
    }
    let a = fields
        .next()
        .ok_or_else(|| "first face index is missing".to_string())?;
    let b = fields
        .next()
        .ok_or_else(|| "second face index is missing".to_string())?;
    let c = fields
        .next()
        .ok_or_else(|| "third face index is missing".to_string())?;
    if fields.next().is_some() {
        return Err("face line does not contain exactly four single-space fields".to_string());
    }
    for (ordinal, spelling) in [("first", a), ("second", b), ("third", c)] {
        if spelling.is_empty()
            || spelling.len() > MAX_CANONICAL_U32_BYTES
            || !spelling.bytes().all(|byte| byte.is_ascii_digit())
        {
            return Err(format!(
                "{ordinal} face index exceeds the {MAX_CANONICAL_U32_BYTES}-byte unsigned token bound"
            ));
        }
        let value = spelling
            .parse::<u32>()
            .map_err(|_| format!("{ordinal} face index is outside u32"))?;
        let value_usize =
            usize::try_from(value).map_err(|_| format!("{ordinal} face index is outside usize"))?;
        if value_usize >= vertices {
            return Err(format!(
                "{ordinal} face index {value} is outside the {vertices}-vertex mesh"
            ));
        }
        if require_writer_spelling && value.to_string() != spelling {
            return Err(format!(
                "{ordinal} face index is not the writer's canonical u32 spelling"
            ));
        }
    }
    Ok(())
}

fn validate_named_group_face_ranges(
    groups: &[NamedFaceGroup],
    mesh_faces: usize,
    work: EvidenceWork<'_>,
    source_index: usize,
) -> Result<(), ImportSummaryError> {
    const FACE_REFERENCES_PER_TILE: usize = EVIDENCE_POLL_BYTES / core::mem::size_of::<u32>();

    let mut inspected_bytes = 0u64;
    work.checkpoint(
        SolveEvidencePhase::NamedGroupFaceRange,
        Some(source_index),
        inspected_bytes,
    )
    .map_err(|_| ImportSummaryError::Cancelled)?;
    for group in groups {
        for tile in group.faces.chunks(FACE_REFERENCES_PER_TILE) {
            work.checkpoint(
                SolveEvidencePhase::NamedGroupFaceRange,
                Some(source_index),
                inspected_bytes,
            )
            .map_err(|_| ImportSummaryError::Cancelled)?;
            let mut range_error = None;
            let mut inspected_faces = 0usize;
            for face in tile {
                inspected_faces += 1;
                let face = match usize::try_from(*face) {
                    Ok(face) => face,
                    Err(_) => {
                        range_error = Some(format!(
                            "source {source_index} named group `{}` carries a face outside usize",
                            group.name
                        ));
                        break;
                    }
                };
                if face >= mesh_faces {
                    range_error = Some(format!(
                        "source {source_index} named group `{}` references face {face} outside the {mesh_faces}-face promoted mesh",
                        group.name
                    ));
                    break;
                }
            }
            let tile_bytes = inspected_faces
                .checked_mul(core::mem::size_of::<u32>())
                .and_then(|bytes| u64::try_from(bytes).ok())
                .ok_or_else(|| {
                    ImportSummaryError::Invalid(
                        "named-group face-range work overflowed u64".to_string(),
                    )
                })?;
            inspected_bytes = inspected_bytes.checked_add(tile_bytes).ok_or_else(|| {
                ImportSummaryError::Invalid(
                    "named-group face-range byte count overflowed u64".to_string(),
                )
            })?;
            work.checkpoint(
                SolveEvidencePhase::NamedGroupFaceRange,
                Some(source_index),
                inspected_bytes,
            )
            .map_err(|_| ImportSummaryError::Cancelled)?;
            if let Some(problem) = range_error {
                return Err(ImportSummaryError::Invalid(problem));
            }
        }
    }
    Ok(())
}

fn parse_assignment_report_counts(
    text: &str,
    expected_source: &str,
    expected_unit: &str,
    expected_assignments: &[(String, bool)],
    mesh_faces: usize,
    work: EvidenceWork<'_>,
    source_index: usize,
) -> Result<usize, ImportSummaryError> {
    let mut cursor = JsonCursor::with_work(
        text,
        work,
        SolveEvidencePhase::AssignmentReportParse,
        Some(source_index),
    );
    let result = parse_assignment_report_counts_cursor(
        &mut cursor,
        expected_source,
        expected_unit,
        expected_assignments,
        mesh_faces,
    );
    let final_poll = cursor.checkpoint_current("after assignment-report parsing");
    if cursor.cancellation_observed() {
        Err(ImportSummaryError::Cancelled)
    } else {
        final_poll.map_err(ImportSummaryError::Invalid)?;
        result.map_err(ImportSummaryError::Invalid)
    }
}

fn parse_assignment_report_counts_cursor(
    cursor: &mut JsonCursor<'_>,
    expected_source: &str,
    expected_unit: &str,
    expected_assignments: &[(String, bool)],
    mesh_faces: usize,
) -> Result<usize, String> {
    cursor.expect("{\"kind\":")?;
    if cursor.parse_string()? != "mesh-assignment-receipt" {
        return Err("kind is not mesh-assignment-receipt".to_string());
    }
    cursor.expect(",\"version\":")?;
    if cursor.parse_string()? != MESH_ASSIGNMENT_SEMANTICS_VERSION {
        return Err("version does not match the fs-io writer".to_string());
    }
    cursor.expect(",\"source_identity\":")?;
    let source_identity = cursor.parse_string()?;
    require_solve_evidence_label(
        &source_identity,
        MAX_SOLVE_EVIDENCE_LABEL_BYTES,
        "assignment-report source identity",
    )?;
    if source_identity != expected_source {
        return Err("source identity differs from the import summary".to_string());
    }
    cursor.expect(",\"length_unit\":")?;
    let length_unit = cursor.parse_string()?;
    require_solve_evidence_label(
        &length_unit,
        MAX_SOLVE_EVIDENCE_LABEL_BYTES,
        "assignment-report length unit",
    )?;
    if length_unit != expected_unit {
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
        let subject = cursor.parse_string()?;
        require_solve_evidence_label(
            &subject,
            MAX_SOLVE_EVIDENCE_LABEL_BYTES,
            "assignment-report subject",
        )?;
        if subject != *expected_subject {
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
        parse_finite_vector3(cursor)?;
        cursor.expect("],\"bounds_max\":[")?;
        parse_finite_vector3(cursor)?;
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

fn validate_import_ir(
    ir: &str,
    spec: &ProjectSpec,
    canonical_project_json: &str,
    work: EvidenceWork<'_>,
) -> Result<ImportIrAttestation, ImportSummaryError> {
    let mut cursor = JsonCursor::with_work(ir, work, SolveEvidencePhase::ImportIrParse, None);
    let result = validate_import_ir_cursor(&mut cursor, spec, canonical_project_json, work);
    let final_poll = cursor.checkpoint_current("after import-IR parsing");
    if let Some(error) = cursor.work_exceeded {
        Err(ImportSummaryError::WorkEnvelope(error))
    } else if cursor.cancellation_observed() {
        Err(ImportSummaryError::Cancelled)
    } else {
        final_poll.map_err(ImportSummaryError::Invalid)?;
        result.map_err(ImportSummaryError::Invalid)
    }
}

fn validate_import_ir_cursor(
    cursor: &mut JsonCursor<'_>,
    spec: &ProjectSpec,
    canonical_project_json: &str,
    work: EvidenceWork<'_>,
) -> Result<ImportIrAttestation, String> {
    cursor.expect("{\"schema\":")?;
    let schema = cursor.parse_string()?;
    if schema != IMPORT_IR_SCHEMA {
        return Err(format!(
            "import IR schema is `{schema}`, not `{IMPORT_IR_SCHEMA}`"
        ));
    }
    cursor.expect(",\"project\":")?;
    let project_json = cursor.take_value()?;
    let canonical_matches = match evidence_bytes_equal(
        project_json.as_bytes(),
        canonical_project_json.as_bytes(),
        work,
        SolveEvidencePhase::ImportIrCanonicalCompare,
        None,
    ) {
        Ok(matches) => matches,
        Err(EvidenceCompareError::Cancelled) => {
            cursor.cancelled = true;
            return Err("solve evidence canonical-project comparison stopped".to_string());
        }
        Err(EvidenceCompareError::WorkEnvelope(error)) => {
            cursor.work_exceeded = Some(error);
            return Err(
                "solve evidence canonical-project comparison exceeded the work envelope"
                    .to_string(),
            );
        }
    };
    if !canonical_matches {
        return Err("import IR does not embed the exact canonical project JSON".to_string());
    }
    let limits = parse_import_ir_limits(cursor)?;
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
        require_solve_evidence_label(
            &source_identity,
            limits.max_label_bytes,
            &format!("import IR source {index} identity"),
        )?;
        let expected_identity = match geometry_source_identity_controlled(artifact, work, index) {
            Ok(identity) => identity,
            Err(EvidenceCompareError::Cancelled) => {
                cursor.cancelled = true;
                return Err("project source-identity derivation stopped".to_string());
            }
            Err(EvidenceCompareError::WorkEnvelope(error)) => {
                cursor.work_exceeded = Some(error);
                return Err(
                    "project source-identity derivation exceeded the work envelope".to_string(),
                );
            }
        };
        if source_identity != expected_identity {
            return Err(format!(
                "import IR source {index} identity `{source_identity}` does not match `{expected_identity}`"
            ));
        }
        cursor.expect(",\"policy\":")?;
        sources.push(parse_import_ir_policy(
            cursor, spec, artifact, index, limits, work,
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
    work: EvidenceWork<'_>,
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
            let named_groups = parse_import_ir_named_groups(cursor, source_index, limits, work)?;
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
            let named_groups = parse_import_ir_named_groups(cursor, source_index, limits, work)?;
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

fn require_solve_evidence_label(
    value: &str,
    declared_max: usize,
    label: &str,
) -> Result<(), String> {
    let effective_max = declared_max.min(MAX_SOLVE_EVIDENCE_LABEL_BYTES);
    if value.len() > effective_max {
        Err(format!(
            "{label} is {} bytes above the solve label ceiling {effective_max}",
            value.len()
        ))
    } else {
        Ok(())
    }
}

fn validate_import_ir_unit(
    spec: &ProjectSpec,
    artifact: &GeometryArtifact,
    source_index: usize,
    unit: &str,
    max_label_bytes: usize,
) -> Result<(), String> {
    require_solve_evidence_label(
        unit,
        max_label_bytes,
        &format!("import IR source {source_index} length unit"),
    )?;
    if unit.is_empty() || unit.trim() != unit || unit.chars().any(char::is_control) {
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

fn import_ir_duplicate_checkpoint(
    cursor: &mut JsonCursor<'_>,
    work: EvidenceWork<'_>,
    source_index: usize,
    units: u64,
    context: &str,
) -> Result<(), String> {
    if work
        .checkpoint(
            SolveEvidencePhase::ImportIrDuplicateCheck,
            Some(source_index),
            units,
        )
        .is_err()
    {
        cursor.cancelled = true;
        Err(format!("solve evidence parsing stopped {context}"))
    } else {
        Ok(())
    }
}

fn parse_import_ir_named_groups(
    cursor: &mut JsonCursor<'_>,
    source_index: usize,
    limits: ImportIrLimits,
    work: EvidenceWork<'_>,
) -> Result<Vec<NamedFaceGroup>, String> {
    cursor.expect("[")?;
    let mut groups = 0usize;
    let mut face_references = 0usize;
    let mut duplicate_units = 0u64;
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
            require_solve_evidence_label(
                &name,
                limits.max_label_bytes,
                &format!("import IR source {source_index} named-group label"),
            )?;
            if name.is_empty() || name.trim() != name || name.chars().any(char::is_control) {
                return Err(format!(
                    "import IR source {source_index} named-group label is not nonempty, trim-canonical, control-free, and bounded"
                ));
            }
            cursor.checkpoint_current("before named-group allocation")?;
            let group_reserve = named_groups.try_reserve(1);
            cursor.checkpoint_current("after named-group allocation")?;
            group_reserve.map_err(|_| "import IR named-group allocation refused".to_string())?;
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
                cursor.checkpoint_current("before named-group face allocation")?;
                let face_reserve = group_faces.try_reserve(1);
                cursor.checkpoint_current("after named-group face allocation")?;
                face_reserve
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
            import_ir_duplicate_checkpoint(
                cursor,
                work,
                source_index,
                duplicate_units,
                "before named-group face-set allocation",
            )?;
            let mut seen_faces = SolveEvidenceSet::default();
            let reserve = seen_faces.try_reserve(group_faces.len());
            import_ir_duplicate_checkpoint(
                cursor,
                work,
                source_index,
                duplicate_units,
                "after named-group face-set allocation",
            )?;
            reserve.map_err(|_| {
                "import IR named-group face duplicate-set allocation refused".to_string()
            })?;
            let faces_per_tile = (EVIDENCE_POLL_BYTES / core::mem::size_of::<u32>()).max(1);
            for tile in group_faces.chunks(faces_per_tile) {
                import_ir_duplicate_checkpoint(
                    cursor,
                    work,
                    source_index,
                    duplicate_units,
                    "before named-group face duplicate tile",
                )?;
                let mut duplicate = None;
                let mut inspected = 0usize;
                for face in tile {
                    inspected += 1;
                    if !seen_faces.insert(*face) {
                        duplicate = Some(*face);
                        break;
                    }
                }
                let inspected_bytes = inspected
                    .checked_mul(core::mem::size_of::<u32>())
                    .and_then(|bytes| u64::try_from(bytes).ok())
                    .ok_or_else(|| {
                        "import IR named-group duplicate-work count overflowed u64".to_string()
                    })?;
                duplicate_units =
                    duplicate_units
                        .checked_add(inspected_bytes)
                        .ok_or_else(|| {
                            "import IR named-group duplicate-work count overflowed u64".to_string()
                        })?;
                import_ir_duplicate_checkpoint(
                    cursor,
                    work,
                    source_index,
                    duplicate_units,
                    "after named-group face duplicate tile",
                )?;
                if let Some(duplicate) = duplicate {
                    return Err(format!(
                        "import IR source {source_index} named group repeats face {duplicate}"
                    ));
                }
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
    import_ir_duplicate_checkpoint(
        cursor,
        work,
        source_index,
        duplicate_units,
        "before named-group name-set allocation",
    )?;
    let mut seen_names = SolveEvidenceSet::default();
    let reserve = seen_names.try_reserve(named_groups.len());
    import_ir_duplicate_checkpoint(
        cursor,
        work,
        source_index,
        duplicate_units,
        "after named-group name-set allocation",
    )?;
    reserve
        .map_err(|_| "import IR named-group name duplicate-set allocation refused".to_string())?;
    for group in &named_groups {
        import_ir_duplicate_checkpoint(
            cursor,
            work,
            source_index,
            duplicate_units,
            "before named-group name duplicate item",
        )?;
        let mut key = String::new();
        let key_reserve = key.try_reserve_exact(group.name.len());
        import_ir_duplicate_checkpoint(
            cursor,
            work,
            source_index,
            duplicate_units,
            "after named-group name-key allocation",
        )?;
        key_reserve
            .map_err(|_| "import IR named-group duplicate-key allocation refused".to_string())?;
        key.push_str(&group.name);
        let inserted = seen_names.insert(key);
        let name_bytes = u64::try_from(group.name.len())
            .map_err(|_| "import IR named-group name length is outside u64".to_string())?;
        duplicate_units = duplicate_units.checked_add(name_bytes).ok_or_else(|| {
            "import IR named-group duplicate-work count overflowed u64".to_string()
        })?;
        import_ir_duplicate_checkpoint(
            cursor,
            work,
            source_index,
            duplicate_units,
            "after named-group name duplicate item",
        )?;
        if !inserted {
            return Err(format!(
                "import IR source {source_index} repeats named group `{}`",
                group.name
            ));
        }
    }
    Ok(named_groups)
}

fn parse_geometry_import_summary(
    text: &str,
    spec: &ProjectSpec,
    project_hash: ContentHash,
    work: EvidenceWork<'_>,
) -> Result<Vec<VerifiedImport>, ImportSummaryError> {
    let mut cursor =
        JsonCursor::with_work(text, work, SolveEvidencePhase::ImportSummaryParse, None);
    let result = parse_geometry_import_summary_cursor(&mut cursor, spec, project_hash, work);
    let final_poll = cursor.checkpoint_current("after import-summary parsing");
    if let Some(error) = cursor.work_exceeded {
        Err(ImportSummaryError::WorkEnvelope(error))
    } else if cursor.cancellation_observed() {
        Err(ImportSummaryError::Cancelled)
    } else {
        final_poll.map_err(ImportSummaryError::Invalid)?;
        result.map_err(ImportSummaryError::Invalid)
    }
}

fn parse_geometry_import_summary_cursor(
    cursor: &mut JsonCursor<'_>,
    spec: &ProjectSpec,
    project_hash: ContentHash,
    work: EvidenceWork<'_>,
) -> Result<Vec<VerifiedImport>, String> {
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
            entries.push(parse_geometry_import_entry(cursor, project_hash)?);
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
        let expected_identity = match geometry_source_identity_controlled(artifact, work, index) {
            Ok(identity) => identity,
            Err(EvidenceCompareError::Cancelled) => {
                cursor.cancelled = true;
                return Err("project source-identity derivation stopped".to_string());
            }
            Err(EvidenceCompareError::WorkEnvelope(error)) => {
                cursor.work_exceeded = Some(error);
                return Err(
                    "project source-identity derivation exceeded the work envelope".to_string(),
                );
            }
        };
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
    require_solve_evidence_label(
        &role,
        MAX_SOLVE_EVIDENCE_LABEL_BYTES,
        "summary artifact role",
    )?;
    cursor.expect(",\"source_label\":")?;
    let source_label = cursor.parse_string()?;
    require_solve_evidence_label(
        &source_label,
        MAX_SOLVE_EVIDENCE_LABEL_BYTES,
        "summary source label",
    )?;
    if source_label.is_empty() || source_label.chars().any(char::is_control) {
        return Err("summary source label violates the import writer bound".to_string());
    }
    cursor.expect(",\"source_label_authority\":")?;
    let source_label_authority = cursor.parse_string()?;
    if source_label_authority != IMPORT_SOURCE_LABEL_AUTHORITY {
        return Err("summary source-label authority is not caller-reported".to_string());
    }
    cursor.expect(",\"source_identity\":")?;
    let source_identity = cursor.parse_string()?;
    require_solve_evidence_label(
        &source_identity,
        MAX_SOLVE_EVIDENCE_LABEL_BYTES,
        "summary source identity",
    )?;
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
    let max_import_record_bytes = MAX_SOLVE_EVIDENCE_LABEL_BYTES
        .checked_add(65)
        .expect("solve label ceiling plus hash separator fits usize");
    if import_record.len() > max_import_record_bytes {
        return Err(format!(
            "summary import record is {} bytes above the solve composite-record ceiling {max_import_record_bytes}",
            import_record.len()
        ));
    }
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
    work: EvidenceWork<'_>,
) -> Result<(i64, Vec<VerifiedImport>), ImportSummaryError> {
    let mut cursor = JsonCursor::with_work(
        text,
        work,
        SolveEvidencePhase::ResumeStageReceiptParse,
        None,
    );
    let result = parse_import_verify_receipt_cursor(&mut cursor, run, project_hash);
    let final_poll = cursor.checkpoint_current("after stage-receipt parsing");
    if cursor.cancellation_observed() {
        Err(ImportSummaryError::Cancelled)
    } else {
        final_poll.map_err(ImportSummaryError::Invalid)?;
        result.map_err(ImportSummaryError::Invalid)
    }
}

fn parse_import_verify_receipt_cursor(
    cursor: &mut JsonCursor<'_>,
    run: SolveRunId,
    project_hash: ContentHash,
) -> Result<(i64, Vec<VerifiedImport>), String> {
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
            entries.push(parse_verified_import_entry(cursor)?);
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
    require_solve_evidence_label(
        &role,
        MAX_SOLVE_EVIDENCE_LABEL_BYTES,
        "stage-receipt artifact role",
    )?;
    cursor.expect(",\"source_identity\":")?;
    let source_identity = cursor.parse_string()?;
    require_solve_evidence_label(
        &source_identity,
        MAX_SOLVE_EVIDENCE_LABEL_BYTES,
        "stage-receipt source identity",
    )?;
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
    poll: Option<JsonPoll<'a>>,
    cancelled: bool,
    work_exceeded: Option<InvocationWorkExceeded>,
}

struct JsonPoll<'a> {
    work: EvidenceWork<'a>,
    phase: SolveEvidencePhase,
    source_index: Option<usize>,
    next_boundary: u64,
}

impl<'a> JsonCursor<'a> {
    fn new(input: &'a str) -> Self {
        Self {
            input,
            pos: 0,
            poll: None,
            cancelled: false,
            work_exceeded: None,
        }
    }

    fn with_work(
        input: &'a str,
        work: EvidenceWork<'a>,
        phase: SolveEvidencePhase,
        source_index: Option<usize>,
    ) -> Self {
        let cancelled = work.checkpoint(phase, source_index, 0).is_err();
        Self {
            input,
            pos: 0,
            poll: Some(JsonPoll {
                work,
                phase,
                source_index,
                next_boundary: u64::try_from(EVIDENCE_POLL_BYTES).expect("stride fits u64"),
            }),
            cancelled,
            work_exceeded: None,
        }
    }

    fn peek(&self) -> Option<u8> {
        self.input.as_bytes().get(self.pos).copied()
    }

    fn problem(&self, what: &str) -> String {
        format!("{what} at byte {}", self.pos)
    }

    fn expect(&mut self, expected: &str) -> Result<(), String> {
        if self.cancelled {
            return Err("solve evidence parsing stopped at phase entry".to_string());
        }
        if self.input[self.pos..].starts_with(expected) {
            self.advance(expected.len())
        } else {
            Err(self.problem(&format!("expected `{expected}`")))
        }
    }

    fn finish(&mut self) -> Result<(), String> {
        if self.cancelled {
            return Err("solve evidence parsing stopped at phase entry".to_string());
        }
        if self.pos == self.input.len() {
            self.poll_completion()?;
            Ok(())
        } else {
            Err(self.problem("trailing bytes after complete JSON value"))
        }
    }

    fn cancellation_observed(&self) -> bool {
        self.cancelled
    }

    fn checkpoint_current(&mut self, context: &str) -> Result<(), String> {
        if self.cancelled {
            return Err(format!("solve evidence parsing stopped {context}"));
        }
        let Some(poll) = &self.poll else {
            return Ok(());
        };
        let units =
            u64::try_from(self.pos).map_err(|_| "JSON cursor is outside u64".to_string())?;
        if poll
            .work
            .checkpoint(poll.phase, poll.source_index, units)
            .is_err()
        {
            self.cancelled = true;
            Err(format!("solve evidence parsing stopped {context}"))
        } else {
            Ok(())
        }
    }

    fn advance(&mut self, bytes: usize) -> Result<(), String> {
        if self.cancelled {
            return Err("solve evidence parsing stopped at phase entry".to_string());
        }
        let end = self
            .pos
            .checked_add(bytes)
            .ok_or_else(|| self.problem("JSON cursor offset overflow"))?;
        if end > self.input.len() {
            return Err(self.problem("JSON cursor advanced past input"));
        }
        if let Some(poll) = &mut self.poll {
            let end_units =
                u64::try_from(end).map_err(|_| "JSON cursor is outside u64".to_string())?;
            while poll.next_boundary <= end_units {
                if poll
                    .work
                    .checkpoint(poll.phase, poll.source_index, poll.next_boundary)
                    .is_err()
                {
                    self.cancelled = true;
                    return Err("solve evidence parsing stopped at a fixed checkpoint".to_string());
                }
                poll.next_boundary = poll
                    .next_boundary
                    .saturating_add(u64::try_from(EVIDENCE_POLL_BYTES).expect("stride fits u64"));
            }
            self.pos = end;
            if poll.work.plan.is_some()
                && poll
                    .work
                    .checkpoint(poll.phase, poll.source_index, end_units)
                    .is_err()
            {
                self.cancelled = true;
                return Err("solve evidence parsing stopped at a planned checkpoint".to_string());
            }
            return Ok(());
        }
        self.pos = end;
        Ok(())
    }

    fn poll_completion(&mut self) -> Result<(), String> {
        self.checkpoint_current("at completion")
    }

    fn parse_string(&mut self) -> Result<String, String> {
        self.expect("\"")?;
        let mut value = String::new();
        loop {
            let character = self.input[self.pos..]
                .chars()
                .next()
                .ok_or_else(|| self.problem("unterminated JSON string"))?;
            self.advance(character.len_utf8())?;
            match character {
                '"' => return Ok(value),
                '\\' => {
                    let escape = self.input[self.pos..]
                        .chars()
                        .next()
                        .ok_or_else(|| self.problem("unterminated JSON escape"))?;
                    self.advance(escape.len_utf8())?;
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
        let value =
            u16::from_str_radix(digits, 16).map_err(|_| self.problem("invalid Unicode escape"))?;
        self.advance(4)?;
        Ok(value)
    }

    fn parse_i64(&mut self) -> Result<i64, String> {
        let start = self.pos;
        if self.peek() == Some(b'-') {
            self.advance(1)?;
        }
        match self.peek() {
            Some(b'0') => {
                self.advance(1)?;
                if self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    return Err(self.problem("leading zero in JSON integer"));
                }
            }
            Some(b'1'..=b'9') => {
                self.advance(1)?;
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.advance(1)?;
                }
            }
            _ => return Err(self.problem("expected JSON integer")),
        }
        let spelling = &self.input[start..self.pos];
        if spelling.len() > MAX_JSON_INTEGER_BYTES {
            return Err(self.problem("JSON integer exceeds the canonical i64 width"));
        }
        spelling
            .parse::<i64>()
            .map_err(|_| self.problem("JSON integer is outside i64"))
    }

    fn parse_u64(&mut self) -> Result<u64, String> {
        let start = self.pos;
        match self.peek() {
            Some(b'0') => {
                self.advance(1)?;
                if self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    return Err(self.problem("leading zero in JSON unsigned integer"));
                }
            }
            Some(b'1'..=b'9') => {
                self.advance(1)?;
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.advance(1)?;
                }
            }
            _ => return Err(self.problem("expected JSON unsigned integer")),
        }
        let spelling = &self.input[start..self.pos];
        if spelling.len() > MAX_JSON_INTEGER_BYTES {
            return Err(self.problem("JSON unsigned integer exceeds the canonical u64 width"));
        }
        spelling
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
        if spelling.len() > MAX_CANONICAL_F64_BYTES {
            return Err(self.problem("JSON number exceeds the canonical finite-f64 width"));
        }
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
                self.advance(1)?;
                if self.peek() == Some(b'}') {
                    self.advance(1)?;
                    return Ok(());
                }
                loop {
                    let _ = self.parse_string()?;
                    self.expect(":")?;
                    self.skip_value(depth + 1)?;
                    match self.peek() {
                        Some(b',') => self.advance(1)?,
                        Some(b'}') => {
                            self.advance(1)?;
                            return Ok(());
                        }
                        _ => return Err(self.problem("expected `,` or `}` in JSON object")),
                    }
                }
            }
            Some(b'[') => {
                self.advance(1)?;
                if self.peek() == Some(b']') {
                    self.advance(1)?;
                    return Ok(());
                }
                loop {
                    self.skip_value(depth + 1)?;
                    match self.peek() {
                        Some(b',') => self.advance(1)?,
                        Some(b']') => {
                            self.advance(1)?;
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
            self.advance(1)?;
        }
        match self.peek() {
            Some(b'0') => {
                self.advance(1)?;
                if self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    return Err(self.problem("leading zero in JSON number"));
                }
            }
            Some(b'1'..=b'9') => {
                self.advance(1)?;
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.advance(1)?;
                }
            }
            _ => return Err(self.problem("invalid JSON number")),
        }
        if self.peek() == Some(b'.') {
            self.advance(1)?;
            let fraction_start = self.pos;
            while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                self.advance(1)?;
            }
            if self.pos == fraction_start {
                return Err(self.problem("empty JSON number fraction"));
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.advance(1)?;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.advance(1)?;
            }
            let exponent_start = self.pos;
            while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                self.advance(1)?;
            }
            if self.pos == exponent_start {
                return Err(self.problem("empty JSON number exponent"));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{
        EVIDENCE_POLL_BYTES, EvidenceReadError, EvidenceUtf8Error, EvidenceWork,
        InvocationWorkExceeded, InvocationWorkLedger, JsonCursor, MAX_CANONICAL_F64_BYTES,
        MAX_CANONICAL_PLY_VERTEX_LINE_BYTES, MAX_SOLVE_EVIDENCE_LABEL_BYTES,
        MAX_SOLVE_INVOCATION_WORK_BYTES, SolveCancellationPlan, SolveEvidencePhase,
        evidence_bytes_equal, evidence_utf8_string, materialize_evidence_artifact,
        parse_canonical_ply_face_line, parse_canonical_ply_vertex_line, preflight_canonical_ply,
        require_solve_evidence_label, validate_named_group_face_ranges, verify_evidence_artifact,
    };
    use fs_exec::CancelGate;
    use fs_io::NamedFaceGroup;
    use fs_ledger::{Ledger, hash_bytes};

    #[test]
    fn invocation_work_ledger_admits_exact_cap_and_refuses_plus_one_without_advancing() {
        let meter = InvocationWorkLedger::default();
        assert_eq!(
            meter
                .charge(MAX_SOLVE_INVOCATION_WORK_BYTES)
                .expect("the exact invocation-work cap is admitted"),
            MAX_SOLVE_INVOCATION_WORK_BYTES
        );
        let attempted = MAX_SOLVE_INVOCATION_WORK_BYTES
            .checked_add(1)
            .expect("1-GiB cap plus one fits u64");
        assert!(matches!(
            meter.charge(1),
            Err(InvocationWorkExceeded::CumulativeBytes {
                attempted: observed
            }) if observed == attempted
        ));
        assert_eq!(
            meter.used.get(),
            MAX_SOLVE_INVOCATION_WORK_BYTES,
            "a refused charge must not advance deterministic work accounting"
        );
        assert_eq!(
            meter.charge(0).expect("zero-byte retry remains admitted"),
            MAX_SOLVE_INVOCATION_WORK_BYTES
        );
    }

    #[test]
    fn utf8_copy_and_direct_comparison_charge_each_accepted_tile_once() {
        let gate = CancelGate::new_clock_free();
        let utf8_meter = InvocationWorkLedger::default();
        let decoded = evidence_utf8_string(
            b"abc",
            EvidenceWork::new(&gate, None, &utf8_meter),
            SolveEvidencePhase::ImportSummaryUtf8,
            None,
            "charged UTF-8",
        )
        .expect("valid UTF-8 is copied");
        assert_eq!(decoded, "abc");
        assert_eq!(utf8_meter.used.get(), 3);

        let compare_meter = InvocationWorkLedger::default();
        assert!(
            evidence_bytes_equal(
                b"abc",
                b"abc",
                EvidenceWork::new(&gate, None, &compare_meter),
                SolveEvidencePhase::ImportIrCanonicalCompare,
                None,
            )
            .expect("equal bytes compare")
        );
        assert_eq!(compare_meter.used.get(), 3);

        let invalid_meter = InvocationWorkLedger::default();
        assert!(matches!(
            evidence_utf8_string(
                &[0x80],
                EvidenceWork::new(&gate, None, &invalid_meter),
                SolveEvidencePhase::ImportSummaryUtf8,
                None,
                "invalid UTF-8",
            ),
            Err(EvidenceUtf8Error::Invalid(_))
        ));
        assert_eq!(
            invalid_meter.used.get(),
            0,
            "a rejected UTF-8 tile is not charged as accepted copy work"
        );
    }

    #[test]
    fn canonical_ply_finite_float_bound_covers_writer_extremes() {
        let smallest_subnormal = f64::from_bits(1);
        let cases = [
            ("maximum", f64::MAX, 309usize),
            ("minimum-positive", f64::MIN_POSITIVE, 326),
            ("smallest-subnormal", smallest_subnormal, 326),
            ("negative-smallest-subnormal", -smallest_subnormal, 327),
            ("negative-zero", -0.0, 2),
        ];
        for (label, value, expected_len) in cases {
            let spelling = value.to_string();
            assert_eq!(
                spelling.len(),
                expected_len,
                "{label} writer spelling changed"
            );
            assert!(
                spelling.len() <= MAX_CANONICAL_F64_BYTES,
                "{label} exceeds the solve finite-f64 token bound"
            );
            let line = format!("{spelling} {spelling} {spelling}");
            assert!(line.len() <= MAX_CANONICAL_PLY_VERTEX_LINE_BYTES);
            parse_canonical_ply_vertex_line(&line, true)
                .unwrap_or_else(|error| panic!("{label} writer line refused: {error}"));
        }
        parse_canonical_ply_face_line("3 0 1 2", 3, true).expect("canonical in-range triangle");
        assert!(
            parse_canonical_ply_face_line("3 00 1 2", 3, true)
                .expect_err("leading zero is not writer output")
                .contains("canonical u32")
        );
        assert!(
            parse_canonical_ply_face_line("3 0 1 3", 3, true)
                .expect_err("index equal to vertex count is out of range")
                .contains("outside")
        );
    }

    #[test]
    fn incremental_utf8_handles_boundaries_errors_empty_and_cancellation() {
        let gate = CancelGate::new_clock_free();
        let mut split_scalar = vec![b'a'; EVIDENCE_POLL_BYTES - 1];
        split_scalar.extend_from_slice("🦀z".as_bytes());
        let decoded = evidence_utf8_string(
            &split_scalar,
            EvidenceWork::unmetered(&gate, None),
            SolveEvidencePhase::ImportSummaryUtf8,
            None,
            "split scalar",
        )
        .expect("four-byte scalar crossing the tile boundary is valid");
        let mut expected = "a".repeat(EVIDENCE_POLL_BYTES - 1);
        expected.push_str("🦀z");
        assert_eq!(decoded, expected);

        let invalid_gate = CancelGate::new_clock_free();
        let mut invalid = vec![b'a'; EVIDENCE_POLL_BYTES];
        invalid.push(0x80);
        invalid.push(b'z');
        let error = evidence_utf8_string(
            &invalid,
            EvidenceWork::unmetered(&invalid_gate, None),
            SolveEvidencePhase::ImportSummaryUtf8,
            None,
            "invalid boundary",
        )
        .expect_err("a standalone continuation byte is invalid");
        assert!(
            matches!(
                &error,
                EvidenceUtf8Error::Invalid(problem)
                    if problem.contains(&format!("byte {EVIDENCE_POLL_BYTES}"))
            ),
            "{error:?}"
        );

        let empty_gate = CancelGate::new_clock_free();
        assert_eq!(
            evidence_utf8_string(
                b"",
                EvidenceWork::unmetered(&empty_gate, None),
                SolveEvidencePhase::ImportSummaryUtf8,
                None,
                "empty",
            )
            .expect("empty UTF-8"),
            ""
        );

        let cancel_gate = CancelGate::new_clock_free();
        let plan = SolveCancellationPlan::new(
            SolveEvidencePhase::ImportSummaryUtf8,
            None,
            EVIDENCE_POLL_BYTES as u64,
        );
        let cancelled = evidence_utf8_string(
            &vec![b'a'; EVIDENCE_POLL_BYTES + 1],
            EvidenceWork::unmetered(&cancel_gate, Some(&plan)),
            SolveEvidencePhase::ImportSummaryUtf8,
            None,
            "cancelled",
        );
        assert!(matches!(cancelled, Err(EvidenceUtf8Error::Cancelled)));
        assert!(plan.fired());

        let invalid_cancel_gate = CancelGate::new_clock_free();
        let invalid_plan =
            SolveCancellationPlan::new(SolveEvidencePhase::ImportSummaryUtf8, None, 1);
        let invalid_cancelled = evidence_utf8_string(
            &[0x80],
            EvidenceWork::unmetered(&invalid_cancel_gate, Some(&invalid_plan)),
            SolveEvidencePhase::ImportSummaryUtf8,
            None,
            "planned invalid",
        );
        assert!(matches!(
            invalid_cancelled,
            Err(EvidenceUtf8Error::Cancelled)
        ));
        assert!(invalid_plan.fired());
    }

    #[test]
    fn json_cursor_phase_entry_prioritizes_immediate_cancellation() {
        let malformed_gate = CancelGate::new_clock_free();
        malformed_gate.request();
        let mut malformed = JsonCursor::with_work(
            "!",
            EvidenceWork::unmetered(&malformed_gate, None),
            SolveEvidencePhase::ImportIrParse,
            None,
        );
        assert!(malformed.cancellation_observed());
        assert!(
            malformed
                .expect("{")
                .expect_err("requested gate stops before inspecting the first invalid byte")
                .contains("phase entry")
        );
        assert_eq!(malformed.pos, 0);

        let valid_gate = CancelGate::new_clock_free();
        valid_gate.request();
        let mut valid = JsonCursor::with_work(
            "{}",
            EvidenceWork::unmetered(&valid_gate, None),
            SolveEvidencePhase::ImportIrParse,
            None,
        );
        assert!(
            valid
                .expect("{")
                .expect_err("requested gate stops valid input at phase entry")
                .contains("phase entry")
        );
        assert_eq!(valid.pos, 0);

        let planned_gate = CancelGate::new_clock_free();
        let plan = SolveCancellationPlan::new(SolveEvidencePhase::ImportIrParse, None, 0);
        let planned = JsonCursor::with_work(
            "{}",
            EvidenceWork::unmetered(&planned_gate, Some(&plan)),
            SolveEvidencePhase::ImportIrParse,
            None,
        );
        assert!(plan.fired());
        assert!(planned.cancellation_observed());
        assert!(planned_gate.is_requested());

        let late_gate = CancelGate::new_clock_free();
        let mut late = JsonCursor::with_work(
            "{!",
            EvidenceWork::unmetered(&late_gate, None),
            SolveEvidencePhase::ImportIrParse,
            None,
        );
        late.expect("{").expect("first token");
        late_gate.request();
        let syntax = late.expect("\"schema\":");
        assert!(syntax.is_err(), "the next token is deliberately invalid");
        assert!(
            late.checkpoint_current("after deliberate syntax error")
                .is_err()
        );
        assert!(late.cancellation_observed());
    }

    #[test]
    fn evidence_reads_observe_phase_entry_before_empty_or_missing_dispatch() {
        let ledger = Ledger::open(":memory:").expect("ledger");
        let empty = ledger
            .put_artifact("empty-evidence", b"", None)
            .expect("empty artifact");
        let empty_gate = CancelGate::new_clock_free();
        let empty_plan = SolveCancellationPlan::new(SolveEvidencePhase::RawSourceRead, Some(0), 0);
        let empty_result = materialize_evidence_artifact(
            &ledger,
            EvidenceWork::unmetered(&empty_gate, Some(&empty_plan)),
            empty.hash,
            0,
            SolveEvidencePhase::RawSourceRead,
            Some(0),
        );
        assert!(matches!(empty_result, Err(EvidenceReadError::Cancelled)));
        assert!(empty_plan.fired());

        let missing_gate = CancelGate::new_clock_free();
        missing_gate.request();
        let missing = hash_bytes(b"missing-evidence-artifact");
        let missing_result = verify_evidence_artifact(
            &ledger,
            EvidenceWork::unmetered(&missing_gate, None),
            missing,
            1,
            SolveEvidencePhase::PromotionReceiptRead,
            Some(0),
        );
        assert!(matches!(missing_result, Err(EvidenceReadError::Cancelled)));
    }

    #[test]
    fn bounded_error_paths_poll_before_classifying_the_result() {
        let compare_gate = CancelGate::new_clock_free();
        let compare_plan =
            SolveCancellationPlan::new(SolveEvidencePhase::ResumeProjectCanonicalCompare, None, 1);
        let compared = evidence_bytes_equal(
            b"a",
            b"b",
            EvidenceWork::unmetered(&compare_gate, Some(&compare_plan)),
            SolveEvidencePhase::ResumeProjectCanonicalCompare,
            None,
        );
        assert!(compared.is_err(), "post-compare checkpoint wins");
        assert!(compare_plan.fired());

        let range_gate = CancelGate::new_clock_free();
        let range_plan =
            SolveCancellationPlan::new(SolveEvidencePhase::NamedGroupFaceRange, Some(0), 4);
        let groups = [NamedFaceGroup {
            name: "bounded-group".to_string(),
            faces: vec![4],
        }];
        let range = validate_named_group_face_ranges(
            &groups,
            4,
            EvidenceWork::unmetered(&range_gate, Some(&range_plan)),
            0,
        );
        assert!(matches!(range, Err(super::ImportSummaryError::Cancelled)));
        assert!(range_plan.fired());

        let preflight_gate = CancelGate::new_clock_free();
        let preflight_plan =
            SolveCancellationPlan::new(SolveEvidencePhase::PromotedMeshPreflight, Some(0), 3);
        let preflight = preflight_canonical_ply(
            b"bad",
            1,
            1,
            EvidenceWork::unmetered(&preflight_gate, Some(&preflight_plan)),
            0,
        );
        assert!(matches!(
            preflight,
            Err(super::ImportSummaryError::Cancelled)
        ));
        assert!(preflight_plan.fired());
    }

    #[test]
    fn solve_evidence_label_ceiling_is_explicit_at_the_boundary() {
        let exact = "l".repeat(MAX_SOLVE_EVIDENCE_LABEL_BYTES);
        require_solve_evidence_label(&exact, usize::MAX, "exact solve-label boundary fixture")
            .expect("exact 4-KiB label is admitted");
        let oversized = format!("{exact}x");
        assert!(
            require_solve_evidence_label(
                &oversized,
                usize::MAX,
                "oversized solve-label boundary fixture",
            )
            .expect_err("limit plus one is refused")
            .contains("solve label ceiling")
        );
    }
}
