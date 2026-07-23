//! HELM-side geometry import orchestration for the CLI product boundary.
//!
//! Canonical `.fsim` geometry rows name an already-imported receipt/content
//! identity, not a machine-local path. This module therefore accepts raw bytes
//! through an explicit caller-owned library, verifies those bytes against the
//! exact project row, runs fs-io quarantine/promotion and fs-project
//! assignment, and retains every successful artifact and receipt in one
//! fs-ledger operation. Physical source labels are diagnostic provenance only.

use core::fmt::Write as _;
use std::collections::BTreeMap;

use fs_evidence::vv::UnitId;
use fs_exec::Cx;
use fs_io::{
    AssignmentLimits, NamedFaceGroup, STEP_FACETED_DECODER_VERSION, STEP_IMPORT_SEMANTICS_VERSION,
    StepFacetedImportRefusal,
};
use fs_ledger::{
    ContentHash, EdgeRole, ExtensionTable, FiveExplicits, Ledger, LedgerError, OpOutcome,
};
use fs_project::{
    GeometryArtifact, GeometryResolution, ImportedMeshLibrary, ProjectSpec,
    geometry_source_identity, resolve_geometry_assignments,
};

const IMPORT_IR_SCHEMA: &str = "frankensim.cli.geometry-import.v1";
const IMPORT_SUMMARY_SCHEMA: &str = "frankensim.cli.geometry-import-receipt.v1";
const IMPORT_REFUSAL_SCHEMA: &str = "frankensim.cli.geometry-import-refusal.v1";
const STEP_IMPORT_RECEIPT_SCHEMA: &str = "frankensim.cli.faceted-step-import-receipt.v1";
const IMPORT_NO_CLAIM: &str = "the ledger binds exact raw bytes, fs-io receipts, one promoted finite tessellation, and the project assignment report; the project row's legacy FNV hook and a caller-supplied source label do not authenticate custody, physical/CAD sameness, continuum coverage, units, or topology beyond the retained lower-layer claims";

/// Explicit resource envelope for one project geometry-import attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GeometryImportLimits {
    /// Maximum number of project geometry rows and supplied sources.
    pub max_sources: usize,
    /// Maximum bytes in one raw source.
    pub max_source_bytes: usize,
    /// Maximum aggregate raw-source bytes.
    pub max_total_source_bytes: usize,
    /// Resource envelope forwarded to persistent assignment resolution.
    pub assignment: AssignmentLimits,
}

impl GeometryImportLimits {
    /// Product checkpoint defaults: 64 sources, 64 MiB each, 256 MiB total.
    pub const DEFAULT: Self = Self {
        max_sources: 64,
        max_source_bytes: 64 * 1024 * 1024,
        max_total_source_bytes: 256 * 1024 * 1024,
        assignment: AssignmentLimits::DEFAULT,
    };
}

impl Default for GeometryImportLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

#[derive(Debug)]
struct RawGeometrySource {
    label: String,
    bytes: Vec<u8>,
    policy: RawGeometryPolicy,
}

#[derive(Debug)]
enum RawGeometryPolicy {
    Mesh {
        length_unit: String,
        max_hole_edges: usize,
        named_groups: Vec<NamedFaceGroup>,
    },
    FacetedStep {
        root_id: u64,
        length_unit: String,
        target_h: f64,
        named_groups: Vec<NamedFaceGroup>,
    },
}

/// Caller-owned raw geometry inputs keyed by exact canonical project rows.
///
/// The library computes the same strong project-row key as fs-project. A
/// machine-local path may be used as `label`, but neither that spelling nor
/// insertion order participates in project identity.
#[derive(Debug, Default)]
pub struct RawGeometryLibrary {
    sources: BTreeMap<String, RawGeometrySource>,
}

impl RawGeometryLibrary {
    /// Construct an empty raw-source library.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Bind one raw STL/OBJ/PLY input to an exact project geometry row.
    ///
    /// Returns `true` when an earlier binding for the same exact row was
    /// replaced. Replacement is explicit at construction time; import itself
    /// observes one unambiguous source per row.
    pub fn insert_mesh(
        &mut self,
        artifact: &GeometryArtifact,
        label: impl Into<String>,
        bytes: Vec<u8>,
        length_unit: impl Into<String>,
        max_hole_edges: usize,
        named_groups: Vec<NamedFaceGroup>,
    ) -> bool {
        self.sources
            .insert(
                geometry_source_identity(artifact),
                RawGeometrySource {
                    label: label.into(),
                    bytes,
                    policy: RawGeometryPolicy::Mesh {
                        length_unit: length_unit.into(),
                        max_hole_edges,
                        named_groups,
                    },
                },
            )
            .is_some()
    }

    /// Bind one strict triangular faceted-STEP input to an exact project row.
    ///
    /// `root_id`, `length_unit`, and `target_h` are explicit replay policy.
    /// The lower-layer decoder ignores STEP representation/unit context and
    /// admits only its pinned root-reachable triangular resource subset.
    pub fn insert_faceted_step(
        &mut self,
        artifact: &GeometryArtifact,
        label: impl Into<String>,
        bytes: Vec<u8>,
        root_id: u64,
        length_unit: impl Into<String>,
        target_h: f64,
        named_groups: Vec<NamedFaceGroup>,
    ) -> bool {
        self.sources
            .insert(
                geometry_source_identity(artifact),
                RawGeometrySource {
                    label: label.into(),
                    bytes,
                    policy: RawGeometryPolicy::FacetedStep {
                        root_id,
                        length_unit: length_unit.into(),
                        target_h,
                        named_groups,
                    },
                },
            )
            .is_some()
    }
}

/// Ledger identities retained for one successfully imported geometry row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RetainedGeometryImport {
    /// Project-local artifact role.
    pub role: String,
    /// Caller-reported physical source label; provenance only, not identity.
    pub source_label: String,
    /// Strong identity of the exact project geometry row.
    pub source_identity: String,
    /// BLAKE3 content address of the exact hostile source bytes.
    pub raw_source: ContentHash,
    /// BLAKE3 content address of fs-io's exact promotion receipt JSON.
    pub promotion_receipt: ContentHash,
    /// BLAKE3 content address of the deterministic lossless PLY mesh payload.
    pub promoted_mesh: ContentHash,
    /// BLAKE3 content address of fs-project's exact retained assignment report.
    pub assignment_report: ContentHash,
    /// Stable imports-table key whose body is the exact fs-io receipt JSON.
    pub import_record: String,
}

/// Successful atomic import/assignment/retention result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeometryImportRun {
    /// Canonical fs-project content identity.
    pub project_hash: ContentHash,
    /// Terminal successful fs-ledger operation.
    pub op_id: i64,
    /// Content address of the complete orchestration summary.
    pub summary_artifact: ContentHash,
    /// Rows in canonical project geometry declaration order.
    pub artifacts: Vec<RetainedGeometryImport>,
    /// Deterministic assignment statistics table for CLI/progress rendering.
    pub assignment_table: String,
}

impl GeometryImportRun {
    /// Exact orchestration no-claim boundary.
    #[must_use]
    pub const fn no_claim() -> &'static str {
        IMPORT_NO_CLAIM
    }
}

/// Durable failure evidence when a raw source reached the ledger boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordedImportRefusal {
    /// Terminal error operation.
    pub op_id: i64,
    /// Retained structured diagnostic artifact.
    pub diagnostic_artifact: ContentHash,
    /// Retained exact raw sources in project declaration order.
    pub raw_sources: Vec<ContentHash>,
    /// Retained fs-io promotion/refusal receipts that were available.
    pub receipt_artifacts: Vec<ContentHash>,
    /// Retained promoted mesh payloads available before a later refusal.
    pub promoted_meshes: Vec<ContentHash>,
    /// Imports-table keys paired with `receipt_artifacts`.
    pub import_records: Vec<String>,
}

/// Stable actionable import refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeometryImportRefusal {
    /// Machine-facing refusal code.
    pub code: &'static str,
    /// Project geometry role, if the refusal is source-specific.
    pub role: Option<String>,
    /// Diagnosis.
    pub what: String,
    /// Actionable next step.
    pub fix: String,
    /// Durable evidence, when the attempt reached a recordable source stage.
    pub recorded: Option<RecordedImportRefusal>,
}

impl core::fmt::Display for GeometryImportRefusal {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if let Some(role) = &self.role {
            write!(
                formatter,
                "{} for geometry role `{role}`: {}",
                self.code, self.what
            )
        } else {
            write!(formatter, "{}: {}", self.code, self.what)
        }
    }
}

impl std::error::Error for GeometryImportRefusal {}

struct PreparedImport {
    role: String,
    source_identity: String,
    label: String,
    raw_bytes: Vec<u8>,
    receipt_json: String,
    promoted_mesh_bytes: Vec<u8>,
}

struct RefusalEvidence<'a> {
    source_identity: &'a str,
    raw_bytes: &'a [u8],
    receipt_json: Option<&'a str>,
    promoted_mesh_bytes: Option<&'a [u8]>,
}

/// Replay every declared raw geometry source through quarantine, promotion, persistent
/// assignment, and one atomic ledger operation.
///
/// This checkpoint admits STL, OBJ, PLY, and the strict triangular faceted-STEP
/// subset. STEP callers must supply an explicit root, length unit, and sampling
/// spacing; all caller policy is frozen into the ledger operation IR.
///
/// # Errors
///
/// [`GeometryImportRefusal`] for project/source mismatch, resource admission,
/// parser/promotion/assignment refusal, cancellation, or ledger failure.
pub fn import_project_geometry(
    project: &ProjectSpec,
    raw: &RawGeometryLibrary,
    ledger: &Ledger,
    limits: GeometryImportLimits,
    cx: &Cx<'_>,
) -> Result<GeometryImportRun, GeometryImportRefusal> {
    if ledger.in_transaction() {
        return Err(refusal(
            "cli-import-ledger-transaction",
            None,
            "geometry import requires ownership of one atomic ledger transaction",
            "commit or roll back the caller transaction, then retry",
        ));
    }
    checkpoint(cx, "cli-import-cancelled", "geometry-import-entry")?;
    validate_limits(limits)?;

    let findings = project.validate();
    if let Some(finding) = findings.first() {
        return Err(refusal(
            finding.code,
            None,
            finding.what.clone(),
            finding.fix.clone(),
        ));
    }
    let geometry = project.geometry.as_ref().ok_or_else(|| {
        refusal(
            "project-geometry-missing",
            None,
            "project has no geometry section",
            "declare imported geometry receipt rows before import replay",
        )
    })?;
    if geometry.len() > limits.max_sources {
        return Err(refusal(
            "cli-import-resource-bound",
            None,
            format!(
                "project declares {} geometry rows; the explicit cap is {}",
                geometry.len(),
                limits.max_sources
            ),
            "raise the import source cap within the run budget or split the project",
        ));
    }
    if raw.sources.len() != geometry.len() {
        return Err(refusal(
            "cli-import-source-set",
            None,
            format!(
                "project declares {} geometry rows but the caller supplied {} raw sources",
                geometry.len(),
                raw.sources.len()
            ),
            "supply exactly one raw source for every exact project geometry row and no extras",
        ));
    }
    for artifact in geometry {
        let source_identity = geometry_source_identity(artifact);
        if !raw.sources.contains_key(&source_identity) {
            return Err(refusal(
                "cli-import-source-missing",
                Some(&artifact.role),
                format!("no raw source is bound to project geometry identity `{source_identity}`"),
                "bind the physical file to this exact project geometry row",
            ));
        }
    }

    let project_json = fs_project::print_json(project).map_err(project_codec_refusal)?;
    let canonical = fs_project::print_sexpr(project).map_err(project_codec_refusal)?;
    let project_hash = fs_project::canonical_hash(canonical.as_bytes());
    let import_request_ir = import_ir(&project_json, geometry, raw, limits)?;
    let mut total_bytes = 0usize;
    let mut mesh_library = ImportedMeshLibrary::new();
    let mut prepared = Vec::with_capacity(geometry.len());

    for artifact in geometry {
        checkpoint(cx, "cli-import-cancelled", "geometry-import-source")?;
        let source_identity = geometry_source_identity(artifact);
        let Some(source) = raw.sources.get(&source_identity) else {
            return Err(refusal(
                "cli-import-source-missing",
                Some(&artifact.role),
                format!("no raw source is bound to project geometry identity `{source_identity}`"),
                "bind the physical file to this exact project geometry row",
            ));
        };
        validate_source_envelope(artifact, source, limits, &mut total_bytes)?;

        let input = match &source.policy {
            RawGeometryPolicy::Mesh {
                length_unit,
                max_hole_edges,
                named_groups,
            } => prepare_mesh_source(
                ledger,
                project,
                &import_request_ir,
                project_hash,
                artifact,
                source,
                &source_identity,
                length_unit,
                *max_hole_edges,
                named_groups,
                &mut mesh_library,
                cx,
            )?,
            RawGeometryPolicy::FacetedStep {
                root_id,
                length_unit,
                target_h,
                named_groups,
            } => prepare_faceted_step_source(
                ledger,
                project,
                &import_request_ir,
                project_hash,
                artifact,
                source,
                &source_identity,
                *root_id,
                length_unit,
                *target_h,
                named_groups,
                &mut mesh_library,
                cx,
            )?,
        };
        prepared.push(input);
    }

    let resolution = resolve_geometry_assignments(project, &mesh_library, limits.assignment, cx);
    if !resolution.admissible() {
        let first = &resolution.violations[0];
        let error = refusal(
            first.code,
            None,
            resolution
                .violations
                .iter()
                .map(|violation| violation.what.as_str())
                .collect::<Vec<_>>()
                .join("; "),
            resolution
                .violations
                .iter()
                .map(|violation| violation.fix.as_str())
                .collect::<Vec<_>>()
                .join("; "),
        );
        return Err(record_prepared_or_ledger_refusal(
            ledger,
            project,
            &import_request_ir,
            project_hash,
            &prepared,
            error,
        ));
    }
    if resolution.artifacts.len() != prepared.len() {
        let error = refusal(
            "cli-import-assignment-result",
            None,
            format!(
                "assignment resolution returned {} artifact reports for {} prepared imports",
                resolution.artifacts.len(),
                prepared.len()
            ),
            "treat this as an fs-project conformance failure and repair the resolver before retrying",
        );
        return Err(record_prepared_or_ledger_refusal(
            ledger,
            project,
            &import_request_ir,
            project_hash,
            &prepared,
            error,
        ));
    }
    checkpoint(cx, "cli-import-cancelled", "geometry-import-ledger")?;
    persist_success(
        ledger,
        project,
        &import_request_ir,
        project_hash,
        &prepared,
        &resolution,
    )
}

#[allow(clippy::too_many_arguments)]
fn prepare_mesh_source(
    ledger: &Ledger,
    project: &ProjectSpec,
    import_request_ir: &str,
    project_hash: ContentHash,
    artifact: &GeometryArtifact,
    source: &RawGeometrySource,
    source_identity: &str,
    length_unit: &str,
    max_hole_edges: usize,
    named_groups: &[NamedFaceGroup],
    mesh_library: &mut ImportedMeshLibrary,
    cx: &Cx<'_>,
) -> Result<PreparedImport, GeometryImportRefusal> {
    let format: &'static str = match artifact.format.as_str() {
        "stl" => "stl",
        "obj" => "obj",
        "ply" => "ply",
        _ => {
            let error = refusal(
                "cli-import-policy-format-mismatch",
                Some(&artifact.role),
                format!(
                    "mesh import policy cannot replay declared format `{}`",
                    artifact.format
                ),
                "bind STL, OBJ, or PLY rows with insert_mesh; bind STEP rows with insert_faceted_step",
            );
            return Err(record_or_ledger_refusal(
                ledger,
                project,
                import_request_ir,
                project_hash,
                artifact,
                source,
                None,
                error,
            ));
        }
    };

    let quarantined = match fs_io::quarantine::import_mesh(&source.bytes, format) {
        Ok(quarantined) => quarantined,
        Err(error) => {
            let error = refusal(
                "cli-import-parse",
                Some(&artifact.role),
                error.to_string(),
                "repair or re-export the named raw file in the declared bounded format",
            );
            return Err(record_or_ledger_refusal(
                ledger,
                project,
                import_request_ir,
                project_hash,
                artifact,
                source,
                None,
                error,
            ));
        }
    };

    if quarantined.source_receipt.source_hash != artifact.source_hash {
        let error = refusal(
            "cli-import-source-hash-mismatch",
            Some(&artifact.role),
            format!(
                "raw source hashes to {:016x}, but the project receipt row pins {:016x}",
                quarantined.source_receipt.source_hash, artifact.source_hash
            ),
            "select the exact imported bytes named by the project or update the project through an explicit import/migration receipt",
        );
        return Err(record_or_ledger_refusal(
            ledger,
            project,
            import_request_ir,
            project_hash,
            artifact,
            source,
            None,
            error,
        ));
    }
    if quarantined.source_receipt.parser_version != artifact.parser_version {
        let error = refusal(
            "cli-import-parser-version-mismatch",
            Some(&artifact.role),
            format!(
                "the current fs-io parser is `{}`, but the project receipt row pins `{}`",
                quarantined.source_receipt.parser_version, artifact.parser_version
            ),
            "replay with the pinned parser constellation or explicitly migrate the import receipt",
        );
        return Err(record_or_ledger_refusal(
            ledger,
            project,
            import_request_ir,
            project_hash,
            artifact,
            source,
            None,
            error,
        ));
    }
    let face_ordinals_will_move = !named_groups.is_empty()
        && quarantined
            .defects
            .iter()
            .any(|defect| matches!(defect.class, "duplicate-face" | "degenerate-face"));

    checkpoint(cx, "cli-import-cancelled", "geometry-import-promotion")?;
    let (promoted, receipt_json) = match fs_io::promote(quarantined, max_hole_edges) {
        Ok(promoted) => promoted,
        Err(promotion) => {
            let error = refusal(
                "cli-import-promotion-refused",
                Some(&artifact.role),
                format!("blocking defects: {}", promotion.blocking.join(", ")),
                promotion.fixes.join("; "),
            );
            return Err(record_or_ledger_refusal(
                ledger,
                project,
                import_request_ir,
                project_hash,
                artifact,
                source,
                Some(&promotion.receipt_json),
                error,
            ));
        }
    };
    checkpoint(cx, "cli-import-cancelled", "geometry-import-post-promotion")?;
    if face_ordinals_will_move {
        let error = refusal(
            "cli-import-group-remap-unavailable",
            Some(&artifact.role),
            "promotion removed duplicate or degenerate faces, so pre-repair named-group ordinals no longer identify the promoted soup",
            "re-export stable named groups against the promoted mesh, omit the groups and use geometric selectors, or supply an adapter with an explicit repair remap",
        );
        return Err(record_or_ledger_refusal(
            ledger,
            project,
            import_request_ir,
            project_hash,
            artifact,
            source,
            Some(&receipt_json),
            error,
        ));
    }

    let soup = promoted.value;
    let promoted_mesh_bytes = fs_io::ply::write_ply(&soup).into_bytes();
    mesh_library.insert(
        artifact,
        soup,
        length_unit.to_string(),
        named_groups.to_vec(),
    );
    Ok(PreparedImport {
        role: artifact.role.clone(),
        source_identity: source_identity.to_string(),
        label: source.label.clone(),
        raw_bytes: source.bytes.clone(),
        receipt_json,
        promoted_mesh_bytes,
    })
}

#[allow(clippy::too_many_arguments)]
fn prepare_faceted_step_source(
    ledger: &Ledger,
    project: &ProjectSpec,
    import_request_ir: &str,
    project_hash: ContentHash,
    artifact: &GeometryArtifact,
    source: &RawGeometrySource,
    source_identity: &str,
    root_id: u64,
    length_unit: &str,
    target_h: f64,
    named_groups: &[NamedFaceGroup],
    mesh_library: &mut ImportedMeshLibrary,
    cx: &Cx<'_>,
) -> Result<PreparedImport, GeometryImportRefusal> {
    if artifact.format != "step" {
        let error = refusal(
            "cli-import-policy-format-mismatch",
            Some(&artifact.role),
            format!(
                "faceted-STEP import policy cannot replay declared format `{}`",
                artifact.format
            ),
            "declare format `step`, or bind mesh formats with insert_mesh",
        );
        return Err(record_or_ledger_refusal(
            ledger,
            project,
            import_request_ir,
            project_hash,
            artifact,
            source,
            None,
            error,
        ));
    }

    let parsed = match fs_io::parse_step(&source.bytes) {
        Ok(parsed) => parsed,
        Err(error) => {
            let error = refusal(
                "cli-import-step-parse",
                Some(&artifact.role),
                error.to_string(),
                "repair or re-export the hostile Part-21 source within fs-io's explicit syntax limits",
            );
            return Err(record_or_ledger_refusal(
                ledger,
                project,
                import_request_ir,
                project_hash,
                artifact,
                source,
                None,
                error,
            ));
        }
    };
    if parsed.receipt().source_fingerprint() != artifact.source_hash {
        let error = refusal(
            "cli-import-source-hash-mismatch",
            Some(&artifact.role),
            format!(
                "raw STEP source hashes to {:016x}, but the project receipt row pins {:016x}",
                parsed.receipt().source_fingerprint(),
                artifact.source_hash
            ),
            "select the exact Part-21 bytes named by the project or explicitly migrate the import receipt",
        );
        return Err(record_or_ledger_refusal(
            ledger,
            project,
            import_request_ir,
            project_hash,
            artifact,
            source,
            None,
            error,
        ));
    }
    if artifact.parser_version != STEP_FACETED_DECODER_VERSION {
        let error = refusal(
            "cli-import-parser-version-mismatch",
            Some(&artifact.role),
            format!(
                "the current strict faceted decoder is `{STEP_FACETED_DECODER_VERSION}`, but the project receipt row pins `{}`",
                artifact.parser_version
            ),
            "replay with the pinned decoder constellation or explicitly migrate the import receipt",
        );
        return Err(record_or_ledger_refusal(
            ledger,
            project,
            import_request_ir,
            project_hash,
            artifact,
            source,
            Some(&parsed.receipt().to_json()),
            error,
        ));
    }

    let unit = match UnitId::try_new(length_unit.to_string()) {
        Ok(unit) => unit,
        Err(error) => {
            let refusal = refusal(
                "cli-import-step-unit",
                Some(&artifact.role),
                format!("STEP length-unit identity was refused: {error}"),
                "supply a bounded non-blank machine unit identity matching the project assignment",
            );
            return Err(record_or_ledger_refusal(
                ledger,
                project,
                import_request_ir,
                project_hash,
                artifact,
                source,
                Some(&parsed.receipt().to_json()),
                refusal,
            ));
        }
    };
    if root_id == 0 {
        let error = refusal(
            "cli-import-step-root",
            Some(&artifact.role),
            "faceted-STEP root instance ID must be positive",
            "select the exact positive FACETED_BREP instance ID from the admitted Part-21 resource",
        );
        return Err(record_or_ledger_refusal(
            ledger,
            project,
            import_request_ir,
            project_hash,
            artifact,
            source,
            Some(&parsed.receipt().to_json()),
            error,
        ));
    }
    if !target_h.is_finite() || target_h <= 0.0 {
        let error = refusal(
            "cli-import-step-spacing",
            Some(&artifact.role),
            format!("STEP sampling spacing must be finite and positive; got {target_h}"),
            "supply a finite positive target_h expressed in the declared length unit",
        );
        return Err(record_or_ledger_refusal(
            ledger,
            project,
            import_request_ir,
            project_hash,
            artifact,
            source,
            Some(&parsed.receipt().to_json()),
            error,
        ));
    }

    checkpoint(cx, "cli-import-cancelled", "geometry-import-step")?;
    let outcome = match fs_io::import_faceted_brep(&parsed, root_id, unit, target_h, cx) {
        Ok(outcome) => outcome,
        Err(StepFacetedImportRefusal::Decode(error)) => {
            let refusal = refusal(
                "cli-import-step-decode",
                Some(&artifact.role),
                error.to_string(),
                "select a supported root-reachable triangular FACETED_BREP closure and repair the named STEP relationship",
            );
            return Err(record_or_ledger_refusal(
                ledger,
                project,
                import_request_ir,
                project_hash,
                artifact,
                source,
                Some(&parsed.receipt().to_json()),
                refusal,
            ));
        }
        Err(StepFacetedImportRefusal::Import {
            decoder_receipt,
            error,
        }) => {
            let receipt_json =
                step_failure_receipt_json(&decoder_receipt.to_json(), &error.to_string());
            let refusal = refusal(
                "cli-import-step-promotion-refused",
                Some(&artifact.role),
                error.to_string(),
                "repair the decoded root-reachable closure or adjust the finite positive sampling policy without weakening the retained no-claim boundary",
            );
            return Err(record_or_ledger_refusal(
                ledger,
                project,
                import_request_ir,
                project_hash,
                artifact,
                source,
                Some(&receipt_json),
                refusal,
            ));
        }
    };
    checkpoint(
        cx,
        "cli-import-cancelled",
        "geometry-import-step-post-promotion",
    )?;

    let receipt_json = step_success_receipt_json(
        &outcome.decoder_receipt().to_json(),
        &outcome.import().receipt().to_json(),
    );
    if !named_groups.is_empty()
        && outcome
            .import()
            .receipt()
            .repairs()
            .iter()
            .any(|repair| repair.action == "removed")
    {
        let error = refusal(
            "cli-import-group-remap-unavailable",
            Some(&artifact.role),
            "the STEP topology handoff removed faces, so decoder-ordinal named groups no longer identify the repaired soup",
            "omit STEP named groups and use geometric selectors, or supply a format adapter with an explicit decoder-to-repaired face remap",
        );
        return Err(record_or_ledger_refusal(
            ledger,
            project,
            import_request_ir,
            project_hash,
            artifact,
            source,
            Some(&receipt_json),
            error,
        ));
    }
    let soup = outcome.import().repaired_soup().clone();
    let promoted_mesh_bytes = fs_io::ply::write_ply(&soup).into_bytes();
    mesh_library.insert(
        artifact,
        soup,
        length_unit.to_string(),
        named_groups.to_vec(),
    );
    Ok(PreparedImport {
        role: artifact.role.clone(),
        source_identity: source_identity.to_string(),
        label: source.label.clone(),
        raw_bytes: source.bytes.clone(),
        receipt_json,
        promoted_mesh_bytes,
    })
}

fn validate_limits(limits: GeometryImportLimits) -> Result<(), GeometryImportRefusal> {
    if limits.max_sources == 0 || limits.max_source_bytes == 0 || limits.max_total_source_bytes == 0
    {
        return Err(refusal(
            "cli-import-resource-bound",
            None,
            "import source count and byte limits must all be positive",
            "provide a nonzero explicit import resource envelope",
        ));
    }
    Ok(())
}

fn validate_source_envelope(
    artifact: &GeometryArtifact,
    source: &RawGeometrySource,
    limits: GeometryImportLimits,
    total_bytes: &mut usize,
) -> Result<(), GeometryImportRefusal> {
    if source.label.is_empty()
        || source.label.len() > 4096
        || source.label.chars().any(char::is_control)
    {
        return Err(refusal(
            "cli-import-source-label",
            Some(&artifact.role),
            "source label must be 1..=4096 UTF-8 bytes without control characters",
            "use a bounded printable path or caller provenance label",
        ));
    }
    if source.bytes.len() > limits.max_source_bytes {
        return Err(refusal(
            "cli-import-resource-bound",
            Some(&artifact.role),
            format!(
                "raw source is {} bytes; the per-source cap is {}",
                source.bytes.len(),
                limits.max_source_bytes
            ),
            "raise the explicit cap within the memory budget or reduce the input",
        ));
    }
    *total_bytes = total_bytes.checked_add(source.bytes.len()).ok_or_else(|| {
        refusal(
            "cli-import-resource-bound",
            Some(&artifact.role),
            "aggregate raw-source byte count overflowed the platform range",
            "reduce the imported source set",
        )
    })?;
    if *total_bytes > limits.max_total_source_bytes {
        return Err(refusal(
            "cli-import-resource-bound",
            Some(&artifact.role),
            format!(
                "aggregate raw sources reached {total_bytes} bytes; the cap is {}",
                limits.max_total_source_bytes
            ),
            "raise the explicit aggregate cap within the memory budget or split the project",
        ));
    }
    Ok(())
}

fn persist_success(
    ledger: &Ledger,
    project: &ProjectSpec,
    import_request_ir: &str,
    project_hash: ContentHash,
    prepared: &[PreparedImport],
    resolution: &GeometryResolution,
) -> Result<GeometryImportRun, GeometryImportRefusal> {
    let (versions, budget, capability, seed) = explicits(project)?;
    ledger.begin().map_err(ledger_refusal)?;
    let result = (|| -> Result<GeometryImportRun, LedgerError> {
        let op = ledger.begin_op(
            None,
            import_request_ir,
            &FiveExplicits {
                seed: &seed,
                versions: &versions,
                budget: &budget,
                capability: &capability,
            },
            0,
        )?;
        let mut retained = Vec::with_capacity(prepared.len());
        for (input, assignment) in prepared.iter().zip(&resolution.artifacts) {
            let raw = ledger.put_artifact("geometry-source", &input.raw_bytes, None)?;
            let receipt = ledger.put_artifact(
                "geometry-import-receipt",
                input.receipt_json.as_bytes(),
                None,
            )?;
            let mesh =
                ledger.put_artifact("geometry-mesh-ply", &input.promoted_mesh_bytes, None)?;
            let report = ledger.put_artifact(
                "geometry-assignment-report",
                &assignment.report_bytes,
                None,
            )?;
            ledger.link(op, &raw.hash, EdgeRole::In)?;
            for hash in [receipt.hash, mesh.hash, report.hash] {
                ledger.link(op, &hash, EdgeRole::Out)?;
            }
            let import_record = import_record_name(&project_hash, &input.source_identity);
            ledger.put_extension(ExtensionTable::Imports, &import_record, &input.receipt_json)?;
            retained.push(RetainedGeometryImport {
                role: input.role.clone(),
                source_label: input.label.clone(),
                source_identity: input.source_identity.clone(),
                raw_source: raw.hash,
                promotion_receipt: receipt.hash,
                promoted_mesh: mesh.hash,
                assignment_report: report.hash,
                import_record,
            });
        }
        let assignment_table = resolution.render_table();
        let summary = summary_json(project_hash, &retained, &assignment_table);
        let summary_artifact =
            ledger.put_artifact("geometry-import-run-receipt", summary.as_bytes(), None)?;
        ledger.link(op, &summary_artifact.hash, EdgeRole::Out)?;
        ledger.finish_op(op, OpOutcome::Ok, None, 1)?;
        Ok(GeometryImportRun {
            project_hash,
            op_id: op,
            summary_artifact: summary_artifact.hash,
            artifacts: retained,
            assignment_table,
        })
    })();
    finish_owned_transaction(ledger, result)
}

#[allow(clippy::too_many_arguments)]
fn record_or_ledger_refusal(
    ledger: &Ledger,
    project: &ProjectSpec,
    import_request_ir: &str,
    project_hash: ContentHash,
    artifact: &GeometryArtifact,
    source: &RawGeometrySource,
    receipt_json: Option<&str>,
    mut error: GeometryImportRefusal,
) -> GeometryImportRefusal {
    let source_identity = geometry_source_identity(artifact);
    let evidence = [RefusalEvidence {
        source_identity: &source_identity,
        raw_bytes: &source.bytes,
        receipt_json,
        promoted_mesh_bytes: None,
    }];
    match persist_refusal(
        ledger,
        project,
        import_request_ir,
        project_hash,
        &evidence,
        &error,
    ) {
        Ok(recorded) => {
            error.recorded = Some(recorded);
            error
        }
        Err(ledger_error) => GeometryImportRefusal {
            code: "cli-import-ledger",
            role: Some(artifact.role.clone()),
            what: format!(
                "the import refused with `{}` but durable refusal recording also failed: {}",
                error.code, ledger_error.what
            ),
            fix: format!(
                "{}; then repair the ledger failure and retry so the refusal is durably visible",
                error.fix
            ),
            recorded: None,
        },
    }
}

fn record_prepared_or_ledger_refusal(
    ledger: &Ledger,
    project: &ProjectSpec,
    import_request_ir: &str,
    project_hash: ContentHash,
    prepared: &[PreparedImport],
    mut error: GeometryImportRefusal,
) -> GeometryImportRefusal {
    let evidence = prepared
        .iter()
        .map(|input| RefusalEvidence {
            source_identity: &input.source_identity,
            raw_bytes: &input.raw_bytes,
            receipt_json: Some(&input.receipt_json),
            promoted_mesh_bytes: Some(&input.promoted_mesh_bytes),
        })
        .collect::<Vec<_>>();
    match persist_refusal(
        ledger,
        project,
        import_request_ir,
        project_hash,
        &evidence,
        &error,
    ) {
        Ok(recorded) => {
            error.recorded = Some(recorded);
            error
        }
        Err(ledger_error) => GeometryImportRefusal {
            code: "cli-import-ledger",
            role: error.role,
            what: format!(
                "the import refused with `{}` but durable refusal recording also failed: {}",
                error.code, ledger_error.what
            ),
            fix: format!(
                "{}; then repair the ledger failure and retry so the refusal is durably visible",
                error.fix
            ),
            recorded: None,
        },
    }
}

fn persist_refusal(
    ledger: &Ledger,
    project: &ProjectSpec,
    import_request_ir: &str,
    project_hash: ContentHash,
    evidence: &[RefusalEvidence<'_>],
    error: &GeometryImportRefusal,
) -> Result<RecordedImportRefusal, GeometryImportRefusal> {
    let (versions, budget, capability, seed) = explicits(project)?;
    let diagnostic = refusal_json(error);
    ledger.begin().map_err(ledger_refusal)?;
    let result = (|| -> Result<RecordedImportRefusal, LedgerError> {
        let op = ledger.begin_op(
            None,
            import_request_ir,
            &FiveExplicits {
                seed: &seed,
                versions: &versions,
                budget: &budget,
                capability: &capability,
            },
            0,
        )?;
        let diagnostic_artifact =
            ledger.put_artifact("geometry-import-refusal", diagnostic.as_bytes(), None)?;
        ledger.link(op, &diagnostic_artifact.hash, EdgeRole::Out)?;
        let mut raw_sources = Vec::with_capacity(evidence.len());
        let mut receipt_artifacts = Vec::with_capacity(evidence.len());
        let mut promoted_meshes = Vec::with_capacity(evidence.len());
        let mut import_records = Vec::with_capacity(evidence.len());
        for item in evidence {
            let raw = ledger.put_artifact("geometry-source", item.raw_bytes, None)?;
            ledger.link(op, &raw.hash, EdgeRole::In)?;
            raw_sources.push(raw.hash);
            if let Some(receipt_json) = item.receipt_json {
                let receipt = ledger.put_artifact(
                    "geometry-import-receipt",
                    receipt_json.as_bytes(),
                    None,
                )?;
                ledger.link(op, &receipt.hash, EdgeRole::Out)?;
                let name = import_record_name(&project_hash, item.source_identity);
                ledger.put_extension(ExtensionTable::Imports, &name, receipt_json)?;
                receipt_artifacts.push(receipt.hash);
                import_records.push(name);
            }
            if let Some(mesh_bytes) = item.promoted_mesh_bytes {
                let mesh = ledger.put_artifact("geometry-mesh-ply", mesh_bytes, None)?;
                ledger.link(op, &mesh.hash, EdgeRole::Out)?;
                promoted_meshes.push(mesh.hash);
            }
        }
        ledger.finish_op(op, OpOutcome::Error, Some(&diagnostic), 1)?;
        Ok(RecordedImportRefusal {
            op_id: op,
            diagnostic_artifact: diagnostic_artifact.hash,
            raw_sources,
            receipt_artifacts,
            promoted_meshes,
            import_records,
        })
    })();
    finish_owned_transaction(ledger, result)
}

fn finish_owned_transaction<T>(
    ledger: &Ledger,
    result: Result<T, LedgerError>,
) -> Result<T, GeometryImportRefusal> {
    match result {
        Ok(value) => match ledger.commit() {
            Ok(()) => Ok(value),
            Err(error) => {
                let rollback = ledger.rollback();
                Err(transaction_refusal(error, rollback.err()))
            }
        },
        Err(error) => {
            let rollback = ledger.rollback();
            Err(transaction_refusal(error, rollback.err()))
        }
    }
}

fn transaction_refusal(
    primary: LedgerError,
    rollback: Option<LedgerError>,
) -> GeometryImportRefusal {
    let mut what = format!("ledger transaction failed: {primary}");
    if let Some(rollback) = rollback {
        let _ = write!(what, "; rollback also failed: {rollback}");
    }
    refusal(
        "cli-import-ledger",
        None,
        what,
        "repair the ledger or contention failure, verify ledger lint/integrity, and retry",
    )
}

fn explicits(
    project: &ProjectSpec,
) -> Result<(String, String, String, [u8; 8]), GeometryImportRefusal> {
    let versions = project.versions.as_ref().ok_or_else(|| {
        refusal(
            "project-versions-missing",
            None,
            "project versions are unavailable at import",
            "run strict project validation before import",
        )
    })?;
    let budgets = project.budgets.as_ref().ok_or_else(|| {
        refusal(
            "project-budgets-missing",
            None,
            "project budgets are unavailable at import",
            "run strict project validation before import",
        )
    })?;
    let capabilities = project.capabilities.as_ref().ok_or_else(|| {
        refusal(
            "project-capabilities-missing",
            None,
            "project capabilities are unavailable at import",
            "run strict project validation before import",
        )
    })?;
    let seeds = project.seeds.as_ref().ok_or_else(|| {
        refusal(
            "project-seeds-missing",
            None,
            "project seed is unavailable at import",
            "run strict project validation before import",
        )
    })?;

    let versions = format!(
        "{{\"schema\":{},\"constellation\":{},\"workspace\":{},\"fs_io\":{},\"step_faceted_decoder\":{},\"step_import_semantics\":{}}}",
        versions.schema,
        json_string(&versions.constellation),
        json_string(&versions.workspace),
        json_string(fs_io::VERSION),
        json_string(STEP_FACETED_DECODER_VERSION),
        json_string(STEP_IMPORT_SEMANTICS_VERSION),
    );
    let budget = format!(
        "{{\"solve_time\":{},\"solve_time_dims\":{:?},\"memory_bytes\":{},\"accuracy_rel\":{}}}",
        budgets.solve_time.value,
        budgets.solve_time.dims.0,
        budgets.memory_bytes,
        budgets.accuracy_rel,
    );
    let mut capability = String::from("{\"requested\":[");
    for (index, item) in capabilities.iter().enumerate() {
        if index > 0 {
            capability.push(',');
        }
        capability.push_str(&json_string(item));
    }
    capability.push_str("]}");
    Ok((versions, budget, capability, seeds.root.to_le_bytes()))
}

fn import_ir(
    project_json: &str,
    geometry: &[GeometryArtifact],
    raw: &RawGeometryLibrary,
    limits: GeometryImportLimits,
) -> Result<String, GeometryImportRefusal> {
    let assignment = limits.assignment;
    let mut out = format!(
        "{{\"schema\":{},\"project\":{project_json},\"limits\":{{\
         \"max_sources\":{},\"max_source_bytes\":{},\"max_total_source_bytes\":{},\
         \"assignment\":{{\"max_mesh_vertices\":{},\"max_mesh_faces\":{},\
         \"max_requests\":{},\"max_named_groups\":{},\"max_group_faces\":{},\
         \"max_selected_faces\":{},\"max_predicate_tests\":{},\"max_label_bytes\":{}}}}},\
         \"sources\":[",
        json_string(IMPORT_IR_SCHEMA),
        limits.max_sources,
        limits.max_source_bytes,
        limits.max_total_source_bytes,
        assignment.max_mesh_vertices,
        assignment.max_mesh_faces,
        assignment.max_requests,
        assignment.max_named_groups,
        assignment.max_group_faces,
        assignment.max_selected_faces,
        assignment.max_predicate_tests,
        assignment.max_label_bytes,
    );
    for (index, artifact) in geometry.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let source_identity = geometry_source_identity(artifact);
        let source = raw.sources.get(&source_identity).ok_or_else(|| {
            refusal(
                "cli-import-source-missing",
                Some(&artifact.role),
                format!("no raw source is bound to project geometry identity `{source_identity}`"),
                "bind the physical file to this exact project geometry row",
            )
        })?;
        let _ = write!(
            out,
            "{{\"source_identity\":{},\"policy\":",
            json_string(&source_identity)
        );
        match &source.policy {
            RawGeometryPolicy::Mesh {
                length_unit,
                max_hole_edges,
                named_groups,
            } => {
                let _ = write!(
                    out,
                    "{{\"kind\":\"mesh\",\"length_unit\":{},\"max_hole_edges\":{},\"named_groups\":",
                    json_string(length_unit),
                    max_hole_edges,
                );
                push_named_groups_json(&mut out, named_groups);
                out.push('}');
            }
            RawGeometryPolicy::FacetedStep {
                root_id,
                length_unit,
                target_h,
                named_groups,
            } => {
                let _ = write!(
                    out,
                    "{{\"kind\":\"faceted-step\",\"root_id\":{},\"length_unit\":{},\"target_h_bits\":{},\"named_groups\":",
                    root_id,
                    json_string(length_unit),
                    json_string(&format!("{:016x}", target_h.to_bits())),
                );
                push_named_groups_json(&mut out, named_groups);
                out.push('}');
            }
        }
        out.push('}');
    }
    out.push_str("]}");
    Ok(out)
}

fn push_named_groups_json(out: &mut String, groups: &[NamedFaceGroup]) {
    out.push('[');
    for (group_index, group) in groups.iter().enumerate() {
        if group_index > 0 {
            out.push(',');
        }
        let _ = write!(out, "{{\"name\":{},\"faces\":[", json_string(&group.name));
        for (face_index, face) in group.faces.iter().enumerate() {
            if face_index > 0 {
                out.push(',');
            }
            let _ = write!(out, "{face}");
        }
        out.push_str("]}");
    }
    out.push(']');
}

fn step_success_receipt_json(decoder_receipt: &str, import_receipt: &str) -> String {
    format!(
        "{{\"schema\":{},\"status\":\"promoted\",\"decoder\":{decoder_receipt},\"tessellation_import\":{import_receipt},\
         \"authority\":\"retained-lower-layer-receipts\",\
         \"no_claim\":\"the wrapper composes retention only; each nested receipt retains its own authority and no-claim boundary\"}}",
        json_string(STEP_IMPORT_RECEIPT_SCHEMA),
    )
}

fn step_failure_receipt_json(decoder_receipt: &str, error: &str) -> String {
    format!(
        "{{\"schema\":{},\"status\":\"refused\",\"decoder\":{decoder_receipt},\
         \"downstream_refusal\":{{\"message\":{}}},\
         \"authority\":\"retained-decoder-receipt-and-diagnostic\",\
         \"no_claim\":\"no promoted tessellation, assignment, SDF certificate, or physical/CAD sameness claim\"}}",
        json_string(STEP_IMPORT_RECEIPT_SCHEMA),
        json_string(error),
    )
}

fn import_record_name(project_hash: &ContentHash, source_identity: &str) -> String {
    format!("{}:{source_identity}", project_hash.to_hex())
}

fn summary_json(
    project_hash: ContentHash,
    retained: &[RetainedGeometryImport],
    assignment_table: &str,
) -> String {
    let mut out = format!(
        "{{\"schema\":{},\"project_hash\":{},\"artifacts\":[",
        json_string(IMPORT_SUMMARY_SCHEMA),
        json_string(&project_hash.to_hex()),
    );
    for (index, artifact) in retained.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"role\":{},\"source_label\":{},\"source_label_authority\":\"caller-reported\",\"source_identity\":{},\"raw_source\":{},\"promotion_receipt\":{},\"promoted_mesh\":{},\"assignment_report\":{},\"import_record\":{}}}",
            json_string(&artifact.role),
            json_string(&artifact.source_label),
            json_string(&artifact.source_identity),
            json_string(&artifact.raw_source.to_hex()),
            json_string(&artifact.promotion_receipt.to_hex()),
            json_string(&artifact.promoted_mesh.to_hex()),
            json_string(&artifact.assignment_report.to_hex()),
            json_string(&artifact.import_record),
        );
    }
    let _ = write!(
        out,
        "],\"assignment_table\":{},\"authority\":\"retained-import-and-assignment-evidence\",\"no_claim\":{}}}",
        json_string(assignment_table),
        json_string(IMPORT_NO_CLAIM),
    );
    out
}

fn refusal_json(error: &GeometryImportRefusal) -> String {
    let role = error
        .role
        .as_ref()
        .map_or_else(|| "null".to_string(), |role| json_string(role));
    format!(
        "{{\"schema\":{},\"code\":{},\"role\":{},\"what\":{},\"fix\":{},\"no_claim\":{}}}",
        json_string(IMPORT_REFUSAL_SCHEMA),
        json_string(error.code),
        role,
        json_string(&error.what),
        json_string(&error.fix),
        json_string(IMPORT_NO_CLAIM),
    )
}

fn json_string(value: &str) -> String {
    let mut out = String::from("\"");
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            character if character <= '\u{1f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(character));
            }
            character => out.push(character),
        }
    }
    out.push('"');
    out
}

fn checkpoint(
    cx: &Cx<'_>,
    code: &'static str,
    stage: &'static str,
) -> Result<(), GeometryImportRefusal> {
    cx.checkpoint().map_err(|_| {
        refusal(
            code,
            None,
            format!("geometry import was cancelled at `{stage}`"),
            "retry under a live deterministic context; no partial success was published",
        )
    })
}

fn project_codec_refusal(error: fs_project::ProjectError) -> GeometryImportRefusal {
    refusal(error.code, None, error.detail, error.hint)
}

fn ledger_refusal(error: LedgerError) -> GeometryImportRefusal {
    refusal(
        "cli-import-ledger",
        None,
        format!("ledger refused geometry import: {error}"),
        "repair the ledger or contention failure, verify ledger lint/integrity, and retry",
    )
}

fn refusal(
    code: &'static str,
    role: Option<&str>,
    what: impl Into<String>,
    fix: impl Into<String>,
) -> GeometryImportRefusal {
    GeometryImportRefusal {
        code,
        role: role.map(str::to_string),
        what: what.into(),
        fix: fix.into(),
        recorded: None,
    }
}
