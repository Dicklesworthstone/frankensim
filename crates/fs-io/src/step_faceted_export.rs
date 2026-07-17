//! Deterministic resource-level re-emission of sealed triangular STEP B-reps.
//!
//! This module accepts only an already decoded bare-`FACE`
//! [`DecodedFacetedBrep`]. It emits a canonical dense-ID
//! `FACETED_BREP -> CLOSED_SHELL -> FACE ->
//! FACE_OUTER_BOUND -> POLY_LOOP -> CARTESIAN_POINT` resource graph, then
//! parses and decodes its own bytes before publication. This is a strict
//! round-trip rung for the native resource subset, not SDF-to-NURBS re-fit,
//! EXPRESS/AP conformance, unit or product context, or topology authority.

use crate::IoError;
use crate::step::{
    StepDocument, StepEntity, StepHeader, StepInstance, StepLimits, StepStructureReceipt,
    StepValue, parse_step_with_limits, write_step_with_limits,
};
use crate::step_faceted::{
    DecodedFacetedBrep, StepFacetedLimits, StepFacetedReceipt, StepFacetedRefusal,
    decode_faceted_brep_with_limits,
};
use fs_exec::Cx;
use fs_rep_mesh::Soup;

/// Versioned meaning of the strict faceted resource re-emitter.
pub const STEP_FACETED_EXPORT_VERSION: &str = "step-triangular-faceted-brep-export-v1";

const CANCELLATION_POLL_STRIDE: usize = 4_096;

/// Caller-supplied deterministic Part-21 `FILE_NAME` metadata.
///
/// Every field is receipt-bound. The syntax writer performs bounded printable
/// ASCII admission; this type does not authenticate any identity or timestamp.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StepFacetedExportMetadata {
    file_name: String,
    timestamp: String,
    author: String,
    organization: String,
    authorization: String,
}

impl StepFacetedExportMetadata {
    /// Retain exact caller metadata for deterministic export.
    #[must_use]
    pub fn new(
        file_name: impl Into<String>,
        timestamp: impl Into<String>,
        author: impl Into<String>,
        organization: impl Into<String>,
        authorization: impl Into<String>,
    ) -> Self {
        Self {
            file_name: file_name.into(),
            timestamp: timestamp.into(),
            author: author.into(),
            organization: organization.into(),
            authorization: authorization.into(),
        }
    }

    /// Exact output file name.
    #[must_use]
    pub fn file_name(&self) -> &str {
        &self.file_name
    }

    /// Exact caller timestamp string.
    #[must_use]
    pub fn timestamp(&self) -> &str {
        &self.timestamp
    }

    /// Exact caller author string.
    #[must_use]
    pub fn author(&self) -> &str {
        &self.author
    }

    /// Exact caller organization string.
    #[must_use]
    pub fn organization(&self) -> &str {
        &self.organization
    }

    /// Exact caller authorization string.
    #[must_use]
    pub fn authorization(&self) -> &str {
        &self.authorization
    }
}

/// Combined syntax and semantic bounds for one strict faceted re-emission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StepFacetedExportLimits {
    /// Part-21 construction, writing, and replay-parse limits.
    pub syntax: StepLimits,
    /// Replay decoder limits.
    pub faceted: StepFacetedLimits,
}

impl Default for StepFacetedExportLimits {
    fn default() -> Self {
        Self {
            syntax: StepLimits::default(),
            faceted: StepFacetedLimits::default(),
        }
    }
}

/// Fail-closed strict faceted resource export refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum StepFacetedExportRefusal {
    /// The sealed source violated an export prerequisite.
    Admission {
        /// Actionable diagnosis.
        what: String,
    },
    /// Checked construction or allocation exceeded a bound.
    Resource {
        /// Stable construction stage.
        stage: &'static str,
        /// Actionable diagnosis.
        what: String,
    },
    /// Caller cancellation was observed before publication.
    Cancelled {
        /// Stable stage that observed cancellation.
        stage: &'static str,
    },
    /// The bounded syntax writer or replay parser refused.
    Syntax {
        /// `write` or `replay-parse`.
        stage: &'static str,
        /// Exact nested syntax refusal.
        source: IoError,
    },
    /// The native semantic replay decoder refused the emitted graph.
    ReplayDecode {
        /// Exact nested decoder refusal.
        source: StepFacetedRefusal,
    },
    /// Replayed coordinates or triangles differed from the sealed source.
    ReplayMismatch {
        /// First deterministic mismatch.
        what: String,
    },
}

impl core::fmt::Display for StepFacetedExportRefusal {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Admission { what } => {
                write!(formatter, "STEP faceted export admission refused: {what}")
            }
            Self::Resource { stage, what } => {
                write!(
                    formatter,
                    "STEP faceted export resource refused at {stage}: {what}"
                )
            }
            Self::Cancelled { stage } => {
                write!(formatter, "STEP faceted export cancelled at {stage}")
            }
            Self::Syntax { stage, source } => {
                write!(formatter, "STEP faceted export {stage} refused: {source}")
            }
            Self::ReplayDecode { source } => {
                write!(
                    formatter,
                    "STEP faceted export replay decode refused: {source}"
                )
            }
            Self::ReplayMismatch { what } => {
                write!(formatter, "STEP faceted export replay mismatch: {what}")
            }
        }
    }
}

impl std::error::Error for StepFacetedExportRefusal {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Syntax { source, .. } => Some(source),
            Self::ReplayDecode { source } => Some(source),
            _ => None,
        }
    }
}

/// Nested source and output receipts for one exact resource replay.
#[derive(Debug, Clone, PartialEq)]
pub struct StepFacetedExportReceipt {
    writer_version: &'static str,
    metadata: StepFacetedExportMetadata,
    limits: StepFacetedExportLimits,
    source_decoder_receipt: StepFacetedReceipt,
    output_structure_receipt: StepStructureReceipt,
    output_decoder_receipt: StepFacetedReceipt,
}

impl StepFacetedExportReceipt {
    /// Resource writer semantics version.
    #[must_use]
    pub const fn writer_version(&self) -> &'static str {
        self.writer_version
    }

    /// Exact caller-supplied header metadata.
    #[must_use]
    pub const fn metadata(&self) -> &StepFacetedExportMetadata {
        &self.metadata
    }

    /// Exact syntax and replay-decoder limits.
    #[must_use]
    pub const fn limits(&self) -> StepFacetedExportLimits {
        self.limits
    }

    /// Complete sealed source decoder receipt.
    #[must_use]
    pub const fn source_decoder_receipt(&self) -> &StepFacetedReceipt {
        &self.source_decoder_receipt
    }

    /// Source-bound syntax receipt from reparsing the emitted bytes.
    #[must_use]
    pub const fn output_structure_receipt(&self) -> &StepStructureReceipt {
        &self.output_structure_receipt
    }

    /// Semantic decoder receipt from replaying the emitted root.
    #[must_use]
    pub const fn output_decoder_receipt(&self) -> &StepFacetedReceipt {
        &self.output_decoder_receipt
    }
}

/// Canonical Part-21 bytes paired with their sealed replay receipt.
#[derive(Debug, PartialEq)]
pub struct ExportedFacetedBrep {
    bytes: Vec<u8>,
    receipt: StepFacetedExportReceipt,
}

impl ExportedFacetedBrep {
    /// Canonical Part-21 bytes.
    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Nested source, syntax, and replay decoder receipt.
    #[must_use]
    pub const fn receipt(&self) -> &StepFacetedExportReceipt {
        &self.receipt
    }

    /// Consume the sealed export into bytes and receipt.
    #[must_use]
    pub fn into_parts(self) -> (Vec<u8>, StepFacetedExportReceipt) {
        (self.bytes, self.receipt)
    }
}

/// Re-emit one sealed strict triangular faceted B-rep with default bounds.
///
/// # Errors
/// [`StepFacetedExportRefusal`] on admission, resource, cancellation, syntax,
/// semantic replay, or bit-exact geometry replay failure.
pub fn export_decoded_faceted_brep(
    source: &DecodedFacetedBrep,
    metadata: StepFacetedExportMetadata,
    cx: &Cx<'_>,
) -> Result<ExportedFacetedBrep, StepFacetedExportRefusal> {
    export_decoded_faceted_brep_with_limits(
        source,
        metadata,
        StepFacetedExportLimits::default(),
        cx,
    )
}

/// Re-emit one sealed strict triangular faceted B-rep with explicit bounds.
///
/// Instance IDs are dense and deterministic: points first, then one
/// `POLY_LOOP`/`FACE_OUTER_BOUND`/`FACE` triple per canonical source triangle,
/// followed by the shell and root. Binary64 coordinates use Rust's stable
/// shortest round-trip decimal spelling. The emitted bytes are reparsed and
/// decoded, and publication refuses unless every coordinate bit and triangle
/// index replays exactly.
///
/// # Errors
/// [`StepFacetedExportRefusal`] on admission, resource, cancellation, syntax,
/// semantic replay, or bit-exact geometry replay failure.
pub fn export_decoded_faceted_brep_with_limits(
    source: &DecodedFacetedBrep,
    metadata: StepFacetedExportMetadata,
    limits: StepFacetedExportLimits,
    cx: &Cx<'_>,
) -> Result<ExportedFacetedBrep, StepFacetedExportRefusal> {
    checkpoint(cx, "entry")?;
    validate_source(source, limits)?;
    let (document, root_id) = build_document(source, &metadata, limits, cx)?;
    checkpoint(cx, "pre-write")?;
    let bytes = write_step_with_limits(&document, limits.syntax).map_err(|source| {
        StepFacetedExportRefusal::Syntax {
            stage: "write",
            source,
        }
    })?;
    checkpoint(cx, "post-write")?;
    let parsed = parse_step_with_limits(&bytes, limits.syntax).map_err(|source| {
        StepFacetedExportRefusal::Syntax {
            stage: "replay-parse",
            source,
        }
    })?;
    checkpoint(cx, "post-replay-parse")?;
    let replay = decode_faceted_brep_with_limits(&parsed, root_id, limits.faceted, cx)
        .map_err(|source| StepFacetedExportRefusal::ReplayDecode { source })?;
    require_exact_replay(source.soup(), replay.soup())?;
    if source.receipt().profile() != replay.receipt().profile() {
        return Err(StepFacetedExportRefusal::ReplayMismatch {
            what: format!(
                "source profile {:?} replayed as {:?}",
                source.receipt().profile(),
                replay.receipt().profile()
            ),
        });
    }
    checkpoint(cx, "publication")?;
    let receipt = StepFacetedExportReceipt {
        writer_version: STEP_FACETED_EXPORT_VERSION,
        metadata,
        limits,
        source_decoder_receipt: source.receipt().clone(),
        output_structure_receipt: parsed.receipt().clone(),
        output_decoder_receipt: replay.receipt().clone(),
    };
    Ok(ExportedFacetedBrep { bytes, receipt })
}

fn validate_source(
    source: &DecodedFacetedBrep,
    limits: StepFacetedExportLimits,
) -> Result<(), StepFacetedExportRefusal> {
    let soup = source.soup();
    if source.receipt().plane_face_count() != 0
        || source.receipt().bare_face_count() != soup.triangles.len()
    {
        return Err(StepFacetedExportRefusal::Admission {
            what: format!(
                "version 1 re-emits only bare FACE resources; source has {} bare and {} plane-backed faces",
                source.receipt().bare_face_count(),
                source.receipt().plane_face_count()
            ),
        });
    }
    if soup.positions.is_empty() || soup.triangles.is_empty() {
        return Err(StepFacetedExportRefusal::Admission {
            what: "sealed source must contain at least one vertex and triangle".to_string(),
        });
    }
    if soup.positions.len() > limits.faceted.max_vertices {
        return Err(StepFacetedExportRefusal::Resource {
            stage: "source-vertices",
            what: format!(
                "source has {} vertices; cap is {}",
                soup.positions.len(),
                limits.faceted.max_vertices
            ),
        });
    }
    if soup.triangles.len() > limits.faceted.max_triangles {
        return Err(StepFacetedExportRefusal::Resource {
            stage: "source-triangles",
            what: format!(
                "source has {} triangles; cap is {}",
                soup.triangles.len(),
                limits.faceted.max_triangles
            ),
        });
    }
    if source_instance_count(soup.positions.len(), soup.triangles.len())?
        > limits.syntax.max_instances
    {
        return Err(StepFacetedExportRefusal::Resource {
            stage: "instance-plan",
            what: format!(
                "dense faceted graph exceeds syntax instance cap {}",
                limits.syntax.max_instances
            ),
        });
    }
    for (index, point) in soup.positions.iter().enumerate() {
        if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
            return Err(StepFacetedExportRefusal::Admission {
                what: format!("source vertex {index} is non-finite"),
            });
        }
    }
    for (face, triangle) in soup.triangles.iter().enumerate() {
        if triangle[0] == triangle[1] || triangle[0] == triangle[2] || triangle[1] == triangle[2] {
            return Err(StepFacetedExportRefusal::Admission {
                what: format!("source face {face} repeats a vertex"),
            });
        }
        for vertex in triangle {
            if usize::try_from(*vertex)
                .ok()
                .is_none_or(|index| index >= soup.positions.len())
            {
                return Err(StepFacetedExportRefusal::Admission {
                    what: format!("source face {face} references missing vertex {vertex}"),
                });
            }
        }
    }
    Ok(())
}

fn source_instance_count(
    vertices: usize,
    triangles: usize,
) -> Result<usize, StepFacetedExportRefusal> {
    triangles
        .checked_mul(3)
        .and_then(|face_instances| vertices.checked_add(face_instances))
        .and_then(|instances| instances.checked_add(2))
        .ok_or_else(|| StepFacetedExportRefusal::Resource {
            stage: "instance-plan",
            what: "dense faceted instance count overflowed".to_string(),
        })
}

#[allow(clippy::too_many_lines)] // One fixed resource graph plus its checked dense-ID plan.
fn build_document(
    source: &DecodedFacetedBrep,
    metadata: &StepFacetedExportMetadata,
    limits: StepFacetedExportLimits,
    cx: &Cx<'_>,
) -> Result<(StepDocument, u64), StepFacetedExportRefusal> {
    let soup = source.soup();
    let instance_count = source_instance_count(soup.positions.len(), soup.triangles.len())?;
    let mut instances = Vec::new();
    instances
        .try_reserve_exact(instance_count)
        .map_err(|error| StepFacetedExportRefusal::Resource {
            stage: "instances",
            what: error.to_string(),
        })?;

    for (index, point) in soup.positions.iter().enumerate() {
        poll(cx, "point-instances", index)?;
        let id = dense_id(index, 1, "point-id")?;
        instances.push(simple_instance(
            id,
            "CARTESIAN_POINT",
            vec![
                StepValue::String(String::new()),
                StepValue::Aggregate(vec![
                    StepValue::Number(step_number(point.x)),
                    StepValue::Number(step_number(point.y)),
                    StepValue::Number(step_number(point.z)),
                ]),
            ],
        ));
    }

    let face_base = soup.positions.len();
    let mut face_references = Vec::new();
    face_references
        .try_reserve_exact(soup.triangles.len())
        .map_err(|error| StepFacetedExportRefusal::Resource {
            stage: "face-references",
            what: error.to_string(),
        })?;
    for (face, triangle) in soup.triangles.iter().enumerate() {
        poll(cx, "face-instances", face)?;
        let offset = face
            .checked_mul(3)
            .and_then(|offset| face_base.checked_add(offset))
            .ok_or_else(|| StepFacetedExportRefusal::Resource {
                stage: "face-id-plan",
                what: "dense face ID offset overflowed".to_string(),
            })?;
        let loop_id = dense_id(offset, 1, "poly-loop-id")?;
        let bound_id = dense_id(offset, 2, "outer-bound-id")?;
        let face_id = dense_id(offset, 3, "face-id")?;
        let point_references = vec![
            StepValue::Reference(u64::from(triangle[0]) + 1),
            StepValue::Reference(u64::from(triangle[1]) + 1),
            StepValue::Reference(u64::from(triangle[2]) + 1),
        ];
        instances.push(simple_instance(
            loop_id,
            "POLY_LOOP",
            vec![
                StepValue::String(String::new()),
                StepValue::Aggregate(point_references),
            ],
        ));
        instances.push(simple_instance(
            bound_id,
            "FACE_OUTER_BOUND",
            vec![
                StepValue::String(String::new()),
                StepValue::Reference(loop_id),
                StepValue::Enumeration("T".to_string()),
            ],
        ));
        instances.push(simple_instance(
            face_id,
            "FACE",
            vec![
                StepValue::String(String::new()),
                StepValue::Aggregate(vec![StepValue::Reference(bound_id)]),
            ],
        ));
        face_references.push(StepValue::Reference(face_id));
    }

    let shell_offset =
        instance_count
            .checked_sub(2)
            .ok_or_else(|| StepFacetedExportRefusal::Resource {
                stage: "root-id-plan",
                what: "dense shell offset underflowed".to_string(),
            })?;
    let shell_id = dense_id(shell_offset, 1, "shell-id")?;
    let root_id = dense_id(shell_offset, 2, "root-id")?;
    instances.push(simple_instance(
        shell_id,
        "CLOSED_SHELL",
        vec![
            StepValue::String(String::new()),
            StepValue::Aggregate(face_references),
        ],
    ));
    instances.push(simple_instance(
        root_id,
        "FACETED_BREP",
        vec![
            StepValue::String(String::new()),
            StepValue::Reference(shell_id),
        ],
    ));
    checkpoint(cx, "document-publication")?;

    let header = StepHeader {
        file_description: vec![
            StepValue::Aggregate(vec![StepValue::String(
                "FrankenSim strict triangular faceted B-rep resource".to_string(),
            )]),
            StepValue::String("2;1".to_string()),
        ],
        file_name: vec![
            StepValue::String(metadata.file_name.clone()),
            StepValue::String(metadata.timestamp.clone()),
            StepValue::Aggregate(vec![StepValue::String(metadata.author.clone())]),
            StepValue::Aggregate(vec![StepValue::String(metadata.organization.clone())]),
            StepValue::String(format!("fs-io@{}", crate::VERSION)),
            StepValue::String("FrankenSim".to_string()),
            StepValue::String(metadata.authorization.clone()),
        ],
        file_schema: vec![StepValue::Aggregate(vec![StepValue::String(
            source.receipt().profile().schema_identifier().to_string(),
        )])],
    };
    let document = StepDocument { header, instances };
    if document.instances.len() != instance_count
        || document.instances.len() > limits.syntax.max_instances
    {
        return Err(StepFacetedExportRefusal::Resource {
            stage: "instance-publication",
            what: format!(
                "constructed {} instances from a plan of {} under cap {}",
                document.instances.len(),
                instance_count,
                limits.syntax.max_instances
            ),
        });
    }
    Ok((document, root_id))
}

fn dense_id(
    zero_based: usize,
    add: usize,
    stage: &'static str,
) -> Result<u64, StepFacetedExportRefusal> {
    zero_based
        .checked_add(add)
        .and_then(|value| u64::try_from(value).ok())
        .ok_or_else(|| StepFacetedExportRefusal::Resource {
            stage,
            what: "dense Part-21 instance ID overflowed".to_string(),
        })
}

fn simple_instance(id: u64, name: &str, parameters: Vec<StepValue>) -> StepInstance {
    StepInstance {
        id,
        components: vec![StepEntity {
            name: name.to_string(),
            parameters,
        }],
    }
}

fn step_number(value: f64) -> String {
    let mut token = value.to_string();
    if let Some(mut exponent) = token.find('e') {
        if !token[..exponent].contains('.') {
            token.insert_str(exponent, ".0");
            exponent += 2;
        }
        token.replace_range(exponent..exponent + 1, "E");
    }
    token
}

fn require_exact_replay(source: &Soup, replay: &Soup) -> Result<(), StepFacetedExportRefusal> {
    if source.positions.len() != replay.positions.len() {
        return Err(StepFacetedExportRefusal::ReplayMismatch {
            what: format!(
                "vertex count changed from {} to {}",
                source.positions.len(),
                replay.positions.len()
            ),
        });
    }
    if source.triangles.len() != replay.triangles.len() {
        return Err(StepFacetedExportRefusal::ReplayMismatch {
            what: format!(
                "triangle count changed from {} to {}",
                source.triangles.len(),
                replay.triangles.len()
            ),
        });
    }
    for (index, (expected, found)) in source.positions.iter().zip(&replay.positions).enumerate() {
        let expected_bits = [
            expected.x.to_bits(),
            expected.y.to_bits(),
            expected.z.to_bits(),
        ];
        let found_bits = [found.x.to_bits(), found.y.to_bits(), found.z.to_bits()];
        if expected_bits != found_bits {
            return Err(StepFacetedExportRefusal::ReplayMismatch {
                what: format!(
                    "vertex {index} bits changed from {expected_bits:x?} to {found_bits:x?}"
                ),
            });
        }
    }
    for (face, (expected, found)) in source.triangles.iter().zip(&replay.triangles).enumerate() {
        if expected != found {
            return Err(StepFacetedExportRefusal::ReplayMismatch {
                what: format!("triangle {face} changed from {expected:?} to {found:?}"),
            });
        }
    }
    Ok(())
}

fn checkpoint(cx: &Cx<'_>, stage: &'static str) -> Result<(), StepFacetedExportRefusal> {
    cx.checkpoint()
        .map_err(|_| StepFacetedExportRefusal::Cancelled { stage })
}

fn poll(cx: &Cx<'_>, stage: &'static str, offset: usize) -> Result<(), StepFacetedExportRefusal> {
    if offset % CANCELLATION_POLL_STRIDE == 0 {
        checkpoint(cx, stage)?;
    }
    Ok(())
}
