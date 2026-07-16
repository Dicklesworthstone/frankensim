//! Caller-tessellated STEP handoff into the existing mesh-to-SDF lane.
//!
//! The Part-21 syntax kernel does not interpret EXPRESS geometry. This module
//! accepts a tessellation produced by an explicitly identified external or
//! future in-house adapter, removes duplicate/degenerate faces and unreferenced
//! vertices, repairs orientation, then refuses leaks/non-manifoldness rather
//! than filling them. A successful result is still only an **estimated** SDF: component
//! nesting, self-intersection freedom, and generalized-winding sign remain
//! uncertified.

use crate::step::{ParsedStep, STEP_SYNTAX_VERSION};
use fs_evidence::{Evidence, NumericalCertificate, NumericalKind, ProvenanceHash, vv::UnitId};
use fs_exec::{Cx, ExecMode};
use fs_obs::ident::IdentityBuilder;
use fs_rep_mesh::{MeshChart, MeshQuality, MeshSdfError, RepairReceipt, Soup, mesh_to_sdf, repair};
use fs_rep_sdf::TiledSdf;
use std::fmt::Write as _;

/// Maximum vertices admitted in one caller-supplied STEP tessellation.
pub const MAX_STEP_TESSELLATION_VERTICES: usize = 1_000_000;
/// Maximum triangles admitted in one caller-supplied STEP tessellation.
pub const MAX_STEP_TESSELLATION_TRIANGLES: usize = 1_000_000;
/// Maximum localized mesh defects retained in a refusal.
pub const MAX_STEP_LOCALIZED_DEFECTS: usize = 256;
/// Maximum bytes in either tessellator identity string.
pub const MAX_STEP_ADAPTER_ID_BYTES: usize = 256;
/// Conservative auxiliary-memory admission ceiling for topology preprocessing.
pub const MAX_STEP_TOPOLOGY_AUXILIARY_BYTES: usize = 512 * 1024 * 1024;
const ESTIMATED_AUX_BYTES_PER_VERTEX: usize = 64;
const ESTIMATED_AUX_BYTES_PER_TRIANGLE: usize = 384;
const CANCELLATION_POLL_STRIDE: usize = 4_096;
/// Versioned semantics label for the complete STEP-tessellation handoff.
pub const STEP_IMPORT_SEMANTICS_VERSION: &str = "step-tessellation-to-sdf-v1";
/// Domain of the bit-exact, non-cryptographic tessellation fingerprint.
pub const STEP_TESSELLATION_FINGERPRINT_DOMAIN: &str = "org.frankensim.fs-io.step-tessellation.v1";

/// Identity of the adapter that materialized a STEP tessellation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StepTessellatorIdentity {
    /// Stable adapter/tool name.
    pub name: String,
    /// Stable adapter/tool version.
    pub version: String,
    /// Non-cryptographic fingerprint of all tessellation settings.
    pub configuration_fingerprint: u64,
}

/// A localized mesh-integrity defect after non-hole-filling repair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepMeshDefectKind {
    /// An edge has only one incident face.
    BoundaryEdge,
    /// An edge has more than two incident faces.
    NonManifoldEdge,
    /// Two incident faces traverse their common edge in the same direction.
    OrientationConflict,
    /// A vertex-link/half-edge manifold check refused the mesh.
    VertexLinkNonManifold,
    /// Aggregate signed orientation is not finite and outward-positive.
    AggregateOrientation,
}

impl StepMeshDefectKind {
    /// Stable receipt label for this defect class.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::BoundaryEdge => "boundary-edge",
            Self::NonManifoldEdge => "non-manifold-edge",
            Self::OrientationConflict => "orientation-conflict",
            Self::VertexLinkNonManifold => "vertex-link-non-manifold",
            Self::AggregateOrientation => "aggregate-orientation",
        }
    }
}

/// Deterministic location for one rejected integrity condition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StepMeshDefect {
    /// Defect class.
    pub kind: StepMeshDefectKind,
    /// Undirected vertex pair when the defect is edge-local.
    pub edge: Option<[u32; 2]>,
    /// First three incident face ordinals, in repaired-soup traversal order.
    pub incident_faces: [Option<usize>; 3],
    /// Vertex ordinal when the defect is vertex-local.
    pub vertex: Option<u32>,
    /// Total number of incident faces for an edge- or vertex-local defect.
    pub total_incidents: u32,
    /// Additional deterministic diagnosis.
    pub detail: String,
}

/// Fail-closed result before an estimated STEP-derived SDF can be published.
#[derive(Debug, Clone, PartialEq)]
pub enum StepImportRefusal {
    /// Request or raw tessellation admission failed before repair.
    Admission {
        /// Source Part-21 fingerprint.
        source_fingerprint: u64,
        /// Actionable diagnosis.
        what: String,
    },
    /// Repaired tessellation is still leaky, non-manifold, or misoriented.
    MeshIntegrity {
        /// Source Part-21 fingerprint.
        source_fingerprint: u64,
        /// Aggregate quality counters.
        quality: MeshQuality,
        /// Deterministic bounded prefix of localized defects.
        defects: Vec<StepMeshDefect>,
        /// Whether more localized defects existed than the receipt retained.
        localized_truncated: bool,
        /// Duplicate/degenerate/orientation repair actions already attempted.
        repairs: Vec<RepairReceipt>,
    },
    /// Mesh-to-SDF construction refused or preprocessing observed cancellation.
    SdfBuild {
        /// Source Part-21 fingerprint.
        source_fingerprint: u64,
        /// Underlying structured sampler refusal.
        error: MeshSdfError,
        /// Repair actions preceding the sampler.
        repairs: Vec<RepairReceipt>,
    },
    /// The mesh-to-SDF evidence could not be honestly composed.
    Evidence {
        /// Source Part-21 fingerprint.
        source_fingerprint: u64,
        /// Actionable diagnosis.
        what: String,
        /// Repair actions preceding evidence composition.
        repairs: Vec<RepairReceipt>,
    },
    /// A bounded preprocessing allocation could not be admitted.
    Resource {
        /// Source Part-21 fingerprint.
        source_fingerprint: u64,
        /// Preprocessing stage that needed the allocation.
        stage: &'static str,
        /// Actionable diagnosis.
        what: String,
        /// Repair actions preceding the refusal.
        repairs: Vec<RepairReceipt>,
    },
}

impl core::fmt::Display for StepImportRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Admission { what, .. } => {
                write!(f, "STEP tessellation admission refused: {what}")
            }
            Self::MeshIntegrity {
                quality, defects, ..
            } => write!(
                f,
                "STEP tessellation integrity refused: {} boundary, {} non-manifold, {} orientation conflicts, outward={}; {} localized diagnostics",
                quality.boundary_edges,
                quality.nonmanifold_edges,
                quality.orientation_conflicts,
                quality.outward_oriented,
                defects.len()
            ),
            Self::SdfBuild { error, .. } => {
                write!(f, "STEP tessellation SDF construction refused: {error}")
            }
            Self::Evidence { what, .. } => {
                write!(f, "STEP tessellation evidence refused: {what}")
            }
            Self::Resource { stage, what, .. } => {
                write!(
                    f,
                    "STEP tessellation resource admission refused at {stage}: {what}"
                )
            }
        }
    }
}

impl std::error::Error for StepImportRefusal {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::SdfBuild { error, .. } => Some(error),
            _ => None,
        }
    }
}

/// Source-bound receipt for the caller-tessellated estimated-SDF handoff.
#[derive(Debug, Clone, PartialEq)]
pub struct StepImportReceipt {
    source_fingerprint: u64,
    canonical_layout_fingerprint: u64,
    source_tessellation_fingerprint: u64,
    repaired_tessellation_fingerprint: u64,
    schemas: Vec<String>,
    tessellator: StepTessellatorIdentity,
    length_unit: UnitId,
    tessellation_deviation: NumericalCertificate,
    source_vertices: usize,
    source_triangles: usize,
    repaired_vertices: usize,
    repaired_triangles: usize,
    repairs: Vec<RepairReceipt>,
    quality: MeshQuality,
    target_h: f64,
    mesh_sdf_numerical: NumericalCertificate,
    combined_numerical: NumericalCertificate,
    output_provenance: ProvenanceHash,
}

impl StepImportReceipt {
    /// Exact-source non-cryptographic fingerprint from the syntax receipt.
    #[must_use]
    pub const fn source_fingerprint(&self) -> u64 {
        self.source_fingerprint
    }

    /// Layout-canonical non-cryptographic fingerprint from the syntax receipt.
    #[must_use]
    pub const fn canonical_layout_fingerprint(&self) -> u64 {
        self.canonical_layout_fingerprint
    }

    /// Bit-exact non-cryptographic fingerprint of the caller's triangle soup.
    #[must_use]
    pub const fn source_tessellation_fingerprint(&self) -> u64 {
        self.source_tessellation_fingerprint
    }

    /// Bit-exact non-cryptographic fingerprint after recorded safe repairs.
    #[must_use]
    pub const fn repaired_tessellation_fingerprint(&self) -> u64 {
        self.repaired_tessellation_fingerprint
    }

    /// Declared schemas copied from the sealed syntax receipt.
    #[must_use]
    pub fn schema_identifiers(&self) -> &[String] {
        &self.schemas
    }

    /// Caller-declared tessellator identity retained by the receipt.
    #[must_use]
    pub const fn tessellator(&self) -> &StepTessellatorIdentity {
        &self.tessellator
    }

    /// Unit shared by tessellation coordinates, deviation, and sample spacing.
    #[must_use]
    pub const fn length_unit(&self) -> &UnitId {
        &self.length_unit
    }

    /// Caller-declared tessellation deviation band.
    #[must_use]
    pub const fn tessellation_deviation(&self) -> NumericalCertificate {
        self.tessellation_deviation
    }

    /// Total estimated SDF error band after deviation composition.
    #[must_use]
    pub const fn combined_numerical(&self) -> NumericalCertificate {
        self.combined_numerical
    }

    /// Mesh-to-SDF numerical evidence before caller-deviation composition.
    #[must_use]
    pub const fn mesh_sdf_numerical(&self) -> NumericalCertificate {
        self.mesh_sdf_numerical
    }

    /// Requested SDF sample spacing.
    #[must_use]
    pub const fn target_h(&self) -> f64 {
        self.target_h
    }

    /// Provenance fingerprint bound to the published evidence.
    #[must_use]
    pub const fn output_provenance(&self) -> ProvenanceHash {
        self.output_provenance
    }

    /// Basic edge-use and aggregate-orientation diagnostics.
    #[must_use]
    pub const fn quality(&self) -> MeshQuality {
        self.quality
    }

    /// Deterministic repair actions applied before integrity admission.
    #[must_use]
    pub fn repairs(&self) -> &[RepairReceipt] {
        &self.repairs
    }

    /// Canonical JSON for a syntax-to-estimated-SDF ledger event.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut out = String::from(
            "{\"kind\":\"step-tessellation-to-sdf-receipt\",\"authority\":\"estimate\",",
        );
        let _ = write!(
            out,
            "\"syntax_version\":\"{}\",\"fs_io_version\":\"{}\",\
             \"fs_rep_mesh_version\":\"{}\",\"fs_rep_sdf_version\":\"{}\",\
             \"step_import_semantics\":\"{}\",\
             \"tessellation_fingerprint_domain\":\"{}\",\
             \"execution_mode\":\"deterministic\",\
             \"source_fingerprint_fnv1a64\":\"{:016x}\",\
             \"canonical_layout_fingerprint_fnv1a64\":\"{:016x}\",\
             \"source_tessellation_fingerprint_fnv1a64\":\"{:016x}\",\
             \"repaired_tessellation_fingerprint_fnv1a64\":\"{:016x}\",\"schemas\":[",
            STEP_SYNTAX_VERSION,
            crate::VERSION,
            fs_rep_mesh::VERSION,
            fs_rep_sdf::VERSION,
            STEP_IMPORT_SEMANTICS_VERSION,
            STEP_TESSELLATION_FINGERPRINT_DOMAIN,
            self.source_fingerprint,
            self.canonical_layout_fingerprint,
            self.source_tessellation_fingerprint,
            self.repaired_tessellation_fingerprint
        );
        for (index, schema) in self.schemas.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            push_json_string(&mut out, schema);
        }
        out.push_str("],\"tessellator\":{\"name\":");
        push_json_string(&mut out, &self.tessellator.name);
        out.push_str(",\"version\":");
        push_json_string(&mut out, &self.tessellator.version);
        let _ = write!(
            out,
            ",\"configuration_fingerprint_noncrypto64\":\"{:016x}\"}},\"length_unit\":",
            self.tessellator.configuration_fingerprint
        );
        push_json_string(&mut out, self.length_unit.as_str());
        out.push_str(",\"tessellation_deviation\":");
        push_numerical_json(&mut out, self.tessellation_deviation);
        let _ = write!(
            out,
            ",\"source_mesh\":{{\"vertices\":{},\"triangles\":{}}},\
             \"repaired_mesh\":{{\"vertices\":{},\"triangles\":{}}},\
             \"admission_limits\":{{\"max_vertices\":{},\"max_triangles\":{},\
             \"max_topology_auxiliary_bytes\":{}}},\"repairs\":[",
            self.source_vertices,
            self.source_triangles,
            self.repaired_vertices,
            self.repaired_triangles,
            MAX_STEP_TESSELLATION_VERTICES,
            MAX_STEP_TESSELLATION_TRIANGLES,
            MAX_STEP_TOPOLOGY_AUXILIARY_BYTES
        );
        for (index, repair) in self.repairs.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            out.push_str(&repair.to_json());
        }
        let _ = write!(
            out,
            "],\"quality\":{{\"invalid_triangles\":{},\"boundary_edges\":{},\
             \"nonmanifold_edges\":{},\"orientation_conflicts\":{},\
             \"aggregate_outward\":{}}},\"target_h\":{},\"mesh_sdf_numerical\":",
            self.quality.invalid_triangles,
            self.quality.boundary_edges,
            self.quality.nonmanifold_edges,
            self.quality.orientation_conflicts,
            self.quality.outward_oriented,
            self.target_h
        );
        push_numerical_json(&mut out, self.mesh_sdf_numerical);
        out.push_str(",\"combined_numerical\":");
        push_numerical_json(&mut out, self.combined_numerical);
        let _ = write!(
            out,
            ",\"output_provenance_fingerprint\":\"{:016x}\",\
             \"sign_confidence\":\"uncertified\",\
             \"no_claim\":\"caller-supplied tessellation; no EXPRESS interpretation, component-nesting proof, self-intersection certificate, or globally certified winding sign\"}}",
            self.output_provenance.0
        );
        out
    }
}

/// Estimated SDF paired immutably with its STEP-tessellation receipt.
#[derive(Debug)]
pub struct StepImportOutcome {
    evidence: Evidence<TiledSdf>,
    receipt: StepImportReceipt,
}

impl StepImportOutcome {
    /// Read-only SDF evidence.
    #[must_use]
    pub const fn evidence(&self) -> &Evidence<TiledSdf> {
        &self.evidence
    }

    /// Read-only source-bound import receipt.
    #[must_use]
    pub const fn receipt(&self) -> &StepImportReceipt {
        &self.receipt
    }
}

/// Convert an explicitly caller-tessellated STEP document into an estimated
/// tiled SDF.
///
/// Duplicate/degenerate faces and unreferenced vertices may be removed, and
/// orientation may be unified. Holes are never filled: any residual boundary, non-manifold edge,
/// vertex-link failure, orientation conflict, or non-outward aggregate refuses
/// publication with localized diagnostics.
/// `length_unit` applies jointly to soup coordinates, the deviation band, and
/// `target_h`; the API therefore cannot compose differently declared units.
///
/// # Errors
/// [`StepImportRefusal`] for invalid adapter/deviation/mesh inputs, residual
/// integrity defects, cancellation, sampler refusal, or non-finite evidence
/// composition.
#[allow(clippy::too_many_lines)] // One ordered quarantine, topology, evidence, and receipt transaction.
pub fn import_step_tessellation(
    parsed: &ParsedStep,
    soup: Soup,
    tessellator: StepTessellatorIdentity,
    tessellation_deviation: NumericalCertificate,
    length_unit: UnitId,
    target_h: f64,
    cx: &Cx<'_>,
) -> Result<StepImportOutcome, StepImportRefusal> {
    let source_fingerprint = parsed.receipt().source_fingerprint();
    observe_cancellation(cx, source_fingerprint, &[])?;
    if cx.mode() != ExecMode::Deterministic {
        return Err(StepImportRefusal::Admission {
            source_fingerprint,
            what: format!(
                "STEP tessellation publication requires deterministic execution mode; got {}",
                cx.mode().name()
            ),
        });
    }
    validate_adapter(&tessellator).map_err(|what| StepImportRefusal::Admission {
        source_fingerprint,
        what,
    })?;
    validate_deviation(tessellation_deviation).map_err(|what| StepImportRefusal::Admission {
        source_fingerprint,
        what,
    })?;
    if !target_h.is_finite() || target_h <= 0.0 {
        return Err(StepImportRefusal::Admission {
            source_fingerprint,
            what: format!("target_h must be finite and positive; got {target_h}"),
        });
    }
    match validate_soup(&soup, cx) {
        Ok(()) => {}
        Err(SoupValidationError::Admission(what)) => {
            return Err(StepImportRefusal::Admission {
                source_fingerprint,
                what,
            });
        }
        Err(SoupValidationError::Cancelled) => {
            return Err(StepImportRefusal::SdfBuild {
                source_fingerprint,
                error: MeshSdfError::Cancelled,
                repairs: Vec::new(),
            });
        }
    }

    let source_vertices = soup.positions.len();
    let source_triangles = soup.triangles.len();
    let source_tessellation_fingerprint =
        tessellation_fingerprint(&soup, cx).map_err(|error| StepImportRefusal::SdfBuild {
            source_fingerprint,
            error,
            repairs: Vec::new(),
        })?;
    let mut repaired = repair(soup, 0);
    observe_cancellation(cx, source_fingerprint, &repaired.receipts)?;
    match compact_unreferenced_vertices(&mut repaired.soup, &mut repaired.receipts, cx) {
        Ok(()) => {}
        Err(PreprocessError::Cancelled) => {
            return Err(StepImportRefusal::SdfBuild {
                source_fingerprint,
                error: MeshSdfError::Cancelled,
                repairs: repaired.receipts,
            });
        }
        Err(PreprocessError::Resource { stage, what }) => {
            return Err(StepImportRefusal::Resource {
                source_fingerprint,
                stage,
                what,
                repairs: repaired.receipts,
            });
        }
    }
    observe_cancellation(cx, source_fingerprint, &repaired.receipts)?;
    let repaired_vertices = repaired.soup.positions.len();
    let repaired_triangles = repaired.soup.triangles.len();
    let repaired_tessellation_fingerprint =
        tessellation_fingerprint(&repaired.soup, cx).map_err(|error| {
            StepImportRefusal::SdfBuild {
                source_fingerprint,
                error,
                repairs: repaired.receipts.clone(),
            }
        })?;

    let (quality, mut defects, mut localized_truncated) =
        match quality_and_localized_defects(&repaired.soup, cx) {
            Ok(result) => result,
            Err(PreprocessError::Cancelled) => {
                return Err(StepImportRefusal::SdfBuild {
                    source_fingerprint,
                    error: MeshSdfError::Cancelled,
                    repairs: repaired.receipts,
                });
            }
            Err(PreprocessError::Resource { stage, what }) => {
                return Err(StepImportRefusal::Resource {
                    source_fingerprint,
                    stage,
                    what,
                    repairs: repaired.receipts,
                });
            }
        };
    let (vertex_defects, vertex_truncated) = match disconnected_vertex_links(&repaired.soup, cx) {
        Ok(result) => result,
        Err(PreprocessError::Cancelled) => {
            return Err(StepImportRefusal::SdfBuild {
                source_fingerprint,
                error: MeshSdfError::Cancelled,
                repairs: repaired.receipts,
            });
        }
        Err(PreprocessError::Resource { stage, what }) => {
            return Err(StepImportRefusal::Resource {
                source_fingerprint,
                stage,
                what,
                repairs: repaired.receipts,
            });
        }
    };
    let has_vertex_defects = !vertex_defects.is_empty();
    for defect in vertex_defects {
        push_defect_bounded(&mut defects, defect, &mut localized_truncated);
    }
    localized_truncated |= vertex_truncated;
    if repaired_triangles == 0 || !quality.passes_basic_orientation_checks() || has_vertex_defects {
        return Err(StepImportRefusal::MeshIntegrity {
            source_fingerprint,
            quality,
            defects,
            localized_truncated,
            repairs: repaired.receipts,
        });
    }

    let chart = MeshChart::new(repaired.soup);
    let mut evidence =
        mesh_to_sdf(&chart, target_h, cx).map_err(|error| StepImportRefusal::SdfBuild {
            source_fingerprint,
            error,
            repairs: repaired.receipts.clone(),
        })?;
    let mesh_sdf_numerical = evidence.numerical;
    let combined_numerical = compose_estimated_error(mesh_sdf_numerical, tessellation_deviation)
        .map_err(|what| StepImportRefusal::Evidence {
            source_fingerprint,
            what,
            repairs: repaired.receipts.clone(),
        })?;

    let adapter_identity = IdentityBuilder::new("fs-io-step-tessellator-v1")
        .str("name", &tessellator.name)
        .str("version", &tessellator.version)
        .u64(
            "configuration-fingerprint",
            tessellator.configuration_fingerprint,
        )
        .finish();
    let parameter_identity = IdentityBuilder::new("fs-io-step-import-parameters-v1")
        .str("step-import-semantics", STEP_IMPORT_SEMANTICS_VERSION)
        .str("fs-io-version", crate::VERSION)
        .str("fs-rep-mesh-version", fs_rep_mesh::VERSION)
        .str("fs-rep-sdf-version", fs_rep_sdf::VERSION)
        .str("length-unit", length_unit.as_str())
        .str("execution-mode", cx.mode().name())
        .u64("target-h-bits", target_h.to_bits())
        .str(
            "tessellation-deviation-kind",
            numerical_kind_label(tessellation_deviation.kind),
        )
        .u64(
            "tessellation-deviation-lo-bits",
            tessellation_deviation.lo.to_bits(),
        )
        .u64(
            "tessellation-deviation-hi-bits",
            tessellation_deviation.hi.to_bits(),
        )
        .finish();
    evidence.provenance = ProvenanceHash::chain(
        STEP_IMPORT_SEMANTICS_VERSION,
        &[
            ProvenanceHash(source_fingerprint),
            ProvenanceHash(parsed.receipt().canonical_layout_fingerprint()),
            ProvenanceHash(source_tessellation_fingerprint),
            ProvenanceHash(repaired_tessellation_fingerprint),
            ProvenanceHash(adapter_identity.root()),
            ProvenanceHash(parameter_identity.root()),
            evidence.provenance,
        ],
    );
    evidence.qoi = combined_numerical.hi;
    evidence.numerical = combined_numerical;
    evidence
        .model
        .cards
        .push("step-caller-tessellation".to_string());
    evidence.model.assumptions.push(format!(
        "STEP geometry was materialized by caller-identified tessellator {}@{} in unit {} with a caller-supplied {:?} deviation upper bound {}; component nesting, self-intersection freedom, semantic correspondence to Part-21 entities, and generalized-winding sign remain uncertified",
        tessellator.name,
        tessellator.version,
        length_unit,
        tessellation_deviation.kind,
        tessellation_deviation.hi
    ));

    let schemas = clone_schemas(parsed.receipt().schema_identifiers()).map_err(|what| {
        StepImportRefusal::Evidence {
            source_fingerprint,
            what,
            repairs: repaired.receipts.clone(),
        }
    })?;
    let receipt = StepImportReceipt {
        source_fingerprint,
        canonical_layout_fingerprint: parsed.receipt().canonical_layout_fingerprint(),
        source_tessellation_fingerprint,
        repaired_tessellation_fingerprint,
        schemas,
        tessellator,
        length_unit,
        tessellation_deviation,
        source_vertices,
        source_triangles,
        repaired_vertices,
        repaired_triangles,
        repairs: repaired.receipts,
        quality,
        target_h,
        mesh_sdf_numerical,
        combined_numerical,
        output_provenance: evidence.provenance,
    };
    Ok(StepImportOutcome { evidence, receipt })
}

fn validate_adapter(identity: &StepTessellatorIdentity) -> Result<(), String> {
    for (field, value) in [("name", &identity.name), ("version", &identity.version)] {
        if value.is_empty() || value.len() > MAX_STEP_ADAPTER_ID_BYTES {
            return Err(format!(
                "tessellator {field} length {} is outside 1..={MAX_STEP_ADAPTER_ID_BYTES}",
                value.len()
            ));
        }
        if value.trim().is_empty()
            || !value.is_ascii()
            || value
                .bytes()
                .any(|byte| byte.is_ascii_control() || byte == b'\0')
        {
            return Err(format!(
                "tessellator {field} must be non-blank, printable, non-NUL ASCII"
            ));
        }
    }
    Ok(())
}

fn validate_deviation(certificate: NumericalCertificate) -> Result<(), String> {
    if certificate.kind == NumericalKind::NoClaim {
        return Err(
            "tessellation deviation must carry Exact, Enclosure, or Estimate authority".to_string(),
        );
    }
    if !certificate.lo.is_finite()
        || !certificate.hi.is_finite()
        || certificate.lo < 0.0
        || certificate.hi < certificate.lo
    {
        return Err(format!(
            "tessellation deviation must be a finite ordered non-negative band; got {:?} [{}, {}]",
            certificate.kind, certificate.lo, certificate.hi
        ));
    }
    if certificate.kind == NumericalKind::Exact && certificate.lo != certificate.hi {
        return Err(format!(
            "an Exact tessellation deviation must be a singleton; got [{}, {}]",
            certificate.lo, certificate.hi
        ));
    }
    Ok(())
}

enum SoupValidationError {
    Admission(String),
    Cancelled,
}

fn validate_soup(soup: &Soup, cx: &Cx<'_>) -> Result<(), SoupValidationError> {
    if soup.positions.is_empty() {
        return Err(SoupValidationError::Admission(
            "STEP tessellation has no vertices".to_string(),
        ));
    }
    if soup.triangles.is_empty() {
        return Err(SoupValidationError::Admission(
            "STEP tessellation has no triangles".to_string(),
        ));
    }
    if soup.positions.len() > MAX_STEP_TESSELLATION_VERTICES {
        return Err(SoupValidationError::Admission(format!(
            "STEP tessellation has {} vertices; cap is {MAX_STEP_TESSELLATION_VERTICES}",
            soup.positions.len()
        )));
    }
    if soup.triangles.len() > MAX_STEP_TESSELLATION_TRIANGLES {
        return Err(SoupValidationError::Admission(format!(
            "STEP tessellation has {} triangles; cap is {MAX_STEP_TESSELLATION_TRIANGLES}",
            soup.triangles.len()
        )));
    }
    let estimated_auxiliary_bytes = soup
        .positions
        .len()
        .checked_mul(ESTIMATED_AUX_BYTES_PER_VERTEX)
        .and_then(|vertices| {
            soup.triangles
                .len()
                .checked_mul(ESTIMATED_AUX_BYTES_PER_TRIANGLE)
                .and_then(|triangles| vertices.checked_add(triangles))
        })
        .ok_or_else(|| {
            SoupValidationError::Admission(
                "STEP tessellation auxiliary-memory estimate overflowed".to_string(),
            )
        })?;
    if estimated_auxiliary_bytes > MAX_STEP_TOPOLOGY_AUXILIARY_BYTES {
        return Err(SoupValidationError::Admission(format!(
            "STEP tessellation needs an estimated {estimated_auxiliary_bytes} auxiliary bytes; cap is {MAX_STEP_TOPOLOGY_AUXILIARY_BYTES}"
        )));
    }
    for (index, point) in soup.positions.iter().enumerate() {
        if index % CANCELLATION_POLL_STRIDE == 0 && cx.checkpoint().is_err() {
            return Err(SoupValidationError::Cancelled);
        }
        if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
            return Err(SoupValidationError::Admission(format!(
                "STEP tessellation vertex {index} has non-finite coordinates"
            )));
        }
    }
    for (face, triangle) in soup.triangles.iter().enumerate() {
        if face % CANCELLATION_POLL_STRIDE == 0 && cx.checkpoint().is_err() {
            return Err(SoupValidationError::Cancelled);
        }
        for vertex in triangle {
            let valid = usize::try_from(*vertex)
                .ok()
                .is_some_and(|vertex| vertex < soup.positions.len());
            if !valid {
                return Err(SoupValidationError::Admission(format!(
                    "STEP tessellation face {face} references missing vertex {vertex}"
                )));
            }
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct EdgeOccurrence {
    edge: [u32; 2],
    face: usize,
    orientation: i8,
}

fn quality_and_localized_defects(
    soup: &Soup,
    cx: &Cx<'_>,
) -> Result<(MeshQuality, Vec<StepMeshDefect>, bool), PreprocessError> {
    let occurrence_count =
        soup.triangles
            .len()
            .checked_mul(3)
            .ok_or_else(|| PreprocessError::Resource {
                stage: "edge-occurrence-count",
                what: "edge-occurrence count overflowed".to_string(),
            })?;
    let mut occurrences = Vec::new();
    occurrences
        .try_reserve_exact(occurrence_count)
        .map_err(|_| PreprocessError::Resource {
            stage: "edge-occurrences",
            what: format!("allocation refused for {occurrence_count} edge occurrences"),
        })?;

    let reference = soup.positions[0];
    let mut signed_volume_six = 0.0f64;
    for (face, triangle) in soup.triangles.iter().enumerate() {
        if face % CANCELLATION_POLL_STRIDE == 0 {
            cx.checkpoint().map_err(|_| PreprocessError::Cancelled)?;
        }
        let [ia, ib, ic] = *triangle;
        let [a, b, c] = [ia, ib, ic].map(|vertex| soup.positions[vertex as usize]);
        let a = a.delta_from(reference);
        let b = b.delta_from(reference);
        let c = c.delta_from(reference);
        signed_volume_six += a.x * (b.y * c.z - b.z * c.y)
            + a.y * (b.z * c.x - b.x * c.z)
            + a.z * (b.x * c.y - b.y * c.x);
        for offset in 0..3 {
            let a = triangle[offset];
            let b = triangle[(offset + 1) % 3];
            let edge = [a.min(b), a.max(b)];
            let orientation = match a.cmp(&b) {
                core::cmp::Ordering::Less => 1,
                core::cmp::Ordering::Greater => -1,
                core::cmp::Ordering::Equal => 0,
            };
            occurrences.push(EdgeOccurrence {
                edge,
                face,
                orientation,
            });
        }
    }
    cx.checkpoint().map_err(|_| PreprocessError::Cancelled)?;
    occurrences.sort_unstable_by_key(|occurrence| {
        (occurrence.edge[0], occurrence.edge[1], occurrence.face)
    });
    cx.checkpoint().map_err(|_| PreprocessError::Cancelled)?;

    let mut defects = Vec::new();
    let mut truncated = false;
    let mut boundary_edges = 0usize;
    let mut nonmanifold_edges = 0usize;
    let mut orientation_conflicts = 0usize;
    let mut index = 0usize;
    while index < occurrences.len() {
        if index % CANCELLATION_POLL_STRIDE == 0 {
            cx.checkpoint().map_err(|_| PreprocessError::Cancelled)?;
        }
        let edge = occurrences[index].edge;
        let mut end = index;
        let mut orientation = 0i32;
        let mut faces = [None; 3];
        while end < occurrences.len() && occurrences[end].edge == edge {
            let local = end - index;
            if let Some(slot) = faces.get_mut(local) {
                *slot = Some(occurrences[end].face);
            }
            orientation += i32::from(occurrences[end].orientation);
            end += 1;
        }
        let count = end - index;
        let kind = if count == 1 {
            boundary_edges += 1;
            Some(StepMeshDefectKind::BoundaryEdge)
        } else if count > 2 {
            nonmanifold_edges += 1;
            Some(StepMeshDefectKind::NonManifoldEdge)
        } else if count == 2 && orientation != 0 {
            orientation_conflicts += 1;
            Some(StepMeshDefectKind::OrientationConflict)
        } else {
            None
        };
        if let Some(kind) = kind {
            push_defect_bounded(
                &mut defects,
                StepMeshDefect {
                    kind,
                    edge: Some(edge),
                    incident_faces: faces,
                    vertex: None,
                    total_incidents: u32::try_from(count).unwrap_or(u32::MAX),
                    detail: format!(
                        "edge [{}, {}] has {} incidents and orientation balance {}",
                        edge[0], edge[1], count, orientation
                    ),
                },
                &mut truncated,
            );
        }
        index = end;
    }
    let quality = MeshQuality {
        invalid_triangles: 0,
        boundary_edges,
        nonmanifold_edges,
        orientation_conflicts,
        outward_oriented: signed_volume_six.is_finite() && signed_volume_six > 0.0,
    };
    if !quality.outward_oriented {
        push_defect_bounded(
            &mut defects,
            StepMeshDefect {
                kind: StepMeshDefectKind::AggregateOrientation,
                edge: None,
                incident_faces: [None; 3],
                vertex: None,
                total_incidents: 0,
                detail: "aggregate signed orientation is not finite and outward-positive"
                    .to_string(),
            },
            &mut truncated,
        );
    }
    Ok((quality, defects, truncated))
}

fn push_defect_bounded(
    defects: &mut Vec<StepMeshDefect>,
    defect: StepMeshDefect,
    truncated: &mut bool,
) {
    if defects.len() < MAX_STEP_LOCALIZED_DEFECTS {
        defects.push(defect);
    } else {
        *truncated = true;
    }
}

enum PreprocessError {
    Cancelled,
    Resource { stage: &'static str, what: String },
}

fn compact_unreferenced_vertices(
    soup: &mut Soup,
    receipts: &mut Vec<RepairReceipt>,
    cx: &Cx<'_>,
) -> Result<(), PreprocessError> {
    if soup.triangles.is_empty() {
        return Ok(());
    }
    let mut referenced = Vec::<bool>::new();
    referenced
        .try_reserve_exact(soup.positions.len())
        .map_err(|_| PreprocessError::Resource {
            stage: "referenced-vertex-census",
            what: format!(
                "allocation refused for {} referenced-vertex flags",
                soup.positions.len()
            ),
        })?;
    referenced.resize(soup.positions.len(), false);
    for (face, triangle) in soup.triangles.iter().enumerate() {
        if face % CANCELLATION_POLL_STRIDE == 0 && cx.checkpoint().is_err() {
            return Err(PreprocessError::Cancelled);
        }
        for vertex in triangle {
            referenced[*vertex as usize] = true;
        }
    }
    let unreferenced = referenced
        .iter()
        .filter(|&&is_referenced| !is_referenced)
        .count();
    if unreferenced == 0 {
        return Ok(());
    }

    let mut remap = Vec::<u32>::new();
    remap
        .try_reserve_exact(soup.positions.len())
        .map_err(|_| PreprocessError::Resource {
            stage: "vertex-compaction-remap",
            what: format!(
                "allocation refused for {} vertex remap entries",
                soup.positions.len()
            ),
        })?;
    remap.resize(soup.positions.len(), u32::MAX);
    let retained = soup.positions.len() - unreferenced;
    let mut compacted = Vec::new();
    compacted
        .try_reserve_exact(retained)
        .map_err(|_| PreprocessError::Resource {
            stage: "vertex-compaction-output",
            what: format!("allocation refused for {retained} retained vertices"),
        })?;
    for (old, (&point, &is_referenced)) in soup.positions.iter().zip(&referenced).enumerate() {
        if old % CANCELLATION_POLL_STRIDE == 0 && cx.checkpoint().is_err() {
            return Err(PreprocessError::Cancelled);
        }
        if is_referenced {
            let new = u32::try_from(compacted.len()).map_err(|_| PreprocessError::Resource {
                stage: "vertex-compaction-index",
                what: "retained vertex index exceeded u32".to_string(),
            })?;
            remap[old] = new;
            compacted.push(point);
        }
    }
    for (face, triangle) in soup.triangles.iter_mut().enumerate() {
        if face % CANCELLATION_POLL_STRIDE == 0 && cx.checkpoint().is_err() {
            return Err(PreprocessError::Cancelled);
        }
        for vertex in triangle {
            *vertex = remap[*vertex as usize];
        }
    }
    receipts
        .try_reserve(1)
        .map_err(|_| PreprocessError::Resource {
            stage: "vertex-compaction-receipt",
            what: "allocation refused for unreferenced-vertex repair receipt".to_string(),
        })?;
    receipts.push(RepairReceipt {
        defect: "unreferenced-vertex",
        location: format!("{unreferenced} of {} source vertices", soup.positions.len()),
        action: "removed and deterministically reindexed retained vertices".to_string(),
    });
    soup.positions = compacted;
    Ok(())
}

fn observe_cancellation(
    cx: &Cx<'_>,
    source_fingerprint: u64,
    repairs: &[RepairReceipt],
) -> Result<(), StepImportRefusal> {
    cx.checkpoint().map_err(|_| StepImportRefusal::SdfBuild {
        source_fingerprint,
        error: MeshSdfError::Cancelled,
        repairs: repairs.to_vec(),
    })
}

#[derive(Debug, Clone, Copy)]
struct LinkOccurrence {
    vertex: u32,
    from: u32,
    to: u32,
    face: usize,
}

fn disconnected_vertex_links(
    soup: &Soup,
    cx: &Cx<'_>,
) -> Result<(Vec<StepMeshDefect>, bool), PreprocessError> {
    let occurrence_count =
        soup.triangles
            .len()
            .checked_mul(3)
            .ok_or_else(|| PreprocessError::Resource {
                stage: "vertex-link-occurrence-count",
                what: "vertex-link occurrence count overflowed".to_string(),
            })?;
    let mut links = Vec::new();
    links
        .try_reserve_exact(occurrence_count)
        .map_err(|_| PreprocessError::Resource {
            stage: "vertex-link-occurrences",
            what: format!("allocation refused for {occurrence_count} vertex-link occurrences"),
        })?;
    for (face, &[a, b, c]) in soup.triangles.iter().enumerate() {
        if face % CANCELLATION_POLL_STRIDE == 0 && cx.checkpoint().is_err() {
            return Err(PreprocessError::Cancelled);
        }
        links.extend_from_slice(&[
            LinkOccurrence {
                vertex: a,
                from: b,
                to: c,
                face,
            },
            LinkOccurrence {
                vertex: b,
                from: c,
                to: a,
                face,
            },
            LinkOccurrence {
                vertex: c,
                from: a,
                to: b,
                face,
            },
        ]);
    }
    cx.checkpoint().map_err(|_| PreprocessError::Cancelled)?;
    links.sort_unstable_by_key(|link| (link.vertex, link.from, link.to, link.face));
    cx.checkpoint().map_err(|_| PreprocessError::Cancelled)?;

    let mut defects = Vec::new();
    defects
        .try_reserve_exact(MAX_STEP_LOCALIZED_DEFECTS)
        .map_err(|_| PreprocessError::Resource {
            stage: "vertex-link-defects",
            what: format!(
                "allocation refused for {MAX_STEP_LOCALIZED_DEFECTS} vertex-link defects"
            ),
        })?;
    let mut truncated = false;
    let mut group_start = 0usize;
    let mut traversal_steps = 0usize;
    while group_start < links.len() {
        if group_start % CANCELLATION_POLL_STRIDE == 0 && cx.checkpoint().is_err() {
            return Err(PreprocessError::Cancelled);
        }
        let vertex = links[group_start].vertex;
        let mut group_end = group_start + 1;
        while group_end < links.len() && links[group_end].vertex == vertex {
            group_end += 1;
        }
        let group = &links[group_start..group_end];
        let expected = group.len();
        let duplicate_from = group
            .windows(2)
            .find(|pair| pair[0].from == pair[1].from)
            .map(|pair| pair[0].from);
        if let Some(from) = duplicate_from {
            push_defect_bounded(
                &mut defects,
                StepMeshDefect {
                    kind: StepMeshDefectKind::VertexLinkNonManifold,
                    edge: None,
                    incident_faces: [Some(group[0].face), None, None],
                    vertex: Some(vertex),
                    total_incidents: u32::try_from(expected).unwrap_or(u32::MAX),
                    detail: format!(
                        "vertex {vertex} link has multiple outgoing arcs from neighbor {from}"
                    ),
                },
                &mut truncated,
            );
            group_start = group_end;
            continue;
        }

        let start = group[0].from;
        let mut current = start;
        let mut visited = 0usize;
        let mut faces = [None; 3];
        let mut missing_arc = None;
        loop {
            if traversal_steps % CANCELLATION_POLL_STRIDE == 0 && cx.checkpoint().is_err() {
                return Err(PreprocessError::Cancelled);
            }
            traversal_steps = traversal_steps.saturating_add(1);
            let Ok(position) = group.binary_search_by_key(&current, |link| link.from) else {
                missing_arc = Some(current);
                break;
            };
            if let Some(face) = faces.get_mut(visited) {
                *face = Some(group[position].face);
            }
            current = group[position].to;
            visited += 1;
            if current == start || visited >= expected {
                break;
            }
        }
        if missing_arc.is_some() || current != start || visited != expected {
            push_defect_bounded(
                &mut defects,
                StepMeshDefect {
                    kind: StepMeshDefectKind::VertexLinkNonManifold,
                    edge: None,
                    incident_faces: faces,
                    vertex: Some(vertex),
                    total_incidents: u32::try_from(expected).unwrap_or(u32::MAX),
                    detail: format!(
                        "vertex {vertex} link is not one closed fan: traversed {visited} of {expected} incident faces{}",
                        missing_arc.map_or_else(String::new, |from| format!(
                            "; missing arc from neighbor {from}"
                        ))
                    ),
                },
                &mut truncated,
            );
        }
        group_start = group_end;
    }
    Ok((defects, truncated))
}

fn compose_estimated_error(
    mesh: NumericalCertificate,
    tessellation: NumericalCertificate,
) -> Result<NumericalCertificate, String> {
    if mesh.kind == NumericalKind::NoClaim
        || !mesh.lo.is_finite()
        || !mesh.hi.is_finite()
        || mesh.lo < 0.0
        || mesh.hi < 0.0
        || mesh.lo > mesh.hi
    {
        return Err(format!(
            "mesh-to-SDF evidence is not a finite ordered claim: {:?} [{}, {}]",
            mesh.kind, mesh.lo, mesh.hi
        ));
    }
    let sum = mesh.hi + tessellation.hi;
    if !sum.is_finite() || sum < 0.0 {
        return Err(format!(
            "mesh/SDF plus tessellation deviation overflowed: {} + {}",
            mesh.hi, tessellation.hi
        ));
    }
    let upper = sum.next_up();
    if !upper.is_finite() {
        return Err(format!(
            "outward rounding overflowed the mesh/SDF plus tessellation deviation sum {sum}"
        ));
    }
    Ok(NumericalCertificate::estimate(0.0, upper))
}

fn tessellation_fingerprint(soup: &Soup, cx: &Cx<'_>) -> Result<u64, MeshSdfError> {
    let mut hash = fnv1a64_update(
        0xcbf2_9ce4_8422_2325,
        STEP_TESSELLATION_FINGERPRINT_DOMAIN.as_bytes(),
    );
    hash = fnv1a64_update(hash, &(soup.positions.len() as u64).to_le_bytes());
    hash = fnv1a64_update(hash, &(soup.triangles.len() as u64).to_le_bytes());
    for (index, point) in soup.positions.iter().enumerate() {
        if index % CANCELLATION_POLL_STRIDE == 0 {
            cx.checkpoint().map_err(|_| MeshSdfError::Cancelled)?;
        }
        hash = fnv1a64_update(hash, &point.x.to_bits().to_le_bytes());
        hash = fnv1a64_update(hash, &point.y.to_bits().to_le_bytes());
        hash = fnv1a64_update(hash, &point.z.to_bits().to_le_bytes());
    }
    for (index, triangle) in soup.triangles.iter().enumerate() {
        if index % CANCELLATION_POLL_STRIDE == 0 {
            cx.checkpoint().map_err(|_| MeshSdfError::Cancelled)?;
        }
        for vertex in triangle {
            hash = fnv1a64_update(hash, &vertex.to_le_bytes());
        }
    }
    Ok(hash)
}

fn fnv1a64_update(mut hash: u64, bytes: &[u8]) -> u64 {
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

fn clone_schemas(schemas: &[String]) -> Result<Vec<String>, String> {
    let mut cloned = Vec::new();
    cloned
        .try_reserve(schemas.len())
        .map_err(|_| "allocation refused for STEP schema receipt".to_string())?;
    cloned.extend(schemas.iter().cloned());
    Ok(cloned)
}

const fn numerical_kind_label(kind: NumericalKind) -> &'static str {
    match kind {
        NumericalKind::Exact => "exact",
        NumericalKind::Enclosure => "enclosure",
        NumericalKind::Estimate => "estimate",
        NumericalKind::NoClaim => "no-claim",
    }
}

fn push_numerical_json(out: &mut String, certificate: NumericalCertificate) {
    let _ = write!(
        out,
        "{{\"kind\":\"{}\",\"lo\":{},\"hi\":{}}}",
        numerical_kind_label(certificate.kind),
        certificate.lo,
        certificate.hi
    );
}

fn push_json_string(out: &mut String, value: &str) {
    out.push('"');
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            character if character.is_control() => {
                let _ = write!(out, "\\u{:04x}", u32::from(character));
            }
            character => out.push(character),
        }
    }
    out.push('"');
}
