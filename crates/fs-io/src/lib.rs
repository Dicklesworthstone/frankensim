//! fs-io (plan patch Rev J): import/export with QUARANTINE — the world
//! boundary. Real workflows ingest dirty geometry and emit useful
//! artifacts, but NO imported artifact becomes a trusted value without a
//! certification receipt: imports land as [`Quarantined`] and promote to
//! `Evidence` only after the repair suite + validity checks pass.
//! Parsers treat every byte as hostile: bounded resources, structured
//! rejection, never a panic.
//!
//! Layer: L2 (MORPH). Runtime deps: `std`, fs-rep-mesh (repair/validity),
//! fs-rep-sdf, fs-exec, fs-evidence, fs-geom, fs-obs, fs-math. PNG/EXR
//! export is fs-img's job (L5); ledger `imports` rows are written HELM-side
//! from the receipts this crate emits (L2 must not call L6).

pub mod catalog;
pub mod export;
pub mod obj;
pub mod ply;
pub mod quarantine;
pub mod selection;
pub mod step;
pub mod step_faceted;
pub mod step_faceted_export;
pub mod step_import;
pub mod stl;

pub use catalog::{
    CATALOG_CANCELLATION_POLL_STRIDE, CATALOG_CSV_PARSER_VERSION, CATALOG_JSON_PARSER_VERSION,
    CATALOG_READ_RECEIPT_VERSION, CATALOG_SCHEMA_VERSION, Catalog, CatalogCancellationEvidence,
    CatalogCsvLimits, CatalogFormat, CatalogInputIdentity, CatalogJsonLimits,
    CatalogProjectionLimits, CatalogRead, CatalogReadAuthority, CatalogReadCounters,
    CatalogReadLimits, CatalogReadReceipt, CatalogSchemaLimits, CatalogSchemaReceipt, ColumnKind,
    ColumnSpec, Schema, SchemaDefinitionRefusal,
};
pub use export::{export_3mf, export_glb, export_vtk};
pub use quarantine::{ImportDefect, ImportReceipt, PromotionRefusal, Quarantined, promote};
pub use selection::{
    AssignmentLimits, AssignmentReceipt, AssignmentRefusal, AssignmentReport, AssignmentRequest,
    AssignmentStats, HalfSpaceSide, MESH_ASSIGNMENT_POLL_STRIDE, MESH_ASSIGNMENT_SEMANTICS_VERSION,
    MeshSelector, NamedFaceGroup, ResolvedAssignment, resolve_mesh_assignments,
};
pub use step::{
    ParsedStep, STEP_SYNTAX_VERSION, StepDocument, StepEntity, StepHeader, StepInstance,
    StepLimits, StepProfileHint, StepStructureReceipt, StepValue, parse_step,
    parse_step_with_limits, write_step, write_step_with_limits,
};
pub use step_faceted::{
    DecodedFacetedBrep, STEP_FACETED_DECODER_VERSION, STEP_FACETED_MATERIALIZER_NAME,
    StepFacetedImportOutcome, StepFacetedImportRefusal, StepFacetedLimits, StepFacetedProfile,
    StepFacetedReceipt, StepFacetedRefusal, decode_faceted_brep_with_limits, import_faceted_brep,
};
pub use step_faceted_export::{
    ExportedFacetedBrep, STEP_FACETED_EXPORT_VERSION, StepFacetedExportLimits,
    StepFacetedExportMetadata, StepFacetedExportReceipt, StepFacetedExportRefusal,
    export_decoded_faceted_brep, export_decoded_faceted_brep_with_limits,
};
pub use step_import::{
    MAX_STEP_ADAPTER_ID_BYTES, MAX_STEP_LOCALIZED_DEFECTS, MAX_STEP_TESSELLATION_TRIANGLES,
    MAX_STEP_TESSELLATION_VERTICES, MAX_STEP_TOPOLOGY_AUXILIARY_BYTES,
    STEP_IMPORT_SEMANTICS_VERSION, STEP_TESSELLATION_FINGERPRINT_DOMAIN, StepImportOutcome,
    StepImportReceipt, StepImportRefusal, StepMeshDefect, StepMeshDefectKind,
    StepTessellatorIdentity, import_step_tessellation,
};

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Hard cap on elements a parser will allocate for (bounded resource
/// consumption on hostile input; ~100M vertices/faces).
pub const MAX_ELEMENTS: usize = 100_000_000;

/// Structured I/O failures (Decalogue P10) — hostile input is REFUSED,
/// never trusted or panicked on.
#[derive(Debug, Clone, PartialEq)]
pub enum IoError {
    /// Structurally invalid bytes for the declared format.
    Malformed {
        /// Byte or line position where parsing failed.
        at: usize,
        /// Diagnosis.
        what: String,
    },
    /// Valid-looking input outside the implemented subset.
    Unsupported {
        /// What was encountered.
        what: String,
    },
    /// Input exceeds declared resource bounds.
    ResourceBound {
        /// What bound.
        what: String,
    },
    /// Catalog value fails its schema.
    Schema {
        /// Row (1-based, excluding the header).
        row: usize,
        /// Column name.
        column: String,
        /// Diagnosis with the offending text.
        what: String,
    },
    /// An fs-exec cancellation request was observed before publication.
    Cancelled {
        /// Stable catalog stage that observed the request.
        stage: &'static str,
        /// Deterministic byte, row, field, or schema position within the stage.
        at: usize,
    },
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IoError::Malformed { at, what } => write!(f, "malformed input at {at}: {what}"),
            IoError::Unsupported { what } => write!(f, "unsupported: {what}"),
            IoError::ResourceBound { what } => write!(f, "resource bound exceeded: {what}"),
            IoError::Schema { row, column, what } => {
                write!(
                    f,
                    "schema violation at row {row}, column {column:?}: {what}"
                )
            }
            IoError::Cancelled { stage, at } => {
                write!(f, "catalog operation cancelled at {stage} position {at}")
            }
        }
    }
}

impl std::error::Error for IoError {}
