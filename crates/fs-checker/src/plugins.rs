//! Public discoverability facade for FrankenSim's package-bound semantic
//! certificate plugins.
//!
//! This module deliberately owns no second registry, raw `(family, version,
//! bytes)` verifier, payload encoder, arithmetic implementation, verdict type,
//! or witness hash. Its exports are the exact authoritative implementation
//! used by [`crate::check_with_capabilities`], strict-JSON checking, and the
//! release gate. A positive plugin decision is therefore reachable only through a complete
//! [`fs_package::EvidencePackage`] whose [`fs_package::SemanticWitness`] is
//! bound into its claim declaration, package root, origin subject, semantic
//! receipt, and release semantic context.
//!
//! The registry is closed and exact-versioned. An unknown family, unsupported
//! schema, malformed payload, resource overrun, arithmetic mismatch, or plugin
//! panic refuses semantic authority. Registry and report fingerprints are
//! domain-separated identities, not signatures or origin attestations.
//!
//! No-claim boundary: a positive semantic report says only that every attached
//! supported witness recomputed the package's declared finite interval under
//! the compiled checker semantics. It does not establish witness provenance,
//! authenticate the producing solver, prove that the model represents reality,
//! or turn a certificate refutation into a proof that the underlying theorem is
//! false. Integrity, semantic recomputation, origin authentication, falsifier
//! authentication, and release approval remain separate fail-closed axes.

pub use crate::semantic::{
    BOUNDED_LINF_RESIDUAL_FAMILY, EXACT_INTERVAL_FAMILY, INITIAL_SEMANTIC_SCHEMA_VERSION,
    MAX_INTERVAL_NODES, MAX_RESIDUAL_DIMENSION, MAX_RESIDUAL_MATRIX_ENTRIES,
    MAX_SEMANTIC_OPERATIONS, MAX_SEMANTIC_PAYLOAD_BYTES, MAX_SEMANTIC_WITNESS_BYTES,
    MAX_SEMANTIC_WITNESSES, SEMANTIC_IMPLEMENTATION_VERSION, SEMANTIC_PLUGIN_IDENTITY_DOMAIN,
    SEMANTIC_PLUGIN_IDENTITY_VERSION, SEMANTIC_REGISTRY_IDENTITY_DOMAIN,
    SEMANTIC_REGISTRY_IDENTITY_VERSION, SEMANTIC_REPORT_IDENTITY_DOMAIN,
    SEMANTIC_REPORT_IDENTITY_VERSION, SemanticClaimReceipt, SemanticClaimStatus, SemanticFailure,
    SemanticFailureKind, SemanticPluginDescriptor, SemanticReport, SemanticStatus,
    admit_retained_semantic_registry_fingerprint, semantic_plugin_registry,
    semantic_registry_fingerprint, verify_portable_semantics,
};
