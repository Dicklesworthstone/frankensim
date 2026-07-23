//! Solver-free verification of portable certificate witnesses.
//!
//! This module is deliberately a CLOSED registry. External callbacks authenticate
//! origins; they cannot add a semantic family or manufacture a semantic pass.
//! Every built-in consumes a bounded canonical byte string iteratively and
//! recomputes its result from primitive inputs.

use std::panic::{AssertUnwindSafe, catch_unwind};

use fs_package::{ContentHash, EvidencePackage, hash_checker_decision};

/// Exact family id for the bounded interval straight-line program.
pub const EXACT_INTERVAL_FAMILY: &str = "frankensim/exact-interval";
/// Exact family id for the dense linear residual certificate.
pub const BOUNDED_LINF_RESIDUAL_FAMILY: &str = "frankensim/bounded-linf-residual";
/// The only schema version implemented by both initial families.
pub const INITIAL_SEMANTIC_SCHEMA_VERSION: u32 = 1;
/// Version of the checker-owned dispatch, arithmetic, and transcript rules.
pub const SEMANTIC_IMPLEMENTATION_VERSION: u32 = 1;

/// Maximum canonical bytes in one built-in witness.
pub const MAX_SEMANTIC_WITNESS_BYTES: usize = fs_package::MAX_SEMANTIC_WITNESS_PAYLOAD_BYTES;
/// Maximum interval-SLP nodes.
pub const MAX_INTERVAL_NODES: usize = 4_096;
/// Maximum dense residual rows or columns.
pub const MAX_RESIDUAL_DIMENSION: usize = 128;
/// Maximum dense residual matrix entries.
pub const MAX_RESIDUAL_MATRIX_ENTRIES: usize = 16_384;
/// Maximum witness-bearing claims in one package.
pub const MAX_SEMANTIC_WITNESSES: usize = fs_package::MAX_SEMANTIC_WITNESSES;
/// Maximum aggregate canonical witness bytes in one package.
pub const MAX_SEMANTIC_PAYLOAD_BYTES: usize = fs_package::MAX_SEMANTIC_WITNESS_TOTAL_BYTES;
/// Maximum charged primitive arithmetic operations in one package.
pub const MAX_SEMANTIC_OPERATIONS: u64 = 1_000_000;

/// Semantic version of one compiled verifier fingerprint.
pub const SEMANTIC_PLUGIN_IDENTITY_VERSION: u32 = 2;
/// Exact logical domain framed into every compiled verifier fingerprint.
pub const SEMANTIC_PLUGIN_IDENTITY_DOMAIN: &str = "fs-checker:semantic-plugin:v2";
/// Semantic version of the closed-registry fingerprint.
pub const SEMANTIC_REGISTRY_IDENTITY_VERSION: u32 = 2;
/// Exact logical domain framed into the closed-registry fingerprint.
pub const SEMANTIC_REGISTRY_IDENTITY_DOMAIN: &str = "fs-checker:semantic-registry:v2";
/// Semantic version of the callback-free semantic transcript digest.
pub const SEMANTIC_REPORT_IDENTITY_VERSION: u32 = 2;
/// Exact logical domain framed into every semantic transcript digest.
pub const SEMANTIC_REPORT_IDENTITY_DOMAIN: &str = "fs-checker:semantic-report:v2";

const _: () = assert!(MAX_SEMANTIC_WITNESS_BYTES == 256 * 1024);
const _: () = assert!(MAX_SEMANTIC_WITNESSES == 4_096);
const _: () = assert!(MAX_SEMANTIC_PAYLOAD_BYTES == 8 * 1024 * 1024);

const REGISTRY_REVISION: &[u8] =
    b"fs-checker-semantic-registry/v1;exact-interval-v1;bounded-linf-residual-v1";
const EXACT_I53_LIMIT: i64 = 9_007_199_254_740_992;

/// Owner-local compiled-plugin declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const SEMANTIC_PLUGIN_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-checker:semantic-plugin",
    "version_const=SEMANTIC_PLUGIN_IDENTITY_VERSION",
    "version=2",
    "domain=fs-checker:semantic-plugin:v2",
    "domain_const=SEMANTIC_PLUGIN_IDENTITY_DOMAIN",
    "encoder=SemanticPluginDescriptor::fingerprint",
    "encoder_helpers=semantic_plugin_fingerprint_with_schema",
    "schema_constants=SEMANTIC_PLUGIN_IDENTITY_VERSION,SEMANTIC_PLUGIN_IDENTITY_DOMAIN,SEMANTIC_IMPLEMENTATION_VERSION,EXACT_INTERVAL_FAMILY,BOUNDED_LINF_RESIDUAL_FAMILY,MAX_INTERVAL_NODES,MAX_RESIDUAL_DIMENSION,MAX_RESIDUAL_MATRIX_ENTRIES,MAX_SEMANTIC_WITNESS_BYTES,MAX_SEMANTIC_WITNESSES,MAX_SEMANTIC_PAYLOAD_BYTES,MAX_SEMANTIC_OPERATIONS,EXACT_I53_LIMIT,REGISTRY_REVISION,crates/fs-package/src/lib.rs#MAX_SEMANTIC_WITNESS_PAYLOAD_BYTES",
    "schema_functions=SemanticPluginDescriptor::admit_retained_fingerprint,verify_portable_semantics_after_integrity,dispatch_plugin,verify_exact_interval,verify_bounded_linf_residual,finite_scalar,exact_binary,node_ref,outward,next_up,next_down,bounded_detail,OperationBudget::new,OperationBudget::used,OperationBudget::charge,Cursor::new,Cursor::take,Cursor::u8,Cursor::u32,Cursor::u64,Cursor::i64,Cursor::finish,LocalInterval::finite,LocalInterval::point,LocalInterval::add,LocalInterval::sub,LocalInterval::mul,LocalInterval::div,LocalInterval::neg,LocalInterval::hull,LocalInterval::abs_upper,SlpValue::exact,SlpValue::interval,PluginError::malformed,PluginError::resource,PluginError::mismatch,PluginError::unknown_family,PluginError::unsupported_version,PluginError::panic,atom,stable_usize,admit_retained_semantic_hash,crates/fs-package/src/lib.rs#hash_checker_decision",
    "schema_dependencies=fs-package:semantic-witness",
    "digest=blake3-derive-key",
    "encoding=typed-binary",
    "sources=SemanticPluginDescriptor",
    "source_fields=SemanticPluginDescriptor.family:semantic,SemanticPluginDescriptor.schema_version:semantic,SemanticPluginDescriptor.maximum_payload_bytes:semantic",
    "source_bindings=SemanticPluginDescriptor.family>plugin-family,SemanticPluginDescriptor.schema_version>plugin-schema-version,SemanticPluginDescriptor.maximum_payload_bytes>maximum-payload-bytes",
    "external_semantic_fields=identity-version,digest-domain,implementation-version,family-specific-limits,registry-revision",
    "semantic_fields=identity-version,digest-domain,implementation-version,plugin-family,plugin-schema-version,maximum-payload-bytes,family-specific-limits,registry-revision",
    "excluded_fields=none",
    "consumers=SemanticPluginDescriptor::fingerprint,semantic_registry_fingerprint,SemanticClaimReceipt::plugin_fingerprint,verify_portable_semantics",
    "mutations=identity-version:crates/fs-checker/src/semantic.rs#semantic_identity_versions_and_transports_fail_closed,digest-domain:crates/fs-checker/src/semantic.rs#semantic_plugin_identity_fields_move_independently,implementation-version:crates/fs-checker/src/semantic.rs#semantic_plugin_identity_fields_move_independently,plugin-family:crates/fs-checker/src/semantic.rs#semantic_plugin_identity_fields_move_independently,plugin-schema-version:crates/fs-checker/src/semantic.rs#semantic_plugin_identity_fields_move_independently,maximum-payload-bytes:crates/fs-checker/src/semantic.rs#semantic_plugin_identity_fields_move_independently,family-specific-limits:crates/fs-checker/src/semantic.rs#semantic_plugin_identity_fields_move_independently,registry-revision:crates/fs-checker/src/semantic.rs#semantic_plugin_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_semantic_plugin_identity_fields",
    "transport_guard=SemanticPluginDescriptor::admit_retained_fingerprint",
    "version_guard=crates/fs-checker/src/semantic.rs#semantic_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-checker:semantic-plugin",
];

/// Owner-local closed-registry declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const SEMANTIC_REGISTRY_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-checker:semantic-registry",
    "version_const=SEMANTIC_REGISTRY_IDENTITY_VERSION",
    "version=2",
    "domain=fs-checker:semantic-registry:v2",
    "domain_const=SEMANTIC_REGISTRY_IDENTITY_DOMAIN",
    "encoder=semantic_registry_fingerprint",
    "encoder_helpers=semantic_registry_fingerprint_with_schema",
    "schema_constants=SEMANTIC_REGISTRY_IDENTITY_VERSION,SEMANTIC_REGISTRY_IDENTITY_DOMAIN,SEMANTIC_IMPLEMENTATION_VERSION,REGISTRY_REVISION,MAX_SEMANTIC_WITNESS_BYTES,MAX_SEMANTIC_WITNESSES,MAX_SEMANTIC_PAYLOAD_BYTES,MAX_SEMANTIC_OPERATIONS,crates/fs-package/src/lib.rs#MAX_SEMANTIC_WITNESS_FAMILY_BYTES,crates/fs-package/src/lib.rs#MAX_SEMANTIC_WITNESS_PAYLOAD_BYTES,crates/fs-package/src/lib.rs#MAX_SEMANTIC_WITNESSES,crates/fs-package/src/lib.rs#MAX_SEMANTIC_WITNESS_TOTAL_BYTES",
    "schema_functions=semantic_registry_schema,semantic_plugin_registry,SemanticPluginDescriptor::fingerprint,admit_retained_semantic_registry_fingerprint,admit_retained_semantic_hash,atom,stable_usize,crates/fs-package/src/lib.rs#hash_checker_decision",
    "schema_dependencies=fs-checker:semantic-plugin",
    "digest=blake3-derive-key",
    "encoding=typed-binary",
    "sources=SemanticRegistrySchema",
    "source_fields=SemanticRegistrySchema.revision:semantic,SemanticRegistrySchema.implementation_version:semantic,SemanticRegistrySchema.maximum_family_bytes:semantic,SemanticRegistrySchema.maximum_witness_bytes:semantic,SemanticRegistrySchema.maximum_witnesses:semantic,SemanticRegistrySchema.maximum_payload_bytes:semantic,SemanticRegistrySchema.maximum_operations:semantic,SemanticRegistrySchema.plugins:semantic",
    "source_bindings=SemanticRegistrySchema.revision>registry-revision,SemanticRegistrySchema.implementation_version>implementation-version,SemanticRegistrySchema.maximum_family_bytes>maximum-family-bytes,SemanticRegistrySchema.maximum_witness_bytes>maximum-witness-bytes,SemanticRegistrySchema.maximum_witnesses>maximum-witnesses,SemanticRegistrySchema.maximum_payload_bytes>maximum-payload-bytes,SemanticRegistrySchema.maximum_operations>maximum-operations,SemanticRegistrySchema.plugins>plugin-count+plugin-order+plugin-family+plugin-schema-version+plugin-fingerprint",
    "external_semantic_fields=identity-version,digest-domain",
    "semantic_fields=identity-version,digest-domain,registry-revision,implementation-version,maximum-family-bytes,maximum-witness-bytes,maximum-witnesses,maximum-payload-bytes,maximum-operations,plugin-count,plugin-order,plugin-family,plugin-schema-version,plugin-fingerprint",
    "excluded_fields=none",
    "consumers=semantic_registry_fingerprint,SemanticReport::registry_fingerprint,verify_portable_semantics,release-approval-signatures",
    "mutations=identity-version:crates/fs-checker/src/semantic.rs#semantic_identity_versions_and_transports_fail_closed,digest-domain:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,registry-revision:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,implementation-version:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,maximum-family-bytes:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,maximum-witness-bytes:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,maximum-witnesses:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,maximum-payload-bytes:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,maximum-operations:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,plugin-count:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,plugin-order:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,plugin-family:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,plugin-schema-version:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently,plugin-fingerprint:crates/fs-checker/src/semantic.rs#semantic_registry_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_semantic_registry_identity_fields",
    "transport_guard=admit_retained_semantic_registry_fingerprint",
    "version_guard=crates/fs-checker/src/semantic.rs#semantic_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-checker:semantic-registry",
];

/// Owner-local semantic-transcript declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const SEMANTIC_REPORT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-checker:semantic-report",
    "version_const=SEMANTIC_REPORT_IDENTITY_VERSION",
    "version=2",
    "domain=fs-checker:semantic-report:v2",
    "domain_const=SEMANTIC_REPORT_IDENTITY_DOMAIN",
    "encoder=semantic_context_hash",
    "encoder_helpers=semantic_context_hash_with_schema,hash_failure",
    "schema_constants=SEMANTIC_REPORT_IDENTITY_VERSION,SEMANTIC_REPORT_IDENTITY_DOMAIN",
    "schema_functions=SemanticReport::admit_retained_context_hash,SemanticReport::validate_context_hash,semantic_status_tag,semantic_claim_status_tag,SemanticFailureKind::tag,optional_hash,optional_u32,optional_usize,optional_atom,atom,stable_usize,admit_retained_semantic_hash,crates/fs-package/src/lib.rs#hash_checker_decision",
    "schema_dependencies=fs-checker:semantic-registry,fs-package:package-root,fs-package:claim-declaration,fs-package:semantic-witness",
    "digest=blake3-derive-key",
    "encoding=typed-binary",
    "sources=SemanticReport,SemanticClaimReceipt,SemanticFailure,SemanticStatus,SemanticClaimStatus,SemanticFailureKind",
    "source_fields=SemanticReport.status:derived:nested-status-variant-classified-separately,SemanticReport.package_root:semantic,SemanticReport.registry_fingerprint:semantic,SemanticReport.claims:semantic,SemanticReport.failures:semantic,SemanticReport.witnesses:semantic,SemanticReport.payload_bytes:semantic,SemanticReport.operations:semantic,SemanticReport.context_hash:derived:recomputed-from-semantic-fields,SemanticClaimReceipt.claim_index:semantic,SemanticClaimReceipt.claim_id:semantic,SemanticClaimReceipt.claim_hash:semantic,SemanticClaimReceipt.family:semantic,SemanticClaimReceipt.schema_version:semantic,SemanticClaimReceipt.witness_hash:semantic,SemanticClaimReceipt.plugin_fingerprint:semantic,SemanticClaimReceipt.status:derived:nested-status-variant-classified-separately,SemanticClaimReceipt.operations:semantic,SemanticClaimReceipt.failure:semantic,SemanticFailure.claim_index:semantic,SemanticFailure.claim_id:semantic,SemanticFailure.family:semantic,SemanticFailure.schema_version:semantic,SemanticFailure.kind:derived:nested-kind-variant-classified-separately,SemanticFailure.detail:semantic,SemanticStatus.variant:semantic,SemanticClaimStatus.variant:semantic,SemanticFailureKind.variant:semantic",
    "source_bindings=SemanticReport.package_root>package-root,SemanticReport.registry_fingerprint>registry-fingerprint,SemanticReport.claims>claim-count+claim-order,SemanticReport.failures>failure-count+failure-order,SemanticReport.witnesses>witness-count,SemanticReport.payload_bytes>payload-bytes,SemanticReport.operations>total-operations,SemanticClaimReceipt.claim_index>claim-index,SemanticClaimReceipt.claim_id>claim-id,SemanticClaimReceipt.claim_hash>claim-hash,SemanticClaimReceipt.family>claim-family-presence+claim-family,SemanticClaimReceipt.schema_version>claim-schema-version-presence+claim-schema-version,SemanticClaimReceipt.witness_hash>witness-hash-presence+witness-hash,SemanticClaimReceipt.plugin_fingerprint>plugin-fingerprint-presence+plugin-fingerprint,SemanticClaimReceipt.operations>claim-operations,SemanticClaimReceipt.failure>receipt-failure-presence,SemanticFailure.claim_index>failure-claim-index-presence+failure-claim-index,SemanticFailure.claim_id>failure-claim-id-presence+failure-claim-id,SemanticFailure.family>failure-family-presence+failure-family,SemanticFailure.schema_version>failure-schema-version-presence+failure-schema-version,SemanticFailure.detail>failure-detail,SemanticStatus.variant>package-status,SemanticClaimStatus.variant>claim-status,SemanticFailureKind.variant>failure-kind",
    "external_semantic_fields=identity-version,digest-domain",
    "semantic_fields=identity-version,digest-domain,package-root,registry-fingerprint,package-status,witness-count,payload-bytes,total-operations,claim-count,claim-order,claim-index,claim-id,claim-hash,claim-family-presence,claim-family,claim-schema-version-presence,claim-schema-version,witness-hash-presence,witness-hash,plugin-fingerprint-presence,plugin-fingerprint,claim-status,claim-operations,receipt-failure-presence,failure-count,failure-order,failure-claim-index-presence,failure-claim-index,failure-claim-id-presence,failure-claim-id,failure-family-presence,failure-family,failure-schema-version-presence,failure-schema-version,failure-kind,failure-detail",
    "excluded_fields=none",
    "consumers=SemanticReport::context_hash,SemanticReport::validate_context_hash,CheckReport::validate_decision_hash,SignaturePurpose::ReleaseApproval,release-approval-auditors",
    "mutations=identity-version:crates/fs-checker/src/semantic.rs#semantic_identity_versions_and_transports_fail_closed,digest-domain:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,package-root:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,registry-fingerprint:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,package-status:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,witness-count:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,payload-bytes:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,total-operations:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-count:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-order:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-index:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-id:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-hash:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-family-presence:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-family:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-schema-version-presence:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-schema-version:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,witness-hash-presence:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,witness-hash:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,plugin-fingerprint-presence:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,plugin-fingerprint:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-status:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,claim-operations:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,receipt-failure-presence:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-count:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-order:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-claim-index-presence:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-claim-index:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-claim-id-presence:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-claim-id:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-family-presence:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-family:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-schema-version-presence:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-schema-version:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-kind:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently,failure-detail:crates/fs-checker/src/semantic.rs#semantic_report_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_semantic_report_identity_fields",
    "transport_guard=SemanticReport::admit_retained_context_hash",
    "version_guard=crates/fs-checker/src/semantic.rs#semantic_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-checker:semantic-report",
];

/// One compiled-in semantic verifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SemanticPluginDescriptor {
    family: &'static str,
    schema_version: u32,
    maximum_payload_bytes: usize,
}

#[allow(dead_code)]
fn classify_semantic_plugin_identity_fields(descriptor: &SemanticPluginDescriptor) {
    let SemanticPluginDescriptor {
        family,
        schema_version,
        maximum_payload_bytes,
    } = descriptor;
    let _ = (family, schema_version, maximum_payload_bytes);
}

impl SemanticPluginDescriptor {
    /// Exact canonical family id.
    #[must_use]
    pub const fn family(&self) -> &'static str {
        self.family
    }

    /// Exact schema version; registry dispatch never performs version fallback.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Maximum accepted canonical payload size.
    #[must_use]
    pub const fn maximum_payload_bytes(&self) -> usize {
        self.maximum_payload_bytes
    }

    /// Stable verifier-semantics identity bound into positive receipts.
    #[must_use]
    pub fn fingerprint(&self) -> ContentHash {
        semantic_plugin_fingerprint_with_schema(
            self,
            SEMANTIC_PLUGIN_IDENTITY_DOMAIN.as_bytes(),
            SEMANTIC_IMPLEMENTATION_VERSION,
            MAX_INTERVAL_NODES,
            EXACT_I53_LIMIT,
            MAX_RESIDUAL_DIMENSION,
            MAX_RESIDUAL_MATRIX_ENTRIES,
            REGISTRY_REVISION,
        )
    }

    /// Admit a retained plugin fingerprint only under the exact v2 schema and
    /// fixed-width digest transport.
    #[must_use]
    pub fn admit_retained_fingerprint(version: u32, bytes: &[u8]) -> Option<ContentHash> {
        admit_retained_semantic_hash(version, SEMANTIC_PLUGIN_IDENTITY_VERSION, bytes)
    }
}

#[allow(clippy::too_many_arguments)]
fn semantic_plugin_fingerprint_with_schema(
    descriptor: &SemanticPluginDescriptor,
    domain: &[u8],
    implementation_version: u32,
    max_interval_nodes: usize,
    exact_i53_limit: i64,
    max_residual_dimension: usize,
    max_residual_matrix_entries: usize,
    registry_revision: &[u8],
) -> ContentHash {
    let mut bytes = Vec::new();
    atom(&mut bytes, domain);
    atom(&mut bytes, &implementation_version.to_le_bytes());
    atom(&mut bytes, descriptor.family.as_bytes());
    atom(&mut bytes, &descriptor.schema_version.to_le_bytes());
    atom(
        &mut bytes,
        &stable_usize(descriptor.maximum_payload_bytes).to_le_bytes(),
    );
    match descriptor.family {
        EXACT_INTERVAL_FAMILY => {
            atom(&mut bytes, &stable_usize(max_interval_nodes).to_le_bytes());
            atom(&mut bytes, &exact_i53_limit.to_le_bytes());
        }
        BOUNDED_LINF_RESIDUAL_FAMILY => {
            atom(
                &mut bytes,
                &stable_usize(max_residual_dimension).to_le_bytes(),
            );
            atom(
                &mut bytes,
                &stable_usize(max_residual_matrix_entries).to_le_bytes(),
            );
        }
        _ => atom(&mut bytes, b"unknown-compiled-plugin"),
    }
    atom(&mut bytes, registry_revision);
    hash_checker_decision(&bytes)
}

#[derive(Debug, Clone, Copy)]
struct SemanticRegistrySchema<'a> {
    revision: &'a [u8],
    implementation_version: u32,
    maximum_family_bytes: usize,
    maximum_witness_bytes: usize,
    maximum_witnesses: usize,
    maximum_payload_bytes: usize,
    maximum_operations: u64,
    plugins: &'a [SemanticPluginDescriptor],
}

#[allow(dead_code)]
fn classify_semantic_registry_identity_fields(schema: &SemanticRegistrySchema<'_>) {
    let SemanticRegistrySchema {
        revision,
        implementation_version,
        maximum_family_bytes,
        maximum_witness_bytes,
        maximum_witnesses,
        maximum_payload_bytes,
        maximum_operations,
        plugins,
    } = schema;
    let _ = (
        revision,
        implementation_version,
        maximum_family_bytes,
        maximum_witness_bytes,
        maximum_witnesses,
        maximum_payload_bytes,
        maximum_operations,
        plugins,
    );
}

const PLUGINS: [SemanticPluginDescriptor; 2] = [
    SemanticPluginDescriptor {
        family: EXACT_INTERVAL_FAMILY,
        schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION,
        maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES,
    },
    SemanticPluginDescriptor {
        family: BOUNDED_LINF_RESIDUAL_FAMILY,
        schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION,
        maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES,
    },
];

/// The fixed registry in deterministic dispatch order.
#[must_use]
pub const fn semantic_plugin_registry() -> &'static [SemanticPluginDescriptor] {
    &PLUGINS
}

fn semantic_registry_schema() -> SemanticRegistrySchema<'static> {
    SemanticRegistrySchema {
        revision: REGISTRY_REVISION,
        implementation_version: SEMANTIC_IMPLEMENTATION_VERSION,
        maximum_family_bytes: fs_package::MAX_SEMANTIC_WITNESS_FAMILY_BYTES,
        maximum_witness_bytes: MAX_SEMANTIC_WITNESS_BYTES,
        maximum_witnesses: MAX_SEMANTIC_WITNESSES,
        maximum_payload_bytes: MAX_SEMANTIC_PAYLOAD_BYTES,
        maximum_operations: MAX_SEMANTIC_OPERATIONS,
        plugins: semantic_plugin_registry(),
    }
}

/// Stable identity of the complete built-in registry.
#[must_use]
pub fn semantic_registry_fingerprint() -> ContentHash {
    semantic_registry_fingerprint_with_schema(
        &semantic_registry_schema(),
        SEMANTIC_REGISTRY_IDENTITY_DOMAIN.as_bytes(),
    )
}

fn semantic_registry_fingerprint_with_schema(
    schema: &SemanticRegistrySchema<'_>,
    domain: &[u8],
) -> ContentHash {
    let mut bytes = Vec::new();
    atom(&mut bytes, domain);
    atom(&mut bytes, schema.revision);
    atom(&mut bytes, &schema.implementation_version.to_le_bytes());
    atom(
        &mut bytes,
        &stable_usize(schema.maximum_family_bytes).to_le_bytes(),
    );
    atom(
        &mut bytes,
        &stable_usize(schema.maximum_witness_bytes).to_le_bytes(),
    );
    atom(
        &mut bytes,
        &stable_usize(schema.maximum_witnesses).to_le_bytes(),
    );
    atom(
        &mut bytes,
        &stable_usize(schema.maximum_payload_bytes).to_le_bytes(),
    );
    atom(&mut bytes, &schema.maximum_operations.to_le_bytes());
    for descriptor in schema.plugins {
        atom(&mut bytes, descriptor.family.as_bytes());
        atom(&mut bytes, &descriptor.schema_version.to_le_bytes());
        atom(&mut bytes, descriptor.fingerprint().as_bytes());
    }
    hash_checker_decision(&bytes)
}

/// Admit a retained closed-registry fingerprint only under the exact v2
/// schema and fixed-width digest transport.
#[must_use]
pub fn admit_retained_semantic_registry_fingerprint(
    version: u32,
    bytes: &[u8],
) -> Option<ContentHash> {
    admit_retained_semantic_hash(version, SEMANTIC_REGISTRY_IDENTITY_VERSION, bytes)
}

/// Package-level outcome of independent portable-witness verification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SemanticStatus {
    /// No claim attached a portable witness.
    NotProvided,
    /// Every attached portable witness was independently recomputed.
    Verified,
    /// At least one attached witness was unsupported, malformed, false, or over budget.
    Refused,
    /// Package integrity failed before semantic bytes were inspected.
    NotRun,
}

/// Per-claim portable-witness outcome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SemanticClaimStatus {
    /// This claim did not attach a witness.
    NotProvided,
    /// The exact compiled-in family/version recomputed the claimed interval.
    Verified,
    /// The witness was not admissible.
    Refused,
}

/// Stable refusal class. Details localize the exact bounded failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SemanticFailureKind {
    /// The package was not structurally safe to inspect.
    StructuralIntegrity,
    /// No compiled-in family had the exact declared id.
    UnknownFamily,
    /// The family exists but the exact schema version does not.
    UnsupportedVersion,
    /// The canonical payload was malformed or used forbidden arithmetic inputs.
    MalformedPayload,
    /// A byte, node, dimension, or operation bound was exceeded.
    ResourceLimit,
    /// Recomputed mathematics did not equal the declared finite Verified interval.
    ClaimMismatch,
    /// A built-in verifier panicked; authority is never retained across a panic.
    VerifierPanic,
}

impl SemanticFailureKind {
    const fn tag(self) -> &'static [u8] {
        match self {
            Self::StructuralIntegrity => b"structural-integrity",
            Self::UnknownFamily => b"unknown-family",
            Self::UnsupportedVersion => b"unsupported-version",
            Self::MalformedPayload => b"malformed-payload",
            Self::ResourceLimit => b"resource-limit",
            Self::ClaimMismatch => b"claim-mismatch",
            Self::VerifierPanic => b"verifier-panic",
        }
    }
}

/// One bounded semantic refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SemanticFailure {
    claim_index: Option<usize>,
    claim_id: Option<String>,
    family: Option<String>,
    schema_version: Option<u32>,
    kind: SemanticFailureKind,
    detail: String,
}

impl SemanticFailure {
    /// Claim index, absent for a package-integrity refusal.
    #[must_use]
    pub const fn claim_index(&self) -> Option<usize> {
        self.claim_index
    }

    /// Claim id, absent for a package-integrity refusal.
    #[must_use]
    pub fn claim_id(&self) -> Option<&str> {
        self.claim_id.as_deref()
    }

    /// Declared family, when a witness was present.
    #[must_use]
    pub fn family(&self) -> Option<&str> {
        self.family.as_deref()
    }

    /// Declared schema version, when a witness was present.
    #[must_use]
    pub const fn schema_version(&self) -> Option<u32> {
        self.schema_version
    }

    /// Stable refusal class.
    #[must_use]
    pub const fn kind(&self) -> SemanticFailureKind {
        self.kind
    }

    /// Bounded, deterministic failure localization.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

/// Sealed per-claim result bound into the semantic context.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SemanticClaimReceipt {
    claim_index: usize,
    claim_id: String,
    claim_hash: ContentHash,
    family: Option<String>,
    schema_version: Option<u32>,
    witness_hash: Option<ContentHash>,
    plugin_fingerprint: Option<ContentHash>,
    status: SemanticClaimStatus,
    operations: u64,
    failure: Option<SemanticFailure>,
}

impl SemanticClaimReceipt {
    /// Stable package position.
    #[must_use]
    pub const fn claim_index(&self) -> usize {
        self.claim_index
    }

    /// Stable claim identity.
    #[must_use]
    pub fn claim_id(&self) -> &str {
        &self.claim_id
    }

    /// Content hash of the complete raw claim declaration.
    #[must_use]
    pub const fn claim_hash(&self) -> ContentHash {
        self.claim_hash
    }

    /// Exact declared family, when present.
    #[must_use]
    pub fn family(&self) -> Option<&str> {
        self.family.as_deref()
    }

    /// Exact declared schema version, when present.
    #[must_use]
    pub const fn schema_version(&self) -> Option<u32> {
        self.schema_version
    }

    /// Domain-separated hash of family, version, and canonical payload.
    #[must_use]
    pub const fn witness_hash(&self) -> Option<ContentHash> {
        self.witness_hash
    }

    /// Compiled verifier identity, only after exact registry dispatch.
    #[must_use]
    pub const fn plugin_fingerprint(&self) -> Option<ContentHash> {
        self.plugin_fingerprint
    }

    /// Per-claim outcome.
    #[must_use]
    pub const fn status(&self) -> SemanticClaimStatus {
        self.status
    }

    /// Primitive operations charged to this witness, including failed work.
    #[must_use]
    pub const fn operations(&self) -> u64 {
        self.operations
    }

    /// Refusal detail, when any.
    #[must_use]
    pub const fn failure(&self) -> Option<&SemanticFailure> {
        self.failure.as_ref()
    }
}

/// Callback-free semantic transcript. Positive fields are sealed and hash-bound.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SemanticReport {
    status: SemanticStatus,
    package_root: ContentHash,
    registry_fingerprint: ContentHash,
    claims: Vec<SemanticClaimReceipt>,
    failures: Vec<SemanticFailure>,
    witnesses: usize,
    payload_bytes: usize,
    operations: u64,
    context_hash: ContentHash,
}

#[allow(dead_code)]
fn classify_semantic_report_identity_fields(
    report: &SemanticReport,
    receipt: &SemanticClaimReceipt,
    failure: &SemanticFailure,
    semantic_status: SemanticStatus,
    claim_status: SemanticClaimStatus,
    failure_kind: SemanticFailureKind,
) {
    let SemanticReport {
        status,
        package_root,
        registry_fingerprint,
        claims,
        failures,
        witnesses,
        payload_bytes,
        operations,
        context_hash,
    } = report;
    let SemanticClaimReceipt {
        claim_index,
        claim_id,
        claim_hash,
        family,
        schema_version,
        witness_hash,
        plugin_fingerprint,
        status: receipt_status,
        operations: receipt_operations,
        failure: receipt_failure,
    } = receipt;
    let SemanticFailure {
        claim_index: failure_claim_index,
        claim_id: failure_claim_id,
        family: failure_family,
        schema_version: failure_schema_version,
        kind,
        detail,
    } = failure;
    let _ = match semantic_status {
        SemanticStatus::NotProvided => 0_u8,
        SemanticStatus::Verified => 1_u8,
        SemanticStatus::Refused => 2_u8,
        SemanticStatus::NotRun => 3_u8,
    };
    let _ = match claim_status {
        SemanticClaimStatus::NotProvided => 0_u8,
        SemanticClaimStatus::Verified => 1_u8,
        SemanticClaimStatus::Refused => 2_u8,
    };
    let _ = match failure_kind {
        SemanticFailureKind::StructuralIntegrity => 0_u8,
        SemanticFailureKind::UnknownFamily => 1_u8,
        SemanticFailureKind::UnsupportedVersion => 2_u8,
        SemanticFailureKind::MalformedPayload => 3_u8,
        SemanticFailureKind::ResourceLimit => 4_u8,
        SemanticFailureKind::ClaimMismatch => 5_u8,
        SemanticFailureKind::VerifierPanic => 6_u8,
    };
    let _ = (
        status,
        package_root,
        registry_fingerprint,
        claims,
        failures,
        witnesses,
        payload_bytes,
        operations,
        context_hash,
        claim_index,
        claim_id,
        claim_hash,
        family,
        schema_version,
        witness_hash,
        plugin_fingerprint,
        receipt_status,
        receipt_operations,
        receipt_failure,
        failure_claim_index,
        failure_claim_id,
        failure_family,
        failure_schema_version,
        kind,
        detail,
    );
}

impl SemanticReport {
    /// Package-level semantic status.
    #[must_use]
    pub const fn status(&self) -> SemanticStatus {
        self.status
    }

    /// Package root to which this transcript applies.
    #[must_use]
    pub const fn package_root(&self) -> ContentHash {
        self.package_root
    }

    /// Complete compiled registry identity.
    #[must_use]
    pub const fn registry_fingerprint(&self) -> ContentHash {
        self.registry_fingerprint
    }

    /// One sealed result for every package claim, including `NotProvided`.
    #[must_use]
    pub fn claims(&self) -> &[SemanticClaimReceipt] {
        &self.claims
    }

    /// Ordered bounded refusals.
    #[must_use]
    pub fn failures(&self) -> &[SemanticFailure] {
        &self.failures
    }

    /// Number of attached witnesses.
    #[must_use]
    pub const fn witnesses(&self) -> usize {
        self.witnesses
    }

    /// Aggregate canonical witness bytes.
    #[must_use]
    pub const fn payload_bytes(&self) -> usize {
        self.payload_bytes
    }

    /// Aggregate charged primitive operations.
    #[must_use]
    pub const fn operations(&self) -> u64 {
        self.operations
    }

    /// Context producers bind into release approval before attaching a signature.
    #[must_use]
    pub const fn context_hash(&self) -> ContentHash {
        self.context_hash
    }

    /// Admit retained semantic-context bytes only under the exact v2 schema
    /// and fixed-width digest transport.
    #[must_use]
    pub fn admit_retained_context_hash(version: u32, bytes: &[u8]) -> Option<ContentHash> {
        admit_retained_semantic_hash(version, SEMANTIC_REPORT_IDENTITY_VERSION, bytes)
    }

    /// Recompute the complete transcript hash.
    #[must_use]
    pub fn recomputed_context_hash(&self) -> ContentHash {
        semantic_context_hash(self)
    }

    /// Whether the stored transcript hash matches every semantic authority field.
    #[must_use]
    pub fn validate_context_hash(&self) -> bool {
        self.context_hash == self.recomputed_context_hash()
    }

    pub(crate) fn not_run(detail: impl Into<String>) -> Self {
        let failure = SemanticFailure {
            claim_index: None,
            claim_id: None,
            family: None,
            schema_version: None,
            kind: SemanticFailureKind::StructuralIntegrity,
            detail: bounded_detail(detail.into()),
        };
        let mut report = Self {
            status: SemanticStatus::NotRun,
            package_root: ContentHash([0; 32]),
            registry_fingerprint: semantic_registry_fingerprint(),
            claims: Vec::new(),
            failures: vec![failure],
            witnesses: 0,
            payload_bytes: 0,
            operations: 0,
            context_hash: ContentHash([0; 32]),
        };
        report.context_hash = semantic_context_hash(&report);
        report
    }
}

/// Verify every attached witness without invoking any external origin or
/// signature callback. Producers use `context_hash()` from this report when
/// constructing release-approval signatures.
#[must_use]
pub fn verify_portable_semantics(package: &EvidencePackage) -> SemanticReport {
    let package_root = match package.verify_structural_integrity() {
        Ok(root) => root,
        Err(error) => return SemanticReport::not_run(format!("{error}")),
    };
    verify_portable_semantics_after_integrity(package, package_root)
}

pub(crate) fn verify_portable_semantics_after_integrity(
    package: &EvidencePackage,
    package_root: ContentHash,
) -> SemanticReport {
    let mut budget = OperationBudget::new(MAX_SEMANTIC_OPERATIONS);
    let mut witnesses = 0usize;
    let mut payload_bytes = 0usize;
    let mut claims = Vec::with_capacity(package.declared_claims_unverified().len());
    let mut failures = Vec::new();

    for (claim_index, claim) in package.declared_claims_unverified().iter().enumerate() {
        let claim_hash = claim.declared_content_hash_unverified();
        let Some(witness) = claim.declared_semantic_witness_unverified() else {
            claims.push(SemanticClaimReceipt {
                claim_index,
                claim_id: claim.id().to_string(),
                claim_hash,
                family: None,
                schema_version: None,
                witness_hash: None,
                plugin_fingerprint: None,
                status: SemanticClaimStatus::NotProvided,
                operations: 0,
                failure: None,
            });
            continue;
        };

        witnesses = match witnesses.checked_add(1) {
            Some(value) => value,
            None => MAX_SEMANTIC_WITNESSES.saturating_add(1),
        };
        let family = witness.family();
        let schema_version = witness.schema_version();
        let payload = witness.canonical_payload();
        let witness_hash = witness.content_hash();
        let before_operations = budget.used();

        let initial_failure = if witnesses > MAX_SEMANTIC_WITNESSES {
            Some(PluginError::resource("semantic witness count exceeds 4096"))
        } else {
            match payload_bytes.checked_add(payload.len()) {
                Some(total) if total <= MAX_SEMANTIC_PAYLOAD_BYTES => {
                    payload_bytes = total;
                    None
                }
                _ => Some(PluginError::resource(
                    "aggregate semantic payload exceeds 8 MiB",
                )),
            }
        };

        let descriptor = PLUGINS
            .iter()
            .find(|plugin| plugin.family == family && plugin.schema_version == schema_version);
        let dispatch = if let Some(error) = initial_failure {
            Err(error)
        } else if payload.len() > MAX_SEMANTIC_WITNESS_BYTES {
            Err(PluginError::resource(
                "witness payload exceeds the 256 KiB family limit",
            ))
        } else if let Some(descriptor) = descriptor {
            catch_unwind(AssertUnwindSafe(|| {
                dispatch_plugin(
                    descriptor,
                    payload,
                    claim.declared_verified_interval_unverified(),
                    &mut budget,
                )
            }))
            .unwrap_or_else(|_| Err(PluginError::panic()))
        } else if PLUGINS.iter().any(|plugin| plugin.family == family) {
            Err(PluginError::unsupported_version())
        } else {
            Err(PluginError::unknown_family())
        };
        let operations = budget.used().saturating_sub(before_operations);

        match dispatch {
            Ok(()) => claims.push(SemanticClaimReceipt {
                claim_index,
                claim_id: claim.id().to_string(),
                claim_hash,
                family: Some(family.to_string()),
                schema_version: Some(schema_version),
                witness_hash: Some(witness_hash),
                plugin_fingerprint: descriptor.map(SemanticPluginDescriptor::fingerprint),
                status: SemanticClaimStatus::Verified,
                operations,
                failure: None,
            }),
            Err(error) => {
                let failure = SemanticFailure {
                    claim_index: Some(claim_index),
                    claim_id: Some(claim.id().to_string()),
                    family: Some(family.to_string()),
                    schema_version: Some(schema_version),
                    kind: error.kind,
                    detail: bounded_detail(error.detail),
                };
                failures.push(failure.clone());
                claims.push(SemanticClaimReceipt {
                    claim_index,
                    claim_id: claim.id().to_string(),
                    claim_hash,
                    family: Some(family.to_string()),
                    schema_version: Some(schema_version),
                    witness_hash: Some(witness_hash),
                    plugin_fingerprint: descriptor.map(SemanticPluginDescriptor::fingerprint),
                    status: SemanticClaimStatus::Refused,
                    operations,
                    failure: Some(failure),
                });
            }
        }
    }

    let status = if !failures.is_empty() {
        SemanticStatus::Refused
    } else if witnesses == 0 {
        SemanticStatus::NotProvided
    } else {
        SemanticStatus::Verified
    };
    let mut report = SemanticReport {
        status,
        package_root,
        registry_fingerprint: semantic_registry_fingerprint(),
        claims,
        failures,
        witnesses,
        payload_bytes,
        operations: budget.used(),
        context_hash: ContentHash([0; 32]),
    };
    report.context_hash = semantic_context_hash(&report);
    report
}

fn semantic_context_hash(report: &SemanticReport) -> ContentHash {
    semantic_context_hash_with_schema(report, SEMANTIC_REPORT_IDENTITY_DOMAIN.as_bytes())
}

fn semantic_context_hash_with_schema(report: &SemanticReport, domain: &[u8]) -> ContentHash {
    let mut bytes = Vec::new();
    atom(&mut bytes, domain);
    atom(&mut bytes, report.package_root.as_bytes());
    atom(&mut bytes, report.registry_fingerprint.as_bytes());
    atom(&mut bytes, semantic_status_tag(report.status));
    atom(&mut bytes, &stable_usize(report.witnesses).to_le_bytes());
    atom(
        &mut bytes,
        &stable_usize(report.payload_bytes).to_le_bytes(),
    );
    atom(&mut bytes, &report.operations.to_le_bytes());
    atom(&mut bytes, &stable_usize(report.claims.len()).to_le_bytes());
    for receipt in &report.claims {
        atom(&mut bytes, &stable_usize(receipt.claim_index).to_le_bytes());
        atom(&mut bytes, receipt.claim_id.as_bytes());
        atom(&mut bytes, receipt.claim_hash.as_bytes());
        optional_atom(&mut bytes, receipt.family.as_deref().map(str::as_bytes));
        optional_u32(&mut bytes, receipt.schema_version);
        optional_hash(&mut bytes, receipt.witness_hash);
        optional_hash(&mut bytes, receipt.plugin_fingerprint);
        atom(&mut bytes, semantic_claim_status_tag(receipt.status));
        atom(&mut bytes, &receipt.operations.to_le_bytes());
        match &receipt.failure {
            Some(failure) => hash_failure(&mut bytes, failure),
            None => atom(&mut bytes, b"no-failure"),
        }
    }
    atom(
        &mut bytes,
        &stable_usize(report.failures.len()).to_le_bytes(),
    );
    for failure in &report.failures {
        hash_failure(&mut bytes, failure);
    }
    hash_checker_decision(&bytes)
}

fn hash_failure(bytes: &mut Vec<u8>, failure: &SemanticFailure) {
    atom(bytes, b"failure");
    optional_usize(bytes, failure.claim_index);
    optional_atom(bytes, failure.claim_id.as_deref().map(str::as_bytes));
    optional_atom(bytes, failure.family.as_deref().map(str::as_bytes));
    optional_u32(bytes, failure.schema_version);
    atom(bytes, failure.kind.tag());
    atom(bytes, failure.detail.as_bytes());
}

fn semantic_status_tag(status: SemanticStatus) -> &'static [u8] {
    match status {
        SemanticStatus::NotProvided => b"not-provided",
        SemanticStatus::Verified => b"verified",
        SemanticStatus::Refused => b"refused",
        SemanticStatus::NotRun => b"not-run",
    }
}

fn semantic_claim_status_tag(status: SemanticClaimStatus) -> &'static [u8] {
    match status {
        SemanticClaimStatus::NotProvided => b"not-provided",
        SemanticClaimStatus::Verified => b"verified",
        SemanticClaimStatus::Refused => b"refused",
    }
}

fn optional_hash(bytes: &mut Vec<u8>, value: Option<ContentHash>) {
    match value {
        Some(value) => {
            atom(bytes, b"some");
            atom(bytes, value.as_bytes());
        }
        None => atom(bytes, b"none"),
    }
}

fn optional_u32(bytes: &mut Vec<u8>, value: Option<u32>) {
    match value {
        Some(value) => {
            atom(bytes, b"some");
            atom(bytes, &value.to_le_bytes());
        }
        None => atom(bytes, b"none"),
    }
}

fn optional_usize(bytes: &mut Vec<u8>, value: Option<usize>) {
    match value {
        Some(value) => {
            atom(bytes, b"some");
            atom(bytes, &stable_usize(value).to_le_bytes());
        }
        None => atom(bytes, b"none"),
    }
}

fn optional_atom(bytes: &mut Vec<u8>, value: Option<&[u8]>) {
    match value {
        Some(value) => {
            atom(bytes, b"some");
            atom(bytes, value);
        }
        None => atom(bytes, b"none"),
    }
}

fn atom(bytes: &mut Vec<u8>, value: &[u8]) {
    bytes.extend_from_slice(&stable_usize(value.len()).to_le_bytes());
    bytes.extend_from_slice(value);
}

fn stable_usize(value: usize) -> u64 {
    u64::try_from(value).unwrap_or(u64::MAX)
}

fn admit_retained_semantic_hash(
    version: u32,
    expected_version: u32,
    bytes: &[u8],
) -> Option<ContentHash> {
    if version != expected_version || bytes.len() != 32 {
        return None;
    }
    let mut exact = [0_u8; 32];
    exact.copy_from_slice(bytes);
    Some(ContentHash(exact))
}

fn bounded_detail(mut detail: String) -> String {
    const MAX_DETAIL_BYTES: usize = 512;
    if detail.len() <= MAX_DETAIL_BYTES {
        return detail;
    }
    let mut end = MAX_DETAIL_BYTES;
    while !detail.is_char_boundary(end) {
        end -= 1;
    }
    detail.truncate(end);
    detail
}

#[derive(Debug)]
struct PluginError {
    kind: SemanticFailureKind,
    detail: String,
}

impl PluginError {
    fn malformed(detail: impl Into<String>) -> Self {
        Self {
            kind: SemanticFailureKind::MalformedPayload,
            detail: detail.into(),
        }
    }

    fn resource(detail: impl Into<String>) -> Self {
        Self {
            kind: SemanticFailureKind::ResourceLimit,
            detail: detail.into(),
        }
    }

    fn mismatch(detail: impl Into<String>) -> Self {
        Self {
            kind: SemanticFailureKind::ClaimMismatch,
            detail: detail.into(),
        }
    }

    fn unknown_family() -> Self {
        Self {
            kind: SemanticFailureKind::UnknownFamily,
            detail: "no built-in plugin has the exact declared family id".to_string(),
        }
    }

    fn unsupported_version() -> Self {
        Self {
            kind: SemanticFailureKind::UnsupportedVersion,
            detail: "the built-in family does not implement the exact declared schema version"
                .to_string(),
        }
    }

    fn panic() -> Self {
        Self {
            kind: SemanticFailureKind::VerifierPanic,
            detail: "built-in verifier panicked; no semantic authority retained".to_string(),
        }
    }
}

struct OperationBudget {
    used: u64,
    limit: u64,
}

impl OperationBudget {
    const fn new(limit: u64) -> Self {
        Self { used: 0, limit }
    }

    const fn used(&self) -> u64 {
        self.used
    }

    fn charge(&mut self, operations: u64) -> Result<(), PluginError> {
        let next = self
            .used
            .checked_add(operations)
            .ok_or_else(|| PluginError::resource("semantic operation count overflow"))?;
        if next > self.limit {
            return Err(PluginError::resource(
                "aggregate semantic operation limit exceeded",
            ));
        }
        self.used = next;
        Ok(())
    }
}

fn dispatch_plugin(
    descriptor: &SemanticPluginDescriptor,
    payload: &[u8],
    claimed: Option<(f64, f64)>,
    budget: &mut OperationBudget,
) -> Result<(), PluginError> {
    let claimed = claimed.ok_or_else(|| {
        PluginError::mismatch("portable semantic witnesses require a Verified interval claim")
    })?;
    if !claimed.0.is_finite() || !claimed.1.is_finite() || claimed.0 > claimed.1 {
        return Err(PluginError::mismatch(
            "portable semantic witnesses require an ordered finite claimed interval",
        ));
    }
    let result = match descriptor.family {
        EXACT_INTERVAL_FAMILY => verify_exact_interval(payload, budget)?,
        BOUNDED_LINF_RESIDUAL_FAMILY => verify_bounded_linf_residual(payload, budget)?,
        _ => return Err(PluginError::unknown_family()),
    };
    if result.lo.to_bits() != claimed.0.to_bits() || result.hi.to_bits() != claimed.1.to_bits() {
        return Err(PluginError::mismatch(format!(
            "recomputed interval bits [{:016x}, {:016x}] do not equal the claim [{:016x}, {:016x}]",
            result.lo.to_bits(),
            result.hi.to_bits(),
            claimed.0.to_bits(),
            claimed.1.to_bits()
        )));
    }
    Ok(())
}

struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, at: 0 }
    }

    fn take<const N: usize>(&mut self, what: &'static str) -> Result<[u8; N], PluginError> {
        let end = self
            .at
            .checked_add(N)
            .ok_or_else(|| PluginError::malformed("payload offset overflow"))?;
        let slice = self
            .bytes
            .get(self.at..end)
            .ok_or_else(|| PluginError::malformed(format!("truncated {what}")))?;
        self.at = end;
        slice
            .try_into()
            .map_err(|_| PluginError::malformed(format!("invalid {what}")))
    }

    fn u8(&mut self, what: &'static str) -> Result<u8, PluginError> {
        Ok(self.take::<1>(what)?[0])
    }

    fn u32(&mut self, what: &'static str) -> Result<u32, PluginError> {
        Ok(u32::from_le_bytes(self.take(what)?))
    }

    fn u64(&mut self, what: &'static str) -> Result<u64, PluginError> {
        Ok(u64::from_le_bytes(self.take(what)?))
    }

    fn i64(&mut self, what: &'static str) -> Result<i64, PluginError> {
        Ok(i64::from_le_bytes(self.take(what)?))
    }

    fn finish(self) -> Result<(), PluginError> {
        if self.at == self.bytes.len() {
            Ok(())
        } else {
            Err(PluginError::malformed(
                "trailing bytes after the canonical witness payload",
            ))
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct LocalInterval {
    lo: f64,
    hi: f64,
}

impl LocalInterval {
    fn finite(lo: f64, hi: f64) -> Result<Self, PluginError> {
        if !lo.is_finite() || !hi.is_finite() {
            return Err(PluginError::malformed(
                "NaN and infinite interval leaves are forbidden",
            ));
        }
        if lo > hi {
            return Err(PluginError::malformed("inverted interval leaf"));
        }
        Ok(Self { lo, hi })
    }

    const fn point(value: f64) -> Self {
        Self {
            lo: value,
            hi: value,
        }
    }

    fn add(self, other: Self) -> Result<Self, PluginError> {
        outward(self.lo + other.lo, self.hi + other.hi)
    }

    fn sub(self, other: Self) -> Result<Self, PluginError> {
        outward(self.lo - other.hi, self.hi - other.lo)
    }

    fn mul(self, other: Self) -> Result<Self, PluginError> {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        if products.iter().any(|value| value.is_nan()) {
            return Err(PluginError::malformed(
                "interval multiplication produced an indeterminate endpoint",
            ));
        }
        let mut lo = products[0];
        let mut hi = products[0];
        for value in &products[1..] {
            lo = lo.min(*value);
            hi = hi.max(*value);
        }
        outward(lo, hi)
    }

    fn div(self, other: Self) -> Result<Self, PluginError> {
        if other.lo <= 0.0 && 0.0 <= other.hi {
            return Err(PluginError::malformed(
                "division by a zero-containing interval is not admitted",
            ));
        }
        let quotients = [
            self.lo / other.lo,
            self.lo / other.hi,
            self.hi / other.lo,
            self.hi / other.hi,
        ];
        if quotients.iter().any(|value| value.is_nan()) {
            return Err(PluginError::malformed(
                "interval division produced an indeterminate endpoint",
            ));
        }
        let mut lo = quotients[0];
        let mut hi = quotients[0];
        for value in &quotients[1..] {
            lo = lo.min(*value);
            hi = hi.max(*value);
        }
        outward(lo, hi)
    }

    const fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    fn hull(self, other: Self) -> Self {
        Self {
            lo: self.lo.min(other.lo),
            hi: self.hi.max(other.hi),
        }
    }

    fn abs_upper(self) -> Result<f64, PluginError> {
        let upper = self.lo.abs().max(self.hi.abs());
        if upper.is_finite() {
            Ok(upper)
        } else {
            Err(PluginError::malformed(
                "residual arithmetic produced a non-finite bound",
            ))
        }
    }
}

fn outward(lo: f64, hi: f64) -> Result<LocalInterval, PluginError> {
    if lo.is_nan() || hi.is_nan() {
        return Err(PluginError::malformed("arithmetic produced a NaN endpoint"));
    }
    Ok(LocalInterval {
        lo: next_down(lo),
        hi: next_up(hi),
    })
}

fn next_up(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        return value;
    }
    if value == 0.0 {
        return f64::from_bits(1);
    }
    if value > 0.0 {
        f64::from_bits(value.to_bits() + 1)
    } else {
        f64::from_bits(value.to_bits() - 1)
    }
}

fn next_down(value: f64) -> f64 {
    if value.is_nan() || value == f64::NEG_INFINITY {
        return value;
    }
    if value == 0.0 {
        return f64::from_bits((1_u64 << 63) | 1);
    }
    if value > 0.0 {
        f64::from_bits(value.to_bits() - 1)
    } else {
        f64::from_bits(value.to_bits() + 1)
    }
}

#[derive(Debug, Clone, Copy)]
struct SlpValue {
    interval: LocalInterval,
    exact_i53: Option<i64>,
}

impl SlpValue {
    fn exact(value: i64) -> Result<Self, PluginError> {
        if !(-EXACT_I53_LIMIT..=EXACT_I53_LIMIT).contains(&value) {
            return Err(PluginError::malformed(
                "exact-i53 leaf is outside the exactly representable range",
            ));
        }
        Ok(Self {
            interval: LocalInterval::point(value as f64),
            exact_i53: Some(value),
        })
    }

    const fn interval(interval: LocalInterval) -> Self {
        Self {
            interval,
            exact_i53: None,
        }
    }
}

fn exact_binary(
    left: SlpValue,
    right: SlpValue,
    exact: impl FnOnce(i64, i64) -> Option<i64>,
    interval: impl FnOnce(LocalInterval, LocalInterval) -> Result<LocalInterval, PluginError>,
) -> Result<SlpValue, PluginError> {
    if let (Some(left), Some(right)) = (left.exact_i53, right.exact_i53)
        && let Some(result) = exact(left, right)
        && (-EXACT_I53_LIMIT..=EXACT_I53_LIMIT).contains(&result)
    {
        return SlpValue::exact(result);
    }
    Ok(SlpValue::interval(interval(left.interval, right.interval)?))
}

fn node_ref(values: &[SlpValue], index: u32, current: usize) -> Result<SlpValue, PluginError> {
    let index = usize::try_from(index)
        .map_err(|_| PluginError::malformed("node reference does not fit usize"))?;
    if index >= current {
        return Err(PluginError::malformed(
            "node references must point strictly backward",
        ));
    }
    values
        .get(index)
        .copied()
        .ok_or_else(|| PluginError::malformed("node reference is out of range"))
}

fn verify_exact_interval(
    payload: &[u8],
    budget: &mut OperationBudget,
) -> Result<LocalInterval, PluginError> {
    let mut cursor = Cursor::new(payload);
    let nodes = usize::try_from(cursor.u32("interval node count")?)
        .map_err(|_| PluginError::resource("interval node count does not fit usize"))?;
    if nodes == 0 || nodes > MAX_INTERVAL_NODES {
        return Err(PluginError::resource(
            "interval node count must be in 1..=4096",
        ));
    }
    let mut values = Vec::with_capacity(nodes);
    for current in 0..nodes {
        let tag = cursor.u8("interval node tag")?;
        let value =
            match tag {
                0 => {
                    budget.charge(1)?;
                    SlpValue::exact(cursor.i64("exact-i53 value")?)?
                }
                1 => {
                    budget.charge(1)?;
                    let lo = f64::from_bits(cursor.u64("interval lower bits")?);
                    let hi = f64::from_bits(cursor.u64("interval upper bits")?);
                    SlpValue::interval(LocalInterval::finite(lo, hi)?)
                }
                2 | 3 | 4 | 5 | 7 => {
                    let left = node_ref(&values, cursor.u32("left node reference")?, current)?;
                    let right = node_ref(&values, cursor.u32("right node reference")?, current)?;
                    match tag {
                        2 => {
                            budget.charge(3)?;
                            exact_binary(left, right, i64::checked_add, LocalInterval::add)?
                        }
                        3 => {
                            budget.charge(3)?;
                            exact_binary(left, right, i64::checked_sub, LocalInterval::sub)?
                        }
                        4 => {
                            budget.charge(8)?;
                            exact_binary(left, right, i64::checked_mul, LocalInterval::mul)?
                        }
                        5 => {
                            budget.charge(8)?;
                            exact_binary(
                                left,
                                right,
                                |left, right| {
                                    (right != 0 && left.checked_rem(right) == Some(0))
                                        .then(|| left.checked_div(right))
                                        .flatten()
                                },
                                LocalInterval::div,
                            )?
                        }
                        7 => {
                            budget.charge(1)?;
                            SlpValue::interval(left.interval.hull(right.interval))
                        }
                        _ => unreachable!("closed interval tag match"),
                    }
                }
                6 => {
                    let input = node_ref(&values, cursor.u32("unary node reference")?, current)?;
                    budget.charge(1)?;
                    if let Some(exact) = input.exact_i53 {
                        SlpValue::exact(exact.checked_neg().ok_or_else(|| {
                            PluginError::malformed("exact-i53 negation overflow")
                        })?)?
                    } else {
                        SlpValue::interval(input.interval.neg())
                    }
                }
                _ => return Err(PluginError::malformed("unknown interval node tag")),
            };
        values.push(value);
    }
    let result = usize::try_from(cursor.u32("interval result index")?)
        .map_err(|_| PluginError::malformed("result index does not fit usize"))?;
    cursor.finish()?;
    values
        .get(result)
        .map(|value| value.interval)
        .ok_or_else(|| PluginError::malformed("result index is out of range"))
}

fn finite_scalar(cursor: &mut Cursor<'_>, what: &'static str) -> Result<f64, PluginError> {
    let value = f64::from_bits(cursor.u64(what)?);
    if value.is_finite() {
        Ok(value)
    } else {
        Err(PluginError::malformed(format!(
            "{what} must be finite (NaN and infinity are forbidden)"
        )))
    }
}

fn verify_bounded_linf_residual(
    payload: &[u8],
    budget: &mut OperationBudget,
) -> Result<LocalInterval, PluginError> {
    let mut cursor = Cursor::new(payload);
    if cursor.u8("residual norm tag")? != 0 {
        return Err(PluginError::malformed(
            "residual schema v1 supports only norm tag 0 (L-infinity)",
        ));
    }
    let rows = usize::try_from(cursor.u32("residual row count")?)
        .map_err(|_| PluginError::resource("row count does not fit usize"))?;
    let cols = usize::try_from(cursor.u32("residual column count")?)
        .map_err(|_| PluginError::resource("column count does not fit usize"))?;
    if rows == 0 || cols == 0 || rows > MAX_RESIDUAL_DIMENSION || cols > MAX_RESIDUAL_DIMENSION {
        return Err(PluginError::resource(
            "residual rows and columns must each be in 1..=128",
        ));
    }
    let entries = rows
        .checked_mul(cols)
        .ok_or_else(|| PluginError::resource("residual matrix size overflow"))?;
    if entries > MAX_RESIDUAL_MATRIX_ENTRIES {
        return Err(PluginError::resource(
            "residual matrix exceeds 16384 entries",
        ));
    }

    let mut matrix = Vec::with_capacity(entries);
    for _ in 0..entries {
        matrix.push(finite_scalar(&mut cursor, "matrix scalar")?);
    }
    let mut candidate = Vec::with_capacity(cols);
    for _ in 0..cols {
        candidate.push(finite_scalar(&mut cursor, "candidate scalar")?);
    }
    let mut right_hand_side = Vec::with_capacity(rows);
    for _ in 0..rows {
        right_hand_side.push(finite_scalar(&mut cursor, "right-hand-side scalar")?);
    }
    cursor.finish()?;

    let per_entry = u64::try_from(entries)
        .ok()
        .and_then(|value| value.checked_mul(11))
        .ok_or_else(|| PluginError::resource("residual operation count overflow"))?;
    let row_cost = u64::try_from(rows)
        .ok()
        .and_then(|value| value.checked_mul(3))
        .ok_or_else(|| PluginError::resource("residual operation count overflow"))?;
    budget.charge(
        per_entry
            .checked_add(row_cost)
            .ok_or_else(|| PluginError::resource("residual operation count overflow"))?,
    )?;

    let mut maximum = 0.0_f64;
    for row in 0..rows {
        let mut product = LocalInterval::point(0.0);
        for column in 0..cols {
            let coefficient = LocalInterval::point(matrix[row * cols + column]);
            let value = LocalInterval::point(candidate[column]);
            product = product.add(coefficient.mul(value)?)?;
        }
        let residual = LocalInterval::point(right_hand_side[row]).sub(product)?;
        maximum = maximum.max(residual.abs_upper()?);
    }
    Ok(LocalInterval {
        lo: 0.0,
        hi: maximum,
    })
}

#[cfg(test)]
mod identity_tests {
    use super::*;

    static ONE_PLUGIN: [SemanticPluginDescriptor; 1] = [SemanticPluginDescriptor {
        family: EXACT_INTERVAL_FAMILY,
        schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION,
        maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES,
    }];
    static REVERSED_PLUGINS: [SemanticPluginDescriptor; 2] = [
        SemanticPluginDescriptor {
            family: BOUNDED_LINF_RESIDUAL_FAMILY,
            schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION,
            maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES,
        },
        SemanticPluginDescriptor {
            family: EXACT_INTERVAL_FAMILY,
            schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION,
            maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES,
        },
    ];
    static FAMILY_CHANGED_PLUGINS: [SemanticPluginDescriptor; 2] = [
        SemanticPluginDescriptor {
            family: "frankensim/mutated-family",
            schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION,
            maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES,
        },
        SemanticPluginDescriptor {
            family: BOUNDED_LINF_RESIDUAL_FAMILY,
            schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION,
            maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES,
        },
    ];
    static SCHEMA_CHANGED_PLUGINS: [SemanticPluginDescriptor; 2] = [
        SemanticPluginDescriptor {
            family: EXACT_INTERVAL_FAMILY,
            schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION + 1,
            maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES,
        },
        SemanticPluginDescriptor {
            family: BOUNDED_LINF_RESIDUAL_FAMILY,
            schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION,
            maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES,
        },
    ];
    static FINGERPRINT_CHANGED_PLUGINS: [SemanticPluginDescriptor; 2] = [
        SemanticPluginDescriptor {
            family: EXACT_INTERVAL_FAMILY,
            schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION,
            maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES + 1,
        },
        SemanticPluginDescriptor {
            family: BOUNDED_LINF_RESIDUAL_FAMILY,
            schema_version: INITIAL_SEMANTIC_SCHEMA_VERSION,
            maximum_payload_bytes: MAX_SEMANTIC_WITNESS_BYTES,
        },
    ];

    fn semantic_failure_fixture() -> SemanticFailure {
        SemanticFailure {
            claim_index: Some(3),
            claim_id: Some("claim-three".to_string()),
            family: Some(EXACT_INTERVAL_FAMILY.to_string()),
            schema_version: Some(INITIAL_SEMANTIC_SCHEMA_VERSION),
            kind: SemanticFailureKind::ClaimMismatch,
            detail: "fixture mismatch".to_string(),
        }
    }

    fn semantic_report_fixture() -> SemanticReport {
        let first_failure = semantic_failure_fixture();
        let second_failure = SemanticFailure {
            claim_index: Some(9),
            claim_id: Some("claim-nine".to_string()),
            family: Some(BOUNDED_LINF_RESIDUAL_FAMILY.to_string()),
            schema_version: Some(INITIAL_SEMANTIC_SCHEMA_VERSION),
            kind: SemanticFailureKind::MalformedPayload,
            detail: "fixture payload".to_string(),
        };
        let first_receipt = SemanticClaimReceipt {
            claim_index: 3,
            claim_id: "claim-three".to_string(),
            claim_hash: ContentHash([3; 32]),
            family: Some(EXACT_INTERVAL_FAMILY.to_string()),
            schema_version: Some(INITIAL_SEMANTIC_SCHEMA_VERSION),
            witness_hash: Some(ContentHash([4; 32])),
            plugin_fingerprint: Some(PLUGINS[0].fingerprint()),
            status: SemanticClaimStatus::Refused,
            operations: 11,
            failure: Some(first_failure.clone()),
        };
        let second_receipt = SemanticClaimReceipt {
            claim_index: 9,
            claim_id: "claim-nine".to_string(),
            claim_hash: ContentHash([9; 32]),
            family: Some(BOUNDED_LINF_RESIDUAL_FAMILY.to_string()),
            schema_version: Some(INITIAL_SEMANTIC_SCHEMA_VERSION),
            witness_hash: Some(ContentHash([10; 32])),
            plugin_fingerprint: Some(PLUGINS[1].fingerprint()),
            status: SemanticClaimStatus::Verified,
            operations: 17,
            failure: None,
        };
        let mut report = SemanticReport {
            status: SemanticStatus::Refused,
            package_root: ContentHash([1; 32]),
            registry_fingerprint: semantic_registry_fingerprint(),
            claims: vec![first_receipt, second_receipt],
            failures: vec![first_failure, second_failure],
            witnesses: 2,
            payload_bytes: 29,
            operations: 28,
            context_hash: ContentHash([0; 32]),
        };
        report.context_hash = semantic_context_hash(&report);
        report
    }

    #[test]
    fn semantic_identity_versions_and_transports_fail_closed() {
        assert_eq!(SEMANTIC_PLUGIN_IDENTITY_VERSION, 2);
        assert_eq!(
            SEMANTIC_PLUGIN_IDENTITY_DOMAIN,
            "fs-checker:semantic-plugin:v2"
        );
        assert_eq!(SEMANTIC_REGISTRY_IDENTITY_VERSION, 2);
        assert_eq!(
            SEMANTIC_REGISTRY_IDENTITY_DOMAIN,
            "fs-checker:semantic-registry:v2"
        );
        assert_eq!(SEMANTIC_REPORT_IDENTITY_VERSION, 2);
        assert_eq!(
            SEMANTIC_REPORT_IDENTITY_DOMAIN,
            "fs-checker:semantic-report:v2"
        );

        let plugin = PLUGINS[0].fingerprint();
        let registry = semantic_registry_fingerprint();
        let report = semantic_report_fixture();
        assert_eq!(
            SemanticPluginDescriptor::admit_retained_fingerprint(
                SEMANTIC_PLUGIN_IDENTITY_VERSION,
                plugin.as_bytes(),
            ),
            Some(plugin)
        );
        assert_eq!(
            admit_retained_semantic_registry_fingerprint(
                SEMANTIC_REGISTRY_IDENTITY_VERSION,
                registry.as_bytes(),
            ),
            Some(registry)
        );
        assert_eq!(
            SemanticReport::admit_retained_context_hash(
                SEMANTIC_REPORT_IDENTITY_VERSION,
                report.context_hash().as_bytes(),
            ),
            Some(report.context_hash())
        );

        for stale in [0, 1, 3, u32::MAX] {
            assert_eq!(
                SemanticPluginDescriptor::admit_retained_fingerprint(stale, plugin.as_bytes()),
                None
            );
            assert_eq!(
                admit_retained_semantic_registry_fingerprint(stale, registry.as_bytes()),
                None
            );
            assert_eq!(
                SemanticReport::admit_retained_context_hash(
                    stale,
                    report.context_hash().as_bytes(),
                ),
                None
            );
        }
        for malformed in [&[0_u8; 31][..], &[0_u8; 33][..]] {
            assert_eq!(
                SemanticPluginDescriptor::admit_retained_fingerprint(
                    SEMANTIC_PLUGIN_IDENTITY_VERSION,
                    malformed,
                ),
                None
            );
            assert_eq!(
                admit_retained_semantic_registry_fingerprint(
                    SEMANTIC_REGISTRY_IDENTITY_VERSION,
                    malformed,
                ),
                None
            );
            assert_eq!(
                SemanticReport::admit_retained_context_hash(
                    SEMANTIC_REPORT_IDENTITY_VERSION,
                    malformed,
                ),
                None
            );
        }
    }

    #[test]
    fn semantic_plugin_identity_fields_move_independently() {
        let descriptor = PLUGINS[0];
        let baseline = descriptor.fingerprint();
        macro_rules! assert_moves {
            ($descriptor:expr, $domain:expr, $implementation:expr, $nodes:expr, $i53:expr, $rows:expr, $entries:expr, $revision:expr) => {
                assert_ne!(
                    semantic_plugin_fingerprint_with_schema(
                        &$descriptor,
                        $domain,
                        $implementation,
                        $nodes,
                        $i53,
                        $rows,
                        $entries,
                        $revision,
                    ),
                    baseline
                );
            };
        }
        assert_moves!(
            descriptor,
            b"fs-checker:semantic-plugin:v2-mutated",
            SEMANTIC_IMPLEMENTATION_VERSION,
            MAX_INTERVAL_NODES,
            EXACT_I53_LIMIT,
            MAX_RESIDUAL_DIMENSION,
            MAX_RESIDUAL_MATRIX_ENTRIES,
            REGISTRY_REVISION
        );
        assert_moves!(
            descriptor,
            SEMANTIC_PLUGIN_IDENTITY_DOMAIN.as_bytes(),
            SEMANTIC_IMPLEMENTATION_VERSION + 1,
            MAX_INTERVAL_NODES,
            EXACT_I53_LIMIT,
            MAX_RESIDUAL_DIMENSION,
            MAX_RESIDUAL_MATRIX_ENTRIES,
            REGISTRY_REVISION
        );
        assert_moves!(
            SemanticPluginDescriptor {
                family: BOUNDED_LINF_RESIDUAL_FAMILY,
                schema_version: descriptor.schema_version,
                maximum_payload_bytes: descriptor.maximum_payload_bytes,
            },
            SEMANTIC_PLUGIN_IDENTITY_DOMAIN.as_bytes(),
            SEMANTIC_IMPLEMENTATION_VERSION,
            MAX_INTERVAL_NODES,
            EXACT_I53_LIMIT,
            MAX_RESIDUAL_DIMENSION,
            MAX_RESIDUAL_MATRIX_ENTRIES,
            REGISTRY_REVISION
        );
        assert_moves!(
            SemanticPluginDescriptor {
                schema_version: descriptor.schema_version + 1,
                ..descriptor
            },
            SEMANTIC_PLUGIN_IDENTITY_DOMAIN.as_bytes(),
            SEMANTIC_IMPLEMENTATION_VERSION,
            MAX_INTERVAL_NODES,
            EXACT_I53_LIMIT,
            MAX_RESIDUAL_DIMENSION,
            MAX_RESIDUAL_MATRIX_ENTRIES,
            REGISTRY_REVISION
        );
        assert_moves!(
            SemanticPluginDescriptor {
                maximum_payload_bytes: descriptor.maximum_payload_bytes + 1,
                ..descriptor
            },
            SEMANTIC_PLUGIN_IDENTITY_DOMAIN.as_bytes(),
            SEMANTIC_IMPLEMENTATION_VERSION,
            MAX_INTERVAL_NODES,
            EXACT_I53_LIMIT,
            MAX_RESIDUAL_DIMENSION,
            MAX_RESIDUAL_MATRIX_ENTRIES,
            REGISTRY_REVISION
        );
        assert_moves!(
            descriptor,
            SEMANTIC_PLUGIN_IDENTITY_DOMAIN.as_bytes(),
            SEMANTIC_IMPLEMENTATION_VERSION,
            MAX_INTERVAL_NODES + 1,
            EXACT_I53_LIMIT,
            MAX_RESIDUAL_DIMENSION,
            MAX_RESIDUAL_MATRIX_ENTRIES,
            REGISTRY_REVISION
        );
        assert_moves!(
            descriptor,
            SEMANTIC_PLUGIN_IDENTITY_DOMAIN.as_bytes(),
            SEMANTIC_IMPLEMENTATION_VERSION,
            MAX_INTERVAL_NODES,
            EXACT_I53_LIMIT - 1,
            MAX_RESIDUAL_DIMENSION,
            MAX_RESIDUAL_MATRIX_ENTRIES,
            REGISTRY_REVISION
        );
        assert_moves!(
            descriptor,
            SEMANTIC_PLUGIN_IDENTITY_DOMAIN.as_bytes(),
            SEMANTIC_IMPLEMENTATION_VERSION,
            MAX_INTERVAL_NODES,
            EXACT_I53_LIMIT,
            MAX_RESIDUAL_DIMENSION,
            MAX_RESIDUAL_MATRIX_ENTRIES,
            b"mutated-registry-revision"
        );

        let residual = PLUGINS[1];
        let residual_baseline = residual.fingerprint();
        assert_ne!(
            semantic_plugin_fingerprint_with_schema(
                &residual,
                SEMANTIC_PLUGIN_IDENTITY_DOMAIN.as_bytes(),
                SEMANTIC_IMPLEMENTATION_VERSION,
                MAX_INTERVAL_NODES,
                EXACT_I53_LIMIT,
                MAX_RESIDUAL_DIMENSION + 1,
                MAX_RESIDUAL_MATRIX_ENTRIES,
                REGISTRY_REVISION,
            ),
            residual_baseline
        );
        assert_ne!(
            semantic_plugin_fingerprint_with_schema(
                &residual,
                SEMANTIC_PLUGIN_IDENTITY_DOMAIN.as_bytes(),
                SEMANTIC_IMPLEMENTATION_VERSION,
                MAX_INTERVAL_NODES,
                EXACT_I53_LIMIT,
                MAX_RESIDUAL_DIMENSION,
                MAX_RESIDUAL_MATRIX_ENTRIES + 1,
                REGISTRY_REVISION,
            ),
            residual_baseline
        );
    }

    #[test]
    fn semantic_registry_identity_fields_move_independently() {
        let schema = semantic_registry_schema();
        let baseline = semantic_registry_fingerprint_with_schema(
            &schema,
            SEMANTIC_REGISTRY_IDENTITY_DOMAIN.as_bytes(),
        );
        macro_rules! assert_moves {
            ($changed:expr) => {
                assert_ne!(
                    semantic_registry_fingerprint_with_schema(
                        &$changed,
                        SEMANTIC_REGISTRY_IDENTITY_DOMAIN.as_bytes(),
                    ),
                    baseline
                );
            };
        }
        assert_ne!(
            semantic_registry_fingerprint_with_schema(
                &schema,
                b"fs-checker:semantic-registry:v2-mutated",
            ),
            baseline
        );
        assert_moves!(SemanticRegistrySchema {
            revision: b"mutated-registry-revision",
            ..schema
        });
        assert_moves!(SemanticRegistrySchema {
            implementation_version: schema.implementation_version + 1,
            ..schema
        });
        assert_moves!(SemanticRegistrySchema {
            maximum_family_bytes: schema.maximum_family_bytes + 1,
            ..schema
        });
        assert_moves!(SemanticRegistrySchema {
            maximum_witness_bytes: schema.maximum_witness_bytes + 1,
            ..schema
        });
        assert_moves!(SemanticRegistrySchema {
            maximum_witnesses: schema.maximum_witnesses + 1,
            ..schema
        });
        assert_moves!(SemanticRegistrySchema {
            maximum_payload_bytes: schema.maximum_payload_bytes + 1,
            ..schema
        });
        assert_moves!(SemanticRegistrySchema {
            maximum_operations: schema.maximum_operations + 1,
            ..schema
        });

        assert_moves!(SemanticRegistrySchema {
            plugins: &ONE_PLUGIN,
            ..schema
        });
        assert_moves!(SemanticRegistrySchema {
            plugins: &REVERSED_PLUGINS,
            ..schema
        });
        assert_moves!(SemanticRegistrySchema {
            plugins: &FAMILY_CHANGED_PLUGINS,
            ..schema
        });
        assert_moves!(SemanticRegistrySchema {
            plugins: &SCHEMA_CHANGED_PLUGINS,
            ..schema
        });
        assert_moves!(SemanticRegistrySchema {
            plugins: &FINGERPRINT_CHANGED_PLUGINS,
            ..schema
        });
    }

    #[test]
    fn semantic_report_identity_fields_move_independently() {
        let report = semantic_report_fixture();
        let baseline = report.recomputed_context_hash();
        assert_eq!(baseline, report.context_hash());
        macro_rules! assert_moves {
            ($mutation:expr) => {{
                let mut changed = report.clone();
                $mutation(&mut changed);
                assert_ne!(changed.recomputed_context_hash(), baseline);
            }};
        }

        assert_ne!(
            semantic_context_hash_with_schema(&report, b"fs-checker:semantic-report:v2-mutated",),
            baseline
        );
        assert_moves!(|changed: &mut SemanticReport| changed.status = SemanticStatus::NotProvided);
        assert_moves!(|changed: &mut SemanticReport| changed.package_root = ContentHash([31; 32]));
        assert_moves!(
            |changed: &mut SemanticReport| changed.registry_fingerprint = ContentHash([30; 32])
        );
        assert_moves!(|changed: &mut SemanticReport| {
            let duplicate = changed.claims[0].clone();
            changed.claims.push(duplicate);
        });
        assert_moves!(|changed: &mut SemanticReport| changed.claims.swap(0, 1));
        assert_moves!(|changed: &mut SemanticReport| {
            let duplicate = changed.failures[0].clone();
            changed.failures.push(duplicate);
        });
        assert_moves!(|changed: &mut SemanticReport| changed.failures.swap(0, 1));
        assert_moves!(|changed: &mut SemanticReport| changed.witnesses += 1);
        assert_moves!(|changed: &mut SemanticReport| changed.payload_bytes += 1);
        assert_moves!(|changed: &mut SemanticReport| changed.operations += 1);

        assert_moves!(|changed: &mut SemanticReport| changed.claims[0].claim_index += 1);
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0].claim_id.push('x'));
        assert_moves!(
            |changed: &mut SemanticReport| changed.claims[0].claim_hash = ContentHash([29; 32])
        );
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0].family = None);
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .family
            .as_mut()
            .expect("fixture family")
            .push('x'));
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0].schema_version = None);
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0].schema_version = Some(2));
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0].witness_hash = None);
        assert_moves!(
            |changed: &mut SemanticReport| changed.claims[0].witness_hash =
                Some(ContentHash([28; 32]))
        );
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0].plugin_fingerprint = None);
        assert_moves!(
            |changed: &mut SemanticReport| changed.claims[0].plugin_fingerprint =
                Some(ContentHash([27; 32]))
        );
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0].status =
            SemanticClaimStatus::NotProvided);
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0].operations += 1);
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0].failure = None);

        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .failure
            .as_mut()
            .expect("fixture failure")
            .claim_index = None);
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .failure
            .as_mut()
            .expect("fixture failure")
            .claim_index = Some(4));
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .failure
            .as_mut()
            .expect("fixture failure")
            .claim_id = None);
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .failure
            .as_mut()
            .expect("fixture failure")
            .claim_id
            .as_mut()
            .expect("fixture claim id")
            .push('x'));
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .failure
            .as_mut()
            .expect("fixture failure")
            .family = None);
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .failure
            .as_mut()
            .expect("fixture failure")
            .family
            .as_mut()
            .expect("fixture failure family")
            .push('x'));
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .failure
            .as_mut()
            .expect("fixture failure")
            .schema_version = None);
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .failure
            .as_mut()
            .expect("fixture failure")
            .schema_version = Some(2));
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .failure
            .as_mut()
            .expect("fixture failure")
            .kind =
            SemanticFailureKind::VerifierPanic);
        assert_moves!(|changed: &mut SemanticReport| changed.claims[0]
            .failure
            .as_mut()
            .expect("fixture failure")
            .detail
            .push('x'));

        let mut stored_hash_changed = report.clone();
        stored_hash_changed.context_hash = ContentHash([26; 32]);
        assert!(!stored_hash_changed.validate_context_hash());
    }
}
