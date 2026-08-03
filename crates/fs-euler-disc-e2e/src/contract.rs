//! Immutable Euler-disc scientific contract and canonical claim graph.
//!
//! Generic V&V artifacts, claim identifiers/dependencies, and no-claim sets
//! remain owned by `fs-evidence`, `fs-ir`, and `fs-govern`.  This module binds
//! exact instances of those lower schemas into an Euler-only addendum; it does
//! not fork their transport or mint scientific authority.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

use fs_blake3::ContentHash;
use fs_evidence::vv::{
    ApplicabilityPolicy, ArtifactKind, ContextOfUse, QoiId, VV_ARTIFACT_FAMILY,
    VV_SCHEMA_ADMISSION_RECEIPT_IDENTITY_DOMAIN, VV_SCHEMA_VERSION,
};
use fs_govern::evidence_contract::{AUTHORITY_ALGEBRA_VERSION, NoClaimBoundary};
use fs_ir::campaign::{
    CampaignClaim, CampaignClaimId, ClaimDependency, EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1,
    EvidenceGap, EvidenceGapId, EvidenceUse,
};

/// First Euler scientific-contract schema.
pub const EULER_CONTRACT_SCHEMA_VERSION: u32 = 1;
/// First separately versioned claim-policy schema.
pub const EULER_CLAIM_POLICY_SCHEMA_VERSION: u32 = 1;
/// Domain-separated identity for the complete Euler scientific contract.
pub const EULER_CONTRACT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.scientific-contract.v1";
/// Canonical transport domain for the Euler-only claim graph.
pub const EULER_CLAIM_GRAPH_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.claim-graph.v1";
/// Domain for immutable hypothesis-only source declarations.
pub const HYPOTHESIS_SOURCE_DECLARATION_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.hypothesis-source-declaration.v1";
/// Identity domain for exact evidence packets.
pub const EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.claim-evidence-packet.v1";
/// Identity domain for one direct claim-dependency receipt.
pub const EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.prerequisite-assessment-receipt.v1";
/// Identity domain for a complete local structural assessment.
pub const EULER_ASSESSMENT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.claim-policy-assessment.v1";
/// Domain for exact frozen-contract checker receipts.
pub const CONTRACT_CHECK_RECEIPT_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.contract-check-receipt.v1";
/// Domain for the deterministic retained log of one claim-policy assessment.
///
/// The campaign-wide evidence-event and retained-log namespace is owned by the
/// later evidence-log contract. This deliberately narrow domain cannot be
/// mistaken for that broader protocol.
pub const CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.claim-policy-assessment-log.v1";
/// Routing address reserved for the separately versioned, Context-bound
/// aggregate-QoI derivation receipt owned by the downstream QoI/event/scoring
/// contract.
///
/// The receipt schema itself is implemented by bead t6314.1.3. That receipt
/// must bind the exact Context, detailed-observable registry, design set, claim,
/// and aggregate-QoI scoring scope; its admission checker must cross-check the
/// design set against the consuming packet. This v1 leaf binds only its exact
/// routing address and a nonzero content identity; it cannot interpret, admit,
/// or mint such a receipt.
pub const EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA: &str =
    "org.frankensim.fs-euler-disc-e2e.aggregate-qoi-derivation-receipt.v1";
/// First schema for the leaf owner/role routing registry.
pub const EULER_OWNER_MATRIX_SCHEMA_VERSION: u32 = 1;
/// Domain-separated identity for the owner/role routing registry.
pub const EULER_OWNER_MATRIX_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.owner-matrix.v1";

/// Owner-local identity declaration for hypothesis-only source rows.
pub const HYPOTHESIS_SOURCE_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-euler-disc-e2e:hypothesis-source",
    "version_const=EULER_CONTRACT_SCHEMA_VERSION",
    "version=1",
    "domain=org.frankensim.fs-euler-disc-e2e.hypothesis-source-declaration.v1",
    "domain_const=HYPOTHESIS_SOURCE_DECLARATION_DOMAIN",
    "encoder=hypothesis_source_declaration_hash",
    "encoder_helpers=hypothesis_source_declaration_hash_with_schema",
    "schema_functions=HypothesisSource::try_new,HypothesisSource::from_canonical_parts,HypothesisSource::verify_identity,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_constants=EULER_CONTRACT_SCHEMA_VERSION,HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,MAX_EULER_TEXT_BYTES",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=HypothesisSource",
    "source_fields=HypothesisSource.id:semantic,HypothesisSource.locator:semantic,HypothesisSource.declaration_hash:derived:cached-root-of-exact-declaration",
    "source_bindings=HypothesisSource.id>source-id,HypothesisSource.locator>source-locator",
    "external_semantic_fields=identity-domain,identity-version,contract-schema-version,canonical-field-order,length-framing,fixed-numeric-little-endian",
    "semantic_fields=identity-domain,identity-version,contract-schema-version,canonical-field-order,length-framing,fixed-numeric-little-endian,source-id,source-locator",
    "excluded_fields=none",
    "consumers=HypothesisSource::try_new,HypothesisSource::from_canonical_parts,HypothesisSource::declaration_hash,HypothesisSource::verify_identity,EulerScientificContract::try_new",
    "mutations=identity-domain:crates/fs-euler-disc-e2e/src/contract.rs#hypothesis_source_identity_semantic_fields_move_independently,identity-version:crates/fs-euler-disc-e2e/src/contract.rs#hypothesis_source_identity_semantic_fields_move_independently,contract-schema-version:crates/fs-euler-disc-e2e/src/contract.rs#hypothesis_source_identity_semantic_fields_move_independently,canonical-field-order:crates/fs-euler-disc-e2e/src/contract.rs#hypothesis_source_identity_semantic_fields_move_independently,length-framing:crates/fs-euler-disc-e2e/src/contract.rs#hypothesis_source_identity_semantic_fields_move_independently,fixed-numeric-little-endian:crates/fs-euler-disc-e2e/src/contract.rs#hypothesis_source_identity_semantic_fields_move_independently,source-id:crates/fs-euler-disc-e2e/src/contract.rs#hypothesis_source_identity_semantic_fields_move_independently,source-locator:crates/fs-euler-disc-e2e/src/contract.rs#hypothesis_source_identity_semantic_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_hypothesis_source_identity_fields",
    "transport_guard=HypothesisSource::verify_identity",
    "version_guard=crates/fs-euler-disc-e2e/tests/scientific_contract.rs#euler_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-euler-disc-e2e:hypothesis-source",
];

/// Owner-local identity declaration for the exact nine-claim DAG.
pub const EULER_CLAIM_GRAPH_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-euler-disc-e2e:claim-graph",
    "version_const=EULER_CLAIM_POLICY_SCHEMA_VERSION",
    "version=1",
    "domain=org.frankensim.fs-euler-disc-e2e.claim-graph.v1",
    "domain_const=EULER_CLAIM_GRAPH_IDENTITY_DOMAIN",
    "encoder=EulerClaimGraph::content_hash",
    "encoder_helpers=EulerClaimGraph::canonical_bytes,EulerClaimSpec::encode",
    "schema_functions=EulerClaimGraph::try_new,EulerClaimGraph::canonical_bytes,EulerClaimGraph::from_canonical_bytes,EulerClaimGraph::content_hash,EulerClaimSpec::try_new,EulerClaimSpec::encode,EulerClaimKind::tag,EulerClaimKind::from_tag,EulerClaimKind::id,EulerClaimKind::forbids_target_fitting,EulerClaimKind::acceptance_family,EulerClaimKind::required_qoi_ids,EulerClaimKind::required_evidence,EulerAcceptanceFamily::tag,EvidenceRequirement::tag,EvidenceRequirement::from_tag,CanonicalWriter::u8,CanonicalWriter::u32,CanonicalWriter::usize,CanonicalWriter::blob,CanonicalWriter::string,CanonicalWriter::strings,CanonicalReader::new,CanonicalReader::take,CanonicalReader::u8,CanonicalReader::u32,CanonicalReader::usize,CanonicalReader::bounded_len,CanonicalReader::blob,CanonicalReader::hash,CanonicalReader::string,CanonicalReader::is_finished,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-evidence/src/vv/model.rs#vv_id,crates/fs-ir/src/campaign.rs#campaign_id",
    "schema_constants=EULER_CLAIM_POLICY_SCHEMA_VERSION,EULER_CLAIM_GRAPH_IDENTITY_DOMAIN,GRAPH_MAGIC,EULER_CLAIM_REGISTRY,EULER_EVIDENCE_REQUIREMENT_REGISTRY,MAX_EULER_CLAIMS,MAX_EULER_GRAPH_BYTES,MAX_EULER_TEXT_BYTES",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=EulerClaimGraph,EulerClaimSpec",
    "source_fields=EulerClaimGraph.claims:semantic,EulerClaimGraph.dependencies:semantic,EulerClaimSpec.kind:derived:transitively-bound-by-claim-registry,EulerClaimSpec.campaign:derived:transitively-bound-by-claim-registry,EulerClaimSpec.requirements:derived:transitively-bound-by-claim-registry",
    "source_bindings=EulerClaimGraph.claims>claim-registry,EulerClaimGraph.dependencies>dependency-registry",
    "external_semantic_fields=identity-domain,identity-version,claim-policy-schema-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian",
    "semantic_fields=identity-domain,identity-version,claim-policy-schema-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian,claim-registry,dependency-registry",
    "excluded_fields=none",
    "consumers=EulerClaimGraph::content_hash,EulerScientificContract::try_new,EulerScientificContract::canonical_bytes,ContractCheckReceipt::new,ContractCheckReceipt::verify_subject,check_frozen_contract",
    "mutations=identity-domain:crates/fs-euler-disc-e2e/src/contract.rs#claim_graph_identity_semantic_fields_move_independently,identity-version:crates/fs-euler-disc-e2e/src/contract.rs#claim_graph_identity_semantic_fields_move_independently,claim-policy-schema-version:crates/fs-euler-disc-e2e/src/contract.rs#claim_graph_identity_semantic_fields_move_independently,transport-magic:crates/fs-euler-disc-e2e/src/contract.rs#claim_graph_identity_semantic_fields_move_independently,canonical-field-order:crates/fs-euler-disc-e2e/src/contract.rs#claim_graph_identity_semantic_fields_move_independently,length-framing:crates/fs-euler-disc-e2e/src/contract.rs#claim_graph_identity_semantic_fields_move_independently,fixed-numeric-little-endian:crates/fs-euler-disc-e2e/src/contract.rs#claim_graph_identity_semantic_fields_move_independently,claim-registry:crates/fs-euler-disc-e2e/src/contract.rs#claim_graph_identity_semantic_fields_move_independently,dependency-registry:crates/fs-euler-disc-e2e/src/contract.rs#claim_graph_identity_semantic_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_euler_claim_graph_identity_fields",
    "transport_guard=EulerClaimGraph::from_canonical_bytes",
    "version_guard=crates/fs-euler-disc-e2e/tests/scientific_contract.rs#euler_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-euler-disc-e2e:claim-graph",
];

/// Owner-local identity declaration for the leaf owner/role routing registry.
///
/// `OwnerRow::source_schema` values are opaque routing addresses owned by this
/// registry. They deliberately do not create identity-dependency edges to the
/// schemas at those addresses; consumers resolve role meaning through this
/// separately versioned registry instead.
pub const OWNER_MATRIX_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-euler-disc-e2e:owner-matrix",
    "version_const=EULER_OWNER_MATRIX_SCHEMA_VERSION",
    "version=1",
    "domain=org.frankensim.fs-euler-disc-e2e.owner-matrix.v1",
    "domain_const=EULER_OWNER_MATRIX_IDENTITY_DOMAIN",
    "encoder=OwnerMatrix::identity",
    "encoder_helpers=OwnerMatrix::canonical_bytes,encode_owner_matrix_components",
    "schema_functions=OwnerMatrix::try_new,OwnerMatrix::canonical_bytes,OwnerMatrix::from_canonical_bytes,OwnerMatrix::identity,encode_owner_matrix_components,OwnerRole::tag,OwnerRole::from_tag,OwnerRole::expected_owner_crate,OwnerRole::expected_source_schema,OwnerRole::expected_authority_ceiling,AuthorityCeiling::tag,AuthorityCeiling::from_tag,OwnerRow::try_new,CanonicalWriter::u8,CanonicalWriter::u32,CanonicalWriter::usize,CanonicalWriter::blob,CanonicalWriter::string,CanonicalWriter::strings,CanonicalReader::new,CanonicalReader::take,CanonicalReader::u8,CanonicalReader::u32,CanonicalReader::usize,CanonicalReader::bounded_len,CanonicalReader::blob,CanonicalReader::hash,CanonicalReader::string,CanonicalReader::is_finished,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_constants=EULER_OWNER_MATRIX_SCHEMA_VERSION,EULER_OWNER_MATRIX_IDENTITY_DOMAIN,OWNER_MATRIX_MAGIC,EULER_OWNER_ROLE_REGISTRY,MAX_OWNER_MATRIX_BYTES,MAX_EULER_TEXT_BYTES,crates/fs-evidence/src/vv/model.rs#VV_ARTIFACT_FAMILY,crates/fs-evidence/src/vv/model.rs#VV_SCHEMA_ADMISSION_RECEIPT_IDENTITY_DOMAIN,FS_IR_CAMPAIGN_SOURCE_SCHEMA,FS_GOVERN_AUTHORITY_SOURCE_SCHEMA,HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,EULER_CLAIM_GRAPH_IDENTITY_DOMAIN,EULER_CONTRACT_IDENTITY_DOMAIN,EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN,EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,EULER_ASSESSMENT_IDENTITY_DOMAIN,CONTRACT_CHECK_RECEIPT_DOMAIN,CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN,EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=OwnerMatrix,OwnerRow",
    "source_fields=OwnerMatrix.schema_version:semantic,OwnerMatrix.rows:semantic,OwnerMatrix.identity:derived:cached-root-of-exact-registry,OwnerRow.role:derived:transitively-bound-by-owner-registry,OwnerRow.owner_crate:derived:transitively-bound-by-owner-registry,OwnerRow.source_schema:derived:opaque-routing-address-bound-by-owner-registry,OwnerRow.authority_ceiling:derived:transitively-bound-by-owner-registry",
    "source_bindings=OwnerMatrix.schema_version>owner-matrix-schema-version,OwnerMatrix.rows>owner-role-registry",
    "external_semantic_fields=identity-domain,identity-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian,owner-role-tags,authority-ceiling-tags",
    "semantic_fields=identity-domain,identity-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian,owner-role-tags,authority-ceiling-tags,owner-matrix-schema-version,owner-role-registry",
    "excluded_fields=none",
    "consumers=OwnerMatrix::identity,OwnerMatrix::canonical_bytes,OwnerMatrix::from_canonical_bytes,EulerScientificContract::try_new,EulerScientificContract::canonical_bytes,crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_log",
    "mutations=identity-domain:crates/fs-euler-disc-e2e/src/contract.rs#owner_matrix_identity_semantic_fields_move_independently,identity-version:crates/fs-euler-disc-e2e/src/contract.rs#owner_matrix_identity_semantic_fields_move_independently,transport-magic:crates/fs-euler-disc-e2e/src/contract.rs#owner_matrix_identity_semantic_fields_move_independently,canonical-field-order:crates/fs-euler-disc-e2e/src/contract.rs#owner_matrix_identity_semantic_fields_move_independently,length-framing:crates/fs-euler-disc-e2e/src/contract.rs#owner_matrix_identity_semantic_fields_move_independently,fixed-numeric-little-endian:crates/fs-euler-disc-e2e/src/contract.rs#owner_matrix_identity_semantic_fields_move_independently,owner-role-tags:crates/fs-euler-disc-e2e/src/contract.rs#owner_matrix_identity_semantic_fields_move_independently,authority-ceiling-tags:crates/fs-euler-disc-e2e/src/contract.rs#owner_matrix_identity_semantic_fields_move_independently,owner-matrix-schema-version:crates/fs-euler-disc-e2e/src/contract.rs#owner_matrix_identity_semantic_fields_move_independently,owner-role-registry:crates/fs-euler-disc-e2e/src/contract.rs#owner_matrix_identity_semantic_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_owner_matrix_identity_fields",
    "transport_guard=OwnerMatrix::from_canonical_bytes",
    "version_guard=crates/fs-euler-disc-e2e/tests/scientific_contract.rs#owner_matrix_transport_is_exact_versioned_and_fail_closed",
    "coupling_surface=fs-euler-disc-e2e:owner-matrix",
];

/// Owner-local identity declaration for the complete frozen scientific
/// contract transport.
pub const EULER_SCIENTIFIC_CONTRACT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-euler-disc-e2e:scientific-contract",
    "version_const=EULER_CONTRACT_SCHEMA_VERSION",
    "version=1",
    "domain=org.frankensim.fs-euler-disc-e2e.scientific-contract.v1",
    "domain_const=EULER_CONTRACT_IDENTITY_DOMAIN",
    "encoder=scientific_contract_identity",
    "encoder_helpers=encode_contract_components,EulerContextExtension::encode,EulerClaimGraph::canonical_bytes,OwnerMatrix::canonical_bytes",
    "schema_functions=EulerScientificContract::try_new,EulerScientificContract::canonical_bytes,EulerScientificContract::from_canonical_bytes,scientific_contract_identity,encode_contract_components,EulerContextExtension::try_new,EulerContextExtension::observation_frame,EulerContextExtension::encode,ScientificRisk::try_new,HypothesisSource::from_canonical_parts,EulerClaimGraph::canonical_bytes,EulerClaimGraph::from_canonical_bytes,OwnerMatrix::canonical_bytes,OwnerMatrix::from_canonical_bytes,OwnerMatrix::identity,OwnerMatrixIdentity::as_hash,EulerClaimKind::from_tag,CanonicalWriter::u8,CanonicalWriter::u32,CanonicalWriter::usize,CanonicalWriter::blob,CanonicalWriter::string,CanonicalWriter::strings,CanonicalReader::new,CanonicalReader::take,CanonicalReader::u8,CanonicalReader::u32,CanonicalReader::usize,CanonicalReader::bounded_len,CanonicalReader::blob,CanonicalReader::hash,CanonicalReader::string,CanonicalReader::is_finished,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes,crates/fs-evidence/src/vv/model.rs#ContextOfUse::header,crates/fs-evidence/src/vv/model.rs#ContextOfUse::qois,crates/fs-evidence/src/vv/model.rs#ContextOfUse::applicability,crates/fs-evidence/src/vv/model.rs#ContextOfUse::applicability_policy,crates/fs-evidence/src/vv/model.rs#ArtifactHeader::units,crates/fs-evidence/src/vv/model.rs#ApplicabilityDomain::numeric,crates/fs-evidence/src/vv/model.rs#ApplicabilityDomain::categorical,crates/fs-evidence/src/vv/model.rs#NumericDomainAxis::unit,crates/fs-evidence/src/vv/model.rs#CategoricalDomainAxis::allowed,crates/fs-evidence/src/vv/model.rs#QoiSpec::id,crates/fs-evidence/src/vv/model.rs#QoiSpec::unit,crates/fs-evidence/src/vv/model.rs#vv_id,crates/fs-evidence/src/vv/codec.rs#canonical_artifact_bytes,crates/fs-evidence/src/vv/codec.rs#VvArtifact::from_canonical_bytes,crates/fs-evidence/src/vv/codec.rs#content_hash_for,crates/fs-evidence/src/vv/codec.rs#encode_context,crates/fs-evidence/src/vv/codec.rs#decode_context,crates/fs-govern/src/evidence_contract.rs#NoClaimBoundary::new,crates/fs-govern/src/evidence_contract.rs#NoClaimBoundary::entries",
    "schema_constants=EULER_CONTRACT_SCHEMA_VERSION,EULER_CONTRACT_IDENTITY_DOMAIN,CONTRACT_MAGIC,MAX_EULER_CONTRACT_BYTES,MAX_EULER_NO_CLAIMS,MAX_EULER_TEXT_BYTES,CORE_NO_CLAIMS,crates/fs-evidence/src/vv/model.rs#VV_SCHEMA_VERSION,crates/fs-evidence/src/vv/model.rs#VV_ARTIFACT_FAMILY,crates/fs-ir/src/campaign.rs#EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1,crates/fs-govern/src/evidence_contract.rs#AUTHORITY_ALGEBRA_VERSION",
    "schema_dependencies=fs-euler-disc-e2e:claim-graph,fs-euler-disc-e2e:hypothesis-source,fs-euler-disc-e2e:owner-matrix,fs-evidence:vv-artifact",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=EulerScientificContract,EulerContextExtension,ScientificRisk",
    "source_fields=EulerScientificContract.schema_version:semantic,EulerScientificContract.context:derived:canonical-decoded-view-of-context-bytes,EulerScientificContract.context_bytes:semantic,EulerScientificContract.context_hash:derived:native-context-artifact-hash,EulerScientificContract.extension:derived:separately-classified-extension,EulerScientificContract.claim_graph:semantic,EulerScientificContract.no_claims:semantic,EulerScientificContract.owner_matrix:semantic,EulerScientificContract.identity:derived:cached-root-of-complete-contract,EulerContextExtension.users:semantic,EulerContextExtension.apparatus_population:semantic,EulerContextExtension.environment_population:semantic,EulerContextExtension.observation_frame:semantic,EulerContextExtension.decision_alternatives:semantic,EulerContextExtension.risks:semantic,EulerContextExtension.hypothesis_sources:semantic,ScientificRisk.code:derived:transitively-bound-by-risk-registry,ScientificRisk.consequence:derived:transitively-bound-by-risk-registry,ScientificRisk.severity:derived:transitively-bound-by-risk-registry,ScientificRisk.affected_claims:derived:transitively-bound-by-risk-registry,ScientificRisk.decision_alternative:derived:transitively-bound-by-risk-registry",
    "source_bindings=EulerScientificContract.schema_version>contract-schema-version,EulerScientificContract.context_bytes>context-canonical-bytes,EulerScientificContract.claim_graph>claim-graph-canonical-bytes,EulerScientificContract.no_claims>no-claim-boundary,EulerScientificContract.owner_matrix>owner-matrix,EulerContextExtension.users>extension-users,EulerContextExtension.apparatus_population>apparatus-population,EulerContextExtension.environment_population>environment-population,EulerContextExtension.observation_frame>observation-frame,EulerContextExtension.decision_alternatives>decision-alternatives,EulerContextExtension.risks>risk-registry,EulerContextExtension.hypothesis_sources>hypothesis-source-declarations",
    "external_semantic_fields=identity-domain,identity-version,transport-magic,embedded-vv-schema-version,embedded-vv-artifact-family,embedded-campaign-schema-version,embedded-authority-version,canonical-field-order,length-framing,fixed-numeric-little-endian",
    "semantic_fields=identity-domain,identity-version,transport-magic,embedded-vv-schema-version,embedded-vv-artifact-family,embedded-campaign-schema-version,embedded-authority-version,canonical-field-order,length-framing,fixed-numeric-little-endian,contract-schema-version,context-canonical-bytes,claim-graph-canonical-bytes,no-claim-boundary,extension-users,apparatus-population,environment-population,observation-frame,decision-alternatives,risk-registry,hypothesis-source-declarations,owner-matrix",
    "excluded_fields=none",
    "consumers=EulerScientificContract::identity,EulerScientificContract::canonical_bytes,EulerScientificContract::from_canonical_bytes,EvidenceRecord::try_new,ClaimEvidencePacket::try_new,check_frozen_contract,admit_frozen_contract,crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_log,ContractCheckReceipt::verify_subject",
    "mutations=identity-domain:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,identity-version:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,transport-magic:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,embedded-vv-schema-version:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,embedded-vv-artifact-family:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,embedded-campaign-schema-version:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,embedded-authority-version:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,canonical-field-order:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,length-framing:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,fixed-numeric-little-endian:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,contract-schema-version:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,context-canonical-bytes:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,claim-graph-canonical-bytes:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,no-claim-boundary:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,extension-users:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,apparatus-population:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,environment-population:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,observation-frame:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,decision-alternatives:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,risk-registry:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,hypothesis-source-declarations:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently,owner-matrix:crates/fs-euler-disc-e2e/src/contract.rs#scientific_contract_identity_semantic_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_euler_contract_identity_fields",
    "transport_guard=EulerScientificContract::from_canonical_bytes",
    "version_guard=crates/fs-euler-disc-e2e/tests/scientific_contract.rs#euler_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-euler-disc-e2e:scientific-contract",
];
/// Exact generic campaign vocabulary schema bound into owner rows.
pub const FS_IR_CAMPAIGN_SOURCE_SCHEMA: &str = "fs-ir:experiment-campaign-schema-v1";
/// Exact generic no-claim authority schema bound into owner rows.
pub const FS_GOVERN_AUTHORITY_SOURCE_SCHEMA: &str = "fs-govern:authority-algebra-v2";
/// The closed v1 taxonomy contains exactly nine claim kinds.
pub const MAX_EULER_CLAIMS: usize = 9;
/// Maximum no-claim rows accepted by the public contract and its decoder.
pub const MAX_EULER_NO_CLAIMS: usize = MAX_EULER_CLAIMS * 8;
/// Bound on one Euler-local descriptive field.
pub const MAX_EULER_TEXT_BYTES: usize = 4_096;
/// Bound on one canonical claim graph.
pub const MAX_EULER_GRAPH_BYTES: usize = 1024 * 1024;
/// Bound on the complete contract transport, including the generic context.
pub const MAX_EULER_CONTRACT_BYTES: usize = 8 * 1024 * 1024;
/// Bound on the complete owner/role routing-registry transport.
pub const MAX_OWNER_MATRIX_BYTES: usize = 128 * 1024;

/// Exact v1 no-claim set. These statements constrain every local assessment.
pub const CORE_NO_CLAIMS: [&str; 7] = [
    "Transcript and publication sources generate hypotheses only; they are not validation evidence.",
    "Fitting or selecting against protected target outcomes is calibrated reproduction, not emergent prediction.",
    "Agreement in an exponent, event time, or stop time does not identify an energy-loss mechanism.",
    "Geometric similarity does not establish dynamic similarity across scale, material, support, or environment.",
    "Deterministic software verification does not establish physical validation.",
    "A successful blind case is local to its exact declared Context of Use and applicability domain.",
    "Negative and inconclusive results are retained terminal outcomes and are never erased or promoted.",
];

const GRAPH_MAGIC: &[u8; 8] = b"FSEDGR01";
const CONTRACT_MAGIC: &[u8; 8] = b"FSEDSC01";
const OWNER_MATRIX_MAGIC: &[u8; 8] = b"FSEDOM01";

// These assertions deliberately turn an upstream generic-schema move into a
// compile-time integration failure.  The human-readable owner schema strings
// must never drift independently from the numeric versions encoded below.
const _: () = assert!(EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1 == 1);
const _: () = assert!(AUTHORITY_ALGEBRA_VERSION == 2);

/// A deterministic contract construction or checking refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContractError {
    code: &'static str,
    detail: String,
}

impl ContractError {
    pub(crate) fn new(code: &'static str, detail: impl Into<String>) -> Self {
        Self {
            code,
            detail: detail.into(),
        }
    }

    /// Stable machine code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.code
    }

    /// Actionable bounded detail.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for ContractError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {}", self.code, self.detail)
    }
}

impl std::error::Error for ContractError {}

fn checked_text(field: &'static str, value: impl Into<String>) -> Result<String, ContractError> {
    let value = value.into();
    // Apply the byte ceiling before any full-input Unicode/whitespace scan so
    // hostile metadata cannot bypass the advertised bounded-work boundary.
    if value.len() > MAX_EULER_TEXT_BYTES
        || value.trim().is_empty()
        || value.trim() != value
        || value.chars().any(char::is_control)
    {
        return Err(ContractError::new(
            "EulerContractInvalidText",
            format!("{field} must be nonblank, trimmed, control-free, and bounded"),
        ));
    }
    Ok(value)
}

fn canonical_text_set(
    field: &'static str,
    values: Vec<String>,
    allow_empty: bool,
) -> Result<Vec<String>, ContractError> {
    if (!allow_empty && values.is_empty()) || values.len() > MAX_EULER_CLAIMS * 8 {
        return Err(ContractError::new(
            "EulerContractCardinality",
            format!("{field} has an invalid item count"),
        ));
    }
    let original_len = values.len();
    let mut values = values
        .into_iter()
        .map(|value| checked_text(field, value))
        .collect::<Result<Vec<_>, _>>()?;
    values.sort();
    values.dedup();
    if values.len() != original_len {
        return Err(ContractError::new(
            "EulerContractDuplicate",
            format!("{field} contains a duplicate value"),
        ));
    }
    Ok(values)
}

/// Content identity of one complete Euler scientific contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ContractIdentity(ContentHash);

impl ContractIdentity {
    /// Wraps an already-derived contract hash without granting admission.
    ///
    /// This constructor is intentionally structural: callers that receive an
    /// identity independently of an [`EulerScientificContract`] can preserve
    /// and compare it, but the value does not prove that the corresponding
    /// contract exists, is canonical, or matches the frozen v1 contract.
    #[must_use]
    pub const fn from_hash(hash: ContentHash) -> Self {
        Self(hash)
    }

    /// Raw domain-separated content hash.
    #[must_use]
    pub const fn as_hash(self) -> ContentHash {
        self.0
    }
}

impl fmt::Display for ContractIdentity {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

/// Closed Euler-disc claim taxonomy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EulerClaimKind {
    /// Deterministic code/solution verification of a trajectory computation.
    NumericalTrajectoryVerification,
    /// Reproduction after fitting against declared calibration data.
    CalibratedReproduction,
    /// Prediction of a held-out trajectory without target fitting.
    BlindTrajectoryPrediction,
    /// Prediction of terminal events or mechanism crossovers.
    EventOrCrossoverPrediction,
    /// Prediction of a preregistered qualitative effect direction.
    QualitativeEffectDirection,
    /// Ranking of preregistered specimen/configuration alternatives.
    Ranking,
    /// A bounded set or interval containing a nonlinear optimum.
    NonlinearOptimumInterval,
    /// Attribution of work or energy to named, nonoverlapping channels.
    EnergyChannelAttribution,
    /// Discrimination among rival physical mechanisms.
    MechanismAttribution,
}

/// Single executable registry for the closed v1 claim taxonomy.
pub(crate) const EULER_CLAIM_REGISTRY: [EulerClaimKind; MAX_EULER_CLAIMS] = [
    EulerClaimKind::NumericalTrajectoryVerification,
    EulerClaimKind::CalibratedReproduction,
    EulerClaimKind::BlindTrajectoryPrediction,
    EulerClaimKind::EventOrCrossoverPrediction,
    EulerClaimKind::QualitativeEffectDirection,
    EulerClaimKind::Ranking,
    EulerClaimKind::NonlinearOptimumInterval,
    EulerClaimKind::EnergyChannelAttribution,
    EulerClaimKind::MechanismAttribution,
];

impl EulerClaimKind {
    /// Every v1 kind in canonical order.
    pub const ALL: [Self; MAX_EULER_CLAIMS] = EULER_CLAIM_REGISTRY;

    #[must_use]
    /// Stable wire tag for the v1 claim-kind schema.
    pub const fn tag(self) -> u8 {
        match self {
            Self::NumericalTrajectoryVerification => 1,
            Self::CalibratedReproduction => 2,
            Self::BlindTrajectoryPrediction => 3,
            Self::EventOrCrossoverPrediction => 4,
            Self::QualitativeEffectDirection => 5,
            Self::Ranking => 6,
            Self::NonlinearOptimumInterval => 7,
            Self::EnergyChannelAttribution => 8,
            Self::MechanismAttribution => 9,
        }
    }

    pub(crate) fn from_tag(tag: u8) -> Result<Self, ContractError> {
        EULER_CLAIM_REGISTRY
            .into_iter()
            .find(|kind| kind.tag() == tag)
            .ok_or_else(|| {
                ContractError::new(
                    "EulerContractUnknownClaimKind",
                    format!("unknown claim-kind tag {tag}"),
                )
            })
    }

    /// Stable generic campaign claim id.
    #[must_use]
    pub const fn id(self) -> &'static str {
        match self {
            Self::NumericalTrajectoryVerification => "numerical-trajectory-verification",
            Self::CalibratedReproduction => "calibrated-reproduction",
            Self::BlindTrajectoryPrediction => "blind-trajectory-prediction",
            Self::EventOrCrossoverPrediction => "event-or-crossover-prediction",
            Self::QualitativeEffectDirection => "qualitative-effect-direction",
            Self::Ranking => "ranking",
            Self::NonlinearOptimumInterval => "nonlinear-optimum-interval",
            Self::EnergyChannelAttribution => "energy-channel-attribution",
            Self::MechanismAttribution => "mechanism-attribution",
        }
    }

    /// True when target fitting would invalidate the requested emergent claim.
    #[must_use]
    pub const fn forbids_target_fitting(self) -> bool {
        !matches!(
            self,
            Self::NumericalTrajectoryVerification | Self::CalibratedReproduction
        )
    }

    /// Separately versioned acceptance family for this claim kind.
    #[must_use]
    pub const fn acceptance_family(self) -> EulerAcceptanceFamily {
        match self {
            Self::NumericalTrajectoryVerification => {
                EulerAcceptanceFamily::NumericalTrajectoryBound
            }
            Self::CalibratedReproduction => EulerAcceptanceFamily::CalibratedTrajectoryBound,
            Self::BlindTrajectoryPrediction => EulerAcceptanceFamily::BlindTrajectoryBound,
            Self::EventOrCrossoverPrediction => EulerAcceptanceFamily::EventOrCrossoverBound,
            Self::QualitativeEffectDirection => EulerAcceptanceFamily::QualitativeDirection,
            Self::Ranking => EulerAcceptanceFamily::RankingWithTieRule,
            Self::NonlinearOptimumInterval => EulerAcceptanceFamily::BoundedOptimumInterval,
            Self::EnergyChannelAttribution => EulerAcceptanceFamily::EnergyClosureAndAllocation,
            Self::MechanismAttribution => EulerAcceptanceFamily::RivalMechanismDiscrimination,
        }
    }

    /// Exact Context QoI identities that the v1 policy permits for this claim.
    #[must_use]
    pub const fn required_qoi_ids(self) -> &'static [&'static str] {
        match self {
            Self::NumericalTrajectoryVerification => &["numerical-trajectory-error"],
            Self::CalibratedReproduction | Self::BlindTrajectoryPrediction => {
                &["normalized-trajectory-discrepancy"]
            }
            Self::EventOrCrossoverPrediction => &["event-class-disposition", "event-time-error"],
            Self::QualitativeEffectDirection => &["qualitative-effect-disposition"],
            Self::Ranking => &["configuration-ranking-disposition"],
            Self::NonlinearOptimumInterval => {
                &["optimum-containment-disposition", "optimum-interval-width"]
            }
            Self::EnergyChannelAttribution => {
                &["energy-balance-residual", "energy-channel-fraction-error"]
            }
            Self::MechanismAttribution => &[
                "energy-channel-fraction-error",
                "rival-mechanism-disposition",
            ],
        }
    }

    /// Exact minimum evidence roles for this claim under policy v1.
    #[must_use]
    pub const fn required_evidence(self) -> &'static [EvidenceRequirement] {
        use EvidenceRequirement as E;
        match self {
            Self::NumericalTrajectoryVerification => &[
                E::CodeVerification,
                E::SolutionVerification,
                E::IndependentReconstruction,
            ],
            Self::CalibratedReproduction => &[
                E::CodeVerification,
                E::SolutionVerification,
                E::CalibrationPartition,
                E::ApplicabilityCheck,
                E::UncertaintyClosure,
                E::IndependentReconstruction,
            ],
            Self::BlindTrajectoryPrediction
            | Self::EventOrCrossoverPrediction
            | Self::QualitativeEffectDirection
            | Self::Ranking
            | Self::NonlinearOptimumInterval => &[
                E::CodeVerification,
                E::SolutionVerification,
                E::PhysicalValidation,
                E::BlindHoldout,
                E::PreregisteredAnalysis,
                E::ApplicabilityCheck,
                E::UncertaintyClosure,
                E::MultiplicityControl,
                E::IndependentReconstruction,
            ],
            Self::EnergyChannelAttribution => &[
                E::CodeVerification,
                E::SolutionVerification,
                E::PhysicalValidation,
                E::BlindHoldout,
                E::PreregisteredAnalysis,
                E::ApplicabilityCheck,
                E::UncertaintyClosure,
                E::MultiplicityControl,
                E::EnergyBalanceClosure,
                E::IndependentReconstruction,
            ],
            Self::MechanismAttribution => &[
                E::CodeVerification,
                E::SolutionVerification,
                E::PhysicalValidation,
                E::BlindHoldout,
                E::PreregisteredAnalysis,
                E::ApplicabilityCheck,
                E::UncertaintyClosure,
                E::MultiplicityControl,
                E::EnergyBalanceClosure,
                E::IndependentReconstruction,
                E::RivalMechanismDiscrimination,
            ],
        }
    }
}

/// Closed semantic family interpreted by later QoI/event protocol artifacts.
///
/// This enum distinguishes decision meaning now without duplicating the exact
/// event, scoring, censoring, or measurement schemas owned by downstream
/// Euler beads.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EulerAcceptanceFamily {
    /// Independent numerical trajectory-discrepancy bound.
    NumericalTrajectoryBound,
    /// Trajectory bound evaluated on declared calibration data.
    CalibratedTrajectoryBound,
    /// Trajectory bound evaluated on a protected blind holdout.
    BlindTrajectoryBound,
    /// Terminal-event or crossover-time bound.
    EventOrCrossoverBound,
    /// Preregistered direction-of-effect decision.
    QualitativeDirection,
    /// Preregistered ordering with an explicit tie rule.
    RankingWithTieRule,
    /// Bounded interval containing a declared nonlinear optimum.
    BoundedOptimumInterval,
    /// Closed work-energy balance and nonoverlapping channel allocation.
    EnergyClosureAndAllocation,
    /// Discrimination among preregistered rival mechanisms.
    RivalMechanismDiscrimination,
}

impl EulerAcceptanceFamily {
    #[must_use]
    /// Stable wire tag for the v1 acceptance-family schema.
    pub const fn tag(self) -> u8 {
        match self {
            Self::NumericalTrajectoryBound => 1,
            Self::CalibratedTrajectoryBound => 2,
            Self::BlindTrajectoryBound => 3,
            Self::EventOrCrossoverBound => 4,
            Self::QualitativeDirection => 5,
            Self::RankingWithTieRule => 6,
            Self::BoundedOptimumInterval => 7,
            Self::EnergyClosureAndAllocation => 8,
            Self::RivalMechanismDiscrimination => 9,
        }
    }
}

/// Evidence roles required by an Euler claim policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EvidenceRequirement {
    /// Code-verification evidence independent of the target trajectory.
    CodeVerification,
    /// Mesh/time/nonlinear/iterative solution-verification evidence.
    SolutionVerification,
    /// Calibration-only physical data under a sealed partition.
    CalibrationPartition,
    /// External physical validation in the declared Context of Use.
    PhysicalValidation,
    /// Sealed blind-holdout prediction evidence.
    BlindHoldout,
    /// Analysis frozen before protected outcomes were accessible.
    PreregisteredAnalysis,
    /// Exact applicability-domain evaluation for the assessed run.
    ApplicabilityCheck,
    /// Explicit uncertainty-source closure without duplicate accounting.
    UncertaintyClosure,
    /// Multiplicity/censoring/tie handling fixed before scoring.
    MultiplicityControl,
    /// Work-energy closure over uniquely owned channels.
    EnergyBalanceClosure,
    /// Artifact-only reconstruction by a distinct checker implementation.
    IndependentReconstruction,
    /// Discriminating evidence among preregistered rival mechanisms.
    RivalMechanismDiscrimination,
}

/// Single executable registry for the closed v1 evidence-role taxonomy.
pub(crate) const EULER_EVIDENCE_REQUIREMENT_REGISTRY: [EvidenceRequirement; 12] = [
    EvidenceRequirement::CodeVerification,
    EvidenceRequirement::SolutionVerification,
    EvidenceRequirement::CalibrationPartition,
    EvidenceRequirement::PhysicalValidation,
    EvidenceRequirement::BlindHoldout,
    EvidenceRequirement::PreregisteredAnalysis,
    EvidenceRequirement::ApplicabilityCheck,
    EvidenceRequirement::UncertaintyClosure,
    EvidenceRequirement::MultiplicityControl,
    EvidenceRequirement::EnergyBalanceClosure,
    EvidenceRequirement::IndependentReconstruction,
    EvidenceRequirement::RivalMechanismDiscrimination,
];

impl EvidenceRequirement {
    /// Every requirement in canonical schema order.
    pub const ALL: [Self; 12] = EULER_EVIDENCE_REQUIREMENT_REGISTRY;

    #[must_use]
    /// Stable wire tag for the v1 evidence-requirement schema.
    pub const fn tag(self) -> u8 {
        match self {
            Self::CodeVerification => 1,
            Self::SolutionVerification => 2,
            Self::CalibrationPartition => 3,
            Self::PhysicalValidation => 4,
            Self::BlindHoldout => 5,
            Self::PreregisteredAnalysis => 6,
            Self::ApplicabilityCheck => 7,
            Self::UncertaintyClosure => 8,
            Self::MultiplicityControl => 9,
            Self::EnergyBalanceClosure => 10,
            Self::IndependentReconstruction => 11,
            Self::RivalMechanismDiscrimination => 12,
        }
    }

    pub(crate) fn from_tag(tag: u8) -> Result<Self, ContractError> {
        EULER_EVIDENCE_REQUIREMENT_REGISTRY
            .into_iter()
            .find(|requirement| requirement.tag() == tag)
            .ok_or_else(|| {
                ContractError::new(
                    "EulerContractUnknownEvidenceRequirement",
                    format!("unknown evidence-requirement tag {tag}"),
                )
            })
    }

    /// Stable machine code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::CodeVerification => "code-verification",
            Self::SolutionVerification => "solution-verification",
            Self::CalibrationPartition => "calibration-partition",
            Self::PhysicalValidation => "physical-validation",
            Self::BlindHoldout => "blind-holdout",
            Self::PreregisteredAnalysis => "preregistered-analysis",
            Self::ApplicabilityCheck => "applicability-check",
            Self::UncertaintyClosure => "uncertainty-closure",
            Self::MultiplicityControl => "multiplicity-control",
            Self::EnergyBalanceClosure => "energy-balance-closure",
            Self::IndependentReconstruction => "independent-reconstruction",
            Self::RivalMechanismDiscrimination => "rival-mechanism-discrimination",
        }
    }
}

/// One typed risk of making the wrong scientific decision.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ScientificRisk {
    code: String,
    consequence: String,
    severity: u8,
    affected_claims: Vec<EulerClaimKind>,
    decision_alternative: String,
}

impl ScientificRisk {
    /// Construct a bounded risk row. Severity is ordinal 1 (low) through 5 (critical).
    pub fn try_new(
        code: impl Into<String>,
        consequence: impl Into<String>,
        severity: u8,
        mut affected_claims: Vec<EulerClaimKind>,
        decision_alternative: impl Into<String>,
    ) -> Result<Self, ContractError> {
        if !(1..=5).contains(&severity) {
            return Err(ContractError::new(
                "EulerContractInvalidRisk",
                "risk severity must be in 1..=5",
            ));
        }
        if affected_claims.is_empty() || affected_claims.len() > MAX_EULER_CLAIMS {
            return Err(ContractError::new(
                "EulerContractInvalidRisk",
                "a risk must name a bounded nonempty affected-claim set",
            ));
        }
        affected_claims.sort_by_key(|kind| kind.tag());
        let original_claim_count = affected_claims.len();
        affected_claims.dedup();
        if affected_claims.len() != original_claim_count {
            return Err(ContractError::new(
                "EulerContractDuplicate",
                "a risk cannot name the same claim twice",
            ));
        }
        Ok(Self {
            code: checked_text("risk.code", code)?,
            consequence: checked_text("risk.consequence", consequence)?,
            severity,
            affected_claims,
            decision_alternative: checked_text("risk.decision_alternative", decision_alternative)?,
        })
    }

    #[must_use]
    /// Stable machine identity of this risk row.
    pub fn code(&self) -> &str {
        &self.code
    }

    #[must_use]
    /// Consequence of making the wrong decision.
    pub fn consequence(&self) -> &str {
        &self.consequence
    }

    #[must_use]
    /// Ordinal severity from one through five.
    pub const fn severity(&self) -> u8 {
        self.severity
    }

    #[must_use]
    /// Claim kinds affected by this risk.
    pub fn affected_claims(&self) -> &[EulerClaimKind] {
        &self.affected_claims
    }

    #[must_use]
    /// Decision alternative associated with this risk response.
    pub fn decision_alternative(&self) -> &str {
        &self.decision_alternative
    }
}

/// A declaration-bound source locator that may generate hypotheses but cannot
/// satisfy an evidence requirement in this contract.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct HypothesisSource {
    id: String,
    locator: String,
    declaration_hash: ContentHash,
}

impl HypothesisSource {
    /// Construct a hypothesis-only source declaration and derive its identity
    /// from the exact bounded identifier and locator.
    pub fn try_new(
        id: impl Into<String>,
        locator: impl Into<String>,
    ) -> Result<Self, ContractError> {
        let id = checked_text("hypothesis_source.id", id)?;
        let locator = checked_text("hypothesis_source.locator", locator)?;
        let declaration_hash = hypothesis_source_declaration_hash(&id, &locator);
        Ok(Self {
            id,
            locator,
            declaration_hash,
        })
    }

    fn from_canonical_parts(
        id: String,
        locator: String,
        encoded_hash: ContentHash,
    ) -> Result<Self, ContractError> {
        let source = Self::try_new(id, locator)?;
        if source.declaration_hash != encoded_hash {
            return Err(ContractError::new(
                "EulerContractHypothesisSourceHashMismatch",
                "hypothesis-source hash does not match its canonical identifier and locator",
            ));
        }
        Ok(source)
    }

    #[must_use]
    /// Stable source-declaration identity.
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    /// Human- or machine-resolvable source locator.
    pub fn locator(&self) -> &str {
        &self.locator
    }

    /// Hash of this source declaration, not a custody or authenticity proof.
    #[must_use]
    pub const fn declaration_hash(&self) -> ContentHash {
        self.declaration_hash
    }

    /// Recompute the exact source-declaration identity. This validates only
    /// the local declaration bytes, never source authenticity or custody.
    pub fn verify_identity(&self) -> Result<(), ContractError> {
        let expected = hypothesis_source_declaration_hash(&self.id, &self.locator);
        if self.declaration_hash != expected {
            return Err(ContractError::new(
                "EulerContractHypothesisSourceHashMismatch",
                "hypothesis-source identity does not match its exact declaration",
            ));
        }
        Ok(())
    }
}

fn hypothesis_source_declaration_hash(id: &str, locator: &str) -> ContentHash {
    hypothesis_source_declaration_hash_with_schema(
        EULER_CONTRACT_SCHEMA_VERSION,
        HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
        id,
        locator,
    )
}

fn hypothesis_source_declaration_hash_with_schema(
    version: u32,
    domain: &str,
    id: &str,
    locator: &str,
) -> ContentHash {
    let mut preimage = Vec::with_capacity(id.len() + locator.len() + 20);
    preimage.extend_from_slice(&version.to_le_bytes());
    let id_len = u64::try_from(id.len()).expect("bounded source identifier length fits u64");
    preimage.extend_from_slice(&id_len.to_le_bytes());
    preimage.extend_from_slice(id.as_bytes());
    let locator_len = u64::try_from(locator.len()).expect("bounded source locator length fits u64");
    preimage.extend_from_slice(&locator_len.to_le_bytes());
    preimage.extend_from_slice(locator.as_bytes());
    fs_blake3::hash_domain(domain, &preimage)
}

/// Euler-only Context-of-Use fields absent from the generic V&V artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EulerContextExtension {
    users: Vec<String>,
    apparatus_population: String,
    environment_population: String,
    observation_frame: String,
    decision_alternatives: Vec<String>,
    risks: Vec<ScientificRisk>,
    hypothesis_sources: Vec<HypothesisSource>,
}

impl EulerContextExtension {
    /// Construct and canonicalize the Euler-specific addendum.
    pub fn try_new(
        users: Vec<String>,
        apparatus_population: impl Into<String>,
        environment_population: impl Into<String>,
        observation_frame: impl Into<String>,
        decision_alternatives: Vec<String>,
        mut risks: Vec<ScientificRisk>,
        mut hypothesis_sources: Vec<HypothesisSource>,
    ) -> Result<Self, ContractError> {
        if risks.is_empty() || risks.len() > MAX_EULER_CLAIMS * 2 {
            return Err(ContractError::new(
                "EulerContractCardinality",
                "context risks must be nonempty and bounded",
            ));
        }
        risks.sort();
        let original_risk_count = risks.len();
        risks.dedup_by(|left, right| left.code == right.code);
        if risks.len() != original_risk_count {
            return Err(ContractError::new(
                "EulerContractDuplicate",
                "risk codes must be unique",
            ));
        }
        let users = canonical_text_set("context.users", users, false)?;
        let decision_alternatives = canonical_text_set(
            "context.decision_alternatives",
            decision_alternatives,
            false,
        )?;
        for risk in &risks {
            if !decision_alternatives.contains(&risk.decision_alternative) {
                return Err(ContractError::new(
                    "EulerContractDanglingDecisionAlternative",
                    format!(
                        "risk {} references an unknown decision alternative",
                        risk.code
                    ),
                ));
            }
        }
        if hypothesis_sources.is_empty() || hypothesis_sources.len() > MAX_EULER_CLAIMS * 2 {
            return Err(ContractError::new(
                "EulerContractCardinality",
                "hypothesis sources must be nonempty and bounded",
            ));
        }
        hypothesis_sources.sort();
        let original_source_count = hypothesis_sources.len();
        hypothesis_sources.dedup_by(|left, right| left.id == right.id);
        if hypothesis_sources.len() != original_source_count {
            return Err(ContractError::new(
                "EulerContractDuplicate",
                "hypothesis-source ids must be unique",
            ));
        }
        Ok(Self {
            users,
            apparatus_population: checked_text(
                "context.apparatus_population",
                apparatus_population,
            )?,
            environment_population: checked_text(
                "context.environment_population",
                environment_population,
            )?,
            observation_frame: checked_text("context.observation_frame", observation_frame)?,
            decision_alternatives,
            risks,
            hypothesis_sources,
        })
    }

    #[must_use]
    /// Intended user identities in canonical order.
    pub fn users(&self) -> &[String] {
        &self.users
    }

    #[must_use]
    /// Exact prose scope for specimens, bases, and support apparatus.
    pub fn apparatus_population(&self) -> &str {
        &self.apparatus_population
    }

    #[must_use]
    /// Exact prose scope for the admitted environmental population.
    pub fn environment_population(&self) -> &str {
        &self.environment_population
    }

    #[must_use]
    /// Exact observation frame that must also appear in applicability.
    pub fn observation_frame(&self) -> &str {
        &self.observation_frame
    }

    #[must_use]
    /// Canonical decision alternatives.
    pub fn decision_alternatives(&self) -> &[String] {
        &self.decision_alternatives
    }

    #[must_use]
    /// Canonical risk register.
    pub fn risks(&self) -> &[ScientificRisk] {
        &self.risks
    }

    #[must_use]
    /// Declaration-bound sources restricted to hypothesis generation.
    pub fn hypothesis_sources(&self) -> &[HypothesisSource] {
        &self.hypothesis_sources
    }

    fn encode(&self, writer: &mut CanonicalWriter) -> Result<(), ContractError> {
        writer.strings(&self.users)?;
        writer.string(&self.apparatus_population)?;
        writer.string(&self.environment_population)?;
        writer.string(&self.observation_frame)?;
        writer.strings(&self.decision_alternatives)?;
        writer.usize(self.risks.len())?;
        for risk in &self.risks {
            writer.string(&risk.code)?;
            writer.string(&risk.consequence)?;
            writer.u8(risk.severity);
            writer.usize(risk.affected_claims.len())?;
            for claim in &risk.affected_claims {
                writer.u8(claim.tag());
            }
            writer.string(&risk.decision_alternative)?;
        }
        writer.usize(self.hypothesis_sources.len())?;
        for source in &self.hypothesis_sources {
            writer.string(&source.id)?;
            writer.string(&source.locator)?;
            writer
                .bytes
                .extend_from_slice(source.declaration_hash.as_bytes());
        }
        Ok(())
    }
}

/// One Euler claim backed by generic campaign vocabulary and local requirements.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EulerClaimSpec {
    kind: EulerClaimKind,
    campaign: CampaignClaim,
    requirements: Vec<EvidenceRequirement>,
}

impl EulerClaimSpec {
    /// Construct one claim while retaining the generic campaign claim schema.
    #[allow(clippy::too_many_lines)] // One closed nine-kind policy audit is easier to review whole.
    pub fn try_new(
        kind: EulerClaimKind,
        mut campaign: CampaignClaim,
        mut requirements: Vec<EvidenceRequirement>,
    ) -> Result<Self, ContractError> {
        if campaign.id.as_str() != kind.id() {
            return Err(ContractError::new(
                "EulerContractClaimIdentityMismatch",
                format!(
                    "submitted claim id {} must equal frozen id {}",
                    campaign.id,
                    kind.id()
                ),
            ));
        }
        if campaign.qois.is_empty() || campaign.qois.len() > MAX_EULER_CLAIMS {
            return Err(ContractError::new(
                "EulerContractCardinality",
                format!("claim {} needs a bounded nonempty QoI set", kind.id()),
            ));
        }
        campaign.qois.sort();
        let original_qoi_count = campaign.qois.len();
        campaign.qois.dedup();
        if campaign.qois.len() != original_qoi_count {
            return Err(ContractError::new(
                "EulerContractDuplicate",
                format!("claim {} contains a duplicate QoI", kind.id()),
            ));
        }
        let actual_qois = campaign.qois.iter().map(QoiId::as_str).collect::<Vec<_>>();
        let mut expected_qois = kind.required_qoi_ids().to_vec();
        expected_qois.sort_unstable();
        if actual_qois != expected_qois {
            return Err(ContractError::new(
                "EulerContractClaimQoiPolicyMismatch",
                format!(
                    "claim {} must reference exactly {:?}, received {:?}",
                    kind.id(),
                    expected_qois,
                    actual_qois
                ),
            ));
        }
        campaign.hypothesis = checked_text("claim.hypothesis", campaign.hypothesis)?;
        campaign.decision_consequence =
            checked_text("claim.decision_consequence", campaign.decision_consequence)?;
        if campaign.evidence_gaps.len() > MAX_EULER_CLAIMS {
            return Err(ContractError::new(
                "EulerContractCardinality",
                format!("claim {} has too many evidence gaps", kind.id()),
            ));
        }
        campaign
            .evidence_gaps
            .sort_by(|left, right| left.id.cmp(&right.id));
        let mut gap_ids = BTreeSet::new();
        for gap in &mut campaign.evidence_gaps {
            if !campaign.qois.contains(&gap.qoi) {
                return Err(ContractError::new(
                    "EulerContractDanglingQoi",
                    format!("gap {} is outside claim {}", gap.id, kind.id()),
                ));
            }
            if !gap_ids.insert(gap.id.clone()) {
                return Err(ContractError::new(
                    "EulerContractDuplicate",
                    format!("claim {} contains a duplicate evidence gap", kind.id()),
                ));
            }
            gap.expected_evidence = checked_text(
                "claim.evidence_gap.expected",
                core::mem::take(&mut gap.expected_evidence),
            )?;
            gap.description = checked_text(
                "claim.evidence_gap.description",
                core::mem::take(&mut gap.description),
            )?;
        }
        if requirements.is_empty() || requirements.len() > EULER_EVIDENCE_REQUIREMENT_REGISTRY.len()
        {
            return Err(ContractError::new(
                "EulerContractCardinality",
                format!(
                    "claim {} needs a bounded nonempty requirement set",
                    kind.id()
                ),
            ));
        }
        requirements.sort_by_key(|requirement| requirement.tag());
        let original_requirement_count = requirements.len();
        requirements.dedup();
        if requirements.len() != original_requirement_count {
            return Err(ContractError::new(
                "EulerContractDuplicate",
                format!(
                    "claim {} contains a duplicate evidence requirement",
                    kind.id()
                ),
            ));
        }
        let mut expected_requirements = kind.required_evidence().to_vec();
        expected_requirements.sort_by_key(|requirement| requirement.tag());
        if requirements != expected_requirements {
            return Err(ContractError::new(
                "EulerContractClaimEvidencePolicyMismatch",
                format!(
                    "claim {} requirements do not match claim-policy schema v{}",
                    kind.id(),
                    EULER_CLAIM_POLICY_SCHEMA_VERSION
                ),
            ));
        }
        Ok(Self {
            kind,
            campaign,
            requirements,
        })
    }

    #[must_use]
    /// Closed local claim kind.
    pub const fn kind(&self) -> EulerClaimKind {
        self.kind
    }

    #[must_use]
    /// Generic fs-ir campaign vocabulary retained by this claim.
    pub const fn campaign(&self) -> &CampaignClaim {
        &self.campaign
    }

    #[must_use]
    /// Exact minimum evidence roles in canonical order.
    pub fn requirements(&self) -> &[EvidenceRequirement] {
        &self.requirements
    }

    #[must_use]
    /// Version of the separately frozen claim policy.
    pub const fn policy_schema_version(&self) -> u32 {
        EULER_CLAIM_POLICY_SCHEMA_VERSION
    }

    #[must_use]
    /// Decision semantics interpreted by later exact scoring artifacts.
    pub const fn acceptance_family(&self) -> EulerAcceptanceFamily {
        self.kind.acceptance_family()
    }

    fn encode(&self, writer: &mut CanonicalWriter) -> Result<(), ContractError> {
        writer.u8(self.kind.tag());
        writer.u32(EULER_CLAIM_POLICY_SCHEMA_VERSION);
        writer.u8(self.kind.acceptance_family().tag());
        writer.string(self.campaign.id.as_str())?;
        writer.usize(self.campaign.qois.len())?;
        for qoi in &self.campaign.qois {
            writer.string(qoi.as_str())?;
        }
        writer.string(&self.campaign.hypothesis)?;
        writer.string(&self.campaign.decision_consequence)?;
        writer.usize(self.campaign.evidence_gaps.len())?;
        for gap in &self.campaign.evidence_gaps {
            writer.string(gap.id.as_str())?;
            writer.string(gap.qoi.as_str())?;
            writer.string(&gap.expected_evidence)?;
            writer.string(&gap.description)?;
        }
        writer.usize(self.requirements.len())?;
        for requirement in &self.requirements {
            writer.u8(requirement.tag());
        }
        Ok(())
    }
}

/// Canonical, acyclic Euler claim dependency graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EulerClaimGraph {
    claims: BTreeMap<EulerClaimKind, EulerClaimSpec>,
    dependencies: Vec<ClaimDependency>,
}

impl EulerClaimGraph {
    /// Canonicalize claims/dependencies and reject dangling or cyclic edges.
    pub fn try_new(
        claims: Vec<EulerClaimSpec>,
        mut dependencies: Vec<ClaimDependency>,
    ) -> Result<Self, ContractError> {
        if claims.is_empty() || claims.len() > MAX_EULER_CLAIMS {
            return Err(ContractError::new(
                "EulerContractCardinality",
                "claim graph must contain between one and nine claims",
            ));
        }
        if dependencies.len() > MAX_EULER_CLAIMS * (MAX_EULER_CLAIMS - 1) {
            return Err(ContractError::new(
                "EulerContractCardinality",
                "claim graph contains too many dependency rows",
            ));
        }
        let mut by_kind = BTreeMap::new();
        let mut by_id = BTreeMap::new();
        let mut global_gap_ids = BTreeSet::new();
        for claim in claims {
            for gap in &claim.campaign.evidence_gaps {
                if !global_gap_ids.insert(gap.id.clone()) {
                    return Err(ContractError::new(
                        "EulerContractDuplicate",
                        format!("evidence-gap id {} is not globally unique", gap.id),
                    ));
                }
            }
            if by_id
                .insert(claim.campaign.id.clone(), claim.kind)
                .is_some()
                || by_kind.insert(claim.kind, claim).is_some()
            {
                return Err(ContractError::new(
                    "EulerContractDuplicate",
                    "claim graph contains a duplicate kind or id",
                ));
            }
        }
        dependencies.sort_by(|left, right| {
            left.prerequisite
                .cmp(&right.prerequisite)
                .then_with(|| left.dependent.cmp(&right.dependent))
                .then_with(|| {
                    evidence_use_tag(left.use_kind).cmp(&evidence_use_tag(right.use_kind))
                })
        });
        let original_dependency_count = dependencies.len();
        dependencies.dedup();
        if dependencies.len() != original_dependency_count {
            return Err(ContractError::new(
                "EulerContractDuplicate",
                "claim graph contains a duplicate dependency",
            ));
        }
        let mut endpoint_pairs = BTreeSet::new();
        for dependency in &dependencies {
            if !by_id.contains_key(&dependency.prerequisite)
                || !by_id.contains_key(&dependency.dependent)
            {
                return Err(ContractError::new(
                    "EulerContractDanglingClaim",
                    format!(
                        "dependency {} -> {} names a missing claim",
                        dependency.prerequisite, dependency.dependent
                    ),
                ));
            }
            if dependency.prerequisite == dependency.dependent {
                return Err(ContractError::new(
                    "EulerContractClaimCycle",
                    "a claim cannot depend on itself",
                ));
            }
            if !endpoint_pairs.insert((
                dependency.prerequisite.clone(),
                dependency.dependent.clone(),
            )) {
                return Err(ContractError::new(
                    "EulerContractDependencyRoleCollision",
                    format!(
                        "dependency {} -> {} cannot be typed as both calibration and validation evidence",
                        dependency.prerequisite, dependency.dependent
                    ),
                ));
            }
        }
        if !is_acyclic(&by_id, &dependencies) {
            return Err(ContractError::new(
                "EulerContractClaimCycle",
                "claim dependencies contain a directed cycle",
            ));
        }
        Ok(Self {
            claims: by_kind,
            dependencies,
        })
    }

    #[must_use]
    /// Claims keyed by their closed local kind.
    pub const fn claims(&self) -> &BTreeMap<EulerClaimKind, EulerClaimSpec> {
        &self.claims
    }

    #[must_use]
    /// Look up one closed claim kind.
    pub fn claim(&self, kind: EulerClaimKind) -> Option<&EulerClaimSpec> {
        self.claims.get(&kind)
    }

    #[must_use]
    /// Canonically ordered directed claim dependencies.
    pub fn dependencies(&self) -> &[ClaimDependency] {
        &self.dependencies
    }

    /// Versioned canonical bytes for the Euler-local graph.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ContractError> {
        let mut writer = CanonicalWriter::default();
        writer.bytes.extend_from_slice(GRAPH_MAGIC);
        writer.u32(EULER_CLAIM_POLICY_SCHEMA_VERSION);
        writer.usize(self.claims.len())?;
        let mut claims = self.claims.values().collect::<Vec<_>>();
        claims.sort_by_key(|claim| claim.kind.tag());
        for claim in claims {
            claim.encode(&mut writer)?;
        }
        writer.usize(self.dependencies.len())?;
        for dependency in &self.dependencies {
            writer.string(dependency.prerequisite.as_str())?;
            writer.string(dependency.dependent.as_str())?;
            writer.u8(evidence_use_tag(dependency.use_kind));
        }
        if writer.bytes.len() > MAX_EULER_GRAPH_BYTES {
            return Err(ContractError::new(
                "EulerContractGraphTooLarge",
                "canonical claim graph exceeds its byte budget",
            ));
        }
        Ok(writer.bytes)
    }

    /// Strictly decode and re-canonicalize the current graph schema.
    #[allow(clippy::too_many_lines)] // Canonical field order is a single auditable state machine.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, ContractError> {
        if bytes.len() > MAX_EULER_GRAPH_BYTES {
            return Err(ContractError::new(
                "EulerContractGraphTooLarge",
                "claim graph transport exceeds its byte budget",
            ));
        }
        let mut reader =
            CanonicalReader::new(bytes, "EulerContractMalformedCanonical", "claim graph");
        if reader.take(GRAPH_MAGIC.len())? != GRAPH_MAGIC {
            return Err(ContractError::new(
                "EulerContractMalformedCanonical",
                "claim graph magic is invalid",
            ));
        }
        let version = reader.u32()?;
        if version != EULER_CLAIM_POLICY_SCHEMA_VERSION {
            return Err(ContractError::new(
                "EulerContractUnsupportedVersion",
                format!(
                    "claim graph schema {version} is unsupported; v1 has no predecessor migration"
                ),
            ));
        }
        let claim_count = reader.usize()?;
        if claim_count == 0 || claim_count > MAX_EULER_CLAIMS {
            return Err(ContractError::new(
                "EulerContractCardinality",
                "claim graph transport has an invalid claim count",
            ));
        }
        let mut claims = Vec::with_capacity(claim_count);
        for _ in 0..claim_count {
            let kind = EulerClaimKind::from_tag(reader.u8()?)?;
            let policy_version = reader.u32()?;
            if policy_version != EULER_CLAIM_POLICY_SCHEMA_VERSION {
                return Err(ContractError::new(
                    "EulerContractUnsupportedClaimPolicyVersion",
                    format!("claim policy schema {policy_version} is unsupported"),
                ));
            }
            let acceptance_family = reader.u8()?;
            if acceptance_family != kind.acceptance_family().tag() {
                return Err(ContractError::new(
                    "EulerContractClaimAcceptanceMismatch",
                    format!(
                        "claim {} carries acceptance-family tag {acceptance_family}, expected {}",
                        kind.id(),
                        kind.acceptance_family().tag()
                    ),
                ));
            }
            let id = CampaignClaimId::try_new(reader.string()?).map_err(|error| {
                ContractError::new("EulerContractMalformedCanonical", error.to_string())
            })?;
            let qoi_count = reader.usize()?;
            if qoi_count == 0 || qoi_count > MAX_EULER_CLAIMS {
                return Err(ContractError::new(
                    "EulerContractCardinality",
                    "claim transport has an invalid QoI count",
                ));
            }
            let mut qois = Vec::with_capacity(qoi_count);
            for _ in 0..qoi_count {
                qois.push(QoiId::try_new(reader.string()?).map_err(|error| {
                    ContractError::new("EulerContractMalformedCanonical", error.to_string())
                })?);
            }
            let hypothesis = reader.string()?;
            let decision_consequence = reader.string()?;
            let gap_count = reader.usize()?;
            if gap_count > MAX_EULER_CLAIMS {
                return Err(ContractError::new(
                    "EulerContractCardinality",
                    "claim transport has too many evidence gaps",
                ));
            }
            let mut evidence_gaps = Vec::with_capacity(gap_count);
            for _ in 0..gap_count {
                let gap_id = EvidenceGapId::try_new(reader.string()?).map_err(|error| {
                    ContractError::new("EulerContractMalformedCanonical", error.to_string())
                })?;
                let qoi = QoiId::try_new(reader.string()?).map_err(|error| {
                    ContractError::new("EulerContractMalformedCanonical", error.to_string())
                })?;
                evidence_gaps.push(EvidenceGap {
                    id: gap_id,
                    qoi,
                    expected_evidence: reader.string()?,
                    description: reader.string()?,
                });
            }
            let requirement_count = reader.usize()?;
            if requirement_count == 0
                || requirement_count > EULER_EVIDENCE_REQUIREMENT_REGISTRY.len()
            {
                return Err(ContractError::new(
                    "EulerContractCardinality",
                    "claim transport has an invalid requirement count",
                ));
            }
            let mut requirements = Vec::with_capacity(requirement_count);
            for _ in 0..requirement_count {
                requirements.push(EvidenceRequirement::from_tag(reader.u8()?)?);
            }
            claims.push(EulerClaimSpec::try_new(
                kind,
                CampaignClaim {
                    id,
                    qois,
                    hypothesis,
                    decision_consequence,
                    evidence_gaps,
                },
                requirements,
            )?);
        }
        let dependency_count = reader.usize()?;
        if dependency_count > MAX_EULER_CLAIMS * (MAX_EULER_CLAIMS - 1) {
            return Err(ContractError::new(
                "EulerContractCardinality",
                "claim transport has too many dependencies",
            ));
        }
        let mut dependencies = Vec::with_capacity(dependency_count);
        for _ in 0..dependency_count {
            let prerequisite = CampaignClaimId::try_new(reader.string()?).map_err(|error| {
                ContractError::new("EulerContractMalformedCanonical", error.to_string())
            })?;
            let dependent = CampaignClaimId::try_new(reader.string()?).map_err(|error| {
                ContractError::new("EulerContractMalformedCanonical", error.to_string())
            })?;
            let use_kind = match reader.u8()? {
                1 => EvidenceUse::CalibrationInput,
                2 => EvidenceUse::ValidationInput,
                tag => {
                    return Err(ContractError::new(
                        "EulerContractMalformedCanonical",
                        format!("unknown evidence-use tag {tag}"),
                    ));
                }
            };
            dependencies.push(ClaimDependency {
                prerequisite,
                dependent,
                use_kind,
            });
        }
        if !reader.is_finished() {
            return Err(ContractError::new(
                "EulerContractMalformedCanonical",
                "claim graph transport has trailing bytes",
            ));
        }
        let decoded = Self::try_new(claims, dependencies)?;
        if decoded.canonical_bytes()?.as_slice() != bytes {
            return Err(ContractError::new(
                "EulerContractNonCanonical",
                "claim graph bytes are not a canonical fixed point",
            ));
        }
        Ok(decoded)
    }

    /// Content hash of the exact graph transport.
    pub fn content_hash(&self) -> Result<ContentHash, ContractError> {
        Ok(fs_blake3::hash_domain(
            EULER_CLAIM_GRAPH_IDENTITY_DOMAIN,
            &self.canonical_bytes()?,
        ))
    }
}

const fn evidence_use_tag(use_kind: EvidenceUse) -> u8 {
    match use_kind {
        EvidenceUse::CalibrationInput => 1,
        EvidenceUse::ValidationInput => 2,
    }
}

fn is_acyclic(
    by_id: &BTreeMap<CampaignClaimId, EulerClaimKind>,
    dependencies: &[ClaimDependency],
) -> bool {
    let mut indegree = by_id
        .keys()
        .cloned()
        .map(|id| (id, 0_usize))
        .collect::<BTreeMap<_, _>>();
    let mut outgoing = BTreeMap::<CampaignClaimId, Vec<CampaignClaimId>>::new();
    for dependency in dependencies {
        let Some(degree) = indegree.get_mut(&dependency.dependent) else {
            return false;
        };
        *degree += 1;
        outgoing
            .entry(dependency.prerequisite.clone())
            .or_default()
            .push(dependency.dependent.clone());
    }
    let mut ready = indegree
        .iter()
        .filter_map(|(id, degree)| (*degree == 0).then_some(id.clone()))
        .collect::<VecDeque<_>>();
    let mut visited = 0_usize;
    while let Some(id) = ready.pop_front() {
        visited += 1;
        if let Some(next) = outgoing.get(&id) {
            for dependent in next {
                let Some(degree) = indegree.get_mut(dependent) else {
                    return false;
                };
                *degree -= 1;
                if *degree == 0 {
                    ready.push_back(dependent.clone());
                }
            }
        }
    }
    visited == by_id.len()
}

/// Owner/source-schema role in the exact artifact matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum OwnerRole {
    /// Generic Context-of-Use artifact.
    ContextOfUse,
    /// Generic typed V&V artifacts referenced as evidence.
    VvEvidenceArtifact,
    /// Generic whole-case structural schema-admission receipt.
    VvSchemaAdmissionReceipt,
    /// Generic claim-id/dependency vocabulary.
    CampaignClaimVocabulary,
    /// Generic canonical no-claim set.
    NoClaimBoundary,
    /// Euler-local hypothesis-only source declaration.
    HypothesisSourceDeclaration,
    /// Euler-local nine-claim graph and minimum-evidence policy.
    EulerClaimGraph,
    /// Euler-only domain extension, graph policy, and composite identity.
    EulerScientificContract,
    /// Euler-local evidence-reference packet.
    ClaimEvidencePacket,
    /// Euler-local direct prerequisite assessment receipt.
    PrerequisiteAssessmentReceipt,
    /// Euler-local structural policy assessment.
    ClaimPolicyAssessment,
    /// Euler-local exact-contract check receipt.
    ContractCheckReceipt,
    /// Euler-local bounded log for one claim-policy assessment.
    ClaimPolicyAssessmentLog,
    /// This separately versioned owner/role routing registry itself.
    OwnerMatrixRegistry,
    /// Downstream Context-bound detailed-observable-to-aggregate-QoI receipt.
    AggregateQoiDerivationReceipt,
}

/// Single executable registry for the complete v1 artifact-owner matrix.
pub(crate) const EULER_OWNER_ROLE_REGISTRY: [OwnerRole; 15] = [
    OwnerRole::ContextOfUse,
    OwnerRole::VvEvidenceArtifact,
    OwnerRole::VvSchemaAdmissionReceipt,
    OwnerRole::CampaignClaimVocabulary,
    OwnerRole::NoClaimBoundary,
    OwnerRole::HypothesisSourceDeclaration,
    OwnerRole::EulerClaimGraph,
    OwnerRole::EulerScientificContract,
    OwnerRole::ClaimEvidencePacket,
    OwnerRole::PrerequisiteAssessmentReceipt,
    OwnerRole::ClaimPolicyAssessment,
    OwnerRole::ContractCheckReceipt,
    OwnerRole::ClaimPolicyAssessmentLog,
    OwnerRole::OwnerMatrixRegistry,
    OwnerRole::AggregateQoiDerivationReceipt,
];

impl OwnerRole {
    /// Complete closed artifact/owner matrix for v1.
    pub const ALL: [Self; 15] = EULER_OWNER_ROLE_REGISTRY;

    #[must_use]
    const fn tag(self) -> u8 {
        match self {
            Self::ContextOfUse => 1,
            Self::VvEvidenceArtifact => 2,
            Self::VvSchemaAdmissionReceipt => 3,
            Self::CampaignClaimVocabulary => 4,
            Self::NoClaimBoundary => 5,
            Self::HypothesisSourceDeclaration => 6,
            Self::EulerClaimGraph => 7,
            Self::EulerScientificContract => 8,
            Self::ClaimEvidencePacket => 9,
            Self::PrerequisiteAssessmentReceipt => 10,
            Self::ClaimPolicyAssessment => 11,
            Self::ContractCheckReceipt => 12,
            Self::ClaimPolicyAssessmentLog => 13,
            // V1 role tags are append-only: adding the downstream receipt must
            // not silently renumber the pre-existing routing-registry role.
            Self::OwnerMatrixRegistry => 14,
            Self::AggregateQoiDerivationReceipt => 15,
        }
    }

    #[must_use]
    /// Required owning crate for this generic or Euler-local role.
    pub const fn expected_owner_crate(self) -> &'static str {
        match self {
            Self::ContextOfUse | Self::VvEvidenceArtifact | Self::VvSchemaAdmissionReceipt => {
                "fs-evidence"
            }
            Self::CampaignClaimVocabulary => "fs-ir",
            Self::NoClaimBoundary => "fs-govern",
            Self::HypothesisSourceDeclaration
            | Self::EulerClaimGraph
            | Self::EulerScientificContract
            | Self::ClaimEvidencePacket
            | Self::PrerequisiteAssessmentReceipt
            | Self::ClaimPolicyAssessment
            | Self::ContractCheckReceipt
            | Self::ClaimPolicyAssessmentLog
            | Self::AggregateQoiDerivationReceipt
            | Self::OwnerMatrixRegistry => "fs-euler-disc-e2e",
        }
    }

    #[must_use]
    /// Required source-schema identity for this role.
    pub const fn expected_source_schema(self) -> &'static str {
        match self {
            Self::ContextOfUse | Self::VvEvidenceArtifact => VV_ARTIFACT_FAMILY,
            Self::VvSchemaAdmissionReceipt => VV_SCHEMA_ADMISSION_RECEIPT_IDENTITY_DOMAIN,
            Self::CampaignClaimVocabulary => FS_IR_CAMPAIGN_SOURCE_SCHEMA,
            Self::NoClaimBoundary => FS_GOVERN_AUTHORITY_SOURCE_SCHEMA,
            Self::HypothesisSourceDeclaration => HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
            Self::EulerClaimGraph => EULER_CLAIM_GRAPH_IDENTITY_DOMAIN,
            Self::EulerScientificContract => EULER_CONTRACT_IDENTITY_DOMAIN,
            Self::ClaimEvidencePacket => EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN,
            Self::PrerequisiteAssessmentReceipt => EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,
            Self::ClaimPolicyAssessment => EULER_ASSESSMENT_IDENTITY_DOMAIN,
            Self::ContractCheckReceipt => CONTRACT_CHECK_RECEIPT_DOMAIN,
            Self::ClaimPolicyAssessmentLog => CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN,
            Self::AggregateQoiDerivationReceipt => EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA,
            Self::OwnerMatrixRegistry => EULER_OWNER_MATRIX_IDENTITY_DOMAIN,
        }
    }

    #[must_use]
    /// Maximum authority the role can contribute locally.
    pub const fn expected_authority_ceiling(self) -> AuthorityCeiling {
        match self {
            Self::ContextOfUse => AuthorityCeiling::StructuralContextDeclaration,
            Self::VvEvidenceArtifact => AuthorityCeiling::StructuralEvidenceReferenceOnly,
            Self::VvSchemaAdmissionReceipt => AuthorityCeiling::StructuralSchemaAdmissionOnly,
            Self::CampaignClaimVocabulary => AuthorityCeiling::CampaignVocabularyOnly,
            Self::NoClaimBoundary => AuthorityCeiling::CanonicalNoClaimBoundary,
            Self::HypothesisSourceDeclaration => AuthorityCeiling::HypothesisOnly,
            Self::EulerClaimGraph => AuthorityCeiling::StructuralClaimPolicyOnly,
            Self::EulerScientificContract
            | Self::ClaimEvidencePacket
            | Self::ClaimPolicyAssessment => AuthorityCeiling::CandidateEligibilityOnly,
            Self::PrerequisiteAssessmentReceipt => AuthorityCeiling::StructuralDependencyOnly,
            Self::ContractCheckReceipt => AuthorityCeiling::StructuralCheckOnly,
            Self::ClaimPolicyAssessmentLog => AuthorityCeiling::DiagnosticRetentionOnly,
            Self::AggregateQoiDerivationReceipt => {
                AuthorityCeiling::StructuralEvidenceReferenceOnly
            }
            Self::OwnerMatrixRegistry => AuthorityCeiling::StructuralRoutingRegistryOnly,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, ContractError> {
        EULER_OWNER_ROLE_REGISTRY
            .into_iter()
            .find(|role| role.tag() == tag)
            .ok_or_else(|| {
                ContractError::new(
                    "EulerOwnerMatrixMalformedCanonical",
                    format!("unknown owner-role tag {tag}"),
                )
            })
    }
}

/// Closed non-widening authority ceiling for one schema owner.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AuthorityCeiling {
    /// Declares a structural Context of Use but no adequacy result.
    StructuralContextDeclaration,
    /// Carries typed generic evidence references without authenticating them.
    StructuralEvidenceReferenceOnly,
    /// Establishes generic schema well-formedness but no scientific adequacy.
    StructuralSchemaAdmissionOnly,
    /// Supplies claim/dependency vocabulary only.
    CampaignVocabularyOnly,
    /// Supplies a canonical set of prohibitions only.
    CanonicalNoClaimBoundary,
    /// Generates hypotheses but cannot satisfy evidence requirements.
    HypothesisOnly,
    /// Declares claim topology and minima without satisfying those minima.
    StructuralClaimPolicyOnly,
    /// Supports local candidate eligibility but no govern grant.
    CandidateEligibilityOnly,
    /// Binds a direct prerequisite relation without transferring physical authority.
    StructuralDependencyOnly,
    /// Establishes exact structural equality only.
    StructuralCheckOnly,
    /// Retains diagnostics without granting authority.
    DiagnosticRetentionOnly,
    /// Routes artifact roles without granting any addressed schema's authority.
    StructuralRoutingRegistryOnly,
}

impl AuthorityCeiling {
    #[must_use]
    /// Stable machine code for this authority ceiling.
    pub const fn code(self) -> &'static str {
        match self {
            Self::StructuralContextDeclaration => "structural-context-declaration",
            Self::StructuralEvidenceReferenceOnly => "structural-evidence-reference-only",
            Self::StructuralSchemaAdmissionOnly => "structural-schema-admission-only",
            Self::CampaignVocabularyOnly => "campaign-vocabulary-only",
            Self::CanonicalNoClaimBoundary => "canonical-no-claim-boundary",
            Self::HypothesisOnly => "hypothesis-only",
            Self::StructuralClaimPolicyOnly => "structural-claim-policy-only",
            Self::CandidateEligibilityOnly => "candidate-eligibility-only",
            Self::StructuralDependencyOnly => "structural-dependency-only",
            Self::StructuralCheckOnly => "structural-check-only",
            Self::DiagnosticRetentionOnly => "diagnostic-retention-only",
            Self::StructuralRoutingRegistryOnly => "structural-routing-registry-only",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::StructuralContextDeclaration => 1,
            Self::StructuralEvidenceReferenceOnly => 2,
            Self::StructuralSchemaAdmissionOnly => 3,
            Self::CampaignVocabularyOnly => 4,
            Self::CanonicalNoClaimBoundary => 5,
            Self::HypothesisOnly => 6,
            Self::StructuralClaimPolicyOnly => 7,
            Self::CandidateEligibilityOnly => 8,
            Self::StructuralDependencyOnly => 9,
            Self::StructuralCheckOnly => 10,
            Self::DiagnosticRetentionOnly => 11,
            Self::StructuralRoutingRegistryOnly => 12,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, ContractError> {
        match tag {
            1 => Ok(Self::StructuralContextDeclaration),
            2 => Ok(Self::StructuralEvidenceReferenceOnly),
            3 => Ok(Self::StructuralSchemaAdmissionOnly),
            4 => Ok(Self::CampaignVocabularyOnly),
            5 => Ok(Self::CanonicalNoClaimBoundary),
            6 => Ok(Self::HypothesisOnly),
            7 => Ok(Self::StructuralClaimPolicyOnly),
            8 => Ok(Self::CandidateEligibilityOnly),
            9 => Ok(Self::StructuralDependencyOnly),
            10 => Ok(Self::StructuralCheckOnly),
            11 => Ok(Self::DiagnosticRetentionOnly),
            12 => Ok(Self::StructuralRoutingRegistryOnly),
            _ => Err(ContractError::new(
                "EulerOwnerMatrixMalformedCanonical",
                format!("unknown authority-ceiling tag {tag}"),
            )),
        }
    }
}

/// One exact artifact/source-schema owner row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnerRow {
    role: OwnerRole,
    owner_crate: String,
    source_schema: String,
    authority_ceiling: AuthorityCeiling,
}

impl OwnerRow {
    /// Construct one row and refuse every owner/schema/ceiling substitution.
    pub fn try_new(
        role: OwnerRole,
        owner_crate: impl Into<String>,
        source_schema: impl Into<String>,
        authority_ceiling: AuthorityCeiling,
    ) -> Result<Self, ContractError> {
        let row = Self {
            role,
            owner_crate: checked_text("owner.owner_crate", owner_crate)?,
            source_schema: checked_text("owner.source_schema", source_schema)?,
            authority_ceiling,
        };
        if row.owner_crate != role.expected_owner_crate()
            || row.source_schema != role.expected_source_schema()
            || row.authority_ceiling != role.expected_authority_ceiling()
        {
            return Err(ContractError::new(
                "EulerContractGenericSchemaFork",
                format!(
                    "owner row {role:?} must use owner {}, schema {}, and ceiling {}",
                    role.expected_owner_crate(),
                    role.expected_source_schema(),
                    role.expected_authority_ceiling().code()
                ),
            ));
        }
        Ok(row)
    }

    #[must_use]
    /// Role described by this row.
    pub const fn role(&self) -> OwnerRole {
        self.role
    }

    #[must_use]
    /// Exact owning crate.
    pub fn owner_crate(&self) -> &str {
        &self.owner_crate
    }

    #[must_use]
    /// Exact source-schema identity.
    pub fn source_schema(&self) -> &str {
        &self.source_schema
    }

    #[must_use]
    /// Non-widening authority ceiling.
    pub const fn authority_ceiling(&self) -> AuthorityCeiling {
        self.authority_ceiling
    }
}

/// Domain-separated identity of one exact owner/role routing registry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OwnerMatrixIdentity(ContentHash);

impl OwnerMatrixIdentity {
    /// Raw domain-separated registry hash.
    #[must_use]
    pub const fn as_hash(self) -> ContentHash {
        self.0
    }
}

impl fmt::Display for OwnerMatrixIdentity {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

/// Exact owner matrix preventing an Euler-local duplicate of a generic schema.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnerMatrix {
    schema_version: u32,
    rows: BTreeMap<OwnerRole, OwnerRow>,
    identity: OwnerMatrixIdentity,
}

impl OwnerMatrix {
    /// Construct the matrix and enforce the generic owner chokepoints.
    pub fn try_new(rows: Vec<OwnerRow>) -> Result<Self, ContractError> {
        if rows.len() != EULER_OWNER_ROLE_REGISTRY.len() {
            return Err(ContractError::new(
                "EulerContractOwnerMatrixIncomplete",
                format!(
                    "owner matrix must contain exactly {} roles",
                    EULER_OWNER_ROLE_REGISTRY.len()
                ),
            ));
        }
        let mut by_role = BTreeMap::new();
        for row in rows {
            let role = row.role;
            if by_role.insert(role, row).is_some() {
                return Err(ContractError::new(
                    "EulerContractDuplicate",
                    "owner matrix contains a duplicate role",
                ));
            }
        }
        for role in EULER_OWNER_ROLE_REGISTRY {
            let Some(row) = by_role.get(&role) else {
                return Err(ContractError::new(
                    "EulerContractOwnerMatrixIncomplete",
                    format!("owner matrix is missing role {role:?}"),
                ));
            };
            if row.owner_crate != role.expected_owner_crate()
                || row.source_schema != role.expected_source_schema()
                || row.authority_ceiling != role.expected_authority_ceiling()
            {
                return Err(ContractError::new(
                    "EulerContractGenericSchemaFork",
                    format!("{role:?} does not match its frozen owner/schema/ceiling row"),
                ));
            }
        }
        let canonical =
            encode_owner_matrix_components(EULER_OWNER_MATRIX_SCHEMA_VERSION, &by_role)?;
        let identity = owner_matrix_identity(&canonical);
        Ok(Self {
            schema_version: EULER_OWNER_MATRIX_SCHEMA_VERSION,
            rows: by_role,
            identity,
        })
    }

    /// Explicit migration policy for the independently versioned routing
    /// registry. V1 has no predecessor and never approximately reinterprets
    /// another version.
    pub fn migration_policy(schema_version: u32) -> Result<(), ContractError> {
        if schema_version == EULER_OWNER_MATRIX_SCHEMA_VERSION {
            Ok(())
        } else {
            Err(ContractError::new(
                "EulerOwnerMatrixUnsupportedVersion",
                format!(
                    "owner-matrix schema {schema_version} is unsupported; v1 has no predecessor migration"
                ),
            ))
        }
    }

    #[must_use]
    /// Independently versioned routing-registry schema.
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    /// Owner rows keyed by their exact role.
    pub const fn rows(&self) -> &BTreeMap<OwnerRole, OwnerRow> {
        &self.rows
    }

    #[must_use]
    /// Domain-separated identity of the exact canonical routing registry.
    pub const fn identity(&self) -> OwnerMatrixIdentity {
        self.identity
    }

    /// Canonical bounded owner-matrix transport.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ContractError> {
        encode_owner_matrix_components(self.schema_version, &self.rows)
    }

    /// Recompute the exact registry identity from its canonical fields.
    pub fn verify_identity(&self) -> Result<(), ContractError> {
        let expected = owner_matrix_identity(&self.canonical_bytes()?);
        if self.schema_version != EULER_OWNER_MATRIX_SCHEMA_VERSION || self.identity != expected {
            return Err(ContractError::new(
                "EulerOwnerMatrixIdentityMismatch",
                "owner-matrix identity or schema version is stale",
            ));
        }
        Ok(())
    }

    /// Strictly decode and re-canonicalize the current routing-registry schema.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, ContractError> {
        if bytes.len() > MAX_OWNER_MATRIX_BYTES {
            return Err(ContractError::new(
                "EulerOwnerMatrixTooLarge",
                "owner-matrix transport exceeds its byte budget",
            ));
        }
        let mut reader =
            CanonicalReader::new(bytes, "EulerOwnerMatrixMalformedCanonical", "owner matrix");
        if reader.take(OWNER_MATRIX_MAGIC.len())? != OWNER_MATRIX_MAGIC {
            return Err(ContractError::new(
                "EulerOwnerMatrixMalformedCanonical",
                "owner-matrix magic is invalid",
            ));
        }
        let schema_version = reader.u32()?;
        Self::migration_policy(schema_version)?;
        let owner_count =
            reader.bounded_len("owner_matrix.rows", EULER_OWNER_ROLE_REGISTRY.len(), false)?;
        if owner_count != EULER_OWNER_ROLE_REGISTRY.len() {
            return Err(ContractError::new(
                "EulerContractOwnerMatrixIncomplete",
                format!(
                    "canonical owner matrix must contain exactly {} rows",
                    EULER_OWNER_ROLE_REGISTRY.len()
                ),
            ));
        }
        let mut rows = Vec::with_capacity(owner_count);
        for _ in 0..owner_count {
            rows.push(OwnerRow::try_new(
                OwnerRole::from_tag(reader.u8()?)?,
                reader.string()?,
                reader.string()?,
                AuthorityCeiling::from_tag(reader.u8()?)?,
            )?);
        }
        if !reader.is_finished() {
            return Err(ContractError::new(
                "EulerOwnerMatrixMalformedCanonical",
                "owner-matrix transport has trailing bytes",
            ));
        }
        let decoded = Self::try_new(rows)?;
        if decoded.canonical_bytes()?.as_slice() != bytes {
            return Err(ContractError::new(
                "EulerOwnerMatrixNonCanonical",
                "owner-matrix bytes are not a canonical fixed point",
            ));
        }
        Ok(decoded)
    }
}

fn owner_matrix_identity(bytes: &[u8]) -> OwnerMatrixIdentity {
    OwnerMatrixIdentity(fs_blake3::hash_domain(
        EULER_OWNER_MATRIX_IDENTITY_DOMAIN,
        bytes,
    ))
}

fn encode_owner_matrix_components(
    schema_version: u32,
    rows: &BTreeMap<OwnerRole, OwnerRow>,
) -> Result<Vec<u8>, ContractError> {
    let mut writer = CanonicalWriter::default();
    writer.bytes.extend_from_slice(OWNER_MATRIX_MAGIC);
    writer.u32(schema_version);
    writer.usize(rows.len())?;
    let mut rows = rows.values().collect::<Vec<_>>();
    rows.sort_by_key(|row| row.role.tag());
    for row in rows {
        writer.u8(row.role.tag());
        writer.string(&row.owner_crate)?;
        writer.string(&row.source_schema)?;
        writer.u8(row.authority_ceiling.tag());
    }
    if writer.bytes.len() > MAX_OWNER_MATRIX_BYTES {
        return Err(ContractError::new(
            "EulerOwnerMatrixTooLarge",
            "canonical owner matrix exceeds its byte budget",
        ));
    }
    Ok(writer.bytes)
}

/// Complete immutable Euler scientific-contract addendum.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerScientificContract {
    schema_version: u32,
    context: ContextOfUse,
    context_bytes: Vec<u8>,
    context_hash: ContentHash,
    extension: EulerContextExtension,
    claim_graph: EulerClaimGraph,
    no_claims: NoClaimBoundary,
    owner_matrix: OwnerMatrix,
    identity: ContractIdentity,
}

fn scientific_contract_identity(bytes: &[u8]) -> ContractIdentity {
    ContractIdentity(fs_blake3::hash_domain(
        EULER_CONTRACT_IDENTITY_DOMAIN,
        bytes,
    ))
}

impl EulerScientificContract {
    /// Bind exact generic artifacts and Euler policy without widening authority.
    ///
    /// A lookalike local context cannot cross the concrete generic boundary:
    ///
    /// ```compile_fail
    /// use fs_evidence::vv::ContextOfUse;
    /// struct CounterfeitContextOfUse;
    /// fn requires_generic_context(_: ContextOfUse) {}
    /// requires_generic_context(CounterfeitContextOfUse);
    /// ```
    #[allow(clippy::too_many_lines)] // The composite invariant audit must precede publication.
    pub fn try_new(
        context: ContextOfUse,
        extension: EulerContextExtension,
        claim_graph: EulerClaimGraph,
        no_claims: NoClaimBoundary,
        owner_matrix: OwnerMatrix,
    ) -> Result<Self, ContractError> {
        if claim_graph.claims.len() != MAX_EULER_CLAIMS
            || !EULER_CLAIM_REGISTRY
                .into_iter()
                .all(|kind| claim_graph.claims.contains_key(&kind))
        {
            return Err(ContractError::new(
                "EulerContractClaimSetIncomplete",
                "the v1 scientific contract requires all nine claim kinds",
            ));
        }
        let context_bytes = context
            .canonical_bytes()
            .map_err(|error| ContractError::new("EulerContractContextCodec", error.to_string()))?;
        let decoded = ContextOfUse::from_canonical_bytes(&context_bytes)
            .map_err(|error| ContractError::new("EulerContractContextCodec", error.to_string()))?;
        if decoded != context {
            return Err(ContractError::new(
                "EulerContractContextNonCanonical",
                "ContextOfUse does not survive its generic canonical round trip",
            ));
        }
        let context_hash = context
            .content_hash()
            .map_err(|error| ContractError::new("EulerContractContextCodec", error.to_string()))?;
        if context.applicability().numeric().is_empty()
            || context.applicability().categorical().is_empty()
        {
            return Err(ContractError::new(
                "EulerContractApplicabilityIncomplete",
                "the physical Context of Use needs numeric and categorical applicability axes",
            ));
        }
        if context.applicability_policy() != ApplicabilityPolicy::Refuse {
            return Err(ContractError::new(
                "EulerContractApplicabilityPolicyMismatch",
                "the v1 Euler Context of Use must refuse every out-of-domain result",
            ));
        }
        let frame_axis = context
            .applicability()
            .categorical()
            .iter()
            .find(|(axis, _)| axis.as_str() == "observation-frame")
            .map(|(_, axis)| axis)
            .ok_or_else(|| {
                ContractError::new(
                    "EulerContractFrameMissing",
                    "the Context of Use must declare an observation-frame axis",
                )
            })?;
        if !frame_axis.allowed().contains(extension.observation_frame()) {
            return Err(ContractError::new(
                "EulerContractFrameMismatch",
                "the Euler observation frame is outside the Context applicability axis",
            ));
        }
        {
            // The generic header constructor retains units in canonical sorted
            // order, so binary search keeps this cross-object audit bounded even
            // for the generic schema's maximum unit/axis cardinalities.
            let declared_units = context.header().units();
            for qoi in context.qois().values() {
                if declared_units.binary_search(qoi.unit()).is_err() {
                    return Err(ContractError::new(
                        "EulerContractUnitCoverage",
                        format!(
                            "QoI {} unit {} is absent from the Context Five-Explicits header",
                            qoi.id(),
                            qoi.unit()
                        ),
                    ));
                }
            }
            for (axis_id, axis) in context.applicability().numeric() {
                if declared_units.binary_search(axis.unit()).is_err() {
                    return Err(ContractError::new(
                        "EulerContractUnitCoverage",
                        format!(
                            "numeric applicability axis {axis_id} unit {} is absent from the Context Five-Explicits header",
                            axis.unit()
                        ),
                    ));
                }
            }
        }
        let mut referenced_qois = BTreeSet::new();
        for claim in claim_graph.claims.values() {
            for qoi in &claim.campaign.qois {
                if !context.qois().contains_key(qoi) {
                    return Err(ContractError::new(
                        "EulerContractDanglingQoi",
                        format!("claim {} references unknown QoI {qoi}", claim.kind.id()),
                    ));
                }
                referenced_qois.insert(qoi.clone());
            }
        }
        let context_qois = context.qois().keys().cloned().collect::<BTreeSet<_>>();
        if referenced_qois != context_qois {
            return Err(ContractError::new(
                "EulerContractQoiClosure",
                "the frozen claim graph must cover every Context QoI exactly by identity",
            ));
        }
        if no_claims.entries().len() > MAX_EULER_NO_CLAIMS {
            return Err(ContractError::new(
                "EulerContractNoClaimCardinality",
                format!("scientific contract accepts at most {MAX_EULER_NO_CLAIMS} no-claim rows"),
            ));
        }
        // `NoClaimBoundary` is owned by fs-govern and deliberately applies its
        // own whitespace canonicalization.  Re-admit every resulting row under
        // this transport's stricter text contract before publishing an Euler
        // identity; otherwise a control character accepted upstream could be
        // encoded here and then refused by our own canonical decoder.
        for boundary in no_claims.entries() {
            checked_text("scientific_contract.no_claim", boundary.clone())?;
        }
        let no_claim_set = no_claims
            .entries()
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>();
        let missing_no_claims = CORE_NO_CLAIMS
            .into_iter()
            .filter(|boundary| !no_claim_set.contains(boundary))
            .collect::<Vec<_>>();
        if !missing_no_claims.is_empty() {
            return Err(ContractError::new(
                "EulerContractNoClaimBoundaryMissing",
                format!("missing binding core no-claims: {missing_no_claims:?}"),
            ));
        }
        let canonical_bytes = encode_contract_components(
            EULER_CONTRACT_SCHEMA_VERSION,
            &context_bytes,
            &extension,
            &claim_graph,
            &no_claims,
            &owner_matrix,
        )?;
        let identity = scientific_contract_identity(&canonical_bytes);
        Ok(Self {
            schema_version: EULER_CONTRACT_SCHEMA_VERSION,
            context,
            context_bytes,
            context_hash,
            extension,
            claim_graph,
            no_claims,
            owner_matrix,
            identity,
        })
    }

    #[must_use]
    /// Euler composite schema version.
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    /// Exact generic fs-evidence Context of Use.
    pub const fn context(&self) -> &ContextOfUse {
        &self.context
    }

    #[must_use]
    /// Exact retained generic Context canonical bytes.
    pub fn context_canonical_bytes(&self) -> &[u8] {
        &self.context_bytes
    }

    #[must_use]
    /// Generic Context content identity under its native domain.
    pub const fn context_hash(&self) -> ContentHash {
        self.context_hash
    }

    #[must_use]
    /// Euler-only Context extension.
    pub const fn extension(&self) -> &EulerContextExtension {
        &self.extension
    }

    #[must_use]
    /// Closed claim dependency graph.
    pub const fn claim_graph(&self) -> &EulerClaimGraph {
        &self.claim_graph
    }

    #[must_use]
    /// Binding generic no-claim set.
    pub const fn no_claims(&self) -> &NoClaimBoundary {
        &self.no_claims
    }

    #[must_use]
    /// Exact generic/Euler artifact owner matrix.
    pub const fn owner_matrix(&self) -> &OwnerMatrix {
        &self.owner_matrix
    }

    #[must_use]
    /// Domain-separated identity of every semantic contract field.
    pub const fn identity(&self) -> ContractIdentity {
        self.identity
    }

    /// V1 has no predecessor. Unknown, prior, and future versions refuse.
    pub fn migration_policy(source_version: u32) -> Result<(), ContractError> {
        if source_version == EULER_CONTRACT_SCHEMA_VERSION {
            Ok(())
        } else {
            Err(ContractError::new(
                "EulerContractUnsupportedVersion",
                format!(
                    "schema {source_version} is unsupported; v1 has no predecessor and no authority-preserving migration"
                ),
            ))
        }
    }

    /// Canonical composite preimage. The context retains its generic bytes and domain.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ContractError> {
        encode_contract_components(
            self.schema_version,
            &self.context_bytes,
            &self.extension,
            &self.claim_graph,
            &self.no_claims,
            &self.owner_matrix,
        )
    }

    /// Decode, validate, and re-encode the entire v1 contract as a fixed point.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, ContractError> {
        decode_contract(bytes)
    }
}

fn encode_contract_components(
    schema_version: u32,
    context_bytes: &[u8],
    extension: &EulerContextExtension,
    claim_graph: &EulerClaimGraph,
    no_claims: &NoClaimBoundary,
    owner_matrix: &OwnerMatrix,
) -> Result<Vec<u8>, ContractError> {
    let mut writer = CanonicalWriter::default();
    writer.bytes.extend_from_slice(CONTRACT_MAGIC);
    writer.u32(schema_version);
    writer.u32(VV_SCHEMA_VERSION);
    writer.string(VV_ARTIFACT_FAMILY)?;
    writer.u32(EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1);
    writer.u32(AUTHORITY_ALGEBRA_VERSION);
    writer.blob(context_bytes)?;
    extension.encode(&mut writer)?;
    writer.blob(&claim_graph.canonical_bytes()?)?;
    writer.usize(no_claims.entries().len())?;
    for boundary in no_claims.entries() {
        writer.string(boundary)?;
    }
    let owner_matrix_bytes = owner_matrix.canonical_bytes()?;
    let owner_matrix_identity = owner_matrix_identity(&owner_matrix_bytes);
    writer
        .bytes
        .extend_from_slice(owner_matrix_identity.as_hash().as_bytes());
    writer.blob(&owner_matrix_bytes)?;
    if writer.bytes.len() > MAX_EULER_CONTRACT_BYTES {
        return Err(ContractError::new(
            "EulerContractTooLarge",
            "canonical scientific contract exceeds its byte budget",
        ));
    }
    Ok(writer.bytes)
}

#[allow(clippy::too_many_lines)]
fn decode_contract(bytes: &[u8]) -> Result<EulerScientificContract, ContractError> {
    if bytes.len() > MAX_EULER_CONTRACT_BYTES {
        return Err(ContractError::new(
            "EulerContractTooLarge",
            "scientific-contract transport exceeds its byte budget",
        ));
    }
    let mut reader = CanonicalReader::new(
        bytes,
        "EulerContractMalformedCanonical",
        "scientific contract",
    );
    if reader.take(CONTRACT_MAGIC.len())? != CONTRACT_MAGIC {
        return Err(ContractError::new(
            "EulerContractMalformedCanonical",
            "scientific-contract magic is invalid",
        ));
    }
    let schema_version = reader.u32()?;
    EulerScientificContract::migration_policy(schema_version)?;
    let vv_version = reader.u32()?;
    if vv_version != VV_SCHEMA_VERSION || reader.string()? != VV_ARTIFACT_FAMILY {
        return Err(ContractError::new(
            "EulerContractGenericSchemaMismatch",
            "embedded fs-evidence V&V schema identity is stale or counterfeit",
        ));
    }
    if reader.u32()? != EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1 {
        return Err(ContractError::new(
            "EulerContractGenericSchemaMismatch",
            "embedded fs-ir campaign schema version is unsupported",
        ));
    }
    if reader.u32()? != AUTHORITY_ALGEBRA_VERSION {
        return Err(ContractError::new(
            "EulerContractGenericSchemaMismatch",
            "embedded fs-govern authority schema version is unsupported",
        ));
    }
    let context_bytes = reader.blob(MAX_EULER_CONTRACT_BYTES)?;
    let context = ContextOfUse::from_canonical_bytes(context_bytes)
        .map_err(|error| ContractError::new("EulerContractContextCodec", error.to_string()))?;

    let user_count = reader.bounded_len("context.users", MAX_EULER_CLAIMS * 8, false)?;
    let mut users = Vec::with_capacity(user_count);
    for _ in 0..user_count {
        users.push(reader.string()?);
    }
    let apparatus_population = reader.string()?;
    let environment_population = reader.string()?;
    let observation_frame = reader.string()?;
    let alternative_count =
        reader.bounded_len("context.decision_alternatives", MAX_EULER_CLAIMS * 8, false)?;
    let mut decision_alternatives = Vec::with_capacity(alternative_count);
    for _ in 0..alternative_count {
        decision_alternatives.push(reader.string()?);
    }
    let risk_count = reader.bounded_len("context.risks", MAX_EULER_CLAIMS * 2, false)?;
    let mut risks = Vec::with_capacity(risk_count);
    for _ in 0..risk_count {
        let code = reader.string()?;
        let consequence = reader.string()?;
        let severity = reader.u8()?;
        let claim_count = reader.bounded_len("risk.affected_claims", MAX_EULER_CLAIMS, false)?;
        let mut affected_claims = Vec::with_capacity(claim_count);
        for _ in 0..claim_count {
            affected_claims.push(EulerClaimKind::from_tag(reader.u8()?)?);
        }
        let decision_alternative = reader.string()?;
        risks.push(ScientificRisk::try_new(
            code,
            consequence,
            severity,
            affected_claims,
            decision_alternative,
        )?);
    }
    let source_count =
        reader.bounded_len("context.hypothesis_sources", MAX_EULER_CLAIMS * 2, false)?;
    let mut hypothesis_sources = Vec::with_capacity(source_count);
    for _ in 0..source_count {
        hypothesis_sources.push(HypothesisSource::from_canonical_parts(
            reader.string()?,
            reader.string()?,
            reader.hash()?,
        )?);
    }
    let extension = EulerContextExtension::try_new(
        users,
        apparatus_population,
        environment_population,
        observation_frame,
        decision_alternatives,
        risks,
        hypothesis_sources,
    )?;

    let graph = EulerClaimGraph::from_canonical_bytes(reader.blob(MAX_EULER_GRAPH_BYTES)?)?;
    let no_claim_count = reader.bounded_len("contract.no_claims", MAX_EULER_NO_CLAIMS, false)?;
    let mut no_claim_text = Vec::with_capacity(no_claim_count);
    for _ in 0..no_claim_count {
        no_claim_text.push(reader.string()?);
    }
    let no_claim_refs = no_claim_text.iter().map(String::as_str).collect::<Vec<_>>();
    let no_claims = NoClaimBoundary::new(&no_claim_refs)
        .map_err(|error| ContractError::new("EulerContractNoClaimCodec", error.to_string()))?;

    let encoded_owner_matrix_identity = OwnerMatrixIdentity(reader.hash()?);
    let owner_matrix = OwnerMatrix::from_canonical_bytes(reader.blob(MAX_OWNER_MATRIX_BYTES)?)?;
    if encoded_owner_matrix_identity != owner_matrix.identity() {
        return Err(ContractError::new(
            "EulerContractOwnerMatrixIdentityMismatch",
            "embedded owner-matrix identity does not match its canonical registry bytes",
        ));
    }
    if !reader.is_finished() {
        return Err(ContractError::new(
            "EulerContractMalformedCanonical",
            "scientific-contract transport has trailing bytes",
        ));
    }
    let decoded =
        EulerScientificContract::try_new(context, extension, graph, no_claims, owner_matrix)?;
    if decoded.canonical_bytes()?.as_slice() != bytes {
        return Err(ContractError::new(
            "EulerContractNonCanonical",
            "scientific-contract bytes are not a canonical fixed point",
        ));
    }
    Ok(decoded)
}

#[derive(Default)]
struct CanonicalWriter {
    bytes: Vec<u8>,
}

impl CanonicalWriter {
    fn u8(&mut self, value: u8) {
        self.bytes.push(value);
    }

    fn u32(&mut self, value: u32) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn usize(&mut self, value: usize) -> Result<(), ContractError> {
        let value = u32::try_from(value).map_err(|_| {
            ContractError::new(
                "EulerContractCardinality",
                "collection length does not fit the canonical transport",
            )
        })?;
        self.u32(value);
        Ok(())
    }

    fn blob(&mut self, value: &[u8]) -> Result<(), ContractError> {
        self.usize(value.len())?;
        self.bytes.extend_from_slice(value);
        Ok(())
    }

    fn string(&mut self, value: &str) -> Result<(), ContractError> {
        self.blob(value.as_bytes())
    }

    fn strings(&mut self, values: &[String]) -> Result<(), ContractError> {
        self.usize(values.len())?;
        for value in values {
            self.string(value)?;
        }
        Ok(())
    }
}

struct CanonicalReader<'a> {
    bytes: &'a [u8],
    offset: usize,
    malformed_code: &'static str,
    artifact: &'static str,
}

impl<'a> CanonicalReader<'a> {
    const fn new(bytes: &'a [u8], malformed_code: &'static str, artifact: &'static str) -> Self {
        Self {
            bytes,
            offset: 0,
            malformed_code,
            artifact,
        }
    }

    fn take(&mut self, len: usize) -> Result<&'a [u8], ContractError> {
        let end = self.offset.checked_add(len).ok_or_else(|| {
            ContractError::new(
                self.malformed_code,
                format!("{} canonical transport offset overflow", self.artifact),
            )
        })?;
        let value = self.bytes.get(self.offset..end).ok_or_else(|| {
            ContractError::new(
                self.malformed_code,
                format!(
                    "truncated {} canonical transport at byte {}",
                    self.artifact, self.offset
                ),
            )
        })?;
        self.offset = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, ContractError> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32, ContractError> {
        let mut bytes = [0_u8; 4];
        bytes.copy_from_slice(self.take(4)?);
        Ok(u32::from_le_bytes(bytes))
    }

    fn usize(&mut self) -> Result<usize, ContractError> {
        usize::try_from(self.u32()?).map_err(|_| {
            ContractError::new(
                "EulerContractCardinality",
                "canonical collection length does not fit this platform",
            )
        })
    }

    fn bounded_len(
        &mut self,
        field: &'static str,
        maximum: usize,
        allow_empty: bool,
    ) -> Result<usize, ContractError> {
        let value = self.usize()?;
        if value > maximum || (!allow_empty && value == 0) {
            return Err(ContractError::new(
                "EulerContractCardinality",
                format!("{field} has an invalid canonical item count {value}"),
            ));
        }
        Ok(value)
    }

    fn blob(&mut self, maximum: usize) -> Result<&'a [u8], ContractError> {
        let len = self.usize()?;
        if len > maximum {
            return Err(ContractError::new(
                "EulerContractCardinality",
                format!("canonical blob length {len} exceeds limit {maximum}"),
            ));
        }
        self.take(len)
    }

    fn hash(&mut self) -> Result<ContentHash, ContractError> {
        let mut bytes = [0_u8; 32];
        bytes.copy_from_slice(self.take(32)?);
        Ok(ContentHash(bytes))
    }

    fn string(&mut self) -> Result<String, ContractError> {
        let len = self.usize()?;
        if len == 0 || len > MAX_EULER_TEXT_BYTES {
            return Err(ContractError::new(
                self.malformed_code,
                format!("{} canonical string has an invalid length", self.artifact),
            ));
        }
        let bytes = self.take(len)?;
        let value = core::str::from_utf8(bytes).map_err(|_| {
            ContractError::new(
                self.malformed_code,
                format!("{} canonical string is not UTF-8", self.artifact),
            )
        })?;
        checked_text("canonical.string", value)
    }

    const fn is_finished(&self) -> bool {
        self.offset == self.bytes.len()
    }
}

#[allow(dead_code)]
fn classify_hypothesis_source_identity_fields(source: &HypothesisSource) {
    let HypothesisSource {
        id: _,
        locator: _,
        declaration_hash: _,
    } = source;
}

#[allow(dead_code)]
fn classify_euler_claim_graph_identity_fields(graph: &EulerClaimGraph, claim: &EulerClaimSpec) {
    let EulerClaimGraph {
        claims: _,
        dependencies: _,
    } = graph;
    let EulerClaimSpec {
        kind: _,
        campaign: _,
        requirements: _,
    } = claim;
}

#[allow(dead_code)]
fn classify_owner_matrix_identity_fields(matrix: &OwnerMatrix, row: &OwnerRow) {
    let OwnerMatrix {
        schema_version: _,
        rows: _,
        identity: _,
    } = matrix;
    let OwnerRow {
        role: _,
        owner_crate: _,
        source_schema: _,
        authority_ceiling: _,
    } = row;
}

#[allow(dead_code)]
fn classify_euler_contract_identity_fields(
    contract: &EulerScientificContract,
    extension: &EulerContextExtension,
    risk: &ScientificRisk,
) {
    let EulerScientificContract {
        schema_version: _,
        context: _,
        context_bytes: _,
        context_hash: _,
        extension: _,
        claim_graph: _,
        no_claims: _,
        owner_matrix: _,
        identity: _,
    } = contract;
    let EulerContextExtension {
        users: _,
        apparatus_population: _,
        environment_population: _,
        observation_frame: _,
        decision_alternatives: _,
        risks: _,
        hypothesis_sources: _,
    } = extension;
    let ScientificRisk {
        code: _,
        consequence: _,
        severity: _,
        affected_claims: _,
        decision_alternative: _,
    } = risk;
}

/// Stable family assertion used by manifest/owner tests.
#[must_use]
pub const fn context_artifact_kind() -> ArtifactKind {
    ArtifactKind::ContextOfUse
}

#[cfg(test)]
mod tests {
    #![allow(clippy::too_many_lines)] // Mutation batteries intentionally keep one field oracle together.

    use super::*;
    use crate::build_frozen_contract;

    // These batteries prove that each declared preimage field is bound by its
    // identity. Some candidates intentionally represent an incompatible codec
    // or inadmissible source object; constructor/refusal tests separately prove
    // which source objects may enter the canonical path.
    fn assert_hash_moved(field: &str, base: ContentHash, candidate: ContentHash) {
        assert_ne!(
            candidate, base,
            "semantic identity field {field} did not move the content hash"
        );
    }

    #[test]
    fn hypothesis_source_identity_semantic_fields_move_independently() {
        fn preimage(
            version: u32,
            id: &str,
            locator: &str,
            big_endian: bool,
            framed: bool,
            locator_first: bool,
        ) -> Vec<u8> {
            let mut bytes = Vec::new();
            let version_bytes = if big_endian {
                version.to_be_bytes()
            } else {
                version.to_le_bytes()
            };
            bytes.extend_from_slice(&version_bytes);
            let mut push = |value: &str| {
                if framed {
                    let length = u64::try_from(value.len()).expect("bounded source field");
                    let length_bytes = if big_endian {
                        length.to_be_bytes()
                    } else {
                        length.to_le_bytes()
                    };
                    bytes.extend_from_slice(&length_bytes);
                }
                bytes.extend_from_slice(value.as_bytes());
            };
            if locator_first {
                push(locator);
                push(id);
            } else {
                push(id);
                push(locator);
            }
            bytes
        }

        let id = "source-ab";
        let locator = "locator-c";
        let base = hypothesis_source_declaration_hash_with_schema(
            EULER_CONTRACT_SCHEMA_VERSION,
            HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
            id,
            locator,
        );
        let canonical = preimage(
            EULER_CONTRACT_SCHEMA_VERSION,
            id,
            locator,
            false,
            true,
            false,
        );
        assert_eq!(
            base,
            fs_blake3::hash_domain(HYPOTHESIS_SOURCE_DECLARATION_DOMAIN, &canonical)
        );
        for (field, candidate) in [
            (
                "identity-domain",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.other-source.v1",
                    &canonical,
                ),
            ),
            (
                "identity-version",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.hypothesis-source-declaration.v2",
                    &canonical,
                ),
            ),
            (
                "contract-schema-version",
                hypothesis_source_declaration_hash_with_schema(
                    EULER_CONTRACT_SCHEMA_VERSION + 1,
                    HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
                    id,
                    locator,
                ),
            ),
            (
                "source-id",
                hypothesis_source_declaration_hash_with_schema(
                    EULER_CONTRACT_SCHEMA_VERSION,
                    HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
                    "source-ac",
                    locator,
                ),
            ),
            (
                "source-locator",
                hypothesis_source_declaration_hash_with_schema(
                    EULER_CONTRACT_SCHEMA_VERSION,
                    HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
                    id,
                    "locator-d",
                ),
            ),
            (
                "canonical-field-order",
                fs_blake3::hash_domain(
                    HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
                    &preimage(
                        EULER_CONTRACT_SCHEMA_VERSION,
                        id,
                        locator,
                        false,
                        true,
                        true,
                    ),
                ),
            ),
            (
                "length-framing",
                fs_blake3::hash_domain(
                    HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
                    &preimage(
                        EULER_CONTRACT_SCHEMA_VERSION,
                        id,
                        locator,
                        false,
                        false,
                        false,
                    ),
                ),
            ),
            (
                "fixed-numeric-little-endian",
                fs_blake3::hash_domain(
                    HYPOTHESIS_SOURCE_DECLARATION_DOMAIN,
                    &preimage(
                        EULER_CONTRACT_SCHEMA_VERSION,
                        id,
                        locator,
                        true,
                        true,
                        false,
                    ),
                ),
            ),
        ] {
            assert_hash_moved(field, base, candidate);
        }
    }

    #[test]
    fn claim_graph_identity_semantic_fields_move_independently() {
        let contract = build_frozen_contract().expect("frozen contract");
        let graph = contract.claim_graph();
        let bytes = graph.canonical_bytes().expect("graph bytes");
        let base = graph.content_hash().expect("graph identity");
        assert_eq!(
            base,
            fs_blake3::hash_domain(EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, &bytes)
        );

        let mut fewer_claims = graph.clone();
        let first_claim = *fewer_claims.claims.keys().next().expect("claim registry");
        fewer_claims.claims.remove(&first_claim);
        let mut fewer_dependencies = graph.clone();
        fewer_dependencies.dependencies.pop();

        let mut wrong_magic = bytes.clone();
        wrong_magic[0] ^= 1;
        let mut wrong_version = bytes.clone();
        wrong_version[8..12].copy_from_slice(&2_u32.to_le_bytes());
        let mut big_endian_version = bytes.clone();
        big_endian_version[8..12].copy_from_slice(&1_u32.to_be_bytes());
        let mut reordered = Vec::with_capacity(bytes.len());
        reordered.extend_from_slice(&bytes[..8]);
        reordered.extend_from_slice(&bytes[12..16]);
        reordered.extend_from_slice(&bytes[8..12]);
        reordered.extend_from_slice(&bytes[16..]);
        let mut unframed = Vec::with_capacity(bytes.len() - 4);
        unframed.extend_from_slice(&bytes[..12]);
        unframed.extend_from_slice(&bytes[16..]);

        for (field, candidate) in [
            (
                "identity-domain",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.other-claim-graph.v1",
                    &bytes,
                ),
            ),
            (
                "identity-version",
                fs_blake3::hash_domain("org.frankensim.fs-euler-disc-e2e.claim-graph.v2", &bytes),
            ),
            (
                "transport-magic",
                fs_blake3::hash_domain(EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, &wrong_magic),
            ),
            (
                "claim-policy-schema-version",
                fs_blake3::hash_domain(EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, &wrong_version),
            ),
            (
                "canonical-field-order",
                fs_blake3::hash_domain(EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, &reordered),
            ),
            (
                "length-framing",
                fs_blake3::hash_domain(EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, &unframed),
            ),
            (
                "fixed-numeric-little-endian",
                fs_blake3::hash_domain(EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, &big_endian_version),
            ),
            (
                "claim-registry",
                fewer_claims.content_hash().expect("changed claim registry"),
            ),
            (
                "dependency-registry",
                fewer_dependencies
                    .content_hash()
                    .expect("changed dependency registry"),
            ),
        ] {
            assert_hash_moved(field, base, candidate);
        }
    }

    #[test]
    fn owner_matrix_identity_semantic_fields_move_independently() {
        assert_eq!(OwnerRole::OwnerMatrixRegistry.tag(), 14);
        assert_eq!(OwnerRole::AggregateQoiDerivationReceipt.tag(), 15);
        let contract = build_frozen_contract().expect("frozen contract");
        let matrix = contract.owner_matrix();
        let bytes = matrix.canonical_bytes().expect("owner-matrix bytes");
        let base = matrix.identity().as_hash();
        assert_eq!(
            base,
            fs_blake3::hash_domain(EULER_OWNER_MATRIX_IDENTITY_DOMAIN, &bytes)
        );
        assert_eq!(
            OwnerMatrix::from_canonical_bytes(&bytes).expect("owner-matrix fixed point"),
            *matrix
        );

        let mut wrong_magic = bytes.clone();
        wrong_magic[0] ^= 1;
        let mut wrong_version = bytes.clone();
        wrong_version[8..12]
            .copy_from_slice(&(EULER_OWNER_MATRIX_SCHEMA_VERSION + 1).to_le_bytes());
        let mut big_endian_version = bytes.clone();
        big_endian_version[8..12].copy_from_slice(&EULER_OWNER_MATRIX_SCHEMA_VERSION.to_be_bytes());
        let mut unframed = Vec::with_capacity(bytes.len() - 4);
        unframed.extend_from_slice(&bytes[..12]);
        unframed.extend_from_slice(&bytes[16..]);

        let mut reader =
            CanonicalReader::new(&bytes, "EulerOwnerMatrixMalformedCanonical", "owner matrix");
        reader.take(OWNER_MATRIX_MAGIC.len()).expect("matrix magic");
        reader.u32().expect("matrix version");
        let row_count = reader.usize().expect("row count");
        let mut row_spans = Vec::with_capacity(row_count);
        for _ in 0..row_count {
            let start = reader.offset;
            reader.u8().expect("role tag");
            reader.string().expect("owner crate");
            reader.string().expect("source schema");
            reader.u8().expect("authority tag");
            row_spans.push((start, reader.offset));
        }
        let mut reordered = Vec::with_capacity(bytes.len());
        reordered.extend_from_slice(&bytes[..row_spans[0].0]);
        reordered.extend_from_slice(&bytes[row_spans[1].0..row_spans[1].1]);
        reordered.extend_from_slice(&bytes[row_spans[0].0..row_spans[0].1]);
        reordered.extend_from_slice(&bytes[row_spans[1].1..]);

        for (field, candidate_bytes) in [
            ("transport-magic", wrong_magic),
            ("owner-matrix-schema-version", wrong_version),
            ("canonical-field-order", reordered),
            ("length-framing", unframed),
            ("fixed-numeric-little-endian", big_endian_version),
        ] {
            assert_hash_moved(
                field,
                base,
                fs_blake3::hash_domain(EULER_OWNER_MATRIX_IDENTITY_DOMAIN, &candidate_bytes),
            );
            assert!(
                OwnerMatrix::from_canonical_bytes(&candidate_bytes).is_err(),
                "{field} mutation must refuse transport admission"
            );
        }
        for (field, candidate) in [
            (
                "identity-domain",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.other-owner-matrix.v1",
                    &bytes,
                ),
            ),
            (
                "identity-version",
                fs_blake3::hash_domain("org.frankensim.fs-euler-disc-e2e.owner-matrix.v2", &bytes),
            ),
        ] {
            assert_hash_moved(field, base, candidate);
        }

        let alternate_ceiling = |expected: AuthorityCeiling| {
            if expected == AuthorityCeiling::HypothesisOnly {
                AuthorityCeiling::StructuralCheckOnly
            } else {
                AuthorityCeiling::HypothesisOnly
            }
        };
        for (index, role) in OwnerRole::ALL.into_iter().enumerate() {
            let next_role = OwnerRole::ALL[(index + 1) % OwnerRole::ALL.len()];
            let mutations = [
                ("owner-role-tags", 0_u8),
                ("owner-crate-routing", 1_u8),
                ("source-schema-routing-address", 2_u8),
                ("authority-ceiling-tags", 3_u8),
            ];
            for (field, mutation) in mutations {
                let mut changed = matrix.clone();
                let row = changed.rows.get_mut(&role).expect("owner row");
                match mutation {
                    0 => row.role = next_role,
                    1 => row.owner_crate.push_str("-changed"),
                    2 => row.source_schema.push_str("-changed"),
                    3 => row.authority_ceiling = alternate_ceiling(row.authority_ceiling),
                    _ => unreachable!(),
                }
                let candidate_bytes = encode_owner_matrix_components(
                    EULER_OWNER_MATRIX_SCHEMA_VERSION,
                    &changed.rows,
                )
                .expect("bounded mutation bytes");
                assert_hash_moved(
                    &format!("{field}:{}", role.tag()),
                    base,
                    fs_blake3::hash_domain(EULER_OWNER_MATRIX_IDENTITY_DOMAIN, &candidate_bytes),
                );
                assert!(
                    OwnerMatrix::from_canonical_bytes(&candidate_bytes).is_err(),
                    "{field} mutation for {role:?} must refuse"
                );
            }
        }
    }

    #[test]
    fn scientific_contract_identity_semantic_fields_move_independently() {
        let contract = build_frozen_contract().expect("frozen contract");
        let bytes = contract.canonical_bytes().expect("contract bytes");
        let base = contract.identity().as_hash();
        assert_eq!(base, scientific_contract_identity(&bytes).as_hash());
        let encode = |version: u32,
                      context_bytes: &[u8],
                      extension: &EulerContextExtension,
                      no_claims: &NoClaimBoundary,
                      owners: &OwnerMatrix| {
            encode_contract_components(
                version,
                context_bytes,
                extension,
                contract.claim_graph(),
                no_claims,
                owners,
            )
        };

        let mut changed_context_bytes = contract.context_bytes.clone();
        changed_context_bytes.push(0);
        let mut changed_users = contract.extension.clone();
        changed_users.users.push("identity-test-user".to_owned());
        let mut changed_apparatus = contract.extension.clone();
        changed_apparatus.apparatus_population.push_str(" changed");
        let mut changed_environment = contract.extension.clone();
        changed_environment
            .environment_population
            .push_str(" changed");
        let mut changed_frame = contract.extension.clone();
        changed_frame.observation_frame.push_str("-changed");
        let mut changed_alternatives = contract.extension.clone();
        changed_alternatives
            .decision_alternatives
            .push("identity-test-alternative".to_owned());
        let mut changed_risks = contract.extension.clone();
        changed_risks.risks[0].consequence.push_str(" changed");
        let mut changed_hypothesis_sources = contract.extension.clone();
        let original_source = changed_hypothesis_sources
            .hypothesis_sources
            .first()
            .expect("frozen hypothesis source")
            .clone();
        changed_hypothesis_sources.hypothesis_sources[0] = HypothesisSource::try_new(
            original_source.id(),
            format!("{}#identity-test", original_source.locator()),
        )
        .expect("changed hypothesis source");
        let mut changed_claim_graph = contract.claim_graph.clone();
        changed_claim_graph
            .claims
            .values_mut()
            .next()
            .expect("frozen graph claim")
            .campaign
            .hypothesis
            .push_str(" Identity-only graph mutation.");
        let mut changed_owners = contract.owner_matrix.clone();
        changed_owners
            .rows
            .values_mut()
            .next()
            .expect("owner row")
            .source_schema
            .push_str("-changed");

        let mut no_claim_entries = contract.no_claims.entries().to_vec();
        no_claim_entries.push("Identity-only additional no-claim.".to_owned());
        let no_claim_refs = no_claim_entries
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>();
        let changed_no_claims = NoClaimBoundary::new(&no_claim_refs).expect("changed no claims");

        let family_end = 20 + VV_ARTIFACT_FAMILY.len();
        let campaign_offset = family_end;
        let authority_offset = campaign_offset + 4;
        let mut wrong_magic = bytes.clone();
        wrong_magic[0] ^= 1;
        let mut wrong_vv = bytes.clone();
        wrong_vv[12..16].copy_from_slice(&(VV_SCHEMA_VERSION + 1).to_le_bytes());
        let mut wrong_family = bytes.clone();
        wrong_family[20] ^= 1;
        let mut wrong_campaign = bytes.clone();
        wrong_campaign[campaign_offset..authority_offset]
            .copy_from_slice(&(EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1 + 1).to_le_bytes());
        let mut wrong_authority = bytes.clone();
        wrong_authority[authority_offset..authority_offset + 4]
            .copy_from_slice(&(AUTHORITY_ALGEBRA_VERSION + 1).to_le_bytes());
        let mut reordered = bytes.clone();
        reordered[campaign_offset..authority_offset]
            .copy_from_slice(&AUTHORITY_ALGEBRA_VERSION.to_le_bytes());
        reordered[authority_offset..authority_offset + 4]
            .copy_from_slice(&EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1.to_le_bytes());
        let mut unframed = Vec::with_capacity(bytes.len() - 4);
        unframed.extend_from_slice(&bytes[..16]);
        unframed.extend_from_slice(&bytes[20..]);
        let mut big_endian = bytes.clone();
        big_endian[authority_offset..authority_offset + 4]
            .copy_from_slice(&AUTHORITY_ALGEBRA_VERSION.to_be_bytes());

        let variants = vec![
            (
                "identity-domain",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.other-contract.v1",
                    &bytes,
                ),
            ),
            (
                "identity-version",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.scientific-contract.v2",
                    &bytes,
                ),
            ),
            (
                "transport-magic",
                fs_blake3::hash_domain(EULER_CONTRACT_IDENTITY_DOMAIN, &wrong_magic),
            ),
            (
                "embedded-vv-schema-version",
                fs_blake3::hash_domain(EULER_CONTRACT_IDENTITY_DOMAIN, &wrong_vv),
            ),
            (
                "embedded-vv-artifact-family",
                fs_blake3::hash_domain(EULER_CONTRACT_IDENTITY_DOMAIN, &wrong_family),
            ),
            (
                "embedded-campaign-schema-version",
                fs_blake3::hash_domain(EULER_CONTRACT_IDENTITY_DOMAIN, &wrong_campaign),
            ),
            (
                "embedded-authority-version",
                fs_blake3::hash_domain(EULER_CONTRACT_IDENTITY_DOMAIN, &wrong_authority),
            ),
            (
                "canonical-field-order",
                fs_blake3::hash_domain(EULER_CONTRACT_IDENTITY_DOMAIN, &reordered),
            ),
            (
                "length-framing",
                fs_blake3::hash_domain(EULER_CONTRACT_IDENTITY_DOMAIN, &unframed),
            ),
            (
                "fixed-numeric-little-endian",
                fs_blake3::hash_domain(EULER_CONTRACT_IDENTITY_DOMAIN, &big_endian),
            ),
            (
                "contract-schema-version",
                scientific_contract_identity(
                    &encode(
                        EULER_CONTRACT_SCHEMA_VERSION + 1,
                        &contract.context_bytes,
                        &contract.extension,
                        &contract.no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed contract version"),
                )
                .as_hash(),
            ),
            (
                "context-canonical-bytes",
                scientific_contract_identity(
                    &encode(
                        contract.schema_version,
                        &changed_context_bytes,
                        &contract.extension,
                        &contract.no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed context bytes"),
                )
                .as_hash(),
            ),
            (
                "extension-users",
                scientific_contract_identity(
                    &encode(
                        contract.schema_version,
                        &contract.context_bytes,
                        &changed_users,
                        &contract.no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed users"),
                )
                .as_hash(),
            ),
            (
                "apparatus-population",
                scientific_contract_identity(
                    &encode(
                        contract.schema_version,
                        &contract.context_bytes,
                        &changed_apparatus,
                        &contract.no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed apparatus"),
                )
                .as_hash(),
            ),
            (
                "environment-population",
                scientific_contract_identity(
                    &encode(
                        contract.schema_version,
                        &contract.context_bytes,
                        &changed_environment,
                        &contract.no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed environment"),
                )
                .as_hash(),
            ),
            (
                "observation-frame",
                scientific_contract_identity(
                    &encode(
                        contract.schema_version,
                        &contract.context_bytes,
                        &changed_frame,
                        &contract.no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed frame"),
                )
                .as_hash(),
            ),
            (
                "decision-alternatives",
                scientific_contract_identity(
                    &encode(
                        contract.schema_version,
                        &contract.context_bytes,
                        &changed_alternatives,
                        &contract.no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed alternatives"),
                )
                .as_hash(),
            ),
            (
                "risk-registry",
                scientific_contract_identity(
                    &encode(
                        contract.schema_version,
                        &contract.context_bytes,
                        &changed_risks,
                        &contract.no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed risks"),
                )
                .as_hash(),
            ),
            (
                "hypothesis-source-declarations",
                scientific_contract_identity(
                    &encode(
                        contract.schema_version,
                        &contract.context_bytes,
                        &changed_hypothesis_sources,
                        &contract.no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed hypothesis sources"),
                )
                .as_hash(),
            ),
            (
                "claim-graph-canonical-bytes",
                scientific_contract_identity(
                    &encode_contract_components(
                        contract.schema_version,
                        &contract.context_bytes,
                        &contract.extension,
                        &changed_claim_graph,
                        &contract.no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed claim graph"),
                )
                .as_hash(),
            ),
            (
                "no-claim-boundary",
                scientific_contract_identity(
                    &encode(
                        contract.schema_version,
                        &contract.context_bytes,
                        &contract.extension,
                        &changed_no_claims,
                        &contract.owner_matrix,
                    )
                    .expect("changed no claims"),
                )
                .as_hash(),
            ),
            (
                "owner-matrix",
                scientific_contract_identity(
                    &encode(
                        contract.schema_version,
                        &contract.context_bytes,
                        &contract.extension,
                        &contract.no_claims,
                        &changed_owners,
                    )
                    .expect("changed owner matrix"),
                )
                .as_hash(),
            ),
        ];
        for (field, candidate) in variants {
            assert_hash_moved(field, base, candidate);
        }
        for (field, hostile) in [
            ("VV schema version", &wrong_vv),
            ("VV artifact family", &wrong_family),
            ("campaign schema version", &wrong_campaign),
            ("authority algebra version", &wrong_authority),
        ] {
            let error = match EulerScientificContract::from_canonical_bytes(hostile) {
                Ok(_) => panic!("wrong embedded {field} must refuse contract admission"),
                Err(error) => error,
            };
            assert_eq!(
                error.code(),
                "EulerContractGenericSchemaMismatch",
                "wrong embedded {field} returned the wrong refusal"
            );
        }
    }

    #[test]
    fn scientific_contract_refuses_a_demoting_applicability_policy_at_both_boundaries() {
        let base = build_frozen_contract().expect("frozen contract");
        assert_eq!(
            base.context().applicability_policy(),
            ApplicabilityPolicy::Refuse
        );
        let demoting_context = ContextOfUse::try_new(
            base.context().header().clone(),
            base.context().decision().to_owned(),
            base.context().qois().values().cloned().collect(),
            base.context().applicability().clone(),
            ApplicabilityPolicy::Demote,
        )
        .expect("generic Context of Use permits an explicit demotion policy");

        let constructor_error = EulerScientificContract::try_new(
            demoting_context.clone(),
            base.extension().clone(),
            base.claim_graph().clone(),
            base.no_claims().clone(),
            base.owner_matrix().clone(),
        )
        .expect_err("the frozen Euler contract must refuse a demoting context");
        assert_eq!(
            constructor_error.code(),
            "EulerContractApplicabilityPolicyMismatch"
        );

        let hostile_bytes = encode_contract_components(
            EULER_CONTRACT_SCHEMA_VERSION,
            &demoting_context
                .canonical_bytes()
                .expect("generic demoting context bytes"),
            base.extension(),
            base.claim_graph(),
            base.no_claims(),
            base.owner_matrix(),
        )
        .expect("bounded hostile contract bytes");
        let decoder_error = EulerScientificContract::from_canonical_bytes(&hostile_bytes)
            .expect_err("the decoder must apply the same frozen applicability policy");
        assert_eq!(
            decoder_error.code(),
            "EulerContractApplicabilityPolicyMismatch"
        );
    }

    #[test]
    fn scientific_contract_no_claim_rows_share_constructor_and_decoder_text_rules() {
        let base = build_frozen_contract().expect("frozen contract");

        let mut valid_rows = CORE_NO_CLAIMS.to_vec();
        valid_rows.push("Additional bounded no-claim for fixed-point coverage.");
        let valid_no_claims =
            NoClaimBoundary::new(&valid_rows).expect("valid generic no-claim boundary");
        let valid = EulerScientificContract::try_new(
            base.context().clone(),
            base.extension().clone(),
            base.claim_graph().clone(),
            valid_no_claims,
            base.owner_matrix().clone(),
        )
        .expect("valid additional no-claim must publish a contract identity");
        let valid_bytes = valid.canonical_bytes().expect("valid contract bytes");
        assert_eq!(
            EulerScientificContract::from_canonical_bytes(&valid_bytes)
                .expect("valid additional no-claim must decode as a fixed point"),
            valid
        );

        let mut hostile_rows = CORE_NO_CLAIMS.to_vec();
        hostile_rows.push("Additional\0no-claim boundary.");
        let hostile_no_claims = NoClaimBoundary::new(&hostile_rows)
            .expect("the generic boundary currently permits an embedded NUL");
        let error = EulerScientificContract::try_new(
            base.context().clone(),
            base.extension().clone(),
            base.claim_graph().clone(),
            hostile_no_claims,
            base.owner_matrix().clone(),
        )
        .expect_err("constructor must refuse bytes that its own decoder rejects");
        assert_eq!(error.code(), "EulerContractInvalidText");
    }

    #[test]
    fn scientific_contract_identity_dependencies_exclude_opaque_owner_routes() {
        let dependency_row = EULER_SCIENTIFIC_CONTRACT_IDENTITY_SCHEMA_DECLARATION
            .iter()
            .find(|row| row.starts_with("schema_dependencies="))
            .expect("scientific-contract identity declaration has dependencies");
        assert_eq!(
            *dependency_row,
            "schema_dependencies=fs-euler-disc-e2e:claim-graph,fs-euler-disc-e2e:hypothesis-source,fs-euler-disc-e2e:owner-matrix,fs-evidence:vv-artifact"
        );
        assert!(
            !dependency_row.contains("fs-evidence:vv-schema-admission-receipt"),
            "owner-matrix routing addresses are opaque and must not create transitive identity dependencies"
        );
    }

    #[test]
    fn canonical_reader_reports_artifact_specific_truncation_diagnostics() {
        let contract = build_frozen_contract().expect("frozen contract");
        let owner_bytes = contract
            .owner_matrix()
            .canonical_bytes()
            .expect("owner-matrix bytes");
        for end in 0..owner_bytes.len() {
            let error = OwnerMatrix::from_canonical_bytes(&owner_bytes[..end])
                .expect_err("every proper owner-matrix prefix must refuse");
            assert_eq!(error.code(), "EulerOwnerMatrixMalformedCanonical");
            assert!(
                error.detail().contains("owner matrix"),
                "truncation at {end} produced a misleading detail: {}",
                error.detail()
            );
        }

        let graph_error = EulerClaimGraph::from_canonical_bytes(&[])
            .expect_err("an empty claim-graph transport must refuse");
        assert_eq!(graph_error.code(), "EulerContractMalformedCanonical");
        assert_eq!(
            graph_error.detail(),
            "truncated claim graph canonical transport at byte 0"
        );

        let contract_error = EulerScientificContract::from_canonical_bytes(&[])
            .expect_err("an empty scientific-contract transport must refuse");
        assert_eq!(contract_error.code(), "EulerContractMalformedCanonical");
        assert_eq!(
            contract_error.detail(),
            "truncated scientific contract canonical transport at byte 0"
        );
    }

    #[test]
    fn owner_matrix_decoder_uses_its_namespace_without_erasing_semantic_errors() {
        let contract = build_frozen_contract().expect("frozen contract");
        let bytes = contract
            .owner_matrix()
            .canonical_bytes()
            .expect("owner-matrix bytes");
        let mut reader =
            CanonicalReader::new(&bytes, "EulerOwnerMatrixMalformedCanonical", "owner matrix");
        reader.take(OWNER_MATRIX_MAGIC.len()).expect("matrix magic");
        reader.u32().expect("matrix version");
        reader.usize().expect("row count");
        let role_offset = reader.offset;
        reader.u8().expect("role tag");
        reader.string().expect("owner crate");
        reader.string().expect("source schema");
        let ceiling_offset = reader.offset;

        for offset in [role_offset, ceiling_offset] {
            let mut unknown_tag = bytes.clone();
            unknown_tag[offset] = u8::MAX;
            let error = OwnerMatrix::from_canonical_bytes(&unknown_tag)
                .expect_err("an unknown owner-matrix tag must refuse");
            assert_eq!(error.code(), "EulerOwnerMatrixMalformedCanonical");
        }

        let mut schema_fork = contract.owner_matrix().clone();
        schema_fork
            .rows
            .values_mut()
            .next()
            .expect("owner row")
            .source_schema
            .push_str("-changed");
        let fork_bytes =
            encode_owner_matrix_components(EULER_OWNER_MATRIX_SCHEMA_VERSION, &schema_fork.rows)
                .expect("bounded fork bytes");
        let error = OwnerMatrix::from_canonical_bytes(&fork_bytes)
            .expect_err("a complete schema fork must refuse semantically");
        assert_eq!(error.code(), "EulerContractGenericSchemaFork");
    }

    #[test]
    fn graph_decoder_refuses_maximum_plus_one_before_allocation() {
        let mut bytes = Vec::from(GRAPH_MAGIC.as_slice());
        bytes.extend_from_slice(&EULER_CLAIM_POLICY_SCHEMA_VERSION.to_le_bytes());
        bytes.extend_from_slice(&u32::try_from(MAX_EULER_CLAIMS + 1).unwrap().to_le_bytes());
        let error = EulerClaimGraph::from_canonical_bytes(&bytes)
            .expect_err("maximum-plus-one claim count must refuse");
        assert_eq!(error.code(), "EulerContractCardinality");

        let graph = build_frozen_contract()
            .expect("frozen contract")
            .claim_graph()
            .clone();
        let mut hostile = graph.clone();
        let dependency = hostile
            .dependencies
            .first()
            .expect("frozen graph dependency")
            .clone();
        hostile.dependencies = vec![dependency; MAX_EULER_CLAIMS * (MAX_EULER_CLAIMS - 1) + 1];
        let bytes = hostile.canonical_bytes().expect("hostile bounded graph");
        let error = EulerClaimGraph::from_canonical_bytes(&bytes)
            .expect_err("constructor maximum-plus-one dependency count must refuse in decoder");
        assert_eq!(error.code(), "EulerContractCardinality");
    }

    #[test]
    fn graph_decoder_refuses_wrong_version_and_trailing_data() {
        let mut wrong_version = Vec::from(GRAPH_MAGIC.as_slice());
        wrong_version.extend_from_slice(&2_u32.to_le_bytes());
        wrong_version.extend_from_slice(&1_u32.to_le_bytes());
        let error = EulerClaimGraph::from_canonical_bytes(&wrong_version)
            .expect_err("future graph version must refuse");
        assert_eq!(error.code(), "EulerContractUnsupportedVersion");
    }
}
