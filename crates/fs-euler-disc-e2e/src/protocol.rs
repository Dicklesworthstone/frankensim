//! Frozen Euler-disc contract constructor, structural checks, and bounded
//! candidate-assessment policy.
//!
//! A successful result from this module is only structurally eligible for a
//! later evidence/governance review. It is not a physical-validation,
//! mechanism, maturity, runtime, or release grant.

use core::{fmt, mem::size_of};
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::ContentHash;
use fs_evidence::vv::{
    AcceptanceCriterion, ApplicabilityDomain, ApplicabilityPoint, ApplicabilityPolicy,
    ArtifactHeader, ArtifactId, ArtifactKind, AxisId, CategoricalDomainAxis, ContextOfUse,
    DeclaredBudget, NumericDomainAxis, QoiId, QoiSpec, SeedDeclaration, UnitId, VV_ARTIFACT_FAMILY,
    VV_SCHEMA_VERSION,
};
use fs_evidence::{Color, ValidityDomain, validate_color_payload};
use fs_govern::evidence_contract::{AUTHORITY_ALGEBRA_VERSION, NoClaimBoundary};
use fs_ir::campaign::{
    CampaignClaim, CampaignClaimId, ClaimDependency, EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1,
    EvidenceGap, EvidenceGapId, EvidenceUse,
};

pub use crate::contract::{CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN, CONTRACT_CHECK_RECEIPT_DOMAIN};
use crate::contract::{
    CORE_NO_CLAIMS, ContractError, ContractIdentity, EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA,
    EULER_ASSESSMENT_IDENTITY_DOMAIN, EULER_CLAIM_REGISTRY, EULER_CONTRACT_SCHEMA_VERSION,
    EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, EULER_EVIDENCE_REQUIREMENT_REGISTRY,
    EULER_OWNER_ROLE_REGISTRY, EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN, EulerClaimGraph,
    EulerClaimKind, EulerClaimSpec, EulerContextExtension, EulerScientificContract,
    EvidenceRequirement, HypothesisSource, MAX_EULER_CLAIMS, MAX_EULER_TEXT_BYTES, OwnerMatrix,
    OwnerRole, OwnerRow, ScientificRisk,
};

/// Version of the local policy evaluator and its exact identity preimages.
pub const EULER_PROTOCOL_SCHEMA_VERSION: u32 = 1;
/// Maximum evidence rows accepted by one packet.
pub const MAX_EVIDENCE_RECORDS: usize = 32;
/// Maximum direct-DAG prerequisite receipts accepted by one assessment.
pub const MAX_PREREQUISITE_RECEIPTS: usize = MAX_EULER_CLAIMS;
/// Maximum bytes in one protocol-local canonical machine identifier.
pub const MAX_PROTOCOL_ID_BYTES: usize = 256;
/// Maximum axes admitted in one declaration-only validated-color regime.
pub const MAX_VALIDITY_DOMAIN_AXES: usize = 32;
/// Maximum canonical bytes for the validity-domain portion of one validated
/// color: the axis-count frame plus every axis/bounds row.
pub const MAX_VALIDITY_DOMAIN_CANONICAL_BYTES: usize = 8 * 1024;
/// Maximum canonical bytes in one evidence-reference packet.
pub const MAX_EVIDENCE_PACKET_BYTES: usize = 1024 * 1024;
/// Maximum bytes in one redacted claim-policy assessment JSON-lines record.
pub const MAX_ASSESSMENT_LOG_BYTES: usize = 32 * 1024;
/// Maximum bytes in one exact contract-check receipt transport.
pub const MAX_CONTRACT_CHECK_RECEIPT_BYTES: usize = 32 * 1024;

/// Explicit compatibility policy for the local protocol identity-preimage and
/// receipt schemas. V1 has no predecessor and never negotiates an unknown
/// version down to the current one.
pub fn protocol_migration_policy(schema_version: u32) -> Result<(), ContractError> {
    if schema_version == EULER_PROTOCOL_SCHEMA_VERSION {
        Ok(())
    } else {
        Err(ContractError::new(
            "EulerProtocolUnsupportedVersion",
            format!(
                "protocol schema {schema_version} is unsupported; v1 has no predecessor migration"
            ),
        ))
    }
}

// The byte preflight below mirrors the v2 shared Color codec. An upstream
// codec move must stop compilation until the wire-size calculation and the
// packet identity declaration are reviewed together.
const _: () = assert!(fs_evidence::COLOR_ALGEBRA_VERSION == 2);
// Keep the local hostile-input preflight locked to the shared Color identity
// grammar. An upstream limit change must stop compilation until the protocol
// schema and its bounded-error contract are reviewed together.
const MAX_PROTOCOL_COLOR_IDENTITY_BYTES: usize = 256;
const _: () = assert!(fs_evidence::MAX_COLOR_IDENTITY_BYTES == MAX_PROTOCOL_COLOR_IDENTITY_BYTES);

/// Owner-local identity declaration for complete evidence-reference packets.
pub const CLAIM_EVIDENCE_PACKET_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-euler-disc-e2e:claim-evidence-packet",
    "version_const=EULER_PROTOCOL_SCHEMA_VERSION",
    "version=1",
    "domain=org.frankensim.fs-euler-disc-e2e.claim-evidence-packet.v1",
    "domain_const=crates/fs-euler-disc-e2e/src/contract.rs#EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN",
    "encoder=claim_evidence_packet_identity",
    "encoder_helpers=packet_canonical_bytes,packet_canonical_len,evidence_authority_canonical_len,append_applicability_point,EvidenceAuthorityDeclaration::canonical_bytes,preflight_validity_domain_canonical_bytes",
    "schema_functions=ClaimEvidencePacket::try_new,ClaimEvidencePacket::canonical_bytes,ClaimEvidencePacket::verify_identity,EvidenceRecord::try_new,ProtocolSeed::not_applicable,ProtocolSeed::validate,ProtocolBudget::try_new,packet_canonical_bytes,packet_canonical_len,evidence_authority_canonical_len,packet_too_large,checked_size_add,checked_packet_len_add,checked_packet_text_len,claim_evidence_packet_identity,append_applicability_point,EvidenceAuthorityDeclaration::try_new,EvidenceAuthorityDeclaration::class,EvidenceAuthorityDeclaration::canonical_bytes,EvidenceRequirement::source_schema,EvidenceRequirement::source_kind,EvidenceRequirement::authority_class,EvidenceAuthorityClass::code,DeclaredEvidenceAccessClass::code,ReportedScientificDisposition::code,AssessmentDisposition::code,preflight_validity_domain_canonical_bytes,preflight_color_identity_bytes,nonzero_hash,checked_protocol_id,validate_protocol_id,crates/fs-euler-disc-e2e/src/contract.rs#ContractIdentity::as_hash,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimKind::id,crates/fs-euler-disc-e2e/src/contract.rs#EvidenceRequirement::code,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes,crates/fs-evidence/src/color.rs#Color::canonical_bytes,crates/fs-evidence/src/color.rs#validate_color_payload,crates/fs-evidence/src/lib.rs#ValidityDomain::bounds,crates/fs-evidence/src/vv/model.rs#ApplicabilityPoint::numeric,crates/fs-evidence/src/vv/model.rs#ApplicabilityPoint::categorical,crates/fs-evidence/src/vv/model.rs#ArtifactKind::canonical_wire_tag",
    "schema_constants=EULER_PROTOCOL_SCHEMA_VERSION,crates/fs-euler-disc-e2e/src/contract.rs#EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN,PACKET_MAGIC,MAX_PROTOCOL_ID_BYTES,MAX_PROTOCOL_COLOR_IDENTITY_BYTES,MAX_EVIDENCE_RECORDS,MAX_VALIDITY_DOMAIN_AXES,MAX_VALIDITY_DOMAIN_CANONICAL_BYTES,MAX_EVIDENCE_PACKET_BYTES,crates/fs-evidence/src/color.rs#COLOR_ALGEBRA_VERSION,crates/fs-evidence/src/color.rs#MAX_COLOR_IDENTITY_BYTES",
    "schema_dependencies=fs-euler-disc-e2e:scientific-contract",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=ClaimEvidencePacket,EvidenceRecord,ProtocolBudget,ProtocolSeed,EvidenceAuthorityDeclaration,DeclaredEvidenceAccessClass,ReportedScientificDisposition,AssessmentDisposition",
    "source_fields=ClaimEvidencePacket.schema_version:semantic,ClaimEvidencePacket.contract_identity:semantic,ClaimEvidencePacket.case_id:semantic,ClaimEvidencePacket.design_set_identity:semantic,ClaimEvidencePacket.aggregate_qoi_derivation_receipt_identity:semantic,ClaimEvidencePacket.claim:semantic,ClaimEvidencePacket.point:semantic,ClaimEvidencePacket.records:semantic,ClaimEvidencePacket.no_claims_accepted:semantic,ClaimEvidencePacket.target_fitted:semantic,ClaimEvidencePacket.reported_scientific_disposition:semantic,ClaimEvidencePacket.expected_disposition:semantic,ClaimEvidencePacket.units:semantic,ClaimEvidencePacket.seed:semantic,ClaimEvidencePacket.budget:semantic,ClaimEvidencePacket.identity:derived:cached-root-of-complete-packet,EvidenceRecord.contract_identity:derived:transitively-bound-by-evidence-registry,EvidenceRecord.claim:derived:transitively-bound-by-evidence-registry,EvidenceRecord.requirement:derived:transitively-bound-by-evidence-registry,EvidenceRecord.qois:derived:transitively-bound-by-evidence-registry,EvidenceRecord.authority:derived:transitively-bound-by-evidence-registry,EvidenceRecord.artifact_hash:derived:transitively-bound-by-evidence-registry,EvidenceRecord.source_id:derived:transitively-bound-by-evidence-registry,EvidenceRecord.source_schema:derived:transitively-bound-by-evidence-registry,EvidenceRecord.source_kind:derived:transitively-bound-by-evidence-registry,EvidenceRecord.schema_admission_receipt_hash:derived:transitively-bound-by-evidence-registry,EvidenceRecord.access_class:derived:transitively-bound-by-evidence-registry,EvidenceRecord.independent:derived:transitively-bound-by-evidence-registry,ProtocolBudget.max_wall_time_ms:derived:transitively-bound-by-protocol-budget,ProtocolBudget.max_memory_bytes:derived:transitively-bound-by-protocol-budget,ProtocolBudget.normalized_accuracy_limit:derived:transitively-bound-by-protocol-budget,ProtocolSeed.variant:derived:transitively-bound-by-seed-declaration,ProtocolSeed.value:derived:transitively-bound-by-seed-declaration,ProtocolSeed.reason:derived:transitively-bound-by-seed-declaration,EvidenceAuthorityDeclaration.variant:derived:transitively-bound-by-evidence-registry,EvidenceAuthorityDeclaration.receipt_hash:derived:transitively-bound-by-evidence-registry,EvidenceAuthorityDeclaration.color:derived:transitively-bound-by-evidence-registry,DeclaredEvidenceAccessClass.variant:derived:transitively-bound-by-evidence-registry,ReportedScientificDisposition.variant:derived:transitively-bound-by-reported-scientific-disposition,AssessmentDisposition.variant:derived:transitively-bound-by-expected-disposition",
    "source_bindings=ClaimEvidencePacket.schema_version>protocol-schema-version,ClaimEvidencePacket.contract_identity>contract-identity,ClaimEvidencePacket.case_id>case-id,ClaimEvidencePacket.design_set_identity>design-set-identity,ClaimEvidencePacket.aggregate_qoi_derivation_receipt_identity>aggregate-qoi-derivation-receipt-identity,ClaimEvidencePacket.claim>claim-kind,ClaimEvidencePacket.point>applicability-point-anchor,ClaimEvidencePacket.records>evidence-registry,ClaimEvidencePacket.no_claims_accepted>no-claim-acceptance,ClaimEvidencePacket.target_fitted>target-fitting-state,ClaimEvidencePacket.reported_scientific_disposition>reported-scientific-disposition,ClaimEvidencePacket.expected_disposition>expected-disposition,ClaimEvidencePacket.units>unit-set,ClaimEvidencePacket.seed>seed-declaration,ClaimEvidencePacket.budget>protocol-budget",
    "external_semantic_fields=identity-domain,identity-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian,color-canonical-codec,artifact-kind-wire-tags",
    "semantic_fields=identity-domain,identity-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian,color-canonical-codec,artifact-kind-wire-tags,protocol-schema-version,contract-identity,case-id,design-set-identity,aggregate-qoi-derivation-receipt-identity,claim-kind,applicability-point-anchor,no-claim-acceptance,target-fitting-state,reported-scientific-disposition,expected-disposition,unit-set,seed-declaration,evidence-registry,protocol-budget",
    "excluded_fields=applicability-point-signed-zero-bit:ieee-754-signed-zero-is-canonicalized-to-positive-zero,protocol-budget-normalized-accuracy-signed-zero-bit:ieee-754-signed-zero-is-canonicalized-to-positive-zero-before-storage",
    "consumers=ClaimEvidencePacket::try_new,ClaimEvidencePacket::canonical_bytes,ClaimEvidencePacket::identity,ClaimEvidencePacket::verify_identity,ClaimEvidencePacket::assess,PrerequisiteAssessmentReceipt::new,claim_policy_assessment_log",
    "mutations=identity-domain:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,identity-version:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,transport-magic:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,canonical-field-order:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,length-framing:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,fixed-numeric-little-endian:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,color-canonical-codec:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,artifact-kind-wire-tags:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,protocol-schema-version:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,contract-identity:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,case-id:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,design-set-identity:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,aggregate-qoi-derivation-receipt-identity:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,claim-kind:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,applicability-point-anchor:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,no-claim-acceptance:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,target-fitting-state:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,reported-scientific-disposition:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,expected-disposition:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,unit-set:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,seed-declaration:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,evidence-registry:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently,protocol-budget:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_identity_semantic_fields_move_independently",
    "nonsemantic_mutations=applicability-point-signed-zero-bit:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_signed_zero_inputs_are_nonsemantic,protocol-budget-normalized-accuracy-signed-zero-bit:crates/fs-euler-disc-e2e/src/protocol.rs#claim_evidence_packet_signed_zero_inputs_are_nonsemantic",
    "field_guard=classify_claim_evidence_packet_identity_fields",
    "transport_guard=ClaimEvidencePacket::verify_identity",
    "version_guard=crates/fs-euler-disc-e2e/tests/scientific_contract.rs#euler_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-euler-disc-e2e:claim-evidence-packet",
];

/// Owner-local identity declaration for direct claim-DAG receipts.
pub const PREREQUISITE_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-euler-disc-e2e:prerequisite-assessment-receipt",
    "version_const=EULER_PROTOCOL_SCHEMA_VERSION",
    "version=1",
    "domain=org.frankensim.fs-euler-disc-e2e.prerequisite-assessment-receipt.v1",
    "domain_const=crates/fs-euler-disc-e2e/src/contract.rs#EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN",
    "encoder=prerequisite_assessment_receipt_identity",
    "encoder_helpers=prerequisite_receipt_bytes,applicability_point_bytes",
    "schema_functions=PrerequisiteAssessmentReceipt::new,PrerequisiteAssessmentReceipt::canonical_bytes,PrerequisiteAssessmentReceipt::verify,ClaimPolicyAssessment::verify_identity,prerequisite_receipt_bytes,prerequisite_assessment_receipt_identity,applicability_point_bytes,append_applicability_point,crates/fs-euler-disc-e2e/src/contract.rs#ContractIdentity::as_hash,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimKind::id,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes,crates/fs-evidence/src/vv/model.rs#ApplicabilityPoint::numeric,crates/fs-evidence/src/vv/model.rs#ApplicabilityPoint::categorical",
    "schema_constants=EULER_PROTOCOL_SCHEMA_VERSION,crates/fs-euler-disc-e2e/src/contract.rs#EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,PREREQUISITE_RECEIPT_MAGIC",
    "schema_dependencies=fs-euler-disc-e2e:claim-evidence-packet,fs-euler-disc-e2e:claim-policy-assessment",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=PrerequisiteAssessmentReceipt",
    "source_fields=PrerequisiteAssessmentReceipt.schema_version:semantic,PrerequisiteAssessmentReceipt.contract_identity:semantic,PrerequisiteAssessmentReceipt.prerequisite:semantic,PrerequisiteAssessmentReceipt.dependent:semantic,PrerequisiteAssessmentReceipt.use_kind:semantic,PrerequisiteAssessmentReceipt.source_packet_identity:semantic,PrerequisiteAssessmentReceipt.source_assessment_identity:semantic,PrerequisiteAssessmentReceipt.source_design_set_identity:semantic,PrerequisiteAssessmentReceipt.source_point_bytes:semantic,PrerequisiteAssessmentReceipt.identity:derived:cached-root-of-complete-receipt",
    "source_bindings=PrerequisiteAssessmentReceipt.schema_version>protocol-schema-version,PrerequisiteAssessmentReceipt.contract_identity>contract-identity,PrerequisiteAssessmentReceipt.prerequisite>prerequisite-claim,PrerequisiteAssessmentReceipt.dependent>dependent-claim,PrerequisiteAssessmentReceipt.use_kind>evidence-use,PrerequisiteAssessmentReceipt.source_packet_identity>source-packet-identity,PrerequisiteAssessmentReceipt.source_assessment_identity>source-assessment-identity,PrerequisiteAssessmentReceipt.source_design_set_identity>source-design-set-identity,PrerequisiteAssessmentReceipt.source_point_bytes>source-applicability-point-anchor",
    "external_semantic_fields=identity-domain,identity-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian",
    "semantic_fields=identity-domain,identity-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian,protocol-schema-version,contract-identity,prerequisite-claim,dependent-claim,evidence-use,source-packet-identity,source-assessment-identity,source-design-set-identity,source-applicability-point-anchor",
    "excluded_fields=none",
    "consumers=PrerequisiteAssessmentReceipt::identity,PrerequisiteAssessmentReceipt::verify,ClaimEvidencePacket::assess,claim_policy_assessment_log",
    "mutations=identity-domain:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,identity-version:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,transport-magic:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,canonical-field-order:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,length-framing:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,fixed-numeric-little-endian:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,protocol-schema-version:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,contract-identity:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,prerequisite-claim:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,dependent-claim:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,evidence-use:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,source-packet-identity:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,source-assessment-identity:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,source-design-set-identity:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently,source-applicability-point-anchor:crates/fs-euler-disc-e2e/src/protocol.rs#prerequisite_receipt_identity_semantic_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_prerequisite_receipt_identity_fields",
    "transport_guard=PrerequisiteAssessmentReceipt::verify",
    "version_guard=crates/fs-euler-disc-e2e/tests/scientific_contract.rs#euler_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-euler-disc-e2e:prerequisite-assessment-receipt",
];

/// Owner-local identity declaration for retained policy assessments.
pub const CLAIM_POLICY_ASSESSMENT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-euler-disc-e2e:claim-policy-assessment",
    "version_const=EULER_PROTOCOL_SCHEMA_VERSION",
    "version=1",
    "domain=org.frankensim.fs-euler-disc-e2e.claim-policy-assessment.v1",
    "domain_const=crates/fs-euler-disc-e2e/src/contract.rs#EULER_ASSESSMENT_IDENTITY_DOMAIN",
    "encoder=claim_policy_assessment_identity",
    "encoder_helpers=assessment_canonical_bytes",
    "schema_functions=ClaimEvidencePacket::assess,ClaimEvidencePacket::verify_identity,StructurallyAdmittedEulerContract::receipt,StructurallyAdmittedEulerContract::contract,ContractCheckReceipt::verify_subject,assess_packet,prerequisite_violations,PrerequisiteAssessmentReceipt::verify,PrerequisiteAssessmentReceipt::canonical_bytes,applicability_point_bytes,append_applicability_point,evidence_use_code,claim_kind_for_id,point_violations,evidence_hash_references,access_class_violation,evidence_weakness,validated_regime_covers_point,EvidenceAuthorityDeclaration::class,EvidenceRequirement::source_schema,EvidenceRequirement::source_kind,EvidenceRequirement::authority_class,EvidenceAuthorityClass::code,DeclaredEvidenceAccessClass::code,AssessmentDisposition::code,ReportedScientificDisposition::code,ClaimPolicyAssessment::build,ClaimPolicyAssessment::verify_identity,ClaimPolicyAssessmentLog::verify_identity,assessment_canonical_bytes,claim_policy_assessment_identity,crates/fs-euler-disc-e2e/src/contract.rs#ContractError::code,crates/fs-euler-disc-e2e/src/contract.rs#ContractIdentity::as_hash,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::identity,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::claim_graph,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::context,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::extension,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimGraph::claim,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimGraph::dependencies,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimSpec::campaign,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimSpec::requirements,crates/fs-euler-disc-e2e/src/contract.rs#EulerContextExtension::hypothesis_sources,crates/fs-euler-disc-e2e/src/contract.rs#HypothesisSource::declaration_hash,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimKind::id,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimKind::forbids_target_fitting,crates/fs-euler-disc-e2e/src/contract.rs#EvidenceRequirement::code,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes,crates/fs-evidence/src/lib.rs#ValidityDomain::bound,crates/fs-evidence/src/lib.rs#ValidityDomain::bounds,crates/fs-evidence/src/vv/model.rs#ApplicabilityPoint::numeric,crates/fs-evidence/src/vv/model.rs#ApplicabilityPoint::categorical,crates/fs-evidence/src/vv/model.rs#ContextOfUse::applicability,crates/fs-evidence/src/vv/model.rs#ContextOfUse::qois,crates/fs-evidence/src/vv/model.rs#ApplicabilityDomain::numeric,crates/fs-evidence/src/vv/model.rs#ApplicabilityDomain::categorical,crates/fs-evidence/src/vv/model.rs#NumericDomainAxis::bounds,crates/fs-evidence/src/vv/model.rs#CategoricalDomainAxis::allowed,crates/fs-evidence/src/vv/model.rs#QoiSpec::unit",
    "schema_constants=EULER_PROTOCOL_SCHEMA_VERSION,MAX_PREREQUISITE_RECEIPTS,crates/fs-euler-disc-e2e/src/contract.rs#EULER_CLAIM_REGISTRY,crates/fs-euler-disc-e2e/src/contract.rs#EULER_ASSESSMENT_IDENTITY_DOMAIN,ASSESSMENT_MAGIC",
    "schema_dependencies=fs-euler-disc-e2e:claim-evidence-packet,fs-euler-disc-e2e:claim-policy-assessment-log",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=ClaimPolicyAssessment",
    "source_fields=ClaimPolicyAssessment.schema_version:semantic,ClaimPolicyAssessment.contract_identity:semantic,ClaimPolicyAssessment.packet_identity:semantic,ClaimPolicyAssessment.design_set_identity:semantic,ClaimPolicyAssessment.aggregate_qoi_derivation_receipt_identity:semantic,ClaimPolicyAssessment.point_bytes:semantic,ClaimPolicyAssessment.case_id:semantic,ClaimPolicyAssessment.claim:semantic,ClaimPolicyAssessment.disposition:semantic,ClaimPolicyAssessment.reported_scientific_disposition:semantic,ClaimPolicyAssessment.reasons:semantic,ClaimPolicyAssessment.log:semantic,ClaimPolicyAssessment.identity:derived:cached-root-of-complete-assessment",
    "source_bindings=ClaimPolicyAssessment.schema_version>protocol-schema-version,ClaimPolicyAssessment.contract_identity>contract-identity,ClaimPolicyAssessment.packet_identity>packet-identity,ClaimPolicyAssessment.design_set_identity>design-set-identity,ClaimPolicyAssessment.aggregate_qoi_derivation_receipt_identity>aggregate-qoi-derivation-receipt-identity,ClaimPolicyAssessment.point_bytes>applicability-point-anchor,ClaimPolicyAssessment.case_id>case-id,ClaimPolicyAssessment.claim>claim-kind,ClaimPolicyAssessment.disposition>assessment-disposition,ClaimPolicyAssessment.reported_scientific_disposition>reported-scientific-disposition,ClaimPolicyAssessment.reasons>reason-registry,ClaimPolicyAssessment.log>assessment-log-identity",
    "external_semantic_fields=identity-domain,identity-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian",
    "semantic_fields=identity-domain,identity-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian,protocol-schema-version,contract-identity,packet-identity,design-set-identity,aggregate-qoi-derivation-receipt-identity,applicability-point-anchor,case-id,claim-kind,assessment-disposition,reported-scientific-disposition,reason-registry,assessment-log-identity",
    "excluded_fields=none",
    "consumers=ClaimPolicyAssessment::identity,ClaimPolicyAssessment::verify_identity,ClaimPolicyAssessment::as_prerequisite_for,PrerequisiteAssessmentReceipt::new",
    "mutations=identity-domain:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,identity-version:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,transport-magic:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,canonical-field-order:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,length-framing:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,fixed-numeric-little-endian:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,protocol-schema-version:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,contract-identity:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,packet-identity:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,design-set-identity:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,aggregate-qoi-derivation-receipt-identity:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,applicability-point-anchor:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,case-id:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,claim-kind:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,assessment-disposition:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,reported-scientific-disposition:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,reason-registry:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently,assessment-log-identity:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_identity_semantic_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_claim_policy_assessment_identity_fields",
    "transport_guard=ClaimPolicyAssessment::verify_identity",
    "version_guard=crates/fs-euler-disc-e2e/tests/scientific_contract.rs#euler_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-euler-disc-e2e:claim-policy-assessment",
];

/// Owner-local identity declaration for exact checker receipts.
pub const CONTRACT_CHECK_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-euler-disc-e2e:contract-check-receipt",
    "version_const=EULER_PROTOCOL_SCHEMA_VERSION",
    "version=1",
    "domain=org.frankensim.fs-euler-disc-e2e.contract-check-receipt.v1",
    "domain_const=crates/fs-euler-disc-e2e/src/contract.rs#CONTRACT_CHECK_RECEIPT_DOMAIN",
    "encoder=contract_check_receipt_identity",
    "encoder_helpers=check_receipt_bytes",
    "schema_functions=ContractCheckReceipt::new,ContractCheckReceipt::canonical_bytes,ContractCheckReceipt::from_canonical_bytes,ContractCheckReceipt::verify_identity,ContractCheckReceipt::verify_subject,ProtocolReader::new,ProtocolReader::take,ProtocolReader::fixed,ProtocolReader::byte,ProtocolReader::u32,ProtocolReader::count,ProtocolReader::text,ProtocolReader::finish,protocol_migration_policy,check_frozen_contract,literal_frozen_hash,check_receipt_bytes,contract_check_receipt_identity,crates/fs-euler-disc-e2e/src/contract.rs#ContractError::code,crates/fs-euler-disc-e2e/src/contract.rs#ContractIdentity::from_hash,crates/fs-euler-disc-e2e/src/contract.rs#ContractIdentity::as_hash,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::canonical_bytes,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::from_canonical_bytes,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::context_canonical_bytes,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::context,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::context_hash,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::claim_graph,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::identity,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimGraph::canonical_bytes,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimGraph::from_canonical_bytes,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimGraph::content_hash,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes,crates/fs-blake3/src/lib.rs#ContentHash::from_hex,crates/fs-evidence/src/vv/codec.rs#canonical_artifact_bytes,crates/fs-evidence/src/vv/codec.rs#VvArtifact::from_canonical_bytes,crates/fs-evidence/src/vv/codec.rs#encode_context,crates/fs-evidence/src/vv/codec.rs#decode_context",
    "schema_constants=EULER_PROTOCOL_SCHEMA_VERSION,crates/fs-euler-disc-e2e/src/contract.rs#CONTRACT_CHECK_RECEIPT_DOMAIN,CONTRACT_CHECK_RECEIPT_MAGIC,CHECKER_ID,FROZEN_CONTEXT_HASH_HEX,FROZEN_CLAIM_GRAPH_HASH_HEX,FROZEN_CONTRACT_IDENTITY_HEX,MAX_CONTRACT_CHECK_RECEIPT_BYTES",
    "schema_dependencies=fs-euler-disc-e2e:claim-graph,fs-euler-disc-e2e:scientific-contract",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=ContractCheckReceipt",
    "source_fields=ContractCheckReceipt.schema_version:semantic,ContractCheckReceipt.checker_id:semantic,ContractCheckReceipt.subject:semantic,ContractCheckReceipt.context_hash:semantic,ContractCheckReceipt.graph_hash:semantic,ContractCheckReceipt.passed:semantic,ContractCheckReceipt.issues:semantic,ContractCheckReceipt.identity:derived:cached-root-of-complete-receipt",
    "source_bindings=ContractCheckReceipt.schema_version>protocol-schema-version,ContractCheckReceipt.checker_id>checker-id,ContractCheckReceipt.subject>subject-identity,ContractCheckReceipt.context_hash>context-hash,ContractCheckReceipt.graph_hash>graph-hash,ContractCheckReceipt.passed>pass-flag,ContractCheckReceipt.issues>issue-registry",
    "external_semantic_fields=identity-domain,identity-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian",
    "semantic_fields=identity-domain,identity-version,transport-magic,canonical-field-order,length-framing,fixed-numeric-little-endian,protocol-schema-version,checker-id,subject-identity,context-hash,graph-hash,pass-flag,issue-registry",
    "excluded_fields=none",
    "consumers=ContractCheckReceipt::identity,ContractCheckReceipt::verify_identity,ContractCheckReceipt::verify_subject,StructurallyAdmittedEulerContract::receipt,admit_frozen_contract",
    "mutations=identity-domain:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,identity-version:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,transport-magic:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,canonical-field-order:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,length-framing:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,fixed-numeric-little-endian:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,protocol-schema-version:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,checker-id:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,subject-identity:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,context-hash:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,graph-hash:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,pass-flag:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently,issue-registry:crates/fs-euler-disc-e2e/src/protocol.rs#contract_check_receipt_identity_semantic_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_contract_check_receipt_identity_fields",
    "transport_guard=ContractCheckReceipt::from_canonical_bytes",
    "version_guard=crates/fs-euler-disc-e2e/tests/scientific_contract.rs#euler_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-euler-disc-e2e:contract-check-receipt",
];

/// Owner-local identity declaration for the exact retained JSON line.
///
/// Prerequisite and aggregate-QoI derivation receipt identities are labeled
/// with their opaque routing addresses resolved through the leaf owner-matrix
/// registry. The log therefore depends on that registry, not recursively on
/// the routed schemas whose instances it records.
pub const CLAIM_POLICY_ASSESSMENT_LOG_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-euler-disc-e2e:claim-policy-assessment-log",
    "version_const=EULER_PROTOCOL_SCHEMA_VERSION",
    "version=1",
    "domain=org.frankensim.fs-euler-disc-e2e.claim-policy-assessment-log.v1",
    "domain_const=crates/fs-euler-disc-e2e/src/contract.rs#CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN",
    "encoder=claim_policy_assessment_log_identity",
    "encoder_helpers=claim_policy_assessment_log,json_string,push_string_array",
    "schema_functions=claim_policy_assessment_log_identity,claim_policy_assessment_log,ClaimPolicyAssessmentLog::from_json_line,ClaimPolicyAssessmentLog::verify_identity,validate_claim_policy_assessment_log_json_line,validate_assessment_log_evidence_sources,validate_assessment_log_relative_artifacts,expected_assessment_log_prerequisites,observed_assessment_log_authority_class,is_closed_assessment_log_reason,assessment_log_requirement_reason,assessment_log_weakness_requirement,assessment_log_requirement,assessment_log_claim_kind,assessment_log_evidence_slot,assessment_log_evidence_slot_code,assessment_log_access_class,expected_assessment_log_access_class,assessment_log_authority_class,assessment_log_evidence_use,expected_assessment_log_units,is_nonzero_lower_content_hash,malformed_assessment_log,AssessmentLogReader::new,AssessmentLogReader::expect,AssessmentLogReader::string,AssessmentLogReader::unsigned,AssessmentLogReader::boolean,AssessmentLogReader::string_array,AssessmentLogReader::finish,is_strictly_sorted,is_exact_lower_hex,validate_protocol_id,build_frozen_contract,claim_kind_for_id,evidence_use_code,EvidenceRequirement::source_schema,EvidenceRequirement::source_kind,EvidenceRequirement::authority_class,EvidenceAuthorityClass::code,AssessmentDisposition::code,ReportedScientificDisposition::code,protocol_migration_policy,evidence_hash_references,crates/fs-euler-disc-e2e/src/contract.rs#ContractIdentity::as_hash,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::identity,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::owner_matrix,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::claim_graph,crates/fs-euler-disc-e2e/src/contract.rs#EulerScientificContract::extension,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimGraph::dependencies,crates/fs-euler-disc-e2e/src/contract.rs#EulerContextExtension::hypothesis_sources,crates/fs-euler-disc-e2e/src/contract.rs#HypothesisSource::declaration_hash,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimKind::id,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimKind::forbids_target_fitting,crates/fs-euler-disc-e2e/src/contract.rs#EulerClaimKind::required_evidence,crates/fs-euler-disc-e2e/src/contract.rs#EvidenceRequirement::code,crates/fs-euler-disc-e2e/src/contract.rs#OwnerMatrix::rows,crates/fs-euler-disc-e2e/src/contract.rs#OwnerRow::source_schema,crates/fs-evidence/src/vv/model.rs#ArtifactKind::slug,json_string,push_string_array,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#ContentHash::from_hex,crates/fs-blake3/src/lib.rs#ContentHash::to_hex",
    "schema_constants=EULER_PROTOCOL_SCHEMA_VERSION,FROZEN_CONTRACT_IDENTITY_HEX,crates/fs-euler-disc-e2e/src/contract.rs#CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN,crates/fs-euler-disc-e2e/src/contract.rs#EULER_OWNER_MATRIX_IDENTITY_DOMAIN,crates/fs-euler-disc-e2e/src/contract.rs#EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN,crates/fs-euler-disc-e2e/src/contract.rs#EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,crates/fs-euler-disc-e2e/src/contract.rs#EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA,crates/fs-euler-disc-e2e/src/contract.rs#EULER_CLAIM_REGISTRY,crates/fs-euler-disc-e2e/src/contract.rs#EULER_EVIDENCE_REQUIREMENT_REGISTRY,PACKET_SCHEMA,MAX_ASSESSMENT_LOG_BYTES,MAX_ASSESSMENT_LOG_EVIDENCE_SOURCE_BYTES,MAX_ASSESSMENT_LOG_UNIT_ROWS,MAX_ASSESSMENT_LOG_UNIT_LIST_BYTES,MAX_ASSESSMENT_LOG_REASON_ROWS,MAX_ASSESSMENT_LOG_REASON_BYTES,MAX_ASSESSMENT_LOG_ARTIFACT_ROWS,MAX_PROTOCOL_ID_BYTES,MAX_EVIDENCE_RECORDS,MAX_PREREQUISITE_RECEIPTS,ASSESSMENT_LOG_ARTIFACT_SLOT,ASSESSMENT_LOG_SCHEMA_RECEIPT_SLOT,ASSESSMENT_LOG_ROLE_RECEIPT_SLOT,ASSESSMENT_LOG_BASE_EVIDENCE_SLOTS,crates/fs-euler-disc-e2e/src/contract.rs#MAX_EULER_CLAIMS,crates/fs-euler-disc-e2e/src/contract.rs#MAX_EULER_TEXT_BYTES,REPRODUCTION_COMMAND",
    "schema_dependencies=fs-euler-disc-e2e:claim-evidence-packet,fs-euler-disc-e2e:scientific-contract,fs-euler-disc-e2e:owner-matrix",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=ClaimPolicyAssessmentLog",
    "source_fields=ClaimPolicyAssessmentLog.json_line:semantic,ClaimPolicyAssessmentLog.identity:derived:cached-root-of-exact-json-line",
    "source_bindings=ClaimPolicyAssessmentLog.json_line>exact-json-line-bytes",
    "external_semantic_fields=identity-domain,identity-version",
    "semantic_fields=identity-domain,identity-version,exact-json-line-bytes",
    "excluded_fields=none",
    "consumers=ClaimPolicyAssessmentLog::from_json_line,ClaimPolicyAssessmentLog::identity,ClaimPolicyAssessmentLog::verify_identity,ClaimPolicyAssessment::log",
    "mutations=identity-domain:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_log_identity_semantic_fields_move_independently,identity-version:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_log_identity_semantic_fields_move_independently,exact-json-line-bytes:crates/fs-euler-disc-e2e/src/protocol.rs#claim_policy_assessment_log_identity_semantic_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_claim_policy_assessment_log_identity_fields",
    "transport_guard=ClaimPolicyAssessmentLog::from_json_line",
    "version_guard=crates/fs-euler-disc-e2e/tests/scientific_contract.rs#euler_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-euler-disc-e2e:claim-policy-assessment-log",
];

const CHECKER_ID: &str = "fs-euler-disc-e2e-structural-contract-checker-v1";
/// Literal-frozen hash of the exact generic Context encoded by contract v1.
///
/// This is a structural drift oracle, not evidence that the Context is
/// scientifically adequate.
pub const FROZEN_CONTEXT_HASH_HEX: &str =
    "c959a7f5920529fffee13d1d58e3e12bf99e3f546a452e45e0a1c4b281186942";
/// Literal-frozen hash of the exact nine-claim graph encoded by policy v1.
///
/// This is a structural drift oracle, not evidence that the graph's scientific
/// choices are correct or complete.
pub const FROZEN_CLAIM_GRAPH_HASH_HEX: &str =
    "53c810e7afa7c150abc2679126ad46a0c2b5ab048510f188938141b1ce0345fd";
/// Literal-frozen identity of the complete frozen scientific-contract v1.
///
/// Matching this identity establishes exact bytes under the declared hash
/// domain only. It does not establish source custody, physical validation, or
/// governance authority.
pub const FROZEN_CONTRACT_IDENTITY_HEX: &str =
    "e95ae98859836b49370bc0a75749f7c6687cd1552a73ae8177fcfafbcb3d5e60";
const PACKET_SCHEMA: &str = EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN;
const PACKET_MAGIC: &[u8; 8] = b"FSEDPK01";
const PREREQUISITE_RECEIPT_MAGIC: &[u8; 8] = b"FSEDPR01";
const ASSESSMENT_MAGIC: &[u8; 8] = b"FSEDAA01";
const CONTRACT_CHECK_RECEIPT_MAGIC: &[u8; 8] = b"FSEDCK01";
const REPRODUCTION_COMMAND: &str = "cargo test --locked -p fs-euler-disc-e2e --test scientific_contract -- g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded --exact --test-threads=1";
/// Exact focused checker-smoke command retained in every v1 protocol log.
///
/// This command reproduces the bounded log/receipt contract test. It is not a
/// replay command for the physical or numerical case named by a particular
/// log, because v1 retains logical artifact identities but no artifact
/// resolver or executable case bundle.
pub const ASSESSMENT_LOG_REPRODUCTION_COMMAND: &str = REPRODUCTION_COMMAND;

fn claim_evidence_packet_identity(bytes: &[u8]) -> ContentHash {
    fs_blake3::hash_domain(EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, bytes)
}

fn prerequisite_assessment_receipt_identity(bytes: &[u8]) -> ContentHash {
    fs_blake3::hash_domain(EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN, bytes)
}

fn claim_policy_assessment_identity(bytes: &[u8]) -> ContentHash {
    fs_blake3::hash_domain(EULER_ASSESSMENT_IDENTITY_DOMAIN, bytes)
}

fn contract_check_receipt_identity(bytes: &[u8]) -> ContentHash {
    fs_blake3::hash_domain(CONTRACT_CHECK_RECEIPT_DOMAIN, bytes)
}

fn claim_policy_assessment_log_identity(json_line: &str) -> ContentHash {
    fs_blake3::hash_domain(CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN, json_line.as_bytes())
}

fn foreign_error(code: &'static str, error: impl fmt::Display) -> ContractError {
    ContractError::new(code, error.to_string())
}

fn artifact_id(value: &str) -> Result<ArtifactId, ContractError> {
    ArtifactId::try_new(value).map_err(|error| foreign_error("EulerContractContextBuild", error))
}

fn qoi_id(value: &str) -> Result<QoiId, ContractError> {
    QoiId::try_new(value).map_err(|error| foreign_error("EulerContractContextBuild", error))
}

fn axis_id(value: &str) -> Result<AxisId, ContractError> {
    AxisId::try_new(value).map_err(|error| foreign_error("EulerContractContextBuild", error))
}

fn unit_id(value: &str) -> Result<UnitId, ContractError> {
    UnitId::try_new(value).map_err(|error| foreign_error("EulerContractContextBuild", error))
}

fn campaign_claim_id(kind: EulerClaimKind) -> Result<CampaignClaimId, ContractError> {
    CampaignClaimId::try_new(kind.id())
        .map_err(|error| foreign_error("EulerContractCampaignBuild", error))
}

fn build_header() -> Result<ArtifactHeader, ContractError> {
    let units = ["1", "g", "j", "k", "mm", "pa", "rad", "rad-per-s", "s"]
        .into_iter()
        .map(unit_id)
        .collect::<Result<Vec<_>, _>>()?;
    ArtifactHeader::try_new(
        artifact_id("euler-disc-context-v1")?,
        units,
        SeedDeclaration::NotApplicable {
            reason: "The immutable Context of Use performs no randomized computation.".to_owned(),
        },
        DeclaredBudget::NotApplicable {
            reason: "Accuracy budgets are declared per QoI rather than for context construction."
                .to_owned(),
        },
        DeclaredBudget::NotApplicable {
            reason: "Context construction contains no admitted long-running work.".to_owned(),
        },
        DeclaredBudget::NotApplicable {
            reason: "Context construction stores only bounded contract metadata.".to_owned(),
        },
        vec![
            (
                "fs-euler-disc-contract".to_owned(),
                EULER_CONTRACT_SCHEMA_VERSION.to_string(),
            ),
            ("fs-evidence-vv".to_owned(), VV_SCHEMA_VERSION.to_string()),
            (
                "fs-govern-authority".to_owned(),
                AUTHORITY_ALGEBRA_VERSION.to_string(),
            ),
            (
                "fs-ir-campaign".to_owned(),
                EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1.to_string(),
            ),
        ],
        vec![
            "declare-claim-graph".to_owned(),
            "declare-context-of-use".to_owned(),
            "declare-no-claim-boundary".to_owned(),
        ],
    )
    .map_err(|error| foreign_error("EulerContractContextBuild", error))
}

fn numeric_axis(
    id: &str,
    unit: &str,
    lo: f64,
    hi: f64,
) -> Result<NumericDomainAxis, ContractError> {
    NumericDomainAxis::try_new(axis_id(id)?, unit_id(unit)?, lo, hi)
        .map_err(|error| foreign_error("EulerContractContextBuild", error))
}

fn categorical_axis(id: &str, allowed: &[&str]) -> Result<CategoricalDomainAxis, ContractError> {
    CategoricalDomainAxis::try_new(
        axis_id(id)?,
        allowed.iter().map(|value| (*value).to_owned()).collect(),
    )
    .map_err(|error| foreign_error("EulerContractContextBuild", error))
}

fn build_applicability() -> Result<ApplicabilityDomain, ContractError> {
    ApplicabilityDomain::try_new(
        vec![
            numeric_axis("outer-radius", "mm", 12.0, 75.0)?,
            numeric_axis("body-thickness", "mm", 2.0, 25.0)?,
            numeric_axis("total-mass", "g", 20.0, 1_500.0)?,
            numeric_axis("edge-radius", "mm", 0.05, 5.0)?,
            numeric_axis("initial-angular-speed", "rad-per-s", 25.0, 250.0)?,
            numeric_axis("ambient-pressure", "pa", 80_000.0, 105_000.0)?,
            numeric_axis("ambient-temperature", "k", 275.0, 310.0)?,
            numeric_axis("base-slope", "rad", -0.01, 0.01)?,
        ],
        vec![
            categorical_axis(
                "disc-material-family",
                &["stainless-steel", "tungsten-alloy"],
            )?,
            categorical_axis(
                "mass-distribution-family",
                &["annular-ring", "center-weighted-cone", "uniform-disc"],
            )?,
            categorical_axis(
                "base-material-family",
                &["borosilicate-glass", "mirror-glass", "steel"],
            )?,
            categorical_axis("support-assembly", &["rigid-three-point-support"])?,
            categorical_axis("environment-regime", &["still-air-laboratory"])?,
            categorical_axis("contact-regime", &["single-rounded-edge-contact"])?,
            categorical_axis("observation-frame", &["laboratory-inertial-frame"])?,
        ],
    )
    .map_err(|error| foreign_error("EulerContractContextBuild", error))
}

fn qoi(
    id: &str,
    name: &str,
    unit: &str,
    acceptance: AcceptanceCriterion,
) -> Result<QoiSpec, ContractError> {
    QoiSpec::try_new(qoi_id(id)?, name, unit_id(unit)?, acceptance)
        .map_err(|error| foreign_error("EulerContractContextBuild", error))
}

fn build_context() -> Result<ContextOfUse, ContractError> {
    ContextOfUse::try_new(
        build_header()?,
        "Decide, claim by claim, whether Euler-disc simulations remain exploratory, are retained as calibrated reproductions, or are eligible candidates for separately governed blind physical review within the exact declared domain.",
        vec![
            qoi(
                "numerical-trajectory-error",
                "Absolute independent numerical trajectory discrepancy",
                "1",
                AcceptanceCriterion::ClosedRange {
                    lo: 0.0,
                    hi: 1.0e-8,
                },
            )?,
            qoi(
                "normalized-trajectory-discrepancy",
                "Normalized physical trajectory discrepancy",
                "1",
                AcceptanceCriterion::ClosedRange { lo: 0.0, hi: 0.05 },
            )?,
            qoi(
                "event-class-disposition",
                "Preregistered terminal-event or crossover class disposition",
                "1",
                AcceptanceCriterion::CategoryEquals {
                    expected: "matches-preregistered-event-class".to_owned(),
                },
            )?,
            qoi(
                "event-time-error",
                "Absolute terminal-event or crossover time error",
                "s",
                AcceptanceCriterion::ClosedRange { lo: 0.0, hi: 0.5 },
            )?,
            qoi(
                "qualitative-effect-disposition",
                "Preregistered qualitative-effect direction disposition",
                "1",
                AcceptanceCriterion::CategoryEquals {
                    expected: "matches-preregistered-direction".to_owned(),
                },
            )?,
            qoi(
                "configuration-ranking-disposition",
                "Preregistered configuration ranking disposition",
                "1",
                AcceptanceCriterion::CategoryEquals {
                    expected: "matches-preregistered-order-or-tie-rule".to_owned(),
                },
            )?,
            qoi(
                "optimum-containment-disposition",
                "Preregistered nonlinear optimum containment disposition",
                "1",
                AcceptanceCriterion::CategoryEquals {
                    expected: "contains-preregistered-optimum-under-exact-score".to_owned(),
                },
            )?,
            qoi(
                "optimum-interval-width",
                "Width of the retained nonlinear optimum interval",
                "mm",
                AcceptanceCriterion::ClosedRange { lo: 0.0, hi: 0.25 },
            )?,
            qoi(
                "energy-balance-residual",
                "Signed closed work-energy residual",
                "j",
                AcceptanceCriterion::ClosedRange {
                    lo: -0.001,
                    hi: 0.001,
                },
            )?,
            qoi(
                "energy-channel-fraction-error",
                "Absolute energy-channel allocation fraction error",
                "1",
                AcceptanceCriterion::ClosedRange { lo: 0.0, hi: 0.05 },
            )?,
            qoi(
                "rival-mechanism-disposition",
                "Preregistered rival-mechanism discrimination disposition",
                "1",
                AcceptanceCriterion::CategoryEquals {
                    expected: "discriminates-preregistered-rival".to_owned(),
                },
            )?,
        ],
        build_applicability()?,
        ApplicabilityPolicy::Refuse,
    )
    .map_err(|error| foreign_error("EulerContractContextBuild", error))
}

fn hypothesis(kind: EulerClaimKind) -> &'static str {
    match kind {
        EulerClaimKind::NumericalTrajectoryVerification => {
            "Independent discretization and implementation checks bound the numerical trajectory error without using physical target outcomes."
        }
        EulerClaimKind::CalibratedReproduction => {
            "A declared fitted model reproduces only its named calibration partition within the declared applicability domain."
        }
        EulerClaimKind::BlindTrajectoryPrediction => {
            "A frozen model predicts a protected physical trajectory without access to that trajectory during fitting or selection."
        }
        EulerClaimKind::EventOrCrossoverPrediction => {
            "A frozen model predicts a preregistered terminal event or crossover under the later exact event protocol."
        }
        EulerClaimKind::QualitativeEffectDirection => {
            "A frozen model predicts the preregistered direction of a controlled geometric, material, base, or environment contrast."
        }
        EulerClaimKind::Ranking => {
            "A frozen model predicts a preregistered ordering or tie among declared configurations."
        }
        EulerClaimKind::NonlinearOptimumInterval => {
            "A frozen model places a nonlinear design optimum inside a bounded preregistered interval."
        }
        EulerClaimKind::EnergyChannelAttribution => {
            "A frozen energy ledger closes and allocates work among uniquely owned, nonoverlapping channels."
        }
        EulerClaimKind::MechanismAttribution => {
            "Preregistered discriminating observables separate named rival mechanisms beyond exponent or stop-time agreement."
        }
    }
}

fn consequence(kind: EulerClaimKind) -> String {
    format!(
        "Retain the {} outcome as positive, negative, or inconclusive under its exact evidence ceiling; never infer a stronger claim.",
        kind.id()
    )
}

fn build_claim(kind: EulerClaimKind) -> Result<EulerClaimSpec, ContractError> {
    let qois = kind
        .required_qoi_ids()
        .iter()
        .map(|id| qoi_id(id))
        .collect::<Result<Vec<_>, _>>()?;
    let evidence_gaps = qois
        .iter()
        .enumerate()
        .map(|(index, qoi)| {
            let id = EvidenceGapId::try_new(format!("{}-gap-{}", kind.id(), index + 1))
                .map_err(|error| foreign_error("EulerContractCampaignBuild", error))?;
            Ok(EvidenceGap {
                id,
                qoi: qoi.clone(),
                expected_evidence: "euler-evidence-bundle-v1".to_owned(),
                description: format!(
                    "No target outcome is embedded here; later campaign artifacts must supply exact evidence for {}.",
                    qoi.as_str()
                ),
            })
        })
        .collect::<Result<Vec<_>, ContractError>>()?;
    EulerClaimSpec::try_new(
        kind,
        CampaignClaim {
            id: campaign_claim_id(kind)?,
            qois,
            hypothesis: hypothesis(kind).to_owned(),
            decision_consequence: consequence(kind),
            evidence_gaps,
        },
        kind.required_evidence().to_vec(),
    )
}

fn dependency(
    prerequisite: EulerClaimKind,
    dependent: EulerClaimKind,
    use_kind: EvidenceUse,
) -> Result<ClaimDependency, ContractError> {
    Ok(ClaimDependency {
        prerequisite: campaign_claim_id(prerequisite)?,
        dependent: campaign_claim_id(dependent)?,
        use_kind,
    })
}

fn build_graph() -> Result<EulerClaimGraph, ContractError> {
    let claims = EULER_CLAIM_REGISTRY
        .into_iter()
        .map(build_claim)
        .collect::<Result<Vec<_>, _>>()?;
    let numerical = EulerClaimKind::NumericalTrajectoryVerification;
    let calibrated = EulerClaimKind::CalibratedReproduction;
    let blind = EulerClaimKind::BlindTrajectoryPrediction;
    EulerClaimGraph::try_new(
        claims,
        vec![
            dependency(numerical, calibrated, EvidenceUse::ValidationInput)?,
            dependency(numerical, blind, EvidenceUse::ValidationInput)?,
            dependency(calibrated, blind, EvidenceUse::CalibrationInput)?,
            dependency(
                blind,
                EulerClaimKind::EventOrCrossoverPrediction,
                EvidenceUse::ValidationInput,
            )?,
            dependency(
                blind,
                EulerClaimKind::QualitativeEffectDirection,
                EvidenceUse::ValidationInput,
            )?,
            dependency(blind, EulerClaimKind::Ranking, EvidenceUse::ValidationInput)?,
            dependency(
                blind,
                EulerClaimKind::NonlinearOptimumInterval,
                EvidenceUse::ValidationInput,
            )?,
            dependency(
                blind,
                EulerClaimKind::EnergyChannelAttribution,
                EvidenceUse::ValidationInput,
            )?,
            dependency(
                blind,
                EulerClaimKind::MechanismAttribution,
                EvidenceUse::ValidationInput,
            )?,
            dependency(
                EulerClaimKind::EnergyChannelAttribution,
                EulerClaimKind::MechanismAttribution,
                EvidenceUse::ValidationInput,
            )?,
        ],
    )
}

fn source_declaration(id: &str, locator: &str) -> Result<HypothesisSource, ContractError> {
    HypothesisSource::try_new(id, locator)
}

fn build_extension() -> Result<EulerContextExtension, ContractError> {
    let exploratory = "retain-as-exploratory-or-calibrated-only";
    let candidate = "advance-to-separate-candidate-review";
    let refuse = "refuse-or-demote-the-claim";
    let retain_terminal = "retain-as-terminal-non-promotion";
    let physical_claims = EULER_CLAIM_REGISTRY
        .into_iter()
        .filter(|kind| *kind != EulerClaimKind::NumericalTrajectoryVerification)
        .collect::<Vec<_>>();
    let emergent_claims = EULER_CLAIM_REGISTRY
        .into_iter()
        .filter(|kind| kind.forbids_target_fitting())
        .collect::<Vec<_>>();
    EulerContextExtension::try_new(
        vec![
            "frankensim-developers".to_owned(),
            "independent-vv-reviewers".to_owned(),
            "research-software-users".to_owned(),
        ],
        "Precision squat cylindrical, center-weighted, and annular Euler-disc specimens plus declared base and three-point support populations; exact geometry and material packs remain separately content-bound.",
        "Indoor still-air laboratory conditions inside the numeric and categorical applicability axes; no vacuum, turbulent cross-flow, contaminated contact, or unmodeled moving support.",
        "laboratory-inertial-frame",
        vec![
            exploratory.to_owned(),
            candidate.to_owned(),
            refuse.to_owned(),
            retain_terminal.to_owned(),
        ],
        vec![
            ScientificRisk::try_new(
                "software-proof-laundering",
                "Treating deterministic code or solution verification as physical validation would overstate model authority.",
                5,
                physical_claims.clone(),
                refuse,
            )?,
            ScientificRisk::try_new(
                "protected-target-leakage",
                "Using protected target outcomes during fitting or model selection would invalidate an emergent-prediction claim.",
                5,
                emergent_claims,
                exploratory,
            )?,
            ScientificRisk::try_new(
                "mechanism-overidentification",
                "Inferring a unique loss mechanism from an exponent or stopping time would suppress viable rival explanations.",
                5,
                vec![EulerClaimKind::MechanismAttribution],
                refuse,
            )?,
            ScientificRisk::try_new(
                "silent-domain-extrapolation",
                "Reusing a successful case outside its apparatus, scale, material, support, frame, or environment could produce an unsupported decision.",
                4,
                physical_claims,
                refuse,
            )?,
            ScientificRisk::try_new(
                "negative-result-erasure",
                "Discarding negative or inconclusive outcomes would bias the retained campaign record.",
                4,
                EULER_CLAIM_REGISTRY.to_vec(),
                retain_terminal,
            )?,
        ],
        vec![source_declaration(
            "steve-mould-euler-disc-video-and-user-transcript",
            "https://www.youtube.com/watch?v=ti2qiU_JTUQ",
        )?],
    )
}

fn build_owner_matrix() -> Result<OwnerMatrix, ContractError> {
    OwnerMatrix::try_new(
        EULER_OWNER_ROLE_REGISTRY
            .into_iter()
            .map(|role| {
                OwnerRow::try_new(
                    role,
                    role.expected_owner_crate(),
                    role.expected_source_schema(),
                    role.expected_authority_ceiling(),
                )
            })
            .collect::<Result<Vec<_>, _>>()?,
    )
}

/// Build the one frozen v1 contract. This performs structural admission only.
pub fn build_frozen_contract() -> Result<EulerScientificContract, ContractError> {
    let no_claims = NoClaimBoundary::new(&CORE_NO_CLAIMS)
        .map_err(|error| foreign_error("EulerContractNoClaimBuild", error))?;
    EulerScientificContract::try_new(
        build_context()?,
        build_extension()?,
        build_graph()?,
        no_claims,
        build_owner_matrix()?,
    )
}

/// Caller-declared access class for one referenced evidence artifact.
///
/// This is policy metadata only. It is deliberately not the generic
/// `fs_evidence::vv::EvidencePartition`, does not prove observation membership,
/// and does not establish that a blind-release receipt was admitted. A later
/// campaign layer must resolve and re-admit the referenced generic artifact,
/// whole-case schema receipt, selection, and release before physical use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DeclaredEvidenceAccessClass {
    /// The evidence is structural and consumes no data partition.
    NotApplicable,
    /// Data may influence parameter fitting or model selection.
    Calibration,
    /// Held-apart physical validation data.
    Validation,
    /// Capability-separated protected holdout data.
    BlindHoldout,
}

impl DeclaredEvidenceAccessClass {
    // Stable machine code for retained logs. Plain comment: attached
    // rustdoc on dependency fragments is refused by the identity gate.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::NotApplicable => "not-applicable",
            Self::Calibration => "calibration",
            Self::Validation => "validation",
            Self::BlindHoldout => "blind-holdout",
        }
    }
}

/// Caller-reported scientific outcome retained independently from evidence
/// sufficiency. V1 does not evaluate an acceptance criterion or admit the
/// preregistered analysis artifact that would establish this result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ReportedScientificDisposition {
    /// The caller reports that the preregistered criterion was satisfied.
    Positive,
    /// The caller reports that the preregistered criterion was not satisfied.
    Negative,
    /// The caller reports that the analysis could not resolve the criterion.
    Inconclusive,
}

impl ReportedScientificDisposition {
    // Stable machine code for retained logs. Plain comment: attached
    // rustdoc on dependency fragments is refused by the identity gate.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Positive => "positive",
            Self::Negative => "negative",
            Self::Inconclusive => "inconclusive",
        }
    }
}

/// Local protocol result. None of these variants is a govern authority grant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AssessmentDisposition {
    /// Every required reference is present and the caller reports a positive
    /// outcome, but no referenced generic artifact or receipt has been
    /// resolved and re-admitted. A separate candidate review must do so before
    /// this result can contribute physical or governance authority.
    ReferenceCompleteCandidate,
    /// Evidence exists but is weaker than the requested local claim.
    DemotedCandidate,
    /// A caller-reported negative or inconclusive outcome is retained without
    /// promotion.
    RetainedTerminal,
    /// A hard structural, leakage, domain, or evidence guard failed.
    Refused,
}

impl AssessmentDisposition {
    // Stable machine code for retained logs. Plain comment: attached
    // rustdoc on dependency fragments is refused by the identity gate.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::ReferenceCompleteCandidate => "reference-complete-candidate-unreadmitted",
            Self::DemotedCandidate => "demoted-candidate",
            Self::RetainedTerminal => "retained-terminal-non-promotion",
            Self::Refused => "refused",
        }
    }
}

fn checked_protocol_id(
    field: &'static str,
    value: impl Into<String>,
) -> Result<String, ContractError> {
    let value = value.into();
    validate_protocol_id(field, &value)?;
    Ok(value)
}

fn validate_protocol_id(field: &'static str, value: &str) -> Result<(), ContractError> {
    let valid = !value.is_empty()
        && value.len() <= MAX_PROTOCOL_ID_BYTES
        && value.trim() == value
        && value.bytes().all(|byte| {
            byte.is_ascii_alphanumeric()
                || matches!(byte, b'-' | b'_' | b'.' | b'/' | b':' | b'@' | b'+' | b'=')
        });
    if !valid {
        return Err(ContractError::new(
            "EulerProtocolInvalidIdentity",
            format!("{field} must be a bounded canonical machine identity"),
        ));
    }
    Ok(())
}

fn nonzero_hash(field: &'static str, hash: ContentHash) -> Result<ContentHash, ContractError> {
    if hash.as_bytes().iter().all(|byte| *byte == 0) {
        return Err(ContractError::new(
            "EulerProtocolZeroIdentity",
            format!("{field} cannot be all zero"),
        ));
    }
    Ok(hash)
}

const fn evidence_use_code(use_kind: EvidenceUse) -> &'static str {
    match use_kind {
        EvidenceUse::CalibrationInput => "calibration-input",
        EvidenceUse::ValidationInput => "validation-input",
    }
}

fn append_applicability_point(
    bytes: &mut Vec<u8>,
    point: &ApplicabilityPoint,
) -> Result<(), ContractError> {
    write_len(bytes, point.numeric().len())?;
    for (axis, value) in point.numeric() {
        write_text(bytes, axis.as_str())?;
        let value = if *value == 0.0 { 0.0 } else { *value };
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    write_len(bytes, point.categorical().len())?;
    for (axis, value) in point.categorical() {
        write_text(bytes, axis.as_str())?;
        write_text(bytes, value)?;
    }
    Ok(())
}

fn applicability_point_bytes(point: &ApplicabilityPoint) -> Result<Vec<u8>, ContractError> {
    let mut bytes = Vec::with_capacity(512);
    append_applicability_point(&mut bytes, point)?;
    Ok(bytes)
}

impl EvidenceRequirement {
    /// Exact owner schema permitted to contain this evidence reference.
    #[must_use]
    pub const fn source_schema(self) -> &'static str {
        VV_ARTIFACT_FAMILY
    }

    /// Exact generic artifact kind that must carry this role in v1.
    ///
    /// A matching kind is only a structural container constraint. The
    /// referenced artifact still needs an exact generic schema-admission
    /// receipt, and neither fact establishes authenticity or adequacy.
    #[must_use]
    pub const fn source_kind(self) -> ArtifactKind {
        match self {
            EvidenceRequirement::SolutionVerification => ArtifactKind::SolutionVerificationReceipt,
            EvidenceRequirement::CalibrationPartition => ArtifactKind::CalibrationSplit,
            EvidenceRequirement::PreregisteredAnalysis
            | EvidenceRequirement::MultiplicityControl => ArtifactKind::ValidationPlan,
            EvidenceRequirement::CodeVerification
            | EvidenceRequirement::PhysicalValidation
            | EvidenceRequirement::BlindHoldout
            | EvidenceRequirement::ApplicabilityCheck
            | EvidenceRequirement::UncertaintyClosure
            | EvidenceRequirement::EnergyBalanceClosure
            | EvidenceRequirement::IndependentReconstruction
            | EvidenceRequirement::RivalMechanismDiscrimination => {
                ArtifactKind::PredictionAssessment
            }
        }
    }

    /// Categorical local declaration class expected for this role.
    #[must_use]
    pub const fn authority_class(self) -> EvidenceAuthorityClass {
        match self {
            Self::CodeVerification
            | Self::CalibrationPartition
            | Self::BlindHoldout
            | Self::PreregisteredAnalysis
            | Self::ApplicabilityCheck
            | Self::MultiplicityControl
            | Self::IndependentReconstruction => EvidenceAuthorityClass::StructuralProcess,
            Self::SolutionVerification | Self::UncertaintyClosure | Self::EnergyBalanceClosure => {
                EvidenceAuthorityClass::VerifiedNumerics
            }
            Self::PhysicalValidation | Self::RivalMechanismDiscrimination => {
                EvidenceAuthorityClass::ValidatedPhysical
            }
        }
    }
}

/// Closed evidence-declaration category. This is categorical rather than a
/// numeric authority score.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EvidenceAuthorityClass {
    /// A process, partition, preregistration, applicability, or independent-
    /// reconstruction receipt. No numerical color is implied.
    StructuralProcess,
    /// A finite Verified numerical enclosure is required.
    VerifiedNumerics,
    /// A Validated physical comparison with an in-domain regime is required.
    ValidatedPhysical,
}

impl EvidenceAuthorityClass {
    #[must_use]
    /// Stable machine code.
    pub const fn code(self) -> &'static str {
        match self {
            Self::StructuralProcess => "structural-process",
            Self::VerifiedNumerics => "verified-numerics",
            Self::ValidatedPhysical => "validated-physical",
        }
    }
}

/// One role-appropriate authority declaration. These values remain untrusted
/// references until their exact owner artifacts and receipts are re-admitted.
#[derive(Debug, Clone, PartialEq)]
pub enum EvidenceAuthorityDeclaration {
    /// Content identity of the role-specific structural/process receipt.
    StructuralProcess {
        /// Referenced receipt identity. This crate does not re-admit the
        /// referenced receipt or authenticate its producer.
        receipt_hash: ContentHash,
    },
    /// Generic Verified numerical payload.
    VerifiedNumerics {
        /// Declaration-only numerical color; this is not an `AdmittedColor`.
        color: Color,
    },
    /// Generic Validated physical payload.
    ValidatedPhysical {
        /// Declaration-only physical color; this is not an `AdmittedColor`.
        color: Color,
    },
}

/// Preflight the exact v2 canonical validity-domain payload without invoking
/// `Color::canonical_bytes` or allocating a second copy of caller-controlled
/// axis data.
///
/// The shared v2 Color codec writes one eight-byte axis count followed by, for
/// each axis, an eight-byte name-length frame, the name bytes, and two framed
/// eight-byte IEEE-754 bounds. Keep this calculation locked to
/// `COLOR_ALGEBRA_VERSION` above.
fn preflight_validity_domain_canonical_bytes(
    regime: &ValidityDomain,
) -> Result<usize, ContractError> {
    if regime.bounds().len() > MAX_VALIDITY_DOMAIN_AXES {
        return Err(ContractError::new(
            "EulerProtocolValidityDomainCardinality",
            format!(
                "validated-color regime has {} axes; the v1 limit is {MAX_VALIDITY_DOMAIN_AXES}",
                regime.bounds().len()
            ),
        ));
    }

    let mut bytes = size_of::<u64>();
    for axis in regime.bounds().keys() {
        if axis.len() > MAX_PROTOCOL_COLOR_IDENTITY_BYTES {
            return Err(ContractError::new(
                "EulerProtocolMalformedColor",
                format!(
                    "validated-color regime axis identity exceeds the v1 byte limit of {MAX_PROTOCOL_COLOR_IDENTITY_BYTES}"
                ),
            ));
        }
        let row_bytes = size_of::<u64>()
            .checked_add(axis.len())
            .and_then(|value| value.checked_add(2 * (size_of::<u64>() + size_of::<f64>())))
            .ok_or_else(|| {
                ContractError::new(
                    "EulerProtocolValidityDomainTooLarge",
                    "validated-color regime canonical length overflowed usize",
                )
            })?;
        bytes = bytes.checked_add(row_bytes).ok_or_else(|| {
            ContractError::new(
                "EulerProtocolValidityDomainTooLarge",
                "validated-color regime canonical length overflowed usize",
            )
        })?;
        if bytes > MAX_VALIDITY_DOMAIN_CANONICAL_BYTES {
            return Err(ContractError::new(
                "EulerProtocolValidityDomainTooLarge",
                format!(
                    "validated-color regime exceeds the v1 canonical-byte limit of {MAX_VALIDITY_DOMAIN_CANONICAL_BYTES}"
                ),
            ));
        }
    }
    Ok(bytes)
}

fn preflight_color_identity_bytes(color: &Color) -> Result<(), ContractError> {
    let oversized_field = match color {
        Color::Validated { dataset, .. } if dataset.len() > MAX_PROTOCOL_COLOR_IDENTITY_BYTES => {
            Some("validated-color dataset identity")
        }
        Color::Estimated { estimator, .. }
            if estimator.len() > MAX_PROTOCOL_COLOR_IDENTITY_BYTES =>
        {
            Some("estimated-color estimator identity")
        }
        Color::Verified { .. } | Color::Validated { .. } | Color::Estimated { .. } => None,
    };
    if let Some(field) = oversized_field {
        return Err(ContractError::new(
            "EulerProtocolMalformedColor",
            format!("{field} exceeds the v1 byte limit of {MAX_PROTOCOL_COLOR_IDENTITY_BYTES}"),
        ));
    }
    Ok(())
}

impl EvidenceAuthorityDeclaration {
    fn try_new(self) -> Result<Self, ContractError> {
        match &self {
            Self::StructuralProcess { receipt_hash } => {
                nonzero_hash("evidence.role_receipt_hash", *receipt_hash)?;
            }
            Self::VerifiedNumerics { color } | Self::ValidatedPhysical { color } => {
                if let Color::Validated { regime, .. } = color {
                    preflight_validity_domain_canonical_bytes(regime)?;
                }
                // `validate_color_payload` retains the offending identity in
                // its rich upstream error. Refuse oversized public strings
                // before that path can clone and format caller-sized input.
                preflight_color_identity_bytes(color)?;
                validate_color_payload(color)
                    .map_err(|error| foreign_error("EulerProtocolMalformedColor", error))?;
            }
        }
        Ok(self)
    }

    #[must_use]
    /// Categorical declaration class.
    pub const fn class(&self) -> EvidenceAuthorityClass {
        match self {
            Self::StructuralProcess { .. } => EvidenceAuthorityClass::StructuralProcess,
            Self::VerifiedNumerics { .. } => EvidenceAuthorityClass::VerifiedNumerics,
            Self::ValidatedPhysical { .. } => EvidenceAuthorityClass::ValidatedPhysical,
        }
    }

    #[must_use]
    /// Generic color when the category legitimately carries one.
    pub const fn color(&self) -> Option<&Color> {
        match self {
            Self::StructuralProcess { .. } => None,
            Self::VerifiedNumerics { color } | Self::ValidatedPhysical { color } => Some(color),
        }
    }

    fn canonical_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(128);
        match self {
            Self::StructuralProcess { receipt_hash } => {
                bytes.push(1);
                bytes.extend_from_slice(receipt_hash.as_bytes());
            }
            Self::VerifiedNumerics { color } => {
                bytes.push(2);
                bytes.extend_from_slice(&color.canonical_bytes());
            }
            Self::ValidatedPhysical { color } => {
                bytes.push(3);
                bytes.extend_from_slice(&color.canonical_bytes());
            }
        }
        bytes
    }
}

/// One content-bound evidence role. Construction validates payload shape but
/// does not establish authenticity or scientific adequacy.
#[derive(Debug, Clone, PartialEq)]
pub struct EvidenceRecord {
    contract_identity: ContractIdentity,
    claim: EulerClaimKind,
    requirement: EvidenceRequirement,
    qois: Vec<QoiId>,
    authority: EvidenceAuthorityDeclaration,
    artifact_hash: ContentHash,
    source_id: String,
    source_schema: String,
    source_kind: ArtifactKind,
    schema_admission_receipt_hash: ContentHash,
    access_class: DeclaredEvidenceAccessClass,
    independent: bool,
}

fn evidence_hash_references(record: &EvidenceRecord) -> [(&'static str, Option<ContentHash>); 3] {
    let role_receipt = match &record.authority {
        EvidenceAuthorityDeclaration::StructuralProcess { receipt_hash } => Some(*receipt_hash),
        EvidenceAuthorityDeclaration::VerifiedNumerics { .. }
        | EvidenceAuthorityDeclaration::ValidatedPhysical { .. } => None,
    };
    [
        ("artifact", Some(record.artifact_hash)),
        (
            "schema-admission-receipt",
            Some(record.schema_admission_receipt_hash),
        ),
        ("role-receipt", role_receipt),
    ]
}

impl EvidenceRecord {
    /// Construct one exact evidence-role binding.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        contract_identity: ContractIdentity,
        claim: EulerClaimKind,
        requirement: EvidenceRequirement,
        mut qois: Vec<QoiId>,
        authority: EvidenceAuthorityDeclaration,
        artifact_hash: ContentHash,
        source_id: impl Into<String>,
        source_schema: impl Into<String>,
        source_kind: ArtifactKind,
        schema_admission_receipt_hash: ContentHash,
        access_class: DeclaredEvidenceAccessClass,
        independent: bool,
    ) -> Result<Self, ContractError> {
        nonzero_hash("evidence.contract_identity", contract_identity.as_hash())?;
        let authority = authority.try_new()?;
        if qois.is_empty() || qois.len() > MAX_EULER_CLAIMS {
            return Err(ContractError::new(
                "EulerProtocolEvidenceQoiCardinality",
                "evidence must bind a bounded nonempty QoI set",
            ));
        }
        qois.sort();
        let original_qoi_count = qois.len();
        qois.dedup();
        if qois.len() != original_qoi_count {
            return Err(ContractError::new(
                "EulerProtocolDuplicateEvidenceQoi",
                "evidence contains a duplicate QoI binding",
            ));
        }
        let source_schema = checked_protocol_id("evidence.source_schema", source_schema)?;
        if source_schema != requirement.source_schema() {
            return Err(ContractError::new(
                "EulerProtocolSourceSchemaMismatch",
                format!(
                    "{} requires source schema {}, received {}",
                    requirement.code(),
                    requirement.source_schema(),
                    source_schema
                ),
            ));
        }
        if source_kind != requirement.source_kind() {
            return Err(ContractError::new(
                "EulerProtocolSourceKindMismatch",
                format!(
                    "{} requires generic artifact kind {}, received {}",
                    requirement.code(),
                    requirement.source_kind().slug(),
                    source_kind.slug()
                ),
            ));
        }
        Ok(Self {
            contract_identity,
            claim,
            requirement,
            qois,
            authority,
            artifact_hash: nonzero_hash("evidence.artifact_hash", artifact_hash)?,
            source_id: checked_protocol_id("evidence.source_id", source_id)?,
            source_schema,
            source_kind,
            schema_admission_receipt_hash: nonzero_hash(
                "evidence.schema_admission_receipt_hash",
                schema_admission_receipt_hash,
            )?,
            access_class,
            independent,
        })
    }

    #[must_use]
    /// Contract identity named by this evidence.
    pub const fn contract_identity(&self) -> ContractIdentity {
        self.contract_identity
    }

    #[must_use]
    /// Claim kind named by this evidence.
    pub const fn claim(&self) -> EulerClaimKind {
        self.claim
    }

    #[must_use]
    /// Single evidence role supplied by this artifact.
    pub const fn requirement(&self) -> EvidenceRequirement {
        self.requirement
    }

    #[must_use]
    /// Exact QoI set bound by the evidence.
    pub fn qois(&self) -> &[QoiId] {
        &self.qois
    }

    #[must_use]
    /// Role-appropriate categorical authority declaration.
    pub const fn authority(&self) -> &EvidenceAuthorityDeclaration {
        &self.authority
    }

    #[must_use]
    /// Nonzero artifact content identity.
    pub const fn artifact_hash(&self) -> ContentHash {
        self.artifact_hash
    }

    #[must_use]
    /// Stable source artifact identity.
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    #[must_use]
    /// Exact permitted source schema.
    pub fn source_schema(&self) -> &str {
        &self.source_schema
    }

    #[must_use]
    /// Exact generic V&V artifact kind containing this role.
    pub const fn source_kind(&self) -> ArtifactKind {
        self.source_kind
    }

    #[must_use]
    /// Content identity of the generic whole-case schema-admission receipt.
    /// The local protocol checks only that the reference is exact and nonzero;
    /// a later campaign must supply and re-verify the concrete receipt/case.
    pub const fn schema_admission_receipt_hash(&self) -> ContentHash {
        self.schema_admission_receipt_hash
    }

    #[must_use]
    /// Untrusted caller-declared calibration/validation/blind access class.
    ///
    /// This label never substitutes for the generic typed partition,
    /// observation selection, or blind-release receipt.
    pub const fn access_class(&self) -> DeclaredEvidenceAccessClass {
        self.access_class
    }

    #[must_use]
    /// Whether the evidence was reconstructed independently.
    pub const fn independent(&self) -> bool {
        self.independent
    }
}

/// Explicit seed state for a retained protocol case.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProtocolSeed {
    /// Explicit counter/randomization seed.
    Fixed {
        /// Exact counter/randomization seed value.
        value: u64,
    },
    /// Randomness is inapplicable for the named machine-readable reason.
    NotApplicable {
        /// Canonical reason code.
        reason: String,
    },
}

impl ProtocolSeed {
    /// Construct an explicit not-applicable seed declaration.
    pub fn not_applicable(reason: impl Into<String>) -> Result<Self, ContractError> {
        Ok(Self::NotApplicable {
            reason: checked_protocol_id("packet.seed.reason", reason)?,
        })
    }

    fn validate(&self) -> Result<(), ContractError> {
        if let Self::NotApplicable { reason } = self {
            validate_protocol_id("packet.seed.reason", reason)?;
        }
        Ok(())
    }
}

/// Explicit per-case resource and normalized computational-accuracy budgets.
///
/// The accuracy value is dimensionless by definition. It is a resource/control
/// declaration for the computation that prepares this packet, not a Context
/// QoI, claim threshold, observed discrepancy, or criterion-evaluation result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProtocolBudget {
    max_wall_time_ms: u64,
    max_memory_bytes: u64,
    normalized_accuracy_limit: f64,
}

impl ProtocolBudget {
    /// Construct bounded time, memory, and normalized-accuracy declarations.
    pub fn try_new(
        max_wall_time_ms: u64,
        max_memory_bytes: u64,
        normalized_accuracy_limit: f64,
    ) -> Result<Self, ContractError> {
        if max_wall_time_ms == 0
            || max_memory_bytes == 0
            || !normalized_accuracy_limit.is_finite()
            || normalized_accuracy_limit < 0.0
        {
            return Err(ContractError::new(
                "EulerProtocolInvalidBudget",
                "time/memory budgets must be positive and normalized dimensionless accuracy finite and nonnegative",
            ));
        }
        Ok(Self {
            max_wall_time_ms,
            max_memory_bytes,
            normalized_accuracy_limit: if normalized_accuracy_limit == 0.0 {
                0.0
            } else {
                normalized_accuracy_limit
            },
        })
    }

    #[must_use]
    /// Maximum wall time in milliseconds.
    pub const fn max_wall_time_ms(self) -> u64 {
        self.max_wall_time_ms
    }

    #[must_use]
    /// Maximum memory in bytes.
    pub const fn max_memory_bytes(self) -> u64 {
        self.max_memory_bytes
    }

    #[must_use]
    /// Nonnegative dimensionless computational-accuracy budget.
    ///
    /// This is never a scientific acceptance threshold.
    pub const fn normalized_accuracy_limit(self) -> f64 {
        self.normalized_accuracy_limit
    }
}

/// Complete candidate-assessment input with no raw observations.
#[derive(Debug, Clone, PartialEq)]
pub struct ClaimEvidencePacket {
    schema_version: u32,
    contract_identity: ContractIdentity,
    case_id: String,
    design_set_identity: ContentHash,
    aggregate_qoi_derivation_receipt_identity: ContentHash,
    claim: EulerClaimKind,
    point: ApplicabilityPoint,
    records: BTreeMap<EvidenceRequirement, EvidenceRecord>,
    no_claims_accepted: bool,
    target_fitted: bool,
    reported_scientific_disposition: ReportedScientificDisposition,
    expected_disposition: AssessmentDisposition,
    units: Vec<String>,
    seed: ProtocolSeed,
    budget: ProtocolBudget,
    identity: ContentHash,
}

fn packet_too_large() -> ContractError {
    ContractError::new(
        "EulerProtocolPacketTooLarge",
        "canonical evidence packet exceeds its byte budget",
    )
}

fn checked_size_add(total: &mut usize, additional: usize) -> Result<(), ContractError> {
    *total = total.checked_add(additional).ok_or_else(packet_too_large)?;
    Ok(())
}

fn checked_packet_len_add(total: &mut usize, additional: usize) -> Result<(), ContractError> {
    checked_size_add(total, additional)?;
    if *total > MAX_EVIDENCE_PACKET_BYTES {
        return Err(packet_too_large());
    }
    Ok(())
}

fn checked_packet_text_len(total: &mut usize, value: &str) -> Result<(), ContractError> {
    checked_packet_len_add(total, size_of::<u32>())?;
    checked_packet_len_add(total, value.len())
}

/// Exact shared Color-v2 bytes nested after the local declaration tag.
///
/// This deliberately mirrors `Color::canonical_bytes` without allocating its
/// temporary buffer. The compile-time Color-v2 assertion above turns an
/// upstream codec-version move into a mandatory review of this calculation.
fn evidence_authority_canonical_len(
    authority: &EvidenceAuthorityDeclaration,
) -> Result<usize, ContractError> {
    let mut authority_len = 1_usize;
    match authority {
        EvidenceAuthorityDeclaration::StructuralProcess { .. } => {
            checked_size_add(&mut authority_len, 32)?;
        }
        EvidenceAuthorityDeclaration::VerifiedNumerics { color }
        | EvidenceAuthorityDeclaration::ValidatedPhysical { color } => {
            let mut color_len = 2_usize; // Color-v2 version and variant tags.
            match color {
                Color::Verified { .. } => {
                    // Two u64-framed IEEE-754 payloads.
                    checked_size_add(&mut color_len, 2 * (size_of::<u64>() + size_of::<f64>()))?;
                }
                Color::Validated {
                    regime, dataset, ..
                } => {
                    checked_size_add(&mut color_len, size_of::<u64>())?;
                    checked_size_add(&mut color_len, dataset.len())?;
                    checked_size_add(
                        &mut color_len,
                        preflight_validity_domain_canonical_bytes(regime)?,
                    )?;
                }
                Color::Estimated { estimator, .. } => {
                    checked_size_add(&mut color_len, size_of::<u64>())?;
                    checked_size_add(&mut color_len, estimator.len())?;
                    checked_size_add(&mut color_len, size_of::<u64>() + size_of::<f64>())?;
                }
            }
            checked_size_add(&mut authority_len, color_len)?;
        }
    }
    Ok(authority_len)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn packet_canonical_len(
    case_id: &str,
    claim: EulerClaimKind,
    point: &ApplicabilityPoint,
    records: &BTreeMap<EvidenceRequirement, EvidenceRecord>,
    reported_scientific_disposition: ReportedScientificDisposition,
    expected_disposition: AssessmentDisposition,
    units: &[String],
    seed: &ProtocolSeed,
) -> Result<usize, ContractError> {
    let mut len = 0_usize;
    checked_packet_len_add(&mut len, PACKET_MAGIC.len())?;
    checked_packet_len_add(&mut len, size_of::<u32>())?;
    checked_packet_len_add(&mut len, 32)?;
    checked_packet_text_len(&mut len, case_id)?;
    // Exact design-set and aggregate-QoI derivation-receipt identities.
    checked_packet_len_add(&mut len, 2 * 32)?;
    checked_packet_text_len(&mut len, claim.id())?;

    checked_packet_len_add(&mut len, size_of::<u32>())?;
    for (axis, _) in point.numeric() {
        checked_packet_text_len(&mut len, axis.as_str())?;
        checked_packet_len_add(&mut len, size_of::<f64>())?;
    }
    checked_packet_len_add(&mut len, size_of::<u32>())?;
    for (axis, value) in point.categorical() {
        checked_packet_text_len(&mut len, axis.as_str())?;
        checked_packet_text_len(&mut len, value)?;
    }

    checked_packet_len_add(&mut len, size_of::<u32>())?;
    for record in records.values() {
        checked_packet_text_len(&mut len, record.requirement.code())?;
        checked_packet_len_add(&mut len, 32)?;
        checked_packet_text_len(&mut len, record.claim.id())?;
        checked_packet_len_add(&mut len, size_of::<u32>())?;
        for qoi in &record.qois {
            checked_packet_text_len(&mut len, qoi.as_str())?;
        }
        checked_packet_len_add(&mut len, size_of::<u32>())?;
        checked_packet_len_add(
            &mut len,
            evidence_authority_canonical_len(&record.authority)?,
        )?;
        checked_packet_len_add(&mut len, 32)?;
        checked_packet_text_len(&mut len, &record.source_id)?;
        checked_packet_text_len(&mut len, &record.source_schema)?;
        checked_packet_len_add(&mut len, 1)?;
        checked_packet_len_add(&mut len, 32)?;
        checked_packet_text_len(&mut len, record.access_class.code())?;
        checked_packet_len_add(&mut len, 1)?;
    }

    checked_packet_len_add(&mut len, 2)?;
    checked_packet_text_len(&mut len, reported_scientific_disposition.code())?;
    checked_packet_text_len(&mut len, expected_disposition.code())?;
    checked_packet_len_add(&mut len, size_of::<u32>())?;
    for unit in units {
        checked_packet_text_len(&mut len, unit)?;
    }
    checked_packet_len_add(&mut len, 1)?;
    match seed {
        ProtocolSeed::Fixed { .. } => checked_packet_len_add(&mut len, size_of::<u64>())?,
        ProtocolSeed::NotApplicable { reason } => checked_packet_text_len(&mut len, reason)?,
    }
    checked_packet_len_add(&mut len, 3 * size_of::<u64>())?;
    Ok(len)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn packet_canonical_bytes(
    schema_version: u32,
    contract_identity: ContractIdentity,
    case_id: &str,
    design_set_identity: ContentHash,
    aggregate_qoi_derivation_receipt_identity: ContentHash,
    claim: EulerClaimKind,
    point: &ApplicabilityPoint,
    records: &BTreeMap<EvidenceRequirement, EvidenceRecord>,
    no_claims_accepted: bool,
    target_fitted: bool,
    reported_scientific_disposition: ReportedScientificDisposition,
    expected_disposition: AssessmentDisposition,
    units: &[String],
    seed: &ProtocolSeed,
    budget: ProtocolBudget,
) -> Result<Vec<u8>, ContractError> {
    let canonical_len = packet_canonical_len(
        case_id,
        claim,
        point,
        records,
        reported_scientific_disposition,
        expected_disposition,
        units,
        seed,
    )?;
    // The exact checked preflight above happens before this packet buffer is
    // allocated. At most the admitted one-mebibyte budget is requested.
    let mut bytes = Vec::with_capacity(canonical_len);
    bytes.extend_from_slice(PACKET_MAGIC);
    bytes.extend_from_slice(&schema_version.to_le_bytes());
    bytes.extend_from_slice(contract_identity.as_hash().as_bytes());
    write_text(&mut bytes, case_id)?;
    bytes.extend_from_slice(design_set_identity.as_bytes());
    bytes.extend_from_slice(aggregate_qoi_derivation_receipt_identity.as_bytes());
    write_text(&mut bytes, claim.id())?;

    append_applicability_point(&mut bytes, point)?;

    write_len(&mut bytes, records.len())?;
    let mut canonical_records = records.values().collect::<Vec<_>>();
    canonical_records.sort_by_key(|record| record.requirement.code());
    for record in canonical_records {
        write_text(&mut bytes, record.requirement.code())?;
        bytes.extend_from_slice(record.contract_identity.as_hash().as_bytes());
        write_text(&mut bytes, record.claim.id())?;
        write_len(&mut bytes, record.qois.len())?;
        for qoi in &record.qois {
            write_text(&mut bytes, qoi.as_str())?;
        }
        let authority = record.authority.canonical_bytes();
        write_len(&mut bytes, authority.len())?;
        bytes.extend_from_slice(&authority);
        bytes.extend_from_slice(record.artifact_hash.as_bytes());
        write_text(&mut bytes, &record.source_id)?;
        write_text(&mut bytes, &record.source_schema)?;
        bytes.push(record.source_kind.canonical_wire_tag());
        bytes.extend_from_slice(record.schema_admission_receipt_hash.as_bytes());
        write_text(&mut bytes, record.access_class.code())?;
        bytes.push(u8::from(record.independent));
    }

    bytes.push(u8::from(no_claims_accepted));
    bytes.push(u8::from(target_fitted));
    write_text(&mut bytes, reported_scientific_disposition.code())?;
    write_text(&mut bytes, expected_disposition.code())?;
    write_len(&mut bytes, units.len())?;
    for unit in units {
        write_text(&mut bytes, unit)?;
    }
    match seed {
        ProtocolSeed::Fixed { value } => {
            bytes.push(1);
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        ProtocolSeed::NotApplicable { reason } => {
            bytes.push(2);
            write_text(&mut bytes, reason)?;
        }
    }
    bytes.extend_from_slice(&budget.max_wall_time_ms.to_le_bytes());
    bytes.extend_from_slice(&budget.max_memory_bytes.to_le_bytes());
    bytes.extend_from_slice(&budget.normalized_accuracy_limit.to_bits().to_le_bytes());
    debug_assert_eq!(
        bytes.len(),
        canonical_len,
        "packet length preflight must mirror the canonical writer exactly"
    );
    if bytes.len() > MAX_EVIDENCE_PACKET_BYTES {
        return Err(packet_too_large());
    }
    Ok(bytes)
}

impl ClaimEvidencePacket {
    /// Construct and canonicalize an evidence packet without raw observations.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        contract_identity: ContractIdentity,
        case_id: impl Into<String>,
        design_set_identity: ContentHash,
        aggregate_qoi_derivation_receipt_identity: ContentHash,
        claim: EulerClaimKind,
        point: ApplicabilityPoint,
        records: Vec<EvidenceRecord>,
        no_claims_accepted: bool,
        target_fitted: bool,
        reported_scientific_disposition: ReportedScientificDisposition,
        expected_disposition: AssessmentDisposition,
        mut units: Vec<String>,
        seed: ProtocolSeed,
        budget: ProtocolBudget,
    ) -> Result<Self, ContractError> {
        nonzero_hash("packet.contract_identity", contract_identity.as_hash())?;
        nonzero_hash("packet.design_set_identity", design_set_identity)?;
        nonzero_hash(
            "packet.aggregate_qoi_derivation_receipt_identity",
            aggregate_qoi_derivation_receipt_identity,
        )?;
        if design_set_identity == aggregate_qoi_derivation_receipt_identity {
            return Err(ContractError::new(
                "EulerProtocolCrossRoleEvidenceAlias",
                "design-set and aggregate-QoI derivation-receipt identities must be distinct",
            ));
        }
        if records.len() > MAX_EVIDENCE_RECORDS {
            return Err(ContractError::new(
                "EulerProtocolEvidenceCardinality",
                "evidence packet exceeds its bounded row count",
            ));
        }
        let mut by_requirement = BTreeMap::new();
        let mut logical_hash_roles =
            BTreeMap::<ContentHash, (EvidenceRequirement, &'static str)>::new();
        for record in records {
            if record.contract_identity != contract_identity || record.claim != claim {
                return Err(ContractError::new(
                    "EulerProtocolEvidenceBindingMismatch",
                    "evidence contract/claim binding differs from its packet",
                ));
            }
            for (slot, reference) in evidence_hash_references(&record) {
                let Some(reference) = reference else {
                    continue;
                };
                if reference == design_set_identity
                    || reference == aggregate_qoi_derivation_receipt_identity
                {
                    return Err(ContractError::new(
                        "EulerProtocolCrossRoleEvidenceAlias",
                        format!(
                            "evidence {slot} for {} cannot alias the design-set or aggregate-QoI derivation-receipt identity",
                            record.requirement.code()
                        ),
                    ));
                }
                if let Some((first_role, first_slot)) =
                    logical_hash_roles.insert(reference, (record.requirement, slot))
                {
                    return Err(ContractError::new(
                        "EulerProtocolCrossRoleEvidenceAlias",
                        format!(
                            "logical hash {reference} cannot occupy {first_slot} for {} and {slot} for {} without an explicit composite receipt",
                            first_role.code(),
                            record.requirement.code()
                        ),
                    ));
                }
            }
            let requirement = record.requirement;
            if by_requirement.insert(requirement, record).is_some() {
                return Err(ContractError::new(
                    "EulerProtocolDuplicateEvidenceRole",
                    format!("duplicate evidence role {}", requirement.code()),
                ));
            }
        }
        seed.validate()?;
        if units.is_empty() || units.len() > MAX_EULER_CLAIMS * 2 {
            return Err(ContractError::new(
                "EulerProtocolUnitCardinality",
                "packet units must be explicit, nonempty, and bounded",
            ));
        }
        units = units
            .into_iter()
            .map(|unit| checked_protocol_id("packet.unit", unit))
            .collect::<Result<Vec<_>, _>>()?;
        units.sort();
        let original_unit_count = units.len();
        units.dedup();
        if units.len() != original_unit_count {
            return Err(ContractError::new(
                "EulerProtocolDuplicateUnit",
                "packet units contain a duplicate",
            ));
        }
        let case_id = checked_protocol_id("packet.case_id", case_id)?;
        let canonical = packet_canonical_bytes(
            EULER_PROTOCOL_SCHEMA_VERSION,
            contract_identity,
            &case_id,
            design_set_identity,
            aggregate_qoi_derivation_receipt_identity,
            claim,
            &point,
            &by_requirement,
            no_claims_accepted,
            target_fitted,
            reported_scientific_disposition,
            expected_disposition,
            &units,
            &seed,
            budget,
        )?;
        let identity = claim_evidence_packet_identity(&canonical);
        Ok(Self {
            schema_version: EULER_PROTOCOL_SCHEMA_VERSION,
            contract_identity,
            case_id,
            design_set_identity,
            aggregate_qoi_derivation_receipt_identity,
            claim,
            point,
            records: by_requirement,
            no_claims_accepted,
            target_fitted,
            reported_scientific_disposition,
            expected_disposition,
            units,
            seed,
            budget,
            identity,
        })
    }

    #[must_use]
    /// Packet schema version.
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    /// Stable case identity.
    pub fn case_id(&self) -> &str {
        &self.case_id
    }

    #[must_use]
    /// Content identity of the complete comparison/configuration/search set.
    ///
    /// The applicability point is only a campaign anchor; relational claims are
    /// scoped by this exact design set.
    pub const fn design_set_identity(&self) -> ContentHash {
        self.design_set_identity
    }

    #[must_use]
    /// Content identity of the downstream Context-bound receipt that derives
    /// this claim's aggregate Context QoIs from detailed observables.
    ///
    /// The downstream receipt must bind the exact Context, detailed-observable
    /// registry, design set, claim, and aggregate-QoI scoring scope. Its future
    /// admission checker must cross-check that design set with this packet.
    /// This crate requires and binds the reference but cannot perform that
    /// admission or execute the downstream scoring artifact.
    pub const fn aggregate_qoi_derivation_receipt_identity(&self) -> ContentHash {
        self.aggregate_qoi_derivation_receipt_identity
    }

    #[must_use]
    /// Requested claim kind.
    pub const fn claim(&self) -> EulerClaimKind {
        self.claim
    }

    #[must_use]
    /// Exact frozen contract identity.
    pub const fn contract_identity(&self) -> ContractIdentity {
        self.contract_identity
    }

    #[must_use]
    /// Campaign-anchor point evaluated against the Context applicability domain.
    /// Relational scope is carried separately by `design_set_identity`.
    pub const fn point(&self) -> &ApplicabilityPoint {
        &self.point
    }

    #[must_use]
    /// Evidence rows keyed by their unique roles.
    pub const fn records(&self) -> &BTreeMap<EvidenceRequirement, EvidenceRecord> {
        &self.records
    }

    /// Exact canonical packet bytes, including every point, role, color,
    /// partition, receipt, seed, budget, and no-claim/target declaration.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ContractError> {
        packet_canonical_bytes(
            self.schema_version,
            self.contract_identity,
            &self.case_id,
            self.design_set_identity,
            self.aggregate_qoi_derivation_receipt_identity,
            self.claim,
            &self.point,
            &self.records,
            self.no_claims_accepted,
            self.target_fitted,
            self.reported_scientific_disposition,
            self.expected_disposition,
            &self.units,
            &self.seed,
            self.budget,
        )
    }

    #[must_use]
    /// Domain-separated identity of the complete canonical packet.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Recompute and verify the stored packet identity.
    pub fn verify_identity(&self) -> Result<(), ContractError> {
        let expected = claim_evidence_packet_identity(&self.canonical_bytes()?);
        if expected != self.identity {
            return Err(ContractError::new(
                "EulerProtocolPacketIdentityMismatch",
                "evidence packet identity does not match its semantic fields",
            ));
        }
        Ok(())
    }

    /// Evaluate the packet under the exact frozen contract.
    pub fn assess(
        &self,
        admitted: &StructurallyAdmittedEulerContract,
        prerequisites: &[PrerequisiteAssessmentReceipt],
    ) -> Result<ClaimPolicyAssessment, ContractError> {
        if prerequisites.len() > MAX_PREREQUISITE_RECEIPTS {
            return Err(ContractError::new(
                "EulerProtocolPrerequisiteCardinality",
                format!(
                    "assessment accepts at most {MAX_PREREQUISITE_RECEIPTS} direct-DAG prerequisite receipts"
                ),
            ));
        }
        admitted.receipt().verify_subject(admitted.contract())?;
        self.verify_identity()?;
        assess_packet(admitted.contract(), self, prerequisites)
    }
}

/// Content-bound proof that one structurally complete caller-reported-positive
/// assessment is being consumed along one exact direct claim-graph edge at the
/// same design set and campaign-anchor applicability point. It transfers no
/// physical authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrerequisiteAssessmentReceipt {
    schema_version: u32,
    contract_identity: ContractIdentity,
    prerequisite: EulerClaimKind,
    dependent: EulerClaimKind,
    use_kind: EvidenceUse,
    source_packet_identity: ContentHash,
    source_assessment_identity: ContentHash,
    source_design_set_identity: ContentHash,
    source_point_bytes: Vec<u8>,
    identity: ContentHash,
}

#[allow(clippy::too_many_arguments)]
fn prerequisite_receipt_bytes(
    schema_version: u32,
    contract_identity: ContractIdentity,
    prerequisite: EulerClaimKind,
    dependent: EulerClaimKind,
    use_kind: EvidenceUse,
    source_packet_identity: ContentHash,
    source_assessment_identity: ContentHash,
    source_design_set_identity: ContentHash,
    source_point_bytes: &[u8],
) -> Result<Vec<u8>, ContractError> {
    let mut bytes = Vec::with_capacity(256 + source_point_bytes.len());
    bytes.extend_from_slice(PREREQUISITE_RECEIPT_MAGIC);
    bytes.extend_from_slice(&schema_version.to_le_bytes());
    bytes.extend_from_slice(contract_identity.as_hash().as_bytes());
    write_text(&mut bytes, prerequisite.id())?;
    write_text(&mut bytes, dependent.id())?;
    write_text(&mut bytes, evidence_use_code(use_kind))?;
    bytes.extend_from_slice(source_packet_identity.as_bytes());
    bytes.extend_from_slice(source_assessment_identity.as_bytes());
    bytes.extend_from_slice(source_design_set_identity.as_bytes());
    write_len(&mut bytes, source_point_bytes.len())?;
    bytes.extend_from_slice(source_point_bytes);
    Ok(bytes)
}

impl PrerequisiteAssessmentReceipt {
    fn new(
        assessment: &ClaimPolicyAssessment,
        dependent: EulerClaimKind,
        use_kind: EvidenceUse,
    ) -> Result<Self, ContractError> {
        if assessment.disposition != AssessmentDisposition::ReferenceCompleteCandidate
            || assessment.reported_scientific_disposition != ReportedScientificDisposition::Positive
        {
            return Err(ContractError::new(
                "EulerProtocolIneligiblePrerequisite",
                "only a structurally complete caller-reported-positive assessment can back a prerequisite receipt",
            ));
        }
        assessment.verify_identity()?;
        let canonical = prerequisite_receipt_bytes(
            EULER_PROTOCOL_SCHEMA_VERSION,
            assessment.contract_identity,
            assessment.claim,
            dependent,
            use_kind,
            assessment.packet_identity,
            assessment.identity,
            assessment.design_set_identity,
            &assessment.point_bytes,
        )?;
        let identity = prerequisite_assessment_receipt_identity(&canonical);
        Ok(Self {
            schema_version: EULER_PROTOCOL_SCHEMA_VERSION,
            contract_identity: assessment.contract_identity,
            prerequisite: assessment.claim,
            dependent,
            use_kind,
            source_packet_identity: assessment.packet_identity,
            source_assessment_identity: assessment.identity,
            source_design_set_identity: assessment.design_set_identity,
            source_point_bytes: assessment.point_bytes.clone(),
            identity,
        })
    }

    /// Exact canonical receipt bytes.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ContractError> {
        prerequisite_receipt_bytes(
            self.schema_version,
            self.contract_identity,
            self.prerequisite,
            self.dependent,
            self.use_kind,
            self.source_packet_identity,
            self.source_assessment_identity,
            self.source_design_set_identity,
            &self.source_point_bytes,
        )
    }

    /// Recompute and verify every stored receipt binding.
    pub fn verify(&self) -> Result<(), ContractError> {
        let expected = prerequisite_assessment_receipt_identity(&self.canonical_bytes()?);
        if expected != self.identity {
            return Err(ContractError::new(
                "EulerProtocolPrerequisiteReceiptIdentityMismatch",
                "prerequisite receipt identity does not match its semantic fields",
            ));
        }
        Ok(())
    }

    #[must_use]
    /// Direct prerequisite claim.
    pub const fn prerequisite(&self) -> EulerClaimKind {
        self.prerequisite
    }

    #[must_use]
    /// Direct dependent claim.
    pub const fn dependent(&self) -> EulerClaimKind {
        self.dependent
    }

    #[must_use]
    /// Generic calibration/validation use on the exact edge.
    pub const fn use_kind(&self) -> EvidenceUse {
        self.use_kind
    }

    #[must_use]
    /// Exact design set covered by the prerequisite assessment.
    pub const fn source_design_set_identity(&self) -> ContentHash {
        self.source_design_set_identity
    }

    #[must_use]
    /// Domain-separated receipt identity.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

fn point_violations(context: &ContextOfUse, point: &ApplicabilityPoint) -> Vec<String> {
    let mut violations = Vec::new();
    for (axis, domain) in context.applicability().numeric() {
        match point.numeric().get(axis) {
            Some(value)
                if {
                    let (lo, hi) = domain.bounds();
                    value.is_finite() && *value >= lo && *value <= hi
                } => {}
            Some(_) => violations.push(format!("out-of-domain-numeric:{}", axis.as_str())),
            None => violations.push(format!("missing-numeric-axis:{}", axis.as_str())),
        }
    }
    for (axis, domain) in context.applicability().categorical() {
        match point.categorical().get(axis) {
            Some(value) if domain.allowed().contains(value) => {}
            Some(_) => violations.push(format!("out-of-domain-category:{}", axis.as_str())),
            None => violations.push(format!("missing-categorical-axis:{}", axis.as_str())),
        }
    }
    violations
}

fn validated_regime_covers_point(
    regime: &ValidityDomain,
    context: &ContextOfUse,
    point: &ApplicabilityPoint,
) -> bool {
    let point_covers_bound = |axis: &str, lo: f64, hi: f64| {
        lo.is_finite()
            && hi.is_finite()
            && point
                .numeric()
                .iter()
                .find(|(point_axis, _)| point_axis.as_str() == axis)
                .is_some_and(|(_, value)| value.is_finite() && *value >= lo && *value <= hi)
    };

    context.applicability().numeric().keys().all(|axis| {
        regime
            .bound(axis.as_str())
            .is_some_and(|(lo, hi)| point_covers_bound(axis.as_str(), lo, hi))
    }) && regime
        .bounds()
        .iter()
        .all(|(axis, (lo, hi))| point_covers_bound(axis, *lo, *hi))
}

fn evidence_weakness(
    record: &EvidenceRecord,
    context: &ContextOfUse,
    point: &ApplicabilityPoint,
) -> Option<String> {
    let expected = record.requirement.authority_class();
    if record.authority.class() != expected {
        return Some(format!(
            "weak-authority:{}:requires-{}:observed-{}",
            record.requirement.code(),
            expected.code(),
            record.authority.class().code()
        ));
    }

    if expected == EvidenceAuthorityClass::VerifiedNumerics {
        let EvidenceAuthorityDeclaration::VerifiedNumerics {
            color: Color::Verified { lo, hi },
        } = &record.authority
        else {
            return Some(format!(
                "weak-authority:{}:requires-finite-verified-color",
                record.requirement.code()
            ));
        };
        if !lo.is_finite() || !hi.is_finite() {
            return Some(format!(
                "weak-authority:{}:verified-enclosure-is-vacuous",
                record.requirement.code()
            ));
        }
    }

    if expected == EvidenceAuthorityClass::ValidatedPhysical {
        let EvidenceAuthorityDeclaration::ValidatedPhysical {
            color: Color::Validated { regime, .. },
        } = &record.authority
        else {
            return Some(format!(
                "weak-authority:{}:requires-validated-color",
                record.requirement.code()
            ));
        };
        if !validated_regime_covers_point(regime, context, point) {
            return Some(format!(
                "weak-validity-domain:{}:does-not-cover-case",
                record.requirement.code()
            ));
        }
    }
    if matches!(
        record.requirement,
        EvidenceRequirement::IndependentReconstruction
            | EvidenceRequirement::PhysicalValidation
            | EvidenceRequirement::BlindHoldout
            | EvidenceRequirement::RivalMechanismDiscrimination
    ) && !record.independent
    {
        return Some(format!(
            "weak-independence:{}:independent-evidence-required",
            record.requirement.code()
        ));
    }
    None
}

fn access_class_violation(record: &EvidenceRecord) -> Option<String> {
    use DeclaredEvidenceAccessClass as P;
    use EvidenceRequirement as E;
    let expected = match record.requirement {
        E::CalibrationPartition => P::Calibration,
        E::PhysicalValidation | E::RivalMechanismDiscrimination => P::Validation,
        E::BlindHoldout => P::BlindHoldout,
        _ => P::NotApplicable,
    };
    (record.access_class != expected).then(|| {
        format!(
            "access-class-mismatch:{}:expected-{}:observed-{}",
            record.requirement.code(),
            expected.code(),
            record.access_class.code()
        )
    })
}

fn claim_kind_for_id(id: &CampaignClaimId) -> Option<EulerClaimKind> {
    EULER_CLAIM_REGISTRY
        .into_iter()
        .find(|kind| kind.id() == id.as_str())
}

fn prerequisite_violations(
    contract: &EulerScientificContract,
    packet: &ClaimEvidencePacket,
    receipts: &[PrerequisiteAssessmentReceipt],
) -> Result<Vec<String>, ContractError> {
    let dependent_id = packet.claim.id();
    let expected = contract
        .claim_graph()
        .dependencies()
        .iter()
        .filter(|dependency| dependency.dependent.as_str() == dependent_id)
        .map(|dependency| {
            claim_kind_for_id(&dependency.prerequisite)
                .map(|kind| (kind, dependency.use_kind))
                .ok_or_else(|| {
                    ContractError::new(
                        "EulerProtocolUnknownDependencyClaim",
                        "frozen claim graph contains an unknown prerequisite identifier",
                    )
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let point_bytes = applicability_point_bytes(&packet.point)?;
    // Caller slice order is not semantic. Canonicalize the complete receipt
    // tuple before the first refusal is observed so both `first_divergence`
    // and the assessment/log identities remain permutation invariant even
    // for multiple malformed or unexpected receipts. The cached identity is
    // deliberately the final tie-breaker because `verify` may reject it.
    let mut ordered_receipts = receipts.iter().collect::<Vec<_>>();
    ordered_receipts.sort_by(|left, right| {
        (
            left.schema_version,
            left.contract_identity.as_hash(),
            left.prerequisite.id(),
            left.dependent.id(),
            evidence_use_code(left.use_kind),
            left.source_packet_identity,
            left.source_assessment_identity,
            left.source_design_set_identity,
            left.source_point_bytes.as_slice(),
            left.identity,
        )
            .cmp(&(
                right.schema_version,
                right.contract_identity.as_hash(),
                right.prerequisite.id(),
                right.dependent.id(),
                evidence_use_code(right.use_kind),
                right.source_packet_identity,
                right.source_assessment_identity,
                right.source_design_set_identity,
                right.source_point_bytes.as_slice(),
                right.identity,
            ))
    });
    let mut observed = BTreeSet::new();
    let mut violations = Vec::new();
    for receipt in ordered_receipts {
        if let Err(error) = receipt.verify() {
            violations.push(format!("malformed-prerequisite-receipt:{}", error.code()));
            continue;
        }
        let key = (receipt.prerequisite, receipt.use_kind);
        let mut exact_binding = true;
        if receipt.schema_version != EULER_PROTOCOL_SCHEMA_VERSION
            || receipt.contract_identity != contract.identity()
        {
            exact_binding = false;
            violations.push(format!(
                "stale-prerequisite-receipt:{}",
                receipt.prerequisite.id()
            ));
        }
        if receipt.dependent != packet.claim || !expected.contains(&key) {
            exact_binding = false;
            violations.push(format!(
                "unexpected-prerequisite-receipt:{}:{}",
                receipt.prerequisite.id(),
                evidence_use_code(receipt.use_kind)
            ));
        }
        if receipt.source_point_bytes != point_bytes {
            exact_binding = false;
            violations.push(format!(
                "prerequisite-applicability-point-mismatch:{}",
                receipt.prerequisite.id()
            ));
        }
        if receipt.source_design_set_identity != packet.design_set_identity {
            exact_binding = false;
            violations.push(format!(
                "prerequisite-design-set-mismatch:{}",
                receipt.prerequisite.id()
            ));
        }
        if exact_binding && !observed.insert(key) {
            violations.push(format!(
                "duplicate-prerequisite-receipt:{}:{}",
                receipt.prerequisite.id(),
                evidence_use_code(receipt.use_kind)
            ));
        }
    }
    for (kind, use_kind) in expected {
        if !observed.contains(&(kind, use_kind)) {
            violations.push(format!(
                "missing-prerequisite-receipt:{}:{}",
                kind.id(),
                evidence_use_code(use_kind)
            ));
        }
    }
    Ok(violations)
}

#[allow(clippy::too_many_lines)] // One pass preserves refusal/demotion precedence for audit.
fn assess_packet(
    contract: &EulerScientificContract,
    packet: &ClaimEvidencePacket,
    prerequisites: &[PrerequisiteAssessmentReceipt],
) -> Result<ClaimPolicyAssessment, ContractError> {
    // Safe constructors always publish the current protocol version. Keep the
    // evaluator fail-closed for internally forged or future decoded packets,
    // but do not mint a retained reason that the v1 log cannot independently
    // bind: the log's own schema version is necessarily current.
    protocol_migration_policy(packet.schema_version)?;
    let mut hard = Vec::new();
    let mut weak = Vec::new();
    if packet.contract_identity != contract.identity() {
        hard.push("contract-identity-mismatch".to_owned());
    }
    hard.extend(prerequisite_violations(contract, packet, prerequisites)?);
    let Some(claim) = contract.claim_graph().claim(packet.claim) else {
        hard.push("claim-not-present-in-contract".to_owned());
        return ClaimPolicyAssessment::build(
            contract,
            packet,
            prerequisites,
            AssessmentDisposition::Refused,
            hard,
        );
    };
    if !packet.no_claims_accepted {
        hard.push("binding-no-claims-not-accepted".to_owned());
    }
    if packet.target_fitted && packet.claim.forbids_target_fitting() {
        hard.push("protected-target-fitting-invalidates-emergent-claim".to_owned());
    }
    hard.extend(point_violations(contract.context(), &packet.point));

    let expected_qois = claim.campaign().qois.as_slice();
    let mut expected_units = expected_qois
        .iter()
        .filter_map(|qoi| contract.context().qois().get(qoi))
        .map(|qoi| qoi.unit().as_str().to_owned())
        .collect::<Vec<_>>();
    expected_units.sort();
    expected_units.dedup();
    if packet.units != expected_units {
        // Comma is deliberately the v1 unit-list separator because
        // `checked_protocol_id` forbids it inside a unit identifier. `+` is
        // legal in identifiers and would make ["1", "j"] collide with
        // ["1+j"] in the retained reason grammar.
        hard.push(format!(
            "claim-unit-set-mismatch:expected-{}:observed-{}",
            expected_units.join(","),
            packet.units.join(",")
        ));
    }
    let hypothesis_hashes = contract
        .extension()
        .hypothesis_sources()
        .iter()
        .map(HypothesisSource::declaration_hash)
        .collect::<BTreeSet<_>>();
    for (role, identity) in [
        ("design-set", packet.design_set_identity),
        (
            "aggregate-qoi-derivation-receipt",
            packet.aggregate_qoi_derivation_receipt_identity,
        ),
    ] {
        if hypothesis_hashes.contains(&identity) {
            hard.push(format!(
                "hypothesis-source-cannot-satisfy-packet-role:{role}"
            ));
        }
    }
    // The no-hypothesis-as-evidence boundary applies to every retained row,
    // including rows whose role is unexpected for this claim.  Scan the
    // packet, rather than only the required-role loop below, so a structurally
    // valid but unexpected collision is represented as a refusal that the
    // strict assessment-log reader can reproduce.
    for (requirement, record) in &packet.records {
        for (slot, reference) in evidence_hash_references(record) {
            if reference.is_some_and(|hash| hypothesis_hashes.contains(&hash)) {
                hard.push(format!(
                    "hypothesis-source-cannot-satisfy-evidence:{}:{slot}",
                    requirement.code()
                ));
            }
        }
    }
    for requirement in claim.requirements() {
        let Some(record) = packet.records.get(requirement) else {
            hard.push(format!("missing-evidence:{}", requirement.code()));
            continue;
        };
        if record.contract_identity != contract.identity() || record.claim != packet.claim {
            hard.push(format!("stale-evidence-binding:{}", requirement.code()));
        }
        if record.qois != expected_qois {
            hard.push(format!("qoi-binding-mismatch:{}", requirement.code()));
        }
        if record.source_schema != requirement.source_schema() {
            hard.push(format!("source-schema-mismatch:{}", requirement.code()));
        }
        if record.source_kind != requirement.source_kind() {
            hard.push(format!("source-kind-mismatch:{}", requirement.code()));
        }
        if let Some(reason) = access_class_violation(record) {
            hard.push(reason);
        }
        if let Some(reason) = evidence_weakness(record, contract.context(), &packet.point) {
            weak.push(reason);
        }
    }
    for (unexpected, record) in packet
        .records
        .iter()
        .filter(|(requirement, _)| !claim.requirements().contains(requirement))
    {
        hard.push(format!("unexpected-evidence:{}", unexpected.code()));
        if let Some(reason) = evidence_weakness(record, contract.context(), &packet.point) {
            weak.push(reason);
        }
    }
    let disposition = if !hard.is_empty() {
        AssessmentDisposition::Refused
    } else if packet.reported_scientific_disposition != ReportedScientificDisposition::Positive {
        // A weakness limits the referenced support attached to a retained
        // caller-reported negative or inconclusive result, but cannot turn
        // that terminal result back into a candidate state. Keep the
        // diagnostics below without erasing or promoting the report.
        AssessmentDisposition::RetainedTerminal
    } else if weak.is_empty() {
        AssessmentDisposition::ReferenceCompleteCandidate
    } else {
        AssessmentDisposition::DemotedCandidate
    };
    hard.extend(weak);
    if packet.expected_disposition != disposition {
        hard.push(format!(
            "expected-disposition-mismatch:expected-{}:observed-{}",
            packet.expected_disposition.code(),
            disposition.code()
        ));
    }
    ClaimPolicyAssessment::build(contract, packet, prerequisites, disposition, hard)
}

/// Deterministic, bounded, redacted JSON-lines record for one local
/// claim-policy assessment.
///
/// This is not the campaign-wide evidence-event or retained-log protocol.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClaimPolicyAssessmentLog {
    json_line: String,
    identity: ContentHash,
}

const MAX_ASSESSMENT_LOG_REASON_ROWS: usize = MAX_EVIDENCE_RECORDS * 8;
// An evidence-source row concatenates four individually bounded identifiers
// and three delimiters. Closed requirement/kind/schema values are shorter,
// but using their public identifier bound keeps the reader coupled to the
// constructor contract rather than to incidental current spellings.
const MAX_ASSESSMENT_LOG_EVIDENCE_SOURCE_BYTES: usize = MAX_PROTOCOL_ID_BYTES * 4 + 3;
const MAX_ASSESSMENT_LOG_UNIT_ROWS: usize = MAX_EULER_CLAIMS * 2;
const MAX_ASSESSMENT_LOG_UNIT_LIST_BYTES: usize =
    MAX_ASSESSMENT_LOG_UNIT_ROWS * MAX_PROTOCOL_ID_BYTES + MAX_ASSESSMENT_LOG_UNIT_ROWS - 1;
// A unit-set mismatch joins the complete bounded observed unit list into one
// reason. The remaining allowance covers the fixed prefix and the small
// frozen expected list. The 32-KiB whole-line ceiling remains the final bound.
const MAX_ASSESSMENT_LOG_REASON_BYTES: usize =
    MAX_ASSESSMENT_LOG_UNIT_LIST_BYTES + MAX_PROTOCOL_ID_BYTES;
const MAX_ASSESSMENT_LOG_ARTIFACT_ROWS: usize =
    MAX_EVIDENCE_RECORDS * 3 + MAX_PREREQUISITE_RECEIPTS + 3;
fn malformed_assessment_log(detail: impl Into<String>) -> ContractError {
    ContractError::new("EulerProtocolMalformedAssessmentLog", detail)
}

struct AssessmentLogReader<'a> {
    input: &'a str,
    offset: usize,
}

impl<'a> AssessmentLogReader<'a> {
    fn new(input: &'a str) -> Self {
        Self { input, offset: 0 }
    }

    fn expect(&mut self, literal: &str, field: &'static str) -> Result<(), ContractError> {
        if !self.input[self.offset..].starts_with(literal) {
            return Err(malformed_assessment_log(format!(
                "assessment log has a missing, reordered, or mistyped {field} field"
            )));
        }
        self.offset += literal.len();
        Ok(())
    }

    fn string(
        &mut self,
        field: &'static str,
        max_decoded_bytes: usize,
    ) -> Result<String, ContractError> {
        let token_start = self.offset;
        self.expect("\"", field)?;
        let mut decoded = String::new();
        loop {
            let Some(byte) = self.input.as_bytes().get(self.offset).copied() else {
                return Err(malformed_assessment_log(format!(
                    "assessment log {field} string is truncated"
                )));
            };
            match byte {
                b'"' => {
                    self.offset += 1;
                    break;
                }
                b'\\' => {
                    self.offset += 1;
                    let Some(escape) = self.input.as_bytes().get(self.offset).copied() else {
                        return Err(malformed_assessment_log(format!(
                            "assessment log {field} escape is truncated"
                        )));
                    };
                    self.offset += 1;
                    match escape {
                        b'"' => decoded.push('"'),
                        b'\\' => decoded.push('\\'),
                        b'n' => decoded.push('\n'),
                        b'r' => decoded.push('\r'),
                        b't' => decoded.push('\t'),
                        b'u' => {
                            let end = self.offset.checked_add(4).ok_or_else(|| {
                                malformed_assessment_log(format!(
                                    "assessment log {field} unicode escape overflows"
                                ))
                            })?;
                            let Some(hex) = self.input.get(self.offset..end) else {
                                return Err(malformed_assessment_log(format!(
                                    "assessment log {field} unicode escape is truncated"
                                )));
                            };
                            if !hex
                                .bytes()
                                .all(|digit| digit.is_ascii_digit() || matches!(digit, b'a'..=b'f'))
                            {
                                return Err(malformed_assessment_log(format!(
                                    "assessment log {field} unicode escape is not canonical lowercase hexadecimal"
                                )));
                            }
                            let scalar = u32::from_str_radix(hex, 16).map_err(|_| {
                                malformed_assessment_log(format!(
                                    "assessment log {field} unicode escape is invalid"
                                ))
                            })?;
                            let character = char::from_u32(scalar).ok_or_else(|| {
                                malformed_assessment_log(format!(
                                    "assessment log {field} unicode escape is not a scalar value"
                                ))
                            })?;
                            decoded.push(character);
                            self.offset = end;
                        }
                        _ => {
                            return Err(malformed_assessment_log(format!(
                                "assessment log {field} uses a noncanonical escape"
                            )));
                        }
                    }
                }
                0x00..=0x1f => {
                    return Err(malformed_assessment_log(format!(
                        "assessment log {field} contains an unescaped control byte"
                    )));
                }
                0x20..=0x7f => {
                    decoded.push(char::from(byte));
                    self.offset += 1;
                }
                _ => {
                    let character = self.input[self.offset..].chars().next().ok_or_else(|| {
                        malformed_assessment_log(format!(
                            "assessment log {field} contains invalid UTF-8"
                        ))
                    })?;
                    decoded.push(character);
                    self.offset += character.len_utf8();
                }
            }
            if decoded.len() > max_decoded_bytes {
                return Err(malformed_assessment_log(format!(
                    "assessment log {field} exceeds its decoded byte bound"
                )));
            }
        }
        let mut canonical = String::new();
        json_string(&mut canonical, &decoded);
        if self.input.get(token_start..self.offset) != Some(canonical.as_str()) {
            return Err(malformed_assessment_log(format!(
                "assessment log {field} string is not canonically escaped"
            )));
        }
        Ok(decoded)
    }

    fn unsigned(&mut self, field: &'static str) -> Result<u64, ContractError> {
        let start = self.offset;
        while self
            .input
            .as_bytes()
            .get(self.offset)
            .is_some_and(u8::is_ascii_digit)
        {
            self.offset += 1;
        }
        let digits = self.input.get(start..self.offset).unwrap_or_default();
        if digits.is_empty() || (digits.len() > 1 && digits.starts_with('0')) {
            return Err(malformed_assessment_log(format!(
                "assessment log {field} is not a canonical unsigned integer"
            )));
        }
        digits.parse::<u64>().map_err(|_| {
            malformed_assessment_log(format!("assessment log {field} unsigned integer overflows"))
        })
    }

    fn boolean(&mut self, field: &'static str) -> Result<bool, ContractError> {
        if self.input[self.offset..].starts_with("true") {
            self.offset += 4;
            Ok(true)
        } else if self.input[self.offset..].starts_with("false") {
            self.offset += 5;
            Ok(false)
        } else {
            Err(malformed_assessment_log(format!(
                "assessment log {field} is not a canonical boolean"
            )))
        }
    }

    fn string_array(
        &mut self,
        field: &'static str,
        max_items: usize,
        max_item_bytes: usize,
    ) -> Result<Vec<String>, ContractError> {
        self.expect("[", field)?;
        let mut values = Vec::new();
        if self.input[self.offset..].starts_with(']') {
            self.offset += 1;
            return Ok(values);
        }
        loop {
            if values.len() == max_items {
                return Err(malformed_assessment_log(format!(
                    "assessment log {field} exceeds its row bound"
                )));
            }
            values.push(self.string(field, max_item_bytes)?);
            if self.input[self.offset..].starts_with(']') {
                self.offset += 1;
                return Ok(values);
            }
            self.expect(",", field)?;
        }
    }

    fn finish(self) -> Result<(), ContractError> {
        if self.offset == self.input.len() {
            Ok(())
        } else {
            Err(malformed_assessment_log(
                "assessment log contains extra, duplicate, or trailing fields",
            ))
        }
    }
}

fn is_strictly_sorted(values: &[String]) -> bool {
    values.windows(2).all(|pair| pair[0] < pair[1])
}

fn is_exact_lower_hex(value: &str, digits: usize) -> bool {
    value.len() == digits
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || matches!(byte, b'a'..=b'f'))
}

fn is_nonzero_lower_content_hash(value: &str) -> bool {
    is_exact_lower_hex(value, 64) && value.bytes().any(|byte| byte != b'0')
}

fn assessment_log_requirement(code: &str) -> Option<EvidenceRequirement> {
    EULER_EVIDENCE_REQUIREMENT_REGISTRY
        .into_iter()
        .find(|requirement| requirement.code() == code)
}

fn assessment_log_claim_kind(id: &str) -> Option<EulerClaimKind> {
    EULER_CLAIM_REGISTRY
        .iter()
        .copied()
        .find(|claim| claim.id() == id)
}

fn validate_assessment_log_evidence_sources(
    values: &[String],
) -> Result<BTreeSet<EvidenceRequirement>, ContractError> {
    let mut requirements = BTreeSet::new();
    for value in values {
        let (requirement_code, remainder) = value.split_once(':').ok_or_else(|| {
            malformed_assessment_log(
                "assessment log evidence source omits its requirement delimiter",
            )
        })?;
        let requirement = assessment_log_requirement(requirement_code).ok_or_else(|| {
            malformed_assessment_log(
                "assessment log evidence source carries an unknown evidence requirement",
            )
        })?;
        let (source_kind, schema_and_id) = remainder.split_once(':').ok_or_else(|| {
            malformed_assessment_log("assessment log evidence source omits its kind delimiter")
        })?;
        if source_kind != requirement.source_kind().slug() {
            return Err(malformed_assessment_log(
                "assessment log evidence source kind is not the exact requirement container",
            ));
        }
        let source_id = schema_and_id
            .strip_prefix(requirement.source_schema())
            .and_then(|remainder| remainder.strip_prefix(':'))
            .ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log evidence source schema is not the exact requirement route",
                )
            })?;
        if validate_protocol_id("assessment_log.evidence_source_id", source_id).is_err() {
            return Err(malformed_assessment_log(
                "assessment log evidence source id is not a canonical machine identity",
            ));
        }
        if !requirements.insert(requirement) {
            return Err(malformed_assessment_log(
                "assessment log carries multiple source rows for one evidence requirement",
            ));
        }
    }
    Ok(requirements)
}

const ASSESSMENT_LOG_ARTIFACT_SLOT: u8 = 1;
const ASSESSMENT_LOG_SCHEMA_RECEIPT_SLOT: u8 = 1 << 1;
const ASSESSMENT_LOG_ROLE_RECEIPT_SLOT: u8 = 1 << 2;
const ASSESSMENT_LOG_BASE_EVIDENCE_SLOTS: u8 =
    ASSESSMENT_LOG_ARTIFACT_SLOT | ASSESSMENT_LOG_SCHEMA_RECEIPT_SLOT;

fn assessment_log_evidence_slot(slot: &str) -> Option<u8> {
    match slot {
        "artifact" => Some(ASSESSMENT_LOG_ARTIFACT_SLOT),
        "schema-admission-receipt" => Some(ASSESSMENT_LOG_SCHEMA_RECEIPT_SLOT),
        "role-receipt" => Some(ASSESSMENT_LOG_ROLE_RECEIPT_SLOT),
        _ => None,
    }
}

fn assessment_log_evidence_slot_code(slot: u8) -> Option<&'static str> {
    match slot {
        ASSESSMENT_LOG_ARTIFACT_SLOT => Some("artifact"),
        ASSESSMENT_LOG_SCHEMA_RECEIPT_SLOT => Some("schema-admission-receipt"),
        ASSESSMENT_LOG_ROLE_RECEIPT_SLOT => Some("role-receipt"),
        _ => None,
    }
}

fn expected_assessment_log_prerequisites(
    contract: &EulerScientificContract,
    claim: EulerClaimKind,
) -> Result<BTreeSet<(EulerClaimKind, EvidenceUse)>, ContractError> {
    contract
        .claim_graph()
        .dependencies()
        .iter()
        .filter(|dependency| dependency.dependent.as_str() == claim.id())
        .map(|dependency| {
            claim_kind_for_id(&dependency.prerequisite)
                .map(|kind| (kind, dependency.use_kind))
                .ok_or_else(|| {
                    malformed_assessment_log(
                        "frozen claim graph contains an unknown prerequisite identifier",
                    )
                })
        })
        .collect()
}

#[allow(clippy::too_many_lines)]
fn validate_assessment_log_relative_artifacts(
    values: &[String],
    packet_identity: &str,
    design_set_identity: &str,
    aggregate_qoi_derivation_receipt_identity: &str,
    claim: EulerClaimKind,
    observed_disposition: &str,
    evidence_sources: &BTreeSet<EvidenceRequirement>,
    reasons: &[String],
    contract: &EulerScientificContract,
) -> Result<(), ContractError> {
    let mut packet_rows = 0_usize;
    let mut design_set_rows = 0_usize;
    let mut aggregate_derivation_rows = 0_usize;
    let mut evidence_slots = BTreeMap::<EvidenceRequirement, u8>::new();
    let mut evidence_references = BTreeMap::<(EvidenceRequirement, u8), ContentHash>::new();
    let mut evidence_hashes = BTreeSet::<ContentHash>::new();
    let mut prerequisite_rows = BTreeMap::<EulerClaimKind, usize>::new();
    let mut prerequisite_hashes = BTreeSet::<ContentHash>::new();
    for value in values {
        if let Some(hash) = value.strip_prefix("packet:") {
            if hash != packet_identity {
                return Err(malformed_assessment_log(
                    "assessment log packet artifact is not its exact packet identity",
                ));
            }
            packet_rows += 1;
            continue;
        }
        if let Some(hash) = value.strip_prefix("design-set:") {
            if hash != design_set_identity {
                return Err(malformed_assessment_log(
                    "assessment log design-set artifact is not its exact top-level identity",
                ));
            }
            design_set_rows += 1;
            continue;
        }
        if let Some(remainder) = value.strip_prefix("aggregate-qoi-derivation:") {
            let hash = remainder
                .strip_prefix(EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA)
                .and_then(|remainder| remainder.strip_prefix(':'))
                .ok_or_else(|| {
                    malformed_assessment_log(
                        "assessment log aggregate-QoI derivation receipt is not labelled with the exact owner route",
                    )
                })?;
            if hash != aggregate_qoi_derivation_receipt_identity {
                return Err(malformed_assessment_log(
                    "assessment log aggregate-QoI derivation artifact is not its exact top-level identity",
                ));
            }
            aggregate_derivation_rows += 1;
            continue;
        }
        if let Some(remainder) = value.strip_prefix("evidence:") {
            let (requirement_code, slot_and_hash) = remainder.split_once(':').ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log evidence artifact omits its requirement delimiter",
                )
            })?;
            let requirement = assessment_log_requirement(requirement_code).ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log evidence artifact carries an unknown requirement",
                )
            })?;
            if !evidence_sources.contains(&requirement) {
                return Err(malformed_assessment_log(
                    "assessment log evidence artifact has no matching evidence source row",
                ));
            }
            let (slot, hash) = slot_and_hash.split_once(':').ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log evidence artifact omits its slot delimiter",
                )
            })?;
            let slot = assessment_log_evidence_slot(slot).ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log evidence artifact carries an unknown identity slot",
                )
            })?;
            if !is_nonzero_lower_content_hash(hash) {
                return Err(malformed_assessment_log(
                    "assessment log evidence artifact is not a nonzero lowercase 256-bit identity",
                ));
            }
            let hash = ContentHash::from_hex(hash).ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log evidence artifact identity cannot be decoded",
                )
            })?;
            let retained_slots = evidence_slots.entry(requirement).or_default();
            if *retained_slots & slot != 0 {
                return Err(malformed_assessment_log(
                    "assessment log carries multiple identities for one evidence slot",
                ));
            }
            *retained_slots |= slot;
            evidence_references.insert((requirement, slot), hash);
            if !evidence_hashes.insert(hash) {
                return Err(malformed_assessment_log(
                    "assessment log reuses one logical hash across evidence slots",
                ));
            }
            continue;
        }
        if let Some(remainder) = value.strip_prefix("prerequisite:") {
            let (prerequisite, route_and_hash) = remainder.split_once(':').ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log prerequisite artifact omits its claim delimiter",
                )
            })?;
            let prerequisite = assessment_log_claim_kind(prerequisite).ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log prerequisite artifact carries an unknown claim",
                )
            })?;
            let hash = route_and_hash
                .strip_prefix(EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN)
                .and_then(|remainder| remainder.strip_prefix(':'))
                .ok_or_else(|| {
                    malformed_assessment_log(
                        "assessment log prerequisite artifact is not labelled with the exact owner route",
                    )
                })?;
            if !is_nonzero_lower_content_hash(hash) {
                return Err(malformed_assessment_log(
                    "assessment log prerequisite artifact is not a nonzero lowercase 256-bit identity",
                ));
            }
            let hash = ContentHash::from_hex(hash).ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log prerequisite artifact identity cannot be decoded",
                )
            })?;
            if !prerequisite_hashes.insert(hash) {
                return Err(malformed_assessment_log(
                    "assessment log reuses one prerequisite receipt identity under multiple labels",
                ));
            }
            *prerequisite_rows.entry(prerequisite).or_default() += 1;
            continue;
        }
        return Err(malformed_assessment_log(
            "assessment log relative artifact carries an unknown row kind",
        ));
    }
    if packet_rows != 1 {
        return Err(malformed_assessment_log(
            "assessment log must retain exactly one packet artifact identity",
        ));
    }
    if design_set_rows != 1 || aggregate_derivation_rows != 1 {
        return Err(malformed_assessment_log(
            "assessment log must retain exactly one design-set and one aggregate-QoI derivation-receipt identity",
        ));
    }
    if design_set_identity == aggregate_qoi_derivation_receipt_identity {
        return Err(malformed_assessment_log(
            "assessment log cannot alias its design set and aggregate-QoI derivation receipt",
        ));
    }
    if evidence_hashes.iter().any(|hash| {
        let hex = hash.to_hex();
        hex == design_set_identity || hex == aggregate_qoi_derivation_receipt_identity
    }) {
        return Err(malformed_assessment_log(
            "assessment log cannot alias design-set or aggregate-QoI derivation identities into evidence slots",
        ));
    }
    let artifact_requirements = evidence_slots.keys().copied().collect::<BTreeSet<_>>();
    if &artifact_requirements != evidence_sources {
        return Err(malformed_assessment_log(
            "assessment log evidence source and artifact role sets differ",
        ));
    }
    for (requirement, slots) in &evidence_slots {
        if *slots & ASSESSMENT_LOG_BASE_EVIDENCE_SLOTS != ASSESSMENT_LOG_BASE_EVIDENCE_SLOTS {
            return Err(malformed_assessment_log(
                "assessment log evidence row omits its artifact or schema-admission receipt identity",
            ));
        }
        let observed_authority = observed_assessment_log_authority_class(*requirement, reasons)?;
        let expected_slots = if observed_authority == EvidenceAuthorityClass::StructuralProcess {
            ASSESSMENT_LOG_BASE_EVIDENCE_SLOTS | ASSESSMENT_LOG_ROLE_RECEIPT_SLOT
        } else {
            ASSESSMENT_LOG_BASE_EVIDENCE_SLOTS
        };
        if *slots != expected_slots {
            return Err(malformed_assessment_log(
                "assessment log evidence row has an impossible role-receipt slot shape",
            ));
        }
    }
    let hypothesis_hashes = contract
        .extension()
        .hypothesis_sources()
        .iter()
        .map(HypothesisSource::declaration_hash)
        .collect::<BTreeSet<_>>();
    for (role, identity) in [
        ("design-set", design_set_identity),
        (
            "aggregate-qoi-derivation-receipt",
            aggregate_qoi_derivation_receipt_identity,
        ),
    ] {
        let reason = format!("hypothesis-source-cannot-satisfy-packet-role:{role}");
        let collision = hypothesis_hashes
            .iter()
            .any(|hypothesis| hypothesis.to_hex() == identity);
        let reason_present = reasons.binary_search(&reason).is_ok();
        if collision != reason_present
            || (collision && observed_disposition != AssessmentDisposition::Refused.code())
        {
            return Err(malformed_assessment_log(
                "assessment log hypothesis-source collision is not bound bidirectionally to its exact packet role and refusal",
            ));
        }
    }
    for ((requirement, slot), hash) in &evidence_references {
        let slot = assessment_log_evidence_slot_code(*slot).ok_or_else(|| {
            malformed_assessment_log(
                "assessment log retained an unrepresentable evidence identity slot",
            )
        })?;
        let reason = format!(
            "hypothesis-source-cannot-satisfy-evidence:{}:{slot}",
            requirement.code()
        );
        let collision = hypothesis_hashes.contains(hash);
        let reason_present = reasons.binary_search(&reason).is_ok();
        if collision != reason_present
            || (collision && observed_disposition != AssessmentDisposition::Refused.code())
        {
            return Err(malformed_assessment_log(
                "assessment log hypothesis-source collision is not bound bidirectionally to its exact evidence slot and refusal",
            ));
        }
    }

    let expected_prerequisites = expected_assessment_log_prerequisites(contract, claim)?;
    let expected_prerequisite_counts = expected_prerequisites.iter().fold(
        BTreeMap::<EulerClaimKind, usize>::new(),
        |mut counts, (prerequisite, _)| {
            *counts.entry(*prerequisite).or_default() += 1;
            counts
        },
    );
    let malformed_prerequisite_reason = reasons
        .iter()
        .any(|reason| reason.starts_with("malformed-prerequisite-receipt:"));
    let prerequisite_reason_claims = |prefix: &str| {
        reasons
            .iter()
            .filter_map(|reason| reason.strip_prefix(prefix))
            .filter_map(|remainder| remainder.split_once(':').map(|(claim, _)| claim))
            .filter_map(assessment_log_claim_kind)
            .collect::<BTreeSet<_>>()
    };
    let unexpected_prerequisite_claims =
        prerequisite_reason_claims("unexpected-prerequisite-receipt:");
    let duplicate_prerequisite_claims =
        prerequisite_reason_claims("duplicate-prerequisite-receipt:");
    let stale_prerequisite_claims = reasons
        .iter()
        .filter_map(|reason| reason.strip_prefix("stale-prerequisite-receipt:"))
        .filter_map(assessment_log_claim_kind)
        .collect::<BTreeSet<_>>();
    let point_mismatch_prerequisite_claims = reasons
        .iter()
        .filter_map(|reason| reason.strip_prefix("prerequisite-applicability-point-mismatch:"))
        .filter_map(assessment_log_claim_kind)
        .collect::<BTreeSet<_>>();
    let design_set_mismatch_prerequisite_claims = reasons
        .iter()
        .filter_map(|reason| reason.strip_prefix("prerequisite-design-set-mismatch:"))
        .filter_map(assessment_log_claim_kind)
        .collect::<BTreeSet<_>>();
    let mut weakness_roles = BTreeSet::new();
    for reason in reasons {
        if let Some(requirement) = assessment_log_weakness_requirement(reason)
            && !weakness_roles.insert(requirement)
        {
            return Err(malformed_assessment_log(
                "assessment log carries multiple mutually exclusive weakness reasons for one evidence role",
            ));
        }
        if let Some(remainder) = reason.strip_prefix("hypothesis-source-cannot-satisfy-evidence:") {
            let (requirement, slot) = remainder.split_once(':').ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log hypothesis-source reason omits its evidence slot",
                )
            })?;
            let requirement = assessment_log_requirement(requirement).ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log hypothesis-source reason carries an unknown evidence role",
                )
            })?;
            let slot = assessment_log_evidence_slot(slot).ok_or_else(|| {
                malformed_assessment_log(
                    "assessment log hypothesis-source reason carries an unknown evidence slot",
                )
            })?;
            if !evidence_references.contains_key(&(requirement, slot)) {
                return Err(malformed_assessment_log(
                    "assessment log hypothesis-source reason names an evidence slot that was not retained",
                ));
            }
        }
        if let Some(prerequisite) = reason
            .strip_prefix("stale-prerequisite-receipt:")
            .or_else(|| reason.strip_prefix("prerequisite-applicability-point-mismatch:"))
            .or_else(|| reason.strip_prefix("prerequisite-design-set-mismatch:"))
            .and_then(assessment_log_claim_kind)
            && !prerequisite_rows.contains_key(&prerequisite)
        {
            return Err(malformed_assessment_log(
                "assessment log prerequisite reason has no retained receipt artifact",
            ));
        }
        if reason.starts_with("malformed-prerequisite-receipt:") && prerequisite_rows.is_empty() {
            return Err(malformed_assessment_log(
                "assessment log malformed-prerequisite reason has no retained receipt artifact",
            ));
        }
        for prefix in [
            "unexpected-prerequisite-receipt:",
            "duplicate-prerequisite-receipt:",
            "missing-prerequisite-receipt:",
        ] {
            let Some(remainder) = reason.strip_prefix(prefix) else {
                continue;
            };
            let Some((prerequisite, use_kind)) = remainder.split_once(':') else {
                continue;
            };
            let (Some(prerequisite), Some(use_kind)) = (
                assessment_log_claim_kind(prerequisite),
                assessment_log_evidence_use(use_kind),
            ) else {
                continue;
            };
            if prefix == "missing-prerequisite-receipt:"
                && !expected_prerequisites.contains(&(prerequisite, use_kind))
            {
                return Err(malformed_assessment_log(
                    "assessment log missing-prerequisite reason is not an incoming frozen-DAG edge",
                ));
            }
            if prefix != "missing-prerequisite-receipt:"
                && !prerequisite_rows.contains_key(&prerequisite)
            {
                return Err(malformed_assessment_log(
                    "assessment log prerequisite reason has no retained receipt artifact",
                ));
            }
            if prefix == "duplicate-prerequisite-receipt:"
                && !expected_prerequisites.contains(&(prerequisite, use_kind))
            {
                return Err(malformed_assessment_log(
                    "assessment log duplicate-prerequisite reason is not an incoming frozen-DAG edge",
                ));
            }
        }
    }
    for (prerequisite, use_kind) in &expected_prerequisites {
        let retained_rows = prerequisite_rows
            .get(prerequisite)
            .copied()
            .unwrap_or_default();
        // With one frozen edge per prerequisite claim, one retained row that
        // is itself labelled stale, point-mismatched, design-set-mismatched,
        // or unexpected cannot
        // satisfy the edge. A malformed diagnostic is intentionally excluded:
        // the writer deduplicates artifact strings, so a valid exact receipt
        // and a malformed transplant carrying its cached identity can collapse
        // to the same retained row while the exact receipt still satisfies the
        // edge. More than one row may likewise include one exact receipt, so
        // the line alone cannot infer absence in that case.
        let sole_row_is_invalid = retained_rows == 1
            && (stale_prerequisite_claims.contains(prerequisite)
                || point_mismatch_prerequisite_claims.contains(prerequisite)
                || design_set_mismatch_prerequisite_claims.contains(prerequisite)
                || unexpected_prerequisite_claims.contains(prerequisite));
        if retained_rows == 0 || sole_row_is_invalid {
            let required_reason = format!(
                "missing-prerequisite-receipt:{}:{}",
                prerequisite.id(),
                evidence_use_code(*use_kind)
            );
            if reasons.binary_search(&required_reason).is_err()
                || observed_disposition != AssessmentDisposition::Refused.code()
            {
                return Err(malformed_assessment_log(
                    "assessment log absent prerequisite artifact is not bound to its exact missing reason and refusal",
                ));
            }
        }
    }
    for (prerequisite, observed_count) in &prerequisite_rows {
        let expected_count = expected_prerequisite_counts
            .get(prerequisite)
            .copied()
            .unwrap_or_default();
        // A global malformed-receipt diagnostic is deliberately not assigned
        // to a claim because a corrupt cached identity cannot be trusted as a
        // semantic binding. A stale, point-mismatch, or design-set-mismatch
        // diagnostic for this claim, however, proves that at least one row
        // verified far enough to
        // expose its semantic fields. For an out-of-DAG claim that verified
        // row must also have produced the exact unexpected-receipt reason;
        // an unrelated malformed row cannot suppress it.
        let malformed_may_explain_unexpected_claim = malformed_prerequisite_reason
            && !stale_prerequisite_claims.contains(prerequisite)
            && !point_mismatch_prerequisite_claims.contains(prerequisite)
            && !design_set_mismatch_prerequisite_claims.contains(prerequisite);
        if expected_count == 0
            && !malformed_may_explain_unexpected_claim
            && !unexpected_prerequisite_claims.contains(prerequisite)
        {
            return Err(malformed_assessment_log(
                "assessment log retains an unexpected prerequisite claim without its evaluator-required refusal reason",
            ));
        }
        if *observed_count > expected_count
            && !malformed_prerequisite_reason
            && !unexpected_prerequisite_claims.contains(prerequisite)
            && !duplicate_prerequisite_claims.contains(prerequisite)
            && !stale_prerequisite_claims.contains(prerequisite)
            && !point_mismatch_prerequisite_claims.contains(prerequisite)
            && !design_set_mismatch_prerequisite_claims.contains(prerequisite)
        {
            return Err(malformed_assessment_log(
                "assessment log retains excess prerequisite artifacts without a line-observable evaluator reason",
            ));
        }
    }
    if observed_disposition != AssessmentDisposition::Refused.code() {
        let required = claim
            .required_evidence()
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        if evidence_sources != &required {
            return Err(malformed_assessment_log(
                "assessment log non-refused disposition does not retain the exact required evidence role set",
            ));
        }
        if prerequisite_rows != expected_prerequisite_counts {
            return Err(malformed_assessment_log(
                "assessment log non-refused disposition does not retain the exact incoming frozen-DAG prerequisite set",
            ));
        }
    }
    Ok(())
}

fn observed_assessment_log_authority_class(
    requirement: EvidenceRequirement,
    reasons: &[String],
) -> Result<EvidenceAuthorityClass, ContractError> {
    let expected = requirement.authority_class();
    let prefix = format!(
        "weak-authority:{}:requires-{}:observed-",
        requirement.code(),
        expected.code()
    );
    let mut observed = None;
    for reason in reasons {
        let Some(code) = reason.strip_prefix(&prefix) else {
            continue;
        };
        let class = assessment_log_authority_class(code).ok_or_else(|| {
            malformed_assessment_log(
                "assessment log weak-authority reason carries an unknown observed class",
            )
        })?;
        if class == expected || observed.replace(class).is_some() {
            return Err(malformed_assessment_log(
                "assessment log carries an impossible weak-authority class binding",
            ));
        }
    }
    Ok(observed.unwrap_or(expected))
}

fn assessment_log_requirement_reason(reason: &str, prefix: &str) -> Option<EvidenceRequirement> {
    let suffix = reason.strip_prefix(prefix)?;
    (!suffix.contains(':'))
        .then(|| assessment_log_requirement(suffix))
        .flatten()
}

fn assessment_log_weakness_requirement(reason: &str) -> Option<EvidenceRequirement> {
    let requirement = if let Some(remainder) = reason.strip_prefix("weak-authority:") {
        remainder.split_once(':')?.0
    } else if let Some(remainder) = reason.strip_prefix("weak-validity-domain:") {
        remainder.strip_suffix(":does-not-cover-case")?
    } else {
        let remainder = reason.strip_prefix("weak-independence:")?;
        remainder.strip_suffix(":independent-evidence-required")?
    };
    assessment_log_requirement(requirement)
}

fn assessment_log_access_class(code: &str) -> bool {
    matches!(
        code,
        "not-applicable" | "calibration" | "validation" | "blind-holdout"
    )
}

fn expected_assessment_log_access_class(requirement: EvidenceRequirement) -> &'static str {
    match requirement {
        EvidenceRequirement::CalibrationPartition => "calibration",
        EvidenceRequirement::PhysicalValidation
        | EvidenceRequirement::RivalMechanismDiscrimination => "validation",
        EvidenceRequirement::BlindHoldout => "blind-holdout",
        _ => "not-applicable",
    }
}

fn assessment_log_authority_class(code: &str) -> Option<EvidenceAuthorityClass> {
    match code {
        "structural-process" => Some(EvidenceAuthorityClass::StructuralProcess),
        "verified-numerics" => Some(EvidenceAuthorityClass::VerifiedNumerics),
        "validated-physical" => Some(EvidenceAuthorityClass::ValidatedPhysical),
        _ => None,
    }
}

fn assessment_log_evidence_use(code: &str) -> Option<EvidenceUse> {
    match code {
        "calibration-input" => Some(EvidenceUse::CalibrationInput),
        "validation-input" => Some(EvidenceUse::ValidationInput),
        _ => None,
    }
}

fn expected_assessment_log_units(claim: EulerClaimKind) -> &'static [&'static str] {
    match claim {
        EulerClaimKind::NumericalTrajectoryVerification
        | EulerClaimKind::CalibratedReproduction
        | EulerClaimKind::BlindTrajectoryPrediction
        | EulerClaimKind::QualitativeEffectDirection
        | EulerClaimKind::Ranking
        | EulerClaimKind::MechanismAttribution => &["1"],
        EulerClaimKind::EventOrCrossoverPrediction => &["1", "s"],
        EulerClaimKind::NonlinearOptimumInterval => &["1", "mm"],
        EulerClaimKind::EnergyChannelAttribution => &["1", "j"],
    }
}

#[allow(clippy::too_many_lines)] // Closed v1 reason grammar is audited in one exhaustive matcher.
fn is_closed_assessment_log_reason(
    reason: &str,
    claim: EulerClaimKind,
    expected_disposition: &str,
    observed_disposition: &str,
    units: &[String],
    evidence_sources: &BTreeSet<EvidenceRequirement>,
) -> bool {
    if matches!(
        reason,
        "contract-identity-mismatch"
            | "binding-no-claims-not-accepted"
            | "protected-target-fitting-invalidates-emergent-claim"
    ) {
        return true;
    }
    let expected_mismatch = format!(
        "expected-disposition-mismatch:expected-{expected_disposition}:observed-{observed_disposition}"
    );
    if reason == expected_mismatch {
        return true;
    }
    for prefix in ["out-of-domain-numeric:", "missing-numeric-axis:"] {
        if let Some(axis) = reason.strip_prefix(prefix) {
            return matches!(
                axis,
                "outer-radius"
                    | "body-thickness"
                    | "total-mass"
                    | "edge-radius"
                    | "initial-angular-speed"
                    | "ambient-pressure"
                    | "ambient-temperature"
                    | "base-slope"
            );
        }
    }
    for prefix in ["out-of-domain-category:", "missing-categorical-axis:"] {
        if let Some(axis) = reason.strip_prefix(prefix) {
            return matches!(
                axis,
                "disc-material-family"
                    | "mass-distribution-family"
                    | "base-material-family"
                    | "support-assembly"
                    | "environment-regime"
                    | "contact-regime"
                    | "observation-frame"
            );
        }
    }
    if let Some(remainder) = reason.strip_prefix("claim-unit-set-mismatch:expected-") {
        let Some((expected, observed)) = remainder.split_once(":observed-") else {
            return false;
        };
        return expected == expected_assessment_log_units(claim).join(",")
            && observed == units.join(",")
            && observed != expected;
    }
    for prefix in [
        "missing-evidence:",
        "stale-evidence-binding:",
        "qoi-binding-mismatch:",
        "unexpected-evidence:",
    ] {
        if let Some(requirement) = assessment_log_requirement_reason(reason, prefix) {
            let required = claim.required_evidence().contains(&requirement);
            return match prefix {
                "missing-evidence:" => required && !evidence_sources.contains(&requirement),
                "unexpected-evidence:" => evidence_sources.contains(&requirement) && !required,
                _ => required && evidence_sources.contains(&requirement),
            };
        }
    }
    if let Some(remainder) = reason.strip_prefix("hypothesis-source-cannot-satisfy-evidence:") {
        let Some((requirement, slot)) = remainder.split_once(':') else {
            return false;
        };
        return assessment_log_requirement(requirement)
            .is_some_and(|requirement| evidence_sources.contains(&requirement))
            && assessment_log_evidence_slot(slot).is_some();
    }
    if let Some(role) = reason.strip_prefix("hypothesis-source-cannot-satisfy-packet-role:") {
        return matches!(role, "design-set" | "aggregate-qoi-derivation-receipt");
    }
    if let Some(remainder) = reason.strip_prefix("access-class-mismatch:") {
        let mut parts = remainder.split(':');
        let (Some(requirement), Some(expected), Some(observed), None) =
            (parts.next(), parts.next(), parts.next(), parts.next())
        else {
            return false;
        };
        let Some(requirement) = assessment_log_requirement(requirement) else {
            return false;
        };
        let (Some(expected), Some(observed)) = (
            expected.strip_prefix("expected-"),
            observed.strip_prefix("observed-"),
        ) else {
            return false;
        };
        return evidence_sources.contains(&requirement)
            && claim.required_evidence().contains(&requirement)
            && expected == expected_assessment_log_access_class(requirement)
            && assessment_log_access_class(observed)
            && observed != expected;
    }
    if let Some(remainder) = reason.strip_prefix("weak-authority:") {
        let Some((requirement, weakness)) = remainder.split_once(':') else {
            return false;
        };
        let Some(requirement) = assessment_log_requirement(requirement) else {
            return false;
        };
        if !evidence_sources.contains(&requirement) {
            return false;
        }
        if matches!(
            weakness,
            "requires-finite-verified-color" | "verified-enclosure-is-vacuous"
        ) {
            return requirement.authority_class() == EvidenceAuthorityClass::VerifiedNumerics;
        }
        if weakness == "requires-validated-color" {
            return requirement.authority_class() == EvidenceAuthorityClass::ValidatedPhysical;
        }
        let Some((expected, observed)) = weakness.split_once(":observed-") else {
            return false;
        };
        let Some(expected) = expected.strip_prefix("requires-") else {
            return false;
        };
        return expected == requirement.authority_class().code()
            && assessment_log_authority_class(observed).is_some()
            && observed != expected;
    }
    if let Some(remainder) = reason.strip_prefix("weak-validity-domain:") {
        let Some(requirement) = remainder.strip_suffix(":does-not-cover-case") else {
            return false;
        };
        return assessment_log_requirement(requirement).is_some_and(|requirement| {
            evidence_sources.contains(&requirement)
                && requirement.authority_class() == EvidenceAuthorityClass::ValidatedPhysical
        });
    }
    if let Some(remainder) = reason.strip_prefix("weak-independence:") {
        let Some(requirement) = remainder.strip_suffix(":independent-evidence-required") else {
            return false;
        };
        return assessment_log_requirement(requirement).is_some_and(|requirement| {
            evidence_sources.contains(&requirement)
                && matches!(
                    requirement,
                    EvidenceRequirement::IndependentReconstruction
                        | EvidenceRequirement::PhysicalValidation
                        | EvidenceRequirement::BlindHoldout
                        | EvidenceRequirement::RivalMechanismDiscrimination
                )
        });
    }
    if let Some(code) = reason.strip_prefix("malformed-prerequisite-receipt:") {
        return code == "EulerProtocolPrerequisiteReceiptIdentityMismatch";
    }
    for prefix in [
        "stale-prerequisite-receipt:",
        "prerequisite-applicability-point-mismatch:",
        "prerequisite-design-set-mismatch:",
    ] {
        if let Some(prerequisite) = reason.strip_prefix(prefix) {
            return assessment_log_claim_kind(prerequisite).is_some();
        }
    }
    for prefix in [
        "unexpected-prerequisite-receipt:",
        "duplicate-prerequisite-receipt:",
        "missing-prerequisite-receipt:",
    ] {
        if let Some(remainder) = reason.strip_prefix(prefix) {
            let Some((prerequisite, use_kind)) = remainder.split_once(':') else {
                return false;
            };
            return assessment_log_claim_kind(prerequisite).is_some()
                && assessment_log_evidence_use(use_kind).is_some();
        }
    }
    false
}

#[allow(clippy::too_many_lines)] // The exact field sequence is the v1 transport contract.
fn validate_claim_policy_assessment_log_json_line(json_line: &str) -> Result<(), ContractError> {
    if json_line.is_empty()
        || json_line.len() > MAX_ASSESSMENT_LOG_BYTES
        || !json_line.ends_with('\n')
        || json_line[..json_line.len().saturating_sub(1)].contains('\n')
        || json_line[..json_line.len().saturating_sub(1)].contains('\r')
    {
        return Err(malformed_assessment_log(
            "assessment log must be one bounded canonical JSON line ending in LF",
        ));
    }
    let line = &json_line[..json_line.len() - 1];
    let mut reader = AssessmentLogReader::new(line);
    reader.expect("{\"schema_version\":1,\"protocol_id\":", "schema_version")?;
    let protocol_id = reader.string("protocol_id", MAX_EULER_TEXT_BYTES)?;
    if protocol_id != CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN {
        return Err(malformed_assessment_log(
            "assessment log protocol_id is not the exact v1 domain",
        ));
    }
    reader.expect(",\"contract_identity\":", "contract_identity")?;
    let contract_identity = reader.string("contract_identity", 64)?;
    reader.expect(",\"packet_contract_identity\":", "packet_contract_identity")?;
    let packet_contract_identity = reader.string("packet_contract_identity", 64)?;
    reader.expect(",\"packet_identity\":", "packet_identity")?;
    let packet_identity = reader.string("packet_identity", 64)?;
    reader.expect(",\"design_set_identity\":", "design_set_identity")?;
    let design_set_identity = reader.string("design_set_identity", 64)?;
    reader.expect(
        ",\"aggregate_qoi_derivation_receipt_identity\":",
        "aggregate_qoi_derivation_receipt_identity",
    )?;
    let aggregate_qoi_derivation_receipt_identity =
        reader.string("aggregate_qoi_derivation_receipt_identity", 64)?;
    for (field, value) in [
        ("contract_identity", contract_identity.as_str()),
        (
            "packet_contract_identity",
            packet_contract_identity.as_str(),
        ),
        ("packet_identity", packet_identity.as_str()),
        ("design_set_identity", design_set_identity.as_str()),
        (
            "aggregate_qoi_derivation_receipt_identity",
            aggregate_qoi_derivation_receipt_identity.as_str(),
        ),
    ] {
        if !is_exact_lower_hex(value, 64) || value.bytes().all(|byte| byte == b'0') {
            return Err(malformed_assessment_log(format!(
                "assessment log {field} is not a nonzero lowercase 256-bit identity"
            )));
        }
    }
    if contract_identity != FROZEN_CONTRACT_IDENTITY_HEX {
        return Err(malformed_assessment_log(
            "assessment log contract_identity is not the exact frozen v1 contract",
        ));
    }
    let frozen_contract = build_frozen_contract().map_err(|_| {
        malformed_assessment_log("the local frozen v1 contract cannot be reconstructed")
    })?;
    if frozen_contract.identity().as_hash().to_hex() != FROZEN_CONTRACT_IDENTITY_HEX {
        return Err(malformed_assessment_log(
            "the reconstructed frozen v1 contract does not match its literal identity",
        ));
    }
    reader.expect(",\"case_id\":", "case_id")?;
    let case_id = reader.string("case_id", MAX_PROTOCOL_ID_BYTES)?;
    if validate_protocol_id("assessment_log.case_id", &case_id).is_err() {
        return Err(malformed_assessment_log(
            "assessment log case_id is not a canonical machine identity",
        ));
    }
    reader.expect(",\"claim\":", "claim")?;
    let claim = reader.string("claim", MAX_PROTOCOL_ID_BYTES)?;
    let claim_kind = assessment_log_claim_kind(&claim).ok_or_else(|| {
        malformed_assessment_log("assessment log claim is not a closed v1 claim identifier")
    })?;
    reader.expect(",\"packet_source_id\":", "packet_source_id")?;
    let packet_source_id = reader.string("packet_source_id", MAX_PROTOCOL_ID_BYTES)?;
    if packet_source_id != case_id {
        return Err(malformed_assessment_log(
            "assessment log packet_source_id does not equal its exact case_id",
        ));
    }
    reader.expect(",\"packet_source_schema\":", "packet_source_schema")?;
    let packet_source_schema = reader.string("packet_source_schema", MAX_EULER_TEXT_BYTES)?;
    if packet_source_schema != PACKET_SCHEMA {
        return Err(malformed_assessment_log(
            "assessment log packet_source_schema is not the exact packet domain",
        ));
    }
    reader.expect(",\"evidence_sources\":", "evidence_sources")?;
    let evidence_sources = reader.string_array(
        "evidence_sources",
        MAX_EVIDENCE_RECORDS,
        MAX_ASSESSMENT_LOG_EVIDENCE_SOURCE_BYTES,
    )?;
    if evidence_sources.iter().any(String::is_empty) || !is_strictly_sorted(&evidence_sources) {
        return Err(malformed_assessment_log(
            "assessment log evidence_sources must be unique and canonically ordered",
        ));
    }
    let evidence_requirements = validate_assessment_log_evidence_sources(&evidence_sources)?;
    reader.expect(",\"units\":", "units")?;
    let units =
        reader.string_array("units", MAX_ASSESSMENT_LOG_UNIT_ROWS, MAX_PROTOCOL_ID_BYTES)?;
    if units.is_empty()
        || !is_strictly_sorted(&units)
        || units
            .iter()
            .any(|unit| validate_protocol_id("assessment_log.unit", unit).is_err())
    {
        return Err(malformed_assessment_log(
            "assessment log units must be nonempty unique canonical machine identities",
        ));
    }
    reader.expect(",\"seed\":{\"kind\":", "seed")?;
    let seed_kind = reader.string("seed.kind", 32)?;
    match seed_kind.as_str() {
        "fixed" => {
            reader.expect(",\"value\":", "seed.value")?;
            let _ = reader.unsigned("seed.value")?;
        }
        "not-applicable" => {
            reader.expect(",\"reason\":", "seed.reason")?;
            let reason = reader.string("seed.reason", MAX_PROTOCOL_ID_BYTES)?;
            if validate_protocol_id("assessment_log.seed.reason", &reason).is_err() {
                return Err(malformed_assessment_log(
                    "assessment log seed.reason is not a canonical machine identity",
                ));
            }
        }
        _ => {
            return Err(malformed_assessment_log(
                "assessment log seed.kind is not a closed v1 variant",
            ));
        }
    }
    reader.expect("},\"budgets\":{\"max_wall_time_ms\":", "budgets")?;
    let max_wall_time_ms = reader.unsigned("budgets.max_wall_time_ms")?;
    reader.expect(",\"max_memory_bytes\":", "budgets.max_memory_bytes")?;
    let max_memory_bytes = reader.unsigned("budgets.max_memory_bytes")?;
    reader.expect(
        ",\"normalized_accuracy_limit_bits\":",
        "budgets.normalized_accuracy_limit_bits",
    )?;
    let normalized_accuracy_bits = reader.string("budgets.normalized_accuracy_limit_bits", 16)?;
    if max_wall_time_ms == 0
        || max_memory_bytes == 0
        || !is_exact_lower_hex(&normalized_accuracy_bits, 16)
    {
        return Err(malformed_assessment_log(
            "assessment log budgets are not the exact positive/canonical v1 representation",
        ));
    }
    let bits = u64::from_str_radix(&normalized_accuracy_bits, 16).map_err(|_| {
        malformed_assessment_log("assessment log normalized_accuracy_limit_bits is invalid")
    })?;
    let normalized_accuracy_limit = f64::from_bits(bits);
    if !normalized_accuracy_limit.is_finite()
        || normalized_accuracy_limit < 0.0
        || (normalized_accuracy_limit == 0.0 && bits != 0)
    {
        return Err(malformed_assessment_log(
            "assessment log normalized_accuracy_limit_bits is nonfinite, negative, or noncanonical zero",
        ));
    }
    reader.expect("},\"expected_disposition\":", "expected_disposition")?;
    let expected_disposition = reader.string("expected_disposition", 64)?;
    reader.expect(",\"observed_disposition\":", "observed_disposition")?;
    let observed_disposition = reader.string("observed_disposition", 64)?;
    let disposition_codes = [
        AssessmentDisposition::ReferenceCompleteCandidate.code(),
        AssessmentDisposition::DemotedCandidate.code(),
        AssessmentDisposition::RetainedTerminal.code(),
        AssessmentDisposition::Refused.code(),
    ];
    if !disposition_codes.contains(&expected_disposition.as_str())
        || !disposition_codes.contains(&observed_disposition.as_str())
    {
        return Err(malformed_assessment_log(
            "assessment log carries an unknown disposition",
        ));
    }
    reader.expect(
        ",\"reported_scientific_disposition\":",
        "reported_scientific_disposition",
    )?;
    let reported = reader.string("reported_scientific_disposition", 32)?;
    if !["positive", "negative", "inconclusive"].contains(&reported.as_str()) {
        return Err(malformed_assessment_log(
            "assessment log carries an unknown reported scientific disposition",
        ));
    }
    reader.expect(",\"target_fitted\":", "target_fitted")?;
    let target_fitted = reader.boolean("target_fitted")?;
    reader.expect(",\"applicability_state\":", "applicability_state")?;
    if reader.string("applicability_state", 128)?
        != "campaign-anchor-point-plus-content-bound-design-set"
    {
        return Err(malformed_assessment_log(
            "assessment log applicability_state is not the exact v1 value",
        ));
    }
    reader.expect(
        ",\"criterion_evaluation_state\":",
        "criterion_evaluation_state",
    )?;
    if reader.string("criterion_evaluation_state", 128)?
        != "aggregate-qoi-derivation-receipt-unreadmitted-evaluation-deferred-to-preregistered-analysis"
    {
        return Err(malformed_assessment_log(
            "assessment log criterion_evaluation_state is not the exact v1 value",
        ));
    }
    reader.expect(",\"first_divergence\":", "first_divergence")?;
    let first_divergence = if reader.input[reader.offset..].starts_with("null") {
        reader.offset += 4;
        None
    } else {
        Some(reader.string("first_divergence", MAX_ASSESSMENT_LOG_REASON_BYTES)?)
    };
    reader.expect(",\"reasons\":", "reasons")?;
    let reasons = reader.string_array(
        "reasons",
        MAX_ASSESSMENT_LOG_REASON_ROWS,
        MAX_ASSESSMENT_LOG_REASON_BYTES,
    )?;
    if reasons.iter().any(String::is_empty) || !is_strictly_sorted(&reasons) {
        return Err(malformed_assessment_log(
            "assessment log reasons must be unique and canonically ordered",
        ));
    }
    if reasons.iter().any(|reason| {
        !is_closed_assessment_log_reason(
            reason,
            claim_kind,
            &expected_disposition,
            &observed_disposition,
            &units,
            &evidence_requirements,
        )
    }) {
        return Err(malformed_assessment_log(
            "assessment log carries a reason outside the closed v1 reason grammar",
        ));
    }
    let contract_mismatch_reason = "contract-identity-mismatch".to_owned();
    let packet_contract_mismatch = packet_contract_identity != contract_identity;
    let contract_mismatch_reason_present = reasons.binary_search(&contract_mismatch_reason).is_ok();
    if packet_contract_mismatch != contract_mismatch_reason_present
        || (packet_contract_mismatch
            && observed_disposition != AssessmentDisposition::Refused.code())
    {
        return Err(malformed_assessment_log(
            "assessment log packet contract identity is not bound bidirectionally to its exact mismatch reason and refusal",
        ));
    }
    // The evaluator visits each context axis once and can observe either an
    // absent value or a present value outside the domain, never both. It also
    // retains at most one evidence record per role, so one role cannot carry
    // two different observed access classes. The reason grammar validates
    // each row in isolation; enforce these writer-level exclusivity rules
    // separately so a hostile line cannot splice individually valid but
    // mutually contradictory diagnostics together.
    let mut point_diagnostic_axes = BTreeSet::<(&str, &str)>::new();
    let mut access_diagnostic_roles = BTreeSet::<EvidenceRequirement>::new();
    for reason in &reasons {
        let point_diagnostic = reason
            .strip_prefix("out-of-domain-numeric:")
            .map(|axis| ("numeric", axis))
            .or_else(|| {
                reason
                    .strip_prefix("missing-numeric-axis:")
                    .map(|axis| ("numeric", axis))
            })
            .or_else(|| {
                reason
                    .strip_prefix("out-of-domain-category:")
                    .map(|axis| ("categorical", axis))
            })
            .or_else(|| {
                reason
                    .strip_prefix("missing-categorical-axis:")
                    .map(|axis| ("categorical", axis))
            });
        if point_diagnostic.is_some_and(|diagnostic| !point_diagnostic_axes.insert(diagnostic)) {
            return Err(malformed_assessment_log(
                "assessment log carries mutually exclusive point diagnostics for one context axis",
            ));
        }

        if let Some(requirement) = reason
            .strip_prefix("access-class-mismatch:")
            .and_then(|remainder| remainder.split_once(':').map(|(role, _)| role))
            .and_then(assessment_log_requirement)
            && !access_diagnostic_roles.insert(requirement)
        {
            return Err(malformed_assessment_log(
                "assessment log carries multiple observed access classes for one evidence role",
            ));
        }
    }
    let expected_unit_set = expected_assessment_log_units(claim_kind);
    let expected_units = expected_unit_set.join(",");
    let observed_units = units.join(",");
    let unit_mismatch_reason =
        format!("claim-unit-set-mismatch:expected-{expected_units}:observed-{observed_units}");
    let unit_mismatch = units.len() != expected_unit_set.len()
        || units
            .iter()
            .map(String::as_str)
            .ne(expected_unit_set.iter().copied());
    let unit_mismatch_retained = reasons.binary_search(&unit_mismatch_reason).is_ok();
    if unit_mismatch != unit_mismatch_retained
        || (unit_mismatch && observed_disposition != AssessmentDisposition::Refused.code())
    {
        return Err(malformed_assessment_log(
            "assessment log unit state is not bound bidirectionally to its exact mismatch reason and refusal",
        ));
    }
    for requirement in EULER_EVIDENCE_REQUIREMENT_REGISTRY {
        let required = claim_kind.required_evidence().contains(&requirement);
        let observed = evidence_requirements.contains(&requirement);
        let missing_reason = format!("missing-evidence:{}", requirement.code());
        let unexpected_reason = format!("unexpected-evidence:{}", requirement.code());
        let missing_retained = reasons.binary_search(&missing_reason).is_ok();
        let unexpected_retained = reasons.binary_search(&unexpected_reason).is_ok();
        if (required && !observed) != missing_retained
            || (!required && observed) != unexpected_retained
            || ((missing_retained || unexpected_retained)
                && observed_disposition != AssessmentDisposition::Refused.code())
        {
            return Err(malformed_assessment_log(
                "assessment log evidence-role presence is not bound bidirectionally to its missing/unexpected reason and refusal",
            ));
        }
    }
    match first_divergence.as_ref() {
        None if reasons.is_empty() => {}
        Some(first) if reasons.binary_search(first).is_ok() => {}
        _ => {
            return Err(malformed_assessment_log(
                "assessment log first_divergence must be exactly one retained reason",
            ));
        }
    }
    let expected_mismatch_reason = format!(
        "expected-disposition-mismatch:expected-{expected_disposition}:observed-{observed_disposition}"
    );
    if reasons.iter().any(|reason| {
        reason.starts_with("expected-disposition-mismatch:") && reason != &expected_mismatch_reason
    }) {
        return Err(malformed_assessment_log(
            "assessment log carries a noncanonical expected-disposition mismatch reason",
        ));
    }
    let expected_mismatch_present = reasons.binary_search(&expected_mismatch_reason).is_ok();
    if (expected_disposition != observed_disposition) != expected_mismatch_present {
        return Err(malformed_assessment_log(
            "assessment log expected/observed disposition mismatch is not bound to its exact reason",
        ));
    }
    let target_fitting_reason = "protected-target-fitting-invalidates-emergent-claim".to_owned();
    let target_fitting_present = reasons.binary_search(&target_fitting_reason).is_ok();
    let target_fitting_required = target_fitted && claim_kind.forbids_target_fitting();
    if target_fitting_required != target_fitting_present
        || (target_fitting_required && observed_disposition != "refused")
    {
        return Err(malformed_assessment_log(
            "assessment log target-fitting state is not bound to its refusal reason/disposition",
        ));
    }
    if (matches!(
        observed_disposition.as_str(),
        "reference-complete-candidate-unreadmitted" | "demoted-candidate"
    ) && reported != ReportedScientificDisposition::Positive.code())
        || (observed_disposition == AssessmentDisposition::RetainedTerminal.code()
            && reported == ReportedScientificDisposition::Positive.code())
    {
        return Err(malformed_assessment_log(
            "assessment log reported outcome is inconsistent with its candidate/terminal disposition",
        ));
    }
    let policy_reasons = reasons
        .iter()
        .filter(|reason| *reason != &expected_mismatch_reason)
        .collect::<Vec<_>>();
    let disposition_reason_shape_is_valid = match observed_disposition.as_str() {
        "reference-complete-candidate-unreadmitted" => policy_reasons.is_empty(),
        "demoted-candidate" => {
            !policy_reasons.is_empty()
                && policy_reasons
                    .iter()
                    .all(|reason| reason.starts_with("weak-"))
        }
        "retained-terminal-non-promotion" => policy_reasons
            .iter()
            .all(|reason| reason.starts_with("weak-")),
        "refused" => {
            !policy_reasons.is_empty()
                && policy_reasons
                    .iter()
                    .any(|reason| !reason.starts_with("weak-"))
        }
        _ => false,
    };
    if !disposition_reason_shape_is_valid {
        return Err(malformed_assessment_log(
            "assessment log reason shape is impossible for its observed disposition",
        ));
    }
    reader.expect(",\"authority_state\":", "authority_state")?;
    let authority_state = reader.string("authority_state", 64)?;
    let expected_authority_state = match observed_disposition.as_str() {
        "reference-complete-candidate-unreadmitted" => "unreadmitted-reference-candidate-only",
        "demoted-candidate" => "demoted-below-requested-claim",
        "retained-terminal-non-promotion" => "terminal-non-promotion",
        "refused" => "no-candidate-authority",
        _ => {
            return Err(malformed_assessment_log(
                "assessment log carries an unknown observed disposition",
            ));
        }
    };
    if authority_state != expected_authority_state {
        return Err(malformed_assessment_log(
            "assessment log authority_state does not match observed_disposition",
        ));
    }
    reader.expect(",\"no_claim_state\":", "no_claim_state")?;
    let no_claim_state = reader.string("no_claim_state", 32)?;
    if !["accepted", "not-accepted"].contains(&no_claim_state.as_str()) {
        return Err(malformed_assessment_log(
            "assessment log no_claim_state is not a closed v1 value",
        ));
    }
    let no_claim_reason = "binding-no-claims-not-accepted".to_owned();
    let no_claim_reason_present = reasons.binary_search(&no_claim_reason).is_ok();
    let no_claims_not_accepted = no_claim_state == "not-accepted";
    if no_claims_not_accepted != no_claim_reason_present
        || (no_claims_not_accepted && observed_disposition != "refused")
    {
        return Err(malformed_assessment_log(
            "assessment log no-claim state is not bound to its refusal reason/disposition",
        ));
    }
    reader.expect(",\"relative_artifacts\":", "relative_artifacts")?;
    let relative_artifacts = reader.string_array(
        "relative_artifacts",
        MAX_ASSESSMENT_LOG_ARTIFACT_ROWS,
        MAX_EULER_TEXT_BYTES,
    )?;
    if relative_artifacts.is_empty()
        || relative_artifacts.iter().any(String::is_empty)
        || !is_strictly_sorted(&relative_artifacts)
    {
        return Err(malformed_assessment_log(
            "assessment log relative_artifacts must be nonempty, unique, and ordered",
        ));
    }
    validate_assessment_log_relative_artifacts(
        &relative_artifacts,
        &packet_identity,
        &design_set_identity,
        &aggregate_qoi_derivation_receipt_identity,
        claim_kind,
        &observed_disposition,
        &evidence_requirements,
        &reasons,
        &frozen_contract,
    )?;
    reader.expect(",\"reproduction_command\":", "reproduction_command")?;
    if reader.string("reproduction_command", MAX_EULER_TEXT_BYTES)? != REPRODUCTION_COMMAND {
        return Err(malformed_assessment_log(
            "assessment log reproduction_command is not the exact focused checker smoke",
        ));
    }
    reader.expect(
        ",\"artifact_resolution_state\":",
        "artifact_resolution_state",
    )?;
    if reader.string("artifact_resolution_state", 128)?
        != "logical-content-identities-only-not-persisted-by-this-crate"
    {
        return Err(malformed_assessment_log(
            "assessment log artifact_resolution_state is not the exact v1 boundary",
        ));
    }
    reader.expect(",\"redaction\":", "redaction")?;
    if reader.string("redaction", 128)?
        != "bounded-structured-protocol-metadata-no-raw-payload-or-artifact-bytes"
    {
        return Err(malformed_assessment_log(
            "assessment log redaction value is not the exact v1 boundary",
        ));
    }
    reader.expect("}", "terminal object")?;
    reader.finish()
}

impl ClaimPolicyAssessmentLog {
    /// Admit one bounded canonical v1 assessment log and derive its exact byte
    /// identity. The strict co-versioned reader enforces the complete field
    /// order, JSON types, closed values, canonical escaping, bounds, and
    /// locally checkable cross-field bindings; it is not a general JSON parser.
    /// It cannot re-evaluate or authenticate the absent packet, prerequisite
    /// receipts, referenced evidence artifacts, producers, or scientific claim.
    pub fn from_json_line(json_line: impl Into<String>) -> Result<Self, ContractError> {
        let json_line = json_line.into();
        validate_claim_policy_assessment_log_json_line(&json_line)?;
        let identity = claim_policy_assessment_log_identity(&json_line);
        Ok(Self {
            json_line,
            identity,
        })
    }

    #[must_use]
    /// Complete deterministic redacted JSON-lines record.
    pub fn json_line(&self) -> &str {
        &self.json_line
    }

    #[must_use]
    /// Domain-separated hash of the exact JSON line.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Recompute the retained-line identity and enforce the output contract.
    pub fn verify_identity(&self) -> Result<(), ContractError> {
        validate_claim_policy_assessment_log_json_line(&self.json_line)?;
        let expected = claim_policy_assessment_log_identity(&self.json_line);
        if expected != self.identity {
            return Err(ContractError::new(
                "EulerProtocolAssessmentLogIdentityMismatch",
                "assessment log identity does not match its exact JSON bytes",
            ));
        }
        Ok(())
    }
}

fn json_string(output: &mut String, value: &str) {
    output.push('"');
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            control if control.is_control() => {
                use core::fmt::Write as _;
                let _ = write!(output, "\\u{:04x}", control as u32);
            }
            other => output.push(other),
        }
    }
    output.push('"');
}

fn push_string_array(output: &mut String, values: &[String]) {
    output.push('[');
    for (index, value) in values.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        json_string(output, value);
    }
    output.push(']');
}

#[allow(clippy::too_many_lines)] // Exact JSON field order is part of the retained identity.
fn claim_policy_assessment_log(
    contract: &EulerScientificContract,
    packet: &ClaimEvidencePacket,
    prerequisites: &[PrerequisiteAssessmentReceipt],
    disposition: AssessmentDisposition,
    reasons: &[String],
    first_divergence: Option<&str>,
) -> Result<ClaimPolicyAssessmentLog, ContractError> {
    use core::fmt::Write as _;
    let prerequisite_route = contract
        .owner_matrix()
        .rows()
        .get(&OwnerRole::PrerequisiteAssessmentReceipt)
        .ok_or_else(|| {
            ContractError::new(
                "EulerProtocolOwnerRouteMissing",
                "owner matrix is missing the prerequisite-receipt route",
            )
        })?;
    let derivation_route = contract
        .owner_matrix()
        .rows()
        .get(&OwnerRole::AggregateQoiDerivationReceipt)
        .ok_or_else(|| {
            ContractError::new(
                "EulerProtocolOwnerRouteMissing",
                "owner matrix is missing the aggregate-QoI derivation-receipt route",
            )
        })?;
    let mut artifacts = packet
        .records
        .values()
        .flat_map(|record| {
            evidence_hash_references(record)
                .into_iter()
                .filter_map(move |(slot, reference)| {
                    reference.map(|hash| {
                        format!(
                            "evidence:{}:{slot}:{}",
                            record.requirement.code(),
                            hash.to_hex()
                        )
                    })
                })
        })
        .collect::<Vec<_>>();
    artifacts.push(format!("packet:{}", packet.identity.to_hex()));
    artifacts.push(format!(
        "design-set:{}",
        packet.design_set_identity.to_hex()
    ));
    artifacts.push(format!(
        "aggregate-qoi-derivation:{}:{}",
        derivation_route.source_schema(),
        packet.aggregate_qoi_derivation_receipt_identity.to_hex()
    ));
    artifacts.extend(prerequisites.iter().map(|receipt| {
        format!(
            "prerequisite:{}:{}:{}",
            receipt.prerequisite.id(),
            prerequisite_route.source_schema(),
            receipt.identity.to_hex()
        )
    }));
    artifacts.sort();
    artifacts.dedup();
    let mut evidence_sources = packet
        .records
        .values()
        .map(|record| {
            format!(
                "{}:{}:{}:{}",
                record.requirement.code(),
                record.source_kind.slug(),
                record.source_schema,
                record.source_id
            )
        })
        .collect::<Vec<_>>();
    evidence_sources.sort();
    evidence_sources.dedup();
    let authority_state = match disposition {
        AssessmentDisposition::ReferenceCompleteCandidate => {
            "unreadmitted-reference-candidate-only"
        }
        AssessmentDisposition::DemotedCandidate => "demoted-below-requested-claim",
        AssessmentDisposition::RetainedTerminal => "terminal-non-promotion",
        AssessmentDisposition::Refused => "no-candidate-authority",
    };
    let no_claim_state = if packet.no_claims_accepted {
        "accepted"
    } else {
        "not-accepted"
    };
    let mut json = String::with_capacity(4_096);
    let _ = write!(
        json,
        "{{\"schema_version\":{EULER_PROTOCOL_SCHEMA_VERSION},\"protocol_id\":"
    );
    json_string(&mut json, CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN);
    json.push_str(",\"contract_identity\":");
    json_string(&mut json, &contract.identity().as_hash().to_hex());
    json.push_str(",\"packet_contract_identity\":");
    json_string(&mut json, &packet.contract_identity.as_hash().to_hex());
    json.push_str(",\"packet_identity\":");
    json_string(&mut json, &packet.identity.to_hex());
    json.push_str(",\"design_set_identity\":");
    json_string(&mut json, &packet.design_set_identity.to_hex());
    json.push_str(",\"aggregate_qoi_derivation_receipt_identity\":");
    json_string(
        &mut json,
        &packet.aggregate_qoi_derivation_receipt_identity.to_hex(),
    );
    json.push_str(",\"case_id\":");
    json_string(&mut json, &packet.case_id);
    json.push_str(",\"claim\":");
    json_string(&mut json, packet.claim.id());
    json.push_str(",\"packet_source_id\":");
    json_string(&mut json, &packet.case_id);
    json.push_str(",\"packet_source_schema\":");
    json_string(&mut json, PACKET_SCHEMA);
    json.push_str(",\"evidence_sources\":");
    push_string_array(&mut json, &evidence_sources);
    json.push_str(",\"units\":");
    push_string_array(&mut json, &packet.units);
    json.push_str(",\"seed\":");
    match &packet.seed {
        ProtocolSeed::Fixed { value: seed } => {
            let _ = write!(json, "{{\"kind\":\"fixed\",\"value\":{seed}}}");
        }
        ProtocolSeed::NotApplicable { reason } => {
            json.push_str("{\"kind\":\"not-applicable\",\"reason\":");
            json_string(&mut json, reason);
            json.push('}');
        }
    }
    let _ = write!(
        json,
        ",\"budgets\":{{\"max_wall_time_ms\":{},\"max_memory_bytes\":{},\"normalized_accuracy_limit_bits\":\"{:016x}\"}}",
        packet.budget.max_wall_time_ms,
        packet.budget.max_memory_bytes,
        packet.budget.normalized_accuracy_limit.to_bits()
    );
    json.push_str(",\"expected_disposition\":");
    json_string(&mut json, packet.expected_disposition.code());
    json.push_str(",\"observed_disposition\":");
    json_string(&mut json, disposition.code());
    json.push_str(",\"reported_scientific_disposition\":");
    json_string(&mut json, packet.reported_scientific_disposition.code());
    let _ = write!(
        json,
        ",\"target_fitted\":{},\"applicability_state\":\"campaign-anchor-point-plus-content-bound-design-set\",\"criterion_evaluation_state\":\"aggregate-qoi-derivation-receipt-unreadmitted-evaluation-deferred-to-preregistered-analysis\"",
        packet.target_fitted
    );
    json.push_str(",\"first_divergence\":");
    if let Some(divergence) = first_divergence {
        json_string(&mut json, divergence);
    } else {
        json.push_str("null");
    }
    json.push_str(",\"reasons\":");
    push_string_array(&mut json, reasons);
    json.push_str(",\"authority_state\":");
    json_string(&mut json, authority_state);
    json.push_str(",\"no_claim_state\":");
    json_string(&mut json, no_claim_state);
    json.push_str(",\"relative_artifacts\":");
    push_string_array(&mut json, &artifacts);
    json.push_str(",\"reproduction_command\":");
    json_string(&mut json, REPRODUCTION_COMMAND);
    json.push_str(",\"artifact_resolution_state\":\"logical-content-identities-only-not-persisted-by-this-crate\",\"redaction\":\"bounded-structured-protocol-metadata-no-raw-payload-or-artifact-bytes\"}\n");
    if json.len() > MAX_ASSESSMENT_LOG_BYTES {
        return Err(ContractError::new(
            "EulerProtocolLogTooLarge",
            "claim-policy assessment log exceeds its byte budget",
        ));
    }
    let identity = claim_policy_assessment_log_identity(&json);
    let log = ClaimPolicyAssessmentLog {
        json_line: json,
        identity,
    };
    log.verify_identity()?;
    Ok(log)
}

/// Policy result with exact reasons and its retained redacted log.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClaimPolicyAssessment {
    schema_version: u32,
    contract_identity: ContractIdentity,
    packet_identity: ContentHash,
    design_set_identity: ContentHash,
    aggregate_qoi_derivation_receipt_identity: ContentHash,
    point_bytes: Vec<u8>,
    case_id: String,
    claim: EulerClaimKind,
    disposition: AssessmentDisposition,
    reported_scientific_disposition: ReportedScientificDisposition,
    reasons: Vec<String>,
    log: ClaimPolicyAssessmentLog,
    identity: ContentHash,
}

#[allow(clippy::too_many_arguments)]
fn assessment_canonical_bytes(
    schema_version: u32,
    contract_identity: ContractIdentity,
    packet_identity: ContentHash,
    design_set_identity: ContentHash,
    aggregate_qoi_derivation_receipt_identity: ContentHash,
    point_bytes: &[u8],
    case_id: &str,
    claim: EulerClaimKind,
    disposition: AssessmentDisposition,
    reported_scientific_disposition: ReportedScientificDisposition,
    reasons: &[String],
    log_identity: ContentHash,
) -> Result<Vec<u8>, ContractError> {
    let mut bytes = Vec::with_capacity(512);
    bytes.extend_from_slice(ASSESSMENT_MAGIC);
    bytes.extend_from_slice(&schema_version.to_le_bytes());
    bytes.extend_from_slice(contract_identity.as_hash().as_bytes());
    bytes.extend_from_slice(packet_identity.as_bytes());
    bytes.extend_from_slice(design_set_identity.as_bytes());
    bytes.extend_from_slice(aggregate_qoi_derivation_receipt_identity.as_bytes());
    write_len(&mut bytes, point_bytes.len())?;
    bytes.extend_from_slice(point_bytes);
    write_text(&mut bytes, case_id)?;
    write_text(&mut bytes, claim.id())?;
    write_text(&mut bytes, disposition.code())?;
    write_text(&mut bytes, reported_scientific_disposition.code())?;
    write_len(&mut bytes, reasons.len())?;
    for reason in reasons {
        write_text(&mut bytes, reason)?;
    }
    bytes.extend_from_slice(log_identity.as_bytes());
    Ok(bytes)
}

impl ClaimPolicyAssessment {
    fn build(
        contract: &EulerScientificContract,
        packet: &ClaimEvidencePacket,
        prerequisites: &[PrerequisiteAssessmentReceipt],
        disposition: AssessmentDisposition,
        mut reasons: Vec<String>,
    ) -> Result<Self, ContractError> {
        let first_divergence = reasons.first().cloned();
        reasons.sort();
        reasons.dedup();
        let log = claim_policy_assessment_log(
            contract,
            packet,
            prerequisites,
            disposition,
            &reasons,
            first_divergence.as_deref(),
        )?;
        let point_bytes = applicability_point_bytes(&packet.point)?;
        let canonical = assessment_canonical_bytes(
            EULER_PROTOCOL_SCHEMA_VERSION,
            contract.identity(),
            packet.identity,
            packet.design_set_identity,
            packet.aggregate_qoi_derivation_receipt_identity,
            &point_bytes,
            &packet.case_id,
            packet.claim,
            disposition,
            packet.reported_scientific_disposition,
            &reasons,
            log.identity,
        )?;
        let identity = claim_policy_assessment_identity(&canonical);
        Ok(Self {
            schema_version: EULER_PROTOCOL_SCHEMA_VERSION,
            contract_identity: contract.identity(),
            packet_identity: packet.identity,
            design_set_identity: packet.design_set_identity,
            aggregate_qoi_derivation_receipt_identity: packet
                .aggregate_qoi_derivation_receipt_identity,
            point_bytes,
            case_id: packet.case_id.clone(),
            claim: packet.claim,
            disposition,
            reported_scientific_disposition: packet.reported_scientific_disposition,
            reasons,
            log,
            identity,
        })
    }

    /// Recompute the complete assessment and nested assessment-log identities.
    pub fn verify_identity(&self) -> Result<(), ContractError> {
        self.log.verify_identity()?;
        let expected = claim_policy_assessment_identity(&assessment_canonical_bytes(
            self.schema_version,
            self.contract_identity,
            self.packet_identity,
            self.design_set_identity,
            self.aggregate_qoi_derivation_receipt_identity,
            &self.point_bytes,
            &self.case_id,
            self.claim,
            self.disposition,
            self.reported_scientific_disposition,
            &self.reasons,
            self.log.identity,
        )?);
        if expected != self.identity {
            return Err(ContractError::new(
                "EulerProtocolAssessmentIdentityMismatch",
                "claim assessment identity does not match its semantic fields",
            ));
        }
        Ok(())
    }

    /// Bind this exact retained-positive structural assessment to a proposed
    /// direct dependent edge. The dependent assessment still checks that the
    /// edge exists in the frozen graph and uses the same exact design set and
    /// campaign-anchor applicability point.
    pub fn as_prerequisite_for(
        &self,
        dependent: EulerClaimKind,
        use_kind: EvidenceUse,
    ) -> Result<PrerequisiteAssessmentReceipt, ContractError> {
        PrerequisiteAssessmentReceipt::new(self, dependent, use_kind)
    }

    #[must_use]
    /// Contract identity assessed.
    pub const fn contract_identity(&self) -> ContractIdentity {
        self.contract_identity
    }

    #[must_use]
    /// Complete source packet identity.
    pub const fn packet_identity(&self) -> ContentHash {
        self.packet_identity
    }

    #[must_use]
    /// Exact comparison/configuration/search-set identity assessed.
    pub const fn design_set_identity(&self) -> ContentHash {
        self.design_set_identity
    }

    #[must_use]
    /// Referenced Context-bound aggregate-QoI derivation receipt.
    pub const fn aggregate_qoi_derivation_receipt_identity(&self) -> ContentHash {
        self.aggregate_qoi_derivation_receipt_identity
    }

    #[must_use]
    /// Case identity assessed.
    pub fn case_id(&self) -> &str {
        &self.case_id
    }

    #[must_use]
    /// Claim kind assessed.
    pub const fn claim(&self) -> EulerClaimKind {
        self.claim
    }

    #[must_use]
    /// Local candidate/demotion/terminal/refusal result.
    pub const fn disposition(&self) -> AssessmentDisposition {
        self.disposition
    }

    #[must_use]
    /// Caller-reported scientific outcome; no criterion evaluation is implied.
    pub const fn reported_scientific_disposition(&self) -> ReportedScientificDisposition {
        self.reported_scientific_disposition
    }

    #[must_use]
    /// Canonically ordered refusal, weakness, terminal-retention, and
    /// expected-disposition-mismatch diagnostics.
    pub fn reasons(&self) -> &[String] {
        &self.reasons
    }

    #[must_use]
    /// Bounded, redacted, content-addressed retained log.
    pub const fn log(&self) -> &ClaimPolicyAssessmentLog {
        &self.log
    }

    #[must_use]
    /// Domain-separated identity of the assessment and nested retained log.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Content-bound structural-check result. Passed means exact structural
/// equality with the frozen v1 contract, never scientific adequacy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContractCheckReceipt {
    schema_version: u32,
    checker_id: String,
    subject: ContractIdentity,
    context_hash: ContentHash,
    graph_hash: ContentHash,
    passed: bool,
    issues: Vec<String>,
    identity: ContentHash,
}

impl ContractCheckReceipt {
    fn new(
        subject: ContractIdentity,
        context_hash: ContentHash,
        graph_hash: ContentHash,
        passed: bool,
        mut issues: Vec<String>,
    ) -> Result<Self, ContractError> {
        issues.sort();
        issues.dedup();
        let checker_id = CHECKER_ID.to_owned();
        let canonical = check_receipt_bytes(
            EULER_PROTOCOL_SCHEMA_VERSION,
            &checker_id,
            subject,
            context_hash,
            graph_hash,
            passed,
            &issues,
        )?;
        let identity = contract_check_receipt_identity(&canonical);
        Ok(Self {
            schema_version: EULER_PROTOCOL_SCHEMA_VERSION,
            checker_id,
            subject,
            context_hash,
            graph_hash,
            passed,
            issues,
            identity,
        })
    }

    #[must_use]
    /// Contract-check receipt schema version.
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    /// Structural checker implementation identity.
    pub fn checker_id(&self) -> &str {
        &self.checker_id
    }

    #[must_use]
    /// Exact Euler contract subject.
    pub const fn subject(&self) -> ContractIdentity {
        self.subject
    }

    #[must_use]
    /// Exact generic Context content hash.
    pub const fn context_hash(&self) -> ContentHash {
        self.context_hash
    }

    #[must_use]
    /// Exact Euler claim-graph content hash.
    pub const fn graph_hash(&self) -> ContentHash {
        self.graph_hash
    }

    #[must_use]
    /// Whether every structural self-consistency check passed.
    pub const fn passed(&self) -> bool {
        self.passed
    }

    #[must_use]
    /// Canonically ordered checker findings.
    pub fn issues(&self) -> &[String] {
        &self.issues
    }

    #[must_use]
    /// Domain-separated receipt content identity.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Canonical bounded receipt bytes.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ContractError> {
        check_receipt_bytes(
            self.schema_version,
            &self.checker_id,
            self.subject,
            self.context_hash,
            self.graph_hash,
            self.passed,
            &self.issues,
        )
    }

    /// Decode and re-admit one exact canonical v1 receipt. The transport does
    /// not carry authority by itself: a consumer must still call
    /// [`Self::verify_subject`] with the exact contract.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, ContractError> {
        if bytes.len() > MAX_CONTRACT_CHECK_RECEIPT_BYTES {
            return Err(ContractError::new(
                "EulerContractCheckReceiptTooLarge",
                "contract-check receipt exceeds its byte budget",
            ));
        }
        let mut reader = ProtocolReader::new(bytes);
        if reader.fixed::<8>("receipt.magic")? != *CONTRACT_CHECK_RECEIPT_MAGIC {
            return Err(ContractError::new(
                "EulerContractCheckReceiptMagic",
                "contract-check receipt magic is not canonical v1",
            ));
        }
        let schema_version = reader.u32("receipt.schema_version")?;
        protocol_migration_policy(schema_version)?;
        let checker_id = reader.text("receipt.checker_id", MAX_EULER_TEXT_BYTES)?;
        if checker_id != CHECKER_ID {
            return Err(ContractError::new(
                "EulerContractCheckReceiptChecker",
                "receipt does not name the exact current checker implementation",
            ));
        }
        let subject =
            ContractIdentity::from_hash(ContentHash(reader.fixed::<32>("receipt.subject")?));
        let context_hash = ContentHash(reader.fixed::<32>("receipt.context_hash")?);
        let graph_hash = ContentHash(reader.fixed::<32>("receipt.graph_hash")?);
        nonzero_hash("receipt.subject", subject.as_hash())?;
        nonzero_hash("receipt.context_hash", context_hash)?;
        nonzero_hash("receipt.graph_hash", graph_hash)?;
        let passed = match reader.byte("receipt.passed")? {
            0 => false,
            1 => true,
            value => {
                return Err(ContractError::new(
                    "EulerContractCheckReceiptBoolean",
                    format!("receipt passed flag has unknown tag {value}"),
                ));
            }
        };
        let issue_count = reader.count("receipt.issues", MAX_EVIDENCE_RECORDS)?;
        let mut issues = Vec::with_capacity(issue_count);
        for _ in 0..issue_count {
            issues.push(reader.text("receipt.issue", MAX_EULER_TEXT_BYTES)?);
        }
        reader.finish("contract-check receipt")?;
        if passed != issues.is_empty() {
            return Err(ContractError::new(
                "EulerContractCheckReceiptInconsistent",
                "passing receipt must have no issues and a failing receipt must name at least one",
            ));
        }
        let decoded = Self::new(subject, context_hash, graph_hash, passed, issues)?;
        if decoded.canonical_bytes()?.as_slice() != bytes {
            return Err(ContractError::new(
                "EulerContractCheckReceiptNonCanonical",
                "receipt bytes do not round-trip to one canonical representation",
            ));
        }
        Ok(decoded)
    }

    /// Recompute the receipt identity and enforce its checker/schema/result
    /// invariants without granting subject authority.
    pub fn verify_identity(&self) -> Result<(), ContractError> {
        let expected = contract_check_receipt_identity(&self.canonical_bytes()?);
        if self.schema_version != EULER_PROTOCOL_SCHEMA_VERSION
            || self.checker_id != CHECKER_ID
            || self.passed != self.issues.is_empty()
            || self.identity != expected
        {
            return Err(ContractError::new(
                "EulerContractCheckReceiptIdentityMismatch",
                "receipt identity, schema, checker, or pass/issues invariant is stale",
            ));
        }
        Ok(())
    }

    /// Re-run the literal-anchor checker and require this receipt to equal its
    /// freshly computed passing result for the supplied subject.
    pub fn verify_subject(&self, contract: &EulerScientificContract) -> Result<(), ContractError> {
        self.verify_identity()?;
        let independently_recomputed = check_frozen_contract(contract)?;
        if !independently_recomputed.passed || self != &independently_recomputed {
            return Err(ContractError::new(
                "EulerContractStaleCheckReceipt",
                "contract-check receipt is not the freshly recomputed passing literal-anchor result for this subject",
            ));
        }
        Ok(())
    }
}

fn write_len(bytes: &mut Vec<u8>, len: usize) -> Result<(), ContractError> {
    let len = u32::try_from(len).map_err(|_| {
        ContractError::new(
            "EulerContractCheckReceiptTooLarge",
            "receipt collection length exceeds u32 transport",
        )
    })?;
    bytes.extend_from_slice(&len.to_le_bytes());
    Ok(())
}

fn write_text(bytes: &mut Vec<u8>, value: &str) -> Result<(), ContractError> {
    write_len(bytes, value.len())?;
    bytes.extend_from_slice(value.as_bytes());
    Ok(())
}

fn check_receipt_bytes(
    schema_version: u32,
    checker_id: &str,
    subject: ContractIdentity,
    context_hash: ContentHash,
    graph_hash: ContentHash,
    passed: bool,
    issues: &[String],
) -> Result<Vec<u8>, ContractError> {
    let mut bytes = Vec::with_capacity(512);
    bytes.extend_from_slice(CONTRACT_CHECK_RECEIPT_MAGIC);
    bytes.extend_from_slice(&schema_version.to_le_bytes());
    write_text(&mut bytes, checker_id)?;
    bytes.extend_from_slice(subject.as_hash().as_bytes());
    bytes.extend_from_slice(context_hash.as_bytes());
    bytes.extend_from_slice(graph_hash.as_bytes());
    bytes.push(u8::from(passed));
    write_len(&mut bytes, issues.len())?;
    for issue in issues {
        write_text(&mut bytes, issue)?;
    }
    if bytes.len() > MAX_CONTRACT_CHECK_RECEIPT_BYTES {
        return Err(ContractError::new(
            "EulerContractCheckReceiptTooLarge",
            "contract-check receipt exceeds its byte budget",
        ));
    }
    Ok(bytes)
}

struct ProtocolReader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> ProtocolReader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn take(&mut self, len: usize, field: &'static str) -> Result<&'a [u8], ContractError> {
        let end = self.offset.checked_add(len).ok_or_else(|| {
            ContractError::new(
                "EulerContractCheckReceiptLengthOverflow",
                format!("{field} length overflows the receipt transport"),
            )
        })?;
        let value = self.bytes.get(self.offset..end).ok_or_else(|| {
            ContractError::new(
                "EulerContractCheckReceiptTruncated",
                format!("{field} is truncated"),
            )
        })?;
        self.offset = end;
        Ok(value)
    }

    fn fixed<const N: usize>(&mut self, field: &'static str) -> Result<[u8; N], ContractError> {
        self.take(N, field)?.try_into().map_err(|_| {
            ContractError::new(
                "EulerContractCheckReceiptTruncated",
                format!("{field} has the wrong fixed width"),
            )
        })
    }

    fn byte(&mut self, field: &'static str) -> Result<u8, ContractError> {
        Ok(self.fixed::<1>(field)?[0])
    }

    fn u32(&mut self, field: &'static str) -> Result<u32, ContractError> {
        Ok(u32::from_le_bytes(self.fixed::<4>(field)?))
    }

    fn count(&mut self, field: &'static str, max: usize) -> Result<usize, ContractError> {
        let value = usize::try_from(self.u32(field)?).map_err(|_| {
            ContractError::new(
                "EulerContractCheckReceiptCardinality",
                format!("{field} cannot be represented on this platform"),
            )
        })?;
        if value > max {
            return Err(ContractError::new(
                "EulerContractCheckReceiptCardinality",
                format!("{field} count {value} exceeds {max}"),
            ));
        }
        Ok(value)
    }

    fn text(&mut self, field: &'static str, max: usize) -> Result<String, ContractError> {
        let len = self.count(field, max)?;
        let bytes = self.take(len, field)?;
        let value = core::str::from_utf8(bytes).map_err(|_| {
            ContractError::new(
                "EulerContractCheckReceiptUtf8",
                format!("{field} is not UTF-8"),
            )
        })?;
        checked_protocol_id(field, value)
    }

    fn finish(self, what: &'static str) -> Result<(), ContractError> {
        if self.offset != self.bytes.len() {
            return Err(ContractError::new(
                "EulerContractCheckReceiptTrailingBytes",
                format!("{what} has trailing bytes"),
            ));
        }
        Ok(())
    }
}

fn literal_frozen_hash(field: &'static str, hex: &str) -> Result<ContentHash, ContractError> {
    ContentHash::from_hex(hex).ok_or_else(|| {
        ContractError::new(
            "EulerContractFrozenDigestInvalid",
            format!("{field} is not an exact 32-byte hexadecimal digest"),
        )
    })
}

/// Re-decode both generic and Euler transports and compare their identities
/// with literal-frozen digests that are not derived by this call.
pub fn check_frozen_contract(
    contract: &EulerScientificContract,
) -> Result<ContractCheckReceipt, ContractError> {
    let mut issues = Vec::new();
    let complete_bytes = contract.canonical_bytes()?;
    match EulerScientificContract::from_canonical_bytes(&complete_bytes) {
        Ok(decoded) if decoded == *contract => {}
        Ok(_) => issues.push("whole-contract-roundtrip-mismatch".to_owned()),
        Err(error) => issues.push(format!("whole-contract-decode:{}", error.code())),
    }
    match ContextOfUse::from_canonical_bytes(contract.context_canonical_bytes()) {
        Ok(decoded) if decoded == *contract.context() => {}
        Ok(_) => issues.push("context-roundtrip-mismatch".to_owned()),
        Err(_) => issues.push("context-decode-refused".to_owned()),
    }
    let graph_bytes = contract.claim_graph().canonical_bytes()?;
    match EulerClaimGraph::from_canonical_bytes(&graph_bytes) {
        Ok(decoded) if decoded == *contract.claim_graph() => {}
        Ok(_) => issues.push("claim-graph-roundtrip-mismatch".to_owned()),
        Err(error) => issues.push(format!("claim-graph-decode:{}", error.code())),
    }
    let graph_hash = contract.claim_graph().content_hash()?;
    let reviewed_context = literal_frozen_hash("frozen context hash", FROZEN_CONTEXT_HASH_HEX)?;
    let reviewed_graph =
        literal_frozen_hash("frozen claim-graph hash", FROZEN_CLAIM_GRAPH_HASH_HEX)?;
    let reviewed_contract =
        literal_frozen_hash("frozen contract identity", FROZEN_CONTRACT_IDENTITY_HEX)?;
    if contract.context_hash() != reviewed_context {
        issues.push("not-the-literal-frozen-v1-context".to_owned());
    }
    if graph_hash != reviewed_graph {
        issues.push("not-the-literal-frozen-v1-claim-graph".to_owned());
    }
    if contract.identity().as_hash() != reviewed_contract {
        issues.push("not-the-literal-frozen-v1-contract".to_owned());
    }
    ContractCheckReceipt::new(
        contract.identity(),
        contract.context_hash(),
        graph_hash,
        issues.is_empty(),
        issues,
    )
}

/// Opaque structural admission wrapper. It exposes no conversion into
/// fs-govern runtime, maturity, physical-validation, or release authority.
///
/// ```compile_fail
/// use fs_euler_disc_e2e::StructurallyAdmittedEulerContract;
/// use fs_govern::evidence_contract::{AuthorityGrant, ProvedAuthority};
/// fn cannot_promote(
///     admitted: StructurallyAdmittedEulerContract,
/// ) -> AuthorityGrant<ProvedAuthority> {
///     admitted.into()
/// }
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct StructurallyAdmittedEulerContract {
    contract: EulerScientificContract,
    receipt: ContractCheckReceipt,
}

impl StructurallyAdmittedEulerContract {
    #[must_use]
    /// Exact structurally admitted frozen contract.
    pub const fn contract(&self) -> &EulerScientificContract {
        &self.contract
    }

    #[must_use]
    /// Passing structural-check receipt.
    pub const fn receipt(&self) -> &ContractCheckReceipt {
        &self.receipt
    }
}

/// Admit only the exact frozen contract after the literal-digest structural
/// check.
pub fn admit_frozen_contract(
    contract: EulerScientificContract,
) -> Result<StructurallyAdmittedEulerContract, ContractError> {
    let receipt = check_frozen_contract(&contract)?;
    if !receipt.passed() {
        return Err(ContractError::new(
            "EulerContractStructuralCheckRefused",
            format!("structural checker issues: {:?}", receipt.issues()),
        ));
    }
    receipt.verify_subject(&contract)?;
    Ok(StructurallyAdmittedEulerContract { contract, receipt })
}

#[allow(dead_code, clippy::too_many_arguments)]
fn classify_claim_evidence_packet_identity_fields(
    packet: &ClaimEvidencePacket,
    record: &EvidenceRecord,
    budget: &ProtocolBudget,
    seed_source: &ProtocolSeed,
    authority_source: &EvidenceAuthorityDeclaration,
    access_class_source: DeclaredEvidenceAccessClass,
    reported_disposition_source: ReportedScientificDisposition,
    assessment_disposition_source: AssessmentDisposition,
) {
    let ClaimEvidencePacket {
        schema_version: _,
        contract_identity: _,
        case_id: _,
        design_set_identity: _,
        aggregate_qoi_derivation_receipt_identity: _,
        claim: _,
        point: _,
        records: _,
        no_claims_accepted: _,
        target_fitted: _,
        reported_scientific_disposition: _,
        expected_disposition: _,
        units: _,
        seed: _,
        budget: _,
        identity: _,
    } = packet;
    let EvidenceRecord {
        contract_identity: _,
        claim: _,
        requirement: _,
        qois: _,
        authority: _,
        artifact_hash: _,
        source_id: _,
        source_schema: _,
        source_kind: _,
        schema_admission_receipt_hash: _,
        access_class: _,
        independent: _,
    } = record;
    let ProtocolBudget {
        max_wall_time_ms: _,
        max_memory_bytes: _,
        normalized_accuracy_limit: _,
    } = budget;
    match seed_source {
        ProtocolSeed::Fixed { value } => {
            let _ = value;
        }
        ProtocolSeed::NotApplicable { reason } => {
            let _ = reason;
        }
    }
    match authority_source {
        EvidenceAuthorityDeclaration::StructuralProcess { receipt_hash } => {
            let _ = receipt_hash;
        }
        EvidenceAuthorityDeclaration::VerifiedNumerics { color }
        | EvidenceAuthorityDeclaration::ValidatedPhysical { color } => {
            let _ = color;
        }
    }
    match access_class_source {
        DeclaredEvidenceAccessClass::NotApplicable
        | DeclaredEvidenceAccessClass::Calibration
        | DeclaredEvidenceAccessClass::Validation
        | DeclaredEvidenceAccessClass::BlindHoldout => {}
    }
    match reported_disposition_source {
        ReportedScientificDisposition::Positive
        | ReportedScientificDisposition::Negative
        | ReportedScientificDisposition::Inconclusive => {}
    }
    match assessment_disposition_source {
        AssessmentDisposition::ReferenceCompleteCandidate
        | AssessmentDisposition::DemotedCandidate
        | AssessmentDisposition::RetainedTerminal
        | AssessmentDisposition::Refused => {}
    }
}

#[allow(dead_code)]
fn classify_prerequisite_receipt_identity_fields(receipt: &PrerequisiteAssessmentReceipt) {
    let PrerequisiteAssessmentReceipt {
        schema_version: _,
        contract_identity: _,
        prerequisite: _,
        dependent: _,
        use_kind: _,
        source_packet_identity: _,
        source_assessment_identity: _,
        source_design_set_identity: _,
        source_point_bytes: _,
        identity: _,
    } = receipt;
}

#[allow(dead_code)]
fn classify_claim_policy_assessment_identity_fields(assessment: &ClaimPolicyAssessment) {
    let ClaimPolicyAssessment {
        schema_version: _,
        contract_identity: _,
        packet_identity: _,
        design_set_identity: _,
        aggregate_qoi_derivation_receipt_identity: _,
        point_bytes: _,
        case_id: _,
        claim: _,
        disposition: _,
        reported_scientific_disposition: _,
        reasons: _,
        log: _,
        identity: _,
    } = assessment;
}

#[allow(dead_code)]
fn classify_contract_check_receipt_identity_fields(receipt: &ContractCheckReceipt) {
    let ContractCheckReceipt {
        schema_version: _,
        checker_id: _,
        subject: _,
        context_hash: _,
        graph_hash: _,
        passed: _,
        issues: _,
        identity: _,
    } = receipt;
}

#[allow(dead_code)]
fn classify_claim_policy_assessment_log_identity_fields(log: &ClaimPolicyAssessmentLog) {
    let ClaimPolicyAssessmentLog {
        json_line: _,
        identity: _,
    } = log;
}

#[cfg(test)]
mod tests {
    #![allow(clippy::too_many_lines)] // Mutation batteries intentionally keep one field oracle together.

    use super::*;

    fn test_hash(label: &str) -> ContentHash {
        fs_blake3::hash_domain(
            "org.frankensim.fs-euler-disc-e2e.protocol-unit-test.v1",
            label.as_bytes(),
        )
    }

    // These batteries prove identity-preimage binding. Alternate codec
    // candidates need not be admissible under the production decoder; the
    // constructor and hostile-transport tests cover admission separately.
    fn assert_protocol_hash_moved(field: &str, base: ContentHash, candidate: ContentHash) {
        assert_ne!(
            candidate, base,
            "semantic protocol identity field {field} did not move the content hash"
        );
    }

    fn test_point(contract: &EulerScientificContract, fraction: f64) -> ApplicabilityPoint {
        let numeric = contract
            .context()
            .applicability()
            .numeric()
            .iter()
            .map(|(axis, domain)| {
                let (lo, hi) = domain.bounds();
                (axis.clone(), lo + (hi - lo) * fraction)
            })
            .collect();
        let categorical = contract
            .context()
            .applicability()
            .categorical()
            .iter()
            .map(|(axis, domain)| {
                (
                    axis.clone(),
                    domain
                        .allowed()
                        .iter()
                        .next()
                        .expect("categorical domain is nonempty")
                        .clone(),
                )
            })
            .collect();
        ApplicabilityPoint::try_new(numeric, categorical).expect("test point")
    }

    fn test_packet(contract: &EulerScientificContract) -> ClaimEvidencePacket {
        ClaimEvidencePacket::try_new(
            contract.identity(),
            "unit-case-a",
            test_hash("unit-design-set"),
            test_hash("unit-aggregate-qoi-derivation-receipt"),
            EulerClaimKind::NumericalTrajectoryVerification,
            test_point(contract, 0.5),
            Vec::new(),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
            vec!["1".to_owned()],
            ProtocolSeed::Fixed { value: 7 },
            ProtocolBudget::try_new(100, 200, 0.25).expect("test budget"),
        )
        .expect("test packet")
    }

    #[test]
    fn evidence_authority_length_preflight_matches_every_canonical_variant() {
        let empty_regime = ValidityDomain::unconstrained();
        assert!(
            empty_regime.bounds().is_empty(),
            "the empty case must exercise the zero-axis validity-domain encoding"
        );
        let boundary_regime = (0..(MAX_VALIDITY_DOMAIN_AXES - 1))
            .fold(ValidityDomain::unconstrained(), |regime, index| {
                regime.with(format!("axis-{index:02}"), -1.0, 1.0)
            });
        let boundary_regime =
            boundary_regime.with("z".repeat(MAX_PROTOCOL_COLOR_IDENTITY_BYTES), -2.0, 2.0);
        assert_eq!(
            boundary_regime.bounds().len(),
            MAX_VALIDITY_DOMAIN_AXES,
            "the boundary case must exercise the exact regime-axis limit"
        );
        assert!(
            boundary_regime
                .bounds()
                .keys()
                .any(|axis| axis.len() == MAX_PROTOCOL_COLOR_IDENTITY_BYTES),
            "the boundary case must exercise the exact regime-axis identity limit"
        );

        let cases = [
            (
                "structural-process",
                EvidenceAuthorityDeclaration::StructuralProcess {
                    receipt_hash: test_hash("authority-length-structural-receipt"),
                },
            ),
            (
                "verified-numerics/verified-color",
                EvidenceAuthorityDeclaration::VerifiedNumerics {
                    color: Color::Verified { lo: 0.25, hi: 0.5 },
                },
            ),
            (
                "validated-physical/verified-color",
                EvidenceAuthorityDeclaration::ValidatedPhysical {
                    color: Color::Verified { lo: 0.5, hi: 0.75 },
                },
            ),
            (
                "verified-numerics/validated-color-empty-regime",
                EvidenceAuthorityDeclaration::VerifiedNumerics {
                    color: Color::Validated {
                        regime: empty_regime.clone(),
                        dataset: "authority-length-empty-regime-a".to_owned(),
                    },
                },
            ),
            (
                "validated-physical/validated-color-empty-regime",
                EvidenceAuthorityDeclaration::ValidatedPhysical {
                    color: Color::Validated {
                        regime: empty_regime,
                        dataset: "authority-length-empty-regime-b".to_owned(),
                    },
                },
            ),
            (
                "verified-numerics/validated-color-boundary-regime",
                EvidenceAuthorityDeclaration::VerifiedNumerics {
                    color: Color::Validated {
                        regime: boundary_regime.clone(),
                        dataset: "authority-length-boundary-regime-a".to_owned(),
                    },
                },
            ),
            (
                "validated-physical/validated-color-boundary-regime",
                EvidenceAuthorityDeclaration::ValidatedPhysical {
                    color: Color::Validated {
                        regime: boundary_regime,
                        dataset: "authority-length-boundary-regime-b".to_owned(),
                    },
                },
            ),
            (
                "verified-numerics/estimated-color",
                EvidenceAuthorityDeclaration::VerifiedNumerics {
                    color: Color::Estimated {
                        estimator: "authority-length-estimator-a".to_owned(),
                        dispersion: 0.125,
                    },
                },
            ),
            (
                "validated-physical/estimated-color",
                EvidenceAuthorityDeclaration::ValidatedPhysical {
                    color: Color::Estimated {
                        estimator: "authority-length-estimator-b".to_owned(),
                        dispersion: 0.25,
                    },
                },
            ),
        ];

        // This is a codec-size mirror test, so the table intentionally includes
        // the canonical zero-axis Validated encoding even though later payload
        // admission rejects that scientifically unusable regime.
        for (case, authority) in cases {
            let preflight_len = evidence_authority_canonical_len(&authority)
                .unwrap_or_else(|error| panic!("{case}: length preflight failed: {error}"));
            assert_eq!(
                preflight_len,
                authority.canonical_bytes().len(),
                "{case}: preflight length diverged from the canonical writer"
            );
        }
    }

    #[test]
    fn zero_contract_identities_refuse_before_identity_publication() {
        let contract = build_frozen_contract().expect("frozen contract");
        let zero_contract = ContractIdentity::from_hash(ContentHash([0; 32]));
        let claim = EulerClaimKind::NumericalTrajectoryVerification;
        let requirement = EvidenceRequirement::CodeVerification;
        let qoi = contract
            .context()
            .qois()
            .keys()
            .next()
            .expect("frozen QoI")
            .clone();
        let error = EvidenceRecord::try_new(
            zero_contract,
            claim,
            requirement,
            vec![qoi],
            EvidenceAuthorityDeclaration::StructuralProcess {
                receipt_hash: test_hash("nonzero-role-receipt"),
            },
            test_hash("nonzero-artifact"),
            "zero-contract-record",
            requirement.source_schema(),
            requirement.source_kind(),
            test_hash("nonzero-schema-receipt"),
            DeclaredEvidenceAccessClass::NotApplicable,
            true,
        )
        .expect_err("zero evidence contract identity must refuse");
        assert_eq!(error.code(), "EulerProtocolZeroIdentity");

        let error = ClaimEvidencePacket::try_new(
            zero_contract,
            "zero-contract-packet",
            test_hash("zero-contract-design-set"),
            test_hash("zero-contract-aggregate-qoi-derivation-receipt"),
            claim,
            test_point(&contract, 0.5),
            Vec::new(),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
            vec!["unit-a".to_owned()],
            ProtocolSeed::Fixed { value: 7 },
            ProtocolBudget::try_new(100, 200, 0.25).expect("test budget"),
        )
        .expect_err("zero packet contract identity must refuse");
        assert_eq!(error.code(), "EulerProtocolZeroIdentity");
    }

    #[test]
    fn packet_requires_distinct_nonzero_design_and_qoi_derivation_identities() {
        let contract = build_frozen_contract().expect("frozen contract");
        let claim = EulerClaimKind::NumericalTrajectoryVerification;
        let construct = |design_set_identity, aggregate_qoi_derivation_receipt_identity| {
            ClaimEvidencePacket::try_new(
                contract.identity(),
                "packet-binding-identity-refusal",
                design_set_identity,
                aggregate_qoi_derivation_receipt_identity,
                claim,
                test_point(&contract, 0.5),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["1".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                ProtocolBudget::try_new(100, 200, 0.25).expect("test budget"),
            )
        };

        let zero = ContentHash([0; 32]);
        let design_set = test_hash("nonzero-design-set");
        let derivation = test_hash("nonzero-aggregate-qoi-derivation-receipt");
        for error in [
            construct(zero, derivation).expect_err("zero design set must refuse"),
            construct(design_set, zero).expect_err("zero derivation receipt must refuse"),
        ] {
            assert_eq!(error.code(), "EulerProtocolZeroIdentity");
        }
        let error = construct(design_set, design_set)
            .expect_err("one identity cannot satisfy two semantic roles");
        assert_eq!(error.code(), "EulerProtocolCrossRoleEvidenceAlias");
    }

    #[test]
    fn color_identity_preflight_preserves_the_boundary_and_bounds_refusal_detail() {
        fn validated_authority(dataset: String) -> EvidenceAuthorityDeclaration {
            EvidenceAuthorityDeclaration::ValidatedPhysical {
                color: Color::Validated {
                    regime: ValidityDomain::unconstrained().with("outer-radius", 0.0, 1.0),
                    dataset,
                },
            }
        }

        fn estimated_authority(estimator: String) -> EvidenceAuthorityDeclaration {
            EvidenceAuthorityDeclaration::VerifiedNumerics {
                color: Color::Estimated {
                    estimator,
                    dispersion: 0.125,
                },
            }
        }

        validated_authority("d".repeat(MAX_PROTOCOL_COLOR_IDENTITY_BYTES))
            .try_new()
            .expect("a dataset identity exactly at the shared byte limit must remain admissible");
        estimated_authority("e".repeat(MAX_PROTOCOL_COLOR_IDENTITY_BYTES))
            .try_new()
            .expect(
                "an estimator identity exactly at the shared byte limit must remain admissible",
            );

        let exact_axis = EvidenceAuthorityDeclaration::ValidatedPhysical {
            color: Color::Validated {
                regime: ValidityDomain::unconstrained().with(
                    "a".repeat(MAX_PROTOCOL_COLOR_IDENTITY_BYTES),
                    0.0,
                    1.0,
                ),
                dataset: "dataset-for-axis-boundary".to_owned(),
            },
        };
        exact_axis.try_new().expect(
            "a regime axis identity exactly at the shared byte limit must remain admissible",
        );

        let dataset_plus_one =
            validated_authority("d".repeat(MAX_PROTOCOL_COLOR_IDENTITY_BYTES + 1))
                .try_new()
                .expect_err("a dataset identity one byte over the limit must refuse");
        let dataset_very_large =
            validated_authority("d".repeat(MAX_PROTOCOL_COLOR_IDENTITY_BYTES * 4_096))
                .try_new()
                .expect_err("a very large dataset identity must refuse before upstream cloning");
        assert_eq!(dataset_plus_one.code(), "EulerProtocolMalformedColor");
        assert_eq!(dataset_very_large.code(), "EulerProtocolMalformedColor");
        assert_eq!(
            dataset_plus_one.detail(),
            "validated-color dataset identity exceeds the v1 byte limit of 256"
        );
        assert_eq!(dataset_very_large.detail(), dataset_plus_one.detail());

        let estimator_plus_one =
            estimated_authority("e".repeat(MAX_PROTOCOL_COLOR_IDENTITY_BYTES + 1))
                .try_new()
                .expect_err("an estimator identity one byte over the limit must refuse");
        let estimator_very_large =
            estimated_authority("e".repeat(MAX_PROTOCOL_COLOR_IDENTITY_BYTES * 4_096))
                .try_new()
                .expect_err("a very large estimator identity must refuse before upstream cloning");
        assert_eq!(estimator_plus_one.code(), "EulerProtocolMalformedColor");
        assert_eq!(estimator_very_large.code(), "EulerProtocolMalformedColor");
        assert_eq!(
            estimator_plus_one.detail(),
            "estimated-color estimator identity exceeds the v1 byte limit of 256"
        );
        assert_eq!(estimator_very_large.detail(), estimator_plus_one.detail());

        let axis_plus_one = EvidenceAuthorityDeclaration::ValidatedPhysical {
            color: Color::Validated {
                regime: ValidityDomain::unconstrained().with(
                    "a".repeat(MAX_PROTOCOL_COLOR_IDENTITY_BYTES + 1),
                    0.0,
                    1.0,
                ),
                dataset: "dataset-for-axis-plus-one".to_owned(),
            },
        }
        .try_new()
        .expect_err("a regime axis identity one byte over the limit must refuse");
        let axis_very_large = EvidenceAuthorityDeclaration::ValidatedPhysical {
            color: Color::Validated {
                regime: ValidityDomain::unconstrained().with(
                    "a".repeat(MAX_PROTOCOL_COLOR_IDENTITY_BYTES * 4_096),
                    0.0,
                    1.0,
                ),
                dataset: "dataset-for-large-axis".to_owned(),
            },
        }
        .try_new()
        .expect_err("a very large regime axis identity must refuse before upstream cloning");
        assert_eq!(axis_plus_one.code(), "EulerProtocolMalformedColor");
        assert_eq!(axis_very_large.code(), "EulerProtocolMalformedColor");
        assert_eq!(
            axis_plus_one.detail(),
            "validated-color regime axis identity exceeds the v1 byte limit of 256"
        );
        assert_eq!(axis_very_large.detail(), axis_plus_one.detail());

        // The exact equality checks above prove that public error detail does
        // not scale with rejected input. Keep an explicit size guard as a
        // backstop against future attempts to make the diagnostic verbose.
        assert!(dataset_plus_one.detail().len() < MAX_PROTOCOL_COLOR_IDENTITY_BYTES);
        assert!(estimator_plus_one.detail().len() < MAX_PROTOCOL_COLOR_IDENTITY_BYTES);
        assert!(axis_plus_one.detail().len() < MAX_PROTOCOL_COLOR_IDENTITY_BYTES);
    }

    #[test]
    fn direct_protocol_seed_variant_validates_borrowed_bounded_reason() {
        let exact = ProtocolSeed::NotApplicable {
            reason: "s".repeat(MAX_PROTOCOL_ID_BYTES),
        };
        exact
            .validate()
            .expect("a direct seed reason exactly at the byte limit must remain admissible");

        let plus_one = ProtocolSeed::NotApplicable {
            reason: "s".repeat(MAX_PROTOCOL_ID_BYTES + 1),
        }
        .validate()
        .expect_err("a direct seed reason one byte over the limit must refuse");
        let very_large = ProtocolSeed::NotApplicable {
            reason: "s".repeat(MAX_PROTOCOL_ID_BYTES * 4_096),
        }
        .validate()
        .expect_err("a very large direct seed reason must refuse without cloning it");

        assert_eq!(plus_one.code(), "EulerProtocolInvalidIdentity");
        assert_eq!(very_large.code(), "EulerProtocolInvalidIdentity");
        assert_eq!(
            plus_one.detail(),
            "packet.seed.reason must be a bounded canonical machine identity"
        );
        assert_eq!(very_large.detail(), plus_one.detail());
        assert!(plus_one.detail().len() < MAX_PROTOCOL_ID_BYTES);
    }

    #[test]
    fn claim_evidence_packet_signed_zero_inputs_are_nonsemantic() {
        fn point_with_zero(contract: &EulerScientificContract, zero: f64) -> ApplicabilityPoint {
            let base = test_point(contract, 0.5);
            let zero_axis = base
                .numeric()
                .keys()
                .next()
                .expect("frozen numeric axis")
                .clone();
            let numeric = base
                .numeric()
                .iter()
                .map(|(axis, value)| (axis.clone(), if axis == &zero_axis { zero } else { *value }))
                .collect();
            let categorical = base
                .categorical()
                .iter()
                .map(|(axis, value)| (axis.clone(), value.clone()))
                .collect();
            ApplicabilityPoint::try_new(numeric, categorical).expect("signed-zero point")
        }

        let contract = build_frozen_contract().expect("frozen contract");
        let make_packet = |zero| {
            ClaimEvidencePacket::try_new(
                contract.identity(),
                "signed-zero-case",
                test_hash("signed-zero-design-set"),
                test_hash("signed-zero-aggregate-qoi-derivation-receipt"),
                EulerClaimKind::NumericalTrajectoryVerification,
                point_with_zero(&contract, zero),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                ProtocolBudget::try_new(100, 200, 0.25).expect("test budget"),
            )
            .expect("signed-zero packet")
        };
        let positive = make_packet(0.0);
        let negative = make_packet(-0.0);
        let positive_zero = positive
            .point()
            .numeric()
            .values()
            .find(|value| **value == 0.0)
            .expect("positive zero retained");
        let negative_zero = negative
            .point()
            .numeric()
            .values()
            .find(|value| **value == 0.0)
            .expect("negative zero retained");
        assert_ne!(positive_zero.to_bits(), negative_zero.to_bits());
        assert_eq!(
            positive.canonical_bytes().expect("positive-zero bytes"),
            negative.canonical_bytes().expect("negative-zero bytes")
        );
        assert_eq!(positive.identity(), negative.identity());

        let make_budget_packet = |normalized_accuracy_limit| {
            ClaimEvidencePacket::try_new(
                contract.identity(),
                "signed-zero-budget-case",
                test_hash("signed-zero-budget-design-set"),
                test_hash("signed-zero-budget-aggregate-qoi-derivation-receipt"),
                EulerClaimKind::NumericalTrajectoryVerification,
                test_point(&contract, 0.5),
                Vec::new(),
                true,
                false,
                ReportedScientificDisposition::Positive,
                AssessmentDisposition::Refused,
                vec!["unit-a".to_owned()],
                ProtocolSeed::Fixed { value: 7 },
                ProtocolBudget::try_new(100, 200, normalized_accuracy_limit)
                    .expect("signed-zero budget"),
            )
            .expect("signed-zero budget packet")
        };
        let positive_budget = make_budget_packet(0.0);
        let negative_budget = make_budget_packet(-0.0);
        assert_eq!(
            positive_budget.budget.normalized_accuracy_limit().to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            negative_budget.budget.normalized_accuracy_limit().to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            positive_budget
                .canonical_bytes()
                .expect("positive-zero budget bytes"),
            negative_budget
                .canonical_bytes()
                .expect("negative-zero budget bytes")
        );
        assert_eq!(positive_budget.identity(), negative_budget.identity());
    }

    fn packet_hash(packet: &ClaimEvidencePacket) -> ContentHash {
        claim_evidence_packet_identity(&packet.canonical_bytes().expect("packet bytes"))
    }

    #[test]
    fn claim_evidence_packet_identity_semantic_fields_move_independently() {
        fn replace_unique_bytes(
            field: &str,
            source: &[u8],
            canonical: &[u8],
            alternate: &[u8],
        ) -> Vec<u8> {
            assert_eq!(
                canonical.len(),
                alternate.len(),
                "{field} alternate codec must preserve transport length"
            );
            let offsets = source
                .windows(canonical.len())
                .enumerate()
                .filter_map(|(offset, window)| (window == canonical).then_some(offset))
                .collect::<Vec<_>>();
            assert_eq!(
                offsets.len(),
                1,
                "{field} canonical field must occur exactly once"
            );
            let mut bytes = source.to_vec();
            bytes[offsets[0]..offsets[0] + alternate.len()].copy_from_slice(alternate);
            bytes
        }

        let contract = build_frozen_contract().expect("frozen contract");
        let base_packet = test_packet(&contract);
        let bytes = base_packet.canonical_bytes().expect("packet bytes");
        let base = packet_hash(&base_packet);
        assert_eq!(base, base_packet.identity);

        let mut variants = Vec::<(&str, ClaimEvidencePacket)>::new();
        let mut candidate = base_packet.clone();
        candidate.schema_version += 1;
        variants.push(("protocol-schema-version", candidate));
        let mut candidate = base_packet.clone();
        candidate.contract_identity = ContractIdentity::from_hash(test_hash("other-contract"));
        variants.push(("contract-identity", candidate));
        let mut candidate = base_packet.clone();
        candidate.case_id = "unit-case-b".to_owned();
        variants.push(("case-id", candidate));
        let mut candidate = base_packet.clone();
        candidate.design_set_identity = test_hash("other-design-set");
        variants.push(("design-set-identity", candidate));
        let mut candidate = base_packet.clone();
        candidate.aggregate_qoi_derivation_receipt_identity =
            test_hash("other-aggregate-qoi-derivation-receipt");
        variants.push(("aggregate-qoi-derivation-receipt-identity", candidate));
        let mut candidate = base_packet.clone();
        candidate.claim = EulerClaimKind::CalibratedReproduction;
        variants.push(("claim-kind", candidate));
        let mut candidate = base_packet.clone();
        candidate.point = test_point(&contract, 0.25);
        variants.push(("applicability-point-anchor", candidate));
        let mut candidate = base_packet.clone();
        candidate.no_claims_accepted = false;
        variants.push(("no-claim-acceptance", candidate));
        let mut candidate = base_packet.clone();
        candidate.target_fitted = true;
        variants.push(("target-fitting-state", candidate));
        let mut candidate = base_packet.clone();
        candidate.reported_scientific_disposition = ReportedScientificDisposition::Negative;
        variants.push(("reported-scientific-disposition", candidate));
        let mut candidate = base_packet.clone();
        candidate.expected_disposition = AssessmentDisposition::DemotedCandidate;
        variants.push(("expected-disposition", candidate));
        let mut candidate = base_packet.clone();
        candidate.units = vec!["unit-b".to_owned()];
        variants.push(("unit-set", candidate));
        let mut candidate = base_packet.clone();
        candidate.seed = ProtocolSeed::Fixed { value: 8 };
        variants.push(("seed-fixed-value", candidate));
        let mut candidate = base_packet.clone();
        candidate.seed = ProtocolSeed::NotApplicable {
            reason: "identity-test-no-seed".to_owned(),
        };
        variants.push(("seed-variant-and-reason", candidate));
        let mut candidate = base_packet.clone();
        candidate.budget.max_wall_time_ms += 1;
        variants.push(("protocol-budget-wall-time", candidate));
        let mut candidate = base_packet.clone();
        candidate.budget.max_memory_bytes += 1;
        variants.push(("protocol-budget-memory", candidate));
        let mut candidate = base_packet.clone();
        candidate.budget.normalized_accuracy_limit = 0.5;
        variants.push(("protocol-budget-normalized-accuracy", candidate));

        let requirement = EvidenceRequirement::SolutionVerification;
        let color = Color::Verified { lo: 1.25, hi: 2.75 };
        let admission_receipt = test_hash("schema-receipt");
        let record = EvidenceRecord::try_new(
            contract.identity(),
            base_packet.claim,
            requirement,
            vec![
                contract
                    .context()
                    .qois()
                    .keys()
                    .next()
                    .expect("frozen QoI registry")
                    .clone(),
            ],
            EvidenceAuthorityDeclaration::VerifiedNumerics {
                color: color.clone(),
            },
            test_hash("artifact"),
            "unit-source",
            requirement.source_schema(),
            requirement.source_kind(),
            admission_receipt,
            DeclaredEvidenceAccessClass::NotApplicable,
            true,
        )
        .expect("valid identity-test evidence row");
        let mut with_record = base_packet.clone();
        with_record.records.insert(requirement, record);
        variants.push(("evidence-registry", with_record.clone()));

        for (field, candidate) in variants {
            assert_protocol_hash_moved(field, base, packet_hash(&candidate));
        }

        let with_record_hash = packet_hash(&with_record);
        let mut record_variants = Vec::<(&str, ClaimEvidencePacket)>::new();
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .contract_identity = ContractIdentity::from_hash(test_hash("record-contract"));
        record_variants.push(("evidence-record-contract-identity", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .claim = EulerClaimKind::CalibratedReproduction;
        record_variants.push(("evidence-record-claim", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .requirement = EvidenceRequirement::UncertaintyClosure;
        record_variants.push(("evidence-record-requirement", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .qois = vec![QoiId::try_new("alternate-qoi").expect("alternate QoI")];
        record_variants.push(("evidence-record-qois", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .authority = EvidenceAuthorityDeclaration::StructuralProcess {
            receipt_hash: test_hash("alternate-role-receipt"),
        };
        record_variants.push(("evidence-record-authority-variant", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .authority = EvidenceAuthorityDeclaration::VerifiedNumerics {
            color: Color::Verified { lo: 1.5, hi: 2.75 },
        };
        record_variants.push(("evidence-record-authority-payload", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .artifact_hash = test_hash("alternate-artifact");
        record_variants.push(("evidence-record-artifact-identity", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .source_id = "alternate-source".to_owned();
        record_variants.push(("evidence-record-source-id", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .source_schema = "org.frankensim.alternate-source-schema.v1".to_owned();
        record_variants.push(("evidence-record-source-schema", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .source_kind = ArtifactKind::ContextOfUse;
        record_variants.push(("evidence-record-source-kind", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .schema_admission_receipt_hash = test_hash("alternate-schema-receipt");
        record_variants.push(("evidence-record-schema-receipt", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .access_class = DeclaredEvidenceAccessClass::Calibration;
        record_variants.push(("evidence-record-access-class", candidate));
        let mut candidate = with_record.clone();
        candidate
            .records
            .get_mut(&requirement)
            .expect("identity record")
            .independent = false;
        record_variants.push(("evidence-record-independence", candidate));
        for (field, candidate) in record_variants {
            assert_protocol_hash_moved(field, with_record_hash, packet_hash(&candidate));
        }

        let with_record_bytes = with_record
            .canonical_bytes()
            .expect("packet with evidence bytes");
        let canonical_color = color.canonical_bytes();
        let mut big_endian_color = canonical_color.clone();
        big_endian_color[10..18].copy_from_slice(&1.25_f64.to_bits().to_be_bytes());
        big_endian_color[26..34].copy_from_slice(&2.75_f64.to_bits().to_be_bytes());
        let changed_color_codec = replace_unique_bytes(
            "color-canonical-codec",
            &with_record_bytes,
            &canonical_color,
            &big_endian_color,
        );
        assert_protocol_hash_moved(
            "color-canonical-codec",
            with_record_hash,
            claim_evidence_packet_identity(&changed_color_codec),
        );
        let mut canonical_kind_and_receipt = Vec::with_capacity(33);
        canonical_kind_and_receipt.push(requirement.source_kind().canonical_wire_tag());
        canonical_kind_and_receipt.extend_from_slice(admission_receipt.as_bytes());
        let mut alternate_kind_and_receipt = canonical_kind_and_receipt.clone();
        alternate_kind_and_receipt[0] = alternate_kind_and_receipt[0].wrapping_add(1);
        let changed_kind_codec = replace_unique_bytes(
            "artifact-kind-wire-tags",
            &with_record_bytes,
            &canonical_kind_and_receipt,
            &alternate_kind_and_receipt,
        );
        assert_protocol_hash_moved(
            "artifact-kind-wire-tags",
            with_record_hash,
            claim_evidence_packet_identity(&changed_kind_codec),
        );

        let mut wrong_magic = bytes.clone();
        wrong_magic[0] ^= 1;
        let mut big_endian = bytes.clone();
        big_endian[8..12].copy_from_slice(&1_u32.to_be_bytes());
        let mut reordered = Vec::with_capacity(bytes.len());
        reordered.extend_from_slice(&bytes[..8]);
        reordered.extend_from_slice(&bytes[12..44]);
        reordered.extend_from_slice(&bytes[8..12]);
        reordered.extend_from_slice(&bytes[44..]);
        let mut unframed = Vec::with_capacity(bytes.len() - 4);
        unframed.extend_from_slice(&bytes[..44]);
        unframed.extend_from_slice(&bytes[48..]);
        for (field, candidate) in [
            (
                "identity-domain",
                fs_blake3::hash_domain("org.frankensim.fs-euler-disc-e2e.other-packet.v1", &bytes),
            ),
            (
                "identity-version",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.claim-evidence-packet.v2",
                    &bytes,
                ),
            ),
            (
                "transport-magic",
                fs_blake3::hash_domain(EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, &wrong_magic),
            ),
            (
                "canonical-field-order",
                fs_blake3::hash_domain(EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, &reordered),
            ),
            (
                "length-framing",
                fs_blake3::hash_domain(EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, &unframed),
            ),
            (
                "fixed-numeric-little-endian",
                fs_blake3::hash_domain(EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, &big_endian),
            ),
        ] {
            assert_protocol_hash_moved(field, base, candidate);
        }
    }

    #[test]
    fn prerequisite_receipt_identity_semantic_fields_move_independently() {
        let contract = build_frozen_contract().expect("frozen contract");
        let packet = test_packet(&contract);
        let point_bytes = applicability_point_bytes(packet.point()).expect("point bytes");
        let mut receipt = PrerequisiteAssessmentReceipt {
            schema_version: EULER_PROTOCOL_SCHEMA_VERSION,
            contract_identity: contract.identity(),
            prerequisite: EulerClaimKind::NumericalTrajectoryVerification,
            dependent: EulerClaimKind::CalibratedReproduction,
            use_kind: EvidenceUse::ValidationInput,
            source_packet_identity: packet.identity(),
            source_assessment_identity: test_hash("assessment"),
            source_design_set_identity: packet.design_set_identity(),
            source_point_bytes: point_bytes,
            identity: test_hash("placeholder"),
        };
        receipt.identity = prerequisite_assessment_receipt_identity(
            &receipt.canonical_bytes().expect("receipt bytes"),
        );
        let bytes = receipt.canonical_bytes().expect("receipt bytes");
        let base = receipt.identity;
        receipt.verify().expect("base receipt identity");

        let mut variants = Vec::<(&str, PrerequisiteAssessmentReceipt)>::new();
        let mut candidate = receipt.clone();
        candidate.schema_version += 1;
        variants.push(("protocol-schema-version", candidate));
        let mut candidate = receipt.clone();
        candidate.contract_identity = ContractIdentity::from_hash(test_hash("other-contract"));
        variants.push(("contract-identity", candidate));
        let mut candidate = receipt.clone();
        candidate.prerequisite = EulerClaimKind::CalibratedReproduction;
        variants.push(("prerequisite-claim", candidate));
        let mut candidate = receipt.clone();
        candidate.dependent = EulerClaimKind::BlindTrajectoryPrediction;
        variants.push(("dependent-claim", candidate));
        let mut candidate = receipt.clone();
        candidate.use_kind = EvidenceUse::CalibrationInput;
        variants.push(("evidence-use", candidate));
        let mut candidate = receipt.clone();
        candidate.source_packet_identity = test_hash("other-packet");
        variants.push(("source-packet-identity", candidate));
        let mut candidate = receipt.clone();
        candidate.source_assessment_identity = test_hash("other-assessment");
        variants.push(("source-assessment-identity", candidate));
        let mut candidate = receipt.clone();
        candidate.source_design_set_identity = test_hash("other-design-set");
        variants.push(("source-design-set-identity", candidate));
        let mut candidate = receipt.clone();
        candidate.source_point_bytes.push(0);
        variants.push(("source-applicability-point-anchor", candidate));
        for (field, candidate) in variants {
            assert_protocol_hash_moved(
                field,
                base,
                prerequisite_assessment_receipt_identity(
                    &candidate
                        .canonical_bytes()
                        .expect("candidate receipt bytes"),
                ),
            );
        }

        let mut wrong_magic = bytes.clone();
        wrong_magic[0] ^= 1;
        let mut big_endian = bytes.clone();
        big_endian[8..12].copy_from_slice(&1_u32.to_be_bytes());
        let mut reordered = Vec::with_capacity(bytes.len());
        reordered.extend_from_slice(&bytes[..8]);
        reordered.extend_from_slice(&bytes[12..44]);
        reordered.extend_from_slice(&bytes[8..12]);
        reordered.extend_from_slice(&bytes[44..]);
        let mut unframed = Vec::with_capacity(bytes.len() - 4);
        unframed.extend_from_slice(&bytes[..44]);
        unframed.extend_from_slice(&bytes[48..]);
        for (field, candidate) in [
            (
                "identity-domain",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.other-prerequisite.v1",
                    &bytes,
                ),
            ),
            (
                "identity-version",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.prerequisite-assessment-receipt.v2",
                    &bytes,
                ),
            ),
            (
                "transport-magic",
                fs_blake3::hash_domain(EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN, &wrong_magic),
            ),
            (
                "canonical-field-order",
                fs_blake3::hash_domain(EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN, &reordered),
            ),
            (
                "length-framing",
                fs_blake3::hash_domain(EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN, &unframed),
            ),
            (
                "fixed-numeric-little-endian",
                fs_blake3::hash_domain(EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN, &big_endian),
            ),
        ] {
            assert_protocol_hash_moved(field, base, candidate);
        }
    }

    #[test]
    fn claim_policy_assessment_identity_semantic_fields_move_independently() {
        let contract = build_frozen_contract().expect("frozen contract");
        let packet = test_packet(&contract);
        let assessment = assess_packet(&contract, &packet, &[])
            .expect("evaluator-reachable baseline assessment");
        let canonical = assessment_canonical_bytes(
            assessment.schema_version,
            assessment.contract_identity,
            assessment.packet_identity,
            assessment.design_set_identity,
            assessment.aggregate_qoi_derivation_receipt_identity,
            &assessment.point_bytes,
            &assessment.case_id,
            assessment.claim,
            assessment.disposition,
            assessment.reported_scientific_disposition,
            &assessment.reasons,
            assessment.log.identity,
        )
        .expect("assessment bytes");
        let base = assessment.identity;
        assessment.verify_identity().expect("assessment identity");

        let mut variants = Vec::<(&str, ClaimPolicyAssessment)>::new();
        let mut candidate = assessment.clone();
        candidate.schema_version += 1;
        variants.push(("protocol-schema-version", candidate));
        let mut candidate = assessment.clone();
        candidate.contract_identity = ContractIdentity::from_hash(test_hash("other-contract"));
        variants.push(("contract-identity", candidate));
        let mut candidate = assessment.clone();
        candidate.packet_identity = test_hash("other-packet");
        variants.push(("packet-identity", candidate));
        let mut candidate = assessment.clone();
        candidate.design_set_identity = test_hash("other-assessment-design-set");
        variants.push(("design-set-identity", candidate));
        let mut candidate = assessment.clone();
        candidate.aggregate_qoi_derivation_receipt_identity =
            test_hash("other-assessment-aggregate-qoi-derivation-receipt");
        variants.push(("aggregate-qoi-derivation-receipt-identity", candidate));
        let mut candidate = assessment.clone();
        candidate.point_bytes.push(0);
        variants.push(("applicability-point-anchor", candidate));
        let mut candidate = assessment.clone();
        candidate.case_id = "unit-case-b".to_owned();
        variants.push(("case-id", candidate));
        let mut candidate = assessment.clone();
        candidate.claim = EulerClaimKind::CalibratedReproduction;
        variants.push(("claim-kind", candidate));
        let mut candidate = assessment.clone();
        candidate.disposition = AssessmentDisposition::DemotedCandidate;
        variants.push(("assessment-disposition", candidate));
        let mut candidate = assessment.clone();
        candidate.reported_scientific_disposition = ReportedScientificDisposition::Negative;
        variants.push(("reported-scientific-disposition", candidate));
        let mut candidate = assessment.clone();
        candidate
            .reasons
            .push("binding-no-claims-not-accepted".to_owned());
        candidate.reasons.sort();
        variants.push(("reason-registry", candidate));

        for (field, candidate) in variants {
            let candidate_bytes = assessment_canonical_bytes(
                candidate.schema_version,
                candidate.contract_identity,
                candidate.packet_identity,
                candidate.design_set_identity,
                candidate.aggregate_qoi_derivation_receipt_identity,
                &candidate.point_bytes,
                &candidate.case_id,
                candidate.claim,
                candidate.disposition,
                candidate.reported_scientific_disposition,
                &candidate.reasons,
                candidate.log.identity,
            )
            .expect("candidate assessment bytes");
            assert_protocol_hash_moved(
                field,
                base,
                claim_policy_assessment_identity(&candidate_bytes),
            );
        }
        let alternate_log_identity = fs_blake3::hash_domain(
            "org.frankensim.fs-euler-disc-e2e.other-assessment-log.v1",
            assessment.log.json_line.as_bytes(),
        );
        let candidate_bytes = assessment_canonical_bytes(
            assessment.schema_version,
            assessment.contract_identity,
            assessment.packet_identity,
            assessment.design_set_identity,
            assessment.aggregate_qoi_derivation_receipt_identity,
            &assessment.point_bytes,
            &assessment.case_id,
            assessment.claim,
            assessment.disposition,
            assessment.reported_scientific_disposition,
            &assessment.reasons,
            alternate_log_identity,
        )
        .expect("alternate log-identity assessment bytes");
        assert_protocol_hash_moved(
            "assessment-log-identity",
            base,
            claim_policy_assessment_identity(&candidate_bytes),
        );

        let bytes = canonical;
        let mut wrong_magic = bytes.clone();
        wrong_magic[0] ^= 1;
        let mut big_endian = bytes.clone();
        big_endian[8..12].copy_from_slice(&1_u32.to_be_bytes());
        let mut reordered = Vec::with_capacity(bytes.len());
        reordered.extend_from_slice(&bytes[..8]);
        reordered.extend_from_slice(&bytes[12..44]);
        reordered.extend_from_slice(&bytes[8..12]);
        reordered.extend_from_slice(&bytes[44..]);
        let point_length_offset = 8 + 4 + 32 + 32 + 32 + 32;
        let mut unframed = Vec::with_capacity(bytes.len() - 4);
        unframed.extend_from_slice(&bytes[..point_length_offset]);
        unframed.extend_from_slice(&bytes[point_length_offset + 4..]);
        for (field, candidate) in [
            (
                "identity-domain",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.other-assessment.v1",
                    &bytes,
                ),
            ),
            (
                "identity-version",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.claim-policy-assessment.v2",
                    &bytes,
                ),
            ),
            (
                "transport-magic",
                fs_blake3::hash_domain(EULER_ASSESSMENT_IDENTITY_DOMAIN, &wrong_magic),
            ),
            (
                "canonical-field-order",
                fs_blake3::hash_domain(EULER_ASSESSMENT_IDENTITY_DOMAIN, &reordered),
            ),
            (
                "length-framing",
                fs_blake3::hash_domain(EULER_ASSESSMENT_IDENTITY_DOMAIN, &unframed),
            ),
            (
                "fixed-numeric-little-endian",
                fs_blake3::hash_domain(EULER_ASSESSMENT_IDENTITY_DOMAIN, &big_endian),
            ),
        ] {
            assert_protocol_hash_moved(field, base, candidate);
        }
    }

    #[test]
    fn contract_check_receipt_identity_semantic_fields_move_independently() {
        let contract = build_frozen_contract().expect("frozen contract");
        let receipt = ContractCheckReceipt::new(
            contract.identity(),
            contract.context_hash(),
            contract.claim_graph().content_hash().expect("graph hash"),
            false,
            vec!["unit-issue-a".to_owned()],
        )
        .expect("check receipt");
        let bytes = receipt.canonical_bytes().expect("receipt bytes");
        let base = receipt.identity;

        let mut variants = Vec::<(&str, ContractCheckReceipt)>::new();
        let mut candidate = receipt.clone();
        candidate.schema_version += 1;
        variants.push(("protocol-schema-version", candidate));
        let mut candidate = receipt.clone();
        candidate.checker_id = "unit-checker-b".to_owned();
        variants.push(("checker-id", candidate));
        let mut candidate = receipt.clone();
        candidate.subject = ContractIdentity::from_hash(test_hash("other-subject"));
        variants.push(("subject-identity", candidate));
        let mut candidate = receipt.clone();
        candidate.context_hash = test_hash("other-context");
        variants.push(("context-hash", candidate));
        let mut candidate = receipt.clone();
        candidate.graph_hash = test_hash("other-graph");
        variants.push(("graph-hash", candidate));
        let mut candidate = receipt.clone();
        candidate.passed = true;
        variants.push(("pass-flag", candidate));
        let mut candidate = receipt.clone();
        candidate.issues.push("unit-issue-b".to_owned());
        variants.push(("issue-registry", candidate));
        for (field, candidate) in variants {
            assert_protocol_hash_moved(
                field,
                base,
                contract_check_receipt_identity(
                    &candidate.canonical_bytes().expect("candidate check bytes"),
                ),
            );
        }

        let mut wrong_magic = bytes.clone();
        wrong_magic[0] ^= 1;
        let mut big_endian = bytes.clone();
        big_endian[8..12].copy_from_slice(&1_u32.to_be_bytes());
        let checker_length = usize::try_from(u32::from_le_bytes(
            bytes[12..16].try_into().expect("checker length bytes"),
        ))
        .expect("checker length");
        let subject_offset = 16 + checker_length;
        let context_offset = subject_offset + 32;
        let graph_offset = context_offset + 32;
        let mut reordered = bytes.clone();
        reordered[context_offset..graph_offset]
            .copy_from_slice(&bytes[graph_offset..graph_offset + 32]);
        reordered[graph_offset..graph_offset + 32]
            .copy_from_slice(&bytes[context_offset..graph_offset]);
        let mut unframed = Vec::with_capacity(bytes.len() - 4);
        unframed.extend_from_slice(&bytes[..12]);
        unframed.extend_from_slice(&bytes[16..]);
        for (field, candidate) in [
            (
                "identity-domain",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.other-check-receipt.v1",
                    &bytes,
                ),
            ),
            (
                "identity-version",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.contract-check-receipt.v2",
                    &bytes,
                ),
            ),
            (
                "transport-magic",
                fs_blake3::hash_domain(CONTRACT_CHECK_RECEIPT_DOMAIN, &wrong_magic),
            ),
            (
                "canonical-field-order",
                fs_blake3::hash_domain(CONTRACT_CHECK_RECEIPT_DOMAIN, &reordered),
            ),
            (
                "length-framing",
                fs_blake3::hash_domain(CONTRACT_CHECK_RECEIPT_DOMAIN, &unframed),
            ),
            (
                "fixed-numeric-little-endian",
                fs_blake3::hash_domain(CONTRACT_CHECK_RECEIPT_DOMAIN, &big_endian),
            ),
        ] {
            assert_protocol_hash_moved(field, base, candidate);
        }
    }

    #[test]
    fn claim_policy_assessment_log_identity_semantic_fields_move_independently() {
        let line = "{\"case\":\"unit-a\"}\n";
        let base = claim_policy_assessment_log_identity(line);
        for (field, candidate) in [
            (
                "identity-domain",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.other-log.v1",
                    line.as_bytes(),
                ),
            ),
            (
                "identity-version",
                fs_blake3::hash_domain(
                    "org.frankensim.fs-euler-disc-e2e.claim-policy-assessment-log.v2",
                    line.as_bytes(),
                ),
            ),
            (
                "exact-json-line-bytes",
                claim_policy_assessment_log_identity("{\"case\":\"unit-b\"}\n"),
            ),
        ] {
            assert_protocol_hash_moved(field, base, candidate);
        }
    }

    #[test]
    fn claim_policy_assessment_log_admits_valid_input_and_preflights_the_byte_ceiling() {
        assert_eq!(
            EULER_CLAIM_REGISTRY,
            EulerClaimKind::ALL,
            "the canonical decoder registry must match the public contract taxonomy"
        );
        assert_eq!(
            EULER_EVIDENCE_REQUIREMENT_REGISTRY,
            EvidenceRequirement::ALL,
            "the canonical decoder evidence registry must match the public contract taxonomy"
        );
        assert_eq!(
            EULER_OWNER_ROLE_REGISTRY,
            OwnerRole::ALL,
            "the canonical owner registry must match the public owner-role taxonomy"
        );
        let contract = build_frozen_contract().expect("frozen contract");
        let packet = test_packet(&contract);
        let mut reasons = packet
            .claim
            .required_evidence()
            .iter()
            .map(|requirement| format!("missing-evidence:{}", requirement.code()))
            .collect::<Vec<_>>();
        reasons.sort();
        let valid = claim_policy_assessment_log(
            &contract,
            &packet,
            &[],
            AssessmentDisposition::Refused,
            &reasons,
            reasons.first().map(String::as_str),
        )
        .expect("closed-grammar refusal log");
        valid
            .verify_identity()
            .expect("a valid retained log must be admitted");
        ClaimPolicyAssessmentLog::from_json_line(valid.json_line.clone())
            .expect("the valid canonical assessment log must decode");
        assert!(valid.json_line.len() < MAX_ASSESSMENT_LOG_BYTES);

        let mut partial_missing_roles = reasons.clone();
        partial_missing_roles
            .pop()
            .expect("multiple required roles");
        let error = claim_policy_assessment_log(
            &contract,
            &packet,
            &[],
            AssessmentDisposition::Refused,
            &partial_missing_roles,
            partial_missing_roles.first().map(String::as_str),
        )
        .expect_err("a partial missing-role diagnostic set must not decode as canonical");
        assert_eq!(error.code(), "EulerProtocolMalformedAssessmentLog");

        let mut false_missing_role = reasons.clone();
        false_missing_role.push("missing-evidence:physical-validation".to_owned());
        false_missing_role.sort();
        let error = claim_policy_assessment_log(
            &contract,
            &packet,
            &[],
            AssessmentDisposition::Refused,
            &false_missing_role,
            false_missing_role.first().map(String::as_str),
        )
        .expect_err("a known role outside the claim policy must not masquerade as missing");
        assert_eq!(error.code(), "EulerProtocolMalformedAssessmentLog");

        let mut oversized_line = valid.json_line.clone();
        oversized_line.insert_str(
            oversized_line.len() - 1,
            &" ".repeat(MAX_ASSESSMENT_LOG_BYTES + 1 - oversized_line.len()),
        );
        assert_eq!(oversized_line.len(), MAX_ASSESSMENT_LOG_BYTES + 1);
        let decode_error = ClaimPolicyAssessmentLog::from_json_line(oversized_line.clone())
            .expect_err("maximum-plus-one retained-log bytes must refuse before hashing");
        assert_eq!(decode_error.code(), "EulerProtocolMalformedAssessmentLog");
        let oversized = ClaimPolicyAssessmentLog {
            identity: claim_policy_assessment_log_identity(&oversized_line),
            json_line: oversized_line,
        };
        let error = oversized
            .verify_identity()
            .expect_err("maximum-plus-one retained-log bytes must refuse");
        assert_eq!(error.code(), "EulerProtocolMalformedAssessmentLog");
    }

    #[test]
    fn assessment_log_artifacts_refuse_a_hypothesis_collision_in_an_unretained_slot() {
        let contract = build_frozen_contract().expect("frozen contract");
        let requirement = EvidenceRequirement::SolutionVerification;
        let packet_identity = test_hash("slot-packet").to_hex();
        let design_set_identity = test_hash("slot-design-set").to_hex();
        let aggregate_qoi_derivation_receipt_identity =
            test_hash("slot-aggregate-qoi-derivation-receipt").to_hex();
        let mut artifacts = vec![
            format!("packet:{packet_identity}"),
            format!("design-set:{design_set_identity}"),
            format!(
                "aggregate-qoi-derivation:{EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA}:{aggregate_qoi_derivation_receipt_identity}"
            ),
            format!(
                "evidence:{}:artifact:{}",
                requirement.code(),
                test_hash("slot-artifact").to_hex()
            ),
            format!(
                "evidence:{}:schema-admission-receipt:{}",
                requirement.code(),
                test_hash("slot-schema-receipt").to_hex()
            ),
        ];
        artifacts.sort();
        let evidence_sources = BTreeSet::from([requirement]);
        let reasons = vec![format!(
            "hypothesis-source-cannot-satisfy-evidence:{}:role-receipt",
            requirement.code()
        )];
        let error = validate_assessment_log_relative_artifacts(
            &artifacts,
            &packet_identity,
            &design_set_identity,
            &aggregate_qoi_derivation_receipt_identity,
            EulerClaimKind::NumericalTrajectoryVerification,
            AssessmentDisposition::Refused.code(),
            &evidence_sources,
            &reasons,
            &contract,
        )
        .expect_err("a collision reason cannot name a role receipt absent from the artifact set");
        assert_eq!(error.code(), "EulerProtocolMalformedAssessmentLog");
    }

    #[test]
    fn assessment_log_artifacts_refuse_multiple_weaknesses_for_one_role() {
        let contract = build_frozen_contract().expect("frozen contract");
        let requirement = EvidenceRequirement::PhysicalValidation;
        let packet_identity = test_hash("weakness-packet").to_hex();
        let design_set_identity = test_hash("weakness-design-set").to_hex();
        let aggregate_qoi_derivation_receipt_identity =
            test_hash("weakness-aggregate-qoi-derivation-receipt").to_hex();
        let mut artifacts = vec![
            format!("packet:{packet_identity}"),
            format!("design-set:{design_set_identity}"),
            format!(
                "aggregate-qoi-derivation:{EULER_AGGREGATE_QOI_DERIVATION_RECEIPT_SCHEMA}:{aggregate_qoi_derivation_receipt_identity}"
            ),
            format!(
                "evidence:{}:artifact:{}",
                requirement.code(),
                test_hash("weakness-artifact").to_hex()
            ),
            format!(
                "evidence:{}:role-receipt:{}",
                requirement.code(),
                test_hash("weakness-role-receipt").to_hex()
            ),
            format!(
                "evidence:{}:schema-admission-receipt:{}",
                requirement.code(),
                test_hash("weakness-schema-receipt").to_hex()
            ),
        ];
        artifacts.sort();
        let evidence_sources = BTreeSet::from([requirement]);
        let mut reasons = vec![
            "weak-authority:physical-validation:requires-validated-physical:observed-structural-process"
                .to_owned(),
            "weak-independence:physical-validation:independent-evidence-required".to_owned(),
        ];
        reasons.sort();
        let error = validate_assessment_log_relative_artifacts(
            &artifacts,
            &packet_identity,
            &design_set_identity,
            &aggregate_qoi_derivation_receipt_identity,
            EulerClaimKind::BlindTrajectoryPrediction,
            AssessmentDisposition::Refused.code(),
            &evidence_sources,
            &reasons,
            &contract,
        )
        .expect_err("one evidence row cannot carry two evaluator-exclusive weakness reasons");
        assert_eq!(error.code(), "EulerProtocolMalformedAssessmentLog");
    }

    #[test]
    fn malformed_prerequisite_permutations_preserve_first_divergence_and_identities() {
        let contract = build_frozen_contract().expect("frozen contract");
        let packet = test_packet(&contract);
        let source_point_bytes = applicability_point_bytes(&packet.point).expect("point bytes");
        let mut receipt = PrerequisiteAssessmentReceipt {
            schema_version: EULER_PROTOCOL_SCHEMA_VERSION,
            contract_identity: contract.identity(),
            prerequisite: EulerClaimKind::CalibratedReproduction,
            dependent: packet.claim,
            use_kind: EvidenceUse::ValidationInput,
            source_packet_identity: test_hash("prerequisite-source-packet"),
            source_assessment_identity: test_hash("prerequisite-source-assessment"),
            source_design_set_identity: packet.design_set_identity,
            source_point_bytes,
            identity: test_hash("temporary-prerequisite-identity"),
        };
        receipt.identity = prerequisite_assessment_receipt_identity(
            &receipt.canonical_bytes().expect("receipt bytes"),
        );
        receipt.verify().expect("baseline receipt identity");

        let mut malformed_identity = receipt.clone();
        malformed_identity.identity = test_hash("corrupt-prerequisite-identity");
        assert!(malformed_identity.verify().is_err());

        let mut stale_schema = receipt;
        stale_schema.schema_version += 1;
        stale_schema.identity = prerequisite_assessment_receipt_identity(
            &stale_schema.canonical_bytes().expect("stale receipt bytes"),
        );
        stale_schema
            .verify()
            .expect("self-consistent stale receipt identity");

        let forward = assess_packet(
            &contract,
            &packet,
            &[malformed_identity.clone(), stale_schema.clone()],
        )
        .expect("forward malformed receipt assessment");
        let reverse = assess_packet(&contract, &packet, &[stale_schema, malformed_identity])
            .expect("reverse malformed receipt assessment");
        assert_eq!(
            forward, reverse,
            "malformed receipt order must not alter first divergence or retained identities"
        );
        assert!(
            forward
                .reasons
                .iter()
                .any(|reason| { reason.starts_with("malformed-prerequisite-receipt:") })
        );
        assert!(
            forward
                .reasons
                .iter()
                .any(|reason| reason.starts_with("stale-prerequisite-receipt:"))
        );
        forward.verify_identity().expect("forward identity");
        reverse.verify_identity().expect("reverse identity");
        assert_eq!(forward.log.identity, reverse.log.identity);
        assert_eq!(forward.identity, reverse.identity);
    }

    #[test]
    fn prerequisite_design_set_mismatch_refuses_and_is_retained_by_the_strict_log() {
        let contract = build_frozen_contract().expect("frozen contract");
        let target = ClaimEvidencePacket::try_new(
            contract.identity(),
            "prerequisite-design-set-mismatch",
            test_hash("target-design-set"),
            test_hash("target-aggregate-qoi-derivation-receipt"),
            EulerClaimKind::CalibratedReproduction,
            test_point(&contract, 0.5),
            Vec::new(),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
            vec!["1".to_owned()],
            ProtocolSeed::Fixed { value: 23 },
            ProtocolBudget::try_new(100, 200, 0.25).expect("test budget"),
        )
        .expect("target packet");
        let source_point_bytes =
            applicability_point_bytes(target.point()).expect("target point bytes");
        let mut receipt = PrerequisiteAssessmentReceipt {
            schema_version: EULER_PROTOCOL_SCHEMA_VERSION,
            contract_identity: contract.identity(),
            prerequisite: EulerClaimKind::NumericalTrajectoryVerification,
            dependent: EulerClaimKind::CalibratedReproduction,
            use_kind: EvidenceUse::ValidationInput,
            source_packet_identity: test_hash("mismatched-design-source-packet"),
            source_assessment_identity: test_hash("mismatched-design-source-assessment"),
            source_design_set_identity: test_hash("different-design-set"),
            source_point_bytes,
            identity: test_hash("temporary-mismatched-design-receipt"),
        };
        receipt.identity = prerequisite_assessment_receipt_identity(
            &receipt.canonical_bytes().expect("receipt bytes"),
        );

        let assessment = assess_packet(&contract, &target, &[receipt])
            .expect("design-set mismatch must remain representable");
        assert_eq!(assessment.disposition, AssessmentDisposition::Refused);
        assert!(assessment.reasons.iter().any(|reason| {
            reason == "prerequisite-design-set-mismatch:numerical-trajectory-verification"
        }));
        assert!(assessment.reasons.iter().any(|reason| {
            reason
                == "missing-prerequisite-receipt:numerical-trajectory-verification:validation-input"
        }));
        assessment
            .log
            .verify_identity()
            .expect("strict log must retain the design-set mismatch grammar");
        assessment.verify_identity().expect("assessment identity");
    }

    #[test]
    fn exact_and_malformed_prerequisites_sharing_one_retained_identity_do_not_invent_missing_edge()
    {
        let contract = build_frozen_contract().expect("frozen contract");
        let target = ClaimEvidencePacket::try_new(
            contract.identity(),
            "deduplicated-prerequisite-artifact",
            test_hash("deduplicated-design-set"),
            test_hash("deduplicated-aggregate-qoi-derivation-receipt"),
            EulerClaimKind::CalibratedReproduction,
            test_point(&contract, 0.5),
            Vec::new(),
            true,
            false,
            ReportedScientificDisposition::Positive,
            AssessmentDisposition::Refused,
            vec!["1".to_owned()],
            ProtocolSeed::Fixed { value: 19 },
            ProtocolBudget::try_new(100, 200, 0.25).expect("test budget"),
        )
        .expect("target packet");
        let source_point_bytes =
            applicability_point_bytes(target.point()).expect("target point bytes");
        let mut exact = PrerequisiteAssessmentReceipt {
            schema_version: EULER_PROTOCOL_SCHEMA_VERSION,
            contract_identity: contract.identity(),
            prerequisite: EulerClaimKind::NumericalTrajectoryVerification,
            dependent: EulerClaimKind::CalibratedReproduction,
            use_kind: EvidenceUse::ValidationInput,
            source_packet_identity: test_hash("deduplicated-source-packet"),
            source_assessment_identity: test_hash("deduplicated-source-assessment"),
            source_design_set_identity: target.design_set_identity,
            source_point_bytes,
            identity: test_hash("temporary-deduplicated-identity"),
        };
        exact.identity = prerequisite_assessment_receipt_identity(
            &exact.canonical_bytes().expect("exact receipt bytes"),
        );
        exact.verify().expect("exact receipt");

        let mut malformed = exact.clone();
        malformed.source_assessment_identity = test_hash("transplanted-source-assessment");
        assert_eq!(malformed.identity, exact.identity);
        assert!(malformed.verify().is_err());

        let assessment = assess_packet(&contract, &target, &[exact, malformed])
            .expect("deduplicated exact plus malformed receipts remain representable");
        assert!(assessment.reasons.iter().any(|reason| {
            reason
                == "malformed-prerequisite-receipt:EulerProtocolPrerequisiteReceiptIdentityMismatch"
        }));
        assert!(!assessment.reasons.iter().any(|reason| {
            reason.starts_with("missing-prerequisite-receipt:numerical-trajectory-verification:")
        }));
        assessment.verify_identity().expect("assessment identity");
        assessment
            .log
            .verify_identity()
            .expect("assessment log identity");
    }

    #[test]
    fn frozen_contract_is_whole_transport_fixed_point() {
        let contract = build_frozen_contract().expect("frozen contract");
        let bytes = contract.canonical_bytes().expect("contract bytes");
        let decoded =
            EulerScientificContract::from_canonical_bytes(&bytes).expect("decode fixed point");
        assert_eq!(decoded, contract);
        assert_eq!(decoded.canonical_bytes().expect("re-encode"), bytes);
    }

    #[test]
    fn independent_check_receipt_verifies_exact_subject() {
        let contract = build_frozen_contract().expect("frozen contract");
        let receipt = check_frozen_contract(&contract).expect("check receipt");
        assert!(receipt.passed(), "{:?}", receipt.issues());
        receipt.verify_subject(&contract).expect("receipt binding");
        assert_ne!(receipt.identity(), contract.identity().as_hash());
    }
}
