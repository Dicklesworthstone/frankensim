//! Euler-disc flagship integration boundary.
//!
//! The crate freezes the scientific Context of Use, claim taxonomy, evidence
//! ceilings, and binding no-claims. Its only executable physics slice is a
//! bounded ideal no-slip rolling baseline with an explicit thin-disc
//! small-angle oracle; it does not implement compliant contact, base, gas,
//! experiment, or target-outcome prediction.

#![forbid(unsafe_code)]

pub mod baseline;
pub mod contract;
pub mod protocol;

pub use baseline::{
    BaselineDynamicsClass, BaselineEnergyLedger, BaselineEquilibriumReceipt, BaselineRefusal,
    BaselineRefusalReason, BaselineRunOutput, BaselineSample, BaselineState,
    BaselineSupportDiagnostic, BaselineTerminal, BaselineTrajectory, SquatDiscInput,
    run_ideal_conservative_baseline,
};

pub use contract::{
    AuthorityCeiling, CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN, CONTRACT_CHECK_RECEIPT_DOMAIN,
    CORE_NO_CLAIMS, ContractError, ContractIdentity, EULER_ASSESSMENT_IDENTITY_DOMAIN,
    EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, EULER_CLAIM_POLICY_SCHEMA_VERSION,
    EULER_CONTRACT_IDENTITY_DOMAIN, EULER_CONTRACT_SCHEMA_VERSION,
    EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, EULER_OWNER_MATRIX_IDENTITY_DOMAIN,
    EULER_OWNER_MATRIX_SCHEMA_VERSION, EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,
    EulerAcceptanceFamily, EulerClaimGraph, EulerClaimKind, EulerClaimSpec, EulerContextExtension,
    EulerScientificContract, EvidenceRequirement, FS_GOVERN_AUTHORITY_SOURCE_SCHEMA,
    FS_IR_CAMPAIGN_SOURCE_SCHEMA, HYPOTHESIS_SOURCE_DECLARATION_DOMAIN, HypothesisSource,
    MAX_EULER_CLAIMS, MAX_EULER_NO_CLAIMS, MAX_OWNER_MATRIX_BYTES, OwnerMatrix,
    OwnerMatrixIdentity, OwnerRole, OwnerRow, ScientificRisk,
};
pub use protocol::{
    AssessmentDisposition, ClaimEvidencePacket, ClaimPolicyAssessment, ClaimPolicyAssessmentLog,
    ContractCheckReceipt, DeclaredEvidenceAccessClass, EULER_PROTOCOL_SCHEMA_VERSION,
    EvidenceAuthorityClass, EvidenceAuthorityDeclaration, EvidenceRecord,
    FROZEN_CLAIM_GRAPH_HASH_HEX, FROZEN_CONTEXT_HASH_HEX, FROZEN_CONTRACT_IDENTITY_HEX,
    MAX_PREREQUISITE_RECEIPTS, MAX_PROTOCOL_ID_BYTES, MAX_VALIDITY_DOMAIN_AXES,
    MAX_VALIDITY_DOMAIN_CANONICAL_BYTES, PrerequisiteAssessmentReceipt, ProtocolBudget,
    ProtocolSeed, ReportedScientificDisposition, StructurallyAdmittedEulerContract,
    admit_frozen_contract, build_frozen_contract, check_frozen_contract,
};
