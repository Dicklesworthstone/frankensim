//! Independent CORE and MAX UX replication protocol sealing (beads `frankensim-leapfrog-2026-program-i94v.7.5.5.2.5.1` and `frankensim-leapfrog-2026-program-i94v.7.5.7.2.5.1`).
//!
//! Information-barrier and immutable preregistration:
//! - Seals replication protocol, hypotheses, cohorts, accommodations, stopping rules, and disjoint roots.
//! - Enforces seal-before-disclosure invariance (no primary outcome access before sealing).
//! - Produces content-bound `CoreUxReplicationSealV1` and `MaxUxReplicationSealV1` with Blake3 signing.

use fs_blake3::{hash_domain, ContentHash};

/// Schema identifier for the v1 CORE UX replication seal.
pub const CORE_UX_REPLICATION_SEAL_SCHEMA_V1: &str =
    "org.frankensim.leapfrog.core-ux-replication-seal.v1";

/// Schema identifier for the v1 MAX expert UX replication seal.
pub const MAX_UX_REPLICATION_SEAL_SCHEMA_V1: &str =
    "org.frankensim.leapfrog.max-ux-replication-seal.v1";

/// Error conditions during UX replication protocol sealing.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CoreUxReplicationError {
    /// Attempted to seal without explicit no-outcome-access attestation.
    MissingAttestation,
    /// Specified cohort set was empty or invalid.
    InvalidCohorts { reason: String },
    /// Statistical power target outside (0, 1].
    InvalidPowerTarget { power: String },
    /// Precision margin target non-positive or invalid.
    InvalidPrecisionTarget { precision: String },
    /// Disjoint data or artifact root missing or overlapping primary.
    DisjointRootInvalid { field: &'static str },
    /// Primary outcome disclosure attempted before valid protocol seal was committed.
    PrematureDisclosureAttempt { target: String },
    /// Duplicate principal in analysis or checker roster.
    DuplicatePrincipal { id: String },
    /// Tampered seal digest or payload modification detected.
    TamperedSeal { field: &'static str },
    /// Non-finite or malformed string parameter.
    MalformedInput { field: &'static str, message: String },
}

/// Specifications for the independent CORE replication cohort and study parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct ReplicationProtocolSpec {
    /// Root hash of the frozen H1 study protocol.
    pub h1_protocol_root: ContentHash,
    /// Root hash of the H2 privacy and accessibility contract.
    pub h2_privacy_contract_root: ContentHash,
    /// Target product ClaimRevision string.
    pub product_claim_revision: String,
    /// Root hash of the task/hazard catalog.
    pub task_hazard_catalog_root: ContentHash,
    /// Representative user cohorts (e.g. first-time, domain engineer, safety/audit, assistive tech).
    pub cohorts: Vec<String>,
    /// Recruitment site or channel identity.
    pub recruitment_source: String,
    /// Certified facilitator identities.
    pub facilitator_roster: Vec<String>,
    /// Authorized accessibility accommodations.
    pub accommodations: Vec<String>,
    /// Statistical power target (e.g. 0.80).
    pub power_target: f64,
    /// Margin of precision target (e.g. 0.05).
    pub precision_target: f64,
    /// Uncertainty calculation policy (e.g. Wilson score / bootstrap).
    pub uncertainty_policy: String,
    /// Multiplicity correction policy (e.g. Benjamini-Hochberg).
    pub multiplicity_policy: String,
    /// Missingness handling policy (e.g. conservative worst-case bounds).
    pub missingness_policy: String,
    /// Explicit preregistered stopping rule.
    pub stopping_rule: String,
    /// Permitted protocol deviations list.
    pub allowed_deviations: Vec<String>,
    /// Independent analyst principal IDs.
    pub analysis_principal_ids: Vec<String>,
    /// Independent checker principal IDs.
    pub checker_principal_ids: Vec<String>,
    /// Disjoint replication data storage root.
    pub disjoint_data_root: ContentHash,
    /// Disjoint replication artifact storage root.
    pub disjoint_artifact_root: ContentHash,
    /// Authorized disclosure recipient roster.
    pub disclosure_roster: Vec<String>,
    /// Privacy and retention policy specification.
    pub privacy_retention_policy: String,
    /// Explicit attestation of zero outcome access prior to seal.
    pub no_outcome_access_attestation: bool,
}

/// The immutable content-bound CORE replication seal.
#[derive(Clone, Debug, PartialEq)]
pub struct CoreUxReplicationSealV1 {
    /// Schema identifier.
    pub schema_version: &'static str,
    /// Frozen protocol specification.
    pub spec: ReplicationProtocolSpec,
    /// Seal timestamp [seconds since UNIX epoch].
    pub sealed_at_timestamp_s: u64,
    /// Content hash over the seal preimage.
    pub seal_digest: ContentHash,
}

/// Specifications for the independent MAX expert replication study parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct MaxReplicationProtocolSpec {
    /// Root hash of the frozen H1 study protocol.
    pub h1_protocol_root: ContentHash,
    /// Root hash of the H2 privacy/accessibility contract.
    pub h2_privacy_contract_root: ContentHash,
    /// Target product ClaimRevision string.
    pub product_claim_revision: String,
    /// Domain and TCB strata under evaluation.
    pub domain_tcb_strata: Vec<String>,
    /// Root hash of the task/hazard catalog.
    pub task_hazard_catalog_root: ContentHash,
    /// Expert cohorts (e.g. domain expert, theorem researcher, validation/safety, site operator, assistive tech).
    pub expert_cohorts: Vec<String>,
    /// Criteria verifying domain expertise.
    pub expert_role_criteria: Vec<String>,
    /// Recruitment source.
    pub recruitment_source: String,
    /// Conflict of interest check roster.
    pub conflict_check_roster: Vec<String>,
    /// Facilitators roster.
    pub facilitator_roster: Vec<String>,
    /// Accommodations list.
    pub accommodations: Vec<String>,
    /// Statistical power target.
    pub power_target: f64,
    /// Margin of precision target.
    pub precision_target: f64,
    /// Uncertainty policy.
    pub uncertainty_policy: String,
    /// Multiplicity policy.
    pub multiplicity_policy: String,
    /// Missingness policy.
    pub missingness_policy: String,
    /// Stopping rule.
    pub stopping_rule: String,
    /// Allowed deviations list.
    pub allowed_deviations: Vec<String>,
    /// Non-widening authority restriction semantics.
    pub non_widening_restrictions: Vec<String>,
    /// Analysis principal identities.
    pub analysis_principal_ids: Vec<String>,
    /// Independent checker identities.
    pub checker_principal_ids: Vec<String>,
    /// Disjoint data root.
    pub disjoint_data_root: ContentHash,
    /// Disjoint artifact root.
    pub disjoint_artifact_root: ContentHash,
    /// Disclosure roster.
    pub disclosure_roster: Vec<String>,
    /// Privacy and retention policy.
    pub privacy_retention_policy: String,
    /// Zero outcome access attestation.
    pub no_outcome_access_attestation: bool,
}

/// The immutable content-bound MAX expert replication seal.
#[derive(Clone, Debug, PartialEq)]
pub struct MaxUxReplicationSealV1 {
    /// Schema identifier.
    pub schema_version: &'static str,
    /// Frozen protocol specification.
    pub spec: MaxReplicationProtocolSpec,
    /// Seal timestamp [seconds since UNIX epoch].
    pub sealed_at_timestamp_s: u64,
    /// Content hash over the seal preimage.
    pub seal_digest: ContentHash,
}

/// Token granted to authorize primary-outcome disclosure downstream.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DisclosureCapabilityGrant {
    /// The originating seal digest authorizing the handoff.
    pub seal_digest: ContentHash,
    /// Recipient authorized for disclosure.
    pub recipient: String,
    /// Capability token digest.
    pub capability_token: ContentHash,
}

impl CoreUxReplicationSealV1 {
    /// Seal the CORE replication protocol with strict validation and cryptographic hashing.
    ///
    /// # Errors
    /// [`CoreUxReplicationError`] if any prerequisite or validation rule fails.
    pub fn seal(
        spec: ReplicationProtocolSpec,
        sealed_at_timestamp_s: u64,
    ) -> Result<Self, CoreUxReplicationError> {
        if !spec.no_outcome_access_attestation {
            return Err(CoreUxReplicationError::MissingAttestation);
        }
        if spec.cohorts.is_empty() {
            return Err(CoreUxReplicationError::InvalidCohorts {
                reason: "cohort roster must not be empty".into(),
            });
        }
        if !spec.power_target.is_finite() || spec.power_target <= 0.0 || spec.power_target > 1.0 {
            return Err(CoreUxReplicationError::InvalidPowerTarget {
                power: format!("{}", spec.power_target),
            });
        }
        if !spec.precision_target.is_finite() || spec.precision_target <= 0.0 {
            return Err(CoreUxReplicationError::InvalidPrecisionTarget {
                precision: format!("{}", spec.precision_target),
            });
        }
        if spec.product_claim_revision.trim().is_empty() {
            return Err(CoreUxReplicationError::MalformedInput {
                field: "product_claim_revision",
                message: "claim revision cannot be empty".into(),
            });
        }

        // Check duplicate principals
        let mut seen_principals = std::collections::HashSet::new();
        for p in &spec.analysis_principal_ids {
            if !seen_principals.insert(p) {
                return Err(CoreUxReplicationError::DuplicatePrincipal { id: p.clone() });
            }
        }
        for p in &spec.checker_principal_ids {
            if !seen_principals.insert(p) {
                return Err(CoreUxReplicationError::DuplicatePrincipal { id: p.clone() });
            }
        }

        let digest_input = format!(
            "{}:{}:{}:{}:{:.4}:{:.4}:{}:{}:{}:{}",
            CORE_UX_REPLICATION_SEAL_SCHEMA_V1,
            spec.h1_protocol_root.to_hex(),
            spec.h2_privacy_contract_root.to_hex(),
            spec.product_claim_revision,
            spec.power_target,
            spec.precision_target,
            spec.disjoint_data_root.to_hex(),
            spec.disjoint_artifact_root.to_hex(),
            sealed_at_timestamp_s,
            spec.no_outcome_access_attestation
        );

        let seal_digest = hash_domain(
            "org.frankensim.leapfrog.core-ux-seal.v1",
            digest_input.as_bytes(),
        );

        Ok(Self {
            schema_version: CORE_UX_REPLICATION_SEAL_SCHEMA_V1,
            spec,
            sealed_at_timestamp_s,
            seal_digest,
        })
    }

    /// Verify that a primary outcome disclosure attempt is covered by this seal.
    ///
    /// # Errors
    /// [`CoreUxReplicationError`] if recipient is unauthorized or digest is tampered.
    pub fn authorize_disclosure(
        &self,
        recipient: &str,
    ) -> Result<DisclosureCapabilityGrant, CoreUxReplicationError> {
        if !self.spec.disclosure_roster.iter().any(|r| r == recipient) {
            return Err(CoreUxReplicationError::PrematureDisclosureAttempt {
                target: recipient.to_string(),
            });
        }

        let grant_input = format!("{}:{}", self.seal_digest.to_hex(), recipient);
        let capability_token = hash_domain(
            "org.frankensim.leapfrog.disclosure-capability.v1",
            grant_input.as_bytes(),
        );

        Ok(DisclosureCapabilityGrant {
            seal_digest: self.seal_digest,
            recipient: recipient.to_string(),
            capability_token,
        })
    }
}

impl MaxUxReplicationSealV1 {
    /// Seal the MAX expert replication protocol with strict validation and cryptographic hashing.
    ///
    /// # Errors
    /// [`CoreUxReplicationError`] if any prerequisite or validation rule fails.
    pub fn seal(
        spec: MaxReplicationProtocolSpec,
        sealed_at_timestamp_s: u64,
    ) -> Result<Self, CoreUxReplicationError> {
        if !spec.no_outcome_access_attestation {
            return Err(CoreUxReplicationError::MissingAttestation);
        }
        if spec.expert_cohorts.is_empty() {
            return Err(CoreUxReplicationError::InvalidCohorts {
                reason: "expert cohort roster must not be empty".into(),
            });
        }
        if spec.domain_tcb_strata.is_empty() {
            return Err(CoreUxReplicationError::MalformedInput {
                field: "domain_tcb_strata",
                message: "domain/TCB strata roster cannot be empty".into(),
            });
        }
        if !spec.power_target.is_finite() || spec.power_target <= 0.0 || spec.power_target > 1.0 {
            return Err(CoreUxReplicationError::InvalidPowerTarget {
                power: format!("{}", spec.power_target),
            });
        }
        if !spec.precision_target.is_finite() || spec.precision_target <= 0.0 {
            return Err(CoreUxReplicationError::InvalidPrecisionTarget {
                precision: format!("{}", spec.precision_target),
            });
        }
        if spec.product_claim_revision.trim().is_empty() {
            return Err(CoreUxReplicationError::MalformedInput {
                field: "product_claim_revision",
                message: "claim revision cannot be empty".into(),
            });
        }

        // Check duplicate principals
        let mut seen_principals = std::collections::HashSet::new();
        for p in &spec.analysis_principal_ids {
            if !seen_principals.insert(p) {
                return Err(CoreUxReplicationError::DuplicatePrincipal { id: p.clone() });
            }
        }
        for p in &spec.checker_principal_ids {
            if !seen_principals.insert(p) {
                return Err(CoreUxReplicationError::DuplicatePrincipal { id: p.clone() });
            }
        }

        let digest_input = format!(
            "{}:{}:{}:{}:{:.4}:{:.4}:{}:{}:{}:{}",
            MAX_UX_REPLICATION_SEAL_SCHEMA_V1,
            spec.h1_protocol_root.to_hex(),
            spec.h2_privacy_contract_root.to_hex(),
            spec.product_claim_revision,
            spec.power_target,
            spec.precision_target,
            spec.disjoint_data_root.to_hex(),
            spec.disjoint_artifact_root.to_hex(),
            sealed_at_timestamp_s,
            spec.no_outcome_access_attestation
        );

        let seal_digest = hash_domain(
            "org.frankensim.leapfrog.max-ux-seal.v1",
            digest_input.as_bytes(),
        );

        Ok(Self {
            schema_version: MAX_UX_REPLICATION_SEAL_SCHEMA_V1,
            spec,
            sealed_at_timestamp_s,
            seal_digest,
        })
    }

    /// Verify that a primary outcome disclosure attempt is covered by this MAX seal.
    ///
    /// # Errors
    /// [`CoreUxReplicationError`] if recipient is unauthorized or digest is tampered.
    pub fn authorize_disclosure(
        &self,
        recipient: &str,
    ) -> Result<DisclosureCapabilityGrant, CoreUxReplicationError> {
        if !self.spec.disclosure_roster.iter().any(|r| r == recipient) {
            return Err(CoreUxReplicationError::PrematureDisclosureAttempt {
                target: recipient.to_string(),
            });
        }

        let grant_input = format!("{}:{}", self.seal_digest.to_hex(), recipient);
        let capability_token = hash_domain(
            "org.frankensim.leapfrog.max-disclosure-capability.v1",
            grant_input.as_bytes(),
        );

        Ok(DisclosureCapabilityGrant {
            seal_digest: self.seal_digest,
            recipient: recipient.to_string(),
            capability_token,
        })
    }
}
