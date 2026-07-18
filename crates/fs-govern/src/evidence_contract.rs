//! Phase 0B-A evidence-contract algebra.
//!
//! This module owns the small, pure vocabulary shared by later schema,
//! checker, graph, ledger, and runtime-admission layers.  It deliberately
//! does not read files, authenticate signatures, or persist receipts.  Its
//! job is narrower and load-bearing: make semantic identity canonical, keep
//! truth and admission axes separate, and make authority-bearing objects
//! impossible to construct by filling public fields.
//!
//! Public values in this Phase 0B-A module are descriptive candidates.  Raw
//! hashes and enum selections cannot mint positive or negative authority:
//! grants, authenticated checker decisions, live heads, and runtime admissions
//! have no public widening constructor until Phase 0B-B verifies durable
//! receipts.
//!
//! ```compile_fail
//! use fs_govern::evidence_contract::{
//!     AuthorityGrant, AuthorityState, ProvedAuthority,
//! };
//!
//! fn cannot_widen(state: AuthorityState) -> AuthorityGrant<ProvedAuthority> {
//!     state.into()
//! }
//! ```
//!
//! ```compile_fail
//! use fs_govern::evidence_contract::{RuntimeAdmission, RuntimeAssessment};
//!
//! fn assessment_is_not_authority(candidate: RuntimeAssessment) -> RuntimeAdmission {
//!     candidate.into()
//! }
//! ```
//!
//! ```compile_fail
//! use fs_govern::evidence_contract::{NonvacuousEvidence, SatisfiableEvidence};
//!
//! fn cannot_substitute(satisfiable: SatisfiableEvidence) -> NonvacuousEvidence {
//!     satisfiable.into()
//! }
//! ```
//!
//! ```compile_fail
//! use fs_govern::evidence_contract::{SatisfiabilityState, UnsatisfiableEvidence};
//!
//! fn opposite_polarity_is_not_substitutable(evidence: UnsatisfiableEvidence) -> SatisfiabilityState {
//!     SatisfiabilityState::Satisfiable(evidence)
//! }
//! ```
//!
//! ```compile_fail
//! use fs_govern::evidence_contract::{NonvacuityState, VacuousEvidence};
//!
//! fn vacuous_is_not_nonvacuous(evidence: VacuousEvidence) -> NonvacuityState {
//!     NonvacuityState::Nonvacuous(evidence)
//! }
//! ```
//!
//! ```compile_fail
//! use fs_govern::evidence_contract::{ReproductionFailedEvidence, ReproductionState};
//!
//! fn failure_is_not_reproduction(evidence: ReproductionFailedEvidence) -> ReproductionState {
//!     ReproductionState::Reproduced(evidence)
//! }
//! ```
//!
//! ```compile_fail
//! use fs_govern::evidence_contract::{InvalidationBinding, InvalidationState};
//!
//! fn raw_fields_cannot_mint_invalidation(binding: InvalidationBinding) -> InvalidationState {
//!     InvalidationState { binding: Some(binding) }
//! }
//! ```

// Public vocabulary is exhaustively documented by the code-owned catalog and
// the drift-checked fs-govern contract; repeating prose on every enum arm and
// typed-ID expansion would make those two truth sources easier to diverge.
#![allow(missing_docs)]
// The private authority-boundary tests stay adjacent to the private grant
// constructors they exercise; later public candidate types remain below.
#![allow(clippy::items_after_test_module)]

use crate::{
    json_escape,
    lanes::{LaneCharter, PROOF_LANE_IDENTITY_DOMAIN, ProofLaneId},
};
use fs_blake3::ContentHash;
use std::{collections::BTreeMap, marker::PhantomData};

/// Current semantic algebra.  Wire/schema work in Phase 0B-B must bind this
/// value instead of inferring meaning from enum ordinals.
pub const AUTHORITY_ALGEBRA_VERSION: u32 = 2;

/// Current default policy vocabulary.  A policy identity also binds all of
/// its explicit requirements, so changing a guard changes the identity.
pub const AUTHORITY_POLICY_VERSION: u32 = 2;

/// Historical baseline accepted by [`migrate_legacy_v0`].
pub const LEGACY_AUTHORITY_SCHEMA_VERSION: u32 = 0;

/// Persisted v1 identities are semantically ambiguous under the v2 quantifier
/// and polarity rules and therefore have no compatibility reinterpretation.
pub const RETIRED_AUTHORITY_SCHEMA_VERSION: u32 = 1;

/// Maximum bytes in one canonical scalar text field before allocation-heavy
/// canonicalization.
pub const MAX_AUTHORITY_TEXT_BYTES: usize = 4096;

/// Maximum members in any canonical set-like field.
pub const MAX_AUTHORITY_SET_MEMBERS: usize = 256;

/// Maximum rendered diagnostic bytes.  Logging is bounded data, not an
/// unbounded dump of attacker-controlled claim text.
pub const MAX_AUTHORITY_LOG_BYTES: usize = 16 * 1024;

/// Canonical identity domains.  These are public so schema/catalog tooling can
/// drift-check the exact namespaces without duplicating strings.
pub const CLAIM_STATEMENT_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.claim-statement.v2";
pub const QUANTIFIED_DOMAIN_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.quantified-domain.v2";
pub const ASSUMPTION_SET_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.assumption-set.v2";
pub const SEMANTIC_CLAIM_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.semantic-claim.v2";
pub const CLAIM_LANE_BINDING_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.claim-lane-binding.v2";
pub const CLAIM_INSTANCE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.claim-instance.v2";
pub const EVIDENCE_REF_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.evidence-ref.v2";
pub const EVIDENCE_STATE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.evidence-state.v2";
pub const SATISFIABLE_EVIDENCE_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.satisfiable-evidence.v2";
pub const UNSATISFIABLE_EVIDENCE_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.unsatisfiable-evidence.v2";
pub const NONVACUOUS_EVIDENCE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.nonvacuous-evidence.v2";
pub const VACUOUS_EVIDENCE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.vacuous-evidence.v2";
pub const REPRODUCTION_FAILED_EVIDENCE_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.reproduction-failed-evidence.v2";
pub const REPRODUCED_EVIDENCE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.reproduced-evidence.v2";
pub const EXACT_INSTANCE_DECISION_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.exact-instance-decision.v2";
pub const INVALIDATION_BINDING_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.invalidation-binding.v2";
pub const AUTHORITY_STATE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.authority-state.v2";
pub const SUPPORT_EDGE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.support-edge.v2";
pub const ATTACK_EDGE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.attack-edge.v2";
pub const COUNTEREXAMPLE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.counterexample.v2";
pub const ADJUDICATION_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.counterexample-adjudication.v2";
pub const REVOCATION_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.revocation-tombstone.v2";
pub const VERIFIED_REVOCATION_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.verified-revocation-tombstone.v2";
pub const CAPABILITY_POLICY_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.capability-policy.v2";
pub const CHECKER_DECISION_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.checker-decision.v2";
pub const AUTHORITY_HEAD_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.authority-head.v2";
pub const RUNTIME_ADMISSION_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.runtime-admission.v2";
pub const INFERENCE_RULE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.inference-rule.v2";
pub const AUTHORITY_MIGRATION_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.authority-migration.v2";

/// Structured refusal for every contract constructor and transition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AuthorityError {
    EmptyField {
        what: &'static str,
    },
    TooLarge {
        what: &'static str,
        observed: usize,
        cap: usize,
    },
    DuplicateMember {
        what: &'static str,
        key: String,
    },
    InvalidValue {
        what: &'static str,
    },
    MissingIdentity {
        what: &'static str,
    },
    SchemaVersionRefused {
        observed: u32,
        supported: u32,
    },
    CatalogSchemaVersionRefused {
        observed: u32,
        supported: u32,
    },
    MigrationUnavailable {
        observed: u32,
        target: u32,
    },
    IdentityMismatch {
        what: &'static str,
    },
    IncompatibleAxes {
        what: &'static str,
    },
    IllegalTransition {
        from: &'static str,
        to: &'static str,
    },
    TerminalState {
        state: &'static str,
    },
    CompositionConflict {
        axis: &'static str,
    },
    CapabilityMissing {
        capability: String,
        version: u32,
    },
    AssumptionNotAccepted {
        assumption: String,
    },
    NoClaimNotAccepted {
        boundary: String,
    },
    CheckerRefused {
        verdict: CheckerVerdict,
    },
    RuntimeRequirementNotMet {
        requirement: &'static str,
    },
    AdjudicationNotRevocable,
    LogCapacityExceeded {
        needed: usize,
        cap: usize,
    },
}

impl AuthorityError {
    /// Stable ranked remediation code for agent-facing diagnostics.
    #[must_use]
    pub fn remedy_code(&self) -> &'static str {
        match self {
            Self::EmptyField { .. } => "supply-required-field",
            Self::TooLarge { .. } | Self::LogCapacityExceeded { .. } => "reduce-bounded-input",
            Self::DuplicateMember { .. } => "deduplicate-semantic-member",
            Self::InvalidValue { .. } => "repair-invalid-value",
            Self::MissingIdentity { .. } => "supply-content-identity",
            Self::SchemaVersionRefused { .. } => "use-supported-schema-or-explicit-migration",
            Self::CatalogSchemaVersionRefused { .. } => {
                "regenerate-authority-catalog-under-supported-schema"
            }
            Self::MigrationUnavailable { .. } => "regenerate-under-current-authority-algebra",
            Self::IdentityMismatch { .. } => "bind-exact-source-identity",
            Self::IncompatibleAxes { .. } | Self::CompositionConflict { .. } => {
                "adjudicate-conflicting-authority"
            }
            Self::IllegalTransition { .. } | Self::TerminalState { .. } => {
                "open-versioned-successor"
            }
            Self::CapabilityMissing { .. } => "supply-required-capability",
            Self::AssumptionNotAccepted { .. } => "accept-or-discharge-assumption",
            Self::NoClaimNotAccepted { .. } => "accept-or-narrow-no-claim-boundary",
            Self::CheckerRefused { .. } => "resolve-checker-decision",
            Self::RuntimeRequirementNotMet { .. } => "satisfy-runtime-policy-axis",
            Self::AdjudicationNotRevocable => "adjudicate-genuine-counterexample-first",
        }
    }
}

impl core::fmt::Display for AuthorityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyField { what } => write!(f, "authority field `{what}` is empty"),
            Self::TooLarge {
                what,
                observed,
                cap,
            } => write!(f, "authority field `{what}` has size {observed}, cap {cap}"),
            Self::DuplicateMember { what, key } => {
                write!(f, "duplicate {what} member `{key}`")
            }
            Self::InvalidValue { what } => write!(f, "invalid authority value: {what}"),
            Self::MissingIdentity { what } => write!(f, "missing identity for {what}"),
            Self::SchemaVersionRefused {
                observed,
                supported,
            } => write!(
                f,
                "authority schema v{observed} refused; this algebra supports v{supported}"
            ),
            Self::CatalogSchemaVersionRefused {
                observed,
                supported,
            } => write!(
                f,
                "authority-catalog serialization v{observed} refused; this reader supports v{supported}"
            ),
            Self::MigrationUnavailable { observed, target } => write!(
                f,
                "authority schema v{observed} has no semantics-preserving migration to v{target}; regenerate from source evidence"
            ),
            Self::IdentityMismatch { what } => write!(f, "identity mismatch: {what}"),
            Self::IncompatibleAxes { what } => write!(f, "incompatible authority axes: {what}"),
            Self::IllegalTransition { from, to } => {
                write!(f, "illegal evidence transition {from} -> {to}")
            }
            Self::TerminalState { state } => write!(f, "authority state {state} is terminal"),
            Self::CompositionConflict { axis } => {
                write!(f, "authority composition conflict on {axis}")
            }
            Self::CapabilityMissing {
                capability,
                version,
            } => write!(
                f,
                "required capability `{capability}` v{version} is missing"
            ),
            Self::AssumptionNotAccepted { assumption } => {
                write!(
                    f,
                    "conditional-proof assumption `{assumption}` is not accepted"
                )
            }
            Self::NoClaimNotAccepted { boundary } => {
                write!(
                    f,
                    "no-claim boundary `{boundary}` is not accepted by policy"
                )
            }
            Self::CheckerRefused { verdict } => {
                write!(f, "checker did not accept: {}", verdict.code())
            }
            Self::RuntimeRequirementNotMet { requirement } => {
                write!(
                    f,
                    "runtime authority requirement `{requirement}` is not met"
                )
            }
            Self::AdjudicationNotRevocable => {
                write!(
                    f,
                    "only a genuine-counterexample adjudication can mint a revocation"
                )
            }
            Self::LogCapacityExceeded { needed, cap } => {
                write!(f, "authority log needs {needed} bytes, cap {cap}")
            }
        }
    }
}

impl std::error::Error for AuthorityError {}

fn preflight_text(what: &'static str, raw: &str) -> Result<(), AuthorityError> {
    if raw.len() > MAX_AUTHORITY_TEXT_BYTES {
        return Err(AuthorityError::TooLarge {
            what,
            observed: raw.len(),
            cap: MAX_AUTHORITY_TEXT_BYTES,
        });
    }
    if raw.split_whitespace().next().is_none() {
        return Err(AuthorityError::EmptyField { what });
    }
    Ok(())
}

fn canonical_text(what: &'static str, raw: &str) -> Result<String, AuthorityError> {
    preflight_text(what, raw)?;
    Ok(raw.split_whitespace().collect::<Vec<_>>().join(" "))
}

fn canonical_set(
    what: &'static str,
    values: &[&str],
    allow_empty: bool,
) -> Result<Vec<String>, AuthorityError> {
    if values.len() > MAX_AUTHORITY_SET_MEMBERS {
        return Err(AuthorityError::TooLarge {
            what,
            observed: values.len(),
            cap: MAX_AUTHORITY_SET_MEMBERS,
        });
    }
    if values.is_empty() && !allow_empty {
        return Err(AuthorityError::EmptyField { what });
    }
    let mut canonical = values
        .iter()
        .map(|value| canonical_text(what, value))
        .collect::<Result<Vec<_>, _>>()?;
    canonical.sort_unstable();
    canonical.dedup();
    Ok(canonical)
}

fn require_hash(what: &'static str, hash: ContentHash) -> Result<ContentHash, AuthorityError> {
    if hash.as_bytes().iter().all(|byte| *byte == 0) {
        Err(AuthorityError::MissingIdentity { what })
    } else {
        Ok(hash)
    }
}

#[derive(Default)]
struct CanonicalBytes(Vec<u8>);

impl CanonicalBytes {
    fn field(&mut self, tag: u8, bytes: &[u8]) {
        self.0.push(tag);
        self.0.extend_from_slice(
            &u64::try_from(bytes.len())
                .expect("bounded canonical field length fits u64")
                .to_le_bytes(),
        );
        self.0.extend_from_slice(bytes);
    }

    fn u8(&mut self, tag: u8, value: u8) {
        self.field(tag, &[value]);
    }

    fn u32(&mut self, tag: u8, value: u32) {
        self.field(tag, &value.to_le_bytes());
    }

    fn u64(&mut self, tag: u8, value: u64) {
        self.field(tag, &value.to_le_bytes());
    }

    fn hash(&mut self, tag: u8, value: ContentHash) {
        self.field(tag, value.as_bytes());
    }
}

macro_rules! typed_id {
    ($name:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(ContentHash);

        impl $name {
            #[must_use]
            pub fn as_hash(&self) -> &ContentHash {
                &self.0
            }
        }

        impl core::fmt::Display for $name {
            fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                write!(f, "{}", self.0)
            }
        }
    };
}

typed_id!(ClaimStatementId);
typed_id!(QuantifiedDomainId);
typed_id!(AssumptionSetId);
typed_id!(SemanticClaimId);
typed_id!(ClaimLaneBindingId);
typed_id!(ClaimInstanceId);
typed_id!(EvidenceId);
typed_id!(EvidenceStateId);
typed_id!(SatisfiableEvidenceId);
typed_id!(UnsatisfiableEvidenceId);
typed_id!(NonvacuousEvidenceId);
typed_id!(VacuousEvidenceId);
typed_id!(ReproductionFailedEvidenceId);
typed_id!(ReproducedEvidenceId);
typed_id!(ExactInstanceDecisionId);
typed_id!(InvalidationBindingId);
typed_id!(AuthorityStateId);
typed_id!(SupportEdgeId);
typed_id!(AttackEdgeId);
typed_id!(CounterexampleId);
typed_id!(AdjudicationId);
typed_id!(RevocationId);
typed_id!(VerifiedRevocationId);
typed_id!(CapabilityPolicyId);
typed_id!(CheckerDecisionId);
typed_id!(AuthorityHeadId);
typed_id!(RuntimeAdmissionId);
typed_id!(InferenceRuleId);
typed_id!(AuthorityMigrationId);

/// Canonical conjunction of statement clauses.  Clause order and cosmetic
/// whitespace do not change identity; adding, removing, or changing a clause
/// does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClaimStatement {
    clauses: Vec<String>,
    identity: ClaimStatementId,
}

impl ClaimStatement {
    pub fn new(clauses: &[&str]) -> Result<Self, AuthorityError> {
        let clauses = canonical_set("statement clause", clauses, false)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.u64(2, clauses.len() as u64);
        for clause in &clauses {
            bytes.field(3, clause.as_bytes());
        }
        let identity = ClaimStatementId(fs_blake3::hash_domain(
            CLAIM_STATEMENT_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self { clauses, identity })
    }

    #[must_use]
    pub fn clauses(&self) -> &[String] {
        &self.clauses
    }

    #[must_use]
    pub fn identity(&self) -> ClaimStatementId {
        self.identity
    }
}

/// Quantifier attached to one independent domain binding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Quantifier {
    ForAll,
    Exists,
}

impl Quantifier {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::ForAll => "forall",
            Self::Exists => "exists",
        }
    }

    fn tag(self) -> u8 {
        match self {
            Self::ForAll => 1,
            Self::Exists => 2,
        }
    }
}

/// One named binding inside an explicit quantifier block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DomainVariable {
    name: String,
    domain: String,
}

impl DomainVariable {
    pub fn new(name: &str, domain: &str) -> Result<Self, AuthorityError> {
        Ok(Self {
            name: canonical_text("domain variable", name)?,
            domain: canonical_text("variable domain", domain)?,
        })
    }

    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[must_use]
    pub fn domain(&self) -> &str {
        &self.domain
    }
}

/// One order-sensitive quantifier block. The caller must explicitly declare
/// whether its variables form a commutative product. Adjacent blocks are never
/// merged, even when their quantifiers agree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuantifierBlock {
    quantifier: Quantifier,
    commutative: bool,
    variables: Vec<DomainVariable>,
}

impl QuantifierBlock {
    pub fn ordered(
        quantifier: Quantifier,
        variables: Vec<DomainVariable>,
    ) -> Result<Self, AuthorityError> {
        Self::new(quantifier, variables, false)
    }

    pub fn commutative(
        quantifier: Quantifier,
        variables: Vec<DomainVariable>,
    ) -> Result<Self, AuthorityError> {
        Self::new(quantifier, variables, true)
    }

    fn new(
        quantifier: Quantifier,
        mut variables: Vec<DomainVariable>,
        commutative: bool,
    ) -> Result<Self, AuthorityError> {
        if variables.is_empty() {
            return Err(AuthorityError::EmptyField {
                what: "quantifier-block variables",
            });
        }
        if variables.len() > MAX_AUTHORITY_SET_MEMBERS {
            return Err(AuthorityError::TooLarge {
                what: "quantifier-block variables",
                observed: variables.len(),
                cap: MAX_AUTHORITY_SET_MEMBERS,
            });
        }
        if commutative {
            variables.sort_by(|left, right| left.name.cmp(&right.name));
        }
        let mut names = variables
            .iter()
            .map(|variable| variable.name.as_str())
            .collect::<Vec<_>>();
        names.sort_unstable();
        for pair in names.windows(2) {
            if pair[0] == pair[1] {
                return Err(AuthorityError::DuplicateMember {
                    what: "domain variable",
                    key: pair[0].to_string(),
                });
            }
        }
        Ok(Self {
            quantifier,
            commutative,
            variables,
        })
    }

    #[must_use]
    pub fn quantifier(&self) -> Quantifier {
        self.quantifier
    }

    #[must_use]
    pub fn is_commutative(&self) -> bool {
        self.commutative
    }

    #[must_use]
    pub fn variables(&self) -> &[DomainVariable] {
        &self.variables
    }
}

/// Canonical ordered quantifier blocks plus commutative predicate clauses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuantifiedDomain {
    blocks: Vec<QuantifierBlock>,
    predicates: Vec<String>,
    identity: QuantifiedDomainId,
}

impl QuantifiedDomain {
    pub fn new(blocks: Vec<QuantifierBlock>, predicates: &[&str]) -> Result<Self, AuthorityError> {
        if blocks.is_empty() {
            return Err(AuthorityError::EmptyField {
                what: "quantifier blocks",
            });
        }
        if blocks.len() > MAX_AUTHORITY_SET_MEMBERS {
            return Err(AuthorityError::TooLarge {
                what: "quantifier blocks",
                observed: blocks.len(),
                cap: MAX_AUTHORITY_SET_MEMBERS,
            });
        }
        let variable_count = blocks.iter().try_fold(0_usize, |count, block| {
            count.checked_add(block.variables.len())
        });
        let variable_count = variable_count.ok_or(AuthorityError::TooLarge {
            what: "quantified variables",
            observed: usize::MAX,
            cap: MAX_AUTHORITY_SET_MEMBERS,
        })?;
        if variable_count > MAX_AUTHORITY_SET_MEMBERS {
            return Err(AuthorityError::TooLarge {
                what: "quantified variables",
                observed: variable_count,
                cap: MAX_AUTHORITY_SET_MEMBERS,
            });
        }
        let mut names = blocks
            .iter()
            .flat_map(|block| block.variables.iter())
            .map(|variable| variable.name.as_str())
            .collect::<Vec<_>>();
        names.sort_unstable();
        for pair in names.windows(2) {
            if pair[0] == pair[1] {
                return Err(AuthorityError::DuplicateMember {
                    what: "domain variable",
                    key: pair[0].to_string(),
                });
            }
        }
        let predicates = canonical_set("domain predicate", predicates, true)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.u64(2, blocks.len() as u64);
        for block in &blocks {
            bytes.u8(3, block.quantifier.tag());
            bytes.u8(4, u8::from(block.commutative));
            bytes.u64(5, block.variables.len() as u64);
            for variable in &block.variables {
                bytes.field(6, variable.name.as_bytes());
                bytes.field(7, variable.domain.as_bytes());
            }
        }
        bytes.u64(8, predicates.len() as u64);
        for predicate in &predicates {
            bytes.field(9, predicate.as_bytes());
        }
        let identity = QuantifiedDomainId(fs_blake3::hash_domain(
            QUANTIFIED_DOMAIN_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            blocks,
            predicates,
            identity,
        })
    }

    #[must_use]
    pub fn blocks(&self) -> &[QuantifierBlock] {
        &self.blocks
    }

    #[must_use]
    pub fn predicates(&self) -> &[String] {
        &self.predicates
    }

    #[must_use]
    pub fn identity(&self) -> QuantifiedDomainId {
        self.identity
    }
}

/// Canonical assumption conjunction.  Empty means explicitly unconditional,
/// not absent metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssumptionSet {
    assumptions: Vec<String>,
    identity: AssumptionSetId,
}

impl AssumptionSet {
    pub fn new(assumptions: &[&str]) -> Result<Self, AuthorityError> {
        let assumptions = canonical_set("assumption", assumptions, true)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.u64(2, assumptions.len() as u64);
        for assumption in &assumptions {
            bytes.field(3, assumption.as_bytes());
        }
        let identity = AssumptionSetId(fs_blake3::hash_domain(
            ASSUMPTION_SET_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            assumptions,
            identity,
        })
    }

    #[must_use]
    pub fn assumptions(&self) -> &[String] {
        &self.assumptions
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.assumptions.is_empty()
    }

    #[must_use]
    pub fn identity(&self) -> AssumptionSetId {
        self.identity
    }
}

/// Exact unit factor.  Equivalence is structural: factor order is ignored and
/// duplicate symbols are exponent-summed; no heuristic alias table is trusted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnitFactor {
    symbol: String,
    exponent: i8,
}

impl UnitFactor {
    pub fn new(symbol: &str, exponent: i8) -> Result<Self, AuthorityError> {
        if exponent == 0 {
            return Err(AuthorityError::InvalidValue {
                what: "unit exponent must be nonzero",
            });
        }
        Ok(Self {
            symbol: canonical_text("unit symbol", symbol)?,
            exponent,
        })
    }

    #[must_use]
    pub fn symbol(&self) -> &str {
        &self.symbol
    }

    #[must_use]
    pub fn exponent(&self) -> i8 {
        self.exponent
    }
}

/// Exact rational unit scale and canonical factors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnitSystem {
    scale_numerator: u64,
    scale_denominator: u64,
    factors: Vec<UnitFactor>,
}

impl UnitSystem {
    pub fn new(
        scale_numerator: u64,
        scale_denominator: u64,
        factors: Vec<UnitFactor>,
    ) -> Result<Self, AuthorityError> {
        if scale_numerator == 0 || scale_denominator == 0 {
            return Err(AuthorityError::InvalidValue {
                what: "unit scale must be a positive rational",
            });
        }
        if factors.len() > MAX_AUTHORITY_SET_MEMBERS {
            return Err(AuthorityError::TooLarge {
                what: "unit factors",
                observed: factors.len(),
                cap: MAX_AUTHORITY_SET_MEMBERS,
            });
        }
        let divisor = gcd(scale_numerator, scale_denominator);
        let mut merged = BTreeMap::<String, i16>::new();
        for factor in factors {
            let exponent = merged.entry(factor.symbol).or_default();
            *exponent = exponent.checked_add(i16::from(factor.exponent)).ok_or(
                AuthorityError::InvalidValue {
                    what: "unit exponent sum overflow",
                },
            )?;
        }
        let mut factors = Vec::with_capacity(merged.len());
        for (symbol, exponent) in merged {
            if exponent == 0 {
                continue;
            }
            let exponent = i8::try_from(exponent).map_err(|_| AuthorityError::InvalidValue {
                what: "final unit exponent sum exceeds i8",
            })?;
            factors.push(UnitFactor { symbol, exponent });
        }
        Ok(Self {
            scale_numerator: scale_numerator / divisor,
            scale_denominator: scale_denominator / divisor,
            factors,
        })
    }

    #[must_use]
    pub fn dimensionless() -> Self {
        Self {
            scale_numerator: 1,
            scale_denominator: 1,
            factors: Vec::new(),
        }
    }

    #[must_use]
    pub fn scale(&self) -> (u64, u64) {
        (self.scale_numerator, self.scale_denominator)
    }

    #[must_use]
    pub fn factors(&self) -> &[UnitFactor] {
        &self.factors
    }

    fn encode(&self, bytes: &mut CanonicalBytes, base: u8) {
        bytes.u64(base, self.scale_numerator);
        bytes.u64(base + 1, self.scale_denominator);
        bytes.u64(base + 2, self.factors.len() as u64);
        for factor in &self.factors {
            bytes.field(base + 3, factor.symbol.as_bytes());
            bytes.u8(base + 4, factor.exponent as u8);
        }
    }
}

const fn gcd(mut left: u64, mut right: u64) -> u64 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

/// Explicit, dimensioned resource budget.  Zero is a meaningful explicit
/// value (for example a no-reviewer dry diagnostic), so it is not treated as
/// missing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AuthorityBudget {
    pub work_units: u64,
    pub memory_bytes: u64,
    pub wall_time_millis: u64,
    pub reviewer_slots: u32,
}

/// Version binding from one named component to its exact version string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VersionBinding {
    component: String,
    version: String,
}

impl VersionBinding {
    pub fn new(component: &str, version: &str) -> Result<Self, AuthorityError> {
        Ok(Self {
            component: canonical_text("version component", component)?,
            version: canonical_text("component version", version)?,
        })
    }

    #[must_use]
    pub fn component(&self) -> &str {
        &self.component
    }

    #[must_use]
    pub fn version(&self) -> &str {
        &self.version
    }
}

/// Exact capability name/version pair.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CapabilityBinding {
    name: String,
    version: u32,
}

impl CapabilityBinding {
    pub fn new(name: &str, version: u32) -> Result<Self, AuthorityError> {
        if version == 0 {
            return Err(AuthorityError::InvalidValue {
                what: "capability version must be nonzero",
            });
        }
        Ok(Self {
            name: canonical_text("capability name", name)?,
            version,
        })
    }

    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[must_use]
    pub fn version(&self) -> u32 {
        self.version
    }
}

/// The Five Explicits: units, seed, budgets, versions, and capabilities.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiveExplicits {
    units: UnitSystem,
    seed: u64,
    budget: AuthorityBudget,
    versions: Vec<VersionBinding>,
    capabilities: Vec<CapabilityBinding>,
}

impl FiveExplicits {
    pub fn new(
        units: UnitSystem,
        seed: u64,
        budget: AuthorityBudget,
        mut versions: Vec<VersionBinding>,
        mut capabilities: Vec<CapabilityBinding>,
    ) -> Result<Self, AuthorityError> {
        if versions.len() > MAX_AUTHORITY_SET_MEMBERS {
            return Err(AuthorityError::TooLarge {
                what: "version bindings",
                observed: versions.len(),
                cap: MAX_AUTHORITY_SET_MEMBERS,
            });
        }
        if capabilities.len() > MAX_AUTHORITY_SET_MEMBERS {
            return Err(AuthorityError::TooLarge {
                what: "capability bindings",
                observed: capabilities.len(),
                cap: MAX_AUTHORITY_SET_MEMBERS,
            });
        }
        versions.sort_by(|left, right| left.component.cmp(&right.component));
        for pair in versions.windows(2) {
            if pair[0].component == pair[1].component {
                return Err(AuthorityError::DuplicateMember {
                    what: "version binding",
                    key: pair[0].component.clone(),
                });
            }
        }
        capabilities.sort();
        for pair in capabilities.windows(2) {
            if pair[0].name == pair[1].name {
                return Err(AuthorityError::DuplicateMember {
                    what: "capability binding",
                    key: pair[0].name.clone(),
                });
            }
        }
        Ok(Self {
            units,
            seed,
            budget,
            versions,
            capabilities,
        })
    }

    #[must_use]
    pub fn units(&self) -> &UnitSystem {
        &self.units
    }

    #[must_use]
    pub fn seed(&self) -> u64 {
        self.seed
    }

    #[must_use]
    pub fn budget(&self) -> AuthorityBudget {
        self.budget
    }

    #[must_use]
    pub fn versions(&self) -> &[VersionBinding] {
        &self.versions
    }

    #[must_use]
    pub fn capabilities(&self) -> &[CapabilityBinding] {
        &self.capabilities
    }

    fn encode(&self, bytes: &mut CanonicalBytes) {
        self.units.encode(bytes, 20);
        bytes.u64(30, self.seed);
        bytes.u64(31, self.budget.work_units);
        bytes.u64(32, self.budget.memory_bytes);
        bytes.u64(33, self.budget.wall_time_millis);
        bytes.u32(34, self.budget.reviewer_slots);
        bytes.u64(35, self.versions.len() as u64);
        for version in &self.versions {
            bytes.field(36, version.component.as_bytes());
            bytes.field(37, version.version.as_bytes());
        }
        bytes.u64(38, self.capabilities.len() as u64);
        for capability in &self.capabilities {
            bytes.field(39, capability.name.as_bytes());
            bytes.u32(40, capability.version);
        }
    }
}

/// Explicit boundaries on what the claim does not establish.  The set is
/// canonical and travels into both semantic and exact-instance identities.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NoClaimBoundary {
    entries: Vec<String>,
}

impl NoClaimBoundary {
    pub fn new(entries: &[&str]) -> Result<Self, AuthorityError> {
        Ok(Self {
            entries: canonical_set("no-claim boundary", entries, true)?,
        })
    }

    #[must_use]
    pub fn none() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    #[must_use]
    pub fn entries(&self) -> &[String] {
        &self.entries
    }
}

/// Explicit binding between one structured claim surface and one validated
/// proof-lane charter.  The binder/artifact identify the external semantic
/// review that asserted correspondence; the algebra never infers it from a
/// bare lane hash.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClaimLaneBinding {
    statement: ClaimStatementId,
    domain: QuantifiedDomainId,
    assumptions: AssumptionSetId,
    lane: ProofLaneId,
    binding_artifact: ContentHash,
    binder: ContentHash,
    identity: ClaimLaneBindingId,
}

impl ClaimLaneBinding {
    pub fn new(
        statement: &ClaimStatement,
        domain: &QuantifiedDomain,
        assumptions: &AssumptionSet,
        lane: &LaneCharter,
        binding_artifact: ContentHash,
        binder: ContentHash,
    ) -> Result<Self, AuthorityError> {
        let binding_artifact = require_hash("claim/lane binding artifact", binding_artifact)?;
        let binder = require_hash("claim/lane binder", binder)?;
        let lane = lane.lane_id();
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *statement.identity.as_hash());
        bytes.hash(3, *domain.identity.as_hash());
        bytes.hash(4, *assumptions.identity.as_hash());
        bytes.hash(5, *lane.as_hash());
        bytes.hash(6, binding_artifact);
        bytes.hash(7, binder);
        let identity = ClaimLaneBindingId(fs_blake3::hash_domain(
            CLAIM_LANE_BINDING_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            statement: statement.identity,
            domain: domain.identity,
            assumptions: assumptions.identity,
            lane,
            binding_artifact,
            binder,
            identity,
        })
    }

    #[must_use]
    pub fn lane(&self) -> ProofLaneId {
        self.lane
    }

    #[must_use]
    pub fn binding_artifact(&self) -> ContentHash {
        self.binding_artifact
    }

    #[must_use]
    pub fn binder(&self) -> ContentHash {
        self.binder
    }

    #[must_use]
    pub fn identity(&self) -> ClaimLaneBindingId {
        self.identity
    }
}

/// Exact semantic claim plus its execution/provenance context.  The semantic
/// identity binds statement, product domain, assumptions, units, and no-claim
/// boundaries.  The instance identity additionally binds the claim/lane
/// binding (including binder and artifact), seed, budget, versions,
/// capabilities, and the algebra version.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClaimInstance {
    statement: ClaimStatement,
    domain: QuantifiedDomain,
    assumptions: AssumptionSet,
    lane_binding: ClaimLaneBinding,
    proof_lane: ProofLaneId,
    explicits: FiveExplicits,
    no_claim: NoClaimBoundary,
    semantic_identity: SemanticClaimId,
    identity: ClaimInstanceId,
}

impl ClaimInstance {
    pub fn new(
        statement: ClaimStatement,
        domain: QuantifiedDomain,
        assumptions: AssumptionSet,
        lane_binding: ClaimLaneBinding,
        explicits: FiveExplicits,
        no_claim: NoClaimBoundary,
    ) -> Result<Self, AuthorityError> {
        if lane_binding.statement != statement.identity
            || lane_binding.domain != domain.identity
            || lane_binding.assumptions != assumptions.identity
        {
            return Err(AuthorityError::IdentityMismatch {
                what: "claim/lane binding does not bind this statement/domain/assumption triple",
            });
        }
        let proof_lane = lane_binding.lane;
        let mut semantic = CanonicalBytes::default();
        semantic.u32(1, AUTHORITY_ALGEBRA_VERSION);
        semantic.hash(2, *statement.identity().as_hash());
        semantic.hash(3, *domain.identity().as_hash());
        semantic.hash(4, *assumptions.identity().as_hash());
        explicits.units.encode(&mut semantic, 5);
        semantic.u64(15, no_claim.entries.len() as u64);
        for boundary in &no_claim.entries {
            semantic.field(16, boundary.as_bytes());
        }
        let semantic_identity = SemanticClaimId(fs_blake3::hash_domain(
            SEMANTIC_CLAIM_IDENTITY_DOMAIN,
            &semantic.0,
        ));

        let mut exact = CanonicalBytes::default();
        exact.u32(1, AUTHORITY_ALGEBRA_VERSION);
        exact.hash(2, *semantic_identity.as_hash());
        exact.hash(3, *lane_binding.identity.as_hash());
        explicits.encode(&mut exact);
        let identity = ClaimInstanceId(fs_blake3::hash_domain(
            CLAIM_INSTANCE_IDENTITY_DOMAIN,
            &exact.0,
        ));
        Ok(Self {
            statement,
            domain,
            assumptions,
            lane_binding,
            proof_lane,
            explicits,
            no_claim,
            semantic_identity,
            identity,
        })
    }

    #[must_use]
    pub fn statement(&self) -> &ClaimStatement {
        &self.statement
    }

    #[must_use]
    pub fn domain(&self) -> &QuantifiedDomain {
        &self.domain
    }

    #[must_use]
    pub fn assumptions(&self) -> &AssumptionSet {
        &self.assumptions
    }

    #[must_use]
    pub fn proof_lane(&self) -> ProofLaneId {
        self.proof_lane
    }

    #[must_use]
    pub fn lane_binding(&self) -> ClaimLaneBinding {
        self.lane_binding
    }

    #[must_use]
    pub fn explicits(&self) -> &FiveExplicits {
        &self.explicits
    }

    #[must_use]
    pub fn no_claim(&self) -> &NoClaimBoundary {
        &self.no_claim
    }

    #[must_use]
    pub fn semantic_identity(&self) -> SemanticClaimId {
        self.semantic_identity
    }

    #[must_use]
    pub fn identity(&self) -> ClaimInstanceId {
        self.identity
    }
}

/// Closed evidence-kind vocabulary.  The wrappers below prevent a
/// satisfiability artifact from being substituted for a nonvacuity artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum EvidenceKind {
    Satisfiable,
    Unsatisfiable,
    Nonvacuous,
    Vacuous,
    ReproductionFailed,
    Reproduced,
    KernelProof,
    ScaleQualification,
    Support,
    Attack,
    Counterexample,
    Adjudication,
    Revocation,
}

impl EvidenceKind {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::Satisfiable => "satisfiable",
            Self::Unsatisfiable => "unsatisfiable",
            Self::Nonvacuous => "nonvacuous",
            Self::Vacuous => "vacuous",
            Self::ReproductionFailed => "reproduction-failed",
            Self::Reproduced => "reproduced",
            Self::KernelProof => "kernel-proof",
            Self::ScaleQualification => "scale-qualification",
            Self::Support => "support",
            Self::Attack => "attack",
            Self::Counterexample => "counterexample",
            Self::Adjudication => "adjudication",
            Self::Revocation => "revocation",
        }
    }

    fn tag(self) -> u8 {
        match self {
            Self::Satisfiable => 1,
            Self::Unsatisfiable => 2,
            Self::Nonvacuous => 3,
            Self::Vacuous => 4,
            Self::ReproductionFailed => 5,
            Self::Reproduced => 6,
            Self::KernelProof => 7,
            Self::ScaleQualification => 8,
            Self::Support => 9,
            Self::Attack => 10,
            Self::Counterexample => 11,
            Self::Adjudication => 12,
            Self::Revocation => 13,
        }
    }
}

/// Content-addressed evidence reference bound to one exact claim and one
/// checker identity.  This is evidence data, not an authority capability.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceRef {
    kind: EvidenceKind,
    claim: ClaimInstanceId,
    artifact: ContentHash,
    checker: ContentHash,
    schema_version: u32,
    identity: EvidenceId,
}

impl EvidenceRef {
    pub fn new(
        kind: EvidenceKind,
        claim: ClaimInstanceId,
        artifact: ContentHash,
        checker: ContentHash,
        schema_version: u32,
    ) -> Result<Self, AuthorityError> {
        if schema_version != AUTHORITY_ALGEBRA_VERSION {
            return Err(AuthorityError::SchemaVersionRefused {
                observed: schema_version,
                supported: AUTHORITY_ALGEBRA_VERSION,
            });
        }
        let artifact = require_hash("evidence artifact", artifact)?;
        let checker = require_hash("evidence checker", checker)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, schema_version);
        bytes.u8(2, kind.tag());
        bytes.hash(3, *claim.as_hash());
        bytes.hash(4, artifact);
        bytes.hash(5, checker);
        let identity = EvidenceId(fs_blake3::hash_domain(
            EVIDENCE_REF_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            kind,
            claim,
            artifact,
            checker,
            schema_version,
            identity,
        })
    }

    #[must_use]
    pub fn kind(&self) -> EvidenceKind {
        self.kind
    }

    #[must_use]
    pub fn claim(&self) -> ClaimInstanceId {
        self.claim
    }

    #[must_use]
    pub fn artifact(&self) -> ContentHash {
        self.artifact
    }

    #[must_use]
    pub fn checker(&self) -> ContentHash {
        self.checker
    }

    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    pub fn identity(&self) -> EvidenceId {
        self.identity
    }
}

macro_rules! conclusion_evidence_wrapper {
    ($name:ident, $id:ident, $kind:ident, $domain:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub struct $name {
            evidence: EvidenceRef,
            identity: $id,
        }

        impl $name {
            pub fn new(
                claim: ClaimInstanceId,
                artifact: ContentHash,
                checker: ContentHash,
            ) -> Result<Self, AuthorityError> {
                let evidence = EvidenceRef::new(
                    EvidenceKind::$kind,
                    claim,
                    artifact,
                    checker,
                    AUTHORITY_ALGEBRA_VERSION,
                )?;
                let mut bytes = CanonicalBytes::default();
                bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
                bytes.hash(2, *evidence.identity.as_hash());
                Ok(Self {
                    evidence,
                    identity: $id(fs_blake3::hash_domain($domain, &bytes.0)),
                })
            }

            #[must_use]
            pub fn evidence(&self) -> EvidenceRef {
                self.evidence
            }

            #[must_use]
            pub fn identity(&self) -> $id {
                self.identity
            }
        }
    };
}

conclusion_evidence_wrapper!(
    SatisfiableEvidence,
    SatisfiableEvidenceId,
    Satisfiable,
    SATISFIABLE_EVIDENCE_IDENTITY_DOMAIN
);
conclusion_evidence_wrapper!(
    UnsatisfiableEvidence,
    UnsatisfiableEvidenceId,
    Unsatisfiable,
    UNSATISFIABLE_EVIDENCE_IDENTITY_DOMAIN
);
conclusion_evidence_wrapper!(
    ReproductionFailedEvidence,
    ReproductionFailedEvidenceId,
    ReproductionFailed,
    REPRODUCTION_FAILED_EVIDENCE_IDENTITY_DOMAIN
);
conclusion_evidence_wrapper!(
    ReproducedEvidence,
    ReproducedEvidenceId,
    Reproduced,
    REPRODUCED_EVIDENCE_IDENTITY_DOMAIN
);

/// Strength class for a nonvacuity witness.  Classes are deliberately
/// incomparable without an explicit versioned inference rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NonvacuityKind {
    Point,
    OpenFamily,
    PositiveMeasureFamily,
    ScaleFamily,
    Custom,
}

impl NonvacuityKind {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::Point => "point",
            Self::OpenFamily => "open-family",
            Self::PositiveMeasureFamily => "positive-measure-family",
            Self::ScaleFamily => "scale-family",
            Self::Custom => "custom",
        }
    }

    fn tag(self) -> u8 {
        match self {
            Self::Point => 1,
            Self::OpenFamily => 2,
            Self::PositiveMeasureFamily => 3,
            Self::ScaleFamily => 4,
            Self::Custom => 5,
        }
    }
}

/// Exact topology/measure/scale/custom context plus fibre identity for one
/// nonvacuity strength class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NonvacuityStrength {
    kind: NonvacuityKind,
    context: ContentHash,
    fibre: ContentHash,
}

impl NonvacuityStrength {
    fn new(
        kind: NonvacuityKind,
        context: ContentHash,
        fibre: ContentHash,
    ) -> Result<Self, AuthorityError> {
        Ok(Self {
            kind,
            context: require_hash("nonvacuity strength context", context)?,
            fibre: require_hash("nonvacuity fibre", fibre)?,
        })
    }

    pub fn point(point: ContentHash, fibre: ContentHash) -> Result<Self, AuthorityError> {
        Self::new(NonvacuityKind::Point, point, fibre)
    }

    pub fn open_family(topology: ContentHash, fibre: ContentHash) -> Result<Self, AuthorityError> {
        Self::new(NonvacuityKind::OpenFamily, topology, fibre)
    }

    pub fn positive_measure_family(
        measure: ContentHash,
        fibre: ContentHash,
    ) -> Result<Self, AuthorityError> {
        Self::new(NonvacuityKind::PositiveMeasureFamily, measure, fibre)
    }

    pub fn scale_family(
        scale_space: ContentHash,
        fibre: ContentHash,
    ) -> Result<Self, AuthorityError> {
        Self::new(NonvacuityKind::ScaleFamily, scale_space, fibre)
    }

    pub fn custom(policy: ContentHash, fibre: ContentHash) -> Result<Self, AuthorityError> {
        Self::new(NonvacuityKind::Custom, policy, fibre)
    }

    #[must_use]
    pub fn kind(&self) -> NonvacuityKind {
        self.kind
    }

    #[must_use]
    pub fn context(&self) -> ContentHash {
        self.context
    }

    #[must_use]
    pub fn fibre(&self) -> ContentHash {
        self.fibre
    }

    #[must_use]
    pub fn satisfies(&self, required: &Self) -> bool {
        self == required
    }

    fn encode(&self, bytes: &mut CanonicalBytes, tag: u8) {
        bytes.u8(tag, self.kind.tag());
        bytes.hash(tag + 1, self.context);
        bytes.hash(tag + 2, self.fibre);
    }
}

macro_rules! strength_bound_conclusion_wrapper {
    ($name:ident, $id:ident, $kind:ident, $domain:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub struct $name {
            evidence: EvidenceRef,
            strength: NonvacuityStrength,
            identity: $id,
        }

        impl $name {
            pub fn new(
                claim: ClaimInstanceId,
                artifact: ContentHash,
                checker: ContentHash,
                strength: NonvacuityStrength,
            ) -> Result<Self, AuthorityError> {
                let evidence = EvidenceRef::new(
                    EvidenceKind::$kind,
                    claim,
                    artifact,
                    checker,
                    AUTHORITY_ALGEBRA_VERSION,
                )?;
                let mut bytes = CanonicalBytes::default();
                bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
                bytes.hash(2, *evidence.identity.as_hash());
                strength.encode(&mut bytes, 3);
                Ok(Self {
                    evidence,
                    strength,
                    identity: $id(fs_blake3::hash_domain($domain, &bytes.0)),
                })
            }

            #[must_use]
            pub fn evidence(&self) -> EvidenceRef {
                self.evidence
            }

            #[must_use]
            pub fn strength(&self) -> NonvacuityStrength {
                self.strength
            }

            #[must_use]
            pub fn identity(&self) -> $id {
                self.identity
            }
        }
    };
}

strength_bound_conclusion_wrapper!(
    NonvacuousEvidence,
    NonvacuousEvidenceId,
    Nonvacuous,
    NONVACUOUS_EVIDENCE_IDENTITY_DOMAIN
);
strength_bound_conclusion_wrapper!(
    VacuousEvidence,
    VacuousEvidenceId,
    Vacuous,
    VACUOUS_EVIDENCE_IDENTITY_DOMAIN
);

/// Request/drain/finalize proof for cancellation.  Every component is
/// content-identified and cancellation is a terminal evidence state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CancellationReceipt {
    request: ContentHash,
    drain: ContentHash,
    finalize: ContentHash,
}

impl CancellationReceipt {
    pub fn new(
        request: ContentHash,
        drain: ContentHash,
        finalize: ContentHash,
    ) -> Result<Self, AuthorityError> {
        Ok(Self {
            request: require_hash("cancellation request", request)?,
            drain: require_hash("cancellation drain", drain)?,
            finalize: require_hash("cancellation finalize", finalize)?,
        })
    }

    #[must_use]
    pub fn request(&self) -> ContentHash {
        self.request
    }

    #[must_use]
    pub fn drain(&self) -> ContentHash {
        self.drain
    }

    #[must_use]
    pub fn finalize(&self) -> ContentHash {
        self.finalize
    }
}

/// Lifecycle of one evidence object.  This is deliberately separate from
/// truth, scale qualification, reproduction, and runtime admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceState {
    Proposed,
    Checked,
    Adjudicated,
    Failed(EvidenceId),
    Cancelled(CancellationReceipt),
}

impl EvidenceState {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::Proposed => "proposed",
            Self::Checked => "checked",
            Self::Adjudicated => "adjudicated",
            Self::Failed(_) => "failed",
            Self::Cancelled(_) => "cancelled",
        }
    }

    #[must_use]
    pub fn is_terminal(self) -> bool {
        matches!(
            self,
            Self::Adjudicated | Self::Failed(_) | Self::Cancelled(_)
        )
    }

    /// Monotone lifecycle transition.  Cancellation may terminate proposed or
    /// checked work only after request/drain/finalize evidence exists.
    pub fn transition(self, next: Self) -> Result<Self, AuthorityError> {
        if self.is_terminal() {
            return Err(AuthorityError::TerminalState { state: self.code() });
        }
        let valid = matches!(
            (self, next),
            (Self::Proposed, Self::Checked)
                | (Self::Checked, Self::Adjudicated)
                | (Self::Proposed | Self::Checked, Self::Failed(_))
                | (Self::Proposed | Self::Checked, Self::Cancelled(_))
        );
        if valid {
            Ok(next)
        } else {
            Err(AuthorityError::IllegalTransition {
                from: self.code(),
                to: next.code(),
            })
        }
    }
}

/// Versioned lifecycle record for one exact evidence reference.
#[derive(Debug, PartialEq, Eq)]
pub struct EvidenceLifecycle {
    evidence: EvidenceRef,
    state: EvidenceState,
    predecessor: Option<EvidenceStateId>,
    identity: EvidenceStateId,
}

impl EvidenceLifecycle {
    #[must_use]
    pub fn proposed(evidence: EvidenceRef) -> Self {
        Self::from_state(evidence, EvidenceState::Proposed, None)
    }

    pub fn transition(&mut self, next: EvidenceState) -> Result<(), AuthorityError> {
        let predecessor = self.identity;
        let state = self.state.transition(next)?;
        let successor = Self::from_state(self.evidence, state, Some(predecessor));
        *self = successor;
        Ok(())
    }

    #[must_use]
    pub fn evidence(&self) -> EvidenceRef {
        self.evidence
    }

    #[must_use]
    pub fn state(&self) -> EvidenceState {
        self.state
    }

    #[must_use]
    pub fn predecessor(&self) -> Option<EvidenceStateId> {
        self.predecessor
    }

    #[must_use]
    pub fn identity(&self) -> EvidenceStateId {
        self.identity
    }

    fn from_state(
        evidence: EvidenceRef,
        state: EvidenceState,
        predecessor: Option<EvidenceStateId>,
    ) -> Self {
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *evidence.identity.as_hash());
        match predecessor {
            None => bytes.u8(3, 0),
            Some(previous) => {
                bytes.u8(3, 1);
                bytes.hash(4, *previous.as_hash());
            }
        }
        match state {
            EvidenceState::Proposed => bytes.u8(5, 0),
            EvidenceState::Checked => bytes.u8(5, 1),
            EvidenceState::Adjudicated => bytes.u8(5, 2),
            EvidenceState::Failed(failure) => {
                bytes.u8(5, 3);
                bytes.hash(6, *failure.as_hash());
            }
            EvidenceState::Cancelled(receipt) => {
                bytes.u8(5, 4);
                bytes.hash(7, receipt.request);
                bytes.hash(8, receipt.drain);
                bytes.hash(9, receipt.finalize);
            }
        }
        Self {
            evidence,
            state,
            predecessor,
            identity: EvidenceStateId(fs_blake3::hash_domain(
                EVIDENCE_STATE_IDENTITY_DOMAIN,
                &bytes.0,
            )),
        }
    }
}

/// Descriptive scientific-truth classification. `Refuted` is incomparable
/// with positive proof states; their conservative meet is the shared
/// `Unknown` bottom rather than an order-dependent choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TruthState {
    Unknown,
    ConditionalProof,
    Proved,
    Refuted,
}

impl TruthState {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::ConditionalProof => "conditional-proof",
            Self::Proved => "proved",
            Self::Refuted => "refuted",
        }
    }

    fn tag(self) -> u8 {
        match self {
            Self::Unknown => 0,
            Self::ConditionalProof => 1,
            Self::Proved => 2,
            Self::Refuted => 3,
        }
    }

    /// Partial order for positive authority.  Unknown is bottom;
    /// conditional proof is below proof; refutation is a separate branch.
    #[must_use]
    pub fn leq(self, other: Self) -> bool {
        self == other
            || self == Self::Unknown
            || (self == Self::ConditionalProof && other == Self::Proved)
    }

    fn meet(self, other: Self) -> Result<Self, AuthorityError> {
        if self.leq(other) {
            Ok(self)
        } else if other.leq(self) {
            Ok(other)
        } else {
            Ok(Self::Unknown)
        }
    }
}

/// Whether the quantified statement has a model.  This is intentionally not
/// inferred from proof or nonvacuity states.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SatisfiabilityState {
    Unknown,
    Satisfiable(SatisfiableEvidence),
    Unsatisfiable(UnsatisfiableEvidence),
}

impl SatisfiabilityState {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::Satisfiable(_) => "satisfiable",
            Self::Unsatisfiable(_) => "unsatisfiable",
        }
    }

    fn encode(self, bytes: &mut CanonicalBytes, tag: u8) {
        match self {
            Self::Unknown => bytes.u8(tag, 0),
            Self::Satisfiable(evidence) => {
                bytes.u8(tag, 1);
                bytes.hash(tag + 1, *evidence.identity().as_hash());
            }
            Self::Unsatisfiable(evidence) => {
                bytes.u8(tag, 2);
                bytes.hash(tag + 1, *evidence.identity().as_hash());
            }
        }
    }

    fn meet(self, other: Self) -> Result<Self, AuthorityError> {
        if self == other {
            return Ok(self);
        }
        Ok(Self::Unknown)
    }

    fn supports(self, required: Self) -> bool {
        required == Self::Unknown || self == required
    }
}

/// Whether the stated domain is known to contain a nontrivial admitted
/// instance.  Kept distinct from satisfiability because a syntactic model and
/// a strength-matched nonvacuity witness are different evidence obligations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NonvacuityState {
    Unknown,
    Nonvacuous(NonvacuousEvidence),
    Vacuous(VacuousEvidence),
}

impl NonvacuityState {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::Nonvacuous(_) => "nonvacuous",
            Self::Vacuous(_) => "vacuous",
        }
    }

    fn encode(self, bytes: &mut CanonicalBytes, tag: u8) {
        match self {
            Self::Unknown => bytes.u8(tag, 0),
            Self::Nonvacuous(evidence) => {
                bytes.u8(tag, 1);
                bytes.hash(tag + 1, *evidence.identity().as_hash());
            }
            Self::Vacuous(evidence) => {
                bytes.u8(tag, 2);
                bytes.hash(tag + 1, *evidence.identity().as_hash());
            }
        }
    }

    fn meet(self, other: Self) -> Result<Self, AuthorityError> {
        if self == other {
            return Ok(self);
        }
        Ok(Self::Unknown)
    }

    fn supports(self, required: Self) -> bool {
        required == Self::Unknown || self == required
    }
}

/// Exact-instance decision polarity.  The verdict is bound into the decision
/// identity and revalidated against the state variant that carries it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExactInstanceVerdict {
    Refused,
    Admitted,
}

impl ExactInstanceVerdict {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::Refused => "refused",
            Self::Admitted => "admitted",
        }
    }

    fn tag(self) -> u8 {
        match self {
            Self::Refused => 1,
            Self::Admitted => 2,
        }
    }
}

/// Descriptive exact-instance decision candidate. It binds the exact claim,
/// policy, checker, verdict, artifact, and schema; it authenticates none of
/// them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExactInstanceDecisionCandidate {
    claim: ClaimInstanceId,
    policy: CapabilityPolicyId,
    checker: ContentHash,
    verdict: ExactInstanceVerdict,
    artifact: ContentHash,
    schema_version: u32,
    identity: ExactInstanceDecisionId,
}

impl ExactInstanceDecisionCandidate {
    pub fn new(
        claim: ClaimInstanceId,
        policy: CapabilityPolicyId,
        checker: ContentHash,
        verdict: ExactInstanceVerdict,
        artifact: ContentHash,
        schema_version: u32,
    ) -> Result<Self, AuthorityError> {
        if schema_version != AUTHORITY_ALGEBRA_VERSION {
            return Err(AuthorityError::SchemaVersionRefused {
                observed: schema_version,
                supported: AUTHORITY_ALGEBRA_VERSION,
            });
        }
        let checker = require_hash("exact-instance checker", checker)?;
        let artifact = require_hash("exact-instance decision artifact", artifact)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, schema_version);
        bytes.hash(2, *claim.as_hash());
        bytes.hash(3, *policy.as_hash());
        bytes.hash(4, checker);
        bytes.u8(5, verdict.tag());
        bytes.hash(6, artifact);
        let identity = ExactInstanceDecisionId(fs_blake3::hash_domain(
            EXACT_INSTANCE_DECISION_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            claim,
            policy,
            checker,
            verdict,
            artifact,
            schema_version,
            identity,
        })
    }

    #[must_use]
    pub fn claim(&self) -> ClaimInstanceId {
        self.claim
    }

    #[must_use]
    pub fn policy(&self) -> CapabilityPolicyId {
        self.policy
    }

    #[must_use]
    pub fn checker(&self) -> ContentHash {
        self.checker
    }

    #[must_use]
    pub fn verdict(&self) -> ExactInstanceVerdict {
        self.verdict
    }

    #[must_use]
    pub fn artifact(&self) -> ContentHash {
        self.artifact
    }

    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    pub fn identity(&self) -> ExactInstanceDecisionId {
        self.identity
    }
}

/// Exact-instance admission axis. Admission and refusal carry a fully bound
/// immutable decision candidate rather than an untyped receipt hash.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExactInstanceAdmission {
    NotEvaluated,
    Refused(ExactInstanceDecisionCandidate),
    Admitted(ExactInstanceDecisionCandidate),
}

impl ExactInstanceAdmission {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::NotEvaluated => "not-evaluated",
            Self::Refused(_) => "refused",
            Self::Admitted(_) => "admitted",
        }
    }

    fn encode(self, bytes: &mut CanonicalBytes, tag: u8) {
        match self {
            Self::NotEvaluated => bytes.u8(tag, 0),
            Self::Refused(decision) => {
                bytes.u8(tag, 1);
                bytes.hash(tag + 1, *decision.identity.as_hash());
            }
            Self::Admitted(decision) => {
                bytes.u8(tag, 2);
                bytes.hash(tag + 1, *decision.identity.as_hash());
            }
        }
    }

    fn meet(self, other: Self) -> Result<Self, AuthorityError> {
        if self == other {
            return Ok(self);
        }
        Ok(Self::NotEvaluated)
    }
}

/// Independent proof-kernel checking axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelState {
    NotChecked,
    KernelChecked(EvidenceRef),
}

impl KernelState {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::NotChecked => "not-checked",
            Self::KernelChecked(_) => "kernel-checked",
        }
    }

    fn meet(self, other: Self) -> Self {
        if self == other {
            self
        } else {
            Self::NotChecked
        }
    }
}

/// Independent empirical/scale qualification axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScaleState {
    NotQualified,
    ScaleQualified(EvidenceRef),
}

impl ScaleState {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::NotQualified => "not-qualified",
            Self::ScaleQualified(_) => "scale-qualified",
        }
    }

    fn meet(self, other: Self) -> Self {
        if self == other {
            self
        } else {
            Self::NotQualified
        }
    }
}

/// Independent reproduction axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReproductionState {
    NotAttempted,
    Failed(ReproductionFailedEvidence),
    Reproduced(ReproducedEvidence),
}

impl ReproductionState {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::NotAttempted => "not-attempted",
            Self::Failed(_) => "failed",
            Self::Reproduced(_) => "reproduced",
        }
    }

    fn meet(self, other: Self) -> Result<Self, AuthorityError> {
        if self == other {
            return Ok(self);
        }
        Ok(Self::NotAttempted)
    }
}

/// Exact live-head binding for an authenticated invalidation transition. The
/// fields and constructor are private so a public hash cannot forge an
/// invalidated authority state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvalidationBinding {
    target_claim: ClaimInstanceId,
    target_state: AuthorityStateId,
    target_head: AuthorityHeadId,
    target_generation: u64,
    tombstone: VerifiedRevocationId,
    identity: InvalidationBindingId,
}

impl InvalidationBinding {
    fn new(
        target_claim: ClaimInstanceId,
        target_state: AuthorityStateId,
        target_head: AuthorityHeadId,
        target_generation: u64,
        tombstone: VerifiedRevocationId,
    ) -> Self {
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *target_claim.as_hash());
        bytes.hash(3, *target_state.as_hash());
        bytes.hash(4, *target_head.as_hash());
        bytes.u64(5, target_generation);
        bytes.hash(6, *tombstone.as_hash());
        Self {
            target_claim,
            target_state,
            target_head,
            target_generation,
            tombstone,
            identity: InvalidationBindingId(fs_blake3::hash_domain(
                INVALIDATION_BINDING_IDENTITY_DOMAIN,
                &bytes.0,
            )),
        }
    }

    #[must_use]
    pub fn target_claim(&self) -> ClaimInstanceId {
        self.target_claim
    }

    #[must_use]
    pub fn target_state(&self) -> AuthorityStateId {
        self.target_state
    }

    #[must_use]
    pub fn target_head(&self) -> AuthorityHeadId {
        self.target_head
    }

    #[must_use]
    pub fn target_generation(&self) -> u64 {
        self.target_generation
    }

    #[must_use]
    pub fn tombstone(&self) -> VerifiedRevocationId {
        self.tombstone
    }

    #[must_use]
    pub fn identity(&self) -> InvalidationBindingId {
        self.identity
    }
}

/// Transitive invalidation axis. Its representation is private: callers can
/// inspect a binding obtained from an authority state but cannot construct an
/// invalidated value from it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvalidationState {
    binding: Option<InvalidationBinding>,
}

impl InvalidationState {
    #[allow(non_upper_case_globals)]
    pub const Clear: Self = Self { binding: None };

    fn invalidated(binding: InvalidationBinding) -> Self {
        Self {
            binding: Some(binding),
        }
    }

    #[must_use]
    pub fn code(self) -> &'static str {
        if self.binding.is_some() {
            "invalidated"
        } else {
            "clear"
        }
    }

    #[must_use]
    pub fn is_clear(self) -> bool {
        self.binding.is_none()
    }

    #[must_use]
    pub fn binding(self) -> Option<InvalidationBinding> {
        self.binding
    }

    fn meet(self, other: Self) -> Result<Self, AuthorityError> {
        match (self.binding, other.binding) {
            (None, None) => Ok(Self::Clear),
            (Some(binding), None) | (None, Some(binding)) => Ok(Self::invalidated(binding)),
            (Some(left), Some(right)) if left == right => Ok(Self::invalidated(left)),
            (Some(_), Some(_)) => Err(AuthorityError::CompositionConflict {
                axis: "invalidation",
            }),
        }
    }
}

/// Descriptive product state for all scientific-authority candidate axes.
/// Fields are private and every constructor rechecks exact evidence-kind and
/// claim bindings, but raw artifacts/checker identities are not authenticated
/// here and this value cannot mint an [`AuthorityGrant`] publicly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuthorityState {
    claim: ClaimInstance,
    truth: TruthState,
    satisfiability: SatisfiabilityState,
    nonvacuity: NonvacuityState,
    exact_admission: ExactInstanceAdmission,
    kernel: KernelState,
    scale: ScaleState,
    reproduction: ReproductionState,
    invalidation: InvalidationState,
    identity: AuthorityStateId,
}

impl AuthorityState {
    /// Construct a clear descriptive state. Passing a previously obtained
    /// invalidation binding is refused: only the authenticated live-head
    /// revocation transition may create an invalidated successor.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        claim: ClaimInstance,
        truth: TruthState,
        satisfiability: SatisfiabilityState,
        nonvacuity: NonvacuityState,
        exact_admission: ExactInstanceAdmission,
        kernel: KernelState,
        scale: ScaleState,
        reproduction: ReproductionState,
        invalidation: InvalidationState,
    ) -> Result<Self, AuthorityError> {
        if !invalidation.is_clear() {
            return Err(AuthorityError::IllegalTransition {
                from: "public-authority-construction",
                to: "invalidated",
            });
        }
        Self::from_axes(
            claim,
            truth,
            satisfiability,
            nonvacuity,
            exact_admission,
            kernel,
            scale,
            reproduction,
            InvalidationState::Clear,
        )
    }

    #[allow(clippy::too_many_lines)]
    #[allow(clippy::too_many_arguments)]
    fn from_axes(
        claim: ClaimInstance,
        truth: TruthState,
        satisfiability: SatisfiabilityState,
        nonvacuity: NonvacuityState,
        exact_admission: ExactInstanceAdmission,
        kernel: KernelState,
        scale: ScaleState,
        reproduction: ReproductionState,
        invalidation: InvalidationState,
    ) -> Result<Self, AuthorityError> {
        let claim_id = claim.identity();
        if truth == TruthState::ConditionalProof && claim.assumptions.is_empty() {
            return Err(AuthorityError::IncompatibleAxes {
                what: "conditional proof requires an explicit nonempty assumption set",
            });
        }
        if truth == TruthState::Proved && !claim.assumptions.is_empty() {
            return Err(AuthorityError::IncompatibleAxes {
                what: "an assumption-bearing proof is conditional, not unqualified Proved",
            });
        }
        if matches!(satisfiability, SatisfiabilityState::Unsatisfiable(_))
            && matches!(nonvacuity, NonvacuityState::Nonvacuous(_))
        {
            return Err(AuthorityError::IncompatibleAxes {
                what: "an unsatisfiable domain cannot be nonvacuous",
            });
        }
        if matches!(exact_admission, ExactInstanceAdmission::Admitted(_))
            && truth == TruthState::Refuted
        {
            return Err(AuthorityError::IncompatibleAxes {
                what: "a refuted exact instance cannot remain admitted",
            });
        }
        if matches!(exact_admission, ExactInstanceAdmission::Admitted(_))
            && matches!(satisfiability, SatisfiabilityState::Unsatisfiable(_))
        {
            return Err(AuthorityError::IncompatibleAxes {
                what: "an unsatisfiable exact instance cannot remain admitted",
            });
        }
        if matches!(exact_admission, ExactInstanceAdmission::Admitted(_))
            && matches!(nonvacuity, NonvacuityState::Vacuous(_))
        {
            return Err(AuthorityError::IncompatibleAxes {
                what: "a vacuous exact instance cannot remain admitted",
            });
        }
        if matches!(exact_admission, ExactInstanceAdmission::Admitted(_))
            && !invalidation.is_clear()
        {
            return Err(AuthorityError::IncompatibleAxes {
                what: "an invalidated exact instance cannot remain admitted",
            });
        }
        if let Some(binding) = invalidation.binding() {
            if binding.target_claim != claim_id {
                return Err(AuthorityError::IdentityMismatch {
                    what: "invalidation binding targets another claim",
                });
            }
        }
        validate_exact_instance_admission(claim_id, exact_admission)?;
        validate_satisfiability_binding(claim_id, satisfiability)?;
        validate_nonvacuity_binding(claim_id, nonvacuity)?;
        validate_evidence_axis(claim_id, kernel_evidence(kernel), EvidenceKind::KernelProof)?;
        validate_evidence_axis(
            claim_id,
            scale_evidence(scale),
            EvidenceKind::ScaleQualification,
        )?;
        validate_reproduction_binding(claim_id, reproduction)?;

        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *claim_id.as_hash());
        bytes.u8(3, truth.tag());
        satisfiability.encode(&mut bytes, 4);
        nonvacuity.encode(&mut bytes, 6);
        exact_admission.encode(&mut bytes, 8);
        encode_optional_evidence(&mut bytes, 10, kernel_evidence(kernel));
        encode_optional_evidence(&mut bytes, 12, scale_evidence(scale));
        encode_reproduction(&mut bytes, 14, reproduction);
        match invalidation.binding() {
            None => bytes.u8(17, 0),
            Some(binding) => {
                bytes.u8(17, 1);
                bytes.hash(18, *binding.identity.as_hash());
            }
        }
        let identity = AuthorityStateId(fs_blake3::hash_domain(
            AUTHORITY_STATE_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            claim,
            truth,
            satisfiability,
            nonvacuity,
            exact_admission,
            kernel,
            scale,
            reproduction,
            invalidation,
            identity,
        })
    }

    pub fn unknown(claim: ClaimInstance) -> Result<Self, AuthorityError> {
        Self::new(
            claim,
            TruthState::Unknown,
            SatisfiabilityState::Unknown,
            NonvacuityState::Unknown,
            ExactInstanceAdmission::NotEvaluated,
            KernelState::NotChecked,
            ScaleState::NotQualified,
            ReproductionState::NotAttempted,
            InvalidationState::Clear,
        )
    }

    #[must_use]
    pub fn claim(&self) -> &ClaimInstance {
        &self.claim
    }

    #[must_use]
    pub fn truth(&self) -> TruthState {
        self.truth
    }

    #[must_use]
    pub fn satisfiability(&self) -> SatisfiabilityState {
        self.satisfiability
    }

    #[must_use]
    pub fn nonvacuity(&self) -> NonvacuityState {
        self.nonvacuity
    }

    #[must_use]
    pub fn exact_admission(&self) -> ExactInstanceAdmission {
        self.exact_admission
    }

    #[must_use]
    pub fn kernel(&self) -> KernelState {
        self.kernel
    }

    #[must_use]
    pub fn scale(&self) -> ScaleState {
        self.scale
    }

    #[must_use]
    pub fn reproduction(&self) -> ReproductionState {
        self.reproduction
    }

    #[must_use]
    pub fn invalidation(&self) -> InvalidationState {
        self.invalidation
    }

    #[must_use]
    pub fn identity(&self) -> AuthorityStateId {
        self.identity
    }

    /// Compare only the clear descriptive scientific-evidence product. Exact
    /// admission decisions are intentionally excluded from this relation.
    #[must_use]
    pub fn scientific_evidence_refines(&self, required: &Self) -> bool {
        self.invalidation.is_clear()
            && required.invalidation.is_clear()
            && self.claim.identity == required.claim.identity
            && required.truth.leq(self.truth)
            && self.satisfiability.supports(required.satisfiability)
            && self.nonvacuity.supports(required.nonvacuity)
            && optional_evidence_dominates(
                kernel_evidence(self.kernel),
                kernel_evidence(required.kernel),
            )
            && optional_evidence_dominates(
                scale_evidence(self.scale),
                scale_evidence(required.scale),
            )
            && reproduction_dominates(self.reproduction, required.reproduction)
    }

    /// Compare only the exact-instance decision axis on clear states.
    #[must_use]
    pub fn exact_decision_refines(&self, required: &Self) -> bool {
        self.invalidation.is_clear()
            && required.invalidation.is_clear()
            && self.claim.identity == required.claim.identity
            && admission_dominates(self.exact_admission, required.exact_admission)
    }

    /// Runtime substitutability requires both states to be clear and requires
    /// refinement on both the scientific and exact-decision relations.
    #[must_use]
    pub fn is_safe_runtime_substitute_for(&self, required: &Self) -> bool {
        self.scientific_evidence_refines(required) && self.exact_decision_refines(required)
    }

    /// Deny-biased evidence intersection. Any invalidation survives the meet,
    /// and distinct invalidation histories conflict rather than being guessed
    /// equivalent. This operation is deliberately **not** a substitutability
    /// greatest lower bound; lattice laws apply only to the clear descriptive
    /// scientific-evidence product.
    pub fn deny_biased_meet(&self, other: &Self) -> Result<Self, AuthorityError> {
        if self.claim.identity != other.claim.identity {
            return Err(AuthorityError::IdentityMismatch {
                what: "authority states refer to different exact claim instances",
            });
        }
        Self::from_axes(
            self.claim.clone(),
            self.truth.meet(other.truth)?,
            self.satisfiability.meet(other.satisfiability)?,
            self.nonvacuity.meet(other.nonvacuity)?,
            self.exact_admission.meet(other.exact_admission)?,
            self.kernel.meet(other.kernel),
            self.scale.meet(other.scale),
            self.reproduction.meet(other.reproduction)?,
            self.invalidation.meet(other.invalidation)?,
        )
    }

    #[must_use]
    pub fn unknown_grant(&self) -> AuthorityGrant<UnknownAuthority> {
        AuthorityGrant::new(self.clone())
    }

    #[allow(dead_code)]
    fn conditional_grant(&self) -> Result<AuthorityGrant<ConditionalAuthority>, AuthorityError> {
        if self.truth != TruthState::ConditionalProof {
            return Err(AuthorityError::RuntimeRequirementNotMet {
                requirement: "conditional-proof",
            });
        }
        Ok(AuthorityGrant::new(self.clone()))
    }

    #[allow(dead_code)]
    fn proved_grant(&self) -> Result<AuthorityGrant<ProvedAuthority>, AuthorityError> {
        if self.truth != TruthState::Proved {
            return Err(AuthorityError::RuntimeRequirementNotMet {
                requirement: "proved",
            });
        }
        Ok(AuthorityGrant::new(self.clone()))
    }

    #[allow(dead_code)]
    fn refuted_grant(&self) -> Result<AuthorityGrant<RefutedAuthority>, AuthorityError> {
        if self.truth != TruthState::Refuted {
            return Err(AuthorityError::RuntimeRequirementNotMet {
                requirement: "refuted",
            });
        }
        Ok(AuthorityGrant::new(self.clone()))
    }
}

#[cfg(test)]
mod authority_boundary_tests {
    use super::*;

    fn hash(label: &str) -> ContentHash {
        fs_blake3::hash_domain(
            "frankensim.fs-govern.evidence-contract-unit-test.v2",
            label.as_bytes(),
        )
    }

    fn claim() -> ClaimInstance {
        claim_with_seed(7)
    }

    fn claim_with_seed(seed: u64) -> ClaimInstance {
        let statement = ClaimStatement::new(&["residual is bounded"]).expect("statement");
        let domain = QuantifiedDomain::new(
            vec![
                QuantifierBlock::commutative(
                    Quantifier::ForAll,
                    vec![DomainVariable::new("mesh", "admitted meshes").expect("domain variable")],
                )
                .expect("quantifier block"),
            ],
            &[],
        )
        .expect("domain");
        let assumptions = AssumptionSet::new(&[]).expect("assumptions");
        let lane = crate::lanes::LaneCharter::new(
            "residual is bounded",
            "admitted meshes",
            &[],
            "checked",
            "unknown",
            "counterexample mesh",
            "authority-boundary-unit",
        )
        .expect("lane");
        let binding = ClaimLaneBinding::new(
            &statement,
            &domain,
            &assumptions,
            &lane,
            hash("lane-binding"),
            hash("lane-binder"),
        )
        .expect("lane binding");
        let explicits = FiveExplicits::new(
            UnitSystem::dimensionless(),
            seed,
            AuthorityBudget {
                work_units: 1,
                memory_bytes: 1,
                wall_time_millis: 1,
                reviewer_slots: 1,
            },
            vec![],
            vec![],
        )
        .expect("explicits");
        ClaimInstance::new(
            statement,
            domain,
            assumptions,
            binding,
            explicits,
            NoClaimBoundary::new(&[]).expect("no claim"),
        )
        .expect("claim")
    }

    fn evidence(claim: &ClaimInstance, kind: EvidenceKind, label: &str) -> EvidenceRef {
        EvidenceRef::new(
            kind,
            claim.identity(),
            hash(&format!("{label}-artifact")),
            hash(&format!("{label}-checker")),
            AUTHORITY_ALGEBRA_VERSION,
        )
        .expect("evidence")
    }

    fn policy() -> CapabilityPolicy {
        CapabilityPolicy::new(
            TruthRequirement::ProvedOnly,
            true,
            Some(
                NonvacuityStrength::scale_family(hash("scales"), hash("fibre")).expect("strength"),
            ),
            true,
            true,
            true,
            vec![],
            &[],
            &[],
        )
        .expect("policy")
    }

    fn proved_state(claim: ClaimInstance, policy: &CapabilityPolicy) -> AuthorityState {
        let satisfiability =
            SatisfiableEvidence::new(claim.identity(), hash("sat-artifact"), hash("sat-checker"))
                .expect("satisfiability");
        let nonvacuity = NonvacuousEvidence::new(
            claim.identity(),
            hash("nonvacuity-artifact"),
            hash("nonvacuity-checker"),
            NonvacuityStrength::scale_family(hash("scales"), hash("fibre")).expect("strength"),
        )
        .expect("nonvacuity");
        let exact_decision = ExactInstanceDecisionCandidate::new(
            claim.identity(),
            policy.identity(),
            hash("authenticated-checker"),
            ExactInstanceVerdict::Admitted,
            hash("decision-receipt"),
            AUTHORITY_ALGEBRA_VERSION,
        )
        .expect("exact decision");
        let reproduction = ReproducedEvidence::new(
            claim.identity(),
            hash("reproduction-artifact"),
            hash("reproduction-checker"),
        )
        .expect("reproduction");
        AuthorityState::new(
            claim.clone(),
            TruthState::Proved,
            SatisfiabilityState::Satisfiable(satisfiability),
            NonvacuityState::Nonvacuous(nonvacuity),
            ExactInstanceAdmission::Admitted(exact_decision),
            KernelState::KernelChecked(evidence(&claim, EvidenceKind::KernelProof, "kernel")),
            ScaleState::ScaleQualified(evidence(&claim, EvidenceKind::ScaleQualification, "scale")),
            ReproductionState::Reproduced(reproduction),
            InvalidationState::Clear,
        )
        .expect("proved state")
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn opaque_admission_requires_internal_authentication_and_current_head() {
        let policy = policy();
        let state = proved_state(claim(), &policy);
        let candidate = CheckerDecisionCandidate::new(
            state.claim().identity(),
            state.identity(),
            policy.identity(),
            hash("authenticated-checker"),
            CheckerVerdict::Accept,
            hash("decision-receipt"),
            None,
        )
        .expect("candidate");
        let decision = CheckerDecision::from_authenticated_candidate(candidate);
        let grant = state.proved_grant().expect("internal grant");
        let mut head = AuthorityHead::initial(&state);
        let other_state = proved_state(claim_with_seed(8), &policy);
        let other_head = AuthorityHead::initial(&other_state);
        assert!(matches!(
            RuntimeAdmission::admit(&grant, &policy, &decision, &other_head),
            Err(AuthorityError::IdentityMismatch { .. })
        ));
        let admission =
            RuntimeAdmission::admit(&grant, &policy, &decision, &head).expect("admission");
        admission.validate_current(&head).expect("current head");
        assert_eq!(admission.authority_head(), head.identity());
        assert_eq!(admission.head_generation(), 0);

        head.advance(&state, &state)
            .expect("same-state head advancement");
        assert_eq!(head.generation(), 1);
        assert_eq!(head.predecessor(), Some(admission.authority_head()));
        assert_eq!(
            admission
                .validate_current(&head)
                .expect_err("stale admission must fail closed"),
            AuthorityError::RuntimeRequirementNotMet {
                requirement: "current-authority-head"
            }
        );
        let refreshed = RuntimeAdmission::admit(&grant, &policy, &decision, &head)
            .expect("refreshed admission");
        refreshed.validate_current(&head).expect("refreshed head");

        let counterexample = CounterexampleCandidate::new(
            state.claim(),
            evidence(&state.claim, EvidenceKind::Counterexample, "counterexample"),
        )
        .expect("counterexample");
        let adjudication = CounterexampleAdjudication::new(
            &counterexample,
            CounterexampleVerdict::GenuineCounterexample,
            evidence(&state.claim, EvidenceKind::Adjudication, "adjudication"),
        )
        .expect("adjudication");
        let tombstone = RevocationTombstone::new(
            &state,
            &adjudication,
            "authenticated revocation fixture",
            evidence(&state.claim, EvidenceKind::Revocation, "revocation"),
        )
        .expect("tombstone");
        let verified = VerifiedRevocationTombstone::from_authenticated_candidate(&tombstone, &head)
            .expect("verified tombstone");
        let invalidated = head
            .revoke(&state, &verified)
            .expect("atomic revocation head advancement");
        assert_eq!(head.generation(), 2);
        assert_eq!(
            refreshed
                .validate_current(&head)
                .expect_err("pre-revocation admission must fail closed"),
            AuthorityError::RuntimeRequirementNotMet {
                requirement: "current-authority-head"
            }
        );
        let invalidated_grant = invalidated
            .proved_grant()
            .expect("historical truth classification remains proved");
        let invalidated_candidate = CheckerDecisionCandidate::new(
            invalidated.claim().identity(),
            invalidated.identity(),
            policy.identity(),
            hash("authenticated-checker"),
            CheckerVerdict::Accept,
            hash("invalidated-decision-receipt"),
            None,
        )
        .expect("invalidated candidate");
        let invalidated_decision =
            CheckerDecision::from_authenticated_candidate(invalidated_candidate);
        assert_eq!(
            RuntimeAdmission::admit(&invalidated_grant, &policy, &invalidated_decision, &head,)
                .expect_err("invalidated head cannot admit"),
            AuthorityError::RuntimeRequirementNotMet {
                requirement: "not-invalidated"
            }
        );

        let mut saturated = AuthorityHead::from_parts(
            state.claim().identity(),
            state.identity(),
            state.invalidation(),
            u64::MAX,
            None,
        );
        let saturated_id = saturated.identity();
        assert_eq!(
            saturated
                .advance(&state, &state)
                .expect_err("generation overflow must refuse"),
            AuthorityError::InvalidValue {
                what: "authority-head generation overflow"
            }
        );
        assert_eq!(saturated.identity(), saturated_id);
        assert_eq!(saturated.generation(), u64::MAX);

        let revoked_id = head.identity();
        assert!(matches!(
            head.advance(&invalidated, &state),
            Err(AuthorityError::IllegalTransition {
                from: "invalidated-authority-head",
                to: "clear"
            })
        ));
        assert_eq!(head.identity(), revoked_id);
        assert_eq!(head.invalidation(), invalidated.invalidation());

        let replacement_tombstone = RevocationTombstone::new(
            &state,
            &adjudication,
            "different authenticated revocation fixture",
            evidence(
                &state.claim,
                EvidenceKind::Revocation,
                "replacement-revocation",
            ),
        )
        .expect("replacement tombstone");
        let replacement_head = AuthorityHead::initial(&state);
        let replacement_verified = VerifiedRevocationTombstone::from_authenticated_candidate(
            &replacement_tombstone,
            &replacement_head,
        )
        .expect("replacement verified tombstone");
        assert!(matches!(
            head.revoke(&state, &replacement_verified),
            Err(AuthorityError::IdentityMismatch { .. })
        ));
        assert_eq!(head.identity(), revoked_id);
        assert_eq!(head.invalidation(), invalidated.invalidation());

        let mut stale_head = AuthorityHead::initial(&state);
        let stale_verified = VerifiedRevocationTombstone::from_authenticated_candidate(
            &replacement_tombstone,
            &stale_head,
        )
        .expect("stale fixture verification");
        stale_head
            .advance(&state, &state)
            .expect("advance after verification");
        assert!(matches!(
            stale_head.revoke(&state, &stale_verified),
            Err(AuthorityError::IdentityMismatch { .. })
        ));

        let mut clean_head = AuthorityHead::initial(&state);
        let other_counterexample = CounterexampleCandidate::new(
            other_state.claim(),
            evidence(
                other_state.claim(),
                EvidenceKind::Counterexample,
                "other-counterexample",
            ),
        )
        .expect("other counterexample");
        let other_adjudication = CounterexampleAdjudication::new(
            &other_counterexample,
            CounterexampleVerdict::GenuineCounterexample,
            evidence(
                other_state.claim(),
                EvidenceKind::Adjudication,
                "other-adjudication",
            ),
        )
        .expect("other adjudication");
        let other_tombstone = RevocationTombstone::new(
            &other_state,
            &other_adjudication,
            "wrong exact target",
            evidence(
                other_state.claim(),
                EvidenceKind::Revocation,
                "other-revocation",
            ),
        )
        .expect("other tombstone");
        let other_verified = VerifiedRevocationTombstone::from_authenticated_candidate(
            &other_tombstone,
            &other_head,
        )
        .expect("other verified tombstone");
        assert!(matches!(
            clean_head.revoke(&state, &other_verified),
            Err(AuthorityError::IdentityMismatch { .. })
        ));
    }
}

/// Ambiguous rank used by the pre-Phase-0B v0 record.  It deliberately mixed
/// proof, admission, and reproduction concepts and therefore cannot retain
/// positive authority during migration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LegacyAuthorityRankV0 {
    Unknown,
    Supported,
    Proved,
    Refuted,
}

impl LegacyAuthorityRankV0 {
    fn tag(self) -> u8 {
        match self {
            Self::Unknown => 0,
            Self::Supported => 1,
            Self::Proved => 2,
            Self::Refuted => 3,
        }
    }
}

/// Retained v0 authority record.  These booleans are intentionally available
/// only as historical input; no current API consumes them as authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegacyAuthorityV0 {
    claim: ClaimInstance,
    rank: LegacyAuthorityRankV0,
    admitted: bool,
    reproduced: bool,
    source_record: ContentHash,
}

impl LegacyAuthorityV0 {
    pub fn new(
        claim: ClaimInstance,
        rank: LegacyAuthorityRankV0,
        admitted: bool,
        reproduced: bool,
        source_record: ContentHash,
    ) -> Result<Self, AuthorityError> {
        Ok(Self {
            claim,
            rank,
            admitted,
            reproduced,
            source_record: require_hash("legacy authority source record", source_record)?,
        })
    }
}

/// Monotone historical migration result.  V0 positive booleans are retained
/// in the migration identity but demoted to orthogonal Unknown/NotEvaluated
/// axes because they lack the exact evidence bindings required by v2.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuthorityMigration {
    source_schema: u32,
    source_record: ContentHash,
    state: AuthorityState,
    demotions: &'static [&'static str],
    identity: AuthorityMigrationId,
}

impl AuthorityMigration {
    #[must_use]
    pub fn source_schema(&self) -> u32 {
        self.source_schema
    }

    #[must_use]
    pub fn source_record(&self) -> ContentHash {
        self.source_record
    }

    #[must_use]
    pub fn state(&self) -> &AuthorityState {
        &self.state
    }

    #[must_use]
    pub fn demotions(&self) -> &'static [&'static str] {
        self.demotions
    }

    #[must_use]
    pub fn identity(&self) -> AuthorityMigrationId {
        self.identity
    }
}

const LEGACY_V0_DEMOTIONS: &[&str] = &[
    "rank->truth:unknown",
    "admitted->exact-instance:not-evaluated",
    "reproduced->reproduction:not-attempted",
    "kernel:unbound",
    "scale:unbound",
];

/// Migrate the sole supported historical baseline without widening it.  Any
/// unknown/future source schema refuses; adding another migration requires a
/// new explicit function and test matrix.
pub fn migrate_legacy_v0(
    source_schema: u32,
    legacy: LegacyAuthorityV0,
) -> Result<AuthorityMigration, AuthorityError> {
    if source_schema == RETIRED_AUTHORITY_SCHEMA_VERSION {
        return Err(AuthorityError::MigrationUnavailable {
            observed: source_schema,
            target: AUTHORITY_ALGEBRA_VERSION,
        });
    }
    if source_schema != LEGACY_AUTHORITY_SCHEMA_VERSION {
        return Err(AuthorityError::SchemaVersionRefused {
            observed: source_schema,
            supported: LEGACY_AUTHORITY_SCHEMA_VERSION,
        });
    }
    let mut bytes = CanonicalBytes::default();
    bytes.u32(1, source_schema);
    bytes.u32(2, AUTHORITY_ALGEBRA_VERSION);
    bytes.hash(3, *legacy.claim.identity.as_hash());
    bytes.u8(4, legacy.rank.tag());
    bytes.u8(5, u8::from(legacy.admitted));
    bytes.u8(6, u8::from(legacy.reproduced));
    bytes.hash(7, legacy.source_record);
    for demotion in LEGACY_V0_DEMOTIONS {
        bytes.field(8, demotion.as_bytes());
    }
    let state = AuthorityState::unknown(legacy.claim)?;
    bytes.hash(9, *state.identity.as_hash());
    let identity = AuthorityMigrationId(fs_blake3::hash_domain(
        AUTHORITY_MIGRATION_IDENTITY_DOMAIN,
        &bytes.0,
    ));
    Ok(AuthorityMigration {
        source_schema,
        source_record: legacy.source_record,
        state,
        demotions: LEGACY_V0_DEMOTIONS,
        identity,
    })
}

/// One code-owned row in the generated authority-object catalog.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AuthorityCatalogRow {
    pub object_kind: &'static str,
    /// Authority algebra under which this object is interpreted.
    pub algebra_version: u32,
    /// Version of the object's own content-identity preimage/domain contract.
    pub identity_schema_version: u32,
    pub identity_domain: &'static str,
    pub identity_sources: &'static str,
    pub binding: &'static str,
    pub no_claim: &'static str,
}

/// Serialization schema for the catalog artifact itself. Catalog format,
/// authority algebra, and per-object identity schema are separate dimensions.
///
/// V3 is an intentional pre-freeze source-format correction. The row-level
/// identity-schema field was added while the outer marker was accidentally
/// left at v2, so that label can denote incompatible historical row shapes.
/// V2 must never be silently interpreted as v3 even though the authority
/// algebra itself remains v2.
pub const AUTHORITY_CATALOG_SCHEMA_VERSION: u32 = 3;
pub const RETIRED_AUTHORITY_CATALOG_SCHEMA_VERSION: u32 = 2;
pub const AUTHORITY_IDENTITY_SCHEMA_VERSION: u32 = 2;
pub const PROOF_LANE_IDENTITY_SCHEMA_VERSION: u32 = 1;

/// Fail-closed catalog-format admission gate for future wire readers.
///
/// This validates only the outer catalog serialization. It deliberately does
/// not compare or rewrite the independent authority-algebra and per-row
/// identity-schema axes.
pub fn validate_authority_catalog_schema_version(observed: u32) -> Result<(), AuthorityError> {
    if observed == AUTHORITY_CATALOG_SCHEMA_VERSION {
        Ok(())
    } else {
        Err(AuthorityError::CatalogSchemaVersionRefused {
            observed,
            supported: AUTHORITY_CATALOG_SCHEMA_VERSION,
        })
    }
}

/// Closed v3 serialization of the authority-algebra-v2 catalog. It is generated
/// from this code table, never hand-built by a dashboard. Contract tests compare
/// every row and every descriptive column against the generated Markdown
/// embedded in `fs-govern/CONTRACT.md`.
pub const AUTHORITY_CATALOG_ROWS: &[AuthorityCatalogRow] = &[
    AuthorityCatalogRow {
        object_kind: "claim-statement",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: CLAIM_STATEMENT_IDENTITY_DOMAIN,
        identity_sources: "canonical conjunction clauses",
        binding: "clause-order invariant; clause mutation moves identity",
        no_claim: "does not prove statement truth",
    },
    AuthorityCatalogRow {
        object_kind: "quantified-domain",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: QUANTIFIED_DOMAIN_IDENTITY_DOMAIN,
        identity_sources: "ordered quantifier blocks, explicit block commutativity, domain predicates",
        binding: "block order is semantic; only declared-commutative intra-block variables sort",
        no_claim: "does not prove satisfiability or nonvacuity",
    },
    AuthorityCatalogRow {
        object_kind: "assumption-set",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: ASSUMPTION_SET_IDENTITY_DOMAIN,
        identity_sources: "canonical assumption conjunction",
        binding: "assumption-order invariant; semantic mutation moves identity",
        no_claim: "does not discharge assumptions",
    },
    AuthorityCatalogRow {
        object_kind: "semantic-claim",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: SEMANTIC_CLAIM_IDENTITY_DOMAIN,
        identity_sources: "statement, domain, assumptions, exact units, no-claim",
        binding: "semantic root excludes execution budget/seed/version/capability context",
        no_claim: "semantic identity is not an exact execution instance",
    },
    AuthorityCatalogRow {
        object_kind: "claim-lane-binding",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: CLAIM_LANE_BINDING_IDENTITY_DOMAIN,
        identity_sources: "statement/domain/assumption roots, validated lane, binder, artifact",
        binding: "claim instances reject a binding minted for another structured claim",
        no_claim: "binding data does not authenticate a dishonest binder",
    },
    AuthorityCatalogRow {
        object_kind: "claim-instance",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: CLAIM_INSTANCE_IDENTITY_DOMAIN,
        identity_sources: "semantic claim, claim-lane binding, Five Explicits",
        binding: "semantic and exact-instance roots are distinct",
        no_claim: "content identity is not admission",
    },
    AuthorityCatalogRow {
        object_kind: "proof-lane",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: PROOF_LANE_IDENTITY_SCHEMA_VERSION,
        identity_domain: PROOF_LANE_IDENTITY_DOMAIN,
        identity_sources: "validated LaneCharter",
        binding: "reuses non-forgeable lanes::ProofLaneId",
        no_claim: "lane identity is not proof",
    },
    AuthorityCatalogRow {
        object_kind: "evidence-ref",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: EVIDENCE_REF_IDENTITY_DOMAIN,
        identity_sources: "kind, exact claim, artifact, checker, schema",
        binding: "conclusion polarity is encoded before typed wrapper construction",
        no_claim: "reference is not authenticated authority",
    },
    AuthorityCatalogRow {
        object_kind: "satisfiable-evidence",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: SATISFIABLE_EVIDENCE_IDENTITY_DOMAIN,
        identity_sources: "satisfiable evidence reference",
        binding: "type and identity both bind positive satisfiability polarity",
        no_claim: "candidate evidence is not authenticated satisfiability authority",
    },
    AuthorityCatalogRow {
        object_kind: "unsatisfiable-evidence",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: UNSATISFIABLE_EVIDENCE_IDENTITY_DOMAIN,
        identity_sources: "unsatisfiable evidence reference",
        binding: "type and identity both bind negative satisfiability polarity",
        no_claim: "candidate evidence is not authenticated unsatisfiability authority",
    },
    AuthorityCatalogRow {
        object_kind: "nonvacuous-evidence",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: NONVACUOUS_EVIDENCE_IDENTITY_DOMAIN,
        identity_sources: "nonvacuous evidence reference, strength kind, context, fibre",
        binding: "positive polarity and exact strength are type- and identity-bound",
        no_claim: "one strength class cannot widen without an inference rule",
    },
    AuthorityCatalogRow {
        object_kind: "vacuous-evidence",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: VACUOUS_EVIDENCE_IDENTITY_DOMAIN,
        identity_sources: "vacuous evidence reference, strength kind, context, fibre",
        binding: "negative polarity and exact strength are type- and identity-bound",
        no_claim: "one strength class cannot widen without an inference rule",
    },
    AuthorityCatalogRow {
        object_kind: "reproduction-failed-evidence",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: REPRODUCTION_FAILED_EVIDENCE_IDENTITY_DOMAIN,
        identity_sources: "reproduction-failed evidence reference",
        binding: "type and identity both bind failed-reproduction polarity",
        no_claim: "failure evidence is not authenticated adjudication",
    },
    AuthorityCatalogRow {
        object_kind: "reproduced-evidence",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: REPRODUCED_EVIDENCE_IDENTITY_DOMAIN,
        identity_sources: "reproduced evidence reference",
        binding: "type and identity both bind successful-reproduction polarity",
        no_claim: "reproduction is distinct from proof and admission",
    },
    AuthorityCatalogRow {
        object_kind: "evidence-state",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: EVIDENCE_STATE_IDENTITY_DOMAIN,
        identity_sources: "exact evidence reference, predecessor, lifecycle/cancellation fields",
        binding: "exclusive transitions replace the token; terminal states cannot revive",
        no_claim: "lifecycle completion does not establish statement truth",
    },
    AuthorityCatalogRow {
        object_kind: "authority-state",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: AUTHORITY_STATE_IDENTITY_DOMAIN,
        identity_sources: "exact claim and all orthogonal authority axes",
        binding: "split scientific/exact refinement, clear runtime substitution, deny-biased meet",
        no_claim: "descriptive classifications are not authenticated authority",
    },
    AuthorityCatalogRow {
        object_kind: "exact-instance-decision",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: EXACT_INSTANCE_DECISION_IDENTITY_DOMAIN,
        identity_sources: "claim, policy, checker, verdict, artifact, schema",
        binding: "runtime validation requires admitted polarity and exact field agreement",
        no_claim: "decision candidate is not an authenticated admission receipt",
    },
    AuthorityCatalogRow {
        object_kind: "invalidation-binding",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: INVALIDATION_BINDING_IDENTITY_DOMAIN,
        identity_sources: "target claim/state/head/generation and verified tombstone",
        binding: "private constructor binds the exact live predecessor atomically",
        no_claim: "binding is not cryptographic receipt authentication",
    },
    AuthorityCatalogRow {
        object_kind: "inference-rule",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: INFERENCE_RULE_IDENTITY_DOMAIN,
        identity_sources: "name, version, definition artifact",
        binding: "default rule set is empty",
        no_claim: "registered rule is not assumed sound",
    },
    AuthorityCatalogRow {
        object_kind: "support-edge",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: SUPPORT_EDGE_IDENTITY_DOMAIN,
        identity_sources: "source state, target claim/lane, rule, evidence",
        binding: "exact endpoint and rule identities",
        no_claim: "candidate edge is neither authenticated nor proof of graph acyclicity",
    },
    AuthorityCatalogRow {
        object_kind: "attack-edge",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: ATTACK_EDGE_IDENTITY_DOMAIN,
        identity_sources: "candidate, target claim/lane, evidence",
        binding: "candidate target/domain must match",
        no_claim: "attack is not adjudication",
    },
    AuthorityCatalogRow {
        object_kind: "counterexample-candidate",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: COUNTEREXAMPLE_IDENTITY_DOMAIN,
        identity_sources: "target claim/domain and counterexample evidence",
        binding: "exact candidate-to-domain identity",
        no_claim: "candidate is not a refutation",
    },
    AuthorityCatalogRow {
        object_kind: "counterexample-adjudication",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: ADJUDICATION_IDENTITY_DOMAIN,
        identity_sources: "candidate, target, verdict, adjudication evidence",
        binding: "only genuine verdict can derive a tombstone candidate",
        no_claim: "candidate adjudication cannot advance an authoritative head",
    },
    AuthorityCatalogRow {
        object_kind: "revocation-tombstone",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: REVOCATION_IDENTITY_DOMAIN,
        identity_sources: "target state, genuine adjudication, reason, evidence",
        binding: "permanent exact-state invalidation",
        no_claim: "descriptive tombstone is not an authenticated revocation receipt",
    },
    AuthorityCatalogRow {
        object_kind: "verified-revocation-tombstone",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: VERIFIED_REVOCATION_IDENTITY_DOMAIN,
        identity_sources: "authenticated tombstone candidate, target head and generation",
        binding: "opaque wrapper has no public minting constructor",
        no_claim: "wrapper alone makes no cryptographic-authority claim",
    },
    AuthorityCatalogRow {
        object_kind: "capability-policy",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: CAPABILITY_POLICY_IDENTITY_DOMAIN,
        identity_sources: "axis/strength requirements, capabilities, accepted assumptions/no-claims",
        binding: "every guard changes policy identity",
        no_claim: "policy data is neither capability possession nor checker authority",
    },
    AuthorityCatalogRow {
        object_kind: "checker-decision",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: CHECKER_DECISION_IDENTITY_DOMAIN,
        identity_sources: "claim, authority, policy, checker, verdict, artifact, cancellation",
        binding: "public candidate is exact data; opaque decision has no public Phase 0B-A mint",
        no_claim: "candidate verdict is neither authentication nor statement truth",
    },
    AuthorityCatalogRow {
        object_kind: "authority-head",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: AUTHORITY_HEAD_IDENTITY_DOMAIN,
        identity_sources: "claim, exact state, invalidation, generation, predecessor head",
        binding: "atomic advancement preserves permanent invalidation and replaces the head token",
        no_claim: "durable single-head authentication is Phase 0B-B scope",
    },
    AuthorityCatalogRow {
        object_kind: "runtime-admission",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: RUNTIME_ADMISSION_IDENTITY_DOMAIN,
        identity_sources: "claim, authority, policy, checker decision, current head identity/generation",
        binding: "positive typestate plus exact product-policy and live-head validation",
        no_claim: "does not widen claim scope or survive revocation",
    },
    AuthorityCatalogRow {
        object_kind: "authority-migration",
        algebra_version: AUTHORITY_ALGEBRA_VERSION,
        identity_schema_version: AUTHORITY_IDENTITY_SCHEMA_VERSION,
        identity_domain: AUTHORITY_MIGRATION_IDENTITY_DOMAIN,
        identity_sources: "legacy schema/record/rank/booleans and explicit demotions",
        binding: "v0 ambiguity demotes to v2 Unknown axes; persisted v1 is refused",
        no_claim: "migration never restores legacy positive authority",
    },
];

#[must_use]
pub fn authority_catalog() -> &'static [AuthorityCatalogRow] {
    AUTHORITY_CATALOG_ROWS
}

/// Deterministic Markdown rows generated from every catalog column. The
/// contract embeds this output verbatim so drift checks are bidirectional.
#[must_use]
pub fn authority_catalog_markdown_rows() -> String {
    use core::fmt::Write as _;

    let mut out = String::new();
    for row in AUTHORITY_CATALOG_ROWS {
        writeln!(
            out,
            "| `{}` | `{}` | `{}` | `{}` | {} | {} | {} |",
            row.object_kind,
            row.algebra_version,
            row.identity_schema_version,
            row.identity_domain,
            row.identity_sources,
            row.binding,
            row.no_claim,
        )
        .expect("writing to a String is infallible");
    }
    out
}

/// Deterministic JSON generated from [`AUTHORITY_CATALOG_ROWS`].
#[must_use]
pub fn authority_catalog_json() -> String {
    use core::fmt::Write as _;

    let mut out = format!(
        "{{\"schema\":\"frankensim-authority-catalog-v{AUTHORITY_CATALOG_SCHEMA_VERSION}\",\"catalog_schema_version\":{AUTHORITY_CATALOG_SCHEMA_VERSION},\"algebra_version\":"
    );
    write!(out, "{AUTHORITY_ALGEBRA_VERSION},\"rows\":[")
        .expect("writing to a String is infallible");
    for (index, row) in AUTHORITY_CATALOG_ROWS.iter().enumerate() {
        if index != 0 {
            out.push(',');
        }
        write!(
            out,
            "{{\"object_kind\":\"{}\",\"algebra_version\":{},\"identity_schema_version\":{},\"identity_domain\":\"{}\",\"identity_sources\":\"{}\",\"binding\":\"{}\",\"no_claim\":\"{}\"}}",
            json_escape(row.object_kind),
            row.algebra_version,
            row.identity_schema_version,
            json_escape(row.identity_domain),
            json_escape(row.identity_sources),
            json_escape(row.binding),
            json_escape(row.no_claim),
        )
        .expect("writing to a String is infallible");
    }
    out.push_str("]}");
    out
}

/// Bounded structured log for one authority state and its optional policy/
/// checker context.  It names the exact object/source identities, versions,
/// every authority axis, no-claim boundaries, and ranked remedy code.
pub fn authority_log_json(
    state: &AuthorityState,
    policy: Option<&CapabilityPolicy>,
    decision: Option<&CheckerDecisionCandidate>,
    error: Option<&AuthorityError>,
) -> Result<String, AuthorityError> {
    use core::fmt::Write as _;

    match (policy, decision) {
        (None, Some(_)) => {
            return Err(AuthorityError::IncompatibleAxes {
                what: "checker-candidate log context requires its policy",
            });
        }
        (Some(policy), Some(decision))
            if decision.claim != state.claim.identity
                || decision.authority != state.identity
                || decision.policy != policy.identity =>
        {
            return Err(AuthorityError::IdentityMismatch {
                what: "authority log context does not bind the exact state and policy",
            });
        }
        _ => {}
    }

    let policy_identity = policy
        .map(|value| value.identity.to_string())
        .unwrap_or_else(|| "none".to_string());
    let decision_identity = decision
        .map(|value| value.identity().to_string())
        .unwrap_or_else(|| "none".to_string());
    let checker_verdict = decision
        .map(|value| value.verdict().code())
        .unwrap_or("not-run");
    let checker_identity = decision
        .map(|value| value.checker().to_string())
        .unwrap_or_else(|| "none".to_string());
    let remedy = error.map_or("none", AuthorityError::remedy_code);
    let mut out = String::new();
    write!(
        out,
        "{{\"object_kind\":\"authority-state\",\"source_identity\":\"{}\",\"authority_identity\":\"{}\",\"schema_version\":{},\"policy_version\":{},\"policy_identity\":\"{}\",\"checker_decision\":\"{}\",\"checker_identity\":\"{}\",\"checker_verdict\":\"{}\",\"axes\":{{\"truth\":\"{}\",\"satisfiability\":\"{}\",\"nonvacuity\":\"{}\",\"exact_instance\":\"{}\",\"kernel\":\"{}\",\"scale\":\"{}\",\"reproduction\":\"{}\",\"invalidation\":\"{}\"}},\"no_claim\":[",
        state.claim.identity,
        state.identity,
        AUTHORITY_ALGEBRA_VERSION,
        AUTHORITY_POLICY_VERSION,
        policy_identity,
        decision_identity,
        checker_identity,
        checker_verdict,
        state.truth.code(),
        state.satisfiability.code(),
        state.nonvacuity.code(),
        state.exact_admission.code(),
        state.kernel.code(),
        state.scale.code(),
        state.reproduction.code(),
        state.invalidation.code(),
    )
    .expect("writing to a String is infallible");
    for (index, boundary) in state.claim.no_claim.entries.iter().enumerate() {
        if index != 0 {
            out.push(',');
        }
        write!(out, "\"{}\"", json_escape(boundary)).expect("writing to a String is infallible");
    }
    write!(out, "],\"remedy\":\"{}\"}}", json_escape(remedy))
        .expect("writing to a String is infallible");
    if out.len() > MAX_AUTHORITY_LOG_BYTES {
        return Err(AuthorityError::LogCapacityExceeded {
            needed: out.len(),
            cap: MAX_AUTHORITY_LOG_BYTES,
        });
    }
    Ok(out)
}

/// Minimum truth authority accepted by a runtime policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TruthRequirement {
    ConditionalOrProved,
    ProvedOnly,
}

impl TruthRequirement {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::ConditionalOrProved => "conditional-or-proved",
            Self::ProvedOnly => "proved-only",
        }
    }

    fn tag(self) -> u8 {
        match self {
            Self::ConditionalOrProved => 1,
            Self::ProvedOnly => 2,
        }
    }
}

/// Exact capability and authority-axis policy for runtime consumption.
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CapabilityPolicy {
    truth: TruthRequirement,
    require_satisfiable: bool,
    nonvacuity_requirement: Option<NonvacuityStrength>,
    require_kernel_checked: bool,
    require_scale_qualified: bool,
    require_reproduced: bool,
    required_capabilities: Vec<CapabilityBinding>,
    accepted_assumptions: Vec<String>,
    accepted_no_claims: Vec<String>,
    identity: CapabilityPolicyId,
}

impl CapabilityPolicy {
    #[allow(clippy::fn_params_excessive_bools)]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        truth: TruthRequirement,
        require_satisfiable: bool,
        nonvacuity_requirement: Option<NonvacuityStrength>,
        require_kernel_checked: bool,
        require_scale_qualified: bool,
        require_reproduced: bool,
        mut required_capabilities: Vec<CapabilityBinding>,
        accepted_assumptions: &[&str],
        accepted_no_claims: &[&str],
    ) -> Result<Self, AuthorityError> {
        if required_capabilities.len() > MAX_AUTHORITY_SET_MEMBERS {
            return Err(AuthorityError::TooLarge {
                what: "required capabilities",
                observed: required_capabilities.len(),
                cap: MAX_AUTHORITY_SET_MEMBERS,
            });
        }
        required_capabilities.sort();
        for pair in required_capabilities.windows(2) {
            if pair[0].name == pair[1].name {
                return Err(AuthorityError::DuplicateMember {
                    what: "required capability",
                    key: pair[0].name.clone(),
                });
            }
        }
        let accepted_assumptions = canonical_set(
            "accepted conditional assumption",
            accepted_assumptions,
            true,
        )?;
        let accepted_no_claims =
            canonical_set("accepted no-claim boundary", accepted_no_claims, true)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_POLICY_VERSION);
        bytes.u8(2, truth.tag());
        bytes.u8(3, u8::from(require_satisfiable));
        match nonvacuity_requirement {
            None => bytes.u8(4, 0),
            Some(required) => {
                bytes.u8(4, 1);
                required.encode(&mut bytes, 5);
            }
        }
        bytes.u8(8, u8::from(require_kernel_checked));
        bytes.u8(9, u8::from(require_scale_qualified));
        bytes.u8(10, u8::from(require_reproduced));
        bytes.u64(11, required_capabilities.len() as u64);
        for capability in &required_capabilities {
            bytes.field(12, capability.name.as_bytes());
            bytes.u32(13, capability.version);
        }
        bytes.u64(14, accepted_assumptions.len() as u64);
        for assumption in &accepted_assumptions {
            bytes.field(15, assumption.as_bytes());
        }
        bytes.u64(16, accepted_no_claims.len() as u64);
        for boundary in &accepted_no_claims {
            bytes.field(17, boundary.as_bytes());
        }
        let identity = CapabilityPolicyId(fs_blake3::hash_domain(
            CAPABILITY_POLICY_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            truth,
            require_satisfiable,
            nonvacuity_requirement,
            require_kernel_checked,
            require_scale_qualified,
            require_reproduced,
            required_capabilities,
            accepted_assumptions,
            accepted_no_claims,
            identity,
        })
    }

    #[must_use]
    pub fn truth(&self) -> TruthRequirement {
        self.truth
    }

    #[must_use]
    pub fn require_satisfiable(&self) -> bool {
        self.require_satisfiable
    }

    #[must_use]
    pub fn nonvacuity_requirement(&self) -> Option<NonvacuityStrength> {
        self.nonvacuity_requirement
    }

    #[must_use]
    pub fn require_kernel_checked(&self) -> bool {
        self.require_kernel_checked
    }

    #[must_use]
    pub fn require_scale_qualified(&self) -> bool {
        self.require_scale_qualified
    }

    #[must_use]
    pub fn require_reproduced(&self) -> bool {
        self.require_reproduced
    }

    #[must_use]
    pub fn required_capabilities(&self) -> &[CapabilityBinding] {
        &self.required_capabilities
    }

    #[must_use]
    pub fn accepted_assumptions(&self) -> &[String] {
        &self.accepted_assumptions
    }

    #[must_use]
    pub fn accepted_no_claims(&self) -> &[String] {
        &self.accepted_no_claims
    }

    #[must_use]
    pub fn identity(&self) -> CapabilityPolicyId {
        self.identity
    }
}

/// Checker result, distinct from truth, reproduction, and admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckerVerdict {
    Accept,
    Refuse,
    Indeterminate,
    Cancelled,
}

impl CheckerVerdict {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::Accept => "accept",
            Self::Refuse => "refuse",
            Self::Indeterminate => "indeterminate",
            Self::Cancelled => "cancelled",
        }
    }

    fn tag(self) -> u8 {
        match self {
            Self::Accept => 1,
            Self::Refuse => 2,
            Self::Indeterminate => 3,
            Self::Cancelled => 4,
        }
    }
}

/// Plain candidate decision. It becomes a [`CheckerDecision`] only through the
/// private receipt-authentication boundary implemented with Phase 0B-B.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CheckerDecisionCandidate {
    claim: ClaimInstanceId,
    authority: AuthorityStateId,
    policy: CapabilityPolicyId,
    checker: ContentHash,
    verdict: CheckerVerdict,
    decision_artifact: ContentHash,
    cancellation: Option<CancellationReceipt>,
    identity: CheckerDecisionId,
}

impl CheckerDecisionCandidate {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        claim: ClaimInstanceId,
        authority: AuthorityStateId,
        policy: CapabilityPolicyId,
        checker: ContentHash,
        verdict: CheckerVerdict,
        decision_artifact: ContentHash,
        cancellation: Option<CancellationReceipt>,
    ) -> Result<Self, AuthorityError> {
        let checker = require_hash("checker identity", checker)?;
        let decision_artifact = require_hash("checker decision artifact", decision_artifact)?;
        match (verdict, cancellation) {
            (CheckerVerdict::Cancelled, None) => {
                return Err(AuthorityError::IncompatibleAxes {
                    what: "cancelled checker decision needs request-drain-finalize evidence",
                });
            }
            (CheckerVerdict::Cancelled, Some(_)) | (_, None) => {}
            (_, Some(_)) => {
                return Err(AuthorityError::IncompatibleAxes {
                    what: "only a cancelled checker decision may carry cancellation evidence",
                });
            }
        }
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *claim.as_hash());
        bytes.hash(3, *authority.as_hash());
        bytes.hash(4, *policy.as_hash());
        bytes.hash(5, checker);
        bytes.u8(6, verdict.tag());
        bytes.hash(7, decision_artifact);
        match cancellation {
            None => bytes.u8(8, 0),
            Some(receipt) => {
                bytes.u8(8, 1);
                bytes.hash(9, receipt.request);
                bytes.hash(10, receipt.drain);
                bytes.hash(11, receipt.finalize);
            }
        }
        let identity = CheckerDecisionId(fs_blake3::hash_domain(
            CHECKER_DECISION_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            claim,
            authority,
            policy,
            checker,
            verdict,
            decision_artifact,
            cancellation,
            identity,
        })
    }

    #[must_use]
    pub fn identity(&self) -> CheckerDecisionId {
        self.identity
    }

    #[must_use]
    pub fn claim(&self) -> ClaimInstanceId {
        self.claim
    }

    #[must_use]
    pub fn authority(&self) -> AuthorityStateId {
        self.authority
    }

    #[must_use]
    pub fn policy(&self) -> CapabilityPolicyId {
        self.policy
    }

    #[must_use]
    pub fn checker(&self) -> ContentHash {
        self.checker
    }

    #[must_use]
    pub fn verdict(&self) -> CheckerVerdict {
        self.verdict
    }

    #[must_use]
    pub fn decision_artifact(&self) -> ContentHash {
        self.decision_artifact
    }

    #[must_use]
    pub fn cancellation(&self) -> Option<CancellationReceipt> {
        self.cancellation
    }
}

/// Pure, non-authoritative compatibility result for a descriptive state,
/// policy, and checker-decision candidate.  This value deliberately has no
/// conversion into [`AuthorityGrant`], [`CheckerDecision`], or
/// [`RuntimeAdmission`].  Phase 0B-B must authenticate receipts and the live
/// authority head before minting any of those opaque values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RuntimeAssessment {
    claim: ClaimInstanceId,
    authority: AuthorityStateId,
    policy: CapabilityPolicyId,
    decision_candidate: CheckerDecisionId,
}

impl RuntimeAssessment {
    #[must_use]
    pub fn claim(&self) -> ClaimInstanceId {
        self.claim
    }

    #[must_use]
    pub fn authority(&self) -> AuthorityStateId {
        self.authority
    }

    #[must_use]
    pub fn policy(&self) -> CapabilityPolicyId {
        self.policy
    }

    #[must_use]
    pub fn decision_candidate(&self) -> CheckerDecisionId {
        self.decision_candidate
    }
}

/// Evaluate the public data algebra without authenticating it.  `Ok` means
/// only "eligible candidate"; it is never runtime authority.
pub fn assess_runtime_candidate(
    state: &AuthorityState,
    policy: &CapabilityPolicy,
    candidate: &CheckerDecisionCandidate,
) -> Result<RuntimeAssessment, AuthorityError> {
    validate_runtime_candidate(state, policy, candidate)?;
    Ok(RuntimeAssessment {
        claim: state.claim.identity,
        authority: state.identity,
        policy: policy.identity,
        decision_candidate: candidate.identity,
    })
}

/// Capability-authenticated immutable checker decision.  Phase 0B-A exposes
/// this object for typed consumption but no public minting API: the durable
/// receipt verifier added by Phase 0B-B owns that trust boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CheckerDecision {
    candidate: CheckerDecisionCandidate,
}

impl CheckerDecision {
    #[allow(dead_code)]
    fn from_authenticated_candidate(candidate: CheckerDecisionCandidate) -> Self {
        Self { candidate }
    }

    #[must_use]
    pub fn claim(&self) -> ClaimInstanceId {
        self.candidate.claim
    }

    #[must_use]
    pub fn authority(&self) -> AuthorityStateId {
        self.candidate.authority
    }

    #[must_use]
    pub fn policy(&self) -> CapabilityPolicyId {
        self.candidate.policy
    }

    #[must_use]
    pub fn checker(&self) -> ContentHash {
        self.candidate.checker
    }

    #[must_use]
    pub fn verdict(&self) -> CheckerVerdict {
        self.candidate.verdict
    }

    #[must_use]
    pub fn decision_artifact(&self) -> ContentHash {
        self.candidate.decision_artifact
    }

    #[must_use]
    pub fn cancellation(&self) -> Option<CancellationReceipt> {
        self.candidate.cancellation
    }

    #[must_use]
    pub fn identity(&self) -> CheckerDecisionId {
        self.candidate.identity
    }
}

/// Opaque live-head token for one claim's authority lineage. It is neither
/// `Clone` nor `Copy`; validated advancement mutably replaces the token and
/// binds its prior identity as the predecessor. Phase 0B-B will persist and
/// authenticate the single current head. No public constructor exists in
/// Phase 0B-A.
#[derive(Debug, PartialEq, Eq)]
pub struct AuthorityHead {
    claim: ClaimInstanceId,
    state: AuthorityStateId,
    invalidation: InvalidationState,
    generation: u64,
    predecessor: Option<AuthorityHeadId>,
    identity: AuthorityHeadId,
}

impl AuthorityHead {
    #[allow(dead_code)]
    fn initial(state: &AuthorityState) -> Self {
        Self::from_parts(
            state.claim.identity,
            state.identity,
            state.invalidation,
            0,
            None,
        )
    }

    #[allow(dead_code)]
    fn advance(
        &mut self,
        current: &AuthorityState,
        successor_state: &AuthorityState,
    ) -> Result<(), AuthorityError> {
        if self.state != current.identity || self.claim != current.claim.identity {
            return Err(AuthorityError::IdentityMismatch {
                what: "authority-head predecessor is stale or belongs to another claim",
            });
        }
        if successor_state.claim.identity != self.claim {
            return Err(AuthorityError::IdentityMismatch {
                what: "authority-head successor belongs to another claim",
            });
        }
        if !self.invalidation.is_clear()
            || !current.invalidation.is_clear()
            || !successor_state.invalidation.is_clear()
        {
            return Err(AuthorityError::IllegalTransition {
                from: self.invalidation.code(),
                to: successor_state.invalidation.code(),
            });
        }
        let generation = self
            .generation
            .checked_add(1)
            .ok_or(AuthorityError::InvalidValue {
                what: "authority-head generation overflow",
            })?;
        let successor = Self::from_parts(
            self.claim,
            successor_state.identity,
            successor_state.invalidation,
            generation,
            Some(self.identity),
        );
        *self = successor;
        Ok(())
    }

    /// Atomically derive the invalidated state and advance this exact live
    /// head. The authenticated wrapper has no public constructor in Phase
    /// 0B-A; durable receipt verification owns that boundary.
    #[allow(dead_code)]
    fn revoke(
        &mut self,
        current: &AuthorityState,
        tombstone: &VerifiedRevocationTombstone,
    ) -> Result<AuthorityState, AuthorityError> {
        if self.state != current.identity || self.claim != current.claim.identity {
            return Err(AuthorityError::IdentityMismatch {
                what: "revocation predecessor is not the current exact authority head",
            });
        }
        if !self.invalidation.is_clear() || !current.invalidation.is_clear() {
            return Err(AuthorityError::IllegalTransition {
                from: "invalidated-authority-head",
                to: "invalidated-authority-head",
            });
        }
        if tombstone.target_claim != self.claim
            || tombstone.target_state != self.state
            || tombstone.target_head != self.identity
            || tombstone.target_generation != self.generation
        {
            return Err(AuthorityError::IdentityMismatch {
                what: "verified tombstone does not target the current exact authority head",
            });
        }
        let generation = self
            .generation
            .checked_add(1)
            .ok_or(AuthorityError::InvalidValue {
                what: "authority-head generation overflow",
            })?;
        let binding = InvalidationBinding::new(
            self.claim,
            self.state,
            self.identity,
            self.generation,
            tombstone.identity,
        );
        let invalidated = AuthorityState::from_axes(
            current.claim.clone(),
            current.truth,
            current.satisfiability,
            current.nonvacuity,
            ExactInstanceAdmission::NotEvaluated,
            current.kernel,
            current.scale,
            current.reproduction,
            InvalidationState::invalidated(binding),
        )?;
        let successor = Self::from_parts(
            self.claim,
            invalidated.identity,
            invalidated.invalidation,
            generation,
            Some(self.identity),
        );
        *self = successor;
        Ok(invalidated)
    }

    fn from_parts(
        claim: ClaimInstanceId,
        state: AuthorityStateId,
        invalidation: InvalidationState,
        generation: u64,
        predecessor: Option<AuthorityHeadId>,
    ) -> Self {
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *claim.as_hash());
        bytes.hash(3, *state.as_hash());
        match invalidation.binding() {
            None => bytes.u8(4, 0),
            Some(binding) => {
                bytes.u8(4, 1);
                bytes.hash(5, *binding.identity.as_hash());
            }
        }
        bytes.u64(6, generation);
        match predecessor {
            None => bytes.u8(7, 0),
            Some(previous) => {
                bytes.u8(7, 1);
                bytes.hash(8, *previous.as_hash());
            }
        }
        Self {
            claim,
            state,
            invalidation,
            generation,
            predecessor,
            identity: AuthorityHeadId(fs_blake3::hash_domain(
                AUTHORITY_HEAD_IDENTITY_DOMAIN,
                &bytes.0,
            )),
        }
    }

    #[must_use]
    pub fn claim(&self) -> ClaimInstanceId {
        self.claim
    }

    #[must_use]
    pub fn state(&self) -> AuthorityStateId {
        self.state
    }

    #[must_use]
    pub fn invalidation(&self) -> InvalidationState {
        self.invalidation
    }

    #[must_use]
    pub fn generation(&self) -> u64 {
        self.generation
    }

    #[must_use]
    pub fn predecessor(&self) -> Option<AuthorityHeadId> {
        self.predecessor
    }

    #[must_use]
    pub fn identity(&self) -> AuthorityHeadId {
        self.identity
    }
}

/// Opaque positive runtime admission.  There is no constructor from hashes or
/// booleans; minting requires a positive type-state grant, an exact policy,
/// a capability-authenticated accepting checker decision, and the current
/// authority-head generation.
#[derive(Debug, PartialEq, Eq)]
pub struct RuntimeAdmission {
    claim: ClaimInstanceId,
    authority: AuthorityStateId,
    policy: CapabilityPolicyId,
    checker_decision: CheckerDecisionId,
    authority_head: AuthorityHeadId,
    head_generation: u64,
    identity: RuntimeAdmissionId,
}

impl RuntimeAdmission {
    #[allow(dead_code)]
    fn admit<S: RuntimeTruthAuthority>(
        grant: &AuthorityGrant<S>,
        policy: &CapabilityPolicy,
        decision: &CheckerDecision,
        head: &AuthorityHead,
    ) -> Result<Self, AuthorityError> {
        let state = grant.state();
        if head.claim != state.claim.identity || head.state != state.identity {
            return Err(AuthorityError::IdentityMismatch {
                what: "runtime admission does not bind the current authority head",
            });
        }
        validate_runtime_candidate(state, policy, &decision.candidate)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.u32(2, AUTHORITY_POLICY_VERSION);
        bytes.hash(3, *state.claim.identity.as_hash());
        bytes.hash(4, *state.identity.as_hash());
        bytes.hash(5, *policy.identity.as_hash());
        bytes.hash(6, *decision.identity().as_hash());
        bytes.hash(7, *head.identity.as_hash());
        bytes.u64(8, head.generation);
        let identity = RuntimeAdmissionId(fs_blake3::hash_domain(
            RUNTIME_ADMISSION_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            claim: state.claim.identity,
            authority: state.identity,
            policy: policy.identity,
            checker_decision: decision.identity(),
            authority_head: head.identity,
            head_generation: head.generation,
            identity,
        })
    }

    #[allow(dead_code)]
    fn validate_current(&self, head: &AuthorityHead) -> Result<(), AuthorityError> {
        if self.claim == head.claim
            && self.authority == head.state
            && self.authority_head == head.identity
            && self.head_generation == head.generation
        {
            Ok(())
        } else {
            Err(AuthorityError::RuntimeRequirementNotMet {
                requirement: "current-authority-head",
            })
        }
    }

    #[must_use]
    pub fn claim(&self) -> ClaimInstanceId {
        self.claim
    }

    #[must_use]
    pub fn authority(&self) -> AuthorityStateId {
        self.authority
    }

    #[must_use]
    pub fn policy(&self) -> CapabilityPolicyId {
        self.policy
    }

    #[must_use]
    pub fn checker_decision(&self) -> CheckerDecisionId {
        self.checker_decision
    }

    #[must_use]
    pub fn authority_head(&self) -> AuthorityHeadId {
        self.authority_head
    }

    #[must_use]
    pub fn head_generation(&self) -> u64 {
        self.head_generation
    }

    #[must_use]
    pub fn identity(&self) -> RuntimeAdmissionId {
        self.identity
    }
}

#[allow(clippy::too_many_lines)]
fn validate_runtime_candidate(
    state: &AuthorityState,
    policy: &CapabilityPolicy,
    candidate: &CheckerDecisionCandidate,
) -> Result<(), AuthorityError> {
    if candidate.verdict != CheckerVerdict::Accept {
        return Err(AuthorityError::CheckerRefused {
            verdict: candidate.verdict,
        });
    }
    if candidate.claim != state.claim.identity
        || candidate.authority != state.identity
        || candidate.policy != policy.identity
    {
        return Err(AuthorityError::IdentityMismatch {
            what: "checker candidate does not bind claim, authority state, and policy exactly",
        });
    }
    if !state.invalidation.is_clear() {
        return Err(AuthorityError::RuntimeRequirementNotMet {
            requirement: "not-invalidated",
        });
    }
    let exact_decision = match state.exact_admission {
        ExactInstanceAdmission::Admitted(decision)
            if decision.verdict == ExactInstanceVerdict::Admitted =>
        {
            decision
        }
        _ => {
            return Err(AuthorityError::RuntimeRequirementNotMet {
                requirement: "exact-instance-admitted",
            });
        }
    };
    if exact_decision.claim != state.claim.identity
        || exact_decision.policy != policy.identity
        || exact_decision.checker != candidate.checker
        || exact_decision.artifact != candidate.decision_artifact
    {
        return Err(AuthorityError::IdentityMismatch {
            what: "exact-instance admission does not bind the runtime claim, policy, checker, and decision artifact exactly",
        });
    }
    match policy.truth {
        TruthRequirement::ConditionalOrProved
            if !matches!(
                state.truth,
                TruthState::ConditionalProof | TruthState::Proved
            ) =>
        {
            return Err(AuthorityError::RuntimeRequirementNotMet {
                requirement: "conditional-or-proved",
            });
        }
        TruthRequirement::ProvedOnly if state.truth != TruthState::Proved => {
            return Err(AuthorityError::RuntimeRequirementNotMet {
                requirement: "proved-only",
            });
        }
        _ => {}
    }
    if state.truth == TruthState::ConditionalProof {
        for assumption in state.claim.assumptions.assumptions() {
            if policy
                .accepted_assumptions
                .binary_search(assumption)
                .is_err()
            {
                return Err(AuthorityError::AssumptionNotAccepted {
                    assumption: assumption.clone(),
                });
            }
        }
    }
    if policy.require_satisfiable
        && !matches!(state.satisfiability, SatisfiabilityState::Satisfiable(_))
    {
        return Err(AuthorityError::RuntimeRequirementNotMet {
            requirement: "satisfiable",
        });
    }
    if let Some(required) = policy.nonvacuity_requirement {
        match state.nonvacuity {
            NonvacuityState::Nonvacuous(evidence) if evidence.strength().satisfies(&required) => {}
            _ => {
                return Err(AuthorityError::RuntimeRequirementNotMet {
                    requirement: "nonvacuity-strength",
                });
            }
        }
    }
    if policy.require_kernel_checked && !matches!(state.kernel, KernelState::KernelChecked(_)) {
        return Err(AuthorityError::RuntimeRequirementNotMet {
            requirement: "kernel-checked",
        });
    }
    if policy.require_scale_qualified && !matches!(state.scale, ScaleState::ScaleQualified(_)) {
        return Err(AuthorityError::RuntimeRequirementNotMet {
            requirement: "scale-qualified",
        });
    }
    if policy.require_reproduced && !matches!(state.reproduction, ReproductionState::Reproduced(_))
    {
        return Err(AuthorityError::RuntimeRequirementNotMet {
            requirement: "reproduced",
        });
    }
    // These are exact Five-Explicits declarations.  Actual possession remains
    // the responsibility of the authenticated Phase 0B-B checker receipt.
    for required in &policy.required_capabilities {
        if !state.claim.explicits.capabilities.contains(required) {
            return Err(AuthorityError::CapabilityMissing {
                capability: required.name.clone(),
                version: required.version,
            });
        }
    }
    for boundary in &state.claim.no_claim.entries {
        if policy.accepted_no_claims.binary_search(boundary).is_err() {
            return Err(AuthorityError::NoClaimNotAccepted {
                boundary: boundary.clone(),
            });
        }
    }
    Ok(())
}

fn validate_exact_instance_admission(
    claim: ClaimInstanceId,
    admission: ExactInstanceAdmission,
) -> Result<(), AuthorityError> {
    match admission {
        ExactInstanceAdmission::NotEvaluated => Ok(()),
        ExactInstanceAdmission::Refused(decision) => {
            if decision.claim != claim {
                return Err(AuthorityError::IdentityMismatch {
                    what: "exact-instance refusal targets another claim",
                });
            }
            if decision.verdict != ExactInstanceVerdict::Refused {
                return Err(AuthorityError::IncompatibleAxes {
                    what: "exact-instance refusal variant carries an admitted verdict",
                });
            }
            Ok(())
        }
        ExactInstanceAdmission::Admitted(decision) => {
            if decision.claim != claim {
                return Err(AuthorityError::IdentityMismatch {
                    what: "exact-instance admission targets another claim",
                });
            }
            if decision.verdict != ExactInstanceVerdict::Admitted {
                return Err(AuthorityError::IncompatibleAxes {
                    what: "exact-instance admission variant carries a refused verdict",
                });
            }
            Ok(())
        }
    }
}

fn validate_satisfiability_binding(
    claim: ClaimInstanceId,
    state: SatisfiabilityState,
) -> Result<(), AuthorityError> {
    match state {
        SatisfiabilityState::Unknown => Ok(()),
        SatisfiabilityState::Satisfiable(evidence) => {
            if evidence.evidence().claim() == claim {
                Ok(())
            } else {
                Err(AuthorityError::IdentityMismatch {
                    what: "satisfiability evidence targets another claim",
                })
            }
        }
        SatisfiabilityState::Unsatisfiable(evidence) => {
            if evidence.evidence().claim() == claim {
                Ok(())
            } else {
                Err(AuthorityError::IdentityMismatch {
                    what: "unsatisfiability evidence targets another claim",
                })
            }
        }
    }
}

fn validate_nonvacuity_binding(
    claim: ClaimInstanceId,
    state: NonvacuityState,
) -> Result<(), AuthorityError> {
    match state {
        NonvacuityState::Unknown => Ok(()),
        NonvacuityState::Nonvacuous(evidence) => {
            if evidence.evidence().claim() == claim {
                Ok(())
            } else {
                Err(AuthorityError::IdentityMismatch {
                    what: "nonvacuity evidence targets another claim",
                })
            }
        }
        NonvacuityState::Vacuous(evidence) => {
            if evidence.evidence().claim() == claim {
                Ok(())
            } else {
                Err(AuthorityError::IdentityMismatch {
                    what: "vacuity evidence targets another claim",
                })
            }
        }
    }
}

fn kernel_evidence(state: KernelState) -> Option<EvidenceRef> {
    match state {
        KernelState::NotChecked => None,
        KernelState::KernelChecked(evidence) => Some(evidence),
    }
}

fn scale_evidence(state: ScaleState) -> Option<EvidenceRef> {
    match state {
        ScaleState::NotQualified => None,
        ScaleState::ScaleQualified(evidence) => Some(evidence),
    }
}

fn reproduction_evidence(state: ReproductionState) -> Option<EvidenceRef> {
    match state {
        ReproductionState::NotAttempted => None,
        ReproductionState::Failed(evidence) => Some(evidence.evidence()),
        ReproductionState::Reproduced(evidence) => Some(evidence.evidence()),
    }
}

fn validate_reproduction_binding(
    claim: ClaimInstanceId,
    state: ReproductionState,
) -> Result<(), AuthorityError> {
    let Some(evidence) = reproduction_evidence(state) else {
        return Ok(());
    };
    if evidence.claim() == claim {
        Ok(())
    } else {
        Err(AuthorityError::IdentityMismatch {
            what: "reproduction conclusion targets another claim",
        })
    }
}

fn validate_evidence_axis(
    claim: ClaimInstanceId,
    evidence: Option<EvidenceRef>,
    required_kind: EvidenceKind,
) -> Result<(), AuthorityError> {
    let Some(evidence) = evidence else {
        return Ok(());
    };
    if evidence.claim != claim {
        return Err(AuthorityError::IdentityMismatch {
            what: "authority-axis evidence targets another claim",
        });
    }
    if evidence.kind != required_kind {
        return Err(AuthorityError::IncompatibleAxes {
            what: "authority-axis evidence has the wrong typed kind",
        });
    }
    Ok(())
}

fn encode_optional_evidence(bytes: &mut CanonicalBytes, tag: u8, evidence: Option<EvidenceRef>) {
    match evidence {
        None => bytes.u8(tag, 0),
        Some(evidence) => {
            bytes.u8(tag, 1);
            bytes.hash(tag + 1, *evidence.identity.as_hash());
        }
    }
}

fn encode_reproduction(bytes: &mut CanonicalBytes, tag: u8, state: ReproductionState) {
    match state {
        ReproductionState::NotAttempted => bytes.u8(tag, 0),
        ReproductionState::Failed(evidence) => {
            bytes.u8(tag, 1);
            bytes.hash(tag + 1, *evidence.identity().as_hash());
        }
        ReproductionState::Reproduced(evidence) => {
            bytes.u8(tag, 2);
            bytes.hash(tag + 1, *evidence.identity().as_hash());
        }
    }
}

fn admission_dominates(current: ExactInstanceAdmission, required: ExactInstanceAdmission) -> bool {
    required == ExactInstanceAdmission::NotEvaluated || current == required
}

fn optional_evidence_dominates(
    current: Option<EvidenceRef>,
    required: Option<EvidenceRef>,
) -> bool {
    required.is_none() || current == required
}

fn reproduction_dominates(current: ReproductionState, required: ReproductionState) -> bool {
    required == ReproductionState::NotAttempted || current == required
}

mod grant_sealed {
    pub trait Sealed {}
}

/// Marker for unclassified/unknown truth authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UnknownAuthority;
/// Marker for conditional-proof authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConditionalAuthority;
/// Marker for unqualified proved authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProvedAuthority;
/// Marker for refuted authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RefutedAuthority;

impl grant_sealed::Sealed for UnknownAuthority {}
impl grant_sealed::Sealed for ConditionalAuthority {}
impl grant_sealed::Sealed for ProvedAuthority {}
impl grant_sealed::Sealed for RefutedAuthority {}

/// Public marker surface with private implementations (sealed).
pub trait AuthorityMarker: grant_sealed::Sealed {}

impl AuthorityMarker for UnknownAuthority {}
impl AuthorityMarker for ConditionalAuthority {}
impl AuthorityMarker for ProvedAuthority {}
impl AuthorityMarker for RefutedAuthority {}

/// Type-state view over one validated [`AuthorityState`].  Marker transitions
/// are available only through checked `AuthorityState::*_grant` methods.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuthorityGrant<S: AuthorityMarker> {
    state: AuthorityState,
    marker: PhantomData<S>,
}

impl<S: AuthorityMarker> AuthorityGrant<S> {
    fn new(state: AuthorityState) -> Self {
        Self {
            state,
            marker: PhantomData,
        }
    }

    #[must_use]
    pub fn state(&self) -> &AuthorityState {
        &self.state
    }
}

/// Sealed marker implemented only by positive proof states.  Unknown and
/// refuted grants therefore cannot satisfy runtime APIs at compile time.
pub trait RuntimeTruthAuthority: AuthorityMarker {}

impl RuntimeTruthAuthority for ConditionalAuthority {}
impl RuntimeTruthAuthority for ProvedAuthority {}

/// Versioned inference rule.  Rule identities are data; the default algebra
/// assumes no extension rules, so a powerful new theorem cannot silently
/// become a built-in widening conversion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InferenceRule {
    name: String,
    version: u32,
    definition_artifact: ContentHash,
    identity: InferenceRuleId,
}

impl InferenceRule {
    pub fn new(
        name: &str,
        version: u32,
        definition_artifact: ContentHash,
    ) -> Result<Self, AuthorityError> {
        if version == 0 {
            return Err(AuthorityError::InvalidValue {
                what: "inference rule version must be nonzero",
            });
        }
        let name = canonical_text("inference rule", name)?;
        let definition_artifact = require_hash("inference rule artifact", definition_artifact)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.field(2, name.as_bytes());
        bytes.u32(3, version);
        bytes.hash(4, definition_artifact);
        let identity = InferenceRuleId(fs_blake3::hash_domain(
            INFERENCE_RULE_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            name,
            version,
            definition_artifact,
            identity,
        })
    }

    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[must_use]
    pub fn version(&self) -> u32 {
        self.version
    }

    #[must_use]
    pub fn definition_artifact(&self) -> ContentHash {
        self.definition_artifact
    }

    #[must_use]
    pub fn identity(&self) -> InferenceRuleId {
        self.identity
    }
}

/// The built-in extension-rule set is deliberately empty.  Consumers must
/// bind an explicit versioned [`InferenceRule`] before using one.
pub const DEFAULT_INFERENCE_RULES: &[InferenceRule] = &[];

/// Descriptive support-edge candidate over exact source-state and target-claim
/// identities. Authentication, graph acyclicity, and allocation are downstream
/// Phase 0B-C responsibilities; this type only fixes the edge algebra.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SupportEdge {
    source: AuthorityStateId,
    target: ClaimInstanceId,
    proof_lane: ProofLaneId,
    rule: InferenceRuleId,
    evidence: EvidenceId,
    identity: SupportEdgeId,
}

impl SupportEdge {
    pub fn new(
        source: &AuthorityState,
        target: &ClaimInstance,
        rule: &InferenceRule,
        evidence: EvidenceRef,
    ) -> Result<Self, AuthorityError> {
        if !source.invalidation.is_clear() {
            return Err(AuthorityError::IncompatibleAxes {
                what: "invalidated authority cannot support downstream authority",
            });
        }
        validate_evidence_axis(target.identity, Some(evidence), EvidenceKind::Support)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *source.identity.as_hash());
        bytes.hash(3, *target.identity.as_hash());
        bytes.hash(4, *target.proof_lane.as_hash());
        bytes.hash(5, *rule.identity.as_hash());
        bytes.hash(6, *evidence.identity.as_hash());
        let identity = SupportEdgeId(fs_blake3::hash_domain(
            SUPPORT_EDGE_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            source: source.identity,
            target: target.identity,
            proof_lane: target.proof_lane,
            rule: rule.identity,
            evidence: evidence.identity,
            identity,
        })
    }

    #[must_use]
    pub fn source(&self) -> AuthorityStateId {
        self.source
    }

    #[must_use]
    pub fn target(&self) -> ClaimInstanceId {
        self.target
    }

    #[must_use]
    pub fn proof_lane(&self) -> ProofLaneId {
        self.proof_lane
    }

    #[must_use]
    pub fn rule(&self) -> InferenceRuleId {
        self.rule
    }

    #[must_use]
    pub fn evidence(&self) -> EvidenceId {
        self.evidence
    }

    #[must_use]
    pub fn identity(&self) -> SupportEdgeId {
        self.identity
    }
}

/// Candidate counterexample bound to the exact claim and domain it attacks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CounterexampleCandidate {
    target: ClaimInstanceId,
    domain: QuantifiedDomainId,
    evidence: EvidenceId,
    identity: CounterexampleId,
}

impl CounterexampleCandidate {
    pub fn new(target: &ClaimInstance, evidence: EvidenceRef) -> Result<Self, AuthorityError> {
        validate_evidence_axis(
            target.identity,
            Some(evidence),
            EvidenceKind::Counterexample,
        )?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *target.identity.as_hash());
        bytes.hash(3, *target.domain.identity.as_hash());
        bytes.hash(4, *evidence.identity.as_hash());
        let identity = CounterexampleId(fs_blake3::hash_domain(
            COUNTEREXAMPLE_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            target: target.identity,
            domain: target.domain.identity,
            evidence: evidence.identity,
            identity,
        })
    }

    #[must_use]
    pub fn target(&self) -> ClaimInstanceId {
        self.target
    }

    #[must_use]
    pub fn domain(&self) -> QuantifiedDomainId {
        self.domain
    }

    #[must_use]
    pub fn evidence(&self) -> EvidenceId {
        self.evidence
    }

    #[must_use]
    pub fn identity(&self) -> CounterexampleId {
        self.identity
    }
}

/// An attack edge connects an exact candidate to its exact target.  It does
/// not itself adjudicate whether the candidate is genuine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AttackEdge {
    candidate: CounterexampleId,
    target: ClaimInstanceId,
    proof_lane: ProofLaneId,
    evidence: EvidenceId,
    identity: AttackEdgeId,
}

impl AttackEdge {
    pub fn new(
        candidate: &CounterexampleCandidate,
        target: &ClaimInstance,
        evidence: EvidenceRef,
    ) -> Result<Self, AuthorityError> {
        if candidate.target != target.identity || candidate.domain != target.domain.identity {
            return Err(AuthorityError::IdentityMismatch {
                what: "attack candidate is not bound to the target claim/domain",
            });
        }
        validate_evidence_axis(target.identity, Some(evidence), EvidenceKind::Attack)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *candidate.identity.as_hash());
        bytes.hash(3, *target.identity.as_hash());
        bytes.hash(4, *target.proof_lane.as_hash());
        bytes.hash(5, *evidence.identity.as_hash());
        let identity = AttackEdgeId(fs_blake3::hash_domain(
            ATTACK_EDGE_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            candidate: candidate.identity,
            target: target.identity,
            proof_lane: target.proof_lane,
            evidence: evidence.identity,
            identity,
        })
    }

    #[must_use]
    pub fn candidate(&self) -> CounterexampleId {
        self.candidate
    }

    #[must_use]
    pub fn target(&self) -> ClaimInstanceId {
        self.target
    }

    #[must_use]
    pub fn proof_lane(&self) -> ProofLaneId {
        self.proof_lane
    }

    #[must_use]
    pub fn evidence(&self) -> EvidenceId {
        self.evidence
    }

    #[must_use]
    pub fn identity(&self) -> AttackEdgeId {
        self.identity
    }
}

/// Adjudicated status of a counterexample candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CounterexampleVerdict {
    GenuineCounterexample,
    OutOfDomain,
    ArtifactDefect,
    Indeterminate,
}

impl CounterexampleVerdict {
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::GenuineCounterexample => "genuine-counterexample",
            Self::OutOfDomain => "out-of-domain",
            Self::ArtifactDefect => "artifact-defect",
            Self::Indeterminate => "indeterminate",
        }
    }

    fn tag(self) -> u8 {
        match self {
            Self::GenuineCounterexample => 1,
            Self::OutOfDomain => 2,
            Self::ArtifactDefect => 3,
            Self::Indeterminate => 4,
        }
    }
}

/// Immutable counterexample-adjudication candidate. A descriptive tombstone
/// can only be derived from the genuine-counterexample variant; authoritative
/// adjudication and head advancement remain receipt-gated downstream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CounterexampleAdjudication {
    candidate: CounterexampleId,
    target: ClaimInstanceId,
    verdict: CounterexampleVerdict,
    evidence: EvidenceId,
    identity: AdjudicationId,
}

impl CounterexampleAdjudication {
    pub fn new(
        candidate: &CounterexampleCandidate,
        verdict: CounterexampleVerdict,
        evidence: EvidenceRef,
    ) -> Result<Self, AuthorityError> {
        validate_evidence_axis(candidate.target, Some(evidence), EvidenceKind::Adjudication)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *candidate.identity.as_hash());
        bytes.hash(3, *candidate.target.as_hash());
        bytes.u8(4, verdict.tag());
        bytes.hash(5, *evidence.identity.as_hash());
        let identity = AdjudicationId(fs_blake3::hash_domain(
            ADJUDICATION_IDENTITY_DOMAIN,
            &bytes.0,
        ));
        Ok(Self {
            candidate: candidate.identity,
            target: candidate.target,
            verdict,
            evidence: evidence.identity,
            identity,
        })
    }

    #[must_use]
    pub fn candidate(&self) -> CounterexampleId {
        self.candidate
    }

    #[must_use]
    pub fn target(&self) -> ClaimInstanceId {
        self.target
    }

    #[must_use]
    pub fn verdict(&self) -> CounterexampleVerdict {
        self.verdict
    }

    #[must_use]
    pub fn evidence(&self) -> EvidenceId {
        self.evidence
    }

    #[must_use]
    pub fn identity(&self) -> AdjudicationId {
        self.identity
    }
}

/// Descriptive revocation/tombstone candidate identity. It binds the exact
/// state, genuine-counterexample adjudication candidate, and typed revocation
/// artifact. It cannot advance the opaque authoritative head publicly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RevocationTombstone {
    target_claim: ClaimInstanceId,
    target_state: AuthorityStateId,
    adjudication: AdjudicationId,
    reason: String,
    evidence: EvidenceId,
    identity: RevocationId,
}

impl RevocationTombstone {
    pub fn new(
        target: &AuthorityState,
        adjudication: &CounterexampleAdjudication,
        reason: &str,
        evidence: EvidenceRef,
    ) -> Result<Self, AuthorityError> {
        if adjudication.verdict != CounterexampleVerdict::GenuineCounterexample {
            return Err(AuthorityError::AdjudicationNotRevocable);
        }
        if adjudication.target != target.claim.identity {
            return Err(AuthorityError::IdentityMismatch {
                what: "counterexample adjudication targets another claim",
            });
        }
        validate_evidence_axis(
            target.claim.identity,
            Some(evidence),
            EvidenceKind::Revocation,
        )?;
        let reason = canonical_text("revocation reason", reason)?;
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *target.claim.identity.as_hash());
        bytes.hash(3, *target.identity.as_hash());
        bytes.hash(4, *adjudication.identity.as_hash());
        bytes.field(5, reason.as_bytes());
        bytes.hash(6, *evidence.identity.as_hash());
        let identity = RevocationId(fs_blake3::hash_domain(REVOCATION_IDENTITY_DOMAIN, &bytes.0));
        Ok(Self {
            target_claim: target.claim.identity,
            target_state: target.identity,
            adjudication: adjudication.identity,
            reason,
            evidence: evidence.identity,
            identity,
        })
    }

    #[must_use]
    pub fn target_claim(&self) -> ClaimInstanceId {
        self.target_claim
    }

    #[must_use]
    pub fn target_state(&self) -> AuthorityStateId {
        self.target_state
    }

    #[must_use]
    pub fn adjudication(&self) -> AdjudicationId {
        self.adjudication
    }

    #[must_use]
    pub fn reason(&self) -> &str {
        &self.reason
    }

    #[must_use]
    pub fn evidence(&self) -> EvidenceId {
        self.evidence
    }

    #[must_use]
    pub fn identity(&self) -> RevocationId {
        self.identity
    }
}

/// Receipt-authenticated revocation token. The public type is inspectable but
/// has no public constructor; Phase 0B-B durable receipt verification is the
/// only intended minting boundary. This wrapper makes no cryptographic claim
/// by itself.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VerifiedRevocationTombstone {
    target_claim: ClaimInstanceId,
    target_state: AuthorityStateId,
    target_head: AuthorityHeadId,
    target_generation: u64,
    candidate: RevocationId,
    identity: VerifiedRevocationId,
}

impl VerifiedRevocationTombstone {
    #[allow(dead_code)]
    fn from_authenticated_candidate(
        candidate: &RevocationTombstone,
        head: &AuthorityHead,
    ) -> Result<Self, AuthorityError> {
        if candidate.target_claim != head.claim || candidate.target_state != head.state {
            return Err(AuthorityError::IdentityMismatch {
                what: "revocation receipt does not target the exact live authority head",
            });
        }
        if !head.invalidation.is_clear() {
            return Err(AuthorityError::IllegalTransition {
                from: "invalidated-authority-head",
                to: "verified-revocation",
            });
        }
        let mut bytes = CanonicalBytes::default();
        bytes.u32(1, AUTHORITY_ALGEBRA_VERSION);
        bytes.hash(2, *candidate.identity.as_hash());
        bytes.hash(3, *head.identity.as_hash());
        bytes.u64(4, head.generation);
        Ok(Self {
            target_claim: candidate.target_claim,
            target_state: candidate.target_state,
            target_head: head.identity,
            target_generation: head.generation,
            candidate: candidate.identity,
            identity: VerifiedRevocationId(fs_blake3::hash_domain(
                VERIFIED_REVOCATION_IDENTITY_DOMAIN,
                &bytes.0,
            )),
        })
    }

    #[must_use]
    pub fn target_claim(&self) -> ClaimInstanceId {
        self.target_claim
    }

    #[must_use]
    pub fn target_state(&self) -> AuthorityStateId {
        self.target_state
    }

    #[must_use]
    pub fn target_head(&self) -> AuthorityHeadId {
        self.target_head
    }

    #[must_use]
    pub fn target_generation(&self) -> u64 {
        self.target_generation
    }

    #[must_use]
    pub fn candidate(&self) -> RevocationId {
        self.candidate
    }

    #[must_use]
    pub fn identity(&self) -> VerifiedRevocationId {
        self.identity
    }
}
