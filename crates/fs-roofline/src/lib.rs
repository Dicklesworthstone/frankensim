//! fs-roofline: the roofline harness (plan §14; Decalogue P6).
//!
//! Performance claims as FALSIFIABLE targets: every registered kernel is
//! benchmarked against its arithmetic-intensity-derived limit on the actual
//! machine — measured axes, never spec-sheet numbers — with dispersion
//! reported and results ledgered under the machine fingerprint. "A target
//! that was never re-measured is a lie waiting to happen."
//!
//! Layer: L6 (consumes fs-substrate probes, fs-simd primitives, and writes
//! fs-ledger records). Reporting-only in v0: attainment verdicts inform;
//! gating bands belong to nightly runs on ledgered reference machines.

pub mod authority;
pub mod axes;
pub mod baseline;
pub mod checkpoint;
pub mod kernels;
pub mod production;
pub mod stats;

#[cfg(test)]
pub use authority::StaticKeyRegistry;
pub use authority::{
    ConfiguredPromotionAuthority, KeyVerdict, NoPromotionAuthority, PromotionAttestation,
    PromotionAuthorityConfigError, PromotionAuthorityDecision, PromotionAuthorityVerifier,
};
pub use axes::MachineAxes;
pub use baseline::{
    AttestedBaselineStore, BASELINE_SCHEMA_VERSION, BaselineAxes, BaselineCandidate,
    BaselineClockError, BaselineIdentity, BaselineProvenance, BaselineStore, BaselineVerdict,
    PromotionError, candidate_axis_admission, days_since_epoch_now, promote_baseline,
};
pub use fs_blake3::ContentHash;

use fs_ledger::{EdgeRole, EventRow, FiveExplicits, Ledger, LedgerError, OpOutcome, now_wall_ns};

pub mod regress;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Semantic version of the finalized registry-run identity.
pub const FINALIZED_RUN_IDENTITY_VERSION: u32 = 3;
/// BLAKE3 derive-key domain for the finalized registry-run identity.
pub const FINALIZED_RUN_DOMAIN: &str = "org.frankensim.fs-roofline.finalized-run.v3";
/// Semantic version of the ordered result-manifest child digest.
pub const RESULT_MANIFEST_IDENTITY_VERSION: u32 = 1;
/// BLAKE3 derive-key domain for the ordered result-manifest child digest.
pub const RESULT_MANIFEST_DOMAIN: &str = "org.frankensim.fs-roofline.run-result-manifest.v1";
const RESULT_MANIFEST_SCHEMA: &str = "fs-roofline-run-manifest-v1";
/// Ledger protocol version for receipt-backed production roofline evidence.
pub const PRODUCTION_PROTOCOL_VERSION: &str = "production-v3";
pub(crate) const PRODUCTION_PROTOCOL_FIELD: &str = "\"protocol\":\"production-v3\"";
/// Ledger protocol version for exploratory/report-only roofline evidence.
pub const CUSTOM_REGISTRY_PROTOCOL_VERSION: &str = "custom-registry";
const CUSTOM_REGISTRY_PROTOCOL_FIELD: &str = "\"protocol\":\"custom-registry\"";
/// Configured durable ledger path used by external FEEC/FFT performance gates.
pub const EXTERNAL_PERF_GATE_LEDGER_ENV: &str = "FRANKENSIM_ROOFLINE_LEDGER";
/// Maximum exact final-gate JSON retained by the external gate recorder.
///
/// The operation IR embeds the JSON value as well as its content hash, so the
/// bound leaves headroom beneath fs-ledger's one-MiB operation-field ceiling.
pub const MAX_EXTERNAL_PERF_GATE_JSON_BYTES: usize = 1_000_000;
const CUSTOM_TUNE_SHAPE_CLASS: &str = "roofline-candidate-v1";
const DEPGRAPH_RECEIPT_ARTIFACT_KIND: &str = "fs-la-depgraph-receipt";
const DEPGRAPH_RECEIPT_ARTIFACT_META: &str =
    "{\"schema\":\"fs-la-depgraph-receipt-v1\",\"trust\":\"operator-observed\"}";
// Exact producer ceiling from fs-la/depgraph_receipt_format.rs. A source-pin
// test detects producer drift without adding a backwards L6 -> L1 edge.
const MAX_DEPGRAPH_RECEIPT_BYTES: u64 = 1_048_576;
const INCONSISTENT_RECEIPT_REFUSAL: &str = "dependency receipt fields compiled into fs-session are incomplete or fail their exported domain digest";

const ROOFLINE_SEED: &[u8] = b"roofline";
const ROOFLINE_BUDGET: &str = "{\"wall_s\":600}";
const ROOFLINE_CAPABILITY: &str = "{\"ops\":[\"perf.roofline\"]}";
/// Semantic version of the current executable-content identity.
pub const EXECUTABLE_BUILD_IDENTITY_VERSION: u32 = 1;
/// BLAKE3 derive-key domain for the current executable-content identity.
pub const ROOFLINE_EXECUTABLE_DOMAIN: &str = "org.frankensim.fs-roofline.executable.v1";

/// Owner-local execution-binding declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const EXECUTION_BINDING_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-roofline:execution-binding",
    "version_const=EXECUTION_BINDING_IDENTITY_VERSION",
    "version=4",
    "domain=org.frankensim.fs-roofline.execution-binding.v4",
    "domain_const=EXECUTION_BINDING_DOMAIN",
    "encoder=KernelExecutionBinding::receipt_identity",
    "encoder_helpers=KernelExecutionBinding::canonical_json,execution_binding_receipt_with_domain,execution_path_json,execution_path_identity",
    "schema_constants=EXECUTION_BINDING_IDENTITY_VERSION,EXECUTION_BINDING_DOMAIN,EXECUTION_BINDING_KIND,crates/fs-session/src/gemm_tune.rs#GEMM_TUNE_ROW_RECEIPT_IDENTITY_VERSION,crates/fs-session/src/gemm_tune.rs#GEMM_TUNE_ROW_RECEIPT_DOMAIN,crates/fs-session/src/gemm_tune.rs#GEMM_EXECUTION_RECEIPT_IDENTITY_VERSION,crates/fs-session/src/gemm_tune.rs#GEMM_EXECUTION_RECEIPT_DOMAIN",
    "schema_functions=KernelExecutionBinding::gemm,KernelExecutionBinding::is_valid_for,KernelExecutionBinding::stable_equivalent,execution_path_is_complete,execution_path_shape_eq,stable_decision_binding,crates/fs-session/src/gemm_tune.rs#ValidatedGemmTuneRow::receipt_identity,crates/fs-session/src/gemm_tune.rs#ValidatedGemmTuneRow::receipt_json,crates/fs-session/src/gemm_tune.rs#GemmExecutionReceipt::receipt_identity,crates/fs-session/src/gemm_tune.rs#GemmExecutionReceipt::canonical_bytes,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:gemm-tune-row-receipt,fs-session:gemm-execution-receipt",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=KernelExecutionBinding,GemmDecisionBinding",
    "source_fields=KernelExecutionBinding.gemm:derived:nested-decision-fields-classified-separately,GemmDecisionBinding.scoped_tune_key:semantic,GemmDecisionBinding.shape_class:semantic,GemmDecisionBinding.canonical_plan:semantic,GemmDecisionBinding.source:semantic,GemmDecisionBinding.operation_tier:semantic,GemmDecisionBinding.build_identity:semantic,GemmDecisionBinding.tune_row_identity:derived:cached-child-receipt-identity,GemmDecisionBinding.validated_row:derived:nested-fs-session-tune-row-receipt,GemmDecisionBinding.execution_path:derived:nested-fs-session-execution-receipt,GemmDecisionBinding.execution_path_identity:derived:cached-child-receipt-identity",
    "source_bindings=GemmDecisionBinding.scoped_tune_key>scoped-tune-key,GemmDecisionBinding.shape_class>shape-class,GemmDecisionBinding.canonical_plan>canonical-plan,GemmDecisionBinding.source>decision-source,GemmDecisionBinding.operation_tier>operation-tier,GemmDecisionBinding.build_identity>build-identity",
    "external_semantic_fields=digest-domain,identity-version,kind-tag,tune-row-receipt-child,execution-path-receipt-child",
    "semantic_fields=digest-domain,identity-version,kind-tag,scoped-tune-key,shape-class,canonical-plan,decision-source,operation-tier,build-identity,tune-row-receipt-child,execution-path-receipt-child",
    "excluded_fields=none",
    "consumers=Attainment::to_jsonl,Attainment::is_citable_against,stable_decision_binding,run_admission_error",
    "mutations=digest-domain:crates/fs-roofline/src/lib.rs#execution_binding_identity_versions_fail_closed,identity-version:crates/fs-roofline/src/lib.rs#execution_binding_identity_versions_fail_closed,kind-tag:crates/fs-roofline/src/lib.rs#execution_binding_identity_versions_fail_closed,scoped-tune-key:crates/fs-roofline/src/kernels.rs#citable_gemm_receipt_rejects_every_bound_field_tamper,shape-class:crates/fs-roofline/src/kernels.rs#citable_gemm_receipt_rejects_every_bound_field_tamper,canonical-plan:crates/fs-roofline/src/kernels.rs#citable_gemm_receipt_rejects_every_bound_field_tamper,decision-source:crates/fs-roofline/src/kernels.rs#citable_gemm_receipt_rejects_every_bound_field_tamper,operation-tier:crates/fs-roofline/src/kernels.rs#citable_gemm_receipt_rejects_every_bound_field_tamper,build-identity:crates/fs-roofline/src/kernels.rs#citable_gemm_receipt_rejects_every_bound_field_tamper,tune-row-receipt-child:crates/fs-roofline/src/kernels.rs#citable_gemm_receipt_rejects_every_bound_field_tamper,execution-path-receipt-child:crates/fs-roofline/src/kernels.rs#citable_gemm_receipt_rejects_every_bound_field_tamper",
    "nonsemantic_mutations=none",
    "field_guard=classify_execution_binding_identity_fields",
    "transport_guard=KernelExecutionBinding::is_valid_for",
    "version_guard=crates/fs-roofline/src/lib.rs#execution_binding_identity_versions_fail_closed",
    "coupling_surface=fs-roofline:execution-binding",
];

/// Owner-local finalized-run declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const FINALIZED_RUN_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-roofline:finalized-run",
    "version_const=FINALIZED_RUN_IDENTITY_VERSION",
    "version=3",
    "domain=org.frankensim.fs-roofline.finalized-run.v3",
    "domain_const=FINALIZED_RUN_DOMAIN",
    "encoder=finalized_run_receipt",
    "encoder_helpers=FinalizedRunIdentityInput::from_run,finalized_run_receipt_from_input,push_receipt_field,run_result_manifest_json,manifest_entry_json",
    "schema_constants=FINALIZED_RUN_IDENTITY_VERSION,FINALIZED_RUN_DOMAIN,RESULT_MANIFEST_IDENTITY_VERSION,RESULT_MANIFEST_DOMAIN,RESULT_MANIFEST_SCHEMA",
    "schema_functions=parse_result_manifest,valid_manifest_identifier,receipt_recomputes_from_stored_rows,Attainment::to_jsonl,AxisAdmissionSnapshot::receipt_json,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#hash_bytes",
    "schema_dependencies=fs-roofline:baseline-record,fs-roofline:promotion-authority-policy,fs-roofline:production-axes-receipt,fs-roofline:execution-binding",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=FinalizedRunIdentityInput",
    "source_fields=FinalizedRunIdentityInput.baseline_receipt:semantic,FinalizedRunIdentityInput.result_payloads:semantic,FinalizedRunIdentityInput.result_manifest:semantic",
    "source_bindings=FinalizedRunIdentityInput.baseline_receipt>baseline-receipt-bytes,FinalizedRunIdentityInput.result_payloads>result-count+ordered-result-payloads,FinalizedRunIdentityInput.result_manifest>result-manifest-bytes",
    "external_semantic_fields=digest-domain,identity-version,receipt-length-prefixes,result-manifest-domain,result-manifest-version",
    "semantic_fields=digest-domain,identity-version,receipt-length-prefixes,baseline-receipt-bytes,result-count,ordered-result-payloads,result-manifest-domain,result-manifest-version,result-manifest-bytes",
    "excluded_fields=none",
    "consumers=FinalizedRegistryRun::receipt_identity,finalize_registry_tuning,record_run_with_protocol,receipt_recomputes_from_stored_rows",
    "mutations=digest-domain:crates/fs-roofline/src/lib.rs#finalized_run_identity_fields_move_independently,identity-version:crates/fs-roofline/src/lib.rs#finalized_run_identity_versions_fail_closed,receipt-length-prefixes:crates/fs-roofline/src/lib.rs#finalized_run_identity_fields_move_independently,baseline-receipt-bytes:crates/fs-roofline/src/lib.rs#finalized_run_identity_fields_move_independently,result-count:crates/fs-roofline/src/lib.rs#finalized_run_identity_fields_move_independently,ordered-result-payloads:crates/fs-roofline/src/lib.rs#finalized_run_identity_fields_move_independently,result-manifest-domain:crates/fs-roofline/src/lib.rs#finalized_run_identity_fields_move_independently,result-manifest-version:crates/fs-roofline/src/lib.rs#finalized_run_identity_versions_fail_closed,result-manifest-bytes:crates/fs-roofline/src/lib.rs#finalized_run_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_finalized_run_identity_fields",
    "transport_guard=receipt_recomputes_from_stored_rows",
    "version_guard=crates/fs-roofline/src/lib.rs#finalized_run_identity_versions_fail_closed",
    "coupling_surface=fs-roofline:finalized-run",
];

/// Owner-local executable-build declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const EXECUTABLE_BUILD_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-roofline:executable-build",
    "version_const=EXECUTABLE_BUILD_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-roofline.executable.v1",
    "domain_const=ROOFLINE_EXECUTABLE_DOMAIN",
    "encoder=read_executable_build_identity",
    "encoder_helpers=executable_build_identity_from_input",
    "schema_constants=EXECUTABLE_BUILD_IDENTITY_VERSION,ROOFLINE_EXECUTABLE_DOMAIN",
    "schema_functions=executable_build_identity,require_stable_executable_identity,crates/fs-blake3/src/lib.rs#Blake3::new,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#Blake3::finalize,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=ExecutableBuildIdentityInput",
    "source_fields=ExecutableBuildIdentityInput.byte_len:semantic,ExecutableBuildIdentityInput.raw_hash:semantic",
    "source_bindings=ExecutableBuildIdentityInput.byte_len>executable-byte-count,ExecutableBuildIdentityInput.raw_hash>raw-executable-content-hash",
    "external_semantic_fields=digest-domain,identity-version,length-prefix-layout",
    "semantic_fields=digest-domain,identity-version,length-prefix-layout,executable-byte-count,raw-executable-content-hash",
    "excluded_fields=executable-path:path-location-is-not-content,read-chunk-size:streaming-implementation-only",
    "consumers=versions_json,record_run_with_protocol,staleness_at,production_op_envelope_is_valid",
    "mutations=digest-domain:crates/fs-roofline/src/lib.rs#executable_build_identity_fields_move_independently,identity-version:crates/fs-roofline/src/lib.rs#executable_build_identity_versions_fail_closed,length-prefix-layout:crates/fs-roofline/src/lib.rs#executable_build_identity_fields_move_independently,executable-byte-count:crates/fs-roofline/src/lib.rs#executable_build_identity_fields_move_independently,raw-executable-content-hash:crates/fs-roofline/src/lib.rs#executable_build_identity_fields_move_independently",
    "nonsemantic_mutations=executable-path:crates/fs-roofline/src/lib.rs#executable_build_identity_excludes_path_and_chunking,read-chunk-size:crates/fs-roofline/src/lib.rs#executable_build_identity_excludes_path_and_chunking",
    "field_guard=classify_executable_build_identity_fields",
    "transport_guard=require_stable_executable_identity",
    "version_guard=crates/fs-roofline/src/lib.rs#executable_build_identity_versions_fail_closed",
    "coupling_surface=fs-roofline:executable-build",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EvidenceNamespace {
    Production,
    Custom,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct DependencyReceiptBinding {
    bytes: &'static str,
    domain_digest: fs_blake3::ContentHash,
    artifact_hash: fs_blake3::ContentHash,
}

impl DependencyReceiptBinding {
    pub(crate) fn from_parts(
        bytes: &'static str,
        domain_digest: fs_blake3::ContentHash,
    ) -> Option<Self> {
        (fs_blake3::hash_domain(fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN, bytes.as_bytes())
            == domain_digest)
            .then(|| Self {
                bytes,
                domain_digest,
                artifact_hash: fs_ledger::hash_bytes(bytes.as_bytes()),
            })
    }

    pub(crate) fn from_build_evidence(
        evidence: fs_session::GemmTuneBuildEvidence,
    ) -> Result<Self, &'static str> {
        if evidence.graph_class != fs_session::GemmGraphEvidenceClass::OperatorObservedReceipt {
            return Err(
                "dependency graph uses the development equivalence salt; production citation requires an exact operator-observed normal/build receipt",
            );
        }
        let receipt = evidence
            .dependency_receipt
            .ok_or(INCONSISTENT_RECEIPT_REFUSAL)?;
        let raw_digest = evidence
            .dependency_receipt_digest
            .ok_or(INCONSISTENT_RECEIPT_REFUSAL)?;
        let digest =
            fs_blake3::ContentHash::from_hex(raw_digest).ok_or(INCONSISTENT_RECEIPT_REFUSAL)?;
        let binding = Self::from_parts(receipt, digest).ok_or(INCONSISTENT_RECEIPT_REFUSAL)?;
        let expected_identity = format!("receipt:{digest}");
        if evidence.graph_class_identity != expected_identity {
            return Err(INCONSISTENT_RECEIPT_REFUSAL);
        }
        Ok(binding)
    }

    pub(crate) fn current() -> Result<Self, &'static str> {
        Self::from_build_evidence(fs_session::gemm_tune_build_evidence())
    }
}

/// One-shot capability proving that every registry kernel observed the exact
/// run's aggregate admission decision and completed its tuning lifecycle.
/// Fields are private and the value is deliberately neither `Clone` nor `Copy`.
#[derive(Debug)]
pub struct FinalizedRegistryRun {
    receipt: fs_blake3::ContentHash,
    admitted: bool,
    consumed: bool,
}

impl FinalizedRegistryRun {
    /// Whether the finalized result set passed aggregate measurement
    /// admission. Citation additionally requires the sealed production
    /// protocol; a custom-registry result can pass this numerical gate without
    /// acquiring publication authority.
    #[must_use]
    pub fn admitted(&self) -> bool {
        self.admitted
    }

    /// Content identity of the exact axes, baseline decision, and result set
    /// seen during finalization.
    #[must_use]
    pub fn receipt_identity(&self) -> fs_blake3::ContentHash {
        self.receipt
    }
}

/// Operator-trusted historical-baseline policy for one REPORT-ONLY run.
///
/// `None` is not a permissive default: it represents a first/unbaselined run.
/// Even when this policy's numerical verdict is [`BaselineVerdict::Trusted`],
/// it carries no promotion-authority proof and therefore cannot enter the
/// citable production protocol. Use an [`AttestedAxisBaselinePolicy`] minted
/// by [`AttestedBaselineStore::policy_for_run`] for that protocol.
#[derive(Clone, Copy)]
pub struct AxisBaselinePolicy<'a> {
    baseline: Option<&'a BaselineAxes>,
    identity: &'a BaselineIdentity,
    now_day: u64,
}

impl<'a> AxisBaselinePolicy<'a> {
    /// Bind the selected baseline (if any), declared current environment, and
    /// observed epoch day used for age admission.
    #[must_use]
    pub const fn new(
        baseline: Option<&'a BaselineAxes>,
        identity: &'a BaselineIdentity,
        now_day: u64,
    ) -> Self {
        Self {
            baseline,
            identity,
            now_day,
        }
    }

    /// Evaluate the pre/probe/post baseline math without promotion authority.
    #[must_use]
    pub fn verdict(&self, pre: &MachineAxes, post: &MachineAxes) -> BaselineVerdict {
        candidate_axis_admission(pre, post, self.baseline, self.identity, self.now_day)
    }

    /// Domain-separated identity of the selected baseline, if one exists.
    #[must_use]
    pub fn baseline_hash(&self) -> Option<fs_blake3::ContentHash> {
        self.baseline.map(BaselineAxes::content_hash)
    }

    /// Canonical, self-contained receipt for the baseline admission decision.
    #[must_use]
    pub fn receipt_json(&self, pre: &MachineAxes, post: &MachineAxes) -> String {
        self.snapshot(pre, post).receipt_json().to_string()
    }

    /// Freeze this operator-trusted decision into the same immutable receipt
    /// shape used by the production protocol, but with tier `candidate`.
    #[must_use]
    pub fn snapshot(&self, pre: &MachineAxes, post: &MachineAxes) -> AxisAdmissionSnapshot {
        AxisAdmissionSnapshot::candidate(
            self.baseline,
            self.identity,
            self.now_day,
            pre,
            post,
            self.verdict(pre, post),
        )
    }
}

/// An owned, authority-attested baseline policy for exactly one production
/// run. Fields are private and there is no public constructor: the attested
/// store checks source retention and captures one atomic authority decision
/// before minting this value.
pub struct AttestedAxisBaselinePolicy {
    baseline: BaselineAxes,
    identity: BaselineIdentity,
    now_day: u64,
    attestation: PromotionAttestation,
    source_receipts: Vec<fs_blake3::ContentHash>,
    authority_decision: PromotionAuthorityDecision,
}

impl AttestedAxisBaselinePolicy {
    /// Crate-internal minting boundary used by `AttestedBaselineStore` after
    /// it has checked attestation structure, exact authority, and source
    /// receipt availability.
    #[must_use]
    pub(crate) fn from_verified(
        baseline: BaselineAxes,
        identity: BaselineIdentity,
        now_day: u64,
        attestation: PromotionAttestation,
        source_receipts: Vec<fs_blake3::ContentHash>,
        authority_decision: PromotionAuthorityDecision,
    ) -> Self {
        Self {
            baseline,
            identity,
            now_day,
            attestation,
            source_receipts,
            authority_decision,
        }
    }

    /// Domain-separated identity of the exact selected baseline.
    #[must_use]
    pub fn baseline_hash(&self) -> fs_blake3::ContentHash {
        self.baseline.content_hash()
    }

    /// Consume the one-run policy and freeze the verifier decision together
    /// with both observed probes. No live verifier survives this boundary.
    #[must_use]
    pub fn decide(self, pre: &MachineAxes, post: &MachineAxes) -> AxisAdmissionSnapshot {
        let decision_day = days_since_epoch_now().ok();
        self.decide_on_day(pre, post, decision_day)
    }

    fn decide_on_day(
        self,
        pre: &MachineAxes,
        post: &MachineAxes,
        decision_day: Option<u64>,
    ) -> AxisAdmissionSnapshot {
        let sources_match =
            self.source_receipts.as_slice() == self.baseline.provenance().source_receipts();
        let authority_verdict = self.authority_decision.verdict();
        let attestation_well_formed = self.attestation.well_formed();
        let verdict = if decision_day.is_none() {
            BaselineVerdict::Unauthorized {
                verdict: "clock-unavailable",
            }
        } else if decision_day != Some(self.now_day) {
            BaselineVerdict::Unauthorized {
                verdict: "policy-day-mismatch",
            }
        } else if !attestation_well_formed {
            BaselineVerdict::Unauthorized {
                verdict: "malformed-attestation",
            }
        } else if authority_verdict != KeyVerdict::Authorized {
            BaselineVerdict::Unauthorized {
                verdict: authority_verdict.name(),
            }
        } else if !sources_match {
            BaselineVerdict::InvalidBaseline {
                reason:
                    "attested policy source receipts differ from the canonical baseline provenance"
                        .to_string(),
            }
        } else {
            candidate_axis_admission(
                pre,
                post,
                Some(&self.baseline),
                &self.identity,
                self.now_day,
            )
        };
        AxisAdmissionSnapshot::attested(
            self.baseline,
            self.identity,
            self.now_day,
            self.attestation,
            self.source_receipts,
            self.authority_decision,
            decision_day,
            pre,
            post,
            verdict,
        )
    }

    #[cfg(test)]
    fn decide_at(
        self,
        pre: &MachineAxes,
        post: &MachineAxes,
        decision_day: u64,
    ) -> AxisAdmissionSnapshot {
        self.decide_on_day(pre, post, Some(decision_day))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AxisAdmissionTier {
    Candidate,
    Attested,
}

/// Immutable, owned admission evidence for one exact pre/post probe pair.
///
/// The canonical bytes are created once. Finalization, eligibility, finalized
/// run identity, and ledger recording all reuse those bytes verbatim.
#[derive(Debug, Clone, PartialEq)]
pub struct AxisAdmissionSnapshot {
    tier: AxisAdmissionTier,
    receipt: String,
    baseline_hash: Option<fs_blake3::ContentHash>,
    verdict: BaselineVerdict,
    authority_admitted: bool,
    decision_day: Option<u64>,
}

impl AxisAdmissionSnapshot {
    fn candidate(
        baseline: Option<&BaselineAxes>,
        identity: &BaselineIdentity,
        now_day: u64,
        pre: &MachineAxes,
        post: &MachineAxes,
        verdict: BaselineVerdict,
    ) -> Self {
        let baseline_json =
            baseline.map_or_else(|| "null".to_string(), BaselineAxes::canonical_json);
        let baseline_hash = baseline.map(BaselineAxes::content_hash);
        let baseline_hash_json =
            baseline_hash.map_or_else(|| "null".to_string(), |hash| format!("\"{hash}\""));
        let required_sources = match baseline {
            Some(record) => record.provenance().source_receipts(),
            None => &[],
        };
        let receipt = format!(
            "{{\"schema\":\"fs-roofline-axis-admission-v2\",\"tier\":\"candidate\",\"now_day\":{now_day},\"decision_day\":null,\"identity\":{},\"pre\":{},\"post\":{},\"baseline_hash\":{baseline_hash_json},\"baseline\":{baseline_json},\"attestation\":null,\"required_source_receipts\":{},\"authority\":{{\"verdict\":\"not-attested\",\"policy_receipt\":null}},\"verdict\":{}}}",
            baseline_identity_json(identity),
            machine_axes_receipt_json(pre),
            machine_axes_receipt_json(post),
            source_receipts_json(required_sources),
            verdict.to_jsonl(),
        );
        Self {
            tier: AxisAdmissionTier::Candidate,
            receipt,
            baseline_hash,
            verdict,
            authority_admitted: false,
            decision_day: None,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn attested(
        baseline: BaselineAxes,
        identity: BaselineIdentity,
        now_day: u64,
        attestation: PromotionAttestation,
        source_receipts: Vec<fs_blake3::ContentHash>,
        authority_decision: PromotionAuthorityDecision,
        decision_day: Option<u64>,
        pre: &MachineAxes,
        post: &MachineAxes,
        verdict: BaselineVerdict,
    ) -> Self {
        let baseline_hash = baseline.content_hash();
        let authority_verdict = authority_decision.verdict();
        let authority_admitted = attestation.well_formed()
            && authority_verdict == KeyVerdict::Authorized
            && decision_day == Some(now_day)
            && source_receipts.as_slice() == baseline.provenance().source_receipts();
        let decision_day_json =
            decision_day.map_or_else(|| "null".to_string(), |day| day.to_string());
        let receipt = format!(
            "{{\"schema\":\"fs-roofline-axis-admission-v2\",\"tier\":\"attested\",\"now_day\":{now_day},\"decision_day\":{decision_day_json},\"identity\":{},\"pre\":{},\"post\":{},\"baseline_hash\":\"{baseline_hash}\",\"baseline\":{},\"attestation\":{{\"key_id\":\"{}\",\"signature\":\"{}\"}},\"required_source_receipts\":{},\"authority\":{{\"verdict\":\"{}\",\"policy_receipt\":\"{}\"}},\"verdict\":{}}}",
            baseline_identity_json(&identity),
            machine_axes_receipt_json(pre),
            machine_axes_receipt_json(post),
            baseline.canonical_json(),
            json_escape(attestation.key_id()),
            json_escape(attestation.signature()),
            source_receipts_json(&source_receipts),
            authority_verdict.name(),
            authority_decision.policy_receipt(),
            verdict.to_jsonl(),
        );
        Self {
            tier: AxisAdmissionTier::Attested,
            receipt,
            baseline_hash: Some(baseline_hash),
            verdict,
            authority_admitted,
            decision_day,
        }
    }

    /// Exact canonical receipt bytes retained at the decision boundary.
    #[must_use]
    pub fn receipt_json(&self) -> &str {
        &self.receipt
    }

    /// Final baseline verdict frozen into this snapshot.
    #[must_use]
    pub fn verdict(&self) -> &BaselineVerdict {
        &self.verdict
    }

    /// Selected baseline identity, if this run selected one.
    #[must_use]
    pub fn baseline_hash(&self) -> Option<fs_blake3::ContentHash> {
        self.baseline_hash
    }

    /// True only when an attested tier carried a well-formed attestation, an
    /// Authorized atomic authority decision, and the retained source list
    /// matched canonical provenance.
    #[must_use]
    pub fn authority_admitted(&self) -> bool {
        self.tier == AxisAdmissionTier::Attested && self.authority_admitted
    }

    /// Whether this frozen baseline decision still supports citation now.
    ///
    /// This combines the attested tier, the trusted numerical verdict, and a
    /// fresh live epoch-day check. It is necessary but not sufficient for a
    /// citable measurement: kernel binding, timing, and recording provenance
    /// remain separate obligations.
    #[must_use]
    pub fn baseline_citation_eligible(&self) -> bool {
        self.baseline_citation_error().is_none()
    }

    /// Structured reason this snapshot cannot support baseline citation now.
    #[must_use]
    pub fn baseline_citation_error(&self) -> Option<String> {
        let live_day_error = self.live_day_error();
        self.baseline_citation_error_with_live_day(live_day_error)
    }

    fn baseline_citation_error_with_live_day(
        &self,
        live_day_error: Option<String>,
    ) -> Option<String> {
        if !self.verdict.trusted() {
            return Some(format!(
                "historical baseline admission refused: {}",
                self.verdict.to_jsonl()
            ));
        }
        if !self.authority_admitted() {
            return Some(
                "baseline admission snapshot lacks authorized promotion authority".to_string(),
            );
        }
        live_day_error
    }

    /// Unix-epoch day observed when the attested policy was consumed. A
    /// candidate snapshot has no decision day.
    #[must_use]
    pub fn decision_day(&self) -> Option<u64> {
        self.decision_day
    }

    pub(crate) fn live_day_error(&self) -> Option<String> {
        if self.tier != AxisAdmissionTier::Attested {
            return None;
        }
        let Some(decision_day) = self.decision_day else {
            return Some(
                "attested baseline admission could not establish a decision day".to_string(),
            );
        };
        match days_since_epoch_now() {
            Ok(current_day) => self.day_error_at(current_day),
            Err(error) => Some(format!(
                "cannot revalidate attested baseline admission day: {error}"
            )),
        }
    }

    fn day_error_at(&self, current_day: u64) -> Option<String> {
        let decision_day = self.decision_day?;
        (current_day != decision_day).then(|| {
            format!(
                "attested baseline admission was frozen on day {decision_day}, but current day is {current_day}"
            )
        })
    }
}

fn source_receipts_json(receipts: &[fs_blake3::ContentHash]) -> String {
    let entries = receipts
        .iter()
        .map(|receipt| format!("\"{receipt}\""))
        .collect::<Vec<_>>()
        .join(",");
    format!("[{entries}]")
}

fn baseline_identity_json(identity: &BaselineIdentity) -> String {
    format!(
        "{{\"fingerprint\":\"{:016x}\",\"cpu_brand\":\"{}\",\"logical_cpus\":{},\"os\":\"{}\",\"arch\":\"{}\",\"firmware\":\"{}\"}}",
        identity.fingerprint(),
        json_escape(identity.cpu_brand()),
        identity.logical_cpus(),
        json_escape(identity.os()),
        json_escape(identity.arch()),
        json_escape(identity.firmware()),
    )
}

fn machine_axes_receipt_json(axes: &MachineAxes) -> String {
    production::machine_axes_receipt_json(axes)
}

/// Shape-class prefix under which versioned roofline rows land in the ledger
/// `tune` table.
pub const TUNE_SHAPE_CLASS: &str = "roofline-v8";

/// Versioned ledger shape-class key for a kernel implementation.
#[must_use]
pub fn tune_shape_class(version: &str) -> String {
    format!("{TUNE_SHAPE_CLASS}:{version}")
}

/// Append-only shape key for one finalized measurement operation.
#[must_use]
pub fn tune_measurement_shape_class(
    version: &str,
    run_receipt: fs_blake3::ContentHash,
    op: i64,
) -> String {
    format!("{}:run={run_receipt}:op={op}", tune_shape_class(version))
}

/// Append-only shape key for exploratory custom-registry evidence.
///
/// This namespace is intentionally disjoint from [`TUNE_SHAPE_CLASS`], so a
/// caller-controlled row can neither satisfy nor poison production staleness.
#[must_use]
pub fn candidate_measurement_shape_class(
    version: &str,
    run_receipt: fs_blake3::ContentHash,
    op: i64,
) -> String {
    format!("{CUSTOM_TUNE_SHAPE_CLASS}:{version}:run={run_receipt}:op={op}")
}

/// Composite machine key for baseline-bound roofline rows.
#[must_use]
pub fn roofline_machine_key(fingerprint: u64, baseline: fs_blake3::ContentHash) -> [u8; 40] {
    let mut key = [0_u8; 40];
    key[..8].copy_from_slice(&fingerprint.to_le_bytes());
    key[8..].copy_from_slice(baseline.as_bytes());
    key
}

/// Which machine axis a kernel is measured against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Threading {
    /// One thread: per-core bandwidth/compute axes.
    SingleThread,
    /// All logical cores: aggregate axes.
    AllCore,
}

/// Machine quantity against which a declared target fraction is evaluated.
///
/// This is independent of [`RoofSide`]: the binding roof remains part of every
/// report, while a plan target such as "75% of measured peak FLOPs" must stay
/// compute-relative even on a machine where memory bandwidth binds first.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TargetAxis {
    /// Fraction of `min(bandwidth limit, compute limit)`.
    BindingRoof,
    /// Fraction of the measured floating-point compute peak.
    ComputePeak,
    /// Fraction of the measured memory-bandwidth peak.
    MemoryBandwidth,
}

impl TargetAxis {
    /// Stable receipt spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::BindingRoof => "binding_roof",
            Self::ComputePeak => "compute_peak",
            Self::MemoryBandwidth => "memory_bandwidth",
        }
    }
}

/// Static description of a benchmarkable kernel: identity plus the
/// arithmetic-intensity model that derives its machine-specific limit.
#[derive(Debug, Clone, Copy)]
pub struct KernelSpec {
    /// Registry name (ledger key; kebab-case).
    pub name: &'static str,
    /// Kernel version (bumped when the implementation changes — attainment
    /// history is only comparable within one version).
    pub version: &'static str,
    /// Bytes moved to/from memory per element processed.
    pub bytes_per_elem: f64,
    /// Floating-point operations per element processed.
    pub flops_per_elem: f64,
    /// Measurement threading model.
    pub threading: Threading,
    /// Axis whose measured value defines `target_fraction`.
    pub target_axis: TargetAxis,
    /// Target as a fraction of [`KernelSpec::target_axis`]. `None` =
    /// report-only, no band claimed.
    pub target_fraction: Option<f64>,
}

/// One benchmarkable kernel: owns its buffers; `run_once` is the timed unit.
pub trait RooflineKernel {
    /// The kernel's spec (identity + intensity model + target).
    fn spec(&self) -> KernelSpec;
    /// Elements processed per `run_once` call.
    fn elements(&self) -> usize;
    /// Execute one timed repetition.
    ///
    /// # Errors
    /// Returns a structured diagnostic when the kernel cannot complete the
    /// repetition. A failed repetition never produces an attainment row.
    fn run_once(&mut self) -> Result<(), String>;
    /// Sealed execution decision made by the most recent `run_once`, when this
    /// kernel has decision provenance beyond its static [`KernelSpec`].
    ///
    /// The default is `None`; the production GEMM route returns a binding that
    /// includes its exact tune key, plan, source, implementation tier, build,
    /// and validated tune-row identity.
    fn execution_binding(&self) -> Option<KernelExecutionBinding> {
        None
    }
    /// Newly measured, validated tune state awaiting the enclosing evidence
    /// transaction. Returning it does not publish it.
    fn pending_tune_publication(&self) -> Option<fs_session::ValidatedGemmTuneRow> {
        None
    }
    /// Complete process-local tuning state after the enclosing registry run's
    /// aggregate admission decision. The measured [`Attainment`] already owns
    /// any pending publication clone; kernels discard their pending marker on
    /// either outcome and additionally invalidate rejected local decisions.
    ///
    /// # Errors
    /// Returns a structured diagnostic if lifecycle cleanup fails. The default
    /// is a no-op for kernels without tune state.
    fn finalize_tuning(&mut self, _admitted: bool) -> Result<(), String> {
        Ok(())
    }
    /// Discard process-local tuning state after an incomplete registry run or
    /// failed lifecycle finalization.
    ///
    /// This hook must be idempotent and must not publish durable evidence. It
    /// is called for every registry member after any peer fails so a later
    /// reuse cannot inherit authority from a run that produced no token.
    ///
    /// # Errors
    /// Returns a structured diagnostic if local cleanup fails. The registry
    /// still attempts every remaining abort hook before reporting errors.
    fn abort_tuning(&mut self) -> Result<(), String> {
        Ok(())
    }
}

/// A sealed execution-decision receipt sampled after one timed repetition.
///
/// Fields are deliberately private. Built-in kernels construct this only from
/// producer-validated values; consumers can inspect the canonical JSON emitted
/// by [`Attainment::to_jsonl`] but cannot mint citable bindings directly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KernelExecutionBinding {
    gemm: GemmDecisionBinding,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct GemmDecisionBinding {
    scoped_tune_key: String,
    shape_class: String,
    canonical_plan: String,
    source: &'static str,
    operation_tier: String,
    build_identity: String,
    tune_row_identity: fs_ledger::ContentHash,
    validated_row: fs_session::ValidatedGemmTuneRow,
    execution_path: fs_session::GemmExecutionReceipt,
    execution_path_identity: fs_ledger::ContentHash,
}

/// Semantic version of the complete GEMM decision/execution binding.
pub const EXECUTION_BINDING_IDENTITY_VERSION: u32 = 4;
/// BLAKE3 derive-key domain for the complete GEMM decision/execution binding.
pub const EXECUTION_BINDING_DOMAIN: &str = "org.frankensim.fs-roofline.execution-binding.v4";
const EXECUTION_BINDING_KIND: &str = "gemm-v4";

#[allow(dead_code)]
fn classify_execution_binding_identity_fields(
    binding: &KernelExecutionBinding,
    gemm_source: &GemmDecisionBinding,
) {
    let KernelExecutionBinding { gemm } = binding;
    let GemmDecisionBinding {
        scoped_tune_key,
        shape_class,
        canonical_plan,
        source,
        operation_tier,
        build_identity,
        tune_row_identity,
        validated_row,
        execution_path,
        execution_path_identity,
    } = gemm_source;
    let _ = (
        gemm,
        scoped_tune_key,
        shape_class,
        canonical_plan,
        source,
        operation_tier,
        build_identity,
        tune_row_identity,
        validated_row,
        execution_path,
        execution_path_identity,
    );
}

fn execution_path_is_complete(path: &fs_session::GemmExecutionReceipt) -> bool {
    path.is_complete()
}

fn execution_path_json(path: &fs_session::GemmExecutionReceipt) -> String {
    let panels = path
        .panels
        .iter()
        .map(|panel| {
            format!(
                "{{\"kernel\":\"{}\",\"mode\":\"{}\",\"declared_run\":{},\"completed\":{},\"total\":{}}}",
                json_escape(&panel.kernel),
                json_escape(&panel.mode),
                panel.declared_run,
                panel.completed,
                panel.total,
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"declared_run\":{},\"completed_tiles\":{},\"total_tiles\":{},\"panel_count\":{},\"panels\":[{}]}}",
        path.declared_run,
        path.completed_tiles,
        path.total_tiles,
        path.panels.len(),
        panels,
    )
}

fn execution_path_identity(
    path: &fs_session::GemmExecutionReceipt,
) -> Result<fs_ledger::ContentHash, String> {
    path.receipt_identity()
        .map_err(|error| format!("cannot bind the complete GEMM execution receipt: {error}"))
}

fn execution_binding_receipt_with_domain(
    canonical_json: &str,
    domain: &str,
) -> fs_ledger::ContentHash {
    fs_blake3::hash_domain(domain, canonical_json.as_bytes())
}

impl KernelExecutionBinding {
    #[allow(clippy::too_many_arguments)] // every independent provenance field stays explicit
    pub(crate) fn gemm(
        scoped_tune_key: String,
        shape_class: String,
        canonical_plan: String,
        source: &'static str,
        operation_tier: String,
        build_identity: String,
        validated_row: fs_session::ValidatedGemmTuneRow,
        machine: u64,
        execution_path: fs_session::GemmExecutionReceipt,
    ) -> Result<Self, String> {
        if !validated_row.matches_decision(&scoped_tune_key, &shape_class, machine, &canonical_plan)
        {
            return Err(
                "validated GEMM tune row does not bind the dispatched decision".to_string(),
            );
        }
        let tune_row_identity = validated_row.receipt_identity();
        let execution_path_identity = execution_path_identity(&execution_path)?;
        let binding = Self {
            gemm: GemmDecisionBinding {
                scoped_tune_key,
                shape_class,
                canonical_plan,
                source,
                operation_tier,
                build_identity,
                tune_row_identity,
                validated_row,
                execution_path,
                execution_path_identity,
            },
        };
        if !binding.is_valid_for(machine) {
            return Err("GEMM execution binding is internally inconsistent".to_string());
        }
        Ok(binding)
    }

    fn is_valid_for(&self, machine: u64) -> bool {
        let gemm = &self.gemm;
        !gemm.scoped_tune_key.is_empty()
            && !gemm.shape_class.is_empty()
            && !gemm.canonical_plan.is_empty()
            && gemm.source == "tuned"
            && !gemm.operation_tier.is_empty()
            && !gemm.build_identity.is_empty()
            && gemm
                .scoped_tune_key
                .contains(&format!("/shape={}/", gemm.shape_class))
            && gemm
                .scoped_tune_key
                .contains(&format!("/tier={}/", gemm.operation_tier))
            && gemm
                .scoped_tune_key
                .ends_with(&format!("/build={}", gemm.build_identity))
            && gemm.tune_row_identity == gemm.validated_row.receipt_identity()
            && execution_path_is_complete(&gemm.execution_path)
            && execution_path_identity(&gemm.execution_path)
                .is_ok_and(|identity| gemm.execution_path_identity == identity)
            && gemm.validated_row.matches_decision(
                &gemm.scoped_tune_key,
                &gemm.shape_class,
                machine,
                &gemm.canonical_plan,
            )
    }

    fn declared_run(&self) -> u64 {
        self.gemm.execution_path.declared_run
    }

    fn stable_equivalent(&self, other: &Self) -> bool {
        let left = &self.gemm;
        let right = &other.gemm;
        left.scoped_tune_key == right.scoped_tune_key
            && left.shape_class == right.shape_class
            && left.canonical_plan == right.canonical_plan
            && left.source == right.source
            && left.operation_tier == right.operation_tier
            && left.build_identity == right.build_identity
            && left.tune_row_identity == right.tune_row_identity
            && left.validated_row == right.validated_row
            && execution_path_shape_eq(&left.execution_path, &right.execution_path)
    }

    fn canonical_json(&self) -> String {
        let gemm = &self.gemm;
        format!(
            "{{\"kind\":\"{EXECUTION_BINDING_KIND}\",\"scoped_tune_key\":\"{}\",\"shape_class\":\"{}\",\"plan\":\"{}\",\"source\":\"{}\",\"operation_tier\":\"{}\",\"build_identity\":\"{}\",\"tune_row_identity\":\"{}\",\"tune_row\":{},\"execution_path_identity\":\"{}\",\"execution_path\":{}}}",
            json_escape(&gemm.scoped_tune_key),
            json_escape(&gemm.shape_class),
            json_escape(&gemm.canonical_plan),
            gemm.source,
            json_escape(&gemm.operation_tier),
            json_escape(&gemm.build_identity),
            gemm.tune_row_identity,
            gemm.validated_row.receipt_json(),
            gemm.execution_path_identity,
            execution_path_json(&gemm.execution_path),
        )
    }

    fn receipt_identity(&self) -> fs_ledger::ContentHash {
        execution_binding_receipt_with_domain(&self.canonical_json(), EXECUTION_BINDING_DOMAIN)
    }

    fn validated_row(&self) -> &fs_session::ValidatedGemmTuneRow {
        &self.gemm.validated_row
    }
}

fn execution_path_shape_eq(
    left: &fs_session::GemmExecutionReceipt,
    right: &fs_session::GemmExecutionReceipt,
) -> bool {
    left.is_complete()
        && right.is_complete()
        && left.completed_tiles == right.completed_tiles
        && left.total_tiles == right.total_tiles
        && left.memory == right.memory
        && left.panels.len() == right.panels.len()
        && left.panels.iter().zip(&right.panels).all(|(left, right)| {
            left.kernel == right.kernel
                && left.mode == right.mode
                && left.completed == right.completed
                && left.total == right.total
        })
}

/// Which side of the roofline binds the limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RoofSide {
    /// Memory-bandwidth-bound at this intensity.
    Bandwidth,
    /// Compute-bound at this intensity.
    Compute,
}

/// Verdict against the kernel's declared band.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// Attainment ≥ target fraction.
    WithinBand,
    /// Attainment < target fraction.
    BelowBand,
    /// No target declared: report-only.
    NoTarget,
    /// The MEASUREMENT ENVIRONMENT is invalid: the probed axes fail
    /// absolute plausibility floors, or the kernel "beat" its roofline
    /// by an impossible margin (stale/contention-crushed axes). A gate
    /// must never pass — or fail — on a machine that was useless while
    /// it was measured (bead 1n61: a load-68 window collapsed both the
    /// STREAM probe and the kernel ~1000× together, and the RATIO
    /// self-normalized to a vacuous within_band).
    EnvironmentInvalid,
}

impl Verdict {
    /// Stable lowercase name for logs and ledger rows.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Verdict::WithinBand => "within_band",
            Verdict::BelowBand => "below_band",
            Verdict::NoTarget => "no_target",
            Verdict::EnvironmentInvalid => "environment_invalid",
        }
    }
}

/// One kernel's measured attainment against the machine roofline.
#[derive(Debug, Clone)]
pub struct Attainment {
    /// Kernel name.
    pub kernel: String,
    /// Kernel version.
    pub version: String,
    /// Median elements/second across repetitions.
    pub elems_per_sec: f64,
    /// Achieved memory traffic, GB/s.
    pub achieved_gbs: f64,
    /// Achieved compute, GFLOP/s.
    pub achieved_gflops: f64,
    /// Roofline limit in elements/second for this machine + intensity.
    pub limit_elems_per_sec: f64,
    /// Which axis binds.
    pub roof: RoofSide,
    /// `elems_per_sec / limit_elems_per_sec` (1.0 = at the roof).
    pub attainment: f64,
    /// Fraction of the spec's declared target axis. This equals `attainment`
    /// only when `target_axis == binding_roof` (or the declared physical axis
    /// happens to bind).
    pub target_attainment: f64,
    /// Relative interquartile dispersion of the repetition times
    /// ((p75 − p25) / median): a benchmark without variance bars is
    /// folklore.
    pub dispersion: f64,
    /// Repetitions measured (after warmup).
    pub reps: usize,
    /// Verdict against the declared band.
    pub verdict: Verdict,
    /// Why this row cannot support a verdict. Present exactly when
    /// `verdict == EnvironmentInvalid`.
    pub invalid_reason: Option<String>,
    axis_binding: AxisBinding,
    spec_binding: SpecBinding,
    measurement_origin: MeasurementOrigin,
    pending_tune_publication: Option<fs_session::ValidatedGemmTuneRow>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MeasurementOrigin {
    Analytic,
    Timed {
        elements: usize,
        warmup_runs: usize,
        sample_seconds_bits: Vec<u64>,
        decision_bindings: Vec<Option<KernelExecutionBinding>>,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AxisBinding {
    fingerprint: u64,
    logical_cpus: u32,
    bandwidth_single_bits: u64,
    bandwidth_all_core_bits: u64,
    peak_single_bits: u64,
    peak_all_core_bits: u64,
}

impl AxisBinding {
    fn new(axes: &MachineAxes) -> Self {
        Self {
            fingerprint: axes.fingerprint,
            logical_cpus: axes.logical_cpus,
            bandwidth_single_bits: axes.bandwidth_single_gbs.to_bits(),
            bandwidth_all_core_bits: axes.bandwidth_all_core_gbs.to_bits(),
            peak_single_bits: axes.peak_single_gflops.to_bits(),
            peak_all_core_bits: axes.peak_all_core_gflops.to_bits(),
        }
    }

    fn matches(&self, axes: &MachineAxes) -> bool {
        self == &Self::new(axes)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SpecBinding {
    kernel: String,
    version: String,
    bytes_per_elem_bits: u64,
    flops_per_elem_bits: u64,
    threading: Threading,
    target_axis: TargetAxis,
    target_fraction_bits: Option<u64>,
}

impl SpecBinding {
    fn new(spec: &KernelSpec) -> Self {
        Self {
            kernel: spec.name.to_string(),
            version: spec.version.to_string(),
            bytes_per_elem_bits: spec.bytes_per_elem.to_bits(),
            flops_per_elem_bits: spec.flops_per_elem.to_bits(),
            threading: spec.threading,
            target_axis: spec.target_axis,
            target_fraction_bits: spec.target_fraction.map(f64::to_bits),
        }
    }
}

fn stable_decision_binding(
    bindings: &[Option<KernelExecutionBinding>],
) -> Option<&KernelExecutionBinding> {
    let first = bindings.first()?.as_ref()?;
    bindings
        .iter()
        .enumerate()
        .all(|(ordinal, candidate)| {
            let Some(candidate) = candidate.as_ref() else {
                return false;
            };
            let Ok(ordinal) = u64::try_from(ordinal) else {
                return false;
            };
            let Some(expected_run) = first.declared_run().checked_add(ordinal) else {
                return false;
            };
            candidate.declared_run() == expected_run && candidate.stable_equivalent(first)
        })
        .then_some(first)
}

impl Attainment {
    /// One JSON line for logs/agents (stable field order).
    #[must_use]
    #[allow(clippy::too_many_lines)] // one canonical receipt schema, kept visibly ordered
    pub fn to_jsonl(&self) -> String {
        let invalid_reason = self.invalid_reason.as_ref().map_or_else(
            || "null".to_string(),
            |reason| format!("\"{}\"", json_escape(reason)),
        );
        let target_bits = self
            .spec_binding
            .target_fraction_bits
            .map_or_else(|| "null".to_string(), |bits| format!("\"{bits:016x}\""));
        let measurement = match &self.measurement_origin {
            MeasurementOrigin::Analytic => "{\"origin\":\"analytic\"}".to_string(),
            MeasurementOrigin::Timed {
                elements,
                warmup_runs,
                sample_seconds_bits,
                decision_bindings,
            } => {
                let sample = stats::sample_from_times(
                    sample_seconds_bits
                        .iter()
                        .copied()
                        .map(f64::from_bits)
                        .collect(),
                );
                let sample_bits = sample_seconds_bits
                    .iter()
                    .map(|bits| format!("\"{bits:016x}\""))
                    .collect::<Vec<_>>()
                    .join(",");
                let decision_bits = decision_bindings
                    .iter()
                    .map(|binding| {
                        binding.as_ref().map_or_else(
                            || "null".to_string(),
                            |binding| format!("\"{}\"", binding.receipt_identity()),
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                format!(
                    "{{\"origin\":\"timed\",\"elements\":{elements},\"warmup_runs\":{warmup_runs},\"sample_seconds_bits\":[{sample_bits}],\"decision_binding_hashes\":[{decision_bits}],\"median_seconds_bits\":\"{:016x}\",\"p25_seconds_bits\":\"{:016x}\",\"p75_seconds_bits\":\"{:016x}\",\"dispersion_bits\":\"{:016x}\"}}",
                    sample.median.to_bits(),
                    sample.p25.to_bits(),
                    sample.p75.to_bits(),
                    sample.dispersion.to_bits(),
                )
            }
        };
        let execution = match &self.measurement_origin {
            MeasurementOrigin::Timed {
                decision_bindings, ..
            } => stable_decision_binding(decision_bindings).map_or_else(
                || "null".to_string(),
                KernelExecutionBinding::canonical_json,
            ),
            MeasurementOrigin::Analytic => "null".to_string(),
        };
        format!(
            "{{\"receipt_version\":3,\"kernel\":\"{}\",\"version\":\"{}\",\"machine\":\"{:016x}\",\
             \"axes\":{{\"logical_cpus\":{},\"bandwidth_single_bits\":\"{:016x}\",\"bandwidth_all_core_bits\":\"{:016x}\",\"peak_single_bits\":\"{:016x}\",\"peak_all_core_bits\":\"{:016x}\"}},\
             \"spec\":{{\"bytes_per_elem_bits\":\"{:016x}\",\"flops_per_elem_bits\":\"{:016x}\",\"threading\":\"{}\",\"target_axis\":\"{}\",\"target_fraction_bits\":{}}},\"measurement\":{},\"execution\":{},\"elems_per_sec_bits\":\"{:016x}\",\"gbs_bits\":\"{:016x}\",\"gflops_bits\":\"{:016x}\",\"limit_elems_per_sec_bits\":\"{:016x}\",\"attainment_bits\":\"{:016x}\",\"target_attainment_bits\":\"{:016x}\",\"dispersion_bits\":\"{:016x}\",\"elems_per_sec_display\":{:.3e},\
             \"gbs\":{:.3},\"gflops\":{:.3},\"limit_elems_per_sec\":{:.3e},\
             \"roof\":\"{}\",\"attainment\":{:.4},\"target_attainment\":{:.4},\"dispersion\":{:.4},\
             \"reps\":{},\"verdict\":\"{}\",\"invalid_reason\":{}}}",
            json_escape(&self.kernel),
            json_escape(&self.version),
            self.axis_binding.fingerprint,
            self.axis_binding.logical_cpus,
            self.axis_binding.bandwidth_single_bits,
            self.axis_binding.bandwidth_all_core_bits,
            self.axis_binding.peak_single_bits,
            self.axis_binding.peak_all_core_bits,
            self.spec_binding.bytes_per_elem_bits,
            self.spec_binding.flops_per_elem_bits,
            match self.spec_binding.threading {
                Threading::SingleThread => "single_thread",
                Threading::AllCore => "all_core",
            },
            self.spec_binding.target_axis.name(),
            target_bits,
            measurement,
            execution,
            self.elems_per_sec.to_bits(),
            self.achieved_gbs.to_bits(),
            self.achieved_gflops.to_bits(),
            self.limit_elems_per_sec.to_bits(),
            self.attainment.to_bits(),
            self.target_attainment.to_bits(),
            self.dispersion.to_bits(),
            self.elems_per_sec,
            self.achieved_gbs,
            self.achieved_gflops,
            self.limit_elems_per_sec,
            match self.roof {
                RoofSide::Bandwidth => "bandwidth",
                RoofSide::Compute => "compute",
            },
            self.attainment,
            self.target_attainment,
            self.dispersion,
            self.reps,
            self.verdict.name(),
            invalid_reason,
        )
    }

    #[allow(clippy::too_many_lines)] // fail-closed rederivation is easier to audit as one flow
    fn is_citable_against(&self, axes: &MachineAxes) -> bool {
        if !self.axis_binding.matches(axes)
            || self.kernel != self.spec_binding.kernel
            || self.version != self.spec_binding.version
            || self.kernel.trim().is_empty()
            || self.version.trim().is_empty()
            || self.reps == 0
            || self.invalid_reason.is_some()
            || self.verdict == Verdict::EnvironmentInvalid
        {
            return false;
        }
        let MeasurementOrigin::Timed {
            elements,
            warmup_runs,
            sample_seconds_bits,
            decision_bindings,
        } = &self.measurement_origin
        else {
            return false;
        };
        let sample_times: Vec<_> = sample_seconds_bits
            .iter()
            .copied()
            .map(f64::from_bits)
            .collect();
        if *elements == 0
            || sample_times.is_empty()
            || sample_times.len() != self.reps
            || decision_bindings.len() != self.reps
            || sample_times
                .iter()
                .any(|seconds| !seconds.is_finite() || *seconds <= 0.0)
        {
            return false;
        }
        let sample = stats::sample_from_times(sample_times);
        if ((*elements as f64) / sample.median).to_bits() != self.elems_per_sec.to_bits()
            || sample.dispersion.to_bits() != self.dispersion.to_bits()
        {
            return false;
        }
        let bytes_per_elem = f64::from_bits(self.spec_binding.bytes_per_elem_bits);
        let flops_per_elem = f64::from_bits(self.spec_binding.flops_per_elem_bits);
        let target = self.spec_binding.target_fraction_bits.map(f64::from_bits);
        if !bytes_per_elem.is_finite()
            || !flops_per_elem.is_finite()
            || bytes_per_elem < 0.0
            || flops_per_elem < 0.0
            || (bytes_per_elem == 0.0 && flops_per_elem == 0.0)
            || target.is_some_and(|value| !value.is_finite() || value <= 0.0 || value > 1.0)
            || !self.elems_per_sec.is_finite()
            || self.elems_per_sec < 0.0
            || !self.dispersion.is_finite()
            || self.dispersion < 0.0
        {
            return false;
        }
        let (bandwidth_gbs, peak_gflops) = match self.spec_binding.threading {
            Threading::SingleThread => (axes.bandwidth_single_gbs, axes.peak_single_gflops),
            Threading::AllCore => (axes.bandwidth_all_core_gbs, axes.peak_all_core_gflops),
        };
        let bandwidth_limit = if bytes_per_elem > 0.0 {
            bandwidth_gbs * 1e9 / bytes_per_elem
        } else {
            f64::INFINITY
        };
        let compute_limit = if flops_per_elem > 0.0 {
            peak_gflops * 1e9 / flops_per_elem
        } else {
            f64::INFINITY
        };
        let (limit, roof) = if bandwidth_limit <= compute_limit {
            (bandwidth_limit, RoofSide::Bandwidth)
        } else {
            (compute_limit, RoofSide::Compute)
        };
        if !limit.is_finite() || limit <= 0.0 {
            return false;
        }
        let attainment = self.elems_per_sec / limit;
        let achieved_gbs = self.elems_per_sec * bytes_per_elem / 1e9;
        let achieved_gflops = self.elems_per_sec * flops_per_elem / 1e9;
        let target_attainment = match self.spec_binding.target_axis {
            TargetAxis::BindingRoof => attainment,
            TargetAxis::ComputePeak => achieved_gflops / peak_gflops,
            TargetAxis::MemoryBandwidth => achieved_gbs / bandwidth_gbs,
        };
        if !attainment.is_finite()
            || attainment > 1.5
            || !target_attainment.is_finite()
            || !achieved_gbs.is_finite()
            || !achieved_gflops.is_finite()
        {
            return false;
        }
        let expected_verdict = match target {
            None => Verdict::NoTarget,
            Some(fraction) if target_attainment >= fraction => Verdict::WithinBand,
            Some(_) => Verdict::BelowBand,
        };
        let decision_valid = if self.kernel == "gemm-f64" {
            *warmup_runs > 0
                && stable_decision_binding(decision_bindings)
                    .is_some_and(|binding| binding.is_valid_for(axes.fingerprint))
        } else {
            decision_bindings.iter().all(Option::is_none)
        };
        let pending_valid = self
            .pending_tune_publication
            .as_ref()
            .is_none_or(|pending| {
                self.kernel == "gemm-f64"
                    && stable_decision_binding(decision_bindings).is_some_and(|binding| {
                        pending.receipt_identity() == binding.validated_row().receipt_identity()
                    })
            });
        self.roof == roof
            && self.limit_elems_per_sec.to_bits() == limit.to_bits()
            && self.attainment.to_bits() == attainment.to_bits()
            && self.target_attainment.to_bits() == target_attainment.to_bits()
            && self.achieved_gbs.to_bits() == achieved_gbs.to_bits()
            && self.achieved_gflops.to_bits() == achieved_gflops.to_bits()
            && self.verdict == expected_verdict
            && decision_valid
            && pending_valid
    }
}

fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                use std::fmt::Write as _;
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out
}

fn spec_error(spec: &KernelSpec) -> Option<&'static str> {
    if spec.name.trim().is_empty() || spec.version.trim().is_empty() {
        return Some("kernel name and version must be non-empty");
    }
    if !spec.bytes_per_elem.is_finite()
        || !spec.flops_per_elem.is_finite()
        || spec.bytes_per_elem < 0.0
        || spec.flops_per_elem < 0.0
        || (spec.bytes_per_elem == 0.0 && spec.flops_per_elem == 0.0)
    {
        return Some("kernel intensity must be finite, non-negative, and exercise an axis");
    }
    if let Some(target) = spec.target_fraction
        && (!target.is_finite() || target <= 0.0 || target > 1.0)
    {
        return Some("target fraction must be finite and in (0, 1]");
    }
    if spec.target_fraction.is_some()
        && match spec.target_axis {
            TargetAxis::BindingRoof => false,
            TargetAxis::ComputePeak => spec.flops_per_elem == 0.0,
            TargetAxis::MemoryBandwidth => spec.bytes_per_elem == 0.0,
        }
    {
        return Some("declared target axis is not exercised by this kernel");
    }
    None
}

/// Compute attainment from a measured rate and the machine axes. Pure
/// arithmetic — meta-tested against hand calculations (`rf_002`).
#[must_use]
pub fn attainment_for(spec: &KernelSpec, elems_per_sec: f64, axes: &MachineAxes) -> Attainment {
    attainment_with_dispersion(spec, elems_per_sec, 0.0, 0, axes)
}

/// [`attainment_for`] with measured dispersion and repetition count.
#[must_use]
pub fn attainment_with_dispersion(
    spec: &KernelSpec,
    elems_per_sec: f64,
    dispersion: f64,
    reps: usize,
    axes: &MachineAxes,
) -> Attainment {
    attainment_with_origin(
        spec,
        elems_per_sec,
        dispersion,
        reps,
        axes,
        MeasurementOrigin::Analytic,
    )
}

#[allow(clippy::too_many_lines)] // one roof decision per branch; splitting obscures the min-roof flow
fn attainment_with_origin(
    spec: &KernelSpec,
    elems_per_sec: f64,
    dispersion: f64,
    reps: usize,
    axes: &MachineAxes,
    measurement_origin: MeasurementOrigin,
) -> Attainment {
    let spec_error = spec_error(spec);
    let measurement_error = if !elems_per_sec.is_finite() || elems_per_sec < 0.0 {
        Some("measured element rate is non-finite or negative")
    } else if !dispersion.is_finite() || dispersion < 0.0 {
        Some("measured dispersion is non-finite or negative")
    } else {
        None
    };
    let safe_rate = if measurement_error.is_none() {
        elems_per_sec
    } else {
        0.0
    };
    let safe_bytes = if spec.bytes_per_elem.is_finite() && spec.bytes_per_elem >= 0.0 {
        spec.bytes_per_elem
    } else {
        0.0
    };
    let safe_flops = if spec.flops_per_elem.is_finite() && spec.flops_per_elem >= 0.0 {
        spec.flops_per_elem
    } else {
        0.0
    };
    let (bandwidth_gbs, peak_gflops) = match spec.threading {
        Threading::SingleThread => (axes.bandwidth_single_gbs, axes.peak_single_gflops),
        Threading::AllCore => (axes.bandwidth_all_core_gbs, axes.peak_all_core_gflops),
    };
    // Limits in elements/second on each axis; +inf when the kernel does not
    // exercise an axis (zero bytes or zero flops per element).
    let bw_limit = if safe_bytes > 0.0 && bandwidth_gbs.is_finite() && bandwidth_gbs > 0.0 {
        bandwidth_gbs * 1e9 / safe_bytes
    } else {
        f64::INFINITY
    };
    let comp_limit = if safe_flops > 0.0 && peak_gflops.is_finite() && peak_gflops > 0.0 {
        peak_gflops * 1e9 / safe_flops
    } else {
        f64::INFINITY
    };
    let (limit, roof) = if bw_limit <= comp_limit {
        (bw_limit, RoofSide::Bandwidth)
    } else {
        (comp_limit, RoofSide::Compute)
    };
    let raw_attainment = if limit.is_finite() && limit > 0.0 {
        safe_rate / limit
    } else {
        0.0
    };
    let raw_achieved_gbs = safe_rate * safe_bytes / 1e9;
    let raw_achieved_gflops = safe_rate * safe_flops / 1e9;
    let raw_target_attainment = match spec.target_axis {
        TargetAxis::BindingRoof => raw_attainment,
        TargetAxis::ComputePeak if peak_gflops.is_finite() && peak_gflops > 0.0 => {
            raw_achieved_gflops / peak_gflops
        }
        TargetAxis::MemoryBandwidth if bandwidth_gbs.is_finite() && bandwidth_gbs > 0.0 => {
            raw_achieved_gbs / bandwidth_gbs
        }
        TargetAxis::ComputePeak | TargetAxis::MemoryBandwidth => f64::NAN,
    };
    let derived_error = if !raw_attainment.is_finite()
        || !raw_target_attainment.is_finite()
        || !raw_achieved_gbs.is_finite()
        || !raw_achieved_gflops.is_finite()
    {
        Some("derived roofline quantities overflowed or became non-finite")
    } else {
        None
    };
    let attainment = if raw_attainment.is_finite() {
        raw_attainment
    } else {
        0.0
    };
    let target_attainment = if raw_target_attainment.is_finite() {
        raw_target_attainment
    } else {
        0.0
    };
    // Environment validity BEFORE band comparison (bead 1n61): the
    // ratio is meaningless when the axes are implausible, and an
    // attainment materially above 1 means the kernel outran its own
    // roofline — the axes were probed under different (crushed)
    // conditions than the kernel run. Refuse to gate either way.
    let invalid_reason = axes
        .plausibility_error()
        .or(spec_error)
        .or(measurement_error)
        .or(derived_error)
        .or_else(|| (attainment > 1.5).then_some("attainment exceeds the credible roofline band"));
    let verdict = if invalid_reason.is_some() {
        Verdict::EnvironmentInvalid
    } else {
        match spec.target_fraction {
            None => Verdict::NoTarget,
            Some(t) if target_attainment >= t => Verdict::WithinBand,
            Some(_) => Verdict::BelowBand,
        }
    };
    Attainment {
        kernel: spec.name.to_string(),
        version: spec.version.to_string(),
        elems_per_sec: safe_rate,
        achieved_gbs: if raw_achieved_gbs.is_finite() {
            raw_achieved_gbs
        } else {
            0.0
        },
        achieved_gflops: if raw_achieved_gflops.is_finite() {
            raw_achieved_gflops
        } else {
            0.0
        },
        limit_elems_per_sec: if limit.is_finite() { limit } else { 0.0 },
        roof,
        attainment,
        target_attainment,
        dispersion: if dispersion.is_finite() && dispersion >= 0.0 {
            dispersion
        } else {
            0.0
        },
        reps,
        verdict,
        invalid_reason: invalid_reason.map(str::to_string),
        axis_binding: AxisBinding::new(axes),
        spec_binding: SpecBinding::new(spec),
        measurement_origin,
        pending_tune_publication: None,
    }
}

/// Largest untimed warmup count accepted by any public roofline measurement.
pub const MAX_MEASUREMENT_WARMUP: usize = 1_000;
/// Largest timed repetition count accepted by any public roofline measurement.
pub const MAX_MEASUREMENT_REPS: usize = 1_000;
/// Largest kernel registry accepted by the public harness.
pub const MAX_REGISTRY_KERNELS: usize = 256;
/// Largest aggregate kernel invocation count accepted by one registry run.
pub const MAX_REGISTRY_KERNEL_INVOCATIONS: usize = 250_000;

fn validate_measurement_work(warmup: usize, reps: usize, kernels: usize) -> Result<usize, String> {
    if warmup > MAX_MEASUREMENT_WARMUP {
        return Err(format!(
            "roofline warmup must be in 0..={MAX_MEASUREMENT_WARMUP}, got {warmup}"
        ));
    }
    if reps == 0 || reps > MAX_MEASUREMENT_REPS {
        return Err(format!(
            "roofline repetitions must be in 1..={MAX_MEASUREMENT_REPS}, got {reps}"
        ));
    }
    if kernels == 0 || kernels > MAX_REGISTRY_KERNELS {
        return Err(format!(
            "roofline registry must contain 1..={MAX_REGISTRY_KERNELS} kernels, got {kernels}"
        ));
    }
    let runs_per_kernel = warmup
        .checked_add(reps)
        .ok_or_else(|| "roofline warmup + repetition count overflowed usize".to_string())?;
    let invocations = runs_per_kernel
        .checked_mul(kernels)
        .ok_or_else(|| "roofline aggregate kernel invocation count overflowed usize".to_string())?;
    if invocations > MAX_REGISTRY_KERNEL_INVOCATIONS {
        return Err(format!(
            "roofline registry requires {invocations} kernel invocations, exceeding the {MAX_REGISTRY_KERNEL_INVOCATIONS}-invocation bound"
        ));
    }
    Ok(runs_per_kernel)
}

/// Measure one kernel (warmup + repetitions) and compute its attainment.
///
/// # Errors
/// Refuses zero or hostile repetition counts before allocation or execution,
/// and reports bounded sample-buffer reservation failures.
pub fn measure(
    kernel: &mut dyn RooflineKernel,
    warmup: usize,
    reps: usize,
    axes: &MachineAxes,
) -> Result<Attainment, String> {
    validate_measurement_work(warmup, reps, 1)?;
    let spec = kernel.spec();
    let elements = kernel.elements();
    let elems = elements as f64;
    let mut times = Vec::new();
    times
        .try_reserve_exact(reps)
        .map_err(|_| format!("cannot reserve {reps} roofline timing samples"))?;
    let mut decision_bindings = Vec::new();
    decision_bindings
        .try_reserve_exact(reps)
        .map_err(|_| format!("cannot reserve {reps} roofline decision bindings"))?;
    for invocation in 0..warmup {
        kernel.run_once().map_err(|error| {
            format!(
                "roofline kernel `{}` failed during warmup invocation {invocation}: {error}",
                spec.name
            )
        })?;
    }
    for invocation in 0..reps {
        let start = std::time::Instant::now();
        kernel.run_once().map_err(|error| {
            format!(
                "roofline kernel `{}` failed during timed invocation {invocation}: {error}",
                spec.name
            )
        })?;
        times.push(start.elapsed().as_secs_f64());
        decision_bindings.push(kernel.execution_binding());
    }
    let sample = stats::sample_from_times(times);
    let elems_per_sec = if sample.median > 0.0 {
        elems / sample.median
    } else {
        0.0
    };
    let mut result = attainment_with_origin(
        &spec,
        elems_per_sec,
        sample.dispersion,
        reps,
        axes,
        MeasurementOrigin::Timed {
            elements,
            warmup_runs: warmup,
            sample_seconds_bits: sample.times.iter().map(|value| value.to_bits()).collect(),
            decision_bindings,
        },
    );
    result.pending_tune_publication = kernel.pending_tune_publication();
    Ok(result)
}

/// Run every kernel in the registry after admitting the complete work shape.
///
/// # Errors
/// Refuses empty/oversized registries and hostile aggregate repetition counts
/// before any kernel executes, and reports bounded result/sample allocation
/// failures. If any kernel fails, every registry member receives an idempotent
/// tuning abort before the error is returned.
pub fn run_registry(
    registry: &mut [Box<dyn RooflineKernel>],
    warmup: usize,
    reps: usize,
    axes: &MachineAxes,
) -> Result<Vec<Attainment>, String> {
    validate_measurement_work(warmup, reps, registry.len())?;
    let mut results = Vec::new();
    results
        .try_reserve_exact(registry.len())
        .map_err(|_| format!("cannot reserve {} roofline results", registry.len()))?;
    let mut measurement_error = None;
    for kernel in registry.iter_mut() {
        match measure(kernel.as_mut(), warmup, reps, axes) {
            Ok(result) => results.push(result),
            Err(error) => {
                measurement_error = Some(error);
                break;
            }
        }
    }
    if let Some(error) = measurement_error {
        let abort_diagnostics = abort_registry_tuning(registry);
        if abort_diagnostics.is_empty() {
            return Err(error);
        }
        return Err(format!(
            "{error}; registry tuning abort also failed: {}",
            abort_diagnostics.join("; ")
        ));
    }
    poison_invalid_run(&mut results);
    Ok(results)
}

fn abort_registry_tuning(registry: &mut [Box<dyn RooflineKernel>]) -> Vec<String> {
    let mut diagnostics = Vec::new();
    for (index, kernel) in registry.iter_mut().enumerate() {
        if let Err(error) = kernel.abort_tuning() {
            diagnostics.push(format!("kernel[{index}]: {error}"));
        }
    }
    diagnostics
}

/// Whether a registry result set passes aggregate measurement admission for
/// these exact measured axes. This checks timing, binding, baseline, and result
/// integrity only. It does not establish production-registry provenance and is
/// therefore not, by itself, citation authority; only
/// [`production::ProductionRun::citation_eligible`] combines both conditions.
#[must_use]
pub fn run_passes_measurement_admission(
    axes: &MachineAxes,
    post_axes: &MachineAxes,
    baseline: AxisBaselinePolicy<'_>,
    results: &[Attainment],
) -> bool {
    let snapshot = baseline.snapshot(axes, post_axes);
    run_admission_error_for_snapshot(axes, &snapshot, results).is_none()
}

fn push_receipt_field(payload: &mut Vec<u8>, bytes: &[u8]) {
    let len = u64::try_from(bytes.len()).expect("receipt field length fits u64");
    payload.extend_from_slice(&len.to_le_bytes());
    payload.extend_from_slice(bytes);
}

/// One row of the op-bound ordered result manifest: which kernel/version sits
/// at which position of the finalized result set, and the content hash of its
/// exact payload bytes.
struct ManifestEntry {
    ordinal: u64,
    kernel: String,
    version: String,
    payload: fs_blake3::ContentHash,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FinalizedRunIdentityInput {
    baseline_receipt: String,
    result_payloads: Vec<String>,
    result_manifest: String,
}

impl FinalizedRunIdentityInput {
    fn from_run(admission: &AxisAdmissionSnapshot, results: &[Attainment]) -> Self {
        Self {
            baseline_receipt: admission.receipt_json().to_string(),
            result_payloads: results.iter().map(Attainment::to_jsonl).collect(),
            result_manifest: run_result_manifest_json(results),
        }
    }
}

#[allow(dead_code)]
fn classify_finalized_run_identity_fields(input: &FinalizedRunIdentityInput) {
    let FinalizedRunIdentityInput {
        baseline_receipt,
        result_payloads,
        result_manifest,
    } = input;
    let _ = (baseline_receipt, result_payloads, result_manifest);
}

fn manifest_entry_json(ordinal: u64, kernel: &str, version: &str, payload: &str) -> String {
    format!(
        "{{\"ordinal\":{ordinal},\"kernel\":\"{}\",\"version\":\"{}\",\"payload\":\"{payload}\"}}",
        json_escape(kernel),
        json_escape(version),
    )
}

fn valid_manifest_identifier(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 128
        && value.bytes().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
        })
}

/// Canonical JSON for the ordered result manifest. Kernel names and versions
/// are static registry identifiers; names needing JSON escaping are refused at
/// parse time, so serialization never escapes.
fn run_result_manifest_json(results: &[Attainment]) -> String {
    let entries = results
        .iter()
        .enumerate()
        .map(|(ordinal, result)| {
            manifest_entry_json(
                ordinal as u64,
                &result.kernel,
                &result.version,
                &fs_ledger::hash_bytes(result.to_jsonl().as_bytes()).to_string(),
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("{{\"schema\":\"{RESULT_MANIFEST_SCHEMA}\",\"entries\":[{entries}]}}")
}

/// Strict parse of a run-result manifest. Ordinals must be exactly
/// `0..entries.len()` in order; kernel/version must be plain identifiers
/// (no quotes, backslashes, or control bytes — escaping would break the
/// byte-exact round-trip this parser enforces).
fn parse_result_manifest(text: &str) -> Option<Vec<ManifestEntry>> {
    let body = text
        .strip_prefix(&format!(
            "{{\"schema\":\"{RESULT_MANIFEST_SCHEMA}\",\"entries\":["
        ))?
        .strip_suffix("]}")?;
    let mut entries = Vec::new();
    if !body.is_empty() {
        for raw in body.split("},") {
            let raw_entry = if raw.ends_with('}') {
                raw.to_string()
            } else {
                format!("{raw}}}")
            };
            let inner = raw_entry
                .strip_prefix("{\"ordinal\":")?
                .strip_suffix("\"}")?;
            let (ordinal_text, rest) = inner.split_once(",\"kernel\":\"")?;
            let (kernel, rest) = rest.split_once("\",\"version\":\"")?;
            let (version, payload_hex) = rest.split_once("\",\"payload\":\"")?;
            if !valid_manifest_identifier(kernel) || !valid_manifest_identifier(version) {
                return None;
            }
            let ordinal: u64 = ordinal_text.parse().ok()?;
            if ordinal != u64::try_from(entries.len()).ok()? {
                return None;
            }
            let payload = fs_blake3::ContentHash::from_hex(payload_hex)?;
            entries.push(ManifestEntry {
                ordinal,
                kernel: kernel.to_string(),
                version: version.to_string(),
                payload,
            });
        }
    }
    // Byte-exact round trip: any formatting the serializer would not produce
    // (whitespace, reordered fields, uppercase hex) is refused.
    let reserialized = format!(
        "{{\"schema\":\"{RESULT_MANIFEST_SCHEMA}\",\"entries\":[{}]}}",
        entries
            .iter()
            .map(|e| manifest_entry_json(e.ordinal, &e.kernel, &e.version, &e.payload.to_string()))
            .collect::<Vec<_>>()
            .join(",")
    );
    (reserialized == text).then_some(entries)
}

fn finalized_run_receipt_from_input(
    input: &FinalizedRunIdentityInput,
    finalized_domain: &str,
    result_manifest_domain: &str,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    push_receipt_field(&mut payload, input.baseline_receipt.as_bytes());
    let result_count = u64::try_from(input.result_payloads.len()).expect("result count fits u64");
    payload.extend_from_slice(&result_count.to_le_bytes());
    for result_payload in &input.result_payloads {
        push_receipt_field(&mut payload, result_payload.as_bytes());
    }
    let manifest_hash =
        fs_blake3::hash_domain(result_manifest_domain, input.result_manifest.as_bytes());
    push_receipt_field(&mut payload, manifest_hash.as_bytes());
    fs_blake3::hash_domain(finalized_domain, &payload)
}

fn finalized_run_receipt(
    admission: &AxisAdmissionSnapshot,
    results: &[Attainment],
) -> fs_blake3::ContentHash {
    let input = FinalizedRunIdentityInput::from_run(admission, results);
    finalized_run_receipt_from_input(&input, FINALIZED_RUN_DOMAIN, RESULT_MANIFEST_DOMAIN)
}

/// Complete the registry's tuning lifecycle at the same admission boundary
/// used for citable performance evidence.
///
/// Every kernel sees the single aggregate decision. Publication belongs to
/// [`record_run`]'s evidence transaction; this hook clears kernel-owned pending
/// markers and invalidates rejected local decisions.
///
/// # Errors
/// Calls every kernel, then returns all lifecycle diagnostics in deterministic
/// registry order. Any diagnostic triggers an idempotent tuning abort across
/// the complete registry before this function returns without a token.
pub fn finalize_registry_tuning(
    registry: &mut [Box<dyn RooflineKernel>],
    axes: &MachineAxes,
    post_axes: &MachineAxes,
    baseline: AxisBaselinePolicy<'_>,
    results: &[Attainment],
) -> Result<FinalizedRegistryRun, String> {
    let admission = baseline.snapshot(axes, post_axes);
    finalize_registry_tuning_with_snapshot(registry, axes, &admission, results)
}

pub(crate) fn finalize_registry_tuning_with_snapshot(
    registry: &mut [Box<dyn RooflineKernel>],
    axes: &MachineAxes,
    admission: &AxisAdmissionSnapshot,
    results: &[Attainment],
) -> Result<FinalizedRegistryRun, String> {
    let mut diagnostics = Vec::new();
    let mut registry_matches = registry.len() == results.len();
    if !registry_matches {
        diagnostics.push(format!(
            "registry/result length mismatch: {} kernels for {} result rows",
            registry.len(),
            results.len()
        ));
    }
    for (index, (kernel, result)) in registry.iter().zip(results).enumerate() {
        let spec = kernel.spec();
        if spec.name != result.kernel || spec.version != result.version {
            registry_matches = false;
            diagnostics.push(format!(
                "kernel[{index}] identity mismatch: registry {}/{} vs result {}/{}",
                spec.name, spec.version, result.kernel, result.version
            ));
        }
        let expected_binding = match &result.measurement_origin {
            MeasurementOrigin::Timed {
                decision_bindings, ..
            } => decision_bindings.last().cloned().flatten(),
            MeasurementOrigin::Analytic => None,
        };
        if kernel.execution_binding() != expected_binding {
            registry_matches = false;
            diagnostics.push(format!(
                "kernel[{index}] execution state changed after this result was measured"
            ));
        }
        if kernel.pending_tune_publication() != result.pending_tune_publication {
            registry_matches = false;
            diagnostics.push(format!(
                "kernel[{index}] pending tune publication changed after this result was measured"
            ));
        }
    }
    let admitted =
        registry_matches && run_admission_error_for_snapshot(axes, admission, results).is_none();
    for (index, kernel) in registry.iter_mut().enumerate() {
        if let Err(error) = kernel.finalize_tuning(admitted) {
            diagnostics.push(format!("kernel[{index}]: {error}"));
        }
    }
    if !diagnostics.is_empty() {
        let abort_diagnostics = abort_registry_tuning(registry);
        diagnostics.extend(
            abort_diagnostics
                .into_iter()
                .map(|diagnostic| format!("abort {diagnostic}")),
        );
        return Err(format!(
            "tuning lifecycle finalization failed with {} issue(s): {}",
            diagnostics.len(),
            diagnostics.join("; ")
        ));
    }
    Ok(FinalizedRegistryRun {
        receipt: finalized_run_receipt(admission, results),
        admitted,
        consumed: false,
    })
}

/// Explain why a result set cannot be admitted as citable performance
/// evidence. `None` means every row is a bound timed measurement.
#[must_use]
pub fn run_admission_error(
    axes: &MachineAxes,
    post_axes: &MachineAxes,
    baseline: AxisBaselinePolicy<'_>,
    results: &[Attainment],
) -> Option<String> {
    let snapshot = baseline.snapshot(axes, post_axes);
    run_admission_error_for_snapshot(axes, &snapshot, results)
}

pub(crate) fn run_admission_error_for_snapshot(
    axes: &MachineAxes,
    admission: &AxisAdmissionSnapshot,
    results: &[Attainment],
) -> Option<String> {
    let live_day_error = admission.live_day_error();
    run_admission_error_for_snapshot_with_live_error(axes, admission, results, live_day_error)
}

pub(crate) fn citable_run_admission_error_for_snapshot(
    axes: &MachineAxes,
    admission: &AxisAdmissionSnapshot,
    results: &[Attainment],
) -> Option<String> {
    let live_day_error = admission.live_day_error();
    run_admission_error_for_snapshot_with_live_error(
        axes,
        admission,
        results,
        live_day_error.clone(),
    )
    .or_else(|| admission.baseline_citation_error_with_live_day(live_day_error))
}

fn run_admission_error_for_snapshot_with_live_error(
    axes: &MachineAxes,
    admission: &AxisAdmissionSnapshot,
    results: &[Attainment],
    live_day_error: Option<String>,
) -> Option<String> {
    if results.is_empty() {
        return Some("registry produced no measured kernels".to_string());
    }
    if let Some(error) = live_day_error {
        return Some(error);
    }
    let baseline_verdict = admission.verdict();
    if !baseline_verdict.trusted() {
        return Some(format!(
            "historical baseline admission refused: {}",
            baseline_verdict.to_jsonl()
        ));
    }
    let mut identities = std::collections::BTreeSet::new();
    for (index, result) in results.iter().enumerate() {
        if !valid_manifest_identifier(&result.kernel) || !valid_manifest_identifier(&result.version)
        {
            return Some(format!(
                "row {index} has a non-canonical kernel/version identifier; expected 1..=128 ASCII alphanumeric, '-', '_', '.', ':', or '/' bytes"
            ));
        }
        if !identities.insert(result.kernel.as_str()) {
            return Some(format!(
                "duplicate kernel identity at row {index}: {}/{}",
                result.kernel, result.version
            ));
        }
        if !result.is_citable_against(axes) {
            return Some(format!(
                "row {index} ({}/{}) is not a bound timed measurement for these axes",
                result.kernel, result.version
            ));
        }
    }
    None
}

fn poison_invalid_run(results: &mut [Attainment]) {
    let Some(origin) = results
        .iter()
        .find(|result| result.verdict == Verdict::EnvironmentInvalid)
    else {
        return;
    };
    let reason = format!(
        "registry invalidated by {}: {}",
        origin.kernel,
        origin
            .invalid_reason
            .as_deref()
            .unwrap_or("invalid evidence")
    );
    for result in results {
        result.verdict = Verdict::EnvironmentInvalid;
        result.invalid_reason = Some(reason.clone());
    }
}

// ---------------------------------------------------------------------------
// §14.1 target table as data
// ---------------------------------------------------------------------------

/// One row of the plan §14.1 target table. `landed = false` rows are
/// visible from day one so nothing is silently uncovered — they flip as the
/// owning kernels register.
#[derive(Debug, Clone, Copy)]
pub struct TargetRow {
    /// Kernel family name.
    pub kernel: &'static str,
    /// What the target means (unit and roof context).
    pub statement: &'static str,
    /// Whether an implementation is registered in this harness yet.
    pub landed: bool,
}

/// The §14.1 table. Statements, not claims: every value must be re-measured
/// on a fingerprinted machine before anyone may cite it.
pub const SECTION_14_1_TARGETS: &[TargetRow] = &[
    TargetRow {
        kernel: "lbm-d3q19-stream-collide",
        statement: "≥1.0 GLUP/s (M-class) / ≥0.6 GLUP/s (TR-class), bandwidth-bound",
        landed: false,
    },
    TargetRow {
        kernel: "gemm-f64",
        statement: "≥75% of measured peak FLOPs for the selected SIMD tier",
        landed: true,
    },
    TargetRow {
        kernel: "spmv-sell-c-sigma",
        statement: "≥85% of measured STREAM-class bandwidth",
        landed: false,
    },
    TargetRow {
        kernel: "feec-apply-p4",
        statement: "≥30% of peak FLOPs, sum-factorized",
        landed: false,
    },
    TargetRow {
        kernel: "batched-small-dense",
        statement: "≥60% of peak FLOPs, SIMD-across-elements",
        landed: false,
    },
    TargetRow {
        kernel: "fft-3d-pencil",
        statement: "≥40% of the memory-bound limit",
        landed: false,
    },
    TargetRow {
        kernel: "sdf-primary-rays",
        statement: "≥80 Mray/s (M-class) / ≥120 Mray/s (TR-class)",
        landed: false,
    },
];

// v4 adds the exact dependency-receipt artifact and domain digest. Earlier
// rows cannot prove which resolved normal/build graph produced the GEMM code,
// so roofline-v8 does not reuse their shape keys. The v8 namespace also
// retires production-v2 rows: their admission proof is not valid under the
// attested production-v3 protocol, and append-only v7 history must not poison
// later v3 evidence.
const ROOFLINE_ROW_SCHEMA: &str = "fs-roofline-ledger-row-v4";
const ROOFLINE_PAYLOAD_ARTIFACT_KIND: &str = "roofline-benchmark-result";
const ROOFLINE_PAYLOAD_ARTIFACT_META: &str = "{\"schema\":\"fs-roofline-benchmark-result-v1\"}";
/// Maximum age of a roofline measurement that can be reported as fresh.
pub const STALENESS_MAX_AGE_NS: i64 = 30 * 24 * 60 * 60 * 1_000_000_000;

fn versions_json(build_identity: fs_blake3::ContentHash) -> String {
    format!("{{\"frankensim_executable\":\"{build_identity}\",\"fs-roofline\":\"{VERSION}\"}}")
}

fn executable_build_identity() -> Result<fs_blake3::ContentHash, LedgerError> {
    static IDENTITY: std::sync::OnceLock<Result<fs_blake3::ContentHash, LedgerError>> =
        std::sync::OnceLock::new();
    IDENTITY.get_or_init(read_executable_build_identity).clone()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ExecutableBuildIdentityInput {
    byte_len: u64,
    raw_hash: fs_blake3::ContentHash,
}

#[allow(dead_code)]
fn classify_executable_build_identity_fields(input: &ExecutableBuildIdentityInput) {
    let ExecutableBuildIdentityInput { byte_len, raw_hash } = input;
    let _ = (byte_len, raw_hash);
}

fn executable_build_identity_from_input(
    input: &ExecutableBuildIdentityInput,
    domain: &str,
) -> fs_blake3::ContentHash {
    let mut preimage = [0_u8; 40];
    preimage[..8].copy_from_slice(&input.byte_len.to_le_bytes());
    preimage[8..].copy_from_slice(input.raw_hash.as_bytes());
    fs_blake3::hash_domain(domain, &preimage)
}

fn read_executable_build_identity() -> Result<fs_blake3::ContentHash, LedgerError> {
    use std::io::Read as _;

    let path = std::env::current_exe().map_err(|error| LedgerError::Invalid {
        field: "executable_identity".to_string(),
        problem: format!("cannot resolve current executable: {error}"),
    })?;
    let mut file = std::fs::File::open(&path).map_err(|error| LedgerError::Invalid {
        field: "executable_identity".to_string(),
        problem: format!("cannot open current executable {}: {error}", path.display()),
    })?;
    let mut hasher = fs_blake3::Blake3::new();
    let mut total = 0_u64;
    let mut chunk = Vec::new();
    chunk
        .try_reserve_exact(64 * 1024)
        .map_err(|error| LedgerError::Invalid {
            field: "executable_identity".to_string(),
            problem: format!("cannot reserve executable hash buffer: {error}"),
        })?;
    chunk.resize(64 * 1024, 0_u8);
    loop {
        let read = file
            .read(&mut chunk)
            .map_err(|error| LedgerError::Invalid {
                field: "executable_identity".to_string(),
                problem: format!("cannot read current executable {}: {error}", path.display()),
            })?;
        if read == 0 {
            break;
        }
        total = total
            .checked_add(u64::try_from(read).expect("read chunk length fits u64"))
            .ok_or_else(|| LedgerError::Invalid {
                field: "executable_identity".to_string(),
                problem: "current executable length exceeds u64".to_string(),
            })?;
        hasher.update(&chunk[..read]);
    }
    let raw = hasher.finalize();
    Ok(executable_build_identity_from_input(
        &ExecutableBuildIdentityInput {
            byte_len: total,
            raw_hash: raw,
        },
        ROOFLINE_EXECUTABLE_DOMAIN,
    ))
}

fn require_stable_executable_identity(
    captured: fs_blake3::ContentHash,
    observed: fs_blake3::ContentHash,
) -> Result<fs_blake3::ContentHash, LedgerError> {
    if observed == captured {
        Ok(captured)
    } else {
        Err(LedgerError::Invalid {
            field: "executable_identity".to_string(),
            problem: format!(
                "current executable drifted between pre-measurement capture {captured} and recording {observed}"
            ),
        })
    }
}

#[derive(Debug)]
struct RooflineRowParams {
    op: i64,
    run_receipt: fs_blake3::ContentHash,
    payload_artifact: fs_blake3::ContentHash,
    dependency_receipt_artifact: Option<fs_blake3::ContentHash>,
    dependency_receipt_digest: Option<fs_blake3::ContentHash>,
    baseline_hash: fs_blake3::ContentHash,
    build_identity: fs_blake3::ContentHash,
    reps: u64,
    post_axis_bits: [u64; 4],
}

impl RooflineRowParams {
    fn to_json(&self) -> String {
        let dependency_artifact = self
            .dependency_receipt_artifact
            .map_or_else(|| "null".to_string(), |hash| format!("\"{hash}\""));
        let dependency_digest = self
            .dependency_receipt_digest
            .map_or_else(|| "null".to_string(), |hash| format!("\"{hash}\""));
        format!(
            "{{\"schema\":\"{ROOFLINE_ROW_SCHEMA}\",\"op\":{},\"run_receipt\":\"{}\",\"payload_artifact\":\"{}\",\"dependency_receipt_artifact\":{dependency_artifact},\"dependency_receipt_digest\":{dependency_digest},\"baseline_hash\":\"{}\",\"build_identity\":\"{}\",\"reps\":{},\"post_bandwidth_single_bits\":\"{:016x}\",\"post_bandwidth_all_core_bits\":\"{:016x}\",\"post_peak_single_bits\":\"{:016x}\",\"post_peak_all_core_bits\":\"{:016x}\"}}",
            self.op,
            self.run_receipt,
            self.payload_artifact,
            self.baseline_hash,
            self.build_identity,
            self.reps,
            self.post_axis_bits[0],
            self.post_axis_bits[1],
            self.post_axis_bits[2],
            self.post_axis_bits[3],
        )
    }
}

fn parse_roofline_row_params(text: &str) -> Option<RooflineRowParams> {
    fn take(rest: &mut &str, prefix: &str) -> Option<()> {
        *rest = rest.strip_prefix(prefix)?;
        Some(())
    }
    fn decimal(rest: &mut &str) -> Option<u64> {
        let end = rest
            .find(|ch: char| !ch.is_ascii_digit())
            .unwrap_or(rest.len());
        if end == 0 || (end > 1 && rest.starts_with('0')) {
            return None;
        }
        let (digits, tail) = rest.split_at(end);
        *rest = tail;
        digits.parse().ok()
    }
    fn hash(rest: &mut &str) -> Option<fs_blake3::ContentHash> {
        let (raw, tail) = rest.split_once('"')?;
        if raw.len() != 64 {
            return None;
        }
        *rest = tail;
        fs_blake3::ContentHash::from_hex(raw)
    }
    #[allow(clippy::option_option)] // parse failure, canonical null, and canonical hash are distinct
    fn optional_hash(rest: &mut &str) -> Option<Option<fs_blake3::ContentHash>> {
        if let Some(tail) = rest.strip_prefix("null") {
            *rest = tail;
            Some(None)
        } else {
            *rest = rest.strip_prefix('"')?;
            hash(rest).map(Some)
        }
    }
    fn bits(rest: &mut &str) -> Option<u64> {
        let (raw, tail) = rest.split_once('"')?;
        if raw.len() != 16
            || !raw
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return None;
        }
        *rest = tail;
        u64::from_str_radix(raw, 16).ok()
    }

    let mut rest = text;
    take(
        &mut rest,
        &format!("{{\"schema\":\"{ROOFLINE_ROW_SCHEMA}\",\"op\":"),
    )?;
    let op = i64::try_from(decimal(&mut rest)?).ok()?;
    if op <= 0 {
        return None;
    }
    take(&mut rest, ",\"run_receipt\":\"")?;
    let run_receipt = hash(&mut rest)?;
    take(&mut rest, ",\"payload_artifact\":\"")?;
    let payload_artifact = hash(&mut rest)?;
    take(&mut rest, ",\"dependency_receipt_artifact\":")?;
    let dependency_receipt_artifact = optional_hash(&mut rest)?;
    take(&mut rest, ",\"dependency_receipt_digest\":")?;
    let dependency_receipt_digest = optional_hash(&mut rest)?;
    take(&mut rest, ",\"baseline_hash\":\"")?;
    let baseline_hash = hash(&mut rest)?;
    take(&mut rest, ",\"build_identity\":\"")?;
    let build_identity = hash(&mut rest)?;
    take(&mut rest, ",\"reps\":")?;
    let reps = decimal(&mut rest)?;
    if reps == 0 {
        return None;
    }
    take(&mut rest, ",\"post_bandwidth_single_bits\":\"")?;
    let bandwidth_single = bits(&mut rest)?;
    take(&mut rest, ",\"post_bandwidth_all_core_bits\":\"")?;
    let bandwidth_all_core = bits(&mut rest)?;
    take(&mut rest, ",\"post_peak_single_bits\":\"")?;
    let peak_single = bits(&mut rest)?;
    take(&mut rest, ",\"post_peak_all_core_bits\":\"")?;
    let peak_all_core = bits(&mut rest)?;
    take(&mut rest, "}")?;
    if !rest.is_empty() {
        return None;
    }
    let params = RooflineRowParams {
        op,
        run_receipt,
        payload_artifact,
        dependency_receipt_artifact,
        dependency_receipt_digest,
        baseline_hash,
        build_identity,
        reps,
        post_axis_bits: [
            bandwidth_single,
            bandwidth_all_core,
            peak_single,
            peak_all_core,
        ],
    };
    (params.to_json() == text).then_some(params)
}

// ---------------------------------------------------------------------------
// Ledger integration and staleness
// ---------------------------------------------------------------------------

const EXTERNAL_PERF_GATE_ADMISSION_ARTIFACT_KIND: &str = "fs-roofline-axis-admission-receipt";
const EXTERNAL_PERF_GATE_ADMISSION_ARTIFACT_META: &str =
    "{\"schema\":\"fs-roofline-axis-admission-v2\",\"role\":\"external-perf-gate-input\"}";
const EXTERNAL_PERF_GATE_RESULT_ARTIFACT_KIND: &str = "fs-roofline-external-perf-gate";
const EXTERNAL_PERF_GATE_RESULT_ARTIFACT_META: &str =
    "{\"schema\":\"fs-roofline-external-perf-gate-v1\"}";
const EXTERNAL_PERF_GATE_SESSION: &[u8] = b"roofline-external-gate";
const EXTERNAL_PERF_GATE_SEED: &[u8] = b"external-perf-gate";
const MAX_EXTERNAL_PERF_GATE_FIELDS: usize = 64;
const MAX_EXTERNAL_PERF_GATE_NESTING: usize = 64;

/// External performance lane admitted by the centralized positive-gate
/// recorder.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExternalPerfGateLane {
    /// High-order FEEC sum-factorization lane.
    Feec,
    /// Memory-resident Stockham FFT lane.
    Fft,
}

impl ExternalPerfGateLane {
    /// Stable lower-case lane identity retained in the ledger receipt.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Feec => "feec",
            Self::Fft => "fft",
        }
    }

    const fn metric(self) -> &'static str {
        match self {
            Self::Feec => "feec-gate",
            Self::Fft => "fft-gate",
        }
    }

    const fn event_kind(self) -> &'static str {
        match self {
            Self::Feec => "external_perf_gate_feec",
            Self::Fft => "external_perf_gate_fft",
        }
    }
}

/// Receipt for one atomically retained external performance gate.
///
/// The operation and event ids are ledger-local. The two content hashes
/// identify the exact admission and final-gate bytes independently of those
/// local row ids.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExternalPerfGateReceipt {
    operation_id: i64,
    event_id: i64,
    lane: ExternalPerfGateLane,
    admission_artifact: fs_blake3::ContentHash,
    final_gate_artifact: fs_blake3::ContentHash,
    receipt: String,
}

impl ExternalPerfGateReceipt {
    fn new(
        operation_id: i64,
        event_id: i64,
        lane: ExternalPerfGateLane,
        admission_artifact: fs_blake3::ContentHash,
        final_gate_artifact: fs_blake3::ContentHash,
    ) -> Self {
        let receipt = format!(
            "{{\"schema\":\"fs-roofline-external-perf-gate-receipt-v1\",\"op\":{operation_id},\"event\":{event_id},\"lane\":\"{}\",\"admission_artifact\":\"{admission_artifact}\",\"final_gate_artifact\":\"{final_gate_artifact}\"}}",
            lane.as_str(),
        );
        Self {
            operation_id,
            event_id,
            lane,
            admission_artifact,
            final_gate_artifact,
            receipt,
        }
    }

    /// Ledger-local operation id linking both artifacts.
    #[must_use]
    pub const fn operation_id(&self) -> i64 {
        self.operation_id
    }

    /// Ledger-local row id of the lane-qualified diagnostic event.
    #[must_use]
    pub const fn event_id(&self) -> i64 {
        self.event_id
    }

    /// Typed external performance lane.
    #[must_use]
    pub const fn lane(&self) -> ExternalPerfGateLane {
        self.lane
    }

    /// Content identity of the exact axis-admission receipt bytes.
    #[must_use]
    pub const fn admission_artifact(&self) -> fs_blake3::ContentHash {
        self.admission_artifact
    }

    /// Content identity of the exact final-gate JSON bytes.
    #[must_use]
    pub const fn final_gate_artifact(&self) -> fs_blake3::ContentHash {
        self.final_gate_artifact
    }

    /// Canonical structured receipt for logging and later lookup.
    #[must_use]
    pub fn receipt_json(&self) -> &str {
        &self.receipt
    }
}

#[derive(Default)]
struct ExternalPerfGateFields<'a> {
    metric: Option<&'a str>,
    citation_eligible: Option<&'a str>,
    recorded: Option<&'a str>,
    report_only: Option<&'a str>,
    admission: Option<&'a str>,
}

fn external_gate_invalid(problem: impl Into<String>) -> LedgerError {
    LedgerError::Invalid {
        field: "external_perf_gate".to_string(),
        problem: problem.into(),
    }
}

fn skip_json_whitespace(bytes: &[u8], mut cursor: usize) -> usize {
    while bytes
        .get(cursor)
        .is_some_and(|byte| matches!(byte, b' ' | b'\n' | b'\r' | b'\t'))
    {
        cursor += 1;
    }
    cursor
}

fn scan_json_string(bytes: &[u8], start: usize) -> Option<usize> {
    if bytes.get(start) != Some(&b'"') {
        return None;
    }
    let mut cursor = start + 1;
    while let Some(byte) = bytes.get(cursor).copied() {
        match byte {
            b'"' => return Some(cursor + 1),
            b'\\' => {
                cursor = cursor.checked_add(2)?;
            }
            0x00..=0x1f => return None,
            _ => cursor += 1,
        }
    }
    None
}

fn scan_json_value(bytes: &[u8], start: usize) -> Option<usize> {
    match bytes.get(start).copied()? {
        b'"' => scan_json_string(bytes, start),
        b'{' | b'[' => {
            let mut stack = [0_u8; MAX_EXTERNAL_PERF_GATE_NESTING];
            let mut depth = 0_usize;
            let mut cursor = start;
            while let Some(byte) = bytes.get(cursor).copied() {
                match byte {
                    b'"' => cursor = scan_json_string(bytes, cursor)?,
                    b'{' | b'[' => {
                        if depth == stack.len() {
                            return None;
                        }
                        stack[depth] = byte;
                        depth = depth.checked_add(1)?;
                        cursor += 1;
                    }
                    b'}' | b']' => {
                        let opener = *stack.get(depth.checked_sub(1)?)?;
                        if !matches!((opener, byte), (b'{', b'}') | (b'[', b']')) {
                            return None;
                        }
                        depth = depth.checked_sub(1)?;
                        cursor += 1;
                        if depth == 0 {
                            return Some(cursor);
                        }
                    }
                    _ => cursor += 1,
                }
            }
            None
        }
        _ => {
            let mut cursor = start;
            while bytes
                .get(cursor)
                .is_some_and(|byte| !matches!(byte, b',' | b'}'))
            {
                cursor += 1;
            }
            let end = (start..cursor)
                .rev()
                .find(|index| !bytes[*index].is_ascii_whitespace())?
                + 1;
            (end > start).then_some(end)
        }
    }
}

fn set_external_gate_field<'a>(
    field: &mut Option<&'a str>,
    key: &str,
    value: &'a str,
) -> Result<(), LedgerError> {
    if field.replace(value).is_some() {
        return Err(external_gate_invalid(format!(
            "duplicate top-level {key:?} field is ambiguous"
        )));
    }
    Ok(())
}

fn parse_external_gate_fields(text: &str) -> Result<ExternalPerfGateFields<'_>, LedgerError> {
    let bytes = text.as_bytes();
    let mut cursor = skip_json_whitespace(bytes, 0);
    if bytes.get(cursor) != Some(&b'{') {
        return Err(external_gate_invalid(
            "final_gate_json must be one top-level JSON object",
        ));
    }
    cursor += 1;
    let mut fields = ExternalPerfGateFields::default();
    let mut field_count = 0_usize;
    loop {
        cursor = skip_json_whitespace(bytes, cursor);
        if bytes.get(cursor) == Some(&b'}') {
            cursor += 1;
            break;
        }
        field_count = field_count
            .checked_add(1)
            .ok_or_else(|| external_gate_invalid("top-level field count overflowed"))?;
        if field_count > MAX_EXTERNAL_PERF_GATE_FIELDS {
            return Err(external_gate_invalid(format!(
                "final_gate_json exceeds the {MAX_EXTERNAL_PERF_GATE_FIELDS}-field top-level bound"
            )));
        }
        let key_start = cursor;
        let key_end = scan_json_string(bytes, key_start)
            .ok_or_else(|| external_gate_invalid("malformed top-level JSON field name"))?;
        let key = &text[key_start + 1..key_end - 1];
        if key.is_empty()
            || !key
                .bytes()
                .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'_')
        {
            return Err(external_gate_invalid(
                "top-level field names must use canonical lower-case ASCII identifiers",
            ));
        }
        cursor = skip_json_whitespace(bytes, key_end);
        if bytes.get(cursor) != Some(&b':') {
            return Err(external_gate_invalid(format!(
                "top-level field {key:?} lacks a value separator"
            )));
        }
        cursor = skip_json_whitespace(bytes, cursor + 1);
        let value_start = cursor;
        let value_end = scan_json_value(bytes, value_start).ok_or_else(|| {
            external_gate_invalid(format!("top-level field {key:?} has no complete value"))
        })?;
        let value = text[value_start..value_end].trim_end();
        match key {
            "metric" => set_external_gate_field(&mut fields.metric, key, value)?,
            "citation_eligible" => {
                set_external_gate_field(&mut fields.citation_eligible, key, value)?;
            }
            "recorded" => set_external_gate_field(&mut fields.recorded, key, value)?,
            "report_only" => set_external_gate_field(&mut fields.report_only, key, value)?,
            "admission" => set_external_gate_field(&mut fields.admission, key, value)?,
            _ => {}
        }
        cursor = skip_json_whitespace(bytes, value_end);
        match bytes.get(cursor) {
            Some(b',') => cursor += 1,
            Some(b'}') => {
                cursor += 1;
                break;
            }
            _ => {
                return Err(external_gate_invalid(format!(
                    "top-level field {key:?} is not followed by a comma or object end"
                )));
            }
        }
    }
    cursor = skip_json_whitespace(bytes, cursor);
    if cursor != bytes.len() {
        return Err(external_gate_invalid(
            "bytes follow the final-gate JSON object",
        ));
    }
    Ok(fields)
}

fn validate_external_perf_gate_inputs(
    lane: ExternalPerfGateLane,
    admission: &AxisAdmissionSnapshot,
    final_gate_json: &str,
) -> Result<(), LedgerError> {
    if let Some(reason) = admission.baseline_citation_error() {
        return Err(external_gate_invalid(format!(
            "authority-admitted baseline snapshot required: {reason}"
        )));
    }
    if final_gate_json.len() > MAX_EXTERNAL_PERF_GATE_JSON_BYTES {
        return Err(external_gate_invalid(format!(
            "final_gate_json is {} bytes, exceeding the {MAX_EXTERNAL_PERF_GATE_JSON_BYTES}-byte bound",
            final_gate_json.len()
        )));
    }
    let fields = parse_external_gate_fields(final_gate_json)?;
    let expected_metric = format!("\"{}\"", lane.metric());
    if fields.metric != Some(expected_metric.as_str()) {
        return Err(external_gate_invalid(format!(
            "lane {:?} requires top-level metric {expected_metric}",
            lane
        )));
    }
    if fields.citation_eligible != Some("true") {
        return Err(external_gate_invalid(
            "top-level citation_eligible must be the literal true",
        ));
    }
    if fields.recorded != Some("true") {
        return Err(external_gate_invalid(
            "top-level recorded must be the literal true",
        ));
    }
    if fields.report_only != Some("false") {
        return Err(external_gate_invalid(
            "top-level report_only must be the literal false",
        ));
    }
    if fields.admission != Some(admission.receipt_json()) {
        return Err(external_gate_invalid(
            "top-level admission must preserve the exact supplied AxisAdmissionSnapshot receipt",
        ));
    }
    Ok(())
}

fn checked_external_gate_len(field: &str, bytes: &[u8]) -> Result<u64, LedgerError> {
    u64::try_from(bytes.len()).map_err(|_| {
        external_gate_invalid(format!("{field} byte length cannot be represented as u64"))
    })
}

fn require_exact_external_gate_artifact(
    ledger: &Ledger,
    field: &str,
    hash: &fs_blake3::ContentHash,
    expected: &[u8],
    expected_kind: &str,
    expected_meta: &str,
) -> Result<(), LedgerError> {
    let bound = checked_external_gate_len(field, expected)?;
    let info = ledger.artifact_info(hash)?.ok_or_else(|| {
        external_gate_invalid(format!("{field} artifact envelope disappeared after write"))
    })?;
    if info.hash != *hash
        || info.len != bound
        || info.kind != expected_kind
        || info.meta.as_deref() != Some(expected_meta)
    {
        return Err(external_gate_invalid(format!(
            "{field} artifact envelope differs on exact re-read"
        )));
    }
    let stored = ledger.get_artifact_bounded(hash, bound)?.ok_or_else(|| {
        external_gate_invalid(format!("{field} artifact disappeared after write"))
    })?;
    if stored.as_slice() != expected {
        return Err(external_gate_invalid(format!(
            "{field} artifact bytes differ from the exact supplied bytes"
        )));
    }
    Ok(())
}

fn rollback_external_gate_error(ledger: &Ledger, error: LedgerError) -> LedgerError {
    if !ledger.in_transaction() {
        return error;
    }
    match ledger.rollback() {
        Ok(()) => error,
        Err(rollback) => external_gate_invalid(format!(
            "external gate write failed ({error}); rollback also failed ({rollback})"
        )),
    }
}

fn external_gate_path_is_nondurable(path: &str) -> bool {
    path.eq_ignore_ascii_case(":memory:")
        || path
            .get(..5)
            .is_some_and(|prefix| prefix.eq_ignore_ascii_case("file:"))
}

/// Open a durable ledger path and atomically retain one authority-admitted
/// external performance gate.
///
/// Empty paths, SQLite's `:memory:` sentinel, and all `file:` URI forms are
/// refused: FEEC/FFT cannot emit recorded/citable evidence unless the receipt
/// survives the process. URI interpretation is deliberately outside this
/// narrow durable-path boundary.
/// The supplied gate must be a bounded JSON object for the selected lane with
/// literal positive citation/recording fields and the exact admission receipt.
///
/// # Errors
/// Invalid paths, non-authoritative admission, report-only or malformed gate
/// JSON, and all ledger write/read failures are returned as `LedgerError`.
pub fn record_external_perf_gate_at_path(
    path: &str,
    lane: ExternalPerfGateLane,
    admission: &AxisAdmissionSnapshot,
    final_gate_json: &str,
) -> Result<ExternalPerfGateReceipt, LedgerError> {
    let trimmed = path.trim();
    if trimmed.is_empty() || external_gate_path_is_nondurable(trimmed) {
        return Err(external_gate_invalid(
            "a non-empty non-URI durable ledger path is required; SQLite memory/URI modes are refused",
        ));
    }
    // Refuse authority, lane, and bounded structural mismatches before opening
    // (and potentially creating) a path. The Ledger-backed core repeats these
    // checks, and begin_op supplies the engine's complete JSON validation.
    validate_external_perf_gate_inputs(lane, admission, final_gate_json)?;
    let ledger = Ledger::open(path)?;
    record_external_perf_gate_in_ledger(&ledger, lane, admission, final_gate_json)
}

/// Atomically retain one authority-admitted external performance gate in an
/// already-open ledger.
///
/// This Ledger-typed entry point is the in-memory test seam. Production FEEC
/// and FFT callers use [`record_external_perf_gate_at_path`] and therefore do
/// not need a direct fs-ledger dependency. The recorder owns its transaction;
/// an already-open caller transaction is refused rather than committed or
/// rolled back. Both exact artifacts and their role-qualified edges, plus the
/// completed operation envelope, are re-read before commit. A lane-qualified
/// non-authoritative diagnostic event projects that op and the same
/// artifacts/gate. Its positive row id and exact one-row count increment prove
/// insertion only: fs-ledger intentionally exposes no general event-payload
/// reader, so the exact artifact and op IR remain the byte-authoritative
/// records. Admission is checked again after those writes so a UTC-day
/// rollover aborts the complete transaction.
///
/// # Errors
/// Non-authoritative or mismatched inputs, caller transaction ownership,
/// write/read disagreement, and ledger failures are all fail-closed errors.
#[allow(clippy::too_many_lines)] // one visible all-or-nothing evidence boundary
pub fn record_external_perf_gate_in_ledger(
    ledger: &Ledger,
    lane: ExternalPerfGateLane,
    admission: &AxisAdmissionSnapshot,
    final_gate_json: &str,
) -> Result<ExternalPerfGateReceipt, LedgerError> {
    validate_external_perf_gate_inputs(lane, admission, final_gate_json)?;
    if ledger.in_transaction() {
        return Err(external_gate_invalid(
            "recorder must own its transaction; commit or roll back the caller transaction first",
        ));
    }

    let admission_bytes = admission.receipt_json().as_bytes();
    let gate_bytes = final_gate_json.as_bytes();
    let admission_hash = fs_ledger::hash_bytes(admission_bytes);
    let final_gate_hash = fs_ledger::hash_bytes(gate_bytes);
    if admission_hash == final_gate_hash {
        return Err(external_gate_invalid(
            "admission and final-gate artifacts must have distinct content identities",
        ));
    }
    let admission_len = checked_external_gate_len("admission", admission_bytes)?;
    let final_gate_len = checked_external_gate_len("final_gate_json", gate_bytes)?;
    let total_bytes = admission_len
        .checked_add(final_gate_len)
        .ok_or_else(|| external_gate_invalid("combined artifact byte count overflowed"))?;
    let versions = format!("{{\"fs-roofline\":\"{VERSION}\",\"external_perf_gate_protocol\":1}}");
    let budget = format!("{{\"artifact_bytes\":{total_bytes}}}");
    let capability = format!(
        "{{\"ops\":[\"perf.external-gate\"],\"lane\":\"{}\"}}",
        lane.as_str()
    );
    // Embedding the exact final gate delegates full JSON validation to the
    // same engine that stores the operation, while the content hash and Out
    // edge bind its independently retained artifact.
    let ir = format!(
        "{{\"op\":\"perf.external-gate\",\"schema\":\"fs-roofline-external-perf-gate-op-v1\",\"lane\":\"{}\",\"admission_artifact\":\"{admission_hash}\",\"final_gate_artifact\":\"{final_gate_hash}\",\"final_gate\":{final_gate_json}}}",
        lane.as_str(),
    );
    let explicits = FiveExplicits {
        seed: EXTERNAL_PERF_GATE_SEED,
        versions: &versions,
        budget: &budget,
        capability: &capability,
    };
    let recorded_at = now_wall_ns();
    let decision_day = admission
        .decision_day()
        .ok_or_else(|| external_gate_invalid("authority-admitted snapshot lacks a decision day"))?;
    if wall_ns_day(recorded_at) != Some(decision_day) {
        return Err(external_gate_invalid(
            "external gate start timestamp is outside the attested admission day",
        ));
    }

    ledger.begin()?;
    if !ledger.in_transaction() {
        return Err(external_gate_invalid(
            "ledger reported a successful begin without an open transaction",
        ));
    }
    let write_result: Result<(i64, i64), LedgerError> = (|| {
        let op = ledger.begin_op(
            Some(EXTERNAL_PERF_GATE_SESSION),
            &ir,
            &explicits,
            recorded_at,
        )?;
        let stored_admission = ledger.put_artifact(
            EXTERNAL_PERF_GATE_ADMISSION_ARTIFACT_KIND,
            admission_bytes,
            Some(EXTERNAL_PERF_GATE_ADMISSION_ARTIFACT_META),
        )?;
        if stored_admission.hash != admission_hash || stored_admission.len != admission_len {
            return Err(external_gate_invalid(
                "stored admission artifact receipt differs from the exact write request",
            ));
        }
        let stored_gate = ledger.put_artifact(
            EXTERNAL_PERF_GATE_RESULT_ARTIFACT_KIND,
            gate_bytes,
            Some(EXTERNAL_PERF_GATE_RESULT_ARTIFACT_META),
        )?;
        if stored_gate.hash != final_gate_hash || stored_gate.len != final_gate_len {
            return Err(external_gate_invalid(
                "stored final-gate artifact receipt differs from the exact write request",
            ));
        }
        ledger.link(op, &admission_hash, EdgeRole::In)?;
        ledger.link(op, &final_gate_hash, EdgeRole::Out)?;
        require_exact_external_gate_artifact(
            ledger,
            "admission",
            &admission_hash,
            admission_bytes,
            EXTERNAL_PERF_GATE_ADMISSION_ARTIFACT_KIND,
            EXTERNAL_PERF_GATE_ADMISSION_ARTIFACT_META,
        )?;
        require_exact_external_gate_artifact(
            ledger,
            "final_gate_json",
            &final_gate_hash,
            gate_bytes,
            EXTERNAL_PERF_GATE_RESULT_ARTIFACT_KIND,
            EXTERNAL_PERF_GATE_RESULT_ARTIFACT_META,
        )?;
        if !ledger.edge_exists(op, &admission_hash, EdgeRole::In)?
            || !ledger.edge_exists(op, &final_gate_hash, EdgeRole::Out)?
        {
            return Err(external_gate_invalid(
                "role-qualified admission/final-gate lineage edge failed exact re-read",
            ));
        }
        let event_count_before = ledger.table_count("events")?;
        let event_payload = format!(
            "{{\"schema\":\"fs-roofline-external-perf-gate-event-v1\",\"op\":{op},\"lane\":\"{}\",\"admission_artifact\":\"{admission_hash}\",\"final_gate_artifact\":\"{final_gate_hash}\",\"final_gate\":{final_gate_json}}}",
            lane.as_str(),
        );
        let event_id = ledger.append_event(&EventRow {
            session: Some(EXTERNAL_PERF_GATE_SESSION),
            t: recorded_at,
            kind: lane.event_kind(),
            payload: Some(&event_payload),
        })?;
        let expected_event_count = event_count_before
            .checked_add(1)
            .ok_or_else(|| external_gate_invalid("event count overflowed"))?;
        if event_id <= 0 || ledger.table_count("events")? != expected_event_count {
            return Err(external_gate_invalid(
                "lane-qualified diagnostic event write receipt failed its precommit existence check",
            ));
        }
        let completed_at = now_wall_ns();
        if completed_at < recorded_at
            || wall_ns_day(completed_at) != Some(decision_day)
            || wall_ns_day(recorded_at) != Some(decision_day)
        {
            return Err(external_gate_invalid(
                "external gate timestamps must be monotone and remain within the attested admission day",
            ));
        }
        ledger.finish_op(op, OpOutcome::Ok, None, completed_at)?;
        let stored_op = ledger.op(op)?.ok_or_else(|| {
            external_gate_invalid("external performance-gate operation disappeared after write")
        })?;
        if stored_op.session.as_deref() != Some(EXTERNAL_PERF_GATE_SESSION)
            || stored_op.ir != ir
            || stored_op.seed.as_slice() != EXTERNAL_PERF_GATE_SEED
            || stored_op.versions != versions
            || stored_op.budget != budget
            || stored_op.capability != capability
            || stored_op.t_start != recorded_at
            || stored_op.t_end != Some(completed_at)
            || stored_op.outcome.as_deref() != Some("ok")
            || stored_op.diag.is_some()
        {
            return Err(external_gate_invalid(
                "completed external performance-gate operation differs on exact re-read",
            ));
        }
        if let Some(reason) = admission.baseline_citation_error() {
            return Err(external_gate_invalid(format!(
                "admission ceased to be citation-eligible before commit: {reason}"
            )));
        }
        if !ledger.in_transaction() {
            return Err(external_gate_invalid(
                "ledger transaction ended before the atomic gate commit",
            ));
        }
        Ok((op, event_id))
    })();

    let (op, event_id) = match write_result {
        Ok(receipt) => receipt,
        Err(error) => return Err(rollback_external_gate_error(ledger, error)),
    };
    if let Err(error) = ledger.commit() {
        return Err(rollback_external_gate_error(ledger, error));
    }
    if ledger.in_transaction() {
        let error = external_gate_invalid(
            "ledger reported a successful commit but retained an open transaction",
        );
        return Err(rollback_external_gate_error(ledger, error));
    }
    Ok(ExternalPerfGateReceipt::new(
        op,
        event_id,
        lane,
        admission_hash,
        final_gate_hash,
    ))
}

/// Record a harness run atomically in the ledger. Admitted timed rows receive
/// metrics, `benchmark_result` events, content-addressed output artifacts, and
/// versioned tune rows. Rejected input receives an exact baseline-admission
/// receipt plus one explicit rejection event and no normal-looking
/// measurements.
/// A successful commit consumes each result's one-shot fresh-row marker;
/// rollback retains it so the same transaction can be retried.
///
/// Evidence recorded through this public entry point is stamped
/// `"protocol":"custom-registry"` in the operation `ir` (bead fz2.5): the
/// caller supplied the registry and both axis probes, so nothing proves the
/// kernels are the production set or that the post-probe was observed after
/// the timed repetitions. Custom-registry evidence is explicitly
/// NON-CITABLE for performance claims; the sealed
/// [`production::ProductionRun`] protocol is the only path that records
/// `"protocol":"production-v3"` together with a retained dependency receipt.
///
/// # Errors
/// Ledger errors propagate and roll back the whole write set.
pub fn record_run(
    ledger: &Ledger,
    axes: &MachineAxes,
    post_axes: &MachineAxes,
    baseline: AxisBaselinePolicy<'_>,
    finalized: &mut FinalizedRegistryRun,
    results: &mut [Attainment],
) -> Result<i64, LedgerError> {
    let admission = baseline.snapshot(axes, post_axes);
    record_run_with_protocol(
        ledger,
        axes,
        post_axes,
        &admission,
        finalized,
        results,
        CUSTOM_REGISTRY_PROTOCOL_FIELD,
        None,
        None,
        EvidenceNamespace::Custom,
        None,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // one auditable all-or-nothing evidence transaction
pub(crate) fn record_run_with_protocol(
    ledger: &Ledger,
    axes: &MachineAxes,
    post_axes: &MachineAxes,
    admission: &AxisAdmissionSnapshot,
    finalized: &mut FinalizedRegistryRun,
    results: &mut [Attainment],
    protocol_ir_fields: &str,
    dependency_receipt: Option<DependencyReceiptBinding>,
    citation_refusal: Option<&str>,
    evidence_namespace: EvidenceNamespace,
    captured_build_identity: Option<fs_blake3::ContentHash>,
) -> Result<i64, LedgerError> {
    let live_day_error = admission.live_day_error();
    let measurement_error = run_admission_error_for_snapshot_with_live_error(
        axes,
        admission,
        results,
        live_day_error.clone(),
    );
    let measurement_admitted = measurement_error.is_none();
    let delayed_day_refusal = finalized.admitted
        && !measurement_admitted
        && live_day_error.is_some()
        && measurement_error == live_day_error;
    let admission_error = match (citation_refusal, measurement_error.as_deref()) {
        (Some(forced), Some(measurement)) if forced != measurement => {
            Some(format!("{forced}; {measurement}"))
        }
        (Some(forced), _) => Some(forced.to_string()),
        (None, Some(measurement)) => Some(measurement.to_string()),
        (None, None) => None,
    };
    let run_valid = admission_error.is_none();
    if finalized.consumed {
        return Err(LedgerError::Invalid {
            field: "finalized_run".to_string(),
            problem: "the finalized roofline run was already recorded".to_string(),
        });
    }
    let expected_receipt = finalized_run_receipt(admission, results);
    if finalized.receipt != expected_receipt
        || (finalized.admitted != measurement_admitted && !delayed_day_refusal)
    {
        return Err(LedgerError::Invalid {
            field: "finalized_run".to_string(),
            problem:
                "axes, baseline decision, results, or admission changed after registry finalization"
                    .to_string(),
        });
    }
    let baseline_receipt = admission.receipt_json().to_string();
    let build_identity = if let Some(captured) = captured_build_identity {
        let observed = read_executable_build_identity()?;
        require_stable_executable_identity(captured, observed)?
    } else {
        executable_build_identity()?
    };
    let versions = versions_json(build_identity);
    let explicits = FiveExplicits {
        seed: ROOFLINE_SEED,
        versions: &versions,
        budget: ROOFLINE_BUDGET,
        capability: ROOFLINE_CAPABILITY,
    };
    // The protocol stamp sits between `admitted` and the receipt/manifest
    // tail; `baseline_admission` must stay the final field (staleness
    // extracts the baseline receipt bytes by stripping the closing brace).
    let ir = format!(
        "{{\"op\":\"perf.roofline\",\"kernels\":{},\"fingerprint\":\"{:016x}\",\"post_fingerprint\":\"{:016x}\",\"measurement_admitted\":{measurement_admitted},\"admitted\":{run_valid},{protocol_ir_fields},\"finalized_run_receipt\":\"{}\",\"result_manifest\":{},\"baseline_admission\":{baseline_receipt}}}",
        results.len(),
        axes.fingerprint,
        post_axes.fingerprint,
        finalized.receipt,
        run_result_manifest_json(results),
    );
    ledger.begin()?;
    let write_result: Result<i64, LedgerError> = (|| {
        let op = ledger.begin_op(Some(b"roofline"), &ir, &explicits, now_wall_ns())?;
        if run_valid && let Some(binding) = dependency_receipt {
            let stored = ledger.put_artifact(
                DEPGRAPH_RECEIPT_ARTIFACT_KIND,
                binding.bytes.as_bytes(),
                Some(DEPGRAPH_RECEIPT_ARTIFACT_META),
            )?;
            if stored.hash != binding.artifact_hash {
                return Err(LedgerError::Invalid {
                    field: "dependency_receipt".to_string(),
                    problem: "stored dependency receipt identity differs from the protocol binding"
                        .to_string(),
                });
            }
            ledger.link(op, &stored.hash, EdgeRole::In)?;
        }
        ledger.append_event(&EventRow {
            session: Some(b"roofline"),
            t: 0,
            kind: "axis_baseline_admission",
            payload: Some(&baseline_receipt),
        })?;
        if run_valid {
            let baseline_hash = admission
                .baseline_hash()
                .ok_or_else(|| LedgerError::Invalid {
                    field: "baseline".to_string(),
                    problem: "trusted roofline admission has no selected baseline".to_string(),
                })?;
            let machine_key = roofline_machine_key(axes.fingerprint, baseline_hash);
            for r in results.iter() {
                if let MeasurementOrigin::Timed {
                    decision_bindings, ..
                } = &r.measurement_origin
                    && let Some(binding) = stable_decision_binding(decision_bindings)
                {
                    // Publish the exact row from the sealed binding even when
                    // this measurement adopted it from another ledger. This
                    // evidence transaction is deliberately insert-only: a
                    // clone or delayed receipt must never replace a newer
                    // cache row. Refresh/replacement is a separate explicit
                    // protocol, never authority carried by Attainment.
                    binding
                        .validated_row()
                        .publish_if_absent_or_identical(ledger)?;
                }
                ledger.record_metric(
                    op,
                    0,
                    &format!("{}.elems_per_sec", r.kernel),
                    r.elems_per_sec,
                )?;
                ledger.record_metric(op, 0, &format!("{}.attainment", r.kernel), r.attainment)?;
                ledger.record_metric(op, 0, &format!("{}.dispersion", r.kernel), r.dispersion)?;
                let payload = r.to_jsonl();
                ledger.append_event(&EventRow {
                    session: Some(b"roofline"),
                    t: 0,
                    kind: "benchmark_result",
                    payload: Some(&payload),
                })?;
                let payload_artifact = ledger.put_artifact(
                    ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                    payload.as_bytes(),
                    Some(ROOFLINE_PAYLOAD_ARTIFACT_META),
                )?;
                ledger.link(op, &payload_artifact.hash, EdgeRole::Out)?;
                let params = RooflineRowParams {
                    op,
                    run_receipt: finalized.receipt,
                    payload_artifact: payload_artifact.hash,
                    dependency_receipt_artifact: dependency_receipt
                        .map(|binding| binding.artifact_hash),
                    dependency_receipt_digest: dependency_receipt
                        .map(|binding| binding.domain_digest),
                    baseline_hash,
                    build_identity,
                    reps: u64::try_from(r.reps).map_err(|_| LedgerError::Invalid {
                        field: "reps".to_string(),
                        problem: "roofline repetition count exceeds u64".to_string(),
                    })?,
                    post_axis_bits: [
                        post_axes.bandwidth_single_gbs.to_bits(),
                        post_axes.bandwidth_all_core_gbs.to_bits(),
                        post_axes.peak_single_gflops.to_bits(),
                        post_axes.peak_all_core_gflops.to_bits(),
                    ],
                }
                .to_json();
                let shape_class = match evidence_namespace {
                    EvidenceNamespace::Production => {
                        tune_measurement_shape_class(&r.version, finalized.receipt, op)
                    }
                    EvidenceNamespace::Custom => {
                        candidate_measurement_shape_class(&r.version, finalized.receipt, op)
                    }
                };
                ledger.tune_put_if_absent(
                    &r.kernel,
                    &shape_class,
                    &machine_key,
                    &params,
                    &payload,
                )?;
                let stored = ledger
                    .tune_get(&r.kernel, &shape_class, &machine_key)?
                    .ok_or_else(|| LedgerError::Invalid {
                        field: "tune".to_string(),
                        problem: "roofline insert-if-absent returned without a stored row"
                            .to_string(),
                    })?;
                if stored.params != params || stored.measured != payload {
                    return Err(LedgerError::Invalid {
                        field: "tune".to_string(),
                        problem: format!(
                            "refusing conflicting roofline evidence for kernel {:?}, shape {:?}",
                            r.kernel, shape_class
                        ),
                    });
                }
            }
            let completed_at = now_wall_ns();
            if evidence_namespace == EvidenceNamespace::Production
                && wall_ns_day(completed_at) != admission.decision_day()
            {
                return Err(LedgerError::Invalid {
                    field: "baseline_admission.decision_day".to_string(),
                    problem: "roofline recording crossed the attested admission day boundary"
                        .to_string(),
                });
            }
            ledger.finish_op(op, OpOutcome::Ok, None, completed_at)?;
        } else {
            let reason = admission_error
                .as_deref()
                .unwrap_or("unknown admission failure");
            let payload = format!(
                "{{\"code\":\"roofline_evidence_rejected\",\"reason\":\"{}\",\"effect\":\"no_measurements_or_tune_rows_published\"}}",
                json_escape(reason)
            );
            ledger.append_event(&EventRow {
                session: Some(b"roofline"),
                t: 0,
                kind: "roofline_run_rejected",
                payload: Some(&payload),
            })?;
            let diagnostic = format!(
                "{{\"code\":\"roofline_evidence_rejected\",\"reason\":\"{}\"}}",
                json_escape(reason)
            );
            ledger.finish_op(op, OpOutcome::Error, Some(&diagnostic), now_wall_ns())?;
        }
        Ok(op)
    })();
    match write_result {
        Ok(op) => match ledger.commit() {
            Ok(()) => {
                finalized.consumed = true;
                for result in results.iter_mut() {
                    result.pending_tune_publication = None;
                }
                Ok(op)
            }
            Err(error) => {
                let _ = ledger.rollback();
                Err(error)
            }
        },
        Err(error) => {
            let _ = ledger.rollback();
            Err(error)
        }
    }
}

/// Staleness state of one kernel's ledgered attainment on this machine.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Staleness {
    /// At least one semantically valid row from this exact executable is no
    /// more than [`STALENESS_MAX_AGE_NS`] old.
    Fresh,
    /// Matching current-build evidence exists, but its newest successful
    /// operation is older than [`STALENESS_MAX_AGE_NS`].
    Expired,
    /// The observation clock predates the newest matching operation receipt.
    ClockRollback,
    /// Rows exist, but none for the current fingerprint — the machine
    /// drifted and every cited number is stale until re-measured.
    FingerprintDrift,
    /// The current machine fingerprint matches historical rows, but no trusted
    /// baseline was selected for comparison.
    BaselineUnavailable,
    /// The machine matches but the selected historical baseline changed.
    BaselineDrift,
    /// Semantically valid rows exist for this machine and baseline, but only
    /// for a different executable-content identity.
    BuildDrift,
    /// A row exists under the exact current version, machine, and baseline key,
    /// but its canonical parameters, payload artifact/edge, operation receipt,
    /// or executable identity no longer agree. Corrupt evidence is never
    /// treated as fresh.
    CorruptEvidence,
    /// No roofline rows at all: never measured.
    NeverMeasured,
}

/// Check one kernel's staleness against the current fingerprint and admitted
/// historical baseline identity.
///
/// # Errors
/// Ledger errors propagate.
pub fn staleness(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
) -> Result<Staleness, LedgerError> {
    staleness_at(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
        now_wall_ns(),
    )
}

/// Deterministic form of [`staleness`] evaluated at an explicit wall-clock
/// nanosecond. Supplying time makes expiry and rollback tests replayable.
///
/// # Errors
/// Ledger and executable-identity errors propagate.
pub fn staleness_at(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
    observed_wall_ns: i64,
) -> Result<Staleness, LedgerError> {
    let current_build = executable_build_identity()?;
    staleness_at_with_build(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
        observed_wall_ns,
        current_build,
    )
}

fn staleness_at_with_build(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
    observed_wall_ns: i64,
    current_build: fs_blake3::ContentHash,
) -> Result<Staleness, LedgerError> {
    staleness_at_with_build_and_dependency(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
        observed_wall_ns,
        current_build,
        DependencyReceiptBinding::current().ok(),
    )
}

#[allow(clippy::too_many_arguments)]
/// Shared row-selection prefix of the staleness lattice: everything decided
/// BEFORE per-row validation (bead vm3i — the checkpointed fast path must
/// classify NeverMeasured/FingerprintDrift/BaselineUnavailable/BaselineDrift
/// identically to the exhaustive path, so both call this).
pub(crate) enum RowSelection {
    Verdict(Staleness),
    Rows(Vec<fs_ledger::TuneRow>),
}

pub(crate) fn select_matching_rows(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
) -> Result<RowSelection, LedgerError> {
    let rows = match ledger.tune_rows(kernel) {
        Ok(rows) => rows,
        Err(LedgerError::TuneCorrupt { .. }) => {
            return Ok(RowSelection::Verdict(Staleness::CorruptEvidence));
        }
        Err(error) => return Err(error),
    };
    let shape_prefix = format!("{}:run=", tune_shape_class(version));
    let roofline_rows: Vec<_> = rows
        .into_iter()
        .filter(|r| r.shape_class.starts_with(&shape_prefix))
        .collect();
    if roofline_rows.is_empty() {
        return Ok(RowSelection::Verdict(Staleness::NeverMeasured));
    }
    let fp = current_fingerprint.to_le_bytes();
    let same_machine: Vec<_> = roofline_rows
        .into_iter()
        .filter(|row| row.machine.get(..8) == Some(fp.as_slice()))
        .collect();
    if same_machine.is_empty() {
        return Ok(RowSelection::Verdict(Staleness::FingerprintDrift));
    }
    let Some(current_baseline) = current_baseline else {
        return Ok(RowSelection::Verdict(Staleness::BaselineUnavailable));
    };
    let key = roofline_machine_key(current_fingerprint, current_baseline);
    let matching: Vec<_> = same_machine
        .into_iter()
        .filter(|row| row.machine == key)
        .collect();
    if matching.is_empty() {
        return Ok(RowSelection::Verdict(Staleness::BaselineDrift));
    }
    Ok(RowSelection::Rows(matching))
}

/// Final age/rollback classification over a completed build scan (bead vm3i:
/// shared by the exhaustive and checkpointed paths).
pub(crate) fn classify_scanned_rows(build_scan: BuildRowScan, observed_wall_ns: i64) -> Staleness {
    let Some(recorded_at_ns) = build_scan.newest_current_build else {
        return if build_scan.saw_foreign_build {
            Staleness::BuildDrift
        } else {
            Staleness::CorruptEvidence
        };
    };
    if observed_wall_ns < recorded_at_ns {
        return Staleness::ClockRollback;
    }
    if observed_wall_ns.saturating_sub(recorded_at_ns) > STALENESS_MAX_AGE_NS {
        return Staleness::Expired;
    }
    Staleness::Fresh
}

#[allow(clippy::too_many_arguments)] // Exact roofline row identity and dependency receipt key.
fn staleness_at_with_build_and_dependency(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
    observed_wall_ns: i64,
    current_build: fs_blake3::ContentHash,
    expected_dependency: Option<DependencyReceiptBinding>,
) -> Result<Staleness, LedgerError> {
    let matching_rows = match select_matching_rows(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
    )? {
        RowSelection::Verdict(v) => return Ok(v),
        RowSelection::Rows(rows) => rows,
    };
    // select_matching_rows only returns Rows after the baseline gate.
    let current_baseline = current_baseline.expect("baseline present when rows match");
    let mut build_scan = BuildRowScan::default();
    for row in &matching_rows {
        let Some(validated) = validate_roofline_row(
            ledger,
            row,
            kernel,
            version,
            current_fingerprint,
            current_baseline,
            expected_dependency,
        )?
        else {
            return Ok(Staleness::CorruptEvidence);
        };
        if !build_scan.observe(&validated, current_build) {
            return Ok(Staleness::CorruptEvidence);
        }
    }
    Ok(classify_scanned_rows(build_scan, observed_wall_ns))
}

pub(crate) struct ValidatedRooflineRow {
    pub(crate) build_identity: fs_blake3::ContentHash,
    pub(crate) recorded_at_ns: i64,
    pub(crate) dependency_matches_current: bool,
    /// The row's op-bound dependency-receipt digests (bead vm3i): the
    /// checkpoint stores these so the fast path can re-derive
    /// `dependency_matches_current` against a FUTURE current binding
    /// without refetching the op.
    pub(crate) dependency_receipt_digest: fs_blake3::ContentHash,
    pub(crate) dependency_receipt_artifact: fs_blake3::ContentHash,
}

#[derive(Default, Clone, Copy)]
pub(crate) struct BuildRowScan {
    pub(crate) newest_current_build: Option<i64>,
    pub(crate) saw_foreign_build: bool,
}

impl BuildRowScan {
    /// Returns false only when a row claims the current executable while its
    /// dependency receipt does not match the receipt compiled into that
    /// executable. Foreign-build rows remain valid history when their own
    /// retained receipt is structurally sound.
    fn observe(
        &mut self,
        row: &ValidatedRooflineRow,
        current_build: fs_blake3::ContentHash,
    ) -> bool {
        if row.build_identity == current_build {
            if !row.dependency_matches_current {
                return false;
            }
            self.newest_current_build = Some(
                self.newest_current_build
                    .map_or(row.recorded_at_ns, |newest| newest.max(row.recorded_at_ns)),
            );
        } else {
            self.saw_foreign_build = true;
        }
        true
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CanonicalProductionOp {
    kernel_count: u64,
    fingerprint: u64,
    post_fingerprint: u64,
    run_nonce: fs_blake3::ContentHash,
    pre_axes_receipt: fs_blake3::ContentHash,
    post_axes_receipt: fs_blake3::ContentHash,
    dependency_receipt_digest: fs_blake3::ContentHash,
    dependency_receipt_artifact: fs_blake3::ContentHash,
    finalized_run_receipt: fs_blake3::ContentHash,
    result_manifest: String,
    baseline_admission: String,
}

impl CanonicalProductionOp {
    fn to_json(&self) -> String {
        format!(
            "{{\"op\":\"perf.roofline\",\"kernels\":{},\"fingerprint\":\"{:016x}\",\"post_fingerprint\":\"{:016x}\",\"measurement_admitted\":true,\"admitted\":true,{PRODUCTION_PROTOCOL_FIELD},\"run_nonce\":\"{}\",\"pre_axes_receipt\":\"{}\",\"post_axes_receipt\":\"{}\",\"dependency_graph_evidence\":\"operator-observed-receipt\",\"dependency_receipt_digest\":\"{}\",\"dependency_receipt_artifact\":\"{}\",\"finalized_run_receipt\":\"{}\",\"result_manifest\":{},\"baseline_admission\":{}}}",
            self.kernel_count,
            self.fingerprint,
            self.post_fingerprint,
            self.run_nonce,
            self.pre_axes_receipt,
            self.post_axes_receipt,
            self.dependency_receipt_digest,
            self.dependency_receipt_artifact,
            self.finalized_run_receipt,
            self.result_manifest,
            self.baseline_admission,
        )
    }
}

fn parse_canonical_production_op(ir: &str) -> Option<CanonicalProductionOp> {
    fn take(rest: &mut &str, prefix: &str) -> Option<()> {
        *rest = rest.strip_prefix(prefix)?;
        Some(())
    }
    fn decimal(rest: &mut &str) -> Option<u64> {
        let end = rest
            .find(|character: char| !character.is_ascii_digit())
            .unwrap_or(rest.len());
        if end == 0 || (end > 1 && rest.starts_with('0')) {
            return None;
        }
        let (digits, tail) = rest.split_at(end);
        *rest = tail;
        digits.parse().ok()
    }
    fn hash(rest: &mut &str) -> Option<fs_blake3::ContentHash> {
        let (hex, tail) = rest.split_once('"')?;
        let hash = fs_blake3::ContentHash::from_hex(hex)?;
        *rest = tail;
        Some(hash)
    }
    fn hex_u64(rest: &mut &str) -> Option<u64> {
        let (raw, tail) = rest.split_once('"')?;
        if raw.len() != 16
            || !raw
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return None;
        }
        *rest = tail;
        u64::from_str_radix(raw, 16).ok()
    }

    if ir.len() > fs_ledger::MAX_OP_IR_BYTES {
        return None;
    }
    let mut rest = ir;
    take(&mut rest, "{\"op\":\"perf.roofline\",\"kernels\":")?;
    let kernel_count = decimal(&mut rest)?;
    if kernel_count == 0 {
        return None;
    }
    take(&mut rest, ",\"fingerprint\":\"")?;
    let fingerprint = hex_u64(&mut rest)?;
    take(&mut rest, ",\"post_fingerprint\":\"")?;
    let post_fingerprint = hex_u64(&mut rest)?;
    take(
        &mut rest,
        &format!(
            ",\"measurement_admitted\":true,\"admitted\":true,{PRODUCTION_PROTOCOL_FIELD},\"run_nonce\":\""
        ),
    )?;
    let run_nonce = hash(&mut rest)?;
    take(&mut rest, ",\"pre_axes_receipt\":\"")?;
    let pre_axes_receipt = hash(&mut rest)?;
    take(&mut rest, ",\"post_axes_receipt\":\"")?;
    let post_axes_receipt = hash(&mut rest)?;
    take(
        &mut rest,
        ",\"dependency_graph_evidence\":\"operator-observed-receipt\",\"dependency_receipt_digest\":\"",
    )?;
    let dependency_receipt_digest = hash(&mut rest)?;
    take(&mut rest, ",\"dependency_receipt_artifact\":\"")?;
    let dependency_receipt_artifact = hash(&mut rest)?;
    take(&mut rest, ",\"finalized_run_receipt\":\"")?;
    let finalized_run_receipt = hash(&mut rest)?;
    take(&mut rest, ",\"result_manifest\":")?;
    let (result_manifest, baseline_and_end) = rest.split_once(",\"baseline_admission\":")?;
    let baseline_admission = baseline_and_end.strip_suffix('}')?;
    let entries = parse_result_manifest(result_manifest)?;
    if u64::try_from(entries.len()).ok()? != kernel_count {
        return None;
    }
    let parsed = CanonicalProductionOp {
        kernel_count,
        fingerprint,
        post_fingerprint,
        run_nonce,
        pre_axes_receipt,
        post_axes_receipt,
        dependency_receipt_digest,
        dependency_receipt_artifact,
        finalized_run_receipt,
        result_manifest: result_manifest.to_string(),
        baseline_admission: baseline_admission.to_string(),
    };
    (parsed.to_json() == ir).then_some(parsed)
}

struct BaselineReceiptView<'a> {
    now_day: u64,
    decision_day: u64,
    identity: &'a str,
    pre: &'a str,
    post: &'a str,
    baseline_hash: fs_blake3::ContentHash,
    baseline: &'a str,
    required_sources: Vec<fs_blake3::ContentHash>,
}

fn parse_baseline_receipt_view(text: &str) -> Option<BaselineReceiptView<'_>> {
    let prefix = "{\"schema\":\"fs-roofline-axis-admission-v2\",\"tier\":\"attested\",\"now_day\":";
    if !text.starts_with(prefix)
        || text.matches(",\"pre\":").count() != 1
        || text.matches(",\"post\":").count() != 1
        || text.matches(",\"baseline_hash\":").count() != 1
        || text.matches(",\"attestation\":").count() != 1
        || text.matches(",\"required_source_receipts\":").count() != 1
        || text.matches(",\"authority\":").count() != 1
        || text.matches(",\"verdict\":").count() != 1
    {
        return None;
    }
    let day_and_later = text.strip_prefix(prefix)?;
    let (now_day_text, decision_and_later) = day_and_later.split_once(",\"decision_day\":")?;
    let (decision_day_text, identity_and_later) =
        decision_and_later.split_once(",\"identity\":")?;
    let now_day = now_day_text.parse::<u64>().ok()?;
    let decision_day = decision_day_text.parse::<u64>().ok()?;
    if now_day.to_string() != now_day_text
        || decision_day.to_string() != decision_day_text
        || now_day != decision_day
    {
        return None;
    }
    let (identity, pre_and_later) = identity_and_later.split_once(",\"pre\":")?;
    if !identity.starts_with('{') || !identity.ends_with('}') {
        return None;
    }
    let (pre, post_and_later) = pre_and_later.split_once(",\"post\":")?;
    let (post, baseline_and_later) = post_and_later.split_once(",\"baseline_hash\":")?;
    if !pre.starts_with('{')
        || !pre.ends_with('}')
        || !post.starts_with('{')
        || !post.ends_with('}')
    {
        return None;
    }
    let (baseline_hash, later) = baseline_and_later.split_once(",\"baseline\":")?;
    let baseline_hash = baseline_hash.strip_prefix('"')?.strip_suffix('"')?;
    let parsed_baseline_hash = fs_blake3::ContentHash::from_hex(baseline_hash)?;
    if parsed_baseline_hash.to_hex() != baseline_hash {
        return None;
    }
    let (baseline, attestation_and_sources) = later.split_once(",\"attestation\":")?;
    if !baseline.starts_with('{') || !baseline.ends_with('}') {
        return None;
    }
    let (attestation, sources_and_authority) =
        attestation_and_sources.split_once(",\"required_source_receipts\":")?;
    let attestation_body = attestation
        .strip_prefix("{\"key_id\":\"")?
        .strip_suffix("\"}")?;
    let (key_id, signature) = attestation_body.split_once("\",\"signature\":\"")?;
    if key_id.is_empty() || signature.is_empty() {
        return None;
    }
    let (required_sources, authority_and_verdict) =
        sources_and_authority.split_once(",\"authority\":")?;
    let source_body = required_sources.strip_prefix('[')?.strip_suffix(']')?;
    let mut previous_source = None;
    let mut parsed_sources = Vec::new();
    for encoded in source_body.split(',') {
        let hex = encoded.strip_prefix('"')?.strip_suffix('"')?;
        let source = fs_blake3::ContentHash::from_hex(hex)?;
        if source.to_string() != hex || previous_source.is_some_and(|previous| previous >= source) {
            return None;
        }
        previous_source = Some(source);
        parsed_sources.push(source);
    }
    if parsed_sources.len() < baseline::MIN_PROMOTION_RUNS {
        return None;
    }
    let (authority, verdict) = authority_and_verdict.split_once(",\"verdict\":")?;
    let policy_receipt = authority
        .strip_prefix("{\"verdict\":\"authorized\",\"policy_receipt\":\"")?
        .strip_suffix("\"}")?;
    let parsed_policy_receipt = fs_blake3::ContentHash::from_hex(policy_receipt)?;
    if parsed_policy_receipt.to_hex() != policy_receipt {
        return None;
    }
    if verdict != "{\"baseline\":\"trusted\"}}" {
        return None;
    }
    Some(BaselineReceiptView {
        now_day,
        decision_day,
        identity,
        pre,
        post,
        baseline_hash: parsed_baseline_hash,
        baseline,
        required_sources: parsed_sources,
    })
}

fn parse_machine_axes_receipt(text: &str) -> Option<(u64, u64, [u64; 4])> {
    fn take(rest: &mut &str, prefix: &str) -> Option<()> {
        *rest = rest.strip_prefix(prefix)?;
        Some(())
    }
    fn decimal(rest: &mut &str) -> Option<u64> {
        let end = rest
            .find(|character: char| !character.is_ascii_digit())
            .unwrap_or(rest.len());
        if end == 0 || (end > 1 && rest.starts_with('0')) {
            return None;
        }
        let (digits, tail) = rest.split_at(end);
        *rest = tail;
        digits.parse().ok()
    }
    fn hex_u64(rest: &mut &str) -> Option<u64> {
        let (raw, tail) = rest.split_once('"')?;
        if raw.len() != 16
            || !raw
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return None;
        }
        *rest = tail;
        u64::from_str_radix(raw, 16).ok()
    }

    let mut rest = text;
    take(&mut rest, "{\"fingerprint\":\"")?;
    let fingerprint = hex_u64(&mut rest)?;
    take(&mut rest, ",\"cpu_brand\":\"")?;
    let (_, tail) = rest.split_once("\",\"logical_cpus\":")?;
    rest = tail;
    let logical_cpus = decimal(&mut rest)?;
    if logical_cpus == 0 {
        return None;
    }
    take(&mut rest, ",\"bandwidth_single_bits\":\"")?;
    let bandwidth_single = hex_u64(&mut rest)?;
    take(&mut rest, ",\"bandwidth_all_core_bits\":\"")?;
    let bandwidth_all_core = hex_u64(&mut rest)?;
    take(&mut rest, ",\"peak_single_bits\":\"")?;
    let peak_single = hex_u64(&mut rest)?;
    take(&mut rest, ",\"peak_all_core_bits\":\"")?;
    let peak_all_core = hex_u64(&mut rest)?;
    take(&mut rest, "}")?;
    if !rest.is_empty() {
        return None;
    }
    Some((
        fingerprint,
        logical_cpus,
        [
            bandwidth_single,
            bandwidth_all_core,
            peak_single,
            peak_all_core,
        ],
    ))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ValidatedProtocolAxes {
    pre_logical_cpus: u64,
    pre_axis_bits: [u64; 4],
    post_axis_bits: [u64; 4],
    decision_day: u64,
}

fn validate_protocol_axes(
    protocol: &CanonicalProductionOp,
    current_fingerprint: u64,
    current_baseline: fs_blake3::ContentHash,
) -> Option<ValidatedProtocolAxes> {
    let baseline = parse_baseline_receipt_view(&protocol.baseline_admission)?;
    let (pre_fingerprint, pre_logical_cpus, pre_axis_bits) =
        parse_machine_axes_receipt(baseline.pre)?;
    let (post_fingerprint, _, post_axis_bits) = parse_machine_axes_receipt(baseline.post)?;
    let parsed_baseline = BaselineStore::from_jsonl(baseline.baseline).ok()?;
    let exact_baseline = parsed_baseline.for_fingerprint(pre_fingerprint)?;
    let exact_identity = baseline_identity_json(exact_baseline.identity());
    if protocol.fingerprint != current_fingerprint
        || pre_fingerprint != current_fingerprint
        || protocol.post_fingerprint != post_fingerprint
        || baseline.baseline_hash != current_baseline
        || exact_baseline.content_hash() != baseline.baseline_hash
        || exact_baseline.canonical_json() != baseline.baseline
        || baseline.identity != exact_identity.as_str()
        || baseline.required_sources.as_slice() != exact_baseline.provenance().source_receipts()
        || baseline.now_day != baseline.decision_day
        || fs_blake3::hash_domain(
            production::PRODUCTION_AXES_RECEIPT_DOMAIN,
            baseline.pre.as_bytes(),
        ) != protocol.pre_axes_receipt
        || fs_blake3::hash_domain(
            production::PRODUCTION_AXES_RECEIPT_DOMAIN,
            baseline.post.as_bytes(),
        ) != protocol.post_axes_receipt
    {
        return None;
    }
    Some(ValidatedProtocolAxes {
        pre_logical_cpus,
        pre_axis_bits,
        post_axis_bits,
        decision_day: baseline.decision_day,
    })
}

const WALL_NS_PER_DAY: u64 = 86_400 * 1_000_000_000;

fn wall_ns_day(wall_ns: i64) -> Option<u64> {
    Some(u64::try_from(wall_ns).ok()? / WALL_NS_PER_DAY)
}

fn artifact_bytes_for_validation(
    ledger: &Ledger,
    artifact: &fs_blake3::ContentHash,
    max_bytes: u64,
) -> Result<Option<Vec<u8>>, LedgerError> {
    match ledger.get_artifact_bounded(artifact, max_bytes) {
        Ok(bytes) => Ok(bytes),
        Err(LedgerError::Corrupt { .. } | LedgerError::ArtifactReadLimit { .. }) => Ok(None),
        Err(error) => Err(error),
    }
}

fn artifact_envelope_is_valid(
    ledger: &Ledger,
    artifact: &fs_blake3::ContentHash,
    expected_kind: &str,
    expected_meta: &str,
) -> Result<bool, LedgerError> {
    let info = match ledger.artifact_info(artifact) {
        Ok(info) => info,
        Err(LedgerError::Corrupt { .. }) => return Ok(false),
        Err(error) => return Err(error),
    };
    Ok(info.is_some_and(|info| {
        info.kind == expected_kind && info.meta.as_deref() == Some(expected_meta)
    }))
}

fn dependency_receipt_is_structurally_valid(
    ledger: &Ledger,
    op: i64,
    protocol: &CanonicalProductionOp,
) -> Result<bool, LedgerError> {
    if !artifact_envelope_is_valid(
        ledger,
        &protocol.dependency_receipt_artifact,
        DEPGRAPH_RECEIPT_ARTIFACT_KIND,
        DEPGRAPH_RECEIPT_ARTIFACT_META,
    )? {
        return Ok(false);
    }
    let Some(bytes) = artifact_bytes_for_validation(
        ledger,
        &protocol.dependency_receipt_artifact,
        MAX_DEPGRAPH_RECEIPT_BYTES,
    )?
    else {
        return Ok(false);
    };
    Ok(
        fs_blake3::hash_domain(fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN, &bytes)
            == protocol.dependency_receipt_digest
            && ledger.edge_exists(op, &protocol.dependency_receipt_artifact, EdgeRole::In)?,
    )
}

fn dependency_receipt_matches_binding(
    protocol: &CanonicalProductionOp,
    expected: Option<DependencyReceiptBinding>,
) -> bool {
    expected.is_some_and(|expected| {
        protocol.dependency_receipt_artifact == expected.artifact_hash
            && protocol.dependency_receipt_digest == expected.domain_digest
    })
}

fn production_op_envelope_is_valid(
    op: &fs_ledger::OpRow,
    build_identity: fs_blake3::ContentHash,
) -> bool {
    op.session.as_deref() == Some(b"roofline".as_slice())
        && op.seed == ROOFLINE_SEED
        && op.versions == versions_json(build_identity)
        && op.budget == ROOFLINE_BUDGET
        && op.capability == ROOFLINE_CAPABILITY
        && op.outcome.as_deref() == Some("ok")
        && op.diag.is_none()
        && op.t_end.is_some_and(|end| end >= op.t_start)
}

fn validate_roofline_row(
    ledger: &Ledger,
    row: &fs_ledger::TuneRow,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: fs_blake3::ContentHash,
    expected_dependency: Option<DependencyReceiptBinding>,
) -> Result<Option<ValidatedRooflineRow>, LedgerError> {
    let Some(params) = parse_roofline_row_params(&row.params) else {
        return Ok(None);
    };
    let op = match ledger.op(params.op) {
        Ok(Some(op)) => op,
        Ok(None) | Err(LedgerError::OpCorrupt { .. }) => return Ok(None),
        Err(error) => return Err(error),
    };
    let Some(recorded_at_ns) = op.t_end else {
        return Ok(None);
    };
    let Some(protocol) = parse_canonical_production_op(&op.ir) else {
        return Ok(None);
    };
    let Some(validated_axes) =
        validate_protocol_axes(&protocol, current_fingerprint, current_baseline)
    else {
        return Ok(None);
    };
    if wall_ns_day(recorded_at_ns) != Some(validated_axes.decision_day)
        || !dependency_receipt_is_structurally_valid(ledger, params.op, &protocol)?
        || op.id != params.op
        || !production_op_envelope_is_valid(&op, params.build_identity)
        || protocol.finalized_run_receipt != params.run_receipt
    {
        return Ok(None);
    }
    let machine_key = roofline_machine_key(current_fingerprint, current_baseline);
    if !stored_manifest_member_is_valid(
        ledger,
        row,
        kernel,
        version,
        &protocol,
        &params,
        machine_key,
        params.op,
        current_baseline,
        params.build_identity,
        validated_axes.pre_logical_cpus,
        validated_axes.pre_axis_bits,
        validated_axes.post_axis_bits,
    )? {
        return Ok(None);
    }

    if !receipt_recomputes_from_stored_rows(
        ledger,
        &protocol,
        &params,
        kernel,
        version,
        machine_key,
        validated_axes.pre_logical_cpus,
        validated_axes.pre_axis_bits,
        validated_axes.post_axis_bits,
    )? {
        return Ok(None);
    }

    Ok(Some(ValidatedRooflineRow {
        build_identity: params.build_identity,
        recorded_at_ns,
        dependency_matches_current: dependency_receipt_matches_binding(
            &protocol,
            expected_dependency,
        ),
        dependency_receipt_digest: protocol.dependency_receipt_digest,
        dependency_receipt_artifact: protocol.dependency_receipt_artifact,
    }))
}

#[allow(clippy::too_many_arguments)]
fn stored_manifest_member_is_valid(
    ledger: &Ledger,
    row: &fs_ledger::TuneRow,
    kernel: &str,
    version: &str,
    protocol: &CanonicalProductionOp,
    params: &RooflineRowParams,
    machine_key: [u8; 40],
    op_id: i64,
    baseline_hash: fs_blake3::ContentHash,
    build_identity: fs_blake3::ContentHash,
    pre_logical_cpus: u64,
    pre_axis_bits: [u64; 4],
    post_axis_bits: [u64; 4],
) -> Result<bool, LedgerError> {
    if row.kernel != kernel
        || params.op != op_id
        || params.run_receipt != protocol.finalized_run_receipt
        || params.dependency_receipt_artifact != Some(protocol.dependency_receipt_artifact)
        || params.dependency_receipt_digest != Some(protocol.dependency_receipt_digest)
        || params.baseline_hash != baseline_hash
        || params.build_identity != build_identity
        || params.post_axis_bits != post_axis_bits
        || row.machine != machine_key
        || row.shape_class
            != tune_measurement_shape_class(version, protocol.finalized_run_receipt, op_id)
        || params.payload_artifact != fs_ledger::hash_bytes(row.measured.as_bytes())
        || !artifact_envelope_is_valid(
            ledger,
            &params.payload_artifact,
            ROOFLINE_PAYLOAD_ARTIFACT_KIND,
            ROOFLINE_PAYLOAD_ARTIFACT_META,
        )?
    {
        return Ok(false);
    }
    let Some(artifact_bytes) = artifact_bytes_for_validation(
        ledger,
        &params.payload_artifact,
        u64::try_from(row.measured.len()).expect("a Rust slice length fits u64"),
    )?
    else {
        return Ok(false);
    };
    if artifact_bytes.as_slice() != row.measured.as_bytes()
        || !ledger.edge_exists(params.op, &params.payload_artifact, EdgeRole::Out)?
    {
        return Ok(false);
    }
    let measured_prefix = format!(
        "{{\"receipt_version\":3,\"kernel\":\"{}\",\"version\":\"{}\",\"machine\":\"{:016x}\",\"axes\":{{\"logical_cpus\":{pre_logical_cpus},\"bandwidth_single_bits\":\"{:016x}\",\"bandwidth_all_core_bits\":\"{:016x}\",\"peak_single_bits\":\"{:016x}\",\"peak_all_core_bits\":\"{:016x}\"}},",
        json_escape(kernel),
        json_escape(version),
        protocol.fingerprint,
        pre_axis_bits[0],
        pre_axis_bits[1],
        pre_axis_bits[2],
        pre_axis_bits[3],
    );
    let Some((_, reps_and_later)) = row.measured.split_once("\"reps\":") else {
        return Ok(false);
    };
    let Some((reps, _)) = reps_and_later.split_once(",\"verdict\":") else {
        return Ok(false);
    };
    Ok(row.measured.starts_with(&measured_prefix)
        && row.measured.matches("\"reps\":").count() == 1
        && reps.parse::<u64>().ok() == Some(params.reps))
}

/// Reconstruct the finalized run receipt from the operation-bound ordered
/// result manifest and the rows actually stored today (bead gp3.15). The
/// manifest lives in the op's `ir`, which no ledger API mutates after
/// `begin_op`; the receipt binds baseline receipt bytes, ordered payload
/// bytes, and the manifest hash. A writer who replaces one payload plus its
/// matching artifact/params while retaining the old run receipt now fails
/// this recomputation instead of classifying as fresh.
#[allow(clippy::too_many_arguments)]
fn receipt_recomputes_from_stored_rows(
    ledger: &Ledger,
    protocol: &CanonicalProductionOp,
    params: &RooflineRowParams,
    kernel: &str,
    version: &str,
    machine_key: [u8; 40],
    pre_logical_cpus: u64,
    pre_axis_bits: [u64; 4],
    post_axis_bits: [u64; 4],
) -> Result<bool, LedgerError> {
    let Some(entries) = parse_result_manifest(&protocol.result_manifest) else {
        return Ok(false);
    };
    if entries.is_empty()
        || !entries.iter().any(|e| {
            e.kernel == kernel && e.version == version && e.payload == params.payload_artifact
        })
    {
        return Ok(false);
    }
    let mut receipt_payload = Vec::new();
    push_receipt_field(&mut receipt_payload, protocol.baseline_admission.as_bytes());
    let entry_count = u64::try_from(entries.len()).expect("manifest entry count fits u64");
    receipt_payload.extend_from_slice(&entry_count.to_le_bytes());
    for entry in &entries {
        let shape =
            tune_measurement_shape_class(&entry.version, protocol.finalized_run_receipt, params.op);
        let Some(stored) = ledger.tune_get(&entry.kernel, &shape, &machine_key)? else {
            return Ok(false);
        };
        let Some(stored_params) = parse_roofline_row_params(&stored.params) else {
            return Ok(false);
        };
        if stored_params.op != params.op
            || stored_params.baseline_hash != params.baseline_hash
            || !stored_manifest_member_is_valid(
                ledger,
                &stored,
                &entry.kernel,
                &entry.version,
                protocol,
                &stored_params,
                machine_key,
                params.op,
                params.baseline_hash,
                params.build_identity,
                pre_logical_cpus,
                pre_axis_bits,
                post_axis_bits,
            )?
            || stored_params.payload_artifact != entry.payload
        {
            return Ok(false);
        }
        push_receipt_field(&mut receipt_payload, stored.measured.as_bytes());
    }
    let manifest_hash =
        fs_blake3::hash_domain(RESULT_MANIFEST_DOMAIN, protocol.result_manifest.as_bytes());
    push_receipt_field(&mut receipt_payload, manifest_hash.as_bytes());
    Ok(
        fs_blake3::hash_domain(FINALIZED_RUN_DOMAIN, &receipt_payload)
            == protocol.finalized_run_receipt,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    struct ReceiptKernel {
        elements: usize,
        value: u64,
    }

    impl RooflineKernel for ReceiptKernel {
        fn spec(&self) -> KernelSpec {
            KernelSpec {
                name: "receipt-kernel",
                version: "1",
                bytes_per_elem: 8.0,
                flops_per_elem: 1.0,
                threading: Threading::SingleThread,
                target_axis: TargetAxis::BindingRoof,
                target_fraction: None,
            }
        }

        fn elements(&self) -> usize {
            self.elements
        }

        fn run_once(&mut self) -> Result<(), String> {
            for _ in 0..1024 {
                self.value = std::hint::black_box(
                    self.value
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1),
                );
            }
            Ok(())
        }
    }

    struct FailingRunKernel {
        calls: usize,
        fail_at: usize,
    }

    impl RooflineKernel for FailingRunKernel {
        fn spec(&self) -> KernelSpec {
            ReceiptKernel {
                elements: 1,
                value: 0,
            }
            .spec()
        }

        fn elements(&self) -> usize {
            1
        }

        fn run_once(&mut self) -> Result<(), String> {
            let invocation = self.calls;
            self.calls = self
                .calls
                .checked_add(1)
                .ok_or_else(|| "test kernel invocation counter overflowed".to_string())?;
            if invocation == self.fail_at {
                return Err("injected kernel refusal".to_string());
            }
            Ok(())
        }
    }

    struct AbortProbeKernel {
        name: &'static str,
        fail: bool,
        pending: std::rc::Rc<std::cell::Cell<bool>>,
        aborts: std::rc::Rc<std::cell::Cell<usize>>,
    }

    impl RooflineKernel for AbortProbeKernel {
        fn spec(&self) -> KernelSpec {
            let mut spec = ReceiptKernel {
                elements: 1,
                value: 0,
            }
            .spec();
            spec.name = self.name;
            spec
        }

        fn elements(&self) -> usize {
            1
        }

        fn run_once(&mut self) -> Result<(), String> {
            if self.fail {
                return Err("injected registry peer refusal".to_string());
            }
            self.pending.set(true);
            Ok(())
        }

        fn abort_tuning(&mut self) -> Result<(), String> {
            self.pending.set(false);
            self.aborts.set(self.aborts.get() + 1);
            Ok(())
        }
    }

    struct AdmissionProbeKernel {
        observed: std::rc::Rc<std::cell::Cell<Option<bool>>>,
    }

    impl RooflineKernel for AdmissionProbeKernel {
        fn spec(&self) -> KernelSpec {
            ReceiptKernel {
                elements: 1,
                value: 0,
            }
            .spec()
        }

        fn elements(&self) -> usize {
            1
        }

        fn run_once(&mut self) -> Result<(), String> {
            Ok(())
        }

        fn finalize_tuning(&mut self, admitted: bool) -> Result<(), String> {
            self.observed.set(Some(admitted));
            Ok(())
        }
    }

    struct FallibleAdmissionProbeKernel {
        id: usize,
        observed: std::rc::Rc<std::cell::RefCell<Vec<(usize, bool)>>>,
        aborted: std::rc::Rc<std::cell::RefCell<Vec<usize>>>,
        failure: Option<&'static str>,
    }

    impl RooflineKernel for FallibleAdmissionProbeKernel {
        fn spec(&self) -> KernelSpec {
            let mut spec = ReceiptKernel {
                elements: 1,
                value: 0,
            }
            .spec();
            spec.name = ["receipt-kernel-0", "receipt-kernel-1", "receipt-kernel-2"]
                .get(self.id)
                .copied()
                .unwrap_or("receipt-kernel-out-of-range");
            spec
        }

        fn elements(&self) -> usize {
            1
        }

        fn run_once(&mut self) -> Result<(), String> {
            let mut value = self.id as u64;
            for _ in 0..1024 {
                value = std::hint::black_box(
                    value
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1),
                );
            }
            std::hint::black_box(value);
            Ok(())
        }

        fn finalize_tuning(&mut self, admitted: bool) -> Result<(), String> {
            self.observed.borrow_mut().push((self.id, admitted));
            self.failure.map_or(Ok(()), |error| Err(error.to_string()))
        }

        fn abort_tuning(&mut self) -> Result<(), String> {
            self.aborted.borrow_mut().push(self.id);
            Ok(())
        }
    }

    fn synthetic_axes() -> MachineAxes {
        MachineAxes {
            fingerprint: 0xABCD,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 100.0,
            bandwidth_all_core_gbs: 400.0,
            peak_single_gflops: 50.0,
            peak_all_core_gflops: 300.0,
        }
    }

    fn trusted_baseline(axes: &MachineAxes) -> (BaselineAxes, BaselineIdentity) {
        trusted_baseline_at(axes, 20_000)
    }

    fn live_trusted_baseline(axes: &MachineAxes) -> (BaselineAxes, BaselineIdentity) {
        trusted_baseline_at(axes, test_day().saturating_sub(10))
    }

    fn trusted_baseline_at(
        axes: &MachineAxes,
        promoted_day: u64,
    ) -> (BaselineAxes, BaselineIdentity) {
        let identity =
            BaselineIdentity::current(axes, "test-firmware").expect("valid synthetic identity");
        let candidates: Vec<_> = (0_u64..3)
            .map(|ordinal| {
                BaselineCandidate::from_receipt(
                    axes.clone(),
                    identity.clone(),
                    fs_blake3::hash_domain(
                        "fs-roofline.lib-test-baseline-source.v1",
                        &ordinal.to_le_bytes(),
                    ),
                )
                .expect("valid synthetic candidate")
            })
            .collect();
        let baseline = promote_baseline(
            &candidates,
            "test-operator",
            "deterministic lib receipt fixture",
            promoted_day,
            90,
        )
        .expect("valid synthetic baseline");
        (baseline, identity)
    }

    fn test_day() -> u64 {
        days_since_epoch_now().expect("unit-test clock after Unix epoch")
    }

    fn attested_policy(
        baseline: &BaselineAxes,
        identity: &BaselineIdentity,
        now_day: u64,
    ) -> AttestedAxisBaselinePolicy {
        AttestedAxisBaselinePolicy::from_verified(
            baseline.clone(),
            identity.clone(),
            now_day,
            PromotionAttestation::new("test-authority", "test-signature"),
            baseline.provenance().source_receipts().to_vec(),
            PromotionAuthorityDecision::new(
                KeyVerdict::Authorized,
                fs_blake3::hash_domain(
                    "fs-roofline.lib-test-policy.v1",
                    baseline.content_hash().as_bytes(),
                ),
            ),
        )
    }

    fn timed_receipt_fixture() -> (MachineAxes, BaselineAxes, BaselineIdentity, Attainment) {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let mut kernel = ReceiptKernel {
            elements: 1,
            value: 0,
        };
        let timed = measure(&mut kernel, 0, 3, &axes).expect("bounded receipt measurement");
        (axes, baseline, identity, timed)
    }

    #[test]
    fn attested_policy_and_snapshot_are_scoped_to_one_epoch_day() {
        let axes = synthetic_axes();
        let (baseline, identity) = live_trusted_baseline(&axes);
        let today = test_day();
        let admitted = attested_policy(&baseline, &identity, today).decide_at(&axes, &axes, today);
        assert!(admitted.authority_admitted());
        assert!(admitted.baseline_citation_eligible());
        assert_eq!(admitted.decision_day(), Some(today));
        assert!(admitted.day_error_at(today).is_none());
        assert!(admitted.day_error_at(today.saturating_add(1)).is_some());

        let delayed = attested_policy(&baseline, &identity, today).decide_at(
            &axes,
            &axes,
            today.saturating_add(1),
        );
        assert!(!delayed.authority_admitted());
        assert!(!delayed.baseline_citation_eligible());
        assert_eq!(
            delayed.verdict(),
            &BaselineVerdict::Unauthorized {
                verdict: "policy-day-mismatch"
            }
        );
        assert!(delayed.receipt_json().contains(&format!(
            "\"now_day\":{today},\"decision_day\":{}",
            today.saturating_add(1)
        )));
    }

    #[test]
    fn attested_baseline_receipt_requires_lowercase_child_hashes() {
        let axes = synthetic_axes();
        let today = test_day();
        let (baseline, identity) = trusted_baseline_at(&axes, today.saturating_sub(10));
        let policy_receipt = fs_blake3::hash_domain(
            "fs-roofline.lib-test-canonical-policy.v1",
            baseline.content_hash().as_bytes(),
        );
        let snapshot = AttestedAxisBaselinePolicy::from_verified(
            baseline.clone(),
            identity,
            today,
            PromotionAttestation::new("test-authority", "test-signature"),
            baseline.provenance().source_receipts().to_vec(),
            PromotionAuthorityDecision::new(KeyVerdict::Authorized, policy_receipt),
        )
        .decide_at(&axes, &axes, today);
        let canonical = snapshot.receipt_json();
        assert!(parse_baseline_receipt_view(canonical).is_some());

        let baseline_hex = baseline.content_hash().to_hex();
        let uppercase_baseline = canonical.replacen(
            &format!("\"baseline_hash\":\"{baseline_hex}\""),
            &format!(
                "\"baseline_hash\":\"{}\"",
                baseline_hex.to_ascii_uppercase()
            ),
            1,
        );
        assert_ne!(uppercase_baseline, canonical);
        assert!(parse_baseline_receipt_view(&uppercase_baseline).is_none());

        let policy_hex = policy_receipt.to_hex();
        let uppercase_policy = canonical.replacen(
            &format!("\"policy_receipt\":\"{policy_hex}\""),
            &format!("\"policy_receipt\":\"{}\"", policy_hex.to_ascii_uppercase()),
            1,
        );
        assert_ne!(uppercase_policy, canonical);
        assert!(parse_baseline_receipt_view(&uppercase_policy).is_none());
    }

    fn external_gate_snapshot(
        axes: &MachineAxes,
        key_id: &str,
    ) -> (BaselineAxes, AxisAdmissionSnapshot, fs_blake3::ContentHash) {
        for _ in 0..3 {
            let today = test_day();
            let (baseline, identity) = trusted_baseline_at(axes, today.saturating_sub(10));
            let policy_receipt = fs_blake3::hash_domain(
                "fs-roofline.lib-test-external-gate-policy.v1",
                key_id.as_bytes(),
            );
            let snapshot = AttestedAxisBaselinePolicy::from_verified(
                baseline.clone(),
                identity,
                today,
                PromotionAttestation::new(key_id, "external-gate-test-signature"),
                baseline.provenance().source_receipts().to_vec(),
                PromotionAuthorityDecision::new(KeyVerdict::Authorized, policy_receipt),
            )
            .decide_at(axes, axes, today);
            if test_day() == today && snapshot.baseline_citation_eligible() {
                return (baseline, snapshot, policy_receipt);
            }
        }
        panic!("could not construct an external-gate fixture within one stable epoch day");
    }

    fn external_gate_json(lane: ExternalPerfGateLane, admission: &AxisAdmissionSnapshot) -> String {
        format!(
            "{{\"metric\":\"{}\",\"target_met\":true,\"citation_eligible\":true,\"recorded\":true,\"report_only\":false,\"reason\":null,\"admission\":{}}}",
            lane.metric(),
            admission.receipt_json(),
        )
    }

    #[test]
    fn external_gate_recorder_preserves_exact_authority_and_gate_receipts() {
        let axes = synthetic_axes();
        let key_id = "rotated-external-gate-key";
        let (baseline, snapshot, policy_receipt) = external_gate_snapshot(&axes, key_id);
        let gate_json = external_gate_json(ExternalPerfGateLane::Feec, &snapshot);
        let expected_admission = snapshot.receipt_json().as_bytes().to_vec();
        let expected_gate = gate_json.as_bytes().to_vec();
        let ledger = Ledger::open(":memory:").expect("in-memory external-gate ledger");

        let receipt = record_external_perf_gate_in_ledger(
            &ledger,
            ExternalPerfGateLane::Feec,
            &snapshot,
            &gate_json,
        )
        .expect("authority-admitted external gate must record atomically");

        assert_eq!(receipt.lane(), ExternalPerfGateLane::Feec);
        assert!(receipt.event_id() > 0);
        assert_eq!(ledger.table_count("events").expect("count events"), 1);
        assert_eq!(
            receipt.admission_artifact(),
            fs_ledger::hash_bytes(&expected_admission)
        );
        assert_eq!(
            receipt.final_gate_artifact(),
            fs_ledger::hash_bytes(&expected_gate)
        );
        assert_eq!(
            ledger
                .get_artifact_bounded(
                    &receipt.admission_artifact(),
                    u64::try_from(expected_admission.len()).expect("bounded fixture"),
                )
                .expect("read admission artifact")
                .expect("admission artifact retained"),
            expected_admission
        );
        assert_eq!(
            ledger
                .get_artifact_bounded(
                    &receipt.final_gate_artifact(),
                    u64::try_from(expected_gate.len()).expect("bounded fixture"),
                )
                .expect("read gate artifact")
                .expect("gate artifact retained"),
            expected_gate
        );
        assert!(
            ledger
                .edge_exists(
                    receipt.operation_id(),
                    &receipt.admission_artifact(),
                    EdgeRole::In,
                )
                .expect("admission edge query")
        );
        assert!(
            ledger
                .edge_exists(
                    receipt.operation_id(),
                    &receipt.final_gate_artifact(),
                    EdgeRole::Out,
                )
                .expect("gate edge query")
        );
        let stored_op = ledger
            .op(receipt.operation_id())
            .expect("operation query")
            .expect("operation retained");
        assert!(
            stored_op
                .ir
                .contains(&format!("\"final_gate\":{gate_json}"))
        );
        assert!(receipt.receipt_json().contains("\"lane\":\"feec\""));
        assert!(
            receipt
                .receipt_json()
                .contains(&format!("\"event\":{}", receipt.event_id()))
        );

        let stored_admission = String::from_utf8(
            ledger
                .get_artifact_bounded(
                    &receipt.admission_artifact(),
                    u64::try_from(snapshot.receipt_json().len()).expect("bounded fixture"),
                )
                .expect("read admission artifact")
                .expect("admission artifact retained"),
        )
        .expect("canonical admission UTF-8");
        assert!(stored_admission.contains(&format!("\"key_id\":\"{key_id}\"")));
        assert!(stored_admission.contains("\"signature\":\"external-gate-test-signature\""));
        assert!(stored_admission.contains(&format!("\"policy_receipt\":\"{policy_receipt}\"")));
        for source in baseline.provenance().source_receipts() {
            assert!(
                stored_admission.contains(&format!("\"{source}\"")),
                "source receipt {source} must survive byte-for-byte"
            );
        }
    }

    #[test]
    fn external_gate_recorder_refuses_candidate_and_leaves_no_rows() {
        let axes = synthetic_axes();
        let today = test_day();
        let (baseline, identity) = trusted_baseline_at(&axes, today.saturating_sub(10));
        let candidate =
            AxisBaselinePolicy::new(Some(&baseline), &identity, today).snapshot(&axes, &axes);
        assert!(!candidate.authority_admitted());
        let gate_json = external_gate_json(ExternalPerfGateLane::Fft, &candidate);
        let ledger = Ledger::open(":memory:").expect("in-memory external-gate ledger");

        let error = record_external_perf_gate_in_ledger(
            &ledger,
            ExternalPerfGateLane::Fft,
            &candidate,
            &gate_json,
        )
        .expect_err("candidate admission must stay report-only");
        assert!(
            error
                .to_string()
                .contains("lacks authorized promotion authority")
        );
        assert_eq!(ledger.table_count("ops").expect("count ops"), 0);
        assert_eq!(ledger.table_count("artifacts").expect("count artifacts"), 0);
        assert_eq!(ledger.table_count("edges").expect("count edges"), 0);
        assert_eq!(ledger.table_count("events").expect("count events"), 0);
        assert!(!ledger.in_transaction());
    }

    #[test]
    fn external_gate_recorder_rejects_nonpositive_or_mismatched_payloads() {
        let axes = synthetic_axes();
        let (_, snapshot, _) = external_gate_snapshot(&axes, "external-gate-key");
        let valid = external_gate_json(ExternalPerfGateLane::Feec, &snapshot);
        let wrong_admission = valid.replacen(snapshot.receipt_json(), "{}", 1);
        let duplicate_citation = valid.replacen('{', "{\"citation_eligible\":true,", 1);
        let trailing_comma = format!(
            "{},}}",
            valid
                .strip_suffix('}')
                .expect("fixture is a top-level object")
        );
        let deeply_nested = format!(
            "{{\"nested\":{}null{},{}",
            "[".repeat(MAX_EXTERNAL_PERF_GATE_NESTING + 1),
            "]".repeat(MAX_EXTERNAL_PERF_GATE_NESTING + 1),
            &valid[1..],
        );
        let invalid = [
            valid.replacen("\"feec-gate\"", "\"fft-gate\"", 1),
            valid.replacen(
                "\"citation_eligible\":true",
                "\"citation_eligible\":false",
                1,
            ),
            valid.replacen("\"recorded\":true", "\"recorded\":false", 1),
            valid.replacen("\"report_only\":false", "\"report_only\":true", 1),
            wrong_admission,
            duplicate_citation,
            trailing_comma,
            deeply_nested,
            format!("{valid} false"),
        ];
        let ledger = Ledger::open(":memory:").expect("in-memory external-gate ledger");
        for payload in invalid {
            record_external_perf_gate_in_ledger(
                &ledger,
                ExternalPerfGateLane::Feec,
                &snapshot,
                &payload,
            )
            .expect_err("nonpositive or mismatched external gate must be refused");
        }
        assert_eq!(ledger.table_count("ops").expect("count ops"), 0);
        assert_eq!(ledger.table_count("artifacts").expect("count artifacts"), 0);
    }

    #[test]
    fn external_gate_recorder_refuses_caller_transactions_and_memory_paths() {
        let axes = synthetic_axes();
        let (_, snapshot, _) = external_gate_snapshot(&axes, "external-gate-key");
        let gate_json = external_gate_json(ExternalPerfGateLane::Fft, &snapshot);
        let ledger = Ledger::open(":memory:").expect("in-memory external-gate ledger");
        ledger.begin().expect("caller-owned transaction");
        record_external_perf_gate_in_ledger(
            &ledger,
            ExternalPerfGateLane::Fft,
            &snapshot,
            &gate_json,
        )
        .expect_err("recorder must not compose with a caller transaction");
        assert!(ledger.in_transaction(), "caller transaction remains owned");
        ledger.rollback().expect("caller rollback");

        for path in [
            "",
            "   ",
            ":memory:",
            "  :memory:  ",
            "file::memory:",
            "file::memory:?cache=shared",
            "file:gate?mode=memory&cache=shared",
            "FILE:gate?MODE=MEMORY",
            "file:persistent.db",
        ] {
            record_external_perf_gate_at_path(
                path,
                ExternalPerfGateLane::Fft,
                &snapshot,
                &gate_json,
            )
            .expect_err("positive evidence requires a durable ledger path");
        }
    }

    #[test]
    fn production_protocol_version_and_field_are_locked_together() {
        assert_eq!(
            PRODUCTION_PROTOCOL_FIELD,
            format!("\"protocol\":\"{PRODUCTION_PROTOCOL_VERSION}\"")
        );
        assert_eq!(
            CUSTOM_REGISTRY_PROTOCOL_FIELD,
            format!("\"protocol\":\"{CUSTOM_REGISTRY_PROTOCOL_VERSION}\"")
        );
    }

    #[test]
    fn execution_binding_identity_versions_fail_closed() {
        assert_eq!(EXECUTION_BINDING_IDENTITY_VERSION, 4);
        assert!(EXECUTION_BINDING_DOMAIN.ends_with(".v4"));
        assert_eq!(EXECUTION_BINDING_KIND, "gemm-v4");
        let current_json = format!("{{\"kind\":\"{EXECUTION_BINDING_KIND}\"}}");
        let current =
            execution_binding_receipt_with_domain(&current_json, EXECUTION_BINDING_DOMAIN);
        assert_ne!(
            current,
            execution_binding_receipt_with_domain(
                &current_json,
                "org.frankensim.fs-roofline.execution-binding.v3",
            ),
            "the old domain cannot alias the complete v4 binding"
        );
        assert_ne!(
            current,
            execution_binding_receipt_with_domain(
                "{\"kind\":\"gemm-v3\"}",
                EXECUTION_BINDING_DOMAIN,
            ),
            "the retained kind tag is versioned independently of formatting"
        );
    }

    #[test]
    fn finalized_run_identity_fields_move_independently() {
        fn identity(input: &FinalizedRunIdentityInput) -> fs_blake3::ContentHash {
            finalized_run_receipt_from_input(input, FINALIZED_RUN_DOMAIN, RESULT_MANIFEST_DOMAIN)
        }

        let input = FinalizedRunIdentityInput {
            baseline_receipt: "{\"baseline\":\"trusted\"}".to_string(),
            result_payloads: vec!["row-a".to_string(), "row-b".to_string()],
            result_manifest: format!("{{\"schema\":\"{RESULT_MANIFEST_SCHEMA}\",\"entries\":[]}}"),
        };
        let original = identity(&input);

        assert_ne!(
            original,
            finalized_run_receipt_from_input(
                &input,
                "org.frankensim.fs-roofline.finalized-run-foreign.v3",
                RESULT_MANIFEST_DOMAIN,
            ),
            "the finalized-run digest domain is semantic"
        );
        assert_ne!(
            original,
            finalized_run_receipt_from_input(
                &input,
                FINALIZED_RUN_DOMAIN,
                "org.frankensim.fs-roofline.run-result-manifest-foreign.v1",
            ),
            "the result-manifest child domain is semantic"
        );

        let mut altered = input.clone();
        altered.baseline_receipt.push('x');
        assert_ne!(
            original,
            identity(&altered),
            "baseline receipt must move identity"
        );
        let mut altered = input.clone();
        altered.result_payloads.push("row-c".to_string());
        assert_ne!(
            original,
            identity(&altered),
            "result count must move identity"
        );
        let mut altered = input.clone();
        altered.result_payloads.swap(0, 1);
        assert_ne!(
            original,
            identity(&altered),
            "ordered result payloads must move identity"
        );
        let mut altered = input.clone();
        altered.result_manifest.push('x');
        assert_ne!(
            original,
            identity(&altered),
            "result-manifest bytes must move identity"
        );

        let mut unframed = Vec::new();
        unframed.extend_from_slice(input.baseline_receipt.as_bytes());
        for row in &input.result_payloads {
            unframed.extend_from_slice(row.as_bytes());
        }
        unframed.extend_from_slice(
            fs_blake3::hash_domain(RESULT_MANIFEST_DOMAIN, input.result_manifest.as_bytes())
                .as_bytes(),
        );
        assert_ne!(
            original,
            fs_blake3::hash_domain(FINALIZED_RUN_DOMAIN, &unframed),
            "removing the count and length prefixes must move identity"
        );
    }

    #[test]
    fn finalized_run_identity_versions_fail_closed() {
        assert_eq!(FINALIZED_RUN_IDENTITY_VERSION, 3);
        assert!(FINALIZED_RUN_DOMAIN.ends_with(".v3"));
        assert_eq!(RESULT_MANIFEST_IDENTITY_VERSION, 1);
        assert!(RESULT_MANIFEST_DOMAIN.ends_with(".v1"));
        let current = format!("{{\"schema\":\"{RESULT_MANIFEST_SCHEMA}\",\"entries\":[]}}");
        assert!(parse_result_manifest(&current).is_some());
        assert!(
            parse_result_manifest(&current.replace("manifest-v1", "manifest-v2")).is_none(),
            "a stale or future child-manifest version must fail closed"
        );
        let input = FinalizedRunIdentityInput {
            baseline_receipt: "baseline".to_string(),
            result_payloads: vec!["result".to_string()],
            result_manifest: current,
        };
        assert_ne!(
            finalized_run_receipt_from_input(&input, FINALIZED_RUN_DOMAIN, RESULT_MANIFEST_DOMAIN,),
            finalized_run_receipt_from_input(
                &input,
                "org.frankensim.fs-roofline.finalized-run.v4",
                RESULT_MANIFEST_DOMAIN,
            ),
            "a version/domain rotation must move the finalized-run identity"
        );
    }

    #[test]
    fn executable_build_identity_fields_move_independently() {
        let input = ExecutableBuildIdentityInput {
            byte_len: 3,
            raw_hash: fs_blake3::hash_bytes(b"abc"),
        };
        let original = executable_build_identity_from_input(&input, ROOFLINE_EXECUTABLE_DOMAIN);
        assert_ne!(
            original,
            executable_build_identity_from_input(
                &input,
                "org.frankensim.fs-roofline.executable-foreign.v1",
            ),
            "the executable digest domain is semantic"
        );
        let mut altered = input;
        altered.byte_len += 1;
        assert_ne!(
            original,
            executable_build_identity_from_input(&altered, ROOFLINE_EXECUTABLE_DOMAIN),
            "the executable byte count is semantic"
        );
        let mut altered = input;
        altered.raw_hash = fs_blake3::hash_bytes(b"abd");
        assert_ne!(
            original,
            executable_build_identity_from_input(&altered, ROOFLINE_EXECUTABLE_DOMAIN),
            "the raw executable content hash is semantic"
        );
        let mut reversed = [0_u8; 40];
        reversed[..32].copy_from_slice(input.raw_hash.as_bytes());
        reversed[32..].copy_from_slice(&input.byte_len.to_le_bytes());
        assert_ne!(
            original,
            fs_blake3::hash_domain(ROOFLINE_EXECUTABLE_DOMAIN, &reversed),
            "the length-prefix layout is semantic"
        );
    }

    #[test]
    fn executable_build_identity_excludes_path_and_chunking() {
        fn streamed_identity(_path: &str, chunks: &[&[u8]]) -> fs_blake3::ContentHash {
            let mut hasher = fs_blake3::Blake3::new();
            let mut byte_len = 0_u64;
            for chunk in chunks {
                byte_len += u64::try_from(chunk.len()).expect("fixture length fits u64");
                hasher.update(chunk);
            }
            executable_build_identity_from_input(
                &ExecutableBuildIdentityInput {
                    byte_len,
                    raw_hash: hasher.finalize(),
                },
                ROOFLINE_EXECUTABLE_DOMAIN,
            )
        }

        assert_eq!(
            streamed_identity("/first/path", &[b"abcdef"]),
            streamed_identity("/other/path", &[b"ab", b"c", b"def"]),
            "path spelling and read chunking are not executable content"
        );
    }

    #[test]
    fn executable_build_identity_versions_fail_closed() {
        assert_eq!(EXECUTABLE_BUILD_IDENTITY_VERSION, 1);
        assert!(ROOFLINE_EXECUTABLE_DOMAIN.ends_with(".v1"));
        let input = ExecutableBuildIdentityInput {
            byte_len: 3,
            raw_hash: fs_blake3::hash_bytes(b"abc"),
        };
        let current = executable_build_identity_from_input(&input, ROOFLINE_EXECUTABLE_DOMAIN);
        let stale = executable_build_identity_from_input(
            &input,
            "org.frankensim.fs-roofline.executable.v2",
        );
        assert!(
            require_stable_executable_identity(current, stale).is_err(),
            "a stale or future executable identity cannot satisfy the current build guard"
        );
    }

    #[test]
    fn public_measurement_resources_refuse_before_kernel_execution() {
        let axes = synthetic_axes();
        let mut kernel = ReceiptKernel {
            elements: 1,
            value: 0,
        };
        for (warmup, reps) in [
            (0, 0),
            (MAX_MEASUREMENT_WARMUP + 1, 1),
            (0, MAX_MEASUREMENT_REPS + 1),
            (usize::MAX, 1),
            (0, usize::MAX),
        ] {
            assert!(measure(&mut kernel, warmup, reps, &axes).is_err());
            assert_eq!(kernel.value, 0, "refused work must not execute the kernel");
        }

        let mut empty: Vec<Box<dyn RooflineKernel>> = Vec::new();
        assert!(run_registry(&mut empty, 0, 1, &axes).is_err());
        let mut oversized: Vec<Box<dyn RooflineKernel>> = (0..=MAX_REGISTRY_KERNELS)
            .map(|_| {
                Box::new(ReceiptKernel {
                    elements: 1,
                    value: 0,
                }) as Box<dyn RooflineKernel>
            })
            .collect();
        assert!(run_registry(&mut oversized, 0, 1, &axes).is_err());

        let mut invocation_heavy: Vec<Box<dyn RooflineKernel>> = (0..MAX_REGISTRY_KERNELS)
            .map(|_| {
                Box::new(ReceiptKernel {
                    elements: 1,
                    value: 0,
                }) as Box<dyn RooflineKernel>
            })
            .collect();
        assert!(
            run_registry(&mut invocation_heavy, MAX_MEASUREMENT_WARMUP, 1, &axes,).is_err(),
            "aggregate invocation admission must precede registry execution"
        );
    }

    #[test]
    fn kernel_failures_propagate_without_constructing_a_timing_row() {
        let axes = synthetic_axes();
        let mut warmup_failure = FailingRunKernel {
            calls: 0,
            fail_at: 0,
        };
        let error = measure(&mut warmup_failure, 1, 1, &axes)
            .expect_err("warmup refusal must fail the measurement");
        assert_eq!(
            error,
            "roofline kernel `receipt-kernel` failed during warmup invocation 0: injected kernel refusal"
        );

        let mut timed_failure = FailingRunKernel {
            calls: 0,
            fail_at: 1,
        };
        let error = measure(&mut timed_failure, 1, 2, &axes)
            .expect_err("timed refusal must fail the measurement");
        assert_eq!(
            error,
            "roofline kernel `receipt-kernel` failed during timed invocation 0: injected kernel refusal"
        );
    }

    #[test]
    fn later_registry_failure_aborts_every_kernel_tuning_state() {
        let axes = synthetic_axes();
        let first_pending = std::rc::Rc::new(std::cell::Cell::new(false));
        let second_pending = std::rc::Rc::new(std::cell::Cell::new(false));
        let first_aborts = std::rc::Rc::new(std::cell::Cell::new(0));
        let second_aborts = std::rc::Rc::new(std::cell::Cell::new(0));
        let mut registry: Vec<Box<dyn RooflineKernel>> = vec![
            Box::new(AbortProbeKernel {
                name: "abort-probe-first",
                fail: false,
                pending: std::rc::Rc::clone(&first_pending),
                aborts: std::rc::Rc::clone(&first_aborts),
            }),
            Box::new(AbortProbeKernel {
                name: "abort-probe-second",
                fail: true,
                pending: std::rc::Rc::clone(&second_pending),
                aborts: std::rc::Rc::clone(&second_aborts),
            }),
        ];

        let error = run_registry(&mut registry, 0, 1, &axes)
            .expect_err("a later kernel refusal must abort the entire registry");
        assert!(error.contains("abort-probe-second"), "{error}");
        assert!(!first_pending.get());
        assert!(!second_pending.get());
        assert_eq!(first_aborts.get(), 1);
        assert_eq!(second_aborts.get(), 1);
    }

    #[test]
    fn crushed_axes_cannot_gate_vacuously() {
        // Bead 1n61 counterexample replay: on a load-68 host both the
        // STREAM probe and the FFT kernel collapsed ~1000× together
        // (axes 0.2 GB/s, kernel 0.156 GB/s on a ~200 GB/s machine) and
        // the RATIO self-normalized to 0.89 = a vacuous within_band.
        let crushed = MachineAxes {
            fingerprint: 0xDEAD,
            cpu_brand: "crushed".to_string(),
            logical_cpus: 128,
            bandwidth_single_gbs: 0.2,
            bandwidth_all_core_gbs: 0.4,
            peak_single_gflops: 0.05,
            peak_all_core_gflops: 0.4,
        };
        assert!(!crushed.plausible());
        let spec = KernelSpec {
            name: "fft-roundtrip",
            version: "1n61",
            bytes_per_elem: 672.0,
            flops_per_elem: 172.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.40),
        };
        // 0.156 GB/s effective — 89% of the crushed axis.
        let a = attainment_for(&spec, 0.156e9 / 672.0, &crushed);
        assert_eq!(
            a.verdict,
            Verdict::EnvironmentInvalid,
            "a crushed environment must refuse to gate (got attainment {:.3})",
            a.attainment
        );
    }

    #[test]
    fn over_roof_attainment_poisons_the_gate() {
        // Healthy axes but the kernel 'beats' its roofline by 2× — the
        // axes are stale relative to the run; refuse.
        let spec = KernelSpec {
            name: "axpy",
            version: "1n61",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.6),
        };
        let a = attainment_for(&spec, 2.0 * 100.0e9 / 24.0, &synthetic_axes());
        assert_eq!(a.verdict, Verdict::EnvironmentInvalid);
        // Slightly over 1 (measurement jitter) still gates normally.
        let b = attainment_for(&spec, 1.2 * 100.0e9 / 24.0, &synthetic_axes());
        assert_eq!(b.verdict, Verdict::WithinBand);
    }

    #[test]
    fn invalid_numeric_inputs_fail_closed_and_remain_json() {
        let base = KernelSpec {
            name: "probe\"escaped",
            version: "1",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.5),
        };
        let bad_rate = attainment_with_dispersion(&base, f64::NAN, 0.0, 3, &synthetic_axes());
        assert_eq!(bad_rate.verdict, Verdict::EnvironmentInvalid);
        assert!(bad_rate.elems_per_sec.is_finite());
        assert!(bad_rate.to_jsonl().contains("probe\\\"escaped"));
        assert!(!bad_rate.to_jsonl().contains("NaN"));

        let bad_dispersion =
            attainment_with_dispersion(&base, 1.0, f64::INFINITY, 3, &synthetic_axes());
        assert_eq!(bad_dispersion.verdict, Verdict::EnvironmentInvalid);
        let bad_target = KernelSpec {
            target_fraction: Some(f64::NAN),
            ..base
        };
        assert_eq!(
            attainment_for(&bad_target, 1.0, &synthetic_axes()).verdict,
            Verdict::EnvironmentInvalid
        );
    }

    #[test]
    fn one_invalid_row_poisons_every_registry_verdict() {
        let normal = KernelSpec {
            name: "normal",
            version: "1",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.5),
        };
        let impossible = KernelSpec {
            name: "impossible",
            ..normal
        };
        let mut results = vec![
            attainment_for(&normal, 100.0e9 / 24.0 * 0.8, &synthetic_axes()),
            attainment_for(&impossible, 100.0e9 / 24.0 * 2.0, &synthetic_axes()),
        ];
        assert_eq!(results[0].verdict, Verdict::WithinBand);
        assert_eq!(results[1].verdict, Verdict::EnvironmentInvalid);
        poison_invalid_run(&mut results);
        assert!(
            results
                .iter()
                .all(|row| row.verdict == Verdict::EnvironmentInvalid)
        );
        assert!(results.iter().all(|row| {
            row.invalid_reason
                .as_deref()
                .is_some_and(|r| r.contains("impossible"))
        }));
    }

    #[test]
    fn only_bound_timed_receipts_are_citable() {
        let (axes, baseline, identity, timed) = timed_receipt_fixture();
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let spec = ReceiptKernel {
            elements: 1,
            value: 0,
        }
        .spec();
        let analytic = attainment_with_dispersion(&spec, 1.0, 0.0, 1, &axes);
        assert!(
            !run_passes_measurement_admission(&axes, &axes, baseline_policy, &[analytic]),
            "an analytic helper result is not measurement evidence"
        );

        assert!(run_passes_measurement_admission(
            &axes,
            &axes,
            baseline_policy,
            std::slice::from_ref(&timed)
        ));
        assert!(timed.to_jsonl().contains("\"sample_seconds_bits\""));
        let mut malformed_identity = timed.clone();
        malformed_identity.kernel = "probe\"},injected".to_string();
        let malformed_manifest =
            run_result_manifest_json(std::slice::from_ref(&malformed_identity));
        assert!(
            malformed_manifest.contains("probe\\\"},injected"),
            "manifest serialization must remain valid JSON before admission refusal"
        );
        assert!(
            run_admission_error(
                &axes,
                &axes,
                baseline_policy,
                std::slice::from_ref(&malformed_identity),
            )
            .is_some_and(|reason| reason.contains("non-canonical kernel/version"))
        );
        let mut drifted_post = axes.clone();
        drifted_post.bandwidth_single_gbs = 60.0;
        assert!(!run_passes_measurement_admission(
            &axes,
            &drifted_post,
            baseline_policy,
            std::slice::from_ref(&timed)
        ));

        let mut tampered = timed.clone();
        tampered.dispersion += 0.01;
        assert!(!run_passes_measurement_admission(
            &axes,
            &axes,
            baseline_policy,
            &[tampered]
        ));
        let mut tampered_target = timed.clone();
        tampered_target.target_attainment += 0.01;
        assert!(!run_passes_measurement_admission(
            &axes,
            &axes,
            baseline_policy,
            &[tampered_target]
        ));
        let mut tampered_axis = timed.clone();
        tampered_axis.spec_binding.target_axis = TargetAxis::ComputePeak;
        assert!(!run_passes_measurement_admission(
            &axes,
            &axes,
            baseline_policy,
            &[tampered_axis]
        ));
        assert!(
            run_admission_error(
                &axes,
                &axes,
                baseline_policy,
                &[timed.clone(), timed.clone()],
            )
            .is_some()
        );

        let mut empty = ReceiptKernel {
            elements: 0,
            value: 0,
        };
        let empty_row = measure(&mut empty, 0, 1, &axes).expect("bounded empty-kernel measurement");
        assert!(!run_passes_measurement_admission(
            &axes,
            &axes,
            baseline_policy,
            &[empty_row]
        ));

        let unbaselined = AxisBaselinePolicy::new(None, &identity, 20_010);
        assert!(
            !run_passes_measurement_admission(
                &axes,
                &axes,
                unbaselined,
                std::slice::from_ref(&timed)
            ),
            "first-run candidate evidence must never authorize itself"
        );
    }

    #[test]
    fn stable_sustained_contention_cannot_cross_the_production_admission_boundary() {
        let quiet = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&quiet);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut crushed = quiet.clone();
        crushed.bandwidth_single_gbs = 10.0;
        crushed.bandwidth_all_core_gbs = 40.0;
        crushed.peak_single_gflops = 10.0;
        crushed.peak_all_core_gflops = 60.0;
        assert!(crushed.plausible());
        assert!(crushed.reprobe_error(&crushed).is_none());

        let mut kernel = ReceiptKernel {
            elements: 1,
            value: 0,
        };
        let timed = measure(&mut kernel, 1, 3, &crushed).expect("bounded contention measurement");
        assert!(!run_passes_measurement_admission(
            &crushed,
            &crushed,
            baseline_policy,
            std::slice::from_ref(&timed)
        ));
        let receipt = baseline_policy.receipt_json(&crushed, &crushed);
        assert!(receipt.contains("\"baseline\":\"degraded\""));
        assert!(receipt.contains(&baseline.content_hash().to_hex()));
        assert!(!receipt.contains("NaN"));
    }

    #[test]
    fn registry_tuning_hook_uses_the_complete_admission_decision() {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut measured = ReceiptKernel {
            elements: 1,
            value: 0,
        };
        let result =
            measure(&mut measured, 0, 3, &axes).expect("bounded admission-hook measurement");
        let observed = std::rc::Rc::new(std::cell::Cell::new(None));
        let mut registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(AdmissionProbeKernel {
            observed: std::rc::Rc::clone(&observed),
        })];

        assert!(
            finalize_registry_tuning(
                &mut registry,
                &axes,
                &axes,
                baseline_policy,
                std::slice::from_ref(&result),
            )
            .expect("admitted hook")
            .admitted()
        );
        assert_eq!(observed.get(), Some(true));

        observed.set(None);
        let mut drifted_post = axes.clone();
        drifted_post.bandwidth_single_gbs *= 0.5;
        assert!(
            !finalize_registry_tuning(
                &mut registry,
                &axes,
                &drifted_post,
                baseline_policy,
                &[result],
            )
            .expect("rejected hook")
            .admitted()
        );
        assert_eq!(observed.get(), Some(false));
    }

    #[test]
    fn registry_tuning_hook_drains_every_kernel_after_a_middle_failure() {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let observed = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));
        let aborted = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));
        let mut registry: Vec<Box<dyn RooflineKernel>> = (0..3)
            .map(|id| {
                Box::new(FallibleAdmissionProbeKernel {
                    id,
                    observed: std::rc::Rc::clone(&observed),
                    aborted: std::rc::Rc::clone(&aborted),
                    failure: (id == 1).then_some("middle cleanup failed"),
                }) as Box<dyn RooflineKernel>
            })
            .collect();
        let results =
            run_registry(&mut registry, 0, 3, &axes).expect("bounded fallible-hook registry run");

        let error =
            finalize_registry_tuning(&mut registry, &axes, &axes, baseline_policy, &results)
                .expect_err("middle failure must be reported after every hook drains");
        assert_eq!(
            observed.borrow().as_slice(),
            &[(0, true), (1, true), (2, true)],
            "first, failing middle, and last kernel must see the same admission decision"
        );
        assert_eq!(
            aborted.borrow().as_slice(),
            &[0, 1, 2],
            "a failed lifecycle must abort process-local state in every kernel"
        );
        assert_eq!(
            error,
            "tuning lifecycle finalization failed with 1 issue(s): kernel[1]: middle cleanup failed"
        );
    }

    #[test]
    fn registry_finalization_refuses_missing_or_mismatched_kernels() {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut measured = ReceiptKernel {
            elements: 1,
            value: 0,
        };
        let result = measure(&mut measured, 0, 1, &axes).expect("bounded finalization measurement");

        let missing = finalize_registry_tuning(
            &mut [],
            &axes,
            &axes,
            baseline_policy,
            std::slice::from_ref(&result),
        )
        .expect_err("an empty registry cannot finalize a nonempty result set");
        assert!(missing.contains("length mismatch"));

        let observed = std::rc::Rc::new(std::cell::Cell::new(None));
        let mut mismatched: Vec<Box<dyn RooflineKernel>> = vec![Box::new(AdmissionProbeKernel {
            observed: std::rc::Rc::clone(&observed),
        })];
        let mut mismatched_result = result;
        mismatched_result.version = "different-version".to_string();
        let mismatch = finalize_registry_tuning(
            &mut mismatched,
            &axes,
            &axes,
            baseline_policy,
            &[mismatched_result],
        )
        .expect_err("a different kernel cannot finalize the measured row");
        assert!(mismatch.contains("identity mismatch"));
        assert_eq!(
            observed.get(),
            Some(false),
            "identity mismatch must be part of the aggregate admission decision"
        );
    }

    #[test]
    fn attainment_matches_hand_calculation_bandwidth_bound() {
        // axpy: 24 B/elem, 2 flop/elem on axes (100 GB/s, 50 GFLOP/s):
        // bw limit = 100e9/24 = 4.1667e9 elem/s; compute = 50e9/2 = 25e9.
        // Bandwidth binds. At 2.0833e9 elem/s attainment = 0.5 exactly.
        let spec = KernelSpec {
            name: "axpy",
            version: "1",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.6),
        };
        let a = attainment_for(&spec, 100.0e9 / 24.0 / 2.0, &synthetic_axes());
        assert_eq!(a.roof, RoofSide::Bandwidth);
        assert!((a.attainment - 0.5).abs() < 1e-12, "got {}", a.attainment);
        assert!((a.achieved_gbs - 50.0).abs() < 1e-9);
        assert_eq!(a.verdict, Verdict::BelowBand);
    }

    #[test]
    fn attainment_matches_hand_calculation_compute_bound() {
        // High-intensity kernel: 1 B/elem, 100 flop/elem.
        // bw limit = 100e9 elem/s; compute = 50e9/100 = 0.5e9 → compute binds.
        let spec = KernelSpec {
            name: "dense",
            version: "1",
            bytes_per_elem: 1.0,
            flops_per_elem: 100.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.5),
        };
        let a = attainment_for(&spec, 0.4e9, &synthetic_axes());
        assert_eq!(a.roof, RoofSide::Compute);
        assert!((a.attainment - 0.8).abs() < 1e-12);
        assert_eq!(a.verdict, Verdict::WithinBand);
        // All-core axes flip the limit.
        let all = KernelSpec {
            threading: Threading::AllCore,
            ..spec
        };
        let b = attainment_for(&all, 0.4e9, &synthetic_axes());
        assert!((b.limit_elems_per_sec - 3.0e9).abs() < 1.0);
    }

    #[test]
    fn compute_target_is_not_laundered_through_the_binding_bandwidth_roof() {
        let bandwidth_starved = MachineAxes {
            fingerprint: 0xC0DE,
            cpu_brand: "high-compute-low-bandwidth".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 10.0,
            bandwidth_all_core_gbs: 20.0,
            peak_single_gflops: 1_000.0,
            peak_all_core_gflops: 2_000.0,
        };
        let spec = KernelSpec {
            name: "gemm-target-model",
            version: "ss0n",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::ComputePeak,
            target_fraction: Some(0.75),
        };
        let above_bandwidth_target = attainment_for(&spec, 1.0e9, &bandwidth_starved);
        assert_eq!(above_bandwidth_target.roof, RoofSide::Bandwidth);
        assert!((above_bandwidth_target.attainment - 0.8).abs() < 1e-12);
        assert!((above_bandwidth_target.target_attainment - 0.001).abs() < 1e-12);
        assert_eq!(above_bandwidth_target.verdict, Verdict::BelowBand);
        assert!(
            above_bandwidth_target
                .to_jsonl()
                .contains("\"target_axis\":\"compute_peak\"")
        );

        let compute_bound = KernelSpec {
            bytes_per_elem: 1.0,
            flops_per_elem: 100.0,
            ..spec
        };
        let exact_boundary = attainment_for(&compute_bound, 0.375e9, &synthetic_axes());
        assert_eq!(exact_boundary.roof, RoofSide::Compute);
        assert_eq!(
            exact_boundary.target_attainment.to_bits(),
            0.75f64.to_bits()
        );
        assert_eq!(exact_boundary.verdict, Verdict::WithinBand);
    }

    #[test]
    fn no_target_reports_without_verdict() {
        let spec = KernelSpec {
            name: "probe",
            version: "1",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: None,
        };
        let a = attainment_for(&spec, 1.0e9, &synthetic_axes());
        assert_eq!(a.verdict, Verdict::NoTarget);
        assert!(a.to_jsonl().contains("\"verdict\":\"no_target\""));
    }

    #[test]
    fn section_14_1_table_is_complete_and_honest() {
        assert_eq!(SECTION_14_1_TARGETS.len(), 7, "all §14.1 families present");
        let registry = kernels::production_registry(1, &synthetic_axes())
            .expect("bounded production registry fixture");
        let registered: std::collections::BTreeSet<_> =
            registry.iter().map(|kernel| kernel.spec().name).collect();
        for row in SECTION_14_1_TARGETS {
            assert_eq!(
                row.landed,
                registered.contains(row.kernel),
                "{} target registration and landed status drifted",
                row.kernel,
            );
        }
    }

    #[test]
    fn versions_json_binds_the_exact_executable_identity() {
        let identity = fs_blake3::hash_domain("fs-roofline.test-executable.v1", b"binary");
        assert_eq!(
            versions_json(identity),
            format!("{{\"frankensim_executable\":\"{identity}\",\"fs-roofline\":\"{VERSION}\"}}")
        );
    }

    #[test]
    fn production_receipt_v3_feeds_exact_scoped_plan_cost_model() {
        let receipt_len = usize::try_from(MAX_DEPGRAPH_RECEIPT_BYTES).unwrap();
        let synthetic_receipt: &'static str =
            Box::leak(format!("{{\"d\":\"{}\"}}", "x".repeat(receipt_len - 8)).into_boxed_str());
        assert_eq!(
            synthetic_receipt.len(),
            receipt_len,
            "fixture must exercise the producer's exact 1 MiB receipt cap"
        );
        let axes = synthetic_axes();
        let (baseline, identity) = live_trusted_baseline(&axes);
        let baseline_policy = attested_policy(&baseline, &identity, test_day());
        let registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(ReceiptKernel {
            elements: 64,
            value: 0,
        })];
        let run = production::ProductionProbe::from_observed(axes.clone())
            .run_with_test_receipt(
                production::ProductionRunConfig {
                    n: 64,
                    warmup: 1,
                    reps: 3,
                },
                baseline_policy,
                registry,
                || axes.clone(),
                synthetic_receipt,
            )
            .expect("seal real receipt-v3 fixture");
        assert!(run.citation_eligible());
        let kernel = run.results()[0].kernel.clone();
        let version = run.results()[0].version.clone();
        let build_identity = read_executable_build_identity()
            .expect("test executable identity remains readable after measurement");
        let dependency = DependencyReceiptBinding::from_parts(
            synthetic_receipt,
            fs_blake3::hash_domain(
                fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN,
                synthetic_receipt.as_bytes(),
            ),
        )
        .expect("exact-cap dependency receipt agrees with its digest");
        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let op = run.record(&ledger).expect("record production evidence");
        let recorded_at = ledger
            .op(op)
            .expect("query exact-cap op")
            .expect("stored exact-cap op")
            .t_end
            .expect("exact-cap op finished");
        let rows = ledger.tune_rows(&kernel).expect("read exact tune row");
        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        let model = fs_plan::cost_model_from_tune(&ledger, &kernel, &row.shape_class, &row.machine)
            .expect("strict planner loader accepts producer-authored row");
        assert_eq!(model.n_obs(), 3, "every timed repetition is retained");
        assert!(model.predict(64.0).is_ok());
        assert_eq!(
            staleness_at_with_build_and_dependency(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                Some(baseline.content_hash()),
                recorded_at,
                build_identity,
                Some(dependency),
            )
            .expect("classify exact-cap retained receipt"),
            Staleness::Fresh
        );
    }

    #[test]
    fn retained_over_cap_dependency_is_refused_by_both_consumers() {
        let receipt_len = usize::try_from(MAX_DEPGRAPH_RECEIPT_BYTES).unwrap() + 1;
        let synthetic_receipt: &'static str =
            Box::leak(format!("{{\"d\":\"{}\"}}", "x".repeat(receipt_len - 8)).into_boxed_str());
        let axes = synthetic_axes();
        let (baseline, identity) = live_trusted_baseline(&axes);
        let baseline_policy = attested_policy(&baseline, &identity, test_day());
        let registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(ReceiptKernel {
            elements: 64,
            value: 0,
        })];
        let run = production::ProductionProbe::from_observed(axes.clone())
            .run_with_test_receipt(
                production::ProductionRunConfig {
                    n: 64,
                    warmup: 1,
                    reps: 3,
                },
                baseline_policy,
                registry,
                || axes.clone(),
                synthetic_receipt,
            )
            .expect("seal hostile over-cap retained receipt");
        let kernel = run.results()[0].kernel.clone();
        let version = run.results()[0].version.clone();
        let build_identity = read_executable_build_identity()
            .expect("test executable identity remains readable after measurement");
        let dependency = DependencyReceiptBinding::from_parts(
            synthetic_receipt,
            fs_blake3::hash_domain(
                fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN,
                synthetic_receipt.as_bytes(),
            ),
        )
        .expect("over-cap fixture still agrees with its retained digest");
        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let op = run
            .record(&ledger)
            .expect("retain hostile dependency receipt");
        let recorded_at = ledger
            .op(op)
            .expect("query hostile receipt op")
            .expect("stored hostile receipt op")
            .t_end
            .expect("hostile receipt op finished");
        let rows = ledger
            .tune_rows(&kernel)
            .expect("read hostile exact tune row");
        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert!(matches!(
            fs_plan::cost_model_from_tune(&ledger, &kernel, &row.shape_class, &row.machine),
            Err(fs_plan::TuneModelError::Ledger(
                LedgerError::ArtifactReadLimit {
                    limit,
                    observed,
                    ..
                }
            )) if limit == MAX_DEPGRAPH_RECEIPT_BYTES
                && observed == MAX_DEPGRAPH_RECEIPT_BYTES + 1
        ));
        assert_eq!(
            staleness_at_with_build_and_dependency(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                Some(baseline.content_hash()),
                recorded_at,
                build_identity,
                Some(dependency),
            )
            .expect("classify hostile retained receipt"),
            Staleness::CorruptEvidence
        );
    }

    #[test]
    fn dependency_receipt_cap_matches_the_fs_la_producer_source() {
        let source = include_str!("../../fs-la/depgraph_receipt_format.rs");
        let declaration = source
            .lines()
            .find_map(|line| {
                line.trim()
                    .strip_prefix("pub const MAX_RECEIPT_BYTES: usize = ")
            })
            .and_then(|value| value.strip_suffix(';'))
            .expect("fs-la producer must declare MAX_RECEIPT_BYTES");
        let producer_cap = declaration
            .replace('_', "")
            .parse::<u64>()
            .expect("fs-la producer cap must remain a decimal byte count");
        assert_eq!(producer_cap, MAX_DEPGRAPH_RECEIPT_BYTES);
    }

    #[test]
    fn staleness_refuses_a_different_current_executable() {
        const SYNTHETIC_RECEIPT: &str = "{\"schema\":\"fs-roofline-synthetic-dependency-receipt-v1\",\"purpose\":\"build-drift-unit-test\"}";
        let axes = synthetic_axes();
        let (baseline, identity) = live_trusted_baseline(&axes);
        let baseline_policy = attested_policy(&baseline, &identity, test_day());
        let registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(ReceiptKernel {
            elements: 1,
            value: 0,
        })];
        let run = production::ProductionProbe::from_observed(axes.clone())
            .run_with_test_receipt(
                production::ProductionRunConfig {
                    n: 1,
                    warmup: 0,
                    reps: 1,
                },
                baseline_policy,
                registry,
                || axes.clone(),
                SYNTHETIC_RECEIPT,
            )
            .expect("seal fixture");
        assert!(
            run.citation_eligible(),
            "fixture must establish production provenance"
        );
        let kernel = run.results()[0].kernel.clone();
        let version = run.results()[0].version.clone();
        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let op = run.record(&ledger).expect("record fixture");
        let recorded_at = ledger
            .op(op)
            .expect("query op")
            .expect("stored op")
            .t_end
            .expect("finished op");
        let foreign_build = fs_blake3::hash_domain(
            "org.frankensim.fs-roofline.foreign-executable.test.v1",
            b"rebuilt binary",
        );
        let dependency = DependencyReceiptBinding::from_parts(
            SYNTHETIC_RECEIPT,
            fs_blake3::hash_domain(
                fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN,
                SYNTHETIC_RECEIPT.as_bytes(),
            ),
        )
        .expect("synthetic receipt digest agrees");
        assert_eq!(
            staleness_at_with_build_and_dependency(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                Some(baseline.content_hash()),
                recorded_at,
                foreign_build,
                Some(dependency),
            )
            .expect("classify"),
            Staleness::BuildDrift
        );
    }

    #[test]
    fn bit_identical_reruns_keep_distinct_operation_bound_rows() {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut first_registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(ReceiptKernel {
            elements: 1,
            value: 0,
        })];
        let mut first_results = run_registry(&mut first_registry, 0, 1, &axes)
            .expect("bounded first identical registry run");
        let mut second_results = first_results.clone();
        let mut first_finalized = finalize_registry_tuning(
            &mut first_registry,
            &axes,
            &axes,
            baseline_policy,
            &first_results,
        )
        .expect("finalize first identical run");
        let mut second_registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(ReceiptKernel {
            elements: 1,
            value: 0,
        })];
        let mut second_finalized = finalize_registry_tuning(
            &mut second_registry,
            &axes,
            &axes,
            baseline_policy,
            &second_results,
        )
        .expect("finalize second identical run");
        assert_eq!(
            first_finalized.receipt_identity(),
            second_finalized.receipt_identity(),
            "fixture must exercise an exact receipt collision"
        );

        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let first_op = record_run(
            &ledger,
            &axes,
            &axes,
            baseline_policy,
            &mut first_finalized,
            &mut first_results,
        )
        .expect("record first identical run");
        let second_op = record_run(
            &ledger,
            &axes,
            &axes,
            baseline_policy,
            &mut second_finalized,
            &mut second_results,
        )
        .expect("record second identical run");
        assert_ne!(first_op, second_op);

        let rows = ledger
            .tune_rows("receipt-kernel")
            .expect("query retained run history");
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().any(|row| {
            row.shape_class
                == candidate_measurement_shape_class(
                    "1",
                    first_finalized.receipt_identity(),
                    first_op,
                )
        }));
        assert!(rows.iter().any(|row| {
            row.shape_class
                == candidate_measurement_shape_class(
                    "1",
                    second_finalized.receipt_identity(),
                    second_op,
                )
        }));
    }

    #[test]
    fn historical_dependency_receipt_does_not_poison_a_current_build_row() {
        let current_build = fs_blake3::hash_domain("fs-roofline.test-build.v1", b"current");
        let historical_build = fs_blake3::hash_domain("fs-roofline.test-build.v1", b"historical");
        let mut scan = BuildRowScan::default();
        let zero = fs_blake3::hash_domain("fs-roofline.test-dep.v1", b"placeholder");
        assert!(scan.observe(
            &ValidatedRooflineRow {
                build_identity: historical_build,
                recorded_at_ns: 10,
                dependency_matches_current: false,
                dependency_receipt_digest: zero,
                dependency_receipt_artifact: zero,
            },
            current_build,
        ));
        assert!(scan.observe(
            &ValidatedRooflineRow {
                build_identity: current_build,
                recorded_at_ns: 20,
                dependency_matches_current: true,
                dependency_receipt_digest: zero,
                dependency_receipt_artifact: zero,
            },
            current_build,
        ));
        assert!(scan.saw_foreign_build);
        assert_eq!(scan.newest_current_build, Some(20));

        assert!(!scan.observe(
            &ValidatedRooflineRow {
                build_identity: current_build,
                recorded_at_ns: 30,
                dependency_matches_current: false,
                dependency_receipt_digest: zero,
                dependency_receipt_artifact: zero,
            },
            current_build,
        ));
    }

    #[test]
    fn production_operation_envelope_binds_all_five_explicit_columns_and_diagnostic() {
        fn assert_mutation_refused(
            valid: &fs_ledger::OpRow,
            build: fs_blake3::ContentHash,
            mutate: impl FnOnce(&mut fs_ledger::OpRow),
        ) {
            let mut altered = valid.clone();
            mutate(&mut altered);
            assert!(!production_op_envelope_is_valid(&altered, build));
        }

        let build = fs_blake3::hash_domain("fs-roofline.test-build.v1", b"envelope");
        let valid = fs_ledger::OpRow {
            id: 7,
            session: Some(b"roofline".to_vec()),
            ir: "{}".to_string(),
            seed: ROOFLINE_SEED.to_vec(),
            versions: versions_json(build),
            budget: ROOFLINE_BUDGET.to_string(),
            capability: ROOFLINE_CAPABILITY.to_string(),
            t_start: 1,
            t_end: Some(2),
            outcome: Some("ok".to_string()),
            diag: None,
        };
        assert!(production_op_envelope_is_valid(&valid, build));

        assert_mutation_refused(&valid, build, |op| op.session = Some(b"other".to_vec()));
        assert_mutation_refused(&valid, build, |op| op.seed = b"other".to_vec());
        assert_mutation_refused(&valid, build, |op| op.versions = "{}".to_string());
        assert_mutation_refused(&valid, build, |op| op.budget = "{}".to_string());
        assert_mutation_refused(&valid, build, |op| op.capability = "{}".to_string());
        assert_mutation_refused(&valid, build, |op| {
            op.outcome = Some("error".to_string());
        });
        assert_mutation_refused(&valid, build, |op| op.diag = Some("{}".to_string()));
        assert_mutation_refused(&valid, build, |op| op.t_end = Some(0));
    }

    #[test]
    fn executable_identity_drift_is_fail_closed() {
        let captured = fs_blake3::hash_domain("fs-roofline.test-build.v1", b"captured");
        let drifted = fs_blake3::hash_domain("fs-roofline.test-build.v1", b"drifted");
        assert_eq!(
            require_stable_executable_identity(captured, captured).expect("stable identity"),
            captured
        );
        let error = require_stable_executable_identity(captured, drifted)
            .expect_err("identity drift must refuse recording");
        assert!(error.to_string().contains("drifted between"));
    }
}
