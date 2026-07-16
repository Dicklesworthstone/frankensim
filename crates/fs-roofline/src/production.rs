//! Sealed production-run protocol (bead fz2.5): the only path to citable
//! roofline evidence.
//!
//! The public `RooflineKernel`/`run_registry`/`record_run` surface treats
//! caller-supplied kernel implementations and `MachineAxes` as a trust root:
//! a caller can clone a valid GEMM `KernelExecutionBinding` into a fake
//! custom kernel named `gemm-f64`, or discard a drifted post-probe and pass
//! the pre-probe twice. That surface stays available for harness tests and
//! exploration, but everything it records is stamped
//! `"protocol":"custom-registry"` and is explicitly NON-CITABLE.
//!
//! The protocol is three opaque stages:
//!
//! 1. [`ProductionProbe::observe`] performs the pre-run axis probe and mints
//!    the per-run nonce. The caller may READ the observed axes (baseline
//!    selection needs them) but can never supply its own.
//! 2. [`ProductionProbe::run`] owns production registry selection, timed
//!    warmup/repetitions, the post-run axis probe (observed strictly after
//!    the timed loop), aggregate admission, and tune finalization, yielding
//!    a [`ProductionRun`].
//! 3. [`ProductionRun::record`] commits atomically and consumes the run,
//!    yielding neutral [`RecordedProductionRun`]. Its operation `ir` carries
//!    `"protocol":"production-v3"`, the nonce, content hashes of both
//!    observed axis receipts, and the retained dependency-receipt binding.
//!    Only [`RecordedProductionRun::revalidate`] can add current authority,
//!    dependency, build, clock, and exact-ledger proof and mint
//!    [`FreshProductionEvidence`].
//!
//! Trust model: the nonce is a process-unique challenge, not cryptographic
//! proof. Type opacity prevents ordinary API consumers from constructing a
//! `ProductionRun`, but `fs-ledger` intentionally exposes general mutation
//! APIs. A trusted ledger writer can therefore mint or replace internally
//! consistent rows. External authentication of the ledger/package is a
//! separate proof obligation; this crate detects corruption inside that
//! trusted-writer boundary and makes no cryptographic-authority claim.

use std::collections::BTreeSet;

use fs_ledger::{Ledger, LedgerError};

use crate::kernels::production_registry_with_ledger;
use crate::{
    Attainment, AttestedAxisBaselinePolicy, AxisAdmissionSnapshot, AxisBaselinePolicy,
    CUSTOM_REGISTRY_PROTOCOL_FIELD, DependencyReceiptBinding, FinalizedRegistryRun, MachineAxes,
    PRODUCTION_PROTOCOL_FIELD, RooflineKernel, citable_run_admission_error_for_snapshot,
    finalize_registry_tuning, finalize_registry_tuning_with_snapshot, json_escape,
    record_run_with_protocol, run_admission_error_for_snapshot, run_registry,
};

const RUN_NONCE_DOMAIN: &str = "org.frankensim.fs-roofline.production-run-nonce.v1";
/// Semantic version of a production machine-axis observation receipt.
pub const PRODUCTION_AXES_RECEIPT_IDENTITY_VERSION: u32 = 1;
/// BLAKE3 derive-key domain for a production machine-axis observation receipt.
pub const PRODUCTION_AXES_RECEIPT_DOMAIN: &str =
    "org.frankensim.fs-roofline.production-axes-receipt.v1";
/// Semantic version of a live dependency-authority policy receipt.
pub const DEPENDENCY_AUTHORITY_POLICY_IDENTITY_VERSION: u32 = 1;
/// Domain for the exact dependency-authority policy sampled by revalidation.
pub const DEPENDENCY_AUTHORITY_POLICY_DOMAIN: &str =
    "frankensim.fs-roofline.dependency-authority-policy.v1";
/// Maximum accepted size of a configured dependency-authority policy.
pub const MAX_DEPENDENCY_AUTHORITY_POLICY_BYTES: usize = 1024 * 1024;

/// Owner-local dependency-authority declaration consumed by
/// `xtask check-identities`.
#[allow(dead_code)]
pub const DEPENDENCY_AUTHORITY_POLICY_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-roofline:dependency-authority-policy",
    "version_const=DEPENDENCY_AUTHORITY_POLICY_IDENTITY_VERSION",
    "version=1",
    "domain=frankensim.fs-roofline.dependency-authority-policy.v1",
    "domain_const=DEPENDENCY_AUTHORITY_POLICY_DOMAIN",
    "encoder=dependency_authority_policy_receipt",
    "encoder_helpers=dependency_authority_policy_receipt_with_domain",
    "schema_constants=DEPENDENCY_AUTHORITY_POLICY_IDENTITY_VERSION,DEPENDENCY_AUTHORITY_POLICY_DOMAIN,MAX_DEPENDENCY_AUTHORITY_POLICY_BYTES",
    "schema_functions=ConfiguredDependencyReceiptAuthority::from_text,ConfiguredDependencyReceiptAuthority::policy_receipt,ConfiguredDependencyReceiptAuthority::verify,DependencyReceiptDecision::new,DependencyReceiptDecision::policy_receipt,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=DependencyAuthorityPolicyIdentityInput",
    "source_fields=DependencyAuthorityPolicyIdentityInput.canonical_bytes:semantic",
    "source_bindings=DependencyAuthorityPolicyIdentityInput.canonical_bytes>canonical-policy-bytes",
    "external_semantic_fields=digest-domain,identity-version",
    "semantic_fields=digest-domain,identity-version,canonical-policy-bytes",
    "excluded_fields=none",
    "consumers=ConfiguredDependencyReceiptAuthority::from_text,ConfiguredDependencyReceiptAuthority::verify,RecordedProductionRun::revalidate,FreshProductionEvidence::dependency_authority_fingerprint,crates/fs-roofline/src/bin/roofline.rs#load_dependency_authority",
    "mutations=digest-domain:crates/fs-roofline/src/production.rs#dependency_authority_policy_identity_fields_move_independently,identity-version:crates/fs-roofline/src/production.rs#dependency_authority_policy_identity_versions_fail_closed,canonical-policy-bytes:crates/fs-roofline/src/production.rs#dependency_authority_policy_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_dependency_authority_policy_identity_fields",
    "transport_guard=ConfiguredDependencyReceiptAuthority::from_text",
    "version_guard=crates/fs-roofline/src/production.rs#dependency_authority_policy_identity_versions_fail_closed",
    "coupling_surface=fs-roofline:dependency-authority-policy",
];

/// Owner-local production-axis declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const PRODUCTION_AXES_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-roofline:production-axes-receipt",
    "version_const=PRODUCTION_AXES_RECEIPT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-roofline.production-axes-receipt.v1",
    "domain_const=PRODUCTION_AXES_RECEIPT_DOMAIN",
    "encoder=axes_receipt",
    "encoder_helpers=ProductionAxesReceiptInput::from_axes,production_axes_receipt_json,machine_axes_receipt_json,axes_receipt_with_domain",
    "schema_constants=PRODUCTION_AXES_RECEIPT_IDENTITY_VERSION,PRODUCTION_AXES_RECEIPT_DOMAIN",
    "schema_functions=production_axes_receipt_is_current,crates/fs-roofline/src/lib.rs#json_escape,crates/fs-roofline/src/lib.rs#parse_machine_axes_receipt,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=ProductionAxesReceiptInput",
    "source_fields=ProductionAxesReceiptInput.fingerprint:semantic,ProductionAxesReceiptInput.cpu_brand:semantic,ProductionAxesReceiptInput.logical_cpus:semantic,ProductionAxesReceiptInput.bandwidth_single_bits:semantic,ProductionAxesReceiptInput.bandwidth_all_core_bits:semantic,ProductionAxesReceiptInput.peak_single_bits:semantic,ProductionAxesReceiptInput.peak_all_core_bits:semantic",
    "source_bindings=ProductionAxesReceiptInput.fingerprint>machine-fingerprint,ProductionAxesReceiptInput.cpu_brand>cpu-brand-utf8,ProductionAxesReceiptInput.logical_cpus>logical-cpus,ProductionAxesReceiptInput.bandwidth_single_bits>bandwidth-single-bits,ProductionAxesReceiptInput.bandwidth_all_core_bits>bandwidth-all-core-bits,ProductionAxesReceiptInput.peak_single_bits>peak-single-bits,ProductionAxesReceiptInput.peak_all_core_bits>peak-all-core-bits",
    "external_semantic_fields=digest-domain,identity-version",
    "semantic_fields=digest-domain,identity-version,machine-fingerprint,cpu-brand-utf8,logical-cpus,bandwidth-single-bits,bandwidth-all-core-bits,peak-single-bits,peak-all-core-bits",
    "excluded_fields=none",
    "consumers=ProductionRun::record,ReportOnlyProductionRun::record,validate_protocol_axes,AxisAdmissionSnapshot::receipt_json",
    "mutations=digest-domain:crates/fs-roofline/src/production.rs#production_axes_receipt_identity_fields_move_independently,identity-version:crates/fs-roofline/src/production.rs#production_axes_receipt_versions_fail_closed,machine-fingerprint:crates/fs-roofline/src/production.rs#production_axes_receipt_identity_fields_move_independently,cpu-brand-utf8:crates/fs-roofline/src/production.rs#production_axes_receipt_identity_fields_move_independently,logical-cpus:crates/fs-roofline/src/production.rs#production_axes_receipt_identity_fields_move_independently,bandwidth-single-bits:crates/fs-roofline/src/production.rs#production_axes_receipt_identity_fields_move_independently,bandwidth-all-core-bits:crates/fs-roofline/src/production.rs#production_axes_receipt_identity_fields_move_independently,peak-single-bits:crates/fs-roofline/src/production.rs#production_axes_receipt_identity_fields_move_independently,peak-all-core-bits:crates/fs-roofline/src/production.rs#production_axes_receipt_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_production_axes_receipt_identity_fields",
    "transport_guard=production_axes_receipt_is_current",
    "version_guard=crates/fs-roofline/src/production.rs#production_axes_receipt_versions_fail_closed",
    "coupling_surface=fs-roofline:production-axes-receipt",
];

#[derive(Debug, Clone, PartialEq, Eq)]
struct ProductionAxesReceiptInput {
    fingerprint: u64,
    cpu_brand: String,
    logical_cpus: u32,
    bandwidth_single_bits: u64,
    bandwidth_all_core_bits: u64,
    peak_single_bits: u64,
    peak_all_core_bits: u64,
}

impl ProductionAxesReceiptInput {
    fn from_axes(axes: &MachineAxes) -> Self {
        Self {
            fingerprint: axes.fingerprint,
            cpu_brand: axes.cpu_brand.clone(),
            logical_cpus: axes.logical_cpus,
            bandwidth_single_bits: axes.bandwidth_single_gbs.to_bits(),
            bandwidth_all_core_bits: axes.bandwidth_all_core_gbs.to_bits(),
            peak_single_bits: axes.peak_single_gflops.to_bits(),
            peak_all_core_bits: axes.peak_all_core_gflops.to_bits(),
        }
    }
}

#[allow(dead_code)]
fn classify_production_axes_receipt_identity_fields(input: &ProductionAxesReceiptInput) {
    let ProductionAxesReceiptInput {
        fingerprint,
        cpu_brand,
        logical_cpus,
        bandwidth_single_bits,
        bandwidth_all_core_bits,
        peak_single_bits,
        peak_all_core_bits,
    } = input;
    let _ = (
        fingerprint,
        cpu_brand,
        logical_cpus,
        bandwidth_single_bits,
        bandwidth_all_core_bits,
        peak_single_bits,
        peak_all_core_bits,
    );
}

struct DependencyAuthorityPolicyIdentityInput<'a> {
    canonical_bytes: &'a [u8],
}

fn dependency_authority_policy_receipt(canonical_bytes: &[u8]) -> fs_blake3::ContentHash {
    dependency_authority_policy_receipt_with_domain(
        DEPENDENCY_AUTHORITY_POLICY_DOMAIN,
        &DependencyAuthorityPolicyIdentityInput { canonical_bytes },
    )
}

fn dependency_authority_policy_receipt_with_domain(
    domain: &str,
    input: &DependencyAuthorityPolicyIdentityInput<'_>,
) -> fs_blake3::ContentHash {
    fs_blake3::hash_domain(domain, input.canonical_bytes)
}

#[allow(dead_code)]
fn classify_dependency_authority_policy_identity_fields(
    input: &DependencyAuthorityPolicyIdentityInput<'_>,
) {
    let DependencyAuthorityPolicyIdentityInput { canonical_bytes: _ } = input;
}
#[cfg(test)]
const DEVELOPMENT_SALT_REFUSAL: &str = "dependency graph uses the development equivalence salt; production citation requires an exact operator-observed normal/build receipt";

static NONCE_COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
const REPORT_ONLY_REFUSAL: &str =
    "operator-trusted candidate baseline has no attested promotion authority";
const MAX_REPORT_ONLY_REFUSAL_BYTES: usize = 4096;
const REPORT_ONLY_REFUSAL_DIGEST_DOMAIN: &str = "org.frankensim.fs-roofline.report-only-refusal.v1";
const SEALED_FINALIZATION_REFUSAL: &str = "sealed registry finalization did not admit this run";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CitationAuthority {
    Receipt(DependencyReceiptBinding),
    Refused(&'static str),
}

impl CitationAuthority {
    fn from_build() -> Self {
        DependencyReceiptBinding::current().map_or_else(Self::Refused, Self::Receipt)
    }

    const fn refusal(self) -> Option<&'static str> {
        match self {
            Self::Receipt(_) => None,
            Self::Refused(reason) => Some(reason),
        }
    }

    const fn receipt(self) -> Option<DependencyReceiptBinding> {
        match self {
            Self::Receipt(binding) => Some(binding),
            Self::Refused(_) => None,
        }
    }
}

/// Explicit live trust roots used to revalidate one recorded production run:
/// the baseline/promotion policy, source-retention inventory, and independent
/// dependency-receipt revocation authority.
/// Construction does no verification; [`RecordedProductionRun::revalidate`]
/// samples the authority exactly once after authenticating the exact ledgered
/// operation and manifest.
pub struct ProductionFreshnessContext<'a> {
    baselines: &'a crate::AttestedBaselineStore,
    authority: &'a dyn crate::PromotionAuthorityVerifier,
    retained_sources: &'a BTreeSet<fs_blake3::ContentHash>,
    dependency_authority: &'a dyn DependencyReceiptAuthority,
}

impl<'a> ProductionFreshnessContext<'a> {
    /// Pin the current baseline store, promotion authority, and retained
    /// source inventory for one named revalidation boundary.
    #[must_use]
    pub const fn new(
        baselines: &'a crate::AttestedBaselineStore,
        authority: &'a dyn crate::PromotionAuthorityVerifier,
        retained_sources: &'a BTreeSet<fs_blake3::ContentHash>,
        dependency_authority: &'a dyn DependencyReceiptAuthority,
    ) -> Self {
        Self {
            baselines,
            authority,
            retained_sources,
            dependency_authority,
        }
    }
}

/// Live operator verdict for the exact dependency receipt compiled into the
/// executable asking for freshness.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DependencyReceiptVerdict {
    /// The receipt remains authorized for production citation.
    Authorized,
    /// The operator revoked the receipt even though its retained bytes remain
    /// structurally valid.
    Revoked,
}

/// One atomic answer from a dependency-receipt authority. The verdict and
/// exact policy receipt travel together so a caller cannot mix a decision
/// sampled under one revocation policy with another policy's identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DependencyReceiptDecision {
    verdict: DependencyReceiptVerdict,
    policy_receipt: fs_blake3::ContentHash,
}

impl DependencyReceiptDecision {
    /// Bind `verdict` to the exact live authority policy that produced it.
    #[must_use]
    pub const fn new(
        verdict: DependencyReceiptVerdict,
        policy_receipt: fs_blake3::ContentHash,
    ) -> Self {
        Self {
            verdict,
            policy_receipt,
        }
    }

    /// Typed live authority verdict.
    #[must_use]
    pub const fn verdict(self) -> DependencyReceiptVerdict {
        self.verdict
    }

    /// Content identity of the exact revocation policy used for the verdict.
    #[must_use]
    pub const fn policy_receipt(self) -> fs_blake3::ContentHash {
        self.policy_receipt
    }
}

/// Injected live authority for dependency-receipt revocation.
pub trait DependencyReceiptAuthority {
    /// Judge the exact domain digest and artifact identity retained by the
    /// recorded operation. Revalidation calls this once per named boundary.
    #[must_use]
    fn verify(
        &self,
        digest: fs_blake3::ContentHash,
        artifact: fs_blake3::ContentHash,
    ) -> DependencyReceiptDecision;
}

/// Bounded operator policy for dependency-receipt revocation.
///
/// The canonical format is a strictly ascending list of 64-character
/// lowercase dependency-receipt digests, one per newline-terminated line. An
/// empty file is the explicit no-revocations policy. The exact canonical bytes
/// are hashed into every [`DependencyReceiptDecision`]; listed digests return
/// [`DependencyReceiptVerdict::Revoked`] and every other structurally verified
/// digest returns [`DependencyReceiptVerdict::Authorized`].
#[derive(Debug)]
pub struct ConfiguredDependencyReceiptAuthority {
    revoked_digests: BTreeSet<fs_blake3::ContentHash>,
    policy_receipt: fs_blake3::ContentHash,
}

impl ConfiguredDependencyReceiptAuthority {
    /// Parse one immutable canonical revocation policy.
    ///
    /// # Errors
    /// Refuses oversized, non-canonical, duplicated, unsorted, or malformed
    /// input.
    pub fn from_text(text: &str) -> Result<Self, String> {
        if text.len() > MAX_DEPENDENCY_AUTHORITY_POLICY_BYTES {
            return Err(format!(
                "dependency-authority policy exceeds the {MAX_DEPENDENCY_AUTHORITY_POLICY_BYTES}-byte bound"
            ));
        }
        let mut revoked_digests = BTreeSet::new();
        if !text.is_empty() {
            let body = text.strip_suffix('\n').ok_or_else(|| {
                "dependency-authority policy must be canonical newline-terminated lowercase hex"
                    .to_string()
            })?;
            if body.is_empty() {
                return Err(
                    "dependency-authority policy uses an empty file, not a blank line, for no revocations"
                        .to_string(),
                );
            }
            let mut previous = None;
            for (index, line) in body.split('\n').enumerate() {
                if line.len() != 64
                    || !line
                        .bytes()
                        .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
                {
                    return Err(format!(
                        "dependency-authority policy line {} must be exactly 64 lowercase hexadecimal bytes",
                        index + 1
                    ));
                }
                let digest = fs_blake3::ContentHash::from_hex(line).ok_or_else(|| {
                    format!(
                        "dependency-authority policy line {} is not a content hash",
                        index + 1
                    )
                })?;
                if previous.is_some_and(|prior| digest <= prior) {
                    return Err(format!(
                        "dependency-authority policy line {} is not in strict ascending order",
                        index + 1
                    ));
                }
                previous = Some(digest);
                let inserted = revoked_digests.insert(digest);
                debug_assert!(inserted);
            }
        }
        Ok(Self {
            revoked_digests,
            policy_receipt: dependency_authority_policy_receipt(text.as_bytes()),
        })
    }

    /// Content identity of the exact canonical policy bytes.
    #[must_use]
    pub const fn policy_receipt(&self) -> fs_blake3::ContentHash {
        self.policy_receipt
    }
}

impl DependencyReceiptAuthority for ConfiguredDependencyReceiptAuthority {
    fn verify(
        &self,
        digest: fs_blake3::ContentHash,
        _artifact: fs_blake3::ContentHash,
    ) -> DependencyReceiptDecision {
        let verdict = if self.revoked_digests.contains(&digest) {
            DependencyReceiptVerdict::Revoked
        } else {
            DependencyReceiptVerdict::Authorized
        };
        DependencyReceiptDecision::new(verdict, self.policy_receipt())
    }
}

/// Why an opaque recorded receipt could not mint current positive evidence.
#[derive(Debug)]
pub enum ProductionFreshnessError {
    /// The sealed run durably recorded an admission refusal, not evidence.
    RecordedRefusal,
    /// The exact operation, manifest, rows, artifacts, or typed receipt no
    /// longer agree.
    CorruptRecordedEvidence,
    /// The observation clock predates the operation completion receipt.
    ClockRollback {
        /// Completion time authenticated from the recorded operation.
        recorded_at_ns: i64,
        /// Live wall-clock observation used for revalidation.
        observed_ns: i64,
    },
    /// The exact evidence is older than the supported freshness window.
    Expired {
        /// Non-negative live age of the recorded operation.
        age_ns: i64,
    },
    /// The current attested store has no baseline for the recorded machine.
    BaselineUnavailable,
    /// The current store replaced the exact baseline record.
    BaselineReplaced,
    /// The current store replaced or removed the recorded attestation.
    PromotionAttestationChanged,
    /// A recorded source receipt is absent from the current retention set.
    SourceReceiptUnavailable(
        /// Missing source receipt identity.
        fs_blake3::ContentHash,
    ),
    /// The current promotion authority refused the recorded attestation.
    PromotionAuthorityRefused {
        /// Exact live authority verdict.
        verdict: crate::KeyVerdict,
    },
    /// The verifier's exact current policy differs from the recorded policy.
    PromotionPolicyChanged,
    /// No current operator-observed dependency receipt is available.
    DependencyAuthorityUnavailable,
    /// The current dependency receipt differs from the recorded one.
    DependencyAuthorityChanged,
    /// The live dependency authority explicitly revoked the exact receipt.
    DependencyAuthorityRevoked,
    /// The executable asking for freshness differs from the recorded build.
    BuildDrift,
    /// Durable ledger access failed.
    Ledger(
        /// Underlying durable-ledger failure.
        LedgerError,
    ),
}

impl core::fmt::Display for ProductionFreshnessError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::RecordedRefusal => write!(f, "recorded roofline run was refused at admission"),
            Self::CorruptRecordedEvidence => {
                write!(f, "recorded roofline operation or typed receipt is corrupt")
            }
            Self::ClockRollback {
                recorded_at_ns,
                observed_ns,
            } => write!(
                f,
                "roofline freshness clock rolled back: recorded at {recorded_at_ns}, observed {observed_ns}"
            ),
            Self::Expired { age_ns } => {
                write!(f, "recorded roofline evidence expired at age {age_ns} ns")
            }
            Self::BaselineUnavailable => {
                write!(f, "current attested baseline is unavailable")
            }
            Self::BaselineReplaced => write!(f, "current baseline replaced the recorded baseline"),
            Self::PromotionAttestationChanged => write!(
                f,
                "current baseline attestation differs from the recorded attestation"
            ),
            Self::SourceReceiptUnavailable(receipt) => {
                write!(
                    f,
                    "recorded baseline source receipt {receipt} is unavailable"
                )
            }
            Self::PromotionAuthorityRefused { verdict } => write!(
                f,
                "current promotion authority refused the recorded attestation: {}",
                verdict.name()
            ),
            Self::PromotionPolicyChanged => {
                write!(
                    f,
                    "current promotion-authority policy differs from the recorded policy"
                )
            }
            Self::DependencyAuthorityUnavailable => {
                write!(f, "current dependency receipt authority is unavailable")
            }
            Self::DependencyAuthorityChanged => {
                write!(
                    f,
                    "current dependency receipt differs from the recorded receipt"
                )
            }
            Self::DependencyAuthorityRevoked => {
                write!(f, "live dependency authority revoked the recorded receipt")
            }
            Self::BuildDrift => write!(f, "current executable differs from the recorded build"),
            Self::Ledger(error) => write!(f, "roofline freshness ledger failure: {error}"),
        }
    }
}

impl core::error::Error for ProductionFreshnessError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Ledger(error) => Some(error),
            _ => None,
        }
    }
}

impl From<LedgerError> for ProductionFreshnessError {
    fn from(error: LedgerError) -> Self {
        Self::Ledger(error)
    }
}

/// Neutral proof that one sealed production operation committed. Recording is
/// not freshness: only [`Self::revalidate`] can mint
/// [`FreshProductionEvidence`]. Fields are private and the type has no
/// conversion from a bare operation id.
#[derive(Debug)]
pub struct RecordedProductionRun {
    op: i64,
    run_receipt: fs_blake3::ContentHash,
    baseline_hash: Option<fs_blake3::ContentHash>,
    promotion_policy_receipt: Option<fs_blake3::ContentHash>,
    dependency_receipt_digest: Option<fs_blake3::ContentHash>,
    dependency_receipt_artifact: Option<fs_blake3::ContentHash>,
    recorded_at_ns: i64,
    admitted: bool,
}

impl RecordedProductionRun {
    /// Ledger operation id, retained for diagnostics and lineage only.
    #[must_use]
    pub const fn op_id(&self) -> i64 {
        self.op
    }

    /// Finalized identity of the exact baseline receipt and ordered results.
    #[must_use]
    pub const fn run_receipt(&self) -> fs_blake3::ContentHash {
        self.run_receipt
    }

    /// Exact baseline selected by the frozen admission, when present.
    #[must_use]
    pub const fn baseline_hash(&self) -> Option<fs_blake3::ContentHash> {
        self.baseline_hash
    }

    /// Exact promotion-policy fingerprint frozen into admission.
    #[must_use]
    pub const fn baseline_authority_fingerprint(&self) -> Option<fs_blake3::ContentHash> {
        self.promotion_policy_receipt
    }

    /// Exact dependency receipt digest retained by the operation. This is
    /// evidence identity, not the live authority-policy identity.
    #[must_use]
    pub const fn dependency_receipt_digest(&self) -> Option<fs_blake3::ContentHash> {
        self.dependency_receipt_digest
    }

    /// Operation completion time committed in the ledger.
    #[must_use]
    pub const fn recorded_at_ns(&self) -> i64 {
        self.recorded_at_ns
    }

    /// Reconstruct a typed receipt for an already-recorded, successfully
    /// admitted sealed production operation. Refusal operations remain
    /// durable diagnostics but are intentionally outside this positive-capable
    /// loader.
    ///
    /// # Errors
    /// Returns a typed corruption refusal when the operation or any exact
    /// manifest member fails authentication, and propagates ledger failures.
    pub fn load(ledger: &Ledger, op: i64) -> Result<Self, ProductionFreshnessError> {
        let Some(validated) = crate::validate_recorded_production_run(ledger, op, None)? else {
            return Err(ProductionFreshnessError::CorruptRecordedEvidence);
        };
        Ok(Self {
            op: validated.op,
            run_receipt: validated.run_receipt,
            baseline_hash: Some(validated.baseline.content_hash()),
            promotion_policy_receipt: Some(validated.promotion_policy_receipt),
            dependency_receipt_digest: Some(validated.dependency_receipt_digest),
            dependency_receipt_artifact: Some(validated.dependency_receipt_artifact),
            recorded_at_ns: validated.recorded_at_ns,
            admitted: true,
        })
    }

    /// Authenticate this exact operation, then perform one named live
    /// promotion-authority recheck and mint the only positive evidence type.
    ///
    /// # Errors
    /// Returns a specific freshness refusal for revocation, authority or
    /// dependency drift, tamper, build drift, rollback, or expiry; ledger
    /// failures are preserved as [`ProductionFreshnessError::Ledger`].
    pub fn revalidate(
        &self,
        ledger: &Ledger,
        current: &ProductionFreshnessContext<'_>,
    ) -> Result<FreshProductionEvidence, ProductionFreshnessError> {
        if !self.admitted {
            return Err(ProductionFreshnessError::RecordedRefusal);
        }
        let dependency = DependencyReceiptBinding::current()
            .map_err(|_| ProductionFreshnessError::DependencyAuthorityUnavailable)?;
        self.revalidate_at_with_dependency(ledger, current, fs_ledger::now_wall_ns(), dependency)
    }

    fn revalidate_at_with_dependency(
        &self,
        ledger: &Ledger,
        current: &ProductionFreshnessContext<'_>,
        observed_ns: i64,
        dependency: DependencyReceiptBinding,
    ) -> Result<FreshProductionEvidence, ProductionFreshnessError> {
        if !self.admitted {
            return Err(ProductionFreshnessError::RecordedRefusal);
        }
        if self.dependency_receipt_digest != Some(dependency.domain_digest)
            || self.dependency_receipt_artifact != Some(dependency.artifact_hash)
        {
            return Err(ProductionFreshnessError::DependencyAuthorityChanged);
        }
        let Some(validated) =
            crate::validate_recorded_production_run(ledger, self.op, Some(dependency))?
        else {
            return Err(ProductionFreshnessError::CorruptRecordedEvidence);
        };
        if validated.op != self.op
            || validated.run_receipt != self.run_receipt
            || Some(validated.baseline.content_hash()) != self.baseline_hash
            || Some(validated.promotion_policy_receipt) != self.promotion_policy_receipt
            || Some(validated.dependency_receipt_digest) != self.dependency_receipt_digest
            || Some(validated.dependency_receipt_artifact) != self.dependency_receipt_artifact
            || validated.recorded_at_ns != self.recorded_at_ns
        {
            return Err(ProductionFreshnessError::CorruptRecordedEvidence);
        }
        let dependency_decision = current.dependency_authority.verify(
            validated.dependency_receipt_digest,
            validated.dependency_receipt_artifact,
        );
        if dependency_decision.verdict() != DependencyReceiptVerdict::Authorized {
            return Err(ProductionFreshnessError::DependencyAuthorityRevoked);
        }
        let current_build = crate::executable_build_identity()?;
        if validated.build_identity != current_build {
            return Err(ProductionFreshnessError::BuildDrift);
        }
        let fingerprint = validated.baseline.identity().fingerprint();
        let Some(current_baseline) = current.baselines.for_fingerprint(fingerprint) else {
            return Err(ProductionFreshnessError::BaselineUnavailable);
        };
        let Some(current_attestation) = current.baselines.attestation_for(fingerprint) else {
            return Err(ProductionFreshnessError::PromotionAttestationChanged);
        };
        if current_baseline != &validated.baseline {
            return Err(ProductionFreshnessError::BaselineReplaced);
        }
        if current_attestation != &validated.attestation {
            return Err(ProductionFreshnessError::PromotionAttestationChanged);
        }
        if let Some(missing) = validated
            .baseline
            .provenance()
            .source_receipts()
            .iter()
            .find(|receipt| !current.retained_sources.contains(receipt))
        {
            return Err(ProductionFreshnessError::SourceReceiptUnavailable(*missing));
        }
        let decision = validated
            .baseline
            .authority_verdict(Some(&validated.attestation), current.authority);
        if decision.verdict() != crate::KeyVerdict::Authorized {
            return Err(ProductionFreshnessError::PromotionAuthorityRefused {
                verdict: decision.verdict(),
            });
        }
        if decision.policy_receipt() != validated.promotion_policy_receipt {
            return Err(ProductionFreshnessError::PromotionPolicyChanged);
        }
        if observed_ns < validated.recorded_at_ns {
            return Err(ProductionFreshnessError::ClockRollback {
                recorded_at_ns: validated.recorded_at_ns,
                observed_ns,
            });
        }
        let age_ns = observed_ns.saturating_sub(validated.recorded_at_ns);
        if age_ns > crate::STALENESS_MAX_AGE_NS {
            return Err(ProductionFreshnessError::Expired { age_ns });
        }
        Ok(FreshProductionEvidence {
            op: validated.op,
            run_receipt: validated.run_receipt,
            baseline_authority_fingerprint: validated.promotion_policy_receipt,
            dependency_authority_fingerprint: dependency_decision.policy_receipt(),
            recorded_at_ns: validated.recorded_at_ns,
            revalidated_at_ns: observed_ns,
        })
    }

    #[cfg(test)]
    fn revalidate_at_for_test(
        &self,
        ledger: &Ledger,
        current: &ProductionFreshnessContext<'_>,
        observed_ns: i64,
        dependency: DependencyReceiptBinding,
    ) -> Result<FreshProductionEvidence, ProductionFreshnessError> {
        self.revalidate_at_with_dependency(ledger, current, observed_ns, dependency)
    }
}

/// Positive exact-operation production evidence that was fresh at
/// [`Self::revalidated_at_ns`]. This is a point-in-time receipt, not a lease:
/// later citation must revalidate again so intervening revocation or expiry is
/// observed. There is no public constructor and the type is deliberately
/// neither `Clone` nor `Copy`.
/// A committed operation id (or even a [`RecordedProductionRun`]) cannot stand
/// in for the positive type:
///
/// ```compile_fail
/// use fs_roofline::production::FreshProductionEvidence;
/// let bare_operation_id: i64 = 7;
/// let _: FreshProductionEvidence = bare_operation_id;
/// ```
///
/// ```compile_fail
/// use fs_roofline::production::{FreshProductionEvidence, RecordedProductionRun};
/// fn cite(_: &FreshProductionEvidence) {}
/// fn recorded_is_not_fresh(recorded: &RecordedProductionRun) {
///     cite(recorded);
/// }
/// ```
#[derive(Debug)]
pub struct FreshProductionEvidence {
    op: i64,
    run_receipt: fs_blake3::ContentHash,
    baseline_authority_fingerprint: fs_blake3::ContentHash,
    dependency_authority_fingerprint: fs_blake3::ContentHash,
    recorded_at_ns: i64,
    revalidated_at_ns: i64,
}

impl FreshProductionEvidence {
    /// Exact successful production operation proved by this value.
    #[must_use]
    pub const fn op_id(&self) -> i64 {
        self.op
    }

    /// Finalized identity of the exact baseline and ordered result manifest.
    #[must_use]
    pub const fn run_receipt(&self) -> fs_blake3::ContentHash {
        self.run_receipt
    }

    /// Promotion-policy fingerprint sampled by the live revalidation.
    #[must_use]
    pub const fn baseline_authority_fingerprint(&self) -> fs_blake3::ContentHash {
        self.baseline_authority_fingerprint
    }

    /// Dependency-policy fingerprint sampled by the live revalidation.
    #[must_use]
    pub const fn dependency_authority_fingerprint(&self) -> fs_blake3::ContentHash {
        self.dependency_authority_fingerprint
    }

    /// Ledgered completion time authenticated from the exact operation.
    #[must_use]
    pub const fn recorded_at_ns(&self) -> i64 {
        self.recorded_at_ns
    }

    /// Live wall-clock observation used for rollback and expiry admission.
    #[must_use]
    pub const fn revalidated_at_ns(&self) -> i64 {
        self.revalidated_at_ns
    }
}

/// Sizing and repetition parameters for one production run.
#[derive(Debug, Clone, Copy)]
pub struct ProductionRunConfig {
    /// Vector-kernel element count (GEMM derives its edge from this).
    pub n: usize,
    /// Untimed warmup repetitions per kernel.
    pub warmup: usize,
    /// Timed repetitions per kernel.
    pub reps: usize,
}

/// Largest vector length accepted by the sealed production runner.
///
/// The shipped registry holds several `f64` vectors plus three approximately
/// `n`-element GEMM matrices, so this bounds aggregate allocation as well as
/// each individual kernel input.
pub const MAX_PRODUCTION_ELEMENTS: usize = crate::kernels::MAX_VECTOR_KERNEL_ELEMENTS;
/// Largest total warmup plus timed invocation count for each production kernel.
///
/// This is intentionally much tighter than the generic measurement API's
/// per-field limits. It bounds dispatch/receipt overhead, cancellation latency,
/// and repeated minimum-shape amplification while leaving ample room above the
/// shipped default's eleven runs.
pub const MAX_PRODUCTION_KERNEL_RUNS: usize = 64;
/// Largest untimed warmup count accepted by the sealed production runner.
///
/// At least one timed repetition is mandatory, so warmups alone can occupy at
/// most 63 of the 64 per-kernel invocation slots.
pub const MAX_PRODUCTION_WARMUP: usize = MAX_PRODUCTION_KERNEL_RUNS - 1;
/// Largest timed repetition count accepted by the sealed production runner.
pub const MAX_PRODUCTION_REPS: usize = MAX_PRODUCTION_KERNEL_RUNS;
/// Largest modeled floating-point work admitted across the complete registry.
///
/// `2^39` admits the shipped profile's approximately 189 billion FLOPs while
/// limiting the maximum vector/GEMM shape to three complete registry runs.
pub const MAX_PRODUCTION_REGISTRY_FLOPS: u128 = 1 << 39;
/// Largest modeled logical byte traffic admitted across the complete registry.
///
/// `2^33` admits the shipped profile's approximately 3.3 GiB of declared
/// traffic while independently bounding bandwidth-heavy repetition profiles.
pub const MAX_PRODUCTION_REGISTRY_BYTES: u128 = 1 << 33;

impl ProductionRunConfig {
    /// Validate resource-driving inputs before registry allocation or timing.
    ///
    /// # Errors
    /// Returns a stable diagnostic for zero or out-of-envelope inputs.
    pub fn validate(self) -> Result<(), String> {
        if self.n == 0 || self.n > MAX_PRODUCTION_ELEMENTS {
            return Err(format!(
                "production n must be in 1..={MAX_PRODUCTION_ELEMENTS}, got {}",
                self.n
            ));
        }
        if self.warmup > MAX_PRODUCTION_WARMUP {
            return Err(format!(
                "production warmup must be in 0..={MAX_PRODUCTION_WARMUP}, got {}",
                self.warmup
            ));
        }
        if self.reps == 0 || self.reps > MAX_PRODUCTION_REPS {
            return Err(format!(
                "production reps must be in 1..={MAX_PRODUCTION_REPS}, got {}",
                self.reps
            ));
        }
        let runs_per_kernel = self
            .warmup
            .checked_add(self.reps)
            .ok_or_else(|| "production warmup + repetition count overflowed usize".to_string())?;
        if runs_per_kernel > MAX_PRODUCTION_KERNEL_RUNS {
            return Err(format!(
                "production warmup + reps must be at most {MAX_PRODUCTION_KERNEL_RUNS} runs per kernel, got {runs_per_kernel}"
            ));
        }
        let work = crate::kernels::production_registry_work(self.n, runs_per_kernel)?;
        debug_assert_eq!(work.runs_per_kernel, runs_per_kernel);
        if work.total_flops > MAX_PRODUCTION_REGISTRY_FLOPS {
            return Err(format!(
                "production registry requires {} modeled FLOPs, exceeding the {MAX_PRODUCTION_REGISTRY_FLOPS}-FLOP bound",
                work.total_flops
            ));
        }
        if work.total_bytes > MAX_PRODUCTION_REGISTRY_BYTES {
            return Err(format!(
                "production registry requires {} modeled logical bytes, exceeding the {MAX_PRODUCTION_REGISTRY_BYTES}-byte bound",
                work.total_bytes
            ));
        }
        Ok(())
    }
}

/// Stage one of the sealed protocol: a pre-run axis probe this crate
/// performed itself, plus the minted per-run nonce.
///
/// No public constructor accepts axes; the only public way in is
/// [`ProductionProbe::observe`], which probes the actual machine.
pub struct ProductionProbe {
    axes: MachineAxes,
    nonce: fs_blake3::ContentHash,
}

impl ProductionProbe {
    /// Probe the machine and mint this run's nonce.
    #[must_use]
    pub fn observe() -> Self {
        Self::from_observed(MachineAxes::probe())
    }

    /// Test seam (`pub(crate)`): inject a synthetic pre-probe. Unreachable
    /// by API consumers, so a forged probe cannot enter the protocol.
    pub(crate) fn from_observed(axes: MachineAxes) -> Self {
        let nonce = mint_nonce(&axes);
        Self { axes, nonce }
    }

    /// The observed pre-run axes (read-only; baseline selection needs them).
    #[must_use]
    pub fn axes(&self) -> &MachineAxes {
        &self.axes
    }

    /// Run the production registry and finalize, consuming the probe.
    ///
    /// The tune ledger (optional) lets the GEMM kernel adopt a previously
    /// validated row; the registry (and with it fsqlite's `!Send`
    /// connection) is dropped before this returns, so the caller may reopen
    /// the same database for [`ProductionRun::record`].
    ///
    /// # Errors
    /// Structured diagnostics from tuning finalization; admission refusal is
    /// NOT an error — the run comes back with `citation_eligible() == false`
    /// and can be recorded as an explicit rejection.
    ///
    /// An operator-trusted candidate policy is rejected by the type system:
    ///
    /// ```compile_fail
    /// use fs_roofline::{AxisBaselinePolicy, BaselineIdentity};
    /// use fs_roofline::production::{ProductionProbe, ProductionRunConfig};
    ///
    /// let probe = ProductionProbe::observe();
    /// let identity = BaselineIdentity::current(probe.axes(), "declared-firmware")
    ///     .expect("valid declared identity");
    /// let candidate = AxisBaselinePolicy::new(None, &identity, 0);
    /// let config = ProductionRunConfig { n: 1, warmup: 0, reps: 1 };
    /// let _ = probe.run(config, candidate, None).expect("type-sealed run");
    /// ```
    pub fn run(
        self,
        config: ProductionRunConfig,
        baseline: AttestedAxisBaselinePolicy,
        tune_ledger: Option<Ledger>,
    ) -> Result<ProductionRun, String> {
        config.validate()?;
        let registry = production_registry_with_ledger(config.n, &self.axes, tune_ledger)?;
        self.run_with_parts(config, baseline, registry, MachineAxes::probe)
    }

    /// Measure the sealed registry with an operator-trusted baseline while
    /// making the non-citable boundary explicit in the return type.
    ///
    /// The returned value deliberately has no `citation_eligible` method and
    /// records only candidate/report-only evidence with a structured refusal.
    pub fn run_report_only(
        self,
        config: ProductionRunConfig,
        baseline: AxisBaselinePolicy<'_>,
        tune_ledger: Option<Ledger>,
    ) -> Result<ReportOnlyProductionRun, String> {
        config.validate()?;
        let registry = production_registry_with_ledger(config.n, &self.axes, tune_ledger)?;
        self.run_report_only_with_parts(config, baseline, registry, MachineAxes::probe)
    }

    /// Protocol core with injected registry and post-probe (`pub(crate)`
    /// test seam: drifted-post and finalizer-failure paths need determinism;
    /// API consumers cannot reach this to forge a run).
    pub(crate) fn run_with_parts(
        self,
        config: ProductionRunConfig,
        baseline: AttestedAxisBaselinePolicy,
        registry: Vec<Box<dyn RooflineKernel>>,
        post_probe: impl FnOnce() -> MachineAxes,
    ) -> Result<ProductionRun, String> {
        self.run_with_parts_and_authority(
            config,
            baseline,
            registry,
            post_probe,
            CitationAuthority::from_build(),
        )
    }

    fn run_with_parts_and_authority(
        self,
        config: ProductionRunConfig,
        baseline: AttestedAxisBaselinePolicy,
        mut registry: Vec<Box<dyn RooflineKernel>>,
        post_probe: impl FnOnce() -> MachineAxes,
        citation_authority: CitationAuthority,
    ) -> Result<ProductionRun, String> {
        config.validate()?;
        let build_identity = crate::read_executable_build_identity().map_err(|error| {
            format!("cannot capture pre-measurement executable identity: {error}")
        })?;
        let results = run_registry(&mut registry, config.warmup, config.reps, &self.axes)?;
        let post_axes = post_probe();
        let admission = baseline.decide(&self.axes, &post_axes);
        let finalized = finalize_registry_tuning_with_snapshot(
            &mut registry,
            &self.axes,
            &admission,
            &results,
        )?;
        drop(registry);
        Ok(ProductionRun {
            axes: self.axes,
            post_axes,
            admission,
            nonce: self.nonce,
            results,
            finalized,
            citation_authority,
            build_identity,
        })
    }

    pub(crate) fn run_report_only_with_parts(
        self,
        config: ProductionRunConfig,
        baseline: AxisBaselinePolicy<'_>,
        mut registry: Vec<Box<dyn RooflineKernel>>,
        post_probe: impl FnOnce() -> MachineAxes,
    ) -> Result<ReportOnlyProductionRun, String> {
        config.validate()?;
        let build_identity = crate::read_executable_build_identity().map_err(|error| {
            format!("cannot capture pre-measurement executable identity: {error}")
        })?;
        let results = run_registry(&mut registry, config.warmup, config.reps, &self.axes)?;
        let post_axes = post_probe();
        let admission = baseline.snapshot(&self.axes, &post_axes);
        let finalized = finalize_registry_tuning_with_snapshot(
            &mut registry,
            &self.axes,
            &admission,
            &results,
        )?;
        drop(registry);
        Ok(ReportOnlyProductionRun {
            axes: self.axes,
            post_axes,
            admission,
            nonce: self.nonce,
            results,
            finalized,
            build_identity,
            report_only_refusal: REPORT_ONLY_REFUSAL.to_string(),
        })
    }

    #[cfg(test)]
    pub(crate) fn run_with_test_receipt(
        self,
        config: ProductionRunConfig,
        baseline: AttestedAxisBaselinePolicy,
        registry: Vec<Box<dyn RooflineKernel>>,
        post_probe: impl FnOnce() -> MachineAxes,
        receipt: &'static str,
    ) -> Result<ProductionRun, String> {
        let digest =
            fs_blake3::hash_domain(fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN, receipt.as_bytes());
        let binding = DependencyReceiptBinding::from_parts(receipt, digest)
            .expect("test receipt digest was computed from the same bytes");
        self.run_with_parts_and_authority(
            config,
            baseline,
            registry,
            post_probe,
            CitationAuthority::Receipt(binding),
        )
    }
}

/// One complete, sealed production registry run.
///
/// Fields are private, there is no public constructor, and the value is
/// neither `Clone` nor `Copy`: the only way to obtain one is
/// [`ProductionProbe::run`], which performed both probes and timed the
/// production registry itself. [`ProductionRun::record`] consumes the run.
pub struct ProductionRun {
    axes: MachineAxes,
    post_axes: MachineAxes,
    admission: AxisAdmissionSnapshot,
    nonce: fs_blake3::ContentHash,
    results: Vec<Attainment>,
    finalized: FinalizedRegistryRun,
    citation_authority: CitationAuthority,
    build_identity: fs_blake3::ContentHash,
}

impl std::fmt::Debug for ProductionRun {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ProductionRun")
            .field(
                "fingerprint",
                &format_args!("{:016x}", self.axes.fingerprint),
            )
            .field("kernels", &self.results.len())
            .field("nonce", &self.nonce)
            .field("citation_eligible", &self.citation_eligible())
            .finish_non_exhaustive()
    }
}

impl ProductionRun {
    /// The pre-run axis probe observed by the protocol.
    #[must_use]
    pub fn axes(&self) -> &MachineAxes {
        &self.axes
    }

    /// The post-run axis probe, observed strictly after the timed loop.
    #[must_use]
    pub fn post_axes(&self) -> &MachineAxes {
        &self.post_axes
    }

    /// The measured result set in registry order.
    #[must_use]
    pub fn results(&self) -> &[Attainment] {
        &self.results
    }

    /// The per-run nonce bound into the recorded operation.
    #[must_use]
    pub fn nonce(&self) -> fs_blake3::ContentHash {
        self.nonce
    }

    /// Whether this run passed aggregate admission and has production receipt
    /// provenance. This is a pre-commit eligibility predicate, not a claim
    /// that evidence was durably recorded or remains fresh.
    #[must_use]
    pub fn citation_eligible(&self) -> bool {
        self.admission_error().is_none()
    }

    /// Why admission refused this run, if it did.
    #[must_use]
    pub fn admission_error(&self) -> Option<String> {
        citable_run_admission_error_for_snapshot(&self.axes, &self.admission, &self.results)
            .or_else(|| {
                (!self.finalized.admitted()).then(|| SEALED_FINALIZATION_REFUSAL.to_string())
            })
            .or_else(|| self.citation_authority.refusal().map(str::to_string))
    }

    /// The baseline-admission receipt for this run's exact probe pair.
    #[must_use]
    pub fn receipt_json(&self) -> &str {
        self.admission.receipt_json()
    }

    /// Immutable admission decision retained by the run.
    #[must_use]
    pub fn admission_snapshot(&self) -> &AxisAdmissionSnapshot {
        &self.admission
    }

    /// Selected baseline identity, when present.
    #[must_use]
    pub fn baseline_hash(&self) -> Option<fs_blake3::ContentHash> {
        self.admission.baseline_hash()
    }

    /// Record the run atomically, consuming it, and return a neutral typed
    /// receipt. The operation `ir` carries
    /// `"protocol":"production-v3"`, the per-run nonce, content hashes of
    /// both observed axis receipts, and dependency-receipt provenance. This
    /// does not claim freshness; callers must explicitly invoke
    /// [`RecordedProductionRun::revalidate`].
    ///
    /// # Errors
    /// Ledger errors propagate and roll back the whole write set; the run is
    /// consumed either way (a failed transaction cannot be replayed into a
    /// different ledger with edited results).
    pub fn record(mut self, ledger: &Ledger) -> Result<RecordedProductionRun, LedgerError> {
        let run_receipt = self.finalized.receipt_identity();
        let baseline_hash = self.admission.baseline_hash();
        let promotion_policy_receipt = self.admission.authority_policy_receipt();
        let dependency = self.citation_authority.receipt();
        let protocol_fields = match self.citation_authority {
            CitationAuthority::Receipt(binding) => format!(
                "{PRODUCTION_PROTOCOL_FIELD},\"run_nonce\":\"{}\",\"pre_axes_receipt\":\"{}\",\"post_axes_receipt\":\"{}\",\"dependency_graph_evidence\":\"operator-observed-receipt\",\"dependency_receipt_digest\":\"{}\",\"dependency_receipt_artifact\":\"{}\"",
                self.nonce,
                axes_receipt(&self.axes),
                axes_receipt(&self.post_axes),
                binding.domain_digest,
                binding.artifact_hash,
            ),
            CitationAuthority::Refused(reason) => format!(
                "{PRODUCTION_PROTOCOL_FIELD},\"run_nonce\":\"{}\",\"pre_axes_receipt\":\"{}\",\"post_axes_receipt\":\"{}\",\"dependency_graph_evidence\":\"development-equivalence-salt\",\"citation_refusal\":\"{}\"",
                self.nonce,
                axes_receipt(&self.axes),
                axes_receipt(&self.post_axes),
                json_escape(reason),
            ),
        };
        let recorded = record_run_with_protocol(
            ledger,
            &self.axes,
            &self.post_axes,
            &self.admission,
            &mut self.finalized,
            &mut self.results,
            &protocol_fields,
            self.citation_authority.receipt(),
            self.citation_authority.refusal(),
            crate::EvidenceNamespace::Production,
            Some(self.build_identity),
        )?;
        Ok(RecordedProductionRun {
            op: recorded.op,
            run_receipt,
            baseline_hash,
            promotion_policy_receipt,
            dependency_receipt_digest: dependency.map(|binding| binding.domain_digest),
            dependency_receipt_artifact: dependency.map(|binding| binding.artifact_hash),
            recorded_at_ns: recorded.recorded_at_ns,
            admitted: recorded.admitted,
        })
    }
}

/// One measured sealed-registry run that intentionally carries only
/// candidate/report-only baseline trust.
///
/// This type has no `citation_eligible` API and cannot record into the
/// production evidence namespace.
///
/// ```compile_fail
/// use fs_roofline::{AxisBaselinePolicy, BaselineIdentity};
/// use fs_roofline::production::{ProductionProbe, ProductionRunConfig};
///
/// let probe = ProductionProbe::observe();
/// let identity = BaselineIdentity::current(probe.axes(), "declared-firmware")
///     .expect("valid declared identity");
/// let candidate = AxisBaselinePolicy::new(None, &identity, 0);
/// let config = ProductionRunConfig { n: 1, warmup: 0, reps: 1 };
/// let run = probe
///     .run_report_only(config, candidate, None)
///     .expect("report-only run");
/// let _ = run.citation_eligible();
/// ```
pub struct ReportOnlyProductionRun {
    axes: MachineAxes,
    post_axes: MachineAxes,
    admission: AxisAdmissionSnapshot,
    nonce: fs_blake3::ContentHash,
    results: Vec<Attainment>,
    finalized: FinalizedRegistryRun,
    build_identity: fs_blake3::ContentHash,
    report_only_refusal: String,
}

impl std::fmt::Debug for ReportOnlyProductionRun {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ReportOnlyProductionRun")
            .field(
                "fingerprint",
                &format_args!("{:016x}", self.axes.fingerprint),
            )
            .field("kernels", &self.results.len())
            .field("nonce", &self.nonce)
            .finish_non_exhaustive()
    }
}

impl ReportOnlyProductionRun {
    /// The pre-run axis probe observed by the sealed runner.
    #[must_use]
    pub fn axes(&self) -> &MachineAxes {
        &self.axes
    }

    /// The post-run axis probe observed after timing.
    #[must_use]
    pub fn post_axes(&self) -> &MachineAxes {
        &self.post_axes
    }

    /// Measured results in registry order.
    #[must_use]
    pub fn results(&self) -> &[Attainment] {
        &self.results
    }

    /// Per-run nonce bound into the candidate operation.
    #[must_use]
    pub fn nonce(&self) -> fs_blake3::ContentHash {
        self.nonce
    }

    /// Structured report-only refusal (or a more specific measurement
    /// refusal when the baseline or rows also fail).
    #[must_use]
    pub fn admission_error(&self) -> Option<String> {
        Some(
            run_admission_error_for_snapshot(&self.axes, &self.admission, &self.results)
                .map_or_else(
                    || self.report_only_refusal.clone(),
                    |measurement| format!("{}; {measurement}", self.report_only_refusal),
                ),
        )
    }

    /// Exact candidate-tier admission snapshot bytes.
    #[must_use]
    pub fn receipt_json(&self) -> &str {
        self.admission.receipt_json()
    }

    /// Immutable admission decision retained by the run.
    #[must_use]
    pub fn admission_snapshot(&self) -> &AxisAdmissionSnapshot {
        &self.admission
    }

    /// Selected baseline identity, when present.
    #[must_use]
    pub fn baseline_hash(&self) -> Option<fs_blake3::ContentHash> {
        self.admission.baseline_hash()
    }

    /// Retain the entrypoint refusal that forced this already-report-only
    /// run out of the attested path. Refusals beyond the durable diagnostic
    /// bound retain a UTF-8-safe prefix plus the original byte length and a
    /// domain-separated digest, so an already completed report-only
    /// measurement cannot be discarded merely because its source/path error
    /// was long. This can only specialize diagnostics on a type that has no
    /// production-evidence writer.
    pub fn with_configuration_refusal(mut self, refusal: String) -> Result<Self, String> {
        if refusal.is_empty() {
            return Err("report-only configuration refusal must be nonempty".to_string());
        }
        self.report_only_refusal = if refusal.len() <= MAX_REPORT_ONLY_REFUSAL_BYTES {
            refusal
        } else {
            let original_bytes = refusal.len();
            let digest =
                fs_blake3::hash_domain(REPORT_ONLY_REFUSAL_DIGEST_DOMAIN, refusal.as_bytes());
            let suffix = format!(
                "...[truncated; original_bytes={original_bytes}; full_refusal_digest={digest}]"
            );
            let prefix_limit = MAX_REPORT_ONLY_REFUSAL_BYTES
                .checked_sub(suffix.len())
                .expect("bounded refusal suffix fits its owner-local ceiling");
            let mut prefix_end = prefix_limit.min(refusal.len());
            while !refusal.is_char_boundary(prefix_end) {
                prefix_end -= 1;
            }
            let mut bounded = refusal[..prefix_end].to_string();
            bounded.push_str(&suffix);
            debug_assert!(bounded.len() <= MAX_REPORT_ONLY_REFUSAL_BYTES);
            bounded
        };
        Ok(self)
    }

    /// Record only a structured candidate/report-only rejection. No metrics,
    /// benchmark-result events, production tune rows, or dependency receipt
    /// artifacts are published.
    pub fn record(mut self, ledger: &Ledger) -> Result<i64, LedgerError> {
        let protocol_fields = format!(
            "{CUSTOM_REGISTRY_PROTOCOL_FIELD},\"run_nonce\":\"{}\",\"pre_axes_receipt\":\"{}\",\"post_axes_receipt\":\"{}\",\"citation_refusal\":\"{}\"",
            self.nonce,
            axes_receipt(&self.axes),
            axes_receipt(&self.post_axes),
            json_escape(&self.report_only_refusal),
        );
        record_run_with_protocol(
            ledger,
            &self.axes,
            &self.post_axes,
            &self.admission,
            &mut self.finalized,
            &mut self.results,
            &protocol_fields,
            None,
            Some(&self.report_only_refusal),
            crate::EvidenceNamespace::Custom,
            Some(self.build_identity),
        )
        .map(|recorded| recorded.op)
    }
}

fn production_axes_receipt_json(input: &ProductionAxesReceiptInput) -> String {
    format!(
        "{{\"fingerprint\":\"{:016x}\",\"cpu_brand\":\"{}\",\"logical_cpus\":{},\"bandwidth_single_bits\":\"{:016x}\",\"bandwidth_all_core_bits\":\"{:016x}\",\"peak_single_bits\":\"{:016x}\",\"peak_all_core_bits\":\"{:016x}\"}}",
        input.fingerprint,
        json_escape(&input.cpu_brand),
        input.logical_cpus,
        input.bandwidth_single_bits,
        input.bandwidth_all_core_bits,
        input.peak_single_bits,
        input.peak_all_core_bits,
    )
}

pub(crate) fn machine_axes_receipt_json(axes: &MachineAxes) -> String {
    production_axes_receipt_json(&ProductionAxesReceiptInput::from_axes(axes))
}

fn axes_receipt_with_domain(
    input: &ProductionAxesReceiptInput,
    domain: &str,
) -> fs_blake3::ContentHash {
    fs_blake3::hash_domain(domain, production_axes_receipt_json(input).as_bytes())
}

/// Content hash of one observed probe's canonical JSONL receipt.
fn axes_receipt(axes: &MachineAxes) -> fs_blake3::ContentHash {
    axes_receipt_with_domain(
        &ProductionAxesReceiptInput::from_axes(axes),
        PRODUCTION_AXES_RECEIPT_DOMAIN,
    )
}

#[allow(dead_code)]
fn production_axes_receipt_is_current(
    axes: &MachineAxes,
    retained: fs_blake3::ContentHash,
) -> bool {
    axes_receipt(axes) == retained
}

/// Process-unique per-run challenge: wall clock, pid, a monotone counter,
/// and the pre-probe receipt. Uniqueness, not secrecy — see the module docs.
fn mint_nonce(axes: &MachineAxes) -> fs_blake3::ContentHash {
    let count = NONCE_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let mut material = Vec::new();
    material.extend_from_slice(&fs_ledger::now_wall_ns().to_le_bytes());
    material.extend_from_slice(&u64::from(std::process::id()).to_le_bytes());
    material.extend_from_slice(&count.to_le_bytes());
    material.extend_from_slice(axes.to_jsonl().as_bytes());
    fs_blake3::hash_domain(RUN_NONCE_DOMAIN, &material)
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::collections::BTreeSet;
    use std::rc::Rc;

    use super::*;
    use crate::kernels::default_registry;
    use crate::{
        BaselineAxes, BaselineCandidate, BaselineIdentity, KernelSpec, STALENESS_MAX_AGE_NS,
        Staleness, TargetAxis, Threading, promote_baseline, roofline_machine_key, staleness,
        staleness_at,
    };

    fn synthetic_axes(fingerprint: u64) -> MachineAxes {
        // Roofs far above any real machine (bead xjhz): cache-resident test
        // kernels must never outrun the fixture roof.
        MachineAxes {
            fingerprint,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 100_000.0,
            bandwidth_all_core_gbs: 400_000.0,
            peak_single_gflops: 50_000.0,
            peak_all_core_gflops: 300_000.0,
        }
    }

    #[test]
    fn production_axes_receipt_identity_fields_move_independently() {
        fn assert_moves(
            original: fs_blake3::ContentHash,
            altered: &ProductionAxesReceiptInput,
            field: &str,
        ) {
            assert_ne!(
                original,
                axes_receipt_with_domain(altered, PRODUCTION_AXES_RECEIPT_DOMAIN),
                "mutating {field} must move the production-axis receipt"
            );
        }

        let axes = synthetic_axes(0xA11CE);
        let input = ProductionAxesReceiptInput::from_axes(&axes);
        let original = axes_receipt(&axes);
        assert!(production_axes_receipt_is_current(&axes, original));
        assert_ne!(
            original,
            axes_receipt_with_domain(
                &input,
                "org.frankensim.fs-roofline.production-axes-receipt-foreign.v1",
            ),
            "the digest domain is semantic"
        );

        let mut altered = input.clone();
        altered.fingerprint += 1;
        assert_moves(original, &altered, "machine-fingerprint");
        let mut altered = input.clone();
        altered.cpu_brand.push('x');
        assert_moves(original, &altered, "cpu-brand-utf8");
        let mut altered = input.clone();
        altered.logical_cpus += 1;
        assert_moves(original, &altered, "logical-cpus");
        let mut altered = input.clone();
        altered.bandwidth_single_bits ^= 1;
        assert_moves(original, &altered, "bandwidth-single-bits");
        let mut altered = input.clone();
        altered.bandwidth_all_core_bits ^= 1;
        assert_moves(original, &altered, "bandwidth-all-core-bits");
        let mut altered = input.clone();
        altered.peak_single_bits ^= 1;
        assert_moves(original, &altered, "peak-single-bits");
        let mut altered = input;
        altered.peak_all_core_bits ^= 1;
        assert_moves(original, &altered, "peak-all-core-bits");
    }

    #[test]
    fn production_axes_receipt_versions_fail_closed() {
        assert_eq!(PRODUCTION_AXES_RECEIPT_IDENTITY_VERSION, 1);
        assert!(PRODUCTION_AXES_RECEIPT_DOMAIN.ends_with(".v1"));
        let axes = synthetic_axes(0xA11CE);
        let stale = axes_receipt_with_domain(
            &ProductionAxesReceiptInput::from_axes(&axes),
            "org.frankensim.fs-roofline.production-axes-receipt.v2",
        );
        assert!(
            !production_axes_receipt_is_current(&axes, stale),
            "a receipt under a stale or future identity domain must not be admitted"
        );
    }

    #[test]
    fn dependency_authority_policy_identity_fields_move_independently() {
        let input = DependencyAuthorityPolicyIdentityInput {
            canonical_bytes: b"",
        };
        let current = dependency_authority_policy_receipt(input.canonical_bytes);
        let changed_domain = dependency_authority_policy_receipt_with_domain(
            "frankensim.fs-roofline.dependency-authority-policy-shadow.v1",
            &input,
        );
        let changed_bytes = dependency_authority_policy_receipt(b"builtin\trevoked\n");
        assert_ne!(current, changed_domain, "the digest domain is semantic");
        assert_ne!(
            current, changed_bytes,
            "the exact policy bytes are semantic"
        );
    }

    #[test]
    fn dependency_authority_policy_identity_versions_fail_closed() {
        assert_eq!(DEPENDENCY_AUTHORITY_POLICY_IDENTITY_VERSION, 1);
        assert!(DEPENDENCY_AUTHORITY_POLICY_DOMAIN.ends_with(".v1"));
        let input = DependencyAuthorityPolicyIdentityInput {
            canonical_bytes: b"",
        };
        let current = dependency_authority_policy_receipt(input.canonical_bytes);
        let future = dependency_authority_policy_receipt_with_domain(
            "frankensim.fs-roofline.dependency-authority-policy.v2",
            &input,
        );
        let configured = ConfiguredDependencyReceiptAuthority::from_text("")
            .expect("empty no-revocations policy is canonical");
        assert_eq!(configured.policy_receipt(), current);
        assert_ne!(configured.policy_receipt(), future);
    }

    fn trusted_baseline(axes: &MachineAxes) -> (BaselineAxes, BaselineIdentity) {
        let identity =
            BaselineIdentity::current(axes, "test-firmware").expect("valid synthetic identity");
        let candidates: Vec<_> = (0_u64..3)
            .map(|ordinal| {
                BaselineCandidate::from_receipt(
                    axes.clone(),
                    identity.clone(),
                    fs_blake3::hash_domain(
                        "fs-roofline.production-baseline-source.v1",
                        &ordinal.to_le_bytes(),
                    ),
                )
                .expect("valid synthetic candidate")
            })
            .collect();
        let baseline = promote_baseline(
            &candidates,
            "test-operator",
            "deterministic production-protocol fixture",
            test_day().saturating_sub(10),
            90,
        )
        .expect("valid synthetic baseline");
        (baseline, identity)
    }

    fn attested_policy(
        baseline: &BaselineAxes,
        identity: &BaselineIdentity,
    ) -> AttestedAxisBaselinePolicy {
        attested_policy_on_day(baseline, identity, test_day())
    }

    fn attested_policy_on_day(
        baseline: &BaselineAxes,
        identity: &BaselineIdentity,
        now_day: u64,
    ) -> AttestedAxisBaselinePolicy {
        let policy_receipt = test_promotion_policy_receipt(baseline);
        AttestedAxisBaselinePolicy::from_verified(
            baseline.clone(),
            identity.clone(),
            now_day,
            crate::PromotionAttestation::new("test-authority", "test-signature"),
            baseline.provenance().source_receipts().to_vec(),
            crate::PromotionAuthorityDecision::new(crate::KeyVerdict::Authorized, policy_receipt),
        )
    }

    fn test_promotion_policy_receipt(baseline: &BaselineAxes) -> fs_blake3::ContentHash {
        fs_blake3::hash_domain(
            "fs-roofline.test-promotion-policy.v1",
            baseline.content_hash().as_bytes(),
        )
    }

    struct LiveTestAuthority {
        calls: Cell<usize>,
        verdict: crate::KeyVerdict,
        policy_receipt: fs_blake3::ContentHash,
    }

    struct FlippingTestAuthority {
        calls: Cell<usize>,
        admitted_policy: fs_blake3::ContentHash,
        revoked_policy: fs_blake3::ContentHash,
    }

    impl crate::PromotionAuthorityVerifier for FlippingTestAuthority {
        fn verify(
            &self,
            _key_id: &str,
            _signature: &str,
            _message: &[u8],
        ) -> crate::PromotionAuthorityDecision {
            let call = self.calls.get();
            self.calls.set(call + 1);
            if call == 0 {
                crate::PromotionAuthorityDecision::new(
                    crate::KeyVerdict::Authorized,
                    self.admitted_policy,
                )
            } else {
                crate::PromotionAuthorityDecision::new(
                    crate::KeyVerdict::RevokedKey,
                    self.revoked_policy,
                )
            }
        }
    }

    impl LiveTestAuthority {
        fn new(verdict: crate::KeyVerdict, policy_receipt: fs_blake3::ContentHash) -> Self {
            Self {
                calls: Cell::new(0),
                verdict,
                policy_receipt,
            }
        }
    }

    impl crate::PromotionAuthorityVerifier for LiveTestAuthority {
        fn verify(
            &self,
            _key_id: &str,
            _signature: &str,
            _message: &[u8],
        ) -> crate::PromotionAuthorityDecision {
            self.calls.set(self.calls.get() + 1);
            crate::PromotionAuthorityDecision::new(self.verdict, self.policy_receipt)
        }
    }

    fn live_baseline_store(
        baseline: &BaselineAxes,
    ) -> (
        crate::AttestedBaselineStore,
        BTreeSet<fs_blake3::ContentHash>,
    ) {
        let retained: BTreeSet<_> = baseline
            .provenance()
            .source_receipts()
            .iter()
            .copied()
            .collect();
        let authority = LiveTestAuthority::new(
            crate::KeyVerdict::Authorized,
            test_promotion_policy_receipt(baseline),
        );
        let mut store = crate::AttestedBaselineStore::new();
        store
            .admit_verified(
                baseline.clone(),
                crate::PromotionAttestation::new("test-authority", "test-signature"),
                &authority,
                &retained,
            )
            .expect("authorized live baseline fixture");
        (store, retained)
    }

    fn test_day() -> u64 {
        crate::days_since_epoch_now().expect("unit-test clock after Unix epoch")
    }

    fn temp_db(tag: &str) -> String {
        static NEXT: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        std::env::temp_dir()
            .join(format!(
                "fs-roofline-prod-{tag}-{}-{n}.db",
                std::process::id()
            ))
            .display()
            .to_string()
    }

    const CONFIG: ProductionRunConfig = ProductionRunConfig {
        n: 1 << 10,
        warmup: 0,
        reps: 1,
    };

    const TEST_DEPGRAPH_RECEIPT: &str = "{\"schema\":\"fs-roofline-synthetic-dependency-receipt-v1\",\"purpose\":\"unit-test-only\"}";
    const SUBSTITUTE_DEPGRAPH_RECEIPT: &str = "{\"schema\":\"fs-roofline-synthetic-dependency-receipt-v1\",\"purpose\":\"substitution-attack\"}";
    static TEST_DEPENDENCY_AUTHORITY: TestDependencyAuthority = TestDependencyAuthority;

    struct TestDependencyAuthority;

    impl DependencyReceiptAuthority for TestDependencyAuthority {
        fn verify(
            &self,
            _digest: fs_blake3::ContentHash,
            _artifact: fs_blake3::ContentHash,
        ) -> DependencyReceiptDecision {
            DependencyReceiptDecision::new(
                DependencyReceiptVerdict::Authorized,
                dependency_authority_policy_receipt(b"test\tallow-all\n"),
            )
        }
    }

    struct LiveTestDependencyAuthority {
        calls: Cell<usize>,
        verdict: DependencyReceiptVerdict,
        policy_receipt: fs_blake3::ContentHash,
    }

    impl LiveTestDependencyAuthority {
        fn new(verdict: DependencyReceiptVerdict, policy_receipt: fs_blake3::ContentHash) -> Self {
            Self {
                calls: Cell::new(0),
                verdict,
                policy_receipt,
            }
        }
    }

    impl DependencyReceiptAuthority for LiveTestDependencyAuthority {
        fn verify(
            &self,
            _digest: fs_blake3::ContentHash,
            _artifact: fs_blake3::ContentHash,
        ) -> DependencyReceiptDecision {
            self.calls.set(self.calls.get() + 1);
            DependencyReceiptDecision::new(self.verdict, self.policy_receipt)
        }
    }

    fn test_receipt_binding() -> DependencyReceiptBinding {
        let digest = fs_blake3::hash_domain(
            fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN,
            TEST_DEPGRAPH_RECEIPT.as_bytes(),
        );
        DependencyReceiptBinding::from_parts(TEST_DEPGRAPH_RECEIPT, digest)
            .expect("test receipt digest agrees")
    }

    fn test_receipt_authority() -> CitationAuthority {
        CitationAuthority::Receipt(test_receipt_binding())
    }

    #[allow(clippy::too_many_arguments)]
    fn receipt_staleness_at(
        ledger: &Ledger,
        kernel: &str,
        version: &str,
        fingerprint: u64,
        baseline: Option<fs_blake3::ContentHash>,
        observed_wall_ns: i64,
        dependency: DependencyReceiptBinding,
    ) -> Result<Staleness, fs_ledger::LedgerError> {
        crate::staleness_at_with_build_and_dependency(
            ledger,
            kernel,
            version,
            fingerprint,
            baseline,
            observed_wall_ns,
            crate::executable_build_identity()?,
            Some(dependency),
        )
    }

    struct CountingKernel {
        runs: Rc<Cell<usize>>,
        value: u64,
    }

    impl crate::RooflineKernel for CountingKernel {
        fn spec(&self) -> KernelSpec {
            KernelSpec {
                name: "counting-kernel",
                version: "1",
                bytes_per_elem: 8.0,
                flops_per_elem: 1.0,
                threading: Threading::SingleThread,
                target_axis: TargetAxis::BindingRoof,
                target_fraction: None,
            }
        }

        fn elements(&self) -> usize {
            64
        }

        fn run_once(&mut self) -> Result<(), String> {
            self.runs.set(self.runs.get() + 1);
            for _ in 0..64 {
                self.value = std::hint::black_box(
                    self.value
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1),
                );
            }
            Ok(())
        }
    }

    struct FailingFinalizeKernel {
        inner: CountingKernel,
    }

    impl crate::RooflineKernel for FailingFinalizeKernel {
        fn spec(&self) -> KernelSpec {
            self.inner.spec()
        }
        fn elements(&self) -> usize {
            self.inner.elements()
        }
        fn run_once(&mut self) -> Result<(), String> {
            self.inner.run_once()
        }
        fn finalize_tuning(&mut self, _admitted: bool) -> Result<(), String> {
            Err("tune ledger unavailable mid-finalize".to_string())
        }
    }

    #[test]
    fn nonces_are_unique_per_probe() {
        let a = ProductionProbe::from_observed(synthetic_axes(0xA));
        let b = ProductionProbe::from_observed(synthetic_axes(0xA));
        assert_ne!(
            a.nonce, b.nonce,
            "identical axes must still mint distinct nonces"
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One exact/limit-plus-one work-admission matrix.
    fn production_config_rejects_zero_and_unbounded_work_before_running() {
        for (config, expected) in [
            (
                ProductionRunConfig {
                    n: 0,
                    warmup: 0,
                    reps: 1,
                },
                "production n",
            ),
            (
                ProductionRunConfig {
                    n: MAX_PRODUCTION_ELEMENTS + 1,
                    warmup: 0,
                    reps: 1,
                },
                "production n",
            ),
            (
                ProductionRunConfig {
                    n: 1,
                    warmup: MAX_PRODUCTION_WARMUP + 1,
                    reps: 1,
                },
                "production warmup",
            ),
            (
                ProductionRunConfig {
                    n: 1,
                    warmup: 0,
                    reps: 0,
                },
                "production reps",
            ),
            (
                ProductionRunConfig {
                    n: 1,
                    warmup: 0,
                    reps: MAX_PRODUCTION_REPS + 1,
                },
                "production reps",
            ),
            (
                ProductionRunConfig {
                    n: MAX_PRODUCTION_ELEMENTS,
                    warmup: MAX_PRODUCTION_WARMUP,
                    reps: 2,
                },
                "warmup + reps",
            ),
            (
                ProductionRunConfig {
                    n: 1,
                    warmup: MAX_PRODUCTION_WARMUP,
                    reps: MAX_PRODUCTION_REPS,
                },
                "warmup + reps",
            ),
        ] {
            let error = config.validate().expect_err("invalid config must fail");
            assert!(error.contains(expected), "unexpected diagnostic: {error}");
        }
        ProductionRunConfig {
            n: MAX_PRODUCTION_ELEMENTS,
            warmup: 0,
            reps: 1,
        }
        .validate()
        .expect("maximum allocation with one pass is admitted without allocating");

        ProductionRunConfig {
            n: 1,
            warmup: 0,
            reps: MAX_PRODUCTION_KERNEL_RUNS,
        }
        .validate()
        .expect("the exact per-kernel run cap is admitted");
        let error = ProductionRunConfig {
            n: 1,
            warmup: 1,
            reps: MAX_PRODUCTION_KERNEL_RUNS,
        }
        .validate()
        .expect_err("one run beyond the per-kernel cap must be refused");
        assert!(
            error.contains("warmup + reps"),
            "unexpected diagnostic: {error}"
        );

        ProductionRunConfig {
            n: 1 << 22,
            warmup: 2,
            reps: 9,
        }
        .validate()
        .expect("the shipped default profile remains inside every derived-work cap");

        ProductionRunConfig {
            n: MAX_PRODUCTION_ELEMENTS,
            warmup: 2,
            reps: 1,
        }
        .validate()
        .expect("three maximum-shape runs remain below the modeled FLOP cap");
        let error = ProductionRunConfig {
            n: MAX_PRODUCTION_ELEMENTS,
            warmup: 3,
            reps: 1,
        }
        .validate()
        .expect_err("four maximum-shape runs must exceed the modeled FLOP cap");
        assert!(
            error.contains("modeled FLOPs"),
            "unexpected diagnostic: {error}"
        );

        ProductionRunConfig {
            n: 1 << 22,
            warmup: 27,
            reps: 1,
        }
        .validate()
        .expect("twenty-eight default-shape runs remain below the byte cap");
        let error = ProductionRunConfig {
            n: 1 << 22,
            warmup: 28,
            reps: 1,
        }
        .validate()
        .expect_err("twenty-nine default-shape runs must exceed the byte cap");
        assert!(
            error.contains("modeled logical bytes"),
            "unexpected diagnostic: {error}"
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one complete retained-receipt lineage attack matrix
    fn dependency_receipt_verifier_requires_input_edge_and_intact_bytes() {
        let CitationAuthority::Receipt(binding) = test_receipt_authority() else {
            unreachable!("test helper is receipt-backed")
        };
        let db = temp_db("dependency-verifier");
        let ledger = Ledger::open(&db).expect("open ledger");
        let stored = ledger
            .put_artifact(
                crate::DEPGRAPH_RECEIPT_ARTIFACT_KIND,
                binding.bytes.as_bytes(),
                Some(crate::DEPGRAPH_RECEIPT_ARTIFACT_META),
            )
            .expect("store receipt");
        assert_eq!(stored.hash, binding.artifact_hash);
        let explicits = fs_ledger::FiveExplicits {
            seed: b"dependency-verifier",
            versions: "{}",
            budget: "{}",
            capability: "{}",
        };
        let op = ledger
            .begin_op(None, "{}", &explicits, 1)
            .expect("begin verifier fixture");
        let placeholder = fs_blake3::hash_domain(
            "fs-roofline.dependency-verifier-placeholder.v1",
            b"placeholder",
        );
        let protocol = crate::CanonicalProductionOp {
            kernel_count: 1,
            fingerprint: 1,
            post_fingerprint: 1,
            run_nonce: placeholder,
            pre_axes_receipt: placeholder,
            post_axes_receipt: placeholder,
            dependency_receipt_digest: binding.domain_digest,
            dependency_receipt_artifact: binding.artifact_hash,
            finalized_run_receipt: placeholder,
            result_manifest: "{\"schema\":\"fs-roofline-run-manifest-v1\",\"entries\":[]}"
                .to_string(),
            baseline_admission: "{}".to_string(),
        };
        let mut missing_protocol = protocol.clone();
        missing_protocol.dependency_receipt_artifact =
            fs_blake3::hash_domain("fs-roofline.missing-dependency-receipt.v1", b"missing");
        assert!(
            !crate::dependency_receipt_is_structurally_valid(&ledger, op, &missing_protocol)
                .expect("missing-artifact verdict")
        );
        assert!(
            !crate::dependency_receipt_is_structurally_valid(&ledger, op, &protocol)
                .expect("missing-edge verdict"),
            "retained bytes without an op input edge are not lineage evidence"
        );
        ledger
            .link(op, &binding.artifact_hash, fs_ledger::EdgeRole::In)
            .expect("link receipt input");
        assert!(
            crate::dependency_receipt_is_structurally_valid(&ledger, op, &protocol)
                .expect("linked receipt verdict")
        );
        let substitute = DependencyReceiptBinding::from_parts(
            SUBSTITUTE_DEPGRAPH_RECEIPT,
            fs_blake3::hash_domain(
                fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN,
                SUBSTITUTE_DEPGRAPH_RECEIPT.as_bytes(),
            ),
        )
        .expect("substitute digest agrees with its own bytes");
        ledger
            .put_artifact(
                crate::DEPGRAPH_RECEIPT_ARTIFACT_KIND,
                substitute.bytes.as_bytes(),
                Some(crate::DEPGRAPH_RECEIPT_ARTIFACT_META),
            )
            .expect("store internally consistent substitute");
        ledger
            .link(op, &substitute.artifact_hash, fs_ledger::EdgeRole::In)
            .expect("link substitute input");
        let mut substituted_protocol = protocol.clone();
        substituted_protocol.dependency_receipt_digest = substitute.domain_digest;
        substituted_protocol.dependency_receipt_artifact = substitute.artifact_hash;
        assert!(
            crate::dependency_receipt_is_structurally_valid(&ledger, op, &substituted_protocol)
                .expect("substitution structure verdict"),
            "a historical receipt is validated against its own retained bytes"
        );
        assert!(
            !crate::dependency_receipt_matches_binding(&substituted_protocol, Some(binding)),
            "a structurally sound historical receipt must not impersonate today's build receipt"
        );
        ledger
            .finish_op(op, fs_ledger::OpOutcome::Ok, None, 2)
            .expect("finish verifier fixture after all lineage is attached");
        ledger
            .corrupt_artifact_for_test(&binding.artifact_hash)
            .expect("tamper receipt bytes");
        assert!(
            !crate::dependency_receipt_is_structurally_valid(&ledger, op, &protocol)
                .expect("tampered receipt verdict"),
            "content corruption must classify as invalid evidence, not a valid lineage edge"
        );
        cleanup_db(&db);
    }

    #[test]
    fn payload_artifact_envelope_requires_exact_kind_and_metadata() {
        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let wrong_kind = ledger
            .put_artifact(
                "untrusted-result",
                b"wrong-kind",
                Some(crate::ROOFLINE_PAYLOAD_ARTIFACT_META),
            )
            .expect("store wrong-kind fixture");
        assert!(
            !crate::artifact_envelope_is_valid(
                &ledger,
                &wrong_kind.hash,
                crate::ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                crate::ROOFLINE_PAYLOAD_ARTIFACT_META,
            )
            .expect("wrong-kind verdict")
        );

        let wrong_meta = ledger
            .put_artifact(
                crate::ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                b"wrong-meta",
                Some("{\"schema\":\"attacker-controlled\"}"),
            )
            .expect("store wrong-metadata fixture");
        assert!(
            !crate::artifact_envelope_is_valid(
                &ledger,
                &wrong_meta.hash,
                crate::ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                crate::ROOFLINE_PAYLOAD_ARTIFACT_META,
            )
            .expect("wrong-metadata verdict")
        );
    }

    #[test]
    fn post_probe_is_observed_strictly_after_every_timed_repetition() {
        let axes = synthetic_axes(0xB);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = attested_policy(&baseline, &identity);
        let runs = Rc::new(Cell::new(0_usize));
        let registry: Vec<Box<dyn crate::RooflineKernel>> = vec![Box::new(CountingKernel {
            runs: Rc::clone(&runs),
            value: 1,
        })];
        let probe = ProductionProbe::from_observed(axes.clone());
        let runs_at_post = Rc::new(Cell::new(usize::MAX));
        let observed = Rc::clone(&runs_at_post);
        let counter = Rc::clone(&runs);
        let config = ProductionRunConfig {
            n: 64,
            warmup: 2,
            reps: 3,
        };
        let run = probe
            .run_with_parts(config, policy, registry, move || {
                observed.set(counter.get());
                axes.clone()
            })
            .expect("protocol run");
        // warmup(2) + reps(3): the post-probe fired only after all five.
        assert_eq!(runs_at_post.get(), 5);
        assert_eq!(run.results().len(), 1);
    }

    #[test]
    fn drifted_post_probe_refuses_citation_and_records_a_rejection() {
        let axes = synthetic_axes(0xC);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = attested_policy(&baseline, &identity);
        let mut drifted = axes.clone();
        drifted.bandwidth_single_gbs *= 0.3;
        drifted.bandwidth_all_core_gbs *= 0.3;
        let probe = ProductionProbe::from_observed(axes);
        let run = probe
            .run_with_parts(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || drifted,
            )
            .expect("protocol run");
        assert!(
            !run.citation_eligible(),
            "drifted post-probe must refuse citation eligibility"
        );
        let reason = run.admission_error().expect("admission diagnostic");
        assert!(
            reason.contains("baseline admission refused"),
            "unexpected diagnostic: {reason}"
        );

        let db = temp_db("drift");
        let ledger = Ledger::open(&db).expect("open ledger");
        let kernel = run.results()[0].kernel.clone();
        let version = run.results()[0].version.clone();
        let fingerprint = run.axes().fingerprint;
        let baseline_hash = Some(baseline.content_hash());
        let recorded = run.record(&ledger).expect("record rejection");
        let op = recorded.op_id();
        let ir = ledger.op(op).unwrap().expect("op row").ir;
        assert!(ir.contains("\"protocol\":\"production-v3\""));
        assert!(ir.contains("\"admitted\":false"));
        let (store, retained) = live_baseline_store(&baseline);
        let authority = LiveTestAuthority::new(
            crate::KeyVerdict::Authorized,
            test_promotion_policy_receipt(&baseline),
        );
        let current = ProductionFreshnessContext::new(
            &store,
            &authority,
            &retained,
            &TEST_DEPENDENCY_AUTHORITY,
        );
        assert!(matches!(
            recorded.revalidate_at_for_test(
                &ledger,
                &current,
                recorded.recorded_at_ns() + 1,
                test_receipt_binding(),
            ),
            Err(ProductionFreshnessError::RecordedRefusal)
        ));
        // A rejected run publishes no tune evidence.
        assert_eq!(
            staleness_at(&ledger, &kernel, &version, fingerprint, baseline_hash, 1)
                .expect("staleness"),
            Staleness::NeverMeasured
        );
        cleanup_db(&db);
    }

    #[test]
    fn partial_finalizer_failure_yields_no_recordable_run() {
        let axes = synthetic_axes(0xD);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = attested_policy(&baseline, &identity);
        let registry: Vec<Box<dyn crate::RooflineKernel>> = vec![Box::new(FailingFinalizeKernel {
            inner: CountingKernel {
                runs: Rc::new(Cell::new(0)),
                value: 1,
            },
        })];
        let probe = ProductionProbe::from_observed(axes);
        let error = probe
            .run_with_parts(CONFIG, policy, registry, || synthetic_axes(0xD))
            .expect_err("finalizer failure must poison the whole run");
        assert!(
            error.contains("tune ledger unavailable mid-finalize"),
            "diagnostic must name the failing kernel's reason: {error}"
        );
        // No ProductionRun exists, so nothing can reach a ledger at all.
    }

    #[test]
    fn development_salt_is_report_only_even_when_measurements_admit() {
        let axes = synthetic_axes(0xD1);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = attested_policy(&baseline, &identity);
        let baseline_hash = Some(baseline.content_hash());
        let probe = ProductionProbe::from_observed(axes.clone());
        let post = axes.clone();
        let run = probe
            .run_with_parts_and_authority(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                CitationAuthority::Refused(DEVELOPMENT_SALT_REFUSAL),
            )
            .expect("report-only run");
        assert!(run.finalized.admitted(), "numeric admission is independent");
        assert!(
            !run.citation_eligible(),
            "a development salt is never citation evidence"
        );
        assert_eq!(
            run.admission_error().as_deref(),
            Some(DEVELOPMENT_SALT_REFUSAL)
        );
        let kernel = run.results()[0].kernel.clone();
        let version = run.results()[0].version.clone();

        let db = temp_db("salt-refusal");
        let ledger = Ledger::open(&db).expect("open ledger");
        let recorded = run.record(&ledger).expect("record structured refusal");
        let op = recorded.op_id();
        let row = ledger.op(op).unwrap().expect("refusal op");
        assert!(row.ir.contains("\"measurement_admitted\":true"));
        assert!(row.ir.contains("\"admitted\":false"));
        assert!(row.ir.contains("\"citation_refusal\":"));
        assert_eq!(row.outcome.as_deref(), Some("error"));
        assert_eq!(
            staleness_at(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                baseline_hash,
                row.t_end.expect("finished refusal"),
            )
            .expect("staleness"),
            Staleness::NeverMeasured,
        );
        cleanup_db(&db);
    }

    #[test]
    fn candidate_policy_records_only_structured_report_only_refusal() {
        const CONFIGURATION_REFUSAL: &str =
            "configured promotion authority key ops/perf-fixture is revoked";
        let axes = synthetic_axes(0xD2);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, test_day());
        let baseline_hash = Some(baseline.content_hash());
        let post = axes.clone();
        let run = ProductionProbe::from_observed(axes.clone())
            .run_report_only_with_parts(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
            )
            .expect("report-only run")
            .with_configuration_refusal(CONFIGURATION_REFUSAL.to_string())
            .expect("bounded configured refusal");
        assert!(run.finalized.admitted(), "the numerical measurement admits");
        assert!(run.receipt_json().contains("\"tier\":\"candidate\""));
        assert_eq!(
            run.admission_error().as_deref(),
            Some(CONFIGURATION_REFUSAL)
        );
        let kernel = run.results()[0].kernel.clone();

        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let op = run.record(&ledger).expect("record report-only refusal");
        let row = ledger.op(op).unwrap().expect("report-only op");
        assert!(row.ir.contains("\"protocol\":\"custom-registry\""));
        assert!(
            row.ir
                .contains(&format!("\"citation_refusal\":\"{CONFIGURATION_REFUSAL}\""))
        );
        assert!(row.ir.contains("\"measurement_admitted\":true"));
        assert!(row.ir.contains("\"admitted\":false"));
        assert_eq!(row.outcome.as_deref(), Some("error"));
        assert!(matches!(
            RecordedProductionRun::load(&ledger, op),
            Err(ProductionFreshnessError::CorruptRecordedEvidence)
        ));
        assert!(
            row.diag
                .as_deref()
                .is_some_and(|diag| diag.contains(CONFIGURATION_REFUSAL)),
            "the durable rejection diagnostic retains the exact configured refusal"
        );
        assert!(
            ledger
                .tune_rows(&kernel)
                .expect("candidate tune rows")
                .is_empty(),
            "report-only recording must not publish candidate or production tune rows"
        );
        assert!(
            ledger
                .op_artifact_hashes(op)
                .expect("report-only artifact edges")
                .is_empty(),
            "report-only recording must not retain result or dependency artifacts"
        );
        assert_eq!(
            staleness_at(
                &ledger,
                &kernel,
                "v1",
                axes.fingerprint,
                baseline_hash,
                row.t_end.expect("finished refusal"),
            )
            .expect("report-only staleness"),
            Staleness::NeverMeasured,
        );
    }

    #[test]
    fn oversized_configuration_refusal_is_bounded_and_remains_recordable() {
        let axes = synthetic_axes(0xD21);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, test_day());
        let post = axes.clone();
        let refusal = format!("configuration path {}", "é".repeat(3_000));
        let original_bytes = refusal.len();
        let digest = fs_blake3::hash_domain(REPORT_ONLY_REFUSAL_DIGEST_DOMAIN, refusal.as_bytes());
        let run = ProductionProbe::from_observed(axes)
            .run_report_only_with_parts(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
            )
            .expect("report-only run")
            .with_configuration_refusal(refusal)
            .expect("an oversized refusal is canonicalized, not rejected");
        assert!(run.report_only_refusal.len() <= MAX_REPORT_ONLY_REFUSAL_BYTES);
        assert!(
            run.report_only_refusal
                .contains(&format!("original_bytes={original_bytes}"))
        );
        assert!(
            run.report_only_refusal
                .contains(&format!("full_refusal_digest={digest}"))
        );

        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let op = run
            .record(&ledger)
            .expect("bounded report-only refusal remains recordable");
        let row = ledger.op(op).unwrap().expect("report-only operation");
        assert_eq!(row.outcome.as_deref(), Some("error"));
        assert!(row.ir.contains("\"citation_refusal\":"));
    }

    #[test]
    fn production_freezes_one_authority_decision_and_reuses_exact_snapshot() {
        struct CountingPromotionAuthority {
            calls: Cell<usize>,
            policy_receipt: fs_blake3::ContentHash,
        }

        impl crate::PromotionAuthorityVerifier for CountingPromotionAuthority {
            fn verify(
                &self,
                _key_id: &str,
                _signature: &str,
                _message: &[u8],
            ) -> crate::PromotionAuthorityDecision {
                self.calls.set(self.calls.get() + 1);
                crate::PromotionAuthorityDecision::new(
                    crate::KeyVerdict::Authorized,
                    self.policy_receipt,
                )
            }
        }

        let axes = synthetic_axes(0xD3);
        let (baseline, identity) = trusted_baseline(&axes);
        let retained: std::collections::BTreeSet<_> = baseline
            .provenance()
            .source_receipts()
            .iter()
            .copied()
            .collect();
        let authority = CountingPromotionAuthority {
            calls: Cell::new(0),
            policy_receipt: fs_blake3::hash_domain(
                "fs-roofline.counting-policy.v1",
                b"one atomic decision",
            ),
        };
        let attestation = crate::PromotionAttestation::new("counting-key", "signature");
        let mut store = crate::AttestedBaselineStore::new();
        store
            .admit_verified(baseline.clone(), attestation, &authority, &retained)
            .expect("admit fixture");
        authority.calls.set(0);
        let policy = store
            .policy_for_run(&identity, &authority, &retained)
            .expect("mint one-run policy");
        assert_eq!(authority.calls.get(), 1, "policy mint verifies once");

        let post = axes.clone();
        let run = ProductionProbe::from_observed(axes)
            .run_with_test_receipt(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                TEST_DEPGRAPH_RECEIPT,
            )
            .expect("production run");
        let exact_snapshot = run.receipt_json().to_string();
        assert_eq!(
            authority.calls.get(),
            1,
            "timing/finalization never reverify"
        );

        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let recorded = run.record(&ledger).expect("record exact snapshot");
        let op = recorded.op_id();
        assert_eq!(authority.calls.get(), 1, "recording never reverify");
        let ir = ledger.op(op).unwrap().expect("production op").ir;
        assert!(
            ir.ends_with(&format!("\"baseline_admission\":{exact_snapshot}}}")),
            "the exact frozen snapshot must be the operation's terminal admission receipt"
        );
    }

    #[test]
    fn delayed_finalized_run_records_only_a_structured_day_refusal() {
        let axes = synthetic_axes(0xD4);
        let (baseline, identity) = trusted_baseline(&axes);
        let today = test_day();
        let yesterday = today.checked_sub(1).expect("current epoch day is positive");
        let admission = attested_policy_on_day(&baseline, &identity, yesterday)
            .decide_at(&axes, &axes, yesterday);
        assert!(admission.authority_admitted());
        assert!(admission.verdict().trusted());
        assert!(!admission.baseline_citation_eligible());

        let mut registry = default_registry(1 << 10).expect("bounded registry fixture");
        let results = run_registry(&mut registry, 0, 1, &axes).expect("bounded delayed run");
        let finalized = FinalizedRegistryRun {
            receipt: crate::finalized_run_receipt(&admission, &results),
            admitted: true,
            consumed: false,
        };
        let kernel = results[0].kernel.clone();
        let run = ProductionRun {
            axes: axes.clone(),
            post_axes: axes,
            admission,
            nonce: fs_blake3::hash_domain("fs-roofline.delayed-run-test-nonce.v1", b"yesterday"),
            results,
            finalized,
            citation_authority: test_receipt_authority(),
            build_identity: crate::read_executable_build_identity()
                .expect("capture test executable identity"),
        };

        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let recorded = run
            .record(&ledger)
            .expect("record delayed structured refusal");
        let op = recorded.op_id();
        let row = ledger.op(op).unwrap().expect("delayed refusal op");
        assert!(row.ir.contains("\"measurement_admitted\":false"));
        assert!(row.ir.contains("\"admitted\":false"));
        assert_eq!(row.outcome.as_deref(), Some("error"));
        assert!(
            row.diag
                .as_deref()
                .is_some_and(|diag| diag.contains("current day")),
            "structured refusal must identify the expired decision day"
        );
        assert!(
            ledger
                .tune_rows(&kernel)
                .expect("delayed tune rows")
                .is_empty()
        );
        assert!(
            ledger
                .op_artifact_hashes(op)
                .expect("delayed artifact edges")
                .is_empty()
        );
    }

    #[test]
    fn frozen_finalization_refusal_controls_eligibility_and_diagnostic_together() {
        for _attempt in 0..2 {
            let today = test_day();
            let axes = synthetic_axes(0xD41);
            let (baseline, identity) = trusted_baseline(&axes);
            let admission =
                attested_policy_on_day(&baseline, &identity, today).decide_at(&axes, &axes, today);
            let mut registry = default_registry(1 << 10).expect("bounded registry fixture");
            let results = run_registry(&mut registry, 0, 1, &axes).expect("bounded frozen run");
            let finalized = FinalizedRegistryRun {
                receipt: crate::finalized_run_receipt(&admission, &results),
                admitted: false,
                consumed: false,
            };
            let run = ProductionRun {
                axes: axes.clone(),
                post_axes: axes,
                admission,
                nonce: fs_blake3::hash_domain(
                    "fs-roofline.finalization-refusal-test-nonce.v1",
                    b"not-admitted",
                ),
                results,
                finalized,
                citation_authority: test_receipt_authority(),
                build_identity: crate::read_executable_build_identity()
                    .expect("capture test executable identity"),
            };

            let refusal = run.admission_error();
            let eligible = run.citation_eligible();
            if test_day() != today {
                continue;
            }
            assert_eq!(refusal.as_deref(), Some(SEALED_FINALIZATION_REFUSAL));
            assert!(!eligible);
            return;
        }
        panic!("test clock crossed UTC midnight twice while checking finalization refusal");
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one end-to-end protocol and staleness state matrix
    fn successful_production_run_records_nonce_and_both_axis_receipts() {
        let axes = synthetic_axes(0xE);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = attested_policy(&baseline, &identity);
        let baseline_hash = Some(baseline.content_hash());
        let probe = ProductionProbe::from_observed(axes.clone());
        let nonce = probe.nonce;
        let post = axes.clone();
        let authority = test_receipt_authority();
        let CitationAuthority::Receipt(dependency_receipt) = authority else {
            unreachable!("test helper is receipt-backed")
        };
        let run = probe
            .run_with_parts_and_authority(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                authority,
            )
            .expect("protocol run");
        assert!(
            run.citation_eligible(),
            "stable synthetic probes must pass the synthetic eligibility fixture"
        );
        assert_eq!(run.nonce(), nonce);

        let db = temp_db("ok");
        let ledger = Ledger::open(&db).expect("open ledger");
        let kernel = run.results()[0].kernel.clone();
        let version = run.results()[0].version.clone();
        let recorded = run.record(&ledger).expect("record production run");
        let op = recorded.op_id();
        let row = ledger.op(op).unwrap().expect("op row");
        let recorded_at = row.t_end.expect("finished op");
        assert!(row.ir.contains("\"protocol\":\"production-v3\""));
        assert!(
            row.ir
                .contains("\"dependency_graph_evidence\":\"operator-observed-receipt\"")
        );
        assert!(row.ir.contains(&format!("\"run_nonce\":\"{nonce}\"")));
        assert!(
            row.ir
                .contains(&format!("\"pre_axes_receipt\":\"{}\"", axes_receipt(&axes)))
        );
        assert!(row.ir.contains(&format!(
            "\"post_axes_receipt\":\"{}\"",
            axes_receipt(&axes)
        )));
        assert!(row.ir.contains("\"admitted\":true"));
        assert!(
            ledger
                .edge_exists(
                    op,
                    &dependency_receipt.artifact_hash,
                    fs_ledger::EdgeRole::In
                )
                .expect("dependency receipt edge")
        );
        let dependency_info = ledger
            .artifact_info(&dependency_receipt.artifact_hash)
            .expect("dependency receipt metadata")
            .expect("retained dependency receipt");
        assert_eq!(dependency_info.kind, crate::DEPGRAPH_RECEIPT_ARTIFACT_KIND);
        assert_eq!(
            dependency_info.meta.as_deref(),
            Some(crate::DEPGRAPH_RECEIPT_ARTIFACT_META)
        );
        assert_eq!(
            ledger
                .get_artifact(&dependency_receipt.artifact_hash)
                .expect("dependency receipt bytes")
                .as_deref(),
            Some(TEST_DEPGRAPH_RECEIPT.as_bytes())
        );
        assert_eq!(
            receipt_staleness_at(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                baseline_hash,
                recorded_at + STALENESS_MAX_AGE_NS,
                dependency_receipt,
            )
            .expect("staleness"),
            Staleness::Fresh
        );
        assert_eq!(
            receipt_staleness_at(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                baseline_hash,
                recorded_at + STALENESS_MAX_AGE_NS + 1,
                dependency_receipt,
            )
            .expect("expired staleness"),
            Staleness::Expired
        );
        assert_eq!(
            receipt_staleness_at(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                baseline_hash,
                recorded_at - 1,
                dependency_receipt,
            )
            .expect("clock rollback staleness"),
            Staleness::ClockRollback
        );
        assert_eq!(
            staleness(&ledger, &kernel, &version, axes.fingerprint, None)
                .expect("missing baseline staleness"),
            Staleness::BaselineUnavailable
        );
        let foreign_baseline = fs_blake3::hash_domain("fs-roofline.foreign-baseline.v1", b"other");
        assert_eq!(
            staleness(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                Some(foreign_baseline),
            )
            .expect("baseline drift staleness"),
            Staleness::BaselineDrift
        );
        assert_eq!(
            staleness(
                &ledger,
                &kernel,
                "different-version",
                axes.fingerprint,
                baseline_hash,
            )
            .expect("version drift staleness"),
            Staleness::NeverMeasured
        );
        assert_eq!(
            staleness(&ledger, &kernel, &version, 0xFFFF, baseline_hash)
                .expect("fingerprint drift staleness"),
            Staleness::FingerprintDrift
        );

        let current_row = ledger
            .tune_rows(&kernel)
            .expect("production tune rows")
            .into_iter()
            .find(|row| {
                row.machine == roofline_machine_key(axes.fingerprint, baseline.content_hash())
            })
            .expect("current production row");
        ledger
            .tune_put(
                &current_row.kernel,
                &current_row.shape_class,
                &current_row.machine,
                &current_row.params,
                "{}",
            )
            .expect("inject valid-JSON payload corruption");
        assert_eq!(
            receipt_staleness_at(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                baseline_hash,
                recorded_at + 1,
                dependency_receipt,
            )
            .expect("corrupt production staleness"),
            Staleness::CorruptEvidence
        );
        cleanup_db(&db);
    }

    struct RecordedManifestRun {
        ledger: Ledger,
        baseline: BaselineAxes,
        recorded: RecordedProductionRun,
        kernels: Vec<(String, String)>,
        recorded_at: i64,
        dependency: DependencyReceiptBinding,
    }

    fn recorded_manifest_run(db: &str) -> RecordedManifestRun {
        let ledger = Ledger::open(db).expect("open ledger");
        let axes = synthetic_axes(0xBEEF);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = attested_policy(&baseline, &identity);
        let probe = ProductionProbe::from_observed(axes.clone());
        let post = axes.clone();
        let run = probe
            .run_with_parts_and_authority(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                test_receipt_authority(),
            )
            .expect("sealed manifest fixture");
        assert!(run.citation_eligible());
        let kernels = run
            .results()
            .iter()
            .map(|result| (result.kernel.clone(), result.version.clone()))
            .collect();
        let recorded = run.record(&ledger).expect("record sealed manifest fixture");
        let op = recorded.op_id();
        let recorded_at = ledger
            .op(op)
            .unwrap()
            .expect("recorded op")
            .t_end
            .expect("finished op");
        RecordedManifestRun {
            ledger,
            baseline,
            recorded,
            kernels,
            recorded_at,
            dependency: test_receipt_binding(),
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One table-like live-authority refusal matrix shares one sealed fixture.
    fn typed_record_requires_exact_live_authority_before_it_becomes_fresh() {
        let db = temp_db("typed-live-authority");
        let run = recorded_manifest_run(&db);
        let (store, retained) = live_baseline_store(&run.baseline);
        let policy_receipt = test_promotion_policy_receipt(&run.baseline);

        let loaded = RecordedProductionRun::load(&run.ledger, run.recorded.op_id())
            .expect("existing sealed operation reconstructs a typed receipt");
        assert_eq!(loaded.op_id(), run.recorded.op_id());
        assert_eq!(loaded.run_receipt(), run.recorded.run_receipt());
        assert_eq!(
            loaded.baseline_authority_fingerprint(),
            Some(policy_receipt)
        );
        assert_eq!(
            loaded.dependency_receipt_digest(),
            Some(run.dependency.domain_digest)
        );
        assert_eq!(loaded.recorded_at_ns(), run.recorded_at);

        let authority = LiveTestAuthority::new(crate::KeyVerdict::Authorized, policy_receipt);
        let dependency_policy_receipt =
            fs_blake3::hash_domain("fs-roofline.test-live-dependency-policy.v1", b"authorized");
        let dependency_authority = LiveTestDependencyAuthority::new(
            DependencyReceiptVerdict::Authorized,
            dependency_policy_receipt,
        );
        let current =
            ProductionFreshnessContext::new(&store, &authority, &retained, &dependency_authority);
        let fresh: FreshProductionEvidence = run
            .recorded
            .revalidate_at_for_test(&run.ledger, &current, run.recorded_at + 1, run.dependency)
            .expect("one exact live recheck mints positive evidence");
        assert_eq!(authority.calls.get(), 1, "one named live authority sample");
        assert_eq!(
            dependency_authority.calls.get(),
            1,
            "one named dependency-authority sample"
        );
        assert_eq!(fresh.op_id(), run.recorded.op_id());
        assert_eq!(fresh.run_receipt(), run.recorded.run_receipt());
        assert_eq!(fresh.baseline_authority_fingerprint(), policy_receipt);
        assert_eq!(
            fresh.dependency_authority_fingerprint(),
            dependency_policy_receipt
        );
        assert_ne!(
            fresh.dependency_authority_fingerprint(),
            run.dependency.domain_digest,
            "live authority identity must not be confused with evidence identity"
        );
        assert_eq!(fresh.recorded_at_ns(), run.recorded_at);
        assert_eq!(fresh.revalidated_at_ns(), run.recorded_at + 1);

        let revoked = LiveTestAuthority::new(crate::KeyVerdict::RevokedKey, policy_receipt);
        let current = ProductionFreshnessContext::new(
            &store,
            &revoked,
            &retained,
            &TEST_DEPENDENCY_AUTHORITY,
        );
        assert!(matches!(
            run.recorded.revalidate_at_for_test(
                &run.ledger,
                &current,
                run.recorded_at + 1,
                run.dependency,
            ),
            Err(ProductionFreshnessError::PromotionAuthorityRefused {
                verdict: crate::KeyVerdict::RevokedKey
            })
        ));
        assert_eq!(revoked.calls.get(), 1, "revocation is sampled once");

        let rotated_policy = fs_blake3::hash_domain(
            "fs-roofline.test-promotion-policy.v2",
            run.baseline.content_hash().as_bytes(),
        );
        let rotated = LiveTestAuthority::new(crate::KeyVerdict::Authorized, rotated_policy);
        let current = ProductionFreshnessContext::new(
            &store,
            &rotated,
            &retained,
            &TEST_DEPENDENCY_AUTHORITY,
        );
        assert!(matches!(
            run.recorded.revalidate_at_for_test(
                &run.ledger,
                &current,
                run.recorded_at + 1,
                run.dependency,
            ),
            Err(ProductionFreshnessError::PromotionPolicyChanged)
        ));
        assert_eq!(rotated.calls.get(), 1, "policy rotation is sampled once");

        let authority = LiveTestAuthority::new(crate::KeyVerdict::Authorized, policy_receipt);
        let current = ProductionFreshnessContext::new(
            &store,
            &authority,
            &retained,
            &TEST_DEPENDENCY_AUTHORITY,
        );
        assert!(matches!(
            run.recorded.revalidate_at_for_test(
                &run.ledger,
                &current,
                run.recorded_at - 1,
                run.dependency,
            ),
            Err(ProductionFreshnessError::ClockRollback { .. })
        ));
        assert!(matches!(
            run.recorded.revalidate_at_for_test(
                &run.ledger,
                &current,
                run.recorded_at + crate::STALENESS_MAX_AGE_NS + 1,
                run.dependency,
            ),
            Err(ProductionFreshnessError::Expired { .. })
        ));

        let substitute_digest = fs_blake3::hash_domain(
            fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN,
            SUBSTITUTE_DEPGRAPH_RECEIPT.as_bytes(),
        );
        let substitute =
            DependencyReceiptBinding::from_parts(SUBSTITUTE_DEPGRAPH_RECEIPT, substitute_digest)
                .expect("substitute dependency receipt is internally consistent");
        assert!(matches!(
            run.recorded.revalidate_at_for_test(
                &run.ledger,
                &current,
                run.recorded_at + 1,
                substitute,
            ),
            Err(ProductionFreshnessError::DependencyAuthorityChanged)
        ));

        let revoked_dependency = LiveTestDependencyAuthority::new(
            DependencyReceiptVerdict::Revoked,
            fs_blake3::hash_domain("fs-roofline.test-live-dependency-policy.v1", b"revoked"),
        );
        let current =
            ProductionFreshnessContext::new(&store, &authority, &retained, &revoked_dependency);
        assert!(matches!(
            run.recorded.revalidate_at_for_test(
                &run.ledger,
                &current,
                run.recorded_at + 1,
                run.dependency,
            ),
            Err(ProductionFreshnessError::DependencyAuthorityRevoked)
        ));
        assert_eq!(
            revoked_dependency.calls.get(),
            1,
            "dependency revocation is sampled exactly once"
        );

        let missing_source = *retained.iter().next().expect("fixture retains sources");
        let mut incomplete_retention = retained.clone();
        assert!(incomplete_retention.remove(&missing_source));
        let authority = LiveTestAuthority::new(crate::KeyVerdict::Authorized, policy_receipt);
        let current = ProductionFreshnessContext::new(
            &store,
            &authority,
            &incomplete_retention,
            &TEST_DEPENDENCY_AUTHORITY,
        );
        assert!(matches!(
            run.recorded.revalidate_at_for_test(
                &run.ledger,
                &current,
                run.recorded_at + 1,
                run.dependency,
            ),
            Err(ProductionFreshnessError::SourceReceiptUnavailable(receipt))
                if receipt == missing_source
        ));

        let mut replacement_axes = synthetic_axes(0xBEEF);
        replacement_axes.bandwidth_single_gbs *= 0.99;
        let (replacement, _) = trusted_baseline(&replacement_axes);
        assert_ne!(replacement.content_hash(), run.baseline.content_hash());
        let (replacement_store, replacement_retained) = live_baseline_store(&replacement);
        let replacement_authority = LiveTestAuthority::new(
            crate::KeyVerdict::Authorized,
            test_promotion_policy_receipt(&replacement),
        );
        let current = ProductionFreshnessContext::new(
            &replacement_store,
            &replacement_authority,
            &replacement_retained,
            &TEST_DEPENDENCY_AUTHORITY,
        );
        assert!(matches!(
            run.recorded.revalidate_at_for_test(
                &run.ledger,
                &current,
                run.recorded_at + 1,
                run.dependency,
            ),
            Err(ProductionFreshnessError::BaselineReplaced)
        ));
        cleanup_db(&db);
    }

    #[test]
    fn exact_typed_receipt_refuses_manifest_tamper_even_with_other_fresh_history() {
        let db = temp_db("typed-manifest-tamper");
        let run = recorded_manifest_run(&db);
        let other = recorded_manifest_run(&db);
        assert_ne!(run.recorded.op_id(), other.recorded.op_id());
        let (store, retained) = live_baseline_store(&run.baseline);
        let authority = LiveTestAuthority::new(
            crate::KeyVerdict::Authorized,
            test_promotion_policy_receipt(&run.baseline),
        );
        let current = ProductionFreshnessContext::new(
            &store,
            &authority,
            &retained,
            &TEST_DEPENDENCY_AUTHORITY,
        );
        let row = run
            .ledger
            .tune_rows(&run.kernels[0].0)
            .expect("read both exact-op rows")
            .into_iter()
            .find(|row| {
                crate::parse_roofline_row_params(&row.params)
                    .is_some_and(|params| params.op == run.recorded.op_id())
            })
            .expect("find the row owned by the receipt under attack");
        run.ledger
            .tune_put(
                &row.kernel,
                &row.shape_class,
                &row.machine,
                &row.params,
                "{}",
            )
            .expect("inject manifest-member payload tamper");
        assert!(matches!(
            run.recorded.revalidate_at_for_test(
                &run.ledger,
                &current,
                run.recorded_at + 1,
                run.dependency,
            ),
            Err(ProductionFreshnessError::CorruptRecordedEvidence)
        ));
        assert_eq!(
            authority.calls.get(),
            0,
            "corrupt exact evidence fails before consulting live authority"
        );
        assert_eq!(
            manifest_probe(&other, &other.kernels[0].0, &other.kernels[0].1),
            Staleness::Fresh,
            "the sibling honest operation remains fresh but cannot cover the damaged receipt"
        );
        cleanup_db(&db);
    }

    #[test]
    fn mutable_authority_is_sampled_once_per_admission_and_named_revalidation() {
        let db = temp_db("flipping-live-authority");
        let ledger = Ledger::open(&db).expect("open flipping-authority ledger");
        let axes = synthetic_axes(0xF11F);
        let (baseline, identity) = trusted_baseline(&axes);
        let (store, retained) = live_baseline_store(&baseline);
        let admitted_policy = test_promotion_policy_receipt(&baseline);
        let authority = FlippingTestAuthority {
            calls: Cell::new(0),
            admitted_policy,
            revoked_policy: fs_blake3::hash_domain(
                "fs-roofline.test-revoked-policy.v1",
                baseline.content_hash().as_bytes(),
            ),
        };
        let policy = store
            .policy_for_run_at(&identity, test_day(), &authority, &retained)
            .expect("first atomic authority sample admits the run");
        assert_eq!(authority.calls.get(), 1);
        let post = axes.clone();
        let run = ProductionProbe::from_observed(axes)
            .run_with_parts_and_authority(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                test_receipt_authority(),
            )
            .expect("sealed run reuses the immutable authority snapshot");
        assert_eq!(
            authority.calls.get(),
            1,
            "measurement/finalization never reverify"
        );
        let recorded = run
            .record(&ledger)
            .expect("record frozen authority receipt");
        assert_eq!(authority.calls.get(), 1, "recording never reverify");
        let current = ProductionFreshnessContext::new(
            &store,
            &authority,
            &retained,
            &TEST_DEPENDENCY_AUTHORITY,
        );
        assert!(matches!(
            recorded.revalidate_at_for_test(
                &ledger,
                &current,
                recorded.recorded_at_ns() + 1,
                test_receipt_binding(),
            ),
            Err(ProductionFreshnessError::PromotionAuthorityRefused {
                verdict: crate::KeyVerdict::RevokedKey
            })
        ));
        assert_eq!(
            authority.calls.get(),
            2,
            "the only second sample is the explicitly named live revalidation"
        );
        cleanup_db(&db);
    }

    #[test]
    fn key_rotation_recovers_only_after_new_attestation_and_new_recording() {
        let db = temp_db("authority-rotation-repromotion");
        let old = recorded_manifest_run(&db);
        let rotated_policy = fs_blake3::hash_domain(
            "fs-roofline.test-rotated-authority-policy.v1",
            old.baseline.content_hash().as_bytes(),
        );
        let rotated_authority =
            LiveTestAuthority::new(crate::KeyVerdict::Authorized, rotated_policy);
        let retained: BTreeSet<_> = old
            .baseline
            .provenance()
            .source_receipts()
            .iter()
            .copied()
            .collect();
        let mut rotated_store = crate::AttestedBaselineStore::new();
        rotated_store
            .admit_verified(
                old.baseline.clone(),
                crate::PromotionAttestation::new("rotated-authority", "rotated-signature"),
                &rotated_authority,
                &retained,
            )
            .expect("same immutable baseline is re-endorsed under the rotated authority");
        let rotated_current = ProductionFreshnessContext::new(
            &rotated_store,
            &rotated_authority,
            &retained,
            &TEST_DEPENDENCY_AUTHORITY,
        );
        assert!(matches!(
            old.recorded.revalidate_at_for_test(
                &old.ledger,
                &rotated_current,
                old.recorded_at + 1,
                old.dependency,
            ),
            Err(ProductionFreshnessError::PromotionAttestationChanged)
        ));

        let identity = old.baseline.identity().clone();
        let policy = rotated_store
            .policy_for_run_at(&identity, test_day(), &rotated_authority, &retained)
            .expect("rotated authority mints a new one-run policy");
        let axes = synthetic_axes(0xBEEF);
        let post = axes.clone();
        let run = ProductionProbe::from_observed(axes)
            .run_with_parts_and_authority(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                test_receipt_authority(),
            )
            .expect("new run freezes the rotated authority decision");
        let recorded = run
            .record(&old.ledger)
            .expect("new run records under the rotated authority");
        let fresh = recorded
            .revalidate_at_for_test(
                &old.ledger,
                &rotated_current,
                recorded.recorded_at_ns() + 1,
                test_receipt_binding(),
            )
            .expect("newly attested and newly recorded evidence becomes Fresh");
        assert_eq!(fresh.op_id(), recorded.op_id());
        assert_eq!(fresh.baseline_authority_fingerprint(), rotated_policy);
        cleanup_db(&db);
    }

    #[test]
    fn custom_registry_history_cannot_poison_a_fresh_production_row() {
        let db = temp_db("custom-history-isolation");
        let production = recorded_manifest_run(&db);
        let axes = synthetic_axes(0xBEEF);
        let identity = BaselineIdentity::current(&axes, "test-firmware")
            .expect("synthetic identity agrees with the retained baseline");
        let policy = AxisBaselinePolicy::new(Some(&production.baseline), &identity, test_day());
        let mut registry = default_registry(1 << 10).expect("bounded registry fixture");
        let mut results =
            run_registry(&mut registry, 0, 1, &axes).expect("bounded exploratory registry run");
        let mut finalized = finalize_registry_tuning(&mut registry, &axes, &axes, policy, &results)
            .expect("finalize exploratory run");
        crate::record_run(
            &production.ledger,
            &axes,
            &axes,
            policy,
            &mut finalized,
            &mut results,
        )
        .expect("record exploratory row in its candidate namespace");

        let (kernel, version) = &production.kernels[0];
        let rows = production
            .ledger
            .tune_rows(kernel)
            .expect("query production and candidate rows");
        assert_eq!(
            rows.iter()
                .filter(|row| row.shape_class.starts_with(crate::TUNE_SHAPE_CLASS))
                .count(),
            1
        );
        assert_eq!(
            rows.iter()
                .filter(|row| row.shape_class.starts_with(crate::CUSTOM_TUNE_SHAPE_CLASS))
                .count(),
            1
        );
        assert_eq!(
            manifest_probe(&production, kernel, version),
            Staleness::Fresh,
            "candidate history must not enter the production staleness scan"
        );
        cleanup_db(&db);
    }

    fn manifest_probe(run: &RecordedManifestRun, kernel: &str, version: &str) -> Staleness {
        receipt_staleness_at(
            &run.ledger,
            kernel,
            version,
            0xBEEF,
            Some(run.baseline.content_hash()),
            run.recorded_at + 1,
            run.dependency,
        )
        .expect("manifest staleness probe")
    }

    fn roofline_row(ledger: &Ledger, kernel: &str) -> fs_ledger::TuneRow {
        let mut rows: Vec<_> = ledger
            .tune_rows(kernel)
            .expect("tune rows")
            .into_iter()
            .filter(|row| row.shape_class.contains(":run="))
            .collect();
        assert_eq!(rows.len(), 1, "expected one roofline row for {kernel}");
        rows.pop().expect("row")
    }

    fn splice_payload(ledger: &Ledger, row: &fs_ledger::TuneRow, new_measured: &str) {
        let old_hash = fs_ledger::hash_bytes(row.measured.as_bytes()).to_string();
        let new_hash = fs_ledger::hash_bytes(new_measured.as_bytes());
        let artifact = ledger
            .put_artifact(
                crate::ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                new_measured.as_bytes(),
                Some("{\"schema\":\"fs-roofline-benchmark-result-v1\"}"),
            )
            .expect("store forged artifact");
        assert_eq!(artifact.hash, new_hash);
        let op: i64 = row
            .params
            .split_once("\"op\":")
            .and_then(|(_, rest)| rest.split_once(','))
            .and_then(|(digits, _)| digits.parse().ok())
            .expect("op id in params");
        assert!(matches!(
            ledger.link(op, &new_hash, fs_ledger::EdgeRole::Out),
            Err(fs_ledger::LedgerError::OpLineageSealed { op: sealed }) if sealed == op
        ));
        assert!(
            !ledger
                .edge_exists(op, &new_hash, fs_ledger::EdgeRole::Out)
                .expect("forged edge absence")
        );
        let forged_params = row.params.replace(&old_hash, &new_hash.to_string());
        assert_ne!(forged_params, row.params);
        ledger
            .tune_put(
                &row.kernel,
                &row.shape_class,
                &row.machine,
                &forged_params,
                new_measured,
            )
            .expect("overwrite row");
    }

    fn altered_measured(measured: &str) -> String {
        let (before, after) = measured
            .split_once("\"dispersion\":")
            .expect("dispersion field");
        let end = after.find([',', '}']).expect("field end");
        let forged = format!("{before}\"dispersion\":9.5e-1{}", &after[end..]);
        assert_ne!(forged, measured);
        forged
    }

    #[test]
    fn manifest_replacement_poisons_the_row_and_its_siblings() {
        let db = temp_db("manifest-splice");
        let run = recorded_manifest_run(&db);
        let (kernel_a, version_a) = run.kernels[0].clone();
        let (kernel_b, version_b) = run.kernels[1].clone();
        assert_eq!(
            manifest_probe(&run, &kernel_a, &version_a),
            Staleness::Fresh
        );
        assert_eq!(
            manifest_probe(&run, &kernel_b, &version_b),
            Staleness::Fresh
        );

        let row = roofline_row(&run.ledger, &kernel_a);
        splice_payload(&run.ledger, &row, &altered_measured(&row.measured));
        assert_eq!(
            manifest_probe(&run, &kernel_a, &version_a),
            Staleness::CorruptEvidence
        );
        assert_eq!(
            manifest_probe(&run, &kernel_b, &version_b),
            Staleness::CorruptEvidence
        );
        cleanup_db(&db);
    }

    #[test]
    fn sibling_parameter_tamper_poisons_every_manifest_member() {
        let db = temp_db("manifest-sibling-params");
        let run = recorded_manifest_run(&db);
        let (kernel_a, _) = run.kernels[0].clone();
        let (kernel_b, version_b) = run.kernels[1].clone();
        assert_eq!(
            manifest_probe(&run, &kernel_b, &version_b),
            Staleness::Fresh
        );

        let row = roofline_row(&run.ledger, &kernel_a);
        let tampered_params = row.params.replace("\"reps\":1,", "\"reps\":2,");
        assert_ne!(
            tampered_params, row.params,
            "fixture must alter sibling params"
        );
        run.ledger
            .tune_put(
                &row.kernel,
                &row.shape_class,
                &row.machine,
                &tampered_params,
                &row.measured,
            )
            .expect("overwrite sibling params");

        assert_eq!(
            manifest_probe(&run, &kernel_b, &version_b),
            Staleness::CorruptEvidence,
            "querying an untouched row must still validate every sibling's canonical params"
        );
        cleanup_db(&db);
    }

    #[test]
    fn sibling_artifact_corruption_poisons_every_manifest_member() {
        let db = temp_db("manifest-sibling-artifact");
        let run = recorded_manifest_run(&db);
        let (kernel_a, _) = run.kernels[0].clone();
        let (kernel_b, version_b) = run.kernels[1].clone();
        let row = roofline_row(&run.ledger, &kernel_a);
        let params = crate::parse_roofline_row_params(&row.params).expect("canonical row params");
        run.ledger
            .corrupt_artifact_for_test(&params.payload_artifact)
            .expect("corrupt sibling artifact bytes");

        assert_eq!(
            manifest_probe(&run, &kernel_b, &version_b),
            Staleness::CorruptEvidence,
            "an untouched row cannot stay Fresh when a sibling artifact is corrupt"
        );
        cleanup_db(&db);
    }

    #[test]
    fn production_operation_parser_rejects_noncanonical_and_ambiguous_ir() {
        let db = temp_db("canonical-operation-parser");
        let run = recorded_manifest_run(&db);
        let (kernel, _) = run.kernels[0].clone();
        let row = roofline_row(&run.ledger, &kernel);
        let params = crate::parse_roofline_row_params(&row.params).expect("canonical row params");
        let ir = run
            .ledger
            .op(params.op)
            .expect("query op")
            .expect("recorded op")
            .ir;
        let parsed = crate::parse_canonical_production_op(&ir).expect("canonical production IR");
        assert_eq!(parsed.to_json(), ir);
        let legacy_protocol = ir.replacen(
            "\"protocol\":\"production-v3\"",
            "\"protocol\":\"production-v2\"",
            1,
        );
        assert!(
            crate::parse_canonical_production_op(&legacy_protocol).is_none(),
            "operator-trusted production-v2 history must never re-enter the v3 citable parser"
        );
        let legacy_admission = ir.replacen(
            "fs-roofline-axis-admission-v2",
            "fs-roofline-axis-admission-v1",
            1,
        );
        assert!(
            crate::parse_canonical_production_op(&legacy_admission)
                .and_then(|legacy| crate::validate_protocol_axes(
                    &legacy,
                    0xBEEF,
                    run.baseline.content_hash(),
                ))
                .is_none(),
            "axis-admission-v1 history must fail the v3 attestation boundary"
        );
        assert!(
            crate::validate_protocol_axes(&parsed, 0xBEEF, run.baseline.content_hash(),).is_some(),
            "canonical pre/post receipts must bind the recorded fingerprints, axes, and baseline"
        );

        let source = run.baseline.provenance().source_receipts()[0];
        let foreign_source = fs_blake3::hash_domain(
            "fs-roofline.substituted-source-receipt.v1",
            source.as_bytes(),
        );
        let substituted_source_ir = ir.replacen(
            &format!("\"required_source_receipts\":[\"{source}\""),
            &format!("\"required_source_receipts\":[\"{foreign_source}\""),
            1,
        );
        let substituted_source = crate::parse_canonical_production_op(&substituted_source_ir)
            .expect("source substitution remains transport-canonical");
        assert!(
            crate::validate_protocol_axes(
                &substituted_source,
                0xBEEF,
                run.baseline.content_hash(),
            )
            .is_none(),
            "required source receipts must exactly equal canonical baseline provenance"
        );

        let substituted_identity_ir = ir.replacen(
            "\"firmware\":\"test-firmware\"",
            "\"firmware\":\"substituted-firmware\"",
            1,
        );
        let substituted_identity = crate::parse_canonical_production_op(&substituted_identity_ir)
            .expect("identity substitution remains transport-canonical");
        assert!(
            crate::validate_protocol_axes(
                &substituted_identity,
                0xBEEF,
                run.baseline.content_hash(),
            )
            .is_none(),
            "snapshot identity must exactly equal the canonical baseline identity"
        );

        let mut substituted_axes_receipt = parsed.clone();
        substituted_axes_receipt.pre_axes_receipt =
            fs_blake3::hash_domain("fs-roofline.substituted-axes-receipt.v1", b"substitute");
        assert!(
            crate::validate_protocol_axes(
                &substituted_axes_receipt,
                0xBEEF,
                run.baseline.content_hash(),
            )
            .is_none(),
            "a canonical operation cannot substitute a different pre-probe receipt"
        );

        let duplicate = ir.replacen(
            "\"measurement_admitted\":true,",
            "\"measurement_admitted\":true,\"measurement_admitted\":true,",
            1,
        );
        assert!(crate::parse_canonical_production_op(&duplicate).is_none());
        assert!(crate::parse_canonical_production_op(&format!("{ir} ")).is_none());
        let reordered = ir.replacen(
            "\"measurement_admitted\":true,\"admitted\":true",
            "\"admitted\":true,\"measurement_admitted\":true",
            1,
        );
        assert!(crate::parse_canonical_production_op(&reordered).is_none());
        cleanup_db(&db);
    }

    #[test]
    fn production_v2_tune_namespace_cannot_poison_v3_history() {
        let db = temp_db("legacy-production-namespace");
        let run = recorded_manifest_run(&db);
        let (kernel, version) = run.kernels[0].clone();
        let current = roofline_row(&run.ledger, &kernel);
        let legacy_shape = current
            .shape_class
            .replacen(crate::TUNE_SHAPE_CLASS, "roofline-v7", 1);
        assert_ne!(legacy_shape, current.shape_class);
        run.ledger
            .tune_put(
                &current.kernel,
                &legacy_shape,
                &current.machine,
                &current.params,
                &current.measured,
            )
            .expect("retain append-only production-v2 row");

        assert_eq!(
            manifest_probe(&run, &kernel, &version),
            Staleness::Fresh,
            "retained roofline-v7 history must be outside the production-v3 scan"
        );
        cleanup_db(&db);
    }

    #[test]
    fn manifest_rejects_rows_added_after_finalization() {
        let db = temp_db("manifest-added");
        let run = recorded_manifest_run(&db);
        let (kernel, version) = run.kernels[0].clone();
        let row = roofline_row(&run.ledger, &kernel);
        let ghost = "ghost-kernel";
        let ghost_measured = row.measured.replace(
            &format!("\"kernel\":\"{kernel}\""),
            &format!("\"kernel\":\"{ghost}\""),
        );
        let ghost_hash = fs_ledger::hash_bytes(ghost_measured.as_bytes());
        run.ledger
            .put_artifact(
                crate::ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                ghost_measured.as_bytes(),
                Some("{\"schema\":\"fs-roofline-benchmark-result-v1\"}"),
            )
            .expect("store ghost artifact");
        let op: i64 = row
            .params
            .split_once("\"op\":")
            .and_then(|(_, rest)| rest.split_once(','))
            .and_then(|(digits, _)| digits.parse().ok())
            .expect("op id in params");
        assert!(matches!(
            run.ledger
                .link(op, &ghost_hash, fs_ledger::EdgeRole::Out),
            Err(fs_ledger::LedgerError::OpLineageSealed { op: sealed }) if sealed == op
        ));
        assert!(
            !run.ledger
                .edge_exists(op, &ghost_hash, fs_ledger::EdgeRole::Out)
                .expect("ghost edge absence")
        );
        let ghost_params = row.params.replace(
            &fs_ledger::hash_bytes(row.measured.as_bytes()).to_string(),
            &ghost_hash.to_string(),
        );
        run.ledger
            .tune_put(
                ghost,
                &row.shape_class,
                &row.machine,
                &ghost_params,
                &ghost_measured,
            )
            .expect("insert ghost row");
        assert_eq!(
            manifest_probe(&run, ghost, &version),
            Staleness::CorruptEvidence
        );
        assert_eq!(manifest_probe(&run, &kernel, &version), Staleness::Fresh);
        cleanup_db(&db);
    }

    #[test]
    fn identical_receipt_backed_rerun_history_stays_fresh() {
        let db = temp_db("manifest-rerun");
        let first = recorded_manifest_run(&db);
        let (first_kernel, first_version) = first.kernels[0].clone();
        assert_eq!(
            manifest_probe(&first, &first_kernel, &first_version),
            Staleness::Fresh
        );

        let axes = synthetic_axes(0xBEEF);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = attested_policy(&baseline, &identity);
        let probe = ProductionProbe::from_observed(axes.clone());
        let post = axes.clone();
        let run = probe
            .run_with_parts_and_authority(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                test_receipt_authority(),
            )
            .expect("second sealed run");
        let second = run.record(&first.ledger).expect("record second run");
        let second_op = second.op_id();
        let rerecorded_at = first
            .ledger
            .op(second_op)
            .unwrap()
            .expect("second op")
            .t_end
            .expect("finished second op");
        for (kernel, version) in &first.kernels {
            assert_eq!(
                receipt_staleness_at(
                    &first.ledger,
                    kernel,
                    version,
                    0xBEEF,
                    Some(first.baseline.content_hash()),
                    rerecorded_at + 1,
                    first.dependency,
                )
                .expect("staleness probe"),
                Staleness::Fresh
            );
        }
        cleanup_db(&db);
    }

    #[test]
    fn pre_dependency_receipt_rows_are_retired_as_corrupt() {
        let db = temp_db("manifest-v3-row");
        let run = recorded_manifest_run(&db);
        let (kernel, version) = run.kernels[0].clone();
        assert_eq!(manifest_probe(&run, &kernel, &version), Staleness::Fresh);
        let row = roofline_row(&run.ledger, &kernel);
        let old_params = row.params.replace(
            "\"schema\":\"fs-roofline-ledger-row-v4\"",
            "\"schema\":\"fs-roofline-ledger-row-v3\"",
        );
        assert_ne!(old_params, row.params);
        run.ledger
            .tune_put(
                &row.kernel,
                &row.shape_class,
                &row.machine,
                &old_params,
                &row.measured,
            )
            .expect("downgrade row schema");
        assert_eq!(
            manifest_probe(&run, &kernel, &version),
            Staleness::CorruptEvidence
        );
        cleanup_db(&db);
    }

    fn cleanup_db(path: &str) {
        for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
            let _ = std::fs::remove_file(format!("{path}{suffix}"));
        }
    }
}
