//! Trusted historical axis baselines (bead dfh3): sustained-contention
//! detection that pre/post agreement cannot provide.
//!
//! The hole this closes: [`MachineAxes::reprobe_error`] only checks that
//! the pre-run and post-run probes AGREE — a host that was already
//! crushed before the first probe and stayed crushed through the second
//! (6 GB/s pre AND post on a normally 100+ GB/s machine) self-normalizes
//! and passes. Absolute floors stay as coarse last-resort sanity, but
//! they cannot be tight enough to catch a 10x degradation on a fast
//! machine without refusing slow-but-honest reference machines.
//!
//! The fix is a separately ADMITTED baseline: a fingerprint-specific
//! record of what this machine's axes measure when quiet, carrying
//! provenance (who promoted it, retained source-receipt identities, why),
//! an age policy,
//! and the environment identity (OS/arch/firmware declaration) it is
//! valid for. Citable gates then require the CURRENT axes to sit inside
//! declared bands around the trusted baseline.
//!
//! Trust discipline (the acceptance's four laws):
//! 1. First-run measurements are CANDIDATE evidence — nothing a probe
//!    measures about itself can authorize itself.
//! 2. Promotion is explicit and governed: at least
//!    [`MIN_PROMOTION_RUNS`] mutually-agreeing candidate runs plus a
//!    named operator annotation and a non-blank justification. The annotation
//!    is not an authenticated signature; protected-store/operator authority is
//!    an explicit trust boundary until a verifier capability is wired.
//! 3. Admission against the baseline refuses degraded, suspiciously
//!    fast, stale, and identity-drifted axes — each with a distinct,
//!    teaching verdict.
//! 4. Baseline updates go through the same promotion gate; there is no
//!    in-place mutation API.

use crate::authority::{KeyVerdict, PromotionAttestation, PromotionAuthorityVerifier};
use crate::axes::{MAX_AXIS_REPROBE_DRIFT, MachineAxes};
use fs_blake3::{ContentHash, hash_domain};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

/// Lower trust band: each current axis must be at least this fraction of
/// its baseline value. Below is SUSTAINED CONTENTION (the dfh3
/// counterexample measured 0.06 on the workbench).
pub const BASELINE_LOW_BAND: f64 = 0.70;

/// Upper trust band: each current axis must be at most this multiple of
/// its baseline value. Above means the machine is no longer the machine
/// the baseline describes (firmware/hardware change, or the baseline was
/// promoted from a degraded window) — re-promotion required either way.
pub const BASELINE_HIGH_BAND: f64 = 1.15;

/// Default and maximum baseline age policies, in days. A baseline older
/// than its policy is STALE: silent firmware/OS updates accumulate.
pub const DEFAULT_BASELINE_AGE_DAYS: u32 = 90;
/// Hard cap any policy must respect.
pub const MAX_BASELINE_AGE_DAYS: u32 = 365;

/// Minimum mutually-agreeing candidate runs behind one promotion.
/// This is a governance floor, not a statistical confidence claim: repeated
/// agreement cannot independently prove that the host was quiet.
pub const MIN_PROMOTION_RUNS: usize = 3;

/// Maximum encoded baseline-store size accepted by parsers and CLI readers.
pub const MAX_BASELINE_STORE_BYTES: usize = 1024 * 1024;
const MAX_BASELINE_LINE_BYTES: usize = 16 * 1024;
const MAX_BASELINE_STRING_BYTES: usize = 4096;

/// Canonical baseline-record schema. A policy or probe-semantics change must
/// use a new schema rather than silently reinterpreting old evidence.
pub const BASELINE_SCHEMA_VERSION: u32 = 1;

/// Domain for the content identity of one admitted baseline record.
pub const BASELINE_HASH_DOMAIN: &str = "frankensim.fs-roofline.baseline.v1";

/// The environment identity a baseline is valid for. `firmware` is a
/// DECLARED string (OS build / kernel release / SMC version — whatever
/// the operator's fleet discipline tracks): declared at promotion,
/// compared verbatim at admission. A mismatch is identity drift, never
/// a band question.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaselineIdentity {
    /// fs-substrate topology fingerprint.
    fingerprint: u64,
    /// CPU brand string.
    cpu_brand: String,
    /// Logical CPU count.
    logical_cpus: u32,
    /// Operating system (`std::env::consts::OS` at promotion).
    os: String,
    /// ISA (`std::env::consts::ARCH` at promotion).
    arch: String,
    /// Declared firmware/OS-build identity.
    firmware: String,
}

impl BaselineIdentity {
    /// The identity of the current process's environment for `axes`,
    /// with the operator's declared firmware string.
    ///
    /// # Errors
    /// Refuses blank, control-bearing, or oversized identity fields and an
    /// axis record with no logical CPUs.
    pub fn current(
        axes: &MachineAxes,
        firmware: impl Into<String>,
    ) -> Result<BaselineIdentity, PromotionError> {
        let identity = BaselineIdentity {
            fingerprint: axes.fingerprint,
            cpu_brand: axes.cpu_brand.clone(),
            logical_cpus: axes.logical_cpus,
            os: std::env::consts::OS.to_string(),
            arch: std::env::consts::ARCH.to_string(),
            firmware: firmware.into(),
        };
        validate_identity(&identity)?;
        Ok(identity)
    }

    /// Topology fingerprint.
    #[must_use]
    pub fn fingerprint(&self) -> u64 {
        self.fingerprint
    }

    /// CPU brand string.
    #[must_use]
    pub fn cpu_brand(&self) -> &str {
        &self.cpu_brand
    }

    /// Logical CPU count.
    #[must_use]
    pub fn logical_cpus(&self) -> u32 {
        self.logical_cpus
    }

    /// Operating-system identity.
    #[must_use]
    pub fn os(&self) -> &str {
        &self.os
    }

    /// Architecture identity.
    #[must_use]
    pub fn arch(&self) -> &str {
        &self.arch
    }

    /// Declared firmware/OS-build identity.
    #[must_use]
    pub fn firmware(&self) -> &str {
        &self.firmware
    }
}

/// One candidate axis probe plus the content identity of its retained source
/// receipt. The hash makes promotions traceable and duplicate-resistant; it
/// does not by itself authenticate who produced the receipt. Production
/// promotion therefore remains a trusted-operator boundary until a signature
/// verifier capability is wired.
#[derive(Debug, Clone)]
pub struct BaselineCandidate {
    axes: MachineAxes,
    identity: BaselineIdentity,
    source_receipt: ContentHash,
}

impl BaselineCandidate {
    /// Bind a plausible axis probe and its environment to a retained receipt.
    ///
    /// # Errors
    /// Refuses implausible axes or an identity that does not describe them.
    pub fn from_receipt(
        axes: MachineAxes,
        identity: BaselineIdentity,
        source_receipt: ContentHash,
    ) -> Result<Self, PromotionError> {
        validate_identity(&identity)?;
        if let Some(reason) = axes.plausibility_error() {
            return Err(PromotionError {
                detail: format!("candidate axes fail plausibility floors: {reason}"),
            });
        }
        validate_identity_matches_axes(&identity, &axes)?;
        Ok(Self {
            axes,
            identity,
            source_receipt,
        })
    }

    /// Candidate axes.
    #[must_use]
    pub fn axes(&self) -> &MachineAxes {
        &self.axes
    }

    /// Candidate environment identity.
    #[must_use]
    pub fn identity(&self) -> &BaselineIdentity {
        &self.identity
    }

    /// Retained source-receipt identity.
    #[must_use]
    pub fn source_receipt(&self) -> ContentHash {
        self.source_receipt
    }
}

/// Who promoted a baseline, from what, and why.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaselineProvenance {
    /// Named operator annotation (non-blank, but not authenticated).
    promoted_by: String,
    /// Non-blank justification recorded with the promotion.
    justification: String,
    /// Promotion day (days since the Unix epoch — see
    /// [`days_since_epoch_now`]).
    promoted_day: u64,
    /// Sorted, unique identities of the retained candidate receipts.
    source_receipts: Vec<ContentHash>,
}

impl BaselineProvenance {
    /// Named promoting operator. This is an annotation, not an authenticated
    /// signature.
    #[must_use]
    pub fn promoted_by(&self) -> &str {
        &self.promoted_by
    }

    /// Promotion justification.
    #[must_use]
    pub fn justification(&self) -> &str {
        &self.justification
    }

    /// Promotion day, as a Unix-epoch day.
    #[must_use]
    pub fn promoted_day(&self) -> u64 {
        self.promoted_day
    }

    /// Sorted, unique source-receipt identities.
    #[must_use]
    pub fn source_receipts(&self) -> &[ContentHash] {
        &self.source_receipts
    }

    /// Number of retained source receipts.
    #[must_use]
    pub fn source_runs(&self) -> usize {
        self.source_receipts.len()
    }
}

/// A trusted, admitted baseline for one machine fingerprint.
#[derive(Debug, Clone, PartialEq)]
pub struct BaselineAxes {
    /// Canonical record schema.
    schema_version: u32,
    /// The environment this baseline describes.
    identity: BaselineIdentity,
    /// Trusted single-thread STREAM bandwidth, GB/s.
    bandwidth_single_gbs: f64,
    /// Trusted all-core STREAM bandwidth, GB/s.
    bandwidth_all_core_gbs: f64,
    /// Trusted single-thread FMA throughput, GFLOP/s.
    peak_single_gflops: f64,
    /// Trusted all-core FMA throughput, GFLOP/s.
    peak_all_core_gflops: f64,
    /// Promotion provenance.
    provenance: BaselineProvenance,
    /// This baseline's age policy in days (≤ [`MAX_BASELINE_AGE_DAYS`]).
    age_policy_days: u32,
}

/// Why a promotion was refused.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PromotionError {
    /// The refusal, in teaching form.
    pub detail: String,
}

impl core::fmt::Display for PromotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "baseline promotion refused: {}", self.detail)
    }
}

impl core::error::Error for PromotionError {}

/// Wall-clock failure while establishing an epoch day for promotion or
/// admission. A citable path must fail closed rather than substitute day zero.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaselineClockError {
    detail: &'static str,
}

impl BaselineClockError {
    /// Stable diagnostic text.
    #[must_use]
    pub fn detail(&self) -> &'static str {
        self.detail
    }
}

impl core::fmt::Display for BaselineClockError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "baseline clock refused: {}", self.detail)
    }
}

impl core::error::Error for BaselineClockError {}

/// An attested baseline store (bead fz2.7): every record travels with
/// its [`PromotionAttestation`], and ADMISSION verifies authorization,
/// signature, and source-receipt availability against injected
/// capabilities — a locally editable store is no longer a silent trust
/// root. Transport lines are `{"record":<canonical>,"attestation":
/// {"key_id":..,"signature":..}}`; the record part is byte-identical to
/// the signed preimage, so replay is deterministic.
#[derive(Debug, Default)]
pub struct AttestedBaselineStore {
    store: BaselineStore,
    attestations: BTreeMap<u64, PromotionAttestation>,
}

impl AttestedBaselineStore {
    /// An empty store.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// The admitted baseline for a fingerprint, if any.
    #[must_use]
    pub fn for_fingerprint(&self, fingerprint: u64) -> Option<&BaselineAxes> {
        self.store.for_fingerprint(fingerprint)
    }

    /// The attestation admitted with a fingerprint's baseline.
    #[must_use]
    pub fn attestation_for(&self, fingerprint: u64) -> Option<&PromotionAttestation> {
        self.attestations.get(&fingerprint)
    }

    /// VERIFIED admission (the only way in): the attestation must be
    /// authorized by the injected authority over the record's exact
    /// content hash, and every source receipt the provenance names must
    /// be AVAILABLE (present in the retained-receipt set) — a promotion
    /// whose evidence has vanished is refused, not trusted.
    ///
    /// # Errors
    /// [`PromotionError`] naming the refusal: unattested records,
    /// forged/edited records (wrong signature), unknown or revoked
    /// keys, missing source receipts, and every structural refusal of
    /// the underlying store.
    pub fn admit_verified(
        &mut self,
        baseline: BaselineAxes,
        attestation: PromotionAttestation,
        authority: &dyn PromotionAuthorityVerifier,
        available_receipts: &BTreeSet<ContentHash>,
    ) -> Result<(), PromotionError> {
        if !attestation.well_formed() {
            return Err(PromotionError {
                detail: "promotion attestation has a blank key id or signature".to_string(),
            });
        }
        match baseline.authority_verdict(Some(&attestation), authority) {
            KeyVerdict::Authorized => {}
            refused => {
                return Err(PromotionError {
                    detail: format!(
                        "promotion authority refused key {:?}: {}",
                        attestation.key_id(),
                        refused.name()
                    ),
                });
            }
        }
        if let Some(missing) = baseline
            .provenance()
            .source_receipts()
            .iter()
            .find(|receipt| !available_receipts.contains(receipt))
        {
            return Err(PromotionError {
                detail: format!(
                    "source receipt {} named by the promotion is not available in the \
                     retained-receipt set — evidence must outlive the promotion it backs",
                    missing.to_hex()
                ),
            });
        }
        let fingerprint = baseline.identity().fingerprint();
        self.store.admit(baseline)?;
        self.attestations.insert(fingerprint, attestation);
        Ok(())
    }

    /// Serialize as attested JSON lines.
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        let mut out = String::new();
        for (fingerprint, attestation) in &self.attestations {
            let Some(baseline) = self.store.for_fingerprint(*fingerprint) else {
                continue;
            };
            out.push_str("{\"record\":");
            out.push_str(&baseline.canonical_json());
            out.push_str(",\"attestation\":{\"key_id\":");
            push_json_string(&mut out, attestation.key_id());
            out.push_str(",\"signature\":");
            push_json_string(&mut out, attestation.signature());
            out.push_str("}}\n");
        }
        out
    }

    /// Parse attested JSON lines STRICTLY. Parsing does NOT verify —
    /// no capability exists at parse time; every CITABLE USE re-verifies
    /// through [`citable_axis_admission_authorized`], so a tampered
    /// store line is caught at the decision point, deterministically.
    ///
    /// # Errors
    /// [`PromotionError`] naming the offending line.
    pub fn from_jsonl(text: &str) -> Result<Self, PromotionError> {
        if text.len() > MAX_BASELINE_STORE_BYTES {
            return Err(PromotionError {
                detail: format!("baseline store exceeds the {MAX_BASELINE_STORE_BYTES}-byte bound"),
            });
        }
        let mut attested = AttestedBaselineStore::new();
        for (line_number, line) in text.lines().enumerate() {
            let refuse = |why: &str| PromotionError {
                detail: format!("attested store line {}: {why}", line_number + 1),
            };
            let rest = line
                .strip_prefix("{\"record\":")
                .ok_or_else(|| refuse("missing the record envelope"))?;
            // The canonical record is a FLAT object: it ends at the first '}'.
            let record_end = rest
                .find('}')
                .ok_or_else(|| refuse("unterminated record object"))?;
            let record_json = &rest[..=record_end];
            let record_store = BaselineStore::from_jsonl(record_json)
                .map_err(|e| refuse(&format!("record part: {e}")))?;
            let tail = &rest[record_end + 1..];
            let tail = tail
                .strip_prefix(",\"attestation\":{\"key_id\":")
                .ok_or_else(|| refuse("missing the attestation envelope"))?;
            let (key_id, tail) =
                take_json_string(tail).ok_or_else(|| refuse("malformed key id"))?;
            let tail = tail
                .strip_prefix(",\"signature\":")
                .ok_or_else(|| refuse("missing the signature field"))?;
            let (signature, tail) =
                take_json_string(tail).ok_or_else(|| refuse("malformed signature"))?;
            if tail != "}}" {
                return Err(refuse("trailing bytes after the attestation"));
            }
            let attestation = PromotionAttestation::new(key_id, signature);
            if !attestation.well_formed() {
                return Err(refuse("blank key id or signature"));
            }
            for (fingerprint, baseline) in record_store.baselines {
                if attested.attestations.contains_key(&fingerprint) {
                    return Err(refuse(&format!("duplicate fingerprint {fingerprint:016x}")));
                }
                attested.store.baselines.insert(fingerprint, baseline);
                attested
                    .attestations
                    .insert(fingerprint, attestation.clone());
            }
        }
        Ok(attested)
    }
}

/// [`citable_axis_admission`] plus MANDATORY promotion-authority
/// verification (bead fz2.7): the citable tier for binding gates. The
/// baseline (when present) must carry an attestation the injected
/// authority accepts over its exact content hash; unattested, forged,
/// edited, unknown-key, and revoked-key records all refuse with a
/// typed [`BaselineVerdict::Unauthorized`] BEFORE any band math.
#[must_use]
#[allow(clippy::too_many_arguments)] // the complete admission context, spelled out
pub fn citable_axis_admission_authorized(
    pre: &MachineAxes,
    post: &MachineAxes,
    baseline: Option<&BaselineAxes>,
    attestation: Option<&PromotionAttestation>,
    identity: &BaselineIdentity,
    now_day: u64,
    authority: &dyn PromotionAuthorityVerifier,
) -> BaselineVerdict {
    if let Some(baseline) = baseline {
        match baseline.authority_verdict(attestation, authority) {
            KeyVerdict::Authorized => {}
            refused => {
                return BaselineVerdict::Unauthorized {
                    verdict: refused.name(),
                };
            }
        }
    }
    citable_axis_admission(pre, post, baseline, identity, now_day)
}

/// Take one JSON string literal off the front of `text`; returns the
/// unescaped value and the remaining tail.
fn take_json_string(text: &str) -> Option<(String, &str)> {
    let rest = text.strip_prefix('"')?;
    let mut out = String::new();
    let mut chars = rest.char_indices();
    loop {
        let (index, c) = chars.next()?;
        match c {
            '"' => return Some((out, &rest[index + 1..])),
            '\\' => {
                let (_, escaped) = chars.next()?;
                match escaped {
                    '"' => out.push('"'),
                    '\\' => out.push('\\'),
                    'n' => out.push('\n'),
                    'r' => out.push('\r'),
                    't' => out.push('\t'),
                    'u' => {
                        let mut code = 0u32;
                        for _ in 0..4 {
                            let (_, hex) = chars.next()?;
                            code = code * 16 + hex.to_digit(16)?;
                        }
                        out.push(char::from_u32(code)?);
                    }
                    _ => return None,
                }
            }
            c if c.is_control() => return None,
            c => out.push(c),
        }
    }
}

/// The verdict of checking current axes against a trusted baseline.
/// Only [`BaselineVerdict::Trusted`] supports citable gates.
#[derive(Debug, Clone, PartialEq)]
pub enum BaselineVerdict {
    /// Every axis sits inside the declared bands: the host is behaving
    /// like its trusted self.
    Trusted,
    /// No admitted baseline exists for this fingerprint: current
    /// measurements are CANDIDATE evidence only.
    Unbaselined,
    /// An axis fell below [`BASELINE_LOW_BAND`] of its baseline —
    /// sustained contention or thermal collapse; pre/post agreement is
    /// irrelevant.
    Degraded {
        /// Which axis.
        axis: &'static str,
        /// current / baseline.
        ratio: f64,
    },
    /// An axis exceeded [`BASELINE_HIGH_BAND`] of its baseline — the
    /// machine is not the machine the baseline describes; re-promote.
    Suspect {
        /// Which axis.
        axis: &'static str,
        /// current / baseline.
        ratio: f64,
    },
    /// The baseline is older than its age policy.
    Stale {
        /// Days since promotion.
        age_days: u64,
        /// The policy that was exceeded.
        limit_days: u32,
    },
    /// Fingerprint/topology/OS/arch/firmware mismatch: the baseline
    /// does not describe this environment at all.
    IdentityDrift {
        /// The first field that differed.
        field: &'static str,
    },
    /// The supplied axes were not credible enough to evaluate.
    InvalidAxes {
        /// Which probe failed.
        probe: &'static str,
        /// Stable plausibility diagnostic.
        reason: &'static str,
    },
    /// The pre/post probes did not corroborate one another.
    ReprobeFailed {
        /// Stable reprobe diagnostic.
        reason: &'static str,
    },
    /// The observed clock precedes the promotion timestamp. Saturating age
    /// arithmetic would make this baseline appear permanently young.
    ClockRollback {
        /// Observed Unix-epoch day.
        now_day: u64,
        /// Baseline promotion Unix-epoch day.
        promoted_day: u64,
    },
    /// The promotion-authority check refused this baseline (bead fz2.7):
    /// unattested, forged/edited (wrong signature), unknown key, or
    /// revoked key. Only an Authorized attestation supports citable use.
    Unauthorized {
        /// The typed authority verdict name (stable).
        verdict: &'static str,
    },
    /// An in-memory record violated the sealed baseline invariants.
    InvalidBaseline {
        /// Structural validation diagnostic.
        reason: String,
    },
}

impl BaselineVerdict {
    /// True only for [`BaselineVerdict::Trusted`].
    #[must_use]
    pub fn trusted(&self) -> bool {
        matches!(self, BaselineVerdict::Trusted)
    }

    /// One-line JSON for reports/ledger.
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        match self {
            BaselineVerdict::Trusted => "{\"baseline\":\"trusted\"}".to_string(),
            BaselineVerdict::Unbaselined => "{\"baseline\":\"unbaselined\"}".to_string(),
            BaselineVerdict::Degraded { axis, ratio } => {
                format!(
                    "{{\"baseline\":\"degraded\",\"axis\":\"{axis}\",\"ratio\":{}}}",
                    json_f64(*ratio)
                )
            }
            BaselineVerdict::Suspect { axis, ratio } => {
                format!(
                    "{{\"baseline\":\"suspect\",\"axis\":\"{axis}\",\"ratio\":{}}}",
                    json_f64(*ratio)
                )
            }
            BaselineVerdict::Stale {
                age_days,
                limit_days,
            } => format!(
                "{{\"baseline\":\"stale\",\"age_days\":{age_days},\"limit_days\":{limit_days}}}"
            ),
            BaselineVerdict::IdentityDrift { field } => {
                format!("{{\"baseline\":\"identity-drift\",\"field\":\"{field}\"}}")
            }
            BaselineVerdict::Unauthorized { verdict } => {
                format!("{{\"baseline\":\"unauthorized\",\"authority\":\"{verdict}\"}}")
            }
            BaselineVerdict::InvalidAxes { probe, reason } => format!(
                "{{\"baseline\":\"invalid-axes\",\"probe\":\"{probe}\",\"reason\":\"{}\"}}",
                json_escaped(reason)
            ),
            BaselineVerdict::ReprobeFailed { reason } => format!(
                "{{\"baseline\":\"reprobe-failed\",\"reason\":\"{}\"}}",
                json_escaped(reason)
            ),
            BaselineVerdict::ClockRollback {
                now_day,
                promoted_day,
            } => format!(
                "{{\"baseline\":\"clock-rollback\",\"now_day\":{now_day},\"promoted_day\":{promoted_day}}}"
            ),
            BaselineVerdict::InvalidBaseline { reason } => format!(
                "{{\"baseline\":\"invalid-record\",\"reason\":\"{}\"}}",
                json_escaped(reason)
            ),
        }
    }
}

/// Days since the Unix epoch, from the system clock. Tests inject their
/// own day; production callers use this.
///
/// # Errors
/// Refuses a wall clock before the Unix epoch. Returning zero here would make
/// future-dated records appear fresh through saturating subtraction.
pub fn days_since_epoch_now() -> Result<u64, BaselineClockError> {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() / 86_400)
        .map_err(|_| BaselineClockError {
            detail: "system wall clock precedes the Unix epoch",
        })
}

fn json_f64(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.3}")
    } else {
        "null".to_string()
    }
}

fn json_escaped(value: &str) -> String {
    let mut out = String::new();
    push_json_string_body(&mut out, value);
    out
}

fn axes_quad(axes: &MachineAxes) -> [(&'static str, f64); 4] {
    [
        ("bandwidth_single_gbs", axes.bandwidth_single_gbs),
        ("bandwidth_all_core_gbs", axes.bandwidth_all_core_gbs),
        ("peak_single_gflops", axes.peak_single_gflops),
        ("peak_all_core_gflops", axes.peak_all_core_gflops),
    ]
}

fn validate_text(field: &str, value: &str) -> Result<(), PromotionError> {
    if value.trim().is_empty() {
        return Err(PromotionError {
            detail: format!("{field} must be non-blank"),
        });
    }
    if value.len() > MAX_BASELINE_STRING_BYTES {
        return Err(PromotionError {
            detail: format!("{field} exceeds the {MAX_BASELINE_STRING_BYTES}-byte string bound"),
        });
    }
    if value.chars().any(char::is_control) {
        return Err(PromotionError {
            detail: format!("{field} must not contain control characters"),
        });
    }
    Ok(())
}

fn validate_identity(identity: &BaselineIdentity) -> Result<(), PromotionError> {
    if identity.logical_cpus == 0 {
        return Err(PromotionError {
            detail: "identity logical CPU count is zero".to_string(),
        });
    }
    for (field, value) in [
        ("cpu_brand", identity.cpu_brand.as_str()),
        ("os", identity.os.as_str()),
        ("arch", identity.arch.as_str()),
        ("firmware", identity.firmware.as_str()),
    ] {
        validate_text(field, value)?;
    }
    Ok(())
}

fn validate_identity_matches_axes(
    identity: &BaselineIdentity,
    axes: &MachineAxes,
) -> Result<(), PromotionError> {
    if identity.fingerprint != axes.fingerprint {
        return Err(PromotionError {
            detail: "candidate identity has a different machine fingerprint".to_string(),
        });
    }
    if identity.logical_cpus != axes.logical_cpus || identity.cpu_brand != axes.cpu_brand {
        return Err(PromotionError {
            detail: "candidate identity has a different topology identity".to_string(),
        });
    }
    Ok(())
}

fn validate_baseline(baseline: &BaselineAxes) -> Result<(), PromotionError> {
    if baseline.schema_version != BASELINE_SCHEMA_VERSION {
        return Err(PromotionError {
            detail: format!(
                "unsupported baseline schema {} (expected {BASELINE_SCHEMA_VERSION})",
                baseline.schema_version
            ),
        });
    }
    validate_identity(&baseline.identity)?;
    validate_text("promoted_by", &baseline.provenance.promoted_by)?;
    validate_text("justification", &baseline.provenance.justification)?;
    if baseline.provenance.source_receipts.len() < MIN_PROMOTION_RUNS {
        return Err(PromotionError {
            detail: format!(
                "baseline requires at least {MIN_PROMOTION_RUNS} source receipts, got {}",
                baseline.provenance.source_receipts.len()
            ),
        });
    }
    if baseline
        .provenance
        .source_receipts
        .windows(2)
        .any(|pair| matches!(pair, [left, right] if left >= right))
    {
        return Err(PromotionError {
            detail: "source receipt identities must be sorted and unique".to_string(),
        });
    }
    if baseline.age_policy_days == 0 || baseline.age_policy_days > MAX_BASELINE_AGE_DAYS {
        return Err(PromotionError {
            detail: format!(
                "age policy {} days is outside 1..={MAX_BASELINE_AGE_DAYS}",
                baseline.age_policy_days
            ),
        });
    }
    let axes = MachineAxes {
        fingerprint: baseline.identity.fingerprint,
        cpu_brand: baseline.identity.cpu_brand.clone(),
        logical_cpus: baseline.identity.logical_cpus,
        bandwidth_single_gbs: baseline.bandwidth_single_gbs,
        bandwidth_all_core_gbs: baseline.bandwidth_all_core_gbs,
        peak_single_gflops: baseline.peak_single_gflops,
        peak_all_core_gflops: baseline.peak_all_core_gflops,
    };
    if let Some(reason) = axes.plausibility_error() {
        return Err(PromotionError {
            detail: format!("baseline axes fail plausibility: {reason}"),
        });
    }
    if baseline.canonical_json_unchecked().len() > MAX_BASELINE_LINE_BYTES {
        return Err(PromotionError {
            detail: format!(
                "canonical baseline exceeds the {MAX_BASELINE_LINE_BYTES}-byte line bound"
            ),
        });
    }
    Ok(())
}

impl BaselineAxes {
    fn baseline_quad(&self) -> [(&'static str, f64); 4] {
        [
            ("bandwidth_single_gbs", self.bandwidth_single_gbs),
            ("bandwidth_all_core_gbs", self.bandwidth_all_core_gbs),
            ("peak_single_gflops", self.peak_single_gflops),
            ("peak_all_core_gflops", self.peak_all_core_gflops),
        ]
    }

    /// Canonical record schema.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Environment identity described by this baseline.
    #[must_use]
    pub fn identity(&self) -> &BaselineIdentity {
        &self.identity
    }

    /// Trusted single-thread memory bandwidth, GB/s.
    #[must_use]
    pub fn bandwidth_single_gbs(&self) -> f64 {
        self.bandwidth_single_gbs
    }

    /// Trusted all-core memory bandwidth, GB/s.
    #[must_use]
    pub fn bandwidth_all_core_gbs(&self) -> f64 {
        self.bandwidth_all_core_gbs
    }

    /// Trusted single-thread compute peak, GFLOP/s.
    #[must_use]
    pub fn peak_single_gflops(&self) -> f64 {
        self.peak_single_gflops
    }

    /// Trusted all-core compute peak, GFLOP/s.
    #[must_use]
    pub fn peak_all_core_gflops(&self) -> f64 {
        self.peak_all_core_gflops
    }

    /// Promotion provenance, including every source-receipt identity.
    #[must_use]
    pub fn provenance(&self) -> &BaselineProvenance {
        &self.provenance
    }

    /// Maximum trusted age in days.
    #[must_use]
    pub fn age_policy_days(&self) -> u32 {
        self.age_policy_days
    }

    /// Check `current` axes (already floor-plausible) against this
    /// baseline on day `now_day`.
    #[must_use]
    pub fn verdict(
        &self,
        current: &MachineAxes,
        identity: &BaselineIdentity,
        now_day: u64,
    ) -> BaselineVerdict {
        if let Err(error) = validate_baseline(self) {
            return BaselineVerdict::InvalidBaseline {
                reason: error.detail,
            };
        }
        if let Some(reason) = current.plausibility_error() {
            return BaselineVerdict::InvalidAxes {
                probe: "current",
                reason,
            };
        }
        if identity.fingerprint != current.fingerprint
            || self.identity.fingerprint != identity.fingerprint
        {
            return BaselineVerdict::IdentityDrift {
                field: "fingerprint",
            };
        }
        if identity.logical_cpus != current.logical_cpus
            || self.identity.logical_cpus != identity.logical_cpus
        {
            return BaselineVerdict::IdentityDrift {
                field: "logical_cpus",
            };
        }
        if identity.cpu_brand != current.cpu_brand || self.identity.cpu_brand != identity.cpu_brand
        {
            return BaselineVerdict::IdentityDrift { field: "cpu_brand" };
        }
        if identity.os != std::env::consts::OS || self.identity.os != identity.os {
            return BaselineVerdict::IdentityDrift { field: "os" };
        }
        if identity.arch != std::env::consts::ARCH || self.identity.arch != identity.arch {
            return BaselineVerdict::IdentityDrift { field: "arch" };
        }
        if self.identity.firmware != identity.firmware {
            return BaselineVerdict::IdentityDrift { field: "firmware" };
        }
        if now_day < self.provenance.promoted_day {
            return BaselineVerdict::ClockRollback {
                now_day,
                promoted_day: self.provenance.promoted_day,
            };
        }
        let age_days = now_day - self.provenance.promoted_day;
        if age_days > u64::from(self.age_policy_days) {
            return BaselineVerdict::Stale {
                age_days,
                limit_days: self.age_policy_days,
            };
        }
        for ((axis, current_value), (_, trusted_value)) in
            axes_quad(current).into_iter().zip(self.baseline_quad())
        {
            let ratio = current_value / trusted_value;
            if !ratio.is_finite() || ratio < BASELINE_LOW_BAND {
                return BaselineVerdict::Degraded { axis, ratio };
            }
            if ratio > BASELINE_HIGH_BAND {
                return BaselineVerdict::Suspect { axis, ratio };
            }
        }
        BaselineVerdict::Trusted
    }

    fn canonical_json_unchecked(&self) -> String {
        let mut s = String::with_capacity(768);
        let _ = write!(
            s,
            "{{\"schema_version\":{},\"low_band_bits\":\"{:016x}\",\"high_band_bits\":\"{:016x}\",\"promotion_drift_bits\":\"{:016x}\",\"fingerprint\":\"{:016x}\",\"cpu_brand\":",
            self.schema_version,
            BASELINE_LOW_BAND.to_bits(),
            BASELINE_HIGH_BAND.to_bits(),
            MAX_AXIS_REPROBE_DRIFT.to_bits(),
            self.identity.fingerprint
        );
        push_json_string(&mut s, &self.identity.cpu_brand);
        let _ = write!(
            s,
            ",\"logical_cpus\":{},\"os\":",
            self.identity.logical_cpus
        );
        push_json_string(&mut s, &self.identity.os);
        s.push_str(",\"arch\":");
        push_json_string(&mut s, &self.identity.arch);
        s.push_str(",\"firmware\":");
        push_json_string(&mut s, &self.identity.firmware);
        let _ = write!(
            s,
            ",\"bandwidth_single_bits\":\"{:016x}\",\"bandwidth_all_core_bits\":\"{:016x}\",\
             \"peak_single_bits\":\"{:016x}\",\"peak_all_core_bits\":\"{:016x}\",\"source_receipts\":[",
            self.bandwidth_single_gbs.to_bits(),
            self.bandwidth_all_core_gbs.to_bits(),
            self.peak_single_gflops.to_bits(),
            self.peak_all_core_gflops.to_bits(),
        );
        for (index, receipt) in self.provenance.source_receipts.iter().enumerate() {
            if index > 0 {
                s.push(',');
            }
            push_json_string(&mut s, &receipt.to_hex());
        }
        s.push_str("],\"promoted_by\":");
        push_json_string(&mut s, &self.provenance.promoted_by);
        s.push_str(",\"justification\":");
        push_json_string(&mut s, &self.provenance.justification);
        let _ = write!(
            s,
            ",\"promoted_day\":{},\"source_runs\":{},\"age_policy_days\":{}}}",
            self.provenance.promoted_day,
            self.provenance.source_receipts.len(),
            self.age_policy_days
        );
        s
    }

    /// Canonical, self-contained JSON preimage. All floats are represented by
    /// their exact bits and source receipts are sorted.
    #[must_use]
    pub fn canonical_json(&self) -> String {
        self.canonical_json_unchecked()
    }

    /// Domain-separated identity of [`Self::canonical_json`].
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        hash_domain(BASELINE_HASH_DOMAIN, self.canonical_json().as_bytes())
    }

    /// The exact bytes a promotion authority signs (bead fz2.7): the
    /// record's content hash, which already binds the canonical
    /// schema/domain hash, sorted source receipt identities, band and
    /// drift policy, machine identity, and promotion time — signing the
    /// hash signs all of them, and editing ANY of them (including the
    /// free-text operator) invalidates the signature.
    #[must_use]
    pub fn promotion_message(&self) -> [u8; 32] {
        *self.content_hash().as_bytes()
    }

    /// Judge this record's attestation against an injected authority.
    /// `None` (unattested) is UnknownKey — never authorized.
    #[must_use]
    pub fn authority_verdict(
        &self,
        attestation: Option<&PromotionAttestation>,
        authority: &dyn PromotionAuthorityVerifier,
    ) -> KeyVerdict {
        match attestation {
            None => KeyVerdict::UnknownKey,
            Some(attestation) if !attestation.well_formed() => KeyVerdict::UnknownKey,
            Some(attestation) => authority.verify(
                attestation.key_id(),
                attestation.signature(),
                &self.promotion_message(),
            ),
        }
    }

    /// Backwards-compatible JSON-lines spelling (without the trailing newline).
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        self.canonical_json()
    }
}

// Compute the conservative per-axis maximum after enforcing mutual agreement.
fn mutually_agreeing_maxima(
    candidates: &[BaselineCandidate],
    first: &BaselineCandidate,
) -> Result<[f64; 4], PromotionError> {
    let mut minima = [f64::INFINITY; 4];
    let mut maxima = [0.0f64; 4];
    for candidate in candidates {
        for ((minimum, maximum), (_, value)) in minima
            .iter_mut()
            .zip(maxima.iter_mut())
            .zip(axes_quad(&candidate.axes))
        {
            *minimum = minimum.min(value);
            *maximum = maximum.max(value);
        }
    }
    let mut promoted = [0.0f64; 4];
    for (((axis, minimum), maximum), promoted_axis) in axes_quad(&first.axes)
        .map(|(axis, _)| axis)
        .into_iter()
        .zip(minima)
        .zip(maxima)
        .zip(promoted.iter_mut())
    {
        if (maximum - minimum) / maximum > MAX_AXIS_REPROBE_DRIFT {
            return Err(PromotionError {
                detail: format!(
                    "candidate runs disagree on {axis} beyond the {MAX_AXIS_REPROBE_DRIFT} drift \
                     band ({minimum:.2} .. {maximum:.2}) — measure on a quiet host"
                ),
            });
        }
        *promoted_axis = maximum;
    }
    Ok(promoted)
}

/// Promote a trusted baseline from candidate runs — THE only way a
/// baseline comes to exist.
///
/// Requirements (each refused with a teaching detail):
/// - at least [`MIN_PROMOTION_RUNS`] candidates with unique retained receipt
///   identities;
/// - every run floor-plausible, same complete environment identity, and every
///   axis pair within [`MAX_AXIS_REPROBE_DRIFT`] of the run minimum
///   (mutual agreement — one quiet run among crushed ones cannot
///   launder the set);
/// - a named operator annotation and non-blank justification;
/// - an age policy within [`MAX_BASELINE_AGE_DAYS`].
///
/// The promoted axes are the per-axis MAXIMUM over the runs: the best
/// mutually-corroborated measurement is the closest to the machine's
/// true quiet capability, and a too-low baseline would inflate every
/// later attainment claim.
///
/// This function makes the record structurally traceable; it cannot prove that
/// a caller-supplied source hash came from an authentic probe. Signature-backed
/// source receipt verification remains a no-claim boundary.
///
/// # Errors
/// [`PromotionError`] naming the failed requirement.
pub fn promote_baseline(
    candidates: &[BaselineCandidate],
    promoted_by: impl Into<String>,
    justification: impl Into<String>,
    promoted_day: u64,
    age_policy_days: u32,
) -> Result<BaselineAxes, PromotionError> {
    let promoted_by = promoted_by.into();
    let justification = justification.into();
    if promoted_by.trim().is_empty() {
        return Err(PromotionError {
            detail: "promotion requires a named operator".to_string(),
        });
    }
    validate_text("promoted_by", &promoted_by)?;
    if justification.trim().is_empty() {
        return Err(PromotionError {
            detail: "promotion requires a non-blank justification".to_string(),
        });
    }
    validate_text("justification", &justification)?;
    if age_policy_days == 0 || age_policy_days > MAX_BASELINE_AGE_DAYS {
        return Err(PromotionError {
            detail: format!(
                "age policy {age_policy_days} days is outside 1..={MAX_BASELINE_AGE_DAYS}"
            ),
        });
    }
    if candidates.len() < MIN_PROMOTION_RUNS {
        return Err(PromotionError {
            detail: format!(
                "promotion requires at least {MIN_PROMOTION_RUNS} candidate runs, got {}",
                candidates.len()
            ),
        });
    }
    let first = candidates.first().ok_or_else(|| PromotionError {
        detail: "promotion has no first candidate".to_string(),
    })?;
    let identity = first.identity.clone();
    validate_identity(&identity)?;
    let mut source_receipts = BTreeSet::new();
    for (index, candidate) in candidates.iter().enumerate() {
        if let Some(reason) = candidate.axes.plausibility_error() {
            return Err(PromotionError {
                detail: format!("candidate run {index} fails plausibility floors: {reason}"),
            });
        }
        validate_identity_matches_axes(&candidate.identity, &candidate.axes)?;
        if candidate.identity != identity {
            return Err(PromotionError {
                detail: format!("candidate run {index} has a different environment identity"),
            });
        }
        if !source_receipts.insert(candidate.source_receipt) {
            return Err(PromotionError {
                detail: format!("candidate run {index} reuses a source receipt identity"),
            });
        }
    }
    // Mutual agreement: for each axis, max/min across runs must sit within the
    // reprobe drift band. The maximum is conservative against inflated claims.
    let [
        bandwidth_single_gbs,
        bandwidth_all_core_gbs,
        peak_single_gflops,
        peak_all_core_gflops,
    ] = mutually_agreeing_maxima(candidates, first)?;
    let baseline = BaselineAxes {
        schema_version: BASELINE_SCHEMA_VERSION,
        identity,
        bandwidth_single_gbs,
        bandwidth_all_core_gbs,
        peak_single_gflops,
        peak_all_core_gflops,
        provenance: BaselineProvenance {
            promoted_by,
            justification,
            promoted_day,
            source_receipts: source_receipts.into_iter().collect(),
        },
        age_policy_days,
    };
    validate_baseline(&baseline)?;
    Ok(baseline)
}

/// The combined citable-axis admission: absolute floors (last-resort
/// sanity), pre/post agreement, AND baseline trust. `baseline = None`
/// yields [`BaselineVerdict::Unbaselined`] — measurements proceed as
/// candidate evidence but nothing citable may be minted from them.
#[must_use]
pub fn citable_axis_admission(
    pre: &MachineAxes,
    post: &MachineAxes,
    baseline: Option<&BaselineAxes>,
    identity: &BaselineIdentity,
    now_day: u64,
) -> BaselineVerdict {
    if let Some(reason) = pre.plausibility_error() {
        return BaselineVerdict::InvalidAxes {
            probe: "pre",
            reason,
        };
    }
    if let Some(reason) = post.plausibility_error() {
        return BaselineVerdict::InvalidAxes {
            probe: "post",
            reason,
        };
    }
    if let Some(reason) = pre.reprobe_error(post) {
        return BaselineVerdict::ReprobeFailed { reason };
    }
    match baseline {
        None => BaselineVerdict::Unbaselined,
        Some(trusted) => {
            let pre_verdict = trusted.verdict(pre, identity, now_day);
            if !pre_verdict.trusted() {
                return pre_verdict;
            }
            trusted.verdict(post, identity, now_day)
        }
    }
}

/// A strict JSON-lines baseline store: one admitted baseline per
/// fingerprint. Duplicate fingerprints, malformed lines, and oversized
/// stores are corruption (fail closed), mirroring the tune store.
#[derive(Debug, Default)]
pub struct BaselineStore {
    baselines: BTreeMap<u64, BaselineAxes>,
}

impl BaselineStore {
    /// An empty store.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// The admitted baseline for a fingerprint, if any.
    #[must_use]
    pub fn for_fingerprint(&self, fingerprint: u64) -> Option<&BaselineAxes> {
        self.baselines.get(&fingerprint)
    }

    /// Admit a promoted baseline, REPLACING any previous baseline for
    /// the same fingerprint (updates go through promotion, so the new
    /// record carries its own provenance).
    ///
    /// # Errors
    /// Refuses structurally invalid records, rollback/non-monotone updates,
    /// and updates that would exceed the bounded store representation.
    pub fn admit(&mut self, baseline: BaselineAxes) -> Result<(), PromotionError> {
        validate_baseline(&baseline)?;
        let fingerprint = baseline.identity.fingerprint;
        if let Some(previous) = self.baselines.get(&fingerprint) {
            if previous.content_hash() == baseline.content_hash() {
                return Ok(());
            }
            if baseline.provenance.promoted_day <= previous.provenance.promoted_day {
                return Err(PromotionError {
                    detail: format!(
                        "baseline update for {fingerprint:016x} is not newer than day {}",
                        previous.provenance.promoted_day
                    ),
                });
            }
        }
        let replaced_bytes = self
            .baselines
            .get(&fingerprint)
            .map_or(0, |old| old.canonical_json().len() + 1);
        let projected =
            self.to_jsonl().len() - replaced_bytes + baseline.canonical_json().len() + 1;
        if projected > MAX_BASELINE_STORE_BYTES {
            return Err(PromotionError {
                detail: format!("baseline store exceeds the {MAX_BASELINE_STORE_BYTES}-byte bound"),
            });
        }
        self.baselines.insert(fingerprint, baseline);
        Ok(())
    }

    /// Serialize as JSON lines.
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        let mut out = String::new();
        for baseline in self.baselines.values() {
            out.push_str(&baseline.to_jsonl());
            out.push('\n');
        }
        out
    }

    /// Parse a JSON-lines store STRICTLY: every line must be a canonical
    /// baseline record; duplicate fingerprints are corruption.
    /// Structural validity and content identity detect substitution only when
    /// a caller also binds an expected hash. This parser does not authenticate
    /// operator authority; the store path remains a protected trust root.
    ///
    /// # Errors
    /// [`PromotionError`] (the store shares the promotion trust domain)
    /// naming the offending line.
    pub fn from_jsonl(text: &str) -> Result<Self, PromotionError> {
        if text.len() > MAX_BASELINE_STORE_BYTES {
            return Err(PromotionError {
                detail: format!("baseline store exceeds the {MAX_BASELINE_STORE_BYTES}-byte bound"),
            });
        }
        let mut store = BaselineStore::new();
        for (line_number, line) in text.lines().enumerate() {
            if line.is_empty() {
                return Err(PromotionError {
                    detail: format!(
                        "baseline store line {} is blank, not canonical",
                        line_number + 1
                    ),
                });
            }
            let baseline = parse_baseline_line(line).ok_or_else(|| PromotionError {
                detail: format!("baseline store line {} is not canonical", line_number + 1),
            })?;
            if store
                .for_fingerprint(baseline.identity.fingerprint)
                .is_some()
            {
                return Err(PromotionError {
                    detail: format!(
                        "baseline store line {} duplicates fingerprint {:016x}",
                        line_number + 1,
                        baseline.identity.fingerprint
                    ),
                });
            }
            validate_baseline(&baseline).map_err(|error| PromotionError {
                detail: format!(
                    "baseline store line {} is invalid: {}",
                    line_number + 1,
                    error.detail
                ),
            })?;
            store
                .baselines
                .insert(baseline.identity.fingerprint, baseline);
        }
        Ok(store)
    }
}

fn push_json_string(out: &mut String, value: &str) {
    out.push('"');
    push_json_string_body(out, value);
    out.push('"');
}

fn push_json_string_body(out: &mut String, value: &str) {
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
}

/// Minimal strict line parser for the canonical writer grammar above.
struct LineParser<'a> {
    rest: &'a str,
}

impl LineParser<'_> {
    fn take(&mut self, token: &str) -> Option<()> {
        self.rest = self.rest.strip_prefix(token)?;
        Some(())
    }

    fn string(&mut self) -> Option<String> {
        self.take("\"")?;
        let mut out = String::new();
        let mut chars = self.rest.char_indices();
        loop {
            let (index, ch) = chars.next()?;
            match ch {
                '"' => {
                    self.rest = self.rest.get(index + 1..)?;
                    if out.len() > MAX_BASELINE_STRING_BYTES {
                        return None;
                    }
                    return Some(out);
                }
                '\\' => {
                    let (_, escaped) = chars.next()?;
                    match escaped {
                        '"' => out.push('"'),
                        '\\' => out.push('\\'),
                        'n' => out.push('\n'),
                        'r' => out.push('\r'),
                        't' => out.push('\t'),
                        'u' => {
                            let mut code = 0u32;
                            for _ in 0..4 {
                                let (_, hex) = chars.next()?;
                                code = code * 16 + hex.to_digit(16)?;
                            }
                            out.push(char::from_u32(code)?);
                        }
                        _ => return None,
                    }
                }
                c if c.is_control() => return None,
                c => out.push(c),
            }
        }
    }

    fn hex_u64(&mut self) -> Option<u64> {
        let raw = self.string()?;
        if raw.len() != 16
            || !raw
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return None;
        }
        u64::from_str_radix(&raw, 16).ok()
    }

    fn content_hash(&mut self) -> Option<ContentHash> {
        let raw = self.string()?;
        if raw.len() != 64
            || !raw
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return None;
        }
        ContentHash::from_hex(&raw)
    }

    fn decimal_u64(&mut self) -> Option<u64> {
        let end = self
            .rest
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(self.rest.len());
        if end == 0 {
            return None;
        }
        let (digits, rest) = self.rest.split_at(end);
        if digits.len() > 1 && digits.starts_with('0') {
            return None; // canonical integers only
        }
        self.rest = rest;
        digits.parse().ok()
    }
}

fn parse_baseline_line(line: &str) -> Option<BaselineAxes> {
    if line.len() > MAX_BASELINE_LINE_BYTES {
        return None;
    }
    let mut p = LineParser { rest: line };
    p.take("{\"schema_version\":")?;
    let schema_version = u32::try_from(p.decimal_u64()?).ok()?;
    p.take(",\"low_band_bits\":")?;
    if p.hex_u64()? != BASELINE_LOW_BAND.to_bits() {
        return None;
    }
    p.take(",\"high_band_bits\":")?;
    if p.hex_u64()? != BASELINE_HIGH_BAND.to_bits() {
        return None;
    }
    p.take(",\"promotion_drift_bits\":")?;
    if p.hex_u64()? != MAX_AXIS_REPROBE_DRIFT.to_bits() {
        return None;
    }
    p.take(",\"fingerprint\":")?;
    let fingerprint = p.hex_u64()?;
    p.take(",\"cpu_brand\":")?;
    let cpu_brand = p.string()?;
    p.take(",\"logical_cpus\":")?;
    let logical_cpus = u32::try_from(p.decimal_u64()?).ok()?;
    p.take(",\"os\":")?;
    let os = p.string()?;
    p.take(",\"arch\":")?;
    let arch = p.string()?;
    p.take(",\"firmware\":")?;
    let firmware = p.string()?;
    p.take(",\"bandwidth_single_bits\":")?;
    let bandwidth_single_gbs = f64::from_bits(p.hex_u64()?);
    p.take(",\"bandwidth_all_core_bits\":")?;
    let bandwidth_all_core_gbs = f64::from_bits(p.hex_u64()?);
    p.take(",\"peak_single_bits\":")?;
    let peak_single_gflops = f64::from_bits(p.hex_u64()?);
    p.take(",\"peak_all_core_bits\":")?;
    let peak_all_core_gflops = f64::from_bits(p.hex_u64()?);
    p.take(",\"source_receipts\":[")?;
    let mut source_receipts = Vec::new();
    if p.rest.starts_with(']') {
        p.take("]")?;
    } else {
        loop {
            source_receipts.push(p.content_hash()?);
            if p.rest.starts_with(']') {
                p.take("]")?;
                break;
            }
            p.take(",")?;
        }
    }
    p.take(",\"promoted_by\":")?;
    let promoted_by = p.string()?;
    p.take(",\"justification\":")?;
    let justification = p.string()?;
    p.take(",\"promoted_day\":")?;
    let promoted_day = p.decimal_u64()?;
    p.take(",\"source_runs\":")?;
    let source_runs = usize::try_from(p.decimal_u64()?).ok()?;
    p.take(",\"age_policy_days\":")?;
    let age_policy_days = u32::try_from(p.decimal_u64()?).ok()?;
    p.take("}")?;
    if !p.rest.is_empty() {
        return None;
    }
    if source_runs != source_receipts.len() {
        return None;
    }
    let baseline = BaselineAxes {
        schema_version,
        identity: BaselineIdentity {
            fingerprint,
            cpu_brand,
            logical_cpus,
            os,
            arch,
            firmware,
        },
        bandwidth_single_gbs,
        bandwidth_all_core_gbs,
        peak_single_gflops,
        peak_all_core_gflops,
        provenance: BaselineProvenance {
            promoted_by,
            justification,
            promoted_day,
            source_receipts,
        },
        age_policy_days,
    };
    validate_baseline(&baseline).ok()?;
    (baseline.canonical_json() == line).then_some(baseline)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quiet_axes() -> MachineAxes {
        MachineAxes {
            fingerprint: 0xF1,
            cpu_brand: "synthetic-m".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 100.0,
            bandwidth_all_core_gbs: 220.0,
            peak_single_gflops: 45.0,
            peak_all_core_gflops: 300.0,
        }
    }

    fn identity() -> BaselineIdentity {
        BaselineIdentity::current(&quiet_axes(), "build-24F74").expect("valid identity")
    }

    fn receipt(label: &str) -> ContentHash {
        hash_domain(
            "frankensim.fs-roofline.baseline.test-source",
            label.as_bytes(),
        )
    }

    fn candidate(axes: MachineAxes, firmware: &str, label: &str) -> BaselineCandidate {
        let identity = BaselineIdentity::current(&axes, firmware).expect("candidate identity");
        BaselineCandidate::from_receipt(axes, identity, receipt(label)).expect("candidate")
    }

    fn candidates_for(axes: &MachineAxes, prefix: &str) -> Vec<BaselineCandidate> {
        (0..3)
            .map(|index| candidate(axes.clone(), "build-24F74", &format!("{prefix}-{index}")))
            .collect()
    }

    fn quiet_candidates(prefix: &str) -> Vec<BaselineCandidate> {
        candidates_for(&quiet_axes(), prefix)
    }

    fn promoted_axes_at(axes: &MachineAxes, day: u64, prefix: &str) -> BaselineAxes {
        promote_baseline(
            &candidates_for(axes, prefix),
            "operator-a",
            "three quiet-window runs on the reference host",
            day,
            90,
        )
        .expect("canonical promotion")
    }

    fn promoted_at(day: u64, prefix: &str) -> BaselineAxes {
        promoted_axes_at(&quiet_axes(), day, prefix)
    }

    fn promoted() -> BaselineAxes {
        promoted_at(20_000, "quiet")
    }

    /// QUIET DRILL: axes at the baseline are trusted end-to-end.
    #[test]
    fn quiet_axes_within_bands_are_trusted() {
        let baseline = promoted();
        let current = quiet_axes();
        assert_eq!(
            baseline.verdict(&current, &identity(), 20_010),
            BaselineVerdict::Trusted
        );
        assert!(
            citable_axis_admission(&current, &current, Some(&baseline), &identity(), 20_010)
                .trusted()
        );
    }

    /// DEGRADED DRILL — the dfh3 counterexample: 6 GB/s pre AND post on
    /// a 100 GB/s host passes floors and pre/post agreement but must be
    /// refused by the baseline.
    #[test]
    fn sustained_contention_is_refused_despite_pre_post_agreement() {
        let baseline = promoted();
        let mut crushed = quiet_axes();
        crushed.bandwidth_single_gbs = 6.0;
        crushed.bandwidth_all_core_gbs = 13.0;
        // Floors pass (6 > 5) and the pair self-agrees...
        assert!(crushed.plausible());
        assert!(crushed.reprobe_error(&crushed).is_none());
        // ...but the baseline sees a 0.06 ratio.
        let verdict =
            citable_axis_admission(&crushed, &crushed, Some(&baseline), &identity(), 20_010);
        assert!(
            matches!(
                verdict,
                BaselineVerdict::Degraded {
                    axis: "bandwidth_single_gbs",
                    ..
                }
            ),
            "{verdict:?}"
        );
        assert!(!verdict.trusted());
    }

    /// A suspiciously FAST axis is not a pass either: this is not the
    /// machine the baseline describes.
    #[test]
    fn faster_than_baseline_beyond_band_is_suspect() {
        let baseline = promoted();
        let mut upgraded = quiet_axes();
        upgraded.peak_single_gflops = 45.0 * 1.4;
        upgraded.peak_all_core_gflops = 300.0 * 1.4;
        let verdict = baseline.verdict(&upgraded, &identity(), 20_010);
        assert!(
            matches!(
                verdict,
                BaselineVerdict::Suspect {
                    axis: "peak_single_gflops",
                    ..
                }
            ),
            "{verdict:?}"
        );
    }

    /// STALE DRILL: a baseline past its age policy refuses.
    #[test]
    fn stale_baseline_is_refused_by_age_policy() {
        let baseline = promoted();
        let verdict = baseline.verdict(&quiet_axes(), &identity(), 20_000 + 91);
        assert_eq!(
            verdict,
            BaselineVerdict::Stale {
                age_days: 91,
                limit_days: 90
            }
        );
        // The day before the boundary still trusts.
        assert!(
            baseline
                .verdict(&quiet_axes(), &identity(), 20_000 + 90)
                .trusted()
        );
        let rollback = baseline.verdict(&quiet_axes(), &identity(), 19_999);
        assert_eq!(
            rollback,
            BaselineVerdict::ClockRollback {
                now_day: 19_999,
                promoted_day: 20_000,
            }
        );
        assert!(!rollback.trusted());
    }

    /// FIRMWARE-DRIFT DRILL: a changed firmware declaration is identity
    /// drift, refused before any band math.
    #[test]
    fn firmware_drift_is_identity_refusal() {
        let baseline = promoted();
        let moved =
            BaselineIdentity::current(&quiet_axes(), "build-25A01").expect("moved identity");
        assert_eq!(
            baseline.verdict(&quiet_axes(), &moved, 20_010),
            BaselineVerdict::IdentityDrift { field: "firmware" }
        );
        let mut other_machine = quiet_axes();
        other_machine.fingerprint = 0xF2;
        assert_eq!(
            baseline.verdict(&other_machine, &identity(), 20_010),
            BaselineVerdict::IdentityDrift {
                field: "fingerprint"
            }
        );
    }

    /// FIRST-RUN LAW: no baseline → Unbaselined, never Trusted; the
    /// measurement is candidate evidence and cannot authorize itself.
    #[test]
    fn first_run_measurements_are_candidates_not_baselines() {
        let current = quiet_axes();
        let verdict = citable_axis_admission(&current, &current, None, &identity(), 20_010);
        assert_eq!(verdict, BaselineVerdict::Unbaselined);
        assert!(!verdict.trusted());
    }

    #[test]
    fn coarse_refusals_are_structured_and_always_valid_json() {
        let baseline = promoted();
        let mut invalid = quiet_axes();
        invalid.bandwidth_single_gbs = f64::NAN;
        let invalid_verdict = citable_axis_admission(
            &invalid,
            &quiet_axes(),
            Some(&baseline),
            &identity(),
            20_010,
        );
        assert!(matches!(
            invalid_verdict,
            BaselineVerdict::InvalidAxes { probe: "pre", .. }
        ));
        let invalid_json = invalid_verdict.to_jsonl();
        assert!(!invalid_json.contains("NaN") && !invalid_json.contains("inf"));

        let mut drifted = quiet_axes();
        drifted.bandwidth_single_gbs *= 0.5;
        let reprobe = citable_axis_admission(
            &quiet_axes(),
            &drifted,
            Some(&baseline),
            &identity(),
            20_010,
        );
        assert!(matches!(reprobe, BaselineVerdict::ReprobeFailed { .. }));
        assert!(reprobe.to_jsonl().starts_with('{'));

        let defensive = BaselineVerdict::Degraded {
            axis: "synthetic",
            ratio: f64::NAN,
        }
        .to_jsonl();
        assert!(defensive.contains("\"ratio\":null"));
    }

    /// GOVERNED PROMOTION: blank operator/justification, too few runs,
    /// disagreeing runs, foreign-fingerprint runs, and out-of-policy age
    /// are each refused with a teaching detail.
    #[test]
    fn promotion_is_governed_and_fails_closed() {
        let runs = quiet_candidates("governed");
        let refuse = |result: Result<BaselineAxes, PromotionError>, needle: &str| {
            let err = result.expect_err(needle);
            assert!(err.detail.contains(needle), "{}", err.detail);
        };
        refuse(
            promote_baseline(&runs, "  ", "why", 1, 90),
            "named operator",
        );
        refuse(promote_baseline(&runs, "op", " ", 1, 90), "justification");
        refuse(
            promote_baseline(&runs[..2], "op", "why", 1, 90),
            "at least 3",
        );
        refuse(promote_baseline(&runs, "op", "why", 1, 0), "age policy");
        refuse(promote_baseline(&runs, "op", "why", 1, 9999), "age policy");
        // One crushed run among quiet ones: mutual agreement refuses.
        let mut crushed = quiet_axes();
        crushed.bandwidth_single_gbs = 6.0;
        crushed.bandwidth_all_core_gbs = 13.0;
        let mixed = vec![
            candidate(quiet_axes(), "build-24F74", "mixed-0"),
            candidate(quiet_axes(), "build-24F74", "mixed-1"),
            candidate(crushed, "build-24F74", "mixed-2"),
        ];
        refuse(promote_baseline(&mixed, "op", "why", 1, 90), "disagree");
        // A foreign full environment identity cannot join a promotion.
        let mut foreign_axes = quiet_axes();
        foreign_axes.fingerprint = 0xF2;
        let foreign = vec![
            candidate(quiet_axes(), "build-24F74", "foreign-0"),
            candidate(foreign_axes, "build-24F74", "foreign-1"),
            candidate(quiet_axes(), "build-24F74", "foreign-2"),
        ];
        refuse(
            promote_baseline(&foreign, "op", "why", 1, 90),
            "environment identity",
        );
        // An implausible run cannot become a candidate.
        let mut implausible = quiet_axes();
        implausible.peak_single_gflops = f64::NAN;
        let implausible_identity =
            BaselineIdentity::current(&implausible, "build-24F74").expect("identity fields valid");
        refuse(
            BaselineCandidate::from_receipt(
                implausible,
                implausible_identity,
                receipt("implausible"),
            )
            .map(|_| promoted()),
            "plausibility floors",
        );
        // Three array entries must name three retained receipts.
        let duplicate_receipt = receipt("duplicate");
        let duplicate = (0..3)
            .map(|_| {
                BaselineCandidate::from_receipt(quiet_axes(), identity(), duplicate_receipt)
                    .expect("candidate")
            })
            .collect::<Vec<_>>();
        refuse(
            promote_baseline(&duplicate, "op", "why", 1, 90),
            "reuses a source receipt",
        );
    }

    #[test]
    fn canonical_hash_binds_schema_axes_policy_and_source_receipts() {
        let baseline = promoted();
        assert_eq!(baseline.schema_version(), BASELINE_SCHEMA_VERSION);
        assert_eq!(baseline.identity().fingerprint(), 0xF1);
        assert_eq!(baseline.provenance().source_runs(), 3);
        assert!(
            baseline
                .provenance()
                .source_receipts()
                .windows(2)
                .all(|pair| pair[0] < pair[1])
        );
        assert_eq!(
            baseline.content_hash(),
            hash_domain(BASELINE_HASH_DOMAIN, baseline.canonical_json().as_bytes())
        );

        let mut reversed = quiet_candidates("order");
        reversed.reverse();
        let reversed = promote_baseline(
            &reversed,
            "operator-a",
            "three quiet-window runs on the reference host",
            20_000,
            90,
        )
        .expect("order-independent promotion");
        let ordered = promote_baseline(
            &quiet_candidates("order"),
            "operator-a",
            "three quiet-window runs on the reference host",
            20_000,
            90,
        )
        .expect("ordered promotion");
        assert_eq!(ordered.content_hash(), reversed.content_hash());
        assert_ne!(baseline.content_hash(), ordered.content_hash());
    }

    /// Store round-trip is lossless; corruption and duplicates refuse.
    #[test]
    fn store_round_trips_and_fails_closed() {
        let mut store = BaselineStore::new();
        store.admit(promoted()).expect("admit");
        let text = store.to_jsonl();
        let back = BaselineStore::from_jsonl(&text).expect("canonical store parses");
        assert_eq!(back.for_fingerprint(0xF1), Some(&promoted()));
        assert!(back.for_fingerprint(0xF2).is_none());
        // Tampered line: refused.
        assert!(BaselineStore::from_jsonl(&text.replace("operator-a", "")).is_err());
        assert!(BaselineStore::from_jsonl("{\"not\":\"a baseline\"}\n").is_err());
        // Duplicate fingerprint: corruption.
        let duplicated = format!("{text}{text}");
        assert!(BaselineStore::from_jsonl(&duplicated).is_err());
        assert!(BaselineStore::from_jsonl(&format!("\n{text}")).is_err());

        // Semantics that promotion could not produce are refused on load.
        let implausible = text.replace(
            &format!("{:016x}", 100.0f64.to_bits()),
            &format!("{:016x}", 1.0f64.to_bits()),
        );
        assert!(BaselineStore::from_jsonl(&implausible).is_err());
        assert!(
            BaselineStore::from_jsonl(&text.replace("\"source_runs\":3", "\"source_runs\":4"))
                .is_err()
        );
        assert!(
            BaselineStore::from_jsonl(&text.replace(
                &format!("{:016x}", BASELINE_LOW_BAND.to_bits()),
                &format!("{:016x}", 0.5f64.to_bits()),
            ))
            .is_err(),
            "serialized admission policy is part of the baseline identity"
        );
        assert!(
            BaselineStore::from_jsonl(&text.replace("operator-a", "oper\\u0061tor-a")).is_err(),
            "alternate JSON encodings are not canonical"
        );

        // Updates are idempotent by hash and otherwise strictly monotone.
        store.admit(promoted()).expect("idempotent admit");
        let older = promoted_at(19_999, "older");
        assert!(store.admit(older).is_err());
        let refreshed = promoted_at(21_000, "refreshed");
        store.admit(refreshed.clone()).expect("newer promotion");
        assert_eq!(store.for_fingerprint(0xF1), Some(&refreshed));
        assert_eq!(store.to_jsonl().lines().count(), 1);

        // The store cannot be bypassed with an invalid in-memory record.
        let mut invalid = promoted_at(22_000, "invalid");
        invalid.bandwidth_single_gbs = 1.0;
        assert!(store.admit(invalid).is_err());

        // Canonical store order is a function of content, not admission order.
        let mut second_axes = quiet_axes();
        second_axes.fingerprint = 0xF2;
        let first = promoted_at(20_000, "deterministic-first");
        let second = promoted_axes_at(&second_axes, 20_000, "deterministic-second");
        let mut forward = BaselineStore::new();
        forward.admit(first.clone()).expect("first");
        forward.admit(second.clone()).expect("second");
        let mut reverse = BaselineStore::new();
        reverse.admit(second).expect("second");
        reverse.admit(first).expect("first");
        assert_eq!(forward.to_jsonl(), reverse.to_jsonl());
    }

    #[test]
    fn identity_and_serialization_bounds_fail_closed() {
        assert!(BaselineIdentity::current(&quiet_axes(), " ").is_err());
        assert!(BaselineIdentity::current(&quiet_axes(), "firmware\0").is_err());
        assert!(BaselineIdentity::current(&quiet_axes(), "x".repeat(4097)).is_err());

        let candidates = quiet_candidates("bounds");
        assert!(promote_baseline(&candidates, "operator-a", "x".repeat(4097), 20_000, 90).is_err());

        let mut huge_brand = quiet_axes();
        huge_brand.cpu_brand = "x".repeat(MAX_BASELINE_STRING_BYTES);
        let candidates = (0..3)
            .map(|index| {
                candidate(
                    huge_brand.clone(),
                    "y".repeat(4096).as_str(),
                    &format!("huge-{index}"),
                )
            })
            .collect::<Vec<_>>();
        let error = promote_baseline(&candidates, "z".repeat(4096), "q".repeat(4096), 20_000, 90)
            .expect_err("aggregate line bound");
        assert!(error.detail.contains("line bound"), "{}", error.detail);
    }
    // ---- fz2.7: promotion-authority drills -------------------------------

    fn authority_with(key: &str) -> crate::StaticKeyRegistry {
        let mut registry = crate::StaticKeyRegistry::new();
        registry.authorize(key);
        registry
    }

    fn attest(baseline: &BaselineAxes, key: &str) -> PromotionAttestation {
        PromotionAttestation::new(
            key,
            crate::StaticKeyRegistry::tag(key, &baseline.promotion_message()),
        )
    }

    fn retained(baseline: &BaselineAxes) -> BTreeSet<ContentHash> {
        baseline
            .provenance()
            .source_receipts()
            .iter()
            .copied()
            .collect()
    }

    /// Authorized admission round-trips through the attested store and
    /// re-verifies identically after reload (deterministic replay).
    #[test]
    fn authorized_admission_round_trips_and_replays() {
        let baseline = promoted();
        let authority = authority_with("ops/2026-q3");
        let attestation = attest(&baseline, "ops/2026-q3");
        let mut store = AttestedBaselineStore::new();
        store
            .admit_verified(
                baseline.clone(),
                attestation.clone(),
                &authority,
                &retained(&baseline),
            )
            .expect("authorized admission");
        let text = store.to_jsonl();
        let back = AttestedBaselineStore::from_jsonl(&text).expect("attested store parses");
        let fingerprint = baseline.identity().fingerprint();
        assert_eq!(back.for_fingerprint(fingerprint), Some(&baseline));
        assert_eq!(back.attestation_for(fingerprint), Some(&attestation));
        // Citable use re-verifies the reloaded record end-to-end.
        let verdict = citable_axis_admission_authorized(
            &quiet_axes(),
            &quiet_axes(),
            back.for_fingerprint(fingerprint),
            back.attestation_for(fingerprint),
            &identity(),
            20_010,
            &authority,
        );
        assert_eq!(verdict, BaselineVerdict::Trusted);
        assert_eq!(back.to_jsonl(), text, "replay is byte-identical");
    }

    /// FORGED OPERATOR / EDITED RECORD: the signature covers the content
    /// hash, so editing the free-text operator (or anything else) makes
    /// the attestation WrongSignature.
    #[test]
    fn forged_operator_and_edited_record_invalidate_the_attestation() {
        let signed = promoted();
        let authority = authority_with("ops/2026-q3");
        let attestation = attest(&signed, "ops/2026-q3");
        // Same candidates, different operator: a forged promoted_by.
        let forged = promote_baseline(
            &quiet_candidates("quiet"),
            "operator-EVIL",
            "three quiet-window runs on the reference host",
            20_000,
            90,
        )
        .expect("structurally valid");
        assert_eq!(
            forged.authority_verdict(Some(&attestation), &authority),
            KeyVerdict::WrongSignature,
            "operator edits move the signed hash"
        );
        // Edited record: different justification, same everything else.
        let edited = promote_baseline(
            &quiet_candidates("quiet"),
            "operator-a",
            "EDITED justification",
            20_000,
            90,
        )
        .expect("structurally valid");
        assert_eq!(
            edited.authority_verdict(Some(&attestation), &authority),
            KeyVerdict::WrongSignature
        );
        let mut store = AttestedBaselineStore::new();
        let refused = store
            .admit_verified(forged, attestation, &authority, &retained(&signed))
            .expect_err("forged record refused");
        assert!(refused.detail.contains("wrong-signature"), "{refused}");
    }

    /// MISSING SOURCE: a promotion whose named receipt is not retained
    /// is refused at admission, naming the hash.
    #[test]
    fn missing_source_receipt_refuses_admission() {
        let baseline = promoted();
        let authority = authority_with("ops/2026-q3");
        let attestation = attest(&baseline, "ops/2026-q3");
        let mut available = retained(&baseline);
        let dropped = *baseline
            .provenance()
            .source_receipts()
            .first()
            .expect("receipts");
        available.remove(&dropped);
        let refused = AttestedBaselineStore::new()
            .admit_verified(baseline, attestation, &authority, &available)
            .expect_err("missing source refused");
        assert!(refused.detail.contains(&dropped.to_hex()), "{refused}");
    }

    /// DUPLICATE SOURCE: two candidates sharing one retained receipt
    /// cannot corroborate each other — refused at promotion.
    #[test]
    fn duplicate_source_receipts_refuse_promotion() {
        let mut candidates = quiet_candidates("dup");
        candidates[1] = candidate(quiet_axes(), "build-24F74", "dup-0"); // same as [0]
        let refused = promote_baseline(&candidates, "operator-a", "why", 20_000, 90)
            .expect_err("duplicate receipts refused");
        assert!(refused.detail.contains("reuses"), "{refused}");
    }

    /// WRONG KEY: a valid tag claimed under a different key id fails —
    /// authorized keys cannot vouch for each other.
    #[test]
    fn wrong_key_claims_are_refused() {
        let baseline = promoted();
        let mut authority = crate::StaticKeyRegistry::new();
        authority.authorize("ops/a");
        authority.authorize("ops/b");
        let tag_a = crate::StaticKeyRegistry::tag("ops/a", &baseline.promotion_message());
        let claimed_as_b = PromotionAttestation::new("ops/b", tag_a);
        assert_eq!(
            baseline.authority_verdict(Some(&claimed_as_b), &authority),
            KeyVerdict::WrongSignature
        );
        let unknown = attest(&baseline, "ops/never-registered");
        assert_eq!(
            baseline.authority_verdict(Some(&unknown), &authority),
            KeyVerdict::UnknownKey
        );
        // Unattested is never authorized, under ANY authority.
        assert_eq!(
            baseline.authority_verdict(None, &authority),
            KeyVerdict::UnknownKey
        );
        assert_eq!(
            baseline.authority_verdict(None, &crate::NoPromotionAuthority),
            KeyVerdict::UnknownKey
        );
    }

    /// REVOKED KEY + VALID ROTATION: revocation retroactively demands
    /// re-promotion; a rotation to a newly authorized key re-verifies.
    #[test]
    fn revocation_demands_repromotion_and_rotation_recovers() {
        let baseline = promoted();
        let mut authority = authority_with("ops/2026-q3");
        let old = attest(&baseline, "ops/2026-q3");
        let mut store = AttestedBaselineStore::new();
        store
            .admit_verified(
                baseline.clone(),
                old.clone(),
                &authority,
                &retained(&baseline),
            )
            .expect("initially authorized");
        authority.revoke("ops/2026-q3");
        let verdict = citable_axis_admission_authorized(
            &quiet_axes(),
            &quiet_axes(),
            store.for_fingerprint(baseline.identity().fingerprint()),
            store.attestation_for(baseline.identity().fingerprint()),
            &identity(),
            20_010,
            &authority,
        );
        assert_eq!(
            verdict,
            BaselineVerdict::Unauthorized {
                verdict: "revoked-key"
            },
            "revoked keys stop citable use BEFORE any band math"
        );
        // Rotation: authorize the new key, RE-PROMOTE (newer day), re-attest.
        authority.authorize("ops/2026-q4");
        let rotated = promoted_at(20_005, "rotated");
        let fresh = attest(&rotated, "ops/2026-q4");
        store
            .admit_verified(rotated.clone(), fresh, &authority, &retained(&rotated))
            .expect("rotated admission");
        let verdict = citable_axis_admission_authorized(
            &quiet_axes(),
            &quiet_axes(),
            store.for_fingerprint(rotated.identity().fingerprint()),
            store.attestation_for(rotated.identity().fingerprint()),
            &identity(),
            20_010,
            &authority,
        );
        assert_eq!(verdict, BaselineVerdict::Trusted);
    }

    /// Tampered attested-store lines refuse at parse or at re-verify.
    #[test]
    fn tampered_attested_lines_fail_closed() {
        let baseline = promoted();
        let authority = authority_with("ops/2026-q3");
        let attestation = attest(&baseline, "ops/2026-q3");
        let mut store = AttestedBaselineStore::new();
        store
            .admit_verified(
                baseline.clone(),
                attestation,
                &authority,
                &retained(&baseline),
            )
            .expect("authorized admission");
        let text = store.to_jsonl();
        // Record tamper: the line still PARSES but citable re-verify refuses.
        let tampered = text.replace("operator-a", "operator-b");
        let back = AttestedBaselineStore::from_jsonl(&tampered).expect("parses structurally");
        let fingerprint = baseline.identity().fingerprint();
        let verdict = citable_axis_admission_authorized(
            &quiet_axes(),
            &quiet_axes(),
            back.for_fingerprint(fingerprint),
            back.attestation_for(fingerprint),
            &identity(),
            20_010,
            &authority,
        );
        assert_eq!(
            verdict,
            BaselineVerdict::Unauthorized {
                verdict: "wrong-signature"
            }
        );
        // Envelope tamper: refused at parse.
        assert!(
            AttestedBaselineStore::from_jsonl(&text.replace("\"attestation\"", "\"x\"")).is_err()
        );
    }
}
