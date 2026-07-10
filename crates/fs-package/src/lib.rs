//! fs-package — machine-checkable evidence packages (plan addendum,
//! Proposal 12). Layer: L6.
//!
//! When FrankenSim asserts "this design meets spec", the assertion travels as a
//! self-contained, CONTENT-ADDRESSED bundle: the color-typed claims, the raw
//! certificate data behind each (carried in the [`fs_evidence::Color`]),
//! provenance (code version + the constellation lockfile), and a Merkle root
//! over the package identity so any tamper is detectable. A standalone,
//! open-source CHECKER re-verifies the package WITHOUT re-running any solver —
//! that structural re-verification is exactly [`EvidencePackage::verify`].
//!
//! Completeness is enforced, not assumed: a validated-color claim that is
//! missing its regime tag OR its anchoring dataset FAILS verification (an
//! unfalsifiable "validated" claim is worse than none). An all-estimated
//! package is still valid and round-trips — honesty about low confidence is
//! not a defect.
//!
//! The Merkle tree uses an in-house FNV-1a content hash (pure Rust, std only —
//! Franken-compliant); a cryptographic signature is DETACHED and OPTIONAL (the
//! bundle is verifiable by content address regardless). Everything is
//! deterministic: the same package yields the same root and JSON.

use fs_evidence::{Color, ColorRank, IntervalOp, compose};

pub mod coverage;
pub use coverage::{ConceptPresence, CoverageStatus, package_coverage, package_presence};

/// A COMPOSITION RECEIPT (schema v3, bead xfxq): this claim's color was
/// derived from earlier claims in the package, and the standalone
/// checker re-runs the derivation — `compose` folded over the parents'
/// colors in order must EQUAL the claimed color exactly. Parents are
/// indices into the package's claim list and must precede this claim
/// (a DAG by construction).
#[derive(Debug, Clone, PartialEq)]
pub struct CompositionReceipt {
    /// Parent claim indices, in fold order (each < this claim's index).
    pub parents: Vec<usize>,
    /// The ledger operation the derivation used.
    pub op: IntervalOp,
}

/// One falsifier's adversarial record against a claim (schema v3):
/// negative results travel WITH the claim; a refuted claim fails
/// verification outright.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FalsifierRecord {
    /// Which falsifier ran.
    pub name: String,
    /// Adversarial attempts executed.
    pub attempts: u64,
    /// Did it refute the claim?
    pub refuted: bool,
    /// Outcome summary.
    pub detail: String,
}

/// An anchoring-dataset identity (schema v3): the reference data behind
/// a validated claim, by stable id and content hash — not just a name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnchorRecord {
    /// Stable dataset identity.
    pub dataset_id: String,
    /// Content hash (hex) of the dataset artifact.
    pub content_hash: String,
}

/// One claim in an evidence package: a statement plus its epistemic color
/// (which carries the certificate data — an interval, a regime+dataset, or an
/// estimator+dispersion), and optionally its composition receipt,
/// falsifier records, and dataset anchors (schema v3).
#[derive(Debug, Clone, PartialEq)]
pub struct Claim {
    /// A stable claim id.
    pub id: String,
    /// The human-readable claim.
    pub statement: String,
    /// The epistemic color + its certificate payload.
    pub color: Color,
    /// The composition receipt, when this claim is derived.
    pub receipt: Option<CompositionReceipt>,
    /// Falsifier records (adversarial attempts + outcomes).
    pub falsifiers: Vec<FalsifierRecord>,
    /// Anchoring datasets.
    pub anchors: Vec<AnchorRecord>,
}

impl Claim {
    /// A claim.
    #[must_use]
    pub fn new(id: impl Into<String>, statement: impl Into<String>, color: Color) -> Claim {
        Claim {
            id: id.into(),
            statement: statement.into(),
            color,
            receipt: None,
            falsifiers: Vec::new(),
            anchors: Vec::new(),
        }
    }

    /// Attach a composition receipt (builder style).
    #[must_use]
    pub fn with_receipt(mut self, parents: Vec<usize>, op: IntervalOp) -> Claim {
        self.receipt = Some(CompositionReceipt { parents, op });
        self
    }

    /// Attach a falsifier record (builder style).
    #[must_use]
    pub fn with_falsifier(mut self, rec: FalsifierRecord) -> Claim {
        self.falsifiers.push(rec);
        self
    }

    /// Attach a dataset anchor (builder style).
    #[must_use]
    pub fn with_anchor(
        mut self,
        dataset_id: impl Into<String>,
        content_hash: impl Into<String>,
    ) -> Claim {
        self.anchors.push(AnchorRecord {
            dataset_id: dataset_id.into(),
            content_hash: content_hash.into(),
        });
        self
    }

    /// A canonical string used for content hashing (bit-exact on floats).
    fn canonical(&self) -> String {
        use core::fmt::Write as _;
        let mut out = String::from("claim|");
        push_atom(&mut out, &self.id);
        push_atom(&mut out, &self.statement);
        match &self.color {
            Color::Verified { lo, hi } => {
                out.push_str("verified|");
                let _ = write!(out, "{}|{}|", lo.to_bits(), hi.to_bits());
            }
            Color::Validated { regime, dataset } => {
                out.push_str("validated|");
                for (k, (lo, hi)) in regime.bounds() {
                    push_atom(&mut out, k);
                    let _ = write!(out, "{}|{}|", lo.to_bits(), hi.to_bits());
                }
                push_atom(&mut out, dataset);
            }
            Color::Estimated {
                estimator,
                dispersion,
            } => {
                out.push_str("estimated|");
                push_atom(&mut out, estimator);
                let _ = write!(out, "{}|", dispersion.to_bits());
            }
        }
        // Schema-v3 fields bind into the content address too.
        match &self.receipt {
            Some(r) => {
                let _ = write!(out, "receipt:{}|", op_name(r.op));
                for &p in &r.parents {
                    let _ = write!(out, "{p}|");
                }
            }
            None => out.push_str("no-receipt|"),
        }
        for fr in &self.falsifiers {
            out.push_str("falsifier|");
            push_atom(&mut out, &fr.name);
            let _ = write!(out, "{}|{}|", fr.attempts, fr.refuted);
            push_atom(&mut out, &fr.detail);
        }
        for a in &self.anchors {
            out.push_str("anchor|");
            push_atom(&mut out, &a.dataset_id);
            push_atom(&mut out, &a.content_hash);
        }
        out
    }
}

/// Stable op name for hashing/JSON.
fn op_name(op: IntervalOp) -> &'static str {
    match op {
        IntervalOp::Add => "add",
        IntervalOp::Mul => "mul",
        IntervalOp::Hull => "hull",
    }
}

fn op_parse(name: &str) -> Option<IntervalOp> {
    match name {
        "add" => Some(IntervalOp::Add),
        "mul" => Some(IntervalOp::Mul),
        "hull" => Some(IntervalOp::Hull),
        _ => None,
    }
}

/// Where a package came from — enough to reproduce it.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Provenance {
    /// The code version / commit that produced the claims.
    pub code_version: String,
    /// The pinned dependency constellation (lockfile digest).
    pub constellation_lock: String,
}

impl Provenance {
    /// Provenance.
    #[must_use]
    pub fn new(
        code_version: impl Into<String>,
        constellation_lock: impl Into<String>,
    ) -> Provenance {
        Provenance {
            code_version: code_version.into(),
            constellation_lock: constellation_lock.into(),
        }
    }
}

/// A self-contained, content-addressed evidence bundle.
#[derive(Debug, Clone, PartialEq)]
pub struct EvidencePackage {
    /// The format version (stability promise for external checkers).
    pub format_version: u32,
    /// The claims, in order.
    pub claims: Vec<Claim>,
    /// Provenance.
    pub provenance: Provenance,
    /// An OPTIONAL detached signature over the Merkle root.
    pub signature: Option<String>,
}

/// The by-color budget pie over a package's claims.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ColorBreakdown {
    /// Verified-color claims.
    pub verified: usize,
    /// Validated-color claims.
    pub validated: usize,
    /// Estimated-color claims.
    pub estimated: usize,
}

/// The result of verifying a package.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackageReport {
    /// The recomputed content address (Merkle root).
    pub merkle_root: u64,
    /// The by-color budget pie.
    pub breakdown: ColorBreakdown,
    /// The number of claims.
    pub claims: usize,
}

/// A structured verification failure.
#[derive(Debug, Clone, PartialEq)]
pub enum PackageError {
    /// Required reproducibility provenance is blank.
    IncompleteProvenance {
        /// What is missing (`"code_version"` or `"constellation_lock"`).
        missing: &'static str,
    },
    /// A claim id is blank or duplicates an earlier claim id.
    InvalidClaimId {
        /// The claim's position in the package.
        index: usize,
        /// The invalid id (blank for the blank-id case).
        id: String,
        /// Why it is invalid (`"blank"` or `"duplicate"`).
        reason: &'static str,
    },
    /// A validated claim is missing part of its evidence.
    IncompleteValidatedClaim {
        /// The claim id.
        claim: String,
        /// What is missing (`"regime"` or `"dataset"`).
        missing: &'static str,
    },
    /// A verified claim's certificate interval is not a finite `[lo <= hi]`.
    IncompleteVerifiedClaim {
        /// The claim id.
        claim: String,
    },
    /// A validated claim has a malformed validity-domain axis.
    InvalidValidatedRegime {
        /// The claim id.
        claim: String,
        /// The malformed axis name (blank for a blank name).
        axis: String,
    },
    /// An estimated claim is missing its estimator identity.
    IncompleteEstimatedClaim {
        /// The claim id.
        claim: String,
        /// What is missing (`"estimator"`).
        missing: &'static str,
    },
    /// An estimated claim's dispersion is NaN or negative. Positive infinity
    /// is the lower-layer algebra's explicit no-quantitative-claim sentinel.
    InvalidEstimatedDispersion {
        /// The claim id.
        claim: String,
    },
    /// Finite claim magnitudes overflowed while deriving the package budget.
    MagnitudeOverflow {
        /// Claim at which the finite subtotal became non-finite.
        claim: String,
        /// Budget component (`"verified_width"` or `"estimated_dispersion"`).
        component: &'static str,
    },
    /// The in-memory package cannot fit the standalone checker's bounded
    /// transport envelope.
    TransportLimit {
        /// Field or container that exceeded its limit.
        what: String,
        /// Configured upper bound.
        limit: usize,
    },
    /// The declared format version is unsupported.
    UnsupportedFormat {
        /// The version found.
        found: u32,
    },
    /// A composition receipt does not re-derive the claimed color: the
    /// checker re-ran `compose` over the parents and got a different
    /// result — a forged or stale derivation (schema v3).
    ReceiptMismatch {
        /// The claim id.
        claim: String,
    },
    /// A receipt references a parent at or after the claim itself (the
    /// derivation DAG must point strictly backwards), or out of range.
    BadReceiptParent {
        /// The claim id.
        claim: String,
        /// The offending parent index.
        parent: usize,
    },
    /// A falsifier REFUTED this claim; a refuted claim cannot verify.
    RefutedClaim {
        /// The claim id.
        claim: String,
        /// The refuting falsifier.
        falsifier: String,
    },
}

/// The one format version this build understands. v2 (bead qmao.6.1)
/// added complete color payloads + the strict parser + root
/// recomputation; v3 (bead xfxq) adds composition receipts (checker
/// re-runs the derivation), falsifier records (refuted claims fail),
/// and dataset anchors — all bound into the content address.
pub const FORMAT_VERSION: u32 = 3;

impl EvidencePackage {
    /// An empty package at the current format version.
    #[must_use]
    pub fn new(provenance: Provenance) -> EvidencePackage {
        EvidencePackage {
            format_version: FORMAT_VERSION,
            claims: Vec::new(),
            provenance,
            signature: None,
        }
    }

    /// Add a claim (builder style).
    #[must_use]
    pub fn with_claim(mut self, claim: Claim) -> EvidencePackage {
        self.claims.push(claim);
        self
    }

    /// Attach a detached signature (builder style).
    #[must_use]
    pub fn signed(mut self, signature: impl Into<String>) -> EvidencePackage {
        self.signature = Some(signature.into());
        self
    }

    /// The content address: an FNV-1a Merkle root over the package identity
    /// (format version, provenance, and ordered claims). Detached signatures are
    /// excluded so signing does not change the address.
    #[must_use]
    pub fn merkle_root(&self) -> u64 {
        let mut level: Vec<u64> = Vec::with_capacity(self.claims.len() + 1);
        level.push(fnv1a(self.package_header().as_bytes()));
        level.extend(self.claims.iter().map(|c| fnv1a(c.canonical().as_bytes())));
        while level.len() > 1 {
            let mut next = Vec::with_capacity(level.len().div_ceil(2));
            for pair in level.chunks(2) {
                match pair {
                    [a, b] => next.push(combine(*a, *b)),
                    [a] => next.push(*a), // odd node carries up
                    _ => {}
                }
            }
            level = next;
        }
        match level.as_slice() {
            [root] => *root,
            [] => fnv1a(b"fs-package:empty-internal-level"),
            _ => fnv1a(b"fs-package:invalid-internal-level"),
        }
    }

    fn package_header(&self) -> String {
        use core::fmt::Write as _;
        let mut out = String::from("package|");
        let _ = write!(
            out,
            "format:{}|claims:{}|",
            self.format_version,
            self.claims.len()
        );
        push_atom(&mut out, &self.provenance.code_version);
        push_atom(&mut out, &self.provenance.constellation_lock);
        out
    }

    /// The by-color budget pie over the claims.
    #[must_use]
    pub fn color_breakdown(&self) -> ColorBreakdown {
        let mut b = ColorBreakdown::default();
        for c in &self.claims {
            match c.color.rank() {
                ColorRank::Verified => b.verified += 1,
                ColorRank::Validated => b.validated += 1,
                ColorRank::Estimated => b.estimated += 1,
            }
        }
        b
    }

    fn verify_claim(&self, index: usize, claim: &Claim) -> Result<(), PackageError> {
        // Schema-v3 semantic re-verification (solver-free): refuted falsifiers
        // fail and composition receipts are independently re-derived.
        if let Some(fr) = claim.falsifiers.iter().find(|f| f.refuted) {
            return Err(PackageError::RefutedClaim {
                claim: claim.id.clone(),
                falsifier: fr.name.clone(),
            });
        }
        if let Some(receipt) = &claim.receipt {
            let mut derived: Option<Color> = None;
            for &parent in &receipt.parents {
                if parent >= index {
                    return Err(PackageError::BadReceiptParent {
                        claim: claim.id.clone(),
                        parent,
                    });
                }
                let parent_color = &self.claims[parent].color;
                derived = Some(match derived {
                    None => parent_color.clone(),
                    Some(current) => compose(&current, parent_color, receipt.op),
                });
            }
            if !matches!(derived, Some(color) if color == claim.color) {
                return Err(PackageError::ReceiptMismatch {
                    claim: claim.id.clone(),
                });
            }
        }
        match &claim.color {
            Color::Verified { lo, hi } => {
                if !(lo.is_finite() && hi.is_finite() && lo <= hi) {
                    return Err(PackageError::IncompleteVerifiedClaim {
                        claim: claim.id.clone(),
                    });
                }
            }
            Color::Validated { regime, dataset } => {
                if regime.bounds().is_empty() {
                    return Err(PackageError::IncompleteValidatedClaim {
                        claim: claim.id.clone(),
                        missing: "regime",
                    });
                }
                if dataset.trim().is_empty() {
                    return Err(PackageError::IncompleteValidatedClaim {
                        claim: claim.id.clone(),
                        missing: "dataset",
                    });
                }
                if let Some((axis, _)) = regime.bounds().iter().find(|(axis, (lo, hi))| {
                    axis.trim().is_empty() || !lo.is_finite() || !hi.is_finite() || lo > hi
                }) {
                    return Err(PackageError::InvalidValidatedRegime {
                        claim: claim.id.clone(),
                        axis: axis.clone(),
                    });
                }
            }
            Color::Estimated {
                estimator,
                dispersion,
            } => {
                if estimator.trim().is_empty() {
                    return Err(PackageError::IncompleteEstimatedClaim {
                        claim: claim.id.clone(),
                        missing: "estimator",
                    });
                }
                if dispersion.is_nan() || *dispersion < 0.0 {
                    return Err(PackageError::InvalidEstimatedDispersion {
                        claim: claim.id.clone(),
                    });
                }
            }
        }
        Ok(())
    }

    /// Re-verify the package WITHOUT any solver: the format must be supported
    /// and every claim's certificate must be complete for its color. Returns
    /// the content address + budget pie on success.
    ///
    /// # Errors
    /// [`PackageError`] on an unsupported format or an incomplete claim.
    pub fn verify(&self) -> Result<PackageReport, PackageError> {
        if self.format_version != FORMAT_VERSION {
            return Err(PackageError::UnsupportedFormat {
                found: self.format_version,
            });
        }
        self.verify_transport_limits()?;
        if self.provenance.code_version.trim().is_empty() {
            return Err(PackageError::IncompleteProvenance {
                missing: "code_version",
            });
        }
        if self.provenance.constellation_lock.trim().is_empty() {
            return Err(PackageError::IncompleteProvenance {
                missing: "constellation_lock",
            });
        }
        let mut claim_ids = std::collections::BTreeSet::new();
        for (index, c) in self.claims.iter().enumerate() {
            if c.id.trim().is_empty() {
                return Err(PackageError::InvalidClaimId {
                    index,
                    id: c.id.clone(),
                    reason: "blank",
                });
            }
            if !claim_ids.insert(c.id.as_str()) {
                return Err(PackageError::InvalidClaimId {
                    index,
                    id: c.id.clone(),
                    reason: "duplicate",
                });
            }
            self.verify_claim(index, c)?;
        }
        self.verify_finite_magnitude_sums()?;
        Ok(PackageReport {
            merkle_root: self.merkle_root(),
            breakdown: self.color_breakdown(),
            claims: self.claims.len(),
        })
    }

    fn verify_transport_limits(&self) -> Result<(), PackageError> {
        fn check_count(what: &str, count: usize) -> Result<(), PackageError> {
            if count > MAX_JSON_CONTAINER_ITEMS {
                return Err(PackageError::TransportLimit {
                    what: what.to_string(),
                    limit: MAX_JSON_CONTAINER_ITEMS,
                });
            }
            Ok(())
        }

        fn escaped_len(value: &str) -> usize {
            value
                .chars()
                .map(|ch| match ch {
                    '"' | '\\' | '\n' | '\r' | '\t' => 2,
                    c if c.is_control() => 6,
                    c => c.len_utf8(),
                })
                .sum()
        }

        fn add_text(total: &mut usize, what: &str, value: &str) -> Result<(), PackageError> {
            if value.len() > MAX_JSON_STRING_BYTES {
                return Err(PackageError::TransportLimit {
                    what: what.to_string(),
                    limit: MAX_JSON_STRING_BYTES,
                });
            }
            *total = total
                .checked_add(escaped_len(value))
                .and_then(|sum| sum.checked_add(2))
                .ok_or_else(|| PackageError::TransportLimit {
                    what: "serialized package size".to_string(),
                    limit: MAX_PACKAGE_BYTES,
                })?;
            Ok(())
        }

        check_count("claims", self.claims.len())?;
        let mut bytes = 512usize;
        let mut nodes = 12usize;
        add_text(
            &mut bytes,
            "provenance.code_version",
            &self.provenance.code_version,
        )?;
        add_text(
            &mut bytes,
            "provenance.constellation_lock",
            &self.provenance.constellation_lock,
        )?;
        if let Some(signature) = &self.signature {
            add_text(&mut bytes, "signature", signature)?;
        }
        for (index, claim) in self.claims.iter().enumerate() {
            bytes = bytes
                .checked_add(256)
                .ok_or_else(|| PackageError::TransportLimit {
                    what: "serialized package size".to_string(),
                    limit: MAX_PACKAGE_BYTES,
                })?;
            nodes = nodes.saturating_add(10);
            add_text(&mut bytes, &format!("claims[{index}].id"), &claim.id)?;
            add_text(
                &mut bytes,
                &format!("claims[{index}].statement"),
                &claim.statement,
            )?;
            match &claim.color {
                Color::Verified { .. } => {}
                Color::Validated { regime, dataset } => {
                    check_count("validated regime axes", regime.bounds().len())?;
                    nodes = nodes.saturating_add(regime.bounds().len() * 3);
                    add_text(&mut bytes, &format!("claims[{index}].dataset"), dataset)?;
                    for axis in regime.bounds().keys() {
                        add_text(&mut bytes, &format!("claims[{index}].regime axis"), axis)?;
                        bytes = bytes.saturating_add(64);
                    }
                }
                Color::Estimated { estimator, .. } => {
                    add_text(&mut bytes, &format!("claims[{index}].estimator"), estimator)?
                }
            }
            if let Some(receipt) = &claim.receipt {
                check_count("receipt parents", receipt.parents.len())?;
                nodes = nodes.saturating_add(receipt.parents.len() + 3);
                bytes = bytes.saturating_add(32 * receipt.parents.len() + 64);
            }
            check_count("falsifiers", claim.falsifiers.len())?;
            check_count("anchors", claim.anchors.len())?;
            nodes = nodes.saturating_add(claim.falsifiers.len() * 6 + claim.anchors.len() * 4);
            bytes = bytes.saturating_add(claim.falsifiers.len() * 160 + claim.anchors.len() * 128);
            for falsifier in &claim.falsifiers {
                add_text(&mut bytes, "falsifier.name", &falsifier.name)?;
                add_text(&mut bytes, "falsifier.detail", &falsifier.detail)?;
            }
            for anchor in &claim.anchors {
                add_text(&mut bytes, "anchor.dataset_id", &anchor.dataset_id)?;
                add_text(&mut bytes, "anchor.content_hash", &anchor.content_hash)?;
            }
            if bytes > MAX_PACKAGE_BYTES {
                return Err(PackageError::TransportLimit {
                    what: "serialized package size".to_string(),
                    limit: MAX_PACKAGE_BYTES,
                });
            }
            if nodes > MAX_JSON_NODES {
                return Err(PackageError::TransportLimit {
                    what: "serialized JSON nodes".to_string(),
                    limit: MAX_JSON_NODES,
                });
            }
        }
        Ok(())
    }

    fn verify_finite_magnitude_sums(&self) -> Result<(), PackageError> {
        let mut verified_width = 0.0f64;
        let mut estimated_finite = 0.0f64;
        for claim in &self.claims {
            match &claim.color {
                Color::Verified { lo, hi } => {
                    let width = hi - lo;
                    let next = verified_width + width;
                    if !width.is_finite() || !next.is_finite() {
                        return Err(PackageError::MagnitudeOverflow {
                            claim: claim.id.clone(),
                            component: "verified_width",
                        });
                    }
                    verified_width = next;
                }
                Color::Estimated { dispersion, .. } if dispersion.is_finite() => {
                    let next = estimated_finite + dispersion;
                    if !next.is_finite() {
                        return Err(PackageError::MagnitudeOverflow {
                            claim: claim.id.clone(),
                            component: "estimated_dispersion",
                        });
                    }
                    estimated_finite = next;
                }
                Color::Estimated { .. } | Color::Validated { .. } => {}
            }
        }
        Ok(())
    }

    /// The per-claim uncertainty MAGNITUDE attribution (bead qmao.6.1):
    /// the budget pie over error magnitudes, not claim counts. Verified
    /// claims contribute their interval width, estimated claims their
    /// dispersion; validated claims carry regional trust with no
    /// numeric bound and are reported as an unquantified COUNT rather
    /// than laundered into a number.
    #[must_use]
    pub fn magnitude_budget(&self) -> MagnitudeBudget {
        let mut b = MagnitudeBudget::default();
        for c in &self.claims {
            match &c.color {
                Color::Verified { lo, hi } => b.verified_width += hi - lo,
                Color::Validated { .. } => b.validated_unquantified += 1,
                Color::Estimated { dispersion, .. } => b.estimated_dispersion += dispersion,
            }
        }
        b.quantified_total = b.verified_width + b.estimated_dispersion;
        b
    }

    /// Emit the package as deterministic, self-describing JSON —
    /// schema v2: COMPLETE color payloads (floats as bit-exact hex),
    /// provenance, signature, the content root, and the magnitude
    /// budget. [`EvidencePackage::from_json`] round-trips this
    /// semantically and refuses anything else.
    #[must_use]
    pub fn to_json(&self) -> String {
        use core::fmt::Write as _;
        let mut out = String::new();
        let _ = write!(
            out,
            "{{\"format_version\":{},\"merkle_root\":\"{:016x}\",\"provenance\":{{\"code_version\":\"{}\",\"constellation_lock\":\"{}\"}},\"signature\":",
            self.format_version,
            self.merkle_root(),
            json_escape(&self.provenance.code_version),
            json_escape(&self.provenance.constellation_lock),
        );
        match &self.signature {
            Some(s) => {
                let _ = write!(out, "\"{}\"", json_escape(s));
            }
            None => out.push_str("null"),
        }
        let mb = self.magnitude_budget();
        let _ = write!(
            out,
            ",\"magnitude_budget\":{{\"verified_width_bits\":\"{:016x}\",\"estimated_dispersion_bits\":\"{:016x}\",\"validated_unquantified\":{}}}",
            mb.verified_width.to_bits(),
            mb.estimated_dispersion.to_bits(),
            mb.validated_unquantified
        );
        out.push_str(",\"claims\":[");
        for (i, c) in self.claims.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            push_claim_json(&mut out, c);
        }
        out.push_str("]}");
        out
    }
}

fn push_claim_json(out: &mut String, claim: &Claim) {
    use core::fmt::Write as _;
    let _ = write!(
        out,
        "{{\"id\":\"{}\",\"statement\":\"{}\",\"color\":",
        json_escape(&claim.id),
        json_escape(&claim.statement),
    );
    match &claim.color {
        Color::Verified { lo, hi } => {
            let _ = write!(
                out,
                "{{\"kind\":\"verified\",\"lo_bits\":\"{:016x}\",\"hi_bits\":\"{:016x}\"}}",
                lo.to_bits(),
                hi.to_bits()
            );
        }
        Color::Validated { regime, dataset } => {
            let _ = write!(out, "{{\"kind\":\"validated\",\"regime\":{{");
            for (index, (axis, (lo, hi))) in regime.bounds().iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                let _ = write!(
                    out,
                    "\"{}\":[\"{:016x}\",\"{:016x}\"]",
                    json_escape(axis),
                    lo.to_bits(),
                    hi.to_bits()
                );
            }
            let _ = write!(out, "}},\"dataset\":\"{}\"}}", json_escape(dataset));
        }
        Color::Estimated {
            estimator,
            dispersion,
        } => {
            let _ = write!(
                out,
                "{{\"kind\":\"estimated\",\"estimator\":\"{}\",\"dispersion_bits\":\"{:016x}\"}}",
                json_escape(estimator),
                dispersion.to_bits()
            );
        }
    }
    match &claim.receipt {
        Some(receipt) => {
            let _ = write!(
                out,
                ",\"receipt\":{{\"op\":\"{}\",\"parents\":{:?}}}",
                op_name(receipt.op),
                receipt.parents
            );
        }
        None => out.push_str(",\"receipt\":null"),
    }
    out.push_str(",\"falsifiers\":[");
    for (index, falsifier) in claim.falsifiers.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"name\":\"{}\",\"attempts\":{},\"refuted\":{},\"detail\":\"{}\"}}",
            json_escape(&falsifier.name),
            falsifier.attempts,
            falsifier.refuted,
            json_escape(&falsifier.detail)
        );
    }
    out.push_str("],\"anchors\":[");
    for (index, anchor) in claim.anchors.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"dataset_id\":\"{}\",\"content_hash\":\"{}\"}}",
            json_escape(&anchor.dataset_id),
            json_escape(&anchor.content_hash)
        );
    }
    out.push_str("]}");
}

/// The magnitude budget (see [`EvidencePackage::magnitude_budget`]).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MagnitudeBudget {
    /// Σ (hi − lo) over verified claims.
    pub verified_width: f64,
    /// Σ dispersion over estimated claims.
    pub estimated_dispersion: f64,
    /// Validated claims (regional trust, no numeric bound — counted,
    /// never converted into a fake magnitude).
    pub validated_unquantified: usize,
    /// verified_width + estimated_dispersion (reconciles with the
    /// parts by construction; the parser re-derives and refuses drift).
    pub quantified_total: f64,
}

/// FNV-1a 64-bit hash (pure Rust, std only).
fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h
}

/// Combine two child hashes into a parent hash.
fn combine(a: u64, b: u64) -> u64 {
    let bytes = a.to_le_bytes().into_iter().chain(b.to_le_bytes());
    let mut buf = [0u8; 16];
    for (slot, byte) in buf.iter_mut().zip(bytes) {
        *slot = byte;
    }
    fnv1a(&buf)
}

fn push_atom(out: &mut String, value: &str) {
    use core::fmt::Write as _;
    let _ = write!(out, "{}:", value.len());
    out.push_str(value);
    out.push('|');
}

/// Minimal JSON string escaping.
fn json_escape(s: &str) -> String {
    use core::fmt::Write as _;
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Strict schema-v3 parser (bead qmao.6.1): the package is a PROOF
// ARTIFACT, so parsing fails closed — unknown fields, missing fields,
// wrong types, bad hex, non-finite certificates, a magnitude budget
// that does not re-derive, or an embedded root that does not recompute
// from the parsed fields are each a structured refusal.
// ---------------------------------------------------------------------------

/// A structured parse failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    /// What was being parsed.
    pub what: String,
    /// Why it refused.
    pub why: String,
}

impl core::fmt::Display for ParseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "package parse refused at {}: {}", self.what, self.why)
    }
}

impl core::error::Error for ParseError {}

/// Maximum serialized package size accepted by the standalone checker.
pub const MAX_PACKAGE_BYTES: usize = 64 * 1024 * 1024;
/// Maximum JSON nesting depth before schema mapping.
pub const MAX_JSON_DEPTH: usize = 64;
/// Maximum total JSON values in one package.
pub const MAX_JSON_NODES: usize = 1_000_000;
/// Maximum decoded bytes in a JSON string or object key.
pub const MAX_JSON_STRING_BYTES: usize = 1024 * 1024;
/// Maximum members in any one object or array.
pub const MAX_JSON_CONTAINER_ITEMS: usize = 100_000;
/// Numeric fields in this schema are bounded integers; longer tokens are
/// hostile or malformed even before exact conversion.
pub const MAX_JSON_NUMBER_BYTES: usize = 128;

/// Minimal JSON value for the strict mapper.
#[derive(Debug, Clone, PartialEq)]
enum Jv {
    Null,
    Bool(bool),
    Str(String),
    /// Raw decimal spelling. Integer-valued schema fields must never pass
    /// through `f64`, which cannot represent every `u64` exactly.
    Num(String),
    Arr(Vec<Jv>),
    Obj(Vec<(String, Jv)>),
}

impl Jv {
    fn kind(&self) -> &'static str {
        match self {
            Self::Null => "null",
            Self::Bool(_) => "boolean",
            Self::Str(_) => "string",
            Self::Num(_) => "number",
            Self::Arr(_) => "array",
            Self::Obj(_) => "object",
        }
    }
}

struct Jp<'a> {
    b: &'a [u8],
    at: usize,
    nodes: usize,
}

impl Jp<'_> {
    fn err(&self, what: &str, why: impl Into<String>) -> ParseError {
        ParseError {
            what: format!("{what} (byte {})", self.at),
            why: why.into(),
        }
    }

    fn ws(&mut self) {
        while self
            .b
            .get(self.at)
            .is_some_and(|c| matches!(c, b' ' | b'\t' | b'\n' | b'\r'))
        {
            self.at += 1;
        }
    }

    fn eat(&mut self, c: u8, what: &str) -> Result<(), ParseError> {
        self.ws();
        if self.b.get(self.at) == Some(&c) {
            self.at += 1;
            Ok(())
        } else {
            Err(self.err(what, format!("expected {:?}", char::from(c))))
        }
    }

    fn value(&mut self) -> Result<Jv, ParseError> {
        self.value_at(0)
    }

    fn value_at(&mut self, depth: usize) -> Result<Jv, ParseError> {
        if depth > MAX_JSON_DEPTH {
            return Err(self.err(
                "value",
                format!("nesting depth exceeds limit {MAX_JSON_DEPTH}"),
            ));
        }
        self.nodes = self.nodes.checked_add(1).ok_or_else(|| {
            self.err(
                "value",
                "JSON node counter overflowed before schema mapping",
            )
        })?;
        if self.nodes > MAX_JSON_NODES {
            return Err(self.err(
                "value",
                format!("JSON node count exceeds limit {MAX_JSON_NODES}"),
            ));
        }
        self.ws();
        match self.b.get(self.at) {
            Some(b'"') => Ok(Jv::Str(self.string()?)),
            Some(b'{') => {
                self.at += 1;
                let mut fields = Vec::new();
                let mut keys = std::collections::BTreeSet::new();
                self.ws();
                if self.b.get(self.at) == Some(&b'}') {
                    self.at += 1;
                    return Ok(Jv::Obj(fields));
                }
                loop {
                    let key = self.string()?;
                    self.eat(b':', "object")?;
                    let v = self.value_at(depth + 1)?;
                    if !keys.insert(key.clone()) {
                        return Err(self.err("object", format!("duplicate key {key:?}")));
                    }
                    fields.push((key, v));
                    if fields.len() > MAX_JSON_CONTAINER_ITEMS {
                        return Err(self.err(
                            "object",
                            format!("object member count exceeds limit {MAX_JSON_CONTAINER_ITEMS}"),
                        ));
                    }
                    self.ws();
                    match self.b.get(self.at) {
                        Some(b',') => {
                            self.at += 1;
                            self.ws();
                        }
                        Some(b'}') => {
                            self.at += 1;
                            return Ok(Jv::Obj(fields));
                        }
                        _ => return Err(self.err("object", "expected ',' or '}'")),
                    }
                }
            }
            Some(b'[') => {
                self.at += 1;
                let mut items = Vec::new();
                self.ws();
                if self.b.get(self.at) == Some(&b']') {
                    self.at += 1;
                    return Ok(Jv::Arr(items));
                }
                loop {
                    items.push(self.value_at(depth + 1)?);
                    if items.len() > MAX_JSON_CONTAINER_ITEMS {
                        return Err(self.err(
                            "array",
                            format!("array element count exceeds limit {MAX_JSON_CONTAINER_ITEMS}"),
                        ));
                    }
                    self.ws();
                    match self.b.get(self.at) {
                        Some(b',') => {
                            self.at += 1;
                        }
                        Some(b']') => {
                            self.at += 1;
                            return Ok(Jv::Arr(items));
                        }
                        _ => return Err(self.err("array", "expected ',' or ']'")),
                    }
                }
            }
            Some(b'n') => {
                if self.b[self.at..].starts_with(b"null") {
                    self.at += 4;
                    Ok(Jv::Null)
                } else {
                    Err(self.err("literal", "unknown literal"))
                }
            }
            Some(b't') => {
                if self.b[self.at..].starts_with(b"true") {
                    self.at += 4;
                    Ok(Jv::Bool(true))
                } else {
                    Err(self.err("literal", "unknown literal"))
                }
            }
            Some(b'f') => {
                if self.b[self.at..].starts_with(b"false") {
                    self.at += 5;
                    Ok(Jv::Bool(false))
                } else {
                    Err(self.err("literal", "unknown literal"))
                }
            }
            Some(c) if c.is_ascii_digit() || *c == b'-' => {
                let start = self.at;
                while self.b.get(self.at).is_some_and(|c| {
                    c.is_ascii_digit() || matches!(c, b'-' | b'+' | b'.' | b'e' | b'E')
                }) {
                    self.at += 1;
                }
                if self.at - start > MAX_JSON_NUMBER_BYTES {
                    return Err(self.err(
                        "number",
                        format!("number token exceeds {MAX_JSON_NUMBER_BYTES} bytes"),
                    ));
                }
                let text = core::str::from_utf8(&self.b[start..self.at]).unwrap_or("");
                text.parse::<f64>()
                    .map(|_| Jv::Num(text.to_string()))
                    .map_err(|_| self.err("number", format!("bad number {text:?}")))
            }
            _ => Err(self.err("value", "unexpected byte or end of input")),
        }
    }

    fn string(&mut self) -> Result<String, ParseError> {
        self.ws();
        if self.b.get(self.at) != Some(&b'"') {
            return Err(self.err("string", "expected '\"'"));
        }
        self.at += 1;
        let mut out = String::new();
        loop {
            match self.b.get(self.at) {
                None => return Err(self.err("string", "unterminated")),
                Some(b'"') => {
                    self.at += 1;
                    return Ok(out);
                }
                Some(b'\\') => {
                    self.at += 1;
                    match self.b.get(self.at) {
                        Some(b'"') => out.push('"'),
                        Some(b'\\') => out.push('\\'),
                        Some(b'n') => out.push('\n'),
                        Some(b'r') => out.push('\r'),
                        Some(b't') => out.push('\t'),
                        Some(b'u') => {
                            let hex = self
                                .b
                                .get(self.at + 1..self.at + 5)
                                .and_then(|h| core::str::from_utf8(h).ok())
                                .and_then(|h| u32::from_str_radix(h, 16).ok())
                                .and_then(char::from_u32)
                                .ok_or_else(|| self.err("string", "bad \\u escape"))?;
                            out.push(hex);
                            self.at += 4;
                        }
                        _ => return Err(self.err("string", "bad escape")),
                    }
                    if out.len() > MAX_JSON_STRING_BYTES {
                        return Err(self.err(
                            "string",
                            format!("decoded string exceeds {MAX_JSON_STRING_BYTES} bytes"),
                        ));
                    }
                    self.at += 1;
                }
                Some(&c) if c < 0x20 => {
                    return Err(self.err("string", "unescaped control character"));
                }
                Some(&c) => {
                    // Multi-byte UTF-8 passes through byte-wise.
                    let len = if c < 0x80 {
                        1
                    } else if c >> 5 == 0b110 {
                        2
                    } else if c >> 4 == 0b1110 {
                        3
                    } else {
                        4
                    };
                    let chunk = self
                        .b
                        .get(self.at..self.at + len)
                        .and_then(|ch| core::str::from_utf8(ch).ok())
                        .ok_or_else(|| self.err("string", "invalid UTF-8"))?;
                    out.push_str(chunk);
                    if out.len() > MAX_JSON_STRING_BYTES {
                        return Err(self.err(
                            "string",
                            format!("decoded string exceeds {MAX_JSON_STRING_BYTES} bytes"),
                        ));
                    }
                    self.at += len;
                }
            }
        }
    }
}

fn obj_fields(v: Jv, what: &str) -> Result<Vec<(String, Jv)>, ParseError> {
    match v {
        Jv::Obj(f) => Ok(f),
        other => Err(ParseError {
            what: what.to_string(),
            why: format!("expected an object, got {}", other.kind()),
        }),
    }
}

/// Take field `key` from `fields`; strict mappers call this for every
/// expected key and then refuse leftovers.
fn take_field(fields: &mut Vec<(String, Jv)>, key: &str, what: &str) -> Result<Jv, ParseError> {
    let idx = fields
        .iter()
        .position(|(k, _)| k == key)
        .ok_or(ParseError {
            what: what.to_string(),
            why: format!("missing required field {key:?}"),
        })?;
    Ok(fields.remove(idx).1)
}

fn no_leftovers(fields: &[(String, Jv)], what: &str) -> Result<(), ParseError> {
    if let Some((k, _)) = fields.first() {
        return Err(ParseError {
            what: what.to_string(),
            why: format!("unknown field {k:?} (schema v3 is closed — fail closed)"),
        });
    }
    Ok(())
}

fn as_str(v: Jv, what: &str) -> Result<String, ParseError> {
    match v {
        Jv::Str(s) => Ok(s),
        other => Err(ParseError {
            what: what.to_string(),
            why: format!("expected a string, got {}", other.kind()),
        }),
    }
}

fn hex_u64(v: Jv, what: &str) -> Result<u64, ParseError> {
    let hex = as_str(v, what)?;
    if hex.len() != 16 || !hex.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err(ParseError {
            what: what.to_string(),
            why: format!("expected 16 hex digits, got {hex:?}"),
        });
    }
    Ok(u64::from_str_radix(&hex, 16).expect("validated hexadecimal u64"))
}

fn bits_f64(v: Jv, what: &str, must_be_finite: bool) -> Result<f64, ParseError> {
    let value = f64::from_bits(hex_u64(v, what)?);
    if must_be_finite && !value.is_finite() {
        return Err(ParseError {
            what: what.to_string(),
            why: format!("non-finite value {value} where a finite certificate is required"),
        });
    }
    Ok(value)
}

fn decimal_u64(v: Jv, what: &str) -> Result<u64, ParseError> {
    match v {
        Jv::Num(text) if !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit()) => {
            text.parse::<u64>().map_err(|_| ParseError {
                what: what.to_string(),
                why: format!("unsigned integer out of range: {text:?}"),
            })
        }
        other => Err(ParseError {
            what: what.to_string(),
            why: format!("expected an unsigned decimal integer, got {}", other.kind()),
        }),
    }
}

fn decimal_usize(v: Jv, what: &str) -> Result<usize, ParseError> {
    let value = decimal_u64(v, what)?;
    usize::try_from(value).map_err(|_| ParseError {
        what: what.to_string(),
        why: format!("unsigned integer {value} does not fit usize"),
    })
}

fn parse_package_fields(text: &str) -> Result<Vec<(String, Jv)>, ParseError> {
    if text.len() > MAX_PACKAGE_BYTES {
        return Err(ParseError {
            what: "package".to_string(),
            why: format!("input exceeds the {MAX_PACKAGE_BYTES}-byte package limit"),
        });
    }
    let mut parser = Jp {
        b: text.as_bytes(),
        at: 0,
        nodes: 0,
    };
    let root = parser.value()?;
    parser.ws();
    if parser.at != parser.b.len() {
        return Err(ParseError {
            what: "package".to_string(),
            why: "trailing bytes after the package object".to_string(),
        });
    }
    obj_fields(root, "package")
}

fn parse_format_version(fields: &mut Vec<(String, Jv)>) -> Result<u32, ParseError> {
    let raw = decimal_u64(
        take_field(fields, "format_version", "package")?,
        "format_version",
    )?;
    let version = u32::try_from(raw).map_err(|_| ParseError {
        what: "format_version".to_string(),
        why: format!("version {raw} does not fit u32"),
    })?;
    if version != FORMAT_VERSION {
        return Err(ParseError {
            what: "format_version".to_string(),
            why: format!("unsupported version {version} (this build reads {FORMAT_VERSION})"),
        });
    }
    Ok(version)
}

fn parse_provenance(fields: &mut Vec<(String, Jv)>) -> Result<Provenance, ParseError> {
    let mut provenance = obj_fields(take_field(fields, "provenance", "package")?, "provenance")?;
    let parsed = Provenance {
        code_version: as_str(
            take_field(&mut provenance, "code_version", "provenance")?,
            "code_version",
        )?,
        constellation_lock: as_str(
            take_field(&mut provenance, "constellation_lock", "provenance")?,
            "constellation_lock",
        )?,
    };
    no_leftovers(&provenance, "provenance")?;
    Ok(parsed)
}

fn parse_signature(fields: &mut Vec<(String, Jv)>) -> Result<Option<String>, ParseError> {
    match take_field(fields, "signature", "package")? {
        Jv::Null => Ok(None),
        Jv::Str(signature) => Ok(Some(signature)),
        other => Err(ParseError {
            what: "signature".to_string(),
            why: format!("expected a string or null, got {}", other.kind()),
        }),
    }
}

fn parse_magnitude_budget(fields: &mut Vec<(String, Jv)>) -> Result<MagnitudeBudget, ParseError> {
    let mut budget = obj_fields(
        take_field(fields, "magnitude_budget", "package")?,
        "magnitude_budget",
    )?;
    let verified_width = bits_f64(
        take_field(&mut budget, "verified_width_bits", "magnitude_budget")?,
        "verified_width_bits",
        false,
    )?;
    let estimated_dispersion = bits_f64(
        take_field(&mut budget, "estimated_dispersion_bits", "magnitude_budget")?,
        "estimated_dispersion_bits",
        false,
    )?;
    let validated_unquantified = decimal_usize(
        take_field(&mut budget, "validated_unquantified", "magnitude_budget")?,
        "validated_unquantified",
    )?;
    no_leftovers(&budget, "magnitude_budget")?;
    Ok(MagnitudeBudget {
        verified_width,
        estimated_dispersion,
        validated_unquantified,
        quantified_total: verified_width + estimated_dispersion,
    })
}

fn parse_claims(fields: &mut Vec<(String, Jv)>) -> Result<Vec<Claim>, ParseError> {
    let values = match take_field(fields, "claims", "package")? {
        Jv::Arr(items) => items,
        other => {
            return Err(ParseError {
                what: "claims".to_string(),
                why: format!("expected an array, got {}", other.kind()),
            });
        }
    };
    values
        .into_iter()
        .enumerate()
        .map(|(index, value)| parse_claim(value, index))
        .collect()
}

fn verify_declarations(
    package: &EvidencePackage,
    declared_root: u64,
    declared_budget: MagnitudeBudget,
) -> Result<(), ParseError> {
    let recomputed_budget = package.magnitude_budget();
    if recomputed_budget.verified_width.to_bits() != declared_budget.verified_width.to_bits()
        || recomputed_budget.estimated_dispersion.to_bits()
            != declared_budget.estimated_dispersion.to_bits()
        || recomputed_budget.validated_unquantified != declared_budget.validated_unquantified
    {
        return Err(ParseError {
            what: "magnitude_budget".to_string(),
            why: "declared budget does not re-derive from the claims (tamper or drift)".to_string(),
        });
    }
    let recomputed_root = package.merkle_root();
    if recomputed_root != declared_root {
        return Err(ParseError {
            what: "merkle_root".to_string(),
            why: format!(
                "embedded root {declared_root:016x} does not recompute from the parsed fields \
                 (got {recomputed_root:016x}) — tampered or forged content"
            ),
        });
    }
    Ok(())
}

impl EvidencePackage {
    /// Parse schema-v3 JSON STRICTLY and semantically: every field
    /// mapped, unknown fields refused, floats reconstructed bit-exactly,
    /// the magnitude budget re-derived and compared, and the embedded
    /// content root recomputed from the parsed fields — a package whose
    /// root does not recompute is tampered or forged, and never loads.
    ///
    /// # Errors
    /// [`ParseError`] naming the field and the refusal.
    pub fn from_json(text: &str) -> Result<EvidencePackage, ParseError> {
        let mut fields = parse_package_fields(text)?;
        let format_version = parse_format_version(&mut fields)?;
        let declared_root = hex_u64(
            take_field(&mut fields, "merkle_root", "package")?,
            "merkle_root",
        )?;
        let provenance = parse_provenance(&mut fields)?;
        let signature = parse_signature(&mut fields)?;
        let declared_budget = parse_magnitude_budget(&mut fields)?;
        let claims = parse_claims(&mut fields)?;
        no_leftovers(&fields, "package")?;
        let pkg = EvidencePackage {
            format_version,
            claims,
            provenance,
            signature,
        };
        verify_declarations(&pkg, declared_root, declared_budget)?;
        pkg.verify().map_err(|error| ParseError {
            what: "package semantics".to_string(),
            why: format!("{error:?}"),
        })?;
        Ok(pkg)
    }
}

fn parse_claim(v: Jv, index: usize) -> Result<Claim, ParseError> {
    let what = format!("claims[{index}]");
    let mut f = obj_fields(v, &what)?;
    let id = as_str(take_field(&mut f, "id", &what)?, &what)?;
    let statement = as_str(take_field(&mut f, "statement", &what)?, &what)?;
    let color = parse_color(take_field(&mut f, "color", &what)?, &what)?;
    let receipt_v = take_field(&mut f, "receipt", &what)?;
    let falsifiers_v = take_field(&mut f, "falsifiers", &what)?;
    let anchors_v = take_field(&mut f, "anchors", &what)?;
    no_leftovers(&f, &what)?;
    Ok(Claim {
        id,
        statement,
        color,
        receipt: parse_receipt(receipt_v, &what)?,
        falsifiers: parse_falsifiers(falsifiers_v, &what)?,
        anchors: parse_anchors(anchors_v, &what)?,
    })
}

fn parse_color(value: Jv, what: &str) -> Result<Color, ParseError> {
    let mut fields = obj_fields(value, what)?;
    let kind = as_str(take_field(&mut fields, "kind", what)?, what)?;
    let color = match kind.as_str() {
        "verified" => {
            let lo = bits_f64(take_field(&mut fields, "lo_bits", what)?, what, true)?;
            let hi = bits_f64(take_field(&mut fields, "hi_bits", what)?, what, true)?;
            if lo > hi {
                return Err(ParseError {
                    what: what.to_string(),
                    why: format!("verified interval inverted: {lo} > {hi}"),
                });
            }
            Color::Verified { lo, hi }
        }
        "validated" => {
            let regime_fields = obj_fields(take_field(&mut fields, "regime", what)?, what)?;
            let mut domain = fs_evidence::ValidityDomain::unconstrained();
            for (param, bounds) in regime_fields {
                let Jv::Arr(pair) = bounds else {
                    return Err(ParseError {
                        what: what.to_string(),
                        why: format!("regime {param:?} must be a [lo_bits, hi_bits] pair"),
                    });
                };
                let [lo_v, hi_v]: [Jv; 2] = pair.try_into().map_err(|_| ParseError {
                    what: what.to_string(),
                    why: format!("regime {param:?} must have exactly two bounds"),
                })?;
                let lo = bits_f64(lo_v, what, true)?;
                let hi = bits_f64(hi_v, what, true)?;
                if param.trim().is_empty() {
                    return Err(ParseError {
                        what: what.to_string(),
                        why: "regime axis name must be non-blank".to_string(),
                    });
                }
                if lo > hi {
                    return Err(ParseError {
                        what: what.to_string(),
                        why: format!("regime axis {param:?} has inverted bounds: {lo} > {hi}"),
                    });
                }
                domain = domain.with(param, lo, hi);
            }
            let dataset = as_str(take_field(&mut fields, "dataset", what)?, what)?;
            Color::Validated {
                regime: domain,
                dataset,
            }
        }
        "estimated" => {
            let estimator = as_str(take_field(&mut fields, "estimator", what)?, what)?;
            let dispersion = bits_f64(
                take_field(&mut fields, "dispersion_bits", what)?,
                what,
                false,
            )?;
            if dispersion.is_nan() || dispersion < 0.0 {
                return Err(ParseError {
                    what: what.to_string(),
                    why: format!("NaN or negative dispersion {dispersion}"),
                });
            }
            Color::Estimated {
                estimator,
                dispersion,
            }
        }
        other => {
            return Err(ParseError {
                what: what.to_string(),
                why: format!("unknown color kind {other:?} — fail closed"),
            });
        }
    };
    no_leftovers(&fields, "claim color")?;
    Ok(color)
}

fn parse_receipt(value: Jv, what: &str) -> Result<Option<CompositionReceipt>, ParseError> {
    let Jv::Obj(mut fields) = value else {
        return match value {
            Jv::Null => Ok(None),
            other => Err(ParseError {
                what: what.to_string(),
                why: format!("receipt must be an object or null, got {}", other.kind()),
            }),
        };
    };
    let op_name = as_str(take_field(&mut fields, "op", what)?, what)?;
    let op = op_parse(&op_name).ok_or_else(|| ParseError {
        what: what.to_string(),
        why: format!("unknown receipt op {op_name:?} — fail closed"),
    })?;
    let parents = match take_field(&mut fields, "parents", what)? {
        Jv::Arr(items) => items
            .into_iter()
            .map(|value| decimal_usize(value, what))
            .collect::<Result<Vec<usize>, ParseError>>()?,
        other => {
            return Err(ParseError {
                what: what.to_string(),
                why: format!("receipt parents must be an array, got {}", other.kind()),
            });
        }
    };
    no_leftovers(&fields, "claim receipt")?;
    Ok(Some(CompositionReceipt { parents, op }))
}

fn parse_falsifiers(value: Jv, what: &str) -> Result<Vec<FalsifierRecord>, ParseError> {
    let Jv::Arr(items) = value else {
        return Err(ParseError {
            what: what.to_string(),
            why: "falsifiers must be an array".to_string(),
        });
    };
    items
        .into_iter()
        .map(|value| {
            let mut fields = obj_fields(value, what)?;
            let name = as_str(take_field(&mut fields, "name", what)?, what)?;
            let attempts = decimal_u64(take_field(&mut fields, "attempts", what)?, what)?;
            let refuted = match take_field(&mut fields, "refuted", what)? {
                Jv::Bool(value) => value,
                other => {
                    return Err(ParseError {
                        what: what.to_string(),
                        why: format!("falsifier refuted must be a bool, got {}", other.kind()),
                    });
                }
            };
            let detail = as_str(take_field(&mut fields, "detail", what)?, what)?;
            no_leftovers(&fields, "falsifier record")?;
            Ok(FalsifierRecord {
                name,
                attempts,
                refuted,
                detail,
            })
        })
        .collect()
}

fn parse_anchors(value: Jv, what: &str) -> Result<Vec<AnchorRecord>, ParseError> {
    let Jv::Arr(items) = value else {
        return Err(ParseError {
            what: what.to_string(),
            why: "anchors must be an array".to_string(),
        });
    };
    items
        .into_iter()
        .map(|value| {
            let mut fields = obj_fields(value, what)?;
            let record = AnchorRecord {
                dataset_id: as_str(take_field(&mut fields, "dataset_id", what)?, what)?,
                content_hash: as_str(take_field(&mut fields, "content_hash", what)?, what)?,
            };
            no_leftovers(&fields, "anchor record")?;
            Ok(record)
        })
        .collect()
}
