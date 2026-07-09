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

use fs_evidence::{Color, ColorRank};

/// One claim in an evidence package: a statement plus its epistemic color
/// (which carries the certificate data — an interval, a regime+dataset, or an
/// estimator+dispersion).
#[derive(Debug, Clone, PartialEq)]
pub struct Claim {
    /// A stable claim id.
    pub id: String,
    /// The human-readable claim.
    pub statement: String,
    /// The epistemic color + its certificate payload.
    pub color: Color,
}

impl Claim {
    /// A claim.
    #[must_use]
    pub fn new(id: impl Into<String>, statement: impl Into<String>, color: Color) -> Claim {
        Claim {
            id: id.into(),
            statement: statement.into(),
            color,
        }
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
        out
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
    /// The declared format version is unsupported.
    UnsupportedFormat {
        /// The version found.
        found: u32,
    },
}

/// The one format version this build understands.
pub const FORMAT_VERSION: u32 = 1;

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
        for c in &self.claims {
            match &c.color {
                Color::Verified { lo, hi } => {
                    if !(lo.is_finite() && hi.is_finite() && lo <= hi) {
                        return Err(PackageError::IncompleteVerifiedClaim {
                            claim: c.id.clone(),
                        });
                    }
                }
                Color::Validated { regime, dataset } => {
                    if regime.bounds().is_empty() {
                        return Err(PackageError::IncompleteValidatedClaim {
                            claim: c.id.clone(),
                            missing: "regime",
                        });
                    }
                    if dataset.trim().is_empty() {
                        return Err(PackageError::IncompleteValidatedClaim {
                            claim: c.id.clone(),
                            missing: "dataset",
                        });
                    }
                }
                Color::Estimated { .. } => {} // estimated needs no certificate
            }
        }
        Ok(PackageReport {
            merkle_root: self.merkle_root(),
            breakdown: self.color_breakdown(),
            claims: self.claims.len(),
        })
    }

    /// Emit the package as deterministic, self-describing JSON.
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
        out.push_str(",\"claims\":[");
        for (i, c) in self.claims.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"id\":\"{}\",\"statement\":\"{}\",\"color\":\"{}\"}}",
                json_escape(&c.id),
                json_escape(&c.statement),
                c.color.rank_str(),
            );
        }
        out.push_str("]}");
        out
    }
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

/// Rank-string helper on [`Color`] for JSON.
trait RankStr {
    fn rank_str(&self) -> &'static str;
}
impl RankStr for Color {
    fn rank_str(&self) -> &'static str {
        match self.rank() {
            ColorRank::Verified => "verified",
            ColorRank::Validated => "validated",
            ColorRank::Estimated => "estimated",
        }
    }
}
