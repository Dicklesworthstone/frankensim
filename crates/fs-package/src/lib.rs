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

/// The one format version this build understands. v2 (bead qmao.6.1):
/// claims carry their COMPLETE color payloads bit-exactly, the JSON
/// round-trips through a strict parser, and the embedded content root
/// must recompute from the parsed fields.
pub const FORMAT_VERSION: u32 = 2;

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
            let _ = write!(
                out,
                "{{\"id\":\"{}\",\"statement\":\"{}\",\"color\":",
                json_escape(&c.id),
                json_escape(&c.statement),
            );
            match &c.color {
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
                    for (j, (k, (lo, hi))) in regime.bounds().iter().enumerate() {
                        if j > 0 {
                            out.push(',');
                        }
                        let _ = write!(
                            out,
                            "\"{}\":[\"{:016x}\",\"{:016x}\"]",
                            json_escape(k),
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
            out.push('}');
        }
        out.push_str("]}");
        out
    }
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

// ---------------------------------------------------------------------------
// Strict schema-v2 parser (bead qmao.6.1): the package is a PROOF
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

/// Minimal JSON value for the strict mapper.
#[derive(Debug, Clone, PartialEq)]
enum Jv {
    Null,
    Str(String),
    Num(f64),
    Arr(Vec<Jv>),
    Obj(Vec<(String, Jv)>),
}

struct Jp<'a> {
    b: &'a [u8],
    at: usize,
}

impl<'a> Jp<'a> {
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
        self.ws();
        match self.b.get(self.at) {
            Some(b'"') => Ok(Jv::Str(self.string()?)),
            Some(b'{') => {
                self.at += 1;
                let mut fields = Vec::new();
                self.ws();
                if self.b.get(self.at) == Some(&b'}') {
                    self.at += 1;
                    return Ok(Jv::Obj(fields));
                }
                loop {
                    let key = self.string()?;
                    self.eat(b':', "object")?;
                    let v = self.value()?;
                    if fields.iter().any(|(k, _)| *k == key) {
                        return Err(self.err("object", format!("duplicate key {key:?}")));
                    }
                    fields.push((key, v));
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
                    items.push(self.value()?);
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
            Some(c) if c.is_ascii_digit() || *c == b'-' => {
                let start = self.at;
                while self
                    .b
                    .get(self.at)
                    .is_some_and(|c| c.is_ascii_digit() || matches!(c, b'-' | b'+' | b'.' | b'e' | b'E'))
                {
                    self.at += 1;
                }
                let text = core::str::from_utf8(&self.b[start..self.at]).unwrap_or("");
                text.parse()
                    .map(Jv::Num)
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
                    self.at += 1;
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
            why: format!("expected an object, got {other:?}"),
        }),
    }
}

/// Take field `key` from `fields`; strict mappers call this for every
/// expected key and then refuse leftovers.
fn take_field(fields: &mut Vec<(String, Jv)>, key: &str, what: &str) -> Result<Jv, ParseError> {
    let idx = fields.iter().position(|(k, _)| k == key).ok_or(ParseError {
        what: what.to_string(),
        why: format!("missing required field {key:?}"),
    })?;
    Ok(fields.remove(idx).1)
}

fn no_leftovers(fields: &[(String, Jv)], what: &str) -> Result<(), ParseError> {
    if let Some((k, _)) = fields.first() {
        return Err(ParseError {
            what: what.to_string(),
            why: format!("unknown field {k:?} (schema v2 is closed — fail closed)"),
        });
    }
    Ok(())
}

fn as_str(v: Jv, what: &str) -> Result<String, ParseError> {
    match v {
        Jv::Str(s) => Ok(s),
        other => Err(ParseError {
            what: what.to_string(),
            why: format!("expected a string, got {other:?}"),
        }),
    }
}

fn bits_f64(v: Jv, what: &str, must_be_finite: bool) -> Result<f64, ParseError> {
    let hex = as_str(v, what)?;
    if hex.len() != 16 || !hex.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err(ParseError {
            what: what.to_string(),
            why: format!("expected 16 hex digits of f64 bits, got {hex:?}"),
        });
    }
    let bits = u64::from_str_radix(&hex, 16).expect("validated hex");
    let value = f64::from_bits(bits);
    if must_be_finite && !value.is_finite() {
        return Err(ParseError {
            what: what.to_string(),
            why: format!("non-finite value {value} where a finite certificate is required"),
        });
    }
    Ok(value)
}

impl EvidencePackage {
    /// Parse schema-v2 JSON STRICTLY and semantically: every field
    /// mapped, unknown fields refused, floats reconstructed bit-exactly,
    /// the magnitude budget re-derived and compared, and the embedded
    /// content root recomputed from the parsed fields — a package whose
    /// root does not recompute is tampered or forged, and never loads.
    ///
    /// # Errors
    /// [`ParseError`] naming the field and the refusal.
    pub fn from_json(text: &str) -> Result<EvidencePackage, ParseError> {
        let mut jp = Jp {
            b: text.as_bytes(),
            at: 0,
        };
        let root_value = jp.value()?;
        jp.ws();
        if jp.at != jp.b.len() {
            return Err(ParseError {
                what: "package".to_string(),
                why: "trailing bytes after the package object".to_string(),
            });
        }
        let mut fields = obj_fields(root_value, "package")?;
        let format_version = match take_field(&mut fields, "format_version", "package")? {
            Jv::Num(n) if n.fract() == 0.0 && (0.0..=f64::from(u32::MAX)).contains(&n) => n as u32,
            other => {
                return Err(ParseError {
                    what: "format_version".to_string(),
                    why: format!("expected a u32, got {other:?}"),
                });
            }
        };
        if format_version != FORMAT_VERSION {
            return Err(ParseError {
                what: "format_version".to_string(),
                why: format!("unsupported version {format_version} (this build reads {FORMAT_VERSION})"),
            });
        }
        let declared_root = {
            let hex = as_str(take_field(&mut fields, "merkle_root", "package")?, "merkle_root")?;
            u64::from_str_radix(&hex, 16).map_err(|_| ParseError {
                what: "merkle_root".to_string(),
                why: format!("expected 16 hex digits, got {hex:?}"),
            })?
        };
        let mut prov = obj_fields(take_field(&mut fields, "provenance", "package")?, "provenance")?;
        let provenance = Provenance {
            code_version: as_str(take_field(&mut prov, "code_version", "provenance")?, "code_version")?,
            constellation_lock: as_str(
                take_field(&mut prov, "constellation_lock", "provenance")?,
                "constellation_lock",
            )?,
        };
        no_leftovers(&prov, "provenance")?;
        let signature = match take_field(&mut fields, "signature", "package")? {
            Jv::Null => None,
            Jv::Str(s) => Some(s),
            other => {
                return Err(ParseError {
                    what: "signature".to_string(),
                    why: format!("expected a string or null, got {other:?}"),
                });
            }
        };
        let mut mb = obj_fields(
            take_field(&mut fields, "magnitude_budget", "package")?,
            "magnitude_budget",
        )?;
        let declared_width = bits_f64(
            take_field(&mut mb, "verified_width_bits", "magnitude_budget")?,
            "verified_width_bits",
            false,
        )?;
        let declared_disp = bits_f64(
            take_field(&mut mb, "estimated_dispersion_bits", "magnitude_budget")?,
            "estimated_dispersion_bits",
            false,
        )?;
        let declared_unq = match take_field(&mut mb, "validated_unquantified", "magnitude_budget")? {
            Jv::Num(n) if n.fract() == 0.0 && n >= 0.0 => n as usize,
            other => {
                return Err(ParseError {
                    what: "validated_unquantified".to_string(),
                    why: format!("expected a count, got {other:?}"),
                });
            }
        };
        no_leftovers(&mb, "magnitude_budget")?;
        let claims_v = match take_field(&mut fields, "claims", "package")? {
            Jv::Arr(items) => items,
            other => {
                return Err(ParseError {
                    what: "claims".to_string(),
                    why: format!("expected an array, got {other:?}"),
                });
            }
        };
        no_leftovers(&fields, "package")?;
        let mut claims = Vec::with_capacity(claims_v.len());
        for (i, cv) in claims_v.into_iter().enumerate() {
            claims.push(parse_claim(cv, i)?);
        }
        let pkg = EvidencePackage {
            format_version,
            claims,
            provenance,
            signature,
        };
        // SEMANTIC re-derivation: the magnitude budget and the content
        // root must recompute from the parsed fields, bit for bit.
        let mb2 = pkg.magnitude_budget();
        if mb2.verified_width.to_bits() != declared_width.to_bits()
            || mb2.estimated_dispersion.to_bits() != declared_disp.to_bits()
            || mb2.validated_unquantified != declared_unq
        {
            return Err(ParseError {
                what: "magnitude_budget".to_string(),
                why: "declared budget does not re-derive from the claims (tamper or drift)"
                    .to_string(),
            });
        }
        let recomputed = pkg.merkle_root();
        if recomputed != declared_root {
            return Err(ParseError {
                what: "merkle_root".to_string(),
                why: format!(
                    "embedded root {declared_root:016x} does not recompute from the parsed \
                     fields (got {recomputed:016x}) — tampered or forged content"
                ),
            });
        }
        Ok(pkg)
    }
}

fn parse_claim(v: Jv, index: usize) -> Result<Claim, ParseError> {
    let what = format!("claims[{index}]");
    let mut f = obj_fields(v, &what)?;
    let id = as_str(take_field(&mut f, "id", &what)?, &what)?;
    let statement = as_str(take_field(&mut f, "statement", &what)?, &what)?;
    let mut cf = obj_fields(take_field(&mut f, "color", &what)?, &what)?;
    no_leftovers(&f, &what)?;
    let kind = as_str(take_field(&mut cf, "kind", &what)?, &what)?;
    let color = match kind.as_str() {
        "verified" => {
            let lo = bits_f64(take_field(&mut cf, "lo_bits", &what)?, &what, true)?;
            let hi = bits_f64(take_field(&mut cf, "hi_bits", &what)?, &what, true)?;
            if lo > hi {
                return Err(ParseError {
                    what,
                    why: format!("verified interval inverted: {lo} > {hi}"),
                });
            }
            Color::Verified { lo, hi }
        }
        "validated" => {
            let regime_fields = obj_fields(take_field(&mut cf, "regime", &what)?, &what)?;
            let mut domain = fs_evidence::ValidityDomain::unconstrained();
            for (param, bounds) in regime_fields {
                let Jv::Arr(pair) = bounds else {
                    return Err(ParseError {
                        what,
                        why: format!("regime {param:?} must be a [lo_bits, hi_bits] pair"),
                    });
                };
                let [lo_v, hi_v]: [Jv; 2] = pair.try_into().map_err(|_| ParseError {
                    what: what.clone(),
                    why: format!("regime {param:?} must have exactly two bounds"),
                })?;
                let lo = bits_f64(lo_v, &what, true)?;
                let hi = bits_f64(hi_v, &what, true)?;
                domain = domain.with(param, lo, hi);
            }
            let dataset = as_str(take_field(&mut cf, "dataset", &what)?, &what)?;
            Color::Validated {
                regime: domain,
                dataset,
            }
        }
        "estimated" => {
            let estimator = as_str(take_field(&mut cf, "estimator", &what)?, &what)?;
            let dispersion = bits_f64(take_field(&mut cf, "dispersion_bits", &what)?, &what, true)?;
            if dispersion < 0.0 {
                return Err(ParseError {
                    what,
                    why: format!("negative dispersion {dispersion}"),
                });
            }
            Color::Estimated {
                estimator,
                dispersion,
            }
        }
        other => {
            return Err(ParseError {
                what,
                why: format!("unknown color kind {other:?} — fail closed"),
            });
        }
    };
    no_leftovers(&cf, "claim color")?;
    Ok(Claim {
        id,
        statement,
        color,
    })
}
