//! Canonical serialization: problems round-trip through the six-base
//! `fsopt v2` line form whose floats are BIT PATTERNS (exact round-trip),
//! and the problem HASH (FNV-1a 64 over the canonical body) is the study
//! identity. The receipt-bearing reader also decodes exact explicit five-base
//! v1 inputs emitted by either known historical v1 string encoding, appends
//! `mol = 0`, and binds the complete old/new artifacts with BLAKE3. Parsing
//! rebuilds THROUGH the validating builder, so a tampered file cannot smuggle
//! in an ill-typed graph — revalidation is free.

use crate::ir::{
    ConstraintKind, Expr, Manifold, NodeId, OptError, Problem, ProblemBuilder, ProblemTag, Sense,
    VarId,
};
use fs_blake3::hash_bytes;
use fs_qty::Dims;
use std::fmt::Write as _;

pub use fs_blake3::ContentHash;

/// Escape a string for single-token embedding. Percent-encodes the token
/// delimiters AND every control / non-ASCII byte, so ANY value (including
/// multibyte UTF-8) round-trips exactly. Graphic ASCII except `%` passes
/// through verbatim, so ASCII fields serialize byte-identically to before.
fn esc(s: &str) -> String {
    use std::fmt::Write;
    let mut out = String::with_capacity(s.len());
    for &b in s.as_bytes() {
        if b.is_ascii_graphic() && b != b'%' {
            out.push(b as char);
        } else {
            let _ = write!(out, "%{b:02X}");
        }
    }
    out
}

/// Original v1 writer encoding from `293b1c1`: only token delimiters were
/// escaped, while all other Unicode scalar values were written literally.
/// This is retained solely to recognize exact artifacts emitted before the
/// byte-percent writer landed in `b44a1a9`; current v2 always uses [`esc`].
fn legacy_literal_esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for character in s.chars() {
        match character {
            '%' => out.push_str("%25"),
            ' ' => out.push_str("%20"),
            '\n' => out.push_str("%0A"),
            _ => out.push(character),
        }
    }
    out
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TokenEncoding {
    PercentBytes,
    LegacyLiteralUtf8,
}

fn escape_for_wire(value: &str, encoding: TokenEncoding) -> String {
    match encoding {
        TokenEncoding::PercentBytes => esc(value),
        TokenEncoding::LegacyLiteralUtf8 => legacy_literal_esc(value),
    }
}

/// Inverse of [`esc`]. Decodes `%XX` on the raw BYTE stream (so a multibyte
/// value reassembles correctly) and rejects malformed escapes or decoded bytes
/// that are not valid UTF-8. A legacy migration receipt must never certify a
/// lossy replacement-character rewrite as an amount-dimension crosswalk.
fn unesc(s: &str) -> Result<String, &'static str> {
    let src = s.as_bytes();
    let mut bytes = Vec::with_capacity(src.len());
    let mut i = 0;
    while i < src.len() {
        if src[i] == b'%' {
            if i + 3 > src.len() {
                return Err("truncated percent escape");
            }
            let (Some(hi), Some(lo)) = (hex_nibble(src[i + 1]), hex_nibble(src[i + 2])) else {
                return Err("percent escape must contain exactly two hexadecimal digits");
            };
            bytes.push((hi << 4) | lo);
            i += 3;
            continue;
        }
        bytes.push(src[i]);
        i += 1;
    }
    String::from_utf8(bytes).map_err(|_| "percent-decoded token is not valid UTF-8")
}

fn hex_nibble(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'a'..=b'f' => Some(b - b'a' + 10),
        b'A'..=b'F' => Some(b - b'A' + 10),
        _ => None,
    }
}

/// Version of the line-oriented `fsopt` representation read from a file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WireVersion {
    /// Legacy five-base dimensions `(m, kg, s, K, A)` with an explicit v1
    /// header.
    V1,
    /// Canonical six-base dimensions `(m, kg, s, K, A, mol)`.
    V2,
}

/// The only admitted semantic rule for converting a v1 dimension vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FiveToSixRule {
    /// Preserve all five exponents and append an exact zero mole exponent.
    AppendMoleZero,
}

/// Immutable evidence that exact legacy bytes were mapped to canonical v2.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DimensionCrosswalkReceipt {
    source_version: WireVersion,
    target_version: WireVersion,
    old_hash: ContentHash,
    new_hash: ContentHash,
    rule: FiveToSixRule,
}

impl DimensionCrosswalkReceipt {
    /// Source schema named by the receipt.
    #[must_use]
    pub const fn source_version(&self) -> WireVersion {
        self.source_version
    }

    /// Target schema named by the receipt.
    #[must_use]
    pub const fn target_version(&self) -> WireVersion {
        self.target_version
    }

    /// BLAKE3 content hash of the exact source artifact, including its hash line.
    #[must_use]
    pub const fn old_hash(&self) -> ContentHash {
        self.old_hash
    }

    /// BLAKE3 content hash of the exact canonical target artifact.
    #[must_use]
    pub const fn new_hash(&self) -> ContentHash {
        self.new_hash
    }

    /// Semantic migration rule applied.
    #[must_use]
    pub const fn rule(&self) -> FiveToSixRule {
        self.rule
    }

    /// Verify this receipt against the exact preserved source and target bytes.
    #[must_use]
    pub fn verifies(&self, old_bytes: &[u8], new_bytes: &[u8]) -> bool {
        self.source_version == WireVersion::V1
            && self.target_version == WireVersion::V2
            && self.rule == FiveToSixRule::AppendMoleZero
            && hash_bytes(old_bytes) == self.old_hash
            && hash_bytes(new_bytes) == self.new_hash
    }
}

/// A parsed problem together with its source provenance and mandatory v1
/// semantic-crosswalk receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedProblem {
    problem: Problem,
    source_version: WireVersion,
    source_hash: u64,
    migration: Option<DimensionCrosswalkReceipt>,
}

impl ParsedProblem {
    /// Revalidated optimization problem. Legacy v1 dimensions have `mol = 0`.
    #[must_use]
    pub const fn problem(&self) -> &Problem {
        &self.problem
    }

    /// Wire version declared by the source header.
    #[must_use]
    pub const fn source_version(&self) -> WireVersion {
        self.source_version
    }

    /// FNV-1a integrity hash embedded in the exact accepted source artifact.
    #[must_use]
    pub const fn source_hash(&self) -> u64 {
        self.source_hash
    }

    /// Five-to-six migration evidence; present exactly for v1 input.
    #[must_use]
    pub const fn migration(&self) -> Option<&DimensionCrosswalkReceipt> {
        self.migration.as_ref()
    }

    /// Consume the parsed result without allowing any provenance component to
    /// disappear implicitly.
    #[must_use]
    pub fn into_parts(self) -> (Problem, WireVersion, u64, Option<DimensionCrosswalkReceipt>) {
        (
            self.problem,
            self.source_version,
            self.source_hash,
            self.migration,
        )
    }
}

fn dims_str(d: Dims, version: WireVersion) -> String {
    match version {
        WireVersion::V1 => format!("({},{},{},{},{})", d.0[0], d.0[1], d.0[2], d.0[3], d.0[4]),
        WireVersion::V2 => format!(
            "({},{},{},{},{},{})",
            d.0[0], d.0[1], d.0[2], d.0[3], d.0[4], d.0[5]
        ),
    }
}

fn parse_dims(s: &str, version: WireVersion) -> Option<Dims> {
    let inner = s.strip_prefix('(')?.strip_suffix(')')?;
    let parts: Vec<&str> = inner.split(',').collect();
    let expected_len = match version {
        WireVersion::V1 => 5,
        WireVersion::V2 => 6,
    };
    if parts.len() != expected_len {
        return None;
    }
    let mut d = [0i8; 6];
    for (slot, p) in d.iter_mut().zip(&parts) {
        *slot = p.parse().ok()?;
    }
    Some(Dims(d))
}

fn f64_hex(v: f64) -> String {
    format!("{:016X}", v.to_bits())
}

fn parse_f64_hex(s: &str) -> Option<f64> {
    u64::from_str_radix(s, 16).ok().map(f64::from_bits)
}

fn manifold_str(m: Manifold) -> String {
    match m {
        Manifold::Rn { dim } => format!("Rn({dim})"),
        Manifold::Sphere { ambient } => format!("Sphere({ambient})"),
        Manifold::So3 => "So3".to_string(),
        Manifold::Stiefel { n, p } => format!("Stiefel({n},{p})"),
    }
}

fn parse_manifold(s: &str) -> Option<Manifold> {
    if s == "So3" {
        return Some(Manifold::So3);
    }
    if let Some(inner) = s.strip_prefix("Rn(").and_then(|r| r.strip_suffix(')')) {
        return Some(Manifold::Rn {
            dim: inner.parse().ok()?,
        });
    }
    if let Some(inner) = s.strip_prefix("Sphere(").and_then(|r| r.strip_suffix(')')) {
        return Some(Manifold::Sphere {
            ambient: inner.parse().ok()?,
        });
    }
    if let Some(inner) = s.strip_prefix("Stiefel(").and_then(|r| r.strip_suffix(')')) {
        let (n, p) = inner.split_once(',')?;
        return Some(Manifold::Stiefel {
            n: n.parse().ok()?,
            p: p.parse().ok()?,
        });
    }
    None
}

/// Canonical single-token form of one expression (the hash-consing key
/// AND the serialized body).
#[must_use]
pub(crate) fn expr_key(e: &Expr) -> String {
    expr_key_for_wire(e, WireVersion::V2, TokenEncoding::PercentBytes)
}

fn expr_key_for_wire(e: &Expr, version: WireVersion, encoding: TokenEncoding) -> String {
    match e {
        Expr::Var(v) => format!("var {}", v.0),
        Expr::Component { of, index } => format!("component {} {index}", of.0),
        Expr::Const { value, dims } => {
            format!("const {} {}", f64_hex(*value), dims_str(*dims, version))
        }
        Expr::Add(a, b) => format!("add {} {}", a.0, b.0),
        Expr::Sub(a, b) => format!("sub {} {}", a.0, b.0),
        Expr::Mul(a, b) => format!("mul {} {}", a.0, b.0),
        Expr::Div(a, b) => format!("div {} {}", a.0, b.0),
        Expr::Neg(a) => format!("neg {}", a.0),
        Expr::Powi { base, exp } => format!("powi {} {exp}", base.0),
        Expr::Sqrt(a) => format!("sqrt {}", a.0),
        Expr::Exp(a) => format!("exp {}", a.0),
        Expr::Ln(a) => format!("ln {}", a.0),
        Expr::Tanh(a) => format!("tanh {}", a.0),
        Expr::Dot(a, b) => format!("dot {} {}", a.0, b.0),
        Expr::NormSq(a) => format!("norm_sq {}", a.0),
        Expr::Min(a, b) => format!("min {} {}", a.0, b.0),
        Expr::Max(a, b) => format!("max {} {}", a.0, b.0),
        Expr::Abs(a) => format!("abs {}", a.0),
        Expr::PdeResidual {
            study,
            over,
            adjoint_available,
            dims,
        } => format!(
            "pde_residual {} {} {} {}",
            escape_for_wire(study, encoding),
            over.0,
            u8::from(*adjoint_available),
            dims_str(*dims, version)
        ),
        Expr::Expectation { of, uq_config } => {
            format!(
                "expectation {} {}",
                of.0,
                escape_for_wire(uq_config, encoding)
            )
        }
        Expr::Cvar {
            of,
            alpha,
            uq_config,
        } => format!(
            "cvar {} {} {}",
            of.0,
            f64_hex(*alpha),
            escape_for_wire(uq_config, encoding)
        ),
        Expr::Quantile { of, q, uq_config } => {
            format!(
                "quantile {} {} {}",
                of.0,
                f64_hex(*q),
                escape_for_wire(uq_config, encoding)
            )
        }
    }
}

/// FNV-1a 64 (in-house; stable across platforms).
fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

fn body_for_wire(problem: &Problem, version: WireVersion, encoding: TokenEncoding) -> String {
    let version_tag = match version {
        WireVersion::V1 => "v1",
        WireVersion::V2 => "v2",
    };
    let mut s = format!("fsopt {version_tag}\n");
    for (i, v) in problem.vars.iter().enumerate() {
        let _ = writeln!(
            s,
            "var {i} {} {} {}",
            escape_for_wire(&v.name, encoding),
            manifold_str(v.manifold),
            dims_str(v.dims, version)
        );
    }
    for (i, e) in problem.exprs.iter().enumerate() {
        let _ = writeln!(s, "expr {i} {}", expr_key_for_wire(e, version, encoding));
    }
    for o in &problem.objectives {
        let sense = match o.sense {
            Sense::Minimize => "min",
            Sense::Maximize => "max",
        };
        let _ = writeln!(s, "objective {sense} {} {}", o.node.0, f64_hex(o.weight));
    }
    for c in &problem.constraints {
        let kind = match c.kind {
            ConstraintKind::EqZero => "eq0",
            ConstraintKind::LeZero => "le0",
        };
        let _ = writeln!(
            s,
            "constraint {kind} {} {}",
            c.node.0,
            escape_for_wire(&c.name, encoding)
        );
    }
    for t in &problem.tags {
        match t {
            ProblemTag::MultiFidelity { levels } => {
                let _ = writeln!(s, "tag multi_fidelity {levels}");
            }
            ProblemTag::ChanceConstrained { prob } => {
                let _ = writeln!(s, "tag chance {}", f64_hex(*prob));
            }
            ProblemTag::Bilevel { inner_hash } => {
                let _ = writeln!(s, "tag bilevel {inner_hash:016X}");
            }
        }
    }
    let _ = writeln!(s, "budget {}", problem.budget.max_evals);
    s
}

fn body(problem: &Problem) -> String {
    body_for_wire(problem, WireVersion::V2, TokenEncoding::PercentBytes)
}

fn serialize_for_wire(problem: &Problem, version: WireVersion, encoding: TokenEncoding) -> String {
    let body = body_for_wire(problem, version, encoding);
    format!("{body}hash {:016X}\n", fnv1a(body.as_bytes()))
}

/// The problem hash: FNV-1a 64 over the canonical body (study
/// identity; two structurally identical builds hash identically).
#[must_use]
pub fn problem_hash(problem: &Problem) -> u64 {
    fnv1a(body(problem).as_bytes())
}

/// Serialize to the canonical text form (hash line appended).
#[must_use]
pub fn serialize(problem: &Problem) -> String {
    serialize_for_wire(problem, WireVersion::V2, TokenEncoding::PercentBytes)
}

fn perr(line: usize, what: impl Into<String>) -> OptError {
    OptError::Parse {
        line,
        what: what.into(),
    }
}

fn require_end<'a, I>(tokens: &mut I, line: usize, directive: &str) -> Result<(), OptError>
where
    I: Iterator<Item = &'a str>,
{
    if tokens.next().is_some() {
        return Err(perr(
            line,
            format!("{directive} directive has trailing fields"),
        ));
    }
    Ok(())
}

fn first_difference_line(actual: &str, canonical: &str) -> usize {
    let mut line = 1usize;
    for (actual_byte, canonical_byte) in actual.bytes().zip(canonical.bytes()) {
        if actual_byte != canonical_byte {
            return line;
        }
        if actual_byte == b'\n' {
            line += 1;
        }
    }
    line
}

fn terminal_hash_line(lines: &[&str]) -> Result<usize, OptError> {
    let mut found = None;
    for (line_index, line) in lines.iter().enumerate() {
        if line.split(' ').next() != Some("hash") {
            continue;
        }
        if found.is_some() {
            return Err(perr(
                line_index + 1,
                "duplicate hash directive; exactly one terminal hash is required",
            ));
        }
        found = Some(line_index);
    }
    let line_index = found.ok_or_else(|| {
        perr(
            lines.len().max(1),
            "missing terminal hash directive; exactly one is required",
        )
    })?;
    if line_index + 1 != lines.len() {
        return Err(perr(
            line_index + 1,
            "hash directive must be the final line of the artifact",
        ));
    }
    Ok(line_index)
}

/// Parse the current canonical v2 form, REBUILDING through the validating
/// builder and verifying its integrity hash.
///
/// Legacy v1 input is rejected here because returning only [`Problem`] would
/// discard mandatory migration evidence. Use [`parse_with_version`] for
/// explicit v1 artifacts.
///
/// # Errors
/// [`OptError::Parse`] (with line numbers) or any builder error.
pub fn parse(text: &str) -> Result<Problem, OptError> {
    let parsed = parse_with_version(text)?;
    if parsed.migration.is_some() {
        return Err(perr(
            1,
            "legacy five-base fsopt requires parse_with_version so its semantic-crosswalk receipt is retained",
        ));
    }
    Ok(parsed.problem)
}

/// Parse the canonical form while retaining its source wire version and
/// embedded hash for provenance/crosswalk recording by the caller.
/// Explicit v1 and v2 are admitted only when the complete source bytes equal a
/// known canonical encoding for that wire version. V1 recognizes both actual
/// historical writers; a headerless artifact has no schema identity and is
/// refused.
///
/// # Errors
/// [`OptError::Parse`] (with line numbers) or any builder error.
#[allow(clippy::too_many_lines)] // one grammar rule per block
pub fn parse_with_version(text: &str) -> Result<ParsedProblem, OptError> {
    let mut b = ProblemBuilder::new();
    let mut expected_hash: Option<u64> = None;
    // A version header is mandatory. Defaulting a headerless artifact to v1
    // would invent provenance for bytes no historical writer emitted.
    let mut source_version = WireVersion::V2;
    let mut explicit_version = false;
    let mut saw_body_directive = false;
    let mut saw_budget = false;
    // Strip only the one canonical terminal LF and split on LF ourselves.
    // `str::lines` would silently remove a preceding CR, but the original v1
    // writer could emit a literal CR inside a final string token. Exact legacy
    // recognition must preserve that byte; ordinary CRLF input still fails
    // the declared-version/canonical-byte checks below.
    let line_text = text.strip_suffix('\n').unwrap_or(text);
    let lines: Vec<&str> = line_text.split('\n').collect();
    let hash_line = terminal_hash_line(&lines)?;
    for (ln0, line) in lines.iter().enumerate() {
        let ln = ln0 + 1;
        let mut tok = line.split(' ');
        let head = tok.next().unwrap_or("");
        if ln == 1 && head != "fsopt" {
            return Err(perr(
                ln,
                "missing leading fsopt version header; headerless artifacts have no schema identity",
            ));
        }
        if !head.is_empty() && head != "fsopt" {
            saw_body_directive = true;
        }
        match head {
            "fsopt" => {
                if explicit_version {
                    return Err(perr(ln, "duplicate fsopt version header"));
                }
                if saw_body_directive {
                    return Err(perr(ln, "fsopt version header must precede the body"));
                }
                source_version = match tok.next() {
                    Some("v1") => WireVersion::V1,
                    Some("v2") => WireVersion::V2,
                    _ => {
                        return Err(perr(
                            ln,
                            "unsupported version (expected `fsopt v1` or `fsopt v2`)",
                        ));
                    }
                };
                require_end(&mut tok, ln, "fsopt version header")?;
                explicit_version = true;
            }
            "var" => {
                let ix: u32 = tok
                    .next()
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "var: missing index"))?;
                let name = unesc(tok.next().ok_or_else(|| perr(ln, "var: missing name"))?)
                    .map_err(|what| perr(ln, format!("var: bad name: {what}")))?;
                let manifold = tok
                    .next()
                    .and_then(parse_manifold)
                    .ok_or_else(|| perr(ln, "var: bad manifold"))?;
                let dims = tok
                    .next()
                    .and_then(|value| parse_dims(value, source_version))
                    .ok_or_else(|| perr(ln, "var: bad dims"))?;
                require_end(&mut tok, ln, "var")?;
                let assigned = b.var(&name, manifold, dims);
                if assigned.0 != ix {
                    return Err(perr(
                        ln,
                        format!(
                            "var id mismatch: file says {ix}, canonical rebuild gives {}",
                            assigned.0
                        ),
                    ));
                }
            }
            "expr" => {
                let ix: u32 = tok
                    .next()
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "expr: missing index"))?;
                let rest: Vec<&str> = tok.collect();
                let id = parse_expr(&mut b, &rest, ln, source_version)?;
                if id.0 != ix {
                    return Err(perr(
                        ln,
                        format!(
                            "expr id mismatch: file says {ix}, canonical rebuild gives {} \
                             (duplicate or reordered nodes break identity)",
                            id.0
                        ),
                    ));
                }
            }
            "objective" => {
                let sense = match tok.next() {
                    Some("min") => Sense::Minimize,
                    Some("max") => Sense::Maximize,
                    _ => return Err(perr(ln, "objective: bad sense")),
                };
                let node: u32 = tok
                    .next()
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "objective: bad node"))?;
                let weight = tok
                    .next()
                    .and_then(parse_f64_hex)
                    .ok_or_else(|| perr(ln, "objective: bad weight"))?;
                require_end(&mut tok, ln, "objective")?;
                b.objective(NodeId(node), sense, weight)?;
            }
            "constraint" => {
                let kind = match tok.next() {
                    Some("eq0") => ConstraintKind::EqZero,
                    Some("le0") => ConstraintKind::LeZero,
                    _ => return Err(perr(ln, "constraint: bad kind")),
                };
                let node: u32 = tok
                    .next()
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "constraint: bad node"))?;
                let name = unesc(
                    tok.next()
                        .ok_or_else(|| perr(ln, "constraint: missing name"))?,
                )
                .map_err(|what| perr(ln, format!("constraint: bad name: {what}")))?;
                require_end(&mut tok, ln, "constraint")?;
                b.constraint(NodeId(node), kind, &name)?;
            }
            "tag" => {
                let tag = match tok.next() {
                    Some("multi_fidelity") => {
                        let levels = tok
                            .next()
                            .and_then(|t| t.parse().ok())
                            .ok_or_else(|| perr(ln, "tag: bad levels"))?;
                        ProblemTag::MultiFidelity { levels }
                    }
                    Some("chance") => {
                        let prob = tok
                            .next()
                            .and_then(parse_f64_hex)
                            .ok_or_else(|| perr(ln, "tag: bad prob"))?;
                        ProblemTag::ChanceConstrained { prob }
                    }
                    Some("bilevel") => {
                        let inner_hash = tok
                            .next()
                            .and_then(|t| u64::from_str_radix(t, 16).ok())
                            .ok_or_else(|| perr(ln, "tag: bad hash"))?;
                        ProblemTag::Bilevel { inner_hash }
                    }
                    _ => return Err(perr(ln, "unknown tag")),
                };
                require_end(&mut tok, ln, "tag")?;
                b.tag(tag);
            }
            "budget" => {
                if saw_budget {
                    return Err(perr(ln, "duplicate budget directive"));
                }
                let n = tok
                    .next()
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "budget: bad count"))?;
                require_end(&mut tok, ln, "budget")?;
                b.set_budget(n);
                saw_budget = true;
            }
            "hash" => {
                let hash = tok
                    .next()
                    .and_then(|t| u64::from_str_radix(t, 16).ok())
                    .ok_or_else(|| perr(ln, "hash: bad value"))?;
                require_end(&mut tok, ln, "hash")?;
                expected_hash = Some(hash);
            }
            "" => return Err(perr(ln, "blank lines are not canonical fsopt syntax")),
            other => return Err(perr(ln, format!("unknown directive `{other}`"))),
        }
    }
    if !saw_budget {
        return Err(perr(hash_line + 1, "missing budget directive"));
    }
    let problem = b.finish();
    let h = expected_hash.ok_or_else(|| perr(hash_line + 1, "hash: bad value"))?;
    let mut body_text = lines[..hash_line].join("\n");
    body_text.push('\n');
    let actual = fnv1a(body_text.as_bytes());
    if actual != h {
        return Err(perr(
            hash_line + 1,
            format!("integrity hash mismatch: file {h:016X}, content {actual:016X}"),
        ));
    }
    let canonical_source =
        serialize_for_wire(&problem, source_version, TokenEncoding::PercentBytes);
    let historical_v1_source = (source_version == WireVersion::V1)
        .then(|| serialize_for_wire(&problem, WireVersion::V1, TokenEncoding::LegacyLiteralUtf8));
    if text != canonical_source
        && historical_v1_source
            .as_deref()
            .is_none_or(|historical| text != historical)
    {
        let line = first_difference_line(text, &canonical_source);
        return Err(perr(
            line,
            "artifact is not an exact encoding emitted by a known writer for its declared wire version; refusing to certify unrelated normalization as AppendMoleZero",
        ));
    }
    let migration = if source_version == WireVersion::V1 {
        let canonical = serialize(&problem);
        Some(DimensionCrosswalkReceipt {
            source_version,
            target_version: WireVersion::V2,
            old_hash: hash_bytes(text.as_bytes()),
            new_hash: hash_bytes(canonical.as_bytes()),
            rule: FiveToSixRule::AppendMoleZero,
        })
    } else {
        None
    };
    Ok(ParsedProblem {
        problem,
        source_version,
        source_hash: h,
        migration,
    })
}

#[allow(clippy::too_many_lines)] // One auditable block keeps every grammar arm and refusal position together.
fn parse_expr(
    b: &mut ProblemBuilder,
    toks: &[&str],
    ln: usize,
    version: WireVersion,
) -> Result<NodeId, OptError> {
    let require_arity = |expected: usize, kind: &str| -> Result<(), OptError> {
        if toks.len() != expected {
            return Err(perr(
                ln,
                format!(
                    "expr {kind}: expected {} operand(s), found {}",
                    expected - 1,
                    toks.len().saturating_sub(1)
                ),
            ));
        }
        Ok(())
    };
    let get = |i: usize| -> Result<&str, OptError> {
        toks.get(i)
            .copied()
            .ok_or_else(|| perr(ln, "expr: missing operand"))
    };
    let node = |i: usize| -> Result<NodeId, OptError> {
        Ok(NodeId(
            get(i)?
                .parse()
                .map_err(|_| perr(ln, "expr: bad node operand"))?,
        ))
    };
    match get(0)? {
        "var" => {
            require_arity(2, "var")?;
            let v: u32 = get(1)?.parse().map_err(|_| perr(ln, "var: bad id"))?;
            b.var_ref(VarId(v))
        }
        "component" => {
            require_arity(3, "component")?;
            let of = node(1)?;
            let index = get(2)?
                .parse()
                .map_err(|_| perr(ln, "component: bad index"))?;
            b.component(of, index)
        }
        "const" => {
            require_arity(3, "const")?;
            let value = parse_f64_hex(get(1)?).ok_or_else(|| perr(ln, "const: bad value"))?;
            let dims = parse_dims(get(2)?, version).ok_or_else(|| perr(ln, "const: bad dims"))?;
            Ok(b.konst(value, dims))
        }
        "add" => {
            require_arity(3, "add")?;
            b.add(node(1)?, node(2)?)
        }
        "sub" => {
            require_arity(3, "sub")?;
            b.sub(node(1)?, node(2)?)
        }
        "mul" => {
            require_arity(3, "mul")?;
            b.mul(node(1)?, node(2)?)
        }
        "div" => {
            require_arity(3, "div")?;
            b.div(node(1)?, node(2)?)
        }
        "neg" => {
            require_arity(2, "neg")?;
            b.neg(node(1)?)
        }
        "powi" => {
            require_arity(3, "powi")?;
            let base = node(1)?;
            let exp = get(2)?
                .parse()
                .map_err(|_| perr(ln, "powi: bad exponent"))?;
            b.powi(base, exp)
        }
        "sqrt" => {
            require_arity(2, "sqrt")?;
            b.sqrt(node(1)?)
        }
        "exp" => {
            require_arity(2, "exp")?;
            b.exp(node(1)?)
        }
        "ln" => {
            require_arity(2, "ln")?;
            b.ln(node(1)?)
        }
        "tanh" => {
            require_arity(2, "tanh")?;
            b.tanh(node(1)?)
        }
        "dot" => {
            require_arity(3, "dot")?;
            b.dot(node(1)?, node(2)?)
        }
        "norm_sq" => {
            require_arity(2, "norm_sq")?;
            b.norm_sq(node(1)?)
        }
        "min" => {
            require_arity(3, "min")?;
            b.min_of(node(1)?, node(2)?)
        }
        "max" => {
            require_arity(3, "max")?;
            b.max_of(node(1)?, node(2)?)
        }
        "abs" => {
            require_arity(2, "abs")?;
            b.abs(node(1)?)
        }
        "pde_residual" => {
            require_arity(5, "pde_residual")?;
            let study =
                unesc(get(1)?).map_err(|what| perr(ln, format!("pde: bad study token: {what}")))?;
            let over: u32 = get(2)?.parse().map_err(|_| perr(ln, "pde: bad var"))?;
            let adj = match get(3)? {
                "0" => false,
                "1" => true,
                _ => return Err(perr(ln, "pde: adjoint flag must be exactly 0 or 1")),
            };
            let dims = parse_dims(get(4)?, version).ok_or_else(|| perr(ln, "pde: bad dims"))?;
            b.pde_residual(&study, VarId(over), adj, dims)
        }
        "expectation" => {
            require_arity(3, "expectation")?;
            let of = node(1)?;
            let cfg = unesc(get(2)?)
                .map_err(|what| perr(ln, format!("expectation: bad config token: {what}")))?;
            b.expectation(of, &cfg)
        }
        "cvar" => {
            require_arity(4, "cvar")?;
            let of = node(1)?;
            let alpha = parse_f64_hex(get(2)?).ok_or_else(|| perr(ln, "cvar: bad alpha"))?;
            let cfg = unesc(get(3)?)
                .map_err(|what| perr(ln, format!("cvar: bad config token: {what}")))?;
            b.cvar(of, alpha, &cfg)
        }
        "quantile" => {
            require_arity(4, "quantile")?;
            let of = node(1)?;
            let q = parse_f64_hex(get(2)?).ok_or_else(|| perr(ln, "quantile: bad q"))?;
            let cfg = unesc(get(3)?)
                .map_err(|what| perr(ln, format!("quantile: bad config token: {what}")))?;
            b.quantile(of, q, &cfg)
        }
        other => Err(perr(ln, format!("unknown expr kind `{other}`"))),
    }
}

#[cfg(test)]
mod tests {
    use super::{
        ContentHash, FiveToSixRule, TokenEncoding, WireVersion, esc, fnv1a, parse,
        parse_with_version, problem_hash, serialize, serialize_for_wire, unesc,
    };
    use crate::ir::{NodeId, ProblemBuilder, Sense};
    use fs_qty::Dims;

    fn with_hash(body: &str) -> String {
        assert!(body.ends_with('\n'));
        format!("{body}hash {:016X}\n", fnv1a(body.as_bytes()))
    }

    #[test]
    fn esc_unesc_round_trips_utf8_and_delimiters() {
        // Non-ASCII, delimiters, and control bytes must round-trip EXACTLY —
        // the old byte-wise `as char` decode was a lossy Latin-1 pass that
        // corrupted every multibyte field (breaking the BITWISE round-trip
        // contract) and could panic slicing a str at a non-char boundary.
        for s in [
            "café",
            "α+β·γ",
            "a b\n%c",
            "π²∇🎯",
            "plain-ascii_123",
            "",
            "100%",
        ] {
            assert_eq!(
                unesc(&esc(s)).expect("writer emits valid escapes"),
                *s,
                "round-trip failed for {s:?}"
            );
        }
        // ASCII fields keep the exact prior wire format (backward compatible).
        assert_eq!(esc("a b%c\n"), "a%20b%25c%0A");
        assert_eq!(unesc("a%20b%25c%0A").expect("canonical escapes"), "a b%c\n");
        for malformed in ["%\u{20AC}x", "%", "%2", "%GG", "%FF"] {
            assert!(
                unesc(malformed).is_err(),
                "malformed or invalid-UTF-8 escape must refuse: {malformed:?}"
            );
        }
    }

    #[test]
    fn exact_v1_bytes_decode_with_zero_amount_and_source_provenance() {
        const LEGACY_HASH: u64 = 0xEA73_E3CB_2B7D_E122;
        let legacy_body = concat!(
            "fsopt v1\n",
            "expr 0 const 3FF0000000000000 (1,2,3,4,5)\n",
            "objective min 0 3FF0000000000000\n",
            "budget 0\n",
        );
        let legacy_hash = fnv1a(legacy_body.as_bytes());
        assert_eq!(
            legacy_hash, LEGACY_HASH,
            "legacy bytes and hash are immutable"
        );
        let legacy_text = format!("{legacy_body}hash {legacy_hash:016X}\n");

        let decoded = parse_with_version(&legacy_text).expect("exact v1 bytes decode");
        assert_eq!(decoded.source_version(), WireVersion::V1);
        assert_eq!(decoded.source_hash(), legacy_hash);
        assert_eq!(
            decoded.problem().node_dims(NodeId(0)),
            Dims([1, 2, 3, 4, 5, 0]),
            "legacy inputs acquire an explicit zero amount exponent"
        );
        assert_eq!(
            serialize_for_wire(
                decoded.problem(),
                WireVersion::V1,
                TokenEncoding::PercentBytes,
            ),
            legacy_text,
            "the explicit-v1 canonical serializer preserves the pinned artifact"
        );

        let canonical = serialize(decoded.problem());
        assert!(canonical.starts_with("fsopt v2\n"));
        assert!(canonical.contains("(1,2,3,4,5,0)"));
        let receipt = decoded.migration().expect("v1 receipt is mandatory");
        assert_eq!(receipt.source_version(), WireVersion::V1);
        assert_eq!(receipt.target_version(), WireVersion::V2);
        assert_eq!(receipt.rule(), FiveToSixRule::AppendMoleZero);
        assert_eq!(
            receipt.old_hash(),
            ContentHash::from_hex(
                "3fb7c9e7e9e1c827030c9dffcc6d5a5139b93e4a13910d2fd7a1f402b384a963"
            )
            .expect("pinned complete-v1 BLAKE3")
        );
        assert_eq!(
            receipt.new_hash(),
            ContentHash::from_hex(
                "96d79826aea409f01491608a789834c2dd90bf70bb8e455347ba04b90b21c7b6"
            )
            .expect("pinned canonical-v2 BLAKE3")
        );
        assert!(receipt.verifies(legacy_text.as_bytes(), canonical.as_bytes()));
        assert!(
            !receipt.verifies(legacy_body.as_bytes(), canonical.as_bytes()),
            "the old hash binds the terminal FNV line too"
        );
        let canonical_hash_line = canonical.rfind("hash ").expect("canonical hash line");
        assert!(
            !receipt.verifies(
                legacy_text.as_bytes(),
                &canonical.as_bytes()[..canonical_hash_line]
            ),
            "the new hash binds the complete canonical artifact"
        );
        assert!(!receipt.verifies(b"tampered", canonical.as_bytes()));
        assert!(
            parse(&legacy_text)
                .expect_err("receipt-discarding parse must reject v1")
                .to_string()
                .contains("semantic-crosswalk receipt")
        );
        let reparsed = parse_with_version(&canonical).expect("canonical v2 bytes decode");
        assert_eq!(reparsed.source_version(), WireVersion::V2);
        assert_eq!(reparsed.source_hash(), problem_hash(decoded.problem()));
        assert!(reparsed.migration().is_none());
        assert_eq!(reparsed.problem(), decoded.problem());
    }

    #[test]
    fn original_literal_utf8_v1_bytes_decode_strictly_with_a_receipt() {
        const HISTORICAL_V1_HASH: u64 = 0xE4D3_9AEE_EEEB_FE36;
        // The original v1 writer in 293b1c1 emitted non-ASCII token bytes
        // literally. The exact historical grammar remains recognizable even
        // though current v2 percent-encodes every non-ASCII byte.
        let historical_body = concat!(
            "fsopt v1\n",
            "var 0 café Rn(1) (0,0,0,0,0)\n",
            "expr 0 const 0000000000000000 (0,0,0,0,0)\n",
            "constraint eq0 0 tail\r\n",
            "budget 0\n",
        );
        assert_eq!(fnv1a(historical_body.as_bytes()), HISTORICAL_V1_HASH);
        let historical = with_hash(historical_body);
        let decoded = parse_with_version(&historical).expect("literal-UTF-8 v1 decodes");
        assert_eq!(decoded.source_version(), WireVersion::V1);
        assert_eq!(decoded.source_hash(), HISTORICAL_V1_HASH);
        assert_eq!(decoded.problem().vars[0].name, "café");
        assert_eq!(decoded.problem().vars[0].dims, Dims::NONE);
        assert_eq!(decoded.problem().constraints[0].name, "tail\r");
        assert_eq!(
            serialize_for_wire(
                decoded.problem(),
                WireVersion::V1,
                TokenEncoding::LegacyLiteralUtf8,
            ),
            historical,
            "the original-v1 writer preserves the pinned artifact"
        );
        let canonical = serialize(decoded.problem());
        assert!(canonical.contains("caf%C3%A9"));
        assert!(canonical.contains("tail%0D"));
        let receipt = decoded.migration().expect("historical v1 receipt");
        assert_eq!(
            receipt.old_hash(),
            ContentHash::from_hex(
                "7e364590e79894541c7af6afcdff09aebec5cd8416e2237ccb2db0223c9ade12"
            )
            .expect("pinned literal-UTF-8-v1 BLAKE3")
        );
        assert_eq!(
            receipt.new_hash(),
            ContentHash::from_hex(
                "7eb0ada0b3eb4c4f8f71880fab46853b34c4094c6ce662dc6307f41b4eb752c4"
            )
            .expect("pinned percent-encoded-v2 BLAKE3")
        );
        assert!(receipt.verifies(historical.as_bytes(), canonical.as_bytes()));
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One test is the complete receipt-eligibility refusal matrix.
    fn noncanonical_sources_cannot_receive_append_mole_only_receipts() {
        const BASE_BODY: &str = concat!(
            "fsopt v1\n",
            "expr 0 const 3FF0000000000000 (0,0,0,0,0)\n",
            "objective min 0 3FF0000000000000\n",
            "budget 0\n",
        );
        let rejects = |artifact: &str, case: &str| {
            let Err(error) = parse_with_version(artifact) else {
                panic!("{case} must refuse before issuing a receipt");
            };
            assert!(!error.to_string().is_empty(), "{case} returns a diagnosis");
        };

        rejects(
            &with_hash(concat!(
                "expr 0 const 3FF0000000000000 (0,0,0,0,0)\n",
                "budget 0\n",
            )),
            "headerless artifact without a schema identity",
        );

        rejects(
            &with_hash(concat!(
                "fsopt v1\n",
                "expr 0 const 3FF0000000000000 (0,0,0,0,0)\n",
                "objective min 0 3FF0000000000000 trailing\n",
                "budget 0\n",
            )),
            "extra directive field",
        );
        rejects(
            &with_hash(concat!(
                "fsopt v1\n",
                "expr 0 const 3FF0000000000000 (0,0,0,0,0) trailing\n",
                "budget 0\n",
            )),
            "extra expression operand",
        );
        rejects(
            &with_hash("fsopt v1\nvar 9 x Rn(1) (0,0,0,0,0)\nbudget 0\n"),
            "wrong variable index",
        );
        rejects(
            &with_hash(concat!(
                "fsopt v1\n",
                "var 0 x Rn(1) (0,0,0,0,0)\n",
                "expr 0 pde_residual study 0 2 (0,0,0,0,0)\n",
                "budget 0\n",
            )),
            "non-Boolean adjoint token",
        );
        for escaped_name in ["%", "%2", "%GG", "%FF"] {
            rejects(
                &with_hash(&format!(
                    "fsopt v1\nvar 0 {escaped_name} Rn(1) (0,0,0,0,0)\nbudget 0\n"
                )),
                "malformed or invalid-UTF-8 escape",
            );
        }
        rejects(
            &with_hash(concat!(
                "fsopt v1\n",
                "expr 0 const 3FF0000000000000 (0,0,0,0,0)\n",
                "\n",
                "budget 0\n",
            )),
            "blank body line",
        );
        rejects(
            &with_hash(concat!(
                "fsopt v1\n",
                "expr 0 const 3FF0000000000000 (0,0,0,0,0)\n",
            )),
            "missing budget",
        );
        rejects(
            &with_hash(concat!(
                "fsopt v1\n",
                "expr 0 const 3FF0000000000000 (0,0,0,0,0)\n",
                "budget 0\n",
                "budget 0\n",
            )),
            "duplicate budget",
        );

        let canonical = with_hash(BASE_BODY);
        rejects(&canonical.replace('\n', "\r\n"), "CRLF-normalized artifact");
        rejects(
            canonical
                .strip_suffix('\n')
                .expect("canonical artifact ends with newline"),
            "missing final newline",
        );

        let noncanonical_v2 = with_hash(concat!(
            "fsopt v2\n",
            "expr 0 const 3ff0000000000000 (0,0,0,0,0,0)\n",
            "budget 0\n",
        ));
        rejects(&noncanonical_v2, "noncanonical current-v2 spelling");
        let canonical_v2_body = concat!(
            "fsopt v2\n",
            "expr 0 const 3FF0000000000000 (0,0,0,0,0,0)\n",
            "budget 0\n",
        );
        rejects(
            &with_hash(&canonical_v2_body.replace("expr 0", "expr  0")),
            "v2 whitespace normalization",
        );
        rejects(
            &with_hash(&canonical_v2_body.replace('\n', "\r\n")),
            "v2 line-ending normalization",
        );
    }

    #[test]
    fn v2_writer_and_parser_preserve_amount_dimensions() {
        let mut builder = ProblemBuilder::new();
        let amount = builder.konst(2.0, Dims([0, 0, 0, 0, 0, 1]));
        builder
            .objective(amount, Sense::Minimize, 1.0)
            .expect("amount-valued scalar objective");
        let problem = builder.finish();

        let text = serialize(&problem);
        assert!(text.starts_with("fsopt v2\n"));
        assert!(text.contains("(0,0,0,0,0,1)"));
        assert_eq!(parse(&text).expect("v2 round-trip"), problem);

        let v1_with_six_dims = with_hash(concat!(
            "fsopt v1\n",
            "expr 0 const 3FF0000000000000 (0,0,0,0,0,1)\n",
        ));
        let v2_with_five_dims = with_hash(concat!(
            "fsopt v2\n",
            "expr 0 const 3FF0000000000000 (0,0,0,0,0)\n",
        ));
        let v1_error = parse_with_version(&v1_with_six_dims)
            .expect_err("explicit v1 arity stays exact")
            .to_string();
        assert!(v1_error.contains("line 2") && v1_error.contains("bad dims"));
        let v2_error = parse_with_version(&v2_with_five_dims)
            .expect_err("v2 arity stays exact")
            .to_string();
        assert!(v2_error.contains("line 2") && v2_error.contains("bad dims"));
    }

    #[test]
    fn exactly_one_terminal_hash_directive_is_required() {
        let mut builder = ProblemBuilder::new();
        let scalar = builder.konst(1.0, Dims::NONE);
        builder
            .objective(scalar, Sense::Minimize, 1.0)
            .expect("scalar objective");
        let canonical = serialize(&builder.finish());
        parse(&canonical).expect("canonical artifact has one terminal hash");

        let hash_at = canonical.rfind("hash ").expect("hash directive");
        let body = &canonical[..hash_at];
        let hash_line = &canonical[hash_at..];
        assert!(parse_with_version(body).is_err(), "missing hash must fail");
        assert!(
            parse_with_version(&format!("{canonical}{hash_line}")).is_err(),
            "duplicate hash must fail"
        );
        assert!(
            parse_with_version(&format!("{canonical}budget 1\n")).is_err(),
            "nonterminal hash must fail"
        );
        assert!(
            parse_with_version(&format!("{body}{} trailing\n", hash_line.trim_end())).is_err(),
            "hash trailing fields must fail"
        );
    }
}
