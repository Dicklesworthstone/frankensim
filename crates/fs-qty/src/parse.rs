//! SI unit-expression parsing: `"0.12Pa*s"`, `"0.5L/s"`, `"65deg"`,
//! `"0.061N/m"`, `"0.03m2/s3"`, `"12mm"`, `"2h"`, `"5772.22"` → [`QtyAny`].
//!
//! This is the literal syntax FrankenScript studies use (plan Appendix C);
//! fs-ir's admission checker parses budgets and BC values through this exact
//! grammar, so unit errors die at admission (plan §11.1).
//!
//! Grammar (whitespace tolerated around the number and between factors):
//!
//! ```text
//! qty      := number unit-expr?
//! unit-expr:= factor ( ('*' | '·' | '/') factor )*      // strict left-to-right
//! factor   := symbol exponent?
//! exponent := '^'? '-'? digits                          // m2, m^2, s^-1
//! symbol   := longest-match named unit, else prefix + named unit
//! ```
//!
//! Policy notes (documented no-claims):
//! - Angles: `rad` is dimensionless; `deg` converts by π/180.
//! - `degC` is AFFINE and therefore only legal as a whole, lone unit with
//!   exponent 1 (`"20degC"`); compounds like `degC/s` are rejected with a
//!   teaching error (differences of Celsius are kelvin — say `K/s`).
//! - Information/monetary units (`GiB`, `B`, …) are rejected here with a
//!   pointer to fs-ir budget syntax; they are not physical dimensions.
//! - `mol` is the sixth admitted base dimension; `cd` remains outside the
//!   vector until photometry is real (no-claim).

use crate::{DIMENSION_COUNT, Dims, QtyAny};
use core::fmt;
use fs_blake3::hash_domain;

pub use fs_blake3::ContentHash;

const PARSE_INPUT_IDENTITY_DOMAIN: &str = "frankensim.fs-qty.parse-input.v1";

/// Explicit work and diagnostic-retention budget for one quantity literal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParseBudget {
    max_input_bytes: usize,
    max_factors: usize,
    max_token_bytes: usize,
    max_diagnostic_bytes: usize,
}

impl ParseBudget {
    /// Conservative compatibility budget used by [`parse_qty`].
    pub const DEFAULT: Self = Self::new(4_096, 256, 64, 256);

    /// Construct a budget. Zero is meaningful and tested for every field.
    #[must_use]
    pub const fn new(
        max_input_bytes: usize,
        max_factors: usize,
        max_token_bytes: usize,
        max_diagnostic_bytes: usize,
    ) -> Self {
        Self {
            max_input_bytes,
            max_factors,
            max_token_bytes,
            max_diagnostic_bytes,
        }
    }

    /// Maximum admitted UTF-8 source bytes.
    #[must_use]
    pub const fn max_input_bytes(self) -> usize {
        self.max_input_bytes
    }

    /// Maximum unit factors in one expression.
    #[must_use]
    pub const fn max_factors(self) -> usize {
        self.max_factors
    }

    /// Maximum bytes in a number, unit, or exponent token.
    #[must_use]
    pub const fn max_token_bytes(self) -> usize {
        self.max_token_bytes
    }

    /// Maximum retained UTF-8 excerpt bytes in an error.
    #[must_use]
    pub const fn max_diagnostic_bytes(self) -> usize {
        self.max_diagnostic_bytes
    }
}

impl Default for ParseBudget {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// The bounded parser resource whose admission failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParseResource {
    /// Total source bytes.
    InputBytes,
    /// Number, unit, or exponent token bytes.
    TokenBytes,
    /// Unit factors evaluated.
    Factors,
}

/// Where and why parsing failed, with a suggested fix (errors-as-guidance).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    /// UTF-8-safe bounded source excerpt containing the failure location.
    pub preview: String,
    /// Byte offset of `preview` within the full input.
    pub preview_start: usize,
    /// Full source length without retaining the source.
    pub input_bytes: usize,
    /// Full-input identity when the input passed byte admission. Oversized
    /// input is deliberately not scanned merely to manufacture a hash.
    pub source_hash: Option<ContentHash>,
    /// Byte offset of the failure.
    pub at: usize,
    /// What went wrong.
    pub kind: ParseErrorKind,
    /// A machine-usable suggestion.
    pub help: String,
}

/// Failure classes for unit-expression parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseErrorKind {
    /// No leading number.
    MissingNumber,
    /// Unit token not recognized.
    UnknownUnit(String),
    /// `degC` used in a compound/exponentiated position.
    AffineUnitInCompound,
    /// Information units are not physical dimensions.
    InformationUnit(String),
    /// Trailing garbage after a valid expression.
    TrailingInput,
    /// Exponent did not parse or exceeded the public dimension cap.
    BadExponent,
    /// A literal or unit conversion produced NaN, infinity, or an
    /// unrepresentable underflow to zero.
    NonFiniteValue,
    /// An explicit parser work/retention limit was exceeded.
    BudgetExceeded {
        /// Bounded resource.
        resource: ParseResource,
        /// Configured maximum.
        limit: usize,
        /// Exact observation when available, otherwise a proven lower bound.
        observed_at_least: usize,
    },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "cannot parse quantity excerpt {:?} (source bytes {}..{} of {}, hash {}) at byte {}: {:?}; {}",
            self.preview,
            self.preview_start,
            self.preview_start + self.preview.len(),
            self.input_bytes,
            self.source_hash.map_or_else(
                || "unavailable-before-byte-admission".to_string(),
                |hash| hash.to_string()
            ),
            self.at,
            self.kind,
            self.help
        )
    }
}

impl core::error::Error for ParseError {}

impl ParseError {
    /// Verify that the retained full-input identity belongs to `source`.
    /// Oversized inputs return `false` because byte admission deliberately
    /// precedes the hash pass.
    #[must_use]
    pub fn verifies_source(&self, source: &str) -> bool {
        let Some(expected) = self.source_hash else {
            return false;
        };
        self.input_bytes == source.len()
            && expected == hash_domain(PARSE_INPUT_IDENTITY_DOMAIN, source.as_bytes())
    }
}

/// A named unit: symbol → (scale-to-SI, dimension). `degC` handled specially.
struct Unit {
    symbol: &'static str,
    scale: f64,
    dims: Dims,
}

const D_NONE: Dims = Dims([0, 0, 0, 0, 0, 0]);
const D_M: Dims = Dims([1, 0, 0, 0, 0, 0]);
const D_KG: Dims = Dims([0, 1, 0, 0, 0, 0]);
const D_S: Dims = Dims([0, 0, 1, 0, 0, 0]);
const D_K: Dims = Dims([0, 0, 0, 1, 0, 0]);
const D_A: Dims = Dims([0, 0, 0, 0, 1, 0]);
const D_MOL: Dims = Dims([0, 0, 0, 0, 0, 1]);
const D_N: Dims = Dims([1, 1, -2, 0, 0, 0]);
const D_PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const D_J: Dims = Dims([2, 1, -2, 0, 0, 0]);
const D_W: Dims = Dims([2, 1, -3, 0, 0, 0]);
const D_HZ: Dims = Dims([0, 0, -1, 0, 0, 0]);
const D_M3: Dims = Dims([3, 0, 0, 0, 0, 0]);
const D_V: Dims = Dims([2, 1, -3, 0, -1, 0]);
const D_C: Dims = Dims([0, 0, 1, 0, 1, 0]);
const D_WB: Dims = Dims([2, 1, -2, 0, -1, 0]);
const D_H: Dims = Dims([2, 1, -2, 0, -2, 0]);
const D_OHM: Dims = Dims([2, 1, -3, 0, -2, 0]);
const D_SIEMENS: Dims = Dims([-2, -1, 3, 0, 2, 0]);
const D_F: Dims = Dims([-2, -1, 4, 0, 2, 0]);
const D_T: Dims = Dims([0, 1, -2, 0, -1, 0]);

/// Longest-match table. Order does not matter (lookup takes the longest
/// symbol that matches the whole token before falling back to prefix+unit).
const UNITS: &[Unit] = &[
    Unit {
        symbol: "m",
        scale: 1.0,
        dims: D_M,
    },
    Unit {
        symbol: "g",
        scale: 1e-3,
        dims: D_KG,
    }, // gram; kg arrives via prefix
    Unit {
        symbol: "s",
        scale: 1.0,
        dims: D_S,
    },
    Unit {
        symbol: "K",
        scale: 1.0,
        dims: D_K,
    },
    Unit {
        symbol: "A",
        scale: 1.0,
        dims: D_A,
    },
    Unit {
        symbol: "mol",
        scale: 1.0,
        dims: D_MOL,
    },
    Unit {
        symbol: "N",
        scale: 1.0,
        dims: D_N,
    },
    Unit {
        symbol: "Pa",
        scale: 1.0,
        dims: D_PA,
    },
    Unit {
        symbol: "J",
        scale: 1.0,
        dims: D_J,
    },
    Unit {
        symbol: "W",
        scale: 1.0,
        dims: D_W,
    },
    Unit {
        symbol: "V",
        scale: 1.0,
        dims: D_V,
    },
    Unit {
        symbol: "C",
        scale: 1.0,
        dims: D_C,
    },
    Unit {
        symbol: "Wb",
        scale: 1.0,
        dims: D_WB,
    },
    Unit {
        symbol: "H",
        scale: 1.0,
        dims: D_H,
    },
    Unit {
        symbol: "Ohm",
        scale: 1.0,
        dims: D_OHM,
    },
    Unit {
        symbol: "S",
        scale: 1.0,
        dims: D_SIEMENS,
    },
    Unit {
        symbol: "F",
        scale: 1.0,
        dims: D_F,
    },
    Unit {
        symbol: "T",
        scale: 1.0,
        dims: D_T,
    },
    Unit {
        symbol: "Hz",
        scale: 1.0,
        dims: D_HZ,
    },
    Unit {
        symbol: "L",
        scale: 1e-3,
        dims: D_M3,
    },
    Unit {
        symbol: "min",
        scale: 60.0,
        dims: D_S,
    },
    Unit {
        symbol: "h",
        scale: 3600.0,
        dims: D_S,
    },
    Unit {
        symbol: "rad",
        scale: 1.0,
        dims: D_NONE,
    },
    Unit {
        symbol: "deg",
        scale: core::f64::consts::PI / 180.0,
        dims: D_NONE,
    },
    Unit {
        symbol: "%",
        scale: 1e-2,
        dims: D_NONE,
    },
];

/// SI prefixes accepted before a named unit.
const PREFIXES: &[(&str, f64)] = &[
    ("p", 1e-12),
    ("n", 1e-9),
    ("u", 1e-6),
    ("µ", 1e-6),
    ("m", 1e-3),
    ("c", 1e-2),
    ("d", 1e-1),
    ("k", 1e3),
    ("M", 1e6),
    ("G", 1e9),
    ("T", 1e12),
];

/// Information-unit symbols we explicitly refuse with a teaching error.
const INFORMATION_UNITS: &[&str] = &["B", "iB", "KiB", "MiB", "GiB", "TiB", "bit"];

/// Public parser cap for every syntactic and accumulated dimension exponent.
/// Keeping this wider than the stored `i8` representation is deliberate: all
/// arithmetic and cap checks happen before the final narrowing conversion.
const MAX_ABS_UNIT_EXPONENT: i32 = 60;

#[derive(Clone, Copy)]
struct ParseContext<'a> {
    input: &'a str,
    budget: ParseBudget,
    source_hash: ContentHash,
}

impl<'a> ParseContext<'a> {
    fn admitted(input: &'a str, budget: ParseBudget) -> Self {
        Self {
            input,
            budget,
            source_hash: hash_domain(PARSE_INPUT_IDENTITY_DOMAIN, input.as_bytes()),
        }
    }
}

fn diagnostic_excerpt(input: &str, at: usize, max_bytes: usize) -> (usize, String) {
    let mut at = at.min(input.len());
    while !input.is_char_boundary(at) {
        at -= 1;
    }
    if max_bytes == 0 {
        return (at, String::new());
    }
    let mut start = at.saturating_sub(max_bytes / 2);
    while start < at && !input.is_char_boundary(start) {
        start += 1;
    }
    let mut end = start.saturating_add(max_bytes).min(input.len());
    while end > start && !input.is_char_boundary(end) {
        end -= 1;
    }
    (start, input[start..end].to_string())
}

fn error_from_source(
    input: &str,
    budget: ParseBudget,
    source_hash: Option<ContentHash>,
    at: usize,
    kind: ParseErrorKind,
    help: &str,
) -> ParseError {
    let at = at.min(input.len());
    let (preview_start, preview) = diagnostic_excerpt(input, at, budget.max_diagnostic_bytes());
    ParseError {
        preview,
        preview_start,
        input_bytes: input.len(),
        source_hash,
        at,
        kind,
        help: help.to_string(),
    }
}

fn err(ctx: &ParseContext<'_>, at: usize, kind: ParseErrorKind, help: &str) -> ParseError {
    error_from_source(ctx.input, ctx.budget, Some(ctx.source_hash), at, kind, help)
}

fn budget_error(
    ctx: &ParseContext<'_>,
    at: usize,
    resource: ParseResource,
    limit: usize,
    observed_at_least: usize,
    help: &str,
) -> ParseError {
    err(
        ctx,
        at,
        ParseErrorKind::BudgetExceeded {
            resource,
            limit,
            observed_at_least,
        },
        help,
    )
}

fn trim_start_at<'a>(text: &'a str, pos: &mut usize) -> &'a str {
    let trimmed = text.trim_start();
    *pos += text.len() - trimmed.len();
    trimmed
}

fn nonfinite(ctx: &ParseContext<'_>, at: usize, help: &str) -> ParseError {
    err(ctx, at, ParseErrorKind::NonFiniteValue, help)
}

fn significand_is_nonzero(number: &str) -> bool {
    let exponent_at = number
        .find('e')
        .or_else(|| number.find('E'))
        .unwrap_or(number.len());
    number[..exponent_at]
        .bytes()
        .any(|digit| matches!(digit, b'1'..=b'9'))
}

/// Resolve one unit token (no exponent) to (scale, dims).
fn resolve_token(ctx: &ParseContext<'_>, at: usize, tok: &str) -> Result<(f64, Dims), ParseError> {
    // Whole-token named unit wins (so `min` is minutes, not milli-inches).
    if let Some(u) = UNITS.iter().find(|u| u.symbol == tok) {
        return Ok((u.scale, u.dims));
    }
    // Information units get a dedicated refusal.
    if INFORMATION_UNITS.iter().any(|s| tok.ends_with(s)) {
        return Err(err(
            ctx,
            at,
            ParseErrorKind::InformationUnit(tok.to_string()),
            "information units (bytes) are not physical dimensions; memory/time budgets \
             use fs-ir budget syntax, e.g. (mem 96GiB)",
        ));
    }
    // Prefix + named unit (prefix is at most one char here; `da` unsupported).
    let mut chars = tok.char_indices();
    if let Some((_, first)) = chars.next() {
        let rest_start = chars.next().map_or(tok.len(), |(i, _)| i);
        let rest = &tok[rest_start..];
        if !rest.is_empty()
            && let Some(&(_, scale)) = PREFIXES.iter().find(|(p, _)| p.starts_with(first))
            && let Some(u) = UNITS.iter().find(|u| u.symbol == rest)
        {
            return Ok((scale * u.scale, u.dims));
        }
    }
    Err(err(
        ctx,
        at,
        ParseErrorKind::UnknownUnit(tok.to_string()),
        "expected an SI unit like m, kg, s, K, A, mol, N, Pa, J, W, V, C, Wb, H, Ohm, S, F, T, \
         Hz, L, min, h, rad, deg, % \
         with an optional prefix (p n u m c d k M G T)",
    ))
}

/// Scan the leading number of a quantity literal; returns its byte length.
fn scan_number(ctx: &ParseContext<'_>, s: &str, at: usize) -> Result<usize, ParseError> {
    let bytes = s.as_bytes();
    let mut end = 0;
    let mut seen_digit = false;
    while end < bytes.len() {
        let c = bytes[end] as char;
        let is_num = c.is_ascii_digit()
            || c == '.'
            || (end == 0 && (c == '+' || c == '-'))
            || ((c == 'e' || c == 'E')
                && seen_digit
                && bytes
                    .get(end + 1)
                    .is_some_and(|&n| (n as char).is_ascii_digit() || n == b'+' || n == b'-'));
        if !is_num {
            break;
        }
        if c.is_ascii_digit() {
            seen_digit = true;
        }
        let advance = if c == 'e' || c == 'E' {
            2 // consume the sign/digit that justified accepting 'e'
        } else {
            1
        };
        let next = end + advance;
        if next > ctx.budget.max_token_bytes() {
            return Err(budget_error(
                ctx,
                at,
                ParseResource::TokenBytes,
                ctx.budget.max_token_bytes(),
                next,
                "numeric token exceeds the parser budget; shorten the literal",
            ));
        }
        end = next;
    }
    Ok(end)
}

fn scan_unit_token(ctx: &ParseContext<'_>, rest: &str, at: usize) -> Result<usize, ParseError> {
    let mut token_bytes = 0;
    for (index, symbol) in rest.char_indices() {
        if !(symbol.is_alphabetic() || symbol == 'µ' || symbol == '%') {
            break;
        }
        let next = index + symbol.len_utf8();
        if next > ctx.budget.max_token_bytes() {
            return Err(budget_error(
                ctx,
                at,
                ParseResource::TokenBytes,
                ctx.budget.max_token_bytes(),
                next,
                "unit token exceeds the parser budget; use a supported short SI symbol",
            ));
        }
        token_bytes = next;
    }
    Ok(token_bytes)
}

/// Parse an optional exponent (`^`? `-`? digits) at the head of `rest`.
/// Returns `(exponent, bytes_consumed)`; the default exponent is 1. Parsing and
/// validation stay in a wider integer so `-128` can never reach `i8::MIN`.
fn parse_exponent(
    ctx: &ParseContext<'_>,
    pos: usize,
    rest: &str,
) -> Result<(i32, usize), ParseError> {
    let mut r = rest;
    let mut consumed = 0;
    if let Some(stripped) = r.strip_prefix('^') {
        r = stripped;
        consumed += 1;
        if consumed > ctx.budget.max_token_bytes() {
            return Err(budget_error(
                ctx,
                pos,
                ParseResource::TokenBytes,
                ctx.budget.max_token_bytes(),
                consumed,
                "unit exponent token exceeds the parser budget",
            ));
        }
    }
    let neg = if let Some(stripped) = r.strip_prefix('-') {
        r = stripped;
        consumed += 1;
        if consumed > ctx.budget.max_token_bytes() {
            return Err(budget_error(
                ctx,
                pos,
                ParseResource::TokenBytes,
                ctx.budget.max_token_bytes(),
                consumed,
                "unit exponent token exceeds the parser budget",
            ));
        }
        true
    } else {
        false
    };
    let mut dig_len = 0;
    for (index, digit) in r.char_indices() {
        if !digit.is_ascii_digit() {
            break;
        }
        let next = index + digit.len_utf8();
        if consumed.saturating_add(next) > ctx.budget.max_token_bytes() {
            return Err(budget_error(
                ctx,
                pos,
                ParseResource::TokenBytes,
                ctx.budget.max_token_bytes(),
                consumed.saturating_add(next),
                "unit exponent token exceeds the parser budget",
            ));
        }
        dig_len = next;
    }
    if dig_len == 0 {
        if consumed > 0 {
            return Err(err(
                ctx,
                pos,
                ParseErrorKind::BadExponent,
                "dangling ^ or - without digits",
            ));
        }
        return Ok((1, 0));
    }
    let magnitude: u64 = r[..dig_len].parse().map_err(|_| {
        err(
            ctx,
            pos,
            ParseErrorKind::BadExponent,
            "unit exponent is too large to parse; use an integer from -60 through 60",
        )
    })?;
    if magnitude > MAX_ABS_UNIT_EXPONENT as u64 {
        return Err(err(
            ctx,
            pos,
            ParseErrorKind::BadExponent,
            "unit exponent exceeds the supported ±60 cap; reduce the exponent",
        ));
    }
    let magnitude = i32::try_from(magnitude).map_err(|_| {
        err(
            ctx,
            pos,
            ParseErrorKind::BadExponent,
            "unit exponent is outside the supported wider-integer domain",
        )
    })?;
    let exp = if neg {
        magnitude.checked_neg().ok_or_else(|| {
            err(
                ctx,
                pos,
                ParseErrorKind::BadExponent,
                "negative unit exponent is outside the supported wider-integer domain",
            )
        })?
    } else {
        magnitude
    };
    Ok((exp, consumed + dig_len))
}

fn checked_dims_after_factor(
    ctx: &ParseContext<'_>,
    at: usize,
    dims: [i32; DIMENSION_COUNT],
    factor_dims: Dims,
    exponent: i32,
    divide: bool,
) -> Result<[i32; DIMENSION_COUNT], ParseError> {
    let mut next = dims;
    for index in 0..DIMENSION_COUNT {
        let scaled = i32::from(factor_dims.0[index])
            .checked_mul(exponent)
            .ok_or_else(|| {
                err(
                    ctx,
                    at,
                    ParseErrorKind::BadExponent,
                    "unit-dimension scaling overflowed; reduce the exponent",
                )
            })?;
        let signed = if divide {
            scaled.checked_neg().ok_or_else(|| {
                err(
                    ctx,
                    at,
                    ParseErrorKind::BadExponent,
                    "unit-dimension division overflowed; reduce the exponent",
                )
            })?
        } else {
            scaled
        };
        let accumulated = dims[index].checked_add(signed).ok_or_else(|| {
            err(
                ctx,
                at,
                ParseErrorKind::BadExponent,
                "accumulated unit dimension overflowed; shorten the unit chain",
            )
        })?;
        if accumulated.unsigned_abs() > MAX_ABS_UNIT_EXPONENT as u32 {
            return Err(err(
                ctx,
                at,
                ParseErrorKind::BadExponent,
                "accumulated unit exponent exceeds the supported ±60 cap; check for a runaway unit chain",
            ));
        }
        next[index] = accumulated;
    }
    Ok(next)
}

fn narrow_dims(
    ctx: &ParseContext<'_>,
    at: usize,
    dims: [i32; DIMENSION_COUNT],
) -> Result<Dims, ParseError> {
    let mut narrowed = [0i8; DIMENSION_COUNT];
    for (target, exponent) in narrowed.iter_mut().zip(dims) {
        *target = i8::try_from(exponent).map_err(|_| {
            err(
                ctx,
                at,
                ParseErrorKind::BadExponent,
                "validated unit exponent could not be represented; reduce the exponent",
            )
        })?;
    }
    Ok(Dims(narrowed))
}

/// Parse a quantity literal into a [`QtyAny`].
///
/// # Errors
/// Returns a [`ParseError`] with position, bounded source diagnostics, kind,
/// and a suggested fix. This compatibility entry point always applies
/// [`ParseBudget::DEFAULT`].
pub fn parse_qty(input: &str) -> Result<QtyAny, ParseError> {
    parse_qty_with_budget(input, ParseBudget::DEFAULT)
}

/// Parse a quantity literal under an explicit work and diagnostic budget.
///
/// Byte admission happens before trimming, hashing, or token scanning. Errors
/// for admitted input carry its full content hash; an oversized input carries
/// its exact length and bounded excerpt but is not scanned merely to hash it.
///
/// # Errors
/// Returns a structured [`ParseError`] for syntax, representability, or budget
/// refusal.
#[allow(clippy::too_many_lines)] // Keeping the strict left-to-right grammar and refusal positions together is auditable.
pub fn parse_qty_with_budget(input: &str, budget: ParseBudget) -> Result<QtyAny, ParseError> {
    if input.len() > budget.max_input_bytes() {
        return Err(error_from_source(
            input,
            budget,
            None,
            budget.max_input_bytes(),
            ParseErrorKind::BudgetExceeded {
                resource: ParseResource::InputBytes,
                limit: budget.max_input_bytes(),
                observed_at_least: input.len(),
            },
            "quantity source exceeds the byte budget; shorten or pre-admit the literal",
        ));
    }
    let ctx = ParseContext::admitted(input, budget);
    let s = input.trim();
    let base = input.len() - input.trim_start().len();

    let end = scan_number(&ctx, s, base)?;
    let num: f64 = s[..end].parse().map_err(|_| {
        err(
            &ctx,
            base,
            ParseErrorKind::MissingNumber,
            "a quantity starts with a number, e.g. 0.12Pa*s",
        )
    })?;
    if !num.is_finite() {
        return Err(nonfinite(
            &ctx,
            base,
            "quantity literals must start with a finite number; reduce the magnitude",
        ));
    }
    if num == 0.0 && significand_is_nonzero(&s[..end]) {
        return Err(nonfinite(
            &ctx,
            base,
            "the nonzero decimal significand underflowed to zero; use a representable magnitude",
        ));
    }
    let mut pos = base + end;
    let mut rest = trim_start_at(&s[end..], &mut pos);

    // --- bare number: dimensionless ---
    if rest.is_empty() {
        return Ok(QtyAny::dimensionless(num));
    }

    // --- special-case lone affine unit degC ---
    if rest == "degC" {
        if budget.max_factors() == 0 {
            return Err(budget_error(
                &ctx,
                pos,
                ParseResource::Factors,
                0,
                1,
                "unit expression has too many factors; simplify the dimensions",
            ));
        }
        if rest.len() > budget.max_token_bytes() {
            return Err(budget_error(
                &ctx,
                pos,
                ParseResource::TokenBytes,
                budget.max_token_bytes(),
                rest.len(),
                "unit token exceeds the parser budget; use a supported short SI symbol",
            ));
        }
        let value = num + 273.15;
        if !value.is_finite() {
            return Err(nonfinite(
                &ctx,
                pos,
                "converting degrees Celsius to coherent kelvin overflowed; reduce the magnitude",
            ));
        }
        return Ok(QtyAny::new(value, D_K));
    }

    // --- unit expression, strict left-to-right ---
    let mut value = num;
    let mut dims = [0i32; DIMENSION_COUNT];
    let mut divide = false;
    let mut factors = 0usize;
    loop {
        let factor_at = pos;
        if factors >= budget.max_factors() {
            return Err(budget_error(
                &ctx,
                factor_at,
                ParseResource::Factors,
                budget.max_factors(),
                factors.saturating_add(1),
                "unit expression has too many factors; simplify the dimensions",
            ));
        }
        factors += 1;
        // token = leading unit letters/µ/%
        let tok_len = scan_unit_token(&ctx, rest, pos)?;
        if tok_len == 0 {
            return Err(err(
                &ctx,
                pos,
                ParseErrorKind::TrailingInput,
                "expected a unit symbol here",
            ));
        }
        let tok = &rest[..tok_len];
        if tok.contains("degC") {
            return Err(err(
                &ctx,
                pos,
                ParseErrorKind::AffineUnitInCompound,
                "degC is affine and only legal alone (e.g. \"20degC\"); temperature \
                 differences and rates are kelvin — write K or K/s",
            ));
        }
        let (scale, tok_dims) = resolve_token(&ctx, pos, tok)?;
        rest = &rest[tok_len..];
        pos += tok_len;

        // optional exponent: '^'? '-'? digits
        let exponent_at = pos;
        let (exp, consumed) = parse_exponent(&ctx, exponent_at, rest)?;
        rest = &rest[consumed..];
        pos += consumed;

        let arithmetic_at = if consumed == 0 {
            factor_at
        } else {
            exponent_at
        };
        dims = checked_dims_after_factor(&ctx, arithmetic_at, dims, tok_dims, exp, divide)?;

        // Apply the numerical scale only after the dimension metadata is known
        // to be admissible. No nonfinite intermediate may enter QtyAny or IR.
        let factor_scale = crate::powi_pinned(scale, exp);
        if !factor_scale.is_finite() || factor_scale <= 0.0 {
            return Err(nonfinite(
                &ctx,
                arithmetic_at,
                "raising this positive unit scale to its exponent overflowed or underflowed to zero; reduce the exponent or prefix",
            ));
        }
        let next_value = if divide {
            value / factor_scale
        } else {
            value * factor_scale
        };
        if !next_value.is_finite() || (value != 0.0 && next_value == 0.0) {
            return Err(nonfinite(
                &ctx,
                arithmetic_at,
                "applying this unit factor produced a nonfinite value or underflowed a nonzero quantity to zero; reduce the literal magnitude, prefix, or exponent",
            ));
        }
        value = next_value;

        // separator or end
        rest = trim_start_at(rest, &mut pos);
        match rest.chars().next() {
            None => return Ok(QtyAny::new(value, narrow_dims(&ctx, pos, dims)?)),
            Some(c @ ('*' | '·')) => {
                divide = false;
                pos += c.len_utf8();
                rest = trim_start_at(&rest[c.len_utf8()..], &mut pos);
            }
            Some('/') => {
                divide = true;
                pos += 1;
                rest = trim_start_at(&rest[1..], &mut pos);
            }
            Some(_) => {
                return Err(err(
                    &ctx,
                    pos,
                    ParseErrorKind::TrailingInput,
                    "expected *, ·, / or end of input after a unit factor",
                ));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{DynViscosity, Pressure, SurfaceTension, Time, VolumetricFlowRate};

    /// The Appendix C literal battery: every unit literal that appears in the
    /// plan's example studies must parse to the right value and dimension.
    #[test]
    fn appendix_c_literals() {
        let cases: &[(&str, f64, Dims)] = &[
            ("0.12Pa*s", 0.12, DynViscosity::DIMS),
            ("0.061N/m", 0.061, SurfaceTension::DIMS),
            ("0.5L/s", 5e-4, VolumetricFlowRate::DIMS),
            ("3mm", 3e-3, Dims([1, 0, 0, 0, 0, 0])),
            ("12mm", 12e-3, Dims([1, 0, 0, 0, 0, 0])),
            ("0deg", 0.0, Dims::NONE),
            ("65deg", 65.0 * core::f64::consts::PI / 180.0, Dims::NONE),
            ("3s", 3.0, Time::DIMS),
            ("2h", 7200.0, Time::DIMS),
            ("0.03m2/s3", 0.03, Dims([2, 0, -3, 0, 0, 0])),
            ("15rad/s", 15.0, Dims([0, 0, -1, 0, 0, 0])),
            ("8m/s", 8.0, Dims([1, 0, -1, 0, 0, 0])),
            ("2e-2", 0.02, Dims::NONE),
            ("5e-3", 5e-3, Dims::NONE),
            ("1e-5", 1e-5, Dims::NONE),
            ("30s", 30.0, Time::DIMS),
            ("24m", 24.0, Dims([1, 0, 0, 0, 0, 0])),
        ];
        for (text, want_value, want_dims) in cases {
            let q = parse_qty(text).unwrap_or_else(|e| panic!("{text}: {e}"));
            assert!(
                (q.value - want_value).abs() <= 1e-12 * want_value.abs().max(1.0),
                "{text}: value {} != {}",
                q.value,
                want_value
            );
            assert_eq!(q.dims, *want_dims, "{text}: dims {:?}", q.dims);
        }
    }

    #[test]
    fn caret_exponents_and_negative_exponents() {
        let q = parse_qty("9.81m/s^2").expect("parses");
        assert_eq!(q.dims, Dims([1, 0, -2, 0, 0, 0]));
        assert!((q.value - 9.81).abs() < 1e-12);
        let q = parse_qty("2.5s^-1").expect("parses");
        assert_eq!(q.dims, Dims([0, 0, -1, 0, 0, 0]));
    }

    #[test]
    fn prefixes_resolve_with_longest_match_first() {
        // `min` must be minutes, not milli-"in".
        let q = parse_qty("2min").expect("parses");
        assert!((q.value - 120.0).abs() < 1e-12);
        // kN, MPa, GPa, um.
        assert!((parse_qty("3kN").unwrap().value - 3000.0).abs() < 1e-9);
        assert!((parse_qty("200MPa").unwrap().value - 2e8).abs() < 1.0);
        assert!((parse_qty("70GPa").unwrap().value - 7e10).abs() < 10.0);
        assert!((parse_qty("5um").unwrap().value - 5e-6).abs() < 1e-18);
        // kg is prefix k + gram.
        let kg = parse_qty("1.2kg").expect("parses");
        assert!((kg.value - 1.2).abs() < 1e-12);
        assert_eq!(kg.dims, Dims([0, 1, 0, 0, 0, 0]));
    }

    #[test]
    fn six_base_and_electrical_tokens_resolve_without_prefix_collisions() {
        let cases = [
            ("2mol", crate::units::moles(2.0).erase()),
            ("3V", crate::units::volts(3.0).erase()),
            ("4C", crate::units::coulombs(4.0).erase()),
            ("5Wb", crate::units::webers(5.0).erase()),
            ("6H", crate::units::henries(6.0).erase()),
            ("7Ohm", crate::units::ohms(7.0).erase()),
            ("8S", crate::units::siemens(8.0).erase()),
            ("9F", crate::units::farads(9.0).erase()),
            ("10T", crate::units::teslas(10.0).erase()),
            ("11mT", crate::units::teslas(0.011).erase()),
            ("12TW", crate::units::watts(12e12).erase()),
            ("13mS", crate::units::siemens(0.013).erase()),
            ("14mH", crate::units::henries(0.014).erase()),
            ("15mWb", crate::units::webers(0.015).erase()),
            ("16kOhm", crate::units::ohms(16e3).erase()),
            ("17uF", crate::units::farads(17e-6).erase()),
        ];
        for (text, expected) in cases {
            let parsed = parse_qty(text).unwrap_or_else(|e| panic!("{text}: {e}"));
            assert_eq!(parsed.dims, expected.dims, "{text}");
            assert!(
                (parsed.value - expected.value).abs() <= 1e-12 * expected.value.abs().max(1.0),
                "{text}: {} != {}",
                parsed.value,
                expected.value
            );
        }

        // `degC` remains the affine temperature token, while lone `C` is
        // coulomb. Whole-token matching must beat prefix interpretation.
        let celsius = parse_qty("20degC").expect("affine temperature");
        assert_eq!(celsius, crate::units::celsius(20.0).erase());
        assert_eq!(
            parse_qty("1THz").expect("tera-hertz"),
            crate::units::hertz(1e12).erase()
        );
    }

    #[test]
    fn compound_chains_apply_strict_left_to_right() {
        // kg/m/s == kg·m⁻¹·s⁻¹ under strict left-to-right division.
        let q = parse_qty("1kg/m/s").expect("parses");
        assert_eq!(
            q.dims,
            Pressure::DIMS
                .checked_plus(Dims([0, 0, 1, 0, 0, 0]))
                .expect("in-range exponents")
        );
        // density: kg/m3
        let d = parse_qty("1000kg/m3").expect("parses");
        assert_eq!(d.dims, Dims([-3, 1, 0, 0, 0, 0]));
    }

    #[test]
    fn celsius_is_affine_and_lone_only() {
        let t = parse_qty("20degC").expect("parses");
        assert_eq!(t.dims, Dims([0, 0, 0, 1, 0, 0]));
        assert!((t.value - 293.15).abs() < 1e-9);
        let e = parse_qty("20degC/s").unwrap_err();
        assert_eq!(e.kind, ParseErrorKind::AffineUnitInCompound);
        assert!(
            e.help.contains("kelvin"),
            "teaching help expected: {}",
            e.help
        );
    }

    #[test]
    fn information_units_are_refused_with_guidance() {
        let e = parse_qty("96GiB").unwrap_err();
        assert!(
            matches!(e.kind, ParseErrorKind::InformationUnit(_)),
            "{e:?}"
        );
        assert!(
            e.help.contains("fs-ir"),
            "help must point at budget syntax: {}",
            e.help
        );
    }

    #[test]
    fn unknown_units_name_the_token_and_suggest() {
        let e = parse_qty("3flurbs").unwrap_err();
        match &e.kind {
            ParseErrorKind::UnknownUnit(t) => assert_eq!(t, "flurbs"),
            k => panic!("wrong kind {k:?}"),
        }
        assert!(e.help.contains("SI unit"));
    }

    #[test]
    fn percent_parses_as_dimensionless_hundredth() {
        let q = parse_qty("15%").expect("parses");
        assert!((q.value - 0.15).abs() < 1e-15);
        assert!(q.dims.is_none());
    }

    #[test]
    fn format_then_reparse_dimensionless_round_trip() {
        // Full format→parse round-trips for arbitrary dims need unit
        // reconstruction (kg^1·m^-1 form is display-only); the dimensionless
        // path must round-trip exactly.
        for i in 0..64u32 {
            let v = f64::from(i).mul_add(0.31, -3.0);
            let s = format!("{v}");
            let q = parse_qty(&s).expect("parses");
            assert_eq!(q.value.to_bits(), v.to_bits());
            assert!(q.dims.is_none());
        }
    }
}

#[cfg(test)]
mod hardening {
    use super::{ParseBudget, ParseErrorKind, ParseResource, parse_qty, parse_qty_with_budget};
    use crate::Dims;

    fn structured_refusal(input: &str, expected_kind: &ParseErrorKind, expected_at: usize) {
        let outcome = std::panic::catch_unwind(|| parse_qty(input));
        assert!(outcome.is_ok(), "parser panicked for {input:?}");
        let error = outcome
            .expect("panic outcome checked")
            .expect_err("hostile literal must refuse");
        assert!(
            error.verifies_source(input),
            "error lost its source identity"
        );
        assert_eq!(error.input_bytes, input.len());
        assert!(error.preview.len() <= ParseBudget::DEFAULT.max_diagnostic_bytes());
        assert!(input.is_char_boundary(error.preview_start));
        assert!(input.is_char_boundary(error.preview_start + error.preview.len()));
        assert_eq!(&error.kind, expected_kind, "wrong refusal for {input:?}");
        assert_eq!(error.at, expected_at, "wrong byte offset for {input:?}");
        assert!(
            error.at <= input.len(),
            "out-of-bounds offset for {input:?}"
        );
        assert!(!error.help.is_empty(), "missing guidance for {input:?}");
    }

    /// The parser sits behind fs-ir admission and will see agent-supplied
    /// text. It must never panic — every outcome is Ok or a structured
    /// ParseError. Deterministic pseudo-random garbage battery (an LCG, not
    /// fs-rand, because UTIL crates cannot depend on L1).
    #[test]
    fn no_panic_on_garbage_input() {
        const ALPHABET: &[&str] = &[
            "0",
            "9",
            ".",
            "-",
            "+",
            "e",
            "E",
            "m",
            "g",
            "s",
            "K",
            "A",
            "N",
            "P",
            "a",
            "z",
            "µ",
            "·",
            "*",
            "/",
            "^",
            "%",
            " ",
            "\t",
            "deg",
            "C",
            "GiB",
            "\u{1F980}",
            "min",
        ];
        let mut state: u64 = 0x5EED_0001;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        for _case in 0..20_000 {
            let len = next() % 24;
            let mut s = String::new();
            for _ in 0..len {
                s.push_str(ALPHABET[next() % ALPHABET.len()]);
            }
            // Must not panic; the Result content is irrelevant here.
            let _ = parse_qty(&s);
        }
    }

    /// Boundary battery: empty, lone signs, huge exponents, deep unit chains.
    #[test]
    fn boundary_inputs_return_structured_errors() {
        for bad in [
            "", " ", "+", "-", ".", "e5", "1m^", "1m^-", "1s^999", "3", "1//s", "1m**s",
        ] {
            match parse_qty(bad) {
                Ok(q) => assert!(
                    q.dims.is_none() || !bad.is_empty(),
                    "unexpected acceptance of {bad:?} -> {q:?}"
                ),
                Err(e) => {
                    assert!(!e.help.is_empty(), "error for {bad:?} must carry guidance");
                }
            }
        }
        // A moderate chain parses fine; a runaway 200-deep chain must return
        // a STRUCTURED error, never panic (this exact case found the i8
        // overflow this battery exists to catch).
        let moderate = format!("1{}", "m/s*".repeat(20) + "m");
        assert!(parse_qty(&moderate).is_ok(), "moderate chain should parse");
        let runaway = format!("1{}", "m/s*".repeat(200) + "m");
        let e = parse_qty(&runaway).expect_err("runaway chain must be rejected");
        assert_eq!(e.kind, ParseErrorKind::BadExponent);
        assert!(
            e.help.contains("runaway"),
            "teaching help expected: {}",
            e.help
        );
    }

    #[test]
    fn exponent_cap_is_checked_wide_before_narrowing() {
        assert_eq!(
            parse_qty("1m^60").expect("positive cap is valid").dims,
            Dims([60, 0, 0, 0, 0, 0])
        );
        assert_eq!(
            parse_qty("1m^-60").expect("negative cap is valid").dims,
            Dims([-60, 0, 0, 0, 0, 0])
        );
        assert_eq!(
            parse_qty("1Pa^30")
                .expect("derived dimension at cap is valid")
                .dims,
            Dims([-30, 30, -60, 0, 0, 0])
        );

        for (input, at) in [
            ("1m^61", 2),
            ("1m^-61", 2),
            ("1m^-128", 2),
            ("1Pa^31", 3),
            ("1Pa^127", 3),
            ("  1m^-128", 4),
        ] {
            structured_refusal(input, &ParseErrorKind::BadExponent, at);
        }
    }

    #[test]
    fn repeated_factor_chains_enforce_every_intermediate_cap() {
        let product_at_cap = format!("1{}", vec!["m"; 60].join("*"));
        let product_over_cap = format!("1{}", vec!["m"; 61].join("*"));
        let product_then_cancel = format!("{product_over_cap}/m");
        let quotient_at_cap = format!("1rad{}", "/m".repeat(60));
        let quotient_over_cap = format!("1rad{}", "/m".repeat(61));

        assert_eq!(
            parse_qty(&product_at_cap)
                .expect("positive repeated chain at cap")
                .dims,
            Dims([60, 0, 0, 0, 0, 0])
        );
        assert_eq!(
            parse_qty(&quotient_at_cap)
                .expect("negative repeated chain at cap")
                .dims,
            Dims([-60, 0, 0, 0, 0, 0])
        );
        for input in [&product_over_cap, &product_then_cancel, &quotient_over_cap] {
            let outcome = std::panic::catch_unwind(|| parse_qty(input));
            assert!(outcome.is_ok(), "parser panicked for repeated chain");
            let error = outcome
                .expect("panic outcome checked")
                .expect_err("cap+1 intermediate must refuse");
            assert_eq!(error.kind, ParseErrorKind::BadExponent);
            assert!(error.help.contains("±60"), "unexpected guidance: {error}");
            assert!(error.at <= input.len());
        }
    }

    #[test]
    fn every_nonfinite_quantity_path_refuses() {
        for (input, at) in [
            ("1e999", 0),
            ("1e999m", 0),
            ("1e-999", 0),
            ("-1e-999", 0),
            ("  1e999degC", 2),
            ("1e308TW", 5),
            ("1Trad^26", 5),
            ("1THz^-30", 4),
            ("1e-308pHz^2", 9),
            ("1e-308rad/THz^25", 13),
            ("1rad/THz^-30", 8),
        ] {
            structured_refusal(input, &ParseErrorKind::NonFiniteValue, at);
        }
    }

    #[test]
    fn textual_zero_remains_valid_with_large_exponents() {
        for input in ["0", "-0.0", "0e999", "-0e999"] {
            let quantity = parse_qty(input).unwrap_or_else(|error| panic!("{input}: {error}"));
            assert_eq!(quantity.value.abs().to_bits(), 0.0_f64.to_bits(), "{input}");
            assert!(quantity.dims.is_none(), "{input}");
        }
    }

    #[test]
    fn whitespace_is_counted_in_error_offsets() {
        let error =
            parse_qty("  1m *   flurbs").expect_err("unknown unit after whitespace must refuse");
        assert!(matches!(error.kind, ParseErrorKind::UnknownUnit(_)));
        assert_eq!(error.at, 9);
    }

    #[test]
    fn byte_factor_and_token_budgets_hold_at_exact_boundaries() {
        let byte_cap = ParseBudget::DEFAULT.max_input_bytes();
        let at_cap = format!("1m{}", " ".repeat(byte_cap - 2));
        assert_eq!(at_cap.len(), byte_cap);
        assert!(parse_qty(&at_cap).is_ok(), "exact byte boundary must admit");
        let over_cap = format!("{at_cap} ");
        let error = parse_qty(&over_cap).expect_err("byte boundary + 1 must refuse");
        assert!(matches!(
            error.kind,
            ParseErrorKind::BudgetExceeded {
                resource: ParseResource::InputBytes,
                limit,
                observed_at_least,
            } if limit == byte_cap && observed_at_least == byte_cap + 1
        ));
        assert_eq!(
            error.source_hash, None,
            "oversized input must not be hash-scanned"
        );
        assert!(error.preview.len() <= ParseBudget::DEFAULT.max_diagnostic_bytes());

        let factor_cap = ParseBudget::DEFAULT.max_factors();
        let factors_at_cap = format!("1{}", vec!["rad"; factor_cap].join("*"));
        assert!(
            parse_qty(&factors_at_cap).is_ok(),
            "exact factor boundary must admit"
        );
        let factors_over_cap = format!("{factors_at_cap}*rad");
        let error = parse_qty(&factors_over_cap).expect_err("factor boundary + 1 must refuse");
        assert!(matches!(
            error.kind,
            ParseErrorKind::BudgetExceeded {
                resource: ParseResource::Factors,
                limit,
                observed_at_least,
            } if limit == factor_cap && observed_at_least == factor_cap + 1
        ));

        let token_cap = ParseBudget::DEFAULT.max_token_bytes();
        let unknown_at_cap = format!("1{}", "x".repeat(token_cap));
        let error = parse_qty(&unknown_at_cap).expect_err("bounded unknown token must refuse");
        assert!(
            matches!(error.kind, ParseErrorKind::UnknownUnit(ref token) if token.len() == token_cap)
        );
        let unknown_over_cap = format!("{unknown_at_cap}x");
        let error = parse_qty(&unknown_over_cap).expect_err("token boundary + 1 must refuse");
        assert!(matches!(
            error.kind,
            ParseErrorKind::BudgetExceeded {
                resource: ParseResource::TokenBytes,
                limit,
                observed_at_least,
            } if limit == token_cap && observed_at_least > token_cap
        ));
    }

    #[test]
    fn exponent_and_zero_budgets_are_explicit() {
        let exponent_budget = ParseBudget::new(64, 1, 3, 32);
        assert!(parse_qty_with_budget("1m^60", exponent_budget).is_ok());
        let error = parse_qty_with_budget("1m^060", exponent_budget)
            .expect_err("exponent token boundary + 1 must refuse");
        assert!(matches!(
            error.kind,
            ParseErrorKind::BudgetExceeded {
                resource: ParseResource::TokenBytes,
                limit: 3,
                observed_at_least: 4,
            }
        ));

        let zero = ParseBudget::new(0, 0, 0, 0);
        let empty = parse_qty_with_budget("", zero).expect_err("empty input is not a quantity");
        assert_eq!(empty.preview, "");
        assert!(empty.verifies_source(""));
        let nonempty = parse_qty_with_budget("1", zero)
            .expect_err("zero byte budget must reject nonempty input");
        assert!(matches!(
            nonempty.kind,
            ParseErrorKind::BudgetExceeded {
                resource: ParseResource::InputBytes,
                limit: 0,
                observed_at_least: 1,
            }
        ));
        assert_eq!(nonempty.preview, "");
        assert_eq!(nonempty.source_hash, None);
    }

    #[test]
    fn multibyte_previews_keep_exact_offsets_and_utf8_boundaries() {
        let input = "  1m * λλλλλλλλ";
        let budget = ParseBudget::new(128, 8, 5, 9);
        let first_lambda = input.find('λ').expect("fixture contains lambda");
        let error = parse_qty_with_budget(input, budget)
            .expect_err("multi-byte token over budget must refuse");
        assert_eq!(error.at, first_lambda);
        assert!(matches!(
            &error.kind,
            ParseErrorKind::BudgetExceeded {
                resource: ParseResource::TokenBytes,
                limit: 5,
                observed_at_least: 6,
            }
        ));
        assert!(error.verifies_source(input));
        assert!(error.preview.len() <= 9);
        assert!(input.is_char_boundary(error.preview_start));
        assert!(input.is_char_boundary(error.preview_start + error.preview.len()));
        assert!(
            error.preview_start <= error.at
                && error.at <= error.preview_start + error.preview.len()
        );
        assert_eq!(error, parse_qty_with_budget(input, budget).unwrap_err());
    }
}
