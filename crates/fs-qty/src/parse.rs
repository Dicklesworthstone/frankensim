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
//! - `mol`/`cd` are outside FrankenSim's 5-dimension vector (no-claim).

use crate::{Dims, QtyAny};
use core::fmt;

/// Where and why parsing failed, with a suggested fix (errors-as-guidance).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    /// The full input.
    pub input: String,
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
    /// Exponent didn't parse or overflowed i8.
    BadExponent,
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "cannot parse quantity {:?} at byte {}: {:?}; {}",
            self.input, self.at, self.kind, self.help
        )
    }
}

impl core::error::Error for ParseError {}

/// A named unit: symbol → (scale-to-SI, dimension). `degC` handled specially.
struct Unit {
    symbol: &'static str,
    scale: f64,
    dims: Dims,
}

const D_NONE: Dims = Dims([0, 0, 0, 0, 0]);
const D_M: Dims = Dims([1, 0, 0, 0, 0]);
const D_KG: Dims = Dims([0, 1, 0, 0, 0]);
const D_S: Dims = Dims([0, 0, 1, 0, 0]);
const D_K: Dims = Dims([0, 0, 0, 1, 0]);
const D_A: Dims = Dims([0, 0, 0, 0, 1]);
const D_N: Dims = Dims([1, 1, -2, 0, 0]);
const D_PA: Dims = Dims([-1, 1, -2, 0, 0]);
const D_J: Dims = Dims([2, 1, -2, 0, 0]);
const D_W: Dims = Dims([2, 1, -3, 0, 0]);
const D_HZ: Dims = Dims([0, 0, -1, 0, 0]);
const D_M3: Dims = Dims([3, 0, 0, 0, 0]);

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

fn err(input: &str, at: usize, kind: ParseErrorKind, help: &str) -> ParseError {
    ParseError {
        input: input.to_string(),
        at,
        kind,
        help: help.to_string(),
    }
}

/// Resolve one unit token (no exponent) to (scale, dims).
fn resolve_token(input: &str, at: usize, tok: &str) -> Result<(f64, Dims), ParseError> {
    // Whole-token named unit wins (so `min` is minutes, not milli-inches).
    if let Some(u) = UNITS.iter().find(|u| u.symbol == tok) {
        return Ok((u.scale, u.dims));
    }
    // Information units get a dedicated refusal.
    if INFORMATION_UNITS.iter().any(|s| tok.ends_with(s)) {
        return Err(err(
            input,
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
        input,
        at,
        ParseErrorKind::UnknownUnit(tok.to_string()),
        "expected an SI unit like m, kg, s, K, A, N, Pa, J, W, Hz, L, min, h, rad, deg, % \
         with an optional prefix (p n u m c d k M G T)",
    ))
}

/// Scan the leading number of a quantity literal; returns its byte length.
fn scan_number(s: &str) -> usize {
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
        if c == 'e' || c == 'E' {
            end += 1; // consume the sign/digit that justified accepting 'e'
        }
        end += 1;
    }
    end
}

/// Parse an optional exponent (`^`? `-`? digits) at the head of `rest`.
/// Returns `(exponent, bytes_consumed)`; the default exponent is 1.
fn parse_exponent(input: &str, pos: usize, rest: &str) -> Result<(i8, usize), ParseError> {
    let mut r = rest;
    let mut consumed = 0;
    if let Some(stripped) = r.strip_prefix('^') {
        r = stripped;
        consumed += 1;
    }
    let neg = if let Some(stripped) = r.strip_prefix('-') {
        r = stripped;
        consumed += 1;
        true
    } else {
        false
    };
    let dig_len = r
        .char_indices()
        .find(|(_, c)| !c.is_ascii_digit())
        .map_or(r.len(), |(i, _)| i);
    if dig_len == 0 {
        if consumed > 0 {
            return Err(err(
                input,
                pos,
                ParseErrorKind::BadExponent,
                "dangling ^ or - without digits",
            ));
        }
        return Ok((1, 0));
    }
    let parsed: i32 = r[..dig_len].parse().map_err(|_| {
        err(
            input,
            pos,
            ParseErrorKind::BadExponent,
            "exponent must be a small integer",
        )
    })?;
    let signed = if neg { -parsed } else { parsed };
    let exp = i8::try_from(signed).map_err(|_| {
        err(
            input,
            pos,
            ParseErrorKind::BadExponent,
            "exponent magnitude too large",
        )
    })?;
    Ok((exp, consumed + dig_len))
}

/// Parse a quantity literal into a [`QtyAny`].
///
/// # Errors
/// Returns a [`ParseError`] with position, kind, and a suggested fix.
pub fn parse_qty(input: &str) -> Result<QtyAny, ParseError> {
    let s = input.trim();
    let base = input.len() - input.trim_start().len();

    let end = scan_number(s);
    let num: f64 = s[..end].parse().map_err(|_| {
        err(
            input,
            base,
            ParseErrorKind::MissingNumber,
            "a quantity starts with a number, e.g. 0.12Pa*s",
        )
    })?;
    let mut rest = s[end..].trim_start();
    let mut pos = base + end;

    // --- bare number: dimensionless ---
    if rest.is_empty() {
        return Ok(QtyAny::dimensionless(num));
    }

    // --- special-case lone affine unit degC ---
    if rest == "degC" {
        return Ok(QtyAny::new(num + 273.15, D_K));
    }

    // --- unit expression, strict left-to-right ---
    let mut value = num;
    let mut dims = Dims::NONE;
    let mut divide = false;
    loop {
        // token = leading unit letters/µ/%
        let tok_len = rest
            .char_indices()
            .find(|(_, c)| !(c.is_alphabetic() || *c == 'µ' || *c == '%'))
            .map_or(rest.len(), |(i, _)| i);
        if tok_len == 0 {
            return Err(err(
                input,
                pos,
                ParseErrorKind::TrailingInput,
                "expected a unit symbol here",
            ));
        }
        let tok = &rest[..tok_len];
        if tok.contains("degC") {
            return Err(err(
                input,
                pos,
                ParseErrorKind::AffineUnitInCompound,
                "degC is affine and only legal alone (e.g. \"20degC\"); temperature \
                 differences and rates are kelvin — write K or K/s",
            ));
        }
        let (scale, tok_dims) = resolve_token(input, pos, tok)?;
        rest = &rest[tok_len..];
        pos += tok_len;

        // optional exponent: '^'? '-'? digits
        let (exp, consumed) = parse_exponent(input, pos, rest)?;
        rest = &rest[consumed..];
        pos += consumed;

        // apply factor
        let factor_scale = scale.powi(i32::from(exp));
        let factor_dims = tok_dims.times(exp);
        if divide {
            value /= factor_scale;
            dims = dims.minus(factor_dims);
        } else {
            value *= factor_scale;
            dims = dims.plus(factor_dims);
        }
        // Semantic sanity: physically meaningful exponents are tiny. Reject
        // accumulations beyond ±60 with a structured error well before the
        // (saturating) i8 arithmetic could mask an adversarial chain.
        if dims.0.iter().any(|&d| d.abs() > 60) {
            return Err(err(
                input,
                pos,
                ParseErrorKind::BadExponent,
                "accumulated unit exponent beyond ±60 is not physically meaningful; \
                 check for a runaway unit chain",
            ));
        }

        // separator or end
        rest = rest.trim_start();
        match rest.chars().next() {
            None => return Ok(QtyAny::new(value, dims)),
            Some('*' | '·') => {
                divide = false;
                let c = rest.chars().next().expect("just matched");
                rest = rest[c.len_utf8()..].trim_start();
                pos += c.len_utf8();
            }
            Some('/') => {
                divide = true;
                rest = rest[1..].trim_start();
                pos += 1;
            }
            Some(_) => {
                return Err(err(
                    input,
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
            ("3mm", 3e-3, Dims([1, 0, 0, 0, 0])),
            ("12mm", 12e-3, Dims([1, 0, 0, 0, 0])),
            ("0deg", 0.0, Dims::NONE),
            ("65deg", 65.0 * core::f64::consts::PI / 180.0, Dims::NONE),
            ("3s", 3.0, Time::DIMS),
            ("2h", 7200.0, Time::DIMS),
            ("0.03m2/s3", 0.03, Dims([2, 0, -3, 0, 0])),
            ("15rad/s", 15.0, Dims([0, 0, -1, 0, 0])),
            ("8m/s", 8.0, Dims([1, 0, -1, 0, 0])),
            ("2e-2", 0.02, Dims::NONE),
            ("5e-3", 5e-3, Dims::NONE),
            ("1e-5", 1e-5, Dims::NONE),
            ("30s", 30.0, Time::DIMS),
            ("24m", 24.0, Dims([1, 0, 0, 0, 0])),
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
        assert_eq!(q.dims, Dims([1, 0, -2, 0, 0]));
        assert!((q.value - 9.81).abs() < 1e-12);
        let q = parse_qty("2.5s^-1").expect("parses");
        assert_eq!(q.dims, Dims([0, 0, -1, 0, 0]));
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
        assert_eq!(kg.dims, Dims([0, 1, 0, 0, 0]));
    }

    #[test]
    fn compound_chains_apply_strict_left_to_right() {
        // kg/m/s == kg·m⁻¹·s⁻¹ under strict left-to-right division.
        let q = parse_qty("1kg/m/s").expect("parses");
        assert_eq!(q.dims, Pressure::DIMS.plus(Dims([0, 0, 1, 0, 0])));
        // density: kg/m3
        let d = parse_qty("1000kg/m3").expect("parses");
        assert_eq!(d.dims, Dims([-3, 1, 0, 0, 0]));
    }

    #[test]
    fn celsius_is_affine_and_lone_only() {
        let t = parse_qty("20degC").expect("parses");
        assert_eq!(t.dims, Dims([0, 0, 0, 1, 0]));
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
    use super::{ParseErrorKind, parse_qty};

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
}
