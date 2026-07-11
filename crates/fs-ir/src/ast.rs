//! The FrankenScript typed AST (plan §11.1): the ONE tree both concrete
//! syntaxes parse into. Atoms are the system's real nouns — dimensioned
//! quantities (fs-qty dims), seeds, counts, strings, symbols, keywords —
//! and every node carries a source span so structured errors point at
//! exact locations in agent-submitted programs.
//!
//! Equality discipline: spans are provenance, not meaning. Use
//! [`Node::same_shape`] for semantic comparison (the isomorphism property
//! s-expr ↔ JSON ↔ AST is stated in terms of it); derived `PartialEq`
//! includes spans and is for same-source comparisons only.

use fs_qty::Dims;

/// Byte range in the source text a node came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Span {
    /// Start byte offset (inclusive).
    pub start: usize,
    /// End byte offset (exclusive).
    pub end: usize,
}

impl Span {
    /// A span covering `start..end`.
    #[must_use]
    pub fn new(start: usize, end: usize) -> Span {
        Span { start, end }
    }
}

/// Byte-count units the IR recognizes for capability/budget grants
/// (information units are deliberately OUTSIDE fs-qty's SI domain).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountUnit {
    /// Bytes.
    B,
    /// 2¹⁰ bytes.
    KiB,
    /// 2²⁰ bytes.
    MiB,
    /// 2³⁰ bytes.
    GiB,
    /// Processor cores.
    Cores,
}

/// A bounded exact decimal used by count literals containing a decimal
/// point or exponent. Its value is `(-1)^negative * significand * 10^exponent10`.
/// Keeping this decimal form exact prevents source text such as
/// `0.99999999999999999B` from rounding into a different resource claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DecimalCount {
    negative: bool,
    significand: u128,
    exponent10: i32,
}

impl DecimalCount {
    pub(crate) fn parse(text: &str) -> Option<Self> {
        let bytes = text.as_bytes();
        let mut index = 0;
        let negative = match bytes.first() {
            Some(b'-') => {
                index = 1;
                true
            }
            Some(b'+') => {
                index = 1;
                false
            }
            _ => false,
        };

        let mut significand = 0_u128;
        let mut fractional_digits = 0_i32;
        let mut saw_digit = false;
        let mut saw_decimal = false;
        while let Some(&byte) = bytes.get(index) {
            match byte {
                b'0'..=b'9' => {
                    saw_digit = true;
                    significand = significand
                        .checked_mul(10)?
                        .checked_add(u128::from(byte - b'0'))?;
                    if saw_decimal {
                        fractional_digits = fractional_digits.checked_add(1)?;
                    }
                    index += 1;
                }
                b'.' if !saw_decimal => {
                    saw_decimal = true;
                    index += 1;
                }
                _ => break,
            }
        }
        if !saw_digit {
            return None;
        }

        let mut written_exponent = 0_i32;
        let mut saw_exponent = false;
        if matches!(bytes.get(index), Some(b'e' | b'E')) {
            saw_exponent = true;
            index += 1;
            let exponent_start = index;
            if matches!(bytes.get(index), Some(b'+' | b'-')) {
                index += 1;
            }
            let digits_start = index;
            while bytes.get(index).is_some_and(u8::is_ascii_digit) {
                index += 1;
            }
            if index == digits_start {
                return None;
            }
            written_exponent = text[exponent_start..index].parse().ok()?;
        }
        if index != bytes.len()
            || (!saw_decimal
                && !saw_exponent
                && text.as_bytes()[0] != b'+'
                && text.as_bytes()[0] != b'-')
        {
            return None;
        }

        let mut exponent10 = written_exponent.checked_sub(fractional_digits)?;
        if significand == 0 {
            return Some(Self {
                negative: false,
                significand: 0,
                exponent10: 0,
            });
        }
        while significand.is_multiple_of(10) {
            significand /= 10;
            exponent10 = exponent10.checked_add(1)?;
        }
        Some(Self {
            negative,
            significand,
            exponent10,
        })
    }

    fn checked_integral_bytes(self, shift: u32) -> Option<u64> {
        if self.negative {
            return None;
        }
        let integral = if self.exponent10 >= 0 {
            let exponent = u32::try_from(self.exponent10).ok()?;
            let scaled = self
                .significand
                .checked_mul(10_u128.checked_pow(exponent)?)?;
            scaled
                .checked_shl(shift)
                .filter(|value| *value >> shift == scaled)?
        } else {
            // Divide the exact decimal denominator 2^n*5^n before any
            // multiplication. This both proves integrality and avoids a
            // transient u128 overflow for values whose reduced byte count is
            // small.
            let denominator_power = self.exponent10.unsigned_abs();
            // A nonzero u128 cannot contain more than 55 factors of five;
            // cap before the loop so hostile exponents remain O(1).
            if denominator_power > 55 {
                return None;
            }
            let mut reduced = self.significand;
            for _ in 0..denominator_power {
                if !reduced.is_multiple_of(5) {
                    return None;
                }
                reduced /= 5;
            }
            if denominator_power > shift {
                let remaining_twos = denominator_power - shift;
                if reduced.trailing_zeros() < remaining_twos {
                    return None;
                }
                reduced >> remaining_twos
            } else {
                let extra_twos = shift - denominator_power;
                reduced
                    .checked_shl(extra_twos)
                    .filter(|value| *value >> extra_twos == reduced)?
            }
        };
        u64::try_from(integral).ok()
    }

    pub(crate) fn canonical(self) -> String {
        format!(
            "{}{}e{}",
            if self.negative { "-" } else { "" },
            self.significand,
            self.exponent10
        )
    }

    fn approx_f64(self) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        // det-ok: documented APPROXIMATE view; identity/enforcement use the exact significand/exponent form, so a 1-ULP profile divergence cannot reach any canonical hash
        let magnitude = (self.significand as f64) * 10_f64.powi(self.exponent10);
        if self.negative { -magnitude } else { magnitude }
    }
}

/// Exact count magnitude (bead gp3.20). Bare integer literals retain a
/// `u128`; decimal/exponent forms retain a bounded exact decimal. Neither
/// identity nor resource enforcement projects through binary floating point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountValue {
    /// Written as a bare integer literal: exact.
    Exact(u128),
    /// Written with a decimal point, exponent, or explicit sign. The name
    /// distinguishes syntax classes; the represented value may be integral.
    Fractional(DecimalCount),
}

impl CountValue {
    /// Magnitude for REPORTING/telemetry only — may round above 2^53;
    /// identity and enforcement never go through this.
    #[must_use]
    pub fn approx_f64(self) -> f64 {
        match self {
            #[allow(clippy::cast_precision_loss)]
            CountValue::Exact(v) => v as f64,
            CountValue::Fractional(decimal) => decimal.approx_f64(),
        }
    }

    /// Exact non-negative integral magnitude in the written unit. This is the
    /// authority-safe path for integral count resources such as concurrent
    /// cores; it refuses negative, fractional, and overflowing values.
    #[must_use]
    pub fn integral_count(self) -> Option<u64> {
        match self {
            CountValue::Exact(value) => u64::try_from(value).ok(),
            CountValue::Fractional(decimal) => decimal.checked_integral_bytes(0),
        }
    }

    /// Exact integral BYTES under `unit`, or `None` when the question
    /// has no exact answer: Cores (not bytes), u64/u128 overflow, a negative
    /// value, or a decimal that does not scale to a whole byte.
    #[must_use]
    pub fn integral_bytes(self, unit: CountUnit) -> Option<u64> {
        let shift = unit.byte_shift()?;
        match self {
            CountValue::Exact(v) => {
                let scaled = v.checked_shl(shift).filter(|s| s >> shift == v)?;
                u64::try_from(scaled).ok()
            }
            CountValue::Fractional(decimal) => decimal.checked_integral_bytes(shift),
        }
    }
}

impl CountUnit {
    /// log2 of the byte factor; `None` for non-byte counts (Cores).
    #[must_use]
    pub fn byte_shift(self) -> Option<u32> {
        match self {
            CountUnit::B => Some(0),
            CountUnit::KiB => Some(10),
            CountUnit::MiB => Some(20),
            CountUnit::GiB => Some(30),
            CountUnit::Cores => None,
        }
    }

    /// Canonical suffix text.
    #[must_use]
    pub fn suffix(self) -> &'static str {
        match self {
            CountUnit::B => "B",
            CountUnit::KiB => "KiB",
            CountUnit::MiB => "MiB",
            CountUnit::GiB => "GiB",
            CountUnit::Cores => "cores",
        }
    }

    /// Parse a suffix.
    #[must_use]
    pub fn from_suffix(s: &str) -> Option<CountUnit> {
        match s {
            "B" => Some(CountUnit::B),
            "KiB" => Some(CountUnit::KiB),
            "MiB" => Some(CountUnit::MiB),
            "GiB" => Some(CountUnit::GiB),
            "cores" => Some(CountUnit::Cores),
            _ => None,
        }
    }
}

/// A node's meaning.
#[derive(Debug, Clone, PartialEq)]
pub enum NodeKind {
    /// Integer literal.
    Int(i64),
    /// Float literal.
    Float(f64),
    /// Dimensioned quantity: SI value + dims (fs-qty), plus the ORIGINAL
    /// literal text. fs-qty normalizes to SI (65deg → 1.134… rad), so the
    /// literal must be preserved verbatim for lossless printing; semantic
    /// equality uses (value, dims) only.
    Qty {
        /// Value in SI base units.
        value: f64,
        /// The five base dimensions.
        dims: Dims,
        /// The literal as written (e.g. `"65deg"`), printed verbatim.
        text: String,
    },
    /// Non-SI count (memory grants, core counts). Integer literals are
    /// EXACT (gp3.20); fractional forms are explicitly fractional.
    Count {
        /// Magnitude in the written unit.
        value: CountValue,
        /// The unit.
        unit: CountUnit,
    },
    /// Seed literal (`0x...`): the Five Explicits' seed pillar.
    Seed(u64),
    /// String literal.
    Str(String),
    /// Symbol (operator names, identifiers).
    Symbol(String),
    /// Keyword (`:name`).
    Keyword(String),
    /// A form: `(head args...)`.
    List(Vec<Node>),
}

/// One AST node: meaning + provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    /// The node's meaning.
    pub kind: NodeKind,
    /// Where it came from in the source text.
    pub span: Span,
}

impl Node {
    /// A node with an empty span (synthesized nodes, e.g. verb lowering).
    #[must_use]
    pub fn synthetic(kind: NodeKind) -> Node {
        Node {
            kind,
            span: Span::default(),
        }
    }

    /// Semantic equality: identical trees ignoring spans and Qty unit
    /// presentation (floats compare by bits; NaN never occurs — both
    /// parsers reject non-finite literals).
    #[must_use]
    pub fn same_shape(&self, other: &Node) -> bool {
        match (&self.kind, &other.kind) {
            (NodeKind::Int(a), NodeKind::Int(b)) => a == b,
            (NodeKind::Float(a), NodeKind::Float(b)) => a.to_bits() == b.to_bits(),
            (
                NodeKind::Qty {
                    value: va,
                    dims: da,
                    ..
                },
                NodeKind::Qty {
                    value: vb,
                    dims: db,
                    ..
                },
            ) => va.to_bits() == vb.to_bits() && da == db,
            (
                NodeKind::Count {
                    value: va,
                    unit: ua,
                },
                NodeKind::Count {
                    value: vb,
                    unit: ub,
                },
            ) => {
                // Mixed syntax classes remain distinct on purpose: `2B`
                // and `2.0B` are different written claims even though both
                // enforce to two bytes.
                va == vb && ua == ub
            }
            (NodeKind::Seed(a), NodeKind::Seed(b)) => a == b,
            (NodeKind::Str(a), NodeKind::Str(b))
            | (NodeKind::Symbol(a), NodeKind::Symbol(b))
            | (NodeKind::Keyword(a), NodeKind::Keyword(b)) => a == b,
            (NodeKind::List(a), NodeKind::List(b)) => {
                a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.same_shape(y))
            }
            _ => false,
        }
    }

    /// The head symbol of a form, if this is a `(symbol ...)` list.
    #[must_use]
    pub fn head(&self) -> Option<&str> {
        if let NodeKind::List(items) = &self.kind
            && let Some(first) = items.first()
            && let NodeKind::Symbol(s) = &first.kind
        {
            return Some(s);
        }
        None
    }

    /// The elements of a list node.
    #[must_use]
    pub fn items(&self) -> Option<&[Node]> {
        match &self.kind {
            NodeKind::List(items) => Some(items),
            _ => None,
        }
    }
}
