//! Wasm32-clean FrankenScript syntax AST (`ast`), isomorphic s-expression
//! (`sexpr`) and JSON (`json`) codecs, and structured [`IrError`].
//!
//! Moved isomorphically out of `fs-ir` (`frankensim-rc-1109-root-3gc12.11.9`,
//! Cut R2a) so `fs-ir`, `fs-cli`, and `fs-runtime` share one parser with zero
//! `fs-ledger`/SQLite dependency on `wasm32-unknown-unknown`.

#![forbid(unsafe_code)]

pub mod ast;
pub mod codec;
pub mod json;
pub mod sexpr;

pub use ast::{CountUnit, CountValue, Node, NodeKind, Span};
pub use codec::{
    AdmittedAssembly, AssemblyAnalysisMode, AssemblyRequirement, BOX_AXES, CodecRefusal,
    DeclaredEnvelopes, GeomEnvelope, MAX_PACKS, PackCatalog, PackEntry, PackKind, ProbeStat,
    TargetEnvelope, admit_assembly_source,
};

use core::fmt;
use fs_qty::Dims;

/// Canonical SI unit labels indexed by base dimension exponents `[L, M, T, Theta, I, N]`.
pub const UNIT_DIMS: &[(&str, Dims)] = &[
    ("m", Dims([1, 0, 0, 0, 0, 0])),
    ("m^3", Dims([3, 0, 0, 0, 0, 0])),
    ("m/s", Dims([1, 0, -1, 0, 0, 0])),
    ("s", Dims([0, 0, 1, 0, 0, 0])),
    ("K", Dims([0, 0, 0, 1, 0, 0])),
    ("W", Dims([2, 1, -3, 0, 0, 0])),
    ("rad/s", Dims([0, 0, -1, 0, 0, 0])),
    ("1", Dims::NONE),
    ("m/N", Dims([0, -1, 2, 0, 0, 0])),
    ("m/(N*s)", Dims([0, -1, 1, 0, 0, 0])),
    ("K/W", Dims([-2, -1, 3, 1, 0, 0])),
    ("kg", Dims([0, 1, 0, 0, 0, 0])),
    ("N", Dims([1, 1, -2, 0, 0, 0])),
    ("N*m", Dims([2, 1, -2, 0, 0, 0])),
    ("Pa", Dims([-1, 1, -2, 0, 0, 0])),
    ("V", Dims([2, 1, -3, 0, -1, 0])),
    ("A", Dims([0, 0, 0, 0, 1, 0])),
];

/// Map SI base dimensions to their canonical display unit label, when known.
#[must_use]
pub fn dims_to_unit(d: Dims) -> Option<&'static str> {
    UNIT_DIMS.iter().find(|(_, x)| *x == d).map(|(u, _)| *u)
}

/// What went wrong while reading a program.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IrErrorKind {
    /// A character that cannot start a form.
    UnexpectedChar,
    /// Input ended mid-form.
    UnexpectedEnd,
    /// Content after the single top-level form.
    TrailingInput,
    /// A `(` without its `)`.
    UnclosedParen,
    /// A `"` without its closing `"`.
    UnclosedString,
    /// An unknown escape sequence.
    BadEscape,
    /// A bare `:`.
    BadKeyword,
    /// A malformed numeric literal.
    BadNumber,
    /// A malformed `0x…` seed.
    BadSeed,
    /// A numeric token that is not an int, float, quantity, or count.
    BadQuantity,
    /// Nesting beyond the depth cap.
    TooDeep,
    /// Malformed JSON structure.
    JsonSyntax,
    /// An unknown atom tag in the JSON mapping.
    JsonUnknownTag,
    /// A tagged literal that does not match its tag.
    JsonTagMismatch,
    /// A versioned artifact uses unsupported language semantics.
    UnsupportedVersion,
    /// Persisted versioned bytes are semantically valid but not canonical.
    NonCanonical,
    /// Expected a study form.
    NotAStudy,
    /// A recognized clause with the wrong shape.
    MalformedClause,
}

impl IrErrorKind {
    /// Stable machine-readable code.
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            IrErrorKind::UnexpectedChar => "IrUnexpectedChar",
            IrErrorKind::UnexpectedEnd => "IrUnexpectedEnd",
            IrErrorKind::TrailingInput => "IrTrailingInput",
            IrErrorKind::UnclosedParen => "IrUnclosedParen",
            IrErrorKind::UnclosedString => "IrUnclosedString",
            IrErrorKind::BadEscape => "IrBadEscape",
            IrErrorKind::BadKeyword => "IrBadKeyword",
            IrErrorKind::BadNumber => "IrBadNumber",
            IrErrorKind::BadSeed => "IrBadSeed",
            IrErrorKind::BadQuantity => "IrBadQuantity",
            IrErrorKind::TooDeep => "IrTooDeep",
            IrErrorKind::JsonSyntax => "IrJsonSyntax",
            IrErrorKind::JsonUnknownTag => "IrJsonUnknownTag",
            IrErrorKind::JsonTagMismatch => "IrJsonTagMismatch",
            IrErrorKind::UnsupportedVersion => "IrUnsupportedVersion",
            IrErrorKind::NonCanonical => "IrNonCanonical",
            IrErrorKind::NotAStudy => "IrNotAStudy",
            IrErrorKind::MalformedClause => "IrMalformedClause",
        }
    }
}

/// A structured parse/recognition failure: span + diagnosis + fix hint
/// (Decalogue P10 — a refusal that teaches).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IrError {
    /// The offending byte range in the source.
    pub span: Span,
    /// The failure class.
    pub kind: IrErrorKind,
    /// What is wrong.
    pub detail: String,
    /// How to fix it.
    pub hint: String,
}

impl fmt::Display for IrError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} at bytes {}..{}: {}; fix: {}",
            self.kind.code(),
            self.span.start,
            self.span.end,
            self.detail,
            self.hint
        )
    }
}

impl std::error::Error for IrError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sexpr_and_json_round_trip_and_dims_to_unit_match() {
        let ast = sexpr::parse("(assembly \"demo\" (version 1) (seed 0x22))").unwrap();
        let json = json::print_checked(&ast).unwrap();
        let back = json::parse(&json).unwrap();
        assert_eq!(
            sexpr::print_checked(&back).unwrap(),
            sexpr::print_checked(&ast).unwrap()
        );
        assert_eq!(dims_to_unit(Dims([1, 1, -2, 0, 0, 0])), Some("N"));
        assert_eq!(dims_to_unit(Dims([0, 0, -1, 0, 0, 0])), Some("rad/s"));
        assert_eq!(dims_to_unit(Dims([1, 2, 3, 0, 0, 0])), None);
        assert!(PackCatalog::from_raw_bytes(&[b"BADMAGIC"]).is_err());
        let steel: &[u8] = include_bytes!(
            "../../../examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
        );
        let dedup = PackCatalog::from_raw_bytes(&[steel, steel]).unwrap();
        assert_eq!(dedup.materials().len(), 1);
    }
}
