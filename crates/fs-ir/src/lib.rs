//! fs-ir: FrankenScript — the system's ONE TRUE INTERFACE (plan §11.1;
//! Decalogue P10). A typed, versioned intermediate representation with two
//! isomorphic concrete syntaxes: canonical s-expressions and a lossless
//! JSON mapping. Agents emit whichever their tooling prefers; BOTH parse
//! to the same typed AST (tested property, not aspiration).
//!
//! - Atoms are the system's real nouns: dimensioned quantities (fs-qty:
//!   `0.12Pa*s`, `65deg`, `0.5L/s`), counts (`384GiB`, `96cores`), seeds
//!   (`0xF00D0002`), strings, symbols, keywords.
//! - Every node carries a byte span; every parse failure is a structured
//!   [`IrError`] with the offending span and a fix hint (refusals teach).
//! - High-level verbs lower to explicit IR with an inspectable trace
//!   ([`lower::lower`]) — progressive disclosure with nothing hidden.
//! - [`study::Study`] recognizes the Appendix C study forms and extracts
//!   the Five Explicits' pillars (validity POLICY is the admission bead's).
//!
//! Layer: L6 (HELM). Runtime deps: `std` + fs-qty.

pub mod admission;
#[cfg(feature = "ladder-planner")]
pub mod anytime;
pub mod ast;
pub mod json;
pub mod lower;
#[cfg(feature = "ladder-planner")]
pub mod planner;
pub mod query;
pub mod sexpr;
pub mod study;

pub use ast::{CountUnit, CountValue, Node, NodeKind, Span};
pub use lower::{LowerStep, Lowered, lower};
pub use query::{FieldRegistry, Qoi, QoiMeta, Query, QueryAdmission, Target};
pub use study::Study;

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The IR language version this build reads and writes. Programs may pin
/// it; readers refuse newer language versions (never guess semantics).
pub const IR_VERSION: u32 = 1;

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
