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

/// The IR language version this build reads and writes. Bare syntax parsers are
/// intentionally syntax-only; use [`VersionedProgram`] for persisted/replayed
/// artifacts whose canonical identity must bind and enforce this version.
pub const IR_VERSION: u32 = 2;

/// A canonical, version-bound FrankenScript artifact envelope.
///
/// Both concrete encodings represent the same envelope AST:
/// `(frankensim-ir :version 2 :program <node>)`. Construction always writes the
/// current version; parsing rejects every unsupported version rather than
/// guessing migration semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct VersionedProgram {
    version: u32,
    program: Node,
}

impl VersionedProgram {
    /// Bind a program to the language version written by this build.
    #[must_use]
    pub fn current(program: Node) -> Self {
        Self {
            version: IR_VERSION,
            program,
        }
    }

    /// Parse and enforce a canonical s-expression envelope.
    ///
    /// # Errors
    /// Syntax/shape errors and every version other than [`IR_VERSION`] refuse
    /// with a structured [`IrError`].
    pub fn parse_sexpr(src: &str) -> Result<Self, IrError> {
        Self::from_envelope(sexpr::parse(src)?)
    }

    /// Parse and enforce a canonical JSON envelope.
    ///
    /// # Errors
    /// Syntax/shape errors and every version other than [`IR_VERSION`] refuse
    /// with a structured [`IrError`].
    pub fn parse_json(src: &str) -> Result<Self, IrError> {
        Self::from_envelope(json::parse(src)?)
    }

    /// Validate an already parsed envelope AST.
    ///
    /// # Errors
    /// Malformed envelopes and unsupported language versions refuse.
    pub fn from_envelope(node: Node) -> Result<Self, IrError> {
        let envelope_span = node.span;
        let NodeKind::List(mut items) = node.kind else {
            return Err(version_envelope_error(envelope_span));
        };
        if items.len() != 5
            || !matches!(&items[0].kind, NodeKind::Symbol(head) if head == "frankensim-ir")
            || !matches!(&items[1].kind, NodeKind::Keyword(key) if key == "version")
            || !matches!(&items[3].kind, NodeKind::Keyword(key) if key == "program")
        {
            return Err(version_envelope_error(envelope_span));
        }
        let version_span = items[2].span;
        let NodeKind::Int(written_version) = &items[2].kind else {
            return Err(version_envelope_error(version_span));
        };
        let Ok(version) = u32::try_from(*written_version) else {
            return Err(version_envelope_error(version_span));
        };
        if version != IR_VERSION {
            return Err(IrError {
                span: version_span,
                kind: IrErrorKind::UnsupportedVersion,
                detail: format!(
                    "IR language version {version} is unsupported by this version-{IR_VERSION} reader"
                ),
                hint: "use an explicit, audited migration before replay; never rewrite version semantics implicitly"
                    .to_string(),
            });
        }
        Ok(Self {
            version,
            program: items.swap_remove(4),
        })
    }

    /// Bound language version.
    #[must_use]
    pub const fn version(&self) -> u32 {
        self.version
    }

    /// The enclosed program AST.
    #[must_use]
    pub const fn program(&self) -> &Node {
        &self.program
    }

    /// Consume the envelope and return its program AST.
    #[must_use]
    pub fn into_program(self) -> Node {
        self.program
    }

    /// Canonical version-bound s-expression.
    #[must_use]
    pub fn print_sexpr(&self) -> String {
        sexpr::print(&self.envelope_node())
    }

    /// Canonical version-bound JSON mapping.
    #[must_use]
    pub fn print_json(&self) -> String {
        json::print(&self.envelope_node())
    }

    fn envelope_node(&self) -> Node {
        Node::synthetic(NodeKind::List(vec![
            Node::synthetic(NodeKind::Symbol("frankensim-ir".to_string())),
            Node::synthetic(NodeKind::Keyword("version".to_string())),
            Node::synthetic(NodeKind::Int(i64::from(self.version))),
            Node::synthetic(NodeKind::Keyword("program".to_string())),
            self.program.clone(),
        ]))
    }
}

fn version_envelope_error(span: Span) -> IrError {
    IrError {
        span,
        kind: IrErrorKind::MalformedClause,
        detail: "expected (frankensim-ir :version <u32> :program <node>)".to_string(),
        hint: "persist and replay programs through VersionedProgram; bare parsers are syntax-only"
            .to_string(),
    }
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
