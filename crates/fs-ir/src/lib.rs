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
//! - [`machine`] provides the default `[S]` durable entity/topology-lineage
//!   kernel, admitted machine graph, and separately identified behavior overlay.
//!
//! Layer: L6 (HELM). Production dependencies and feature-gated deltas are
//! declared explicitly in `Cargo.toml`.

pub mod admission;
#[cfg(feature = "ladder-planner")]
pub mod anytime;
pub mod ast;
pub mod campaign;
pub mod catalog;
pub mod conformance;
pub mod cross_contract;
#[cfg(feature = "derived-crosswalk")]
pub mod derived_crosswalk;
pub mod json;
pub mod lower;
pub mod machine;
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
pub const IR_VERSION: u32 = 3;

/// A canonical, version-bound FrankenScript artifact envelope.
///
/// Both concrete encodings represent the same envelope AST:
/// `(frankensim-ir :version 3 :program <node>)`. Construction always writes the
/// current version; parsing rejects every unsupported version rather than
/// guessing migration semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct VersionedProgram {
    version: u32,
    program: Node,
}

impl VersionedProgram {
    /// Bind a trusted program to the language version written by this build.
    ///
    /// This compatibility constructor panics for a caller-forged invalid AST.
    /// New boundary code should use [`VersionedProgram::try_current`].
    #[must_use]
    pub fn current(program: Node) -> Self {
        Self::try_current(program).expect("invalid AST passed to VersionedProgram::current")
    }

    /// Validate a program and bind it to the current language version.
    ///
    /// # Errors
    /// Rejects the first invalid AST atom with its exact tree path and span,
    /// including a program too deep to fit inside the persisted envelope.
    pub fn try_current(program: Node) -> Result<Self, IrError> {
        let artifact = Self {
            version: IR_VERSION,
            program,
        };
        artifact.envelope_node().validate()?;
        Ok(artifact)
    }

    /// Parse and enforce a canonical s-expression envelope.
    ///
    /// # Errors
    /// Syntax/shape errors, noncanonical persisted bytes, and every version
    /// other than [`IR_VERSION`] refuse with a structured [`IrError`].
    pub fn parse_sexpr(src: &str) -> Result<Self, IrError> {
        let artifact = Self::from_envelope(sexpr::parse(src)?)?;
        require_canonical_input(src, &artifact.print_sexpr_checked()?, "s-expression")?;
        Ok(artifact)
    }

    /// Parse and enforce a canonical JSON envelope.
    ///
    /// # Errors
    /// Syntax/shape errors, noncanonical persisted bytes, and every version
    /// other than [`IR_VERSION`] refuse with a structured [`IrError`].
    pub fn parse_json(src: &str) -> Result<Self, IrError> {
        let artifact = Self::from_envelope(json::parse(src)?)?;
        require_canonical_input(src, &artifact.print_json_checked()?, "JSON")?;
        Ok(artifact)
    }

    /// Validate an already parsed envelope AST.
    ///
    /// # Errors
    /// Malformed envelopes and unsupported language versions refuse.
    pub fn from_envelope(node: Node) -> Result<Self, IrError> {
        node.validate()?;
        let envelope_span = node.span;
        let NodeKind::List(mut items) = node.kind else {
            return Err(version_envelope_error(
                envelope_span,
                "version envelope root $ must be a list",
            ));
        };
        let head = items.first().ok_or_else(|| {
            version_envelope_error(envelope_span, "version envelope is missing head $[0]")
        })?;
        if !matches!(&head.kind, NodeKind::Symbol(value) if value == "frankensim-ir") {
            return Err(version_envelope_error(
                head.span,
                "expected symbol frankensim-ir at version envelope $[0]",
            ));
        }
        let version_key = items.get(1).ok_or_else(|| {
            version_envelope_error(
                Span::new(envelope_span.end, envelope_span.end),
                "version envelope is missing keyword :version at $[1]",
            )
        })?;
        if !matches!(&version_key.kind, NodeKind::Keyword(key) if key == "version") {
            let detail = match &version_key.kind {
                NodeKind::Keyword(key) => format!(
                    "unknown or out-of-order version envelope keyword :{key} at $[1]; expected :version"
                ),
                _ => "expected keyword :version at version envelope $[1]".to_string(),
            };
            return Err(version_envelope_error(version_key.span, &detail));
        }
        let version_node = items.get(2).ok_or_else(|| {
            version_envelope_error(
                Span::new(envelope_span.end, envelope_span.end),
                "version envelope is missing u32 value at $[2]",
            )
        })?;
        let version_span = version_node.span;
        let NodeKind::Int(written_version) = &version_node.kind else {
            return Err(version_envelope_error(
                version_span,
                "expected an integer language version at version envelope $[2]",
            ));
        };
        let Ok(version) = u32::try_from(*written_version) else {
            return Err(version_envelope_error(
                version_span,
                "language version at version envelope $[2] is outside u32",
            ));
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
        let program_key = items.get(3).ok_or_else(|| {
            version_envelope_error(
                Span::new(envelope_span.end, envelope_span.end),
                "version envelope is missing keyword :program at $[3]",
            )
        })?;
        if !matches!(&program_key.kind, NodeKind::Keyword(key) if key == "program") {
            let detail = match &program_key.kind {
                NodeKind::Keyword(key) if key == "version" => {
                    "duplicate version envelope keyword :version at $[3]".to_string()
                }
                NodeKind::Keyword(key) => {
                    format!("unknown version envelope keyword :{key} at $[3]")
                }
                _ => "expected keyword :program at version envelope $[3]".to_string(),
            };
            return Err(version_envelope_error(program_key.span, &detail));
        }
        if items.len() < 5 {
            return Err(version_envelope_error(
                Span::new(envelope_span.end, envelope_span.end),
                "version envelope is missing program value at $[4]",
            ));
        }
        if let Some(trailing) = items.get(5) {
            return Err(version_envelope_error(
                trailing.span,
                "unexpected trailing version envelope value at $[5]",
            ));
        }
        Ok(Self {
            version,
            program: items.pop().expect("exact envelope length checked"),
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
        self.print_sexpr_checked()
            .expect("VersionedProgram contains an invalid AST")
    }

    /// Checked canonical version-bound s-expression.
    ///
    /// # Errors
    /// Refuses an invalid program AST rather than emitting ambiguous bytes.
    pub fn print_sexpr_checked(&self) -> Result<String, IrError> {
        sexpr::print_checked(&self.envelope_node())
    }

    /// Canonical version-bound JSON mapping.
    #[must_use]
    pub fn print_json(&self) -> String {
        self.print_json_checked()
            .expect("VersionedProgram contains an invalid AST")
    }

    /// Checked canonical version-bound JSON mapping.
    ///
    /// # Errors
    /// Refuses an invalid program AST rather than emitting ambiguous bytes.
    pub fn print_json_checked(&self) -> Result<String, IrError> {
        json::print_checked(&self.envelope_node())
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

fn version_envelope_error(span: Span, detail: &str) -> IrError {
    IrError {
        span,
        kind: IrErrorKind::MalformedClause,
        detail: detail.to_string(),
        hint: "persist and replay programs through VersionedProgram; bare parsers are syntax-only"
            .to_string(),
    }
}

fn require_canonical_input(src: &str, canonical: &str, syntax: &str) -> Result<(), IrError> {
    if src == canonical {
        return Ok(());
    }
    Err(IrError {
        span: first_text_difference_span(src, canonical),
        kind: IrErrorKind::NonCanonical,
        detail: format!("versioned {syntax} artifact is not in its one canonical byte encoding"),
        hint: "parse untrusted bare syntax, validate/migrate it explicitly, then persist the checked VersionedProgram rendering"
            .to_string(),
    })
}

fn first_text_difference_span(input: &str, canonical: &str) -> Span {
    let mut input_chars = input.char_indices();
    let mut canonical_chars = canonical.chars();
    loop {
        match (input_chars.next(), canonical_chars.next()) {
            (Some((_offset, input_char)), Some(canonical_char)) if input_char == canonical_char => {
            }
            (Some((offset, input_char)), _) => {
                return Span::new(offset, offset + input_char.len_utf8());
            }
            (None, Some(_)) => return Span::new(input.len(), input.len()),
            (None, None) => unreachable!("unequal strings must have a first difference"),
        }
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
