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

impl CountUnit {
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
    /// Non-SI count (memory grants, core counts).
    Count {
        /// Magnitude in the written unit.
        value: f64,
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
        Node { kind, span: Span::default() }
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
                NodeKind::Qty { value: va, dims: da, .. },
                NodeKind::Qty { value: vb, dims: db, .. },
            ) => va.to_bits() == vb.to_bits() && da == db,
            (
                NodeKind::Count { value: va, unit: ua },
                NodeKind::Count { value: vb, unit: ub },
            ) => va.to_bits() == vb.to_bits() && ua == ub,
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
