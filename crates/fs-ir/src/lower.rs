//! Verb lowering (plan §11.1): high-level verbs exist for terseness but
//! LOWER TO EXPLICIT IR, and the lowering is inspectable — progressive
//! disclosure with nothing hidden. Every injected default is named in the
//! trace, so an agent can see exactly what its shorthand committed it to.

use crate::ast::{Node, NodeKind, Span};
use crate::{IrError, IrErrorKind, sexpr};

/// One lowering event.
#[derive(Debug, Clone)]
pub struct LowerStep {
    /// The verb that was expanded.
    pub verb: String,
    /// Where the verb appeared in the source.
    pub span: Span,
    /// Defaults the lowering injected (named explicitly — P10).
    pub injected: Vec<String>,
    /// The explicit form it became (canonical s-expr rendering).
    pub expansion: String,
}

/// A lowered program with its inspection trace.
#[derive(Debug, Clone)]
pub struct Lowered {
    /// The fully explicit program.
    pub node: Node,
    /// Every expansion performed, in traversal order.
    pub trace: Vec<LowerStep>,
}

/// Lower every known verb in the tree (recursively; innermost first).
///
/// # Errors
/// Structured [`IrError`] pointing at a malformed verb usage.
pub fn lower(node: &Node) -> Result<Lowered, IrError> {
    let mut trace = Vec::new();
    let node = lower_inner(node, &mut trace)?;
    Ok(Lowered { node, trace })
}

fn lower_inner(node: &Node, trace: &mut Vec<LowerStep>) -> Result<Node, IrError> {
    let lowered_children = match &node.kind {
        NodeKind::List(items) => {
            let children: Result<Vec<Node>, IrError> =
                items.iter().map(|n| lower_inner(n, trace)).collect();
            Node {
                kind: NodeKind::List(children?),
                span: node.span,
            }
        }
        _ => node.clone(),
    };
    match lowered_children.head() {
        Some("optimize-shape") => lower_optimize_shape(&lowered_children, trace),
        Some("simulate-pour") => lower_simulate_pour(&lowered_children, trace),
        _ => Ok(lowered_children),
    }
}

fn sym(s: &str) -> Node {
    Node::synthetic(NodeKind::Symbol(s.to_string()))
}

fn kw(s: &str) -> Node {
    Node::synthetic(NodeKind::Keyword(s.to_string()))
}

fn list(items: Vec<Node>) -> Node {
    Node::synthetic(NodeKind::List(items))
}

/// `(optimize-shape :min J :over X [:method M] [:until U])` →
/// `(ascent.optimize (min J) :over X :method M :until U :emit (ledger report))`
fn lower_optimize_shape(node: &Node, trace: &mut Vec<LowerStep>) -> Result<Node, IrError> {
    let items = node.items().expect("verb head implies list");
    let missing = |what: &str| IrError {
        span: node.span,
        kind: IrErrorKind::MalformedClause,
        detail: format!("optimize-shape needs {what}"),
        hint: "(optimize-shape :min <objective> :over <levers> [:method M] [:until U])".to_string(),
    };
    let args = &items[1..];
    if !args.len().is_multiple_of(2) {
        let trailing = args.last().unwrap_or(node);
        let detail = match &trailing.kind {
            NodeKind::Keyword(name) => format!("optimize-shape has dangling :{name} with no value"),
            _ => "optimize-shape has a trailing argument outside a keyword/value pair".to_string(),
        };
        return Err(IrError {
            span: trailing.span,
            kind: IrErrorKind::MalformedClause,
            detail,
            hint: "use exact :min/:over/:method/:until keyword/value pairs".to_string(),
        });
    }

    let mut objective = None;
    let mut over = None;
    let mut method = None;
    let mut until = None;
    for pair in args.as_chunks::<2>().0 {
        let NodeKind::Keyword(name) = &pair[0].kind else {
            return Err(IrError {
                span: pair[0].span,
                kind: IrErrorKind::MalformedClause,
                detail: "optimize-shape argument names must be keywords".to_string(),
                hint: "use exact :min/:over/:method/:until keyword/value pairs".to_string(),
            });
        };
        let slot = match name.as_str() {
            "min" => &mut objective,
            "over" => &mut over,
            "method" => &mut method,
            "until" => &mut until,
            _ => {
                return Err(IrError {
                    span: pair[0].span,
                    kind: IrErrorKind::MalformedClause,
                    detail: format!("unknown optimize-shape argument :{name}"),
                    hint: "use only :min, :over, :method, and :until".to_string(),
                });
            }
        };
        if slot.replace(&pair[1]).is_some() {
            return Err(IrError {
                span: pair[0].span,
                kind: IrErrorKind::MalformedClause,
                detail: format!("duplicate optimize-shape argument :{name} is ambiguous"),
                hint: format!("retain exactly one :{name} value"),
            });
        }
    }

    let objective = objective.ok_or_else(|| missing(":min <objective>"))?;
    let over = over.ok_or_else(|| missing(":over <levers>"))?;
    let mut injected = Vec::new();
    let method = method.cloned().unwrap_or_else(|| {
        injected.push(":method (lbfgs :m 17)".to_string());
        list(vec![
            sym("lbfgs"),
            kw("m"),
            Node::synthetic(NodeKind::Int(17)),
        ])
    });
    let until = until.cloned().unwrap_or_else(|| {
        injected.push(":until (grad-norm 1e-5)".to_string());
        list(vec![
            sym("grad-norm"),
            Node::synthetic(NodeKind::Float(1e-5)),
        ])
    });
    injected.push(":emit (ledger report)".to_string());
    let explicit = list(vec![
        sym("ascent.optimize"),
        list(vec![sym("min"), objective.clone()]),
        kw("over"),
        over.clone(),
        kw("method"),
        method,
        kw("until"),
        until,
        kw("emit"),
        list(vec![sym("ledger"), sym("report")]),
    ]);
    trace.push(LowerStep {
        verb: "optimize-shape".to_string(),
        span: node.span,
        injected,
        expansion: sexpr::print(&explicit)?,
    });
    Ok(explicit)
}

/// `(simulate-pour vessel fluid schedule)` →
/// `(flux.free-surface-lbm vessel fluid schedule)`
fn lower_simulate_pour(node: &Node, trace: &mut Vec<LowerStep>) -> Result<Node, IrError> {
    let items = node.items().expect("verb head implies list");
    if items.len() != 4 {
        return Err(IrError {
            span: node.span,
            kind: IrErrorKind::MalformedClause,
            detail: format!("simulate-pour takes 3 arguments, got {}", items.len() - 1),
            hint: "(simulate-pour <vessel> <fluid> <schedule>)".to_string(),
        });
    }
    let mut explicit_items = vec![sym("flux.free-surface-lbm")];
    explicit_items.extend(items[1..].iter().cloned());
    let explicit = list(explicit_items);
    trace.push(LowerStep {
        verb: "simulate-pour".to_string(),
        span: node.span,
        injected: Vec::new(),
        expansion: sexpr::print(&explicit)?,
    });
    Ok(explicit)
}
