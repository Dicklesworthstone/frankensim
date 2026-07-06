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
            Node { kind: NodeKind::List(children?), span: node.span }
        }
        _ => node.clone(),
    };
    match lowered_children.head() {
        Some("optimize-shape") => lower_optimize_shape(&lowered_children, trace),
        Some("simulate-pour") => lower_simulate_pour(&lowered_children, trace),
        _ => Ok(lowered_children),
    }
}

/// Fetch the node following keyword `:name` in an argument list.
fn kw_get<'a>(items: &'a [Node], name: &str) -> Option<&'a Node> {
    items.windows(2).find_map(|pair| match &pair[0].kind {
        NodeKind::Keyword(k) if k == name => Some(&pair[1]),
        _ => None,
    })
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
        hint: "(optimize-shape :min <objective> :over <levers> [:method M] [:until U])"
            .to_string(),
    };
    let objective = kw_get(items, "min").ok_or_else(|| missing(":min <objective>"))?;
    let over = kw_get(items, "over").ok_or_else(|| missing(":over <levers>"))?;
    let mut injected = Vec::new();
    let method = kw_get(items, "method").cloned().unwrap_or_else(|| {
        injected.push(":method (lbfgs :m 17)".to_string());
        list(vec![sym("lbfgs"), kw("m"), Node::synthetic(NodeKind::Int(17))])
    });
    let until = kw_get(items, "until").cloned().unwrap_or_else(|| {
        injected.push(":until (grad-norm 1e-5)".to_string());
        list(vec![sym("grad-norm"), Node::synthetic(NodeKind::Float(1e-5))])
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
        expansion: sexpr::print(&explicit),
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
        expansion: sexpr::print(&explicit),
    });
    Ok(explicit)
}
