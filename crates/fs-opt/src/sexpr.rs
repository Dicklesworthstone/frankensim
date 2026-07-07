//! Deterministic s-expression serialization (the plan's isomorphic-s-expr
//! direction; in-house — no serde, Franken-only). The parser rebuilds
//! through the TYPED constructors, so a tampered file with dimension
//! violations is REJECTED at parse time (re-validation is the point, not
//! an accident). Problem hash = FNV-64 over the canonical form.

use crate::graph::{Node, NodeId, Problem, ValidationError};
use crate::manifold::Manifold;
use fs_qty::Dims;

/// Serialize to the canonical s-expr form (deterministic: node-index
/// order, `{:?}` shortest-round-trip floats).
#[must_use]
pub fn to_sexpr(p: &Problem) -> String {
    use core::fmt::Write as _;
    let mut s = String::from("(problem\n (vars\n");
    for (name, m) in p.vars() {
        let mtxt = match m {
            Manifold::Euclidean(n) => format!("(euclidean {n})"),
            Manifold::Sphere(n) => format!("(sphere {n})"),
            Manifold::So3 => "(so3)".to_string(),
            Manifold::Stiefel(n, q) => format!("(stiefel {n} {q})"),
            Manifold::FixedVolumeLevelSet(n) => format!("(fvls {n})"),
        };
        let _ = writeln!(s, "  (var \"{name}\" {mtxt})");
    }
    s.push_str(" )\n (nodes\n");
    for node in p.nodes() {
        let line = match node {
            Node::Const(v, d) => {
                let dd = d.0;
                format!(
                    "(const {v:?} ({} {} {} {} {}))",
                    dd[0], dd[1], dd[2], dd[3], dd[4]
                )
            }
            Node::Component(v, i) => format!("(comp {} {})", v.0, i),
            Node::NormSq(v) => format!("(normsq {})", v.0),
            Node::Dot(u, w) => format!("(dot {} {})", u.0, w.0),
            Node::Add(a, b) => format!("(add {} {})", a.0, b.0),
            Node::Sub(a, b) => format!("(sub {} {})", a.0, b.0),
            Node::Mul(a, b) => format!("(mul {} {})", a.0, b.0),
            Node::Div(a, b) => format!("(div {} {})", a.0, b.0),
            Node::Neg(a) => format!("(neg {})", a.0),
            Node::Powi(a, n) => format!("(powi {} {n})", a.0),
            Node::Min(a, b) => format!("(min {} {})", a.0, b.0),
            Node::Max(a, b) => format!("(max {} {})", a.0, b.0),
            Node::Abs(a) => format!("(abs {})", a.0),
            Node::Sin(a) => format!("(sin {})", a.0),
            Node::Cos(a) => format!("(cos {})", a.0),
            Node::Exp(a) => format!("(exp {})", a.0),
            Node::Ln(a) => format!("(ln {})", a.0),
            Node::Sqrt(a) => format!("(sqrt {})", a.0),
            Node::Tanh(a) => format!("(tanh {})", a.0),
            Node::PdeResidual { study, adjoint_available } => {
                format!("(pde {study} {adjoint_available})")
            }
            Node::Expectation { inner, uq_config } => {
                format!("(expect {} {uq_config})", inner.0)
            }
        };
        let _ = writeln!(s, "  {line}");
    }
    let _ = writeln!(s, " )\n (objective {}))", p.objective().0);
    s
}

/// FNV-64 problem hash over the canonical serialization — study identity.
#[must_use]
pub fn problem_hash(p: &Problem) -> u64 {
    let s = to_sexpr(p);
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    for b in s.bytes() {
        acc ^= u64::from(b);
        acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
    }
    acc
}

/// Parse error: token position + message.
#[derive(Debug, Clone, PartialEq)]
pub enum ParseError {
    /// Malformed syntax.
    Syntax(String),
    /// Structurally valid but rejected by typed re-validation.
    Invalid(ValidationError),
}

impl From<ValidationError> for ParseError {
    fn from(e: ValidationError) -> ParseError {
        ParseError::Invalid(e)
    }
}

/// Parse the canonical form back into a [`Problem`], re-running every
/// typed constructor (dimension checks and class propagation included).
///
/// # Errors
/// [`ParseError`] on syntax or re-validation failure.
#[allow(clippy::too_many_lines)]
pub fn from_sexpr(text: &str) -> Result<Problem, ParseError> {
    let toks = tokenize(text);
    let mut p = Problem::new();
    let mut i = 0usize;
    let expect = |toks: &[String], i: &mut usize, t: &str| -> Result<(), ParseError> {
        if toks.get(*i).map(String::as_str) == Some(t) {
            *i += 1;
            Ok(())
        } else {
            Err(ParseError::Syntax(format!(
                "expected '{t}' at token {} (got {:?})",
                *i,
                toks.get(*i)
            )))
        }
    };
    let word = |toks: &[String], i: &mut usize| -> Result<String, ParseError> {
        let w = toks
            .get(*i)
            .ok_or_else(|| ParseError::Syntax("unexpected end".into()))?
            .clone();
        *i += 1;
        Ok(w)
    };
    expect(&toks, &mut i, "(")?;
    expect(&toks, &mut i, "problem")?;
    expect(&toks, &mut i, "(")?;
    expect(&toks, &mut i, "vars")?;
    while toks.get(i).map(String::as_str) == Some("(") {
        i += 1;
        expect(&toks, &mut i, "var")?;
        let name = word(&toks, &mut i)?;
        let name = name.trim_matches('"').to_string();
        expect(&toks, &mut i, "(")?;
        let kind = word(&toks, &mut i)?;
        let m = match kind.as_str() {
            "euclidean" => Manifold::Euclidean(parse_num(&word(&toks, &mut i)?)?),
            "sphere" => Manifold::Sphere(parse_num(&word(&toks, &mut i)?)?),
            "so3" => Manifold::So3,
            "stiefel" => {
                let n = parse_num(&word(&toks, &mut i)?)?;
                let q = parse_num(&word(&toks, &mut i)?)?;
                Manifold::Stiefel(n, q)
            }
            "fvls" => Manifold::FixedVolumeLevelSet(parse_num(&word(&toks, &mut i)?)?),
            other => return Err(ParseError::Syntax(format!("unknown manifold {other}"))),
        };
        expect(&toks, &mut i, ")")?;
        expect(&toks, &mut i, ")")?;
        p.variable(&name, m);
    }
    expect(&toks, &mut i, ")")?;
    expect(&toks, &mut i, "(")?;
    expect(&toks, &mut i, "nodes")?;
    let mut ids: Vec<NodeId> = Vec::new();
    while toks.get(i).map(String::as_str) == Some("(") {
        i += 1;
        let op = word(&toks, &mut i)?;
        let id = rebuild_node(&mut p, &op, &toks, &mut i, &ids)?;
        ids.push(id);
        expect(&toks, &mut i, ")")?;
    }
    expect(&toks, &mut i, ")")?;
    expect(&toks, &mut i, "(")?;
    expect(&toks, &mut i, "objective")?;
    let obj: usize = parse_num(&word(&toks, &mut i)?)?;
    p.set_objective(ids[obj]);
    Ok(p)
}

fn rebuild_node(
    p: &mut Problem,
    op: &str,
    toks: &[String],
    i: &mut usize,
    ids: &[NodeId],
) -> Result<NodeId, ParseError> {
    use crate::graph::VarId;
    let mut num = |i: &mut usize| -> Result<usize, ParseError> {
        let w = toks
            .get(*i)
            .ok_or_else(|| ParseError::Syntax("unexpected end".into()))?;
        *i += 1;
        parse_num(w)
    };
    let node_ref = |k: usize| -> Result<NodeId, ParseError> {
        ids.get(k)
            .copied()
            .ok_or_else(|| ParseError::Syntax(format!("forward node reference {k}")))
    };
    Ok(match op {
        "const" => {
            let w = toks
                .get(*i)
                .ok_or_else(|| ParseError::Syntax("unexpected end".into()))?
                .clone();
            *i += 1;
            let v: f64 = w
                .parse()
                .map_err(|_| ParseError::Syntax(format!("bad float {w}")))?;
            if toks.get(*i).map(String::as_str) != Some("(") {
                return Err(ParseError::Syntax("expected dims".into()));
            }
            *i += 1;
            let mut d = [0i8; 5];
            for slot in &mut d {
                let w = toks
                    .get(*i)
                    .ok_or_else(|| ParseError::Syntax("unexpected end".into()))?;
                *i += 1;
                *slot = w
                    .parse()
                    .map_err(|_| ParseError::Syntax(format!("bad dim {w}")))?;
            }
            if toks.get(*i).map(String::as_str) != Some(")") {
                return Err(ParseError::Syntax("expected close of dims".into()));
            }
            *i += 1;
            p.constant(v, Dims(d))
        }
        "comp" => {
            let v = num(i)?;
            let idx = num(i)?;
            p.component(VarId(u32::try_from(v).expect("var id")), u32::try_from(idx).expect("idx"))
        }
        "normsq" => {
            let v = num(i)?;
            p.norm_sq(VarId(u32::try_from(v).expect("var id")))
        }
        "dot" => {
            let u = num(i)?;
            let w = num(i)?;
            p.dot(
                VarId(u32::try_from(u).expect("var id")),
                VarId(u32::try_from(w).expect("var id")),
            )
        }
        "add" => p.add(node_ref(num(i)?)?, node_ref(num(i)?)?)?,
        "sub" => p.sub(node_ref(num(i)?)?, node_ref(num(i)?)?)?,
        "mul" => {
            let a = node_ref(num(i)?)?;
            let b = node_ref(num(i)?)?;
            p.mul(a, b)
        }
        "div" => {
            let a = node_ref(num(i)?)?;
            let b = node_ref(num(i)?)?;
            p.div(a, b)
        }
        "neg" => {
            let a = node_ref(num(i)?)?;
            p.neg(a)
        }
        "powi" => {
            let a = node_ref(num(i)?)?;
            let w = toks
                .get(*i)
                .ok_or_else(|| ParseError::Syntax("unexpected end".into()))?;
            *i += 1;
            let n: i8 = w
                .parse()
                .map_err(|_| ParseError::Syntax(format!("bad power {w}")))?;
            p.powi(a, n)
        }
        "min" => p.min(node_ref(num(i)?)?, node_ref(num(i)?)?)?,
        "max" => p.max(node_ref(num(i)?)?, node_ref(num(i)?)?)?,
        "abs" => {
            let a = node_ref(num(i)?)?;
            p.abs(a)
        }
        "sin" | "cos" | "exp" | "ln" | "sqrt" | "tanh" => {
            let a = node_ref(num(i)?)?;
            let op_static: &'static str = match op {
                "sin" => "sin",
                "cos" => "cos",
                "exp" => "exp",
                "ln" => "ln",
                "sqrt" => "sqrt",
                _ => "tanh",
            };
            p.unary(op_static, a)?
        }
        "pde" => {
            let study = num(i)? as u64;
            let w = toks
                .get(*i)
                .ok_or_else(|| ParseError::Syntax("unexpected end".into()))?;
            *i += 1;
            let adj = w == "true";
            p.pde_residual(study, adj)
        }
        "expect" => {
            let inner = node_ref(num(i)?)?;
            let cfg = num(i)? as u64;
            p.expectation(inner, cfg)
        }
        other => return Err(ParseError::Syntax(format!("unknown op {other}"))),
    })
}

fn parse_num<T: core::str::FromStr>(w: &str) -> Result<T, ParseError> {
    w.parse()
        .map_err(|_| ParseError::Syntax(format!("bad number {w}")))
}

fn tokenize(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    for ch in text.chars() {
        match ch {
            '(' | ')' => {
                if !cur.is_empty() {
                    out.push(core::mem::take(&mut cur));
                }
                out.push(ch.to_string());
            }
            c if c.is_whitespace() => {
                if !cur.is_empty() {
                    out.push(core::mem::take(&mut cur));
                }
            }
            c => cur.push(c),
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}
