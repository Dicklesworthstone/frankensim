//! Canonical serialization: problems round-trip through a line-based
//! text form whose floats are BIT PATTERNS (exact round-trip), and the
//! problem HASH (FNV-1a 64 over the canonical body) is the study
//! identity. Parsing rebuilds THROUGH the validating builder, so a
//! tampered file cannot smuggle in an ill-typed graph — revalidation is
//! free.

use crate::ir::{
    ConstraintKind, Expr, Manifold, NodeId, OptError, Problem, ProblemBuilder, ProblemTag, Sense,
    VarId,
};
use fs_qty::Dims;
use std::fmt::Write as _;

/// Escape a string for single-token embedding. Percent-encodes the token
/// delimiters AND every control / non-ASCII byte, so ANY value (including
/// multibyte UTF-8) round-trips exactly. Graphic ASCII except `%` passes
/// through verbatim, so ASCII fields serialize byte-identically to before.
fn esc(s: &str) -> String {
    use std::fmt::Write;
    let mut out = String::with_capacity(s.len());
    for &b in s.as_bytes() {
        if b.is_ascii_graphic() && b != b'%' {
            out.push(b as char);
        } else {
            let _ = write!(out, "%{b:02X}");
        }
    }
    out
}

/// Inverse of [`esc`]. Decodes `%XX` on the raw BYTE stream (so a multibyte
/// value reassembles correctly — the old byte-wise `as char` was a lossy
/// Latin-1 decode that corrupted every non-ASCII field), and never slices a
/// `str` at a non-char boundary, so a tampered token cannot panic.
fn unesc(s: &str) -> String {
    let src = s.as_bytes();
    let mut bytes = Vec::with_capacity(src.len());
    let mut i = 0;
    while i < src.len() {
        if src[i] == b'%'
            && i + 3 <= src.len()
            && let (Some(hi), Some(lo)) = (hex_nibble(src[i + 1]), hex_nibble(src[i + 2]))
        {
            bytes.push((hi << 4) | lo);
            i += 3;
            continue;
        }
        bytes.push(src[i]);
        i += 1;
    }
    String::from_utf8_lossy(&bytes).into_owned()
}

fn hex_nibble(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'a'..=b'f' => Some(b - b'a' + 10),
        b'A'..=b'F' => Some(b - b'A' + 10),
        _ => None,
    }
}

fn dims_str(d: Dims) -> String {
    format!("({},{},{},{},{})", d.0[0], d.0[1], d.0[2], d.0[3], d.0[4])
}

fn parse_dims(s: &str) -> Option<Dims> {
    let inner = s.strip_prefix('(')?.strip_suffix(')')?;
    let parts: Vec<&str> = inner.split(',').collect();
    if parts.len() != 5 {
        return None;
    }
    let mut d = [0i8; 5];
    for (slot, p) in d.iter_mut().zip(&parts) {
        *slot = p.parse().ok()?;
    }
    Some(Dims(d))
}

fn f64_hex(v: f64) -> String {
    format!("{:016X}", v.to_bits())
}

fn parse_f64_hex(s: &str) -> Option<f64> {
    u64::from_str_radix(s, 16).ok().map(f64::from_bits)
}

fn manifold_str(m: Manifold) -> String {
    match m {
        Manifold::Rn { dim } => format!("Rn({dim})"),
        Manifold::Sphere { ambient } => format!("Sphere({ambient})"),
        Manifold::So3 => "So3".to_string(),
        Manifold::Stiefel { n, p } => format!("Stiefel({n},{p})"),
    }
}

fn parse_manifold(s: &str) -> Option<Manifold> {
    if s == "So3" {
        return Some(Manifold::So3);
    }
    if let Some(inner) = s.strip_prefix("Rn(").and_then(|r| r.strip_suffix(')')) {
        return Some(Manifold::Rn {
            dim: inner.parse().ok()?,
        });
    }
    if let Some(inner) = s.strip_prefix("Sphere(").and_then(|r| r.strip_suffix(')')) {
        return Some(Manifold::Sphere {
            ambient: inner.parse().ok()?,
        });
    }
    if let Some(inner) = s.strip_prefix("Stiefel(").and_then(|r| r.strip_suffix(')')) {
        let (n, p) = inner.split_once(',')?;
        return Some(Manifold::Stiefel {
            n: n.parse().ok()?,
            p: p.parse().ok()?,
        });
    }
    None
}

/// Canonical single-token form of one expression (the hash-consing key
/// AND the serialized body).
#[must_use]
pub(crate) fn expr_key(e: &Expr) -> String {
    match e {
        Expr::Var(v) => format!("var {}", v.0),
        Expr::Component { of, index } => format!("component {} {index}", of.0),
        Expr::Const { value, dims } => {
            format!("const {} {}", f64_hex(*value), dims_str(*dims))
        }
        Expr::Add(a, b) => format!("add {} {}", a.0, b.0),
        Expr::Sub(a, b) => format!("sub {} {}", a.0, b.0),
        Expr::Mul(a, b) => format!("mul {} {}", a.0, b.0),
        Expr::Div(a, b) => format!("div {} {}", a.0, b.0),
        Expr::Neg(a) => format!("neg {}", a.0),
        Expr::Powi { base, exp } => format!("powi {} {exp}", base.0),
        Expr::Sqrt(a) => format!("sqrt {}", a.0),
        Expr::Exp(a) => format!("exp {}", a.0),
        Expr::Ln(a) => format!("ln {}", a.0),
        Expr::Tanh(a) => format!("tanh {}", a.0),
        Expr::Dot(a, b) => format!("dot {} {}", a.0, b.0),
        Expr::NormSq(a) => format!("norm_sq {}", a.0),
        Expr::Min(a, b) => format!("min {} {}", a.0, b.0),
        Expr::Max(a, b) => format!("max {} {}", a.0, b.0),
        Expr::Abs(a) => format!("abs {}", a.0),
        Expr::PdeResidual {
            study,
            over,
            adjoint_available,
            dims,
        } => format!(
            "pde_residual {} {} {} {}",
            esc(study),
            over.0,
            u8::from(*adjoint_available),
            dims_str(*dims)
        ),
        Expr::Expectation { of, uq_config } => {
            format!("expectation {} {}", of.0, esc(uq_config))
        }
        Expr::Cvar {
            of,
            alpha,
            uq_config,
        } => format!("cvar {} {} {}", of.0, f64_hex(*alpha), esc(uq_config)),
        Expr::Quantile { of, q, uq_config } => {
            format!("quantile {} {} {}", of.0, f64_hex(*q), esc(uq_config))
        }
    }
}

/// FNV-1a 64 (in-house; stable across platforms).
fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

fn body(problem: &Problem) -> String {
    let mut s = String::from("fsopt v1\n");
    for (i, v) in problem.vars.iter().enumerate() {
        let _ = writeln!(
            s,
            "var {i} {} {} {}",
            esc(&v.name),
            manifold_str(v.manifold),
            dims_str(v.dims)
        );
    }
    for (i, e) in problem.exprs.iter().enumerate() {
        let _ = writeln!(s, "expr {i} {}", expr_key(e));
    }
    for o in &problem.objectives {
        let sense = match o.sense {
            Sense::Minimize => "min",
            Sense::Maximize => "max",
        };
        let _ = writeln!(s, "objective {sense} {} {}", o.node.0, f64_hex(o.weight));
    }
    for c in &problem.constraints {
        let kind = match c.kind {
            ConstraintKind::EqZero => "eq0",
            ConstraintKind::LeZero => "le0",
        };
        let _ = writeln!(s, "constraint {kind} {} {}", c.node.0, esc(&c.name));
    }
    for t in &problem.tags {
        match t {
            ProblemTag::MultiFidelity { levels } => {
                let _ = writeln!(s, "tag multi_fidelity {levels}");
            }
            ProblemTag::ChanceConstrained { prob } => {
                let _ = writeln!(s, "tag chance {}", f64_hex(*prob));
            }
            ProblemTag::Bilevel { inner_hash } => {
                let _ = writeln!(s, "tag bilevel {inner_hash:016X}");
            }
        }
    }
    let _ = writeln!(s, "budget {}", problem.budget.max_evals);
    s
}

/// The problem hash: FNV-1a 64 over the canonical body (study
/// identity; two structurally identical builds hash identically).
#[must_use]
pub fn problem_hash(problem: &Problem) -> u64 {
    fnv1a(body(problem).as_bytes())
}

/// Serialize to the canonical text form (hash line appended).
#[must_use]
pub fn serialize(problem: &Problem) -> String {
    let b = body(problem);
    format!("{b}hash {:016X}\n", fnv1a(b.as_bytes()))
}

fn perr(line: usize, what: impl Into<String>) -> OptError {
    OptError::Parse {
        line,
        what: what.into(),
    }
}

/// Parse the canonical form, REBUILDING through the validating builder
/// (ill-typed files are rejected with the builder's teaching errors)
/// and verifying the integrity hash.
///
/// # Errors
/// [`OptError::Parse`] (with line numbers) or any builder error.
#[allow(clippy::too_many_lines)] // one grammar rule per block
pub fn parse(text: &str) -> Result<Problem, OptError> {
    let mut b = ProblemBuilder::new();
    let mut expected_hash: Option<u64> = None;
    let mut body_end = 0usize;
    let lines: Vec<&str> = text.lines().collect();
    for (ln0, line) in lines.iter().enumerate() {
        let ln = ln0 + 1;
        let mut tok = line.split(' ');
        let head = tok.next().unwrap_or("");
        match head {
            "fsopt" => {
                if tok.next() != Some("v1") {
                    return Err(perr(ln, "unsupported version (expected `fsopt v1`)"));
                }
            }
            "var" => {
                let _ix: u32 = tok
                    .next()
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "var: missing index"))?;
                let name = unesc(tok.next().ok_or_else(|| perr(ln, "var: missing name"))?);
                let manifold = tok
                    .next()
                    .and_then(parse_manifold)
                    .ok_or_else(|| perr(ln, "var: bad manifold"))?;
                let dims = tok
                    .next()
                    .and_then(parse_dims)
                    .ok_or_else(|| perr(ln, "var: bad dims"))?;
                b.var(&name, manifold, dims);
            }
            "expr" => {
                let ix: u32 = tok
                    .next()
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "expr: missing index"))?;
                let rest: Vec<&str> = tok.collect();
                let id = parse_expr(&mut b, &rest, ln)?;
                if id.0 != ix {
                    return Err(perr(
                        ln,
                        format!(
                            "expr id mismatch: file says {ix}, canonical rebuild gives {} \
                             (duplicate or reordered nodes break identity)",
                            id.0
                        ),
                    ));
                }
            }
            "objective" => {
                let sense = match tok.next() {
                    Some("min") => Sense::Minimize,
                    Some("max") => Sense::Maximize,
                    _ => return Err(perr(ln, "objective: bad sense")),
                };
                let node: u32 = tok
                    .next()
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "objective: bad node"))?;
                let weight = tok
                    .next()
                    .and_then(parse_f64_hex)
                    .ok_or_else(|| perr(ln, "objective: bad weight"))?;
                b.objective(NodeId(node), sense, weight)?;
            }
            "constraint" => {
                let kind = match tok.next() {
                    Some("eq0") => ConstraintKind::EqZero,
                    Some("le0") => ConstraintKind::LeZero,
                    _ => return Err(perr(ln, "constraint: bad kind")),
                };
                let node: u32 = tok
                    .next()
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "constraint: bad node"))?;
                let name = unesc(tok.next().unwrap_or(""));
                b.constraint(NodeId(node), kind, &name)?;
            }
            "tag" => match tok.next() {
                Some("multi_fidelity") => {
                    let levels = tok
                        .next()
                        .and_then(|t| t.parse().ok())
                        .ok_or_else(|| perr(ln, "tag: bad levels"))?;
                    b.tag(ProblemTag::MultiFidelity { levels });
                }
                Some("chance") => {
                    let prob = tok
                        .next()
                        .and_then(parse_f64_hex)
                        .ok_or_else(|| perr(ln, "tag: bad prob"))?;
                    b.tag(ProblemTag::ChanceConstrained { prob });
                }
                Some("bilevel") => {
                    let inner_hash = tok
                        .next()
                        .and_then(|t| u64::from_str_radix(t, 16).ok())
                        .ok_or_else(|| perr(ln, "tag: bad hash"))?;
                    b.tag(ProblemTag::Bilevel { inner_hash });
                }
                _ => return Err(perr(ln, "unknown tag")),
            },
            "budget" => {
                let n = tok
                    .next()
                    .and_then(|t| t.parse().ok())
                    .ok_or_else(|| perr(ln, "budget: bad count"))?;
                b.set_budget(n);
            }
            "hash" => {
                expected_hash = Some(
                    tok.next()
                        .and_then(|t| u64::from_str_radix(t, 16).ok())
                        .ok_or_else(|| perr(ln, "hash: bad value"))?,
                );
                body_end = ln0;
            }
            "" => {}
            other => return Err(perr(ln, format!("unknown directive `{other}`"))),
        }
    }
    let problem = b.finish();
    if let Some(h) = expected_hash {
        let mut body_text = lines[..body_end].join("\n");
        body_text.push('\n');
        let actual = fnv1a(body_text.as_bytes());
        if actual != h {
            return Err(perr(
                body_end + 1,
                format!("integrity hash mismatch: file {h:016X}, content {actual:016X}"),
            ));
        }
    }
    Ok(problem)
}

fn parse_expr(b: &mut ProblemBuilder, toks: &[&str], ln: usize) -> Result<NodeId, OptError> {
    let get = |i: usize| -> Result<&str, OptError> {
        toks.get(i)
            .copied()
            .ok_or_else(|| perr(ln, "expr: missing operand"))
    };
    let node = |i: usize| -> Result<NodeId, OptError> {
        Ok(NodeId(
            get(i)?
                .parse()
                .map_err(|_| perr(ln, "expr: bad node operand"))?,
        ))
    };
    match get(0)? {
        "var" => {
            let v: u32 = get(1)?.parse().map_err(|_| perr(ln, "var: bad id"))?;
            b.var_ref(VarId(v))
        }
        "component" => {
            let of = node(1)?;
            let index = get(2)?
                .parse()
                .map_err(|_| perr(ln, "component: bad index"))?;
            b.component(of, index)
        }
        "const" => {
            let value = parse_f64_hex(get(1)?).ok_or_else(|| perr(ln, "const: bad value"))?;
            let dims = parse_dims(get(2)?).ok_or_else(|| perr(ln, "const: bad dims"))?;
            Ok(b.konst(value, dims))
        }
        "add" => b.add(node(1)?, node(2)?),
        "sub" => b.sub(node(1)?, node(2)?),
        "mul" => b.mul(node(1)?, node(2)?),
        "div" => b.div(node(1)?, node(2)?),
        "neg" => b.neg(node(1)?),
        "powi" => {
            let base = node(1)?;
            let exp = get(2)?
                .parse()
                .map_err(|_| perr(ln, "powi: bad exponent"))?;
            b.powi(base, exp)
        }
        "sqrt" => b.sqrt(node(1)?),
        "exp" => b.exp(node(1)?),
        "ln" => b.ln(node(1)?),
        "tanh" => b.tanh(node(1)?),
        "dot" => b.dot(node(1)?, node(2)?),
        "norm_sq" => b.norm_sq(node(1)?),
        "min" => b.min_of(node(1)?, node(2)?),
        "max" => b.max_of(node(1)?, node(2)?),
        "abs" => b.abs(node(1)?),
        "pde_residual" => {
            let study = unesc(get(1)?);
            let over: u32 = get(2)?.parse().map_err(|_| perr(ln, "pde: bad var"))?;
            let adj = get(3)? == "1";
            let dims = parse_dims(get(4)?).ok_or_else(|| perr(ln, "pde: bad dims"))?;
            b.pde_residual(&study, VarId(over), adj, dims)
        }
        "expectation" => {
            let of = node(1)?;
            let cfg = unesc(get(2)?);
            b.expectation(of, &cfg)
        }
        "cvar" => {
            let of = node(1)?;
            let alpha = parse_f64_hex(get(2)?).ok_or_else(|| perr(ln, "cvar: bad alpha"))?;
            let cfg = unesc(get(3)?);
            b.cvar(of, alpha, &cfg)
        }
        "quantile" => {
            let of = node(1)?;
            let q = parse_f64_hex(get(2)?).ok_or_else(|| perr(ln, "quantile: bad q"))?;
            let cfg = unesc(get(3)?);
            b.quantile(of, q, &cfg)
        }
        other => Err(perr(ln, format!("unknown expr kind `{other}`"))),
    }
}

#[cfg(test)]
mod tests {
    use super::{esc, unesc};

    #[test]
    fn esc_unesc_round_trips_utf8_and_delimiters() {
        // Non-ASCII, delimiters, and control bytes must round-trip EXACTLY —
        // the old byte-wise `as char` decode was a lossy Latin-1 pass that
        // corrupted every multibyte field (breaking the BITWISE round-trip
        // contract) and could panic slicing a str at a non-char boundary.
        for s in ["café", "α+β·γ", "a b\n%c", "π²∇🎯", "plain-ascii_123", "", "100%"] {
            assert_eq!(unesc(&esc(s)), *s, "round-trip failed for {s:?}");
        }
        // ASCII fields keep the exact prior wire format (backward compatible).
        assert_eq!(esc("a b%c\n"), "a%20b%25c%0A");
        assert_eq!(unesc("a%20b%25c%0A"), "a b%c\n");
        // Crafted/truncated tokens must NOT panic (a literal `%` before a
        // 3-byte char used to split a str mid-character).
        let _ = unesc("%\u{20AC}x");
        let _ = unesc("%");
        let _ = unesc("%2");
    }
}
