//! fs-shapeprog — generative geometry program synthesis. Layer: L2.
//!
//! Continuous optimization refines a topology; it rarely INVENTS a grammar.
//! This is the discrete-invention medium: a typed constructive-geometry DSL
//! with SDF semantics, a rewrite engine that SIMPLIFIES and CANONICALIZES
//! programs under geometric identities, and seeded shape-grammar derivation.
//!
//! The load-bearing safety property (the acceptance criterion): a rewrite
//! PRESERVES GEOMETRY within its declared certificate — verified by SDF
//! sampling ([`max_sdf_discrepancy`]). Exact identities (offset composition
//! `offset(offset(a,r₁),r₂) = offset(a, r₁+r₂)`, union commutativity, transform
//! distribution, empty-identity) leave the SDF unchanged; a certified-
//! approximate one (dropping an offset below tolerance) changes it by at most
//! the stated bound. [`canonical_hash`] gives equivalent programs one identity
//! for archive/ledger dedup. Deterministic; no dependencies.

/// A primitive shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shape {
    /// A sphere of the given radius.
    Sphere,
    /// A cube of the given half-extent.
    Cube,
}

/// A constructive-geometry program (an SDF-valued expression).
#[derive(Debug, Clone, PartialEq)]
pub enum Geom {
    /// The empty set (SDF `+∞`).
    Empty,
    /// A primitive with a size parameter (radius / half-extent).
    Primitive {
        /// The shape.
        shape: Shape,
        /// The size (radius / half-extent).
        size: f64,
    },
    /// Boolean union (SDF `min`).
    Union(Box<Geom>, Box<Geom>),
    /// Boolean intersection (SDF `max`).
    Intersect(Box<Geom>, Box<Geom>),
    /// Boolean difference `a \ b` (SDF `max(a, −b)`).
    Difference(Box<Geom>, Box<Geom>),
    /// Grow by `radius` (SDF `child − radius`).
    Offset {
        /// The child.
        child: Box<Geom>,
        /// The offset radius.
        radius: f64,
    },
    /// Translate the child.
    Translate {
        /// The child.
        child: Box<Geom>,
        /// The translation.
        t: [f64; 3],
    },
}

fn sphere(radius: f64) -> Geom {
    Geom::Primitive {
        shape: Shape::Sphere,
        size: radius,
    }
}

impl Geom {
    /// A sphere primitive.
    #[must_use]
    pub fn sphere(radius: f64) -> Geom {
        sphere(radius)
    }
    /// A cube primitive (half-extent).
    #[must_use]
    pub fn cube(half: f64) -> Geom {
        Geom::Primitive {
            shape: Shape::Cube,
            size: half,
        }
    }
    /// Union (owning builder).
    #[must_use]
    pub fn union(self, other: Geom) -> Geom {
        Geom::Union(Box::new(self), Box::new(other))
    }
    /// Offset (owning builder).
    #[must_use]
    pub fn offset(self, radius: f64) -> Geom {
        Geom::Offset {
            child: Box::new(self),
            radius,
        }
    }
    /// Translate (owning builder).
    #[must_use]
    pub fn translate(self, t: [f64; 3]) -> Geom {
        Geom::Translate {
            child: Box::new(self),
            t,
        }
    }

    /// The signed distance at point `p`.
    #[must_use]
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        match self {
            Geom::Empty => f64::INFINITY,
            Geom::Primitive { shape, size } => match shape {
                Shape::Sphere => norm(p) - size,
                Shape::Cube => cube_sdf(p, *size),
            },
            Geom::Union(a, b) => a.sdf(p).min(b.sdf(p)),
            Geom::Intersect(a, b) => a.sdf(p).max(b.sdf(p)),
            Geom::Difference(a, b) => a.sdf(p).max(-b.sdf(p)),
            Geom::Offset { child, radius } => child.sdf(p) - radius,
            Geom::Translate { child, t } => child.sdf([p[0] - t[0], p[1] - t[1], p[2] - t[2]]),
        }
    }

    /// Print as an s-expression.
    #[must_use]
    pub fn to_sexpr(&self) -> String {
        match self {
            Geom::Empty => "(empty)".to_string(),
            Geom::Primitive { shape, size } => {
                let s = match shape {
                    Shape::Sphere => "sphere",
                    Shape::Cube => "cube",
                };
                format!("({s} {})", fmt(*size))
            }
            Geom::Union(a, b) => format!("(union {} {})", a.to_sexpr(), b.to_sexpr()),
            Geom::Intersect(a, b) => format!("(intersect {} {})", a.to_sexpr(), b.to_sexpr()),
            Geom::Difference(a, b) => format!("(difference {} {})", a.to_sexpr(), b.to_sexpr()),
            Geom::Offset { child, radius } => {
                format!("(offset {} {})", child.to_sexpr(), fmt(*radius))
            }
            Geom::Translate { child, t } => format!(
                "(translate {} {} {} {})",
                child.to_sexpr(),
                fmt(t[0]),
                fmt(t[1]),
                fmt(t[2])
            ),
        }
    }

    /// The canonical form: commutative operands (union/intersect) sorted by
    /// their canonical printing.
    #[must_use]
    pub fn canonical(&self) -> Geom {
        match self {
            Geom::Union(a, b) => {
                let (ca, cb) = order(a.canonical(), b.canonical());
                Geom::Union(Box::new(ca), Box::new(cb))
            }
            Geom::Intersect(a, b) => {
                let (ca, cb) = order(a.canonical(), b.canonical());
                Geom::Intersect(Box::new(ca), Box::new(cb))
            }
            Geom::Difference(a, b) => {
                Geom::Difference(Box::new(a.canonical()), Box::new(b.canonical()))
            }
            Geom::Offset { child, radius } => Geom::Offset {
                child: Box::new(child.canonical()),
                radius: *radius,
            },
            Geom::Translate { child, t } => Geom::Translate {
                child: Box::new(child.canonical()),
                t: *t,
            },
            leaf => leaf.clone(),
        }
    }

    /// A content hash of the canonical form — equivalent programs share it
    /// (archive/ledger dedup).
    #[must_use]
    pub fn canonical_hash(&self) -> u64 {
        fnv1a(self.canonical().to_sexpr().as_bytes())
    }

    /// Node count (program size).
    #[must_use]
    pub fn size(&self) -> usize {
        match self {
            Geom::Empty | Geom::Primitive { .. } => 1,
            Geom::Union(a, b) | Geom::Intersect(a, b) | Geom::Difference(a, b) => {
                1 + a.size() + b.size()
            }
            Geom::Offset { child, .. } | Geom::Translate { child, .. } => 1 + child.size(),
        }
    }
}

fn order(a: Geom, b: Geom) -> (Geom, Geom) {
    if a.to_sexpr() <= b.to_sexpr() {
        (a, b)
    } else {
        (b, a)
    }
}

/// The maximum `|SDF_a − SDF_b|` over the sample points — the rewrite-safety
/// check (both non-finite at a point count as agreement).
#[must_use]
pub fn max_sdf_discrepancy(a: &Geom, b: &Geom, samples: &[[f64; 3]]) -> f64 {
    let mut worst = 0.0_f64;
    for &p in samples {
        let (da, db) = (a.sdf(p), b.sdf(p));
        if da.is_infinite() && db.is_infinite() && da.signum() == db.signum() {
            continue;
        }
        worst = worst.max((da - db).abs());
    }
    worst
}

// -- Rewrite engine ---------------------------------------------------------

/// A rewrite's fidelity certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Certificate {
    /// The SDF is preserved exactly.
    Exact,
    /// The SDF changes by at most `bound`.
    Approximate {
        /// The certified error bound.
        bound: f64,
    },
}

/// One applied rewrite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rewrite {
    /// The rule name.
    pub rule: &'static str,
    /// Its certificate.
    pub certificate: Certificate,
}

/// The result of simplifying a program.
#[derive(Debug, Clone, PartialEq)]
pub struct Simplified {
    /// The simplified program.
    pub program: Geom,
    /// The rewrites applied (in order).
    pub rewrites: Vec<Rewrite>,
    /// The maximum certified SDF error introduced (`0` if all exact).
    pub max_error: f64,
}

/// Simplify a program to a fixpoint under the geometric-identity rewrites.
/// Offsets smaller than `tiny_offset_tol` are dropped with a certified bound;
/// every other rewrite is exact.
#[must_use]
pub fn simplify(g: &Geom, tiny_offset_tol: f64) -> Simplified {
    let mut current = g.clone();
    let mut rewrites = Vec::new();
    for _ in 0..64 {
        let before = current.to_sexpr();
        current = rewrite_pass(&current, tiny_offset_tol, &mut rewrites);
        if current.to_sexpr() == before {
            break;
        }
    }
    let max_error = rewrites
        .iter()
        .map(|r| match r.certificate {
            Certificate::Exact => 0.0,
            Certificate::Approximate { bound } => bound,
        })
        .fold(0.0_f64, f64::max);
    Simplified {
        program: current,
        rewrites,
        max_error,
    }
}

fn rewrite_pass(g: &Geom, tol: f64, log: &mut Vec<Rewrite>) -> Geom {
    // simplify children first (bottom-up).
    let g = match g {
        Geom::Union(a, b) => Geom::Union(
            Box::new(rewrite_pass(a, tol, log)),
            Box::new(rewrite_pass(b, tol, log)),
        ),
        Geom::Intersect(a, b) => Geom::Intersect(
            Box::new(rewrite_pass(a, tol, log)),
            Box::new(rewrite_pass(b, tol, log)),
        ),
        Geom::Difference(a, b) => Geom::Difference(
            Box::new(rewrite_pass(a, tol, log)),
            Box::new(rewrite_pass(b, tol, log)),
        ),
        Geom::Offset { child, radius } => Geom::Offset {
            child: Box::new(rewrite_pass(child, tol, log)),
            radius: *radius,
        },
        Geom::Translate { child, t } => Geom::Translate {
            child: Box::new(rewrite_pass(child, tol, log)),
            t: *t,
        },
        leaf => leaf.clone(),
    };
    // apply a root rule.
    match g {
        // offset composition (EXACT for signed distance fields).
        Geom::Offset { child, radius: r2 } => {
            if let Geom::Offset {
                child: inner,
                radius: r1,
            } = *child
            {
                log.push(Rewrite {
                    rule: "offset-compose",
                    certificate: Certificate::Exact,
                });
                Geom::Offset {
                    child: inner,
                    radius: r1 + r2,
                }
            } else if r2.abs() < tol {
                // drop a tiny offset — certified approximate within |r|.
                log.push(Rewrite {
                    rule: "drop-tiny-offset",
                    certificate: Certificate::Approximate { bound: r2.abs() },
                });
                *child
            } else {
                Geom::Offset { child, radius: r2 }
            }
        }
        // empty identities (EXACT).
        Geom::Union(a, b) => match (*a, *b) {
            (Geom::Empty, x) | (x, Geom::Empty) => {
                log.push(Rewrite {
                    rule: "union-identity",
                    certificate: Certificate::Exact,
                });
                x
            }
            (a, b) => Geom::Union(Box::new(a), Box::new(b)),
        },
        Geom::Difference(a, b) => match *b {
            Geom::Empty => {
                log.push(Rewrite {
                    rule: "difference-identity",
                    certificate: Certificate::Exact,
                });
                *a
            }
            b => Geom::Difference(a, Box::new(b)),
        },
        Geom::Intersect(a, b) => match (*a, *b) {
            (Geom::Empty, _) | (_, Geom::Empty) => {
                log.push(Rewrite {
                    rule: "intersect-empty",
                    certificate: Certificate::Exact,
                });
                Geom::Empty
            }
            (a, b) => Geom::Intersect(Box::new(a), Box::new(b)),
        },
        // transform distributes over union (EXACT).
        Geom::Translate { child, t } => match *child {
            Geom::Union(a, b) => {
                log.push(Rewrite {
                    rule: "translate-distributes",
                    certificate: Certificate::Exact,
                });
                Geom::Union(
                    Box::new(Geom::Translate { child: a, t }),
                    Box::new(Geom::Translate { child: b, t }),
                )
            }
            Geom::Empty => {
                log.push(Rewrite {
                    rule: "translate-empty",
                    certificate: Certificate::Exact,
                });
                Geom::Empty
            }
            child => Geom::Translate {
                child: Box::new(child),
                t,
            },
        },
        other => other,
    }
}

// -- Shape grammar ----------------------------------------------------------

/// A shape-grammar production: `count` copies of `unit` spaced by `spacing`,
/// unioned (a rib / module pattern).
#[must_use]
pub fn linear_repeat(unit: &Geom, count: usize, spacing: [f64; 3]) -> Geom {
    if count == 0 {
        return Geom::Empty;
    }
    let mut acc = Geom::Empty;
    for i in 0..count {
        let f = i as f64;
        let copy = unit
            .clone()
            .translate([spacing[0] * f, spacing[1] * f, spacing[2] * f]);
        acc = if i == 0 {
            copy
        } else {
            Geom::Union(Box::new(acc), Box::new(copy))
        };
    }
    acc
}

/// A seeded stochastic derivation: a repeat of `1..=max_count` units, chosen
/// reproducibly from `seed`.
#[must_use]
pub fn stochastic_repeat(unit: &Geom, max_count: usize, spacing: [f64; 3], seed: u64) -> Geom {
    let count = if max_count == 0 {
        0
    } else {
        (seed % max_count as u64) as usize + 1
    };
    linear_repeat(unit, count, spacing)
}

// -- Parser (round-trip) ----------------------------------------------------

/// A parse failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    /// Unexpected end of input.
    UnexpectedEnd,
    /// An unexpected token.
    Unexpected(String),
    /// A malformed number.
    BadNumber(String),
}

/// Parse an s-expression program (round-trips with [`Geom::to_sexpr`]).
///
/// # Errors
/// [`ParseError`] on malformed input.
pub fn parse(s: &str) -> Result<Geom, ParseError> {
    let mut tokens = tokenize(s);
    tokens.reverse();
    let g = parse_expr(&mut tokens)?;
    if tokens.is_empty() {
        Ok(g)
    } else {
        Err(ParseError::Unexpected(tokens.pop().unwrap()))
    }
}

fn tokenize(s: &str) -> Vec<String> {
    s.replace('(', " ( ")
        .replace(')', " ) ")
        .split_whitespace()
        .map(str::to_string)
        .collect()
}

fn num(t: &str) -> Result<f64, ParseError> {
    t.parse::<f64>()
        .map_err(|_| ParseError::BadNumber(t.to_string()))
}

fn parse_expr(tokens: &mut Vec<String>) -> Result<Geom, ParseError> {
    let open = tokens.pop().ok_or(ParseError::UnexpectedEnd)?;
    if open != "(" {
        return Err(ParseError::Unexpected(open));
    }
    let head = tokens.pop().ok_or(ParseError::UnexpectedEnd)?;
    let g = match head.as_str() {
        "empty" => Geom::Empty,
        "sphere" => Geom::sphere(num(&pop(tokens)?)?),
        "cube" => Geom::cube(num(&pop(tokens)?)?),
        "union" => Geom::Union(Box::new(parse_expr(tokens)?), Box::new(parse_expr(tokens)?)),
        "intersect" => {
            Geom::Intersect(Box::new(parse_expr(tokens)?), Box::new(parse_expr(tokens)?))
        }
        "difference" => {
            Geom::Difference(Box::new(parse_expr(tokens)?), Box::new(parse_expr(tokens)?))
        }
        "offset" => {
            let child = Box::new(parse_expr(tokens)?);
            Geom::Offset {
                child,
                radius: num(&pop(tokens)?)?,
            }
        }
        "translate" => {
            let child = Box::new(parse_expr(tokens)?);
            let t = [
                num(&pop(tokens)?)?,
                num(&pop(tokens)?)?,
                num(&pop(tokens)?)?,
            ];
            Geom::Translate { child, t }
        }
        other => return Err(ParseError::Unexpected(other.to_string())),
    };
    let close = tokens.pop().ok_or(ParseError::UnexpectedEnd)?;
    if close != ")" {
        return Err(ParseError::Unexpected(close));
    }
    Ok(g)
}

fn pop(tokens: &mut Vec<String>) -> Result<String, ParseError> {
    tokens.pop().ok_or(ParseError::UnexpectedEnd)
}

// -- helpers ----------------------------------------------------------------

fn norm(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
}

fn cube_sdf(p: [f64; 3], half: f64) -> f64 {
    let q = [p[0].abs() - half, p[1].abs() - half, p[2].abs() - half];
    let outside = norm([q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)]);
    let inside = q[0].max(q[1]).max(q[2]).min(0.0);
    outside + inside
}

fn fmt(x: f64) -> String {
    // stable, round-trippable numeric print.
    let s = format!("{x}");
    if s.contains('.') || s.contains('e') || s.contains("inf") || s.contains("NaN") {
        s
    } else {
        format!("{s}.0")
    }
}

fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h
}
