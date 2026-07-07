//! Rigorous interval evaluation over fs-opt expression graphs: the
//! in-house prover behind Robust and Certification(Interval)
//! constraints. Conservative inclusion per node; anything it cannot
//! bound rigorously (division through zero, ln at nonpositive lo,
//! PDE/stochastic nodes) REFUSES with a teaching reason rather than
//! guessing. Rounding is to-nearest (outward-rounded arithmetic joins
//! with fs-ivl — CONTRACT no-claim); the conformance battery separates
//! the conservativeness law from that caveat with an fp-slack check.

use fs_opt::{Expr, NodeId, Problem};

/// A closed interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Iv {
    /// Lower end.
    pub lo: f64,
    /// Upper end.
    pub hi: f64,
}

impl Iv {
    fn point(v: f64) -> Iv {
        Iv { lo: v, hi: v }
    }

    fn add(self, o: Iv) -> Iv {
        Iv {
            lo: self.lo + o.lo,
            hi: self.hi + o.hi,
        }
    }

    fn sub(self, o: Iv) -> Iv {
        Iv {
            lo: self.lo - o.hi,
            hi: self.hi - o.lo,
        }
    }

    fn neg(self) -> Iv {
        Iv {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    fn mul(self, o: Iv) -> Iv {
        let c = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Iv {
            lo: c.iter().copied().fold(f64::INFINITY, f64::min),
            hi: c.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        }
    }

    fn min_iv(self, o: Iv) -> Iv {
        Iv {
            lo: self.lo.min(o.lo),
            hi: self.hi.min(o.hi),
        }
    }

    fn max_iv(self, o: Iv) -> Iv {
        Iv {
            lo: self.lo.max(o.lo),
            hi: self.hi.max(o.hi),
        }
    }

    fn abs(self) -> Iv {
        if self.lo >= 0.0 {
            self
        } else if self.hi <= 0.0 {
            self.neg()
        } else {
            Iv {
                lo: 0.0,
                hi: self.hi.max(-self.lo),
            }
        }
    }

    fn monotone(self, f: impl Fn(f64) -> f64) -> Iv {
        Iv {
            lo: f(self.lo),
            hi: f(self.hi),
        }
    }
}

/// Why the engine refused (honest gaps, not failures).
#[derive(Debug, Clone, PartialEq)]
pub enum IvalError {
    /// Division by an interval containing zero.
    DivThroughZero,
    /// `ln`/`sqrt` domain violated over the box.
    Domain {
        /// Which op.
        op: &'static str,
    },
    /// PDE/stochastic nodes have no interval semantics here.
    Unevaluable {
        /// Node id.
        node: u32,
    },
    /// Negative `powi` exponents are deferred.
    NegativePow,
    /// Binding count does not cover the variable's components.
    BadBindings,
}

impl core::fmt::Display for IvalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            IvalError::DivThroughZero => write!(
                f,
                "division by an interval containing zero cannot be bounded"
            ),
            IvalError::Domain { op } => {
                write!(f, "`{op}` leaves its domain somewhere in the box")
            }
            IvalError::Unevaluable { node } => write!(
                f,
                "node {node} (PDE/stochastic) has no interval semantics in this engine"
            ),
            IvalError::NegativePow => write!(f, "negative integer powers are deferred"),
            IvalError::BadBindings => {
                write!(f, "the box bindings do not cover the variable components")
            }
        }
    }
}

impl std::error::Error for IvalError {}

enum IvVal {
    S(Iv),
    V(Vec<Iv>),
}

/// Rigorously bound a scalar node over per-component boxes for the
/// problem's FIRST variable (the v1 host shape).
///
/// # Errors
/// [`IvalError`] naming the refusal.
pub fn interval_eval(
    problem: &Problem,
    node: NodeId,
    boxes: &[(f64, f64)],
) -> Result<Iv, IvalError> {
    match ival_at(problem, node, boxes)? {
        IvVal::S(iv) => Ok(iv),
        IvVal::V(_) => Err(IvalError::BadBindings),
    }
}

#[allow(clippy::too_many_lines)] // one inclusion rule per node kind
fn ival_at(problem: &Problem, node: NodeId, boxes: &[(f64, f64)]) -> Result<IvVal, IvalError> {
    let s = |v: IvVal| -> Iv {
        match v {
            IvVal::S(x) => x,
            IvVal::V(_) => unreachable!("builder enforced scalar shapes"),
        }
    };
    let out = match &problem.exprs[node.0 as usize] {
        Expr::Var(_) => IvVal::V(
            boxes
                .iter()
                .map(|&(lo, hi)| Iv { lo, hi })
                .collect::<Vec<Iv>>(),
        ),
        Expr::Component { of, index } => {
            let v = ival_at(problem, *of, boxes)?;
            match v {
                IvVal::V(xs) => IvVal::S(*xs.get(*index as usize).ok_or(IvalError::BadBindings)?),
                IvVal::S(_) => unreachable!("builder enforced vector shape"),
            }
        }
        Expr::Const { value, .. } => IvVal::S(Iv::point(*value)),
        Expr::Add(a, b) => {
            let (x, y) = (ival_at(problem, *a, boxes)?, ival_at(problem, *b, boxes)?);
            match (x, y) {
                (IvVal::S(p), IvVal::S(q)) => IvVal::S(p.add(q)),
                (IvVal::V(p), IvVal::V(q)) => {
                    IvVal::V(p.iter().zip(&q).map(|(u, w)| u.add(*w)).collect())
                }
                _ => unreachable!("builder enforced matching shapes"),
            }
        }
        Expr::Sub(a, b) => {
            let (x, y) = (ival_at(problem, *a, boxes)?, ival_at(problem, *b, boxes)?);
            match (x, y) {
                (IvVal::S(p), IvVal::S(q)) => IvVal::S(p.sub(q)),
                (IvVal::V(p), IvVal::V(q)) => {
                    IvVal::V(p.iter().zip(&q).map(|(u, w)| u.sub(*w)).collect())
                }
                _ => unreachable!("builder enforced matching shapes"),
            }
        }
        Expr::Mul(a, b) => {
            let (x, y) = (ival_at(problem, *a, boxes)?, ival_at(problem, *b, boxes)?);
            match (x, y) {
                (IvVal::S(p), IvVal::S(q)) => IvVal::S(p.mul(q)),
                (IvVal::S(p), IvVal::V(q)) | (IvVal::V(q), IvVal::S(p)) => {
                    IvVal::V(q.iter().map(|w| p.mul(*w)).collect())
                }
                _ => unreachable!("builder rejected vector*vector"),
            }
        }
        Expr::Div(a, b) => {
            let p = s(ival_at(problem, *a, boxes)?);
            let q = s(ival_at(problem, *b, boxes)?);
            if q.lo <= 0.0 && q.hi >= 0.0 {
                return Err(IvalError::DivThroughZero);
            }
            IvVal::S(p.mul(Iv {
                lo: 1.0 / q.hi,
                hi: 1.0 / q.lo,
            }))
        }
        Expr::Neg(a) => match ival_at(problem, *a, boxes)? {
            IvVal::S(p) => IvVal::S(p.neg()),
            IvVal::V(p) => IvVal::V(p.iter().map(|u| u.neg()).collect()),
        },
        Expr::Powi { base, exp } => {
            if *exp < 0 {
                return Err(IvalError::NegativePow);
            }
            let b = s(ival_at(problem, *base, boxes)?);
            let mut acc = Iv::point(1.0);
            for _ in 0..*exp {
                acc = acc.mul(b);
            }
            // Tighten even powers (dependency-aware square).
            if *exp % 2 == 0 && b.lo < 0.0 && b.hi > 0.0 {
                acc.lo = acc.lo.max(0.0);
            }
            IvVal::S(acc)
        }
        Expr::Sqrt(a) => {
            let p = s(ival_at(problem, *a, boxes)?);
            if p.lo < 0.0 {
                return Err(IvalError::Domain { op: "sqrt" });
            }
            IvVal::S(p.monotone(f64::sqrt))
        }
        Expr::Exp(a) => IvVal::S(s(ival_at(problem, *a, boxes)?).monotone(f64::exp)),
        Expr::Ln(a) => {
            let p = s(ival_at(problem, *a, boxes)?);
            if p.lo <= 0.0 {
                return Err(IvalError::Domain { op: "ln" });
            }
            IvVal::S(p.monotone(f64::ln))
        }
        Expr::Tanh(a) => IvVal::S(s(ival_at(problem, *a, boxes)?).monotone(f64::tanh)),
        Expr::Dot(a, b) => {
            let (x, y) = (ival_at(problem, *a, boxes)?, ival_at(problem, *b, boxes)?);
            match (x, y) {
                (IvVal::V(p), IvVal::V(q)) => {
                    let mut acc = Iv::point(0.0);
                    for (u, w) in p.iter().zip(&q) {
                        acc = acc.add(u.mul(*w));
                    }
                    IvVal::S(acc)
                }
                _ => unreachable!("builder enforced vectors"),
            }
        }
        Expr::NormSq(a) => match ival_at(problem, *a, boxes)? {
            IvVal::V(p) => {
                let mut acc = Iv::point(0.0);
                for u in &p {
                    let sq = u.mul(*u);
                    // x·x is nonnegative even when the box straddles 0.
                    acc = acc.add(Iv {
                        lo: sq.lo.max(0.0),
                        hi: sq.hi,
                    });
                }
                IvVal::S(acc)
            }
            IvVal::S(_) => unreachable!("builder enforced vector"),
        },
        Expr::Min(a, b) => {
            IvVal::S(s(ival_at(problem, *a, boxes)?).min_iv(s(ival_at(problem, *b, boxes)?)))
        }
        Expr::Max(a, b) => {
            IvVal::S(s(ival_at(problem, *a, boxes)?).max_iv(s(ival_at(problem, *b, boxes)?)))
        }
        Expr::Abs(a) => IvVal::S(s(ival_at(problem, *a, boxes)?).abs()),
        Expr::PdeResidual { .. }
        | Expr::Expectation { .. }
        | Expr::Cvar { .. }
        | Expr::Quantile { .. } => {
            return Err(IvalError::Unevaluable { node: node.0 });
        }
    };
    Ok(out)
}
