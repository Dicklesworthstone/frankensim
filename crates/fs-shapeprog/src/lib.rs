//! fs-shapeprog — generative geometry program synthesis. Layer: L2.
//!
//! Continuous optimization refines a topology; it rarely INVENTS a grammar.
//! This is the discrete-invention medium: a typed constructive-geometry DSL
//! with SDF semantics, a rewrite engine that SIMPLIFIES and CANONICALIZES
//! programs under geometric identities, and seeded shape-grammar derivation.
//!
//! The load-bearing safety property (the acceptance criterion): a rewrite
//! PRESERVES GEOMETRY within its declared compositional certificate, with
//! [`max_sdf_discrepancy`] as an independent finite-sample falsifier. Exact
//! identities are bit-equivalent under the interpreter's finite-input policy;
//! in particular, consecutive offsets are not reassociated because two rounded
//! subtractions need not equal one subtraction by a rounded sum.
//! Certified-approximate offset drops use a compositional outward bound.
//! [`canonical_hash`] gives equivalent programs one identity for archive/ledger
//! dedup. Deterministic; no dependencies.

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

    fn has_finite_parameters(&self) -> bool {
        match self {
            Geom::Empty => true,
            Geom::Primitive { size, .. } => size.is_finite(),
            Geom::Union(a, b) | Geom::Intersect(a, b) | Geom::Difference(a, b) => {
                a.has_finite_parameters() && b.has_finite_parameters()
            }
            Geom::Offset { child, radius } => radius.is_finite() && child.has_finite_parameters(),
            Geom::Translate { child, t } => {
                t.iter().all(|value| value.is_finite()) && child.has_finite_parameters()
            }
        }
    }

    fn admissible_sdf_evaluation(&self, p: [f64; 3]) -> AdmissibleSdf {
        let mut ignore_visit = || {};
        self.admissible_sdf_evaluation_with(p, &mut ignore_visit)
    }

    fn admissible_sdf_evaluation_with<F>(&self, p: [f64; 3], visit: &mut F) -> AdmissibleSdf
    where
        F: FnMut(),
    {
        visit();
        match self {
            Geom::Empty => AdmissibleSdf::StructuralEmpty,
            Geom::Primitive { .. } => AdmissibleSdf::finite(self.sdf(p)),
            Geom::Union(a, b) => match (
                a.admissible_sdf_evaluation_with(p, visit),
                b.admissible_sdf_evaluation_with(p, visit),
            ) {
                (AdmissibleSdf::StructuralEmpty, AdmissibleSdf::StructuralEmpty) => {
                    AdmissibleSdf::StructuralEmpty
                }
                (AdmissibleSdf::StructuralEmpty, AdmissibleSdf::Finite(value))
                | (AdmissibleSdf::Finite(value), AdmissibleSdf::StructuralEmpty) => {
                    AdmissibleSdf::Finite(value)
                }
                (AdmissibleSdf::Finite(left), AdmissibleSdf::Finite(right)) => {
                    AdmissibleSdf::Finite(left.min(right))
                }
                _ => AdmissibleSdf::Invalid,
            },
            Geom::Intersect(a, b) => match (
                a.admissible_sdf_evaluation_with(p, visit),
                b.admissible_sdf_evaluation_with(p, visit),
            ) {
                (AdmissibleSdf::StructuralEmpty, _) | (_, AdmissibleSdf::StructuralEmpty) => {
                    AdmissibleSdf::StructuralEmpty
                }
                (AdmissibleSdf::Finite(left), AdmissibleSdf::Finite(right)) => {
                    AdmissibleSdf::Finite(left.max(right))
                }
                _ => AdmissibleSdf::Invalid,
            },
            Geom::Difference(a, b) => match (
                a.admissible_sdf_evaluation_with(p, visit),
                b.admissible_sdf_evaluation_with(p, visit),
            ) {
                (AdmissibleSdf::StructuralEmpty, _) => AdmissibleSdf::StructuralEmpty,
                (AdmissibleSdf::Finite(left), AdmissibleSdf::StructuralEmpty) => {
                    AdmissibleSdf::Finite(left)
                }
                (AdmissibleSdf::Finite(left), AdmissibleSdf::Finite(right)) => {
                    AdmissibleSdf::Finite(left.max(-right))
                }
                _ => AdmissibleSdf::Invalid,
            },
            Geom::Offset { child, radius } => {
                match child.admissible_sdf_evaluation_with(p, visit) {
                    AdmissibleSdf::StructuralEmpty => AdmissibleSdf::StructuralEmpty,
                    AdmissibleSdf::Finite(value) => AdmissibleSdf::finite(value - radius),
                    AdmissibleSdf::Invalid => AdmissibleSdf::Invalid,
                }
            }
            Geom::Translate { child, t } => {
                let translated = [p[0] - t[0], p[1] - t[1], p[2] - t[2]];
                match child.admissible_sdf_evaluation_with(translated, visit) {
                    AdmissibleSdf::StructuralEmpty => AdmissibleSdf::StructuralEmpty,
                    AdmissibleSdf::Finite(value)
                        if translated.iter().all(|coordinate| coordinate.is_finite()) =>
                    {
                        AdmissibleSdf::Finite(value)
                    }
                    AdmissibleSdf::Finite(_) | AdmissibleSdf::Invalid => AdmissibleSdf::Invalid,
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum AdmissibleSdf {
    StructuralEmpty,
    Finite(f64),
    Invalid,
}

impl AdmissibleSdf {
    fn finite(value: f64) -> Self {
        if value.is_finite() {
            Self::Finite(value)
        } else {
            Self::Invalid
        }
    }
}

#[cfg(test)]
mod admissible_evaluation_tests {
    use super::{AdmissibleSdf, Geom};

    #[test]
    fn admissibility_walk_visits_each_node_exactly_once() {
        let mut program = Geom::sphere(1.0);
        for _ in 0..256 {
            program = program.translate([0.0, 0.0, 0.0]);
        }

        let mut visits = 0_usize;
        let result = program.admissible_sdf_evaluation_with([0.0, 0.0, 0.0], &mut || {
            visits += 1;
        });

        assert!(matches!(result, AdmissibleSdf::Finite(value) if value == -1.0));
        assert_eq!(
            visits,
            program.size(),
            "the evidence walk must do one local evaluation per AST node"
        );
    }
}

fn order(a: Geom, b: Geom) -> (Geom, Geom) {
    if a.to_sexpr() <= b.to_sexpr() {
        (a, b)
    } else {
        (b, a)
    }
}

/// An outward upper bound on the maximum `|SDF_a − SDF_b|` over the sample
/// points — the rewrite-safety check. Structurally empty SDFs agree at `+∞`;
/// invalid evidence or unrepresentable arithmetic returns `+∞` as a fail-closed
/// sentinel. Every non-structural intermediate branch value and translated
/// coordinate must be finite, so a finite selected Boolean root cannot mask
/// invalid evidence.
#[must_use]
pub fn max_sdf_discrepancy(a: &Geom, b: &Geom, samples: &[[f64; 3]]) -> f64 {
    if samples.is_empty() || !a.has_finite_parameters() || !b.has_finite_parameters() {
        return f64::INFINITY;
    }
    let mut worst = 0.0_f64;
    for &p in samples {
        if !p.iter().all(|value| value.is_finite()) {
            return f64::INFINITY;
        }
        let (da, db) = match (
            a.admissible_sdf_evaluation(p),
            b.admissible_sdf_evaluation(p),
        ) {
            (AdmissibleSdf::StructuralEmpty, AdmissibleSdf::StructuralEmpty) => continue,
            (AdmissibleSdf::Finite(da), AdmissibleSdf::Finite(db)) => (da, db),
            _ => return f64::INFINITY,
        };
        let Some(delta) = outward_abs_difference(da, db) else {
            return f64::INFINITY;
        };
        worst = worst.max(delta);
    }
    worst
}

fn outward_abs_difference(left: f64, right: f64) -> Option<f64> {
    let negated_right = -right;
    let difference = left + negated_right;
    if !difference.is_finite() {
        return None;
    }

    // Knuth two-sum recovers the exact residual of the rounded subtraction.
    // If that residual points away from zero, the rounded magnitude is a lower
    // bound and must advance one lattice point before it can serve as evidence.
    let right_virtual = difference - left;
    let residual = (left - (difference - right_virtual)) + (negated_right - right_virtual);
    let magnitude = difference.abs();
    let rounded_inward =
        (difference > 0.0 && residual > 0.0) || (difference < 0.0 && residual < 0.0);
    if rounded_inward {
        next_up_nonnegative(magnitude)
    } else {
        Some(magnitude)
    }
}

// -- Rewrite engine ---------------------------------------------------------

/// A rewrite's fidelity certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Certificate {
    /// The interpreter result is preserved bit-for-bit for every admitted
    /// finite evaluation.
    Exact,
    /// The SDF changes by at most the local nonnegative `bound`.
    Approximate {
        /// The certified error bound.
        bound: f64,
    },
}

/// One deterministic edge in the path from the program root to a rewrite.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RewritePathStep {
    /// Left child of a union.
    UnionLeft,
    /// Right child of a union.
    UnionRight,
    /// Left child of an intersection.
    IntersectLeft,
    /// Right child of an intersection.
    IntersectRight,
    /// Minuend (left child) of a difference.
    DifferenceLeft,
    /// Subtrahend (right child) of a difference.
    DifferenceRight,
    /// Child of an offset.
    OffsetChild,
    /// Child of a translation.
    TranslateChild,
}

/// The bound-algebra operation whose finite outward result was unavailable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundOperation {
    /// Sequential approximation errors must be added.
    Sequential,
    /// Independent Boolean alternatives use their maximum.
    Alternative,
    /// A parent operator scales a child bound by an absolute Lipschitz factor.
    Scale,
    /// The local rounded-subtraction envelope for dropping an offset.
    LocalOffset,
}

/// Why simplification transactionally returned the original program.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SimplifyRefusal {
    /// A programmatic caller supplied a non-finite geometry parameter.
    NonFiniteProgram,
    /// The tiny-offset tolerance was negative or non-finite.
    InvalidTolerance {
        /// Exact IEEE-754 bits of the rejected tolerance, preserving signed
        /// zero and NaN payload diagnostics without NaN equality ambiguity.
        tolerance_bits: u64,
    },
    /// A required outward bound had no finite `f64` representation.
    UnrepresentableBound {
        /// Zero-based rewrite pass.
        pass: usize,
        /// Program path at which propagation failed.
        path: Vec<RewritePathStep>,
        /// The failed composition operation.
        operation: BoundOperation,
    },
    /// The rewrite system did not reach a syntactic fixed point within the
    /// deterministic pass budget.
    PassLimitExceeded {
        /// Number of passes attempted.
        limit: usize,
    },
}

/// One applied rewrite.
#[derive(Debug, Clone, PartialEq)]
pub struct Rewrite {
    /// The rule name.
    pub rule: &'static str,
    /// Its certificate.
    pub certificate: Certificate,
    /// Zero-based pass in which the rewrite occurred.
    pub pass: usize,
    /// Deterministic root-to-node path before this rewrite.
    pub path: Vec<RewritePathStep>,
    /// Bound for the rewritten subtree after composing the work known at this
    /// log point. The pass-level trace is the final authority for root bounds.
    pub accumulated_bound: f64,
}

/// The outward bound contributed by one rewrite pass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PassBound {
    /// Zero-based pass index.
    pub pass: usize,
    /// Root-level uniform SDF error bound for this pass.
    pub bound: f64,
}

/// The result of simplifying a program.
#[derive(Debug, Clone, PartialEq)]
pub struct Simplified {
    /// The simplified program.
    pub program: Geom,
    /// The rewrites applied (in order).
    pub rewrites: Vec<Rewrite>,
    /// The composed uniform SDF error bound (`0` if all rewrites are exact).
    /// Sequential pass bounds are outward-added; this is not a maximum over
    /// local rewrites despite the legacy field name.
    pub max_error: f64,
    /// Index of the first approximate rewrite in [`Self::rewrites`].
    pub first_lossy_rewrite: Option<usize>,
    /// Per-pass root bounds from which [`Self::max_error`] is replayed.
    pub pass_bounds: Vec<PassBound>,
    /// A transactional refusal. On refusal, `program` is the original input,
    /// `rewrites` and `pass_bounds` are empty, and `max_error` is exactly zero.
    pub refusal: Option<SimplifyRefusal>,
}

impl Simplified {
    /// Whether simplification was refused and transactionally rolled back.
    #[must_use]
    pub fn is_refused(&self) -> bool {
        self.refusal.is_some()
    }
}

#[derive(Debug, Clone, Copy)]
struct ErrorBound(f64);

impl ErrorBound {
    const ZERO: Self = Self(0.0);

    fn checked(value: f64, operation: BoundOperation) -> Result<Self, BoundOperation> {
        if value.is_finite() && value >= 0.0 {
            Ok(Self(value))
        } else {
            Err(operation)
        }
    }

    fn sequential(self, next: Self) -> Result<Self, BoundOperation> {
        let sum = self.0 + next.0;
        if !sum.is_finite() {
            return Err(BoundOperation::Sequential);
        }

        // Knuth two-sum: `sum + residual` is the exact real addition when the
        // finite addition does not overflow. Advance one lattice point only
        // when round-to-nearest landed below the exact nonnegative sum.
        let next_virtual = sum - self.0;
        let residual = (self.0 - (sum - next_virtual)) + (next.0 - next_virtual);
        let upper = if residual > 0.0 {
            next_up_nonnegative(sum).ok_or(BoundOperation::Sequential)?
        } else {
            sum
        };
        Self::checked(upper, BoundOperation::Sequential)
    }

    fn alternative(self, other: Self) -> Self {
        Self(self.0.max(other.0))
    }

    fn scale(
        self,
        absolute_factor: f64,
        operation: BoundOperation,
    ) -> Result<Self, BoundOperation> {
        if !absolute_factor.is_finite() || absolute_factor < 0.0 {
            return Err(operation);
        }
        if self.0 == 0.0 || absolute_factor == 0.0 {
            return Ok(Self::ZERO);
        }
        if absolute_factor == 1.0 {
            return Ok(self);
        }

        let product = self.0 * absolute_factor;
        if !product.is_finite() {
            return Err(operation);
        }
        // Multiplication by two is an exact binary-exponent shift whenever the
        // result is finite (including the subnormal ladder). Other factors get
        // one conservative outward lattice step.
        let upper = if absolute_factor == 2.0 {
            product
        } else {
            next_up_nonnegative(product).ok_or(operation)?
        };
        Self::checked(upper, operation)
    }
}

#[derive(Debug)]
struct PassResult {
    program: Geom,
    error: ErrorBound,
}

#[derive(Debug)]
struct BoundFailure {
    path: Vec<RewritePathStep>,
    operation: BoundOperation,
}

const REWRITE_PASS_LIMIT: usize = 64;

/// Simplify a finite-parameter program to a true fixpoint under the registered
/// rewrites. A finite, nonnegative tolerance admits offsets with
/// `|radius| < tiny_offset_tol` for a certified approximate drop. Bounds are
/// composed according to their actual AST context: sequential effects add
/// outwardly, independent Boolean branches take a maximum, and parent
/// operators apply their declared absolute Lipschitz factor. The tolerance is
/// a local radius-admission threshold, not a promise that the composed
/// [`Simplified::max_error`] is no greater than that scalar.
///
/// Consecutive offsets are deliberately preserved. Reassociating their rounded
/// subtractions is not an exact interpreter identity. Any invalid input, bound
/// overflow, or pass-limit exhaustion returns the original program with an
/// explicit [`SimplifyRefusal`].
#[must_use]
pub fn simplify(g: &Geom, tiny_offset_tol: f64) -> Simplified {
    if !g.has_finite_parameters() {
        return refused(g, SimplifyRefusal::NonFiniteProgram);
    }
    if !tiny_offset_tol.is_finite() || tiny_offset_tol < 0.0 {
        return refused(
            g,
            SimplifyRefusal::InvalidTolerance {
                tolerance_bits: tiny_offset_tol.to_bits(),
            },
        );
    }

    let mut current = g.clone();
    let mut rewrites = Vec::new();
    let mut pass_bounds = Vec::new();
    let mut total = ErrorBound::ZERO;
    let mut converged = false;

    for pass in 0..REWRITE_PASS_LIMIT {
        let before = current.to_sexpr();
        let mut pass_rewrites = Vec::new();
        let result = match rewrite_pass(&current, tiny_offset_tol, pass, &[], &mut pass_rewrites) {
            Ok(result) => result,
            Err(failure) => {
                return refused(
                    g,
                    SimplifyRefusal::UnrepresentableBound {
                        pass,
                        path: failure.path,
                        operation: failure.operation,
                    },
                );
            }
        };

        if result.program.to_sexpr() == before {
            converged = true;
            break;
        }

        total = match total.sequential(result.error) {
            Ok(total) => total,
            Err(operation) => {
                return refused(
                    g,
                    SimplifyRefusal::UnrepresentableBound {
                        pass,
                        path: Vec::new(),
                        operation,
                    },
                );
            }
        };
        pass_bounds.push(PassBound {
            pass,
            bound: result.error.0,
        });
        rewrites.extend(pass_rewrites);
        current = result.program;
    }

    if !converged {
        return refused(
            g,
            SimplifyRefusal::PassLimitExceeded {
                limit: REWRITE_PASS_LIMIT,
            },
        );
    }

    let first_lossy_rewrite = rewrites
        .iter()
        .position(|rewrite| matches!(rewrite.certificate, Certificate::Approximate { .. }));
    Simplified {
        program: current,
        rewrites,
        max_error: total.0,
        first_lossy_rewrite,
        pass_bounds,
        refusal: None,
    }
}

fn refused(g: &Geom, refusal: SimplifyRefusal) -> Simplified {
    Simplified {
        program: g.clone(),
        rewrites: Vec::new(),
        max_error: 0.0,
        first_lossy_rewrite: None,
        pass_bounds: Vec::new(),
        refusal: Some(refusal),
    }
}

fn rewrite_pass(
    g: &Geom,
    tol: f64,
    pass: usize,
    path: &[RewritePathStep],
    log: &mut Vec<Rewrite>,
) -> Result<PassResult, BoundFailure> {
    match g {
        Geom::Empty | Geom::Primitive { .. } => Ok(PassResult {
            program: g.clone(),
            error: ErrorBound::ZERO,
        }),
        Geom::Union(a, b) => {
            if matches!(a.as_ref(), Geom::Empty) {
                let child = rewrite_child(b, tol, pass, path, RewritePathStep::UnionRight, log)?;
                record_exact(log, pass, path, "union-identity", child.error);
                return Ok(child);
            }
            if matches!(b.as_ref(), Geom::Empty) {
                let child = rewrite_child(a, tol, pass, path, RewritePathStep::UnionLeft, log)?;
                record_exact(log, pass, path, "union-identity", child.error);
                return Ok(child);
            }

            let left = rewrite_child(a, tol, pass, path, RewritePathStep::UnionLeft, log)?;
            let right = rewrite_child(b, tol, pass, path, RewritePathStep::UnionRight, log)?;
            let inherited = left.error.alternative(right.error);
            Ok(apply_root_rule(
                Geom::Union(Box::new(left.program), Box::new(right.program)),
                inherited,
                pass,
                path,
                log,
            ))
        }
        Geom::Intersect(a, b) => {
            if matches!(a.as_ref(), Geom::Empty) || matches!(b.as_ref(), Geom::Empty) {
                record_exact(log, pass, path, "intersect-empty", ErrorBound::ZERO);
                return Ok(PassResult {
                    program: Geom::Empty,
                    error: ErrorBound::ZERO,
                });
            }

            let left = rewrite_child(a, tol, pass, path, RewritePathStep::IntersectLeft, log)?;
            let right = rewrite_child(b, tol, pass, path, RewritePathStep::IntersectRight, log)?;
            let inherited = left.error.alternative(right.error);
            Ok(apply_root_rule(
                Geom::Intersect(Box::new(left.program), Box::new(right.program)),
                inherited,
                pass,
                path,
                log,
            ))
        }
        Geom::Difference(a, b) => {
            if matches!(b.as_ref(), Geom::Empty) {
                let child =
                    rewrite_child(a, tol, pass, path, RewritePathStep::DifferenceLeft, log)?;
                record_exact(log, pass, path, "difference-identity", child.error);
                return Ok(child);
            }

            let left = rewrite_child(a, tol, pass, path, RewritePathStep::DifferenceLeft, log)?;
            let right = rewrite_child(b, tol, pass, path, RewritePathStep::DifferenceRight, log)?;
            let inherited = left.error.alternative(right.error);
            Ok(apply_root_rule(
                Geom::Difference(Box::new(left.program), Box::new(right.program)),
                inherited,
                pass,
                path,
                log,
            ))
        }
        Geom::Offset { child, radius } => {
            if matches!(child.as_ref(), Geom::Empty) {
                record_exact(log, pass, path, "offset-empty", ErrorBound::ZERO);
                return Ok(PassResult {
                    program: Geom::Empty,
                    error: ErrorBound::ZERO,
                });
            }

            if radius.abs() < tol {
                // For x and correctly-rounded z = RN(x-r), x itself is a
                // representable candidate for the exact x-r. Nearest rounding
                // gives |z-(x-r)| <= |r|, hence |z-x| <= 2|r|. The factor two
                // is necessary at rounding-cell threshold neighbours.
                let local = ErrorBound::checked(radius.abs(), BoundOperation::LocalOffset)
                    .and_then(|bound| bound.scale(2.0, BoundOperation::LocalOffset))
                    .map_err(|operation| BoundFailure {
                        path: path.to_vec(),
                        operation,
                    })?;
                let log_index = log.len();
                log.push(Rewrite {
                    rule: "drop-tiny-offset",
                    certificate: Certificate::Approximate { bound: local.0 },
                    pass,
                    path: path.to_vec(),
                    accumulated_bound: local.0,
                });

                // Apply the root drop conceptually before simplifying its
                // child. The two errors are sequential on one evaluation path.
                let child =
                    rewrite_child(child, tol, pass, path, RewritePathStep::OffsetChild, log)?;
                let error = local
                    .sequential(child.error)
                    .map_err(|operation| BoundFailure {
                        path: path.to_vec(),
                        operation,
                    })?;
                log[log_index].accumulated_bound = error.0;
                Ok(PassResult {
                    program: child.program,
                    error,
                })
            } else {
                let child =
                    rewrite_child(child, tol, pass, path, RewritePathStep::OffsetChild, log)?;
                // Rounded translation is not globally Lipschitz on the float
                // lattice: tiny adjacent inputs can straddle a much larger
                // rounding-cell boundary after the same shift. Propagate the
                // real affine map with factor one, then add two nearest-rounding
                // envelopes. For z=RN(x-r), x is a representable candidate, so
                // |z-(x-r)| <= |r|; applying this to both child values gives
                // |RN(x-r)-RN(y-r)| <= |x-y| + 2|r|. When the child bound is
                // zero, both evaluations are identical and no envelope is
                // needed. Crucially, offset nodes are never reassociated.
                let error = if child.error.0 == 0.0 {
                    ErrorBound::ZERO
                } else {
                    let propagated =
                        child
                            .error
                            .scale(1.0, BoundOperation::Scale)
                            .map_err(|operation| BoundFailure {
                                path: path.to_vec(),
                                operation,
                            })?;
                    let rounding = ErrorBound::checked(radius.abs(), BoundOperation::Scale)
                        .and_then(|bound| bound.scale(2.0, BoundOperation::Scale))
                        .map_err(|operation| BoundFailure {
                            path: path.to_vec(),
                            operation,
                        })?;
                    propagated
                        .sequential(rounding)
                        .map_err(|operation| BoundFailure {
                            path: path.to_vec(),
                            operation,
                        })?
                };
                Ok(PassResult {
                    program: Geom::Offset {
                        child: Box::new(child.program),
                        radius: *radius,
                    },
                    error,
                })
            }
        }
        Geom::Translate { child, t } => {
            if matches!(child.as_ref(), Geom::Empty) {
                record_exact(log, pass, path, "translate-empty", ErrorBound::ZERO);
                return Ok(PassResult {
                    program: Geom::Empty,
                    error: ErrorBound::ZERO,
                });
            }
            let child =
                rewrite_child(child, tol, pass, path, RewritePathStep::TranslateChild, log)?;
            let error = child
                .error
                .scale(1.0, BoundOperation::Scale)
                .map_err(|operation| BoundFailure {
                    path: path.to_vec(),
                    operation,
                })?;
            Ok(apply_root_rule(
                Geom::Translate {
                    child: Box::new(child.program),
                    t: *t,
                },
                error,
                pass,
                path,
                log,
            ))
        }
    }
}

fn rewrite_child(
    child: &Geom,
    tol: f64,
    pass: usize,
    path: &[RewritePathStep],
    step: RewritePathStep,
    log: &mut Vec<Rewrite>,
) -> Result<PassResult, BoundFailure> {
    let mut child_path = Vec::with_capacity(path.len() + 1);
    child_path.extend_from_slice(path);
    child_path.push(step);
    rewrite_pass(child, tol, pass, &child_path, log)
}

/// Apply a single exact root-level rewrite after its children are simplified.
fn apply_root_rule(
    g: Geom,
    inherited: ErrorBound,
    pass: usize,
    path: &[RewritePathStep],
    log: &mut Vec<Rewrite>,
) -> PassResult {
    let program = match g {
        Geom::Union(a, b) => match (*a, *b) {
            (Geom::Empty, x) | (x, Geom::Empty) => {
                record_exact(log, pass, path, "union-identity", inherited);
                x
            }
            (a, b) => Geom::Union(Box::new(a), Box::new(b)),
        },
        Geom::Difference(a, b) => match *b {
            Geom::Empty => {
                record_exact(log, pass, path, "difference-identity", inherited);
                *a
            }
            b => Geom::Difference(a, Box::new(b)),
        },
        Geom::Intersect(a, b) => match (*a, *b) {
            (Geom::Empty, _) | (_, Geom::Empty) => {
                record_exact(log, pass, path, "intersect-empty", inherited);
                Geom::Empty
            }
            (a, b) => Geom::Intersect(Box::new(a), Box::new(b)),
        },
        Geom::Translate { child, t } => match *child {
            Geom::Union(a, b) => {
                record_exact(log, pass, path, "translate-distributes", inherited);
                Geom::Union(
                    Box::new(Geom::Translate { child: a, t }),
                    Box::new(Geom::Translate { child: b, t }),
                )
            }
            Geom::Empty => {
                record_exact(log, pass, path, "translate-empty", inherited);
                Geom::Empty
            }
            child => Geom::Translate {
                child: Box::new(child),
                t,
            },
        },
        other => other,
    };
    PassResult {
        program,
        error: inherited,
    }
}

fn record_exact(
    log: &mut Vec<Rewrite>,
    pass: usize,
    path: &[RewritePathStep],
    rule: &'static str,
    accumulated: ErrorBound,
) {
    log.push(Rewrite {
        rule,
        certificate: Certificate::Exact,
        pass,
        path: path.to_vec(),
        accumulated_bound: accumulated.0,
    });
}

fn next_up_nonnegative(value: f64) -> Option<f64> {
    debug_assert!(value.is_finite() && value >= 0.0);
    let next = if value == 0.0 {
        f64::from_bits(1)
    } else {
        f64::from_bits(value.to_bits() + 1)
    };
    next.is_finite().then_some(next)
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

/// Parse an s-expression program (round-trips with [`Geom::to_sexpr`] for
/// finite-parameter programs).
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
    let value = t
        .parse::<f64>()
        .map_err(|_| ParseError::BadNumber(t.to_string()))?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ParseError::BadNumber(t.to_string()))
    }
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
