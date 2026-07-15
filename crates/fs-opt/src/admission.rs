//! The single versioned problem-admission validator (bead sj31i.48).
//!
//! Every path that can put structure into a [`Problem`] — the
//! incremental [`crate::ProblemBuilder`], the wire parser (which
//! rebuilds through the builder), and any future migration — validates
//! through the SAME leaf rules in this module, so builder acceptance
//! and [`Problem::admit`] re-acceptance cannot drift. Admission is the
//! defense-in-depth chokepoint: it re-derives every node's
//! shape/dimension/class from the expression list, re-checks every
//! leaf policy and cap, proves reference validity and acyclicity from
//! the arena ordering, and only then mints the problem's
//! [`ProblemSemanticId`]. Violations come back as a COMPLETE report in
//! deterministic order (section by section, ascending index), never a
//! first-error-only refusal.

use crate::ir::Class;
use crate::ir::{
    BilevelRef, Expr, Manifold, NodeId, OptError, Problem, ProblemTag, Shape, Variable, children,
    own_class,
};
use crate::serial::{LegacyProblemHash, ProblemSemanticId};
use fs_qty::Dims;

/// Version of the admission schema: bump when a rule, cap default, or
/// the semantic-id preimage changes meaning.
pub const ADMISSION_SCHEMA_VERSION: u32 = 2;

/// Versioned per-item and aggregate admission caps. Construct via
/// [`AdmissionCaps::default`] and override individual fields; the
/// struct is non-exhaustive so new caps can join without breaking
/// callers.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmissionCaps {
    /// Maximum declared variables.
    pub max_vars: u32,
    /// Maximum expression nodes (also bounds every `NodeId`).
    pub max_nodes: u32,
    /// Maximum objectives.
    pub max_objectives: u32,
    /// Maximum constraints.
    pub max_constraints: u32,
    /// Maximum structure tags.
    pub max_tags: u32,
    /// Maximum bytes in a variable/constraint name.
    pub max_name_bytes: u64,
    /// Maximum bytes in a study/UQ-config identifier.
    pub max_string_bytes: u64,
    /// Maximum multi-fidelity levels.
    pub max_fidelity_levels: u32,
    /// Maximum point storage of ONE variable.
    pub max_point_dim: u32,
    /// Maximum SUMMED point storage across all variables.
    pub max_total_point_storage: u64,
    /// Maximum expression depth. This bounds recursive consumers even
    /// when the graph is a valid, acyclic chain.
    pub max_graph_depth: u32,
    /// Conservative aggregate byte budget for retained names, expression
    /// strings, intern keys, and their canonical wire spellings.
    pub max_total_retained_bytes: u64,
    /// Maximum exact input artifact size accepted by the wire parser.
    pub max_wire_bytes: u64,
}

impl AdmissionCaps {
    /// Compatibility spelling for callers pinned to the original cap
    /// constant. Validation itself is governed by schema v2.
    pub const V1: AdmissionCaps = AdmissionCaps {
        max_vars: 4096,
        max_nodes: 1 << 20,
        max_objectives: 1024,
        max_constraints: 1 << 16,
        max_tags: 1024,
        max_name_bytes: 4096,
        max_string_bytes: 4096,
        max_fidelity_levels: 1024,
        max_point_dim: 1 << 24,
        max_total_point_storage: 1 << 32,
        max_graph_depth: 256,
        max_total_retained_bytes: 1 << 28,
        max_wire_bytes: 1 << 26,
    };

    /// The v2 cap schedule.
    pub const V2: AdmissionCaps = AdmissionCaps {
        max_vars: 4096,
        max_nodes: 1 << 20,
        max_objectives: 1024,
        max_constraints: 1 << 16,
        max_tags: 1024,
        max_name_bytes: 4096,
        max_string_bytes: 4096,
        max_fidelity_levels: 1024,
        max_point_dim: 1 << 24,
        max_total_point_storage: 1 << 32,
        max_graph_depth: 256,
        max_total_retained_bytes: 1 << 28,
        max_wire_bytes: 1 << 26,
    };
}

impl Default for AdmissionCaps {
    fn default() -> Self {
        AdmissionCaps::V2
    }
}

/// One admission violation, locating the offending section and index.
#[derive(Debug, Clone, PartialEq)]
pub enum AdmissionViolation {
    /// An aggregate cap or metadata-vector alignment failure.
    Aggregate {
        /// What failed.
        what: &'static str,
        /// Observed count.
        count: u64,
        /// The cap or expected count.
        cap: u64,
    },
    /// A variable failed leaf validation.
    Var {
        /// Variable index.
        index: u32,
        /// The teaching error.
        error: OptError,
    },
    /// An expression node failed re-derivation or leaf validation.
    Node {
        /// Node index.
        index: u32,
        /// The teaching error.
        error: OptError,
    },
    /// A node references a child at or after itself — the arena order
    /// proof of acyclicity/reference validity fails.
    ChildOrder {
        /// The referencing node.
        node: u32,
        /// The out-of-order child id.
        child: u32,
    },
    /// A cached shape/dimension/class disagrees with re-derivation.
    CacheMismatch {
        /// The node.
        node: u32,
        /// Which cache.
        what: &'static str,
    },
    /// An objective failed validation.
    Objective {
        /// Objective index.
        index: u32,
        /// The teaching error.
        error: OptError,
    },
    /// A constraint failed validation.
    Constraint {
        /// Constraint index.
        index: u32,
        /// The teaching error.
        error: OptError,
    },
    /// A tag failed validation.
    Tag {
        /// Tag index.
        index: u32,
        /// The teaching error.
        error: OptError,
    },
}

impl core::fmt::Display for AdmissionViolation {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            AdmissionViolation::Aggregate { what, count, cap } => {
                write!(f, "aggregate: {what} = {count} (limit {cap})")
            }
            AdmissionViolation::Var { index, error } => write!(f, "var {index}: {error}"),
            AdmissionViolation::Node { index, error } => write!(f, "node {index}: {error}"),
            AdmissionViolation::ChildOrder { node, child } => write!(
                f,
                "node {node} references child {child} at or after itself; the arena \
                 ordering proof of acyclicity fails"
            ),
            AdmissionViolation::CacheMismatch { node, what } => write!(
                f,
                "node {node}: cached {what} disagrees with re-derivation from the \
                 expression list"
            ),
            AdmissionViolation::Objective { index, error } => {
                write!(f, "objective {index}: {error}")
            }
            AdmissionViolation::Constraint { index, error } => {
                write!(f, "constraint {index}: {error}")
            }
            AdmissionViolation::Tag { index, error } => write!(f, "tag {index}: {error}"),
        }
    }
}

/// The complete, deterministically ordered rejection report.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmissionReport {
    schema_version: u32,
    violations: Vec<AdmissionViolation>,
}

impl AdmissionReport {
    /// Admission schema the report was produced under.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// The violations, in deterministic section/index order.
    #[must_use]
    pub fn violations(&self) -> &[AdmissionViolation] {
        &self.violations
    }
}

impl core::fmt::Display for AdmissionReport {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "problem admission (schema v{}) refused with {} violation(s):",
            self.schema_version,
            self.violations.len()
        )?;
        for v in &self.violations {
            write!(f, "\n  - {v}")?;
        }
        Ok(())
    }
}

impl std::error::Error for AdmissionReport {}

/// Evidence that a [`Problem`] passed the common admission validator.
/// Fields are sealed; the only constructor is [`admit_with_caps`], so
/// holding a `ProblemAdmission` means the checks actually ran.
#[derive(Debug, Clone, PartialEq)]
pub struct ProblemAdmission {
    schema_version: u32,
    semantic_id: ProblemSemanticId,
    var_count: u32,
    node_count: u32,
    objective_count: u32,
    constraint_count: u32,
    total_point_storage: u64,
    total_retained_bytes: u64,
    quarantined_legacy_identities: Vec<LegacyProblemHash>,
}

impl ProblemAdmission {
    /// Admission schema version the checks ran under.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// The problem's semantic identity over its normalized admitted
    /// meaning (domain-separated BLAKE3 of the canonical v3 body).
    #[must_use]
    pub fn semantic_id(&self) -> ProblemSemanticId {
        self.semantic_id
    }

    /// Declared variable count.
    #[must_use]
    pub fn var_count(&self) -> u32 {
        self.var_count
    }

    /// Expression node count.
    #[must_use]
    pub fn node_count(&self) -> u32 {
        self.node_count
    }

    /// Objective count.
    #[must_use]
    pub fn objective_count(&self) -> u32 {
        self.objective_count
    }

    /// Constraint count.
    #[must_use]
    pub fn constraint_count(&self) -> u32 {
        self.constraint_count
    }

    /// Summed point storage across all variables.
    #[must_use]
    pub fn total_point_storage(&self) -> u64 {
        self.total_point_storage
    }

    /// Conservative retained/canonical byte accounting for the admitted
    /// problem. The same accounting is enforced incrementally by the
    /// builder before it retains a new item.
    #[must_use]
    pub fn total_retained_bytes(&self) -> u64 {
        self.total_retained_bytes
    }

    /// Legacy FNV-64 identities carried by this problem (bilevel
    /// references parsed from historical artifacts). Inspectable
    /// provenance only — they confer NO execution or certificate
    /// authority, and admission never upgrades them.
    #[must_use]
    pub fn quarantined_legacy_identities(&self) -> &[LegacyProblemHash] {
        &self.quarantined_legacy_identities
    }
}

/// Validate a manifold descriptor (leaf rule shared with the builder).
///
/// # Errors
/// [`OptError::ManifoldInvalid`] / [`OptError::CapExceeded`].
pub(crate) fn validate_manifold(m: &Manifold, caps: &AdmissionCaps) -> Result<(), OptError> {
    match *m {
        Manifold::Rn { dim } => {
            if dim == 0 {
                return Err(OptError::ManifoldInvalid {
                    what: "Rn { dim: 0 } declares a variable with no point storage; a \
                           zero-dimensional variable cannot bind a point (declare a \
                           constant instead)"
                        .to_string(),
                });
            }
        }
        Manifold::Sphere { ambient } => {
            if ambient < 2 {
                return Err(OptError::ManifoldInvalid {
                    what: format!(
                        "Sphere {{ ambient: {ambient} }} is degenerate; the unit sphere \
                         needs ambient >= 2 (its tangent space is ambient - 1 \
                         dimensional)"
                    ),
                });
            }
        }
        Manifold::So3 => {}
        Manifold::Stiefel { n, p } => {
            if p == 0 || p > n {
                return Err(OptError::ManifoldInvalid {
                    what: format!(
                        "Stiefel {{ n: {n}, p: {p} }} needs 1 <= p <= n (orthonormal \
                         p-frames in n dimensions)"
                    ),
                });
            }
        }
    }
    let point = m.point_dim().ok_or_else(|| OptError::ManifoldInvalid {
        what: format!("{m:?} point storage overflows the u32 domain"),
    })?;
    m.tangent_dim().ok_or_else(|| OptError::ManifoldInvalid {
        what: format!("{m:?} tangent dimension is not representable"),
    })?;
    if point > caps.max_point_dim {
        return Err(OptError::CapExceeded {
            what: "variable point storage",
            count: u64::from(point),
            cap: u64::from(caps.max_point_dim),
        });
    }
    Ok(())
}

/// Validate a diagnostic name (leaf rule shared with the builder).
pub(crate) fn validate_name(
    what: &'static str,
    s: &str,
    caps: &AdmissionCaps,
) -> Result<(), OptError> {
    if s.len() as u64 > caps.max_name_bytes {
        return Err(OptError::CapExceeded {
            what,
            count: s.len() as u64,
            cap: caps.max_name_bytes,
        });
    }
    Ok(())
}

/// Validate a study/UQ-config identifier (leaf rule).
pub(crate) fn validate_string(
    what: &'static str,
    s: &str,
    caps: &AdmissionCaps,
) -> Result<(), OptError> {
    if s.len() as u64 > caps.max_string_bytes {
        return Err(OptError::CapExceeded {
            what,
            count: s.len() as u64,
            cap: caps.max_string_bytes,
        });
    }
    Ok(())
}

fn escaped_wire_bytes(s: &str) -> u64 {
    s.as_bytes().iter().fold(0u64, |count, byte| {
        count.saturating_add(if byte.is_ascii_graphic() && *byte != b'%' {
            1
        } else {
            3
        })
    })
}

/// Conservative retained/wire byte charge for a diagnostic name. The
/// fixed allowance covers its directive, ids, separators, and container
/// metadata; raw and escaped bytes cover both the owned value and its
/// canonical spelling.
pub(crate) fn name_retained_bytes(s: &str) -> u64 {
    64u64
        .saturating_add(s.len() as u64)
        .saturating_add(escaped_wire_bytes(s))
}

/// Conservative retained/wire byte charge for one expression, including
/// the owned expression, its intern key, and canonical directive.
pub(crate) fn expr_retained_bytes(e: &Expr) -> u64 {
    let string = match e {
        Expr::PdeResidual { study, .. } => study.as_str(),
        Expr::Expectation { uq_config, .. }
        | Expr::Cvar { uq_config, .. }
        | Expr::Quantile { uq_config, .. } => uq_config.as_str(),
        _ => "",
    };
    192u64
        .saturating_add(string.len() as u64)
        .saturating_add(escaped_wire_bytes(string))
}

pub(crate) fn checked_retained_total(
    current: u64,
    addition: u64,
    caps: &AdmissionCaps,
) -> Result<u64, OptError> {
    let total = current.saturating_add(addition);
    if total > caps.max_total_retained_bytes {
        return Err(OptError::CapExceeded {
            what: "total retained/canonical bytes",
            count: total,
            cap: caps.max_total_retained_bytes,
        });
    }
    Ok(total)
}

/// Derive and enforce the depth of one expression from previously
/// admitted child depths. Arena order makes this a complete acyclicity-
/// compatible depth calculation.
pub(crate) fn derive_depth(
    e: &Expr,
    lookup: &dyn Fn(NodeId) -> Option<u32>,
    caps: &AdmissionCaps,
) -> Result<u32, OptError> {
    let mut depth = 1u32;
    for child in children(e) {
        let child_depth = lookup(child).ok_or(OptError::UnknownNode { id: child.0 })?;
        depth = depth.max(child_depth.checked_add(1).unwrap_or(u32::MAX));
    }
    if depth > caps.max_graph_depth {
        return Err(OptError::CapExceeded {
            what: "graph depth",
            count: u64::from(depth),
            cap: u64::from(caps.max_graph_depth),
        });
    }
    Ok(depth)
}

/// Validate an objective weight: FINITE and NONNEGATIVE (`Sense`
/// already carries direction), refusing `-0.0` so one meaning cannot
/// have two bit-pattern wire identities.
pub(crate) fn validate_weight(weight: f64) -> Result<(), OptError> {
    if !weight.is_finite() {
        return Err(OptError::NonFinite {
            what: "objective weight",
            bits: weight.to_bits(),
        });
    }
    if weight < 0.0 || (weight == 0.0 && weight.is_sign_negative()) {
        return Err(OptError::BadParam {
            what: "objective weight (finite nonnegative; Sense carries direction)",
            value: weight,
        });
    }
    Ok(())
}

/// Validate a structure tag (leaf rule shared with the builder).
pub(crate) fn validate_tag(tag: &ProblemTag, caps: &AdmissionCaps) -> Result<(), OptError> {
    match tag {
        ProblemTag::MultiFidelity { levels } => {
            if *levels == 0 {
                return Err(OptError::BadParam {
                    what: "multi_fidelity levels (must be nonzero)",
                    value: 0.0,
                });
            }
            if *levels > caps.max_fidelity_levels {
                return Err(OptError::CapExceeded {
                    what: "multi_fidelity levels",
                    count: u64::from(*levels),
                    cap: u64::from(caps.max_fidelity_levels),
                });
            }
        }
        ProblemTag::ChanceConstrained { prob } => {
            if !(prob.is_finite() && *prob > 0.0 && *prob < 1.0) {
                return Err(OptError::BadParam {
                    what: "chance probability (finite, strictly inside (0, 1))",
                    value: *prob,
                });
            }
        }
        ProblemTag::Bilevel { .. } => {}
    }
    Ok(())
}

fn dims_add_checked(a: Dims, b: Dims, op: &'static str) -> Result<Dims, OptError> {
    let mut out = [0i8; 6];
    for (i, slot) in out.iter_mut().enumerate() {
        let sum = i16::from(a.0[i]) + i16::from(b.0[i]);
        *slot = i8::try_from(sum).map_err(|_| OptError::DimSumOverflow {
            op,
            left: a.0,
            right: b.0,
        })?;
    }
    Ok(Dims(out))
}

fn dims_sub_checked(a: Dims, b: Dims, op: &'static str) -> Result<Dims, OptError> {
    let mut out = [0i8; 6];
    for (i, slot) in out.iter_mut().enumerate() {
        let diff = i16::from(a.0[i]) - i16::from(b.0[i]);
        *slot = i8::try_from(diff).map_err(|_| OptError::DimSumOverflow {
            op,
            left: a.0,
            right: b.0,
        })?;
    }
    Ok(Dims(out))
}

type NodeInfo = (Shape, Dims, Class);

/// THE expression rule: derive (shape, dims, class) for a candidate
/// expression against already-admitted context, refusing with the same
/// teaching errors the incremental builder reports. Both the builder
/// and [`admit_with_caps`] call this — one validator, two entry
/// points.
#[allow(clippy::too_many_lines)] // one auditable arm per node kind — this IS the rule table
pub(crate) fn derive_expr(
    e: &Expr,
    lookup: &dyn Fn(NodeId) -> Option<NodeInfo>,
    vars: &[Variable],
    caps: &AdmissionCaps,
) -> Result<NodeInfo, OptError> {
    let get = |n: NodeId| lookup(n).ok_or(OptError::UnknownNode { id: n.0 });
    let scalar = |op: &'static str, info: NodeInfo| -> Result<NodeInfo, OptError> {
        match info.0 {
            Shape::Scalar => Ok(info),
            v @ Shape::Vector(_) => Err(OptError::ShapeMismatch {
                op,
                left: v,
                right: Shape::Scalar,
            }),
        }
    };
    let same_shape_dims = |op: &'static str, a: NodeId, b: NodeId| -> Result<NodeInfo, OptError> {
        let (ia, ib) = (get(a)?, get(b)?);
        if ia.0 != ib.0 {
            return Err(OptError::ShapeMismatch {
                op,
                left: ia.0,
                right: ib.0,
            });
        }
        if ia.1 != ib.1 {
            return Err(OptError::DimMismatch {
                op,
                left: ia.1.0,
                right: ib.1.0,
            });
        }
        Ok(ia)
    };
    let (shape, dims) = match e {
        Expr::Var(v) => {
            let var = vars
                .get(v.0 as usize)
                .ok_or(OptError::UnknownVar { id: v.0 })?;
            let point = var
                .manifold
                .point_dim()
                .ok_or_else(|| OptError::ManifoldInvalid {
                    what: format!("{:?} point storage overflows the u32 domain", var.manifold),
                })?;
            (Shape::Vector(point), var.dims)
        }
        Expr::Component { of, index } => {
            let info = get(*of)?;
            match info.0 {
                Shape::Vector(len) if *index < len => (Shape::Scalar, info.1),
                Shape::Vector(len) => return Err(OptError::IndexOut { index: *index, len }),
                Shape::Scalar => {
                    // Diagnostic-only arithmetic: the "required" shape
                    // saturates so component(scalar, u32::MAX) reports
                    // instead of wrapping/panicking (fail-closed on
                    // arbitrary indices).
                    return Err(OptError::ShapeMismatch {
                        op: "component",
                        left: Shape::Scalar,
                        right: Shape::Vector(index.saturating_add(1)),
                    });
                }
            }
        }
        Expr::Const { value, dims } => {
            if !value.is_finite() {
                return Err(OptError::NonFinite {
                    what: "constant value",
                    bits: value.to_bits(),
                });
            }
            (Shape::Scalar, *dims)
        }
        Expr::Add(a, b) => {
            let info = same_shape_dims("add", *a, *b)?;
            (info.0, info.1)
        }
        Expr::Sub(a, b) => {
            let info = same_shape_dims("sub", *a, *b)?;
            (info.0, info.1)
        }
        Expr::Min(a, b) => {
            let info = same_shape_dims("min", *a, *b)?;
            // Min/Max are SCALAR-only end to end: both the strict and
            // the interval evaluator implement only the scalar case, so
            // admitting vectors here would admit a program that panics
            // downstream. Enforce the evaluators' contract at the
            // shared leaf rule (batch-verify High #3).
            if info.0 != Shape::Scalar {
                return Err(OptError::ShapeMismatch {
                    op: "min",
                    left: info.0,
                    right: Shape::Scalar,
                });
            }
            (Shape::Scalar, info.1)
        }
        Expr::Max(a, b) => {
            let info = same_shape_dims("max", *a, *b)?;
            if info.0 != Shape::Scalar {
                return Err(OptError::ShapeMismatch {
                    op: "max",
                    left: info.0,
                    right: Shape::Scalar,
                });
            }
            (Shape::Scalar, info.1)
        }
        Expr::Mul(a, b) => {
            let (ia, ib) = (get(*a)?, get(*b)?);
            let dims = dims_add_checked(ia.1, ib.1, "mul")?;
            let shape = match (ia.0, ib.0) {
                (Shape::Scalar, Shape::Scalar) => Shape::Scalar,
                (Shape::Scalar, Shape::Vector(n)) | (Shape::Vector(n), Shape::Scalar) => {
                    Shape::Vector(n)
                }
                (l, r) => {
                    return Err(OptError::ShapeMismatch {
                        op: "mul",
                        left: l,
                        right: r,
                    });
                }
            };
            (shape, dims)
        }
        Expr::Div(a, b) => {
            let ia = scalar("div", get(*a)?)?;
            let ib = scalar("div", get(*b)?)?;
            (Shape::Scalar, dims_sub_checked(ia.1, ib.1, "div")?)
        }
        Expr::Neg(a) => {
            let info = get(*a)?;
            (info.0, info.1)
        }
        Expr::Powi { base, exp } => {
            let info = scalar("powi", get(*base)?)?;
            let d = info.1.0;
            let mut scaled = [0i8; 6];
            for (out, &b) in scaled.iter_mut().zip(&d) {
                let product = i32::from(b)
                    .checked_mul(*exp)
                    .ok_or(OptError::DimOverflow {
                        op: "powi",
                        dims: d,
                        exponent: *exp,
                    })?;
                *out = i8::try_from(product).map_err(|_| OptError::DimOverflow {
                    op: "powi",
                    dims: d,
                    exponent: *exp,
                })?;
            }
            (Shape::Scalar, Dims(scaled))
        }
        Expr::Sqrt(a) => {
            let info = scalar("sqrt", get(*a)?)?;
            let d = info.1.0;
            if d.iter().any(|x| x % 2 != 0) {
                return Err(OptError::OddDims { dims: d });
            }
            (Shape::Scalar, Dims(d.map(|x| x / 2)))
        }
        Expr::Exp(a) | Expr::Ln(a) | Expr::Tanh(a) => {
            let op = match e {
                Expr::Exp(_) => "exp",
                Expr::Ln(_) => "ln",
                _ => "tanh",
            };
            let info = scalar(op, get(*a)?)?;
            if info.1 != Dims::NONE {
                return Err(OptError::NonDimensionless { op, dims: info.1.0 });
            }
            (Shape::Scalar, Dims::NONE)
        }
        Expr::Dot(a, b) => {
            let (ia, ib) = (get(*a)?, get(*b)?);
            match (ia.0, ib.0) {
                (Shape::Vector(n), Shape::Vector(m)) if n == m => {
                    (Shape::Scalar, dims_add_checked(ia.1, ib.1, "dot")?)
                }
                (l, r) => {
                    return Err(OptError::ShapeMismatch {
                        op: "dot",
                        left: l,
                        right: r,
                    });
                }
            }
        }
        Expr::NormSq(a) => {
            let info = get(*a)?;
            match info.0 {
                Shape::Vector(_) => (Shape::Scalar, dims_add_checked(info.1, info.1, "norm_sq")?),
                s @ Shape::Scalar => {
                    return Err(OptError::ShapeMismatch {
                        op: "norm_sq",
                        left: s,
                        right: Shape::Vector(1),
                    });
                }
            }
        }
        Expr::Abs(a) => {
            let info = scalar("abs", get(*a)?)?;
            (Shape::Scalar, info.1)
        }
        Expr::PdeResidual {
            study, over, dims, ..
        } => {
            if (over.0 as usize) >= vars.len() {
                return Err(OptError::UnknownVar { id: over.0 });
            }
            validate_string("pde study identifier", study, caps)?;
            (Shape::Scalar, *dims)
        }
        Expr::Expectation { of, uq_config } => {
            validate_string("uq config identifier", uq_config, caps)?;
            let info = scalar("expectation", get(*of)?)?;
            (Shape::Scalar, info.1)
        }
        Expr::Cvar {
            of,
            alpha,
            uq_config,
        } => {
            validate_string("uq config identifier", uq_config, caps)?;
            if !(alpha.is_finite() && *alpha > 0.0 && *alpha < 1.0) {
                return Err(OptError::BadParam {
                    what: "cvar alpha",
                    value: *alpha,
                });
            }
            let info = scalar("cvar", get(*of)?)?;
            (Shape::Scalar, info.1)
        }
        Expr::Quantile { of, q, uq_config } => {
            validate_string("uq config identifier", uq_config, caps)?;
            if !(q.is_finite() && *q > 0.0 && *q < 1.0) {
                return Err(OptError::BadParam {
                    what: "quantile q",
                    value: *q,
                });
            }
            let info = scalar("quantile", get(*of)?)?;
            (Shape::Scalar, info.1)
        }
    };
    let class = children(e)
        .iter()
        .filter_map(|c| lookup(*c).map(|i| i.2))
        .chain([own_class(e)])
        .min()
        .unwrap_or(Class::Smooth);
    Ok((shape, dims, class))
}

/// Re-validate the problem and mint its semantic identity. Aggregate
/// count/alignment/byte caps are preflighted first and refuse before
/// proportional allocation. Inside that bounded envelope, every section
/// runs in deterministic order and returns all discovered violations.
///
/// # Errors
/// A complete [`AdmissionReport`].
#[allow(clippy::too_many_lines)] // one section per admission surface, deterministic order
pub(crate) fn admit_with_caps(
    problem: &Problem,
    caps: &AdmissionCaps,
) -> Result<ProblemAdmission, AdmissionReport> {
    let mut violations = Vec::new();

    // Section 1: aggregate caps and metadata-vector alignment.
    let n_nodes = problem.exprs.len() as u64;
    let aggregate: [(&'static str, u64, u64); 5] = [
        (
            "variables",
            problem.vars.len() as u64,
            u64::from(caps.max_vars),
        ),
        ("expression nodes", n_nodes, u64::from(caps.max_nodes)),
        (
            "objectives",
            problem.objectives.len() as u64,
            u64::from(caps.max_objectives),
        ),
        (
            "constraints",
            problem.constraints.len() as u64,
            u64::from(caps.max_constraints),
        ),
        ("tags", problem.tags.len() as u64, u64::from(caps.max_tags)),
    ];
    for (what, count, cap) in aggregate {
        if count > cap {
            violations.push(AdmissionViolation::Aggregate { what, count, cap });
        }
    }
    for (what, len) in [
        ("shape cache length", problem.shapes.len() as u64),
        ("dimension cache length", problem.dims.len() as u64),
        ("class cache length", problem.classes.len() as u64),
    ] {
        if len != n_nodes {
            violations.push(AdmissionViolation::Aggregate {
                what,
                count: len,
                cap: n_nodes,
            });
        }
    }
    if !violations.is_empty() {
        return Err(AdmissionReport {
            schema_version: ADMISSION_SCHEMA_VERSION,
            violations,
        });
    }

    let total_retained_bytes = problem
        .vars
        .iter()
        .map(|var| name_retained_bytes(&var.name))
        .chain(problem.exprs.iter().map(expr_retained_bytes))
        .chain(
            problem
                .constraints
                .iter()
                .map(|constraint| name_retained_bytes(&constraint.name)),
        )
        .fold(0u64, u64::saturating_add);
    if total_retained_bytes > caps.max_total_retained_bytes {
        return Err(AdmissionReport {
            schema_version: ADMISSION_SCHEMA_VERSION,
            violations: vec![AdmissionViolation::Aggregate {
                what: "total retained/canonical bytes",
                count: total_retained_bytes,
                cap: caps.max_total_retained_bytes,
            }],
        });
    }

    // Section 2: variables (manifold policy, name policy, storage sum).
    let mut total_point_storage: u64 = 0;
    for (i, var) in problem.vars.iter().enumerate() {
        let index = i as u32;
        if let Err(error) = validate_name("variable name", &var.name, caps) {
            violations.push(AdmissionViolation::Var { index, error });
        }
        match validate_manifold(&var.manifold, caps) {
            Ok(()) => {
                let storage = u64::from(var.manifold.point_dim().unwrap_or(u32::MAX));
                total_point_storage = total_point_storage.saturating_add(storage);
            }
            Err(error) => violations.push(AdmissionViolation::Var { index, error }),
        }
    }
    if total_point_storage > caps.max_total_point_storage {
        violations.push(AdmissionViolation::Aggregate {
            what: "total point storage",
            count: total_point_storage,
            cap: caps.max_total_point_storage,
        });
    }
    if usize::try_from(total_point_storage).is_err() {
        violations.push(AdmissionViolation::Aggregate {
            what: "target packed point storage",
            count: total_point_storage,
            cap: usize::MAX as u64,
        });
    }
    if violations.iter().any(|violation| {
        matches!(
            violation,
            AdmissionViolation::Aggregate { what, .. }
                if *what == "total point storage" || *what == "target packed point storage"
        )
    }) {
        return Err(AdmissionReport {
            schema_version: ADMISSION_SCHEMA_VERSION,
            violations,
        });
    }

    // Section 3: nodes — arena-order references (acyclicity proof),
    // full re-derivation, and cache agreement. Re-derivation uses the
    // RECOMPUTED tables so one corrupt cache cannot vouch for another;
    // when a node fails, its cached triple is used to continue the
    // scan so every downstream violation is still reported.
    let mut derived: Vec<NodeInfo> = Vec::with_capacity(problem.exprs.len());
    let mut derived_depths: Vec<u32> = Vec::with_capacity(problem.exprs.len());
    let aligned = problem.shapes.len() == problem.exprs.len()
        && problem.dims.len() == problem.exprs.len()
        && problem.classes.len() == problem.exprs.len();
    for (i, e) in problem.exprs.iter().enumerate() {
        let index = i as u32;
        let mut order_ok = true;
        for child in children(e) {
            if child.0 as usize >= i {
                violations.push(AdmissionViolation::ChildOrder {
                    node: index,
                    child: child.0,
                });
                order_ok = false;
            }
        }
        let cached = if aligned {
            Some((problem.shapes[i], problem.dims[i], problem.classes[i]))
        } else {
            None
        };
        let re_derived = if order_ok {
            let lookup = |n: NodeId| derived.get(n.0 as usize).copied();
            derive_expr(e, &lookup, &problem.vars, caps)
        } else {
            Err(OptError::UnknownNode { id: index })
        };
        let re_derived_depth = if order_ok {
            let lookup = |n: NodeId| derived_depths.get(n.0 as usize).copied();
            derive_depth(e, &lookup, caps)
        } else {
            Err(OptError::UnknownNode { id: index })
        };
        match re_derived {
            Ok(info) => {
                if let Some(c) = cached {
                    if c.0 != info.0 {
                        violations.push(AdmissionViolation::CacheMismatch {
                            node: index,
                            what: "shape",
                        });
                    }
                    if c.1 != info.1 {
                        violations.push(AdmissionViolation::CacheMismatch {
                            node: index,
                            what: "dimensions",
                        });
                    }
                    if c.2 != info.2 {
                        violations.push(AdmissionViolation::CacheMismatch {
                            node: index,
                            what: "class",
                        });
                    }
                }
                derived.push(info);
            }
            Err(error) => {
                if order_ok {
                    violations.push(AdmissionViolation::Node { index, error });
                }
                derived.push(cached.unwrap_or((Shape::Scalar, Dims::NONE, Class::Smooth)));
            }
        }
        match re_derived_depth {
            Ok(depth) => derived_depths.push(depth),
            Err(error) => {
                if order_ok {
                    violations.push(AdmissionViolation::Node { index, error });
                }
                derived_depths.push(1);
            }
        }
    }

    // Section 4: objectives (existence, scalar root, weight policy).
    for (i, o) in problem.objectives.iter().enumerate() {
        let index = i as u32;
        match derived.get(o.node.0 as usize) {
            None => violations.push(AdmissionViolation::Objective {
                index,
                error: OptError::UnknownNode { id: o.node.0 },
            }),
            Some(info) if info.0 != Shape::Scalar => {
                violations.push(AdmissionViolation::Objective {
                    index,
                    error: OptError::NotScalar { node: o.node.0 },
                });
            }
            Some(_) => {}
        }
        if let Err(error) = validate_weight(o.weight) {
            violations.push(AdmissionViolation::Objective { index, error });
        }
    }

    // Section 5: constraints (existence, scalar root, name policy).
    for (i, c) in problem.constraints.iter().enumerate() {
        let index = i as u32;
        match derived.get(c.node.0 as usize) {
            None => violations.push(AdmissionViolation::Constraint {
                index,
                error: OptError::UnknownNode { id: c.node.0 },
            }),
            Some(info) if info.0 != Shape::Scalar => {
                violations.push(AdmissionViolation::Constraint {
                    index,
                    error: OptError::NotScalar { node: c.node.0 },
                });
            }
            Some(_) => {}
        }
        if let Err(error) = validate_name("constraint name", &c.name, caps) {
            violations.push(AdmissionViolation::Constraint { index, error });
        }
    }

    // Section 6: tags (policies + legacy-identity quarantine).
    let mut quarantined = Vec::new();
    for (i, t) in problem.tags.iter().enumerate() {
        if let Err(error) = validate_tag(t, caps) {
            violations.push(AdmissionViolation::Tag {
                index: i as u32,
                error,
            });
        }
        if let ProblemTag::Bilevel {
            inner: BilevelRef::LegacyFnv(h),
        } = t
        {
            quarantined.push(*h);
        }
    }

    if !violations.is_empty() {
        return Err(AdmissionReport {
            schema_version: ADMISSION_SCHEMA_VERSION,
            violations,
        });
    }

    let semantic_id = ProblemSemanticId::mint(&crate::serial::canonical_body_v3(problem));
    Ok(ProblemAdmission {
        schema_version: ADMISSION_SCHEMA_VERSION,
        semantic_id,
        var_count: problem.vars.len() as u32,
        node_count: problem.exprs.len() as u32,
        objective_count: problem.objectives.len() as u32,
        constraint_count: problem.constraints.len() as u32,
        total_point_storage,
        total_retained_bytes,
        quarantined_legacy_identities: quarantined,
    })
}
