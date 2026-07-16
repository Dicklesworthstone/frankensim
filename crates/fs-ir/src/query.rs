//! Declarative query language v0 (plan addendum, Proposal 8): the surface
//! that replaces "run a simulation" with "ask a question, spend a budget".
//!
//! An operator poses a QUESTION with a confidence/tolerance requirement and
//! a budget — *"is max von Mises stress under 180 MPa at 95% confidence,
//! answer for under $50?"* — instead of specifying mesh, solver, and
//! timestep and receiving whatever accuracy falls out. This is the
//! imperative→declarative inversion relational databases made in the 1970s,
//! available here for the same reason: a closed, priced algebra over the
//! operations.
//!
//! This module is v0 and DELIBERATELY SCOPED. It defines:
//! - a fixed MENU of QoI functionals ([`Qoi`]) — no general programs — each
//!   carrying the metadata flags the planner reads ([`QoiMeta`]);
//! - a [`Target`] (an absolute tolerance, a confidence level, or both);
//! - a [`Query`] object = `(QoI, target, budget, deadline)`;
//! - admission ([`Query::admit`]) that type-checks a query against the
//!   design's typed fields ([`FieldRegistry`]) and REJECTS ill-posed
//!   queries with ranked, teaching [`Finding`]s (reusing the admission
//!   machinery), in the same milliseconds-class discipline as study
//!   admission;
//! - a concrete IR surface: [`Query::from_node`] / [`Query::to_node`] make a
//!   query an admissible, versioned IR object (`(query …)` form), not a
//!   stringly-typed request — round-tripping under the AST's `same_shape`
//!   isomorphism.
//!
//! The query LANGUAGE is the durable contract; the planner (a separate
//! bead) is swappable underneath it. Anytime/refusal semantics live with the
//! query RESULT (another bead) — this module owns only posing and admitting.
//!
//! Determinism: [`Query::admit`] runs its checks in a fixed order and emits
//! findings deterministically, so a replayed query reproduces the same
//! verdict (the addendum's determinism-as-contract requirement).

use std::collections::BTreeMap;

use crate::admission::{Finding, RankedFix, Severity};
use crate::ast::{Node, NodeKind, Span};
use crate::{IrError, IrErrorKind};
use fs_qty::Dims;

/// Volume dimensions (`m³`) — the measure a spatial integral multiplies by.
const VOLUME_DIMS: Dims = Dims([3, 0, 0, 0, 0, 0]);
/// Time dimensions (`s`) — the deadline's expected dimension.
const TIME_DIMS: Dims = Dims([0, 0, 1, 0, 0, 0]);

/// The QoI functional menu (v0). A fixed set of forms covering the wedge
/// vertical's real questions — NOT a general programming surface.
#[derive(Debug, Clone, PartialEq)]
pub enum Qoi {
    /// `max` of a named scalar field over a named region.
    MaxOverRegion {
        /// The field being interrogated.
        field: String,
        /// The region the max is taken over.
        region: String,
    },
    /// Spatial integral `∫ f dV` of a named field over a named region.
    Integral {
        /// The integrand field.
        field: String,
        /// The region of integration.
        region: String,
    },
    /// Exceedance probability `P(max over region f ≥ threshold)` under a
    /// declared environment distribution (Proposal F). The result is a
    /// dimensionless probability; the threshold carries the field's dims.
    Exceedance {
        /// The field whose exceedance is measured.
        field: String,
        /// The region over which the field is reduced.
        region: String,
        /// Threshold value in SI base units.
        threshold: f64,
        /// The threshold's dimensions (must match the field's).
        threshold_dims: Dims,
        /// The declared environment/hazard distribution the probability is
        /// taken under (a Proposal F artifact, referenced by name).
        environment: String,
    },
}

/// Planner-facing metadata every QoI advertises. (Whether the QoI is
/// inherently probabilistic is determined by the variant, not stored here —
/// see [`Qoi::is_probabilistic`].)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QoiMeta {
    /// Is the functional linear in the field? (linear QoIs admit cheap
    /// adjoint-weighted error estimates).
    pub linear: bool,
    /// Is an adjoint available for goal-oriented (DWR) refinement?
    pub adjoint_available: bool,
    /// Does the fidelity ladder apply (can this QoI be evaluated on a
    /// coarser rung and prolongated)?
    pub ladder_applicable: bool,
}

impl Qoi {
    /// Stable form name (for logging/verdicts).
    #[must_use]
    pub fn kind_name(&self) -> &'static str {
        match self {
            Qoi::MaxOverRegion { .. } => "max-over-region",
            Qoi::Integral { .. } => "integral",
            Qoi::Exceedance { .. } => "exceedance",
        }
    }

    /// The interrogated field's name.
    #[must_use]
    pub fn field(&self) -> &str {
        match self {
            Qoi::MaxOverRegion { field, .. }
            | Qoi::Integral { field, .. }
            | Qoi::Exceedance { field, .. } => field,
        }
    }

    /// The region's name.
    #[must_use]
    pub fn region(&self) -> &str {
        match self {
            Qoi::MaxOverRegion { region, .. }
            | Qoi::Integral { region, .. }
            | Qoi::Exceedance { region, .. } => region,
        }
    }

    /// Planner metadata for this functional.
    #[must_use]
    pub fn meta(&self) -> QoiMeta {
        match self {
            // max is nonlinear (a pointwise sup) but adjoint-estimable via a
            // smoothed-max surrogate; ladder applies.
            Qoi::MaxOverRegion { .. } => QoiMeta {
                linear: false,
                adjoint_available: true,
                ladder_applicable: true,
            },
            // a spatial integral is linear in the field — the DWR sweet spot.
            Qoi::Integral { .. } => QoiMeta {
                linear: true,
                adjoint_available: true,
                ladder_applicable: true,
            },
            // exceedance is a probability under an environment ensemble.
            Qoi::Exceedance { .. } => QoiMeta {
                linear: false,
                adjoint_available: false,
                ladder_applicable: true,
            },
        }
    }

    /// Is the QoI inherently probabilistic (needs an environment
    /// distribution), rather than a deterministic field functional?
    #[must_use]
    pub fn is_probabilistic(&self) -> bool {
        matches!(self, Qoi::Exceedance { .. })
    }

    /// The dimensions of the QoI's VALUE, given the interrogated field's
    /// dimensions: `max` inherits the field's dims; `integral` multiplies by
    /// volume; `exceedance` is a dimensionless probability. `None` when the
    /// integral's volume factor overflows the i8 exponent domain (bead
    /// sj31i.11: a clamped exponent could alias false physics, so the
    /// admission below rejects instead).
    #[must_use]
    pub fn value_dims(&self, field_dims: Dims) -> Option<Dims> {
        match self {
            Qoi::MaxOverRegion { .. } => Some(field_dims),
            Qoi::Integral { .. } => field_dims.checked_plus(VOLUME_DIMS),
            Qoi::Exceedance { .. } => Some(Dims::NONE),
        }
    }
}

/// What the operator requires of the answer.
#[derive(Debug, Clone, PartialEq)]
pub enum Target {
    /// An absolute tolerance on the QoI value (SI value + dims). The answer
    /// is good enough when the certified interval is narrower than this.
    Tolerance {
        /// Half-width, in SI base units.
        value: f64,
        /// The tolerance's dimensions (must match the QoI value dims).
        dims: Dims,
    },
    /// A statistical confidence level in `(0, 1)` that the QoI meets the
    /// implied predicate. `1.0` is deliberately uncertifiable.
    Confidence(f64),
    /// Both an absolute tolerance AND a confidence level — the robust ask.
    ToleranceAndConfidence {
        /// Tolerance half-width in SI base units.
        value: f64,
        /// The tolerance's dimensions.
        dims: Dims,
        /// The confidence level in `(0, 1)`.
        confidence: f64,
    },
}

impl Target {
    /// The tolerance `(value, dims)` if this target carries one.
    #[must_use]
    pub fn tolerance(&self) -> Option<(f64, Dims)> {
        match self {
            Target::Tolerance { value, dims }
            | Target::ToleranceAndConfidence { value, dims, .. } => Some((*value, *dims)),
            Target::Confidence(_) => None,
        }
    }

    /// The confidence level if this target carries one.
    #[must_use]
    pub fn confidence(&self) -> Option<f64> {
        match self {
            Target::Confidence(c) | Target::ToleranceAndConfidence { confidence: c, .. } => {
                Some(*c)
            }
            Target::Tolerance { .. } => None,
        }
    }
}

/// A declarative query v0: a question plus a budget.
#[derive(Debug, Clone, PartialEq)]
pub struct Query {
    /// The quantity of interest.
    pub qoi: Qoi,
    /// What the answer must satisfy.
    pub target: Target,
    /// The compute budget in dollars (the priced "$50").
    pub budget_usd: f64,
    /// The wall-clock deadline in seconds.
    pub deadline_s: f64,
    /// Source provenance (default for programmatically-built queries).
    pub span: Span,
}

/// The design's typed fields: name → SI dimensions. In production this is
/// supplied by the design's function-space interface types (Proposal 13);
/// admission consults it to reject a QoI over a field that does not exist
/// and to check tolerance dimensions.
#[derive(Debug, Clone, Default)]
pub struct FieldRegistry {
    fields: std::collections::BTreeMap<String, Dims>,
}

impl FieldRegistry {
    /// An empty registry.
    #[must_use]
    pub fn new() -> FieldRegistry {
        FieldRegistry {
            fields: std::collections::BTreeMap::new(),
        }
    }

    /// Register a field's dimensions (builder style).
    #[must_use]
    pub fn with(mut self, name: &str, dims: Dims) -> FieldRegistry {
        self.fields.insert(name.to_string(), dims);
        self
    }

    /// The dimensions of a named field, if it exists.
    #[must_use]
    pub fn dims_of(&self, name: &str) -> Option<Dims> {
        self.fields.get(name).copied()
    }

    /// Field names in deterministic (sorted) order — used to teach the
    /// operator what IS available when they name a missing field.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.fields.keys().map(String::as_str).collect()
    }
}

/// The verdict of admitting a query.
#[derive(Debug, Clone, PartialEq)]
pub struct QueryAdmission {
    /// The QoI form name (logging).
    pub qoi_kind: &'static str,
    /// The QoI's value dimensions under the registry (`None` if the field is
    /// unknown, so the dims cannot be derived).
    pub value_dims: Option<Dims>,
    /// True iff no finding has `Reject` severity.
    pub admitted: bool,
    /// Every finding, deterministically ordered (best-first within a check).
    pub findings: Vec<Finding>,
}

impl QueryAdmission {
    /// A one-line, machine-parseable diagnosis for structured logging (never
    /// printed to stdout by library code — the caller decides).
    #[must_use]
    pub fn diagnosis(&self) -> String {
        if self.admitted {
            format!("query admitted: {} (no rejects)", self.qoi_kind)
        } else {
            let rejects = self
                .findings
                .iter()
                .filter(|f| f.severity == Severity::Reject)
                .count();
            format!(
                "query rejected: {} ({rejects} blocking finding(s))",
                self.qoi_kind
            )
        }
    }
}

impl Query {
    /// A query with a default (unset) span.
    #[must_use]
    pub fn new(qoi: Qoi, target: Target, budget_usd: f64, deadline_s: f64) -> Query {
        Query {
            qoi,
            target,
            budget_usd,
            deadline_s,
            span: Span::default(),
        }
    }

    /// Admit the query against a field registry: type-check it and return a
    /// verdict with ranked teaching fixes for anything ill-posed. Pure and
    /// deterministic; runs in constant time (no solves). The checks, in fixed
    /// order: field existence, budget, deadline, confidence, tolerance value,
    /// and dimensional consistency (see the per-check helpers).
    #[must_use]
    pub fn admit(&self, fields: &FieldRegistry) -> QueryAdmission {
        let field_dims = fields.dims_of(self.qoi.field());
        let mut findings = Vec::new();
        findings.extend(self.check_field(fields, field_dims));
        findings.extend(self.check_budget());
        findings.extend(self.check_deadline());
        findings.extend(self.check_confidence());
        findings.extend(self.check_tolerance_value());
        findings.extend(self.check_dimensions(field_dims));
        let admitted = !findings.iter().any(|f| f.severity == Severity::Reject);
        QueryAdmission {
            qoi_kind: self.qoi.kind_name(),
            value_dims: field_dims.and_then(|fd| self.qoi.value_dims(fd)),
            admitted,
            findings,
        }
    }

    /// Build a single-fix rejection finding pointing at this query's span.
    fn reject(
        &self,
        check: &'static str,
        what: String,
        action: String,
        qoi_impact: &str,
    ) -> Finding {
        Finding {
            check,
            severity: Severity::Reject,
            span: self.span,
            what,
            fixes: vec![RankedFix {
                action,
                predicted_wall_s: None,
                qoi_impact: qoi_impact.to_string(),
            }],
        }
    }

    /// The QoI's field must exist in the design.
    fn check_field(&self, fields: &FieldRegistry, field_dims: Option<Dims>) -> Option<Finding> {
        if field_dims.is_some() {
            return None;
        }
        let available = fields.names().join(", ");
        let action = if available.is_empty() {
            "declare the field on the design before querying it".to_string()
        } else {
            format!("use one of the design's fields: {available}")
        };
        Some(self.reject(
            "query.field",
            format!(
                "QoI names field '{}', which is not a field of this design",
                self.qoi.field()
            ),
            action,
            "query cannot be planned against a nonexistent field",
        ))
    }

    /// The dollar budget must be finite and positive.
    fn check_budget(&self) -> Option<Finding> {
        if self.budget_usd.is_finite() && self.budget_usd > 0.0 {
            return None;
        }
        Some(self.reject(
            "query.budget",
            format!(
                "budget must be a finite positive dollar amount, got {}",
                self.budget_usd
            ),
            "grant a positive compute budget, e.g. (budget 50)".to_string(),
            "a zero budget can discharge no query",
        ))
    }

    /// The wall deadline must be finite and positive.
    fn check_deadline(&self) -> Option<Finding> {
        if self.deadline_s.is_finite() && self.deadline_s > 0.0 {
            return None;
        }
        Some(self.reject(
            "query.deadline",
            format!(
                "deadline must be a finite positive number of seconds, got {}",
                self.deadline_s
            ),
            "give a future deadline, e.g. (deadline 30s)".to_string(),
            "a past/zero deadline leaves no time to answer",
        ))
    }

    /// Any confidence must lie strictly in `(0, 1)` — 100% is uncertifiable.
    fn check_confidence(&self) -> Option<Finding> {
        let c = self.target.confidence()?;
        if c.is_finite() && c > 0.0 && c < 1.0 {
            return None;
        }
        let (what, action) = if c >= 1.0 {
            (
                format!("confidence {c} is uncertifiable — no finite evidence proves 100%"),
                "request a confidence strictly below 1.0, e.g. (confidence 0.95)".to_string(),
            )
        } else {
            (
                format!("confidence {c} must be strictly greater than 0"),
                "request a confidence in (0, 1), e.g. (confidence 0.95)".to_string(),
            )
        };
        Some(self.reject(
            "query.confidence",
            what,
            action,
            "an uncertifiable target can never be met",
        ))
    }

    /// Any tolerance half-width must be finite and positive.
    fn check_tolerance_value(&self) -> Option<Finding> {
        let (value, _) = self.target.tolerance()?;
        if value.is_finite() && value > 0.0 {
            return None;
        }
        Some(self.reject(
            "query.target",
            format!("tolerance half-width must be finite and positive, got {value}"),
            "request a positive tolerance, e.g. (tolerance 5MPa)".to_string(),
            "a zero/negative tolerance demands exactness no solve can certify",
        ))
    }

    /// A tolerance's dims must match the QoI value dims, and an exceedance
    /// threshold's dims must match the field's (a stress tolerance on a
    /// probability, or a "5 second" tolerance on a stress, is
    /// self-contradictory). Only meaningful once the field exists.
    fn check_dimensions(&self, field_dims: Option<Dims>) -> Vec<Finding> {
        let Some(fd) = field_dims else {
            return Vec::new();
        };
        let mut out = Vec::new();
        let Some(value_dims) = self.qoi.value_dims(fd) else {
            return vec![
                self.reject(
                    "query.dimensions",
                    format!(
                        "the {} of field '{}' (dims {:?}) overflows the supported i8 exponent \
                     domain",
                        self.qoi.kind_name(),
                        self.qoi.field(),
                        fd.0
                    ),
                    "interrogate a field whose dimensions stay within the i8 exponent domain \
                 after the QoI's volume factor"
                        .to_string(),
                    "dimension overflow; the QoI value cannot carry admissible physics",
                ),
            ];
        };
        if let Some((_, tol_dims)) = self.target.tolerance()
            && tol_dims != value_dims
        {
            out.push(self.reject(
                "query.dimensions",
                format!(
                    "tolerance dims {:?} do not match the QoI value dims {:?} (a {} of field '{}')",
                    tol_dims,
                    value_dims,
                    self.qoi.kind_name(),
                    self.qoi.field()
                ),
                format!("state the tolerance in units matching the QoI value dims {value_dims:?}"),
                "a dimensionally-inconsistent tolerance is meaningless",
            ));
        }
        if let Qoi::Exceedance { threshold_dims, .. } = &self.qoi
            && *threshold_dims != fd
        {
            out.push(self.reject(
                "query.dimensions",
                format!(
                    "exceedance threshold dims {:?} do not match field '{}' dims {:?}",
                    threshold_dims,
                    self.qoi.field(),
                    fd
                ),
                format!("state the threshold in the field's units (dims {fd:?})"),
                "an off-dimension threshold compares incomparable quantities",
            ));
        }
        out
    }

    /// Recognize a query from its concrete IR form:
    /// `(query <qoi> <target> (budget N) (deadline T))`, where `<qoi>` is one
    /// of `(max :field "…" :region "…")`, `(integral …)`, or
    /// `(exceedance :field "…" :region "…" :threshold Q :env "…")`, and
    /// `<target>` is `(tolerance Q)`, `(confidence F)`, or
    /// `(tolerance Q :confidence F)`.
    ///
    /// # Errors
    /// A structured [`IrError`] pointing at the malformed clause.
    pub fn from_node(node: &Node) -> Result<Query, IrError> {
        node.validate()?;
        let items = match node.head() {
            Some("query") => node.items().expect("head implies list"),
            _ => {
                return Err(malformed(
                    node.span,
                    "expected a (query …) form",
                    "wrap it as (query <qoi> <target> (budget N) (deadline T))",
                ));
            }
        };
        let qoi = parse_qoi(items.get(1).ok_or_else(|| {
            malformed(
                node.span,
                "query has no QoI",
                "add a QoI form, e.g. (max :field \"vm\" :region \"bracket\")",
            )
        })?)?;
        let target = parse_target(items.get(2).ok_or_else(|| {
            malformed(
                node.span,
                "query has no target",
                "add a target, e.g. (confidence 0.95)",
            )
        })?)?;
        let mut budget_usd = None;
        let mut deadline_s = None;
        for clause in &items[3..] {
            match clause.head() {
                Some("budget") => {
                    if budget_usd.is_some() {
                        return Err(malformed(
                            clause.span,
                            "duplicate query.budget clause",
                            "provide exactly one (budget N) clause",
                        ));
                    }
                    budget_usd = Some(clause_number(clause)?);
                }
                Some("deadline") => {
                    if deadline_s.is_some() {
                        return Err(malformed(
                            clause.span,
                            "duplicate query.deadline clause",
                            "provide exactly one (deadline T) clause",
                        ));
                    }
                    deadline_s = Some(clause_seconds(clause)?);
                }
                _ => {
                    return Err(malformed(
                        clause.span,
                        "unknown query clause",
                        "queries take (budget N) and (deadline T) after the QoI and target",
                    ));
                }
            }
        }
        let budget_usd = budget_usd.ok_or_else(|| {
            malformed(
                node.span,
                "missing query.budget clause",
                "provide exactly one (budget N) clause",
            )
        })?;
        let deadline_s = deadline_s.ok_or_else(|| {
            malformed(
                node.span,
                "missing query.deadline clause",
                "provide exactly one (deadline T) clause",
            )
        })?;
        if items.len() != 5
            || items[3].head() != Some("budget")
            || items[4].head() != Some("deadline")
        {
            return Err(malformed(
                node.span,
                "query clauses are not in the exact canonical schema",
                "use (query <qoi> <target> (budget N) (deadline T))",
            ));
        }
        Ok(Query {
            qoi,
            target,
            budget_usd,
            deadline_s,
            span: node.span,
        })
    }

    /// Emit the query as its concrete IR form (synthetic spans). Round-trips
    /// with [`Query::from_node`] under `same_shape` (semantic equality
    /// ignores spans and quantity presentation text).
    pub fn to_node(&self) -> Result<Node, IrError> {
        let node = list(vec![
            sym("query"),
            self.qoi_node()?,
            self.target_node()?,
            list(vec![
                sym("budget"),
                Node::synthetic(NodeKind::Float(self.budget_usd)),
            ]),
            list(vec![
                sym("deadline"),
                Node::quantity(self.deadline_s, TIME_DIMS)?,
            ]),
        ]);
        node.validate()?;
        Ok(node)
    }

    fn qoi_node(&self) -> Result<Node, IrError> {
        Ok(match &self.qoi {
            Qoi::MaxOverRegion { field, region } => list(vec![
                sym("max"),
                kw("field"),
                str_node(field),
                kw("region"),
                str_node(region),
            ]),
            Qoi::Integral { field, region } => list(vec![
                sym("integral"),
                kw("field"),
                str_node(field),
                kw("region"),
                str_node(region),
            ]),
            Qoi::Exceedance {
                field,
                region,
                threshold,
                threshold_dims,
                environment,
            } => list(vec![
                sym("exceedance"),
                kw("field"),
                str_node(field),
                kw("region"),
                str_node(region),
                kw("threshold"),
                Node::quantity(*threshold, *threshold_dims)?,
                kw("env"),
                str_node(environment),
            ]),
        })
    }

    fn target_node(&self) -> Result<Node, IrError> {
        Ok(match &self.target {
            Target::Tolerance { value, dims } => {
                list(vec![sym("tolerance"), Node::quantity(*value, *dims)?])
            }
            Target::Confidence(c) => list(vec![
                sym("confidence"),
                Node::synthetic(NodeKind::Float(*c)),
            ]),
            Target::ToleranceAndConfidence {
                value,
                dims,
                confidence,
            } => list(vec![
                sym("tolerance"),
                Node::quantity(*value, *dims)?,
                kw("confidence"),
                Node::synthetic(NodeKind::Float(*confidence)),
            ]),
        })
    }
}

// ---- s-expr helpers -------------------------------------------------------

fn malformed(span: Span, detail: &str, hint: &str) -> IrError {
    IrError {
        span,
        kind: IrErrorKind::MalformedClause,
        detail: detail.to_string(),
        hint: hint.to_string(),
    }
}

fn sym(s: &str) -> Node {
    Node::synthetic(NodeKind::Symbol(s.to_string()))
}
fn kw(s: &str) -> Node {
    Node::synthetic(NodeKind::Keyword(s.to_string()))
}
fn str_node(s: &str) -> Node {
    Node::synthetic(NodeKind::Str(s.to_string()))
}
fn list(items: Vec<Node>) -> Node {
    Node::synthetic(NodeKind::List(items))
}
/// Parse an exact `:key value` schema without silently skipping positional,
/// unknown, duplicate, or dangling arguments.
fn keyword_args<'a>(
    items: &'a [Node],
    start: usize,
    form: &str,
    allowed: &[&str],
) -> Result<BTreeMap<&'a str, &'a Node>, IrError> {
    let mut out = BTreeMap::new();
    let mut index = start;
    while index < items.len() {
        let NodeKind::Keyword(key) = &items[index].kind else {
            return Err(malformed(
                items[index].span,
                &format!("unexpected positional argument at {form}[{index}]"),
                "use only the documented :key value pairs",
            ));
        };
        if !allowed.contains(&key.as_str()) {
            return Err(malformed(
                items[index].span,
                &format!("unknown {form}.{key} keyword"),
                &format!("known {form} keywords: {}", allowed.join(", ")),
            ));
        }
        let value = items.get(index + 1).ok_or_else(|| {
            malformed(
                items[index].span,
                &format!("dangling {form}.{key} keyword"),
                "add its value or remove the keyword",
            )
        })?;
        if out.insert(key.as_str(), value).is_some() {
            return Err(malformed(
                items[index].span,
                &format!("duplicate {form}.{key} keyword"),
                "provide each keyword exactly once",
            ));
        }
        index += 2;
    }
    Ok(out)
}

fn kw_str(
    args: &BTreeMap<&str, &Node>,
    key: &str,
    form: &str,
    span: Span,
) -> Result<String, IrError> {
    if let Some(value) = args.get(key) {
        if let NodeKind::Str(text) = &value.kind {
            return Ok(text.clone());
        }
        return Err(malformed(
            value.span,
            &format!("expected a string at {form}.{key}"),
            "e.g. :field \"vm\"",
        ));
    }
    Err(malformed(
        span,
        &format!("missing required {form}.{key} keyword"),
        "provide every required keyword exactly once",
    ))
}

fn kw_qty(
    args: &BTreeMap<&str, &Node>,
    key: &str,
    form: &str,
    span: Span,
) -> Result<(f64, Dims), IrError> {
    if let Some(value) = args.get(key) {
        if let NodeKind::Qty { value, dims, .. } = &value.kind {
            return Ok((*value, *dims));
        }
        return Err(malformed(
            value.span,
            &format!("expected a dimensioned quantity at {form}.{key}"),
            "e.g. :threshold 180MPa",
        ));
    }
    Err(malformed(
        span,
        &format!("missing required {form}.{key} keyword"),
        "provide every required keyword exactly once",
    ))
}

fn parse_qoi(node: &Node) -> Result<Qoi, IrError> {
    let items = node.items().ok_or_else(|| {
        malformed(
            node.span,
            "QoI must be a form",
            "e.g. (max :field \"vm\" :region \"bracket\")",
        )
    })?;
    match node.head() {
        Some("max") => {
            let args = keyword_args(items, 1, "query.qoi.max", &["field", "region"])?;
            Ok(Qoi::MaxOverRegion {
                field: kw_str(&args, "field", "query.qoi.max", node.span)?,
                region: kw_str(&args, "region", "query.qoi.max", node.span)?,
            })
        }
        Some("integral") => {
            let args = keyword_args(items, 1, "query.qoi.integral", &["field", "region"])?;
            Ok(Qoi::Integral {
                field: kw_str(&args, "field", "query.qoi.integral", node.span)?,
                region: kw_str(&args, "region", "query.qoi.integral", node.span)?,
            })
        }
        Some("exceedance") => {
            let args = keyword_args(
                items,
                1,
                "query.qoi.exceedance",
                &["field", "region", "threshold", "env"],
            )?;
            let (threshold, threshold_dims) =
                kw_qty(&args, "threshold", "query.qoi.exceedance", node.span)?;
            Ok(Qoi::Exceedance {
                field: kw_str(&args, "field", "query.qoi.exceedance", node.span)?,
                region: kw_str(&args, "region", "query.qoi.exceedance", node.span)?,
                threshold,
                threshold_dims,
                environment: kw_str(&args, "env", "query.qoi.exceedance", node.span)?,
            })
        }
        _ => Err(malformed(
            node.span,
            "unknown QoI form",
            "v0 supports (max …), (integral …), (exceedance …)",
        )),
    }
}

fn parse_target(node: &Node) -> Result<Target, IrError> {
    let items = node
        .items()
        .ok_or_else(|| malformed(node.span, "target must be a form", "e.g. (confidence 0.95)"))?;
    match node.head() {
        Some("tolerance") => {
            let (value, dims) = match items.get(1).map(|n| &n.kind) {
                Some(NodeKind::Qty { value, dims, .. }) => (*value, *dims),
                _ => {
                    return Err(malformed(
                        node.span,
                        "tolerance needs a dimensioned quantity",
                        "e.g. (tolerance 5MPa)",
                    ));
                }
            };
            let args = keyword_args(items, 2, "query.target.tolerance", &["confidence"])?;
            if let Some(confidence_node) = args.get("confidence") {
                let confidence = float_of(confidence_node)?;
                Ok(Target::ToleranceAndConfidence {
                    value,
                    dims,
                    confidence,
                })
            } else {
                Ok(Target::Tolerance { value, dims })
            }
        }
        Some("confidence") => {
            if items.len() != 2 {
                let offense = items.get(2).map_or(node.span, |item| item.span);
                return Err(malformed(
                    offense,
                    "query.target.confidence takes exactly one number",
                    "use (confidence F)",
                ));
            }
            let n = items.get(1).ok_or_else(|| {
                malformed(
                    node.span,
                    "confidence needs a number",
                    "e.g. (confidence 0.95)",
                )
            })?;
            Ok(Target::Confidence(float_of(n)?))
        }
        _ => Err(malformed(
            node.span,
            "unknown target form",
            "v0 supports (tolerance Q), (confidence F), (tolerance Q :confidence F)",
        )),
    }
}

fn float_of(node: &Node) -> Result<f64, IrError> {
    match &node.kind {
        NodeKind::Float(f) => Ok(*f),
        NodeKind::Int(i) if integer_is_exact_f64(*i) => {
            #[allow(clippy::cast_precision_loss)]
            let value = *i as f64;
            Ok(value)
        }
        NodeKind::Int(_) => Err(malformed(
            node.span,
            "integer cannot be represented exactly as an f64 query number",
            "use an exactly representable integer or an explicit finite float",
        )),
        _ => Err(malformed(node.span, "expected a number", "e.g. 0.95")),
    }
}

fn integer_is_exact_f64(value: i64) -> bool {
    let magnitude = value.unsigned_abs();
    if magnitude == 0 {
        return true;
    }
    let significant = magnitude >> magnitude.trailing_zeros();
    (u64::BITS - significant.leading_zeros()) <= 53
}

fn clause_number(clause: &Node) -> Result<f64, IrError> {
    let items = clause
        .items()
        .ok_or_else(|| malformed(clause.span, "malformed clause", "e.g. (budget 50)"))?;
    if items.len() != 2 {
        let offense = items.get(2).map_or(clause.span, |item| item.span);
        return Err(malformed(
            offense,
            "query.budget takes exactly one number",
            "use (budget N)",
        ));
    }
    let n = items
        .get(1)
        .ok_or_else(|| malformed(clause.span, "clause needs a value", "e.g. (budget 50)"))?;
    float_of(n)
}

fn clause_seconds(clause: &Node) -> Result<f64, IrError> {
    let items = clause
        .items()
        .ok_or_else(|| malformed(clause.span, "malformed clause", "e.g. (deadline 30s)"))?;
    if items.len() != 2 {
        let offense = items.get(2).map_or(clause.span, |item| item.span);
        return Err(malformed(
            offense,
            "query.deadline takes exactly one time",
            "use (deadline 30s)",
        ));
    }
    match items.get(1).map(|n| &n.kind) {
        // accept a time quantity (dims = seconds) or a bare number of seconds.
        Some(NodeKind::Qty { value, dims, .. }) if *dims == TIME_DIMS => Ok(*value),
        Some(NodeKind::Qty { .. }) => Err(malformed(
            items[1].span,
            "deadline must have time dimensions",
            "e.g. (deadline 30s) or (deadline 5min)",
        )),
        Some(NodeKind::Float(_) | NodeKind::Int(_)) => float_of(&items[1]),
        _ => Err(malformed(
            clause.span,
            "deadline needs a time",
            "e.g. (deadline 30s)",
        )),
    }
}
