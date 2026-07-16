//! STATIC ADMISSION (plan §11.1, the gp3.5 bead): before anything
//! executes, the study passes admission — an ill-typed study is rejected
//! in MILLISECONDS with a structured, spans-attached diagnosis and RANKED
//! FIXES, not discovered at hour six.
//!
//! Raw syntax is first bound to an IR version and fully lowered; malformed
//! shorthand refuses before any authority decision, and the report binds the
//! exact raw and lowered canonical identities. Six dimensions are then timed
//! over only the explicit semantics: the Five Explicits (structure),
//! dimensional analysis (fs-qty dims through the expression graph), chart
//! routability (the Rep Router as an admission predicate), budget feasibility
//! (learned fs-plan cost models with cost-derived fix estimates), capability
//! sufficiency (session token globs), and regime gating (fs-regime reports;
//! `(assert (regime.allows …))` is enforced, and `flux.*` verbs are checked
//! against the report's model verdicts).

use crate::ast::{CountUnit, Node, NodeKind, Span};
use crate::lower::lower;
use crate::study::Study;
use crate::{IR_VERSION, IrError, VersionedProgram};
use fs_geom::{CostOracle, RoutePlanError, RouteRequest, Router};
use fs_plan::{CostEvidenceClass, FreshnessPolicy, SealedCostModel, StalenessVerdict};
use fs_qty::Dims;
use fs_regime::RegimeReport;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::time::Instant;

/// Finding severity: `Reject` blocks admission, `Warn` admits with notice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// Advisory; the study still admits.
    Warn,
    /// Blocks admission.
    Reject,
}

/// One ranked, actionable fix.
#[derive(Debug, Clone, PartialEq)]
pub struct RankedFix {
    /// What to change, concretely.
    pub action: String,
    /// Predicted total wall after the fix (cost-model-derived), if known.
    pub predicted_wall_s: Option<f64>,
    /// The QoI/accuracy consequence, stated honestly.
    pub qoi_impact: String,
}

/// One admission finding.
#[derive(Debug, Clone, PartialEq)]
pub struct Finding {
    /// Which admission dimension produced it.
    pub check: &'static str,
    /// Reject or Warn.
    pub severity: Severity,
    /// The source span the finding points at.
    pub span: Span,
    /// The diagnosis.
    pub what: String,
    /// Ranked fixes (best first).
    pub fixes: Vec<RankedFix>,
}

/// The session's capability token (fs-session owns issuance; admission
/// only checks sufficiency).
#[derive(Debug, Clone, PartialEq)]
pub struct SessionCapability {
    /// Granted operator globs (`flux.*`, `ascent.optimize`, …).
    pub ops: Vec<String>,
    /// Core grant.
    pub cores: u64,
    /// Memory grant in bytes.
    pub mem_bytes: u64,
    /// Wall grant in seconds.
    pub wall_s: f64,
}

/// What to do when a solver choice violates the regime report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegimePolicy {
    /// Violations reject.
    Reject,
    /// Violations warn (exploratory sessions).
    Warn,
}

/// A chart-conversion requirement extracted by lowering (or declared by
/// the caller): admission asks the Router whether a route exists.
#[derive(Debug, Clone, PartialEq)]
pub struct ChartRequirement {
    /// Source chart kind.
    pub from: String,
    /// Destination chart kind.
    pub to: String,
    /// Reference magnitude for error grounding.
    pub scale: f64,
    /// Absolute error budget.
    pub max_abs_error: f64,
    /// Cost budget, seconds.
    pub max_cost_s: f64,
}

/// Everything admission checks against.
pub struct AdmissionContext<'a> {
    /// The Rep Router and its cost oracle (chart feasibility).
    pub router: Option<(&'a Router, &'a dyn CostOracle)>,
    /// Freshness contract for assessing sealed cost evidence at the
    /// decision point (bead l2k92). `None` keeps mint-time classes
    /// (no clock is ever read here); `Some` degrades receipt-backed
    /// models that are aged out, machine-drifted, or future-stamped to
    /// `StaleRooflineReceipt` before the admission notice is written.
    pub cost_freshness: Option<CostFreshnessContext>,
    /// Conversion requirements to verify.
    pub chart_requirements: Vec<ChartRequirement>,
    /// Learned wall-cost models keyed by verb head. Sealed (bead 2pmb):
    /// receipt-backed authority is mintable only by fs-plan's exact
    /// roofline loader; caller-fitted models enter as explicitly
    /// provisional evidence and every admission that costs through one
    /// carries a Warn finding naming the class.
    pub cost_models: BTreeMap<String, SealedCostModel>,
    /// The session capability token.
    pub capability: Option<SessionCapability>,
    /// The regime report for the study's physics, when computed.
    pub regime: Option<&'a RegimeReport>,
    /// Violation policy for regime gating.
    pub regime_policy: RegimePolicy,
}

/// Caller-supplied freshness contract for cost-evidence assessment
/// (bead l2k92): the observation time, the current machine
/// fingerprint, and the preregistered horizon. Pure data — admission
/// never reads a clock or host identity itself.
#[derive(Debug, Clone, PartialEq)]
pub struct CostFreshnessContext {
    /// Observation time in ledger nanoseconds.
    pub now_ns: i64,
    /// The current machine fingerprint (exact roofline machine key).
    pub machine: Vec<u8>,
    /// The preregistered freshness horizon.
    pub policy: FreshnessPolicy,
}

/// Per-check timing (the milliseconds-class evidence).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckTiming {
    /// Check name.
    pub check: &'static str,
    /// Elapsed microseconds.
    pub micros: u128,
}

/// Exact semantic identity binding used by one admission decision.
///
/// The canonical strings are complete versioned envelopes, not weak hashes:
/// equality therefore means byte-identical canonical versioned programs.
/// A malformed raw AST has no canonical raw identity. A refusal after raw
/// binding retains that identity but has no lowered identity; no invalid or
/// partially lowered tree can reach an authority check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoweringReceipt {
    ir_version: u32,
    raw_canonical: Option<String>,
    lowered_canonical: Option<String>,
}

impl LoweringReceipt {
    fn unbound_current() -> Self {
        Self {
            ir_version: IR_VERSION,
            raw_canonical: None,
            lowered_canonical: None,
        }
    }

    fn for_program(program: &VersionedProgram) -> Self {
        Self {
            ir_version: program.version(),
            raw_canonical: Some(program.print_sexpr()),
            lowered_canonical: None,
        }
    }

    fn bind_lowered(&mut self, node: &Node) -> Result<(), IrError> {
        let lowered = VersionedProgram::try_current(node.clone())?;
        self.lowered_canonical = Some(lowered.print_sexpr());
        Ok(())
    }

    /// IR language version governing both identities.
    #[must_use]
    pub const fn ir_version(&self) -> u32 {
        self.ir_version
    }

    /// Exact canonical versioned identity submitted for admission.
    ///
    /// Returns the empty string only when a malformed caller-forged AST could
    /// not be bound to any canonical versioned identity. New code that needs
    /// to distinguish that refusal should use [`Self::raw_canonical_opt`].
    #[must_use]
    pub fn raw_canonical(&self) -> &str {
        self.raw_canonical.as_deref().unwrap_or("")
    }

    /// Exact submitted identity, or `None` when raw binding itself refused.
    #[must_use]
    pub fn raw_canonical_opt(&self) -> Option<&str> {
        self.raw_canonical.as_deref()
    }

    /// Exact canonical versioned identity actually inspected by authority
    /// checks, or `None` when lowering refused.
    #[must_use]
    pub fn lowered_canonical(&self) -> Option<&str> {
        self.lowered_canonical.as_deref()
    }
}

/// The admission verdict.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmissionReport {
    /// Study name, or a phase-specific placeholder when admission refuses
    /// before study recognition.
    pub study: String,
    /// True iff no Reject findings.
    pub admitted: bool,
    /// All findings, deterministic order (check, then span).
    pub findings: Vec<Finding>,
    /// Per-check wall timings.
    pub timings: Vec<CheckTiming>,
    /// Exact binding between submitted and authority-checked semantics.
    pub lowering: LoweringReceipt,
}

impl AdmissionReport {
    /// A canonical, deterministic diagnosis rendering (logs + G0 test).
    #[must_use]
    pub fn diagnosis(&self) -> String {
        let mut out = format!(
            "admission {} for {:?}: {} finding(s)\n",
            if self.admitted { "PASSED" } else { "REJECTED" },
            self.study,
            self.findings.len()
        );
        for f in &self.findings {
            let _ = writeln!(
                out,
                "[{}] {} @{}..{}: {}",
                f.check,
                match f.severity {
                    Severity::Warn => "warn",
                    Severity::Reject => "REJECT",
                },
                f.span.start,
                f.span.end,
                f.what
            );
            for (rank, fix) in f.fixes.iter().enumerate() {
                let wall = fix
                    .predicted_wall_s
                    .map_or(String::new(), |w| format!(" (predicted wall {w:.1}s)"));
                let _ = writeln!(
                    out,
                    "  fix#{rank}: {}{wall} — {}",
                    fix.action, fix.qoi_impact
                );
            }
        }
        out
    }
}

/// Bind a parsed source program to the current IR version, lower it, then run
/// admission over only the explicit lowered semantics.
///
/// This compatibility entry point is for newly parsed source. Persisted and
/// replayed artifacts should call [`admit_versioned`] so unsupported versions
/// are refused at [`VersionedProgram`] construction rather than inferred.
#[must_use]
pub fn admit(node: &Node, cx: &AdmissionContext<'_>) -> AdmissionReport {
    let program = match VersionedProgram::try_current(node.clone()) {
        Ok(program) => program,
        Err(error) => return lowering_refusal(LoweringReceipt::unbound_current(), error),
    };
    admit_versioned(&program, cx)
}

/// Lower and admit an explicitly version-bound FrankenScript program.
///
/// Lowering is the first decision boundary. A malformed shorthand returns one
/// `lowering` rejection without running capability, budget, chart, regime, or
/// other semantic checks over raw or partially lowered syntax.
#[must_use]
pub fn admit_versioned(program: &VersionedProgram, cx: &AdmissionContext<'_>) -> AdmissionReport {
    let mut receipt = LoweringReceipt::for_program(program);
    let lowered = match lower(program.program()) {
        Ok(lowered) => lowered,
        Err(error) => return lowering_refusal(receipt, error),
    };
    if let Err(error) = receipt.bind_lowered(&lowered.node) {
        return lowering_refusal(receipt, error);
    }
    admit_lowered(&lowered.node, cx, receipt)
}

fn lowering_refusal(receipt: LoweringReceipt, error: IrError) -> AdmissionReport {
    AdmissionReport {
        study: "<lowering-refused>".to_string(),
        admitted: false,
        findings: vec![Finding {
            check: "lowering",
            severity: Severity::Reject,
            span: error.span,
            what: error.detail,
            fixes: vec![RankedFix {
                action: error.hint,
                predicted_wall_s: None,
                qoi_impact: "structural fix; no QoI impact".to_string(),
            }],
        }],
        timings: Vec::new(),
        lowering: receipt,
    }
}

fn admit_lowered(
    node: &Node,
    cx: &AdmissionContext<'_>,
    receipt: LoweringReceipt,
) -> AdmissionReport {
    let mut findings: Vec<Finding> = Vec::new();
    let mut timings = Vec::new();
    let study = match Study::from_node(node) {
        Ok(s) => s,
        Err(e) => {
            return AdmissionReport {
                study: "<unparsed>".to_string(),
                admitted: false,
                findings: vec![Finding {
                    check: "structure",
                    severity: Severity::Reject,
                    span: e.span,
                    what: e.detail,
                    fixes: vec![RankedFix {
                        action: e.hint,
                        predicted_wall_s: None,
                        qoi_impact: "structural fix; no QoI impact".to_string(),
                    }],
                }],
                timings,
                lowering: receipt,
            };
        }
    };
    let mut run =
        |name: &'static str, f: &mut dyn FnMut(&mut Vec<Finding>), findings: &mut Vec<Finding>| {
            let t0 = Instant::now();
            f(findings);
            timings.push(CheckTiming {
                check: name,
                micros: t0.elapsed().as_micros(),
            });
        };
    run(
        "explicits",
        &mut |f| check_explicits(node, &study, cx, f),
        &mut findings,
    );
    run(
        "dimensional",
        &mut |f| check_dimensional(&study, f),
        &mut findings,
    );
    run(
        "budget",
        &mut |f| check_budget(&study, cx, f),
        &mut findings,
    );
    run(
        "capability",
        &mut |f| check_capability(&study, cx, f),
        &mut findings,
    );
    run("charts", &mut |f| check_charts(cx, f), &mut findings);
    run(
        "regime",
        &mut |f| check_regime(&study, cx, f),
        &mut findings,
    );
    findings.sort_by(|a, b| a.check.cmp(b.check).then(a.span.start.cmp(&b.span.start)));
    AdmissionReport {
        study: study.name.to_string(),
        admitted: !findings.iter().any(|f| f.severity == Severity::Reject),
        findings,
        timings,
        lowering: receipt,
    }
}

fn reject(check: &'static str, span: Span, what: String, action: String) -> Finding {
    Finding {
        check,
        severity: Severity::Reject,
        span,
        what,
        fixes: vec![RankedFix {
            action,
            predicted_wall_s: None,
            qoi_impact: "structural fix; no QoI impact".to_string(),
        }],
    }
}

// ------------------------------------------------------- Five Explicits

fn check_explicits(
    node: &Node,
    study: &Study<'_>,
    cx: &AdmissionContext<'_>,
    out: &mut Vec<Finding>,
) {
    if study.seed.is_none() {
        out.push(reject(
            "explicits",
            node.span,
            "the seed pillar is missing".to_string(),
            "add (seed 0x…) — every stochastic draw keys off it".to_string(),
        ));
    }
    match study.versions {
        None => out.push(reject(
            "explicits",
            node.span,
            "the versions pillar is missing".to_string(),
            "add (versions (constellation :lock \"…\"))".to_string(),
        )),
        Some(v) if study.constellation_lock().is_none() => out.push(reject(
            "explicits",
            v.span,
            "versions clause carries no constellation lock".to_string(),
            "pin (constellation :lock \"YYYY-MM\") inside (versions …)".to_string(),
        )),
        Some(_) => {}
    }
    if study.budget.is_none() {
        out.push(reject(
            "explicits",
            node.span,
            "the budget pillar is missing".to_string(),
            "add (budget (wall …) (mem …) …) or a qoi-precision budget".to_string(),
        ));
    }
    // Capabilities: satisfied by a session token OR an explicit clause.
    if study.capability.is_none() && cx.capability.is_none() {
        out.push(reject(
            "explicits",
            node.span,
            "no capability grant: neither a session token nor a (capability …) clause".to_string(),
            "attach the session's capability token to admission, or declare \
             (capability :cores … :mem … :wall … :ops (…))"
                .to_string(),
        ));
    }
}

// ---------------------------------------------------------- dimensional

const ARITH_SAME_DIMS: &[&str] = &["+", "-", "min", "max"];
const COMPARE: &[&str] = &["=", "<", ">", "<=", ">="];

fn infer_dims(
    node: &Node,
    env: &BTreeMap<&str, Option<Dims>>,
    out: &mut Vec<Finding>,
) -> Option<Dims> {
    match &node.kind {
        NodeKind::Int(_) | NodeKind::Float(_) => Some(Dims::NONE),
        NodeKind::Qty { dims, .. } => Some(*dims),
        NodeKind::Symbol(s) => env.get(s.as_str()).copied().flatten(),
        NodeKind::Count { .. } | NodeKind::Seed(_) | NodeKind::Str(_) | NodeKind::Keyword(_) => {
            None
        }
        NodeKind::List(items) => {
            let head = node.head();
            // A bare `()` has no head and no operands. `check_dimensional`
            // recurses `infer_dims` into every body clause and let-expression,
            // so an empty list is reachable from parseable input such as
            // `(study "x" ())`; the old `&items[1..]` was a usize range panic
            // (start 1 > len 0). Fail closed: treat it as dimensionless-unknown,
            // exactly like the other operand-free atoms above (Str/Keyword/Count
            // → None). `items.get(1..)` is `None` precisely when `items` is
            // empty, so `?` returns None (unknown dims) instead of panicking.
            let args = items.get(1..)?;
            if let Some(h) = head
                && (ARITH_SAME_DIMS.contains(&h) || COMPARE.contains(&h))
            {
                let mut known: Option<(Dims, Span)> = None;
                for a in args {
                    let d = infer_dims(a, env, out);
                    if let Some(d) = d {
                        match known {
                            None => known = Some((d, a.span)),
                            Some((expect, first_span)) if expect != d => {
                                out.push(Finding {
                                    check: "dimensional",
                                    severity: Severity::Reject,
                                    span: a.span,
                                    what: format!(
                                        "({h} …): operand dims {:?} disagree with {:?} \
                                         (first operand at {}..{})",
                                        d.0, expect.0, first_span.start, first_span.end
                                    ),
                                    fixes: vec![RankedFix {
                                        action: format!(
                                            "express both operands of ({h} …) in the same \
                                             physical dimensions"
                                        ),
                                        predicted_wall_s: None,
                                        qoi_impact: "unit error; the study is meaningless \
                                                     until fixed"
                                            .to_string(),
                                    }],
                                });
                            }
                            Some(_) => {}
                        }
                    }
                }
                return if COMPARE.contains(&head.unwrap_or_default()) {
                    Some(Dims::NONE)
                } else {
                    known.map(|(d, _)| d)
                };
            }
            if head == Some("*") || head == Some("/") {
                let mut acc = Dims::NONE;
                let mut all_known = true;
                for (i, a) in args.iter().enumerate() {
                    match infer_dims(a, env, out) {
                        Some(d) => {
                            // AUTHORITATIVE checked composition (bead
                            // sj31i.11): a saturated exponent could
                            // cancel back into range as false physics,
                            // so overflow REJECTS instead of clamping.
                            let composed = if head == Some("/") && i > 0 {
                                acc.checked_minus(d)
                            } else {
                                acc.checked_plus(d)
                            };
                            match composed {
                                Some(next) => acc = next,
                                None => {
                                    out.push(Finding {
                                        check: "dimensional",
                                        severity: Severity::Reject,
                                        span: a.span,
                                        what: format!(
                                            "({} …): operand dims {:?} overflow the supported \
                                             i8 exponent domain when composed with {:?}",
                                            head.unwrap_or_default(),
                                            d.0,
                                            acc.0
                                        ),
                                        fixes: vec![RankedFix {
                                            action: format!(
                                                "reduce the exponent magnitude of the operands \
                                                 of ({} …); admitted exponents must stay within \
                                                 i8",
                                                head.unwrap_or_default()
                                            ),
                                            predicted_wall_s: None,
                                            qoi_impact: "dimension overflow; the expression \
                                                         cannot carry admissible physics"
                                                .to_string(),
                                        }],
                                    });
                                    return None;
                                }
                            }
                        }
                        None => all_known = false,
                    }
                }
                return all_known.then_some(acc);
            }
            // Unknown verb: recurse for nested errors, result unknown.
            for a in args {
                let _ = infer_dims(a, env, out);
            }
            None
        }
    }
}

fn check_dimensional(study: &Study<'_>, out: &mut Vec<Finding>) {
    let mut env: BTreeMap<&str, Option<Dims>> = BTreeMap::new();
    for (name, expr) in &study.lets {
        let d = infer_dims(expr, &env, out);
        env.insert(name, d);
    }
    for clause in &study.body {
        let _ = infer_dims(clause, &env, out);
    }
}

// --------------------------------------------------------------- budget

fn qty_seconds(node: &Node) -> Option<f64> {
    if let NodeKind::Qty { value, dims, .. } = &node.kind
        && *dims == Dims([0, 0, 1, 0, 0, 0])
    {
        return Some(*value);
    }
    None
}

/// EXACT integral bytes for a Count node (gp3.20): u128-checked
/// scaling for exact literals (overflow refuses BEFORE rounding),
/// separately-defined fractional semantics for decimal forms. `None`
/// for non-Count nodes, Cores, fractional/overflowing results.
fn count_bytes_exact(node: &Node) -> Option<u64> {
    if let NodeKind::Count { value, unit } = &node.kind {
        return value.integral_bytes(*unit);
    }
    None
}

fn invalid_resource_value(
    check: &'static str,
    span: Span,
    resource: &str,
    requirement: &str,
) -> Finding {
    reject(
        check,
        span,
        format!("{resource} must be {requirement}"),
        format!("supply {resource} as {requirement}"),
    )
}

fn check_budget_resource_domains(study: &Study<'_>, out: &mut Vec<Finding>) {
    let Some(budget) = study.budget else {
        return;
    };
    let Some(items) = budget.items() else {
        return;
    };
    if items.len() == 1 {
        out.push(reject(
            "budget",
            budget.span,
            "the budget pillar is empty".to_string(),
            "declare at least one wall, memory, or QoI accuracy budget".to_string(),
        ));
        return;
    }

    let mut seen = BTreeSet::new();
    for clause in &items[1..] {
        let Some(values) = clause.items() else {
            out.push(reject(
                "budget",
                clause.span,
                "budget entries must be parenthesized clauses".to_string(),
                "write entries such as (wall 10s), (mem 1GiB), or (qoi ...)".to_string(),
            ));
            continue;
        };
        let Some(resource) = clause.head() else {
            out.push(reject(
                "budget",
                clause.span,
                "budget entries must have a symbolic name and a value".to_string(),
                "remove the empty entry or name its budget dimension".to_string(),
            ));
            continue;
        };
        if values.len() < 2 {
            out.push(reject(
                "budget",
                clause.span,
                format!("{resource} budget entry has no value"),
                format!("supply a value in ({resource} ...)"),
            ));
            continue;
        }
        if !matches!(resource, "wall" | "mem") {
            // Operator-specific and QoI budgets remain extensible until the
            // catalog lands, but they must still be named, structured clauses
            // with at least one value.
            continue;
        }
        if !seen.insert(resource) {
            out.push(reject(
                "budget",
                clause.span,
                format!("duplicate {resource} budget is ambiguous"),
                format!("retain exactly one ({resource} ...) budget"),
            ));
            continue;
        }
        if values.len() != 2 {
            out.push(reject(
                "budget",
                clause.span,
                format!("{resource} budget takes exactly one value"),
                format!("write ({resource} <value>) with no extra operands"),
            ));
            continue;
        }
        let value = values.get(1);
        match resource {
            "wall" => {
                let seconds = value.and_then(qty_seconds);
                if !seconds.is_some_and(|seconds| seconds.is_finite() && seconds > 0.0) {
                    out.push(invalid_resource_value(
                        "budget",
                        clause.span,
                        "wall budget",
                        "a finite positive time quantity",
                    ));
                }
            }
            "mem" => {
                let bytes = value.and_then(count_bytes_exact);
                if bytes.is_none_or(|bytes| bytes == 0) {
                    out.push(invalid_resource_value(
                        "budget",
                        clause.span,
                        "memory budget",
                        "a finite positive whole-byte count below 2^64",
                    ));
                }
            }
            _ => unreachable!("resource filtered above"),
        }
    }
}

/// The declared wall budget: `(budget (wall 2h) …)`, else the capability
/// clause's `:wall`, else the session token's grant.
fn wall_budget_s(study: &Study<'_>, cx: &AdmissionContext<'_>) -> Option<f64> {
    if let Some(budget) = study.budget
        && let Some(items) = budget.items()
    {
        for clause in items {
            if clause.head() == Some("wall")
                && let Some(w) = clause.items().and_then(|i| i.get(1)).and_then(qty_seconds)
            {
                return Some(w);
            }
        }
    }
    if let Some(cap) = study.capability
        && let Some(items) = cap.items()
    {
        for pair in items.windows(2) {
            if let NodeKind::Keyword(k) = &pair[0].kind
                && k == "wall"
                && let Some(w) = qty_seconds(&pair[1])
            {
                return Some(w);
            }
        }
    }
    cx.capability.as_ref().map(|c| c.wall_s)
}

/// Size feature for a verb call: exactly one numeric
/// `:dof`/`:size`/`:modes` argument, else the unit-size default.
fn size_of_call(items: &[Node], verb: &str) -> Result<f64, String> {
    let mut size = None;
    for (index, item) in items.iter().enumerate() {
        let NodeKind::Keyword(keyword) = &item.kind else {
            continue;
        };
        if !matches!(keyword.as_str(), "dof" | "size" | "modes") {
            continue;
        }
        if size.is_some() {
            return Err(format!(
                "operation {verb:?} declares more than one :dof/:size/:modes feature"
            ));
        }
        let value = items
            .get(index + 1)
            .ok_or_else(|| format!("operation {verb:?} has no value after :{keyword}"))?;
        let value = match &value.kind {
            #[allow(clippy::cast_precision_loss)]
            NodeKind::Int(value) => *value as f64,
            NodeKind::Float(value) => *value,
            _ => {
                return Err(format!(
                    "operation {verb:?} requires a numeric value after :{keyword}"
                ));
            }
        };
        if !value.is_finite() || value < 0.0 {
            return Err(format!(
                "operation {verb:?} requires a finite non-negative value after :{keyword}"
            ));
        }
        size = Some(value);
    }
    Ok(size.unwrap_or(1.0))
}

/// Collect (verb, size, span) for every modeled call in the tree.
type ModeledCall<'a> = (&'a str, f64, Span);
type CostedCall<'a> = (&'a str, f64, f64, Span);

fn modeled_calls<'a>(
    node: &'a Node,
    models: &BTreeMap<String, SealedCostModel>,
    out: &mut Vec<ModeledCall<'a>>,
    findings: &mut Vec<Finding>,
) {
    if let NodeKind::List(items) = &node.kind {
        if let Some(h) = node.head()
            && (h.contains('.') || models.contains_key(h))
        {
            match size_of_call(items, h) {
                Ok(size) if models.contains_key(h) => out.push((h, size, node.span)),
                Ok(_) => {}
                Err(what) => findings.push(reject(
                    "budget",
                    node.span,
                    what,
                    "supply at most one numeric :dof, :size, or :modes value".to_string(),
                )),
            }
        }
        for child in items {
            modeled_calls(child, models, out, findings);
        }
    }
}

fn cost_model_rejection(span: Span, what: String, action: String) -> Finding {
    Finding {
        check: "budget",
        severity: Severity::Reject,
        span,
        what,
        fixes: vec![RankedFix {
            action,
            predicted_wall_s: None,
            qoi_impact: "evidence acquisition only; no modeled QoI change".to_string(),
        }],
    }
}

/// One admitting-with-notice finding for weakened cost evidence
/// (beads 2pmb + l2k92): provisional fits and stale receipts both
/// inform the budget arithmetic, but the admission record names the
/// class — and for stale evidence the exact staleness verdict — at
/// the decision point, never silently.
fn weakened_evidence_notice(
    verb: &str,
    kernel: &str,
    class: CostEvidenceClass,
    verdict: StalenessVerdict,
    span: Span,
) -> Option<Finding> {
    if class.rank() >= CostEvidenceClass::ExactRooflineReceipt.rank() {
        return None;
    }
    let (label, cause, action) = match verdict {
        StalenessVerdict::AgedOut { age_ns, horizon_ns } => (
            "CostModelStale",
            format!("aged {age_ns} ns against a {horizon_ns} ns freshness horizon"),
            format!(
                "re-record {verb} tune evidence through the exact roofline lane; the retained \
                 receipt is older than the preregistered horizon"
            ),
        ),
        StalenessVerdict::MachineDrift => (
            "CostModelStale",
            "recorded under a different machine fingerprint".to_string(),
            format!("re-record {verb} tune evidence on the current machine"),
        ),
        StalenessVerdict::FutureRecording { ahead_ns } => (
            "CostModelStale",
            format!("recorded {ahead_ns} ns in the future of the observation time"),
            format!(
                "repair the {verb} evidence timeline; a future recording is clock skew or \
                 forged provenance"
            ),
        ),
        StalenessVerdict::Fresh | StalenessVerdict::NotApplicable => (
            "CostModelProvisional",
            "not a validated roofline receipt".to_string(),
            format!(
                "record {verb} tune evidence through the exact roofline lane so \
                 fs-plan::cost_model_from_tune can mint receipt-backed authority"
            ),
        ),
    };
    Some(Finding {
        check: "budget",
        severity: Severity::Warn,
        span,
        what: format!(
            "{label}: the wall-cost model for operation {verb:?} is {} evidence \
             (scope {kernel:?}): {cause}",
            class.name(),
        ),
        fixes: vec![RankedFix {
            action,
            predicted_wall_s: None,
            qoi_impact: "evidence provenance only; no modeled QoI change".to_string(),
        }],
    })
}

fn cost_registered_calls<'a>(
    calls: &[ModeledCall<'a>],
    models: &BTreeMap<String, SealedCostModel>,
    freshness: Option<&CostFreshnessContext>,
    out: &mut Vec<Finding>,
) -> Option<(f64, Vec<CostedCall<'a>>)> {
    let mut total = 0.0f64;
    let mut costed = Vec::with_capacity(calls.len());
    let mut refused = false;
    let mut evidence_warned: Vec<&str> = Vec::new();
    for &(verb, size, span) in calls {
        let model = models
            .get(verb)
            .expect("modeled_calls only returns registered models");
        // Sealed authority (beads 2pmb + l2k92): weakened evidence may
        // inform the budget arithmetic, but the admission record must
        // say so — once per verb, admitting with notice, never
        // silently. With a freshness contract attached, receipt-backed
        // evidence is ASSESSED (age, machine drift, future stamps)
        // before the notice decision.
        let (class, verdict) = match freshness {
            Some(contract) => model.assess(contract.now_ns, &contract.machine, contract.policy),
            None => (model.evidence_class(), StalenessVerdict::NotApplicable),
        };
        if !evidence_warned.contains(&verb)
            && let Some(finding) =
                weakened_evidence_notice(verb, model.scope().kernel(), class, verdict, span)
        {
            evidence_warned.push(verb);
            out.push(finding);
        }
        match model.predict(size).map(|sealed| sealed.prediction) {
            Ok(prediction) if prediction.p90.is_finite() && prediction.p90 >= 0.0 => {
                let next_total = total + prediction.p90;
                if next_total.is_finite() {
                    total = next_total;
                    costed.push((verb, size, prediction.p90, span));
                } else {
                    refused = true;
                    out.push(cost_model_rejection(
                        span,
                        format!(
                            "CostModelInvalidPrediction: adding the registered model's p90 for operation {verb:?} at size {size} overflows finite wall time"
                        ),
                        format!(
                            "re-fit {verb} from finite tune observations at representative sizes"
                        ),
                    ));
                }
            }
            Ok(prediction) => {
                refused = true;
                out.push(cost_model_rejection(
                    span,
                    format!(
                        "CostModelInvalidPrediction: registered model for operation {verb:?} returned non-finite or negative p90 {:?} at size {size}",
                        prediction.p90
                    ),
                    format!("re-fit {verb} from valid finite tune observations"),
                ));
            }
            Err(reason) => {
                refused = true;
                out.push(cost_model_rejection(
                    span,
                    format!(
                        "CostModelRefused: registered model for operation {verb:?} at size {size} did not provide a wall-cost prediction: {reason}"
                    ),
                    format!(
                        "run {verb} at representative sizes and record enough tune observations before admission"
                    ),
                ));
            }
        }
    }
    (!refused).then_some((total, costed))
}

fn check_budget(study: &Study<'_>, cx: &AdmissionContext<'_>, out: &mut Vec<Finding>) {
    check_budget_resource_domains(study, out);
    let mut calls = Vec::new();
    for (_, expr) in &study.lets {
        modeled_calls(expr, &cx.cost_models, &mut calls, out);
    }
    for clause in &study.body {
        modeled_calls(clause, &cx.cost_models, &mut calls, out);
    }
    let Some(wall) = wall_budget_s(study, cx) else {
        return; // qoi-precision-only budgets carry no wall bound to screen
    };
    if calls.is_empty() {
        return;
    }
    let Some((total, mut costed)) =
        cost_registered_calls(&calls, &cx.cost_models, cx.cost_freshness.as_ref(), out)
    else {
        return;
    };
    if total <= wall {
        return;
    }
    // BudgetInfeasible with RANKED, cost-model-derived fixes (§11.3).
    costed.sort_by(|a, b| b.2.total_cmp(&a.2));
    let mut fixes = Vec::new();
    if let Some((verb, size, p90, _)) = costed.first() {
        let halved = cx
            .cost_models
            .get(*verb)
            .and_then(|model| model.predict(size / 2.0).ok())
            .map(|sealed| sealed.prediction.p90)
            .filter(|prediction| prediction.is_finite() && *prediction >= 0.0);
        if let Some(halved) = halved {
            fixes.push(RankedFix {
                action: format!(
                    "coarsen {verb}: halve its size feature ({size} -> {})",
                    size / 2.0
                ),
                predicted_wall_s: Some(total - p90 + halved),
                qoi_impact: "resolution halves; the verb's error model governs the \
                             QoI degradation"
                    .to_string(),
            });
        }
        fixes.push(RankedFix {
            action: format!(
                "surrogate-screen {verb} (evaluate candidates on the surrogate, \
                             verify winners at full fidelity)"
            ),
            predicted_wall_s: Some(total - p90 + 0.2 * p90),
            qoi_impact: "screening decisions carry surrogate error; final verdicts are \
                         re-verified"
                .to_string(),
        });
    }
    fixes.push(RankedFix {
        action: format!("relax the wall budget to {:.0}s", (total * 1.2).ceil()),
        predicted_wall_s: Some(total),
        qoi_impact: "no QoI impact; costs more wall".to_string(),
    });
    fixes.sort_by(|a, b| {
        a.predicted_wall_s
            .unwrap_or(f64::INFINITY)
            .total_cmp(&b.predicted_wall_s.unwrap_or(f64::INFINITY))
    });
    let span = costed.first().map_or(Span::default(), |c| c.3);
    out.push(Finding {
        check: "budget",
        severity: Severity::Reject,
        span,
        what: format!(
            "BudgetInfeasible: predicted p90 wall {total:.1}s exceeds the {wall:.1}s bound \
             ({} modeled op(s))",
            costed.len()
        ),
        fixes,
    });
}

// ----------------------------------------------------------- capability

fn glob_matches(pattern: &str, verb: &str) -> bool {
    pattern
        .strip_suffix('*')
        .map_or(pattern == verb, |prefix| verb.starts_with(prefix))
}

/// Whether an operator grant is one canonical exact name or a namespace
/// wildcard ending in `.*`.
///
/// Capability issuers use this same predicate as admission so malformed
/// wildcard spellings cannot acquire broader semantics before IR validation.
#[must_use]
pub fn valid_operator_pattern(pattern: &str) -> bool {
    if pattern.trim() != pattern
        || pattern.is_empty()
        || pattern.starts_with('.')
        || pattern.ends_with('.')
        || pattern.contains("..")
        || !pattern
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-' | b'*'))
    {
        return false;
    }
    match pattern.find('*') {
        None => true,
        Some(index) => {
            index >= 2
                && index + 1 == pattern.len()
                && pattern[..index].ends_with('.')
                && !pattern[..index - 1].contains('*')
        }
    }
}

fn glob_covers(grant: &str, requested: &str) -> bool {
    if let Some(requested_prefix) = requested.strip_suffix('*') {
        return grant
            .strip_suffix('*')
            .is_some_and(|grant_prefix| requested_prefix.starts_with(grant_prefix));
    }
    glob_matches(grant, requested)
}

fn namespaced_verbs<'a>(node: &'a Node, out: &mut Vec<(&'a str, Span)>) {
    if let NodeKind::List(items) = &node.kind {
        if let Some(h) = node.head()
            && h.contains('.')
        {
            out.push((h, node.span));
        }
        for child in items {
            namespaced_verbs(child, out);
        }
    }
}

#[allow(clippy::too_many_lines)] // One ordered, fail-closed capability admission matrix.
fn check_capability(study: &Study<'_>, cx: &AdmissionContext<'_>, out: &mut Vec<Finding>) {
    let token = cx.capability.as_ref();
    if let Some(token) = token {
        if !token.wall_s.is_finite() || token.wall_s < 0.0 {
            out.push(invalid_resource_value(
                "capability",
                Span::default(),
                "session wall grant",
                "finite and non-negative",
            ));
        }
        if token
            .ops
            .iter()
            .any(|pattern| !valid_operator_pattern(pattern))
        {
            out.push(reject(
                "capability",
                Span::default(),
                "the session token contains an invalid operator grant pattern".to_string(),
                "use exact operator names or namespace wildcards such as flux.*".to_string(),
            ));
        }
    }

    let mut declared_ops: Option<Vec<&str>> = None;
    let mut declared_field_count = 0usize;
    if let Some(capability) = study.capability
        && let Some(items) = capability.items()
    {
        let mut seen = BTreeSet::new();
        let fields = &items[1..];
        if !fields.len().is_multiple_of(2) {
            out.push(reject(
                "capability",
                capability.span,
                "capability fields must be exact keyword/value pairs".to_string(),
                "remove the dangling field or supply its value".to_string(),
            ));
        }
        for pair in fields.as_chunks::<2>().0 {
            let field_node = &pair[0];
            let value_node = &pair[1];
            let NodeKind::Keyword(field) = &field_node.kind else {
                out.push(reject(
                    "capability",
                    field_node.span,
                    "capability field names must be keywords".to_string(),
                    "use :cores, :mem, :wall, or :ops before each value".to_string(),
                ));
                continue;
            };
            if !matches!(field.as_str(), "cores" | "mem" | "wall" | "ops") {
                out.push(reject(
                    "capability",
                    field_node.span,
                    format!("unknown capability field :{field}"),
                    "use only :cores, :mem, :wall, and :ops".to_string(),
                ));
                continue;
            }
            declared_field_count += 1;
            if !seen.insert(field.as_str()) {
                out.push(reject(
                    "capability",
                    field_node.span,
                    format!("duplicate :{field} capability field is ambiguous"),
                    format!("retain exactly one :{field} field"),
                ));
                continue;
            }
            match field.as_str() {
                "cores" => {
                    let value = match &value_node.kind {
                        NodeKind::Int(value) => u64::try_from(*value).ok(),
                        NodeKind::Count {
                            value,
                            unit: CountUnit::Cores,
                        } => value.integral_count(),
                        _ => None,
                    };
                    if value.is_none() {
                        out.push(invalid_resource_value(
                            "capability",
                            value_node.span,
                            "core ask",
                            "a non-negative whole core count fitting u64",
                        ));
                    }
                }
                "mem" => {
                    if count_bytes_exact(value_node).is_none() {
                        out.push(invalid_resource_value(
                            "capability",
                            value_node.span,
                            "memory ask",
                            "a finite non-negative whole-byte count below 2^64",
                        ));
                    }
                }
                "wall" => {
                    let value = qty_seconds(value_node);
                    if !value.is_some_and(|value| value.is_finite() && value >= 0.0) {
                        out.push(invalid_resource_value(
                            "capability",
                            value_node.span,
                            "wall ask",
                            "a finite non-negative time quantity",
                        ));
                    }
                }
                "ops" => match &value_node.kind {
                    NodeKind::List(nodes)
                        if !nodes.is_empty()
                            && nodes.iter().all(|node| {
                                matches!(&node.kind, NodeKind::Symbol(pattern) if valid_operator_pattern(pattern))
                            }) =>
                    {
                        declared_ops = Some(
                            nodes
                                .iter()
                                .filter_map(|node| match &node.kind {
                                    NodeKind::Symbol(pattern) => Some(pattern.as_str()),
                                    _ => None,
                                })
                                .collect(),
                        );
                    }
                    _ => out.push(reject(
                        "capability",
                        value_node.span,
                        ":ops must be a non-empty list of exact names or namespace wildcards"
                            .to_string(),
                        "declare operator grants such as :ops (flux.* ascent.optimize)"
                            .to_string(),
                    )),
                },
                _ => unreachable!("capability field filtered above"),
            }
        }
    }

    if token.is_none() {
        if study.capability.is_some() && declared_field_count == 0 {
            out.push(reject(
                "capability",
                study
                    .capability
                    .map_or(Span::default(), |capability| capability.span),
                "the explicit capability pillar contains no recognized grants".to_string(),
                "declare resource limits and a non-empty :ops grant".to_string(),
            ));
        }
        if study.capability.is_some() && declared_ops.is_none() {
            out.push(reject(
                "capability",
                study
                    .capability
                    .map_or(Span::default(), |capability| capability.span),
                "a self-contained capability pillar has no valid :ops grant".to_string(),
                "add a non-empty :ops list or attach a session capability token".to_string(),
            ));
        }
    }

    let mut verbs = Vec::new();
    for (_, expr) in &study.lets {
        namespaced_verbs(expr, &mut verbs);
    }
    for clause in &study.body {
        namespaced_verbs(clause, &mut verbs);
    }
    for (verb, span) in verbs {
        if let Some(token) = token
            && !token.ops.iter().any(|pattern| glob_matches(pattern, verb))
        {
            out.push(Finding {
                check: "capability",
                severity: Severity::Reject,
                span,
                what: format!("operator {verb} is outside the session token's grants"),
                fixes: vec![RankedFix {
                    action: format!(
                        "request a token covering {}.* or remove the {verb} op",
                        verb.split('.').next().unwrap_or(verb)
                    ),
                    predicted_wall_s: None,
                    qoi_impact: "capability change; no QoI impact".to_string(),
                }],
            });
        }
        if let Some(patterns) = &declared_ops
            && !patterns.iter().any(|pattern| glob_matches(pattern, verb))
        {
            out.push(Finding {
                check: "capability",
                severity: Severity::Reject,
                span,
                what: format!("operator {verb} is outside the study's explicit capability grants"),
                fixes: vec![RankedFix {
                    action: format!(
                        "add {}.* to :ops or remove the {verb} op",
                        verb.split('.').next().unwrap_or(verb)
                    ),
                    predicted_wall_s: None,
                    qoi_impact: "capability change; no QoI impact".to_string(),
                }],
            });
        }
    }

    if let (Some(token), Some(patterns)) = (token, &declared_ops) {
        for requested in patterns {
            if !token.ops.iter().any(|grant| glob_covers(grant, requested)) {
                out.push(reject(
                    "capability",
                    study
                        .capability
                        .map_or(Span::default(), |capability| capability.span),
                    format!("declared operator grant {requested:?} exceeds the session token"),
                    "remove the overbroad declared grant or obtain a covering token".to_string(),
                ));
            }
        }
    }

    // Declared asks must fit inside the token.
    if let Some(token) = token
        && let Some(cap) = study.capability
        && let Some(items) = cap.items()
    {
        for pair in items.windows(2) {
            let NodeKind::Keyword(k) = &pair[0].kind else {
                continue;
            };
            let over = match k.as_str() {
                "cores" => match &pair[1].kind {
                    NodeKind::Int(value) => u64::try_from(*value).ok().and_then(|cores| {
                        (cores > token.cores)
                            .then(|| format!("{cores} cores asked, {} granted", token.cores))
                    }),
                    NodeKind::Count {
                        value,
                        unit: CountUnit::Cores,
                    } => value.integral_count().and_then(|cores| {
                        (cores > token.cores)
                            .then(|| format!("{cores} cores asked, {} granted", token.cores))
                    }),
                    _ => None,
                },
                "mem" => count_bytes_exact(&pair[1]).and_then(|bytes| {
                    (bytes > token.mem_bytes)
                        .then(|| format!("{bytes} bytes asked, {} bytes granted", token.mem_bytes))
                }),
                "wall" => qty_seconds(&pair[1]).and_then(|w| {
                    (w > token.wall_s)
                        .then(|| format!("{w:.0}s wall asked, {:.0}s granted", token.wall_s))
                }),
                _ => None,
            };
            if let Some(what) = over {
                out.push(Finding {
                    check: "capability",
                    severity: Severity::Reject,
                    span: pair[1].span,
                    what: format!("capability ask exceeds the session token: {what}"),
                    fixes: vec![RankedFix {
                        action: "lower the ask or obtain a larger session token".to_string(),
                        predicted_wall_s: None,
                        qoi_impact: "resource change; no QoI impact".to_string(),
                    }],
                });
            }
        }
    }
}

// --------------------------------------------------------------- charts

fn check_charts(cx: &AdmissionContext<'_>, out: &mut Vec<Finding>) {
    if cx.chart_requirements.is_empty() {
        return;
    }
    let Some((router, oracle)) = cx.router else {
        out.push(Finding {
            check: "charts",
            severity: Severity::Warn,
            span: Span::default(),
            what: format!(
                "{} chart requirement(s) declared but no Router attached — feasibility \
                 unverified",
                cx.chart_requirements.len()
            ),
            fixes: vec![RankedFix {
                action: "attach the Rep Router (and cost oracle) to the admission context"
                    .to_string(),
                predicted_wall_s: None,
                qoi_impact: "verification gap only".to_string(),
            }],
        });
        return;
    };
    for req in &cx.chart_requirements {
        let request = RouteRequest {
            from: req.from.clone(),
            to: req.to.clone(),
            scale: req.scale,
            max_abs_error: req.max_abs_error,
            max_cost_s: req.max_cost_s,
        };
        match router.plan(&request, oracle) {
            Ok(plan) if !plan.all_certified() => out.push(Finding {
                check: "charts",
                severity: Severity::Reject,
                span: Span::default(),
                what: format!(
                    "conversion route {} -> {} is only estimated; at least one converter is not \
                     declared certificate-backed",
                    req.from, req.to
                ),
                fixes: vec![RankedFix {
                    action: "register and validate a certificate-backed converter chain, or keep \
                             this study outside authoritative admission"
                        .to_string(),
                    predicted_wall_s: Some(plan.predicted_cost_s()),
                    qoi_impact: format!(
                        "estimated route has declared composed error {:.3e}",
                        plan.composed_abs_error()
                    ),
                }],
            }),
            Ok(_) => {}
            Err(error) => {
                let (what, fixes) = match error {
                    RoutePlanError::Infeasible(refusal) => {
                        let fixes = refusal
                            .fixes
                            .iter()
                            .map(|fix| RankedFix {
                                action: fix.clone(),
                                predicted_wall_s: refusal.best_cost_s,
                                qoi_impact: refusal.best_abs_error.map_or_else(
                                    || "route feasibility change".to_string(),
                                    |value| format!("best achievable composed error {value:.3e}"),
                                ),
                            })
                            .collect();
                        (
                            format!("no conversion route {} -> {}: {refusal}", req.from, req.to),
                            fixes,
                        )
                    }
                    invalid => (
                        format!(
                            "conversion route authority {} -> {} is invalid: {invalid}",
                            req.from, req.to
                        ),
                        vec![RankedFix {
                            action:
                                "repair the route request or measured cost/error oracle evidence"
                                    .to_string(),
                            predicted_wall_s: None,
                            qoi_impact: "invalid routing authority cannot support admission"
                                .to_string(),
                        }],
                    ),
                };
                out.push(Finding {
                    check: "charts",
                    severity: Severity::Reject,
                    span: Span::default(),
                    what,
                    fixes,
                });
            }
        }
    }
}

// --------------------------------------------------------------- regime

fn regime_asserts<'a>(node: &'a Node, out: &mut Vec<(&'a str, Span)>) {
    if let NodeKind::List(items) = &node.kind {
        if node.head() == Some("assert")
            && let Some(inner) = items.get(1)
            && inner.head() == Some("regime.allows")
            && let Some(model) = inner.items().and_then(|i| i.get(1))
            && let NodeKind::Symbol(name) = &model.kind
        {
            out.push((name, inner.span));
        }
        for child in items {
            regime_asserts(child, out);
        }
    }
}

fn check_regime(study: &Study<'_>, cx: &AdmissionContext<'_>, out: &mut Vec<Finding>) {
    let severity = match cx.regime_policy {
        RegimePolicy::Reject => Severity::Reject,
        RegimePolicy::Warn => Severity::Warn,
    };
    // Explicit asserts + flux.* solver verbs both face the report.
    let mut asks: Vec<(&str, Span)> = Vec::new();
    let mut verbs = Vec::new();
    for (_, expr) in &study.lets {
        regime_asserts(expr, &mut asks);
        namespaced_verbs(expr, &mut verbs);
    }
    for clause in &study.body {
        regime_asserts(clause, &mut asks);
        namespaced_verbs(clause, &mut verbs);
    }
    asks.extend(verbs.into_iter().filter(|(v, _)| v.starts_with("flux.")));
    if asks.is_empty() {
        return;
    }
    let Some(report) = cx.regime else {
        out.push(Finding {
            check: "regime",
            severity: Severity::Warn,
            span: asks[0].1,
            what: "solver choices present but no RegimeReport attached — regime validity \
                   unverified"
                .to_string(),
            fixes: vec![RankedFix {
                action: "run fs_regime::assess over the study's physical inputs and attach \
                         the report"
                    .to_string(),
                predicted_wall_s: None,
                qoi_impact: "verification gap only".to_string(),
            }],
        });
        return;
    };
    for (model, span) in asks {
        if report.valid_models.iter().any(|m| m == model) {
            continue;
        }
        if let Some((_, reason)) = report.invalid_models.iter().find(|(m, _)| m == model) {
            let mut fixes: Vec<RankedFix> = report
                .valid_models
                .iter()
                .map(|alt| RankedFix {
                    action: format!("switch to {alt} (valid in this regime)"),
                    predicted_wall_s: None,
                    qoi_impact: "model change inside the validity domain".to_string(),
                })
                .collect();
            fixes.truncate(3);
            out.push(Finding {
                check: "regime",
                severity,
                span,
                what: format!(
                    "{model} is outside its validity domain here ({reason}); dominant \
                     balance: {}",
                    report.dominant_balance
                ),
                fixes,
            });
        }
        // Unknown to the registry: silence (cards land with their solvers).
    }
}

#[cfg(test)]
mod evidence_notice_tests {
    use super::*;

    fn notice(class: CostEvidenceClass, verdict: StalenessVerdict) -> Option<Finding> {
        weakened_evidence_notice("gemm", "gemm-f64", class, verdict, Span::default())
    }

    #[test]
    fn fresh_exact_evidence_is_silent() {
        assert!(
            notice(
                CostEvidenceClass::ExactRooflineReceipt,
                StalenessVerdict::Fresh
            )
            .is_none()
        );
    }

    #[test]
    fn provisional_evidence_admits_with_the_2pmb_notice() {
        let finding = notice(
            CostEvidenceClass::ProvisionalUnaudited,
            StalenessVerdict::NotApplicable,
        )
        .expect("provisional warns");
        assert_eq!(finding.severity, Severity::Warn);
        assert!(finding.what.starts_with("CostModelProvisional"));
        assert!(finding.what.contains("provisional-unaudited"));
    }

    #[test]
    fn stale_verdicts_name_their_exact_cause_at_the_decision_point() {
        let aged = notice(
            CostEvidenceClass::StaleRooflineReceipt,
            StalenessVerdict::AgedOut {
                age_ns: 900,
                horizon_ns: 500,
            },
        )
        .expect("aged-out warns");
        assert!(aged.what.starts_with("CostModelStale"));
        assert!(
            aged.what
                .contains("aged 900 ns against a 500 ns freshness horizon")
        );
        assert!(aged.what.contains("stale-roofline-receipt"));

        let drift = notice(
            CostEvidenceClass::StaleRooflineReceipt,
            StalenessVerdict::MachineDrift,
        )
        .expect("drift warns");
        assert!(drift.what.contains("different machine fingerprint"));
        assert!(drift.fixes[0].action.contains("current machine"));

        let future = notice(
            CostEvidenceClass::StaleRooflineReceipt,
            StalenessVerdict::FutureRecording { ahead_ns: 42 },
        )
        .expect("future warns");
        assert!(future.what.contains("recorded 42 ns in the future"));
        assert!(future.fixes[0].action.contains("clock skew or"));
    }
}
