//! Infeasibility DIAGNOSIS: elastic-relaxation solves, minimal unsat
//! cores by deletion filtering, and RANKED repairs with feasibility
//! estimates — the machinery that turns "optimizer failed" into a
//! design conversation.
//!
//! The elastic solve minimizes total hinge violation `Σ max(gᵢ, 0)`
//! over a domain box with multi-start projected subgradient descent
//! (deterministic LCG starts). Feasible ⟺ the elastic optimum's total
//! violation is ~0. The unsat core starts from the elastic support
//! (violated constraints at the optimum) and is refined by the
//! DELETION FILTER. The support is verified jointly infeasible before
//! filtering; otherwise the seed expands to the full, already-proven
//! infeasible set. The result is MINIMAL: dropping ANY member restores
//! feasibility — a property the conformance battery checks against
//! brute-force enumeration.
//!
//! RESOURCE CONTRACT (bead frankensim-constraint-restoration-budget-
//! receipts-x5sev): every entry point first builds a checked
//! [`RestorationWorkPlan`] and admits it against the caller's `Cx`
//! budget ([`AdmittedBudget`]) BEFORE the first charged tile. All work
//! then flows through one shared accountant: canonical skip-mask
//! builds, total-violation passes, finite-difference probes, descent
//! steps and starts, deletion-filter subset solves, and Monte-Carlo
//! repair sampling are checkpointed at deterministic bounded tiles
//! (`CANCELLATION_STRIDE`) and charged with contract unit weights. A
//! cancellation, deadline, poll-quota, or cost-quota stop drains
//! immediately with a typed [`RestorationError::Refused`] carrying the
//! exact [`RestorationWorkReceipt`] — never a partial report. Input
//! faults (malformed domains, duplicate/out-of-range skip indices,
//! inconsistent plans) refuse as [`RestorationError::Invalid`] before
//! any budget authority is consumed. Buffers proportional to
//! dimensions/constraints are allocated without lease admission on
//! lease-less contexts; the receipt records
//! [`RestorationMemoryAuthority::NoLeaseNoClaim`] rather than inventing
//! private memory authority.
//!
//! Receipt granularity: `work_units_charged`, `consumption`, and
//! `memory` describe the WHOLE admitted run (the shared accountant);
//! `starts_completed` reports the PRIMARY elastic solve only, matching
//! the embedded [`ElasticReport`]'s solver-scoped `evals`.

use crate::restoration::{
    RestorationError, RestorationMemoryAuthority, RestorationWorkLimits, RestorationWorkPlan,
    RestorationWorkReceipt, RestorationWorkShape, RESTORATION_UNITS_CONSTRAINT_EVALUATION,
    RESTORATION_UNITS_SKIP_MASK_ENTRY,
};
use crate::{ConError, ConstraintSpec, DomainError, DomainRangeError, push_json_string, scalar_at};
use fs_exec::{AdmittedBudget, BudgetConsumption, BudgetRefusal, Cx};
use fs_opt::{Manifold, Problem};

/// Per-component design-domain box.
#[derive(Debug, Clone, PartialEq)]
pub struct DomainBox {
    /// `(lo, hi)` per component of the sole `Rn` design variable. Admission
    /// requires exact dimension, finite ordered endpoints, and finite spans;
    /// `lo == hi` denotes a valid fixed coordinate.
    pub ranges: Vec<(f64, f64)>,
}

/// The elastic-relaxation solve's outcome.
#[derive(Debug, Clone)]
pub struct ElasticReport {
    /// The minimizer of total violation.
    pub x: Vec<f64>,
    /// Total hinge violation at the optimum (~0 ⟺ feasible).
    pub total_violation: f64,
    /// Per-constraint violations at the optimum.
    pub violations: Vec<f64>,
    /// Objective evaluations spent by THIS solve.
    pub evals: u64,
    /// Retained resource contract: which plan ran, what the shared
    /// accountant enforced, and how much work completed. On success the
    /// consumption carries NO refusal — success is claimable only when
    /// the admitted work ran to completion under its budget.
    pub work: RestorationWorkReceipt,
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }
}

/// Feasibility tolerance for the elastic optimum.
pub(crate) const FEAS_TOL: f64 = 1e-6;

/// Defensive penalty for an evaluator backend that carries a NON-FINITE scalar
/// instead of returning fs-opt's current typed `EvalNonFinite` refusal. It is
/// finite so a raw `NaN.max(0.0)` can never be dropped into false feasibility.
const NONFINITE_PENALTY: f64 = 1e30;

const CANCELLATION_STRIDE: usize = 64;

const PHASE_SKIP_MASK: &str = "fs-constraint:restoration-skip-mask";
const PHASE_TOTAL_PASS: &str = "fs-constraint:restoration-total-pass";
const PHASE_FD_PROBE: &str = "fs-constraint:restoration-fd-probe";
const PHASE_STEP: &str = "fs-constraint:restoration-step";
const PHASE_START: &str = "fs-constraint:restoration-start";
const PHASE_FINAL_VIOLATIONS: &str = "fs-constraint:restoration-final-violations";
const PHASE_FILTER_SOLVE: &str = "fs-constraint:restoration-filter-solve";
const PHASE_REPAIR_CANDIDATE: &str = "fs-constraint:restoration-repair-candidate";
const PHASE_REPAIR_SAMPLE: &str = "fs-constraint:restoration-repair-sample";

/// Why a planned restoration run stopped mid-flight: the budget
/// authority refused (typed; the caller attaches the retained receipt),
/// or an ordinary constraint-calculus fault surfaced.
enum Stop {
    Refused(BudgetRefusal),
    Fault(ConError),
}

impl From<ConError> for Stop {
    fn from(error: ConError) -> Self {
        Self::Fault(error)
    }
}

type Stopped<T> = Result<T, Stop>;

fn cp(budget: &mut AdmittedBudget<'_>, cx: &Cx<'_>, phase: &'static str) -> Stopped<()> {
    budget.checkpoint(phase, cx).map_err(Stop::Refused)
}

fn charge(budget: &mut AdmittedBudget<'_>, phase: &'static str, units: u64) -> Stopped<()> {
    budget.charge_cost(phase, units).map_err(Stop::Refused)
}

fn violation_contribution(value: f64) -> f64 {
    if value.is_finite() {
        value.max(0.0)
    } else {
        NONFINITE_PENALTY
    }
}

fn checked_total_violation(total: f64, contribution: f64) -> Result<f64, ConError> {
    let next = total + contribution;
    if next.is_finite() {
        Ok(next)
    } else {
        Err(ConError::BadParam {
            what: "elastic total violation",
            value: next,
        })
    }
}

fn validate_domain(problem: &Problem, domain: &DomainBox) -> Result<usize, ConError> {
    if problem.vars().len() != 1 {
        return Err(ConError::InvalidDomain(DomainError::HostVariableCount {
            got: problem.vars().len(),
        }));
    }
    let variable = &problem.vars()[0];
    let Manifold::Rn { dim } = variable.manifold else {
        return Err(ConError::InvalidDomain(DomainError::HostVariableManifold {
            got: variable.manifold,
        }));
    };
    let expected = usize::try_from(dim).map_err(|_| {
        ConError::InvalidDomain(DomainError::PointDimensionUnrepresentable { declared: dim })
    })?;
    if domain.ranges.len() != expected {
        return Err(ConError::InvalidDomain(DomainError::DimensionMismatch {
            expected,
            got: domain.ranges.len(),
        }));
    }
    // Pure bounded arithmetic over caller-provided ranges: NO polling
    // here. Every cancellation must surface through the admitted-budget
    // accountant as a typed refusal, never as a bare teaching error.
    for (axis, &(lo, hi)) in domain.ranges.iter().enumerate() {
        if !lo.is_finite() || !hi.is_finite() {
            return Err(ConError::InvalidDomain(DomainError::InvalidRange {
                axis,
                lo,
                hi,
                reason: DomainRangeError::NonFiniteEndpoint,
            }));
        }
        if lo > hi {
            return Err(ConError::InvalidDomain(DomainError::InvalidRange {
                axis,
                lo,
                hi,
                reason: DomainRangeError::Reversed,
            }));
        }
        if !(hi - lo).is_finite() {
            return Err(ConError::InvalidDomain(DomainError::InvalidRange {
                axis,
                lo,
                hi,
                reason: DomainRangeError::UnrepresentableSpan,
            }));
        }
    }
    Ok(expected)
}

/// Canonical admitted skip mask. Duplicate and out-of-range indices are
/// typed refusals — silent tolerance would let a malformed request change
/// which constraints were relaxed without leaving a stable input identity.
fn canonical_skip_mask(specs_len: usize, skip: &[usize]) -> Result<SkipMask, ConError> {
    let mut skipped = vec![false; specs_len];
    let mut skipped_count = 0u32;
    for &index in skip {
        if index >= specs_len {
            return Err(ConError::BadParam {
                what: "restoration skip index out of range",
                value: index as f64,
            });
        }
        if skipped[index] {
            return Err(ConError::BadParam {
                what: "restoration skip indices must be unique",
                value: index as f64,
            });
        }
        skipped[index] = true;
        skipped_count += 1;
    }
    Ok(SkipMask {
        skipped,
        skipped_count,
    })
}

struct SkipMask {
    /// `true` where the constraint is excluded from the elastic objective.
    skipped: Vec<bool>,
    /// Cardinality; the plan-facing count. Ordering of the caller's skip
    /// list never reaches the plan or the solve.
    skipped_count: u32,
}

impl SkipMask {
    fn keeping_only(specs_len: usize, members: &[usize]) -> Self {
        let mut skipped = vec![true; specs_len];
        let mut kept = 0u32;
        for &member in members {
            if !skipped[member] {
                continue; // defensive: duplicated members stay one mask entry
            }
            skipped[member] = false;
            kept += 1;
        }
        Self {
            skipped_count: u32::try_from(specs_len)
                .unwrap_or(u32::MAX)
                .saturating_sub(kept),
            skipped,
        }
    }
}

/// Internal completed-solution payload (no receipt yet — the receipt
/// belongs to the PUBLIC run that admits the plan).
struct ElasticSolution {
    x: Vec<f64>,
    total_violation: f64,
    violations: Vec<f64>,
    evals: u64,
}

/// Progress visible to the caller even when a run refuses mid-flight.
#[derive(Default)]
struct RunTally {
    starts_completed: u32,
}

/// One total-violation pass over the active set: checkpointed at the
/// deterministic stride, charged as completed work only.
fn run_total_pass(
    problem: &Problem,
    specs: &[ConstraintSpec],
    active_indices: &[usize],
    x: &[f64],
    evals: &mut u64,
    budget: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Stopped<f64> {
    cp(budget, cx, PHASE_TOTAL_PASS)?;
    let mut total = 0.0;
    for (ordinal, &i) in active_indices.iter().enumerate() {
        if ordinal % CANCELLATION_STRIDE == 0 {
            cp(budget, cx, PHASE_TOTAL_PASS)?;
        }
        let gi = scalar_at(problem, specs[i].node, x)?;
        total = checked_total_violation(total, violation_contribution(gi))?;
        *evals += 1;
    }
    charge(
        budget,
        PHASE_TOTAL_PASS,
        u64::try_from(active_indices.len()).unwrap_or(u64::MAX)
            * RESTORATION_UNITS_CONSTRAINT_EVALUATION,
    )?;
    Ok(total)
}

/// One multi-start descent iteration: initial pass plus capped
/// subgradient steps with finite-difference probes, all through the
/// shared accountant. Evaluation ORDER is byte-identical to the retired
/// un-budgeted loop, so ambient-budget runs replay legacy results
/// bit-for-bit.
#[allow(clippy::too_many_lines)] // admission ordering, probes, and step updates stay interleaved on purpose
fn run_descent_start(
    problem: &Problem,
    specs: &[ConstraintSpec],
    active_indices: &[usize],
    domain: &DomainBox,
    seed_x: Vec<f64>,
    limits_steps_per_start: u32,
    evals: &mut u64,
    budget: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Stopped<(Vec<f64>, f64)> {
    let clamp = |x: &mut [f64]| {
        for (value, &(lo, hi)) in x.iter_mut().zip(&domain.ranges) {
            *value = value.clamp(lo, hi);
        }
    };
    let mut x = seed_x;
    let mut v = run_total_pass(problem, specs, active_indices, &x, evals, budget, cx)?;
    let diam: f64 = domain
        .ranges
        .iter()
        .map(|&(lo, hi)| hi - lo)
        .fold(0.0, f64::max);
    for step in 0..limits_steps_per_start {
        cp(budget, cx, PHASE_STEP)?;
        if v <= FEAS_TOL {
            break;
        }
        // FD subgradient of the hinge sum.
        let h = 1e-6 * diam.max(1.0);
        let dimension = x.len();
        let mut gradient = vec![0.0; dimension];
        for (k, gk) in gradient.iter_mut().enumerate() {
            if k % CANCELLATION_STRIDE == 0 {
                cp(budget, cx, PHASE_FD_PROBE)?;
            }
            let mut xp = x.clone();
            xp[k] += h;
            clamp(&mut xp);
            let mut xm = x.clone();
            xm[k] -= h;
            clamp(&mut xm);
            let plus = run_total_pass(problem, specs, active_indices, &xp, evals, budget, cx)?;
            let minus = run_total_pass(problem, specs, active_indices, &xm, evals, budget, cx)?;
            *gk = (plus - minus) / (xp[k] - xm[k]).max(1e-300);
        }
        let norm = gradient.iter().map(|value| value * value).sum::<f64>().sqrt();
        if norm < 1e-14 {
            break;
        }
        let lr = 0.3 * diam / (1.0 + f64::from(step) * 0.05) / norm;
        for (value, gv) in x.iter_mut().zip(&gradient) {
            *value -= lr * gv;
        }
        clamp(&mut x);
        v = run_total_pass(problem, specs, active_indices, &x, evals, budget, cx)?;
    }
    Ok((x, v))
}

/// Full elastic solve under the shared accountant: canonical mask build,
/// midpoint seed, multi-start descent, and the canonical final violation
/// pass whose component sum is the published authority.
#[allow(clippy::too_many_arguments)]
fn run_elastic(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    mask: &SkipMask,
    limits: RestorationWorkLimits,
    tally: &mut RunTally,
    budget: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Stopped<ElasticSolution> {
    let constraint_count = specs.len();
    cp(budget, cx, PHASE_SKIP_MASK)?;
    let mut active_indices = Vec::new();
    for i in 0..constraint_count {
        if i % CANCELLATION_STRIDE == 0 {
            cp(budget, cx, PHASE_SKIP_MASK)?;
        }
        if !mask.skipped[i] {
            active_indices.push(i);
        }
    }
    charge(
        budget,
        PHASE_SKIP_MASK,
        u64::try_from(constraint_count).unwrap_or(u64::MAX) * RESTORATION_UNITS_SKIP_MASK_ENTRY,
    )?;

    let mut evals = 0u64;
    let mut rng = Lcg(0x1001_2026_0707_0001);
    let mut best_x: Vec<f64> = domain
        .ranges
        .iter()
        .map(|&(lo, hi)| f64::midpoint(lo, hi))
        .collect();
    let mut best_v =
        run_total_pass(problem, specs, &active_indices, &best_x, &mut evals, budget, cx)?;
    for start in 0..limits.starts {
        cp(budget, cx, PHASE_START)?;
        let seed = if start == 0 {
            best_x.clone()
        } else {
            domain
                .ranges
                .iter()
                .map(|&(lo, hi)| lo + (hi - lo) * rng.unit())
                .collect()
        };
        let (x, v) = run_descent_start(
            problem,
            specs,
            &active_indices,
            domain,
            seed,
            limits.steps_per_start,
            &mut evals,
            budget,
            cx,
        )?;
        tally.starts_completed = tally.starts_completed.saturating_add(1);
        if v < best_v {
            best_v = v;
            best_x = x;
        }
    }
    let mut violations = Vec::with_capacity(constraint_count);
    for (i, spec) in specs.iter().enumerate() {
        if i % CANCELLATION_STRIDE == 0 {
            cp(budget, cx, PHASE_FINAL_VIOLATIONS)?;
        }
        if mask.skipped[i] {
            violations.push(0.0);
        } else {
            violations.push(violation_contribution(scalar_at(
                problem, spec.node, &best_x,
            )?));
            evals += 1;
        }
    }
    charge(
        budget,
        PHASE_FINAL_VIOLATIONS,
        u64::try_from(active_indices.len()).unwrap_or(u64::MAX)
            * RESTORATION_UNITS_CONSTRAINT_EVALUATION,
    )?;
    // The published component vector is the authority source for the published
    // total. Recompute it in the same canonical order instead of trusting the
    // optimizer-carried `best_v`, which may be stale if final evidence
    // evaluation evolves independently of the search loop.
    let total_violation = violations.iter().try_fold(0.0, |total, &violation| {
        checked_total_violation(total, violation)
    })?;
    Ok(ElasticSolution {
        x: best_x,
        total_violation,
        violations,
        evals,
    })
}

fn memory_authority(cx: &Cx<'_>) -> RestorationMemoryAuthority {
    if cx.lease().is_some() {
        RestorationMemoryAuthority::LeaseAdmitted
    } else {
        RestorationMemoryAuthority::NoLeaseNoClaim
    }
}

fn receipt(
    plan: &RestorationWorkPlan,
    budget: &AdmittedBudget<'_>,
    memory: RestorationMemoryAuthority,
    starts_completed: u32,
) -> RestorationWorkReceipt {
    let consumption: BudgetConsumption = budget.consumption();
    RestorationWorkReceipt {
        plan_identity: plan.identity(),
        schema_version: plan.schema_version,
        consumption: Some(consumption),
        work_units_charged: consumption.cost_charged,
        starts_completed,
        memory,
    }
}

/// Checked plan construction shared by both entry points. Domain
/// admission runs first (raw pre-admission polls), then the canonical
/// skip mask fixes the input identity, then the plan states its cost.
fn prepare(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    skip: &[usize],
    limits: RestorationWorkLimits,
    cx: &Cx<'_>,
) -> Result<(SkipMask, RestorationWorkPlan), RestorationError> {
    let invalid = RestorationError::Invalid;
    let dimension = validate_domain(problem, domain).map_err(invalid)?;
    let mask = canonical_skip_mask(specs.len(), skip).map_err(invalid)?;
    let dimensions = u32::try_from(dimension).map_err(|_| {
        invalid(ConError::BadParam {
            what: "restoration host dimension exceeds the plan's u32 identity space",
            value: dimension as f64,
        })
    })?;
    let constraints_total = u32::try_from(specs.len()).map_err(|_| {
        invalid(ConError::BadParam {
            what: "restoration constraint count exceeds the plan's u32 identity space",
            value: specs.len() as f64,
        })
    })?;
    let plan = RestorationWorkPlan::plan(RestorationWorkShape {
        dimensions,
        constraints_total,
        skipped_count: mask.skipped_count,
        limits,
    })
    .map_err(invalid)?;
    Ok((mask, plan))
}

/// Minimize `Σ max(gᵢ(x), 0)` over the box: multi-start projected
/// subgradient descent (deterministic). Small-fixture machinery — the
/// production restoration solver is a later ASCENT bead.
///
/// Builds the default (historical schedule) [`RestorationWorkPlan`],
/// admits it against the caller's `Cx` budget, and returns the report
/// bound to its [`RestorationWorkReceipt`].
///
/// # Errors
/// [`RestorationError::Invalid`] for malformed domains/skip lists before
/// any budget authority is consumed; evaluation teaching errors carried
/// through; [`RestorationError::Refused`] with the retained receipt when
/// admitted work is stopped by cancellation, deadline, poll, or cost
/// authority.
#[allow(clippy::too_many_lines)] // one elastic solve keeps domain admission, planning, and admission ordered; splitting would interleave them.
pub fn elastic_solve(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    skip: &[usize],
    cx: &Cx<'_>,
) -> Result<ElasticReport, RestorationError> {
    let (mask, plan) = prepare(
        problem,
        specs,
        domain,
        skip,
        RestorationWorkLimits::default(),
        cx,
    )?;
    admit_and_run_elastic(problem, specs, domain, mask, plan, cx)
}

/// Run the elastic solve under a CALLER-SUPPLIED [`RestorationWorkPlan`]
/// (bead frankensim-constraint-restoration-budget-receipts-x5sev). The
/// plan must be schema-current, internally consistent, and describe
/// exactly this host: dimension, constraint count, and active-set
/// cardinality are re-derived here and mismatches refuse before
/// admission.
///
/// # Errors
/// [`RestorationError::Invalid`] for malformed domains/skip lists,
/// inconsistent plans, and evaluator faults — all before or outside the
/// budget contract; [`RestorationError::Refused`] with the retained
/// receipt when admitted work stops at a checkpoint or charge boundary.
pub fn elastic_solve_with_plan(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    skip: &[usize],
    plan: RestorationWorkPlan,
    cx: &Cx<'_>,
) -> Result<ElasticReport, RestorationError> {
    let invalid = RestorationError::Invalid;
    let dimension = validate_domain(problem, domain).map_err(invalid)?;
    let mask = canonical_skip_mask(specs.len(), skip).map_err(invalid)?;
    if plan.schema_version != RESTORATION_WORK_PLAN_SCHEMA_VERSION_ID {
        return Err(invalid(ConError::BadParam {
            what: "unknown restoration work plan schema version",
            value: f64::from(plan.schema_version),
        }));
    }
    plan.verify_consistency().map_err(invalid)?;
    let dimension_matches = u32::try_from(dimension)
        .map_or(false, |admitted| admitted == plan.dimensions);
    if !dimension_matches {
        return Err(invalid(ConError::BadParam {
            what: "restoration work plan dimensions must equal the admitted host dimension",
            value: f64::from(plan.dimensions),
        }));
    }
    if usize::try_from(plan.constraints_total) != Ok(specs.len())
        || plan.active_constraints != specs.len() as u32 - mask.skipped_count
    {
        return Err(invalid(ConError::BadParam {
            what: "restoration work plan constraint counts must equal the evaluated set",
            value: f64::from(plan.constraints_total),
        }));
    }
    admit_and_run_elastic(problem, specs, domain, mask, plan, cx)
}

/// Schema-current marker used by [`elastic_solve_with_plan`]; a plain
/// alias keeps the comparison site honest about what it checks.
const RESTORATION_WORK_PLAN_SCHEMA_VERSION_ID: u32 =
    crate::RESTORATION_WORK_PLAN_SCHEMA_VERSION;

/// Shared admission + execution tail for the two public elastic entry
/// points. `mask` and `plan` arrive pre-validated.
fn admit_and_run_elastic(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    mask: SkipMask,
    plan: RestorationWorkPlan,
    cx: &Cx<'_>,
) -> Result<ElasticReport, RestorationError> {
    let memory = memory_authority(cx);
    let mut budget = match AdmittedBudget::admit_ambient(cx, plan.total_work_units) {
        Ok(budget) => budget,
        Err(refusal) => {
            return Err(RestorationError::Refused {
                refusal,
                receipt: RestorationWorkReceipt::refused_admission(&plan),
            });
        }
    };
    let mut tally = RunTally::default();
    match run_elastic(
        problem,
        specs,
        domain,
        &mask,
        plan.limits,
        &mut tally,
        &mut budget,
        cx,
    ) {
        Ok(solution) => Ok(ElasticReport {
            x: solution.x,
            total_violation: solution.total_violation,
            violations: solution.violations,
            evals: solution.evals,
            work: receipt(&plan, &budget, memory, tally.starts_completed),
        }),
        Err(Stop::Fault(error)) => Err(RestorationError::Invalid(error)),
        Err(Stop::Refused(refusal)) => Err(RestorationError::Refused {
            refusal,
            receipt: receipt(&plan, &budget, memory, tally.starts_completed),
        }),
    }
}

/// One suggested repair.
#[derive(Debug, Clone, PartialEq)]
pub struct RepairAction {
    /// What to do, in words (agent-facing).
    pub description: String,
    /// Structured form.
    pub kind: RepairKind,
    /// Estimated probability the repaired space is feasible
    /// (Monte-Carlo over the domain; calibrated in the battery).
    pub feasibility_estimate: f64,
}

/// Structured repair kinds.
#[derive(Debug, Clone, PartialEq)]
pub enum RepairKind {
    /// Relax constraint `index` by adding `slack` to its bound
    /// (`g ≤ 0` becomes `g ≤ slack`).
    RelaxBound {
        /// Which constraint.
        index: usize,
        /// Added slack.
        slack: f64,
    },
    /// Drop a soft constraint entirely.
    DropSoft {
        /// Which constraint.
        index: usize,
    },
}

/// The full diagnosis payload (the agent-facing artifact).
#[derive(Debug, Clone)]
pub struct Diagnosis {
    /// Whether the space is feasible as stated.
    pub feasible: bool,
    /// A feasible point when one exists.
    pub witness: Option<Vec<f64>>,
    /// MINIMAL unsat core (constraint indices), empty when feasible.
    pub core: Vec<usize>,
    /// Ranked repairs (best first), empty when feasible.
    pub repairs: Vec<RepairAction>,
    /// Elastic-solve evidence.
    pub elastic: ElasticReport,
    /// Retained resource contract of the WHOLE diagnosis run: the base
    /// solve, support verification, deletion filtering, verification,
    /// and repair sampling all charged this one admitted budget.
    pub work: RestorationWorkReceipt,
}

impl Diagnosis {
    fn invalid_reason(&self, specs: &[ConstraintSpec]) -> Option<&'static str> {
        if !self.elastic.total_violation.is_finite() {
            return Some("nonfinite-total-violation");
        }
        if self.elastic.total_violation < 0.0 {
            return Some("negative-total-violation");
        }
        if self.elastic.x.iter().any(|value| !value.is_finite()) {
            return Some("nonfinite-elastic-point");
        }
        if self.elastic.violations.len() != specs.len() {
            return Some("component-violation-count-mismatch");
        }
        if self
            .elastic
            .violations
            .iter()
            .any(|value| !value.is_finite())
        {
            return Some("nonfinite-component-violation");
        }
        if self.elastic.violations.iter().any(|&value| value < 0.0) {
            return Some("negative-component-violation");
        }
        let component_total = self.elastic.violations.iter().sum::<f64>();
        if !component_total.is_finite() {
            return Some("nonfinite-component-violation-total");
        }
        #[allow(clippy::float_cmp)]
        // The recomputed component sum must equal the stored total bitwise; that equality IS the determinism invariant.
        if component_total != self.elastic.total_violation {
            return Some("total-component-violation-mismatch");
        }
        if self.core.iter().any(|&index| index >= specs.len()) {
            return Some("unknown-core-constraint");
        }
        if self.core.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Some("noncanonical-core-order");
        }
        if self
            .witness
            .as_ref()
            .is_some_and(|point| point.iter().any(|value| !value.is_finite()))
        {
            return Some("nonfinite-witness");
        }
        for repair in &self.repairs {
            if !repair.feasibility_estimate.is_finite()
                || !(0.0..=1.0).contains(&repair.feasibility_estimate)
            {
                return Some("invalid-feasibility-estimate");
            }
            match &repair.kind {
                RepairKind::RelaxBound { index, slack } => {
                    if *index >= specs.len() {
                        return Some("unknown-repair-constraint");
                    }
                    if !slack.is_finite() || *slack < 0.0 {
                        return Some("invalid-repair-slack");
                    }
                }
                RepairKind::DropSoft { index } => {
                    if *index >= specs.len() {
                        return Some("unknown-repair-constraint");
                    }
                }
            }
        }
        if self.feasible {
            if self.elastic.total_violation > FEAS_TOL {
                return Some("feasible-claim-has-positive-violation");
            }
            if self.witness.is_none() {
                return Some("feasible-claim-missing-witness");
            }
            if self.witness.as_deref() != Some(self.elastic.x.as_slice()) {
                return Some("witness-does-not-match-elastic-point");
            }
            if self
                .elastic
                .violations
                .iter()
                .any(|&violation| violation > FEAS_TOL)
            {
                return Some("feasible-claim-has-component-violation");
            }
            if !self.core.is_empty() || !self.repairs.is_empty() {
                return Some("feasible-claim-has-conflict-evidence");
            }
        } else {
            if self.elastic.total_violation <= FEAS_TOL {
                return Some("infeasible-claim-lacks-positive-violation");
            }
            if self.witness.is_some() {
                return Some("infeasible-claim-has-witness");
            }
            if self.core.is_empty() {
                return Some("infeasible-claim-missing-core");
            }
        }
        None
    }

    /// Canonical JSON payload for the ledger/session surface. Dynamic text is
    /// escaped. A publicly forged inconsistent or non-finite diagnosis emits a
    /// deterministic invalid/no-claim object; it never retains `feasible:true`
    /// while silently replacing required evidence with `null`.
    #[must_use]
    pub fn to_json(&self, specs: &[ConstraintSpec]) -> String {
        use std::fmt::Write as _;

        if let Some(reason) = self.invalid_reason(specs) {
            let mut invalid = "{\"valid\":false,\"reason\":".to_string();
            push_json_string(&mut invalid, reason);
            invalid.push_str(
                ",\"feasible\":false,\"total_violation\":null,\"core\":[],\"repairs\":[]}",
            );
            return invalid;
        }

        let mut s = format!("{{\"feasible\":{},\"total_violation\":", self.feasible);
        let _ = write!(s, "{:.3e}", self.elastic.total_violation);
        s.push_str(",\"core\":[");
        for (k, &i) in self.core.iter().enumerate() {
            if k > 0 {
                s.push(',');
            }
            push_json_string(&mut s, &specs[i].name);
        }
        s.push_str("],\"repairs\":[");
        for (k, r) in self.repairs.iter().enumerate() {
            if k > 0 {
                s.push(',');
            }
            s.push_str("{\"action\":");
            push_json_string(&mut s, &r.description);
            s.push_str(",\"est_feasible\":");
            let _ = write!(s, "{:.2}", r.feasibility_estimate);
            s.push('}');
        }
        s.push_str("]}");
        s
    }
}

/// Monte-Carlo feasible-volume fraction with constraint `relax[i]`
/// slack applied (the repair feasibility estimator; deterministic).
/// Every COMPLETED sample is charged its performed evaluations against
/// the shared accountant before the next sample starts.
fn feasible_fraction(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    relax: &[(usize, f64)],
    drop: Option<usize>,
    samples: u32,
    budget: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Stopped<f64> {
    let mut rng = Lcg(0x1001_2026_0707_0002);
    let mut hits = 0u32;
    for sample in 0..samples {
        if sample % u32::try_from(CANCELLATION_STRIDE).expect("small stride") == 0 {
            cp(budget, cx, PHASE_REPAIR_SAMPLE)?;
        }
        let x: Vec<f64> = domain
            .ranges
            .iter()
            .map(|&(lo, hi)| lo + (hi - lo) * rng.unit())
            .collect();
        let mut ok = true;
        let mut performed = 0u64;
        for (i, spec) in specs.iter().enumerate() {
            if i % CANCELLATION_STRIDE == 0 {
                cp(budget, cx, PHASE_REPAIR_SAMPLE)?;
            }
            if Some(i) == drop {
                continue;
            }
            let slack = relax.iter().find(|(j, _)| *j == i).map_or(0.0, |(_, s)| *s);
            // A non-finite constraint value is undefined here, hence NOT feasible
            // — `NaN > slack` is false, which would otherwise count the sample as
            // feasible and inflate the feasibility estimate.
            let gi = scalar_at(problem, spec.node, &x)?;
            performed += 1;
            if !gi.is_finite() || gi > slack {
                ok = false;
                break;
            }
        }
        charge(
            budget,
            PHASE_REPAIR_SAMPLE,
            performed * RESTORATION_UNITS_CONSTRAINT_EVALUATION,
        )?;
        if ok {
            hits += 1;
        }
    }
    Ok(f64::from(hits) / f64::from(samples))
}

/// Diagnose a constraint set over a domain: feasibility, MINIMAL unsat
/// core (deletion-filtered), and ranked repairs with feasibility
/// estimates.
///
/// Builds one aggregate [`RestorationWorkPlan`] covering the base solve,
/// the structural worst case of deletion filtering (`2N + 2` subset
/// solves) and repair estimation (`3N` Monte-Carlo estimates), admits it
/// ONCE against the caller's `Cx` budget, and binds the retained
/// [`RestorationWorkReceipt`] to the returned [`Diagnosis`] and to its
/// embedded [`ElasticReport`].
///
/// # Errors
/// [`RestorationError::Invalid`] for malformed domains and evaluation
/// teaching errors; [`RestorationError::Refused`] with the retained
/// receipt when the admitted run is stopped by cancellation, deadline,
/// poll, or cost authority at any phase.
#[allow(clippy::too_many_lines)] // the deletion filter's invariant chain stays linear and commented
pub fn diagnose_infeasibility(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    cx: &Cx<'_>,
) -> Result<Diagnosis, RestorationError> {
    let invalid = RestorationError::Invalid;
    let (mask, plan) = prepare(
        problem,
        specs,
        domain,
        &[],
        RestorationWorkLimits::default(),
        cx,
    )?;
    let memory = memory_authority(cx);
    let mut budget = match AdmittedBudget::admit_ambient(cx, plan.total_work_units) {
        Ok(budget) => budget,
        Err(refusal) => {
            return Err(RestorationError::Refused {
                refusal,
                receipt: RestorationWorkReceipt::refused_admission(&plan),
            });
        }
    };
    let stopped =
        |budget: &AdmittedBudget<'_>, refusal: BudgetRefusal| RestorationError::Refused {
            refusal,
            receipt: receipt(&plan, budget, memory, 0),
        };

    let mut tally = RunTally::default();
    let base = match run_elastic(
        problem,
        specs,
        domain,
        &mask,
        plan.limits,
        &mut tally,
        &mut budget,
        cx,
    ) {
        Ok(solution) => solution,
        Err(Stop::Fault(error)) => return Err(invalid(error)),
        Err(Stop::Refused(refusal)) => return Err(stopped(&budget, refusal)),
    };
    if base.total_violation <= FEAS_TOL {
        let work = receipt(&plan, &budget, memory, tally.starts_completed);
        let elastic = ElasticReport {
            x: base.x.clone(),
            total_violation: base.total_violation,
            violations: base.violations,
            evals: base.evals,
            work,
        };
        return Ok(Diagnosis {
            feasible: true,
            witness: Some(base.x),
            core: Vec::new(),
            repairs: Vec::new(),
            elastic,
            work,
        });
    }

    // Candidate core: the elastic support (violated at the optimum). A
    // support identifies the sum-optimum's active trade-off, but need not be
    // jointly infeasible by itself. Verify it before deletion filtering and
    // deterministically expand to the full, already-proven infeasible set when
    // the support is feasible.
    let mut core: Vec<usize> = base
        .violations
        .iter()
        .enumerate()
        .filter(|&(_, &v)| v > FEAS_TOL)
        .map(|(i, _)| i)
        .collect();

    let verify_subset =
        |core_members: &[usize], budget: &mut AdmittedBudget<'_>| -> Stopped<ElasticSolution> {
            cp(budget, cx, PHASE_FILTER_SOLVE)?;
            let subset_mask = SkipMask::keeping_only(specs.len(), core_members);
            let mut subset_tally = RunTally::default();
            run_elastic(
                problem,
                specs,
                domain,
                &subset_mask,
                plan.limits,
                &mut subset_tally,
                budget,
                cx,
            )
        };

    let support = match verify_subset(&core, &mut budget) {
        Ok(solution) => solution,
        Err(Stop::Fault(error)) => return Err(invalid(error)),
        Err(Stop::Refused(refusal)) => return Err(stopped(&budget, refusal)),
    };
    if support.total_violation <= FEAS_TOL {
        core = (0..specs.len()).collect();
    }

    // Deletion filter for MINIMALITY. The current core is jointly infeasible
    // on entry. A removal is installed only when the resulting subset is also
    // jointly infeasible, so that invariant is preserved at every step.
    let mut k = 0;
    while k < core.len() {
        cp(&mut budget, cx, PHASE_FILTER_SOLVE).map_err(|stop| match stop {
            Stop::Refused(refusal) => stopped(&budget, refusal),
            Stop::Fault(error) => invalid(error),
        })?;
        let mut without_members = core.clone();
        without_members.remove(k);
        let without = match verify_subset(&without_members, &mut budget) {
            Ok(solution) => solution,
            Err(Stop::Fault(error)) => return Err(invalid(error)),
            Err(Stop::Refused(refusal)) => return Err(stopped(&budget, refusal)),
        };
        if without.total_violation <= FEAS_TOL {
            k += 1; // necessary: dropping it restores feasibility
        } else {
            core = without_members; // redundant: still infeasible without it
        }
    }
    let verified_core = match verify_subset(&core, &mut budget) {
        Ok(solution) => solution,
        Err(Stop::Fault(error)) => return Err(invalid(error)),
        Err(Stop::Refused(refusal)) => return Err(stopped(&budget, refusal)),
    };
    assert!(
        verified_core.total_violation > FEAS_TOL,
        "deletion filtering must not publish a jointly feasible unsat core"
    );

    // Repairs: relax each core member by graded slacks, or drop it if
    // soft; estimate feasibility by Monte-Carlo volume; rank.
    let mut repairs = Vec::new();
    for &i in &core {
        cp(&mut budget, cx, PHASE_REPAIR_CANDIDATE)
            .map_err(|stop| match stop {
                Stop::Refused(refusal) => stopped(&budget, refusal),
                Stop::Fault(error) => invalid(error),
            })?;
        let scale = base.violations[i].max(FEAS_TOL);
        for factor in [1.1, 1.5] {
            let slack = scale * factor;
            if !slack.is_finite() {
                return Err(invalid(ConError::BadParam {
                    what: "repair slack",
                    value: slack,
                }));
            }
            let est = match feasible_fraction(
                problem,
                specs,
                domain,
                &[(i, slack)],
                None,
                plan.limits.feasibility_samples,
                &mut budget,
                cx,
            ) {
                Ok(est) => est,
                Err(Stop::Fault(error)) => return Err(invalid(error)),
                Err(Stop::Refused(refusal)) => return Err(stopped(&budget, refusal)),
            };
            repairs.push(RepairAction {
                description: format!("relax `{}` by {slack:.3} (g <= {slack:.3})", specs[i].name),
                kind: RepairKind::RelaxBound { index: i, slack },
                feasibility_estimate: est,
            });
        }
        if matches!(specs[i].kind, crate::ConstraintKind::Soft(_)) {
            let est = match feasible_fraction(
                problem,
                specs,
                domain,
                &[],
                Some(i),
                plan.limits.feasibility_samples,
                &mut budget,
                cx,
            ) {
                Ok(est) => est,
                Err(Stop::Fault(error)) => return Err(invalid(error)),
                Err(Stop::Refused(refusal)) => return Err(stopped(&budget, refusal)),
            };
            repairs.push(RepairAction {
                description: format!("drop soft constraint `{}`", specs[i].name),
                kind: RepairKind::DropSoft { index: i },
                feasibility_estimate: est,
            });
        }
    }
    repairs.sort_by(|a, b| {
        b.feasibility_estimate
            .partial_cmp(&a.feasibility_estimate)
            .expect("estimates are finite")
            .then_with(|| a.description.cmp(&b.description))
    });

    let work = receipt(&plan, &budget, memory, tally.starts_completed);
    let elastic = ElasticReport {
        x: base.x,
        total_violation: base.total_violation,
        violations: base.violations,
        evals: base.evals,
        work,
    };
    Ok(Diagnosis {
        feasible: false,
        witness: None,
        core,
        repairs,
        elastic,
        work,
    })
}
