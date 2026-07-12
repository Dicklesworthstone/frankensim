//! VALUE-OF-INFORMATION QUERIES (addendum Proposal C, bead knh1.6;
//! \[F\] — behind the `voi-queries` feature): THE IGNORANCE MARKET, v0
//! as a RANKED LIST. Across everything the ledger is uncertain about,
//! where does one dollar of evidence most change the downstream
//! decision? Decision sensitivity (does the decision FLIP inside the
//! node's interval?) is computed from CACHED surrogate sweeps, crossed
//! with a PRICED PROBE MENU that unifies computational and physical
//! experiments, ranked by SAMPLED FLIP-FRACTION REDUCTION PER DOLLAR.
//!
//! MYOPIC one-step VoI ONLY (the proposal's own discipline: full
//! sequential VoI is intractable and myopic captures most of the
//! value). The output surfaces as (i) the hint on query results — the
//! upgrade the fs-ir anytime module's CONTRACT reserved for Proposal C
//! — and (ii) the scheduler for discrepancy probes.
//!
//! THE KILL CRITERION AS CODE: [`VoiScheduler`] owns one append-only
//! chronological pairwise e-process, the remaining spend budget, and consumed
//! decision snapshots. Scheduling consults the CURRENT audit verdict under
//! exclusive mutation; stale reports carry no spending authority.

use std::collections::BTreeSet;

use fs_blake3::{ContentHash, hash_domain};
use fs_eproc::{LossSpan, PairwiseRace};

/// Maximum uncertainty nodes in one myopic VoI request.
pub const MAX_VOI_NODES: usize = 256;
/// Maximum probe menu entries (and scheduled ranked entries).
pub const MAX_VOI_PROBES: usize = 1024;
/// Maximum visible-ASCII byte length of node, probe, target, and audit names.
pub const MAX_VOI_NAME_BYTES: usize = 128;
/// Maximum interval-sweep grid size.
pub const MAX_VOI_GRID: usize = 1024;
/// Maximum surrogate evaluations admitted by one public VoI call.
pub const MAX_VOI_EVALUATIONS: usize = 4096;
/// Maximum matched-cost observations admitted by one prospective audit.
pub const MAX_VOI_AUDIT_RECORDS: usize = 4096;
/// Maximum distinct ranking snapshots retained by one live scheduler.
pub const MAX_VOI_SCHEDULED_CONTEXTS: usize = 4096;
/// Fixed anytime-valid false-activation level for VoI scheduling authority.
pub const VOI_AUDIT_ALPHA: f64 = 0.05;

const RANKED_MENU_SOURCE_DOMAIN: &str = "frankensim.fs-plan.voi-ranked-source.v1";
const RANKED_MENU_CONTEXT_DOMAIN: &str = "frankensim.fs-plan.voi-ranked-menu.v2";
const AUDIT_CONTEXT_DOMAIN: &str = "frankensim.fs-plan.voi-audit.v2";

/// Why a VoI query, audit, or schedule was refused.
#[derive(Debug, Clone, PartialEq)]
pub enum VoiError {
    /// A bounded collection falls outside its admitted size range.
    SizeLimit {
        /// Collection being validated.
        collection: &'static str,
        /// Supplied element count.
        count: usize,
        /// Inclusive lower bound.
        min: usize,
        /// Inclusive upper bound.
        max: usize,
    },
    /// The surrogate's declared arity differs from the node vector.
    ArityMismatch {
        /// Declared surrogate arity.
        arity: usize,
        /// Supplied node count.
        node_count: usize,
    },
    /// A node/probe/target/audit identity is not bounded visible ASCII.
    InvalidName {
        /// Name category.
        kind: &'static str,
        /// Position in its collection.
        index: usize,
        /// Supplied UTF-8 byte length.
        bytes: usize,
        /// Inclusive byte limit.
        max_bytes: usize,
    },
    /// A supposedly unique name occurs more than once.
    DuplicateName {
        /// Name category.
        kind: &'static str,
        /// Duplicate value (already bounded by [`MAX_VOI_NAME_BYTES`]).
        name: String,
    },
    /// An interval is nonfinite, unordered, too wide for finite
    /// arithmetic, or excludes its nominal value.
    InvalidInterval {
        /// Node name.
        node: String,
        /// Lower endpoint.
        lo: f64,
        /// Nominal value.
        nominal: f64,
        /// Upper endpoint.
        hi: f64,
    },
    /// A surrogate returned a nonfinite decision margin.
    NonFiniteMargin {
        /// Returned margin.
        value: f64,
    },
    /// A sensitivity request names a missing node.
    NodeIndexOutOfRange {
        /// Supplied node index.
        node_idx: usize,
        /// Supplied node count.
        node_count: usize,
    },
    /// The sweep grid is zero or exceeds the declared cap.
    InvalidGrid {
        /// Supplied grid size.
        grid: usize,
        /// Inclusive upper bound.
        max: usize,
    },
    /// A request would exceed the surrogate-evaluation budget.
    EvaluationLimitExceeded {
        /// Required evaluations.
        requested: usize,
        /// Inclusive limit.
        max: usize,
    },
    /// A probe has an invalid numeric field.
    InvalidProbeValue {
        /// Probe name.
        probe: String,
        /// Invalid field (`cost` or `shrink`).
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
    /// A probe target resolves to zero or multiple nodes.
    TargetResolution {
        /// Probe name.
        probe: String,
        /// Requested target.
        target: String,
        /// Number of matching nodes.
        matches: usize,
    },
    /// The scheduling budget is nonfinite or negative.
    InvalidBudget {
        /// Supplied budget.
        budget: f64,
    },
    /// An audit observation has a malformed finite matched-cost pair.
    InvalidAuditCost {
        /// Observation identity.
        observation: String,
        /// Recommended-purchase cost.
        recommended_cost: f64,
        /// Alternative-purchase cost.
        alternative_cost: f64,
    },
    /// An audit compares a purchase with itself.
    InvalidAuditPair {
        /// Observation identity.
        observation: String,
    },
    /// An audit repeats an observation identity and could double-count evidence.
    DuplicateAuditObservation {
        /// Repeated observation identity.
        observation: String,
    },
    /// A ranked menu belongs to a different scheduling policy.
    PolicyScopeMismatch {
        /// Policy scope owned by the scheduler.
        expected: String,
        /// Policy scope bound into the ranked menu.
        actual: String,
    },
    /// One decision snapshot was already evaluated by this scheduler.
    RankingSnapshotAlreadyConsumed {
        /// Duplicate caller-declared decision/ledger snapshot.
        snapshot_id: String,
    },
    /// Scheduling was requested before the live audit crossed its threshold.
    MissingSchedulingAuthority,
    /// Finite inputs could not produce a finite, monotone result.
    ArithmeticRefusal {
        /// Operation that failed.
        operation: &'static str,
        /// Bounded subject name.
        subject: String,
    },
}

impl core::fmt::Display for VoiError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SizeLimit {
                collection,
                count,
                min,
                max,
            } => write!(
                f,
                "{collection} has {count} entries; expected an inclusive range of {min}..={max}"
            ),
            Self::ArityMismatch { arity, node_count } => write!(
                f,
                "surrogate declares arity {arity}, but {node_count} uncertainty node(s) were supplied"
            ),
            Self::InvalidName {
                kind,
                index,
                bytes,
                max_bytes,
            } => write!(
                f,
                "{kind} name at index {index} is not nonempty visible ASCII or is {bytes} bytes long (limit {max_bytes})"
            ),
            Self::DuplicateName { kind, name } => {
                write!(f, "duplicate {kind} name {name:?}")
            }
            Self::InvalidInterval {
                node,
                lo,
                nominal,
                hi,
            } => write!(
                f,
                "node {node:?} has invalid interval [{lo:?}, {hi:?}] with nominal {nominal:?}"
            ),
            Self::NonFiniteMargin { value } => {
                write!(f, "surrogate returned nonfinite margin {value:?}")
            }
            Self::NodeIndexOutOfRange {
                node_idx,
                node_count,
            } => write!(
                f,
                "node index {node_idx} is out of range for {node_count} node(s)"
            ),
            Self::InvalidGrid { grid, max } => {
                write!(f, "sweep grid {grid} is outside 1..={max}")
            }
            Self::EvaluationLimitExceeded { requested, max } => write!(
                f,
                "VoI request needs {requested} surrogate evaluations; the limit is {max}"
            ),
            Self::InvalidProbeValue {
                probe,
                field,
                value,
            } => write!(f, "probe {probe:?} has invalid {field} {value:?}"),
            Self::TargetResolution {
                probe,
                target,
                matches,
            } => write!(
                f,
                "probe {probe:?} target {target:?} resolves to {matches} uncertainty node(s), expected exactly one"
            ),
            Self::InvalidBudget { budget } => {
                write!(
                    f,
                    "probe budget must be finite and non-negative, got {budget:?}"
                )
            }
            Self::InvalidAuditCost {
                observation,
                recommended_cost,
                alternative_cost,
            } => write!(
                f,
                "audit observation {observation:?} requires equal finite positive matched costs, got {recommended_cost:?} and {alternative_cost:?}"
            ),
            Self::InvalidAuditPair { observation } => write!(
                f,
                "audit observation {observation:?} compares a purchase with itself"
            ),
            Self::DuplicateAuditObservation { observation } => write!(
                f,
                "audit observation {observation:?} appears more than once"
            ),
            Self::PolicyScopeMismatch { expected, actual } => write!(
                f,
                "ranked menu policy scope {actual:?} does not match scheduler scope {expected:?}"
            ),
            Self::RankingSnapshotAlreadyConsumed { snapshot_id } => write!(
                f,
                "decision snapshot {snapshot_id:?} was already consumed by this scheduler"
            ),
            Self::MissingSchedulingAuthority => write!(
                f,
                "the live VoI audit has not crossed the anytime-valid scheduling threshold"
            ),
            Self::ArithmeticRefusal { operation, subject } => {
                write!(
                    f,
                    "{operation} for {subject:?} did not remain finite and monotone"
                )
            }
        }
    }
}

impl std::error::Error for VoiError {}

/// One uncertainty node touching a live decision: a named quantity the
/// ledger only knows to an interval.
#[derive(Debug, Clone, PartialEq)]
pub struct UncertaintyNode {
    /// Ledger name.
    pub name: String,
    /// Current uncertainty interval.
    pub lo: f64,
    /// Upper end.
    pub hi: f64,
    /// Nominal (decision-time) value.
    pub nominal: f64,
}

/// A live decision over the uncertain quantities: v0 is a threshold
/// verdict on a cheap cached surrogate `margin(values) > 0`.
pub struct LiveDecision<'a> {
    /// The cached surrogate margin (cheap by Proposals 9/2).
    pub margin: &'a dyn Fn(&[f64]) -> f64,
    /// Node count.
    pub arity: usize,
}

fn validate_size(
    collection: &'static str,
    count: usize,
    min: usize,
    max: usize,
) -> Result<(), VoiError> {
    if (min..=max).contains(&count) {
        Ok(())
    } else {
        Err(VoiError::SizeLimit {
            collection,
            count,
            min,
            max,
        })
    }
}

fn validate_name(kind: &'static str, index: usize, name: &str) -> Result<(), VoiError> {
    if name.is_empty()
        || name.len() > MAX_VOI_NAME_BYTES
        || !name.bytes().all(|byte| byte.is_ascii_graphic())
    {
        Err(VoiError::InvalidName {
            kind,
            index,
            bytes: name.len(),
            max_bytes: MAX_VOI_NAME_BYTES,
        })
    } else {
        Ok(())
    }
}

fn validate_nodes(decision: &LiveDecision<'_>, nodes: &[UncertaintyNode]) -> Result<(), VoiError> {
    validate_size("uncertainty nodes", nodes.len(), 1, MAX_VOI_NODES)?;
    if decision.arity != nodes.len() {
        return Err(VoiError::ArityMismatch {
            arity: decision.arity,
            node_count: nodes.len(),
        });
    }
    let mut names = BTreeSet::new();
    for (index, node) in nodes.iter().enumerate() {
        validate_name("uncertainty node", index, &node.name)?;
        if !names.insert(node.name.as_str()) {
            return Err(VoiError::DuplicateName {
                kind: "uncertainty node",
                name: node.name.clone(),
            });
        }
        let width = node.hi - node.lo;
        if !node.lo.is_finite()
            || !node.hi.is_finite()
            || !node.nominal.is_finite()
            || node.lo > node.hi
            || node.nominal < node.lo
            || node.nominal > node.hi
            || !width.is_finite()
        {
            return Err(VoiError::InvalidInterval {
                node: node.name.clone(),
                lo: node.lo,
                nominal: node.nominal,
                hi: node.hi,
            });
        }
    }
    Ok(())
}

fn validate_grid(grid: usize) -> Result<(), VoiError> {
    if (1..=MAX_VOI_GRID).contains(&grid) {
        Ok(())
    } else {
        Err(VoiError::InvalidGrid {
            grid,
            max: MAX_VOI_GRID,
        })
    }
}

fn validate_evaluations(requested: usize) -> Result<(), VoiError> {
    if requested <= MAX_VOI_EVALUATIONS {
        Ok(())
    } else {
        Err(VoiError::EvaluationLimitExceeded {
            requested,
            max: MAX_VOI_EVALUATIONS,
        })
    }
}

fn evaluate_margin(decision: &LiveDecision<'_>, values: &[f64]) -> Result<f64, VoiError> {
    let margin = (decision.margin)(values);
    if margin.is_finite() {
        Ok(margin)
    } else {
        Err(VoiError::NonFiniteMargin { value: margin })
    }
}

fn nominal_values(nodes: &[UncertaintyNode]) -> Vec<f64> {
    nodes.iter().map(|node| node.nominal).collect()
}

fn node_at(nodes: &[UncertaintyNode], node_idx: usize) -> Result<&UncertaintyNode, VoiError> {
    nodes.get(node_idx).ok_or(VoiError::NodeIndexOutOfRange {
        node_idx,
        node_count: nodes.len(),
    })
}

fn nominal_verdict_validated(
    decision: &LiveDecision<'_>,
    nodes: &[UncertaintyNode],
) -> Result<bool, VoiError> {
    Ok(evaluate_margin(decision, &nominal_values(nodes))? > 0.0)
}

fn flip_probability_validated(
    decision: &LiveDecision<'_>,
    nodes: &[UncertaintyNode],
    base: bool,
    node_idx: usize,
    lo: f64,
    hi: f64,
    grid: usize,
) -> Result<f64, VoiError> {
    let mut values = nominal_values(nodes);
    let node = node_at(nodes, node_idx)?;
    let width = hi - lo;
    let mut flips = 0usize;
    for k in 0..grid {
        #[allow(clippy::cast_precision_loss)]
        let t = (k as f64 + 0.5) / grid as f64;
        let sample = lo + t * width;
        if !sample.is_finite() {
            return Err(VoiError::ArithmeticRefusal {
                operation: "interval sweep",
                subject: node.name.clone(),
            });
        }
        *values
            .get_mut(node_idx)
            .ok_or(VoiError::NodeIndexOutOfRange {
                node_idx,
                node_count: nodes.len(),
            })? = sample;
        if (evaluate_margin(decision, &values)? > 0.0) != base {
            flips += 1;
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let probability = flips as f64 / grid as f64;
    Ok(probability)
}

impl LiveDecision<'_> {
    /// The nominal verdict.
    ///
    /// # Errors
    /// [`VoiError`] when node/arity/interval invariants fail or the
    /// cached surrogate returns a nonfinite margin.
    pub fn nominal_verdict(&self, nodes: &[UncertaintyNode]) -> Result<bool, VoiError> {
        validate_nodes(self, nodes)?;
        validate_evaluations(1)?;
        nominal_verdict_validated(self, nodes)
    }

    /// DECISION SENSITIVITY of one node: sweep the node's interval on
    /// the cached surrogate (others at nominal, `grid` points) and
    /// return the fraction of MIDPOINT GRID SAMPLES where the verdict differs
    /// from nominal. This is a myopic estimate under the uniform interval
    /// measure, not a certified probability.
    ///
    /// # Errors
    /// [`VoiError`] when the request is malformed, exceeds the declared
    /// sweep/evaluation limits, or the surrogate returns a nonfinite
    /// margin.
    pub fn flip_probability(
        &self,
        nodes: &[UncertaintyNode],
        node_idx: usize,
        grid: usize,
    ) -> Result<f64, VoiError> {
        validate_nodes(self, nodes)?;
        let node = node_at(nodes, node_idx)?;
        validate_grid(grid)?;
        let evaluations = grid
            .checked_add(1)
            .ok_or_else(|| VoiError::ArithmeticRefusal {
                operation: "sweep evaluation count",
                subject: node.name.clone(),
            })?;
        validate_evaluations(evaluations)?;
        let base = nominal_verdict_validated(self, nodes)?;
        flip_probability_validated(self, nodes, base, node_idx, node.lo, node.hi, grid)
    }
}

/// The kind of evidence purchase — the menu UNIFIES computational and
/// physical experiments (the epistemic-engine identity made concrete).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProbeKind {
    /// Climb a fidelity rung / refine / add solver accuracy.
    Computational,
    /// Wind-tunnel anchor, CT scan, strain gauge — reality as evidence.
    Physical,
}

/// One priced probe: buying it SHRINKS a node's interval around its
/// nominal by `shrink` (0 < shrink < 1).
#[derive(Debug, Clone, PartialEq)]
pub struct Probe {
    /// Menu name ("climb-to-rung-96", "wind-tunnel-anchor", ...).
    pub name: String,
    /// Which node it tightens.
    pub target: String,
    /// Price in dollars.
    pub cost: f64,
    /// Post-probe interval width as a fraction of the current width.
    pub shrink: f64,
    /// Computational or physical.
    pub kind: ProbeKind,
}

/// One ranked purchase: the myopic VoI score.
#[derive(Debug, Clone, PartialEq)]
pub struct RankedPurchase {
    /// The probe.
    probe: Probe,
    /// Grid-sampled flip fraction before the purchase.
    flip_before: f64,
    /// Grid-sampled flip fraction after the declared contraction.
    flip_after: f64,
    /// THE SCORE: sampled flip-fraction reduction per dollar.
    score: f64,
}

impl RankedPurchase {
    /// The validated probe purchase.
    #[must_use]
    pub fn probe(&self) -> &Probe {
        &self.probe
    }

    /// Grid-sampled flip fraction before the purchase.
    #[must_use]
    pub fn flip_before(&self) -> f64 {
        self.flip_before
    }

    /// Grid-sampled flip fraction after the declared contraction.
    #[must_use]
    pub fn flip_after(&self) -> f64 {
        self.flip_after
    }

    /// Grid-sampled flip-fraction reduction per dollar.
    #[must_use]
    pub fn score(&self) -> f64 {
        self.score
    }
}

/// A complete, canonical ranking for one validated supplied
/// uncertainty/menu/grid snapshot. Rows and context are private so safe callers
/// cannot omit, splice, or reorder scheduling authority after ranking.
#[derive(Debug, PartialEq)]
pub struct RankedMenu {
    rows: Vec<RankedPurchase>,
    source_context_id: ContentHash,
    context_id: ContentHash,
    policy_scope: String,
    snapshot_id: String,
    grid: usize,
}

impl RankedMenu {
    /// Number of ranked purchases in the complete menu.
    #[must_use]
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// A ranked menu produced by [`rank_purchases`] is never empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Inspect one canonical row without exposing mutable membership.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<&RankedPurchase> {
        self.rows.get(index)
    }

    /// Inspect the highest-ranked purchase.
    #[must_use]
    pub fn top(&self) -> Option<&RankedPurchase> {
        self.rows.first()
    }

    /// Iterate over canonical rows for reporting only.
    #[must_use]
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &RankedPurchase> {
        self.rows.iter()
    }

    /// Midpoint grid used for every sampled flip estimate in this menu.
    #[must_use]
    pub fn grid(&self) -> usize {
        self.grid
    }

    /// BLAKE3 identity of the validated node/menu/grid snapshot.
    ///
    /// This binds supplied content but does not identify callback code, prove
    /// catalog completeness, or prove that the snapshot remains current;
    /// callers must compare it with their ledger/session snapshot before use.
    #[must_use]
    pub fn context_id(&self) -> ContentHash {
        self.context_id
    }

    /// Identity of policy, snapshot, nodes, source menu, and grid before
    /// evaluating the decision surrogate.
    #[must_use]
    pub fn source_context_id(&self) -> ContentHash {
        self.source_context_id
    }

    /// Caller-declared policy/version scope bound into this ranking.
    #[must_use]
    pub fn policy_scope(&self) -> &str {
        &self.policy_scope
    }

    /// Caller-declared decision/ledger snapshot bound into this ranking.
    #[must_use]
    pub fn snapshot_id(&self) -> &str {
        &self.snapshot_id
    }
}

/// Structured, grid-qualified query hint. Its private optional purchase keeps
/// the no-sampled-change state distinct from an authoritative zero claim.
#[derive(Debug, Clone, PartialEq)]
pub struct QueryHint {
    context_id: ContentHash,
    grid: usize,
    purchase: Option<RankedPurchase>,
}

impl QueryHint {
    /// Ranked snapshot identity supporting this estimate.
    #[must_use]
    pub fn context_id(&self) -> ContentHash {
        self.context_id
    }

    /// Midpoint grid supporting this estimate.
    #[must_use]
    pub fn grid(&self) -> usize {
        self.grid
    }

    /// Estimated top purchase, or `None` when no sampled row changed the
    /// decision on this grid. `None` is not a proof that no purchase can help.
    #[must_use]
    pub fn purchase(&self) -> Option<&RankedPurchase> {
        self.purchase.as_ref()
    }

    /// Safe deterministic text. Identifiers are escaped and every finite
    /// scalar uses Rust's shortest round-tripping representation.
    #[must_use]
    pub fn render_text(&self) -> String {
        match &self.purchase {
            Some(top) => format!(
                "estimated top evidence on the supplied menu from a {}-point midpoint sweep: {} (${}) - sampled flip fraction {} -> {} on {} ({}/$)",
                self.grid,
                escape_text(&top.probe.name),
                top.probe.cost,
                top.flip_before,
                top.flip_after,
                escape_text(&top.probe.target),
                top.score,
            ),
            None => format!(
                "no sampled purchase changed the decision on the {}-point midpoint sweep; this estimate does not prove that further evidence has zero value",
                self.grid
            ),
        }
    }

    /// Strict JSON rendering for logs and evidence payloads.
    #[must_use]
    pub fn to_json(&self) -> String {
        let context = self.context_id.to_hex();
        match &self.purchase {
            Some(top) => format!(
                "{{\"schema\":\"fs-plan.voi-hint.v1\",\"kind\":\"estimated_purchase\",\"context\":\"{context}\",\"grid\":{},\"probe\":{},\"target\":{},\"cost_dollars\":{},\"sampled_flip_before\":{},\"sampled_flip_after\":{},\"score_per_dollar\":{}}}",
                self.grid,
                json_string(&top.probe.name),
                json_string(&top.probe.target),
                top.probe.cost,
                top.flip_before,
                top.flip_after,
                top.score,
            ),
            None => format!(
                "{{\"schema\":\"fs-plan.voi-hint.v1\",\"kind\":\"no_sampled_change\",\"context\":\"{context}\",\"grid\":{},\"authoritative_zero\":false}}",
                self.grid
            ),
        }
    }
}

impl core::fmt::Display for QueryHint {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(&self.render_text())
    }
}

fn escape_text(value: &str) -> String {
    value.chars().flat_map(char::escape_default).collect()
}

fn json_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for byte in value.bytes() {
        match byte {
            b'"' => out.push_str("\\\""),
            b'\\' => out.push_str("\\\\"),
            _ => out.push(char::from(byte)),
        }
    }
    out.push('"');
    out
}

fn compare_ranked(a: &RankedPurchase, b: &RankedPurchase) -> core::cmp::Ordering {
    b.score
        .total_cmp(&a.score)
        .then(a.probe.cost.total_cmp(&b.probe.cost))
        .then(a.probe.name.cmp(&b.probe.name))
}

fn validate_probe(index: usize, probe: &Probe) -> Result<(), VoiError> {
    validate_name("probe", index, &probe.name)?;
    validate_name("probe target", index, &probe.target)?;
    if !probe.cost.is_finite() || probe.cost <= 0.0 {
        return Err(VoiError::InvalidProbeValue {
            probe: probe.name.clone(),
            field: "cost",
            value: probe.cost,
        });
    }
    if !probe.shrink.is_finite() || probe.shrink <= 0.0 || probe.shrink >= 1.0 {
        return Err(VoiError::InvalidProbeValue {
            probe: probe.name.clone(),
            field: "shrink",
            value: probe.shrink,
        });
    }
    Ok(())
}

fn validate_menu(nodes: &[UncertaintyNode], menu: &[Probe]) -> Result<Vec<usize>, VoiError> {
    validate_size("probe menu", menu.len(), 1, MAX_VOI_PROBES)?;
    let mut names = BTreeSet::new();
    let mut targets = Vec::with_capacity(menu.len());
    for (index, probe) in menu.iter().enumerate() {
        validate_probe(index, probe)?;
        if !names.insert(probe.name.as_str()) {
            return Err(VoiError::DuplicateName {
                kind: "probe",
                name: probe.name.clone(),
            });
        }
        let mut matched = None;
        let mut matches = 0usize;
        for (node_idx, node) in nodes.iter().enumerate() {
            if node.name == probe.target {
                matches += 1;
                matched = Some(node_idx);
            }
        }
        let Some(node_idx) = matched.filter(|_| matches == 1) else {
            return Err(VoiError::TargetResolution {
                probe: probe.name.clone(),
                target: probe.target.clone(),
                matches,
            });
        };
        targets.push(node_idx);
    }
    Ok(targets)
}

fn push_u32(out: &mut Vec<u8>, value: usize, subject: &'static str) -> Result<(), VoiError> {
    let value = u32::try_from(value).map_err(|_| VoiError::ArithmeticRefusal {
        operation: "VoI context length",
        subject: subject.to_string(),
    })?;
    out.extend_from_slice(&value.to_le_bytes());
    Ok(())
}

fn push_text(out: &mut Vec<u8>, value: &str, subject: &'static str) -> Result<(), VoiError> {
    push_u32(out, value.len(), subject)?;
    out.extend_from_slice(value.as_bytes());
    Ok(())
}

fn ranked_source_context(
    nodes: &[UncertaintyNode],
    menu: &[Probe],
    grid: usize,
    policy_scope: &str,
    snapshot_id: &str,
) -> Result<ContentHash, VoiError> {
    let mut canonical = Vec::new();
    canonical.extend_from_slice(&1u32.to_le_bytes());
    push_text(&mut canonical, policy_scope, "VoI policy scope")?;
    push_text(&mut canonical, snapshot_id, "VoI snapshot identity")?;
    push_u32(&mut canonical, grid, "grid")?;
    push_u32(&mut canonical, nodes.len(), "uncertainty nodes")?;
    for node in nodes {
        push_text(&mut canonical, &node.name, "uncertainty node name")?;
        canonical.extend_from_slice(&node.lo.to_bits().to_le_bytes());
        canonical.extend_from_slice(&node.nominal.to_bits().to_le_bytes());
        canonical.extend_from_slice(&node.hi.to_bits().to_le_bytes());
    }
    let mut canonical_menu: Vec<&Probe> = menu.iter().collect();
    canonical_menu.sort_by(|left, right| left.name.cmp(&right.name));
    push_u32(&mut canonical, canonical_menu.len(), "probe menu")?;
    for probe in canonical_menu {
        push_text(&mut canonical, &probe.name, "probe name")?;
        push_text(&mut canonical, &probe.target, "probe target")?;
        canonical.extend_from_slice(&probe.cost.to_bits().to_le_bytes());
        canonical.extend_from_slice(&probe.shrink.to_bits().to_le_bytes());
        canonical.push(match probe.kind {
            ProbeKind::Computational => 0,
            ProbeKind::Physical => 1,
        });
    }
    Ok(hash_domain(RANKED_MENU_SOURCE_DOMAIN, &canonical))
}

fn ranked_output_context(
    source_context_id: ContentHash,
    rows: &[RankedPurchase],
) -> Result<ContentHash, VoiError> {
    let mut canonical = Vec::new();
    canonical.extend_from_slice(&2u32.to_le_bytes());
    canonical.extend_from_slice(source_context_id.as_bytes());
    push_u32(&mut canonical, rows.len(), "ranked output rows")?;
    for row in rows {
        push_text(&mut canonical, &row.probe.name, "ranked probe name")?;
        canonical.extend_from_slice(&row.flip_before.to_bits().to_le_bytes());
        canonical.extend_from_slice(&row.flip_after.to_bits().to_le_bytes());
        canonical.extend_from_slice(&row.score.to_bits().to_le_bytes());
    }
    Ok(hash_domain(RANKED_MENU_CONTEXT_DOMAIN, &canonical))
}

#[derive(Debug, Clone, Copy)]
struct PreparedProbe {
    node_idx: usize,
    post_lo: f64,
    post_hi: f64,
}

fn prepare_probes(
    nodes: &[UncertaintyNode],
    menu: &[Probe],
    targets: &[usize],
) -> Result<Vec<PreparedProbe>, VoiError> {
    let mut prepared = Vec::with_capacity(menu.len());
    for (probe, &node_idx) in menu.iter().zip(targets) {
        let node = node_at(nodes, node_idx)?;
        let contracted_left = (node.nominal - node.lo) * probe.shrink;
        let contracted_right = (node.hi - node.nominal) * probe.shrink;
        let post_lo = node.nominal - contracted_left;
        let post_hi = node.nominal + contracted_right;
        let post_width = post_hi - post_lo;
        let expected_width = (node.hi - node.lo) * probe.shrink;
        if !contracted_left.is_finite()
            || !contracted_right.is_finite()
            || !post_lo.is_finite()
            || !post_hi.is_finite()
            || !post_width.is_finite()
            || !expected_width.is_finite()
            || (node.nominal > node.lo && contracted_left == 0.0)
            || (node.hi > node.nominal && contracted_right == 0.0)
            || post_lo < node.lo
            || post_lo > node.nominal
            || post_hi < node.nominal
            || post_hi > node.hi
            || (node.hi > node.lo && post_width == 0.0)
        {
            return Err(VoiError::ArithmeticRefusal {
                operation: "post-probe interval contraction",
                subject: probe.name.clone(),
            });
        }
        prepared.push(PreparedProbe {
            node_idx,
            post_lo,
            post_hi,
        });
    }
    Ok(prepared)
}

/// Rank the probe menu by sampled flip-fraction reduction per dollar for the live
/// decision — MYOPIC one-step VoI (each probe is evaluated against the
/// CURRENT state only; no sequential tree).
///
/// # Errors
/// [`VoiError`] when the decision, node set, menu, targets, grid, probe
/// economics, callback margins, or derived arithmetic are invalid.
pub fn rank_purchases(
    decision: &LiveDecision<'_>,
    nodes: &[UncertaintyNode],
    menu: &[Probe],
    grid: usize,
    policy_scope: &str,
    snapshot_id: &str,
) -> Result<RankedMenu, VoiError> {
    validate_name("VoI policy scope", 0, policy_scope)?;
    validate_name("VoI snapshot identity", 0, snapshot_id)?;
    validate_nodes(decision, nodes)?;
    validate_grid(grid)?;
    let targets = validate_menu(nodes, menu)?;
    let unique_targets = targets.iter().copied().collect::<BTreeSet<_>>();
    let sweep_count = unique_targets
        .len()
        .checked_add(menu.len())
        .ok_or_else(|| VoiError::ArithmeticRefusal {
            operation: "ranking sweep count",
            subject: "probe menu".to_string(),
        })?;
    let evaluations = grid
        .checked_mul(sweep_count)
        .and_then(|sweeps| sweeps.checked_add(1))
        .ok_or_else(|| VoiError::ArithmeticRefusal {
            operation: "ranking evaluation count",
            subject: "probe menu".to_string(),
        })?;
    validate_evaluations(evaluations)?;
    // All input-derived intervals are prepared before the first callback, so
    // a malformed later probe cannot leave observable partial callback work.
    let prepared = prepare_probes(nodes, menu, &targets)?;
    let source_context_id = ranked_source_context(nodes, menu, grid, policy_scope, snapshot_id)?;
    let base = nominal_verdict_validated(decision, nodes)?;
    let mut flip_before = vec![None; nodes.len()];
    for node_idx in unique_targets {
        let node = node_at(nodes, node_idx)?;
        flip_before[node_idx] = Some(flip_probability_validated(
            decision, nodes, base, node_idx, node.lo, node.hi, grid,
        )?);
    }

    let mut ranked = Vec::with_capacity(menu.len());
    let mut evaluation_order = (0..menu.len()).collect::<Vec<_>>();
    evaluation_order.sort_by(|left, right| menu[*left].name.cmp(&menu[*right].name));
    for index in evaluation_order {
        let probe = &menu[index];
        let prepared = prepared[index];
        let before = flip_before[prepared.node_idx].ok_or_else(|| VoiError::ArithmeticRefusal {
            operation: "pre-probe sweep lookup",
            subject: probe.name.clone(),
        })?;
        let flip_after = flip_probability_validated(
            decision,
            nodes,
            base,
            prepared.node_idx,
            prepared.post_lo,
            prepared.post_hi,
            grid,
        )?;
        let score = (before - flip_after).max(0.0) / probe.cost;
        if !score.is_finite() || score < 0.0 {
            return Err(VoiError::ArithmeticRefusal {
                operation: "sampled flip-fraction-per-dollar score",
                subject: probe.name.clone(),
            });
        }
        ranked.push(RankedPurchase {
            probe: probe.clone(),
            flip_before: before,
            flip_after,
            score,
        });
    }
    ranked.sort_by(compare_ranked);
    let context_id = ranked_output_context(source_context_id, &ranked)?;
    Ok(RankedMenu {
        rows: ranked,
        source_context_id,
        context_id,
        policy_scope: policy_scope.to_string(),
        snapshot_id: snapshot_id.to_string(),
        grid,
    })
}

/// Surface a structured QUERY-RESULT HINT. Every scalar is explicitly a
/// grid-sampled estimate; a sampled zero is never rendered as proof that no
/// evidence could change the decision.
#[must_use]
pub fn hint_for_query(ranked: &RankedMenu) -> QueryHint {
    QueryHint {
        context_id: ranked.context_id,
        grid: ranked.grid,
        purchase: ranked.rows.iter().find(|row| row.score > 0.0).cloned(),
    }
}

/// The audit verdict for reporting. This enum is not scheduling authority;
/// only the live [`VoiScheduler`] owns executable state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AuditVerdict {
    /// Anytime-valid evidence crossed the fixed activation threshold.
    KeepScheduling,
    /// Evidence is absent, insufficient, or has not crossed the threshold.
    DemoteToReporting,
}

/// One validated matched-cost prospective-audit observation.
///
/// Fields are private so raw booleans and unmatched prices cannot enter the
/// e-process without identity, provenance, and economic validation.
#[derive(Debug, Clone, PartialEq)]
pub struct MatchedAuditRecord {
    observation_id: String,
    recommended_id: String,
    alternative_id: String,
    provenance: String,
    matched_cost: f64,
    recommended_changed_decision: bool,
    alternative_changed_decision: bool,
}

impl MatchedAuditRecord {
    /// Construct one matched-cost comparison.
    ///
    /// # Errors
    /// [`VoiError`] unless identities/provenance are bounded visible ASCII,
    /// candidates differ, and both finite positive costs are bit-identical.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        observation_id: impl Into<String>,
        recommended_id: impl Into<String>,
        alternative_id: impl Into<String>,
        provenance: impl Into<String>,
        recommended_cost: f64,
        alternative_cost: f64,
        recommended_changed_decision: bool,
        alternative_changed_decision: bool,
    ) -> Result<Self, VoiError> {
        let observation_id = observation_id.into();
        let recommended_id = recommended_id.into();
        let alternative_id = alternative_id.into();
        let provenance = provenance.into();
        for (kind, value) in [
            ("audit observation", observation_id.as_str()),
            ("recommended purchase", recommended_id.as_str()),
            ("alternative purchase", alternative_id.as_str()),
            ("audit provenance", provenance.as_str()),
        ] {
            validate_name(kind, 0, value)?;
        }
        if recommended_id == alternative_id {
            return Err(VoiError::InvalidAuditPair {
                observation: observation_id,
            });
        }
        if !recommended_cost.is_finite()
            || recommended_cost <= 0.0
            || recommended_cost.to_bits() != alternative_cost.to_bits()
        {
            return Err(VoiError::InvalidAuditCost {
                observation: observation_id,
                recommended_cost,
                alternative_cost,
            });
        }
        // Caller-owned Strings can carry arbitrary spare capacity despite
        // bounded content. Rebuild every retained identity from its validated
        // slice so the audit-record cap is also an operational memory bound.
        let observation_id = observation_id.as_str().to_owned();
        let recommended_id = recommended_id.as_str().to_owned();
        let alternative_id = alternative_id.as_str().to_owned();
        let provenance = provenance.as_str().to_owned();
        Ok(Self {
            observation_id,
            recommended_id,
            alternative_id,
            provenance,
            matched_cost: recommended_cost,
            recommended_changed_decision,
            alternative_changed_decision,
        })
    }

    /// Stable observation identity used to prevent duplicate evidence.
    #[must_use]
    pub fn observation_id(&self) -> &str {
        &self.observation_id
    }

    /// Recommended-purchase identity.
    #[must_use]
    pub fn recommended_id(&self) -> &str {
        &self.recommended_id
    }

    /// Matched alternative-purchase identity.
    #[must_use]
    pub fn alternative_id(&self) -> &str {
        &self.alternative_id
    }

    /// Caller-supplied provenance identity.
    #[must_use]
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// Exact matched cost.
    #[must_use]
    pub fn matched_cost(&self) -> f64 {
        self.matched_cost
    }

    /// Whether the recommended purchase changed the realized decision.
    #[must_use]
    pub fn recommended_changed_decision(&self) -> bool {
        self.recommended_changed_decision
    }

    /// Whether the matched alternative changed the realized decision.
    #[must_use]
    pub fn alternative_changed_decision(&self) -> bool {
        self.alternative_changed_decision
    }
}

/// One authorized, single-epoch scheduling decision. The receipt retains the
/// ranked snapshot, audit evidence root, and exact budget transition instead of
/// returning a provenance-free probe row.
#[derive(Debug, PartialEq)]
pub struct ScheduledPurchase {
    purchase: RankedPurchase,
    ranked_context_id: ContentHash,
    ranked_source_context_id: ContentHash,
    policy_scope: String,
    snapshot_id: String,
    ranked_grid: usize,
    audit_context_id: ContentHash,
    audit_observations: usize,
    audit_log_e_value: f64,
    budget_dollars: f64,
    remaining_budget_dollars: f64,
}

impl ScheduledPurchase {
    /// Authorized purchase.
    #[must_use]
    pub fn purchase(&self) -> &RankedPurchase {
        &self.purchase
    }

    /// Ranked node/menu/grid snapshot identity.
    #[must_use]
    pub fn ranked_context_id(&self) -> ContentHash {
        self.ranked_context_id
    }

    /// Source policy/snapshot/node/menu/grid identity before evaluation.
    #[must_use]
    pub fn ranked_source_context_id(&self) -> ContentHash {
        self.ranked_source_context_id
    }

    /// Scheduling policy/version scope shared by the audit and ranking.
    #[must_use]
    pub fn policy_scope(&self) -> &str {
        &self.policy_scope
    }

    /// Decision/ledger snapshot whose ranking was consumed.
    #[must_use]
    pub fn snapshot_id(&self) -> &str {
        &self.snapshot_id
    }

    /// Midpoint grid supporting the sampled purchase score.
    #[must_use]
    pub fn ranked_grid(&self) -> usize {
        self.ranked_grid
    }

    /// Anytime-valid matched-audit evidence identity.
    #[must_use]
    pub fn audit_context_id(&self) -> ContentHash {
        self.audit_context_id
    }

    /// Matched-cost observation count supporting authority.
    #[must_use]
    pub fn audit_observations(&self) -> usize {
        self.audit_observations
    }

    /// Final log e-value supporting authority.
    #[must_use]
    pub fn audit_log_e_value(&self) -> f64 {
        self.audit_log_e_value
    }

    /// Admitted scheduling budget in dollars.
    #[must_use]
    pub fn budget_dollars(&self) -> f64 {
        self.budget_dollars
    }

    /// Exact remaining budget in dollars after this one purchase.
    #[must_use]
    pub fn remaining_budget_dollars(&self) -> f64 {
        self.remaining_budget_dollars
    }
}

/// Immutable reporting snapshot of one live prospective-audit prefix. This is
/// never scheduling authority; only the owning [`VoiScheduler`] can spend.
#[derive(Debug, Clone, PartialEq)]
pub struct AuditReport {
    policy_scope: String,
    audit_context_id: ContentHash,
    observations: usize,
    log_e_value: f64,
    verdict: AuditVerdict,
}

impl AuditReport {
    /// Reporting verdict.
    #[must_use]
    pub fn verdict(&self) -> AuditVerdict {
        self.verdict
    }

    /// Caller-declared audit policy/version scope.
    #[must_use]
    pub fn policy_scope(&self) -> &str {
        &self.policy_scope
    }

    /// Content identity of the canonical evidence prefix.
    #[must_use]
    pub fn audit_context_id(&self) -> ContentHash {
        self.audit_context_id
    }

    /// Number of matched-cost observations evaluated.
    #[must_use]
    pub fn observations(&self) -> usize {
        self.observations
    }

    /// Final log e-value, useful for reporting progress before activation.
    #[must_use]
    pub fn log_e_value(&self) -> f64 {
        self.log_e_value
    }
}

fn audit_context(
    policy_scope: &str,
    records: &[MatchedAuditRecord],
) -> Result<ContentHash, VoiError> {
    let mut canonical = Vec::new();
    canonical.extend_from_slice(&2u32.to_le_bytes());
    push_text(&mut canonical, policy_scope, "VoI audit policy scope")?;
    canonical.extend_from_slice(&VOI_AUDIT_ALPHA.to_bits().to_le_bytes());
    push_u32(
        &mut canonical,
        MAX_VOI_AUDIT_RECORDS,
        "maximum audit records",
    )?;
    push_u32(&mut canonical, records.len(), "audit records")?;
    for record in records {
        push_text(&mut canonical, &record.observation_id, "audit observation")?;
        push_text(
            &mut canonical,
            &record.recommended_id,
            "recommended purchase",
        )?;
        push_text(
            &mut canonical,
            &record.alternative_id,
            "alternative purchase",
        )?;
        push_text(&mut canonical, &record.provenance, "audit provenance")?;
        canonical.extend_from_slice(&record.matched_cost.to_bits().to_le_bytes());
        canonical.push(u8::from(record.recommended_changed_decision));
        canonical.push(u8::from(record.alternative_changed_decision));
    }
    Ok(hash_domain(AUDIT_CONTEXT_DOMAIN, &canonical))
}

/// Single-owner live VoI audit and purchase scheduler.
///
/// Audit observations enter one append-only [`PairwiseRace`] in prospective
/// order. The scheduler owns the remaining budget and the bounded set of
/// decision snapshots it has already evaluated, so one process cannot reuse a
/// stale ranking or reset the budget through the safe API. Caller-declared
/// identities and outcomes remain unauthenticated until `frankensim-wk4m`.
#[derive(Debug)]
pub struct VoiScheduler {
    policy_scope: String,
    remaining_budget_dollars: f64,
    audit_records: Vec<MatchedAuditRecord>,
    observation_ids: BTreeSet<String>,
    race: PairwiseRace,
    consumed_snapshots: BTreeSet<String>,
}

impl VoiScheduler {
    /// Create one live scheduler for a fixed policy/version and total budget.
    ///
    /// # Errors
    /// [`VoiError`] when the policy identity or budget is malformed.
    pub fn new(policy_scope: impl Into<String>, budget: f64) -> Result<Self, VoiError> {
        let policy_scope = policy_scope.into();
        validate_name("VoI policy scope", 0, &policy_scope)?;
        if !budget.is_finite() || budget < 0.0 {
            return Err(VoiError::InvalidBudget { budget });
        }
        Ok(Self {
            policy_scope: policy_scope.as_str().to_owned(),
            remaining_budget_dollars: budget,
            audit_records: Vec::new(),
            observation_ids: BTreeSet::new(),
            race: PairwiseRace::new(LossSpan::ONE),
            consumed_snapshots: BTreeSet::new(),
        })
    }

    /// Append one prospectively ordered matched-cost result to the one live
    /// e-process. Duplicate identities and limit+1 refuse before wealth or
    /// retained state changes.
    ///
    /// # Errors
    /// [`VoiError`] on a duplicate observation, the audit cap, or invalid
    /// e-process arithmetic.
    pub fn observe_audit(&mut self, record: MatchedAuditRecord) -> Result<(), VoiError> {
        if self.audit_records.len() >= MAX_VOI_AUDIT_RECORDS {
            return Err(VoiError::SizeLimit {
                collection: "VoI audit records",
                count: self.audit_records.len().saturating_add(1),
                min: 0,
                max: MAX_VOI_AUDIT_RECORDS,
            });
        }
        if self.observation_ids.contains(&record.observation_id) {
            return Err(VoiError::DuplicateAuditObservation {
                observation: record.observation_id.clone(),
            });
        }
        let recommended_loss = f64::from(u8::from(!record.recommended_changed_decision));
        let alternative_loss = f64::from(u8::from(!record.alternative_changed_decision));
        let mut next_race = self.race.clone();
        next_race
            .observe(recommended_loss, alternative_loss)
            .map_err(|_| VoiError::ArithmeticRefusal {
                operation: "VoI matched-cost e-process",
                subject: record.observation_id.clone(),
            })?;
        self.observation_ids.insert(record.observation_id.clone());
        self.audit_records.push(record);
        self.race = next_race;
        Ok(())
    }

    /// Immutable reporting snapshot of the current chronological audit prefix.
    /// The report carries no scheduling capability.
    ///
    /// # Errors
    /// [`VoiError`] if the e-process or content identity leaves its finite
    /// bounded domain.
    pub fn audit_report(&self) -> Result<AuditReport, VoiError> {
        let audit_context_id = audit_context(&self.policy_scope, &self.audit_records)?;
        let log_e_value = self.race.log_e_value();
        if !log_e_value.is_finite() {
            return Err(VoiError::ArithmeticRefusal {
                operation: "VoI audit log e-value",
                subject: audit_context_id.to_hex(),
            });
        }
        Ok(AuditReport {
            policy_scope: self.policy_scope.clone(),
            audit_context_id,
            observations: self.audit_records.len(),
            log_e_value,
            verdict: if self.race.a_beats_b(VOI_AUDIT_ALPHA) {
                AuditVerdict::KeepScheduling
            } else {
                AuditVerdict::DemoteToReporting
            },
        })
    }

    /// Fixed policy/version scope owned by this scheduler.
    #[must_use]
    pub fn policy_scope(&self) -> &str {
        &self.policy_scope
    }

    /// Current unspent scheduler budget.
    #[must_use]
    pub fn remaining_budget_dollars(&self) -> f64 {
        self.remaining_budget_dollars
    }

    /// Number of decision snapshots already evaluated by this scheduler.
    #[must_use]
    pub fn consumed_snapshots(&self) -> usize {
        self.consumed_snapshots.len()
    }

    /// Evaluate at most one purchase from one previously unseen decision
    /// snapshot. The current live audit must still be above threshold. All
    /// validation and arithmetic precede mutation; success (including a
    /// no-affordable-purchase result) consumes the snapshot atomically.
    ///
    /// # Errors
    /// [`VoiError`] when audit authority is currently absent, policy scopes
    /// differ, the snapshot was already consumed, retained snapshot capacity
    /// is exhausted, or budget arithmetic cannot decrease monotonically.
    pub fn schedule(&mut self, ranked: RankedMenu) -> Result<Option<ScheduledPurchase>, VoiError> {
        if ranked.policy_scope != self.policy_scope {
            return Err(VoiError::PolicyScopeMismatch {
                expected: self.policy_scope.clone(),
                actual: ranked.policy_scope,
            });
        }
        if self.consumed_snapshots.contains(&ranked.snapshot_id) {
            return Err(VoiError::RankingSnapshotAlreadyConsumed {
                snapshot_id: ranked.snapshot_id,
            });
        }
        if self.consumed_snapshots.len() >= MAX_VOI_SCHEDULED_CONTEXTS {
            return Err(VoiError::SizeLimit {
                collection: "VoI consumed ranking snapshots",
                count: self.consumed_snapshots.len().saturating_add(1),
                min: 0,
                max: MAX_VOI_SCHEDULED_CONTEXTS,
            });
        }
        if !self.race.a_beats_b(VOI_AUDIT_ALPHA) {
            return Err(VoiError::MissingSchedulingAuthority);
        }
        let audit = self.audit_report()?;
        let budget = self.remaining_budget_dollars;
        let purchase = ranked
            .rows
            .iter()
            .find(|row| row.score > 0.0 && row.probe.cost <= budget)
            .cloned();
        let remaining = if let Some(purchase) = &purchase {
            let remaining = budget - purchase.probe.cost;
            if !remaining.is_finite() || remaining < 0.0 || remaining >= budget {
                return Err(VoiError::ArithmeticRefusal {
                    operation: "remaining-budget subtraction",
                    subject: purchase.probe.name.clone(),
                });
            }
            remaining
        } else {
            budget
        };

        self.consumed_snapshots.insert(ranked.snapshot_id.clone());
        self.remaining_budget_dollars = remaining;
        Ok(purchase.map(|purchase| ScheduledPurchase {
            purchase,
            ranked_context_id: ranked.context_id,
            ranked_source_context_id: ranked.source_context_id,
            policy_scope: ranked.policy_scope,
            snapshot_id: ranked.snapshot_id,
            ranked_grid: ranked.grid,
            audit_context_id: audit.audit_context_id,
            audit_observations: audit.observations,
            audit_log_e_value: audit.log_e_value,
            budget_dollars: budget,
            remaining_budget_dollars: remaining,
        }))
    }
}

/// Recompute a bounded reporting-only audit from a supplied chronological
/// prefix. This helper never returns scheduling authority; executable callers
/// must retain one live [`VoiScheduler`] and append observations to it.
///
/// # Errors
/// [`VoiError`] for malformed policy identity, duplicate/oversized records, or
/// invalid e-process arithmetic.
pub fn audit_scheduling(
    policy_scope: &str,
    records: &[MatchedAuditRecord],
) -> Result<AuditReport, VoiError> {
    let mut scheduler = VoiScheduler::new(policy_scope, 0.0)?;
    for record in records {
        scheduler.observe_audit(record.clone())?;
    }
    scheduler.audit_report()
}

#[cfg(test)]
mod retained_capacity_tests {
    use super::*;

    fn spare_capacity(value: &str) -> String {
        let mut out = String::with_capacity(4096);
        out.push_str(value);
        out
    }

    #[test]
    fn audit_authority_rebuilds_caller_owned_string_capacity() {
        let record = MatchedAuditRecord::new(
            spare_capacity("obs-1"),
            spare_capacity("recommended"),
            spare_capacity("alternative"),
            spare_capacity("ledger-row-1"),
            1.0,
            1.0,
            true,
            false,
        )
        .expect("bounded matched audit record");
        for value in [
            &record.observation_id,
            &record.recommended_id,
            &record.alternative_id,
            &record.provenance,
        ] {
            assert!(value.capacity() <= MAX_VOI_NAME_BYTES);
        }

        let scheduler =
            VoiScheduler::new(spare_capacity("policy-v1"), 1.0).expect("bounded scheduler policy");
        assert!(scheduler.policy_scope.capacity() <= MAX_VOI_NAME_BYTES);
    }
}
