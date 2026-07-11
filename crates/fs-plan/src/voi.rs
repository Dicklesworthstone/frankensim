//! VALUE-OF-INFORMATION QUERIES (addendum Proposal C, bead knh1.6;
//! [F] — behind the `voi-queries` feature): THE IGNORANCE MARKET, v0
//! as a RANKED LIST. Across everything the ledger is uncertain about,
//! where does one dollar of evidence most change the downstream
//! decision? Decision sensitivity (does the decision FLIP inside the
//! node's interval?) is computed from CACHED surrogate sweeps, crossed
//! with a PRICED PROBE MENU that unifies computational and physical
//! experiments, ranked by FLIP-PROBABILITY-PER-DOLLAR.
//!
//! MYOPIC one-step VoI ONLY (the proposal's own discipline: full
//! sequential VoI is intractable and myopic captures most of the
//! value). The output surfaces as (i) the hint on query results — the
//! upgrade the fs-ir anytime module's CONTRACT reserved for Proposal C
//! — and (ii) the scheduler for discrepancy probes.
//!
//! THE KILL CRITERION AS CODE: [`audit_verdict`] compares
//! VoI-recommended purchases against agent-chosen alternatives at
//! matched cost; if recommendations do not measurably outperform on
//! realized decision changes, VoI DEMOTES ITSELF to a reporting
//! feature.

use std::collections::BTreeSet;

/// Maximum uncertainty nodes in one myopic VoI request.
pub const MAX_VOI_NODES: usize = 256;
/// Maximum probe menu entries (and scheduled ranked entries).
pub const MAX_VOI_PROBES: usize = 1024;
/// Maximum UTF-8 byte length of node, probe, and target names.
pub const MAX_VOI_NAME_BYTES: usize = 128;
/// Maximum interval-sweep grid size.
pub const MAX_VOI_GRID: usize = 1024;
/// Maximum surrogate evaluations admitted by one public VoI call.
pub const MAX_VOI_EVALUATIONS: usize = 4096;

/// Why a VoI query or schedule was refused.
#[derive(Debug, Clone, PartialEq)]
pub enum VoiError {
    /// A bounded collection is empty or oversized.
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
    /// A node/probe/target name is blank, padded, or oversized.
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
    /// A ranked purchase contains a malformed derived scalar.
    InvalidRankedValue {
        /// Probe name.
        probe: String,
        /// Invalid field.
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
    /// A forged ranked menu repeats a purchase identity.
    DuplicateRankedProbe {
        /// Repeated probe name.
        name: String,
    },
    /// The scheduling budget is nonfinite or negative.
    InvalidBudget {
        /// Supplied budget.
        budget: f64,
    },
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
                "{kind} name at index {index} is blank, padded, or {bytes} bytes long (limit {max_bytes})"
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
            Self::InvalidRankedValue {
                probe,
                field,
                value,
            } => write!(f, "ranked probe {probe:?} has invalid {field} {value:?}"),
            Self::DuplicateRankedProbe { name } => {
                write!(f, "ranked probe {name:?} appears more than once")
            }
            Self::InvalidBudget { budget } => {
                write!(
                    f,
                    "probe budget must be finite and non-negative, got {budget:?}"
                )
            }
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
    if name.is_empty() || name.trim() != name || name.len() > MAX_VOI_NAME_BYTES {
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
    node_idx: usize,
    lo: f64,
    hi: f64,
    grid: usize,
) -> Result<f64, VoiError> {
    let base = nominal_verdict_validated(decision, nodes)?;
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
    /// return the fraction of the interval where the verdict differs
    /// from nominal — the myopic flip probability under the uniform
    /// interval measure (v0's declared prior).
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
        flip_probability_validated(self, nodes, node_idx, node.lo, node.hi, grid)
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
    /// Flip probability before the purchase.
    flip_before: f64,
    /// Expected flip probability after the purchase.
    flip_after: f64,
    /// THE SCORE: expected flip-probability reduction per dollar.
    score: f64,
}

impl RankedPurchase {
    /// The validated probe purchase.
    #[must_use]
    pub fn probe(&self) -> &Probe {
        &self.probe
    }

    /// Flip probability before the purchase.
    #[must_use]
    pub fn flip_before(&self) -> f64 {
        self.flip_before
    }

    /// Expected flip probability after the purchase.
    #[must_use]
    pub fn flip_after(&self) -> f64 {
        self.flip_after
    }

    /// Expected flip-probability reduction per dollar.
    #[must_use]
    pub fn score(&self) -> f64 {
        self.score
    }
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

/// Rank the probe menu by flip-probability-per-dollar for the live
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
) -> Result<Vec<RankedPurchase>, VoiError> {
    validate_nodes(decision, nodes)?;
    validate_grid(grid)?;
    let targets = validate_menu(nodes, menu)?;
    let evaluations = grid
        .checked_add(1)
        .and_then(|per_sweep| per_sweep.checked_mul(2))
        .and_then(|per_probe| menu.len().checked_mul(per_probe))
        .ok_or_else(|| VoiError::ArithmeticRefusal {
            operation: "ranking evaluation count",
            subject: "probe menu".to_string(),
        })?;
    validate_evaluations(evaluations)?;

    let mut ranked = Vec::with_capacity(menu.len());
    for (probe, &node_idx) in menu.iter().zip(&targets) {
        let node = node_at(nodes, node_idx)?;
        let flip_before =
            flip_probability_validated(decision, nodes, node_idx, node.lo, node.hi, grid)?;
        // Myopic post-probe state: the interval shrinks around the
        // nominal by the probe's factor.
        let half = (node.hi - node.lo) / 2.0 * probe.shrink;
        let post_lo = node.nominal - half;
        let post_hi = node.nominal + half;
        if !half.is_finite()
            || !post_lo.is_finite()
            || !post_hi.is_finite()
            || !(post_hi - post_lo).is_finite()
            || post_lo > node.nominal
            || post_hi < node.nominal
        {
            return Err(VoiError::ArithmeticRefusal {
                operation: "post-probe interval",
                subject: probe.name.clone(),
            });
        }
        let flip_after =
            flip_probability_validated(decision, nodes, node_idx, post_lo, post_hi, grid)?;
        let score = (flip_before - flip_after).max(0.0) / probe.cost;
        if !score.is_finite() || score < 0.0 {
            return Err(VoiError::ArithmeticRefusal {
                operation: "flip-probability-per-dollar score",
                subject: probe.name.clone(),
            });
        }
        ranked.push(RankedPurchase {
            probe: probe.clone(),
            flip_before,
            flip_after,
            score,
        });
    }
    ranked.sort_by(compare_ranked);
    Ok(ranked)
}

fn canonical_ranked(ranked: &[RankedPurchase]) -> Result<Vec<&RankedPurchase>, VoiError> {
    validate_size("ranked probe menu", ranked.len(), 1, MAX_VOI_PROBES)?;
    let mut names = BTreeSet::new();
    let mut canonical = Vec::with_capacity(ranked.len());
    for (index, purchase) in ranked.iter().enumerate() {
        validate_probe(index, &purchase.probe)?;
        if !names.insert(purchase.probe.name.as_str()) {
            return Err(VoiError::DuplicateRankedProbe {
                name: purchase.probe.name.clone(),
            });
        }
        for (field, value, valid) in [
            (
                "flip_before",
                purchase.flip_before,
                purchase.flip_before.is_finite() && (0.0..=1.0).contains(&purchase.flip_before),
            ),
            (
                "flip_after",
                purchase.flip_after,
                purchase.flip_after.is_finite() && (0.0..=1.0).contains(&purchase.flip_after),
            ),
            (
                "score",
                purchase.score,
                purchase.score.is_finite() && purchase.score >= 0.0,
            ),
        ] {
            if !valid {
                return Err(VoiError::InvalidRankedValue {
                    probe: purchase.probe.name.clone(),
                    field,
                    value,
                });
            }
        }
        canonical.push(purchase);
    }
    canonical.sort_by(|a, b| compare_ranked(a, b));
    Ok(canonical)
}

/// Surface the top purchase as the QUERY-RESULT HINT — the Proposal-8
/// anytime-hint shape, now priced by decision impact instead of the
/// O(h) extrapolation the fs-ir CONTRACT flagged as interim.
///
/// Input order is not authoritative: the validated score/cost/name
/// comparator is re-applied before selecting a hint.
///
/// # Errors
/// [`VoiError`] when the ranked menu is empty, oversized, duplicated,
/// or internally malformed.
pub fn hint_for_query(ranked: &[RankedPurchase]) -> Result<String, VoiError> {
    let canonical = canonical_ranked(ranked)?;
    let hint = match canonical.first() {
        Some(top) if top.score > 0.0 => format!(
            "highest-value evidence: {} (${:.0}) — expected flip-probability drop \
             {:.3} -> {:.3} on '{}' ({:.4}/$)",
            top.probe.name,
            top.probe.cost,
            top.flip_before,
            top.flip_after,
            top.probe.target,
            top.score
        ),
        _ => "no purchase on the menu changes the decision — spend nothing".to_string(),
    };
    Ok(hint)
}

/// THE PROBE SCHEDULER: greedy top-k purchases under a dollar budget
/// (the discrepancy-probe scheduler surface for color-probes).
///
/// # Errors
/// [`VoiError`] when the budget or a ranked purchase is malformed,
/// purchase identities repeat, or finite positive cost cannot decrease
/// the remaining budget monotonically.
pub fn schedule_probes(
    ranked: &[RankedPurchase],
    budget: f64,
) -> Result<Vec<RankedPurchase>, VoiError> {
    if !budget.is_finite() || budget < 0.0 {
        return Err(VoiError::InvalidBudget { budget });
    }
    let canonical = canonical_ranked(ranked)?;

    let mut remaining = budget;
    let mut out = Vec::new();
    for r in canonical {
        if r.score > 0.0 && r.probe.cost <= remaining {
            let next = remaining - r.probe.cost;
            if !next.is_finite() || next < 0.0 || next >= remaining {
                return Err(VoiError::ArithmeticRefusal {
                    operation: "remaining-budget subtraction",
                    subject: r.probe.name.clone(),
                });
            }
            remaining = next.max(0.0);
            out.push((*r).clone());
        }
    }
    Ok(out)
}

/// One prospective-audit record: a VoI-recommended purchase vs the
/// agent-chosen alternative at MATCHED COST, with the realized outcome
/// (did the evidence change the decision?).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuditRecord {
    /// Did the RECOMMENDED purchase realize a decision change?
    pub recommended_changed_decision: bool,
    /// Did the agent-chosen alternative?
    pub alternative_changed_decision: bool,
}

/// The audit verdict — the kill criterion as code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AuditVerdict {
    /// Recommendations measurably outperform: keep scheduling by VoI.
    KeepScheduling,
    /// They do not: DEMOTE VoI to a reporting feature.
    DemoteToReporting,
}

/// Compare realized decision-change rates at matched cost. VoI keeps
/// its scheduling authority only while it MEASURABLY outperforms.
#[must_use]
pub fn audit_verdict(records: &[AuditRecord]) -> AuditVerdict {
    if records.is_empty() {
        return AuditVerdict::DemoteToReporting; // no evidence, no authority
    }
    let rec = records
        .iter()
        .filter(|r| r.recommended_changed_decision)
        .count();
    let alt = records
        .iter()
        .filter(|r| r.alternative_changed_decision)
        .count();
    if rec > alt {
        AuditVerdict::KeepScheduling
    } else {
        AuditVerdict::DemoteToReporting
    }
}
