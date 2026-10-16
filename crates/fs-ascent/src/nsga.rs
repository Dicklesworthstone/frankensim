//! Baseline NSGA-III (epic bead `frankensim-epic-ascent-7tv.24.13`): a
//! deterministic gradient-free population optimizer whose complete
//! DEFECT-MATRIX battery is frozen BEFORE any adaptive-quality claim
//! rides on it.
//!
//! Design contract obeyed here (see crate CONTRACT.md house rules):
//! determinism first — every ordering decision flows through
//! `f64::total_cmp` plus an explicit index tie-break, so results are a
//! pure function of the input SET (enumeration-order independence is
//! asserted by the battery); refusals are typed and loud; duplicate
//! objective vectors collapse stably; budget accounting is integral
//! (evaluations never exceed the authored ceiling; a partial final
//! generation is discarded rather than half-reported); survivors keep
//! their non-dominated front membership in the receipt for audit.
//!
//! AUTHORITY: Estimate-class. This module provides no convergence
//! certificate and no quality claim about optima; those live with the
//! KKT-bearing tracing engines in [`crate::pareto`].

use fs_exec::{AdmittedBudget, BudgetConsumption, BudgetRefusal, Cx};

/// Authored baseline ceiling on reference-direction family size.
const REFERENCE_DIRECTION_CAP: usize = 10_000;

/// Typed refusals from NSGA-III setup and execution.
#[derive(Debug, Clone, PartialEq)]
pub enum NsgaError {
    /// The initial population is empty.
    PopulationEmpty,
    /// Individuals disagree on objective dimensionality.
    ObjectiveCountMismatch {
        /// Dimension of the first individual's objective vector.
        expected: usize,
        /// Offending individual position (0-based).
        at: usize,
        /// Its observed dimension.
        found: usize,
    },
    /// Decision vectors disagree in dimensionality.
    DecisionCountMismatch {
        /// Dimension of the first individual's decision vector.
        expected: usize,
        /// Offending individual position (0-based).
        at: usize,
        /// Its observed dimension.
        found: usize,
    },
    /// A non-finite objective or decision value entered the pipeline.
    NonFinite {
        /// Individual position (population position or evaluation id).
        individual: usize,
        /// Component position inside that vector.
        component: usize,
        /// Which vector tripped.
        kind: NonFiniteKind,
    },
    /// A box bound pair is malformed (min > max or non-finite) or the
    /// bounds vector length mismatches the decision axis count.
    BoxBoundsInvalid {
        /// Axis index of the offending bound.
        axis: usize,
    },
    /// The reference-direction request is malformed.
    ReferenceInvalid {
        /// Authored refusal reason.
        what: String,
    },
    /// An admitted reference-geometry operation was refused before it could
    /// publish a revised geometry or report.
    ReferenceAdmissionRefused {
        /// Stable fs-exec budget/cancellation refusal rendered for callers
        /// that do not otherwise depend on the executor error vocabulary.
        what: String,
    },
    /// A prepared reference snapshot was completed against a different
    /// geometry revision. Re-admit from the current geometry instead.
    ReferenceSnapshotStale,
    /// The authored evaluation budget cannot host even one complete
    /// generation, so nothing could ever converge honestly.
    BudgetInfeasible {
        /// Minimum individuals per generation.
        minimum_population: usize,
        /// Authored budget.
        budget: usize,
    },
}

/// Which vector carried the non-finite sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NonFiniteKind {
    /// Objective vector component.
    Objective,
    /// Decision vector component.
    Decision,
}

impl core::fmt::Display for NsgaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::PopulationEmpty => write!(f, "initial population is empty"),
            Self::ObjectiveCountMismatch {
                expected,
                at,
                found,
            } => write!(
                f,
                "individual {at} has {found} objectives; expected {expected}"
            ),
            Self::DecisionCountMismatch {
                expected,
                at,
                found,
            } => write!(
                f,
                "individual {at} decision dim {found}; expected {expected}"
            ),
            Self::NonFinite {
                individual,
                component,
                kind,
            } => match kind {
                NonFiniteKind::Objective => write!(
                    f,
                    "individual {individual} objective {component} is non-finite"
                ),
                NonFiniteKind::Decision => write!(
                    f,
                    "individual {individual} decision {component} is non-finite"
                ),
            },
            Self::BoxBoundsInvalid { axis } => {
                write!(f, "box bounds invalid on axis {axis} (finite min<=max)")
            }
            Self::ReferenceInvalid { what } => write!(f, "reference invalid: {what}"),
            Self::ReferenceAdmissionRefused { what } => {
                write!(f, "reference adaptation was not admitted: {what}")
            }
            Self::ReferenceSnapshotStale => write!(
                f,
                "reference adaptation snapshot no longer matches the supplied geometry"
            ),
            Self::BudgetInfeasible {
                minimum_population,
                budget,
            } => write!(
                f,
                "budget {budget} cannot host one generation of \
                 {minimum_population} individuals"
            ),
        }
    }
}

impl core::error::Error for NsgaError {}

fn vec_ordering(a: &[f64], b: &[f64]) -> core::cmp::Ordering {
    let n = a.len().min(b.len());
    for i in 0..n {
        match a[i].total_cmp(&b[i]) {
            core::cmp::Ordering::Equal => {}
            other => return other,
        }
    }
    a.len().cmp(&b.len())
}

/// One candidate: an uninterpreted decision vector plus its objective
/// vector. `x` exists so tie-breaks and duplicate policy key on
/// decision identity independent of objective bits.
#[derive(Debug, Clone, PartialEq)]
pub struct NsgaIndividual {
    /// Decision variables (box-checked when bounds are configured).
    pub x: Vec<f64>,
    /// Minimization objectives.
    pub f: Vec<f64>,
}

impl NsgaIndividual {
    /// Strict Pareto dominance under minimization: no worse everywhere,
    /// strictly better somewhere.
    #[must_use]
    pub fn dominates(&self, other: &Self) -> bool {
        let mut strictly_better = false;
        for (fi, fj) in self.f.iter().zip(other.f.iter()) {
            if fi > fj {
                return false;
            }
            if fi < fj {
                strictly_better = true;
            }
        }
        strictly_better
    }

    /// Canonical total order: objectives lexicographically by
    /// [`f64::total_cmp`], then decisions. Pure function of content.
    #[must_use]
    pub fn canonical_ordering(&self, other: &Self) -> core::cmp::Ordering {
        vec_ordering(&self.f, &other.f).then_with(|| vec_ordering(&self.x, &other.x))
    }
}

/// Fast non-dominated sorting O(N²·m). Fronts are ascending index
/// lists; two equal-objective individuals ride the SAME front (the
/// selection stage owns duplicate collapse, never the sort).
///
/// # Panics
/// Never: counter updates saturate and the loop terminates when every
/// rank-0 batch is drained; the internal invariant (acyclic dominance
/// assigns all indices) is asserted under debug builds only.
#[must_use]
pub fn fast_nondominated_sort(pop: &[NsgaIndividual]) -> Vec<Vec<usize>> {
    let n = pop.len();
    if n == 0 {
        return Vec::new();
    }
    let mut dominates_list = vec![Vec::<usize>::new(); n];
    let mut dominated_count = vec![0usize; n];
    for i in 0..n {
        for j in i + 1..n {
            let i_dom_j = pop[i].dominates(&pop[j]);
            let j_dom_i = pop[j].dominates(&pop[i]);
            if i_dom_j {
                dominates_list[i].push(j);
                dominated_count[j] += 1;
            } else if j_dom_i {
                dominates_list[j].push(i);
                dominated_count[i] += 1;
            }
        }
    }
    let mut fronts: Vec<Vec<usize>> = Vec::new();
    let mut assigned = 0usize;
    let mut current: Vec<usize> = (0..n).filter(|&i| dominated_count[i] == 0).collect();
    while !current.is_empty() {
        assigned += current.len();
        let mut next = Vec::new();
        for &i in &current {
            for &j in &dominates_list[i] {
                if dominated_count[j] > 0 {
                    dominated_count[j] -= 1;
                    if dominated_count[j] == 0 {
                        next.push(j);
                    }
                }
            }
        }
        fronts.push(current);
        current = next;
    }
    debug_assert_eq!(assigned, n, "acyclic dominance must cover all");
    fronts
}

/// Exhaustive m-tuples of non-negative integers summing to `p`, lex
/// ascending. Exact over integers; no accumulation order dependence.
#[must_use]
pub fn partition_tuples(p: usize, m: usize) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    if m == 0 {
        return out;
    }
    let mut prefix = vec![0usize; m];
    walk_partitions(m, 0, p, &mut prefix, &mut out);
    out
}

fn walk_partitions(
    m: usize,
    i: usize,
    left: usize,
    prefix: &mut [usize],
    out: &mut Vec<Vec<usize>>,
) {
    if i + 1 == m {
        prefix[i] = left;
        out.push(prefix.to_vec());
        return;
    }
    for v in 0..=left {
        prefix[i] = v;
        walk_partitions(m, i + 1, left - v, prefix, out);
    }
    prefix[i] = 0;
}

/// Closed-form family size C(p+m−1, m−1): the permutation oracle's
/// independent oracle. Saturates to `usize::MAX` past machine range
/// (the cap check treats that as refusal, never as overflow silence).
#[must_use]
pub fn das_dennis_cardinality(divisions: usize, m: usize) -> usize {
    let nn = divisions.saturating_add(m.saturating_sub(1));
    let k = m.saturating_sub(1);
    let mut acc = 1usize;
    for i in 0..k {
        acc = match acc.checked_mul(nn - i).and_then(|v| v.checked_div(i + 1)) {
            Some(v) => v,
            None => return usize::MAX,
        };
    }
    acc
}

fn das_dennis(divisions: usize, m: usize) -> Vec<Vec<f64>> {
    partition_tuples(divisions, m)
        .into_iter()
        .map(|t| t.into_iter().map(|v| v as f64 / divisions as f64).collect())
        .collect()
}

/// Validate and materialize the reference family for `objectives`
/// axes.
///
/// # Errors
/// [`NsgaError::ReferenceInvalid`] when `divisions == 0`,
/// `objectives < 2`, or the family exceeds
/// [`REFERENCE_DIRECTION_CAP`].
#[must_use]
pub fn build_references(divisions: usize, objectives: usize) -> Result<Vec<Vec<f64>>, NsgaError> {
    if divisions == 0 {
        return Err(NsgaError::ReferenceInvalid {
            what: "divisions must be >= 1".to_owned(),
        });
    }
    if objectives < 2 {
        return Err(NsgaError::ReferenceInvalid {
            what: "needs at least two objectives".to_owned(),
        });
    }
    if das_dennis_cardinality(divisions, objectives) > REFERENCE_DIRECTION_CAP {
        return Err(NsgaError::ReferenceInvalid {
            what: format!(
                "reference family exceeds the {}-direction baseline cap",
                REFERENCE_DIRECTION_CAP
            ),
        });
    }
    Ok(das_dennis(divisions, objectives))
}

/// Schema for deterministic adaptive reference-geometry receipts.
pub const NSGA_REFERENCE_GEOMETRY_SCHEMA_VERSION: u32 = 1;

/// A versioned, canonical NSGA-III reference geometry.
///
/// Directions and observed front points live on the non-negative unit simplex.
/// The geometry is immutable: an accepted adaptation returns a new revision, so
/// callers retain the previous value as their exact rollback point.
#[derive(Debug, Clone, PartialEq)]
pub struct NsgaReferenceGeometry {
    /// Receipt schema version.
    pub schema_version: u32,
    /// Number of objective axes.
    pub objectives: usize,
    /// Lexicographically canonical simplex directions.
    pub directions: Vec<Vec<f64>>,
    /// Monotone revision number; fixed geometry starts at zero.
    pub revision: u64,
}

/// Complete, independently comparable identity of a reference geometry.
///
/// This intentionally retains direction bits rather than claiming a short hash
/// is an authority anchor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NsgaReferenceGeometryIdentity {
    /// Receipt schema version.
    pub schema_version: u32,
    /// Number of objective axes.
    pub objectives: usize,
    /// Geometry revision.
    pub revision: u64,
    /// Canonical IEEE-754 direction bits.
    pub direction_bits: Vec<Vec<u64>>,
}

/// Immutable identity of an admitted normalized-front snapshot.
///
/// The complete direction and front bits are retained deliberately: this is
/// evidence for replay and stale-snapshot refusal, not a short-hash claim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NsgaReferenceGeometrySnapshotIdentity {
    /// Schema of the snapshot encoding.
    pub schema_version: u32,
    /// Geometry against which this front was admitted.
    pub geometry: NsgaReferenceGeometryIdentity,
    /// Canonical policy that governed this decision. Distinct caps or
    /// hysteresis thresholds are distinct replay inputs even on one front.
    pub policy: NsgaReferenceGeometryPolicyIdentity,
    /// Sorted, de-duplicated normalized front-point IEEE-754 bits.
    pub front_bits: Vec<Vec<u64>>,
}

/// Canonical identity of the explicit adaptive-geometry policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NsgaReferenceGeometryPolicyIdentity {
    /// Inclusive cap on the direction family.
    pub max_directions: usize,
    /// IEEE-754 bits of the non-negative cover trigger.
    pub cover_trigger_bits: u64,
}

/// Evidence retained only after an admitted adaptation has completed.
///
/// `budget` is the exact ambient fs-exec budget admitted and charged.  Its
/// deadline, work quota, and cancellation result are therefore not advisory
/// metadata. A refusal returns no adaptation and cannot publish this evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NsgaReferenceGeometryEvidence {
    /// Input snapshot identity used for the decision.
    pub snapshot: NsgaReferenceGeometrySnapshotIdentity,
    /// Exact planned and charged work under the admitted context budget.
    pub budget: BudgetConsumption,
}

/// A bounded adaptation policy. It is explicit so a caller cannot silently
/// spend additional reference-direction capacity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NsgaReferenceGeometryPolicy {
    /// Inclusive maximum direction count for all revisions.
    pub max_directions: usize,
    /// Refine only when the observed-front cover radius is strictly larger.
    /// This is the hysteresis boundary.
    pub cover_trigger: f64,
}

impl NsgaReferenceGeometryPolicy {
    /// Complete canonical policy identity retained in every run receipt.
    #[must_use]
    pub fn identity(self) -> NsgaReferenceGeometryPolicyIdentity {
        NsgaReferenceGeometryPolicyIdentity {
            max_directions: self.max_directions,
            cover_trigger_bits: self.cover_trigger.to_bits(),
        }
    }
}

/// One empty interior reference interval bracketed by observed 2D directions.
///
/// This is a sentinel, not a topology certificate: it reports an occupancy gap
/// on the ordered two-objective reference family only.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NsgaDisconnectedFrontSentinel {
    /// Occupied reference immediately before the gap.
    pub left_reference: usize,
    /// First empty reference in the gap.
    pub first_empty_reference: usize,
    /// Occupied reference immediately after the gap.
    pub right_reference: usize,
}

/// Checkable geometry measurements for one normalized front snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct NsgaReferenceGeometryMetrics {
    /// Number of admitted front points.
    pub front_points: usize,
    /// Number of reference directions associated with at least one point.
    pub covered_directions: usize,
    /// Maximum point-to-nearest-reference Euclidean distance, if a front exists.
    pub front_cover_radius: Option<f64>,
    /// Maximum reference-to-nearest-point Euclidean distance, if a front exists.
    pub reference_cover_radius: Option<f64>,
    /// Maximum difference between observed association share and uniform
    /// reference share, if a front exists.
    pub occupancy_discrepancy: Option<f64>,
    /// Ordered 2D occupancy gaps; higher-dimensional fronts deliberately emit
    /// no topology-shaped sentinel.
    pub disconnected_front_sentinels: Vec<NsgaDisconnectedFrontSentinel>,
}

/// The deterministic result of attempting one bounded refinement.
#[derive(Debug, Clone, PartialEq)]
pub struct NsgaReferenceGeometryAdaptation {
    /// Geometry to use after this decision (identical to the input on hold).
    pub geometry: NsgaReferenceGeometry,
    /// Measurements computed from the admitted front snapshot.
    pub metrics: NsgaReferenceGeometryMetrics,
    /// Explicit route selected by the policy.
    pub decision: NsgaReferenceGeometryDecision,
    /// Present only for the context-admitted path. Context-free inspection via
    /// [`NsgaReferenceGeometry::adapt`] intentionally cannot mint budget or
    /// cancellation evidence.
    pub evidence: Option<NsgaReferenceGeometryEvidence>,
}

/// Why an adaptive reference-geometry attempt appended or held.
#[derive(Debug, Clone, PartialEq)]
pub enum NsgaReferenceGeometryDecision {
    /// A canonical observed point was appended as one new direction.
    Appended {
        /// The simplex point appended to the new geometry revision.
        direction: Vec<f64>,
    },
    /// No admitted point requires refinement beyond the hysteresis threshold.
    HeldWithinCoverTrigger,
    /// Every admitted point was already represented exactly.
    HeldNoNovelDirection,
    /// The front snapshot contained no points.
    HeldEmptyFront,
    /// The explicit direction budget prevented refinement.
    HeldDirectionBudget,
}

/// A prepared, context-admitted reference-geometry adaptation.
///
/// This is the resume/fork boundary for geometry work: completing it against
/// the original geometry is deterministic; completing it against a changed
/// geometry fails closed as stale. The object owns normalized input and the
/// admitted executor accountant, so an untrusted caller cannot substitute a
/// different front after admission.
pub struct NsgaReferenceGeometryAdmission<'clock> {
    snapshot: NsgaReferenceGeometrySnapshotIdentity,
    points: Vec<Vec<f64>>,
    policy: NsgaReferenceGeometryPolicy,
    budget: AdmittedBudget<'clock>,
}

impl<'clock> NsgaReferenceGeometryAdmission<'clock> {
    /// The complete snapshot identity bound at admission.
    #[must_use]
    pub fn snapshot(&self) -> &NsgaReferenceGeometrySnapshotIdentity {
        &self.snapshot
    }

    /// Complete a previously admitted operation without partial publication.
    ///
    /// Cancellation, deadline, and work-quota checks run before the result is
    /// materialized and again immediately before returning it. A refusal is an
    /// error, never a held geometry masquerading as an accepted receipt.
    pub fn complete(
        mut self,
        geometry: &NsgaReferenceGeometry,
        cx: &Cx<'_>,
    ) -> Result<NsgaReferenceGeometryAdaptation, NsgaError> {
        if geometry.identity() != self.snapshot.geometry {
            return Err(NsgaError::ReferenceSnapshotStale);
        }
        self.budget
            .checkpoint("nsga-reference-geometry/entry", cx)
            .map_err(reference_admission_refusal)?;
        let planned = self.budget.consumption().planned_cost;
        self.budget
            .charge_cost("nsga-reference-geometry/work", planned)
            .map_err(reference_admission_refusal)?;
        let mut adaptation = geometry.adapt_normalized(&self.points, self.policy)?;
        self.budget
            .checkpoint("nsga-reference-geometry/publication", cx)
            .map_err(reference_admission_refusal)?;
        adaptation.evidence = Some(NsgaReferenceGeometryEvidence {
            snapshot: self.snapshot,
            budget: self.budget.consumption(),
        });
        Ok(adaptation)
    }
}

impl NsgaReferenceGeometry {
    /// Build the canonical fixed Das--Dennis geometry at revision zero.
    pub fn fixed(divisions: usize, objectives: usize) -> Result<Self, NsgaError> {
        Ok(Self {
            schema_version: NSGA_REFERENCE_GEOMETRY_SCHEMA_VERSION,
            objectives,
            directions: build_references(divisions, objectives)?,
            revision: 0,
        })
    }

    /// Return the complete canonical identity used by metrics and route receipts.
    #[must_use]
    pub fn identity(&self) -> NsgaReferenceGeometryIdentity {
        NsgaReferenceGeometryIdentity {
            schema_version: self.schema_version,
            objectives: self.objectives,
            revision: self.revision,
            direction_bits: self
                .directions
                .iter()
                .map(|direction| direction.iter().map(|value| value.to_bits()).collect())
                .collect(),
        }
    }

    /// Measure one normalized front against this geometry.
    pub fn assess(
        &self,
        normalized_front: &[Vec<f64>],
    ) -> Result<NsgaReferenceGeometryMetrics, NsgaError> {
        validate_reference_geometry(self)?;
        let points = normalize_front_points(normalized_front, self.objectives)?;
        Ok(reference_geometry_metrics(
            &self.directions,
            &points,
            self.objectives,
        ))
    }

    /// Adapt by appending at most one canonical observed point under `policy`.
    ///
    /// The input geometry is never mutated. Callers can therefore roll back by
    /// retaining `self`; no hidden evaluation, seed, or comparison work occurs.
    pub fn adapt(
        &self,
        normalized_front: &[Vec<f64>],
        policy: NsgaReferenceGeometryPolicy,
    ) -> Result<NsgaReferenceGeometryAdaptation, NsgaError> {
        validate_reference_geometry(self)?;
        validate_reference_geometry_policy(policy, self.directions.len())?;
        let points = normalize_front_points(normalized_front, self.objectives)?;
        self.adapt_normalized(&points, policy)
    }

    /// Admit a bounded adaptation under the caller's executor context.
    ///
    /// The exact checked work plan is charged to `cx.budget()`; a finite
    /// deadline without an attached deterministic clock refuses at admission.
    /// The returned object can be retained for replay, but only completes when
    /// its geometry identity still matches.
    pub fn prepare_adaptation<'clock>(
        &self,
        normalized_front: &[Vec<f64>],
        policy: NsgaReferenceGeometryPolicy,
        cx: &Cx<'clock>,
    ) -> Result<NsgaReferenceGeometryAdmission<'clock>, NsgaError> {
        validate_reference_geometry(self)?;
        validate_reference_geometry_policy(policy, self.directions.len())?;
        let points = normalize_front_points(normalized_front, self.objectives)?;
        self.prepare_canonical_adaptation(points, policy, cx)
    }

    fn prepare_canonical_adaptation<'clock>(
        &self,
        points: Vec<Vec<f64>>,
        policy: NsgaReferenceGeometryPolicy,
        cx: &Cx<'clock>,
    ) -> Result<NsgaReferenceGeometryAdmission<'clock>, NsgaError> {
        let planned_work = reference_geometry_work_units(points.len(), self.directions.len())?;
        let budget =
            AdmittedBudget::admit_ambient(cx, planned_work).map_err(reference_admission_refusal)?;
        Ok(NsgaReferenceGeometryAdmission {
            snapshot: reference_snapshot_identity(self.identity(), policy, &points),
            points,
            policy,
            budget,
        })
    }

    fn adapt_normalized(
        &self,
        points: &[Vec<f64>],
        policy: NsgaReferenceGeometryPolicy,
    ) -> Result<NsgaReferenceGeometryAdaptation, NsgaError> {
        let metrics = reference_geometry_metrics(&self.directions, &points, self.objectives);
        if points.is_empty() {
            return Ok(NsgaReferenceGeometryAdaptation {
                geometry: self.clone(),
                metrics,
                decision: NsgaReferenceGeometryDecision::HeldEmptyFront,
                evidence: None,
            });
        }
        if metrics
            .front_cover_radius
            .is_some_and(|radius| radius <= policy.cover_trigger)
        {
            return Ok(NsgaReferenceGeometryAdaptation {
                geometry: self.clone(),
                metrics,
                decision: NsgaReferenceGeometryDecision::HeldWithinCoverTrigger,
                evidence: None,
            });
        }
        if self.directions.len() == policy.max_directions {
            return Ok(NsgaReferenceGeometryAdaptation {
                geometry: self.clone(),
                metrics,
                decision: NsgaReferenceGeometryDecision::HeldDirectionBudget,
                evidence: None,
            });
        }

        let candidate = farthest_novel_point(&points, &self.directions);
        let Some(direction) = candidate else {
            return Ok(NsgaReferenceGeometryAdaptation {
                geometry: self.clone(),
                metrics,
                decision: NsgaReferenceGeometryDecision::HeldNoNovelDirection,
                evidence: None,
            });
        };
        let mut directions = self.directions.clone();
        directions.push(direction.clone());
        directions.sort_by(|left, right| vec_ordering(left, right));
        let geometry = Self {
            schema_version: self.schema_version,
            objectives: self.objectives,
            directions,
            revision: self
                .revision
                .checked_add(1)
                .ok_or_else(|| NsgaError::ReferenceInvalid {
                    what: "reference geometry revision overflow".to_owned(),
                })?,
        };
        Ok(NsgaReferenceGeometryAdaptation {
            geometry,
            metrics,
            decision: NsgaReferenceGeometryDecision::Appended { direction },
            evidence: None,
        })
    }
}

fn reference_admission_refusal(refusal: BudgetRefusal) -> NsgaError {
    NsgaError::ReferenceAdmissionRefused {
        what: refusal.to_string(),
    }
}

fn reference_geometry_work_units(points: usize, directions: usize) -> Result<u64, NsgaError> {
    // Association, reverse cover, and novel-point selection each inspect the
    // same bounded point/direction grid. The extra direction-plus-one charge
    // covers materializing and ordering a possible appended direction, even
    // when the decision later holds. Charging this checked upper bound before
    // materialization makes a later refusal unable to publish partial output.
    let grid = u64::try_from(points)
        .ok()
        .and_then(|points| {
            u64::try_from(directions)
                .ok()
                .and_then(|dirs| points.checked_mul(dirs))
        })
        .and_then(|value| value.checked_mul(3))
        .and_then(|grid| {
            u64::try_from(directions)
                .ok()
                .and_then(|directions| directions.checked_add(1))
                .and_then(|append_work| grid.checked_add(append_work))
        })
        .ok_or_else(|| NsgaError::ReferenceInvalid {
            what: "reference geometry work plan overflow".to_owned(),
        })?;
    Ok(grid)
}

fn reference_snapshot_identity(
    geometry: NsgaReferenceGeometryIdentity,
    policy: NsgaReferenceGeometryPolicy,
    points: &[Vec<f64>],
) -> NsgaReferenceGeometrySnapshotIdentity {
    let mut front_bits: Vec<Vec<u64>> = points
        .iter()
        .map(|point| point.iter().map(|value| value.to_bits()).collect())
        .collect();
    front_bits.sort();
    front_bits.dedup();
    NsgaReferenceGeometrySnapshotIdentity {
        schema_version: NSGA_REFERENCE_GEOMETRY_SCHEMA_VERSION,
        geometry,
        policy: policy.identity(),
        front_bits,
    }
}

fn validate_reference_geometry(geometry: &NsgaReferenceGeometry) -> Result<(), NsgaError> {
    if geometry.schema_version != NSGA_REFERENCE_GEOMETRY_SCHEMA_VERSION {
        return Err(NsgaError::ReferenceInvalid {
            what: format!(
                "unsupported reference geometry schema {}",
                geometry.schema_version
            ),
        });
    }
    if geometry.objectives < 2 || geometry.directions.is_empty() {
        return Err(NsgaError::ReferenceInvalid {
            what: "reference geometry needs at least two objectives and one direction".to_owned(),
        });
    }
    if geometry.directions.len() > REFERENCE_DIRECTION_CAP {
        return Err(NsgaError::ReferenceInvalid {
            what: "reference geometry exceeds the baseline direction cap".to_owned(),
        });
    }
    for (index, direction) in geometry.directions.iter().enumerate() {
        if direction.len() != geometry.objectives
            || direction.iter().any(|value| {
                !value.is_finite() || *value < 0.0 || value.to_bits() == (-0.0f64).to_bits()
            })
            || (direction.iter().sum::<f64>() - 1.0).abs() > 32.0 * f64::EPSILON
        {
            return Err(NsgaError::ReferenceInvalid {
                what: format!("direction {index} is not a finite unit-simplex point"),
            });
        }
    }
    for pair in geometry.directions.windows(2) {
        if vec_ordering(&pair[0], &pair[1]) != core::cmp::Ordering::Less {
            return Err(NsgaError::ReferenceInvalid {
                what: "reference geometry directions must be unique and lexicographically ordered"
                    .to_owned(),
            });
        }
    }
    Ok(())
}

fn validate_reference_geometry_policy(
    policy: NsgaReferenceGeometryPolicy,
    current_directions: usize,
) -> Result<(), NsgaError> {
    if policy.max_directions < current_directions
        || policy.max_directions > REFERENCE_DIRECTION_CAP
        || !policy.cover_trigger.is_finite()
        || policy.cover_trigger < 0.0
        || policy.cover_trigger.to_bits() == (-0.0f64).to_bits()
    {
        return Err(NsgaError::ReferenceInvalid {
            what: "adaptive reference policy needs a finite non-negative trigger and a bounded non-shrinking direction cap".to_owned(),
        });
    }
    Ok(())
}

fn normalize_front_points(
    front: &[Vec<f64>],
    objectives: usize,
) -> Result<Vec<Vec<f64>>, NsgaError> {
    let mut normalized: Vec<Vec<f64>> = front
        .iter()
        .enumerate()
        .map(|(index, point)| {
            if point.len() != objectives
                || point.iter().any(|value| {
                    !value.is_finite() || *value < 0.0 || value.to_bits() == (-0.0f64).to_bits()
                })
            {
                return Err(NsgaError::ReferenceInvalid {
                    what: format!(
                        "front point {index} is not finite non-negative objective geometry"
                    ),
                });
            }
            let sum: f64 = point.iter().sum();
            if !(sum > 0.0) || !sum.is_finite() {
                return Err(NsgaError::ReferenceInvalid {
                    what: format!("front point {index} has zero simplex scale"),
                });
            }
            Ok(point
                .iter()
                .map(|value| canonical_positive_zero(value / sum))
                .collect())
        })
        .collect::<Result<_, _>>()?;
    normalized.sort_by(|left, right| vec_ordering(left, right));
    normalized.dedup_by(|left, right| {
        left.iter()
            .zip(right)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    });
    Ok(normalized)
}

fn canonical_positive_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn squared_distance(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(a, b)| {
            let delta = a - b;
            delta * delta
        })
        .sum()
}

fn nearest_reference(point: &[f64], directions: &[Vec<f64>]) -> (usize, f64) {
    let mut best = (0usize, f64::INFINITY);
    for (index, direction) in directions.iter().enumerate() {
        let distance = squared_distance(point, direction);
        if distance < best.1 || (distance == best.1 && index < best.0) {
            best = (index, distance);
        }
    }
    best
}

fn reference_geometry_metrics(
    directions: &[Vec<f64>],
    points: &[Vec<f64>],
    objectives: usize,
) -> NsgaReferenceGeometryMetrics {
    if points.is_empty() {
        return NsgaReferenceGeometryMetrics {
            front_points: 0,
            covered_directions: 0,
            front_cover_radius: None,
            reference_cover_radius: None,
            occupancy_discrepancy: None,
            disconnected_front_sentinels: Vec::new(),
        };
    }
    let mut counts = vec![0usize; directions.len()];
    let mut front_cover_sq = 0.0f64;
    for point in points {
        let (reference, _) = associate_one(point, directions);
        counts[reference] += 1;
        let (_, distance) = nearest_reference(point, directions);
        front_cover_sq = front_cover_sq.max(distance);
    }
    let reference_cover_sq = directions
        .iter()
        .map(|direction| {
            points
                .iter()
                .map(|point| squared_distance(direction, point))
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0f64, f64::max);
    let uniform = 1.0 / directions.len() as f64;
    let discrepancy = counts
        .iter()
        .map(|count| ((*count as f64 / points.len() as f64) - uniform).abs())
        .fold(0.0f64, f64::max);
    let disconnected_front_sentinels = if objectives == 2 {
        disconnected_sentinels(&counts)
    } else {
        Vec::new()
    };
    NsgaReferenceGeometryMetrics {
        front_points: points.len(),
        covered_directions: counts.iter().filter(|&&count| count > 0).count(),
        front_cover_radius: Some(front_cover_sq.sqrt()),
        reference_cover_radius: Some(reference_cover_sq.sqrt()),
        occupancy_discrepancy: Some(discrepancy),
        disconnected_front_sentinels,
    }
}

fn disconnected_sentinels(counts: &[usize]) -> Vec<NsgaDisconnectedFrontSentinel> {
    let mut sentinels = Vec::new();
    let mut left = None;
    let mut index = 0usize;
    while index < counts.len() {
        if counts[index] > 0 {
            left = Some(index);
            index += 1;
            continue;
        }
        let first_empty = index;
        while index < counts.len() && counts[index] == 0 {
            index += 1;
        }
        if let (Some(left_reference), Some(right_reference)) =
            (left, (index < counts.len()).then_some(index))
        {
            sentinels.push(NsgaDisconnectedFrontSentinel {
                left_reference,
                first_empty_reference: first_empty,
                right_reference,
            });
        }
    }
    sentinels
}

fn farthest_novel_point(points: &[Vec<f64>], directions: &[Vec<f64>]) -> Option<Vec<f64>> {
    let mut candidates = points.to_vec();
    candidates.sort_by(|left, right| vec_ordering(left, right));
    candidates.dedup_by(|left, right| {
        left.iter()
            .zip(right)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    });
    candidates
        .into_iter()
        .filter(|point| {
            !directions.iter().any(|direction| {
                point
                    .iter()
                    .zip(direction)
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            })
        })
        .max_by(|left, right| {
            let left_distance = nearest_reference(left, directions).1;
            let right_distance = nearest_reference(right, directions).1;
            left_distance
                .total_cmp(&right_distance)
                .then_with(|| vec_ordering(right, left))
        })
}

/// Outcome of ideal/extreme/intercept normalization.
#[derive(Debug, Clone)]
struct Normalization {
    /// Per-objective translation subtracted before scaling.
    ideal: Vec<f64>,
    /// Per-objective divisors; positive and finite either way.
    scales: Vec<f64>,
    /// Whether the ASF extreme system was singular and the nadir
    /// fallback engaged (surfaced through the receipt).
    singular_fallback: bool,
}

fn compute_normalization(pop: &[NsgaIndividual], m: usize) -> Normalization {
    let mut ideal = vec![f64::INFINITY; m];
    for ind in pop {
        for j in 0..m {
            if ind.f[j] < ideal[j] {
                ideal[j] = ind.f[j];
            }
        }
    }
    let translated: Vec<Vec<f64>> = pop
        .iter()
        .map(|ind| (0..m).map(|j| ind.f[j] - ideal[j]).collect())
        .collect();

    // Axis extremes via the Achievement Scalarizing Function with
    // canonical tie resolution (value equal => lower index wins).
    let eps_non_axis = 1.0e-6f64;
    let mut extremes: Vec<usize> = Vec::with_capacity(m);
    'extremes: for j in 0..m {
        let mut best_val = f64::INFINITY;
        let mut best_idx = usize::MAX;
        for (idx, t) in translated.iter().enumerate() {
            let mut worst = f64::NEG_INFINITY;
            for l in 0..m {
                let v = if l == j { t[l] } else { t[l] / eps_non_axis };
                worst = worst.max(v);
            }
            if worst < best_val || (worst == best_val && idx < best_idx) {
                best_val = worst;
                best_idx = idx;
            }
        }
        if best_idx == usize::MAX || !best_val.is_finite() {
            break 'extremes;
        }
        extremes.push(best_idx);
    }

    let rows_ok = extremes.len() == m
        && extremes
            .iter()
            .all(|&e| e != usize::MAX && translated[e].iter().any(|&v| v.is_finite()))
        && (0..m).all(|j| (j + 1..m).all(|k| translated[extremes[j]] != translated[extremes[k]]));

    let mut singular_fallback = false;
    let mut solved_scales: Option<Vec<f64>> = None;
    if rows_ok {
        let mut at = vec![vec![0.0f64; m]; m];
        for j in 0..m {
            for l in 0..m {
                at[l][j] = translated[extremes[j]][l];
            }
        }
        let col_max: Vec<f64> = (0..m)
            .map(|l| (0..m).map(|r| at[r][l].abs()).fold(0.0f64, f64::max))
            .collect();
        let mut singular = false;
        for col in 0..m {
            let mut piv = col;
            for r in col..m {
                if at[r][col].abs() > at[piv][col].abs() {
                    piv = r;
                }
            }
            if !(at[piv][col].abs() > 1.0e-12 * col_max[col].max(1.0e-300)) {
                singular = true;
                break;
            }
            at.swap(col, piv);
            for r in col + 1..m {
                let factor = at[r][col] / at[col][col];
                // Indexed form required: row `r` is mutated while
                // pivot row `col` is read; iterator adapters would
                // double-borrow (`clippy::needless_range_loop` wants a
                // transform that cannot hold both loans).
                #[allow(clippy::needless_range_loop)]
                for cc in col..m {
                    at[r][cc] -= factor * at[col][cc];
                }
            }
        }
        if singular {
            singular_fallback = true;
        } else {
            let mut z = vec![1.0f64; m];
            for r in (0..m).rev() {
                let mut s = z[r];
                for cc in r + 1..m {
                    s -= at[r][cc] * z[cc];
                }
                z[r] = s / at[r][r];
            }
            let mut scales = vec![0.0f64; m];
            let mut ok = true;
            for (l, zl) in z.iter().enumerate() {
                if !zl.is_finite() || *zl <= 0.0 {
                    ok = false;
                    break;
                }
                scales[l] = 1.0 / zl;
            }
            if ok {
                solved_scales = Some(scales);
            } else {
                singular_fallback = true;
            }
        }
    } else {
        singular_fallback = true;
    }

    let scales = solved_scales.unwrap_or_else(|| {
        singular_fallback = true;
        (0..m)
            .map(|l| {
                let span = translated.iter().map(|t| t[l]).fold(0.0f64, f64::max);
                if span.is_finite() && span > 0.0 {
                    span
                } else {
                    1.0
                }
            })
            .collect()
    });

    Normalization {
        ideal,
        scales,
        singular_fallback,
    }
}

fn perpendicular_distance(t_norm: &[f64], dir: &[f64]) -> f64 {
    let dot: f64 = t_norm.iter().zip(dir.iter()).map(|(t, d)| t * d).sum();
    let direction_norm_sq: f64 = dir.iter().map(|value| value * value).sum();
    debug_assert!(direction_norm_sq.is_finite() && direction_norm_sq > 0.0);
    let projection = dot / direction_norm_sq;
    t_norm
        .iter()
        .zip(dir.iter())
        .map(|(t, d)| {
            let v = t - projection * d;
            v * v
        })
        .sum::<f64>()
        .sqrt()
}

/// Association of one normalized individual to its closest reference
/// direction. Ties resolve to the LOWER reference index, always.
fn associate_one(t_norm: &[f64], dirs: &[Vec<f64>]) -> (usize, f64) {
    let mut best_idx = 0usize;
    let mut best_dist = f64::INFINITY;
    for (j, d) in dirs.iter().enumerate() {
        let dist = perpendicular_distance(t_norm, d);
        if dist < best_dist || (dist == best_dist && j < best_idx) {
            best_dist = dist;
            best_idx = j;
        }
    }
    (best_idx, best_dist)
}

/// Baseline configuration. Box bounds are REQUIRED: the fixture
/// families the matrix targets are box-bounded, and mutation without a
/// box would silently drift units.
#[derive(Debug, Clone)]
pub struct NsgaConfig {
    /// Das–Dennis divisions p >= 1.
    pub reference_divisions: usize,
    /// Fixed per-generation offspring count (population invariant).
    pub population_size: usize,
    /// Hard generation ceiling.
    pub max_generations: usize,
    /// Maximum objective evaluations across the run; integral and
    /// audited. Partial generations never report.
    pub eval_budget: usize,
    /// Seed for the xorshift64* variation stream.
    pub seed: u64,
    /// Inclusive box bounds per decision axis (order matches x).
    pub bounds: Vec<(f64, f64)>,
}

impl NsgaConfig {
    fn validated(self, m_decisions: usize) -> Result<Self, NsgaError> {
        if self.population_size == 0 {
            return Err(NsgaError::PopulationEmpty);
        }
        if self.bounds.len() != m_decisions {
            return Err(NsgaError::BoxBoundsInvalid {
                axis: self.bounds.len(),
            });
        }
        for (i, (lo, hi)) in self.bounds.iter().enumerate() {
            if !(lo <= hi) || !lo.is_finite() || !hi.is_finite() {
                return Err(NsgaError::BoxBoundsInvalid { axis: i });
            }
        }
        if self.eval_budget < self.population_size {
            return Err(NsgaError::BudgetInfeasible {
                minimum_population: self.population_size,
                budget: self.eval_budget,
            });
        }
        Ok(self)
    }
}

/// Why the run stopped after its last completed generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NsgaStop {
    /// Authored generation ceiling reached.
    MaxGenerations,
    /// The next generation could not fit inside the evaluation budget.
    BudgetBoundary,
}

/// Execution receipt.
#[derive(Debug, Clone)]
pub struct NsgaReport {
    /// Generations completed (environmental selection applied).
    pub generations: usize,
    /// Total objective evaluations charged against the budget.
    pub evaluations: usize,
    /// Why the run stopped.
    pub stop: NsgaStop,
    /// Final survivors (deduplicated, frontier-led).
    pub population: Vec<NsgaIndividual>,
    /// Front membership over [`Self::population`].
    pub fronts: Vec<Vec<usize>>,
    /// Whether ANY normalization step fell back to nadir scaling.
    pub normalization_singular_fallback: bool,
    /// Identity of the geometry admitted before the first generation.
    pub reference_geometry_initial: NsgaReferenceGeometryIdentity,
    /// Identity of the geometry used after the final completed generation.
    pub reference_geometry_final: NsgaReferenceGeometryIdentity,
    /// Policy identity retained even when no generation can run.
    pub reference_geometry_policy: NsgaReferenceGeometryPolicyIdentity,
    /// One checkable metric snapshot per completed generation.
    pub reference_geometry_metrics: Vec<NsgaReferenceGeometryMetrics>,
    /// One explicit geometry route decision per completed generation.
    pub reference_geometry_decisions: Vec<NsgaReferenceGeometryDecision>,
    /// Complete normalized-front identities used by each generation's decision.
    pub reference_geometry_snapshots: Vec<NsgaReferenceGeometrySnapshotIdentity>,
    /// Budget/cancellation evidence for context-admitted generations. The
    /// context-free compatibility route leaves this empty rather than claiming
    /// executor authority it did not receive.
    pub reference_geometry_evidence: Vec<NsgaReferenceGeometryEvidence>,
}

/// Deterministic xorshift64* stream: exact bit-sequence spec so future
/// cross-platform batteries reproduce trajectories bitwise.
struct XorShift(u64);

impl XorShift {
    const fn new(seed: u64) -> Self {
        Self(if seed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            seed
        })
    }

    const fn next_u64(&mut self) -> u64 {
        let s = self.0;
        let s1 = s ^ (s >> 12);
        let s2 = s1 ^ (s1 << 25);
        let s3 = s2 ^ (s2 >> 27);
        self.0 = s3;
        s3.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn unit(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) / ((1u128 << 53) as f64)
    }

    fn below(&mut self, modulus: usize) -> usize {
        (self.next_u64() % modulus as u64) as usize
    }
}

fn clamp_axis(v: f64, axis: usize, bounds: &[(f64, f64)]) -> f64 {
    match bounds.get(axis) {
        Some(&(lo, hi)) => v.clamp(lo, hi),
        None => v,
    }
}

fn pick_parent(pop: &[NsgaIndividual], rng: &mut XorShift) -> usize {
    let a = rng.below(pop.len());
    let b = rng.below(pop.len());
    match pop[a].canonical_ordering(&pop[b]) {
        core::cmp::Ordering::Greater => b,
        _ => a,
    }
}

/// Niching fill for the partial last front. Candidate key at decision
/// time: (niche_count, reference_index, distance, original_index);
/// duplicates inside the front were collapsed beforehand, keeping the
/// canonically smallest twin.
fn niche_fill(
    front: &[usize],
    slots: usize,
    assoc: &[(usize, f64)],
    niche_counts: &mut [usize],
) -> Vec<usize> {
    let mut chosen: Vec<usize> = Vec::with_capacity(slots);
    let mut pending: Vec<usize> = front.to_vec();
    while chosen.len() < slots && !pending.is_empty() {
        let mut best_pos = pending[0];
        let mut best_key = (
            niche_counts[assoc[pending[0]].0],
            assoc[pending[0]].0,
            assoc[pending[0]].1,
            pending[0],
        );
        let mut best_slot = 0usize;
        for (slot, &idx) in pending.iter().enumerate().skip(1) {
            let key = (niche_counts[assoc[idx].0], assoc[idx].0, assoc[idx].1, idx);
            if key < best_key {
                best_key = key;
                best_pos = idx;
                best_slot = slot;
            }
        }
        pending.remove(best_slot);
        niche_counts[assoc[best_pos].0] += 1;
        chosen.push(best_pos);
    }
    chosen
}

fn dedup_keep_canonical(idxs: &[usize], union: &[NsgaIndividual]) -> Vec<usize> {
    let mut seen: Vec<&NsgaIndividual> = Vec::new();
    let mut kept: Vec<usize> = Vec::new();
    // Front indices ascend already; twins keep the FIRST occurrence,
    // which is also the canonical (objective-then-decision smaller one
    // only differs when x ordering disagrees with insertion order — we
    // explicitly re-check canonicity to stay enumeration-invariant).
    for &i in idxs {
        let ind = &union[i];
        let dup = seen.iter().any(|s| s.f == ind.f);
        if dup {
            continue;
        }
        seen.push(ind);
        kept.push(i);
    }
    kept
}

/// Run the baseline engine from `initial_pop`. Objective vectors given
/// IN `initial_pop` are consumed as-is; `eval` is invoked only for
/// newly bred decisions, once each, deterministically ordered.
///
/// # Errors
/// Typed [`NsgaError`] on the frozen refusal classes.
pub fn nsga3_run(
    initial_pop: &[NsgaIndividual],
    cfg: &NsgaConfig,
    eval: &mut dyn FnMut(&[f64]) -> Vec<f64>,
) -> Result<NsgaReport, NsgaError> {
    if initial_pop.is_empty() {
        return Err(NsgaError::PopulationEmpty);
    }
    let geometry = NsgaReferenceGeometry::fixed(cfg.reference_divisions, initial_pop[0].f.len())?;
    let policy = NsgaReferenceGeometryPolicy {
        max_directions: geometry.directions.len(),
        cover_trigger: 0.0,
    };
    nsga3_run_with_reference_geometry(initial_pop, cfg, geometry, policy, eval)
}

/// Run NSGA-III with explicit immutable adaptive reference geometry.
///
/// This is the production environmental-selection route: each completed
/// generation receives exactly one policy decision from its normalized union,
/// and the report retains all pre/post identities, metrics, and decisions for
/// replay. The default [`nsga3_run`] supplies a zero-growth policy to preserve
/// its historical fixed Das--Dennis behavior.
pub fn nsga3_run_with_reference_geometry(
    initial_pop: &[NsgaIndividual],
    cfg: &NsgaConfig,
    geometry: NsgaReferenceGeometry,
    geometry_policy: NsgaReferenceGeometryPolicy,
    eval: &mut dyn FnMut(&[f64]) -> Vec<f64>,
) -> Result<NsgaReport, NsgaError> {
    nsga3_run_with_reference_geometry_inner(initial_pop, cfg, geometry, geometry_policy, eval, None)
}

/// Run production NSGA-III with context-admitted adaptive geometry.
///
/// A cancellation, deadline, or work-budget refusal returns no report, so a
/// caller cannot accidentally publish a partial selection as a completed run.
/// Each completed generation retains its complete snapshot and exact fs-exec
/// budget consumption in [`NsgaReport`].
pub fn nsga3_run_with_reference_geometry_cancellable(
    initial_pop: &[NsgaIndividual],
    cfg: &NsgaConfig,
    geometry: NsgaReferenceGeometry,
    geometry_policy: NsgaReferenceGeometryPolicy,
    eval: &mut dyn FnMut(&[f64]) -> Vec<f64>,
    cx: &Cx<'_>,
) -> Result<NsgaReport, NsgaError> {
    nsga3_run_with_reference_geometry_inner(
        initial_pop,
        cfg,
        geometry,
        geometry_policy,
        eval,
        Some(cx),
    )
}

fn nsga3_run_with_reference_geometry_inner(
    initial_pop: &[NsgaIndividual],
    cfg: &NsgaConfig,
    mut geometry: NsgaReferenceGeometry,
    geometry_policy: NsgaReferenceGeometryPolicy,
    eval: &mut dyn FnMut(&[f64]) -> Vec<f64>,
    cx: Option<&Cx<'_>>,
) -> Result<NsgaReport, NsgaError> {
    if initial_pop.is_empty() {
        return Err(NsgaError::PopulationEmpty);
    }
    let m = initial_pop[0].f.len();
    let dx = initial_pop[0].x.len();
    for (i, ind) in initial_pop.iter().enumerate() {
        if ind.f.len() != m {
            return Err(NsgaError::ObjectiveCountMismatch {
                expected: m,
                at: i,
                found: ind.f.len(),
            });
        }
        if ind.x.len() != dx {
            return Err(NsgaError::DecisionCountMismatch {
                expected: dx,
                at: i,
                found: ind.x.len(),
            });
        }
        for (c, v) in ind.f.iter().enumerate() {
            if !v.is_finite() {
                return Err(NsgaError::NonFinite {
                    individual: i,
                    component: c,
                    kind: NonFiniteKind::Objective,
                });
            }
        }
        for (c, v) in ind.x.iter().enumerate() {
            if !v.is_finite() {
                return Err(NsgaError::NonFinite {
                    individual: i,
                    component: c,
                    kind: NonFiniteKind::Decision,
                });
            }
        }
    }
    let cfg = cfg.clone().validated(dx)?;
    validate_reference_geometry(&geometry)?;
    if geometry.objectives != m {
        return Err(NsgaError::ReferenceInvalid {
            what: format!(
                "reference geometry has {} objectives; population has {m}",
                geometry.objectives
            ),
        });
    }
    validate_reference_geometry_policy(geometry_policy, geometry.directions.len())?;
    let reference_geometry_initial = geometry.identity();
    let reference_geometry_policy = geometry_policy.identity();

    let mut rng = XorShift::new(cfg.seed);
    let mut pop: Vec<NsgaIndividual> = initial_pop.to_vec();
    let mut evaluations = pop.len();
    let mut generations_done = 0usize;
    let mut stop = NsgaStop::MaxGenerations;
    let mut normalization_fallback = false;
    let mut reference_geometry_metrics = Vec::new();
    let mut reference_geometry_decisions = Vec::new();
    let mut reference_geometry_snapshots = Vec::new();
    let mut reference_geometry_evidence = Vec::new();

    for _g in 0..cfg.max_generations {
        if cx.is_some_and(Cx::is_cancel_requested) {
            return Err(NsgaError::ReferenceAdmissionRefused {
                what: "cancellation requested before NSGA generation".to_owned(),
            });
        }
        if evaluations + cfg.population_size > cfg.eval_budget {
            stop = NsgaStop::BudgetBoundary;
            break;
        }
        // Variation + evaluation (integral charging per child).
        let mut children: Vec<NsgaIndividual> = Vec::with_capacity(cfg.population_size);
        while children.len() < cfg.population_size {
            let pa = pick_parent(&pop, &mut rng);
            let pb = pick_parent(&pop, &mut rng);
            let t = 0.25 + 0.5 * rng.unit();
            let mut cx = Vec::with_capacity(dx);
            let mut cy = Vec::with_capacity(dx);
            for k in 0..dx {
                let a = pop[pa].x[k];
                let b = pop[pb].x[k];
                let lo = a.min(b);
                let span = (a - b).abs();
                cx.push(clamp_axis(lo - 0.25 * span * t, k, &cfg.bounds));
                cy.push(clamp_axis(lo + span + 0.25 * span * t, k, &cfg.bounds));
            }
            for child_x in [cx, cy] {
                if children.len() == cfg.population_size {
                    break;
                }
                let f = eval(&child_x);
                evaluations += 1;
                if f.len() != m {
                    return Err(NsgaError::ObjectiveCountMismatch {
                        expected: m,
                        at: evaluations,
                        found: f.len(),
                    });
                }
                for (c, v) in f.iter().enumerate() {
                    if !v.is_finite() {
                        return Err(NsgaError::NonFinite {
                            individual: evaluations,
                            component: c,
                            kind: NonFiniteKind::Objective,
                        });
                    }
                }
                children.push(NsgaIndividual { x: child_x, f });
            }
        }
        debug_assert_eq!(children.len(), cfg.population_size);

        // Environmental selection over the union.
        let mut union = pop.clone();
        union.append(&mut children);
        let norm = compute_normalization(&union, m);
        normalization_fallback |= norm.singular_fallback;
        let fronts_union = fast_nondominated_sort(&union);
        let geometry_front: Vec<Vec<f64>> = fronts_union
            .first()
            .into_iter()
            .flatten()
            .filter_map(|&index| {
                let ind = &union[index];
                let point: Vec<f64> = (0..m)
                    .map(|j| canonical_positive_zero((ind.f[j] - norm.ideal[j]) / norm.scales[j]))
                    .collect();
                // The ideal itself has no ray direction. It remains part of
                // environmental selection, but is intentionally absent from
                // reference-geometry adaptation rather than turning a valid
                // run into a zero-simplex-scale refusal.
                (point.iter().sum::<f64>() > 0.0).then_some(point)
            })
            .collect();
        let canonical_geometry_front = normalize_front_points(&geometry_front, m)?;
        let snapshot = reference_snapshot_identity(
            geometry.identity(),
            geometry_policy,
            &canonical_geometry_front,
        );
        let adaptation = if let Some(cx) = cx {
            geometry
                .prepare_canonical_adaptation(canonical_geometry_front, geometry_policy, cx)?
                .complete(&geometry, cx)?
        } else {
            geometry.adapt_normalized(&canonical_geometry_front, geometry_policy)?
        };
        if let Some(evidence) = adaptation.evidence.as_ref()
            && evidence.snapshot != snapshot
        {
            return Err(NsgaError::ReferenceInvalid {
                what: "admitted reference evidence did not bind the production snapshot".to_owned(),
            });
        }
        geometry = adaptation.geometry;
        reference_geometry_metrics.push(adaptation.metrics);
        reference_geometry_decisions.push(adaptation.decision);
        reference_geometry_snapshots.push(snapshot);
        if let Some(evidence) = adaptation.evidence {
            reference_geometry_evidence.push(evidence);
        }
        let dirs = &geometry.directions;
        let assoc: Vec<(usize, f64)> = union
            .iter()
            .map(|ind| {
                let tn: Vec<f64> = (0..m)
                    .map(|j| (ind.f[j] - norm.ideal[j]) / norm.scales[j])
                    .collect();
                associate_one(&tn, &dirs)
            })
            .collect();

        let mut survivors: Vec<usize> = Vec::with_capacity(cfg.population_size);
        let mut niche_counts = vec![0usize; dirs.len()];
        for front in &fronts_union {
            if survivors.len() == cfg.population_size {
                break;
            }
            let slots = cfg.population_size - survivors.len();
            if front.len() <= slots {
                let kept = dedup_keep_canonical(front, &union);
                for idx in &kept {
                    niche_counts[assoc[*idx].0] += 1;
                    survivors.push(*idx);
                    if survivors.len() == cfg.population_size {
                        break;
                    }
                }
            } else {
                let deduped_front = dedup_keep_canonical(front, &union);
                let needed = slots.min(deduped_front.len());
                let picked = niche_fill(&deduped_front, needed, &assoc, &mut niche_counts);
                survivors.extend_from_slice(&picked);
                break;
            }
        }
        // Exact-size guarantee: dedup in the frontier-led phase may
        // legitimately shrink the candidate set below the configured
        // population (degenerate all-twin universes). Refill CANONICALLY
        // from remaining union members first, allowing objective twins
        // whose decisions differ, and only pad with literal clones as
        // the documented last resort. The population never shrinks.
        if survivors.len() < cfg.population_size {
            let mut ordered: Vec<usize> = (0..union.len()).collect();
            ordered.sort_by(|&a, &b| union[a].canonical_ordering(&union[b]).then(a.cmp(&b)));
            for idx in ordered {
                if survivors.len() == cfg.population_size {
                    break;
                }
                let already = survivors.contains(&idx);
                let twin_of_selected = survivors
                    .iter()
                    .any(|&s| union[s].f == union[idx].f && union[s].x == union[idx].x);
                if already || twin_of_selected {
                    continue;
                }
                survivors.push(idx);
            }
            if survivors.len() < cfg.population_size {
                let base_idx: Vec<usize> = survivors.clone();
                let base_len = base_idx.len();
                while survivors.len() < cfg.population_size && base_len > 0 {
                    survivors.push(base_idx[survivors.len() % base_len]);
                }
            }
        }
        debug_assert_eq!(survivors.len(), cfg.population_size);
        generations_done += 1;
        pop = survivors.into_iter().map(|i| union[i].clone()).collect();
    }

    let fronts_final = fast_nondominated_sort(&pop);
    Ok(NsgaReport {
        generations: generations_done,
        evaluations,
        stop,
        population: pop,
        fronts: fronts_final,
        normalization_singular_fallback: normalization_fallback,
        reference_geometry_initial,
        reference_geometry_final: geometry.identity(),
        reference_geometry_policy,
        reference_geometry_metrics,
        reference_geometry_decisions,
        reference_geometry_snapshots,
        reference_geometry_evidence,
    })
}
