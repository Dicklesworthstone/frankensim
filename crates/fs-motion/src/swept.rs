//! Certified swept-implicit minimization and fail-closed envelope tracing.
//!
//! For a sign-correct exact-distance base chart moving by `M(t)`, the scalar
//!
//! `psi(x) = inf_t phi(M(t)^-1 x)`
//!
//! has the sign of the swept union.  It is not, in general, the swept
//! region's signed distance.  This module keeps that distinction in the type
//! and evidence surface: [`SweptChart::evaluate`] returns a certified interval
//! for `psi`, while the `Chart` implementation deliberately makes no distance
//! or ray-step claim.

use crate::{MotionError, SpacetimeChart, WankelParams};
use fs_evidence::NumericalCertificate;
use fs_exec::Cx;
use fs_ga::Point as GaPoint;
use fs_geom::{Aabb, BettiBounds, Chart, ChartSample, Differentiability, Point3, TraceStepClaim};
use fs_ivl::Interval;
use fs_math::det;

/// Accuracy and work budget for deterministic swept-field minimization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SweptConfig {
    /// Stop with [`SweepDecision::Enclosure`] when the certified infimum
    /// interval is no wider than this value.
    pub value_tolerance: f64,
    /// Do not subdivide a time cell narrower than this value.
    pub time_tolerance: f64,
    /// Maximum number of binary time-cell subdivisions.  Zero is a valid
    /// inspection-only budget and commonly produces `Unknown`.
    pub max_subdivisions: usize,
}

impl Default for SweptConfig {
    fn default() -> Self {
        Self {
            value_tolerance: 1.0e-8,
            time_tolerance: 1.0e-8,
            max_subdivisions: 4_096,
        }
    }
}

impl SweptConfig {
    fn validate(self) -> Result<Self, MotionError> {
        if !self.value_tolerance.is_finite() || self.value_tolerance < 0.0 {
            return Err(MotionError::InvalidConfiguration {
                what: "swept value_tolerance must be finite and nonnegative",
            });
        }
        if !self.time_tolerance.is_finite() || self.time_tolerance <= 0.0 {
            return Err(MotionError::InvalidConfiguration {
                what: "swept time_tolerance must be finite and positive",
            });
        }
        Ok(self)
    }
}

/// Whether the requested swept-field accuracy was certified within budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SweepDecision {
    /// The returned interval is certified and meets the requested width.
    Enclosure,
    /// The returned interval remains sound, but the requested width was not
    /// proved before the time/work budget ended.
    Unknown,
}

/// Receipt for one deterministic branch-and-bound swept-field query.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SweepReceipt {
    /// Certified enclosure of `inf_t phi(M(t)^-1 x)`.
    pub implicit: Interval,
    /// Accuracy outcome; `Unknown` never weakens `implicit`'s containment.
    pub decision: SweepDecision,
    /// A feasible time whose point enclosure supplied the global upper bound.
    pub witness_time: f64,
    /// Number of binary subdivisions performed.
    pub subdivisions: usize,
    /// Number of time-cell enclosures evaluated.
    pub evaluated_cells: usize,
    /// Number of feasible point witnesses evaluated.
    pub evaluated_witnesses: usize,
    /// Largest motor versor-defect bound observed by the query.
    pub max_defect: f64,
}

#[derive(Debug, Clone, Copy)]
struct SweepCell {
    span: Interval,
    lower: f64,
}

fn split_point(span: Interval) -> Option<f64> {
    let mid = span.midpoint();
    (mid > span.lo() && mid < span.hi()).then_some(mid)
}

fn same_float(left: f64, right: f64) -> bool {
    left.total_cmp(&right).is_eq()
}

/// A chart-like view of a moving region's swept implicit union.
#[derive(Debug)]
pub struct SweptChart<C> {
    moving: SpacetimeChart<C>,
    support: Aabb,
    config: SweptConfig,
}

impl<C: Chart> SweptChart<C> {
    /// Construct a swept chart and certify a finite support enclosure for the
    /// whole motion.  Unbounded base support refuses: blindly transforming an
    /// infinite box would manufacture NaNs and void containment.
    pub fn new(
        moving: SpacetimeChart<C>,
        config: SweptConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, MotionError> {
        let config = config.validate()?;
        let base_support = moving.base().support();
        if !base_support.is_finite() {
            return Err(MotionError::UnboundedSupport);
        }
        let support = moving
            .tube()
            .box_action_over(&base_support, moving.tube().domain(), cx)?
            .bounds;
        if !support.is_finite() {
            return Err(MotionError::UnboundedSupport);
        }
        Ok(Self {
            moving,
            support,
            config,
        })
    }

    /// The underlying time-dependent chart.
    #[must_use]
    pub fn moving(&self) -> &SpacetimeChart<C> {
        &self.moving
    }

    /// The minimization configuration.
    #[must_use]
    pub fn config(&self) -> SweptConfig {
        self.config
    }

    /// Evaluate a certified interval containing the swept implicit infimum.
    ///
    /// Every active time cell contributes a rigorous lower bound.  Point
    /// evaluation at a feasible time contributes a rigorous upper bound on
    /// the global infimum.  Binary subdivision preserves a complete cover of
    /// the time domain, so sampling is not used as proof.
    #[allow(clippy::too_many_lines)] // One auditable global-bound refinement transaction.
    pub fn evaluate(&self, x: Point3, cx: &Cx<'_>) -> Result<SweepReceipt, MotionError> {
        if !(x.x.is_finite() && x.y.is_finite() && x.z.is_finite()) {
            return Err(MotionError::NonFiniteInput { what: "point" });
        }
        let domain = self.moving.tube().domain();
        let root = self.moving.eval_over(x, domain, cx)?;
        let witness_time = domain.midpoint();
        let witness = self
            .moving
            .eval_over(x, Interval::point(witness_time), cx)?;
        let mut cells = vec![SweepCell {
            span: domain,
            lower: root.value.lo(),
        }];
        let mut best_upper = witness.value.hi();
        let mut best_time = witness_time;
        let mut subdivisions = 0usize;
        let mut evaluated_cells = 1usize;
        let mut evaluated_witnesses = 1usize;
        let mut max_defect = root.defect.max(witness.defect);

        loop {
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            let lower = cells
                .iter()
                .map(|cell| cell.lower)
                .fold(f64::INFINITY, f64::min);
            if lower > best_upper {
                return Err(MotionError::InconsistentEnclosure {
                    lower,
                    upper: best_upper,
                });
            }
            let width = best_upper - lower;
            if width <= self.config.value_tolerance {
                return Ok(SweepReceipt {
                    implicit: Interval::new(lower, best_upper),
                    decision: SweepDecision::Enclosure,
                    witness_time: best_time,
                    subdivisions,
                    evaluated_cells,
                    evaluated_witnesses,
                    max_defect,
                });
            }
            if subdivisions >= self.config.max_subdivisions {
                return Ok(SweepReceipt {
                    implicit: Interval::new(lower, best_upper),
                    decision: SweepDecision::Unknown,
                    witness_time: best_time,
                    subdivisions,
                    evaluated_cells,
                    evaluated_witnesses,
                    max_defect,
                });
            }

            // Refine the splittable cell with the smallest lower bound; ties
            // use interval endpoints.  This is fixed-order and scheduler-free.
            let mut selected: Option<usize> = None;
            for (index, cell) in cells.iter().enumerate() {
                if cell.span.width() <= self.config.time_tolerance
                    || split_point(cell.span).is_none()
                {
                    continue;
                }
                let replace = selected.is_none_or(|current| {
                    let incumbent = &cells[current];
                    cell.lower
                        .total_cmp(&incumbent.lower)
                        .then_with(|| cell.span.lo().total_cmp(&incumbent.span.lo()))
                        .then_with(|| cell.span.hi().total_cmp(&incumbent.span.hi()))
                        .is_lt()
                });
                if replace {
                    selected = Some(index);
                }
            }
            let Some(selected) = selected else {
                return Ok(SweepReceipt {
                    implicit: Interval::new(lower, best_upper),
                    decision: SweepDecision::Unknown,
                    witness_time: best_time,
                    subdivisions,
                    evaluated_cells,
                    evaluated_witnesses,
                    max_defect,
                });
            };

            let parent = cells.swap_remove(selected);
            let mid = split_point(parent.span).expect("selected cells are splittable");
            for span in [
                Interval::new(parent.span.lo(), mid),
                Interval::new(mid, parent.span.hi()),
            ] {
                let bound = self.moving.eval_over(x, span, cx)?;
                let time = span.midpoint();
                let point = self.moving.eval_over(x, Interval::point(time), cx)?;
                max_defect = max_defect.max(bound.defect).max(point.defect);
                cells.push(SweepCell {
                    span,
                    lower: bound.value.lo(),
                });
                evaluated_cells += 1;
                evaluated_witnesses += 1;
                let candidate = point.value.hi();
                if candidate < best_upper
                    || (same_float(candidate, best_upper) && time.total_cmp(&best_time).is_lt())
                {
                    best_upper = candidate;
                    best_time = time;
                }
            }
            subdivisions += 1;
        }
    }
}

impl<C: Chart> Chart for SweptChart<C> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        match self.evaluate(x, cx) {
            Ok(receipt) => ChartSample {
                // Representative only.  `error` is no-claim because this is
                // not generally distance from the abstract swept region.
                signed_distance: receipt.implicit.midpoint(),
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::no_claim(),
            },
            Err(_) => ChartSample {
                signed_distance: f64::NAN,
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::no_claim(),
            },
        }
    }

    fn support(&self) -> Aabb {
        self.support
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::NoClaim
    }

    fn topology_hint(&self) -> BettiBounds {
        BettiBounds::unknown()
    }

    fn name(&self) -> &'static str {
        "motion/swept-implicit"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::Unknown
    }

    fn inside(&self, x: Point3, cx: &Cx<'_>) -> bool {
        // A straddling or failed enclosure is not silently rounded into an
        // inside claim.  Callers needing three-valued logic use `evaluate`.
        self.evaluate(x, cx)
            .is_ok_and(|receipt| receipt.implicit.hi() < 0.0)
    }
}

/// Three-valued proof state used by envelope oracles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProofState {
    /// The stated predicate is certified over the full time cell.
    Proven,
    /// Its negation is certified over the full time cell.
    Refuted,
    /// Neither direction is certified.
    Unknown,
}

/// Interval evidence needed before an envelope branch may be admitted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnvelopeEvidence {
    /// Enclosure of `F(x,t)` over the time cell.
    pub field: Interval,
    /// Enclosure of `dF/dt(x,t)` over the time cell.
    pub time_derivative: Interval,
    /// Certified margin for the characteristic-system rank test.  A strictly
    /// positive lower endpoint is required for regularity.
    pub rank_margin: Interval,
    /// Validated existence of a characteristic root in this cell (for
    /// example, interval Newton or implicit-manifold continuation).
    pub characteristic_exists: ProofState,
    /// Whether the branch lies inside the declared parameter/trim domain.
    pub within_trim: ProofState,
    /// Whether the branch is visible on the actual swept boundary.
    pub visible: ProofState,
}

/// Classification of one candidate envelope time cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnvelopeBranchClass {
    /// A regular, visible, in-trim interior characteristic branch.
    RegularInterior,
    /// A separately classified contribution at the path domain endpoint.
    ParameterEndpoint,
    /// Certified not to contain the required characteristic root.
    NotCharacteristic,
    /// Certified rank loss; no regular branch is admitted.
    RankSingular,
    /// Outside the declared trimming domain.
    Trimmed,
    /// Hidden from the actual swept boundary.
    Occluded,
    /// Evidence is incomplete, straddles a degeneracy, or mixes an endpoint
    /// with an interior cell.
    Unknown,
}

/// Classify one envelope branch without upgrading interval overlap into root
/// existence.  A zero-containing interval is necessary, never sufficient.
#[must_use]
pub fn classify_envelope_branch(
    domain: Interval,
    span: Interval,
    evidence: EnvelopeEvidence,
) -> EnvelopeBranchClass {
    if !domain.encloses(span) {
        return EnvelopeBranchClass::NotCharacteristic;
    }
    if !evidence.field.contains_zero() {
        return if evidence.characteristic_exists == ProofState::Proven {
            EnvelopeBranchClass::Unknown
        } else {
            EnvelopeBranchClass::NotCharacteristic
        };
    }
    if evidence.within_trim == ProofState::Refuted {
        return EnvelopeBranchClass::Trimmed;
    }
    if evidence.visible == ProofState::Refuted {
        return EnvelopeBranchClass::Occluded;
    }
    if evidence.characteristic_exists == ProofState::Refuted {
        return EnvelopeBranchClass::NotCharacteristic;
    }
    if evidence.characteristic_exists == ProofState::Unknown
        || evidence.within_trim == ProofState::Unknown
        || evidence.visible == ProofState::Unknown
    {
        return EnvelopeBranchClass::Unknown;
    }

    let point_endpoint = same_float(span.lo(), span.hi())
        && (same_float(span.lo(), domain.lo()) || same_float(span.hi(), domain.hi()));
    if point_endpoint {
        return EnvelopeBranchClass::ParameterEndpoint;
    }
    if !evidence.time_derivative.contains_zero() {
        return EnvelopeBranchClass::Unknown;
    }
    if evidence.rank_margin.hi() <= 0.0 {
        return EnvelopeBranchClass::RankSingular;
    }
    if evidence.rank_margin.lo() <= 0.0 {
        return EnvelopeBranchClass::Unknown;
    }
    if same_float(span.lo(), domain.lo()) || same_float(span.hi(), domain.hi()) {
        return EnvelopeBranchClass::Unknown;
    }
    EnvelopeBranchClass::RegularInterior
}

/// A validated provider for characteristic-set evidence.
///
/// Implementations own the derivative theorem and root-existence proof.  This
/// crate never finite-differences an opaque [`Chart`] and calls the result a
/// certificate.
pub trait EnvelopeOracle<C: Chart>: Send + Sync {
    /// Produce evidence valid over the complete closed time `span`, explicitly
    /// bound to the supplied moving chart and its motor tube.
    fn evidence(
        &self,
        moving: &SpacetimeChart<C>,
        point: Point3,
        span: Interval,
        cx: &Cx<'_>,
    ) -> Result<EnvelopeEvidence, MotionError>;
}

/// Accuracy and work budget for envelope characteristic tracing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnvelopeConfig {
    /// Maximum width of an admitted regular branch's time enclosure.
    pub time_tolerance: f64,
    /// Maximum binary subdivisions of the path domain.
    pub max_subdivisions: usize,
}

impl Default for EnvelopeConfig {
    fn default() -> Self {
        Self {
            time_tolerance: 1.0e-8,
            max_subdivisions: 4_096,
        }
    }
}

impl EnvelopeConfig {
    fn validate(self) -> Result<Self, MotionError> {
        if !self.time_tolerance.is_finite() || self.time_tolerance <= 0.0 {
            return Err(MotionError::InvalidConfiguration {
                what: "envelope time_tolerance must be finite and positive",
            });
        }
        Ok(self)
    }
}

/// One retained envelope branch and the evidence that classified it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnvelopeBranch {
    /// Certified time enclosure for the branch.
    pub span: Interval,
    /// Oracle evidence over `span`.
    pub evidence: EnvelopeEvidence,
    /// Fail-closed branch class.
    pub class: EnvelopeBranchClass,
}

/// Reason-count telemetry for deterministic envelope tracing.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EnvelopeTraceStats {
    /// Binary subdivisions performed.
    pub subdivisions: usize,
    /// Cells passed to the oracle, including exact endpoint cells.
    pub examined_cells: usize,
    /// Cells excluded because no characteristic root can occur.
    pub rejected_not_characteristic: usize,
    /// Cells excluded by the declared trim domain.
    pub rejected_trimmed: usize,
    /// Cells excluded by visibility evidence.
    pub rejected_occluded: usize,
    /// Retained regular interior branches.
    pub regular_branches: usize,
    /// Separately retained endpoint contributions.
    pub endpoint_branches: usize,
    /// Retained cells whose classification remains unresolved.
    pub unresolved_branches: usize,
}

/// Whether every possible characteristic cell was either certified, classified
/// as an endpoint, or rigorously rejected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnvelopeDecision {
    /// No unresolved characteristic branch remains.
    Enclosure,
    /// At least one possible branch remains unresolved.
    Unknown,
}

/// Receipt for an envelope trace at one spatial query point.
#[derive(Debug, Clone, PartialEq)]
pub struct EnvelopeTraceReceipt {
    /// Swept-implicit enclosure proving every admitted envelope point remains
    /// inside the containing swept band.
    pub swept: SweepReceipt,
    /// Regular, endpoint, and unresolved branches; rejected cells are counted
    /// in `stats` rather than retained.
    pub branches: Vec<EnvelopeBranch>,
    /// `Enclosure` only when no unresolved branch remains.
    pub decision: EnvelopeDecision,
    /// Deterministic classification/rejection counters.
    pub stats: EnvelopeTraceStats,
}

/// A swept chart plus a validated characteristic-set oracle.
#[derive(Debug)]
pub struct EnvelopeChart<C, O> {
    swept: SweptChart<C>,
    oracle: O,
    config: EnvelopeConfig,
}

impl<C: Chart, O: EnvelopeOracle<C>> EnvelopeChart<C, O> {
    /// Bind an oracle to a swept chart.
    pub fn new(
        swept: SweptChart<C>,
        oracle: O,
        config: EnvelopeConfig,
    ) -> Result<Self, MotionError> {
        Ok(Self {
            swept,
            oracle,
            config: config.validate()?,
        })
    }

    /// The containing swept chart.
    #[must_use]
    pub fn swept(&self) -> &SweptChart<C> {
        &self.swept
    }

    fn record_rejection(stats: &mut EnvelopeTraceStats, class: EnvelopeBranchClass) -> bool {
        match class {
            EnvelopeBranchClass::NotCharacteristic => stats.rejected_not_characteristic += 1,
            EnvelopeBranchClass::Trimmed => stats.rejected_trimmed += 1,
            EnvelopeBranchClass::Occluded => stats.rejected_occluded += 1,
            _ => return false,
        }
        true
    }

    /// Trace and classify all candidate time cells at `point` under the
    /// oracle's interval exclusion and existence theorem.
    #[allow(clippy::too_many_lines)] // One complete partition/classification receipt.
    pub fn trace(&self, point: Point3, cx: &Cx<'_>) -> Result<EnvelopeTraceReceipt, MotionError> {
        let swept = self.swept.evaluate(point, cx)?;
        let domain = self.swept.moving().tube().domain();
        let mut stats = EnvelopeTraceStats::default();
        let mut branches = Vec::new();

        // Endpoint surfaces are not interior characteristics and therefore do
        // not require dF/dt = 0.  They are evaluated as exact parameter cells.
        for endpoint in [domain.lo(), domain.hi()] {
            let span = Interval::point(endpoint);
            let evidence = self.oracle.evidence(self.swept.moving(), point, span, cx)?;
            stats.examined_cells += 1;
            let class = classify_envelope_branch(domain, span, evidence);
            if Self::record_rejection(&mut stats, class) {
                continue;
            }
            if class == EnvelopeBranchClass::ParameterEndpoint {
                stats.endpoint_branches += 1;
            } else {
                stats.unresolved_branches += 1;
            }
            branches.push(EnvelopeBranch {
                span,
                evidence,
                class,
            });
        }

        let mut pending = vec![domain];
        while !pending.is_empty() {
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            // Widest-first subdivision avoids a depth-first budget bias; ties
            // select the earliest time interval.
            let mut selected = 0usize;
            for index in 1..pending.len() {
                let candidate = pending[index];
                let incumbent = pending[selected];
                if candidate
                    .width()
                    .total_cmp(&incumbent.width())
                    .then_with(|| incumbent.lo().total_cmp(&candidate.lo()))
                    .is_gt()
                {
                    selected = index;
                }
            }
            let span = pending.swap_remove(selected);
            let evidence = self.oracle.evidence(self.swept.moving(), point, span, cx)?;
            stats.examined_cells += 1;
            let class = classify_envelope_branch(domain, span, evidence);
            if Self::record_rejection(&mut stats, class) {
                continue;
            }

            let can_split = span.width() > self.config.time_tolerance
                && stats.subdivisions < self.config.max_subdivisions;
            if can_split {
                if let Some(mid) = split_point(span) {
                    pending.push(Interval::new(span.lo(), mid));
                    pending.push(Interval::new(mid, span.hi()));
                    stats.subdivisions += 1;
                    continue;
                }
            }

            let precise_regular = class == EnvelopeBranchClass::RegularInterior
                && span.width() <= self.config.time_tolerance;
            let retained_class =
                if class == EnvelopeBranchClass::RegularInterior && !precise_regular {
                    EnvelopeBranchClass::Unknown
                } else {
                    class
                };
            if precise_regular {
                stats.regular_branches += 1;
            } else {
                stats.unresolved_branches += 1;
            }
            branches.push(EnvelopeBranch {
                span,
                evidence,
                class: retained_class,
            });
        }

        let decision = if stats.unresolved_branches == 0 {
            EnvelopeDecision::Enclosure
        } else {
            EnvelopeDecision::Unknown
        };
        Ok(EnvelopeTraceReceipt {
            swept,
            branches,
            decision,
            stats,
        })
    }
}

impl<C: Chart, O: EnvelopeOracle<C>> Chart for EnvelopeChart<C, O> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        self.swept.eval(x, cx)
    }

    fn support(&self) -> Aabb {
        self.swept.support()
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::NoClaim
    }

    fn topology_hint(&self) -> BettiBounds {
        BettiBounds::unknown()
    }

    fn name(&self) -> &'static str {
        "motion/envelope-candidate"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::Unknown
    }

    fn inside(&self, x: Point3, cx: &Cx<'_>) -> bool {
        self.swept.inside(x, cx)
    }
}

/// Construct an [`EnvelopeChart`] directly from a moving chart and validated
/// oracle.  The motor path is already content-bound inside `moving`'s tube.
pub fn envelope<C: Chart, O: EnvelopeOracle<C>>(
    moving: SpacetimeChart<C>,
    oracle: O,
    swept_config: SweptConfig,
    envelope_config: EnvelopeConfig,
    cx: &Cx<'_>,
) -> Result<EnvelopeChart<C, O>, MotionError> {
    let swept = SweptChart::new(moving, swept_config, cx)?;
    EnvelopeChart::new(swept, oracle, envelope_config)
}

/// Body-fixed ideal Wankel apex point.
///
/// This is a mathematical point, not a finite apex seal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WankelApexPoint {
    /// Distance from rotor center to the ideal apex point.
    pub body_radius: f64,
    /// Body-frame angular phase of this apex.
    pub body_phase: f64,
}

/// Declared finite circular seal-tip geometry.
///
/// Its center locus is distinct from its contact locus; the actual bore is the
/// envelope of the circle family, not either locus alone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WankelSealCircle {
    /// Seal-center distance from the rotor center.
    pub body_center_radius: f64,
    /// Seal-center phase in the rotor body frame.
    pub body_phase: f64,
    /// Physical seal-tip radius.
    pub tip_radius: f64,
    /// Declared radial clearance added to the contact offset.
    pub clearance: f64,
}

fn validate_wankel_scalar(value: f64, what: &'static str) -> Result<(), MotionError> {
    if !value.is_finite() {
        return Err(MotionError::NonFiniteInput { what });
    }
    Ok(())
}

fn wankel_body_locus(
    params: &WankelParams,
    body_radius: f64,
    body_phase: f64,
    time: f64,
) -> Result<Point3, MotionError> {
    for (value, what) in [
        (params.eccentricity, "wankel eccentricity"),
        (params.omega, "wankel omega"),
        (params.crank_phase, "wankel crank phase"),
        (params.rotor_phase, "wankel rotor phase"),
        (body_radius, "wankel body radius"),
        (body_phase, "wankel body phase"),
        (time, "wankel time"),
    ] {
        validate_wankel_scalar(value, what)?;
    }
    if body_radius < 0.0 {
        return Err(MotionError::InvalidGeometry {
            what: "wankel body radius must be nonnegative",
        });
    }

    // Derivation from the constructed pose, not a mnemonic trochoid:
    // M(t) = base * T(e cos(alpha), e sin(alpha)) * R(beta),
    // alpha = omega*t + crank_phase, beta = alpha/3 + rotor_phase.
    // Acting on r(cos(gamma), sin(gamma)) gives the two vector terms below.
    let alpha = params.omega * time + params.crank_phase;
    let beta = alpha / 3.0 + params.rotor_phase;
    let angle = beta + body_phase;
    if !(alpha.is_finite() && beta.is_finite() && angle.is_finite()) {
        return Err(MotionError::NonFiniteInput {
            what: "derived wankel phase",
        });
    }
    let local = GaPoint {
        x: params.eccentricity * det::cos(alpha) + body_radius * det::cos(angle),
        y: params.eccentricity * det::sin(alpha) + body_radius * det::sin(angle),
        z: 0.0,
    };
    let world = params
        .base_pose
        .transform_point(local)
        .map_err(|_| MotionError::PointActionFailed)?;
    if !(world.x.is_finite() && world.y.is_finite() && world.z.is_finite()) {
        return Err(MotionError::NonFiniteInput {
            what: "wankel base-pose action",
        });
    }
    Ok(Point3::new(world.x, world.y, world.z))
}

impl WankelApexPoint {
    /// Derived ideal-apex epitrochoid point at `time`.
    pub fn at(self, params: &WankelParams, time: f64) -> Result<Point3, MotionError> {
        wankel_body_locus(params, self.body_radius, self.body_phase, time)
    }
}

impl WankelSealCircle {
    fn validate(self) -> Result<Self, MotionError> {
        for (value, what) in [
            (self.body_center_radius, "wankel seal-center radius"),
            (self.body_phase, "wankel seal phase"),
            (self.tip_radius, "wankel seal-tip radius"),
            (self.clearance, "wankel seal clearance"),
        ] {
            validate_wankel_scalar(value, what)?;
        }
        if self.body_center_radius < 0.0 || self.tip_radius < 0.0 || self.clearance < 0.0 {
            return Err(MotionError::InvalidGeometry {
                what: "wankel seal radii and clearance must be nonnegative",
            });
        }
        Ok(self)
    }

    /// Derived seal-center locus at `time`.
    pub fn center_at(self, params: &WankelParams, time: f64) -> Result<Point3, MotionError> {
        let seal = self.validate()?;
        wankel_body_locus(params, seal.body_center_radius, seal.body_phase, time)
    }

    /// One declared seal-contact point for a supplied body-frame unit normal.
    ///
    /// Supplying a normal does not prove that it is the visible bore-envelope
    /// normal; that proof belongs to [`EnvelopeOracle`].
    pub fn contact_at(
        self,
        params: &WankelParams,
        time: f64,
        body_unit_normal: [f64; 2],
    ) -> Result<Point3, MotionError> {
        let seal = self.validate()?;
        let [nx, ny] = body_unit_normal;
        validate_wankel_scalar(nx, "wankel seal normal")?;
        validate_wankel_scalar(ny, "wankel seal normal")?;
        let norm = det::sqrt(nx * nx + ny * ny);
        if (norm - 1.0).abs() > 1.0e-12 {
            return Err(MotionError::InvalidGeometry {
                what: "wankel seal contact normal must be unit length",
            });
        }
        let offset = seal.tip_radius + seal.clearance;
        let center_x = seal.body_center_radius * det::cos(seal.body_phase);
        let center_y = seal.body_center_radius * det::sin(seal.body_phase);
        let radius = det::sqrt((center_x + offset * nx).mul_add(
            center_x + offset * nx,
            (center_y + offset * ny) * (center_y + offset * ny),
        ));
        let phase = det::atan2(center_y + offset * ny, center_x + offset * nx);
        wankel_body_locus(params, radius, phase, time)
    }
}
