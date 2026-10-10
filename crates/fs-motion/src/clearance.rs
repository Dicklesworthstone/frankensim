//! Certified two-sided clearance over prescribed rigid motion.
//!
//! A lower bound alone may prove that bodies remain disjoint, but it does not
//! enclose their minimum separation.  [`separation_over`] pairs a complete
//! time-cell cover (global lower bound) with an admissible fixed-time witness
//! (global upper bound).  Every conversion, spatial, motion-model, and
//! optimization error is an explicit additive length budget.

use crate::{MotionError, SpacetimeChart};
use fs_exec::Cx;
use fs_geom::{Chart, Point3};
use fs_ivl::Interval;

/// Additive uncertainty contributions for one static-clearance claim, in
/// metres.  Each component is a certified nonnegative absolute radius.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ClearanceErrors {
    /// Error introduced by chart/representation conversion.
    pub chart_conversion_m: f64,
    /// Static spatial discretization/query error.
    pub spatial_discretization_m: f64,
    /// Difference between the certified constructed tube and the intended
    /// physical motion model.
    pub motion_tube_m: f64,
    /// Incomplete static/time optimization error.
    pub optimization_m: f64,
}

fn add_nonnegative_upper(left: f64, right: f64) -> f64 {
    let sum = left + right;
    if sum.total_cmp(&0.0).is_eq() {
        0.0
    } else {
        sum.next_up()
    }
}

impl ClearanceErrors {
    fn validate(self) -> Result<Self, MotionError> {
        for (value, what) in [
            (self.chart_conversion_m, "clearance chart-conversion error"),
            (
                self.spatial_discretization_m,
                "clearance spatial-discretization error",
            ),
            (self.motion_tube_m, "clearance motion-tube error"),
            (self.optimization_m, "clearance optimization error"),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(MotionError::InvalidEvidence { what });
            }
        }
        Ok(self)
    }

    /// Outward-rounded total additive error radius.
    pub fn total_upper(self) -> Result<f64, MotionError> {
        let valid = self.validate()?;
        let total = add_nonnegative_upper(
            add_nonnegative_upper(valid.chart_conversion_m, valid.spatial_discretization_m),
            add_nonnegative_upper(valid.motion_tube_m, valid.optimization_m),
        );
        if !total.is_finite() {
            return Err(MotionError::InvalidEvidence {
                what: "clearance total error must be finite",
            });
        }
        Ok(total)
    }

    fn pair(left: Self, right: Self, optimization_m: f64) -> Result<Self, MotionError> {
        let left = left.validate()?;
        let right = right.validate()?;
        let paired = Self {
            chart_conversion_m: add_nonnegative_upper(
                left.chart_conversion_m,
                right.chart_conversion_m,
            ),
            spatial_discretization_m: add_nonnegative_upper(
                left.spatial_discretization_m,
                right.spatial_discretization_m,
            ),
            motion_tube_m: add_nonnegative_upper(left.motion_tube_m, right.motion_tube_m),
            optimization_m: add_nonnegative_upper(
                add_nonnegative_upper(left.optimization_m, right.optimization_m),
                optimization_m,
            ),
        };
        paired.validate()
    }
}

/// Certified static lower-bound evidence over one complete time cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClearanceLowerEvidence {
    /// Lower bound before abstract-region error inflation.
    pub raw_lower_m: f64,
    /// Additive error radius that must be subtracted.
    pub errors: ClearanceErrors,
}

/// A feasible fixed-time configuration supplying an upper bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClearanceWitnessEvidence {
    /// Exact time at which the configuration is feasible.
    pub time: f64,
    /// Upper bound on the static gap at that time before error inflation.
    pub raw_upper_m: f64,
    /// Additive error radius that must be added.
    pub errors: ClearanceErrors,
    /// Stable oracle-specific witness kind for reports and replay.
    pub kind: &'static str,
}

/// Chart-bound provider of rigorous static gap evidence.
///
/// Implementations may use `fs-query` convex separation, exact-distance
/// grids, certified converted meshes, or analytic primitives.  They must bind
/// every returned error term to the two supplied moving charts.
pub trait ClearanceOracle<A: Chart, B: Chart>: Send + Sync {
    /// Lower bound valid for every time in `span`.
    fn lower_bound_over(
        &self,
        a: &SpacetimeChart<A>,
        b: &SpacetimeChart<B>,
        span: Interval,
        cx: &Cx<'_>,
    ) -> Result<ClearanceLowerEvidence, MotionError>;

    /// Optional upper-bound witness at the requested feasible time.  Returning
    /// `None` is honest and forces a lower-only receipt.
    fn witness_at(
        &self,
        a: &SpacetimeChart<A>,
        b: &SpacetimeChart<B>,
        time: f64,
        cx: &Cx<'_>,
    ) -> Result<Option<ClearanceWitnessEvidence>, MotionError>;
}

/// Work and accuracy budget for clearance minimization in time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClearanceConfig {
    /// Requested maximum width of a two-sided range.
    pub value_tolerance_m: f64,
    /// Minimum splittable time-cell width.
    pub time_tolerance: f64,
    /// Maximum binary subdivisions.
    pub max_subdivisions: usize,
}

impl Default for ClearanceConfig {
    fn default() -> Self {
        Self {
            value_tolerance_m: 1.0e-8,
            time_tolerance: 1.0e-8,
            max_subdivisions: 4_096,
        }
    }
}

impl ClearanceConfig {
    fn validate(self) -> Result<Self, MotionError> {
        if !self.value_tolerance_m.is_finite() || self.value_tolerance_m < 0.0 {
            return Err(MotionError::InvalidConfiguration {
                what: "clearance value_tolerance_m must be finite and nonnegative",
            });
        }
        if !self.time_tolerance.is_finite() || self.time_tolerance <= 0.0 {
            return Err(MotionError::InvalidConfiguration {
                what: "clearance time_tolerance must be finite and positive",
            });
        }
        Ok(self)
    }
}

/// Availability of the upper side of a clearance range.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClearanceSidedness {
    /// Both lower and feasible-configuration upper bounds are present.
    TwoSided,
    /// Only a global lower bound is available.
    LowerOnly,
}

/// Whether a two-sided clearance interval met its requested width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClearanceDecision {
    /// The two-sided range meets the requested tolerance.
    Enclosure,
    /// Work/width exhausted, or no upper witness exists.  Every bound still
    /// retains its stated one- or two-sided validity.
    Unknown,
}

/// Error budgets attached to the active lower bound and optional best upper
/// witness.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClearanceRangeErrors {
    /// Error radius subtracted from the winning raw lower bound.
    pub lower: ClearanceErrors,
    /// Error radius added to the best raw upper witness.
    pub upper: Option<ClearanceErrors>,
}

/// Certified minimum-clearance receipt over a time span.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClearanceRange {
    /// Global lower bound on minimum signed surface clearance, in metres.
    pub lower_m: f64,
    /// Feasible-configuration upper bound, when one exists.
    pub upper_m: Option<f64>,
    /// Time of the best feasible upper witness.
    pub witness_time: Option<f64>,
    /// Stable witness kind.
    pub witness_kind: Option<&'static str>,
    /// One- versus two-sided evidence.
    pub sidedness: ClearanceSidedness,
    /// Accuracy result.
    pub decision: ClearanceDecision,
    /// Error accounting for both sides.
    pub errors: ClearanceRangeErrors,
    /// Binary subdivisions performed.
    pub subdivisions: usize,
    /// Complete-cover cell evaluations.
    pub evaluated_cells: usize,
    /// Feasible-time witness attempts.
    pub witness_attempts: usize,
    /// Witnesses actually supplied by the oracle.
    pub admitted_witnesses: usize,
}

impl ClearanceRange {
    /// A strictly positive global lower bound proves separation over the full
    /// requested time span.
    #[must_use]
    pub fn separation_proven(self) -> bool {
        self.lower_m > 0.0
    }
}

#[derive(Debug, Clone, Copy)]
struct ClearanceCell {
    span: Interval,
    lower_m: f64,
    errors: ClearanceErrors,
}

#[derive(Debug, Clone, Copy)]
struct AcceptedWitness {
    upper_m: f64,
    evidence: ClearanceWitnessEvidence,
}

fn split_point(span: Interval) -> Option<f64> {
    let mid = span.midpoint();
    (mid > span.lo() && mid < span.hi()).then_some(mid)
}

fn same_float(left: f64, right: f64) -> bool {
    left.total_cmp(&right).is_eq()
}

fn inflate_lower(
    span: Interval,
    evidence: ClearanceLowerEvidence,
) -> Result<ClearanceCell, MotionError> {
    if !evidence.raw_lower_m.is_finite() {
        return Err(MotionError::InvalidEvidence {
            what: "clearance lower bound must be finite",
        });
    }
    let errors = evidence.errors.validate()?;
    let lower_m = (evidence.raw_lower_m - errors.total_upper()?).next_down();
    if !lower_m.is_finite() {
        return Err(MotionError::InvalidEvidence {
            what: "inflated clearance lower bound must be finite",
        });
    }
    Ok(ClearanceCell {
        span,
        lower_m,
        errors,
    })
}

fn inflate_witness(
    requested_time: f64,
    evidence: ClearanceWitnessEvidence,
) -> Result<AcceptedWitness, MotionError> {
    if !same_float(requested_time, evidence.time) {
        return Err(MotionError::InvalidEvidence {
            what: "clearance witness time does not match the requested feasible time",
        });
    }
    if !evidence.raw_upper_m.is_finite() || evidence.kind.is_empty() {
        return Err(MotionError::InvalidEvidence {
            what: "clearance witness must have a finite upper bound and nonempty kind",
        });
    }
    let errors = evidence.errors.validate()?;
    let upper_m = (evidence.raw_upper_m + errors.total_upper()?).next_up();
    if !upper_m.is_finite() {
        return Err(MotionError::InvalidEvidence {
            what: "inflated clearance upper bound must be finite",
        });
    }
    Ok(AcceptedWitness {
        upper_m,
        evidence: ClearanceWitnessEvidence { errors, ..evidence },
    })
}

fn current_lower(cells: &[ClearanceCell]) -> (f64, ClearanceErrors) {
    let mut best = &cells[0];
    for cell in &cells[1..] {
        if cell.lower_m < best.lower_m
            || (same_float(cell.lower_m, best.lower_m)
                && cell.span.lo().total_cmp(&best.span.lo()).is_lt())
        {
            best = cell;
        }
    }
    (best.lower_m, best.errors)
}

#[allow(clippy::too_many_arguments)]
fn build_receipt(
    cells: &[ClearanceCell],
    witness: Option<AcceptedWitness>,
    decision: ClearanceDecision,
    subdivisions: usize,
    evaluated_cells: usize,
    witness_attempts: usize,
    admitted_witnesses: usize,
) -> Result<ClearanceRange, MotionError> {
    let (lower_m, lower_errors) = current_lower(cells);
    let upper_m = witness.map(|item| item.upper_m);
    if let Some(upper) = upper_m
        && lower_m > upper
    {
        return Err(MotionError::InconsistentEnclosure {
            lower: lower_m,
            upper,
        });
    }
    Ok(ClearanceRange {
        lower_m,
        upper_m,
        witness_time: witness.map(|item| item.evidence.time),
        witness_kind: witness.map(|item| item.evidence.kind),
        sidedness: if witness.is_some() {
            ClearanceSidedness::TwoSided
        } else {
            ClearanceSidedness::LowerOnly
        },
        decision,
        errors: ClearanceRangeErrors {
            lower: lower_errors,
            upper: witness.map(|item| item.evidence.errors),
        },
        subdivisions,
        evaluated_cells,
        witness_attempts,
        admitted_witnesses,
    })
}

/// Certify the minimum signed surface clearance over `span`.
///
/// The oracle lower bounds cover every active time cell.  Witnesses are tried
/// at both endpoints, the root midpoint, and every child midpoint.  No witness
/// is invented when the static provider cannot furnish one.
#[allow(clippy::too_many_lines)] // One complete lower/upper branch-and-bound receipt.
pub fn separation_over<A: Chart, B: Chart, O: ClearanceOracle<A, B>>(
    a: &SpacetimeChart<A>,
    b: &SpacetimeChart<B>,
    span: Interval,
    oracle: &O,
    config: ClearanceConfig,
    cx: &Cx<'_>,
) -> Result<ClearanceRange, MotionError> {
    let config = config.validate()?;
    for domain in [a.tube().domain(), b.tube().domain()] {
        if !domain.encloses(span) {
            return Err(MotionError::OutOfDomain {
                lo: span.lo(),
                hi: span.hi(),
                domain_lo: domain.lo(),
                domain_hi: domain.hi(),
            });
        }
    }

    cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
    let root_evidence = oracle.lower_bound_over(a, b, span, cx)?;
    cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
    let root = inflate_lower(span, root_evidence)?;
    let mut cells = vec![root];
    let mut best_witness: Option<AcceptedWitness> = None;
    let mut witness_attempts = 0usize;
    let mut admitted_witnesses = 0usize;
    for time in [span.lo(), span.midpoint(), span.hi()] {
        cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
        witness_attempts += 1;
        let evidence = oracle.witness_at(a, b, time, cx)?;
        cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
        if let Some(evidence) = evidence {
            admitted_witnesses += 1;
            let candidate = inflate_witness(time, evidence)?;
            if best_witness.is_none_or(|best| {
                candidate.upper_m < best.upper_m
                    || (same_float(candidate.upper_m, best.upper_m)
                        && candidate
                            .evidence
                            .time
                            .total_cmp(&best.evidence.time)
                            .is_lt())
            }) {
                best_witness = Some(candidate);
            }
        }
    }

    let mut subdivisions = 0usize;
    let mut evaluated_cells = 1usize;
    loop {
        cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
        let (lower_m, _) = current_lower(&cells);
        if let Some(witness) = best_witness {
            if lower_m > witness.upper_m {
                return Err(MotionError::InconsistentEnclosure {
                    lower: lower_m,
                    upper: witness.upper_m,
                });
            }
            if witness.upper_m - lower_m <= config.value_tolerance_m {
                return build_receipt(
                    &cells,
                    best_witness,
                    ClearanceDecision::Enclosure,
                    subdivisions,
                    evaluated_cells,
                    witness_attempts,
                    admitted_witnesses,
                );
            }
        }
        if subdivisions >= config.max_subdivisions {
            return build_receipt(
                &cells,
                best_witness,
                ClearanceDecision::Unknown,
                subdivisions,
                evaluated_cells,
                witness_attempts,
                admitted_witnesses,
            );
        }

        let mut selected: Option<usize> = None;
        for (index, cell) in cells.iter().enumerate() {
            if cell.span.width() <= config.time_tolerance || split_point(cell.span).is_none() {
                continue;
            }
            if selected.is_none_or(|current| {
                let incumbent = &cells[current];
                cell.lower_m
                    .total_cmp(&incumbent.lower_m)
                    .then_with(|| cell.span.lo().total_cmp(&incumbent.span.lo()))
                    .then_with(|| cell.span.hi().total_cmp(&incumbent.span.hi()))
                    .is_lt()
            }) {
                selected = Some(index);
            }
        }
        let Some(selected) = selected else {
            return build_receipt(
                &cells,
                best_witness,
                ClearanceDecision::Unknown,
                subdivisions,
                evaluated_cells,
                witness_attempts,
                admitted_witnesses,
            );
        };

        let parent = cells.swap_remove(selected);
        let mid = split_point(parent.span).expect("selected clearance cells are splittable");
        for child_span in [
            Interval::new(parent.span.lo(), mid),
            Interval::new(mid, parent.span.hi()),
        ] {
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            let child_evidence = oracle.lower_bound_over(a, b, child_span, cx)?;
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            let child = inflate_lower(child_span, child_evidence)?;
            cells.push(child);
            evaluated_cells += 1;
            let time = child_span.midpoint();
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            witness_attempts += 1;
            let evidence = oracle.witness_at(a, b, time, cx)?;
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            if let Some(evidence) = evidence {
                admitted_witnesses += 1;
                let candidate = inflate_witness(time, evidence)?;
                if best_witness.is_none_or(|best| {
                    candidate.upper_m < best.upper_m
                        || (same_float(candidate.upper_m, best.upper_m)
                            && candidate
                                .evidence
                                .time
                                .total_cmp(&best.evidence.time)
                                .is_lt())
                }) {
                    best_witness = Some(candidate);
                }
            }
        }
        subdivisions += 1;
    }
}

/// Certified proxy for a spherical body in its own moving frame.
///
/// Nonzero errors bind an approximate chart conversion back to the abstract
/// source body; all-zero errors mean the source itself is this exact sphere.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereClearanceProxy {
    /// Body-frame center.
    pub center: Point3,
    /// Sphere radius in metres.
    pub radius_m: f64,
    /// Representation/motion error budget for this body.
    pub errors: ClearanceErrors,
}

impl SphereClearanceProxy {
    fn validate(self) -> Result<Self, MotionError> {
        if !(self.center.x.is_finite()
            && self.center.y.is_finite()
            && self.center.z.is_finite()
            && self.radius_m.is_finite()
            && self.radius_m > 0.0)
        {
            return Err(MotionError::InvalidGeometry {
                what: "sphere clearance proxy needs a finite center and positive finite radius",
            });
        }
        self.errors.validate()?;
        Ok(self)
    }
}

/// Analytic certified clearance oracle for two spherical proxies.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpherePairClearanceOracle {
    a: SphereClearanceProxy,
    b: SphereClearanceProxy,
    optimization_error_m: f64,
    witnesses_enabled: bool,
}

impl SpherePairClearanceOracle {
    /// Construct a sphere-pair oracle.  Disabling witnesses deliberately
    /// exercises/uses the lower-only evidence lane.
    pub fn new(
        a: SphereClearanceProxy,
        b: SphereClearanceProxy,
        optimization_error_m: f64,
        witnesses_enabled: bool,
    ) -> Result<Self, MotionError> {
        let a = a.validate()?;
        let b = b.validate()?;
        if !optimization_error_m.is_finite() || optimization_error_m < 0.0 {
            return Err(MotionError::InvalidEvidence {
                what: "sphere-pair optimization error must be finite and nonnegative",
            });
        }
        if !(a.radius_m + b.radius_m).is_finite() {
            return Err(MotionError::InvalidGeometry {
                what: "sphere-pair radius sum must be finite",
            });
        }
        Ok(Self {
            a,
            b,
            optimization_error_m,
            witnesses_enabled,
        })
    }

    fn errors(self) -> Result<ClearanceErrors, MotionError> {
        ClearanceErrors::pair(self.a.errors, self.b.errors, self.optimization_error_m)
    }

    fn center_distance(
        &self,
        a: &SpacetimeChart<impl Chart>,
        b: &SpacetimeChart<impl Chart>,
        span: Interval,
        cx: &Cx<'_>,
    ) -> Result<Interval, MotionError> {
        let a_center = a.tube().point_action_over(self.a.center, span, cx)?.coords;
        let b_center = b.tube().point_action_over(self.b.center, span, cx)?.coords;
        let mut norm_squared = Interval::point(0.0);
        for axis in 0..3 {
            let magnitude = (a_center[axis] - b_center[axis]).abs();
            norm_squared = norm_squared + magnitude * magnitude;
        }
        Ok(norm_squared.sqrt())
    }

    fn radius_sum_bounds(self) -> (f64, f64) {
        let sum = self.a.radius_m + self.b.radius_m;
        (sum.next_down(), sum.next_up())
    }
}

impl<A: Chart, B: Chart> ClearanceOracle<A, B> for SpherePairClearanceOracle {
    fn lower_bound_over(
        &self,
        a: &SpacetimeChart<A>,
        b: &SpacetimeChart<B>,
        span: Interval,
        cx: &Cx<'_>,
    ) -> Result<ClearanceLowerEvidence, MotionError> {
        let distance = self.center_distance(a, b, span, cx)?;
        let (_, radius_hi) = self.radius_sum_bounds();
        Ok(ClearanceLowerEvidence {
            raw_lower_m: (distance.lo() - radius_hi).next_down(),
            errors: self.errors()?,
        })
    }

    fn witness_at(
        &self,
        a: &SpacetimeChart<A>,
        b: &SpacetimeChart<B>,
        time: f64,
        cx: &Cx<'_>,
    ) -> Result<Option<ClearanceWitnessEvidence>, MotionError> {
        if !self.witnesses_enabled {
            return Ok(None);
        }
        let distance = self.center_distance(a, b, Interval::point(time), cx)?;
        let (radius_lo, _) = self.radius_sum_bounds();
        Ok(Some(ClearanceWitnessEvidence {
            time,
            raw_upper_m: (distance.hi() - radius_lo).next_up(),
            errors: self.errors()?,
            kind: "sphere-center-distance",
        }))
    }
}

/// A pointwise deepest-common-interior witness for two exact SDFs.
///
/// This is an inradius witness at one point.  It is **not** minimum
/// translation distance, rigid-pose separating displacement, EPA depth, or a
/// contact-response gap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OverlapInradiusWitness {
    /// Common-interior point.
    pub point: Point3,
    /// Certified radius of a ball around `point` inside both exact-SDF bodies.
    pub inradius_lower_m: f64,
}

/// Construct a common-interior witness from rigorous exact-SDF value bands.
pub fn overlap_inradius_witness(
    point: Point3,
    field_a: Interval,
    field_b: Interval,
) -> Result<Option<OverlapInradiusWitness>, MotionError> {
    if !(point.x.is_finite()
        && point.y.is_finite()
        && point.z.is_finite()
        && field_a.lo().is_finite()
        && field_a.hi().is_finite()
        && field_b.lo().is_finite()
        && field_b.hi().is_finite())
    {
        return Err(MotionError::InvalidEvidence {
            what: "overlap witness needs finite point and exact-SDF enclosures",
        });
    }
    let max_hi = field_a.hi().max(field_b.hi());
    if max_hi < 0.0 {
        Ok(Some(OverlapInradiusWitness {
            point,
            inradius_lower_m: (-max_hi).next_down().max(0.0),
        }))
    } else {
        Ok(None)
    }
}
