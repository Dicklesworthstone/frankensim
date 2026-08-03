//! Validated line/arc meridian charts revolved about the z axis.
//!
//! This module deliberately admits a narrow, closed representation: a simple,
//! counter-clockwise loop in the `rho >= 0` half-plane.  Its revolution is a
//! compact solid and its geometric signed-distance formula is the planar
//! signed distance at `(hypot(x, y), z)`. Construction establishes that real
//! geometry, but v1 has no directed-rounding proof for binary64 evaluation;
//! `Chart` therefore exposes no certified trace or topology claim.

use fs_evidence::NumericalCertificate;
use fs_exec::Cx;
use fs_geom::{
    Aabb, BettiBounds, Chart, ChartSample, Differentiability, Point3, TraceStepClaim, Vec3,
};

/// Hard work bound for both construction and a single exhaustive query.
pub const MAX_AXISYMMETRIC_SEGMENTS: usize = 1024;
const TAU: f64 = core::f64::consts::TAU;
/// Binary64 relative slack.  Individual predicates turn this into a quantity
/// with the same dimension as the quantity being compared.
const ROUNDING_ULPS: f64 = 256.0 * f64::EPSILON;

/// One point in the meridian half-plane. `radius` is cylindrical radius and
/// must be non-negative; `axial` is the coordinate along the revolution axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeridianPoint {
    /// Cylindrical radius.
    pub radius: f64,
    /// Coordinate along the declared z axis.
    pub axial: f64,
}

impl MeridianPoint {
    /// Construct a meridian point. Validation happens at chart construction so
    /// an entire profile is refused transactionally.
    #[must_use]
    pub const fn new(radius: f64, axial: f64) -> Self {
        Self { radius, axial }
    }
}

/// A bounded meridian boundary feature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MeridianSegment {
    /// A straight meridian edge. Constant-radius lines generate cylindrical
    /// faces; sloped lines generate admissible conical chamfers.
    Line {
        /// Oriented start.
        start: MeridianPoint,
        /// Oriented end.
        end: MeridianPoint,
    },
    /// A proper circular arc. These are genuine circular fillets, not a
    /// polygonal approximation. `clockwise` selects the oriented sweep from
    /// `start` to `end` around `center`.
    Arc {
        /// Oriented start. Construction canonicalizes it onto the retained
        /// circle from the center, start angle, and one canonical radius.
        start: MeridianPoint,
        /// Oriented end. Input may differ from the declared circle by the
        /// scale-aware admission residual; construction canonicalizes the
        /// stored loop vertex so all consumers share one circle.
        end: MeridianPoint,
        /// Circle center in the meridian plane.
        center: MeridianPoint,
        /// True for a clockwise sweep; false for counter-clockwise.
        clockwise: bool,
    },
}

/// User-facing outer-edge treatment for [`squat_disc`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SquatDiscEdgeTreatment {
    /// A C0 rim where the planar caps meet the cylindrical outer face.
    Sharp,
    /// Equal, true circular fillets at the upper and lower outer rims.
    /// `radius == 0.0` is exactly the sharp profile; it does not create a
    /// degenerate zero-radius arc.
    CircularFillet {
        /// Radius of each meridian circular fillet.
        radius: f64,
    },
}

/// Authority attached to an [`AxisymmetricSupportPoint`].
///
/// The feature optimization is analytic, but this v1 result is a binary64
/// mechanics input without directed-rounding certification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxisymmetricSupportAuthority {
    /// Deterministic analytic evaluation, not a numerical certificate.
    Estimate,
}

/// A deterministic minimizer of a normalized body-frame support functional.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricSupportPoint {
    /// Unique body-frame point selected by the requested direction.
    pub point: Point3,
    /// `unit_direction dot point`; minimizing the original nonzero direction
    /// yields the same point.
    pub support_value: f64,
    /// Retained meridian feature that supplied the deterministic minimizer.
    pub source_feature: usize,
    /// Explicit non-certificate authority label.
    pub authority: AxisymmetricSupportAuthority,
}

/// Refusal from [`AxisymmetricChart::minimum_support_point`].
#[derive(Debug, Clone, PartialEq)]
pub enum AxisymmetricSupportError {
    /// Input direction has a non-finite component.
    NonFiniteDirection {
        /// Rejected body-frame direction.
        direction: Vec3,
    },
    /// A zero direction has no support ordering.
    ZeroDirection,
    /// Retained construction evidence no longer verifies.
    InvalidChart(AxisymmetricError),
    /// Cancellation was observed before every retained feature was considered.
    Cancelled,
    /// A co-minimizing feature leaves the support point non-unique: it is
    /// flat under the functional or selects a distinct meridian point.
    NonUniqueFeatureSupport {
        /// Source feature that establishes the non-uniqueness.
        source_feature: usize,
    },
    /// An axial direction selected a positive-radius point, whose azimuthal
    /// revolution is a non-unique ring or face support.
    NonUniqueAzimuthalSupport {
        /// Feature supplying the non-unique support.
        source_feature: usize,
    },
    /// A finite input produced a non-finite candidate or support value.
    NonFiniteResult,
}

impl core::fmt::Display for AxisymmetricSupportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFiniteDirection { direction } => write!(
                f,
                "axisymmetric support direction must be finite, got ({}, {}, {})",
                direction.x, direction.y, direction.z
            ),
            Self::ZeroDirection => write!(f, "axisymmetric support direction must be nonzero"),
            Self::InvalidChart(error) => write!(
                f,
                "axisymmetric support requires a still-valid meridian: {error}"
            ),
            Self::Cancelled => write!(
                f,
                "axisymmetric support search was cancelled before publication"
            ),
            Self::NonUniqueFeatureSupport { source_feature } => write!(
                f,
                "axisymmetric support feature {source_feature} leaves no unique support point"
            ),
            Self::NonUniqueAzimuthalSupport { source_feature } => write!(
                f,
                "axisymmetric support feature {source_feature} selects a non-unique axial ring or face"
            ),
            Self::NonFiniteResult => write!(f, "axisymmetric support evaluation became non-finite"),
        }
    }
}

impl std::error::Error for AxisymmetricSupportError {}

impl MeridianSegment {
    fn start(self) -> MeridianPoint {
        match self {
            Self::Line { start, .. } | Self::Arc { start, .. } => start,
        }
    }
    fn end(self) -> MeridianPoint {
        match self {
            Self::Line { end, .. } | Self::Arc { end, .. } => end,
        }
    }
}

/// Stable v1 semantic fingerprint of a validated profile.
///
/// This identifies exact input bits and segment ordering for cache/provenance
/// routing; it is not an admission or scientific-authority token.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AxisymmetricIdentity(pub u64);

/// Construction evidence retained by every admitted chart.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricConstructionCertificate {
    /// Version of the closed input encoding used for [`AxisymmetricIdentity`].
    pub schema_version: u32,
    /// Fingerprint over every semantic profile field.
    pub identity: AxisymmetricIdentity,
    /// Number of input boundary features, including collapsed axis closures.
    pub input_feature_count: usize,
    /// Number of non-degenerate revolved boundary features searched at query.
    pub surfaced_feature_count: usize,
    /// Positive oriented meridian area used to establish the inside rule.
    pub signed_meridian_area: f64,
    /// Whether the validated loop reaches the axis. This is descriptive
    /// construction data only; v1 publishes no topology certificate.
    pub touches_axis: bool,
}

/// Structured refusal from axisymmetric construction.
#[derive(Debug, Clone, PartialEq)]
pub enum AxisymmetricError {
    /// Profile has fewer than three or more than the fixed maximum features.
    SegmentCount {
        /// Observed input feature count.
        count: usize,
    },
    /// A semantic coordinate was not finite.
    NonFinite {
        /// Named semantic input field.
        field: &'static str,
        /// Observed non-finite value.
        value: f64,
    },
    /// A feature left the admissible cylindrical half-plane.
    NegativeRadius {
        /// Observed negative radius.
        value: f64,
    },
    /// A public squat-disc dimension was non-finite or non-positive.
    NonPositiveDimension {
        /// Named public dimension.
        field: &'static str,
        /// Observed invalid value.
        value: f64,
    },
    /// A public squat-disc fillet radius lay outside its closed admissible
    /// interval `[0, min(outer_radius, thickness / 2)]`.
    InvalidEdgeRadius {
        /// Requested fillet radius.
        radius: f64,
        /// Largest admissible radius for this disc.
        maximum: f64,
    },
    /// Consecutive features do not share a literal endpoint.
    OpenLoop {
        /// Index of the feature whose end fails to join its successor.
        index: usize,
    },
    /// A line or arc had zero length, a zero radius, or a nearly-full sweep.
    DegenerateFeature {
        /// Index of the degenerate feature.
        index: usize,
    },
    /// A cancellable validation consumer observed cancellation before all
    /// construction predicates completed.
    Cancelled,
    /// Arc endpoints are not on the declared circle within a conservative
    /// machine-scale construction tolerance.
    InvalidArc {
        /// Index of the arc with inconsistent circle data.
        index: usize,
    },
    /// The arc travels into negative radius, which revolution cannot admit.
    ArcLeavesHalfPlane {
        /// Index of the arc that crosses into negative radius.
        index: usize,
    },
    /// A proper arc has an interior radial minimum tangent to the revolution
    /// axis. Its revolution pinches, so neither of the v1 topology hints is
    /// sound for the resulting singular set.
    AxisTangentArc {
        /// Index of the singular arc.
        index: usize,
    },
    /// The loop does not orient its material to the left (CCW).
    NonPositiveOrientation,
    /// Non-adjacent features meet, overlap, or cross.
    SelfIntersection {
        /// First intersecting feature index.
        first: usize,
        /// Second intersecting feature index.
        second: usize,
    },
}

impl core::fmt::Display for AxisymmetricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SegmentCount { count } => write!(
                f,
                "axisymmetric profile needs 3..={MAX_AXISYMMETRIC_SEGMENTS} segments, got {count}"
            ),
            Self::NonFinite { field, value } => write!(f, "{field} must be finite, got {value}"),
            Self::NegativeRadius { value } => {
                write!(f, "meridian radius must be non-negative, got {value}")
            }
            Self::NonPositiveDimension { field, value } => write!(
                f,
                "squat-disc {field} must be finite and strictly positive, got {value}"
            ),
            Self::InvalidEdgeRadius { radius, maximum } => write!(
                f,
                "squat-disc edge radius must lie in [0, {maximum}], got {radius}"
            ),
            Self::OpenLoop { index } => write!(
                f,
                "feature {index} does not share its endpoint with the next feature"
            ),
            Self::DegenerateFeature { index } => write!(
                f,
                "feature {index} is degenerate or has an unsupported full-circle sweep"
            ),
            Self::Cancelled => write!(
                f,
                "axisymmetric construction validation was cancelled before completion"
            ),
            Self::InvalidArc { index } => write!(
                f,
                "arc {index} endpoints are not equidistant from its center"
            ),
            Self::ArcLeavesHalfPlane { index } => {
                write!(f, "arc {index} enters negative cylindrical radius")
            }
            Self::AxisTangentArc { index } => write!(
                f,
                "arc {index} is tangent to the revolution axis in its interior; singular pinches are unsupported"
            ),
            Self::NonPositiveOrientation => {
                write!(f, "meridian loop must be simple and counter-clockwise")
            }
            Self::SelfIntersection { first, second } => write!(
                f,
                "features {first} and {second} intersect outside their shared join"
            ),
        }
    }
}

impl std::error::Error for AxisymmetricError {}

/// A globally validated line/arc chart of revolution.
#[derive(Debug, Clone)]
pub struct AxisymmetricChart {
    segments: Vec<MeridianSegment>,
    certificate: AxisymmetricConstructionCertificate,
    support: Aabb,
}

impl AxisymmetricChart {
    /// Build a centered compact disc of revolution with an exact sharp or
    /// circularly filleted outer rim.
    ///
    /// `outer_radius` and `thickness` are strictly positive. A circular
    /// `radius` is admitted exactly when
    /// `0 <= radius <= min(outer_radius, thickness / 2)`. The generated
    /// profile uses literal circular arcs and omits collapsed tangent lines at
    /// either equality boundary, so no polygonization or zero-length feature
    /// can enter the generic validator.
    pub fn squat_disc(
        outer_radius: f64,
        thickness: f64,
        edge_treatment: SquatDiscEdgeTreatment,
    ) -> Result<Self, AxisymmetricError> {
        if !outer_radius.is_finite() || outer_radius <= 0.0 {
            return Err(AxisymmetricError::NonPositiveDimension {
                field: "outer_radius",
                value: outer_radius,
            });
        }
        if !thickness.is_finite() || thickness <= 0.0 {
            return Err(AxisymmetricError::NonPositiveDimension {
                field: "thickness",
                value: thickness,
            });
        }
        let half_thickness = thickness * 0.5;
        let edge_radius = match edge_treatment {
            SquatDiscEdgeTreatment::Sharp => 0.0,
            SquatDiscEdgeTreatment::CircularFillet { radius } => {
                let maximum = outer_radius.min(half_thickness);
                if !radius.is_finite() || radius < 0.0 || radius > maximum {
                    return Err(AxisymmetricError::InvalidEdgeRadius { radius, maximum });
                }
                radius
            }
        };
        if edge_radius == 0.0 {
            return Self::try_new(vec![
                line(
                    point(0.0, -half_thickness),
                    point(outer_radius, -half_thickness),
                ),
                line(
                    point(outer_radius, -half_thickness),
                    point(outer_radius, half_thickness),
                ),
                line(
                    point(outer_radius, half_thickness),
                    point(0.0, half_thickness),
                ),
                line(point(0.0, half_thickness), point(0.0, -half_thickness)),
            ]);
        }

        let lower_axis = point(0.0, -half_thickness);
        let lower_tangent = point(outer_radius - edge_radius, -half_thickness);
        let lower_outer = point(outer_radius, -half_thickness + edge_radius);
        let upper_outer = point(outer_radius, half_thickness - edge_radius);
        let upper_tangent = point(outer_radius - edge_radius, half_thickness);
        let upper_axis = point(0.0, half_thickness);
        let mut segments = Vec::with_capacity(5);
        if edge_radius < outer_radius {
            segments.push(line(lower_axis, lower_tangent));
        }
        segments.push(MeridianSegment::Arc {
            start: lower_tangent,
            end: lower_outer,
            center: point(outer_radius - edge_radius, -half_thickness + edge_radius),
            clockwise: false,
        });
        if edge_radius < half_thickness {
            segments.push(line(lower_outer, upper_outer));
        }
        segments.push(MeridianSegment::Arc {
            start: upper_outer,
            end: upper_tangent,
            center: point(outer_radius - edge_radius, half_thickness - edge_radius),
            clockwise: false,
        });
        if edge_radius < outer_radius {
            segments.push(line(upper_tangent, upper_axis));
        }
        segments.push(line(upper_axis, lower_axis));
        Self::try_new(segments)
    }

    /// Validate a closed oriented line/arc loop. Unsupported topology or
    /// malformed geometry is a typed refusal; binary64 trace/topology
    /// authority remains explicitly unavailable in v1.
    pub fn try_new(mut segments: Vec<MeridianSegment>) -> Result<Self, AxisymmetricError> {
        // First validate the caller's literal loop and bounded endpoint
        // residuals. Canonicalizing afterwards cannot repair an open or
        // malformed profile.
        validate_profile(&segments)?;
        canonicalize_arc_endpoints(&mut segments)?;
        let certificate = validate_profile(&segments)?;
        let mut r_max: f64 = 0.0;
        let mut z_min = f64::INFINITY;
        let mut z_max = f64::NEG_INFINITY;
        for segment in &segments {
            for point in segment_extrema(*segment) {
                r_max = r_max.max(point.radius);
                z_min = z_min.min(point.axial);
                z_max = z_max.max(point.axial);
            }
        }
        Ok(Self {
            segments,
            certificate,
            support: Aabb::new(
                Point3::new(-r_max, -r_max, z_min),
                Point3::new(r_max, r_max, z_max),
            ),
        })
    }

    /// Independently re-run all construction obligations against retained
    /// semantic input. This is useful to consumers that deserialize profiles.
    pub fn verify_construction(
        &self,
    ) -> Result<AxisymmetricConstructionCertificate, AxisymmetricError> {
        validate_profile(&self.segments)
    }

    pub(crate) fn verify_construction_with_cx(
        &self,
        cx: &Cx<'_>,
    ) -> Result<AxisymmetricConstructionCertificate, AxisymmetricError> {
        validate_profile_with_cx(&self.segments, Some(cx))
    }

    /// Retained construction certificate.
    #[must_use]
    pub const fn construction_certificate(&self) -> AxisymmetricConstructionCertificate {
        self.certificate
    }

    /// Exact ordered semantic input, for independent reconstruction/checking.
    #[must_use]
    pub fn segments(&self) -> &[MeridianSegment] {
        &self.segments
    }

    /// Analytically minimize the body-frame linear support functional over
    /// every retained line/arc feature.
    ///
    /// The input is normalized internally, so scaling a finite nonzero
    /// direction cannot change the selected point. A nonzero radial component
    /// selects one azimuthal point. Purely axial requests are refused when
    /// their minimizer has positive radius, rather than inventing a rim point
    /// for a non-unique ring or face support.
    pub fn minimum_support_point(
        &self,
        direction: Vec3,
        cx: &Cx<'_>,
    ) -> Result<AxisymmetricSupportPoint, AxisymmetricSupportError> {
        let unit = normalized_support_direction(direction)?;
        match self.verify_construction_with_cx(cx) {
            Ok(_) => {}
            Err(AxisymmetricError::Cancelled) => return Err(AxisymmetricSupportError::Cancelled),
            Err(error) => return Err(AxisymmetricSupportError::InvalidChart(error)),
        }
        let radial = unit.x.hypot(unit.y);
        let radial_coefficient = -radial;
        let mut best: Option<SupportCandidate> = None;
        let mut non_unique_feature = None;
        for (index, segment) in self.segments.iter().copied().enumerate() {
            if index % 16 == 0 && cx.checkpoint().is_err() {
                return Err(AxisymmetricSupportError::Cancelled);
            }
            let candidate = support_candidate(segment, index, radial_coefficient, unit.z);
            if !candidate.radius.is_finite()
                || !candidate.axial.is_finite()
                || !candidate.value.is_finite()
            {
                return Err(AxisymmetricSupportError::NonFiniteResult);
            }
            match best {
                None => {
                    non_unique_feature = if candidate.flat_feature {
                        Some(candidate.source_feature)
                    } else {
                        None
                    };
                    best = Some(candidate);
                }
                Some(current) => {
                    if support_values_tied(candidate.value, current.value) {
                        if candidate.flat_feature || !same_support_candidate(candidate, current) {
                            non_unique_feature = Some(candidate.source_feature);
                        }
                    } else if candidate.value.total_cmp(&current.value).is_lt() {
                        non_unique_feature = if candidate.flat_feature {
                            Some(candidate.source_feature)
                        } else {
                            None
                        };
                        best = Some(candidate);
                    }
                }
            }
        }
        if cx.checkpoint().is_err() {
            return Err(AxisymmetricSupportError::Cancelled);
        }
        let Some(best) = best else {
            return Err(AxisymmetricSupportError::InvalidChart(
                AxisymmetricError::DegenerateFeature { index: 0 },
            ));
        };
        if let Some(source_feature) = non_unique_feature {
            return Err(AxisymmetricSupportError::NonUniqueFeatureSupport { source_feature });
        }
        if radial == 0.0 && best.radius > 0.0 {
            return Err(AxisymmetricSupportError::NonUniqueAzimuthalSupport {
                source_feature: best.source_feature,
            });
        }
        let point = if radial == 0.0 {
            Point3::new(0.0, 0.0, best.axial)
        } else {
            let azimuth_x = unit.x / radial;
            let azimuth_y = unit.y / radial;
            Point3::new(
                -best.radius * azimuth_x,
                -best.radius * azimuth_y,
                best.axial,
            )
        };
        if !finite3(point) {
            return Err(AxisymmetricSupportError::NonFiniteResult);
        }
        Ok(AxisymmetricSupportPoint {
            point,
            support_value: best.value,
            source_feature: best.source_feature,
            authority: AxisymmetricSupportAuthority::Estimate,
        })
    }

    fn query(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        if !finite3(x) {
            return refused_sample();
        }
        let q = MeridianPoint::new(x.x.hypot(x.y), x.z);
        let mut best = f64::INFINITY;
        let mut second = f64::INFINITY;
        let mut closest = None;
        for (index, segment) in self.segments.iter().copied().enumerate() {
            if index % 16 == 0 && cx.checkpoint().is_err() {
                return refused_sample();
            }
            if axis_closure(segment) {
                continue;
            }
            let candidate = nearest_on_segment(segment, q);
            let d2 = squared_delta(q, candidate);
            if d2 < best {
                second = best;
                best = d2;
                closest = Some(candidate);
            } else if d2 < second {
                second = d2;
            }
        }
        let Some(closest) = closest else {
            return refused_sample();
        };
        let distance = best.sqrt();
        let sign = if distance <= query_tolerance(q, closest) {
            0.0
        } else if inside_even_odd(&self.segments, q) {
            -1.0
        } else {
            1.0
        };
        let signed = sign * distance;
        let unique = (second - best).abs() > tie_tolerance(best, second);
        let gradient = if sign == 0.0
            || !unique
            || q.radius <= query_tolerance(q, closest)
            || distance <= query_tolerance(q, closest)
        {
            None
        } else {
            let dr = (q.radius - closest.radius) / distance * sign;
            let dz = (q.axial - closest.axial) / distance * sign;
            Some(Vec3::new(dr * x.x / q.radius, dr * x.y / q.radius, dz))
        };
        ChartSample {
            signed_distance: signed,
            gradient,
            lipschitz: Some(1.0),
            // The closed-form real geometry has no directed-rounding proof
            // for this binary64 evaluation, so it must not expose an exact
            // trace theorem or a one-ULP enclosure.
            error: NumericalCertificate::estimate(signed, signed),
        }
    }
}

/// Build a centered compact disc of revolution. See
/// [`AxisymmetricChart::squat_disc`] for its exact geometry and validation
/// contract.
pub fn squat_disc(
    outer_radius: f64,
    thickness: f64,
    edge_treatment: SquatDiscEdgeTreatment,
) -> Result<AxisymmetricChart, AxisymmetricError> {
    AxisymmetricChart::squat_disc(outer_radius, thickness, edge_treatment)
}

impl Chart for AxisymmetricChart {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        self.query(x, cx)
    }
    fn support(&self) -> Aabb {
        self.support
    }
    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::NoClaim
    }
    fn trace_value_enclosure(
        &self,
        _x: Point3,
        _sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        NumericalCertificate::no_claim()
    }
    fn topology_hint(&self) -> BettiBounds {
        BettiBounds::unknown()
    }
    fn name(&self) -> &'static str {
        "frep/axisymmetric-line-arc"
    }
    fn differentiability(&self) -> Differentiability {
        Differentiability::C0
    }
}

fn refused_sample() -> ChartSample {
    ChartSample {
        signed_distance: f64::NAN,
        gradient: None,
        lipschitz: None,
        error: NumericalCertificate::no_claim(),
    }
}

fn finite3(p: Point3) -> bool {
    p.x.is_finite() && p.y.is_finite() && p.z.is_finite()
}

#[derive(Debug, Clone, Copy)]
struct SupportCandidate {
    radius: f64,
    axial: f64,
    value: f64,
    source_feature: usize,
    flat_feature: bool,
}

fn same_support_candidate(a: SupportCandidate, b: SupportCandidate) -> bool {
    a.radius == b.radius && a.axial == b.axial
}

fn support_values_tied(a: f64, b: f64) -> bool {
    (a - b).abs() <= tie_tolerance(a, b)
}

fn normalized_support_direction(direction: Vec3) -> Result<Vec3, AxisymmetricSupportError> {
    if !direction.x.is_finite() || !direction.y.is_finite() || !direction.z.is_finite() {
        return Err(AxisymmetricSupportError::NonFiniteDirection { direction });
    }
    let scale = direction
        .x
        .abs()
        .max(direction.y.abs())
        .max(direction.z.abs());
    if scale == 0.0 {
        return Err(AxisymmetricSupportError::ZeroDirection);
    }
    let scaled = Vec3::new(
        direction.x / scale,
        direction.y / scale,
        direction.z / scale,
    );
    let norm = scaled.x.hypot(scaled.y).hypot(scaled.z);
    if !norm.is_finite() || norm == 0.0 {
        return Err(AxisymmetricSupportError::NonFiniteDirection { direction });
    }
    Ok(Vec3::new(scaled.x / norm, scaled.y / norm, scaled.z / norm))
}

fn support_candidate(
    segment: MeridianSegment,
    source_feature: usize,
    radial_coefficient: f64,
    axial_coefficient: f64,
) -> SupportCandidate {
    match segment {
        MeridianSegment::Line { start, end } => {
            let slope = radial_coefficient * (end.radius - start.radius)
                + axial_coefficient * (end.axial - start.axial);
            let point = if slope < 0.0 { end } else { start };
            SupportCandidate {
                radius: point.radius,
                axial: point.axial,
                value: radial_coefficient * point.radius + axial_coefficient * point.axial,
                source_feature,
                flat_feature: slope == 0.0,
            }
        }
        MeridianSegment::Arc { start, end, .. } => {
            let mut best = support_point_candidate(
                start,
                source_feature,
                radial_coefficient,
                axial_coefficient,
            );
            let end_candidate =
                support_point_candidate(end, source_feature, radial_coefficient, axial_coefficient);
            if end_candidate.value.total_cmp(&best.value).is_lt() {
                best = end_candidate;
            }
            let interior_angle = (-axial_coefficient).atan2(-radial_coefficient);
            if arc_contains_angle(segment, interior_angle) {
                let interior = support_point_candidate(
                    arc_point_at_angle(segment, interior_angle),
                    source_feature,
                    radial_coefficient,
                    axial_coefficient,
                );
                if interior.value.total_cmp(&best.value).is_lt() {
                    best = interior;
                }
            }
            best
        }
    }
}

fn support_point_candidate(
    point: MeridianPoint,
    source_feature: usize,
    radial_coefficient: f64,
    axial_coefficient: f64,
) -> SupportCandidate {
    SupportCandidate {
        radius: point.radius,
        axial: point.axial,
        value: radial_coefficient * point.radius + axial_coefficient * point.axial,
        source_feature,
        flat_feature: false,
    }
}
fn point(radius: f64, axial: f64) -> MeridianPoint {
    MeridianPoint::new(radius, axial)
}
fn line(start: MeridianPoint, end: MeridianPoint) -> MeridianSegment {
    MeridianSegment::Line { start, end }
}
fn finite_point(p: MeridianPoint) -> bool {
    p.radius.is_finite() && p.axial.is_finite()
}
fn squared_delta(a: MeridianPoint, b: MeridianPoint) -> f64 {
    let dr = a.radius - b.radius;
    let dz = a.axial - b.axial;
    dr.mul_add(dr, dz * dz)
}
fn scale_of(points: &[MeridianPoint]) -> f64 {
    // Predicate slack must follow local geometric differences, not the world
    // axial origin. A large rigid z translation cannot make a small fillet
    // look less simple or less circular.
    let mut scale = f64::MIN_POSITIVE;
    for (index, point) in points.iter().copied().enumerate() {
        for other in points.iter().copied().skip(index + 1) {
            scale = scale
                .max((point.radius - other.radius).abs())
                .max((point.axial - other.axial).abs());
        }
    }
    scale
}
fn length_tolerance(points: &[MeridianPoint]) -> f64 {
    (ROUNDING_ULPS * scale_of(points)).max(f64::MIN_POSITIVE)
}
fn query_tolerance(a: MeridianPoint, b: MeridianPoint) -> f64 {
    length_tolerance(&[a, b])
}
fn area_tolerance(points: &[MeridianPoint]) -> f64 {
    let scale = scale_of(points);
    ROUNDING_ULPS * scale * scale
}
fn angle_tolerance() -> f64 {
    ROUNDING_ULPS
}
fn discriminant_tolerance(aa: f64, bb: f64, cc: f64) -> f64 {
    // `b^2 - 4ac` has units L^4.  Scaling by the competing L^4 terms,
    // rather than a length tolerance, keeps tangent classification invariant
    // under a uniform unit change.
    ROUNDING_ULPS
        * (bb * bb)
            .abs()
            .max((4.0 * aa * cc).abs())
            .max(f64::MIN_POSITIVE)
}
fn tie_tolerance(a: f64, b: f64) -> f64 {
    (ROUNDING_ULPS * a.abs().max(b.abs()).max(f64::MIN_POSITIVE)).max(f64::MIN_POSITIVE)
}
fn same_point(a: MeridianPoint, b: MeridianPoint) -> bool {
    a.radius == b.radius && a.axial == b.axial
}

/// Replace each admitted arc endpoint with a point on one retained canonical
/// circle, then propagate that point to the adjoining line/arc vertex. The
/// public input still has to close literally before this step; this only
/// removes sub-ULP subtractive disagreement between a fillet's tangent paths
/// and the circle consumed by distance, support, and Green integrals.
fn canonicalize_arc_endpoints(segments: &mut [MeridianSegment]) -> Result<(), AxisymmetricError> {
    let mut vertices: Vec<_> = segments.iter().map(|segment| segment.start()).collect();
    let mut canonical_vertices: Vec<Option<(MeridianPoint, usize)>> = vec![None; segments.len()];

    for (index, segment) in segments.iter().copied().enumerate() {
        let MeridianSegment::Arc {
            start, end, center, ..
        } = segment
        else {
            continue;
        };
        let radius = squared_delta(start, center).sqrt();
        let start_angle = (start.axial - center.axial).atan2(start.radius - center.radius);
        let end_angle = (end.axial - center.axial).atan2(end.radius - center.radius);
        let start = canonical_circle_point(center, radius, start_angle);
        let end = canonical_circle_point(center, radius, end_angle);
        for (vertex, point) in [(index, start), ((index + 1) % segments.len(), end)] {
            match canonical_vertices[vertex] {
                None => canonical_vertices[vertex] = Some((point, index)),
                Some((prior, prior_index)) => {
                    if !nearby_point(prior, point) {
                        return Err(AxisymmetricError::InvalidArc {
                            index: prior_index.min(index),
                        });
                    }
                }
            }
        }
    }

    for (index, replacement) in canonical_vertices.into_iter().enumerate() {
        if let Some((point, _)) = replacement {
            vertices[index] = point;
        }
    }
    for (index, segment) in segments.iter_mut().enumerate() {
        let start = vertices[index];
        let end = vertices[(index + 1) % vertices.len()];
        match segment {
            MeridianSegment::Line {
                start: line_start,
                end: line_end,
            }
            | MeridianSegment::Arc {
                start: line_start,
                end: line_end,
                ..
            } => {
                *line_start = start;
                *line_end = end;
            }
        }
    }
    Ok(())
}

fn canonical_circle_point(center: MeridianPoint, radius: f64, angle: f64) -> MeridianPoint {
    MeridianPoint::new(
        center.radius + radius * angle.cos(),
        center.axial + radius * angle.sin(),
    )
}

fn validate_profile(
    segments: &[MeridianSegment],
) -> Result<AxisymmetricConstructionCertificate, AxisymmetricError> {
    validate_profile_with_cx(segments, None)
}

fn validate_profile_with_cx(
    segments: &[MeridianSegment],
    cx: Option<&Cx<'_>>,
) -> Result<AxisymmetricConstructionCertificate, AxisymmetricError> {
    if !(3..=MAX_AXISYMMETRIC_SEGMENTS).contains(&segments.len()) {
        return Err(AxisymmetricError::SegmentCount {
            count: segments.len(),
        });
    }
    for (index, segment) in segments.iter().copied().enumerate() {
        if index % 16 == 0 && validation_checkpoint(cx)? {
            return Err(AxisymmetricError::Cancelled);
        }
        let points = match segment {
            MeridianSegment::Line { start, end } => [start, end, start],
            MeridianSegment::Arc {
                start, end, center, ..
            } => [start, end, center],
        };
        for point in points {
            if !finite_point(point) {
                return Err(AxisymmetricError::NonFinite {
                    field: "meridian coordinate",
                    value: if !point.radius.is_finite() {
                        point.radius
                    } else {
                        point.axial
                    },
                });
            }
            if point.radius < 0.0 {
                return Err(AxisymmetricError::NegativeRadius {
                    value: point.radius,
                });
            }
        }
        if same_point(segment.start(), segment.end()) {
            return Err(AxisymmetricError::DegenerateFeature { index });
        }
        if matches!(segment, MeridianSegment::Arc { .. }) {
            validate_arc(index, segment)?;
        }
    }
    for index in 0..segments.len() {
        if !same_point(
            segments[index].end(),
            segments[(index + 1) % segments.len()].start(),
        ) {
            return Err(AxisymmetricError::OpenLoop { index });
        }
    }
    let area = signed_area(segments);
    let points: Vec<_> = segments.iter().map(|s| s.start()).collect();
    if !area.is_finite() || area <= area_tolerance(&points) {
        return Err(AxisymmetricError::NonPositiveOrientation);
    }
    for first in 0..segments.len() {
        for second in first + 1..segments.len() {
            if (first * segments.len() + second) % 16 == 0 && validation_checkpoint(cx)? {
                return Err(AxisymmetricError::Cancelled);
            }
            if adjacent(first, second, segments.len()) {
                let join = shared_adjacent_join(segments, first, second);
                if adjacent_features_reintersect(segments[first], segments[second], join) {
                    return Err(AxisymmetricError::SelfIntersection { first, second });
                }
                continue;
            }
            if segments_intersect(segments[first], segments[second]) {
                return Err(AxisymmetricError::SelfIntersection { first, second });
            }
        }
    }
    let surfaced_feature_count = segments
        .iter()
        .copied()
        .filter(|s| !axis_closure(*s))
        .count();
    if surfaced_feature_count == 0 {
        return Err(AxisymmetricError::DegenerateFeature { index: 0 });
    }
    let touches_axis = segments
        .iter()
        .any(|s| s.start().radius == 0.0 || s.end().radius == 0.0);
    Ok(AxisymmetricConstructionCertificate {
        schema_version: 1,
        identity: profile_identity(segments),
        input_feature_count: segments.len(),
        surfaced_feature_count,
        signed_meridian_area: area,
        touches_axis,
    })
}

fn validation_checkpoint(cx: Option<&Cx<'_>>) -> Result<bool, AxisymmetricError> {
    Ok(cx.is_some_and(|context| context.checkpoint().is_err()))
}

fn validate_arc(index: usize, segment: MeridianSegment) -> Result<(), AxisymmetricError> {
    let MeridianSegment::Arc {
        start, end, center, ..
    } = segment
    else {
        unreachable!();
    };
    let rs = squared_delta(start, center).sqrt();
    let re = squared_delta(end, center).sqrt();
    let tolerance = length_tolerance(&[start, end, center]);
    if rs <= tolerance || (rs - re).abs() > tolerance {
        return Err(AxisymmetricError::InvalidArc { index });
    }
    if arc_sweep(segment).abs() <= angle_tolerance()
        || TAU - arc_sweep(segment).abs() <= angle_tolerance()
    {
        return Err(AxisymmetricError::DegenerateFeature { index });
    }
    let leftmost = if arc_contains_angle(segment, core::f64::consts::PI) {
        center.radius - rs
    } else {
        start.radius.min(end.radius)
    };
    if leftmost < -tolerance {
        return Err(AxisymmetricError::ArcLeavesHalfPlane { index });
    }
    if arc_has_interior_axis_tangent(segment, tolerance) {
        return Err(AxisymmetricError::AxisTangentArc { index });
    }
    Ok(())
}

fn arc_has_interior_axis_tangent(segment: MeridianSegment, tolerance: f64) -> bool {
    let MeridianSegment::Arc { center, .. } = segment else {
        return false;
    };
    let radius = arc_radius(segment);
    if (center.radius - radius).abs() > tolerance
        || !arc_contains_angle(segment, core::f64::consts::PI)
    {
        return false;
    }
    let sweep = arc_sweep(segment);
    let start = arc_start_angle(segment);
    let travel = if sweep >= 0.0 {
        (core::f64::consts::PI - start).rem_euclid(TAU)
    } else {
        (start - core::f64::consts::PI).rem_euclid(TAU)
    };
    travel > angle_tolerance() && travel < sweep.abs() - angle_tolerance()
}

fn arc_radius(segment: MeridianSegment) -> f64 {
    match segment {
        MeridianSegment::Arc { start, center, .. } => squared_delta(start, center).sqrt(),
        _ => 0.0,
    }
}
fn arc_start_angle(segment: MeridianSegment) -> f64 {
    match segment {
        MeridianSegment::Arc { start, center, .. } => {
            (start.axial - center.axial).atan2(start.radius - center.radius)
        }
        _ => 0.0,
    }
}
fn arc_sweep(segment: MeridianSegment) -> f64 {
    match segment {
        MeridianSegment::Arc {
            start,
            end,
            center,
            clockwise,
        } => {
            let a = (start.axial - center.axial).atan2(start.radius - center.radius);
            let b = (end.axial - center.axial).atan2(end.radius - center.radius);
            let d = if clockwise {
                (a - b).rem_euclid(TAU)
            } else {
                (b - a).rem_euclid(TAU)
            };
            if clockwise { -d } else { d }
        }
        _ => 0.0,
    }
}
fn arc_contains_angle(segment: MeridianSegment, angle: f64) -> bool {
    let sweep = arc_sweep(segment);
    let start = arc_start_angle(segment);
    let delta = if sweep >= 0.0 {
        (angle - start).rem_euclid(TAU)
    } else {
        (start - angle).rem_euclid(TAU)
    };
    delta <= sweep.abs() + angle_tolerance()
}
fn segment_extrema(segment: MeridianSegment) -> Vec<MeridianPoint> {
    let mut points = vec![segment.start(), segment.end()];
    if matches!(segment, MeridianSegment::Arc { .. }) {
        for angle in [
            0.0,
            core::f64::consts::FRAC_PI_2,
            core::f64::consts::PI,
            3.0 * core::f64::consts::FRAC_PI_2,
        ] {
            if arc_contains_angle(segment, angle) {
                points.push(arc_point_at_angle(segment, angle));
            }
        }
    }
    points
}
fn arc_point_at_angle(segment: MeridianSegment, angle: f64) -> MeridianPoint {
    let MeridianSegment::Arc { center, .. } = segment else {
        unreachable!()
    };
    let radius = arc_radius(segment);
    MeridianPoint::new(
        center.radius + radius * angle.cos(),
        center.axial + radius * angle.sin(),
    )
}
fn axis_closure(segment: MeridianSegment) -> bool {
    matches!(segment, MeridianSegment::Line { start, end } if start.radius == 0.0 && end.radius == 0.0)
}

fn nearest_on_segment(segment: MeridianSegment, q: MeridianPoint) -> MeridianPoint {
    match segment {
        MeridianSegment::Line { start, end } => {
            let dr = end.radius - start.radius;
            let dz = end.axial - start.axial;
            let denom = dr.mul_add(dr, dz * dz);
            let u = (((q.radius - start.radius) * dr + (q.axial - start.axial) * dz) / denom)
                .clamp(0.0, 1.0);
            MeridianPoint::new(start.radius + u * dr, start.axial + u * dz)
        }
        MeridianSegment::Arc { center, .. } => {
            let vr = q.radius - center.radius;
            let vz = q.axial - center.axial;
            if vr == 0.0 && vz == 0.0 {
                return segment.start();
            }
            let angle = vz.atan2(vr);
            if arc_contains_angle(segment, angle) {
                arc_point_at_angle(segment, angle)
            } else {
                let start = segment.start();
                let end = segment.end();
                if squared_delta(q, start) <= squared_delta(q, end) {
                    start
                } else {
                    end
                }
            }
        }
    }
}

fn inside_even_odd(segments: &[MeridianSegment], q: MeridianPoint) -> bool {
    let mut crossings = 0_u32;
    for segment in segments.iter().copied() {
        match segment {
            MeridianSegment::Line { start, end } => {
                if (start.axial <= q.axial && q.axial < end.axial)
                    || (end.axial <= q.axial && q.axial < start.axial)
                {
                    let u = (q.axial - start.axial) / (end.axial - start.axial);
                    if start.radius + u * (end.radius - start.radius) > q.radius {
                        crossings += 1;
                    }
                }
            }
            MeridianSegment::Arc { center, .. } => {
                let radius = arc_radius(segment);
                let v = (q.axial - center.axial) / radius;
                if !(-1.0..=1.0).contains(&v) {
                    continue;
                }
                let a = v.asin();
                let angles = if a.abs() == core::f64::consts::FRAC_PI_2 {
                    [Some(a), None]
                } else {
                    [Some(a), Some(core::f64::consts::PI - a)]
                };
                for angle in angles.into_iter().flatten() {
                    if !arc_contains_angle(segment, angle) {
                        continue;
                    }
                    if !arc_half_open_crossing(segment, angle) {
                        continue;
                    }
                    if center.radius + radius * angle.cos() > q.radius {
                        crossings += 1;
                    }
                }
            }
        }
    }
    crossings % 2 == 1
}

/// Match line crossings' lower-inclusive/upper-exclusive axial convention at
/// every line/arc join. Endpoint inclusion follows the local axial derivative
/// under the oriented sweep, not the endpoints' global axial ordering; major
/// arcs can reverse direction between their endpoints. Interior horizontal
/// extrema are tangent contacts and contribute no parity crossing.
fn arc_half_open_crossing(segment: MeridianSegment, angle: f64) -> bool {
    let sweep = arc_sweep(segment);
    let start_angle = arc_start_angle(segment);
    let travel = if sweep >= 0.0 {
        (angle - start_angle).rem_euclid(TAU)
    } else {
        (start_angle - angle).rem_euclid(TAU)
    };
    let end_distance = sweep.abs() - travel;
    let axial_derivative = sweep.signum() * angle.cos();
    if travel <= angle_tolerance() {
        return axial_derivative > angle_tolerance();
    }
    if end_distance <= angle_tolerance() {
        return axial_derivative < -angle_tolerance();
    }
    angle.cos().abs() > angle_tolerance()
}

fn signed_area(segments: &[MeridianSegment]) -> f64 {
    let axial_origin = segments[0].start().axial;
    0.5 * segments
        .iter()
        .copied()
        .map(|segment| match segment {
            MeridianSegment::Line { start, end } => {
                start.radius * (end.axial - axial_origin)
                    - (start.axial - axial_origin) * end.radius
            }
            MeridianSegment::Arc { center, .. } => {
                let r = arc_radius(segment);
                let a = arc_start_angle(segment);
                let b = a + arc_sweep(segment);
                r * center.radius * (b.sin() - a.sin())
                    - r * (center.axial - axial_origin) * (b.cos() - a.cos())
                    + r * r * (b - a)
            }
        })
        .sum::<f64>()
}

fn adjacent(first: usize, second: usize, len: usize) -> bool {
    second == first + 1 || (first == 0 && second + 1 == len)
}

fn shared_adjacent_join(
    segments: &[MeridianSegment],
    first: usize,
    second: usize,
) -> MeridianPoint {
    if second == first + 1 {
        segments[first].end()
    } else {
        segments[first].start()
    }
}

fn nearby_point(a: MeridianPoint, b: MeridianPoint) -> bool {
    squared_delta(a, b) <= query_tolerance(a, b).powi(2)
}

/// Compare a numerically reconstructed intersection with the literal shared
/// join using the adjacent features' local extent.  A tolerance derived from
/// just the two almost-identical points collapses to zero after a rigid
/// translation and would falsely turn a literal line/arc join into a crossing.
fn nearby_feature_point(
    point: MeridianPoint,
    join: MeridianPoint,
    first: MeridianSegment,
    second: MeridianSegment,
) -> bool {
    let mut points = vec![
        point,
        join,
        first.start(),
        first.end(),
        second.start(),
        second.end(),
    ];
    if let MeridianSegment::Arc { center, .. } = first {
        points.push(center);
    }
    if let MeridianSegment::Arc { center, .. } = second {
        points.push(center);
    }
    squared_delta(point, join) <= length_tolerance(&points).powi(2)
}

/// Adjacent features may meet only at their one literal loop join. Any other
/// crossing, tangency, or overlapping interval makes the meridian invalid.
fn adjacent_features_reintersect(
    a: MeridianSegment,
    b: MeridianSegment,
    join: MeridianPoint,
) -> bool {
    match (a, b) {
        (MeridianSegment::Line { .. }, MeridianSegment::Line { .. }) => {
            adjacent_line_line_reintersects(a, b, join)
        }
        (MeridianSegment::Line { start, end }, arc @ MeridianSegment::Arc { .. })
        | (arc @ MeridianSegment::Arc { .. }, MeridianSegment::Line { start, end }) => {
            line_arc_intersection_points(start, end, arc)
                .into_iter()
                .any(|point| !nearby_feature_point(point, join, a, b))
        }
        (a @ MeridianSegment::Arc { .. }, b @ MeridianSegment::Arc { .. }) => {
            adjacent_arc_arc_reintersects(a, b, join)
        }
    }
}

fn adjacent_line_line_reintersects(
    a: MeridianSegment,
    b: MeridianSegment,
    join: MeridianPoint,
) -> bool {
    let MeridianSegment::Line { start: a0, end: a1 } = a else {
        unreachable!()
    };
    let MeridianSegment::Line { start: b0, end: b1 } = b else {
        unreachable!()
    };
    for (point, start, end) in [(a0, b0, b1), (a1, b0, b1), (b0, a0, a1), (b1, a0, a1)] {
        if !nearby_feature_point(point, join, a, b) && on_line(start, end, point) {
            return true;
        }
    }
    let o1 = cross(a0, a1, b0);
    let o2 = cross(a0, a1, b1);
    let o3 = cross(b0, b1, a0);
    let o4 = cross(b0, b1, a1);
    o1 != 0.0
        && o2 != 0.0
        && o3 != 0.0
        && o4 != 0.0
        && o1.signum() != o2.signum()
        && o3.signum() != o4.signum()
}

fn adjacent_arc_arc_reintersects(
    a: MeridianSegment,
    b: MeridianSegment,
    join: MeridianPoint,
) -> bool {
    let (MeridianSegment::Arc { center: ca, .. }, MeridianSegment::Arc { center: cb, .. }) = (a, b)
    else {
        unreachable!()
    };
    let coincident = nearby_feature_point(ca, cb, a, b)
        && (arc_radius(a) - arc_radius(b)).abs()
            <= length_tolerance(&[ca, cb, a.start(), b.start()]);
    if coincident {
        for (point, other) in [(a.start(), b), (a.end(), b), (b.start(), a), (b.end(), a)] {
            if !nearby_feature_point(point, join, a, b) && point_on_arc(other, point) {
                return true;
            }
        }
        return false;
    }
    arc_arc_intersection_points(a, b)
        .into_iter()
        .any(|point| !nearby_feature_point(point, join, a, b))
}

fn point_on_arc(arc: MeridianSegment, point: MeridianPoint) -> bool {
    let MeridianSegment::Arc { center, .. } = arc else {
        unreachable!();
    };
    let radius = arc_radius(arc);
    if ((squared_delta(point, center).sqrt()) - radius).abs() > query_tolerance(point, center) {
        return false;
    }
    arc_contains_angle(
        arc,
        (point.axial - center.axial).atan2(point.radius - center.radius),
    )
}

fn segments_intersect(a: MeridianSegment, b: MeridianSegment) -> bool {
    match (a, b) {
        (
            MeridianSegment::Line { start: a0, end: a1 },
            MeridianSegment::Line { start: b0, end: b1 },
        ) => line_line_intersect(a0, a1, b0, b1),
        (MeridianSegment::Line { start, end }, arc @ MeridianSegment::Arc { .. })
        | (arc @ MeridianSegment::Arc { .. }, MeridianSegment::Line { start, end }) => {
            line_arc_intersect(start, end, arc)
        }
        (a @ MeridianSegment::Arc { .. }, b @ MeridianSegment::Arc { .. }) => {
            arc_arc_intersect(a, b)
        }
    }
}
fn cross(a: MeridianPoint, b: MeridianPoint, c: MeridianPoint) -> f64 {
    (b.radius - a.radius) * (c.axial - a.axial) - (b.axial - a.axial) * (c.radius - a.radius)
}
fn on_line(a: MeridianPoint, b: MeridianPoint, p: MeridianPoint) -> bool {
    cross(a, b, p).abs() <= area_tolerance(&[a, b, p])
        && p.radius >= a.radius.min(b.radius) - query_tolerance(a, b)
        && p.radius <= a.radius.max(b.radius) + query_tolerance(a, b)
        && p.axial >= a.axial.min(b.axial) - query_tolerance(a, b)
        && p.axial <= a.axial.max(b.axial) + query_tolerance(a, b)
}
fn line_line_intersect(
    a0: MeridianPoint,
    a1: MeridianPoint,
    b0: MeridianPoint,
    b1: MeridianPoint,
) -> bool {
    let o1 = cross(a0, a1, b0);
    let o2 = cross(a0, a1, b1);
    let o3 = cross(b0, b1, a0);
    let o4 = cross(b0, b1, a1);
    (o1.signum() != o2.signum() && o3.signum() != o4.signum())
        || on_line(a0, a1, b0)
        || on_line(a0, a1, b1)
        || on_line(b0, b1, a0)
        || on_line(b0, b1, a1)
}
fn line_arc_intersect(start: MeridianPoint, end: MeridianPoint, arc: MeridianSegment) -> bool {
    !line_arc_intersection_points(start, end, arc).is_empty()
}

fn line_arc_intersection_points(
    start: MeridianPoint,
    end: MeridianPoint,
    arc: MeridianSegment,
) -> Vec<MeridianPoint> {
    let MeridianSegment::Arc { center, .. } = arc else {
        unreachable!()
    };
    let dr = end.radius - start.radius;
    let dz = end.axial - start.axial;
    let fr = start.radius - center.radius;
    let fz = start.axial - center.axial;
    let aa = dr.mul_add(dr, dz * dz);
    let bb = 2.0 * (fr * dr + fz * dz);
    let cc = fr.mul_add(fr, fz * fz) - arc_radius(arc).powi(2);
    let disc = bb.mul_add(bb, -4.0 * aa * cc);
    if disc < -discriminant_tolerance(aa, bb, cc) {
        return Vec::new();
    }
    let root = disc.max(0.0).sqrt();
    let mut intersections = Vec::with_capacity(2);
    for u in [(-bb - root) / (2.0 * aa), (-bb + root) / (2.0 * aa)] {
        if (0.0..=1.0).contains(&u) {
            let p = MeridianPoint::new(start.radius + u * dr, start.axial + u * dz);
            let angle = (p.axial - center.axial).atan2(p.radius - center.radius);
            if arc_contains_angle(arc, angle)
                && !intersections
                    .iter()
                    .copied()
                    .any(|prior| nearby_point(prior, p))
            {
                intersections.push(p);
            }
        }
    }
    intersections
}
fn arc_arc_intersect(a: MeridianSegment, b: MeridianSegment) -> bool {
    let (MeridianSegment::Arc { center: ca, .. }, MeridianSegment::Arc { center: cb, .. }) = (a, b)
    else {
        unreachable!()
    };
    let ra = arc_radius(a);
    let rb = arc_radius(b);
    let dx = cb.radius - ca.radius;
    let dy = cb.axial - ca.axial;
    let d = dx.hypot(dy);
    let tolerance = length_tolerance(&[ca, cb, a.start(), b.start()]);
    if d <= tolerance {
        return (ra - rb).abs() <= tolerance;
    }
    !arc_arc_intersection_points(a, b).is_empty()
}

fn arc_arc_intersection_points(a: MeridianSegment, b: MeridianSegment) -> Vec<MeridianPoint> {
    let (MeridianSegment::Arc { center: ca, .. }, MeridianSegment::Arc { center: cb, .. }) = (a, b)
    else {
        unreachable!()
    };
    let ra = arc_radius(a);
    let rb = arc_radius(b);
    let dx = cb.radius - ca.radius;
    let dy = cb.axial - ca.axial;
    let d = dx.hypot(dy);
    let tolerance = length_tolerance(&[ca, cb, a.start(), b.start()]);
    if d <= tolerance || d > ra + rb + tolerance || d < (ra - rb).abs() - tolerance {
        return Vec::new();
    }
    let x = (ra * ra - rb * rb + d * d) / (2.0 * d);
    let h = (ra * ra - x * x).max(0.0).sqrt();
    let ux = dx / d;
    let uy = dy / d;
    let mut intersections = Vec::with_capacity(2);
    for p in [
        MeridianPoint::new(ca.radius + x * ux - h * uy, ca.axial + x * uy + h * ux),
        MeridianPoint::new(ca.radius + x * ux + h * uy, ca.axial + x * uy - h * ux),
    ] {
        let aa = (p.axial - ca.axial).atan2(p.radius - ca.radius);
        let ab = (p.axial - cb.axial).atan2(p.radius - cb.radius);
        if arc_contains_angle(a, aa)
            && arc_contains_angle(b, ab)
            && !intersections
                .iter()
                .copied()
                .any(|prior| nearby_point(prior, p))
        {
            intersections.push(p);
        }
    }
    intersections
}

fn profile_identity(segments: &[MeridianSegment]) -> AxisymmetricIdentity {
    let mut bytes = Vec::with_capacity(16 + segments.len() * 57);
    bytes.extend_from_slice(b"fs-rep-frep-axisymmetric-v1");
    bytes.extend_from_slice(&(segments.len() as u64).to_le_bytes());
    for segment in segments {
        match *segment {
            MeridianSegment::Line { start, end } => {
                bytes.push(0);
                append_point(&mut bytes, start);
                append_point(&mut bytes, end);
            }
            MeridianSegment::Arc {
                start,
                end,
                center,
                clockwise,
            } => {
                bytes.push(1);
                append_point(&mut bytes, start);
                append_point(&mut bytes, end);
                append_point(&mut bytes, center);
                bytes.push(u8::from(clockwise));
            }
        }
    }
    AxisymmetricIdentity(fs_obs::fnv1a64(&bytes))
}
fn append_point(bytes: &mut Vec<u8>, point: MeridianPoint) {
    bytes.extend_from_slice(&point.radius.to_bits().to_le_bytes());
    bytes.extend_from_slice(&point.axial.to_bits().to_le_bytes());
}
