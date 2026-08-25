//! Axisymmetric support, normal, curvature, reach-estimate, and gap adapters
//! (bead frankensim-b8bxd.2).
//!
//! This module provides L2 query adapters over validated axisymmetric charts
//! ([`fs_rep_frep::axisymmetric::AxisymmetricChart`]):
//!
//! - [`AxisymmetricSupportMap`]: implements [`ConvexSupportMap`] for convex
//!   axisymmetric charts via global analytic maximization over all meridian
//!   features, with certified support slack and contact inflation;
//! - [`axisymmetric_normal`]: evaluates unique Smooth/Axis surface normals and
//!   refuses ambiguous feature boundaries;
//! - [`axisymmetric_curvature`]: evaluates Estimate-authority meridional and
//!   azimuthal principal-curvature magnitudes on one smooth feature;
//! - [`axisymmetric_reach`]: exposes a local-feature-scale reach heuristic with
//!   an explicit Estimate-only no-certificate boundary;
//! - [`AxisymmetricGapOracle`]: pointwise gap oracle pairing an axisymmetric chart
//!   with another chart under explicit pointwise-only no-claim boundaries.

use crate::{ContactInflation, ConvexOverlapWitness, ConvexSupportMap, GapSample, QueryError};
use fs_evidence::NumericalKind;
use fs_exec::Cx;
use fs_geom::{Chart, Point3, Vec3};
use fs_rep_frep::axisymmetric::{
    AxisymmetricChart, AxisymmetricConstructionCertificate, AxisymmetricCurvatureError,
    MeridianPoint, MeridianSegment,
};

/// Binary64 relative slack for support computation.
const AXISYMMETRIC_SUPPORT_SLACK_FACTOR: f64 = 256.0 * f64::EPSILON;

/// Normalize without squaring the original components. Scaling first keeps
/// finite nonzero directions usable at both ends of the binary64 range.
fn normalized_direction(direction: Vec3) -> Option<Vec3> {
    if !direction.x.is_finite() || !direction.y.is_finite() || !direction.z.is_finite() {
        return None;
    }
    let scale = direction
        .x
        .abs()
        .max(direction.y.abs())
        .max(direction.z.abs());
    if scale == 0.0 {
        return None;
    }
    let scaled = Vec3::new(
        direction.x / scale,
        direction.y / scale,
        direction.z / scale,
    );
    let norm = scaled.x.hypot(scaled.y).hypot(scaled.z);
    if !norm.is_finite() || norm == 0.0 {
        return None;
    }
    Some(Vec3::new(scaled.x / norm, scaled.y / norm, scaled.z / norm))
}

fn normalized_meridian_tangent(dr: f64, dz: f64) -> Option<(f64, f64)> {
    let scale = dr.abs().max(dz.abs());
    if !scale.is_finite() || scale == 0.0 {
        return None;
    }
    let scaled_dr = dr / scale;
    let scaled_dz = dz / scale;
    let norm = scaled_dr.hypot(scaled_dz);
    if !norm.is_finite() || norm == 0.0 {
        return None;
    }
    Some((scaled_dr / norm, scaled_dz / norm))
}

fn segment_tangent(segment: MeridianSegment, at_end: bool) -> Option<(f64, f64)> {
    match segment {
        MeridianSegment::Line { start, end } => {
            normalized_meridian_tangent(end.radius - start.radius, end.axial - start.axial)
        }
        MeridianSegment::Arc {
            start,
            end,
            center,
            clockwise,
        } => {
            if clockwise {
                return None;
            }
            let point = if at_end { end } else { start };
            normalized_meridian_tangent(-(point.axial - center.axial), point.radius - center.radius)
        }
    }
}

/// A simple CCW profile with only CCW arcs is convex exactly when every join
/// makes a nonnegative turn. The chart constructor already proves closure and
/// rejects self-intersections; this check supplies the missing local-convexity
/// obligation required by `ConvexSupportMap`.
fn profile_is_convex(segments: &[MeridianSegment]) -> bool {
    let turn_tolerance = AXISYMMETRIC_SUPPORT_SLACK_FACTOR;
    segments.iter().enumerate().all(|(index, segment)| {
        let next = segments[(index + 1) % segments.len()];
        let Some((in_r, in_z)) = segment_tangent(*segment, true) else {
            return false;
        };
        let Some((out_r, out_z)) = segment_tangent(next, false) else {
            return false;
        };
        in_r * out_z - in_z * out_r >= -turn_tolerance
    })
}

fn ccw_arc_contains_angle(start: f64, end: f64, candidate: f64) -> bool {
    let sweep = (end - start).rem_euclid(core::f64::consts::TAU);
    let travel = (candidate - start).rem_euclid(core::f64::consts::TAU);
    travel <= sweep + AXISYMMETRIC_SUPPORT_SLACK_FACTOR
}

fn oriented_arc_contains_angle(start: f64, end: f64, candidate: f64, clockwise: bool) -> bool {
    if clockwise {
        ccw_arc_contains_angle(end, start, candidate)
    } else {
        ccw_arc_contains_angle(start, end, candidate)
    }
}

#[derive(Debug, Clone, Copy)]
struct NearestMeridianFeature {
    closest: MeridianPoint,
    distance: f64,
    normal_r: f64,
    normal_z: f64,
    at_boundary: bool,
    feature_scale: f64,
}

fn local_point_tolerance(query: MeridianPoint, closest: MeridianPoint, feature_scale: f64) -> f64 {
    AXISYMMETRIC_SUPPORT_SLACK_FACTOR
        * feature_scale
            .max(query.radius.abs())
            .max(closest.radius.abs())
            .max(f64::MIN_POSITIVE)
}

fn nearest_meridian_feature(
    segment: MeridianSegment,
    query: MeridianPoint,
) -> Option<NearestMeridianFeature> {
    match segment {
        MeridianSegment::Line { start, end } => {
            // The collapsed r=0 profile edge closes the meridian loop but does
            // not generate a surface under revolution.
            if start.radius == 0.0 && end.radius == 0.0 {
                return None;
            }
            let dr = end.radius - start.radius;
            let dz = end.axial - start.axial;
            let length = dr.hypot(dz);
            let (tangent_r, tangent_z) = normalized_meridian_tangent(dr, dz)?;
            if !length.is_finite() || length == 0.0 {
                return None;
            }
            let along =
                (query.radius - start.radius) * tangent_r + (query.axial - start.axial) * tangent_z;
            let fraction = (along / length).clamp(0.0, 1.0);
            let closest =
                MeridianPoint::new(start.radius + fraction * dr, start.axial + fraction * dz);
            let distance = (query.radius - closest.radius).hypot(query.axial - closest.axial);
            let tolerance = local_point_tolerance(query, closest, length);
            let at_boundary = (closest.radius - start.radius).hypot(closest.axial - start.axial)
                <= tolerance
                || (closest.radius - end.radius).hypot(closest.axial - end.axial) <= tolerance;
            Some(NearestMeridianFeature {
                closest,
                distance,
                normal_r: tangent_z,
                normal_z: -tangent_r,
                at_boundary,
                feature_scale: length,
            })
        }
        MeridianSegment::Arc {
            start,
            end,
            center,
            clockwise,
        } => {
            let radius = (start.radius - center.radius).hypot(start.axial - center.axial);
            if !radius.is_finite() || radius == 0.0 {
                return None;
            }
            let vr = query.radius - center.radius;
            let vz = query.axial - center.axial;
            let center_distance = vr.hypot(vz);
            let candidate_angle = vz.atan2(vr);
            let start_angle = (start.axial - center.axial).atan2(start.radius - center.radius);
            let end_angle = (end.axial - center.axial).atan2(end.radius - center.radius);
            let closest = if center_distance > 0.0
                && oriented_arc_contains_angle(start_angle, end_angle, candidate_angle, clockwise)
            {
                MeridianPoint::new(
                    center.radius + radius * candidate_angle.cos(),
                    center.axial + radius * candidate_angle.sin(),
                )
            } else if (query.radius - start.radius).hypot(query.axial - start.axial)
                <= (query.radius - end.radius).hypot(query.axial - end.axial)
            {
                start
            } else {
                end
            };
            let distance = (query.radius - closest.radius).hypot(query.axial - closest.axial);
            let tolerance = local_point_tolerance(query, closest, radius);
            let at_boundary = (closest.radius - start.radius).hypot(closest.axial - start.axial)
                <= tolerance
                || (closest.radius - end.radius).hypot(closest.axial - end.axial) <= tolerance;
            let orientation = if clockwise { -1.0 } else { 1.0 };
            Some(NearestMeridianFeature {
                closest,
                distance,
                normal_r: orientation * (closest.radius - center.radius) / radius,
                normal_z: orientation * (closest.axial - center.axial) / radius,
                at_boundary,
                feature_scale: radius,
            })
        }
    }
}

fn arc_projection_maximum(
    query_r: f64,
    query_z: f64,
    start: MeridianPoint,
    end: MeridianPoint,
    center: MeridianPoint,
) -> Option<f64> {
    let start_angle = (start.axial - center.axial).atan2(start.radius - center.radius);
    let end_angle = (end.axial - center.axial).atan2(end.radius - center.radius);
    let vr = query_r - center.radius;
    let vz = query_z - center.axial;
    if !vr.is_finite() || !vz.is_finite() {
        return None;
    }

    let mut maximum = (vr * start_angle.cos() + vz * start_angle.sin())
        .max(vr * end_angle.cos() + vz * end_angle.sin());
    if vr != 0.0 || vz != 0.0 {
        let candidate = vz.atan2(vr);
        if ccw_arc_contains_angle(start_angle, end_angle, candidate) {
            maximum = maximum.max(vr.hypot(vz));
        }
    }
    maximum.is_finite().then_some(maximum)
}

/// Convex support map adapter over a validated convex axisymmetric chart.
#[derive(Debug, Clone)]
pub struct AxisymmetricSupportMap {
    chart: AxisymmetricChart,
    certificate: AxisymmetricConstructionCertificate,
    inflation: ContactInflation,
    slack: f64,
    interior: Point3,
    z_min: f64,
    z_max: f64,
    max_radius: f64,
}

impl AxisymmetricSupportMap {
    /// Access the construction certificate of the underlying chart.
    #[must_use]
    pub const fn certificate(&self) -> AxisymmetricConstructionCertificate {
        self.certificate
    }
    /// Construct a convex support map from a validated axisymmetric chart.
    ///
    /// # Errors
    /// [`QueryError::ConvexInvalidShape`] if the chart has a central bore/hole
    /// or if its meridian contour is non-convex.
    pub fn try_new(chart: AxisymmetricChart) -> Result<Self, QueryError> {
        Self::try_new_with_inflation(chart, ContactInflation::exact_zero())
    }

    /// Construct a convex support map with certified contact inflation.
    ///
    /// # Errors
    /// [`QueryError::ConvexInvalidShape`] if the chart has a central bore/hole
    /// or if its meridian contour is non-convex.
    pub fn try_new_with_inflation(
        chart: AxisymmetricChart,
        inflation: ContactInflation,
    ) -> Result<Self, QueryError> {
        let cert = chart.construction_certificate();
        if !cert.touches_axis {
            return Err(QueryError::ConvexInvalidShape {
                reason: "axisymmetric chart with a central bore/annulus is non-convex",
            });
        }

        let segments = chart.segments();
        if segments.is_empty() {
            return Err(QueryError::ConvexInvalidShape {
                reason: "axisymmetric chart has no segments",
            });
        }

        if !profile_is_convex(segments) {
            return Err(QueryError::ConvexInvalidShape {
                reason: "axisymmetric meridian profile is not convex",
            });
        }

        let mut z_min = f64::INFINITY;
        let mut z_max = f64::NEG_INFINITY;
        let mut max_r = 0.0f64;

        for seg in segments {
            match *seg {
                MeridianSegment::Line { start, end } => {
                    z_min = z_min.min(start.axial).min(end.axial);
                    z_max = z_max.max(start.axial).max(end.axial);
                    max_r = max_r.max(start.radius).max(end.radius);
                }
                MeridianSegment::Arc {
                    start,
                    end,
                    center,
                    clockwise,
                } => {
                    if clockwise {
                        // Clockwise in CCW meridian loop means reentrant/concave notch!
                        return Err(QueryError::ConvexInvalidShape {
                            reason: "axisymmetric profile has concave (clockwise) arc",
                        });
                    }
                    let r_arc = (start.radius - center.radius).hypot(start.axial - center.axial);
                    z_min = z_min
                        .min(start.axial)
                        .min(end.axial)
                        .min(center.axial - r_arc);
                    z_max = z_max
                        .max(start.axial)
                        .max(end.axial)
                        .max(center.axial + r_arc);
                    max_r = max_r
                        .max(start.radius)
                        .max(end.radius)
                        .max(center.radius + r_arc);
                }
            }
        }

        if !z_min.is_finite()
            || !z_max.is_finite()
            || !max_r.is_finite()
            || z_max <= z_min
            || max_r <= 0.0
        {
            return Err(QueryError::ConvexInvalidShape {
                reason: "axisymmetric chart has invalid bounding extents",
            });
        }

        let axial_span = z_max - z_min;
        let scale = max_r.max(z_min.abs()).max(z_max.abs()).max(axial_span);
        let base_slack = (scale * AXISYMMETRIC_SUPPORT_SLACK_FACTOR).next_up();
        let total_slack = base_slack + inflation.radius();
        let inflated_extent = scale + inflation.radius();
        if !base_slack.is_finite() || !total_slack.is_finite() || !inflated_extent.is_finite() {
            return Err(QueryError::ConvexInvalidShape {
                reason: "axisymmetric chart exceeds finite support arithmetic",
            });
        }

        let interior = Point3::new(0.0, 0.0, f64::midpoint(z_min, z_max));

        Ok(Self {
            chart,
            certificate: cert,
            inflation,
            slack: total_slack,
            interior,
            z_min,
            z_max,
            max_radius: max_r,
        })
    }

    /// Access the underlying axisymmetric chart.
    #[must_use]
    pub const fn chart(&self) -> &AxisymmetricChart {
        &self.chart
    }

    /// Access the contact inflation applied to this support map.
    #[must_use]
    pub const fn inflation(&self) -> ContactInflation {
        self.inflation
    }

    /// Maximum radius of the solid.
    #[must_use]
    pub const fn max_radius(&self) -> f64 {
        self.max_radius
    }

    /// Axial span [z_min, z_max].
    #[must_use]
    pub const fn axial_bounds(&self) -> (f64, f64) {
        (self.z_min, self.z_max)
    }
}

impl ConvexSupportMap for AxisymmetricSupportMap {
    fn support_point(&self, direction: Vec3) -> Point3 {
        let Some(unit) = normalized_direction(direction) else {
            return Point3::new(0.0, 0.0, self.z_max);
        };
        let unit_x = unit.x;
        let unit_y = unit.y;
        let unit_z = unit.z;

        let radial_dir = unit_x.hypot(unit_y);

        // Maximize unit_radial * r + unit_z * z over all meridian segments
        let mut best_val = f64::NEG_INFINITY;
        let mut best_r = 0.0;
        let mut best_z = self.z_max;

        for seg in self.chart.segments() {
            match *seg {
                MeridianSegment::Line { start, end } => {
                    let v1 = radial_dir * start.radius + unit_z * start.axial;
                    let v2 = radial_dir * end.radius + unit_z * end.axial;
                    if v1 > best_val {
                        best_val = v1;
                        best_r = start.radius;
                        best_z = start.axial;
                    }
                    if v2 > best_val {
                        best_val = v2;
                        best_r = end.radius;
                        best_z = end.axial;
                    }
                }
                MeridianSegment::Arc {
                    start,
                    end,
                    center,
                    clockwise: _,
                } => {
                    // Check endpoints
                    let v1 = radial_dir * start.radius + unit_z * start.axial;
                    let v2 = radial_dir * end.radius + unit_z * end.axial;
                    if v1 > best_val {
                        best_val = v1;
                        best_r = start.radius;
                        best_z = start.axial;
                    }
                    if v2 > best_val {
                        best_val = v2;
                        best_r = end.radius;
                        best_z = end.axial;
                    }

                    // Check extremal normal direction on arc
                    let arc_r = (start.radius - center.radius).hypot(start.axial - center.axial);
                    let cand_r = center.radius + arc_r * radial_dir;
                    let cand_z = center.axial + arc_r * unit_z;
                    if cand_r >= 0.0 {
                        // Check if (cand_r, cand_z) lies in the arc's angular interval
                        let angle_start =
                            (start.axial - center.axial).atan2(start.radius - center.radius);
                        let angle_end =
                            (end.axial - center.axial).atan2(end.radius - center.radius);
                        let angle_cand = (cand_z - center.axial).atan2(cand_r - center.radius);

                        if ccw_arc_contains_angle(angle_start, angle_end, angle_cand) {
                            let v_arc = radial_dir * cand_r + unit_z * cand_z;
                            if v_arc > best_val {
                                best_val = v_arc;
                                best_r = cand_r;
                                best_z = cand_z;
                            }
                        }
                    }
                }
            }
        }

        let (px, py) = if radial_dir > 0.0 {
            let cos_theta = unit_x / radial_dir;
            let sin_theta = unit_y / radial_dir;
            (best_r * cos_theta, best_r * sin_theta)
        } else {
            (0.0, 0.0)
        };

        let mut pt = Point3::new(px, py, best_z);
        let infl = self.inflation.radius();
        if infl > 0.0 {
            pt.x += infl * unit_x;
            pt.y += infl * unit_y;
            pt.z += infl * unit_z;
        }
        pt
    }

    fn interior_point(&self) -> Point3 {
        self.interior
    }

    fn support_slack(&self) -> f64 {
        self.slack
    }

    fn contained_ball_radius(&self, center: Point3) -> Option<f64> {
        if !center.x.is_finite() || !center.y.is_finite() || !center.z.is_finite() {
            return None;
        }
        let cr = center.x.hypot(center.y);
        let cz = center.z;

        // A convex profile is the intersection of its oriented supporting
        // half-spaces. The minimum inward margin therefore proves an entire
        // ball, whereas a bounding-cylinder check proves only a coarse box.
        let mut min_margin = f64::INFINITY;
        for seg in self.chart.segments() {
            match *seg {
                MeridianSegment::Line { start, end } => {
                    // The r=0 closure is an artifact of the meridian half-plane,
                    // not a boundary of the revolved three-dimensional solid.
                    if start.radius == 0.0 && end.radius == 0.0 {
                        continue;
                    }
                    let dr = end.radius - start.radius;
                    let dz = end.axial - start.axial;
                    let (tangent_r, tangent_z) = normalized_meridian_tangent(dr, dz)?;
                    let margin = -tangent_z * (cr - start.radius) + tangent_r * (cz - start.axial);
                    if !margin.is_finite() {
                        return None;
                    }
                    min_margin = min_margin.min(margin.next_down());
                }
                MeridianSegment::Arc {
                    start,
                    end,
                    center: c,
                    clockwise,
                } => {
                    if clockwise {
                        return None;
                    }
                    let arc_r = (start.radius - c.radius).hypot(start.axial - c.axial);
                    let projection = arc_projection_maximum(cr, cz, start, end, c)?;
                    let margin = (arc_r - projection).next_down();
                    if !margin.is_finite() {
                        return None;
                    }
                    min_margin = min_margin.min(margin);
                }
            }
        }

        // `slack` includes both the support arithmetic guard and retained
        // contact inflation, so the proof cannot launder either uncertainty.
        let inradius = (min_margin - self.slack).next_down();
        if inradius > 0.0 && inradius.is_finite() {
            Some(inradius)
        } else {
            None
        }
    }

    fn name(&self) -> &'static str {
        "frep/axisymmetric-convex-support"
    }
}

/// Normal classification at a surface point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalClassification {
    /// Smooth C1 interior of a line or arc segment.
    Smooth,
    /// On the axis of revolution (r = 0).
    Axis,
}

/// Normal vector result from [`axisymmetric_normal`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricNormal {
    /// Outward unit normal.
    pub normal: Vec3,
    /// Query surface point.
    pub point: Point3,
    /// Topological/differential classification.
    pub classification: NormalClassification,
    /// Index of the active meridian feature.
    pub feature_index: usize,
}

/// Evaluate surface normal on an axisymmetric chart.
///
/// Feature boundaries with no unique retained normal refuse instead of
/// publishing one arbitrarily selected incident-feature normal.
///
/// # Errors
/// [`QueryError::InvalidPointSample`] if the point is non-finite;
/// [`QueryError::NotOnBoundary`] if it is not on the retained surface;
/// [`QueryError::InvalidPointArithmetic`] at an ambiguous feature boundary;
/// [`QueryError::Cancelled`] when cancellation is observed.
pub fn axisymmetric_normal(
    chart: &AxisymmetricChart,
    point: Point3,
    cx: &Cx<'_>,
) -> Result<AxisymmetricNormal, QueryError> {
    if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
        return Err(QueryError::InvalidPointSample {
            at: [point.x, point.y, point.z],
        });
    }

    if cx.checkpoint().is_err() {
        return Err(QueryError::Cancelled);
    }

    let r = point.x.hypot(point.y);
    let query = MeridianPoint::new(r, point.z);
    let mut best: Option<(usize, NearestMeridianFeature)> = None;
    let mut tied_feature = false;
    let mut ambiguous_tied_normal = false;

    for (index, segment) in chart.segments().iter().copied().enumerate() {
        if index % 16 == 0 && cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let Some(candidate) = nearest_meridian_feature(segment, query) else {
            continue;
        };
        if !candidate.distance.is_finite() {
            return Err(QueryError::InvalidPointArithmetic {
                reason: "axisymmetric surface residual became non-finite",
            });
        }
        match best {
            None => best = Some((index, candidate)),
            Some((_, incumbent)) => {
                let tie_tolerance = AXISYMMETRIC_SUPPORT_SLACK_FACTOR
                    * incumbent
                        .feature_scale
                        .max(candidate.feature_scale)
                        .max(incumbent.distance)
                        .max(candidate.distance)
                        .max(f64::MIN_POSITIVE);
                if candidate.distance + tie_tolerance < incumbent.distance {
                    best = Some((index, candidate));
                    tied_feature = false;
                    ambiguous_tied_normal = false;
                } else if (candidate.distance - incumbent.distance).abs() <= tie_tolerance {
                    tied_feature = true;
                    let normal_delta = (candidate.normal_r - incumbent.normal_r)
                        .hypot(candidate.normal_z - incumbent.normal_z);
                    if normal_delta > AXISYMMETRIC_SUPPORT_SLACK_FACTOR {
                        ambiguous_tied_normal = true;
                    }
                }
            }
        }
    }

    let Some((feature_index, nearest)) = best else {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "axisymmetric chart has no surfaced feature",
        });
    };
    let tolerance = local_point_tolerance(query, nearest.closest, nearest.feature_scale);
    if nearest.distance > tolerance {
        return Err(QueryError::NotOnBoundary {
            sd: nearest.distance,
        });
    }

    if r <= tolerance {
        if ambiguous_tied_normal
            || nearest.normal_r.abs() > AXISYMMETRIC_SUPPORT_SLACK_FACTOR
            || nearest.normal_z.abs() <= AXISYMMETRIC_SUPPORT_SLACK_FACTOR
        {
            return Err(QueryError::InvalidPointArithmetic {
                reason: "axisymmetric axis point has no unique retained normal",
            });
        }
        return Ok(AxisymmetricNormal {
            normal: Vec3::new(0.0, 0.0, nearest.normal_z.signum()),
            point,
            classification: NormalClassification::Axis,
            feature_index,
        });
    }

    if ambiguous_tied_normal || (nearest.at_boundary && !tied_feature) {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "axisymmetric feature boundary has no unique retained normal",
        });
    }

    let cos_theta = point.x / r;
    let sin_theta = point.y / r;
    let normal = normalized_direction(Vec3::new(
        nearest.normal_r * cos_theta,
        nearest.normal_r * sin_theta,
        nearest.normal_z,
    ));
    let Some(normal) = normal else {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "axisymmetric normal became non-finite",
        });
    };

    Ok(AxisymmetricNormal {
        normal,
        point,
        classification: NormalClassification::Smooth,
        feature_index,
    })
}

/// Principal curvature result from [`axisymmetric_curvature`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricCurvature {
    /// Meridional principal-curvature magnitude [m^-1].
    pub meridional_curvature: f64,
    /// Azimuthal principal-curvature magnitude [m^-1].
    pub azimuthal_curvature: f64,
    /// Mean of the two magnitude estimates [m^-1], not oriented mean curvature.
    pub mean_curvature: f64,
    /// Product of the two magnitude estimates [m^-2], not oriented Gaussian curvature.
    pub gaussian_curvature: f64,
    /// Conservative absolute binary64 rounding/conditioning estimate [m^-1].
    pub uncertainty_m_inverse: f64,
    /// Numerical authority, always [`NumericalKind::Estimate`].
    pub authority: NumericalKind,
    /// Active feature index.
    pub feature_index: usize,
}

/// Evaluate principal-curvature magnitudes on one smooth axisymmetric feature.
///
/// This adapter deliberately inherits the underlying chart's Estimate-only
/// authority. It does not turn local binary64 curvature into a reach
/// certificate.
///
/// # Errors
/// The boundary, ambiguity, and cancellation refusals from
/// [`axisymmetric_normal`], plus [`QueryError::InvalidPointArithmetic`] when
/// the selected feature has no stable local curvature estimate.
pub fn axisymmetric_curvature(
    chart: &AxisymmetricChart,
    point: Point3,
    cx: &Cx<'_>,
) -> Result<AxisymmetricCurvature, QueryError> {
    let norm = axisymmetric_normal(chart, point, cx)?;
    if norm.classification != NormalClassification::Smooth {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "axisymmetric curvature is unavailable at the revolution axis",
        });
    }
    let estimate = chart
        .principal_curvatures_at_feature_point(norm.feature_index, point, cx)
        .map_err(|error| match error {
            AxisymmetricCurvatureError::Cancelled => QueryError::Cancelled,
            AxisymmetricCurvatureError::NonFinitePoint { point } => {
                QueryError::InvalidPointSample {
                    at: [point.x, point.y, point.z],
                }
            }
            AxisymmetricCurvatureError::PointNotOnSelectedFeature { residual_m, .. } => {
                QueryError::NotOnBoundary { sd: residual_m }
            }
            _ => QueryError::InvalidPointArithmetic {
                reason: "axisymmetric curvature is unavailable at this feature point",
            },
        })?;
    let kappa_m = estimate.meridional_m_inverse;
    let kappa_theta = estimate.azimuthal_m_inverse;
    let mean_curvature = (kappa_m + kappa_theta) * 0.5;
    let gaussian_curvature = kappa_m * kappa_theta;

    Ok(AxisymmetricCurvature {
        meridional_curvature: kappa_m,
        azimuthal_curvature: kappa_theta,
        mean_curvature,
        gaussian_curvature,
        uncertainty_m_inverse: estimate.uncertainty_m_inverse,
        authority: NumericalKind::Estimate,
        feature_index: norm.feature_index,
    })
}

/// Estimate-only reach-related feature scales from [`axisymmetric_reach`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricReach {
    /// Heuristic global reach scale [m]; not a lower bound or certificate.
    pub global_reach_estimate: f64,
    /// Smallest retained circular-arc radius [m].
    pub min_meridional_radius_estimate: f64,
    /// Smallest positive sampled meridian radius [m].
    pub min_sampled_positive_radius: f64,
    /// Numerical authority, always [`NumericalKind::Estimate`].
    pub authority: NumericalKind,
}

/// Estimate a reach-related local feature scale over an axisymmetric chart.
///
/// This endpoint/arc-radius diagnostic does not analyze the medial axis or
/// nonlocal self-approach, so it must never be consumed as a reach lower bound.
///
/// # Errors
/// [`QueryError::ConvexInvalidShape`] if the chart is empty or degenerate.
pub fn axisymmetric_reach(chart: &AxisymmetricChart) -> Result<AxisymmetricReach, QueryError> {
    let segments = chart.segments();
    if segments.is_empty() {
        return Err(QueryError::ConvexInvalidShape {
            reason: "empty axisymmetric chart",
        });
    }

    let mut min_arc_r = f64::INFINITY;
    let mut min_r_nonzero = f64::INFINITY;

    for seg in segments {
        match *seg {
            MeridianSegment::Line { start, end } => {
                if start.radius > 1e-6 {
                    min_r_nonzero = min_r_nonzero.min(start.radius);
                }
                if end.radius > 1e-6 {
                    min_r_nonzero = min_r_nonzero.min(end.radius);
                }
            }
            MeridianSegment::Arc {
                start, end, center, ..
            } => {
                let arc_r = (start.radius - center.radius).hypot(start.axial - center.axial);
                if arc_r > 0.0 {
                    min_arc_r = min_arc_r.min(arc_r);
                }
                if start.radius > 1e-6 {
                    min_r_nonzero = min_r_nonzero.min(start.radius);
                }
                if end.radius > 1e-6 {
                    min_r_nonzero = min_r_nonzero.min(end.radius);
                }
            }
        }
    }

    let global_reach = min_arc_r.min(min_r_nonzero);

    Ok(AxisymmetricReach {
        global_reach_estimate: global_reach,
        min_meridional_radius_estimate: min_arc_r,
        min_sampled_positive_radius: min_r_nonzero,
        authority: NumericalKind::Estimate,
    })
}

/// Pointwise gap oracle pairing an axisymmetric chart with another chart.
pub struct AxisymmetricGapOracle<'a> {
    chart_a: &'a AxisymmetricChart,
    chart_b: &'a dyn Chart,
    inflation_a: ContactInflation,
    inflation_b: ContactInflation,
    total_inflation: ContactInflation,
}

impl<'a> AxisymmetricGapOracle<'a> {
    /// Construct a gap oracle between an axisymmetric chart and another chart.
    ///
    /// # Errors
    /// [`QueryError::InvalidContactInflation`] if inflation composition fails.
    pub fn new(chart_a: &'a AxisymmetricChart, chart_b: &'a dyn Chart) -> Result<Self, QueryError> {
        Self::new_with_inflation(
            chart_a,
            chart_b,
            ContactInflation::exact_zero(),
            ContactInflation::exact_zero(),
        )
    }

    /// Construct a gap oracle with certified contact inflation.
    ///
    /// # Errors
    /// [`QueryError::InvalidContactInflation`] if inflation composition fails.
    pub fn new_with_inflation(
        chart_a: &'a AxisymmetricChart,
        chart_b: &'a dyn Chart,
        inflation_a: ContactInflation,
        inflation_b: ContactInflation,
    ) -> Result<Self, QueryError> {
        let total_inflation = inflation_a
            .compose(inflation_b)
            .map_err(QueryError::InvalidContactInflation)?;

        Ok(Self {
            chart_a,
            chart_b,
            inflation_a,
            inflation_b,
            total_inflation,
        })
    }

    /// Pointwise gap query at probe `p`.
    ///
    /// Note: this is a strictly POINTWISE query. It does not certify continuous
    /// collision detection, no-tunneling along a trajectory, or moving-pair separation.
    ///
    /// # Errors
    /// [`QueryError::InvalidPointSample`] if the probe is non-finite.
    pub fn gap_at(&self, p: Point3, cx: &Cx<'_>) -> Result<GapSample, QueryError> {
        if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
            return Err(QueryError::InvalidPointSample {
                at: [p.x, p.y, p.z],
            });
        }

        let sa = self.chart_a.eval(p, cx);
        let sb = self.chart_b.eval(p, cx);

        let phi_a = sa.signed_distance;
        let phi_b = sb.signed_distance;

        let infl_a = self.inflation_a.radius();
        let infl_b = self.inflation_b.radius();
        let infl_total = self.total_inflation.radius();

        let sum_lo = (phi_a + phi_b - infl_total).next_down();
        let sum_hi = (phi_a + phi_b + infl_total).next_up();

        let outside_a = phi_a - infl_a > 0.0;
        let outside_b = phi_b - infl_b > 0.0;
        let separation_upper = if outside_a && outside_b {
            Some(sum_hi)
        } else {
            None
        };

        let inside_a_bound = phi_a + infl_a;
        let inside_b_bound = phi_b + infl_b;
        let max_inside = inside_a_bound.max(inside_b_bound);
        let (overlap_inradius, overlap_witness) = if max_inside < 0.0 {
            let r = (-max_inside).next_down();
            (Some(r), ConvexOverlapWitness::from_common_ball(p, r))
        } else {
            (None, None)
        };

        let normal = match (sa.gradient, sb.gradient) {
            (Some(ga), Some(gb)) => {
                let dx = ga.x - gb.x;
                let dy = ga.y - gb.y;
                let dz = ga.z - gb.z;
                let norm = (dx * dx + dy * dy + dz * dz).sqrt();
                if norm.is_finite() && norm > 1e-12 {
                    Some([dx / norm, dy / norm, dz / norm])
                } else {
                    None
                }
            }
            _ => None,
        };

        Ok(GapSample {
            sum_lo,
            sum_hi,
            separation_upper,
            overlap_inradius,
            normal,
            overlap_witness,
        })
    }
}
