//! Generic certified axisymmetric support, normal, curvature, reach, and gap adapters
//! (bead frankensim-b8bxd.2).
//!
//! This module provides L2 query adapters over validated axisymmetric charts
//! ([`fs_rep_frep::axisymmetric::AxisymmetricChart`]):
//!
//! - [`AxisymmetricSupportMap`]: implements [`ConvexSupportMap`] for convex
//!   axisymmetric charts via global analytic maximization over all meridian
//!   features, with certified support slack and contact inflation;
//! - [`axisymmetric_normal`]: evaluates surface normals with explicit classification
//!   (Smooth, Axis, JoinSetValued);
//! - [`axisymmetric_curvature`]: evaluates meridional and azimuthal principal
//!   curvatures, mean/Gaussian curvatures, and local reach bounds;
//! - [`axisymmetric_reach`]: certified global reach lower bounds;
//! - [`AxisymmetricGapOracle`]: pointwise gap oracle pairing an axisymmetric chart
//!   with another chart under explicit pointwise-only no-claim boundaries.

use crate::{ContactInflation, ConvexOverlapWitness, ConvexSupportMap, GapSample, QueryError};
use fs_exec::Cx;
use fs_geom::{Chart, Point3, Vec3};
use fs_rep_frep::axisymmetric::{
    AxisymmetricChart, AxisymmetricConstructionCertificate, MeridianSegment,
};

/// Binary64 relative slack for support computation.
const AXISYMMETRIC_SUPPORT_SLACK_FACTOR: f64 = 256.0 * f64::EPSILON;

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

        // Validate convexity of the meridian profile
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
                    z_min = z_min.min(start.axial).min(end.axial).min(center.axial - r_arc);
                    z_max = z_max.max(start.axial).max(end.axial).max(center.axial + r_arc);
                    max_r = max_r.max(start.radius).max(end.radius).max(center.radius + r_arc);
                }
            }
        }

        if !z_min.is_finite() || !z_max.is_finite() || !max_r.is_finite() || z_max <= z_min || max_r <= 0.0 {
            return Err(QueryError::ConvexInvalidShape {
                reason: "axisymmetric chart has invalid bounding extents",
            });
        }

        let scale = max_r.max(z_max - z_min);
        let base_slack = (scale * AXISYMMETRIC_SUPPORT_SLACK_FACTOR).next_up();
        let total_slack = base_slack + inflation.radius();

        let interior = Point3::new(0.0, 0.0, (z_min + z_max) * 0.5);

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
        let finite = direction.x.is_finite() && direction.y.is_finite() && direction.z.is_finite();
        if !finite {
            return Point3::new(0.0, 0.0, self.z_max);
        }

        let norm_sq = direction.x * direction.x + direction.y * direction.y + direction.z * direction.z;
        if norm_sq == 0.0 {
            return Point3::new(0.0, 0.0, self.z_max);
        }

        let inv_norm = 1.0 / norm_sq.sqrt();
        let unit_x = direction.x * inv_norm;
        let unit_y = direction.y * inv_norm;
        let unit_z = direction.z * inv_norm;

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
                        let angle_start = (start.axial - center.axial).atan2(start.radius - center.radius);
                        let angle_end = (end.axial - center.axial).atan2(end.radius - center.radius);
                        let angle_cand = (cand_z - center.axial).atan2(cand_r - center.radius);

                        let in_span = if angle_end >= angle_start {
                            angle_cand >= angle_start && angle_cand <= angle_end
                        } else {
                            angle_cand >= angle_start || angle_cand <= angle_end
                        };

                        if in_span {
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

        let (px, py) = if radial_dir > 1e-15 {
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

        // Check if inside axial bounds
        if cz <= self.z_min || cz >= self.z_max || cr >= self.max_radius {
            return None;
        }

        // Compute conservative inradius to boundary features
        let mut min_dist = f64::INFINITY;
        for seg in self.chart.segments() {
            match *seg {
                MeridianSegment::Line { start, end } => {
                    if start.radius < 1e-12 && end.radius < 1e-12 {
                        continue;
                    }
                    let dr = end.radius - start.radius;
                    let dz = end.axial - start.axial;
                    let len_sq = dr * dr + dz * dz;
                    if len_sq > 0.0 {
                        let t = (((cr - start.radius) * dr + (cz - start.axial) * dz) / len_sq).clamp(0.0, 1.0);
                        let proj_r = start.radius + t * dr;
                        let proj_z = start.axial + t * dz;
                        let d = (cr - proj_r).hypot(cz - proj_z);
                        min_dist = min_dist.min(d);
                    }
                }
                MeridianSegment::Arc {
                    start,
                    end: _,
                    center: c,
                    clockwise: _,
                } => {
                    let arc_r = (start.radius - c.radius).hypot(start.axial - c.axial);
                    let d_center = (cr - c.radius).hypot(cz - c.axial);
                    let d = (d_center - arc_r).abs();
                    min_dist = min_dist.min(d);
                }
            }
        }

        let inradius = min_dist - self.inflation.radius();
        if inradius > 0.0 {
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
    /// At a C0 join or corner between features.
    JoinSetValued,
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
/// # Errors
/// [`QueryError::InvalidPointSample`] if the point is non-finite.
pub fn axisymmetric_normal(
    chart: &AxisymmetricChart,
    point: Point3,
    _cx: &Cx<'_>,
) -> Result<AxisymmetricNormal, QueryError> {
    if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
        return Err(QueryError::InvalidPointSample {
            at: [point.x, point.y, point.z],
        });
    }

    let r = point.x.hypot(point.y);
    let z = point.z;

    if r < 1e-12 {
        // Axis point: normal is along z axis
        let normal = if z >= 0.0 {
            Vec3::new(0.0, 0.0, 1.0)
        } else {
            Vec3::new(0.0, 0.0, -1.0)
        };
        return Ok(AxisymmetricNormal {
            normal,
            point,
            classification: NormalClassification::Axis,
            feature_index: 0,
        });
    }

    let cos_theta = point.x / r;
    let sin_theta = point.y / r;

    // Find closest meridian segment
    let mut min_dist_sq = f64::INFINITY;
    let mut best_feature = 0;
    let mut best_normal_rz = (1.0, 0.0);
    let mut is_join = false;

    let segments = chart.segments();
    for (idx, seg) in segments.iter().enumerate() {
        match *seg {
            MeridianSegment::Line { start, end } => {
                if start.radius < 1e-12 && end.radius < 1e-12 {
                    continue;
                }
                let dr = end.radius - start.radius;
                let dz = end.axial - start.axial;
                let len_sq = dr * dr + dz * dz;
                if len_sq > 0.0 {
                    let len = len_sq.sqrt();
                    let t = (((r - start.radius) * dr + (z - start.axial) * dz) / len_sq).clamp(0.0, 1.0);
                    let proj_r = start.radius + t * dr;
                    let proj_z = start.axial + t * dz;
                    let dist_sq = (r - proj_r).powi(2) + (z - proj_z).powi(2);
                    if dist_sq < min_dist_sq {
                        min_dist_sq = dist_sq;
                        best_feature = idx;
                        // Outward normal in (r, z) for CCW segment is (dz / len, -dr / len)
                        best_normal_rz = (dz / len, -dr / len);
                        is_join = t <= 1e-6 || t >= 1.0 - 1e-6;
                    }
                }
            }
            MeridianSegment::Arc {
                start,
                end: _,
                center,
                clockwise: _,
            } => {
                let arc_r = (start.radius - center.radius).hypot(start.axial - center.axial);
                let d_center = (r - center.radius).hypot(z - center.axial);
                let dist_sq = (d_center - arc_r).powi(2);
                if dist_sq < min_dist_sq {
                    min_dist_sq = dist_sq;
                    best_feature = idx;
                    if d_center > 1e-15 {
                        best_normal_rz = (
                            (r - center.radius) / d_center,
                            (z - center.axial) / d_center,
                        );
                    }
                    is_join = false;
                }
            }
        }
    }

    let (nr, nz) = best_normal_rz;
    let normal = Vec3::new(nr * cos_theta, nr * sin_theta, nz);
    let classification = if is_join {
        NormalClassification::JoinSetValued
    } else {
        NormalClassification::Smooth
    };

    Ok(AxisymmetricNormal {
        normal,
        point,
        classification,
        feature_index: best_feature,
    })
}

/// Principal curvature result from [`axisymmetric_curvature`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricCurvature {
    /// Meridional principal curvature (along the meridian profile).
    pub meridional_curvature: f64,
    /// Azimuthal principal curvature (along the circle of revolution).
    pub azimuthal_curvature: f64,
    /// Mean curvature `H = (kappa_m + kappa_theta) / 2`.
    pub mean_curvature: f64,
    /// Gaussian curvature `K = kappa_m * kappa_theta`.
    pub gaussian_curvature: f64,
    /// Certified local reach lower bound.
    pub reach_bound: f64,
    /// Active feature index.
    pub feature_index: usize,
}

/// Evaluate principal curvatures and local reach bound on an axisymmetric chart.
///
/// # Errors
/// [`QueryError::InvalidPointSample`] if the point is non-finite or at the axis.
pub fn axisymmetric_curvature(
    chart: &AxisymmetricChart,
    point: Point3,
    cx: &Cx<'_>,
) -> Result<AxisymmetricCurvature, QueryError> {
    let norm = axisymmetric_normal(chart, point, cx)?;
    let r = point.x.hypot(point.y);

    if r < 1e-12 {
        return Err(QueryError::InvalidPointSample {
            at: [point.x, point.y, point.z],
        });
    }

    let segments = chart.segments();
    let seg = segments.get(norm.feature_index).copied().ok_or_else(|| {
        QueryError::InvalidPointArithmetic {
            reason: "feature index out of bounds",
        }
    })?;

    let (kappa_m, reach_m) = match seg {
        MeridianSegment::Line { .. } => (0.0, f64::INFINITY),
        MeridianSegment::Arc {
            start,
            center,
            clockwise,
            ..
        } => {
            let radius = (start.radius - center.radius).hypot(start.axial - center.axial);
            if radius <= 0.0 {
                (0.0, f64::INFINITY)
            } else {
                let sign = if clockwise { -1.0 } else { 1.0 };
                (sign / radius, radius)
            }
        }
    };

    // Azimuthal curvature kappa_theta = n_r / r (Meusnier theorem)
    let n_r = norm.normal.x.hypot(norm.normal.y);
    let kappa_theta = n_r / r;
    let reach_theta = if kappa_theta.abs() > 1e-12 {
        1.0 / kappa_theta.abs()
    } else {
        f64::INFINITY
    };

    let mean_curvature = (kappa_m + kappa_theta) * 0.5;
    let gaussian_curvature = kappa_m * kappa_theta;
    let reach_bound = reach_m.min(reach_theta).min(r);

    Ok(AxisymmetricCurvature {
        meridional_curvature: kappa_m,
        azimuthal_curvature: kappa_theta,
        mean_curvature,
        gaussian_curvature,
        reach_bound,
        feature_index: norm.feature_index,
    })
}

/// Reach lower bound result from [`axisymmetric_reach`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricReach {
    /// Certified global reach lower bound over the entire revolved shape.
    pub global_reach_lower_bound: f64,
    /// Minimum meridional radius of curvature across filleted features.
    pub min_meridional_radius_of_curvature: f64,
    /// Minimum azimuthal radius of curvature.
    pub min_azimuthal_radius_of_curvature: f64,
}

/// Compute certified reach bounds over an axisymmetric chart.
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
                start,
                end,
                center,
                ..
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
        global_reach_lower_bound: global_reach,
        min_meridional_radius_of_curvature: min_arc_r,
        min_azimuthal_radius_of_curvature: min_r_nonzero,
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
    pub fn new(
        chart_a: &'a AxisymmetricChart,
        chart_b: &'a dyn Chart,
    ) -> Result<Self, QueryError> {
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
