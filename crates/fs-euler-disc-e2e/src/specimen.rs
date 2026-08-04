//! Resolved physical specimens for the Euler-disc campaign.
//!
//! A [`DiscProfileSpec`] is deliberately a small *description* of an
//! axisymmetric solid.  Resolving it produces one `AxisymmetricChart` and then
//! obtains mass, centroid, and principal inertia from that exact same line/arc
//! profile.  A dynamics caller must use the resolved chart for support queries
//! as well; hand-entering cylinder or cone inertia alongside a different
//! contact shape is specifically what this boundary prevents.
//!
//! The geometry and mass integrations are analytic over the represented
//! line/arc profile, but their binary64 evaluations carry `Estimate`/roundoff
//! telemetry only.  This module is therefore an input-consistency layer, not a
//! contact-patch, material-calibration, experimental-validation, or
//! configuration-ranking claim.

use core::fmt;

use fs_exec::Cx;
use fs_rep_frep::{
    AxisymmetricChart, AxisymmetricError, AxisymmetricIdentity, AxisymmetricMassError,
    AxisymmetricMassProperties, MeridianPoint, MeridianSegment, SquatDiscEdgeTreatment,
};

/// The bounded, user-facing profile families admitted by the Euler campaign.
///
/// Every variant revolves a single simple meridian about its local `z` axis.
/// The local origin is a construction coordinate, not necessarily the center
/// of mass: [`ResolvedDiscProfile::mass_properties`] is authoritative for the
/// latter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DiscProfileSpec {
    /// A solid squat cylinder with either sharp or true circular-filleted rims.
    SolidCylinder {
        /// Maximum cylindrical radius [m].
        outer_radius_m: f64,
        /// Distance between the two cap planes [m].
        thickness_m: f64,
        /// Exact outer-rim treatment.
        edge_treatment: SquatDiscEdgeTreatment,
    },
    /// A homogeneous annular cylinder.  Its inner bore is physical geometry,
    /// rather than a mass-only adjustment to a solid-cylinder contact shape.
    AnnularCylinder {
        /// Outer cylindrical radius [m].
        outer_radius_m: f64,
        /// Radius of the through bore [m].
        inner_radius_m: f64,
        /// Distance between cap planes [m].
        thickness_m: f64,
    },
    /// A symmetric double-conical or double-frustum profile.
    ///
    /// `face_radius_m == 0` is a true bicone whose two points lie on the
    /// revolution axis.  Positive values create equal planar end faces and
    /// true conical flanks.  This parameterization keeps the material
    /// symmetric about `z = 0` while making contact geometry differ from a
    /// cylinder.
    SymmetricTapered {
        /// Maximum radius at the equatorial plane [m].
        outer_radius_m: f64,
        /// Radius of each planar end face [m].
        face_radius_m: f64,
        /// Axial tip-to-tip/end-face separation [m].
        thickness_m: f64,
    },
    /// A solid cylinder with equal straight conical chamfers at both rims.
    ///
    /// Both chamfer distances must be positive.  Use `SolidCylinder::Sharp`
    /// to request no chamfer; accepting a half-zero chamfer would conceal a
    /// caller-unit or topology mistake.
    ChamferedCylinder {
        /// Maximum radius on the retained cylindrical band [m].
        outer_radius_m: f64,
        /// Distance between cap planes [m].
        thickness_m: f64,
        /// Radial inset from the cylindrical band to either cap edge [m].
        chamfer_radial_m: f64,
        /// Axial run of either conical chamfer [m].
        chamfer_axial_m: f64,
    },
}

/// Geometry extents declared by a resolved profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscProfileDimensions {
    /// Largest radial extent of the profile [m].
    pub outer_radius_m: f64,
    /// Difference between maximum and minimum meridian axial coordinates [m].
    pub thickness_m: f64,
}

/// A profile resolved into one validated geometry and its matching properties.
#[derive(Clone, Debug)]
pub struct ResolvedDiscProfile {
    /// Original bounded parameterization for provenance and JSON reporting.
    pub spec: DiscProfileSpec,
    /// Validated line/arc solid of revolution used for every support query.
    pub chart: AxisymmetricChart,
    /// Homogeneous volumetric density used for the mass integration [kg/m³].
    pub density_kg_per_m3: f64,
    /// Deterministic identity of the exact retained meridian input.
    pub identity: AxisymmetricIdentity,
    /// Stable outer-radius and axial-thickness dimensions [m].
    pub dimensions: DiscProfileDimensions,
    /// Analytic line/arc mass, center of mass, and centroidal inertia.
    pub mass_properties: AxisymmetricMassProperties,
}

/// Refusal from resolving a profile specification.
#[derive(Clone, Debug, PartialEq)]
pub enum DiscProfileError {
    /// A named profile parameter was non-finite or outside its documented domain.
    InvalidParameter { field: &'static str, value: f64 },
    /// A parameter relationship does not describe the named profile family.
    InvalidRelationship { detail: &'static str },
    /// The generic line/arc chart refused the constructed meridian.
    Geometry(AxisymmetricError),
    /// The exact line/arc mass integration refused to publish properties.
    Mass(AxisymmetricMassError),
}

impl fmt::Display for DiscProfileError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParameter { field, value } => {
                write!(
                    formatter,
                    "invalid Euler-disc profile parameter {field}={value}"
                )
            }
            Self::InvalidRelationship { detail } => {
                write!(
                    formatter,
                    "invalid Euler-disc profile relationship: {detail}"
                )
            }
            Self::Geometry(error) => {
                write!(formatter, "Euler-disc profile geometry refused: {error}")
            }
            Self::Mass(error) => write!(formatter, "Euler-disc profile mass refused: {error}"),
        }
    }
}

impl std::error::Error for DiscProfileError {}

impl DiscProfileSpec {
    /// Resolve the specification into a validated profile and matching mass
    /// properties.  The same `Cx` controls chart mass integration and is
    /// retained by callers for later support queries.
    pub fn resolve(
        self,
        density_kg_per_m3: f64,
        cx: &Cx<'_>,
    ) -> Result<ResolvedDiscProfile, DiscProfileError> {
        if !density_kg_per_m3.is_finite() || density_kg_per_m3 <= 0.0 {
            return Err(DiscProfileError::InvalidParameter {
                field: "density_kg_per_m3",
                value: density_kg_per_m3,
            });
        }
        let (chart, dimensions) = self.chart_and_dimensions()?;
        let identity = chart.construction_certificate().identity;
        let mass_properties = chart
            .mass_properties(density_kg_per_m3, cx)
            .map_err(DiscProfileError::Mass)?;
        Ok(ResolvedDiscProfile {
            spec: self,
            chart,
            density_kg_per_m3,
            identity,
            dimensions,
            mass_properties,
        })
    }

    fn chart_and_dimensions(
        self,
    ) -> Result<(AxisymmetricChart, DiscProfileDimensions), DiscProfileError> {
        match self {
            Self::SolidCylinder {
                outer_radius_m,
                thickness_m,
                edge_treatment,
            } => {
                positive("outer_radius_m", outer_radius_m)?;
                positive("thickness_m", thickness_m)?;
                let chart =
                    AxisymmetricChart::squat_disc(outer_radius_m, thickness_m, edge_treatment)
                        .map_err(DiscProfileError::Geometry)?;
                Ok((
                    chart,
                    DiscProfileDimensions {
                        outer_radius_m,
                        thickness_m,
                    },
                ))
            }
            Self::AnnularCylinder {
                outer_radius_m,
                inner_radius_m,
                thickness_m,
            } => {
                positive("outer_radius_m", outer_radius_m)?;
                positive("inner_radius_m", inner_radius_m)?;
                positive("thickness_m", thickness_m)?;
                if inner_radius_m >= outer_radius_m {
                    return Err(DiscProfileError::InvalidRelationship {
                        detail: "annular inner_radius_m must be smaller than outer_radius_m",
                    });
                }
                let half = 0.5 * thickness_m;
                let chart = chart_from_segments(vec![
                    line(inner_radius_m, -half, outer_radius_m, -half),
                    line(outer_radius_m, -half, outer_radius_m, half),
                    line(outer_radius_m, half, inner_radius_m, half),
                    line(inner_radius_m, half, inner_radius_m, -half),
                ])?;
                Ok((
                    chart,
                    DiscProfileDimensions {
                        outer_radius_m,
                        thickness_m,
                    },
                ))
            }
            Self::SymmetricTapered {
                outer_radius_m,
                face_radius_m,
                thickness_m,
            } => {
                positive("outer_radius_m", outer_radius_m)?;
                nonnegative("face_radius_m", face_radius_m)?;
                positive("thickness_m", thickness_m)?;
                if face_radius_m >= outer_radius_m {
                    return Err(DiscProfileError::InvalidRelationship {
                        detail: "symmetric tapered face_radius_m must be smaller than outer_radius_m",
                    });
                }
                let half = 0.5 * thickness_m;
                let segments = if face_radius_m == 0.0 {
                    vec![
                        line(0.0, -half, outer_radius_m, 0.0),
                        line(outer_radius_m, 0.0, 0.0, half),
                        line(0.0, half, 0.0, -half),
                    ]
                } else {
                    vec![
                        line(0.0, -half, face_radius_m, -half),
                        line(face_radius_m, -half, outer_radius_m, 0.0),
                        line(outer_radius_m, 0.0, face_radius_m, half),
                        line(face_radius_m, half, 0.0, half),
                        line(0.0, half, 0.0, -half),
                    ]
                };
                let chart = chart_from_segments(segments)?;
                Ok((
                    chart,
                    DiscProfileDimensions {
                        outer_radius_m,
                        thickness_m,
                    },
                ))
            }
            Self::ChamferedCylinder {
                outer_radius_m,
                thickness_m,
                chamfer_radial_m,
                chamfer_axial_m,
            } => {
                positive("outer_radius_m", outer_radius_m)?;
                positive("thickness_m", thickness_m)?;
                positive("chamfer_radial_m", chamfer_radial_m)?;
                positive("chamfer_axial_m", chamfer_axial_m)?;
                if chamfer_radial_m >= outer_radius_m {
                    return Err(DiscProfileError::InvalidRelationship {
                        detail: "chamfer_radial_m must be smaller than outer_radius_m",
                    });
                }
                let half = 0.5 * thickness_m;
                if chamfer_axial_m > half {
                    return Err(DiscProfileError::InvalidRelationship {
                        detail: "chamfer_axial_m must not exceed thickness_m / 2",
                    });
                }
                let lower_cap_radius = outer_radius_m - chamfer_radial_m;
                let lower_outer_z = -half + chamfer_axial_m;
                let upper_outer_z = half - chamfer_axial_m;
                let mut segments = vec![
                    line(0.0, -half, lower_cap_radius, -half),
                    line(lower_cap_radius, -half, outer_radius_m, lower_outer_z),
                ];
                if chamfer_axial_m < half {
                    segments.push(line(
                        outer_radius_m,
                        lower_outer_z,
                        outer_radius_m,
                        upper_outer_z,
                    ));
                }
                segments.extend([
                    line(outer_radius_m, upper_outer_z, lower_cap_radius, half),
                    line(lower_cap_radius, half, 0.0, half),
                    line(0.0, half, 0.0, -half),
                ]);
                let chart = chart_from_segments(segments)?;
                Ok((
                    chart,
                    DiscProfileDimensions {
                        outer_radius_m,
                        thickness_m,
                    },
                ))
            }
        }
    }
}

fn point(radius: f64, axial: f64) -> MeridianPoint {
    MeridianPoint::new(radius, axial)
}

fn line(start_radius: f64, start_axial: f64, end_radius: f64, end_axial: f64) -> MeridianSegment {
    MeridianSegment::Line {
        start: point(start_radius, start_axial),
        end: point(end_radius, end_axial),
    }
}

fn chart_from_segments(
    segments: Vec<MeridianSegment>,
) -> Result<AxisymmetricChart, DiscProfileError> {
    AxisymmetricChart::try_new(segments).map_err(DiscProfileError::Geometry)
}

fn positive(field: &'static str, value: f64) -> Result<(), DiscProfileError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(DiscProfileError::InvalidParameter { field, value })
    }
}

fn nonnegative(field: &'static str, value: f64) -> Result<(), DiscProfileError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(DiscProfileError::InvalidParameter { field, value })
    }
}
