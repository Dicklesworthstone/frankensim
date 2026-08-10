//! Validated physical and animated cameras for cinematic rendering.
//!
//! The optical model is an ideal thin lens. Focus distance is the positive
//! axial distance from the lens plane, not distance along each raster ray.
//! Lens distortion, chromatic aberration, diffraction, rolling shutter, and
//! autofocus heuristics are deliberately outside this v1 contract.

use core::fmt;

use fs_exec::{Cancelled, Cx};
use fs_geom::{Point3, Vec3};
use fs_math::det;

use crate::charts::Ray;
use crate::motion::ShutterInterval;

const PI: f64 = core::f64::consts::PI;
const LOOK_AT_SINE_MIN: f64 = 1.0e-10;
const UNIT_TOLERANCE: f64 = 1.0e-12;
const MAX_APERTURE_BLADES: u8 = 64;

/// A ranked remediation for a rejected camera input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CameraFix {
    /// Rank one is the most direct suggested correction.
    pub rank: u8,
    /// Stable machine-readable remediation code.
    pub code: &'static str,
    /// Concise actionable guidance.
    pub message: &'static str,
}

/// Fail-closed camera diagnostics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CameraError {
    /// The execution scope requested cancellation.
    Cancelled,
    /// A point, vector, time, or scalar was not finite.
    NonFiniteInput,
    /// Eye and target coincide or a direction has no usable magnitude.
    DegenerateDirection,
    /// The declared up reference is too nearly parallel to the view axis.
    NearlyCollinearUp,
    /// A legacy-compatible basis is not unit and mutually orthogonal.
    InvalidOrientationBasis,
    /// Focal length, sensor height, FOV, or tangent projection is invalid.
    InvalidProjection,
    /// Focus lies on or behind the lens plane.
    InvalidFocusDistance,
    /// Aperture radius, f-number, blade count, or rotation is invalid.
    InvalidAperture,
    /// ISO sensitivity or exposure compensation metadata is invalid.
    InvalidExposureMetadata,
    /// A lens coordinate was outside the half-open unit square.
    InvalidLensSample,
    /// A shot contains no keyframes.
    EmptyShot,
    /// Shot bounds or shot identity are invalid.
    InvalidShot,
    /// Keyframes are not finite, ordered, or contained by their shot.
    InvalidKeyframes,
    /// Optical or focus interpolation policy changes inside a continuous shot.
    IncompatibleKeyframes,
    /// Evaluation would extrapolate outside the admitted camera timeline.
    Extrapolation,
    /// An exposure is not wholly contained by one shot and would cross a cut.
    ShutterCrossesCut,
    /// An exposure token was used outside its admitted time interval.
    InvalidExposure,
    /// Interpolation produced an invalid pose or focus distance.
    InvalidInterpolation,
}

impl CameraError {
    /// Ranked fixes intended for structured admission diagnostics.
    #[must_use]
    #[allow(clippy::too_many_lines)] // one exhaustive error-to-remediation table
    pub fn ranked_fixes(self) -> &'static [CameraFix] {
        const FINITE: &[CameraFix] = &[CameraFix {
            rank: 1,
            code: "replace_non_finite",
            message: "replace NaN or infinity with explicit finite SI values",
        }];
        const DIRECTION: &[CameraFix] = &[
            CameraFix {
                rank: 1,
                code: "separate_eye_target",
                message: "move the target away from the camera eye",
            },
            CameraFix {
                rank: 2,
                code: "supply_view_direction",
                message: "supply a finite nonzero view direction",
            },
        ];
        const UP: &[CameraFix] = &[
            CameraFix {
                rank: 1,
                code: "choose_least_aligned_axis",
                message: "use the world axis least aligned with the view direction as up",
            },
            CameraFix {
                rank: 2,
                code: "declare_roll",
                message: "supply an explicit non-collinear up vector to declare camera roll",
            },
        ];
        const BASIS: &[CameraFix] = &[
            CameraFix {
                rank: 1,
                code: "orthonormalize_basis",
                message: "supply unit forward/up vectors with zero dot product",
            },
            CameraFix {
                rank: 2,
                code: "use_look_at",
                message: "construct the camera from eye, target, and an explicit up reference",
            },
        ];
        const PROJECTION: &[CameraFix] = &[
            CameraFix {
                rank: 1,
                code: "positive_focal_sensor",
                message: "supply positive focal length and sensor height in metres",
            },
            CameraFix {
                rank: 2,
                code: "bounded_vertical_fov",
                message: "supply a vertical FOV strictly between zero and pi radians",
            },
        ];
        const FOCUS: &[CameraFix] = &[
            CameraFix {
                rank: 1,
                code: "focus_in_front",
                message: "place the focus plane at a positive axial distance",
            },
            CameraFix {
                rank: 2,
                code: "move_focus_target",
                message: "move the tracked focus target in front of the lens plane",
            },
        ];
        const APERTURE: &[CameraFix] = &[
            CameraFix {
                rank: 1,
                code: "valid_aperture",
                message: "use a finite nonnegative radius or a finite positive f-number",
            },
            CameraFix {
                rank: 2,
                code: "valid_blade_count",
                message: "use between 3 and 64 aperture blades",
            },
        ];
        const LENS: &[CameraFix] = &[CameraFix {
            rank: 1,
            code: "unit_lens_sample",
            message: "supply both lens coordinates in the half-open interval [0,1)",
        }];
        const EXPOSURE_METADATA: &[CameraFix] = &[
            CameraFix {
                rank: 1,
                code: "positive_iso",
                message: "supply a finite positive ISO sensitivity",
            },
            CameraFix {
                rank: 2,
                code: "finite_exposure_compensation",
                message: "supply finite exposure compensation in EV",
            },
        ];
        const EMPTY_SHOT: &[CameraFix] = &[CameraFix {
            rank: 1,
            code: "add_camera_keyframe",
            message: "add at least one finite camera keyframe inside the shot bounds",
        }];
        const INVALID_SHOT: &[CameraFix] = &[
            CameraFix {
                rank: 1,
                code: "unique_shot_ids",
                message: "assign every shot a unique nonzero identity",
            },
            CameraFix {
                rank: 2,
                code: "partition_shot_timeline",
                message: "order non-overlapping shot bounds with no instant shared by three shots",
            },
        ];
        const SHOT: &[CameraFix] = &[
            CameraFix {
                rank: 1,
                code: "partition_at_cut",
                message: "clip the shutter to one camera shot and render the adjacent shot separately",
            },
            CameraFix {
                rank: 2,
                code: "extend_shot_bounds",
                message: "extend explicit shot bounds without interpolating across a hard cut",
            },
        ];
        const KEYFRAMES: &[CameraFix] = &[
            CameraFix {
                rank: 1,
                code: "strict_keyframe_order",
                message: "sort camera keyframes by unique finite absolute time",
            },
            CameraFix {
                rank: 2,
                code: "split_optical_change",
                message: "split projection or aperture changes into separate shots",
            },
        ];
        const EXPOSURE: &[CameraFix] = &[CameraFix {
            rank: 1,
            code: "readmit_exposure",
            message: "obtain a new exposure token for this camera timeline and shutter",
        }];
        match self {
            Self::Cancelled => &[],
            Self::NonFiniteInput => FINITE,
            Self::DegenerateDirection => DIRECTION,
            Self::NearlyCollinearUp => UP,
            Self::InvalidOrientationBasis => BASIS,
            Self::InvalidProjection => PROJECTION,
            Self::InvalidFocusDistance => FOCUS,
            Self::InvalidAperture => APERTURE,
            Self::InvalidExposureMetadata => EXPOSURE_METADATA,
            Self::InvalidLensSample => LENS,
            Self::EmptyShot => EMPTY_SHOT,
            Self::InvalidShot => INVALID_SHOT,
            Self::InvalidKeyframes | Self::IncompatibleKeyframes | Self::InvalidInterpolation => {
                KEYFRAMES
            }
            Self::Extrapolation | Self::ShutterCrossesCut => SHOT,
            Self::InvalidExposure => EXPOSURE,
        }
    }
}

impl fmt::Display for CameraError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Cancelled => "camera evaluation cancelled",
            Self::NonFiniteInput => "camera input must be finite",
            Self::DegenerateDirection => "camera view direction is degenerate",
            Self::NearlyCollinearUp => "camera up reference is nearly collinear with view",
            Self::InvalidOrientationBasis => "camera basis is not orthonormal",
            Self::InvalidProjection => "camera projection is invalid",
            Self::InvalidFocusDistance => "camera focus must be in front of the lens",
            Self::InvalidAperture => "camera aperture is invalid",
            Self::InvalidExposureMetadata => "camera exposure metadata is invalid",
            Self::InvalidLensSample => "lens sample is outside [0,1)^2",
            Self::EmptyShot => "camera shot has no keyframes",
            Self::InvalidShot => "camera shot bounds or identity are invalid",
            Self::InvalidKeyframes => "camera keyframes are invalid",
            Self::IncompatibleKeyframes => "camera keyframes cannot be interpolated continuously",
            Self::Extrapolation => "camera evaluation would extrapolate",
            Self::ShutterCrossesCut => "camera shutter crosses a hard cut",
            Self::InvalidExposure => "camera exposure token does not cover the requested time",
            Self::InvalidInterpolation => "camera interpolation produced an invalid state",
        })
    }
}

impl core::error::Error for CameraError {}

impl From<Cancelled> for CameraError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

/// Validated vertical projection. The film aspect ratio remains an explicit
/// render setting; no full-frame sensor size is assumed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CameraProjection {
    vertical_half_tan: f64,
    kind: ProjectionKind,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum ProjectionKind {
    FocalSensor {
        focal_length_m: f64,
        sensor_height_m: f64,
    },
    VerticalFov {
        radians: f64,
    },
    HalfTangent,
}

impl CameraProjection {
    /// Admit physical focal length and vertical sensor height in metres.
    pub fn try_focal_sensor(
        focal_length_m: f64,
        sensor_height_m: f64,
    ) -> Result<Self, CameraError> {
        if !focal_length_m.is_finite()
            || !sensor_height_m.is_finite()
            || focal_length_m <= 0.0
            || sensor_height_m <= 0.0
        {
            return Err(CameraError::InvalidProjection);
        }
        let vertical_half_tan = half_ratio(sensor_height_m, focal_length_m);
        if !vertical_half_tan.is_finite() || vertical_half_tan <= 0.0 {
            return Err(CameraError::InvalidProjection);
        }
        Ok(Self {
            vertical_half_tan,
            kind: ProjectionKind::FocalSensor {
                focal_length_m,
                sensor_height_m,
            },
        })
    }

    /// Admit a vertical field of view strictly between zero and pi radians.
    pub fn try_vertical_fov(radians: f64) -> Result<Self, CameraError> {
        if !radians.is_finite() || radians <= 0.0 || radians >= PI {
            return Err(CameraError::InvalidProjection);
        }
        let vertical_half_tan = det::tan(0.5 * radians);
        if !vertical_half_tan.is_finite() || vertical_half_tan <= 0.0 {
            return Err(CameraError::InvalidProjection);
        }
        Ok(Self {
            vertical_half_tan,
            kind: ProjectionKind::VerticalFov { radians },
        })
    }

    /// Admit a direct vertical half tangent. Zero is valid for exact center-ray
    /// fixtures retained by the legacy tracer.
    pub fn try_half_tangent(vertical_half_tan: f64) -> Result<Self, CameraError> {
        if !vertical_half_tan.is_finite() || vertical_half_tan < 0.0 {
            return Err(CameraError::InvalidProjection);
        }
        Ok(Self {
            vertical_half_tan,
            kind: ProjectionKind::HalfTangent,
        })
    }

    /// `tan(vertical_fov/2)` used by raster projection.
    #[must_use]
    pub const fn vertical_half_tan(self) -> f64 {
        self.vertical_half_tan
    }

    /// Physical focal length when the projection declares one.
    #[must_use]
    pub const fn focal_length_m(self) -> Option<f64> {
        match self.kind {
            ProjectionKind::FocalSensor { focal_length_m, .. } => Some(focal_length_m),
            ProjectionKind::VerticalFov { .. } | ProjectionKind::HalfTangent => None,
        }
    }

    /// Active vertical sensor height (m) when physically declared.
    #[must_use]
    pub const fn sensor_height_m(self) -> Option<f64> {
        match self.kind {
            ProjectionKind::FocalSensor {
                sensor_height_m, ..
            } => Some(sensor_height_m),
            ProjectionKind::VerticalFov { .. } | ProjectionKind::HalfTangent => None,
        }
    }

    /// Explicit full vertical FOV (rad) when that parameterization was used.
    #[must_use]
    pub const fn vertical_fov_rad(self) -> Option<f64> {
        match self.kind {
            ProjectionKind::VerticalFov { radians } => Some(radians),
            ProjectionKind::FocalSensor { .. } | ProjectionKind::HalfTangent => None,
        }
    }
}

/// Deterministic lens coordinate in the half-open unit square.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LensSample {
    u: f64,
    v: f64,
}

impl LensSample {
    /// The aperture centre under both supported sampling maps.
    pub const CENTER: Self = Self { u: 0.5, v: 0.5 };

    /// Admit two finite coordinates in `[0,1)`.
    pub fn try_new(u: f64, v: f64) -> Result<Self, CameraError> {
        if !u.is_finite() || !v.is_finite() || !(0.0..1.0).contains(&u) || !(0.0..1.0).contains(&v)
        {
            return Err(CameraError::InvalidLensSample);
        }
        Ok(Self { u, v })
    }

    /// First unit-square coordinate.
    #[must_use]
    pub const fn u(self) -> f64 {
        self.u
    }

    /// Second unit-square coordinate.
    #[must_use]
    pub const fn v(self) -> f64 {
        self.v
    }
}

/// Validated ideal aperture. Polygon radius is the circumradius. Its internal
/// representation is private so callers cannot bypass the constructors and
/// forge invalid blade tables or non-finite radii.
#[derive(Clone, Debug, PartialEq)]
pub struct Aperture {
    kind: ApertureKind,
}

#[derive(Clone, Debug, PartialEq)]
enum ApertureKind {
    Pinhole,
    Circular {
        radius_m: f64,
    },
    RegularPolygon {
        radius_m: f64,
        blades: u8,
        rotation_rad: f64,
        vertices: Vec<[f64; 2]>,
    },
}

impl Aperture {
    /// Admit a circular aperture. A radius of exactly zero is canonical pinhole.
    pub fn try_circular(radius_m: f64) -> Result<Self, CameraError> {
        if !radius_m.is_finite() || radius_m < 0.0 {
            return Err(CameraError::InvalidAperture);
        }
        Ok(if radius_m == 0.0 {
            Self {
                kind: ApertureKind::Pinhole,
            }
        } else {
            Self {
                kind: ApertureKind::Circular { radius_m },
            }
        })
    }

    /// Admit a circular aperture from focal length and f-number.
    pub fn try_from_f_number(focal_length_m: f64, f_number: f64) -> Result<Self, CameraError> {
        if !focal_length_m.is_finite()
            || !f_number.is_finite()
            || focal_length_m <= 0.0
            || f_number <= 0.0
        {
            return Err(CameraError::InvalidAperture);
        }
        let radius_m = half_ratio(focal_length_m, f_number);
        if !radius_m.is_finite() || radius_m == 0.0 {
            return Err(CameraError::InvalidAperture);
        }
        Self::try_circular(radius_m)
    }

    /// Admit a regular bladed aperture. A zero radius is canonical pinhole.
    pub fn try_regular_polygon(
        radius_m: f64,
        blades: u8,
        rotation_rad: f64,
    ) -> Result<Self, CameraError> {
        if !radius_m.is_finite()
            || radius_m < 0.0
            || !rotation_rad.is_finite()
            || !(3..=MAX_APERTURE_BLADES).contains(&blades)
        {
            return Err(CameraError::InvalidAperture);
        }
        if radius_m == 0.0 {
            return Ok(Self {
                kind: ApertureKind::Pinhole,
            });
        }
        let rotation_rad = canonical_zero(rotation_rad.rem_euclid(2.0 * PI));
        let step = 2.0 * PI / f64::from(blades);
        let vertices = (0..blades)
            .map(|index| {
                let theta = rotation_rad + f64::from(index) * step;
                [radius_m * det::cos(theta), radius_m * det::sin(theta)]
            })
            .collect();
        Ok(Self {
            kind: ApertureKind::RegularPolygon {
                radius_m,
                blades,
                rotation_rad,
                vertices,
            },
        })
    }

    /// Whether this is the exact pinhole limit.
    #[must_use]
    pub const fn is_pinhole(&self) -> bool {
        matches!(&self.kind, ApertureKind::Pinhole)
    }

    /// Aperture radius (m).
    #[must_use]
    pub const fn radius_m(&self) -> f64 {
        match &self.kind {
            ApertureKind::Pinhole => 0.0,
            ApertureKind::Circular { radius_m } | ApertureKind::RegularPolygon { radius_m, .. } => {
                *radius_m
            }
        }
    }

    /// Blade count for a regular polygon aperture.
    #[must_use]
    pub const fn blades(&self) -> Option<u8> {
        match &self.kind {
            ApertureKind::RegularPolygon { blades, .. } => Some(*blades),
            ApertureKind::Pinhole | ApertureKind::Circular { .. } => None,
        }
    }

    /// Canonical blade rotation (rad) for a regular polygon aperture.
    #[must_use]
    pub const fn rotation_rad(&self) -> Option<f64> {
        match &self.kind {
            ApertureKind::RegularPolygon { rotation_rad, .. } => Some(*rotation_rad),
            ApertureKind::Pinhole | ApertureKind::Circular { .. } => None,
        }
    }

    fn sample(&self, sample: LensSample) -> [f64; 2] {
        match &self.kind {
            ApertureKind::Pinhole => [0.0, 0.0],
            ApertureKind::Circular { radius_m } => {
                let sx = 2.0 * sample.u - 1.0;
                let sy = 2.0 * sample.v - 1.0;
                if sx == 0.0 && sy == 0.0 {
                    return [0.0, 0.0];
                }
                let (radius, theta) = if sx.abs() > sy.abs() {
                    (sx, (PI / 4.0) * (sy / sx))
                } else {
                    (sy, PI / 2.0 - (PI / 4.0) * (sx / sy))
                };
                [
                    radius_m * radius * det::cos(theta),
                    radius_m * radius * det::sin(theta),
                ]
            }
            ApertureKind::RegularPolygon {
                blades, vertices, ..
            } => {
                let blade_coordinate = f64::from(*blades) * sample.u;
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let index = blade_coordinate as usize;
                let local = blade_coordinate - index as f64;
                let radial = det::sqrt(local);
                let left = vertices[index];
                let right = vertices[(index + 1) % vertices.len()];
                let along = sample.v;
                [
                    radial * ((1.0 - along) * left[0] + along * right[0]),
                    radial * ((1.0 - along) * left[1] + along * right[1]),
                ]
            }
        }
    }
}

/// Metadata retained for later exposure/color pipelines. The current tracer
/// averages radiance samples and does not claim sensor irradiance calibration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureMetadata {
    /// Nominal sensor sensitivity (ISO).
    sensitivity_iso: f64,
    /// Declared artistic exposure compensation (EV).
    compensation_ev: f64,
}

impl ExposureMetadata {
    /// Neutral, explicit metadata used when no photographic exposure transform
    /// is requested.
    pub const NEUTRAL: Self = Self {
        sensitivity_iso: 100.0,
        compensation_ev: 0.0,
    };

    /// Admit finite positive ISO and finite EV compensation.
    pub fn try_new(sensitivity_iso: f64, compensation_ev: f64) -> Result<Self, CameraError> {
        if !sensitivity_iso.is_finite() || sensitivity_iso <= 0.0 || !compensation_ev.is_finite() {
            return Err(CameraError::InvalidExposureMetadata);
        }
        Ok(Self {
            sensitivity_iso,
            compensation_ev: canonical_zero(compensation_ev),
        })
    }

    /// Nominal sensor sensitivity (ISO).
    #[must_use]
    pub const fn sensitivity_iso(self) -> f64 {
        self.sensitivity_iso
    }

    /// Declared artistic exposure compensation (EV).
    #[must_use]
    pub const fn compensation_ev(self) -> f64 {
        self.compensation_ev
    }
}

/// One fully evaluated ideal camera.
#[derive(Clone, Debug, PartialEq)]
pub struct PhysicalCamera {
    eye: Point3,
    forward: Vec3,
    up: Vec3,
    right: Vec3,
    orientation_xyzw: [f64; 4],
    projection: CameraProjection,
    focus_distance_m: f64,
    aperture: Aperture,
    exposure: ExposureMetadata,
}

/// Projection of a world point through the deterministic optical centre.
///
/// This is the geometric reprojection model used by motion vectors. It is
/// intentionally independent of a stochastic thin-lens sample: aperture blur
/// is not object motion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OpticalCenterProjection {
    /// The point is in front of the camera. NDC uses `+x` right and `+y` up;
    /// values outside `[-1, 1]` are valid off-screen projections.
    InFront {
        /// Normalized device coordinates.
        ndc_xy: [f64; 2],
        /// Positive axial camera depth in metres.
        depth_m: f64,
    },
    /// The point lies on or behind the lens plane and has no perspective
    /// raster correspondence.
    BehindCamera {
        /// Nonpositive signed axial camera depth in metres.
        signed_depth_m: f64,
    },
}

/// Pinhole sensor response for a world point that projects inside one raster
/// pixel.
///
/// The density is the exact solid-angle density induced by uniform sampling
/// over that pixel's tangent-plane footprint. Light-subpath estimators use it
/// as the camera-endpoint importance factor so their splats estimate the same
/// per-pixel radiance average as camera-ray sampling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PinholeRasterSample {
    /// Row-major pixel index in `width * height`.
    pub pixel: u32,
    /// Unit direction from the optical centre to the projected world point.
    pub direction_from_camera: Vec3,
    /// Camera-ray probability density with respect to solid angle [sr^-1].
    pub pdf_solid_angle: f64,
    /// Positive axial camera depth [m].
    pub depth_m: f64,
}

#[derive(Clone, Copy)]
struct CameraBasis {
    forward: Vec3,
    up: Vec3,
    right: Vec3,
}

impl PhysicalCamera {
    /// Deterministically choose the world axis least aligned with a finite
    /// nonzero view direction. This is an explicit remediation helper; look-at
    /// construction never invokes it silently because doing so would change
    /// camera roll.
    pub fn suggested_up_reference(view_direction: Vec3) -> Result<Vec3, CameraError> {
        let forward = scaled_unit(view_direction)?;
        let [x, y, z] = [forward.x.abs(), forward.y.abs(), forward.z.abs()];
        Ok(if x <= y && x <= z {
            Vec3::new(1.0, 0.0, 0.0)
        } else if y <= z {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        })
    }

    /// Construct a robust look-at camera. A nearly collinear explicit up vector
    /// refuses rather than silently changing roll.
    pub fn try_look_at(
        eye: Point3,
        target: Point3,
        up_reference: Vec3,
        projection: CameraProjection,
        focus_distance_m: f64,
        aperture: Aperture,
    ) -> Result<Self, CameraError> {
        ensure_point_finite(eye)?;
        ensure_point_finite(target)?;
        let forward = scaled_direction_between(eye, target)?;
        let up_reference = scaled_unit(up_reference)?;
        let right_cross = cross(forward, up_reference);
        if right_cross.x == 0.0 && right_cross.y == 0.0 && right_cross.z == 0.0 {
            return Err(CameraError::NearlyCollinearUp);
        }
        let cross_magnitude = scaled_norm(right_cross)?;
        if cross_magnitude < LOOK_AT_SINE_MIN {
            return Err(CameraError::NearlyCollinearUp);
        }
        let right = scaled_unit(right_cross)?;
        let up = scaled_unit(cross(right, forward))?;
        Self::from_basis(
            eye,
            CameraBasis { forward, up, right },
            projection,
            focus_distance_m,
            aperture,
            ExposureMetadata::NEUTRAL,
        )
    }

    /// Construct an exact legacy-compatible camera from already-orthonormal
    /// forward/up vectors. Their bits are retained for the pinhole branch.
    pub fn try_legacy_compatible(
        eye: Point3,
        forward: Vec3,
        up: Vec3,
        vertical_half_tan: f64,
        focus_distance_m: f64,
        aperture: Aperture,
    ) -> Result<Self, CameraError> {
        ensure_point_finite(eye)?;
        ensure_vec_finite(forward)?;
        ensure_vec_finite(up)?;
        let forward_norm = scaled_norm(forward)?;
        let up_norm = scaled_norm(up)?;
        if (forward_norm - 1.0).abs() > UNIT_TOLERANCE
            || (up_norm - 1.0).abs() > UNIT_TOLERANCE
            || forward.dot(up).abs() > UNIT_TOLERANCE
        {
            return Err(CameraError::InvalidOrientationBasis);
        }
        let right = cross(forward, up);
        if (scaled_norm(right)? - 1.0).abs() > UNIT_TOLERANCE {
            return Err(CameraError::InvalidOrientationBasis);
        }
        Self::from_basis(
            eye,
            CameraBasis { forward, up, right },
            CameraProjection::try_half_tangent(vertical_half_tan)?,
            focus_distance_m,
            aperture,
            ExposureMetadata::NEUTRAL,
        )
    }

    fn from_basis(
        eye: Point3,
        basis: CameraBasis,
        projection: CameraProjection,
        focus_distance_m: f64,
        aperture: Aperture,
        exposure: ExposureMetadata,
    ) -> Result<Self, CameraError> {
        ensure_point_finite(eye)?;
        ensure_vec_finite(basis.forward)?;
        ensure_vec_finite(basis.up)?;
        ensure_vec_finite(basis.right)?;
        if !focus_distance_m.is_finite() || focus_distance_m <= 0.0 {
            return Err(CameraError::InvalidFocusDistance);
        }
        let orientation_xyzw = quaternion_from_basis(basis.right, basis.up, basis.forward)?;
        Ok(Self {
            eye,
            forward: basis.forward,
            up: basis.up,
            right: basis.right,
            orientation_xyzw,
            projection,
            focus_distance_m,
            aperture,
            exposure,
        })
    }

    /// Attach validated exposure metadata without changing ray geometry.
    #[must_use]
    pub fn with_exposure_metadata(mut self, exposure: ExposureMetadata) -> Self {
        self.exposure = exposure;
        self
    }

    /// Camera eye (m).
    #[must_use]
    pub const fn eye(&self) -> Point3 {
        self.eye
    }

    /// Unit optical axis.
    #[must_use]
    pub const fn forward(&self) -> Vec3 {
        self.forward
    }

    /// Unit image-up axis.
    #[must_use]
    pub const fn up(&self) -> Vec3 {
        self.up
    }

    /// Unit image-right axis.
    #[must_use]
    pub const fn right(&self) -> Vec3 {
        self.right
    }

    /// Vertical projection.
    #[must_use]
    pub const fn projection(&self) -> CameraProjection {
        self.projection
    }

    /// Positive axial focus distance (m).
    #[must_use]
    pub const fn focus_distance_m(&self) -> f64 {
        self.focus_distance_m
    }

    /// Ideal aperture geometry.
    #[must_use]
    pub const fn aperture(&self) -> &Aperture {
        &self.aperture
    }

    /// Derived f-number when both a physical focal length and nonzero aperture
    /// radius are declared.
    #[must_use]
    pub fn f_number(&self) -> Option<f64> {
        self.projection
            .focal_length_m()
            .filter(|_| !self.aperture.is_pinhole())
            .and_then(|focal_length_m| {
                let f_number = half_ratio(focal_length_m, self.aperture.radius_m());
                (f_number.is_finite() && f_number > 0.0).then_some(f_number)
            })
    }

    /// Exposure metadata; not applied by the current radiance integrator.
    #[must_use]
    pub const fn exposure_metadata(&self) -> ExposureMetadata {
        self.exposure
    }

    /// Project a finite world point through the optical centre into NDC.
    ///
    /// `aspect_ratio = width / height` is explicit. The legacy zero-half-tan
    /// center-ray fixture is not an invertible image projection and therefore
    /// refuses here. Thin-lens aperture state does not affect this mapping.
    pub fn project_from_optical_center(
        &self,
        world_point: Point3,
        aspect_ratio: f64,
    ) -> Result<OpticalCenterProjection, CameraError> {
        ensure_point_finite(world_point)?;
        let half_tan = self.projection.vertical_half_tan();
        if !aspect_ratio.is_finite()
            || aspect_ratio <= 0.0
            || !half_tan.is_finite()
            || half_tan <= 0.0
        {
            return Err(CameraError::InvalidProjection);
        }
        let from_eye = world_point.delta_from(self.eye);
        ensure_vec_finite(from_eye)?;
        let depth_m = from_eye.dot(self.forward);
        if !depth_m.is_finite() {
            return Err(CameraError::InvalidProjection);
        }
        if depth_m <= 0.0 {
            return Ok(OpticalCenterProjection::BehindCamera {
                signed_depth_m: canonical_zero(depth_m),
            });
        }
        let ndc_xy = [
            (from_eye.dot(self.right) / depth_m) / aspect_ratio / half_tan,
            (from_eye.dot(self.up) / depth_m) / half_tan,
        ];
        if ndc_xy.iter().any(|coordinate| !coordinate.is_finite()) {
            return Err(CameraError::InvalidProjection);
        }
        Ok(OpticalCenterProjection::InFront { ndc_xy, depth_m })
    }

    /// Evaluate the pinhole camera endpoint corresponding to a world point.
    ///
    /// `Ok(None)` means that the point is behind the camera or outside the
    /// half-open raster. A finite aperture refuses because connecting to the
    /// optical centre would not represent that camera's lens-area integral.
    pub fn pinhole_raster_sample(
        &self,
        world_point: Point3,
        width: u32,
        height: u32,
    ) -> Result<Option<PinholeRasterSample>, CameraError> {
        if !self.aperture.is_pinhole() {
            return Err(CameraError::InvalidAperture);
        }
        if width == 0 || height == 0 {
            return Err(CameraError::InvalidProjection);
        }
        let aspect_ratio = f64::from(width) / f64::from(height);
        let OpticalCenterProjection::InFront { ndc_xy, depth_m } =
            self.project_from_optical_center(world_point, aspect_ratio)?
        else {
            return Ok(None);
        };
        // Raster x grows with NDC x, while raster y grows opposite NDC y.
        // Preserve the exact half-open pixel domain in each orientation.
        let x_inside = (-1.0..1.0).contains(&ndc_xy[0]);
        let y_inside = (-1.0..=1.0).contains(&ndc_xy[1]) && ndc_xy[1] != -1.0;
        if !x_inside || !y_inside {
            return Ok(None);
        }
        let raster_x = 0.5 * (ndc_xy[0] + 1.0) * f64::from(width);
        let raster_y = 0.5 * (1.0 - ndc_xy[1]) * f64::from(height);
        if !(raster_x >= 0.0
            && raster_x < f64::from(width)
            && raster_y >= 0.0
            && raster_y < f64::from(height))
        {
            return Ok(None);
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let pixel_x = raster_x as u32;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let pixel_y = raster_y as u32;
        let pixel = pixel_y
            .checked_mul(width)
            .and_then(|row| row.checked_add(pixel_x))
            .ok_or(CameraError::InvalidProjection)?;

        let from_camera = world_point.delta_from(self.eye);
        let distance_m = scaled_norm(from_camera)?;
        let direction_from_camera = from_camera.scale(1.0 / distance_m);
        let optical_cosine = direction_from_camera.dot(self.forward);
        let half_tan = self.projection.vertical_half_tan();
        let tangent_pixel_area =
            4.0 * aspect_ratio * half_tan * half_tan / (f64::from(width) * f64::from(height));
        let optical_cosine_cubed = optical_cosine * optical_cosine * optical_cosine;
        let pdf_solid_angle = 1.0 / (tangent_pixel_area * optical_cosine_cubed);
        if !(optical_cosine > 0.0
            && tangent_pixel_area > 0.0
            && tangent_pixel_area.is_finite()
            && pdf_solid_angle > 0.0
            && pdf_solid_angle.is_finite())
        {
            return Err(CameraError::InvalidProjection);
        }
        Ok(Some(PinholeRasterSample {
            pixel,
            direction_from_camera,
            pdf_solid_angle,
            depth_m,
        }))
    }

    /// Generate a ray from horizontal and vertical tangent offsets. A pinhole
    /// executes the same operation order as the legacy tracer.
    pub fn generate_ray_from_tangent_offsets(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
    ) -> Result<Ray, CameraError> {
        cx.checkpoint()?;
        if !x_tan.is_finite() || !y_tan.is_finite() {
            return Err(CameraError::NonFiniteInput);
        }
        let raster_vector = Vec3::new(
            self.forward.x + x_tan * self.right.x + y_tan * self.up.x,
            self.forward.y + x_tan * self.right.y + y_tan * self.up.y,
            self.forward.z + x_tan * self.right.z + y_tan * self.up.z,
        );
        let ray = if self.aperture.is_pinhole() {
            Ray {
                origin: self.eye,
                dir: direct_unit(raster_vector)?,
            }
        } else {
            // Because dot(forward, raster_vector)=1, multiplying the unnormalised
            // raster vector by axial focus distance lands on the focus plane.
            let focus_point = self.eye.offset(raster_vector.scale(self.focus_distance_m));
            let [lens_x, lens_y] = self.aperture.sample(lens_sample);
            let lens_origin = self
                .eye
                .offset(self.right.scale(lens_x))
                .offset(self.up.scale(lens_y));
            Ray {
                origin: lens_origin,
                dir: scaled_unit(focus_point.delta_from(lens_origin))?,
            }
        };
        cx.checkpoint()?;
        Ok(ray)
    }

    fn with_interpolated_state(
        &self,
        eye: Point3,
        orientation_xyzw: [f64; 4],
        focus_distance_m: f64,
    ) -> Result<Self, CameraError> {
        let (right, up, forward) = basis_from_quaternion(orientation_xyzw)?;
        Self::from_basis(
            eye,
            CameraBasis { forward, up, right },
            self.projection,
            focus_distance_m,
            self.aperture.clone(),
            self.exposure,
        )
    }

    fn optical_state_matches(&self, other: &Self) -> bool {
        self.projection == other.projection
            && self.aperture == other.aperture
            && self.exposure == other.exposure
    }
}

/// Focus interpolation policy declared by each keyframe.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum KeyframeFocus {
    /// Interpolate positive axial focus distances linearly (m).
    AxialDistance,
    /// Interpolate an explicit identity-resolved world point linearly, then
    /// project it onto the evaluated optical axis. This is tracking, not
    /// heuristic autofocus.
    WorldPoint(Point3),
}

/// One exact camera state in a continuous shot.
#[derive(Clone, Debug, PartialEq)]
pub struct CameraKeyframe {
    absolute_time_s: f64,
    camera: PhysicalCamera,
    focus: KeyframeFocus,
}

impl CameraKeyframe {
    /// Admit a keyframe whose scalar focus distance is interpolated directly.
    pub fn try_new(absolute_time_s: f64, camera: PhysicalCamera) -> Result<Self, CameraError> {
        Self::try_with_focus(absolute_time_s, camera, KeyframeFocus::AxialDistance)
    }

    /// Admit a keyframe that tracks an explicit world-space focus point.
    pub fn try_with_world_focus(
        absolute_time_s: f64,
        camera: PhysicalCamera,
        world_point: Point3,
    ) -> Result<Self, CameraError> {
        ensure_point_finite(world_point)?;
        let distance = camera.forward.dot(world_point.delta_from(camera.eye));
        if !distance.is_finite() || distance <= 0.0 {
            return Err(CameraError::InvalidFocusDistance);
        }
        // Only the focus plane changes. Preserve the admitted endpoint pose
        // bits rather than round-tripping the basis through a quaternion.
        let mut camera = camera;
        camera.focus_distance_m = distance;
        Self::try_with_focus(
            absolute_time_s,
            camera,
            KeyframeFocus::WorldPoint(world_point),
        )
    }

    fn try_with_focus(
        absolute_time_s: f64,
        camera: PhysicalCamera,
        focus: KeyframeFocus,
    ) -> Result<Self, CameraError> {
        if !absolute_time_s.is_finite() {
            return Err(CameraError::NonFiniteInput);
        }
        Ok(Self {
            absolute_time_s: canonical_zero(absolute_time_s),
            camera,
            focus,
        })
    }

    /// Absolute keyframe time (s).
    #[must_use]
    pub const fn absolute_time_s(&self) -> f64 {
        self.absolute_time_s
    }

    /// Exact admitted camera state.
    #[must_use]
    pub const fn camera(&self) -> &PhysicalCamera {
        &self.camera
    }

    /// Declared focus interpolation policy.
    #[must_use]
    pub const fn focus(&self) -> KeyframeFocus {
        self.focus
    }
}

/// One continuous camera shot. Pose and focus hold at their first/last
/// keyframe inside explicit shot bounds; there is never extrapolation outside.
#[derive(Clone, Debug, PartialEq)]
pub struct CameraShot {
    shot_id: u64,
    start_s: f64,
    end_s: f64,
    keyframes: Vec<CameraKeyframe>,
}

impl CameraShot {
    /// Admit a continuous shot with strictly ordered keyframes.
    pub fn try_new(
        shot_id: u64,
        start_s: f64,
        end_s: f64,
        keyframes: Vec<CameraKeyframe>,
    ) -> Result<Self, CameraError> {
        let start_s = canonical_zero(start_s);
        let end_s = canonical_zero(end_s);
        if shot_id == 0 || !start_s.is_finite() || !end_s.is_finite() || start_s > end_s {
            return Err(CameraError::InvalidShot);
        }
        if keyframes.is_empty() {
            return Err(CameraError::EmptyShot);
        }
        if keyframes
            .iter()
            .any(|keyframe| keyframe.absolute_time_s < start_s || keyframe.absolute_time_s > end_s)
            || keyframes
                .windows(2)
                .any(|pair| pair[0].absolute_time_s >= pair[1].absolute_time_s)
        {
            return Err(CameraError::InvalidKeyframes);
        }
        for pair in keyframes.windows(2) {
            if !pair[0].camera.optical_state_matches(&pair[1].camera)
                || !focus_modes_match(pair[0].focus, pair[1].focus)
            {
                return Err(CameraError::IncompatibleKeyframes);
            }
        }
        Ok(Self {
            shot_id,
            start_s,
            end_s,
            keyframes,
        })
    }

    /// Stable nonzero shot identity.
    #[must_use]
    pub const fn shot_id(&self) -> u64 {
        self.shot_id
    }

    /// Inclusive shot start (s).
    #[must_use]
    pub const fn start_s(&self) -> f64 {
        self.start_s
    }

    /// Inclusive shot end (s).
    #[must_use]
    pub const fn end_s(&self) -> f64 {
        self.end_s
    }

    /// Ordered exact keyframes.
    #[must_use]
    pub fn keyframes(&self) -> &[CameraKeyframe] {
        &self.keyframes
    }

    fn evaluate(&self, absolute_time_s: f64) -> Result<PhysicalCamera, CameraError> {
        if !absolute_time_s.is_finite() {
            return Err(CameraError::NonFiniteInput);
        }
        let absolute_time_s = canonical_zero(absolute_time_s);
        if absolute_time_s < self.start_s || absolute_time_s > self.end_s {
            return Err(CameraError::Extrapolation);
        }
        if absolute_time_s <= self.keyframes[0].absolute_time_s {
            return Ok(self.keyframes[0].camera.clone());
        }
        let last = self.keyframes.len() - 1;
        if absolute_time_s >= self.keyframes[last].absolute_time_s {
            return Ok(self.keyframes[last].camera.clone());
        }
        match self
            .keyframes
            .binary_search_by(|keyframe| keyframe.absolute_time_s.total_cmp(&absolute_time_s))
        {
            Ok(index) => Ok(self.keyframes[index].camera.clone()),
            Err(right_index) => interpolate_camera_keyframes(
                &self.keyframes[right_index - 1],
                &self.keyframes[right_index],
                absolute_time_s,
            ),
        }
    }
}

/// Which shot owns an exact zero-width exposure at a hard cut.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CutSide {
    /// The outgoing shot immediately before the cut.
    Before,
    /// The entering shot immediately after the cut.
    After,
}

/// Validated exposure-to-shot binding. This prevents a rare exact boundary
/// sample from selecting the opposite side of a cut.
#[derive(Clone, Copy, Debug)]
pub struct CameraExposure<'a> {
    camera: &'a AnimatedCamera,
    shot_index: usize,
    shot_id: u64,
    open_s: f64,
    close_s: f64,
}

/// One evaluated camera plus the continuous-shot identity that owns it.
/// Temporal reprojection compares this identity to refuse vectors across hard
/// cuts rather than interpreting a cut as extreme camera velocity.
#[derive(Clone, Debug, PartialEq)]
pub struct EvaluatedCamera {
    shot_id: u64,
    camera: PhysicalCamera,
}

impl EvaluatedCamera {
    /// Stable nonzero continuous-shot identity.
    #[must_use]
    pub const fn shot_id(&self) -> u64 {
        self.shot_id
    }

    /// Exact camera evaluated at the requested time.
    #[must_use]
    pub const fn camera(&self) -> &PhysicalCamera {
        &self.camera
    }

    /// Consume the tag and return its evaluated camera.
    #[must_use]
    pub fn into_camera(self) -> PhysicalCamera {
        self.camera
    }
}

impl CameraExposure<'_> {
    /// Shot identity owning the complete exposure.
    #[must_use]
    pub const fn shot_id(self) -> u64 {
        self.shot_id
    }
}

/// Ordered camera shots separated by explicit hard cuts or gaps.
#[derive(Clone, Debug, PartialEq)]
pub struct AnimatedCamera {
    shots: Vec<CameraShot>,
}

impl AnimatedCamera {
    /// Admit non-overlapping shots in time order. Equal adjacent endpoints are
    /// an explicit cut, resolved by [`CutSide`].
    pub fn try_new(shots: Vec<CameraShot>) -> Result<Self, CameraError> {
        if shots.is_empty() {
            return Err(CameraError::EmptyShot);
        }
        if shots
            .windows(2)
            .any(|pair| pair[0].end_s > pair[1].start_s || pair[0].shot_id == pair[1].shot_id)
            || shots.windows(3).any(|triple| {
                // Cut ownership is defined at one exact canonical instant,
                // not by an epsilon neighborhood.
                triple[0].end_s.to_bits() == triple[1].start_s.to_bits()
                    && triple[1].start_s.to_bits() == triple[1].end_s.to_bits()
                    && triple[1].end_s.to_bits() == triple[2].start_s.to_bits()
            })
        {
            return Err(CameraError::InvalidShot);
        }
        let mut ids = std::collections::BTreeSet::new();
        if shots.iter().any(|shot| !ids.insert(shot.shot_id)) {
            return Err(CameraError::InvalidShot);
        }
        Ok(Self { shots })
    }

    /// Construct a static camera held across explicit shot bounds.
    pub fn try_static(
        shot_id: u64,
        start_s: f64,
        end_s: f64,
        camera: PhysicalCamera,
    ) -> Result<Self, CameraError> {
        let keyframe_time = canonical_zero(start_s);
        Self::try_new(vec![CameraShot::try_new(
            shot_id,
            start_s,
            end_s,
            vec![CameraKeyframe::try_new(keyframe_time, camera)?],
        )?])
    }

    /// Ordered camera shots.
    #[must_use]
    pub fn shots(&self) -> &[CameraShot] {
        &self.shots
    }

    /// Bind a complete shutter to exactly one shot. A positive exposure may
    /// touch a cut at one endpoint but may not cross it.
    pub fn admit_shutter(
        &self,
        cx: &Cx<'_>,
        shutter: ShutterInterval,
        cut_side: CutSide,
    ) -> Result<CameraExposure<'_>, CameraError> {
        cx.checkpoint()?;
        let mut first = None;
        let mut second = None;
        for (index, shot) in self.shots.iter().enumerate() {
            cx.checkpoint()?;
            if shutter.open_s() >= shot.start_s && shutter.close_s() <= shot.end_s {
                if first.is_none() {
                    first = Some((index, shot));
                } else if second.is_none() {
                    second = Some((index, shot));
                } else {
                    return Err(CameraError::InvalidShot);
                }
            }
        }
        let first = first.ok_or(CameraError::ShutterCrossesCut)?;
        let chosen = match second {
            None => first,
            Some(second) if shutter.duration_s() == 0.0 => match cut_side {
                CutSide::Before => first,
                CutSide::After => second,
            },
            Some(_) => return Err(CameraError::ShutterCrossesCut),
        };
        cx.checkpoint()?;
        Ok(CameraExposure {
            camera: self,
            shot_index: chosen.0,
            shot_id: chosen.1.shot_id,
            open_s: shutter.open_s(),
            close_s: shutter.close_s(),
        })
    }

    /// Evaluate one exact time using a previously admitted exposure binding.
    pub fn evaluate_exposure(
        &self,
        cx: &Cx<'_>,
        exposure: CameraExposure<'_>,
        absolute_time_s: f64,
    ) -> Result<PhysicalCamera, CameraError> {
        cx.checkpoint()?;
        if !core::ptr::eq(self, exposure.camera)
            || !absolute_time_s.is_finite()
            || absolute_time_s < exposure.open_s
            || absolute_time_s > exposure.close_s
        {
            return Err(CameraError::InvalidExposure);
        }
        let shot = self
            .shots
            .get(exposure.shot_index)
            .filter(|shot| shot.shot_id == exposure.shot_id)
            .ok_or(CameraError::InvalidExposure)?;
        let camera = shot.evaluate(absolute_time_s)?;
        cx.checkpoint()?;
        Ok(camera)
    }

    /// Evaluate without an exposure token. At an exact cut, `cut_side` chooses
    /// the outgoing or entering shot; no pose is blended across the cut.
    pub fn evaluate(
        &self,
        cx: &Cx<'_>,
        absolute_time_s: f64,
        cut_side: CutSide,
    ) -> Result<PhysicalCamera, CameraError> {
        self.evaluate_with_shot(cx, absolute_time_s, cut_side)
            .map(EvaluatedCamera::into_camera)
    }

    /// Evaluate one exact time and retain the owning continuous-shot identity.
    /// At an exact cut, `cut_side` selects both the pose and returned identity.
    pub fn evaluate_with_shot(
        &self,
        cx: &Cx<'_>,
        absolute_time_s: f64,
        cut_side: CutSide,
    ) -> Result<EvaluatedCamera, CameraError> {
        cx.checkpoint()?;
        if !absolute_time_s.is_finite() {
            return Err(CameraError::NonFiniteInput);
        }
        let mut matches = self
            .shots
            .iter()
            .filter(|shot| absolute_time_s >= shot.start_s && absolute_time_s <= shot.end_s);
        let first = matches.next().ok_or(CameraError::Extrapolation)?;
        let chosen = match matches.next() {
            None => first,
            Some(second) => match cut_side {
                CutSide::Before => first,
                CutSide::After => second,
            },
        };
        if matches.next().is_some() {
            return Err(CameraError::InvalidShot);
        }
        let camera = chosen.evaluate(absolute_time_s)?;
        cx.checkpoint()?;
        Ok(EvaluatedCamera {
            shot_id: chosen.shot_id,
            camera,
        })
    }
}

fn interpolate_camera_keyframes(
    left: &CameraKeyframe,
    right: &CameraKeyframe,
    absolute_time_s: f64,
) -> Result<PhysicalCamera, CameraError> {
    let alpha =
        unit_interval_fraction(left.absolute_time_s, right.absolute_time_s, absolute_time_s)?;
    if !alpha.is_finite() || !(0.0..=1.0).contains(&alpha) {
        return Err(CameraError::InvalidInterpolation);
    }
    let eye = Point3::new(
        lerp(left.camera.eye.x, right.camera.eye.x, alpha),
        lerp(left.camera.eye.y, right.camera.eye.y, alpha),
        lerp(left.camera.eye.z, right.camera.eye.z, alpha),
    );
    let orientation = slerp_shortest(
        left.camera.orientation_xyzw,
        right.camera.orientation_xyzw,
        alpha,
    )?;
    let (_, _, forward) = basis_from_quaternion(orientation)?;
    let focus_distance_m = match (left.focus, right.focus) {
        (KeyframeFocus::AxialDistance, KeyframeFocus::AxialDistance) => lerp(
            left.camera.focus_distance_m,
            right.camera.focus_distance_m,
            alpha,
        ),
        (KeyframeFocus::WorldPoint(left_target), KeyframeFocus::WorldPoint(right_target)) => {
            let target = Point3::new(
                lerp(left_target.x, right_target.x, alpha),
                lerp(left_target.y, right_target.y, alpha),
                lerp(left_target.z, right_target.z, alpha),
            );
            forward.dot(target.delta_from(eye))
        }
        _ => return Err(CameraError::IncompatibleKeyframes),
    };
    if !focus_distance_m.is_finite() || focus_distance_m <= 0.0 {
        return Err(CameraError::InvalidFocusDistance);
    }
    left.camera
        .with_interpolated_state(eye, orientation, focus_distance_m)
}

fn focus_modes_match(left: KeyframeFocus, right: KeyframeFocus) -> bool {
    matches!(
        (left, right),
        (KeyframeFocus::AxialDistance, KeyframeFocus::AxialDistance)
            | (KeyframeFocus::WorldPoint(_), KeyframeFocus::WorldPoint(_))
    )
}

fn slerp_shortest(
    left: [f64; 4],
    mut right: [f64; 4],
    alpha: f64,
) -> Result<[f64; 4], CameraError> {
    let mut dot = left
        .iter()
        .zip(right)
        .map(|(first, second)| first * second)
        .sum::<f64>();
    if dot < 0.0 {
        for component in &mut right {
            *component = -*component;
        }
        dot = -dot;
    }
    dot = dot.clamp(-1.0, 1.0);
    let mut interpolated = if dot > 1.0 - 1.0e-12 {
        core::array::from_fn(|index| left[index] + alpha * (right[index] - left[index]))
    } else {
        let angle = det::acos(dot);
        let denominator = det::sin(angle);
        if denominator == 0.0 || !denominator.is_finite() {
            return Err(CameraError::InvalidInterpolation);
        }
        let left_weight = det::sin((1.0 - alpha) * angle) / denominator;
        let right_weight = det::sin(alpha * angle) / denominator;
        core::array::from_fn(|index| left_weight.mul_add(left[index], right_weight * right[index]))
    };
    normalize_quaternion(&mut interpolated)?;
    Ok(interpolated)
}

fn quaternion_from_basis(right: Vec3, up: Vec3, forward: Vec3) -> Result<[f64; 4], CameraError> {
    // Camera local axes are +X=right, +Y=up, -Z=forward. Using +Z=forward
    // here would create a reflection rather than a proper rotation.
    let m00 = right.x;
    let m01 = up.x;
    let m02 = -forward.x;
    let m10 = right.y;
    let m11 = up.y;
    let m12 = -forward.y;
    let m20 = right.z;
    let m21 = up.z;
    let m22 = -forward.z;
    let trace = m00 + m11 + m22;
    let mut quaternion = if trace > 0.0 {
        let s = 2.0 * det::sqrt(trace + 1.0);
        [(m21 - m12) / s, (m02 - m20) / s, (m10 - m01) / s, 0.25 * s]
    } else if m00 > m11 && m00 > m22 {
        let s = 2.0 * det::sqrt(1.0 + m00 - m11 - m22);
        [0.25 * s, (m01 + m10) / s, (m02 + m20) / s, (m21 - m12) / s]
    } else if m11 > m22 {
        let s = 2.0 * det::sqrt(1.0 + m11 - m00 - m22);
        [(m01 + m10) / s, 0.25 * s, (m12 + m21) / s, (m02 - m20) / s]
    } else {
        let s = 2.0 * det::sqrt(1.0 + m22 - m00 - m11);
        [(m02 + m20) / s, (m12 + m21) / s, 0.25 * s, (m10 - m01) / s]
    };
    normalize_quaternion(&mut quaternion)?;
    if quaternion[3] < 0.0 {
        for component in &mut quaternion {
            *component = -*component;
        }
    }
    Ok(quaternion)
}

fn basis_from_quaternion(quaternion: [f64; 4]) -> Result<(Vec3, Vec3, Vec3), CameraError> {
    let [x, y, z, w] = quaternion;
    if quaternion.iter().any(|value| !value.is_finite()) {
        return Err(CameraError::InvalidInterpolation);
    }
    let xx = x * x;
    let yy = y * y;
    let zz = z * z;
    let xy = x * y;
    let xz = x * z;
    let yz = y * z;
    let wx = w * x;
    let wy = w * y;
    let wz = w * z;
    let right = Vec3::new(1.0 - 2.0 * (yy + zz), 2.0 * (xy + wz), 2.0 * (xz - wy));
    let up = Vec3::new(2.0 * (xy - wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz + wx));
    let back = Vec3::new(2.0 * (xz + wy), 2.0 * (yz - wx), 1.0 - 2.0 * (xx + yy));
    Ok((
        scaled_unit(right)?,
        scaled_unit(up)?,
        scaled_unit(back.scale(-1.0))?,
    ))
}

fn normalize_quaternion(quaternion: &mut [f64; 4]) -> Result<(), CameraError> {
    let scale = quaternion
        .iter()
        .map(|component| component.abs())
        .fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale == 0.0 {
        return Err(CameraError::InvalidInterpolation);
    }
    let scaled_squared = quaternion
        .iter()
        .map(|component| {
            let scaled = component / scale;
            scaled * scaled
        })
        .sum::<f64>();
    let norm = scale * det::sqrt(scaled_squared);
    if !norm.is_finite() || norm == 0.0 {
        return Err(CameraError::InvalidInterpolation);
    }
    for component in quaternion {
        *component /= norm;
    }
    Ok(())
}

fn direct_unit(vector: Vec3) -> Result<Vec3, CameraError> {
    let norm = vector.norm();
    if !norm.is_finite() || norm == 0.0 {
        return Err(CameraError::DegenerateDirection);
    }
    Ok(vector.scale(1.0 / norm))
}

fn scaled_norm(vector: Vec3) -> Result<f64, CameraError> {
    ensure_vec_finite(vector)?;
    let scale = vector.x.abs().max(vector.y.abs()).max(vector.z.abs());
    if scale == 0.0 {
        return Err(CameraError::DegenerateDirection);
    }
    let scaled = vector.scale(1.0 / scale);
    let norm = scale * det::sqrt(scaled.dot(scaled));
    if !norm.is_finite() || norm == 0.0 {
        return Err(CameraError::DegenerateDirection);
    }
    Ok(norm)
}

fn scaled_unit(vector: Vec3) -> Result<Vec3, CameraError> {
    ensure_vec_finite(vector)?;
    let scale = vector.x.abs().max(vector.y.abs()).max(vector.z.abs());
    if scale == 0.0 {
        return Err(CameraError::DegenerateDirection);
    }
    let scaled = vector.scale(1.0 / scale);
    let norm = det::sqrt(scaled.dot(scaled));
    if !norm.is_finite() || norm == 0.0 {
        return Err(CameraError::DegenerateDirection);
    }
    Ok(scaled.scale(1.0 / norm))
}

fn scaled_direction_between(origin: Point3, target: Point3) -> Result<Vec3, CameraError> {
    ensure_point_finite(origin)?;
    ensure_point_finite(target)?;
    let direct = target.delta_from(origin);
    if direct.x.is_finite() && direct.y.is_finite() && direct.z.is_finite() {
        return scaled_unit(direct);
    }
    let scale = origin
        .x
        .abs()
        .max(origin.y.abs())
        .max(origin.z.abs())
        .max(target.x.abs())
        .max(target.y.abs())
        .max(target.z.abs());
    if scale == 0.0 || !scale.is_finite() {
        return Err(CameraError::DegenerateDirection);
    }
    scaled_unit(Vec3::new(
        target.x / scale - origin.x / scale,
        target.y / scale - origin.y / scale,
        target.z / scale - origin.z / scale,
    ))
}

fn ensure_point_finite(point: Point3) -> Result<(), CameraError> {
    if point.x.is_finite() && point.y.is_finite() && point.z.is_finite() {
        Ok(())
    } else {
        Err(CameraError::NonFiniteInput)
    }
}

fn ensure_vec_finite(vector: Vec3) -> Result<(), CameraError> {
    if vector.x.is_finite() && vector.y.is_finite() && vector.z.is_finite() {
        Ok(())
    } else {
        Err(CameraError::NonFiniteInput)
    }
}

fn cross(left: Vec3, right: Vec3) -> Vec3 {
    Vec3::new(
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x,
    )
}

fn lerp(left: f64, right: f64, alpha: f64) -> f64 {
    (1.0 - alpha).mul_add(left, alpha * right)
}

fn unit_interval_fraction(left: f64, right: f64, value: f64) -> Result<f64, CameraError> {
    let width = right - left;
    let alpha = if width.is_finite() {
        (value - left) / width
    } else {
        // Finite endpoints can still have an infinite direct difference, for
        // example [-f64::MAX, f64::MAX]. Scaling before subtraction preserves
        // the position within that interval without overflow.
        let scale = left.abs().max(right.abs()).max(value.abs());
        let scaled_left = left / scale;
        let scaled_right = right / scale;
        (value / scale - scaled_left) / (scaled_right - scaled_left)
    };
    if !alpha.is_finite() || !(0.0..=1.0).contains(&alpha) {
        return Err(CameraError::InvalidInterpolation);
    }
    Ok(alpha)
}

fn half_ratio(numerator: f64, denominator: f64) -> f64 {
    if numerator >= denominator {
        0.5 / (denominator / numerator)
    } else {
        0.5 * (numerator / denominator)
    }
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}
