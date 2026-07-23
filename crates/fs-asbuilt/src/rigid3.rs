//! 3-D rigid and similarity registration with known correspondences.
//!
//! The rigid fit is the closed-form weighted Kabsch solution: the global
//! minimizer of the declared scalar-weighted least-squares objective over
//! `SO(3) x R^3`, computed through the deterministic one-sided Jacobi SVD in
//! `fs-la`. Right-handed canonicalization of both singular frames (the third
//! column is rebuilt as the cross product of the first two) makes the optimal
//! rotation exactly `V * U^T` in every admitted case, including coplanar
//! rank-2 cross-covariances, without a separate determinant branch.
//!
//! The similarity variant reports its scale estimate WITH a first-order
//! standard error and a units-suspicion diagnostic instead of silently
//! rescaling: a scale far from 1.0 is more often a unit error than a real
//! manufacturing effect, so the diagnostic names the nearest common
//! unit-conversion factor.
//!
//! The calibrated path extends the crate's 2-D metrology machinery to the
//! 6-dof pose. The estimator stays the scalar-weighted closed-form Kabsch fit
//! (base weights `3 / trace(Sigma_i)`, optional deterministic Huber
//! multipliers); its covariance is the correct first-order sandwich for that
//! estimator under the full declared per-fiducial 3x3 covariance model. The
//! estimator is therefore globally optimal in its declared class but is NOT
//! the generalized-least-squares optimum under anisotropic models; that
//! efficiency gap is a documented no-claim, not a hidden assumption.

#![allow(clippy::needless_range_loop)] // Fixed 3x3/3x6/6x6 indices expose the parameter ordering.
#![allow(clippy::float_cmp)] // Exact zeros distinguish structural refusals from IEEE noise.

use crate::uncertainty::HuberPolicy;
use crate::{MAX_AS_BUILT_POINTS, MIN_FIDUCIALS, ScaledSumSquares, canonical_zero};

/// Identity schema for calibrated 3-D registration models.
pub const RIGID3_SCHEMA_VERSION: u32 = 1;

const IDENTITY_DOMAIN: &str = "org.frankensim.fs-asbuilt.rigid3.v1";
const IDENTITY_KIND: &[u8] = b"rigid3-calibrated-registration-v1";
const POLL_STRIDE: usize = 256;
/// Relative spectral-rank gate shared by the design/measured/cross scatter
/// classification, mirroring the 2-D relative rank gate.
const RANK_RELATIVE_TOLERANCE: f64 = 1e-12;

type Vec3 = [f64; 3];
type Mat3 = [[f64; 3]; 3];
type Mat6 = [[f64; 6]; 6];
type Mat36 = [[f64; 6]; 3];

/// Row-major covariance for pose parameters `(tx, ty, tz, rx, ry, rz)`.
///
/// The rotation block is a left rotation-vector perturbation applied about
/// the image of the weighted design centroid, so the translation block is the
/// covariance of the predicted centroid image and the translation/rotation
/// cross blocks vanish in exact arithmetic.
pub type PoseCovariance6 = [[f64; 6]; 6];

/// A structured 3-D registration failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rigid3Error {
    /// Fewer correspondences than the rigid problem needs.
    TooFewFiducials {
        /// Supplied fiducial count.
        have: usize,
        /// Required minimum.
        need: usize,
    },
    /// The bounded point cap was exceeded.
    TooManyPoints {
        /// Supplied count.
        have: usize,
        /// Maximum accepted count.
        max: usize,
    },
    /// Parallel arrays differ in length.
    LengthMismatch {
        /// Array being checked.
        field: &'static str,
        /// Required length.
        expected: usize,
        /// Supplied length.
        found: usize,
    },
    /// A scalar violates a named finite/range invariant.
    InvalidScalar {
        /// Stable field name.
        field: &'static str,
        /// Stable required domain.
        requirement: &'static str,
    },
    /// The design configuration cannot observe a 3-D rotation.
    DegenerateDesign {
        /// Geometric diagnosis of the collapse.
        diagnosis: DegeneracyDiagnosis,
    },
    /// The measured configuration cannot observe a 3-D rotation.
    DegenerateMeasured {
        /// Geometric diagnosis of the collapse.
        diagnosis: DegeneracyDiagnosis,
    },
    /// The centered cross-covariance is rank deficient even though both
    /// configurations have rank two or more; no rotation axis pairing is
    /// observable from these correspondences.
    UnobservableRotation,
    /// The data prefer a reflection and the two smallest cross-covariance
    /// singular values coincide at the stated relative gate, so more than one
    /// proper rotation attains the least-squares optimum. Refused fail-closed,
    /// mirroring the 2-D ambiguous-global-minimum refusal.
    AmbiguousRotation,
    /// A covariance is not strictly positive definite.
    NonPositiveDefiniteCovariance {
        /// Covariance ordinal in its input family.
        index: usize,
    },
    /// The calibration identity is not a valid evidence leaf.
    InvalidCalibrationIdentity {
        /// Stable grammar reason.
        reason: &'static str,
    },
    /// Cross-fiducial dependence was declared but not quantified.
    UnknownDependence,
    /// The pose information matrix is singular or numerically unresolved.
    SingularInformation,
    /// A finite input produced an unrepresentable aggregate.
    ArithmeticOverflow {
        /// Stable aggregate name.
        field: &'static str,
    },
    /// A bounded output allocation could not be reserved.
    AllocationFailed,
    /// Cancellation was observed at a bounded scan boundary.
    Cancelled {
        /// Stable phase name.
        phase: &'static str,
    },
}

impl core::fmt::Display for Rigid3Error {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::TooFewFiducials { have, need } => {
                write!(
                    formatter,
                    "3-D registration needs at least {need} fiducials, got {have}"
                )
            }
            Self::TooManyPoints { have, max } => {
                write!(formatter, "3-D point count {have} exceeds {max}")
            }
            Self::LengthMismatch {
                field,
                expected,
                found,
            } => write!(formatter, "{field} has length {found}; expected {expected}"),
            Self::InvalidScalar { field, requirement } => {
                write!(formatter, "{field} must be {requirement}")
            }
            Self::DegenerateDesign { diagnosis } => {
                write!(formatter, "design configuration is degenerate: {diagnosis}")
            }
            Self::DegenerateMeasured { diagnosis } => {
                write!(
                    formatter,
                    "measured configuration is degenerate: {diagnosis}"
                )
            }
            Self::UnobservableRotation => {
                write!(
                    formatter,
                    "cross-covariance is rank deficient; rotation unobservable"
                )
            }
            Self::AmbiguousRotation => write!(
                formatter,
                "reflection-preferring data with coincident trailing singular values; the optimal proper rotation is not unique"
            ),
            Self::NonPositiveDefiniteCovariance { index } => {
                write!(
                    formatter,
                    "covariance {index} is not strictly positive definite"
                )
            }
            Self::InvalidCalibrationIdentity { reason } => {
                write!(formatter, "calibration identity rejected: {reason}")
            }
            Self::UnknownDependence => {
                write!(
                    formatter,
                    "cross-fiducial dependence declared but not quantified"
                )
            }
            Self::SingularInformation => {
                write!(
                    formatter,
                    "pose information matrix is singular or unresolved"
                )
            }
            Self::ArithmeticOverflow { field } => {
                write!(formatter, "{field} overflowed the finite binary64 range")
            }
            Self::AllocationFailed => write!(formatter, "bounded output allocation failed"),
            Self::Cancelled { phase } => write!(formatter, "cancelled during {phase}"),
        }
    }
}

impl core::error::Error for Rigid3Error {}

/// Geometric diagnosis of a degenerate point configuration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegeneracyDiagnosis {
    /// All points coincide (zero spatial extent).
    Coincident,
    /// The points span only a line; rotation about that line is unobservable.
    Collinear,
}

impl core::fmt::Display for DegeneracyDiagnosis {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Coincident => write!(formatter, "all points coincide"),
            Self::Collinear => write!(formatter, "points are collinear"),
        }
    }
}

/// A finite 3-D point (design or measured coordinate).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Point3 {
    /// Construct a finite point, canonicalizing signed zero.
    ///
    /// # Errors
    /// Refuses NaN or infinite coordinates.
    pub fn new(x: f64, y: f64, z: f64) -> Result<Point3, Rigid3Error> {
        for (field, value) in [("point.x", x), ("point.y", y), ("point.z", z)] {
            if !value.is_finite() {
                return Err(Rigid3Error::InvalidScalar {
                    field,
                    requirement: "finite",
                });
            }
        }
        Ok(Point3 {
            x: canonical_zero(x),
            y: canonical_zero(y),
            z: canonical_zero(z),
        })
    }

    /// x coordinate.
    #[must_use]
    pub const fn x(self) -> f64 {
        self.x
    }

    /// y coordinate.
    #[must_use]
    pub const fn y(self) -> f64 {
        self.y
    }

    /// z coordinate.
    #[must_use]
    pub const fn z(self) -> f64 {
        self.z
    }

    pub(crate) const fn coords(self) -> Vec3 {
        [self.x, self.y, self.z]
    }
}

/// A 3-D fiducial correspondence: a design reference point and where the scan
/// measured it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fiducial3 {
    design: Point3,
    measured: Point3,
}

impl Fiducial3 {
    /// A fiducial correspondence from already-valid typed points.
    #[must_use]
    pub const fn new(design: Point3, measured: Point3) -> Fiducial3 {
        Fiducial3 { design, measured }
    }

    /// Design-time reference location.
    #[must_use]
    pub const fn design(self) -> Point3 {
        self.design
    }

    /// Measured location.
    #[must_use]
    pub const fn measured(self) -> Point3 {
        self.measured
    }
}

// ---------------------------------------------------------------------------
// Small fixed-size linear algebra (module-local; fs-la owns the SVD/eigen).
// ---------------------------------------------------------------------------

pub(crate) fn sub3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(crate) fn add3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub(crate) fn scale3(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

pub(crate) fn dot3(a: Vec3, b: Vec3) -> f64 {
    a[2].mul_add(b[2], a[0].mul_add(b[0], a[1] * b[1]))
}

pub(crate) fn cross3(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1].mul_add(b[2], -(a[2] * b[1])),
        a[2].mul_add(b[0], -(a[0] * b[2])),
        a[0].mul_add(b[1], -(a[1] * b[0])),
    ]
}

pub(crate) fn norm3(a: Vec3) -> f64 {
    dot3(a, a).sqrt()
}

pub(crate) fn mat3_vec(matrix: &Mat3, vector: Vec3) -> Vec3 {
    [
        dot3(matrix[0], vector),
        dot3(matrix[1], vector),
        dot3(matrix[2], vector),
    ]
}

pub(crate) fn mat3_transpose(matrix: &Mat3) -> Mat3 {
    let mut out = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            out[row][column] = matrix[column][row];
        }
    }
    out
}

pub(crate) fn mat3_mul(left: &Mat3, right: &Mat3) -> Mat3 {
    let mut out = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            let mut value = 0.0;
            for inner in 0..3 {
                value = left[row][inner].mul_add(right[inner][column], value);
            }
            out[row][column] = value;
        }
    }
    out
}

pub(crate) fn det3(matrix: &Mat3) -> f64 {
    dot3(matrix[0], cross3(matrix[1], matrix[2]))
}

/// Rotation angle in radians recovered from a rotation matrix trace, clamped
/// into the representable arccosine domain.
pub(crate) fn rotation_angle_from_matrix(rotation: &Mat3) -> f64 {
    let trace = rotation[0][0] + rotation[1][1] + rotation[2][2];
    let cosine = ((trace - 1.0) / 2.0).clamp(-1.0, 1.0);
    cosine.acos()
}

/// Invert a symmetric positive-definite matrix through a diagonally
/// equilibrated Cholesky factorization, generalizing the crate's 2-D
/// `invert_spd3` to any fixed dimension.
fn invert_spd<const N: usize>(matrix: &[[f64; N]; N]) -> Result<[[f64; N]; N], Rigid3Error> {
    if !matrix.iter().flatten().all(|value| value.is_finite()) {
        return Err(Rigid3Error::ArithmeticOverflow {
            field: "information matrix",
        });
    }
    let mut scales = [0.0; N];
    for index in 0..N {
        let diagonal = matrix[index][index];
        if !diagonal.is_finite() || diagonal <= 0.0 {
            return Err(Rigid3Error::SingularInformation);
        }
        scales[index] = diagonal.sqrt();
    }
    let mut normalized = [[0.0; N]; N];
    for row in 0..N {
        for column in 0..N {
            // Successive division avoids forming a possibly overflowing
            // product of two dimensional scales.
            normalized[row][column] = matrix[row][column] / scales[row] / scales[column];
            if !normalized[row][column].is_finite() {
                return Err(Rigid3Error::SingularInformation);
            }
        }
    }
    let mut lower = [[0.0; N]; N];
    let pivot_floor = 256.0 * f64::EPSILON;
    for row in 0..N {
        for column in 0..=row {
            let mut value = normalized[row][column];
            for prior in 0..column {
                value -= lower[row][prior] * lower[column][prior];
            }
            if row == column {
                if !value.is_finite() || value <= pivot_floor {
                    return Err(Rigid3Error::SingularInformation);
                }
                lower[row][column] = value.sqrt();
            } else {
                lower[row][column] = value / lower[column][column];
            }
        }
    }
    let mut inverse_normalized = [[0.0; N]; N];
    for basis in 0..N {
        let mut forward = [0.0; N];
        for row in 0..N {
            let mut value = if row == basis { 1.0 } else { 0.0 };
            for prior in 0..row {
                value -= lower[row][prior] * forward[prior];
            }
            forward[row] = value / lower[row][row];
        }
        let mut backward = [0.0; N];
        for row in (0..N).rev() {
            let mut value = forward[row];
            for later in row + 1..N {
                value -= lower[later][row] * backward[later];
            }
            backward[row] = value / lower[row][row];
        }
        for row in 0..N {
            inverse_normalized[row][basis] = backward[row];
        }
    }
    let mut inverse = [[0.0; N]; N];
    for row in 0..N {
        for column in 0..N {
            inverse[row][column] = inverse_normalized[row][column] / scales[row] / scales[column];
        }
    }
    if inverse.iter().flatten().all(|value| value.is_finite()) {
        Ok(inverse)
    } else {
        Err(Rigid3Error::ArithmeticOverflow {
            field: "inverse information matrix",
        })
    }
}

/// Midpoint-symmetrize and revalidate a covariance as strictly positive
/// definite, generalizing the crate's 2-D publication gate.
fn symmetrize_and_validate_spd<const N: usize>(
    mut covariance: [[f64; N]; N],
) -> Result<[[f64; N]; N], Rigid3Error> {
    if !covariance.iter().flatten().all(|value| value.is_finite()) {
        return Err(Rigid3Error::ArithmeticOverflow {
            field: "pose covariance",
        });
    }
    for row in 0..N {
        if covariance[row][row] <= 0.0 {
            return Err(Rigid3Error::SingularInformation);
        }
        for column in row + 1..N {
            let symmetric = f64::midpoint(covariance[row][column], covariance[column][row]);
            covariance[row][column] = symmetric;
            covariance[column][row] = symmetric;
        }
    }
    let mut scales = [0.0; N];
    for index in 0..N {
        scales[index] = covariance[index][index].sqrt();
    }
    let mut lower = [[0.0; N]; N];
    let pivot_floor = 256.0 * f64::EPSILON;
    for row in 0..N {
        for column in 0..=row {
            let mut value = covariance[row][column] / scales[row] / scales[column];
            for prior in 0..column {
                value -= lower[row][prior] * lower[column][prior];
            }
            if row == column {
                if !value.is_finite() || value <= pivot_floor {
                    return Err(Rigid3Error::SingularInformation);
                }
                lower[row][column] = value.sqrt();
            } else {
                lower[row][column] = value / lower[column][column];
            }
        }
    }
    Ok(covariance)
}

fn checkpoint(
    cx: &fs_exec::Cx<'_>,
    ordinal: usize,
    phase: &'static str,
) -> Result<(), Rigid3Error> {
    if ordinal.is_multiple_of(POLL_STRIDE) {
        cx.checkpoint()
            .map_err(|_| Rigid3Error::Cancelled { phase })?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Weighted Kabsch core.
// ---------------------------------------------------------------------------

/// Spectral conditioning of the admitted registration, published so callers
/// can see how close the configuration sat to each refusal gate instead of
/// trusting a boolean.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RegistrationCondition {
    design_spectrum: Vec3,
    measured_spectrum: Vec3,
    cross_singular_values: Vec3,
    coplanar_design: bool,
    coplanar_cross: bool,
    reflection_preferred: bool,
}

impl RegistrationCondition {
    /// Descending eigenvalues of the extent-normalized, weighted design
    /// scatter matrix.
    #[must_use]
    pub const fn design_spectrum(&self) -> Vec3 {
        self.design_spectrum
    }

    /// Descending eigenvalues of the extent-normalized, weighted measured
    /// scatter matrix.
    #[must_use]
    pub const fn measured_spectrum(&self) -> Vec3 {
        self.measured_spectrum
    }

    /// Descending singular values of the extent-normalized, weighted
    /// design/measured cross-covariance.
    #[must_use]
    pub const fn cross_singular_values(&self) -> Vec3 {
        self.cross_singular_values
    }

    /// Whether the design configuration is coplanar at the relative rank
    /// gate. Coplanar configurations are accepted: three non-collinear points
    /// determine the rigid pose.
    #[must_use]
    pub const fn coplanar_design(&self) -> bool {
        self.coplanar_design
    }

    /// Whether the cross-covariance is rank two at the relative rank gate.
    #[must_use]
    pub const fn coplanar_cross(&self) -> bool {
        self.coplanar_cross
    }

    /// True when the unconstrained orthogonal optimum is a reflection with a
    /// clearly nonzero trailing singular value: the proper rotation returned
    /// is still the constrained optimum, but mirrored data (a left-handed
    /// scan frame or a mirrored part) is the usual physical cause and is
    /// worth investigating before trusting residuals.
    #[must_use]
    pub const fn reflection_preferred(&self) -> bool {
        self.reflection_preferred
    }
}

struct KabschCore {
    rotation: Mat3,
    centroid_design: Vec3,
    centroid_measured: Vec3,
    /// `sigma1 + sigma2 + d * sigma3` in normalized units (Umeyama numerator).
    signed_sigma_sum: f64,
    /// `sum_i w_i ||dp_i||^2` in normalized units (Umeyama denominator).
    design_scatter: f64,
    scale_design: f64,
    scale_measured: f64,
    condition: RegistrationCondition,
}

struct CloudFrame {
    anchor: Vec3,
    scale: f64,
}

fn cloud_frame(
    fiducials: &[Fiducial3],
    select_measured: bool,
    cx: &fs_exec::Cx<'_>,
    phase: &'static str,
) -> Result<CloudFrame, Rigid3Error> {
    let mut minimum = [f64::INFINITY; 3];
    let mut maximum = [f64::NEG_INFINITY; 3];
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, phase)?;
        let point = if select_measured {
            fiducial.measured.coords()
        } else {
            fiducial.design.coords()
        };
        for axis in 0..3 {
            minimum[axis] = minimum[axis].min(point[axis]);
            maximum[axis] = maximum[axis].max(point[axis]);
        }
    }
    let mut anchor = [0.0; 3];
    let mut scale = 0.0f64;
    for axis in 0..3 {
        anchor[axis] = f64::midpoint(minimum[axis], maximum[axis]);
        scale = scale.max(maximum[axis] - minimum[axis]);
    }
    if !scale.is_finite() {
        return Err(Rigid3Error::ArithmeticOverflow {
            field: "point-set extent",
        });
    }
    Ok(CloudFrame { anchor, scale })
}

fn classify_spectrum(spectrum: Vec3) -> Result<bool, DegeneracyDiagnosis> {
    if spectrum[0] <= 0.0 {
        return Err(DegeneracyDiagnosis::Coincident);
    }
    if spectrum[1] <= RANK_RELATIVE_TOLERANCE * spectrum[0] {
        return Err(DegeneracyDiagnosis::Collinear);
    }
    Ok(spectrum[2] <= RANK_RELATIVE_TOLERANCE * spectrum[0])
}

fn spectrum_of_scatter(scatter: &Mat3) -> Vec3 {
    let flat = [
        scatter[0][0],
        scatter[0][1],
        scatter[0][2],
        scatter[1][0],
        scatter[1][1],
        scatter[1][2],
        scatter[2][0],
        scatter[2][1],
        scatter[2][2],
    ];
    let svd = fs_la::factor::svd_jacobi(&flat, 3, 3);
    [svd.sigma[0], svd.sigma[1], svd.sigma[2]]
}

fn column(flat: &[f64], index: usize) -> Vec3 {
    [flat[index], flat[3 + index], flat[6 + index]]
}

#[allow(clippy::too_many_lines)]
fn weighted_kabsch(
    fiducials: &[Fiducial3],
    weights: &[f64],
    cx: &fs_exec::Cx<'_>,
) -> Result<KabschCore, Rigid3Error> {
    if fiducials.len() < MIN_FIDUCIALS {
        return Err(Rigid3Error::TooFewFiducials {
            have: fiducials.len(),
            need: MIN_FIDUCIALS,
        });
    }
    if fiducials.len() > MAX_AS_BUILT_POINTS {
        return Err(Rigid3Error::TooManyPoints {
            have: fiducials.len(),
            max: MAX_AS_BUILT_POINTS,
        });
    }
    if weights.len() != fiducials.len() {
        return Err(Rigid3Error::LengthMismatch {
            field: "weights",
            expected: fiducials.len(),
            found: weights.len(),
        });
    }
    for weight in weights {
        if !weight.is_finite() || *weight <= 0.0 {
            return Err(Rigid3Error::InvalidScalar {
                field: "weight",
                requirement: "finite and positive",
            });
        }
    }

    let design_frame = cloud_frame(fiducials, false, cx, "rigid3.design-extent")?;
    let measured_frame = cloud_frame(fiducials, true, cx, "rigid3.measured-extent")?;
    if design_frame.scale == 0.0 {
        return Err(Rigid3Error::DegenerateDesign {
            diagnosis: DegeneracyDiagnosis::Coincident,
        });
    }
    if measured_frame.scale == 0.0 {
        return Err(Rigid3Error::DegenerateMeasured {
            diagnosis: DegeneracyDiagnosis::Coincident,
        });
    }

    // Weighted running means of the extent-normalized anchored deviations.
    let mut weight_total = 0.0f64;
    let mut mean_design = [0.0f64; 3];
    let mut mean_measured = [0.0f64; 3];
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "rigid3.centroids")?;
        let weight = weights[index];
        weight_total += weight;
        let gain = weight / weight_total;
        let design = scale3(
            sub3(fiducial.design.coords(), design_frame.anchor),
            1.0 / design_frame.scale,
        );
        let measured = scale3(
            sub3(fiducial.measured.coords(), measured_frame.anchor),
            1.0 / measured_frame.scale,
        );
        for axis in 0..3 {
            mean_design[axis] += gain * (design[axis] - mean_design[axis]);
            mean_measured[axis] += gain * (measured[axis] - mean_measured[axis]);
        }
    }
    if !weight_total.is_finite() || weight_total <= 0.0 {
        return Err(Rigid3Error::ArithmeticOverflow {
            field: "total weight",
        });
    }

    // Weighted scatter and cross-covariance accumulation, normalized frame.
    let mut scatter_design = [[0.0f64; 3]; 3];
    let mut scatter_measured = [[0.0f64; 3]; 3];
    let mut cross = [[0.0f64; 3]; 3];
    let mut design_scatter_trace = 0.0f64;
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "rigid3.scatter")?;
        let weight = weights[index];
        let design = sub3(
            scale3(
                sub3(fiducial.design.coords(), design_frame.anchor),
                1.0 / design_frame.scale,
            ),
            mean_design,
        );
        let measured = sub3(
            scale3(
                sub3(fiducial.measured.coords(), measured_frame.anchor),
                1.0 / measured_frame.scale,
            ),
            mean_measured,
        );
        for row in 0..3 {
            for col in 0..3 {
                scatter_design[row][col] =
                    (weight * design[row]).mul_add(design[col], scatter_design[row][col]);
                scatter_measured[row][col] =
                    (weight * measured[row]).mul_add(measured[col], scatter_measured[row][col]);
                cross[row][col] = (weight * design[row]).mul_add(measured[col], cross[row][col]);
            }
        }
        design_scatter_trace = weight.mul_add(dot3(design, design), design_scatter_trace);
    }
    for aggregate in [&scatter_design, &scatter_measured, &cross] {
        if !aggregate.iter().flatten().all(|value| value.is_finite()) {
            return Err(Rigid3Error::ArithmeticOverflow {
                field: "scatter accumulation",
            });
        }
    }

    let design_spectrum = spectrum_of_scatter(&scatter_design);
    let coplanar_design = classify_spectrum(design_spectrum)
        .map_err(|diagnosis| Rigid3Error::DegenerateDesign { diagnosis })?;
    let measured_spectrum = spectrum_of_scatter(&scatter_measured);
    let _coplanar_measured = classify_spectrum(measured_spectrum)
        .map_err(|diagnosis| Rigid3Error::DegenerateMeasured { diagnosis })?;

    // Deterministic one-sided Jacobi SVD of the 3x3 cross-covariance.
    let cross_flat = [
        cross[0][0],
        cross[0][1],
        cross[0][2],
        cross[1][0],
        cross[1][1],
        cross[1][2],
        cross[2][0],
        cross[2][1],
        cross[2][2],
    ];
    let svd = fs_la::factor::svd_jacobi(&cross_flat, 3, 3);
    let sigma = [svd.sigma[0], svd.sigma[1], svd.sigma[2]];
    if sigma[0] <= 0.0 || sigma[1] <= RANK_RELATIVE_TOLERANCE * sigma[0] {
        return Err(Rigid3Error::UnobservableRotation);
    }
    let coplanar_cross = sigma[2] <= RANK_RELATIVE_TOLERANCE * sigma[0];

    // Reflection handling. For a full-rank cross-covariance the sign of
    // det(U) * det(V) says whether the unconstrained orthogonal optimum is a
    // reflection; when it is and the two trailing singular values coincide,
    // the constrained proper optimum is non-unique and we refuse.
    let u1 = column(&svd.u, 0);
    let u2 = column(&svd.u, 1);
    let v1 = column(&svd.v, 0);
    let v2 = column(&svd.v, 1);
    let mut reflection = false;
    if !coplanar_cross {
        let det_u = det3(&[column(&svd.u, 0), column(&svd.u, 1), column(&svd.u, 2)]);
        let det_v = det3(&[column(&svd.v, 0), column(&svd.v, 1), column(&svd.v, 2)]);
        reflection = det_u * det_v < 0.0;
        if reflection && (sigma[1] - sigma[2]) <= RANK_RELATIVE_TOLERANCE * sigma[0] {
            return Err(Rigid3Error::AmbiguousRotation);
        }
    }

    // Right-handed canonicalization: rebuilding each third column as the
    // cross product of the first two forces det(U) = det(V) = +1, after which
    // the Kabsch optimum is exactly V * U^T in both the full-rank and the
    // coplanar rank-2 case — the classical diag(1, 1, det V U^T) correction
    // is absorbed by the canonicalization.
    let u3 = cross3(u1, u2);
    let v3 = cross3(v1, v2);
    let u3_norm = norm3(u3);
    let v3_norm = norm3(v3);
    if u3_norm <= 0.0 || v3_norm <= 0.0 || !u3_norm.is_finite() || !v3_norm.is_finite() {
        return Err(Rigid3Error::UnobservableRotation);
    }
    let u3 = scale3(u3, 1.0 / u3_norm);
    let v3 = scale3(v3, 1.0 / v3_norm);

    // Columns of U and V as rows here; R = V * U^T means
    // R[r][c] = sum_k V[r][k] U[c][k] over the canonical columns.
    let u_columns = [u1, u2, u3];
    let v_columns = [v1, v2, v3];
    let mut rotation = [[0.0f64; 3]; 3];
    for row in 0..3 {
        for col in 0..3 {
            let mut value = 0.0;
            for k in 0..3 {
                value = v_columns[k][row].mul_add(u_columns[k][col], value);
            }
            rotation[row][col] = value;
        }
    }
    if !rotation.iter().flatten().all(|value| value.is_finite()) {
        return Err(Rigid3Error::ArithmeticOverflow { field: "rotation" });
    }

    let signed_sigma_sum = if reflection {
        sigma[0] + sigma[1] - sigma[2]
    } else {
        sigma[0] + sigma[1] + sigma[2]
    };

    let centroid_design = add3(design_frame.anchor, scale3(mean_design, design_frame.scale));
    let centroid_measured = add3(
        measured_frame.anchor,
        scale3(mean_measured, measured_frame.scale),
    );

    Ok(KabschCore {
        rotation,
        centroid_design,
        centroid_measured,
        signed_sigma_sum,
        design_scatter: design_scatter_trace,
        scale_design: design_frame.scale,
        scale_measured: measured_frame.scale,
        condition: RegistrationCondition {
            design_spectrum,
            measured_spectrum,
            cross_singular_values: sigma,
            coplanar_design,
            coplanar_cross,
            reflection_preferred: reflection,
        },
    })
}

// ---------------------------------------------------------------------------
// Public rigid registration.
// ---------------------------------------------------------------------------

/// A rigid 3-D registration (rotation + translation) mapping design →
/// measured, with the advisory residual it carries forward.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rigid3Registration {
    rotation: Mat3,
    translation: Vec3,
    residual_rms: f64,
    condition: RegistrationCondition,
}

impl Rigid3Registration {
    /// Row-major rotation matrix (design → measured).
    #[must_use]
    pub const fn rotation(&self) -> &Mat3 {
        &self.rotation
    }

    /// Translation vector (design → measured).
    #[must_use]
    pub const fn translation(&self) -> Vec3 {
        self.translation
    }

    /// Unweighted root-mean-square fiducial residual. A global fit
    /// diagnostic, not transform covariance.
    #[must_use]
    pub const fn residual_rms(&self) -> f64 {
        self.residual_rms
    }

    /// Spectral conditioning diagnostics of the admitted fit.
    #[must_use]
    pub const fn condition(&self) -> &RegistrationCondition {
        &self.condition
    }

    /// Rotation angle in radians, recovered from the matrix trace.
    #[must_use]
    pub fn rotation_angle_rad(&self) -> f64 {
        rotation_angle_from_matrix(&self.rotation)
    }

    /// Apply the registration to a design point.
    ///
    /// # Errors
    /// Refuses non-finite arithmetic overflow.
    pub fn apply(&self, point: Point3) -> Result<Point3, Rigid3Error> {
        let rotated = mat3_vec(&self.rotation, point.coords());
        let mapped = add3(rotated, self.translation);
        Point3::new(mapped[0], mapped[1], mapped[2]).map_err(|_| Rigid3Error::ArithmeticOverflow {
            field: "applied point",
        })
    }
}

fn residual_rms_for(
    fiducials: &[Fiducial3],
    rotation: &Mat3,
    translation: Vec3,
    scale: f64,
    cx: &fs_exec::Cx<'_>,
) -> Result<f64, Rigid3Error> {
    let mut residuals = ScaledSumSquares::default();
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "rigid3.residual")?;
        let mapped = add3(
            scale3(mat3_vec(rotation, fiducial.design.coords()), scale),
            translation,
        );
        let delta = sub3(fiducial.measured.coords(), mapped);
        let distance = norm3(delta);
        if !distance.is_finite() {
            return Err(Rigid3Error::ArithmeticOverflow {
                field: "residual distance",
            });
        }
        residuals
            .add(distance)
            .map_err(|_| Rigid3Error::ArithmeticOverflow {
                field: "residual sum of squares",
            })?;
    }
    residuals
        .root_mean_square(fiducials.len() as f64)
        .map_err(|_| Rigid3Error::ArithmeticOverflow {
            field: "residual RMS",
        })
}

fn rigid_from_core(
    fiducials: &[Fiducial3],
    core: &KabschCore,
    cx: &fs_exec::Cx<'_>,
) -> Result<Rigid3Registration, Rigid3Error> {
    let translation = sub3(
        core.centroid_measured,
        mat3_vec(&core.rotation, core.centroid_design),
    );
    if !translation.iter().all(|value| value.is_finite()) {
        return Err(Rigid3Error::ArithmeticOverflow {
            field: "translation",
        });
    }
    let residual_rms = residual_rms_for(fiducials, &core.rotation, translation, 1.0, cx)?;
    Ok(Rigid3Registration {
        rotation: core.rotation,
        translation,
        residual_rms,
        condition: core.condition,
    })
}

/// Solve the rigid 3-D registration that best maps the fiducials' design
/// points onto their measured points (closed-form Kabsch via the
/// deterministic Jacobi SVD). Requires at least [`MIN_FIDUCIALS`]
/// non-collinear correspondences; coplanar configurations are admitted and
/// flagged in the condition payload.
///
/// # Errors
/// Degenerate design/measured configurations with a geometric diagnosis,
/// rank-deficient cross-covariance, the ambiguous reflection hard case,
/// non-finite aggregates, or a structured [`Rigid3Error::Cancelled`].
pub fn register3(
    fiducials: &[Fiducial3],
    cx: &fs_exec::Cx<'_>,
) -> Result<Rigid3Registration, Rigid3Error> {
    let mut weights = Vec::new();
    weights
        .try_reserve_exact(fiducials.len())
        .map_err(|_| Rigid3Error::AllocationFailed)?;
    weights.resize(fiducials.len(), 1.0);
    let core = weighted_kabsch(fiducials, &weights, cx)?;
    let registration = rigid_from_core(fiducials, &core, cx)?;
    cx.checkpoint().map_err(|_| Rigid3Error::Cancelled {
        phase: "rigid3.publish",
    })?;
    Ok(registration)
}

// ---------------------------------------------------------------------------
// Similarity registration (scale reported, never silent).
// ---------------------------------------------------------------------------

/// Common unit-conversion ratios named by the units-suspicion diagnostic.
const KNOWN_CONVERSIONS: [(f64, &str); 12] = [
    (25.4, "25.4 (inch to millimetre)"),
    (0.039_370_078_740_157_48, "1/25.4 (millimetre to inch)"),
    (1000.0, "1000 (metre to millimetre)"),
    (0.001, "1/1000 (millimetre to metre)"),
    (100.0, "100 (metre to centimetre)"),
    (0.01, "1/100 (centimetre to metre)"),
    (2.54, "2.54 (inch to centimetre)"),
    (0.393_700_787_401_574_8, "1/2.54 (centimetre to inch)"),
    (0.0254, "0.0254 (inch to metre)"),
    (39.370_078_740_157_48, "39.37 (metre to inch)"),
    (0.3048, "0.3048 (foot to metre)"),
    (3.280_839_895_013_123, "3.281 (metre to foot)"),
];

/// A named unit-conversion factor near a suspicious scale estimate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnitSuspicion {
    conversion_name: &'static str,
    factor: f64,
    log_gap: f64,
}

impl UnitSuspicion {
    /// Stable description of the nearest known conversion.
    #[must_use]
    pub const fn conversion_name(&self) -> &'static str {
        self.conversion_name
    }

    /// The nearest known conversion ratio.
    #[must_use]
    pub const fn factor(&self) -> f64 {
        self.factor
    }

    /// `|ln(scale / factor)|`: how close the estimate sits to that ratio.
    #[must_use]
    pub const fn log_gap(&self) -> f64 {
        self.log_gap
    }
}

/// The similarity fit's scale estimate with its first-order uncertainty and
/// units-suspicion diagnostic.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScaleAssessment {
    estimate: f64,
    standard_error: f64,
    declared_tolerance: f64,
    suspicion: Option<UnitSuspicion>,
}

impl ScaleAssessment {
    /// Umeyama scale estimate (design → measured).
    #[must_use]
    pub const fn estimate(&self) -> f64 {
        self.estimate
    }

    /// First-order standard error of the scale under an isotropic
    /// homoscedastic noise model estimated from the fit residuals. A
    /// reported diagnostic, not a calibrated bound.
    #[must_use]
    pub const fn standard_error(&self) -> f64 {
        self.standard_error
    }

    /// The caller-declared relative tolerance around 1.0.
    #[must_use]
    pub const fn declared_tolerance(&self) -> f64 {
        self.declared_tolerance
    }

    /// Fires when `|estimate - 1|` exceeds the declared tolerance; names the
    /// nearest common unit-conversion ratio so a unit error is investigated
    /// before the scale is believed.
    #[must_use]
    pub const fn suspicion(&self) -> Option<&UnitSuspicion> {
        self.suspicion.as_ref()
    }
}

/// A similarity 3-D registration (scale, rotation, translation) mapping
/// design → measured, with its scale assessment and advisory residual.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Similarity3Registration {
    rotation: Mat3,
    translation: Vec3,
    scale: ScaleAssessment,
    residual_rms: f64,
    condition: RegistrationCondition,
}

impl Similarity3Registration {
    /// Row-major rotation matrix (design → measured).
    #[must_use]
    pub const fn rotation(&self) -> &Mat3 {
        &self.rotation
    }

    /// Translation vector (design → measured).
    #[must_use]
    pub const fn translation(&self) -> Vec3 {
        self.translation
    }

    /// Scale estimate, uncertainty, and units-suspicion diagnostic.
    #[must_use]
    pub const fn scale(&self) -> &ScaleAssessment {
        &self.scale
    }

    /// Unweighted root-mean-square fiducial residual of the similarity fit.
    #[must_use]
    pub const fn residual_rms(&self) -> f64 {
        self.residual_rms
    }

    /// Spectral conditioning diagnostics of the admitted fit.
    #[must_use]
    pub const fn condition(&self) -> &RegistrationCondition {
        &self.condition
    }

    /// Rotation angle in radians, recovered from the matrix trace.
    #[must_use]
    pub fn rotation_angle_rad(&self) -> f64 {
        rotation_angle_from_matrix(&self.rotation)
    }

    /// Apply the similarity to a design point.
    ///
    /// # Errors
    /// Refuses non-finite arithmetic overflow.
    pub fn apply(&self, point: Point3) -> Result<Point3, Rigid3Error> {
        let mapped = add3(
            scale3(
                mat3_vec(&self.rotation, point.coords()),
                self.scale.estimate,
            ),
            self.translation,
        );
        Point3::new(mapped[0], mapped[1], mapped[2]).map_err(|_| Rigid3Error::ArithmeticOverflow {
            field: "applied point",
        })
    }
}

/// Solve the similarity 3-D registration (rotation + translation + one
/// global scale). The scale is REPORTED with a first-order standard error and
/// a units-suspicion diagnostic: `scale_tolerance` is the caller-declared
/// relative band around 1.0 outside which suspicion fires. There is no
/// default tolerance; silent rescaling is the failure mode this API exists to
/// prevent.
///
/// # Errors
/// Everything [`register3`] refuses, plus an invalid tolerance or a
/// non-positive/non-finite scale aggregate.
pub fn register3_similarity(
    fiducials: &[Fiducial3],
    scale_tolerance: f64,
    cx: &fs_exec::Cx<'_>,
) -> Result<Similarity3Registration, Rigid3Error> {
    if !scale_tolerance.is_finite() || scale_tolerance < 0.0 {
        return Err(Rigid3Error::InvalidScalar {
            field: "scale_tolerance",
            requirement: "finite and non-negative",
        });
    }
    let mut weights = Vec::new();
    weights
        .try_reserve_exact(fiducials.len())
        .map_err(|_| Rigid3Error::AllocationFailed)?;
    weights.resize(fiducials.len(), 1.0);
    let core = weighted_kabsch(fiducials, &weights, cx)?;

    if core.design_scatter <= 0.0 {
        return Err(Rigid3Error::DegenerateDesign {
            diagnosis: DegeneracyDiagnosis::Coincident,
        });
    }
    // Un-normalize the Umeyama ratio: cross-covariance entries carry
    // scale_design * scale_measured, the design scatter carries
    // scale_design^2.
    let scale_estimate =
        core.signed_sigma_sum * core.scale_measured / (core.design_scatter * core.scale_design);
    if !scale_estimate.is_finite() || scale_estimate <= 0.0 {
        return Err(Rigid3Error::InvalidScalar {
            field: "scale estimate",
            requirement: "finite and positive",
        });
    }

    let translation = sub3(
        core.centroid_measured,
        scale3(
            mat3_vec(&core.rotation, core.centroid_design),
            scale_estimate,
        ),
    );
    if !translation.iter().all(|value| value.is_finite()) {
        return Err(Rigid3Error::ArithmeticOverflow {
            field: "translation",
        });
    }
    let residual_rms =
        residual_rms_for(fiducials, &core.rotation, translation, scale_estimate, cx)?;

    // First-order scale standard error under isotropic homoscedastic noise:
    // var(s) ~= sigma^2 / sum w ||dp||^2, with sigma^2 estimated from the
    // similarity residuals at 3n - 7 degrees of freedom.
    let n = fiducials.len() as f64;
    let dof = 3.0f64.mul_add(n, -7.0).max(1.0);
    let noise_variance = residual_rms * residual_rms * n / dof;
    let design_scatter_unnormalized = core.design_scatter * core.scale_design * core.scale_design;
    let standard_error = if design_scatter_unnormalized > 0.0 {
        (noise_variance / design_scatter_unnormalized).sqrt()
    } else {
        return Err(Rigid3Error::DegenerateDesign {
            diagnosis: DegeneracyDiagnosis::Coincident,
        });
    };
    if !standard_error.is_finite() {
        return Err(Rigid3Error::ArithmeticOverflow {
            field: "scale standard error",
        });
    }

    let suspicion = if (scale_estimate - 1.0).abs() > scale_tolerance {
        let mut best: Option<UnitSuspicion> = None;
        for (factor, conversion_name) in KNOWN_CONVERSIONS {
            let log_gap = (scale_estimate / factor).ln().abs();
            let better = match &best {
                None => true,
                Some(current) => log_gap < current.log_gap,
            };
            if better {
                best = Some(UnitSuspicion {
                    conversion_name,
                    factor,
                    log_gap,
                });
            }
        }
        best
    } else {
        None
    };

    cx.checkpoint().map_err(|_| Rigid3Error::Cancelled {
        phase: "rigid3.publish",
    })?;
    Ok(Similarity3Registration {
        rotation: core.rotation,
        translation,
        scale: ScaleAssessment {
            estimate: scale_estimate,
            standard_error,
            declared_tolerance: scale_tolerance,
            suspicion,
        },
        residual_rms,
        condition: core.condition,
    })
}

// ---------------------------------------------------------------------------
// Calibrated 6-dof covariance.
// ---------------------------------------------------------------------------

/// A symmetric, strictly positive-definite 3x3 measurement covariance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Covariance3 {
    xx: f64,
    xy: f64,
    xz: f64,
    yy: f64,
    yz: f64,
    zz: f64,
}

impl Covariance3 {
    /// Construct a finite covariance. Strict positive definiteness is checked
    /// after scale normalization so tiny and huge valid matrices are treated
    /// consistently.
    ///
    /// # Errors
    /// Non-finite entries, non-positive marginal variances, or a normalized
    /// Cholesky pivot at or below zero.
    pub fn new(xx: f64, xy: f64, xz: f64, yy: f64, yz: f64, zz: f64) -> Result<Self, Rigid3Error> {
        for (field, value) in [
            ("covariance.xx", xx),
            ("covariance.xy", xy),
            ("covariance.xz", xz),
            ("covariance.yy", yy),
            ("covariance.yz", yz),
            ("covariance.zz", zz),
        ] {
            if !value.is_finite() {
                return Err(Rigid3Error::InvalidScalar {
                    field,
                    requirement: "finite",
                });
            }
        }
        let scale = xx
            .abs()
            .max(xy.abs())
            .max(xz.abs())
            .max(yy.abs())
            .max(yz.abs())
            .max(zz.abs());
        if xx <= 0.0 || yy <= 0.0 || zz <= 0.0 || scale == 0.0 {
            return Err(Rigid3Error::NonPositiveDefiniteCovariance { index: 0 });
        }
        let normalized = [
            [xx / scale, xy / scale, xz / scale],
            [xy / scale, yy / scale, yz / scale],
            [xz / scale, yz / scale, zz / scale],
        ];
        // Strict normalized Cholesky admission: every pivot must be finite
        // and strictly positive.
        let mut lower = [[0.0f64; 3]; 3];
        for row in 0..3 {
            for column in 0..=row {
                let mut value = normalized[row][column];
                for prior in 0..column {
                    value -= lower[row][prior] * lower[column][prior];
                }
                if row == column {
                    if !value.is_finite() || value <= 0.0 {
                        return Err(Rigid3Error::NonPositiveDefiniteCovariance { index: 0 });
                    }
                    lower[row][column] = value.sqrt();
                } else {
                    lower[row][column] = value / lower[column][column];
                }
            }
        }
        Ok(Self {
            xx,
            xy,
            xz,
            yy,
            yz,
            zz,
        })
    }

    /// x variance.
    #[must_use]
    pub const fn xx(self) -> f64 {
        self.xx
    }

    /// x/y covariance.
    #[must_use]
    pub const fn xy(self) -> f64 {
        self.xy
    }

    /// x/z covariance.
    #[must_use]
    pub const fn xz(self) -> f64 {
        self.xz
    }

    /// y variance.
    #[must_use]
    pub const fn yy(self) -> f64 {
        self.yy
    }

    /// y/z covariance.
    #[must_use]
    pub const fn yz(self) -> f64 {
        self.yz
    }

    /// z variance.
    #[must_use]
    pub const fn zz(self) -> f64 {
        self.zz
    }

    /// Trace.
    #[must_use]
    pub fn trace(self) -> f64 {
        self.xx + self.yy + self.zz
    }

    fn as_matrix(self) -> Mat3 {
        [
            [self.xx, self.xy, self.xz],
            [self.xy, self.yy, self.yz],
            [self.xz, self.yz, self.zz],
        ]
    }

    /// Symmetric principal inverse square root through the deterministic
    /// Jacobi eigendecomposition. The symmetric factor keeps standardized
    /// residual norms equivariant under rigid coordinate-frame rotations.
    fn principal_inverse_sqrt(self, index: usize) -> Result<Mat3, Rigid3Error> {
        let matrix = self.as_matrix();
        let flat = [
            matrix[0][0],
            matrix[0][1],
            matrix[0][2],
            matrix[1][0],
            matrix[1][1],
            matrix[1][2],
            matrix[2][0],
            matrix[2][1],
            matrix[2][2],
        ];
        let (eigenvalues, eigenvectors) = fs_la::eigen::jacobi_eigh(&flat, 3);
        if eigenvalues[0] <= 0.0 || !eigenvalues.iter().all(|value| value.is_finite()) {
            return Err(Rigid3Error::NonPositiveDefiniteCovariance { index });
        }
        let mut result = [[0.0f64; 3]; 3];
        for k in 0..3 {
            let inv_root = 1.0 / eigenvalues[k].sqrt();
            let q = [eigenvectors[k], eigenvectors[3 + k], eigenvectors[6 + k]];
            for row in 0..3 {
                for column in 0..3 {
                    result[row][column] =
                        (inv_root * q[row]).mul_add(q[column], result[row][column]);
                }
            }
        }
        if result.iter().flatten().all(|value| value.is_finite()) {
            Ok(result)
        } else {
            Err(Rigid3Error::NonPositiveDefiniteCovariance { index })
        }
    }
}

/// Declared correlation between distinct 3-D fiducial measurement errors.
///
/// The 2-D module's standardized equicorrelation shortcut is deliberately not
/// carried into 3-D yet; only independence is quantifiable here, and unknown
/// dependence refuses rather than silently assuming independence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossFiducialModel3 {
    /// Distinct fiducial errors are independent; each may still carry a full
    /// 3x3 within-point covariance.
    Independent,
    /// Dependence exists but no cross-covariance model was supplied.
    /// Estimation refuses rather than silently using independence.
    Unknown,
}

/// Complete calibration/noise model for one 3-D fiducial family.
#[derive(Debug, Clone, PartialEq)]
pub struct MetrologyModel3 {
    fiducial_covariances: Vec<Covariance3>,
    cross_fiducial: CrossFiducialModel3,
    huber: HuberPolicy,
    calibration_identity: String,
}

impl MetrologyModel3 {
    /// Construct a model whose per-fiducial covariance order matches the
    /// correspondence order supplied to [`estimate_calibrated_rigid3`].
    ///
    /// # Errors
    /// Invalid sizes, Huber domain, or calibration identity grammar.
    pub fn new(
        fiducial_covariances: Vec<Covariance3>,
        cross_fiducial: CrossFiducialModel3,
        huber: HuberPolicy,
        calibration_identity: impl Into<String>,
    ) -> Result<Self, Rigid3Error> {
        let count = fiducial_covariances.len();
        if count < MIN_FIDUCIALS {
            return Err(Rigid3Error::TooFewFiducials {
                have: count,
                need: MIN_FIDUCIALS,
            });
        }
        if count > MAX_AS_BUILT_POINTS {
            return Err(Rigid3Error::TooManyPoints {
                have: count,
                max: MAX_AS_BUILT_POINTS,
            });
        }
        if let HuberPolicy::Enabled {
            threshold,
            iterations,
        } = huber
        {
            if !threshold.is_finite() || threshold <= 0.0 {
                return Err(Rigid3Error::InvalidScalar {
                    field: "huber.threshold",
                    requirement: "finite and positive",
                });
            }
            if iterations == 0 || iterations > crate::uncertainty::MAX_HUBER_ITERATIONS {
                return Err(Rigid3Error::InvalidScalar {
                    field: "huber.iterations",
                    requirement: "between 1 and MAX_HUBER_ITERATIONS",
                });
            }
        }
        let calibration_identity = calibration_identity.into();
        if let Some(reason) = fs_evidence::color_leaf_identity_reason(&calibration_identity) {
            return Err(Rigid3Error::InvalidCalibrationIdentity { reason });
        }
        Ok(Self {
            fiducial_covariances,
            cross_fiducial,
            huber,
            calibration_identity,
        })
    }

    /// Ordered calibrated fiducial covariances.
    #[must_use]
    pub fn fiducial_covariances(&self) -> &[Covariance3] {
        &self.fiducial_covariances
    }

    /// Declared cross-fiducial model.
    #[must_use]
    pub const fn cross_fiducial(&self) -> CrossFiducialModel3 {
        self.cross_fiducial
    }

    /// Robust policy.
    #[must_use]
    pub const fn huber(&self) -> HuberPolicy {
        self.huber
    }

    /// Calibration artifact identity. Bound into the model root but not
    /// self-authenticating.
    #[must_use]
    pub fn calibration_identity(&self) -> &str {
        &self.calibration_identity
    }
}

/// Per-fiducial standardized residual and robust classification for the 3-D
/// calibrated fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rigid3OutlierDiagnostic {
    standardized_residual_norm: f64,
    robust_weight: f64,
    disposition: crate::uncertainty::OutlierDisposition,
}

impl Rigid3OutlierDiagnostic {
    /// Residual norm after whitening by the declared covariance's symmetric
    /// principal inverse square root.
    #[must_use]
    pub const fn standardized_residual_norm(&self) -> f64 {
        self.standardized_residual_norm
    }

    /// Final Huber multiplier used for the published point estimate.
    #[must_use]
    pub const fn robust_weight(&self) -> f64 {
        self.robust_weight
    }

    /// Stable robust-screen classification.
    #[must_use]
    pub const fn disposition(&self) -> crate::uncertainty::OutlierDisposition {
        self.disposition
    }
}

/// Robustly refitted 3-D registration plus the full first-order sandwich
/// covariance of `(tx, ty, tz, rx, ry, rz)`. For Huber fits the covariance is
/// a frozen-weight sandwich: conditional, first-order, and silent about
/// weight-selection uncertainty.
#[derive(Debug, Clone, PartialEq)]
pub struct CalibratedRigid3Registration {
    registration: Rigid3Registration,
    covariance: Mat6,
    degrees_of_freedom: usize,
    robust_weights: Vec<f64>,
    leverage: Vec<f64>,
    outlier_diagnostics: Vec<Rigid3OutlierDiagnostic>,
    model_identity: fs_blake3::ContentHash,
    robust_conditional: bool,
}

impl CalibratedRigid3Registration {
    /// Robust weighted-Kabsch point estimate.
    #[must_use]
    pub const fn registration(&self) -> &Rigid3Registration {
        &self.registration
    }

    /// Row-major covariance for `(tx, ty, tz, rx, ry, rz)`.
    #[must_use]
    pub const fn covariance(&self) -> &PoseCovariance6 {
        &self.covariance
    }

    /// Residual degrees of freedom, exactly `3n - 6`.
    #[must_use]
    pub const fn degrees_of_freedom(&self) -> usize {
        self.degrees_of_freedom
    }

    /// Final deterministic per-fiducial Huber multipliers. The complete
    /// estimator weight is this multiplier times the model base weight
    /// `3 / trace(Sigma_i)`.
    #[must_use]
    pub fn robust_weights(&self) -> &[f64] {
        &self.robust_weights
    }

    /// Per-fiducial hat-block leverage traces; their sum is the fitted
    /// parameter dimension (6).
    #[must_use]
    pub fn leverage(&self) -> &[f64] {
        &self.leverage
    }

    /// Per-fiducial standardized residual and robust classification.
    #[must_use]
    pub fn outlier_diagnostics(&self) -> &[Rigid3OutlierDiagnostic] {
        &self.outlier_diagnostics
    }

    /// Domain-separated identity binding model, inputs, fit, covariance, and
    /// diagnostics. An integrity address, not authentication.
    #[must_use]
    pub const fn model_identity(&self) -> fs_blake3::ContentHash {
        self.model_identity
    }

    /// Whether the covariance is conditional on adaptive robust weights.
    #[must_use]
    pub const fn robust_conditional(&self) -> bool {
        self.robust_conditional
    }
}

/// The 3x6 pose Jacobian block for one fiducial: `[I3 | -[y]x]` with
/// `y = R (p - centroid)`, the left rotation-vector perturbation about the
/// weighted design centroid image.
fn pose_jacobian(y: Vec3) -> Mat36 {
    [
        [1.0, 0.0, 0.0, 0.0, y[2], -y[1]],
        [0.0, 1.0, 0.0, -y[2], 0.0, y[0]],
        [0.0, 0.0, 1.0, y[1], -y[0], 0.0],
    ]
}

fn residual_vector(
    fiducial: Fiducial3,
    rotation: &Mat3,
    translation: Vec3,
) -> Result<Vec3, Rigid3Error> {
    let mapped = add3(mat3_vec(rotation, fiducial.design.coords()), translation);
    let residual = sub3(fiducial.measured.coords(), mapped);
    if residual.iter().all(|value| value.is_finite()) {
        Ok(residual)
    } else {
        Err(Rigid3Error::ArithmeticOverflow {
            field: "residual vector",
        })
    }
}

/// Estimate a robust 3-D rigid transform and its complete first-order
/// fixed-weight pose covariance. The point estimate is the closed-form
/// weighted Kabsch fit — the exact global minimizer of the declared
/// scalar-weighted objective — with base weights `3 / trace(Sigma_i)` and
/// deterministic Huber multipliers refreshed against declared-covariance
/// standardized residual norms, re-solving after every refresh including the
/// last. The covariance is the sandwich for that estimator under the full
/// declared per-fiducial covariance model; for isotropic models it reduces
/// exactly to the classical generalized-least-squares covariance.
///
/// # Errors
/// Invalid lengths/domains, unknown dependence, degenerate geometry,
/// unresolved information, cancellation, allocation failure, or non-finite
/// arithmetic.
#[allow(clippy::too_many_lines)]
pub fn estimate_calibrated_rigid3(
    fiducials: &[Fiducial3],
    model: &MetrologyModel3,
    cx: &fs_exec::Cx<'_>,
) -> Result<CalibratedRigid3Registration, Rigid3Error> {
    if fiducials.len() != model.fiducial_covariances.len() {
        return Err(Rigid3Error::LengthMismatch {
            field: "fiducial_covariances",
            expected: fiducials.len(),
            found: model.fiducial_covariances.len(),
        });
    }
    if matches!(model.cross_fiducial, CrossFiducialModel3::Unknown) {
        return Err(Rigid3Error::UnknownDependence);
    }

    // Whitening factors and base weights from the declared covariances.
    let count = fiducials.len();
    let mut whiteners = Vec::new();
    whiteners
        .try_reserve_exact(count)
        .map_err(|_| Rigid3Error::AllocationFailed)?;
    let mut base_weights = Vec::new();
    base_weights
        .try_reserve_exact(count)
        .map_err(|_| Rigid3Error::AllocationFailed)?;
    for (index, covariance) in model.fiducial_covariances.iter().enumerate() {
        checkpoint(cx, index, "rigid3.whitening")?;
        whiteners.push(covariance.principal_inverse_sqrt(index)?);
        let trace = covariance.trace();
        if !trace.is_finite() || trace <= 0.0 {
            return Err(Rigid3Error::NonPositiveDefiniteCovariance { index });
        }
        base_weights.push(3.0 / trace);
    }

    let mut huber_multipliers = Vec::new();
    huber_multipliers
        .try_reserve_exact(count)
        .map_err(|_| Rigid3Error::AllocationFailed)?;
    huber_multipliers.resize(count, 1.0);

    let mut estimator_weights = Vec::new();
    estimator_weights
        .try_reserve_exact(count)
        .map_err(|_| Rigid3Error::AllocationFailed)?;
    estimator_weights.extend_from_slice(&base_weights);

    let mut core = weighted_kabsch(fiducials, &estimator_weights, cx)?;
    let mut registration = rigid_from_core(fiducials, &core, cx)?;

    let standardized_norms = |rotation: &Mat3,
                              translation: Vec3,
                              phase: &'static str|
     -> Result<Vec<f64>, Rigid3Error> {
        let mut norms = Vec::new();
        norms
            .try_reserve_exact(count)
            .map_err(|_| Rigid3Error::AllocationFailed)?;
        for (index, fiducial) in fiducials.iter().enumerate() {
            checkpoint(cx, index, phase)?;
            let residual = residual_vector(*fiducial, rotation, translation)?;
            let whitened = mat3_vec(&whiteners[index], residual);
            let norm = norm3(whitened);
            if !norm.is_finite() {
                return Err(Rigid3Error::ArithmeticOverflow {
                    field: "standardized residual norm",
                });
            }
            norms.push(norm);
        }
        Ok(norms)
    };

    if let HuberPolicy::Enabled {
        threshold,
        iterations,
    } = model.huber
    {
        for _ in 0..iterations {
            let norms = standardized_norms(
                &registration.rotation,
                registration.translation,
                "rigid3.robust-norms",
            )?;
            for (index, norm) in norms.iter().enumerate() {
                checkpoint(cx, index, "rigid3.robust-weights")?;
                huber_multipliers[index] = if *norm <= threshold || *norm == 0.0 {
                    1.0
                } else {
                    threshold / norm
                };
                estimator_weights[index] = base_weights[index] * huber_multipliers[index];
            }
            // Re-solve after every refresh, including the last: the
            // published transform and covariance use identical weights.
            core = weighted_kabsch(fiducials, &estimator_weights, cx)?;
            registration = rigid_from_core(fiducials, &core, cx)?;
        }
    }

    // Weighted design centroid for the rotation pivot: the same weights the
    // final solve used, so the translation/rotation cross information
    // vanishes in exact arithmetic.
    let mut weight_total = 0.0f64;
    let mut centroid = [0.0f64; 3];
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "rigid3.covariance-centroid")?;
        let weight = estimator_weights[index];
        weight_total += weight;
        let gain = weight / weight_total;
        let design = fiducial.design.coords();
        for axis in 0..3 {
            centroid[axis] += gain * (design[axis] - centroid[axis]);
        }
    }

    // Information (bread) and sandwich meat under the declared model.
    let mut information = [[0.0f64; 6]; 6];
    let mut meat = [[0.0f64; 6]; 6];
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "rigid3.information")?;
        let weight = estimator_weights[index];
        let pivoted = sub3(fiducial.design.coords(), centroid);
        let y = mat3_vec(&registration.rotation, pivoted);
        let jacobian = pose_jacobian(y);
        let sigma = model.fiducial_covariances[index].as_matrix();
        // sigma_j = Sigma_i * J (3x6).
        let mut sigma_j = [[0.0f64; 6]; 3];
        for row in 0..3 {
            for column in 0..6 {
                let mut value = 0.0;
                for inner in 0..3 {
                    value = sigma[row][inner].mul_add(jacobian[inner][column], value);
                }
                sigma_j[row][column] = value;
            }
        }
        let weight_squared = weight * weight;
        for row in 0..6 {
            for column in 0..6 {
                let mut jt_j = 0.0;
                let mut jt_sigma_j = 0.0;
                for inner in 0..3 {
                    jt_j = jacobian[inner][row].mul_add(jacobian[inner][column], jt_j);
                    jt_sigma_j = jacobian[inner][row].mul_add(sigma_j[inner][column], jt_sigma_j);
                }
                information[row][column] = (weight * jt_j).mul_add(1.0, information[row][column]);
                meat[row][column] = (weight_squared * jt_sigma_j).mul_add(1.0, meat[row][column]);
            }
        }
    }

    let bread = invert_spd::<6>(&information)?;
    let mut raw_covariance = [[0.0f64; 6]; 6];
    // raw = bread * meat * bread^T (bread is symmetric).
    let mut meat_bread = [[0.0f64; 6]; 6];
    for row in 0..6 {
        for column in 0..6 {
            let mut value = 0.0;
            for inner in 0..6 {
                value = meat[row][inner].mul_add(bread[inner][column], value);
            }
            meat_bread[row][column] = value;
        }
    }
    for row in 0..6 {
        for column in 0..6 {
            let mut value = 0.0;
            for inner in 0..6 {
                value = bread[row][inner].mul_add(meat_bread[inner][column], value);
            }
            raw_covariance[row][column] = value;
        }
    }
    let covariance = symmetrize_and_validate_spd::<6>(raw_covariance)?;

    // Leverage: w_i * trace(J_i * bread * J_i^T); the traces sum to 6.
    let mut leverage = Vec::new();
    leverage
        .try_reserve_exact(count)
        .map_err(|_| Rigid3Error::AllocationFailed)?;
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "rigid3.leverage")?;
        let pivoted = sub3(fiducial.design.coords(), centroid);
        let y = mat3_vec(&registration.rotation, pivoted);
        let jacobian = pose_jacobian(y);
        let mut trace = 0.0f64;
        for row in 0..3 {
            for a in 0..6 {
                let mut bread_jt = 0.0;
                for b in 0..6 {
                    bread_jt = bread[a][b].mul_add(jacobian[row][b], bread_jt);
                }
                trace = jacobian[row][a].mul_add(bread_jt, trace);
            }
        }
        let value = estimator_weights[index] * trace;
        if !value.is_finite() {
            return Err(Rigid3Error::ArithmeticOverflow { field: "leverage" });
        }
        leverage.push(value);
    }

    let final_norms = standardized_norms(
        &registration.rotation,
        registration.translation,
        "rigid3.final-norms",
    )?;
    let mut outlier_diagnostics = Vec::new();
    outlier_diagnostics
        .try_reserve_exact(count)
        .map_err(|_| Rigid3Error::AllocationFailed)?;
    for (index, norm) in final_norms.iter().enumerate() {
        checkpoint(cx, index, "rigid3.outlier-diagnostics")?;
        let multiplier = huber_multipliers[index];
        let disposition = match model.huber {
            HuberPolicy::Disabled => crate::uncertainty::OutlierDisposition::NotEvaluated,
            HuberPolicy::Enabled { .. } if multiplier < 1.0 => {
                crate::uncertainty::OutlierDisposition::Downweighted
            }
            HuberPolicy::Enabled { .. } => crate::uncertainty::OutlierDisposition::Retained,
        };
        outlier_diagnostics.push(Rigid3OutlierDiagnostic {
            standardized_residual_norm: *norm,
            robust_weight: multiplier,
            disposition,
        });
    }

    let model_identity = calibrated_identity(
        fiducials,
        model,
        &registration,
        &covariance,
        &huber_multipliers,
        &leverage,
        &outlier_diagnostics,
        cx,
    )?;
    cx.checkpoint().map_err(|_| Rigid3Error::Cancelled {
        phase: "rigid3.calibrated-publish",
    })?;
    Ok(CalibratedRigid3Registration {
        registration,
        covariance,
        degrees_of_freedom: count * 3 - 6,
        robust_weights: huber_multipliers,
        leverage,
        outlier_diagnostics,
        model_identity,
        robust_conditional: !matches!(model.huber, HuberPolicy::Disabled),
    })
}

struct IdentityEncoder {
    hasher: fs_blake3::Blake3,
}

impl IdentityEncoder {
    fn new() -> Self {
        let mut encoder = Self {
            hasher: fs_blake3::Blake3::new(),
        };
        encoder.bytes(IDENTITY_DOMAIN.as_bytes());
        encoder.bytes(IDENTITY_KIND);
        encoder
    }

    fn bytes(&mut self, bytes: &[u8]) {
        self.hasher.update(&(bytes.len() as u64).to_le_bytes());
        self.hasher.update(bytes);
    }

    fn u8(&mut self, value: u8) {
        self.bytes(&[value]);
    }

    fn u32(&mut self, value: u32) {
        self.bytes(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.bytes(&value.to_le_bytes());
    }

    fn f64(&mut self, value: f64) {
        let value = canonical_zero(value);
        self.bytes(&value.to_bits().to_le_bytes());
    }

    fn point(&mut self, point: Point3) {
        self.f64(point.x());
        self.f64(point.y());
        self.f64(point.z());
    }

    fn finish(self) -> fs_blake3::ContentHash {
        let preimage = self.hasher.finalize();
        fs_blake3::hash_domain(IDENTITY_DOMAIN, preimage.as_bytes())
    }
}

#[allow(clippy::too_many_arguments)]
fn calibrated_identity(
    fiducials: &[Fiducial3],
    model: &MetrologyModel3,
    registration: &Rigid3Registration,
    covariance: &Mat6,
    huber_multipliers: &[f64],
    leverage: &[f64],
    outlier_diagnostics: &[Rigid3OutlierDiagnostic],
    cx: &fs_exec::Cx<'_>,
) -> Result<fs_blake3::ContentHash, Rigid3Error> {
    let mut encoder = IdentityEncoder::new();
    encoder.u32(RIGID3_SCHEMA_VERSION);
    encoder.u64(fiducials.len() as u64);
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "rigid3.identity")?;
        encoder.point(fiducial.design());
        encoder.point(fiducial.measured());
        let covariance3 = model.fiducial_covariances[index];
        encoder.f64(covariance3.xx());
        encoder.f64(covariance3.xy());
        encoder.f64(covariance3.xz());
        encoder.f64(covariance3.yy());
        encoder.f64(covariance3.yz());
        encoder.f64(covariance3.zz());
    }
    encoder.u8(match model.cross_fiducial {
        CrossFiducialModel3::Independent => 0,
        CrossFiducialModel3::Unknown => 1,
    });
    match model.huber {
        HuberPolicy::Disabled => {
            encoder.u8(0);
        }
        HuberPolicy::Enabled {
            threshold,
            iterations,
        } => {
            encoder.u8(1);
            encoder.f64(threshold);
            encoder.u8(iterations);
        }
    }
    encoder.bytes(model.calibration_identity.as_bytes());
    for row in registration.rotation() {
        for value in row {
            encoder.f64(*value);
        }
    }
    for value in registration.translation() {
        encoder.f64(value);
    }
    encoder.f64(registration.residual_rms());
    let condition = registration.condition();
    for value in condition.design_spectrum() {
        encoder.f64(value);
    }
    for value in condition.measured_spectrum() {
        encoder.f64(value);
    }
    for value in condition.cross_singular_values() {
        encoder.f64(value);
    }
    encoder.u8(u8::from(condition.coplanar_design()));
    encoder.u8(u8::from(condition.coplanar_cross()));
    encoder.u8(u8::from(condition.reflection_preferred()));
    for row in covariance {
        for value in row {
            encoder.f64(*value);
        }
    }
    for (index, multiplier) in huber_multipliers.iter().enumerate() {
        checkpoint(cx, index, "rigid3.identity-weights")?;
        encoder.f64(*multiplier);
        encoder.f64(leverage[index]);
        let diagnostic = &outlier_diagnostics[index];
        encoder.f64(diagnostic.standardized_residual_norm());
        encoder.u8(match diagnostic.disposition() {
            crate::uncertainty::OutlierDisposition::NotEvaluated => 0,
            crate::uncertainty::OutlierDisposition::Retained => 1,
            crate::uncertainty::OutlierDisposition::Downweighted => 2,
        });
    }
    Ok(encoder.finish())
}
