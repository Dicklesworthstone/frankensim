//! Calibrated spatial uncertainty for 2-D as-built registration.
//!
//! This module is deliberately separate from the legacy residual-RMS screens
//! in the crate root. A registration residual is a fit diagnostic; it is not a
//! transform covariance and is never added to the pointwise bounds here.
//!
//! The implemented model is a first-order, calibrated-noise model. Each
//! fiducial carries a positive-definite 2x2 measurement covariance. Fiducials
//! may be independent or may share one declared equicorrelation after their
//! individual covariance factors whiten them. Deterministic Huber reweighting
//! is available only for the independent model; combining data-dependent
//! weights with the equicorrelation shortcut refuses rather than silently
//! claiming generalized-least-squares coverage.
//!
//! Spatial bounds use the full translation/rotation covariance and a
//! distribution-free Chebyshev plus union bound. The result is simultaneous
//! over the complete inspected family and requires neither independence nor a
//! Gaussian claim for inspection points. Its coverage remains conditional on
//! the supplied covariance and bias bounds being calibrated upper models and
//! on inspection data being disjoint from registration data. Unknown overlap
//! or an unavailable bound on the total registered-inspection bias yields an
//! explicit indeterminate no-claim result.

#![allow(clippy::needless_range_loop)] // Fixed 2x2/3x3/4x4 indices expose the parameter ordering.
#![allow(clippy::float_cmp)] // Exact zeros distinguish structural trust-region and IEEE cases.

use super::{Fiducial, MAX_AS_BUILT_POINTS, Point2, Registration};
use fs_ivl::Interval;

/// Identity schema for calibrated registration models.
pub const REGISTRATION_UNCERTAINTY_SCHEMA_VERSION: u32 = 1;
/// Identity schema for simultaneous as-built decisions.
pub const SPATIAL_EVIDENCE_SCHEMA_VERSION: u32 = 1;
/// Maximum deterministic Huber reweighting passes.
pub const MAX_HUBER_ITERATIONS: u8 = 32;

const IDENTITY_DOMAIN: &str = "org.frankensim.fs-asbuilt.spatial-uncertainty.v1";
const DENY_POLICY_DOMAIN: &str = "org.frankensim.fs-asbuilt.spatial-evidence.deny.v1";
const POLL_STRIDE: usize = 256;

/// A malformed or scientifically unsupported spatial-uncertainty request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpatialUncertaintyError {
    /// Parallel arrays differ in length.
    LengthMismatch {
        /// Array being checked.
        field: &'static str,
        /// Required length.
        expected: usize,
        /// Supplied length.
        found: usize,
    },
    /// Too few fiducials were supplied for three rigid parameters.
    TooFewFiducials {
        /// Supplied fiducial count.
        have: usize,
    },
    /// A bounded point family exceeds the crate-wide cap.
    TooManyPoints {
        /// Supplied point count.
        have: usize,
        /// Maximum accepted count.
        max: usize,
    },
    /// A scalar violates a named finite/range invariant.
    InvalidScalar {
        /// Stable field name.
        field: &'static str,
        /// Stable required domain.
        requirement: &'static str,
    },
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
    /// Robust weights cannot be combined with the selected correlation model.
    RobustCorrelationUnsupported,
    /// Cross-fiducial dependence was declared but not quantified.
    UnknownDependence,
    /// The rigid information matrix is singular or numerically unresolved.
    SingularInformation,
    /// More than one globally minimizing rigid rotation is compatible with the
    /// calibrated fixed-weight objective.
    AmbiguousGlobalMinimum,
    /// A finite input produced an unrepresentable finite-model aggregate.
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

impl core::fmt::Display for SpatialUncertaintyError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::LengthMismatch {
                field,
                expected,
                found,
            } => write!(formatter, "{field} has length {found}; expected {expected}"),
            Self::TooFewFiducials { have } => {
                write!(formatter, "spatial covariance needs at least 3 fiducials, got {have}")
            }
            Self::TooManyPoints { have, max } => {
                write!(formatter, "spatial point count {have} exceeds {max}")
            }
            Self::InvalidScalar { field, requirement } => {
                write!(formatter, "{field} must be {requirement}")
            }
            Self::NonPositiveDefiniteCovariance { index } => {
                write!(formatter, "covariance {index} is not strictly positive definite")
            }
            Self::InvalidCalibrationIdentity { reason } => {
                write!(formatter, "calibration identity is invalid: {reason}")
            }
            Self::RobustCorrelationUnsupported => formatter.write_str(
                "Huber reweighting with cross-fiducial equicorrelation has no implemented covariance claim",
            ),
            Self::UnknownDependence => formatter.write_str(
                "cross-fiducial dependence is unknown; calibrated covariance is unavailable",
            ),
            Self::SingularInformation => {
                formatter.write_str("rigid transform information is singular or unresolved")
            }
            Self::AmbiguousGlobalMinimum => formatter.write_str(
                "calibrated registration has multiple globally minimizing rotations",
            ),
            Self::ArithmeticOverflow { field } => {
                write!(formatter, "spatial uncertainty aggregate {field} overflowed")
            }
            Self::AllocationFailed => formatter.write_str("spatial output allocation failed"),
            Self::Cancelled { phase } => {
                write!(formatter, "spatial uncertainty cancelled during {phase}")
            }
        }
    }
}

impl std::error::Error for SpatialUncertaintyError {}

fn finite(field: &'static str, value: f64) -> Result<f64, SpatialUncertaintyError> {
    if value.is_finite() {
        Ok(if value == 0.0 { 0.0 } else { value })
    } else {
        Err(SpatialUncertaintyError::InvalidScalar {
            field,
            requirement: "finite",
        })
    }
}

fn finite_non_negative(field: &'static str, value: f64) -> Result<f64, SpatialUncertaintyError> {
    let value = finite(field, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(SpatialUncertaintyError::InvalidScalar {
            field,
            requirement: "finite and non-negative",
        })
    }
}

fn checkpoint(
    cx: &fs_exec::Cx<'_>,
    ordinal: usize,
    phase: &'static str,
) -> Result<(), SpatialUncertaintyError> {
    if ordinal.is_multiple_of(POLL_STRIDE) {
        cx.checkpoint()
            .map_err(|_| SpatialUncertaintyError::Cancelled { phase })?;
    }
    Ok(())
}

/// A symmetric, strictly positive-definite 2x2 covariance matrix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Covariance2 {
    xx: f64,
    xy: f64,
    yy: f64,
}

impl Covariance2 {
    /// Construct a finite covariance. Strict positive definiteness is checked
    /// after scale normalization so tiny and huge valid matrices are treated
    /// consistently.
    ///
    /// # Errors
    /// Non-finite entries, non-positive marginal variances, overflowed trace,
    /// or a non-positive normalized determinant.
    pub fn new(xx: f64, xy: f64, yy: f64) -> Result<Self, SpatialUncertaintyError> {
        let xx = finite("covariance.xx", xx)?;
        let xy = finite("covariance.xy", xy)?;
        let yy = finite("covariance.yy", yy)?;
        let scale = xx.abs().max(xy.abs()).max(yy.abs());
        if xx <= 0.0 || yy <= 0.0 || scale == 0.0 {
            return Err(SpatialUncertaintyError::NonPositiveDefiniteCovariance { index: 0 });
        }
        let determinant = (xx / scale) * (yy / scale) - (xy / scale) * (xy / scale);
        if !determinant.is_finite() || determinant <= 0.0 {
            return Err(SpatialUncertaintyError::NonPositiveDefiniteCovariance { index: 0 });
        }
        if !(xx + yy).is_finite() {
            return Err(SpatialUncertaintyError::ArithmeticOverflow {
                field: "covariance trace",
            });
        }
        Ok(Self { xx, xy, yy })
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

    /// y variance.
    #[must_use]
    pub const fn yy(self) -> f64 {
        self.yy
    }

    /// Trace, used by distribution-free radial coverage.
    #[must_use]
    pub fn trace(self) -> f64 {
        self.xx + self.yy
    }

    /// Symmetric principal square root and its inverse. The symmetric factor
    /// makes the standardized equicorrelation model equivariant under rigid
    /// coordinate-frame rotations; a triangular Cholesky factor would not.
    fn principal_sqrt_and_inverse(
        self,
        index: usize,
    ) -> Result<([[f64; 2]; 2], [[f64; 2]; 2]), SpatialUncertaintyError> {
        let scale = self.xx.max(self.yy).max(self.xy.abs());
        let a = self.xx / scale;
        let b = self.xy / scale;
        let d = self.yy / scale;
        let determinant = a * d - b * b;
        if determinant <= 0.0 || !determinant.is_finite() {
            return Err(SpatialUncertaintyError::NonPositiveDefiniteCovariance { index });
        }
        let root_det = determinant.sqrt();
        let denominator = (a + d + 2.0 * root_det).sqrt();
        let root_scale = scale.sqrt();
        let factor = [
            [
                root_scale * (a + root_det) / denominator,
                root_scale * b / denominator,
            ],
            [
                root_scale * b / denominator,
                root_scale * (d + root_det) / denominator,
            ],
        ];
        let factor_det = factor[0][0] * factor[1][1] - factor[0][1] * factor[1][0];
        if !factor.iter().flatten().all(|value| value.is_finite())
            || !factor_det.is_finite()
            || factor_det <= 0.0
        {
            return Err(SpatialUncertaintyError::NonPositiveDefiniteCovariance { index });
        }
        let inverse = [
            [factor[1][1] / factor_det, -factor[0][1] / factor_det],
            [-factor[1][0] / factor_det, factor[0][0] / factor_det],
        ];
        if !inverse.iter().flatten().all(|value| value.is_finite()) {
            return Err(SpatialUncertaintyError::NonPositiveDefiniteCovariance { index });
        }
        Ok((factor, inverse))
    }
}

/// Declared correlation between distinct fiducial measurement pairs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CrossFiducialModel {
    /// Distinct fiducial errors are independent; each pair may still carry
    /// x/y covariance through [`Covariance2`].
    Independent,
    /// After each pair is whitened by its symmetric principal square root,
    /// every distinct pair has correlation `rho * I2`.
    EquicorrelatedStandardized {
        /// Common standardized pair correlation.
        rho: f64,
    },
    /// Dependence exists but no cross-covariance model was supplied. Estimation
    /// refuses rather than silently using independence.
    Unknown,
}

/// Deterministic robust weighting policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HuberPolicy {
    /// Model-based generalized least squares without data-dependent weights.
    Disabled,
    /// Fixed-pass Huber reweighting of each whitened residual norm.
    Enabled {
        /// Positive standardized Huber threshold.
        threshold: f64,
        /// Exact deterministic number of reweighting passes.
        iterations: u8,
    },
}

impl HuberPolicy {
    /// Construct a bounded Huber policy.
    ///
    /// # Errors
    /// A non-positive/non-finite threshold or iteration count outside
    /// `1..=MAX_HUBER_ITERATIONS`.
    pub fn new(threshold: f64, iterations: u8) -> Result<Self, SpatialUncertaintyError> {
        if !threshold.is_finite() || threshold <= 0.0 {
            return Err(SpatialUncertaintyError::InvalidScalar {
                field: "huber.threshold",
                requirement: "finite and positive",
            });
        }
        if iterations == 0 || iterations > MAX_HUBER_ITERATIONS {
            return Err(SpatialUncertaintyError::InvalidScalar {
                field: "huber.iterations",
                requirement: "between 1 and MAX_HUBER_ITERATIONS",
            });
        }
        Ok(Self::Enabled {
            threshold,
            iterations,
        })
    }
}

/// Bound on the total non-stochastic registered-inspection error shared by the
/// complete queried family.
///
/// This is deliberately not a raw fiducial or scanner calibration error. The
/// caller must propagate every systematic registration and inspection effect
/// over the declared query domain and supply one radial upper bound on their
/// combined vector error. Applying that already-propagated bound once avoids
/// both rotational under-counting and double counting.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BiasBound {
    /// Finite radial upper bound in the same length unit as the points.
    Bounded(f64),
    /// No finite bias bound is available; decisions remain indeterminate.
    Unbounded,
}

/// Complete calibration/noise model for one fiducial family.
#[derive(Debug, Clone, PartialEq)]
pub struct MetrologyModel {
    fiducial_covariances: Vec<Covariance2>,
    cross_fiducial: CrossFiducialModel,
    huber: HuberPolicy,
    registered_inspection_bias: BiasBound,
    calibration_identity: String,
}

impl MetrologyModel {
    /// Construct a model whose per-fiducial covariance order matches the
    /// correspondence order supplied to [`estimate_calibrated_registration`].
    ///
    /// # Errors
    /// Invalid sizes, correlation domain, bias bound, or calibration identity.
    pub fn new(
        fiducial_covariances: Vec<Covariance2>,
        cross_fiducial: CrossFiducialModel,
        huber: HuberPolicy,
        registered_inspection_bias: BiasBound,
        calibration_identity: impl Into<String>,
    ) -> Result<Self, SpatialUncertaintyError> {
        let count = fiducial_covariances.len();
        if count < 3 {
            return Err(SpatialUncertaintyError::TooFewFiducials { have: count });
        }
        if count > MAX_AS_BUILT_POINTS {
            return Err(SpatialUncertaintyError::TooManyPoints {
                have: count,
                max: MAX_AS_BUILT_POINTS,
            });
        }
        match huber {
            HuberPolicy::Disabled => {}
            HuberPolicy::Enabled {
                threshold,
                iterations,
            } => {
                if !threshold.is_finite() || threshold <= 0.0 {
                    return Err(SpatialUncertaintyError::InvalidScalar {
                        field: "huber.threshold",
                        requirement: "finite and positive",
                    });
                }
                if iterations == 0 || iterations > MAX_HUBER_ITERATIONS {
                    return Err(SpatialUncertaintyError::InvalidScalar {
                        field: "huber.iterations",
                        requirement: "between 1 and MAX_HUBER_ITERATIONS",
                    });
                }
            }
        }
        match cross_fiducial {
            CrossFiducialModel::Independent => {}
            CrossFiducialModel::Unknown => {}
            CrossFiducialModel::EquicorrelatedStandardized { rho } => {
                let lower = -1.0 / ((count - 1) as f64);
                if !rho.is_finite() || rho <= lower || rho >= 1.0 {
                    return Err(SpatialUncertaintyError::InvalidScalar {
                        field: "cross_fiducial.rho",
                        requirement: "strictly inside (-1/(n-1), 1)",
                    });
                }
                if !matches!(huber, HuberPolicy::Disabled) {
                    return Err(SpatialUncertaintyError::RobustCorrelationUnsupported);
                }
            }
        }
        let registered_inspection_bias = match registered_inspection_bias {
            BiasBound::Bounded(value) => {
                BiasBound::Bounded(finite_non_negative("registered_inspection_bias", value)?)
            }
            BiasBound::Unbounded => BiasBound::Unbounded,
        };
        let calibration_identity = calibration_identity.into();
        if let Some(reason) = fs_evidence::color_leaf_identity_reason(&calibration_identity) {
            return Err(SpatialUncertaintyError::InvalidCalibrationIdentity { reason });
        }
        Ok(Self {
            fiducial_covariances,
            cross_fiducial,
            huber,
            registered_inspection_bias,
            calibration_identity,
        })
    }

    /// Ordered calibrated fiducial covariances.
    #[must_use]
    pub fn fiducial_covariances(&self) -> &[Covariance2] {
        &self.fiducial_covariances
    }

    /// Declared cross-fiducial model.
    #[must_use]
    pub const fn cross_fiducial(&self) -> CrossFiducialModel {
        self.cross_fiducial
    }

    /// Robust policy.
    #[must_use]
    pub const fn huber(&self) -> HuberPolicy {
        self.huber
    }

    /// Total registered-inspection radial bias bound over the query domain.
    #[must_use]
    pub const fn registered_inspection_bias(&self) -> BiasBound {
        self.registered_inspection_bias
    }

    /// Calibration artifact identity. It is bound into the model root but is
    /// not self-authenticating.
    #[must_use]
    pub fn calibration_identity(&self) -> &str {
        &self.calibration_identity
    }
}

/// Row-major covariance matrix for rigid parameters `(tx, ty, theta)`.
pub type TransformCovariance = [[f64; 3]; 3];

type Mat3 = TransformCovariance;
type Mat23 = [[f64; 3]; 2];
type Mat4 = [[f64; 4]; 4];
type Mat24 = [[f64; 4]; 2];

const ZERO3: Mat3 = [[0.0; 3]; 3];
const ZERO4: Mat4 = [[0.0; 4]; 4];

fn add_outer3(target: &mut Mat3, rows: &Mat23, scale: f64) {
    for (left, target_row) in target.iter_mut().enumerate() {
        for (right, value) in target_row.iter_mut().enumerate() {
            *value += scale * (rows[0][left] * rows[0][right] + rows[1][left] * rows[1][right]);
        }
    }
}

fn add_outer4(target: &mut Mat4, rows: &Mat24, scale: f64) {
    for (left, target_row) in target.iter_mut().enumerate() {
        for (right, value) in target_row.iter_mut().enumerate() {
            *value += scale * (rows[0][left] * rows[0][right] + rows[1][left] * rows[1][right]);
        }
    }
}

fn add_transpose_product4(target: &mut [f64; 4], rows: &Mat24, vector: [f64; 2], scale: f64) {
    for component in 0..4 {
        target[component] +=
            scale * (rows[0][component] * vector[0] + rows[1][component] * vector[1]);
    }
}

fn mat3_vec(matrix: &Mat3, vector: [f64; 3]) -> [f64; 3] {
    let mut result = [0.0; 3];
    for (row, output) in result.iter_mut().enumerate() {
        *output =
            matrix[row][0] * vector[0] + matrix[row][1] * vector[1] + matrix[row][2] * vector[2];
    }
    result
}

fn mat3_mul(left: &Mat3, right: &Mat3) -> Mat3 {
    let mut product = ZERO3;
    for (row, product_row) in product.iter_mut().enumerate() {
        for (column, value) in product_row.iter_mut().enumerate() {
            *value = (0..3)
                .map(|inner| left[row][inner] * right[inner][column])
                .sum();
        }
    }
    product
}

fn transpose3(matrix: &Mat3) -> Mat3 {
    let mut transpose = ZERO3;
    for (row, values) in matrix.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            transpose[column][row] = *value;
        }
    }
    transpose
}

fn invert_spd2(matrix: [[f64; 2]; 2]) -> Result<[[f64; 2]; 2], SpatialUncertaintyError> {
    let a = matrix[0][0];
    let d = matrix[1][1];
    if !a.is_finite() || !d.is_finite() || a <= 0.0 || d <= 0.0 {
        return Err(SpatialUncertaintyError::SingularInformation);
    }
    let a_scale = a.sqrt();
    let d_scale = d.sqrt();
    let correlation = 0.5 * (matrix[0][1] + matrix[1][0]) / a_scale / d_scale;
    let determinant = 1.0 - correlation * correlation;
    if !correlation.is_finite() || !determinant.is_finite() || determinant <= 256.0 * f64::EPSILON {
        return Err(SpatialUncertaintyError::SingularInformation);
    }
    let inverse = [
        [
            1.0 / determinant / a,
            -correlation / determinant / a_scale / d_scale,
        ],
        [
            -correlation / determinant / a_scale / d_scale,
            1.0 / determinant / d,
        ],
    ];
    if inverse.iter().flatten().all(|value| value.is_finite()) {
        Ok(inverse)
    } else {
        Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "inverse translation information",
        })
    }
}

fn mat2_vec(matrix: [[f64; 2]; 2], vector: [f64; 2]) -> [f64; 2] {
    [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    ]
}

/// Invert a symmetric positive-definite 3x3 matrix through a diagonally
/// equilibrated Cholesky factor. Translation and rotation have different
/// units, so one global scale would make admission depend on the chosen length
/// unit. The relative pivot guard is part of the geometry/no-claim boundary;
/// unresolved leverage refuses rather than being clamped.
fn invert_spd3(matrix: &Mat3) -> Result<Mat3, SpatialUncertaintyError> {
    if !matrix.iter().flatten().all(|value| value.is_finite()) {
        return Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "information matrix",
        });
    }
    let mut scales = [0.0; 3];
    for index in 0..3 {
        let diagonal = matrix[index][index];
        if !diagonal.is_finite() || diagonal <= 0.0 {
            return Err(SpatialUncertaintyError::SingularInformation);
        }
        scales[index] = diagonal.sqrt();
    }
    let mut normalized = ZERO3;
    for row in 0..3 {
        for column in 0..3 {
            // Dividing successively avoids forming a possibly overflowing
            // product of two dimensional scales.
            normalized[row][column] = matrix[row][column] / scales[row] / scales[column];
            if !normalized[row][column].is_finite() {
                return Err(SpatialUncertaintyError::SingularInformation);
            }
        }
    }
    let mut lower = ZERO3;
    let pivot_floor = 256.0 * f64::EPSILON;
    for row in 0..3 {
        for column in 0..=row {
            let mut value = normalized[row][column];
            for prior in 0..column {
                value -= lower[row][prior] * lower[column][prior];
            }
            if row == column {
                if !value.is_finite() || value <= pivot_floor {
                    return Err(SpatialUncertaintyError::SingularInformation);
                }
                lower[row][column] = value.sqrt();
            } else {
                lower[row][column] = value / lower[column][column];
            }
        }
    }
    let mut inverse_normalized = ZERO3;
    for basis in 0..3 {
        let mut forward = [0.0; 3];
        for row in 0..3 {
            let mut value = if row == basis { 1.0 } else { 0.0 };
            for prior in 0..row {
                value -= lower[row][prior] * forward[prior];
            }
            forward[row] = value / lower[row][row];
        }
        let mut backward = [0.0; 3];
        for row in (0..3).rev() {
            let mut value = forward[row];
            for later in row + 1..3 {
                value -= lower[later][row] * backward[later];
            }
            backward[row] = value / lower[row][row];
        }
        for row in 0..3 {
            inverse_normalized[row][basis] = backward[row];
        }
    }
    let mut inverse = ZERO3;
    for row in 0..3 {
        for column in 0..3 {
            inverse[row][column] = inverse_normalized[row][column] / scales[row] / scales[column];
        }
    }
    if inverse.iter().flatten().all(|value| value.is_finite()) {
        Ok(inverse)
    } else {
        Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "inverse information matrix",
        })
    }
}

fn symmetrize_and_validate_covariance(
    mut covariance: Mat3,
) -> Result<Mat3, SpatialUncertaintyError> {
    if !covariance.iter().flatten().all(|value| value.is_finite()) {
        return Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "registration covariance",
        });
    }
    for row in 0..3 {
        if covariance[row][row] <= 0.0 {
            return Err(SpatialUncertaintyError::SingularInformation);
        }
        for column in row + 1..3 {
            let symmetric = f64::midpoint(covariance[row][column], covariance[column][row]);
            covariance[row][column] = symmetric;
            covariance[column][row] = symmetric;
        }
    }
    let scales = [
        covariance[0][0].sqrt(),
        covariance[1][1].sqrt(),
        covariance[2][2].sqrt(),
    ];
    let mut lower = ZERO3;
    let pivot_floor = 256.0 * f64::EPSILON;
    for row in 0..3 {
        for column in 0..=row {
            let mut value = covariance[row][column] / scales[row] / scales[column];
            for prior in 0..column {
                value -= lower[row][prior] * lower[column][prior];
            }
            if row == column {
                if !value.is_finite() || value <= pivot_floor {
                    return Err(SpatialUncertaintyError::SingularInformation);
                }
                lower[row][column] = value.sqrt();
            } else {
                lower[row][column] = value / lower[column][column];
            }
        }
    }
    Ok(covariance)
}

fn whiten(inverse_factor: [[f64; 2]; 2], rows: Mat23, vector: [f64; 2]) -> (Mat23, [f64; 2]) {
    let mut whitened_rows = [[0.0; 3]; 2];
    for row in 0..2 {
        for column in 0..3 {
            whitened_rows[row][column] =
                inverse_factor[row][0] * rows[0][column] + inverse_factor[row][1] * rows[1][column];
        }
    }
    let whitened_vector = [
        inverse_factor[0][0] * vector[0] + inverse_factor[0][1] * vector[1],
        inverse_factor[1][0] * vector[0] + inverse_factor[1][1] * vector[1],
    ];
    (whitened_rows, whitened_vector)
}

fn whiten4(inverse_factor: [[f64; 2]; 2], rows: Mat24, vector: [f64; 2]) -> (Mat24, [f64; 2]) {
    let mut whitened_rows = [[0.0; 4]; 2];
    for row in 0..2 {
        for column in 0..4 {
            whitened_rows[row][column] =
                inverse_factor[row][0] * rows[0][column] + inverse_factor[row][1] * rows[1][column];
        }
    }
    let whitened_vector = [
        inverse_factor[0][0] * vector[0] + inverse_factor[0][1] * vector[1],
        inverse_factor[1][0] * vector[0] + inverse_factor[1][1] * vector[1],
    ];
    (whitened_rows, whitened_vector)
}

#[derive(Debug)]
struct LinearFitSystem {
    information: Mat4,
    target: [f64; 4],
    design_pivot: Point2,
}

#[allow(clippy::too_many_lines)]
fn linear_fit_system(
    fiducials: &[Fiducial],
    model: &MetrologyModel,
    weights: &[f64],
    cx: &fs_exec::Cx<'_>,
) -> Result<LinearFitSystem, SpatialUncertaintyError> {
    let design_pivot = fiducials
        .first()
        .ok_or(SpatialUncertaintyError::TooFewFiducials { have: 0 })?
        .design();
    let mut information = ZERO4;
    let mut target = [0.0; 4];
    let mut row_sum = [[0.0; 4]; 2];
    let mut measured_sum = [0.0; 2];
    for (index, ((fiducial, covariance), weight)) in fiducials
        .iter()
        .zip(&model.fiducial_covariances)
        .zip(weights)
        .enumerate()
    {
        checkpoint(cx, index, "global registration fit")?;
        let (_, inverse_factor) = covariance.principal_sqrt_and_inverse(index)?;
        let design = fiducial.design();
        let measured = fiducial.measured();
        let design_x = design.x() - design_pivot.x();
        let design_y = design.y() - design_pivot.y();
        if !design_x.is_finite() || !design_y.is_finite() {
            return Err(SpatialUncertaintyError::ArithmeticOverflow {
                field: "pivoted design coordinate",
            });
        }
        let rows = [
            [1.0, 0.0, design_x, -design_y],
            [0.0, 1.0, design_y, design_x],
        ];
        let (rows, measured) = whiten4(inverse_factor, rows, [measured.x(), measured.y()]);
        if !rows.iter().flatten().all(|value| value.is_finite())
            || !measured.iter().all(|value| value.is_finite())
            || !weight.is_finite()
            || *weight <= 0.0
        {
            return Err(SpatialUncertaintyError::ArithmeticOverflow {
                field: "global registration system",
            });
        }
        match model.cross_fiducial {
            CrossFiducialModel::Independent => {
                add_outer4(&mut information, &rows, *weight);
                add_transpose_product4(&mut target, &rows, measured, *weight);
            }
            CrossFiducialModel::EquicorrelatedStandardized { .. } => {
                add_outer4(&mut information, &rows, 1.0);
                add_transpose_product4(&mut target, &rows, measured, 1.0);
                for row in 0..2 {
                    measured_sum[row] += measured[row];
                    for column in 0..4 {
                        row_sum[row][column] += rows[row][column];
                    }
                }
            }
            CrossFiducialModel::Unknown => {
                return Err(SpatialUncertaintyError::UnknownDependence);
            }
        }
    }
    cx.checkpoint()
        .map_err(|_| SpatialUncertaintyError::Cancelled {
            phase: "global registration fit",
        })?;
    if let CrossFiducialModel::EquicorrelatedStandardized { rho } = model.cross_fiducial {
        let count = fiducials.len() as f64;
        let denominator = 1.0 + (count - 1.0) * rho;
        let diagonal = 1.0 / (1.0 - rho);
        let rank_one = -rho / ((1.0 - rho) * denominator);
        for row in 0..4 {
            for column in 0..4 {
                let aggregate =
                    row_sum[0][row] * row_sum[0][column] + row_sum[1][row] * row_sum[1][column];
                information[row][column] =
                    diagonal * information[row][column] + rank_one * aggregate;
            }
            let aggregate = row_sum[0][row] * measured_sum[0] + row_sum[1][row] * measured_sum[1];
            target[row] = diagonal * target[row] + rank_one * aggregate;
        }
    }
    if !information
        .iter()
        .flatten()
        .chain(target.iter())
        .all(|value| value.is_finite())
    {
        return Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "global registration normal equations",
        });
    }
    Ok(LinearFitSystem {
        information,
        target,
        design_pivot,
    })
}

/// Solve `min u' H u - 2 g' u` over the complete unit circle. This is the
/// two-dimensional trust-region subproblem obtained after eliminating rigid
/// translation. The secular equation is monotone above the smallest
/// eigenvalue; its hard case is handled explicitly.
#[derive(Debug, Clone, Copy, PartialEq)]
struct UnitCircleMinimum {
    unit: [f64; 2],
    multiplier: f64,
}

#[allow(clippy::too_many_lines)]
fn minimize_on_unit_circle(
    hessian: [[f64; 2]; 2],
    target: [f64; 2],
) -> Result<UnitCircleMinimum, SpatialUncertaintyError> {
    let a = hessian[0][0];
    let b = f64::midpoint(hessian[0][1], hessian[1][0]);
    let d = hessian[1][1];
    if ![a, b, d, target[0], target[1]]
        .into_iter()
        .all(f64::is_finite)
    {
        return Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "profiled rotation system",
        });
    }
    let midpoint = f64::midpoint(a, d);
    let radius = (0.5 * (a - d)).hypot(b);
    if !midpoint.is_finite() || !radius.is_finite() {
        return Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "profiled rotation eigensystem",
        });
    }
    let smallest = midpoint - radius;
    let largest = midpoint + radius;
    if smallest < 0.0 {
        return Err(SpatialUncertaintyError::SingularInformation);
    }
    let target_norm = target[0].hypot(target[1]);
    if target_norm == 0.0 {
        // A purely quadratic anisotropic profile has the antipodal pair of
        // minimum-axis rotations; an isotropic profile has every rotation.
        return Err(SpatialUncertaintyError::AmbiguousGlobalMinimum);
    }
    if radius == 0.0 {
        return Ok(UnitCircleMinimum {
            unit: [target[0] / target_norm, target[1] / target_norm],
            multiplier: target_norm - smallest,
        });
    }

    // Classify eigen-alignment before constructing trigonometric eigenvectors.
    // An exact hard case can otherwise acquire a tiny spurious minimum-axis
    // projection from `sin_cos` and enter the regular secular solve. Outward
    // arithmetic makes any unresolved alignment fail closed.
    let interval_a = Interval::point(a);
    let interval_b = Interval::point(b);
    let interval_d = Interval::point(d);
    let interval_h0 = Interval::point(target[0]);
    let interval_h1 = Interval::point(target[1]);
    let applied_h0 = interval_a * interval_h0 + interval_b * interval_h1;
    let applied_h1 = interval_b * interval_h0 + interval_d * interval_h1;
    let eigen_cross = interval_h0 * applied_h1 - interval_h1 * applied_h0;
    if eigen_cross.contains_zero() {
        let interval_midpoint = (interval_a + interval_d) * Interval::point(0.5);
        let orientation = interval_h0
            * ((interval_a - interval_midpoint) * interval_h0 + interval_b * interval_h1)
            + interval_h1
                * (interval_b * interval_h0 + (interval_d - interval_midpoint) * interval_h1);
        let half_difference = (interval_a - interval_d) * Interval::point(0.5);
        let interval_gap = Interval::point(2.0)
            * (half_difference * half_difference + interval_b * interval_b).sqrt();
        let interval_target_norm = (interval_h0 * interval_h0 + interval_h1 * interval_h1).sqrt();
        if orientation.lo() > 0.0 {
            if interval_target_norm.hi() < interval_gap.lo() {
                return Err(SpatialUncertaintyError::AmbiguousGlobalMinimum);
            }
            if interval_target_norm.lo() <= interval_gap.hi() {
                // Equality is a unique but quartically flat minimizer; an
                // overlap also cannot certify positive tangent curvature.
                return Err(SpatialUncertaintyError::SingularInformation);
            }
            // A maximum-axis target beyond the eigengap is regular and safely
            // separated from the hard boundary; continue to bisection.
        } else if orientation.hi() >= 0.0 {
            return Err(SpatialUncertaintyError::SingularInformation);
        }
        // Minimum-axis alignment is also regular; bisection remains separated
        // from the singular multiplier and the final KKT gate certifies it.
    }

    // `maximum_axis` is the eigenvector of the largest eigenvalue. The
    // half-angle construction fixes its sign deterministically.
    let angle = 0.5 * (2.0 * b).atan2(a - d);
    let (axis_sine, axis_cosine) = angle.sin_cos();
    let maximum_axis = [axis_cosine, axis_sine];
    let minimum_axis = [-axis_sine, axis_cosine];
    let projected = [
        minimum_axis[0] * target[0] + minimum_axis[1] * target[1],
        maximum_axis[0] * target[0] + maximum_axis[1] * target[1],
    ];
    let mut lower_multiplier = -smallest;
    let mut upper_multiplier = target_norm - smallest;
    if !lower_multiplier.is_finite()
        || !upper_multiplier.is_finite()
        || upper_multiplier <= lower_multiplier
    {
        return Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "profiled rotation multiplier",
        });
    }
    // Fixed bisection is deterministic. The secular residual decreases
    // strictly from positive infinity to at most zero on this bracket.
    for _ in 0..128 {
        let multiplier = f64::midpoint(lower_multiplier, upper_multiplier);
        let first = projected[0] / (smallest + multiplier);
        let second = projected[1] / (largest + multiplier);
        let residual = first * first + second * second - 1.0;
        if !residual.is_finite() || residual > 0.0 {
            lower_multiplier = multiplier;
        } else {
            upper_multiplier = multiplier;
        }
    }
    let first = projected[0] / (smallest + upper_multiplier);
    let second = projected[1] / (largest + upper_multiplier);
    let coordinate_norm = first.hypot(second);
    if !coordinate_norm.is_finite() || coordinate_norm == 0.0 {
        return Err(SpatialUncertaintyError::SingularInformation);
    }
    let coordinates = [first / coordinate_norm, second / coordinate_norm];
    let multiplier = upper_multiplier;
    let unit = [
        minimum_axis[0] * coordinates[0] + maximum_axis[0] * coordinates[1],
        minimum_axis[1] * coordinates[0] + maximum_axis[1] * coordinates[1],
    ];
    let norm = unit[0].hypot(unit[1]);
    if !norm.is_finite() || norm == 0.0 {
        return Err(SpatialUncertaintyError::SingularInformation);
    }
    let multiplier = if multiplier == 0.0 { 0.0 } else { multiplier };
    if !multiplier.is_finite() || multiplier < -smallest {
        return Err(SpatialUncertaintyError::SingularInformation);
    }
    let unit = [unit[0] / norm, unit[1] / norm];
    let kkt_residual = [
        (a + multiplier) * unit[0] + b * unit[1] - target[0],
        b * unit[0] + (d + multiplier) * unit[1] - target[1],
    ];
    let kkt_scale = target_norm
        .max(a.abs())
        .max(b.abs())
        .max(d.abs())
        .max(multiplier.abs())
        .max(f64::MIN_POSITIVE);
    let kkt_tolerance = 4096.0 * f64::EPSILON * kkt_scale;
    let kkt_norm = kkt_residual[0].hypot(kkt_residual[1]);
    if !kkt_tolerance.is_finite() || !kkt_norm.is_finite() || kkt_norm > kkt_tolerance {
        return Err(SpatialUncertaintyError::SingularInformation);
    }
    Ok(UnitCircleMinimum { unit, multiplier })
}

#[derive(Debug, Clone, Copy)]
struct FixedWeightSolution {
    parameters: [f64; 3],
    trust_multiplier: f64,
}

fn solve_global_fixed_weights(
    fiducials: &[Fiducial],
    model: &MetrologyModel,
    weights: &[f64],
    cx: &fs_exec::Cx<'_>,
) -> Result<FixedWeightSolution, SpatialUncertaintyError> {
    let system = linear_fit_system(fiducials, model, weights, cx)?;
    let translation_information = [
        [system.information[0][0], system.information[0][1]],
        [system.information[1][0], system.information[1][1]],
    ];
    let translation_inverse = invert_spd2(translation_information)?;
    let coupling = [
        [system.information[0][2], system.information[0][3]],
        [system.information[1][2], system.information[1][3]],
    ];
    let rotation_information = [
        [system.information[2][2], system.information[2][3]],
        [system.information[3][2], system.information[3][3]],
    ];
    let translation_target = [system.target[0], system.target[1]];
    let inverse_target = mat2_vec(translation_inverse, translation_target);
    let mut profiled_information = rotation_information;
    let mut profiled_target = [system.target[2], system.target[3]];
    for left in 0..2 {
        profiled_target[left] -=
            coupling[0][left] * inverse_target[0] + coupling[1][left] * inverse_target[1];
        for right in 0..2 {
            let inverse_coupling = [
                translation_inverse[0][0] * coupling[0][right]
                    + translation_inverse[0][1] * coupling[1][right],
                translation_inverse[1][0] * coupling[0][right]
                    + translation_inverse[1][1] * coupling[1][right],
            ];
            profiled_information[left][right] -=
                coupling[0][left] * inverse_coupling[0] + coupling[1][left] * inverse_coupling[1];
        }
    }
    let minimum = minimize_on_unit_circle(profiled_information, profiled_target)?;
    let rotation = minimum.unit;
    let translation_rhs = [
        translation_target[0] - coupling[0][0] * rotation[0] - coupling[0][1] * rotation[1],
        translation_target[1] - coupling[1][0] * rotation[0] - coupling[1][1] * rotation[1],
    ];
    let translation = mat2_vec(translation_inverse, translation_rhs);
    let pivot = system.design_pivot;
    let origin_translation = [
        translation[0] - (rotation[0] * pivot.x() - rotation[1] * pivot.y()),
        translation[1] - (rotation[1] * pivot.x() + rotation[0] * pivot.y()),
    ];
    let parameters = [
        origin_translation[0],
        origin_translation[1],
        rotation[1].atan2(rotation[0]),
    ];
    if parameters.iter().all(|value| value.is_finite()) {
        Ok(FixedWeightSolution {
            parameters,
            trust_multiplier: minimum.multiplier,
        })
    } else {
        Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "global registration solution",
        })
    }
}

fn predicted_and_jacobian(
    fiducial: Fiducial,
    parameters: [f64; 3],
) -> Result<(Point2, Mat23), SpatialUncertaintyError> {
    let (sine, cosine) = parameters[2].sin_cos();
    let design = fiducial.design();
    let rotated_x = cosine * design.x() - sine * design.y();
    let rotated_y = sine * design.x() + cosine * design.y();
    let predicted =
        Point2::new(rotated_x + parameters[0], rotated_y + parameters[1]).map_err(|_| {
            SpatialUncertaintyError::ArithmeticOverflow {
                field: "predicted fiducial",
            }
        })?;
    let jacobian = [[1.0, 0.0, -rotated_y], [0.0, 1.0, rotated_x]];
    Ok((predicted, jacobian))
}

#[derive(Debug)]
struct Equations {
    information: Mat3,
    standardized_norms: Vec<f64>,
    whitened_rows: Vec<Mat23>,
}

fn equations(
    fiducials: &[Fiducial],
    model: &MetrologyModel,
    parameters: [f64; 3],
    weights: &[f64],
    cx: &fs_exec::Cx<'_>,
) -> Result<Equations, SpatialUncertaintyError> {
    let mut information_sum = ZERO3;
    let mut row_sum = [[0.0; 3]; 2];
    let mut standardized_norms = Vec::new();
    let mut whitened_rows = Vec::new();
    standardized_norms
        .try_reserve_exact(fiducials.len())
        .map_err(|_| SpatialUncertaintyError::AllocationFailed)?;
    whitened_rows
        .try_reserve_exact(fiducials.len())
        .map_err(|_| SpatialUncertaintyError::AllocationFailed)?;

    for (index, ((fiducial, covariance), weight)) in fiducials
        .iter()
        .zip(&model.fiducial_covariances)
        .zip(weights)
        .enumerate()
    {
        checkpoint(cx, index, "registration-covariance equations")?;
        let (_, inverse_factor) = covariance.principal_sqrt_and_inverse(index)?;
        let (predicted, jacobian) = predicted_and_jacobian(*fiducial, parameters)?;
        let measured = fiducial.measured();
        let residual = [measured.x() - predicted.x(), measured.y() - predicted.y()];
        let (rows, standardized) = whiten(inverse_factor, jacobian, residual);
        if !rows.iter().flatten().all(|value| value.is_finite())
            || !standardized.iter().all(|value| value.is_finite())
        {
            return Err(SpatialUncertaintyError::ArithmeticOverflow {
                field: "whitened registration system",
            });
        }
        standardized_norms.push(standardized[0].hypot(standardized[1]));
        whitened_rows.push(rows);
        match model.cross_fiducial {
            CrossFiducialModel::Independent => {
                add_outer3(&mut information_sum, &rows, *weight);
            }
            CrossFiducialModel::EquicorrelatedStandardized { .. } => {
                add_outer3(&mut information_sum, &rows, 1.0);
                for row in 0..2 {
                    for column in 0..3 {
                        row_sum[row][column] += rows[row][column];
                    }
                }
            }
            CrossFiducialModel::Unknown => {
                return Err(SpatialUncertaintyError::UnknownDependence);
            }
        }
    }
    cx.checkpoint()
        .map_err(|_| SpatialUncertaintyError::Cancelled {
            phase: "registration-covariance equations",
        })?;

    if let CrossFiducialModel::EquicorrelatedStandardized { rho } = model.cross_fiducial {
        let count = fiducials.len() as f64;
        let denominator = 1.0 + (count - 1.0) * rho;
        let diagonal = 1.0 / (1.0 - rho);
        let rank_one = -rho / ((1.0 - rho) * denominator);
        for row in 0..3 {
            for column in 0..3 {
                let aggregate =
                    row_sum[0][row] * row_sum[0][column] + row_sum[1][row] * row_sum[1][column];
                information_sum[row][column] =
                    diagonal * information_sum[row][column] + rank_one * aggregate;
            }
        }
    }
    if !information_sum
        .iter()
        .flatten()
        .all(|value| value.is_finite())
    {
        return Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "registration normal equations",
        });
    }
    Ok(Equations {
        information: information_sum,
        standardized_norms,
        whitened_rows,
    })
}

fn robust_sandwich(
    bread: &Mat3,
    rows: &[Mat23],
    weights: &[f64],
    cx: &fs_exec::Cx<'_>,
) -> Result<Mat3, SpatialUncertaintyError> {
    let mut meat = ZERO3;
    for (index, (row, weight)) in rows.iter().zip(weights).enumerate() {
        checkpoint(cx, index, "robust covariance sandwich")?;
        add_outer3(&mut meat, row, weight * weight);
    }
    cx.checkpoint()
        .map_err(|_| SpatialUncertaintyError::Cancelled {
            phase: "robust covariance sandwich",
        })?;
    Ok(mat3_mul(&mat3_mul(bread, &meat), &transpose3(bread)))
}

fn full_model_leverage(
    bread: &Mat3,
    rows: &[Mat23],
    weights: &[f64],
    cross_fiducial: CrossFiducialModel,
    cx: &fs_exec::Cx<'_>,
) -> Result<Vec<f64>, SpatialUncertaintyError> {
    let mut leverage = Vec::new();
    leverage
        .try_reserve_exact(rows.len())
        .map_err(|_| SpatialUncertaintyError::AllocationFailed)?;
    let (precision_diagonal, precision_rank_one, row_sum) = match cross_fiducial {
        CrossFiducialModel::Independent => (1.0, 0.0, [[0.0; 3]; 2]),
        CrossFiducialModel::EquicorrelatedStandardized { rho } => {
            let mut sum = [[0.0; 3]; 2];
            for (index, pair) in rows.iter().enumerate() {
                checkpoint(cx, index, "correlated leverage preparation")?;
                for row in 0..2 {
                    for column in 0..3 {
                        sum[row][column] += pair[row][column];
                    }
                }
            }
            let count = rows.len() as f64;
            let denominator = 1.0 + (count - 1.0) * rho;
            (1.0 / (1.0 - rho), -rho / ((1.0 - rho) * denominator), sum)
        }
        CrossFiducialModel::Unknown => return Err(SpatialUncertaintyError::UnknownDependence),
    };
    for (index, (pair, weight)) in rows.iter().zip(weights).enumerate() {
        checkpoint(cx, index, "fiducial leverage")?;
        let mut precision_rows = *pair;
        match cross_fiducial {
            CrossFiducialModel::Independent => {
                for value in precision_rows.iter_mut().flatten() {
                    *value *= *weight;
                }
            }
            CrossFiducialModel::EquicorrelatedStandardized { .. } => {
                for row in 0..2 {
                    for column in 0..3 {
                        precision_rows[row][column] = precision_diagonal * pair[row][column]
                            + precision_rank_one * row_sum[row][column];
                    }
                }
            }
            CrossFiducialModel::Unknown => {
                return Err(SpatialUncertaintyError::UnknownDependence);
            }
        }
        let mut trace = 0.0;
        for row in 0..2 {
            let projected = mat3_vec(bread, precision_rows[row]);
            trace += pair[row][0] * projected[0]
                + pair[row][1] * projected[1]
                + pair[row][2] * projected[2];
        }
        if !trace.is_finite()
            || (trace < 0.0 && matches!(cross_fiducial, CrossFiducialModel::Independent))
        {
            return Err(SpatialUncertaintyError::ArithmeticOverflow {
                field: "fiducial leverage",
            });
        }
        leverage.push(trace);
    }
    cx.checkpoint()
        .map_err(|_| SpatialUncertaintyError::Cancelled {
            phase: "fiducial leverage",
        })?;
    Ok(leverage)
}

fn require_rank_two_design(
    fiducials: &[Fiducial],
    cx: &fs_exec::Cx<'_>,
) -> Result<(), SpatialUncertaintyError> {
    let mut scale = 0.0f64;
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "registration rank scale")?;
        scale = scale
            .max(fiducial.design().x().abs())
            .max(fiducial.design().y().abs());
    }
    if scale == 0.0 || !scale.is_finite() {
        return Err(SpatialUncertaintyError::SingularInformation);
    }
    let count = fiducials.len() as f64;
    let (mut mean_x, mut mean_y) = (0.0, 0.0);
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "registration rank centroid")?;
        mean_x += fiducial.design().x() / scale / count;
        mean_y += fiducial.design().y() / scale / count;
    }
    let (mut xx, mut xy, mut yy) = (0.0, 0.0, 0.0);
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "registration rank scatter")?;
        let x = fiducial.design().x() / scale - mean_x;
        let y = fiducial.design().y() / scale - mean_y;
        xx += x * x;
        xy += x * y;
        yy += y * y;
    }
    cx.checkpoint()
        .map_err(|_| SpatialUncertaintyError::Cancelled {
            phase: "registration rank scatter",
        })?;
    let trace = xx + yy;
    let determinant = xx * yy - xy * xy;
    if !trace.is_finite()
        || !determinant.is_finite()
        || trace <= 0.0
        || determinant <= 256.0 * f64::EPSILON * trace * trace
    {
        return Err(SpatialUncertaintyError::SingularInformation);
    }
    Ok(())
}

#[allow(clippy::too_many_lines)]
fn model_identity(
    fiducials: &[Fiducial],
    model: &MetrologyModel,
    registration: &Registration,
    covariance: &Mat3,
    diagnostics: &[OutlierDiagnostic],
    leverage: &[f64],
    cx: &fs_exec::Cx<'_>,
) -> Result<fs_blake3::ContentHash, SpatialUncertaintyError> {
    let mut encoder = IdentityEncoder::new(b"registration-model");
    encoder.u32(REGISTRATION_UNCERTAINTY_SCHEMA_VERSION);
    encoder.bytes(b"global-unit-circle-gls-v1");
    encoder.bytes(model.calibration_identity.as_bytes());
    encoder.u64(fiducials.len() as u64);
    encoder.u64((fiducials.len() * 2 - 3) as u64);
    match model.cross_fiducial {
        CrossFiducialModel::Independent => encoder.u8(0),
        CrossFiducialModel::EquicorrelatedStandardized { rho } => {
            encoder.u8(1);
            encoder.f64(rho);
        }
        CrossFiducialModel::Unknown => encoder.u8(2),
    }
    match model.huber {
        HuberPolicy::Disabled => encoder.u8(0),
        HuberPolicy::Enabled {
            threshold,
            iterations,
        } => {
            encoder.u8(1);
            encoder.f64(threshold);
            encoder.u8(iterations);
        }
    }
    match model.registered_inspection_bias {
        BiasBound::Bounded(bound) => {
            encoder.u8(0);
            encoder.f64(bound);
        }
        BiasBound::Unbounded => encoder.u8(1),
    }
    for (index, ((fiducial, covariance), (diagnostic, leverage))) in fiducials
        .iter()
        .zip(&model.fiducial_covariances)
        .zip(diagnostics.iter().zip(leverage))
        .enumerate()
    {
        checkpoint(cx, index, "registration model identity")?;
        encoder.point(fiducial.design());
        encoder.point(fiducial.measured());
        encoder.covariance(*covariance);
        encoder.f64(diagnostic.standardized_residual_norm);
        encoder.f64(diagnostic.robust_weight);
        encoder.u8(match diagnostic.disposition {
            OutlierDisposition::NotEvaluated => 0,
            OutlierDisposition::Retained => 1,
            OutlierDisposition::Downweighted => 2,
        });
        encoder.f64(*leverage);
    }
    encoder.f64(registration.tx());
    encoder.f64(registration.ty());
    encoder.f64(registration.rotation_rad());
    encoder.f64(registration.residual_rms());
    for value in covariance.iter().flatten() {
        encoder.f64(*value);
    }
    cx.checkpoint()
        .map_err(|_| SpatialUncertaintyError::Cancelled {
            phase: "registration model identity",
        })?;
    Ok(encoder.finish())
}

/// Interpretation of one fiducial's robust outlier screen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutlierDisposition {
    /// Robust screening was disabled for this calibrated model.
    NotEvaluated,
    /// The robust screen retained this fiducial at unit weight.
    Retained,
    /// The robust screen reduced this fiducial's influence.
    Downweighted,
}

/// Read-only diagnostic for one fiducial in input order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OutlierDiagnostic {
    standardized_residual_norm: f64,
    robust_weight: f64,
    disposition: OutlierDisposition,
}

impl OutlierDiagnostic {
    /// Residual norm after the declared covariance whitening.
    #[must_use]
    pub const fn standardized_residual_norm(&self) -> f64 {
        self.standardized_residual_norm
    }

    /// Final fixed robust weight used for the published point estimate.
    #[must_use]
    pub const fn robust_weight(&self) -> f64 {
        self.robust_weight
    }

    /// Stable robust-screen classification.
    #[must_use]
    pub const fn disposition(&self) -> OutlierDisposition {
        self.disposition
    }
}

/// Robustly refitted registration plus the full first-order covariance of
/// `(tx, ty, theta)`. For Huber fits the covariance is a frozen-weight
/// sandwich: it is conditional/first-order and does not claim to cover weight
/// selection uncertainty.
#[derive(Debug, Clone, PartialEq)]
pub struct CalibratedRegistration {
    registration: Registration,
    covariance: Mat3,
    degrees_of_freedom: usize,
    weights: Vec<f64>,
    leverage: Vec<f64>,
    outlier_diagnostics: Vec<OutlierDiagnostic>,
    model_identity: fs_blake3::ContentHash,
    registered_inspection_bias: BiasBound,
    robust_conditional: bool,
}

impl CalibratedRegistration {
    /// Robust/GLS point estimate.
    #[must_use]
    pub const fn registration(&self) -> &Registration {
        &self.registration
    }

    /// Row-major covariance for `(tx, ty, theta)`.
    #[must_use]
    pub const fn covariance(&self) -> &TransformCovariance {
        &self.covariance
    }

    /// Residual degrees of freedom, exactly `2n - 3`.
    #[must_use]
    pub const fn degrees_of_freedom(&self) -> usize {
        self.degrees_of_freedom
    }

    /// Final deterministic per-fiducial robust weights.
    #[must_use]
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// Full-model observation-block leverage traces. Correlated GLS blocks may
    /// be signed because the hat operator is not symmetric in observation
    /// coordinates; their sum is the effective fitted parameter dimension.
    #[must_use]
    pub fn leverage(&self) -> &[f64] {
        &self.leverage
    }

    /// Per-fiducial standardized residual and robust classification.
    #[must_use]
    pub fn outlier_diagnostics(&self) -> &[OutlierDiagnostic] {
        &self.outlier_diagnostics
    }

    /// Domain-separated identity binding model, inputs, fit, covariance, and
    /// diagnostics. It is an integrity address, not authentication.
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

/// Estimate a robust/GLS rigid transform and its complete first-order
/// fixed-weight parameter covariance. Every fixed-weight fit is the global
/// solution of a constrained `(tx, ty, cos(theta), sin(theta))` problem; there
/// is no caller-provided local starting angle.
///
/// # Errors
/// Invalid lengths/domains, unknown dependence, unresolved information,
/// cancellation, allocation failure, or non-finite arithmetic.
#[allow(clippy::too_many_lines)]
pub fn estimate_calibrated_registration(
    fiducials: &[Fiducial],
    model: &MetrologyModel,
    cx: &fs_exec::Cx<'_>,
) -> Result<CalibratedRegistration, SpatialUncertaintyError> {
    if fiducials.len() != model.fiducial_covariances.len() {
        return Err(SpatialUncertaintyError::LengthMismatch {
            field: "fiducial_covariances",
            expected: fiducials.len(),
            found: model.fiducial_covariances.len(),
        });
    }
    if fiducials.len() < 3 {
        return Err(SpatialUncertaintyError::TooFewFiducials {
            have: fiducials.len(),
        });
    }
    if matches!(model.cross_fiducial, CrossFiducialModel::Unknown) {
        return Err(SpatialUncertaintyError::UnknownDependence);
    }
    require_rank_two_design(fiducials, cx)?;
    let mut weights = Vec::new();
    weights
        .try_reserve_exact(fiducials.len())
        .map_err(|_| SpatialUncertaintyError::AllocationFailed)?;
    for index in 0..fiducials.len() {
        checkpoint(cx, index, "registration weight initialization")?;
        weights.push(1.0);
    }
    let mut solution = solve_global_fixed_weights(fiducials, model, &weights, cx)?;
    if let HuberPolicy::Enabled {
        threshold,
        iterations,
    } = model.huber
    {
        for _ in 0..iterations {
            let residuals = equations(fiducials, model, solution.parameters, &weights, cx)?;
            for (index, (weight, norm)) in weights
                .iter_mut()
                .zip(residuals.standardized_norms)
                .enumerate()
            {
                checkpoint(cx, index, "robust weight refresh")?;
                *weight = if norm <= threshold || norm == 0.0 {
                    1.0
                } else {
                    threshold / norm
                };
            }
            // This solve occurs after every refresh, including the last one:
            // the published transform and covariance therefore use the exact
            // same final fixed weights.
            solution = solve_global_fixed_weights(fiducials, model, &weights, cx)?;
        }
    }

    let parameters = solution.parameters;
    let final_system = equations(fiducials, model, parameters, &weights, cx)?;
    let mut sensitivity_information = final_system.information;
    sensitivity_information[2][2] += solution.trust_multiplier;
    let bread = invert_spd3(&sensitivity_information)?;
    let raw_covariance = if matches!(model.huber, HuberPolicy::Disabled) {
        mat3_mul(
            &mat3_mul(&bread, &final_system.information),
            &transpose3(&bread),
        )
    } else {
        robust_sandwich(&bread, &final_system.whitened_rows, &weights, cx)?
    };
    let covariance = symmetrize_and_validate_covariance(raw_covariance)?;
    let leverage = full_model_leverage(
        &bread,
        &final_system.whitened_rows,
        &weights,
        model.cross_fiducial,
        cx,
    )?;
    let mut outlier_diagnostics = Vec::new();
    outlier_diagnostics
        .try_reserve_exact(fiducials.len())
        .map_err(|_| SpatialUncertaintyError::AllocationFailed)?;
    for (index, (standardized_residual_norm, robust_weight)) in final_system
        .standardized_norms
        .iter()
        .zip(&weights)
        .enumerate()
    {
        checkpoint(cx, index, "outlier diagnostics")?;
        let disposition = match model.huber {
            HuberPolicy::Disabled => OutlierDisposition::NotEvaluated,
            HuberPolicy::Enabled { .. } if *robust_weight < 1.0 => OutlierDisposition::Downweighted,
            HuberPolicy::Enabled { .. } => OutlierDisposition::Retained,
        };
        outlier_diagnostics.push(OutlierDiagnostic {
            standardized_residual_norm: *standardized_residual_norm,
            robust_weight: *robust_weight,
            disposition,
        });
    }

    let mut residual_scale = 0.0f64;
    let mut residual_squares = 0.0f64;
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "registration residual diagnostic")?;
        let (predicted, _) = predicted_and_jacobian(*fiducial, parameters)?;
        let measured = fiducial.measured();
        let residual = (measured.x() - predicted.x()).hypot(measured.y() - predicted.y());
        if residual != 0.0 {
            if residual_scale < residual {
                let ratio = residual_scale / residual;
                residual_squares = 1.0 + residual_squares * ratio * ratio;
                residual_scale = residual;
            } else {
                let ratio = residual / residual_scale;
                residual_squares += ratio * ratio;
            }
        }
    }
    let residual_rms = if residual_scale == 0.0 {
        0.0
    } else {
        residual_scale * (residual_squares / fiducials.len() as f64).sqrt()
    };
    let registration = Registration::new(parameters[2], parameters[0], parameters[1], residual_rms)
        .map_err(|_| SpatialUncertaintyError::ArithmeticOverflow {
            field: "refitted registration",
        })?;
    let identity = model_identity(
        fiducials,
        model,
        &registration,
        &covariance,
        &outlier_diagnostics,
        &leverage,
        cx,
    )?;
    cx.checkpoint()
        .map_err(|_| SpatialUncertaintyError::Cancelled {
            phase: "registration publication",
        })?;
    Ok(CalibratedRegistration {
        registration,
        covariance,
        degrees_of_freedom: fiducials.len() * 2 - 3,
        weights,
        leverage,
        outlier_diagnostics,
        model_identity: identity,
        registered_inspection_bias: model.registered_inspection_bias,
        robust_conditional: !matches!(model.huber, HuberPolicy::Disabled),
    })
}

/// Relationship between the inspection measurements and the measurements that
/// fitted the registration. Reuse requires influence cross-covariance terms;
/// this v1 API therefore supports only a declared disjoint family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InspectionRelation {
    /// Inspection measurement errors are independent of the registration
    /// fiducial errors. Inspection points may remain mutually correlated: the
    /// simultaneous union bound does not assume otherwise.
    DisjointFromRegistration,
    /// A measurement is reused or otherwise dependent, but the required cross
    /// covariance was not supplied. No calibrated point bound is emitted.
    UnknownOrOverlapping,
}

/// Calibrated tolerance decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecisionState {
    /// The simultaneous upper bound fits the design tolerance.
    WithinTolerance,
    /// The simultaneous lower bound exceeds the design tolerance.
    ExceedsTolerance,
    /// Bounds overlap the tolerance or a required model component is absent.
    Indeterminate,
}

/// Why the decision has its current authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecisionReason {
    /// Finite simultaneous bounds establish the state.
    SimultaneousBounds,
    /// Finite simultaneous bounds straddle the tolerance.
    ToleranceOverlap,
    /// No finite systematic bias bound was supplied.
    UnboundedBias,
    /// Registration and inspection measurements may be dependent.
    RegistrationInspectionDependence,
    /// Robust covariance is conditional on selected weights and is not used as
    /// a finite-sample coverage certificate.
    AdaptiveWeightsConditional,
}

/// Per-point first-order propagation and simultaneous radial bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointDecisionBound {
    observed_deviation: f64,
    total_covariance: Covariance2,
    simultaneous_radius: f64,
    lower: f64,
    upper: f64,
}

impl PointDecisionBound {
    /// Observed registered-design to scan deviation.
    #[must_use]
    pub const fn observed_deviation(&self) -> f64 {
        self.observed_deviation
    }

    /// `G Cov(tx,ty,theta) G^T + Cov(inspection)`, with no residual-RMS term.
    #[must_use]
    pub const fn total_covariance(&self) -> Covariance2 {
        self.total_covariance
    }

    /// Familywise Chebyshev-union radial uncertainty plus bounded bias.
    #[must_use]
    pub const fn simultaneous_radius(&self) -> f64 {
        self.simultaneous_radius
    }

    /// Non-negative lower bound on this point's true deviation.
    #[must_use]
    pub const fn lower(&self) -> f64 {
        self.lower
    }

    /// Upper bound on this point's true deviation.
    #[must_use]
    pub const fn upper(&self) -> f64 {
        self.upper
    }
}

/// Evidence-bearing familywise decision. Bounds are absent exactly when a
/// required dependence, bias, or robust-selection model is unavailable.
#[derive(Debug, Clone, PartialEq)]
pub struct DecisionEvidence {
    state: DecisionState,
    reason: DecisionReason,
    lower: Option<f64>,
    upper: Option<f64>,
    tolerance: f64,
    confidence: f64,
    family_size: usize,
}

impl DecisionEvidence {
    /// Tri-state decision.
    #[must_use]
    pub const fn state(&self) -> DecisionState {
        self.state
    }

    /// Stable authority/no-claim reason.
    #[must_use]
    pub const fn reason(&self) -> DecisionReason {
        self.reason
    }

    /// Simultaneous maximum-deviation lower bound, when available.
    #[must_use]
    pub const fn lower(&self) -> Option<f64> {
        self.lower
    }

    /// Simultaneous maximum-deviation upper bound, when available.
    #[must_use]
    pub const fn upper(&self) -> Option<f64> {
        self.upper
    }

    /// Design tolerance used for the decision.
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Requested familywise confidence in `(0,1)`.
    #[must_use]
    pub const fn confidence(&self) -> f64 {
        self.confidence
    }

    /// Number of simultaneously covered inspection points.
    #[must_use]
    pub const fn family_size(&self) -> usize {
        self.family_size
    }
}

/// Candidate spatial evidence. Private fields prevent post-construction
/// mutation; its identity is content integrity only until admitted through an
/// injected [`EvidenceVerifier`]. Authentication never changes its scientific
/// rank or converts first-order statistical assumptions into validation.
#[derive(Debug, Clone, PartialEq)]
pub struct CalibratedAsBuiltEvidence {
    point_bounds: Vec<PointDecisionBound>,
    decision: DecisionEvidence,
    model_identity: fs_blake3::ContentHash,
    evidence_identity: fs_blake3::ContentHash,
}

impl CalibratedAsBuiltEvidence {
    /// Pointwise bounds in inspection input order. Empty exactly for an
    /// unavailable/no-claim dependence, bias, or robust-selection model.
    #[must_use]
    pub fn point_bounds(&self) -> &[PointDecisionBound] {
        &self.point_bounds
    }

    /// Familywise tri-state evidence.
    #[must_use]
    pub const fn decision(&self) -> &DecisionEvidence {
        &self.decision
    }

    /// Registration-model identity consumed by the assessment.
    #[must_use]
    pub const fn model_identity(&self) -> fs_blake3::ContentHash {
        self.model_identity
    }

    /// Full content identity of model, confidence, inputs, and outputs.
    #[must_use]
    pub const fn evidence_identity(&self) -> fs_blake3::ContentHash {
        self.evidence_identity
    }

    /// Ask an injected authority to authenticate this exact candidate and
    /// receipt. The candidate remains an estimated/conditional scientific
    /// result; this wrapper authenticates lineage, not truth.
    pub fn authenticate(
        self,
        receipt: EvidenceReceipt,
        verifier: &dyn EvidenceVerifier,
    ) -> Result<AuthenticatedAsBuiltEvidence, EvidenceAuthenticationError> {
        if receipt.schema_version != SPATIAL_EVIDENCE_SCHEMA_VERSION {
            return Err(EvidenceAuthenticationError::SchemaMismatch {
                receipt: receipt.schema_version,
                current: SPATIAL_EVIDENCE_SCHEMA_VERSION,
            });
        }
        if receipt.evidence_identity != self.evidence_identity {
            return Err(EvidenceAuthenticationError::EvidenceMismatch);
        }
        let verification = verifier.verify(&self, &receipt);
        if !verification.accepted {
            return Err(EvidenceAuthenticationError::Refused {
                policy: verification.policy,
            });
        }
        if verification.policy != receipt.policy_fingerprint {
            return Err(EvidenceAuthenticationError::PolicyMismatch {
                receipt: receipt.policy_fingerprint,
                decision: verification.policy,
            });
        }
        Ok(AuthenticatedAsBuiltEvidence {
            evidence: self,
            receipt,
        })
    }
}

#[derive(Debug, Clone, Copy)]
struct PosePropagation {
    covariance: Covariance2,
    trace_enclosure: Interval,
    rotated_design: [Interval; 2],
}

fn pose_covariance(
    registration: &CalibratedRegistration,
    design: Point2,
) -> Result<PosePropagation, SpatialUncertaintyError> {
    let (sine, cosine) = registration.registration.rotation_rad().sin_cos();
    let rotated_x = cosine * design.x() - sine * design.y();
    let rotated_y = sine * design.x() + cosine * design.y();
    let rows = [[1.0, 0.0, -rotated_y], [0.0, 1.0, rotated_x]];
    let angle = Interval::point(registration.registration.rotation_rad());
    let interval_sine = angle.sin();
    let interval_cosine = angle.cos();
    let design_x = Interval::point(design.x());
    let design_y = Interval::point(design.y());
    let interval_rotated_x = interval_cosine * design_x - interval_sine * design_y;
    let interval_rotated_y = interval_sine * design_x + interval_cosine * design_y;
    let interval_rows = [
        [
            Interval::point(1.0),
            Interval::point(0.0),
            -interval_rotated_y,
        ],
        [
            Interval::point(0.0),
            Interval::point(1.0),
            interval_rotated_x,
        ],
    ];
    let covariance = &registration.covariance;
    let quadratic = |left: [f64; 3], right: [f64; 3]| {
        let projected = mat3_vec(covariance, right);
        left[0] * projected[0] + left[1] * projected[1] + left[2] * projected[2]
    };
    let quadratic_interval = |vector: [Interval; 3]| {
        let mut sum = Interval::point(0.0);
        for row in 0..3 {
            for column in 0..3 {
                sum = sum + vector[row] * Interval::point(covariance[row][column]) * vector[column];
            }
        }
        sum
    };
    let pose = Covariance2::new(
        quadratic(rows[0], rows[0]),
        0.5 * (quadratic(rows[0], rows[1]) + quadratic(rows[1], rows[0])),
        quadratic(rows[1], rows[1]),
    )?;
    let trace = quadratic_interval(interval_rows[0]) + quadratic_interval(interval_rows[1]);
    if !trace.lo().is_finite() || !trace.hi().is_finite() || trace.hi() < 0.0 {
        return Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "pose covariance trace enclosure",
        });
    }
    Ok(PosePropagation {
        covariance: pose,
        trace_enclosure: Interval::new(trace.lo().max(0.0), trace.hi()),
        rotated_design: [interval_rotated_x, interval_rotated_y],
    })
}

fn unavailable_decision(
    reason: DecisionReason,
    tolerance: f64,
    confidence: f64,
    family_size: usize,
) -> DecisionEvidence {
    DecisionEvidence {
        state: DecisionState::Indeterminate,
        reason,
        lower: None,
        upper: None,
        tolerance,
        confidence,
        family_size,
    }
}

#[allow(clippy::too_many_arguments)]
fn assessment_identity(
    registration: &CalibratedRegistration,
    design: &[Point2],
    scanned: &[Point2],
    inspection_covariances: &[Covariance2],
    relation: InspectionRelation,
    tolerance: f64,
    confidence: f64,
    point_bounds: &[PointDecisionBound],
    decision: &DecisionEvidence,
    cx: &fs_exec::Cx<'_>,
) -> Result<fs_blake3::ContentHash, SpatialUncertaintyError> {
    let mut encoder = IdentityEncoder::new(b"spatial-decision");
    encoder.u32(SPATIAL_EVIDENCE_SCHEMA_VERSION);
    encoder.bytes(registration.model_identity.as_bytes());
    encoder.u64(design.len() as u64);
    encoder.u8(match relation {
        InspectionRelation::DisjointFromRegistration => 0,
        InspectionRelation::UnknownOrOverlapping => 1,
    });
    encoder.f64(tolerance);
    encoder.f64(confidence);
    for (index, ((design, scanned), covariance)) in design
        .iter()
        .zip(scanned)
        .zip(inspection_covariances)
        .enumerate()
    {
        checkpoint(cx, index, "spatial evidence input identity")?;
        encoder.point(*design);
        encoder.point(*scanned);
        encoder.covariance(*covariance);
    }
    encoder.u64(point_bounds.len() as u64);
    for (index, bound) in point_bounds.iter().enumerate() {
        checkpoint(cx, index, "spatial evidence output identity")?;
        encoder.f64(bound.observed_deviation);
        encoder.covariance(bound.total_covariance);
        encoder.f64(bound.simultaneous_radius);
        encoder.f64(bound.lower);
        encoder.f64(bound.upper);
    }
    encoder.u8(match decision.state {
        DecisionState::WithinTolerance => 0,
        DecisionState::ExceedsTolerance => 1,
        DecisionState::Indeterminate => 2,
    });
    encoder.u8(match decision.reason {
        DecisionReason::SimultaneousBounds => 0,
        DecisionReason::ToleranceOverlap => 1,
        DecisionReason::UnboundedBias => 2,
        DecisionReason::RegistrationInspectionDependence => 3,
        DecisionReason::AdaptiveWeightsConditional => 4,
    });
    match decision.lower {
        Some(value) => {
            encoder.u8(1);
            encoder.f64(value);
        }
        None => encoder.u8(0),
    }
    match decision.upper {
        Some(value) => {
            encoder.u8(1);
            encoder.f64(value);
        }
        None => encoder.u8(0),
    }
    cx.checkpoint()
        .map_err(|_| SpatialUncertaintyError::Cancelled {
            phase: "spatial evidence identity",
        })?;
    Ok(encoder.finish())
}

/// Propagate pose covariance to every inspected point and compute one
/// familywise maximum-deviation decision.
///
/// # Errors
/// Empty/mismatched/oversized families, invalid tolerance/confidence,
/// cancellation, allocation failure, or non-finite propagation.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn assess_calibrated_as_built(
    registration: &CalibratedRegistration,
    design: &[Point2],
    scanned: &[Point2],
    inspection_covariances: &[Covariance2],
    relation: InspectionRelation,
    tolerance: f64,
    confidence: f64,
    cx: &fs_exec::Cx<'_>,
) -> Result<CalibratedAsBuiltEvidence, SpatialUncertaintyError> {
    if design.is_empty() {
        return Err(SpatialUncertaintyError::InvalidScalar {
            field: "inspection.family_size",
            requirement: "positive",
        });
    }
    if design.len() > MAX_AS_BUILT_POINTS {
        return Err(SpatialUncertaintyError::TooManyPoints {
            have: design.len(),
            max: MAX_AS_BUILT_POINTS,
        });
    }
    for (field, found) in [
        ("scanned", scanned.len()),
        ("inspection_covariances", inspection_covariances.len()),
    ] {
        if found != design.len() {
            return Err(SpatialUncertaintyError::LengthMismatch {
                field,
                expected: design.len(),
                found,
            });
        }
    }
    let tolerance = finite_non_negative("tolerance", tolerance)?;
    if !confidence.is_finite() || confidence <= 0.0 || confidence >= 1.0 {
        return Err(SpatialUncertaintyError::InvalidScalar {
            field: "confidence",
            requirement: "strictly between zero and one",
        });
    }

    let unavailable_reason = match (
        relation,
        registration.registered_inspection_bias,
        registration.robust_conditional,
    ) {
        (InspectionRelation::UnknownOrOverlapping, _, _) => {
            Some(DecisionReason::RegistrationInspectionDependence)
        }
        (_, BiasBound::Unbounded, _) => Some(DecisionReason::UnboundedBias),
        (_, _, true) => Some(DecisionReason::AdaptiveWeightsConditional),
        _ => None,
    };
    if let Some(reason) = unavailable_reason {
        let decision = unavailable_decision(reason, tolerance, confidence, design.len());
        let identity = assessment_identity(
            registration,
            design,
            scanned,
            inspection_covariances,
            relation,
            tolerance,
            confidence,
            &[],
            &decision,
            cx,
        )?;
        cx.checkpoint()
            .map_err(|_| SpatialUncertaintyError::Cancelled {
                phase: "spatial evidence publication",
            })?;
        return Ok(CalibratedAsBuiltEvidence {
            point_bounds: Vec::new(),
            decision,
            model_identity: registration.model_identity,
            evidence_identity: identity,
        });
    }
    let BiasBound::Bounded(bias) = registration.registered_inspection_bias else {
        unreachable!("unbounded bias returned above")
    };
    let alpha = Interval::point(1.0) - Interval::point(confidence);
    let family_multiplier = Interval::point(design.len() as f64) / alpha;
    if !family_multiplier.lo().is_finite() || !family_multiplier.hi().is_finite() {
        return Err(SpatialUncertaintyError::ArithmeticOverflow {
            field: "simultaneous family multiplier",
        });
    }
    let mut bounds = Vec::new();
    bounds
        .try_reserve_exact(design.len())
        .map_err(|_| SpatialUncertaintyError::AllocationFailed)?;
    let mut maximum_lower = 0.0f64;
    let mut maximum_upper = 0.0f64;
    for (index, ((design_point, scanned_point), inspection_covariance)) in design
        .iter()
        .zip(scanned)
        .zip(inspection_covariances)
        .enumerate()
    {
        checkpoint(cx, index, "spatial propagation")?;
        let pose = pose_covariance(registration, *design_point)?;
        let total = Covariance2::new(
            pose.covariance.xx + inspection_covariance.xx,
            pose.covariance.xy + inspection_covariance.xy,
            pose.covariance.yy + inspection_covariance.yy,
        )?;
        let mapped = registration
            .registration
            .apply(*design_point)
            .map_err(|_| SpatialUncertaintyError::ArithmeticOverflow {
                field: "registered inspection point",
            })?;
        let mapped_x_enclosure =
            pose.rotated_design[0] + Interval::point(registration.registration.tx());
        let mapped_y_enclosure =
            pose.rotated_design[1] + Interval::point(registration.registration.ty());
        let dx = mapped_x_enclosure - Interval::point(scanned_point.x());
        let dy = mapped_y_enclosure - Interval::point(scanned_point.y());
        let deviation_enclosure = (dx * dx + dy * dy).sqrt();
        let deviation = (mapped.x() - scanned_point.x()).hypot(mapped.y() - scanned_point.y());
        let inspection_trace =
            Interval::point(inspection_covariance.xx) + Interval::point(inspection_covariance.yy);
        let stochastic_square = (pose.trace_enclosure + inspection_trace) * family_multiplier;
        if !deviation.is_finite()
            || !deviation_enclosure.lo().is_finite()
            || !deviation_enclosure.hi().is_finite()
            || !stochastic_square.lo().is_finite()
            || !stochastic_square.hi().is_finite()
            || stochastic_square.hi() < 0.0
        {
            return Err(SpatialUncertaintyError::ArithmeticOverflow {
                field: "simultaneous point radius",
            });
        }
        let stochastic_square =
            Interval::new(stochastic_square.lo().max(0.0), stochastic_square.hi());
        let radius_enclosure = Interval::point(bias) + stochastic_square.sqrt();
        let radius = radius_enclosure.hi();
        let lower = (deviation_enclosure - radius_enclosure).lo().max(0.0);
        let upper = (deviation_enclosure + radius_enclosure).hi();
        if !radius.is_finite() || !upper.is_finite() {
            return Err(SpatialUncertaintyError::ArithmeticOverflow {
                field: "simultaneous point bound",
            });
        }
        maximum_lower = maximum_lower.max(lower);
        maximum_upper = maximum_upper.max(upper);
        bounds.push(PointDecisionBound {
            observed_deviation: deviation,
            total_covariance: total,
            simultaneous_radius: radius,
            lower,
            upper,
        });
    }
    cx.checkpoint()
        .map_err(|_| SpatialUncertaintyError::Cancelled {
            phase: "spatial propagation",
        })?;
    let (state, reason) = if maximum_upper <= tolerance {
        (
            DecisionState::WithinTolerance,
            DecisionReason::SimultaneousBounds,
        )
    } else if maximum_lower > tolerance {
        (
            DecisionState::ExceedsTolerance,
            DecisionReason::SimultaneousBounds,
        )
    } else {
        (
            DecisionState::Indeterminate,
            DecisionReason::ToleranceOverlap,
        )
    };
    let decision = DecisionEvidence {
        state,
        reason,
        lower: Some(maximum_lower),
        upper: Some(maximum_upper),
        tolerance,
        confidence,
        family_size: design.len(),
    };
    let identity = assessment_identity(
        registration,
        design,
        scanned,
        inspection_covariances,
        relation,
        tolerance,
        confidence,
        &bounds,
        &decision,
        cx,
    )?;
    cx.checkpoint()
        .map_err(|_| SpatialUncertaintyError::Cancelled {
            phase: "spatial evidence publication",
        })?;
    Ok(CalibratedAsBuiltEvidence {
        point_bounds: bounds,
        decision,
        model_identity: registration.model_identity,
        evidence_identity: identity,
    })
}

/// Plain receipt data. It becomes authoritative only when an injected
/// [`EvidenceVerifier`] authenticates it against the exact candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceReceipt {
    evidence_identity: fs_blake3::ContentHash,
    schema_version: u32,
    policy_fingerprint: fs_blake3::ContentHash,
}

impl EvidenceReceipt {
    /// Assemble a receipt from retained ledger/policy fields.
    #[must_use]
    pub fn from_parts(
        evidence_identity: fs_blake3::ContentHash,
        schema_version: u32,
        policy_fingerprint: fs_blake3::ContentHash,
    ) -> Self {
        Self {
            evidence_identity,
            schema_version,
            policy_fingerprint,
        }
    }

    /// Bound evidence root.
    #[must_use]
    pub const fn evidence_identity(&self) -> fs_blake3::ContentHash {
        self.evidence_identity
    }

    /// Bound evidence schema.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Bound admitting-policy identity.
    #[must_use]
    pub const fn policy_fingerprint(&self) -> fs_blake3::ContentHash {
        self.policy_fingerprint
    }
}

/// Verifier decision with the exact policy identity that decided.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceVerification {
    accepted: bool,
    policy: fs_blake3::ContentHash,
}

impl EvidenceVerification {
    /// Accept under a named policy.
    #[must_use]
    pub fn accept(policy: fs_blake3::ContentHash) -> Self {
        Self {
            accepted: true,
            policy,
        }
    }

    /// Reject under a named policy.
    #[must_use]
    pub fn reject(policy: fs_blake3::ContentHash) -> Self {
        Self {
            accepted: false,
            policy,
        }
    }
}

/// Injected capability that authenticates one candidate/receipt pair against
/// an external ledger or policy. Implementations are explicit trust roots.
pub trait EvidenceVerifier {
    /// Authenticate the pair and name the deciding policy.
    fn verify(
        &self,
        evidence: &CalibratedAsBuiltEvidence,
        receipt: &EvidenceReceipt,
    ) -> EvidenceVerification;
}

/// Deny-all default: an integrity hash never authenticates itself.
#[derive(Debug, Clone, Copy, Default)]
pub struct NoEvidenceVerifier;

/// Stable deny-all policy identity.
#[must_use]
pub fn no_evidence_policy() -> fs_blake3::ContentHash {
    fs_blake3::hash_domain(DENY_POLICY_DOMAIN, b"deny-all")
}

impl EvidenceVerifier for NoEvidenceVerifier {
    fn verify(
        &self,
        _evidence: &CalibratedAsBuiltEvidence,
        _receipt: &EvidenceReceipt,
    ) -> EvidenceVerification {
        EvidenceVerification::reject(no_evidence_policy())
    }
}

/// Authentication refusal at the local capability boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvidenceAuthenticationError {
    /// Receipt schema is stale or unknown.
    SchemaMismatch {
        /// Receipt version.
        receipt: u32,
        /// Current version.
        current: u32,
    },
    /// Receipt names a different evidence root.
    EvidenceMismatch,
    /// Verifier refused the pair.
    Refused {
        /// Deciding policy.
        policy: fs_blake3::ContentHash,
    },
    /// Accepting verifier named a policy different from the receipt.
    PolicyMismatch {
        /// Receipt policy.
        receipt: fs_blake3::ContentHash,
        /// Verifier decision policy.
        decision: fs_blake3::ContentHash,
    },
}

impl core::fmt::Display for EvidenceAuthenticationError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SchemaMismatch { receipt, current } => write!(
                formatter,
                "spatial evidence schema {receipt} does not match current {current}"
            ),
            Self::EvidenceMismatch => formatter.write_str("spatial evidence root mismatch"),
            Self::Refused { policy } => {
                write!(formatter, "spatial evidence refused by policy {policy}")
            }
            Self::PolicyMismatch { receipt, decision } => write!(
                formatter,
                "spatial evidence receipt policy {receipt} differs from decision {decision}"
            ),
        }
    }
}

impl std::error::Error for EvidenceAuthenticationError {}

/// Opaque lineage-authenticated candidate. Private fields prevent callers from
/// minting it without an accepting verifier.
#[derive(Debug, Clone, PartialEq)]
pub struct AuthenticatedAsBuiltEvidence {
    evidence: CalibratedAsBuiltEvidence,
    receipt: EvidenceReceipt,
}

impl AuthenticatedAsBuiltEvidence {
    /// Exact candidate authenticated by the external policy.
    #[must_use]
    pub const fn evidence(&self) -> &CalibratedAsBuiltEvidence {
        &self.evidence
    }

    /// Exact receipt authenticated with it.
    #[must_use]
    pub const fn receipt(&self) -> &EvidenceReceipt {
        &self.receipt
    }
}

struct IdentityEncoder {
    hasher: fs_blake3::Blake3,
}

impl IdentityEncoder {
    fn new(kind: &[u8]) -> Self {
        let mut encoder = Self {
            hasher: fs_blake3::Blake3::new(),
        };
        encoder.bytes(IDENTITY_DOMAIN.as_bytes());
        encoder.bytes(kind);
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
        let value = if value == 0.0 { 0.0 } else { value };
        self.bytes(&value.to_bits().to_le_bytes());
    }

    fn point(&mut self, point: Point2) {
        self.f64(point.x());
        self.f64(point.y());
    }

    fn covariance(&mut self, covariance: Covariance2) {
        self.f64(covariance.xx);
        self.f64(covariance.xy);
        self.f64(covariance.yy);
    }

    fn finish(self) -> fs_blake3::ContentHash {
        let preimage = self.hasher.finalize();
        fs_blake3::hash_domain(IDENTITY_DOMAIN, preimage.as_bytes())
    }
}

#[cfg(test)]
mod numerical_admission_tests {
    use super::{
        SpatialUncertaintyError, minimize_on_unit_circle, symmetrize_and_validate_covariance,
    };

    #[test]
    fn outward_alignment_gate_refuses_rotated_trust_region_hard_case() {
        assert_eq!(
            minimize_on_unit_circle([[8.0, 6.0], [6.0, 17.0]], [1.0, 2.0]),
            Err(SpatialUncertaintyError::AmbiguousGlobalMinimum)
        );
        assert_eq!(
            minimize_on_unit_circle([[1.0, 0.0], [0.0, 2.0]], [0.0, 0.0]),
            Err(SpatialUncertaintyError::AmbiguousGlobalMinimum)
        );
    }

    #[test]
    fn covariance_publication_is_bit_symmetric_and_spd_only() {
        let covariance = symmetrize_and_validate_covariance([
            [4.0e12, 2.0, -3.0e-3],
            [2.0 + 4.0 * f64::EPSILON, 9.0, 1.0e-6],
            [-3.0e-3, 1.0e-6, 2.5e-8],
        ])
        .expect("anisotropic SPD covariance remains admissible");
        for row in 0..3 {
            for column in 0..3 {
                assert_eq!(
                    covariance[row][column].to_bits(),
                    covariance[column][row].to_bits()
                );
            }
        }
        assert_eq!(
            symmetrize_and_validate_covariance(
                [[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0],]
            ),
            Err(SpatialUncertaintyError::SingularInformation)
        );
    }
}
