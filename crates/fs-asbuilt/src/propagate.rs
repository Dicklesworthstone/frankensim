//! Registration pose covariance propagated into quantity-of-interest
//! uncertainty: the geometry term of the engineering error budget.
//!
//! The crate CONTRACT has always warned that `residual_rms` must not be
//! misused as transform covariance. This module is the other half of that
//! bargain: the CALIBRATED 6-dof pose covariance flows, typed, into
//! fs-evidence's eight-term engineering budget as the `Geometry` source.
//!
//! One transform moves everything, so pose-induced QoI errors are correlated
//! across QoIs. The propagation therefore emits a single content-addressed
//! record carrying the full cross-QoI covariance; every per-QoI geometry
//! term cites that shared record as its provenance/replay artifact, so the
//! correlation structure travels with the budget instead of degrading into
//! independent per-QoI noise.
//!
//! Method honesty: propagation is FIRST-ORDER (linearized through
//! caller-supplied pose sensitivities). When the caller can evaluate the
//! true QoI map at perturbed poses, a deterministic sigma-point spot-check
//! compares the nonlinear deltas against the linear prediction; if the
//! declared tolerance is exceeded the geometry terms DOWNGRADE to
//! `TermValue::Unknown` with a named reason rather than publishing a
//! confidently wrong half-width. Without an evaluator, the method records
//! `LinearizedUnchecked` and the declared reason — silence is not an option.

#![allow(clippy::needless_range_loop)] // Fixed 6-dof indices expose the pose ordering.
#![allow(clippy::float_cmp)] // Exact zeros distinguish structural refusals from IEEE noise.

use crate::canonical_zero;
use crate::rigid3::{CalibratedRigid3Registration, PoseCovariance6};
use fs_evidence::uncertainty::{
    DistributionTerm, EngineeringUncertaintyKind, EngineeringUncertaintyTerm, TermValue,
    UncertaintyArtifactRef,
};

/// Identity schema for geometry-propagation records.
pub const GEOMETRY_PROPAGATION_SCHEMA_VERSION: u32 = 1;
/// Maximum quantities of interest in one propagation.
pub const MAX_PROPAGATION_QOIS: usize = 64;
/// Number of deterministic sigma-point samples used by the spot-check
/// (one +/- pair per pose degree of freedom).
pub const SPOT_CHECK_SAMPLES: u32 = 12;

const IDENTITY_DOMAIN: &str = "org.frankensim.fs-asbuilt.geometry-propagation.v1";
const IDENTITY_KIND: &[u8] = b"geometry-propagation-record-v1";
const PROVENANCE_ROLE: &str = "as-built-pose-propagation";
const POLL_STRIDE: usize = 256;

/// A structured propagation failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PropagateError {
    /// No quantities of interest were declared.
    EmptyQois,
    /// The bounded QoI cap was exceeded.
    TooManyQois {
        /// Supplied count.
        have: usize,
        /// Maximum accepted count.
        max: usize,
    },
    /// Two QoIs share a name; the record would be ambiguous.
    DuplicateQoi {
        /// Ordinal of the second occurrence.
        index: usize,
    },
    /// A scalar violates a named finite/range invariant.
    InvalidScalar {
        /// Stable field name.
        field: &'static str,
        /// Stable required domain.
        requirement: &'static str,
    },
    /// A QoI name or unit failed the bounded-identifier grammar.
    InvalidName {
        /// Ordinal of the offending sensitivity.
        index: usize,
    },
    /// The evaluator returned the wrong number of QoI deltas.
    EvaluatorLengthMismatch {
        /// Required length.
        expected: usize,
        /// Returned length.
        found: usize,
    },
    /// The evaluator refused or returned a non-finite delta.
    EvaluatorRefused {
        /// Sigma-point ordinal (0-based).
        sample: usize,
    },
    /// The pose covariance could not be factored for sampling.
    NonPositiveDefinitePoseCovariance,
    /// Term admission failed inside fs-evidence (should be unreachable for
    /// validated inputs; surfaced instead of panicking).
    TermAdmission,
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

impl core::fmt::Display for PropagateError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyQois => write!(formatter, "no quantities of interest declared"),
            Self::TooManyQois { have, max } => {
                write!(
                    formatter,
                    "{have} quantities of interest exceed the cap {max}"
                )
            }
            Self::DuplicateQoi { index } => {
                write!(
                    formatter,
                    "duplicate quantity-of-interest name at ordinal {index}"
                )
            }
            Self::InvalidScalar { field, requirement } => {
                write!(formatter, "{field} must be {requirement}")
            }
            Self::InvalidName { index } => {
                write!(formatter, "sensitivity {index} has an invalid name or unit")
            }
            Self::EvaluatorLengthMismatch { expected, found } => {
                write!(
                    formatter,
                    "evaluator returned {found} deltas; expected {expected}"
                )
            }
            Self::EvaluatorRefused { sample } => {
                write!(formatter, "evaluator refused sigma point {sample}")
            }
            Self::NonPositiveDefinitePoseCovariance => {
                write!(
                    formatter,
                    "pose covariance is not positive definite for sampling"
                )
            }
            Self::TermAdmission => {
                write!(
                    formatter,
                    "fs-evidence refused a term built from validated inputs"
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

impl core::error::Error for PropagateError {}

fn checkpoint(
    cx: &fs_exec::Cx<'_>,
    ordinal: usize,
    phase: &'static str,
) -> Result<(), PropagateError> {
    if ordinal.is_multiple_of(POLL_STRIDE) {
        cx.checkpoint()
            .map_err(|_| PropagateError::Cancelled { phase })?;
    }
    Ok(())
}

fn bounded_name(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 128
        && value == value.trim()
        && value.bytes().all(|byte| {
            byte.is_ascii_alphanumeric()
                || matches!(byte, b'-' | b'_' | b'.' | b'/' | b':' | b'@' | b'+' | b'=')
        })
}

/// One quantity of interest's first-order pose sensitivity.
///
/// The gradient is `d QoI / d (tx, ty, tz, rx, ry, rz)` in the SAME pose
/// parameterization the calibrated registration's covariance uses: a left
/// rotation-vector perturbation about the weighted design centroid image.
/// Sensitivities computed about a different pivot are silently wrong; the
/// pivot convention is part of this type's contract.
#[derive(Debug, Clone, PartialEq)]
pub struct QoiSensitivity {
    name: String,
    unit: String,
    gradient: [f64; 6],
}

impl QoiSensitivity {
    /// Declare one QoI sensitivity.
    ///
    /// # Errors
    /// Non-finite gradient entries or an invalid name/unit grammar.
    pub fn new(
        name: impl Into<String>,
        unit: impl Into<String>,
        gradient: [f64; 6],
    ) -> Result<Self, PropagateError> {
        let name = name.into();
        let unit = unit.into();
        if !bounded_name(&name) || !bounded_name(&unit) {
            return Err(PropagateError::InvalidName { index: 0 });
        }
        for value in gradient {
            if !value.is_finite() {
                return Err(PropagateError::InvalidScalar {
                    field: "sensitivity gradient",
                    requirement: "finite",
                });
            }
        }
        let mut canonical = [0.0f64; 6];
        for axis in 0..6 {
            canonical[axis] = canonical_zero(gradient[axis]);
        }
        Ok(Self {
            name,
            unit,
            gradient: canonical,
        })
    }

    /// QoI name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// QoI unit label.
    #[must_use]
    pub fn unit(&self) -> &str {
        &self.unit
    }

    /// Pose gradient in `(tx, ty, tz, rx, ry, rz)` order.
    #[must_use]
    pub const fn gradient(&self) -> [f64; 6] {
        self.gradient
    }
}

/// Caller-supplied evaluation of the true QoI map at a perturbed pose,
/// used only by the linearization spot-check.
pub trait QoiEvaluator {
    /// Return the QoI DELTAS from the nominal-pose values under the given
    /// pose perturbation, one per declared sensitivity in declaration order.
    ///
    /// # Errors
    /// Any refusal aborts the propagation as [`PropagateError::EvaluatorRefused`].
    fn evaluate(&self, pose_delta: &[f64; 6]) -> Result<Vec<f64>, String>;
}

/// Caller-declared coverage policy for the published half-widths.
///
/// There is no default: the factor is the caller's declared engineering
/// policy (for example `k = 2` at `level = 0.95`), not a distributional
/// claim this module could certify.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoveragePolicy {
    level: f64,
    factor: f64,
}

impl CoveragePolicy {
    /// Declare a coverage policy.
    ///
    /// # Errors
    /// `level` outside `(0, 1)` or a non-positive/non-finite factor.
    pub fn new(level: f64, factor: f64) -> Result<Self, PropagateError> {
        if !level.is_finite() || level <= 0.0 || level >= 1.0 {
            return Err(PropagateError::InvalidScalar {
                field: "coverage.level",
                requirement: "strictly inside (0, 1)",
            });
        }
        if !factor.is_finite() || factor <= 0.0 {
            return Err(PropagateError::InvalidScalar {
                field: "coverage.factor",
                requirement: "finite and positive",
            });
        }
        Ok(Self { level, factor })
    }

    /// Declared coverage level.
    #[must_use]
    pub const fn level(&self) -> f64 {
        self.level
    }

    /// Declared half-width multiplier on the standard deviation.
    #[must_use]
    pub const fn factor(&self) -> f64 {
        self.factor
    }
}

/// How the propagation was performed and whether its linearization earned
/// a spot-check.
#[derive(Debug, Clone, PartialEq)]
pub enum PropagationMethod {
    /// Linearized propagation with no evaluator available; the declared
    /// reason travels with the record.
    LinearizedUnchecked {
        /// Why no spot-check was possible (never empty).
        reason: String,
    },
    /// Linearized propagation whose sigma-point spot-check stayed inside
    /// the declared tolerance.
    LinearizedSpotChecked {
        /// Deterministic sigma-point count.
        samples: u32,
        /// Worst relative gap observed across QoIs and sigma points.
        max_relative_gap: f64,
        /// The declared acceptance tolerance.
        tolerance: f64,
    },
    /// The spot-check exceeded tolerance: linearization is NOT valid at
    /// this covariance scale, and every geometry term downgrades to
    /// `Unknown` instead of publishing a half-width.
    LinearizationRejected {
        /// Deterministic sigma-point count.
        samples: u32,
        /// Worst relative gap observed.
        max_relative_gap: f64,
        /// The declared acceptance tolerance it exceeded.
        tolerance: f64,
        /// Ordinal of the worst-offending QoI.
        worst_qoi: usize,
    },
}

/// The content-addressed cross-QoI geometry propagation record.
#[derive(Debug, Clone, PartialEq)]
pub struct GeometryPropagation {
    qois: Vec<QoiSensitivity>,
    pose_covariance: PoseCovariance6,
    registration_identity: fs_blake3::ContentHash,
    cross_covariance: Vec<f64>,
    standard_deviations: Vec<f64>,
    coverage: CoveragePolicy,
    method: PropagationMethod,
    record_identity: fs_blake3::ContentHash,
}

impl GeometryPropagation {
    /// Declared sensitivities in record order.
    #[must_use]
    pub fn qois(&self) -> &[QoiSensitivity] {
        &self.qois
    }

    /// The pose covariance this record consumed.
    #[must_use]
    pub const fn pose_covariance(&self) -> &PoseCovariance6 {
        &self.pose_covariance
    }

    /// Model identity of the calibrated registration this record consumed.
    #[must_use]
    pub const fn registration_identity(&self) -> fs_blake3::ContentHash {
        self.registration_identity
    }

    /// Row-major cross-QoI covariance `G Sigma G^T` in record QoI order.
    /// This is the correlation structure the per-QoI budget terms share.
    #[must_use]
    pub fn cross_covariance(&self) -> &[f64] {
        &self.cross_covariance
    }

    /// Per-QoI standard deviations (square roots of the covariance
    /// diagonal), in record QoI order.
    #[must_use]
    pub fn standard_deviations(&self) -> &[f64] {
        &self.standard_deviations
    }

    /// Correlation coefficient between two record QoIs, or `None` when a
    /// marginal is exactly zero (correlation undefined, not assumed).
    #[must_use]
    pub fn correlation(&self, first: usize, second: usize) -> Option<f64> {
        let count = self.qois.len();
        if first >= count || second >= count {
            return None;
        }
        let denominator = self.standard_deviations[first] * self.standard_deviations[second];
        if denominator == 0.0 {
            return None;
        }
        Some(self.cross_covariance[first * count + second] / denominator)
    }

    /// The declared coverage policy.
    #[must_use]
    pub const fn coverage(&self) -> &CoveragePolicy {
        &self.coverage
    }

    /// How the propagation was performed.
    #[must_use]
    pub const fn method(&self) -> &PropagationMethod {
        &self.method
    }

    /// Content identity of the complete record. Every emitted geometry term
    /// cites this identity, which is how the shared correlation structure
    /// travels with per-QoI budgets.
    #[must_use]
    pub const fn record_identity(&self) -> fs_blake3::ContentHash {
        self.record_identity
    }

    /// Build the geometry budget term for one record QoI.
    ///
    /// A validated linearization yields a `Distribution` term (zero mean,
    /// propagated standard deviation, declared coverage half-width). A
    /// rejected or unchecked-but-declared-invalid state never reaches here
    /// as a distribution: when the method is `LinearizationRejected` the
    /// term is `Unknown` with a named reason. Both carry this record's
    /// identity as provenance and replay authority.
    ///
    /// # Errors
    /// An out-of-range ordinal or an internal admission refusal.
    pub fn geometry_term(
        &self,
        ordinal: usize,
    ) -> Result<EngineeringUncertaintyTerm, PropagateError> {
        if ordinal >= self.qois.len() {
            return Err(PropagateError::InvalidScalar {
                field: "qoi ordinal",
                requirement: "within the record's QoI count",
            });
        }
        let artifact = UncertaintyArtifactRef::new(PROVENANCE_ROLE, self.record_identity)
            .map_err(|_| PropagateError::TermAdmission)?;
        let value = match &self.method {
            PropagationMethod::LinearizationRejected {
                max_relative_gap,
                tolerance,
                worst_qoi,
                ..
            } => TermValue::unknown(format!(
                "linearization spot-check exceeded tolerance {tolerance} (worst relative gap {max_relative_gap} at qoi {}); a sampling-based propagation is required at this covariance scale",
                self.qois[*worst_qoi].name()
            ))
            .map_err(|_| PropagateError::TermAdmission)?,
            PropagationMethod::LinearizedUnchecked { .. }
            | PropagationMethod::LinearizedSpotChecked { .. } => {
                let deviation = self.standard_deviations[ordinal];
                TermValue::Distribution(DistributionTerm {
                    mean: 0.0,
                    standard_deviation: deviation,
                    conservative_half_width: self.coverage.factor * deviation,
                    level: self.coverage.level,
                    replay: artifact.clone(),
                })
            }
        };
        EngineeringUncertaintyTerm::try_new(EngineeringUncertaintyKind::Geometry, value, artifact)
            .map_err(|_| PropagateError::TermAdmission)
    }
}

/// Lower-triangular Cholesky factor of the pose covariance after diagonal
/// equilibration, used to build deterministic sigma points.
fn pose_cholesky(covariance: &PoseCovariance6) -> Result<[[f64; 6]; 6], PropagateError> {
    let mut scales = [0.0f64; 6];
    for index in 0..6 {
        let diagonal = covariance[index][index];
        if !diagonal.is_finite() || diagonal <= 0.0 {
            return Err(PropagateError::NonPositiveDefinitePoseCovariance);
        }
        scales[index] = diagonal.sqrt();
    }
    let mut lower = [[0.0f64; 6]; 6];
    let pivot_floor = 256.0 * f64::EPSILON;
    for row in 0..6 {
        for column in 0..=row {
            let mut value = covariance[row][column] / scales[row] / scales[column];
            for prior in 0..column {
                value -= lower[row][prior] * lower[column][prior];
            }
            if row == column {
                if !value.is_finite() || value <= pivot_floor {
                    return Err(PropagateError::NonPositiveDefinitePoseCovariance);
                }
                lower[row][column] = value.sqrt();
            } else {
                lower[row][column] = value / lower[column][column];
            }
        }
    }
    // Undo the equilibration so columns are true covariance factors.
    for row in 0..6 {
        for column in 0..6 {
            lower[row][column] *= scales[row];
        }
    }
    Ok(lower)
}

/// Propagate a calibrated registration's pose covariance into cross-QoI
/// geometry uncertainty.
///
/// Linearized: the cross-QoI covariance is `G Sigma G^T` for the declared
/// gradient matrix `G`. When `evaluator` is supplied, deterministic one-sigma
/// sigma points (`+/- L e_i` for the equilibrated Cholesky factor `L`) probe
/// the true QoI map; the worst relative gap against the linear prediction is
/// compared to `spot_check_tolerance` (relative to each QoI's propagated
/// standard deviation). Exceeding tolerance does not fail the call: the
/// record's method becomes `LinearizationRejected` and every geometry term
/// downgrades to `Unknown` with a named reason. Without an evaluator the
/// caller must declare why (`unchecked_reason`), which travels in the record.
///
/// # Errors
/// Invalid declarations, evaluator refusals/shape mismatches, a non-positive
/// definite pose covariance when sampling is requested, non-finite
/// aggregates, allocation failure, or structured cancellation.
#[allow(clippy::too_many_lines)]
pub fn propagate_pose_covariance(
    registration: &CalibratedRigid3Registration,
    sensitivities: &[QoiSensitivity],
    coverage: CoveragePolicy,
    evaluator: Option<&dyn QoiEvaluator>,
    spot_check_tolerance: f64,
    unchecked_reason: &str,
    cx: &fs_exec::Cx<'_>,
) -> Result<GeometryPropagation, PropagateError> {
    if sensitivities.is_empty() {
        return Err(PropagateError::EmptyQois);
    }
    if sensitivities.len() > MAX_PROPAGATION_QOIS {
        return Err(PropagateError::TooManyQois {
            have: sensitivities.len(),
            max: MAX_PROPAGATION_QOIS,
        });
    }
    for (index, sensitivity) in sensitivities.iter().enumerate() {
        checkpoint(cx, index, "propagate.declarations")?;
        for other in &sensitivities[..index] {
            if other.name == sensitivity.name {
                return Err(PropagateError::DuplicateQoi { index });
            }
        }
    }
    if evaluator.is_some() {
        if !spot_check_tolerance.is_finite() || spot_check_tolerance <= 0.0 {
            return Err(PropagateError::InvalidScalar {
                field: "spot_check_tolerance",
                requirement: "finite and positive when an evaluator is supplied",
            });
        }
    } else if unchecked_reason.trim().is_empty() {
        return Err(PropagateError::InvalidScalar {
            field: "unchecked_reason",
            requirement: "non-empty when no evaluator is supplied",
        });
    }

    let count = sensitivities.len();
    let covariance = registration.covariance();

    // Cross-QoI covariance: C = G Sigma G^T.
    let mut cross = Vec::new();
    cross
        .try_reserve_exact(count * count)
        .map_err(|_| PropagateError::AllocationFailed)?;
    // Precompute Sigma * g_j for each QoI.
    let mut sigma_gradients = Vec::new();
    sigma_gradients
        .try_reserve_exact(count)
        .map_err(|_| PropagateError::AllocationFailed)?;
    for (index, sensitivity) in sensitivities.iter().enumerate() {
        checkpoint(cx, index, "propagate.sigma-gradients")?;
        let gradient = sensitivity.gradient;
        let mut product = [0.0f64; 6];
        for row in 0..6 {
            let mut value = 0.0;
            for column in 0..6 {
                value = covariance[row][column].mul_add(gradient[column], value);
            }
            product[row] = value;
        }
        sigma_gradients.push(product);
    }
    for row in 0..count {
        checkpoint(cx, row, "propagate.cross-covariance")?;
        for column in 0..count {
            let mut value = 0.0;
            for axis in 0..6 {
                value =
                    sensitivities[row].gradient[axis].mul_add(sigma_gradients[column][axis], value);
            }
            if !value.is_finite() {
                return Err(PropagateError::ArithmeticOverflow {
                    field: "cross-QoI covariance",
                });
            }
            cross.push(value);
        }
    }
    // Midpoint-symmetrize so the record is exactly symmetric.
    for row in 0..count {
        for column in row + 1..count {
            let symmetric = f64::midpoint(cross[row * count + column], cross[column * count + row]);
            cross[row * count + column] = symmetric;
            cross[column * count + row] = symmetric;
        }
    }
    let mut deviations = Vec::new();
    deviations
        .try_reserve_exact(count)
        .map_err(|_| PropagateError::AllocationFailed)?;
    for index in 0..count {
        let variance = cross[index * count + index];
        if variance < 0.0 || !variance.is_finite() {
            return Err(PropagateError::ArithmeticOverflow {
                field: "propagated variance",
            });
        }
        deviations.push(variance.sqrt());
    }

    // Deterministic sigma-point spot-check of the linearization.
    let method = if let Some(evaluator) = evaluator {
        let lower = pose_cholesky(covariance)?;
        let mut max_relative_gap = 0.0f64;
        let mut worst_qoi = 0usize;
        let mut sample_ordinal = 0usize;
        for axis in 0..6 {
            for sign in [1.0f64, -1.0] {
                checkpoint(cx, sample_ordinal, "propagate.spot-check")?;
                let mut delta = [0.0f64; 6];
                for row in 0..6 {
                    delta[row] = sign * lower[row][axis];
                }
                let observed =
                    evaluator
                        .evaluate(&delta)
                        .map_err(|_| PropagateError::EvaluatorRefused {
                            sample: sample_ordinal,
                        })?;
                if observed.len() != count {
                    return Err(PropagateError::EvaluatorLengthMismatch {
                        expected: count,
                        found: observed.len(),
                    });
                }
                for (qoi, observed_delta) in observed.iter().enumerate() {
                    if !observed_delta.is_finite() {
                        return Err(PropagateError::EvaluatorRefused {
                            sample: sample_ordinal,
                        });
                    }
                    let mut predicted = 0.0f64;
                    for axis_inner in 0..6 {
                        predicted = sensitivities[qoi].gradient[axis_inner]
                            .mul_add(delta[axis_inner], predicted);
                    }
                    // Relative to the QoI's propagated one-sigma scale; a
                    // zero-scale QoI only accepts an exactly zero gap.
                    let scale = deviations[qoi];
                    let gap = (observed_delta - predicted).abs();
                    let relative = if scale > 0.0 {
                        gap / scale
                    } else if gap == 0.0 {
                        0.0
                    } else {
                        f64::INFINITY
                    };
                    if relative > max_relative_gap {
                        max_relative_gap = relative;
                        worst_qoi = qoi;
                    }
                }
                sample_ordinal += 1;
            }
        }
        if max_relative_gap > spot_check_tolerance {
            PropagationMethod::LinearizationRejected {
                samples: SPOT_CHECK_SAMPLES,
                max_relative_gap,
                tolerance: spot_check_tolerance,
                worst_qoi,
            }
        } else {
            PropagationMethod::LinearizedSpotChecked {
                samples: SPOT_CHECK_SAMPLES,
                max_relative_gap,
                tolerance: spot_check_tolerance,
            }
        }
    } else {
        PropagationMethod::LinearizedUnchecked {
            reason: unchecked_reason.trim().to_owned(),
        }
    };

    let record_identity = record_identity(
        sensitivities,
        covariance,
        registration.model_identity(),
        &cross,
        &deviations,
        coverage,
        &method,
        cx,
    )?;
    cx.checkpoint().map_err(|_| PropagateError::Cancelled {
        phase: "propagate.publish",
    })?;
    Ok(GeometryPropagation {
        qois: sensitivities.to_vec(),
        pose_covariance: *covariance,
        registration_identity: registration.model_identity(),
        cross_covariance: cross,
        standard_deviations: deviations,
        coverage,
        method,
        record_identity,
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

    fn finish(self) -> fs_blake3::ContentHash {
        let preimage = self.hasher.finalize();
        fs_blake3::hash_domain(IDENTITY_DOMAIN, preimage.as_bytes())
    }
}

#[allow(clippy::too_many_arguments)]
fn record_identity(
    sensitivities: &[QoiSensitivity],
    pose_covariance: &PoseCovariance6,
    registration_identity: fs_blake3::ContentHash,
    cross: &[f64],
    deviations: &[f64],
    coverage: CoveragePolicy,
    method: &PropagationMethod,
    cx: &fs_exec::Cx<'_>,
) -> Result<fs_blake3::ContentHash, PropagateError> {
    let mut encoder = IdentityEncoder::new();
    encoder.u32(GEOMETRY_PROPAGATION_SCHEMA_VERSION);
    encoder.bytes(registration_identity.as_bytes());
    encoder.u64(sensitivities.len() as u64);
    for (index, sensitivity) in sensitivities.iter().enumerate() {
        checkpoint(cx, index, "propagate.identity")?;
        encoder.bytes(sensitivity.name.as_bytes());
        encoder.bytes(sensitivity.unit.as_bytes());
        for value in sensitivity.gradient {
            encoder.f64(value);
        }
    }
    for row in pose_covariance {
        for value in row {
            encoder.f64(*value);
        }
    }
    for (index, value) in cross.iter().enumerate() {
        checkpoint(cx, index, "propagate.identity-covariance")?;
        encoder.f64(*value);
    }
    for value in deviations {
        encoder.f64(*value);
    }
    encoder.f64(coverage.level);
    encoder.f64(coverage.factor);
    match method {
        PropagationMethod::LinearizedUnchecked { reason } => {
            encoder.u8(0);
            encoder.bytes(reason.as_bytes());
        }
        PropagationMethod::LinearizedSpotChecked {
            samples,
            max_relative_gap,
            tolerance,
        } => {
            encoder.u8(1);
            encoder.u32(*samples);
            encoder.f64(*max_relative_gap);
            encoder.f64(*tolerance);
        }
        PropagationMethod::LinearizationRejected {
            samples,
            max_relative_gap,
            tolerance,
            worst_qoi,
        } => {
            encoder.u8(2);
            encoder.u32(*samples);
            encoder.f64(*max_relative_gap);
            encoder.f64(*tolerance);
            encoder.u64(*worst_qoi as u64);
        }
    }
    Ok(encoder.finish())
}
