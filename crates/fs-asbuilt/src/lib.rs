//! fs-asbuilt — reality is just another chart (plan addendum, Proposal 11).
//! Layer: L2.
//!
//! A CT scan (or laser point cloud) of the manufactured part is one more
//! REPRESENTATION with its own restriction maps, so "as-built vs as-designed"
//! becomes a δ between two sections — computed by the same watertightness
//! machinery, closing the loop through the physical world. The "validated"
//! color stops being a static stamp and becomes a living, regime-tagged belief.
//!
//! Two facts make this honest:
//! - REGISTRATION (aligning scan to design) is an OPTIMIZATION WITH ERROR.
//!   [`register`] solves the rigid
//!   2-D fit in closed form (no SVD) and is made WELL-POSED by fiducials/datums
//!   specified at design time (≥ 3 non-collinear points) — the
//!   design-for-verification requirement pushed upstream. Its retained RMS is
//!   only a global fit diagnostic, not transform or spatial uncertainty.
//! - The R8 screen is explicit: if that residual exceeds the geometric
//!   deviation under review, the signal is below this advisory noise screen
//!   ([`well_posed`]); no certification claim follows from passing it.
//!
//! The as-built δ ([`as_built_diff`]) is measurement-noise-aware and emits
//! an **estimated candidate** with a proposed regime. A caller-supplied
//! calibration identity is provenance, not authority: this crate exposes no
//! validated-promotion API until an authenticated verifier and retained
//! calibration artifact are available. Both resource-driving entry points
//! require an [`fs_exec::Cx`], poll at fixed point strides, and publish no
//! partial result when cancellation is observed. Deterministic; pure Rust.

use fs_evidence::color_leaf_identity_reason;
pub use fs_evidence::{Color, ValidityDomain};
use fs_ivl::Interval;

const AS_BUILT_ESTIMATOR_DOMAIN: &str = "org.frankensim.fs-asbuilt.diff-estimator.v4";
const AS_BUILT_ESTIMATOR_SCHEMA: &[u8] = b"fs-asbuilt-diff-estimator-v4";
/// Identity-bound work-plan version for resource-driving scans and hashing.
pub const AS_BUILT_WORK_PLAN_VERSION: u32 = 2;
/// Identity-bound cancellation policy version for resource-driving scans.
pub const AS_BUILT_POLL_POLICY_VERSION: u32 = 2;
/// Maximum point visits between cancellation polls inside a complete scan.
pub const AS_BUILT_POLL_STRIDE_POINTS: usize = 256;
/// Maximum calibration-identity bytes between cancellation polls.
pub const AS_BUILT_POLL_STRIDE_BYTES: usize = 256;
/// Maximum points accepted by registration or one as-built comparison.
pub const MAX_AS_BUILT_POINTS: usize = 1_000_000;

const REGISTER_INITIAL_PHASE: &str = "register.initial";
const REGISTER_DESIGN_CENTROID_PHASE: &str = "register.design-centroid";
const REGISTER_MEASURED_CENTROID_PHASE: &str = "register.measured-centroid";
const REGISTER_SCATTER_PHASE: &str = "register.scatter";
const REGISTER_RESIDUAL_PHASE: &str = "register.residual";
const REGISTER_PUBLISH_PHASE: &str = "register.publish";
const DIFF_INITIAL_PHASE: &str = "as-built-diff.initial";
const DIFF_CALIBRATION_VALIDATION_PHASE: &str = "as-built-diff.calibration-validation";
const DIFF_DEVIATIONS_PHASE: &str = "as-built-diff.deviations";
const DIFF_MAXIMUM_PHASE: &str = "as-built-diff.maximum";
const DIFF_IDENTITY_PHASE: &str = "as-built-diff.identity";
const DIFF_CALIBRATION_HASH_PHASE: &str = "as-built-diff.calibration-hash";
const DIFF_PUBLISH_PHASE: &str = "as-built-diff.publish";

/// A 2-D point (design or measured coordinate).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2 {
    /// x coordinate.
    x: f64,
    /// y coordinate.
    y: f64,
}

impl Point2 {
    /// Construct a finite point.
    ///
    /// # Errors
    /// Refuses NaN or infinite coordinates.
    pub fn new(x: f64, y: f64) -> Result<Point2, RegError> {
        require_finite("point.x", x)?;
        require_finite("point.y", y)?;
        Ok(Point2 {
            x: canonical_zero(x),
            y: canonical_zero(y),
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

    fn dist(self, other: Point2) -> Result<f64, RegError> {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let distance = dx.hypot(dy);
        require_finite("point distance", distance)?;
        Ok(distance)
    }
}

/// A fiducial/datum correspondence: a design reference point and where the scan
/// measured it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fiducial {
    /// The design-time reference location.
    design: Point2,
    /// The location the scan measured for it.
    measured: Point2,
}

impl Fiducial {
    /// A fiducial correspondence.
    #[must_use]
    pub fn new(design: Point2, measured: Point2) -> Fiducial {
        Fiducial { design, measured }
    }

    /// Design-time reference location.
    #[must_use]
    pub const fn design(self) -> Point2 {
        self.design
    }

    /// Measured location.
    #[must_use]
    pub const fn measured(self) -> Point2 {
        self.measured
    }
}

/// A structured registration/ingestion failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RegError {
    /// Cancellation was observed at a deterministic work boundary. No partial
    /// registration or diff is published.
    Cancelled {
        /// Stable operation phase at the observing checkpoint.
        phase: &'static str,
        /// Exact logical work units completed before the checkpoint.
        completed_work: u128,
        /// Exact logical work units planned by the constant-time preflight.
        planned_work: u128,
    },
    /// Fewer fiducials than needed for a well-posed fit.
    TooFewFiducials {
        /// Supplied.
        have: usize,
        /// Required.
        need: usize,
    },
    /// The fiducials are (near-)collinear — the fit is ill-posed.
    CollinearFiducials,
    /// The design points have rank, but their measured correspondence does not
    /// contain any measured spread and therefore cannot determine a rotation.
    UnobservableRotation,
    /// The measured correspondence has spread, but outward arithmetic cannot
    /// certify a nonzero rotation objective at binary64 precision.
    RotationCertificationUnresolved,
    /// Two point sets have mismatched lengths.
    LengthMismatch {
        /// Expected.
        expected: usize,
        /// Found.
        found: usize,
    },
    /// An empty point set.
    Empty,
    /// A resource-driving point set exceeds the public bound.
    TooManyPoints {
        /// Supplied point count.
        have: usize,
        /// Maximum accepted point count.
        max: usize,
    },
    /// A public numeric input is NaN or infinite.
    NonFinite {
        /// Stable field or computation name.
        field: &'static str,
    },
    /// A quantity that must be non-negative was negative.
    Negative {
        /// Stable field name.
        field: &'static str,
    },
    /// A calibration candidate identity is not an admissible provenance leaf.
    InvalidCalibrationIdentity {
        /// Stable structural reason from the shared evidence grammar.
        reason: &'static str,
        /// Input byte length, retained without cloning hostile input.
        bytes: usize,
    },
    /// The bounded deviations vector could not reserve memory.
    AllocationFailed,
    /// A canonical identity field length could not be represented as `u64`.
    IdentityEncodingOverflow,
    /// The complete logical-work plan could not be represented exactly.
    WorkPlanOverflow {
        /// Stable operation name.
        operation: &'static str,
    },
    /// Runtime work accounting did not reconcile with its preflight plan.
    WorkPlanMismatch {
        /// Stable operation phase where reconciliation failed.
        phase: &'static str,
        /// Logical work units observed at reconciliation.
        completed_work: u128,
        /// Logical work units expected at reconciliation.
        planned_work: u128,
    },
}

impl core::fmt::Display for RegError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled {
                phase,
                completed_work,
                planned_work,
            } => write!(
                formatter,
                "as-built operation cancelled during {phase} after {completed_work}/{planned_work} logical work units"
            ),
            Self::TooFewFiducials { have, need } => {
                write!(formatter, "need at least {need} fiducials, got {have}")
            }
            Self::CollinearFiducials => formatter.write_str("fiducials are collinear"),
            Self::UnobservableRotation => {
                formatter.write_str("fiducial correspondence does not determine a rotation")
            }
            Self::RotationCertificationUnresolved => formatter.write_str(
                "fiducial rotation objective is unresolved by the outward-rounded certificate",
            ),
            Self::LengthMismatch { expected, found } => {
                write!(formatter, "expected {expected} scanned points, got {found}")
            }
            Self::Empty => formatter.write_str("point set is empty"),
            Self::TooManyPoints { have, max } => {
                write!(formatter, "point count {have} exceeds bound {max}")
            }
            Self::NonFinite { field } => write!(formatter, "{field} must be finite"),
            Self::Negative { field } => write!(formatter, "{field} must be non-negative"),
            Self::InvalidCalibrationIdentity { reason, bytes } => write!(
                formatter,
                "calibration candidate identity is invalid ({reason}, {bytes} bytes)"
            ),
            Self::AllocationFailed => {
                formatter.write_str("could not reserve the bounded deviations vector")
            }
            Self::IdentityEncodingOverflow => {
                formatter.write_str("canonical identity field length exceeds u64")
            }
            Self::WorkPlanOverflow { operation } => {
                write!(formatter, "{operation} logical-work plan exceeds u128")
            }
            Self::WorkPlanMismatch {
                phase,
                completed_work,
                planned_work,
            } => write!(
                formatter,
                "as-built work accounting mismatch during {phase}: completed {completed_work}, expected {planned_work} logical work units"
            ),
        }
    }
}

impl std::error::Error for RegError {}

fn require_finite(field: &'static str, value: f64) -> Result<(), RegError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(RegError::NonFinite { field })
    }
}

fn require_non_negative(field: &'static str, value: f64) -> Result<(), RegError> {
    require_finite(field, value)?;
    if value >= 0.0 {
        Ok(())
    } else {
        Err(RegError::Negative { field })
    }
}

const fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn affine_component_with_scaled_fallback(
    direct: f64,
    first_coefficient: f64,
    first_value: f64,
    second_coefficient: f64,
    second_value: f64,
    translation: f64,
    field: &'static str,
) -> Result<f64, RegError> {
    if direct.is_finite() {
        return Ok(direct);
    }

    // Preserve the ordinary evaluation bits whenever it succeeds. If an
    // intermediate rotation sum overflowed before a finite translation could
    // cancel it, normalize all three addends by their largest finite source
    // magnitude so cancellation occurs before rescaling. The outward interval
    // is essential: scaling alone can round a truly overflowing exact sum back
    // to `f64::MAX`, so recovery is admitted only when the original real
    // three-term affine sum is certified to lie inside the finite range.
    let scale = first_value
        .abs()
        .max(second_value.abs())
        .max(translation.abs());
    if scale == 0.0 {
        return Ok(0.0);
    }
    let scale_enclosure = Interval::point(scale);
    let normalized_enclosure = Interval::point(first_coefficient)
        * (Interval::point(first_value) / scale_enclosure)
        + Interval::point(second_coefficient) * (Interval::point(second_value) / scale_enclosure)
        + Interval::point(translation) / scale_enclosure;
    let affine_enclosure = normalized_enclosure * scale_enclosure;
    if !(affine_enclosure.lo().is_finite() && affine_enclosure.hi().is_finite()) {
        return Err(RegError::NonFinite { field });
    }
    let normalized = first_coefficient.mul_add(
        first_value / scale,
        second_coefficient.mul_add(second_value / scale, translation / scale),
    );
    let recovered = normalized * scale;
    if !recovered.is_finite() || !affine_enclosure.contains(recovered) {
        return Err(RegError::NonFinite { field });
    }
    Ok(recovered)
}

#[derive(Clone, Copy, Debug, Default)]
struct ScaledSumSquares {
    scale: f64,
    scaled_square_sum: f64,
}

impl ScaledSumSquares {
    fn add(&mut self, value: f64) -> Result<(), RegError> {
        debug_assert!(value.is_finite() && value >= 0.0);
        if value == 0.0 {
            return Ok(());
        }
        if self.scale < value {
            let ratio = self.scale / value;
            self.scaled_square_sum = 1.0 + self.scaled_square_sum * ratio * ratio;
            self.scale = value;
        } else {
            let ratio = value / self.scale;
            self.scaled_square_sum += ratio * ratio;
        }
        require_finite(
            "registration scaled residual sum of squares",
            self.scaled_square_sum,
        )
    }

    fn root_mean_square(self, count: f64) -> Result<f64, RegError> {
        if self.scale == 0.0 {
            return Ok(0.0);
        }
        let residual = self.scale * (self.scaled_square_sum / count).sqrt();
        require_finite("registration residual RMS", residual)?;
        Ok(residual)
    }
}

/// A rigid registration (rotation + translation) mapping design → measured,
/// with the residual it carries forward.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Registration {
    /// Rotation angle (radians).
    rotation_rad: f64,
    /// Translation x.
    tx: f64,
    /// Translation y.
    ty: f64,
    /// Root-mean-square residual of the fit (a global fit diagnostic only).
    residual_rms: f64,
}

impl Registration {
    /// Construct a finite rigid registration with a non-negative residual.
    ///
    /// # Errors
    /// Refuses non-finite transform components or a negative residual.
    pub fn new(rotation_rad: f64, tx: f64, ty: f64, residual_rms: f64) -> Result<Self, RegError> {
        require_finite("registration.rotation_rad", rotation_rad)?;
        require_finite("registration.tx", tx)?;
        require_finite("registration.ty", ty)?;
        require_non_negative("registration.residual_rms", residual_rms)?;
        Ok(Self {
            rotation_rad: canonical_zero(rotation_rad),
            tx: canonical_zero(tx),
            ty: canonical_zero(ty),
            residual_rms: canonical_zero(residual_rms),
        })
    }

    /// Rotation angle in radians.
    #[must_use]
    pub const fn rotation_rad(&self) -> f64 {
        self.rotation_rad
    }

    /// x translation.
    #[must_use]
    pub const fn tx(&self) -> f64 {
        self.tx
    }

    /// y translation.
    #[must_use]
    pub const fn ty(&self) -> f64 {
        self.ty
    }

    /// Registration residual RMS.
    #[must_use]
    pub const fn residual_rms(&self) -> f64 {
        self.residual_rms
    }

    /// Map a design point into measured coordinates.
    ///
    /// # Errors
    /// Refuses arithmetic overflow to a non-finite mapped point.
    pub fn apply(&self, point: Point2) -> Result<Point2, RegError> {
        let (s, c) = self.rotation_rad.sin_cos();
        let direct_x = c * point.x - s * point.y + self.tx;
        let direct_y = s * point.x + c * point.y + self.ty;
        Point2::new(
            affine_component_with_scaled_fallback(
                direct_x, c, point.x, -s, point.y, self.tx, "point.x",
            )?,
            affine_component_with_scaled_fallback(
                direct_y, s, point.x, c, point.y, self.ty, "point.y",
            )?,
        )
    }
}

/// The minimum fiducials for a well-posed 2-D rigid fit.
pub const MIN_FIDUCIALS: usize = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RegistrationWorkPlan {
    points_per_scan: u128,
    total: u128,
}

impl RegistrationWorkPlan {
    fn preflight(point_count: usize) -> Result<Self, RegError> {
        if point_count < MIN_FIDUCIALS {
            return Err(RegError::TooFewFiducials {
                have: point_count,
                need: MIN_FIDUCIALS,
            });
        }
        if point_count > MAX_AS_BUILT_POINTS {
            return Err(RegError::TooManyPoints {
                have: point_count,
                max: MAX_AS_BUILT_POINTS,
            });
        }
        let points_per_scan =
            u128::try_from(point_count).map_err(|_| RegError::WorkPlanOverflow {
                operation: "register",
            })?;
        let total = points_per_scan
            .checked_mul(6)
            .ok_or(RegError::WorkPlanOverflow {
                operation: "register",
            })?;
        Ok(Self {
            points_per_scan,
            total,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DiffWorkPlan {
    points_per_scan: u128,
    calibration_validation_byte_units: u128,
    calibration_hash_byte_units: u128,
    total: u128,
}

impl DiffWorkPlan {
    fn preflight(
        design_len: usize,
        scanned_len: usize,
        design_tolerance: f64,
        measurement_noise: f64,
        calibration_candidate: &str,
    ) -> Result<Self, RegError> {
        if design_len == 0 {
            return Err(RegError::Empty);
        }
        if design_len != scanned_len {
            return Err(RegError::LengthMismatch {
                expected: design_len,
                found: scanned_len,
            });
        }
        if design_len > MAX_AS_BUILT_POINTS {
            return Err(RegError::TooManyPoints {
                have: design_len,
                max: MAX_AS_BUILT_POINTS,
            });
        }
        require_non_negative("design_tolerance", design_tolerance)?;
        require_non_negative("measurement_noise", measurement_noise)?;
        if calibration_candidate.len() > fs_evidence::MAX_COLOR_IDENTITY_BYTES {
            return Err(RegError::InvalidCalibrationIdentity {
                reason: "too-long",
                bytes: calibration_candidate.len(),
            });
        }
        let points_per_scan =
            u128::try_from(design_len).map_err(|_| RegError::WorkPlanOverflow {
                operation: "as-built-diff",
            })?;
        let calibration_validation_byte_units = u128::try_from(calibration_candidate.len())
            .map_err(|_| RegError::WorkPlanOverflow {
                operation: "as-built-diff",
            })?;
        let calibration_hash_byte_units = calibration_validation_byte_units;
        let point_work = points_per_scan
            .checked_mul(3)
            .ok_or(RegError::WorkPlanOverflow {
                operation: "as-built-diff",
            })?;
        let byte_work = calibration_validation_byte_units
            .checked_add(calibration_hash_byte_units)
            .ok_or(RegError::WorkPlanOverflow {
                operation: "as-built-diff",
            })?;
        let total = point_work
            .checked_add(byte_work)
            .ok_or(RegError::WorkPlanOverflow {
                operation: "as-built-diff",
            })?;
        Ok(Self {
            points_per_scan,
            calibration_validation_byte_units,
            calibration_hash_byte_units,
            total,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct WorkProgress {
    completed: u128,
    planned: u128,
    operation: &'static str,
}

impl WorkProgress {
    const fn new(planned: u128, operation: &'static str) -> Self {
        Self {
            completed: 0,
            planned,
            operation,
        }
    }

    fn advance(&mut self, units: u128) -> Result<(), RegError> {
        let completed = self
            .completed
            .checked_add(units)
            .ok_or(RegError::WorkPlanOverflow {
                operation: self.operation,
            })?;
        if completed > self.planned {
            return Err(RegError::WorkPlanMismatch {
                phase: self.operation,
                completed_work: completed,
                planned_work: self.planned,
            });
        }
        self.completed = completed;
        Ok(())
    }

    fn complete_point(&mut self) -> Result<(), RegError> {
        self.advance(1)
    }

    fn require_completed(self, phase: &'static str, planned_work: u128) -> Result<(), RegError> {
        if self.completed == planned_work {
            Ok(())
        } else {
            Err(RegError::WorkPlanMismatch {
                phase,
                completed_work: self.completed,
                planned_work,
            })
        }
    }
}

fn operation_checkpoint(
    phase: &'static str,
    progress: WorkProgress,
    poll: &mut impl FnMut(&'static str, u128, u128) -> Result<(), fs_exec::Cancelled>,
) -> Result<(), RegError> {
    poll(phase, progress.completed, progress.planned).map_err(|_| RegError::Cancelled {
        phase,
        completed_work: progress.completed,
        planned_work: progress.planned,
    })
}

fn scan_checkpoint(
    index: usize,
    stride_points: usize,
    phase: &'static str,
    progress: WorkProgress,
    poll: &mut impl FnMut(&'static str, u128, u128) -> Result<(), fs_exec::Cancelled>,
) -> Result<(), RegError> {
    debug_assert!(stride_points != 0);
    if index != 0 && index.is_multiple_of(stride_points) {
        operation_checkpoint(phase, progress, poll)?;
    }
    Ok(())
}

/// Solve the rigid 2-D registration that best maps `fiducials`' design points
/// onto their measured points (closed-form least squares — the 2-D Umeyama /
/// Procrustes rotation). Requires ≥ 3 non-collinear fiducials.
///
/// # Errors
/// [`RegError::TooFewFiducials`], [`RegError::CollinearFiducials`],
/// [`RegError::UnobservableRotation`],
/// [`RegError::RotationCertificationUnresolved`], oversized or non-finite
/// inputs/intermediates, work-plan overflow, or a structured
/// [`RegError::Cancelled`] with exact point-visit progress. The complete work
/// plan is computed before the initial cancellation checkpoint.
pub fn register(fiducials: &[Fiducial], cx: &fs_exec::Cx<'_>) -> Result<Registration, RegError> {
    let mut poll = |_: &'static str, _: u128, _: u128| cx.checkpoint();
    register_with_poll(fiducials, &mut poll)
}

fn register_with_poll(
    fiducials: &[Fiducial],
    poll: &mut impl FnMut(&'static str, u128, u128) -> Result<(), fs_exec::Cancelled>,
) -> Result<Registration, RegError> {
    let plan = RegistrationWorkPlan::preflight(fiducials.len())?;
    let mut progress = WorkProgress::new(plan.total, "register");
    operation_checkpoint(REGISTER_INITIAL_PHASE, progress, poll)?;

    let n = fiducials.len();
    let nf = f64::from(u32::try_from(n).map_err(|_| RegError::TooManyPoints {
        have: n,
        max: MAX_AS_BUILT_POINTS,
    })?);
    operation_checkpoint(REGISTER_DESIGN_CENTROID_PHASE, progress, poll)?;
    let cp = centroid(
        fiducials,
        |fiducial| fiducial.design,
        REGISTER_DESIGN_CENTROID_PHASE,
        &mut progress,
        poll,
    )?;
    progress.require_completed(REGISTER_DESIGN_CENTROID_PHASE, plan.points_per_scan * 2)?;
    operation_checkpoint(REGISTER_MEASURED_CENTROID_PHASE, progress, poll)?;
    let cq = centroid(
        fiducials,
        |fiducial| fiducial.measured,
        REGISTER_MEASURED_CENTROID_PHASE,
        &mut progress,
        poll,
    )?;
    progress.require_completed(REGISTER_MEASURED_CENTROID_PHASE, plan.points_per_scan * 4)?;
    if cp.scale <= 0.0 {
        return Err(RegError::CollinearFiducials);
    }
    if cq.scale <= 0.0 {
        return Err(RegError::UnobservableRotation);
    }

    // scatter of the centered DESIGN points — collinear iff it is rank-deficient.
    let (mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0);
    // cross-covariance terms for the optimal rotation.
    let (mut s_dot, mut s_cross) = (0.0, 0.0);
    let mut s_dot_enclosure = Interval::point(0.0);
    let mut s_cross_enclosure = Interval::point(0.0);
    operation_checkpoint(REGISTER_SCATTER_PHASE, progress, poll)?;
    for (index, f) in fiducials.iter().enumerate() {
        scan_checkpoint(
            index,
            AS_BUILT_POLL_STRIDE_POINTS,
            REGISTER_SCATTER_PHASE,
            progress,
            poll,
        )?;
        // Positive common scales cancel from the rank determinant and the
        // Procrustes rotation objective. Normalize before products so finite
        // tiny or huge geometries do not underflow/overflow fourth-order rank
        // expressions or second-order cross sums.
        let dpx = normalized_offset(f.design.x, cp.anchor.x, cp.scale, "design x")?
            - cp.normalized_mean.x;
        let dpy = normalized_offset(f.design.y, cp.anchor.y, cp.scale, "design y")?
            - cp.normalized_mean.y;
        let dqx = normalized_offset(f.measured.x, cq.anchor.x, cq.scale, "measured x")?
            - cq.normalized_mean.x;
        let dqy = normalized_offset(f.measured.y, cq.anchor.y, cq.scale, "measured y")?
            - cq.normalized_mean.y;
        let dpx_enclosure =
            normalized_offset_enclosure(f.design.x, cp.anchor.x, cp.scale) - cp.normalized_x;
        let dpy_enclosure =
            normalized_offset_enclosure(f.design.y, cp.anchor.y, cp.scale) - cp.normalized_y;
        let dqx_enclosure =
            normalized_offset_enclosure(f.measured.x, cq.anchor.x, cq.scale) - cq.normalized_x;
        let dqy_enclosure =
            normalized_offset_enclosure(f.measured.y, cq.anchor.y, cq.scale) - cq.normalized_y;
        sxx += dpx * dpx;
        syy += dpy * dpy;
        sxy += dpx * dpy;
        s_dot += dpx * dqx + dpy * dqy;
        s_cross += dpx * dqy - dpy * dqx;
        s_dot_enclosure =
            s_dot_enclosure + dpx_enclosure * dqx_enclosure + dpy_enclosure * dqy_enclosure;
        s_cross_enclosure =
            s_cross_enclosure + dpx_enclosure * dqy_enclosure - dpy_enclosure * dqx_enclosure;
        progress.complete_point()?;
    }
    progress.require_completed(REGISTER_SCATTER_PHASE, plan.points_per_scan * 5)?;
    for (field, value) in [
        ("registration design scatter xx", sxx),
        ("registration design scatter yy", syy),
        ("registration design scatter xy", sxy),
        ("registration cross-covariance dot", s_dot),
        ("registration cross-covariance cross", s_cross),
    ] {
        require_finite(field, value)?;
    }
    let det = sxx * syy - sxy * sxy;
    let trace = sxx + syy;
    require_finite("registration scatter determinant", det)?;
    require_finite("registration scatter trace", trace)?;
    let trace_squared = trace * trace;
    require_finite("registration squared scatter trace", trace_squared)?;
    if trace <= 0.0 || det <= 1e-12 * trace_squared {
        return Err(RegError::CollinearFiducials);
    }

    // atan2(0, 0) returns zero, but that is a library convention rather than
    // an inferred rotation. The measured-extent gate above catches a fully
    // collapsed scan. The interval sums then require at least one exact
    // cross-covariance component to be separated from zero; if both can be
    // zero, the orientation objective is unproved and admission fails closed
    // without an arbitrary epsilon.
    let finite_observability_enclosures = [s_dot_enclosure, s_cross_enclosure]
        .into_iter()
        .all(|enclosure| enclosure.lo().is_finite() && enclosure.hi().is_finite());
    if !finite_observability_enclosures
        || (s_dot_enclosure.contains_zero() && s_cross_enclosure.contains_zero())
    {
        return Err(RegError::RotationCertificationUnresolved);
    }

    let rotation_rad = s_cross.atan2(s_dot);
    let (s, c) = rotation_rad.sin_cos();
    let direct_tx = cq.point.x - (c * cp.point.x - s * cp.point.y);
    let direct_ty = cq.point.y - (s * cp.point.x + c * cp.point.y);
    let tx = affine_component_with_scaled_fallback(
        direct_tx,
        -c,
        cp.point.x,
        s,
        cp.point.y,
        cq.point.x,
        "registration.tx",
    )?;
    let ty = affine_component_with_scaled_fallback(
        direct_ty,
        -s,
        cp.point.x,
        -c,
        cp.point.y,
        cq.point.y,
        "registration.ty",
    )?;
    let reg = Registration::new(rotation_rad, tx, ty, 0.0)?;
    // Residual RMS is retained as a global fit diagnostic. It is not transform
    // covariance or a pointwise spatial uncertainty bound.
    let mut residuals = ScaledSumSquares::default();
    operation_checkpoint(REGISTER_RESIDUAL_PHASE, progress, poll)?;
    for (index, fiducial) in fiducials.iter().enumerate() {
        scan_checkpoint(
            index,
            AS_BUILT_POLL_STRIDE_POINTS,
            REGISTER_RESIDUAL_PHASE,
            progress,
            poll,
        )?;
        let distance = reg.apply(fiducial.design)?.dist(fiducial.measured)?;
        residuals.add(distance)?;
        progress.complete_point()?;
    }
    progress.require_completed(REGISTER_RESIDUAL_PHASE, plan.total)?;
    operation_checkpoint(REGISTER_PUBLISH_PHASE, progress, poll)?;
    Registration::new(rotation_rad, tx, ty, residuals.root_mean_square(nf)?)
}

#[derive(Clone, Copy)]
struct CertifiedCentroid {
    point: Point2,
    anchor: Point2,
    normalized_mean: Point2,
    normalized_x: Interval,
    normalized_y: Interval,
    scale: f64,
}

fn stable_running_mean(current: f64, value: f64, count: u32) -> f64 {
    if count == 1 {
        return value;
    }
    let divisor = f64::from(count);
    let difference = value - current;
    if difference.is_finite() {
        current + difference / divisor
    } else {
        // A finite subtraction can overflow only for opposite-sign extreme
        // operands. This convex form keeps both products finite and their sum
        // cannot overflow because the terms have opposite signs.
        let next_weight = 1.0 / divisor;
        current.mul_add(1.0 - next_weight, value * next_weight)
    }
}

fn normalized_offset(
    value: f64,
    anchor: f64,
    scale: f64,
    field: &'static str,
) -> Result<f64, RegError> {
    let direct = (value - anchor) / scale;
    let normalized = if direct.is_finite() {
        direct
    } else {
        // The midpoint anchor normally makes the direct difference finite.
        // Retain a fail-safe scaled subtraction for the extreme opposite-sign
        // endpoint where the raw subtraction rounds to infinity.
        value / scale - anchor / scale
    };
    require_finite(field, normalized)?;
    Ok(normalized)
}

fn normalized_offset_enclosure(value: f64, anchor: f64, scale: f64) -> Interval {
    let scale = Interval::point(scale);
    let direct = (Interval::point(value) - Interval::point(anchor)) / scale;
    if direct.lo().is_finite() && direct.hi().is_finite() {
        direct
    } else {
        Interval::point(value) / scale - Interval::point(anchor) / scale
    }
}

fn centroid(
    fiducials: &[Fiducial],
    select: impl Fn(&Fiducial) -> Point2 + Copy,
    phase: &'static str,
    progress: &mut WorkProgress,
    poll: &mut impl FnMut(&'static str, u128, u128) -> Result<(), fs_exec::Cancelled>,
) -> Result<CertifiedCentroid, RegError> {
    let mut count = 0_u32;
    let (mut mean_x, mut mean_y) = (0.0, 0.0);
    let (mut min_x, mut max_x) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut min_y, mut max_y) = (f64::INFINITY, f64::NEG_INFINITY);
    for (index, fiducial) in fiducials.iter().enumerate() {
        scan_checkpoint(index, AS_BUILT_POLL_STRIDE_POINTS, phase, *progress, poll)?;
        let point = select(fiducial);
        count = count.checked_add(1).ok_or(RegError::WorkPlanOverflow {
            operation: "register",
        })?;
        mean_x = stable_running_mean(mean_x, point.x, count);
        mean_y = stable_running_mean(mean_y, point.y, count);
        min_x = min_x.min(point.x);
        max_x = max_x.max(point.x);
        min_y = min_y.min(point.y);
        max_y = max_y.max(point.y);
        require_finite("point centroid x", mean_x)?;
        require_finite("point centroid y", mean_y)?;
        progress.complete_point()?;
    }
    let anchor = Point2::new(f64::midpoint(min_x, max_x), f64::midpoint(min_y, max_y))?;
    let scale = (max_x - anchor.x)
        .abs()
        .max((min_x - anchor.x).abs())
        .max((max_y - anchor.y).abs())
        .max((min_y - anchor.y).abs());
    require_finite("point-set normalization scale", scale)?;
    operation_checkpoint(phase, *progress, poll)?;

    let mut normalized_sum_x = 0.0;
    let mut normalized_sum_y = 0.0;
    let mut normalized_sum_x_enclosure = Interval::point(0.0);
    let mut normalized_sum_y_enclosure = Interval::point(0.0);
    for (index, fiducial) in fiducials.iter().enumerate() {
        scan_checkpoint(index, AS_BUILT_POLL_STRIDE_POINTS, phase, *progress, poll)?;
        let point = select(fiducial);
        if scale > 0.0 {
            let normalized_x = normalized_offset(point.x, anchor.x, scale, "normalized point x")?;
            let normalized_y = normalized_offset(point.y, anchor.y, scale, "normalized point y")?;
            normalized_sum_x += normalized_x;
            normalized_sum_y += normalized_y;
            normalized_sum_x_enclosure =
                normalized_sum_x_enclosure + normalized_offset_enclosure(point.x, anchor.x, scale);
            normalized_sum_y_enclosure =
                normalized_sum_y_enclosure + normalized_offset_enclosure(point.y, anchor.y, scale);
            require_finite("normalized centroid x sum", normalized_sum_x)?;
            require_finite("normalized centroid y sum", normalized_sum_y)?;
        }
        progress.complete_point()?;
    }
    let divisor = f64::from(count);
    let normalized_mean = if scale > 0.0 {
        Point2::new(normalized_sum_x / divisor, normalized_sum_y / divisor)?
    } else {
        Point2::new(0.0, 0.0)?
    };
    let interval_divisor = Interval::point(divisor);
    let reconstruct = |anchor: f64, normalized: f64, fallback: f64, lo: f64, hi: f64| {
        let anchored = scale.mul_add(normalized, anchor);
        if anchored.is_finite() {
            anchored.clamp(lo, hi)
        } else {
            // The stable online mean is an independently finite convex estimate
            // for the rare case where rounded reconstruction crosses MAX.
            fallback
        }
    };
    let point = Point2::new(
        reconstruct(anchor.x, normalized_mean.x, mean_x, min_x, max_x),
        reconstruct(anchor.y, normalized_mean.y, mean_y, min_y, max_y),
    )?;
    Ok(CertifiedCentroid {
        point,
        anchor,
        normalized_mean,
        normalized_x: normalized_sum_x_enclosure / interval_divisor,
        normalized_y: normalized_sum_y_enclosure / interval_divisor,
        scale,
    })
}

/// The R8 residual-proxy screen: the global registration fit residual is below
/// the supplied geometric-deviation signal. This is an advisory screen, not a
/// proof that registration is trustworthy: `residual_rms` is neither a
/// pointwise uncertainty bound nor a calibrated confidence statement. If the
/// residual meets or exceeds the signal, the as-built loop is premature for
/// that part class (defer to point-sensor assimilation).
#[must_use]
pub fn well_posed(reg: &Registration, certified_deviation: f64) -> bool {
    certified_deviation.is_finite()
        && certified_deviation > 0.0
        && reg.residual_rms < certified_deviation
}

/// The as-built δ between design and scanned sections.
#[derive(Debug, Clone, PartialEq)]
pub struct AsBuiltDiff {
    /// Per-point deviation `||registered(design) − scanned||`.
    deviations: Vec<f64>,
    /// The largest deviation.
    max_deviation: f64,
    /// Advisory one-dispersion screen for the design tolerance.
    within_tolerance: bool,
    /// Advisory one-dispersion screen for whether the maximum deviation rises
    /// above the conservatively combined estimated dispersion.
    above_noise_floor: bool,
    /// Proposed regime for later calibration-authority review.
    proposed_regime: ValidityDomain,
    /// The δ's honest candidate color. This API never emits `Validated`.
    color: Color,
}

impl AsBuiltDiff {
    /// Per-point deviations in the input order.
    #[must_use]
    pub fn deviations(&self) -> &[f64] {
        &self.deviations
    }

    /// Largest point deviation.
    #[must_use]
    pub const fn max_deviation(&self) -> f64 {
        self.max_deviation
    }

    /// Whether the maximum deviation plus one conservatively combined
    /// estimated dispersion fits the supplied design tolerance. This advisory
    /// screen is not a tolerance certificate.
    #[must_use]
    pub const fn within_tolerance(&self) -> bool {
        self.within_tolerance
    }

    /// Whether the largest deviation exceeds one conservatively combined
    /// estimated dispersion. This advisory screen is not a statistical
    /// significance test.
    #[must_use]
    pub const fn above_noise_floor(&self) -> bool {
        self.above_noise_floor
    }

    /// Proposed, unauthenticated validity regime for later review.
    #[must_use]
    pub const fn proposed_regime(&self) -> &ValidityDomain {
        &self.proposed_regime
    }

    /// Honest candidate color produced by [`as_built_diff`].
    #[must_use]
    pub const fn color(&self) -> &Color {
        &self.color
    }
}

/// Compute the as-built δ after registration: apply the registration to each
/// design point and measure its deviation from the corresponding scanned point.
/// The δ is colored ESTIMATED. Its bounded, domain-separated identity binds
/// every point, registration component, tolerance, noise value, and the
/// structurally valid calibration candidate identity, plus the execution mode,
/// every budget field, the exact checked work plan, and the versioned poll
/// policy. The proposed regime is carried separately for later authenticated
/// calibration review. The returned decision booleans are conservative,
/// advisory one-dispersion screens: registration residual RMS is a global fit
/// diagnostic rather than a pointwise uncertainty bound, so neither boolean is
/// a tolerance certificate or statistical-significance claim.
///
/// # Errors
/// Refuses empty/mismatched/oversized point sets, malformed calibration
/// identities, negative or non-finite tolerances/noise, and non-finite
/// arithmetic results. Cancellation returns [`RegError::Cancelled`] with exact
/// progress and never publishes a partial diff.
pub fn as_built_diff(
    reg: &Registration,
    design: &[Point2],
    scanned: &[Point2],
    design_tolerance: f64,
    measurement_noise: f64,
    calibration_candidate: &str,
    cx: &fs_exec::Cx<'_>,
) -> Result<AsBuiltDiff, RegError> {
    let execution = ExecutionIdentity::from_cx(cx);
    let mut poll = |_: &'static str, _: u128, _: u128| cx.checkpoint();
    as_built_diff_with_poll(
        reg,
        design,
        scanned,
        design_tolerance,
        measurement_noise,
        calibration_candidate,
        execution,
        CURRENT_POLL_POLICY,
        &mut poll,
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ExecutionIdentity {
    mode: fs_exec::ExecMode,
    budget: fs_exec::Budget,
}

impl ExecutionIdentity {
    fn from_cx(cx: &fs_exec::Cx<'_>) -> Self {
        Self {
            mode: cx.mode(),
            budget: cx.budget(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PollPolicy {
    version: u32,
    stride_points: usize,
    stride_bytes: usize,
}

const CURRENT_POLL_POLICY: PollPolicy = PollPolicy {
    version: AS_BUILT_POLL_POLICY_VERSION,
    stride_points: AS_BUILT_POLL_STRIDE_POINTS,
    stride_bytes: AS_BUILT_POLL_STRIDE_BYTES,
};

#[allow(clippy::too_many_arguments)]
fn as_built_diff_with_poll(
    reg: &Registration,
    design: &[Point2],
    scanned: &[Point2],
    design_tolerance: f64,
    measurement_noise: f64,
    calibration_candidate: &str,
    execution: ExecutionIdentity,
    poll_policy: PollPolicy,
    poll: &mut impl FnMut(&'static str, u128, u128) -> Result<(), fs_exec::Cancelled>,
) -> Result<AsBuiltDiff, RegError> {
    let plan = DiffWorkPlan::preflight(
        design.len(),
        scanned.len(),
        design_tolerance,
        measurement_noise,
        calibration_candidate,
    )?;
    let mut progress = WorkProgress::new(plan.total, "as-built-diff");
    operation_checkpoint(DIFF_INITIAL_PHASE, progress, poll)?;

    operation_checkpoint(DIFF_CALIBRATION_VALIDATION_PHASE, progress, poll)?;
    let invalid_calibration_reason = color_leaf_identity_reason(calibration_candidate);
    progress.advance(plan.calibration_validation_byte_units)?;
    operation_checkpoint(DIFF_CALIBRATION_VALIDATION_PHASE, progress, poll)?;
    progress.require_completed(
        DIFF_CALIBRATION_VALIDATION_PHASE,
        plan.calibration_validation_byte_units,
    )?;
    if let Some(reason) = invalid_calibration_reason {
        return Err(RegError::InvalidCalibrationIdentity {
            reason,
            bytes: calibration_candidate.len(),
        });
    }

    let mut deviations = Vec::new();
    deviations
        .try_reserve_exact(design.len())
        .map_err(|_| RegError::AllocationFailed)?;
    operation_checkpoint(DIFF_DEVIATIONS_PHASE, progress, poll)?;
    for (index, (design_point, scanned_point)) in design.iter().zip(scanned).enumerate() {
        scan_checkpoint(
            index,
            poll_policy.stride_points,
            DIFF_DEVIATIONS_PHASE,
            progress,
            poll,
        )?;
        deviations.push(reg.apply(*design_point)?.dist(*scanned_point)?);
        progress.complete_point()?;
    }
    progress.require_completed(
        DIFF_DEVIATIONS_PHASE,
        plan.calibration_validation_byte_units + plan.points_per_scan,
    )?;

    operation_checkpoint(DIFF_MAXIMUM_PHASE, progress, poll)?;
    let mut max_deviation = 0.0_f64;
    for (index, deviation) in deviations.iter().copied().enumerate() {
        scan_checkpoint(
            index,
            poll_policy.stride_points,
            DIFF_MAXIMUM_PHASE,
            progress,
            poll,
        )?;
        max_deviation = max_deviation.max(deviation);
        progress.complete_point()?;
    }
    progress.require_completed(
        DIFF_MAXIMUM_PHASE,
        plan.calibration_validation_byte_units + plan.points_per_scan * 2,
    )?;
    require_finite("maximum as-built deviation", max_deviation)?;
    let proposed_regime = ValidityDomain::unconstrained()
        .with(
            "registration_residual",
            0.0,
            reg.residual_rms.max(f64::MIN_POSITIVE),
        )
        .with(
            "measurement_noise",
            0.0,
            measurement_noise.max(f64::MIN_POSITIVE),
        )
        .with(
            "design_tolerance",
            0.0,
            design_tolerance.max(f64::MIN_POSITIVE),
        );
    // `Estimated` dispersions compose additively unless a calibrated
    // independence model establishes a sharper rule. The registration RMS is
    // only a global fit diagnostic, not a pointwise uncertainty bound.
    let dispersion = reg.residual_rms + measurement_noise;
    require_finite("combined as-built dispersion", dispersion)?;
    let estimator = estimator_identity(
        reg,
        design,
        scanned,
        design_tolerance,
        measurement_noise,
        calibration_candidate,
        execution,
        plan,
        poll_policy,
        &mut progress,
        poll,
    )?;
    let within_tolerance =
        max_deviation <= design_tolerance && dispersion <= design_tolerance - max_deviation;
    let above_noise_floor = max_deviation > dispersion;
    progress.require_completed(DIFF_IDENTITY_PHASE, plan.total)?;
    let output = AsBuiltDiff {
        deviations,
        max_deviation,
        within_tolerance,
        above_noise_floor,
        proposed_regime,
        color: Color::Estimated {
            estimator,
            dispersion,
        },
    };
    operation_checkpoint(DIFF_PUBLISH_PHASE, progress, poll)?;
    Ok(output)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn estimator_identity(
    registration: &Registration,
    design: &[Point2],
    scanned: &[Point2],
    design_tolerance: f64,
    measurement_noise: f64,
    calibration_candidate: &str,
    execution: ExecutionIdentity,
    work_plan: DiffWorkPlan,
    poll_policy: PollPolicy,
    progress: &mut WorkProgress,
    poll: &mut impl FnMut(&'static str, u128, u128) -> Result<(), fs_exec::Cancelled>,
) -> Result<String, RegError> {
    fn field(hasher: &mut fs_blake3::Blake3, bytes: &[u8]) -> Result<(), RegError> {
        let length = u64::try_from(bytes.len()).map_err(|_| RegError::IdentityEncodingOverflow)?;
        hasher.update(&length.to_le_bytes());
        hasher.update(bytes);
        Ok(())
    }

    fn number(hasher: &mut fs_blake3::Blake3, label: &[u8], value: f64) -> Result<(), RegError> {
        field(hasher, label)?;
        field(hasher, &canonical_zero(value).to_bits().to_le_bytes())
    }

    fn unsigned(
        hasher: &mut fs_blake3::Blake3,
        label: &[u8],
        bytes: &[u8],
    ) -> Result<(), RegError> {
        field(hasher, label)?;
        field(hasher, bytes)
    }

    operation_checkpoint(DIFF_IDENTITY_PHASE, *progress, poll)?;
    let mut hasher = fs_blake3::Blake3::new();
    field(&mut hasher, AS_BUILT_ESTIMATOR_SCHEMA)?;
    field(&mut hasher, b"execution.mode")?;
    field(&mut hasher, execution.mode.name().as_bytes())?;
    field(&mut hasher, b"execution.budget.deadline-present")?;
    field(
        &mut hasher,
        &[u8::from(execution.budget.deadline.is_some())],
    )?;
    if let Some(deadline) = execution.budget.deadline {
        unsigned(
            &mut hasher,
            b"execution.budget.deadline-nanos",
            &deadline.as_nanos().to_le_bytes(),
        )?;
    }
    unsigned(
        &mut hasher,
        b"execution.budget.poll-quota",
        &execution.budget.poll_quota.to_le_bytes(),
    )?;
    field(&mut hasher, b"execution.budget.cost-quota-present")?;
    field(
        &mut hasher,
        &[u8::from(execution.budget.cost_quota.is_some())],
    )?;
    if let Some(cost_quota) = execution.budget.cost_quota {
        unsigned(
            &mut hasher,
            b"execution.budget.cost-quota",
            &cost_quota.to_le_bytes(),
        )?;
    }
    unsigned(
        &mut hasher,
        b"execution.budget.priority",
        &[execution.budget.priority],
    )?;
    unsigned(
        &mut hasher,
        b"work-plan.version",
        &AS_BUILT_WORK_PLAN_VERSION.to_le_bytes(),
    )?;
    unsigned(
        &mut hasher,
        b"work-plan.deviation-point-visits",
        &work_plan.points_per_scan.to_le_bytes(),
    )?;
    unsigned(
        &mut hasher,
        b"work-plan.maximum-point-visits",
        &work_plan.points_per_scan.to_le_bytes(),
    )?;
    unsigned(
        &mut hasher,
        b"work-plan.identity-point-pair-visits",
        &work_plan.points_per_scan.to_le_bytes(),
    )?;
    unsigned(
        &mut hasher,
        b"work-plan.calibration-validation-byte-units",
        &work_plan.calibration_validation_byte_units.to_le_bytes(),
    )?;
    unsigned(
        &mut hasher,
        b"work-plan.calibration-hash-byte-units",
        &work_plan.calibration_hash_byte_units.to_le_bytes(),
    )?;
    unsigned(
        &mut hasher,
        b"work-plan.total-work-units",
        &work_plan.total.to_le_bytes(),
    )?;
    unsigned(
        &mut hasher,
        b"poll-policy.version",
        &poll_policy.version.to_le_bytes(),
    )?;
    let stride_points =
        u128::try_from(poll_policy.stride_points).map_err(|_| RegError::WorkPlanOverflow {
            operation: "as-built-diff",
        })?;
    unsigned(
        &mut hasher,
        b"poll-policy.stride-points",
        &stride_points.to_le_bytes(),
    )?;
    let stride_bytes =
        u128::try_from(poll_policy.stride_bytes).map_err(|_| RegError::WorkPlanOverflow {
            operation: "as-built-diff",
        })?;
    unsigned(
        &mut hasher,
        b"poll-policy.stride-bytes",
        &stride_bytes.to_le_bytes(),
    )?;
    field(&mut hasher, b"calibration-candidate")?;
    operation_checkpoint(DIFF_CALIBRATION_HASH_PHASE, *progress, poll)?;
    field(&mut hasher, calibration_candidate.as_bytes())?;
    progress.advance(work_plan.calibration_hash_byte_units)?;
    operation_checkpoint(DIFF_CALIBRATION_HASH_PHASE, *progress, poll)?;
    progress.require_completed(
        DIFF_CALIBRATION_HASH_PHASE,
        work_plan.calibration_validation_byte_units
            + work_plan.points_per_scan * 2
            + work_plan.calibration_hash_byte_units,
    )?;
    number(
        &mut hasher,
        b"registration.rotation_rad",
        registration.rotation_rad,
    )?;
    number(&mut hasher, b"registration.tx", registration.tx)?;
    number(&mut hasher, b"registration.ty", registration.ty)?;
    number(
        &mut hasher,
        b"registration.residual_rms",
        registration.residual_rms,
    )?;
    number(&mut hasher, b"design_tolerance", design_tolerance)?;
    number(&mut hasher, b"measurement_noise", measurement_noise)?;
    field(&mut hasher, b"point-count")?;
    let point_count =
        u64::try_from(design.len()).map_err(|_| RegError::IdentityEncodingOverflow)?;
    field(&mut hasher, &point_count.to_le_bytes())?;
    for (ordinal, (design_point, scanned_point)) in design.iter().zip(scanned).enumerate() {
        scan_checkpoint(
            ordinal,
            poll_policy.stride_points,
            DIFF_IDENTITY_PHASE,
            *progress,
            poll,
        )?;
        field(&mut hasher, b"point-pair")?;
        let ordinal = u64::try_from(ordinal).map_err(|_| RegError::IdentityEncodingOverflow)?;
        field(&mut hasher, &ordinal.to_le_bytes())?;
        number(&mut hasher, b"design.x", design_point.x)?;
        number(&mut hasher, b"design.y", design_point.y)?;
        number(&mut hasher, b"scanned.x", scanned_point.x)?;
        number(&mut hasher, b"scanned.y", scanned_point.y)?;
        progress.complete_point()?;
    }
    progress.require_completed(DIFF_IDENTITY_PHASE, work_plan.total)?;
    let preimage_hash = hasher.finalize();
    Ok(format!(
        "asbuilt-diff-v4:{}",
        fs_blake3::hash_domain(AS_BUILT_ESTIMATOR_DOMAIN, preimage_hash.as_bytes())
    ))
}

#[cfg(test)]
mod cancellation_tests {
    use super::*;

    #[test]
    fn g0_scaled_residual_rms_avoids_intermediate_square_overflow() {
        let mut residuals = ScaledSumSquares::default();
        residuals.add(f64::MAX).unwrap();
        residuals.add(f64::MAX).unwrap();
        let rms = residuals.root_mean_square(2.0).unwrap();
        assert_eq!(rms.to_bits(), f64::MAX.to_bits());
    }

    fn points(count: usize) -> Vec<Point2> {
        (0..count)
            .map(|index| {
                let coordinate = f64::from(u32::try_from(index).expect("small test index"));
                Point2::new(coordinate, coordinate.mul_add(0.5, 1.0)).expect("finite test point")
            })
            .collect()
    }

    fn fiducials(count: usize) -> Vec<Fiducial> {
        points(count)
            .into_iter()
            .map(|point| Fiducial::new(point, point))
            .collect()
    }

    fn execution() -> ExecutionIdentity {
        ExecutionIdentity {
            mode: fs_exec::ExecMode::Deterministic,
            budget: fs_exec::Budget::INFINITE,
        }
    }

    fn identity(diff: &AsBuiltDiff) -> &str {
        match diff.color() {
            Color::Estimated { estimator, .. } => estimator,
            other => panic!("expected estimated diff, got {other:?}"),
        }
    }

    #[test]
    fn g4_stride_boundary_and_plus_one_have_exact_phase_progress() {
        let boundary = fiducials(AS_BUILT_POLL_STRIDE_POINTS);
        let boundary_error = register_with_poll(&boundary, &mut |phase, completed, _| {
            if completed == u128::try_from(AS_BUILT_POLL_STRIDE_POINTS).unwrap() {
                Err(fs_exec::Cancelled)
            } else {
                let _ = phase;
                Ok(())
            }
        })
        .expect_err("cancellation at the phase boundary must suppress publication");
        assert_eq!(
            boundary_error,
            RegError::Cancelled {
                phase: REGISTER_DESIGN_CENTROID_PHASE,
                completed_work: 256,
                planned_work: 1_536,
            }
        );

        let plus_one = fiducials(AS_BUILT_POLL_STRIDE_POINTS + 1);
        let plus_one_error = register_with_poll(&plus_one, &mut |phase, completed, _| {
            if phase == REGISTER_DESIGN_CENTROID_PHASE && completed == 256 {
                Err(fs_exec::Cancelled)
            } else {
                Ok(())
            }
        })
        .expect_err("the stride-plus-one scan must poll before its last point");
        assert_eq!(
            plus_one_error,
            RegError::Cancelled {
                phase: REGISTER_DESIGN_CENTROID_PHASE,
                completed_work: 256,
                planned_work: 1_542,
            }
        );

        let second_pass_error = register_with_poll(&plus_one, &mut |phase, completed, _| {
            let second_pass_stride =
                u128::try_from(plus_one.len() + AS_BUILT_POLL_STRIDE_POINTS).unwrap();
            if phase == REGISTER_DESIGN_CENTROID_PHASE && completed == second_pass_stride {
                Err(fs_exec::Cancelled)
            } else {
                Ok(())
            }
        })
        .expect_err("the anchored-normalized pass must poll at its own stride boundary");
        assert_eq!(
            second_pass_error,
            RegError::Cancelled {
                phase: REGISTER_DESIGN_CENTROID_PHASE,
                completed_work: 513,
                planned_work: 1_542,
            }
        );
    }

    #[test]
    fn g4_hostile_maximum_shapes_preflight_before_initial_cancellation() {
        let point = Point2::new(0.0, 0.0).expect("finite hostile-maximum point");
        let maximum_fiducials = vec![Fiducial::new(point, point); MAX_AS_BUILT_POINTS];
        let registration_error = register_with_poll(
            &maximum_fiducials,
            &mut |phase, completed_work, planned_work| {
                assert_eq!(phase, REGISTER_INITIAL_PHASE);
                assert_eq!(completed_work, 0);
                assert_eq!(planned_work, 6_000_000);
                Err(fs_exec::Cancelled)
            },
        )
        .expect_err("maximum registration must remain cancellable before scalar work");
        assert_eq!(
            registration_error,
            RegError::Cancelled {
                phase: REGISTER_INITIAL_PHASE,
                completed_work: 0,
                planned_work: 6_000_000,
            }
        );

        let maximum_points = vec![point; MAX_AS_BUILT_POINTS];
        let registration = Registration::new(0.0, 0.0, 0.0, 0.0).unwrap();
        let diff_error = as_built_diff_with_poll(
            &registration,
            &maximum_points,
            &maximum_points,
            1.0,
            0.1,
            "hostile-maximum-cancel-fixture",
            execution(),
            CURRENT_POLL_POLICY,
            &mut |phase, completed_work, planned_work| {
                assert_eq!(phase, DIFF_INITIAL_PHASE);
                assert_eq!(completed_work, 0);
                assert_eq!(planned_work, 3_000_060);
                Err(fs_exec::Cancelled)
            },
        )
        .expect_err("maximum diff must remain cancellable before allocation or scalar work");
        assert_eq!(
            diff_error,
            RegError::Cancelled {
                phase: DIFF_INITIAL_PHASE,
                completed_work: 0,
                planned_work: 3_000_060,
            }
        );
    }

    #[test]
    fn g4_mid_diff_and_final_publication_cancellation_are_transactional() {
        let registration_fiducials = [
            Fiducial::new(
                Point2::new(0.0, 0.0).unwrap(),
                Point2::new(1.0, 1.0).unwrap(),
            ),
            Fiducial::new(
                Point2::new(2.0, 0.0).unwrap(),
                Point2::new(3.0, 1.0).unwrap(),
            ),
            Fiducial::new(
                Point2::new(0.0, 2.0).unwrap(),
                Point2::new(1.0, 3.0).unwrap(),
            ),
        ];
        let registration_final_error =
            register_with_poll(&registration_fiducials, &mut |phase, _, _| {
                if phase == REGISTER_PUBLISH_PHASE {
                    Err(fs_exec::Cancelled)
                } else {
                    Ok(())
                }
            })
            .expect_err("the registration final checkpoint must precede publication");
        assert_eq!(
            registration_final_error,
            RegError::Cancelled {
                phase: REGISTER_PUBLISH_PHASE,
                completed_work: 18,
                planned_work: 18,
            }
        );

        let reg = Registration::new(0.0, 0.0, 0.0, 0.0).unwrap();
        let design = points(AS_BUILT_POLL_STRIDE_POINTS + 1);
        let mid_error = as_built_diff_with_poll(
            &reg,
            &design,
            &design,
            1.0,
            0.1,
            "mid-cancel-fixture",
            execution(),
            CURRENT_POLL_POLICY,
            &mut |phase, completed, _| {
                if phase == DIFF_DEVIATIONS_PHASE && completed == 274 {
                    Err(fs_exec::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
        .expect_err("mid-scan cancellation must return no partial normal output");
        assert_eq!(
            mid_error,
            RegError::Cancelled {
                phase: DIFF_DEVIATIONS_PHASE,
                completed_work: 274,
                planned_work: 807,
            }
        );

        let one = [Point2::new(0.0, 0.0).unwrap()];
        let final_error = as_built_diff_with_poll(
            &reg,
            &one,
            &one,
            1.0,
            0.1,
            "final-cancel-fixture",
            execution(),
            CURRENT_POLL_POLICY,
            &mut |phase, _, _| {
                if phase == DIFF_PUBLISH_PHASE {
                    Err(fs_exec::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
        .expect_err("the final checkpoint must precede authoritative publication");
        assert_eq!(
            final_error,
            RegError::Cancelled {
                phase: DIFF_PUBLISH_PHASE,
                completed_work: 43,
                planned_work: 43,
            }
        );
    }

    #[test]
    fn g4_calibration_byte_boundaries_have_exact_phase_progress() {
        let reg = Registration::new(0.0, 0.0, 0.0, 0.0).unwrap();
        let one = [Point2::new(0.0, 0.0).unwrap()];
        let calibration = "x".repeat(AS_BUILT_POLL_STRIDE_BYTES);

        let validation_error = as_built_diff_with_poll(
            &reg,
            &one,
            &one,
            1.0,
            0.1,
            &calibration,
            execution(),
            CURRENT_POLL_POLICY,
            &mut |phase, completed, _| {
                if phase == DIFF_CALIBRATION_VALIDATION_PHASE && completed == 256 {
                    Err(fs_exec::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
        .expect_err("validation must poll after the maximum bounded identity scan");
        assert_eq!(
            validation_error,
            RegError::Cancelled {
                phase: DIFF_CALIBRATION_VALIDATION_PHASE,
                completed_work: 256,
                planned_work: 515,
            }
        );

        let hash_error = as_built_diff_with_poll(
            &reg,
            &one,
            &one,
            1.0,
            0.1,
            &calibration,
            execution(),
            CURRENT_POLL_POLICY,
            &mut |phase, completed, _| {
                if phase == DIFF_CALIBRATION_HASH_PHASE && completed == 514 {
                    Err(fs_exec::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
        .expect_err("identity hashing must poll after the maximum bounded identity scan");
        assert_eq!(
            hash_error,
            RegError::Cancelled {
                phase: DIFF_CALIBRATION_HASH_PHASE,
                completed_work: 514,
                planned_work: 515,
            }
        );
    }

    #[test]
    fn g0_work_progress_fails_closed_when_runtime_exceeds_preflight() {
        let mut progress = WorkProgress::new(0, "short-plan-fixture");
        assert_eq!(
            progress.complete_point(),
            Err(RegError::WorkPlanMismatch {
                phase: "short-plan-fixture",
                completed_work: 1,
                planned_work: 0,
            })
        );
    }

    #[test]
    fn g5_poll_policy_version_and_stride_change_identity_not_numerics() {
        let reg = Registration::new(0.0, 0.0, 0.0, 0.0).unwrap();
        let design = [
            Point2::new(0.0, 0.0).unwrap(),
            Point2::new(1.0, 2.0).unwrap(),
        ];
        let scanned = [
            Point2::new(0.1, 0.0).unwrap(),
            Point2::new(1.0, 2.0).unwrap(),
        ];
        let run = |poll_policy| {
            as_built_diff_with_poll(
                &reg,
                &design,
                &scanned,
                0.2,
                0.05,
                "policy-identity-fixture",
                execution(),
                poll_policy,
                &mut |_, _, _| Ok(()),
            )
            .expect("non-cancelled policy fixture")
        };

        let baseline = run(CURRENT_POLL_POLICY);
        let next_version = run(PollPolicy {
            version: CURRENT_POLL_POLICY.version + 1,
            ..CURRENT_POLL_POLICY
        });
        let next_stride = run(PollPolicy {
            stride_points: CURRENT_POLL_POLICY.stride_points + 1,
            ..CURRENT_POLL_POLICY
        });
        let next_byte_stride = run(PollPolicy {
            stride_bytes: CURRENT_POLL_POLICY.stride_bytes + 1,
            ..CURRENT_POLL_POLICY
        });
        assert_ne!(identity(&baseline), identity(&next_version));
        assert_ne!(identity(&baseline), identity(&next_stride));
        assert_ne!(identity(&baseline), identity(&next_byte_stride));
        assert_eq!(baseline.deviations, next_version.deviations);
        assert_eq!(baseline.deviations, next_stride.deviations);
        assert_eq!(baseline.deviations, next_byte_stride.deviations);
        assert_eq!(
            baseline.max_deviation.to_bits(),
            next_version.max_deviation.to_bits()
        );
        assert_eq!(
            baseline.max_deviation.to_bits(),
            next_stride.max_deviation.to_bits()
        );
        assert_eq!(
            baseline.max_deviation.to_bits(),
            next_byte_stride.max_deviation.to_bits()
        );
    }
}
