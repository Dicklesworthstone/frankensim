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
//! - REGISTRATION (aligning scan to design) is an OPTIMIZATION WITH ERROR, so
//!   its error is carried forward, not discarded. [`register`] solves the rigid
//!   2-D fit in closed form (no SVD) and is made WELL-POSED by fiducials/datums
//!   specified at design time (≥ 3 non-collinear points) — the
//!   design-for-verification requirement pushed upstream.
//! - The R8 kill criterion is explicit: if registration uncertainty exceeds the
//!   geometric deviation being certified, the signal is below the noise floor
//!   ([`well_posed`]).
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

const AS_BUILT_ESTIMATOR_DOMAIN: &str = "org.frankensim.fs-asbuilt.diff-estimator.v3";
const AS_BUILT_ESTIMATOR_SCHEMA: &[u8] = b"fs-asbuilt-diff-estimator-v3";
const AS_BUILT_WORK_PLAN_VERSION: u32 = 1;
/// Identity-bound cancellation policy version for resource-driving point scans.
pub const AS_BUILT_POLL_POLICY_VERSION: u32 = 1;
/// Maximum point visits between cancellation polls inside a complete scan.
pub const AS_BUILT_POLL_STRIDE_POINTS: usize = 256;
/// Maximum points accepted by registration or one as-built comparison.
pub const MAX_AS_BUILT_POINTS: usize = 1_000_000;

const REGISTER_INITIAL_PHASE: &str = "register.initial";
const REGISTER_DESIGN_CENTROID_PHASE: &str = "register.design-centroid";
const REGISTER_MEASURED_CENTROID_PHASE: &str = "register.measured-centroid";
const REGISTER_SCATTER_PHASE: &str = "register.scatter";
const REGISTER_RESIDUAL_PHASE: &str = "register.residual";
const REGISTER_PUBLISH_PHASE: &str = "register.publish";
const DIFF_INITIAL_PHASE: &str = "as-built-diff.initial";
const DIFF_DEVIATIONS_PHASE: &str = "as-built-diff.deviations";
const DIFF_MAXIMUM_PHASE: &str = "as-built-diff.maximum";
const DIFF_IDENTITY_PHASE: &str = "as-built-diff.identity";
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
        /// Exact point-visit work completed before the checkpoint.
        completed_work: u128,
        /// Exact point-visit work planned by the constant-time preflight.
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
    /// The complete point-visit work plan could not be represented exactly.
    WorkPlanOverflow {
        /// Stable operation name.
        operation: &'static str,
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
                "as-built operation cancelled during {phase} after {completed_work}/{planned_work} point visits"
            ),
            Self::TooFewFiducials { have, need } => {
                write!(formatter, "need at least {need} fiducials, got {have}")
            }
            Self::CollinearFiducials => formatter.write_str("fiducials are collinear"),
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
                write!(formatter, "{operation} point-visit work plan exceeds u128")
            }
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
    /// Root-mean-square residual of the fit (the registration uncertainty).
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
        Point2::new(
            c * point.x - s * point.y + self.tx,
            s * point.x + c * point.y + self.ty,
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
            .checked_mul(4)
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
        let total = points_per_scan
            .checked_mul(3)
            .ok_or(RegError::WorkPlanOverflow {
                operation: "as-built-diff",
            })?;
        Ok(Self {
            points_per_scan,
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

    fn complete_point(&mut self) -> Result<(), RegError> {
        self.completed = self
            .completed
            .checked_add(1)
            .filter(|completed| *completed <= self.planned)
            .ok_or(RegError::WorkPlanOverflow {
                operation: self.operation,
            })?;
        Ok(())
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
/// [`RegError::TooFewFiducials`], [`RegError::CollinearFiducials`], or a
/// structured [`RegError::Cancelled`] with exact point-visit progress. The
/// complete work plan is computed before the initial cancellation checkpoint.
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
        fiducials.iter().map(|f| f.design),
        REGISTER_DESIGN_CENTROID_PHASE,
        &mut progress,
        poll,
    )?;
    debug_assert_eq!(progress.completed, plan.points_per_scan);
    operation_checkpoint(REGISTER_MEASURED_CENTROID_PHASE, progress, poll)?;
    let cq = centroid(
        fiducials.iter().map(|f| f.measured),
        REGISTER_MEASURED_CENTROID_PHASE,
        &mut progress,
        poll,
    )?;
    debug_assert_eq!(progress.completed, plan.points_per_scan * 2);

    // scatter of the centered DESIGN points — collinear iff it is rank-deficient.
    let (mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0);
    // cross-covariance terms for the optimal rotation.
    let (mut s_dot, mut s_cross) = (0.0, 0.0);
    operation_checkpoint(REGISTER_SCATTER_PHASE, progress, poll)?;
    for (index, f) in fiducials.iter().enumerate() {
        scan_checkpoint(
            index,
            AS_BUILT_POLL_STRIDE_POINTS,
            REGISTER_SCATTER_PHASE,
            progress,
            poll,
        )?;
        let (dpx, dpy) = (f.design.x - cp.x, f.design.y - cp.y);
        let (dqx, dqy) = (f.measured.x - cq.x, f.measured.y - cq.y);
        sxx += dpx * dpx;
        syy += dpy * dpy;
        sxy += dpx * dpy;
        s_dot += dpx * dqx + dpy * dqy;
        s_cross += dpx * dqy - dpy * dqx;
        progress.complete_point()?;
    }
    debug_assert_eq!(progress.completed, plan.points_per_scan * 3);
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

    let rotation_rad = s_cross.atan2(s_dot);
    let (s, c) = rotation_rad.sin_cos();
    let tx = cq.x - (c * cp.x - s * cp.y);
    let ty = cq.y - (s * cp.x + c * cp.y);
    let reg = Registration::new(rotation_rad, tx, ty, 0.0)?;
    // residual RMS = the carried-forward registration uncertainty.
    let mut ss = 0.0;
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
        ss += distance * distance;
        require_finite("registration residual sum of squares", ss)?;
        progress.complete_point()?;
    }
    debug_assert_eq!(progress.completed, plan.total);
    operation_checkpoint(REGISTER_PUBLISH_PHASE, progress, poll)?;
    Registration::new(rotation_rad, tx, ty, (ss / nf).sqrt())
}

fn centroid(
    points: impl Iterator<Item = Point2>,
    phase: &'static str,
    progress: &mut WorkProgress,
    poll: &mut impl FnMut(&'static str, u128, u128) -> Result<(), fs_exec::Cancelled>,
) -> Result<Point2, RegError> {
    let mut n = 0.0;
    let (mut sx, mut sy) = (0.0, 0.0);
    for (index, point) in points.enumerate() {
        scan_checkpoint(index, AS_BUILT_POLL_STRIDE_POINTS, phase, *progress, poll)?;
        sx += point.x;
        sy += point.y;
        n += 1.0;
        require_finite("point centroid x sum", sx)?;
        require_finite("point centroid y sum", sy)?;
        progress.complete_point()?;
    }
    Point2::new(sx / n, sy / n)
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
}

const CURRENT_POLL_POLICY: PollPolicy = PollPolicy {
    version: AS_BUILT_POLL_POLICY_VERSION,
    stride_points: AS_BUILT_POLL_STRIDE_POINTS,
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

    if let Some(reason) = color_leaf_identity_reason(calibration_candidate) {
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
    debug_assert_eq!(progress.completed, plan.total);
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
        b"work-plan.total-point-visits",
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
    field(&mut hasher, b"calibration-candidate")?;
    field(&mut hasher, calibration_candidate.as_bytes())?;
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
    let preimage_hash = hasher.finalize();
    Ok(format!(
        "asbuilt-diff-v3:{}",
        fs_blake3::hash_domain(AS_BUILT_ESTIMATOR_DOMAIN, preimage_hash.as_bytes())
    ))
}

#[cfg(test)]
mod cancellation_tests {
    use super::*;

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
                phase: REGISTER_MEASURED_CENTROID_PHASE,
                completed_work: 256,
                planned_work: 1_024,
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
                planned_work: 1_028,
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
                completed_work: 12,
                planned_work: 12,
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
                if phase == DIFF_DEVIATIONS_PHASE && completed == 256 {
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
                completed_work: 256,
                planned_work: 771,
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
                completed_work: 3,
                planned_work: 3,
            }
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
        assert_ne!(identity(&baseline), identity(&next_version));
        assert_ne!(identity(&baseline), identity(&next_stride));
        assert_eq!(baseline.deviations, next_version.deviations);
        assert_eq!(baseline.deviations, next_stride.deviations);
        assert_eq!(baseline.max_deviation, next_version.max_deviation);
        assert_eq!(baseline.max_deviation, next_stride.max_deviation);
    }
}
