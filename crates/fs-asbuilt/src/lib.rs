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
//! The as-built δ ([`as_built_diff`]) is measurement-noise-aware, colored
//! **validated**, and anchored to the metrology instrument's CALIBRATION
//! CERTIFICATE. Deterministic; depends only on `fs-evidence`.

pub use fs_evidence::{Color, ValidityDomain};

/// A 2-D point (design or measured coordinate).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2 {
    /// x coordinate.
    pub x: f64,
    /// y coordinate.
    pub y: f64,
}

impl Point2 {
    /// A point.
    #[must_use]
    pub fn new(x: f64, y: f64) -> Point2 {
        Point2 { x, y }
    }
    fn dist(self, o: Point2) -> f64 {
        ((self.x - o.x).powi(2) + (self.y - o.y).powi(2)).sqrt()
    }
}

/// A fiducial/datum correspondence: a design reference point and where the scan
/// measured it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fiducial {
    /// The design-time reference location.
    pub design: Point2,
    /// The location the scan measured for it.
    pub measured: Point2,
}

impl Fiducial {
    /// A fiducial correspondence.
    #[must_use]
    pub fn new(design: Point2, measured: Point2) -> Fiducial {
        Fiducial { design, measured }
    }
}

/// A structured registration/ingestion failure.
#[derive(Debug, Clone, PartialEq)]
pub enum RegError {
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
}

/// A rigid registration (rotation + translation) mapping design → measured,
/// with the residual it carries forward.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Registration {
    /// Rotation angle (radians).
    pub rotation_rad: f64,
    /// Translation x.
    pub tx: f64,
    /// Translation y.
    pub ty: f64,
    /// Root-mean-square residual of the fit (the registration uncertainty).
    pub residual_rms: f64,
}

impl Registration {
    /// Map a design point into measured coordinates.
    #[must_use]
    pub fn apply(&self, p: Point2) -> Point2 {
        let (s, c) = self.rotation_rad.sin_cos();
        Point2 {
            x: c * p.x - s * p.y + self.tx,
            y: s * p.x + c * p.y + self.ty,
        }
    }
}

/// The minimum fiducials for a well-posed 2-D rigid fit.
pub const MIN_FIDUCIALS: usize = 3;

/// Solve the rigid 2-D registration that best maps `fiducials`' design points
/// onto their measured points (closed-form least squares — the 2-D Umeyama /
/// Procrustes rotation). Requires ≥ 3 non-collinear fiducials.
///
/// # Errors
/// [`RegError::TooFewFiducials`] or [`RegError::CollinearFiducials`].
pub fn register(fiducials: &[Fiducial]) -> Result<Registration, RegError> {
    let n = fiducials.len();
    if n < MIN_FIDUCIALS {
        return Err(RegError::TooFewFiducials {
            have: n,
            need: MIN_FIDUCIALS,
        });
    }
    let nf = n as f64;
    let cp = centroid(fiducials.iter().map(|f| f.design));
    let cq = centroid(fiducials.iter().map(|f| f.measured));

    // scatter of the centered DESIGN points — collinear iff it is rank-deficient.
    let (mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0);
    // cross-covariance terms for the optimal rotation.
    let (mut s_dot, mut s_cross) = (0.0, 0.0);
    for f in fiducials {
        let (dpx, dpy) = (f.design.x - cp.x, f.design.y - cp.y);
        let (dqx, dqy) = (f.measured.x - cq.x, f.measured.y - cq.y);
        sxx += dpx * dpx;
        syy += dpy * dpy;
        sxy += dpx * dpy;
        s_dot += dpx * dqx + dpy * dqy;
        s_cross += dpx * dqy - dpy * dqx;
    }
    let det = sxx * syy - sxy * sxy;
    let trace = sxx + syy;
    if trace <= f64::EPSILON || det <= 1e-12 * trace * trace {
        return Err(RegError::CollinearFiducials);
    }

    let rotation_rad = s_cross.atan2(s_dot);
    let (s, c) = rotation_rad.sin_cos();
    let tx = cq.x - (c * cp.x - s * cp.y);
    let ty = cq.y - (s * cp.x + c * cp.y);
    let reg = Registration {
        rotation_rad,
        tx,
        ty,
        residual_rms: 0.0,
    };
    // residual RMS = the carried-forward registration uncertainty.
    let ss: f64 = fiducials
        .iter()
        .map(|f| reg.apply(f.design).dist(f.measured).powi(2))
        .sum();
    Ok(Registration {
        residual_rms: (ss / nf).sqrt(),
        ..reg
    })
}

fn centroid(pts: impl Iterator<Item = Point2>) -> Point2 {
    let mut n = 0.0;
    let (mut sx, mut sy) = (0.0, 0.0);
    for p in pts {
        sx += p.x;
        sy += p.y;
        n += 1.0;
    }
    Point2 {
        x: sx / n,
        y: sy / n,
    }
}

/// The R8 well-posedness gate: registration is trustworthy only when its
/// uncertainty (`residual_rms`) is BELOW the geometric deviation being
/// certified. If the residual meets or exceeds the signal, the as-built loop is
/// premature for that part class (defer to point-sensor assimilation).
#[must_use]
pub fn well_posed(reg: &Registration, certified_deviation: f64) -> bool {
    certified_deviation > 0.0 && reg.residual_rms < certified_deviation
}

/// The as-built δ between design and scanned sections.
#[derive(Debug, Clone, PartialEq)]
pub struct AsBuiltDiff {
    /// Per-point deviation `||registered(design) − scanned||`.
    pub deviations: Vec<f64>,
    /// The largest deviation.
    pub max_deviation: f64,
    /// Is the whole part within the design tolerance?
    pub within_tolerance: bool,
    /// Is the max deviation ABOVE the measurement noise floor (distinguishable
    /// from noise)?
    pub above_noise_floor: bool,
    /// The δ's color — validated, anchored to the calibration certificate.
    pub color: Color,
}

/// Compute the as-built δ after registration: apply the registration to each
/// design point and measure its deviation from the corresponding scanned point.
/// The δ is colored VALIDATED, its regime tagged with the registration residual
/// and measurement noise, anchored to the metrology `calibration_cert`.
///
/// # Errors
/// [`RegError::Empty`] / [`RegError::LengthMismatch`].
pub fn as_built_diff(
    reg: &Registration,
    design: &[Point2],
    scanned: &[Point2],
    design_tolerance: f64,
    measurement_noise: f64,
    calibration_cert: &str,
) -> Result<AsBuiltDiff, RegError> {
    if design.is_empty() {
        return Err(RegError::Empty);
    }
    if design.len() != scanned.len() {
        return Err(RegError::LengthMismatch {
            expected: design.len(),
            found: scanned.len(),
        });
    }
    let deviations: Vec<f64> = design
        .iter()
        .zip(scanned)
        .map(|(d, s)| reg.apply(*d).dist(*s))
        .collect();
    let max_deviation = deviations.iter().copied().fold(0.0_f64, f64::max);
    let regime = ValidityDomain::unconstrained()
        .with(
            "registration_residual",
            0.0,
            reg.residual_rms.max(f64::MIN_POSITIVE),
        )
        .with(
            "measurement_noise",
            0.0,
            measurement_noise.max(f64::MIN_POSITIVE),
        );
    Ok(AsBuiltDiff {
        deviations,
        max_deviation,
        within_tolerance: max_deviation <= design_tolerance,
        above_noise_floor: max_deviation > measurement_noise,
        color: Color::Validated {
            regime,
            dataset: calibration_cert.to_string(),
        },
    })
}
