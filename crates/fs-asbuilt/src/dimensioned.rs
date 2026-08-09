//! Dimensional adoption for the as-built lane (bead sj31i.7.2): explicit
//! frames, typed length units, and affine/linear distinctions for
//! coordinates, registration transforms, residuals, tolerances, and
//! uncertainty evidence.
//!
//! The raw [`crate::Point2`]/[`crate::Fiducial`]/[`crate::register`] path
//! stays untouched. This module is the dimensionally typed front door:
//!
//! - a [`FramePoint`] is an AFFINE position bound to a named [`FrameId`]
//!   and [`LengthUnit`]; a [`Displacement`] is a LINEAR quantity. Point
//!   minus point is a displacement, displacement plus point is a point,
//!   and point plus point refuses — there is no implicit unit equivalence
//!   and no silent frame crossing;
//! - millimetres and metres mix only through explicit typed conversion;
//! - a [`DimensionedRegistration`] binds the design/measured frames, the
//!   registration unit, the transform, and the residual into one receipt
//!   identity (`as-built-registration:v1:<hex>`);
//! - uncertainty covariance over the planar state is derived mechanically
//!   through `fs_qty::inference::CovarianceSchema` (length-squared
//!   entries), never caller-declared.
//!
//! Work is admitted and polled through the ambient `fs_exec::Cx` exactly
//! like the raw registration path.

use fs_qty::inference::{CovarianceSchema, SlotSchema, StateSchema};
use fs_qty::semantic::QuantitySpec;

use crate::{Fiducial, Point2, RegError, Registration, register};

/// Length unit for as-built coordinates. The scale factor maps one unit to
/// coherent SI metres; conversion happens ONLY at explicit typed calls.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthUnit {
    /// SI metres.
    Meters,
    /// Millimetres (0.001 m).
    Millimeters,
}

impl LengthUnit {
    /// Coherent-SI scale factor for one unit.
    #[must_use]
    pub const fn scale_to_si(self) -> f64 {
        match self {
            Self::Meters => 1.0,
            Self::Millimeters => 1.0e-3,
        }
    }

    /// Stable wire name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Meters => "m",
            Self::Millimeters => "mm",
        }
    }
}

/// A bounded coordinate-frame identity token.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct FrameId(String);

impl FrameId {
    /// Admit a frame identity: 1..=64 ASCII alphanumeric, `-`, or `_`.
    ///
    /// # Errors
    /// Returns [`RegError::InvalidIdentity`] for empty, oversized, or
    /// non-token identities.
    pub fn try_new(identity: impl Into<String>) -> Result<Self, RegError> {
        let identity = identity.into();
        let valid = !identity.is_empty()
            && identity.len() <= 64
            && identity
                .bytes()
                .all(|byte| byte.is_ascii_alphanumeric() || byte == b'-' || byte == b'_');
        if !valid {
            return Err(RegError::InvalidFrameIdentity {
                reason: "expected 1..=64 ASCII alphanumeric, '-', or '_' bytes",
                bytes: identity.len(),
            });
        }
        Ok(Self(identity))
    }

    /// The identity string.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// An affine position bound to a frame and a length unit. Points support
/// exactly the affine operations: difference (a displacement), translation
/// by a displacement, and explicit unit conversion.
#[derive(Debug, Clone, PartialEq)]
pub struct FramePoint {
    frame: FrameId,
    unit: LengthUnit,
    x: f64,
    y: f64,
}

impl FramePoint {
    /// Admit a finite frame point.
    ///
    /// # Errors
    /// Returns [`RegError`] for non-finite coordinates.
    pub fn new(frame: FrameId, unit: LengthUnit, x: f64, y: f64) -> Result<Self, RegError> {
        Point2::new(x, y)?;
        Ok(Self {
            frame,
            unit,
            x: if x == 0.0 { 0.0 } else { x },
            y: if y == 0.0 { 0.0 } else { y },
        })
    }

    /// The bound frame.
    #[must_use]
    pub const fn frame(&self) -> &FrameId {
        &self.frame
    }

    /// The declared unit.
    #[must_use]
    pub const fn unit(&self) -> LengthUnit {
        self.unit
    }

    /// x coordinate in the declared unit.
    #[must_use]
    pub const fn x(&self) -> f64 {
        self.x
    }

    /// y coordinate in the declared unit.
    #[must_use]
    pub const fn y(&self) -> f64 {
        self.y
    }

    /// Explicit typed unit conversion; the frame is preserved.
    #[must_use]
    pub fn to_unit(&self, unit: LengthUnit) -> Self {
        if unit == self.unit {
            return self.clone();
        }
        let factor = self.unit.scale_to_si() / unit.scale_to_si();
        Self {
            frame: self.frame.clone(),
            unit,
            x: self.x * factor,
            y: self.y * factor,
        }
    }

    /// Affine difference `self - other` as a displacement in this point's
    /// unit. Frames must match; units are converted explicitly.
    ///
    /// # Errors
    /// Returns [`RegError::FrameMismatch`] when the frames differ, or a
    /// finite-computation refusal when the difference overflows.
    pub fn difference(&self, other: &Self) -> Result<Displacement, RegError> {
        if self.frame != other.frame {
            return Err(RegError::FrameMismatch {
                expected: self.frame.as_str().to_string(),
                actual: other.frame.as_str().to_string(),
            });
        }
        let other = other.to_unit(self.unit);
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        if !dx.is_finite() || !dy.is_finite() {
            return Err(RegError::NonFinite { field: "point difference" });
        }
        Ok(Displacement {
            unit: self.unit,
            dx,
            dy,
        })
    }

    /// Affine translation by a displacement; the result stays in this
    /// point's frame and unit.
    ///
    /// # Errors
    /// Returns [`RegError`] on non-finite results.
    pub fn translate(&self, displacement: Displacement) -> Result<Self, RegError> {
        let displacement = displacement.to_unit(self.unit);
        let x = self.x + displacement.dx;
        let y = self.y + displacement.dy;
        Self::new(self.frame.clone(), self.unit, x, y)
    }
}

/// A linear planar quantity in a named unit. Displacements add, subtract,
/// and scale freely; they never carry a frame or an origin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Displacement {
    unit: LengthUnit,
    dx: f64,
    dy: f64,
}

impl Displacement {
    /// Admit a finite displacement.
    ///
    /// # Errors
    /// Returns [`RegError`] for non-finite components.
    pub fn new(unit: LengthUnit, dx: f64, dy: f64) -> Result<Self, RegError> {
        if !dx.is_finite() || !dy.is_finite() {
            return Err(RegError::NonFinite {
                field: "displacement component",
            });
        }
        Ok(Self { unit, dx, dy })
    }

    /// The declared unit.
    #[must_use]
    pub const fn unit(&self) -> LengthUnit {
        self.unit
    }

    /// x component in the declared unit.
    #[must_use]
    pub const fn dx(&self) -> f64 {
        self.dx
    }

    /// y component in the declared unit.
    #[must_use]
    pub const fn dy(&self) -> f64 {
        self.dy
    }

    /// Explicit typed unit conversion.
    #[must_use]
    pub fn to_unit(self, unit: LengthUnit) -> Self {
        if unit == self.unit {
            return self;
        }
        let factor = self.unit.scale_to_si() / unit.scale_to_si();
        Self {
            unit,
            dx: self.dx * factor,
            dy: self.dy * factor,
        }
    }

    /// Linear addition; units are converted explicitly.
    #[must_use]
    pub fn plus(self, other: Self) -> Self {
        let other = other.to_unit(self.unit);
        Self {
            unit: self.unit,
            dx: self.dx + other.dx,
            dy: self.dy + other.dy,
        }
    }

    /// Euclidean magnitude in the declared unit.
    #[must_use]
    pub fn magnitude(self) -> f64 {
        self.dx.hypot(self.dy)
    }
}

/// A length-dimensioned tolerance bound in a named unit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tolerance {
    value: f64,
    unit: LengthUnit,
}

impl Tolerance {
    /// Admit a finite, non-negative tolerance.
    ///
    /// # Errors
    /// Returns [`RegError`] for negative or non-finite values.
    pub fn new(value: f64, unit: LengthUnit) -> Result<Self, RegError> {
        if !value.is_finite() || value < 0.0 {
            return Err(RegError::NonFinite {
                field: "tolerance value",
            });
        }
        Ok(Self { value, unit })
    }

    /// The tolerance value in the declared unit.
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }

    /// The declared unit.
    #[must_use]
    pub const fn unit(&self) -> LengthUnit {
        self.unit
    }

    /// True when the displacement's magnitude fits the tolerance after
    /// explicit unit conversion.
    #[must_use]
    pub fn contains(&self, displacement: Displacement) -> bool {
        let displacement = displacement.to_unit(self.unit);
        displacement.magnitude() <= self.value
    }
}

/// A fiducial correspondence with typed endpoints: a design-frame point and
/// a measured-frame point. The two frames may differ (registration maps
/// between them); units convert explicitly at the registration boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct DimensionedFiducial {
    design: FramePoint,
    measured: FramePoint,
}

impl DimensionedFiducial {
    /// Admit one correspondence.
    #[must_use]
    pub const fn new(design: FramePoint, measured: FramePoint) -> Self {
        Self { design, measured }
    }

    /// Design-frame endpoint.
    #[must_use]
    pub const fn design(&self) -> &FramePoint {
        &self.design
    }

    /// Measured-frame endpoint.
    #[must_use]
    pub const fn measured(&self) -> &FramePoint {
        &self.measured
    }
}

/// A dimensionally bound registration product: the transform, the design
/// and measured frames, the registration unit, the residual, and the
/// receipt identity.
#[derive(Debug, Clone, PartialEq)]
pub struct DimensionedRegistration {
    registration: Registration,
    design_frame: FrameId,
    measured_frame: FrameId,
    unit: LengthUnit,
    identity: String,
}

impl DimensionedRegistration {
    /// The raw registration transform (rotation, translation) expressed in
    /// the registration unit.
    #[must_use]
    pub const fn registration(&self) -> &Registration {
        &self.registration
    }

    /// The design frame.
    #[must_use]
    pub const fn design_frame(&self) -> &FrameId {
        &self.design_frame
    }

    /// The measured frame.
    #[must_use]
    pub const fn measured_frame(&self) -> &FrameId {
        &self.measured_frame
    }

    /// The registration unit.
    #[must_use]
    pub const fn unit(&self) -> LengthUnit {
        self.unit
    }

    /// The residual RMS as a typed tolerance-compatible length.
    #[must_use]
    pub fn residual_length(&self) -> f64 {
        self.registration.residual_rms()
    }

    /// The receipt identity `as-built-registration:v1:<64 lowercase hex>`.
    #[must_use]
    pub fn identity(&self) -> &str {
        &self.identity
    }
}

/// Wire prefix of the dimensioned registration receipt identity.
pub const REGISTRATION_RECEIPT_PREFIX: &str = "as-built-registration:v1:";

const RECEIPT_DOMAIN: &str = "org.frankensim.fs-asbuilt.registration.v1";

fn registration_identity(
    design_frame: &FrameId,
    measured_frame: &FrameId,
    unit: LengthUnit,
    fiducial_count: usize,
    registration: &Registration,
) -> String {
    let mut hasher = fs_blake3::DomainHasher::new(RECEIPT_DOMAIN);
    hasher.update(&(design_frame.as_str().len() as u64).to_le_bytes());
    hasher.update(design_frame.as_str().as_bytes());
    hasher.update(&(measured_frame.as_str().len() as u64).to_le_bytes());
    hasher.update(measured_frame.as_str().as_bytes());
    hasher.update(unit.name().as_bytes());
    hasher.update(&(fiducial_count as u64).to_le_bytes());
    hasher.update(&registration.rotation_rad().to_le_bytes());
    hasher.update(&registration.tx().to_le_bytes());
    hasher.update(&registration.ty().to_le_bytes());
    hasher.update(&registration.residual_rms().to_le_bytes());
    format!("{REGISTRATION_RECEIPT_PREFIX}{}", hasher.finalize())
}

/// Register dimensioned fiducials through the production rigid fit:
/// every endpoint converts explicitly into the declared registration unit,
/// design frames must agree and measured frames must agree (mixed-frame
/// batches refuse), and the raw [`register`] runs on the converted points.
/// The returned product binds frames, unit, transform, residual, and the
/// receipt identity.
///
/// # Errors
/// Returns [`RegError::FrameMismatch`] for mixed design or measured frames,
/// and every refusal of the production registration (too few fiducials,
/// collinearity, cancellation, budgets).
pub fn register_dimensioned(
    fiducials: &[DimensionedFiducial],
    unit: LengthUnit,
    cx: &fs_exec::Cx<'_>,
) -> Result<DimensionedRegistration, RegError> {
    let Some(first) = fiducials.first() else {
        return Err(RegError::TooFewFiducials { have: 0, need: 2 });
    };
    let design_frame = first.design().frame().clone();
    let measured_frame = first.measured().frame().clone();
    let mut raw = Vec::with_capacity(fiducials.len());
    for fiducial in fiducials {
        if fiducial.design().frame() != &design_frame {
            return Err(RegError::FrameMismatch {
                expected: design_frame.as_str().to_string(),
                actual: fiducial.design().frame().as_str().to_string(),
            });
        }
        if fiducial.measured().frame() != &measured_frame {
            return Err(RegError::FrameMismatch {
                expected: measured_frame.as_str().to_string(),
                actual: fiducial.measured().frame().as_str().to_string(),
            });
        }
        let design = fiducial.design().to_unit(unit);
        let measured = fiducial.measured().to_unit(unit);
        raw.push(Fiducial::new(
            Point2::new(design.x(), design.y())?,
            Point2::new(measured.x(), measured.y())?,
        ));
    }
    let registration = register(&raw, cx)?;
    let identity = registration_identity(
        &design_frame,
        &measured_frame,
        unit,
        fiducials.len(),
        &registration,
    );
    Ok(DimensionedRegistration {
        registration,
        design_frame,
        measured_frame,
        unit,
        identity,
    })
}

/// The mechanical covariance schema over the planar as-built state
/// `[x, y]`: both slots are length, so every covariance entry carries
/// length-squared dimensions and every information entry carries
/// inverse-length-squared dimensions, derived by the shared core.
///
/// # Errors
/// Returns the typed dimensional error when the schema cannot be admitted.
pub fn planar_covariance_schema() -> Result<CovarianceSchema, RegError> {
    let length = SlotSchema::new(QuantitySpec::dimensional(fs_qty::Dims([1, 0, 0, 0, 0, 0])));
    let state = StateSchema::try_new(vec![length, length]).map_err(RegError::Dimensional)?;
    Ok(CovarianceSchema::over(state))
}
