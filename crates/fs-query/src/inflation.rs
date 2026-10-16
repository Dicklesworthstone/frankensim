//! Receipt-backed absolute-error inflation for contact-adjacent queries.
//!
//! Contact bounds may only be widened from rigorous conversion, router, or
//! motion-error evidence. The private radius prevents a raw scalar from being
//! mistaken for an admitted receipt, while the arithmetic helpers preserve
//! outward rounding and keep an explicit exact/native zero bit-neutral.

use core::fmt;

use fs_evidence::{Certified, NumericalKind};
use fs_geom::ChainOutcome;

/// A validated nonnegative absolute radius carried into contact bounds.
///
/// Positive values can only be obtained from rigorous conversion or Rep Router
/// receipts. Exact/native geometry makes the explicit [`Self::exact_zero`]
/// assertion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactInflation {
    radius: f64,
}

impl ContactInflation {
    /// Bit-exact zero for geometry used without representation conversion.
    const ZERO: Self = Self { radius: 0.0 };

    /// Explicit exact-zero assertion for a query input with no external
    /// conversion or motion uncertainty.
    #[must_use]
    pub const fn exact_zero() -> Self {
        Self::ZERO
    }

    fn from_certified<T>(receipt: &Certified<T>) -> Result<Self, ContactInflationError> {
        Self::from_receipt(
            receipt.qoi,
            receipt.numerical.kind,
            receipt.numerical.lo,
            receipt.numerical.hi,
        )
    }

    /// Admit a direct converter's unforgeable certified receipt.
    ///
    /// The numerical certificate must describe a finite, ordered,
    /// nonnegative absolute error and contain its finite, nonnegative QoI.
    /// The admitted radius is the certificate's upper endpoint.
    ///
    /// # Errors
    /// [`ContactInflationError`] if the receipt cannot authorize a finite
    /// nonnegative absolute-error radius.
    pub fn from_conversion<T>(receipt: &Certified<T>) -> Result<Self, ContactInflationError> {
        Self::from_certified(receipt)
    }

    /// Admit an unforgeable certified absolute motion-error receipt.
    ///
    /// This is distinct from the deterministic swept radius already supplied
    /// to CCD: it represents uncertainty in that motion model and composes
    /// with any conversion radius before reaching a query bound.
    ///
    /// # Errors
    /// [`ContactInflationError`] if the receipt cannot authorize a finite
    /// nonnegative absolute-error radius.
    pub fn from_motion<T>(receipt: &Certified<T>) -> Result<Self, ContactInflationError> {
        Self::from_certified(receipt)
    }

    fn from_chain_outcome(outcome: &ChainOutcome) -> Result<Self, ContactInflationError> {
        let receipt = outcome.receipt();
        if receipt.value.to_bits() != receipt.qoi.to_bits() {
            return Err(ContactInflationError::ScalarValueMismatch {
                value_bits: receipt.value.to_bits(),
                qoi_bits: receipt.qoi.to_bits(),
            });
        }
        Self::from_receipt(
            receipt.qoi,
            receipt.numerical.kind,
            receipt.numerical.lo,
            receipt.numerical.hi,
        )
    }

    /// Admit an executed Rep Router chain's composed receipt.
    ///
    /// The composed receipt must be Exact/Enclosure-class, finite, ordered,
    /// nonnegative, and contain its finite nonnegative scalar QoI. Estimated,
    /// no-claim, or internally inconsistent routes cannot authorize contact
    /// evidence.
    ///
    /// # Errors
    /// [`ContactInflationError`] if the scalar route receipt is internally
    /// inconsistent, non-rigorous, or malformed.
    pub fn from_route(outcome: &ChainOutcome) -> Result<Self, ContactInflationError> {
        Self::from_chain_outcome(outcome)
    }

    /// Conservatively compose two independent absolute-error radii.
    ///
    /// Either zero operand is returned unchanged, including its exact bit
    /// pattern. Two positive radii are added with upward rounding.
    ///
    /// # Errors
    /// [`ContactInflationError::ArithmeticOverflow`] if the finite sum cannot
    /// be represented.
    pub fn compose(self, other: Self) -> Result<Self, ContactInflationError> {
        if self.radius == 0.0 {
            return Ok(other);
        }
        if other.radius == 0.0 {
            return Ok(self);
        }
        Ok(Self {
            radius: add_upper(self.radius, other.radius)?,
        })
    }

    /// Conservative absolute radius to apply to a bound.
    #[must_use]
    pub const fn radius(self) -> f64 {
        self.radius
    }

    /// Widen an upper endpoint by this radius with upward rounding.
    pub(crate) fn inflate_upper(self, value: f64) -> Result<f64, ContactInflationError> {
        add_upper(value, self.radius)
    }

    /// Widen a lower endpoint by this radius with downward rounding.
    pub(crate) fn deflate_lower(self, value: f64) -> Result<f64, ContactInflationError> {
        subtract_lower(value, self.radius, false)
    }

    /// Widen and clamp a nonnegative lower endpoint at zero.
    pub(crate) fn deflate_nonnegative(self, value: f64) -> Result<f64, ContactInflationError> {
        subtract_lower(value, self.radius, true)
    }

    /// Widen a nonnegative motion or classification radius.
    pub(crate) fn inflate_nonnegative(self, value: f64) -> Result<f64, ContactInflationError> {
        validate_radius(value)?;
        self.inflate_upper(value)
    }

    fn from_receipt(
        qoi: f64,
        kind: NumericalKind,
        lo: f64,
        hi: f64,
    ) -> Result<Self, ContactInflationError> {
        if !matches!(kind, NumericalKind::Exact | NumericalKind::Enclosure) {
            return Err(ContactInflationError::NonRigorousReceipt { kind });
        }
        if !(qoi.is_finite()
            && qoi >= 0.0
            && lo.is_finite()
            && hi.is_finite()
            && lo >= 0.0
            && lo <= qoi
            && qoi <= hi)
        {
            return Err(ContactInflationError::MalformedReceipt {
                qoi_bits: qoi.to_bits(),
                lo_bits: lo.to_bits(),
                hi_bits: hi.to_bits(),
            });
        }
        if kind == NumericalKind::Exact
            && (lo.to_bits() != qoi.to_bits() || hi.to_bits() != qoi.to_bits())
        {
            return Err(ContactInflationError::MalformedReceipt {
                qoi_bits: qoi.to_bits(),
                lo_bits: lo.to_bits(),
                hi_bits: hi.to_bits(),
            });
        }
        Ok(Self {
            // Absolute zero has one canonical representation. This keeps a
            // rigorous `-0.0` endpoint from leaking a distinct zero bit
            // pattern into otherwise bit-neutral native queries.
            radius: if hi == 0.0 { 0.0 } else { hi },
        })
    }
}

impl<T> TryFrom<&Certified<T>> for ContactInflation {
    type Error = ContactInflationError;

    fn try_from(receipt: &Certified<T>) -> Result<Self, Self::Error> {
        Self::from_certified(receipt)
    }
}

impl TryFrom<&ChainOutcome> for ContactInflation {
    type Error = ContactInflationError;

    fn try_from(outcome: &ChainOutcome) -> Result<Self, Self::Error> {
        Self::from_chain_outcome(outcome)
    }
}

/// Why a conversion or motion-error inflation could not be admitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactInflationError {
    /// The receipt carries only estimate/no-claim numerical authority.
    NonRigorousReceipt {
        /// Numerical authority supplied by the receipt.
        kind: NumericalKind,
    },
    /// The receipt's absolute-error enclosure or QoI is structurally invalid.
    MalformedReceipt {
        /// Exact bits of the receipt QoI.
        qoi_bits: u64,
        /// Exact bits of the lower endpoint.
        lo_bits: u64,
        /// Exact bits of the upper endpoint.
        hi_bits: u64,
    },
    /// A scalar router receipt's carried value differs from its stated QoI.
    ScalarValueMismatch {
        /// Exact bits of the carried scalar.
        value_bits: u64,
        /// Exact bits of the certified QoI.
        qoi_bits: u64,
    },
    /// A raw motion/classification radius was not finite and nonnegative.
    InvalidRadius {
        /// Exact bits of the rejected radius.
        radius_bits: u64,
    },
    /// A bound participating in outward arithmetic was non-finite.
    InvalidBound {
        /// Exact bits of the rejected bound.
        bound_bits: u64,
    },
    /// Outward-rounded finite arithmetic overflowed.
    ArithmeticOverflow {
        /// Exact bits of the left operand.
        lhs_bits: u64,
        /// Exact bits of the right operand.
        rhs_bits: u64,
    },
}

impl fmt::Display for ContactInflationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonRigorousReceipt { kind } => write!(
                f,
                "contact inflation requires Exact/Enclosure conversion evidence; got {kind:?}"
            ),
            Self::MalformedReceipt {
                qoi_bits,
                lo_bits,
                hi_bits,
            } => write!(
                f,
                "contact inflation receipt must be finite, ordered, nonnegative, and contain its QoI (qoi={qoi_bits:#018x}, lo={lo_bits:#018x}, hi={hi_bits:#018x})"
            ),
            Self::ScalarValueMismatch {
                value_bits,
                qoi_bits,
            } => write!(
                f,
                "contact inflation route receipt carries scalar {value_bits:#018x} but certifies QoI {qoi_bits:#018x}"
            ),
            Self::InvalidRadius { radius_bits } => write!(
                f,
                "contact motion/classification radius must be finite and nonnegative (bits {radius_bits:#018x})"
            ),
            Self::InvalidBound { bound_bits } => write!(
                f,
                "contact bound arithmetic requires a finite bound (bits {bound_bits:#018x})"
            ),
            Self::ArithmeticOverflow { lhs_bits, rhs_bits } => write!(
                f,
                "contact inflation arithmetic overflowed (lhs={lhs_bits:#018x}, rhs={rhs_bits:#018x})"
            ),
        }
    }
}

impl core::error::Error for ContactInflationError {}

/// Add two finite values and round toward positive infinity.
///
/// A zero right operand is a bit-neutral identity. Callers that use `rhs` as
/// an inflation radius must validate it before this lower-level helper.
fn add_upper(lhs: f64, rhs: f64) -> Result<f64, ContactInflationError> {
    for value in [lhs, rhs] {
        if !value.is_finite() {
            return Err(ContactInflationError::InvalidBound {
                bound_bits: value.to_bits(),
            });
        }
    }
    if rhs == 0.0 {
        return Ok(lhs);
    }
    let sum = lhs + rhs;
    if !sum.is_finite() {
        return Err(ContactInflationError::ArithmeticOverflow {
            lhs_bits: lhs.to_bits(),
            rhs_bits: rhs.to_bits(),
        });
    }
    let upper = sum.next_up();
    if upper.is_finite() {
        Ok(upper)
    } else {
        Err(ContactInflationError::ArithmeticOverflow {
            lhs_bits: lhs.to_bits(),
            rhs_bits: rhs.to_bits(),
        })
    }
}

/// Subtract a nonnegative radius and round toward negative infinity.
///
/// When `clamp_nonnegative` is true, negative results become zero. A zero
/// radius adds no rounding step.
fn subtract_lower(
    lhs: f64,
    radius: f64,
    clamp_nonnegative: bool,
) -> Result<f64, ContactInflationError> {
    if !lhs.is_finite() {
        return Err(ContactInflationError::InvalidBound {
            bound_bits: lhs.to_bits(),
        });
    }
    validate_radius(radius)?;
    let lower = if radius == 0.0 {
        lhs
    } else {
        let difference = lhs - radius;
        if !difference.is_finite() {
            return Err(ContactInflationError::ArithmeticOverflow {
                lhs_bits: lhs.to_bits(),
                rhs_bits: radius.to_bits(),
            });
        }
        difference.next_down()
    };
    Ok(if clamp_nonnegative {
        lower.max(0.0)
    } else {
        lower
    })
}

fn validate_radius(radius: f64) -> Result<(), ContactInflationError> {
    if radius.is_finite() && radius >= 0.0 {
        Ok(())
    } else {
        Err(ContactInflationError::InvalidRadius {
            radius_bits: radius.to_bits(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enclosed(qoi: f64, lo: f64, hi: f64) -> Result<ContactInflation, ContactInflationError> {
        ContactInflation::from_receipt(qoi, NumericalKind::Enclosure, lo, hi)
    }

    #[test]
    fn explicit_exact_zero_is_bit_neutral() {
        let zero = ContactInflation::exact_zero();
        assert_eq!(zero.radius().to_bits(), 0.0f64.to_bits());
        assert_eq!(
            zero.inflate_nonnegative(0.25).unwrap().to_bits(),
            0.25f64.to_bits()
        );

        let admitted = enclosed(0.125, 0.0, 0.125).expect("rigorous receipt");
        assert_eq!(
            zero.compose(admitted).unwrap().radius().to_bits(),
            admitted.radius().to_bits()
        );
        assert_eq!(
            admitted.compose(zero).unwrap().radius().to_bits(),
            admitted.radius().to_bits()
        );
    }

    #[test]
    fn receipt_upper_endpoint_is_admitted_and_composition_rounds_outward() {
        let a = enclosed(0.125, 0.0, 0.25).expect("rigorous receipt");
        let b = enclosed(0.5, 0.25, 0.5).expect("rigorous receipt");
        assert_eq!(a.radius().to_bits(), 0.25f64.to_bits());
        assert!(a.compose(b).unwrap().radius() > 0.75);
    }

    #[test]
    fn malformed_and_non_rigorous_route_receipts_refuse() {
        assert!(matches!(
            ContactInflation::from_receipt(0.25, NumericalKind::Estimate, 0.0, 0.25),
            Err(ContactInflationError::NonRigorousReceipt {
                kind: NumericalKind::Estimate
            })
        ));
        assert!(matches!(
            ContactInflation::from_receipt(
                0.25,
                NumericalKind::NoClaim,
                f64::NEG_INFINITY,
                f64::INFINITY,
            ),
            Err(ContactInflationError::NonRigorousReceipt {
                kind: NumericalKind::NoClaim
            })
        ));

        assert!(matches!(
            enclosed(0.5, 0.0, 0.25),
            Err(ContactInflationError::MalformedReceipt { .. })
        ));

        assert!(matches!(
            enclosed(0.25, 0.0, f64::INFINITY),
            Err(ContactInflationError::MalformedReceipt { .. })
        ));

        assert!(matches!(
            ContactInflation::from_receipt(0.25, NumericalKind::Exact, 0.0, 0.5),
            Err(ContactInflationError::MalformedReceipt { .. })
        ));
    }

    #[test]
    fn outward_helpers_refuse_bad_and_overflowing_inputs() {
        assert!(matches!(
            ContactInflation::exact_zero().inflate_nonnegative(-1.0),
            Err(ContactInflationError::InvalidRadius { .. })
        ));
        assert!(matches!(
            add_upper(f64::MAX, f64::MAX),
            Err(ContactInflationError::ArithmeticOverflow { .. })
        ));
        let half = enclosed(0.5, 0.0, 0.5).expect("rigorous receipt");
        assert_eq!(half.deflate_nonnegative(0.25).unwrap(), 0.0);
    }
}
