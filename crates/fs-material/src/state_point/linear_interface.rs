//! Source-owned scalar bilateral stiffness and viscosity integrated over area.
//! Geometry, force/displacement maps and time integration remain with consumers.

use super::{MaterialStatePointError, ResolvedInterfaceStatePoint, ScalarAdmissibility};
use fs_matdb::{InterfaceSystemCard, MaterialStateId, QueryPoint};
use fs_qty::{Dims, DynViscosity, Length, Pressure};

/// Dimensions of interface stiffness per area [N/m³], not bulk modulus [Pa].
pub fn interface_stiffness_per_area_dims() -> Dims {
    Pressure::DIMS
        .checked_minus(Length::DIMS)
        .expect("fixed SI interface exponents")
}

/// Dimensions of interface resistance per area [N s/m³], not bulk viscosity.
pub fn interface_viscosity_per_area_dims() -> Dims {
    DynViscosity::DIMS
        .checked_minus(Length::DIMS)
        .expect("fixed SI interface exponents")
}

/// Optional source property: the angular frequency [rad/s] at which the
/// viscosity was written as `eta k / omega_ref` from a measured constant loss
/// factor `eta`. A consumer that honours it applies `eta` at every frequency
/// instead of the Kelvin-Voigt reading `eta omega / omega_ref`. It is resolved
/// only when the consumer requests it, with these dimensions (1/s).
pub const LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA: &str = "constant_loss_reference_omega";

/// Dimensions of [`LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA`]: an angular frequency.
pub fn constant_loss_reference_omega_dims() -> Dims {
    Dims([0, 0, -1, 0, 0, 0])
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}

/// One passive scalar law F = k x + c v at a fixed resolved interface state.
///
/// Source order identifies the actual surface roles and their full receipts;
/// it does not make this bilateral law nonreciprocal. The geometry owner must
/// use the same relative-displacement row and its transpose for force, with
/// opposite reactions on the two sides. No bulk material name, missing loss,
/// contact area, prestress, orientation or rest displacement is inferred here.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLinearInterface {
    resolved: ResolvedInterfaceStatePoint,
    stiffness_property: String,
    viscosity_property: String,
    area_m2: f64,
    stiffness_n_m: f64,
    damping_n_s_m: f64,
    omega_band_rad_s: (f64, f64),
    constant_loss_reference_omega_rad_s: Option<f64>,
}

impl ResolvedLinearInterface {
    /// Integrate original stiffness/viscosity-per-area claims over an explicit
    /// positive physical area. Both properties are required; an explicitly
    /// sourced zero permits a pure spring or damper, while two zeros refuse.
    ///
    /// Validity-wide constants and original source curves/surfaces at held,
    /// source-qualified external conditions are admitted. Temperature, moisture,
    /// composition and imposed pressure remain in the original query receipts.
    /// Frequency-dependent samples are not converted into a constant time-domain
    /// law. Both claims need explicit
    /// frequency applicability; consumers must check their actual modes/band.
    /// The complete original state, ordered materials and usage receipts remain
    /// inspectable. Changing state requires a new source resolution.
    ///
    /// # Errors
    /// Missing/ambiguous properties, wrong dimensions, unsupported source curves,
    /// negative/nonfinite coefficients, absent/disjoint frequency support, or
    /// nonpositive/nonfinite area and unrepresentable integrated coefficients.
    pub fn new(
        source: ResolvedInterfaceStatePoint,
        stiffness_property: &str,
        viscosity_property: &str,
        area_m2: f64,
    ) -> Result<Self, MaterialStatePointError> {
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(invalid(
                "linear interface requires a finite positive physical area",
            ));
        }
        let (k, kb) = coefficient(
            &source,
            stiffness_property,
            interface_stiffness_per_area_dims(),
        )?;
        let (c, cb) = coefficient(
            &source,
            viscosity_property,
            interface_viscosity_per_area_dims(),
        )?;
        let mut band = (kb.0.max(cb.0), kb.1.min(cb.1));
        // The anchor is read only when the resolved bundle carries it, and
        // narrows the band to its own applicability like any coefficient.
        let constant_loss_reference_omega_rad_s =
            match source.property(LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA) {
                None => None,
                Some(_) => {
                    let (w, wb) = coefficient(
                        &source,
                        LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA,
                        constant_loss_reference_omega_dims(),
                    )?;
                    if w <= 0.0 || c == 0.0 {
                        return Err(invalid(
                            "a constant-loss anchor is a positive angular frequency on a nonzero viscosity",
                        ));
                    }
                    band = (band.0.max(wb.0), band.1.min(wb.1));
                    Some(w)
                }
            };
        let stiffness_n_m = k * area_m2;
        let damping_n_s_m = c * area_m2;
        if !stiffness_n_m.is_finite()
            || !damping_n_s_m.is_finite()
            || (k > 0.0 && stiffness_n_m == 0.0)
            || (c > 0.0 && damping_n_s_m == 0.0)
            || (stiffness_n_m == 0.0 && damping_n_s_m == 0.0)
            || !band.0.is_finite()
            || !band.1.is_finite()
            || band.0 < 0.0
            || band.1 < band.0
        {
            return Err(invalid(
                "linear interface needs representable active coefficients and intersected frequency support",
            ));
        }
        Ok(Self {
            resolved: source,
            stiffness_property: stiffness_property.into(),
            viscosity_property: viscosity_property.into(),
            area_m2,
            stiffness_n_m,
            damping_n_s_m,
            omega_band_rad_s: band,
            constant_loss_reference_omega_rad_s,
        })
    }

    /// Original complete ordered interface state, property evidence and receipts.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedInterfaceStatePoint {
        &self.resolved
    }

    /// Rebuild stiffness, damping and supported frequency band at the new
    /// interface state, retaining geometry, property roles and exact sources.
    /// Every property in the original bundle refreshes atomically, including
    /// coefficients consumed by other domains sharing this physical interface.
    pub fn at_query_point(
        &self,
        card: &InterfaceSystemCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::new(
            self.resolved.at_query_point(card, point)?,
            &self.stiffness_property,
            &self.viscosity_property,
            self.area_m2,
        )
    }
    /// Original selected stiffness-per-area property name.
    #[must_use]
    pub fn stiffness_property(&self) -> &str {
        &self.stiffness_property
    }
    /// Original selected viscosity-per-area property name.
    #[must_use]
    pub fn viscosity_property(&self) -> &str {
        &self.viscosity_property
    }
    /// Actual area used once in both coefficient integrals [m²].
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Integrated bilateral stiffness [N/m].
    #[must_use]
    pub const fn stiffness_n_m(&self) -> f64 {
        self.stiffness_n_m
    }
    /// Integrated passive resistance [N s/m].
    #[must_use]
    pub const fn damping_n_s_m(&self) -> f64 {
        self.damping_n_s_m
    }
    /// Intersected claim applicability, with original typed axes converted to rad/s.
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.omega_band_rad_s
    }
    /// The source's [`LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA`] [rad/s]
    /// when resolved, else `None` (Kelvin-Voigt).
    #[must_use]
    pub const fn constant_loss_reference_omega_rad_s(&self) -> Option<f64> {
        self.constant_loss_reference_omega_rad_s
    }

    /// Match exact material states at the actual receiving surfaces, in order.
    /// Texture/orientation and geometric co-location remain consumer obligations.
    pub fn require_surface_materials(
        &self,
        first: &MaterialStateId,
        second: &MaterialStateId,
    ) -> Result<(), MaterialStatePointError> {
        let pair = self.resolved.surface_materials();
        if &pair[0] != first || &pair[1] != second {
            return Err(invalid(
                "linear interface surface material assignment mismatch",
            ));
        }
        Ok(())
    }

    /// Require a physical angular frequency inside both original claim domains.
    /// Calling this on a query frequency alone does not validate a full response.
    pub fn require_angular_frequency_rad_s(
        &self,
        omega: f64,
    ) -> Result<(), MaterialStatePointError> {
        if !omega.is_finite() || omega < self.omega_band_rad_s.0 || omega > self.omega_band_rad_s.1
        {
            return Err(invalid(
                "linear interface frequency is outside original source applicability",
            ));
        }
        Ok(())
    }
}

fn coefficient(
    source: &ResolvedInterfaceStatePoint,
    name: &str,
    dims: Dims,
) -> Result<(f64, (f64, f64)), MaterialStatePointError> {
    let p = source
        .property(name)
        .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
            property: name.into(),
            reason: "missing or ambiguous original linear interface coefficient",
        })?;
    if p.requirement().dims() != dims {
        return Err(MaterialStatePointError::DimensionMismatch {
            property: name.into(),
            expected: dims,
            found: p.requirement().dims(),
        });
    }
    let value = p.value_si();
    if !value.is_finite() || value < 0.0 {
        return Err(MaterialStatePointError::OutsideConsumerDomain {
            property: name.into(),
            value,
            admissibility: ScalarAdmissibility::NonNegative,
        });
    }
    if !p.is_constant_at_fixed_environment() {
        return Err(invalid(
            "linear interface coefficients must be constant at the resolved external conditions",
        ));
    }
    let band = p.angular_frequency_band_rad_s()?.ok_or_else(|| {
        invalid("each linear interface coefficient needs explicit frequency applicability")
    })?;
    Ok((value, band))
}
