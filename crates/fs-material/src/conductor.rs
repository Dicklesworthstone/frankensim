//! Source-resolved uniform DC conductor primitive.
//!
//! The adapter consumes exactly one bulk electrical-resistivity claim and a
//! caller-owned uniform length and cross-section. It does not choose a metal,
//! solve a circuit, add contact resistance, or claim temperature evolution.

use core::fmt;

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::Dims;

use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};

/// Canonical source property for bulk electrical resistivity.
pub const ELECTRICAL_RESISTIVITY_PROPERTY: &str = "electrical_resistivity";

/// Electrical resistivity dimensions `[L, M, t, T, I, N]` for ohm metre.
pub const ELECTRICAL_RESISTIVITY_DIMS: Dims = Dims([3, 1, -3, 0, -2, 0]);

/// A uniform conductor resolved from one source-supported material state.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedConductor {
    material: ResolvedMaterialStatePoint,
    length_m: f64,
    area_m2: f64,
    resistance_ohm: f64,
}

impl ResolvedConductor {
    /// Exact source-resolved material bundle, including its usage receipt.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Uniform DC resistance [ohm].
    #[must_use]
    pub const fn resistance_ohm(&self) -> f64 {
        self.resistance_ohm
    }

    /// Caller-owned uniform conductor length [m].
    #[must_use]
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }

    /// Caller-owned uniform conductor cross-sectional area [m²].
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Resistive Joule power `I²R` [W] for a finite signed current [A].
    ///
    /// Current direction does not affect the dissipated power. This is a
    /// scalar constitutive evaluation, not a circuit or energy-balance solve.
    pub fn joule_power_w(&self, current_a: f64) -> Result<f64, ConductorError> {
        joule_power_w(self.resistance_ohm, current_a)
    }
}

/// Refusal from source-resolved uniform conductor construction or evaluation.
#[derive(Clone, Debug, PartialEq)]
pub enum ConductorError {
    /// The material card could not supply one admissible resistivity claim.
    MaterialState(MaterialStatePointError),
    /// A caller-owned geometric or electrical input was non-finite or invalid.
    InvalidInput {
        /// Stable input name.
        quantity: &'static str,
    },
    /// A finite input combination could not produce a finite physical result.
    InvalidDerived {
        /// Stable result name.
        quantity: &'static str,
    },
}

impl fmt::Display for ConductorError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaterialState(source) => {
                write!(formatter, "conductor material refusal: {source}")
            }
            Self::InvalidInput { quantity } => {
                write!(formatter, "invalid conductor input: {quantity}")
            }
            Self::InvalidDerived { quantity } => {
                write!(formatter, "unrepresentable conductor result: {quantity}")
            }
        }
    }
}

impl std::error::Error for ConductorError {}

impl From<MaterialStatePointError> for ConductorError {
    fn from(source: MaterialStatePointError) -> Self {
        Self::MaterialState(source)
    }
}

fn uniform_resistance_ohm(
    resistivity_ohm_m: f64,
    length_m: f64,
    area_m2: f64,
) -> Result<f64, ConductorError> {
    let direct = resistivity_ohm_m * length_m / area_m2;
    let resistance_ohm = if direct.is_finite() && direct > 0.0 {
        direct
    } else {
        resistivity_ohm_m * (length_m / area_m2)
    };
    if !(resistance_ohm.is_finite() && resistance_ohm > 0.0) {
        return Err(ConductorError::InvalidDerived {
            quantity: "resistance_ohm",
        });
    }
    Ok(resistance_ohm)
}

fn joule_power_w(resistance_ohm: f64, current_a: f64) -> Result<f64, ConductorError> {
    if !current_a.is_finite() {
        return Err(ConductorError::InvalidInput {
            quantity: "current_a",
        });
    }
    let direct = current_a * current_a * resistance_ohm;
    let power_w = if direct.is_finite() && (direct > 0.0 || current_a == 0.0) {
        direct
    } else {
        (current_a * resistance_ohm) * current_a
    };
    if !(power_w.is_finite() && (power_w > 0.0 || current_a == 0.0)) {
        return Err(ConductorError::InvalidDerived {
            quantity: "joule_power_w",
        });
    }
    Ok(power_w)
}

/// Resolve a uniform isotropic DC conductor from one resistivity claim.
///
/// The card must provide `electrical_resistivity` in ohm metre dimensions at
/// the complete caller-supplied query point. Geometry is the actual length and
/// area at that state, not reference geometry to be thermally expanded here;
/// no contact resistance, current distribution, circuit topology, temperature
/// update, or transient energy balance is inferred.
///
/// # Errors
/// Refuses invalid geometry, missing/ambiguous/wrong-unit/out-of-domain source
/// claims, incomplete pin plans, and non-finite resistance.
pub fn resolve_uniform_conductor(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    length_m: f64,
    area_m2: f64,
) -> Result<ResolvedConductor, ConductorError> {
    if !(length_m.is_finite() && length_m > 0.0) {
        return Err(ConductorError::InvalidInput {
            quantity: "length_m",
        });
    }
    if !(area_m2.is_finite() && area_m2 > 0.0) {
        return Err(ConductorError::InvalidInput {
            quantity: "area_m2",
        });
    }
    let requirement = ScalarPropertyRequirement::try_new(
        ELECTRICAL_RESISTIVITY_PROPERTY,
        ELECTRICAL_RESISTIVITY_DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?;
    let material = resolve_material_state_point(card, point, &[requirement], selection)?;
    let resistivity_ohm_m = material
        .property(ELECTRICAL_RESISTIVITY_PROPERTY)
        .expect("the sole required conductor property was resolved")
        .value_si();
    let resistance_ohm = uniform_resistance_ohm(resistivity_ohm_m, length_m, area_m2)?;
    Ok(ResolvedConductor {
        material,
        length_m,
        area_m2,
        resistance_ohm,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_resistance_recovers_after_intermediate_product_overflow() {
        let resistance = uniform_resistance_ohm(1.0e200, 1.0e200, 1.0e200).unwrap();
        assert_eq!(resistance, 1.0e200);
    }

    #[test]
    fn g0_joule_power_recovers_after_current_square_underflow() {
        let power = joule_power_w(1.0e200, 1.0e-200).unwrap();
        assert_eq!(power, 1.0e-200);
    }

    #[test]
    fn g0_joule_power_refuses_unrepresentable_nonzero_result() {
        assert_eq!(
            joule_power_w(f64::MIN_POSITIVE, f64::MIN_POSITIVE),
            Err(ConductorError::InvalidDerived {
                quantity: "joule_power_w"
            })
        );
    }
}
