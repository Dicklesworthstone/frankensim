//! Original source coefficients for coupled Hertz sliding and torsional spin.
//!
//! The thermal Coulomb owner retains the ordered interface, pinned claims,
//! explicit temperature/pressure axes and all fixed external conditions. This
//! bridge supplies their actual values to the original `HertzWrenchContact`.
//! Normal force and circular radius remain caller-normal-owner inputs; nominal
//! source pressure is never inferred from them. No new solver, history, elastic
//! storage, heat partition or material authority is introduced.

use core::fmt;

use fs_blake3::ContentHash;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_tribo::FrictionLaw;
use fs_tribo::sliding_patch::contact::{
    HertzWrenchContact, HertzWrenchError, HertzWrenchParameters,
};

use super::{SourceThermalCoulomb, SourceThermalCoulombError, coefficients};

/// Immutable source query and its original native Hertz law for one trial.
/// This is a coefficient/geometry binding, not an accepted thermal checkpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHertzWrench {
    owner: ContentHash,
    source: ResolvedInterfaceStatePoint,
    law: HertzWrenchContact,
}

impl SourceHertzWrench {
    /// Original thermal-friction identity, including ordered source, selection,
    /// fixed query coordinates, coefficient roles, axis mapping and authority.
    #[must_use]
    pub const fn owner_identity(&self) -> ContentHash {
        self.owner
    }

    /// Exact temperature/pressure query, selected claims, uncertainty and usage
    /// receipts. The original interface owner can verify every receipt.
    #[must_use]
    pub const fn friction_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }

    /// Original full translational/torsional law. Its candidate solve and all
    /// mechanical/work/cancellation gates remain unchanged.
    #[must_use]
    pub const fn law(&self) -> &HertzWrenchContact {
        &self.law
    }
}

/// Original source or native Hertz refusal, without any published state.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceHertzWrenchError {
    /// Pinned claim, physical query, unit, applicability or coefficient refusal.
    Source(SourceThermalCoulombError),
    /// Original native admission of normal geometry and Coulomb coefficients.
    Hertz(HertzWrenchError),
}

impl From<SourceThermalCoulombError> for SourceHertzWrenchError {
    fn from(error: SourceThermalCoulombError) -> Self {
        Self::Source(error)
    }
}
impl From<HertzWrenchError> for SourceHertzWrenchError {
    fn from(error: HertzWrenchError) -> Self {
        Self::Hertz(error)
    }
}
impl fmt::Display for SourceHertzWrenchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Source(error) => write!(f, "source-backed Hertz coefficients: {error}"),
            Self::Hertz(error) => write!(f, "source-backed Hertz normal geometry: {error}"),
        }
    }
}
impl std::error::Error for SourceHertzWrenchError {}

impl SourceThermalCoulomb {
    /// Resolve the same selected coefficient claims at actual thermal conditions
    /// and bind them to the original coupled Hertz wrench law.
    ///
    /// The supplied compressive normal force [N] and positive circular radius
    /// [m] must come from the caller's admitted normal patch. Native admission
    /// checks their finite physical domain; this bridge does not prove a Hertz
    /// geometry, derive normal elasticity, or promote that geometry's authority.
    /// The optional nominal pressure [Pa] must exactly match the originally
    /// declared source-axis ownership, even at zero normal load. It is not the
    /// pressure implicitly obtained from `N/(pi*a²)`. Temperature is absolute [K].
    /// Sources with additional evolving coordinates refuse this T/p-only entry;
    /// surface temperatures, speed and path cannot be frozen implicitly.
    ///
    /// The source owner is only borrowed. Every trial preserves its pinned claim
    /// identities and fixed external conditions; all returned receipts and the
    /// native law belong to that exact query. An enclosing thermal/mechanical
    /// owner must accept the solved temperature, mechanics and receipts together.
    pub fn hertz_at(
        &self,
        temperature_k: f64,
        nominal_pressure_pa: Option<f64>,
        normal_force_n: f64,
        contact_radius_m: f64,
    ) -> Result<SourceHertzWrench, SourceHertzWrenchError> {
        let source = self.source_at_thermal_conditions(
            &self.initial_source,
            temperature_k,
            nominal_pressure_pa,
        )?;
        let coefficients = coefficients(&source, &self.static_property, &self.kinetic_property)?;
        let FrictionLaw::Coulomb {
            static_mu,
            kinetic_mu,
        } = coefficients
        else {
            return Err(SourceThermalCoulombError::Invalid(
                "the original thermal owner must supply Coulomb coefficients",
            )
            .into());
        };
        let law = HertzWrenchContact::new(HertzWrenchParameters {
            normal_force_n,
            contact_radius_m,
            static_mu,
            kinetic_mu,
        })?;
        Ok(SourceHertzWrench {
            owner: self.identity,
            source,
            law,
        })
    }
}
