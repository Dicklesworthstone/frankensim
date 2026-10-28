//! Gas substitution in the existing rigid-frame JCA(L) material law.
//!
//! Pore geometry and intrinsic permeability belong to the solid skeleton.
//! Static flow resistivity does not: `sigma = mu / k0` must follow the actual
//! saturating gas. The characterization viscosity is explicit source data;
//! neither the card's name nor a room-air default supplies it.
//!
//! The existing [`PorousCard::melamine_uf_nguyen2024`] owner records this same
//! relation from Nguyen et al., Acta Acustica 8:54 (2024),
//! doi:10.1051/aacus/2024046, Table 1. All dynamic-density, thermal-relaxation
//! and rigid-backing formulas remain solely in [`PorousCard`].

use super::{PorousAmbient, PorousCard};
use super::permeability::IntrinsicPermeability;
use crate::{MaterialError, gas::GasState};
use fs_math::c64::C64;

/// An unchanged porous skeleton with explicitly characterized gas drag.
///
/// This retains the original card and derives intrinsic viscous permeability
/// from its measured/identified resistivity and the supplied characterization
/// viscosity. Each query re-evaluates gas drag, inertia and thermal relaxation
/// from one coherent [`GasState`], through the original JCA(L) owner.
///
/// The frame remains rigid and the pores open and gas-filled. Source support
/// for that unchanged skeleton, continuum flow and gas phase remains the
/// caller's responsibility; gas substitution does not qualify swelling,
/// liquid saturation, condensation, rarefaction or elastic-frame waves.
#[derive(Debug, Clone)]
pub struct PorousGasLaw {
    card: PorousCard,
    reference_dynamic_viscosity_pa_s: f64,
    viscous_permeability: IntrinsicPermeability,
}

impl PorousAmbient {
    /// Map the canonical coherent gas snapshot without independent coefficients.
    ///
    /// # Errors
    /// Propagates [`GasState::validate_coherence`] failures.
    pub fn from_gas(gas: &GasState) -> Result<Self, MaterialError> {
        gas.validate_coherence()?;
        Ok(Self {
            rho0: gas.density,
            c0: gas.sound_speed,
            eta: gas.dynamic_viscosity,
            gamma: gas.gamma,
            p0: gas.pressure,
            prandtl: gas.prandtl,
        })
    }
}

impl PorousGasLaw {
    /// Bind the exact viscosity used to characterize this card's resistivity.
    /// No characterization temperature or gas identity is guessed.
    ///
    /// # Errors
    /// Refuses invalid cards, absent provenance, nonfinite support, nonpositive
    /// viscosity, and unrepresentable intrinsic permeability.
    pub fn try_new(
        card: PorousCard,
        reference_dynamic_viscosity_pa_s: f64,
    ) -> Result<Self, MaterialError> {
        card.validate()?;
        let permeability = reference_dynamic_viscosity_pa_s / card.flow_resistivity;
        if card.provenance.trim().is_empty()
            || !card.valid_hz.1.is_finite()
            || !reference_dynamic_viscosity_pa_s.is_finite()
            || reference_dynamic_viscosity_pa_s <= 0.0
            || !permeability.is_finite()
            || permeability <= 0.0
        {
            return Err(MaterialError::Parameters {
                what: "porous gas law needs source provenance, finite support and positive characterization viscosity/permeability".into(),
            });
        }
        Ok(Self {
            card,
            reference_dynamic_viscosity_pa_s,
            viscous_permeability: IntrinsicPermeability::isotropic(permeability)?,
        })
    }

    /// Original characterization, including its source and frequency support.
    #[must_use]
    pub const fn reference_card(&self) -> &PorousCard {
        &self.card
    }

    /// Supplied characterization-gas viscosity [Pa s].
    #[must_use]
    pub const fn reference_dynamic_viscosity_pa_s(&self) -> f64 {
        self.reference_dynamic_viscosity_pa_s
    }

    /// Intrinsic viscous permeability [m2], unchanged by gas substitution.
    #[must_use]
    pub const fn viscous_permeability_m2(&self) -> f64 {
        self.viscous_permeability.principal_m2()[0]
    }

    /// The same intrinsic pore geometry used by continuum Darcy consumers.
    /// Acoustic frequency support alone does not qualify a steady-flow regime.
    #[must_use]
    pub const fn intrinsic_permeability(&self) -> IntrinsicPermeability {
        self.viscous_permeability
    }

    /// Surface impedance [Pa s/m] of the gas-filled, rigid-backed layer.
    ///
    /// Uses [`PorousCard::surface_impedance_rigid`]'s `exp(+i omega t)`
    /// convention. A consumer with the opposite convention must conjugate;
    /// a volume-flow load must also divide by the explicitly exposed area.
    /// No ambient value is cached.
    ///
    /// # Errors
    /// Refuses an incoherent gas, frequency outside the original source band,
    /// invalid thickness, and nonfinite or nonpassive resulting impedance.
    pub fn surface_impedance_rigid(
        &self,
        gas: &GasState,
        omega: f64,
        thickness_m: f64,
    ) -> Result<C64, MaterialError> {
        let ambient = PorousAmbient::from_gas(gas)?;
        let band = self.card.valid_hz;
        if !omega.is_finite()
            || omega < core::f64::consts::TAU * band.0
            || omega > core::f64::consts::TAU * band.1
        {
            return Err(MaterialError::State {
                what: "porous gas query is outside the card's frequency support".into(),
            });
        }
        let mut current = self.card.clone();
        if gas.dynamic_viscosity != self.reference_dynamic_viscosity_pa_s {
            current.flow_resistivity = gas.dynamic_viscosity / self.viscous_permeability_m2();
        }
        let impedance = current.surface_impedance_rigid(ambient, omega, thickness_m)?;
        if !impedance.re.is_finite() || !impedance.im.is_finite() || impedance.re < 0.0 {
            return Err(MaterialError::State {
                what: "gas-filled porous surface impedance must be finite and passive".into(),
            });
        }
        Ok(impedance)
    }
}

#[cfg(test)]
mod tests;
