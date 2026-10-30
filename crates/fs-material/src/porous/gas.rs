//! Fluid substitution in the existing rigid-frame JCA(L) material law.
//!
//! Pore geometry and intrinsic permeability belong to the solid skeleton.
//! Static flow resistivity does not: `sigma = mu / k0` must follow the actual
//! saturating fluid. The characterization viscosity is explicit source data;
//! neither the card's name nor a room-air default supplies it.
//!
//! The existing [`PorousCard::melamine_uf_nguyen2024`] owner records this same
//! relation from Nguyen et al., Acta Acustica 8:54 (2024),
//! doi:10.1051/aacus/2024046, Table 1. All dynamic-density, thermal-relaxation
//! and rigid-backing formulas remain solely in [`PorousCard`].
//!
//! General fluids replace ideal-gas gamma*p by their original isentropic
//! modulus K_s, while Cp/Cv controls the thermal relaxation. The general-fluid
//! derivation is given in COMSOL's Acoustics Module theory, Eqs. 2-46--2-56:
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.182.html>.
//! This is a constitutive model transfer, not experimental qualification of an
//! air-characterized skeleton in water or of a freely moving porous frame.

use super::permeability::IntrinsicPermeability;
use super::{PorousAmbient, PorousCard};
use super::{characteristic_from_response, rigid_surface_from_characteristic};
use crate::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use crate::{MaterialError, gas::GasState};
use fs_math::c64::C64;

/// An unchanged porous skeleton with explicitly characterized viscous drag.
///
/// This retains the original card and derives intrinsic viscous permeability
/// from its measured/identified resistivity and the supplied characterization
/// viscosity. Each query re-evaluates drag, inertia and thermal relaxation
/// from one coherent fluid state, through the original JCA(L) owner.
///
/// The frame remains rigid and the open pores are fully filled by ONE fluid.
/// Source support for that unchanged skeleton, continuum flow, frequency and
/// phase remains the caller's responsibility. Substitution does not qualify
/// swelling, partial saturation, condensation, rarefaction or frame waves.
#[derive(Debug, Clone)]
pub struct PorousFluidLaw {
    card: PorousCard,
    reference_dynamic_viscosity_pa_s: f64,
    viscous_permeability: IntrinsicPermeability,
}

/// Existing gas-facing spelling of the same source and numerical owner.
pub type PorousGasLaw = PorousFluidLaw;

/// One evaluated rigid-frame response retaining the actual coherent fluid.
/// The original material and pore calibration retain source authority; this
/// numerical sample neither selects new claims nor extends their domains.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PorousFluidResponse {
    fluid: ThermoviscousFluidState,
    omega: f64,
    flow_resistivity: f64,
    density: C64,
    bulk: C64,
    impedance: C64,
    wavenumber: C64,
}

impl PorousFluidResponse {
    /// Actual acoustic, caloric and transport snapshot used together.
    pub const fn fluid(&self) -> &ThermoviscousFluidState {
        &self.fluid
    }
    /// Angular frequency [rad/s].
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    /// Current static resistivity mu/k0 [Pa s/m²], not a new measurement.
    pub const fn flow_resistivity_pa_s_m2(&self) -> f64 {
        self.flow_resistivity
    }
    /// Effective dynamic density [kg/m³], including the porosity convention.
    pub const fn dynamic_density_kg_m3(&self) -> C64 {
        self.density
    }
    /// Effective dynamic bulk modulus [Pa], including the porosity convention.
    pub const fn dynamic_bulk_modulus_pa(&self) -> C64 {
        self.bulk
    }
    /// Specific characteristic impedance [Pa s/m] and wavenumber [rad/m],
    /// with the original exp(+i omega t) convention.
    pub const fn characteristic(&self) -> (C64, C64) {
        (self.impedance, self.wavenumber)
    }
    /// Rigid-backed surface impedance [Pa s/m] using the same sampled fluid.
    pub fn surface_impedance_rigid(&self, thickness_m: f64) -> Result<C64, MaterialError> {
        passive_surface(rigid_surface_from_characteristic(
            self.impedance,
            self.wavenumber,
            thickness_m,
        )?)
    }
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

fn passive_surface(impedance: C64) -> Result<C64, MaterialError> {
    if !finite(impedance) || impedance.re < 0.0 {
        return Err(MaterialError::State {
            what: "fluid-filled porous surface impedance must be finite and passive".into(),
        });
    }
    Ok(impedance)
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

impl PorousFluidLaw {
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

    /// Supplied characterization-fluid viscosity [Pa s].
    #[must_use]
    pub const fn reference_dynamic_viscosity_pa_s(&self) -> f64 {
        self.reference_dynamic_viscosity_pa_s
    }

    /// Intrinsic viscous permeability [m2], unchanged by fluid substitution.
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

    fn admit_frequency(&self, omega: f64) -> Result<(), MaterialError> {
        let band = self.card.valid_hz;
        if !omega.is_finite()
            || omega < core::f64::consts::TAU * band.0
            || omega > core::f64::consts::TAU * band.1
        {
            return Err(MaterialError::State {
                what: "porous fluid query is outside the card's frequency support".into(),
            });
        }
        Ok(())
    }

    fn current_card(&self, viscosity: f64) -> PorousCard {
        let mut current = self.card.clone();
        if viscosity != self.reference_dynamic_viscosity_pa_s {
            current.flow_resistivity = viscosity / self.viscous_permeability_m2();
        }
        current
    }

    /// Re-evaluate this same skeleton in an original gas, liquid or nonideal
    /// fluid owner. Density, adiabatic modulus, viscosity, Cp/Cv and thermal
    /// diffusion come from ONE [`ThermoviscousFluid`] call. No pressure or
    /// ideal-gas closure is invented for a liquid. Cp=Cv gives exactly K_s/phi
    /// at every frequency and removes thermal relaxation, but viscous drag
    /// remains active. Missing transport or Cv refuses at the original owner.
    pub fn effective_fluid<M: ThermoviscousFluid + ?Sized>(
        &self,
        fluid: &M,
        omega: f64,
    ) -> Result<PorousFluidResponse, MaterialError> {
        self.admit_frequency(omega)?;
        let fluid = fluid.thermoviscous_state()?;
        let acoustic = fluid.acoustic();
        let transport = fluid.transport();
        let current = self.current_card(transport.dynamic_viscosity_pa_s());
        let density = current.effective_density_from(
            acoustic.density_kg_m3(),
            transport.dynamic_viscosity_pa_s(),
            omega,
        )?;
        let bulk = current.effective_bulk_from(
            acoustic.density_kg_m3(),
            transport.dynamic_viscosity_pa_s(),
            transport.prandtl(),
            fluid.heat_capacity_ratio(),
            acoustic.adiabatic_bulk_modulus_pa(),
            omega,
        )?;
        let (impedance, wavenumber) = characteristic_from_response(density, bulk, omega);
        if !finite(density)
            || density.re <= 0.0
            || density.im > 0.0
            || !finite(bulk)
            || bulk.re <= 0.0
            || bulk.im < 0.0
            || !finite(impedance)
            || impedance.re <= 0.0
            || !finite(wavenumber)
            || wavenumber.re <= 0.0
            || wavenumber.im > 0.0
        {
            return Err(MaterialError::State {
                what: "fluid-filled porous response must be finite and passive".into(),
            });
        }
        Ok(PorousFluidResponse {
            fluid,
            omega,
            flow_resistivity: current.flow_resistivity,
            density,
            bulk,
            impedance,
            wavenumber,
        })
    }

    /// Rigid-backed surface impedance [Pa s/m] from the actual saturating
    /// fluid. Phasor and exposed-area conversions belong to the consumer.
    /// The same bounded JCA(L) owner supplies gas and liquid results.
    pub fn surface_impedance_rigid_fluid<M: ThermoviscousFluid + ?Sized>(
        &self,
        fluid: &M,
        omega: f64,
        thickness_m: f64,
    ) -> Result<C64, MaterialError> {
        self.effective_fluid(fluid, omega)?
            .surface_impedance_rigid(thickness_m)
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
        self.admit_frequency(omega)?;
        let current = self.current_card(gas.dynamic_viscosity);
        passive_surface(current.surface_impedance_rigid(ambient, omega, thickness_m)?)
    }
}

#[cfg(test)]
mod tests;
