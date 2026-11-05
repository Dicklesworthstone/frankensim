//! Original equilibrium derivatives needed by pressure/velocity/temperature
//! acoustics. Cp/Cv fixes the magnitude of thermal coupling, but cannot supply
//! its sign or the physical temperature/heat ports. Those come from the SAME
//! source's thermal state. See the general-fluid constitutive definitions in
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.aco/aco_ug_thermo.09.04.html>.

use super::{ClassicalBulkAbsorption, invalid};
use crate::MaterialError;
use crate::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use crate::fluid_transport::{ThermalFluidTransport, ThermalFluidTransportState};

/// Frozen, thermodynamically coherent Newtonian medium for both acoustic and
/// thermal diffusion modes. Source receipts remain on the original material;
/// this numerical snapshot adds no equation of state or calibration support.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkThermoviscousEquilibrium {
    absorption: ClassicalBulkAbsorption,
    temperature_k: f64,
    expansion_per_k: f64,
    isothermal_compressibility_pa_inv: f64,
}

impl BulkThermoviscousEquilibrium {
    /// Read the original acoustic/caloric and thermal views once each. Unknown
    /// signed expansion refuses; it is never guessed from Cp/Cv. Bulk
    /// viscosity [Pa s] remains explicit and independent of shear viscosity.
    pub fn from_material<M: ThermoviscousFluid + ThermalFluidTransport + ?Sized>(
        material: &M,
        bulk_viscosity_pa_s: f64,
    ) -> Result<Self, MaterialError> {
        Self::from_states(
            material.thermoviscous_state()?,
            material.thermal_transport_state()?,
            bulk_viscosity_pa_s,
        )
    }

    /// Bind explicitly supplied checked views of ONE equilibrium. This also
    /// supports authored numerical fixtures without claiming material data.
    /// Requires identical transport and
    /// `beta_S = beta_T - alpha_p² T0/(rho Cp)`, with
    /// `beta_T = (Cp/Cv)/(rho c²)`, to floating-point accuracy. Negative and
    /// exactly zero expansion are supported. Inconsistent tuples refuse.
    pub fn from_states(
        fluid: ThermoviscousFluidState,
        thermal: ThermalFluidTransportState,
        bulk_viscosity_pa_s: f64,
    ) -> Result<Self, MaterialError> {
        if fluid.transport() != thermal.transport() {
            return Err(invalid(
                "bulk acoustic and thermal views require identical transport",
            ));
        }
        let alpha = thermal.volumetric_expansion_per_k().ok_or_else(|| {
            invalid("full bulk acoustics requires known signed isobaric expansion")
        })?;
        let acoustic = fluid.acoustic();
        let beta_s = acoustic.adiabatic_compressibility_per_pa();
        let beta_t = fluid.heat_capacity_ratio() * beta_s;
        // Same general-fluid storage identity and evaluation order as the
        // source-bound mixed thermoviscous FE owner. No ideal-gas inference.
        let minimum_cp =
            alpha * alpha * thermal.temperature_k() / acoustic.density_kg_m3() / beta_t;
        let actual = beta_t * (1.0 - minimum_cp / fluid.transport().specific_heat_j_kg_k());
        if !beta_t.is_finite()
            || beta_t <= 0.0
            || !minimum_cp.is_finite()
            || !actual.is_finite()
            || actual <= 0.0
            || (actual - beta_s).abs() > 256.0 * f64::EPSILON * beta_s
        {
            return Err(invalid(
                "bulk Cp/Cv, sound speed and expansion violate one thermodynamic closure",
            ));
        }
        Ok(Self {
            absorption: ClassicalBulkAbsorption::from_material(&fluid, bulk_viscosity_pa_s)?,
            temperature_k: thermal.temperature_k(),
            expansion_per_k: alpha,
            isothermal_compressibility_pa_inv: beta_t,
        })
    }

    /// The same checked medium, including the independently declared bulk viscosity.
    pub const fn absorption(&self) -> &ClassicalBulkAbsorption {
        &self.absorption
    }
    /// Original absolute reference temperature [K].
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    /// Original signed isobaric volumetric expansion [1/K].
    pub const fn isobaric_expansion_per_k(&self) -> f64 {
        self.expansion_per_k
    }
    /// Isothermal compressibility [1/Pa], from the same caloric/acoustic state.
    pub const fn isothermal_compressibility_pa_inv(&self) -> f64 {
        self.isothermal_compressibility_pa_inv
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid_acoustics::FluidAcousticState;
    use crate::fluid_transport::FluidTransportState;

    fn states(alpha: Option<f64>) -> (ThermoviscousFluidState, ThermalFluidTransportState) {
        let transport = FluidTransportState::try_new(2.0, 8.0, 0.6, 0.003).unwrap();
        let beta_s = 1.0 / (2.0 * 100.0);
        let beta_t = beta_s + alpha.unwrap_or(0.08).powi(2) * 12.5 / (2.0 * 8.0);
        (
            ThermoviscousFluidState::try_new(
                FluidAcousticState::try_new(2.0, 10.0).unwrap(),
                transport,
                8.0 * beta_s / beta_t,
            )
            .unwrap(),
            ThermalFluidTransportState::try_new(transport, 12.5, alpha).unwrap(),
        )
    }

    #[test]
    fn equilibrium_keeps_signed_expansion_and_actual_transport() {
        for alpha in [-0.08, 0.0, 0.08] {
            let (fluid, thermal) = states(Some(alpha));
            let bound = BulkThermoviscousEquilibrium::from_states(fluid, thermal, 0.002).unwrap();
            assert_eq!(bound.absorption().state(), fluid);
            assert_eq!(bound.absorption().bulk_viscosity_pa_s(), 0.002);
            assert_eq!(bound.temperature_k(), 12.5);
            assert_eq!(bound.isobaric_expansion_per_k(), alpha);
            let beta_t = bound.isothermal_compressibility_pa_inv();
            let beta_s = beta_t - alpha * alpha * 12.5 / (2.0 * 8.0);
            assert!((beta_s * 200.0 - 1.0).abs() < 2e-15);
        }
    }

    #[test]
    fn equilibrium_refuses_unknown_expansion_and_inconsistent_views() {
        let (fluid, unknown) = states(None);
        assert!(BulkThermoviscousEquilibrium::from_states(fluid, unknown, 0.0).is_err());
        for (temperature, alpha) in [(25.0, 0.08), (12.5, 0.04), (12.5, 0.0)] {
            let wrong =
                ThermalFluidTransportState::try_new(fluid.transport(), temperature, Some(alpha))
                    .unwrap();
            assert!(BulkThermoviscousEquilibrium::from_states(fluid, wrong, 0.0).is_err());
        }
        let wrong_transport = FluidTransportState::try_new(2.0, 8.0, 0.6, 0.006).unwrap();
        let wrong = ThermalFluidTransportState::try_new(wrong_transport, 12.5, Some(0.08)).unwrap();
        assert!(BulkThermoviscousEquilibrium::from_states(fluid, wrong, 0.0).is_err());
        let (_, thermal) = states(Some(0.08));
        for bulk in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(BulkThermoviscousEquilibrium::from_states(fluid, thermal, bulk).is_err());
        }
    }
}
