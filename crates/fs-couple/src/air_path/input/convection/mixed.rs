//! Aiding mixed convection using the original gas and shared group admission.

use super::*;
use fs_material::fluid_transport::ThermalFluidTransport;

/// An unconfined isothermal vertical plate in prescribed laminar aiding flow.
///
/// Signed upward speed, gravity and a caller-admitted Boussinesq bound are
/// explicit. For world-frame vectors, construct settings through
/// [`AidingVerticalPlateFrame::settings`]; oblique flow cannot be replaced by
/// its vertical projection. Exposed area excludes contact-covered patches.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AidingVerticalPlateGasBoundary {
    /// Plate height along the external stream [m], not the exposed patch size.
    pub characteristic_length_m: f64,
    pub exposed_area_m2: f64,
    pub settings: AidingVerticalPlateSettings,
}

impl AidingVerticalPlateGasBoundary {
    pub fn validate(self) -> Result<(), GasConvectionError> {
        if [self.characteristic_length_m, self.exposed_area_m2]
            .into_iter().any(|value| !value.is_finite() || value <= 0.0)
        {
            return Err(GasConvectionError::Invalid(
                "aiding mixed convection needs positive finite plate height and exposed area",
            ));
        }
        self.settings.validate()?;
        Ok(())
    }

    /// Re-query the original source at the arithmetic wall/bulk film mean.
    /// Re, Ra, Pr, expansion and k use that ONE state at retained pressure and
    /// composition; buoyancy and heat use the FULL wall/bulk contrast.
    ///
    /// Manning-Qureshi's original relation assumes constant properties. Film
    /// sampling is this adapter's frozen-property approximation, not a rule
    /// asserted by the paper. The density bound does not establish that other
    /// properties vary weakly. The shared correlation retains all laminar and
    /// Prandtl limits, including the pure forced and natural branches.
    ///
    /// Exact thermal equilibrium at rest gives zero heat without evaluating
    /// the card's invalid Re=Ra=0 point. A resting non-isothermal state with
    /// neither forced nor buoyant driving still refuses; no conduction floor
    /// or separate stagnant-gas path is invented.
    pub fn evaluate(
        self, gas: &IsobaricGasLaw, wall_temperature_k: f64,
    ) -> Result<GasConvectionSample, GasConvectionError> {
        self.validate()?;
        if !wall_temperature_k.is_finite() || wall_temperature_k <= 0.0 {
            return Err(GasConvectionError::Invalid(
                "wall temperature must be positive finite kelvin",
            ));
        }
        let reservoir_temperature_k = gas.reference().state().temperature;
        let property_gas = gas.at_temperature(Temperature::new(
            0.5 * wall_temperature_k + 0.5 * reservoir_temperature_k,
        ))?;
        let thermal = property_gas.thermal_transport_state()?;
        let transport = thermal.transport();
        let expansion = thermal.volumetric_expansion_per_k().ok_or(
            GasConvectionError::Invalid("mixed convection requires source-owned thermal expansion"),
        )?;
        let groups = self.settings.groups(transport, expansion,
            Temperature::new(wall_temperature_k), Temperature::new(reservoir_temperature_k),
            Length::new(self.characteristic_length_m))?;
        let nusselt = if wall_temperature_k == reservoir_temperature_k && groups.reynolds == 0.0 {
            None
        } else {
            Some(fs_convection::evaluate(CorrelationId::ManningQureshiAidingVerticalPlate,
                CorrelationInputs::mixed(groups.reynolds, groups.rayleigh, groups.prandtl))?)
        };
        finish_sample(property_gas, None, nusselt,
            transport.thermal_conductivity_w_m_k(), self.characteristic_length_m,
            self.exposed_area_m2, wall_temperature_k, reservoir_temperature_k)
    }
}
