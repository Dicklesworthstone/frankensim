//! Quiescent-reservoir convection with source-owned thermal expansion.

use super::*;
use fs_material::fluid_transport::ThermalFluidTransport;

/// Exposed patch of one explicitly oriented isothermal body in quiescent gas.
/// The reservoir is infinite and prescribed; heat entering it is external heat
/// leaving the finite body. No free stream or mixed-convection law is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NaturalGasBoundary {
    pub geometry: NaturalConvectionGeometry,
    /// Plate height, cylinder/sphere diameter, or FULL horizontal planform
    /// area/perimeter. A contact-masked exposed patch cannot determine A/P.
    pub characteristic_length_m: f64,
    pub exposed_area_m2: f64,
    /// Nonnegative magnitude. Horizontal face orientation is explicit above.
    pub gravity_m_s2: f64,
}

impl NaturalGasBoundary {
    pub fn validate(self) -> Result<(), GasConvectionError> {
        if [self.characteristic_length_m, self.exposed_area_m2]
            .into_iter().any(|x| !x.is_finite() || x <= 0.0)
            || !self.gravity_m_s2.is_finite() || self.gravity_m_s2 < 0.0
        {
            return Err(GasConvectionError::Invalid(
                "natural convection needs positive finite length/area and nonnegative gravity",
            ));
        }
        if self.geometry == NaturalConvectionGeometry::Sphere
            && self.exposed_area_m2 > core::f64::consts::PI
                * self.characteristic_length_m * self.characteristic_length_m
        {
            return Err(GasConvectionError::Invalid(
                "sphere exposed area cannot exceed pi times diameter squared",
            ));
        }
        Ok(())
    }

    /// Re-query the original gas at mean film temperature and retained bulk
    /// pressure/composition. Its thermal transport owner supplies beta; the
    /// signed full wall-to-bulk contrast chooses horizontal-face stability.
    /// Existing source, phase, Rayleigh and Prandtl restrictions remain active.
    ///
    /// Exact thermal equilibrium has exactly zero heat and entropy production.
    /// Except for the sphere's sourced Nu=2 limit, this does not evaluate or
    /// extrapolate an out-of-range Ra=0 correlation. A nonzero temperature
    /// difference below a card's lower bound still refuses.
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
            GasConvectionError::Invalid("natural convection requires source-owned thermal expansion"),
        )?;
        let difference = wall_temperature_k - reservoir_temperature_k;
        let nusselt = if difference == 0.0
            && self.geometry != NaturalConvectionGeometry::Sphere
        {
            None
        } else {
            let rayleigh = transport.rayleigh(self.gravity_m_s2, expansion,
                difference, self.characteristic_length_m)?;
            Some(self.geometry.evaluate(rayleigh, transport.prandtl())?)
        };
        finish_sample(property_gas, None, nusselt,
            transport.thermal_conductivity_w_m_k(), self.characteristic_length_m,
            self.exposed_area_m2, wall_temperature_k, reservoir_temperature_k)
    }
}
