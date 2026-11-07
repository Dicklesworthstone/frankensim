//! Original isobaric gas properties on a native conducting exterior surface.

use super::*;
use fs_material::gas_species::input::temperature::IsobaricGasLaw;

/// Gas convection on `IntegratedCaloricGrayBoundary`, including the actual
/// temperature drop across the body's original half-cell/coating resistance.
///
/// Pressure, composition, source choices and phase admission remain with the
/// retained `IsobaricGasLaw`. Mean-film laws refresh density, Cp, conductivity
/// and viscosity together at `(Ts + Tbulk)/2`. Their analytic source tangents
/// enter the SAME correlation chain rule used by the liquid exterior owner.
/// Natural convection uses the original ideal-gas isobaric expansion `1/T`;
/// no liquid expansion or dry-air replacement is introduced.
///
/// Sutherland/USSA/Eucken, named mixtures and supported original thermochemical
/// sources supply analytic transport tangents. Thermochemical NASA calorics,
/// selected transport claims and Wilke/WMS denominator derivatives stay with
/// their source owners; unsupported derivatives and genuine source kinks refuse.
/// Bulk-held laws without a wall-viscosity correction retain exactly zero
/// wall-property derivatives and need no active transport tangent.
///
/// Gas and liquid sources share physical geometry, equilibrium heat limits,
/// correlation-domain admission and explicit property-temperature conventions.
/// Rectangular four-wall sections remain unsupported by a pointwise exterior
/// interface. The native owner retains solid storage, surface elimination,
/// radiation and energy/entropy accounting; this is not another heat stepper.
pub struct ImplicitGasConvection<'a, 'g> {
    caloric: &'a dyn ThermalLaw,
    energy_bounds_j_kg: [f64; 2],
    gas: &'g IsobaricGasLaw,
    settings: ExteriorConvectionSettings,
}

impl<'a, 'g> ImplicitGasConvection<'a, 'g> {
    /// Bind the SAME solid caloric object used by the native conducting boundary
    /// and the original gas law. The declared fluid interval is admitted before
    /// constructing the boundary; each actual film/bulk/wall query rechecks it.
    pub fn try_new(
        caloric: &'a dyn ThermalLaw,
        energy_bounds_j_kg: [f64; 2],
        gas: &'g IsobaricGasLaw,
        settings: ExteriorConvectionSettings,
    ) -> Result<Self, ExteriorConvectionError> {
        evaluator::admit(caloric, energy_bounds_j_kg, settings)?;
        let [lower, upper] = settings.fluid_temperature_bounds_k;
        gas.admit_caloric_interval(Temperature::new(lower), Temperature::new(upper))?;
        if settings.property_temperature == ConvectionPropertyTemperature::FilmMean
            || settings.law.needs_surface_viscosity()
        {
            // Ask the original source for analytic derivatives, including at
            // the admitted interval endpoints. No secant or fitted transport
            // model stands in for an unavailable thermochemical tangent.
            for temperature in [lower, upper] {
                gas.transport_temperature_tangent(Temperature::new(temperature))?;
            }
        }
        Ok(Self {
            caloric,
            energy_bounds_j_kg,
            gas,
            settings,
        })
    }

    /// Retained original gas, including its reference pressure, composition and
    /// source authority. Its reference is never mutated by a surface trial.
    pub const fn gas(&self) -> &IsobaricGasLaw {
        self.gas
    }

    /// Physical surface/bath evaluation with original source and correlation
    /// failures intact. Material uncertainty is not promoted to a tangent bound.
    pub fn evaluate(
        &self,
        surface_temperature_k: f64,
        bath_temperature_k: f64,
    ) -> Result<ExteriorConvectionResponse, ExteriorConvectionError> {
        evaluator::evaluate(self, surface_temperature_k, bath_temperature_k)
    }

    fn admit_temperature(&self, temperature: f64) -> Result<(), ExteriorConvectionError> {
        let [lower, upper] = self.settings.fluid_temperature_bounds_k;
        if !temperature.is_finite() || !(lower..=upper).contains(&temperature) {
            return Err(ExteriorConvectionError::Invalid(
                "gas property temperature leaves its original declared interval",
            ));
        }
        Ok(())
    }
}

impl evaluator::ExteriorFluid for ImplicitGasConvection<'_, '_> {
    fn caloric(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn energy_bounds(&self) -> [f64; 2] {
        self.energy_bounds_j_kg
    }
    fn settings(&self) -> ExteriorConvectionSettings {
        self.settings
    }

    fn sample_transport(
        &self,
        temperature: f64,
        derivatives: bool,
    ) -> Result<TransportSample, ExteriorConvectionError> {
        self.admit_temperature(temperature)?;
        let temperature = Temperature::new(temperature);
        if derivatives {
            let sample = self.gas.transport_temperature_tangent(temperature)?;
            Ok(TransportSample {
                state: sample.gas.transport_state()?,
                derivatives: [
                    sample.density_derivative_kg_m3_k,
                    sample.specific_heat_derivative_j_kg_k2,
                    sample.conductivity_derivative_w_m_k2,
                    sample.dynamic_viscosity_derivative_pa_s_k,
                ],
            })
        } else {
            Ok(TransportSample {
                state: self.gas.at_temperature(temperature)?.transport_state()?,
                derivatives: [0.0; 4],
            })
        }
    }

    fn wall_viscosity(
        &self,
        temperature: f64,
    ) -> Result<(f64, f64, bool), ExteriorConvectionError> {
        self.admit_temperature(temperature)?;
        let sample = self
            .gas
            .transport_temperature_tangent(Temperature::new(temperature))?;
        let slope = sample.dynamic_viscosity_derivative_pa_s_k;
        Ok((sample.gas.state().dynamic_viscosity, slope, slope != 0.0))
    }

    fn expansion(
        &self,
        temperature: f64,
        _sample: &TransportSample,
        derivative: bool,
    ) -> Result<(f64, f64), ExteriorConvectionError> {
        self.admit_temperature(temperature)?;
        // Every retained IsobaricGasLaw owner has the same frozen-composition
        // ideal-gas EOS rho=p/(R T). Its isobaric beta is 1/T, independently of
        // the transport fit or NASA heat-capacity variation.
        let beta = finite(1.0 / temperature)?;
        Ok((
            beta,
            if derivative {
                finite(-beta * beta)?
            } else {
                0.0
            },
        ))
    }
}

impl ImplicitGrayExteriorLaw for ImplicitGasConvection<'_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn energy_bounds_j_kg(&self) -> [f64; 2] {
        self.energy_bounds_j_kg
    }
    fn surface_temperature_bounds_k(&self) -> [f64; 2] {
        self.settings.surface_temperature_bounds_k
    }
    fn flux_state_at_surface(
        &self,
        energy: f64,
        surface: f64,
        bath: f64,
    ) -> Result<ImplicitGrayExteriorFluxState, DepositionError> {
        evaluator::flux_state_at_surface(self, energy, surface, bath)
    }
    fn state_at_surface(
        &self,
        energy: f64,
        surface: f64,
        bath: f64,
    ) -> Result<ImplicitGrayExteriorState, DepositionError> {
        evaluator::state_at_surface(self, energy, surface, bath)
    }
}
