//! Original isobaric gas properties on a native conducting exterior surface.

use super::*;
use fs_material::gas_species::input::GasRecord;
use fs_material::gas_species::input::temperature::{
    IsobaricGasLaw, IsobaricGasTransportTangent, PreparedIsobaricGasTransport,
};

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
/// Original temperature/pressure surfaces prepare their cell geometry once and
/// retain fixed-pressure, inward thermal tangents at admitted source endpoints.
/// Bulk-held laws without a selected wall correction retain exactly zero
/// wall-property derivatives and need no active transport tangent.
/// [`Self::with_gnielinski_air_wall_correction`] explicitly selects the
/// original air temperature-ratio correction for developing turbulent tubes.
/// Its actual wall-temperature tangent needs no fitted transport derivative.
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
    transport: Option<PreparedIsobaricGasTransport<'g>>,
    settings: ExteriorConvectionSettings,
    gnielinski_air_wall_correction: bool,
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
        let transport = if settings.property_temperature == ConvectionPropertyTemperature::FilmMean
            || settings.law.needs_surface_viscosity()
        {
            // Ask the original source for analytic derivatives, including at
            // the admitted interval endpoints. No secant or fitted transport
            // model stands in for an unavailable thermochemical tangent.
            let transport = gas.prepare_transport_temperature_tangent()?;
            for temperature in [lower, upper] {
                transport.evaluate(Temperature::new(temperature))?;
            }
            Some(transport)
        } else {
            None
        };
        Ok(Self {
            caloric,
            energy_bounds_j_kg,
            gas,
            transport,
            settings,
            gnielinski_air_wall_correction: false,
        })
    }

    /// Explicitly select Gnielinski's air correction `(Tbulk/Twall)^0.45`
    /// for this developing turbulent circular-tube section.
    ///
    /// Only the retained original dry-air source is admitted. Humid air,
    /// arbitrary named mixtures and thermochemical gases are not silently
    /// declared eligible for the air exponent. Bulk Re, Pr and conductivity
    /// still come from the same source at the actual bath temperature.
    ///
    /// Both temperatures are absolute kelvin. The native surface solve supplies
    /// the actual exterior temperature after material/coating resistance;
    /// its analytic partial is `-0.45 h/Twall`. The coefficient also retains
    /// the explicit `T_bulk` partial at fixed transport groups, not a total
    /// bath-temperature derivative. Prescribed mass flow retains the correction.
    /// Original source intervals and the correlation's ratio admission remain
    /// active. No air-property substitution or second thermal solve is added.
    pub fn with_gnielinski_air_wall_correction(mut self) -> Result<Self, ExteriorConvectionError> {
        if !matches!(
            self.settings.law,
            SourceCellConvectionLaw::ForcedTubeSection { .. }
        ) || self.settings.property_temperature != ConvectionPropertyTemperature::Reservoir
        {
            return Err(ExteriorConvectionError::Invalid(
                "air Gnielinski wall correction requires a bulk-property turbulent circular-tube section",
            ));
        }
        if !matches!(self.gas.reference(), GasRecord::Legacy(state)
            if state.water_mole_fraction == 0.0)
        {
            return Err(ExteriorConvectionError::Invalid(
                "air Gnielinski wall correction requires the retained original dry-air source",
            ));
        }
        let [lo, hi] = self.settings.fluid_temperature_bounds_k;
        let [wall_lo, wall_hi] = self.settings.surface_temperature_bounds_k;
        if wall_lo < lo || wall_hi > hi {
            return Err(ExteriorConvectionError::Invalid(
                "air wall correction requires original gas support over the whole surface interval",
            ));
        }
        self.gnielinski_air_wall_correction = true;
        Ok(self)
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

    fn transport_tangent(
        &self,
        temperature: Temperature,
    ) -> Result<IsobaricGasTransportTangent, ExteriorConvectionError> {
        Ok(if let Some(prepared) = &self.transport {
            prepared.evaluate(temperature)?
        } else {
            self.gas.transport_temperature_tangent(temperature)?
        })
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

    fn turbulent_wall_correction(
        &self,
        temperature: f64,
    ) -> Result<Option<evaluator::TurbulentWallCorrection>, ExteriorConvectionError> {
        if self.gnielinski_air_wall_correction {
            self.admit_temperature(temperature)?;
            Ok(Some(evaluator::TurbulentWallCorrection::Air))
        } else {
            Ok(None)
        }
    }

    fn sample_transport(
        &self,
        temperature: f64,
        derivatives: bool,
    ) -> Result<TransportSample, ExteriorConvectionError> {
        self.admit_temperature(temperature)?;
        let temperature = Temperature::new(temperature);
        if derivatives {
            let sample = self.transport_tangent(temperature)?;
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
        let sample = self.transport_tangent(Temperature::new(temperature))?;
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
