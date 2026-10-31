//! Original IAPWS water calorics for finite-rate evaporation and condensation.
//!
//! This adapter retains an existing `IsobaricWater` chart: liquid enthalpy,
//! density, inverse temperature and thermal entropy all come from that SAME
//! IF97 source. Region 4 supplies vapor pressure and its analytic derivative;
//! Region 2 supplies saturated-vapor enthalpy. The phase-transfer enthalpy is
//! `h_v(T,p_sat(T))-h_l(T,p_liquid)`, so the compressed-liquid pressure term is
//! retained rather than substituting an unrelated constant latent heat.
//!
//! The finite-rate gas film belongs to `fs_thermal::evaporation`. This bounded
//! engineering adapter admits 273.16..=373.15 K, liquid pressure at most 0.2 MPa,
//! and a strictly sub-boiling liquid interval. Its liquid pressure remains the
//! chart's prescribed mechanical load. The caller must reconcile this load
//! with its gas/interface mechanics. No Poynting or humidity enhancement factor,
//! gas fugacity, supercooling, solute activity, curvature or nucleation model is
//! inferred. Source calorics do not certify the ideal binary-film approximation.
//! Above these limits, use an explicitly admitted nonideal phase/transport law.
//!
//! Source: IAPWS R7-97(2012), unchanged equations and coefficients already owned
//! by `fs_thermochem::water::if97`:
//! <https://iapws.org/technical-guidance/release/IF97-Rev.download>.

/// Conservative evaporation/condensation with an evolving original gas parcel.
pub mod gas;

use crate::phase::water::{IsobaricWater, WaterState};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, EvaporationProperties};
use fs_thermochem::water::if97::{If97SaturationPair, saturation_pressure_derivative};

/// Low-pressure pure-water interface on an ORIGINAL immutable caloric chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterEvaporation {
    liquid: IsobaricWater,
    temperature_range_k: [f64; 2],
    enthalpy_range_j_kg: [f64; 2],
}
impl WaterEvaporation {
    /// Pin a caller-selected, strictly liquid part of an existing IF97 chart.
    /// Neither the chart's source references nor its original pressure is changed.
    pub fn try_new(
        liquid: IsobaricWater,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, EvaporationError> {
        let [lo, hi] = temperature_range_k;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo < 273.16
            || hi > 373.15
            || lo >= hi
            || liquid.pressure_pa() > 200_000.0
            || hi >= liquid.saturation().temperature().value()
        {
            return Err(EvaporationError::Invalid(
                "water evaporation requires a strictly liquid 273.16..=373.15 K interval below boiling and liquid pressure <=0.2 MPa",
            ));
        }
        let mut enthalpy_range_j_kg = [0.0; 2];
        for (i, t) in temperature_range_k.into_iter().enumerate() {
            let state = liquid.state_at_temperature(Temperature::new(t),|| false)
                .map_err(|_| EvaporationError::Invalid("original liquid caloric chart does not support the selected evaporation interval"))?;
            if state.liquid_mass_fraction() != 1.0 {
                return Err(EvaporationError::Invalid(
                    "water evaporation requires an all-liquid caloric branch",
                ));
            }
            enthalpy_range_j_kg[i] = state.specific_enthalpy_j_kg();
        }
        if enthalpy_range_j_kg[0] >= enthalpy_range_j_kg[1] {
            return Err(EvaporationError::Unrepresentable(
                "liquid enthalpy interval collapsed",
            ));
        }
        let law = Self {
            liquid,
            temperature_range_k,
            enthalpy_range_j_kg,
        };
        // Admit both source endpoints before publishing this immutable adapter.
        fs_thermal::evaporation::properties(&law, lo)?;
        fs_thermal::evaporation::properties(&law, hi)?;
        Ok(law)
    }

    /// Unmodified pressure-controlled liquid/heating chart.
    pub const fn liquid_chart(&self) -> &IsobaricWater {
        &self.liquid
    }

    pub const fn enthalpy_range_j_kg(&self) -> [f64; 2] {
        self.enthalpy_range_j_kg
    }

    /// Original fixed liquid pressure [Pa], separately from gas partial pressure.
    pub fn liquid_pressure_pa(&self) -> f64 {
        self.liquid.pressure_pa()
    }

    fn admit_temperature(&self, t: f64) -> Result<(), EvaporationError> {
        let [minimum_k, maximum_k] = self.temperature_range_k;
        if !t.is_finite() || !(minimum_k..=maximum_k).contains(&t) {
            return Err(EvaporationError::TemperatureOutOfRange {
                temperature_k: t,
                minimum_k,
                maximum_k,
            });
        }
        Ok(())
    }

    /// Source liquid enthalpy/density at the chart's fixed pressure.
    pub fn liquid_state(&self, temperature_k: f64) -> Result<WaterState, EvaporationError> {
        self.admit_temperature(temperature_k)?;
        self.liquid
            .state_at_temperature(Temperature::new(temperature_k), || false)
            .map_err(|_| EvaporationError::Invalid("original IF97 liquid state refused"))
    }

    /// Both unmodified IF97 saturation branches at this interface temperature.
    pub fn saturation(&self, temperature_k: f64) -> Result<If97SaturationPair, EvaporationError> {
        self.admit_temperature(temperature_k)?;
        If97SaturationPair::at_temperature(Temperature::new(temperature_k))
            .map_err(|_| EvaporationError::Invalid("original IF97 saturation source refused"))
    }

    /// Original saturated vapor-stream enthalpy in the SAME IF97 datum [J/kg].
    pub fn vapor_specific_enthalpy_j_kg(
        &self,
        temperature_k: f64,
    ) -> Result<f64, EvaporationError> {
        Ok(self
            .saturation(temperature_k)?
            .vapor()
            .specific_enthalpy()
            .value())
    }

    /// Original bounded IF97 h-to-T inverse, restricted to the admitted liquid.
    pub fn temperature_at_specific_enthalpy_k(&self, h: f64) -> Result<f64, EvaporationError> {
        let [lo, hi] = self.enthalpy_range_j_kg;
        if !h.is_finite() || !(lo..=hi).contains(&h) {
            return Err(EvaporationError::Invalid(
                "enthalpy leaves the admitted liquid evaporation branch",
            ));
        }
        // Preserve exact source endpoint temperatures despite inverse rounding.
        if h == lo {
            return Ok(self.temperature_range_k[0]);
        }
        if h == hi {
            return Ok(self.temperature_range_k[1]);
        }
        let state = self
            .liquid
            .state_at_specific_enthalpy(h, || false)
            .map_err(|_| {
                EvaporationError::Invalid("original IF97 liquid enthalpy inverse refused")
            })?
            .into_state();
        let t = state.temperature_k();
        self.admit_temperature(t)?;
        if state.liquid_mass_fraction() != 1.0 {
            return Err(EvaporationError::Invalid(
                "liquid evaporation cannot consume a boiling mixture",
            ));
        }
        Ok(t)
    }
}

impl EvaporationLaw for WaterEvaporation {
    fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }

    fn properties(&self, temperature_k: f64) -> Result<EvaporationProperties, EvaporationError> {
        let liquid = self.liquid_state(temperature_k)?;
        let pair = self.saturation(temperature_k)?;
        let pressure_derivative = saturation_pressure_derivative(Temperature::new(temperature_k))
            .map_err(|_| {
            EvaporationError::Invalid("original IF97 vapor-pressure derivative refused")
        })?;
        Ok(EvaporationProperties {
            saturation_pressure_pa: pair.pressure().value(),
            saturation_pressure_derivative_pa_k: pressure_derivative.value(),
            latent_heat_j_kg: pair.vapor().specific_enthalpy().value()
                - liquid.specific_enthalpy_j_kg(),
            liquid_density_kg_m3: liquid.bulk_density_kg_m3(),
            // Same water molecular-mass convention as the original gas species.
            molar_mass_kg_mol: 0.018_015_28,
        })
    }
}

impl ThermalLaw for WaterEvaporation {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_at_specific_enthalpy_k(h)
            .map_err(|_| DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.liquid.maximum_temperature_slope()
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.temperature_k(a)?;
        self.temperature_k(b)?;
        self.liquid
            .specific_entropy_change_j_kg_k(a, b, || false)
            .map(Some)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

#[cfg(test)]
mod tests;
