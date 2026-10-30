//! Named thermal-conductivity correlations for near-ambient water and ice.
//!
//! Liquid water uses IAPWS SR6-08(2011), equation (8) and Table 6, without
//! changing its coefficients. Section 6 permits the temperature-only relation
//! between the vapor/liquid saturation pressure and 0.3 MPa. The existing IF97
//! Region 4 source supplies that admission boundary. The admitted liquid
//! temperature interval is 273.15--383.15 K; the release's optional subcooled
//! extrapolation is deliberately not exposed here. No stability against ice is
//! inferred from a liquid conductivity query.
//!
//! Ice uses the explicitly empirical scalar regression in Y.-C. Yen, *Review
//! of thermal properties of snow, ice, and sea ice*, CRREL Report 81-10 (1981),
//! equation (33), p. 15. It represents ordinary polycrystalline freshwater ice
//! at near-atmospheric conditions, not porous snow, saline ice, pressure effects
//! or a crystal-axis tensor. Its selected near-freezing operating interval is
//! 243.15--273.16 K. The report describes substantial sample/preparation scatter
//! and supplies no numerical uncertainty bound for this regression. In
//! particular, its value near melting differs from the rounded 2.2 W/(m K)
//! discussed in the same report; the original regression is not retuned.
//!
//! These are transport laws, separate from the caloric IAPWS equations. They
//! supply neither heat capacity nor latent heat, and cannot alter an enthalpy
//! reference. Their analytic temperature derivatives are derivatives of these
//! correlations, not certified material-response bounds. Evaluation is native,
//! allocation-free, and uses the repository's deterministic elementary math.
//!
//! Attribution: International Association for the Properties of Water and Steam.
//! <https://iapws.org/technical-guidance/release/LiquidWater.download>.
//! Yen's original report: <https://usace.contentdm.oclc.org/digital/collection/p266001coll1/id/7366>.

use super::if97;
use core::fmt;
use fs_math::det;
use fs_qty::{Pressure, Qty, Temperature};

/// Thermal conductivity in coherent SI [W/(m K)].
pub type ThermalConductivity = Qty<1, 1, -3, -1, 0, 0>;
/// Temperature derivative of thermal conductivity [W/(m K²)].
pub type ConductivityTemperatureDerivative = Qty<1, 1, -3, -2, 0, 0>;

pub const LIQUID_SOURCE_RELEASE: &str = "IAPWS SR6-08(2011), equation (8), Table 6";
pub const LIQUID_SOURCE_URL: &str =
    "https://iapws.org/technical-guidance/release/LiquidWater.download";
pub const ICE_SOURCE_RELEASE: &str = "Yen (1981), CRREL Report 81-10, equation (33)";
pub const ICE_SOURCE_URL: &str =
    "https://usace.contentdm.oclc.org/digital/collection/p266001coll1/id/7366";
pub const LIQUID_MINIMUM_TEMPERATURE_K: f64 = 273.15;
pub const LIQUID_MAXIMUM_TEMPERATURE_K: f64 = 383.15;
pub const LIQUID_MAXIMUM_PRESSURE_PA: f64 = 300_000.0;
/// Selected near-freezing scope of the empirical ice transport law.
pub const ICE_MINIMUM_TEMPERATURE_K: f64 = 243.15;
pub const ICE_MAXIMUM_TEMPERATURE_K: f64 = 273.16;

/// The correlation actually evaluated, including the ice model's empirical status.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaterConductivitySource {
    LiquidIapwsSr6082011,
    IceYen1981Empirical,
}

/// Source-domain or finite-representation refusal; no input is clipped.
#[derive(Clone, Debug, PartialEq)]
pub enum WaterConductivityError {
    TemperatureOutOfRange {
        source: WaterConductivitySource,
        temperature_k: f64,
        minimum_k: f64,
        maximum_k: f64,
    },
    LiquidPressureOutOfRange {
        pressure_pa: f64,
        minimum_pa: f64,
        maximum_pa: f64,
    },
    Saturation(if97::If97Error),
    Unrepresentable,
}
impl fmt::Display for WaterConductivityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "water thermal conductivity: {self:?}")
    }
}
impl std::error::Error for WaterConductivityError {}

/// Immutable value and analytic derivative from one explicitly named source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterConductivityState {
    source: WaterConductivitySource,
    temperature_k: f64,
    conductivity_w_m_k: f64,
    derivative_w_m_k2: f64,
}
impl WaterConductivityState {
    pub const fn source(self) -> WaterConductivitySource {
        self.source
    }
    pub const fn temperature(self) -> Temperature {
        Temperature::new(self.temperature_k)
    }
    pub const fn conductivity(self) -> ThermalConductivity {
        ThermalConductivity::new(self.conductivity_w_m_k)
    }
    pub const fn temperature_derivative(self) -> ConductivityTemperatureDerivative {
        ConductivityTemperatureDerivative::new(self.derivative_w_m_k2)
    }
}

/// Original SR6-08(2011) liquid relation, within its stated pressure applicability.
///
/// Pressure is used for domain admission; the admitted correlation itself is
/// temperature-only. At the lower temperature edge liquid may be metastable
/// relative to ice. Phase selection belongs to the caloric material owner.
pub fn liquid_at_pressure_temperature(
    pressure: Pressure,
    temperature: Temperature,
) -> Result<WaterConductivityState, WaterConductivityError> {
    let source = WaterConductivitySource::LiquidIapwsSr6082011;
    let t = temperature.value();
    admit_temperature(
        source,
        t,
        LIQUID_MINIMUM_TEMPERATURE_K,
        LIQUID_MAXIMUM_TEMPERATURE_K,
    )?;
    let minimum_pa = if97::saturation_pressure(temperature)
        .map_err(WaterConductivityError::Saturation)?
        .value();
    let p = pressure.value();
    if !p.is_finite() || !(minimum_pa..=LIQUID_MAXIMUM_PRESSURE_PA).contains(&p) {
        return Err(WaterConductivityError::LiquidPressureOutOfRange {
            pressure_pa: p,
            minimum_pa,
            maximum_pa: LIQUID_MAXIMUM_PRESSURE_PA,
        });
    }
    let log_reduced_temperature = det::ln(t / 300.0);
    let coefficients = [1.663_0, -1.778_1, 1.156_7, -0.432_115];
    let exponents = [-1.15, -3.4, -6.0, -7.6];
    let mut k = 0.0;
    let mut weighted = 0.0;
    for (coefficient, exponent) in coefficients.into_iter().zip(exponents) {
        let term = coefficient * det::exp(exponent * log_reduced_temperature);
        k += term;
        weighted += exponent * term;
    }
    state(source, t, k, weighted / t)
}

/// Yen's unretuned equation (33) for near-atmospheric polycrystalline fresh ice.
///
/// This temperature-only empirical source has no pressure response. Consumers
/// must explicitly retain its atmospheric-condition assumption; passing an ice
/// EOS state at arbitrary pressure does not extend this transport law's scope.
pub fn ice_at_temperature(
    temperature: Temperature,
) -> Result<WaterConductivityState, WaterConductivityError> {
    let source = WaterConductivitySource::IceYen1981Empirical;
    let t = temperature.value();
    admit_temperature(
        source,
        t,
        ICE_MINIMUM_TEMPERATURE_K,
        ICE_MAXIMUM_TEMPERATURE_K,
    )?;
    let k = 9.828 * det::exp(-0.0057 * t);
    state(source, t, k, -0.0057 * k)
}

fn admit_temperature(
    source: WaterConductivitySource,
    temperature_k: f64,
    minimum_k: f64,
    maximum_k: f64,
) -> Result<(), WaterConductivityError> {
    if temperature_k.is_finite() && (minimum_k..=maximum_k).contains(&temperature_k) {
        Ok(())
    } else {
        Err(WaterConductivityError::TemperatureOutOfRange {
            source,
            temperature_k,
            minimum_k,
            maximum_k,
        })
    }
}

fn state(
    source: WaterConductivitySource,
    temperature_k: f64,
    conductivity_w_m_k: f64,
    derivative_w_m_k2: f64,
) -> Result<WaterConductivityState, WaterConductivityError> {
    if !conductivity_w_m_k.is_finite()
        || conductivity_w_m_k <= 0.0
        || !derivative_w_m_k2.is_finite()
    {
        return Err(WaterConductivityError::Unrepresentable);
    }
    Ok(WaterConductivityState {
        source,
        temperature_k,
        conductivity_w_m_k,
        derivative_w_m_k2,
    })
}

#[cfg(test)]
mod tests;
