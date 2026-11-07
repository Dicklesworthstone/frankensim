//! Local temperature derivatives on the original NASA-9 source regions.

use crate::{
    Nasa9RegionV1, Nasa9StandardStateModelV1, ThermochemErrorV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, source_heat_capacity,
};
use fs_qty::{Qty, Temperature};

/// Molar heat-capacity temperature derivative [J/(mol K²)].
pub type MolarHeatCapacityTemperatureDerivativeV1 = Qty<2, 1, -2, -2, 0, -1>;

/// Source support and numerical refusals specific to a local Cp tangent.
#[derive(Clone, Debug, PartialEq)]
pub enum Nasa9TemperatureDerivativeErrorV1 {
    /// The original Cp evaluator refused the temperature or source arithmetic.
    Source(ThermochemErrorV1),
    /// The original adjacent regions do not define one Cp tangent at this seam.
    NonDifferentiable {
        /// Original region immediately below the shared temperature.
        lower_region: usize,
        /// Original region selected at and above the shared temperature.
        upper_region: usize,
        /// Exact requested absolute-temperature bits.
        temperature_bits: u64,
    },
    /// Differentiating an admitted source produced a non-finite result.
    Unrepresentable {
        /// Original source region.
        region: usize,
        /// Exact rejected derivative bits.
        bits: u64,
    },
}

impl core::fmt::Display for Nasa9TemperatureDerivativeErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(error) => error.fmt(f),
            Self::NonDifferentiable {
                lower_region,
                upper_region,
                temperature_bits,
            } => write!(
                f,
                "NASA-9 Cp has no unique tangent between regions {lower_region} and {upper_region} at temperature bits {temperature_bits:#018x}",
            ),
            Self::Unrepresentable { region, bits } => write!(
                f,
                "NASA-9 Cp temperature derivative is not finite in region {region} (bits {bits:#018x})",
            ),
        }
    }
}
impl core::error::Error for Nasa9TemperatureDerivativeErrorV1 {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Source(error) => Some(error),
            _ => None,
        }
    }
}
impl From<ThermochemErrorV1> for Nasa9TemperatureDerivativeErrorV1 {
    fn from(error: ThermochemErrorV1) -> Self {
        Self::Source(error)
    }
}

impl Nasa9StandardStateModelV1 {
    /// Differentiate the original molar Cp polynomial at this temperature.
    ///
    /// Region selection and finite Cp admission are the original evaluator's.
    /// The query allocates nothing, evaluates no perturbed temperatures, and
    /// does not change caloric references. A source endpoint uses its inward
    /// derivative. At a shared seam, BOTH original Cp values and derivatives
    /// must agree to rounding tolerance; an authored jump or corner refuses.
    /// This supplies a constitutive derivative, not a derivative uncertainty
    /// bound, positive-capacity certificate or source-domain extension.
    pub fn evaluate_heat_capacity_temperature_derivative(
        &self,
        temperature: Temperature,
    ) -> Result<MolarHeatCapacityTemperatureDerivativeV1, Nasa9TemperatureDerivativeErrorV1> {
        let cp = self.evaluate_heat_capacity(temperature)?.value();
        let t = temperature.value();
        let (index, region) = self
            .region_for_temperature(t)
            .expect("original Cp admitted its region");
        let derivative = slope(region, index, t)?;
        if index > 0
            && self.regions[index - 1].temperature_max().value() == t
            && region.temperature_min().value() == t
        {
            let previous = &self.regions[index - 1];
            let previous_cp = source_heat_capacity(previous, t);
            let previous_slope = slope(previous, index - 1, t)?;
            if !rounding_equal(cp, previous_cp) || !rounding_equal(derivative, previous_slope) {
                return Err(Nasa9TemperatureDerivativeErrorV1::NonDifferentiable {
                    lower_region: index - 1,
                    upper_region: index,
                    temperature_bits: t.to_bits(),
                });
            }
        }
        Ok(MolarHeatCapacityTemperatureDerivativeV1::new(derivative))
    }
}

fn slope(
    region: &Nasa9RegionV1,
    index: usize,
    t: f64,
) -> Result<f64, Nasa9TemperatureDerivativeErrorV1> {
    let [a0, a1, _, a3, a4, a5, a6, _, _] = *region.coefficients();
    let inverse = 1.0 / t;
    let inverse2 = inverse * inverse;
    let t2 = t * t;
    let derivative = R
        * (-2.0 * ((a0 * inverse2) * inverse) - (a1 * inverse) * inverse
            + a3
            + 2.0 * (a4 * t)
            + 3.0 * (a5 * t2)
            + 4.0 * (a6 * (t2 * t)));
    if !derivative.is_finite() {
        return Err(Nasa9TemperatureDerivativeErrorV1::Unrepresentable {
            region: index,
            bits: derivative.to_bits(),
        });
    }
    Ok(derivative)
}

fn rounding_equal(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}

#[cfg(test)]
mod tests;
