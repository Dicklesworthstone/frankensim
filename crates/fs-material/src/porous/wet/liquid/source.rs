//! Selected material-card liquid retention parameters and applicability.

use super::{
    Arc, MaterialError, WetPoreLiquidLaw, WetPoreLiquidSupport, fraction, invalid, positive,
};
use crate::porous::permeability::{PorousGasLayer, SourceAxis, source_axis};
use crate::state_point::{ResolvedMaterialStatePoint, is_fixed_external_coordinate};
use fs_qty::{Dims, QuantitySpec};

/// Exact selected property names. Saturations are liquid/gas pore-volume
/// fractions, not relative humidity or wet/dry mass moisture content. P0 is Pa,
/// never pressure head or an inverse-length/pressure alpha coefficient.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WetPoreLiquidProperties<'a> {
    pub residual_liquid_saturation: &'a str,
    pub residual_gas_saturation: &'a str,
    pub van_genuchten_m: &'a str,
    pub connectivity_exponent: &'a str,
    pub pressure_scale_pa: &'a str,
}
impl<'a> WetPoreLiquidProperties<'a> {
    fn names(self) -> [&'a str; 5] {
        [
            self.residual_liquid_saturation,
            self.residual_gas_saturation,
            self.van_genuchten_m,
            self.connectivity_exponent,
            self.pressure_scale_pa,
        ]
    }
    fn from_names([slr, sgr, m, ell, p0]: [&'a str; 5]) -> Self {
        Self {
            residual_liquid_saturation: slr,
            residual_gas_saturation: sgr,
            van_genuchten_m: m,
            connectivity_exponent: ell,
            pressure_scale_pa: p0,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct LiquidSource {
    material: ResolvedMaterialStatePoint,
    properties: [String; 5],
}

impl WetPoreLiquidLaw {
    /// Explicit dimensional-one coordinate for liquid pore-volume saturation.
    pub const SOURCE_SATURATION_AXIS: &'static str = "liquid_saturation";
    /// Magnitude of superficial velocity through bulk area, explicitly in m/s.
    pub const SOURCE_SPEED_AXIS: &'static str = "superficial_speed";

    /// Bind five selected validity-wide constant parameters to their complete
    /// original material bundle, observations, uncertainty and query context.
    /// Actual T/absolute liquid-p/speed/saturation limits intersect all selected
    /// source domains. Unsupported evolving axes and sampled curves refuse.
    /// The caller explicitly identifies the saturation basis through these keys;
    /// equal dimensions do not convert humidity or mass moisture into saturation.
    /// On continuum attachment, this calibration must share the intrinsic
    /// permeability layer's original material card and complete query context.
    pub fn from_material(
        source: ResolvedMaterialStatePoint,
        properties: WetPoreLiquidProperties<'_>,
        mut support: WetPoreLiquidSupport,
    ) -> Result<Self, MaterialError> {
        support.validate()?;
        let names = properties.names();
        let mut values = [0.0; 5];
        for (i, name) in names.iter().enumerate() {
            let property = source
                .property(name)
                .ok_or_else(|| invalid("missing selected liquid retention parameter"))?;
            let dims = if i == 4 {
                fs_qty::Pressure::DIMS
            } else {
                Dims::NONE
            };
            if property.requirement().quantity() != QuantitySpec::dimensional(dims)
                || !property.source_axes().is_empty()
                || property.answer().receipt.decision
                    != fs_matdb::EvaluationDecision::ConstantWithinValidity
                || !property.value_si().is_finite()
                || property.value_si() < 0.0
            {
                return Err(invalid(
                    "liquid retention parameters need nonnegative validity-wide constants in exact dimensionless pore-fraction/shape or Pa bases",
                ));
            }
            values[i] = property.value_si();
        }
        for property in source.properties() {
            let answer = property.answer();
            for (axis, &bounds) in answer.evidence.model.validity.bounds() {
                let quantity = answer.receipt.axis_quantities.get(axis).copied();
                let value = source.state_coordinate(axis);
                admit_coordinate(&mut support, axis, quantity, value, bounds)?;
            }
        }
        let [slr, sgr, m, ell, p0] = values;
        let provenance = format!(
            "source-bound liquid retention on material {:?}",
            source.card_identity()
        );
        let mut law = Self::try_new(slr, sgr, m, ell, p0, support, provenance)?;
        law.1 = Some(Arc::new(LiquidSource {
            material: source,
            properties: names.map(str::to_owned),
        }));
        Ok(law)
    }
    /// Complete original immutable claims, query, observations and uncertainty.
    pub fn source(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.1.as_ref().map(|source| &source.material)
    }
    pub fn source_properties(&self) -> Option<WetPoreLiquidProperties<'_>> {
        self.1.as_ref().map(|source| {
            WetPoreLiquidProperties::from_names(source.properties.each_ref().map(String::as_str))
        })
    }
    pub(crate) fn admit_layer_source(&self, layer: &PorousGasLayer) -> Result<(), MaterialError> {
        let Some(source) = &self.1 else { return Ok(()) };
        let original = layer.source().ok_or_else(|| {
            invalid("source-backed liquid retention needs source-backed intrinsic permeability")
        })?;
        if original.card_identity() != source.material.card_identity()
            || original.query_point() != source.material.query_point()
            || original.axis_quantities() != source.material.axis_quantities()
        {
            return Err(invalid(
                "liquid retention and intrinsic permeability must share the original material card and query context",
            ));
        }
        Ok(())
    }
}

fn admit_coordinate(
    support: &mut WetPoreLiquidSupport,
    axis: &str,
    quantity: Option<QuantitySpec>,
    value: Option<f64>,
    bounds: (f64, f64),
) -> Result<(), MaterialError> {
    if axis == WetPoreLiquidLaw::SOURCE_SATURATION_AXIS {
        if quantity != Some(QuantitySpec::dimensional(Dims::NONE))
            || !value.is_some_and(fraction)
            || !fraction(bounds.0)
            || !fraction(bounds.1)
            || bounds.0 > bounds.1
        {
            return Err(invalid(
                "liquid_saturation source axis requires explicit dimensional-one pore-volume fraction and support in [0,1]",
            ));
        }
        support.liquid_saturation = (
            support.liquid_saturation.0.max(bounds.0),
            support.liquid_saturation.1.min(bounds.1),
        );
    } else if axis == WetPoreLiquidLaw::SOURCE_SPEED_AXIS {
        if quantity != Some(QuantitySpec::dimensional(fs_qty::Velocity::DIMS))
            || !value.is_some_and(|v| v.is_finite() && v >= 0.0)
            || bounds.0 != 0.0
            || !positive(bounds.1)
        {
            return Err(invalid(
                "liquid source superficial_speed requires m/s support from zero to a finite positive maximum",
            ));
        }
        support.maximum_superficial_speed_m_s = support.maximum_superficial_speed_m_s.min(bounds.1);
    } else {
        match source_axis(axis, quantity)? {
            Some(SourceAxis::Temperature) => {
                if !value.is_some_and(positive) {
                    return Err(invalid(
                        "liquid retention source needs positive absolute temperature",
                    ));
                }
                support.temperature_k = (
                    support.temperature_k.0.max(bounds.0),
                    support.temperature_k.1.min(bounds.1),
                );
            }
            Some(SourceAxis::Pressure) => {
                if !value.is_some_and(positive) {
                    return Err(invalid(
                        "liquid retention source needs positive absolute liquid pressure",
                    ));
                }
                support.liquid_pressure_pa = (
                    support.liquid_pressure_pa.0.max(bounds.0),
                    support.liquid_pressure_pa.1.min(bounds.1),
                );
            }
            None if is_fixed_external_coordinate(axis, quantity, value) => {}
            None => {
                return Err(invalid(
                    "liquid retention source has an unsupported evolving coordinate",
                ));
            }
        }
    }
    support.validate()
}

#[cfg(test)]
mod tests;
