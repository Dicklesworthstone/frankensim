//! One material/phase condition update for electrical and thermal properties.
//! The canonical QueryPoint, property resolver and source claims remain owners.

use crate::phase::EquilibriumPhaseState;
use crate::state_point::{MaterialStatePointError, ResolvedScalarProperty};
use fs_matdb::{MaterialCard, PropertyValue, QueryPoint};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// Original material and externally specified conditions of a finite thermal body.
#[derive(Clone, Copy, Debug)]
pub struct ThermalMaterialPoint<'a> {
    /// Same material card as the body's enthalpy chart.
    pub material: &'a MaterialCard,
    /// All external state coordinates; only the explicitly named axes change.
    pub context: &'a QueryPoint,
    /// Existing absolute-temperature coordinate [K].
    pub temperature_axis: &'a str,
    /// Optional existing dimensional scalar liquid MASS fraction coordinate.
    /// None means no evolving fraction coordinate, not phase independence.
    pub liquid_fraction_axis: Option<&'a str>,
}

/// Original material refusal stays visible; no representative value is inserted.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalMaterialError {
    /// Incompatible identity/quantity, absent coordinate or bad geometry.
    Invalid(&'static str),
    /// The original condition-point owner refused the update.
    Coordinate(String),
    /// The canonical property owner refused the exact source query.
    Material(MaterialStatePointError),
}
impl core::fmt::Display for ThermalMaterialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "thermal material: {s}"),
            Self::Coordinate(s) => write!(f, "thermal material coordinate: {s}"),
            Self::Material(e) => write!(f, "thermal material source: {e}"),
        }
    }
}
impl std::error::Error for ThermalMaterialError {}
fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

impl ThermalMaterialPoint<'_> {
    /// Update only T and the declared phase coordinate while retaining the
    /// original card, pressure/rate/etc. conditions and quantity descriptors.
    pub fn at(&self, body: EquilibriumPhaseState) -> Result<QueryPoint, ThermalMaterialError> {
        let invalid = ThermalMaterialError::Invalid;
        if self.material.content_hash() != body.material_card_identity() {
            return Err(invalid("property and enthalpy must use the same material card"));
        }
        let tq = self.context.axis_quantities().get(self.temperature_axis);
        if !self.context.axes().contains_key(self.temperature_axis)
            || tq.is_some_and(|q| *q != QuantitySpec::dimensional(Temperature::DIMS)
                && *q != semantic(QuantityKind::AbsoluteTemperature)) {
            return Err(invalid("an existing absolute-temperature axis is required"));
        }
        let mut point = self.context.clone().with(self.temperature_axis, body.temperature_k())
            .map_err(|e| ThermalMaterialError::Coordinate(e.to_string()))?;
        if let Some(axis) = self.liquid_fraction_axis {
            if axis == self.temperature_axis || !point.axes().contains_key(axis)
                || point.axis_quantities().get(axis).is_some_and(|q|
                    *q != QuantitySpec::dimensional(Dims([0;6]))) {
                return Err(invalid("liquid mass fraction needs a distinct existing dimensional scalar axis"));
            }
            point = point.with(axis, body.liquid_mass_fraction())
                .map_err(|e| ThermalMaterialError::Coordinate(e.to_string()))?;
        }
        Ok(point)
    }

    /// Admit a scalar static property already resolved from this binding's card.
    /// Curves may use the declared temperature or liquid-mass-fraction axis.
    /// This checks source form and frequency support, not the property's units,
    /// physical range, or continuous support between sampled states. Call `at`
    /// and the common state-point resolver before admitting the returned property.
    pub fn admit_static(&self, property: &ResolvedScalarProperty) -> Result<(), ThermalMaterialError> {
        let invalid = ThermalMaterialError::Invalid;
        let claim = self.material.claims().claim(property.answer().receipt.selected)
            .ok_or(invalid("selected thermal property claim is absent"))?;
        match &claim.value {
            PropertyValue::Scalar { .. } => {},
            PropertyValue::Curve { abscissa, abscissa_dims, .. }
                if (abscissa == self.temperature_axis && *abscissa_dims == Temperature::DIMS)
                    || (Some(abscissa.as_str()) == self.liquid_fraction_axis
                        && *abscissa_dims == Dims([0;6])) => {},
            _ => return Err(invalid("static property must be scalar or depend on a declared thermal coordinate")),
        }
        if property.angular_frequency_band_rad_s()
            .map_err(|e| ThermalMaterialError::Coordinate(format!("{e:?}")))?
            .is_some_and(|(low, _)| low > 0.0) {
            return Err(invalid("material source does not admit a static response"));
        }
        Ok(())
    }
}
