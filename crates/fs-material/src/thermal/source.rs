//! One material/phase condition update for electrical and thermal properties.
//! The canonical QueryPoint, property resolver and source claims remain owners.

use crate::phase::EquilibriumPhaseState;
use crate::state_point::{MaterialStatePointError, ResolvedScalarProperty};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

/// Total hemispherical, dimensionless surface emissivity. This is not a
/// wavelength sample, directional absorptance or a bulk optical coefficient.
pub const SURFACE_EMISSIVITY_PROPERTY: &str = "hemispherical-total-emissivity";

/// Total hemispherical emissivity is dimensionless.
pub const EMISSIVITY_DIMS: Dims = Dims::NONE;

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
    /// The original fluid transport owner refused its current physical state.
    Fluid(crate::MaterialError),
}
impl core::fmt::Display for ThermalMaterialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "thermal material: {s}"),
            Self::Coordinate(s) => write!(f, "thermal material coordinate: {s}"),
            Self::Material(e) => write!(f, "thermal material source: {e}"),
            Self::Fluid(e) => write!(f, "thermal fluid source: {e}"),
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
            return Err(invalid(
                "property and enthalpy must use the same material card",
            ));
        }
        let tq = self.context.axis_quantities().get(self.temperature_axis);
        if !self.context.axes().contains_key(self.temperature_axis)
            || tq.is_some_and(|q| {
                *q != QuantitySpec::dimensional(Temperature::DIMS)
                    && *q != semantic(QuantityKind::AbsoluteTemperature)
            })
        {
            return Err(invalid("an existing absolute-temperature axis is required"));
        }
        let mut point = self
            .context
            .clone()
            .with(self.temperature_axis, body.temperature_k())
            .map_err(|e| ThermalMaterialError::Coordinate(e.to_string()))?;
        if let Some(axis) = self.liquid_fraction_axis {
            if axis == self.temperature_axis
                || !point.axes().contains_key(axis)
                || point
                    .axis_quantities()
                    .get(axis)
                    .is_some_and(|q| *q != QuantitySpec::dimensional(Dims([0; 6])))
            {
                return Err(invalid(
                    "liquid mass fraction needs a distinct existing dimensional scalar axis",
                ));
            }
            point = point
                .with(axis, body.liquid_mass_fraction())
                .map_err(|e| ThermalMaterialError::Coordinate(e.to_string()))?;
        }
        Ok(point)
    }

    /// Admit a scalar static property already resolved from this binding's card.
    /// Source axes may be the declared temperature/liquid-mass-fraction axes or
    /// typed fixed pressure, composition, moisture and concentration conditions.
    /// This checks source form and frequency support, not the property's units,
    /// physical range, or continuous support between sampled states. Call `at`
    /// and the common state-point resolver before admitting the returned property.
    pub fn admit_static(
        &self,
        property: &ResolvedScalarProperty,
    ) -> Result<(), ThermalMaterialError> {
        self.admit_static_with_axes(property, &[])
    }

    /// Admit the same instantaneous scalar law with explicitly coupled extra
    /// coordinates. Every added axis must already have its exact quantity on
    /// the original context; thermal axes cannot be replaced or retyped.
    /// This does not solve those coordinates: the calling physical owner must
    /// close them with its own state equation. Frequency support stays static.
    pub fn admit_static_with_axes(
        &self,
        property: &ResolvedScalarProperty,
        additional_axes: &[(&str, QuantitySpec)],
    ) -> Result<(), ThermalMaterialError> {
        let invalid = ThermalMaterialError::Invalid;
        for (index, (axis, quantity)) in additional_axes.iter().enumerate() {
            if *axis == self.temperature_axis
                || Some(*axis) == self.liquid_fraction_axis
                || quantity.dims() == Temperature::DIMS
                || !self.context.axes().contains_key(*axis)
                || self.context.axis_quantities().get(*axis) != Some(quantity)
                || additional_axes[..index]
                    .iter()
                    .any(|(other, _)| other == axis)
            {
                return Err(invalid(
                    "additional static coordinates need unique existing typed nonthermal axes",
                ));
            }
        }
        let claim = self
            .material
            .claims()
            .claim(property.answer().receipt.selected)
            .ok_or(invalid("selected thermal property claim is absent"))?;
        if claim.validity.bounds().keys().any(|axis| {
            axis != self.temperature_axis
                && (self
                    .context
                    .axis_quantities()
                    .get(axis)
                    .is_some_and(|q| q.dims() == Temperature::DIMS)
                    || axis == "T" && !self.context.axis_quantities().contains_key(axis))
        }) {
            return Err(invalid(
                "a second source temperature requires an explicitly coupled thermal owner",
            ));
        }
        if !claim.value.axes().all(|(axis, dims)| {
            (axis == self.temperature_axis && dims == Temperature::DIMS)
                || (Some(axis) == self.liquid_fraction_axis && dims == Dims([0; 6]))
                || additional_axes.iter().any(|(declared_axis, quantity)| {
                    *declared_axis == axis && quantity.dims() == dims
                })
                || (axis != self.temperature_axis
                    && dims != Temperature::DIMS
                    && self
                        .context
                        .axis_quantities()
                        .get(axis)
                        .is_some_and(|q| q.dims() == dims)
                    && crate::state_point::is_fixed_external_coordinate(
                        axis,
                        self.context.axis_quantities().get(axis).copied(),
                        self.context.axes().get(axis).copied(),
                    ))
        }) {
            return Err(invalid(if additional_axes.is_empty() {
                "static property must depend only on owned thermal coordinates and typed fixed external conditions"
            } else {
                "static property must depend only on explicitly coupled coordinates"
            }));
        }
        if property
            .angular_frequency_band_rad_s()
            .map_err(|e| ThermalMaterialError::Coordinate(format!("{e:?}")))?
            .is_some_and(|(low, _)| low > 0.0)
        {
            return Err(invalid("material source does not admit a static response"));
        }
        Ok(())
    }
}
