//! Material-bound series conduction between the two energy midpoint states.
//!
//! G = A / (d_gas/k_gas + d_body/k_body). Each layer is represented by its
//! uniform midpoint coefficient; this is not a resolved temperature profile,
//! convection, contact resistance or a phase-front transport model.
//! Reference: MIT Unified Engineering, Thermal Resistance Circuits:
//! https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node118.html

use super::ConductanceSample;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_matdb::{MaterialCard, PropertyKey, PropertyValue, QueryPoint};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// Original body source, explicit condition axes, and caller-owned layer geometry.
/// No independent gas conductivity is accepted: it comes from the same resolved
/// midpoint that drove gas energy and piston momentum.
pub struct MaterialConductionLayers<'a> {
    /// Must be the SAME card that owns the body's enthalpy/phase chart.
    pub material: &'a MaterialCard,
    /// All external material conditions. Only the declared thermal axes change.
    pub context: &'a QueryPoint,
    /// Existing absolute-temperature coordinate, in kelvin.
    pub temperature_axis: &'a str,
    /// Optional existing dimensional scalar liquid MASS fraction coordinate.
    /// Declare it when source conductivity depends on phase fraction. Absent
    /// means a source model without an evolving fraction coordinate, not a
    /// mixture rule or evidence of phase-independent conductivity.
    pub liquid_fraction_axis: Option<&'a str>,
    /// Exact scalar isotropic conductivity key, dimensional or static semantic.
    pub conductivity: &'a PropertyKey,
    /// Existing explicit claim-selection/pinning policy.
    pub selection: MaterialPropertySelection,
    /// Common transfer area [m2]; zero declares isolation, never missing data.
    pub area_m2: f64,
    /// Positive effective gas conduction distance [m].
    pub gas_distance_m: f64,
    /// Positive effective body conduction distance [m].
    pub body_distance_m: f64,
}

/// Original material refusal stays visible; no representative value is inserted.
#[derive(Clone, Debug, PartialEq)]
pub enum LayerConductionError {
    /// Incompatible identity/quantity, absent coordinate or bad geometry.
    Invalid(&'static str),
    /// The original condition-point owner refused the update.
    Coordinate(String),
    /// The canonical property owner refused the exact source query.
    Material(MaterialStatePointError),
}
impl core::fmt::Display for LayerConductionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "material conduction layers: {s}"),
            Self::Coordinate(s) => write!(f, "material conduction coordinate: {s}"),
            Self::Material(e) => write!(f, "material conduction source: {e}"),
        }
    }
}
impl std::error::Error for LayerConductionError {}
fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

impl MaterialConductionLayers<'_> {
    /// Resolve the exact source at the actual midpoint T and, when declared,
    /// liquid fraction. Preserve all other axes, quantity descriptors and source
    /// receipts. In particular a latent plateau can change conductivity via its
    /// fraction coordinate without inventing a temperature change.
    ///
    /// Both coefficients are resolved even for zero area, so missing material
    /// data cannot silently become insulation. Explicit zero body conductivity
    /// means insulation. Spatial homogenization, phase/pressure validity and the
    /// effective layer lengths remain source/caller obligations.
    pub fn resolve(
        &self, gas: &ResolvedThermochemicalGasState, body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<ResolvedMaterialStatePoint>, LayerConductionError> {
        let invalid = LayerConductionError::Invalid;
        if self.material.content_hash() != body.material_card_identity() {
            return Err(invalid("conductivity and enthalpy must use the same material card"));
        }
        if !self.area_m2.is_finite() || self.area_m2 < 0.0
            || !positive(self.gas_distance_m) || !positive(self.body_distance_m) {
            return Err(invalid("finite nonnegative area and positive layer distances required"));
        }
        let q = self.conductivity.quantity();
        if q != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
            && q != semantic(QuantityKind::ThermalConductivity) {
            return Err(invalid("scalar static thermal conductivity quantity required"));
        }
        let tq = self.context.axis_quantities().get(self.temperature_axis);
        if !self.context.axes().contains_key(self.temperature_axis)
            || tq.is_some_and(|q| *q != QuantitySpec::dimensional(Temperature::DIMS)
                && *q != semantic(QuantityKind::AbsoluteTemperature)) {
            return Err(invalid("an existing absolute-temperature axis is required"));
        }
        let mut point = self.context.clone().with(self.temperature_axis, body.temperature_k())
            .map_err(|e| LayerConductionError::Coordinate(e.to_string()))?;
        if let Some(axis) = self.liquid_fraction_axis {
            if axis == self.temperature_axis || !point.axes().contains_key(axis)
                || point.axis_quantities().get(axis).is_some_and(|q|
                    *q != QuantitySpec::dimensional(Dims([0;6]))) {
                return Err(invalid("liquid mass fraction needs a distinct existing dimensional scalar axis"));
            }
            point = point.with(axis, body.liquid_mass_fraction())
                .map_err(|e| LayerConductionError::Coordinate(e.to_string()))?;
        }
        let requirement = ScalarPropertyRequirement::try_with_key(self.conductivity, ScalarAdmissibility::NonNegative)
            .map_err(LayerConductionError::Material)?;
        let source = resolve_material_state_point(self.material, &point, &[requirement], self.selection.clone())
            .map_err(LayerConductionError::Material)?;
        let property = source.property(self.conductivity.name())
            .ok_or(invalid("resolved conductivity is absent"))?;
        let claim = self.material.claims().claim(property.answer().receipt.selected)
            .ok_or(invalid("selected conductivity claim is absent"))?;
        match &claim.value {
            PropertyValue::Scalar { .. } => {},
            PropertyValue::Curve { abscissa, abscissa_dims, .. }
                if (abscissa == self.temperature_axis && *abscissa_dims == Temperature::DIMS)
                    || (Some(abscissa.as_str()) == self.liquid_fraction_axis
                        && *abscissa_dims == Dims([0;6])) => {},
            _ => return Err(invalid("conductivity must be scalar or depend on a declared thermal coordinate")),
        }
        if property.angular_frequency_band_rad_s()
            .map_err(|e| LayerConductionError::Coordinate(format!("{e:?}")))?
            .is_some_and(|(low, _)| low > 0.0) {
            return Err(invalid("conductivity source does not admit static heat transport"));
        }
        let kb = property.value_si();
        let kg = gas.state().thermal_conductivity;
        if !positive(kg) { return Err(invalid("gas midpoint conductivity is not positive and finite")); }
        let conductance = if self.area_m2 == 0.0 || kb == 0.0 { 0.0 } else {
            let gas_r = self.gas_distance_m/kg;
            let body_r = self.body_distance_m/kb;
            let total_r = gas_r+body_r;
            let value = self.area_m2/total_r;
            if !positive(gas_r) || !positive(body_r) || !positive(total_r) || !positive(value) {
                return Err(invalid("layer resistance or conductance is unrepresentable"));
            }
            value
        };
        Ok(ConductanceSample { conductance_w_k: conductance, source })
    }
}

#[cfg(test)]
mod tests;
