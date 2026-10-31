//! Material-bound series conduction between the two energy midpoint states.
//!
//! G = A / (d_gas/k_gas + d_body/k_body). Each layer is represented by its
//! uniform midpoint coefficient; this is not a resolved temperature profile,
//! convection, contact resistance or a phase-front transport model.
//! Reference: MIT Unified Engineering, Thermal Resistance Circuits:
//! https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node118.html

use super::ConductanceSample;
/// Compatibility name for the shared thermal/electrical material-source refusal.
pub use super::source::ThermalMaterialError as LayerConductionError;
use super::source::ThermalMaterialPoint;
#[cfg(test)]
use fs_matdb::PropertyValue;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
#[cfg(test)]
use fs_qty::{Dims, Temperature};

/// Original body source, explicit condition axes, and caller-owned layer geometry.
/// No independent gas conductivity is accepted: it comes from the same resolved
/// midpoint that drove gas energy and piston momentum. The existing fluid
/// transport interface retains the actual ideal or real-fluid state.
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

fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

impl MaterialConductionLayers<'_> {
    fn binding(&self) -> ThermalMaterialPoint<'_> {
        ThermalMaterialPoint { material: self.material, context: self.context,
            temperature_axis: self.temperature_axis, liquid_fraction_axis: self.liquid_fraction_axis }
    }

    fn validate_geometry_and_key(&self) -> Result<(), LayerConductionError> {
        if !self.area_m2.is_finite() || self.area_m2 < 0.0
            || !positive(self.gas_distance_m) || !positive(self.body_distance_m) {
            return Err(LayerConductionError::Invalid("finite nonnegative area and positive layer distances required"));
        }
        let q = self.conductivity.quantity();
        if q != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
            && q != semantic(QuantityKind::ThermalConductivity) {
            return Err(LayerConductionError::Invalid("scalar static thermal conductivity quantity required"));
        }
        Ok(())
    }

    /// Resolve the exact source at the actual midpoint T and, when declared,
    /// liquid fraction. Preserve all other axes, quantity descriptors and source
    /// receipts. In particular a latent plateau can change conductivity via its
    /// fraction coordinate without inventing a temperature change.
    ///
    /// Both coefficients are resolved even for zero area, so missing material
    /// data cannot silently become insulation. Explicit zero body conductivity
    /// means insulation. Spatial homogenization, phase/pressure validity and the
    /// effective layer lengths remain source/caller obligations.
    pub fn resolve<G: FluidTransport + ?Sized>(
        &self,
        gas: &G,
        body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<ResolvedMaterialStatePoint>, LayerConductionError> {
        let point = self.binding().at(body)?;
        self.validate_geometry_and_key()?;
        let requirement = ScalarPropertyRequirement::try_with_key(self.conductivity, ScalarAdmissibility::NonNegative)
            .map_err(LayerConductionError::Material)?;
        let source = resolve_material_state_point(self.material, &point, &[requirement], self.selection.clone())
            .map_err(LayerConductionError::Material)?;
        self.resolve_from_material(gas, body, source)
    }

    /// Consume an existing complete material snapshot without querying a second
    /// property bundle. Electrical/thermal network consumers retain the SAME
    /// source identity and receipts, including auxiliary electrical properties.
    ///
    /// The card, every coordinate and its quantity convention must match the
    /// actual body phase state and this binding's external conditions. The exact
    /// conductivity key, static-law admission, nonnegative value and any explicit
    /// conductivity pin are checked even at zero area. A pre-resolved snapshot
    /// owns its policy selection; policies are not run again against another
    /// set of requirements. This does not refresh a stale material state.
    pub fn resolve_from_material<G: FluidTransport + ?Sized>(
        &self,
        gas: &G,
        body: EquilibriumPhaseState,
        source: ResolvedMaterialStatePoint,
    ) -> Result<ConductanceSample<ResolvedMaterialStatePoint>, LayerConductionError> {
        let invalid = LayerConductionError::Invalid;
        let binding = self.binding();
        let point = binding.at(body)?;
        self.validate_geometry_and_key()?;
        if source.card_identity() != self.material.content_hash()
            || source.axis_quantities() != point.axis_quantities()
            || source.query_point().len() != point.axes().len()
            || source.query_point().iter().zip(point.axes()).any(|((a, x), (b, y))|
                a != b || x.to_bits() != y.to_bits()) {
            return Err(invalid("layer snapshot must retain the exact material and physical condition point"));
        }
        let property = source.property_by_key(self.conductivity)
            .ok_or(invalid("resolved conductivity key or quantity is absent"))?;
        if let MaterialPropertySelection::PinnedByProperty(pins) = &self.selection {
            if pins.len() != 1 || pins[0].0 != self.conductivity.name()
                || pins[0].1 != property.answer().receipt.selected {
                return Err(invalid("layer snapshot does not retain the explicit conductivity pin"));
            }
        }
        binding.admit_static(property)?;
        let kb = property.value_si();
        let kg = gas
            .transport_state()
            .map_err(LayerConductionError::Fluid)?
            .thermal_conductivity_w_m_k();
        if !kb.is_finite() || kb < 0.0 {
            return Err(invalid("body conductivity is not nonnegative and finite"));
        }
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
pub(in crate::cavity_phs::piston) mod tests;
