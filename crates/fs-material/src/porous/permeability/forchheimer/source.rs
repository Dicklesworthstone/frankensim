//! Original material-card beta and applicability in the shared porous owner.
//!
//! Only a validity-wide constant beta [1/m] is admitted. A measured pressure
//! loss, dimensionless drag coefficient, apparent permeability, or a sample
//! from an evolving beta curve is not silently converted into that quantity.

use super::super::{SourceAxis, source_axis};
use super::{
    DarcyFlowConstraint, ForchheimerCalibration, MaterialError, PorousGasLayer, PorousGasSupport,
    invalid, nonnegative, positive,
};
use crate::state_point::{ResolvedMaterialStatePoint, is_fixed_external_coordinate};
use fs_qty::{Dims, QuantitySpec};

const SOURCE_SPEED_AXIS: &str = "superficial_speed";

#[derive(Clone, Debug, PartialEq)]
pub(super) struct ForchheimerSource {
    material: ResolvedMaterialStatePoint,
    property: String,
}

impl ForchheimerCalibration {
    /// Exact beta basis in rho*beta*|q|q: inverse metres, with no factor of 1/2.
    pub const BETA_DIMS: Dims = Dims([-1, 0, 0, 0, 0, 0]);
    /// Optional source applicability coordinate: the magnitude of superficial
    /// velocity, explicitly described in dimensional m/s. Pore velocity and
    /// Reynolds number are different coordinates and are not inferred.
    pub const SOURCE_SPEED_AXIS: &'static str = SOURCE_SPEED_AXIS;

    /// Bind beta from one exact selected original material property.
    ///
    /// Retains the complete source bundle, query, uncertainty and selected claim.
    /// T/p limits and any explicit superficial-speed limit are intersected with
    /// the caller's calibration support. All other admitted source conditions
    /// remain fixed. Normal and transverse constraint are explicit specimen
    /// applicability, not inferred from a property name or a scalar value.
    ///
    /// This source must be attached to a permeability layer with the same
    /// immutable material card and complete original query. Selecting k and
    /// beta together, then using `PorousGasLayer::with_material_forchheimer`,
    /// supplies that binding directly. No fluid identity, Ergun relationship,
    /// porosity, pore speed or measurement-to-model qualification is invented.
    ///
    /// # Errors
    /// Refuses wrong quantity, a nonconstant beta, unsupported source axes,
    /// disjoint support, and a source speed range that excludes zero: current
    /// shared consumers enforce a maximum speed, not a nonzero lower bound.
    pub fn from_material(
        source: ResolvedMaterialStatePoint,
        beta_property: &str,
        normal: [f64; 3],
        constraint: DarcyFlowConstraint,
        mut support: PorousGasSupport,
    ) -> Result<Self, MaterialError> {
        support.validate()?;
        let property = source
            .property(beta_property)
            .ok_or_else(|| invalid("missing selected Forchheimer beta source property"))?;
        if property.requirement().quantity() != QuantitySpec::dimensional(Self::BETA_DIMS)
            || !property.source_axes().is_empty()
            || property.answer().receipt.decision
                != fs_matdb::EvaluationDecision::ConstantWithinValidity
            || !nonnegative(property.value_si())
        {
            return Err(invalid(
                "Forchheimer beta source must be a nonnegative validity-wide constant in exact inverse metres",
            ));
        }
        let beta = property.value_si();
        // Retain the whole resolved bundle's original applicability, just as
        // the intrinsic-permeability owner does. A companion claim cannot be
        // silently detached when both were selected for this specimen state.
        for property in source.properties() {
            let answer = property.answer();
            for (axis, &bounds) in answer.evidence.model.validity.bounds() {
                let quantity = answer.receipt.axis_quantities.get(axis).copied();
                let value = source
                    .query_point()
                    .iter()
                    .find(|(name, _)| name == axis)
                    .map(|(_, value)| *value);
                admit_source_coordinate(&mut support, axis, quantity, value, bounds)?;
            }
        }
        let provenance = format!(
            "source-bound Forchheimer beta {beta_property} on material {:?}",
            source.card_identity()
        );
        let mut calibration = Self::try_new(beta, normal, constraint, support, provenance)?;
        calibration.source = Some(std::sync::Arc::new(ForchheimerSource {
            material: source,
            property: beta_property.into(),
        }));
        Ok(calibration)
    }

    /// Unchanged original source claims, units, query coordinates and uncertainty.
    pub fn source(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.source.as_ref().map(|source| &source.material)
    }

    /// Exact selected beta property name, alongside its complete original bundle.
    pub fn source_property(&self) -> Option<&str> {
        self.source.as_ref().map(|source| source.property.as_str())
    }

    pub(super) fn admit_layer_source(&self, layer: &PorousGasLayer) -> Result<(), MaterialError> {
        let Some(source) = &self.source else {
            return Ok(());
        };
        let original = layer.source().ok_or_else(|| invalid(
            "source-backed Forchheimer beta requires the original source-backed permeability layer"))?;
        if original.card_identity() != source.material.card_identity()
            || original.query_point() != source.material.query_point()
            || original.axis_quantities() != source.material.axis_quantities()
        {
            return Err(invalid(
                "Forchheimer beta and intrinsic permeability must retain the same material card and source context",
            ));
        }
        Ok(())
    }
}

impl PorousGasLayer {
    /// Resolve beta from this layer's unchanged original material bundle and
    /// attach its directional calibration. Source-selected k and beta therefore
    /// reach every existing gas and continuum consumer through the same layer.
    pub fn with_material_forchheimer(
        &self,
        beta_property: &str,
        normal: [f64; 3],
        constraint: DarcyFlowConstraint,
        support: PorousGasSupport,
    ) -> Result<Self, MaterialError> {
        let source = self.source().ok_or_else(|| {
            invalid(
                "material Forchheimer construction requires this layer's original source bundle",
            )
        })?;
        let calibration = ForchheimerCalibration::from_material(
            source.clone(),
            beta_property,
            normal,
            constraint,
            support,
        )?;
        self.with_forchheimer_calibration(calibration)
    }
}

// One admitted applicability coordinate. The same maximum-speed field is
// already enforced by compressible gas, frozen graph and native liquid owners.
fn admit_source_coordinate(
    support: &mut PorousGasSupport,
    axis: &str,
    quantity: Option<QuantitySpec>,
    value: Option<f64>,
    bounds: (f64, f64),
) -> Result<(), MaterialError> {
    if axis == SOURCE_SPEED_AXIS {
        if quantity != Some(QuantitySpec::dimensional(fs_qty::Velocity::DIMS))
            || !value.is_some_and(nonnegative)
            || bounds.0 != 0.0
            || !positive(bounds.1)
        {
            return Err(invalid(
                "Forchheimer source superficial_speed needs explicit m/s and support from zero to a positive maximum",
            ));
        }
        support.maximum_superficial_speed_m_s = support.maximum_superficial_speed_m_s.min(bounds.1);
        return support.validate();
    }
    match source_axis(axis, quantity)? {
        Some(SourceAxis::Temperature) => {
            if !value.is_some_and(positive) {
                return Err(invalid(
                    "Forchheimer source needs positive absolute temperature",
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
                    "Forchheimer source needs positive absolute pressure",
                ));
            }
            support.pressure_pa = (
                support.pressure_pa.0.max(bounds.0),
                support.pressure_pa.1.min(bounds.1),
            );
        }
        None if is_fixed_external_coordinate(axis, quantity, value) => {}
        None => {
            return Err(invalid(
                "Forchheimer source has an unsupported evolving coordinate",
            ));
        }
    }
    support.validate()
}

#[cfg(test)]
mod domain_tests;
#[cfg(test)]
mod tests;
