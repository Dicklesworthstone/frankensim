//! Sourced solid conductivity at the same temperature epoch as conjugate flow.
//!
//! The canonical material owner selects and refreshes the complete original
//! bundle. Only the declared absolute-temperature coordinate changes; geometry,
//! pressure, composition, source identity and all other conditions stay fixed.

use fs_matdb::MaterialCard;
use fs_material::state_point::{
    MaterialStatePointError, ResolvedMaterialStatePoint, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

/// One original solid source with an explicit grid-frame conductivity mapping.
///
/// Repeating a property for all three axes declares isotropy. Three independent
/// names declare principal axes ALREADY aligned to the voxel grid. This does
/// not rotate tensors, infer missing components or reinterpret a scalar as an
/// arbitrary anisotropic law. A new material needs only a different source
/// bundle and mapping, not a new thermal or flow solver.
#[derive(Clone, Debug)]
pub struct SolidConductivityLaw<'a> {
    card: &'a MaterialCard,
    reference: ResolvedMaterialStatePoint,
    temperature_axis: String,
    properties: [String; 3],
}

/// Original source evaluation and its actual conductivity in the geometry frame.
#[derive(Clone, Debug, PartialEq)]
pub struct SolidConductivitySample {
    /// Original voxel index, never a newly invented material index.
    pub cell: usize,
    /// Original solid assignment in the voxel domain.
    pub material_index: u16,
    /// Actual trial or accepted endpoint absolute temperature [K].
    pub temperature_k: f64,
    /// Complete original material bundle, selected claims and usage receipts.
    pub material: ResolvedMaterialStatePoint,
    /// Positive grid-aligned principal conductivities [W/(m K)].
    pub conductivity_w_m_k: [f64; 3],
}

impl<'a> SolidConductivityLaw<'a> {
    /// Bind an already resolved material without discarding its other properties.
    /// Every future query pins all original source claims. Unsupported states
    /// therefore refuse rather than switch tables or copy reference constants.
    ///
    /// `properties` is in x/y/z grid order; repeated names explicitly share the
    /// same scalar. Conductivity must have static W/(m K) quantity semantics and
    /// a positive value. Field/frequency/strain/history-dependent sources need
    /// their own state owners and cannot be frozen by this steady heat adapter.
    pub fn try_new(
        card: &'a MaterialCard,
        reference: ResolvedMaterialStatePoint,
        temperature_axis: &str,
        properties: [&str; 3],
    ) -> Result<Self, MaterialStatePointError> {
        if reference.card_identity() != card.content_hash() {
            return Err(invalid("solid conductivity requires its original material card"));
        }
        let point = reference.query();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static,
        ));
        let typed = match point.axis_quantities().get(temperature_axis) {
            Some(q) => *q == absolute || *q == QuantitySpec::dimensional(Temperature::DIMS),
            None => temperature_axis == "T",
        };
        if !typed || !point.axes().get(temperature_axis)
            .is_some_and(|t| t.is_finite() && *t > 0.0)
        {
            return Err(invalid("solid conductivity needs a positive absolute-temperature coordinate"));
        }
        let law = Self { card, reference, temperature_axis: temperature_axis.into(),
            properties: properties.map(str::to_owned) };
        law.conductivity(&law.reference)?;
        Ok(law)
    }

    /// Complete original snapshot, including any other coupled requirements.
    pub const fn reference(&self) -> &ResolvedMaterialStatePoint { &self.reference }

    /// Source coefficients at the original query, for native topology metadata.
    pub fn reference_conductivity_w_m_k(&self) -> [f64; 3] {
        self.conductivity(&self.reference).expect("immutable admitted solid source")
    }

    /// Refresh this same source at one actual solver temperature. No mutable
    /// history is committed by an equilibrium query.
    pub fn at_temperature(
        &self, cell: usize, material_index: u16, temperature_k: f64,
    ) -> Result<SolidConductivitySample, MaterialStatePointError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(invalid("solid conductivity temperature must be positive and finite"));
        }
        let point = self.reference.query().with(&self.temperature_axis, temperature_k)
            .map_err(|_| invalid("solid temperature violates its retained quantity convention"))?;
        let material = self.reference.at_query_point(self.card, &point)?;
        let conductivity_w_m_k = self.conductivity(&material)?;
        Ok(SolidConductivitySample {
            cell, material_index, temperature_k, material, conductivity_w_m_k,
        })
    }

    fn conductivity(
        &self, material: &ResolvedMaterialStatePoint,
    ) -> Result<[f64; 3], MaterialStatePointError> {
        let semantic = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity, ValueForm::Static,
        ));
        let mut values = [0.0; 3];
        for (value, name) in values.iter_mut().zip(&self.properties) {
            let property = material.property(name)
                .ok_or_else(|| invalid("solid conductivity mapping names an absent source property"))?;
            let quantity = property.requirement().quantity();
            if quantity != semantic
                && quantity != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
                || !property.is_static_at_fixed_external_conditions()
                || !property.value_si().is_finite() || property.value_si() <= 0.0
            {
                return Err(invalid("solid conductivity needs positive static thermal sources at fixed external conditions"));
            }
            *value = property.value_si();
        }
        Ok(values)
    }
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}
