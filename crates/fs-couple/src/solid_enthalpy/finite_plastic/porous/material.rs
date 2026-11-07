//! Original-card caloric binding for retained finite porous thermomechanics.
//!
//! The full non-isochoric constitutive history remains with the porous owner.
//! Caloric admission shares the finite solid's original zero-pressure chart,
//! pinned receipts, thermal domain and positive combined heat-capacity gate.
//! Reference density fixes mass; it is not replaced by temperature-dependent
//! chart metadata or inferred from the changing plastic volume.

use fs_material::phase::{EquilibriumPhaseState, ResolvedHeatCapacityCurve};
use fs_material::plastic::finite::porous::thermal::source::ResolvedAffineThermalFinitePorous;
use fs_material::plastic::finite::porous::thermal::{
    AffineThermalFinitePorous, FiniteThermalPorousMaterial,
};

use super::super::FiniteThermoplasticError;
use super::super::material::{
    FiniteThermoplasticQuery, ReferenceCaloricMaterial, ReferenceCaloricQuery,
};

/// One original material card, thermal porous family and exact caloric chart.
/// All energetic coefficients and capacities refer to the original volume;
/// current physical density is the retained reference density divided by J.
#[derive(Clone, Debug)]
pub struct FiniteThermalPorousMaterialBinding<'a> {
    query: FiniteThermoplasticQuery<'a>,
    family: ResolvedAffineThermalFinitePorous,
    caloric: ReferenceCaloricMaterial<'a>,
    minimum_total_heat_capacity_j_m3_k: f64,
}

impl<'a> FiniteThermalPorousMaterialBinding<'a> {
    /// Retain exact original source identities after complete cold-path admission.
    pub fn try_new(
        query: FiniteThermoplasticQuery<'a>,
        initial_specific_enthalpy_j_kg: f64,
    ) -> Result<Self, FiniteThermoplasticError> {
        let model =
            query
                .card
                .models()
                .get(query.model_index)
                .ok_or(FiniteThermoplasticError::Invalid(
                    "finite thermal model index is outside the material card",
                ))?;
        let family = ResolvedAffineThermalFinitePorous::try_from_model_card_at(
            model,
            &query.context,
            &query.temperature_axis,
        )
        .map_err(|error| FiniteThermoplasticError::Source(error.to_string()))?;
        let parameters = family.law().parameters();
        let caloric = ReferenceCaloricMaterial::try_new(
            ReferenceCaloricQuery {
                card: query.card,
                caloric: query.caloric,
                context: query.context.clone(),
                temperature_axis: query.temperature_axis.clone(),
                caloric_reference: query.caloric_reference,
            },
            initial_specific_enthalpy_j_kg,
            parameters.reference_temperature_k,
            parameters.temperature_range_k,
        )?;
        let minimum_total_heat_capacity_j_m3_k = caloric.minimum_heat_capacity_j_m3_k()
            + family.law().minimum_mechanical_heat_capacity_j_m3_k();
        if !minimum_total_heat_capacity_j_m3_k.is_finite()
            || minimum_total_heat_capacity_j_m3_k <= 0.0
        {
            return Err(FiniteThermoplasticError::Invalid(
                "caloric plus mechanical heat capacity must remain positive over the full admitted state domain",
            ));
        }
        Ok(Self {
            query,
            family,
            caloric,
            minimum_total_heat_capacity_j_m3_k,
        })
    }

    pub const fn query(&self) -> &FiniteThermoplasticQuery<'a> {
        &self.query
    }
    pub const fn family(&self) -> &ResolvedAffineThermalFinitePorous {
        &self.family
    }
    pub const fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.query.caloric
    }
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.caloric.reference_density_kg_m3()
    }
    /// Conservative lower bound including the exact caloric chart's chords.
    pub const fn minimum_total_heat_capacity_j_m3_k(&self) -> f64 {
        self.minimum_total_heat_capacity_j_m3_k
    }
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<EquilibriumPhaseState, fs_material::phase::PhaseStateError> {
        self.query.caloric.curve().state_at_specific_enthalpy(h)
    }
}

impl FiniteThermalPorousMaterial for FiniteThermalPorousMaterialBinding<'_> {
    fn thermal_finite_porous(&self) -> &AffineThermalFinitePorous {
        self.family.law()
    }
}
