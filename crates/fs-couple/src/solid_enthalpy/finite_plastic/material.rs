//! One original material card, finite constitutive family and caloric chart.

mod caloric;
pub use caloric::{ReferenceCaloricMaterial, ReferenceCaloricQuery};

use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::phase::{EquilibriumPhaseState, ResolvedHeatCapacityCurve};
use fs_material::plastic::finite::thermal::source::ResolvedAffineThermalFiniteJ2;
use fs_material::plastic::finite::thermal::{AffineThermalFiniteJ2, FiniteThermalJ2Material};

use super::FiniteThermoplasticError;

/// Thermodynamic reference of the supplied caloric enthalpy chart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CaloricReference {
    /// Stress-free caloric reference at zero pressure: h0(T) = u_caloric(T).
    /// Mechanical energy and pressure work belong to the finite solid. This is
    /// not a conversion of an arbitrary finite-pressure enthalpy chart into u.
    StressFreeZeroPressure,
}

/// Explicit same-card finite mechanical/caloric source request.
#[derive(Clone, Debug)]
pub struct FiniteThermoplasticQuery<'a> {
    pub card: &'a MaterialCard,
    pub model_index: usize,
    pub caloric: &'a ResolvedHeatCapacityCurve,
    pub context: QueryPoint,
    pub temperature_axis: String,
    pub caloric_reference: CaloricReference,
}

/// Immutable source family, preserving the exact original caloric object.
/// Density at the stress-free reference temperature supplies reference mass;
/// current physical density is rho_reference / det(F). The caloric chart's
/// density metadata at another temperature is never substituted for that mass.
#[derive(Clone, Debug)]
pub struct FiniteThermoplasticMaterial<'a> {
    query: FiniteThermoplasticQuery<'a>,
    family: ResolvedAffineThermalFiniteJ2,
    caloric: ReferenceCaloricMaterial<'a>,
    minimum_total_heat_capacity_j_m3_k: f64,
}

impl<'a> FiniteThermoplasticMaterial<'a> {
    pub fn try_new(
        query: FiniteThermoplasticQuery<'a>,
        initial_specific_enthalpy_j_kg: f64,
    ) -> Result<Self, FiniteThermoplasticError> {
        let invalid = FiniteThermoplasticError::Invalid;
        let model = query.card.models().get(query.model_index).ok_or(invalid(
            "finite thermal model index is outside the material card",
        ))?;
        let family = ResolvedAffineThermalFiniteJ2::try_from_model_card_at(
            model,
            &query.context,
            &query.temperature_axis,
        )
        .map_err(|e| FiniteThermoplasticError::Source(e.to_string()))?;
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
            return Err(invalid(
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
    pub const fn family(&self) -> &ResolvedAffineThermalFiniteJ2 {
        &self.family
    }
    pub const fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.query.caloric
    }
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.caloric.reference_density_kg_m3()
    }
    /// Conservative lower bound, including the exact chart's chord capacities.
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

impl FiniteThermalJ2Material for FiniteThermoplasticMaterial<'_> {
    fn thermal_finite_j2(&self) -> &AffineThermalFiniteJ2 {
        self.family.law()
    }
}
