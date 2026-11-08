//! One original material card, finite constitutive family and caloric chart.

mod caloric;
pub use caloric::{ReferenceCaloricMaterial, ReferenceCaloricQuery};

use fs_matdb::{ClaimSelection, MaterialCard, QueryPoint};
use fs_material::phase::{EquilibriumPhaseState, ResolvedHeatCapacityCurve};
use fs_material::plastic::finite::thermal::source::ResolvedAffineThermalFiniteJ2;
use fs_material::plastic::finite::thermal::source::recovery::ResolvedAffineThermalRecoveringCreep;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2ThermalError, FiniteJ2ThermalIncrement,
    FiniteThermalInelasticMaterial, FiniteThermalJ2Material,
};

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
pub struct FiniteThermoplasticMaterial<'a, F = ResolvedAffineThermalFiniteJ2> {
    query: FiniteThermoplasticQuery<'a>,
    family: F,
    caloric: ReferenceCaloricMaterial<'a>,
    minimum_total_heat_capacity_j_m3_k: f64,
}

/// Source-bound recovering finite creep on the original caloric chart and mass.
pub type FiniteRecoveringThermoplasticMaterial<'a> =
    FiniteThermoplasticMaterial<'a, ResolvedAffineThermalRecoveringCreep>;

/// Complete typed finite material plus its exact reference caloric inventory.
/// The coupled solver uses this seam without copying a constitutive law or
/// replacing a material's internal history with a plasticity-only state.
pub trait FiniteThermoplasticMaterialFamily: FiniteThermalInelasticMaterial {
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve;
    fn reference_density_kg_m3(&self) -> f64;
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
        Self::from_family(query, initial_specific_enthalpy_j_kg, family)
    }

    /// Admit original same-card temperature/plastic-strain observations and the
    /// original caloric chart on the existing reciprocal energy owner. The
    /// tabulated thermal law must be selected explicitly by its model card.
    pub fn try_new_with_tabulated_hardening(
        query: FiniteThermoplasticQuery<'a>,
        initial_specific_enthalpy_j_kg: f64,
        curve_selection: ClaimSelection,
        maximum_source_samples: usize,
    ) -> Result<Self, FiniteThermoplasticError> {
        let family = ResolvedAffineThermalFiniteJ2::try_from_material_card_at(
            query.card,
            query.model_index,
            &query.context,
            &query.temperature_axis,
            curve_selection,
            maximum_source_samples,
        )
        .map_err(|e| FiniteThermoplasticError::Source(e.to_string()))?;
        Self::from_family(query, initial_specific_enthalpy_j_kg, family)
    }
}

impl<'a> FiniteRecoveringThermoplasticMaterial<'a> {
    /// Select an explicit recovering finite model and retain its original card.
    pub fn try_new_recovering(
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
        let family = ResolvedAffineThermalRecoveringCreep::try_from_model_card_at(
            model,
            &query.context,
            &query.temperature_axis,
        )
        .map_err(|error| FiniteThermoplasticError::Source(error.to_string()))?;
        Self::from_family(query, initial_specific_enthalpy_j_kg, family)
    }
}

impl<'a, F> FiniteThermoplasticMaterial<'a, F> {
    fn from_family(
        query: FiniteThermoplasticQuery<'a>,
        initial_specific_enthalpy_j_kg: f64,
        family: F,
    ) -> Result<Self, FiniteThermoplasticError>
    where
        F: FiniteThermalInelasticMaterial,
    {
        let invalid = FiniteThermoplasticError::Invalid;
        let caloric = ReferenceCaloricMaterial::try_new(
            ReferenceCaloricQuery {
                card: query.card,
                caloric: query.caloric,
                context: query.context.clone(),
                temperature_axis: query.temperature_axis.clone(),
                caloric_reference: query.caloric_reference,
            },
            initial_specific_enthalpy_j_kg,
            family.reference_temperature_k(),
            family.temperature_range_k(),
        )?;
        let minimum_total_heat_capacity_j_m3_k = caloric.minimum_heat_capacity_j_m3_k()
            + family.minimum_mechanical_heat_capacity_j_m3_k();
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
    pub const fn family(&self) -> &F {
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

impl<F> FiniteThermoplasticMaterialFamily for FiniteThermoplasticMaterial<'_, F>
where
    Self: FiniteThermalInelasticMaterial,
{
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve {
        FiniteThermoplasticMaterial::caloric_source(self)
    }
    fn reference_density_kg_m3(&self) -> f64 {
        FiniteThermoplasticMaterial::reference_density_kg_m3(self)
    }
}

impl FiniteThermalInelasticMaterial for FiniteRecoveringThermoplasticMaterial<'_> {
    type State = <ResolvedAffineThermalRecoveringCreep as FiniteThermalInelasticMaterial>::State;
    type Frozen = <ResolvedAffineThermalRecoveringCreep as FiniteThermalInelasticMaterial>::Frozen;

    fn initial_state(&self) -> Self::State {
        self.family.initial_state()
    }
    fn reference_temperature_k(&self) -> f64 {
        self.family.reference_temperature_k()
    }
    fn temperature_range_k(&self) -> [f64; 2] {
        self.family.temperature_range_k()
    }
    fn minimum_mechanical_heat_capacity_j_m3_k(&self) -> f64 {
        self.family.minimum_mechanical_heat_capacity_j_m3_k()
    }
    fn at_temperature_for_inspection(
        &self,
        temperature_k: f64,
    ) -> Result<Self::Frozen, FiniteJ2ThermalError> {
        self.family.at_temperature_for_inspection(temperature_k)
    }
    fn at_temperature_for_step(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self::Frozen, FiniteJ2ThermalError> {
        self.family
            .at_temperature_for_step(temperature_k, duration_s)
    }
    fn checked_increment_for_step(
        &self,
        old_f: &[f64; 9],
        old: &Self::State,
        old_temperature_k: f64,
        new_f: &[f64; 9],
        new_temperature_k: f64,
        yield_tolerance_pa: f64,
        duration_s: f64,
    ) -> Result<FiniteJ2ThermalIncrement<Self::State>, FiniteJ2ThermalError> {
        self.family.checked_increment_for_step(
            old_f,
            old,
            old_temperature_k,
            new_f,
            new_temperature_k,
            yield_tolerance_pa,
            duration_s,
        )
    }
}
