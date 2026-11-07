//! Original reference calorics shared by finite mechanical families.
//!
//! This admission owns no mechanical constitutive update. It retains the exact
//! original solid caloric chart, typed fixed context, pinned Cp/density receipts
//! and reference-density query. Current volume and temperature never rebind mass.

use super::CaloricReference;
use crate::solid_enthalpy::finite_plastic::FiniteThermoplasticError;
use fs_matdb::{ClaimSelection, MaterialCard, QueryPoint};
use fs_material::phase::{EquilibriumPhaseState, ResolvedHeatCapacityCurve};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

/// Explicit original caloric source and its fixed external conditions.
#[derive(Clone, Debug)]
pub struct ReferenceCaloricQuery<'a> {
    pub card: &'a MaterialCard,
    pub caloric: &'a ResolvedHeatCapacityCurve,
    pub context: QueryPoint,
    pub temperature_axis: String,
    pub caloric_reference: CaloricReference,
}

/// Source-admitted caloric storage per unit original mass. Mechanical energy
/// and any mechanical heat-capacity correction remain with the selected law.
#[derive(Clone, Debug)]
pub struct ReferenceCaloricMaterial<'a> {
    query: ReferenceCaloricQuery<'a>,
    reference_state: ResolvedMaterialStatePoint,
    reference_temperature_k: f64,
    reference_density_kg_m3: f64,
    minimum_heat_capacity_j_m3_k: f64,
}

impl<'a> ReferenceCaloricMaterial<'a> {
    /// Admit the complete solid caloric chart inside the caller's explicitly
    /// selected temperature span. This grants no mechanical or fracture-law
    /// validity outside that law's own independently checked source domain.
    pub fn try_new(
        query: ReferenceCaloricQuery<'a>,
        initial_specific_enthalpy_j_kg: f64,
        reference_temperature_k: f64,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, FiniteThermoplasticError> {
        let invalid = FiniteThermoplasticError::Invalid;
        let [minimum_k, maximum_k] = temperature_range_k;
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        if !minimum_k.is_finite()
            || !maximum_k.is_finite()
            || minimum_k <= 0.0
            || maximum_k <= minimum_k
            || !reference_temperature_k.is_finite()
            || !(minimum_k..=maximum_k).contains(&reference_temperature_k)
            || query.temperature_axis.is_empty()
            || query.context.axis_quantities().get(&query.temperature_axis) != Some(&absolute)
            || query.context.axes().get(&query.temperature_axis) != Some(&reference_temperature_k)
        {
            return Err(invalid(
                "positive caloric temperature interval and exact typed absolute reference temperature required",
            ));
        }
        match query.caloric_reference {
            CaloricReference::StressFreeZeroPressure => {}
        }
        if query
            .context
            .axis_quantities()
            .iter()
            .any(|(axis, quantity)| {
                quantity
                    .semantic_type()
                    .is_some_and(|ty| ty.kind() == QuantityKind::Pressure)
                    && query
                        .context
                        .axes()
                        .get(axis)
                        .is_some_and(|value| *value != 0.0)
            })
        {
            return Err(invalid(
                "the zero-pressure caloric reference cannot use a nonzero pressure-conditioned context",
            ));
        }
        let curve = query.caloric.curve();
        let knots = curve.knots();
        if curve.material_card_identity() != query.card.content_hash()
            || knots.len() < 2
            || knots.len() > fs_material::phase::MAX_HEAT_CAPACITY_ENTHALPY_KNOTS
            || knots.iter().any(|k| k.liquid_mass_fraction != 0.0)
            || knots
                .windows(2)
                .any(|p| p[0].temperature_k >= p[1].temperature_k)
        {
            return Err(invalid(
                "one same-card strictly increasing fully solid caloric chart required",
            ));
        }
        let state = curve
            .state_at_specific_enthalpy(initial_specific_enthalpy_j_kg)
            .map_err(FiniteThermoplasticError::Caloric)?;
        if state.temperature_k() != reference_temperature_k
            || knots[0].temperature_k < minimum_k
            || knots[knots.len() - 1].temperature_k > maximum_k
        {
            return Err(invalid(
                "initial caloric state must equal the original reference temperature and the complete chart stay inside the declared thermal domain",
            ));
        }
        let receipts = query.caloric.receipts();
        let context: Vec<_> = query
            .context
            .axes()
            .iter()
            .filter(|(axis, _)| *axis != &query.temperature_axis)
            .map(|(axis, value)| (axis.clone(), *value))
            .collect();
        if receipts.len() < 4
            || receipts.len() % 2 != 0
            || receipts.iter().any(|r| {
                r.axis_quantities != *query.context.axis_quantities()
                    || r.query_point
                        .iter()
                        .filter(|(axis, _)| axis != &query.temperature_axis)
                        .cloned()
                        .collect::<Vec<_>>()
                        != context
            })
            || receipts.chunks_exact(2).any(|pair| {
                pair[0].selected != receipts[0].selected || pair[1].selected != receipts[1].selected
            })
        {
            return Err(invalid(
                "original pinned Cp/density receipts and one fixed external context required",
            ));
        }
        let mut requirements = Vec::with_capacity(2);
        let mut pins = Vec::with_capacity(2);
        for (i, dims) in [Dims([2, 0, -2, -1, 0, 0]), Dims([-3, 1, 0, 0, 0, 0])]
            .into_iter()
            .enumerate()
        {
            let claim = query
                .card
                .claims()
                .claim(receipts[i].selected)
                .ok_or(invalid(
                    "caloric source claim is absent from the shared material card",
                ))?;
            if claim.key.quantity() != QuantitySpec::dimensional(dims) {
                return Err(invalid(
                    "original caloric Cp/density claims require their exact SI quantities",
                ));
            }
            requirements.push(
                ScalarPropertyRequirement::try_with_key(
                    &claim.key,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .map_err(|e| FiniteThermoplasticError::Source(e.to_string()))?,
            );
            pins.push((claim.key.name().to_owned(), receipts[i].selected));
        }
        let reference_state = resolve_material_state_point(
            query.card,
            &query.context,
            &requirements,
            MaterialPropertySelection::PinnedByProperty(pins),
        )
        .map_err(|e| FiniteThermoplasticError::Source(e.to_string()))?;
        if reference_state
            .properties()
            .iter()
            .any(|p| !p.has_only_fixed_external_conditions(Some(&query.temperature_axis), false))
        {
            return Err(invalid(
                "caloric sources cannot freeze evolving mechanical or rate coordinates",
            ));
        }
        let reference_density_kg_m3 = reference_state
            .property(&receipts[1].property)
            .ok_or(invalid(
                "original reference density role is absent from the resolved source",
            ))?
            .value_si();
        // Reference density comes from the original claim at T_ref. Chart
        // density metadata does not replace the reference mass inventory.
        let lower = query
            .context
            .clone()
            .with(&query.temperature_axis, knots[0].temperature_k)
            .map_err(|e| FiniteThermoplasticError::Source(e.to_string()))?;
        let upper = query
            .context
            .clone()
            .with(
                &query.temperature_axis,
                knots[knots.len() - 1].temperature_k,
            )
            .map_err(|e| FiniteThermoplasticError::Source(e.to_string()))?;
        let cp_claim = query
            .card
            .claims()
            .claim(receipts[0].selected)
            .ok_or(invalid("original heat-capacity claim is absent"))?;
        let cp = query
            .card
            .claims()
            .query_line_typed(
                &cp_claim.key,
                &query.temperature_axis,
                &lower,
                &upper,
                ClaimSelection::Pinned(receipts[0].selected),
            )
            .map_err(|e| FiniteThermoplasticError::Source(e.to_string()))?;
        let mut minimum_cp = cp
            .knots
            .iter()
            .map(|(_, sample)| sample.evidence.value.value)
            .fold(f64::INFINITY, f64::min);
        for index in 0..knots.len() - 1 {
            minimum_cp = minimum_cp.min(
                query
                    .caloric
                    .curve()
                    .sensible_heat_capacity_bounds_j_kg_k(index)
                    .map_err(|_| invalid("solid caloric heat capacity is unrepresentable"))?[0],
            );
        }
        let minimum_heat_capacity_j_m3_k = reference_density_kg_m3 * minimum_cp;
        if !minimum_cp.is_finite()
            || minimum_cp <= 0.0
            || !minimum_heat_capacity_j_m3_k.is_finite()
            || minimum_heat_capacity_j_m3_k <= 0.0
        {
            return Err(invalid(
                "reference caloric heat capacity must remain positive over the full admitted temperature domain",
            ));
        }
        Ok(Self {
            query,
            reference_state,
            reference_temperature_k,
            reference_density_kg_m3,
            minimum_heat_capacity_j_m3_k,
        })
    }

    pub const fn query(&self) -> &ReferenceCaloricQuery<'a> {
        &self.query
    }
    pub const fn reference_state(&self) -> &ResolvedMaterialStatePoint {
        &self.reference_state
    }
    pub const fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.query.caloric
    }
    pub const fn reference_temperature_k(&self) -> f64 {
        self.reference_temperature_k
    }
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.reference_density_kg_m3
    }
    pub const fn minimum_heat_capacity_j_m3_k(&self) -> f64 {
        self.minimum_heat_capacity_j_m3_k
    }
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<EquilibriumPhaseState, fs_material::phase::PhaseStateError> {
        self.query.caloric.curve().state_at_specific_enthalpy(h)
    }
}
