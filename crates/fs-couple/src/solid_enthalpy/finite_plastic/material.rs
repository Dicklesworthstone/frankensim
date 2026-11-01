//! One original material card, finite constitutive family and caloric chart.

use fs_matdb::{ClaimSelection, MaterialCard, QueryPoint};
use fs_material::phase::{EquilibriumPhaseState, ResolvedHeatCapacityCurve};
use fs_material::plastic::finite::thermal::source::ResolvedAffineThermalFiniteJ2;
use fs_material::plastic::finite::thermal::{AffineThermalFiniteJ2, FiniteThermalJ2Material};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::semantic::QuantityKind;
use fs_qty::{Dims, QuantitySpec};

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
    reference_density_kg_m3: f64,
    minimum_total_heat_capacity_j_m3_k: f64,
}

impl<'a> FiniteThermoplasticMaterial<'a> {
    pub fn try_new(
        query: FiniteThermoplasticQuery<'a>,
        initial_specific_enthalpy_j_kg: f64,
    ) -> Result<Self, FiniteThermoplasticError> {
        let invalid = FiniteThermoplasticError::Invalid;
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
        let model = query.card.models().get(query.model_index).ok_or(invalid(
            "finite thermal model index is outside the material card",
        ))?;
        let family = ResolvedAffineThermalFiniteJ2::try_from_model_card_at(
            model,
            &query.context,
            &query.temperature_axis,
        )
        .map_err(|e| FiniteThermoplasticError::Source(e.to_string()))?;
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
        let parameters = family.law().parameters();
        if state.temperature_k() != parameters.reference_temperature_k
            || query.context.axes().get(&query.temperature_axis).copied()
                != Some(state.temperature_k())
            || knots[0].temperature_k < parameters.temperature_range_k[0]
            || knots[knots.len() - 1].temperature_k > parameters.temperature_range_k[1]
        {
            return Err(invalid(
                "initial caloric state must be at the model reference temperature and the complete chart inside its thermal domain",
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
        let bundle = resolve_material_state_point(
            query.card,
            &query.context,
            &requirements,
            MaterialPropertySelection::PinnedByProperty(pins),
        )
        .map_err(|e| FiniteThermoplasticError::Source(e.to_string()))?;
        if bundle
            .properties()
            .iter()
            .any(|p| !p.has_only_fixed_external_conditions(Some(&query.temperature_axis), false))
        {
            return Err(invalid(
                "caloric sources cannot freeze evolving mechanical or rate coordinates",
            ));
        }
        let reference_density_kg_m3 = bundle
            .property(&receipts[1].property)
            .ok_or(invalid(
                "original reference density role is absent from the resolved source",
            ))?
            .value_si();
        // Reference density comes from the original claim at T_ref. The
        // caloric owner's specific-volume interpolation may approximate its
        // density metadata; it neither rescales mass nor supplies current rho.
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
        for pair in knots.windows(2) {
            minimum_cp = minimum_cp.min(
                (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                    / (pair[1].temperature_k - pair[0].temperature_k),
            );
        }
        let minimum_total_heat_capacity_j_m3_k = reference_density_kg_m3 * minimum_cp
            + family.law().minimum_mechanical_heat_capacity_j_m3_k();
        if !minimum_cp.is_finite()
            || minimum_cp <= 0.0
            || !minimum_total_heat_capacity_j_m3_k.is_finite()
            || minimum_total_heat_capacity_j_m3_k <= 0.0
        {
            return Err(invalid(
                "caloric plus mechanical heat capacity must remain positive over the full admitted state domain",
            ));
        }
        Ok(Self {
            query,
            family,
            reference_density_kg_m3,
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
        self.reference_density_kg_m3
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
