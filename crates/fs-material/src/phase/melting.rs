//! Source-backed sharp solid/liquid melting on the original enthalpy chart.
//!
//! The caller declares an isothermal equilibrium transition and supplies its
//! source keys, both sensible branches and complete source context. This is
//! not a phase-diagram, kinetic, hysteresis or supercooling model. Independently
//! tabulated enthalpy is deliberately not added to integrated Cp and fusion heat.

use super::{
    EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery, MAX_HEAT_CAPACITY_ENTHALPY_KNOTS,
    PhaseStateError, ResolvedHeatCapacityCurve, SolidLiquidPhase, heat_capacity_source_endpoints,
};
use fs_blake3::DomainHasher;
use fs_matdb::{MaterialCard, PropertyKey, PropertyUsageReceipt, QueryPoint, SelectionPolicy};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const SPECIFIC_ENERGY: Dims = Dims([2, 0, -2, 0, 0, 0]);

/// Complete source selection for one declared sharp equilibrium transition.
#[derive(Clone, Copy, Debug)]
pub struct MeltingSourceQuery<'a> {
    /// Solid Cp/density support ending exactly at the source melting point.
    pub solid: HeatCapacitySourceQuery<'a>,
    /// Liquid Cp/density support beginning exactly at that same melting point.
    pub liquid: HeatCapacitySourceQuery<'a>,
    /// Source equilibrium melting temperature, typed as absolute temperature.
    pub melting_point: &'a PropertyKey,
    /// Source mass-specific fusion enthalpy, J/kg, at the transition context.
    pub latent_heat_fusion: &'a PropertyKey,
    /// Complete transition point, including its absolute-temperature coordinate.
    /// Every non-temperature coordinate must match both sensible branches,
    /// except for the explicitly declared phase selector.
    pub transition: &'a QueryPoint,
    pub transition_selection: SelectionPolicy,
    /// Optional source table selector: zero for solid, one for liquid. It is
    /// dimensionless and categorical, not an evolving liquid mass fraction.
    /// The transition query omits it. With None, branches must share all
    /// non-temperature context and may use independently named property keys.
    pub phase_selector: Option<&'a str>,
}

/// One original h/T/fraction/density chart with every source usage retained.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedMeltingCurve {
    curve: EquilibriumEnthalpyPhaseCurve,
    melting_temperature_k: f64,
    latent_heat_fusion_j_kg: f64,
    plateau_enthalpy_bounds_j_kg: [f64; 2],
}

impl ResolvedMeltingCurve {
    /// Integrate sourced solid/liquid Cp and insert fusion energy exactly once.
    ///
    /// Both branches use the existing heat-capacity owner, including its
    /// continuous-support checks, retained breakpoints and interpolation bound.
    /// Their temperatures must meet the independently queried melting point
    /// exactly; no endpoint extension, gap filling or transition smoothing is
    /// performed. Pressure, composition and other source context cannot change
    /// between the branches and transition query.
    ///
    /// The declared enthalpy datum applies at the solid branch's lower endpoint.
    /// `maximum_interpolation_error_j_kg` bounds sensible-branch interpolation
    /// relative to the admitted piecewise-linear Cp laws, not source uncertainty
    /// or floating-point roundoff. The combined chart retains the original
    /// generated-knot budget. Fusion is an exact-temperature plateau in that
    /// chart; its width must remain representable at the chosen enthalpy datum.
    ///
    /// Density is independently sourced on each branch. The original chart
    /// interpolates specific volume through coexistence. A consumer keeps
    /// specimen mass invariant and handles the corresponding volume change;
    /// this constructor neither creates mass nor asserts fixed-volume geometry.
    pub fn try_new(
        card: &MaterialCard,
        query: MeltingSourceQuery<'_>,
        reference_specific_enthalpy_j_kg: f64,
        maximum_interpolation_error_j_kg: f64,
    ) -> Result<Self, PhaseStateError> {
        let axis = query.solid.temperature_axis;
        if axis.is_empty() || query.liquid.temperature_axis != axis {
            return Err(invalid(
                "melting branches require the same temperature axis",
            ));
        }
        if query.melting_point.quantity() != absolute_temperature()
            || query.latent_heat_fusion.quantity() != QuantitySpec::dimensional(SPECIFIC_ENERGY)
        {
            return Err(invalid(
                "melting requires absolute-temperature and mass-specific fusion-energy keys",
            ));
        }
        let solid_span = heat_capacity_source_endpoints(query.solid)?;
        let liquid_span = heat_capacity_source_endpoints(query.liquid)?;
        validate_context(query)?;
        let answers = [query.melting_point, query.latent_heat_fusion].map(|key| {
            card.claims()
                .query_typed(key, query.transition, query.transition_selection)
                .map_err(|error| PhaseStateError::SourceQuery {
                    property: key.name().into(),
                    reason: error.to_string(),
                })
        });
        let [melting, fusion] = answers;
        let melting = melting?;
        let fusion = fusion?;
        let melting_temperature_k = melting.evidence.value.value;
        let latent_heat_fusion_j_kg = fusion.evidence.value.value;
        if !melting_temperature_k.is_finite()
            || melting_temperature_k <= 0.0
            || !latent_heat_fusion_j_kg.is_finite()
            || latent_heat_fusion_j_kg <= 0.0
        {
            return Err(invalid(
                "melting temperature and fusion energy must be finite and positive",
            ));
        }
        if solid_span[1] != melting_temperature_k
            || liquid_span[0] != melting_temperature_k
            || query.transition.axes().get(axis) != Some(&melting_temperature_k)
        {
            return Err(invalid(
                "solid end, liquid start and transition query must equal the sourced melting point",
            ));
        }

        let solid = ResolvedHeatCapacityCurve::try_new(
            card,
            query.solid,
            SolidLiquidPhase::Solid,
            reference_specific_enthalpy_j_kg,
            maximum_interpolation_error_j_kg,
        )?;
        let solid_end_h = solid
            .curve()
            .knots()
            .last()
            .expect("admitted solid branch")
            .specific_enthalpy_j_kg;
        let liquid_start_h = solid_end_h + latent_heat_fusion_j_kg;
        let represented_latent_heat = liquid_start_h - solid_end_h;
        if !liquid_start_h.is_finite()
            || !represented_latent_heat.is_finite()
            || represented_latent_heat <= 0.0
        {
            return Err(invalid(
                "fusion energy is not representable at the selected enthalpy datum",
            ));
        }
        let liquid = ResolvedHeatCapacityCurve::try_new(
            card,
            query.liquid,
            SolidLiquidPhase::Liquid,
            liquid_start_h,
            maximum_interpolation_error_j_kg,
        )?;
        let count = solid.curve().knots().len() + liquid.curve().knots().len();
        if count > MAX_HEAT_CAPACITY_ENTHALPY_KNOTS {
            return Err(invalid(
                "combined melting chart exceeds the generated-knot budget",
            ));
        }
        let mut knots = Vec::with_capacity(count);
        knots.extend_from_slice(solid.curve().knots());
        knots.extend_from_slice(liquid.curve().knots());
        let mut curve = EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), knots)?;
        let mut receipts = solid.receipts().to_vec();
        receipts.push(melting.receipt);
        receipts.push(fusion.receipt);
        receipts.extend_from_slice(liquid.receipts());
        let mut identity = DomainHasher::new("org.frankensim.fs-material.sourced-melting-curve.v1");
        identity.update(curve.identity().as_bytes());
        identity.update(solid.curve().identity().as_bytes());
        identity.update(liquid.curve().identity().as_bytes());
        identity.update(&(receipts.len() as u64).to_le_bytes());
        for receipt in &receipts {
            identity.update(receipt.content_hash().as_bytes());
        }
        curve.identity = identity.finalize();
        curve.source_receipts = receipts;
        Ok(Self {
            curve,
            melting_temperature_k,
            latent_heat_fusion_j_kg,
            plateau_enthalpy_bounds_j_kg: [solid_end_h, liquid_start_h],
        })
    }

    /// Original equilibrium chart consumed by existing thermal/phase owners.
    #[must_use]
    pub const fn curve(&self) -> &EquilibriumEnthalpyPhaseCurve {
        &self.curve
    }

    /// Solid Cp/density receipts, transition temperature/fusion receipts, then
    /// liquid Cp/density receipts. Conditions and selected claims remain intact.
    #[must_use]
    pub fn receipts(&self) -> &[PropertyUsageReceipt] {
        self.curve.source_receipts()
    }

    #[must_use]
    pub const fn melting_temperature_k(&self) -> f64 {
        self.melting_temperature_k
    }

    /// Admitted source value, before finite-precision addition to the datum.
    #[must_use]
    pub const fn latent_heat_fusion_j_kg(&self) -> f64 {
        self.latent_heat_fusion_j_kg
    }

    /// Actual representable enthalpies at fully solid and fully liquid melting.
    #[must_use]
    pub const fn plateau_enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.plateau_enthalpy_bounds_j_kg
    }
}

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn invalid(what: &'static str) -> PhaseStateError {
    PhaseStateError::InvalidCurve { what }
}

fn validate_context(query: MeltingSourceQuery<'_>) -> Result<(), PhaseStateError> {
    let axis = query.solid.temperature_axis;
    if let Some(selector) = query.phase_selector {
        if selector.is_empty() || selector == axis {
            return Err(invalid(
                "phase selector must be a distinct nonempty source axis",
            ));
        }
        for (branch, value) in [(query.solid, 0.0), (query.liquid, 1.0)] {
            for point in [branch.lower, branch.upper] {
                if point.axes().get(selector) != Some(&value)
                    || point.axis_quantities().get(selector)
                        != Some(&QuantitySpec::dimensional(Dims::NONE))
                {
                    return Err(invalid(
                        "phase selector must be typed solid zero or liquid one",
                    ));
                }
            }
        }
        if query.transition.axes().contains_key(selector)
            || query.transition.axis_quantities().contains_key(selector)
        {
            return Err(invalid(
                "transition query must not select a single-phase table",
            ));
        }
    }
    let shared = |name: &str| name != axis && Some(name) != query.phase_selector;
    for point in [
        query.solid.lower,
        query.solid.upper,
        query.liquid.lower,
        query.liquid.upper,
        query.transition,
    ] {
        if point.axis_quantities().get(axis) != Some(&absolute_temperature())
            || !point
                .axes()
                .get(axis)
                .is_some_and(|t| t.is_finite() && *t > 0.0)
        {
            return Err(invalid(
                "melting queries require explicitly typed positive absolute temperatures",
            ));
        }
        if !point
            .axes()
            .iter()
            .filter(|(name, _)| shared(name))
            .eq(query
                .transition
                .axes()
                .iter()
                .filter(|(name, _)| shared(name)))
            || !point
                .axis_quantities()
                .iter()
                .filter(|(name, _)| shared(name))
                .eq(query
                    .transition
                    .axis_quantities()
                    .iter()
                    .filter(|(name, _)| shared(name)))
        {
            return Err(invalid(
                "melting branches and transition require identical external source context",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
