//! Original caloric increments with one complete source-bound liquid endpoint.
//!
//! This prepares an immutable endpoint for a physical owner. It does not solve
//! heat transfer, advance a fluid clock, supply mechanical work or move volume.

use super::LiquidTemperatureLaw;
use crate::liquid::{
    ResolvedLiquidState,
    surface::{LiquidCapillarity, ResolvedLiquidSurfaceState},
};
use crate::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError, SolidLiquidPhase,
};
use crate::state_point::MaterialStatePointError;
use fs_qty::Temperature;

/// Source and caloric data for one proposed liquid endpoint.
///
/// Owns its source bundle and reanchored temperature law independently. A
/// consumer may borrow both into its continuing native state while this
/// endpoint stays alive; there is no self-reference or replacement caloric
/// chart. Preparing or dropping this value does not change the accepted state.
#[derive(Clone, Debug)]
pub struct LiquidCaloricEndpoint<'a> {
    initial: EquilibriumPhaseState,
    state: EquilibriumPhaseState,
    source: ResolvedLiquidSurfaceState,
    law: LiquidTemperatureLaw<'a>,
    requested_increment_j_kg: f64,
    represented_increment_j_kg: f64,
    maximum_density_defect_kg_m3: f64,
}

impl<'a> LiquidCaloricEndpoint<'a> {
    /// Original accepted state on the unchanged caloric chart.
    #[must_use]
    pub const fn initial(&self) -> EquilibriumPhaseState {
        self.initial
    }

    /// Proposed state from that same original chart and enthalpy datum.
    #[must_use]
    pub const fn state(&self) -> EquilibriumPhaseState {
        self.state
    }

    /// Noncapillary momentum source at the actual proposed temperature.
    ///
    /// The complete original bundle, including optional additional claims, is
    /// retained. Capillary stress is explicitly neglected by this endpoint
    /// API; a coupled surface owner must retain its own capillary treatment.
    #[must_use]
    pub const fn source(&self) -> &ResolvedLiquidSurfaceState {
        &self.source
    }

    /// Original fixed-context law reanchored to the same resolved endpoint.
    #[must_use]
    pub const fn law(&self) -> &LiquidTemperatureLaw<'a> {
        &self.law
    }

    /// Requested NET specific enthalpy increment; no work is added or removed.
    #[must_use]
    pub const fn net_specific_enthalpy_increment_j_kg(&self) -> f64 {
        self.requested_increment_j_kg
    }

    /// Actual endpoint-minus-initial enthalpy, before multiplying by mass.
    ///
    /// Difference from the request is floating-point storage quantization.
    /// The consuming thermal owner checks its own extensive energy budget.
    #[must_use]
    pub const fn represented_specific_enthalpy_increment_j_kg(&self) -> f64 {
        self.represented_increment_j_kg
    }

    /// Largest actual source/chart density defect across the two endpoints.
    ///
    /// This is an observed numerical/model discrepancy, not source uncertainty.
    #[must_use]
    pub const fn maximum_source_density_defect_kg_m3(&self) -> f64 {
        self.maximum_density_defect_kg_m3
    }
}

/// Original material/caloric refusal while preparing an owned liquid endpoint.
#[derive(Clone, Debug, PartialEq)]
pub enum LiquidCaloricError {
    /// The original source declined a complete liquid or surface state.
    Material(MaterialStatePointError),
    /// The original caloric chart declined the state or enthalpy increment.
    Caloric(PhaseStateError),
    /// The source and caloric objects do not describe one admitted liquid path.
    InvalidBinding {
        /// Failed original-state or source-context invariant.
        what: &'static str,
    },
    /// Original source and chart densities exceed the caller's explicit budget.
    DensityAgreement {
        /// Maximum actual initial/final source-chart density difference.
        observed_kg_m3: f64,
        /// Explicit numerical/model accuracy allowance supplied by the caller.
        maximum_allowed_kg_m3: f64,
    },
}

impl core::fmt::Display for LiquidCaloricError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for LiquidCaloricError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Caloric(error) => Some(error),
            _ => None,
        }
    }
}

impl From<MaterialStatePointError> for LiquidCaloricError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}

impl From<PhaseStateError> for LiquidCaloricError {
    fn from(error: PhaseStateError) -> Self {
        Self::Caloric(error)
    }
}

impl<'a> LiquidTemperatureLaw<'a> {
    /// Reanchor this SAME source law at a new supported absolute temperature.
    ///
    /// Refreshes the original complete bundle once. Original requirements,
    /// selected claims, viscosity convention, external conditions and any
    /// density-derived expansion path survive; reconstructing a law from only
    /// its material bundle would lose that separately retained density path.
    /// This changes no accepted state, enthalpy datum or physical inventory.
    ///
    /// # Errors
    /// Retains every original temperature, source-support and transport refusal.
    pub fn reanchored_at_temperature(
        &self,
        temperature: Temperature,
    ) -> Result<Self, MaterialStatePointError> {
        let state = self.at_temperature(temperature)?;
        Ok(self.reanchor_resolved(state))
    }

    /// Prepare a noncapillary liquid endpoint from original NET enthalpy.
    ///
    /// The original chart alone integrates/inverts the enthalpy change. Its
    /// state and material identities must match, and both endpoints must be
    /// liquid. Every retained caloric receipt must keep its selected source
    /// and complete nonthermal query context. Authored charts without source
    /// receipts retain their explicit original constitutive assumptions; this
    /// method does not promote them to source-qualified caloric data.
    ///
    /// Initial source data are reused at the reference temperature, otherwise
    /// refreshed once. The final source is resolved once after the caloric
    /// update; equal-temperature endpoints reuse the same resolved state.
    /// Both source/chart density discrepancies must fit the supplied explicit
    /// accuracy allowance. Neither density is overwritten to force agreement.
    /// The surface state and reanchored law adopt that SAME endpoint bundle.
    ///
    /// Use the original chart during any scalar heat/work solve and call this
    /// method only for its accepted candidate. No table is rebuilt or sampled
    /// into another chart, and no intermediate Newton trial needs a material
    /// query. The endpoint remains a proposal until its physical owner accepts
    /// mass, volume, transport, source work, cancellation and time together.
    ///
    /// For an explicitly isobaric Cp chart, ambient p*dV is already inside h.
    /// The supplied increment is NET thermal enthalpy after any separately
    /// owned mechanical-energy partition. This method neither subtracts that
    /// ambient work again nor turns pressure-projection/remap losses into heat.
    /// It does not infer pressure qualification or an absolute pressure datum.
    ///
    /// # Errors
    /// Refuses a mismatched chart/card/source context, nonliquid endpoint,
    /// unsupported source refresh, invalid density budget, excessive density
    /// disagreement or an unrepresentable original caloric increment.
    pub fn prepare_caloric_endpoint(
        &self,
        caloric: &EquilibriumEnthalpyPhaseCurve,
        initial: EquilibriumPhaseState,
        net_specific_enthalpy_increment_j_kg: f64,
        maximum_source_density_defect_kg_m3: f64,
    ) -> Result<LiquidCaloricEndpoint<'a>, LiquidCaloricError> {
        if !maximum_source_density_defect_kg_m3.is_finite()
            || maximum_source_density_defect_kg_m3 < 0.0
            || caloric.material_card_identity() != self.reference_state.material().card_identity()
        {
            return Err(invalid(
                "liquid caloric endpoint needs its original card and a finite density-accuracy budget",
            ));
        }
        let state =
            caloric.advance_specific_energy(initial, net_specific_enthalpy_increment_j_kg)?;
        self.finish_caloric_endpoint(
            caloric,
            initial,
            state,
            net_specific_enthalpy_increment_j_kg,
            maximum_source_density_defect_kg_m3,
        )
    }

    /// Refresh sources at an EXACT endpoint already accepted on the original chart.
    ///
    /// Native thermal owners supply their actual final state, so subtracting and
    /// re-adding h cannot move its last bit or select a different source epoch.
    /// Both states must belong to this same chart. The full original claim,
    /// fixed-context, liquid-phase and source-density checks remain unchanged.
    /// The reported requested increment equals the represented endpoint change;
    /// no heat, mechanical work or independently integrated state is inferred.
    pub fn prepare_caloric_state_endpoint(
        &self,
        caloric: &EquilibriumEnthalpyPhaseCurve,
        initial: EquilibriumPhaseState,
        state: EquilibriumPhaseState,
        maximum_source_density_defect_kg_m3: f64,
    ) -> Result<LiquidCaloricEndpoint<'a>, LiquidCaloricError> {
        self.finish_caloric_endpoint(
            caloric,
            initial,
            state,
            state.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg(),
            maximum_source_density_defect_kg_m3,
        )
    }

    fn finish_caloric_endpoint(
        &self,
        caloric: &EquilibriumEnthalpyPhaseCurve,
        initial: EquilibriumPhaseState,
        state: EquilibriumPhaseState,
        net_specific_enthalpy_increment_j_kg: f64,
        maximum_source_density_defect_kg_m3: f64,
    ) -> Result<LiquidCaloricEndpoint<'a>, LiquidCaloricError> {
        if !maximum_source_density_defect_kg_m3.is_finite()
            || maximum_source_density_defect_kg_m3 < 0.0
            || caloric.material_card_identity() != self.reference_state.material().card_identity()
        {
            return Err(invalid(
                "liquid caloric endpoint needs its original card and a finite density-accuracy budget",
            ));
        }
        if initial.phase_curve_identity() != caloric.identity()
            || state.phase_curve_identity() != caloric.identity()
        {
            return Err(PhaseStateError::StateCurveMismatch.into());
        }
        if initial.phase() != SolidLiquidPhase::Liquid || state.phase() != SolidLiquidPhase::Liquid
        {
            return Err(invalid(
                "liquid caloric endpoint cannot discard a solid or mixed phase",
            ));
        }
        for receipt in caloric.source_receipts() {
            let property = self
                .reference_state
                .material()
                .property(&receipt.property)
                .ok_or_else(|| {
                    invalid("liquid caloric source is absent from the original transport bundle")
                })?;
            if property.answer().receipt.selected != receipt.selected
                || receipt.axis_quantities != *self.reference_point.axis_quantities()
                || !receipt
                    .query_point
                    .iter()
                    .filter(|(axis, _)| axis != &self.temperature_axis)
                    .map(|(axis, value)| (axis, value))
                    .eq(self
                        .reference_point
                        .axes()
                        .iter()
                        .filter(|(axis, _)| *axis != &self.temperature_axis))
            {
                return Err(invalid(
                    "liquid caloric source changed its original claim or fixed query context",
                ));
            }
        }
        let reference_temperature = self.reference_point.axes()[&self.temperature_axis];
        let refreshed_initial = if initial.temperature_k() == reference_temperature {
            None
        } else {
            Some(self.at_temperature(Temperature::new(initial.temperature_k()))?)
        };
        let initial_source = refreshed_initial.as_ref().unwrap_or(&self.reference_state);
        let maximum_initial_defect = density_defect(initial, initial_source)?;
        let endpoint_source = if state.temperature_k() == initial.temperature_k() {
            refreshed_initial.unwrap_or_else(|| self.reference_state.clone())
        } else if state.temperature_k() == reference_temperature {
            self.reference_state.clone()
        } else {
            self.at_temperature(Temperature::new(state.temperature_k()))?
        };
        let maximum_density_defect_kg_m3 =
            maximum_initial_defect.max(density_defect(state, &endpoint_source)?);
        if maximum_density_defect_kg_m3 > maximum_source_density_defect_kg_m3 {
            return Err(LiquidCaloricError::DensityAgreement {
                observed_kg_m3: maximum_density_defect_kg_m3,
                maximum_allowed_kg_m3: maximum_source_density_defect_kg_m3,
            });
        }
        let represented_increment_j_kg =
            state.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg();
        if !represented_increment_j_kg.is_finite() {
            return Err(invalid(
                "represented liquid enthalpy increment is unrepresentable",
            ));
        }
        let source = ResolvedLiquidSurfaceState::try_from_material(
            endpoint_source.material().clone(),
            self.viscosity_basis,
            LiquidCapillarity::Neglected,
        )?;
        let law = self.reanchor_resolved(endpoint_source);
        Ok(LiquidCaloricEndpoint {
            initial,
            state,
            source,
            law,
            requested_increment_j_kg: net_specific_enthalpy_increment_j_kg,
            represented_increment_j_kg,
            maximum_density_defect_kg_m3,
        })
    }

    fn reanchor_resolved(&self, reference_state: ResolvedLiquidState) -> Self {
        Self {
            card: self.card,
            reference_point: reference_state.material().query(),
            reference_state,
            temperature_axis: self.temperature_axis.clone(),
            viscosity_basis: self.viscosity_basis,
        }
    }
}

fn density_defect(
    phase: EquilibriumPhaseState,
    source: &ResolvedLiquidState,
) -> Result<f64, LiquidCaloricError> {
    let defect = (phase.bulk_density_kg_m3() - source.density_kg_m3()).abs();
    if !defect.is_finite() {
        return Err(invalid(
            "liquid source/chart density difference is unrepresentable",
        ));
    }
    Ok(defect)
}

fn invalid(what: &'static str) -> LiquidCaloricError {
    LiquidCaloricError::InvalidBinding { what }
}

#[cfg(test)]
mod tests;
