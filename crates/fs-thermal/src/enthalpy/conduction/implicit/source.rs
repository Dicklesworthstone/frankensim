//! Local heat-generation feedback on the original implicit conduction graph.

use super::{ImplicitConductionChain, ImplicitConductionError};
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, ThermalLaw};
use std::collections::BTreeMap;

/// Signed power INTO the original cell and its derivative with respect to that
/// cell's specific stored energy. Positive power heats, negative power cools.
/// This is total cell power [W], not a volumetric density [W/m³].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitHeatSourceState {
    pub power_w: f64,
    pub power_enthalpy_derivative_w_kg_j: f64,
}

impl ImplicitHeatSourceState {
    pub(super) fn admit(self) -> Result<(), DepositionError> {
        if !self.power_w.is_finite() || !self.power_enthalpy_derivative_w_kg_j.is_finite() {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
}

/// Immutable local source law evaluated at the SAME state as heat conduction.
///
/// A physical source owner supplies its complete signed power and analytic h
/// derivative. Temperature feedback includes its own dT/dh chain rule; a
/// phase-dependent source can change across a latent plateau even at constant
/// temperature. Every evaluation must enforce the source's validity domain.
/// Speculative calls may neither advance kinetics nor mutate physical history.
///
/// This autonomous local closure does not infer an electrical circuit, reaction
/// inventory, source entropy, or cross-cell feedback. A coupled source needs its
/// corresponding physical owner. Independent authored and observed source watts
/// remain additive; the caller must not supply the same power through both.
pub trait ImplicitHeatSourceLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw;
    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError>;
}

impl<'a> ImplicitConductionChain<'a> {
    /// Attach a bounded subset of local endpoint heat sources.
    ///
    /// Original caloric objects, cell identities and source domains are checked
    /// before returning. Empty, duplicate, unknown or repeated attachments
    /// refuse. Sources join the existing joule residual and diagonal Newton
    /// tangent in backward Euler, both SDIRK2 stages and adaptive advancement.
    /// They do not introduce another timestepper or an extra energy credit.
    /// Positive feedback can destroy equilibrium or make a chosen interval
    /// inadmissible; the existing nonlinear retry and work limits still apply.
    pub fn with_heat_sources(
        mut self,
        thermal: &DepositionSimulation<'_>,
        sources: &[(CellId, &'a dyn ImplicitHeatSourceLaw)],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if !self.heat_sources.is_empty() || sources.is_empty() || sources.len() > self.ids.len() {
            return Err(ImplicitConductionError::Invalid(
                "one bounded nonempty implicit heat-source attachment is required",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        let mut bindings = BTreeMap::new();
        for &(id, law) in sources {
            checkpoint()?;
            let cell = thermal
                .cells
                .get(&id)
                .ok_or(ImplicitConductionError::Invalid(
                    "implicit heat source requires an active original cell",
                ))?;
            if !thermal.uses_material_law(cell.material_id, law.caloric_law())
                || bindings.insert(id, law).is_some()
            {
                return Err(ImplicitConductionError::Invalid(
                    "implicit heat source requires a unique cell and its original caloric object",
                ));
            }
        }
        self.heat_sources = bindings;
        self.endpoint(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }
}
