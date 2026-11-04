//! Complete two-sided conductive fluxes on the original implicit heat graph.

use super::{ImplicitConductionChain, ImplicitConductionError};
use crate::enthalpy::conduction::ConductiveLink;
use crate::enthalpy::{ContactPatch, DepositionError, DepositionSimulation, ThermalLaw};
use std::collections::BTreeMap;

/// One passive link flux and its full derivatives at the same trial state.
/// Positive power is from A to B: `Q = G * (T_A - T_B)`. The derivatives
/// differentiate this complete Q, not G, with the other enthalpy held fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitLinkState {
    pub conductance_w_k: f64,
    pub power_enthalpy_derivatives_w_kg_j: [f64; 2],
}

impl ImplicitLinkState {
    pub(super) fn admit(self) -> Result<(), DepositionError> {
        if !self.conductance_w_k.is_finite()
            || self.conductance_w_k < 0.0
            || self
                .power_enthalpy_derivatives_w_kg_j
                .iter()
                .any(|value| !value.is_finite())
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
}

/// Immutable constitutive face solve using the original ordered caloric laws.
///
/// This replaces the bulk-frozen series approximation on selected links. A
/// material owner may integrate k(T) through both half cells and their contact,
/// then supply the secant and analytic endpoint tangents. It must admit every
/// trial against its source domains, and must not advance or mutate history.
/// The supplied patch retains the original ordered cells, area and coverage.
/// The supplied geometry is the original link's complete resistance model;
/// independent contact laws cannot also own this link. No internal half-cell
/// storage, additional heat source or different time integration is introduced.
pub trait ImplicitLinkLaw {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2];
    fn state_at_specific_enthalpies(
        &self,
        enthalpies_j_kg: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError>;
}

impl<'a> ImplicitConductionChain<'a> {
    /// Attach one bounded subset of complete conductive link laws.
    ///
    /// Unselected links preserve their original coefficients and arithmetic.
    /// Caloric objects must match the actual ordered A/B cells. Duplicate IDs,
    /// overlap with source contacts and repeated attachment refuse. Source and
    /// geometric admission is probed before returning, without changing any
    /// thermal state; the same callback is used at every Newton trial.
    pub fn with_link_laws(
        mut self,
        thermal: &DepositionSimulation<'_>,
        links: &[(u64, &'a dyn ImplicitLinkLaw)],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if !self.link_laws.is_empty()
            || links.is_empty()
            || links.len() > self.conduction.links.len()
        {
            return Err(ImplicitConductionError::Invalid(
                "one bounded nonempty implicit link attachment is required",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        let mut bindings = BTreeMap::new();
        for &(id, law) in links {
            checkpoint()?;
            let link = self
                .conduction
                .links
                .get(&id)
                .ok_or(ImplicitConductionError::Invalid("inactive implicit flux link"))?;
            if self.contacts.contains_key(&id) || bindings.insert(id, law).is_some() {
                return Err(ImplicitConductionError::Invalid(
                    "implicit link flux needs a unique owner without a separate contact law",
                ));
            }
            let caloric = law.caloric_laws();
            for (side, cell) in [link.a, link.b].into_iter().enumerate() {
                if !thermal.uses_material_law(thermal.cells[&cell].material_id, caloric[side]) {
                    return Err(ImplicitConductionError::Invalid(
                        "implicit link flux requires original ordered caloric objects",
                    ));
                }
            }
        }
        self.link_laws = bindings;
        self.endpoint(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }
}
