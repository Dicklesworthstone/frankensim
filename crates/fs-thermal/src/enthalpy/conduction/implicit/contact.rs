//! Original contact constitutive tangents in the existing enthalpy residual.
//!
//! With `D=d_a/k_a+R''+d_b/k_b`, the native flux remains `G*(T_b-T_a)`
//! into A, where `G=A/D`. Contact adds `-G²/A * dR''/dh_i` to each bulk
//! conductance derivative. The original Jacobian converts h derivatives to
//! energy-increment derivatives using the corresponding invariant cell mass.
//! Opposite heat rows receive the same transfer with opposite signs.

use super::{ImplicitConductionChain, ImplicitConductionError, finite};
use crate::enthalpy::{DepositionError, DepositionSimulation, ThermalLaw};
use std::collections::BTreeMap;

/// One source-owned area-specific resistance and its two enthalpy derivatives.
/// Derivatives hold the other side's enthalpy and all external conditions fixed.
/// A latent phase derivative can be nonzero even when both dT/dh values vanish.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitContactState {
    pub resistance_m2_k_w: f64,
    pub resistance_enthalpy_derivatives_m2_kg_k_w_j: [f64; 2],
}
impl ImplicitContactState {
    pub(super) fn admit(self) -> Result<(), DepositionError> {
        if !self.resistance_m2_k_w.is_finite()
            || self.resistance_m2_k_w < 0.0
            || self
                .resistance_enthalpy_derivatives_m2_kg_k_w_j
                .iter()
                .any(|v| !v.is_finite())
            || (self.resistance_m2_k_w == 0.0
                && self.resistance_enthalpy_derivatives_m2_kg_k_w_j != [0.0; 2])
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }

    pub(super) fn conductance_derivative(
        self,
        g: f64,
        area: f64,
        side: usize,
    ) -> Result<f64, DepositionError> {
        let derivative = self.resistance_enthalpy_derivatives_m2_kg_k_w_j[side];
        if g == 0.0 || derivative == 0.0 {
            return Ok(0.0);
        }
        // The resistance fraction R''/D is bounded. Do not form G² directly.
        finite(
            -g * finite((g / area) * self.resistance_m2_k_w)?
                * (derivative / self.resistance_m2_k_w),
        )
    }
}

/// Immutable two-sided law on the SAME caloric objects as the original cells.
/// Material owners select sources and close temperature/phase chain rules.
/// Geometry, area, pressure-work and interfacial heat storage are not supplied
/// by this interface. Preparation and trial evaluations must not mutate history.
pub trait ImplicitContactLaw {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2];
    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2];
    fn state_at_specific_enthalpies(
        &self,
        enthalpies_j_kg: [f64; 2],
    ) -> Result<ImplicitContactState, DepositionError>;
}

impl<'a> ImplicitConductionChain<'a> {
    /// Attach a bounded subset of original link contact laws before advancing.
    ///
    /// Ordered caloric objects must match the actual link's A/B cells. Every
    /// selected geometric link must declare zero fixed contact resistance:
    /// the attached source then owns that entire physical resistance once.
    /// Unselected links retain their original coefficients and arithmetic.
    /// Replacing or accumulating an already-attached source set is refused.
    /// Attachment probes actual source states without changing thermal history.
    pub fn with_contacts(
        mut self,
        thermal: &DepositionSimulation<'_>,
        contacts: &[(u64, &'a dyn ImplicitContactLaw)],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if !self.contacts.is_empty()
            || contacts.is_empty()
            || contacts.len() > self.conduction.links.len()
        {
            return Err(ImplicitConductionError::Invalid(
                "one bounded nonempty implicit contact attachment is required",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        let mut bindings = BTreeMap::new();
        for &(id, law) in contacts {
            checkpoint()?;
            let link = self
                .conduction
                .links
                .get(&id)
                .ok_or(ImplicitConductionError::Invalid(
                    "inactive implicit contact link",
                ))?;
            if link.geometry.contact_resistance_m2_k_w != 0.0 || bindings.insert(id, law).is_some()
            {
                return Err(ImplicitConductionError::Invalid(
                    "implicit source contact must replace zero authored resistance on a unique link",
                ));
            }
            let caloric = law.caloric_laws();
            let bounds = law.enthalpy_bounds_j_kg();
            for (side, cell) in [link.a, link.b].into_iter().enumerate() {
                let [lo, hi] = bounds[side];
                if !lo.is_finite()
                    || !hi.is_finite()
                    || lo >= hi
                    || !thermal.uses_material_law(thermal.cells[&cell].material_id, caloric[side])
                {
                    return Err(ImplicitConductionError::Invalid(
                        "implicit contact requires finite enthalpy bounds and original ordered caloric objects",
                    ));
                }
            }
        }
        self.contacts = bindings;
        self.endpoint(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }
}
