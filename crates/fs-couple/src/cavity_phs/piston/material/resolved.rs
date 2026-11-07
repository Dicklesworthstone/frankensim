//! Source-derived rigid inertia on the original resolved shock/piston owners.
//!
//! The retained body supplies total translational mass; the gas mesh alone owns
//! area and position. Native Euler steps own waves, moving volumes, midpoint
//! mechanics, acceptance and joint publication. No lumped chamber substitutes
//! for the resolved field. Duration arguments are maxima; each native report
//! fixes the actual accepted interval. The original density source stays fixed:
//! gas heating does not re-query solid mass or deposit heat in the specimen.
//! This guided axial wall does not model deformation, rotation, thermal
//! expansion or general multidimensional FSI.

mod shared;

use fs_gas::dynamics::{Euler1d, StepControl};
use fs_gas::dynamics::piston::{
    Piston, PistonControl, PistonError, PistonLoad, PistonStep, PistonSuspension,
    SuspendedPistonStep,
};
use fs_gas::multicomponent::{MixtureEuler1d, StepReport as MixtureStep};
use fs_gas::multicomponent::pressure_equilibrium::{
    ContactEquilibriumControl, ContactEquilibriumStepReport, ContactWorkspace,
};
use fs_mbd_material::ResolvedBody;
use fs_qty::Time;

/// A source-bound rigid specimen constrained to translate with one axial wall.
///
/// All regions translate together, so the canonical body's total mass applies;
/// rotational inertia and centre-of-mass offsets are not an effective axial
/// mass. Geometry coverage, rigidity, guide reactions and the association with
/// a continuing gas field remain the caller's physical-model obligations.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialResolvedPiston {
    body: ResolvedBody,
    piston: Piston,
}

impl MaterialResolvedPiston {
    /// Bind original geometry/density receipts to native inertia, with no
    /// independent numeric mass. This creates a new wall, not a live mass swap.
    pub fn new(body: ResolvedBody, velocity_m_s: f64) -> Result<Self, PistonError> {
        let piston = Piston::new(body.inertia().mass(), velocity_m_s)?;
        Ok(Self { body, piston })
    }

    /// Original immutable geometry, density queries and canonical inertia.
    #[must_use]
    pub const fn body(&self) -> &ResolvedBody {
        &self.body
    }

    /// Accepted native mass and velocity; no mutable mass-replacement escape.
    #[must_use]
    pub const fn piston(&self) -> &Piston {
        &self.piston
    }

    /// Advance frozen-composition pressure waves and the material wall together.
    /// Native source, geometry, CFL, balance and cancellation refusals preserve
    /// the accepted gas, mesh, clock, piston and original material source.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self, gas: &mut Euler1d, duration: Time, gas_controls: StepControl,
        load: PistonLoad, controls: PistonControl, cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep, PistonError> {
        gas.step_with_piston(
            &mut self.piston, duration, gas_controls, load, controls, cancelled,
        )
    }

    /// Add the native stationary spring/damper to the same pressure trials.
    /// The native report separates spring storage, damping loss and drive work;
    /// changing the suspension between calls remains external constitutive work.
    #[allow(clippy::too_many_arguments)]
    pub fn step_suspended(
        &mut self, gas: &mut Euler1d, duration: Time, gas_controls: StepControl,
        load: PistonLoad, suspension: PistonSuspension, controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep, PistonError> {
        gas.step_with_suspended_piston(
            &mut self.piston, duration, gas_controls, load, suspension, controls, cancelled,
        )
    }

    /// Retain every transported species and its original caloric reference.
    /// This does not separately advance chemistry, molecular transport or heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_mixture(
        &mut self, gas: &mut MixtureEuler1d, duration: Time, gas_controls: StepControl,
        load: PistonLoad, controls: PistonControl, cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<MixtureStep>, PistonError> {
        gas.step_with_piston(
            &mut self.piston, duration, gas_controls, load, controls, cancelled,
        )
    }

    /// Compose native species transport with the original spring/damper owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_mixture_suspended(
        &mut self, gas: &mut MixtureEuler1d, duration: Time, gas_controls: StepControl,
        load: PistonLoad, suspension: PistonSuspension, controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep<MixtureStep>, PistonError> {
        gas.step_with_suspended_piston(
            &mut self.piston, duration, gas_controls, load, suspension, controls, cancelled,
        )
    }

    /// Use the native opt-in caloric-contact correction and its separate energy
    /// defect account. Caller-owned scratch, pressure work and physical balance
    /// budgets retain their original owners; no defect is relabeled as heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_mixture_with_contact_equilibrium(
        &mut self, gas: &mut MixtureEuler1d, duration: Time, gas_controls: StepControl,
        load: PistonLoad, controls: PistonControl, equilibrium: ContactEquilibriumControl,
        workspace: &mut ContactWorkspace, cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<ContactEquilibriumStepReport<MixtureStep>>, PistonError> {
        gas.step_with_piston_and_contact_equilibrium(
            &mut self.piston, duration, gas_controls, load, controls,
            equilibrium, workspace, cancelled,
        )
    }

    /// Preserve the same contact admission with native spring storage and loss.
    /// Rejection or cancellation publishes none of the accepted physical owners.
    #[allow(clippy::too_many_arguments)]
    pub fn step_mixture_suspended_with_contact_equilibrium(
        &mut self, gas: &mut MixtureEuler1d, duration: Time, gas_controls: StepControl,
        load: PistonLoad, suspension: PistonSuspension, controls: PistonControl,
        equilibrium: ContactEquilibriumControl, workspace: &mut ContactWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep<ContactEquilibriumStepReport<MixtureStep>>, PistonError> {
        gas.step_with_suspended_piston_and_contact_equilibrium(
            &mut self.piston, duration, gas_controls, load, suspension, controls,
            equilibrium, workspace, cancelled,
        )
    }
}

#[cfg(test)]
mod tests;
