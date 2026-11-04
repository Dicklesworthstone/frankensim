//! Pressure-dependent thermal support exchange with retained solid inertia.
//!
//! The original contact Newton residual owns every thermal proposal. Its same
//! mechanical endpoint now includes the accepted velocity and consistent mass;
//! no trial advances the live history. Physical plastic/reversible heat enters
//! each original cell once. Kinetic change and backward-Euler loss close the
//! combined work balance, but the numerical loss never becomes material heat.
//!
//! Contacts keep their original fixed, complete, prescribed reference faces.
//! This does not infer separation, impact, friction or evolving contact area.

use super::super::contact::{
    ContactCandidate, ContactError, ContactReport, ContactSolve, ThermalSupportContact,
};
use super::{ThermoplasticDynamics, ThermoplasticDynamicsReport};
use fs_exec::Cx;
use fs_solid::plastic3::IncrementSettings;
use fs_thermal::enthalpy::DepositionSimulation;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;

/// Publish the finite thermal bodies and their inertial solid together.
pub struct DynamicContactCandidate<'m, 't, S> {
    thermal: DepositionSimulation<'t>,
    dynamics: ThermoplasticDynamics<'m>,
    report: ThermoplasticDynamicsReport<ContactReport<S>>,
}
impl<'m, 't, S> DynamicContactCandidate<'m, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn dynamics(&self) -> &ThermoplasticDynamics<'m> {
        &self.dynamics
    }
    pub const fn report(&self) -> &ThermoplasticDynamicsReport<ContactReport<S>> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermoplasticDynamics<'m>,
        ThermoplasticDynamicsReport<ContactReport<S>>,
    ) {
        (self.thermal, self.dynamics, self.report)
    }
}

fn finish<'m, 't, S, E>(
    candidate: ContactCandidate<'m, 't, S>,
    motion: Option<super::Motion>,
    cx: &Cx<'_>,
) -> Result<DynamicContactCandidate<'m, 't, S>, ContactError<E>> {
    let motion = motion.ok_or(ContactError::Invalid(
        "inertial contact candidate must retain native velocity and kinetic account",
    ))?;
    let (thermal, solid, material) = candidate.into_parts();
    cx.checkpoint().map_err(|_| ContactError::Cancelled)?;
    Ok(DynamicContactCandidate {
        thermal,
        dynamics: ThermoplasticDynamics {
            solid,
            velocity_mps: motion.velocity_mps,
        },
        report: ThermoplasticDynamicsReport {
            material,
            inertia: motion.energy,
        },
    })
}

impl<'m> ThermoplasticDynamics<'m> {
    /// Close contact pressure, both temperatures and actual material motion in
    /// one endpoint solve. Prescribed displacement and nodal force are endpoint
    /// values. Every proposal restarts from the same accepted u/v/J2/enthalpy.
    /// The callback supplies unrelated material transport and external heat;
    /// contact conductance and mechanical heat must not be supplied twice.
    /// Resume with BOTH returned owners. Failure leaves every accepted input intact.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_thermal_contacts<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        contacts: &[ThermalSupportContact<'_>],
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: ContactSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cx: &Cx<'_>,
    ) -> Result<DynamicContactCandidate<'m, 't, S>, ContactError<E>> {
        let (candidate, motion) = self.solid.advance_contacts_candidate(
            thermal,
            contacts,
            target_s,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
            guess_j,
            resolve,
            false,
            Some(&self.velocity_mps),
            cx,
        )?;
        finish(candidate, motion, cx)
    }

    /// Same inertial/contact residual with source-resolved bulk conduction in
    /// series. The callback supplies fresh BULK-ONLY conductance for each named
    /// support link. The original interface owner combines that coefficient
    /// with the pressure/temperature-resolved resistance, including insulation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_series_thermal_contacts<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        contacts: &[ThermalSupportContact<'_>],
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: ContactSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cx: &Cx<'_>,
    ) -> Result<DynamicContactCandidate<'m, 't, S>, ContactError<E>> {
        let (candidate, motion) = self.solid.advance_contacts_candidate(
            thermal,
            contacts,
            target_s,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
            guess_j,
            resolve,
            true,
            Some(&self.velocity_mps),
            cx,
        )?;
        finish(candidate, motion, cx)
    }
}
