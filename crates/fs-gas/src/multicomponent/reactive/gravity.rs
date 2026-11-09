//! Chemistry and conservative gravity share the original accepted transport.

use super::*;
use crate::multicomponent::{IsothermalMixtureGravity, MixtureGravityStep};

/// One complete reacting step in a prescribed gravitational potential.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveGravityStep {
    /// Native chemistry/transport/chemistry transaction, including gravity work
    /// exactly once in its mechanical source and gas-energy accounts.
    pub reactive: ReactiveStepReport,
    /// Native cell-center potential-energy quadrature before the whole step [J].
    pub potential_energy_before_j: f64,
    /// The same quadrature after the final chemistry half step [J].
    pub potential_energy_after_j: f64,
    /// Compensated local mass changes times the original potential [J].
    pub potential_energy_change_j: f64,
    /// Gas plus potential-energy change, including physical formation energy,
    /// less boundary transfer [J]. Numerical defects are not added as heat.
    pub total_energy_residual_j: f64,
    /// Physical gravity integrated against the actual transport-stage masses.
    pub physical_gravity_impulse_kg_m_s: f64,
    /// Represented reference-traction correction, separate from physical force.
    pub reference_traction_impulse_defect_kg_m_s: f64,
}

impl ReactingEuler1d {
    /// Advance finite-rate chemistry in the native gravity-balanced atmosphere.
    ///
    /// C(h/2), A_gravity(h), C(h/2) uses the existing conservative species,
    /// momentum and energy update. Its same Riemann mass flux supplies gravity
    /// work and potential transfer; an independent rho*u source is never added.
    /// Reaction changes the actual composition, calorics and pressure while the
    /// prepared atmosphere remains an immutable reconstruction reference.
    ///
    /// A shortened acoustic step restarts both chemistry halves from the
    /// accepted field. The complete chemical-plus-gravitational balance is
    /// checked for representability before any field, clock or report publishes.
    /// Source and geometry admission belong to `IsothermalMixtureGravity`.
    /// Zero gravity follows the original unforced numerical path.
    ///
    /// This closed, stationary, inviscid composition does not supply a molecular
    /// gravitational-diffusion law, reaction equilibrium, moving boundaries or
    /// an entropy-stability certificate. An atmosphere is mechanically at rest
    /// only while its actual chemical/thermal state also remains stationary.
    pub fn step_with_isothermal_gravity(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        reference: &IsothermalMixtureGravity,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveGravityStep, ReactiveError> {
        let mut account = GravityAccount::new(reference, &self.accepted, &mut cancelled)?;
        self.step_inner_forced::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            stabilized::MolecularMethod::Explicit(DissipationControl::default()),
            [None; 2],
            None,
            None,
            None,
            Some(&mut account),
            &mut cancelled,
        )?;
        Ok(account
            .complete
            .expect("accepted gravity composition has its balance"))
    }
}

pub(super) struct GravityAccount<'a> {
    reference: &'a IsothermalMixtureGravity,
    potential_before: f64,
    transport: Option<MixtureGravityStep>,
    complete: Option<ReactiveGravityStep>,
}

impl<'a> GravityAccount<'a> {
    fn new(
        reference: &'a IsothermalMixtureGravity,
        accepted: &MixtureEuler1d,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ReactiveError> {
        // This also checks the exact source, mesh and reflecting boundaries
        // before paying for a chemistry solve or preparing a trial state.
        let potential_before = reference.potential_energy(accepted, cancelled)?;
        Ok(Self {
            reference,
            potential_before,
            transport: None,
            complete: None,
        })
    }

    pub(super) fn begin_attempt(&mut self) {
        self.transport = None;
        self.complete = None;
    }

    pub(super) fn advance_transport(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StepControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportStepReport, ReactiveError> {
        let step =
            flow.step_with_isothermal_gravity(duration, controls, self.reference, cancelled)?;
        self.transport = Some(step);
        Ok(step.transport)
    }

    pub(super) fn finish(
        &mut self,
        before: &MixtureEuler1d,
        after: &MixtureEuler1d,
        reactive: ReactiveStepReport,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        let transport = self.transport.ok_or(ReactiveError::Invalid(
            "missing accepted reacting gravity transport",
        ))?;
        let potential_after = self.reference.potential_energy(after, &mut *cancelled)?;
        let potential_change = self
            .reference
            .potential_energy_change(before, after, cancelled)?;
        // The original reactive gas balance already subtracts gravity work.
        // Restore that one term and add potential change to audit the closed
        // gas+gravity system, including any chemical mass roundoff at its true
        // cell potential. No missing mass or energy is silently repaired.
        let total_residual = sum([
            reactive.absolute_energy_residual_j,
            reactive.transport.source_change.work,
            potential_change,
        ]);
        if !total_residual.is_finite() {
            return Err(ReactiveError::Invalid(
                "unrepresentable reacting gas plus gravitational energy balance",
            ));
        }
        self.complete = Some(ReactiveGravityStep {
            reactive,
            potential_energy_before_j: self.potential_before,
            potential_energy_after_j: potential_after,
            potential_energy_change_j: potential_change,
            total_energy_residual_j: total_residual,
            physical_gravity_impulse_kg_m_s: transport.physical_gravity_impulse_kg_m_s,
            reference_traction_impulse_defect_kg_m_s: transport
                .reference_traction_impulse_defect_kg_m_s,
        });
        Ok(())
    }
}
