//! Contact-only geometry changes at the original implicit thermal endpoint.

use super::*;
use std::cell::RefCell;

/// A private thermal trial whose only editable state is its contact patches.
///
/// Cell identities, material objects, masses, enthalpies, geometry, boundary
/// data, clock and event cursor remain owned by the original thermal step.
/// Read the actual changed coverage before resolving endpoint transport.
/// Each batch uses the existing contact-update budgets and admission. Any
/// refused batch poisons the trial, even if the source callback catches it.
pub struct ContactPatchEndpoint<'e, 'm> {
    thermal: &'e mut DepositionSimulation<'m>,
    checkpoint: &'e mut dyn FnMut() -> Result<(), DepositionError>,
    failure: Option<DepositionError>,
}

impl<'m> ContactPatchEndpoint<'_, 'm> {
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        self.thermal
    }

    /// Replace only complete, expected original patches on this private trial.
    /// Removal restores actual ambient exposure. Replacement coverage is checked
    /// across the complete batch; transient overlapping partitions are allowed.
    pub fn replace_contact_patches(
        &mut self,
        updates: &[ContactPatchUpdate],
    ) -> Result<(), DepositionError> {
        if let Some(error) = &self.failure {
            return Err(error.clone());
        }
        let result = self
            .thermal
            .apply_contact_patch_updates(updates, &mut self.checkpoint);
        if let Err(error) = &result {
            self.failure = Some(error.clone());
        }
        result
    }
}

impl<'s, 'm> BackwardEulerStep<'s, 'm> {
    /// Resolve contact geometry and transport at the same thermal endpoint.
    ///
    /// Updated patches persist only if the ordinary local and total energy
    /// gates accept this trial. Their actual coverage enters the ORIGINAL heat
    /// rate routine, including ambient cooling and prescribed reservoirs. The
    /// callback cannot change storage or advance any scheduled event. Events
    /// exactly at the target remain pending in the accepted left-limit state.
    /// Endpoint contact geometry is a backward-Euler approximation; no swept
    /// contact residence, topology-event location or convergence is inferred.
    pub fn evaluate_with_contact_patches<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(&mut ContactPatchEndpoint<'_, 'm>) -> Result<StepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_contacts_option(energy_increments_j, resolve, false, checkpoint)
    }

    /// The same endpoint contact update with the native exact caloric entropy.
    /// Contact coverage and fluxes used by both energy and entropy are identical.
    pub fn evaluate_with_contact_patches_and_entropy<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(&mut ContactPatchEndpoint<'_, 'm>) -> Result<StepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_contacts_option(energy_increments_j, resolve, true, checkpoint)
    }

    fn evaluate_contacts_option<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(&mut ContactPatchEndpoint<'_, 'm>) -> Result<StepCoefficients<T>, E>,
        with_entropy: bool,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        let checkpoint = RefCell::new(checkpoint);
        let check = || (*checkpoint.borrow_mut())();
        self.evaluate_endpoint_option(
            energy_increments_j,
            |state| {
                let mut poll = || check();
                let mut endpoint = ContactPatchEndpoint {
                    thermal: state,
                    checkpoint: &mut poll,
                    failure: None,
                };
                let result = resolve(&mut endpoint);
                if let Some(error) = endpoint.failure {
                    return Err(error.into());
                }
                result.map(|heat| {
                    (
                        TransportStepCoefficients::from(BoundaryStepCoefficients::from(heat)),
                        None,
                        BoundaryExchanges::default(),
                        SourcePowers::default(),
                    )
                })
            },
            false,
            with_entropy,
            None,
            check,
        )
    }
}
