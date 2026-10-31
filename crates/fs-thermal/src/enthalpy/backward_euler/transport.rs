//! Endpoint-evaluated conservative transfers on the original implicit residual.
//!
//! A transfer debits and credits the same power, independently of conductive
//! surface coverage. The supplying physical model, not this carrier, owns its
//! closure and accuracy. In particular directed mass-carried enthalpy must not
//! be converted into a symmetric Fourier conductance or external heater pair.
use super::*;
use super::super::advection::Advected;

/// Complete endpoint transport. Existing source/boundary coefficients retain
/// their meaning. `transfers` use the SAME native pairwise power owner as the
/// explicit observer path and never contribute to external energy supply.
#[derive(Debug)]
pub struct TransportStepCoefficients<T> {
    pub heat: BoundaryStepCoefficients<T>,
    /// Positive from a to b. Bounded by the original thermal max_links budget.
    pub transfers: Vec<HeatTransfer>,
}
impl<T> From<BoundaryStepCoefficients<T>> for TransportStepCoefficients<T> {
    fn from(heat: BoundaryStepCoefficients<T>) -> Self {
        Self { heat, transfers: Vec::new() }
    }
}

impl Advected<()> {
    /// Query the ORIGINAL upwind mass-flow graph at a proposed endpoint.
    ///
    /// This reuses its source and conservative-transfer evaluation, including
    /// per-cell fixed-mass admission, source caloric/datum checks, component
    /// reference subtraction and prescribed stationary-species admission. The
    /// original observer, thermal state, inventory and audit never change.
    /// No accepted step, explicit CFL restriction or entropy-flow integration
    /// is implied by this instantaneous endpoint query.
    ///
    /// Feed this result into the native implicit residual (or combine its
    /// external sources and transfers with another physical source). Recirculating
    /// energy stays in `transfers`; inlet/outlet enthalpy stays in `sources_w`.
    /// The fixed graph is cloned as bounded scratch; no nonlinear source law or
    /// transport equation is duplicated. Refining time and the upwind control
    /// volumes remains necessary even when the algebraic residual is small.
    pub fn endpoint_transport(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<TransportStepCoefficients<()>, DepositionError> {
        checkpoint()?;
        let mut proposal = self.clone();
        proposal.begin()?;
        let sources_w = proposal.sources(thermal, &mut checkpoint)?;
        let transfers = proposal.heat_transfers(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(TransportStepCoefficients {
            heat: StepCoefficients { sources_w, conductances_w_k: BTreeMap::new(),
                payload: () }.into(),
            transfers,
        })
    }
}

impl<'m> BackwardEulerStep<'_, 'm> {
    /// Evaluate external sources, surface/link conductances and directed internal
    /// energy transfers at ONE proposed endpoint. The existing cell-energy
    /// residual and local/global acceptance gates remain authoritative.
    /// No integrator, caloric law, accepted history or numerical heat is added.
    /// Fixed topology/inventory, source repeatability and time refinement remain
    /// caller obligations. A transfer's entropy is not assumed to be conduction.
    pub fn evaluate_with_transport<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(&DepositionSimulation<'m>) -> Result<TransportStepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_transport_option(energy_increments_j, resolve, false, checkpoint)
    }

    /// Same conservative residual with the original caloric entropy account.
    /// Signed transfer Q/T contributions remain unclassified until their owner
    /// supplies entropy-flow/mixing semantics; they are not passive-conduction
    /// generation merely because they exchange energy internally.
    pub fn evaluate_with_transport_and_entropy<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(&DepositionSimulation<'m>) -> Result<TransportStepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_transport_option(energy_increments_j, resolve, true, checkpoint)
    }
}

#[cfg(test)]
#[path = "transport/tests.rs"]
mod tests;
