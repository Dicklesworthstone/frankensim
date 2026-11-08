//! Second-order original-enthalpy transport on an admitted frozen momentum field.
//!
//! The native thermal owner supplies SSPRK2 stages, local error control and
//! accepted heat/mass ledgers. This adapter retains the solved face mass flux,
//! selected spatial reconstruction and original material/geometry admission.
//! Only `MaterialFlowAdvection<_, ()>` is supported: a live momentum solve or
//! arbitrary inner observer cannot be advanced as a repeatable stage query.

use super::{FluidTransport, MaterialFlowAdvection};
pub use fs_conduction::transient::deposition::advection::explicit::{
    AdaptiveSsprk2Controls, AdaptiveSsprk2Error, AdaptiveSsprk2Report, Ssprk2Controls,
    Ssprk2Report, Ssprk2TimeError,
};
use fs_conduction::transient::deposition::backward_euler::{
    BoundaryStepCoefficients, StepCoefficients,
};
use fs_conduction::transient::deposition::{DepositionError, DepositionSimulation};
use fs_exec::CancelGate;
use std::collections::BTreeMap;

fn checkpoint(gate: &CancelGate) -> Result<(), DepositionError> {
    if gate.is_requested() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}

fn unchanged(
    _: &DepositionSimulation<'_>,
) -> Result<BoundaryStepCoefficients<()>, DepositionError> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    }
    .into())
}

impl<F: FluidTransport + ?Sized> MaterialFlowAdvection<'_, F, ()> {
    /// Advance the original thermal enthalpy with SSPRK2 and the already
    /// solved mass flux. Donor or selected limited traces, including latent
    /// storage, remain with the native thermal owner. Each stage observes the
    /// original combined conduction/advection CFL and caloric domain.
    ///
    /// Returns the accepted material-flow wrapper with its actual transport
    /// ledger and reconstruction retained. Frozen momentum is unchanged; this
    /// does not add phase mobility or infer momentum/pressure-work heating.
    pub fn advance_ssprk2_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: Ssprk2Controls,
        gate: &CancelGate,
    ) -> Result<(DepositionSimulation<'m>, Self, Ssprk2Report), DepositionError> {
        self.advance_ssprk2_with_coefficients(thermal, target_s, controls, unchanged, gate)
    }

    /// Stage-local source watts and link/ambient/reservoir conductances through
    /// the original coefficient carrier. Queries see actual trial enthalpy and
    /// time, may be repeated/discarded, and must not advance physical history.
    /// The callback must honor authored boundary state at native events; an
    /// endpoint query uses the pre-event state even at the event's exact time.
    ///
    /// The original fluid cell map, mass, geometry and caloric-law identity are
    /// admitted before queries and before final publication. Native transfers,
    /// source heat and accepted flow throughput are each credited once.
    pub fn advance_ssprk2_with_coefficients<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: Ssprk2Controls,
        mut coefficients: impl FnMut(
            &DepositionSimulation<'m>,
        ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
        gate: &CancelGate,
    ) -> Result<(DepositionSimulation<'m>, Self, Ssprk2Report), DepositionError> {
        self.admit(thermal, &mut || checkpoint(gate))?;
        let (next, advected, report) = self.advected.advance_ssprk2_with_coefficients(
            thermal,
            target_s,
            controls,
            |stage| {
                self.admit(stage, &mut || checkpoint(gate))?;
                coefficients(stage)
            },
            || checkpoint(gate),
        )?;
        self.admit(&next, &mut || checkpoint(gate))?;
        let accepted = Self {
            material: self.material,
            transport: self.transport,
            domain: self.domain,
            field: self.field,
            input: self.input,
            advected,
        };
        Ok((next, accepted, report))
    }

    /// Adaptive SSPRK2 using exactly one coarse step versus two fine steps.
    /// Accepts only the conservative fine trajectory and its physical ledgers;
    /// no extrapolated enthalpy or duplicated mass throughput is published.
    /// The mandatory h norm remains sensitive on latent plateaus and scales
    /// physical changes, not absolute formation enthalpy. Optional temperature
    /// control queries the same original T(h). This estimates local smooth
    /// time error, not spatial error or a certified global error bound.
    pub fn advance_ssprk2_adaptive_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: AdaptiveSsprk2Controls,
        gate: &CancelGate,
    ) -> Result<(DepositionSimulation<'m>, Self, AdaptiveSsprk2Report), AdaptiveSsprk2Error> {
        self.advance_ssprk2_adaptive_with_coefficients(thermal, target_s, controls, unchanged, gate)
    }

    /// Adaptive transport with repeatable stage coefficients under the same
    /// contract as [`Self::advance_ssprk2_with_coefficients`]. All discarded
    /// stages share the caller's native stage/trial budgets. Native events
    /// execute once on the accepted fine timeline; failure/cancellation returns
    /// neither a partial thermal state nor a partially accumulated flow ledger.
    pub fn advance_ssprk2_adaptive_with_coefficients<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: AdaptiveSsprk2Controls,
        mut coefficients: impl FnMut(
            &DepositionSimulation<'m>,
        ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
        gate: &CancelGate,
    ) -> Result<(DepositionSimulation<'m>, Self, AdaptiveSsprk2Report), AdaptiveSsprk2Error> {
        self.admit(thermal, &mut || checkpoint(gate))?;
        let (next, advected, report) = self.advected.advance_ssprk2_adaptive_with_coefficients(
            thermal,
            target_s,
            controls,
            |stage| {
                self.admit(stage, &mut || checkpoint(gate))?;
                coefficients(stage)
            },
            || checkpoint(gate),
        )?;
        self.admit(&next, &mut || checkpoint(gate))?;
        let accepted = Self {
            material: self.material,
            transport: self.transport,
            domain: self.domain,
            field: self.field,
            input: self.input,
            advected,
        };
        Ok((next, accepted, report))
    }
}
