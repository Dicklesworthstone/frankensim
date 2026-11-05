//! Whole-assembly accuracy for source-owned molecular species and heat transport.
//!
//! The original coupled controller compares complete coarse and fine physical
//! trajectories. The reusable molecular workspace is the only mutable capture;
//! no accepted inventory, thermal state or mechanical state lives in that cache.

use super::super::diffusion::coupled::{CoupledDiffusionReport, TransportLedger};
use super::*;
use fs_gas::material::transport::parcel::ParcelMolecularError;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkStep, transport::GasTransportStep,
};

pub use super::super::diffusion::coupled::{
    CoupledDiffusionBalance as CoupledMolecularBalance,
    CoupledDiffusionControls as CoupledMolecularControls,
    CoupledDiffusionError as CoupledMolecularError,
    CoupledDiffusionEstimate as CoupledMolecularEstimate,
    CoupledDiffusionScales as CoupledMolecularScales,
};

/// Accepted native molecular face reports, reactions and complete storage accounts.
pub type CoupledMolecularReport = CoupledDiffusionReport<ParcelMolecularStep>;

impl TransportLedger for ParcelMolecularStep {
    fn species_exchange(&self) -> &SpeciesNetworkStep {
        &self.transport.exchange
    }

    fn refinable(error: &MaterialFlowError) -> bool {
        match error {
            MaterialFlowError::Molecular(ParcelMolecularError::Transport(error)) => {
                // Use the existing native-path classification unchanged. Source,
                // geometry and constitutive errors never become timestep retries.
                <GasTransportStep as TransportLedger>::refinable(&MaterialFlowError::Transport(
                    error.clone(),
                ))
            }
            // Feedback can also refuse a native path at the complete final
            // endpoint, outside the molecular owner's nested transport error.
            other => <GasTransportStep as TransportLedger>::refinable(other),
        }
    }
}

impl<'m> MaterialExchangeState<'m> {
    /// Resolve molecular diffusion, optional other gas paths, reactions, wall
    /// heat, circuit storage and moving-chamber feedback on one adaptive window.
    ///
    /// Every complete fine step queries the original molecular law at the
    /// current composition, temperature and pressure. Soret separation, reciprocal
    /// Dufour heat, Fourier conduction and partial enthalpies therefore respond
    /// to the preceding accepted heat/reaction/motion step. The same prepared
    /// face workspace is reused through all coarse, fine and refused trials.
    ///
    /// Only the original two-half-step states and transfers are accepted. Error
    /// scales independently cover each species inventory, gas T/p/volume, piston
    /// velocity, thermal-cell energy, charge and flux. Accepted donor enthalpy
    /// and the molecular heat correction enter the original whole-window balance
    /// once; a nonzero heat-only path also invalidates the closed two-port audit.
    ///
    /// This controls local error in the existing first-order physical split; it
    /// is not a global error certificate or a simultaneous implicit PDE solve.
    /// External drives, path geometry and topology stay fixed during the window.
    /// Source receipts in `molecular` retain the native last-query semantics on
    /// refusal; all live physical owners remain unchanged until publication.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_molecular_transport_candidate<L: TransportLaw>(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        molecular: &mut ParcelMolecularTransport<L>,
        links: &[ParcelMolecularLink],
        paths: GasTransportPaths<'_>,
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        time: CoupledMolecularControls,
    ) -> Result<(Self, CoupledMolecularReport), CoupledMolecularError> {
        self.coupled_transport_window_candidate(
            cx,
            target_s,
            ports,
            drives,
            reactions,
            molecular_controls.physical.exchange,
            controls,
            guess_j,
            time,
            |state, end, subdrives, guess| {
                state.molecular_transport_step_candidate(
                    prepared,
                    cx,
                    end,
                    ports,
                    subdrives,
                    reactions,
                    gas_patches,
                    molecular,
                    links,
                    paths,
                    molecular_controls,
                    controls,
                    guess,
                )
            },
        )
    }

    /// Publish once after complete-window accuracy, physical balances and the
    /// final cancellation checkpoint. A failed later interval publishes no prefix.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_coupled_molecular_transport<L: TransportLaw>(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        molecular: &mut ParcelMolecularTransport<L>,
        links: &[ParcelMolecularLink],
        paths: GasTransportPaths<'_>,
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        time: CoupledMolecularControls,
    ) -> Result<CoupledMolecularReport, CoupledMolecularError> {
        let (next, report) = self.coupled_molecular_transport_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            gas_patches,
            molecular,
            links,
            paths,
            molecular_controls,
            controls,
            guess_j,
            time,
        )?;
        cx.checkpoint()
            .map_err(|_| CoupledMolecularError::Cancelled)?;
        *self = next;
        Ok(report)
    }
}
