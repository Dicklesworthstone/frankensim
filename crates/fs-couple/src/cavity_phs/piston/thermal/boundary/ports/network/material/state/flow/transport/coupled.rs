//! Complete-assembly time accuracy for the original combined transport paths.
//!
//! Every path retains its source, force convention and common finite-inventory
//! owner. Only the ingress differs from molecular diffusion; the complete
//! assembly's adaptive loop, physical error coordinates and accounts are shared.
use super::super::diffusion::coupled::{CoupledDiffusionReport, TransportLedger};
use super::*;
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::SnapshotNozzleError,
    species::{
        diffusion::network::DiffusionNetworkError,
        network::{
            SpeciesNetworkError, SpeciesNetworkStep,
            transport::{membrane::MembraneNetworkError, porous::PorousNetworkError},
        },
    },
};

pub use super::super::diffusion::coupled::{
    CoupledDiffusionBalance as CoupledTransportBalance,
    CoupledDiffusionControls as CoupledTransportControls,
    CoupledDiffusionError as CoupledTransportError,
    CoupledDiffusionEstimate as CoupledTransportEstimate,
    CoupledDiffusionScales as CoupledTransportScales,
};

/// The accepted original combined-path receipts, including layered porous fluxes.
pub type CoupledTransportReport = CoupledDiffusionReport<GasTransportStep>;

impl TransportLedger for GasTransportStep {
    fn species_exchange(&self) -> &SpeciesNetworkStep {
        &self.exchange
    }

    fn refinable(error: &MaterialFlowError) -> bool {
        // Match the original finite-path trajectory's typed, timestep-sensitive
        // refusals. A source/domain/geometry failure remains fatal even when it
        // occurred during a trial; no message matching or tolerance relaxation.
        matches!(
            error,
            MaterialFlowError::Transport(
                GasTransportError::Exchange(
                    SpeciesNetworkError::Withdrawal { .. }
                        | SpeciesNetworkError::Parcel {
                            error: GasEnergyError::OutsideBracket { .. }
                                | GasEnergyError::NotConverged { .. },
                            ..
                        }
                ) | GasTransportError::Nozzle(SnapshotNozzleError::EquilibriumCrossed { .. })
                    | GasTransportError::Porous(PorousNetworkError::EquilibriumCrossed { .. })
                    | GasTransportError::Membrane(MembraneNetworkError::EquilibriumCrossed { .. })
                    | GasTransportError::Diffusion(DiffusionNetworkError::EntropyDecrease { .. })
            )
        )
    }
}

impl<'m> MaterialExchangeState<'m> {
    /// Adapt all declared gas transport paths together with chemistry, heat,
    /// electrical storage and moving-chamber feedback over one held-input window.
    ///
    /// One full coarse step is compared with two complete fine steps through the
    /// ORIGINAL material assembly. Each fine step resamples all nozzle, molecular,
    /// Darcy, selective membrane, homogeneous/layered dusty-gas and finite effusion
    /// paths from its current gas state before their ONE shared inventory update.
    /// Arrivals, reactions, temperature, pressure and volume therefore influence
    /// the following transport solve. Only original fine states are accepted.
    ///
    /// The dimensioned estimate covers every species inventory, gas T/p/volume,
    /// piston velocity, thermal-cell energy, capacitor charge and inductor flux.
    /// The report retains each accepted native path receipt, reaction extent and
    /// whole-window species/energy/momentum balances at the original tolerances.
    /// Coarse, failed and rejected work all spend the same bounded call budget.
    ///
    /// This is local accuracy control of the existing first-order split. Path
    /// geometry, cell masses/topology and external loads remain authored/fixed;
    /// events through the target and external reservoir/vacuum streams require
    /// their own finite owners. A dusty-gas path already contains molecular,
    /// Knudsen and viscous mechanisms for its pore area; do not add them twice.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_transport_paths_candidate(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        paths: GasTransportPaths<'_>,
        transport_controls: GasTransportSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        time: CoupledTransportControls,
    ) -> Result<(Self, CoupledTransportReport), CoupledTransportError> {
        self.coupled_transport_window_candidate(
            cx,
            target_s,
            ports,
            drives,
            reactions,
            transport_controls.exchange,
            controls,
            guess_j,
            time,
            |state, end, subdrives, guess| {
                state.transport_paths_step_candidate(
                    prepared,
                    cx,
                    end,
                    ports,
                    subdrives,
                    reactions,
                    gas_patches,
                    paths,
                    transport_controls,
                    controls,
                    guess,
                )
            },
        )
    }

    /// Publish the complete accepted assembly once, after time accuracy,
    /// whole-window physical balances and the last cancellation checkpoint.
    /// Every failure leaves all live gas, heat, electrical and piston owners intact.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_coupled_transport_paths(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        paths: GasTransportPaths<'_>,
        transport_controls: GasTransportSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        time: CoupledTransportControls,
    ) -> Result<CoupledTransportReport, CoupledTransportError> {
        let (next, report) = self.coupled_transport_paths_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            gas_patches,
            paths,
            transport_controls,
            controls,
            guess_j,
            time,
        )?;
        cx.checkpoint()
            .map_err(|_| CoupledTransportError::Cancelled)?;
        *self = next;
        Ok(report)
    }
}
