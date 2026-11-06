//! Exponential binary transport on the original complete material assembly.
//!
//! The finite gas owner resolves source-dependent mixing and partial enthalpy.
//! The existing flow admission then encloses diffusion, reactions, thermal and
//! electrical feedback, and piston motion in one original-to-final account.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::network::binary::{
    self as transport_binary, BinaryDiffusionNetworkAdvanceControls,
    BinaryDiffusionNetworkAdvanceReport,
};

impl<'m> MaterialExchangeState<'m> {
    /// Preview adaptive exponential binary diffusion and the full live assembly.
    ///
    /// Every open edge requires the original binary molar-frame source closure.
    /// Actual finite volumes, compositions, temperatures and pair diffusivities
    /// refresh at each pair update. The gas owner checks species and temperature
    /// error, charges rejected/failed pair work and accumulates actual accepted
    /// species enthalpy transfers. Native source and withdrawal limits still
    /// apply; removing the explicit diffusion stability limit does not remove
    /// the physical model's admission gates.
    ///
    /// All accepted gross transfers mark the corresponding chamber audits open,
    /// including equimolar counterdiffusion with zero net mass. Diffusion and
    /// reaction enthalpy already belong to the actual gas inventories and are
    /// never deposited again as external heat. The existing assembly balance
    /// starts before diffusion, and no owner or accepted history changes until
    /// the complete feedback solve and whole-window mass/energy audit succeed.
    ///
    /// This remains a first-order diffusion/reaction/heat-motion split. The gas
    /// error estimate controls transport alone, not the complete coupled split.
    /// Its binary molar frame does not supply pressure/thermal diffusion, nozzle
    /// flow or simultaneous moving-volume transport. Cut intervals at authored
    /// forcing changes and refine the complete assembly interval as needed.
    #[allow(clippy::too_many_arguments)]
    pub fn adaptive_binary_diffusion_step_candidate(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        links: &[DiffusionNetworkLink<'_>],
        diffusion_controls: DiffusionNetworkSolve,
        time_controls: BinaryDiffusionNetworkAdvanceControls,
        controls: NetworkSolve,
        guess_j: &[f64],
    ) -> Result<
        (
            Self,
            MaterialFlowAdvance<BinaryDiffusionNetworkAdvanceReport>,
        ),
        MaterialFlowError,
    > {
        let snapshot = self.diffusion_snapshot(cx, target_s, drives, controls)?;
        let (parcels, flow) = transport_binary::advance_adaptive_candidate(
            &snapshot.parcels,
            links,
            snapshot.duration_s,
            &snapshot.endpoints,
            diffusion_controls,
            time_controls,
            || cx.checkpoint().is_err(),
        )?
        .into_parts();
        let mut moved = BTreeSet::new();
        for transfer in &flow.transfers {
            checkpoint(cx)?;
            if transfer.mass_kg > 0.0 {
                moved.insert(transfer.donor);
                moved.insert(transfer.receiver);
            }
        }
        self.finish_flow_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            gas_patches,
            controls,
            guess_j,
            FlowStage {
                parcels,
                moved,
                report: flow,
                mass_tolerance_kg: diffusion_controls.exchange.mass_tolerance_kg,
            },
        )
    }

    /// Publish adaptive binary transport and every coupled owner atomically.
    ///
    /// Pair-work exhaustion, missing sources, a late feedback refusal or
    /// cancellation leaves the accepted gas, chemistry, circuit, thermal field,
    /// piston state and physical clock unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_adaptive_binary_diffusion(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        links: &[DiffusionNetworkLink<'_>],
        diffusion_controls: DiffusionNetworkSolve,
        time_controls: BinaryDiffusionNetworkAdvanceControls,
        controls: NetworkSolve,
        guess_j: &[f64],
    ) -> Result<MaterialFlowAdvance<BinaryDiffusionNetworkAdvanceReport>, MaterialFlowError> {
        let (next, report) = self.adaptive_binary_diffusion_step_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            gas_patches,
            links,
            diffusion_controls,
            time_controls,
            controls,
            guess_j,
        )?;
        checkpoint(cx)?;
        *self = next;
        Ok(report)
    }
}
