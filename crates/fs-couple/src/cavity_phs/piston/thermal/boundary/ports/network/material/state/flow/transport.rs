//! Both existing transport laws feed ONE inventory before the original coupled step.
//! No accepted gas, circuit, thermal or mechanical history changes during a trial.
use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
    nozzle::network::snapshot::SnapshotNozzleLink,
    species::{diffusion::network::DiffusionNetworkLink,
        network::transport::{self as shared, GasTransportError, GasTransportSolve, GasTransportStep}},
};

impl From<GasTransportError> for MaterialFlowError {
    fn from(error: GasTransportError) -> Self {
        match error { GasTransportError::Cancelled => Self::Cancelled, other => Self::Transport(other) }
    }
}

impl<'m> MaterialExchangeState<'m> {
    /// Preview pressure restrictions and molecular diffusion on one accepted
    /// snapshot, then react the combined inventories and run the ORIGINAL joint
    /// electrical/thermal/moving-wall solve. All stages use ONE physical interval.
    ///
    /// Lists declare distinct physical paths, with their own source-backed laws;
    /// the adapter does not add a nozzle velocity to a molar-average diffusion
    /// velocity on the same aperture. Diffusion retains its explicit admitted
    /// temperature/pressure differences and source domains. Heating or motion
    /// can invalidate that admission on the next interval; nothing is clamped.
    ///
    /// All outgoing species across BOTH mechanisms share the same original
    /// inventory allowance. Incoming material cannot replenish that allowance
    /// or traverse a second path in this explicit stage. The combined transport
    /// state is reconstructed once per volume and chemistry runs once per bound
    /// chamber. Nonlinear heat/electrical/mechanical trials never repeat either.
    ///
    /// Partial enthalpy and formation energy are already included in the gas;
    /// neither is added again as heat. Existing whole-assembly mass/energy checks
    /// enclose the original pre-transport gas, circuit, thermal and kinetic state.
    /// Any GROSS species transfer makes the closed two-port chamber audit
    /// unavailable even at zero net mass, while retaining the enclosing account.
    ///
    /// This is explicit first-order transport then reaction and joint feedback,
    /// not implicit CFD, stiff kinetics, shared crank motion or jet thrust.
    #[allow(clippy::too_many_arguments)]
    pub fn transport_step_candidate(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], nozzles: &[SnapshotNozzleLink],
        diffusion: &[DiffusionNetworkLink<'_>], transport_controls: GasTransportSolve,
        controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<(Self, MaterialFlowAdvance<GasTransportStep>), MaterialFlowError> {
        checkpoint(cx)?;
        let dt = target_s - self.time_s();
        if !dt.is_finite() || dt <= 0.0 || self.chambers.is_empty()
            || self.chambers.len() > controls.maximum_chambers
            || self.chambers.len() > MAX_SNAPSHOT_PARCELS
            || drives.len() != self.chambers.len()
            || nozzles.len().checked_add(diffusion.len()).is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS) {
            return Err(MaterialFlowError::Invalid(
                "one forward clock, bounded combined paths and one drive per live chamber required"));
        }
        let mut endpoints = BTreeMap::new();
        for drive in drives {
            checkpoint(cx)?;
            if !self.chambers.contains_key(&drive.chamber)
                || endpoints.insert(drive.chamber, drive.piston_solve.gas).is_some() {
                return Err(MaterialFlowError::Invalid("duplicate drive or unknown live chamber"));
            }
        }
        let mut parcels = BTreeMap::new();
        for (&id, chamber) in &self.chambers {
            checkpoint(cx)?;
            parcels.insert(id, chamber.gas().clone());
        }
        let (parcels, flow) = shared::step_candidate(&parcels, nozzles, diffusion, dt,
            &endpoints, transport_controls, || cx.checkpoint().is_err())?.into_parts();
        let moved = flow.exchange.parcels.iter()
            .filter(|p| p.incoming_mass_kg > 0.0 || p.outgoing_mass_kg > 0.0)
            .map(|p| p.id).collect();
        self.finish_flow_candidate(prepared, cx, target_s, ports, drives, reactions,
            gas_patches, controls, guess_j, FlowStage { parcels, moved, report: flow,
                mass_tolerance_kg: transport_controls.exchange.mass_tolerance_kg })
    }

    /// Publish all gas, charge/flux, enthalpy, wall motion and port histories
    /// together, once, after combined transport and all coupled balances pass.
    /// Failed source, budget, nonlinear or final-cancellation checks leave the
    /// original time and every physical owner intact for a retry.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_transport(
        &mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], nozzles: &[SnapshotNozzleLink],
        diffusion: &[DiffusionNetworkLink<'_>], transport_controls: GasTransportSolve,
        controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<MaterialFlowAdvance<GasTransportStep>, MaterialFlowError> {
        let (next, report) = self.transport_step_candidate(prepared, cx, target_s, ports,
            drives, reactions, gas_patches, nozzles, diffusion, transport_controls, controls, guess_j)?;
        checkpoint(cx)?;
        *self = next;
        Ok(report)
    }
}
