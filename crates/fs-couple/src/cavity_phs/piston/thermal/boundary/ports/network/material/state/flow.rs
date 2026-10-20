//! Gas connectivity on the same live electrical/thermal/reactive assembly.
//! Stage the original conservative nozzle graph once, then call the original
//! assembled reaction/exchange advance. All stages use ONE physical interval.
/// Shared-state nozzle exchange on the same complete physical assembly.
pub mod snapshot;

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::{
    self, GasNetworkError, GasNetworkLink, GasNetworkSolve, GasNetworkStep,
};
use std::collections::BTreeSet;

#[derive(Debug)]
pub enum MaterialFlowError {
    Invalid(&'static str),
    Flow(GasNetworkError),
    SnapshotFlow(network::snapshot::SnapshotNozzleError),
    Exchange(ReactiveExchangeError),
    Gas { chamber: u64, error: fs_material::MaterialError },
    Cancelled,
}
impl core::fmt::Display for MaterialFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material flow assembly: {self:?}")
    }
}
impl std::error::Error for MaterialFlowError {}
impl From<GasNetworkError> for MaterialFlowError {
    fn from(error: GasNetworkError) -> Self {
        match error { GasNetworkError::Cancelled => Self::Cancelled, other => Self::Flow(other) }
    }
}
impl From<network::snapshot::SnapshotNozzleError> for MaterialFlowError {
    fn from(error: network::snapshot::SnapshotNozzleError) -> Self {
        match error {
            network::snapshot::SnapshotNozzleError::Cancelled => Self::Cancelled,
            other => Self::SnapshotFlow(other),
        }
    }
}
impl From<ReactiveExchangeError> for MaterialFlowError {
    fn from(error: ReactiveExchangeError) -> Self {
        match error { ReactiveExchangeError::Cancelled => Self::Cancelled, other => Self::Exchange(other) }
    }
}
/// Full-interval extensive changes, including flow and any subsequent reaction.
#[derive(Clone, Debug, PartialEq)]
pub struct FlowChamberStep {
    pub chamber: u64,
    pub mass_change_kg: f64,
    pub gas_energy_change_j: f64,
    pub entropy_change_j_k: f64,
}
/// Stage reports retain their original reference states. The total balance
/// below encloses ALL live owners from before flow until after exchange/motion.
pub struct MaterialFlowAdvance<F = GasNetworkStep> {
    /// The actual selected transport owner's report, before reaction/heat/motion.
    pub flow: F,
    /// This unchanged report starts AFTER the nozzle graph stage.
    pub exchange: MaterialExchangeAdvance,
    pub chambers: Vec<FlowChamberStep>,
    pub total_mass_residual_kg: f64,
    pub total_energy_residual_j: f64,
}

// Both transport discretizations feed the same whole-interval balance and
// atomic state publication. No second chemistry, circuit or thermal solve.
struct FlowStage<F> {
    parcels: BTreeMap<u64, fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel>,
    moved: BTreeSet<u64>,
    report: F,
    mass_tolerance_kg: f64,
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), MaterialFlowError> {
    cx.checkpoint().map_err(|_| MaterialFlowError::Cancelled)
}
impl<'m> MaterialExchangeState<'m> {
    /// Preview nozzle exchange -> optional reactions -> coupled wall heat and
    /// circuit/momentum evolution. Every stage is the original physical owner;
    /// no material or species registry, caloric law or flow-rate formula is added.
    ///
    /// Link endpoints are live chamber IDs. Branching/cyclic paths use the
    /// original graph's explicit stable-ID Lie splitting; this is first order,
    /// not simultaneous flow/chemistry/heat integration. Both sides of every
    /// path are finite gas stores. To model an outside reservoir it must be an
    /// explicit finite owner; no ambient pressure or composition is invented.
    /// The currently supported stores all have the existing piston dynamics.
    ///
    /// Transfer occurs only once before nonlinear feedback. Carried enthalpy
    /// and formation energy are never deposited again as heat. Source catalogs
    /// must cover any arriving species when reactions are enabled. No vacuum,
    /// flame, pipe storage, jet thrust, check-valve or actuator model is inferred.
    /// Existing pair pressure-crossing, source and withdrawal checks remain.
    ///
    /// Changed-mass chamber port accounts are marked unavailable for closed
    /// two-port auditing because a species-stream port is missing. Actual heat
    /// and mechanical transfers remain observable; use the enclosing balance.
    #[allow(clippy::too_many_arguments)]
    pub fn flow_step_candidate(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], links: &[GasNetworkLink], flow_controls: GasNetworkSolve,
        controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<(Self, MaterialFlowAdvance), MaterialFlowError> {
        use MaterialFlowError as Error;
        checkpoint(cx)?;
        if self.chambers.len() > controls.maximum_chambers || self.chambers.len() > flow_controls.maximum_parcels {
            return Err(Error::Invalid("live chambers exceed flow or coupled work bounds"));
        }
        let dt = target_s - self.time_s();
        let mut parcels = BTreeMap::new();
        for (&id, chamber) in &self.chambers {
            checkpoint(cx)?;
            parcels.insert(id, chamber.gas().clone());
        }
        let (parcels, flow) = network::step_candidate(&parcels, links, dt, flow_controls,
            || cx.checkpoint().is_err())?.into_parts();
        let mut moved = BTreeSet::new();
        for link in &flow.links {
            if link.transfer.transferred_mass_kg != 0.0 {
                moved.insert(link.donor); moved.insert(link.receiver);
            }
        }
        self.finish_flow_candidate(prepared, cx, target_s, ports, drives, reactions,
            gas_patches, controls, guess_j, FlowStage { parcels, moved, report: flow,
                mass_tolerance_kg: flow_controls.mass_tolerance_kg })
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_flow_candidate<F>(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], controls: NetworkSolve, guess_j: &[f64],
        stage: FlowStage<F>,
    ) -> Result<(Self, MaterialFlowAdvance<F>), MaterialFlowError> {
        use MaterialFlowError as Error;
        checkpoint(cx)?;
        let dt = target_s - self.time_s();
        let FlowStage { parcels, moved, report: flow, mass_tolerance_kg } = stage;
        let mut chambers = BTreeMap::new();
        for (id, parcel) in parcels {
            checkpoint(cx)?;
            let mut chamber = self.chambers[&id].clone();
            chamber.gas = parcel;
            chambers.insert(id, chamber);
        }
        let (mut next, exchange) = self.step_from_chambers(&chambers, prepared, cx, target_s, ports, drives,
            reactions, gas_patches, controls, guess_j)?;
        let mut mass = 0.0;
        let mut total = (exchange.report.electrical.stored_after_j - exchange.report.electrical.stored_before_j)
            + exchange.report.thermal.stored_energy_change_j - exchange.report.electrical.native.supplied
            - exchange.step.external_network_heat_j;
        let mut accounts = Vec::with_capacity(self.chambers.len());
        for boundary in &exchange.report.source.chambers {
            checkpoint(cx)?;
            let id = boundary.id;
            let before = &self.chambers[&id];
            let after = next.chambers.get_mut(&id).ok_or(Error::Invalid("flow advance omitted a live chamber"))?;
            let dm = after.gas().mass_kg() - before.gas().mass_kg();
            let u0 = before.gas().state().specific_internal_energy_j_kg();
            let du = after.gas().mass_kg() * (after.gas().state().specific_internal_energy_j_kg() - u0) + dm * u0;
            let s0 = before.gas().state().specific_entropy_j_kg_k().map_err(|error| Error::Gas { chamber: id, error })?;
            let s1 = after.gas().state().specific_entropy_j_kg_k().map_err(|error| Error::Gas { chamber: id, error })?;
            let ds = after.gas().mass_kg() * (s1 - s0) + dm * s0;
            let mechanical = &boundary.boundary.piston;
            let dk = (0.5 * before.mass_kg() * (after.velocity_m_s() - before.velocity_m_s()))
                * (after.velocity_m_s() + before.velocity_m_s());
            let drive = drives.iter().find(|d| d.chamber == id).ok_or(Error::Invalid("missing accepted drive"))?;
            mass += dm;
            total += du + dk + mechanical.exterior_pressure_work_j - mechanical.external_force_work_j - drive.heat_j;
            if ![dm, du, ds, dk, mass, total].iter().all(|v| v.is_finite()) {
                return Err(Error::Invalid("whole-interval mass, entropy or energy is unrepresentable"));
            }
            if moved.contains(&id) {
                // Replace the post-transfer account with ONE whole-interval tick.
                // Even net-zero mass can carry species/enthalpy through this node.
                after.ports = before.ports.clone();
                after.ports.accept(crate::port_adoption::StepAccount {
                    tick: before.ports.next_tick(), duration_s: dt,
                    storage_change_j: Some(du + dk), loss_j: Some(0.0), included_work_j: 0.0,
                    entropy_change_j_k: ds, tolerance_j: Some(controls.total_energy_tolerance_j),
                    unavailable: Some("gas network mass exchange needs stream ports; use the whole-assembly flow balance"),
                }, [crate::port_adoption::PortTransfer::work(
                    mechanical.external_force_work_j - mechanical.exterior_pressure_work_j, true),
                    crate::port_adoption::PortTransfer::heat(mechanical.gas.supplied_heat_j,
                        mechanical.midpoint.state().temperature)]);
            }
            accounts.push(FlowChamberStep { chamber: id, mass_change_kg: dm, gas_energy_change_j: du, entropy_change_j_k: ds });
        }
        if !mass.is_finite() || !total.is_finite() || mass.abs() > mass_tolerance_kg
            || total.abs() > controls.total_energy_tolerance_j {
            return Err(Error::Invalid("flow, reactions, finite bodies, circuit and motion fail whole-interval conservation"));
        }
        checkpoint(cx)?;
        Ok((next, MaterialFlowAdvance { flow, exchange, chambers: accounts,
            total_mass_residual_kg: mass, total_energy_residual_j: total }))
    }

    /// Publish the complete flow/reactive/electrical/thermal/mechanical state
    /// together. Failure at any earlier edge, source or nonlinear trial leaves
    /// all live gas inventories, charge/flux, enthalpy, port history and time intact.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_flow(
        &mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], links: &[GasNetworkLink], flow_controls: GasNetworkSolve,
        controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<MaterialFlowAdvance, MaterialFlowError> {
        let (next, report) = self.flow_step_candidate(prepared, cx, target_s, ports, drives, reactions,
            gas_patches, links, flow_controls, controls, guess_j)?;
        checkpoint(cx)?;
        *self = next;
        Ok(report)
    }
}
