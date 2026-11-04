//! Molecular mixing enters the original reacting/electrical/thermal/moving
//! assembly before its joint feedback solve. Every explicit substep samples one
//! accepted gas snapshot; the original owner transports partial enthalpy.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    snapshot::MAX_SNAPSHOT_PARCELS,
    species::diffusion::{GasDiffusionDriving, network::{
        self as transport, DiffusionNetworkLink, DiffusionNetworkSolve, DiffusionNetworkStep,
        adaptive::{self, DiffusionNetworkAdvanceControls, DiffusionNetworkAdvanceReport},
    }},
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};

struct DiffusionSnapshot {
    duration_s: f64,
    parcels: BTreeMap<u64, ClosedGasParcel>,
    endpoints: BTreeMap<u64, GasEnergySolve>,
}

impl MaterialExchangeState<'_> {
    fn diffusion_snapshot(
        &self, cx: &Cx<'_>, target_s: f64, drives: &[ChamberDrive], controls: NetworkSolve,
    ) -> Result<DiffusionSnapshot, MaterialFlowError> {
        checkpoint(cx)?;
        let duration_s = target_s - self.time_s();
        if !duration_s.is_finite() || duration_s <= 0.0 || self.chambers.is_empty()
            || self.chambers.len() > controls.maximum_chambers
            || self.chambers.len() > MAX_SNAPSHOT_PARCELS
            || drives.len() != self.chambers.len() {
            return Err(MaterialFlowError::Invalid(
                "one forward clock and one gas solver per bounded live chamber required"));
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
        Ok(DiffusionSnapshot { duration_s, parcels, endpoints })
    }
}

impl<'m> MaterialExchangeState<'m> {
    /// Preview source-backed Maxwell–Stefan diffusion, then original chemistry,
    /// wall heat, circuit and piston motion over the SAME physical interval.
    ///
    /// This admits counterdiffusion even at equal total pressure. Binary pair
    /// cards and actual face composition/T/p determine the rates; names and gas
    /// viscosity do not substitute for diffusivity. Every face sees the original
    /// accepted inventories. All outgoing species limits apply jointly, before
    /// any incoming material can replenish them. Each endpoint uses its existing
    /// drive's caloric solver, independent of the number or order of faces.
    ///
    /// The original isothermal/isobaric ideal-gas, molar-average-frame closure
    /// remains in force, with explicit admitted T/p differences. It excludes
    /// pressure/thermal diffusion, membrane forces and resolved bulk flow. A
    /// molar-average counterflow need not have zero MASS flow. Do not overlay a
    /// separate nozzle step and assume the combination is a simultaneous solve.
    ///
    /// The split is first order: diffusion, optional reactions, then the existing
    /// coupled heat/motion solve. Arrivals may react in this interval but cannot
    /// diffuse through a second face until the next interval. Partial enthalpy
    /// and chemical formation energy already belong to the gas inventories and
    /// are never deposited again as external heat. Sources must cover arriving
    /// species. Refinement and all original entropy/source/withdrawal gates apply.
    ///
    /// The whole-assembly mass/energy balance encloses the original pre-diffusion
    /// state. Any gross inlet/outlet marks the chamber's closed two-port audit
    /// unavailable, even for zero net mass. Gas inventory, reactions, charge,
    /// enthalpy, motion and port history publish together or remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn diffusion_step_candidate(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], links: &[DiffusionNetworkLink<'_>],
        diffusion_controls: DiffusionNetworkSolve, controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<(Self, MaterialFlowAdvance<DiffusionNetworkStep>), MaterialFlowError> {
        let snapshot = self.diffusion_snapshot(cx, target_s, drives, controls)?;
        let (parcels, flow) = transport::step_candidate(&snapshot.parcels, links, snapshot.duration_s,
            &snapshot.endpoints, diffusion_controls, || cx.checkpoint().is_err())?.into_parts();
        let moved = flow.exchange.parcels.iter()
            .filter(|p| p.incoming_mass_kg > 0.0 || p.outgoing_mass_kg > 0.0)
            .map(|p| p.id).collect();
        self.finish_flow_candidate(prepared, cx, target_s, ports, drives, reactions,
            gas_patches, controls, guess_j, FlowStage { parcels, moved, report: flow,
                mass_tolerance_kg: diffusion_controls.exchange.mass_tolerance_kg })
    }

    /// Advance all live owners once after original diffusion, reaction, coupled
    /// feedback and full-interval conservation have succeeded. Late cancellation
    /// or a source/nonlinear refusal preserves the complete accepted assembly.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_diffusion(
        &mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], links: &[DiffusionNetworkLink<'_>],
        diffusion_controls: DiffusionNetworkSolve, controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<MaterialFlowAdvance<DiffusionNetworkStep>, MaterialFlowError> {
        let (next, report) = self.diffusion_step_candidate(prepared, cx, target_s, ports,
            drives, reactions, gas_patches, links, diffusion_controls, controls, guess_j)?;
        checkpoint(cx)?;
        *self = next;
        Ok(report)
    }

    /// Resolve molecular transport through the complete window with bounded
    /// species-inventory/temperature error, then react and solve heat/motion once.
    ///
    /// Select composition-only molar diffusion or composition-plus-pressure
    /// mass-frame diffusion explicitly. Each accepted half step resamples actual
    /// composition, T/p and the original binary material laws. Arrivals can
    /// continue to another face in later substeps, while every withdrawal gate
    /// applies to one simultaneous snapshot. Rejected probes spend work only.
    ///
    /// The returned report retains gross transfers over ALL accepted substeps,
    /// including reversals. These mark open chamber audits even when the final
    /// net mass or the last substep's transfer vanishes. The existing complete
    /// assembly audit starts before transport, and every owner publishes once.
    ///
    /// Diffusion/reaction/heat-motion splitting is still first order. The local
    /// diffusion error estimate is not a bound on the complete coupled split;
    /// no Soret/Dufour law, bulk flow or simultaneous mechanical residual is
    /// introduced. Original physical/source/entropy budgets remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn adaptive_diffusion_step_candidate(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], links: &[DiffusionNetworkLink<'_>],
        diffusion_controls: DiffusionNetworkSolve, time_controls: DiffusionNetworkAdvanceControls,
        driving: GasDiffusionDriving, controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<(Self, MaterialFlowAdvance<DiffusionNetworkAdvanceReport>), MaterialFlowError> {
        let snapshot = self.diffusion_snapshot(cx, target_s, drives, controls)?;
        let (parcels, flow) = adaptive::advance_candidate_with_driving(
            &snapshot.parcels, links, snapshot.duration_s, &snapshot.endpoints,
            diffusion_controls, time_controls, driving, || cx.checkpoint().is_err(),
        )?.into_parts();
        let mut moved = BTreeSet::new();
        for transfer in &flow.transfers {
            checkpoint(cx)?;
            if transfer.mass_kg > 0.0 {
                moved.insert(transfer.donor);
                moved.insert(transfer.receiver);
            }
        }
        self.finish_flow_candidate(prepared, cx, target_s, ports, drives, reactions,
            gas_patches, controls, guess_j, FlowStage { parcels, moved, report: flow,
                mass_tolerance_kg: diffusion_controls.exchange.mass_tolerance_kg })
    }

    /// Atomically publish error-controlled molecular transport followed by the
    /// original complete material assembly. Work exhaustion, late source failure
    /// or cancellation leaves gas, thermal, electrical and mechanical state intact.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_adaptive_diffusion(
        &mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], links: &[DiffusionNetworkLink<'_>],
        diffusion_controls: DiffusionNetworkSolve, time_controls: DiffusionNetworkAdvanceControls,
        driving: GasDiffusionDriving, controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<MaterialFlowAdvance<DiffusionNetworkAdvanceReport>, MaterialFlowError> {
        let (next, report) = self.adaptive_diffusion_step_candidate(prepared, cx, target_s,
            ports, drives, reactions, gas_patches, links, diffusion_controls, time_controls,
            driving, controls, guess_j)?;
        checkpoint(cx)?;
        *self = next;
        Ok(report)
    }
}
