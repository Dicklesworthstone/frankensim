//! Molecular mixing enters the original reacting/electrical/thermal/moving
//! assembly before its joint feedback solve. All faces read one accepted gas
//! snapshot; the original material owner transports partial species enthalpy.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    snapshot::MAX_SNAPSHOT_PARCELS,
    species::diffusion::network::{
        self as transport, DiffusionNetworkLink, DiffusionNetworkSolve, DiffusionNetworkStep,
    },
};

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
        checkpoint(cx)?;
        let dt = target_s - self.time_s();
        if !dt.is_finite() || dt <= 0.0 || self.chambers.is_empty()
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
        let (parcels, flow) = transport::step_candidate(&parcels, links, dt,
            &endpoints, diffusion_controls, || cx.checkpoint().is_err())?.into_parts();
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
}
