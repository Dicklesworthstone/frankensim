//! All restrictions sample the same accepted gas state before the existing
//! reaction/electrical/thermal/mechanical assembly advances. Material mixing
//! stays with the shared gas owner; no arrival can leave again in this stage.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::{self as transport, SnapshotNozzleLink, SnapshotNozzleStep},
    snapshot::GasSnapshotSolve,
};

impl<'m> MaterialExchangeState<'m> {
    /// Preview simultaneous nozzle sampling, then the original coupled step.
    ///
    /// Every inlet/outlet reads the ORIGINAL accepted donor's pressure,
    /// temperature, composition and source-bearing enthalpy. Shared transport
    /// aggregates all species and extensive fluxes before reconstructing any
    /// endpoint. Link IDs only label reports; they cannot choose mixing order.
    /// Each chamber's endpoint gas inversion uses its existing drive's
    /// `piston_solve.gas`, so a volume has one caloric admission irrespective of
    /// how many inlets/outlets meet there. Every drive must occur exactly once.
    ///
    /// This is explicit first-order transport followed by reaction and coupled
    /// heat/motion, all over ONE interval. It is not a monolithic flow/reaction
    /// solve. The original finite-volume, quasi-steady nozzle, gross-withdrawal
    /// and joint pressure-crossing restrictions still apply. No pressure-flow
    /// work or formation energy is added again as a heat source.
    ///
    /// Whole-assembly mass/energy and changed-inventory port semantics are the
    /// same as `flow_step_candidate`. All owners remain unchanged on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn snapshot_flow_step_candidate(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        links: &[SnapshotNozzleLink],
        flow_controls: GasSnapshotSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
    ) -> Result<(Self, MaterialFlowAdvance<SnapshotNozzleStep>), MaterialFlowError> {
        checkpoint(cx)?;
        let dt = target_s - self.time_s();
        if !dt.is_finite()
            || dt <= 0.0
            || self.chambers.len() > controls.maximum_chambers
            || drives.len() != self.chambers.len()
        {
            return Err(MaterialFlowError::Invalid(
                "one forward clock and one gas solver per bounded live chamber required",
            ));
        }
        let mut endpoints = BTreeMap::new();
        for drive in drives {
            checkpoint(cx)?;
            if !self.chambers.contains_key(&drive.chamber)
                || endpoints
                    .insert(drive.chamber, drive.piston_solve.gas)
                    .is_some()
            {
                return Err(MaterialFlowError::Invalid(
                    "duplicate drive or unknown live chamber",
                ));
            }
        }
        let mut parcels = BTreeMap::new();
        for (&id, chamber) in &self.chambers {
            checkpoint(cx)?;
            parcels.insert(id, chamber.gas().clone());
        }
        let (parcels, flow) =
            transport::step_candidate(&parcels, links, dt, &endpoints, flow_controls, || {
                cx.checkpoint().is_err()
            })?
            .into_parts();
        let moved = flow
            .exchange
            .parcels
            .iter()
            .filter(|p| p.incoming_mass_kg > 0.0 || p.outgoing_mass_kg > 0.0)
            .map(|p| p.id)
            .collect();
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
                mass_tolerance_kg: flow_controls.mass_tolerance_kg,
            },
        )
    }

    /// Publish all gas, charge/flux, enthalpy, wall motion and port histories
    /// together, after the shared transport and coupled balances both pass.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_snapshot_flow(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        links: &[SnapshotNozzleLink],
        flow_controls: GasSnapshotSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
    ) -> Result<MaterialFlowAdvance<SnapshotNozzleStep>, MaterialFlowError> {
        let (next, report) = self.snapshot_flow_step_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            gas_patches,
            links,
            flow_controls,
            controls,
            guess_j,
        )?;
        checkpoint(cx)?;
        *self = next;
        Ok(report)
    }
}
