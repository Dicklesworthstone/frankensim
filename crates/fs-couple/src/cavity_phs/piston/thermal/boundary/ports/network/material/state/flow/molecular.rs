//! Source-owned continuum molecular transport in the existing complete assembly.
//!
//! The original gas face supplies mass-frame Maxwell--Stefan, pressure and
//! thermal diffusion, Fourier conduction and reciprocal Dufour heat. One native
//! inventory update combines these fluxes with every separately declared path;
//! the original chemistry, circuit, thermal and moving-wall owners then advance.

use super::*;
use fs_gas::material::transport::parcel::{
    ParcelMolecularLink, ParcelMolecularSolve, ParcelMolecularStep, ParcelMolecularTransport,
};
use fs_gas::multicomponent::dissipation::TransportLaw;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    snapshot::MAX_SNAPSHOT_PARCELS, species::network::transport::GasTransportPaths,
};

impl<'m> MaterialExchangeState<'m> {
    /// Preview molecular transport followed by the existing complete assembly.
    ///
    /// Each face resolves its original law at the actual accepted composition,
    /// temperature and pressure. Soret separation and Dufour heat participate
    /// together only when that law supplies explicit cross coefficients. The
    /// molecular mass frame has zero net diffusive mass, even for unequal molar
    /// masses; each original species still carries its own partial enthalpy.
    /// The face's additional heat correction is included exactly once, in the
    /// same caloric reconstruction as the species and all `paths` transfers.
    ///
    /// `links` and `paths` declare distinct physical connections. Do not count
    /// one aperture twice or interpret this zero-bulk-velocity molecular flux
    /// as a resolved jet. Every outgoing species shares one original inventory
    /// allowance. Incoming material can react in this interval but cannot cross
    /// another face until the next interval. Source, entropy, withdrawal and
    /// final native-path driving-force gates remain with their original owners.
    ///
    /// The split is first order: molecular/native transport, reaction, then
    /// joint heat, circuit and piston feedback. It needs time refinement. The
    /// next interval resolves the changed physical state again; no initial
    /// material coefficient or transport temperature is frozen by this adapter.
    ///
    /// Gas inventories, wall motion, electrical/thermal stores, histories and
    /// clock publish together. A refused trial may update the molecular law's
    /// private last-query receipts, as documented by `ParcelMolecularTransport`;
    /// those caches are never accepted gas state. Any species stream or extra
    /// heat path makes the chamber's closed two-port audit unavailable; the
    /// enclosing whole-assembly energy and mass balances remain authoritative.
    #[allow(clippy::too_many_arguments)]
    pub fn molecular_transport_step_candidate<L: TransportLaw>(
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
    ) -> Result<(Self, MaterialFlowAdvance<ParcelMolecularStep>), MaterialFlowError> {
        checkpoint(cx)?;
        let duration_s = target_s - self.time_s();
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || self.chambers.is_empty()
            || self.chambers.len() > controls.maximum_chambers
            || self.chambers.len() > MAX_SNAPSHOT_PARCELS
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
        let (parcels, flow) = molecular
            .step_candidate(
                &parcels,
                links,
                paths,
                duration_s,
                &endpoints,
                molecular_controls,
                || cx.checkpoint().is_err(),
            )?
            .into_parts();
        let mut moved: BTreeSet<_> = flow
            .transport
            .exchange
            .parcels
            .iter()
            .filter(|p| p.incoming_mass_kg > 0.0 || p.outgoing_mass_kg > 0.0)
            .map(|p| p.id)
            .collect();
        // A cell's net heat can vanish between incident links. Retain the
        // actual edges, including cancelling molecular corrections, rather
        // than inferring a closed two-port history from that net balance.
        for link in &flow.links {
            checkpoint(cx)?;
            if link.heat_correction_left_to_right_j != 0.0 {
                moved.extend(link.parcels);
            }
        }
        for heat in &flow.transport.exchange.heat {
            checkpoint(cx)?;
            if heat.heat_left_to_right_j != 0.0 {
                moved.extend(heat.parcels);
            }
        }
        let (next, report) = self.finish_flow_candidate(
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
                mass_tolerance_kg: molecular_controls.physical.exchange.mass_tolerance_kg,
            },
        )?;
        let mut endpoints = BTreeMap::new();
        for (&id, chamber) in &next.chambers {
            checkpoint(cx)?;
            endpoints.insert(id, chamber.gas().clone());
        }
        report
            .flow
            .transport
            .admit_endpoint_driving_forces(&endpoints, || cx.checkpoint().is_err())?;
        Ok((next, report))
    }

    /// Publish all physical owners once after molecular, reaction, feedback,
    /// final driving-force, conservation and cancellation checks have passed.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_molecular_transport<L: TransportLaw>(
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
    ) -> Result<MaterialFlowAdvance<ParcelMolecularStep>, MaterialFlowError> {
        let (next, report) = self.molecular_transport_step_candidate(
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
        )?;
        checkpoint(cx)?;
        *self = next;
        Ok(report)
    }
}
