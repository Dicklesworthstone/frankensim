//! Source-owned continuum molecular transport in the existing complete assembly.
//!
//! The original gas face supplies mass-frame Maxwell--Stefan, pressure and
//! thermal diffusion, Fourier conduction and reciprocal Dufour heat. One native
//! inventory update combines these fluxes with every separately declared path;
//! the original chemistry, circuit, thermal and moving-wall owners then advance.

/// Adaptive accuracy over the complete molecular/thermal/electromechanical assembly.
#[path = "molecular_coupled.rs"]
pub mod coupled;
/// Molecular transport with caller-selected gas/body conductance and receipts.
#[path = "molecular_gas_law.rs"]
pub mod gas_law;

/// Explicit held material baths and internal molecular exchange on the live assembly.
#[path = "molecular_reservoir.rs"]
pub mod reservoir;
#[path = "molecular_stage.rs"]
mod stage;

use super::*;
use fs_gas::material::transport::parcel::{
    ParcelMolecularLink, ParcelMolecularSolve, ParcelMolecularStep, ParcelMolecularTransport,
};
use fs_gas::multicomponent::dissipation::TransportLaw;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::GasTransportPaths;

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
    /// Every declared molecular endpoint is requalified after complete feedback,
    /// including an initially zero-flux face; unrelated chambers are not admitted
    /// against a molecular law that does not own their connection.
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
        let stage = self.molecular_stage(
            cx,
            target_s,
            drives,
            molecular,
            links,
            paths,
            molecular_controls,
            controls,
        )?;
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
            stage,
        )?;
        next.admit_molecular_endpoint(cx, molecular, &report.flow, molecular_controls)?;
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
