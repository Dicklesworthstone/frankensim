//! Held material baths and internal molecular faces on the live assembly.
//!
//! The original gas owner supplies every species/enthalpy/heat transfer and
//! its open entropy audit. The existing reaction, finite-wall, circuit and
//! piston owners then consume exactly those gas inventories. Their independently
//! checked energy increments compose without a new subtraction of gross gas
//! storage and bath enthalpy. Transport retains its original source-reference
//! cancellation; feedback retains its original caloric/reaction limits.

use super::*;
use fs_gas::material::transport::parcel::reservoir::{
    ReservoirMolecularLink, ReservoirMolecularStep,
};

#[path = "molecular_reservoir_balance.rs"]
mod balance;
pub use balance::{
    MolecularReservoirBalance, MolecularReservoirChamberBalance, MolecularReservoirSpeciesBalance,
};

/// Original open transport and feedback receipts, plus their shared physical budget.
pub struct MolecularReservoirAdvance {
    pub flow: ReservoirMolecularStep,
    pub exchange: MaterialExchangeAdvance,
    pub balance: MolecularReservoirBalance,
}

impl<'m> MaterialExchangeState<'m> {
    /// Preview simultaneous internal/bath molecular transport, then full feedback.
    ///
    /// Every declared face samples one original inventory and spends one native
    /// gross species allowance. Each bath has explicit held composition, T/p and
    /// original material sources. Both stream signs carry that bath's enthalpy;
    /// Fourier/Dufour heat and the native sensible correction enter one caloric
    /// reconstruction. No bath species or heat is applied again during feedback.
    ///
    /// The returned balance combines the native reference-cancelled transport
    /// residual and original feedback residual under one unchanged absolute budget.
    /// Actual original/final species inventories are checked separately against
    /// internal transfers, bath streams and original reaction stoichiometry.
    /// A large transported formation reference does not erase the separately
    /// retained feedback increment. Reaction calorics retain their original limits.
    ///
    /// This is the original first-order split with fixed authored faces and held
    /// baths. The native bath owner admits internal molecular faces; bulk-flow,
    /// porous, membrane and effusion paths need their own simultaneous owner.
    /// Geometry and bath states are not advanced. Refine the timestep to resolve
    /// coupled feedback; the closed-system adaptive API does not admit bath streams.
    ///
    /// Final molecular source qualification covers every connected endpoint even
    /// when its initial flux was zero. Live gas, motion, thermal/circuit stores,
    /// port histories and clock publish together only after the final checkpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn molecular_reservoir_step_candidate<L: TransportLaw>(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        molecular: &mut ParcelMolecularTransport<L>,
        internal_links: &[ParcelMolecularLink],
        reservoirs: &[ReservoirMolecularLink<'_>],
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
    ) -> Result<(Self, MolecularReservoirAdvance), MaterialFlowError> {
        let inputs = super::stage::inputs(self, cx, target_s, drives, controls)?;
        let (parcels, flow) = molecular
            .network_candidate(
                &inputs.parcels,
                internal_links,
                reservoirs,
                inputs.duration_s,
                &inputs.endpoints,
                molecular_controls,
                || cx.checkpoint().is_err(),
            )?
            .into_parts();
        let mut staged = BTreeMap::new();
        for (id, parcel) in parcels {
            checkpoint(cx)?;
            let mut chamber = self
                .chambers
                .get(&id)
                .ok_or(MaterialFlowError::Invalid(
                    "reservoir transport returned an unknown chamber",
                ))?
                .clone();
            chamber.gas = parcel;
            staged.insert(id, chamber);
        }
        let (mut next, exchange) = self.step_from_chambers(
            &staged,
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            gas_patches,
            controls,
            guess_j,
        )?;
        let balance = balance::finish(
            self,
            &staged,
            &mut next,
            cx,
            &flow,
            &exchange,
            reactions,
            molecular_controls,
            controls,
        )?;
        let connected: BTreeSet<_> = flow
            .internal_links
            .iter()
            .flat_map(|link| link.parcels)
            .chain(flow.links.iter().map(|link| link.parcel))
            .collect();
        for id in connected {
            checkpoint(cx)?;
            let gas = next.chambers.get(&id).ok_or(MaterialFlowError::Invalid(
                "missing final reservoir-connected molecular owner",
            ))?;
            molecular.qualify_parcel(gas.gas(), molecular_controls.face, || {
                cx.checkpoint().is_err()
            })?;
        }
        checkpoint(cx)?;
        Ok((
            next,
            MolecularReservoirAdvance {
                flow,
                exchange,
                balance,
            },
        ))
    }

    /// Publish one complete held-bath molecular/thermal/electromechanical interval.
    /// Refusal at any transport or later feedback stage publishes no physical prefix.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_molecular_reservoir<L: TransportLaw>(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        molecular: &mut ParcelMolecularTransport<L>,
        internal_links: &[ParcelMolecularLink],
        reservoirs: &[ReservoirMolecularLink<'_>],
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
    ) -> Result<MolecularReservoirAdvance, MaterialFlowError> {
        let (next, report) = self.molecular_reservoir_step_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            gas_patches,
            molecular,
            internal_links,
            reservoirs,
            molecular_controls,
            controls,
            guess_j,
        )?;
        checkpoint(cx)?;
        *self = next;
        Ok(report)
    }
}
