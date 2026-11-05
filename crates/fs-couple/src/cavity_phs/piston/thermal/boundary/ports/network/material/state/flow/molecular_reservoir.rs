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
/// Held-bath molecular transport with caller-selected gas/body heat laws.
#[path = "molecular_reservoir_gas_law.rs"]
pub mod gas_law;
#[path = "molecular_reservoir_stage.rs"]
mod stage;
pub use balance::{
    MolecularReservoirBalance, MolecularReservoirChamberBalance, MolecularReservoirSpeciesBalance,
};

/// Original open transport and feedback receipts, plus their shared physical budget.
pub struct MolecularReservoirAdvance<I = ResolvedMaterialStatePoint> {
    pub flow: ReservoirMolecularStep,
    pub exchange: MaterialExchangeAdvance<I>,
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
        let (staged, flow) = self.molecular_reservoir_stage(
            cx,
            target_s,
            drives,
            molecular,
            internal_links,
            reservoirs,
            molecular_controls,
            controls,
        )?;
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
        next.admit_molecular_reservoir_endpoint(cx, molecular, &flow, molecular_controls)?;
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
