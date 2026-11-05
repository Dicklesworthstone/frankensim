//! Held molecular baths with the original selected gas-surface heat assembly.

use super::*;

impl<'m> MaterialExchangeState<'m> {
    /// Preview held-bath/internal molecular transport and selected gas/body heat.
    ///
    /// Every molecular face shares the same original finite inventory and native
    /// species allowance. The native open owner includes both bath stream signs,
    /// original partial enthalpies, Fourier/Dufour heat and the sensible correction
    /// in one caloric reconstruction. The existing reaction owner then runs once.
    /// The selected law samples the actual transported/reacted composition and gas
    /// midpoint on each joint heat, finite-body, circuit and moving-wall trial.
    /// Its typed source receipts are retained unchanged in `exchange`.
    ///
    /// `resolve` returns total conductance in W/K for an explicitly exposed
    /// `surface`. It must be a pure constitutive evaluation: do not advance
    /// transport or chemistry, deposit heat, or mutate physical state inside it.
    /// The original surface owner admits aggregate area and accounts for each
    /// surface heat transfer once. Bath/internal molecular faces and gas/body
    /// surfaces must denote distinct physical connections.
    ///
    /// This is the same first-order, fixed-geometry, held-bath split as
    /// [`Self::molecular_reservoir_step_candidate`]. Source-relative transport and
    /// original feedback energy defects share one unchanged absolute budget;
    /// actual endpoint species are checked against transport and stoichiometry.
    /// Final source qualification covers every connected endpoint, including
    /// initially zero-flux faces. Refusal publishes no physical owner or clock.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn molecular_reservoir_with_gas_law_candidate<L: TransportLaw, I, E>(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface],
        molecular: &mut ParcelMolecularTransport<L>,
        internal_links: &[ParcelMolecularLink],
        reservoirs: &[ReservoirMolecularLink<'_>],
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<(Self, MolecularReservoirAdvance<I>), MaterialGasLawFlowError<E>> {
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
        let (mut next, exchange) = self.step_with_gas_law_from_chambers(
            &staged, prepared, cx, target_s, ports, drives, reactions, surfaces, controls, guess_j,
            resolve,
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

    /// Publish the complete held-bath, chemistry and selected-heat interval after
    /// all native source, conservation and cancellation gates have succeeded.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_molecular_reservoir_and_gas_law<L: TransportLaw, I, E>(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface],
        molecular: &mut ParcelMolecularTransport<L>,
        internal_links: &[ParcelMolecularLink],
        reservoirs: &[ReservoirMolecularLink<'_>],
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<MolecularReservoirAdvance<I>, MaterialGasLawFlowError<E>> {
        let (next, report) = self.molecular_reservoir_with_gas_law_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            surfaces,
            molecular,
            internal_links,
            reservoirs,
            molecular_controls,
            controls,
            guess_j,
            resolve,
        )?;
        cx.checkpoint()
            .map_err(|_| MaterialGasLawFlowError::Cancelled)?;
        *self = next;
        Ok(report)
    }
}
