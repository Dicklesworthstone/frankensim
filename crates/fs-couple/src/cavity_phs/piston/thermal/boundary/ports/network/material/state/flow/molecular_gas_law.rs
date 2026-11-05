//! Molecular species/heat transport with the original selected gas-surface law.

use super::*;

impl<'m> MaterialExchangeState<'m> {
    /// Preview molecular transport and the caller-selected gas heat assembly.
    ///
    /// All molecular and independently declared native paths share the same
    /// original finite inventory and one caloric reconstruction. The existing
    /// reaction owner then runs once. The selected heat law receives the actual
    /// transported/reacted gas composition and energy midpoint on each joint
    /// heat, circuit, finite-body and wall-motion trial. Its typed source receipt
    /// is retained unchanged in the original exchange report.
    ///
    /// `resolve` returns total conductance in W/K for an explicitly exposed
    /// `surface`. It must be a pure constitutive evaluation: it must not advance
    /// transport or chemistry, deposit heat, or mutate physical state. The
    /// original surface owner admits aggregate area and accounts for each heat
    /// transfer once. Molecular Dufour/Fourier heat belongs to the declared gas
    /// links; selected surface heat belongs to the distinct gas/body patches.
    ///
    /// This remains the first-order split described by
    /// [`Self::molecular_transport_step_candidate`]. Native final driving forces
    /// and every declared molecular endpoint's original source/law are admitted
    /// again after complete feedback. A refusal publishes no physical owner.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn molecular_transport_with_gas_law_candidate<L: TransportLaw, I, E>(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface],
        molecular: &mut ParcelMolecularTransport<L>,
        links: &[ParcelMolecularLink],
        paths: GasTransportPaths<'_>,
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<(Self, MaterialFlowAdvance<ParcelMolecularStep, I>), MaterialGasLawFlowError<E>>
    {
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
        let (next, report) = self.finish_flow_with_gas_law_candidate(
            prepared, cx, target_s, ports, drives, reactions, surfaces, controls, guess_j, stage,
            resolve,
        )?;
        next.admit_molecular_endpoint(cx, molecular, &report.flow, molecular_controls)?;
        Ok((next, report))
    }

    /// Publish molecular, chemistry, selected surface heat and all feedback
    /// owners together after original source, balance and cancellation gates.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_molecular_transport_and_gas_law<L: TransportLaw, I, E>(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface],
        molecular: &mut ParcelMolecularTransport<L>,
        links: &[ParcelMolecularLink],
        paths: GasTransportPaths<'_>,
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<MaterialFlowAdvance<ParcelMolecularStep, I>, MaterialGasLawFlowError<E>> {
        let (next, report) = self.molecular_transport_with_gas_law_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            surfaces,
            molecular,
            links,
            paths,
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
