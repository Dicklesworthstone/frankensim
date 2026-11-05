//! Adaptive held-bath transport with original caller-selected gas/body heat.
//!
//! The constitutive callback participates in every complete coarse and fine
//! physical trial. It neither advances transport nor owns a second heat deposit.
use super::*;

/// Preserve the selected law's source error inside the failed physical interval.
pub type CoupledMolecularReservoirGasLawError<E> =
    CoupledPhysicalError<CoupledMolecularReservoirEstimate, MaterialGasLawFlowError<E>>;

/// Accepted native transport and original typed gas-surface source receipts.
pub type CoupledMolecularReservoirGasLawReport<I = ResolvedMaterialStatePoint> =
    CoupledMolecularReservoirReport<I>;

impl<'m> MaterialExchangeState<'m> {
    /// Adapt internal/bath molecular transport, chemistry, selected surface heat,
    /// finite-body/circuit storage and piston motion on one physical window.
    ///
    /// Every original fine step evaluates `resolve` at the actual transported
    /// and reacted gas composition and the gas state of the joint heat/work
    /// trial. Its conductance therefore responds to evolving species, T/p,
    /// finite-body temperatures and electrical heating. The complete state and
    /// each directional molecular/bath stream use the same coarse/fine window.
    /// Only native fine states and their original typed receipts are retained.
    ///
    /// `resolve` returns total conductance [W/K] for the declared surface. It
    /// must be a pure constitutive evaluation, without advancing chemistry,
    /// transport or physical stores. The surface owner deposits its heat once;
    /// molecular Fourier/Dufour heat belongs to distinct gas/gas or bath faces.
    /// Neither source receipts nor callback errors require `Clone` or erasure.
    ///
    /// Original native timestep-sensitive physical failures may refine. A
    /// selected constitutive-source refusal is fatal and remains typed, even
    /// after an earlier interval was privately accepted. A failure publishes
    /// none of the gas, mechanical, thermal, circuit, clock or port history.
    ///
    /// State/stream estimates remain local to the existing first-order split.
    /// Bath states, authored faces/surfaces and external mechanical/electrical
    /// inputs stay held; prescribed gas heat is apportioned at constant power.
    /// Original whole-window conservation allowances are never multiplied by
    /// the number of accepted fine steps or relaxed for a selected heat law.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn coupled_molecular_reservoir_with_gas_law_candidate<L: TransportLaw, I, E>(
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
        time: CoupledMolecularReservoirControls,
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<
        (Self, CoupledMolecularReservoirGasLawReport<I>),
        CoupledMolecularReservoirGasLawError<E>,
    > {
        self.coupled_molecular_reservoir_window_candidate(
            cx,
            target_s,
            ports,
            drives,
            reactions,
            molecular_controls,
            controls,
            guess_j,
            time,
            |state, end, subdrives, guess| {
                state.molecular_reservoir_with_gas_law_candidate(
                    prepared,
                    cx,
                    end,
                    ports,
                    subdrives,
                    reactions,
                    surfaces,
                    molecular,
                    internal_links,
                    reservoirs,
                    molecular_controls,
                    controls,
                    guess,
                    &resolve,
                )
            },
            // The typed selected-law failure owner already classifies the
            // native transport and complete-feedback timestep refusals.
            |_| false,
        )
    }

    /// Publish one complete accepted reservoir window after original physical
    /// balances, directional accuracy and the final cancellation checkpoint.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_coupled_molecular_reservoir_and_gas_law<L: TransportLaw, I, E>(
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
        time: CoupledMolecularReservoirControls,
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<CoupledMolecularReservoirGasLawReport<I>, CoupledMolecularReservoirGasLawError<E>>
    {
        let (next, report) = self.coupled_molecular_reservoir_with_gas_law_candidate(
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
            time,
            resolve,
        )?;
        cx.checkpoint()
            .map_err(|_| CoupledMolecularReservoirGasLawError::Cancelled)?;
        *self = next;
        Ok(report)
    }
}
