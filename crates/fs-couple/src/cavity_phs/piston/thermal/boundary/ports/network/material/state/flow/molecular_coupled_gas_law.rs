//! Adaptive molecular transport through the original selected gas/body heat law.

use super::super::super::diffusion::coupled::CoupledPhysicalError;
use super::*;

/// Complete physical refusals preserve the selected heat law's original error.
pub type CoupledMolecularGasLawError<E> =
    CoupledPhysicalError<CoupledMolecularEstimate, MaterialGasLawFlowError<E>>;

/// Original fine-step molecular, reaction and typed heat-source receipts.
/// Source receipts need not implement `Clone` or any serialization trait.
pub type CoupledMolecularGasLawReport<I = ResolvedMaterialStatePoint> =
    CoupledDiffusionReport<ParcelMolecularStep, I>;

impl<'m> MaterialExchangeState<'m> {
    /// Adapt molecular transport, reactions and caller-selected gas/body heat
    /// through the same complete thermal, electrical and moving-wall trajectory.
    ///
    /// Each physical probe calls the original selected-law interval candidate.
    /// The constitutive callback therefore sees the actual transported/reacted
    /// composition and current heat/work midpoint on every fine step and inner
    /// heat trial. Its typed source receipt travels with that original step;
    /// neither the callback nor an earlier conductance is frozen for the window.
    ///
    /// The shared controller compares every species inventory, gas T/p/volume,
    /// piston velocity, thermal-cell energy, capacitor charge and inductor flux.
    /// Only the two original fine states and receipts are accepted. Independent
    /// whole-window species, energy and momentum accounts use the original
    /// conservation tolerances; refining the time grid does not enlarge them.
    /// The prepared molecular workspace is reused through all physical probes.
    ///
    /// `resolve` is a pure constitutive evaluation returning total conductance
    /// [W/K] on each explicitly declared `surface`. It must not evolve transport
    /// or chemistry, deposit heat, or mutate physical owners. Surface heat and
    /// molecular Fourier/Dufour heat belong to distinct authored connections.
    /// Typed source/domain failures remain fatal; only original native numerical
    /// step-size refusals can retry. No source error is stringified or erased.
    ///
    /// This controls local error of the existing first-order physical split.
    /// Geometry, topology, electrical and mechanical inputs remain held; gas
    /// heat drives are apportioned at constant power over the requested window.
    /// Refused probes may update the molecular cache's last-query receipts, but
    /// publish no live gas, thermal/circuit state, motion, clock or port history.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn coupled_molecular_transport_with_gas_law_candidate<L: TransportLaw, I, E>(
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
        time: CoupledMolecularControls,
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<(Self, CoupledMolecularGasLawReport<I>), CoupledMolecularGasLawError<E>> {
        self.coupled_transport_window_with_error_candidate(
            cx,
            target_s,
            ports,
            drives,
            reactions,
            molecular_controls.physical.exchange,
            controls,
            guess_j,
            time,
            |state, end, subdrives, guess| {
                state.molecular_transport_with_gas_law_candidate(
                    prepared,
                    cx,
                    end,
                    ports,
                    subdrives,
                    reactions,
                    surfaces,
                    molecular,
                    links,
                    paths,
                    molecular_controls,
                    controls,
                    guess,
                    &resolve,
                )
            },
            // The selected-error owner classifies both molecular and feedback
            // refusals. No adapter-specific retry or source fallback is added.
            |_| false,
        )
    }

    /// Publish the complete accepted molecular and selected-heat trajectory
    /// once, after window accuracy, physical balances and final cancellation.
    /// A refusal after an accepted private interval publishes no prefix.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_coupled_molecular_transport_and_gas_law<L: TransportLaw, I, E>(
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
        time: CoupledMolecularControls,
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<CoupledMolecularGasLawReport<I>, CoupledMolecularGasLawError<E>> {
        let (next, report) = self.coupled_molecular_transport_with_gas_law_candidate(
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
            time,
            resolve,
        )?;
        cx.checkpoint()
            .map_err(|_| CoupledMolecularGasLawError::Cancelled)?;
        *self = next;
        Ok(report)
    }
}
