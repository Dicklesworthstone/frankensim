//! Adaptive combined transport through the original selected gas heat assembly.
use super::super::super::diffusion::coupled::CoupledPhysicalError;
use super::*;

/// Original native path receipts with the selected heat law's accepted sources.
pub type CoupledTransportGasLawReport<I> = CoupledDiffusionReport<GasTransportStep, I>;
/// Complete-state accuracy and physical refusals retain the caller's source type.
pub type CoupledTransportGasLawError<E> =
    CoupledPhysicalError<CoupledTransportEstimate, MaterialGasLawFlowError<E>>;

impl<'m> MaterialExchangeState<'m> {
    /// Adapt the complete native transport path set and a selected gas heat law.
    ///
    /// Every coarse and fine trial samples nozzle, molecular, Darcy, selective
    /// membrane, homogeneous/layered dusty-gas and finite effusion paths from ONE
    /// initial inventory. The original shared species owner performs one update,
    /// followed by chemistry, finite-body heat, electrical and moving-wall feedback.
    /// The selected law sees the transported/reacted mixture and actual gas
    /// midpoint during every native heat/work trial; accepted typed sources are
    /// retained without cloning or replacing them with a preset coefficient.
    ///
    /// `resolve` returns total conductance [W/K] for the exposed `surface`. It must
    /// be a pure constitutive evaluation and must not evolve or deposit physics.
    /// Rejected trials may query it repeatedly. All authored path and surface
    /// geometries stay fixed, and each physical connection must be represented
    /// once: a dusty-gas area already contains its molecular/Knudsen/viscous paths.
    ///
    /// The shared complete-state comparison controls every species inventory,
    /// gas T/p/volume, piston velocity, cell energy, charge and flux. Only the
    /// two original fine steps are accepted. Whole-window mass/species, energy
    /// and momentum tolerances retain their original absolute values. Typed
    /// native step-sensitive failures may refine; a selected source refusal is
    /// fatal. This is local control of the original first-order split, with held
    /// external loads and fixed topology; open bath/vacuum streams require their
    /// own complete-window accounting API.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn coupled_transport_paths_with_gas_law_candidate<I, E>(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface],
        paths: GasTransportPaths<'_>,
        transport_controls: GasTransportSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        time: CoupledTransportControls,
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<(Self, CoupledTransportGasLawReport<I>), CoupledTransportGasLawError<E>> {
        self.coupled_transport_window_with_error_candidate(
            cx,
            target_s,
            ports,
            drives,
            reactions,
            transport_controls.exchange,
            controls,
            guess_j,
            time,
            |state, end, subdrives, guess| {
                state.transport_paths_step_with_gas_law_candidate(
                    prepared,
                    cx,
                    end,
                    ports,
                    subdrives,
                    reactions,
                    surfaces,
                    paths,
                    transport_controls,
                    controls,
                    guess,
                    &resolve,
                )
            },
            |_| false,
        )
    }

    /// Publish all accepted gas, thermal, circuit and wall owners together after
    /// whole-window accuracy, physical balances and the final cancellation gate.
    /// A later source or work refusal cannot publish an accepted private prefix.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_coupled_transport_paths_and_gas_law<I, E>(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface],
        paths: GasTransportPaths<'_>,
        transport_controls: GasTransportSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        time: CoupledTransportControls,
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<CoupledTransportGasLawReport<I>, CoupledTransportGasLawError<E>> {
        let (next, report) = self.coupled_transport_paths_with_gas_law_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            surfaces,
            paths,
            transport_controls,
            controls,
            guess_j,
            time,
            resolve,
        )?;
        cx.checkpoint()
            .map_err(|_| CoupledTransportGasLawError::Cancelled)?;
        *self = next;
        Ok(report)
    }
}
