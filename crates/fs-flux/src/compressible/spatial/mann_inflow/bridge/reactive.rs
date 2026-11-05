//! One original OU realization through finite-rate chemistry and transport.

use super::*;
use crate::compressible::spatial::mixture::{self as gas, molecular, reactive};
use crate::compressible::spatial::{BodyForce3d, ForceControl3d};
use fs_gas::multicomponent::dissipation::TransportLaw;

/// Optional original held volume force, authored Mach controls and reusable
/// forcing workspace. This carries native data without an additional source,
/// clock, allocation or random draw. Scratch remains owned by the caller.
pub type MannMixtureForce3<'force, 'scratch> = (
    BodyForce3d<'force>,
    ForceControl3d,
    &'scratch mut gas::ForcingWorkspace3d,
);

impl MannReservoirInflow3 {
    /// Evolve chemistry and waves with the same conditional OU boundary path.
    /// Original source-owned kinetics, species references, heat release and
    /// chemical/physical energy accounting remain in `ReactingEuler3d`.
    /// Heat-driven CFL retries replay the retained path at their new times.
    ///
    /// `forcing` optionally composes the original held mechanical source.
    /// Prepare profiles on the Euler field before constructing the reactor.
    /// Bridge work, immutable source and scratch-frame semantics follow
    /// [`Self::step_mixture_bridge`]. This adds velocity fluctuations only.
    #[allow(clippy::too_many_arguments)]
    pub fn step_reacting_bridge(
        &self,
        flow: &mut reactive::ReactingEuler3d,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: reactive::ReactiveStepControl,
        forcing: Option<MannMixtureForce3<'_, '_>>,
        maximum_normal_draws: usize,
        baseline: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<MixtureReservoirFaceState, MixtureError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirBridgeReport3<reactive::ReactiveStepReport3d>, MannReservoirError3>
    {
        self.drive_mixture_bridge(
            bridge,
            frame,
            NativeGeometry::mixture(flow.flow()),
            maximum_duration,
            maximum_normal_draws,
            baseline,
            cancelled,
            |duration, budget, sample, cancel| match forcing {
                Some((force, force_controls, workspace)) => flow
                    .step_with_forcing_and_reservoir_profiles(
                        duration,
                        controls,
                        force,
                        force_controls,
                        workspace,
                        budget,
                        sample,
                        cancel,
                    )
                    .map_err(Into::into),
                None => flow
                    .step_with_reservoir_profiles(duration, controls, budget, sample, cancel)
                    .map_err(Into::into),
            },
        )
    }

    /// Original molecular/chemistry/advection/chemistry/molecular evolution
    /// observes the same path through both SSP diffusion halves and central
    /// waves. Maxwell--Stefan diffusion, enthalpies, full stress and Fourier
    /// heat remain with the supplied original transport operator.
    ///
    /// `maximum_workspace_bytes` covers the native reacting owner, operator
    /// and optional force workspace. Separate atmosphere/profile budgets span
    /// every phase and discarded whole trial. Both physical owners publish
    /// reports/loads only when the complete interval is accepted.
    #[allow(clippy::too_many_arguments)]
    pub fn step_reacting_with_dissipation_bridge<L: TransportLaw>(
        &self,
        flow: &mut reactive::ReactingEuler3d,
        operator: &mut molecular::MolecularTransport3d<L>,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: reactive::ReactiveStepControl,
        dissipation_controls: molecular::DissipationControl,
        maximum_workspace_bytes: usize,
        forcing: Option<MannMixtureForce3<'_, '_>>,
        maximum_normal_draws: usize,
        baseline: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<MixtureReservoirFaceState, MixtureError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirBridgeReport3<reactive::ReactiveStepReport3d>, MannReservoirError3>
    {
        self.drive_mixture_bridge(
            bridge,
            frame,
            NativeGeometry::mixture(flow.flow()),
            maximum_duration,
            maximum_normal_draws,
            baseline,
            cancelled,
            |duration, budget, sample, cancel| match forcing {
                Some((force, force_controls, workspace)) => flow
                    .step_with_dissipation_forcing_and_reservoir_profiles(
                        operator,
                        duration,
                        controls,
                        dissipation_controls,
                        maximum_workspace_bytes,
                        force,
                        force_controls,
                        workspace,
                        budget,
                        sample,
                        cancel,
                    )
                    .map_err(Into::into),
                None => flow
                    .step_with_dissipation_and_reservoir_profiles(
                        operator,
                        duration,
                        controls,
                        dissipation_controls,
                        maximum_workspace_bytes,
                        budget,
                        sample,
                        cancel,
                    )
                    .map_err(Into::into),
            },
        )
    }

    /// RKL2 counterpart with original Legendre stage times and a native
    /// residual budget shared by both molecular halves and chemistry retries.
    /// Enlarged diffusion stability does not establish stochastic temporal
    /// accuracy; native physical step caps remain caller-controlled.
    #[allow(clippy::too_many_arguments)]
    pub fn step_reacting_with_stabilized_dissipation_bridge<L: TransportLaw>(
        &self,
        flow: &mut reactive::ReactingEuler3d,
        operator: &mut molecular::MolecularTransport3d<L>,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: reactive::ReactiveStepControl,
        dissipation_controls: molecular::StabilizedDissipationControl,
        maximum_workspace_bytes: usize,
        forcing: Option<MannMixtureForce3<'_, '_>>,
        maximum_normal_draws: usize,
        baseline: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<MixtureReservoirFaceState, MixtureError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        MannReservoirBridgeReport3<reactive::StabilizedReactiveStepReport3d>,
        MannReservoirError3,
    > {
        self.drive_mixture_bridge(
            bridge,
            frame,
            NativeGeometry::mixture(flow.flow()),
            maximum_duration,
            maximum_normal_draws,
            baseline,
            cancelled,
            |duration, budget, sample, cancel| match forcing {
                Some((force, force_controls, workspace)) => flow
                    .step_with_stabilized_dissipation_forcing_and_reservoir_profiles(
                        operator,
                        duration,
                        controls,
                        dissipation_controls,
                        maximum_workspace_bytes,
                        force,
                        force_controls,
                        workspace,
                        budget,
                        sample,
                        cancel,
                    )
                    .map_err(Into::into),
                None => flow
                    .step_with_stabilized_dissipation_and_reservoir_profiles(
                        operator,
                        duration,
                        controls,
                        dissipation_controls,
                        maximum_workspace_bytes,
                        budget,
                        sample,
                        cancel,
                    )
                    .map_err(Into::into),
            },
        )
    }
}
