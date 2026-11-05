//! Native mechanical forcing and evolving inlet observations share one step.

use super::*;
use crate::compressible::spatial::mixture::{self as gas, molecular as mixture};
use crate::compressible::spatial::molecular as frozen;
use crate::compressible::spatial::{BodyForce3d, ForceControl3d, ForcingWorkspace3d};
use fs_gas::multicomponent::dissipation::TransportLaw;

impl MannReservoirInflow3 {
    /// Compose the original held volume force with one evolving Mann boundary
    /// realization. Native half-kicks retain actual-density impulse/work and
    /// authored Mach/CFL admission. The force has no invented random component.
    /// All bridge, profile, replay and scratch-frame semantics follow
    /// [`Self::step_frozen_bridge`]; the source stays at its accepted endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn step_forced_frozen_bridge(
        &self,
        flow: &mut Euler3d,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: StepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        maximum_normal_draws: usize,
        baseline: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirBridgeReport3<StepReport>, MannReservoirError3> {
        self.drive_frozen_bridge(
            bridge,
            frame,
            NativeGeometry::frozen(flow),
            maximum_duration,
            maximum_normal_draws,
            baseline,
            cancelled,
            |duration, budget, sample, cancel| {
                flow.step_with_forcing_and_reservoir_profiles(
                    duration,
                    controls,
                    force,
                    force_controls,
                    workspace,
                    budget,
                    sample,
                    cancel,
                )
                .map_err(Into::into)
            },
        )
    }

    /// The same original force and OU path on the evolving species owner.
    /// Canonical species and referenced/absolute energy retain native accounts;
    /// the supplied baseline owns composition, pressure and mean velocity.
    #[allow(clippy::too_many_arguments)]
    pub fn step_forced_mixture_bridge(
        &self,
        flow: &mut MixtureEuler3d,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: StepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut gas::ForcingWorkspace3d,
        maximum_normal_draws: usize,
        baseline: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<MixtureReservoirFaceState, MixtureError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirBridgeReport3<gas::StepReport>, MannReservoirError3> {
        self.drive_mixture_bridge(
            bridge,
            frame,
            NativeGeometry::mixture(flow),
            maximum_duration,
            maximum_normal_draws,
            baseline,
            cancelled,
            |duration, budget, sample, cancel| {
                flow.step_with_forcing_and_reservoir_profiles(
                    duration,
                    controls,
                    force,
                    force_controls,
                    workspace,
                    budget,
                    sample,
                    cancel,
                )
                .map_err(Into::into)
            },
        )
    }

    /// Sustained forced viscous flow with the same OU inlet through both native
    /// molecular halves and central force/advection/force. Mechanical source
    /// work remains separate from Fourier, traction and advective exchange.
    #[allow(clippy::too_many_arguments)]
    pub fn step_forced_molecular_frozen_bridge<L: TransportLaw>(
        &self,
        flow: &mut frozen::MolecularEuler3d<L>,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: frozen::MolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        maximum_normal_draws: usize,
        baseline: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirBridgeReport3<frozen::MolecularStepReport3d>, MannReservoirError3>
    {
        self.drive_frozen_bridge(
            bridge,
            frame,
            NativeGeometry::frozen(flow.flow()),
            maximum_duration,
            maximum_normal_draws,
            baseline,
            cancelled,
            |duration, budget, sample, cancel| {
                flow.step_with_forcing_and_reservoir_profiles(
                    duration,
                    controls,
                    force,
                    force_controls,
                    workspace,
                    budget,
                    sample,
                    cancel,
                )
                .map_err(Into::into)
            },
        )
    }

    /// RKL2 molecular counterpart with independent whole-call residual, face
    /// and conditional-normal work bounds. All native physical operators and
    /// endpoint publication remain those of the supplied coupled flow owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_forced_molecular_frozen_bridge<L: TransportLaw>(
        &self,
        flow: &mut frozen::MolecularEuler3d<L>,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: frozen::MolecularStepControl,
        stabilization: frozen::MolecularStabilization,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        maximum_normal_draws: usize,
        baseline: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        MannReservoirBridgeReport3<frozen::StabilizedMolecularStepReport3d>,
        MannReservoirError3,
    > {
        self.drive_frozen_bridge(
            bridge,
            frame,
            NativeGeometry::frozen(flow.flow()),
            maximum_duration,
            maximum_normal_draws,
            baseline,
            cancelled,
            |duration, budget, sample, cancel| {
                flow.step_stabilized_with_forcing_and_reservoir_profiles(
                    duration,
                    controls,
                    stabilization,
                    force,
                    force_controls,
                    workspace,
                    budget,
                    sample,
                    cancel,
                )
                .map_err(Into::into)
            },
        )
    }

    /// Forced Maxwell--Stefan/Fourier/full-stress composition with the original
    /// complete species and enthalpy source, sharing one bounded OU observation
    /// path through accepted and rejected molecular/Euler split attempts.
    #[allow(clippy::too_many_arguments)]
    pub fn step_forced_molecular_mixture_bridge<L: TransportLaw>(
        &self,
        flow: &mut mixture::MolecularEuler3d<L>,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: mixture::MolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut gas::ForcingWorkspace3d,
        maximum_normal_draws: usize,
        baseline: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<MixtureReservoirFaceState, MixtureError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirBridgeReport3<mixture::MolecularStepReport3d>, MannReservoirError3>
    {
        self.drive_mixture_bridge(
            bridge,
            frame,
            NativeGeometry::mixture(flow.flow()),
            maximum_duration,
            maximum_normal_draws,
            baseline,
            cancelled,
            |duration, budget, sample, cancel| {
                flow.step_with_forcing_and_reservoir_profiles(
                    duration,
                    controls,
                    force,
                    force_controls,
                    workspace,
                    budget,
                    sample,
                    cancel,
                )
                .map_err(Into::into)
            },
        )
    }

    /// Original RKL2 diffusion stages and held mechanical forcing consume the
    /// same turbulent bath, with original atomicity and source admission.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_forced_molecular_mixture_bridge<L: TransportLaw>(
        &self,
        flow: &mut mixture::MolecularEuler3d<L>,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: mixture::StabilizedMolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut gas::ForcingWorkspace3d,
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
        MannReservoirBridgeReport3<mixture::StabilizedMolecularStepReport3d>,
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
            |duration, budget, sample, cancel| {
                flow.step_stabilized_with_forcing_and_reservoir_profiles(
                    duration,
                    controls,
                    force,
                    force_controls,
                    workspace,
                    budget,
                    sample,
                    cancel,
                )
                .map_err(Into::into)
            },
        )
    }
}
