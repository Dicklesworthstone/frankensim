//! Native molecular composition consumes the same bounded bridge sampler.

use super::*;
use crate::compressible::spatial::mixture::molecular as mixture;
use crate::compressible::spatial::molecular as frozen;
use fs_gas::multicomponent::dissipation::TransportLaw;

impl MannReservoirInflow3 {
    /// Observe one OU path through both molecular halves and central Euler,
    /// including all whole-step retries. Molecular bath heat/traction/work use
    /// these same faces when the frozen owner enables reservoir exchange.
    /// The original owner retains its wall, transport law and optional SGS.
    /// Preparation, scratch-frame and budget semantics match
    /// [`Self::step_frozen_bridge`].
    #[allow(clippy::too_many_arguments)]
    pub fn step_molecular_frozen_bridge<L: TransportLaw>(
        &self,
        flow: &mut frozen::MolecularEuler3d<L>,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: frozen::MolecularStepControl,
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
                flow.step_with_reservoir_profiles(duration, controls, budget, sample, cancel)
                    .map_err(Into::into)
            },
        )
    }

    /// The same evolving inlet at the original RKL2 molecular abscissae and
    /// central SSP Euler times. Source normal draws, face/mode samples and
    /// native residual evaluations have independent cumulative bounds.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_molecular_frozen_bridge<L: TransportLaw>(
        &self,
        flow: &mut frozen::MolecularEuler3d<L>,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: frozen::MolecularStepControl,
        stabilization: frozen::MolecularStabilization,
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
                flow.step_stabilized_with_reservoir_profiles(
                    duration,
                    controls,
                    stabilization,
                    budget,
                    sample,
                    cancel,
                )
                .map_err(Into::into)
            },
        )
    }

    /// Original Maxwell--Stefan species transport, partial enthalpy, Fourier
    /// heat and full Newtonian stress observe the same path as central Euler.
    /// The baseline fills the complete canonical composition at each actual
    /// stage; no species, temperature or pressure noise is synthesized here.
    /// Existing native source and wall-inertia refusals remain in force.
    #[allow(clippy::too_many_arguments)]
    pub fn step_molecular_mixture_bridge<L: TransportLaw>(
        &self,
        flow: &mut mixture::MolecularEuler3d<L>,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: mixture::MolecularStepControl,
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
                flow.step_with_reservoir_profiles(duration, controls, budget, sample, cancel)
                    .map_err(Into::into)
            },
        )
    }

    /// RKL2 counterpart of [`Self::step_molecular_mixture_bridge`], retaining
    /// the native operator's boundary quadrature, source laws and atomic split.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_molecular_mixture_bridge<L: TransportLaw>(
        &self,
        flow: &mut mixture::MolecularEuler3d<L>,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: mixture::StabilizedMolecularStepControl,
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
                flow.step_stabilized_with_reservoir_profiles(
                    duration, controls, budget, sample, cancel,
                )
                .map_err(Into::into)
            },
        )
    }
}
