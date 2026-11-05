//! Prescribed exterior histories on the original reacting split transaction.

use super::*;
use crate::compressible::spatial::mixture::{Error, ReservoirFaceState};

impl ReactingEuler3d {
    /// Advance chemistry/waves/chemistry with a prescribed exterior history.
    ///
    /// Prepare reservoir profiles before constructing this owner. The provider
    /// fills partial densities in the original gas's canonical species order
    /// and returns velocity and pressure at each requested physical stage time.
    /// It must replay the same prescribed history after a shortened trial;
    /// sampling must not advance its accepted stochastic/source state.
    ///
    /// One face-evaluation limit covers initial exterior CFL admission, all
    /// original wave stages and every discarded complete chemistry trial.
    /// Provider, exterior admission and work-limit errors are terminal. Cells,
    /// endpoint profiles, chemistry/source/boundary balances, loads and physical
    /// clock publish together after final cancellation admission.
    pub fn step_with_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport3d, ReactiveError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        self.step_impl::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            DissipationControl::default(),
            super::super::MAX_MIXTURE_WORKSPACE_BYTES,
            None,
            None,
            &[],
            Some(&mut context),
            cancelled,
        )
    }

    /// Compose original held body forcing with stage-time reacting reservoirs.
    ///
    /// Force kicks keep the original specific-work, Mach and final acoustic
    /// admission. Boundary transport uses the time-resolved prescribed species,
    /// pressure and velocity; it does not add a second stochastic volume force.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_forcing_and_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport3d, ReactiveError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        self.step_impl::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            DissipationControl::default(),
            super::super::MAX_MIXTURE_WORKSPACE_BYTES,
            Some((force, force_controls, workspace)),
            None,
            &[],
            Some(&mut context),
            cancelled,
        )
    }

    /// Sample one prescribed history through original D/C/A/C/D transport.
    ///
    /// Molecular residuals use physical times in [t,t+dt/2] and [t+dt/2,t+dt];
    /// the central wave residuals use t and t+dt. Source-owned chemistry keeps
    /// both original fixed-volume half intervals and energy references. One
    /// provider budget spans both molecular halves, Euler and all heat-release
    /// retries. Source-owned diffusion, species enthalpy, heat and full stress
    /// remain authoritative. Finite wall inertia and wall acceleration are not
    /// supported by this entry point; prescribed constant walls are retained.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dissipation_and_reservoir_profiles<L: TransportLaw>(
        &mut self,
        operator: &mut MolecularTransport3d<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        maximum_workspace_bytes: usize,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport3d, ReactiveError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        self.step_impl(
            maximum_duration,
            controls,
            Some(operator),
            dissipation_controls,
            maximum_workspace_bytes,
            None,
            None,
            &[],
            Some(&mut context),
            cancelled,
        )
    }

    /// Add original held body work to the timed D/C/A/C/D transaction.
    /// The complete byte allowance includes supplied force workspace; all
    /// physical owners publish only after the entire accepted split succeeds.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dissipation_forcing_and_reservoir_profiles<L: TransportLaw>(
        &mut self,
        operator: &mut MolecularTransport3d<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        maximum_workspace_bytes: usize,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport3d, ReactiveError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        self.step_impl(
            maximum_duration,
            controls,
            Some(operator),
            dissipation_controls,
            maximum_workspace_bytes,
            Some((force, force_controls, workspace)),
            None,
            &[],
            Some(&mut context),
            cancelled,
        )
    }

    /// Use original RKL2 molecular stages with time-resolved reacting inflow.
    ///
    /// Actual molecular evaluations and prescribed-face evaluations each retain
    /// one whole-step budget through both halves and all discarded trials.
    /// Heat-release shortening replays the same physical source interval. The
    /// explicit macrostep bound remains an accuracy choice; stabilization adds
    /// neither a positivity guarantee nor a new chemical/transport model.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_stabilized_dissipation_and_reservoir_profiles<L: TransportLaw>(
        &mut self,
        operator: &mut MolecularTransport3d<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        maximum_workspace_bytes: usize,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveStepReport3d, ReactiveError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        let mut work = StepWork::new(dissipation_controls);
        let coupled = self.step_impl(
            maximum_duration,
            controls,
            Some(operator),
            dissipation_controls.transport,
            maximum_workspace_bytes,
            None,
            Some(&mut work),
            &[],
            Some(&mut context),
            cancelled,
        )?;
        Ok(StabilizedReactiveStepReport3d {
            coupled,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    /// Compose held forcing, finite-rate chemistry and original RKL2 stages
    /// with one prescribed time history and one atomic physical publication.
    /// Finite wall inertia and wall acceleration remain explicitly unsupported.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_stabilized_dissipation_forcing_and_reservoir_profiles<L: TransportLaw>(
        &mut self,
        operator: &mut MolecularTransport3d<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        maximum_workspace_bytes: usize,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveStepReport3d, ReactiveError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        let mut work = StepWork::new(dissipation_controls);
        let coupled = self.step_impl(
            maximum_duration,
            controls,
            Some(operator),
            dissipation_controls.transport,
            maximum_workspace_bytes,
            Some((force, force_controls, workspace)),
            Some(&mut work),
            &[],
            Some(&mut context),
            cancelled,
        )?;
        Ok(StabilizedReactiveStepReport3d {
            coupled,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }
}
