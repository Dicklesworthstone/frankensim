//! Prescribed time histories on the original molecular stage and split owners.

use super::super::{Error, ReservoirFaceState, reservoir_stages::Context};
use super::*;
use crate::compressible::mixture::molecular::stabilized::IntervalFailure;

/// Keep a split phase relative to the original clock. Adding an already
/// rounded half-step origin to another half-step can exceed `t + dt` by one
/// ulp, which would query outside a bounded prescribed stochastic interval.
#[derive(Clone, Copy)]
pub(in crate::compressible::spatial::mixture) struct Clock {
    pub(in crate::compressible::spatial::mixture) origin_s: f64,
    pub(in crate::compressible::spatial::mixture) phase_start_s: f64,
    pub(in crate::compressible::spatial::mixture) phase_end_s: f64,
}

impl Clock {
    pub(in crate::compressible::spatial::mixture) fn direct(origin_s: f64, duration: f64) -> Self {
        Self {
            origin_s,
            phase_start_s: 0.0,
            phase_end_s: duration,
        }
    }

    fn at(self, elapsed: f64, stage_offset: f64, duration: f64) -> f64 {
        let phase_time = if elapsed == duration || stage_offset == duration - elapsed {
            // The controller's exact remaining interval names the accepted
            // phase endpoint; this is not a clamp of external sample times.
            self.phase_end_s
        } else {
            self.phase_start_s + elapsed + stage_offset
        };
        self.origin_s + phase_time
    }
}

/// Borrow the existing profile banks; a stage changes only the private bank.
pub(super) enum Sampling<'a, 'b> {
    Held(&'a Reservoirs),
    Timed {
        reservoirs: &'a mut Reservoirs,
        context: &'a mut Context<'b>,
        clock: Clock,
        duration: f64,
    },
}

impl Sampling<'_, '_> {
    pub(super) fn active(&self) -> &Reservoirs {
        match self {
            Self::Held(reservoirs) => reservoirs,
            Self::Timed { reservoirs, .. } => reservoirs,
        }
    }

    pub(super) fn at(
        &mut self,
        gas: &MixtureGas,
        grid: &CartesianGrid,
        elapsed: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<reservoir::View<'_>, DissipationError> {
        self.at_stage(gas, grid, elapsed, 0.0, cancelled)
    }

    pub(super) fn at_stage(
        &mut self,
        gas: &MixtureGas,
        grid: &CartesianGrid,
        elapsed: f64,
        stage_offset: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<reservoir::View<'_>, DissipationError> {
        match self {
            Self::Held(reservoirs) => Ok(reservoir::View::from(*reservoirs)),
            Self::Timed {
                reservoirs,
                context,
                clock,
                duration,
            } => {
                context.sample_at(
                    clock.at(elapsed, stage_offset, *duration),
                    gas,
                    grid,
                    reservoirs,
                    cancelled,
                )?;
                Ok(reservoir::View {
                    reservoirs,
                    trial: true,
                })
            }
        }
    }

    pub(super) fn failure(&self, error: DissipationError) -> Failure {
        let terminal = match self {
            Self::Held(_) => false,
            Self::Timed { context, .. } => context.failed(),
        };
        Failure { error, terminal }
    }
}

/// Original numerical refusals still halve; an external source refusal does not.
pub(super) struct Failure {
    pub(super) error: DissipationError,
    terminal: bool,
}

impl From<DissipationError> for Failure {
    fn from(error: DissipationError) -> Self {
        Self {
            error,
            terminal: false,
        }
    }
}

impl IntervalFailure for Failure {
    fn is_terminal(&self) -> bool {
        self.terminal || self.error.is_terminal()
    }

    fn step_failure(self, elapsed: f64, dt: f64) -> Self {
        Self {
            error: self.error.step_failure(elapsed, dt),
            terminal: self.terminal,
        }
    }
}

pub(in crate::compressible::spatial::mixture) fn admit(
    flow: &MixtureEuler3d,
    clock: Clock,
    duration: f64,
    finite_inertia: bool,
) -> Result<(), DissipationError> {
    if finite_inertia {
        return Err(DissipationError::Invalid(
            "time-resolved molecular reservoirs require prescribed walls without finite inertia",
        ));
    }
    if flow.reservoirs.profiles.is_none() {
        return Err(DissipationError::Invalid(
            "mixture reservoir profiles are not prepared",
        ));
    }
    let start = clock.origin_s + clock.phase_start_s;
    let end = clock.origin_s + clock.phase_end_s;
    if !clock.origin_s.is_finite()
        || !clock.phase_start_s.is_finite()
        || !clock.phase_end_s.is_finite()
        || clock.phase_start_s < 0.0
        || clock.phase_end_s < clock.phase_start_s
        || !start.is_finite()
        || !duration.is_finite()
        || duration < 0.0
        || !end.is_finite()
        || (duration > 0.0 && end <= start)
    {
        return Err(DissipationError::Invalid(
            "unrepresentable molecular reservoir time interval",
        ));
    }
    Ok(())
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Advance molecular transport with complete reservoirs at each SSP stage.
    ///
    /// `physical_start` supplies the absolute clock because this operator leaves
    /// Euler time unchanged. Every substep and shortened retry samples its own
    /// start and endpoint. The provider fills every canonical partial density,
    /// then returns Cartesian velocity and pressure, admitted by the original
    /// gas source. Prepare profile storage on `flow` before constructing this
    /// operator. The existing active/trial banks are reused without allocation.
    ///
    /// `maximum_face_evaluations` counts every callback, including discarded
    /// trials. Provider, source-admission and budget errors are terminal. The
    /// provider must be pure and replayable at absolute time; external side
    /// effects cannot be rolled back. Smooth data retain the SSP temporal order.
    ///
    /// Cells, endpoint profiles, molecular report and boundary exchange publish
    /// together after the whole interval succeeds. Finite wall inertia and
    /// reactive/adaptive compositions are not provided by this entry point.
    /// A zero duration preserves profiles and makes no provider calls.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_reservoir_profiles(
        &mut self,
        flow: &mut MixtureEuler3d,
        physical_start: Time,
        duration: Time,
        controls: DissipationControl,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        let report = self.advance_driven_with_profiles_unpublished(
            flow,
            duration,
            controls,
            self.maximum_workspace_bytes,
            None,
            None,
            &WallWork::default(),
            Some((
                &mut context,
                Clock::direct(physical_start.value(), duration.value()),
            )),
            &mut cancelled,
        )?;
        self.publish(report, None);
        Ok(report)
    }

    /// Original RKL2 transport with prescribed traces at its Legendre abscissae.
    ///
    /// The origin, provider, work bound and atomic publication have the same
    /// semantics as `advance_with_reservoir_profiles`. Each recurrence residual
    /// samples `physical_start + elapsed + c_j * dt`, including endpoint
    /// admission. The native macrostep cap remains the physical accuracy
    /// control; the enlarged parabolic stability interval is not an error bound.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stabilized_with_reservoir_profiles(
        &mut self,
        flow: &mut MixtureEuler3d,
        physical_start: Time,
        duration: Time,
        controls: StabilizedDissipationControl,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularReport, DissipationError> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        let report = self.advance_stabilized_with_profiles_unpublished(
            flow,
            duration,
            controls,
            self.maximum_workspace_bytes,
            &mut 0,
            Some((
                &mut context,
                Clock::direct(physical_start.value(), duration.value()),
            )),
            &mut cancelled,
        )?;
        self.publish(report.transport, None);
        Ok(report)
    }
}
