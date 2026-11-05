//! Nonautonomous observations on the original molecular transport owners.

use super::*;
use crate::compressible::mixture::molecular::stabilized::IntervalFailure;
use fs_gas::multicomponent::dissipation::DissipationError;

pub(super) use super::super::reservoir::TimedSource;

/// Keep the original physical origin through split phases. Adding two rounded
/// half-step clocks can leave the admitted outer interval by one ulp; offsets
/// are combined before the sole addition to the original origin instead.
#[derive(Clone, Copy)]
pub(super) struct Clock {
    origin: f64,
    phase_start: f64,
    phase_end: f64,
}

impl Clock {
    pub(super) fn interval(origin: f64, duration: f64) -> Self {
        Self {
            origin,
            phase_start: 0.0,
            phase_end: duration,
        }
    }

    pub(super) fn split(origin: f64, duration: f64, second_half: bool) -> Self {
        Self {
            origin,
            phase_start: if second_half { 0.5 * duration } else { 0.0 },
            phase_end: if second_half {
                duration
            } else {
                0.5 * duration
            },
        }
    }

    pub(super) fn at(self, elapsed: f64, duration: f64) -> f64 {
        self.origin
            + if elapsed == duration {
                self.phase_end
            } else {
                self.phase_start + elapsed
            }
    }

    pub(super) fn admit(self, duration: f64) -> Result<(), MolecularError> {
        let start = self.at(0.0, duration);
        let end = self.at(duration, duration);
        if !start.is_finite() || !end.is_finite() || end <= start {
            return Err(MolecularError::Invalid(
                "unrepresentable timed molecular interval",
            ));
        }
        Ok(())
    }
}

pub(super) fn admit_interval(start: f64, duration: f64) -> Result<(), MolecularError> {
    if !start.is_finite() || !(start + duration).is_finite() || start + duration <= start {
        return Err(MolecularError::Invalid(
            "unrepresentable timed molecular interval",
        ));
    }
    Ok(())
}

/// One complete source admission before an original spatial evaluation. The
/// returned immutable view cannot outlive the borrow of its reusable cache.
#[allow(clippy::too_many_arguments)]
pub(super) fn stage<'a>(
    reservoirs: &'a mut super::super::reservoir::Reservoirs,
    source: Option<&mut TimedSource<'_>>,
    gas: &FrozenGas,
    grid: &CartesianGrid,
    time: f64,
    exchange: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<super::super::reservoir::View<'a>>, MolecularError> {
    let timed = source.is_some();
    if let Some(source) = source {
        source.stage(&mut reservoirs.profiles, gas, grid, time, cancelled)?;
    }
    Ok(exchange.then(|| {
        if timed {
            reservoirs.staged()
        } else {
            reservoirs.accepted()
        }
    }))
}

/// The generic RKL controller retains its numerical retry policy while a
/// provider, exterior-source or profile-budget refusal keeps its original
/// error and terminates immediately. No public error category is reclassified.
pub(super) struct AttemptFailure {
    pub(super) error: MolecularError,
    pub(super) source_failed: bool,
}

impl From<DissipationError> for AttemptFailure {
    fn from(error: DissipationError) -> Self {
        Self {
            error: error.into(),
            source_failed: false,
        }
    }
}

impl IntervalFailure for AttemptFailure {
    fn is_terminal(&self) -> bool {
        self.source_failed || self.error.is_terminal()
    }

    fn step_failure(self, elapsed: f64, dt: f64) -> Self {
        Self {
            error: self.error.step_failure(elapsed, dt),
            ..self
        }
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Integrate prescribed bath traces at the actual nonautonomous SSP times.
    ///
    /// `physical_start` is explicit because a molecular operator does not own
    /// or advance the Euler clock. Each substep observes its start and endpoint;
    /// rejected trials resample the reduced interval. Molecular exchange still
    /// requires [`Self::with_reservoir_exchange`]. Original caloric admission,
    /// full Newtonian stress/work, Fourier heat and selected SGS are unchanged.
    ///
    /// Prepare the flow's profile geometry first. The provider has the same
    /// replayable, read-only contract as [`Euler3d::step_with_reservoir_profiles`]:
    /// identical arguments must return identical physical data, including when
    /// retry times move backwards. One face budget covers the complete duration
    /// and every discarded trial; source and budget refusals are terminal.
    ///
    /// Reusable active/trial caches avoid per-step allocations. Only a complete
    /// interval publishes cells, molecular loads and its endpoint profile.
    /// Failure preserves all accepted data; Euler time and step metadata remain
    /// unchanged even on success. A following held step uses the new endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_reservoir_profiles(
        &mut self,
        flow: &mut Euler3d,
        physical_start: Time,
        duration: Time,
        controls: MolecularControl,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularReport3d, MolecularError> {
        let mut source = TimedSource::new(&mut sample, maximum_face_evaluations);
        self.advance_impl(
            flow,
            duration,
            controls,
            Clock::interval(physical_start.value(), duration.value()),
            Some(&mut source),
            cancelled,
        )
    }

    /// Integrate a replayable prescribed bath at the original RKL2 abscissae.
    ///
    /// The anchor is observed once at the macrostep start; each Legendre stage
    /// observes its own physical time and the final admission observes the
    /// endpoint. A rejected macrostep rebuilds its anchor and all visited bath
    /// states. This retains native RKL2 coefficients and signed field/boundary
    /// weights, with separate cumulative operator and provider work budgets.
    /// Physical origins, source ownership and atomic publication follow
    /// [`Self::advance_with_reservoir_profiles`]. `maximum_step_s` still controls
    /// time resolution; stability acceleration is not a time-error estimate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stabilized_with_reservoir_profiles(
        &mut self,
        flow: &mut Euler3d,
        physical_start: Time,
        duration: Time,
        controls: MolecularControl,
        stabilization: MolecularStabilization,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularReport3d, MolecularError> {
        let mut source = TimedSource::new(&mut sample, maximum_face_evaluations);
        self.advance_stabilized_with_budget(
            flow,
            duration,
            controls,
            stabilization,
            &mut 0,
            Clock::interval(physical_start.value(), duration.value()),
            Some(&mut source),
            &mut cancelled,
        )
    }
}
