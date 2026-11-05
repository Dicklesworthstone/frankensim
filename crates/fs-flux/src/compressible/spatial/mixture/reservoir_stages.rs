//! Complete prescribed exterior gas evaluated at the native SSP stage times.

use super::*;

impl MixtureEuler3d {
    /// Advance with complete reservoir profiles prescribed at the RK stage times.
    ///
    /// Prepare profile ownership with `prepare_reservoir_profiles` first. The
    /// provider receives absolute fluid time [s], the original side/face and its
    /// physical centre [m], and every canonical partial density [kg/m³]. It must
    /// fill all species, including zeros, and return Cartesian velocity [m/s]
    /// and absolute pressure [Pa]. The original gas source admits each state.
    ///
    /// The first residual uses `t`; the second uses `t + dt` for the actual
    /// attempt. A shortened retry resamples both times. The initial `t` sample
    /// also bounds the first interval. Every stage's exterior acoustic rate and
    /// native Riemann wave speeds must pass the combined three-axis CFL gate.
    /// `maximum_face_evaluations` counts every provider call, including rejected
    /// attempts; it is never reset during retries. Budget exhaustion, provider
    /// refusal and exterior-source admission errors are terminal.
    ///
    /// The provider must be pure and replayable at an absolute time. Mutable
    /// closure access supports counters or observation scratch, not advancing a
    /// stochastic realization or external inventory. Such external side effects
    /// cannot be rolled back by this owner. Piecewise/discontinuous time data do
    /// not inherit the smooth-data temporal accuracy of SSP RK2.
    ///
    /// Existing admitted active/trial profile storage is reused without
    /// allocating numerical buffers. Only a completely accepted step publishes
    /// endpoint profiles, cells, time, loads and report together. Ordinary
    /// `step`, forcing/contact options and composed molecular/reacting owners
    /// retain their explicit held-profile semantics.
    pub fn step_with_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        let mut context = Context::new(&mut sample, maximum_face_evaluations);
        self.step_impl_with_contact(
            maximum_duration,
            controls,
            None,
            None,
            Some(&mut context),
            cancelled,
        )
    }
}

type Sample<'a> = dyn FnMut(f64, BoundarySide, usize, [f64; 3], &mut [f64]) -> Result<ReservoirFaceState, Error>
    + 'a;

pub(super) struct Context<'a> {
    sample: &'a mut Sample<'a>,
    maximum_face_evaluations: usize,
    evaluations: usize,
    start_ready: bool,
    failed: bool,
}

impl<'a> Context<'a> {
    pub(super) fn new(sample: &'a mut Sample<'a>, maximum_face_evaluations: usize) -> Self {
        Self {
            sample,
            maximum_face_evaluations,
            evaluations: 0,
            start_ready: false,
            failed: false,
        }
    }

    pub(super) fn failed(&self) -> bool {
        self.failed
    }

    pub(super) fn sample_at(
        &mut self,
        time_s: f64,
        gas: &MixtureGas,
        grid: &CartesianGrid,
        reservoirs: &mut reservoir::Reservoirs,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.start_ready = false;
        let result = (|| {
            if !time_s.is_finite() {
                return Err(Error::Invalid(
                    "finite mixture reservoir sample time required",
                ));
            }
            let profiles = reservoirs.profiles.as_mut().ok_or(Error::Invalid(
                "mixture reservoir profiles are not prepared",
            ))?;
            profiles.stage(
                gas,
                grid,
                &mut |side, face, position, densities| {
                    if self.evaluations >= self.maximum_face_evaluations {
                        return Err(Error::Invalid(
                            "mixture reservoir profile face evaluation budget",
                        ));
                    }
                    self.evaluations += 1;
                    (self.sample)(time_s, side, face, position, densities)
                },
                cancelled,
            )
        })();
        self.failed |= result.is_err();
        result
    }

    pub(super) fn prepare_start(
        &mut self,
        time_s: f64,
        gas: &MixtureGas,
        grid: &CartesianGrid,
        reservoirs: &mut reservoir::Reservoirs,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.sample_at(time_s, gas, grid, reservoirs, cancelled)?;
        self.start_ready = true;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn attempt(
        &mut self,
        time_s: f64,
        gas: &MixtureGas,
        grid: &CartesianGrid,
        boundaries: Boundaries3d,
        reservoirs: &mut reservoir::Reservoirs,
        original: &[f64],
        work: &mut Workspace,
        dt: f64,
        cfl: f64,
        reconstruction: Reconstruction,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<stepping::StageStats, AttemptError> {
        stepping::rk2(cancelled, |final_stage, cancelled| {
            // The initial sample served interval selection as well as the first
            // residual. Once consumed, every retry samples its own two stages.
            if final_stage || !std::mem::take(&mut self.start_ready) {
                self.sample_at(
                    time_s + if final_stage { dt } else { 0.0 },
                    gas,
                    grid,
                    reservoirs,
                    cancelled,
                )
                .map_err(AttemptError::terminal)?;
            }
            let exterior_rate = reservoirs.maximum_rate_at(grid, true);
            stepping::check_cfl::<Error>(exterior_rate * dt, cfl)
                .map_err(AttemptError::retryable)?;
            let mut stats = rk_stage(
                gas,
                grid,
                boundaries,
                reservoirs,
                true,
                original,
                work,
                dt,
                cfl,
                reconstruction,
                None,
                None,
                final_stage,
                cancelled,
            )
            .map_err(AttemptError::retryable)?;
            stats.rate = stats.rate.max(exterior_rate);
            Ok(stats)
        })
    }
}

/// Preserve the native retry policy for numerical trials while making an
/// external provider's refusal terminal. The original error is returned intact.
pub(super) struct AttemptError {
    pub(super) error: Error,
    terminal: bool,
}

impl AttemptError {
    pub(super) fn retryable(error: Error) -> Self {
        Self {
            error,
            terminal: false,
        }
    }

    fn terminal(error: Error) -> Self {
        Self {
            error,
            terminal: true,
        }
    }
}

impl stepping::Failure for AttemptError {
    fn invalid(reason: &'static str) -> Self {
        Self::retryable(Error::Invalid(reason))
    }

    fn cancelled() -> Self {
        Self::terminal(Error::Cancelled)
    }

    fn cfl(observed: f64, maximum: f64) -> Self {
        Self::retryable(Error::CflExceeded { observed, maximum })
    }

    fn retries(attempts: usize, last: Self) -> Self {
        Self::retryable(Error::StepRetries {
            attempts,
            last: Box::new(last.error),
        })
    }

    fn is_cancelled(&self) -> bool {
        // `advance` uses this predicate for immediate return without halving.
        self.terminal || matches!(self.error, Error::Cancelled)
    }

    fn is_cfl(&self) -> bool {
        <Error as stepping::Failure>::is_cfl(&self.error)
    }
}

#[cfg(test)]
#[path = "tests/reservoir_stages.rs"]
mod tests;
