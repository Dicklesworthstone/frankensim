//! Replayable prescribed exterior states at the two nonautonomous SSP times.

use super::*;

type Sample<'a> =
    dyn FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error> + 'a;

/// Only the provider and its remaining work budget are retained on the stack.
/// The admitted active/trial arrays remain the sole numerical face storage.
pub(in crate::compressible::spatial) struct TimedSource<'a> {
    sample: &'a mut Sample<'a>,
    remaining: usize,
    terminal_failure: bool,
}

impl TimedSource<'_> {
    pub(in crate::compressible::spatial) fn failed(&self) -> bool {
        self.terminal_failure
    }

    pub(in crate::compressible::spatial) fn stage(
        &mut self,
        profiles: &mut Profiles,
        gas: &FrozenGas,
        grid: &CartesianGrid,
        time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let result = (|| {
            poll(cancelled)?;
            if !time.is_finite() {
                return Err(Error::Invalid("nonfinite reservoir profile sample time"));
            }
            if profiles.values.is_empty() {
                return Err(Error::Invalid(
                    "reservoir profiles must be prepared before timed stepping",
                ));
            }
            // Refuse a scan which cannot fit before invoking any of its faces.
            // Earlier scans, including all rejected trials, consume the same
            // budget. A callback which refuses still counts as an evaluation.
            if profiles.values.len() > self.remaining {
                return Err(Error::Invalid(
                    "timed reservoir profile face evaluation budget exhausted",
                ));
            }
            profiles.stage(
                gas,
                grid,
                |side, face, point| {
                    self.remaining -= 1;
                    (self.sample)(time, side, face, point)
                },
                cancelled,
            )
        })();
        if result.is_err() {
            // Provider, exterior-source and work-budget refusals are terminal.
            // They are not a numerical cell/CFL failure that halving can cure.
            self.terminal_failure = true;
        }
        result
    }
}

impl Euler3d {
    /// Advance with complete exterior face data at the actual SSP stage times.
    ///
    /// Prepare profile geometry with [`Self::prepare_reservoir_profiles`] first.
    /// `sample(time_s, side, face, position_m)` supplies all prepared faces in
    /// canonical side/face order at the accepted start time and at `time + dt`.
    /// Source-admitted endpoint data participate in the unsplit CFL check and
    /// the second residual. A smaller retry samples its own start and endpoint
    /// times again. Integrated boundary exchanges retain both RK weights.
    ///
    /// The callback must be a deterministic, pure observation of a read-only,
    /// replayable source at the requested physical time. It must not advance a
    /// random stream, an OU clock, an external solver or irreversible state.
    /// Interior mutability may be used for diagnostics, but must not change the
    /// returned physical state for identical arguments. Neither a sample-hold
    /// approximation nor an interpolation of previously accepted profiles is
    /// substituted for a requested stage observation.
    ///
    /// `maximum_face_evaluations` bounds all callback invocations across this
    /// call, including rejected attempts; each full face scan must fit before
    /// it starts. Provider, gas-source and budget refusals terminate the call.
    /// Numerical cell/CFL failures retain the native bounded retry policy.
    /// The existing two profile arrays are reused without allocation or cloning.
    /// A failed or cancelled call preserves accepted profiles, cells, clock,
    /// report and loads. A successful call publishes its endpoint profile with
    /// the physical step, so subsequent held stepping uses that endpoint.
    ///
    /// This method advances native Euler only. It does not implicitly retime
    /// body-force, molecular, reacting or evolving stochastic-source owners.
    pub fn step_with_reservoir_profiles(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        maximum_face_evaluations: usize,
        mut sample: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        let mut source = TimedSource {
            sample: &mut sample,
            remaining: maximum_face_evaluations,
            terminal_failure: false,
        };
        self.step_impl(
            maximum_duration,
            controls,
            None,
            Some(&mut source),
            cancelled,
        )
    }
}

#[cfg(test)]
mod tests;
