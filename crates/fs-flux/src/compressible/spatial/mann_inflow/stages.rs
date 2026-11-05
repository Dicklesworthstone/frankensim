//! Observe the same frozen spatial realization at native SSP RK2 abscissae.

use super::*;
use crate::compressible::spatial::{StepControl, StepReport};
use fs_qty::Time;

/// Native accepted-step evidence together with bounded atmosphere work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannReservoirStepReport3<R> {
    pub flow: R,
    /// Includes observations in rejected attempts, not just the accepted pair.
    pub evaluated_faces: usize,
    /// `evaluated_faces * represented_wavevectors`, across the complete call.
    pub mode_samples: usize,
    /// Requested duration clipped by the frozen field's advective phase bound.
    /// The native solver can shorten it further for its actual face CFL.
    pub maximum_duration: Time,
}

// One stack-owned cached trace per stage, shared by every face. Source geometry,
// topology, phase and frame checks do not repeat for each point on a plane.
struct StageSampler<'a> {
    inlet: &'a MannReservoirInflow3,
    source: MannReservoirSource3<'a>,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    maximum_duration: Time,
    maximum_face_evaluations: usize,
    trace: MannReservoirTrace3<'a>,
    evaluated_faces: usize,
    mode_count: usize,
}

impl<'a> StageSampler<'a> {
    fn new(
        inlet: &'a MannReservoirInflow3,
        source: MannReservoirSource3<'a>,
        grid: CartesianGrid,
        boundaries: Boundaries3d,
        time_s: f64,
        maximum_duration: Time,
    ) -> Result<Self, MannReservoirError3> {
        let convection = match source {
            MannReservoirSource3::Frozen {
                convection_velocity_m_s,
                ..
            } => convection_velocity_m_s,
            MannReservoirSource3::Evolving(_) | MannReservoirSource3::BridgeFrame(_) => {
                return Err(MannReservoirError3::Input(
                    "use step_*_bridge with a prepared path for evolving stage-time Mann forcing",
                ));
            }
        };
        let requested = maximum_duration.value();
        if !requested.is_finite() || requested <= 0.0 {
            return Err(MannReservoirError3::Input(
                "positive finite Mann step duration required",
            ));
        }
        let phase_rate = inlet.k_max_per_m * speed(convection);
        if !phase_rate.is_finite() {
            return Err(MannReservoirError3::AdvectiveNyquist {
                phase_bound_rad: phase_rate,
            });
        }
        // Leave a few ulps for the subsequent product checked by trace(). This
        // limits advected spatial modes, not the unbounded frequency of OU noise.
        let duration = if phase_rate > 0.0 {
            requested.min(core::f64::consts::PI * (1.0 - 8.0 * f64::EPSILON) / phase_rate)
        } else {
            requested
        };
        let trace = inlet.trace(&grid, boundaries, time_s, duration, source)?;
        let required = inlet
            .mode_samples
            .checked_mul(2)
            .ok_or(MannReservoirError3::Input("Mann stage work overflow"))?;
        if required > inlet.config.maximum_mode_samples {
            return Err(MannReservoirError3::Work {
                required,
                maximum: inlet.config.maximum_mode_samples,
            });
        }
        let mode_count = inlet.mode_samples / inlet.faces;
        Ok(Self {
            inlet,
            source,
            grid,
            boundaries,
            maximum_duration: Time::new(duration),
            maximum_face_evaluations: inlet.config.maximum_mode_samples / mode_count,
            trace,
            evaluated_faces: 0,
            mode_count,
        })
    }

    fn fluctuation(
        &mut self,
        time_s: f64,
        side: BoundarySide,
        face: usize,
    ) -> Result<[f64; 3], MannReservoirError3> {
        if self.trace.report().sampled_time_s != time_s {
            self.trace = self.inlet.trace(
                &self.grid,
                self.boundaries,
                time_s,
                self.maximum_duration.value(),
                self.source,
            )?;
        }
        // The native profile owner charges every callback before invoking it,
        // including EOS refusals and rejected numerical retries.
        self.evaluated_faces += 1;
        self.trace.fluctuation(side, face)
    }

    fn report<R>(&self, flow: R) -> MannReservoirStepReport3<R> {
        MannReservoirStepReport3 {
            flow,
            evaluated_faces: self.evaluated_faces,
            mode_samples: self.evaluated_faces * self.mode_count,
            maximum_duration: self.maximum_duration,
        }
    }
}

impl MannReservoirInflow3 {
    /// Advance frozen-composition gas with the frozen Mann realization observed
    /// at `t` and each candidate `t + dt`, using native SSP RK2 boundary weights.
    ///
    /// First prepare exactly this inlet's sides with `apply_frozen(Prepare)`.
    /// `baseline` supplies time-dependent mean velocity, density and pressure;
    /// it must be deterministic and replayable, because failed trials resample
    /// earlier times. The atmosphere is borrowed and never advanced. All mode
    /// evaluations, including rejected trials, share `maximum_mode_samples`.
    /// Cells, time, loads, reports and endpoint profiles publish together.
    ///
    /// Only `MannReservoirSource3::Frozen` is admitted. A fixed-step OU endpoint
    /// cannot define values at arbitrary adaptive RK times. The held interface
    /// remains available when the caller owns that stochastic clock.
    pub fn step_frozen<'a>(
        &'a self,
        flow: &mut Euler3d,
        source: MannReservoirSource3<'a>,
        maximum_duration: Time,
        controls: StepControl,
        mut baseline: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
        ) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirStepReport3<StepReport>, MannReservoirError3> {
        poll(&mut cancelled)?;
        if flow.prepared_reservoir_profile_sides() != self.selected {
            return Err(MannReservoirError3::Input(
                "stage-time step requires exactly the prepared Mann sides",
            ));
        }
        let mut sampler = StageSampler::new(
            self,
            source,
            *flow.grid(),
            flow.boundaries(),
            flow.time().value(),
            maximum_duration,
        )?;
        let mut sample_error = None;
        let result = flow.step_with_reservoir_profiles(
            sampler.maximum_duration,
            controls,
            sampler.maximum_face_evaluations,
            |time, side, face, point| {
                let u = sampler.fluctuation(time, side, face).map_err(|error| {
                    sample_error = Some(error);
                    Error::Invalid("Mann stage observation refused")
                })?;
                let mut mean = baseline(time, side, face, point)?;
                for (value, fluctuation) in mean.velocity.iter_mut().zip(u) {
                    *value += fluctuation;
                }
                Ok(mean)
            },
            &mut cancelled,
        );
        match result {
            Ok(report) => Ok(sampler.report(report)),
            Err(error) => Err(sample_error.unwrap_or_else(|| error.into())),
        }
    }

    /// The same stage-time observation for the original multicomponent gas.
    /// `baseline` fills every canonical partial density (including zeros) at
    /// the requested physical face and time. No temperature, composition or
    /// pressure fluctuations are invented by the turbulence adapter.
    ///
    /// Preparation, replay, atomic publication, source-clock and cumulative-work
    /// requirements are identical to [`Self::step_frozen`].
    pub fn step_mixture<'a>(
        &'a self,
        flow: &mut MixtureEuler3d,
        source: MannReservoirSource3<'a>,
        maximum_duration: Time,
        controls: StepControl,
        mut baseline: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<MixtureReservoirFaceState, MixtureError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirStepReport3<super::super::mixture::StepReport>, MannReservoirError3>
    {
        poll(&mut cancelled)?;
        if flow.prepared_reservoir_profile_sides() != self.selected {
            return Err(MannReservoirError3::Input(
                "stage-time step requires exactly the prepared Mann sides",
            ));
        }
        let mut sampler = StageSampler::new(
            self,
            source,
            *flow.grid(),
            flow.boundaries(),
            flow.time().value(),
            maximum_duration,
        )?;
        let mut sample_error = None;
        let result = flow.step_with_reservoir_profiles(
            sampler.maximum_duration,
            controls,
            sampler.maximum_face_evaluations,
            |time, side, face, point, densities: &mut [f64]| {
                let u = sampler.fluctuation(time, side, face).map_err(|error| {
                    sample_error = Some(error);
                    MixtureError::Invalid("Mann stage observation refused")
                })?;
                let mut mean = baseline(time, side, face, point, densities)?;
                for (value, fluctuation) in mean.velocity.iter_mut().zip(u) {
                    *value += fluctuation;
                }
                Ok(mean)
            },
            &mut cancelled,
        );
        match result {
            Ok(report) => Ok(sampler.report(report)),
            Err(error) => Err(sample_error.unwrap_or_else(|| error.into())),
        }
    }
}
