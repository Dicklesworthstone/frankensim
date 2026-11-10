//! Physical-time error control on the original source-owned mixed SDIRK2 owner.

use super::{
    MaterialDissipatedHeat, MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxSystem,
    MaterialFluxViscosityTangent, MaterialNavierStokesBoundary, MaterialNavierStokesSdirk2Step,
};
use crate::ns::boundary::OpenBoundaryPolicy;
use crate::ns::material::heat::DissipatedHeatAccumulator;
use crate::ns::newton::{FluxNewtonControls, FluxNewtonError, FluxSdirk2Error, FluxSdirk2Stage};
use fs_exec::CancelGate;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_time::PiController;

#[path = "adaptive/heat_control.rs"]
mod heat_control;
use heat_control::HeatControl;
pub use heat_control::{
    MaterialNavierStokesHeatAdaptiveAdvance, MaterialSdirk2HeatControls, MaterialSdirk2HeatReport,
    MaterialSdirk2HeatTolerance,
};

/// Original physical edge laws evaluated at the actual absolute stage time [s].
/// Topology and policy are fixed; data may vary smoothly in time and position.
/// Velocities are m/s, total tractions Pa and positions metres.
#[derive(Clone, Copy)]
pub enum MaterialNavierStokesTimeBoundary<'a> {
    Velocity {
        edge: usize,
        value: &'a dyn Fn(f64, [f64; 2]) -> [f64; 2],
    },
    NavierSlip {
        edge: usize,
        value: &'a dyn Fn(f64, [f64; 2]) -> [f64; 2],
        /// Spatially constant on this edge, in Pa s/m; finite and nonnegative.
        tangential_resistance_pa_s_per_m: &'a dyn Fn(f64) -> f64,
    },
    Open {
        edge: usize,
        traction: &'a dyn Fn(f64, [f64; 2], [f64; 2]) -> [f64; 2],
        backflow: OpenBoundaryPolicy,
    },
    OpenReservoir {
        edge: usize,
        traction: &'a dyn Fn(f64, [f64; 2], [f64; 2]) -> [f64; 2],
        velocity: &'a dyn Fn(f64, [f64; 2]) -> [f64; 2],
    },
}

/// Bounded step-doubling controls. Step sizes refer to a complete comparison;
/// the retained fine solution uses two half intervals.
#[derive(Clone, Debug)]
pub struct MaterialSdirk2AdaptiveControls {
    pub initial_step_s: f64,
    /// Floor for an unclipped comparison. A shorter final interval is allowed
    /// only when its actual stages are representable and its error is admitted.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    /// Positive absolute tolerance on physical, area-weighted RMS velocity.
    pub velocity_absolute_tolerance_m_s: f64,
    pub velocity_relative_tolerance: f64,
    /// Includes ALL rejected and accepted comparisons, at most six Newton
    /// stages each. Every stage also retains its original Newton/Krylov limits.
    pub maximum_attempts: usize,
    pub newton: FluxNewtonControls,
    /// Shared controller settings; defaults use order-three local-error gains.
    pub controller: PiController,
    /// Previous accepted normalized error, in (0,1], for resumed PI proposals.
    pub previous_error_ratio: f64,
}

impl Default for MaterialSdirk2AdaptiveControls {
    fn default() -> Self {
        Self {
            initial_step_s: 1e-3,
            minimum_step_s: 1e-12,
            maximum_step_s: 0.1,
            velocity_absolute_tolerance_m_s: 1e-6,
            velocity_relative_tolerance: 1e-3,
            maximum_attempts: 1000,
            newton: FluxNewtonControls::default(),
            controller: PiController {
                k_p: 0.7 / 3.0,
                k_i: 0.4 / 3.0,
                safety: 0.9,
                max_growth: 2.0,
                max_shrink: 0.2,
            },
            previous_error_ratio: 1.0,
        }
    }
}

/// Which original SDIRK2 solve supplied a refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaterialSdirk2Trial {
    Coarse,
    FirstHalf,
    SecondHalf,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MaterialSdirk2StageFailure {
    pub trial: MaterialSdirk2Trial,
    /// Actual absolute time of the refused stage, or the subinterval start
    /// when the original SDIRK2 owner rejects input before either stage.
    pub time_s: f64,
    pub error: FluxSdirk2Error<MaterialFluxError>,
}

/// No partial advancement or altered input history accompanies an error.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialSdirk2AdaptiveError {
    Invalid(&'static str),
    Cancelled,
    Stage(MaterialSdirk2StageFailure),
    /// Accepted-history heat could not be represented; no partial advancement
    /// or accumulated heat is published. This is terminal, not a retry signal.
    Heat {
        time_s: f64,
        error: MaterialFluxError,
    },
    WorkLimit {
        attempts: usize,
        last_numerical_failure: Option<MaterialSdirk2StageFailure>,
    },
    StepUnderflow {
        time_s: f64,
        timestep_s: f64,
    },
    Tolerance {
        time_s: f64,
        timestep_s: f64,
        normalized_error: f64,
    },
}
impl core::fmt::Display for MaterialSdirk2AdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive mixed material SDIRK2: {self:?}")
    }
}
impl std::error::Error for MaterialSdirk2AdaptiveError {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialSdirk2AdaptiveReport {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub attempts: usize,
    pub accepted_comparisons: usize,
    pub rejected_comparisons: usize,
    pub numerical_rejections: usize,
    pub stage_solves: usize,
    pub maximum_accepted_error: f64,
    pub last_step_start_time_s: f64,
    pub last_error_ratio: f64,
    /// Unclipped PI proposal, suitable for the next call's initial step.
    pub next_step_s: f64,
}

/// Only the last ACCEPTED fine comparison is retained, with its original
/// material owners and physical reactions. Earlier states are not archived;
/// their accepted original-stage heat is accumulated in bounded storage.
pub struct MaterialNavierStokesAdaptiveAdvance<G> {
    pub first_half: MaterialNavierStokesSdirk2Step<G>,
    /// Its endpoint is the final result at report.end_time_s.
    pub second_half: MaterialNavierStokesSdirk2Step<G>,
    pub report: MaterialSdirk2AdaptiveReport,
    /// Whole-interval bulk and slip production [J/m], summed from accepted
    /// fine half steps only. Coarse and rejected candidates contribute nothing.
    pub dissipated_heat: MaterialDissipatedHeat,
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Advance to an exact physical endpoint using source-owned mixed SDIRK2
    /// step doubling. Force, boundary and material callbacks receive actual
    /// absolute stage times, including retries and both half intervals. They
    /// must be pure functions of time and their declared spatial/state inputs.
    /// Force density is N/m³, velocity m/s, time s and material tangents Pa s².
    ///
    /// Each comparison solves one full step and two half steps from the SAME
    /// accepted history. The fine endpoint is accepted without Richardson
    /// extrapolation. The error indicator is the physical RMS velocity norm of
    /// (u_fine-u_coarse)/3, divided by atol + rtol*max(RMS(u_old),RMS(u_fine)).
    /// RMS uses the original finite-element kinetic-energy quadrature divided
    /// by mesh area, never the norm of representation-dependent edge moments.
    /// Pressure is not part of this temporal error estimate.
    ///
    /// Original material sources, fixed density, mixed boundary policies,
    /// analytic Newton operators and physical reactions remain authoritative.
    /// Source/boundary/input refusals terminate the whole call. Only numerical
    /// Newton line-search, correction-solve and convergence failures may retry
    /// with a smaller step. The attempt budget bounds all six stages per trial,
    /// including rejected trials; per-stage Newton/Krylov limits remain active.
    ///
    /// Cancellation is checked before each two-stage native solve and before
    /// publication. The synchronous native Newton solver is not interrupted
    /// inside that pair. Caller history is immutable; only complete successful
    /// advancement returns a result. Storage is O(mesh size), with one reused
    /// difference vector and only the last accepted fine pair retained.
    ///
    /// Positive-duration forward intervals are required. This is local time
    /// error control under resolved spatial/nonlinear error and smooth data,
    /// not a global error certificate or adaptive adjoint. Accepted original
    /// stage powers are integrated into separate bulk and slip heat over the
    /// WHOLE requested interval. Coarse/rejected heat is discarded, and no
    /// source is queried again or past source archive retained. Velocity error
    /// control does not certify heat accuracy: two strongly damped endpoints
    /// can agree while underresolving fast-mode dissipation history. Stiff
    /// accuracy does not imply nonlinear positivity, passivity, caloric
    /// evolution or an integrated energy balance.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_navier_stokes_mixed_rheology_sdirk2<G: FluidTransport>(
        &self,
        start_time_s: f64,
        end_time_s: f64,
        previous_velocity: &[f64],
        force_density: &dyn Fn(f64, [f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesTimeBoundary<'_>],
        controls: &MaterialSdirk2AdaptiveControls,
        resolve: impl FnMut(
            f64,
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
        gate: &CancelGate,
    ) -> Result<MaterialNavierStokesAdaptiveAdvance<G>, MaterialSdirk2AdaptiveError> {
        self.advance_mixed_sdirk2_inner(
            start_time_s,
            end_time_s,
            previous_velocity,
            force_density,
            boundaries,
            controls,
            None,
            resolve,
            gate,
        )
        .map(|(flow, _)| flow)
    }

    /// Advance with additional original cell/edge heat error control and an
    /// accepted-history budget on absolute signed SDIRK2 kinetic time defects.
    ///
    /// This uses the SAME six native Newton stages per comparison, source
    /// callbacks, absolute stage times, accepted fine history and work bounds
    /// as the velocity-only entry. Physical heat is never corrected or replaced
    /// by the time defect. The original PI settings see the largest enabled
    /// indicator; conservative order-three gains also bound heat-relative and
    /// duration-allocated defect criteria whose ratios can scale as dt².
    ///
    /// Absolute defect credit is allocated by comparison dt/requested duration;
    /// relative credit uses only that comparison's original generated heat.
    /// This prevents an early stiff transient consuming the whole absolute
    /// allowance, at the cost of additional resolution in long requested windows.
    /// These are local temporal indicators and a time-quadrature account, not
    /// a pressure/stress, spatial-error, thermal-evolution or passivity certificate.
    ///
    /// # Errors
    /// All velocity-only refusals remain transactional. Invalid heat controls
    /// refuse before callbacks. Failed heat/defect tolerance uses the same
    /// bounded rejection path; no partial state, heat or defect account returns.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_navier_stokes_mixed_rheology_sdirk2_with_heat_control<G: FluidTransport>(
        &self,
        start_time_s: f64,
        end_time_s: f64,
        previous_velocity: &[f64],
        force_density: &dyn Fn(f64, [f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesTimeBoundary<'_>],
        controls: &MaterialSdirk2AdaptiveControls,
        heat_controls: &MaterialSdirk2HeatControls,
        resolve: impl FnMut(
            f64,
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
        gate: &CancelGate,
    ) -> Result<MaterialNavierStokesHeatAdaptiveAdvance<G>, MaterialSdirk2AdaptiveError> {
        let (flow, heat_control) = self.advance_mixed_sdirk2_inner(
            start_time_s,
            end_time_s,
            previous_velocity,
            force_density,
            boundaries,
            controls,
            Some(heat_controls),
            resolve,
            gate,
        )?;
        Ok(MaterialNavierStokesHeatAdaptiveAdvance {
            flow,
            heat_control: heat_control.ok_or(MaterialSdirk2AdaptiveError::Invalid(
                "missing admitted heat-control account",
            ))?,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_mixed_sdirk2_inner<G: FluidTransport>(
        &self,
        start_time_s: f64,
        end_time_s: f64,
        previous_velocity: &[f64],
        force_density: &dyn Fn(f64, [f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesTimeBoundary<'_>],
        controls: &MaterialSdirk2AdaptiveControls,
        heat_controls: Option<&MaterialSdirk2HeatControls>,
        mut resolve: impl FnMut(
            f64,
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
        gate: &CancelGate,
    ) -> Result<
        (
            MaterialNavierStokesAdaptiveAdvance<G>,
            Option<MaterialSdirk2HeatReport>,
        ),
        MaterialSdirk2AdaptiveError,
    > {
        poll(gate)?;
        controls.admit(start_time_s, end_time_s)?;
        let mut heat_control = heat_controls
            .map(|c| HeatControl::new(*c, end_time_s - start_time_s))
            .transpose()?;
        self.system
            .sdirk2_stage_timestep(previous_velocity, controls.initial_step_s)
            .map_err(|error| {
                MaterialSdirk2AdaptiveError::Stage(MaterialSdirk2StageFailure {
                    trial: MaterialSdirk2Trial::Coarse,
                    time_s: start_time_s,
                    error: FluxSdirk2Error {
                        stage: None,
                        error: MaterialFluxError::Newton(error),
                    },
                })
            })?;
        let area: f64 = self.system.mesh.areas.iter().sum();
        if !area.is_finite() || area <= 0.0 {
            return Err(MaterialSdirk2AdaptiveError::Invalid(
                "finite positive mesh area required",
            ));
        }
        let mut history = previous_velocity.to_vec();
        let mut history_rms = self.physical_rms_velocity(&history, area)?;
        let mut difference = vec![0.0; self.system.n_u];
        let mut time = start_time_s;
        let mut proposal = controls.initial_step_s;
        let mut error_previous = controls.previous_error_ratio;
        let mut last_failure = None;
        let mut accepted = None;
        let mut heat = DissipatedHeatAccumulator::new(self.system.mesh.areas.len());
        let mut report = MaterialSdirk2AdaptiveReport {
            start_time_s,
            end_time_s,
            attempts: 0,
            accepted_comparisons: 0,
            rejected_comparisons: 0,
            numerical_rejections: 0,
            stage_solves: 0,
            maximum_accepted_error: 0.0,
            last_step_start_time_s: start_time_s,
            last_error_ratio: error_previous,
            next_step_s: proposal,
        };
        while time < end_time_s {
            poll(gate)?;
            if report.attempts == controls.maximum_attempts {
                return Err(MaterialSdirk2AdaptiveError::WorkLimit {
                    attempts: report.attempts,
                    last_numerical_failure: last_failure,
                });
            }
            let remaining = end_time_s - time;
            let clipped = proposal > remaining;
            let end = if proposal >= remaining {
                end_time_s
            } else {
                time + proposal
            };
            let dt = end - time;
            let midpoint = time + 0.5 * dt;
            stage_times(time, midpoint)?;
            stage_times(midpoint, end)?;
            stage_times(time, end)?;
            report.attempts += 1;
            let comparison = (|| {
                let mut step = |trial, start, end, previous: &[f64]| {
                    self.timed_mixed_sdirk2(
                        trial,
                        start,
                        end,
                        previous,
                        force_density,
                        boundaries,
                        controls.newton,
                        &mut resolve,
                        &mut report.stage_solves,
                        gate,
                    )
                };
                let coarse = step(MaterialSdirk2Trial::Coarse, time, end, &history)?;
                let first = step(MaterialSdirk2Trial::FirstHalf, time, midpoint, &history)?;
                let middle = &first.endpoint.flow.step.solution.discrete().x[..self.system.n_u];
                let second = step(MaterialSdirk2Trial::SecondHalf, midpoint, end, middle)?;
                Ok::<_, MaterialSdirk2AdaptiveError>((coarse, first, second))
            })();
            let (coarse, first, second) = match comparison {
                Ok(value) => value,
                Err(MaterialSdirk2AdaptiveError::Stage(failure)) if failure.retryable() => {
                    report.rejected_comparisons += 1;
                    report.numerical_rejections += 1;
                    let next = 0.5 * dt;
                    if next < controls.minimum_step_s || next >= dt || time + next <= time {
                        return Err(MaterialSdirk2AdaptiveError::Stage(failure));
                    }
                    last_failure = Some(failure);
                    proposal = next;
                    continue;
                }
                Err(error) => return Err(error),
            };
            let fine = &second.endpoint.flow.step.solution.discrete().x[..self.system.n_u];
            let coarse_velocity = &coarse.endpoint.flow.step.solution.discrete().x;
            for i in 0..difference.len() {
                difference[i] = (fine[i] - coarse_velocity[i]) / 3.0;
            }
            let fine_rms = self.physical_rms_velocity(fine, area)?;
            let scale = controls.velocity_absolute_tolerance_m_s
                + controls.velocity_relative_tolerance * history_rms.max(fine_rms);
            let error = self.physical_rms_velocity(&difference, area)? / scale;
            if !scale.is_finite() || !error.is_finite() {
                return Err(MaterialSdirk2AdaptiveError::Invalid(
                    "unrepresentable velocity error scale",
                ));
            }
            let heat_trial = heat_control
                .as_ref()
                .map(|control| control.trial(&coarse, &first, &second, dt, gate))
                .transpose()?;
            let error = heat_trial
                .as_ref()
                .map_or(error, |heat| error.max(heat.error_ratio()));
            let factor = if error.is_infinite() {
                controls.controller.max_shrink
            } else {
                controls.factor(error, error_previous)?
            };
            let suggested = (dt * factor).clamp(controls.minimum_step_s, controls.maximum_step_s);
            if error > 1.0 {
                report.rejected_comparisons += 1;
                if dt <= controls.minimum_step_s || suggested >= dt || time + suggested <= time {
                    return Err(MaterialSdirk2AdaptiveError::Tolerance {
                        time_s: time,
                        timestep_s: dt,
                        normalized_error: error,
                    });
                }
                proposal = suggested;
                last_failure = None;
                continue;
            }
            poll(gate)?;
            for half in [&first, &second] {
                heat.add_integral(&half.dissipated_heat)
                    .map_err(|error| MaterialSdirk2AdaptiveError::Heat { time_s: end, error })?;
            }
            if let (Some(control), Some(trial)) = (&mut heat_control, heat_trial) {
                control.accept(trial)?;
            }
            history.copy_from_slice(fine);
            history_rms = fine_rms;
            report.accepted_comparisons += 1;
            report.maximum_accepted_error = report.maximum_accepted_error.max(error);
            report.last_step_start_time_s = time;
            error_previous = error.max(1e-300);
            report.last_error_ratio = error_previous;
            time = end;
            if !clipped {
                proposal = suggested;
            }
            report.next_step_s = proposal;
            last_failure = None;
            accepted = Some((first, second));
        }
        poll(gate)?;
        let (first_half, second_half) = accepted.ok_or(MaterialSdirk2AdaptiveError::Invalid(
            "positive duration did not produce an accepted interval",
        ))?;
        let dissipated_heat = heat.finish(end_time_s - start_time_s).map_err(|error| {
            MaterialSdirk2AdaptiveError::Heat {
                time_s: end_time_s,
                error,
            }
        })?;
        let heat_control = heat_control.map(HeatControl::finish).transpose()?;
        poll(gate)?;
        Ok((
            MaterialNavierStokesAdaptiveAdvance {
                first_half,
                second_half,
                report,
                dissipated_heat,
            },
            heat_control,
        ))
    }

    fn physical_rms_velocity(
        &self,
        velocity: &[f64],
        area: f64,
    ) -> Result<f64, MaterialSdirk2AdaptiveError> {
        let squared = 2.0 * (self.system.kinetic_energy(velocity) / area);
        if !squared.is_finite() || squared < 0.0 {
            return Err(MaterialSdirk2AdaptiveError::Invalid(
                "unrepresentable physical velocity norm",
            ));
        }
        Ok(squared.sqrt())
    }

    #[allow(clippy::too_many_arguments)]
    fn timed_mixed_sdirk2<G: FluidTransport>(
        &self,
        trial: MaterialSdirk2Trial,
        start: f64,
        end: f64,
        previous: &[f64],
        force: &dyn Fn(f64, [f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesTimeBoundary<'_>],
        controls: FluxNewtonControls,
        resolve: &mut impl FnMut(
            f64,
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<MaterialFluxViscosityTangent<G>, MaterialError>,
        stage_solves: &mut usize,
        gate: &CancelGate,
    ) -> Result<MaterialNavierStokesSdirk2Step<G>, MaterialSdirk2AdaptiveError> {
        poll(gate)?;
        let [first_time, endpoint_time] = stage_times(start, end)?;
        let time = |stage| match stage {
            FluxSdirk2Stage::First => first_time,
            FluxSdirk2Stage::Endpoint => endpoint_time,
        };
        let result = with_time_boundaries(boundaries, first_time, |first| {
            with_time_boundaries(boundaries, endpoint_time, |endpoint| {
                self.sdirk2_navier_stokes_mixed_rheology_newton_step(
                    &|stage, point| force(time(stage), point),
                    first,
                    endpoint,
                    previous,
                    end - start,
                    controls,
                    |stage, point, reference| resolve(time(stage), point, reference),
                )
            })
        });
        *stage_solves += match &result {
            Ok(_) => 2,
            Err(error) => match error.stage {
                Some(FluxSdirk2Stage::First) => 1,
                Some(FluxSdirk2Stage::Endpoint) => 2,
                None => 0,
            },
        };
        result.map_err(|error| {
            MaterialSdirk2AdaptiveError::Stage(MaterialSdirk2StageFailure {
                trial,
                time_s: error.stage.map_or(start, time),
                error,
            })
        })
    }
}

impl MaterialSdirk2StageFailure {
    fn retryable(&self) -> bool {
        matches!(
            self.error.error,
            MaterialFluxError::Newton(
                FluxNewtonError::LinearSolve { .. }
                    | FluxNewtonError::LineSearch { .. }
                    | FluxNewtonError::NonConvergence { .. }
            )
        )
    }
}

impl MaterialSdirk2AdaptiveControls {
    fn admit(&self, start: f64, end: f64) -> Result<(), MaterialSdirk2AdaptiveError> {
        let pi = &self.controller;
        if !start.is_finite()
            || !end.is_finite()
            || end <= start
            || !(end - start).is_finite()
            || !self.initial_step_s.is_finite()
            || !self.minimum_step_s.is_finite()
            || !self.maximum_step_s.is_finite()
            || self.minimum_step_s <= 0.0
            || self.initial_step_s < self.minimum_step_s
            || self.initial_step_s > self.maximum_step_s
            || !self.velocity_absolute_tolerance_m_s.is_finite()
            || self.velocity_absolute_tolerance_m_s <= 0.0
            || !self.velocity_relative_tolerance.is_finite()
            || self.velocity_relative_tolerance < 0.0
            || self.maximum_attempts == 0
            || self.maximum_attempts > usize::MAX / 6
            || self.newton.maximum_iterations == 0
            || self.newton.maximum_backtracks > 52
            || !self.newton.residual_tolerance.is_finite()
            || self.newton.residual_tolerance <= 0.0
            || !self.previous_error_ratio.is_finite()
            || self.previous_error_ratio <= 0.0
            || self.previous_error_ratio > 1.0
            || !pi.k_p.is_finite()
            || pi.k_p < 0.0
            || !pi.k_i.is_finite()
            || pi.k_i < 0.0
            || !pi.safety.is_finite()
            || pi.safety <= 0.0
            || pi.safety > 1.0
            || !pi.max_growth.is_finite()
            || pi.max_growth < 1.0
            || !pi.max_shrink.is_finite()
            || pi.max_shrink <= 0.0
            || pi.max_shrink >= 1.0
        {
            return Err(MaterialSdirk2AdaptiveError::Invalid(
                "forward interval, step, tolerance, controller or work bounds",
            ));
        }
        Ok(())
    }

    fn factor(&self, error: f64, previous: f64) -> Result<f64, MaterialSdirk2AdaptiveError> {
        let pi = &self.controller;
        let factor = if error <= 1.0 {
            pi.safety
                * fs_math::det::pow(error.max(1e-300), -pi.k_p)
                * fs_math::det::pow(previous, pi.k_i)
        } else {
            pi.safety * fs_math::det::pow(error, -1.0 / 3.0)
        };
        if factor.is_nan() {
            return Err(MaterialSdirk2AdaptiveError::Invalid(
                "unrepresentable PI controller factor",
            ));
        }
        Ok(factor.clamp(
            pi.max_shrink,
            if error <= 1.0 { pi.max_growth } else { 1.0 },
        ))
    }
}

fn stage_times(start: f64, end: f64) -> Result<[f64; 2], MaterialSdirk2AdaptiveError> {
    let dt = end - start;
    let first = start + FluxSdirk2Stage::First.time_fraction() * dt;
    if !dt.is_finite() || dt <= 0.0 || first <= start || first >= end {
        return Err(MaterialSdirk2AdaptiveError::StepUnderflow {
            time_s: start,
            timestep_s: dt,
        });
    }
    Ok([first, end])
}

fn poll(gate: &CancelGate) -> Result<(), MaterialSdirk2AdaptiveError> {
    if gate.is_requested() {
        Err(MaterialSdirk2AdaptiveError::Cancelled)
    } else {
        Ok(())
    }
}

// Homogeneous closure vectors avoid a box per edge and reuse native admission;
// the source owns every boundary datum and no pressure/velocity law is invented.
fn with_time_boundaries<T>(
    boundaries: &[MaterialNavierStokesTimeBoundary<'_>],
    time: f64,
    solve: impl FnOnce(&[MaterialNavierStokesBoundary<'_>]) -> T,
) -> T {
    let velocities: Vec<_> = boundaries
        .iter()
        .map(|boundary| {
            move |point| match boundary {
                MaterialNavierStokesTimeBoundary::Velocity { value, .. }
                | MaterialNavierStokesTimeBoundary::NavierSlip { value, .. } => value(time, point),
                MaterialNavierStokesTimeBoundary::OpenReservoir { velocity, .. } => {
                    velocity(time, point)
                }
                MaterialNavierStokesTimeBoundary::Open { .. } => {
                    unreachable!("open stress has no prescribed velocity")
                }
            }
        })
        .collect();
    let tractions: Vec<_> = boundaries
        .iter()
        .map(|boundary| {
            move |point, normal| match boundary {
                MaterialNavierStokesTimeBoundary::Open { traction, .. }
                | MaterialNavierStokesTimeBoundary::OpenReservoir { traction, .. } => {
                    traction(time, point, normal)
                }
                _ => unreachable!("velocity/slip edge has no prescribed traction"),
            }
        })
        .collect();
    let sampled: Vec<_> = boundaries
        .iter()
        .zip(&velocities)
        .zip(&tractions)
        .map(|((boundary, velocity), traction)| match *boundary {
            MaterialNavierStokesTimeBoundary::Velocity { edge, .. } => {
                MaterialNavierStokesBoundary::Velocity {
                    edge,
                    value: velocity,
                }
            }
            MaterialNavierStokesTimeBoundary::NavierSlip {
                edge,
                tangential_resistance_pa_s_per_m,
                ..
            } => MaterialNavierStokesBoundary::NavierSlip {
                edge,
                value: velocity,
                tangential_resistance_pa_s_per_m: tangential_resistance_pa_s_per_m(time),
            },
            MaterialNavierStokesTimeBoundary::Open { edge, backflow, .. } => {
                MaterialNavierStokesBoundary::Open {
                    edge,
                    traction,
                    backflow,
                }
            }
            MaterialNavierStokesTimeBoundary::OpenReservoir { edge, .. } => {
                MaterialNavierStokesBoundary::OpenReservoir {
                    edge,
                    traction,
                    velocity,
                }
            }
        })
        .collect();
    solve(&sampled)
}
