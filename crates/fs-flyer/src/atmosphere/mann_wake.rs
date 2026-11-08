//! Full-field Mann turbulence in the native connected material wake dynamics.
//!
//! Every vertex and every rejected/adaptive stage observes one conditioned OU
//! history over the complete shedding interval. Stage frames are shared by all
//! vertices; no interpolated amplitudes or independent per-vertex noise enters
//! the trajectory. The original VPM owners retain circulation, topology, core
//! evolution, images, error estimates and transactional wake publication.
//!
//! Positions, velocities and gradients use the same right-handed z-up Mann/VPM
//! frame, not the flight integrator's forward/right/down world chart. The
//! uniform mean remains `SelfInducedShed::ambient_velocity_mps`, unless the
//! optional original log law supplies it. This prescribed atmospheric forcing
//! does not add a wake SGS closure or resolve viscous vortex-tube stretching.
//! OU paths are continuous but not differentiable: the native integrator's
//! formal order for smooth forcing is not claimed for these stochastic paths.

use std::cell::RefCell;

use fs_atmo::mann_field::affine::{MannAffineAccuracy, MannAffineError, MannAffineField};
use fs_atmo::mann_field::fixed_step::{
    FixedStepMannCheckpoint, FixedStepMannField, MAX_MANN_OU_BRIDGE_DEPTH, MannOuBridge,
    MannOuBridgeFrame,
};
use fs_atmo::mann_field::wall::MannPlaneWall;
use fs_atmo::{FlatSiteLogLaw, FlowSample};
use fs_vpm::farfield3d::FarField;
use fs_vpm::filament3d::{
    FilamentWake, Refusal, SelfInducedShed, WakeAdaptiveReport, WakeAdvectionReport, WakeBoundary,
    WakeTimeSpec, WakeTreeSpec,
};

/// Physical observation and shared conditional-path identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannWakeConfig {
    pub bridge_seed: u64,
    pub maximum_bridge_depth: usize,
    /// Optional slip plane in the common z-up frame. Native wake images must
    /// declare this same plane; an incompatible normal mean is refused.
    pub wall: Option<MannPlaneWall>,
    /// Optional original neutral log law, directed along positive x. Requires
    /// `wall` for the height origin and zero native uniform ambient velocity.
    pub mean: Option<FlatSiteLogLaw>,
}

/// Additional atmospheric work over one complete attempted wake interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MannWakeBudget {
    /// Includes bridge/endpoint preparation and every conditional path visit,
    /// including revisited stage times after rejected wake trials.
    pub maximum_normal_draws: usize,
    /// Counts all native incident-field queries. Exact queries cost O(modes);
    /// affine queries cost O(1), with O(modes) preparation per observed stage.
    pub maximum_field_samples: usize,
}

/// Explicit spatial observation for one attempted physical wake interval.
/// The anchor remains fixed across vertices, adaptive trials and stage order.
/// Repeat the same policy when replaying the native wake/wind checkpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MannWakeSampling {
    Exact,
    Affine {
        /// In the common Mann/VPM frame; must lie on the declared wall, if any.
        anchor_m: [f64; 3],
        /// Bounds the represented spatial field, not wake trajectory error.
        accuracy: MannAffineAccuracy,
    },
}

/// Original integrator and spatial-approximation choice, with no hidden default.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MannWakeIntegration {
    Midpoint {
        tree: Option<WakeTreeSpec>,
    },
    Adaptive {
        time: WakeTimeSpec,
        tree: Option<WakeTreeSpec>,
    },
}

/// Original native temporal/spatial evidence from the successful interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MannWakeAdvection {
    /// The native direct midpoint entry point has no field-work report.
    DirectMidpoint,
    TreeMidpoint(WakeAdvectionReport),
    Adaptive(WakeAdaptiveReport),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannWakeStepReport {
    pub advection: MannWakeAdvection,
    pub normal_draws: usize,
    pub field_samples: usize,
    /// Conditional frame updates; consecutive vertices at one stage share it.
    pub stage_frames: usize,
    pub maximum_bridge_depth_used: usize,
    /// Once per visited OU stage in affine mode, including revisited stages.
    /// All vertices at the same stage reuse this fixed-size observation.
    pub affine_frames: usize,
    /// Largest admitted analytical spatial truncation estimates. These use
    /// ordinary f64 arithmetic, not interval-certified error bounds.
    pub maximum_affine_error: MannAffineError,
    pub start_time_s: f64,
    pub end_time_s: f64,
}

/// Store alongside the accepted native near/far wake, including its tick.
#[derive(Clone, Debug, PartialEq)]
pub struct MannWakeCheckpoint {
    pub version: u32,
    pub field: FixedStepMannCheckpoint,
    pub config: MannWakeConfig,
    pub next_wake_tick: u64,
}

/// Atmospheric owner for successive native shedding intervals. Retained
/// endpoints reuse their storage. One bounded bridge and frame are prepared
/// per interval; successful per-vertex/stage observations allocate nothing.
#[derive(Clone, Debug)]
pub struct MannWakeWind {
    accepted: FixedStepMannField,
    endpoint: FixedStepMannField,
    config: MannWakeConfig,
    next_wake_tick: u64,
    maximum_workspace_bytes: usize,
}

fn error(code: &'static str, message: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![
            "use the same admitted wind/wake clock, physical frame and explicit budgets".into(),
        ],
    }
}

fn atmosphere(error: fs_atmo::Refusal) -> Refusal {
    Refusal {
        code: error.code,
        message: error.message,
        ranked_repairs: error.ranked_repairs,
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(error(
            "mann-wake-cancelled",
            "atmospheric wake interval cancelled",
        ))
    } else {
        Ok(())
    }
}

impl MannWakeWind {
    /// Bind this snapshot to the supplied native wake's current tick. The
    /// source time origin is retained; it need not equal tick / shedding rate.
    /// The source dt must equal the native interval exactly. The byte budget
    /// covers two retained endpoints plus one transient bridge and frame,
    /// conservatively counting shared temporal arrays once per owner. Original
    /// spectral geometry, fixed metadata and native wake storage are separate.
    pub fn prepare(
        field: FixedStepMannField,
        wake: &FilamentWake,
        config: MannWakeConfig,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, Refusal> {
        if field.dt_s() != 1.0 / wake.cert.shed_hz {
            return Err(error(
                "mann-wake-clock",
                "wind dt differs from the native shedding interval",
            ));
        }
        Self::from_parts(field, config, wake.ticks, maximum_workspace_bytes)
    }

    fn from_parts(
        field: FixedStepMannField,
        config: MannWakeConfig,
        next_wake_tick: u64,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, Refusal> {
        if !(1..=MAX_MANN_OU_BRIDGE_DEPTH).contains(&config.maximum_bridge_depth)
            || (config.mean.is_some() && config.wall.is_none())
        {
            return Err(error(
                "mann-wake-invalid",
                "invalid bridge depth or log-law height origin",
            ));
        }
        if let Some(mean) = config.mean {
            mean.admit().map_err(atmosphere)?;
        }
        let peak = field
            .workspace_bytes()
            .checked_mul(4)
            .ok_or_else(|| error("mann-wake-budget", "wind workspace size overflows"))?;
        if peak > maximum_workspace_bytes {
            return Err(error(
                "mann-wake-budget",
                "wind endpoint/bridge workspace exceeds its byte budget",
            ));
        }
        let result = Self {
            endpoint: field.clone(),
            accepted: field,
            config,
            next_wake_tick,
            maximum_workspace_bytes,
        };
        if result
            .workspace_bytes()
            .checked_mul(2)
            .is_none_or(|v| v > maximum_workspace_bytes)
        {
            return Err(error(
                "mann-wake-budget",
                "allocated wind capacities exceed the byte budget",
            ));
        }
        Ok(result)
    }

    /// Retained numerical capacities; a step admits bridge/frame capacities too.
    pub fn workspace_bytes(&self) -> usize {
        self.accepted.workspace_bytes() + self.endpoint.workspace_bytes()
    }
    pub fn time_s(&self) -> f64 {
        self.accepted.time_s()
    }
    pub fn dt_s(&self) -> f64 {
        self.accepted.dt_s()
    }

    /// Accepted atmospheric velocity and analytic gradient in the common frame.
    /// The separate native uniform ambient velocity is not included here.
    pub fn sample(&self, point_m: [f64; 3]) -> Result<FlowSample, Refusal> {
        let sample = match self.config.wall {
            Some(wall) => self.accepted.sample_wall_blocked(point_m, wall),
            None => self.accepted.sample(point_m),
        }
        .map_err(atmosphere)?;
        add_mean(sample, point_m, self.config)
    }

    pub fn checkpoint(&self) -> MannWakeCheckpoint {
        MannWakeCheckpoint {
            version: 1,
            field: self.accepted.checkpoint(),
            config: self.config,
            next_wake_tick: self.next_wake_tick,
        }
    }

    pub fn resume(
        checkpoint: MannWakeCheckpoint,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, Refusal> {
        if checkpoint.version != 1 {
            return Err(error(
                "mann-wake-checkpoint",
                "unsupported atmospheric wake checkpoint version",
            ));
        }
        let field = FixedStepMannField::resume(checkpoint.field, maximum_workspace_bytes / 4)
            .map_err(atmosphere)?;
        Self::from_parts(
            field,
            checkpoint.config,
            checkpoint.next_wake_tick,
            maximum_workspace_bytes,
        )
    }

    /// Advance the actual near wake and atmospheric endpoint together. Native
    /// rejected trials observe the same bridge; failure changes neither owner.
    #[allow(clippy::too_many_arguments)]
    pub fn shed_near(
        &mut self,
        near: &mut FilamentWake,
        gamma: &[f64],
        config: SelfInducedShed,
        boundary: WakeBoundary<'_>,
        integration: MannWakeIntegration,
        budget: MannWakeBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannWakeStepReport, Refusal> {
        self.shed_near_with_sampling(
            near,
            gamma,
            config,
            boundary,
            integration,
            MannWakeSampling::Exact,
            budget,
            cancelled,
        )
    }

    /// Native near-wake advection with an explicitly budgeted affine option.
    /// Every new OU stage rebuilds the affine value, gradient and error bounds;
    /// a spatial-budget refusal rolls back the complete native transaction.
    #[allow(clippy::too_many_arguments)]
    pub fn shed_near_with_sampling(
        &mut self,
        near: &mut FilamentWake,
        gamma: &[f64],
        config: SelfInducedShed,
        boundary: WakeBoundary<'_>,
        integration: MannWakeIntegration,
        sampling: MannWakeSampling,
        budget: MannWakeBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannWakeStepReport, Refusal> {
        self.advance_interval(
            near.ticks,
            near.cert.shed_hz,
            config,
            boundary,
            sampling,
            budget,
            cancelled,
            |variation, cancelled| match integration {
                MannWakeIntegration::Midpoint { tree: None } => near
                    .shed_self_induced_with_environment(
                        gamma, config, boundary, variation, cancelled,
                    )
                    .map(|()| MannWakeAdvection::DirectMidpoint),
                MannWakeIntegration::Midpoint { tree: Some(tree) } => near
                    .shed_self_induced_with_tree(
                        gamma, config, boundary, tree, variation, cancelled,
                    )
                    .map(MannWakeAdvection::TreeMidpoint),
                MannWakeIntegration::Adaptive { time, tree } => near
                    .shed_self_induced_adaptive(
                        gamma, config, time, boundary, tree, variation, cancelled,
                    )
                    .map(MannWakeAdvection::Adaptive),
            },
        )
    }

    /// The same single history drives retained near/far sources, including
    /// their shared vertices. Publication includes both physical partitions.
    #[allow(clippy::too_many_arguments)]
    pub fn shed_hybrid(
        &mut self,
        near: &mut FilamentWake,
        far: &mut FarField,
        gamma: &[f64],
        config: SelfInducedShed,
        boundary: WakeBoundary<'_>,
        integration: MannWakeIntegration,
        budget: MannWakeBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannWakeStepReport, Refusal> {
        self.shed_hybrid_with_sampling(
            near,
            far,
            gamma,
            config,
            boundary,
            integration,
            MannWakeSampling::Exact,
            budget,
            cancelled,
        )
    }

    /// The same authored affine field drives both native near/far partitions.
    /// No partition-specific fit or independent turbulent sample is generated.
    #[allow(clippy::too_many_arguments)]
    pub fn shed_hybrid_with_sampling(
        &mut self,
        near: &mut FilamentWake,
        far: &mut FarField,
        gamma: &[f64],
        config: SelfInducedShed,
        boundary: WakeBoundary<'_>,
        integration: MannWakeIntegration,
        sampling: MannWakeSampling,
        budget: MannWakeBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannWakeStepReport, Refusal> {
        self.advance_interval(
            near.ticks,
            near.cert.shed_hz,
            config,
            boundary,
            sampling,
            budget,
            cancelled,
            |variation, cancelled| match integration {
                MannWakeIntegration::Midpoint { tree: None } => far
                    .shed_self_induced_with_environment(
                        near, gamma, config, boundary, variation, cancelled,
                    )
                    .map(|()| MannWakeAdvection::DirectMidpoint),
                MannWakeIntegration::Midpoint { tree: Some(tree) } => far
                    .shed_self_induced_with_tree(
                        near, gamma, config, boundary, tree, variation, cancelled,
                    )
                    .map(MannWakeAdvection::TreeMidpoint),
                MannWakeIntegration::Adaptive { time, tree } => far
                    .shed_self_induced_adaptive(
                        near, gamma, config, time, boundary, tree, variation, cancelled,
                    )
                    .map(MannWakeAdvection::Adaptive),
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_interval(
        &mut self,
        wake_tick: u64,
        shed_hz: f64,
        config: SelfInducedShed,
        boundary: WakeBoundary<'_>,
        sampling: MannWakeSampling,
        budget: MannWakeBudget,
        mut cancelled: impl FnMut() -> bool,
        advance: impl FnOnce(
            &mut dyn FnMut([f64; 3], f64) -> Result<[f64; 3], Refusal>,
            &mut dyn FnMut() -> bool,
        ) -> Result<MannWakeAdvection, Refusal>,
    ) -> Result<MannWakeStepReport, Refusal> {
        poll(&mut cancelled)?;
        if let MannWakeSampling::Affine { anchor_m, accuracy } = sampling {
            accuracy.admit().map_err(atmosphere)?;
            if !anchor_m.iter().all(|x| x.is_finite())
                || self
                    .config
                    .wall
                    .is_some_and(|wall| anchor_m[2] != wall.height_m())
            {
                return Err(error(
                    "mann-affine-anchor",
                    "affine wake anchor must be finite and lie on the declared wall",
                ));
            }
        }
        if wake_tick != self.next_wake_tick || self.dt_s() != 1.0 / shed_hz {
            return Err(error(
                "mann-wake-clock",
                "wake and wind are not at the bound shedding interval",
            ));
        }
        let next_tick = wake_tick
            .checked_add(1)
            .ok_or_else(|| error("mann-wake-clock", "wake tick overflows"))?;
        let next_step = self
            .accepted
            .step()
            .checked_add(1)
            .ok_or_else(|| error("mann-wake-clock", "wind step overflows"))?;
        let matching_wall = match (self.config.wall, boundary.image_plane) {
            (None, None) => true,
            (Some(wall), Some(plane)) => {
                wall.height_m() == plane.z0_m && config.ambient_velocity_mps[2] == 0.0
            }
            _ => false,
        };
        if !matching_wall || (self.config.mean.is_some() && config.ambient_velocity_mps != [0.0; 3])
        {
            return Err(error(
                "mann-wake-boundary",
                "incident wind, mean and native wake images are incompatible",
            ));
        }
        // The original bridge and original accepted-endpoint transition each
        // evaluate the same 4n counter-addressed innovations. Count both costs;
        // neither introduces an independent endpoint realization.
        let endpoint_draws = 8 * self.accepted.model().mode_count();
        if endpoint_draws > budget.maximum_normal_draws {
            return Err(error(
                "mann-wake-work",
                "wind endpoints exceed the interval normal-draw budget",
            ));
        }
        self.endpoint
            .copy_temporal_state_from(&self.accepted)
            .map_err(atmosphere)?;
        self.endpoint
            .advance_to(next_step, endpoint_draws / 2, &mut cancelled)
            .map_err(atmosphere)?;
        let remaining_bytes = self.maximum_workspace_bytes - self.workspace_bytes();
        let bridge = MannOuBridge::prepare(
            &self.accepted,
            self.config.bridge_seed,
            self.config.maximum_bridge_depth,
            remaining_bytes / 2,
            endpoint_draws / 2,
            &mut cancelled,
        )
        .map_err(atmosphere)?;
        let frame = bridge
            .frame(remaining_bytes - bridge.workspace_bytes(), &mut cancelled)
            .map_err(atmosphere)?;
        let mut stage = StageWind {
            bridge,
            frame,
            config: self.config,
            sampling,
            affine: None,
            budget,
            normal_draws: endpoint_draws,
            field_samples: 0,
            stage_frames: 0,
            maximum_depth: 0,
            affine_frames: 0,
            maximum_affine_error: MannAffineError::default(),
        };
        let cancelled = RefCell::new(cancelled);
        let advection = advance(
            &mut |point, elapsed| {
                stage
                    .sample(point, elapsed, || (cancelled.borrow_mut())())
                    .map(|flow| flow.u)
            },
            &mut || (cancelled.borrow_mut())(),
        )?;
        let report = MannWakeStepReport {
            advection,
            normal_draws: stage.normal_draws,
            field_samples: stage.field_samples,
            stage_frames: stage.stage_frames,
            maximum_bridge_depth_used: stage.maximum_depth,
            affine_frames: stage.affine_frames,
            maximum_affine_error: stage.maximum_affine_error,
            start_time_s: stage.bridge.start_time_s(),
            end_time_s: stage.bridge.end_time_s(),
        };
        // Native VPM has performed its final cancellation check and published
        // the wake. Only infallible endpoint publication may follow it.
        core::mem::swap(&mut self.accepted, &mut self.endpoint);
        self.next_wake_tick = next_tick;
        Ok(report)
    }
}

struct StageWind {
    bridge: MannOuBridge,
    frame: MannOuBridgeFrame,
    config: MannWakeConfig,
    sampling: MannWakeSampling,
    affine: Option<MannAffineField>,
    budget: MannWakeBudget,
    normal_draws: usize,
    field_samples: usize,
    stage_frames: usize,
    maximum_depth: usize,
    affine_frames: usize,
    maximum_affine_error: MannAffineError,
}

impl StageWind {
    fn sample(
        &mut self,
        point: [f64; 3],
        elapsed: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FlowSample, Refusal> {
        poll(&mut cancelled)?;
        if !elapsed.is_finite() || elapsed < 0.0 || elapsed > self.bridge.dt_s() {
            return Err(error(
                "mann-wake-clock",
                "native stage is outside its complete bridge interval",
            ));
        }
        if self.field_samples >= self.budget.maximum_field_samples {
            return Err(error(
                "mann-wake-work",
                "native incident-field queries exceed their interval budget",
            ));
        }
        self.field_samples += 1;
        // Match both original physical endpoints exactly. The bridge already
        // maps their represented gap to the nominal OU interval. This maps the
        // native elapsed stage clock the same way, without interpolating noise.
        let time = if elapsed == 0.0 {
            self.bridge.start_time_s()
        } else if elapsed == self.bridge.dt_s() {
            self.bridge.end_time_s()
        } else {
            let time = self.bridge.start_time_s()
                + (self.bridge.end_time_s() - self.bridge.start_time_s())
                    * (elapsed / self.bridge.dt_s());
            if time <= self.bridge.start_time_s() || time >= self.bridge.end_time_s() {
                return Err(error(
                    "mann-wake-clock",
                    "interior wake stage has no distinct physical time",
                ));
            }
            time
        };
        if time != self.frame.time_s() {
            let report = self
                .bridge
                .sample_frame(
                    time,
                    &mut self.frame,
                    self.budget.maximum_normal_draws - self.normal_draws,
                    &mut cancelled,
                )
                .map_err(atmosphere)?;
            self.normal_draws += report.normal_draws;
            self.stage_frames += 1;
            self.maximum_depth = self.maximum_depth.max(report.tree_depth);
            // A returning adaptive stage rebuilds from its original bridge
            // realization; an affine field from another time is never reused.
            self.affine = None;
        }
        let sample = match self.sampling {
            MannWakeSampling::Exact => match self.config.wall {
                Some(wall) => self.frame.sample_wall_blocked(point, wall),
                None => self.frame.sample(point),
            }
            .map_err(atmosphere)?,
            MannWakeSampling::Affine { anchor_m, accuracy } => {
                if self.affine.is_none() {
                    self.affine = Some(
                        self.frame
                            .affine_observation(anchor_m, self.config.wall, &mut cancelled)
                            .map_err(atmosphere)?,
                    );
                    self.affine_frames += 1;
                }
                let sample = self
                    .affine
                    .as_ref()
                    .expect("prepared current stage")
                    .sample(point, accuracy)
                    .map_err(atmosphere)?;
                self.maximum_affine_error.velocity_mps = self
                    .maximum_affine_error
                    .velocity_mps
                    .max(sample.error.velocity_mps);
                self.maximum_affine_error.gradient_per_s = self
                    .maximum_affine_error
                    .gradient_per_s
                    .max(sample.error.gradient_per_s);
                sample.flow
            }
        };
        add_mean(sample, point, self.config)
    }
}

fn add_mean(
    mut sample: FlowSample,
    point: [f64; 3],
    config: MannWakeConfig,
) -> Result<FlowSample, Refusal> {
    if let Some(mean) = config.mean {
        let height = point[2]
            - config
                .wall
                .expect("admitted log-law height origin")
                .height_m();
        sample.u[0] += mean.speed(height);
        sample.grad[0][2] += mean.dspeed_dh(height);
    }
    if !sample
        .u
        .iter()
        .chain(sample.grad.iter().flatten())
        .all(|x| x.is_finite())
    {
        return Err(error(
            "mann-wake-numerical",
            "incident wind or analytic gradient is not finite",
        ));
    }
    Ok(sample)
}

#[cfg(test)]
mod tests;
