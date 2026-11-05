//! Fixed-clock evolution of the complete original nonlinear plate and material history.
use super::{AcousticRealizeError, VkBody, render_step_error};
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplit, LinearSplitBudget, PreparedExactZoh,
    PreparedNonlinearTones, admit_midpoint_band,
};
use fs_phs::{
    PortHamiltonian, PreparedStepError, PreparedStepRecord, ResistanceWitness, StepWorkspace,
};

/// Accepted compact acoustic observation and the original mechanical/material work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PreparedVkStep {
    /// Compact baffled-monopole pressure from the interval-mean acceleration [Pa].
    pub pressure_pa: f64,
    /// Accepted endpoint volume velocity [m^3/s].
    pub volume_velocity_m3_s: f64,
    /// Complete storage, material/acoustic losses and physical input work [J].
    /// Acoustic loss must not also be routed as solid material heat.
    pub work: PreparedStepRecord,
}

struct NonlinearEvolution {
    system: PortHamiltonian,
    workspace: StepWorkspace,
}

enum Evolution {
    Exact(Box<PreparedExactZoh<'static>>),
    Nonlinear(Box<NonlinearEvolution>),
}
impl Evolution {
    fn system(&self) -> &PortHamiltonian {
        match self {
            Self::Exact(flow) => flow.system(),
            Self::Nonlinear(flow) => &flow.system,
        }
    }
}

/// The same plate, damping and material histories on a fixed physical clock.
///
/// Quadratic bodies cache the complete exponential and independent heat/work
/// integrals. Nonlinear bodies retain the canonical Gonzalez equation and its
/// reusable Newton workspace. Both keep all physical states and cross terms.
/// Successful warm steps allocate no heap memory for these plate storage owners.
/// This is a numerical reuse guarantee, not a hard real-time timing certificate.
pub struct PreparedVkBody {
    evolution: Evolution,
    x: Vec<f64>,
    candidate: Vec<f64>,
    output: Vec<f64>,
    areas: Vec<f64>,
    dt_s: f64,
}

impl VkBody {
    /// Prepare this complete body for repeated steps of `dt_s`, preserving its
    /// accepted displacement, momentum and constitutive memory exactly.
    ///
    /// Preparation consumes the original body and allocates caches/scratch.
    /// An already moving body is allowed: this changes numerical preparation,
    /// not its material law or physical state. The original storage owner
    /// chooses exact quadratic propagation or nonlinear discrete gradients.
    ///
    /// # Errors
    /// Invalid clock/projections, original render capacity, unresolved nonlinear
    /// pole-rate band, or refusal by the existing preparation/solver owners.
    pub fn into_prepared(self, dt_s: f64) -> Result<PreparedVkBody, AcousticRealizeError> {
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || self.sys.port_dim() != 1
            || self.areas.is_empty()
            || self
                .areas
                .len()
                .checked_mul(2)
                .is_none_or(|n| n > self.x.len())
            || self.areas.iter().chain(&self.x).any(|x| !x.is_finite())
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "prepared nonlinear plate needs a positive fixed clock and finite original state/projections",
            });
        }
        let n = self.sys.state_dim();
        let evolution = if self.linear.is_exact_quadratic() {
            Evolution::Exact(Box::new(
                PreparedExactZoh::new_owned(
                    self.sys,
                    dt_s,
                    LinearSplitBudget {
                        max_modes: n,
                        max_omega_dt: core::f64::consts::PI,
                    },
                )
                .map_err(render_step_error)?,
            ))
        } else {
            admit_midpoint_band(self.linear.pole_rate_bound_rad_s(), dt_s)
                .map_err(|e| render_step_error(fs_phs::PhsError::TimeResolution(e)))?;
            let workspace = StepWorkspace::new(&self.sys).map_err(render_step_error)?;
            Evolution::Nonlinear(Box::new(NonlinearEvolution {
                system: self.sys,
                workspace,
            }))
        };
        Ok(PreparedVkBody {
            evolution,
            x: self.x,
            candidate: vec![0.0; n],
            output: vec![0.0; 1],
            areas: self.areas,
            dt_s,
        })
    }
}

impl PreparedVkBody {
    /// The fixed physical duration of every accepted step [s].
    #[must_use]
    pub fn time_step_s(&self) -> f64 {
        self.dt_s
    }

    /// Original `[q0,p0,q1,p1,...]` followed by constitutive and acoustic/cavity
    /// histories in their attachment order, retaining each load's coordinates.
    #[must_use]
    pub fn state(&self) -> &[f64] {
        &self.x
    }

    /// Complete mechanical, constitutive and attached acoustic storage [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.evolution.system().hamiltonian(&self.x)
    }

    /// Endpoint volume velocity in the original surface projection [m^3/s].
    #[must_use]
    pub fn volume_velocity(&self) -> f64 {
        self.areas
            .iter()
            .enumerate()
            .map(|(k, a)| a * self.x[2 * k + 1])
            .sum()
    }

    /// Advance the unchanged body under one held physical force [N].
    /// Density and distance affect only the compact acoustic observation.
    /// No physical state is published on solver or observation refusal.
    pub fn drive_and_radiate(
        &mut self,
        force_n: f64,
        rho: f64,
        listener_m: f64,
    ) -> Result<PreparedVkStep, AcousticRealizeError> {
        self.drive_and_radiate_controlled(force_n, rho, listener_m, || false)
    }

    /// Same physical step with cooperative cancellation. Nonlinear Newton polls
    /// its original boundaries; a bounded cached linear contraction is atomic.
    /// Cancellation before publication leaves every physical history unchanged.
    pub fn drive_and_radiate_controlled(
        &mut self,
        force_n: f64,
        rho: f64,
        listener_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedVkStep, AcousticRealizeError> {
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        if !force_n.is_finite()
            || !rho.is_finite()
            || rho <= 0.0
            || !listener_m.is_finite()
            || listener_m <= 0.0
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "prepared plate step requires finite force and positive finite observer density/distance",
            });
        }
        let work = match &mut self.evolution {
            Evolution::Exact(flow) => flow
                .step_into(&self.x, &[force_n], &mut self.candidate, &mut self.output)
                .map_err(render_step_error)?,
            Evolution::Nonlinear(flow) => flow
                .workspace
                .step_into_controlled(
                    &flow.system,
                    &self.x,
                    &[force_n],
                    self.dt_s,
                    &mut self.candidate,
                    &mut self.output,
                    &mut cancelled,
                )
                .map_err(|e| match e {
                    PreparedStepError::Cancelled => AcousticRealizeError::Cancelled,
                    PreparedStepError::Solver(error) => render_step_error(error),
                })?,
        };
        let frame = compact_step(
            &self.areas,
            &self.x,
            &self.candidate,
            self.dt_s,
            rho,
            listener_m,
            work,
        )?;
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        core::mem::swap(&mut self.x, &mut self.candidate);
        Ok(frame)
    }
}

/// Fixed mechanical/observation clock and complete continuous-source resources.
#[derive(Clone, Copy, Debug)]
pub struct VkToneConfig {
    /// Duration of every accepted physical interval [s].
    pub time_step_s: f64,
    /// Fraction of Nyquist allowed for every retained linear mode and source
    /// carrier. This does not bound nonlinear or envelope-generated harmonics;
    /// those still require physical-clock refinement and output antialiasing.
    pub nyquist_guard_fraction: f64,
    /// Complete pHS cache and adapter numerical payload/work admission.
    /// The caller's coefficient arrays and allocation metadata are excluded.
    pub dense_budget: DenseExponentialBudget,
}

/// Original nonlinear plate motion under continuously evolving port efforts.
///
/// Retains every displacement, momentum and material/acoustic history. The pHS
/// owner advances the original quadratic block and source quadratures exactly,
/// with its existing nonlinear discrete-gradient endpoint solve. It does not
/// replace the nonlinear storage with independent linear oscillators.
///
/// Pressure is the same compact baffled-monopole interval-mean acceleration
/// observation as [`PreparedVkBody`]. It is source-time pressure at the declared
/// distance, with no propagation delay or new acoustic reaction inferred here.
/// Physical state, channel heat and the returned observation publish together.
/// Warm stepping and cached phase rotation allocate no adapter scratch.
pub struct PreparedVkTones {
    flow: PreparedNonlinearTones,
    x: Vec<f64>,
    candidate: Vec<f64>,
    output: Vec<f64>,
    areas: Vec<f64>,
    channel_heat: Vec<f64>,
    candidate_heat: Vec<f64>,
}

impl VkBody {
    /// Prepare constant-amplitude physical tones on the original force port.
    /// `channels` supplies optional original resistance observations, in caller
    /// order. The source owner supplies their physical interpretation; channel
    /// overlap is not turned into a partition of total material/acoustic heat.
    pub fn into_prepared_tones(
        self,
        frequencies: &[f64],
        config: VkToneConfig,
        channels: &[ResistanceWitness],
    ) -> Result<PreparedVkTones, AcousticRealizeError> {
        self.into_prepared_tones_controlled(frequencies, config, channels, || false)
    }

    /// Prepare constant-amplitude tones with cancellation during cold work.
    pub fn into_prepared_tones_controlled(
        self,
        frequencies: &[f64],
        config: VkToneConfig,
        channels: &[ResistanceWitness],
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedVkTones, AcousticRealizeError> {
        PreparedVkTones::prepare(
            self.sys,
            &self.linear,
            self.x,
            self.areas,
            frequencies,
            config,
            channels,
            false,
            cancelled,
        )
    }

    /// Prepare affine cosine/sine coefficient envelopes over each interval.
    /// Both endpoint pairs use the interval-start carrier basis. The caller
    /// rotates the accepted end pair to continue its source phase.
    pub fn into_prepared_affine_tones(
        self,
        frequencies: &[f64],
        config: VkToneConfig,
        channels: &[ResistanceWitness],
    ) -> Result<PreparedVkTones, AcousticRealizeError> {
        self.into_prepared_affine_tones_controlled(frequencies, config, channels, || false)
    }

    /// Prepare affine tone envelopes with cancellation during cold work.
    pub fn into_prepared_affine_tones_controlled(
        self,
        frequencies: &[f64],
        config: VkToneConfig,
        channels: &[ResistanceWitness],
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedVkTones, AcousticRealizeError> {
        PreparedVkTones::prepare(
            self.sys,
            &self.linear,
            self.x,
            self.areas,
            frequencies,
            config,
            channels,
            true,
            cancelled,
        )
    }
}

impl PreparedVkTones {
    // PlateBankImage uses this same owner for its complete original ports and
    // histories. The source supplies its matching split and signed projections;
    // the pHS preparation independently admits the complete physical system.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare(
        system: PortHamiltonian,
        linear: &LinearSplit,
        initial: Vec<f64>,
        areas: Vec<f64>,
        frequencies: &[f64],
        config: VkToneConfig,
        channels: &[ResistanceWitness],
        affine: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AcousticRealizeError> {
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let (n, m) = (system.state_dim(), system.port_dim());
        let dt = config.time_step_s;
        if !dt.is_finite()
            || dt <= 0.0
            || m == 0
            || initial.len() != n
            || areas.is_empty()
            || areas.len().checked_mul(2).is_none_or(|d| d > n)
            || areas.iter().chain(&initial).any(|x| !x.is_finite())
            || !config.nyquist_guard_fraction.is_finite()
            || !(0.0 < config.nyquist_guard_fraction && config.nyquist_guard_fraction < 1.0)
        {
            return Err(tone_invalid(
                "plate tones require complete finite state/projections and an explicit positive clock/Nyquist guard",
            ));
        }
        let maximum_phase = core::f64::consts::PI * config.nyquist_guard_fraction;
        if frequencies.is_empty()
            || frequencies
                .iter()
                .chain(linear.angular_frequencies_rad_s())
                .any(|w| !w.is_finite() || *w < 0.0 || *w * dt > maximum_phase)
        {
            return Err(tone_invalid(
                "all original plate modes and source carriers must fit the declared observation band",
            ));
        }
        let adapter_bytes = n
            .checked_mul(2)
            .and_then(|v| v.checked_add(areas.len()))
            .and_then(|v| v.checked_add(m))
            .and_then(|v| v.checked_add(channels.len().checked_mul(2)?))
            .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()))
            .ok_or_else(|| tone_invalid("plate tone state/channel payload overflows"))?;
        let mut dense = config.dense_budget;
        dense.max_dense_bytes = dense
            .max_dense_bytes
            .checked_sub(adapter_bytes)
            .ok_or_else(|| {
                tone_invalid("plate tone state/channel payload exceeds its byte budget")
            })?;
        let split = LinearSplitBudget {
            max_modes: n,
            max_omega_dt: core::f64::consts::PI,
        };
        let flow = if affine {
            PreparedNonlinearTones::new_owned_affine_with_resistance_channels_controlled(
                system,
                frequencies,
                dt,
                split,
                dense,
                channels,
                &mut cancelled,
            )
        } else {
            PreparedNonlinearTones::new_owned_with_resistance_channels_controlled(
                system,
                frequencies,
                dt,
                split,
                dense,
                channels,
                &mut cancelled,
            )
        }
        .map_err(tone_step_error)?;
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        Ok(Self {
            flow,
            x: initial,
            candidate: vec![0.0; n],
            output: vec![0.0; m],
            areas,
            channel_heat: vec![0.0; channels.len()],
            candidate_heat: vec![0.0; channels.len()],
        })
    }

    pub fn time_step_s(&self) -> f64 {
        self.flow.duration()
    }
    pub fn state(&self) -> &[f64] {
        &self.x
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.flow.system().hamiltonian(&self.x)
    }
    pub fn volume_velocity(&self) -> f64 {
        self.areas
            .iter()
            .enumerate()
            .map(|(k, area)| area * self.x[2 * k + 1])
            .sum()
    }
    pub fn port_count(&self) -> usize {
        self.flow.port_dim()
    }
    pub fn tone_count(&self) -> usize {
        self.flow.tone_count()
    }
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        self.flow.angular_frequencies_rad_s()
    }
    pub fn is_affine(&self) -> bool {
        self.flow.is_affine()
    }

    /// Last accepted interval's original resistance observations [J], in the
    /// constructor's channel order. Refused trials leave this slice unchanged.
    pub fn channel_heat_j(&self) -> &[f64] {
        &self.channel_heat
    }

    /// Rotate a caller-owned coefficient pair only after accepting its interval.
    pub fn advance_force_phase(
        &self,
        tone: usize,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), AcousticRealizeError> {
        self.flow
            .advance_force_phase(tone, cosine, sine)
            .map_err(render_step_error)
    }

    /// Continuous efforts in original port-major/tone-major order. The body's
    /// drive port uses force [N]; a complete bank's face port uses pressure [Pa].
    pub fn drive_and_radiate(
        &mut self,
        cosine_n: &[f64],
        sine_n: &[f64],
        rho: f64,
        listener_m: f64,
    ) -> Result<PreparedVkStep, AcousticRealizeError> {
        self.drive_and_radiate_controlled(cosine_n, sine_n, rho, listener_m, || false)
    }

    /// Cancellable constant-amplitude step; no state, heat or pressure publishes
    /// until the original nonlinear solve and finite acoustic observation pass.
    pub fn drive_and_radiate_controlled(
        &mut self,
        cosine_n: &[f64],
        sine_n: &[f64],
        rho: f64,
        listener_m: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedVkStep, AcousticRealizeError> {
        self.drive(cosine_n, sine_n, None, rho, listener_m, cancelled)
    }

    /// Affine effort endpoints in the original port units (force [N] or face
    /// pressure [Pa]), all in the interval-start carrier basis.
    #[allow(clippy::too_many_arguments)]
    pub fn drive_affine_and_radiate(
        &mut self,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        rho: f64,
        listener_m: f64,
    ) -> Result<PreparedVkStep, AcousticRealizeError> {
        self.drive_affine_and_radiate_controlled(
            cosine_start,
            sine_start,
            cosine_end,
            sine_end,
            rho,
            listener_m,
            || false,
        )
    }

    /// Cancellable affine step with the same accepted state/heat boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn drive_affine_and_radiate_controlled(
        &mut self,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        rho: f64,
        listener_m: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedVkStep, AcousticRealizeError> {
        self.drive(
            cosine_start,
            sine_start,
            Some((cosine_end, sine_end)),
            rho,
            listener_m,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn drive(
        &mut self,
        cosine: &[f64],
        sine: &[f64],
        end: Option<(&[f64], &[f64])>,
        rho: f64,
        listener_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedVkStep, AcousticRealizeError> {
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        if !rho.is_finite() || rho <= 0.0 || !listener_m.is_finite() || listener_m <= 0.0 {
            return Err(tone_invalid(
                "plate tone observation needs positive finite density and distance",
            ));
        }
        let work = if let Some((c, s)) = end {
            self.flow.step_affine_observed_into_controlled(
                &self.x,
                cosine,
                sine,
                c,
                s,
                &mut self.candidate,
                &mut self.output,
                &mut self.candidate_heat,
                &mut cancelled,
            )
        } else {
            self.flow.step_observed_into_controlled(
                &self.x,
                cosine,
                sine,
                &mut self.candidate,
                &mut self.output,
                &mut self.candidate_heat,
                &mut cancelled,
            )
        }
        .map_err(tone_step_error)?;
        let frame = compact_step(
            &self.areas,
            &self.x,
            &self.candidate,
            self.time_step_s(),
            rho,
            listener_m,
            work,
        )?;
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        core::mem::swap(&mut self.x, &mut self.candidate);
        core::mem::swap(&mut self.channel_heat, &mut self.candidate_heat);
        Ok(frame)
    }
}

#[allow(clippy::too_many_arguments)]
fn compact_step(
    areas: &[f64],
    initial: &[f64],
    candidate: &[f64],
    dt: f64,
    rho: f64,
    listener_m: f64,
    work: PreparedStepRecord,
) -> Result<PreparedVkStep, AcousticRealizeError> {
    let mut pressure_pa = 0.0;
    let mut volume_velocity_m3_s = 0.0;
    for (k, area) in areas.iter().enumerate() {
        let velocity = candidate[2 * k + 1];
        let acceleration = (velocity - initial[2 * k + 1]) / dt;
        pressure_pa += rho * area * acceleration / (core::f64::consts::TAU * listener_m);
        volume_velocity_m3_s += area * velocity;
    }
    if !pressure_pa.is_finite() || !volume_velocity_m3_s.is_finite() {
        return Err(tone_invalid(
            "prepared plate acoustic observation is not finite",
        ));
    }
    Ok(PreparedVkStep {
        pressure_pa,
        volume_velocity_m3_s,
        work,
    })
}

fn tone_invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn tone_step_error(error: PreparedStepError) -> AcousticRealizeError {
    match error {
        PreparedStepError::Cancelled => AcousticRealizeError::Cancelled,
        PreparedStepError::Solver(error) => render_step_error(error),
    }
}
