//! Fixed-clock evolution of the complete original nonlinear plate and material history.
use super::{AcousticRealizeError, VkBody, render_step_error};
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh, admit_midpoint_band};
use fs_phs::{PortHamiltonian, PreparedStepError, PreparedStepRecord, StepWorkspace};

/// Accepted compact acoustic observation and the original mechanical/material work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PreparedVkStep {
    /// Compact baffled-monopole pressure from the interval-mean acceleration [Pa].
    pub pressure_pa: f64,
    /// Accepted endpoint volume velocity [m^3/s].
    pub volume_velocity_m3_s: f64,
    /// Complete nonlinear storage, Maxwell/Kelvin heat and physical input work [J].
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

    /// Complete original `[q0,p0,q1,p1,...,material histories...]` state.
    #[must_use]
    pub fn state(&self) -> &[f64] {
        &self.x
    }

    /// Complete mechanical and material-memory storage [J].
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
        let mut pressure_pa = 0.0;
        let mut volume_velocity_m3_s = 0.0;
        for (k, area) in self.areas.iter().enumerate() {
            let velocity = self.candidate[2 * k + 1];
            let acceleration = (velocity - self.x[2 * k + 1]) / self.dt_s;
            pressure_pa += rho * area * acceleration / (core::f64::consts::TAU * listener_m);
            volume_velocity_m3_s += area * velocity;
        }
        if !pressure_pa.is_finite() || !volume_velocity_m3_s.is_finite() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "prepared plate acoustic observation is not finite",
            });
        }
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        core::mem::swap(&mut self.x, &mut self.candidate);
        Ok(PreparedVkStep {
            pressure_pa,
            volume_velocity_m3_s,
            work,
        })
    }
}
