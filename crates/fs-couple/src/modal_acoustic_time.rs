//! Time-domain realization of mass-normalized structural modes observed
//! through physical acoustic pressure transfers.
//!
//! This module is deliberately geometry- and material-agnostic. An upstream
//! structural solver supplies `omega` and a damping ratio, a force projector
//! supplies generalized force in `N / sqrt(kg)`, and an acoustic solver
//! supplies complex pressure per generalized modal velocity. The state then
//! has the standard mass-normalized units
//! `q [m sqrt(kg)]`, `qdot [m sqrt(kg) / s]`, and modal energy [J].
//!
//! Each sample advances the viscously damped oscillator exactly for a
//! zero-order-held force. The pressure observation
//! `p = Re(H) qdot + Im(H) omega q` is the exact real-valued realization of a
//! complex transfer `H` at that mode's natural frequency under the shared
//! `exp(-i omega t)` convention. It is therefore a narrow-band modal
//! radiation realization, not a broadband replacement for vector fitting an
//! acoustic frequency response.

use fs_math::{c64::C64, det};

/// Maximum number of modes admitted by one runtime model.
pub const MAX_TIME_DOMAIN_ACOUSTIC_MODES: usize = 4_096;

/// One mass-normalized damped mode and its observer transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalAcousticMode {
    /// Natural angular frequency [rad/s].
    pub angular_frequency_rad_s: f64,
    /// Viscous modal damping ratio `zeta`.
    pub damping_ratio: f64,
    /// Complex observer pressure per generalized modal velocity
    /// `[Pa s / (m sqrt(kg))]`, under `exp(-i omega t)`.
    pub pressure_per_modal_velocity: C64,
}

/// Explicit numerical and physical state limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalAcousticTimeBudget {
    /// Fraction of Nyquist below which every natural frequency must lie.
    pub nyquist_guard_fraction: f64,
    /// Maximum absolute mass-normalized displacement [m sqrt(kg)].
    pub maximum_abs_displacement_m_sqrt_kg: f64,
    /// Maximum absolute mass-normalized velocity [m sqrt(kg)/s].
    pub maximum_abs_velocity_m_sqrt_kg_per_s: f64,
    /// Maximum total modal energy [J].
    pub maximum_total_energy_j: f64,
    /// Maximum absolute physical observer pressure [Pa].
    pub maximum_abs_pressure_pa: f64,
}

impl ModalAcousticTimeBudget {
    /// Bounded reference values for audible structural modes. These are safety
    /// ceilings, not normalizations or mastering controls.
    #[must_use]
    pub const fn audible_reference() -> Self {
        Self {
            nyquist_guard_fraction: 0.9,
            maximum_abs_displacement_m_sqrt_kg: 1.0,
            maximum_abs_velocity_m_sqrt_kg_per_s: 100_000.0,
            maximum_total_energy_j: 1.0e9,
            maximum_abs_pressure_pa: 1.0e7,
        }
    }
}

/// State of one mass-normalized modal coordinate at a sample boundary.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModalAcousticState {
    /// Generalized modal displacement [m sqrt(kg)].
    pub displacement_m_sqrt_kg: f64,
    /// Generalized modal velocity [m sqrt(kg)/s].
    pub velocity_m_sqrt_kg_per_s: f64,
}

/// One exact zero-order-hold sample transition.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalAcousticFrame {
    /// Physical pressure at the observer after the transition [Pa].
    pub observer_pressure_pa: f64,
    /// Energy of each retained mode after the transition [J].
    pub modal_energy_j: Vec<f64>,
    /// Sum of retained modal energies after the transition [J].
    pub total_modal_energy_j: f64,
    /// Work done by the held generalized forces during this sample [J].
    pub input_work_j: f64,
    /// Energy removed by viscous damping during this sample [J].
    ///
    /// Values within `dissipation_roundoff_tolerance_j` of zero may be
    /// slightly negative because this is computed as `work - delta_energy`.
    pub viscous_dissipation_j: f64,
    /// Scale-aware roundoff allowance attached to the dissipation value [J].
    pub dissipation_roundoff_tolerance_j: f64,
}

/// Typed refusal from model admission or a transactional sample step.
#[derive(Clone, Debug, PartialEq)]
pub enum ModalAcousticTimeError {
    /// An input scalar or array violates its physical domain.
    InvalidInput {
        /// Failed physical or numerical invariant.
        what: &'static str,
    },
    /// A natural frequency violates the caller's Nyquist guard.
    ModeAboveNyquistGuard {
        /// Zero-based mode index.
        mode: usize,
        /// Natural frequency [Hz].
        frequency_hz: f64,
        /// Maximum admitted frequency [Hz].
        maximum_hz: f64,
    },
    /// The force vector cardinality differs from the admitted mode count.
    ForceCountMismatch {
        /// Admitted number of modes.
        expected: usize,
        /// Number of generalized forces supplied.
        found: usize,
    },
    /// A restored checkpoint supplied one state per the wrong number of modes.
    StateCountMismatch {
        /// Admitted number of modes.
        expected: usize,
        /// Number of supplied states.
        found: usize,
    },
    /// A state-dependent observer supplied one transfer per the wrong number
    /// of admitted modes.
    TransferCountMismatch {
        /// Admitted number of modes.
        expected: usize,
        /// Number of supplied complex transfers.
        found: usize,
    },
    /// A candidate step crossed a caller-owned state/energy/pressure ceiling.
    BudgetExceeded {
        /// Quantity whose absolute magnitude crossed its ceiling.
        what: &'static str,
        /// Candidate magnitude.
        value: f64,
        /// Caller-authored ceiling.
        limit: f64,
    },
    /// A supposedly dissipative oscillator produced energy beyond roundoff.
    NegativeDissipation {
        /// Computed work-minus-energy difference [J].
        dissipation_j: f64,
        /// Scale-aware floating-point allowance [J].
        tolerance_j: f64,
    },
}

impl core::fmt::Display for ModalAcousticTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { what } => {
                write!(f, "FS-COUPLE-MODAL-TIME-INPUT: {what}")
            }
            Self::ModeAboveNyquistGuard {
                mode,
                frequency_hz,
                maximum_hz,
            } => write!(
                f,
                "FS-COUPLE-MODAL-TIME-NYQUIST: mode {mode} frequency {frequency_hz:.6e} Hz exceeds {maximum_hz:.6e} Hz"
            ),
            Self::ForceCountMismatch { expected, found } => write!(
                f,
                "FS-COUPLE-MODAL-TIME-SHAPE: expected {expected} modal forces, found {found}"
            ),
            Self::StateCountMismatch { expected, found } => write!(
                f,
                "FS-COUPLE-MODAL-TIME-STATE-SHAPE: expected {expected} modal states, found {found}"
            ),
            Self::TransferCountMismatch { expected, found } => write!(
                f,
                "FS-COUPLE-MODAL-TIME-TRANSFER-SHAPE: expected {expected} observer transfers, found {found}"
            ),
            Self::BudgetExceeded { what, value, limit } => write!(
                f,
                "FS-COUPLE-MODAL-TIME-BUDGET: {what} {value:.6e} exceeds {limit:.6e}"
            ),
            Self::NegativeDissipation {
                dissipation_j,
                tolerance_j,
            } => write!(
                f,
                "FS-COUPLE-MODAL-TIME-PASSIVITY: dissipation {dissipation_j:.6e} J is below roundoff allowance {:.6e} J",
                -*tolerance_j
            ),
        }
    }
}

impl std::error::Error for ModalAcousticTimeError {}

/// Exact-ZOH runtime for independent mass-normalized modes.
#[derive(Clone, Debug)]
pub struct ModalAcousticTimeModel {
    sample_period_s: f64,
    modes: Vec<ModalAcousticMode>,
    states: Vec<ModalAcousticState>,
    budget: ModalAcousticTimeBudget,
}

impl ModalAcousticTimeModel {
    /// Admit one fixed-rate runtime model with all modes initially at rest.
    ///
    /// # Errors
    /// Refuses invalid sample rate, empty/oversized mode sets, non-finite or
    /// nonphysical mode parameters, unsafe Nyquist placement, or bad budgets.
    pub fn try_new(
        sample_rate_hz: u32,
        modes: Vec<ModalAcousticMode>,
        budget: ModalAcousticTimeBudget,
    ) -> Result<Self, ModalAcousticTimeError> {
        if sample_rate_hz == 0 {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "sample rate must be positive",
            });
        }
        if modes.is_empty() || modes.len() > MAX_TIME_DOMAIN_ACOUSTIC_MODES {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "mode count must be in 1..=4096",
            });
        }
        validate_budget(budget)?;
        let maximum_hz = 0.5 * f64::from(sample_rate_hz) * budget.nyquist_guard_fraction;
        for (mode_index, mode) in modes.iter().enumerate() {
            if !(mode.angular_frequency_rad_s > 0.0
                && mode.angular_frequency_rad_s.is_finite()
                && mode.damping_ratio >= 0.0
                && mode.damping_ratio.is_finite()
                && mode.pressure_per_modal_velocity.re.is_finite()
                && mode.pressure_per_modal_velocity.im.is_finite())
            {
                return Err(ModalAcousticTimeError::InvalidInput {
                    what: "modes need positive finite frequency, non-negative finite damping, and finite pressure transfer",
                });
            }
            let frequency_hz = mode.angular_frequency_rad_s / core::f64::consts::TAU;
            if frequency_hz > maximum_hz {
                return Err(ModalAcousticTimeError::ModeAboveNyquistGuard {
                    mode: mode_index,
                    frequency_hz,
                    maximum_hz,
                });
            }
        }
        let states = vec![ModalAcousticState::default(); modes.len()];
        Ok(Self {
            sample_period_s: f64::from(sample_rate_hz).recip(),
            modes,
            states,
            budget,
        })
    }

    /// Borrow the exact admitted modes.
    #[must_use]
    pub fn modes(&self) -> &[ModalAcousticMode] {
        &self.modes
    }

    /// Borrow the current sample-boundary states.
    #[must_use]
    pub fn states(&self) -> &[ModalAcousticState] {
        &self.states
    }

    /// Replace the complete modal state from an accepted external checkpoint.
    ///
    /// This is the transactional restart seam used when another coupled
    /// physics owner, rather than the acoustic runtime itself, owns commit and
    /// rollback. The candidate must match the admitted modal basis and the
    /// same displacement, velocity, and total-energy budgets as a normal
    /// advance. A refusal leaves the runtime unchanged.
    pub fn restore_states(
        &mut self,
        states: &[ModalAcousticState],
    ) -> Result<(), ModalAcousticTimeError> {
        if states.len() != self.modes.len() {
            return Err(ModalAcousticTimeError::StateCountMismatch {
                expected: self.modes.len(),
                found: states.len(),
            });
        }
        let mut total_energy_j = 0.0;
        for (mode, state) in self.modes.iter().zip(states) {
            if !(state.displacement_m_sqrt_kg.is_finite()
                && state.velocity_m_sqrt_kg_per_s.is_finite())
            {
                return Err(ModalAcousticTimeError::InvalidInput {
                    what: "restored modal states must be finite",
                });
            }
            check_limit(
                "absolute modal displacement",
                state.displacement_m_sqrt_kg.abs(),
                self.budget.maximum_abs_displacement_m_sqrt_kg,
            )?;
            check_limit(
                "absolute modal velocity",
                state.velocity_m_sqrt_kg_per_s.abs(),
                self.budget.maximum_abs_velocity_m_sqrt_kg_per_s,
            )?;
            total_energy_j += modal_energy(*mode, *state);
        }
        if !total_energy_j.is_finite() {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "restored modal energy must be finite",
            });
        }
        check_limit(
            "total modal energy",
            total_energy_j,
            self.budget.maximum_total_energy_j,
        )?;
        self.states.clone_from_slice(states);
        Ok(())
    }

    /// Nominal output sample period [s].
    #[must_use]
    pub const fn sample_period_s(&self) -> f64 {
        self.sample_period_s
    }

    /// Initialize every mode at the static equilibrium of one held force.
    ///
    /// This is the causal initial condition for a simulation window that
    /// begins after a load has already been applied for long enough that its
    /// transient vibration has decayed. Each mass-normalized coordinate is set
    /// to `q = force / omega^2` with zero velocity. It must not be used to hide
    /// a real load-on event; callers that model application of a new force
    /// should retain the zero state and advance that force normally.
    ///
    /// The update is transactional: malformed forces or a state/energy budget
    /// refusal leave every modal coordinate unchanged.
    pub fn initialize_static_equilibrium(
        &mut self,
        generalized_force_n_per_sqrt_kg: &[f64],
    ) -> Result<(), ModalAcousticTimeError> {
        if generalized_force_n_per_sqrt_kg.len() != self.modes.len() {
            return Err(ModalAcousticTimeError::ForceCountMismatch {
                expected: self.modes.len(),
                found: generalized_force_n_per_sqrt_kg.len(),
            });
        }
        let mut candidate = Vec::with_capacity(self.states.len());
        let mut total_energy_j = 0.0;
        for (mode, force) in self.modes.iter().zip(generalized_force_n_per_sqrt_kg) {
            if !force.is_finite() {
                return Err(ModalAcousticTimeError::InvalidInput {
                    what: "static-equilibrium modal forces must be finite",
                });
            }
            let omega = mode.angular_frequency_rad_s;
            let state = ModalAcousticState {
                displacement_m_sqrt_kg: force / (omega * omega),
                velocity_m_sqrt_kg_per_s: 0.0,
            };
            check_limit(
                "absolute modal displacement",
                state.displacement_m_sqrt_kg.abs(),
                self.budget.maximum_abs_displacement_m_sqrt_kg,
            )?;
            total_energy_j += modal_energy(*mode, state);
            candidate.push(state);
        }
        if !total_energy_j.is_finite() {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "static-equilibrium modal energy must be finite",
            });
        }
        check_limit(
            "total modal energy",
            total_energy_j,
            self.budget.maximum_total_energy_j,
        )?;
        self.states = candidate;
        Ok(())
    }

    /// Observe the current modal state through caller-supplied complex
    /// pressure transfers.
    ///
    /// This is the generic seam for a moving source/listener or a changing
    /// acoustic medium. Each transfer has the same units and
    /// `exp(-i omega t)` convention as [`ModalAcousticMode::pressure_per_modal_velocity`].
    /// The method does not advance or mutate oscillator state.
    ///
    /// # Errors
    /// Refuses wrong cardinality, nonfinite transfers/pressure, or the
    /// admitted absolute-pressure ceiling.
    pub fn observer_pressure_with_transfers(
        &self,
        pressure_per_modal_velocity: &[C64],
    ) -> Result<f64, ModalAcousticTimeError> {
        if pressure_per_modal_velocity.len() != self.modes.len() {
            return Err(ModalAcousticTimeError::TransferCountMismatch {
                expected: self.modes.len(),
                found: pressure_per_modal_velocity.len(),
            });
        }
        let mut pressure_pa = 0.0;
        for ((mode, state), transfer) in self
            .modes
            .iter()
            .zip(&self.states)
            .zip(pressure_per_modal_velocity)
        {
            if !(transfer.re.is_finite() && transfer.im.is_finite()) {
                return Err(ModalAcousticTimeError::InvalidInput {
                    what: "observer pressure transfers must be finite",
                });
            }
            pressure_pa += transfer.re * state.velocity_m_sqrt_kg_per_s
                + transfer.im * mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg;
        }
        if !pressure_pa.is_finite() {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "observer pressure evaluation produced a non-finite result",
            });
        }
        check_limit(
            "absolute observer pressure",
            pressure_pa.abs(),
            self.budget.maximum_abs_pressure_pa,
        )?;
        Ok(pressure_pa)
    }

    /// Observe oscillatory pressure about the static equilibrium of the
    /// currently held generalized force.
    ///
    /// A complex narrow-band transfer represents the quadrature pair of a
    /// mode oscillating near its natural frequency. Applying its imaginary
    /// part directly to the total displacement would incorrectly turn the
    /// static compliance `q_static = force / omega^2` into DC acoustic
    /// pressure. This variant removes that equilibrium displacement before
    /// applying the quadrature term while preserving the actual modal
    /// velocity. It is appropriate when the structural step is driven by the
    /// same zero-order-held force supplied here.
    ///
    /// The method remains a narrow-band modal radiation realization. It does
    /// not turn one-frequency transfers into a broadband retarded-time solve.
    ///
    /// # Errors
    /// Refuses wrong transfer/force cardinality, nonfinite inputs or pressure,
    /// or the admitted absolute-pressure ceiling.
    pub fn observer_pressure_with_transfers_about_static_equilibrium(
        &self,
        pressure_per_modal_velocity: &[C64],
        held_generalized_force_n_per_sqrt_kg: &[f64],
    ) -> Result<f64, ModalAcousticTimeError> {
        if pressure_per_modal_velocity.len() != self.modes.len() {
            return Err(ModalAcousticTimeError::TransferCountMismatch {
                expected: self.modes.len(),
                found: pressure_per_modal_velocity.len(),
            });
        }
        if held_generalized_force_n_per_sqrt_kg.len() != self.modes.len() {
            return Err(ModalAcousticTimeError::ForceCountMismatch {
                expected: self.modes.len(),
                found: held_generalized_force_n_per_sqrt_kg.len(),
            });
        }
        let mut pressure_pa = 0.0;
        for (((mode, state), transfer), force) in self
            .modes
            .iter()
            .zip(&self.states)
            .zip(pressure_per_modal_velocity)
            .zip(held_generalized_force_n_per_sqrt_kg)
        {
            if !(transfer.re.is_finite() && transfer.im.is_finite() && force.is_finite()) {
                return Err(ModalAcousticTimeError::InvalidInput {
                    what: "observer transfers and held modal forces must be finite",
                });
            }
            let omega = mode.angular_frequency_rad_s;
            let dynamic_displacement = state.displacement_m_sqrt_kg - force / (omega * omega);
            pressure_pa += transfer.re * state.velocity_m_sqrt_kg_per_s
                + transfer.im * omega * dynamic_displacement;
        }
        if !pressure_pa.is_finite() {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "equilibrium-relative observer pressure became non-finite",
            });
        }
        check_limit(
            "absolute observer pressure",
            pressure_pa.abs(),
            self.budget.maximum_abs_pressure_pa,
        )?;
        Ok(pressure_pa)
    }

    /// Advance all modes by one exact zero-order-held sample.
    ///
    /// The step is transactional: a refusal leaves every state unchanged.
    ///
    /// # Errors
    /// Refuses a wrong-sized/non-finite force vector, a budget crossing, or a
    /// negative dissipation defect beyond scale-aware roundoff.
    pub fn step(
        &mut self,
        generalized_force_n_per_sqrt_kg: &[f64],
    ) -> Result<ModalAcousticFrame, ModalAcousticTimeError> {
        self.step_duration(generalized_force_n_per_sqrt_kg, self.sample_period_s)
    }

    /// Advance all modes for one positive caller-specified subinterval under
    /// exact zero-order-held generalized forces.
    ///
    /// This seam lets a fixed-rate audio producer split one sample at exact
    /// mechanics-control boundaries without smearing force changes. The
    /// nominal sample rate still owns the Nyquist admission gate.
    ///
    /// # Errors
    /// Refuses nonpositive/nonfinite duration and every condition documented
    /// by [`Self::step`]. A refusal leaves every modal state unchanged.
    pub fn step_duration(
        &mut self,
        generalized_force_n_per_sqrt_kg: &[f64],
        duration_s: f64,
    ) -> Result<ModalAcousticFrame, ModalAcousticTimeError> {
        if !(duration_s > 0.0 && duration_s.is_finite()) {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "step duration must be positive and finite",
            });
        }
        if generalized_force_n_per_sqrt_kg.len() != self.modes.len() {
            return Err(ModalAcousticTimeError::ForceCountMismatch {
                expected: self.modes.len(),
                found: generalized_force_n_per_sqrt_kg.len(),
            });
        }
        if generalized_force_n_per_sqrt_kg
            .iter()
            .any(|force| !force.is_finite())
        {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "generalized forces must be finite",
            });
        }

        let mut candidate = Vec::with_capacity(self.states.len());
        let mut modal_energy_j = Vec::with_capacity(self.states.len());
        let mut pressure_pa = 0.0;
        let mut input_work_j = 0.0;
        let mut energy_before_j = 0.0;
        let mut energy_after_j = 0.0;
        for ((mode, state), force) in self
            .modes
            .iter()
            .zip(&self.states)
            .zip(generalized_force_n_per_sqrt_kg)
        {
            let next = advance_exact_zoh(*mode, *state, *force, duration_s);
            check_limit(
                "absolute modal displacement",
                next.displacement_m_sqrt_kg.abs(),
                self.budget.maximum_abs_displacement_m_sqrt_kg,
            )?;
            check_limit(
                "absolute modal velocity",
                next.velocity_m_sqrt_kg_per_s.abs(),
                self.budget.maximum_abs_velocity_m_sqrt_kg_per_s,
            )?;
            let before = modal_energy(*mode, *state);
            let after = modal_energy(*mode, next);
            energy_before_j += before;
            energy_after_j += after;
            input_work_j += force * (next.displacement_m_sqrt_kg - state.displacement_m_sqrt_kg);
            pressure_pa += mode.pressure_per_modal_velocity.re * next.velocity_m_sqrt_kg_per_s
                + mode.pressure_per_modal_velocity.im
                    * mode.angular_frequency_rad_s
                    * next.displacement_m_sqrt_kg;
            modal_energy_j.push(after);
            candidate.push(next);
        }
        check_limit(
            "total modal energy",
            energy_after_j,
            self.budget.maximum_total_energy_j,
        )?;
        check_limit(
            "absolute observer pressure",
            pressure_pa.abs(),
            self.budget.maximum_abs_pressure_pa,
        )?;
        if !pressure_pa.is_finite()
            || !energy_before_j.is_finite()
            || !energy_after_j.is_finite()
            || !input_work_j.is_finite()
        {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "sample transition produced a non-finite result",
            });
        }
        let viscous_dissipation_j = input_work_j - (energy_after_j - energy_before_j);
        let scale = input_work_j
            .abs()
            .max(energy_before_j)
            .max(energy_after_j)
            .max(f64::MIN_POSITIVE);
        let dissipation_roundoff_tolerance_j = 256.0 * f64::EPSILON * scale;
        if viscous_dissipation_j < -dissipation_roundoff_tolerance_j {
            return Err(ModalAcousticTimeError::NegativeDissipation {
                dissipation_j: viscous_dissipation_j,
                tolerance_j: dissipation_roundoff_tolerance_j,
            });
        }
        self.states = candidate;
        Ok(ModalAcousticFrame {
            observer_pressure_pa: pressure_pa,
            modal_energy_j,
            total_modal_energy_j: energy_after_j,
            input_work_j,
            viscous_dissipation_j,
            dissipation_roundoff_tolerance_j,
        })
    }
}

fn validate_budget(budget: ModalAcousticTimeBudget) -> Result<(), ModalAcousticTimeError> {
    if !(budget.nyquist_guard_fraction > 0.0
        && budget.nyquist_guard_fraction < 1.0
        && budget.nyquist_guard_fraction.is_finite())
    {
        return Err(ModalAcousticTimeError::InvalidInput {
            what: "Nyquist guard fraction must be finite and in (0,1)",
        });
    }
    for (value, what) in [
        (
            budget.maximum_abs_displacement_m_sqrt_kg,
            "maximum displacement must be positive and finite",
        ),
        (
            budget.maximum_abs_velocity_m_sqrt_kg_per_s,
            "maximum velocity must be positive and finite",
        ),
        (
            budget.maximum_total_energy_j,
            "maximum energy must be positive and finite",
        ),
        (
            budget.maximum_abs_pressure_pa,
            "maximum pressure must be positive and finite",
        ),
    ] {
        if !(value > 0.0 && value.is_finite()) {
            return Err(ModalAcousticTimeError::InvalidInput { what });
        }
    }
    Ok(())
}

fn check_limit(what: &'static str, value: f64, limit: f64) -> Result<(), ModalAcousticTimeError> {
    if value > limit {
        return Err(ModalAcousticTimeError::BudgetExceeded { what, value, limit });
    }
    Ok(())
}

fn modal_energy(mode: ModalAcousticMode, state: ModalAcousticState) -> f64 {
    0.5 * (state.velocity_m_sqrt_kg_per_s * state.velocity_m_sqrt_kg_per_s
        + mode.angular_frequency_rad_s
            * mode.angular_frequency_rad_s
            * state.displacement_m_sqrt_kg
            * state.displacement_m_sqrt_kg)
}

fn advance_exact_zoh(
    mode: ModalAcousticMode,
    state: ModalAcousticState,
    force: f64,
    dt: f64,
) -> ModalAcousticState {
    let omega = mode.angular_frequency_rad_s;
    let zeta = mode.damping_ratio;
    let equilibrium = force / (omega * omega);
    let q0 = state.displacement_m_sqrt_kg - equilibrium;
    let v0 = state.velocity_m_sqrt_kg_per_s;
    let critical_delta = zeta - 1.0;
    let (q, v) = if critical_delta.abs() <= 1.0e-8 {
        let decay = det::exp(-omega * dt);
        let omega_dt = omega * dt;
        (
            decay * ((1.0 + omega_dt) * q0 + dt * v0),
            decay * (-omega * omega * dt * q0 + (1.0 - omega_dt) * v0),
        )
    } else if zeta < 1.0 {
        let root = det::sqrt(1.0 - zeta * zeta);
        let damped_omega = omega * root;
        let angle = damped_omega * dt;
        let sine = det::sin(angle);
        let cosine = det::cos(angle);
        let decay = det::exp(-zeta * omega * dt);
        let damping_over_damped = zeta / root;
        (
            decay * ((cosine + damping_over_damped * sine) * q0 + sine / damped_omega * v0),
            decay * (-omega / root * sine * q0 + (cosine - damping_over_damped * sine) * v0),
        )
    } else {
        let root = det::sqrt(zeta * zeta - 1.0);
        let slow = -omega * (zeta - root);
        let fast = -omega * (zeta + root);
        let denominator = slow - fast;
        let slow_amplitude = (v0 - fast * q0) / denominator;
        let fast_amplitude = (slow * q0 - v0) / denominator;
        let slow_e = det::exp(slow * dt);
        let fast_e = det::exp(fast * dt);
        (
            slow_amplitude * slow_e + fast_amplitude * fast_e,
            slow * slow_amplitude * slow_e + fast * fast_amplitude * fast_e,
        )
    };
    ModalAcousticState {
        displacement_m_sqrt_kg: q + equilibrium,
        velocity_m_sqrt_kg_per_s: v,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn budget() -> ModalAcousticTimeBudget {
        ModalAcousticTimeBudget {
            maximum_abs_displacement_m_sqrt_kg: 100.0,
            maximum_abs_velocity_m_sqrt_kg_per_s: 100_000.0,
            maximum_total_energy_j: 1.0e12,
            maximum_abs_pressure_pa: 1.0e12,
            ..ModalAcousticTimeBudget::audible_reference()
        }
    }

    #[test]
    fn g1_undamped_constant_force_matches_closed_form_and_conserves_work_energy() {
        let omega = 2.0 * core::f64::consts::PI * 1_000.0;
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::new(2.0, -0.25),
        };
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
        let force = 3.0;
        let frame = model.step(&[force]).unwrap();
        let dt = 1.0 / 48_000.0;
        let expected_q = force / (omega * omega) * (1.0 - det::cos(omega * dt));
        let expected_v = force / omega * det::sin(omega * dt);
        assert!((model.states()[0].displacement_m_sqrt_kg - expected_q).abs() < 1.0e-18);
        assert!((model.states()[0].velocity_m_sqrt_kg_per_s - expected_v).abs() < 1.0e-15);
        assert!((frame.input_work_j - frame.total_modal_energy_j).abs() < 1.0e-15);
        assert!(frame.viscous_dissipation_j.abs() <= frame.dissipation_roundoff_tolerance_j);
        let expected_pressure = 2.0 * expected_v - 0.25 * omega * expected_q;
        assert!((frame.observer_pressure_pa - expected_pressure).abs() < 1.0e-15);
    }

    #[test]
    fn g0_damped_free_decay_is_passive_in_all_three_damping_regimes() {
        for damping_ratio in [0.02, 1.0, 2.0] {
            let mode = ModalAcousticMode {
                angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 400.0,
                damping_ratio,
                pressure_per_modal_velocity: C64::from_re(1.0),
            };
            let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
            model.states[0].displacement_m_sqrt_kg = 1.0e-3;
            let before = modal_energy(mode, model.states[0]);
            let frame = model.step(&[0.0]).unwrap();
            assert!(frame.total_modal_energy_j < before);
            assert!(frame.viscous_dissipation_j >= -frame.dissipation_roundoff_tolerance_j);
        }
    }

    #[test]
    fn transactional_budget_refusal_leaves_state_unchanged() {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 1_000.0,
            damping_ratio: 0.01,
            pressure_per_modal_velocity: C64::from_re(1.0),
        };
        let mut tight = budget();
        tight.maximum_abs_velocity_m_sqrt_kg_per_s = 1.0e-12;
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], tight).unwrap();
        let before = model.states().to_vec();
        assert!(matches!(
            model.step(&[1.0]),
            Err(ModalAcousticTimeError::BudgetExceeded { .. })
        ));
        assert_eq!(model.states(), before);
    }

    #[test]
    fn g3_splitting_a_held_force_at_a_control_boundary_preserves_the_solution() {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 2_400.0,
            damping_ratio: 0.07,
            pressure_per_modal_velocity: C64::new(1.5, -0.4),
        };
        let mut whole = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
        let mut split = whole.clone();
        let whole_frame = whole.step(&[2.5]).unwrap();
        split.step_duration(&[2.5], 0.5 / 48_000.0).unwrap();
        let split_frame = split.step_duration(&[2.5], 0.5 / 48_000.0).unwrap();
        assert!(
            (whole.states()[0].displacement_m_sqrt_kg - split.states()[0].displacement_m_sqrt_kg)
                .abs()
                < 1.0e-18
        );
        assert!(
            (whole.states()[0].velocity_m_sqrt_kg_per_s
                - split.states()[0].velocity_m_sqrt_kg_per_s)
                .abs()
                < 1.0e-13
        );
        assert!(
            (whole_frame.observer_pressure_pa - split_frame.observer_pressure_pa).abs() < 1.0e-13
        );
    }

    #[test]
    fn g0_state_dependent_observer_transfer_is_exact_and_nonmutating() {
        let omega = 2.0 * core::f64::consts::PI * 1_200.0;
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: 0.03,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
        model.step(&[4.0]).unwrap();
        let before = model.states().to_vec();
        let transfer = C64::new(1.75, -0.3);
        let pressure = model.observer_pressure_with_transfers(&[transfer]).unwrap();
        let expected = transfer.re * before[0].velocity_m_sqrt_kg_per_s
            + transfer.im * omega * before[0].displacement_m_sqrt_kg;
        assert!((pressure - expected).abs() < 1.0e-15);
        assert_eq!(model.states(), before);

        assert!(matches!(
            model.observer_pressure_with_transfers(&[]),
            Err(ModalAcousticTimeError::TransferCountMismatch {
                expected: 1,
                found: 0
            })
        ));
        assert!(matches!(
            model.observer_pressure_with_transfers(&[C64::new(f64::NAN, 0.0)]),
            Err(ModalAcousticTimeError::InvalidInput { .. })
        ));
        assert_eq!(model.states(), before);
    }

    #[test]
    fn g0_static_equilibrium_does_not_radiate_through_narrow_band_quadrature() {
        let omega = 2.0 * core::f64::consts::PI * 1_500.0;
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: 0.02,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
        let force = 7.0;
        model.states[0] = ModalAcousticState {
            displacement_m_sqrt_kg: force / (omega * omega),
            velocity_m_sqrt_kg_per_s: 0.0,
        };
        let transfer = C64::new(2.0, -0.75);
        let legacy = model.observer_pressure_with_transfers(&[transfer]).unwrap();
        let pressure = model
            .observer_pressure_with_transfers_about_static_equilibrium(&[transfer], &[force])
            .unwrap();
        assert!(legacy.abs() > 0.0);
        assert_eq!(pressure.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            model
                .observer_pressure_with_transfers_about_static_equilibrium(&[transfer], &[0.0])
                .unwrap()
                .to_bits(),
            legacy.to_bits()
        );
    }

    #[test]
    fn g0_static_equilibrium_initializer_sets_exact_compliance_transactionally() {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 1_500.0,
            damping_ratio: 0.02,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
        let force = 3.25;
        model.initialize_static_equilibrium(&[force]).unwrap();
        let omega = model.modes()[0].angular_frequency_rad_s;
        assert_eq!(
            model.states()[0],
            ModalAcousticState {
                displacement_m_sqrt_kg: force / (omega * omega),
                velocity_m_sqrt_kg_per_s: 0.0,
            }
        );

        let before = model.states().to_vec();
        assert!(matches!(
            model.initialize_static_equilibrium(&[f64::NAN]),
            Err(ModalAcousticTimeError::InvalidInput { .. })
        ));
        assert_eq!(model.states(), before);
    }

    #[test]
    fn g0_external_checkpoint_restore_is_exact_and_transactional() {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 900.0,
            damping_ratio: 0.015,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
        let restored = [ModalAcousticState {
            displacement_m_sqrt_kg: 2.0e-5,
            velocity_m_sqrt_kg_per_s: -0.03,
        }];
        model.restore_states(&restored).unwrap();
        assert_eq!(model.states(), restored);

        let before = model.states().to_vec();
        assert!(matches!(
            model.restore_states(&[]),
            Err(ModalAcousticTimeError::StateCountMismatch {
                expected: 1,
                found: 0
            })
        ));
        assert_eq!(model.states(), before);
        assert!(matches!(
            model.restore_states(&[ModalAcousticState {
                displacement_m_sqrt_kg: f64::NAN,
                velocity_m_sqrt_kg_per_s: 0.0,
            }]),
            Err(ModalAcousticTimeError::InvalidInput { .. })
        ));
        assert_eq!(model.states(), before);
    }
}
