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
//!
//! Explicit free coordinates use the exact zero-stiffness held-force limit,
//! not a small artificial natural frequency. They carry kinetic energy and
//! exchange physical forces through the same modal ports, but have no direct
//! narrow-band pressure transfer. This is translation in a supplied fixed
//! coordinate, not a rotating six-degree-of-freedom rigid-body model.
//!
//! Owner: `fs-modal` (layer L1), moved here from `fs-couple` under owner decision OD10 (IRP M1.1-W1a) so that strings, solids
//! and couplers share one acyclic owner. `fs_couple::modal_acoustic_time` re-exports this module under every former name and
//! keeps the opt-in `broadband` observer, which needs `fs-vfit`, `fs-phs` and the broadband radiation owner.

use std::sync::Arc;

use fs_math::{c64::C64, det};

#[path = "time/heat.rs"]
mod heat;
use heat::{DampedTransition, ViscousIntegral};

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
    /// Independently integrates `2*zeta*omega*v(t)^2` along the exact held-force
    /// trajectory. This is nonnegative physical heat, including losses smaller
    /// than subtraction roundoff in the stored mechanical energy.
    pub viscous_dissipation_j: f64,
    /// Independent mode-local viscous integrals for physical loss routing.
    /// Their sum is the aggregate heat; undamped and free modes report exact zero.
    pub modal_viscous_losses: Vec<ModalViscousLoss>,
    /// Scale-aware roundoff allowance attached to the dissipation value [J].
    pub dissipation_roundoff_tolerance_j: f64,
}

/// Viscous energy removal for one retained mode, with numerical uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalViscousLoss {
    /// Nonnegative integral of the physical viscous power [J].
    /// No subtraction of stored mechanical energies enters this value.
    pub energy_j: f64,
    /// Scale-aware floating-point allowance for the integral contraction [J].
    pub roundoff_tolerance_j: f64,
}

/// Reusable sample-transition storage, allocated outside the audio callback.
///
/// One workspace serves any model with the same mode count; it stores no
/// constitutive coefficients or history. Successful `step_into` calls return
/// its borrowed diagnostics. A refusal may overwrite scratch entries but never
/// publishes a frame or changes the model's accepted state.
#[derive(Clone, Debug)]
pub struct ModalAcousticWorkspace {
    candidate: Vec<ModalAcousticState>,
    frame: ModalAcousticFrame,
}

impl ModalAcousticWorkspace {
    /// Allocate all candidate-state and diagnostic buffers for an admitted model.
    #[must_use]
    pub fn new(model: &ModalAcousticTimeModel) -> Self {
        let count = model.modes.len();
        Self {
            candidate: vec![ModalAcousticState::default(); count],
            frame: ModalAcousticFrame {
                observer_pressure_pa: 0.0,
                modal_energy_j: vec![0.0; count],
                total_modal_energy_j: 0.0,
                input_work_j: 0.0,
                viscous_dissipation_j: 0.0,
                modal_viscous_losses: vec![
                    ModalViscousLoss {
                        energy_j: 0.0,
                        roundoff_tolerance_j: 0.0,
                    };
                    count
                ],
                dissipation_roundoff_tolerance_j: 0.0,
            },
        }
    }

    /// Diagnostics published by the last successful `step_into` (all zero before the first).
    #[must_use]
    pub fn frame(&self) -> &ModalAcousticFrame {
        &self.frame
    }

    /// Number of modes this workspace can stage without resizing.
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.candidate.len()
    }
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
    transition_kernel: Arc<ExactZohKernel>,
}

#[derive(Debug)]
struct ExactZohKernel {
    sample_period_s: f64,
    full_duration_bits: u64,
    half_duration_bits: u64,
    full: Box<[PreparedZohCoefficients]>,
    half: std::sync::OnceLock<Box<[PreparedZohCoefficients]>>,
}

impl ExactZohKernel {
    fn new(modes: &[ModalAcousticMode], sample_period_s: f64) -> Self {
        let half_period_s = 0.5 * sample_period_s;
        let full = modes
            .iter()
            .copied()
            .map(|mode| PreparedZohCoefficients::new(mode, sample_period_s))
            .collect::<Vec<_>>()
            .into_boxed_slice();
        Self {
            sample_period_s,
            full_duration_bits: sample_period_s.to_bits(),
            half_duration_bits: half_period_s.to_bits(),
            full,
            half: std::sync::OnceLock::new(),
        }
    }

    fn for_duration(
        &self,
        modes: &[ModalAcousticMode],
        duration_s: f64,
    ) -> Option<&[PreparedZohCoefficients]> {
        let duration_bits = duration_s.to_bits();
        if duration_bits == self.full_duration_bits {
            Some(&self.full)
        } else if duration_bits == self.half_duration_bits {
            let half_period_s = 0.5 * self.sample_period_s;
            Some(self.half.get_or_init(|| {
                modes
                    .iter()
                    .copied()
                    .map(|mode| PreparedZohCoefficients::new(mode, half_period_s))
                    .collect::<Vec<_>>()
                    .into_boxed_slice()
            }))
        } else {
            None
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct PreparedZohCoefficients {
    transition: ExactZohCoefficients,
    viscous_integral: ViscousIntegral,
}
impl PreparedZohCoefficients {
    fn new(mode: ModalAcousticMode, dt: f64) -> Self {
        Self {
            transition: exact_zoh_coefficients(mode, dt),
            viscous_integral: ViscousIntegral::new(mode, dt),
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum ExactZohCoefficients {
    Free {
        duration_s: f64,
    },
    Undamped {
        sine: f64,
        one_minus_cosine: f64,
    },
    Damped {
        position_retention: f64,
        velocity_retention: f64,
        velocity_gain: f64,
        displacement_gain: f64,
        position_decrement: f64,
        velocity_decrement: f64,
    },
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
        Self::try_new_impl(sample_rate_hz, modes, budget, false)
    }

    /// Explicitly admit zero-stiffness, mass-normalized coordinates alongside
    /// ordinary modes. A zero frequency requires zero damping ratio and zero
    /// pressure transfer: neither a drag law nor free-body radiation is inferred.
    /// Positive-frequency admission and transition arithmetic are unchanged.
    /// Initial coordinates are zero; use `restore_states` for supplied motion.
    /// Static-compliance initialization refuses models containing free coordinates.
    ///
    /// # Errors
    /// All `try_new` limits, or a free coordinate with damping/acoustic transfer.
    pub fn try_new_with_free_coordinates(
        sample_rate_hz: u32,
        modes: Vec<ModalAcousticMode>,
        budget: ModalAcousticTimeBudget,
    ) -> Result<Self, ModalAcousticTimeError> {
        Self::try_new_impl(sample_rate_hz, modes, budget, true)
    }

    /// Construct one freely translating mass from kilograms, metres and m/s.
    /// Its coordinate is q=sqrt(m)*x, velocity=sqrt(m)*dx/dt, so a force F [N]
    /// enters as F/sqrt(m), and the conjugate displacement/force-port shape is
    /// 1/sqrt(m). The same transformation must be used in attachment maps.
    /// There is no tether, drag, imposed restitution or direct acoustic output.
    /// The existing budgets retain their mass-normalized units and joules.
    ///
    /// # Errors
    /// Nonpositive/nonfinite mass, nonfinite motion, or clock/state/budget limits.
    pub fn try_free_mass(
        sample_rate_hz: u32,
        mass_kg: f64,
        displacement_m: f64,
        velocity_m_s: f64,
        budget: ModalAcousticTimeBudget,
    ) -> Result<Self, ModalAcousticTimeError> {
        if !mass_kg.is_finite()
            || mass_kg <= 0.0
            || !displacement_m.is_finite()
            || !velocity_m_s.is_finite()
        {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "free mass requires positive finite kilograms and finite physical motion",
            });
        }
        let root_mass = det::sqrt(mass_kg);
        let mut model = Self::try_new_with_free_coordinates(
            sample_rate_hz,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            }],
            budget,
        )?;
        model.restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: root_mass * displacement_m,
            velocity_m_sqrt_kg_per_s: root_mass * velocity_m_s,
        }])?;
        Ok(model)
    }

    fn try_new_impl(
        sample_rate_hz: u32,
        modes: Vec<ModalAcousticMode>,
        budget: ModalAcousticTimeBudget,
        allow_free_coordinates: bool,
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
            if allow_free_coordinates && mode.angular_frequency_rad_s == 0.0 {
                if mode.damping_ratio != 0.0 || mode.pressure_per_modal_velocity != C64::ZERO {
                    return Err(ModalAcousticTimeError::InvalidInput {
                        what: "free coordinates require zero damping ratio and zero pressure transfer",
                    });
                }
            } else if !(mode.angular_frequency_rad_s > 0.0
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
        let sample_period_s = f64::from(sample_rate_hz).recip();
        let transition_kernel = Arc::new(ExactZohKernel::new(&modes, sample_period_s));
        Ok(Self {
            sample_period_s,
            modes,
            states,
            budget,
            transition_kernel,
        })
    }

    /// Borrow the exact admitted modes.
    #[must_use]
    pub fn modes(&self) -> &[ModalAcousticMode] {
        &self.modes
    }

    /// A copy of this model with every narrow-band pressure transfer set to zero. This read-only transfer edit changes no
    /// mechanical coefficient, cached transition or state; it exists so that a broadband receiver can own the pressure ceiling.
    #[must_use]
    pub fn without_pressure_transfers(&self) -> Self {
        let mut copy = self.clone();
        for mode in &mut copy.modes {
            mode.pressure_per_modal_velocity = C64::ZERO;
        }
        copy
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

    #[inline]
    pub(crate) fn copy_states_from(&mut self, source: &Self) {
        debug_assert_eq!(self.states.len(), source.states.len());
        self.states.copy_from_slice(&source.states);
    }

    /// One exact zero-order-hold step of mode `mode_index` from `state` under the held generalized `force`, with this model's
    /// cached full-step coefficients. Public so that single-mode consumers (the bow contact, compact radiators and the
    /// broadband observer in `fs-couple`) take the very same step the model takes.
    #[inline]
    pub fn advance_exact_zoh_step(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        force: f64,
    ) -> ModalAcousticState {
        advance_exact_zoh_cached(
            self.modes[mode_index],
            state,
            force,
            self.transition_kernel.full[mode_index].transition,
        )
    }

    /// Exact zero-order-hold one-sample linear map `(qq, qv, vq, vv, bq, bv)`
    /// of `self.modes()[mode_index]`, evaluated on the three basis states
    /// `(q = 1/omega, v = 0, F = 0)` (or `q = 1` when `omega == 0`),
    /// `(q = 0, v = 1, F = 0)`, and `(q = 0, v = 0, F = 1)` using the admitted
    /// model's cached full-step transition coefficients.
    #[must_use]
    pub fn probe_linear_map(&self, mode_index: usize) -> (f64, f64, f64, f64, f64, f64) {
        let omega = self.modes[mode_index].angular_frequency_rad_s;
        let (qq, vq) = if omega == 0.0 {
            let s_q = self.advance_exact_zoh_step(
                mode_index,
                ModalAcousticState {
                    displacement_m_sqrt_kg: 1.0,
                    velocity_m_sqrt_kg_per_s: 0.0,
                },
                0.0,
            );
            (s_q.displacement_m_sqrt_kg, s_q.velocity_m_sqrt_kg_per_s)
        } else {
            let s_q = self.advance_exact_zoh_step(
                mode_index,
                ModalAcousticState {
                    displacement_m_sqrt_kg: 1.0 / omega,
                    velocity_m_sqrt_kg_per_s: 0.0,
                },
                0.0,
            );
            (
                s_q.displacement_m_sqrt_kg * omega,
                s_q.velocity_m_sqrt_kg_per_s * omega,
            )
        };
        let s_v = self.advance_exact_zoh_step(
            mode_index,
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: 1.0,
            },
            0.0,
        );
        let s_f = self.advance_exact_zoh_step(mode_index, ModalAcousticState::default(), 1.0);
        (
            qq,
            s_v.displacement_m_sqrt_kg,
            vq,
            s_v.velocity_m_sqrt_kg_per_s,
            s_f.displacement_m_sqrt_kg,
            s_f.velocity_m_sqrt_kg_per_s,
        )
    }

    /// Nominal output sample period [s].
    #[must_use]
    pub const fn sample_period_s(&self) -> f64 {
        self.sample_period_s
    }

    /// Original immutable physical/observer budgets, also enforced by coupled
    /// owners that restore a complete externally propagated candidate state.
    #[must_use]
    pub const fn budget(&self) -> ModalAcousticTimeBudget {
        self.budget
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
            if omega == 0.0 {
                return Err(ModalAcousticTimeError::InvalidInput {
                    what: "a free coordinate has no unique isolated static equilibrium; supply its motion",
                });
            }
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
            if mode.angular_frequency_rad_s == 0.0 && *transfer != C64::ZERO {
                return Err(ModalAcousticTimeError::InvalidInput {
                    what: "free coordinates have no narrow-band pressure transfer",
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
            if omega == 0.0 {
                if *transfer != C64::ZERO {
                    return Err(ModalAcousticTimeError::InvalidInput {
                        what: "free coordinates have no narrow-band pressure transfer",
                    });
                }
                // No acoustic contribution, and no fictitious F/omega^2 baseline.
                continue;
            }
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
        let mut workspace = ModalAcousticWorkspace::new(self);
        self.step_duration_into(generalized_force_n_per_sqrt_kg, duration_s, &mut workspace)?;
        Ok(workspace.frame)
    }

    /// Advance one nominal sample using caller-owned, preallocated storage.
    ///
    /// This is the same transition and all the same gates as [`Self::step`],
    /// with no buffer construction, resize, clone or allocation during the step.
    /// The borrowed frame remains valid until the workspace is reused; callers
    /// needing an owned history may explicitly clone it outside the callback.
    ///
    /// # Errors
    /// Every refusal from [`Self::step`], plus an incompatible workspace size.
    /// On refusal the model remains unchanged; scratch diagnostics are not published.
    pub fn step_into<'w>(
        &mut self,
        generalized_force_n_per_sqrt_kg: &[f64],
        workspace: &'w mut ModalAcousticWorkspace,
    ) -> Result<&'w ModalAcousticFrame, ModalAcousticTimeError> {
        self.step_duration_into(
            generalized_force_n_per_sqrt_kg,
            self.sample_period_s,
            workspace,
        )
    }

    /// Allocation-free exact-ZOH advance for a positive caller-specified duration.
    /// Full/half samples reuse the existing cached transition coefficients;
    /// other durations retain the existing exact fallback. No new integrator or
    /// relaxed budget/passivity test is introduced.
    ///
    /// # Errors
    /// Invalid duration, force, workspace size, or any state/passivity refusal
    /// from [`Self::step`]. No accepted model state changes on failure.
    #[allow(clippy::too_many_lines)] // one shared, transactional physics loop
    pub fn step_duration_into<'w>(
        &mut self,
        generalized_force_n_per_sqrt_kg: &[f64],
        duration_s: f64,
        workspace: &'w mut ModalAcousticWorkspace,
    ) -> Result<&'w ModalAcousticFrame, ModalAcousticTimeError> {
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

        if workspace.mode_count() != self.modes.len() {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "modal workspace must match the admitted mode count",
            });
        }
        let mut pressure_pa = 0.0;
        let mut input_work_j = 0.0;
        let mut energy_before_j = 0.0;
        let mut energy_after_j = 0.0;
        let mut viscous_dissipation_j = 0.0;
        let mut dissipation_roundoff_tolerance_j = 0.0;
        let transition_coefficients = self.transition_kernel.for_duration(&self.modes, duration_s);
        for (mode_index, ((mode, state), force)) in self
            .modes
            .iter()
            .zip(&self.states)
            .zip(generalized_force_n_per_sqrt_kg)
            .enumerate()
        {
            let prepared = match transition_coefficients {
                Some(coefficients) => coefficients[mode_index],
                None => PreparedZohCoefficients::new(*mode, duration_s),
            };
            let (next, displacement_increment) =
                advance_exact_zoh_increment(*mode, *state, *force, prepared.transition);
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
            // The held force acts over the displacement increment the transition
            // itself added to q0, never over `next - state`. Rounding q0 + increment
            // onto the grid of q costs up to eps*|q|/2, and a force times that is
            // the whole balance error whenever one step moves q by a small fraction
            // of itself: the stored positions are then the wrong operands. This is
            // also the free coordinate's reason (it stores no positional energy, and
            // subtracting absolute positions can fail near velocity reversal).
            let work = force * displacement_increment;
            input_work_j += work;
            let (loss, allowance) = prepared
                .viscous_integral
                .evaluate(*mode, *state, *force, duration_s);
            let balance_loss = work - (after - before);
            let balance_allowance =
                256.0 * f64::EPSILON * work.abs().max(before).max(after).max(f64::MIN_POSITIVE)
                    + allowance;
            if !loss.is_finite() || !allowance.is_finite() {
                return Err(ModalAcousticTimeError::InvalidInput {
                    what: "mode-local viscous integral must be finite",
                });
            }
            if balance_loss < -balance_allowance {
                return Err(ModalAcousticTimeError::NegativeDissipation {
                    dissipation_j: balance_loss,
                    tolerance_j: balance_allowance,
                });
            }
            if (balance_loss - loss).abs() > balance_allowance {
                return Err(ModalAcousticTimeError::InvalidInput {
                    what: "independent viscous heat and held-force energy balance disagree",
                });
            }
            viscous_dissipation_j += loss;
            dissipation_roundoff_tolerance_j += allowance;
            workspace.frame.modal_viscous_losses[mode_index] = ModalViscousLoss {
                energy_j: loss,
                roundoff_tolerance_j: allowance,
            };
            pressure_pa += mode.pressure_per_modal_velocity.re * next.velocity_m_sqrt_kg_per_s
                + mode.pressure_per_modal_velocity.im
                    * mode.angular_frequency_rad_s
                    * next.displacement_m_sqrt_kg;
            workspace.frame.modal_energy_j[mode_index] = after;
            workspace.candidate[mode_index] = next;
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
        let balance_loss = input_work_j - (energy_after_j - energy_before_j);
        let scale = input_work_j
            .abs()
            .max(energy_before_j)
            .max(energy_after_j)
            .max(f64::MIN_POSITIVE);
        dissipation_roundoff_tolerance_j +=
            (16.0 + self.modes.len() as f64) * f64::EPSILON * viscous_dissipation_j;
        let balance_allowance = 256.0 * f64::EPSILON * scale + dissipation_roundoff_tolerance_j;
        if !viscous_dissipation_j.is_finite()
            || !dissipation_roundoff_tolerance_j.is_finite()
            || (balance_loss - viscous_dissipation_j).abs() > balance_allowance
        {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "aggregate viscous integral and held-force energy balance disagree",
            });
        }
        self.states.copy_from_slice(&workspace.candidate);
        workspace.frame.observer_pressure_pa = pressure_pa;
        workspace.frame.total_modal_energy_j = energy_after_j;
        workspace.frame.input_work_j = input_work_j;
        workspace.frame.viscous_dissipation_j = viscous_dissipation_j;
        workspace.frame.dissipation_roundoff_tolerance_j = dissipation_roundoff_tolerance_j;
        Ok(&workspace.frame)
    }

    /// Solve the exact-ZOH held generalized force that drives one physical
    /// port to `target_port_velocity_m_s` after `duration_s`, without
    /// advancing state.
    ///
    /// `port_shapes` are the same conjugate mode-shape weights used both to
    /// project a physical port force into each mode's generalized force
    /// (`generalized_force_k = port_shapes[k] * port_force`) and to recover
    /// physical port velocity from modal velocities
    /// (`port_velocity = sum_k port_shapes[k] * velocity_k`). The
    /// zero-order-held transition is affine in the applied force, so this
    /// evaluates the same transition used by [`Self::step_duration`] at
    /// force `0` (free modal motion) and force `1` (unit-port-force
    /// mobility) per mode, then inverts the resulting scalar affine equation
    /// for the port force. Use that force as `port_shapes[k] * force` when
    /// actually stepping, to reach the same endpoint this solved for.
    ///
    /// # Errors
    /// Refuses a wrong-length shape vector, a nonpositive/nonfinite
    /// duration, nonfinite shapes or target velocity, or a nonpositive or
    /// nonfinite port mobility (for example an undamped port sampled at a
    /// duration where a held force does not change the endpoint velocity at
    /// all). Never mutates state.
    pub fn held_force_for_port_velocity(
        &self,
        port_shapes: &[f64],
        target_port_velocity_m_s: f64,
        duration_s: f64,
    ) -> Result<f64, ModalAcousticTimeError> {
        if !(duration_s > 0.0 && duration_s.is_finite()) {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "held-force duration must be positive and finite",
            });
        }
        if port_shapes.len() != self.modes.len() {
            return Err(ModalAcousticTimeError::ForceCountMismatch {
                expected: self.modes.len(),
                found: port_shapes.len(),
            });
        }
        if !target_port_velocity_m_s.is_finite() || port_shapes.iter().any(|phi| !phi.is_finite()) {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "held-force port shapes and target velocity must be finite",
            });
        }
        let transition_coefficients = self.transition_kernel.for_duration(&self.modes, duration_s);
        let mut free_port_velocity_m_s = 0.0;
        let mut port_mobility_m_s_per_n = 0.0;
        for (mode_index, (mode, state)) in self.modes.iter().zip(&self.states).enumerate() {
            let phi = port_shapes[mode_index];
            let (free, unit) = match transition_coefficients {
                Some(coefficients) => (
                    advance_exact_zoh_cached(
                        *mode,
                        *state,
                        0.0,
                        coefficients[mode_index].transition,
                    ),
                    advance_exact_zoh_cached(
                        *mode,
                        ModalAcousticState::default(),
                        1.0,
                        coefficients[mode_index].transition,
                    ),
                ),
                None => (
                    advance_exact_zoh(*mode, *state, 0.0, duration_s),
                    advance_exact_zoh(*mode, ModalAcousticState::default(), 1.0, duration_s),
                ),
            };
            // Mobility is the force response from rest, independent of the
            // accepted state. Subtracting two large free velocities can erase it.
            let mode_mobility = unit.velocity_m_sqrt_kg_per_s;
            free_port_velocity_m_s += phi * free.velocity_m_sqrt_kg_per_s;
            port_mobility_m_s_per_n += phi * phi * mode_mobility;
        }
        if !free_port_velocity_m_s.is_finite() {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "free port velocity evaluation produced a non-finite result",
            });
        }
        if !(port_mobility_m_s_per_n > 0.0 && port_mobility_m_s_per_n.is_finite()) {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "port mobility must be positive and finite",
            });
        }
        let held_force_n =
            (target_port_velocity_m_s - free_port_velocity_m_s) / port_mobility_m_s_per_n;
        if !held_force_n.is_finite() {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "held-force solve produced a non-finite result",
            });
        }
        Ok(held_force_n)
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
    f64::midpoint(
        state.velocity_m_sqrt_kg_per_s * state.velocity_m_sqrt_kg_per_s,
        mode.angular_frequency_rad_s
            * mode.angular_frequency_rad_s
            * state.displacement_m_sqrt_kg
            * state.displacement_m_sqrt_kg,
    )
}

fn exact_zoh_coefficients(mode: ModalAcousticMode, dt: f64) -> ExactZohCoefficients {
    let omega = mode.angular_frequency_rad_s;
    if omega == 0.0 {
        return ExactZohCoefficients::Free { duration_s: dt };
    }
    let zeta = mode.damping_ratio;
    if zeta == 0.0 {
        let angle = omega * dt;
        let half_sine = det::sin(0.5 * angle);
        return ExactZohCoefficients::Undamped {
            sine: det::sin(angle),
            one_minus_cosine: 2.0 * half_sine * half_sine,
        };
    }
    // For q'' + c q' + omega^2 q = F, B is the velocity response to a
    // held unit force and G is its displacement response. E11 = 1-omega^2 G,
    // E12 = B, E21 = -omega^2 B, E22 = E11-c B. Compute the convolution
    // itself: subtracting and restoring F/omega^2 loses the entire attack
    // when the equilibrium is large compared with one sample's motion.
    //
    // A small-interval entire series followed by transition doubling covers
    // under-, critical- and overdamping without root subtraction or a band
    // snapped to critical damping. Only preparation pays this cost. B/h and
    // G/h^2 stay scaled during doubling, including strongly overdamped modes.
    let mut angle = omega * dt;
    let mut drag = 2.0 * (zeta * angle);
    if !angle.is_finite() || !drag.is_finite() {
        return ExactZohCoefficients::Damped {
            position_retention: f64::NAN,
            velocity_retention: f64::NAN,
            velocity_gain: f64::NAN,
            displacement_gain: f64::NAN,
            position_decrement: f64::NAN,
            velocity_decrement: f64::NAN,
        };
    }
    let mut doublings = 0;
    while angle.abs().max(drag.abs()) > 0.25 {
        angle *= 0.5;
        drag *= 0.5;
        doublings += 1;
    }
    let mut transition = DampedTransition::small(angle, drag);
    for _ in 0..doublings {
        transition.double();
    }
    transition.coefficients(dt)
}

// The exact integral of velocity for a held force on a free coordinate.
// Shared by state evolution and work evaluation; no absolute-position subtraction.
fn free_displacement_increment(velocity: f64, force: f64, duration_s: f64) -> f64 {
    duration_s * (velocity + 0.5 * (duration_s * force))
}

fn advance_exact_zoh_cached(
    mode: ModalAcousticMode,
    state: ModalAcousticState,
    force: f64,
    coefficients: ExactZohCoefficients,
) -> ModalAcousticState {
    advance_exact_zoh_increment(mode, state, force, coefficients).0
}

// The same transition, bit for bit in q and v, that also returns the displacement
// increment it added to q0. Free, undamped and damped short steps form q = q0 + increment,
// so the increment is that exact addend and keeps the digits the sum rounds away; `q - q0`
// cannot recover them. The damped long-step branch is a sum of three comparable terms with
// no separate addend, so it returns `q - q0`.
fn advance_exact_zoh_increment(
    mode: ModalAcousticMode,
    state: ModalAcousticState,
    force: f64,
    coefficients: ExactZohCoefficients,
) -> (ModalAcousticState, f64) {
    let omega = mode.angular_frequency_rad_s;
    let q0 = state.displacement_m_sqrt_kg;
    let v0 = state.velocity_m_sqrt_kg_per_s;
    let (q, v, displacement_increment) = match coefficients {
        ExactZohCoefficients::Free { duration_s } => {
            // Exact polynomial solution of q''=force under a held force.
            // Full, half and arbitrary durations use this same arithmetic.
            let delta_v = duration_s * force;
            let displacement_increment = free_displacement_increment(v0, force, duration_s);
            return (
                ModalAcousticState {
                    displacement_m_sqrt_kg: state.displacement_m_sqrt_kg + displacement_increment,
                    velocity_m_sqrt_kg_per_s: v0 + delta_v,
                },
                displacement_increment,
            );
        }
        ExactZohCoefficients::Undamped {
            sine,
            one_minus_cosine,
        } => {
            return advance_undamped_zoh(mode, state, force, sine, one_minus_cosine);
        }
        ExactZohCoefficients::Damped {
            position_retention,
            velocity_retention,
            velocity_gain,
            displacement_gain,
            position_decrement,
            velocity_decrement,
        } => {
            let restoring_force = force - omega * omega * q0;
            let (q, displacement_increment) = if position_decrement.abs() < 0.5 {
                let increment = velocity_gain * v0 + displacement_gain * restoring_force;
                (q0 + increment, increment)
            } else {
                let q = position_retention * q0 + velocity_gain * v0 + displacement_gain * force;
                (q, q - q0)
            };
            let v = if velocity_decrement.abs() < 0.5 {
                v0 + (velocity_gain * restoring_force - velocity_decrement * v0)
            } else {
                velocity_retention * v0 + velocity_gain * restoring_force
            };
            (q, v, displacement_increment)
        }
    };
    (
        ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v,
        },
        displacement_increment,
    )
}

// Advance about the current displacement rather than subtracting and adding
// the possibly much larger force equilibrium. In a short step, 1-cos(theta)
// is evaluated by its half-angle identity, so force * delta_q retains the
// precision needed by the unchanged work/energy check.
fn advance_undamped_zoh(
    mode: ModalAcousticMode,
    state: ModalAcousticState,
    force: f64,
    sine: f64,
    one_minus_cosine: f64,
) -> (ModalAcousticState, f64) {
    let omega = mode.angular_frequency_rad_s;
    let q = state.displacement_m_sqrt_kg;
    let v = state.velocity_m_sqrt_kg_per_s;
    let displacement_increment =
        (v / omega) * sine + (force / (omega * omega) - q) * one_minus_cosine;
    let velocity_increment = (force / omega - omega * q) * sine - v * one_minus_cosine;
    (
        ModalAcousticState {
            displacement_m_sqrt_kg: q + displacement_increment,
            velocity_m_sqrt_kg_per_s: v + velocity_increment,
        },
        displacement_increment,
    )
}

/// One exact zero-order-hold step of a single mode from `state` under the held generalized `force` over `dt`.
pub fn advance_exact_zoh(
    mode: ModalAcousticMode,
    state: ModalAcousticState,
    force: f64,
    dt: f64,
) -> ModalAcousticState {
    advance_exact_zoh_cached(mode, state, force, exact_zoh_coefficients(mode, dt))
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

    fn assert_state_bits_eq(actual: ModalAcousticState, expected: ModalAcousticState) {
        assert_eq!(
            actual.displacement_m_sqrt_kg.to_bits(),
            expected.displacement_m_sqrt_kg.to_bits()
        );
        assert_eq!(
            actual.velocity_m_sqrt_kg_per_s.to_bits(),
            expected.velocity_m_sqrt_kg_per_s.to_bits()
        );
    }

    #[test]
    fn g1_damped_held_force_matches_high_precision_exponential_across_scales() {
        // Independent 90-decimal-digit exponential of the augmented matrix
        // [[0,1,0],[-omega²,-2*zeta*omega,F],[0,0,0]], acting on [q,v,1].
        // Inputs were lifted from their exact f64 values before exponentiation.
        // Cover soft modes, tiny physical loss, both sides of critical damping,
        // and a viscous mode whose slow pole used to cancel to zero.
        let cases = [
            (
                1.0e-6,
                0.05,
                100.0,
                1.0 / 48_000.0,
                0.0,
                0.0,
                2.170_138_888_887_381_5e-8,
                0.002_083_333_333_331_163,
            ),
            (
                1000.0,
                1.0e-12,
                24_200.0,
                1.0 / 48_000.0,
                0.0002,
                -0.3,
                0.000_198_958_597_058_404_05,
                0.200_028_933_615_378_96,
            ),
            (
                1000.0,
                0.999_999_995,
                3.0,
                0.002,
                0.0002,
                -0.3,
                1.781_981_824_719_25e-6,
                -0.012_721_516_531_311_369,
            ),
            (
                1000.0,
                1.000_000_005,
                3.0,
                0.002,
                0.0002,
                -0.3,
                1.781_983_077_021_73e-6,
                -0.012_721_516_717_171_822,
            ),
            (
                100.0,
                1.0e12,
                2.0e14,
                0.1,
                0.0,
                0.0,
                0.099_999_999_999_745,
                0.999_999_999_995,
            ),
        ];
        for (omega, zeta, force, dt, q, v, expected_q, expected_v) in cases {
            let mode = ModalAcousticMode {
                angular_frequency_rad_s: omega,
                damping_ratio: zeta,
                pressure_per_modal_velocity: C64::ZERO,
            };
            let initial = ModalAcousticState {
                displacement_m_sqrt_kg: q,
                velocity_m_sqrt_kg_per_s: v,
            };
            let actual = advance_exact_zoh(mode, initial, force, dt);
            assert!(
                (actual.displacement_m_sqrt_kg - expected_q).abs() <= 2.0e-13 * expected_q.abs(),
                "omega={omega}, zeta={zeta}: q={} expected={expected_q}",
                actual.displacement_m_sqrt_kg
            );
            assert!(
                (actual.velocity_m_sqrt_kg_per_s - expected_v).abs() <= 2.0e-13 * expected_v.abs(),
                "omega={omega}, zeta={zeta}: v={} expected={expected_v}",
                actual.velocity_m_sqrt_kg_per_s
            );
            let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
            model.restore_states(&[initial]).unwrap();
            let frame = model.step_duration(&[force], dt).unwrap();
            assert_state_bits_eq(model.states()[0], actual);
            assert!(frame.viscous_dissipation_j >= -frame.dissipation_roundoff_tolerance_j);
        }
    }

    #[test]
    fn g1_port_mobility_does_not_subtract_large_free_velocities() {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: 1000.0,
            damping_ratio: 0.05,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let mut limits = budget();
        limits.maximum_abs_velocity_m_sqrt_kg_per_s = 1.0e15;
        limits.maximum_total_energy_j = 1.0e30;
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], limits).unwrap();
        let state = ModalAcousticState {
            displacement_m_sqrt_kg: 0.0,
            velocity_m_sqrt_kg_per_s: 1.0e14,
        };
        model.restore_states(&[state]).unwrap();
        let dt = model.sample_period_s();
        let free = advance_exact_zoh(mode, state, 0.0, dt).velocity_m_sqrt_kg_per_s;
        let response = advance_exact_zoh(mode, ModalAcousticState::default(), 1.0, dt)
            .velocity_m_sqrt_kg_per_s;
        let target = free + 1.0e12 * response;
        let force = model
            .held_force_for_port_velocity(&[1.0], target, dt)
            .unwrap();
        let actual = advance_exact_zoh(mode, state, force, dt).velocity_m_sqrt_kg_per_s;
        assert!((actual - target).abs() <= 2.0 * f64::EPSILON * target.abs());
        assert_eq!(
            model.states(),
            &[state],
            "solving the port must not advance its state"
        );
    }

    #[test]
    fn g1_forced_undamped_steps_preserve_work_without_equilibrium_cancellation() {
        let state = ModalAcousticState {
            displacement_m_sqrt_kg: 0.0002,
            velocity_m_sqrt_kg_per_s: -0.3,
        };
        // Independent Taylor solution of q'' = F - omega² q. Each parity
        // advances by two derivative orders; no equilibrium subtraction or
        // trigonometric implementation is shared with production.
        let reference = |omega: f64, force: f64, dt: f64| {
            let acceleration = force - omega * omega * state.displacement_m_sqrt_kg;
            let theta_squared = (omega * dt).powi(2);
            let mut q_odd = state.velocity_m_sqrt_kg_per_s * dt;
            let mut q_even = 0.5 * acceleration * dt * dt;
            let mut v_odd = acceleration * dt;
            let mut v_even = -0.5 * omega * omega * state.velocity_m_sqrt_kg_per_s * dt * dt;
            let mut q = state.displacement_m_sqrt_kg;
            let mut v = state.velocity_m_sqrt_kg_per_s;
            for pair in 0..6 {
                q += q_odd + q_even;
                v += v_odd + v_even;
                let n = f64::from(2 * pair + 1);
                q_odd *= -theta_squared / ((n + 1.0) * (n + 2.0));
                v_odd *= -theta_squared / ((n + 1.0) * (n + 2.0));
                q_even *= -theta_squared / ((n + 2.0) * (n + 3.0));
                v_even *= -theta_squared / ((n + 2.0) * (n + 3.0));
            }
            (q, v)
        };
        for omega in [800.0, 1000.0, 3200.0] {
            for force in [-24_200.0, 24_200.0] {
                for fraction in [1.0, 0.5, 0.0625] {
                    let mode = ModalAcousticMode {
                        angular_frequency_rad_s: omega,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO,
                    };
                    let mut model =
                        ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
                    model.restore_states(&[state]).unwrap();
                    let dt = fraction * model.sample_period_s();
                    let frame = model.step_duration(&[force], dt).unwrap();
                    let (q, v) = reference(omega, force, dt);
                    let actual = model.states()[0];
                    assert!((actual.displacement_m_sqrt_kg - q).abs() < 1.0e-18);
                    assert!((actual.velocity_m_sqrt_kg_per_s - v).abs() < 2.0e-15);
                    assert!(
                        frame.viscous_dissipation_j.abs() <= frame.dissipation_roundoff_tolerance_j
                    );
                }
            }
        }
        // Reproduce the old formula's false energy rejection at the exact
        // parameters exposed by the attempted contact constraint.
        let omega = 1000.0;
        let force = 24_200.0;
        let dt = 1.0 / 48_000.0;
        let equilibrium = force / (omega * omega);
        let centered_q = state.displacement_m_sqrt_kg - equilibrium;
        let (sine, cosine) = (det::sin(omega * dt), det::cos(omega * dt));
        let old = ModalAcousticState {
            displacement_m_sqrt_kg: cosine * centered_q
                + sine / omega * state.velocity_m_sqrt_kg_per_s
                + equilibrium,
            velocity_m_sqrt_kg_per_s: -omega * sine * centered_q
                + cosine * state.velocity_m_sqrt_kg_per_s,
        };
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let before = modal_energy(mode, state);
        let after = modal_energy(mode, old);
        let work = force * (old.displacement_m_sqrt_kg - state.displacement_m_sqrt_kg);
        let allowance = 256.0 * f64::EPSILON * work.abs().max(before).max(after);
        assert!(work - (after - before) < -allowance);
    }

    #[test]
    fn g1_damped_held_force_retains_motion_far_below_force_equilibrium() {
        // F/omega^2 = 3e6, but the first displacement is only 1.5e-12.
        // Subtracting/restoring equilibrium erased that physical attack.
        let dt = 1.0e-6;
        let force = 3.0;
        for damping_ratio in [1.0e-9, 0.2, 1.0, 2.0] {
            let mode = ModalAcousticMode {
                angular_frequency_rad_s: 1.0e-3,
                damping_ratio,
                pressure_per_modal_velocity: C64::ZERO,
            };
            let next = advance_exact_zoh(mode, ModalAcousticState::default(), force, dt);
            // Independent short-time expansion through the first drag term;
            // the omitted relative term is O((c dt)^2 + (omega dt)^2).
            let drag_dt = 2.0 * damping_ratio * mode.angular_frequency_rad_s * dt;
            let expected_q = 0.5 * force * dt * dt * (1.0 - drag_dt / 3.0);
            let expected_v = force * dt * (1.0 - 0.5 * drag_dt);
            assert!(
                (next.displacement_m_sqrt_kg / expected_q - 1.0).abs() < 2.0e-14,
                "zeta={damping_ratio} q={} expected={expected_q}",
                next.displacement_m_sqrt_kg
            );
            assert!((next.velocity_m_sqrt_kg_per_s / expected_v - 1.0).abs() < 2.0e-14);
        }
    }

    #[test]
    fn g1_near_critical_modes_keep_their_actual_decay_law() {
        let omega = 7.0;
        let dt = 0.8 / omega;
        let initial = ModalAcousticState {
            displacement_m_sqrt_kg: 1.0,
            velocity_m_sqrt_kg_per_s: 0.0,
        };
        let mut positions = Vec::new();
        for zeta in [1.0_f64 - 5.0e-9, 1.0, 1.0 + 5.0e-9] {
            let mode = ModalAcousticMode {
                angular_frequency_rad_s: omega,
                damping_ratio: zeta,
                pressure_per_modal_velocity: C64::ZERO,
            };
            let delta = zeta * omega * dt;
            let angle = omega * dt * ((1.0 - zeta) * (1.0 + zeta)).abs().sqrt();
            let (cosine, sinc) = if zeta < 1.0 {
                (angle.cos(), angle.sin() / angle)
            } else if zeta > 1.0 {
                (angle.cosh(), angle.sinh() / angle)
            } else {
                (1.0, 1.0)
            };
            let expected_q = (-delta).exp() * (cosine + delta * sinc);
            let expected_v = -omega * omega * dt * (-delta).exp() * sinc;
            let next = advance_exact_zoh(mode, initial, 0.0, dt);
            assert!((next.displacement_m_sqrt_kg - expected_q).abs() < 4.0e-14);
            assert!((next.velocity_m_sqrt_kg_per_s - expected_v).abs() < 4.0e-14);
            positions.push(next.displacement_m_sqrt_kg);
        }
        // These physically different damping ratios previously all used zeta=1.
        assert!(positions[0] < positions[1] && positions[1] < positions[2]);
    }

    #[test]
    fn g1_strongly_overdamped_mode_retains_its_slow_relaxation() {
        let omega = 2.0;
        let zeta = 1.0e10;
        // Stable product-of-roots expression; zeta - sqrt(zeta^2-1) rounds
        // to zero here even though the slow decay is physically nonzero.
        let root = ((zeta - 1.0_f64) * (zeta + 1.0)).sqrt();
        let slow = -omega / (zeta + root);
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: zeta,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let initial = ModalAcousticState {
            displacement_m_sqrt_kg: 0.75,
            velocity_m_sqrt_kg_per_s: 0.75 * slow,
        };
        let dt = -1.0 / slow;
        let expected_q = 0.75 * (-1.0_f64).exp();
        let expected_v = slow * expected_q;
        let next = advance_exact_zoh(mode, initial, 0.0, dt);
        assert!((next.displacement_m_sqrt_kg / expected_q - 1.0).abs() < 2.0e-12);
        assert!((next.velocity_m_sqrt_kg_per_s / expected_v - 1.0).abs() < 2.0e-12);
        let half = advance_exact_zoh(mode, initial, 0.0, 0.5 * dt);
        let split = advance_exact_zoh(mode, half, 0.0, 0.5 * dt);
        assert!((split.displacement_m_sqrt_kg / expected_q - 1.0).abs() < 2.0e-12);
        assert!((split.velocity_m_sqrt_kg_per_s / expected_v - 1.0).abs() < 2.0e-12);
    }

    #[test]
    fn g1_critical_ring_down_preserves_a_tail_smaller_than_machine_epsilon() {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: 1.0,
            damping_ratio: 1.0,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let initial = ModalAcousticState {
            displacement_m_sqrt_kg: 1.0,
            velocity_m_sqrt_kg_per_s: 0.0,
        };
        let next = advance_exact_zoh(mode, initial, 0.0, 40.0);
        assert!((next.displacement_m_sqrt_kg / (41.0 * (-40.0_f64).exp()) - 1.0).abs() < 2.0e-12);
        assert!(
            (next.velocity_m_sqrt_kg_per_s / (-40.0 * (-40.0_f64).exp()) - 1.0).abs() < 2.0e-12
        );
    }

    #[test]
    fn g1_cached_high_q_ring_down_preserves_physical_phase_and_decay() {
        let omega = core::f64::consts::TAU * 440.0;
        let zeta = 1.0e-4;
        let dt = 1.0 / 48_000.0;
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: zeta,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let coefficients = exact_zoh_coefficients(mode, dt);
        let mut state = ModalAcousticState {
            displacement_m_sqrt_kg: 1.0,
            velocity_m_sqrt_kg_per_s: 0.0,
        };
        for _ in 0..48_000 {
            state = advance_exact_zoh_cached(mode, state, 0.0, coefficients);
        }
        // One second / 440 cycles: independently evaluate the continuous
        // free-decay solution, including the damping-induced frequency shift.
        let root = (1.0 - zeta * zeta).sqrt();
        let angle = omega * root;
        let decay = (-zeta * omega).exp();
        let expected_q = decay * (angle.cos() + zeta / root * angle.sin());
        let expected_v = -decay * omega / root * angle.sin();
        assert!((state.displacement_m_sqrt_kg - expected_q).abs() < 2.0e-10);
        assert!((state.velocity_m_sqrt_kg_per_s - expected_v).abs() / omega < 2.0e-10);
    }

    #[test]
    fn g0_cached_full_and_half_transitions_are_bit_exact_in_every_damping_regime() {
        let sample_period_s = f64::from(48_000_u32).recip();
        let state = ModalAcousticState {
            displacement_m_sqrt_kg: -2.75e-4,
            velocity_m_sqrt_kg_per_s: 0.1875,
        };
        for damping_ratio in [
            0.0,
            0.03,
            1.0 - 2.0e-8,
            1.0 - 0.5e-8,
            1.0,
            1.0 + 0.5e-8,
            1.0 + 2.0e-8,
            2.0,
        ] {
            let mode = ModalAcousticMode {
                angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 1_375.0,
                damping_ratio,
                pressure_per_modal_velocity: C64::new(1.25, -0.375),
            };
            for dt in [sample_period_s, 0.5 * sample_period_s] {
                let expected = advance_exact_zoh(mode, state, 3.125, dt);
                let actual =
                    advance_exact_zoh_cached(mode, state, 3.125, exact_zoh_coefficients(mode, dt));
                assert_state_bits_eq(actual, expected);
            }
        }

        let model = ModalAcousticTimeModel::try_new(
            48_000,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 1_375.0,
                damping_ratio: 0.03,
                pressure_per_modal_velocity: C64::ZERO,
            }],
            budget(),
        )
        .unwrap();
        let clone = model.clone();
        assert!(Arc::ptr_eq(
            &model.transition_kernel,
            &clone.transition_kernel
        ));
        assert!(
            model
                .transition_kernel
                .for_duration(&model.modes, sample_period_s)
                .is_some()
        );
        assert!(
            model
                .transition_kernel
                .for_duration(&model.modes, 0.5 * sample_period_s)
                .is_some()
        );
    }

    #[test]
    fn g0_arbitrary_duration_falls_back_exactly_and_refusals_remain_transactional() {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 925.0,
            damping_ratio: 0.08,
            pressure_per_modal_velocity: C64::new(0.75, -0.125),
        };
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
        let state = ModalAcousticState {
            displacement_m_sqrt_kg: 1.25e-4,
            velocity_m_sqrt_kg_per_s: -0.0625,
        };
        model.restore_states(&[state]).unwrap();
        let dt = 0.37 / 48_000.0;
        assert!(
            model
                .transition_kernel
                .for_duration(&model.modes, dt)
                .is_none()
        );
        let expected = advance_exact_zoh(mode, state, -1.75, dt);
        model.step_duration(&[-1.75], dt).unwrap();
        assert_state_bits_eq(model.states()[0], expected);

        let before_refusals = model.states().to_vec();
        assert_eq!(
            model.step_duration(&[], 0.0),
            Err(ModalAcousticTimeError::InvalidInput {
                what: "step duration must be positive and finite",
            })
        );
        assert!(matches!(
            model.step_duration(&[f64::NAN], dt),
            Err(ModalAcousticTimeError::InvalidInput {
                what: "generalized forces must be finite"
            })
        ));
        assert_eq!(model.states(), before_refusals);
    }

    #[test]
    #[ignore = "temporary Pass 18 exact-ZOH performance A/B"]
    fn pass18_optimization_micro_ab_exact_zoh() {
        use std::hint::black_box;
        use std::time::Instant;

        const MODES: usize = 32;
        const STEPS: usize = 32_768;
        const ROUNDS: usize = 6;

        let dt = f64::from(6_144_000_u32).recip();
        let modes = (0..MODES)
            .map(|index| ModalAcousticMode {
                angular_frequency_rad_s: 2.0
                    * core::f64::consts::PI
                    * (646.0 + 113.0 * index as f64),
                damping_ratio: 0.000_2 + 0.000_07 * index as f64,
                pressure_per_modal_velocity: C64::ZERO,
            })
            .collect::<Vec<_>>();
        let coefficients = modes
            .iter()
            .copied()
            .map(|mode| exact_zoh_coefficients(mode, dt))
            .collect::<Vec<_>>();

        let run = |cached: bool| {
            let mut states = (0..MODES)
                .map(|index| ModalAcousticState {
                    displacement_m_sqrt_kg: 1.0e-8 * (index as f64 + 1.0),
                    velocity_m_sqrt_kg_per_s: -2.0e-5 * (index as f64 + 1.0),
                })
                .collect::<Vec<_>>();
            let start = Instant::now();
            for step in 0..STEPS {
                for index in 0..MODES {
                    let mode = black_box(modes[index]);
                    let state = black_box(states[index]);
                    let force =
                        black_box((index as f64 + 1.0) * 2.0e-4 + (step & 31) as f64 * 1.0e-7);
                    states[index] = if cached {
                        advance_exact_zoh_cached(mode, state, force, coefficients[index])
                    } else {
                        advance_exact_zoh(mode, state, force, dt)
                    };
                }
            }
            let elapsed = start.elapsed();
            let digest = states.iter().fold(0_u64, |bits, state| {
                bits ^ state.displacement_m_sqrt_kg.to_bits().rotate_left(7)
                    ^ state.velocity_m_sqrt_kg_per_s.to_bits()
            });
            black_box((elapsed, digest))
        };

        let mut uncached_total = 0.0;
        let mut cached_total = 0.0;
        for round in 0..ROUNDS {
            let (uncached, uncached_digest, cached, cached_digest) = if round % 2 == 0 {
                let (uncached, uncached_digest) = run(false);
                let (cached, cached_digest) = run(true);
                (uncached, uncached_digest, cached, cached_digest)
            } else {
                let (cached, cached_digest) = run(true);
                let (uncached, uncached_digest) = run(false);
                (uncached, uncached_digest, cached, cached_digest)
            };
            assert_eq!(cached_digest, uncached_digest);
            let uncached_s = uncached.as_secs_f64();
            let cached_s = cached.as_secs_f64();
            uncached_total += uncached_s;
            cached_total += cached_s;
            eprintln!(
                "PASS18_AB round={round} uncached_s={uncached_s:.6} cached_s={cached_s:.6} speedup={:.4}",
                uncached_s / cached_s
            );
        }
        eprintln!(
            "PASS18_SUMMARY uncached_s={uncached_total:.6} cached_s={cached_total:.6} speedup={:.4}",
            uncached_total / cached_total
        );
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
    fn g1_mode_local_loss_survives_a_large_undamped_energy_background() {
        let large = ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 400.0,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let small = ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 1000.0,
            damping_ratio: 0.2,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let small_state = ModalAcousticState {
            displacement_m_sqrt_kg: 1e-9,
            velocity_m_sqrt_kg_per_s: 0.0,
        };
        let mut combined =
            ModalAcousticTimeModel::try_new(48_000, vec![large, small], budget()).unwrap();
        combined
            .restore_states(&[
                ModalAcousticState {
                    displacement_m_sqrt_kg: 10.0,
                    velocity_m_sqrt_kg_per_s: 0.0,
                },
                small_state,
            ])
            .unwrap();
        let mut isolated = ModalAcousticTimeModel::try_new(48_000, vec![small], budget()).unwrap();
        isolated.restore_states(&[small_state]).unwrap();
        let combined = combined.step(&[0.0, 0.0]).unwrap();
        let isolated = isolated.step(&[0.0]).unwrap();
        let loss = combined.modal_viscous_losses[1];
        assert_eq!(loss, isolated.modal_viscous_losses[0]);
        assert!(loss.energy_j > 100.0 * loss.roundoff_tolerance_j);
        assert!(loss.energy_j < f64::EPSILON * combined.total_modal_energy_j);
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
    fn g0_mode_cap_4096_admitted_4097_refused() {
        // 192 kHz keeps every mode below the Nyquist guard (0.9 * 96 kHz), so only the mode-count cap is exercised.
        let ladder = |count: usize| -> Vec<ModalAcousticMode> {
            (0..count)
                .map(|k| ModalAcousticMode {
                    angular_frequency_rad_s: 2.0
                        * core::f64::consts::PI
                        * (100.0 + 10.0 * k as f64),
                    damping_ratio: 0.01,
                    pressure_per_modal_velocity: C64::from_re(1.0),
                })
                .collect()
        };
        assert_eq!(MAX_TIME_DOMAIN_ACOUSTIC_MODES, 4_096);
        let admitted = ModalAcousticTimeModel::try_new(192_000, ladder(4_096), budget())
            .expect("exactly the cap is admitted");
        assert_eq!(admitted.modes().len(), 4_096);
        let refused = ModalAcousticTimeError::InvalidInput {
            what: "mode count must be in 1..=4096",
        };
        assert_eq!(
            ModalAcousticTimeModel::try_new(192_000, ladder(4_097), budget()).err(),
            Some(refused.clone())
        );
        assert_eq!(
            ModalAcousticTimeModel::try_new(192_000, Vec::new(), budget()).err(),
            Some(refused)
        );
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

    #[test]
    fn g3_held_force_for_port_velocity_reaches_target_in_every_damping_regime() {
        // Mixed-sign multimode port (CONTRACT.md's own G3 scope), every
        // canonical damping regime, and full/half/subsample durations.
        for damping_ratio in [0.0_f64, 1.0e-9, 0.12, 1.0, 1.7] {
            let modes = vec![
                ModalAcousticMode {
                    angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 220.0,
                    damping_ratio,
                    pressure_per_modal_velocity: C64::ZERO,
                },
                ModalAcousticMode {
                    angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 880.0,
                    damping_ratio,
                    pressure_per_modal_velocity: C64::ZERO,
                },
            ];
            let mut model = ModalAcousticTimeModel::try_new(48_000, modes, budget()).unwrap();
            model
                .restore_states(&[
                    ModalAcousticState {
                        displacement_m_sqrt_kg: 3.0e-4,
                        velocity_m_sqrt_kg_per_s: -0.2,
                    },
                    ModalAcousticState {
                        displacement_m_sqrt_kg: -1.5e-4,
                        velocity_m_sqrt_kg_per_s: 0.05,
                    },
                ])
                .unwrap();
            let shapes = [0.8_f64, -0.6];
            let before = model.states().to_vec();
            let dt = model.sample_period_s();
            for duration_s in [dt, 0.5 * dt, 0.37 * dt] {
                let target_port_velocity_m_s = 0.75;
                let force = model
                    .held_force_for_port_velocity(&shapes, target_port_velocity_m_s, duration_s)
                    .unwrap();
                let mut probe = model.clone();
                let generalized_forces: Vec<f64> = shapes.iter().map(|phi| phi * force).collect();
                probe
                    .step_duration(&generalized_forces, duration_s)
                    .unwrap();
                let achieved_port_velocity_m_s: f64 = probe
                    .states()
                    .iter()
                    .zip(&shapes)
                    .map(|(s, phi)| phi * s.velocity_m_sqrt_kg_per_s)
                    .sum();
                assert!(
                    (achieved_port_velocity_m_s - target_port_velocity_m_s).abs() < 1.0e-9,
                    "damping_ratio={damping_ratio} duration_s={duration_s} achieved={achieved_port_velocity_m_s}"
                );
            }
            assert_eq!(model.states(), before, "solve must never mutate state");
        }
    }

    #[test]
    fn g1_held_force_for_port_velocity_beats_the_old_acceleration_cancellation_heuristic() {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 300.0,
            damping_ratio: 0.05,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
        model
            .restore_states(&[ModalAcousticState {
                displacement_m_sqrt_kg: 2.0e-3,
                velocity_m_sqrt_kg_per_s: -0.4,
            }])
            .unwrap();
        let shapes = [1.0_f64];
        // A deliberately non-small sub-step (omega * duration_s = 1 rad),
        // where the old small-angle "acceleration cancellation" estimate is
        // no longer a good approximation of the exact-ZOH endpoint
        // constraint this method solves exactly.
        let duration_s = 1.0 / mode.angular_frequency_rad_s;
        let target_port_velocity_m_s = 0.6;

        let exact_force = model
            .held_force_for_port_velocity(&shapes, target_port_velocity_m_s, duration_s)
            .unwrap();
        let mut exact_probe = model.clone();
        exact_probe
            .step_duration(&[shapes[0] * exact_force], duration_s)
            .unwrap();
        let exact_velocity = exact_probe.states()[0].velocity_m_sqrt_kg_per_s;
        assert!((exact_velocity - target_port_velocity_m_s).abs() < 1.0e-9);

        // The formula bowed_string.rs used before this bead: cancel the
        // current contact acceleration and correct the sub-step velocity
        // error linearly. Reproduced verbatim from the pre-fix source (see
        // frankensim commit e8589b00) to demonstrate its pinning error.
        let w_point: f64 = shapes.iter().map(|phi| phi * phi).sum();
        let v_rel = target_port_velocity_m_s
            - model
                .states()
                .iter()
                .zip(&shapes)
                .map(|(s, phi)| phi * s.velocity_m_sqrt_kg_per_s)
                .sum::<f64>();
        let accel_hold: f64 = model
            .modes()
            .iter()
            .zip(model.states())
            .zip(&shapes)
            .map(|((m, s), phi)| {
                phi * (2.0
                    * m.damping_ratio
                    * m.angular_frequency_rad_s
                    * s.velocity_m_sqrt_kg_per_s
                    + m.angular_frequency_rad_s
                        * m.angular_frequency_rad_s
                        * s.displacement_m_sqrt_kg)
            })
            .sum::<f64>()
            / w_point;
        let old_force = accel_hold + v_rel / (w_point * duration_s);
        let mut old_probe = model.clone();
        old_probe
            .step_duration(&[shapes[0] * old_force], duration_s)
            .unwrap();
        let old_velocity = old_probe.states()[0].velocity_m_sqrt_kg_per_s;

        // The old heuristic materially misses the endpoint target over this
        // non-small step; the new exact solve does not.
        assert!(
            (old_velocity - target_port_velocity_m_s).abs() > 1.0e-4,
            "old_velocity={old_velocity} target={target_port_velocity_m_s}"
        );
    }

    #[test]
    fn g0_held_force_for_port_velocity_refuses_malformed_inputs_and_never_mutates() {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 500.0,
            damping_ratio: 0.1,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
        model
            .restore_states(&[ModalAcousticState {
                displacement_m_sqrt_kg: 1.0e-4,
                velocity_m_sqrt_kg_per_s: 0.02,
            }])
            .unwrap();
        let before = model.states().to_vec();
        let dt = model.sample_period_s();

        assert!(matches!(
            model.held_force_for_port_velocity(&[1.0], 0.5, 0.0),
            Err(ModalAcousticTimeError::InvalidInput {
                what: "held-force duration must be positive and finite"
            })
        ));
        assert!(matches!(
            model.held_force_for_port_velocity(&[1.0], 0.5, -dt),
            Err(ModalAcousticTimeError::InvalidInput {
                what: "held-force duration must be positive and finite"
            })
        ));
        assert!(matches!(
            model.held_force_for_port_velocity(&[1.0], 0.5, f64::NAN),
            Err(ModalAcousticTimeError::InvalidInput {
                what: "held-force duration must be positive and finite"
            })
        ));
        assert!(matches!(
            model.held_force_for_port_velocity(&[], 0.5, dt),
            Err(ModalAcousticTimeError::ForceCountMismatch {
                expected: 1,
                found: 0
            })
        ));
        assert!(matches!(
            model.held_force_for_port_velocity(&[f64::NAN], 0.5, dt),
            Err(ModalAcousticTimeError::InvalidInput { .. })
        ));
        assert!(matches!(
            model.held_force_for_port_velocity(&[1.0], f64::NAN, dt),
            Err(ModalAcousticTimeError::InvalidInput { .. })
        ));

        // An undamped port sampled after exactly 3/4 of one period has
        // strictly negative unit-force mobility (sin(3*pi/2) = -1): a held
        // force there does not move the endpoint velocity the way the port
        // needs, and the solve must refuse rather than return a poisoned
        // answer.
        let undamped = ModalAcousticMode {
            angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 440.0,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        };
        let negative_mobility_model =
            ModalAcousticTimeModel::try_new(48_000, vec![undamped], budget()).unwrap();
        let three_quarter_period_s = 1.5 * core::f64::consts::PI / undamped.angular_frequency_rad_s;
        assert!(matches!(
            negative_mobility_model.held_force_for_port_velocity(
                &[1.0],
                0.5,
                three_quarter_period_s
            ),
            Err(ModalAcousticTimeError::InvalidInput {
                what: "port mobility must be positive and finite"
            })
        ));

        assert_eq!(
            model.states(),
            before,
            "every refusal must leave state untouched"
        );
    }

    #[test]
    fn probe_linear_map_matches_step_for_oscillatory_and_free_coordinates() {
        let modes = vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            },
            ModalAcousticMode {
                angular_frequency_rad_s: 2.0 * core::f64::consts::PI * 880.0,
                damping_ratio: 0.015,
                pressure_per_modal_velocity: C64::ZERO,
            },
        ];
        let mut model =
            ModalAcousticTimeModel::try_new_with_free_coordinates(48_000, modes, budget()).unwrap();
        let initial = [
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.002,
                velocity_m_sqrt_kg_per_s: -0.35,
            },
            ModalAcousticState {
                displacement_m_sqrt_kg: -0.0004,
                velocity_m_sqrt_kg_per_s: 0.8,
            },
        ];
        model.restore_states(&initial).unwrap();
        let forces = [1.25, -3.5];
        let m0 = model.probe_linear_map(0);
        let m1 = model.probe_linear_map(1);
        model.step(&forces).unwrap();
        let after = model.states();
        for (idx, (q0, f, (qq, qv, vq, vv, bq, bv), s1)) in [
            (0usize, initial[0], forces[0], m0, after[0]),
            (1usize, initial[1], forces[1], m1, after[1]),
        ]
        .into_iter()
        .map(|(i, q0, f, m, s1)| (i, (q0, f, m, s1)))
        {
            let pred_q = qq * q0.displacement_m_sqrt_kg + qv * q0.velocity_m_sqrt_kg_per_s + bq * f;
            let pred_v = vq * q0.displacement_m_sqrt_kg + vv * q0.velocity_m_sqrt_kg_per_s + bv * f;
            assert!(
                (pred_q - s1.displacement_m_sqrt_kg).abs() <= 1e-15,
                "mode {idx} displacement mismatch: {pred_q} vs {}",
                s1.displacement_m_sqrt_kg
            );
            assert!(
                (pred_v - s1.velocity_m_sqrt_kg_per_s).abs() <= 1e-12,
                "mode {idx} velocity mismatch: {pred_v} vs {}",
                s1.velocity_m_sqrt_kg_per_s
            );
        }
    }

    // Held-force work over the displacement increment the transition itself added to q0.
    //
    // The four steps that first refused "independent viscous heat and held-force energy balance
    // disagree" at 95b9c9e5a, when the work was formed as force * (next - state): the old-law cello
    // C2 and double-bass E1 open-string cards and, under the hyperbolic bow law, cello C2 and
    // double-bass A1 (partial 3 or 5 of a 65-205 Hz string at the 768 kHz bow substep). Every input
    // is the exact f64 captured at the refusal, as bits. `exact_work` and `exact_heat` come from an
    // independent 80-digit exponential of [[0,1,0],[-w^2,-2 zeta w,F],[0,0,0]] acting on [q0,v0,1],
    // with those inputs lifted exactly: work = F (q1 - q0), heat = integral of 2 zeta w v^2.
    struct RefusedStep {
        omega: u64,
        force: u64,
        q0: u64,
        v0: u64,
        exact_work: f64,
        exact_heat: f64,
    }
    const REFUSED_ZETA_BITS: u64 = 0x3f53_fafe_e44c_1a60;
    const BOW_SUBSTEP_BITS: u64 = 0x3eb5_d867_c3ec_e2a5;
    const REFUSED_STEPS: [RefusedStep; 4] = [
        RefusedStep {
            omega: 0x4093_45aa_98f2_4490,
            force: 0x403f_024c_8db4_1037,
            q0: 0x3e51_0e23_7118_2ec8,
            v0: 0xbf01_481e_91e8_8b32,
            exact_work: -5.164_046_616_300_688_5e-10,
            exact_heat: 1.172_078_936_764_100_5e-15,
        },
        RefusedStep {
            omega: 0x4094_36bc_910f_c724,
            force: 0x4043_f4ce_ee39_e398,
            q0: 0x3e64_e504_1662_8a0c,
            v0: 0xbef7_368d_b8f6_02b8,
            exact_work: 1.977_238_145_884_529_9e-10,
            exact_heat: 9.811_487_601_742_133_9e-16,
        },
        RefusedStep {
            omega: 0x4093_45aa_98f2_4490,
            force: 0x4034_5047_8abb_7a73,
            q0: 0xbe62_e8d0_d272_3d93,
            v0: 0xbeed_8580_96f1_6e38,
            exact_work: -2.161_101_970_446_351_2e-11,
            exact_heat: 2.321_886_795_679_498_0e-16,
        },
        RefusedStep {
            omega: 0x409a_ff80_fd79_77af,
            force: 0x4042_c7f3_5ffe_981f,
            q0: 0x3e52_aac7_9e51_f943,
            v0: 0xbefb_3424_831b_5b3a,
            exact_work: -7.447_397_818_591_547_4e-11,
            exact_heat: 1.103_567_335_447_121_4e-15,
        },
    ];

    fn refused_step_model(step: &RefusedStep) -> (ModalAcousticTimeModel, f64) {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: f64::from_bits(step.omega),
            damping_ratio: f64::from_bits(REFUSED_ZETA_BITS),
            pressure_per_modal_velocity: C64::ZERO,
        };
        let mut model = ModalAcousticTimeModel::try_new(768_000, vec![mode], budget()).unwrap();
        model
            .restore_states(&[ModalAcousticState {
                displacement_m_sqrt_kg: f64::from_bits(step.q0),
                velocity_m_sqrt_kg_per_s: f64::from_bits(step.v0),
            }])
            .unwrap();
        assert_eq!(model.sample_period_s().to_bits(), BOW_SUBSTEP_BITS);
        (model, f64::from_bits(step.force))
    }

    #[test]
    fn g1_held_force_work_is_the_transition_increment_on_the_steps_that_refused() {
        // Red before the change: each of these steps was refused at the unchanged tolerance.
        for (index, step) in REFUSED_STEPS.iter().enumerate() {
            let (mut model, force) = refused_step_model(step);
            let dt = model.sample_period_s();
            let before = modal_energy(model.modes[0], model.states()[0]);
            let frame = model
                .step_duration(&[force], dt)
                .unwrap_or_else(|error| panic!("step {index} refused: {error:?}"));
            // The balance compares the work with the stored energies, so the work is judged on the
            // check's own scale, eps * max(|work|, E), not on |work|, which can be far below E.
            let scale = f64::EPSILON
                * step.exact_work.abs().max(before).max(frame.total_modal_energy_j);
            let work_error = (frame.input_work_j - step.exact_work).abs();
            assert!(
                work_error <= 4.0 * scale,
                "step {index}: work {:e} against exact {:e}, error {:.3} x eps*max(|work|, E)",
                frame.input_work_j,
                step.exact_work,
                work_error / scale
            );
            println!(
                "{{\"test\":\"work_on_refused_steps\",\"step\":{index},\"work_error_over_eps_scale\":{:.3}}}",
                work_error / scale
            );
            let heat_error = (frame.viscous_dissipation_j - step.exact_heat).abs() / step.exact_heat;
            assert!(
                heat_error <= 8.0 * f64::EPSILON,
                "step {index}: heat {:e} against exact {:e}",
                frame.viscous_dissipation_j,
                step.exact_heat
            );
        }
    }

    #[test]
    fn g1_a_planted_heat_error_of_one_part_in_ten_thousand_still_refuses() {
        for (index, step) in REFUSED_STEPS.iter().enumerate() {
            let (mut model, force) = refused_step_model(step);
            let dt = model.sample_period_s();
            // The same cached transition, but a viscous integral prepared for a damping ratio 1e-4
            // larger: the heat the work/energy balance is compared against is wrong by 1e-4 of itself.
            let mut planted_mode = model.modes[0];
            planted_mode.damping_ratio *= 1.0 + 1.0e-4;
            let mut kernel = ExactZohKernel::new(&model.modes, dt);
            let mut full = kernel.full.to_vec();
            full[0].viscous_integral = ViscousIntegral::new(planted_mode, dt);
            kernel.full = full.into_boxed_slice();
            model.transition_kernel = Arc::new(kernel);
            match model.step_duration(&[force], dt) {
                Err(ModalAcousticTimeError::InvalidInput { what }) => assert_eq!(
                    what, "independent viscous heat and held-force energy balance disagree",
                    "step {index}"
                ),
                other => panic!("step {index}: a heat error of 1e-4 must refuse, got {other:?}"),
            }
        }
    }

    // The per-mode balance test of `step_duration_into`, written the same way with the same
    // tolerance: |balance_loss - loss| / allowance, so a value above 1 refuses. It exists to measure
    // how large a displacement error the balance can see. The two tests above run the real check;
    // the assertions below tie this restatement to it on the same four steps, in both directions.
    fn balance_margin(
        mode: ModalAcousticMode,
        state: ModalAcousticState,
        next: ModalAcousticState,
        force: f64,
        work: f64,
        prepared: &PreparedZohCoefficients,
        duration_s: f64,
    ) -> f64 {
        let before = modal_energy(mode, state);
        let after = modal_energy(mode, next);
        let (loss, allowance) = prepared
            .viscous_integral
            .evaluate(mode, state, force, duration_s);
        let balance_loss = work - (after - before);
        let balance_allowance =
            256.0 * f64::EPSILON * work.abs().max(before).max(after).max(f64::MIN_POSITIVE)
                + allowance;
        (balance_loss - loss).abs() / balance_allowance
    }

    #[test]
    fn g1_a_displacement_error_far_above_roundoff_still_refuses() {
        let mut smallest = Vec::new();
        for (index, step) in REFUSED_STEPS.iter().enumerate() {
            let (model, force) = refused_step_model(step);
            let (mode, state, dt) = (model.modes[0], model.states()[0], model.sample_period_s());
            let prepared = PreparedZohCoefficients::new(mode, dt);
            let (next, increment) =
                advance_exact_zoh_increment(mode, state, force, prepared.transition);
            // Work from the stored positions, as before the change: refused, as the real check did.
            let from_positions =
                force * (next.displacement_m_sqrt_kg - state.displacement_m_sqrt_kg);
            assert!(
                balance_margin(mode, state, next, force, from_positions, &prepared, dt) > 1.0,
                "step {index}"
            );
            // Work over the carried increment: accepted, as the real check now does.
            let carried = force * increment;
            let margin = |relative: f64| {
                let planted = ModalAcousticState {
                    displacement_m_sqrt_kg: next.displacement_m_sqrt_kg * (1.0 + relative),
                    ..next
                };
                balance_margin(mode, state, planted, force, carried, &prepared, dt)
            };
            assert!(margin(0.0) < 1.0, "step {index}");
            // The work no longer reads the stored displacement, so the balance sees an error in it
            // only through the stored energy, w^2 |q| delta, where it used to see |F - w^2 q| delta.
            // An error of 1e-12 of q still refuses, in either direction.
            assert!(margin(1.0e-12) > 1.0 && margin(-1.0e-12) > 1.0, "step {index}");
            let mut relative = 1.0e-17;
            while relative < 1.0e-3 && !(margin(relative) > 1.0 && margin(-relative) > 1.0) {
                relative *= 1.25;
            }
            smallest.push(relative);
        }
        println!(
            "{{\"test\":\"displacement_error_sensitivity\",\"smallest_refused_relative_error\":{smallest:?}}}"
        );
        assert!(smallest.iter().all(|relative| *relative <= 2.0e-13), "{smallest:?}");
    }

    #[test]
    fn g1_the_returned_increment_is_the_addend_of_the_displacement_it_forms() {
        let dt = 1.0 / 48_000.0;
        let state = ModalAcousticState {
            displacement_m_sqrt_kg: 3.0e-6,
            velocity_m_sqrt_kg_per_s: -2.0e-3,
        };
        // Free, undamped, damped short step, damped long step (omega dt = 2.5).
        for (omega, zeta, long_step) in [
            (0.0, 0.0, false),
            (2000.0, 0.0, false),
            (2000.0, 1.0e-3, false),
            (120_000.0, 0.02, true),
        ] {
            let mode = ModalAcousticMode {
                angular_frequency_rad_s: omega,
                damping_ratio: zeta,
                pressure_per_modal_velocity: C64::ZERO,
            };
            let coefficients = exact_zoh_coefficients(mode, dt);
            let (next, increment) = advance_exact_zoh_increment(mode, state, 3.0, coefficients);
            assert_state_bits_eq(next, advance_exact_zoh_cached(mode, state, 3.0, coefficients));
            let q0 = state.displacement_m_sqrt_kg;
            if long_step {
                assert_eq!(increment.to_bits(), (next.displacement_m_sqrt_kg - q0).to_bits());
            } else {
                assert_eq!(next.displacement_m_sqrt_kg.to_bits(), (q0 + increment).to_bits());
            }
        }
    }
}
