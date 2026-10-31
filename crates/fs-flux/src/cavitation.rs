//! Spherical Rayleigh--Plesset bubble dynamics in a frozen Newtonian liquid.
//!
//! `R R'' + 3 R'^2/2 = (p_v + p_g - p_inf - 2 sigma/R - 4 mu R'/R)/rho`.
//! The noncondensable gas follows the explicitly supplied polytropic law;
//! vapor partial pressure, density, viscosity, surface tension and sound speed
//! come from one [`CavitationState`]. The imposed liquid pressure uses an
//! absolute reference and may be tensile. It is never clamped to vapor pressure.
//!
//! See [Brennen, chapters 2--4](https://media.library.caltech.edu/CaltechBOOK:1995.001/chap2.htm).
//! This incompressible, spherical model excludes liquid acoustic back-reaction,
//! shocks, walls, bubble interactions, gas diffusion and evolving interface
//! temperature. Vapor remains at the supplied saturation pressure; phase mass
//! and latent heat are not silently inferred from geometric volume. The small
//! vapor/liquid density-ratio approximation to interface kinematics is checked
//! when vapor density is supplied; otherwise it remains an explicit caller
//! applicability condition. A radius or wall-Mach boundary refuses the
//! advance instead of manufacturing a minimum-radius rebound.

use core::f64::consts::PI;
use core::fmt;
use fs_material::cavitation::CavitationState;

/// Invalid input, an explicit model boundary, or an unsuccessful integration.
#[derive(Clone, Debug, PartialEq)]
pub enum CavitationError {
    /// A named input is nonfinite or outside its admitted range.
    InvalidInput(&'static str),
    /// The radius left the caller's spherical-continuum domain.
    RadiusOutsideDomain {
        radius_m: f64,
        minimum_radius_m: f64,
        maximum_radius_m: f64,
    },
    /// The incompressible liquid approximation no longer applies.
    WallMachExceeded {
        wall_mach_number: f64,
        maximum_wall_mach_number: f64,
    },
    /// Neglecting interfacial liquid/vapor mass-flux kinematics is inadmissible.
    VaporDensityRatioExceeded {
        vapor_liquid_density_ratio: f64,
        maximum_vapor_density_ratio: f64,
    },
    /// A physical or numerical derived quantity is not representable.
    NonFiniteDerived(&'static str),
    /// All allowed trial steps were consumed before the requested endpoint.
    StepBudgetExceeded {
        attempts: usize,
        accepted_steps: usize,
    },
    /// The requested local error cannot be met at the minimum step.
    ToleranceNotMet { normalized_error: f64, step_s: f64 },
    /// Advancing the floating-point clock is no longer possible.
    StepUnderflow,
    /// Cancellation was observed before publishing the complete interval.
    Cancelled,
}

impl fmt::Display for CavitationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "spherical cavitation calculation refused: {self:?}")
    }
}
impl std::error::Error for CavitationError {}

/// Fixed noncondensable-gas inventory represented by `p_g R^(3 kappa)`.
/// This is separate from the liquid's saturation vapor pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleGas {
    reference_radius_m: f64,
    reference_pressure_pa: f64,
    polytropic_exponent: f64,
}

impl BubbleGas {
    /// The gas pressure may be zero for a vapor-only Rayleigh limit. The
    /// positive exponent is caller supplied: 1 is isothermal, and an adiabatic
    /// choice requires the gas's own heat-capacity ratio.
    pub fn try_new(
        reference_radius_m: f64,
        reference_pressure_pa: f64,
        polytropic_exponent: f64,
    ) -> Result<Self, CavitationError> {
        if !positive(reference_radius_m)
            || !nonnegative(reference_pressure_pa)
            || !positive(polytropic_exponent)
            || !(3.0 * polytropic_exponent).is_finite()
        {
            return Err(CavitationError::InvalidInput("polytropic gas reference"));
        }
        Ok(Self {
            reference_radius_m,
            reference_pressure_pa,
            polytropic_exponent,
        })
    }

    /// Construct the gas pressure needed for a motionless reference bubble:
    /// `p_g0 = p_inf0 - p_v + 2 sigma/R0`. A negative result cannot represent
    /// this fixed-gas equilibrium and is refused.
    pub fn equilibrium(
        liquid: CavitationState,
        reference_radius_m: f64,
        absolute_pressure_pa: f64,
        polytropic_exponent: f64,
    ) -> Result<Self, CavitationError> {
        if !positive(reference_radius_m) || !absolute_pressure_pa.is_finite() {
            return Err(CavitationError::InvalidInput("bubble equilibrium"));
        }
        Self::try_new(
            reference_radius_m,
            absolute_pressure_pa - liquid.saturation_pressure_pa()
                + 2.0 * liquid.surface_tension_n_m() / reference_radius_m,
            polytropic_exponent,
        )
    }

    pub const fn reference_radius_m(self) -> f64 {
        self.reference_radius_m
    }
    pub const fn reference_pressure_pa(self) -> f64 {
        self.reference_pressure_pa
    }
    pub const fn polytropic_exponent(self) -> f64 {
        self.polytropic_exponent
    }

    fn pressure(self, radius_m: f64) -> Result<f64, CavitationError> {
        // A vapor-only bubble must not evaluate an overflowing gas power and
        // then multiply it by zero during collapse.
        let pressure = if self.reference_pressure_pa == 0.0 {
            0.0
        } else {
            self.reference_pressure_pa
                * (self.reference_radius_m / radius_m).powf(3.0 * self.polytropic_exponent)
        };
        finite(pressure, "gas partial pressure")
    }
}

/// Caller-declared spherical-continuum and incompressibility boundaries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleLimits {
    /// Smallest admitted radius [m], strictly positive; no state is clamped.
    pub minimum_radius_m: f64,
    /// Largest admitted radius [m], chosen for the surrounding geometry.
    pub maximum_radius_m: f64,
    /// Largest `abs(R')/c_l`; must lie in `(0, 0.1]` for this model.
    pub maximum_wall_mach_number: f64,
    /// Largest admitted `rho_v/rho_l`, in `(0, 0.1]`. Checked whenever the
    /// material supplies vapor density. Without that optional datum, the caller
    /// explicitly assumes the neglected interfacial mass flux is small.
    pub maximum_vapor_density_ratio: f64,
}

impl BubbleLimits {
    fn validate(self) -> Result<(), CavitationError> {
        if !positive(self.minimum_radius_m)
            || !positive(self.maximum_radius_m)
            || self.maximum_radius_m <= self.minimum_radius_m
            || !positive(self.maximum_wall_mach_number)
            || self.maximum_wall_mach_number > 0.1
            || !positive(self.maximum_vapor_density_ratio)
            || self.maximum_vapor_density_ratio > 0.1
        {
            return Err(CavitationError::InvalidInput("bubble applicability limits"));
        }
        Ok(())
    }
}

/// Accepted physical bubble state. Cloning it is a complete physical checkpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleState {
    time_s: f64,
    radius_m: f64,
    wall_velocity_m_s: f64,
}

impl BubbleState {
    pub fn try_new(
        time_s: f64,
        radius_m: f64,
        wall_velocity_m_s: f64,
    ) -> Result<Self, CavitationError> {
        if !time_s.is_finite() || !positive(radius_m) || !wall_velocity_m_s.is_finite() {
            return Err(CavitationError::InvalidInput("bubble state"));
        }
        Ok(Self {
            time_s,
            radius_m,
            wall_velocity_m_s,
        })
    }
    pub const fn time_s(self) -> f64 {
        self.time_s
    }
    pub const fn radius_m(self) -> f64 {
        self.radius_m
    }
    pub const fn wall_velocity_m_s(self) -> f64 {
        self.wall_velocity_m_s
    }
}

/// Local embedded-error controls and a strict trial-step work budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleNumerics {
    /// Dimensionless relative tolerance in `[0, 1]`.
    pub relative_tolerance: f64,
    /// Positive radius error scale [m].
    pub radius_absolute_tolerance_m: f64,
    /// Positive wall-velocity error scale [m/s].
    pub velocity_absolute_tolerance_m_s: f64,
    pub initial_step_s: f64,
    pub minimum_step_s: f64,
    /// Also bounds pressure-waveform sampling intervals. The caller must resolve
    /// its drive; sampled local error cannot discover an arbitrary missed pulse.
    pub maximum_step_s: f64,
    /// Counts both accepted and rejected trial steps.
    pub maximum_attempts: usize,
}

impl BubbleNumerics {
    fn validate(self) -> Result<(), CavitationError> {
        if !nonnegative(self.relative_tolerance)
            || self.relative_tolerance > 1.0
            || !positive(self.radius_absolute_tolerance_m)
            || !positive(self.velocity_absolute_tolerance_m_s)
            || !positive(self.minimum_step_s)
            || !positive(self.maximum_step_s)
            || !positive(self.initial_step_s)
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || self.maximum_attempts == 0
            || self.maximum_attempts > (usize::MAX - 2) / 7
        {
            return Err(CavitationError::InvalidInput("bubble integration controls"));
        }
        Ok(())
    }
}

/// Instantaneous physical source and mechanical quantities at an accepted state.
/// `V''` supplies a compact monopole source; it is not microphone pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleObservation {
    pub state: BubbleState,
    pub imposed_absolute_pressure_pa: f64,
    pub gas_partial_pressure_pa: f64,
    pub wall_acceleration_m_s2: f64,
    pub wall_mach_number: f64,
    pub volume_m3: f64,
    pub volume_rate_m3_s: f64,
    pub volume_acceleration_m3_s2: f64,
    /// Kinetic energy of the infinite incompressible exterior liquid [J].
    pub liquid_kinetic_energy_j: f64,
    pub surface_energy_j: f64,
    /// Nonnegative instantaneous `16 pi mu R R'^2` [W], not an integrated heat.
    pub viscous_dissipation_w: f64,
}

/// Successful complete-interval advance; errors never modify the caller state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleStep {
    pub observation: BubbleObservation,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub rhs_evaluations: usize,
    /// Largest accepted embedded-error ratio, at most one. This is a local
    /// numerical indicator, not a bound on global or physical-model error.
    pub maximum_normalized_local_error: f64,
}

/// Immutable physical model; runtime radius/velocity live in [`BubbleState`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighPlessetBubble {
    liquid: CavitationState,
    gas: BubbleGas,
    limits: BubbleLimits,
}

impl RayleighPlessetBubble {
    pub fn try_new(
        liquid: CavitationState,
        gas: BubbleGas,
        limits: BubbleLimits,
    ) -> Result<Self, CavitationError> {
        limits.validate()?;
        if let Some(vapor) = liquid.saturated_vapor() {
            let ratio = vapor.density_kg_m3() / liquid.liquid_density_kg_m3();
            if ratio > limits.maximum_vapor_density_ratio {
                return Err(CavitationError::VaporDensityRatioExceeded {
                    vapor_liquid_density_ratio: ratio,
                    maximum_vapor_density_ratio: limits.maximum_vapor_density_ratio,
                });
            }
        }
        Ok(Self {
            liquid,
            gas,
            limits,
        })
    }
    pub const fn liquid(self) -> CavitationState {
        self.liquid
    }
    pub const fn gas(self) -> BubbleGas {
        self.gas
    }
    pub const fn limits(self) -> BubbleLimits {
        self.limits
    }

    /// Evaluate the physical acceleration and source without advancing state.
    pub fn observe(
        &self,
        state: &BubbleState,
        absolute_pressure_pa: f64,
    ) -> Result<BubbleObservation, CavitationError> {
        let (gas_partial_pressure_pa, acceleration) = self.dynamics(state, absolute_pressure_pa)?;
        let radius = state.radius_m;
        let velocity = state.wall_velocity_m_s;
        let area = 4.0 * PI * radius * radius;
        let volume = area * (radius / 3.0);
        if !positive(area) || !positive(volume) {
            return Err(CavitationError::NonFiniteDerived("bubble area or volume"));
        }
        Ok(BubbleObservation {
            state: *state,
            imposed_absolute_pressure_pa: absolute_pressure_pa,
            gas_partial_pressure_pa,
            wall_acceleration_m_s2: acceleration,
            wall_mach_number: velocity.abs() / self.liquid.liquid_sound_speed_m_s(),
            volume_m3: volume,
            volume_rate_m3_s: finite(area * velocity, "volume rate")?,
            volume_acceleration_m3_s2: finite(
                4.0 * PI * radius * (2.0 * velocity * velocity + radius * acceleration),
                "volume acceleration",
            )?,
            liquid_kinetic_energy_j: finite(
                1.5 * self.liquid.liquid_density_kg_m3() * volume * velocity * velocity,
                "exterior kinetic energy",
            )?,
            surface_energy_j: finite(area * self.liquid.surface_tension_n_m(), "surface energy")?,
            viscous_dissipation_w: finite(
                16.0 * PI * self.liquid.dynamic_viscosity_pa_s() * radius * velocity * velocity,
                "viscous dissipation",
            )?,
        })
    }

    /// Advance under a deterministic, bounded-work pressure callback [Pa].
    /// Split at known drive discontinuities and use `maximum_step_s` to resolve
    /// its bandwidth. The pressure may be negative under metastable tension.
    pub fn advance<Pressure: Fn(f64) -> f64>(
        &self,
        state: &mut BubbleState,
        end_time_s: f64,
        numerics: &BubbleNumerics,
        pressure: &Pressure,
    ) -> Result<BubbleStep, CavitationError> {
        self.advance_controlled(state, end_time_s, numerics, pressure, &mut || false)
    }

    /// Dormand--Prince 5(4), fixed two-component stack storage and reuse of the
    /// accepted endpoint derivative. Cancellation is polled before each RHS and
    /// final publication. Failed stages can refine the step, but cannot relax
    /// physical limits. The complete requested interval publishes atomically.
    /// Caller-owned pressure callbacks must themselves bound their work.
    pub fn advance_controlled<Pressure, Cancel>(
        &self,
        state: &mut BubbleState,
        end_time_s: f64,
        numerics: &BubbleNumerics,
        pressure: &Pressure,
        cancelled: &mut Cancel,
    ) -> Result<BubbleStep, CavitationError>
    where
        Pressure: Fn(f64) -> f64,
        Cancel: FnMut() -> bool,
    {
        poll(cancelled)?;
        numerics.validate()?;
        if !end_time_s.is_finite()
            || end_time_s < state.time_s
            || !(end_time_s - state.time_s).is_finite()
        {
            return Err(CavitationError::InvalidInput("forward time interval"));
        }
        let mut current = *state;
        let mut rhs_evaluations = 0;
        let mut derivative = self.rhs(current, pressure, cancelled, &mut rhs_evaluations)?;
        let mut step = numerics.initial_step_s;
        let mut accepted_steps = 0;
        let mut rejected_steps = 0;
        let mut maximum_normalized_local_error: f64 = 0.0;
        while current.time_s < end_time_s {
            if accepted_steps + rejected_steps == numerics.maximum_attempts {
                return Err(CavitationError::StepBudgetExceeded {
                    attempts: numerics.maximum_attempts,
                    accepted_steps,
                });
            }
            poll(cancelled)?;
            let remaining = end_time_s - current.time_s;
            step = step.min(remaining);
            let next_time = if step == remaining {
                end_time_s
            } else {
                current.time_s + step
            };
            if next_time <= current.time_s || next_time > end_time_s {
                return Err(CavitationError::StepUnderflow);
            }
            step = next_time - current.time_s;
            // Compare represented endpoint times. Subtracting a rounded clock
            // can make an already minimum-size step look slightly larger than
            // the minimum and otherwise retry the same failed trial forever.
            let can_refine = next_time > current.time_s + numerics.minimum_step_s;
            let trial = self.trial(
                current,
                derivative,
                next_time,
                step,
                numerics,
                pressure,
                cancelled,
                &mut rhs_evaluations,
            );
            let (next, endpoint_derivative, error) = match trial {
                Ok(trial) => trial,
                Err(error) if refinable(&error) && can_refine => {
                    rejected_steps += 1;
                    step = (0.5 * step).max(numerics.minimum_step_s);
                    continue;
                }
                Err(error) => return Err(error),
            };
            if error > 1.0 {
                rejected_steps += 1;
                if !can_refine {
                    return Err(CavitationError::ToleranceNotMet {
                        normalized_error: error,
                        step_s: step,
                    });
                }
            } else {
                current = next;
                derivative = endpoint_derivative;
                accepted_steps += 1;
                maximum_normalized_local_error = maximum_normalized_local_error.max(error);
            }
            let factor = if error == 0.0 {
                5.0
            } else {
                (0.9 * error.powf(-0.2)).clamp(0.2, 5.0)
            };
            step = (step * factor).clamp(numerics.minimum_step_s, numerics.maximum_step_s);
        }
        poll(cancelled)?;
        let observation = self.observe(&current, pressure(current.time_s))?;
        poll(cancelled)?;
        *state = current;
        Ok(BubbleStep {
            observation,
            accepted_steps,
            rejected_steps,
            rhs_evaluations,
            maximum_normalized_local_error,
        })
    }

    fn dynamics(&self, state: &BubbleState, pressure: f64) -> Result<(f64, f64), CavitationError> {
        let radius = state.radius_m;
        if !radius.is_finite()
            || radius < self.limits.minimum_radius_m
            || radius > self.limits.maximum_radius_m
        {
            return Err(CavitationError::RadiusOutsideDomain {
                radius_m: radius,
                minimum_radius_m: self.limits.minimum_radius_m,
                maximum_radius_m: self.limits.maximum_radius_m,
            });
        }
        let velocity = state.wall_velocity_m_s;
        let mach = velocity.abs() / self.liquid.liquid_sound_speed_m_s();
        if !mach.is_finite() || mach > self.limits.maximum_wall_mach_number {
            return Err(CavitationError::WallMachExceeded {
                wall_mach_number: mach,
                maximum_wall_mach_number: self.limits.maximum_wall_mach_number,
            });
        }
        if !pressure.is_finite() {
            return Err(CavitationError::InvalidInput("imposed absolute pressure"));
        }
        let gas_pressure = self.gas.pressure(radius)?;
        // This ordering exactly preserves the equilibrium constructor's
        // pressure balance at R0, and includes both capillary and viscous stress.
        let driving_pressure = gas_pressure
            - (pressure - self.liquid.saturation_pressure_pa()
                + 2.0 * self.liquid.surface_tension_n_m() / radius)
            - 4.0 * self.liquid.dynamic_viscosity_pa_s() * velocity / radius;
        let acceleration = finite(
            (driving_pressure / self.liquid.liquid_density_kg_m3() - 1.5 * velocity * velocity)
                / radius,
            "wall acceleration",
        )?;
        Ok((gas_pressure, acceleration))
    }

    fn rhs<Pressure: Fn(f64) -> f64, Cancel: FnMut() -> bool>(
        &self,
        state: BubbleState,
        pressure: &Pressure,
        cancelled: &mut Cancel,
        evaluations: &mut usize,
    ) -> Result<[f64; 2], CavitationError> {
        poll(cancelled)?;
        *evaluations += 1;
        let (_, acceleration) = self.dynamics(&state, pressure(state.time_s))?;
        Ok([state.wall_velocity_m_s, acceleration])
    }

    #[allow(clippy::too_many_arguments)]
    fn trial<Pressure: Fn(f64) -> f64, Cancel: FnMut() -> bool>(
        &self,
        state: BubbleState,
        first_derivative: [f64; 2],
        next_time: f64,
        step: f64,
        numerics: &BubbleNumerics,
        pressure: &Pressure,
        cancelled: &mut Cancel,
        evaluations: &mut usize,
    ) -> Result<(BubbleState, [f64; 2], f64), CavitationError> {
        let mut stages = [[0.0; 2]; 7];
        stages[0] = first_derivative;
        let initial = [state.radius_m, state.wall_velocity_m_s];
        let mut values = initial;
        for i in 0..6 {
            for component in 0..2 {
                let mut sum = 0.0;
                for (j, stage) in stages.iter().enumerate().take(i + 1) {
                    sum += A[i][j] * stage[component];
                }
                values[component] = initial[component] + step * sum;
            }
            let stage_time = if C[i] == 1.0 {
                next_time
            } else {
                state.time_s + C[i] * step
            };
            stages[i + 1] = self.rhs(
                BubbleState {
                    time_s: stage_time,
                    radius_m: values[0],
                    wall_velocity_m_s: values[1],
                },
                pressure,
                cancelled,
                evaluations,
            )?;
        }
        let absolute = [
            numerics.radius_absolute_tolerance_m,
            numerics.velocity_absolute_tolerance_m_s,
        ];
        let mut error: f64 = 0.0;
        for component in 0..2 {
            let mut difference = 0.0;
            for (stage, weight) in stages.iter().zip(ERROR_WEIGHT) {
                difference += weight * stage[component];
            }
            let scale = finite(
                absolute[component]
                    + numerics.relative_tolerance
                        * initial[component].abs().max(values[component].abs()),
                "local error scale",
            )?;
            error = error.max(finite(
                (step * difference).abs() / scale,
                "embedded local error",
            )?);
        }
        Ok((
            BubbleState {
                time_s: next_time,
                radius_m: values[0],
                wall_velocity_m_s: values[1],
            },
            stages[6],
            error,
        ))
    }
}

// Standard Dormand--Prince coefficients. The final A row is the fifth-order
// endpoint; ERROR_WEIGHT is its difference from the embedded fourth-order rule.
const C: [f64; 6] = [0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0];
const A: [[f64; 6]; 6] = [
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
    ],
];
const ERROR_WEIGHT: [f64; 7] = [
    35.0 / 384.0 - 5179.0 / 57600.0,
    0.0,
    500.0 / 1113.0 - 7571.0 / 16695.0,
    125.0 / 192.0 - 393.0 / 640.0,
    -2187.0 / 6784.0 + 92097.0 / 339200.0,
    11.0 / 84.0 - 187.0 / 2100.0,
    -1.0 / 40.0,
];

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64, name: &'static str) -> Result<f64, CavitationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CavitationError::NonFiniteDerived(name))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CavitationError> {
    if cancelled() {
        Err(CavitationError::Cancelled)
    } else {
        Ok(())
    }
}
fn refinable(error: &CavitationError) -> bool {
    matches!(
        error,
        CavitationError::RadiusOutsideDomain { .. }
            | CavitationError::WallMachExceeded { .. }
            | CavitationError::NonFiniteDerived(_)
    )
}

#[cfg(test)]
#[path = "cavitation/tests.rs"]
mod tests;
