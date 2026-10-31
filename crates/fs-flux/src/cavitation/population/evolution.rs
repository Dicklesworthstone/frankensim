//! Nonlinear radial evolution of an original, fixed-number bubble population.
//!
//! Each class uses its original Keller--Miksis equation and adaptive stepper.
//! Absolute number density remains per fixed reference mixture volume; only
//! radius and wall velocity evolve. The aggregate alpha, alpha' and alpha''
//! are sums of the actual native volumes and their geometric derivatives.
//! Liquid displacement/transport, nucleation, breakup, coalescence, diffusion,
//! thermal exchange and inter-bubble pressure feedback are not supplied here.
//! The common prescribed pressure history is owned by the caller.
//!
//! Dilution is admitted for the entire interval by the conservative envelope
//! sum_j n_j (4 pi/3) R_max,j^3, using each original model's authored radius
//! limit. The native stepper checks those radii at every physical RHS stage;
//! checking only asynchronous class endpoints would miss transient growth.
//! This envelope does not certify weak cloud interaction or spatial compactness.
//!
//! A proposal uses preallocated scratch states. Inspect it before commit; on
//! any error, cancellation, or dropped proposal, all accepted states and their
//! common clock remain unchanged. No gas inventory or number weight is altered.

use super::super::{BubbleNumerics, BubblePressureSample, BubbleState, CavitationError};
use super::{BubblePopulation, Sum};

/// Population bounds in addition to every original class model's own domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PopulationEvolutionLimits {
    /// Largest aggregate volume fraction, in (0, 0.01]. The entire authored
    /// radius envelope must fit, not only the initial or final actual volume.
    pub maximum_void_fraction: f64,
}

/// Native numerical controls and one shared budget for all classes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PopulationEvolutionNumerics {
    pub per_class: BubbleNumerics,
    /// Maximum SUM of accepted and rejected native trials in this interval.
    /// Zero admits an empty population or a zero-duration observation only.
    pub maximum_total_attempts: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum PopulationEvolutionError {
    InvalidInput(&'static str),
    Input(CavitationError),
    Class {
        class_index: usize,
        cause: CavitationError,
    },
    StateCountMismatch {
        expected: usize,
        actual: usize,
    },
    ClockMismatch {
        class_index: usize,
        expected_time_s: f64,
        actual_time_s: f64,
    },
    DilutionEnvelopeExceeded {
        maximum_void_fraction_envelope: f64,
        maximum_void_fraction: f64,
    },
    DilutionExceeded {
        void_fraction: f64,
        maximum_void_fraction: f64,
    },
    WorkBudgetExceeded {
        class_index: usize,
        consumed_attempts: usize,
        maximum_total_attempts: usize,
    },
    Cancelled,
}
impl core::fmt::Display for PopulationEvolutionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "bubble population evolution refused: {self:?}")
    }
}
impl std::error::Error for PopulationEvolutionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Input(error) | Self::Class { cause: error, .. } => Some(error),
            _ => None,
        }
    }
}
impl From<CavitationError> for PopulationEvolutionError {
    fn from(error: CavitationError) -> Self {
        if error == CavitationError::Cancelled {
            Self::Cancelled
        } else {
            Self::Input(error)
        }
    }
}

/// Actual geometric source per unit fixed reference mixture volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PopulationObservation {
    pub time_s: f64,
    pub imposed_pressure: BubblePressureSample,
    pub occupied_classes: usize,
    pub number_density_per_m3: f64,
    pub void_fraction: f64,
    pub void_fraction_rate_per_s: f64,
    pub void_fraction_acceleration_per_s2: f64,
    /// sum_j n_j R_j [1/m²], useful for a separately admitted cloud model.
    pub number_radius_moment_per_m2: f64,
    pub maximum_radius_m: f64,
    pub maximum_wall_mach_number: f64,
}

/// Native work totals and one common accepted endpoint, without a class clone.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PopulationEvolutionStep {
    pub observation: PopulationObservation,
    pub advanced_classes: usize,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    /// Sum of the original steppers' RHS counters. Endpoint observation work
    /// is an additional bounded traversal, not included in those counters.
    pub rhs_evaluations: usize,
    pub maximum_normalized_local_error: f64,
}

/// Original models plus one accepted common clock and reusable candidate states.
#[derive(Debug)]
pub struct EvolvingBubblePopulation {
    population: BubblePopulation,
    states: Vec<BubbleState>,
    scratch: Vec<BubbleState>,
    time_s: f64,
    limits: PopulationEvolutionLimits,
    maximum_void_fraction_envelope: f64,
}

impl EvolvingBubblePopulation {
    pub fn try_new(
        population: BubblePopulation,
        time_s: f64,
        states: Vec<BubbleState>,
        limits: PopulationEvolutionLimits,
        initial_pressure: BubblePressureSample,
    ) -> Result<Self, PopulationEvolutionError> {
        if !time_s.is_finite()
            || !limits.maximum_void_fraction.is_finite()
            || limits.maximum_void_fraction <= 0.0
            || limits.maximum_void_fraction > 0.01
        {
            return Err(PopulationEvolutionError::InvalidInput(
                "initial clock or dilution limit",
            ));
        }
        if states.len() != population.classes().len() {
            return Err(PopulationEvolutionError::StateCountMismatch {
                expected: population.classes().len(),
                actual: states.len(),
            });
        }
        let mut envelope = Sum::default();
        for class in population.classes() {
            if class.number_density_per_m3() == 0.0 {
                continue;
            }
            let radius = class.bubble().limits().maximum_radius_m;
            let area = 4.0 * core::f64::consts::PI * radius * radius;
            let maximum_fraction = class.number_density_per_m3() * (area * (radius / 3.0));
            if !maximum_fraction.is_finite() || maximum_fraction <= 0.0 {
                return Err(PopulationEvolutionError::InvalidInput(
                    "unrepresentable radius envelope",
                ));
            }
            envelope.add(maximum_fraction);
        }
        let maximum_void_fraction_envelope = envelope.value();
        if !maximum_void_fraction_envelope.is_finite()
            || maximum_void_fraction_envelope > limits.maximum_void_fraction
        {
            return Err(PopulationEvolutionError::DilutionEnvelopeExceeded {
                maximum_void_fraction_envelope,
                maximum_void_fraction: limits.maximum_void_fraction,
            });
        }
        let owner = Self {
            population,
            scratch: states.clone(),
            states,
            time_s,
            limits,
            maximum_void_fraction_envelope,
        };
        owner.observe(initial_pressure)?;
        Ok(owner)
    }

    /// Initialize all slots at their authored reference radii and zero velocity.
    /// The supplied pressure is still checked; equilibrium is not manufactured.
    pub fn at_reference(
        population: BubblePopulation,
        time_s: f64,
        limits: PopulationEvolutionLimits,
        initial_pressure: BubblePressureSample,
    ) -> Result<Self, PopulationEvolutionError> {
        let states = population
            .classes()
            .iter()
            .map(|class| {
                BubbleState::try_new(time_s, class.bubble().gas().reference_radius_m(), 0.0)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::try_new(population, time_s, states, limits, initial_pressure)
    }

    pub const fn population(&self) -> &BubblePopulation {
        &self.population
    }
    pub fn states(&self) -> &[BubbleState] {
        &self.states
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn limits(&self) -> PopulationEvolutionLimits {
        self.limits
    }
    pub const fn maximum_void_fraction_envelope(&self) -> f64 {
        self.maximum_void_fraction_envelope
    }

    /// Read actual native sources at the current accepted common clock.
    pub fn observe(
        &self,
        pressure: BubblePressureSample,
    ) -> Result<PopulationObservation, PopulationEvolutionError> {
        self.observe_states(&self.states, self.time_s, pressure, &mut || false)
    }

    pub fn advance<Pressure: Fn(f64) -> BubblePressureSample>(
        &mut self,
        end_time_s: f64,
        numerics: &PopulationEvolutionNumerics,
        pressure: &Pressure,
    ) -> Result<PopulationEvolutionStep, PopulationEvolutionError> {
        self.advance_controlled(end_time_s, numerics, pressure, &mut || false)
    }

    pub fn advance_controlled<Pressure, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &PopulationEvolutionNumerics,
        pressure: &Pressure,
        cancelled: &mut Cancel,
    ) -> Result<PopulationEvolutionStep, PopulationEvolutionError>
    where
        Pressure: Fn(f64) -> BubblePressureSample,
        Cancel: FnMut() -> bool,
    {
        Ok(self
            .propose_advance_controlled(end_time_s, numerics, pressure, cancelled)?
            .commit())
    }

    /// Prepare all classes under the same smooth authored drive and derivative.
    /// No accepted state is changed until the returned proposal is committed.
    /// Per-class adaptive clocks can differ internally because the conservative
    /// radius envelope has already admitted dilution throughout the interval.
    ///
    /// The pressure callback must be deterministic in time, with the derivative
    /// of that same smooth history. Calls may revisit earlier times between
    /// classes and rejected trials. The native maximum step must resolve the
    /// authored drive; sampled local error cannot detect an arbitrary missed pulse.
    pub fn propose_advance_controlled<Pressure, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &PopulationEvolutionNumerics,
        pressure: &Pressure,
        cancelled: &mut Cancel,
    ) -> Result<PopulationAdvance<'_>, PopulationEvolutionError>
    where
        Pressure: Fn(f64) -> BubblePressureSample,
        Cancel: FnMut() -> bool,
    {
        poll(cancelled)?;
        numerics.per_class.validate()?;
        if !end_time_s.is_finite()
            || end_time_s < self.time_s
            || !(end_time_s - self.time_s).is_finite()
        {
            return Err(PopulationEvolutionError::InvalidInput(
                "forward population time interval",
            ));
        }
        // This also bounds the sum of native RHS counters and prevents a work
        // budget itself from overflowing the aggregate counters.
        if numerics.maximum_total_attempts > (usize::MAX - self.population.occupied_classes()) / 7 {
            return Err(PopulationEvolutionError::InvalidInput(
                "population attempt budget",
            ));
        }
        self.scratch.copy_from_slice(&self.states);
        let mut accepted_steps = 0;
        let mut rejected_steps = 0;
        let mut rhs_evaluations = 0;
        let mut advanced_classes = 0;
        let mut maximum_normalized_local_error: f64 = 0.0;
        if end_time_s > self.time_s {
            for (class_index, (class, state)) in self
                .population
                .classes()
                .iter()
                .zip(self.scratch.iter_mut())
                .enumerate()
            {
                poll(cancelled)?;
                if class.number_density_per_m3() == 0.0 {
                    // An absent class has no radial equation to advance.
                    *state = BubbleState::try_new(
                        end_time_s,
                        state.radius_m(),
                        state.wall_velocity_m_s(),
                    )?;
                    continue;
                }
                let consumed = accepted_steps + rejected_steps;
                let remaining = numerics.maximum_total_attempts - consumed;
                if remaining == 0 {
                    return Err(PopulationEvolutionError::WorkBudgetExceeded {
                        class_index,
                        consumed_attempts: consumed,
                        maximum_total_attempts: numerics.maximum_total_attempts,
                    });
                }
                let mut controls = numerics.per_class;
                controls.maximum_attempts = controls.maximum_attempts.min(remaining);
                let report = class
                    .bubble()
                    .advance_controlled(state, end_time_s, &controls, pressure, cancelled)
                    .map_err(|cause| match cause {
                        CavitationError::Cancelled => PopulationEvolutionError::Cancelled,
                        CavitationError::StepBudgetExceeded { attempts, .. }
                            if controls.maximum_attempts == remaining =>
                        {
                            PopulationEvolutionError::WorkBudgetExceeded {
                                class_index,
                                consumed_attempts: consumed + attempts,
                                maximum_total_attempts: numerics.maximum_total_attempts,
                            }
                        }
                        cause => PopulationEvolutionError::Class { class_index, cause },
                    })?;
                accepted_steps += report.accepted_steps;
                rejected_steps += report.rejected_steps;
                rhs_evaluations += report.rhs_evaluations;
                advanced_classes += 1;
                maximum_normalized_local_error =
                    maximum_normalized_local_error.max(report.maximum_normalized_local_error);
            }
        }
        poll(cancelled)?;
        let observation =
            self.observe_states(&self.scratch, end_time_s, pressure(end_time_s), cancelled)?;
        poll(cancelled)?;
        let step = PopulationEvolutionStep {
            observation,
            advanced_classes,
            accepted_steps,
            rejected_steps,
            rhs_evaluations,
            maximum_normalized_local_error,
        };
        Ok(PopulationAdvance { owner: self, step })
    }

    fn observe_states(
        &self,
        states: &[BubbleState],
        time_s: f64,
        pressure: BubblePressureSample,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PopulationObservation, PopulationEvolutionError> {
        BubblePressureSample::try_new(
            pressure.absolute_pressure_pa,
            pressure.time_derivative_pa_s,
        )?;
        let mut fraction = Sum::default();
        let mut rate = Sum::default();
        let mut acceleration = Sum::default();
        let mut number_radius = Sum::default();
        let mut maximum_radius_m: f64 = 0.0;
        let mut maximum_wall_mach_number: f64 = 0.0;
        for (class_index, (class, state)) in
            self.population.classes().iter().zip(states).enumerate()
        {
            poll(cancelled)?;
            if state.time_s() != time_s {
                return Err(PopulationEvolutionError::ClockMismatch {
                    class_index,
                    expected_time_s: time_s,
                    actual_time_s: state.time_s(),
                });
            }
            let number = class.number_density_per_m3();
            if number == 0.0 {
                continue;
            }
            let native = class
                .bubble()
                .observe(state, pressure)
                .map_err(|cause| PopulationEvolutionError::Class { class_index, cause })?;
            let weighted_volume = number * native.volume_m3;
            let weighted_rate = number * native.volume_rate_m3_s;
            let weighted_acceleration = number * native.volume_acceleration_m3_s2;
            let weighted_radius = number * state.radius_m();
            // An occupied class cannot silently disappear in the mixture.
            // Exact native zero derivatives remain valid, as does cancellation
            // between different representable class contributions.
            if ![
                weighted_volume,
                weighted_rate,
                weighted_acceleration,
                weighted_radius,
            ]
            .into_iter()
            .all(f64::is_finite)
                || weighted_volume <= 0.0
                || weighted_radius <= 0.0
                || (native.volume_rate_m3_s != 0.0 && weighted_rate == 0.0)
                || (native.volume_acceleration_m3_s2 != 0.0 && weighted_acceleration == 0.0)
            {
                return Err(PopulationEvolutionError::InvalidInput(
                    "unrepresentable weighted bubble source",
                ));
            }
            fraction.add(weighted_volume);
            rate.add(weighted_rate);
            acceleration.add(weighted_acceleration);
            number_radius.add(weighted_radius);
            maximum_radius_m = maximum_radius_m.max(state.radius_m());
            maximum_wall_mach_number = maximum_wall_mach_number.max(native.wall_mach_number);
        }
        let void_fraction = fraction.value();
        let void_fraction_rate_per_s = rate.value();
        let void_fraction_acceleration_per_s2 = acceleration.value();
        let number_radius_moment_per_m2 = number_radius.value();
        if ![
            void_fraction,
            void_fraction_rate_per_s,
            void_fraction_acceleration_per_s2,
            number_radius_moment_per_m2,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(PopulationEvolutionError::InvalidInput(
                "unrepresentable aggregate bubble source",
            ));
        }
        if void_fraction > self.limits.maximum_void_fraction {
            return Err(PopulationEvolutionError::DilutionExceeded {
                void_fraction,
                maximum_void_fraction: self.limits.maximum_void_fraction,
            });
        }
        Ok(PopulationObservation {
            time_s,
            imposed_pressure: pressure,
            occupied_classes: self.population.occupied_classes(),
            number_density_per_m3: self.population.number_density_per_m3(),
            void_fraction,
            void_fraction_rate_per_s,
            void_fraction_acceleration_per_s2,
            number_radius_moment_per_m2,
            maximum_radius_m,
            maximum_wall_mach_number,
        })
    }
}

/// A complete proposed physical interval borrowing its original owner.
/// Dropping it preserves the accepted clock and every accepted class state.
#[must_use = "inspect and commit the proposal to publish the physical interval"]
#[derive(Debug)]
pub struct PopulationAdvance<'a> {
    owner: &'a mut EvolvingBubblePopulation,
    step: PopulationEvolutionStep,
}
impl PopulationAdvance<'_> {
    pub const fn observation(&self) -> &PopulationObservation {
        &self.step.observation
    }
    pub const fn step(&self) -> &PopulationEvolutionStep {
        &self.step
    }
    pub fn states(&self) -> &[BubbleState] {
        &self.owner.scratch
    }
    pub const fn population(&self) -> &BubblePopulation {
        &self.owner.population
    }

    /// Infallible publication of this already admitted proposal. A downstream
    /// observer should validate and poll its cancellation before calling this.
    pub fn commit(self) -> PopulationEvolutionStep {
        let Self { owner, step } = self;
        core::mem::swap(&mut owner.states, &mut owner.scratch);
        owner.time_s = step.observation.time_s;
        step
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PopulationEvolutionError> {
    if cancelled() {
        Err(PopulationEvolutionError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "evolution/tests.rs"]
mod tests;
