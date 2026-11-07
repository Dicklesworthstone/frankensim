//! Inertial material worldlines intersected with a receiver's past light cone.
//!
//! Positions, coordinate clocks and velocities share one inertial world frame.
//! The source clock is proper time along its material worldline; it is not the
//! world coordinate clock. DD clock pairs retain short flight times at large
//! epochs. The published emission event is checked against `|x_o-x_s|=c(t_o-t_s)`
//! before its source history can be sampled. This is the retarded-event geometry
//! used in, for example, Fitzpatrick's relativistic moving-source derivation:
//! <https://farside.ph.utexas.edu/teaching/jk1/Electromagnetism/node150.html>.
//!
//! The uniform owner describes one material point. [`trajectory`] composes
//! finite nonuniform histories from continuous inertial segments with their
//! accumulated proper clocks. Neither owner supplies a surface/visibility
//! intersection, apparent-area map, smooth accelerated-worldline integration,
//! dielectric path or curved-spacetime solver. Its numerical residual is an
//! admission check, not an interval certificate for uncertain input geometry.

use super::{DopplerError, FrameVelocity, VACUUM_LIGHT_SPEED_M_S, VacuumDoppler, dd_dot, unit};
use fs_math::dd::Dd;

/// Finite nonuniform motion represented by continuous inertial segments.
pub mod trajectory;
pub use trajectory::{PiecewiseVacuumWorldline, VacuumTrajectoryKnot};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RetardedVacuumError {
    InvalidClock,
    InvalidPosition,
    InvalidHistory,
    /// At least two finite, strictly time-ordered position knots are required.
    InvalidTrajectory,
    /// The authored trajectory exceeds its explicit segment-storage cap.
    TrajectoryBudget,
    AllocationRefused,
    InvalidTolerance,
    /// The sole light-cone intersection has zero flight time and no ray.
    Coincident,
    /// The unique past emission lies outside the authored proper-time history.
    OutsideHistory,
    Residual {
        residual_m: f64,
        tolerance_m: f64,
    },
    /// Lost material-clock offset, converted back to coordinate seconds * c.
    ClockResolution {
        equivalent_distance_error_m: f64,
        tolerance_m: f64,
    },
    Doppler(DopplerError),
    NumericalRange,
}

impl core::fmt::Display for RetardedVacuumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "retarded vacuum observation: {self:?}")
    }
}
impl core::error::Error for RetardedVacuumError {}
impl From<DopplerError> for RetardedVacuumError {
    fn from(value: DopplerError) -> Self {
        Self::Doppler(value)
    }
}

/// A receiver event and its instantaneous SI velocity, all in the world frame.
/// The receiver need not move uniformly away from this event.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VacuumReception {
    coordinate_time_s: Dd,
    position_m: [f64; 3],
    frame: FrameVelocity,
}

impl VacuumReception {
    pub fn try_new(
        coordinate_time_s: Dd,
        position_m: [f64; 3],
        velocity_m_s: [f64; 3],
    ) -> Result<Self, RetardedVacuumError> {
        position(position_m)?;
        Ok(Self {
            coordinate_time_s: clock(coordinate_time_s)?,
            position_m,
            frame: FrameVelocity::try_new(velocity_m_s)?,
        })
    }
}

/// An affine source trajectory with an explicitly bounded material history.
/// The paired reference clocks identify the SAME event on this worldline.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformVacuumWorldline {
    reference_coordinate_time_s: Dd,
    reference_proper_time_s: Dd,
    reference_position_m: [Dd; 3],
    velocity_m_s: [Dd; 3],
    proper_history_s: [Dd; 2],
    frame: FrameVelocity,
}

impl UniformVacuumWorldline {
    pub fn try_new(
        reference_coordinate_time_s: Dd,
        reference_proper_time_s: Dd,
        reference_position_m: [f64; 3],
        velocity_m_s: [f64; 3],
        proper_history_s: [Dd; 2],
    ) -> Result<Self, RetardedVacuumError> {
        position(reference_position_m)?;
        let proper_history_s = [clock(proper_history_s[0])?, clock(proper_history_s[1])?];
        if proper_history_s[1].lt(proper_history_s[0]) {
            return Err(RetardedVacuumError::InvalidHistory);
        }
        Ok(Self {
            reference_coordinate_time_s: clock(reference_coordinate_time_s)?,
            reference_proper_time_s: clock(reference_proper_time_s)?,
            reference_position_m: reference_position_m.map(Dd::from_f64),
            velocity_m_s: velocity_m_s.map(Dd::from_f64),
            proper_history_s,
            frame: FrameVelocity::try_new(velocity_m_s)?,
        })
    }

    /// Closed support of the actual source history in material proper seconds.
    pub const fn proper_history_s(&self) -> [Dd; 2] {
        self.proper_history_s
    }

    /// Unique positive-flight-time null ray, with fixed work and no allocation.
    ///
    /// Let q=x_o-x_s(t_o), beta=v_s/c and ell=c(t_o-t_s). The null equation is
    /// `(1-beta²) ell² - 2(q.beta) ell - q² = 0`. Scaling q before evaluating
    /// the discriminant avoids range² overflow; rationalizing the negative-b
    /// branch avoids cancellation for a source receding along the ray.
    ///
    /// `max_residual_m` is a finite nonnegative absolute distance tolerance.
    /// The residual is recomputed from the PUBLISHED DD emission clock and
    /// position, so an unresolved delay at a very large epoch is refused. The
    /// same tolerance separately bounds c times the coordinate-clock equivalent
    /// of precision lost when publishing the material proper-clock offset.
    pub fn observe(
        &self,
        receiver: VacuumReception,
        max_residual_m: f64,
    ) -> Result<RetardedVacuumRay, RetardedVacuumError> {
        if !max_residual_m.is_finite() || max_residual_m < 0.0 {
            return Err(RetardedVacuumError::InvalidTolerance);
        }
        let offset = checked(receiver.coordinate_time_s - self.reference_coordinate_time_s)?;
        let q = finite_vector(core::array::from_fn(|i| {
            (Dd::from_f64(receiver.position_m[i]) - self.reference_position_m[i])
                - self.velocity_m_s[i] * offset
        }))?;
        let scale = q.iter().fold(0.0_f64, |a, x| a.max(x.hi.abs()));
        if scale == 0.0 {
            return Err(RetardedVacuumError::Coincident);
        }
        let scale = Dd::from_f64(scale);
        let normalized = q.map(|x| x / scale);
        let squared = dd_dot(normalized, normalized);
        let b = dd_dot(normalized, self.frame.beta);
        let a = self.frame.inverse_gamma * self.frame.inverse_gamma;
        let discriminant = (b * b + a * squared).sqrt();
        let length = if b.lt(Dd::ZERO) {
            scale * (squared / (discriminant - b))
        } else {
            scale * ((discriminant + b) / a)
        };
        let c = Dd::from_f64(VACUUM_LIGHT_SPEED_M_S);
        let emission_time_s = checked(receiver.coordinate_time_s - checked(length)? / c)?;

        self.observe_emission(receiver, emission_time_s, max_residual_m)
    }

    /// Common publication checks for an analytic root or an exact knot event.
    fn observe_emission(
        &self,
        receiver: VacuumReception,
        emission_time_s: Dd,
        max_residual_m: f64,
    ) -> Result<RetardedVacuumRay, RetardedVacuumError> {
        let c = Dd::from_f64(VACUUM_LIGHT_SPEED_M_S);
        // Reconstruct every dependent quantity from the represented event,
        // not from an unpublishable intermediate difference of large clocks.
        let source_offset = checked(emission_time_s - self.reference_coordinate_time_s)?;
        let proper_offset = checked(self.frame.inverse_gamma * source_offset)?;
        let emission_proper_time_s = checked(self.reference_proper_time_s + proper_offset)?;
        let proper_clock_residual_s = checked(
            ((emission_proper_time_s - self.reference_proper_time_s) - proper_offset).abs(),
        )?;
        let clock_error_m = checked(c * (proper_clock_residual_s / self.frame.inverse_gamma))?;
        if Dd::from_f64(max_residual_m).lt(clock_error_m) {
            return Err(RetardedVacuumError::ClockResolution {
                equivalent_distance_error_m: clock_error_m.to_f64(),
                tolerance_m: max_residual_m,
            });
        }
        if emission_proper_time_s.lt(self.proper_history_s[0])
            || self.proper_history_s[1].lt(emission_proper_time_s)
        {
            return Err(RetardedVacuumError::OutsideHistory);
        }
        let emission_position_m = finite_vector(core::array::from_fn(|i| {
            self.reference_position_m[i] + self.velocity_m_s[i] * source_offset
        }))?;
        let displacement = finite_vector(core::array::from_fn(|i| {
            Dd::from_f64(receiver.position_m[i]) - emission_position_m[i]
        }))?;
        let range = norm(displacement)?;
        let flight_time_s = checked(receiver.coordinate_time_s - emission_time_s)?;
        if flight_time_s == Dd::ZERO && range == Dd::ZERO {
            return Err(RetardedVacuumError::Coincident);
        }
        if !Dd::ZERO.lt(flight_time_s) || !Dd::ZERO.lt(range) {
            return Err(RetardedVacuumError::NumericalRange);
        }
        let residual = checked((range - c * flight_time_s).abs())?;
        if Dd::from_f64(max_residual_m).lt(residual) {
            return Err(RetardedVacuumError::Residual {
                residual_m: residual.to_f64(),
                tolerance_m: max_residual_m,
            });
        }
        let direction = unit(displacement.map(|x| (x / range).to_f64()))?;
        let doppler = VacuumDoppler::with_frames(direction, self.frame, receiver.frame)?;
        let coordinate_time_rate = checked(
            Dd::from_f64(doppler.frequency_ratio()) * receiver.frame.inverse_gamma
                / self.frame.inverse_gamma,
        )?
        .to_f64();
        if coordinate_time_rate <= 0.0 {
            return Err(RetardedVacuumError::NumericalRange);
        }
        Ok(RetardedVacuumRay {
            emission_time_s,
            emission_proper_time_s,
            emission_position_m,
            reception_time_s: receiver.coordinate_time_s,
            flight_time_s,
            range_m: range.to_f64(),
            light_cone_residual_m: residual.to_f64(),
            proper_clock_residual_s: proper_clock_residual_s.to_f64(),
            world_direction: direction,
            coordinate_time_rate,
            doppler,
        })
    }
}

/// One admitted physical ray, retaining its geometry, clocks and Lorentz view.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RetardedVacuumRay {
    emission_time_s: Dd,
    emission_proper_time_s: Dd,
    emission_position_m: [Dd; 3],
    reception_time_s: Dd,
    flight_time_s: Dd,
    range_m: f64,
    light_cone_residual_m: f64,
    proper_clock_residual_s: f64,
    world_direction: [f64; 3],
    coordinate_time_rate: f64,
    doppler: VacuumDoppler,
}

impl RetardedVacuumRay {
    pub const fn emission_coordinate_time_s(&self) -> Dd {
        self.emission_time_s
    }

    /// The original material history must be sampled on this proper clock.
    pub const fn emission_proper_time_s(&self) -> Dd {
        self.emission_proper_time_s
    }

    /// DD coordinates retain small displacements from large world origins.
    pub const fn emission_position_m(&self) -> [Dd; 3] {
        self.emission_position_m
    }

    pub const fn reception_coordinate_time_s(&self) -> Dd {
        self.reception_time_s
    }

    pub const fn flight_time_s(&self) -> Dd {
        self.flight_time_s
    }

    pub const fn range_m(&self) -> f64 {
        self.range_m
    }

    pub const fn light_cone_residual_m(&self) -> f64 {
        self.light_cone_residual_m
    }

    /// Precision lost when adding the proper-time offset to its reference.
    pub const fn proper_clock_residual_s(&self) -> f64 {
        self.proper_clock_residual_s
    }

    /// Points from the actual emission position toward the receiver event.
    pub const fn world_propagation_direction(&self) -> [f64; 3] {
        self.world_direction
    }

    pub const fn doppler(&self) -> VacuumDoppler {
        self.doppler
    }

    /// dt_emission / dt_reception for the declared receiver event velocity.
    pub const fn coordinate_time_rate(&self) -> f64 {
        self.coordinate_time_rate
    }

    /// d(tau_source) / d(tau_receiver), equal to nu_receiver/nu_source.
    /// Do not apply this again after sampling the retarded source history.
    pub const fn proper_time_rate(&self) -> f64 {
        self.doppler.frequency_ratio()
    }
}

fn clock(value: Dd) -> Result<Dd, RetardedVacuumError> {
    if !value.hi.is_finite() || !value.lo.is_finite() {
        return Err(RetardedVacuumError::InvalidClock);
    }
    checked(Dd::from_pair(value.hi, value.lo)).map_err(|_| RetardedVacuumError::InvalidClock)
}

fn position(value: [f64; 3]) -> Result<(), RetardedVacuumError> {
    if value.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(RetardedVacuumError::InvalidPosition)
    }
}

fn checked(value: Dd) -> Result<Dd, RetardedVacuumError> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err(RetardedVacuumError::NumericalRange)
    }
}

fn finite_vector(value: [Dd; 3]) -> Result<[Dd; 3], RetardedVacuumError> {
    for x in value {
        checked(x)?;
    }
    Ok(value)
}

fn norm(value: [Dd; 3]) -> Result<Dd, RetardedVacuumError> {
    let scale = value.iter().fold(0.0_f64, |a, x| a.max(x.hi.abs()));
    if scale == 0.0 {
        return Ok(Dd::ZERO);
    }
    let scale = Dd::from_f64(scale);
    let normalized = value.map(|x| x / scale);
    checked(scale * dd_dot(normalized, normalized).sqrt())
}

#[cfg(test)]
mod tests;
