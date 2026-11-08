//! Continuous piecewise-inertial source or receiver motion on a finite history.
//!
//! Positions interpolate linearly in world coordinate time. Each segment's
//! material proper time advances by `dt * sqrt(1-v²/c²)`, accumulated in DD;
//! a later segment never reinterprets the earlier clock with its own gamma.
//! The shared local Lorentz owner supplies the Doppler factor at emission.
//! See UCSD Physics 110B, sections 19.5.2 and 19.6.1, for the proper-time
//! integral and retarded reception of successive inertial legs:
//! <https://courses.physics.ucsd.edu/2019/Winter/physics110b/LECTURES/110B_ALL.pdf>.
//!
//! A velocity discontinuity is an ideal instantaneous change, not a resolved
//! acceleration. At an exact interior knot the outgoing segment supplies the
//! right-hand frequency and clock derivative. The last knot uses the incoming
//! segment. Refining sampled motion is the caller's responsibility; this is
//! neither a smooth-trajectory error bound nor a source-visibility solver.

use super::{
    RetardedVacuumError, RetardedVacuumRay, UniformVacuumWorldline, VacuumReception, checked,
    clock, emission_input, finite_vector, norm, position,
};
use crate::spectral::doppler::{FrameVelocity, VACUUM_LIGHT_SPEED_M_S};
use fs_math::dd::Dd;

/// An authored event on the source's world-coordinate position history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VacuumTrajectoryKnot {
    pub coordinate_time_s: Dd,
    pub position_m: [f64; 3],
}

#[derive(Clone, Copy, Debug)]
struct Segment {
    motion: UniformVacuumWorldline,
    /// Lost preceding proper-clock increments, without hiding epoch rounding.
    prefix_clock_residual_s: Dd,
}

/// Prepared finite trajectory with one allocation and O(log N) ray queries.
///
/// All segments are admitted as strictly subluminal before publication. Knot
/// positions, coordinate clocks and derived SI velocities share DD arithmetic.
/// The proper-time origin belongs to the first knot. Queries neither allocate
/// nor extrapolate outside the closed authored history.
#[derive(Clone, Debug)]
pub struct PiecewiseVacuumWorldline {
    segments: Vec<Segment>,
    last_coordinate_time_s: Dd,
    last_position_m: [Dd; 3],
    proper_history_s: [Dd; 2],
}

impl PiecewiseVacuumWorldline {
    /// Prepare at most `max_segments` intervals from strictly ordered knots.
    /// Storage and numerical admission finish before this trajectory is usable.
    /// The actual source temperature/history remains with its original owner.
    pub fn try_new(
        knots: &[VacuumTrajectoryKnot],
        first_proper_time_s: Dd,
        max_segments: usize,
    ) -> Result<Self, RetardedVacuumError> {
        if knots.len() < 2 {
            return Err(RetardedVacuumError::InvalidTrajectory);
        }
        if knots.len() - 1 > max_segments {
            return Err(RetardedVacuumError::TrajectoryBudget);
        }
        let first_proper_time_s = clock(first_proper_time_s)?;
        let mut time = clock(knots[0].coordinate_time_s)?;
        position(knots[0].position_m)?;
        let mut location = knots[0].position_m.map(Dd::from_f64);
        let mut proper_offset = Dd::ZERO;
        let mut prefix_clock_residual_s = Dd::ZERO;
        let mut segments = Vec::new();
        segments
            .try_reserve_exact(knots.len() - 1)
            .map_err(|_| RetardedVacuumError::AllocationRefused)?;
        for knot in &knots[1..] {
            let next_time = clock(knot.coordinate_time_s)?;
            position(knot.position_m)?;
            let next_location = knot.position_m.map(Dd::from_f64);
            let duration = checked(next_time - time)?;
            if !Dd::ZERO.lt(duration) {
                return Err(RetardedVacuumError::InvalidTrajectory);
            }
            let velocity_m_s = finite_vector(core::array::from_fn(|i| {
                (next_location[i] - location[i]) / duration
            }))?;
            let frame = FrameVelocity::from_dd_velocity(velocity_m_s)?;
            let proper_duration = checked(duration * frame.inverse_gamma)?;
            let next_proper_offset = checked(proper_offset + proper_duration)?;
            if !Dd::ZERO.lt(proper_duration) || !proper_offset.lt(next_proper_offset) {
                return Err(RetardedVacuumError::NumericalRange);
            }
            segments.push(Segment {
                motion: UniformVacuumWorldline {
                    reference_coordinate_time_s: time,
                    // Keep the epoch separate until the final event is known.
                    reference_proper_time_s: proper_offset,
                    reference_position_m: location,
                    velocity_m_s,
                    proper_history_s: [proper_offset, next_proper_offset],
                    frame,
                },
                prefix_clock_residual_s,
            });
            prefix_clock_residual_s = checked(
                prefix_clock_residual_s
                    + checked(((next_proper_offset - proper_offset) - proper_duration).abs())?,
            )?;
            time = next_time;
            location = next_location;
            proper_offset = next_proper_offset;
        }
        let last_proper_time_s = checked(first_proper_time_s + proper_offset)?;
        if !first_proper_time_s.lt(last_proper_time_s) {
            return Err(RetardedVacuumError::NumericalRange);
        }
        Ok(Self {
            segments,
            last_coordinate_time_s: time,
            last_position_m: location,
            proper_history_s: [first_proper_time_s, last_proper_time_s],
        })
    }

    /// Number of admitted inertial intervals.
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    /// Closed authored source support in world-coordinate seconds.
    pub fn coordinate_history_s(&self) -> [Dd; 2] {
        [
            self.segments[0].motion.reference_coordinate_time_s,
            self.last_coordinate_time_s,
        ]
    }

    /// The same closed source support in accumulated material proper seconds.
    pub const fn proper_history_s(&self) -> [Dd; 2] {
        self.proper_history_s
    }

    /// Locate a moving detector event on its accumulated proper clock.
    ///
    /// A binary search selects the proper-time interval before inverting that
    /// interval's own clock rate. All preceding time-dilation increments remain
    /// intact. Interior knots use the outgoing velocity and the last endpoint
    /// the incoming velocity, just as source observations do. Queries allocate
    /// nothing, never extrapolate and retain DD receiver positions.
    pub fn reception_at_proper_time(
        &self,
        proper_time_s: Dd,
        max_residual_m: f64,
    ) -> Result<VacuumReception, RetardedVacuumError> {
        if !max_residual_m.is_finite() || max_residual_m < 0.0 {
            return Err(RetardedVacuumError::InvalidTolerance);
        }
        let proper_time_s = clock(proper_time_s)?;
        if proper_time_s.lt(self.proper_history_s[0]) || self.proper_history_s[1].lt(proper_time_s)
        {
            return Err(RetardedVacuumError::OutsideHistory);
        }
        let count = self.segments.len();
        let mut lower = 0;
        let mut upper = count;
        while upper - lower > 1 {
            let middle = lower + (upper - lower) / 2;
            if proper_time_s.lt(self.proper_knot_time(middle)?) {
                upper = middle;
            } else {
                lower = middle;
            }
        }
        let start = self.proper_knot_time(lower)?;
        if !start.lt(self.proper_knot_time(lower + 1)?)
            || (proper_time_s == start && lower > 0 && self.proper_knot_time(lower - 1)? == start)
        {
            // A collapsed receiver epoch cannot choose between distinct knot
            // velocities, even if a permissive distance tolerance was supplied.
            return Err(RetardedVacuumError::NumericalRange);
        }
        let exact_knot = if proper_time_s == start {
            Some(lower)
        } else if proper_time_s == self.proper_history_s[1] {
            Some(count)
        } else {
            None
        };
        let offset = if let Some(knot) = exact_knot {
            self.proper_knot_offset(knot)
        } else {
            checked(proper_time_s - self.proper_history_s[0])?
        };
        let segment = &self.segments[lower];
        let receiver = if let Some(knot) = exact_knot {
            segment
                .motion
                .reception_event(self.knot(knot).0, offset, max_residual_m)?
        } else {
            segment
                .motion
                .reception_at_proper_time(offset, max_residual_m)?
        };
        self.publish_reception(segment, receiver, proper_time_s, offset, max_residual_m)
    }

    fn publish_reception(
        &self,
        segment: &Segment,
        mut receiver: VacuumReception,
        proper_time_s: Dd,
        offset: Dd,
        max_residual_m: f64,
    ) -> Result<VacuumReception, RetardedVacuumError> {
        let epoch_loss = checked(
            checked(((proper_time_s - self.proper_history_s[0]) - offset).abs())?
                + checked((self.proper_history_s[0] + offset - proper_time_s).abs())?,
        )?;
        receiver.proper_clock_residual_s = checked(
            receiver.proper_clock_residual_s + segment.prefix_clock_residual_s + epoch_loss,
        )?;
        receiver.proper_time_s = Some(proper_time_s);
        receiver.admitted(max_residual_m)
    }

    /// Receive a fixed emission event on this finite detector trajectory.
    /// The future-cone residual decreases monotonically over subluminal motion,
    /// so a binary search selects one segment before the analytic root. Exact
    /// interior knots use the outgoing velocity; the final knot uses incoming.
    pub fn receive_emission(
        &self,
        emission_coordinate_time_s: Dd,
        emission_position_m: [Dd; 3],
        max_residual_m: f64,
    ) -> Result<VacuumReception, RetardedVacuumError> {
        let (time, location) = emission_input(
            emission_coordinate_time_s,
            emission_position_m,
            max_residual_m,
        )?;
        let residual_at = |index| {
            let (arrival, position) = self.knot(index);
            let displacement = finite_vector(core::array::from_fn(|i| position[i] - location[i]))?;
            checked(
                norm(displacement)?
                    - checked(Dd::from_f64(VACUUM_LIGHT_SPEED_M_S) * checked(arrival - time)?)?,
            )
        };
        let count = self.segments.len();
        let first = residual_at(0)?;
        let last = residual_at(count)?;
        if first.lt(Dd::ZERO) || Dd::ZERO.lt(last) {
            return Err(RetardedVacuumError::OutsideHistory);
        }
        let mut lower = 0;
        let mut upper = count;
        let mut exact = if first == Dd::ZERO {
            Some(0)
        } else if last == Dd::ZERO {
            lower = count - 1;
            Some(count)
        } else {
            None
        };
        while exact.is_none() && upper - lower > 1 {
            let middle = lower + (upper - lower) / 2;
            let residual = residual_at(middle)?;
            if !residual.lt(Dd::ZERO) {
                lower = middle;
                if residual == Dd::ZERO {
                    exact = Some(middle);
                }
            } else {
                upper = middle;
            }
        }
        let segment = &self.segments[lower];
        let receiver = if let Some(index) = exact {
            segment
                .motion
                .receive_at(time, location, self.knot(index).0, max_residual_m)?
        } else {
            segment
                .motion
                .receive_emission(time, location, max_residual_m)?
        };
        let offset = receiver
            .proper_time_s
            .ok_or(RetardedVacuumError::NumericalRange)?;
        let proper = checked(self.proper_history_s[0] + offset)?;
        self.publish_reception(segment, receiver, proper, offset, max_residual_m)
    }

    /// Select the unique causal segment, then use the original analytic solve.
    ///
    /// At a fixed reception position, `t + |x_receiver-x_source(t)|/c`
    /// strictly increases along every admitted subluminal segment. Knot
    /// residuals therefore locate the past-light-cone intersection with a
    /// binary search; there is no iteration count depending on source speed.
    /// `max_residual_m` retains the original geometry/clock admission meaning.
    pub fn observe(
        &self,
        receiver: VacuumReception,
        max_residual_m: f64,
    ) -> Result<RetardedVacuumRay, RetardedVacuumError> {
        let receiver = receiver.admitted(max_residual_m)?;
        let count = self.segments.len();
        let first = self.knot_residual(0, receiver)?;
        let last = self.knot_residual(count, receiver)?;
        if Dd::ZERO.lt(first) || last.lt(Dd::ZERO) {
            return Err(RetardedVacuumError::OutsideHistory);
        }
        if first == Dd::ZERO {
            return self.observe_segment(0, Some(0), receiver, max_residual_m);
        }
        if last == Dd::ZERO {
            return self.observe_segment(count - 1, Some(count), receiver, max_residual_m);
        }
        let mut lower = 0;
        let mut upper = count;
        while upper - lower > 1 {
            let middle = lower + (upper - lower) / 2;
            let residual = self.knot_residual(middle, receiver)?;
            if residual == Dd::ZERO {
                return self.observe_segment(middle, Some(middle), receiver, max_residual_m);
            }
            if residual.lt(Dd::ZERO) {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        self.observe_segment(lower, None, receiver, max_residual_m)
    }

    pub(super) fn knot(&self, index: usize) -> (Dd, [Dd; 3]) {
        if index == self.segments.len() {
            (self.last_coordinate_time_s, self.last_position_m)
        } else {
            let motion = &self.segments[index].motion;
            (
                motion.reference_coordinate_time_s,
                motion.reference_position_m,
            )
        }
    }

    fn proper_knot_offset(&self, index: usize) -> Dd {
        if index == self.segments.len() {
            self.segments[index - 1].motion.proper_history_s[1]
        } else {
            self.segments[index].motion.reference_proper_time_s
        }
    }

    pub(super) fn proper_knot_time(&self, index: usize) -> Result<Dd, RetardedVacuumError> {
        checked(self.proper_history_s[0] + self.proper_knot_offset(index))
    }

    fn knot_residual(
        &self,
        index: usize,
        receiver: VacuumReception,
    ) -> Result<Dd, RetardedVacuumError> {
        let (time, location) = self.knot(index);
        let displacement = finite_vector(core::array::from_fn(|i| {
            receiver.position_m[i] - location[i]
        }))?;
        checked(
            norm(displacement)?
                - checked(
                    Dd::from_f64(VACUUM_LIGHT_SPEED_M_S)
                        * checked(receiver.coordinate_time_s - time)?,
                )?,
        )
    }

    fn observe_segment(
        &self,
        segment_index: usize,
        exact_knot: Option<usize>,
        receiver: VacuumReception,
        max_residual_m: f64,
    ) -> Result<RetardedVacuumRay, RetardedVacuumError> {
        let segment = &self.segments[segment_index];
        let mut ray = if let Some(knot) = exact_knot {
            segment
                .motion
                .observe_emission(receiver, self.knot(knot).0, max_residual_m)?
        } else {
            segment.motion.observe(receiver, max_residual_m)?
        };
        let proper_offset = ray.emission_proper_time_s;
        let proper_time = checked(self.proper_history_s[0] + proper_offset)?;
        let epoch_loss = checked(((proper_time - self.proper_history_s[0]) - proper_offset).abs())?;
        let clock_residual_s = checked(
            segment.prefix_clock_residual_s
                + Dd::from_f64(ray.proper_clock_residual_s)
                + epoch_loss,
        )?;
        let error_m = checked(
            Dd::from_f64(VACUUM_LIGHT_SPEED_M_S)
                * (clock_residual_s / segment.motion.frame.inverse_gamma),
        )?;
        if Dd::from_f64(max_residual_m).lt(error_m) {
            return Err(RetardedVacuumError::ClockResolution {
                equivalent_distance_error_m: error_m.to_f64(),
                tolerance_m: max_residual_m,
            });
        }
        ray.emission_proper_time_s = proper_time;
        ray.proper_clock_residual_s = clock_residual_s.to_f64();
        Ok(ray)
    }
}

#[cfg(test)]
mod tests;
