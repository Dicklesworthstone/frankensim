//! Continuous piecewise-inertial source motion on a finite authored history.
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
    clock, finite_vector, norm, position,
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
        if !max_residual_m.is_finite() || max_residual_m < 0.0 {
            return Err(RetardedVacuumError::InvalidTolerance);
        }
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

    fn knot(&self, index: usize) -> (Dd, [Dd; 3]) {
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

    fn knot_residual(
        &self,
        index: usize,
        receiver: VacuumReception,
    ) -> Result<Dd, RetardedVacuumError> {
        let (time, location) = self.knot(index);
        let displacement = finite_vector(core::array::from_fn(|i| {
            Dd::from_f64(receiver.position_m[i]) - location[i]
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
