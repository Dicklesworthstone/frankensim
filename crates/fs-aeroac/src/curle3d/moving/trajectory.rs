//! Retarded observation of finite, smoothly accelerating source trajectories.
//!
//! Each accepted position/velocity/acceleration pair defines one C2 quintic
//! Hermite segment. Position and BOTH derivatives are evaluated from that same
//! polynomial; the interpolant is an observation model, not a second mechanics
//! integrator or a bound on the original solver's error.
//!
//! Preparation uses outward `fs_ivl` Bernstein hulls to admit speed and receiver
//! clearance over every complete segment. Inconclusive hulls subdivide under
//! explicit work/depth caps. The admitted speed is strictly subsonic, so
//! `d(tau + R(tau)/c)/d tau = 1 - n.v/c > 0`: every covered receiver time has
//! one retarded root. Safeguarded Newton/bisection then queries the original
//! moving compact-pressure owner. No trajectory, strength, or tail extrapolates.

use super::{
    MovingCompactObserver, PreparedMovingCompactSample, SourceKinematics, dot, norm, residual,
};
use crate::AeroacError;
use fs_ivl::Interval;
use fs_math::dd::Dd;

/// Shared, observer-independent C2 path preparation.
pub mod path;
/// Causal observation of two independently moving, accepted paths.
pub mod receiver;
/// Exact ideal-wall paths on the original finite source trajectory and clock.
pub mod reflection;

pub use path::PreparedKinematicTrajectory;
pub use receiver::{
    MovingReceiverBudget, MovingReceiverConfig, PreparedMovingReceiverSample,
    PreparedMovingReceiverSource,
};
pub use reflection::{
    IdealReflectionObserver, PreparedReflectedSourceTrajectory, ReflectedTrajectoryBudget,
};

const MAX_BOUND_DEPTH: usize = 32;
const MAX_ROOT_ITERATIONS: usize = 256;

/// One accepted physical source event in the medium's fixed coordinate frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TimedSourceKinematics {
    /// Physical emission time [s].
    pub source_time_s: f64,
    /// Position, velocity and acceleration at exactly that time.
    pub kinematics: SourceKinematics,
}

/// Retained storage, preparation work and per-query refinement bounds.
#[derive(Clone, Copy, Debug)]
pub struct TrajectoryBudget {
    /// Maximum original quintic segments; the number of accepted knots minus one.
    pub maximum_segments: usize,
    /// Logical owner plus segment payload bytes. Allocator metadata/rounding and
    /// the fixed preparation stack (at most 33 bound nodes) are excluded.
    pub maximum_payload_bytes: usize,
    /// Total Bernstein-bound nodes visited across all source segments.
    pub maximum_bound_nodes: usize,
    /// Maximum binary subdivision depth, inclusive, in 0..=32.
    pub maximum_subdivision_depth: usize,
    /// Maximum interior root refinements per query, in 1..=256. Endpoint checks
    /// and final pressure preparation add bounded constant work.
    pub maximum_root_iterations: usize,
}

#[derive(Clone, Copy, Debug)]
struct QuinticSegment {
    first: TimedSourceKinematics,
    last: TimedSourceKinematics,
    duration: Dd,
    // Position controls are offsets from first.kinematics.position_m.
    position: [[Dd; 3]; 6],
    velocity: [[Dd; 3]; 5],
    acceleration: [[Dd; 3]; 4],
    first_arrival: f64,
    last_arrival: f64,
}

#[derive(Clone, Copy)]
struct BoundNode {
    position: [[Interval; 3]; 6],
    velocity: [[Interval; 3]; 5],
    depth: usize,
}

/// A finite source path admitted for ONE original microphone and medium.
///
/// Constructor-only subdivision changes no source knots. Runtime queries use
/// the original quintics, allocate nothing and depend only on the requested
/// timestamp. Bounds concern the represented interpolation; material physics,
/// finite-source compactness, interpolation convergence and bandwidth remain
/// caller-owned. Force vectors and derivatives must still use the fixed frame.
#[derive(Debug)]
pub struct PreparedSourceTrajectory {
    observer: MovingCompactObserver,
    segments: Vec<QuinticSegment>,
    maximum_source_time_rate: f64,
    maximum_root_iterations: usize,
}

impl PreparedSourceTrajectory {
    /// Compile accepted kinematics and admit the complete finite trajectory.
    /// Refusal or cancellation publishes no partially prepared owner.
    pub fn new(
        observer: MovingCompactObserver,
        samples: &[TimedSourceKinematics],
        budget: TrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let count = samples
            .len()
            .checked_sub(1)
            .ok_or_else(|| invalid("trajectory needs two knots"))?;
        if count == 0
            || count > budget.maximum_segments
            || budget.maximum_bound_nodes < count
            || budget.maximum_subdivision_depth > MAX_BOUND_DEPTH
            || !(1..=MAX_ROOT_ITERATIONS).contains(&budget.maximum_root_iterations)
        {
            return Err(invalid("finite trajectory segment, bound and root budgets"));
        }
        let bytes = count
            .checked_mul(size_of::<QuinticSegment>())
            .and_then(|n| n.checked_add(size_of::<Self>()))
            .ok_or_else(|| invalid("trajectory retained payload overflow"))?;
        if bytes > budget.maximum_payload_bytes {
            return Err(invalid("trajectory exceeds retained payload budget"));
        }
        let first = observer.prepare_emission(samples[0].source_time_s, samples[0].kinematics)?;
        let mut arrival = first.receiver_time_s();
        let mut segments = Vec::new();
        segments
            .try_reserve_exact(count)
            .map_err(|_| invalid("cannot reserve trajectory segments"))?;
        let mut visited = 0_usize;
        let mut maximum_mach = 0.0_f64;
        for pair in samples.windows(2) {
            poll(cancelled)?;
            if !pair[1].source_time_s.is_finite() || pair[1].source_time_s <= pair[0].source_time_s
            {
                return Err(invalid(
                    "trajectory source times must be finite and strictly increasing",
                ));
            }
            let end = observer.prepare_emission(pair[1].source_time_s, pair[1].kinematics)?;
            if end.receiver_time_s() <= arrival {
                return Err(invalid(
                    "trajectory knot arrivals are not representably increasing",
                ));
            }
            let (segment, bounds, stationary) =
                QuinticSegment::compile(pair[0], pair[1], arrival, end.receiver_time_s())?;
            let admitted_mach = admit_segment(
                &observer,
                &segment,
                bounds,
                stationary,
                budget,
                &mut visited,
                cancelled,
            )?;
            maximum_mach = maximum_mach.max(admitted_mach);
            arrival = end.receiver_time_s();
            segments.push(segment);
        }
        let maximum_source_time_rate = if maximum_mach == 0.0 {
            1.0
        } else {
            (Interval::point(1.0) / (Interval::point(1.0) - Interval::point(maximum_mach))).hi()
        };
        if !maximum_source_time_rate.is_finite() || maximum_source_time_rate <= 0.0 {
            return Err(invalid("trajectory Doppler bound is unrepresentable"));
        }
        poll(cancelled)?;
        Ok(Self {
            observer,
            segments,
            maximum_source_time_rate,
            maximum_root_iterations: budget.maximum_root_iterations,
        })
    }

    /// Original fixed microphone, medium snapshot and admission limits.
    pub const fn observer(&self) -> &MovingCompactObserver {
        &self.observer
    }
    /// Complete retained source-time domain [s], inclusive.
    pub fn source_time_interval_s(&self) -> [f64; 2] {
        [
            self.segments[0].first.source_time_s,
            self.segments[self.segments.len() - 1].last.source_time_s,
        ]
    }
    /// Rounded arrivals of the first and last retained events [s], inclusive.
    pub fn receiver_time_interval_s(&self) -> [f64; 2] {
        [
            self.segments[0].first_arrival,
            self.segments[self.segments.len() - 1].last_arrival,
        ]
    }
    /// Conservative whole-trajectory upper bound on d tau/d t.
    /// Geometric amplitude modulation still needs its own bandwidth declaration.
    pub const fn maximum_source_time_rate(&self) -> f64 {
        self.maximum_source_time_rate
    }
    /// Per-query interior refinement cap retained from construction.
    pub const fn maximum_root_iterations(&self) -> usize {
        self.maximum_root_iterations
    }
    /// Number of original interpolation spans; bound subdivision adds no spans.
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    /// Read the same C2 position polynomial and its first two derivatives.
    /// Original knots return the accepted kinematics verbatim.
    pub fn kinematics_at(&self, source_time: f64) -> Result<SourceKinematics, AeroacError> {
        let [first, last] = self.source_time_interval_s();
        if !source_time.is_finite() || source_time < first || source_time > last {
            return Err(invalid("kinematics query is outside retained trajectory"));
        }
        let index = self
            .segments
            .partition_point(|segment| segment.last.source_time_s < source_time);
        self.segments[index].evaluate(source_time, true)
    }

    /// Invert the monotone arrival map and prepare the original pressure law.
    ///
    /// The requested receiver timestamp is retained exactly. Admission measures
    /// `|t - tau - R/c|` after rounding tau to its PUBLIC f64 timestamp and
    /// re-evaluating kinematics there. An unrepresentable root therefore refuses
    /// instead of returning strengths and geometry on different clocks.
    pub fn prepare_receiver(
        &self,
        receiver_time: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingCompactSample, AeroacError> {
        poll(cancelled)?;
        let [first, last] = self.receiver_time_interval_s();
        if !receiver_time.is_finite()
            || receiver_time < first
            || receiver_time > last
            || !maximum_time_residual_s.is_finite()
            || maximum_time_residual_s < 0.0
        {
            return Err(invalid(
                "receiver query is outside finite trajectory or residual domain",
            ));
        }
        let index = self
            .segments
            .partition_point(|segment| segment.last_arrival < receiver_time);
        let segment = &self.segments[index];
        self.prepare_segment(
            &self.observer,
            segment,
            receiver_time,
            [segment.first_arrival, segment.last_arrival],
            maximum_time_residual_s,
            cancelled,
        )
    }

    // A reflected path uses the original polynomial with a virtual receiver.
    // Only its two endpoint arrivals change; root refinement and pressure
    // evaluation have one owner for direct and reflected propagation.
    fn prepare_segment(
        &self,
        observer: &MovingCompactObserver,
        segment: &QuinticSegment,
        receiver_time: f64,
        [first_arrival, last_arrival]: [f64; 2],
        maximum_time_residual_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingCompactSample, AeroacError> {
        if receiver_time == first_arrival {
            return Self::finish(
                observer,
                segment,
                segment.first.source_time_s,
                receiver_time,
                maximum_time_residual_s,
                cancelled,
            );
        }
        if receiver_time == last_arrival {
            return Self::finish(
                observer,
                segment,
                segment.last.source_time_s,
                receiver_time,
                maximum_time_residual_s,
                cancelled,
            );
        }
        let mut lo = segment.first.source_time_s;
        let mut hi = segment.last.source_time_s;
        let mut flo = Self::arrival_residual(observer, segment, lo, receiver_time)?.0;
        let mut fhi = Self::arrival_residual(observer, segment, hi, receiver_time)?.0;
        if !flo.lt(Dd::ZERO) || !Dd::ZERO.lt(fhi) {
            return Err(invalid(
                "rounded trajectory arrivals do not bracket the receiver query",
            ));
        }
        let fraction = (Dd::from_f64(receiver_time) - Dd::from_f64(first_arrival))
            / (Dd::from_f64(last_arrival) - Dd::from_f64(first_arrival));
        let mut time = (Dd::from_f64(lo) + segment.duration * fraction).to_f64();
        if time <= lo || time >= hi {
            time = lo.midpoint(hi);
        }
        let mut best_error = flo.abs().to_f64().min(fhi.abs().to_f64());
        for iteration in 0..self.maximum_root_iterations {
            poll(cancelled)?;
            if time <= lo || time >= hi {
                let endpoint = if flo.abs().lt(fhi.abs()) { lo } else { hi };
                return Self::finish(
                    observer,
                    segment,
                    endpoint,
                    receiver_time,
                    maximum_time_residual_s,
                    cancelled,
                );
            }
            let (value, derivative) =
                Self::arrival_residual(observer, segment, time, receiver_time)?;
            let error = value.abs().to_f64();
            best_error = best_error.min(error);
            if error <= maximum_time_residual_s {
                return Self::finish(
                    observer,
                    segment,
                    time,
                    receiver_time,
                    maximum_time_residual_s,
                    cancelled,
                );
            }
            if value.lt(Dd::ZERO) {
                lo = time;
                flo = value;
            } else {
                hi = time;
                fhi = value;
            }
            let newton = (Dd::from_f64(time) - value / Dd::from_f64(derivative)).to_f64();
            // An occasional forced bisection bounds stalled boundary-hugging
            // Newton sequences while retaining rapid local convergence.
            time = if iteration % 4 != 3 && newton.is_finite() && lo < newton && newton < hi {
                newton
            } else {
                lo.midpoint(hi)
            };
        }
        Err(AeroacError::NotConverged {
            what: "trajectory retarded-time iteration budget",
            residual: best_error,
        })
    }

    fn arrival_residual(
        observer: &MovingCompactObserver,
        segment: &QuinticSegment,
        time: f64,
        receiver: f64,
    ) -> Result<(Dd, f64), AeroacError> {
        let state = segment.evaluate(time, false)?;
        let displacement = core::array::from_fn(|i| observer.observer_m[i] - state.position_m[i]);
        let radius = norm(displacement);
        let c = observer.medium.sound_speed_m_s();
        let delay = radius / c;
        let mach = state.velocity_m_s.map(|v| v / c);
        observer.admit_speed(norm(mach))?;
        if !radius.is_finite()
            || radius < observer.limits.minimum_distance_m
            || !delay.is_finite()
            || delay <= 0.0
        {
            return Err(invalid("trajectory retarded geometry is unrepresentable"));
        }
        let derivative = 1.0 - dot(displacement.map(|v| v / radius), mach);
        let value = Dd::from_f64(time) - Dd::from_f64(receiver) + Dd::from_f64(delay);
        if !finite_dd(value) || !derivative.is_finite() || derivative <= 0.0 {
            return Err(invalid(
                "trajectory retarded residual or derivative is unrepresentable",
            ));
        }
        Ok((value, derivative))
    }

    fn finish(
        observer: &MovingCompactObserver,
        segment: &QuinticSegment,
        time: f64,
        receiver: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingCompactSample, AeroacError> {
        let mut prepared = observer.prepare_emission(time, segment.evaluate(time, true)?)?;
        let error = residual(receiver, time, prepared.propagation_delay_s());
        if !error.is_finite() || error > tolerance {
            return Err(AeroacError::NotConverged {
                what: "trajectory retarded-time residual at represented source time",
                residual: error,
            });
        }
        if receiver <= time {
            return Err(invalid("trajectory arrival is not representably causal"));
        }
        prepared.sample.receiver_time_s = receiver;
        prepared.sample.retarded_time_residual_s = error;
        poll(cancelled)?;
        Ok(prepared)
    }
}

impl QuinticSegment {
    fn compile(
        first: TimedSourceKinematics,
        last: TimedSourceKinematics,
        first_arrival: f64,
        last_arrival: f64,
    ) -> Result<(Self, BoundNode, bool), AeroacError> {
        let duration = Dd::from_f64(last.source_time_s) - Dd::from_f64(first.source_time_s);
        if !finite_dd(duration) || duration.to_f64() <= 0.0 {
            return Err(invalid("quintic source duration is unrepresentable"));
        }
        let stationary = first.kinematics.position_m == last.kinematics.position_m
            && first.kinematics.velocity_m_s == [0.0; 3]
            && last.kinematics.velocity_m_s == [0.0; 3]
            && first.kinematics.acceleration_m_s2 == [0.0; 3]
            && last.kinematics.acceleration_m_s2 == [0.0; 3];
        let mut position = [[Dd::ZERO; 3]; 6];
        let mut velocity = [[Dd::ZERO; 3]; 5];
        let mut acceleration = [[Dd::ZERO; 3]; 4];
        let zero = Interval::point(0.0);
        let mut bounds = BoundNode {
            position: [[zero; 3]; 6],
            velocity: [[zero; 3]; 5],
            depth: 0,
        };
        if !stationary {
            let h = Interval::point(last.source_time_s) - Interval::point(first.source_time_s);
            if h.lo() <= 0.0 || !h.hi().is_finite() {
                return Err(invalid("quintic duration enclosure"));
            }
            for axis in 0..3 {
                let a = first.kinematics;
                let b = last.kinematics;
                let delta = Dd::from_f64(b.position_m[axis]) - Dd::from_f64(a.position_m[axis]);
                let v0 = scaled_velocity(duration, a.velocity_m_s[axis]);
                let v1 = scaled_velocity(duration, b.velocity_m_s[axis]);
                let a0 = scaled_acceleration(duration, a.acceleration_m_s2[axis]);
                let a1 = scaled_acceleration(duration, b.acceleration_m_s2[axis]);
                if (a.velocity_m_s[axis] != 0.0 && v0 == Dd::ZERO)
                    || (b.velocity_m_s[axis] != 0.0 && v1 == Dd::ZERO)
                    || (a.acceleration_m_s2[axis] != 0.0 && a0 == Dd::ZERO)
                    || (b.acceleration_m_s2[axis] != 0.0 && a1 == Dd::ZERO)
                {
                    return Err(invalid("quintic endpoint derivative underflow"));
                }
                let controls = [
                    Dd::ZERO,
                    v0,
                    Dd::from_f64(2.0) * v0 + a0,
                    delta - Dd::from_f64(2.0) * v1 + a1,
                    delta - v1,
                    delta,
                ];
                let delta_bound =
                    Interval::point(b.position_m[axis]) - Interval::point(a.position_m[axis]);
                let v0b = velocity_bound(h, a.velocity_m_s[axis]);
                let v1b = velocity_bound(h, b.velocity_m_s[axis]);
                let a0b = acceleration_bound(h, a.acceleration_m_s2[axis]);
                let a1b = acceleration_bound(h, b.acceleration_m_s2[axis]);
                let enclosure = [
                    zero,
                    v0b,
                    Interval::point(2.0) * v0b + a0b,
                    delta_bound - Interval::point(2.0) * v1b + a1b,
                    delta_bound - v1b,
                    delta_bound,
                ];
                for i in 0..6 {
                    position[i][axis] = controls[i];
                    bounds.position[i][axis] = enclosure[i].hull(enclose_dd(controls[i])?);
                }
                for i in 0..5 {
                    velocity[i][axis] =
                        ((controls[i + 1] - controls[i]) / duration) * Dd::from_f64(5.0);
                    bounds.velocity[i][axis] =
                        (((bounds.position[i + 1][axis] - bounds.position[i][axis]) / h)
                            * Interval::point(5.0))
                        .hull(enclose_dd(velocity[i][axis])?);
                }
                for i in 0..4 {
                    acceleration[i][axis] = ((velocity[i + 1][axis] - velocity[i][axis])
                        / duration)
                        * Dd::from_f64(4.0);
                }
            }
        }
        if position
            .iter()
            .chain(velocity.iter())
            .chain(acceleration.iter())
            .flatten()
            .any(|v| !finite_dd(*v))
        {
            return Err(invalid(
                "quintic position or derivative coefficients are unrepresentable",
            ));
        }
        Ok((
            Self {
                first,
                last,
                duration,
                position,
                velocity,
                acceleration,
                first_arrival,
                last_arrival,
            },
            bounds,
            stationary,
        ))
    }

    fn evaluate(&self, time: f64, acceleration: bool) -> Result<SourceKinematics, AeroacError> {
        if time == self.first.source_time_s {
            return Ok(self.first.kinematics);
        }
        if time == self.last.source_time_s {
            return Ok(self.last.kinematics);
        }
        if !time.is_finite() || time < self.first.source_time_s || time > self.last.source_time_s {
            return Err(invalid(
                "quintic evaluation would extrapolate source motion",
            ));
        }
        let u = (Dd::from_f64(time) - Dd::from_f64(self.first.source_time_s)) / self.duration;
        let offset = evaluate_controls(self.position, u);
        let velocity = evaluate_controls(self.velocity, u);
        let acceleration = if acceleration {
            evaluate_controls(self.acceleration, u)
        } else {
            [Dd::ZERO; 3]
        };
        let result = SourceKinematics {
            position_m: core::array::from_fn(|i| {
                (Dd::from_f64(self.first.kinematics.position_m[i]) + offset[i]).to_f64()
            }),
            velocity_m_s: velocity.map(Dd::to_f64),
            acceleration_m_s2: acceleration.map(Dd::to_f64),
        };
        if result
            .position_m
            .iter()
            .chain(&result.velocity_m_s)
            .chain(&result.acceleration_m_s2)
            .any(|v| !v.is_finite())
        {
            return Err(invalid("quintic kinematics is unrepresentable"));
        }
        Ok(result)
    }
}

fn admit_segment(
    observer: &MovingCompactObserver,
    segment: &QuinticSegment,
    initial: BoundNode,
    stationary: bool,
    budget: TrajectoryBudget,
    visited: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, AeroacError> {
    let relative_receiver = core::array::from_fn(|i| {
        Interval::point(observer.observer_m[i])
            - Interval::point(segment.first.kinematics.position_m[i])
    });
    let mut stack = [initial; MAX_BOUND_DEPTH + 1];
    let mut pending = 1_usize;
    let mut maximum_mach = 0.0_f64;
    while pending > 0 {
        poll(cancelled)?;
        if *visited >= budget.maximum_bound_nodes {
            return Err(AeroacError::NotConverged {
                what: "trajectory whole-segment bound-node budget",
                residual: 1.0,
            });
        }
        *visited += 1;
        pending -= 1;
        let node = stack[pending];
        let mach = if stationary {
            0.0
        } else {
            velocity_hull_mach(&node.velocity, observer.medium.sound_speed_m_s())
        };
        let clear = position_hull_clear(
            &node.position,
            relative_receiver,
            observer.limits.minimum_distance_m,
        );
        if mach.is_finite() && mach <= observer.limits.maximum_mach && clear {
            maximum_mach = maximum_mach.max(mach);
            continue;
        }
        if node.depth >= budget.maximum_subdivision_depth {
            return Err(AeroacError::NotConverged {
                what: "trajectory speed/clearance hull remains inadmissible at depth limit",
                residual: if mach.is_finite() {
                    (mach - observer.limits.maximum_mach).max(if clear { 0.0 } else { 1.0 })
                } else {
                    f64::INFINITY
                },
            });
        }
        let (left_position, right_position) = split_controls(node.position);
        let (left_velocity, right_velocity) = split_controls(node.velocity);
        stack[pending] = BoundNode {
            position: right_position,
            velocity: right_velocity,
            depth: node.depth + 1,
        };
        stack[pending + 1] = BoundNode {
            position: left_position,
            velocity: left_velocity,
            depth: node.depth + 1,
        };
        pending += 2;
    }
    Ok(maximum_mach)
}

fn velocity_hull_mach(controls: &[[Interval; 3]; 5], sound_speed: f64) -> f64 {
    let c = Interval::point(sound_speed);
    let mut maximum = 0.0_f64;
    for control in controls {
        let mut square = Interval::point(0.0);
        for component in control {
            let scaled = (*component / c).abs();
            square = square + scaled * scaled;
        }
        maximum = maximum.max(square.sqrt().hi());
    }
    maximum
}

fn position_hull_clear(
    controls: &[[Interval; 3]; 6],
    receiver: [Interval; 3],
    clearance: f64,
) -> bool {
    let mut square = Interval::point(0.0);
    for axis in 0..3 {
        let mut hull = controls[0][axis];
        for control in &controls[1..] {
            hull = hull.hull(control[axis]);
        }
        let separation = (receiver[axis] - hull).abs();
        let minimum = separation.lo();
        if minimum >= clearance {
            return true;
        }
        if minimum > 0.0 {
            // All components reaching this branch are below the clearance;
            // normalization avoids squaring physical distances that overflow.
            let scaled = Interval::point(minimum) / Interval::point(clearance);
            square = square + scaled * scaled;
        }
    }
    square.lo() >= 1.0
}

fn evaluate_controls<const N: usize>(mut controls: [[Dd; 3]; N], u: Dd) -> [Dd; 3] {
    let complement = Dd::ONE - u;
    for level in 1..N {
        for i in 0..N - level {
            for axis in 0..3 {
                controls[i][axis] = complement * controls[i][axis] + u * controls[i + 1][axis];
            }
        }
    }
    controls[0]
}

fn split_controls<const N: usize>(
    mut controls: [[Interval; 3]; N],
) -> ([[Interval; 3]; N], [[Interval; 3]; N]) {
    let mut left = controls;
    let mut right = controls;
    let half = Interval::point(0.5);
    for level in 1..N {
        for i in 0..N - level {
            for axis in 0..3 {
                // Scaling first also avoids overflowing a finite convex mean.
                controls[i][axis] = half * controls[i][axis] + half * controls[i + 1][axis];
            }
        }
        left[level] = controls[0];
        right[N - 1 - level] = controls[N - 1 - level];
    }
    (left, right)
}

fn scaled_velocity(h: Dd, value: f64) -> Dd {
    if value == 0.0 {
        return Dd::ZERO;
    }
    if h.to_f64() >= 1.0 {
        (h / Dd::from_f64(5.0)) * Dd::from_f64(value)
    } else {
        (h * Dd::from_f64(value)) / Dd::from_f64(5.0)
    }
}
fn scaled_acceleration(h: Dd, value: f64) -> Dd {
    if value == 0.0 {
        return Dd::ZERO;
    }
    if h.to_f64() >= 1.0 {
        ((h / Dd::from_f64(20.0)) * Dd::from_f64(value)) * h
    } else {
        ((Dd::from_f64(value) * h) * h) / Dd::from_f64(20.0)
    }
}
fn velocity_bound(h: Interval, value: f64) -> Interval {
    if value == 0.0 {
        Interval::point(0.0)
    } else {
        (h / Interval::point(5.0)) * Interval::point(value)
    }
}
fn acceleration_bound(h: Interval, value: f64) -> Interval {
    if value == 0.0 {
        Interval::point(0.0)
    } else if h.hi() >= 1.0 {
        ((h / Interval::point(20.0)) * Interval::point(value)) * h
    } else {
        ((Interval::point(value) * h) * h) / Interval::point(20.0)
    }
}
fn enclose_dd(value: Dd) -> Result<Interval, AeroacError> {
    if !finite_dd(value) {
        return Err(invalid("quintic coefficient overflow"));
    }
    if value.lo == 0.0 {
        Ok(Interval::point(value.hi))
    } else {
        Ok(Interval::point(value.hi) + Interval::point(value.lo))
    }
}
fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
fn invalid(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled {
            what: "finite moving source trajectory",
        })
    } else {
        Ok(())
    }
}
