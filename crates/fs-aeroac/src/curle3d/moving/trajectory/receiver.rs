//! Pressure sampled by a moving microphone on its original physical clock.
//!
//! At each receiver time `t`, solve `t = tau + |x(t)-y(tau)|/c`. The pressure
//! law is still the Eulerian source field evaluated at `x(t)`: source velocity
//! remains relative to the medium, and receiver motion multiplies no pressure
//! coefficient. It changes the observed source clock by
//! `d tau/d t = (1-n.v_receiver/c)/(1-n.v_source/c)`.

use super::path::PreparedKinematicTrajectory;
use super::{
    MAX_BOUND_DEPTH, MAX_ROOT_ITERATIONS, MovingCompactObserver, PreparedMovingCompactSample,
    QuinticSegment, SourceKinematics, dot, finite_dd, invalid, norm, poll, residual,
    split_controls,
};
use crate::AeroacError;
use crate::curle3d::moving::MovingCompactLimits;
use fs_ivl::Interval;
use fs_math::dd::Dd;

mod doppler;

/// The finite observation window and its physical error limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingReceiverConfig {
    /// Original receiver times [s], inclusive; a single timestamp is allowed.
    pub receiver_time_interval_s: [f64; 2],
    /// Minimum separation of admitted retarded source/receiver events [m].
    pub minimum_distance_m: f64,
    /// Absolute `|t-tau-R/c|` tolerance at the represented source timestamp [s].
    pub maximum_time_residual_s: f64,
}

/// Finite work limits for pair preparation and every retarded query.
#[derive(Clone, Copy, Debug)]
pub struct MovingReceiverBudget {
    /// Original source/receiver span pairs visited during clearance admission.
    pub maximum_segment_pairs: usize,
    /// Total causal Bernstein pair-bound nodes across those span pairs.
    pub maximum_bound_nodes: usize,
    /// Total alternating source/receiver subdivisions per node, in 0..=32.
    pub maximum_subdivision_depth: usize,
    /// Interior root refinements per query, in 1..=256; also limited by the
    /// source path's original root cap. Knot search adds O(log(source spans)).
    pub maximum_root_iterations: usize,
}

/// One moving-microphone sample, with pressure and clock metadata kept distinct.
#[derive(Clone, Copy, Debug)]
pub struct PreparedMovingReceiverSample {
    /// Original Eulerian pressure coefficients at the microphone's position.
    /// Its `source_time_rate()` remains the fixed-position factor `1/kappa`.
    pub pressure: PreparedMovingCompactSample,
    /// Accepted-knot microphone kinematics in the shared fluid rest frame at
    /// exactly the requested time. Its path also exposes laboratory kinematics.
    pub receiver_kinematics: SourceKinematics,
    /// Actual `d tau/d t` along this microphone path, including both motions.
    pub observed_source_time_rate: f64,
}

/// Borrowed source/receiver motion admitted for one finite observation window.
///
/// Both paths retain one original homogeneous, lossless medium and the same
/// explicit rest frame. Laboratory paths in a uniform flow share one original
/// `UniformFlowFrame`; mixing independently shifted frames refuses.
/// Endpoint roots prove continuous history coverage because both paths are
/// strictly subsonic. Clearance is admitted with conservative interval bounds
/// over source/receiver span pairs, pruning events outside the sound cone.
/// Inconclusive bounds refuse under explicit work caps; no knot is changed.
/// Queries allocate nothing and never advance either mechanical trajectory.
#[derive(Clone, Copy, Debug)]
pub struct PreparedMovingReceiverSource<'a> {
    source: &'a PreparedKinematicTrajectory,
    receiver: &'a PreparedKinematicTrajectory,
    config: MovingReceiverConfig,
    maximum_source_time_rate: f64,
    maximum_root_iterations: usize,
}

impl<'a> PreparedMovingReceiverSource<'a> {
    /// Admit shared material, finite clocks, retarded coverage and clearance.
    /// Two endpoint root solves and the pair-bound work are preparation-only.
    pub fn new(
        source: &'a PreparedKinematicTrajectory,
        receiver: &'a PreparedKinematicTrajectory,
        config: MovingReceiverConfig,
        budget: MovingReceiverBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let [first, last] = config.receiver_time_interval_s;
        let [receiver_first, receiver_last] = receiver.time_interval_s();
        if source.medium() != receiver.medium()
            || source.flow_frame() != receiver.flow_frame()
            || !first.is_finite()
            || !last.is_finite()
            || first > last
            || first < receiver_first
            || last > receiver_last
            || !config.minimum_distance_m.is_finite()
            || config.minimum_distance_m <= 0.0
            || !config.maximum_time_residual_s.is_finite()
            || config.maximum_time_residual_s < 0.0
            || budget.maximum_segment_pairs == 0
            || budget.maximum_bound_nodes == 0
            || budget.maximum_subdivision_depth > MAX_BOUND_DEPTH
            || !(1..=MAX_ROOT_ITERATIONS).contains(&budget.maximum_root_iterations)
        {
            return Err(invalid(
                "moving receiver medium, finite window, clearance or work budget",
            ));
        }
        let maximum_source_time_rate = doppler::maximum_rate(source, receiver);
        if !maximum_source_time_rate.is_finite() || maximum_source_time_rate <= 0.0 {
            return Err(invalid("moving receiver Doppler bound is unrepresentable"));
        }
        let owner = Self {
            source,
            receiver,
            config,
            maximum_source_time_rate,
            maximum_root_iterations: budget
                .maximum_root_iterations
                .min(source.maximum_root_iterations()),
        };
        // Since d(tau)/dt is positive, coverage at both receiver endpoints
        // covers every intermediate receiver time. These queries do not need
        // an already-proved clearance hull: only their final roots gate it.
        owner.prepare_receiver(first, cancelled)?;
        if last != first {
            owner.prepare_receiver(last, cancelled)?;
        }
        admit_clearance(&owner, budget, cancelled)?;
        poll(cancelled)?;
        Ok(owner)
    }

    /// Original finite source path; never resampled or mechanically advanced.
    pub const fn source(&self) -> &'a PreparedKinematicTrajectory {
        self.source
    }
    /// Original finite microphone path, on its own accepted physical times.
    pub const fn receiver(&self) -> &'a PreparedKinematicTrajectory {
        self.receiver
    }
    /// The admitted finite receiver window and residual/clearance limits.
    pub const fn config(&self) -> MovingReceiverConfig {
        self.config
    }
    /// Continuous upper observed source-clock rate. Intersects the original
    /// independent-speed bound with a signed relative-velocity hull bound;
    /// common translation therefore does not imply opposing-motion Doppler.
    /// The bound includes velocities at DIFFERENT emission/receiver times.
    /// Additional geometric amplitude modulation still needs a bandwidth limit.
    pub const fn maximum_source_time_rate(&self) -> f64 {
        self.maximum_source_time_rate
    }
    /// Effective per-query root refinement cap, preserving the source's cap.
    pub const fn maximum_root_iterations(&self) -> usize {
        self.maximum_root_iterations
    }

    /// Observe the source at one original microphone timestamp.
    ///
    /// Every trial fixes the microphone at `x(t)` while evaluating the source
    /// at its trial emission time. Trial positions may cross the microphone;
    /// only the final causal event must satisfy the distance limit.
    pub fn prepare_receiver(
        &self,
        receiver_time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingReceiverSample, AeroacError> {
        poll(cancelled)?;
        let [first, last] = self.config.receiver_time_interval_s;
        if !receiver_time.is_finite() || receiver_time < first || receiver_time > last {
            return Err(invalid(
                "query is outside the admitted moving receiver window",
            ));
        }
        let receiver_state = self.receiver.kinematics_at(receiver_time)?;
        let first_segment = &self.source.segments[0].polynomial;
        let last_segment = &self.source.segments[self.source.segments.len() - 1].polynomial;
        let first_time = first_segment.first.source_time_s;
        let last_time = last_segment.last.source_time_s;
        let first_value = self
            .arrival_residual(first_segment, first_time, receiver_time, receiver_state)?
            .0;
        let last_value = self
            .arrival_residual(last_segment, last_time, receiver_time, receiver_state)?
            .0;
        // Reconcile only an exactly matching rounded knot arrival. A residual
        // tolerance alone cannot authorize extrapolation beyond the history.
        if rounded_arrival_matches(first_time, receiver_time, first_value) {
            return self.finish(
                first_segment,
                first_time,
                receiver_time,
                receiver_state,
                cancelled,
            );
        }
        if rounded_arrival_matches(last_time, receiver_time, last_value) {
            return self.finish(
                last_segment,
                last_time,
                receiver_time,
                receiver_state,
                cancelled,
            );
        }
        if Dd::ZERO.lt(first_value) || last_value.lt(Dd::ZERO) {
            return Err(invalid(
                "moving receiver has no retarded event in the retained source history",
            ));
        }
        // Arrival is monotone in emission time at this fixed x(t). Search
        // original source knots without a per-microphone arrival cache.
        let mut begin = 0;
        let mut end = self.source.segments.len();
        while begin < end {
            poll(cancelled)?;
            let middle = begin + (end - begin) / 2;
            let segment = &self.source.segments[middle].polynomial;
            let value = self
                .arrival_residual(
                    segment,
                    segment.last.source_time_s,
                    receiver_time,
                    receiver_state,
                )?
                .0;
            if value.lt(Dd::ZERO) {
                begin = middle + 1;
            } else {
                end = middle;
            }
        }
        let index = begin.min(self.source.segments.len() - 1);
        let segment = &self.source.segments[index].polynomial;
        let mut lo = segment.first.source_time_s;
        let mut hi = segment.last.source_time_s;
        let mut flo = self
            .arrival_residual(segment, lo, receiver_time, receiver_state)?
            .0;
        let mut fhi = self
            .arrival_residual(segment, hi, receiver_time, receiver_state)?
            .0;
        for (time, value) in [(lo, flo), (hi, fhi)] {
            if rounded_arrival_matches(time, receiver_time, value)
                || value.abs().to_f64() <= self.config.maximum_time_residual_s
            {
                return self.finish(segment, time, receiver_time, receiver_state, cancelled);
            }
        }
        if !flo.lt(Dd::ZERO) || !Dd::ZERO.lt(fhi) {
            return Err(invalid("moving receiver retarded segment is not bracketed"));
        }
        let fraction = (Dd::ZERO - flo) / (fhi - flo);
        let mut time = (Dd::from_f64(lo) + segment.duration * fraction).to_f64();
        if time <= lo || time >= hi {
            time = lo.midpoint(hi);
        }
        let mut best_error = flo.abs().to_f64().min(fhi.abs().to_f64());
        for iteration in 0..self.maximum_root_iterations {
            poll(cancelled)?;
            if time <= lo || time >= hi {
                let endpoint = if flo.abs().lt(fhi.abs()) { lo } else { hi };
                return self.finish(segment, endpoint, receiver_time, receiver_state, cancelled);
            }
            let (value, derivative) =
                self.arrival_residual(segment, time, receiver_time, receiver_state)?;
            let error = value.abs().to_f64();
            best_error = best_error.min(error);
            if error <= self.config.maximum_time_residual_s {
                return self.finish(segment, time, receiver_time, receiver_state, cancelled);
            }
            if value.lt(Dd::ZERO) {
                lo = time;
                flo = value;
            } else {
                hi = time;
                fhi = value;
            }
            let newton = if derivative > 0.0 {
                (Dd::from_f64(time) - value / Dd::from_f64(derivative)).to_f64()
            } else {
                f64::NAN
            };
            time = if iteration % 4 != 3 && newton.is_finite() && lo < newton && newton < hi {
                newton
            } else {
                lo.midpoint(hi)
            };
        }
        Err(AeroacError::NotConverged {
            what: "moving receiver retarded-time iteration budget",
            residual: best_error,
        })
    }

    fn arrival_residual(
        &self,
        segment: &QuinticSegment,
        time: f64,
        receiver_time: f64,
        receiver: SourceKinematics,
    ) -> Result<(Dd, f64), AeroacError> {
        let source = segment.evaluate(time, false)?;
        let displacement = core::array::from_fn(|i| receiver.position_m[i] - source.position_m[i]);
        let distance = norm(displacement);
        let c = self.source.medium().sound_speed_m_s();
        let delay = distance / c;
        if !distance.is_finite() || !delay.is_finite() {
            return Err(invalid("moving receiver trial distance is unrepresentable"));
        }
        // At R=0 the norm derivative is undefined. The residual is still
        // finite and monotone; use a bracket step instead of a Newton step.
        let derivative = if distance == 0.0 {
            0.0
        } else {
            1.0 - dot(
                displacement.map(|v| v / distance),
                source.velocity_m_s.map(|v| v / c),
            )
        };
        let value = Dd::from_f64(time) - Dd::from_f64(receiver_time) + Dd::from_f64(delay);
        if !finite_dd(value) || !derivative.is_finite() || derivative < 0.0 {
            return Err(invalid("moving receiver trial residual is unrepresentable"));
        }
        Ok((value, derivative))
    }

    fn finish(
        &self,
        segment: &QuinticSegment,
        time: f64,
        receiver_time: f64,
        receiver_state: SourceKinematics,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingReceiverSample, AeroacError> {
        let source_state = segment.evaluate(time, true)?;
        let observer = MovingCompactObserver::new(
            self.source.medium(),
            receiver_state.position_m,
            MovingCompactLimits {
                maximum_mach: self.source.declared_maximum_mach,
                minimum_distance_m: self.config.minimum_distance_m,
            },
        )?;
        let mut pressure = observer.prepare_emission(time, source_state)?;
        let error = residual(receiver_time, time, pressure.propagation_delay_s());
        if !error.is_finite() || error > self.config.maximum_time_residual_s {
            return Err(AeroacError::NotConverged {
                what: "moving receiver retarded residual at represented source time",
                residual: error,
            });
        }
        if receiver_time <= time {
            return Err(invalid(
                "moving receiver arrival is not representably causal",
            ));
        }
        let direction = core::array::from_fn(|i| {
            (receiver_state.position_m[i] - source_state.position_m[i]) / pressure.distance_m()
        });
        let c = self.source.medium().sound_speed_m_s();
        let numerator = 1.0 - dot(direction, receiver_state.velocity_m_s.map(|v| v / c));
        let observed_source_time_rate = numerator * pressure.source_time_rate();
        if !observed_source_time_rate.is_finite() || observed_source_time_rate <= 0.0 {
            return Err(invalid(
                "moving receiver observed clock rate is unrepresentable",
            ));
        }
        pressure.sample.receiver_time_s = receiver_time;
        pressure.sample.retarded_time_residual_s = error;
        poll(cancelled)?;
        Ok(PreparedMovingReceiverSample {
            pressure,
            receiver_kinematics: receiver_state,
            observed_source_time_rate,
        })
    }
}

fn rounded_arrival_matches(source_time: f64, receiver_time: f64, signed_residual: Dd) -> bool {
    // signed_residual = source_time - receiver_time + delay; DD cancellation
    // recovers the original rounded arrival without subtracting a large epoch.
    let arrival = (Dd::from_f64(receiver_time) + signed_residual).to_f64();
    source_time < receiver_time && arrival == receiver_time
}

#[derive(Clone, Copy)]
struct PairNode {
    source_position: [[Interval; 3]; 6],
    receiver_position: [[Interval; 3]; 6],
    source_first: Interval,
    source_last: Interval,
    receiver_first: Interval,
    receiver_last: Interval,
    depth: usize,
}

fn admit_clearance(
    owner: &PreparedMovingReceiverSource<'_>,
    budget: MovingReceiverBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AeroacError> {
    let [first, last] = owner.config.receiver_time_interval_s;
    let window = Interval::new(first, last);
    let c = Interval::point(owner.source.medium().sound_speed_m_s());
    let tolerance = owner.config.maximum_time_residual_s;
    let residual_interval = Interval::new(-tolerance, tolerance);
    let mut visited = 0_usize;
    let mut pairs = 0_usize;
    let first_segment = owner
        .receiver
        .segments
        .partition_point(|segment| segment.polynomial.last.source_time_s < first);
    for receiver in &owner.receiver.segments[first_segment..] {
        poll(cancelled)?;
        let r = &receiver.polynomial;
        if r.first.source_time_s > last {
            break;
        }
        for source in &owner.source.segments {
            poll(cancelled)?;
            if pairs >= budget.maximum_segment_pairs {
                return Err(invalid(
                    "moving receiver clearance exceeds segment-pair budget",
                ));
            }
            pairs += 1;
            let s = &source.polynomial;
            let origin_delta = core::array::from_fn(|i| {
                Interval::point(r.first.kinematics.position_m[i])
                    - Interval::point(s.first.kinematics.position_m[i])
            });
            let initial = PairNode {
                source_position: source.position_bounds,
                receiver_position: receiver.position_bounds,
                source_first: Interval::point(s.first.source_time_s),
                source_last: Interval::point(s.last.source_time_s),
                receiver_first: Interval::point(r.first.source_time_s),
                receiver_last: Interval::point(r.last.source_time_s),
                depth: 0,
            };
            let mut stack = [initial; MAX_BOUND_DEPTH + 1];
            let mut pending = 1;
            while pending > 0 {
                poll(cancelled)?;
                if visited >= budget.maximum_bound_nodes {
                    return Err(AeroacError::NotConverged {
                        what: "moving receiver clearance bound-node budget",
                        residual: 1.0,
                    });
                }
                visited += 1;
                pending -= 1;
                let node = stack[pending];
                let Some(receiver_time) =
                    Interval::new(node.receiver_first.lo(), node.receiver_last.hi())
                        .intersect(window)
                else {
                    continue;
                };
                let source_time = Interval::new(node.source_first.lo(), node.source_last.hi());
                let travel_time = receiver_time - source_time + residual_interval;
                if travel_time.hi() < 0.0 {
                    continue;
                }
                let causal_distance =
                    Interval::new(travel_time.lo().max(0.0), travel_time.hi().max(0.0)) * c;
                if causal_distance.lo() >= owner.config.minimum_distance_m {
                    continue;
                }
                let distance = distance_hull(&node, origin_delta);
                if distance.lo() >= owner.config.minimum_distance_m
                    || causal_distance.hi() < distance.lo()
                    || distance.hi() < causal_distance.lo()
                {
                    continue;
                }
                if node.depth >= budget.maximum_subdivision_depth {
                    return Err(AeroacError::NotConverged {
                        what: "moving receiver retarded-clearance hull remains inconclusive",
                        residual: (owner.config.minimum_distance_m - distance.lo()).max(0.0),
                    });
                }
                let half = Interval::point(0.5);
                let (left, right) = if node.depth % 2 == 0 {
                    let (a, b) = split_controls(node.source_position);
                    let middle = half * node.source_first + half * node.source_last;
                    (
                        PairNode {
                            source_position: a,
                            source_last: middle,
                            depth: node.depth + 1,
                            ..node
                        },
                        PairNode {
                            source_position: b,
                            source_first: middle,
                            depth: node.depth + 1,
                            ..node
                        },
                    )
                } else {
                    let (a, b) = split_controls(node.receiver_position);
                    let middle = half * node.receiver_first + half * node.receiver_last;
                    (
                        PairNode {
                            receiver_position: a,
                            receiver_last: middle,
                            depth: node.depth + 1,
                            ..node
                        },
                        PairNode {
                            receiver_position: b,
                            receiver_first: middle,
                            depth: node.depth + 1,
                            ..node
                        },
                    )
                };
                stack[pending] = right;
                stack[pending + 1] = left;
                pending += 2;
            }
        }
    }
    Ok(())
}

fn distance_hull(node: &PairNode, origin_delta: [Interval; 3]) -> Interval {
    let components: [Interval; 3] = core::array::from_fn(|axis| {
        let mut source = node.source_position[0][axis];
        let mut receiver = node.receiver_position[0][axis];
        for i in 1..6 {
            source = source.hull(node.source_position[i][axis]);
            receiver = receiver.hull(node.receiver_position[i][axis]);
        }
        (origin_delta[axis] + receiver - source).abs()
    });
    let scale = components.iter().fold(0.0_f64, |a, b| a.max(b.hi()));
    if !scale.is_finite() {
        return Interval::new(0.0, f64::INFINITY);
    }
    if scale == 0.0 {
        return Interval::point(0.0);
    }
    let scale = Interval::point(scale);
    let mut square = Interval::point(0.0);
    for component in components {
        let normalized = component / scale;
        square = square + normalized * normalized;
    }
    square.sqrt() * scale
}
