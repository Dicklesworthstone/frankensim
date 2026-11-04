//! Observer-independent accepted paths for sources and microphones.
//!
//! These paths use the original quintic compiler, derivative evaluator and
//! outward Bernstein subdivision. They prescribe motion; no mechanical state
//! is integrated, inferred from audio, or advanced by observation.

use super::{
    BoundNode, MAX_BOUND_DEPTH, MAX_ROOT_ITERATIONS, QuinticSegment, SourceKinematics,
    TimedSourceKinematics, TrajectoryBudget, invalid, norm, poll, split_controls,
    velocity_hull_mach,
};
use crate::AeroacError;
use crate::curle3d::moving::flow::UniformFlowFrame;
use fs_ivl::Interval;
use fs_material::fluid_acoustics::FluidAcousticState;

#[derive(Debug)]
pub(super) struct KinematicSegment {
    pub(super) polynomial: QuinticSegment,
    // Relative to polynomial.first.kinematics.position_m, as in the compiler.
    pub(super) position_bounds: [[Interval; 3]; 6],
}

/// Finite C2 motion with a whole-path, strictly subsonic speed bound.
///
/// Source and microphone paths have the same interpolation owner. The legacy
/// `TimedSourceKinematics::source_time_s` field holds the physical event time
/// when used for either role; it is not shifted by a propagation delay.
/// No receiver, source strength, acoustic clearance or radiation law is bound
/// here. Those belong to the prepared source/receiver pair.
#[derive(Debug)]
pub struct PreparedKinematicTrajectory {
    pub(super) medium: FluidAcousticState,
    flow_frame: Option<UniformFlowFrame>,
    pub(super) segments: Vec<KinematicSegment>,
    pub(super) declared_maximum_mach: f64,
    maximum_mach: f64,
    maximum_root_iterations: usize,
}

impl PreparedKinematicTrajectory {
    /// Compile the original accepted knots and admit speed over every span.
    ///
    /// `maximum_payload_bytes` counts this owner and its retained segment
    /// payload, including the original interval position controls. Allocator
    /// metadata and the fixed, at most 33-node preparation stack are excluded.
    /// The root cap is retained for use when this path is a source.
    pub fn new(
        medium: FluidAcousticState,
        samples: &[TimedSourceKinematics],
        maximum_mach: f64,
        budget: TrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::compile(medium, samples, None, maximum_mach, budget, cancelled)
    }

    /// Prepare accepted laboratory motion in one homogeneous moving fluid.
    ///
    /// Position and velocity are transformed at each accepted event, before
    /// the SAME original C2 quintic is compiled and bounded. Acceleration and
    /// physical timestamps are unchanged. No intermediate knot array is made;
    /// preparation storage/work and allocation-free queries retain the original
    /// budgets. `maximum_mach` bounds speed relative to the fluid, including
    /// every interior point of every segment.
    ///
    /// Prepare the microphone with the same `frame`, then use the existing
    /// `PreparedMovingReceiverSource` and accepted-history renderer. Their
    /// causal solver, pressure coefficients and observed source-clock rate
    /// automatically use both paths' fluid-relative motion. The general
    /// whole-window Doppler bound remains conservative.
    pub fn in_uniform_flow(
        medium: FluidAcousticState,
        samples: &[TimedSourceKinematics],
        frame: UniformFlowFrame,
        maximum_mach: f64,
        budget: TrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::compile(
            medium,
            samples,
            Some(frame),
            maximum_mach,
            budget,
            cancelled,
        )
    }

    fn compile(
        medium: FluidAcousticState,
        samples: &[TimedSourceKinematics],
        flow_frame: Option<UniformFlowFrame>,
        maximum_mach: f64,
        budget: TrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let count = samples.len().checked_sub(1).unwrap_or(0);
        if !maximum_mach.is_finite()
            || !(0.0..1.0).contains(&maximum_mach)
            || count == 0
            || count > budget.maximum_segments
            || budget.maximum_bound_nodes < count
            || budget.maximum_subdivision_depth > MAX_BOUND_DEPTH
            || !(1..=MAX_ROOT_ITERATIONS).contains(&budget.maximum_root_iterations)
        {
            return Err(invalid(
                "kinematic trajectory Mach, segment and work budgets",
            ));
        }
        let bytes = count
            .checked_mul(size_of::<KinematicSegment>())
            .and_then(|n| n.checked_add(size_of::<Self>()))
            .ok_or_else(|| invalid("kinematic trajectory payload overflow"))?;
        if bytes > budget.maximum_payload_bytes {
            return Err(invalid(
                "kinematic trajectory exceeds retained payload budget",
            ));
        }
        let c = medium.sound_speed_m_s();
        let transform = |mut sample: TimedSourceKinematics| -> Result<_, AeroacError> {
            if let Some(frame) = flow_frame {
                sample.kinematics =
                    frame.to_fluid_frame(sample.source_time_s, sample.kinematics)?;
            }
            Ok(sample)
        };
        let mut first = transform(samples[0])?;
        admit_knot(first, c, maximum_mach)?;
        let mut segments = Vec::new();
        segments
            .try_reserve_exact(count)
            .map_err(|_| invalid("cannot reserve kinematic trajectory segments"))?;
        let mut visited = 0_usize;
        let mut admitted_mach = 0.0_f64;
        for pair in samples.windows(2) {
            poll(cancelled)?;
            let last = transform(pair[1])?;
            admit_knot(last, c, maximum_mach)?;
            if last.source_time_s <= first.source_time_s {
                return Err(invalid("kinematic trajectory times must increase strictly"));
            }
            // Arrival fields belong only to the fixed-receiver adapter. The
            // shared polynomial evaluator never reads these placeholders.
            let (polynomial, bounds, stationary) = QuinticSegment::compile(first, last, 0.0, 0.0)?;
            let bound = admit_speed(
                bounds,
                stationary,
                c,
                maximum_mach,
                budget,
                &mut visited,
                cancelled,
            )?;
            admitted_mach = admitted_mach.max(bound);
            segments.push(KinematicSegment {
                polynomial,
                position_bounds: bounds.position,
            });
            first = last;
        }
        poll(cancelled)?;
        Ok(Self {
            medium,
            flow_frame,
            segments,
            declared_maximum_mach: maximum_mach,
            maximum_mach: admitted_mach,
            maximum_root_iterations: budget.maximum_root_iterations,
        })
    }

    /// Original material snapshot, shared by both paths of a prepared pair.
    pub const fn medium(&self) -> FluidAcousticState {
        self.medium
    }
    /// Original explicit laboratory-to-fluid frame, if motion was authored in flow.
    pub const fn flow_frame(&self) -> Option<UniformFlowFrame> {
        self.flow_frame
    }
    /// Complete physical event-time domain [s], inclusive.
    pub fn time_interval_s(&self) -> [f64; 2] {
        [
            self.segments[0].polynomial.first.source_time_s,
            self.segments[self.segments.len() - 1]
                .polynomial
                .last
                .source_time_s,
        ]
    }
    /// Actual outward upper speed/c bound retained by preparation.
    pub const fn maximum_mach(&self) -> f64 {
        self.maximum_mach
    }
    /// Original cap available for each retarded query when used as a source.
    pub const fn maximum_root_iterations(&self) -> usize {
        self.maximum_root_iterations
    }
    /// Original quintic span count; bound subdivision changes no interpolation.
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }
    /// Read the accepted-knot interpolant and its two derivatives in the fluid
    /// rest frame used by the original acoustic solver.
    /// An out-of-domain query refuses instead of extending a trajectory.
    pub fn kinematics_at(&self, time_s: f64) -> Result<SourceKinematics, AeroacError> {
        self.segment_at(time_s)?.evaluate(time_s, true)
    }

    /// Observe the same represented motion in its original laboratory frame.
    /// Quiescent trajectories return their original interpolant directly.
    pub fn laboratory_kinematics_at(&self, time_s: f64) -> Result<SourceKinematics, AeroacError> {
        let state = self.kinematics_at(time_s)?;
        self.flow_frame
            .map_or(Ok(state), |frame| frame.to_laboratory_frame(time_s, state))
    }

    pub(super) fn segment_at(&self, time_s: f64) -> Result<&QuinticSegment, AeroacError> {
        let [first, last] = self.time_interval_s();
        if !time_s.is_finite() || time_s < first || time_s > last {
            return Err(invalid("query is outside retained kinematic trajectory"));
        }
        let index = self
            .segments
            .partition_point(|s| s.polynomial.last.source_time_s < time_s);
        Ok(&self.segments[index].polynomial)
    }
}

fn admit_knot(knot: TimedSourceKinematics, c: f64, maximum_mach: f64) -> Result<(), AeroacError> {
    let state = knot.kinematics;
    if !knot.source_time_s.is_finite()
        || state
            .position_m
            .iter()
            .chain(&state.velocity_m_s)
            .chain(&state.acceleration_m_s2)
            .any(|v| !v.is_finite())
        || norm(state.velocity_m_s.map(|v| v / c)) > maximum_mach
    {
        return Err(invalid(
            "kinematic trajectory knot is nonfinite or exceeds Mach limit",
        ));
    }
    Ok(())
}

fn admit_speed(
    initial: BoundNode,
    stationary: bool,
    c: f64,
    maximum_mach: f64,
    budget: TrajectoryBudget,
    visited: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, AeroacError> {
    let mut stack = [initial; MAX_BOUND_DEPTH + 1];
    let mut pending = 1;
    let mut admitted = 0.0_f64;
    while pending > 0 {
        poll(cancelled)?;
        if *visited >= budget.maximum_bound_nodes {
            return Err(AeroacError::NotConverged {
                what: "kinematic trajectory speed bound-node budget",
                residual: 1.0,
            });
        }
        *visited += 1;
        pending -= 1;
        let node = stack[pending];
        let mach = if stationary {
            0.0
        } else {
            velocity_hull_mach(&node.velocity, c)
        };
        if mach.is_finite() && mach <= maximum_mach {
            admitted = admitted.max(mach);
            continue;
        }
        if node.depth >= budget.maximum_subdivision_depth {
            return Err(AeroacError::NotConverged {
                what: "kinematic trajectory speed hull remains inadmissible",
                residual: (mach - maximum_mach).max(0.0),
            });
        }
        let (left_velocity, right_velocity) = split_controls(node.velocity);
        // Position bounds are retained separately in their original form;
        // this speed-only stack needs only the velocity subdivision.
        stack[pending] = BoundNode {
            velocity: right_velocity,
            depth: node.depth + 1,
            ..node
        };
        stack[pending + 1] = BoundNode {
            velocity: left_velocity,
            depth: node.depth + 1,
            ..node
        };
        pending += 2;
    }
    Ok(admitted)
}
