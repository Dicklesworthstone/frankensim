//! Ideal-wall images at a microphone's original physical arrival event.
//!
//! A stationary reflection H changes the microphone to H^-1 x(t), including
//! its velocity in the observed-clock derivative. The source interpolant,
//! force frame, source clock and Eulerian pressure law remain unchanged.
//! Every path uses the original moving-pair causal solver and interval bounds.

use super::super::receiver::{ReceiverImage, ReceiverPreparationWork};
use super::super::{
    MAX_BOUND_DEPTH, MovingReceiverBudget, MovingReceiverConfig, PreparedKinematicTrajectory,
    PreparedMovingReceiverSample, PreparedMovingReceiverSource, invalid, poll, split_controls,
};
use super::{PositionNode, ReflectedTrajectoryBudget, admit_planes, apply_parity};
use crate::AeroacError;
use crate::curle3d::IdealAcousticPlane;
use crate::curle3d::moving::{
    CompactSourceStrength, MovingCompactObserver, PreparedMovingCompactSample,
};
use fs_ivl::Interval;
use fs_math::dd::Dd;

/// One original C2 microphone path and one complete orthogonal ideal-wall
/// group. Every source in a coherent render borrows this same owner.
///
/// The walls are stationary in the homogeneous lossless fluid-rest frame.
/// Uniform-flow paths refuse: their laboratory walls require a separately
/// admitted frame transform. The receiver may reach a wall when its complete
/// interpolant can be proved inside; an inconclusive rounded hull refuses.
#[derive(Debug)]
pub struct MovingReflectionObserver<'a> {
    receiver: &'a PreparedKinematicTrajectory,
    planes: [IdealAcousticPlane; 3],
    plane_count: usize,
    signs: [f64; 8],
}

impl<'a> MovingReflectionObserver<'a> {
    /// Admit all 2^m paths and the original receiver's complete half-space
    /// confinement once. The bound budget counts ORIGINAL receiver spans and
    /// their subdivisions; no virtual trajectory is allocated or interpolated.
    pub fn new(
        receiver: &'a PreparedKinematicTrajectory,
        planes: &[IdealAcousticPlane],
        maximum_source_paths: usize,
        budget: ReflectedTrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let path_count = admit_planes(planes, maximum_source_paths)?;
        if receiver.flow_frame().is_some() {
            return Err(invalid(
                "moving wall images require paths in the stationary fluid frame",
            ));
        }
        admit_halfspaces(receiver, planes, false, budget, cancelled)?;
        let mut retained = [planes[0]; 3];
        retained[..planes.len()].copy_from_slice(planes);
        let mut signs = [1.0; 8];
        for (mask, sign) in signs.iter_mut().enumerate().take(path_count) {
            for (index, plane) in planes.iter().enumerate() {
                if mask & (1 << index) != 0 {
                    *sign *= plane.boundary().sign();
                }
            }
        }
        poll(cancelled)?;
        Ok(Self {
            receiver,
            planes: retained,
            plane_count: planes.len(),
            signs,
        })
    }

    /// Original microphone trajectory, without image copies or clock shifts.
    pub const fn receiver(&self) -> &'a PreparedKinematicTrajectory {
        self.receiver
    }
    pub fn planes(&self) -> &[IdealAcousticPlane] {
        &self.planes[..self.plane_count]
    }
    /// Complete path count per physical source, including direct path zero.
    pub const fn path_count(&self) -> usize {
        1 << self.plane_count
    }
}

/// All ideal-wall paths from one original source to one moving microphone.
///
/// The small fixed array retains eight borrowed pair views at most. Source
/// and receiver histories remain singular; preparation and pressure queries
/// allocate nothing. Each path has its own retarded time and observed Doppler
/// rate. Missing source coverage refuses the COMPLETE observation window.
#[derive(Debug)]
pub struct PreparedReflectedMovingReceiverSource<'a> {
    observer: &'a MovingReflectionObserver<'a>,
    paths: [PreparedMovingReceiverSource<'a>; 8],
    maximum_source_time_rate: f64,
}

impl<'a> PreparedReflectedMovingReceiverSource<'a> {
    /// Admit the original source's whole half-space confinement and every
    /// image's causal coverage and clearance. `pair_budget` is AGGREGATE across
    /// all paths: source/receiver span pairs and bound nodes are never reset
    /// for a new image. Its root cap remains a per-path query cap.
    pub fn new(
        source: &'a PreparedKinematicTrajectory,
        observer: &'a MovingReflectionObserver<'a>,
        config: MovingReceiverConfig,
        source_bounds: ReflectedTrajectoryBudget,
        pair_budget: MovingReceiverBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        if source.flow_frame().is_some() {
            return Err(invalid(
                "reflected moving sources require the stationary fluid frame",
            ));
        }
        admit_halfspaces(source, observer.planes(), true, source_bounds, cancelled)?;
        let mut work = ReceiverPreparationWork::default();
        let direct = PreparedMovingReceiverSource::new_image(
            source,
            observer.receiver,
            None,
            config,
            pair_budget,
            &mut work,
            cancelled,
        )?;
        let mut paths = [direct; 8];
        let mut maximum_source_time_rate = direct.maximum_source_time_rate();
        for (mask, path) in paths
            .iter_mut()
            .enumerate()
            .take(observer.path_count())
            .skip(1)
        {
            poll(cancelled)?;
            *path = PreparedMovingReceiverSource::new_image(
                source,
                observer.receiver,
                Some(ReceiverImage {
                    planes: observer.planes(),
                    mask,
                }),
                config,
                pair_budget,
                &mut work,
                cancelled,
            )?;
            maximum_source_time_rate =
                maximum_source_time_rate.max(path.maximum_source_time_rate());
        }
        poll(cancelled)?;
        Ok(Self {
            observer,
            paths,
            maximum_source_time_rate,
        })
    }

    pub const fn reflection_observer(&self) -> &'a MovingReflectionObserver<'a> {
        self.observer
    }
    pub const fn source(&self) -> &'a PreparedKinematicTrajectory {
        self.paths[0].source()
    }
    pub const fn receiver(&self) -> &'a PreparedKinematicTrajectory {
        self.observer.receiver
    }
    pub const fn config(&self) -> MovingReceiverConfig {
        self.paths[0].config()
    }
    pub const fn path_count(&self) -> usize {
        self.observer.path_count()
    }
    /// Upper source-clock rate covering EVERY retained path and both original
    /// motions; receiver motion is a sampling rate, never a pressure gain.
    pub const fn maximum_source_time_rate(&self) -> f64 {
        self.maximum_source_time_rate
    }
    pub const fn maximum_root_iterations(&self) -> usize {
        self.paths[0].maximum_root_iterations()
    }
    pub fn path_maximum_source_time_rate(&self, path: usize) -> Result<f64, AeroacError> {
        Ok(self.path(path)?.maximum_source_time_rate())
    }

    fn path(&self, index: usize) -> Result<&PreparedMovingReceiverSource<'a>, AeroacError> {
        if index >= self.path_count() {
            return Err(invalid(
                "moving reflection path index exceeds the complete image group",
            ));
        }
        Ok(&self.paths[index])
    }

    /// Physical microphone position with this source's original material and
    /// clearance limits. The direct path is preserved bit-for-bit.
    pub fn observer_at(&self, time: f64) -> Result<MovingCompactObserver, AeroacError> {
        self.paths[0].observer_at(time)
    }

    /// Prepare one image on the ORIGINAL source history. The pressure already
    /// includes rigid/release parity. Receiver kinematics remain physical;
    /// the observed source-clock rate includes the virtual receiver's velocity.
    pub fn prepare_path(
        &self,
        path: usize,
        receiver_time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingReceiverSample, AeroacError> {
        poll(cancelled)?;
        let mut sample = self
            .path(path)?
            .prepare_receiver(receiver_time, cancelled)?;
        sample.pressure = apply_parity(sample.pressure, self.observer.signs[path]);
        Ok(sample)
    }

    /// Original emission endpoint at the specified physical microphone event.
    /// This supports exact rounded-arrival reconciliation; it never supplies
    /// extrapolation, a root clamp or missing source-strength history.
    pub fn prepare_path_emission(
        &self,
        path: usize,
        source_time: f64,
        receiver_time: f64,
    ) -> Result<PreparedMovingCompactSample, AeroacError> {
        let observer = self.path(path)?.observer_at(receiver_time)?;
        Ok(apply_parity(
            observer.prepare_emission(source_time, self.source().kinematics_at(source_time)?)?,
            self.observer.signs[path],
        ))
    }

    /// Coherent pressure from one physical Q/F history, queried on its own
    /// clock at each path's retarded time. A failed history/path aborts the sum.
    pub fn pressure_at(
        &self,
        receiver_time: f64,
        source_at: &mut impl FnMut(f64) -> Result<CompactSourceStrength, AeroacError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, AeroacError> {
        let mut pressure = Dd::ZERO;
        for path in 0..self.path_count() {
            let sample = self.prepare_path(path, receiver_time, cancelled)?.pressure;
            pressure = pressure
                + Dd::from_f64(
                    sample
                        .evaluate(source_at(sample.source_time_s())?)?
                        .pressure_pa,
                );
        }
        poll(cancelled)?;
        let value = pressure.to_f64();
        if !value.is_finite() {
            return Err(invalid("moving reflected pressure is unrepresentable"));
        }
        Ok(value)
    }
}

fn admit_halfspaces(
    path: &PreparedKinematicTrajectory,
    planes: &[IdealAcousticPlane],
    strict: bool,
    budget: ReflectedTrajectoryBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AeroacError> {
    if budget.maximum_bound_nodes < path.segment_count()
        || budget.maximum_subdivision_depth > MAX_BOUND_DEPTH
    {
        return Err(invalid("moving reflection half-space bound budget"));
    }
    let mut visited = 0;
    for segment in &path.segments {
        let polynomial = &segment.polynomial;
        for knot in [polynomial.first, polynomial.last] {
            for plane in planes {
                let distance = plane.distance(knot.kinematics.position_m)?.to_f64();
                if distance < 0.0 || (strict && distance == 0.0) {
                    return Err(invalid(
                        "moving reflection path leaves its physical half-space",
                    ));
                }
            }
        }
        let origin = polynomial.first.kinematics.position_m;
        let initial = PositionNode {
            controls: segment.position_bounds,
            depth: 0,
        };
        let mut stack = [initial; MAX_BOUND_DEPTH + 1];
        let mut pending = 1;
        while pending > 0 {
            poll(cancelled)?;
            if visited >= budget.maximum_bound_nodes {
                return Err(AeroacError::NotConverged {
                    what: "moving reflection half-space bound-node budget",
                    residual: 1.0,
                });
            }
            visited += 1;
            pending -= 1;
            let node = stack[pending];
            let inside = planes.iter().all(|plane| {
                // The compiler's mathematical polynomial is exactly constant
                // in this coordinate even if a generic interval subtraction
                // widened 0-0. This admits tangential motion ON an axis wall.
                let exactly_on_axis_wall = !strict
                    && (0..3).any(|axis| {
                        plane.inward_normal()[axis].abs() == 1.0
                            && (0..3)
                                .all(|other| other == axis || plane.inward_normal()[other] == 0.0)
                            && origin[axis] == plane.point_m()[axis]
                            && polynomial.last.kinematics.position_m[axis] == origin[axis]
                            && polynomial
                                .position
                                .iter()
                                .all(|control| control[axis] == Dd::ZERO)
                    });
                exactly_on_axis_wall
                    || node.controls.iter().all(|control| {
                        let mut distance = Interval::point(0.0);
                        for axis in 0..3 {
                            let normal = plane.inward_normal()[axis];
                            if normal == 0.0 {
                                continue;
                            }
                            // Exact zero components preserve an exactly represented
                            // axis-aligned boundary trajectory. No near-zero distance
                            // or nonzero normal component is snapped to the wall.
                            let offset = control[axis];
                            if offset.lo() == 0.0
                                && offset.hi() == 0.0
                                && origin[axis] == plane.point_m()[axis]
                            {
                                continue;
                            }
                            let position = if offset.lo() == 0.0 && offset.hi() == 0.0 {
                                Interval::point(origin[axis])
                            } else {
                                Interval::point(origin[axis]) + offset
                            };
                            let term = (position - Interval::point(plane.point_m()[axis]))
                                * Interval::point(normal);
                            distance = if distance.lo() == 0.0 && distance.hi() == 0.0 {
                                term
                            } else {
                                distance + term
                            };
                        }
                        if strict {
                            distance.lo() > 0.0
                        } else {
                            distance.lo() >= 0.0
                        }
                    })
            });
            if inside {
                continue;
            }
            if node.depth >= budget.maximum_subdivision_depth {
                return Err(AeroacError::NotConverged {
                    what: "moving reflection whole-path half-space hull remains inconclusive",
                    residual: 1.0,
                });
            }
            let (left, right) = split_controls(node.controls);
            stack[pending] = PositionNode {
                controls: right,
                depth: node.depth + 1,
            };
            stack[pending + 1] = PositionNode {
                controls: left,
                depth: node.depth + 1,
            };
            pending += 2;
        }
    }
    Ok(())
}
