//! Retarded ideal-wall reflection on an original moving source history.
//!
//! For a stationary orthogonal image transform H, the Euclidean Green function
//! obeys G(x,H y)=G(H^-1 x,y). Using a virtual receiver H^-1 x therefore keeps
//! the ORIGINAL source position, velocity, acceleration, force and force rate.
//! It is exactly equivalent to transforming every image-source vector, without
//! cloning or reinterpolating the trajectory. Rigid/release parity multiplies
//! all eight original pressure coefficients. Each path keeps its own retarded
//! time; a reflected transient is not a gain applied to the direct waveform.
//!
//! The medium is homogeneous, lossless and at rest in the stationary walls'
//! frame. This is the finite ideal-plane image construction of Allen & Berkley,
//! JASA 65 (1979), 943–950, doi:10.1121/1.382599, with the moving Green-function
//! derivatives owned by the original compact source operator. Finite impedance,
//! general wedges, moving walls, dispersion and feedback are absent. Moving
//! microphones use the same image law through [`moving`].

use super::{
    MAX_BOUND_DEPTH, PreparedSourceTrajectory, QuinticSegment, dot, enclose_dd, invalid, poll,
    position_hull_clear, split_controls,
};
use crate::AeroacError;
use crate::curle3d::IdealAcousticPlane;
use crate::curle3d::moving::{
    CompactSourceStrength, MovingCompactObserver, PreparedMovingCompactSample,
};
use fs_ivl::Interval;
use fs_math::dd::Dd;

/// Reflected observation on one original moving microphone trajectory.
pub mod moving;

pub use moving::{MovingReflectionObserver, PreparedReflectedMovingReceiverSource};

#[derive(Clone, Copy, Debug)]
struct ImagePath {
    observer: MovingCompactObserver,
    sign: f64,
}

/// A fixed microphone and the COMPLETE image group of one to three orthogonal
/// ideal planes. This small fixed-size value is shared by all physical sources.
/// Walls and the microphone are stationary in the original fluid rest frame.
#[derive(Clone, Debug)]
pub struct IdealReflectionObserver {
    observer: MovingCompactObserver,
    planes: [IdealAcousticPlane; 3],
    plane_count: usize,
    paths: [ImagePath; 8],
}

impl IdealReflectionObserver {
    /// Admit the original microphone and all 2^m paths, including direct sound.
    /// Normals must be perpendicular within 64 machine epsilons; no geometry
    /// is snapped. The microphone may lie on any of the selected ideal walls.
    /// `maximum_source_paths` limits the number of paths PER physical source.
    /// Construction is bounded constant work with no allocation.
    pub fn new(
        observer: MovingCompactObserver,
        planes: &[IdealAcousticPlane],
        maximum_source_paths: usize,
    ) -> Result<Self, AeroacError> {
        let path_count = admit_planes(planes, maximum_source_paths)?;
        for plane in planes {
            if plane.distance(observer.observer_position_m())?.to_f64() < 0.0 {
                return Err(invalid(
                    "transient microphone is outside the fluid half-space",
                ));
            }
        }
        let direct = ImagePath {
            observer,
            sign: 1.0,
        };
        let mut paths = [direct; 8];
        for (mask, path) in paths.iter_mut().enumerate().take(path_count).skip(1) {
            let mut position = observer.observer_position_m();
            let mut sign = 1.0;
            // Invert the source-image composition, without relying on rounded
            // reflection matrices commuting bit-for-bit in an arbitrary frame.
            for (index, plane) in planes.iter().enumerate().rev() {
                if mask & (1 << index) != 0 {
                    let before = plane.distance(position)?.to_f64();
                    position = plane.reflected_position(position)?;
                    let after = plane.distance(position)?.to_f64();
                    if (before > 0.0 && after >= 0.0) || (before < 0.0 && after <= 0.0) {
                        return Err(invalid("virtual microphone is not representable"));
                    }
                    sign *= plane.boundary().sign();
                }
            }
            *path = ImagePath {
                observer: MovingCompactObserver::new(
                    observer.medium(),
                    position,
                    observer.limits(),
                )?,
                sign,
            };
        }
        let mut retained_planes = [planes[0]; 3];
        retained_planes[..planes.len()].copy_from_slice(planes);
        Ok(Self {
            observer,
            planes: retained_planes,
            plane_count: planes.len(),
            paths,
        })
    }

    /// Original physical microphone, fluid snapshot and motion/clearance limits.
    pub const fn observer(&self) -> &MovingCompactObserver {
        &self.observer
    }
    /// Original ideal boundaries, in their deterministic image-mask order.
    pub fn planes(&self) -> &[IdealAcousticPlane] {
        &self.planes[..self.plane_count]
    }
    /// All paths per physical source; index zero is the original direct path.
    pub const fn path_count(&self) -> usize {
        1 << self.plane_count
    }

    fn path(&self, index: usize) -> Result<ImagePath, AeroacError> {
        if index >= self.path_count() {
            return Err(invalid(
                "transient reflection path index is outside the image group",
            ));
        }
        Ok(self.paths[index])
    }
}

fn admit_planes(
    planes: &[IdealAcousticPlane],
    maximum_source_paths: usize,
) -> Result<usize, AeroacError> {
    if !(1..=3).contains(&planes.len()) {
        return Err(invalid(
            "transient reflection needs one to three ideal planes",
        ));
    }
    let path_count = 1_usize << planes.len();
    if path_count > maximum_source_paths {
        return Err(invalid(
            "transient reflection exceeds its complete path budget",
        ));
    }
    for (index, plane) in planes.iter().enumerate() {
        for other in &planes[..index] {
            if dot(plane.inward_normal(), other.inward_normal()).abs() > 64.0 * f64::EPSILON {
                return Err(invalid("transient reflection planes must be perpendicular"));
            }
        }
    }
    Ok(path_count)
}

/// Additional whole-trajectory admission work. Source storage and root caps
/// remain with the original trajectory. Reflection retains no heap allocation;
/// its preparation stack has at most 33 six-control position hulls.
#[derive(Clone, Copy, Debug)]
pub struct ReflectedTrajectoryBudget {
    /// Total Bernstein nodes across ALL original source spans. Each checks at
    /// most three half-spaces and eight virtual-microphone clearances.
    pub maximum_bound_nodes: usize,
    /// Maximum binary subdivision depth, inclusive, in 0..=32.
    pub maximum_subdivision_depth: usize,
}

/// One original finite C2 source observed through all admitted ideal-wall paths.
///
/// The original trajectory and shared reflection observer are borrowed, never
/// copied or advanced. Each path has its own finite arrival interval. Their
/// intersection is the complete receiver window: missing prehistory or tails
/// are refused, never treated as silence or extrapolated.
#[derive(Debug)]
pub struct PreparedReflectedSourceTrajectory<'a> {
    trajectory: &'a PreparedSourceTrajectory,
    observer: &'a IdealReflectionObserver,
    arrival_intervals: [[f64; 2]; 8],
    receiver_interval: [f64; 2],
}

impl<'a> PreparedReflectedSourceTrajectory<'a> {
    /// Admit every source span strictly inside all physical half-spaces and
    /// clear of every virtual microphone. Outward Bernstein hulls subdivide
    /// under the explicit work budget; inconclusive bounds refuse.
    pub fn new(
        trajectory: &'a PreparedSourceTrajectory,
        observer: &'a IdealReflectionObserver,
        budget: ReflectedTrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let physical = observer.observer();
        if trajectory.observer().medium() != physical.medium()
            || trajectory.observer().observer_position_m() != physical.observer_position_m()
            || trajectory.observer().limits() != physical.limits()
            || budget.maximum_bound_nodes < trajectory.segment_count()
            || budget.maximum_subdivision_depth > MAX_BOUND_DEPTH
        {
            return Err(invalid(
                "reflected trajectory observer, material or bound budget mismatch",
            ));
        }
        let first = trajectory.segments[0].first;
        let last = trajectory.segments[trajectory.segments.len() - 1].last;
        let mut arrival_intervals = [[0.0; 2]; 8];
        let mut receiver_interval = [f64::NEG_INFINITY, f64::INFINITY];
        for (index, interval) in arrival_intervals
            .iter_mut()
            .enumerate()
            .take(observer.path_count())
        {
            poll(cancelled)?;
            let microphone = observer.paths[index].observer;
            *interval = if index == 0 {
                trajectory.receiver_time_interval_s()
            } else {
                [
                    microphone
                        .prepare_emission(first.source_time_s, first.kinematics)?
                        .receiver_time_s(),
                    microphone
                        .prepare_emission(last.source_time_s, last.kinematics)?
                        .receiver_time_s(),
                ]
            };
            if interval[1] <= interval[0] {
                return Err(invalid(
                    "reflected trajectory arrivals are not representably increasing",
                ));
            }
            receiver_interval[0] = receiver_interval[0].max(interval[0]);
            receiver_interval[1] = receiver_interval[1].min(interval[1]);
        }
        if receiver_interval[0] > receiver_interval[1] {
            return Err(invalid(
                "retained trajectory has no complete reflected receiver window",
            ));
        }
        let mut visited = 0;
        for segment in &trajectory.segments {
            admit_span(observer, segment, budget, &mut visited, cancelled)?;
        }
        poll(cancelled)?;
        Ok(Self {
            trajectory,
            observer,
            arrival_intervals,
            receiver_interval,
        })
    }

    pub const fn trajectory(&self) -> &'a PreparedSourceTrajectory {
        self.trajectory
    }
    pub const fn reflection_observer(&self) -> &'a IdealReflectionObserver {
        self.observer
    }
    pub const fn path_count(&self) -> usize {
        self.observer.path_count()
    }
    /// Finite arrival window covered by EVERY retained image path.
    pub const fn receiver_time_interval_s(&self) -> [f64; 2] {
        self.receiver_interval
    }
    /// Path-specific coverage, useful for diagnosing missing reflected tails.
    pub fn path_receiver_time_interval_s(&self, path: usize) -> Result<[f64; 2], AeroacError> {
        self.observer.path(path)?;
        Ok(self.arrival_intervals[path])
    }
    /// The original whole-path bound covers every virtual receiver direction.
    pub const fn maximum_source_time_rate(&self) -> f64 {
        self.trajectory.maximum_source_time_rate()
    }
    /// Per-path interior refinements. Complete work multiplies by path count.
    pub const fn maximum_root_iterations(&self) -> usize {
        self.trajectory.maximum_root_iterations()
    }

    /// Prepare one physical propagation path at its OWN retarded source time.
    /// Apply ORIGINAL Q/F strengths; ideal-wall parity is already in the eight
    /// coefficients. Paths are ordered by plane bit mask, direct path first.
    /// Queries allocate nothing and share the original safeguarded root solver.
    pub fn prepare_path(
        &self,
        path: usize,
        receiver_time: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingCompactSample, AeroacError> {
        poll(cancelled)?;
        let image = self.observer.path(path)?;
        let [first, last] = self.arrival_intervals[path];
        if !receiver_time.is_finite()
            || receiver_time < first
            || receiver_time > last
            || !maximum_time_residual_s.is_finite()
            || maximum_time_residual_s < 0.0
        {
            return Err(invalid(
                "reflected query exceeds finite path history or residual domain",
            ));
        }
        if path == 0 {
            return self.trajectory.prepare_receiver(
                receiver_time,
                maximum_time_residual_s,
                cancelled,
            );
        }
        // Search original knots with this path's arrival map. A per-image
        // trajectory or full arrival array would multiply retained history.
        let mut begin = 0;
        let mut end = self.trajectory.segments.len();
        while begin < end {
            poll(cancelled)?;
            let middle = begin + (end - begin) / 2;
            let knot = self.trajectory.segments[middle].last;
            let arrival = image
                .observer
                .prepare_emission(knot.source_time_s, knot.kinematics)?
                .receiver_time_s();
            if arrival < receiver_time {
                begin = middle + 1;
            } else {
                end = middle;
            }
        }
        let segment = &self.trajectory.segments[begin];
        let first_arrival = image
            .observer
            .prepare_emission(segment.first.source_time_s, segment.first.kinematics)?
            .receiver_time_s();
        let last_arrival = image
            .observer
            .prepare_emission(segment.last.source_time_s, segment.last.kinematics)?
            .receiver_time_s();
        let prepared = self.trajectory.prepare_segment(
            &image.observer,
            segment,
            receiver_time,
            [first_arrival, last_arrival],
            maximum_time_residual_s,
            cancelled,
        )?;
        Ok(apply_parity(prepared, image.sign))
    }

    /// Prepare an original emission event for exact rounded-endpoint
    /// reconciliation. This does not clamp a retarded root or extend a history.
    pub fn prepare_path_emission(
        &self,
        path: usize,
        source_time: f64,
    ) -> Result<PreparedMovingCompactSample, AeroacError> {
        let image = self.observer.path(path)?;
        let motion = self.trajectory.kinematics_at(source_time)?;
        Ok(apply_parity(
            image.observer.prepare_emission(source_time, motion)?,
            image.sign,
        ))
    }

    /// Coherent pressure [Pa] from ONE physical history, queried separately at
    /// each path's retarded time. A history error aborts the complete result;
    /// no unsupported image is silently dropped. The callback must return the
    /// original fixed-frame strengths and derivatives at the requested time.
    /// Unrounded path pressures are retained until the final coherent sum.
    pub fn pressure_at(
        &self,
        receiver_time: f64,
        maximum_time_residual_s: f64,
        source_at: &mut impl FnMut(f64) -> Result<CompactSourceStrength, AeroacError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, AeroacError> {
        let [first, last] = self.receiver_interval;
        if !receiver_time.is_finite() || receiver_time < first || receiver_time > last {
            return Err(invalid("receiver time lacks complete reflected history"));
        }
        let mut pressure = Dd::ZERO;
        for path in 0..self.path_count() {
            let prepared =
                self.prepare_path(path, receiver_time, maximum_time_residual_s, cancelled)?;
            let strength = source_at(prepared.source_time_s())?;
            let (_, path_pressure) = prepared.evaluate_with_pressure_dd(strength)?;
            pressure = pressure + path_pressure;
        }
        poll(cancelled)?;
        let pressure = pressure.to_f64();
        if !pressure.is_finite() {
            return Err(invalid("coherent reflected pressure is unrepresentable"));
        }
        Ok(pressure)
    }
}

fn apply_parity(mut sample: PreparedMovingCompactSample, sign: f64) -> PreparedMovingCompactSample {
    if sign < 0.0 {
        sample.coefficients = sample.coefficients.map(|value| -value);
    }
    sample
}

#[derive(Clone, Copy)]
struct PositionNode {
    controls: [[Interval; 3]; 6],
    depth: usize,
}

fn admit_span(
    observer: &IdealReflectionObserver,
    segment: &QuinticSegment,
    budget: ReflectedTrajectoryBudget,
    visited: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AeroacError> {
    let zero = Interval::point(0.0);
    let mut controls = [[zero; 3]; 6];
    for (control, original) in controls.iter_mut().zip(segment.position) {
        for axis in 0..3 {
            control[axis] = enclose_dd(original[axis])?;
        }
    }
    // Exact accepted knots are returned verbatim by the original evaluator.
    for axis in 0..3 {
        controls[5][axis] = controls[5][axis].hull(
            Interval::point(segment.last.kinematics.position_m[axis])
                - Interval::point(segment.first.kinematics.position_m[axis]),
        );
    }
    for knot in [segment.first, segment.last] {
        for plane in observer.planes() {
            if plane.distance(knot.kinematics.position_m)?.to_f64() <= 0.0 {
                return Err(invalid(
                    "source trajectory must remain strictly inside every half-space",
                ));
            }
        }
    }
    let origin = segment.first.kinematics.position_m;
    let mut stack = [PositionNode { controls, depth: 0 }; MAX_BOUND_DEPTH + 1];
    let mut pending = 1;
    while pending > 0 {
        poll(cancelled)?;
        if *visited >= budget.maximum_bound_nodes {
            return Err(AeroacError::NotConverged {
                what: "reflected trajectory bound-node budget",
                residual: 1.0,
            });
        }
        *visited += 1;
        pending -= 1;
        let node = stack[pending];
        let inside = observer.planes().iter().all(|plane| {
            node.controls.iter().all(|control| {
                let mut distance = zero;
                for axis in 0..3 {
                    // Enclose the representable world position before taking
                    // a wall distance. Cancelling a large common origin first
                    // could hide a trajectory that rounds onto the wall.
                    let position = Interval::point(origin[axis]) + control[axis];
                    distance = distance
                        + (position - Interval::point(plane.point_m()[axis]))
                            * Interval::point(plane.inward_normal()[axis]);
                }
                distance.lo() > 0.0
            })
        });
        let clear = observer.paths[..observer.path_count()].iter().all(|path| {
            let relative = core::array::from_fn(|i| {
                Interval::point(path.observer.observer_position_m()[i]) - Interval::point(origin[i])
            });
            position_hull_clear(
                &node.controls,
                relative,
                observer.observer.limits().minimum_distance_m,
            )
        });
        if inside && clear {
            continue;
        }
        if node.depth >= budget.maximum_subdivision_depth {
            return Err(AeroacError::NotConverged {
                what: "reflected source half-space/clearance hull remains inadmissible",
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
    Ok(())
}
