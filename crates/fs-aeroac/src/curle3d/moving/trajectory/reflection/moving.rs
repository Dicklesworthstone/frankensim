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
use super::{
    PositionNode, ReflectedCompactField, ReflectedFieldSum, ReflectedTrajectoryBudget,
    admit_planes, apply_image, image_velocity,
};
use crate::AeroacError;
use crate::curle3d::IdealAcousticPlane;
use crate::curle3d::array::wedge::{WedgeImages, WedgeReceiverImage};
use crate::curle3d::moving::{
    CompactSourceStrength, MovingCompactObserver, PreparedMovingCompactSample, dot, norm,
};
use fs_ivl::Interval;
use fs_math::dd::Dd;

/// One original C2 microphone path and one complete orthogonal or pi/m wedge
/// ideal-wall group. Every source in a coherent render borrows this same owner.
///
/// The walls are stationary in the homogeneous lossless fluid-rest frame.
/// Uniform-flow constructors admit laboratory walls only when the original
/// mean flow is tangent to every wall. The receiver may reach a wall when its complete
/// interpolant can be proved inside; an inconclusive rounded hull refuses.
#[derive(Debug)]
pub struct MovingReflectionObserver<'a> {
    receiver: &'a PreparedKinematicTrajectory,
    planes: [IdealAcousticPlane; 3],
    plane_count: usize,
    signs: Vec<f64>,
    wedge_images: Vec<WedgeReceiverImage>,
    wedge_order: Option<usize>,
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
        Self::orthogonal(receiver, planes, path_count, budget, cancelled)
    }

    /// Reflect accepted laboratory motion in a homogeneous tangential flow.
    ///
    /// `planes` are stationary LABORATORY walls. They are translated into the
    /// original trajectory's fluid frame once, at preparation. For n.U = 0,
    /// n.(x-a) = n.(x'-a'), so the same stationary image group applies at every
    /// time. A normal flow component refuses instead of freezing a moving wall.
    /// Source and microphone must retain the exact same `UniformFlowFrame`.
    pub fn in_uniform_flow(
        receiver: &'a PreparedKinematicTrajectory,
        planes: &[IdealAcousticPlane],
        maximum_source_paths: usize,
        budget: ReflectedTrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let path_count = admit_planes(planes, maximum_source_paths)?;
        let transformed = fluid_frame_planes(receiver, planes)?;
        Self::orthogonal(
            receiver,
            &transformed[..planes.len()],
            path_count,
            budget,
            cancelled,
        )
    }

    fn orthogonal(
        receiver: &'a PreparedKinematicTrajectory,
        planes: &[IdealAcousticPlane],
        path_count: usize,
        budget: ReflectedTrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        budget.admit_payload(
            path_count.checked_mul(size_of::<f64>())
                .and_then(|bytes| bytes.checked_add(size_of::<Self>())),
        )?;
        admit_halfspaces(receiver, planes, false, budget, cancelled)?;
        let mut retained = [planes[0]; 3];
        retained[..planes.len()].copy_from_slice(planes);
        let mut signs = Vec::new();
        signs
            .try_reserve_exact(path_count)
            .map_err(|_| invalid("moving image allocation exceeds available memory"))?;
        signs.resize(path_count, 1.0);
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
            wedge_images: Vec::new(),
            wedge_order: None,
        })
    }

    /// Prepare all 2m inverse images of the SAME moving microphone in an
    /// admitted ideal pi/m wedge. Geometry, boundary signs and mixed-wall
    /// admissibility share the harmonic owner. Each affine image is cached
    /// once; no query repeats a sequence of m rotations or copies a history.
    /// The complete physical receiver trajectory must remain in both walls.
    pub fn new_wedge(
        receiver: &'a PreparedKinematicTrajectory,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
        budget: ReflectedTrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let images = WedgeImages::new(planes, order)?;
        if images.path_count() > maximum_source_paths {
            return Err(invalid("moving wedge exceeds its complete path budget"));
        }
        if receiver.flow_frame().is_some() {
            return Err(invalid(
                "moving wall images require paths in the stationary fluid frame",
            ));
        }
        Self::wedge(receiver, planes, images, order, budget, cancelled)
    }

    /// Complete pi/m wedge images in a flow parallel to the wedge edge.
    /// Walls are authored in laboratory coordinates; all original source and
    /// receiver histories, budgets and per-image retarded times are reused.
    pub fn wedge_in_uniform_flow(
        receiver: &'a PreparedKinematicTrajectory,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
        budget: ReflectedTrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let transformed = fluid_frame_planes(receiver, &planes)?;
        let planes = [transformed[0], transformed[1]];
        let images = WedgeImages::new(planes, order)?;
        if images.path_count() > maximum_source_paths {
            return Err(invalid("moving wedge exceeds its complete path budget"));
        }
        Self::wedge(receiver, planes, images, order, budget, cancelled)
    }

    fn wedge(
        receiver: &'a PreparedKinematicTrajectory,
        planes: [IdealAcousticPlane; 2],
        images: WedgeImages,
        order: usize,
        budget: ReflectedTrajectoryBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        budget.admit_payload(
            images.path_count().checked_mul(size_of::<f64>() + size_of::<WedgeReceiverImage>())
                .and_then(|bytes| bytes.checked_add(size_of::<Self>())),
        )?;
        admit_halfspaces(receiver, &planes, false, budget, cancelled)?;
        let mut signs = Vec::new();
        let mut wedge_images = Vec::new();
        signs
            .try_reserve_exact(images.path_count())
            .map_err(|_| invalid("moving wedge allocation exceeds available memory"))?;
        wedge_images
            .try_reserve_exact(images.path_count())
            .map_err(|_| invalid("moving wedge allocation exceeds available memory"))?;
        images.visit_inverse(|_, image, sign| {
            poll(cancelled)?;
            wedge_images.push(image);
            signs.push(sign);
            Ok(())
        })?;
        poll(cancelled)?;
        Ok(Self {
            receiver,
            planes: [planes[0], planes[1], planes[0]],
            plane_count: 2,
            signs,
            wedge_images,
            wedge_order: Some(order),
        })
    }

    /// Original microphone trajectory, without image copies or clock shifts.
    pub const fn receiver(&self) -> &'a PreparedKinematicTrajectory {
        self.receiver
    }
    /// Admitted walls in the original trajectories' FLUID frame.
    pub fn planes(&self) -> &[IdealAcousticPlane] {
        &self.planes[..self.plane_count]
    }
    /// Complete path count per physical source, including direct path zero.
    pub const fn path_count(&self) -> usize {
        self.signs.len()
    }
    /// Dihedral order when prepared by `new_wedge`, otherwise orthogonal planes.
    pub const fn wedge_order(&self) -> Option<usize> {
        self.wedge_order
    }
}

// A plane's tangential displacement under x' = x - origin - U*(t-epoch)
// does not change its locus. Retain the exact stored normal; normalizing it a
// second time could change a wall's rounded geometry or wedge admission.
fn fluid_frame_planes(
    receiver: &PreparedKinematicTrajectory,
    planes: &[IdealAcousticPlane],
) -> Result<[IdealAcousticPlane; 3], AeroacError> {
    let frame = receiver
        .flow_frame()
        .ok_or_else(|| invalid("uniform-flow reflections need a declared trajectory frame"))?;
    if !(1..=3).contains(&planes.len()) {
        return Err(invalid("uniform-flow reflection plane count"));
    }
    let velocity = frame.velocity_m_s();
    let scale = velocity.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
    let direction = if scale == 0.0 {
        [0.0; 3]
    } else {
        let scaled = velocity.map(|v| v / scale);
        let length = norm(scaled);
        scaled.map(|v| v / length)
    };
    let mut transformed = [planes[0]; 3];
    for (target, original) in transformed.iter_mut().zip(planes) {
        if dot(direction, original.inward_normal()).abs() > 64.0 * f64::EPSILON {
            return Err(invalid(
                "uniform flow must be tangent to every reflecting wall",
            ));
        }
        let point = frame.to_fluid_frame(
            frame.reference_time_s(),
            crate::curle3d::moving::SourceKinematics {
                position_m: original.point_m(),
                ..crate::curle3d::moving::SourceKinematics::default()
            },
        )?;
        // Reuse the exact normal via a translated copy owned by its plane type.
        *target = original.translated_to(point.position_m)?;
    }
    Ok(transformed)
}

/// All ideal-wall paths from one original source to one moving microphone.
///
/// Bounded storage retains one borrowed pair view per image. Source
/// and receiver histories remain singular; pressure queries allocate nothing.
/// Each path has its own retarded time and observed Doppler
/// rate. Missing source coverage refuses the COMPLETE observation window.
#[derive(Debug)]
pub struct PreparedReflectedMovingReceiverSource<'a> {
    observer: &'a MovingReflectionObserver<'a>,
    paths: Vec<PreparedMovingReceiverSource<'a>>,
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
        if source.flow_frame() != observer.receiver.flow_frame() {
            return Err(invalid(
                "reflected source and microphone must share their original fluid frame",
            ));
        }
        source_bounds.admit_payload(
            observer.path_count().checked_mul(size_of::<PreparedMovingReceiverSource<'a>>())
                .and_then(|bytes| bytes.checked_add(size_of::<Self>())),
        )?;
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
        let mut paths = Vec::new();
        paths
            .try_reserve_exact(observer.path_count())
            .map_err(|_| invalid("reflected pair allocation exceeds available memory"))?;
        paths.push(direct);
        let mut maximum_source_time_rate = direct.maximum_source_time_rate();
        for mask in 1..observer.path_count() {
            poll(cancelled)?;
            let image = if observer.wedge_order.is_some() {
                ReceiverImage::Wedge(&observer.wedge_images[mask])
            } else {
                ReceiverImage::Orthogonal {
                    planes: observer.planes(),
                    mask,
                }
            };
            let path = PreparedMovingReceiverSource::new_image(
                source,
                observer.receiver,
                Some(image),
                config,
                pair_budget,
                &mut work,
                cancelled,
            )?;
            maximum_source_time_rate =
                maximum_source_time_rate.max(path.maximum_source_time_rate());
            paths.push(path);
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
    pub fn source(&self) -> &'a PreparedKinematicTrajectory {
        self.paths[0].source()
    }
    pub const fn receiver(&self) -> &'a PreparedKinematicTrajectory {
        self.observer.receiver
    }
    pub fn config(&self) -> MovingReceiverConfig {
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
    pub fn maximum_root_iterations(&self) -> usize {
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
    /// includes rigid/release parity; monopole velocity follows the transpose
    /// image map into the physical receiver axes. Receiver kinematics remain physical;
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
        sample.pressure = apply_image(
            sample.pressure,
            self.observer.signs[path],
            self.observer.planes(),
            path,
            self.observer
                .wedge_images
                .get(path)
                .filter(|_| path != 0)
                .map(|image| &image.matrix),
        );
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
        Ok(apply_image(
            observer.prepare_emission(source_time, self.source().kinematics_at(source_time)?)?,
            self.observer.signs[path],
            self.observer.planes(),
            path,
            self.observer
                .wedge_images
                .get(path)
                .filter(|_| path != 0)
                .map(|image| &image.matrix),
        ))
    }

    /// Coherent pressure from one physical Q/F history, queried on its own
    /// clock at each path's retarded time. A failed history/path aborts the sum.
    /// Unrounded path pressures are retained until the final coherent sum.
    pub fn pressure_at(
        &self,
        receiver_time: f64,
        source_at: &mut impl FnMut(f64) -> Result<CompactSourceStrength, AeroacError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, AeroacError> {
        let mut pressure = Dd::ZERO;
        for path in 0..self.path_count() {
            let sample = self.prepare_path(path, receiver_time, cancelled)?.pressure;
            let strength = source_at(sample.source_time_s())?;
            let (_, path_pressure) = sample.evaluate_with_pressure_dd(strength)?;
            pressure = pressure + path_pressure;
        }
        poll(cancelled)?;
        let value = pressure.to_f64();
        if !value.is_finite() {
            return Err(invalid("moving reflected pressure is unrepresentable"));
        }
        Ok(value)
    }

    /// Full coherent Q/F field at the microphone's original arrival event.
    ///
    /// Each image uses its own causal source time and its virtual microphone
    /// position at `receiver_time`. The history callback receives that path
    /// index, virtual observer and source time. Supply the ORIGINAL force
    /// integral `(3 n (n.F) - F)/R^3` [N s/m³] at that fixed virtual position,
    /// using the original source history in the fluid-rest frame. A history
    /// computed at another path or a previous microphone event is invalid.
    /// The observer's `force_history_integrand` supplies the native kernel.
    ///
    /// Velocity is the fluid's acoustic perturbation in physical receiver
    /// axes; microphone velocity is not added. The result retains near fields
    /// and cross-path pressure-work flux, including in a tangential uniform
    /// flow. Mean-flow advection of acoustic energy is not part of `p*u`.
    /// The caller owns explicit prehistory, quadrature accuracy and continuous
    /// past-worldline clearance. Any failure aborts the complete field.
    pub fn force_field_at(
        &self,
        receiver_time: f64,
        source_at: &mut impl FnMut(f64) -> Result<CompactSourceStrength, AeroacError>,
        force_history_at: &mut impl FnMut(
            usize,
            &MovingCompactObserver,
            f64,
        ) -> Result<[f64; 3], AeroacError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReflectedCompactField, AeroacError> {
        poll(cancelled)?;
        let mut sum = ReflectedFieldSum::default();
        for path in 0..self.path_count() {
            let original = self
                .path(path)?
                .prepare_receiver(receiver_time, cancelled)?
                .pressure;
            let source_time = original.source_time_s();
            let residual = original.retarded_time_residual_s();
            let virtual_observer = self.path(path)?.observer_at(receiver_time)?;
            let prepared = virtual_observer
                .prepare_force_sample(original, self.source().kinematics_at(source_time)?)?;
            poll(cancelled)?;
            let strength = source_at(source_time)?;
            poll(cancelled)?;
            let history = force_history_at(path, &virtual_observer, source_time)?;
            let (_, pressure, velocity) = prepared.evaluate_with_fields_dd(strength, history)?;
            let sign = self.observer.signs[path];
            sum.add(
                pressure * Dd::from_f64(sign),
                image_velocity(
                    velocity,
                    sign,
                    self.observer.planes(),
                    path,
                    self.observer
                        .wedge_images
                        .get(path)
                        .filter(|_| path != 0)
                        .map(|image| &image.matrix),
                ),
                residual,
            );
        }
        poll(cancelled)?;
        sum.finish(receiver_time)
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
