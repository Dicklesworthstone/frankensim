//! Retarded ideal-wall reflection on an original moving source history.
//!
//! For a stationary isometric image transform H, the Euclidean Green function
//! obeys G(x,H y)=G(H^-1 x,y). Using a virtual receiver H^-1 x therefore keeps
//! the ORIGINAL source position, velocity, acceleration, force and force rate.
//! It is exactly equivalent to transforming every image-source vector, without
//! cloning or reinterpolating the trajectory. Rigid/release parity multiplies
//! all eight original pressure coefficients. The monopole velocity gradient
//! also follows the transpose image map back to the physical receiver axes.
//! Each path keeps its own retarded time; a reflected transient is not a gain
//! applied to the direct waveform.
//!
//! The medium is homogeneous, lossless and at rest in the stationary walls'
//! frame. This is the finite ideal-plane image construction of Allen & Berkley,
//! JASA 65 (1979), 943–950, doi:10.1121/1.382599, with the moving Green-function
//! derivatives owned by the original compact source operator. Finite impedance,
//! general wedge angles, moving walls, dispersion and feedback are absent. Moving
//! microphones use the same image law through [`moving`].

use super::{
    BoundNode, MAX_BOUND_DEPTH, PreparedSourceTrajectory, QuinticSegment, doppler, dot, enclose_dd,
    invalid, poll, position_hull_clear, split_controls,
};
use crate::AeroacError;
use crate::curle3d::IdealAcousticPlane;
use crate::curle3d::array::wedge::WedgeImages;
use crate::curle3d::moving::{
    CompactSourceStrength, MovingCompactObserver, PreparedMovingCompactSample,
};
use fs_ivl::Interval;
use fs_math::dd::Dd;

/// Reflected observation on one original moving microphone trajectory.
pub mod moving;

pub use moving::{MovingReflectionObserver, PreparedReflectedMovingReceiverSource};

/// Coherent instantaneous Q/F field of one complete ideal-wall image group.
///
/// Velocity is the acoustic perturbation in the original fluid-rest axes,
/// including force induction and past forcing. It is not microphone velocity.
/// The signed pressure-work flux is formed AFTER summing all path fields;
/// neither a harmonic half-product nor a plane-wave approximation is used.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReflectedCompactField {
    /// Common physical microphone event [s].
    pub receiver_time_s: f64,
    /// Coherent total pressure [Pa].
    pub pressure_pa: f64,
    /// Coherent total Cartesian acoustic particle velocity [m/s].
    pub particle_velocity_m_s: [f64; 3],
    /// Signed instantaneous pressure-work flux [W/m²].
    pub instantaneous_intensity_w_m2: [f64; 3],
    /// Largest absolute retarded-time residual over the complete image group.
    pub maximum_time_residual_s: f64,
}

#[derive(Default)]
struct ReflectedFieldSum {
    pressure: Dd,
    velocity: [Dd; 3],
    maximum_time_residual_s: f64,
}

impl ReflectedFieldSum {
    fn add(&mut self, pressure: Dd, velocity: [Dd; 3], residual: f64) {
        self.pressure = self.pressure + pressure;
        for (sum, value) in self.velocity.iter_mut().zip(velocity) {
            *sum = *sum + value;
        }
        self.maximum_time_residual_s = self.maximum_time_residual_s.max(residual.abs());
    }

    fn finish(self, receiver_time_s: f64) -> Result<ReflectedCompactField, AeroacError> {
        let pressure_pa = self.pressure.to_f64();
        let particle_velocity_m_s = self.velocity.map(Dd::to_f64);
        let instantaneous_intensity_w_m2 = self.velocity.map(|u| (self.pressure * u).to_f64());
        if !pressure_pa.is_finite()
            || particle_velocity_m_s
                .iter()
                .chain(&instantaneous_intensity_w_m2)
                .any(|value| !value.is_finite())
        {
            return Err(invalid("coherent reflected field is unrepresentable"));
        }
        Ok(ReflectedCompactField {
            receiver_time_s,
            pressure_pa,
            particle_velocity_m_s,
            instantaneous_intensity_w_m2,
            maximum_time_residual_s: self.maximum_time_residual_s,
        })
    }
}

#[derive(Clone, Copy, Debug)]
struct ImagePath {
    observer: MovingCompactObserver,
    sign: f64,
    inverse_jacobian: Option<[[f64; 3]; 3]>,
}

/// A fixed microphone and the COMPLETE image group of one to three orthogonal
/// ideal planes, or the two walls of an admitted pi/m wedge.
/// Prepared path geometry is shared by all physical sources.
/// Walls and the microphone are stationary in the original fluid rest frame.
#[derive(Clone, Debug)]
pub struct IdealReflectionObserver {
    observer: MovingCompactObserver,
    planes: [IdealAcousticPlane; 3],
    plane_count: usize,
    paths: Vec<ImagePath>,
    wedge_order: Option<usize>,
}

impl IdealReflectionObserver {
    /// Admit the original microphone and all 2^m paths, including direct sound.
    /// Normals must be perpendicular within 64 machine epsilons; no geometry
    /// is snapped. The microphone may lie on any of the selected ideal walls.
    /// `maximum_source_paths` limits the number of paths PER physical source.
    /// Construction retains only the bounded image geometry, not trajectories.
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
            inverse_jacobian: None,
        };
        let mut paths = Vec::new();
        paths
            .try_reserve_exact(path_count)
            .map_err(|_| invalid("transient image allocation exceeds available memory"))?;
        paths.resize(path_count, direct);
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
                inverse_jacobian: None,
            };
        }
        let mut retained_planes = [planes[0]; 3];
        retained_planes[..planes.len()].copy_from_slice(planes);
        Ok(Self {
            observer,
            planes: retained_planes,
            plane_count: planes.len(),
            paths,
            wedge_order: None,
        })
    }

    /// Retain every exact dihedral image of an ideal pi/m wedge, m in 2..=128.
    /// Equal wall conditions admit every supported m; mixed rigid/release walls
    /// require even m. The harmonic owner's geometry and group admission are
    /// reused, with no angle snapping or finite reflection-order truncation.
    /// Direct path zero preserves its original arithmetic. Every other path
    /// has its own retarded time on the same source history.
    pub fn new_wedge(
        observer: MovingCompactObserver,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
    ) -> Result<Self, AeroacError> {
        let images = WedgeImages::new(planes, order)?;
        if images.path_count() > maximum_source_paths {
            return Err(invalid("transient wedge exceeds its complete path budget"));
        }
        for plane in planes {
            if plane.distance(observer.observer_position_m())?.to_f64() < 0.0 {
                return Err(invalid(
                    "transient microphone is outside the acoustic wedge",
                ));
            }
        }
        let mut paths = Vec::new();
        paths
            .try_reserve_exact(images.path_count())
            .map_err(|_| invalid("transient wedge allocation exceeds available memory"))?;
        images.visit_inverse(|index, image, sign| {
            let microphone = if index == 0 {
                observer
            } else {
                MovingCompactObserver::new(
                    observer.medium(),
                    image.position(observer.observer_position_m())?,
                    observer.limits(),
                )?
            };
            paths.push(ImagePath {
                observer: microphone,
                sign,
                inverse_jacobian: (index != 0).then_some(image.matrix),
            });
            Ok(())
        })?;
        Ok(Self {
            observer,
            planes: [planes[0], planes[1], planes[0]],
            plane_count: 2,
            paths,
            wedge_order: Some(order),
        })
    }

    /// Original physical microphone, fluid snapshot and motion/clearance limits.
    pub const fn observer(&self) -> &MovingCompactObserver {
        &self.observer
    }
    /// Original ideal boundaries, in the supplied physical order.
    pub fn planes(&self) -> &[IdealAcousticPlane] {
        &self.planes[..self.plane_count]
    }
    /// All paths per physical source; index zero is the original direct path.
    pub const fn path_count(&self) -> usize {
        self.paths.len()
    }
    /// Dihedral order when prepared by `new_wedge`, otherwise orthogonal planes.
    pub const fn wedge_order(&self) -> Option<usize> {
        self.wedge_order
    }

    /// Geometry and vector frame for one stationary source's reflected field.
    ///
    /// The original source must lie strictly inside every physical half-space.
    /// Returns the existing virtual receiver, its pressure parity, and a
    /// row-major matrix mapping virtual particle velocity to the physical frame.
    /// The velocity matrix already includes parity. It is the transpose of the
    /// receiver-image Jacobian. Orthogonal paths apply vector reflections in
    /// FORWARD plane order; wedge paths use the cached inverse map's transpose.
    /// This retains the actual rounded transforms without assuming they
    /// commute. No source history, propagation phase or field is evaluated.
    pub fn stationary_field_path(
        &self,
        path: usize,
        source_m: [f64; 3],
    ) -> Result<(MovingCompactObserver, f64, [[f64; 3]; 3]), AeroacError> {
        let image = self.path(path)?;
        for plane in self.planes() {
            if plane.distance(source_m)?.to_f64() <= 0.0 {
                return Err(invalid(
                    "stationary source must lie strictly inside the fluid half-spaces",
                ));
            }
        }
        if let Some(jacobian) = image.inverse_jacobian {
            let velocity = core::array::from_fn(|row| {
                core::array::from_fn(|column| {
                    (Dd::from_f64(image.sign) * Dd::from_f64(jacobian[column][row])).to_f64()
                })
            });
            return Ok((image.observer, image.sign, velocity));
        }
        let columns: [[Dd; 3]; 3] = core::array::from_fn(|column| {
            let mut vector = [Dd::ZERO; 3];
            vector[column] = Dd::ONE;
            for (index, plane) in self.planes().iter().enumerate() {
                if path & (1 << index) == 0 {
                    continue;
                }
                let normal = plane.inward_normal();
                let projection = (0..3).fold(Dd::ZERO, |sum, axis| {
                    sum + vector[axis] * Dd::from_f64(normal[axis])
                });
                vector = core::array::from_fn(|axis| {
                    vector[axis] - projection * Dd::from_f64(2.0 * normal[axis])
                });
            }
            vector
        });
        let velocity = core::array::from_fn(|row| {
            core::array::from_fn(|column| {
                (columns[column][row] * Dd::from_f64(image.sign)).to_f64()
            })
        });
        Ok((image.observer, image.sign, velocity))
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

/// Additional whole-trajectory admission work and retained image bounds.
/// Source polynomials and root caps remain with the original trajectory.
/// Each fixed-microphone source caches one rate per ORIGINAL span, taking the
/// maximum over its image group rather than allocating a span-by-image matrix.
#[derive(Clone, Copy, Debug)]
pub struct ReflectedTrajectoryBudget {
    /// Logical owner and owned arrays: image metadata, borrowed pair records or
    /// fixed-microphone arrival/rate caches. Borrowed trajectories, allocator
    /// metadata/rounding and fixed-depth preparation stacks are excluded.
    pub maximum_payload_bytes: usize,
    /// Total Bernstein nodes across ALL original source spans. Each checks at
    /// most three half-spaces and the complete admitted virtual-microphone set.
    pub maximum_bound_nodes: usize,
    /// Maximum binary subdivision depth, inclusive, in 0..=32.
    pub maximum_subdivision_depth: usize,
}

impl ReflectedTrajectoryBudget {
    fn admit_payload(self, bytes: Option<usize>) -> Result<(), AeroacError> {
        let bytes = bytes.ok_or_else(|| invalid("reflected trajectory payload overflows"))?;
        if bytes > self.maximum_payload_bytes {
            return Err(invalid("reflected trajectory exceeds retained payload budget"));
        }
        Ok(())
    }
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
    arrival_intervals: Vec<[f64; 2]>,
    span_maximum_rates: Vec<f64>,
    maximum_source_time_rate: f64,
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
        budget.admit_payload(
            observer.path_count().checked_mul(size_of::<[f64; 2]>())
                .and_then(|bytes| trajectory.segment_count().checked_mul(size_of::<f64>())
                    .and_then(|rates| bytes.checked_add(rates)))
                .and_then(|bytes| bytes.checked_add(size_of::<Self>())),
        )?;
        let first = trajectory.segments[0].polynomial.first;
        let last = trajectory.segments[trajectory.segments.len() - 1].polynomial.last;
        let mut arrival_intervals = Vec::new();
        arrival_intervals
            .try_reserve_exact(observer.path_count())
            .map_err(|_| invalid("reflected arrival allocation exceeds available memory"))?;
        arrival_intervals.resize(observer.path_count(), [0.0; 2]);
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
        let mut span_maximum_rates = Vec::new();
        span_maximum_rates.try_reserve_exact(trajectory.segment_count())
            .map_err(|_| invalid("reflected Doppler cache allocation exceeds available memory"))?;
        let mut visited = 0;
        let mut maximum_source_time_rate = 0.0_f64;
        for segment in &trajectory.segments {
            let rate = admit_span(
                observer,
                &segment.polynomial,
                trajectory.maximum_mach,
                budget,
                &mut visited,
                cancelled,
            )?;
            maximum_source_time_rate = maximum_source_time_rate.max(rate);
            span_maximum_rates.push(rate);
        }
        poll(cancelled)?;
        Ok(Self {
            trajectory,
            observer,
            arrival_intervals,
            span_maximum_rates,
            maximum_source_time_rate,
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
    /// Directional ceiling covering every ORIGINAL span and virtual receiver.
    /// A direct-path redshift is never substituted for an image's direction.
    pub const fn maximum_source_time_rate(&self) -> f64 {
        self.maximum_source_time_rate
    }

    /// Directional ceiling for a finite receiver window shared by all images.
    ///
    /// Outward endpoint searches bracket the union of causally possible source
    /// spans. One cached max-over-images value per original span bounds that
    /// union without a span-by-image matrix. Partial spans keep their complete
    /// bound. Work is O(paths * log N + retained spans), with no allocation.
    pub fn maximum_source_time_rate_in_window(
        &self,
        receiver_time_interval_s: [f64; 2],
        maximum_time_residual_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, AeroacError> {
        poll(cancelled)?;
        doppler::admit_window(
            receiver_time_interval_s,
            self.receiver_interval,
            maximum_time_residual_s,
        )?;
        let mut first = self.trajectory.segment_count();
        let mut last = 0;
        for path in &self.observer.paths {
            let range = doppler::source_spans(
                self.trajectory,
                &path.observer,
                receiver_time_interval_s,
                maximum_time_residual_s,
                cancelled,
            )?;
            if !range.is_empty() {
                first = first.min(range.start);
                last = last.max(range.end);
            }
        }
        let mut maximum = 0.0_f64;
        if first < last {
            for rate in &self.span_maximum_rates[first..last] {
                poll(cancelled)?;
                maximum = maximum.max(*rate);
            }
        }
        Ok(if maximum > 0.0 { maximum } else { self.maximum_source_time_rate })
    }
    /// Per-path interior refinements. Complete work multiplies by path count.
    pub const fn maximum_root_iterations(&self) -> usize {
        self.trajectory.maximum_root_iterations()
    }

    /// Prepare one physical propagation path at its OWN retarded source time.
    /// Apply ORIGINAL Q/F strengths; ideal-wall parity is already in the eight
    /// pressure coefficients and the monopole velocity uses physical receiver
    /// axes. Orthogonal paths use plane bit masks; wedge paths use
    /// R^k, R^k S1 with R=S2 S1. The direct path is always first.
    /// Queries allocate nothing and share the original safeguarded root solver.
    pub fn prepare_path(
        &self,
        path: usize,
        receiver_time: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingCompactSample, AeroacError> {
        let prepared =
            self.prepare_original_path(path, receiver_time, maximum_time_residual_s, cancelled)?;
        let image = self.observer.path(path)?;
        Ok(apply_image(
            prepared,
            image.sign,
            self.observer.planes(),
            path,
            image.inverse_jacobian.as_ref(),
        ))
    }

    // Original source axes at the virtual receiver, before image parity or
    // vector rotation. Pressure and force-field queries share this exact root.
    fn prepare_original_path(
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
            let knot = self.trajectory.segments[middle].polynomial.last;
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
        let segment = &self.trajectory.segments[begin].polynomial;
        let first_arrival = image
            .observer
            .prepare_emission(segment.first.source_time_s, segment.first.kinematics)?
            .receiver_time_s();
        let last_arrival = image
            .observer
            .prepare_emission(segment.last.source_time_s, segment.last.kinematics)?
            .receiver_time_s();
        self.trajectory.prepare_segment(
            &image.observer,
            segment,
            receiver_time,
            [first_arrival, last_arrival],
            maximum_time_residual_s,
            cancelled,
        )
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
        Ok(apply_image(
            image.observer.prepare_emission(source_time, motion)?,
            image.sign,
            self.observer.planes(),
            path,
            image.inverse_jacobian.as_ref(),
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

    /// Full coherent volume/force field at one physical microphone event.
    ///
    /// `source_at` supplies ORIGINAL fixed-frame Q/F strengths and derivatives
    /// at each path's own retarded time. `force_history_at` receives the path
    /// index, its VIRTUAL observer and that same source time. It must return
    /// `integral (3 n (n.F) - F)/R^3 ds` [N s/m³] on the original source clock,
    /// using original source positions and forces at this virtual observer.
    /// [`MovingCompactObserver::force_history_integrand`] supplies the kernel.
    /// A direct-path integral cannot be reused for a different virtual observer.
    ///
    /// Explicit prehistory, integration accuracy and continuous separation of
    /// every virtual observer from the complete past worldline remain the
    /// history owner's obligations. Zero history declares zero prior forcing.
    /// A path, history or cancellation error aborts the entire result. Queries
    /// allocate nothing; each path solves its original causal root only once.
    pub fn force_field_at(
        &self,
        receiver_time: f64,
        maximum_time_residual_s: f64,
        source_at: &mut impl FnMut(f64) -> Result<CompactSourceStrength, AeroacError>,
        force_history_at: &mut impl FnMut(
            usize,
            &MovingCompactObserver,
            f64,
        ) -> Result<[f64; 3], AeroacError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReflectedCompactField, AeroacError> {
        poll(cancelled)?;
        let [first, last] = self.receiver_interval;
        if !receiver_time.is_finite() || receiver_time < first || receiver_time > last {
            return Err(invalid("receiver time lacks complete reflected history"));
        }
        let mut sum = ReflectedFieldSum::default();
        for path in 0..self.path_count() {
            let original = self.prepare_original_path(
                path,
                receiver_time,
                maximum_time_residual_s,
                cancelled,
            )?;
            let source_time = original.source_time_s();
            let residual = original.retarded_time_residual_s();
            let image = self.observer.path(path)?;
            let prepared = image
                .observer
                .prepare_force_sample(original, self.trajectory.kinematics_at(source_time)?)?;
            poll(cancelled)?;
            let strength = source_at(source_time)?;
            poll(cancelled)?;
            let history = force_history_at(path, &image.observer, source_time)?;
            let (_, pressure, velocity) = prepared.evaluate_with_fields_dd(strength, history)?;
            sum.add(
                pressure * Dd::from_f64(image.sign),
                image_velocity(
                    velocity,
                    image.sign,
                    self.observer.planes(),
                    path,
                    image.inverse_jacobian.as_ref(),
                ),
                residual,
            );
        }
        poll(cancelled)?;
        sum.finish(receiver_time)
    }
}

// Velocity is a polar vector: parity times the transpose of the virtual
// receiver's image Jacobian. Keep DD fields until AFTER coherent summation.
fn image_velocity(
    mut velocity: [Dd; 3],
    sign: f64,
    planes: &[IdealAcousticPlane],
    mask: usize,
    inverse_jacobian: Option<&[[f64; 3]; 3]>,
) -> [Dd; 3] {
    if sign < 0.0 {
        velocity = velocity.map(|value| -value);
    }
    if let Some(jacobian) = inverse_jacobian {
        return core::array::from_fn(|row| {
            (0..3).fold(Dd::ZERO, |sum, axis| {
                sum + Dd::from_f64(jacobian[axis][row]) * velocity[axis]
            })
        });
    }
    for (index, plane) in planes.iter().enumerate() {
        if mask & (1 << index) == 0 {
            continue;
        }
        let normal = plane.inward_normal();
        let projection = (0..3).fold(Dd::ZERO, |sum, axis| {
            sum + Dd::from_f64(normal[axis]) * velocity[axis]
        });
        for axis in 0..3 {
            velocity[axis] = velocity[axis] - Dd::from_f64(2.0 * normal[axis]) * projection;
        }
    }
    velocity
}

fn apply_image(
    mut sample: PreparedMovingCompactSample,
    sign: f64,
    planes: &[IdealAcousticPlane],
    mask: usize,
    inverse_jacobian: Option<&[[f64; 3]; 3]>,
) -> PreparedMovingCompactSample {
    if sign < 0.0 {
        sample.coefficients = sample.coefficients.map(|value| -value);
    }
    if let Some(mut coefficients) = sample.monopole_velocity {
        if sign < 0.0 {
            coefficients = coefficients.map(|row| row.map(|value| -value));
        }
        // A receiver gradient follows the TRANSPOSE of its inverse image map.
        // Force-source pressure columns stay in the original source frame.
        if let Some(jacobian) = inverse_jacobian {
            let original = coefficients;
            coefficients = core::array::from_fn(|row| {
                core::array::from_fn(|column| {
                    (0..3)
                        .fold(Dd::ZERO, |sum, axis| {
                            sum + Dd::from_f64(jacobian[axis][row])
                                * Dd::from_f64(original[axis][column])
                        })
                        .to_f64()
                })
            });
        } else {
            // Orthogonal virtual receivers apply planes in reverse order;
            // retain the original forward reflection arithmetic here.
            for (index, plane) in planes.iter().enumerate() {
                if mask & (1 << index) == 0 {
                    continue;
                }
                let normal = plane.inward_normal();
                for column in 0..2 {
                    let projection = (0..3).fold(Dd::ZERO, |sum, axis| {
                        sum + Dd::from_f64(normal[axis]) * Dd::from_f64(coefficients[axis][column])
                    });
                    for axis in 0..3 {
                        coefficients[axis][column] = (Dd::from_f64(coefficients[axis][column])
                            - Dd::from_f64(2.0 * normal[axis]) * projection)
                            .to_f64();
                    }
                }
            }
        }
        sample.monopole_velocity = coefficients
            .iter()
            .flatten()
            .all(|value| value.is_finite())
            .then_some(coefficients);
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
    maximum_mach: f64,
    budget: ReflectedTrajectoryBudget,
    visited: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, AeroacError> {
    let zero = Interval::point(0.0);
    let mut controls = [[zero; 3]; 6];
    let mut velocity = [[zero; 3]; 5];
    for (control, original) in velocity.iter_mut().zip(segment.velocity) {
        for axis in 0..3 {
            control[axis] = enclose_dd(original[axis])?;
        }
    }
    for axis in 0..3 {
        velocity[0][axis] = velocity[0][axis]
            .hull(Interval::point(segment.first.kinematics.velocity_m_s[axis]));
        velocity[4][axis] = velocity[4][axis]
            .hull(Interval::point(segment.last.kinematics.velocity_m_s[axis]));
    }
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
    let mut stack = [BoundNode { position: controls, velocity, depth: 0 }; MAX_BOUND_DEPTH + 1];
    let mut pending = 1;
    let mut maximum_rate = 0.0_f64;
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
            node.position.iter().all(|control| {
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
                &node.position,
                relative,
                observer.observer.limits().minimum_distance_m,
            )
        });
        if inside && clear {
            for path in &observer.paths {
                maximum_rate = maximum_rate.max(doppler::maximum_rate(
                    &path.observer, origin, &node, maximum_mach,
                ));
            }
            continue;
        }
        if node.depth >= budget.maximum_subdivision_depth {
            return Err(AeroacError::NotConverged {
                what: "reflected source half-space/clearance hull remains inadmissible",
                residual: 1.0,
            });
        }
        let (left, right) = split_controls(node.position);
        let (left_velocity, right_velocity) = split_controls(node.velocity);
        stack[pending] = BoundNode {
            position: right,
            velocity: right_velocity,
            depth: node.depth + 1,
        };
        stack[pending + 1] = BoundNode {
            position: left,
            velocity: left_velocity,
            depth: node.depth + 1,
        };
        pending += 2;
    }
    Ok(maximum_rate)
}
