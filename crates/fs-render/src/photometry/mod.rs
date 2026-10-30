//! Calibrated, headless fixture transport in metres, radians, lumens and lux.
//!
//! The scene borrows the spectral tracer's actual [`Primitive`] geometry and
//! materials. Source intensity is flux-normalized, rather than interpreting
//! display RGB or the legacy relative-XYZ film as a photometric measurement.
//! Direct illuminance evaluates `I(direction) * cos(theta) * T / distance²`.
//! It is the unscattered, unreflected component, not total room illuminance.
//! [`path::PathIlluminanceJob`] adds sampled reflected and volume transport to
//! that direct calculation using the existing spectral path integrator.
//! Calibrated finite mesh sources in [`area`] share those same paths and report
//! their incident illuminance in the sampled continuous component.
//!
//! The tabulated CIE 1931 observer and 80-bin spectral discretization are
//! declared numerical models. Neither manufacturer calibration nor model error
//! has a certified bound here. [`NumericalCertificate::no_claim`] deliberately
//! distinguishes an analytic binary64 calculation from a rigorous enclosure.
//! No Monte Carlo sampling error is invented for a deterministic direct sum.
//! Source order is retained; identical inputs and progressive partitions replay
//! bitwise on the same ISA/toolchain. Cross-ISA equality is not asserted.

use crate::charts::Ray;
use crate::instances::SharedGeometry;
use crate::tracer::photometric::geometry::PreparedGeometry;
use crate::tracer::{Material, Primitive, Shape, TracerError};
use fs_evidence::NumericalCertificate;
use fs_exec::{Cancelled, Cx};
use fs_geom::{Point3, Vec3};
use std::sync::Arc;

pub mod area;
pub mod image;
pub mod medium;
pub mod observer;
pub mod path;
pub mod source;
pub mod transport;

use medium::{HomogeneousMedium, MediumError};
use source::{PhotometricSource, SourceError};
use transport::{PhotometricMedium, photopic_transmittance};

/// Version of the calibrated transport model and output interpretation.
/// Version 5 admits passive tabulated diffuse reflectance. Version 4's
/// shared tabulated haze phase laws, version 3's fixture
/// angular partitions and version 2's published CIE observer remain in use.
pub const PHOTOMETRY_SEMANTICS_VERSION: u32 = 5;

/// Explicit geometric and output budgets for an offline calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransportOptions {
    /// Receiver endpoint exclusion for visibility, in metres. This is a
    /// geometric tolerance, not a claimed bound on illuminance error.
    pub ray_epsilon_m: f64,
    /// Maximum supported connection or camera-ray length in metres. Longer
    /// connections refuse rather than silently losing distant blockers.
    pub max_distance_m: f64,
    /// Maximum retained grid/image/probe values per job.
    pub max_outputs: usize,
}

impl Default for TransportOptions {
    fn default() -> Self {
        Self {
            ray_epsilon_m: 1.0e-6,
            max_distance_m: 10_000.0,
            max_outputs: 16_777_216,
        }
    }
}

/// Refusals preserve invalid input, unsupported physics and cancellation.
#[derive(Debug)]
pub enum PhotometryError {
    /// Cancellation leaves the previously completed output prefix intact.
    Cancelled,
    /// Probe position or normal is non-finite or has no direction.
    InvalidProbe,
    /// Grid axes, dimensions or coordinates are invalid.
    InvalidGrid,
    /// The requested grid index is outside the declared extent.
    GridIndexOutOfRange,
    /// Geometric tolerance, distance or memory budget is invalid.
    InvalidOptions,
    /// Camera ray, distance, image dimensions or aperture are unsupported.
    InvalidCameraRay,
    /// Single-scattering quadrature settings are invalid.
    InvalidScatterSettings,
    /// The declared fixture angular intervals need more panels than the
    /// caller permits. No unresolved interval is silently discarded.
    ScatteringPartitionBudgetExceeded {
        /// Increase this total panel cap to admit the known angular intervals.
        max_panels: u32,
    },
    /// The output count exceeds the explicit job budget.
    OutputBudgetExceeded,
    /// Output storage could not be reserved.
    AllocationRefused,
    /// Point-source approximation is singular at this receiver.
    SourceTooClose {
        /// Source in the caller's retained order.
        source: usize,
    },
    /// A connection is outside the caller's explicit tracing range.
    DistanceLimitExceeded,
    /// An admitted calculation overflowed or produced a non-finite value.
    NonFiniteResult,
    /// A scene object requests physics not represented by this transport rung.
    UnsupportedPrimitive {
        /// Object in the borrowed tracer primitive slice.
        primitive: usize,
        /// Specific unsupported input; never silently treated as opaque.
        reason: &'static str,
    },
    /// Invalid or stale mesh/appearance input.
    InvalidPrimitive {
        /// Object in the borrowed tracer primitive slice.
        primitive: usize,
    },
    /// Calibrated-source validation/evaluation refusal.
    Source(SourceError),
    /// Participating-medium validation/evaluation refusal.
    Medium(MediumError),
    /// The original physical-grid transport refused a field, ray or budget.
    PhysicalMedium(crate::volumes::transport::VolumeTransportError),
    /// The requested image requires medium terms this estimator cannot supply.
    UnsupportedMedium {
        /// Use a component-specific or full spectral path entry point instead.
        reason: &'static str,
    },
    /// The original renderer's geometry backend refused the ray.
    Geometry(TracerError),
    /// The existing physical camera refused ray generation.
    Camera(crate::camera::CameraError),
    /// Bounded parallel execution refused memory, configuration or scheduling.
    Execution(crate::tracer::RenderExecutionError),
}

impl core::fmt::Display for PhotometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "photometric transport refused: {self:?}")
    }
}

impl std::error::Error for PhotometryError {}

impl From<Cancelled> for PhotometryError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

impl From<SourceError> for PhotometryError {
    fn from(value: SourceError) -> Self {
        Self::Source(value)
    }
}

impl From<MediumError> for PhotometryError {
    fn from(value: MediumError) -> Self {
        Self::Medium(value)
    }
}

impl From<crate::volumes::transport::VolumeTransportError> for PhotometryError {
    fn from(value: crate::volumes::transport::VolumeTransportError) -> Self {
        match value {
            crate::volumes::transport::VolumeTransportError::Cancelled => Self::Cancelled,
            other => Self::PhysicalMedium(other),
        }
    }
}

impl From<TracerError> for PhotometryError {
    fn from(value: TracerError) -> Self {
        match value {
            TracerError::Cancelled => Self::Cancelled,
            other => Self::Geometry(other),
        }
    }
}

impl From<crate::camera::CameraError> for PhotometryError {
    fn from(value: crate::camera::CameraError) -> Self {
        match value {
            crate::camera::CameraError::Cancelled => Self::Cancelled,
            other => Self::Camera(other),
        }
    }
}

impl From<crate::tracer::RenderExecutionError> for PhotometryError {
    fn from(value: crate::tracer::RenderExecutionError) -> Self {
        match value {
            crate::tracer::RenderExecutionError::Tracer(error) => error.into(),
            crate::tracer::RenderExecutionError::Executor(fs_exec::RunError::Cancelled {
                ..
            }) => Self::Cancelled,
            other => Self::Execution(other),
        }
    }
}

/// A one-sided infinitesimal receiver, including a free-space calculation point.
/// Normalization never changes its declared position in metres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IlluminanceProbe {
    position_m: Point3,
    normal: Vec3,
}

impl IlluminanceProbe {
    /// Admit a finite position and finite nonzero surface normal.
    pub fn try_new(position_m: Point3, normal: Vec3) -> Result<Self, PhotometryError> {
        if !finite_point(position_m) {
            return Err(PhotometryError::InvalidProbe);
        }
        let normal = unit(normal).ok_or(PhotometryError::InvalidProbe)?;
        Ok(Self { position_m, normal })
    }

    /// World-space receiver position in metres.
    #[must_use]
    pub const fn position_m(self) -> Point3 {
        self.position_m
    }

    /// Unit front-side normal, in the caller's world frame.
    #[must_use]
    pub const fn normal(self) -> Vec3 {
        self.normal
    }
}

/// Point-sampled calculation plane. `step_u` and `step_v` are physical vectors
/// between samples; their cross product defines the receiving side. Values are
/// row-major, `origin + x*step_u + y*step_v`, and are not cell-area averages.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalculationGrid {
    origin: Point3,
    step_u: Vec3,
    step_v: Vec3,
    normal: Vec3,
    width: u32,
    height: u32,
    len: usize,
}

impl CalculationGrid {
    /// Construct a plane without allocating its probes.
    pub fn try_new(
        origin: Point3,
        step_u: Vec3,
        step_v: Vec3,
        width: u32,
        height: u32,
    ) -> Result<Self, PhotometryError> {
        if !finite_point(origin) || width == 0 || height == 0 {
            return Err(PhotometryError::InvalidGrid);
        }
        let u = unit(step_u).ok_or(PhotometryError::InvalidGrid)?;
        let v = unit(step_v).ok_or(PhotometryError::InvalidGrid)?;
        let normal = unit(cross(u, v)).ok_or(PhotometryError::InvalidGrid)?;
        let len = usize::try_from(u64::from(width) * u64::from(height))
            .map_err(|_| PhotometryError::InvalidGrid)?;
        let grid = Self {
            origin,
            step_u,
            step_v,
            normal,
            width,
            height,
            len,
        };
        // An affine plane is finite throughout when all four corners are.
        for index in [0, width as usize - 1, len - width as usize, len - 1] {
            grid.probe(index)?;
        }
        Ok(grid)
    }

    /// Number of requested receiver samples.
    #[must_use]
    pub const fn len(self) -> usize {
        self.len
    }

    /// An admitted grid always contains at least one sample.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        false
    }

    /// Raster dimensions `(width, height)`.
    #[must_use]
    pub const fn dimensions(self) -> (u32, u32) {
        (self.width, self.height)
    }

    /// Generate one calculation point; no geometry conversion is performed.
    pub fn probe(self, index: usize) -> Result<IlluminanceProbe, PhotometryError> {
        if index >= self.len {
            return Err(PhotometryError::GridIndexOutOfRange);
        }
        let x = (index % self.width as usize) as f64;
        let y = (index / self.width as usize) as f64;
        IlluminanceProbe::try_new(
            self.origin
                .offset(self.step_u.scale(x))
                .offset(self.step_v.scale(y)),
            self.normal,
        )
    }
}

/// The calibrated direct component of illuminance at one receiver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IlluminanceEstimate {
    /// Unscattered direct illuminance in lux (lm/m²).
    pub lux: f64,
    /// No rigorous floating-point enclosure is currently supplied. The
    /// analytic formula has no sampling variance, which does not make it exact.
    pub numerical: NumericalCertificate,
    /// Bound on calibration, geometry, observer and omitted reflected/scattered
    /// contributions. `None` means unknown, not zero or a confidence interval.
    pub total_error_bound_lux: Option<f64>,
}

/// Immutable borrowed scene for a snapshot of fixture transforms/state.
///
/// Geometry and opaque materials are exactly those used by the existing
/// spectral renderer, including material-state bindings. Static snapshots are
/// intentional: advance fixture transforms by constructing a new scene/job.
/// Point emitters can represent every emitter in a multi-emitter fixture.
/// Material optics are used by the image path; direct lux depends on visibility.
/// A prepared scene index reuses the native mesh/instance intersections and is
/// shared by calibrated path jobs. Its broad phase preserves original leaf
/// order; it changes neither the physical estimator nor its uncertainty.
pub struct PhotometricScene<'a> {
    primitives: &'a [Primitive],
    geometry: Arc<PreparedGeometry>,
    sources: &'a [PhotometricSource],
    medium: PhotometricMedium<'a>,
    options: TransportOptions,
}

impl<'a> PhotometricScene<'a> {
    /// Admit native triangle meshes and rigid mesh instances. Chart proximity
    /// hits do not establish whether a surface lies before a finite connection
    /// endpoint, so chart visibility explicitly refuses. Transparent objects,
    /// animated geometry and legacy uncalibrated emission also refuse instead
    /// of disappearing from transport or changing units. Meshes must retain
    /// their original BVHs; do not mutate `TriMesh`'s public vertices after
    /// construction.
    pub fn try_new(
        primitives: &'a [Primitive],
        sources: &'a [PhotometricSource],
        medium: HomogeneousMedium,
        options: TransportOptions,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::try_with_medium(
            primitives,
            sources,
            PhotometricMedium::Homogeneous(medium),
            options,
            cx,
        )
    }

    /// Admit the same static scene with an original physical medium snapshot.
    /// Direct illuminance includes only fixture light which has neither
    /// scattered nor reflected. Temperature-dependent medium emission is not
    /// part of that direct component. Image methods explicitly distinguish the
    /// direct surface component from unsupported grid haze/thermal emission.
    pub fn try_with_medium(
        primitives: &'a [Primitive],
        sources: &'a [PhotometricSource],
        medium: PhotometricMedium<'a>,
        options: TransportOptions,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::admit(primitives, sources, medium, options, false, cx)
    }

    /// The spectral path adapter supplies native tangent frames and dielectric
    /// boundary state. The direct-only image adapter does not, so its public
    /// admission must retain the narrower material contract.
    pub(crate) fn try_with_path_medium(
        primitives: &'a [Primitive],
        sources: &'a [PhotometricSource],
        medium: PhotometricMedium<'a>,
        options: TransportOptions,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::admit(primitives, sources, medium, options, true, cx)
    }

    fn admit(
        primitives: &'a [Primitive],
        sources: &'a [PhotometricSource],
        medium: PhotometricMedium<'a>,
        options: TransportOptions,
        native_paths: bool,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        if !options.ray_epsilon_m.is_finite()
            || options.ray_epsilon_m <= 0.0
            || !options.max_distance_m.is_finite()
            || options.max_distance_m <= options.ray_epsilon_m
            || options.max_outputs == 0
        {
            return Err(PhotometryError::InvalidOptions);
        }
        let has_point_emission = sources.iter().any(|source| source.flux_lm() > 0.0);
        for (index, primitive) in primitives.iter().enumerate() {
            cx.checkpoint()?;
            let unsupported = match &primitive.shape {
                Shape::AnimatedInstance(_) => Some("supply a static fixture/geometry snapshot"),
                Shape::Chart(_) => Some("chart segment endpoint visibility is not yet supported"),
                Shape::Instance(instance)
                    if matches!(instance.geometry(), SharedGeometry::Chart(_)) =>
                {
                    Some("chart segment endpoint visibility is not yet supported")
                }
                _ if primitive.emission.is_some() => {
                    Some("legacy emission has no luminous calibration")
                }
                _ => match primitive.material {
                    Material::Dielectric { .. } if !native_paths => {
                        Some("refractive connections are not yet supported")
                    }
                    Material::Dielectric { .. } if has_point_emission => {
                        // A continuous refractive walk can hit finite emitters,
                        // but cannot find a delta-position point fixture. The
                        // missing specular point connections must not become a
                        // silent zero inside a result labelled total lux.
                        Some(
                            "active point fixtures with glass require specular point-source connections",
                        )
                    }
                    Material::BrushedConductor { .. } if !native_paths => {
                        Some("anisotropic tangent-frame shading is not yet supported")
                    }
                    Material::GrayDiffuseThermal { .. }
                    | Material::SpectralDiffuseThermal { .. } => {
                        Some("thermal sources require a separate calibrated emission integral")
                    }
                    _ => None,
                },
            };
            if let Some(reason) = unsupported {
                return Err(PhotometryError::UnsupportedPrimitive {
                    primitive: index,
                    reason,
                });
            }
            match primitive.material {
                Material::Lambertian { reflectance } | Material::Ggx { reflectance, .. } => {
                    if reflectance.c.iter().any(|c| !c.is_finite()) {
                        return Err(PhotometryError::InvalidPrimitive { primitive: index });
                    }
                }
                _ => {}
            }
            if let Material::Ggx { alpha, .. } = primitive.material {
                if !alpha.is_finite() || alpha <= 0.0 {
                    return Err(PhotometryError::InvalidPrimitive { primitive: index });
                }
            }
            let mesh = match &primitive.shape {
                Shape::Mesh(mesh) => Some(mesh),
                Shape::Instance(instance) => match instance.geometry() {
                    SharedGeometry::Mesh(mesh) => Some(mesh.as_ref()),
                    SharedGeometry::Chart(_) => None,
                },
                _ => None,
            };
            if let Some(mesh) = mesh {
                for vertex in &mesh.vertices {
                    cx.checkpoint()?;
                    if vertex.iter().any(|v| !v.is_finite()) {
                        return Err(PhotometryError::InvalidPrimitive { primitive: index });
                    }
                }
                for triangle in &mesh.triangles {
                    cx.checkpoint()?;
                    if triangle.iter().any(|&v| v as usize >= mesh.vertices.len()) {
                        return Err(PhotometryError::InvalidPrimitive { primitive: index });
                    }
                }
            }
        }
        Ok(Self {
            primitives,
            geometry: Arc::new(PreparedGeometry::try_new(primitives, cx)?),
            sources,
            medium,
            options,
        })
    }

    pub(crate) fn prepared_geometry(&self) -> Arc<PreparedGeometry> {
        Arc::clone(&self.geometry)
    }

    #[cfg(test)]
    pub(crate) fn disable_acceleration_for_test(&mut self) {
        self.geometry = Arc::new(PreparedGeometry::linear_for_test());
    }

    /// Direct lux from every source, in caller order with compensated summation.
    /// Receiver endpoint exclusion affects visibility only, never the inverse
    /// square or extinction distance. Every source and BVH/chart polls `Cx`.
    pub fn illuminance(
        &self,
        cx: &Cx<'_>,
        probe: &IlluminanceProbe,
    ) -> Result<IlluminanceEstimate, PhotometryError> {
        cx.checkpoint()?;
        let mut total = CompensatedSum::default();
        for (index, source) in self.sources.iter().enumerate() {
            cx.checkpoint()?;
            if source.flux_lm() == 0.0 {
                continue;
            }
            let connection = self.connection(cx, source, index, probe.position_m)?;
            let cosine = probe.normal.dot(connection.towards_source).max(0.0);
            if cosine == 0.0 || connection.illuminance_normal_lux == 0.0 {
                continue;
            }
            if self.visible(cx, source.position_m(), &connection)? {
                total.add(connection.illuminance_normal_lux * cosine);
            }
        }
        let lux = total.finish()?;
        Ok(IlluminanceEstimate {
            lux,
            numerical: NumericalCertificate::no_claim(),
            total_error_bound_lux: None,
        })
    }

    /// Evaluate a full row-major grid with private staging. On refusal there is
    /// no completed-grid result; use [`IlluminanceJob::for_grid`] for progress.
    pub fn illuminance_grid(
        &self,
        grid: &CalculationGrid,
        cx: &Cx<'_>,
    ) -> Result<Vec<IlluminanceEstimate>, PhotometryError> {
        let mut job = IlluminanceJob::for_grid(self, *grid)?;
        job.advance(cx, grid.len())?;
        Ok(job.values)
    }

    fn connection(
        &self,
        cx: &Cx<'_>,
        source: &PhotometricSource,
        index: usize,
        receiver: Point3,
    ) -> Result<LightConnection, PhotometryError> {
        let mut connection = self.unattenuated_connection(source, index, receiver)?;
        if connection.illuminance_normal_lux == 0.0 {
            return Ok(connection);
        }
        let transmission = match self.medium.homogeneous() {
            Some(medium) => medium.transmittance(connection.distance_m)?,
            None => {
                let values = self.medium.transmittance(
                    cx,
                    source.position_m(),
                    connection.towards_source.scale(-1.0),
                    connection.distance_m,
                )?;
                photopic_transmittance(source.spectrum(), &values)?
            }
        };
        connection.illuminance_normal_lux *= transmission;
        if !connection.illuminance_normal_lux.is_finite() {
            return Err(PhotometryError::NonFiniteResult);
        }
        Ok(connection)
    }

    fn unattenuated_connection(
        &self,
        source: &PhotometricSource,
        index: usize,
        receiver: Point3,
    ) -> Result<LightConnection, PhotometryError> {
        let delta = source.position_m().delta_from(receiver);
        let towards_source =
            unit(delta).ok_or(PhotometryError::SourceTooClose { source: index })?;
        let distance_m = scaled_norm(delta).ok_or(PhotometryError::NonFiniteResult)?;
        if distance_m <= self.options.ray_epsilon_m {
            return Err(PhotometryError::SourceTooClose { source: index });
        }
        if distance_m > self.options.max_distance_m {
            return Err(PhotometryError::DistanceLimitExceeded);
        }
        let intensity_cd = source.intensity_cd(towards_source.scale(-1.0))?;
        // Divide sequentially to avoid overflow of distance² for distant rigs.
        let illuminance_normal_lux = intensity_cd / distance_m / distance_m;
        if !illuminance_normal_lux.is_finite() || illuminance_normal_lux < 0.0 {
            return Err(PhotometryError::NonFiniteResult);
        }
        Ok(LightConnection {
            towards_source,
            distance_m,
            illuminance_normal_lux,
        })
    }

    fn visible(
        &self,
        cx: &Cx<'_>,
        emitter: Point3,
        connection: &LightConnection,
    ) -> Result<bool, PhotometryError> {
        let ray = Ray {
            origin: emitter,
            dir: connection.towards_source.scale(-1.0),
        };
        let end = connection.distance_m - self.options.ray_epsilon_m;
        Ok(self
            .geometry
            .intersect(self.primitives, cx, &ray, end)?
            .is_none())
    }
}

struct LightConnection {
    towards_source: Vec3,
    distance_m: f64,
    illuminance_normal_lux: f64,
}

enum Receivers<'a> {
    Probes(&'a [IlluminanceProbe]),
    Grid(CalculationGrid),
}

/// Resumable direct calculation bound by immutable borrows to one scene and one
/// receiver set. Each completed receiver commits atomically. Cancellation or a
/// backend refusal leaves the completed prefix available and the current
/// receiver uncommitted; resume using a fresh, uncancelled context.
pub struct IlluminanceJob<'a, 'scene> {
    scene: &'a PhotometricScene<'scene>,
    receivers: Receivers<'a>,
    count: usize,
    values: Vec<IlluminanceEstimate>,
}

impl<'a, 'scene> IlluminanceJob<'a, 'scene> {
    /// Begin a calculation on caller-supplied surface points.
    pub fn new(
        scene: &'a PhotometricScene<'scene>,
        probes: &'a [IlluminanceProbe],
    ) -> Result<Self, PhotometryError> {
        Self::begin(scene, Receivers::Probes(probes), probes.len())
    }

    /// Begin a grid calculation without allocating a separate probe array.
    pub fn for_grid(
        scene: &'a PhotometricScene<'scene>,
        grid: CalculationGrid,
    ) -> Result<Self, PhotometryError> {
        Self::begin(scene, Receivers::Grid(grid), grid.len())
    }

    fn begin(
        scene: &'a PhotometricScene<'scene>,
        receivers: Receivers<'a>,
        count: usize,
    ) -> Result<Self, PhotometryError> {
        if count > scene.options.max_outputs {
            return Err(PhotometryError::OutputBudgetExceeded);
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(count)
            .map_err(|_| PhotometryError::AllocationRefused)?;
        Ok(Self {
            scene,
            receivers,
            count,
            values,
        })
    }

    /// Compute at most `max_points` complete receivers, returning the number
    /// appended by this call. Partition size cannot alter arithmetic order.
    pub fn advance(&mut self, cx: &Cx<'_>, max_points: usize) -> Result<usize, PhotometryError> {
        cx.checkpoint()?;
        let start = self.values.len();
        let end = start.saturating_add(max_points).min(self.count);
        for index in start..end {
            let probe = match self.receivers {
                Receivers::Probes(probes) => probes[index],
                Receivers::Grid(grid) => grid.probe(index)?,
            };
            let value = self.scene.illuminance(cx, &probe)?;
            cx.checkpoint()?;
            self.values.push(value);
        }
        Ok(self.values.len() - start)
    }

    /// Only the successfully completed prefix; unfinished receivers are absent.
    #[must_use]
    pub fn results(&self) -> &[IlluminanceEstimate] {
        &self.values
    }

    /// Whether every declared receiver has completed.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.values.len() == self.count
    }
}

#[derive(Default)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }

    fn finish(self) -> Result<f64, PhotometryError> {
        let value = self.sum + self.correction;
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(PhotometryError::NonFiniteResult)
        }
    }
}

fn finite_point(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

fn unit(v: Vec3) -> Option<Vec3> {
    if !v.x.is_finite() || !v.y.is_finite() || !v.z.is_finite() {
        return None;
    }
    let scale = v.x.abs().max(v.y.abs()).max(v.z.abs());
    if scale == 0.0 {
        return None;
    }
    let scaled = Vec3::new(v.x / scale, v.y / scale, v.z / scale);
    Some(scaled.scale(1.0 / scaled.norm()))
}

fn scaled_norm(v: Vec3) -> Option<f64> {
    let scale = v.x.abs().max(v.y.abs()).max(v.z.abs());
    if !scale.is_finite() || scale == 0.0 {
        return None;
    }
    let norm = Vec3::new(v.x / scale, v.y / scale, v.z / scale).norm() * scale;
    (norm.is_finite() && norm > 0.0).then_some(norm)
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}
