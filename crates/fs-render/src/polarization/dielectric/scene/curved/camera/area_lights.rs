//! Finite-area source irradiance through original native refractive images.
//!
//! A prepared source borrows an original surface emitter and native mesh faces.
//! It samples area without tessellating charts or constructing surrogate lamps.
//! For one area sample, Le*cos(source)/pdf_area is a differential radiant
//! intensity. The native polarized bundle applies Fresnel transfer, absorption,
//! focusing and receiver obliquity exactly once. No inverse-square or additional
//! refractive-index factor is introduced.
//!
//! The standalone irradiance estimates are composed with camera radiance by
//! [`camera`], using matching light-sampling and emitter-hit MIS weights.
//! Finite authored image seeds do not establish global image completeness.

use super::{CameraRadianceError, CurvedCameraScene, EmissionSides, SurfaceEmitter};
use crate::charts::Ray;
#[cfg(test)]
use crate::charts::TriMesh;
use crate::charts::area_source::MeshAreaError;
use crate::gradient_index::bundle::ProjectionOptions;
use crate::gradient_index::event::BoundarySide;
use crate::instances::{GeometryInstance, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::Interaction;
use crate::polarization::dielectric::scene::curved::CurvedPolarizedScene;
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSceneError;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleConnection, BundleConnectionError, BundleConnectionLaunch, BundleConnectionOptions,
    BundleConnectionTarget, BundleSurfaceSource,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageCatalog, BundleImageError, BundleImageOptions, BundleImageSeed,
    BundleImageSeedOutcome, BundleImageWorkspace,
};
use crate::polarization::dielectric::transport::bundle::PolarizedBundleError;
#[cfg(feature = "volumes")]
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighScene;
#[cfg(feature = "volumes")]
use crate::polarization::dielectric::scene::curved::participating::transmittance::CurvedTransmittanceOptions;
#[cfg(feature = "volumes")]
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::{
    CurvedTransmittanceImageError, CurvedTransmittanceImageWorkspace,
};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// Reciprocal optical-image densities and complementary camera MIS weights.
pub mod mis;

/// Finite-area next-event sampling and complementary native emitter hits.
pub mod camera;

/// Original molecular camera paths with finite-area surface and volume lighting.
#[cfg(feature = "volumes")]
pub mod participating;

#[derive(Clone, Debug, PartialEq)]
pub enum AreaEmitterError {
    InvalidEmitter { object_id: u64 },
    UnsupportedGeometry,
    TriangleBudget { required: usize, limit: usize },
    InvalidTriangle { triangle: usize },
    Mesh(MeshAreaError),
    AllocationFailure,
    InvalidSample,
    SampleOwnerMismatch,
    UnsupportedEmissionSide,
    SourceContactMismatch,
    ScatteringOwnerMismatch,
    EmptyImages,
    UnresolvedSeed { seed: usize, error: BundleConnectionError },
    Camera(CameraRadianceError),
    Source(CurvedBundleSceneError),
    Images(BundleImageError),
    Refractive(PolarizedBundleError),
    #[cfg(feature = "volumes")]
    Survival(CurvedTransmittanceImageError),
    NumericalRange,
    Cancelled,
}

impl From<fs_exec::Cancelled> for AreaEmitterError {
    fn from(_: fs_exec::Cancelled) -> Self {
        Self::Cancelled
    }
}
impl From<MeshAreaError> for AreaEmitterError {
    fn from(error: MeshAreaError) -> Self { Self::Mesh(error) }
}
impl core::fmt::Display for AreaEmitterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native finite-area emitter: {self:?}")
    }
}
impl core::error::Error for AreaEmitterError {}

/// Bounds the original mesh, preparation and each native hit lookup. Subsets
/// still inspect excluded original faces for overlapping source support.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AreaEmitterLimits {
    /// Whole original mesh bound, also bounding selected CDF construction.
    pub max_triangles: usize,
    /// Native BVH visits, candidate pairs and source faces during admission.
    pub max_preparation_work: u64,
    /// Native BVH visits and candidate faces in one contact membership query.
    pub max_lookup_work: u64,
}

impl Default for AreaEmitterLimits {
    fn default() -> Self {
        Self { max_triangles: 65_536, max_preparation_work: 4_194_304, max_lookup_work: 65_536 }
    }
}

#[derive(Debug)]
struct TriangleArea {
    original_triangle: usize,
    area_m2: f64,
    cumulative_m2: f64,
}

/// Original source geometry and radiance, with bounded sampling metadata.
/// All original mesh faces or an explicit sorted subset may emit. Selected
/// faces have disjoint interiors from every coplanar original face; no mesh
/// union or multiplicity correction is inferred. Faces need not be connected
/// or coplanar. Their original windings and per-side radiance are retained.
/// The mesh is borrowed, never copied or converted. Image seeds remain world
/// directions and must support the selected face's emitting hemisphere.
pub struct CurvedAreaEmitter<'s, 'm, 'e> {
    scene: &'s CurvedPolarizedScene<'m>,
    emitter: &'e SurfaceEmitter,
    instance: &'m GeometryInstance,
    triangles: Vec<TriangleArea>,
    total_area_m2: f64,
    limits: AreaEmitterLimits,
}

/// A sealed area sample from the original emitter and geometry instance.
/// The density includes triangle selection and is measured in inverse m².
/// It does not include a side, direction, wavelength or optical-image PDF.
#[derive(Clone, Copy, Debug)]
pub struct AreaEmitterSample<'a> {
    emitter: &'a SurfaceEmitter,
    instance: &'a GeometryInstance,
    // This exact preparation owns the conditional subset PDF.
    support: &'a [TriangleArea],
    triangle: usize,
    position_m: Point3,
    normal: Vec3,
    pdf_area_per_m2: f64,
}

impl AreaEmitterSample<'_> {
    pub const fn object_id(&self) -> u64 { self.emitter.object_id() }
    pub const fn triangle_index(&self) -> usize { self.triangle }
    pub const fn position_m(&self) -> Point3 { self.position_m }
    pub const fn normal(&self) -> Vec3 { self.normal }
    pub const fn pdf_area_per_m2(&self) -> f64 { self.pdf_area_per_m2 }
}

/// One Monte Carlo estimate of the area integral, summing distinct images at
/// the sampled source point. Average independent returned estimates; do not
/// divide by the number of images. No partial estimate is returned on failure.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AreaIrradianceSample {
    irradiance_w_m2_nm: f64,
    images: u32,
    traces: u32,
    fields: u64,
    geometry: u64,
    source_draws: u32,
    source_candidates: u32,
}

impl AreaIrradianceSample {
    pub const fn spectral_irradiance_w_m2_nm(&self) -> f64 { self.irradiance_w_m2_nm }
    pub const fn image_count(&self) -> u32 { self.images }
    pub const fn trace_evaluations(&self) -> u32 { self.traces }
    pub const fn field_evaluations(&self) -> u64 { self.fields }
    pub const fn query_work(&self) -> u64 { self.geometry }
    pub const fn source_free_flight_draws(&self) -> u32 { self.source_draws }
    pub const fn source_candidates(&self) -> u32 { self.source_candidates }
}

impl<'s, 'm, 'e> CurvedAreaEmitter<'s, 'm, 'e> {
    /// Bind the existing calibrated spectrum to one original native source
    /// object. Detector receivers are supported; the source itself cannot be
    /// a detector. No camera scene or replacement photometric source is built.
    pub fn try_new(
        cx: &Cx<'_>,
        scene: &'s CurvedPolarizedScene<'m>,
        emitter: &'e SurfaceEmitter,
        max_triangles: usize,
    ) -> Result<Self, AreaEmitterError> {
        Self::try_new_with_limits(cx, scene, emitter,
            AreaEmitterLimits { max_triangles, ..AreaEmitterLimits::default() })
    }

    /// Prepare all original faces with explicit geometry and work limits.
    pub fn try_new_with_limits(
        cx: &Cx<'_>, scene: &'s CurvedPolarizedScene<'m>, emitter: &'e SurfaceEmitter,
        limits: AreaEmitterLimits,
    ) -> Result<Self, AreaEmitterError> {
        Self::prepare(cx, scene, emitter, None, limits)
    }

    /// Prepare a strictly increasing, nonempty subset of original triangle
    /// indices. Area density is conditional on this support; another disjoint
    /// source preparation can use its own world-space image seeds. Every
    /// selected face is checked against excluded faces for coplanar overlap.
    pub fn try_new_faces(
        cx: &Cx<'_>, scene: &'s CurvedPolarizedScene<'m>, emitter: &'e SurfaceEmitter,
        faces: &[usize], limits: AreaEmitterLimits,
    ) -> Result<Self, AreaEmitterError> {
        Self::prepare(cx, scene, emitter, Some(faces), limits)
    }

    fn prepare(
        cx: &Cx<'_>, scene: &'s CurvedPolarizedScene<'m>, emitter: &'e SurfaceEmitter,
        faces: Option<&[usize]>, limits: AreaEmitterLimits,
    ) -> Result<Self, AreaEmitterError> {
        cx.checkpoint()?;
        let object_id = emitter.object_id();
        let invalid = || AreaEmitterError::InvalidEmitter { object_id };
        let binding = scene.bindings
            .binary_search_by_key(&object_id, |binding| binding.object_id)
            .map_err(|_| invalid())?;
        if matches!(scene.bindings[binding].interaction, Interaction::Detector { .. }) {
            return Err(invalid());
        }
        let geometry = scene.geometry;
        let index = geometry.instances()
            .binary_search_by_key(&object_id, GeometryInstance::object_id)
            .map_err(|_| invalid())?;
        let instance = &geometry.instances()[index];
        let SharedGeometry::Mesh(mesh) = instance.geometry() else {
            return Err(AreaEmitterError::UnsupportedGeometry);
        };
        if mesh.triangles.is_empty() || limits.max_triangles == 0
            || mesh.triangles.len() > limits.max_triangles {
            return Err(AreaEmitterError::TriangleBudget {
                required: mesh.triangles.len(), limit: limits.max_triangles,
            });
        }
        let count = faces.map_or(mesh.triangles.len(), <[usize]>::len);
        if count == 0 || count > mesh.triangles.len() {
            return Err(MeshAreaError::InvalidSupport.into());
        }
        let mut selected = Vec::new();
        selected.try_reserve_exact(count).map_err(|_| AreaEmitterError::AllocationFailure)?;
        if let Some(faces) = faces { selected.extend_from_slice(faces); }
        else { selected.extend(0..mesh.triangles.len()); }
        mesh.validate_area_source_faces(cx, &selected, limits.max_preparation_work)?;
        let mut triangles = Vec::new();
        triangles.try_reserve_exact(count)
            .map_err(|_| AreaEmitterError::AllocationFailure)?;
        let mut sum = Sum::default();
        let mut previous = 0.0;
        for triangle in selected {
            cx.checkpoint()?;
            let indices = mesh.triangles[triangle];
            let points = indices.map(|index| mesh.vertices.get(index as usize).copied());
            let [Some(a), Some(b), Some(c)] = points else {
                return Err(AreaEmitterError::InvalidTriangle { triangle });
            };
            if a.into_iter().chain(b).chain(c).any(|x| !x.is_finite()) {
                return Err(AreaEmitterError::InvalidTriangle { triangle });
            }
            let cross = cross(delta(b, a), delta(c, a));
            let area = 0.5 * norm(cross);
            if !area.is_finite() || area <= 0.0 {
                return Err(AreaEmitterError::InvalidTriangle { triangle });
            }
            sum.add(area);
            let cumulative = sum.value();
            if !cumulative.is_finite() || cumulative <= previous {
                return Err(AreaEmitterError::NumericalRange);
            }
            triangles.push(TriangleArea { original_triangle: triangle,
                area_m2: area, cumulative_m2: cumulative });
            previous = cumulative;
        }
        cx.checkpoint()?;
        Ok(Self { scene, emitter, instance, triangles, total_area_m2: previous, limits })
    }

    pub const fn source_emitter(&self) -> &'e SurfaceEmitter { self.emitter }
    pub const fn source_instance(&self) -> &'m GeometryInstance { self.instance }
    pub const fn total_area_m2(&self) -> f64 { self.total_area_m2 }

    /// Disjoint conditional source supports can share one original emitter.
    /// Preparation already excludes positive-area overlap between distinct
    /// original faces, so sorted native face identity determines overlap here.
    pub(super) fn support_overlaps(&self, other: &Self) -> bool {
        if !core::ptr::eq(self.instance, other.instance) { return false; }
        let (mut a, mut b) = (0, 0);
        while a < self.triangles.len() && b < other.triangles.len() {
            let left = self.triangles[a].original_triangle;
            let right = other.triangles[b].original_triangle;
            if left == right { return true; }
            if left < right { a += 1; } else { b += 1; }
        }
        false
    }

    /// Sample the original mesh with two variates in [0,1). An exact triangle
    /// edge/vertex is unresolved for the source departure and refuses, without
    /// jitter, snapping, or success-conditioned resampling.
    pub fn sample_area(&self, samples: [f64; 2]) -> Result<AreaEmitterSample<'_>, AreaEmitterError> {
        if samples.into_iter().any(|u| !u.is_finite() || !(0.0..1.0).contains(&u)) {
            return Err(AreaEmitterError::InvalidSample);
        }
        let target = samples[0] * self.total_area_m2;
        let slot = self.triangles.partition_point(|entry| entry.cumulative_m2 <= target);
        let entry = self.triangles.get(slot).ok_or(AreaEmitterError::InvalidSample)?;
        let triangle = entry.original_triangle;
        let lower = if slot == 0 { 0.0 } else { self.triangles[slot - 1].cumulative_m2 };
        let interval = entry.cumulative_m2 - lower;
        let u = (target - lower) / interval;
        let root = det::sqrt(u);
        let barycentric = [1.0 - root, root * (1.0 - samples[1]), root * samples[1]];
        if barycentric.into_iter().any(|w| !w.is_finite() || w <= 0.0 || w >= 1.0) {
            return Err(AreaEmitterError::InvalidSample);
        }
        let SharedGeometry::Mesh(mesh) = self.instance.geometry() else {
            return Err(AreaEmitterError::UnsupportedGeometry);
        };
        let [a, b, c] = mesh.triangles[triangle].map(|index| mesh.vertices[index as usize]);
        let first = delta(b, a);
        let second = delta(c, a);
        let local = Point3::new(
            a[0] + first.x * barycentric[1] + second.x * barycentric[2],
            a[1] + first.y * barycentric[1] + second.y * barycentric[2],
            a[2] + first.z * barycentric[1] + second.z * barycentric[2],
        );
        let position_m = self.instance.transform().transform_point(local);
        let raw_normal = self.instance.transform().transform_vector(cross(first, second));
        let length = norm(raw_normal);
        let normal = Vec3::new(raw_normal.x / length, raw_normal.y / length, raw_normal.z / length);
        // Use the actual represented CDF interval so sampling and its density
        // agree even when accumulated-area rounding changes a face's interval.
        let pdf_area_per_m2 = (interval / self.total_area_m2) / entry.area_m2;
        if ![position_m.x, position_m.y, position_m.z, normal.x, normal.y, normal.z,
            pdf_area_per_m2].into_iter().all(f64::is_finite) || pdf_area_per_m2 <= 0.0 {
            return Err(AreaEmitterError::NumericalRange);
        }
        Ok(AreaEmitterSample { emitter: self.emitter, instance: self.instance,
            support: self.triangles.as_slice(), triangle,
            position_m, normal, pdf_area_per_m2 })
    }

    /// Incident spectral irradiance at an original detector or scattering
    /// surface, or transverse irradiance at a native volume target. Continuous
    /// surface targets use their actual native tangent, including curved
    /// receiving geometry. Each supplied seed must
    /// depart the explicitly selected source side. Both-sided radiance is per
    /// side: selecting a side introduces no factor of two or one-half.
    #[allow(clippy::too_many_arguments)]
    pub fn incident_irradiance(
        &self, cx: &Cx<'_>, sample: AreaEmitterSample<'_>, side: BoundarySide,
        initial_medium: usize, transverse_axis: Vec3, wavelength_nm: f64,
        target: BundleConnectionTarget, connection: BundleConnectionOptions,
        seeds: &[BundleImageSeed<'_>], images: BundleImageOptions,
        workspace: &mut BundleImageWorkspace,
    ) -> Result<AreaIrradianceSample, AreaEmitterError> {
        cx.checkpoint()?;
        let (launch, source, radiance) = self.launch(sample, side, initial_medium,
            transverse_axis, wavelength_nm, seeds)?;
        let catalog = self.scene.discover_bundle_images_from_surface_with_workspace(
            cx, launch, source, target, connection, seeds, images, workspace,
        ).map_err(AreaEmitterError::Images)?;
        complete(catalog)?;
        self.sum_images(sample, side, radiance, target, connection.projection,
            catalog.images().iter().map(|image| (image.connection(), 0.0, false)),
            AreaIrradianceSample { traces: catalog.trace_evaluations(),
                fields: catalog.field_evaluations(), geometry: catalog.query_work(),
                ..AreaIrradianceSample::default() })
    }

    /// The molecular owner must borrow this exact original native scene.
    /// Discovery precedes independent source tracking; duplicate seeds consume
    /// no duplicate survival estimate. Image, trace and tracking bounds remain
    /// cumulative, and every required native seed must resolve.
    #[cfg(feature = "volumes")]
    #[allow(clippy::too_many_arguments)]
    pub fn incident_irradiance_with_transmittance_samples(
        &self, cx: &Cx<'_>, sample: AreaEmitterSample<'_>, side: BoundarySide,
        initial_medium: usize, transverse_axis: Vec3, wavelength_nm: f64,
        target: BundleConnectionTarget, connection: BundleConnectionOptions,
        seeds: &[BundleImageSeed<'_>], images: BundleImageOptions,
        scattering: &CurvedRayleighScene<'_, '_, '_>,
        tracking: CurvedTransmittanceOptions, workspace: &mut CurvedTransmittanceImageWorkspace,
        samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<AreaIrradianceSample, AreaEmitterError> {
        cx.checkpoint()?;
        if !core::ptr::eq(self.scene, scattering.surface_scene()) {
            return Err(AreaEmitterError::ScatteringOwnerMismatch);
        }
        let (launch, source, radiance) = self.launch(sample, side, initial_medium,
            transverse_axis, wavelength_nm, seeds)?;
        let catalog = scattering.discover_bundle_images_from_surface_with_transmittance_samples(
            cx, launch, source, target, connection, seeds, images, tracking, workspace, samples,
        ).map_err(AreaEmitterError::Survival)?;
        complete(catalog.geometric_catalog())?;
        self.sum_images(sample, side, radiance, target, connection.projection,
            catalog.images().iter().map(|image| (image.connection(),
                image.log_scattering_survival(), image.ratio_tracking_zero())),
            AreaIrradianceSample { traces: catalog.trace_evaluations(),
                fields: catalog.field_evaluations(), geometry: catalog.query_work(),
                source_draws: catalog.free_flight_draws(),
                source_candidates: catalog.tracking_candidates(),
                ..AreaIrradianceSample::default() })
    }

    fn launch(
        &self, sample: AreaEmitterSample<'_>, side: BoundarySide, initial_medium: usize,
        transverse_axis: Vec3, wavelength_nm: f64, seeds: &[BundleImageSeed<'_>],
    ) -> Result<(BundleConnectionLaunch, BundleSurfaceSource, f64), AreaEmitterError> {
        if !core::ptr::eq(sample.emitter, self.emitter)
            || !core::ptr::eq(sample.instance, self.instance)
            || !core::ptr::eq(sample.support, self.triangles.as_slice()) {
            return Err(AreaEmitterError::SampleOwnerMismatch);
        }
        if side == BoundarySide::OnSurface
            || (self.emitter.sides() == EmissionSides::Front && side != BoundarySide::Positive) {
            return Err(AreaEmitterError::UnsupportedEmissionSide);
        }
        let seed = seeds.first().ok_or(AreaEmitterError::Images(BundleImageError::EmptySeeds))?;
        let source = BundleSurfaceSource::try_new(self.emitter.object_id(), side)
            .map_err(AreaEmitterError::Source)?;
        let radiance = self.emitter.spectrum().evaluate(wavelength_nm).map_err(AreaEmitterError::Camera)?;
        Ok((BundleConnectionLaunch { ray: Ray { origin: sample.position_m,
            dir: seed.initial_direction }, transverse_axis, wavelength_nm,
            normalized_source_stokes: Stokes::unpolarized(1.0), initial_medium }, source, radiance))
    }

    #[allow(clippy::too_many_arguments)]
    fn sum_images<'p>(
        &self, sample: AreaEmitterSample<'_>, side: BoundarySide, radiance: f64,
        target: BundleConnectionTarget, projection: ProjectionOptions,
        images: impl Iterator<Item = (&'p BundleConnection, f64, bool)>,
        mut result: AreaIrradianceSample,
    ) -> Result<AreaIrradianceSample, AreaEmitterError> {
        let mut sum = Sum::default();
        for (connection, log_survival, ratio_zero) in images {
            let path = connection.path();
            let contact = path.source_contact().ok_or(AreaEmitterError::SourceContactMismatch)?;
            let bundle = path.bundle();
            let initial = bundle.propagation().initial_ray();
            if contact.object_id != self.emitter.object_id()
                || contact.geometry_identity != self.instance.geometry_identity()
                || contact.frame_identity != self.instance.frame_identity()
                || initial.origin != sample.position_m
                || contact.witness.normal.dot(sample.normal) < 1.0 - 1e-8 {
                return Err(AreaEmitterError::SourceContactMismatch);
            }
            let signed_cosine = contact.witness.normal.dot(initial.dir);
            if !signed_cosine.is_finite() || match side {
                BoundarySide::Positive => signed_cosine <= 0.0,
                BoundarySide::Negative => signed_cosine >= 0.0,
                BoundarySide::OnSurface => true,
            } {
                return Err(AreaEmitterError::SourceContactMismatch);
            }
            let cosine = signed_cosine.abs().min(1.0);
            // A volume target's plane only localizes the endpoint. Its physical
            // receiving area is perpendicular to the actual arriving ray.
            let normal = if target.is_volume_receiver() {
                bundle.propagation().propagation().ray.dir.scale(-1.0)
            } else if target.is_surface_receiver() {
                let contact = path.last_contact().ok_or(AreaEmitterError::SourceContactMismatch)?;
                let normal = contact.witness.normal;
                if normal.dot(bundle.propagation().propagation().ray.dir) > 0.0 {
                    normal.scale(-1.0)
                } else { normal }
            } else { target.normal() };
            let value = weighted_irradiance(bundle, normal, projection, radiance,
                cosine, sample.pdf_area_per_m2, log_survival, ratio_zero)?;
            sum.add(value);
            result.images += 1;
        }
        let value = sum.value();
        if !value.is_finite() || value < 0.0 {
            return Err(AreaEmitterError::NumericalRange);
        }
        result.irradiance_w_m2_nm = value;
        Ok(result)
    }
}

impl<'s, 'm> CurvedAreaEmitter<'s, 'm, 's> {
    /// Reuse the exact already-admitted camera emitter when no detector
    /// receiver is needed. This delegates to the same original-source owner.
    pub fn from_camera(
        cx: &Cx<'_>, camera: &CurvedCameraScene<'s, 'm>, object_id: u64, max_triangles: usize,
    ) -> Result<Self, AreaEmitterError> {
        let invalid = || AreaEmitterError::InvalidEmitter { object_id };
        let slot = camera.scene.bindings
            .binary_search_by_key(&object_id, |binding| binding.object_id)
            .map_err(|_| invalid())?;
        let emitter = camera.emitters.get(slot).and_then(|entry| *entry).ok_or_else(invalid)?;
        Self::try_new(cx, camera.scene, emitter, max_triangles)
    }

    /// Conditional original-face support of an already admitted camera source.
    pub fn from_camera_faces(
        cx: &Cx<'_>, camera: &CurvedCameraScene<'s, 'm>, object_id: u64,
        faces: &[usize], limits: AreaEmitterLimits,
    ) -> Result<Self, AreaEmitterError> {
        let invalid = || AreaEmitterError::InvalidEmitter { object_id };
        let slot = camera.scene.bindings
            .binary_search_by_key(&object_id, |binding| binding.object_id)
            .map_err(|_| invalid())?;
        let emitter = camera.emitters.get(slot).and_then(|entry| *entry).ok_or_else(invalid)?;
        Self::try_new_faces(cx, camera.scene, emitter, faces, limits)
    }
}

fn complete(catalog: &BundleImageCatalog) -> Result<(), AreaEmitterError> {
    for (seed, result) in catalog.seed_results().iter().enumerate() {
        if let BundleImageSeedOutcome::Miss(error) = result.outcome() {
            return Err(AreaEmitterError::UnresolvedSeed { seed, error: error.clone() });
        }
    }
    if catalog.images().is_empty() { return Err(AreaEmitterError::EmptyImages); }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn weighted_irradiance(
    bundle: &crate::polarization::dielectric::transport::bundle::PolarizedRayBundle,
    normal: Vec3, projection: ProjectionOptions, radiance: f64, cosine: f64,
    pdf_area: f64, log_survival: f64, ratio_zero: bool,
) -> Result<f64, AreaEmitterError> {
    if ratio_zero {
        return if log_survival == f64::NEG_INFINITY { Ok(0.0) }
            else { Err(AreaEmitterError::NumericalRange) };
    }
    if !log_survival.is_finite() || log_survival > 0.0 {
        return Err(AreaEmitterError::NumericalRange);
    }
    if radiance == 0.0 || cosine == 0.0 || bundle.interface_log_power() == f64::NEG_INFINITY {
        return Ok(0.0);
    }
    let intensity = (radiance * cosine) / pdf_area;
    if intensity.is_finite() && intensity > 0.0 {
        match bundle.point_source_irradiance(normal, projection, intensity) {
            Ok(irradiance) => {
                let value = irradiance * det::exp(log_survival);
                if value.is_finite() && value > 0.0 { return Ok(value); }
            }
            Err(PolarizedBundleError::NumericalRange) => {}
            Err(error) => return Err(AreaEmitterError::Refractive(error)),
        }
    }
    // Combine the same native factors before exponentiation when intermediate
    // intensity or survival cannot be represented. Positive underflow is never
    // silently reported as exact darkness.
    let area = bundle.propagation().project_at_endpoint(normal, projection)
        .map_err(|error| AreaEmitterError::Refractive(error.into()))?
        .area_per_solid_angle_m2();
    let logarithm = det::ln(radiance) + det::ln(cosine) - det::ln(pdf_area)
        + bundle.interface_log_power() - bundle.propagation().propagation().optical_depth
        - det::ln(area) + log_survival;
    let value = det::exp(logarithm);
    if value.is_finite() && value > 0.0 { Ok(value) }
    else { Err(AreaEmitterError::NumericalRange) }
}

#[derive(Default)]
struct Sum { total: f64, correction: f64 }
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.total + value;
        self.correction += if self.total.abs() >= value.abs() {
            (self.total - next) + value
        } else { (value - next) + self.total };
        self.total = next;
    }
    fn value(&self) -> f64 { self.total + self.correction }
}
fn delta(a: [f64; 3], b: [f64; 3]) -> Vec3 {
    Vec3::new(a[0] - b[0], a[1] - b[1], a[2] - b[2])
}
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
}
fn norm(a: Vec3) -> f64 { det::hypot(det::hypot(a.x, a.y), a.z) }

#[cfg(test)]
mod tests;

