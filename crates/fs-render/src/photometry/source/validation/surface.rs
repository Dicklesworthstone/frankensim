//! Direct diffuse luminance at a barycentric point of an ORIGINAL native mesh.
//!
//! This is a surface-sample reference, not an uncertified camera first-hit.
//! Relative barycentric weights locate the point; the original winding and
//! retained rigid placement determine its normal. The observer selects the
//! same two-sided diffuse convention as the native tracer. Both visibility
//! connections exclude only the last ray_epsilon_m at the sampled surface;
//! attenuation always covers their full lengths.
//!
//! Each source contributes I*cos/r² times the spectral sum of source power,
//! rho(lambda)/pi, both transmissions, and the retained photopic weights.
//! Colored reflection is never a gray multiplier on already integrated lux.
//! The Lambertian equation is PBRT 4e, Reflection Models / Diffuse Reflection;
//! extinction of outgoing radiance is its Equation of Transfer, Eq. 14.2.
//!
//! Bound only this direct, once-reflected, unscattered component. Omitted
//! interreflection/in-scattering, emission, calibration, interpolation and
//! physical geometry/material fidelity, and the physical effect of endpoint
//! exclusions remain unknown. No pixel footprint or first-hit claim is made.

use super::{
    ContentHash, Cx, DomainHasher, EnclosureUnavailable, Interval, PhotometricScene,
    PhotometryError, Prepared, cross, dot, finite, hash_float, intensity_bound,
    nonnegative, norm, prepare, spectral_medium::PreparedMedium,
    visibility::{Visibility, world_triangle},
};
use crate::photometry::source::{
    LUMINOUS_EFFICACY_LM_W, SPECTRAL_BINS, photopic_weight_nm, wavelength_nm,
};
use crate::spectral::SpectralReflectance;
use crate::tracer::Material;
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::enclosure::{
    GRID_BUNDLE_ENCLOSURE_SEMANTICS_VERSION, GridEnclosureBudget,
};
use fs_evidence::NumericalCertificate;
use fs_geom::Point3;
use fs_ivl::decay::exp_neg;
use std::{io::{self, Write}, sync::Arc};

type Vector = [Interval; 3];
/// Version of mesh sampling, spectral reflection, visibility and output scope.
pub const DIFFUSE_SURFACE_ENCLOSURE_SEMANTICS_VERSION: u32 = 1;
/// Maximum prepared scalar source responses (4 MiB of initialized intervals).
/// Spectral grids use only a fixed material packet; scene assets remain borrowed.
pub const MAX_DIFFUSE_SURFACE_SOURCES: usize = 262_144;

/// A point on one original triangle, observed from an explicit SI position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiffuseSurfaceSample {
    triangle: usize,
    weights: [f64; 3],
    observer: Point3,
}
impl DiffuseSurfaceSample {
    /// Weights are finite, nonnegative and not all zero. Normalize their exact
    /// real sum, not a rounded sum assumed to equal one. [1,1,1] is a centroid.
    /// Triangle indices refer to the primitive selected by the reference.
    pub fn try_new(triangle: usize, relative_weights: [f64; 3], observer_m: Point3)
        -> Result<Self, SurfaceEnclosureError> {
        if relative_weights.iter().any(|v| !v.is_finite() || *v < 0.0)
            || relative_weights.iter().all(|v| *v == 0.0)
            || [observer_m.x, observer_m.y, observer_m.z].iter().any(|v| !v.is_finite()) {
            return Err(SurfaceEnclosureError::InvalidSample);
        }
        Ok(Self { triangle, weights: relative_weights.map(|v| if v == 0.0 { 0.0 } else { v }), observer: observer_m })
    }
    /// Index in the original primitive's triangle array.
    #[must_use]
    pub const fn triangle(self) -> usize { self.triangle }
    /// Original relative weights; no rounded point replaces this declaration.
    #[must_use]
    pub const fn relative_weights(self) -> [f64; 3] { self.weights }
    /// Position at which outgoing radiance is observed, in metres.
    #[must_use]
    pub const fn observer_m(self) -> Point3 { self.observer }
}

/// Explicit refusal; a failed material, visibility or grid bound is not black.
#[derive(Debug)]
pub enum SurfaceEnclosureError {
    /// Invalid primitive/triangle/weights/observer or degenerate surface geometry.
    InvalidSample,
    /// The selected original material is not passive tabulated SpectralDiffuse.
    UnsupportedMaterial,
    /// The observer's side of the surface plane could not be proved.
    UnresolvedObserverSide,
    /// Storage/output limit exceeded or zero grid-node budget.
    Budget,
    /// Bounded storage allocation failed.
    Allocation,
    /// Intermediate or final outward bounds are not finite/representable.
    NumericalRange,
    /// Original optional source/mesh preparation could not establish an enclosure.
    Unavailable(EnclosureUnavailable),
    /// Original native medium refused the full connection or its work budget.
    PhysicalMedium(VolumeTransportError),
    /// Cancellation or original distance/geometry refusal.
    Transport(PhotometryError),
}
impl core::fmt::Display for SurfaceEnclosureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "diffuse surface enclosure refused: {self:?}")
    }
}
impl std::error::Error for SurfaceEnclosureError {}
impl From<PhotometryError> for SurfaceEnclosureError {
    fn from(value: PhotometryError) -> Self { Self::Transport(value) }
}
impl From<fs_exec::Cancelled> for SurfaceEnclosureError {
    fn from(_: fs_exec::Cancelled) -> Self { Self::Transport(PhotometryError::Cancelled) }
}
impl From<VolumeTransportError> for SurfaceEnclosureError {
    fn from(value: VolumeTransportError) -> Self {
        match value {
            VolumeTransportError::Cancelled => Self::Transport(PhotometryError::Cancelled),
            other => Self::PhysicalMedium(other),
        }
    }
}

/// Component-scoped reference at an original mesh sample, before exposure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnclosedDiffuseLuminance {
    /// Midpoint of the numerical enclosure, in cd/m²; not sampled image radiance.
    pub luminance_cd_m2: f64,
    /// Encloses the retained direct diffuse numerical model, including visibility.
    pub numerical: NumericalCertificate,
    /// Outward absolute error about the reported midpoint, in cd/m².
    pub numerical_error_bound_cd_m2: f64,
    /// Native material/mesh primitive selected by the reference.
    pub primitive: usize,
    /// Native triangle selected by this observation.
    pub triangle: usize,
    /// Potentially nonzero source connections whose visibility is unresolved.
    pub uncertain_light_connections: usize,
    /// True only when a potentially nonzero observer connection is unresolved.
    pub uncertain_observer_visibility: bool,
    /// Actual geometry, source bins, material knots, medium, observer and sample.
    pub numerical_model_identity: ContentHash,
    /// Calibration, model discrepancy and omitted transport remain unknown.
    pub total_error_bound_cd_m2: Option<f64>,
}
impl EnclosedDiffuseLuminance {
    /// Write a separate surface-component CSV header, not a full-image certificate.
    pub fn write_csv_header(writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(b"index,scope,primitive,triangle,luminance_cd_m2,numerical_lower_cd_m2,numerical_upper_cd_m2,numerical_error_bound_cd_m2,uncertain_light_connections,uncertain_observer_visibility,numerical_model_identity,total_physical_error_bound_cd_m2,version\n")
    }
    /// Stream round-trip values and explicitly unknown total physical error.
    pub fn write_csv_row(self, index: usize, writer: &mut impl Write) -> io::Result<()> {
        writeln!(writer, "{index},direct-diffuse-surface,{},{},{},{},{},{},{},{},{},unknown,{DIFFUSE_SURFACE_ENCLOSURE_SEMANTICS_VERSION}",
            self.primitive, self.triangle, self.luminance_cd_m2, self.numerical.lo,
            self.numerical.hi, self.numerical_error_bound_cd_m2, self.uncertain_light_connections,
            self.uncertain_observer_visibility, self.numerical_model_identity)
    }
}

/// Prepared response of one native material over shared scene/medium metadata.
/// Reflections and attenuation retain the ORIGINAL source spectra and material.
/// Queries allocate nothing. A sibling primitive can share the same prepared
/// geometry and medium through `for_primitive`, without a second scene copy.
pub struct DiffuseSurfaceReference<'a, 'scene> {
    scene: &'a PhotometricScene<'scene>,
    prepared: Arc<Prepared<'scene>>,
    primitive: usize,
    weights: [Interval; SPECTRAL_BINS],
    gray_response: Option<Vec<Interval>>,
    grid_budget: GridEnclosureBudget,
    identity: ContentHash,
}
impl<'a, 'scene> DiffuseSurfaceReference<'a, 'scene> {
    /// Select an actual native SpectralDiffuse material and prepare its scene.
    /// Other opaque materials may occlude; refractive geometry explicitly refuses.
    pub fn new(scene: &'a PhotometricScene<'scene>, primitive: usize, cx: &Cx<'_>)
        -> Result<Self, SurfaceEnclosureError> {
        cx.checkpoint()?;
        material(scene, primitive)?;
        if scene.sources.len() > MAX_DIFFUSE_SURFACE_SOURCES { return Err(SurfaceEnclosureError::Budget); }
        for object in scene.primitives {
            cx.checkpoint()?;
            if matches!(object.material, Material::Dielectric { .. }) {
                return Err(SurfaceEnclosureError::Unavailable(EnclosureUnavailable::GeometryVisibility));
            }
        }
        let prepared = Arc::new(prepare(scene, cx)?);
        Self::bind(scene, prepared, primitive, GridEnclosureBudget::default(), cx)
    }

    /// Select another primitive while sharing original visibility/extinction metadata.
    pub fn for_primitive(&self, primitive: usize, cx: &Cx<'_>) -> Result<Self, SurfaceEnclosureError> {
        Self::bind(self.scene, Arc::clone(&self.prepared), primitive, self.grid_budget, cx)
    }

    /// Set additional native-grid enclosure work. The original cell cap still applies.
    /// Immutable references/jobs retain their own prior budget.
    pub fn with_grid_budget(mut self, budget: GridEnclosureBudget) -> Result<Self, SurfaceEnclosureError> {
        if budget.max_nodes_per_segment == 0 { return Err(SurfaceEnclosureError::Budget); }
        self.grid_budget = budget;
        Ok(self)
    }

    fn bind(scene: &'a PhotometricScene<'scene>, prepared: Arc<Prepared<'scene>>, primitive: usize,
        grid_budget: GridEnclosureBudget, cx: &Cx<'_>) -> Result<Self, SurfaceEnclosureError> {
        cx.checkpoint()?;
        let reflectance = material(scene, primitive)?;
        let (model, medium) = match prepared.as_ref() {
            Prepared::Available { identity, medium, .. } => (*identity, medium),
            Prepared::Unavailable(reason) => return Err(SurfaceEnclosureError::Unavailable(*reason)),
        };
        let pi = Interval::new(core::f64::consts::PI.next_down(), core::f64::consts::PI.next_up());
        let mut weights = [zero(); SPECTRAL_BINS];
        let mut hash = DomainHasher::new("org.frankensim.photometry.diffuse-surface.v1");
        hash.update(&model.0);
        hash.update(&DIFFUSE_SURFACE_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&GRID_BUNDLE_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&(primitive as u64).to_le_bytes());
        hash.update(&(reflectance.knots().len() as u64).to_le_bytes());
        for &(lambda, rho) in reflectance.knots() { hash_float(&mut hash, lambda); hash_float(&mut hash, rho); }
        let observer = crate::photometry::observer::OBSERVER_ID.as_bytes();
        hash.update(&(observer.len() as u64).to_le_bytes()); hash.update(observer);
        hash_float(&mut hash, LUMINOUS_EFFICACY_LM_W);
        for (bin, weight) in weights.iter_mut().enumerate() {
            cx.checkpoint()?;
            let lambda = wavelength_nm(bin);
            let rho = reflectance_enclosure(reflectance, lambda);
            *weight = if rho.hi() == 0.0 { zero() } else {
                nonnegative(mul(mul(rho, point(photopic_weight_nm(bin))), point(LUMINOUS_EFFICACY_LM_W)) / pi)
            };
            if !finite(*weight) { return Err(SurfaceEnclosureError::NumericalRange); }
            hash_float(&mut hash, lambda); hash_float(&mut hash, photopic_weight_nm(bin));
        }
        // Gray transport still needs the ACTUAL spectral material response.
        // Cache only one scalar per source, not a copied source x material packet.
        let mut gray_response = if medium.needs_spectrum() { None } else { Some(Vec::new()) };
        if let Some(values) = &mut gray_response {
            values.try_reserve_exact(scene.sources.len()).map_err(|_| SurfaceEnclosureError::Allocation)?;
        }
        for source in scene.sources {
            cx.checkpoint()?;
            let mut response = zero();
            for (bin, &power) in source.spectrum().bins_w_per_lumen_nm().iter().enumerate() {
                hash_float(&mut hash, power);
                if gray_response.is_some() { response = add(response, mul(point(power), weights[bin])); }
            }
            if !finite(response) { return Err(SurfaceEnclosureError::NumericalRange); }
            if let Some(values) = &mut gray_response { values.push(response); }
        }
        cx.checkpoint()?;
        Ok(Self { scene, prepared, primitive, weights, gray_response, grid_budget, identity: hash.finalize() })
    }

    /// Bound direct diffuse radiance at an original barycentric mesh location.
    /// Observer distance changes extinction/visibility, never a fictitious 1/r²
    /// radiance loss. A side unresolved by outward predicates refuses explicitly.
    pub fn enclose(&self, cx: &Cx<'_>, sample: DiffuseSurfaceSample)
        -> Result<EnclosedDiffuseLuminance, SurfaceEnclosureError> {
        cx.checkpoint()?;
        let (position, mut normal) = sample_geometry(&self.scene.primitives[self.primitive], sample)?;
        let observer = [sample.observer.x, sample.observer.y, sample.observer.z].map(point);
        let to_observer = subtract(observer, position);
        let observer_distance = distance(to_observer, self.scene.options.ray_epsilon_m, self.scene.options.max_distance_m)?;
        let side = dot(normal, to_observer);
        if side.lo() > 0.0 {} else if side.hi() < 0.0 { normal = normal.map(|v| -v); }
        else { return Err(SurfaceEnclosureError::UnresolvedObserverSide); }
        let (visibility, medium) = match self.prepared.as_ref() {
            Prepared::Available { visibility, medium, .. } => (visibility, medium),
            Prepared::Unavailable(reason) => return Err(SurfaceEnclosureError::Unavailable(*reason)),
        };
        let camera_visibility = visibility.visibility(self.scene.primitives, observer, to_observer.map(|v| -v),
            visibility_end(observer_distance, self.scene.options.ray_epsilon_m), cx)?;
        let mut value = zero();
        let mut uncertain = 0;
        if camera_visibility != Visibility::Blocked && self.weights.iter().any(|v| v.hi() > 0.0)
            && self.scene.sources.iter().any(|s| s.flux_lm() > 0.0) {
            let camera_packet = match medium {
                PreparedMedium::Grid { packet, .. } => Some(packet.transmittance_box(cx, observer, position, self.grid_budget)?),
                PreparedMedium::Gray { .. } => None,
            };
            for (index, source) in self.scene.sources.iter().enumerate() {
                cx.checkpoint()?;
                if source.flux_lm() == 0.0 { continue; }
                let origin = source.transform().translation_m().map(point);
                let delta = subtract(origin, position);
                let length = distance(delta, self.scene.options.ray_epsilon_m, self.scene.options.max_distance_m)?;
                let direction = delta.map(|v| v / length);
                let cosine = nonnegative(dot(normal, direction));
                if cosine.hi() == 0.0 { continue; }
                let intensity = intensity_bound(source, direction.map(|v| -v), cx)?
                    .ok_or(SurfaceEnclosureError::NumericalRange)?;
                if intensity.hi() == 0.0 { continue; }
                let visible = visibility.visibility(self.scene.primitives, origin, delta.map(|v| -v),
                    visibility_end(length, self.scene.options.ray_epsilon_m), cx)?;
                if visible == Visibility::Blocked { continue; }
                let response = match medium {
                    PreparedMedium::Gray { extinction, vacuum } => {
                        let transmission = if *vacuum { point(1.0) } else {
                            exp_neg(mul(*extinction, add(length, observer_distance)))
                                .ok_or(SurfaceEnclosureError::NumericalRange)?
                        };
                        mul(self.gray_response.as_ref().ok_or(SurfaceEnclosureError::NumericalRange)?[index], transmission)
                    }
                    PreparedMedium::Grid { packet, .. } => {
                        let light = packet.transmittance_box(cx, origin, position, self.grid_budget)?;
                        let camera = camera_packet.as_ref().ok_or(SurfaceEnclosureError::NumericalRange)?;
                        let mut response = zero();
                        for (bin, &power) in source.spectrum().bins_w_per_lumen_nm().iter().enumerate() {
                            if bin % 16 == 0 { cx.checkpoint()?; }
                            response = add(response, mul(mul(point(power), self.weights[bin]), mul(light[bin], camera[bin])));
                        }
                        response
                    }
                };
                let mut contribution = mul(mul(nonnegative(intensity / length / length), cosine), response);
                if !finite(contribution) { return Err(SurfaceEnclosureError::NumericalRange); }
                if visible == Visibility::Uncertain && contribution.hi() > 0.0 {
                    contribution = Interval::new(0.0, contribution.hi()); uncertain += 1;
                }
                value = add(value, contribution);
                if !finite(value) { return Err(SurfaceEnclosureError::NumericalRange); }
            }
        }
        let observer_uncertain = camera_visibility == Visibility::Uncertain && value.hi() > 0.0;
        if observer_uncertain { value = Interval::new(0.0, value.hi()); }
        let midpoint = value.lo() * 0.5 + value.hi() * 0.5;
        let error = if value.lo() == value.hi() { 0.0 } else {
            (point(midpoint) - point(value.lo())).hi().max((point(value.hi()) - point(midpoint)).hi())
        };
        if !error.is_finite() { return Err(SurfaceEnclosureError::NumericalRange); }
        let mut hash = DomainHasher::new("org.frankensim.photometry.diffuse-surface.sample.v1");
        hash.update(&self.identity.0); hash.update(&(sample.triangle as u64).to_le_bytes());
        for x in sample.weights.into_iter().chain([sample.observer.x, sample.observer.y, sample.observer.z]) {
            hash_float(&mut hash, x);
        }
        cx.checkpoint()?;
        Ok(EnclosedDiffuseLuminance {
            luminance_cd_m2: midpoint, numerical: NumericalCertificate::enclosure(value.lo(), value.hi()),
            numerical_error_bound_cd_m2: error, primitive: self.primitive, triangle: sample.triangle,
            uncertain_light_connections: uncertain, uncertain_observer_visibility: observer_uncertain,
            numerical_model_identity: hash.finalize(), total_error_bound_cd_m2: None,
        })
    }
}

/// Borrowed surface observations with atomic publication of each completed value.
pub struct DiffuseSurfaceJob<'r, 'a, 'scene> {
    reference: &'r DiffuseSurfaceReference<'a, 'scene>,
    samples: &'r [DiffuseSurfaceSample],
    results: Vec<EnclosedDiffuseLuminance>,
}
impl<'r, 'a, 'scene> DiffuseSurfaceJob<'r, 'a, 'scene> {
    /// Admit output storage; neither mesh data nor observations are copied.
    pub fn new(reference: &'r DiffuseSurfaceReference<'a, 'scene>, samples: &'r [DiffuseSurfaceSample], cx: &Cx<'_>)
        -> Result<Self, SurfaceEnclosureError> {
        cx.checkpoint()?;
        if samples.len() > reference.scene.options.max_outputs { return Err(SurfaceEnclosureError::Budget); }
        let mut results = Vec::new();
        results.try_reserve_exact(samples.len()).map_err(|_| SurfaceEnclosureError::Allocation)?;
        cx.checkpoint()?;
        Ok(Self { reference, samples, results })
    }
    /// Chunking does not enter the arithmetic. Failure retains the accepted prefix.
    pub fn advance(&mut self, cx: &Cx<'_>, max_samples: usize) -> Result<usize, SurfaceEnclosureError> {
        cx.checkpoint()?;
        let start = self.results.len();
        let end = start.saturating_add(max_samples).min(self.samples.len());
        for &sample in &self.samples[start..end] {
            let result = self.reference.enclose(cx, sample)?;
            cx.checkpoint()?; self.results.push(result);
        }
        Ok(self.results.len() - start)
    }
    /// Complete accepted prefix; unfinished observations are absent, not black.
    #[must_use]
    pub fn results(&self) -> &[EnclosedDiffuseLuminance] { &self.results }
    /// Whether every declared observation completed.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.results.len() == self.samples.len() }
}

fn material<'a>(scene: &'a PhotometricScene<'_>, primitive: usize) -> Result<&'a SpectralReflectance, SurfaceEnclosureError> {
    match &scene.primitives.get(primitive).ok_or(SurfaceEnclosureError::InvalidSample)?.material {
        Material::SpectralDiffuse { reflectance } => Ok(reflectance),
        _ => Err(SurfaceEnclosureError::UnsupportedMaterial),
    }
}
fn reflectance_enclosure(curve: &SpectralReflectance, lambda: f64) -> Interval {
    let knots = curve.knots();
    let right = knots.partition_point(|&(wavelength, _)| wavelength < lambda);
    if knots[right].0 == lambda { return point(knots[right].1); }
    let (x, a) = knots[right - 1]; let (y, b) = knots[right];
    if a == b { return point(a); }
    let t = (point(lambda) - point(x)) / (point(y) - point(x));
    let t = Interval::new(t.lo().max(0.0), t.hi().min(1.0));
    let u = point(1.0) - t;
    let u = Interval::new(u.lo().max(0.0), u.hi().min(1.0));
    let result = add(mul(point(a), u), mul(point(b), t));
    Interval::new(result.lo().max(a.min(b)), result.hi().min(a.max(b)))
}
fn sample_geometry(primitive: &crate::tracer::Primitive, sample: DiffuseSurfaceSample)
    -> Result<(Vector, Vector), SurfaceEnclosureError> {
    let [a, b, c] = world_triangle(primitive, sample.triangle).ok_or(SurfaceEnclosureError::InvalidSample)?;
    let edges = [subtract(b, a), subtract(c, a)];
    let scale = edges.iter().flatten().map(|v| v.lo().abs().max(v.hi().abs())).fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale == 0.0 { return Err(SurfaceEnclosureError::InvalidSample); }
    let [e, f] = edges.map(|edge| edge.map(|v| v / point(scale)));
    let raw = cross(e, f);
    let size = norm(raw).ok_or(SurfaceEnclosureError::NumericalRange)?;
    if size.lo() <= 0.0 { return Err(SurfaceEnclosureError::InvalidSample); }
    let normal = raw.map(|v| v / size);
    let peak = sample.weights.into_iter().fold(0.0_f64, f64::max);
    let weights = sample.weights.map(|v| if v == 0.0 { zero() } else { nonnegative(point(v) / point(peak)) });
    let total = add(add(weights[0], weights[1]), weights[2]);
    if total.lo() <= 0.0 || !finite(total) { return Err(SurfaceEnclosureError::NumericalRange); }
    let [_, u, v] = weights.map(|w| if w.hi() == 0.0 { zero() } else { nonnegative(w / total) });
    let position = std::array::from_fn(|i| a[i] + edges[0][i] * u + edges[1][i] * v);
    if position.iter().chain(&normal).any(|v| !finite(*v)) { return Err(SurfaceEnclosureError::NumericalRange); }
    Ok((position, normal))
}
fn distance(vector: Vector, minimum: f64, maximum: f64) -> Result<Interval, SurfaceEnclosureError> {
    let length = norm(vector).ok_or(SurfaceEnclosureError::NumericalRange)?;
    if length.lo() <= minimum { return Err(SurfaceEnclosureError::InvalidSample); }
    if length.hi() > maximum { return Err(PhotometryError::DistanceLimitExceeded.into()); }
    Ok(length)
}
fn visibility_end(length: Interval, epsilon: f64) -> Interval {
    let end = point(1.0) - point(epsilon) / length;
    Interval::new(end.lo().clamp(0.0, 1.0), end.hi().clamp(0.0, 1.0))
}
fn subtract(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|i| if a[i].lo() == a[i].hi() && a[i] == b[i] { zero() } else { a[i] - b[i] })
}
fn point(x: f64) -> Interval { Interval::point(x) }
fn zero() -> Interval { point(0.0) }
fn add(a: Interval, b: Interval) -> Interval {
    if a.hi() == 0.0 { b } else if b.hi() == 0.0 { a } else { nonnegative(a + b) }
}
fn mul(a: Interval, b: Interval) -> Interval {
    if a.hi() == 0.0 || b.hi() == 0.0 { zero() } else { nonnegative(a * b) }
}

#[cfg(test)]
mod tests;
