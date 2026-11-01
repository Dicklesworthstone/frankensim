//! Bounded single-scattered haze radiance along an explicit camera segment.
//!
//! For q(t)=camera+t*(end-camera), integrate, in cd/m²,
//! sigma_s * |end-camera| * sum_i I_i(q-source_i) p(mu_i)
//! * exp(-sigma_t*(|q-source_i|+t*|end-camera|)) / |q-source_i|² dt.
//! Both directions in mu describe photon propagation, as in the original phase
//! owner. The shared source-law and opaque-mesh interval evaluators bound EVERY
//! point of a panel. A narrow beam or shadow missed by midpoint samples cannot
//! disappear. Both camera and light connections retain uncertain visibility.
//!
//! The endpoints are original SI points; no rounded camera normalization changes
//! the segment. The camera leg has no endpoint offset. The light leg retains the
//! scene's receiver-only visibility exclusion, which never shortens attenuation.
//! Occluded camera intervals contribute zero; no nominal first-hit distance is
//! promoted into an exact boundary. Refine the largest radiance uncertainty.
//!
//! This is the single-scattering term of the original homogeneous gray model,
//! not a total image or a replacement path tracer. Native nominal calculations
//! are untouched. Higher scattering orders, surface terms, thermal emission,
//! finite-source optics, calibration, spectral/field/model discrepancy and the
//! physical effect of endpoint exclusions remain unbounded. The transport
//! equation is described in PBRT 4e, Light Transport II, Equation of Transfer;
//! use the native propagation-direction convention, not PBRT's away vectors.

use super::{
    ContentHash, Cx, DomainHasher, EnclosureUnavailable, Interval, PhotometricScene,
    PhotometryError, Prepared, dot, finite, hash_float, intensity_bound,
    nonnegative, norm, prepare, visibility::{Visibility, VisibilityIndex},
};
use crate::photometry::medium::PhaseFunction;
use fs_evidence::NumericalCertificate;
use fs_geom::Point3;
use fs_ivl::decay::exp_neg;
use std::{cmp::Ordering, collections::BinaryHeap, io::{self, Write}};

mod phase;

/// Version of the target integral, phase enclosure, visibility and output scope.
pub const SCATTER_ENCLOSURE_SEMANTICS_VERSION: u32 = 1;
/// Hard ceiling on source-panel evaluations and retained frontier entries.
pub const MAX_SCATTER_ENCLOSURE_PANELS: usize = 262_144;

type Vector = [Interval; 3];

/// Original camera and far endpoint in metres, defining a finite straight ray.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScatterSegment { camera: Point3, end: Point3 }
impl ScatterSegment {
    /// Admit finite distinct endpoints. Numeric length and source separation
    /// are checked with outward arithmetic when this segment is evaluated.
    pub fn try_new(camera: Point3, end: Point3) -> Result<Self, ScatterEnclosureError> {
        if [camera.x, camera.y, camera.z, end.x, end.y, end.z].into_iter().any(|v| !v.is_finite())
            || camera == end { return Err(ScatterEnclosureError::InvalidSegment); }
        Ok(Self { camera, end })
    }
    /// Camera position, in metres.
    #[must_use]
    pub const fn camera_m(self) -> Point3 { self.camera }
    /// Far endpoint, in metres; surfaces may occlude part of this segment.
    #[must_use]
    pub const fn end_m(self) -> Point3 { self.end }
}

/// One common numerical-error/work budget for phase, beam, shadows and extinction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScatterEnclosureBudget {
    /// Finite nonnegative requested error about the returned midpoint, in cd/m².
    pub absolute_error_cd_m2: f64,
    /// Root sources plus two children per bisection; never exceeds the hard cap.
    pub max_evaluations: usize,
    /// Binary subdivision depth, at most 48. Zero retains whole-segment bounds.
    pub max_depth: u32,
}
impl Default for ScatterEnclosureBudget {
    fn default() -> Self {
        Self { absolute_error_cd_m2: 1e-4, max_evaluations: 16_384, max_depth: 40 }
    }
}

/// Explicit inability to enclose a segment; no failed value becomes black.
#[derive(Debug)]
pub enum ScatterEnclosureError {
    /// Nonfinite/coincident segment endpoints or unresolved segment length.
    InvalidSegment,
    /// Invalid accuracy/work limits or insufficient room for root sources/outputs.
    Budget,
    /// Allocation of bounded frontier/output storage failed.
    Allocation,
    /// Physical-grid media need a separate heterogeneous in-scattering enclosure.
    HomogeneousOnly,
    /// The phase law is not a wavelength-independent admitted native law.
    UnsupportedPhase,
    /// Separation from a point-source singularity was not proved for the whole ray.
    UnresolvedSourceSeparation {
        /// Source in the original scene order.
        source: usize,
    },
    /// A finite useful arithmetic enclosure could not be retained.
    NumericalRange,
    /// Original optional source/visibility preparation was unavailable.
    Unavailable(EnclosureUnavailable),
    /// Original cancellation, geometry or declared-distance refusal.
    Transport(PhotometryError),
}
impl core::fmt::Display for ScatterEnclosureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "single-scattering enclosure refused: {self:?}")
    }
}
impl std::error::Error for ScatterEnclosureError {}
impl From<PhotometryError> for ScatterEnclosureError {
    fn from(e: PhotometryError) -> Self { Self::Transport(e) }
}
impl From<fs_exec::Cancelled> for ScatterEnclosureError {
    fn from(_: fs_exec::Cancelled) -> Self { Self::Transport(PhotometryError::Cancelled) }
}

/// Every termination preserves the bound; only the first establishes tolerance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScatterEnclosureTermination {
    /// Outward midpoint error is no greater than the requested absolute error.
    ToleranceMet,
    /// The next two children would exceed the common evaluation budget.
    PanelBudget,
    /// No refinable panel remains at the admitted depth or arithmetic resolution.
    ResolutionLimit,
}

/// A numerical reference for the single-scattered component, before exposure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnclosedSingleScatter {
    /// Deterministic enclosure midpoint, not a Monte Carlo estimate.
    pub luminance_cd_m2: f64,
    /// Outward enclosure of this explicitly scoped numerical integral.
    pub numerical: NumericalCertificate,
    /// Error about the midpoint, including unresolved shadows and beam boundaries.
    pub numerical_error_bound_cd_m2: f64,
    /// Requested accuracy, distinct from the achieved bound.
    pub requested_error_cd_m2: f64,
    /// Why work stopped; budget exhaustion is not a tolerance pass.
    pub termination: ScatterEnclosureTermination,
    /// Evaluated source-panels, including discarded parents and analytic zeros.
    pub evaluations: usize,
    /// Final possibly nonzero panels with unresolved camera or light visibility.
    pub uncertain_visibility_panels: usize,
    /// Original source/geometry/gray coefficients/phase, endpoints and semantics.
    pub numerical_model_identity: ContentHash,
    /// Higher scattering orders remain unknown, not zero.
    pub omitted_multiple_scattering_bound_cd_m2: Option<f64>,
    /// Calibration, construction, geometry/model discrepancy and total error unknown.
    pub total_error_bound_cd_m2: Option<f64>,
}
impl EnclosedSingleScatter {
    /// Separate component-scoped CSV, not a certificate for a full image.
    pub fn write_csv_header(writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(b"index,scope,luminance_cd_m2,numerical_lower_cd_m2,numerical_upper_cd_m2,numerical_error_bound_cd_m2,requested_error_cd_m2,termination,evaluations,uncertain_visibility_panels,numerical_model_identity,omitted_multiple_scattering_bound_cd_m2,total_physical_error_bound_cd_m2,version\n")
    }
    /// Stream round-trip values with explicit unknown omitted and physical errors.
    pub fn write_csv_row(self, index: usize, writer: &mut impl Write) -> io::Result<()> {
        let status = match self.termination {
            ScatterEnclosureTermination::ToleranceMet => "tolerance-met",
            ScatterEnclosureTermination::PanelBudget => "panel-budget",
            ScatterEnclosureTermination::ResolutionLimit => "resolution-limit",
        };
        writeln!(writer, "{index},homogeneous-gray-single-scatter,{},{},{},{},{},{status},{},{},{},unknown,unknown,{SCATTER_ENCLOSURE_SEMANTICS_VERSION}",
            self.luminance_cd_m2, self.numerical.lo, self.numerical.hi,
            self.numerical_error_bound_cd_m2, self.requested_error_cd_m2,
            self.evaluations, self.uncertain_visibility_panels, self.numerical_model_identity)
    }
}

/// Prepared native source and visibility metadata, shared by all camera segments.
/// No geometry, source table, spectrum or phase table is independently reconstructed.
pub struct HomogeneousScatterReference<'a, 'scene> {
    scene: &'a PhotometricScene<'scene>,
    visibility: VisibilityIndex,
    extinction: Interval,
    scattering: f64,
    phase: &'a PhaseFunction,
    identity: ContentHash,
}
impl<'a, 'scene> HomogeneousScatterReference<'a, 'scene> {
    /// Admit the existing gray medium and supported opaque geometry/source laws.
    pub fn new(scene: &'a PhotometricScene<'scene>, cx: &Cx<'_>) -> Result<Self, ScatterEnclosureError> {
        cx.checkpoint()?;
        let medium = scene.medium.homogeneous().ok_or(ScatterEnclosureError::HomogeneousOnly)?;
        let phase = medium.phase_ref();
        if matches!(phase, PhaseFunction::SpectralTable(_) | PhaseFunction::SpectralMixture(_)) {
            return Err(ScatterEnclosureError::UnsupportedPhase);
        }
        for primitive in scene.primitives {
            cx.checkpoint()?;
            if matches!(primitive.material, crate::tracer::Material::Dielectric { .. }) {
                return Err(ScatterEnclosureError::Unavailable(EnclosureUnavailable::GeometryVisibility));
            }
        }
        let (identity, visibility) = match prepare(scene, cx)? {
            Prepared::Available { identity, visibility, .. } => (identity, visibility),
            Prepared::Unavailable(reason) => return Err(ScatterEnclosureError::Unavailable(reason)),
        };
        let mut hash = DomainHasher::new("org.frankensim.photometry.single-scatter.v1");
        hash.update(&SCATTER_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&identity.0);
        phase.hash_into(&mut hash);
        cx.checkpoint()?;
        Ok(Self { scene, visibility,
            extinction: add(point(medium.absorption_per_m()), point(medium.scattering_per_m())),
            scattering: medium.scattering_per_m(), phase, identity: hash.finalize() })
    }

    /// Enclose one camera segment. Jobs amortize frontier allocation across rays.
    pub fn enclose(&self, cx: &Cx<'_>, segment: ScatterSegment, budget: ScatterEnclosureBudget)
        -> Result<EnclosedSingleScatter, ScatterEnclosureError> {
        cx.checkpoint()?;
        let mut workspace = Workspace::new(budget, self.scene.sources.len())?;
        self.evaluate(cx, segment, budget, &mut workspace)
    }

    fn evaluate(&self, cx: &Cx<'_>, segment: ScatterSegment, budget: ScatterEnclosureBudget,
        work: &mut Workspace) -> Result<EnclosedSingleScatter, ScatterEnclosureError> {
        cx.checkpoint()?;
        work.frontier.clear();
        let p = segment.camera;
        let camera = coordinates(p);
        let delta = difference(coordinates(segment.end), camera);
        let length = norm(delta).ok_or(ScatterEnclosureError::InvalidSegment)?;
        if length.lo() <= 0.0 { return Err(ScatterEnclosureError::InvalidSegment); }
        if length.hi() > self.scene.options.max_distance_m {
            return Err(PhotometryError::DistanceLimitExceeded.into());
        }
        let ray = SegmentGeometry { camera, delta, length };
        let mut raw = zero();
        let mut evaluations = 0;
        let mut uncertain = 0;
        for (source, fixture) in self.scene.sources.iter().enumerate() {
            cx.checkpoint()?;
            let mut panel = Panel { source, lo: 0.0, hi: 1.0, minimum_distance: 0.0,
                depth: 0, order: evaluations, range: zero(), uncertain: false };
            if self.scattering > 0.0 && fixture.flux_lm() > 0.0 {
                let offset = difference(coordinates(fixture.position_m()), camera);
                // Convex squared distance has its minimum at the clamped
                // projection. Enclose that projection, not a sampled near miss.
                let projection = dot(offset, delta) / (length * length);
                if !finite(projection) { return Err(ScatterEnclosureError::NumericalRange); }
                let t = Interval::new(projection.lo().clamp(0.0, 1.0), projection.hi().clamp(0.0, 1.0));
                let closest = difference(offset, delta.map(|v| v * t));
                panel.minimum_distance = norm(closest).ok_or(ScatterEnclosureError::NumericalRange)?.lo();
                if panel.minimum_distance <= self.scene.options.ray_epsilon_m {
                    return Err(ScatterEnclosureError::UnresolvedSourceSeparation { source });
                }
                let far = norm(offset).ok_or(ScatterEnclosureError::NumericalRange)?.hi()
                    .max(norm(difference(offset, delta)).ok_or(ScatterEnclosureError::NumericalRange)?.hi());
                if far > self.scene.options.max_distance_m { return Err(PhotometryError::DistanceLimitExceeded.into()); }
                self.panel(cx, &ray, &mut panel)?;
            }
            evaluations += 1;
            raw = add(raw, panel.range);
            if !finite(raw) { return Err(ScatterEnclosureError::NumericalRange); }
            uncertain += usize::from(panel.uncertain);
            work.retain(panel, budget);
        }
        let mut enclosed = raw;
        let termination = loop {
            cx.checkpoint()?;
            if midpoint_error(enclosed).1 <= budget.absolute_error_cd_m2 { break ScatterEnclosureTermination::ToleranceMet; }
            if work.frontier.is_empty() { break ScatterEnclosureTermination::ResolutionLimit; }
            if evaluations + 2 > budget.max_evaluations { break ScatterEnclosureTermination::PanelBudget; }
            let parent = work.frontier.pop().ok_or(ScatterEnclosureError::NumericalRange)?;
            let middle = parent.lo + (parent.hi - parent.lo) * 0.5;
            let mut left = Panel { hi: middle, depth: parent.depth+1, order: evaluations,
                range: zero(), uncertain: false, ..parent };
            let mut right = Panel { lo: middle, depth: parent.depth+1, order: evaluations+1,
                range: zero(), uncertain: false, ..parent };
            self.panel(cx, &ray, &mut left)?;
            self.panel(cx, &ray, &mut right)?;
            evaluations += 2;
            // Keep bounds on the sums of LEAF ENDPOINTS separate from the
            // intersected truth enclosure; only the former permits replacement.
            let lo = ((point(raw.lo()) - point(parent.range.lo())) + point(left.range.lo()) + point(right.range.lo())).lo().max(0.0);
            let hi = ((point(raw.hi()) - point(parent.range.hi())) + point(left.range.hi()) + point(right.range.hi())).hi();
            if !lo.is_finite() || !hi.is_finite() || lo > hi { return Err(ScatterEnclosureError::NumericalRange); }
            raw = Interval::new(lo, hi);
            enclosed = enclosed.intersect(raw).ok_or(ScatterEnclosureError::NumericalRange)?;
            uncertain -= usize::from(parent.uncertain);
            for child in [left, right] {
                uncertain += usize::from(child.uncertain);
                work.retain(child, budget);
            }
        };
        let (value, error) = midpoint_error(enclosed);
        if !error.is_finite() { return Err(ScatterEnclosureError::NumericalRange); }
        let mut hash = DomainHasher::new("org.frankensim.photometry.single-scatter.segment.v1");
        hash.update(&self.identity.0);
        for v in [segment.camera.x, segment.camera.y, segment.camera.z, segment.end.x, segment.end.y, segment.end.z] {
            hash_float(&mut hash, v);
        }
        cx.checkpoint()?;
        Ok(EnclosedSingleScatter { luminance_cd_m2: value,
            numerical: NumericalCertificate::enclosure(enclosed.lo(), enclosed.hi()),
            numerical_error_bound_cd_m2: error, requested_error_cd_m2: budget.absolute_error_cd_m2,
            termination, evaluations, uncertain_visibility_panels: uncertain,
            numerical_model_identity: hash.finalize(), omitted_multiple_scattering_bound_cd_m2: None,
            total_error_bound_cd_m2: None })
    }

    fn panel(&self, cx: &Cx<'_>, ray: &SegmentGeometry, panel: &mut Panel) -> Result<(), ScatterEnclosureError> {
        cx.checkpoint()?;
        let t = Interval::new(panel.lo, panel.hi);
        let camera_path = ray.delta.map(|v| v * t);
        let position = std::array::from_fn(|i| ray.camera[i] + camera_path[i]);
        let camera_state = self.visibility.visibility(self.scene.primitives, ray.camera, camera_path, point(1.0), cx)?;
        if camera_state == Visibility::Blocked { return Ok(()); }
        let source = &self.scene.sources[panel.source];
        let origin = coordinates(source.position_m());
        let towards = difference(origin, position);
        let radius = norm(towards).ok_or(ScatterEnclosureError::NumericalRange)?;
        let lower = radius.lo().max(panel.minimum_distance);
        if lower > radius.hi() { return Err(ScatterEnclosureError::NumericalRange); }
        let radius = Interval::new(lower, radius.hi());
        let towards_unit = towards.map(|v| v / radius);
        // A wide panel's direction box may contain zero even though the real
        // line stays separated. Preserve the source's global intensity ceiling
        // until refinement resolves that dependency, rather than dropping light.
        let intensity = intensity_bound(source, towards_unit.map(|v| -v), cx)?
            .unwrap_or_else(|| Interval::new(0.0,
                (point(source.flux_lm()) * point(source.distribution.peak_per_sr)).hi()));
        if !finite(intensity) { return Err(ScatterEnclosureError::NumericalRange); }
        if intensity.hi() == 0.0 { return Ok(()); }
        let end = point(1.0) - point(self.scene.options.ray_epsilon_m) / radius;
        let end = Interval::new(end.lo().clamp(0.0, 1.0), end.hi().clamp(0.0, 1.0));
        let light_state = self.visibility.visibility(self.scene.primitives, origin, towards.map(|v| -v), end, cx)?;
        if light_state == Visibility::Blocked { return Ok(()); }
        let cosine = dot(towards_unit, ray.delta.map(|v| v / ray.length));
        if !finite(cosine) { return Err(ScatterEnclosureError::NumericalRange); }
        let mu = Interval::new(cosine.lo().clamp(-1.0, 1.0), cosine.hi().clamp(-1.0, 1.0));
        let phase = phase::enclose(self.phase, mu, cx)?;
        // The sum of both leg lengths is convex along the camera segment:
        // its maximum is at a panel endpoint. Triangle inequality supplies an
        // independent lower bound, including the exact collinear path limit.
        let offset = difference(origin, ray.camera);
        let straight = norm(offset).ok_or(ScatterEnclosureError::NumericalRange)?;
        let mut upper = 0.0_f64;
        for endpoint in [panel.lo, panel.hi] {
            let radius = norm(difference(offset, ray.delta.map(|v| v * point(endpoint))))
                .ok_or(ScatterEnclosureError::NumericalRange)?;
            upper = upper.max(add(radius, mul(point(endpoint), ray.length)).hi());
        }
        let path = add(radius, mul(t, ray.length));
        let lower = path.lo().max(straight.lo());
        let upper = path.hi().min(upper);
        if lower > upper { return Err(ScatterEnclosureError::NumericalRange); }
        let transmission = exp_neg(mul(self.extinction, Interval::new(lower, upper)))
            .ok_or(ScatterEnclosureError::NumericalRange)?;
        let value = mul(mul(mul(nonnegative(intensity / radius / radius), transmission), phase),
            mul(mul(point(self.scattering), ray.length), nonnegative(point(panel.hi)-point(panel.lo))));
        if !finite(value) { return Err(ScatterEnclosureError::NumericalRange); }
        panel.uncertain = value.hi() > 0.0 && (camera_state == Visibility::Uncertain || light_state == Visibility::Uncertain);
        panel.range = if panel.uncertain { Interval::new(0.0, value.hi()) } else { value };
        Ok(())
    }
}

struct SegmentGeometry { camera: Vector, delta: Vector, length: Interval }
#[derive(Clone, Copy)]
struct Panel {
    source: usize, lo: f64, hi: f64, minimum_distance: f64,
    depth: u32, order: usize, range: Interval, uncertain: bool,
}
impl PartialEq for Panel { fn eq(&self, other: &Self) -> bool { self.cmp(other) == Ordering::Equal } }
impl Eq for Panel {}
impl PartialOrd for Panel { fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) } }
impl Ord for Panel {
    fn cmp(&self, other: &Self) -> Ordering {
        self.range.width().total_cmp(&other.range.width()).then_with(|| other.order.cmp(&self.order))
    }
}
struct Workspace { frontier: BinaryHeap<Panel> }
impl Workspace {
    fn new(budget: ScatterEnclosureBudget, sources: usize) -> Result<Self, ScatterEnclosureError> {
        if !budget.absolute_error_cd_m2.is_finite() || budget.absolute_error_cd_m2 < 0.0
            || budget.max_depth > 48 || budget.max_evaluations == 0
            || budget.max_evaluations > MAX_SCATTER_ENCLOSURE_PANELS || sources > budget.max_evaluations {
            return Err(ScatterEnclosureError::Budget);
        }
        let mut frontier = BinaryHeap::new();
        // Each split spends two evaluations but increases the leaf count by one.
        let capacity = sources + (budget.max_evaluations - sources) / 2;
        frontier.try_reserve_exact(capacity).map_err(|_| ScatterEnclosureError::Allocation)?;
        Ok(Self { frontier })
    }
    fn retain(&mut self, panel: Panel, budget: ScatterEnclosureBudget) {
        if panel.depth < budget.max_depth && panel.range.hi() > panel.range.lo() { self.frontier.push(panel); }
    }
}

/// Progressive original camera segments. Reference preparation and bounded heap
/// are reused. A failed/cancelled segment is unpublished and can be retried;
/// completed segments remain bit-identical under different advance partitions.
pub struct ScatterEnclosureJob<'a, 'scene, 'field> {
    reference: &'a HomogeneousScatterReference<'scene, 'field>,
    segments: &'a [ScatterSegment], budget: ScatterEnclosureBudget,
    workspace: Workspace, values: Vec<EnclosedSingleScatter>,
}
impl<'a, 'scene, 'field> ScatterEnclosureJob<'a, 'scene, 'field> {
    /// Borrow finite camera segments without copying geometry or phase data.
    pub fn new(reference: &'a HomogeneousScatterReference<'scene, 'field>, segments: &'a [ScatterSegment],
        budget: ScatterEnclosureBudget, cx: &Cx<'_>) -> Result<Self, ScatterEnclosureError> {
        cx.checkpoint()?;
        if segments.len() > reference.scene.options.max_outputs { return Err(ScatterEnclosureError::Budget); }
        let workspace = Workspace::new(budget, reference.scene.sources.len())?;
        let mut values = Vec::new();
        values.try_reserve_exact(segments.len()).map_err(|_| ScatterEnclosureError::Allocation)?;
        cx.checkpoint()?;
        Ok(Self { reference, segments, budget, workspace, values })
    }
    /// Append at most `max_segments` whole component results; no partial ray leaks.
    pub fn advance(&mut self, cx: &Cx<'_>, max_segments: usize) -> Result<usize, ScatterEnclosureError> {
        cx.checkpoint()?;
        let start = self.values.len();
        let end = start.saturating_add(max_segments).min(self.segments.len());
        for index in start..end {
            let result = self.reference.evaluate(cx, self.segments[index], self.budget, &mut self.workspace)?;
            cx.checkpoint()?;
            self.values.push(result);
        }
        Ok(self.values.len()-start)
    }
    /// Completed prefix only; absent rays are not zero radiance.
    #[must_use]
    pub fn results(&self) -> &[EnclosedSingleScatter] { &self.values }
    /// Every declared segment has a completed enclosure.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.values.len() == self.segments.len() }
}

fn point(v: f64) -> Interval { Interval::point(v) }
fn zero() -> Interval { point(0.0) }
fn coordinates(p: Point3) -> Vector { [p.x, p.y, p.z].map(point) }
fn difference(a: Vector, b: Vector) -> Vector { std::array::from_fn(|i| a[i]-b[i]) }
fn add(a: Interval, b: Interval) -> Interval {
    if a.hi() == 0.0 { b } else if b.hi() == 0.0 { a } else { nonnegative(a+b) }
}
fn mul(a: Interval, b: Interval) -> Interval {
    if a.hi() == 0.0 || b.hi() == 0.0 { zero() } else { nonnegative(a*b) }
}
fn midpoint_error(range: Interval) -> (f64, f64) {
    let mid = range.lo() + (range.hi()-range.lo())*0.5;
    let error = if range.lo() == range.hi() { 0.0 } else {
        (point(mid)-point(range.lo())).hi().max((point(range.hi())-point(mid)).hi())
    };
    (mid, error)
}

#[cfg(test)]
mod tests;
