//! Direct finite-area illuminance with bounded opaque-mesh shadows in vacuum.
//!
//! Reuse the native emitter preparation and analytical projected-triangle
//! integral. The original interval visibility hierarchy classifies ALL rays
//! from a source subtriangle, not a few sampled shadow rays. Clear patches use
//! the analytical integral, blocked patches are zero, and unresolved patches
//! retain [0, unoccluded upper]. Bisect the largest uncertain lux contribution
//! first. A work/depth limit returns the remaining enclosure, never a fabricated
//! tolerance pass. No source, primitive, material or path estimator is modified.
//!
//! This encloses selected-source DIRECT emission only. Opaque static meshes and
//! rigid instances may occlude one another, including faces of folded sources.
//! Reflection, refraction, haze, projector optics, calibration and physical mesh
//! fidelity remain outside the model. As in the existing direct connector, the
//! source endpoint and a declared receiver-end distance are excluded. The
//! physical effect of that exclusion and total physical error remain unknown.

use super::*;
use crate::photometry::PhotometryError;
use crate::photometry::source::validation::{EnclosureUnavailable, visibility::{Visibility, VisibilityIndex}};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Version of adaptive patch partitioning, shadow predicates and output scope.
pub const AREA_OCCLUSION_SEMANTICS_VERSION: u32 = 1;
/// Hard ceiling on patch evaluations and reserved frontier slots per receiver.
pub const MAX_AREA_SHADOW_PATCHES: usize = 262_144;

/// Per-receiver accuracy/work limits. Exhaustion is reported in the result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AreaOcclusionBudget {
    /// Requested maximum midpoint numerical error, in lux; finite and >= 0.
    pub absolute_error_lux: f64,
    /// Includes initial source triangles and both children of every bisection.
    /// Must cover all prepared source triangles and not exceed the hard ceiling.
    pub max_patch_evaluations: usize,
    /// Maximum binary subdivision depth, in 0..=32; zero retains root bounds.
    pub max_depth: u32,
    /// Positive finite distance excluded at the RECEIVER only, in metres.
    pub ray_epsilon_m: f64,
}
impl Default for AreaOcclusionBudget {
    fn default() -> Self {
        Self { absolute_error_lux: 0.1, max_patch_evaluations: 16_384, max_depth: 24, ray_epsilon_m: 1e-6 }
    }
}

/// Refusal before a complete receiver enclosure is available.
#[derive(Debug)]
pub enum AreaOcclusionError {
    /// Cancellation preserves all previously published receivers.
    Cancelled,
    /// Invalid accuracy, subdivision or endpoint settings.
    InvalidOptions,
    /// Initial triangles cannot fit the work limit, or output storage is capped.
    Budget,
    /// Bounded frontier/output reservation failed.
    Allocation,
    /// Original area-source admission or analytical evaluation refusal.
    Reference(AreaReferenceError),
    /// Unsupported or over-budget native visibility representation.
    GeometryUnavailable(EnclosureUnavailable),
    /// Transparent geometry cannot be treated as an opaque shadow caster.
    TransmissiveGeometry,
    /// Original visibility preparation/evaluation failure.
    Transport(PhotometryError),
}
impl core::fmt::Display for AreaOcclusionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "occluded area reference refused: {self:?}")
    }
}
impl std::error::Error for AreaOcclusionError {}
impl From<fs_exec::Cancelled> for AreaOcclusionError {
    fn from(_: fs_exec::Cancelled) -> Self { Self::Cancelled }
}
impl From<AreaReferenceError> for AreaOcclusionError {
    fn from(error: AreaReferenceError) -> Self {
        match error { AreaReferenceError::Cancelled => Self::Cancelled, other => Self::Reference(other) }
    }
}
impl From<PhotometryError> for AreaOcclusionError {
    fn from(error: PhotometryError) -> Self {
        match error { PhotometryError::Cancelled => Self::Cancelled, other => Self::Transport(other) }
    }
}

/// A tolerance pass is distinct from a bounded but insufficiently tight result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AreaOcclusionTermination {
    /// The returned midpoint's outward numerical error meets the requested lux.
    ToleranceMet,
    /// No more patch evaluations were available; the current bound is retained.
    PatchBudget,
    /// No refinable patches remain (depth limit or irreducible arithmetic width).
    ResolutionLimit,
}

/// Complete direct-area result including unresolved shadow uncertainty.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OccludedAreaIlluminance {
    /// Deterministic interval midpoint, not a sampled/path-tracer estimate.
    pub lux: f64,
    /// Selected-source direct vacuum illuminance with the declared mesh shadows.
    pub numerical: NumericalCertificate,
    /// Outward midpoint error, INCLUDING remaining penumbra ambiguity, in lux.
    pub numerical_error_bound_lux: f64,
    /// Requested numerical accuracy, kept separate from achieved accuracy.
    pub requested_error_lux: f64,
    /// Why adaptive evaluation stopped; every variant retains a valid bound.
    pub termination: AreaOcclusionTermination,
    /// Actual patch evaluations, including initial triangles and discarded parents.
    pub evaluated_patches: usize,
    /// Final leaf patches whose possibly nonzero contribution has unknown visibility.
    pub uncertain_visibility_patches: usize,
    /// Final leaf patches using a geometric fallback for the angular integral.
    pub coarse_integral_patches: usize,
    /// Selected-source model, original occluders, probe, and endpoint convention.
    /// This is not authenticated calibration evidence.
    pub numerical_model_identity: ContentHash,
    /// Unknown calibration, physical geometry and omitted-transport error.
    pub total_error_bound_lux: Option<f64>,
}
impl OccludedAreaIlluminance {
    /// CSV labels its direct/vacuum scope and actual numerical termination.
    pub fn write_csv_header(writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(b"index,scope,lux,numerical_lower_lux,numerical_upper_lux,numerical_error_bound_lux,requested_error_lux,termination,evaluated_patches,uncertain_visibility_patches,coarse_integral_patches,numerical_model_identity,total_physical_error_bound_lux,version\n")
    }
    /// Missing total error remains unknown; I/O failures propagate.
    pub fn write_csv_row(self, index: usize, writer: &mut impl Write) -> io::Result<()> {
        let reason = match self.termination {
            AreaOcclusionTermination::ToleranceMet => "tolerance-met",
            AreaOcclusionTermination::PatchBudget => "patch-budget",
            AreaOcclusionTermination::ResolutionLimit => "resolution-limit",
        };
        writeln!(writer, "{index},opaque-mesh-vacuum-direct-area,{},{},{},{},{},{reason},{},{},{},{},unknown,{AREA_OCCLUSION_SEMANTICS_VERSION}",
            self.lux, self.numerical.lo, self.numerical.hi, self.numerical_error_bound_lux,
            self.requested_error_lux, self.evaluated_patches, self.uncertain_visibility_patches,
            self.coarse_integral_patches, self.numerical_model_identity)
    }
}

/// Shared original-source preparation plus the existing interval visibility index.
/// Source meshes are enclosed once; all occluder geometry remains borrowed.
pub struct OccludedAreaReference<'scene> {
    sources: UnoccludedAreaReference,
    primitives: &'scene [Primitive],
    visibility: VisibilityIndex,
    identity: ContentHash,
}
impl<'scene> OccludedAreaReference<'scene> {
    /// Admit uniform native mesh emitters and the same scene's opaque blockers.
    /// All scene geometry participates, including non-emitting and back-facing
    /// source faces. Dielectrics, animated/chart shapes and projectors refuse.
    pub fn new(
        primitives: &'scene [Primitive], sources: &[PhotometricAreaEmitter],
        preparation: AreaReferenceBudget, cx: &Cx<'_>,
    ) -> Result<Self, AreaOcclusionError> {
        cx.checkpoint()?;
        for primitive in primitives {
            cx.checkpoint()?;
            if matches!(primitive.material, Material::Dielectric { .. }) {
                return Err(AreaOcclusionError::TransmissiveGeometry);
            }
        }
        let sources = UnoccludedAreaReference::new(primitives, sources, preparation, cx)?;
        let mut hash = DomainHasher::new("org.frankensim.photometry.occluded-area.v1");
        hash.update(&AREA_OCCLUSION_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&sources.identity.0);
        let visibility = VisibilityIndex::prepare(primitives, &mut hash, cx)?
            .map_err(AreaOcclusionError::GeometryUnavailable)?;
        cx.checkpoint()?;
        Ok(Self { sources, primitives, visibility, identity: hash.finalize() })
    }

    /// Bound one receiver. Progressive jobs reuse the bounded frontier allocation.
    pub fn illuminance(
        &self, cx: &Cx<'_>, probe: &IlluminanceProbe, budget: AreaOcclusionBudget,
    ) -> Result<OccludedAreaIlluminance, AreaOcclusionError> {
        cx.checkpoint()?;
        let mut workspace = Workspace::new(budget, self.sources.faces.len())?;
        self.evaluate(cx, probe, budget, &mut workspace)
    }

    fn evaluate(
        &self, cx: &Cx<'_>, probe: &IlluminanceProbe, budget: AreaOcclusionBudget, work: &mut Workspace,
    ) -> Result<OccludedAreaIlluminance, AreaOcclusionError> {
        cx.checkpoint()?;
        work.frontier.clear();
        let p = probe.position_m();
        let receiver = [p.x, p.y, p.z].map(point);
        let n = probe.normal();
        let normal = normalized([n.x, n.y, n.z].map(point)).ok_or(AreaReferenceError::NumericalRange)?;
        let mut raw = zero();
        let mut evaluations = 0;
        let mut uncertain = 0;
        let mut coarse = 0;
        for (root, face) in self.sources.faces.iter().enumerate() {
            cx.checkpoint()?;
            let mut patch = Patch { root, vertices: face.vertices, area: face.area,
                depth: 0, range: zero(), uncertain: false, coarse: false, order: evaluations };
            if face.luminance.hi() > 0.0 {
                // Establish source-plane separation ONCE for the original face.
                // Its exact subdivisions stay in that plane, including the
                // source-face exclusion used by the shared visibility owner.
                let first = normalized(difference(face.vertices[0], receiver))
                    .ok_or(AreaReferenceError::UnresolvedSourcePlane { triangle: root })?;
                let side = dot(face.perpendicular, first);
                if !finite(side) { return Err(AreaReferenceError::NumericalRange.into()); }
                if side.contains_zero() { return Err(AreaReferenceError::UnresolvedSourcePlane { triangle: root }.into()); }
                if face.both_sides || side.hi() < 0.0 {
                    self.evaluate_patch(&mut patch, receiver, normal, budget, cx)?;
                }
            }
            evaluations += 1;
            raw = nonnegative(add(raw, patch.range));
            if !finite(raw) { return Err(AreaReferenceError::NumericalRange.into()); }
            uncertain += usize::from(patch.uncertain);
            coarse += usize::from(patch.coarse);
            work.retain(patch, budget);
        }
        let mut enclosed = raw;
        let termination = loop {
            cx.checkpoint()?;
            if midpoint_and_error(enclosed).1 <= budget.absolute_error_lux {
                break AreaOcclusionTermination::ToleranceMet;
            }
            if work.frontier.is_empty() { break AreaOcclusionTermination::ResolutionLimit; }
            if evaluations + 2 > budget.max_patch_evaluations { break AreaOcclusionTermination::PatchBudget; }
            let parent = work.frontier.pop().ok_or(AreaReferenceError::NumericalRange)?;
            let mut children = parent.split(evaluations);
            for child in &mut children {
                self.evaluate_patch(child, receiver, normal, budget, cx)?;
            }
            evaluations += 2;
            // Raw endpoints bound the SUM OF LEAF ENDPOINTS, not just the true
            // total. This invariant makes parent replacement sound. Never use
            // the separately intersected enclosure as this accumulator.
            raw = replace_parent(raw, parent.range, children[0].range, children[1].range)?;
            enclosed = enclosed.intersect(raw).ok_or(AreaReferenceError::NumericalRange)?;
            uncertain -= usize::from(parent.uncertain);
            coarse -= usize::from(parent.coarse);
            for child in children {
                uncertain += usize::from(child.uncertain);
                coarse += usize::from(child.coarse);
                work.retain(child, budget);
            }
        };
        let (lux, error) = midpoint_and_error(enclosed);
        if !error.is_finite() { return Err(AreaReferenceError::NumericalRange.into()); }
        let mut hash = DomainHasher::new("org.frankensim.photometry.occluded-area.probe.v1");
        hash.update(&self.identity.0);
        for value in [p.x, p.y, p.z, n.x, n.y, n.z, budget.ray_epsilon_m] { hash_float(&mut hash, value); }
        cx.checkpoint()?;
        Ok(OccludedAreaIlluminance { lux, numerical: NumericalCertificate::enclosure(enclosed.lo(), enclosed.hi()),
            numerical_error_bound_lux: error, requested_error_lux: budget.absolute_error_lux,
            termination, evaluated_patches: evaluations, uncertain_visibility_patches: uncertain,
            coarse_integral_patches: coarse, numerical_model_identity: hash.finalize(), total_error_bound_lux: None })
    }

    fn evaluate_patch(
        &self, patch: &mut Patch, receiver: Vector, normal: Vector,
        budget: AreaOcclusionBudget, cx: &Cx<'_>,
    ) -> Result<(), AreaOcclusionError> {
        cx.checkpoint()?;
        let face = &self.sources.faces[patch.root];
        let vectors = patch.vertices.map(|v| difference(v, receiver));
        let (factor, coarse) = projected_triangle(vectors, normal, patch.area);
        if factor.hi() == 0.0 { return Ok(()); }
        let origin = std::array::from_fn(|i| patch.vertices[0][i].hull(patch.vertices[1][i]).hull(patch.vertices[2][i]));
        let direction = difference(receiver, origin);
        let distance = norm(direction).ok_or(AreaReferenceError::NumericalRange)?;
        let end = if distance.lo() > 0.0 {
            let range = point(1.0) - point(budget.ray_epsilon_m) / distance;
            Interval::new(range.lo().clamp(0.0, 1.0), range.hi().clamp(0.0, 1.0))
        } else {
            // The box may contain the receiver even when no point of the
            // actual triangle does. Keep every possible endpoint convention.
            Interval::new(0.0, 1.0)
        };
        let state = self.visibility.visibility_from_face(self.primitives, patch.vertices, receiver, end, face.origin_face, cx)?;
        if state == Visibility::Blocked { return Ok(()); }
        let value = nonnegative(mul(face.luminance, factor));
        if !finite(value) { return Err(AreaReferenceError::NumericalRange.into()); }
        patch.uncertain = state == Visibility::Uncertain && value.hi() > 0.0;
        patch.coarse = coarse;
        patch.range = if patch.uncertain { Interval::new(0.0, value.hi()) } else { value };
        Ok(())
    }
}

fn midpoint_and_error(range: Interval) -> (f64, f64) {
    let mid = range.lo() + (range.hi() - range.lo()) * 0.5;
    (mid, sub(point(mid), point(range.lo())).hi().max(sub(point(range.hi()), point(mid)).hi()))
}

fn replace_parent(total: Interval, parent: Interval, a: Interval, b: Interval) -> Result<Interval, AreaReferenceError> {
    let lower = add(add(sub(point(total.lo()), point(parent.lo())), point(a.lo())), point(b.lo())).lo().max(0.0);
    let upper = add(add(sub(point(total.hi()), point(parent.hi())), point(a.hi())), point(b.hi())).hi();
    if !lower.is_finite() || !upper.is_finite() || lower > upper { return Err(AreaReferenceError::NumericalRange); }
    Ok(Interval::new(lower, upper))
}

#[derive(Clone, Copy)]
struct Patch {
    root: usize,
    vertices: [Vector; 3],
    area: Interval,
    range: Interval,
    depth: u32,
    uncertain: bool,
    coarse: bool,
    order: usize,
}
impl Patch {
    fn split(self, order: usize) -> [Self; 2] {
        // Longest bounding edge is a heuristic only; any edge bisection exactly
        // partitions the original triangle and halves its geometric area.
        let edge = (0..3).max_by(|&a, &b| {
            let length = |i: usize| norm(difference(self.vertices[(i+1)%3], self.vertices[i]))
                .map_or(f64::INFINITY, |v| v.hi());
            length(a).total_cmp(&length(b)).then_with(|| b.cmp(&a))
        }).unwrap_or(0);
        let a = edge; let b = (edge+1)%3; let c = (edge+2)%3;
        let midpoint = std::array::from_fn(|i| add(mul(point(0.5), self.vertices[a][i]), mul(point(0.5), self.vertices[b][i])));
        let triangles = [[self.vertices[a], midpoint, self.vertices[c]], [midpoint, self.vertices[b], self.vertices[c]]];
        std::array::from_fn(|i| Self { vertices: triangles[i], area: mul(self.area, point(0.5)),
            depth: self.depth+1, range: zero(), uncertain: false, coarse: false, order: order+i, ..self })
    }
}
impl PartialEq for Patch { fn eq(&self, other: &Self) -> bool { self.cmp(other) == Ordering::Equal } }
impl Eq for Patch {}
impl PartialOrd for Patch { fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) } }
impl Ord for Patch {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.range.hi()-self.range.lo()).total_cmp(&(other.range.hi()-other.range.lo()))
            .then_with(|| other.order.cmp(&self.order))
    }
}

struct Workspace { frontier: BinaryHeap<Patch> }
impl Workspace {
    fn new(budget: AreaOcclusionBudget, roots: usize) -> Result<Self, AreaOcclusionError> {
        if !budget.absolute_error_lux.is_finite() || budget.absolute_error_lux < 0.0
            || budget.max_depth > 32 || !budget.ray_epsilon_m.is_finite() || budget.ray_epsilon_m <= 0.0
            || budget.max_patch_evaluations == 0 || budget.max_patch_evaluations > MAX_AREA_SHADOW_PATCHES {
            return Err(AreaOcclusionError::InvalidOptions);
        }
        if roots > budget.max_patch_evaluations { return Err(AreaOcclusionError::Budget); }
        let mut frontier = BinaryHeap::new();
        // Frontier leaves never outnumber evaluated patches. All later pushes
        // fit this reservation, even if cancellation interrupts a split.
        frontier.try_reserve_exact(budget.max_patch_evaluations).map_err(|_| AreaOcclusionError::Allocation)?;
        Ok(Self { frontier })
    }
    fn retain(&mut self, patch: Patch, budget: AreaOcclusionBudget) {
        if (patch.uncertain || patch.coarse) && patch.depth < budget.max_depth && patch.range.hi() > patch.range.lo() {
            self.frontier.push(patch);
        }
    }
}

/// Progressive receivers sharing immutable geometry and one bounded frontier.
/// Refinement is internal to each receiver; cancellation discards that receiver
/// and retains the published prefix. Restarting recomputes only its unpublished
/// receiver. Chunk sizes do not change bounds, ordering or identities.
pub struct OccludedAreaJob<'a, 'scene> {
    reference: &'a OccludedAreaReference<'scene>,
    receivers: Receivers<'a>,
    count: usize,
    budget: AreaOcclusionBudget,
    work: Workspace,
    values: Vec<OccludedAreaIlluminance>,
}
impl<'a, 'scene> OccludedAreaJob<'a, 'scene> {
    /// Borrow original probe positions and normals.
    pub fn new(reference: &'a OccludedAreaReference<'scene>, probes: &'a [IlluminanceProbe], budget: AreaOcclusionBudget, cx: &Cx<'_>) -> Result<Self, AreaOcclusionError> {
        Self::begin(reference, Receivers::Probes(probes), probes.len(), budget, cx)
    }
    /// Use the original point-sampled calculation plane, not cell-area averages.
    pub fn for_grid(reference: &'a OccludedAreaReference<'scene>, grid: CalculationGrid, budget: AreaOcclusionBudget, cx: &Cx<'_>) -> Result<Self, AreaOcclusionError> {
        Self::begin(reference, Receivers::Grid(grid), grid.len(), budget, cx)
    }
    fn begin(reference: &'a OccludedAreaReference<'scene>, receivers: Receivers<'a>, count: usize, budget: AreaOcclusionBudget, cx: &Cx<'_>) -> Result<Self, AreaOcclusionError> {
        cx.checkpoint()?;
        if count > reference.sources.max_outputs { return Err(AreaOcclusionError::Budget); }
        let work = Workspace::new(budget, reference.sources.faces.len())?;
        let mut values = Vec::new();
        values.try_reserve_exact(count).map_err(|_| AreaOcclusionError::Allocation)?;
        cx.checkpoint()?;
        Ok(Self { reference, receivers, count, budget, work, values })
    }
    /// Append complete numerical enclosures, including explicitly budget-limited
    /// ones. No failed or cancelled receiver is published as zero or a tolerance pass.
    pub fn advance(&mut self, cx: &Cx<'_>, max_points: usize) -> Result<usize, AreaOcclusionError> {
        cx.checkpoint()?;
        let start = self.values.len();
        let end = start.saturating_add(max_points).min(self.count);
        for index in start..end {
            let probe = match &self.receivers {
                Receivers::Probes(probes) => probes[index],
                Receivers::Grid(grid) => grid.probe(index).map_err(|_| AreaReferenceError::NumericalRange)?,
            };
            let value = self.reference.evaluate(cx, &probe, self.budget, &mut self.work)?;
            cx.checkpoint()?;
            self.values.push(value);
        }
        Ok(self.values.len()-start)
    }
    /// Completed receiver prefix only.
    #[must_use]
    pub fn results(&self) -> &[OccludedAreaIlluminance] { &self.values }
    /// Every receiver has a result; inspect termination for requested accuracy.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.values.len() == self.count }
}

#[cfg(test)]
mod tests;
