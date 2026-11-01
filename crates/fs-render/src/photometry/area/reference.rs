//! Deterministic, unoccluded vacuum references for native uniform mesh emitters.
//!
//! Integrate the receiving cosine over each triangle's spherical projection by
//! a boundary integral, after clipping to the receiving hemisphere. There is no
//! point-source replacement, random sampling, or numerical area quadrature.
//! Geometry, rigid placements, spectra, flux and emitting sides come from the
//! original `Primitive` / `PhotometricAreaEmitter` owners. Outward arithmetic
//! includes geometric area normalization and the retained spectral/observer sum.
//!
//! This is an explicitly unoccluded, zero-bounce, vacuum reference, not a scene
//! renderer: other primitives, mutual source occlusion, reflection, refraction,
//! haze and projector optics are not evaluated. Use it to validate the original
//! path estimator on scenes satisfying those assumptions. A bound here never
//! strengthens an unrelated path estimate or establishes total physical accuracy.
//! Source-construction history, physical fixture/geometry calibration, spectral
//! discretization and wavelengths outside the retained band remain unbounded.

use super::{EmissionSidedness, PhotometricAreaEmitter};
use crate::charts::TriMesh;
use crate::instances::{RigidTransform, SharedGeometry};
use crate::photometry::{CalculationGrid, IlluminanceProbe, Receivers};
use crate::photometry::observer;
use crate::photometry::source::{BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W};
use crate::tracer::{Material, Primitive, Shape};
use fs_blake3::{ContentHash, DomainHasher};
use fs_evidence::NumericalCertificate;
use fs_exec::Cx;
use fs_ivl::{Interval, angle::atan2_upper_half_plane};
use std::io::{self, Write};

pub mod occlusion;

/// Version of the reference's geometry, clipping and numerical interpretation.
pub const AREA_REFERENCE_SEMANTICS_VERSION: u32 = 1;

/// Explicit preparation storage and progressive output limits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AreaReferenceBudget {
    /// Maximum supplied emitter bindings, including zero-flux bindings.
    pub max_emitters: usize,
    /// Maximum total bound triangles; preparation also uses one emitter's keys.
    pub max_triangles: usize,
    /// Maximum retained receiver results in one progressive job.
    pub max_outputs: usize,
}

impl Default for AreaReferenceBudget {
    fn default() -> Self {
        Self { max_emitters: 4096, max_triangles: 65_536, max_outputs: 1_048_576 }
    }
}

/// Explicit reference refusal; a failed receiver is never published as zero.
#[derive(Debug)]
pub enum AreaReferenceError {
    /// Preparation or evaluation was cancelled.
    Cancelled,
    /// A declared emitter, triangle or output limit was exceeded.
    Budget,
    /// Bounded workspace allocation failed.
    Allocation,
    /// Missing, degenerate or repeated triangles, or duplicate source bindings.
    InvalidGeometry,
    /// Only static opaque mesh/rigid-instance, uniform Lambertian sources qualify.
    UnsupportedSource,
    /// A required finite positive geometric/radiometric quantity was unresolved.
    NumericalRange,
    /// Receiver separation from this source triangle's plane was not established.
    UnresolvedSourcePlane {
        /// Triangle in the prepared, source-order traversal.
        triangle: usize,
    },
}

impl core::fmt::Display for AreaReferenceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "unoccluded area reference refused: {self:?}")
    }
}
impl std::error::Error for AreaReferenceError {}
impl From<fs_exec::Cancelled> for AreaReferenceError {
    fn from(_: fs_exec::Cancelled) -> Self { Self::Cancelled }
}

type Vector = [Interval; 3];

struct Face {
    // Original face identity for source-endpoint exclusion in occlusion queries.
    origin_face: (usize, usize),
    vertices: [Vector; 3],
    // Positively rescaled winding normal, used only for source-side admission.
    perpendicular: Vector,
    area: Interval,
    luminance: Interval,
    both_sides: bool,
}

/// Immutable prepared numerical view of the original selected emitters.
///
/// World-coordinate enclosures and source normalizations are prepared once.
/// Receiver queries use O(triangles) work and fixed local storage, with no
/// per-receiver allocation. Only the selected source geometry is retained; the
/// identity is not a complete scene hash and does not imply unobstructed rays.
pub struct UnoccludedAreaReference {
    faces: Vec<Face>,
    identity: ContentHash,
    max_outputs: usize,
}

impl UnoccludedAreaReference {
    /// Prepare original native meshes and source bindings under explicit limits.
    /// No source/material object is modified. Duplicate/coincident faces within
    /// a source refuse as in native area admission. Projectors refuse explicitly.
    pub fn new(
        primitives: &[Primitive], sources: &[PhotometricAreaEmitter],
        budget: AreaReferenceBudget, cx: &Cx<'_>,
    ) -> Result<Self, AreaReferenceError> {
        cx.checkpoint()?;
        if sources.len() > budget.max_emitters { return Err(AreaReferenceError::Budget); }
        let mut total = 0_usize;
        let mut indices = Vec::new();
        indices.try_reserve_exact(sources.len()).map_err(|_| AreaReferenceError::Allocation)?;
        for source in sources {
            cx.checkpoint()?;
            let (mesh, _) = mesh(primitives, source)?;
            if mesh.triangles.is_empty() { return Err(AreaReferenceError::InvalidGeometry); }
            total = total.checked_add(mesh.triangles.len()).ok_or(AreaReferenceError::Budget)?;
            if total > budget.max_triangles { return Err(AreaReferenceError::Budget); }
            indices.push(source.primitive_index());
        }
        indices.sort_unstable();
        if indices.windows(2).any(|w| w[0] == w[1]) { return Err(AreaReferenceError::InvalidGeometry); }
        let mut faces = Vec::new();
        faces.try_reserve_exact(total).map_err(|_| AreaReferenceError::Allocation)?;
        let mut hash = DomainHasher::new("org.frankensim.photometry.unoccluded-area.v1");
        hash.update(&AREA_REFERENCE_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&fs_ivl::INTERVAL_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&fs_ivl::angle::ANGLE_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&(sources.len() as u64).to_le_bytes());
        for source in sources {
            cx.checkpoint()?;
            let (mesh, placement) = mesh(primitives, source)?;
            let first = faces.len();
            let mut area = zero();
            let mut keys = Vec::new();
            keys.try_reserve_exact(mesh.triangles.len()).map_err(|_| AreaReferenceError::Allocation)?;
            hash.update(&(source.primitive_index() as u64).to_le_bytes());
            hash.update(&(mesh.triangles.len() as u64).to_le_bytes());
            let both = source.sidedness() == EmissionSidedness::Both;
            hash.update(&[u8::from(both)]);
            hash_float(&mut hash, source.flux_lm());
            for (triangle_index, triangle) in mesh.triangles.iter().enumerate() {
                cx.checkpoint()?;
                let mut vertices = [[zero(); 3]; 3];
                let mut key = [[0_u64; 3]; 3];
                for corner in 0..3 {
                    let point = *mesh.vertices.get(triangle[corner] as usize)
                        .ok_or(AreaReferenceError::InvalidGeometry)?;
                    if point.iter().any(|v| !v.is_finite()) { return Err(AreaReferenceError::InvalidGeometry); }
                    // Proper rigid placement is injective: coincident local
                    // triangles are coincident world triangles, regardless of winding.
                    key[corner] = point.map(|v| if v == 0.0 { 0 } else { v.to_bits() });
                    vertices[corner] = transform(point, placement);
                    for coordinate in vertices[corner] { hash_range(&mut hash, coordinate); }
                }
                key.sort_unstable();
                keys.push(key);
                let edge_a = difference(vertices[1], vertices[0]);
                let edge_b = difference(vertices[2], vertices[0]);
                let scale = magnitude_scale(&[edge_a, edge_b])
                    .ok_or(AreaReferenceError::NumericalRange)?;
                if scale == 0.0 { return Err(AreaReferenceError::InvalidGeometry); }
                let a = edge_a.map(|x| x / point(scale));
                let b = edge_b.map(|x| x / point(scale));
                let perpendicular = cross(a, b);
                let length = norm(perpendicular).ok_or(AreaReferenceError::NumericalRange)?;
                if length.lo() <= 0.0 { return Err(AreaReferenceError::InvalidGeometry); }
                let face_area = mul(mul(mul(point(0.5), length), point(scale)), point(scale));
                if !finite(face_area) || face_area.lo() <= 0.0 { return Err(AreaReferenceError::NumericalRange); }
                area = add(area, face_area);
                faces.push(Face { origin_face: (source.primitive_index(), triangle_index), vertices, perpendicular, area: face_area, luminance: zero(), both_sides: both });
            }
            keys.sort_unstable();
            if keys.windows(2).any(|w| w[0] == w[1]) { return Err(AreaReferenceError::InvalidGeometry); }
            if !finite(area) { return Err(AreaReferenceError::NumericalRange); }
            // Retain the actual spectral normalization; do not assume its
            // binary64 round trip is exactly one photopic lumen.
            let mut photopic = zero();
            for (bin, &power) in source.spectrum().bins_w_per_lumen_nm().iter().enumerate() {
                cx.checkpoint()?;
                let observer_y = observer::xyz_at_bin(bin)[1];
                hash_float(&mut hash, power);
                hash_float(&mut hash, observer_y);
                photopic = add(photopic, mul(point(power), point(observer_y)));
            }
            photopic = mul(mul(photopic, point(BIN_WIDTH_NM)), point(LUMINOUS_EFFICACY_LM_W));
            let luminance = if source.flux_lm() == 0.0 { zero() } else {
                nonnegative(mul(point(source.flux_lm()) / area / pi()
                    / point(if both { 2.0 } else { 1.0 }), photopic))
            };
            if !finite(luminance) { return Err(AreaReferenceError::NumericalRange); }
            hash_range(&mut hash, luminance);
            for face in &mut faces[first..] { cx.checkpoint()?; face.luminance = luminance; }
        }
        cx.checkpoint()?;
        Ok(Self { faces, identity: hash.finalize(), max_outputs: budget.max_outputs })
    }

    /// Enclose selected-source direct illuminance under the named assumptions.
    /// A difficult horizon/near-collinear edge falls back to a positive geometric
    /// upper bound, not a guessed branch. Such triangles are counted in the result.
    pub fn illuminance(
        &self, cx: &Cx<'_>, probe: &IlluminanceProbe,
    ) -> Result<AreaReferenceIlluminance, AreaReferenceError> {
        cx.checkpoint()?;
        let p = probe.position_m();
        let receiver = [p.x, p.y, p.z].map(point);
        let n = probe.normal();
        let normal = normalized([n.x, n.y, n.z].map(point))
            .ok_or(AreaReferenceError::NumericalRange)?;
        let mut range = zero();
        let mut coarse_triangles = 0;
        for (index, face) in self.faces.iter().enumerate() {
            cx.checkpoint()?;
            if face.luminance.hi() == 0.0 { continue; }
            let vectors = face.vertices.map(|v| difference(v, receiver));
            let first = normalized(vectors[0]).ok_or(AreaReferenceError::UnresolvedSourcePlane { triangle: index })?;
            let side = dot(face.perpendicular, first);
            if !finite(side) { return Err(AreaReferenceError::NumericalRange); }
            if side.lo() <= 0.0 && side.hi() >= 0.0 {
                return Err(AreaReferenceError::UnresolvedSourcePlane { triangle: index });
            }
            if !face.both_sides && side.lo() > 0.0 { continue; }
            let (factor, coarse) = projected_triangle(vectors, normal, face.area);
            coarse_triangles += usize::from(coarse);
            range = nonnegative(add(range, mul(face.luminance, factor)));
            if !finite(range) { return Err(AreaReferenceError::NumericalRange); }
        }
        let lux = range.lo() + (range.hi() - range.lo()) * 0.5;
        let error = sub(point(lux), point(range.lo())).hi()
            .max(sub(point(range.hi()), point(lux)).hi());
        if !error.is_finite() { return Err(AreaReferenceError::NumericalRange); }
        let mut hash = DomainHasher::new("org.frankensim.photometry.unoccluded-area.probe.v1");
        hash.update(&self.identity.0);
        for value in [p.x, p.y, p.z, n.x, n.y, n.z] { hash_float(&mut hash, value); }
        cx.checkpoint()?;
        Ok(AreaReferenceIlluminance {
            lux, numerical: NumericalCertificate::enclosure(range.lo(), range.hi()),
            numerical_error_bound_lux: error, coarse_triangles,
            numerical_model_identity: hash.finalize(), total_error_bound_lux: None,
        })
    }
}

/// An enclosure of the ideal unoccluded vacuum reference, not of a full scene.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AreaReferenceIlluminance {
    /// Deterministic midpoint in lux; not a native path sample or sample mean.
    pub lux: f64,
    /// Outward enclosure including area normalization, clipping and spectral sum.
    pub numerical: NumericalCertificate,
    /// Bound on midpoint error with respect to this conditional numerical model.
    pub numerical_error_bound_lux: f64,
    /// Triangles using a conservative geometric range instead of a resolved edge sum.
    pub coarse_triangles: usize,
    /// Selected source geometry, spectrum, numerical versions and original probe.
    pub numerical_model_identity: ContentHash,
    /// Unknown physical/model error; this implementation always returns `None`.
    pub total_error_bound_lux: Option<f64>,
}

impl AreaReferenceIlluminance {
    /// Output explicitly names the unoccluded, vacuum, direct-area scope.
    pub fn write_csv_header(writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(b"index,scope,lux,numerical_lower_lux,numerical_upper_lux,numerical_error_bound_lux,coarse_triangles,numerical_model_identity,calibration_error_bound_lux,spectral_discretization_error_bound_lux,physical_model_error_bound_lux,total_error_bound_lux,version\n")
    }
    /// Stream one conditional reference; I/O errors propagate to the caller.
    pub fn write_csv_row(self, index: usize, writer: &mut impl Write) -> io::Result<()> {
        writeln!(writer, "{index},unoccluded-vacuum-direct-area,{},{},{},{},{},{},unknown,unknown,unknown,unknown,{AREA_REFERENCE_SEMANTICS_VERSION}",
            self.lux, self.numerical.lo, self.numerical.hi, self.numerical_error_bound_lux,
            self.coarse_triangles, self.numerical_model_identity)
    }
}

/// Progressive receiver queries sharing one immutable source preparation.
pub struct UnoccludedAreaJob<'a> {
    reference: &'a UnoccludedAreaReference,
    receivers: Receivers<'a>,
    count: usize,
    values: Vec<AreaReferenceIlluminance>,
}
impl<'a> UnoccludedAreaJob<'a> {
    /// Borrow explicit probes without copying them.
    pub fn new(reference: &'a UnoccludedAreaReference, probes: &'a [IlluminanceProbe], cx: &Cx<'_>) -> Result<Self, AreaReferenceError> {
        Self::begin(reference, Receivers::Probes(probes), probes.len(), cx)
    }
    /// Use the original calculation plane and its row-major probe convention.
    pub fn for_grid(reference: &'a UnoccludedAreaReference, grid: CalculationGrid, cx: &Cx<'_>) -> Result<Self, AreaReferenceError> {
        Self::begin(reference, Receivers::Grid(grid), grid.len(), cx)
    }
    fn begin(reference: &'a UnoccludedAreaReference, receivers: Receivers<'a>, count: usize, cx: &Cx<'_>) -> Result<Self, AreaReferenceError> {
        cx.checkpoint()?;
        if count > reference.max_outputs { return Err(AreaReferenceError::Budget); }
        let mut values = Vec::new();
        values.try_reserve_exact(count).map_err(|_| AreaReferenceError::Allocation)?;
        cx.checkpoint()?;
        Ok(Self { reference, receivers, count, values })
    }
    /// Publish complete receivers only; errors preserve the accepted prefix.
    pub fn advance(&mut self, cx: &Cx<'_>, max_points: usize) -> Result<usize, AreaReferenceError> {
        cx.checkpoint()?;
        let start = self.values.len();
        let end = start.saturating_add(max_points).min(self.count);
        for index in start..end {
            let probe = match &self.receivers {
                Receivers::Probes(probes) => probes[index],
                Receivers::Grid(grid) => grid.probe(index).map_err(|_| AreaReferenceError::NumericalRange)?,
            };
            let value = self.reference.illuminance(cx, &probe)?;
            cx.checkpoint()?;
            self.values.push(value);
        }
        Ok(self.values.len() - start)
    }
    /// Only completed receiver results; uncomputed probes are not reported as zero.
    #[must_use]
    pub fn results(&self) -> &[AreaReferenceIlluminance] { &self.values }
    /// Every admitted receiver has been published.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.values.len() == self.count }
}

fn mesh<'a>(primitives: &'a [Primitive], source: &PhotometricAreaEmitter) -> Result<(&'a TriMesh, Option<RigidTransform>), AreaReferenceError> {
    if source.projector().is_some() { return Err(AreaReferenceError::UnsupportedSource); }
    let primitive = primitives.get(source.primitive_index()).ok_or(AreaReferenceError::InvalidGeometry)?;
    if primitive.emission.is_some() || matches!(primitive.material, Material::Dielectric { .. }) {
        return Err(AreaReferenceError::UnsupportedSource);
    }
    match &primitive.shape {
        Shape::Mesh(mesh) => Ok((mesh, None)),
        Shape::Instance(instance) => match instance.geometry() {
            SharedGeometry::Mesh(mesh) => Ok((mesh.as_ref(), Some(instance.transform()))),
            SharedGeometry::Chart(_) => Err(AreaReferenceError::UnsupportedSource),
        },
        _ => Err(AreaReferenceError::UnsupportedSource),
    }
}

fn transform(v: [f64; 3], placement: Option<RigidTransform>) -> Vector {
    let v = v.map(point);
    let Some(placement) = placement else { return v; };
    let [x, y, z, w] = placement.rotation_xyzw().map(point);
    let q = [x, y, z];
    let twice = cross(q, v).map(|x| mul(point(2.0), x));
    let double = cross(q, twice);
    let translation = placement.translation_m().map(point);
    std::array::from_fn(|i| add(add(add(v[i], mul(w, twice[i])), double[i]), translation[i]))
}

/// Stokes' theorem: integral_Omega omega dOmega is half the sum over edges of
/// theta * unit(a cross b). Clip the original planar triangle to n.r >= 0
/// BEFORE its spherical projection, so the receiving cosine is never negative.
fn projected_triangle(v: [Vector; 3], normal: Vector, area: Interval) -> (Interval, bool) {
    let upper = geometric_upper(v, area);
    let fallback = || (Interval::new(0.0, upper), true);
    let d = v.map(|p| dot(normal, p));
    if d.iter().all(|x| x.hi() <= 0.0) { return (zero(), false); }
    if d.iter().any(|x| x.lo() <= 0.0 && x.hi() >= 0.0 && *x != zero()) { return fallback(); }
    let mut polygon = [[zero(); 3]; 4];
    let mut count = 0;
    // Clipping a triangle against one half-plane creates at most four vertices.
    for end in 0..3 {
        let start = (end + 2) % 3;
        let inside_start = d[start].lo() >= 0.0;
        let inside_end = d[end].lo() >= 0.0;
        if inside_start != inside_end {
            let p = if d[start] == zero() { v[start] }
                else if d[end] == zero() { v[end] }
                else {
                    let t = d[start] / sub(d[start], d[end]);
                    if !finite(t) { return fallback(); }
                    let t = Interval::new(t.lo().max(0.0), t.hi().min(1.0));
                    std::array::from_fn(|i| add(mul(sub(point(1.0), t), v[start][i]), mul(t, v[end][i])))
                };
            if count == 0 || polygon[count - 1] != p {
                if count == 4 { return fallback(); }
                polygon[count] = p; count += 1;
            }
        }
        if inside_end && (count == 0 || polygon[count - 1] != v[end]) {
            if count == 4 { return fallback(); }
            polygon[count] = v[end]; count += 1;
        }
    }
    if count > 1 && polygon[0] == polygon[count - 1] { count -= 1; }
    if count < 3 { return (zero(), false); }
    for vertex in &mut polygon[..count] {
        let Some(direction) = normalized(*vertex) else { return fallback(); };
        *vertex = direction;
    }
    let mut boundary = [zero(); 3];
    for edge in 0..count {
        let a = polygon[edge];
        let b = polygon[(edge + 1) % count];
        let perpendicular = cross(a, b);
        let Some(sine) = norm(perpendicular) else { return fallback(); };
        if sine.lo() <= 0.0 { return fallback(); }
        let Some(angle) = atan2_upper_half_plane(sine, dot(a, b)) else { return fallback(); };
        for axis in 0..3 { boundary[axis] = add(boundary[axis], mul(perpendicular[axis] / sine, angle)); }
    }
    let projected = mul(point(0.5), dot(normal, boundary)).abs();
    match projected.intersect(Interval::new(0.0, upper)) {
        Some(bound) => (bound, false),
        // No narrowing on a failed numerical intersection.
        None => fallback(),
    }
}

fn geometric_upper(v: [Vector; 3], area: Interval) -> f64 {
    // Cosine-weighted solid angle cannot exceed pi; also dOmega <= dA/r².
    let box_vector = std::array::from_fn(|axis| v[0][axis].hull(v[1][axis]).hull(v[2][axis]));
    match norm(box_vector) {
        Some(distance) if distance.lo() > 0.0 => {
            let bound = (area / point(distance.lo()) / point(distance.lo())).hi();
            bound.min(pi().hi())
        }
        _ => pi().hi(),
    }
}

fn point(x: f64) -> Interval { Interval::point(x) }
fn zero() -> Interval { point(0.0) }
fn pi() -> Interval { Interval::new(fs_math::next_down(core::f64::consts::PI), fs_math::next_up(core::f64::consts::PI)) }
fn finite(x: Interval) -> bool { x.lo().is_finite() && x.hi().is_finite() }
fn nonnegative(x: Interval) -> Interval { Interval::new(x.lo().max(0.0), x.hi().max(0.0)) }
fn add(a: Interval, b: Interval) -> Interval { if a == zero() { b } else if b == zero() { a } else { a + b } }
fn sub(a: Interval, b: Interval) -> Interval {
    if b == zero() { a } else if a == b && a.lo() == a.hi() { zero() } else { a - b }
}
fn mul(a: Interval, b: Interval) -> Interval { if a == zero() || b == zero() { zero() } else { a * b } }
fn difference(a: Vector, b: Vector) -> Vector { std::array::from_fn(|i| sub(a[i], b[i])) }
fn dot(a: Vector, b: Vector) -> Interval { add(add(mul(a[0], b[0]), mul(a[1], b[1])), mul(a[2], b[2])) }
fn cross(a: Vector, b: Vector) -> Vector {
    [sub(mul(a[1], b[2]), mul(a[2], b[1])), sub(mul(a[2], b[0]), mul(a[0], b[2])), sub(mul(a[0], b[1]), mul(a[1], b[0]))]
}
fn magnitude_scale(vectors: &[Vector]) -> Option<f64> {
    let mut scale = 0.0_f64;
    for v in vectors { for &x in v { if !finite(x) { return None; } scale = scale.max(x.lo().abs()).max(x.hi().abs()); } }
    Some(scale)
}
fn norm(v: Vector) -> Option<Interval> {
    let scale = magnitude_scale(&[v])?;
    if scale == 0.0 { return Some(zero()); }
    let squares = v.map(|x| { let a = (x / point(scale)).abs(); nonnegative(mul(a, a)) });
    let result = mul(nonnegative(add(add(squares[0], squares[1]), squares[2])).sqrt(), point(scale));
    finite(result).then_some(nonnegative(result))
}
fn normalized(v: Vector) -> Option<Vector> {
    let length = norm(v)?;
    if length.lo() <= 0.0 { return None; }
    let unit = v.map(|x| if x == zero() { zero() } else { x / length });
    unit.iter().all(|x| finite(*x)).then_some(unit)
}
fn hash_float(hash: &mut DomainHasher, value: f64) { hash.update(&value.to_bits().to_le_bytes()); }
fn hash_range(hash: &mut DomainHasher, value: Interval) { hash_float(hash, value.lo()); hash_float(hash, value.hi()); }

#[cfg(test)]
mod tests;
