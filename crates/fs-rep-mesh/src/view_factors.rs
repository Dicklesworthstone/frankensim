//! Geometry-derived diffuse enclosure factors, with explicit sampling and
//! reciprocal/completeness projection. The caller owns the five-dimensional
//! unit-cube sampler (normally Owen-scrambled Sobol); this L2 kernel owns area
//! sampling, cosine-weighted rays, nearest opaque hits and geometric accounting.
//!
//! Triangle normals point INTO the radiating void. Unassigned triangles are
//! still opaque: hitting one refuses, it does not become an invisible obstacle.
//! No escaping ray is renormalized away. The symmetric exchange-area projection
//! is explicit, bounded, and reported alongside the original ray counts.
//!
//! These are floating-point ESTIMATES. Prefix agreement is not an error bound,
//! unseen openings/occlusion are not excluded, and projected conservation is not
//! proof of geometric accuracy. No surface may be called certified from this
//! result. Fixed inputs and samples replay on the same ISA/math profile; no
//! cross-ISA or invariance under remeshing is claimed.

use fs_blake3::{ContentHash, hash_domain};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

use crate::ray_triangle_watertight;

/// Resource and numerical policy, supplied before any ray work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViewFactorConfig {
    /// Even power of two in 256..=1,048,576, for EACH emitting patch.
    pub rays_per_surface: u32,
    /// Cap on ray/triangle intersection visits, including opaque unassigned faces.
    /// The complete all-pairs visit count is checked before tracing.
    pub max_triangle_tests: u64,
    /// At most 4096 simultaneous symmetric-scaling steps.
    pub max_balance_iterations: usize,
    /// Maximum absolute change between N/2-prefix and N raw factors; not a CI.
    pub max_sampling_change: f64,
    /// Maximum absolute change of any factor due to reciprocity/closure projection.
    pub max_factor_adjustment: f64,
}

/// Named refusals; failure never returns a partially usable factor matrix.
#[derive(Debug, Clone, PartialEq)]
pub enum ViewFactorError {
    /// An input or representable arithmetic requirement failed.
    Invalid(&'static str),
    /// The complete bounded trace exceeds the declared work allowance.
    Budget { required: u64, limit: u64 },
    /// Cancellation/deadline was observed at a bounded checkpoint.
    Interrupted,
    /// A ray left the represented enclosure; no implicit ambient row is added.
    Escape { surface: usize, ray: u32 },
    /// An opaque first hit has no radiative patch owner.
    UnassignedHit { surface: usize, ray: u32, face: usize },
    /// The receiving face points away from the void seen by the emitter.
    BackfaceHit { surface: usize, ray: u32, face: usize },
    /// The two nested sample counts disagree beyond the explicit allowance.
    SamplingChange { observed: f64, allowed: f64 },
    /// Symmetric area scaling did not close within the numerical work budget.
    Balance { iterations: usize, residual: f64 },
    /// Conservation projection would change the raw estimate too much.
    Adjustment { observed: f64, allowed: f64 },
}
impl core::fmt::Display for ViewFactorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "enclosure view factors: {self:?}")
    }
}
impl std::error::Error for ViewFactorError {}

type Result<T> = core::result::Result<T, ViewFactorError>;
const MAX_FACES: usize = 8192;
const MAX_VERTICES: usize = 32768;
const MAX_SURFACES: usize = 256;
const MAX_TESTS: u64 = 200_000_000;
const BALANCE_TOLERANCE: f64 = 1e-13;

/// An immutable estimate and the actual corrections which produced it.
#[derive(Debug, Clone)]
pub struct ViewFactorEstimate {
    geometry_identity: ContentHash,
    areas: Vec<f64>,
    counts: Vec<Vec<u32>>,
    raw: Vec<Vec<f64>>,
    factors: Vec<Vec<f64>>,
    config: ViewFactorConfig,
    sampling_change: f64,
    adjustment: f64,
    balance_iterations: usize,
    triangle_tests: u64,
}
impl ViewFactorEstimate {
    /// Exact input geometry, winding, owner order and patch-count identity.
    pub const fn geometry_identity(&self) -> ContentHash { self.geometry_identity }
    /// Integrated areas in the caller's patch order.
    pub fn areas(&self) -> &[f64] { &self.areas }
    /// Nearest-hit counts, before ANY reciprocity correction.
    pub fn counts(&self) -> &[Vec<u32>] { &self.counts }
    /// Count/N factors before ANY reciprocity correction.
    pub fn raw_factors(&self) -> &[Vec<f64>] { &self.raw }
    /// Symmetric-area-scaled factors; numerical conservation, not geometric proof.
    pub fn factors(&self) -> &[Vec<f64>] { &self.factors }
    /// Exact policy consumed by generation.
    pub const fn config(&self) -> ViewFactorConfig { self.config }
    /// Observed maximum difference between N and N/2 ray-prefix factors.
    pub const fn sampling_change(&self) -> f64 { self.sampling_change }
    /// Observed maximum absolute projected-minus-raw factor change.
    pub const fn adjustment(&self) -> f64 { self.adjustment }
    /// Number of completed symmetric-scaling updates.
    pub const fn balance_iterations(&self) -> usize { self.balance_iterations }
    /// Intersection visits actually performed.
    pub const fn triangle_tests(&self) -> u64 { self.triangle_tests }
}

fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| ViewFactorError::Interrupted)
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() { Ok(x) } else { Err(ViewFactorError::Invalid("nonfinite geometry/projection arithmetic")) }
}
fn dot(a: [f64;3], b: [f64;3]) -> f64 { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
fn sub(a: [f64;3], b: [f64;3]) -> [f64;3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }
fn cross(a: [f64;3], b: [f64;3]) -> [f64;3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn unit(a: [f64;3]) -> Result<[f64;3]> {
    let length = finite(dot(a,a).sqrt())?;
    if length <= 0.0 { return Err(ViewFactorError::Invalid("degenerate or unrepresentable triangle/frame")); }
    Ok([a[0]/length,a[1]/length,a[2]/length])
}
fn point(a: [f64;3]) -> Point3 { Point3::new(a[0],a[1],a[2]) }
struct Face { points: [Point3;3], xyz: [[f64;3];3], normal: [f64;3], tangent: [f64;3], bitangent: [f64;3] }

/// Hash the complete ordered geometry and ownership, without interpreting it
/// as a topology certificate. The same defensive shape caps as tracing apply.
pub fn geometry_identity(positions: &[[f64;3]], triangles: &[[u32;3]],
    owners: &[Option<usize>], surfaces: usize) -> Result<ContentHash> {
    if positions.is_empty() || positions.len()>MAX_VERTICES || triangles.is_empty()
        || triangles.len()>MAX_FACES || triangles.len()!=owners.len()
        || !(2..=MAX_SURFACES).contains(&surfaces) {
        return Err(ViewFactorError::Invalid("geometry shape/cap or patch count is invalid"));
    }
    let mut bytes=Vec::new();
    for n in [positions.len(),triangles.len(),surfaces] { bytes.extend_from_slice(&(n as u64).to_le_bytes()); }
    for xyz in positions { for &x in xyz { finite(x)?; bytes.extend_from_slice(&x.to_bits().to_le_bytes()); } }
    for (triangle,owner) in triangles.iter().zip(owners) {
        let mut unique=*triangle; unique.sort_unstable();
        if unique[0]==unique[1] || unique[1]==unique[2] || unique[2] as usize>=positions.len()
            || owner.is_some_and(|i|i>=surfaces) {
            return Err(ViewFactorError::Invalid("invalid triangle indices or patch owner"));
        }
        for v in triangle { bytes.extend_from_slice(&v.to_le_bytes()); }
        bytes.extend_from_slice(&owner.map_or(u64::MAX, |i|i as u64).to_le_bytes());
    }
    Ok(hash_domain("org.frankensim.fs-rep-mesh.view-factor-geometry.v1", &bytes))
}

/// Trace a complete closed diffuse enclosure using a caller-owned unit-cube
/// sequence. Draws are indexed by (patch, ray), not worker or execution order.
/// The five coordinates must lie STRICTLY in (0,1): patch-area choice, two
/// triangle-area coordinates, and two cosine-hemisphere coordinates.
///
/// Every triangle, including an unassigned one, is an opaque occluder. Normals
/// must face the radiating void. The full opaque geometry, not just emitters,
/// must be supplied. Standard binary64 elementary functions are used; this is
/// a same-profile numerical model, not a cross-ISA deterministic certificate.
#[allow(clippy::too_many_arguments)]
pub fn estimate_view_factors(cx: &Cx<'_>, positions: &[[f64;3]], triangles: &[[u32;3]],
    owners: &[Option<usize>], surfaces: usize, config: ViewFactorConfig,
    mut sample: impl FnMut(usize,u32)->[f64;5]) -> Result<ViewFactorEstimate> {
    poll(cx)?;
    let identity=geometry_identity(positions,triangles,owners,surfaces)?;
    let rays=config.rays_per_surface;
    if !(256..=1_048_576).contains(&rays) || !rays.is_power_of_two()
        || !(1..=4096).contains(&config.max_balance_iterations)
        || !(config.max_sampling_change.is_finite() && config.max_sampling_change>0.0 && config.max_sampling_change<1.0)
        || !(config.max_factor_adjustment.is_finite() && config.max_factor_adjustment>0.0 && config.max_factor_adjustment<1.0)
        || config.max_triangle_tests==0 || config.max_triangle_tests>MAX_TESTS {
        return Err(ViewFactorError::Invalid("invalid ray, projection or work policy"));
    }
    let required=(surfaces as u64).checked_mul(u64::from(rays))
        .and_then(|n|n.checked_mul(triangles.len().saturating_sub(1) as u64))
        .ok_or(ViewFactorError::Invalid("ray work overflow"))?;
    if required>config.max_triangle_tests { return Err(ViewFactorError::Budget {required,limit:config.max_triangle_tests}); }
    let mut faces=Vec::with_capacity(triangles.len());
    let mut areas=vec![0.0;surfaces];
    let mut selections=vec![Vec::new();surfaces];
    for (slot,vertices) in triangles.iter().enumerate() {
        poll(cx)?;
        let xyz=vertices.map(|v|positions[v as usize]);
        let edge=sub(xyz[1],xyz[0]);
        let normal_raw=cross(edge,sub(xyz[2],xyz[0]));
        let area=finite(0.5*dot(normal_raw,normal_raw).sqrt())?;
        if area<=0.0 { return Err(ViewFactorError::Invalid("zero or unrepresentable triangle area")); }
        let normal=unit(normal_raw)?;
        let tangent=unit(edge)?;
        let bitangent=unit(cross(normal,tangent))?;
        faces.push(Face {points:xyz.map(point),xyz,normal,tangent,bitangent});
        if let Some(owner)=owners[slot] {
            let next=finite(areas[owner]+area)?;
            if next<=areas[owner] {return Err(ViewFactorError::Invalid("patch area accumulation lost a face"));}
            areas[owner]=next; selections[owner].push((slot,next));
        }
    }
    if areas.iter().any(|a|*a<=0.0) {return Err(ViewFactorError::Invalid("empty emitting patch"));}
    let mut counts=vec![vec![0_u32;surfaces];surfaces];
    let mut prefix=vec![vec![0_u32;surfaces];surfaces];
    let mut visits=0_u64;
    for surface in 0..surfaces {
        for ray in 0..rays {
            poll(cx)?;
            let u=sample(surface,ray);
            if u.iter().any(|v|!v.is_finite() || *v<=0.0 || *v>=1.0) {
                return Err(ViewFactorError::Invalid("sampler must return five coordinates strictly in (0,1)"));
            }
            let selected=selections[surface].partition_point(|(_,end)|*end<=u[0]*areas[surface]);
            let source=selections[surface].get(selected).ok_or(ViewFactorError::Invalid("area sample rounded past patch"))?.0;
            let face=&faces[source];
            let root=u[1].sqrt(); let bary=[1.0-root,root*(1.0-u[2]),root*u[2]];
            let mut xyz=[0.0;3];
            for (k,x) in xyz.iter_mut().enumerate() {
                *x=finite(bary[0]*face.xyz[0][k]+bary[1]*face.xyz[1][k]+bary[2]*face.xyz[2][k])?;
            }
            let radial=u[3].sqrt(); let axial=(1.0-u[3]).sqrt();
            let (sine,cosine)=(core::f64::consts::TAU*u[4]).sin_cos();
            let local=[radial*cosine,radial*sine,axial];
            let direction=unit(core::array::from_fn(|k|local[0]*face.tangent[k]+local[1]*face.bitangent[k]+local[2]*face.normal[k]))?;
            let origin=point(xyz); let dir=Vec3::new(direction[0],direction[1],direction[2]);
            let mut nearest=None; let mut distance=f64::INFINITY;
            for (slot,other) in faces.iter().enumerate() {
                if slot==source {continue;}
                if visits%256==0 {poll(cx)?;} visits+=1;
                if let Some(t)=ray_triangle_watertight(origin,dir,other.points[0],other.points[1],other.points[2])
                    && t>0.0 && t.is_finite() && t<distance {
                    distance=t; nearest=Some(slot);
                }
            }
            let hit=nearest.ok_or(ViewFactorError::Escape {surface,ray})?;
            let owner=owners[hit].ok_or(ViewFactorError::UnassignedHit {surface,ray,face:hit})?;
            if dot(faces[hit].normal,direction)>=0.0 {return Err(ViewFactorError::BackfaceHit {surface,ray,face:hit});}
            counts[surface][owner]+=1;
            if ray<rays/2 {prefix[surface][owner]+=1;}
        }
    }
    let raw:Vec<Vec<f64>>=counts.iter().map(|row|row.iter().map(|&n|f64::from(n)/f64::from(rays)).collect()).collect();
    let mut change=0.0_f64;
    for i in 0..surfaces {for j in 0..surfaces {change=change.max((raw[i][j]-2.0*f64::from(prefix[i][j])/f64::from(rays)).abs());}}
    if change>config.max_sampling_change {return Err(ViewFactorError::SamplingChange {observed:change,allowed:config.max_sampling_change});}
    let (factors,iterations)=balance(cx,&areas,&raw,config.max_balance_iterations)?;
    let mut adjustment=0.0_f64;
    for i in 0..surfaces {for j in 0..surfaces {adjustment=adjustment.max((factors[i][j]-raw[i][j]).abs());}}
    if adjustment>config.max_factor_adjustment {return Err(ViewFactorError::Adjustment {observed:adjustment,allowed:config.max_factor_adjustment});}
    poll(cx)?;
    Ok(ViewFactorEstimate {geometry_identity:identity,areas,counts,raw,factors,config,
        sampling_change:change,adjustment,balance_iterations:iterations,triangle_tests:visits})
}

// Symmetric nonnegative exchange area G, followed by simultaneous diagonal
// scaling X G X to the measured area margins. Unlike row normalization, this
// preserves reciprocity AND the observed zero pattern. It can fail for an
// unresolved/incompatible sampled support; no diagonal or environment is added.
fn balance(cx:&Cx<'_>,areas:&[f64],raw:&[Vec<f64>],max_iterations:usize)->Result<(Vec<Vec<f64>>,usize)> {
    let n=areas.len(); let scale=areas.iter().copied().fold(0.0,f64::max);
    let target:Vec<_>=areas.iter().map(|a|a/scale).collect();
    if target.iter().any(|&a|a<=0.0 || !a.is_finite()) {return Err(ViewFactorError::Invalid("unrepresentable area ratios"));}
    let mut exchange=vec![vec![0.0;n];n];
    for i in 0..n {for j in i..n {
        let g=finite(0.5*target[i]*raw[i][j]+0.5*target[j]*raw[j][i])?;
        exchange[i][j]=g; exchange[j][i]=g;
    }}
    let mut diagonal=vec![1.0;n];
    for iteration in 0..=max_iterations {
        poll(cx)?;
        let mut ratios=vec![0.0;n]; let mut residual=0.0_f64;
        for i in 0..n {
            let mut row=0.0;
            for j in 0..n {row=finite(row+finite(exchange[i][j]*diagonal[j])?)?;}
            row=finite(row*diagonal[i])?;
            if row<=0.0 {return Err(ViewFactorError::Invalid("empty reciprocal support or underflow"));}
            ratios[i]=finite(target[i]/row)?;
            residual=residual.max(finite(row/target[i]-1.0)?.abs());
        }
        if residual<=BALANCE_TOLERANCE {
            let mut result=vec![vec![0.0;n];n];
            for i in 0..n {for j in i..n {
                let g=finite(diagonal[i]*exchange[i][j]*diagonal[j])?;
                result[i][j]=finite(g/target[i])?;
                result[j][i]=finite(g/target[j])?;
                if !(0.0..=1.0).contains(&result[i][j]) || !(0.0..=1.0).contains(&result[j][i]) {
                    return Err(ViewFactorError::Invalid("projected factor outside [0,1]"));
                }
            }}
            return Ok((result,iteration));
        }
        if iteration==max_iterations {return Err(ViewFactorError::Balance {iterations:iteration,residual});}
        for i in 0..n {diagonal[i]=finite(diagonal[i]*ratios[i].sqrt())?;}
    }
    unreachable!("bounded projection exits at its final iteration")
}

#[cfg(test)]
mod tests;
