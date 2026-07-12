//! SPECTRAL PATH TRACER v1 (bead 872c, WS3; [F] — behind the `tracer`
//! feature): hero-wavelength Monte-Carlo transport with next-event
//! estimation + BSDF sampling combined by the crate's MIS balance
//! heuristic, over the certified chart/BVH backends, producing CIE XYZ
//! film and byte-exact EXR through fs-img.
//!
//! DETERMINISM: every random draw comes from a counter-based stream
//! keyed by `(pixel, sample, bounce-dimension)` — Philox 4×32-10 for
//! path decisions, optionally Owen-scrambled Sobol' for the pixel/
//! wavelength dimensions ([`Sampler::OwenSobol`], decorrelated across
//! pixels by a Philox-derived scramble seed). No draw depends on
//! scheduling, so images are bitwise invariant to tile traversal and
//! worker count, and a render RESUMED from an `spp` checkpoint equals
//! the straight-through render bitwise (the pause–serialize–resume
//! doctrine applied to images). All transcendentals in the radiance
//! path go through `fs_math::det` (goldens hash these bits — no
//! platform libm), and Fresnel/roughness powers are explicit
//! multiplications, never `powi` (the a55x/4xnt hazard class).
//!
//! v1 scope (documented, falsifiable): single rectangular area light
//! per scene for NEE; lights are also scene geometry so BSDF paths
//! find them (MIS-weighted both ways); materials are Lambertian and
//! GGX (Smith separable G, Schlick Fresnel with the spectral
//! reflectance as F0); no volumetric media (the `volumes` module is
//! separate); no environment light; no Russian roulette (fixed depth
//! keeps work deterministic).

use crate::charts::{Hit, Ray, TraceTermination, TriMesh, sphere_trace};
use crate::spectral::{
    LAMBDA_MAX, LAMBDA_MIN, LiftedSpectrum, cie_x, cie_y, cie_z, xyz_e_to_d65, xyz_to_linear_srgb,
    y_integral,
};
use crate::{balance_heuristic, hero_wavelengths};
use fs_exec::{Cancelled, Cx};
use fs_geom::{Chart, Point3, Vec3};
use fs_math::det;
use fs_rand::philox::philox4x32_10;
use fs_rand::qmc::Sobol;

/// Bit-affecting semantic surface version of the tracer (see
/// golden-couplings.json): the path-integrator estimator shape, the
/// Philox/Sobol stream keying, the BSDF forms, the CMF/adaptation
/// constants it inherits from `spectral`, and the EXR channel layout.
/// Bump ONLY with a semantic justification per docs/GOLDEN_POLICY.md.
pub const TRACER_BIT_SEMANTICS_VERSION: u32 = 1;

const PI: f64 = core::f64::consts::PI;
/// Hero-wavelength packet width (the bead's 4-wavelength packets).
pub const PACKET: usize = 4;
/// Self-intersection offset along the normal when spawning rays.
const RAY_EPS: f64 = 1e-6;
/// Sphere-trace surface tolerance.
const TRACE_EPS: f64 = 1e-7;

/// The per-draw uniform stream: Philox keyed by (pixel, sample,
/// dimension). Counter-based — random access, no state shared between
/// pixels/samples/workers.
struct PathRng {
    pixel: u32,
    sample: u32,
    dim: u32,
    key: [u32; 2],
}

impl PathRng {
    fn next2(&mut self) -> (f64, f64) {
        let out = philox4x32_10([self.pixel, self.sample, self.dim, 0x7261_7972], self.key);
        self.dim += 1;
        (u32_unit(out[0]), u32_unit(out[1]))
    }
}

fn u32_unit(x: u32) -> f64 {
    f64::from(x) / 4_294_967_296.0
}

/// Pixel-space sampler for the (jitter-x, jitter-y, hero-λ) dimensions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sampler {
    /// Independent Philox draws for every dimension.
    Iid,
    /// Owen-scrambled Sobol' over the three pixel dimensions,
    /// decorrelated across pixels by a Philox-derived scramble seed
    /// (the ambition-round upgrade; its equal-spp variance claim is
    /// measured in the battery, not assumed).
    OwenSobol,
}

/// How direct lighting is estimated — [`DirectStrategy::Mis`] is the
/// product setting; the single-technique modes exist so the battery
/// can MEASURE that MIS beats either alone (the bead's acceptance).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirectStrategy {
    /// Next-event estimation only.
    NeeOnly,
    /// BSDF sampling only.
    BsdfOnly,
    /// Both, combined with the balance heuristic.
    Mis,
}

/// A surface material.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Material {
    /// Ideal diffuse with a spectral reflectance.
    Lambertian {
        /// Reflectance spectrum (bounded (0,1) by construction).
        reflectance: LiftedSpectrum,
    },
    /// GGX microfacet (Smith separable shadowing, Schlick Fresnel with
    /// the spectral reflectance as F0).
    Ggx {
        /// F0 reflectance spectrum.
        reflectance: LiftedSpectrum,
        /// Roughness α (GGX convention, > 0).
        alpha: f64,
    },
}

/// Scene geometry: a triangle mesh (BVH) or any certified chart
/// (sphere-traced SDF/F-rep through the default [S] backend surface hardened by
/// bead 8ll9).
pub enum Shape {
    /// Triangle mesh over the deterministic median-split BVH.
    Mesh(TriMesh),
    /// A certified-Lipschitz chart, sphere-traced.
    Chart(Box<dyn Chart>),
}

/// One scene object.
pub struct Primitive {
    /// Geometry.
    pub shape: Shape,
    /// Material (ignored for pure emitters in v1: lights do not
    /// reflect).
    pub material: Material,
    /// Emitted radiance: spectrum × scale (None = non-emissive).
    pub emission: Option<(LiftedSpectrum, f64)>,
}

/// The single rectangular area light (v1) used by next-event
/// estimation. The SAME rectangle must also be present as an emissive
/// mesh primitive (index `prim`) so BSDF-sampled paths hit it.
pub struct RectLight {
    /// One corner.
    pub corner: Point3,
    /// First edge.
    pub edge_u: Vec3,
    /// Second edge.
    pub edge_v: Vec3,
    /// Index of the emissive primitive this light corresponds to.
    pub prim: usize,
    /// Emitted radiance spectrum × scale (must match the primitive's).
    pub emission: (LiftedSpectrum, f64),
}

impl RectLight {
    fn area(&self) -> f64 {
        cross(self.edge_u, self.edge_v).norm()
    }

    fn normal(&self) -> Vec3 {
        let n = cross(self.edge_u, self.edge_v);
        n.scale(1.0 / n.norm())
    }
}

/// Pinhole camera. `half_tan` is tan(fov/2) supplied directly — the
/// library takes no trig on its API surface.
pub struct Camera {
    /// Eye point.
    pub eye: Point3,
    /// Unit view direction.
    pub forward: Vec3,
    /// Unit up (orthogonal to forward).
    pub up: Vec3,
    /// tan(vertical fov / 2).
    pub half_tan: f64,
}

/// A renderable scene.
pub struct Scene {
    /// Objects (lights included as emissive primitives).
    pub primitives: Vec<Primitive>,
    /// The NEE light (v1: exactly one).
    pub light: RectLight,
    /// Camera.
    pub camera: Camera,
}

/// Render settings.
#[derive(Debug, Clone, Copy)]
pub struct Settings {
    /// Image width (pixels).
    pub width: u32,
    /// Image height (pixels).
    pub height: u32,
    /// Samples per pixel.
    pub spp: u32,
    /// Maximum path depth (bounces).
    pub max_depth: u32,
    /// Pixel-dimension sampler.
    pub sampler: Sampler,
    /// Direct-lighting strategy.
    pub strategy: DirectStrategy,
    /// Stream seed (the replay identity).
    pub seed: u64,
}

/// Fail-closed spectral-tracer diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TracerError {
    /// The supplied execution context requested cancellation.
    Cancelled,
    /// A chart backend stopped in a state other than a clean miss or certified
    /// residual hit.
    BackendFailure(TraceTermination),
    /// A chart returned a terminal result without retaining its typed
    /// no-tunneling claim. Uncertified misses are not geometry absence.
    UncertifiedTrace,
    /// A progressive sample range had its exclusive end before its start.
    InvalidRange { from: u32, to: u32 },
    /// Render dimensions or a public film buffer were structurally invalid or
    /// could not be allocated without exceeding address-space bounds.
    InvalidInput,
    /// Shading requires a finite surface normal; no arbitrary fallback normal
    /// may be minted.
    MissingNormal,
}

impl core::fmt::Display for TracerError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => formatter.write_str("spectral render cancelled"),
            Self::BackendFailure(termination) => {
                write!(formatter, "chart backend stopped with {termination:?}")
            }
            Self::UncertifiedTrace => {
                formatter.write_str("chart backend produced an uncertified trace result")
            }
            Self::InvalidRange { from, to } => {
                write!(formatter, "invalid progressive sample range {from}..{to}")
            }
            Self::InvalidInput => formatter.write_str("invalid spectral render input"),
            Self::MissingNormal => formatter.write_str("surface hit has no finite normal"),
        }
    }
}

impl core::error::Error for TracerError {}

impl From<Cancelled> for TracerError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

/// Accumulated CIE XYZ film: `spp` samples summed per pixel (divide on
/// output). Checkpointable: rendering samples `[a, b)` then `[b, c)`
/// into the same film equals rendering `[a, c)` bitwise.
#[derive(Debug, Clone, PartialEq)]
pub struct Film {
    /// Width.
    pub width: u32,
    /// Height.
    pub height: u32,
    /// Row-major XYZ sums.
    pub xyz: Vec<[f64; 3]>,
    /// Samples accumulated so far.
    pub spp_done: u32,
}

impl Film {
    /// An empty film.
    #[must_use]
    pub fn new(width: u32, height: u32) -> Film {
        Self::try_new(width, height).expect("film dimensions must fit the address space")
    }

    /// Allocate an empty film without panicking on malformed or excessive
    /// dimensions.
    pub fn try_new(width: u32, height: u32) -> Result<Film, TracerError> {
        if width == 0 || height == 0 {
            return Err(TracerError::InvalidInput);
        }
        let len = (width as usize)
            .checked_mul(height as usize)
            .ok_or(TracerError::InvalidInput)?;
        let mut xyz = Vec::new();
        xyz.try_reserve_exact(len)
            .map_err(|_| TracerError::InvalidInput)?;
        xyz.resize(len, [0.0; 3]);
        Ok(Film {
            width,
            height,
            xyz,
            spp_done: 0,
        })
    }

    /// Linear-sRGB planes (R, G, B row-major), Bradford-adapted like
    /// the rest of the spectral pipeline; sums divided by `spp_done`.
    #[must_use]
    pub fn to_linear_srgb(&self) -> [Vec<f32>; 3] {
        let n = self.xyz.len();
        let mut planes = [vec![0.0f32; n], vec![0.0f32; n], vec![0.0f32; n]];
        let inv = if self.spp_done == 0 {
            0.0
        } else {
            1.0 / f64::from(self.spp_done)
        };
        for (i, xyz) in self.xyz.iter().enumerate() {
            let rgb = xyz_to_linear_srgb(xyz_e_to_d65([xyz[0] * inv, xyz[1] * inv, xyz[2] * inv]));
            for (p, v) in planes.iter_mut().zip(rgb) {
                #[allow(clippy::cast_possible_truncation)]
                {
                    p[i] = v as f32;
                }
            }
        }
        planes
    }
}

/// Render samples `[from, to)` for every pixel into `film` (progressive
/// accumulation; `film.spp_done` must equal `from`).
///
/// # Panics
/// If the film dimensions or checkpoint do not match. A malformed public XYZ
/// buffer or invalid range is returned as a structured error.
pub fn render_range(
    scene: &Scene,
    cx: &Cx<'_>,
    s: &Settings,
    film: &mut Film,
    from: u32,
    to: u32,
) -> Result<(), TracerError> {
    assert_eq!((film.width, film.height), (s.width, s.height), "film shape");
    if film.width == 0 || film.height == 0 || s.width == 0 || s.height == 0 {
        return Err(TracerError::InvalidInput);
    }
    let expected_len = (film.width as usize)
        .checked_mul(film.height as usize)
        .ok_or(TracerError::InvalidInput)?;
    if film.xyz.len() != expected_len {
        return Err(TracerError::InvalidInput);
    }
    assert_eq!(film.spp_done, from, "progressive checkpoint mismatch");
    if to < from {
        return Err(TracerError::InvalidRange { from, to });
    }
    cx.checkpoint()?;
    if to == from {
        return Ok(());
    }
    let key = [(s.seed & 0xffff_ffff) as u32, (s.seed >> 32) as u32];
    let sobol = Sobol::scrambled(3, s.seed);
    let kn = 1.0 / y_integral();
    // Cancellation and backend refusals are transactional: a failed range
    // leaves both the accumulated sums and checkpoint unchanged, so retrying
    // cannot double-count a partially completed range.
    let mut staged_xyz = Vec::new();
    staged_xyz
        .try_reserve_exact(film.xyz.len())
        .map_err(|_| TracerError::InvalidInput)?;
    for chunk in film.xyz.chunks(4096) {
        cx.checkpoint()?;
        staged_xyz.extend_from_slice(chunk);
    }
    for py in 0..s.height {
        cx.checkpoint()?;
        for px in 0..s.width {
            let pixel = py * s.width + px;
            let slot = &mut staged_xyz[pixel as usize];
            for sample in from..to {
                cx.checkpoint()?;
                let (jx, jy, ul) = pixel_dims(s, &sobol, key, pixel, sample);
                let xyz = trace_path(scene, cx, s, kn, pixel, sample, jx, jy, ul)?;
                slot[0] += xyz[0];
                slot[1] += xyz[1];
                slot[2] += xyz[2];
            }
        }
    }
    cx.checkpoint()?;
    film.xyz = staged_xyz;
    film.spp_done = to;
    Ok(())
}

/// Render the full image (fresh film, samples `[0, spp)`).
pub fn render(scene: &Scene, cx: &Cx<'_>, s: &Settings) -> Result<Film, TracerError> {
    cx.checkpoint()?;
    let mut film = Film::try_new(s.width, s.height)?;
    render_range(scene, cx, s, &mut film, 0, s.spp)?;
    Ok(film)
}

/// The (jitter-x, jitter-y, hero-λ) dimensions for one (pixel, sample).
fn pixel_dims(
    s: &Settings,
    sobol: &Sobol,
    key: [u32; 2],
    pixel: u32,
    sample: u32,
) -> (f64, f64, f64) {
    match s.sampler {
        Sampler::Iid => {
            let a = philox4x32_10([pixel, sample, 0xdead_0001, 0], key);
            (u32_unit(a[0]), u32_unit(a[1]), u32_unit(a[2]))
        }
        Sampler::OwenSobol => {
            // One Sobol' point per sample index; Cranley–Patterson-free
            // decorrelation across pixels via a per-pixel Philox shift
            // of the SAMPLE INDEX ordering is not net-preserving, so
            // instead the scramble seed is shared and the pixel enters
            // through a Philox-derived toroidal shift of the point —
            // net-preserving per pixel, decorrelated across pixels.
            let mut pt = [0.0f64; 3];
            sobol.point(sample, &mut pt);
            let shift = philox4x32_10([pixel, 0x50b0_1000, 0, 0], key);
            let wrap = |x: f64, u: u32| {
                let v = x + u32_unit(u);
                if v >= 1.0 { v - 1.0 } else { v }
            };
            (
                wrap(pt[0], shift[0]),
                wrap(pt[1], shift[1]),
                wrap(pt[2], shift[2]),
            )
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // one integrator, one story
fn trace_path(
    scene: &Scene,
    cx: &Cx<'_>,
    s: &Settings,
    kn: f64,
    pixel: u32,
    sample: u32,
    jx: f64,
    jy: f64,
    ul: f64,
) -> Result<[f64; 3], TracerError> {
    let key = [(s.seed & 0xffff_ffff) as u32, (s.seed >> 32) as u32];
    let mut rng = PathRng {
        pixel,
        sample,
        dim: 1,
        key,
    };
    // Hero wavelengths: one stratified draw covers the packet.
    let hero = LAMBDA_MIN + ul * (LAMBDA_MAX - LAMBDA_MIN);
    let lambdas = hero_wavelengths(hero, PACKET, LAMBDA_MIN, LAMBDA_MAX);
    // Camera ray.
    let px = pixel % s.width;
    let py = pixel / s.width;
    let (w, h) = (f64::from(s.width), f64::from(s.height));
    let ndc_x = (2.0 * (f64::from(px) + jx) / w - 1.0) * s.camera_aspect() * scene.camera.half_tan;
    let ndc_y = (1.0 - 2.0 * (f64::from(py) + jy) / h) * scene.camera.half_tan;
    let right = cross(scene.camera.forward, scene.camera.up);
    let dir = unit(Vec3::new(
        scene.camera.forward.x + ndc_x * right.x + ndc_y * scene.camera.up.x,
        scene.camera.forward.y + ndc_x * right.y + ndc_y * scene.camera.up.y,
        scene.camera.forward.z + ndc_x * right.z + ndc_y * scene.camera.up.z,
    ));
    let mut ray = Ray {
        origin: scene.camera.eye,
        dir,
    };
    let mut throughput = [1.0f64; PACKET];
    let mut radiance = [0.0f64; PACKET];
    let mut prev_bsdf_pdf: Option<f64> = None;
    let mut prev_origin = ray.origin;
    for _depth in 0..s.max_depth {
        cx.checkpoint()?;
        let Some((prim_idx, hit)) = intersect(scene, cx, &ray)? else {
            break;
        };
        let prim = &scene.primitives[prim_idx];
        let n = oriented_normal(&hit, &ray)?;
        if let Some((spec, scale)) = &prim.emission {
            // MIS weight against NEE for this light, seen from the
            // previous vertex.
            #[allow(clippy::match_same_arms)] // distinct reasons for the same weight
            let weight = match (s.strategy, prev_bsdf_pdf) {
                (_, None) => 1.0, // camera ray: only BSDF "found" it
                (DirectStrategy::BsdfOnly, _) => 1.0,
                (DirectStrategy::NeeOnly, _) => 0.0,
                (DirectStrategy::Mis, Some(bp)) => {
                    let d = hit.point.delta_from(prev_origin);
                    let d2 = d.dot(d);
                    let cos_l = scene.light.normal().dot(unit(d)).abs();
                    let pdf_nee = if cos_l > 1e-12 {
                        d2 / (cos_l * scene.light.area())
                    } else {
                        0.0
                    };
                    balance_heuristic(1, bp, 1, pdf_nee)
                }
            };
            for (k, &l) in lambdas.iter().enumerate() {
                radiance[k] += throughput[k] * spec.eval(l) * scale * weight;
            }
            break; // v1: emitters do not reflect
        }
        // Next-event estimation.
        if s.strategy != DirectStrategy::BsdfOnly {
            let (u1, u2) = rng.next2();
            let q = scene
                .light
                .corner
                .offset(scene.light.edge_u.scale(u1))
                .offset(scene.light.edge_v.scale(u2));
            let to_light = q.delta_from(hit.point);
            let d2 = to_light.dot(to_light);
            let dist = d2.sqrt();
            let wi = to_light.scale(1.0 / dist);
            let cos_s = n.dot(wi);
            let cos_l = scene.light.normal().dot(wi).abs();
            if cos_s > 0.0 && cos_l > 1e-9 {
                let shadow = Ray {
                    origin: hit.point.offset(n.scale(RAY_EPS)),
                    dir: wi,
                };
                let vis = match intersect(scene, cx, &shadow)? {
                    Some((i, h)) => i == scene.light.prim && h.t > dist - 1e-4,
                    None => false,
                };
                if vis {
                    let pdf_nee = d2 / (cos_l * scene.light.area());
                    let wo = ray.dir.scale(-1.0);
                    let bsdf_pdf = bsdf_pdf(&prim.material, n, wo, wi);
                    let weight = match s.strategy {
                        DirectStrategy::Mis => balance_heuristic(1, pdf_nee, 1, bsdf_pdf),
                        _ => 1.0,
                    };
                    let (espec, escale) = &scene.light.emission;
                    for (k, &l) in lambdas.iter().enumerate() {
                        let f = bsdf_eval(&prim.material, n, wo, wi, l);
                        radiance[k] +=
                            throughput[k] * f * cos_s * espec.eval(l) * escale / pdf_nee * weight;
                    }
                }
            }
        }
        // BSDF sampling for the next bounce.
        let (u1, u2) = rng.next2();
        let wo = ray.dir.scale(-1.0);
        let Some((wi, pdf)) = bsdf_sample(&prim.material, n, wo, u1, u2) else {
            break;
        };
        let cos_s = n.dot(wi).max(0.0);
        if pdf <= 0.0 || cos_s <= 0.0 {
            break;
        }
        for (k, &l) in lambdas.iter().enumerate() {
            throughput[k] *= bsdf_eval(&prim.material, n, wo, wi, l) * cos_s / pdf;
        }
        prev_bsdf_pdf = Some(pdf);
        prev_origin = hit.point;
        ray = Ray {
            origin: hit.point.offset(n.scale(RAY_EPS)),
            dir: wi,
        };
    }
    // Hero-wavelength estimator → XYZ (same normalization convention
    // as `spectral::xyz_of_spectrum`: Y of unit radiance is 1).
    let range = LAMBDA_MAX - LAMBDA_MIN;
    let mut xyz = [0.0f64; 3];
    for (k, &l) in lambdas.iter().enumerate() {
        let w = radiance[k] * range / PACKET as f64 * kn;
        xyz[0] += w * cie_x(l);
        xyz[1] += w * cie_y(l);
        xyz[2] += w * cie_z(l);
    }
    Ok(xyz)
}

impl Settings {
    fn camera_aspect(&self) -> f64 {
        f64::from(self.width) / f64::from(self.height)
    }
}

fn intersect(scene: &Scene, cx: &Cx<'_>, ray: &Ray) -> Result<Option<(usize, Hit)>, TracerError> {
    let mut best: Option<(usize, Hit)> = None;
    for (i, prim) in scene.primitives.iter().enumerate() {
        cx.checkpoint()?;
        let hit = match &prim.shape {
            Shape::Mesh(mesh) => mesh.intersect_with_cx(cx, ray)?,
            Shape::Chart(chart) => {
                let (hit, audit) = sphere_trace(chart.as_ref(), cx, ray, 1e4, TRACE_EPS, 1.0);
                if matches!(
                    audit.termination,
                    TraceTermination::Hit | TraceTermination::Miss
                ) && !audit.certified
                {
                    return Err(TracerError::UncertifiedTrace);
                }
                match audit.termination {
                    TraceTermination::Cancelled => return Err(TracerError::Cancelled),
                    TraceTermination::Miss => None,
                    TraceTermination::Hit => {
                        Some(hit.ok_or(TracerError::BackendFailure(TraceTermination::Hit))?)
                    }
                    termination => return Err(TracerError::BackendFailure(termination)),
                }
            }
        };
        if let Some(h) = hit
            && best.as_ref().is_none_or(|(_, bh)| h.t < bh.t)
        {
            best = Some((i, h));
        }
    }
    Ok(best)
}

fn oriented_normal(hit: &Hit, ray: &Ray) -> Result<Vec3, TracerError> {
    let n = hit.normal.ok_or(TracerError::MissingNormal)?;
    if !n.x.is_finite() || !n.y.is_finite() || !n.z.is_finite() || n.norm() <= 0.0 {
        return Err(TracerError::MissingNormal);
    }
    Ok(if n.dot(ray.dir) > 0.0 {
        n.scale(-1.0)
    } else {
        n
    })
}

// ---- BSDF machinery --------------------------------------------------

/// Deterministic orthonormal basis from a unit normal (Frisvad's
/// branch on the pole, fixed arithmetic).
fn basis(n: Vec3) -> (Vec3, Vec3) {
    if n.z < -0.999_999_9 {
        return (Vec3::new(0.0, -1.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
    }
    let a = 1.0 / (1.0 + n.z);
    let b = -n.x * n.y * a;
    (
        Vec3::new(1.0 - n.x * n.x * a, b, -n.x),
        Vec3::new(b, 1.0 - n.y * n.y * a, -n.y),
    )
}

fn to_world(n: Vec3, local: [f64; 3]) -> Vec3 {
    let (t, b) = basis(n);
    Vec3::new(
        t.x * local[0] + b.x * local[1] + n.x * local[2],
        t.y * local[0] + b.y * local[1] + n.y * local[2],
        t.z * local[0] + b.z * local[1] + n.z * local[2],
    )
}

/// Cosine-weighted hemisphere sample around `n` using `det` trig (this
/// path feeds the frozen goldens; the crate-root helper uses platform
/// trig and stays for the un-hashed v0 batteries).
fn cosine_sample(n: Vec3, u1: f64, u2: f64) -> (Vec3, f64) {
    let r = u1.sqrt();
    let phi = 2.0 * PI * u2;
    let (sp, cp) = (det::sin(phi), det::cos(phi));
    let z = (1.0 - u1).max(0.0).sqrt();
    (to_world(n, [r * cp, r * sp, z]), z / PI)
}

fn ggx_d(alpha: f64, cos_m: f64) -> f64 {
    if cos_m <= 0.0 {
        return 0.0;
    }
    let a2 = alpha * alpha;
    let c2 = cos_m * cos_m;
    let t = c2 * (a2 - 1.0) + 1.0;
    a2 / (PI * t * t)
}

fn smith_g1(alpha: f64, cos_v: f64) -> f64 {
    let a2 = alpha * alpha;
    2.0 * cos_v / (cos_v + (a2 + (1.0 - a2) * cos_v * cos_v).sqrt())
}

fn schlick(f0: f64, cos_i: f64) -> f64 {
    let m = (1.0 - cos_i).clamp(0.0, 1.0);
    let m2 = m * m;
    let m5 = m2 * m2 * m; // explicit powers — never powi (hazard class)
    f0 + (1.0 - f0) * m5
}

fn bsdf_eval(mat: &Material, n: Vec3, wo: Vec3, wi: Vec3, lambda: f64) -> f64 {
    let (cos_o, cos_i) = (n.dot(wo), n.dot(wi));
    if cos_o <= 0.0 || cos_i <= 0.0 {
        return 0.0;
    }
    match mat {
        Material::Lambertian { reflectance } => reflectance.eval(lambda) / PI,
        Material::Ggx { reflectance, alpha } => {
            let hsum = Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z);
            let hn = hsum.norm();
            if hn < 1e-12 {
                return 0.0;
            }
            let m = hsum.scale(1.0 / hn);
            let d = ggx_d(*alpha, n.dot(m));
            let g = smith_g1(*alpha, cos_o) * smith_g1(*alpha, cos_i);
            let f = schlick(reflectance.eval(lambda), wo.dot(m).max(0.0));
            d * g * f / (4.0 * cos_o * cos_i)
        }
    }
}

fn bsdf_pdf(mat: &Material, n: Vec3, wo: Vec3, wi: Vec3) -> f64 {
    let cos_i = n.dot(wi);
    if cos_i <= 0.0 || n.dot(wo) <= 0.0 {
        return 0.0;
    }
    match mat {
        Material::Lambertian { .. } => cos_i / PI,
        Material::Ggx { alpha, .. } => {
            let hsum = Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z);
            let hn = hsum.norm();
            if hn < 1e-12 {
                return 0.0;
            }
            let m = hsum.scale(1.0 / hn);
            let wom = wo.dot(m);
            if wom <= 0.0 {
                return 0.0;
            }
            ggx_d(*alpha, n.dot(m)) * n.dot(m).max(0.0) / (4.0 * wom)
        }
    }
}

fn bsdf_sample(mat: &Material, n: Vec3, wo: Vec3, u1: f64, u2: f64) -> Option<(Vec3, f64)> {
    match mat {
        Material::Lambertian { .. } => {
            let (wi, pdf) = cosine_sample(n, u1, u2);
            (pdf > 0.0).then_some((wi, pdf))
        }
        Material::Ggx { alpha, .. } => {
            // Sample the half-vector from D(m)·cos m (standard GGX NDF
            // sampling; VNDF is a recorded follow-up).
            let a2 = alpha * alpha;
            let cos_m2 = ((1.0 - u1) / (u1 * (a2 - 1.0) + 1.0)).clamp(0.0, 1.0);
            let cos_m = cos_m2.sqrt();
            let sin_m = (1.0 - cos_m2).max(0.0).sqrt();
            let phi = 2.0 * PI * u2;
            let m = to_world(n, [sin_m * det::cos(phi), sin_m * det::sin(phi), cos_m]);
            let wom = wo.dot(m);
            if wom <= 0.0 {
                return None;
            }
            let wi = Vec3::new(
                2.0 * wom * m.x - wo.x,
                2.0 * wom * m.y - wo.y,
                2.0 * wom * m.z - wo.z,
            );
            if n.dot(wi) <= 0.0 {
                return None;
            }
            let pdf = ggx_d(*alpha, n.dot(m)) * n.dot(m).max(0.0) / (4.0 * wom);
            (pdf > 0.0).then_some((wi, pdf))
        }
    }
}

// ---- vector helpers (fs-geom's Vec3 has no cross) ---------------------

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.norm())
}

/// Encode a film as a linear-sRGB float EXR (channels R, G, B) —
/// byte-exact through fs-img's writer.
///
/// # Errors
/// Propagates [`fs_img::ImgError`] on shape defects.
pub fn film_to_exr(film: &Film) -> Result<Vec<u8>, fs_img::ImgError> {
    let [r, g, b] = film.to_linear_srgb();
    let ch = |name: &str, data: Vec<f32>| fs_img::Channel {
        name: name.to_string(),
        ty: fs_img::PixelType::Float,
        data,
    };
    fs_img::write_exr(
        film.width,
        film.height,
        &[ch("R", r), ch("G", g), ch("B", b)],
    )
}
