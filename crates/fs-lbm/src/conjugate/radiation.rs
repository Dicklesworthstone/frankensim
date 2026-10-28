//! Gray-diffuse radiation from exposed solid surfaces, between them and to
//! the surroundings, coupled to the conjugate energy equation.
//!
//! # Model
//!
//! Every solid voxel face that borders a fluid voxel is an emitting patch of
//! its own area with the emissivity of its material. Its escape factor `F` is
//! the fraction of its diffuse (cosine-weighted) hemisphere that leaves the
//! domain through a face declared as surroundings (an opening or inlet,
//! typically) without first hitting a solid voxel or a non-surroundings
//! domain face; it is estimated by `rays_per_face` deterministic rays
//! marched voxel by voxel (Amanatides–Woo traversal) with counter-based
//! random streams keyed by `(seed, cell, face, ray)`, so the estimate is
//! bit-reproducible and its standard error is `sqrt(F (1 - F) / rays)`.
//! The patch then loses
//!
//! ```text
//! q = eps sigma F A (T^4 - T_amb^4)
//! ```
//!
//! to the effective surroundings temperature `T_amb = <T_surroundings^4>^(1/4)`
//! over its escaping rays, preserving the sum of incident radiant powers.
//! The energy equation carries `q` as a cell sink Newton-linearized about the
//! previous iterate (`G = 4 eps sigma F A T^3`), iterated to convergence.
//!
//! # Surface-to-surface exchange
//!
//! With `RadiationConfig::surface_exchange` (the default) the emitting faces
//! also exchange with each other: they are grouped into planar tiles of up
//! to `patch_size x patch_size` faces of one material and orientation
//! (patches), and every ray now ends where it is absorbed: on an emitting
//! patch, or in the surroundings. Non-emitting solids (emissivity 0) and
//! opaque domain faces are perfect diffuse reflectors (re-radiating
//! surfaces): a ray hitting one is re-emitted cosine-weighted from the hit
//! point, up to `max_bounces` times (a ray that exhausts the budget counts
//! as returning to its own patch). The Monte Carlo exchange factors are
//! made reciprocal (`A_i F_ij = A_j F_ji`, averaged) and each row's
//! remainder becomes the patch's self-view, so the exchange conserves energy
//! exactly: the net heat all surfaces lose equals the heat the surroundings
//! receive. The irradiation `H` of each patch solves the gray-diffuse
//! radiosity system
//!
//! ```text
//! H_i = sum_j F_ij J_j + F_is sigma T_s^4,   J_j = eps_j sigma <T^4>_j + (1 - eps_j) H_j
//! ```
//!
//! by Gauss–Seidel, and each face loses `eps A (sigma T^4 - H)`, Newton-
//! linearized in its own temperature (`G = 4 eps sigma A T^3`) with `H`
//! lagged, iterated to convergence.
//!
//! # No-claim boundaries
//!
//! Without `surface_exchange`, a ray that hits a solid or a wall is simply
//! not escaping (exact when the obstructing surfaces are at the emitter's
//! temperature). With it, a patch shares one irradiation over its tile, and
//! non-emitting surfaces and opaque domain faces reflect perfectly (an
//! adiabatic re-radiating wall; declare radiating enclosure walls as solid
//! cells with an emissivity). Air is transparent; gray, diffuse, opaque
//! surfaces; all factors are Monte Carlo estimates.

use fs_exec::CancelGate;

use super::domain::{FluidProperties, SolidMaterial, Voxel, VoxelDomain};
use super::energy::{
    CellSink, EnergyConfig, EnergySolution, ThermalFace, ThermalSetup, solve_energy,
};
use super::flow::FlowField;
use super::{ChtError, finite, poll};

/// Stefan–Boltzmann constant, W/(m^2 K^4) (CODATA 2018, exact).
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;

const IRRADIATION_ITERATIONS: usize = 100_000;
const IRRADIATION_TOLERANCE: f64 = 1e-13;

/// Radiation controls.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiationConfig {
    /// Emissivity of each solid material in `solids` order (`0` disables).
    pub emissivity: Vec<f64>,
    /// Surroundings temperature seen through each domain face, K, in
    /// `Face3::ALL` order; `None` is an opaque, non-exchanging face.
    pub surroundings_k: [Option<f64>; 6],
    /// Rays per exposed face.
    pub rays_per_face: usize,
    /// Seed of the counter-based ray streams.
    pub seed: u64,
    /// Nonlinear iteration budget for steady radiation and each transient step.
    pub max_iterations: usize,
    /// Largest temperature change between Picard iterates at convergence,
    /// relative to the temperature span. Transient steps also require each
    /// cell's physical radiation power to agree with its solved tangent at
    /// this relative tolerance (with a floating-point roundoff allowance).
    pub tolerance: f64,
    /// Exchange between the emitting surfaces (see the module docs); when
    /// false, only the escape to the surroundings is modelled.
    pub surface_exchange: bool,
    /// Patch tile edge, in faces, for `surface_exchange`.
    pub patch_size: usize,
    /// Diffuse reflections a ray may take off non-emitting surfaces and
    /// opaque domain faces, for `surface_exchange`.
    pub max_bounces: usize,
}

impl RadiationConfig {
    /// Defaults for the given emissivities and surroundings.
    #[must_use]
    pub fn new(emissivity: Vec<f64>, surroundings_k: [Option<f64>; 6]) -> Self {
        Self {
            emissivity,
            surroundings_k,
            rays_per_face: 256,
            seed: 0x5EED_0FA1,
            max_iterations: 50,
            tolerance: 1e-9,
            surface_exchange: true,
            patch_size: 4,
            max_bounces: 64,
        }
    }
}

/// One emitting solid face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExposedFace {
    /// Solid cell.
    pub cell: usize,
    /// Local face of the cell (`Face3` order) bordering fluid.
    pub face: usize,
    /// Emissivity of the cell's material.
    pub emissivity: f64,
    /// Escape factor to the surroundings.
    pub escape: f64,
    /// Fourth-power mean surroundings temperature of the escaping rays, K:
    /// `<T_surroundings^4>^(1/4)`, preserving their incident radiant power.
    pub surroundings_k: f64,
}

/// Radiation evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiationReport {
    /// Emitting faces with a non-zero escape factor.
    pub exposed_faces: usize,
    /// Rays traced from the faces with a non-zero escape factor.
    pub rays: usize,
    /// Picard iterations.
    pub iterations: usize,
    /// Net heat radiated to the surroundings at the final temperatures
    /// (nonlinear law), W.
    pub radiated_w: f64,
    /// Largest temperature change of the last Picard iterate over the span.
    pub temperature_change: f64,
    /// Exchange patches (0 without `surface_exchange`).
    pub patches: usize,
}

fn splitmix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

fn unit(x: u64) -> f64 {
    // 53 random bits in [0, 1).
    (x >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

/// Ray parameter at which a ray from `origin` along `dir` leaves cell index
/// `i` across its face in the `step` direction on `axis` (non-uniform
/// Amanatides–Woo: each boundary from its own face coordinate).
fn exit_time(domain: &VoxelDomain, axis: usize, i: i64, step: i64, origin: f64, dir: f64) -> f64 {
    let face = if step > 0 { i as usize + 1 } else { i as usize };
    (domain.face_coord(axis, face) - origin) / dir
}

/// March one ray from `origin` (inside fluid cell `start`) along `dir`;
/// returns the domain-face index if it leaves through a declared
/// surroundings face before hitting a solid or an opaque face.
fn march(
    domain: &VoxelDomain,
    surroundings: &[Option<f64>; 6],
    start: [usize; 3],
    origin: [f64; 3],
    dir: [f64; 3],
) -> Option<usize> {
    let dims = domain.dims();
    let mut cell = start.map(|v| v as i64);
    let mut t_max = [f64::INFINITY; 3];
    let mut step = [0i64; 3];
    for a in 0..3 {
        if dir[a] != 0.0 {
            step[a] = if dir[a] > 0.0 { 1 } else { -1 };
            t_max[a] = exit_time(domain, a, cell[a], step[a], origin[a], dir[a]);
        }
    }
    let budget = 3 * (dims[0] + dims[1] + dims[2]) + 8;
    for _ in 0..budget {
        let axis = if t_max[0] <= t_max[1] && t_max[0] <= t_max[2] {
            0
        } else if t_max[1] <= t_max[2] {
            1
        } else {
            2
        };
        cell[axis] += step[axis];
        if cell[axis] >= 0 && cell[axis] < dims[axis] as i64 {
            t_max[axis] = exit_time(
                domain,
                axis,
                cell[axis],
                step[axis],
                origin[axis],
                dir[axis],
            );
        }
        if cell[axis] < 0 || cell[axis] >= dims[axis] as i64 {
            let side = 2 * axis + usize::from(step[axis] > 0);
            return surroundings[side].map(|_| side);
        }
        let index = domain.index(cell[0] as usize, cell[1] as usize, cell[2] as usize);
        if !domain.is_fluid(index) {
            return None;
        }
    }
    None
}

/// Admission shared by both radiation geometries.
fn admit_radiation(solids: &[SolidMaterial], config: &RadiationConfig) -> Result<(), ChtError> {
    if config.emissivity.len() != solids.len() {
        return Err(ChtError::InvalidInput {
            field: "radiation.emissivity",
            reason: format!(
                "expected {} emissivities (one per solid), got {}",
                solids.len(),
                config.emissivity.len()
            ),
        });
    }
    for &e in &config.emissivity {
        finite("radiation.emissivity", e)?;
        if !(0.0..=1.0).contains(&e) {
            return Err(ChtError::InvalidInput {
                field: "radiation.emissivity",
                reason: format!("{e} is outside [0, 1]"),
            });
        }
    }
    for t in config.surroundings_k.iter().flatten() {
        finite("radiation.surroundings_k", *t)?;
        if *t <= 0.0 {
            return Err(ChtError::InvalidInput {
                field: "radiation.surroundings_k",
                reason: "absolute temperatures must be positive".into(),
            });
        }
    }
    if config.rays_per_face == 0 {
        return Err(ChtError::InvalidInput {
            field: "radiation.rays_per_face",
            reason: "must be at least one".into(),
        });
    }
    Ok(())
}

/// Escape factors of every emitting solid face (non-zero emissivity,
/// bordering fluid).
///
/// # Errors
/// Input refusals or [`ChtError::Cancelled`].
pub fn escape_factors(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    config: &RadiationConfig,
    gate: &CancelGate,
) -> Result<Vec<ExposedFace>, ChtError> {
    admit_radiation(solids, config)?;
    let surroundings = config.surroundings_k.map(|t| t.unwrap_or(0.0));
    let mut out = Vec::new();
    for c in 0..domain.cell_count() {
        if c.is_multiple_of(1024) {
            poll(gate)?;
        }
        let Voxel::Solid(material) = domain.voxel_at(c) else {
            continue;
        };
        let emissivity = config.emissivity[usize::from(material)];
        if emissivity == 0.0 {
            continue;
        }
        let at = domain.coords(c);
        for face in 0..6 {
            let Some(n) = domain.neighbor(c, face) else {
                continue;
            };
            if !domain.is_fluid(n) {
                continue;
            }
            let axis = face / 2;
            let sign = if face % 2 == 1 { 1.0 } else { -1.0 };
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
            let plane = if face % 2 == 1 {
                domain.face_coord(axis, at[axis] + 1)
            } else {
                domain.face_coord(axis, at[axis])
            };
            let nudge = 1e-9 * domain.width(axis, at[axis]);
            let start = domain.coords(n);
            let mut escapes_by_side = [0usize; 6];
            let stream = config
                .seed
                .wrapping_add((c as u64).wrapping_mul(6).wrapping_add(face as u64) << 20);
            for ray in 0..config.rays_per_face {
                let base = splitmix(stream.wrapping_add(ray as u64));
                let (r1, r2, r3, r4) = (
                    unit(splitmix(base ^ 1)),
                    unit(splitmix(base ^ 2)),
                    unit(splitmix(base ^ 3)),
                    unit(splitmix(base ^ 4)),
                );
                // Cosine-weighted direction about the outward normal.
                let radius = r1.sqrt();
                let phi = std::f64::consts::TAU * r2;
                let mut dir = [0.0; 3];
                dir[axis] = sign * (1.0 - r1).max(0.0).sqrt();
                dir[u] = radius * phi.cos();
                dir[v] = radius * phi.sin();
                // Uniform origin on the face, nudged into the fluid cell.
                let mut origin = [0.0; 3];
                origin[axis] = sign.mul_add(nudge, plane);
                origin[u] = r3.mul_add(domain.width(u, at[u]), domain.face_coord(u, at[u]));
                origin[v] = r4.mul_add(domain.width(v, at[v]), domain.face_coord(v, at[v]));
                if let Some(side) = march(domain, &config.surroundings_k, start, origin, dir) {
                    escapes_by_side[side] += 1;
                }
            }
            let escaped = escapes_by_side.iter().sum::<usize>();
            if escaped > 0 {
                // Radiant power adds in T^4, not T. Count exits first so
                // this reduction costs at most six terms per emitting face.
                // Scale by a reached temperature to avoid overflow in T^4
                // and underflow caused by a hotter, unseen domain face.
                let scale = escapes_by_side
                    .iter()
                    .zip(surroundings)
                    .filter_map(|(&count, t)| (count > 0).then_some(t))
                    .fold(0.0, f64::max);
                let scaled_fourth_sum = escapes_by_side
                    .iter()
                    .zip(surroundings)
                    .filter(|(count, _)| **count > 0)
                    .map(|(&count, t)| count as f64 * (t / scale).powi(4))
                    .sum::<f64>();
                let effective_k = scale * (scaled_fourth_sum / escaped as f64).sqrt().sqrt();
                out.push(ExposedFace {
                    cell: c,
                    face,
                    emissivity,
                    escape: escaped as f64 / config.rays_per_face as f64,
                    surroundings_k: effective_k,
                });
            }
        }
    }
    Ok(out)
}

/// Radiative sinks of `faces` linearized about `temperature` by Newton:
/// `q(T) ~ q(T_k) + q'(T_k) (T - T_k)` with `q' = 4 eps sigma F A T_k^3`,
/// written as the sink `G (T - T_s)` with `G = q'` and
/// `T_s = T_k - q(T_k) / q'` (exact at `T = T_k`, quadratic convergence).
#[must_use]
pub fn radiative_sinks(
    domain: &VoxelDomain,
    faces: &[ExposedFace],
    temperature: &[f64],
) -> Vec<CellSink> {
    faces
        .iter()
        .map(|f| {
            let (t, ta) = (temperature[f.cell], f.surroundings_k);
            let area = domain.face_area(f.cell, f.face / 2);
            let scale = f.emissivity * STEFAN_BOLTZMANN * f.escape * area;
            let slope = 4.0 * scale * t * t * t;
            let q = scale * (t.powi(4) - ta.powi(4));
            CellSink {
                cell: f.cell,
                conductance_w_k: slope,
                temperature_k: t - q / slope.max(f64::MIN_POSITIVE),
            }
        })
        .collect()
}

/// Net radiated power at `temperature` by the nonlinear law, W.
#[must_use]
pub fn radiated_power(domain: &VoxelDomain, faces: &[ExposedFace], temperature: &[f64]) -> f64 {
    faces
        .iter()
        .map(|f| {
            let (t, ta) = (temperature[f.cell], f.surroundings_k);
            let area = domain.face_area(f.cell, f.face / 2);
            f.emissivity * STEFAN_BOLTZMANN * f.escape * area * (t.powi(4) - ta.powi(4))
        })
        .sum()
}

/// Where a traced ray ends.
enum Hit {
    /// Left through a surroundings face at this temperature.
    Surroundings(f64),
    /// Entered solid `cell` through its local `face` at `point`, from the
    /// fluid cell `from`.
    Solid {
        cell: usize,
        face: usize,
        point: [f64; 3],
        from: [usize; 3],
    },
    /// Reached opaque domain face `side` at `point` from fluid cell `from`.
    Opaque {
        side: usize,
        point: [f64; 3],
        from: [usize; 3],
    },
    /// Traversal budget exhausted.
    Lost,
}

fn trace(
    domain: &VoxelDomain,
    surroundings: &[Option<f64>; 6],
    start: [usize; 3],
    origin: [f64; 3],
    dir: [f64; 3],
) -> Hit {
    let dims = domain.dims();
    let mut cell = start.map(|v| v as i64);
    let mut t_max = [f64::INFINITY; 3];
    let mut step = [0i64; 3];
    for a in 0..3 {
        if dir[a] != 0.0 {
            step[a] = if dir[a] > 0.0 { 1 } else { -1 };
            t_max[a] = exit_time(domain, a, cell[a], step[a], origin[a], dir[a]);
        }
    }
    let budget = 3 * (dims[0] + dims[1] + dims[2]) + 8;
    for _ in 0..budget {
        let axis = if t_max[0] <= t_max[1] && t_max[0] <= t_max[2] {
            0
        } else if t_max[1] <= t_max[2] {
            1
        } else {
            2
        };
        let from = cell.map(|v| v as usize);
        let t = t_max[axis];
        let point = [0, 1, 2].map(|a| t.mul_add(dir[a], origin[a]));
        cell[axis] += step[axis];
        if cell[axis] >= 0 && cell[axis] < dims[axis] as i64 {
            t_max[axis] = exit_time(
                domain,
                axis,
                cell[axis],
                step[axis],
                origin[axis],
                dir[axis],
            );
        }
        if cell[axis] < 0 || cell[axis] >= dims[axis] as i64 {
            let side = 2 * axis + usize::from(step[axis] > 0);
            return match surroundings[side] {
                Some(t) => Hit::Surroundings(t),
                None => Hit::Opaque { side, point, from },
            };
        }
        let index = domain.index(cell[0] as usize, cell[1] as usize, cell[2] as usize);
        if !domain.is_fluid(index) {
            // Entered through the solid's face facing the ray's origin.
            let face = 2 * axis + usize::from(step[axis] < 0);
            return Hit::Solid {
                cell: index,
                face,
                point,
                from,
            };
        }
    }
    Hit::Lost
}

/// A cosine-weighted direction about the unit normal `sign * e_axis`.
fn cosine_direction(axis: usize, sign: f64, r1: f64, r2: f64) -> [f64; 3] {
    let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
    let radius = r1.sqrt();
    let phi = std::f64::consts::TAU * r2;
    let mut dir = [0.0; 3];
    dir[axis] = sign * (1.0 - r1).max(0.0).sqrt();
    dir[u] = radius * phi.cos();
    dir[v] = radius * phi.sin();
    dir
}

/// Surface-to-surface radiation geometry: emitting faces grouped into
/// patches with their Monte Carlo exchange and escape fractions (see the
/// module docs).
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceExchange {
    /// Emitting faces `(cell, local face, patch)`.
    faces: Vec<(usize, usize, usize)>,
    /// Per patch: emissivity.
    emissivity: Vec<f64>,
    /// Per patch: area, m^2.
    patch_area: Vec<f64>,
    /// Per patch: reciprocal exchange fractions `F_ij` to other patches.
    exchange: Vec<Vec<(usize, f64)>>,
    /// Per patch: self-view fraction `F_ii`.
    self_view: Vec<f64>,
    /// Per patch: escaping fraction `F_is` and `F_is <sigma T_s^4>`, W/m^2.
    escape: Vec<(f64, f64)>,
    /// Area of each emitting face (as `faces`), m^2.
    face_areas: Vec<f64>,
    /// Rays traced.
    rays: usize,
}

impl SurfaceExchange {
    /// Trace the exchange geometry of every emitting face.
    ///
    /// # Errors
    /// Input refusals (as [`escape_factors`], plus `patch_size` and
    /// `max_bounces` of at least one) or [`ChtError::Cancelled`].
    #[allow(clippy::too_many_lines)] // patching, tracing, reciprocity
    pub fn build(
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        config: &RadiationConfig,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        admit_radiation(solids, config)?;
        if config.patch_size == 0 || config.max_bounces == 0 {
            return Err(ChtError::InvalidInput {
                field: "radiation.patch_size",
                reason: "patch_size and max_bounces must be at least one".into(),
            });
        }
        let nudge = 1e-9 * domain.min_width();
        // Emitting faces, keyed by material, orientation, plane and tile.
        let mut keyed: Vec<([usize; 5], usize, usize, f64)> = Vec::new();
        for c in 0..domain.cell_count() {
            if c.is_multiple_of(4096) {
                poll(gate)?;
            }
            let Voxel::Solid(material) = domain.voxel_at(c) else {
                continue;
            };
            let emissivity = config.emissivity[usize::from(material)];
            if emissivity == 0.0 {
                continue;
            }
            let at = domain.coords(c);
            for face in 0..6 {
                let Some(n) = domain.neighbor(c, face) else {
                    continue;
                };
                if !domain.is_fluid(n) {
                    continue;
                }
                let axis = face / 2;
                let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
                let key = [
                    usize::from(material),
                    face,
                    at[axis],
                    at[u] / config.patch_size,
                    at[v] / config.patch_size,
                ];
                keyed.push((key, c, face, emissivity));
            }
        }
        keyed.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        let mut faces = Vec::with_capacity(keyed.len());
        let mut face_areas = Vec::with_capacity(keyed.len());
        let (mut emissivity, mut patch_area) = (Vec::new(), Vec::new());
        for (i, &(key, cell, face, eps)) in keyed.iter().enumerate() {
            if i == 0 || keyed[i - 1].0 != key {
                emissivity.push(eps);
                patch_area.push(0.0f64);
            }
            let patch = emissivity.len() - 1;
            let area = domain.face_area(cell, face / 2);
            patch_area[patch] += area;
            face_areas.push(area);
            faces.push((cell, face, patch));
        }
        let patches = emissivity.len();
        let mut lookup: Vec<(usize, usize)> = faces
            .iter()
            .map(|&(cell, face, patch)| (cell * 6 + face, patch))
            .collect();
        lookup.sort_unstable();
        let patch_of = |cell: usize, face: usize| {
            lookup
                .binary_search_by_key(&(cell * 6 + face), |&(k, _)| k)
                .ok()
                .map(|i| lookup[i].1)
        };
        // Ray counts per patch: absorbed by other patches, escaped (with
        // sigma T^4), and returned to the emitter.
        let mut hits: Vec<std::collections::BTreeMap<usize, usize>> =
            vec![std::collections::BTreeMap::new(); patches];
        let mut escaped = vec![(0usize, 0.0f64); patches];
        let mut traced = vec![0usize; patches];
        for (index, &(c, face, patch)) in faces.iter().enumerate() {
            if index.is_multiple_of(256) {
                poll(gate)?;
            }
            let at = domain.coords(c);
            let axis = face / 2;
            let sign = if face % 2 == 1 { 1.0 } else { -1.0 };
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
            let plane = if face % 2 == 1 {
                domain.face_coord(axis, at[axis] + 1)
            } else {
                domain.face_coord(axis, at[axis])
            };
            let first = domain.coords(domain.neighbor(c, face).expect("borders fluid"));
            let stream = config
                .seed
                .wrapping_add((c as u64).wrapping_mul(6).wrapping_add(face as u64) << 20);
            for ray in 0..config.rays_per_face {
                traced[patch] += 1;
                let base = splitmix(stream.wrapping_add(ray as u64));
                let draw = |k: u64| unit(splitmix(base ^ k));
                let mut dir = cosine_direction(axis, sign, draw(1), draw(2));
                let mut origin = [0.0; 3];
                origin[axis] = sign.mul_add(nudge, plane);
                origin[u] = draw(3).mul_add(domain.width(u, at[u]), domain.face_coord(u, at[u]));
                origin[v] = draw(4).mul_add(domain.width(v, at[v]), domain.face_coord(v, at[v]));
                let mut start = first;
                for bounce in 0..config.max_bounces {
                    // A reflector re-emits about its inward normal.
                    let reflection = match trace(domain, &config.surroundings_k, start, origin, dir)
                    {
                        Hit::Surroundings(t) => {
                            escaped[patch].0 += 1;
                            escaped[patch].1 += STEFAN_BOLTZMANN * t.powi(4);
                            None
                        }
                        Hit::Solid {
                            cell,
                            face: entered,
                            point,
                            from,
                        } => {
                            if let Some(target) = patch_of(cell, entered) {
                                *hits[patch].entry(target).or_insert(0) += 1;
                                None
                            } else {
                                // Non-emitting solid: the outward normal of
                                // the entered face points back into `from`.
                                let normal = if entered % 2 == 1 { 1.0 } else { -1.0 };
                                Some((entered / 2, normal, point, from))
                            }
                        }
                        Hit::Opaque { side, point, from } => {
                            let normal = if side % 2 == 1 { -1.0 } else { 1.0 };
                            Some((side / 2, normal, point, from))
                        }
                        Hit::Lost => None,
                    };
                    let Some((normal_axis, normal, point, from)) = reflection else {
                        break;
                    };
                    let k = 8 + 2 * bounce as u64;
                    dir = cosine_direction(normal_axis, normal, draw(k), draw(k + 1));
                    origin = point;
                    origin[normal_axis] = normal.mul_add(nudge, point[normal_axis]);
                    start = from;
                }
                // Lost rays and exhausted bounces return to the emitter (no
                // net exchange): they stay in the self-view remainder.
            }
        }
        // Reciprocity AND exact row sums: average the raw exchange A_i F_ij
        // with A_j F_ji into a symmetric S (self-hits on the diagonal), then
        // balance it symmetrically, G = D S D, so every row sums to the
        // emitter's absorbed share A_i (1 - F_is) (symmetric Sinkhorn
        // iteration d <- sqrt(d r / (S d))). Averaging alone breaks the row
        // sums, and clipping the resulting negative self-views biased a
        // black-enclosure test by 1 % in T^4.
        let area = patch_area.clone();
        let raw = |i: usize, j: usize| {
            hits[i]
                .get(&j)
                .map_or(0.0, |&n| n as f64 / traced[i] as f64)
        };
        let escape: Vec<(f64, f64)> = (0..patches)
            .map(|i| {
                let n = traced[i] as f64;
                (escaped[i].0 as f64 / n, escaped[i].1 / n)
            })
            .collect();
        // Self-hits: rays returning to their own patch, lost rays and
        // exhausted bounce budgets.
        let returned: Vec<f64> = (0..patches)
            .map(|i| {
                let out: f64 = hits[i]
                    .iter()
                    .filter(|&(&j, _)| j != i)
                    .map(|(_, &n)| n as f64)
                    .sum();
                (1.0 - out / traced[i] as f64 - escape[i].0).max(0.0)
            })
            .collect();
        let mut symmetric: Vec<Vec<(usize, f64)>> = vec![Vec::new(); patches];
        for i in 0..patches {
            symmetric[i].push((i, area[i] * returned[i]));
            for &j in hits[i].keys() {
                if j == i {
                    continue;
                }
                let shared = 0.5 * area[i].mul_add(raw(i, j), area[j] * raw(j, i));
                symmetric[i].push((j, shared));
                if !hits[j].contains_key(&i) {
                    symmetric[j].push((i, shared));
                }
            }
        }
        for row in &mut symmetric {
            row.sort_by_key(|&(j, _)| j);
        }
        let target: Vec<f64> = (0..patches)
            .map(|i| area[i] * (1.0 - escape[i].0))
            .collect();
        let mut d = vec![1.0f64; patches];
        for sweep in 0..10_000 {
            if sweep % 64 == 0 {
                poll(gate)?;
            }
            let mut worst = 0.0f64;
            for i in 0..patches {
                let sum: f64 = symmetric[i].iter().map(|&(j, g)| g * d[j]).sum();
                if sum > 0.0 && target[i] > 0.0 {
                    worst = worst.max((d[i] * sum / target[i] - 1.0).abs());
                    d[i] = (d[i] * target[i] / sum).sqrt();
                }
            }
            if worst <= 1e-13 {
                break;
            }
        }
        let mut exchange = vec![Vec::new(); patches];
        let mut self_view = vec![0.0f64; patches];
        for i in 0..patches {
            for &(j, g) in &symmetric[i] {
                let f = d[i] * g * d[j] / area[i];
                if j == i {
                    self_view[i] = f;
                } else if f > 0.0 {
                    exchange[i].push((j, f));
                }
            }
        }
        Ok(Self {
            faces,
            emissivity,
            patch_area,
            exchange,
            self_view,
            escape,
            face_areas,
            rays: traced.iter().sum(),
        })
    }

    /// Patches.
    #[must_use]
    pub fn patches(&self) -> usize {
        self.emissivity.len()
    }

    /// Rays traced.
    #[must_use]
    pub fn rays(&self) -> usize {
        self.rays
    }

    /// Area-mean `sigma T^4` of each patch's faces, W/m^2.
    fn blackbody(&self, temperature: &[f64]) -> Vec<f64> {
        let mut sum = vec![0.0; self.patches()];
        for (&(cell, _, patch), area) in self.faces.iter().zip(&self.face_areas) {
            sum[patch] += area * STEFAN_BOLTZMANN * temperature[cell].powi(4);
        }
        sum.iter()
            .zip(&self.patch_area)
            .map(|(s, a)| s / a)
            .collect()
    }

    /// Irradiation of every patch at `temperature`, W/m^2 (Gauss–Seidel on
    /// the radiosity system to 1e-13 of the largest value).
    #[must_use]
    pub fn irradiation(&self, temperature: &[f64]) -> Vec<f64> {
        // Retain the public sampling API and its original change criterion.
        // Physical solvers use the fallible checked path below.
        match self.irradiation_iterations(
            temperature,
            || Ok::<(), std::convert::Infallible>(()),
            |_, _| Ok(true),
        ) {
            Ok((irradiation, _, _)) => irradiation,
            Err(never) => match never {},
        }
    }

    fn irradiation_iterations<E>(
        &self,
        temperature: &[f64],
        mut checkpoint: impl FnMut() -> Result<(), E>,
        mut accept: impl FnMut(&[f64], &[f64]) -> Result<bool, E>,
    ) -> Result<(Vec<f64>, Vec<f64>, bool), E> {
        checkpoint()?;
        let e = self.blackbody(temperature);
        let n = self.patches();
        let mut h: Vec<f64> = (0..n).map(|i| e[i]).collect();
        for _ in 0..IRRADIATION_ITERATIONS {
            checkpoint()?;
            let mut change = 0.0f64;
            let mut scale = 0.0f64;
            for i in 0..n {
                let mut incoming = self.escape[i].1;
                for &(j, f) in &self.exchange[i] {
                    incoming = f.mul_add(
                        (1.0 - self.emissivity[j]).mul_add(h[j], self.emissivity[j] * e[j]),
                        incoming,
                    );
                }
                // The self-view solved in place.
                let (fs, eps) = (self.self_view[i], self.emissivity[i]);
                let next = fs.mul_add(eps * e[i], incoming) / (1.0 - fs * (1.0 - eps));
                change = change.max((next - h[i]).abs());
                scale = scale.max(next.abs());
                h[i] = next;
            }
            if change <= IRRADIATION_TOLERANCE * scale && accept(&e, &h)? {
                return Ok((h, e, true));
            }
        }
        Ok((h, e, false))
    }

    fn checked_irradiation(
        &self, temperature: &[f64], gate: &CancelGate,
    ) -> Result<Vec<f64>, ChtError> {
        let (h, e, converged) = self.irradiation_iterations(
            temperature,
            || poll(gate),
            |e, h| Ok(self.irradiation_residual(e, h)? <= IRRADIATION_TOLERANCE),
        )?;
        if !converged {
            return Err(ChtError::SolverNotConverged {
                system: "radiation irradiation",
                iterations: IRRADIATION_ITERATIONS,
                relative_residual: self.irradiation_residual(&e, &h)?,
                tolerance: IRRADIATION_TOLERANCE,
            });
        }
        poll(gate)?;
        Ok(h)
    }

    // Reconstruct the ORIGINAL irradiation equations and the independent
    // face-emission/escape energy account. A small iterate change alone can
    // accept a badly wrong heat transfer when emissivity is close to zero.
    fn irradiation_residual(&self, e: &[f64], h: &[f64]) -> Result<f64, ChtError> {
        let mut radiosity = Vec::with_capacity(self.patches());
        for (i, (&emission, &incoming)) in e.iter().zip(h).enumerate() {
            finite("radiation.blackbody_w_m2", emission)?;
            finite("radiation.irradiation_w_m2", incoming)?;
            let value = (1.0 - self.emissivity[i]).mul_add(
                incoming, self.emissivity[i] * emission,
            );
            finite("radiation.radiosity_w_m2", value)?;
            if emission < 0.0 || incoming < 0.0 || value < 0.0 {
                return Err(ChtError::InvalidInput {
                    field: "radiation.irradiation",
                    reason: "radiative intensities must be non-negative".into(),
                });
            }
            radiosity.push(value);
        }
        let mut worst = 0.0_f64;
        let (mut net, mut external, mut heat_scale, mut roundoff_scale) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..self.patches() {
            let mut incoming = self.escape[i].1;
            for &(j, factor) in &self.exchange[i] {
                incoming = factor.mul_add(radiosity[j], incoming);
            }
            incoming = self.self_view[i].mul_add(radiosity[i], incoming);
            finite("radiation.irradiation_equation_w_m2", incoming)?;
            let scale = incoming.abs().max(h[i].abs()).max(f64::MIN_POSITIVE);
            worst = worst.max(
                ((incoming - h[i]).abs() - 64.0 * f64::EPSILON * scale).max(0.0) / scale,
            );
            let eps_area = self.emissivity[i] * self.patch_area[i];
            let emitted = eps_area * (e[i] - h[i]);
            let escaped = self.patch_area[i]
                * self.escape[i].0.mul_add(radiosity[i], -self.escape[i].1);
            finite("radiation.patch_emission_w", emitted)?;
            finite("radiation.patch_escape_w", escaped)?;
            net += emitted;
            external += escaped;
            heat_scale += emitted.abs() + escaped.abs();
            roundoff_scale += eps_area * (e[i].abs() + h[i].abs())
                + self.patch_area[i]
                    * ((self.escape[i].0 * radiosity[i]).abs() + self.escape[i].1.abs());
        }
        let closure = net - external;
        finite("radiation.enclosure_closure_w", closure)?;
        finite("radiation.enclosure_heat_scale_w", heat_scale)?;
        finite("radiation.enclosure_roundoff_scale_w", roundoff_scale)?;
        Ok(worst.max(
            (closure.abs() - 64.0 * f64::EPSILON * roundoff_scale).max(0.0)
                / heat_scale.max(f64::MIN_POSITIVE),
        ))
    }

    /// Radiative sinks at `temperature` with the patch irradiation
    /// `irradiation`: each face loses `eps A (sigma T^4 - H)`, Newton-
    /// linearized in its own temperature.
    #[must_use]
    pub fn sinks(&self, temperature: &[f64], irradiation: &[f64]) -> Vec<CellSink> {
        self.faces
            .iter()
            .zip(&self.face_areas)
            .map(|(&(cell, _, patch), face_area)| {
                let t = temperature[cell];
                let eps_a = self.emissivity[patch] * face_area;
                let slope = 4.0 * eps_a * STEFAN_BOLTZMANN * t * t * t;
                let q = eps_a * STEFAN_BOLTZMANN.mul_add(t.powi(4), -irradiation[patch]);
                CellSink {
                    cell,
                    conductance_w_k: slope,
                    temperature_k: t - q / slope.max(f64::MIN_POSITIVE),
                }
            })
            .collect()
    }

    /// Net heat the surroundings receive at `temperature` and `irradiation`,
    /// W (equal to the net heat all patches lose: the exchange is
    /// reciprocal).
    #[must_use]
    pub fn radiated_to_surroundings(&self, temperature: &[f64], irradiation: &[f64]) -> f64 {
        let e = self.blackbody(temperature);
        (0..self.patches())
            .map(|i| {
                let eps = self.emissivity[i];
                let radiosity = (1.0 - eps).mul_add(irradiation[i], eps * e[i]);
                let area = self.patch_area[i];
                area * self.escape[i].0.mul_add(radiosity, -self.escape[i].1)
            })
            .sum()
    }

    /// Net heat all emitting faces lose at `temperature` and `irradiation`,
    /// W, by the face law.
    #[must_use]
    pub fn net_emission(&self, temperature: &[f64], irradiation: &[f64]) -> f64 {
        self.faces
            .iter()
            .zip(&self.face_areas)
            .map(|(&(cell, _, patch), area)| {
                self.emissivity[patch]
                    * area
                    * STEFAN_BOLTZMANN.mul_add(temperature[cell].powi(4), -irradiation[patch])
            })
            .sum()
    }
}

/// The radiation geometry in force: escape factors only, or the
/// surface-to-surface exchange.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Radiators {
    Escape(Vec<ExposedFace>),
    Exchange(SurfaceExchange),
}

/// Radiation evaluated once at one admissible temperature field. The same
/// irradiation supplies the next tangent, the physical cell powers, and
/// the heat actually escaping to the surroundings.
pub(crate) struct RadiationEndpoint {
    pub(crate) sinks: Vec<CellSink>,
    pub(crate) cell_power_w: Vec<f64>,
    pub(crate) cell_scale_w: Vec<f64>,
    pub(crate) radiated_w: f64,
}

impl Radiators {
    pub(crate) fn build(
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        config: &RadiationConfig,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        Ok(if config.surface_exchange {
            Self::Exchange(SurfaceExchange::build(domain, solids, config, gate)?)
        } else {
            Self::Escape(escape_factors(domain, solids, config, gate)?)
        })
    }

    /// Evaluate the nonlinear face law and its next Newton tangent together.
    pub(crate) fn endpoint(
        &self,
        domain: &VoxelDomain,
        temperature: &[f64],
        gate: &CancelGate,
    ) -> Result<RadiationEndpoint, ChtError> {
        if temperature.len() != domain.cell_count() {
            return Err(ChtError::InvalidInput {
                field: "radiation.temperature",
                reason: "expected one absolute temperature per cell".into(),
            });
        }
        for (cell, &t) in temperature.iter().enumerate() {
            if cell % 1024 == 0 {
                poll(gate)?;
            }
            super::finite_positive("radiation.temperature", t)?;
        }
        let mut cell_power_w = vec![0.0; domain.cell_count()];
        let mut cell_scale_w = vec![0.0; domain.cell_count()];
        let (sinks, radiated_w) = match self {
            Self::Escape(faces) => {
                for face in faces {
                    let scale = face.emissivity
                        * STEFAN_BOLTZMANN
                        * face.escape
                        * domain.face_area(face.cell, face.face / 2);
                    let emission = scale * temperature[face.cell].powi(4);
                    let absorption = scale * face.surroundings_k.powi(4);
                    cell_power_w[face.cell] += emission - absorption;
                    cell_scale_w[face.cell] += emission.abs() + absorption.abs();
                }
                (
                    radiative_sinks(domain, faces, temperature),
                    cell_power_w.iter().sum(),
                )
            }
            Self::Exchange(model) => {
                let irradiation = model.checked_irradiation(temperature, gate)?;
                for (&(cell, _, patch), &area) in model.faces.iter().zip(&model.face_areas) {
                    let eps_a = model.emissivity[patch] * area;
                    let emission = STEFAN_BOLTZMANN * temperature[cell].powi(4);
                    cell_power_w[cell] +=
                        eps_a * STEFAN_BOLTZMANN.mul_add(temperature[cell].powi(4), -irradiation[patch]);
                    cell_scale_w[cell] += eps_a * (emission.abs() + irradiation[patch].abs());
                }
                (
                    model.sinks(temperature, &irradiation),
                    model.radiated_to_surroundings(temperature, &irradiation),
                )
            }
        };
        for (&power, &scale) in cell_power_w.iter().zip(&cell_scale_w) {
            finite("radiation.power_w", power)?;
            finite("radiation.power_scale_w", scale)?;
        }
        finite("radiation.radiated_w", radiated_w)?;
        for sink in &sinks {
            finite("radiation.conductance_w_k", sink.conductance_w_k)?;
            finite("radiation.sink_temperature_k", sink.temperature_k)?;
        }
        Ok(RadiationEndpoint {
            sinks,
            cell_power_w,
            cell_scale_w,
            radiated_w,
        })
    }


    /// A reference temperature for the first linearization: the mean
    /// surroundings of escaping faces, else `fallback`.
    fn reference(&self, fallback: f64) -> f64 {
        match self {
            Self::Escape(faces) if !faces.is_empty() => {
                faces.iter().map(|f| f.surroundings_k).sum::<f64>() / faces.len() as f64
            }
            _ => fallback,
        }
    }

    fn counts(&self, rays_per_face: usize) -> (usize, usize, usize) {
        match self {
            Self::Escape(faces) => (faces.len(), faces.len() * rays_per_face, 0),
            Self::Exchange(model) => (model.faces.len(), model.rays(), model.patches()),
        }
    }
}

/// Steady conjugate energy with surface radiation to the surroundings on a
/// given flow, by Picard iteration of the linearized radiative sinks.
///
/// # Errors
/// Input and solver refusals, [`ChtError::SolverNotConverged`] (system
/// `"radiation"`) when the Picard budget ends first, or
/// [`ChtError::Cancelled`].
#[allow(clippy::too_many_arguments)] // physics inputs + both configurations
pub fn solve_energy_radiating(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    energy: &EnergyConfig,
    radiation: &RadiationConfig,
    gate: &CancelGate,
) -> Result<(EnergySolution, RadiationReport), ChtError> {
    super::finite_positive("radiation.tolerance", radiation.tolerance)?;
    if radiation.max_iterations == 0 {
        return Err(ChtError::InvalidInput {
            field: "radiation.max_iterations", reason: "must be at least one".into(),
        });
    }
    let radiators = Radiators::build(domain, solids, radiation, gate)?;
    // Start from the sinks linearized about a reference temperature (the
    // surroundings, else a declared face or fixed temperature): without
    // radiation the first iterate can lack the dominant heat path.
    let reference = radiators.reference(reference_temperature(setup, radiation));
    let mut first = setup.clone();
    first.cell_sinks.extend(
        radiators.endpoint(domain, &vec![reference; domain.cell_count()], gate)?.sinks,
    );
    let mut solution = solve_energy(domain, fluid, solids, flow, &first, energy, gate)?;
    let mut endpoint = radiators.endpoint(domain, &solution.temperature, gate)?;
    let mut change = f64::INFINITY;
    let mut iterations = 0usize;
    while iterations < radiation.max_iterations {
        poll(gate)?;
        iterations += 1;
        let mut step = setup.clone();
        step.cell_sinks.extend_from_slice(&endpoint.sinks);
        let next = solve_energy(domain, fluid, solids, flow, &step, energy, gate)?;
        let next_endpoint = radiators.endpoint(domain, &next.temperature, gate)?;
        let (lo, hi) = next
            .temperature
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), t| {
                (lo.min(*t), hi.max(*t))
            });
        let largest = next
            .temperature
            .iter()
            .zip(&solution.temperature)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        change = largest / (hi - lo).max(f64::MIN_POSITIVE);
        solution = next;
        endpoint = next_endpoint;
        if change <= radiation.tolerance {
            break;
        }
    }
    if change > radiation.tolerance {
        return Err(ChtError::SolverNotConverged {
            system: "radiation",
            iterations,
            relative_residual: change,
            tolerance: radiation.tolerance,
        });
    }
    let (exposed_faces, rays, patches) = radiators.counts(radiation.rays_per_face);
    let report = RadiationReport {
        exposed_faces,
        rays,
        iterations,
        radiated_w: endpoint.radiated_w,
        temperature_change: change,
        patches,
    };
    Ok((solution, report))
}

/// The first linearization's temperature: the mean declared surroundings,
/// else the first face or fixed-cell temperature, else 300 K.
pub(crate) fn reference_temperature(setup: &ThermalSetup, radiation: &RadiationConfig) -> f64 {
    let declared: Vec<f64> = radiation.surroundings_k.iter().flatten().copied().collect();
    if !declared.is_empty() {
        return declared.iter().sum::<f64>() / declared.len() as f64;
    }
    setup
        .faces
        .iter()
        .find_map(|rule| match *rule {
            ThermalFace::Temperature(t) | ThermalFace::Inflow { temperature: t } => Some(t),
            ThermalFace::Convective { ambient, .. } => Some(ambient),
            _ => None,
        })
        .or_else(|| setup.fixed_temperature.first().map(|&(_, t)| t))
        .unwrap_or(300.0)
}

#[cfg(test)]
mod checked_irradiation_tests {
    use super::*;

    // Exact equal-area parallel plates isolate the thermal radiosity solve
    // from Monte Carlo geometry error. Every emitted ray reaches its neighbor.
    fn closed_plates(emissivity: f64) -> SurfaceExchange {
        SurfaceExchange {
            faces: vec![(0, 0, 0), (1, 1, 1)],
            emissivity: vec![emissivity; 2],
            patch_area: vec![1.0; 2],
            exchange: vec![vec![(1, 1.0)], vec![(0, 1.0)]],
            self_view: vec![0.0; 2],
            escape: vec![(0.0, 0.0); 2],
            face_areas: vec![1.0; 2],
            rays: 0,
        }
    }

    #[test]
    fn checked_steady_consumers_refuse_unresolved_closed_cavity() {
        use crate::conjugate::{
            FvBoundary, FvBuoyancyConfig, SimpleConfig, fv_natural_convection,
        };
        let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |p| {
            if p[0] > 1.0 && p[0] < 2.0 { Voxel::Fluid } else { Voxel::Solid(0) }
        }).unwrap();
        let fluid = FluidProperties {
            density_kg_m3: 1.0, specific_heat_j_kg_k: 1.0,
            conductivity_w_m_k: 1.0, kinematic_viscosity_m2_s: 1.0,
        };
        let solids = [SolidMaterial::new("reflective plate", 1.0)];
        let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
        setup.fixed_temperature = vec![(0, 400.0), (1, 350.0), (2, 300.0)];
        let mut radiation = RadiationConfig::new(vec![1e-9], [None; 6]);
        radiation.rays_per_face = 64;
        radiation.patch_size = 1;
        let gate = CancelGate::new();
        assert!(matches!(
            solve_energy_radiating(&domain, &fluid, &solids, &FlowField::quiescent(&domain),
                &setup, &EnergyConfig::default(), &radiation, &gate),
            Err(ChtError::SolverNotConverged { system: "radiation irradiation", .. })
        ));
        let mut natural = FvBuoyancyConfig::new([0.0; 3], 0.0, 350.0,
            SimpleConfig::new([FvBoundary::wall(); 6]));
        natural.max_couplings = 4;
        natural.sweeps_per_coupling = 1;
        natural.radiation = Some(radiation);
        assert!(matches!(
            fv_natural_convection(&domain, &fluid, &solids, &setup, &natural, &gate),
            Err(ChtError::SolverNotConverged { system: "radiation irradiation", .. })
        ));
        natural.radiation.as_mut().unwrap().emissivity[0] = 0.8;
        let (steady, report) = solve_energy_radiating(
            &domain, &fluid, &solids, &FlowField::quiescent(&domain),
            &setup, &EnergyConfig::default(), natural.radiation.as_ref().unwrap(), &gate,
        ).unwrap();
        assert_eq!(steady.temperature, vec![400.0, 350.0, 300.0]);
        assert_eq!(report.radiated_w, 0.0);
        let coupled = fv_natural_convection(
            &domain, &fluid, &solids, &setup, &natural, &gate,
        ).unwrap();
        assert_eq!(coupled.energy.temperature, vec![400.0, 350.0, 300.0]);
        assert_eq!(coupled.report.radiated_w, Some(0.0));
        for balance in [steady.report.balance, coupled.energy.report.balance] {
            for value in [
                balance.source_w, balance.fixed_cell_injection_w, balance.boundary_outflow_w,
                balance.advective_outflow_w, balance.sink_outflow_w, balance.residual_w,
                balance.relative_residual,
            ] {
                assert!(value.is_finite(), "{balance:?}");
            }
        }
    }

    #[test]
    fn checked_closed_cavity_matches_gray_parallel_plate_heat() {
        let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
        let gate = CancelGate::new();
        for emissivity in [1.0, 0.8, 0.1, 0.01] {
            let endpoint = Radiators::Exchange(closed_plates(emissivity))
                .endpoint(&domain, &[400.0, 300.0], &gate).unwrap();
            let expected = STEFAN_BOLTZMANN * (400.0_f64.powi(4) - 300.0_f64.powi(4))
                / (2.0 / emissivity - 1.0);
            assert!((endpoint.cell_power_w[0] - expected).abs() < 1e-11 * expected);
            assert!((endpoint.cell_power_w[1] + expected).abs() < 1e-11 * expected);
            assert_eq!(endpoint.radiated_w, 0.0);
            let closure = endpoint.cell_power_w.iter().sum::<f64>();
            assert!(closure.abs() < 1e-11 * expected, "eps={emissivity}: {closure}");
        }
    }

    #[test]
    fn checked_closed_cavity_refuses_unresolved_reflection_and_observes_cancellation() {
        let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
        let gate = CancelGate::new();
        // 1e-9 exhausts the sweep budget; 1e-15 tricks the legacy iterate-change
        // threshold almost immediately. Neither may publish fictitious net heat.
        for emissivity in [1e-9, 1e-15] {
            let model = Radiators::Exchange(closed_plates(emissivity));
            let equilibrium = model.endpoint(&domain, &[350.0; 2], &gate).unwrap();
            assert!(equilibrium.cell_power_w.iter().all(|q| q.abs() < 1e-20));
            assert!(matches!(
                model.endpoint(&domain, &[400.0, 300.0], &gate),
                Err(ChtError::SolverNotConverged {
                    system: "radiation irradiation", iterations: IRRADIATION_ITERATIONS, ..
                })
            ));
        }
        let cancelled = CancelGate::new();
        cancelled.request();
        assert!(matches!(
            closed_plates(0.8).checked_irradiation(&[400.0, 300.0], &cancelled),
            Err(ChtError::Cancelled)
        ));
    }
}
