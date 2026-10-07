//! LVEL algebraic turbulence model (Agonafer, Gan-Li & Spalding, ASME HTD
//! 324, 1996) for the SIMPLEC voxel flow: the electronics-cooling model of
//! choice where fan-driven channel flows are transitional or turbulent and
//! meshes are too coarse for two-equation models.
//!
//! # Model
//!
//! At each fluid cell, from the wall distance `L` and the local speed `V`,
//! the local Reynolds number `Re_L = V L / nu` fixes the law-of-the-wall
//! coordinates through Spalding's unified profile
//!
//! ```text
//! y+ = u+ + (1/E) [exp(k u+) - 1 - k u+ - (k u+)^2 / 2 - (k u+)^3 / 6],
//! u+ y+ = Re_L,
//! ```
//!
//! solved for `u+` by Newton, and the effective viscosity is its slope
//!
//! ```text
//! nu_eff / nu = dy+/du+ = 1 + (k/E) [exp(k u+) - 1 - k u+ - (k u+)^2 / 2],
//! ```
//!
//! with `k = 0.41`, `E = 8.6`. That tangent viscosity drives interior
//! diffusion; the discrete wall flux over the half cell to a wall uses the
//! secant `nu_w / nu = y+ / u+` of the same profile, so the wall shear
//! `mu_w u_P / L` reproduces `rho u_tau^2` of the unified law of the wall
//! (using the tangent there over-predicts the friction: measured +61 % at
//! Re_Dh = 2e4 before this distinction). Both degrade to laminar flow as
//! `Re_L -> 0`. The energy equation takes the turbulent conductivity
//! `rho c_p nu_t / Pr_t` (`Pr_t = 0.9`).
//!
//! The wall distance runs from each fluid voxel centre to the box of the
//! nearest solid voxel or wall-type domain face (exact for axis-aligned
//! walls), the nearest seed from the separable linear-time feature
//! transform of Felzenszwalb & Huttenlocher (Theory of Computing 8, 2012) on
//! the cell centres of a uniform or graded grid padded by one ghost layer
//! behind each wall face.
//!
//! # No-claim boundaries
//!
//! LVEL is an algebraic, equilibrium model: no transport of turbulence,
//! history, separation physics, or transition prediction; its friction and
//! heat-transfer accuracy is that of the law of the wall (measured bands
//! are recorded on the validation fixtures, not promised elsewhere).

use super::domain::VoxelDomain;

/// von Karman constant of the LVEL profile.
pub const KAPPA: f64 = 0.41;
/// Log-law constant `E` of Spalding's profile.
pub const SPALDING_E: f64 = 8.6;
/// Turbulent Prandtl number of the energy equation.
pub const TURBULENT_PRANDTL: f64 = 0.9;

/// Spalding's `y+(u+)`.
fn spalding_y_plus(u: f64) -> f64 {
    let ku = KAPPA * u;
    u + (ku.exp() - 1.0 - ku - 0.5 * ku * ku - ku * ku * ku / 6.0) / SPALDING_E
}

/// `dy+/du+`, the effective-to-molecular viscosity ratio.
fn spalding_slope(u: f64) -> f64 {
    let ku = KAPPA * u;
    1.0 + KAPPA * (ku.exp() - 1.0 - ku - 0.5 * ku * ku) / SPALDING_E
}

/// Spalding's `u+` at local Reynolds number `Re_L = u+ y+(u+)` (Newton,
/// monotone, from the laminar sublayer guess `u+ = sqrt(Re)`).
fn spalding_u_plus(reynolds: f64) -> f64 {
    let mut u = reynolds.sqrt().min(60.0);
    for _ in 0..60 {
        let g = u * spalding_y_plus(u) - reynolds;
        let dg = spalding_y_plus(u) + u * spalding_slope(u);
        let next = (u - g / dg).max(0.5 * u);
        let done = (next - u).abs() <= 1e-12 * u.max(1.0);
        u = next;
        if done {
            break;
        }
    }
    u
}

/// The two viscosity ratios of the unified law of the wall at local
/// Reynolds number `Re_L = V L / nu`:
///
/// - `tangent = nu_eff / nu - 1 = dy+/du+ - 1`, the local eddy viscosity
///   ratio (interior diffusion);
/// - `secant = y+ / u+ = Re_L / u+^2`, the ratio that makes the discrete
///   wall flux over the distance `L`, `nu_w u / L`, reproduce the wall
///   shear `u_tau^2` of the law of the wall (wall-adjacent cells).
///
/// Both tend to the laminar values (0 and 1) as `Re_L -> 0`.
#[must_use]
pub fn law_of_the_wall_ratios(reynolds: f64) -> (f64, f64) {
    if !(reynolds > 0.0) {
        return (0.0, 1.0);
    }
    let u = spalding_u_plus(reynolds);
    (
        (spalding_slope(u) - 1.0).max(0.0),
        (reynolds / (u * u)).max(1.0),
    )
}

/// `nu_t / nu` (the tangent ratio of [`law_of_the_wall_ratios`]).
#[must_use]
pub fn eddy_viscosity_ratio(reynolds: f64) -> f64 {
    law_of_the_wall_ratios(reynolds).0
}

/// One-dimensional squared Euclidean distance transform with features
/// (Felzenszwalb & Huttenlocher, on arbitrary sample positions `x`):
/// `out[q] = min_p ((x_q - x_p)^2 + f[p])` over finite `f[p]`, and
/// `feature_out[q] = feature[p]` of the minimizing `p` (all infinite gives
/// infinity). `v` and `z` are scratch of length `f.len()` and
/// `f.len() + 1`.
#[allow(clippy::too_many_arguments)] // the transform's inputs and scratch
fn edt_1d(
    x: &[f64],
    f: &[f64],
    feature: &[[usize; 3]],
    out: &mut [f64],
    feature_out: &mut [[usize; 3]],
    v: &mut [usize],
    z: &mut [f64],
) {
    let n = f.len();
    // Lower envelope of the parabolas rooted at finite samples.
    let mut k: Option<usize> = None;
    for q in 0..n {
        if !f[q].is_finite() {
            continue;
        }
        let fq = x[q].mul_add(x[q], f[q]);
        loop {
            match k {
                None => {
                    k = Some(0);
                    v[0] = q;
                    z[0] = f64::NEG_INFINITY;
                    z[1] = f64::INFINITY;
                    break;
                }
                Some(top) => {
                    let p = v[top];
                    let s = (fq - x[p].mul_add(x[p], f[p])) / (2.0 * (x[q] - x[p]));
                    if s <= z[top] {
                        k = top.checked_sub(1);
                        continue;
                    }
                    let next = top + 1;
                    v[next] = q;
                    z[next] = s;
                    z[next + 1] = f64::INFINITY;
                    k = Some(next);
                    break;
                }
            }
        }
    }
    if k.is_none() {
        out[..n].fill(f64::INFINITY);
        return;
    }
    let mut j = 0usize;
    for q in 0..n {
        while z[j + 1] < x[q] {
            j += 1;
        }
        let p = v[j];
        let d = x[q] - x[p];
        out[q] = d.mul_add(d, f[p]);
        feature_out[q] = feature[p];
    }
}

/// Distance (m) from each fluid cell centre to the nearest wall: a solid
/// voxel or a domain face whose `wall[face]` is true (openings, inlets and
/// symmetry planes are not walls). The nearest seed cell comes from an
/// exact separable feature transform over the cell centres (uniform or
/// graded), with one ghost layer behind each wall face; the distance is
/// then to that seed cell's box (exact for axis-aligned walls), at least
/// half the cell's smallest width. Solid cells get zero.
#[must_use]
#[allow(clippy::too_many_lines)] // padding, three passes, box distances
pub fn wall_distance(domain: &VoxelDomain, wall: [bool; 6]) -> Vec<f64> {
    let [nx, ny, nz] = domain.dims();
    let n = [nx, ny, nz];
    // Padded grid: one ghost layer on every side; ghost cells are seeds on
    // wall faces, far away otherwise. Ghost centres mirror the edge cells.
    let dims = [nx + 2, ny + 2, nz + 2];
    let index = |x: usize, y: usize, z: usize| (z * dims[1] + y) * dims[0] + x;
    let width = |a: usize, i: usize| domain.width(a, i.clamp(1, n[a]) - 1);
    let position: [Vec<f64>; 3] = [0, 1, 2].map(|a| {
        (0..dims[a])
            .map(|i| {
                if i == 0 {
                    domain.face_coord(a, 0) - 0.5 * domain.width(a, 0)
                } else if i == dims[a] - 1 {
                    domain.face_coord(a, n[a]) + 0.5 * domain.width(a, n[a] - 1)
                } else {
                    0.5 * (domain.face_coord(a, i - 1) + domain.face_coord(a, i))
                }
            })
            .collect()
    });
    let mut grid = vec![f64::INFINITY; dims[0] * dims[1] * dims[2]];
    let mut feature = vec![[0usize; 3]; grid.len()];
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let inside = [x, y, z]
                    .iter()
                    .zip(&dims)
                    .all(|(&c, &d)| c > 0 && c + 1 < d);
                let seed = if inside {
                    !domain.is_fluid(domain.index(x - 1, y - 1, z - 1))
                } else {
                    // A ghost cell seeds when it lies behind exactly a wall
                    // face of the domain (any axis where it is outside).
                    let mut on_wall = false;
                    let mut on_open = false;
                    for (axis, &c) in [x, y, z].iter().enumerate() {
                        if c == 0 {
                            if wall[2 * axis] {
                                on_wall = true;
                            } else {
                                on_open = true;
                            }
                        } else if c + 1 == dims[axis] {
                            if wall[2 * axis + 1] {
                                on_wall = true;
                            } else {
                                on_open = true;
                            }
                        }
                    }
                    on_wall && !on_open
                };
                if seed {
                    grid[index(x, y, z)] = 0.0;
                    feature[index(x, y, z)] = [x, y, z];
                }
            }
        }
    }
    // Separable transform along x, then y, then z (squared metres).
    let longest = dims[0].max(dims[1]).max(dims[2]);
    let (mut f, mut out) = (vec![0.0; longest], vec![0.0; longest]);
    let (mut feat, mut feat_out) = (vec![[0usize; 3]; longest], vec![[0usize; 3]; longest]);
    let (mut v, mut zz) = (vec![0usize; longest], vec![0.0; longest + 2]);
    for axis in 0..3 {
        let len = dims[axis];
        let (u, w) = ((axis + 1) % 3, (axis + 2) % 3);
        for b in 0..dims[w] {
            for a in 0..dims[u] {
                let at = |i: usize| {
                    let mut c = [0usize; 3];
                    c[axis] = i;
                    c[u] = a;
                    c[w] = b;
                    index(c[0], c[1], c[2])
                };
                for i in 0..len {
                    f[i] = grid[at(i)];
                    feat[i] = feature[at(i)];
                }
                if f[..len].iter().all(|x| x.is_infinite()) {
                    continue;
                }
                edt_1d(
                    &position[axis],
                    &f[..len],
                    &feat[..len],
                    &mut out[..len],
                    &mut feat_out[..len],
                    &mut v,
                    &mut zz,
                );
                for i in 0..len {
                    grid[at(i)] = out[i];
                    feature[at(i)] = feat_out[i];
                }
            }
        }
    }
    (0..domain.cell_count())
        .map(|c| {
            if !domain.is_fluid(c) {
                return 0.0;
            }
            let [x, y, z] = domain.coords(c);
            let cell = index(x + 1, y + 1, z + 1);
            if !grid[cell].is_finite() {
                return f64::INFINITY;
            }
            // Distance from this centre to the nearest seed's box.
            let seed = feature[cell];
            let own = domain.widths(c);
            let mut squared = 0.0;
            for a in 0..3 {
                let gap = ((position[a][[x, y, z][a] + 1] - position[a][seed[a]]).abs()
                    - 0.5 * width(a, seed[a]))
                .max(0.0);
                squared = gap.mul_add(gap, squared);
            }
            let floor = 0.5 * own.iter().fold(f64::INFINITY, |m, w| m.min(*w));
            squared.sqrt().max(floor)
        })
        .collect()
}
