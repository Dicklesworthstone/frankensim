//! LVEL algebraic turbulence model (Agonafer, Gan-Li & Spalding, ASME HTD
//! 324, 1996) for the SIMPLEC voxel flow: the electronics-cooling model of
//! choice where fan-driven channel flows are transitional or turbulent and
//! meshes are too coarse for two-equation models.
//!
//! # Model
//!
//! At each fluid cell, the wall distance `L` and the local speed `V` relative
//! to that wall give `Re_L = V L / nu`, which fixes the law-of-the-wall
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
//! `mu_w |u_P - u_wall| / L` reproduces `rho u_tau^2` of the unified wall law
//! (using the tangent there over-predicts the friction: measured +61 % at
//! Re_Dh = 2e4 before this distinction). Both degrade to laminar flow as
//! `Re_L -> 0`. The energy equation takes the turbulent conductivity
//! `rho c_p nu_t / Pr_t` (`Pr_t = 0.9`).
//!
//! The wall distance runs from each fluid voxel centre to the box of the
//! nearest solid voxel or wall-type domain face (exact for axis-aligned
//! walls). A separable linear-time feature transform minimizes distance to
//! each seed's actual cell interval, including on graded grids. Minimizing
//! distance to solid centres and then measuring only the selected box would
//! choose the wrong surface near unequal cells and some obstacle corners.
//!
//! # No-claim boundaries
//!
//! LVEL is an algebraic, equilibrium model: no transport of turbulence,
//! history, separation physics, or transition prediction; its friction and
//! heat-transfer accuracy is that of the law of the wall (measured bands
//! are recorded on the validation fixtures, not promised elsewhere).

use super::domain::VoxelDomain;

/// Prandtl-aware equilibrium thermal wall resistance and coupling.
pub mod thermal_wall;

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

/// Intersection of the costs `f_p + dist(x, [a,b])^2` and
/// `f_q + dist(x, [c,d])^2`, where `a < b <= c < d` and
/// `offset = f_q - f_p`. Their difference increases monotonically, so
/// there is one switch from p to q. Depending on the offset, it lies in
/// the gap, inside one interval, or beyond both. The latter two cases are
/// essential on graded grids; a centre-parabola intersection is not valid.
fn interval_crossing([a, b]: [f64; 2], [c, d]: [f64; 2], offset: f64) -> f64 {
    let gap = c - b;
    let root = offset.abs().sqrt();
    if root <= gap {
        if gap == 0.0 {
            b
        } else {
            b + 0.5 * (gap + offset / gap)
        }
    } else if offset > 0.0 {
        if root <= d - b {
            b + root
        } else {
            b + 0.5 * ((d - b) + offset / (d - b))
        }
    } else if root <= c - a {
        c - root
    } else {
        a + 0.5 * ((c - a) + offset / (c - a))
    }
}

/// Squared distance transform to ordered, nonoverlapping cell intervals:
/// `out[q] = min_p (dist(x[q], [edges[p], edges[p+1]])^2 + f[p])`.
/// It retains the minimizing source feature. Each interval is pushed and
/// popped at most once, and each sample advances through the lower envelope
/// once, giving linear work as for the point-parabola transform. Exact
/// switch ties prefer the earlier interval. `v` and `z` are reusable scratch.
#[allow(clippy::too_many_arguments)] // the transform's inputs and scratch
fn edt_1d(
    x: &[f64],
    edges: &[f64],
    f: &[f64],
    feature: &[[usize; 3]],
    out: &mut [f64],
    feature_out: &mut [[usize; 3]],
    v: &mut [usize],
    z: &mut [f64],
) {
    let n = f.len();
    // Lower envelope of the squared-distance functions of finite intervals.
    let mut k: Option<usize> = None;
    for q in 0..n {
        if !f[q].is_finite() {
            continue;
        }
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
                    let s = interval_crossing(
                        [edges[p], edges[p + 1]],
                        [edges[q], edges[q + 1]],
                        f[q] - f[p],
                    );
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
        let d = (edges[p] - x[q]).max(x[q] - edges[p + 1]).max(0.0);
        out[q] = d.mul_add(d, f[p]);
        feature_out[q] = feature[p];
    }
}

/// Distance (m) from each fluid cell centre to the nearest wall: a solid
/// voxel or a domain face whose `wall[face]` is true (openings, inlets and
/// symmetry planes are not walls). The separable transform minimizes distance
/// to the actual solid boxes, including corners and unequal cell widths.
/// Ghost boxes behind wall faces represent those domain planes. The distance
/// is at least half the cell's smallest width. Solid cells get zero.
#[must_use]
pub fn wall_distance(domain: &VoxelDomain, wall: [bool; 6]) -> Vec<f64> {
    let mut distance = Vec::with_capacity(domain.cell_count());
    wall_features(domain, wall, |d, _| distance.push(d));
    distance
}

/// The same distance transform also retains the identity of a domain wall.
/// Internal solid voxels have no boundary face and remain stationary. Keeping
/// the feature selected by the transform ties the velocity reference to the
/// distance actually used, including deterministic ties between walls.
pub(super) fn wall_distance_and_boundary(
    domain: &VoxelDomain,
    wall: [bool; 6],
) -> (Vec<f64>, Vec<Option<u8>>) {
    let mut distance = Vec::with_capacity(domain.cell_count());
    let mut boundary = Vec::with_capacity(domain.cell_count());
    wall_features(domain, wall, |d, face| {
        distance.push(d);
        boundary.push(face);
    });
    (distance, boundary)
}

/// Share the distance/feature transform without a second grid scan or an
/// intermediate per-cell tuple allocation. The additional small boundary
/// cache is absent from laminar/SGS solvers.
#[allow(clippy::too_many_lines)] // padding and three interval-transform passes
fn wall_features(domain: &VoxelDomain, wall: [bool; 6], mut emit: impl FnMut(f64, Option<u8>)) {
    let [nx, ny, nz] = domain.dims();
    let n = [nx, ny, nz];
    // Padded grid: one ghost layer on every side; ghost cells are seeds on
    // wall faces, far away otherwise. Ghost intervals mirror the edge cells.
    let dims = [nx + 2, ny + 2, nz + 2];
    let index = |x: usize, y: usize, z: usize| (z * dims[1] + y) * dims[0] + x;
    let edges: [Vec<f64>; 3] = [0, 1, 2].map(|a| {
        (0..=dims[a])
            .map(|i| {
                if i == 0 {
                    domain.face_coord(a, 0) - domain.width(a, 0)
                } else if i == dims[a] {
                    domain.face_coord(a, n[a]) + domain.width(a, n[a] - 1)
                } else {
                    domain.face_coord(a, i - 1)
                }
            })
            .collect()
    });
    let position: [Vec<f64>; 3] =
        core::array::from_fn(|a| edges[a].windows(2).map(|e| 0.5 * (e[0] + e[1])).collect());
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
                    &edges[axis],
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
    for c in 0..domain.cell_count() {
        if !domain.is_fluid(c) {
            emit(0.0, None);
            continue;
        }
        let [x, y, z] = domain.coords(c);
        let cell = index(x + 1, y + 1, z + 1);
        if !grid[cell].is_finite() {
            emit(f64::INFINITY, None);
            continue;
        }
        // The interval transform already minimized the complete box distance;
        // this same feature selects the moving-wall velocity reference.
        let seed = feature[cell];
        let own = domain.widths(c);
        let floor = 0.5 * own.iter().fold(f64::INFINITY, |m, w| m.min(*w));
        let boundary = (0..3).find_map(|a| {
            if seed[a] == 0 {
                Some((2 * a) as u8)
            } else if seed[a] == n[a] + 1 {
                Some((2 * a + 1) as u8)
            } else {
                None
            }
        });
        emit(grid[cell].sqrt().max(floor), boundary);
    }
}

#[cfg(test)]
mod distance_tests {
    use super::*;
    use crate::conjugate::Voxel;

    #[test]
    fn interval_envelope_crosses_in_gaps_interiors_and_exterior_tails() {
        // G0: check the five distinct intersection cases against the defining
        // distances, independently of the three-dimensional transform.
        let left = [0.0, 2.0];
        let right = [3.0, 8.0];
        for offset in [-100.0_f64, -4.0, 0.0, 4.0, 100.0] {
            let x = interval_crossing(left, right, offset);
            let squared = |x: f64, interval: [f64; 2]| {
                (interval[0] - x).max(x - interval[1]).max(0.0).powi(2)
            };
            let residual = squared(x, left) - squared(x, right) - offset;
            assert!(
                residual.abs() < 1e-12,
                "offset={offset}, x={x}, residual={residual}"
            );
            assert!(squared(x - 0.01, left) - squared(x - 0.01, right) < offset);
            assert!(squared(x + 0.01, left) - squared(x + 0.01, right) > offset);
        }
        assert_eq!(interval_crossing([0.0, 2.0], [2.0, 3.0], 0.0), 2.0);
    }

    #[test]
    fn nearest_solid_centre_is_not_the_nearest_surface() {
        // The narrower solid has the nearest centre (2.75 vs 5.5), but
        // its surface is farther away (2.5 vs 1.5).
        let mut graded = VoxelDomain::graded([
            vec![1.0, 1.0, 8.0, 0.5],
            vec![1.0, 2.0, 0.5, 5.0],
            vec![1.0],
        ])
        .unwrap();
        graded.set(0, 2, 0, Voxel::Solid(0));
        graded.set(2, 0, 0, Voxel::Solid(0));
        let distances = wall_distance(&graded, [false; 6]);
        assert_eq!(distances[graded.index(0, 0, 0)], 1.5);

        // This also fails on a uniform grid: the farther centre's corner
        // lies closer than the other voxel's face (sqrt(40.5) vs 6.5).
        let mut uniform = VoxelDomain::new(8, 8, 1, 1.0).unwrap();
        uniform.set(0, 7, 0, Voxel::Solid(0));
        uniform.set(5, 5, 0, Voxel::Solid(0));
        assert_eq!(wall_distance(&uniform, [false; 6])[0], 40.5_f64.sqrt());
    }

    #[test]
    fn wall_distance_and_feature_match_direct_box_and_plane_geometry() {
        // G0/G3: the small direct oracle enumerates actual occupied boxes
        // and domain planes. It does not use centres, seeds or the transform.
        for widths in [
            [vec![1.0; 4], vec![1.0; 3], vec![1.0; 2]],
            [
                vec![0.25, 2.0, 0.5, 4.0],
                vec![3.0, 0.125, 1.0],
                vec![0.25, 4.0],
            ],
        ] {
            for scale in [1e-4, 1.0, 1e4] {
                let mut domain = VoxelDomain::graded(
                    widths
                        .clone()
                        .map(|axis| axis.into_iter().map(|w| w * scale).collect()),
                )
                .unwrap();
                for [x, y, z] in [[0, 2, 0], [2, 0, 1], [3, 2, 1]] {
                    domain.set(x, y, z, Voxel::Solid(0));
                }
                for mask in 0..64 {
                    let walls = core::array::from_fn(|side| mask & (1 << side) != 0);
                    let (distance, boundary) = wall_distance_and_boundary(&domain, walls);
                    for cell in 0..domain.cell_count() {
                        if !domain.is_fluid(cell) {
                            assert_eq!(distance[cell], 0.0);
                            assert_eq!(boundary[cell], None);
                            continue;
                        }
                        let [x, y, z] = domain.coords(cell);
                        let point = domain.center(x, y, z);
                        let mut solid_distance = f64::INFINITY;
                        for other in 0..domain.cell_count() {
                            if domain.is_fluid(other) {
                                continue;
                            }
                            let at = domain.coords(other);
                            let squared: f64 = (0..3)
                                .map(|axis| {
                                    (domain.face_coord(axis, at[axis]) - point[axis])
                                        .max(point[axis] - domain.face_coord(axis, at[axis] + 1))
                                        .max(0.0)
                                        .powi(2)
                                })
                                .sum();
                            solid_distance = solid_distance.min(squared.sqrt());
                        }
                        let plane_distance = |side: usize| {
                            let axis = side / 2;
                            if side % 2 == 0 {
                                point[axis]
                            } else {
                                domain.face_coord(axis, domain.dims()[axis]) - point[axis]
                            }
                        };
                        let expected = (0..6)
                            .filter(|&side| walls[side])
                            .map(plane_distance)
                            .fold(solid_distance, f64::min);
                        let selected = boundary[cell].map_or(solid_distance, |side| {
                            assert!(walls[usize::from(side)]);
                            plane_distance(usize::from(side))
                        });
                        assert!(
                            (distance[cell] - expected).abs() < 2e-14 * scale,
                            "mask={mask}, cell={cell}, got={}, expected={expected}",
                            distance[cell]
                        );
                        assert!(
                            (selected - expected).abs() < 2e-14 * scale,
                            "wrong wall frame: mask={mask}, cell={cell}, boundary={:?}",
                            boundary[cell]
                        );
                    }
                }
            }
        }
        let empty = VoxelDomain::new(2, 2, 2, 1.0).unwrap();
        let (distance, boundary) = wall_distance_and_boundary(&empty, [false; 6]);
        assert!(distance.iter().all(|d| d.is_infinite()));
        assert!(boundary.iter().all(Option::is_none));
    }
}
