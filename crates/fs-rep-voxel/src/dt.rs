//! EXACT Euclidean distance transform (Felzenszwalb–Huttenlocher lower
//! envelopes, separable over the three axes) over the active set's
//! bounding box. Squared distances are computed in integer-exact voxel
//! units — the conformance suite checks EQUALITY against the O(n²)
//! brute force, not a tolerance. Deterministic: fixed pass order.

use crate::field::OccupancyField;

/// A dense distance field over the active set's bounding box.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceField {
    /// Bounding-box min voxel coordinate.
    pub min: [i32; 3],
    /// Grid dimensions (voxels).
    pub dims: [usize; 3],
    /// Squared distance (voxel units) to the nearest ACTIVE voxel,
    /// row-major x-fastest; `f64::INFINITY` when there are no seeds.
    pub sq: Vec<f64>,
    /// Voxel edge length (m) — converts voxel distances to world.
    pub voxel_size: f64,
}

impl DistanceField {
    /// Euclidean distance (world units) at a voxel inside the box.
    #[must_use]
    pub fn distance(&self, coord: [i32; 3]) -> Option<f64> {
        let idx = self.index(coord)?;
        Some(fs_math::det::sqrt(self.sq[idx]) * self.voxel_size)
    }

    fn index(&self, coord: [i32; 3]) -> Option<usize> {
        let mut idx = 0usize;
        let mut stride = 1usize;
        for ((&c, &lo), &dim) in coord.iter().zip(&self.min).zip(&self.dims) {
            let rel = c.checked_sub(lo)?;
            if rel < 0 {
                return None;
            }
            let rel = usize::try_from(rel).ok()?;
            if rel >= dim {
                return None;
            }
            idx += rel * stride;
            stride *= dim;
        }
        Some(idx)
    }
}

/// One-dimensional exact squared-distance transform (lower envelope of
/// parabolas). `f` holds squared distances; output overwrites it.
fn dt_1d(f: &mut [f64], v: &mut [usize], z: &mut [f64]) {
    let n = f.len();
    // The envelope is built from FINITE parabolas only (+inf sources can
    // never be nearest); an all-infinite line stays infinite.
    let Some(first) = (0..n).find(|&i| f[i].is_finite()) else {
        return;
    };
    let mut k = 0usize;
    v[0] = first;
    z[0] = f64::NEG_INFINITY;
    z[1] = f64::INFINITY;
    #[allow(clippy::cast_precision_loss)] // line lengths are small
    let sq = |x: usize| (x * x) as f64;
    for q in (first + 1)..n {
        if !f[q].is_finite() {
            continue;
        }
        loop {
            let p = v[k];
            let s = ((f[q] + sq(q)) - (f[p] + sq(p))) / (2.0 * (q as f64 - p as f64));
            if s <= z[k] {
                if k == 0 {
                    break;
                }
                k -= 1;
            } else {
                k += 1;
                v[k] = q;
                z[k] = s;
                z[k + 1] = f64::INFINITY;
                break;
            }
        }
    }
    k = 0;
    let out: Vec<f64> = (0..n)
        .map(|q| {
            while z[k + 1] < q as f64 {
                k += 1;
            }
            let p = v[k];
            let d = q as f64 - p as f64;
            d * d + f[p]
        })
        .collect();
    f.copy_from_slice(&out);
}

/// Exact Euclidean DT of an occupancy field over its active bounding box
/// (distance TO the active set; active voxels get 0). Returns `None` for
/// an empty active set.
#[must_use]
pub fn euclidean_dt(field: &OccupancyField) -> Option<DistanceField> {
    let mut min = [i32::MAX; 3];
    let mut max = [i32::MIN; 3];
    for (c, _) in field.grid.iter_active() {
        for k in 0..3 {
            min[k] = min[k].min(c[k]);
            max[k] = max[k].max(c[k]);
        }
    }
    if min[0] == i32::MAX {
        return None;
    }
    let dims: [usize; 3] =
        core::array::from_fn(|k| usize::try_from(max[k] - min[k] + 1).expect("bounded box"));
    let total = dims[0] * dims[1] * dims[2];
    let mut sq = vec![f64::INFINITY; total];
    let at = |x: usize, y: usize, z: usize| x + dims[0] * (y + dims[1] * z);
    for (c, _) in field.grid.iter_active() {
        let x = usize::try_from(c[0] - min[0]).expect("in box");
        let y = usize::try_from(c[1] - min[1]).expect("in box");
        let z = usize::try_from(c[2] - min[2]).expect("in box");
        sq[at(x, y, z)] = 0.0;
    }
    let max_dim = dims.iter().copied().max().unwrap_or(1);
    let mut line = vec![0.0f64; max_dim];
    let mut v = vec![0usize; max_dim];
    let mut zbuf = vec![0.0f64; max_dim + 1];
    // Pass 1: x lines.
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for (i, slot) in line.iter_mut().take(dims[0]).enumerate() {
                *slot = sq[at(i, y, z)];
            }
            dt_1d(&mut line[..dims[0]], &mut v, &mut zbuf);
            for i in 0..dims[0] {
                sq[at(i, y, z)] = line[i];
            }
        }
    }
    // Pass 2: y lines.
    for z in 0..dims[2] {
        for x in 0..dims[0] {
            for (i, slot) in line.iter_mut().take(dims[1]).enumerate() {
                *slot = sq[at(x, i, z)];
            }
            dt_1d(&mut line[..dims[1]], &mut v, &mut zbuf);
            for i in 0..dims[1] {
                sq[at(x, i, z)] = line[i];
            }
        }
    }
    // Pass 3: z lines.
    for y in 0..dims[1] {
        for x in 0..dims[0] {
            for (i, slot) in line.iter_mut().take(dims[2]).enumerate() {
                *slot = sq[at(x, y, i)];
            }
            dt_1d(&mut line[..dims[2]], &mut v, &mut zbuf);
            for i in 0..dims[2] {
                sq[at(x, y, i)] = line[i];
            }
        }
    }
    Some(DistanceField {
        min,
        dims,
        sq,
        voxel_size: field.voxel_size,
    })
}
