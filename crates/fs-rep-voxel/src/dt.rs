//! EXACT Euclidean distance transform (Felzenszwalb–Huttenlocher lower
//! envelopes, separable over the three axes) over the active set's
//! bounding box. Squared distances are computed in integer-exact voxel
//! units — the conformance suite checks EQUALITY against the O(n²)
//! brute force, not a tolerance. Deterministic: fixed pass order.

use crate::{VoxelError, field::OccupancyField};

#[derive(Debug, Clone, Copy)]
pub(crate) struct CheckedBox {
    pub(crate) min: [i32; 3],
    pub(crate) max: [i32; 3],
    pub(crate) dims: [usize; 3],
    pub(crate) total: usize,
}

pub(crate) fn active_bounds(field: &OccupancyField) -> Option<([i32; 3], [i32; 3])> {
    let mut min = [i32::MAX; 3];
    let mut max = [i32::MIN; 3];
    let mut found = false;
    for (c, _) in field.grid.iter_active() {
        found = true;
        for axis in 0..3 {
            min[axis] = min[axis].min(c[axis]);
            max[axis] = max[axis].max(c[axis]);
        }
    }
    found.then_some((min, max))
}

pub(crate) fn checked_dense_box(
    min: [i32; 3],
    max: [i32; 3],
    halo: u32,
    max_voxels: usize,
    operation: &'static str,
) -> Result<CheckedBox, VoxelError> {
    let mut padded_min = [0i32; 3];
    let mut padded_max = [0i32; 3];
    let mut dims_u64 = [0u64; 3];
    let halo_i64 = i64::from(halo);
    for axis in 0..3 {
        let lo = i64::from(min[axis]) - halo_i64;
        let hi = i64::from(max[axis]) + halo_i64;
        padded_min[axis] = i32::try_from(lo).map_err(|_| VoxelError::CoordinateRange {
            operation,
            axis,
            min: min[axis],
            max: max[axis],
            halo,
        })?;
        padded_max[axis] = i32::try_from(hi).map_err(|_| VoxelError::CoordinateRange {
            operation,
            axis,
            min: min[axis],
            max: max[axis],
            halo,
        })?;
        let span = hi - lo + 1;
        dims_u64[axis] = u64::try_from(span).map_err(|_| VoxelError::CoordinateRange {
            operation,
            axis,
            min: min[axis],
            max: max[axis],
            halo,
        })?;
    }

    let required = dims_u64
        .iter()
        .try_fold(1u128, |volume, &dim| volume.checked_mul(u128::from(dim)))
        .ok_or(VoxelError::DenseVolumeOverflow {
            operation,
            dims: dims_u64,
        })?;
    if required > max_voxels as u128 {
        return Err(VoxelError::VoxelBudgetExceeded {
            operation,
            required,
            maximum: max_voxels,
        });
    }
    let to_usize = |dim| {
        usize::try_from(dim).map_err(|_| VoxelError::DenseVolumeOverflow {
            operation,
            dims: dims_u64,
        })
    };
    let dims = [
        to_usize(dims_u64[0])?,
        to_usize(dims_u64[1])?,
        to_usize(dims_u64[2])?,
    ];
    let total = usize::try_from(required).map_err(|_| VoxelError::DenseVolumeOverflow {
        operation,
        dims: dims_u64,
    })?;
    Ok(CheckedBox {
        min: padded_min,
        max: padded_max,
        dims,
        total,
    })
}

fn ensure_vec_capacity<T>(
    len: usize,
    operation: &'static str,
    dims: [usize; 3],
) -> Result<(), VoxelError> {
    let max_len = (isize::MAX as usize) / core::mem::size_of::<T>().max(1);
    if len <= max_len {
        return Ok(());
    }
    Err(VoxelError::DenseVolumeOverflow {
        operation,
        dims: dims.map(|dim| u64::try_from(dim).unwrap_or(u64::MAX)),
    })
}

/// A dense distance field over the active set's bounding box.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceField {
    /// Bounding-box min voxel coordinate.
    min: [i32; 3],
    /// Grid dimensions (voxels).
    dims: [usize; 3],
    /// Squared distance (voxel units) to the nearest ACTIVE voxel,
    /// row-major x-fastest; `f64::INFINITY` when there are no seeds.
    sq: Vec<f64>,
    /// Voxel edge length (m) — converts voxel distances to world.
    voxel_size: f64,
}

impl DistanceField {
    /// Bounding-box minimum voxel coordinate.
    #[must_use]
    pub fn min(&self) -> [i32; 3] {
        self.min
    }

    /// Dense bounding-box dimensions.
    #[must_use]
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// Checked row-major squared distances in voxel units.
    #[must_use]
    pub fn squared_distances(&self) -> &[f64] {
        &self.sq
    }

    /// Voxel edge length used to convert distances to world units.
    #[must_use]
    pub fn voxel_size(&self) -> f64 {
        self.voxel_size
    }

    /// Euclidean distance (world units) at a voxel inside the box.
    #[must_use]
    pub fn distance(&self, coord: [i32; 3]) -> Option<f64> {
        let idx = self.index(coord)?;
        Some(fs_math::det::sqrt(*self.sq.get(idx)?) * self.voxel_size)
    }

    fn index(&self, coord: [i32; 3]) -> Option<usize> {
        let mut idx = 0usize;
        let mut stride = 1usize;
        for ((&c, &lo), &dim) in coord.iter().zip(&self.min).zip(&self.dims) {
            let rel = i64::from(c) - i64::from(lo);
            if rel.is_negative() {
                return None;
            }
            let rel = usize::try_from(rel).ok()?;
            if rel >= dim {
                return None;
            }
            idx = idx.checked_add(rel.checked_mul(stride)?)?;
            stride = stride.checked_mul(dim)?;
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
    #[allow(clippy::cast_precision_loss)] // authorized dense line length
    let sq = |x: usize| {
        let x = x as f64;
        x * x
    };
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
///
/// `max_voxels` is an explicit authorization for the dense bounding-box
/// allocation; sparse active count is not used as an allocation proxy.
///
/// # Errors
/// Returns a structured coordinate, volume, or budget error before
/// allocation when the active bounding box is not admissible.
pub fn euclidean_dt(
    field: &OccupancyField,
    max_voxels: usize,
) -> Result<Option<DistanceField>, VoxelError> {
    let Some((min, max)) = active_bounds(field) else {
        return Ok(None);
    };
    let checked = checked_dense_box(min, max, 0, max_voxels, "euclidean distance transform")?;
    let dims = checked.dims;
    let total = checked.total;
    ensure_vec_capacity::<f64>(total, "euclidean distance transform", dims)?;
    let mut sq = vec![f64::INFINITY; total];
    let at = |x: usize, y: usize, z: usize| x + dims[0] * (y + dims[1] * z);
    for (c, _) in field.grid.iter_active() {
        let x = usize::try_from(i64::from(c[0]) - i64::from(min[0])).map_err(|_| {
            VoxelError::CoordinateRange {
                operation: "euclidean distance transform",
                axis: 0,
                min: min[0],
                max: max[0],
                halo: 0,
            }
        })?;
        let y = usize::try_from(i64::from(c[1]) - i64::from(min[1])).map_err(|_| {
            VoxelError::CoordinateRange {
                operation: "euclidean distance transform",
                axis: 1,
                min: min[1],
                max: max[1],
                halo: 0,
            }
        })?;
        let z = usize::try_from(i64::from(c[2]) - i64::from(min[2])).map_err(|_| {
            VoxelError::CoordinateRange {
                operation: "euclidean distance transform",
                axis: 2,
                min: min[2],
                max: max[2],
                halo: 0,
            }
        })?;
        sq[at(x, y, z)] = 0.0;
    }
    let max_dim = dims.iter().copied().max().unwrap_or(1);
    ensure_vec_capacity::<f64>(max_dim, "euclidean distance transform scratch", dims)?;
    ensure_vec_capacity::<usize>(max_dim, "euclidean distance transform scratch", dims)?;
    let z_len = max_dim
        .checked_add(1)
        .ok_or(VoxelError::DenseVolumeOverflow {
            operation: "euclidean distance transform scratch",
            dims: dims.map(|dim| u64::try_from(dim).unwrap_or(u64::MAX)),
        })?;
    ensure_vec_capacity::<f64>(z_len, "euclidean distance transform scratch", dims)?;
    let mut line = vec![0.0f64; max_dim];
    let mut v = vec![0usize; max_dim];
    let mut zbuf = vec![0.0f64; z_len];
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
    Ok(Some(DistanceField {
        min: checked.min,
        dims,
        sq,
        voxel_size: field.voxel_size(),
    }))
}
