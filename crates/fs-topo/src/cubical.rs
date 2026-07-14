//! Cubical topology: exact Betti numbers of voxel solids (union-find
//! components + exact Euler characteristic + complement cavities, with
//! `b1` closed by duality), true 0-dimensional persistence over the
//! voxel filtration (elder rule), persistence-thresholded feature
//! counting, and chart-level verification with HONEST resolution
//! caveats. Sequential and deterministic (P2); the chunked-parallel
//! reduction for 10⁸⁺-voxel fields is a CONTRACT no-claim routed to
//! the perf lane.

use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ClippedChart, Point3, SamplingDomain, SamplingDomainError};

/// Deterministic upper bound on cells admitted by chart voxelization.
pub const MAX_VOXELIZE_CELLS: u64 = 1_000_000;

/// Structured chart-voxelization failure.
#[derive(Debug, Clone, PartialEq)]
pub enum VoxelizeError {
    /// The source support or explicit clip is not an admissible finite
    /// three-dimensional sampling domain.
    SamplingDomain(SamplingDomainError),
    /// The longest-axis cell resolution was zero.
    InvalidResolution {
        /// The offending resolution.
        n: u32,
    },
    /// The requested resolution produced a non-finite or zero cell size.
    InvalidCellSize {
        /// The requested longest-axis resolution.
        n: u32,
        /// Longest admitted domain span.
        longest: f64,
    },
    /// Per-axis dimension derivation or checked multiplication overflowed.
    VoxelCountOverflow {
        /// Per-axis dimensions involved in the overflow.
        dims: [u64; 3],
    },
    /// The requested field exceeds the deterministic voxel-work cap.
    VoxelLimit {
        /// Per-axis cell dimensions.
        dims: [u64; 3],
        /// Total cells required.
        need: u128,
        /// Deterministic cell cap.
        cap: u64,
    },
    /// A chart returned a non-finite nominal field value.
    InvalidSample {
        /// Point at which evaluation failed.
        point: Point3,
        /// Raw bits of the rejected value.
        value_bits: u64,
    },
    /// Voxelization observed cancellation at a bounded polling point.
    Cancelled {
        /// Voxels fully evaluated before cancellation was observed.
        completed_voxels: u64,
    },
}

impl core::fmt::Display for VoxelizeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SamplingDomain(error) => write!(f, "{error}"),
            Self::InvalidResolution { n } => write!(
                f,
                "voxelization refused: longest-axis resolution must be positive, got {n}"
            ),
            Self::InvalidCellSize { n, longest } => write!(
                f,
                "voxelization refused: span {longest} at resolution {n} has no finite positive cell size"
            ),
            Self::VoxelCountOverflow { dims } => write!(
                f,
                "voxelization refused: cell dimensions {dims:?} overflow the addressable field"
            ),
            Self::VoxelLimit { dims, need, cap } => write!(
                f,
                "voxelization refused: field {dims:?} requires {need} cells, exceeding the {cap} cell cap; lower n or shrink the clip"
            ),
            Self::InvalidSample { point, value_bits } => write!(
                f,
                "voxelization refused: chart sample at {point:?} is non-finite (f64 bits {value_bits:#018x})"
            ),
            Self::Cancelled { completed_voxels } => write!(
                f,
                "voxelization cancelled after {completed_voxels} completed voxels"
            ),
        }
    }
}

impl core::error::Error for VoxelizeError {}

impl From<SamplingDomainError> for VoxelizeError {
    fn from(error: SamplingDomainError) -> Self {
        Self::SamplingDomain(error)
    }
}

/// A dense voxel field (values at cell centers, `x` fastest).
#[derive(Debug, Clone)]
pub struct VoxelField {
    /// Cells per axis.
    pub dims: [u32; 3],
    /// Cell-center values (signed distance or density).
    pub values: Vec<f64>,
    /// Requested maximum cell width (metadata). Each axis is partitioned over
    /// its exact span, so shorter-axis widths may be smaller.
    pub h: f64,
}

impl VoxelField {
    fn idx(&self, x: u32, y: u32, z: u32) -> usize {
        (usize::try_from(z).expect("u32 fits usize")
            * usize::try_from(self.dims[1]).expect("u32 fits usize")
            + usize::try_from(y).expect("u32 fits usize"))
            * usize::try_from(self.dims[0]).expect("u32 fits usize")
            + usize::try_from(x).expect("u32 fits usize")
    }

    /// Occupancy at threshold: `value < level`.
    fn filled(&self, x: u32, y: u32, z: u32, level: f64) -> bool {
        self.values[self.idx(x, y, z)] < level
    }
}

/// Voxelize a chart's signed distance over its support at `n` cells on
/// the longest axis. Every axis is partitioned over its exact admitted span;
/// cell widths are at most the longest-axis `h` and centers remain inside the
/// support even for extreme aspect ratios.
///
/// # Errors
/// [`VoxelizeError`] when the resolution, finite sampling domain, checked
/// cell count, or deterministic work cap is inadmissible, or when
/// cancellation is observed at a bounded polling point.
pub fn voxelize(chart: &dyn Chart, n: u32, cx: &Cx<'_>) -> Result<VoxelField, VoxelizeError> {
    if n == 0 {
        return Err(VoxelizeError::InvalidResolution { n });
    }
    let domain = SamplingDomain::admit(chart.support(), None)?;
    let b = domain.bounds();
    let spans = domain.spans();
    let ext = [spans.x, spans.y, spans.z];
    let longest = domain.max_span();
    let h = longest / f64::from(n);
    if !h.is_finite() || h <= 0.0 {
        return Err(VoxelizeError::InvalidCellSize { n, longest });
    }
    let mut dims_u64 = [0u64; 3];
    for axis in 0..3 {
        let cells = (ext[axis] / h).ceil().max(1.0);
        if !cells.is_finite() || cells > f64::from(u32::MAX) {
            dims_u64[axis] = u64::MAX;
            return Err(VoxelizeError::VoxelCountOverflow { dims: dims_u64 });
        }
        let mut cells = cells as u64;
        // A nearest-rounded quotient can land exactly on an integer even when
        // the realized stored width is one ulp larger than the public `h`.
        // Validate the actual per-axis width and increment until the metadata
        // promise is true in the same f64 arithmetic consumers observe.
        loop {
            let width_cells = u32::try_from(cells)
                .map_err(|_| VoxelizeError::VoxelCountOverflow { dims: dims_u64 })?;
            let realized_width = ext[axis] / f64::from(width_cells);
            if !realized_width.is_finite() {
                dims_u64[axis] = u64::MAX;
                return Err(VoxelizeError::VoxelCountOverflow { dims: dims_u64 });
            }
            if realized_width <= h {
                break;
            }
            cells = cells
                .checked_add(1)
                .ok_or(VoxelizeError::VoxelCountOverflow { dims: dims_u64 })?;
            if cells > u64::from(u32::MAX) {
                dims_u64[axis] = cells;
                return Err(VoxelizeError::VoxelCountOverflow { dims: dims_u64 });
            }
        }
        dims_u64[axis] = cells;
    }
    let need = dims_u64
        .iter()
        .try_fold(1u128, |product, &dim| product.checked_mul(u128::from(dim)))
        .ok_or(VoxelizeError::VoxelCountOverflow { dims: dims_u64 })?;
    if need > u128::from(MAX_VOXELIZE_CELLS) {
        return Err(VoxelizeError::VoxelLimit {
            dims: dims_u64,
            need,
            cap: MAX_VOXELIZE_CELLS,
        });
    }
    let total =
        usize::try_from(need).map_err(|_| VoxelizeError::VoxelCountOverflow { dims: dims_u64 })?;
    let dims = dims_u64.map(|dim| u32::try_from(dim).expect("u32 dimension checked"));
    let mut values = Vec::with_capacity(total);
    let coordinate = |axis: usize, index: u32| {
        let lo = match axis {
            0 => b.min.x,
            1 => b.min.y,
            _ => b.min.z,
        };
        let fraction = (f64::from(index) + 0.5) / f64::from(dims[axis]);
        lo + ext[axis] * fraction
    };
    let mut completed_voxels = 0u64;
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                if completed_voxels.is_multiple_of(256) {
                    cx.checkpoint()
                        .map_err(|_| VoxelizeError::Cancelled { completed_voxels })?;
                }
                let p = Point3::new(coordinate(0, x), coordinate(1, y), coordinate(2, z));
                let value = chart.eval(p, cx).signed_distance;
                if !value.is_finite() {
                    return Err(VoxelizeError::InvalidSample {
                        point: p,
                        value_bits: value.to_bits(),
                    });
                }
                values.push(value);
                completed_voxels += 1;
            }
        }
    }
    cx.checkpoint()
        .map_err(|_| VoxelizeError::Cancelled { completed_voxels })?;
    Ok(VoxelField { dims, values, h })
}

/// Voxelize the geometric intersection `chart ∩ clip` at `n` cells on the
/// longest axis. The explicit clip is part of the sampled chart field.
///
/// # Errors
/// [`VoxelizeError`] under the same conditions as [`voxelize`], plus an
/// invalid, empty, or degenerate explicit clip.
pub fn voxelize_clipped(
    chart: &dyn Chart,
    n: u32,
    clip: Aabb,
    cx: &Cx<'_>,
) -> Result<VoxelField, VoxelizeError> {
    if n == 0 {
        return Err(VoxelizeError::InvalidResolution { n });
    }
    let clipped = ClippedChart::new(chart, clip)?;
    voxelize(&clipped, n, cx)
}

struct UnionFind {
    parent: Vec<u32>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n as u32).collect(),
        }
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            self.parent[x as usize] = self.parent[self.parent[x as usize] as usize];
            x = self.parent[x as usize];
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        // Deterministic: smaller root wins (index order, not rank).
        let (lo, hi) = (ra.min(rb), ra.max(rb));
        self.parent[hi as usize] = lo;
        true
    }
}

/// Exact Betti numbers of the voxel SOLID `{value < level}` (union of
/// closed unit cubes): `b0` by 26-CONNECTED union-find (closed cubes touch
/// through faces, edges, and corners), `b2` as bounded complement components
/// (6-connected empty cells against a virtual outside — the (26, 6)
/// digital-topology pair), and `b1 = b0 + b2 − χ` with χ counted EXACTLY on
/// the cubical complex, so the three are mutually consistent (bead rcvl).
#[must_use]
#[allow(clippy::too_many_lines)] // b0, b2, and the exact Euler count are one derivation
pub fn betti(field: &VoxelField, level: f64) -> (u32, u32, u32) {
    let [nx, ny, nz] = field.dims;
    let total = field.values.len();
    // b0: components of filled cells under 26-CONNECTIVITY. The solid is a
    // union of CLOSED unit cubes, so cubes sharing a face, EDGE, or corner are
    // topologically connected — the digital-topology partner of the 6-connected
    // VOID model used for b2 and the closed-cube χ counted below (the (26, 6)
    // pair). A 6-connected b0 was internally inconsistent with χ, inflating
    // b1 = b0 + b2 − χ into PHANTOM TUNNELS on any diagonal contact (bead rcvl).
    // The 13 already-visited neighbours in the 3×3×3 stencil (z→y→x order).
    #[rustfmt::skip]
    #[allow(clippy::items_after_statements)] // the stencil belongs with the b0 loop
    const BACKWARD_26: [(isize, isize, isize); 13] = [
        (-1, -1, -1), (-1, -1, 0), (-1, -1, 1),
        (-1,  0, -1), (-1,  0, 0), (-1,  0, 1),
        (-1,  1, -1), (-1,  1, 0), (-1,  1, 1),
        ( 0, -1, -1), ( 0, -1, 0), ( 0, -1, 1),
        ( 0,  0, -1),
    ];
    let mut uf = UnionFind::new(total);
    let mut filled_count = 0u64;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if !field.filled(x, y, z, level) {
                    continue;
                }
                filled_count += 1;
                let i = field.idx(x, y, z) as u32;
                for &(dz, dy, dx) in &BACKWARD_26 {
                    let (nxp, nyp, nzp) = (
                        i64::from(x) + dx as i64,
                        i64::from(y) + dy as i64,
                        i64::from(z) + dz as i64,
                    );
                    if nxp < 0 || nyp < 0 || nzp < 0 {
                        continue;
                    }
                    let (nxp, nyp, nzp) = (nxp as u32, nyp as u32, nzp as u32);
                    if nxp >= nx || nyp >= ny || nzp >= nz {
                        continue;
                    }
                    if field.filled(nxp, nyp, nzp, level) {
                        uf.union(i, field.idx(nxp, nyp, nzp) as u32);
                    }
                }
            }
        }
    }
    if filled_count == 0 {
        return (0, 0, 0);
    }
    let mut b0 = 0u32;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = field.idx(x, y, z) as u32;
                if field.filled(x, y, z, level) && uf.find(i) == i {
                    b0 += 1;
                }
            }
        }
    }
    // b2: bounded components of EMPTY cells (union with the virtual
    // outside node; unbounded = touching the boundary).
    let outside = total as u32;
    let mut ufc = UnionFind::new(total + 1);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if field.filled(x, y, z, level) {
                    continue;
                }
                let i = field.idx(x, y, z) as u32;
                if x == 0 || y == 0 || z == 0 || x == nx - 1 || y == ny - 1 || z == nz - 1 {
                    ufc.union(i, outside);
                }
                if x > 0 && !field.filled(x - 1, y, z, level) {
                    ufc.union(i, field.idx(x - 1, y, z) as u32);
                }
                if y > 0 && !field.filled(x, y - 1, z, level) {
                    ufc.union(i, field.idx(x, y - 1, z) as u32);
                }
                if z > 0 && !field.filled(x, y, z - 1, level) {
                    ufc.union(i, field.idx(x, y, z - 1) as u32);
                }
            }
        }
    }
    let out_root = ufc.find(outside);
    let mut b2 = 0u32;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = field.idx(x, y, z) as u32;
                if !field.filled(x, y, z, level) && ufc.find(i) == i && i != out_root {
                    b2 += 1;
                }
            }
        }
    }
    // χ = V − E + F − C on the closed cubical complex: a k-cell is
    // present iff ANY incident voxel is filled.
    let filled_at = |x: i64, y: i64, z: i64| -> bool {
        x >= 0
            && y >= 0
            && z >= 0
            && (x as u32) < nx
            && (y as u32) < ny
            && (z as u32) < nz
            && field.filled(x as u32, y as u32, z as u32, level)
    };
    let mut v_count = 0i64;
    let mut e_count = 0i64;
    let mut f_count = 0i64;
    for z in 0..=i64::from(nz) {
        for y in 0..=i64::from(ny) {
            for x in 0..=i64::from(nx) {
                // Vertex at lattice point (x,y,z): 8 incident voxels.
                let mut any = false;
                for dz in -1..=0 {
                    for dy in -1..=0 {
                        for dx in -1..=0 {
                            any |= filled_at(x + dx, y + dy, z + dz);
                        }
                    }
                }
                if any {
                    v_count += 1;
                }
                // Edges along +x/+y/+z from this lattice point: 4
                // incident voxels each.
                // Edges along +x/+y/+z: 4 incident voxels each, with
                // offsets ONLY in the two perpendicular axes.
                for axis in 0..3usize {
                    if (axis == 0 && x >= i64::from(nx))
                        || (axis == 1 && y >= i64::from(ny))
                        || (axis == 2 && z >= i64::from(nz))
                    {
                        continue;
                    }
                    let mut any_e = false;
                    for u in -1..=0i64 {
                        for w in -1..=0i64 {
                            let (ex, ey, ez) = match axis {
                                0 => (x, y + u, z + w),
                                1 => (x + u, y, z + w),
                                _ => (x + u, y + w, z),
                            };
                            any_e |= filled_at(ex, ey, ez);
                        }
                    }
                    if any_e {
                        e_count += 1;
                    }
                }
                // Faces normal to each axis at this lattice corner: 2
                // incident voxels each.
                let face_ok: [bool; 3] = [
                    y < i64::from(ny) && z < i64::from(nz),
                    x < i64::from(nx) && z < i64::from(nz),
                    x < i64::from(nx) && y < i64::from(ny),
                ];
                for (axis, &ok) in face_ok.iter().enumerate() {
                    if !ok {
                        continue;
                    }
                    let (c1, c2) = match axis {
                        0 => ((x - 1, y, z), (x, y, z)),
                        1 => ((x, y - 1, z), (x, y, z)),
                        _ => ((x, y, z - 1), (x, y, z)),
                    };
                    if filled_at(c1.0, c1.1, c1.2) || filled_at(c2.0, c2.1, c2.2) {
                        f_count += 1;
                    }
                }
            }
        }
    }
    let filled_i64 = i64::try_from(filled_count).expect("voxel counts fit i64");
    let chi = v_count - e_count + f_count - filled_i64;
    let b1 = (i64::from(b0) + i64::from(b2) - chi).max(0) as u32;
    (b0, b1, b2)
}

/// One 0-dimensional persistence bar.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bar {
    /// Birth value (component appears).
    pub birth: f64,
    /// Death value (`f64::INFINITY` for essential components).
    pub death: f64,
}

impl Bar {
    /// Persistence (lifetime).
    #[must_use]
    pub fn persistence(&self) -> f64 {
        self.death - self.birth
    }
}

/// True 0-dimensional persistence of the sublevel filtration (elder
/// rule: the younger component dies at each merge). Deterministic:
/// voxels sorted by (value, index).
#[must_use]
pub fn persistence0(field: &VoxelField) -> Vec<Bar> {
    let n = field.values.len();
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_by(|&a, &b| {
        field.values[a as usize]
            .partial_cmp(&field.values[b as usize])
            .expect("finite values")
            .then(a.cmp(&b))
    });
    let [nx, ny, nz] = field.dims;
    let mut uf = UnionFind::new(n);
    let mut birth: Vec<f64> = vec![f64::NAN; n]; // birth of the ROOT's component
    let mut active = vec![false; n];
    let mut bars = Vec::new();
    for &i in &order {
        let v = field.values[i as usize];
        active[i as usize] = true;
        birth[i as usize] = v;
        let x = i % nx;
        let y = (i / nx) % ny;
        let z = i / (nx * ny);
        let neighbors: [Option<u32>; 6] = [
            (x > 0).then(|| i - 1),
            (x + 1 < nx).then(|| i + 1),
            (y > 0).then(|| i - nx),
            (y + 1 < ny).then(|| i + nx),
            (z > 0).then(|| i - nx * ny),
            (z + 1 < nz).then(|| i + nx * ny),
        ];
        for nb in neighbors.into_iter().flatten() {
            if !active[nb as usize] {
                continue;
            }
            let ri = uf.find(i);
            let rn = uf.find(nb);
            if ri == rn {
                continue;
            }
            // Elder rule: the younger birth dies now.
            let (bi, bn) = (birth[ri as usize], birth[rn as usize]);
            let (survivor_birth, dying_birth) = if bi <= bn { (bi, bn) } else { (bn, bi) };
            if v > dying_birth {
                bars.push(Bar {
                    birth: dying_birth,
                    death: v,
                });
            }
            uf.union(ri, rn);
            let root = uf.find(ri);
            birth[root as usize] = survivor_birth;
        }
    }
    // Essential components.
    for i in 0..n as u32 {
        if uf.find(i) == i && active[i as usize] {
            bars.push(Bar {
                birth: birth[i as usize],
                death: f64::INFINITY,
            });
        }
    }
    bars.sort_by(|a, b| {
        a.birth
            .partial_cmp(&b.birth)
            .expect("finite births")
            .then(a.death.partial_cmp(&b.death).expect("ordered"))
    });
    bars
}

/// Count features whose persistence exceeds `tau` (noise-robust
/// component counting — the ASCENT consumption hook).
#[must_use]
pub fn count_persistent(bars: &[Bar], tau: f64) -> usize {
    bars.iter().filter(|b| b.persistence() > tau).count()
}

/// Verified-at-resolution Betti numbers of a chart's region: voxelize
/// at `n` and compute. HONEST framing: this is exact for the VOXEL
/// solid; features thinner than the cell size can be missed, so the
/// result is an Estimate-grade refinement of `topology_hint()`, not a
/// proof (interval-certified topology is the sheaf bead's).
///
/// # Errors
/// [`VoxelizeError`] under the same conditions as [`voxelize`].
pub fn verify_topology(
    chart: &dyn Chart,
    n: u32,
    cx: &Cx<'_>,
) -> Result<(u32, u32, u32), VoxelizeError> {
    let field = voxelize(chart, n, cx)?;
    Ok(betti(&field, 0.0))
}

/// Verified-at-resolution Betti numbers for the geometric intersection
/// `chart ∩ clip`.
///
/// # Errors
/// [`VoxelizeError`] under the same conditions as [`voxelize_clipped`].
pub fn verify_topology_clipped(
    chart: &dyn Chart,
    n: u32,
    clip: Aabb,
    cx: &Cx<'_>,
) -> Result<(u32, u32, u32), VoxelizeError> {
    let field = voxelize_clipped(chart, n, clip, cx)?;
    Ok(betti(&field, 0.0))
}
