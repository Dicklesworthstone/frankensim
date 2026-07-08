//! Volumetric rendering (bead qfx.3): heterogeneous media by WOODCOCK
//! (delta/null-collision) tracking — the unbiased workhorse pointed
//! directly at live simulation fields. [`VolumeGrid`] BORROWS its
//! density buffer (zero-copy: a running LBM simulation's field renders
//! without a byte moving); majorants come either as a global bound or
//! as a tiled per-block maximum grid (the FrankenVDB tile-maxima
//! wiring is a recorded successor — no fvdb crate in-workspace yet).
//! Beer–Lambert is the analytic homogeneous fast path; HG and Rayleigh
//! phase functions sample exactly; emission uses the collision
//! estimator E[B·1_real] = ∫ σ T B ds with Planck spectral weights.
//! All sampling is per-stream Philox (fs-rand) — images replay
//! bitwise, tile-order independent.

use fs_rand::StreamKey;

/// A borrowed piecewise-constant density field on a regular grid —
/// ZERO-COPY by construction: the buffer belongs to whoever simulates.
pub struct VolumeGrid<'a> {
    /// Cells per axis.
    pub dims: [usize; 3],
    /// Cell-centered densities (x-major: `i + nx·(j + ny·k)`).
    pub data: &'a [f64],
    /// World origin of the grid's min corner.
    pub origin: [f64; 3],
    /// Cell size per axis.
    pub cell: [f64; 3],
}

impl<'a> VolumeGrid<'a> {
    /// Wrap a borrowed buffer.
    ///
    /// # Panics
    /// If `data.len() != nx·ny·nz`.
    #[must_use]
    pub fn new(
        dims: [usize; 3],
        data: &'a [f64],
        origin: [f64; 3],
        cell: [f64; 3],
    ) -> VolumeGrid<'a> {
        assert_eq!(
            data.len(),
            dims[0] * dims[1] * dims[2],
            "field buffer length must match dims"
        );
        VolumeGrid {
            dims,
            data,
            origin,
            cell,
        }
    }

    /// Density at world point `p` (nearest cell; zero outside).
    #[must_use]
    pub fn sigma_at(&self, p: [f64; 3]) -> f64 {
        let mut idx = [0usize; 3];
        for a in 0..3 {
            let u = (p[a] - self.origin[a]) / self.cell[a];
            if u < 0.0 {
                return 0.0;
            }
            let i = u as usize;
            if i >= self.dims[a] {
                return 0.0;
            }
            idx[a] = i;
        }
        self.data[idx[0] + self.dims[0] * (idx[1] + self.dims[1] * idx[2])]
    }

    /// Global majorant (max over cells; ≥ every σ by construction).
    #[must_use]
    pub fn global_majorant(&self) -> f64 {
        self.data.iter().fold(0.0f64, |m, &v| m.max(v))
    }

    /// World-space bounds.
    #[must_use]
    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        let hi = [
            (self.dims[0] as f64).mul_add(self.cell[0], self.origin[0]),
            (self.dims[1] as f64).mul_add(self.cell[1], self.origin[1]),
            (self.dims[2] as f64).mul_add(self.cell[2], self.origin[2]),
        ];
        (self.origin, hi)
    }
}

/// A tiled majorant: per-block maxima over `block³` cells — the
/// free hierarchy sparse fields give (tile maxima), built here
/// deterministically from the dense grid.
pub struct MajorantGrid {
    /// Blocks per axis.
    pub bdims: [usize; 3],
    /// Cells per block edge.
    pub block: usize,
    /// Per-block maxima.
    pub maxima: Vec<f64>,
}

impl MajorantGrid {
    /// Build from a grid.
    #[must_use]
    pub fn build(grid: &VolumeGrid<'_>, block: usize) -> MajorantGrid {
        let block = block.max(1);
        let bdims = [
            grid.dims[0].div_ceil(block),
            grid.dims[1].div_ceil(block),
            grid.dims[2].div_ceil(block),
        ];
        let mut maxima = vec![0.0f64; bdims[0] * bdims[1] * bdims[2]];
        for k in 0..grid.dims[2] {
            for j in 0..grid.dims[1] {
                for i in 0..grid.dims[0] {
                    let v = grid.data[i + grid.dims[0] * (j + grid.dims[1] * k)];
                    let b = (i / block) + bdims[0] * ((j / block) + bdims[1] * (k / block));
                    maxima[b] = maxima[b].max(v);
                }
            }
        }
        MajorantGrid {
            bdims,
            block,
            maxima,
        }
    }

    /// Majorant at world point `p` for the given grid geometry (global
    /// max outside — conservative).
    #[must_use]
    pub fn at(&self, grid: &VolumeGrid<'_>, p: [f64; 3]) -> f64 {
        let mut idx = [0usize; 3];
        for a in 0..3 {
            let u = (p[a] - grid.origin[a]) / grid.cell[a];
            if u < 0.0 {
                return 0.0;
            }
            let i = u as usize;
            if i >= grid.dims[a] {
                return 0.0;
            }
            idx[a] = i / self.block;
        }
        self.maxima[idx[0] + self.bdims[0] * (idx[1] + self.bdims[1] * idx[2])]
    }
}

/// A ray segment: `origin + t·dir` for `t ∈ [t0, t1]`.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin.
    pub origin: [f64; 3],
    /// Direction (unit for metric lengths).
    pub dir: [f64; 3],
    /// Segment start.
    pub t0: f64,
    /// Segment end.
    pub t1: f64,
}

impl Ray {
    /// Point at parameter `t`.
    #[must_use]
    pub fn at(&self, t: f64) -> [f64; 3] {
        [
            self.dir[0].mul_add(t, self.origin[0]),
            self.dir[1].mul_add(t, self.origin[1]),
            self.dir[2].mul_add(t, self.origin[2]),
        ]
    }
}

/// Beer–Lambert transmittance (the homogeneous analytic fast path).
#[must_use]
pub fn beer_lambert(sigma: f64, length: f64) -> f64 {
    (-sigma * length).exp()
}

/// One Woodcock (delta-tracking) transmittance sample along
/// `origin + t·dir`, `t ∈ [t0, t1]`: unbiased {0, 1} estimator whose
/// mean is exp(−∫σ) for ANY `majorant_bound ≥ max σ` — looseness
/// changes only the null-collision count (the battery's unbiasedness
/// gate). `majorant_at` (tile maxima) is a LOOKUP-THINNING stage: a
/// candidate above the local tile bound skips the field lookup
/// entirely; the per-tile-rate DDA traversal is a recorded successor.
#[must_use]
pub fn woodcock_transmittance(
    grid: &VolumeGrid<'_>,
    majorant_at: &dyn Fn([f64; 3]) -> f64,
    majorant_bound: f64,
    ray: Ray,
    stream: &mut fs_rand::Stream,
) -> (f64, u32) {
    let mut t = ray.t0;
    let mut nulls = 0u32;
    if majorant_bound <= 0.0 {
        return (1.0, 0);
    }
    loop {
        // Free flight against the GLOBAL bound (a constant-rate
        // Poisson process thinned twice: once against the local tile
        // majorant, once against the true σ — both thinnings keep the
        // estimator unbiased and the tile stage is where loose global
        // bounds stop costing σ-lookups).
        let u = stream.next_f64().max(f64::MIN_POSITIVE);
        t -= u.ln() / majorant_bound;
        if t >= ray.t1 {
            return (1.0, nulls);
        }
        let p = ray.at(t);
        let m_local = majorant_at(p).min(majorant_bound);
        let v = stream.next_f64() * majorant_bound;
        if v < m_local {
            // Candidate real collision inside the tile bound.
            if v < grid.sigma_at(p) {
                return (0.0, nulls);
            }
            nulls += 1;
        } else {
            nulls += 1;
        }
    }
}

/// Collision-based emission estimator along a ray segment: returns
/// `source(x)` at a REAL collision, 0 on escape — its mean is
/// ∫ σ(s)·T(s)·source(s) ds (for a constant source B and constant σ:
/// B·(1 − exp(−σL)), the closed form the battery gates).
#[must_use]
pub fn woodcock_emission(
    grid: &VolumeGrid<'_>,
    majorant_bound: f64,
    source: &dyn Fn([f64; 3]) -> f64,
    ray: Ray,
    stream: &mut fs_rand::Stream,
) -> f64 {
    if majorant_bound <= 0.0 {
        return 0.0;
    }
    let mut t = ray.t0;
    loop {
        let u = stream.next_f64().max(f64::MIN_POSITIVE);
        t -= u.ln() / majorant_bound;
        if t >= ray.t1 {
            return 0.0;
        }
        let p = ray.at(t);
        if stream.next_f64() * majorant_bound < grid.sigma_at(p) {
            return source(p);
        }
    }
}

/// Planck spectral radiance (unnormalized shape) at wavelength
/// `lambda_nm` and temperature `t_kelvin` — the blackbody weight for
/// emissive media. Uses the standard c₂ = hc/k in nm·K.
#[must_use]
pub fn planck(lambda_nm: f64, t_kelvin: f64) -> f64 {
    const C2_NM_K: f64 = 1.438_776_877e7;
    let x = C2_NM_K / (lambda_nm * t_kelvin);
    let l5 = lambda_nm.powi(5);
    1.0 / (l5 * x.exp_m1())
}

/// Sample the Henyey–Greenstein phase function: returns cosθ with
/// pdf ∝ (1−g²)/(1 + g² − 2g·cosθ)^{3/2}; `E[cosθ] = g`.
#[must_use]
pub fn hg_sample_cos(g: f64, u: f64) -> f64 {
    if g.abs() < 1e-6 {
        2.0f64.mul_add(u, -1.0)
    } else {
        let s = (1.0 - g * g) / 2.0f64.mul_add(g * u, 1.0 - g);
        (g * g + 1.0 - s * s) / (2.0 * g)
    }
}

/// HG phase pdf in cosθ (normalized over cosθ ∈ [−1, 1] with the
/// azimuthal 1/2π folded out).
#[must_use]
pub fn hg_pdf_cos(g: f64, cos_theta: f64) -> f64 {
    let denom = (2.0 * g).mul_add(-cos_theta, 1.0 + g * g);
    0.5 * (1.0 - g * g) / denom.powf(1.5)
}

/// Sample the Rayleigh phase function (pdf ∝ 1 + cos²θ): exact
/// inversion via the depressed cubic; `E[cosθ] = 0`, `E[cos²θ] = 2/5`.
#[must_use]
pub fn rayleigh_sample_cos(u: f64) -> f64 {
    // CDF: (3/8)(c + c³/3) + 1/2 = u  ⇒  c³ + 3c − (8u − 4) = 0.
    let rhs = 8.0f64.mul_add(u, -4.0);
    // Cardano (single real root: discriminant > 0 for all u).
    let d = (rhs * rhs / 4.0 + 1.0).sqrt();
    let c1 = (rhs / 2.0 + d).cbrt();
    let c2 = (rhs / 2.0 - d).cbrt();
    (c1 + c2).clamp(-1.0, 1.0)
}

/// A deterministic per-pixel stream for image rendering.
#[must_use]
pub fn pixel_stream(seed: u64, pixel: u32) -> fs_rand::Stream {
    StreamKey {
        seed,
        kernel: 0x0301,
        tile: pixel,
    }
    .stream()
}

/// Orthographic transmittance image of a grid, rays along −z through
/// the full slab: `res × res` pixels over the grid's xy bounds,
/// `spp` Woodcock samples each, tiled majorant. Deterministic:
/// per-pixel streams, so tile ORDER cannot matter (gated).
#[must_use]
pub fn render_transmittance(
    grid: &VolumeGrid<'_>,
    majorant: &MajorantGrid,
    res: usize,
    spp: u32,
    seed: u64,
) -> Vec<f64> {
    let (lo, hi) = grid.bounds();
    let bound = grid.global_majorant();
    let mut img = vec![0.0f64; res * res];
    for py in 0..res {
        for px in 0..res {
            let pixel = u32::try_from(py * res + px).expect("resolution fits u32");
            let mut stream = pixel_stream(seed, pixel);
            let x = lo[0] + (hi[0] - lo[0]) * ((px as f64 + 0.5) / res as f64);
            let y = lo[1] + (hi[1] - lo[1]) * ((py as f64 + 0.5) / res as f64);
            let ray = Ray {
                origin: [x, y, hi[2] + 1.0],
                dir: [0.0, 0.0, -1.0],
                t0: 1.0,
                t1: 1.0 + (hi[2] - lo[2]),
            };
            let mut acc = 0.0;
            for _ in 0..spp {
                let m_at = |p: [f64; 3]| majorant.at(grid, p);
                let (tr, _) = woodcock_transmittance(grid, &m_at, bound, ray, &mut stream);
                acc += tr;
            }
            img[py * res + px] = acc / f64::from(spp);
        }
    }
    img
}
