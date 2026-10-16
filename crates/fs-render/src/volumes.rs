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

/// Bit-semantics version for scalar-transfer emission/absorption rendering.
///
/// Bump this when transfer interpolation, stream keying, ray placement, or
/// collision-estimator arithmetic changes in a way that can move image bits.
pub const TRANSFER_RENDER_SEMANTICS_VERSION: u32 = 1;

/// One knot in a scalar-to-optical-properties transfer function.
///
/// `source_rgb` is linear source radiance at a real extinction event, not an
/// opacity or a display-encoded color.  A homogeneous segment therefore has
/// expected radiance `source_rgb * (1 - exp(-extinction * length))`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransferPoint {
    /// Scalar-field value at this knot.
    pub scalar: f64,
    /// Extinction coefficient, in inverse world-length units.
    pub extinction: f64,
    /// Linear-RGB source radiance associated with a real collision.
    pub source_rgb: [f64; 3],
}

/// Interpolated optical properties returned by [`TransferFunction::sample`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransferSample {
    /// Extinction coefficient, in inverse world-length units.
    pub extinction: f64,
    /// Linear-RGB source radiance associated with a real collision.
    pub source_rgb: [f64; 3],
}

/// A validated, piecewise-linear scalar transfer function.
///
/// Values outside the knot domain clamp to the nearest endpoint.  Knot
/// scalars are strictly increasing, so interpolation and tie behavior are
/// deterministic.
#[derive(Debug, Clone, PartialEq)]
pub struct TransferFunction {
    points: Vec<TransferPoint>,
}

/// Fail-closed diagnostics for transfer construction and bounded rendering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DvrError {
    /// At least one transfer knot is required.
    EmptyTransferFunction,
    /// A transfer knot scalar is not finite.
    NonFiniteTransferScalar {
        /// Index of the invalid knot.
        index: usize,
    },
    /// Adjacent knot scalars are not strictly increasing.
    NonIncreasingTransferScalar {
        /// Index of the latter knot.
        index: usize,
    },
    /// An adjacent scalar interval overflows finite arithmetic.
    UnrepresentableTransferInterval {
        /// Index of the latter knot.
        index: usize,
    },
    /// A transfer extinction is negative or non-finite.
    InvalidTransferExtinction {
        /// Index of the invalid knot.
        index: usize,
    },
    /// A transfer source-radiance channel is negative or non-finite.
    InvalidTransferSource {
        /// Index of the invalid knot.
        index: usize,
        /// RGB channel index.
        channel: usize,
    },
    /// A queried scalar is not finite.
    NonFiniteScalarSample,
    /// Grid dimensions, length, origin, cell size, or bounds are invalid.
    InvalidGrid,
    /// A field element is not finite.
    NonFiniteGridValue {
        /// Linear x-major field index.
        index: usize,
    },
    /// Ray coordinates or parameter bounds are invalid.
    InvalidRay,
    /// The Woodcock majorant is negative or non-finite.
    InvalidMajorant,
    /// A sampled extinction exceeded the caller-supplied majorant.
    MajorantViolation,
    /// Image width or height is zero or their product overflowed.
    InvalidImageSize,
    /// Samples per pixel or collision steps per sample is zero.
    InvalidSamplingBudget,
    /// The requested image exceeds its explicit pixel budget.
    PixelBudgetExceeded {
        /// Requested pixels.
        requested: usize,
        /// Admitted pixels.
        limit: usize,
    },
    /// The borrowed field exceeds its explicit majorant-scan budget.
    GridCellBudgetExceeded {
        /// Borrowed field elements.
        requested: usize,
        /// Admitted field elements.
        limit: usize,
    },
    /// The requested primary samples exceed their explicit budget.
    SampleBudgetExceeded {
        /// Requested primary samples.
        requested: u64,
        /// Admitted primary samples.
        limit: u64,
    },
    /// Pixel-count to primary-sample-count arithmetic overflowed.
    SampleCountOverflow,
    /// A pixel identity cannot be represented by the keyed stream.
    PixelIdentityOverflow,
    /// The image buffer could not be reserved in full.
    ImageAllocationRefused {
        /// Requested pixels.
        pixels: usize,
    },
    /// Null-collision tracking exhausted its per-sample hard limit.
    TrackingStepBudgetExceeded,
}

impl TransferFunction {
    /// Validate and take ownership of a transfer-function knot vector.
    pub fn new(points: Vec<TransferPoint>) -> Result<Self, DvrError> {
        if points.is_empty() {
            return Err(DvrError::EmptyTransferFunction);
        }
        for (index, point) in points.iter().enumerate() {
            if !point.scalar.is_finite() {
                return Err(DvrError::NonFiniteTransferScalar { index });
            }
            if !point.extinction.is_finite() || point.extinction < 0.0 {
                return Err(DvrError::InvalidTransferExtinction { index });
            }
            for (channel, value) in point.source_rgb.iter().enumerate() {
                if !value.is_finite() || *value < 0.0 {
                    return Err(DvrError::InvalidTransferSource { index, channel });
                }
            }
            if index > 0 {
                let span = point.scalar - points[index - 1].scalar;
                if span <= 0.0 {
                    return Err(DvrError::NonIncreasingTransferScalar { index });
                }
                if !span.is_finite() {
                    return Err(DvrError::UnrepresentableTransferInterval { index });
                }
            }
        }
        Ok(Self { points })
    }

    /// Return the validated knots in deterministic scalar order.
    #[must_use]
    pub fn points(&self) -> &[TransferPoint] {
        &self.points
    }

    /// Sample the transfer function, clamping outside its knot domain.
    pub fn sample(&self, scalar: f64) -> Result<TransferSample, DvrError> {
        if !scalar.is_finite() {
            return Err(DvrError::NonFiniteScalarSample);
        }
        Ok(self.sample_finite(scalar))
    }

    fn sample_finite(&self, scalar: f64) -> TransferSample {
        let upper = self.points.partition_point(|point| point.scalar <= scalar);
        let point = if upper == 0 {
            self.points[0]
        } else if upper == self.points.len() {
            self.points[self.points.len() - 1]
        } else {
            let lo = self.points[upper - 1];
            let hi = self.points[upper];
            let weight = (scalar - lo.scalar) / (hi.scalar - lo.scalar);
            let lerp = |a: f64, b: f64| (b - a).mul_add(weight, a);
            TransferPoint {
                scalar,
                extinction: lerp(lo.extinction, hi.extinction),
                source_rgb: [
                    lerp(lo.source_rgb[0], hi.source_rgb[0]),
                    lerp(lo.source_rgb[1], hi.source_rgb[1]),
                    lerp(lo.source_rgb[2], hi.source_rgb[2]),
                ],
            }
        };
        TransferSample {
            extinction: point.extinction,
            source_rgb: point.source_rgb,
        }
    }

    /// Maximum mapped extinction over the grid's piecewise-constant cells.
    ///
    /// This admits `max_grid_cells`, then validates the borrowed field and
    /// geometry before scanning it, so the returned value is a conservative
    /// global Woodcock majorant.
    pub fn extinction_majorant(
        &self,
        grid: &VolumeGrid<'_>,
        max_grid_cells: usize,
    ) -> Result<f64, DvrError> {
        if grid.data.len() > max_grid_cells {
            return Err(DvrError::GridCellBudgetExceeded {
                requested: grid.data.len(),
                limit: max_grid_cells,
            });
        }
        validate_transfer_grid(grid)?;
        let mut majorant = 0.0f64;
        for &scalar in grid.data {
            majorant = majorant.max(self.sample_finite(scalar).extinction);
        }
        Ok(majorant)
    }
}

/// Admission settings for an orthographic transfer-function render.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DvrSettings {
    /// `[width, height]` in pixels.
    pub resolution: [usize; 2],
    /// Independent Woodcock estimates averaged per pixel.
    pub samples_per_pixel: u32,
    /// Deterministic image seed.
    pub seed: u64,
    /// Hard image-allocation budget.
    pub max_pixels: usize,
    /// Hard field elements scanned to derive the mapped majorant.
    pub max_grid_cells: usize,
    /// Hard primary-sample budget across the image.
    pub max_samples: u64,
    /// Hard null-collision steps per primary sample.
    pub max_tracking_steps_per_sample: u32,
}

/// A borrowed piecewise-constant scalar field on a regular grid — ZERO-COPY
/// by construction: the buffer belongs to whoever simulates.  Legacy
/// Woodcock entry points interpret the scalar as density/extinction directly;
/// transfer rendering maps it through [`TransferFunction`].
pub struct VolumeGrid<'a> {
    /// Cells per axis.
    pub dims: [usize; 3],
    /// Cell-centered scalar values (x-major: `i + nx·(j + ny·k)`).
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
        self.value_at(p).unwrap_or(0.0)
    }

    /// Scalar value at world point `p`, or `None` outside the grid.
    ///
    /// Unlike [`Self::sigma_at`], this distinguishes an out-of-domain point
    /// from an in-domain zero.  Transfer functions need that distinction
    /// because scalar zero need not map to vacuum.
    #[must_use]
    pub fn value_at(&self, p: [f64; 3]) -> Option<f64> {
        if p.iter().any(|value| !value.is_finite()) {
            return None;
        }
        let mut idx = [0usize; 3];
        for a in 0..3 {
            let u = (p[a] - self.origin[a]) / self.cell[a];
            if u < 0.0 {
                return None;
            }
            let i = u as usize;
            if i >= self.dims[a] {
                return None;
            }
            idx[a] = i;
        }
        Some(self.data[idx[0] + self.dims[0] * (idx[1] + self.dims[1] * idx[2])])
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

fn validate_transfer_grid(grid: &VolumeGrid<'_>) -> Result<(), DvrError> {
    if grid.dims.contains(&0)
        || grid
            .dims
            .iter()
            .try_fold(1usize, |count, &dim| count.checked_mul(dim))
            != Some(grid.data.len())
        || grid.origin.iter().any(|value| !value.is_finite())
        || grid
            .cell
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(DvrError::InvalidGrid);
    }
    let (_, hi) = grid.bounds();
    if hi.iter().zip(grid.origin).any(|(upper, lower)| {
        !upper.is_finite() || *upper <= lower || !(*upper - lower).is_finite()
    }) {
        return Err(DvrError::InvalidGrid);
    }
    if let Some(index) = grid.data.iter().position(|value| !value.is_finite()) {
        return Err(DvrError::NonFiniteGridValue { index });
    }
    Ok(())
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

/// One bounded Woodcock estimate of transfer-function emission/absorption.
///
/// The transfer maps each in-grid scalar to extinction and linear-RGB source
/// radiance.  The supplied global `majorant_bound` is checked against every
/// cell before tracking; an insufficient bound fails instead of biasing the
/// estimate.  `max_tracking_steps` is a hard work limit, and exhaustion
/// returns no radiance claim.
pub fn woodcock_transfer_emission(
    grid: &VolumeGrid<'_>,
    transfer: &TransferFunction,
    majorant_bound: f64,
    ray: Ray,
    stream: &mut fs_rand::Stream,
    max_tracking_steps: u32,
    max_grid_cells: usize,
) -> Result<[f64; 3], DvrError> {
    if max_tracking_steps == 0 {
        return Err(DvrError::InvalidSamplingBudget);
    }
    validate_transfer_ray(ray)?;
    if !majorant_bound.is_finite() || majorant_bound < 0.0 {
        return Err(DvrError::InvalidMajorant);
    }
    let required_majorant = transfer.extinction_majorant(grid, max_grid_cells)?;
    if majorant_bound < required_majorant {
        return Err(DvrError::MajorantViolation);
    }
    track_transfer_emission(
        grid,
        transfer,
        majorant_bound,
        ray,
        stream,
        max_tracking_steps,
    )
}

fn validate_transfer_ray(ray: Ray) -> Result<(), DvrError> {
    let direction_norm_sq = ray.dir.iter().map(|value| value * value).sum::<f64>();
    if ray.origin.iter().any(|value| !value.is_finite())
        || ray.dir.iter().any(|value| !value.is_finite())
        || !ray.t0.is_finite()
        || !ray.t1.is_finite()
        || ray.t1 < ray.t0
        || !direction_norm_sq.is_finite()
        || direction_norm_sq <= 0.0
    {
        return Err(DvrError::InvalidRay);
    }
    Ok(())
}

fn track_transfer_emission(
    grid: &VolumeGrid<'_>,
    transfer: &TransferFunction,
    majorant_bound: f64,
    ray: Ray,
    stream: &mut fs_rand::Stream,
    max_tracking_steps: u32,
) -> Result<[f64; 3], DvrError> {
    if majorant_bound == 0.0 || ray.t0 == ray.t1 {
        return Ok([0.0; 3]);
    }
    let mut t = ray.t0;
    for _ in 0..max_tracking_steps {
        let u = stream.next_f64().max(f64::MIN_POSITIVE);
        t -= u.ln() / majorant_bound;
        if t >= ray.t1 {
            return Ok([0.0; 3]);
        }
        let Some(scalar) = grid.value_at(ray.at(t)) else {
            continue;
        };
        let optical = transfer.sample_finite(scalar);
        if optical.extinction > majorant_bound {
            return Err(DvrError::MajorantViolation);
        }
        if stream.next_f64() * majorant_bound < optical.extinction {
            return Ok(optical.source_rgb);
        }
    }
    Err(DvrError::TrackingStepBudgetExceeded)
}

/// Planck spectral radiance (unnormalized shape) at wavelength
/// `lambda_nm` and temperature `t_kelvin` — the blackbody weight for
/// emissive media. Uses the standard c₂ = hc/k in nm·K.
#[must_use]
pub fn planck(lambda_nm: f64, t_kelvin: f64) -> f64 {
    const C2_NM_K: f64 = 1.438_776_877e7;
    let x = C2_NM_K / (lambda_nm * t_kelvin);
    let l5 = fs_math::det::powi(lambda_nm, 5);
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

/// Render a bounded orthographic linear-RGB emission/absorption image.
///
/// Rays travel along `-z` through the full grid slab.  Each pixel owns a
/// Philox stream keyed independently from the transmittance renderer, so
/// scheduling and tile order cannot alter its samples.  Admission checks the
/// field-scan, image, primary-sample, and per-sample tracking budgets before
/// allocation or rendering.  Any later tracking refusal drops the private
/// staging buffer and returns an error; callers never receive a partial image.
pub fn render_transfer_emission(
    grid: &VolumeGrid<'_>,
    transfer: &TransferFunction,
    settings: DvrSettings,
) -> Result<Vec<[f64; 3]>, DvrError> {
    let [width, height] = settings.resolution;
    let pixels = width
        .checked_mul(height)
        .filter(|_| width > 0 && height > 0)
        .ok_or(DvrError::InvalidImageSize)?;
    if settings.samples_per_pixel == 0 || settings.max_tracking_steps_per_sample == 0 {
        return Err(DvrError::InvalidSamplingBudget);
    }
    if pixels > settings.max_pixels {
        return Err(DvrError::PixelBudgetExceeded {
            requested: pixels,
            limit: settings.max_pixels,
        });
    }
    let pixel_count = u64::try_from(pixels).map_err(|_| DvrError::SampleCountOverflow)?;
    let samples = pixel_count
        .checked_mul(u64::from(settings.samples_per_pixel))
        .ok_or(DvrError::SampleCountOverflow)?;
    if samples > settings.max_samples {
        return Err(DvrError::SampleBudgetExceeded {
            requested: samples,
            limit: settings.max_samples,
        });
    }
    if u32::try_from(pixels - 1).is_err() {
        return Err(DvrError::PixelIdentityOverflow);
    }

    let majorant = transfer.extinction_majorant(grid, settings.max_grid_cells)?;
    let (lo, hi) = grid.bounds();
    let ray_length = hi[2] - lo[2];
    let mut image = Vec::new();
    image
        .try_reserve_exact(pixels)
        .map_err(|_| DvrError::ImageAllocationRefused { pixels })?;
    image.resize(pixels, [0.0; 3]);
    if majorant == 0.0 {
        return Ok(image);
    }

    let sample_scale = 1.0 / f64::from(settings.samples_per_pixel);
    for py in 0..height {
        for px in 0..width {
            let pixel_index = py * width + px;
            let pixel = u32::try_from(pixel_index).map_err(|_| DvrError::PixelIdentityOverflow)?;
            let mut stream = StreamKey {
                seed: settings.seed,
                kernel: 0x0302,
                tile: pixel,
            }
            .stream();
            let x = (hi[0] - lo[0]).mul_add((px as f64 + 0.5) / width as f64, lo[0]);
            let y = (hi[1] - lo[1]).mul_add((py as f64 + 0.5) / height as f64, lo[1]);
            let ray = Ray {
                origin: [x, y, hi[2]],
                dir: [0.0, 0.0, -1.0],
                t0: 0.0,
                t1: ray_length,
            };
            let mut sum = [0.0; 3];
            for _ in 0..settings.samples_per_pixel {
                let sample = track_transfer_emission(
                    grid,
                    transfer,
                    majorant,
                    ray,
                    &mut stream,
                    settings.max_tracking_steps_per_sample,
                )?;
                sum[0] += sample[0];
                sum[1] += sample[1];
                sum[2] += sample[2];
            }
            image[pixel_index] = [
                sum[0] * sample_scale,
                sum[1] * sample_scale,
                sum[2] * sample_scale,
            ];
        }
    }
    Ok(image)
}
