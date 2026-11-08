//! Preintegrated emission and extinction of a continuous scalar field.
//!
//! The scalar varies linearly along each segment. Splitting at every crossed
//! transfer knot makes extinction and source radiance linear on each piece;
//! its optical depth is then analytic. Source transport integrates the same
//! `sigma * source * exp(-optical_depth)` model as the cell-centered renderer.
//! It does not interpolate already-classified opacity or miss a narrow
//! transfer-function peak between two transparent scalar samples.
//!
//! See Engel, Kraus and Ertl (2001), *High-Quality Pre-Integrated Volume
//! Rendering*, section 3, <https://doi.org/10.1145/383507.383515>. Here source
//! radiance is multiplied by extinction, matching [`super::TransferPoint`].
//! There is no scattering solve, refraction, display gamma or measured-optics
//! claim. Quadrature errors are estimates, not interval certificates.

use super::{DvrError, TransferFunction, TransferSample};
use fs_math::det;

#[cfg(feature = "tracer")]
mod ray;
#[cfg(feature = "tracer")]
pub use ray::{
    NodeRayBudget, PreintegratedRayError, PreintegratedRayGrid, PreintegratedRayPixel,
    PreintegratedRaySettings,
};

/// Accuracy and hard work limits for one preintegration call.
///
/// For image rendering the work limits apply to the complete image, rather
/// than resetting for every pixel or segment. Analytic constant-source and
/// constant-extinction pieces consume no quadrature evaluations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreintegrationSettings {
    /// Absolute error target for the normalized mean source position in each
    /// transfer interval, in `[64 * f64::EPSILON, 1]`. The resulting radiance
    /// error estimate is returned by [`PreintegratedPixel::source_error_rgb`].
    pub source_tolerance: f64,
    /// Maximum processed transfer intervals, including vacuum intervals.
    pub max_intervals: usize,
    /// Maximum integrand evaluations across adaptive 8/16-point Gauss rules.
    pub max_quadrature_evaluations: u64,
}

/// Borrowed node-centered scalar data with trilinear reconstruction.
///
/// Dimensions count nodes, not cells; every axis must have at least two nodes.
/// The occupied extent is `(dims - 1) * spacing`. Extinction uses inverse
/// coordinate units, so the caller must use the same length unit for spacing
/// and the transfer function. No density calibration is inferred.
pub struct NodeVolumeGrid<'a> {
    /// Nodes per axis.
    pub dims: [usize; 3],
    /// Node values, indexed by `x + nx * (y + ny * z)`.
    pub data: &'a [f64],
    /// World position of the minimum-corner node.
    pub origin: [f64; 3],
    /// Positive node spacing in world-coordinate units.
    pub spacing: [f64; 3],
}

/// Allocation and work admission for an orthographic node-centered render.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreintegratedImageSettings {
    /// `[width, height]`; rays pass through pixel centers and face `-z`.
    pub resolution: [usize; 2],
    /// Maximum output pixels.
    pub max_pixels: usize,
    /// Maximum borrowed scalar nodes inspected during validation.
    pub max_grid_nodes: usize,
    /// Maximum ray segments, equal to `width * height * (nz - 1)`.
    pub max_segments: u64,
    /// Source accuracy and whole-image transfer/quadrature work limits.
    pub integration: PreintegrationSettings,
}

/// Premultiplied radiance and physical opacity of an integrated segment/ray.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreintegratedPixel {
    emission_rgb: [f64; 3],
    optical_depth: f64,
    source_error_rgb: [f64; 3],
}

impl PreintegratedPixel {
    /// Accumulated linear source radiance; already weighted by absorption.
    #[must_use]
    pub const fn emission_rgb(self) -> [f64; 3] {
        self.emission_rgb
    }

    /// Dimensionless line integral of the mapped extinction coefficient.
    #[must_use]
    pub const fn optical_depth(self) -> f64 {
        self.optical_depth
    }

    /// Fraction of a direct background beam transmitted along this ray.
    #[must_use]
    pub fn transmittance(self) -> f64 {
        det::exp(-self.optical_depth)
    }

    /// Extinguished fraction, without thin-layer cancellation in `1 - T`.
    #[must_use]
    pub fn opacity(self) -> f64 {
        -det::expm1(-self.optical_depth)
    }

    /// Accumulated attenuated quadrature error estimates, in radiance units.
    ///
    /// These exclude scalar reconstruction, transfer-data and roundoff errors;
    /// they are not certified bounds. Optical depth itself is analytic for the
    /// declared piecewise-linear reconstruction.
    #[must_use]
    pub const fn source_error_rgb(self) -> [f64; 3] {
        self.source_error_rgb
    }

    /// Composite over nonnegative linear background radiance as `E + T * B`.
    pub fn composite(self, background_rgb: [f64; 3]) -> Result<[f64; 3], DvrError> {
        if background_rgb.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(DvrError::InvalidBackgroundRadiance);
        }
        let result = std::array::from_fn(|channel| {
            self.emission_rgb[channel]
                + attenuated_product(background_rgb[channel], 1.0, self.optical_depth)
        });
        if result.iter().any(|v| !v.is_finite()) {
            return Err(DvrError::UnrepresentableRadiance);
        }
        Ok(result)
    }
}

/// Refusals leave no partially rendered image in caller-owned storage.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PreintegrationError {
    /// An existing transfer/image admission invariant was violated.
    Transfer(DvrError),
    /// Source accuracy is non-finite or outside its supported range.
    InvalidSourceTolerance,
    /// Endpoints, their difference, or the nonnegative length are non-finite.
    InvalidSegment,
    /// A nonzero interval length or optical depth exceeded binary64 range.
    UnrepresentableOpticalDepth,
    /// Processing another transfer interval would exceed the explicit limit.
    TransferIntervalBudgetExceeded,
    /// Another quadrature rule would exceed the explicit evaluation limit.
    QuadratureBudgetExceeded,
    /// Adaptive integration exhausted its bounded binary subdivision depth.
    QuadratureResolutionLimit,
    /// A cancellable continuous-ray operation was stopped at a work boundary.
    Cancelled,
}

impl From<DvrError> for PreintegrationError {
    fn from(error: DvrError) -> Self {
        Self::Transfer(error)
    }
}

/// Integrate a linear scalar segment, ordered from observer to background.
///
/// Knot crossings are integrated even if both endpoint samples map to vacuum.
/// Reversing endpoints preserves opacity but can change the attenuated color.
/// The segment length is physical world length, not an opacity sampling step.
/// Work uses fixed-size stack storage; it does not allocate a lookup table.
pub fn preintegrate_transfer_segment(
    transfer: &TransferFunction,
    scalar_front: f64,
    scalar_back: f64,
    length: f64,
    settings: PreintegrationSettings,
) -> Result<PreintegratedPixel, PreintegrationError> {
    let mut work = Work::new(settings)?;
    let mut ray = RayIntegral::default();
    integrate_segment(
        transfer,
        scalar_front,
        scalar_back,
        length,
        &mut work,
        &mut ray,
    )?;
    Ok(ray.finish())
}

/// Render borrowed node-centered data with continuous trilinear reconstruction.
///
/// Bilinear scalar interpolation locates each pixel's sample on every z plane;
/// interpolation along each resulting ray segment is linear and is integrated
/// before compositing. This avoids treating nodes as cells, inventing extra
/// thickness, or interpolating classified RGBA instead of the physical field.
/// Work is bounded by the declared node/segment/transfer/quadrature budgets.
/// This scalar CPU path has no per-call cancellation hook or lateral antialiasing.
pub fn render_preintegrated_transfer(
    grid: &NodeVolumeGrid<'_>,
    transfer: &TransferFunction,
    settings: PreintegratedImageSettings,
) -> Result<Vec<PreintegratedPixel>, PreintegrationError> {
    let mut work = Work::new(settings.integration)?;
    let pixels = admit_image(grid, settings)?;
    let mut image = Vec::new();
    image
        .try_reserve_exact(pixels)
        .map_err(|_| DvrError::ImageAllocationRefused { pixels })?;
    let [width, height] = settings.resolution;
    for py in 0..height {
        let (iy, fy) = pixel_node_coordinate(py, height, grid.dims[1]);
        for px in 0..width {
            let (ix, fx) = pixel_node_coordinate(px, width, grid.dims[0]);
            let mut front = grid.plane_scalar(ix, iy, grid.dims[2] - 1, fx, fy);
            let mut ray = RayIntegral::default();
            for z in (0..grid.dims[2] - 1).rev() {
                let back = grid.plane_scalar(ix, iy, z, fx, fy);
                integrate_segment(transfer, front, back, grid.spacing[2], &mut work, &mut ray)?;
                front = back;
            }
            image.push(ray.finish());
        }
    }
    Ok(image)
}

impl NodeVolumeGrid<'_> {
    fn plane_scalar(&self, x: usize, y: usize, z: usize, fx: f64, fy: f64) -> f64 {
        let base = x + self.dims[0] * (y + self.dims[1] * z);
        let lower = lerp(self.data[base], self.data[base + 1], fx);
        let upper = lerp(
            self.data[base + self.dims[0]],
            self.data[base + self.dims[0] + 1],
            fx,
        );
        lerp(lower, upper, fy)
    }
}

fn admit_image(
    grid: &NodeVolumeGrid<'_>,
    settings: PreintegratedImageSettings,
) -> Result<usize, DvrError> {
    let [width, height] = settings.resolution;
    let pixels = width
        .checked_mul(height)
        .filter(|_| width > 0 && height > 0)
        .ok_or(DvrError::InvalidImageSize)?;
    if pixels > settings.max_pixels {
        return Err(DvrError::PixelBudgetExceeded {
            requested: pixels,
            limit: settings.max_pixels,
        });
    }
    if grid.data.len() > settings.max_grid_nodes {
        return Err(DvrError::GridCellBudgetExceeded {
            requested: grid.data.len(),
            limit: settings.max_grid_nodes,
        });
    }
    if grid.dims.iter().any(|n| *n < 2)
        || grid.dims.iter().try_fold(1usize, |n, &d| n.checked_mul(d)) != Some(grid.data.len())
    {
        return Err(DvrError::InvalidGrid);
    }
    let segments = u64::try_from(pixels)
        .ok()
        .zip(u64::try_from(grid.dims[2] - 1).ok())
        .and_then(|(n, z)| n.checked_mul(z))
        .ok_or(DvrError::CellVisitCountOverflow)?;
    if segments > settings.max_segments {
        return Err(DvrError::CellVisitBudgetExceeded {
            requested: segments,
            limit: settings.max_segments,
        });
    }
    for axis in 0..3 {
        let length = (grid.dims[axis] - 1) as f64 * grid.spacing[axis];
        let upper = grid.origin[axis] + length;
        if !grid.origin[axis].is_finite()
            || !grid.spacing[axis].is_finite()
            || grid.spacing[axis] <= 0.0
            || !length.is_finite()
            || length <= 0.0
            || !upper.is_finite()
            || upper <= grid.origin[axis]
        {
            return Err(DvrError::InvalidGrid);
        }
    }
    if let Some(index) = grid.data.iter().position(|v| !v.is_finite()) {
        return Err(DvrError::NonFiniteGridValue { index });
    }
    Ok(pixels)
}

// Exact integer cell ownership avoids subtracting a large translated origin.
// Valid dimensions match a real f64 slice, bounding this u128 product even on
// 64-bit targets; pixel allocation is admitted before these coordinates run.
fn pixel_node_coordinate(pixel: usize, resolution: usize, nodes: usize) -> (usize, f64) {
    let numerator = (2 * pixel as u128 + 1) * (nodes as u128 - 1);
    let denominator = 2 * resolution as u128;
    (
        (numerator / denominator) as usize,
        (numerator % denominator) as f64 / denominator as f64,
    )
}

struct Work<'a> {
    settings: PreintegrationSettings,
    intervals: usize,
    evaluations: u64,
    checkpoint: Option<&'a mut dyn FnMut() -> Result<(), PreintegrationError>>,
}

impl<'a> Work<'a> {
    fn new(settings: PreintegrationSettings) -> Result<Self, PreintegrationError> {
        if !settings.source_tolerance.is_finite()
            || !(64.0 * f64::EPSILON..=1.0).contains(&settings.source_tolerance)
        {
            return Err(PreintegrationError::InvalidSourceTolerance);
        }
        Ok(Self {
            settings,
            intervals: 0,
            evaluations: 0,
            checkpoint: None,
        })
    }

    #[cfg(feature = "tracer")]
    fn with_checkpoint(
        settings: PreintegrationSettings,
        checkpoint: &'a mut dyn FnMut() -> Result<(), PreintegrationError>,
    ) -> Result<Self, PreintegrationError> {
        let mut work = Self::new(settings)?;
        work.checkpoint = Some(checkpoint);
        Ok(work)
    }

    fn interval(&mut self) -> Result<(), PreintegrationError> {
        if let Some(check) = self.checkpoint.as_mut() {
            check()?;
        }
        if self.intervals >= self.settings.max_intervals {
            return Err(PreintegrationError::TransferIntervalBudgetExceeded);
        }
        self.intervals += 1;
        Ok(())
    }

    fn quadrature(&mut self) -> Result<(), PreintegrationError> {
        if let Some(check) = self.checkpoint.as_mut() {
            check()?;
        }
        if self.settings.max_quadrature_evaluations - self.evaluations < 24 {
            return Err(PreintegrationError::QuadratureBudgetExceeded);
        }
        self.evaluations += 24;
        Ok(())
    }
}

fn integrate_segment(
    transfer: &TransferFunction,
    front: f64,
    back: f64,
    length: f64,
    work: &mut Work<'_>,
    ray: &mut RayIntegral,
) -> Result<(), PreintegrationError> {
    let span = back - front;
    if !front.is_finite()
        || !back.is_finite()
        || !span.is_finite()
        || !length.is_finite()
        || length < 0.0
    {
        return Err(PreintegrationError::InvalidSegment);
    }
    if length == 0.0 {
        return Ok(());
    }
    if front == back {
        work.interval()?;
        let optical = transfer.sample_finite(front);
        return ray.add_piece(optical, optical, length, work);
    }
    let points = transfer.points();
    let increasing = span > 0.0;
    let mut knot = if increasing {
        points.partition_point(|p| p.scalar <= front)
    } else {
        points.partition_point(|p| p.scalar < front)
    };
    let mut scalar = front;
    let mut optical = transfer.sample_finite(scalar);
    loop {
        let next = if increasing && knot < points.len() && points[knot].scalar < back {
            let value = points[knot].scalar;
            knot += 1;
            value
        } else if !increasing && knot > 0 && points[knot - 1].scalar > back {
            knot -= 1;
            points[knot].scalar
        } else {
            back
        };
        work.interval()?;
        let delta = (next - scalar).abs();
        let mut piece_length = length * (delta / span.abs());
        if piece_length == 0.0 && delta > 0.0 {
            piece_length = det::exp(det::ln(length) + det::ln(delta) - det::ln(span.abs()));
        }
        if piece_length == 0.0 || !piece_length.is_finite() {
            return Err(PreintegrationError::UnrepresentableOpticalDepth);
        }
        let next_optical = transfer.sample_finite(next);
        ray.add_piece(optical, next_optical, piece_length, work)?;
        if next == back {
            return Ok(());
        }
        scalar = next;
        optical = next_optical;
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, value: f64) {
        if self.value.is_infinite() {
            return;
        }
        let corrected = value - self.correction;
        let next = self.value + corrected;
        self.correction = if next.is_finite() {
            (next - self.value) - corrected
        } else {
            0.0
        };
        self.value = next;
    }
}

#[derive(Default)]
struct RayIntegral {
    depth: Sum,
    emission: [Sum; 3],
    error: [Sum; 3],
}

impl RayIntegral {
    fn add_piece(
        &mut self,
        front: TransferSample,
        back: TransferSample,
        length: f64,
        work: &mut Work<'_>,
    ) -> Result<(), PreintegrationError> {
        let low = front.extinction.min(back.extinction);
        let high = front.extinction.max(back.extinction);
        if high == 0.0 {
            return Ok(());
        }
        let mean_extinction = low + 0.5 * (high - low);
        let mut tau = mean_extinction * length;
        if tau == 0.0 {
            // Averaging subnormal coefficients can underflow even when their
            // product with a long segment has a representable optical depth.
            tau = det::exp(det::ln(high) + det::ln(length) + det::ln(0.5 * (1.0 + low / high)));
        }
        if !tau.is_finite() || tau == 0.0 {
            return Err(PreintegrationError::UnrepresentableOpticalDepth);
        }
        let alpha = -det::expm1(-tau);
        let (position, error) = if front.source_rgb == back.source_rgb {
            (0.0, 0.0)
        } else if front.extinction == back.extinction {
            (constant_extinction_position(tau, alpha), 0.0)
        } else {
            let quantile = SourceQuantile {
                tau,
                alpha,
                transmission: det::exp(-tau),
                front: front.extinction / high,
                back: back.extinction / high,
            };
            integrate_mean(&quantile, 0.0, 1.0, work.settings.source_tolerance, 0, work)?
        };
        for channel in 0..3 {
            let source = lerp(
                front.source_rgb[channel],
                back.source_rgb[channel],
                position,
            );
            self.emission[channel].add(attenuated_product(source, alpha, self.depth.value));
            let variation = (back.source_rgb[channel] - front.source_rgb[channel]).abs();
            self.error[channel].add(attenuated_product(
                variation * error,
                alpha,
                self.depth.value,
            ));
            if !self.emission[channel].value.is_finite() || !self.error[channel].value.is_finite() {
                return Err(DvrError::UnrepresentableRadiance.into());
            }
        }
        self.depth.add(tau);
        Ok(())
    }

    fn finish(self) -> PreintegratedPixel {
        PreintegratedPixel {
            emission_rgb: self.emission.map(|s| s.value),
            optical_depth: self.depth.value,
            source_error_rgb: self.error.map(|s| s.value),
        }
    }
}

fn lerp(a: f64, b: f64, weight: f64) -> f64 {
    if a.is_sign_negative() == b.is_sign_negative() {
        (b - a).mul_add(weight, a)
    } else {
        (1.0 - weight).mul_add(a, weight * b)
    }
}

fn attenuated_product(source: f64, weight: f64, depth: f64) -> f64 {
    if source == 0.0 || weight == 0.0 {
        return 0.0;
    }
    let transmission = det::exp(-depth);
    let product = source * weight * transmission;
    if transmission < f64::MIN_POSITIVE || product == 0.0 {
        det::exp(det::ln(source) + det::ln(weight) - depth)
    } else {
        product
    }
}

// Mean position of the first interaction conditional on interaction, for a
// homogeneous slab: 1/tau - 1/expm1(tau). The thin-limit series avoids loss of
// the leading 1/2. Its first omitted term is below 3e-17 at tau = 0.1.
fn constant_extinction_position(tau: f64, alpha: f64) -> f64 {
    if tau < 0.1 {
        let square = tau * tau;
        0.5 - tau
            * (1.0 / 12.0
                - square * (1.0 / 720.0 - square * (1.0 / 30_240.0 - square / 1_209_600.0)))
    } else {
        1.0 / tau - det::exp(-tau) / alpha
    }
}

struct SourceQuantile {
    tau: f64,
    alpha: f64,
    transmission: f64,
    front: f64,
    back: f64,
}

impl SourceQuantile {
    fn at(&self, t: f64) -> f64 {
        // Integrate over cumulative absorption probability rather than distance:
        // y = -log(1 - alpha*p). This resolves thick slabs without marching a
        // tiny mean free path. Smoothstep p(t) regularizes both endpoint limits.
        let p = t * t * (3.0 - 2.0 * t);
        let opposite = 1.0 - t;
        let q = opposite * opposite * (1.0 + 2.0 * t);
        let removed = self.alpha * p;
        let fraction = if removed < 0.5 {
            let remainder = 1.0 - removed;
            let logarithmic_ratio = if remainder == 1.0 {
                1.0
            } else {
                -det::ln(remainder) / (1.0 - remainder)
            };
            // Normalize before multiplying by p so thin subnormal optical
            // depths retain the full conditional interaction distribution.
            logarithmic_ratio * (self.alpha / self.tau) * p
        } else {
            -det::ln(q + p * self.transmission) / self.tau
        };
        let fraction = fraction.clamp(0.0, 1.0);
        // Stable inverse of the quadratic integrated linear extinction. The
        // discriminant is a positive interpolation, even for falling sigma.
        let discriminant =
            (1.0 - fraction) * self.front * self.front + fraction * self.back * self.back;
        let position = if fraction == 0.0 {
            0.0
        } else {
            fraction * (self.front + self.back) / (self.front + discriminant.sqrt())
        };
        position * (6.0 * t * opposite)
    }
}

// Positive half of Gauss-Legendre rules on [-1, 1]. Evaluating both signs
// preserves symmetry and never samples a singular interval endpoint.
const GAUSS_8: [(f64, f64); 4] = [
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];
const GAUSS_16: [(f64, f64); 8] = [
    (0.095_012_509_837_637_44, 0.189_450_610_455_068_5),
    (0.281_603_550_779_258_9, 0.182_603_415_044_923_6),
    (0.458_016_777_657_227_4, 0.169_156_519_395_002_5),
    (0.617_876_244_402_643_8, 0.149_595_988_816_576_7),
    (0.755_404_408_355_003, 0.124_628_971_255_533_9),
    (0.865_631_202_387_831_8, 0.095_158_511_682_492_78),
    (0.944_575_023_073_232_6, 0.062_253_523_938_647_89),
    (0.989_400_934_991_649_9, 0.027_152_459_411_754_1),
];

fn gauss(rule: &[(f64, f64)], f: &SourceQuantile, center: f64, half: f64) -> f64 {
    let mut sum = 0.0;
    for &(x, weight) in rule {
        sum += weight * (f.at(center - half * x) + f.at(center + half * x));
    }
    half * sum
}

fn integrate_mean(
    f: &SourceQuantile,
    left: f64,
    right: f64,
    tolerance: f64,
    depth: u32,
    work: &mut Work<'_>,
) -> Result<(f64, f64), PreintegrationError> {
    work.quadrature()?;
    let half = 0.5 * (right - left);
    let center = left + half;
    let coarse = gauss(&GAUSS_8, f, center, half);
    let fine = gauss(&GAUSS_16, f, center, half);
    let error = (fine - coarse).abs();
    if error <= tolerance {
        return Ok((fine, error));
    }
    // Fixed recursion depth bounds stack use independently of caller budgets.
    if depth == 40 || center == left || center == right {
        return Err(PreintegrationError::QuadratureResolutionLimit);
    }
    let (a, a_error) = integrate_mean(f, left, center, 0.5 * tolerance, depth + 1, work)?;
    let (b, b_error) = integrate_mean(f, center, right, 0.5 * tolerance, depth + 1, work)?;
    Ok((a + b, a_error + b_error))
}

#[cfg(test)]
mod tests;
