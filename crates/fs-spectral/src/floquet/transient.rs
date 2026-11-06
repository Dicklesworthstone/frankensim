//! Finite-horizon singular growth of a complete monodromy at one fixed section.
//!
//! These are powers of the supplied map in its disclosed dimensionless Euclidean
//! state chart. They are not physical energy, a phase quotient, within-period
//! growth, a bound between requested horizons, or a nonlinear basin estimate.

use crate::dense::{DenseSpectrumError, eigenpair_row_residuals};
use fs_la::{
    eigen_complex::{ComplexSvdError, svd_complex_with_cancel},
    gemm::{GemmScalar, gemm_scalar_checked},
};
use fs_math::{c64::C64, dd::Dd, det};

/// TOTAL limits for one complete set of requested integer-period horizons.
#[derive(Clone, Copy, Debug)]
pub struct StroboscopicGrowthControls {
    /// Additionally bounded by the native SVD dimension limit of 1024.
    pub maximum_dimension: usize,
    /// Additionally bounded by 4096 returned horizons.
    pub maximum_horizons: usize,
    pub maximum_periods: u64,
    /// All new live numerical arrays and returned directions are included.
    /// Borrowed inputs, allocator metadata and caller storage are excluded.
    pub maximum_payload_bytes: usize,
    /// Conservative TOTAL arithmetic admission score, not measured FLOPs.
    pub maximum_work: u64,
    pub svd_max_sweeps: usize,
    pub svd_tolerance: f64,
    /// Original-row forward/adjoint residual and unit-direction norm allowance.
    pub maximum_relative_residual: f64,
}

/// Largest singular gain and its actual initial/final directions at one horizon.
#[derive(Clone, Debug)]
pub struct StroboscopicGrowthSample {
    pub periods: u64,
    /// Some(0) only for a zero computed power. For a nonzero power, None means
    /// the positive ordinary gain overflows or underflows; consult log_gain.
    pub gain: Option<f64>,
    /// Natural logarithm of the gain, finite for every nonzero computed power.
    /// None denotes zero gain, without using a nonfinite logarithm sentinel.
    pub log_gain: Option<f64>,
    /// Complete unit right singular direction in the original supplied chart.
    pub initial_direction: Vec<C64>,
    /// Complete unit left singular direction, with the SAME complex phase as
    /// initial_direction. None for zero gain: its final perturbation is zero.
    pub final_direction: Option<Vec<C64>>,
    /// Original-row backward residual of A*v=sigma*u for the normalized power.
    pub forward_relative_residual: f64,
    /// Original-row backward residual of A^H*u=sigma*v, without Gram squaring.
    pub adjoint_relative_residual: f64,
    pub svd_sweeps: usize,
}

/// Numerical evidence at exactly the requested horizons, in their input order.
#[derive(Clone, Debug)]
pub struct StroboscopicGrowthReport {
    pub dimension: usize,
    pub samples: Vec<StroboscopicGrowthSample>,
    /// Number of binary matrix products, including exact-zero fast products.
    pub matrix_products: usize,
    pub admitted_payload_bytes: usize,
    pub admitted_work: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum StroboscopicGrowthError {
    Invalid(&'static str),
    Budget {
        quantity: &'static str,
        requested: u128,
        cap: u128,
    },
    SingularVectors(ComplexSvdError),
    Arithmetic(&'static str),
    DirectionNorm {
        periods: u64,
        defect: f64,
    },
    Residual {
        periods: u64,
        forward: f64,
        adjoint: f64,
    },
    Cancelled,
}
impl core::fmt::Display for StroboscopicGrowthError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "stroboscopic growth refused: {self:?}")
    }
}
impl std::error::Error for StroboscopicGrowthError {}
type Result<T> = core::result::Result<T, StroboscopicGrowthError>;

/// Measure ||M^k||_2 and optimal COMPLETE state perturbations at selected k.
///
/// Counts must be nonempty, sorted, unique nonnegative integers. Zero requests
/// the identity map. No mode, phase direction or defective Jordan component is
/// removed, and no eigendecomposition is required. The original matrix is a
/// complete monodromy supplied by its physical/integration owner.
///
/// A positive global normalization and a double-double accumulated logarithmic
/// scale keep binary powers in range. Rescaling never silently erases a nonzero
/// component; unrepresentable matrix arithmetic refuses. Nonzero gains outside
/// f64's ordinary range remain finite log_gain values, with gain=None.
///
/// Each gap g from the preceding horizon uses popcount(g)+floor(log2(g)) matrix
/// products (zero products for g=0), through the shared scalar-generic GEMM.
/// Only a constant number of full matrices is live. The shared direct complex
/// SVD supplies the largest singular pair; independent original-row checks of
/// [0 A; A^H 0][u;v]=sigma[u;v] verify both equations without forming A^H*A.
/// Singular-pair residuals concern the COMPUTED power, not accumulated powering
/// error, an exact norm enclosure or convergence of the original monodromy.
///
/// With dimension n, H horizons, S native sweeps and B binary products, the
/// admitted work is 128*B*n^3 + H*(128*(S+1)*n^3+1024*n^2+512*n)
/// +64*n^2+256*H. The payload bound is
/// 16*(16*n^2+(4*H+32)*n)+H*sizeof(sample)+4096 bytes, including the native SVD,
/// augmented residual overlap and every returned direction. All limits are
/// admitted before allocation. Cancellation is polled per matrix row, binary
/// product, native SVD boundary and residual row; no partial report escapes.
///
/// # Errors
/// Invalid complete data/controls, TOTAL resource exhaustion, unrepresentable
/// arithmetic, native SVD refusal, failed independent residuals or cancellation.
pub fn stroboscopic_growth(
    matrix: &[C64],
    dimension: usize,
    periods: &[u64],
    controls: StroboscopicGrowthControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<StroboscopicGrowthReport> {
    poll(&mut cancelled)?;
    let (bytes, work) = admit(matrix, dimension, periods, controls, &mut cancelled)?;
    let n = dimension;
    let original = Power::from_matrix(matrix, n, &mut cancelled)?;
    let mut current = Power::identity(n);
    let mut binary = Power::identity(n);
    let mut scratch = vec![Scalar(C64::ZERO); n * n];
    let mut samples = Vec::with_capacity(periods.len());
    let mut previous = 0;
    let mut products = 0;
    for &horizon in periods {
        poll(&mut cancelled)?;
        let mut gap = horizon - previous;
        if gap != 0 {
            for (out, input) in binary
                .values
                .chunks_exact_mut(n)
                .zip(original.values.chunks_exact(n))
            {
                poll(&mut cancelled)?;
                out.copy_from_slice(input);
            }
            binary.log_scale = original.log_scale;
            binary.zero = original.zero;
        }
        while gap != 0 {
            poll(&mut cancelled)?;
            if gap & 1 != 0 {
                let (log_scale, zero) =
                    multiply(&current, &binary, &mut scratch, n, &mut cancelled)?;
                core::mem::swap(&mut current.values, &mut scratch);
                current.log_scale = log_scale;
                current.zero = zero;
                products += 1;
            }
            gap >>= 1;
            if gap != 0 {
                let (log_scale, zero) =
                    multiply(&binary, &binary, &mut scratch, n, &mut cancelled)?;
                core::mem::swap(&mut binary.values, &mut scratch);
                binary.log_scale = log_scale;
                binary.zero = zero;
                products += 1;
            }
        }
        samples.push(assess(&current, n, horizon, controls, &mut cancelled)?);
        previous = horizon;
    }
    poll(&mut cancelled)?;
    Ok(StroboscopicGrowthReport {
        dimension: n,
        samples,
        matrix_products: products,
        admitted_payload_bytes: bytes,
        admitted_work: work,
    })
}

// Only the scalar adapter is local; every matrix product remains in fs-la.
#[derive(Clone, Copy)]
#[repr(transparent)]
struct Scalar(C64);
impl GemmScalar for Scalar {
    fn zero() -> Self {
        Self(C64::ZERO)
    }
    fn is_exact_zero(self) -> bool {
        self.0 == C64::ZERO
    }
    fn mul_add(self, multiplier: Self, addend: Self) -> Self {
        let product = self.0 * multiplier.0;
        // A product of two nonzero complex numbers cannot be zero. Refuse a
        // total product underflow instead of inventing an exact nilpotent map.
        if self.0 != C64::ZERO && multiplier.0 != C64::ZERO && product == C64::ZERO {
            Self(C64::new(f64::NAN, f64::NAN))
        } else {
            Self(product + addend.0)
        }
    }
}

struct Power {
    values: Vec<Scalar>,
    log_scale: Dd,
    zero: bool,
}
impl Power {
    fn from_matrix(matrix: &[C64], n: usize, cancelled: &mut impl FnMut() -> bool) -> Result<Self> {
        let mut values = Vec::with_capacity(matrix.len());
        for row in matrix.chunks_exact(n) {
            poll(cancelled)?;
            values.extend(row.iter().copied().map(Scalar));
        }
        let logarithm = normalize(&mut values, n, cancelled)?;
        Ok(Self {
            values,
            log_scale: Dd::from_f64(logarithm.unwrap_or(0.0)),
            zero: logarithm.is_none(),
        })
    }
    fn identity(n: usize) -> Self {
        let mut values = vec![Scalar(C64::ZERO); n * n];
        for i in 0..n {
            values[i * n + i] = Scalar(C64::ONE);
        }
        Self {
            values,
            log_scale: Dd::ZERO,
            zero: false,
        }
    }
}

fn normalize(
    values: &mut [Scalar],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<f64>> {
    let mut scale = 0.0_f64;
    for row in values.chunks_exact(n) {
        poll(cancelled)?;
        for &Scalar(z) in row {
            if !finite(z) {
                return Err(arithmetic("matrix power is not representable"));
            }
            scale = scale.max(z.re.abs()).max(z.im.abs());
        }
    }
    if scale == 0.0 {
        return Ok(None);
    }
    for row in values.chunks_exact_mut(n) {
        poll(cancelled)?;
        for value in row {
            let z = value.0;
            let normalized = C64::new(z.re / scale, z.im / scale);
            if !finite(normalized)
                || (z.re != 0.0 && normalized.re == 0.0)
                || (z.im != 0.0 && normalized.im == 0.0)
            {
                return Err(arithmetic(
                    "matrix rescaling would erase a nonzero component",
                ));
            }
            value.0 = normalized;
        }
    }
    let logarithm = det::ln(scale);
    if !logarithm.is_finite() {
        return Err(arithmetic("matrix scale logarithm is not representable"));
    }
    Ok(Some(logarithm))
}

fn multiply(
    left: &Power,
    right: &Power,
    output: &mut [Scalar],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Dd, bool)> {
    for row in 0..n {
        poll(cancelled)?;
        let out = &mut output[row * n..(row + 1) * n];
        if left.zero || right.zero {
            out.fill(Scalar(C64::ZERO));
        } else {
            gemm_scalar_checked(
                1,
                n,
                n,
                Scalar(C64::ONE),
                &left.values[row * n..(row + 1) * n],
                &right.values,
                Scalar(C64::ZERO),
                out,
            )
            .map_err(|_| arithmetic("complete binary matrix product shape required"))?;
        }
    }
    let Some(logarithm) = normalize(output, n, cancelled)? else {
        return Ok((Dd::ZERO, true));
    };
    let scale = left.log_scale + right.log_scale + Dd::from_f64(logarithm);
    if !scale.to_f64().is_finite() {
        return Err(arithmetic(
            "accumulated power logarithm is not representable",
        ));
    }
    Ok((scale, false))
}

fn assess(
    power: &Power,
    n: usize,
    periods: u64,
    controls: StroboscopicGrowthControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StroboscopicGrowthSample> {
    let mut matrix = Vec::with_capacity(n * n);
    for row in power.values.chunks_exact(n) {
        poll(cancelled)?;
        matrix.extend(row.iter().map(|value| value.0));
    }
    let svd = svd_complex_with_cancel(
        &matrix,
        n,
        n,
        controls.svd_tolerance,
        controls.svd_max_sweeps,
        &mut *cancelled,
    )
    .map_err(|error| {
        if error == ComplexSvdError::Cancelled {
            StroboscopicGrowthError::Cancelled
        } else {
            StroboscopicGrowthError::SingularVectors(error)
        }
    })?;
    if svd.singular_values.len() != n
        || svd.u.len() != n * n
        || svd.v.len() != n * n
        || svd
            .singular_values
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        || svd.singular_values.windows(2).any(|pair| pair[0] < pair[1])
        || svd.u.iter().chain(&svd.v).any(|&z| !finite(z))
    {
        return Err(arithmetic(
            "complete finite descending native singular decomposition required",
        ));
    }
    let sigma = svd.singular_values[0];
    if (sigma == 0.0) != power.zero {
        return Err(arithmetic(
            "zero singular gain does not match the complete matrix power",
        ));
    }
    let mut initial = Vec::with_capacity(n);
    let mut final_direction = Vec::with_capacity(n);
    let (mut initial_norm, mut final_norm) = (0.0, 0.0);
    for row in 0..n {
        poll(cancelled)?;
        initial.push(svd.v[row * n]);
        final_direction.push(svd.u[row * n]);
        initial_norm = det::hypot(initial_norm, svd.v[row * n].abs());
        final_norm = det::hypot(final_norm, svd.u[row * n].abs());
    }
    let defect = (initial_norm - 1.0).abs().max((final_norm - 1.0).abs());
    if !defect.is_finite() || defect > controls.maximum_relative_residual {
        return Err(StroboscopicGrowthError::DirectionNorm { periods, defect });
    }
    let size = 2 * n;
    let mut augmented = vec![C64::ZERO; size * size];
    for row in 0..n {
        poll(cancelled)?;
        for column in 0..n {
            let z = matrix[row * n + column];
            augmented[row * size + n + column] = z;
            augmented[(n + column) * size + row] = z.conj();
        }
    }
    let mut pair = Vec::with_capacity(size);
    pair.extend_from_slice(&final_direction);
    pair.extend_from_slice(&initial);
    let residuals = eigenpair_row_residuals(
        &augmented,
        size,
        C64::from_re(sigma),
        &pair,
        &mut *cancelled,
    )
    .map_err(|error| match error {
        DenseSpectrumError::Cancelled => StroboscopicGrowthError::Cancelled,
        DenseSpectrumError::Invalid { what } => arithmetic(what),
        _ => arithmetic("original singular-pair residual evaluation failed"),
    })?;
    let forward = residuals[..n].iter().copied().fold(0.0, f64::max);
    let adjoint = residuals[n..].iter().copied().fold(0.0, f64::max);
    if forward > controls.maximum_relative_residual || adjoint > controls.maximum_relative_residual
    {
        return Err(StroboscopicGrowthError::Residual {
            periods,
            forward,
            adjoint,
        });
    }
    let (gain, log_gain) = if power.zero {
        (Some(0.0), None)
    } else {
        let logarithm = (power.log_scale + Dd::from_f64(det::ln(sigma))).to_f64();
        if !logarithm.is_finite() {
            return Err(arithmetic("singular gain logarithm is not representable"));
        }
        let ordinary = det::exp(logarithm);
        (
            ordinary
                .is_finite()
                .then_some(ordinary)
                .filter(|&value| value > 0.0),
            Some(logarithm),
        )
    };
    poll(cancelled)?;
    Ok(StroboscopicGrowthSample {
        periods,
        gain,
        log_gain,
        initial_direction: initial,
        final_direction: (!power.zero).then_some(final_direction),
        forward_relative_residual: forward,
        adjoint_relative_residual: adjoint,
        svd_sweeps: svd.sweeps,
    })
}

fn admit(
    matrix: &[C64],
    n: usize,
    periods: &[u64],
    controls: StroboscopicGrowthControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(usize, u64)> {
    if n == 0
        || n > 1024
        || n.checked_mul(n) != Some(matrix.len())
        || periods.is_empty()
        || periods.len() > 4096
        || controls.svd_max_sweeps == 0
        || !controls.svd_tolerance.is_finite()
        || controls.svd_tolerance <= f64::EPSILON
        || controls.svd_tolerance >= 1.0
        || !controls.maximum_relative_residual.is_finite()
        || controls.maximum_relative_residual <= 0.0
        || controls.maximum_relative_residual >= 1.0
    {
        return Err(StroboscopicGrowthError::Invalid(
            "complete bounded matrix, nonempty horizons and positive native tolerances required",
        ));
    }
    bound(
        "state dimension",
        n as u128,
        controls.maximum_dimension as u128,
    )?;
    bound(
        "requested horizons",
        periods.len() as u128,
        controls.maximum_horizons as u128,
    )?;
    let mut prior = None;
    let mut products = 0_u128;
    for &horizon in periods {
        poll(cancelled)?;
        if prior.is_some_and(|old| horizon <= old) {
            return Err(StroboscopicGrowthError::Invalid(
                "sorted unique nonnegative period counts required",
            ));
        }
        bound(
            "period horizon",
            horizon as u128,
            controls.maximum_periods as u128,
        )?;
        let gap = horizon - prior.unwrap_or(0);
        if gap > 0 {
            products += gap.count_ones() as u128 + (63 - gap.leading_zeros()) as u128;
        }
        prior = Some(horizon);
    }
    for row in matrix.chunks_exact(n) {
        poll(cancelled)?;
        if row.iter().any(|&z| !finite(z)) {
            return Err(StroboscopicGrowthError::Invalid(
                "finite original monodromy entries required",
            ));
        }
    }
    let (d, h, sweeps) = (
        n as u128,
        periods.len() as u128,
        controls.svd_max_sweeps as u128,
    );
    let square = d * d;
    let cube = square * d;
    let per_horizon = checked_add(
        checked_mul(checked_mul(128, checked_add(sweeps, 1)?)?, cube)?,
        1024 * square + 512 * d,
    )?;
    let work = checked_add(
        checked_add(
            checked_mul(128 * cube, products)?,
            checked_mul(h, per_horizon)?,
        )?,
        64 * square + 256 * h,
    )?;
    let bytes = checked_add(
        checked_mul(16, 16 * square + (4 * h + 32) * d)?,
        checked_add(
            checked_mul(h, core::mem::size_of::<StroboscopicGrowthSample>() as u128)?,
            4096,
        )?,
    )?;
    bound("total numerical work", work, controls.maximum_work as u128)?;
    bound(
        "total payload bytes",
        bytes,
        controls.maximum_payload_bytes as u128,
    )?;
    Ok((bytes as usize, work as u64))
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn arithmetic(message: &'static str) -> StroboscopicGrowthError {
    StroboscopicGrowthError::Arithmetic(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(StroboscopicGrowthError::Cancelled)
    } else {
        Ok(())
    }
}
fn bound(quantity: &'static str, requested: u128, cap: u128) -> Result<()> {
    if requested > cap {
        Err(StroboscopicGrowthError::Budget {
            quantity,
            requested,
            cap,
        })
    } else {
        Ok(())
    }
}
fn checked_add(a: u128, b: u128) -> Result<u128> {
    a.checked_add(b)
        .ok_or(arithmetic("total resource envelope overflows"))
}
fn checked_mul(a: u128, b: u128) -> Result<u128> {
    a.checked_mul(b)
        .ok_or(arithmetic("total resource envelope overflows"))
}
