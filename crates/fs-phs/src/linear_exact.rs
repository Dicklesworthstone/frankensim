//! Opt-in linear propagation; the reference Gonzalez path is unchanged.
//!
//! Exact ZOH retains the full state and all damping/coupling entries. A Jacobi
//! basis balances the quadratic energy before a dense exponential; it is not
//! modal truncation. Held-force work and loss are independently integrated by
//! a jointly integrated finite-horizon Gramian, not inferred from the energy
//! difference. Stable doubling never forms an exponentially growing inverse
//! propagator for a physically decaying mode.
//!
//! Frequency matching replaces Q by Q + delta_Q in the Gonzalez Hamiltonian.
//! In energy coordinates K = sqrt(Q) J sqrt(Q), each eigenvalue omega^2 of
//! -K^2 gets the factor tan(omega*dt/2)/(omega*dt/2). Thus the lossless linear
//! Cayley map has exactly the original phase below pi. Nonlinear storage is
//! retained verbatim. The ledger is for this MODIFIED Hamiltonian; static
//! compliance and damped rates are not exact in this mode. Use ZOH for held
//! linear input. Neither mode certifies a nonlinear trajectory or a render.

use crate::{PhsError, PortHamiltonian, StepRecord, Storage, discrete_gradient};
use fs_la::eigen::{admit_jacobi_eigh, jacobi_eigh};
use fs_la::factor::lu;
use fs_math::det;

mod gramian;
use gramian::{accumulate, identity, product};
/// Cold caching of the same complete quadratic trajectory and its integrals.
pub mod prepared;

/// Conservative midpoint angular-rate band. For a damped oscillatory pole
/// with |lambda|*dt <= 0.04, Cayley-map phase error is below one cent even in
/// the nearly critically damped limit (the lossless error is smaller).
pub const MIDPOINT_ONE_CENT_ANGULAR_STEP: f64 = 0.04;

/// Dense Van Loan block ceiling, checked before matrix allocation by callers
/// preparing complete coupled render operators. No state truncation is allowed.
pub const MAX_LINEAR_EXPONENTIAL_BLOCK_DIM: usize = 256;

/// A non-exact linear block exceeded its declared resolved temporal band.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ResolvedBandError {
    /// No finite nonnegative rate or positive finite clock was supplied.
    Invalid { rate_rad_s: f64, dt_s: f64 },
    /// Raising the clock is necessary; coefficients are never retuned/clamped.
    Outside {
        /// Supplied upper bound on the linear pole magnitude [1/s].
        rate_rad_s: f64,
        /// Requested timestep [s].
        dt_s: f64,
        /// Maximum admitted rate times timestep [rad].
        maximum_angular_step: f64,
        /// Largest admitted angular-rate bound divided by 2*pi [Hz].
        maximum_frequency_hz: f64,
    },
}

impl core::fmt::Display for ResolvedBandError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { rate_rad_s, dt_s } => write!(
                f,
                "midpoint phase admission requires a finite nonnegative pole-rate bound and positive dt: rate={rate_rad_s}, dt={dt_s}"
            ),
            Self::Outside {
                rate_rad_s,
                dt_s,
                maximum_angular_step,
                maximum_frequency_hz,
            } => write!(
                f,
                "midpoint phase refuses rate*dt={} above {maximum_angular_step}; this clock resolves rate/(2*pi) <= {maximum_frequency_hz} Hz within one cent",
                rate_rad_s * dt_s
            ),
        }
    }
}
impl std::error::Error for ResolvedBandError {}

/// Admit a midpoint/Cayley linear block against a <=1-cent phase budget.
/// Supply a bound on ALL its linear pole magnitudes, including damping and
/// coupled stiffness. A single bare component frequency is not such a bound.
pub fn admit_midpoint_band(rate_rad_s: f64, dt_s: f64) -> Result<(), ResolvedBandError> {
    if !rate_rad_s.is_finite() || rate_rad_s < 0.0 || !dt_s.is_finite() || dt_s <= 0.0 {
        return Err(ResolvedBandError::Invalid { rate_rad_s, dt_s });
    }
    if rate_rad_s * dt_s > MIDPOINT_ONE_CENT_ANGULAR_STEP {
        return Err(ResolvedBandError::Outside {
            rate_rad_s,
            dt_s,
            maximum_angular_step: MIDPOINT_ONE_CENT_ANGULAR_STEP,
            maximum_frequency_hz: MIDPOINT_ONE_CENT_ANGULAR_STEP / (core::f64::consts::TAU * dt_s),
        });
    }
    Ok(())
}

/// Explicit decomposition H(x) = x^T Q x / 2 + H_remainder(x).
#[derive(Clone, Debug, PartialEq)]
pub struct QuadraticPart {
    /// Symmetric PSD Q, row-major, in the owner's unchanged state order.
    pub matrix: Vec<f64>,
    /// True only if the remainder is identically zero, not merely small.
    pub exact: bool,
}

/// The caller chooses the linear propagator explicitly.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LinearPropagator {
    /// Full coupled exponential, with zero-order-held input and integrated work.
    ExactZoh,
    /// Lossless-frequency-matched Gonzalez solve, with modified energy.
    FrequencyMatched,
}

/// Dense preparation capacity and the frequency-matching band.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearSplitBudget {
    /// Maximum scalar basis coordinates, including zero-energy coordinates.
    /// A second-order oscillator uses two; no coordinates are discarded.
    pub max_modes: usize,
    /// Maximum omega*dt for frequency matching; additionally must remain < pi.
    pub max_omega_dt: f64,
}

/// Full-dimensional energy basis and lossless frequencies, bound to J/R/G/Q.
#[derive(Clone, Debug)]
pub struct LinearSplit {
    n: usize,
    m: usize,
    j: Vec<f64>,
    r: Vec<f64>,
    g: Vec<f64>,
    quadratic: QuadraticPart,
    values: Vec<f64>,
    basis: Vec<f64>,
    frequencies: Vec<f64>,
    frequency_basis: Vec<f64>,
    budget: LinearSplitBudget,
    pole_rate_bound: f64,
}

fn refuse(what: &'static str) -> PhsError {
    PhsError::LinearSplitBudget { what }
}

fn finite(values: &[f64]) -> bool {
    values.iter().all(|v| v.is_finite())
}

fn quadratic(a: &[f64], x: &[f64]) -> f64 {
    let n = x.len();
    let mut result = 0.0;
    for i in 0..n {
        for j in 0..n {
            result += x[i] * a[i * n + j] * x[j];
        }
    }
    result
}

// Bounded scaling/squaring with a converged Taylor series at ||A||_inf <= 1/2.
// Dense work is deliberately capped; no sparse, cancellation-latency or
// allocation-free claim. These caps apply before forming the block matrix.
/// Exponential of a finite real square operator with at most 256 rows.
/// Reuses the exact-ZOH owner's deterministic scaling/squaring algorithm.
/// The caller includes the requested duration in `a`.
///
/// # Errors
/// Refuses invalid shape, nonfinite entries, work-cap overflow or failure of
/// the bounded Taylor/scaling iteration. No stability certificate is inferred.
pub fn exponential(a: &[f64], n: usize) -> Result<Vec<f64>, PhsError> {
    exponential_with_budget(a, n, DenseExponentialBudget {
        max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
        max_dense_bytes: usize::MAX,
        max_multiply_terms: u64::MAX,
    })
}

/// Explicit cold dense-exponential admission for owners that prepare and cache
/// a complete operator. This does not enlarge the ordinary split/ZOH capacity.
/// The checked byte count covers input, output and all simultaneous matrix
/// workspace; the term count reserves 32 Taylor products and 48 squarings.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DenseExponentialBudget {
    /// Matrix dimension ceiling, at most MAX_PREPARED_EXPONENTIAL_DIM (6144).
    pub max_dimension: usize,
    /// Maximum simultaneous dense matrix bytes, including the supplied input.
    pub max_dense_bytes: usize,
    /// Maximum scalar multiply-add terms in dense products.
    pub max_multiply_terms: u64,
}

/// Hard capacity of the explicit budgeted cold path; defaults remain 256. 6144 is the 2d x 2d Van Loan operator of a
/// 3072-dimensional generator (fs-couple's 1024 prepared modes, three coordinates each).
pub const MAX_PREPARED_EXPONENTIAL_DIM: usize = 6144;

/// The same deterministic Taylor/scaling exponential under an explicit larger
/// cold work/byte budget. No approximate backend, state removal, stability
/// certificate or hot-path timing guarantee is selected by this entry point.
///
/// # Errors
/// Refuses capacity/workspace/work overflow before allocating a matrix, and
/// the same nonfinite or unconverged operators as [`exponential`].
pub fn exponential_with_budget(
    a: &[f64],
    n: usize,
    budget: DenseExponentialBudget,
) -> Result<Vec<f64>, PhsError> {
    admit_exponential(n, budget)?;
    if n == 0
        || n > budget.max_dimension
        || n.checked_mul(n) != Some(a.len())
        || a.iter().any(|x| !x.is_finite())
    {
        return Err(refuse("invalid or oversized exponential operator"));
    }
    let mut norm = 0.0_f64;
    for row in a.chunks_exact(n) {
        norm = norm.max(row.iter().map(|v| v.abs()).sum());
    }
    if !norm.is_finite() {
        return Err(refuse("nonfinite exponential operator"));
    }
    let mut scale = 1.0;
    let mut squarings = 0;
    while norm * scale > 0.5 {
        if squarings == 48 {
            return Err(refuse("exponential scaling exceeds 48 squarings"));
        }
        scale *= 0.5;
        squarings += 1;
    }
    let scaled: Vec<_> = a.iter().map(|x| x * scale).collect();
    let mut result = identity(n);
    let mut term = result.clone();
    let mut converged = false;
    for order in 1..=32 {
        term = product(&term, &scaled, n);
        let mut largest = 0.0_f64;
        for (sum, value) in result.iter_mut().zip(&mut term) {
            *value /= f64::from(order);
            *sum += *value;
            largest = largest.max(value.abs());
        }
        if largest <= 0.25 * f64::EPSILON {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(refuse("exponential series did not converge"));
    }
    for _ in 0..squarings {
        result = product(&result, &result, n);
    }
    if !finite(&result) {
        return Err(refuse("exponential is unrepresentable"));
    }
    Ok(result)
}

/// The dense work and byte admission of an n x n exponential, before any
/// operator is read or matrix allocated.
fn admit_exponential(n: usize, budget: DenseExponentialBudget) -> Result<(), PhsError> {
    let extent = n.checked_mul(n).ok_or_else(|| refuse("exponential extent overflow"))?;
    let bytes = extent.checked_mul(5).and_then(|x| x.checked_mul(core::mem::size_of::<f64>()))
        .ok_or_else(|| refuse("exponential workspace overflow"))?;
    let terms = (n as u64).checked_pow(3).and_then(|x| x.checked_mul(80))
        .ok_or_else(|| refuse("exponential work overflow"))?;
    if budget.max_dimension == 0 || budget.max_dimension > MAX_PREPARED_EXPONENTIAL_DIM
        || bytes > budget.max_dense_bytes || terms > budget.max_multiply_terms {
        return Err(refuse("exponential dense byte/work admission"));
    }
    Ok(())
}

/// The three blocks of a block upper-triangular 2d x 2d operator or of its
/// exponential, `[[first, cross], [0, last]]`, each d x d and row-major.
#[derive(Clone, Debug, PartialEq)]
pub struct UpperBlocks {
    /// Upper-left block.
    pub first: Vec<f64>,
    /// Upper-right block.
    pub cross: Vec<f64>,
    /// Lower-right block.
    pub last: Vec<f64>,
}

impl UpperBlocks {
    fn entries(&self) -> impl Iterator<Item = &f64> {
        self.first.iter().chain(&self.cross).chain(&self.last)
    }

    /// self x right. Each entry of the 2d x 2d product adds its k terms in the
    /// textbook order; the terms that meet the zero block are exact zeros and
    /// add nothing, so they are not formed: four d x d products instead of
    /// eight.
    fn product(&self, right: &Self, d: usize) -> Self {
        let mut cross = vec![0.0; d * d];
        accumulate(&mut cross, &self.first, &right.cross, d);
        accumulate(&mut cross, &self.cross, &right.last, d);
        Self { first: product(&self.first, &right.first, d), cross, last: product(&self.last, &right.last, d) }
    }
}

/// [`exponential_with_budget`] of the block upper-triangular operator
/// `[[first, cross], [0, last]]` (each block d x d), returned as the three
/// blocks of the exponential. Its lower-left block is exactly zero.
///
/// It is the same computation bit for bit, with the same admission of the
/// 2d x 2d operator and the same refusals. The norm, scaling, Taylor series
/// and squarings are identical; only the products are formed by blocks. That
/// is exact: every product of block upper-triangular matrices keeps the zero
/// block (its entries are sums of exact zeros, which stay +0), and a term
/// that meets it adds an exact zero to a sum that never holds -0. Van Loan's
/// integral operator `[[-A^T h, R h], [0, A h]]` has this form.
///
/// # Errors
/// The refusals of [`exponential_with_budget`] for the assembled operator.
pub fn upper_block_exponential_with_budget(
    operator: &UpperBlocks,
    d: usize,
    budget: DenseExponentialBudget,
) -> Result<UpperBlocks, PhsError> {
    let n = d.checked_mul(2).ok_or_else(|| refuse("exponential extent overflow"))?;
    admit_exponential(n, budget)?;
    if n == 0
        || n > budget.max_dimension
        || [&operator.first, &operator.cross, &operator.last].iter().any(|b| d.checked_mul(d) != Some(b.len()))
        || operator.entries().any(|x| !x.is_finite())
    {
        return Err(refuse("invalid or oversized exponential operator"));
    }
    // The row sums of the assembled operator, each in its own column order.
    let mut norm = 0.0_f64;
    for (first, cross) in operator.first.chunks_exact(d).zip(operator.cross.chunks_exact(d)) {
        norm = norm.max(first.iter().chain(cross).map(|v| v.abs()).sum());
    }
    for last in operator.last.chunks_exact(d) {
        norm = norm.max(core::iter::repeat(&0.0).take(d).chain(last).map(|v: &f64| v.abs()).sum());
    }
    if !norm.is_finite() {
        return Err(refuse("nonfinite exponential operator"));
    }
    let mut scale = 1.0;
    let mut squarings = 0;
    while norm * scale > 0.5 {
        if squarings == 48 {
            return Err(refuse("exponential scaling exceeds 48 squarings"));
        }
        scale *= 0.5;
        squarings += 1;
    }
    let scale_block = |block: &[f64]| -> Vec<f64> { block.iter().map(|x| x * scale).collect() };
    let scaled = UpperBlocks {
        first: scale_block(&operator.first), cross: scale_block(&operator.cross), last: scale_block(&operator.last),
    };
    let mut result = UpperBlocks { first: identity(d), cross: vec![0.0; d * d], last: identity(d) };
    let mut term = result.clone();
    let mut converged = false;
    for order in 1..=32 {
        term = term.product(&scaled, d);
        // Entry by entry, as the assembled loop; the largest term is the same
        // maximum, since the zero block's terms are zero.
        let mut largest = 0.0_f64;
        for (sum, value) in [&mut result.first, &mut result.cross, &mut result.last].into_iter()
            .zip([&mut term.first, &mut term.cross, &mut term.last])
            .flat_map(|(sums, values)| sums.iter_mut().zip(values.iter_mut()))
        {
            *value /= f64::from(order);
            *sum += *value;
            largest = largest.max(value.abs());
        }
        if largest <= 0.25 * f64::EPSILON {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(refuse("exponential series did not converge"));
    }
    for _ in 0..squarings {
        result = result.product(&result, d);
    }
    if result.entries().any(|x| !x.is_finite()) {
        return Err(refuse("exponential is unrepresentable"));
    }
    Ok(result)
}

/// Prepare without probing H or dropping any state, damping, or port terms.
///
/// # Errors
/// Refuses absent explicit quadratic data, capacity, malformed raw systems,
/// nonfinite/negative energy eigenvalues, or an inaccurate eigendecomposition.
pub fn split_quadratic(
    sys: &PortHamiltonian,
    budget: LinearSplitBudget,
) -> Result<LinearSplit, PhsError> {
    let (n, m) = (sys.n, sys.m);
    let extent = n.checked_mul(n).ok_or_else(|| refuse("state extent"))?;
    let block = n
        .checked_add(m)
        .and_then(|x| x.checked_add(1))
        .and_then(|x| x.checked_mul(2))
        .ok_or_else(|| refuse("ZOH block extent"))?;
    if n == 0
        || n > budget.max_modes
        || block > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
        || !budget.max_omega_dt.is_finite()
        || budget.max_omega_dt <= 0.0
    {
        return Err(refuse(
            "state capacity, ZOH block <= 256, and positive finite band required",
        ));
    }
    if sys.j.len() != extent
        || sys.r.len() != extent
        || sys.g.len() != n.checked_mul(m).ok_or_else(|| refuse("port extent"))?
        || !finite(&sys.j)
        || !finite(&sys.r)
        || !finite(&sys.g)
    {
        return Err(refuse("finite J/R/G shapes required"));
    }
    crate::require_skew(&sys.j, n, "split J")?;
    crate::require_symmetric(&sys.r, n, "split R")?;
    admit_jacobi_eigh(n).map_err(|_| refuse("Jacobi work admission"))?;
    crate::require_psd(&sys.r, n, "split R")?;
    let part = sys
        .storage
        .quadratic_part()
        .ok_or_else(|| refuse("storage has no explicit quadratic part"))?;
    if part.matrix.len() != extent || !finite(&part.matrix) {
        return Err(refuse("finite quadratic matrix shape required"));
    }
    crate::require_symmetric(&part.matrix, n, "split Q")?;
    let (mut values, basis) = checked_eigh(&part.matrix, n)?;
    // The relative floor PortHamiltonian admission uses: a zero-energy
    // coordinate (a free mass, a rigid mode) can return from Jacobi as -eps.
    let floor = crate::PSD_TOL * crate::matrix_scale(&part.matrix);
    if values.iter().any(|v| *v < -floor) {
        return Err(PhsError::NotPsd { what: "split Q" });
    }
    for value in &mut values {
        *value = value.max(0.0);
    }
    let mut skew = vec![0.0; extent];
    for i in 0..n {
        for j in 0..n {
            for a in 0..n {
                for b in 0..n {
                    skew[i * n + j] += basis[a * n + i] * sys.j[a * n + b] * basis[b * n + j];
                }
            }
            skew[i * n + j] *= det::sqrt(values[i]) * det::sqrt(values[j]);
        }
    }
    // K^T K = -K^2; forming the Gram matrix preserves exact symmetry.
    let mut squared = vec![0.0; extent];
    for i in 0..n {
        for j in i..n {
            let value = (0..n).map(|k| skew[k * n + i] * skew[k * n + j]).sum();
            squared[i * n + j] = value;
            squared[j * n + i] = value;
        }
    }
    if !finite(&squared) {
        return Err(refuse("frequency operator overflow"));
    }
    let (_, frequency_basis) = checked_eigh(&squared, n)?;
    // Evaluate the Gram Rayleigh quotient as ||K v||². An odd-dimensional
    // skew operator has a zero mode; subtractive Jacobi updates can otherwise
    // report a tiny negative eigenvalue there. This is a sum of squares, not
    // a clamp or a dropped low-frequency coordinate.
    let mut frequencies = Vec::with_capacity(n);
    for k in 0..n {
        let mut omega2 = 0.0;
        for i in 0..n {
            let component: f64 = (0..n)
                .map(|j| skew[i * n + j] * frequency_basis[j * n + k])
                .sum();
            omega2 += component * component;
        }
        if !omega2.is_finite() {
            return Err(refuse("squared frequency overflow"));
        }
        frequencies.push(det::sqrt(omega2));
    }
    // Bound the FULL linear generator in balanced energy coordinates. This
    // includes non-diagonal damping and gyroscopic/interconnection entries.
    let mut infinity_norm = 0.0_f64;
    let mut columns = vec![0.0; n];
    for i in 0..n {
        let wi = if values[i] > 0.0 {
            det::sqrt(values[i])
        } else {
            1.0
        };
        let mut row = 0.0;
        for j in 0..n {
            let wj = if values[j] > 0.0 {
                det::sqrt(values[j])
            } else {
                1.0
            };
            let mut jr = 0.0;
            for a in 0..n {
                for b in 0..n {
                    jr +=
                        basis[a * n + i] * (sys.j[a * n + b] - sys.r[a * n + b]) * basis[b * n + j];
                }
            }
            let entry = (wi * jr * values[j] / wj).abs();
            row += entry;
            columns[j] += entry;
        }
        if !row.is_finite() {
            return Err(refuse("linear pole-rate norm overflow"));
        }
        infinity_norm = infinity_norm.max(row);
    }
    let one_norm = columns.into_iter().fold(0.0_f64, f64::max);
    let pole_rate_bound = det::sqrt(infinity_norm) * det::sqrt(one_norm);
    if !pole_rate_bound.is_finite() {
        return Err(refuse("linear spectral norm bound overflow"));
    }
    Ok(LinearSplit {
        n,
        m,
        j: sys.j.clone(),
        r: sys.r.clone(),
        g: sys.g.clone(),
        quadratic: part,
        values,
        basis,
        frequencies,
        frequency_basis,
        budget,
        pole_rate_bound,
    })
}

fn checked_eigh(a: &[f64], n: usize) -> Result<(Vec<f64>, Vec<f64>), PhsError> {
    let (values, vectors) = jacobi_eigh(a, n);
    let scale = a.iter().fold(f64::MIN_POSITIVE, |s, v| s.max(v.abs()));
    for k in 0..n {
        for i in 0..n {
            let av: f64 = (0..n).map(|j| a[i * n + j] * vectors[j * n + k]).sum();
            if !av.is_finite()
                || (av - values[k] * vectors[i * n + k]).abs() > 1e-12 * n as f64 * scale
            {
                return Err(PhsError::Eigen);
            }
        }
    }
    if !finite(&values) || !finite(&vectors) {
        return Err(PhsError::Eigen);
    }
    Ok((values, vectors))
}

impl LinearSplit {
    /// Conservative upper bound on magnitudes of ALL linear poles [1/s].
    #[must_use]
    pub fn pole_rate_bound_rad_s(&self) -> f64 {
        self.pole_rate_bound
    }
    /// Frequencies of all energy-basis coordinates [rad/s], ascending.
    /// Each oscillatory pair occurs twice; zero modes remain present.
    #[must_use]
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.frequencies
    }

    /// The matched Hamiltonian's quadratic correction in original coordinates.
    /// Its static compliance change must be disclosed by a consumer.
    pub fn quadratic_correction(&self, dt: f64) -> Result<Vec<f64>, PhsError> {
        if !dt.is_finite() || dt < 0.0 {
            return Err(refuse("finite nonnegative dt required"));
        }
        let mut factors = Vec::with_capacity(self.n);
        for &omega in &self.frequencies {
            let z = omega * dt;
            if !z.is_finite() || z >= core::f64::consts::PI || z > self.budget.max_omega_dt {
                return Err(PhsError::FrequencyMatchUnadmitted { omega_dt: z });
            }
            factors.push(if z == 0.0 {
                0.0
            } else {
                det::tan(0.5 * z) / (0.5 * z) - 1.0
            });
        }
        let n = self.n;
        let mut directions = vec![0.0; n * n];
        for i in 0..n {
            for k in 0..n {
                for j in 0..n {
                    directions[i * n + k] += self.basis[i * n + j]
                        * det::sqrt(self.values[j])
                        * self.frequency_basis[j * n + k];
                }
            }
        }
        let mut correction = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let value = (0..n)
                    .map(|k| directions[i * n + k] * factors[k] * directions[j * n + k])
                    .sum();
                correction[i * n + j] = value;
                correction[j * n + i] = value;
            }
        }
        if !finite(&correction) {
            return Err(refuse("matched quadratic energy overflow"));
        }
        Ok(correction)
    }
}

/// Render stepping with exact quadratic propagation and a resolved-band
/// Gonzalez fallback for an explicitly declared nonlinear remainder. No
/// storage probe, modified Hamiltonian, silent fallback or modal truncation.
pub fn step_render(
    sys: &PortHamiltonian,
    split: &LinearSplit,
    x0: &[f64],
    u: &[f64],
    dt: f64,
) -> Result<StepRecord, PhsError> {
    if split.quadratic.exact {
        return step_split(sys, split, LinearPropagator::ExactZoh, x0, u, dt);
    }
    if sys.n != split.n
        || sys.m != split.m
        || sys.j != split.j
        || sys.r != split.r
        || sys.g != split.g
        || sys.storage.quadratic_part().as_ref() != Some(&split.quadratic)
    {
        return Err(refuse("render split does not match J/R/G/storage"));
    }
    admit_midpoint_band(split.pole_rate_bound, dt).map_err(PhsError::TimeResolution)?;
    crate::step(sys, x0, u, dt)
}

/// Upper bound on the linear oscillatory pole magnitudes of an effort-joined
/// descriptor. The algebraic velocity constraint is an orthogonal projection
/// in the quadratic energy coordinates; its norm is one. Thus the full
/// differential generator's spectral-norm bound also bounds its restriction
/// to that constraint. No constraint stiffness or multiplier mass is invented.
pub fn descriptor_pole_rate_bound(sys: &crate::DescriptorPortHamiltonian) -> Result<f64, PhsError> {
    let part = sys
        .storage
        .quadratic_part()
        .ok_or_else(|| refuse("descriptor storage has no explicit quadratic part"))?;
    let n = sys.n_diff;
    if part.matrix.len() != sys.n * sys.n {
        return Err(refuse("descriptor quadratic shape"));
    }
    let mut q = vec![0.0; n * n];
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    for a in 0..n {
        for b in 0..n {
            q[a * n + b] = part.matrix[a * sys.n + b];
            j[a * n + b] = sys.j[a * sys.n + b];
            r[a * n + b] = sys.r[a * sys.n + b];
        }
    }
    let quadratic = crate::QuadraticStorage::new(q, n)?;
    let differential = PortHamiltonian::new(n, 0, j, r, vec![], Box::new(quadratic))?;
    let split = split_quadratic(
        &differential,
        LinearSplitBudget {
            max_modes: n,
            max_omega_dt: core::f64::consts::PI,
        },
    )?;
    // A zero-energy coordinate has no oscillatory restoring force, but can
    // receive motion. The bound above retains its generator entries as well.
    Ok(split.pole_rate_bound)
}

/// Step the complete state. Errors publish no state and never switch propagator.
///
/// ZOH is admitted only for explicitly quadratic storage. Frequency matching
/// keeps the full nonlinear remainder in a strict Gonzalez/Newton solve.
pub fn step_split(
    sys: &PortHamiltonian,
    split: &LinearSplit,
    prop: LinearPropagator,
    x0: &[f64],
    u: &[f64],
    dt: f64,
) -> Result<StepRecord, PhsError> {
    if sys.n != split.n
        || sys.m != split.m
        || sys.j != split.j
        || sys.r != split.r
        || sys.g != split.g
        || sys.storage.quadratic_part().as_ref() != Some(&split.quadratic)
    {
        return Err(refuse("split does not match J/R/G/storage"));
    }
    if x0.len() != split.n
        || u.len() != split.m
        || !finite(x0)
        || !finite(u)
        || !dt.is_finite()
        || dt < 0.0
    {
        return Err(refuse("finite state/input and nonnegative dt required"));
    }
    match prop {
        LinearPropagator::ExactZoh => exact_zoh(sys, split, x0, u, dt),
        LinearPropagator::FrequencyMatched => {
            let correction = split.quadratic_correction(dt)?;
            let storage = MatchedStorage {
                original: sys.storage.as_ref(),
                correction: &correction,
            };
            gonzalez(sys, &storage, x0, u, dt)
        }
    }
}

fn exact_zoh(
    sys: &PortHamiltonian,
    split: &LinearSplit,
    x0: &[f64],
    u: &[f64],
    dt: f64,
) -> Result<StepRecord, PhsError> {
    exact_zoh_loss(sys, split, x0, u, dt, &sys.r, true)
}

/// Integrate one declared PSD resistance channel along the SAME complete
/// coupled trajectory. This does not replace the dynamics' R by that channel.
/// Independent channel integrals let coupled owners retain their original
/// component/connection loss accounts after migrating to a full exponential.
pub fn integrated_resistance_loss(
    sys: &PortHamiltonian,
    split: &LinearSplit,
    x0: &[f64],
    u: &[f64],
    dt: f64,
    resistance: &[f64],
) -> Result<f64, PhsError> {
    let n = split.n;
    if sys.n != n
        || sys.m != split.m
        || sys.j != split.j
        || sys.r != split.r
        || sys.g != split.g
        || sys.storage.quadratic_part().as_ref() != Some(&split.quadratic)
        || x0.len() != n
        || u.len() != split.m
        || !finite(x0)
        || !finite(u)
        || !dt.is_finite()
        || dt < 0.0
        || resistance.len() != n * n
        || !finite(resistance)
    {
        return Err(refuse(
            "finite matching channel-loss split/state/input/clock required",
        ));
    }
    crate::require_symmetric(resistance, n, "loss channel R")?;
    crate::require_psd(resistance, n, "loss channel R")?;
    Ok(exact_zoh_loss(sys, split, x0, u, dt, resistance, false)?.dissipated)
}

fn exact_zoh_loss(
    sys: &PortHamiltonian,
    split: &LinearSplit,
    x0: &[f64],
    u: &[f64],
    dt: f64,
    loss_resistance: &[f64],
    check_balance: bool,
) -> Result<StepRecord, PhsError> {
    if !split.quadratic.exact {
        return Err(refuse("ZOH requires identically quadratic storage"));
    }
    if dt == 0.0 {
        return Ok(StepRecord {
            x: x0.to_vec(),
            y: sys.output(x0),
            delta_h: 0.0,
            supplied: 0.0,
            dissipated: 0.0,
            newton_iters: 0,
            solver_residual: 0.0,
            solver_certificate: None,
        });
    }
    let n = split.n;
    let d = n + 1 + split.m;
    let weights: Vec<_> = split
        .values
        .iter()
        .map(|&v| if v > 0.0 { det::sqrt(v) } else { 1.0 })
        .collect();
    // z = W V^T x, x = V W^-1 z; zero-energy directions retain unit scale.
    let mut forward = vec![0.0; n * n];
    let mut inverse = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            forward[i * n + j] = weights[i] * split.basis[j * n + i];
            inverse[j * n + i] = split.basis[j * n + i] / weights[i];
        }
    }
    let mut flow = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                flow[i * n + j] +=
                    (sys.j[i * n + k] - sys.r[i * n + k]) * split.quadratic.matrix[k * n + j];
            }
        }
    }
    let balanced = product(&product(&forward, &flow, n), &inverse, n);
    let effort_map = product(&split.quadratic.matrix, &inverse, n);
    let mut b = vec![0.0; d * d];
    let mut initial = vec![0.0; d];
    initial[n] = 1.0;
    for i in 0..n {
        initial[i] = (0..n).map(|j| forward[i * n + j] * x0[j]).sum();
        for j in 0..n {
            b[i * d + j] = balanced[i * n + j];
        }
        for a in 0..n {
            for p in 0..split.m {
                b[i * d + n] += forward[i * n + a] * sys.g[a * split.m + p] * u[p];
            }
        }
    }
    // Integrators for EVERY port output, independently of the held input.
    for p in 0..split.m {
        for j in 0..n {
            for a in 0..n {
                b[(n + 1 + p) * d + j] += sys.g[a * split.m + p] * effort_map[a * n + j];
            }
        }
    }
    let mut loss = vec![0.0; d * d];
    let r_effort = product(loss_resistance, &effort_map, n);
    for i in 0..n {
        for j in 0..n {
            loss[i * d + j] = (0..n)
                .map(|a| effort_map[a * n + i] * r_effort[a * n + j])
                .sum();
        }
    }
    // Preserve the original dense capacity gate. Joint Taylor/Lyapunov
    // propagation works on d x d matrices and doubles W as W + E^T W E.
    // No exp(-B^T dt) is formed: stiff physical decay retains finite heat.
    admit_exponential(2 * d, DenseExponentialBudget {
        max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
        max_dense_bytes: usize::MAX,
        max_multiply_terms: u64::MAX,
    })?;
    let flow = gramian::flow_integral(&b, &loss, d, dt).map_err(refuse)?;
    let mut end = vec![0.0; d];
    for i in 0..d {
        for j in 0..d {
            end[i] += flow.transition[i * d + j] * initial[j];
        }
    }
    let dissipated = quadratic(&flow.loss_integral, &initial);
    let x: Vec<_> = (0..n)
        .map(|i| (0..n).map(|j| inverse[i * n + j] * end[j]).sum())
        .collect();
    let y: Vec<_> = end[n + 1..].iter().map(|v| v / dt).collect();
    let supplied = end[n + 1..].iter().zip(u).map(|(a, b)| a * b).sum::<f64>();
    let delta_h = sys.hamiltonian(&x) - sys.hamiltonian(x0);
    let energy_scale = sys
        .hamiltonian(x0)
        .abs()
        .max(sys.hamiltonian(&x).abs())
        .max(dissipated.abs())
        .max(supplied.abs())
        .max(f64::MIN_POSITIVE);
    let defect = delta_h + dissipated - supplied;
    if !finite(&x)
        || !finite(&y)
        || !finite(&[delta_h, dissipated, supplied])
        || dissipated < -1e-12 * energy_scale
        || (check_balance && defect.abs() > 1e-12 * energy_scale)
    {
        return Err(refuse(
            "independent exponential energy balance exceeds 1e-12 relative",
        ));
    }
    Ok(StepRecord {
        x,
        y,
        delta_h,
        dissipated,
        supplied,
        newton_iters: 0,
        solver_residual: 0.0,
        solver_certificate: None,
    })
}

struct MatchedStorage<'a> {
    original: &'a dyn Storage,
    correction: &'a [f64],
}
impl Storage for MatchedStorage<'_> {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.original.hamiltonian(x) + 0.5 * quadratic(self.correction, x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.original.gradient(x, out);
        for i in 0..x.len() {
            for j in 0..x.len() {
                out[i] += self.correction[i * x.len() + j] * x[j];
            }
        }
    }
}

fn gonzalez(
    sys: &PortHamiltonian,
    storage: &dyn Storage,
    x0: &[f64],
    u: &[f64],
    dt: f64,
) -> Result<StepRecord, PhsError> {
    let n = sys.n;
    let forcing: Vec<f64> = (0..n)
        .map(|i| (0..sys.m).map(|p| sys.g[i * sys.m + p] * u[p]).sum())
        .collect();
    let residual = |x: &[f64]| {
        let dg = discrete_gradient(storage, x0, x);
        (0..n)
            .map(|i| {
                x[i] - x0[i]
                    - dt * (forcing[i]
                        + (0..n)
                            .map(|j| (sys.j[i * n + j] - sys.r[i * n + j]) * dg[j])
                            .sum::<f64>())
            })
            .collect::<Vec<_>>()
    };
    let mut x = x0.to_vec();
    for iterations in 0..=50 {
        let res = residual(&x);
        if !finite(&res) || !finite(&x) {
            return Err(refuse("nonfinite Gonzalez trial"));
        }
        let norm = res.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
        let scale = x
            .iter()
            .chain(x0)
            .fold(f64::MIN_POSITIVE, |s, v| s.max(v.abs()))
            .max(forcing.iter().fold(0.0_f64, |s, v| s.max(dt * v.abs())));
        if norm <= 1e-13 * scale {
            let dg = discrete_gradient(storage, x0, &x);
            let dissipated = dt * quadratic(&sys.r, &dg);
            let y = sys.output_from_effort(&dg);
            let supplied = dt * u.iter().zip(&y).map(|(a, b)| a * b).sum::<f64>();
            let delta_h = storage.hamiltonian(&x) - storage.hamiltonian(x0);
            if !finite(&[delta_h, dissipated, supplied]) {
                return Err(refuse("nonfinite matched ledger"));
            }
            return Ok(StepRecord {
                x,
                y,
                dissipated,
                supplied,
                delta_h,
                newton_iters: iterations,
                solver_residual: norm,
                solver_certificate: None,
            });
        }
        if iterations == 50 {
            return Err(PhsError::NewtonStalled { residual: norm });
        }
        let mut jac = vec![0.0; n * n];
        for j in 0..n {
            let h = 1e-5 * x[j].abs().max(scale);
            if h == 0.0 || !h.is_finite() {
                return Err(refuse("Gonzalez difference increment"));
            }
            let mut plus = x.clone();
            let mut minus = x.clone();
            plus[j] += h;
            minus[j] -= h;
            let rp = residual(&plus);
            let rm = residual(&minus);
            for i in 0..n {
                jac[i * n + j] = (rp[i] - rm[i]) / (2.0 * h);
            }
        }
        if !finite(&jac) {
            return Err(refuse("nonfinite Gonzalez Jacobian"));
        }
        let factor = lu(&jac, n).map_err(|_| PhsError::NewtonStalled { residual: norm })?;
        let mut delta = res;
        factor.solve(&mut delta);
        for i in 0..n {
            x[i] -= delta[i];
        }
    }
    unreachable!("bounded Newton loop always accepts or refuses")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{QuadraticStorage, ScalarLaw, SeparableStorage, mass_spring_damper};

    #[test]
    fn explicit_prepared_exponential_preserves_default_and_checks_dense_resources() {
        let n = MAX_LINEAR_EXPONENTIAL_BLOCK_DIM + 1;
        let input = vec![0.0; n * n];
        let budget = DenseExponentialBudget {
            max_dimension: n,
            max_dense_bytes: 5 * n * n * core::mem::size_of::<f64>(),
            max_multiply_terms: 80 * (n as u64).pow(3),
        };
        assert!(exponential(&input, n).is_err());
        assert_eq!(exponential_with_budget(&input, n, budget).unwrap(), identity(n));
        assert!(exponential_with_budget(&input, n, DenseExponentialBudget {
            max_dense_bytes: budget.max_dense_bytes - 1, ..budget }).is_err());
        assert!(exponential_with_budget(&input, n, DenseExponentialBudget {
            max_multiply_terms: budget.max_multiply_terms - 1, ..budget }).is_err());
        assert!(exponential_with_budget(&input, n, DenseExponentialBudget {
            max_dimension: MAX_PREPARED_EXPONENTIAL_DIM + 1, ..budget }).is_err());
        let ordinary = [-0.37, 0.2, -0.1, -0.8];
        assert_eq!(exponential(&ordinary, 2).unwrap(),
            exponential_with_budget(&ordinary, 2, budget).unwrap());
    }

    #[test]
    fn public_exponential_shape_cap_and_scalar_oracle() {
        assert_eq!(
            exponential(&vec![0.0; 256 * 256], 256).unwrap(),
            identity(256)
        );
        assert!(exponential(&[], 257).is_err());
        assert!(exponential(&[], 0).is_err());
        assert!(exponential(&[0.0; 3], 2).is_err());
        assert!(exponential(&[f64::NAN], 1).is_err());
        assert!((exponential(&[-0.37], 1).unwrap()[0] - det::exp(-0.37)).abs() < 2e-15);
    }

    fn budget(n: usize) -> LinearSplitBudget {
        LinearSplitBudget {
            max_modes: n,
            max_omega_dt: core::f64::consts::PI,
        }
    }

    #[test]
    fn midpoint_resolved_band_cap_cap_plus_one_and_damped_poles() {
        let cap = MIDPOINT_ONE_CENT_ANGULAR_STEP;
        assert!(admit_midpoint_band(1.0, cap).is_ok());
        assert!(matches!(
            admit_midpoint_band(1.0, cap + 1e-12),
            Err(ResolvedBandError::Outside { .. })
        ));
        for (rate, dt) in [
            (f64::NAN, 1.0),
            (1.0, 0.0),
            (-1.0, 1.0),
            (1.0, f64::INFINITY),
        ] {
            assert!(matches!(
                admit_midpoint_band(rate, dt),
                Err(ResolvedBandError::Invalid { .. })
            ));
        }
        // Analytic complex Cayley map, including the near-critical limit.
        for zeta in [0.0, 0.3, 0.9, 0.999_999] {
            let beta = det::sqrt(1.0 - zeta * zeta);
            let phase = det::atan2(beta * cap, 1.0 - 0.25 * cap * cap);
            let cents = 1200.0 * det::ln(phase / (beta * cap)) / core::f64::consts::LN_2;
            assert!(cents.abs() <= 1.0, "{zeta}: {cents} cents");
        }
    }

    #[test]
    fn independent_resistance_channels_sum_along_the_same_coupled_trajectory() {
        let q = QuadraticStorage::new(vec![4.0, 0.0, 0.0, 1.0], 2).unwrap();
        let sys = PortHamiltonian::new(
            2,
            1,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.7],
            vec![0.0, 1.0],
            Box::new(q),
        )
        .unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        for dt in [0.05, 0.5, 1.25] {
            let whole = step_render(&sys, &split, &[0.3, 0.4], &[0.2], dt).unwrap();
            let a = integrated_resistance_loss(
                &sys,
                &split,
                &[0.3, 0.4],
                &[0.2],
                dt,
                &[0.0, 0.0, 0.0, 0.2],
            )
            .unwrap();
            let b = integrated_resistance_loss(
                &sys,
                &split,
                &[0.3, 0.4],
                &[0.2],
                dt,
                &[0.0, 0.0, 0.0, 0.5],
            )
            .unwrap();
            assert!(a > 0.0 && b > 0.0);
            assert!((a + b - whole.dissipated).abs() <= 1e-12 * whole.dissipated);
        }
        assert!(
            integrated_resistance_loss(
                &sys,
                &split,
                &[0.3, 0.4],
                &[0.2],
                0.1,
                &[0.0, 0.0, 0.0, -0.1]
            )
            .is_err()
        );
    }

    #[test]
    fn descriptor_band_bounds_the_actual_common_velocity_mode() {
        let a = mass_spring_damper(1.0, 1.0, 0.0).unwrap();
        let b = mass_spring_damper(1.0, 4.0, 0.0).unwrap();
        let sys = crate::join_port(a, b, 0, 0).unwrap();
        let bound = descriptor_pole_rate_bound(&sys).unwrap();
        let dt = MIDPOINT_ONE_CENT_ANGULAR_STEP / bound;
        assert!(admit_midpoint_band(bound, dt).is_ok());
        assert!(admit_midpoint_band(bound, dt * (1.0 + 1e-10)).is_err());
        let omega = det::sqrt(2.5);
        for scale in [0.25, 0.5, 1.0] {
            let h = dt * scale;
            let rec = crate::step_descriptor(&sys, &[1e-4, 0.0, 1e-4, 0.0, 0.0], &[], h).unwrap();
            let q = 0.5 * (rec.x[0] + rec.x[2]);
            let v = 0.5 * (rec.x[1] + rec.x[3]);
            let phase = det::atan2(-v / omega, q);
            let cents = 1200.0 * det::ln(phase / (omega * h)) / core::f64::consts::LN_2;
            assert!(cents.abs() <= 1.0, "{cents}");
        }
    }

    #[test]
    fn quadratic_metadata_survives_relaxation_and_full_galerkin_transport() {
        let base = mass_spring_damper(2.0, 3.0, 0.2).unwrap();
        let expanded = base
            .with_relaxation_branches(vec![crate::RelaxationBranch {
                projection: vec![1.0, 0.0],
                stiffness: 4.0,
                relaxation_time_s: 0.3,
            }])
            .unwrap();
        let part = expanded.storage.quadratic_part().unwrap();
        assert!(part.exact);
        assert_eq!(
            part.matrix,
            vec![7.0, 0.0, -2.0, 0.0, 0.5, 0.0, -2.0, 0.0, 1.0]
        );
        let basis = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0];
        let reduced = crate::reduce_galerkin(&expanded, &basis, 3).unwrap();
        let transported = reduced.storage.quadratic_part().unwrap();
        assert!(transported.exact);
        assert_eq!(
            transported.matrix,
            vec![7.0, -2.0, 0.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.5]
        );
        let split = split_quadratic(&reduced, budget(3)).unwrap();
        let record = step_render(&reduced, &split, &[0.1, 0.02, 0.3], &[0.5], 0.2).unwrap();
        assert!(record.balance_residual().abs() < 1e-12);
        let nonlinear = crate::duffing_oscillator(1.0, 1.0, 2.0, 0.0).unwrap();
        // The public v1 reducer admits quadratic storage only. A quadratic
        // metadata hook does not authorize replacing a live quartic remainder.
        assert!(matches!(
            crate::reduce_galerkin(&nonlinear, &[1.0, 0.0, 0.0, 1.0], 2),
            Err(PhsError::Dimension { .. })
        ));
        // Exercise the storage transport itself without claiming that the
        // public quadratic-only reducer accepts this nonlinear system.
        let reduced_nonlinear = crate::ReducedStorage {
            v: vec![1.0, 0.0, 0.0, 1.0],
            n: 2,
            k: 2,
            inner: nonlinear.storage,
        };
        let nonlinear_part = reduced_nonlinear.quadratic_part().unwrap();
        assert!(!nonlinear_part.exact);
        let point = [0.3, 0.2];
        let quadratic_energy = 0.5 * (point[0] * point[0] + point[1] * point[1]);
        assert!(reduced_nonlinear.hamiltonian(&point) > quadratic_energy);
    }

    #[test]
    fn render_step_keeps_nonlinear_storage_and_refuses_unresolved_clocks() {
        let nonlinear = crate::duffing_oscillator(1.0, 1.0, 7.0, 0.0).unwrap();
        let split = split_quadratic(&nonlinear, budget(2)).unwrap();
        assert!(step_render(&nonlinear, &split, &[0.3, 0.2], &[0.0], 0.04).is_ok());
        assert!(matches!(
            step_render(&nonlinear, &split, &[0.3, 0.2], &[0.0], 0.040_001),
            Err(PhsError::TimeResolution(ResolvedBandError::Outside { .. }))
        ));
        let record = step_render(&nonlinear, &split, &[0.3, 0.2], &[0.0], 0.04).unwrap();
        assert!(record.balance_residual().abs() <= 1e-12 * nonlinear.hamiltonian(&[0.3, 0.2]));
        let linear = mass_spring_damper(1.0, 1.0, 0.0).unwrap();
        let linear_split = split_quadratic(&linear, budget(2)).unwrap();
        let exact = step_render(&linear, &linear_split, &[1.0, 0.0], &[0.0], 2.5).unwrap();
        relative(exact.x[0], det::cos(2.5), 1.0, 1e-12);
    }

    // G1: independent closed-form underdamped oscillator, mass = 1.
    fn analytic(omega: f64, zeta: f64, t: f64, force: f64, initial: [f64; 2]) -> [f64; 2] {
        let a = zeta * omega;
        let w = omega * det::sqrt(1.0 - zeta * zeta);
        let equilibrium = force / (omega * omega);
        let q = initial[0] - equilibrium;
        let b = (initial[1] + a * q) / w;
        let (c, s, e) = (det::cos(w * t), det::sin(w * t), det::exp(-a * t));
        [
            equilibrium + e * (q * c + b * s),
            e * ((-a * q + w * b) * c + (-a * b - w * q) * s),
        ]
    }

    fn relative(actual: f64, expected: f64, scale: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * scale,
            "actual={actual:.17e} expected={expected:.17e} scale={scale:.3e}"
        );
    }

    fn phase(map: [[f64; 2]; 2]) -> f64 {
        let re = 0.5 * (map[0][0] + map[1][1]);
        let determinant = map[0][0] * map[1][1] - map[0][1] * map[1][0];
        let im = det::sqrt(determinant - re * re);
        det::atan2(im, re).abs()
    }

    fn map(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        prop: LinearPropagator,
        dt: f64,
    ) -> [[f64; 2]; 2] {
        let a = step_split(sys, split, prop, &[1.0, 0.0], &[0.0], dt)
            .unwrap()
            .x;
        let b = step_split(sys, split, prop, &[0.0, 1.0], &[0.0], dt)
            .unwrap()
            .x;
        [[a[0], b[0]], [a[1], b[1]]]
    }

    #[test]
    fn exact_zoh_matches_closed_form_damped_oscillator() {
        let (omega, zeta, dt, force) = (2.0, 0.01, 0.0005, 0.3);
        let sys = mass_spring_damper(1.0, omega * omega, 2.0 * zeta * omega).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        let initial = [0.7, -0.2];
        let mut x = initial.to_vec();
        for _ in 0..10_000 {
            x = step_split(&sys, &split, LinearPropagator::ExactZoh, &x, &[force], dt)
                .unwrap()
                .x;
        }
        let expected = analytic(omega, zeta, 10_000.0 * dt, force, initial);
        for i in 0..2 {
            relative(x[i], expected[i], 1.0, 1e-12);
        }
    }

    #[test]
    fn quadratic_storage_steps_without_frequency_warping() {
        let omega = 17.0;
        let sys = mass_spring_damper(1.0, omega * omega, 0.0).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        assert_eq!(split.angular_frequencies_rad_s().len(), 2);
        for z in [0.1, 1.0, 2.5] {
            for prop in [
                LinearPropagator::ExactZoh,
                LinearPropagator::FrequencyMatched,
            ] {
                let observed = phase(map(&sys, &split, prop, z / omega));
                let cents = 1200.0 * det::ln(observed / z) / core::f64::consts::LN_2;
                assert!(
                    cents.abs() <= 0.01,
                    "{prop:?}: {cents} cents at omega*dt={z}"
                );
            }
        }
    }

    #[test]
    fn damped_mode_decays_at_exact_rate() {
        let (omega, zeta) = (19.0, 0.13);
        let sys = mass_spring_damper(1.0, omega * omega, 2.0 * zeta * omega).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        for z in [0.1, 1.0, 2.5] {
            let a = map(&sys, &split, LinearPropagator::ExactZoh, z / omega);
            let radius = det::sqrt(a[0][0] * a[1][1] - a[0][1] * a[1][0]);
            let expected = det::exp(-zeta * z);
            relative(radius, expected, expected, 1e-12);
        }
    }

    #[test]
    fn frequency_matched_refuses_at_pi() {
        let sys = mass_spring_damper(1.0, 1.0, 0.0).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        let pi = core::f64::consts::PI;
        // Band admission is separate from representability/Newton admission:
        // an arbitrarily near-pole step need not have a representable residual.
        assert!(split.quadratic_correction(pi * (1.0 - 1e-9)).is_ok());
        assert!(
            matches!(step_split(&sys,&split,LinearPropagator::FrequencyMatched,&[1.0,0.0],&[0.0],pi),
            Err(PhsError::FrequencyMatchUnadmitted { omega_dt }) if omega_dt == pi)
        );
        let limited = split_quadratic(
            &sys,
            LinearSplitBudget {
                max_modes: 2,
                max_omega_dt: 1.0,
            },
        )
        .unwrap();
        assert!(limited.quadratic_correction(1.0).is_ok());
        assert!(matches!(
            limited.quadratic_correction(1.0 + 1e-12),
            Err(PhsError::FrequencyMatchUnadmitted { .. })
        ));
    }

    #[test]
    fn split_step_keeps_discrete_energy_identity() {
        let sys = mass_spring_damper(1.0, 25.0, 0.4).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        let (x, u, dt) = ([0.3, -0.7], [0.8], 0.1);
        for prop in [
            LinearPropagator::ExactZoh,
            LinearPropagator::FrequencyMatched,
        ] {
            let record = step_split(&sys, &split, prop, &x, &u, dt).unwrap();
            assert!(record.balance_residual().abs() <= 1e-12 * sys.hamiltonian(&x));
            assert!(record.dissipated > 0.0);
            if prop == LinearPropagator::FrequencyMatched {
                let correction = split.quadratic_correction(dt).unwrap();
                let storage = MatchedStorage {
                    original: sys.storage.as_ref(),
                    correction: &correction,
                };
                let dg = discrete_gradient(&storage, &x, &record.x);
                relative(
                    record.dissipated,
                    dt * quadratic(&sys.r, &dg),
                    record.dissipated,
                    1e-12,
                );
            } else {
                // For the force port y = velocity: integral y dt = q1-q0.
                relative(record.supplied, u[0] * (record.x[0] - x[0]), 1.0, 1e-12);
                let expected = analytic(5.0, 0.04, dt, u[0], x);
                for i in 0..2 {
                    relative(record.x[i], expected[i], 1.0, 1e-12);
                }
            }
        }
    }

    #[test]
    fn nonquadratic_remainder_still_uses_gonzalez() {
        let sys = crate::duffing_oscillator(1.0, 9.0, 7.0, 0.0).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        assert!(
            step_split(
                &sys,
                &split,
                LinearPropagator::ExactZoh,
                &[0.3, 0.2],
                &[0.0],
                0.1
            )
            .is_err()
        );
        let correction = split.quadratic_correction(0.1).unwrap();
        let storage = MatchedStorage {
            original: sys.storage.as_ref(),
            correction: &correction,
        };
        let mut x = vec![0.3, 0.2];
        let h0 = storage.hamiltonian(&x);
        for _ in 0..100 {
            let record = step_split(
                &sys,
                &split,
                LinearPropagator::FrequencyMatched,
                &x,
                &[0.0],
                0.1,
            )
            .unwrap();
            assert!(record.balance_residual().abs() <= 1e-12 * h0);
            x = record.x;
        }
        relative(storage.hamiltonian(&x), h0, h0, 1e-12);
        let linear = mass_spring_damper(1.0, 9.0, 0.0).unwrap();
        let linear_split = split_quadratic(&linear, budget(2)).unwrap();
        let linear_step = step_split(
            &linear,
            &linear_split,
            LinearPropagator::FrequencyMatched,
            &[0.3, 0.2],
            &[0.0],
            0.1,
        )
        .unwrap();
        let nonlinear_step = step_split(
            &sys,
            &split,
            LinearPropagator::FrequencyMatched,
            &[0.3, 0.2],
            &[0.0],
            0.1,
        )
        .unwrap();
        assert!(
            (nonlinear_step.x[1] - linear_step.x[1]).abs() > 1e-3,
            "quartic force must be live"
        );
    }

    #[test]
    fn split_budget_cap_and_cap_plus_one() {
        let one = mass_spring_damper(1.0, 1.0, 0.0).unwrap();
        assert!(split_quadratic(&one, budget(2)).is_ok());
        assert!(matches!(
            split_quadratic(&one, budget(1)),
            Err(PhsError::LinearSplitBudget { .. })
        ));
        let three = PortHamiltonian::new(
            3,
            0,
            vec![0.0; 9],
            vec![0.0; 9],
            vec![],
            Box::new(QuadraticStorage::new(identity(3), 3).unwrap()),
        )
        .unwrap();
        assert!(matches!(
            split_quadratic(&three, budget(2)),
            Err(PhsError::LinearSplitBudget { .. })
        ));
        assert!(split_quadratic(&three, budget(3)).is_ok());
    }

    #[test]
    fn exact_zoh_preserves_static_compliance_and_free_mass() {
        let sys = mass_spring_damper(2.0, 25.0, 0.8).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        let record = step_split(
            &sys,
            &split,
            LinearPropagator::ExactZoh,
            &[2.0 / 25.0, 0.0],
            &[2.0],
            0.3,
        )
        .unwrap();
        relative(record.x[0], 2.0 / 25.0, 1.0, 1e-12);
        relative(record.x[1], 0.0, 1.0, 1e-12);
        let free = mass_spring_damper(2.0, 0.0, 0.0).unwrap();
        let prepared = split_quadratic(&free, budget(2)).unwrap();
        let result = step_split(
            &free,
            &prepared,
            LinearPropagator::ExactZoh,
            &[0.3, 0.4],
            &[0.7],
            0.2,
        )
        .unwrap();
        relative(
            result.x[0],
            0.3 + 0.4 / 2.0 * 0.2 + 0.5 * 0.7 / 2.0 * 0.2 * 0.2,
            1.0,
            1e-12,
        );
        relative(result.x[1], 0.4 + 0.7 * 0.2, 1.0, 1e-12);
    }

    #[test]
    fn unsupported_storage_stale_split_and_nonfinite_inputs_refuse() {
        let sys = mass_spring_damper(1.0, 1.0, 0.0).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        let changed = mass_spring_damper(1.0, 2.0, 0.0).unwrap();
        assert!(
            step_split(
                &changed,
                &split,
                LinearPropagator::ExactZoh,
                &[1.0, 0.0],
                &[0.0],
                0.1
            )
            .is_err()
        );
        for dt in [f64::NAN, f64::INFINITY, -1.0] {
            assert!(
                step_split(
                    &sys,
                    &split,
                    LinearPropagator::ExactZoh,
                    &[1.0, 0.0],
                    &[0.0],
                    dt
                )
                .is_err()
            );
        }
        assert!(
            step_split(
                &sys,
                &split,
                LinearPropagator::ExactZoh,
                &[f64::NAN, 0.0],
                &[0.0],
                0.1
            )
            .is_err()
        );
        let unsupported = PortHamiltonian::new(
            1,
            0,
            vec![0.0],
            vec![0.0],
            vec![],
            Box::new(SeparableStorage::new(vec![ScalarLaw {
                h: |x| 0.5 * x * x,
                dh: |x| x,
            }])),
        )
        .unwrap();
        assert!(matches!(
            split_quadratic(&unsupported, budget(1)),
            Err(PhsError::LinearSplitBudget { .. })
        ));
    }

    #[test]
    fn coupled_quadratic_block_retains_off_diagonal_storage_and_loss() {
        let sys = PortHamiltonian::new(
            2,
            0,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.2, 0.1, 0.1, 0.3],
            vec![],
            Box::new(QuadraticStorage::new(vec![5.0, 4.0, 4.0, 5.0], 2).unwrap()),
        )
        .unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        let (dt, initial) = (0.17, [0.3, -0.2]);
        let record =
            step_split(&sys, &split, LinearPropagator::ExactZoh, &initial, &[], dt).unwrap();
        // Independently multiply (J-R)Q: [[2.6,3.7],[-6.7,-5.9]].
        // For A=alpha I+B, trace B=0 and B²=-beta² I,
        // exp(Ah)=exp(alpha h)[cos(beta h)I+sin(beta h)B/beta].
        let a = [2.6, 3.7, -6.7, -5.9];
        let alpha = 0.5 * (a[0] + a[3]);
        let beta = det::sqrt(a[0] * a[3] - a[1] * a[2] - alpha * alpha);
        let envelope = det::exp(alpha * dt);
        let s = det::sin(beta * dt) / beta;
        let c = det::cos(beta * dt);
        for i in 0..2 {
            let expected = envelope
                * (c * initial[i]
                    + s * ((a[i * 2] - if i == 0 { alpha } else { 0.0 }) * initial[0]
                        + (a[i * 2 + 1] - if i == 1 { alpha } else { 0.0 }) * initial[1]));
            relative(record.x[i], expected, 1.0, 1e-12);
        }
        assert!(record.balance_residual().abs() <= 1e-12 * sys.hamiltonian(&initial));
    }

    #[test]
    fn odd_skew_block_keeps_its_zero_frequency_coordinate() {
        let sys = PortHamiltonian::new(
            3,
            0,
            vec![0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, 0.0],
            vec![0.0; 9],
            vec![],
            Box::new(
                QuadraticStorage::new(
                    vec![1.0, 0.0, 0.0, 0.0, 99_999_999.0, 0.0, 0.0, 0.0, 1.0],
                    3,
                )
                .unwrap(),
            ),
        )
        .unwrap();
        let split = split_quadratic(&sys, budget(3)).unwrap();
        assert_eq!(split.angular_frequencies_rad_s().len(), 3);
        assert!(split.angular_frequencies_rad_s()[0] < 1e-10);
        let x = [1.0, 1e-8, 0.0]; // Q*x = [1,1,0] up to declared roundoff.
        let record = step_split(&sys, &split, LinearPropagator::ExactZoh, &x, &[], 1e-5).unwrap();
        assert!(record.balance_residual().abs() <= 1e-12 * sys.hamiltonian(&x));
    }

    #[test]
    fn blocked_product_equals_the_textbook_loop_bit_for_bit() {
        // Sign-mixed entries over forty binades, with exact zeros of both
        // signs; sizes below, at and between multiples of the four-row block.
        let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
        let mut entries = |count: usize| -> Vec<f64> {
            (0..count).map(|i| {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
                match i % 11 {
                    3 => 0.0,
                    7 => -0.0,
                    _ => ((seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 2f64.powi((seed % 40) as i32 - 20),
                }
            }).collect()
        };
        assert!(product(&[], &[], 0).is_empty());
        for n in [1, 2, 3, 4, 5, 7, 8, 9, 13, 64, 99] {
            let (a, b) = (entries(n * n), entries(n * n));
            let mut textbook = vec![0.0; n * n];
            for i in 0..n {
                for k in 0..n {
                    for j in 0..n {
                        textbook[i * n + j] += a[i * n + k] * b[k * n + j];
                    }
                }
            }
            let blocked = product(&a, &b, n);
            assert!(blocked.iter().zip(&textbook).all(|(x, y)| x.to_bits() == y.to_bits()), "n = {n}");
        }
    }

    #[test]
    fn upper_block_exponential_equals_the_assembled_exponential_bit_for_bit() {
        let budget = DenseExponentialBudget {
            max_dimension: 64, max_dense_bytes: usize::MAX, max_multiply_terms: u64::MAX,
        };
        let mut seed = 0x6a09_e667_f3bc_c909_u64;
        let mut entries = |count: usize, size: f64| -> Vec<f64> {
            (0..count).map(|i| {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
                if i % 7 == 3 { 0.0 } else { ((seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * size }
            }).collect()
        };
        let assemble = |blocks: &UpperBlocks, d: usize| -> Vec<f64> {
            let n = 2 * d;
            let mut a = vec![0.0; n * n];
            for i in 0..d {
                for j in 0..d {
                    a[i * n + j] = blocks.first[i * d + j];
                    a[i * n + d + j] = blocks.cross[i * d + j];
                    a[(d + i) * n + d + j] = blocks.last[i * d + j];
                }
            }
            a
        };
        // Norms from inside the series radius (no squaring) to dozens of squarings.
        for d in [1, 2, 3, 5, 8, 13] {
            for size in [0.05, 1.0, 40.0, 3.0e5] {
                let operator = UpperBlocks {
                    first: entries(d * d, size), cross: entries(d * d, size), last: entries(d * d, size),
                };
                let full = exponential_with_budget(&assemble(&operator, d), 2 * d, budget);
                let blocks = upper_block_exponential_with_budget(&operator, d, budget);
                match (full, blocks) {
                    (Ok(full), Ok(blocks)) => {
                        let n = 2 * d;
                        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
                        assert_eq!(bits(&assemble(&blocks, d)), bits(&full), "d = {d}, size = {size}");
                        // The assembled exponential's lower-left block is exactly +0.
                        assert!((0..d * d).all(|k| full[(d + k / d) * n + k % d].to_bits() == 0), "d = {d}");
                    },
                    (Err(full), Err(blocks)) => assert_eq!(full.to_string(), blocks.to_string(), "d = {d}, size = {size}"),
                    (full, blocks) => panic!("d = {d}, size = {size}: {full:?} vs {blocks:?}"),
                }
            }
        }
        // The same refusals: a nonfinite entry, a malformed block, and the 2d admission.
        let mut bad = UpperBlocks { first: vec![1.0; 4], cross: vec![1.0; 4], last: vec![1.0; 4] };
        bad.cross[3] = f64::NAN;
        let mut assembled = assemble(&bad, 2);
        assembled[7] = f64::NAN;
        assert_eq!(upper_block_exponential_with_budget(&bad, 2, budget).unwrap_err().to_string(),
            exponential_with_budget(&assembled, 4, budget).unwrap_err().to_string());
        bad.cross = vec![1.0; 3];
        assert!(upper_block_exponential_with_budget(&bad, 2, budget).is_err());
        let small = DenseExponentialBudget { max_dimension: 3, ..budget };
        let ok = UpperBlocks { first: vec![0.1; 4], cross: vec![0.1; 4], last: vec![0.1; 4] };
        assert_eq!(upper_block_exponential_with_budget(&ok, 2, small).unwrap_err().to_string(),
            exponential_with_budget(&assemble(&ok, 2), 4, small).unwrap_err().to_string());
    }
}
