//! Small-signal response of the original periodically varying tangent.
//!
//! For a base orbit of angular frequency `omega`, solve
//! `sum_l [delta_kl P(i*(nu+k*omega)) - A_(k-l)] X_l = F_k`, where
//! `A(t) = d island(t,x_base(t))/dx`. Both signs of every retained sideband
//! participate. Replacing A(t) by its mean would erase parametric conversion.
//!
//! The physical perturbation is `Re sum_k X_k exp(i*(nu+k*omega)*t)`;
//! coefficients are peak complex amplitudes, not the half-amplitudes in HbOrbit.
//! At coincident positive/negative physical frequencies their contributions
//! must be combined coherently. No sideband is folded, clipped or discarded.
//!
//! This is a finite harmonic-transfer approximation (Sandberg, Moellerstedt
//! and Bernhardsson, IEEE TAC 50(12), 2005, doi:10.1109/TAC.2005.860294).
//! Independent sample grids check the base balance and tangent quadrature;
//! an omitted-band equation check detects unresolved response coupling. These
//! are numerical checks, not infinite-tail or forward-error enclosures. Repeat
//! sideband orders and base-orbit resolution. A finite particular response
//! does not establish Floquet stability or decay of startup transients.
//! The source and its tangent must repeat at the authored base period; finite
//! quadrature does not establish periodicity of an arbitrary time callback.

use crate::{AftSum, HbOrbit, IslandStateJacobian, OrbitProblem, TAU, normalized_difference};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_blocked_with_cancel};
use fs_math::{c64::C64, dd::Dd, det};

/// Independent original state and equation units; numerical conditioning only.
#[derive(Clone, Debug)]
pub struct SidebandScales {
    /// Positive characteristic magnitude of each original state coordinate.
    pub states: Vec<f64>,
    /// Positive characteristic magnitude of each original equation row.
    pub equations: Vec<f64>,
}

/// Explicit finite spectral, sampling, numerical payload and accuracy limits.
#[derive(Clone, Copy, Debug)]
pub struct SidebandControls {
    /// Retain sideband indices -harmonics through +harmonics, at most 64.
    pub harmonics: usize,
    /// Coarse quadrature; a separate 2*samples+1 grid is always evaluated.
    /// Must resolve the authored base order and tangent orders through 3*N+2.
    pub samples_per_period: usize,
    /// Full complex dimension d*(2*N+1), also hard-capped at 1024.
    pub maximum_unknowns: usize,
    /// Conservative peak numerical payload of preparation, one factor and one
    /// solve. Excludes source/caller storage, allocator metadata and stacks.
    pub maximum_payload_bytes: usize,
    /// Bounds numerical preparation plus one factor/solve; source callbacks'
    /// own internal work is excluded. Each later factor/solve has the same cap.
    pub maximum_work: u64,
    /// Sampled RMS of the complete original base equations divided by equation
    /// scales, including defects outside the base's retained Fourier order.
    pub maximum_base_residual: f64,
    /// Absolute Frobenius norm of the two normalized tangent Fourier arrays'
    /// difference. This does not certify absence of aliases on both grids.
    pub maximum_tangent_sampling_defect: f64,
    /// Componentwise backward-error ceiling for the retained complex solve.
    pub maximum_relative_residual: f64,
    /// Norm of normalized omitted equation coefficients for N<|k|<=2*N+2.
    /// Outside input coefficients are exactly zero by the input's definition.
    /// This is an absolute equation-scale allowance, not a forward state bound.
    pub maximum_omitted_residual: f64,
}

/// No partial preparation, factor or response is published on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SidebandError {
    /// Invalid shape, physical domain, or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The original source does not declare an exact complete state tangent.
    TangentUnavailable,
    /// A declared tangent refused or left a nonfinite/unwritten entry.
    TangentRefused {
        /// Zero-based failed quadrature sample.
        sample: usize,
        /// Total samples on the failing grid.
        samples: usize,
    },
    /// The original physical equation left a nonfinite/unwritten entry.
    SourceRefused {
        /// Zero-based failed quadrature sample.
        sample: usize,
        /// Total samples on the failing grid.
        samples: usize,
    },
    /// An explicit count, numerical payload or operation limit was exceeded.
    Budget {
        /// Resource being admitted.
        quantity: &'static str,
        /// Complete requested numerical resource.
        requested: u128,
        /// Authored resource ceiling.
        cap: u128,
    },
    /// The supplied base waveform fails its original sampled equations.
    BaseUnresolved {
        /// Measured original-equation RMS in the declared numerical units.
        residual: f64,
        /// Authored upper limit.
        limit: f64,
    },
    /// The independently sampled periodic tangent has not converged.
    SamplingUnresolved {
        /// Measured normalized Fourier-matrix difference.
        defect: f64,
        /// Authored upper limit.
        limit: f64,
    },
    /// No unique finite retained response exists; no regularization is added.
    Singular {
        /// Zero-based failed elimination pivot in sideband-major ordering.
        index: usize,
    },
    /// The retained complex solve exceeds its backward-error limit.
    Residual {
        /// Measured componentwise backward error.
        measured: f64,
        /// Authored upper limit.
        limit: f64,
    },
    /// The solved waveform leaves excessive forcing in checked omitted bands.
    OmittedSidebands {
        /// Measured norm in the declared original-equation numerical units.
        measured: f64,
        /// Authored upper limit.
        limit: f64,
    },
    /// Cancellation observed before publication.
    Cancelled,
}
impl core::fmt::Display for SidebandError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-ORBIT-SIDEBAND: {self:?}")
    }
}
impl std::error::Error for SidebandError {}
impl From<LuComplexError> for SidebandError {
    fn from(value: LuComplexError) -> Self {
        match value {
            LuComplexError::Cancelled => Self::Cancelled,
            LuComplexError::Singular { index } => Self::Singular { index },
        }
    }
}
type Result<T> = core::result::Result<T, SidebandError>;
fn invalid(what: &'static str) -> SidebandError {
    SidebandError::Invalid(what)
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SidebandError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn bound(quantity: &'static str, requested: u128, cap: u128) -> Result<()> {
    if requested > cap {
        Err(SidebandError::Budget {
            quantity,
            requested,
            cap,
        })
    } else {
        Ok(())
    }
}

/// Original source and its cached periodic tangent, shared by probe frequencies.
/// The source must remain deterministic and unchanged for this borrow's lifetime.
pub struct PreparedPeriodicLinearization<'a, P> {
    problem: &'a P,
    dimension: usize,
    omega: f64,
    controls: SidebandControls,
    scales: SidebandScales,
    tangent: Vec<C64>,
    /// Maximum independently recomputed original base balance on the two grids.
    pub base_residual: f64,
    /// Measured normalized tangent difference between the two sample grids.
    pub tangent_sampling_defect: f64,
}

impl<'a, P: OrbitProblem> PreparedPeriodicLinearization<'a, P> {
    /// Prepare the full physical tangent, without inferring one from finite
    /// differences or trusting the supplied orbit's previous convergence fields.
    /// Cancellation is checked per sample, Fourier order and source/port call.
    pub fn new(
        problem: &'a P,
        orbit: &HbOrbit,
        controls: SidebandControls,
        scales: &SidebandScales,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let d = problem.dim();
        admit(d, orbit, controls, scales)?;
        let callback = problem
            .island_state_jacobian()
            .ok_or(SidebandError::TangentUnavailable)?;
        let coarse = sample(
            problem,
            callback,
            orbit,
            controls.harmonics,
            controls.samples_per_period,
            scales,
            &mut cancelled,
        )?;
        let fine = sample(
            problem,
            callback,
            orbit,
            controls.harmonics,
            2 * controls.samples_per_period + 1,
            scales,
            &mut cancelled,
        )?;
        let base_residual = coarse.1.max(fine.1);
        if base_residual > controls.maximum_base_residual {
            return Err(SidebandError::BaseUnresolved {
                residual: base_residual,
                limit: controls.maximum_base_residual,
            });
        }
        let mut tangent_sampling_defect = 0.0;
        for (&a, &b) in coarse.0.iter().zip(&fine.0) {
            tangent_sampling_defect = det::hypot(tangent_sampling_defect, (a - b).abs());
        }
        if !tangent_sampling_defect.is_finite()
            || tangent_sampling_defect > controls.maximum_tangent_sampling_defect
        {
            return Err(SidebandError::SamplingUnresolved {
                defect: tangent_sampling_defect,
                limit: controls.maximum_tangent_sampling_defect,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            problem,
            dimension: d,
            omega: orbit.omega,
            controls,
            scales: scales.clone(),
            tangent: fine.0,
            base_residual,
            tangent_sampling_defect,
        })
    }

    /// Number of complete original state coordinates.
    pub fn dimension(&self) -> usize {
        self.dimension
    }
    /// Largest retained signed sideband index.
    pub fn harmonics(&self) -> usize {
        self.controls.harmonics
    }
    /// Authored base-orbit angular frequency in radians per unit time.
    pub fn orbital_angular_frequency(&self) -> f64 {
        self.omega
    }

    /// Factor the complete coupled sideband system once for many excitations.
    /// Negative shifted frequencies are passed to the original linear port;
    /// no time realization, conjugacy or out-of-domain extrapolation is inferred.
    /// No damping floor, phase pin or diagonal decoupling resolves singularities.
    pub fn factor(
        &self,
        carrier_angular_frequency: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedSidebandResponse<'_>> {
        poll(&mut cancelled)?;
        if !carrier_angular_frequency.is_finite() {
            return Err(invalid("finite signed carrier angular frequency required"));
        }
        let d = self.dimension;
        let h = self.controls.harmonics;
        let n = d * (2 * h + 1);
        let mut operator = vec![C64::ZERO; n * n];
        for band in 0..2 * h + 1 {
            poll(&mut cancelled)?;
            let k = band as isize - h as isize;
            let frequency = (k as f64).mul_add(self.omega, carrier_angular_frequency);
            if !frequency.is_finite() {
                return Err(invalid("shifted frequency is not representable"));
            }
            let port = self.problem.port(C64::new(0.0, frequency));
            poll(&mut cancelled)?;
            check_port(&port, d)?;
            for row in 0..d {
                poll(&mut cancelled)?;
                for other in 0..2 * h + 1 {
                    let order = band as isize - other as isize;
                    for col in 0..d {
                        let mut value = -coefficient(&self.tangent, d, order, row, col);
                        if band == other {
                            value =
                                value + scale_entry(port[row * d + col], row, col, &self.scales)?;
                        }
                        if !finite(value) {
                            return Err(invalid("sideband operator is not representable"));
                        }
                        operator[(band * d + row) * n + other * d + col] = value;
                    }
                }
            }
        }
        // One explicitly selected worker keeps this numerical payload bound
        // independent of an outer thread-local LU scheduling declaration.
        let factor = lu_complex_blocked_with_cancel(&operator, n, 1, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(PreparedSidebandResponse {
            dimension: d,
            omega: self.omega,
            carrier: carrier_angular_frequency,
            controls: self.controls,
            scales: &self.scales,
            tangent: &self.tangent,
            operator,
            factor,
        })
    }
}

/// One reusable complex factor in the original source's periodic state chart.
pub struct PreparedSidebandResponse<'a> {
    dimension: usize,
    omega: f64,
    carrier: f64,
    controls: SidebandControls,
    scales: &'a SidebandScales,
    tangent: &'a [C64],
    operator: Vec<C64>,
    factor: LuComplex,
}

/// Peak complex perturbations in sideband-major, original-state order.
#[derive(Clone, Debug)]
pub struct SidebandResponse {
    /// Original physical state units, indexed by `(k+N)*dimension+state`.
    pub coefficients: Vec<C64>,
    /// Maximum componentwise backward error of the retained complex equations.
    pub relative_residual: f64,
    /// Norm of original-equation defects, divided by declared equation scales,
    /// in the explicitly checked omitted sidebands N<|k|<=2*N+2.
    pub omitted_residual: f64,
}

impl PreparedSidebandResponse<'_> {
    /// Number of complete original state coordinates.
    pub fn dimension(&self) -> usize {
        self.dimension
    }
    /// Largest retained signed sideband index.
    pub fn harmonics(&self) -> usize {
        self.controls.harmonics
    }
    /// Authored base-orbit angular frequency in radians per unit time.
    pub fn orbital_angular_frequency(&self) -> f64 {
        self.omega
    }
    /// Signed probe carrier angular frequency in radians per unit time.
    pub fn carrier_angular_frequency(&self) -> f64 {
        self.carrier
    }

    /// Solve complete original-equation forcing phasors. Index `(k+N)*d+j`
    /// addresses state/equation j at carrier+k*omega. Unlisted sidebands are
    /// zero by definition, not an implicit truncation of a supplied input.
    /// The cached source tangent and factor serve every RHS without callbacks.
    pub fn solve(
        &self,
        forcing: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SidebandResponse> {
        poll(&mut cancelled)?;
        let d = self.dimension;
        let h = self.controls.harmonics;
        let n = d * (2 * h + 1);
        if forcing.len() != n || forcing.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "complete finite sideband-major equation forcing required",
            ));
        }
        let rhs: Vec<_> = forcing
            .iter()
            .enumerate()
            .map(|(i, z)| {
                C64::new(
                    z.re / self.scales.equations[i % d],
                    z.im / self.scales.equations[i % d],
                )
            })
            .collect();
        if rhs.iter().zip(forcing).any(|(&a, &b)| {
            !finite(a) || (b.re != 0.0 && a.re == 0.0) || (b.im != 0.0 && a.im == 0.0)
        }) {
            return Err(invalid("normalized forcing is not representable"));
        }
        let mut solution = rhs.clone();
        self.factor
            .solve_with_cancel(&mut solution, &mut cancelled)?;
        if solution.iter().any(|&z| !finite(z)) {
            return Err(invalid("sideband solution is not representable"));
        }
        let relative_residual = residual(&self.operator, &solution, &rhs, &mut cancelled)?;
        if relative_residual > self.controls.maximum_relative_residual {
            return Err(SidebandError::Residual {
                measured: relative_residual,
                limit: self.controls.maximum_relative_residual,
            });
        }
        let mut omitted_residual = 0.0;
        for magnitude in h + 1..=2 * h + 2 {
            for k in [-(magnitude as isize), magnitude as isize] {
                for row in 0..d {
                    poll(&mut cancelled)?;
                    let mut sum = ComplexSum::new();
                    for band in 0..2 * h + 1 {
                        let order = k - (band as isize - h as isize);
                        for col in 0..d {
                            sum.product(
                                coefficient(self.tangent, d, order, row, col),
                                solution[band * d + col],
                            );
                        }
                    }
                    omitted_residual = det::hypot(omitted_residual, sum.value().abs());
                }
            }
        }
        if !omitted_residual.is_finite()
            || omitted_residual > self.controls.maximum_omitted_residual
        {
            return Err(SidebandError::OmittedSidebands {
                measured: omitted_residual,
                limit: self.controls.maximum_omitted_residual,
            });
        }
        for (i, value) in solution.iter_mut().enumerate() {
            let old = *value;
            *value = value.scale(self.scales.states[i % d]);
            if !finite(*value)
                || (old.re != 0.0 && value.re == 0.0)
                || (old.im != 0.0 && value.im == 0.0)
            {
                return Err(invalid("physical sideband response is not representable"));
            }
        }
        poll(&mut cancelled)?;
        Ok(SidebandResponse {
            coefficients: solution,
            relative_residual,
            omitted_residual,
        })
    }
}

fn admit(d: usize, orbit: &HbOrbit, c: SidebandControls, s: &SidebandScales) -> Result<()> {
    let width = orbit.coeffs.first().map_or(0, Vec::len);
    if d == 0
        || width == 0
        || width > 513
        || c.harmonics > 64
        || !positive(orbit.omega)
        || !positive(TAU / orbit.omega)
        || orbit.coeffs.len() != d
        || orbit
            .coeffs
            .iter()
            .any(|r| r.len() != width || r[0].im != 0.0 || r.iter().any(|&z| !finite(z)))
        || s.states.len() != d
        || s.equations.len() != d
        || s.states.iter().chain(&s.equations).any(|&v| !positive(v))
        || !positive(c.maximum_base_residual)
        || !positive(c.maximum_tangent_sampling_defect)
        || !positive(c.maximum_relative_residual)
        || c.maximum_relative_residual >= 1.0
        || !positive(c.maximum_omitted_residual)
        || !(8..=65_536).contains(&c.samples_per_period)
        || c.samples_per_period < 2 * (3 * c.harmonics + 2).max(width - 1) + 1
    {
        return Err(invalid(
            "finite full orbit, physical scales, resolved quadrature and positive limits required",
        ));
    }
    let dimension = d as u128;
    let n = dimension * (2 * c.harmonics + 1) as u128;
    bound("complex unknowns", n, c.maximum_unknowns.min(1024) as u128)?;
    let entries = dimension * dimension * (3 * c.harmonics + 3) as u128;
    let base = dimension * width as u128;
    // Includes both Fourier grids, their compensated sums, original callbacks,
    // retained operator/factor, LU tile descriptors, solve and residual scratch.
    let payload = 80 * entries
        + 96 * base
        + 128 * dimension * dimension
        + 128 * dimension
        + 96 * n * n
        + 512 * n
        + 32 * width.max(3 * c.harmonics + 3) as u128;
    bound(
        "numerical payload bytes",
        payload,
        c.maximum_payload_bytes as u128,
    )?;
    let samples = (3 * c.samples_per_period + 1) as u128;
    let work = samples * (16 * entries + 32 * base + 32 * dimension * dimension)
        + 16 * n * n * n
        + 128 * n * n;
    bound("numerical work", work, c.maximum_work as u128)
}

fn check_port(port: &[C64], d: usize) -> Result<()> {
    if port.len() != d * d || port.iter().any(|&z| !finite(z)) {
        Err(invalid(
            "original port must provide a complete finite matrix at every queried frequency",
        ))
    } else {
        Ok(())
    }
}

fn scale_entry(z: C64, row: usize, col: usize, s: &SidebandScales) -> Result<C64> {
    let value = C64::new(
        (z.re / s.equations[row]) * s.states[col],
        (z.im / s.equations[row]) * s.states[col],
    );
    if !finite(value) || (z.re != 0.0 && value.re == 0.0) || (z.im != 0.0 && value.im == 0.0) {
        Err(invalid("normalized original operator is not representable"))
    } else {
        Ok(value)
    }
}

fn coefficient(a: &[C64], d: usize, order: isize, row: usize, col: usize) -> C64 {
    let value = a[order.unsigned_abs() * d * d + row * d + col];
    if order < 0 { value.conj() } else { value }
}

/// Shared double-double owner for independent original-row contractions.
struct ComplexSum {
    re: Dd,
    im: Dd,
}
impl ComplexSum {
    fn new() -> Self {
        Self {
            re: Dd::ZERO,
            im: Dd::ZERO,
        }
    }
    fn product(&mut self, a: C64, b: C64) {
        let [ar, ai, br, bi] = [a.re, a.im, b.re, b.im].map(Dd::from_f64);
        self.re = self.re + ar * br - ai * bi;
        self.im = self.im + ar * bi + ai * br;
    }
    fn value(&self) -> C64 {
        C64::new(self.re.to_f64(), self.im.to_f64())
    }
}

fn residual(
    a: &[C64],
    x: &[C64],
    rhs: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let n = x.len();
    let mut maximum: f64 = 0.0;
    for row in 0..n {
        poll(cancelled)?;
        let mut sum = ComplexSum::new();
        let mut scale = AftSum::new(rhs[row].abs());
        for col in 0..n {
            sum.product(a[row * n + col], x[col]);
            scale.add(a[row * n + col].abs() * x[col].abs());
        }
        sum.re = sum.re - Dd::from_f64(rhs[row].re);
        sum.im = sum.im - Dd::from_f64(rhs[row].im);
        let defect = sum.value().abs();
        let scale = scale.total();
        if !defect.is_finite() || !scale.is_finite() {
            return Err(invalid("sideband residual is not representable"));
        }
        maximum = maximum.max(if scale == 0.0 {
            if defect == 0.0 { 0.0 } else { f64::INFINITY }
        } else {
            defect / scale
        });
    }
    Ok(maximum)
}

#[allow(clippy::too_many_arguments)]
fn sample<P: OrbitProblem>(
    problem: &P,
    callback: IslandStateJacobian<P>,
    orbit: &HbOrbit,
    h: usize,
    samples: usize,
    scales: &SidebandScales,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, f64)> {
    let d = problem.dim();
    let square = d * d;
    let base_order = orbit.coeffs[0].len() - 1;
    let tangent_order = 3 * h + 2;
    let max_order = base_order.max(tangent_order);
    let mut tangent_re: Vec<_> = (0..(tangent_order + 1) * square)
        .map(|_| AftSum::new(0.0))
        .collect();
    let mut tangent_im: Vec<_> = (0..tangent_re.len()).map(|_| AftSum::new(0.0)).collect();
    // Evaluate the original linear port on EVERY supplied base coefficient.
    // The time-domain comparison below measures the full sampled equation,
    // not only a converged projection which can hide generated harmonics.
    let mut linear = vec![C64::ZERO; (base_order + 1) * d];
    for order in 0..=base_order {
        poll(cancelled)?;
        let frequency = order as f64 * orbit.omega;
        if !frequency.is_finite() {
            return Err(invalid("base harmonic frequency is not representable"));
        }
        let port = problem.port(C64::new(0.0, frequency));
        poll(cancelled)?;
        check_port(&port, d)?;
        for row in 0..d {
            let mut sum = ComplexSum::new();
            for col in 0..d {
                let z = orbit.coeffs[col][order];
                let value = C64::new(z.re / scales.states[col], z.im / scales.states[col]);
                if !finite(value)
                    || (z.re != 0.0 && value.re == 0.0)
                    || (z.im != 0.0 && value.im == 0.0)
                {
                    return Err(invalid("normalized base coefficient is not representable"));
                }
                sum.product(scale_entry(port[row * d + col], row, col, scales)?, value);
            }
            linear[order * d + row] = sum.value();
        }
    }
    if linear.iter().any(|&z| !finite(z)) {
        return Err(invalid("base linear response is not representable"));
    }
    let mut reference = vec![0.0; square];
    let mut x = vec![0.0; d];
    let mut g = vec![0.0; d];
    let mut jacobian = vec![0.0; square];
    let mut cosine = vec![0.0; max_order + 1];
    let mut sine = vec![0.0; max_order + 1];
    let period = TAU / orbit.omega;
    let sample_weight = 1.0 / det::sqrt(samples as f64);
    let mut base_residual = 0.0;
    for value in &linear[..d] {
        base_residual = det::hypot(base_residual, value.im);
    }
    for sample in 0..samples {
        poll(cancelled)?;
        let phase = TAU * (sample as f64 / samples as f64);
        let time = period * (sample as f64 / samples as f64);
        for order in 0..=max_order {
            (sine[order], cosine[order]) = det::sin_cos(order as f64 * phase);
        }
        for (row, state) in x.iter_mut().enumerate() {
            let mut sum = AftSum::new(orbit.coeffs[row][0].re);
            for order in 1..=base_order {
                let z = orbit.coeffs[row][order];
                sum.add(2.0 * (z.re * cosine[order] - z.im * sine[order]));
            }
            *state = sum.total();
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Err(invalid("base waveform is not representable"));
        }
        g.fill(f64::NAN);
        problem.island(time, &x, &mut g);
        poll(cancelled)?;
        if g.iter().any(|v| !v.is_finite()) {
            return Err(SidebandError::SourceRefused { sample, samples });
        }
        jacobian.fill(f64::NAN);
        if !callback(problem, time, &x, &mut jacobian) || jacobian.iter().any(|v| !v.is_finite()) {
            return Err(SidebandError::TangentRefused { sample, samples });
        }
        poll(cancelled)?;
        for row in 0..d {
            g[row] /= scales.equations[row];
            if !g[row].is_finite() {
                return Err(invalid("normalized source is not representable"));
            }
            let mut balance = Dd::from_f64(linear[row].re) - Dd::from_f64(g[row]);
            for order in 1..=base_order {
                let z = linear[order * d + row];
                balance = balance
                    + Dd::from_f64(2.0)
                        * (Dd::from_f64(z.re) * Dd::from_f64(cosine[order])
                            - Dd::from_f64(z.im) * Dd::from_f64(sine[order]));
            }
            base_residual = det::hypot(base_residual, sample_weight * balance.to_f64());
            for col in 0..d {
                jacobian[row * d + col] =
                    scale_entry(C64::from_re(jacobian[row * d + col]), row, col, scales)?.re;
            }
        }
        if sample == 0 {
            reference.copy_from_slice(&jacobian);
        }
        for (a, &b) in jacobian.iter_mut().zip(&reference) {
            *a = normalized_difference(*a, b, samples);
        }
        for order in 0..=tangent_order {
            poll(cancelled)?;
            for entry in 0..square {
                tangent_re[order * square + entry].add(cosine[order] * jacobian[entry]);
                tangent_im[order * square + entry].add(-sine[order] * jacobian[entry]);
            }
        }
    }
    for (sum, &value) in tangent_re.iter_mut().zip(&reference) {
        sum.add(value);
    }
    let tangent: Vec<_> = tangent_re
        .iter()
        .zip(&tangent_im)
        .map(|(a, b)| C64::new(a.total(), b.total()))
        .collect();
    if tangent.iter().any(|&z| !finite(z)) {
        return Err(invalid("periodic tangent is not representable"));
    }
    if !base_residual.is_finite() {
        return Err(invalid("base balance is not representable"));
    }
    Ok((tangent, base_residual))
}
