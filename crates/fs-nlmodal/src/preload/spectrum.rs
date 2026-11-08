//! Conservative harmonics about the held preload, from the original energy.

use super::PreloadedModalStorage;
use crate::duffing_kinematic_spectrum_exact;
use fs_math::{dd::Dd, det, qd::Qd};

/// Bounded cold quadrature; no warm simulation state is changed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PreloadedSpectrumBudget {
    /// Maximum even theta intervals; a power of two in 64..=65536.
    pub max_intervals: usize,
    /// Requested relative frequency and normalized kinematic refinement change.
    /// This is a numerical convergence tolerance, not an interval enclosure.
    pub tolerance: f64,
}

impl Default for PreloadedSpectrumBudget {
    fn default() -> Self {
        Self {
            max_intervals: 16384,
            tolerance: 1e-8,
        }
    }
}

/// Measured finite omitted band, never a bound on the infinite Fourier tail.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PreloadedSpectrumGuard {
    /// Last evaluated order; the guard comprises H+1 through this order (2H).
    pub last_harmonic: usize,
    /// Sum of absolute displacement coefficients in this finite band.
    pub displacement_l1: f64,
    /// Corresponding physical velocity sum, in displacement units per second.
    pub velocity_l1: f64,
    /// Corresponding physical acceleration sum, in displacement units per second².
    pub acceleration_l1: f64,
}

/// Conservative Fourier seed in the source's incremental displacement units.
#[derive(Debug, Clone, PartialEq)]
pub struct PreloadedSpectrum {
    /// Fundamental physical angular frequency, in radians per second.
    pub angular_frequency: f64,
    /// Authored positive turning displacement at phase zero.
    pub positive_turning_displacement: f64,
    /// Opposite turning displacement, generally not minus the positive one.
    pub negative_turning_displacement: f64,
    /// Mean incremental displacement; the static reference is not added again.
    pub displacement_dc: f64,
    /// Signed peak cosine coefficients for EVERY order h=1..=H.
    /// Positive-frequency complex displacement coefficients are half these.
    pub cosine_amplitudes: Vec<f64>,
    /// Signed peak sine coefficients -h*omega*a_h in physical velocity units.
    pub velocity_sine_amplitudes: Vec<f64>,
    /// Signed peak cosine coefficients -(h*omega)²*a_h in acceleration units.
    pub acceleration_cosine_amplitudes: Vec<f64>,
    /// Directly computed finite guard band. Small guard content alone does not
    /// certify all higher harmonics; increase H to check spectral resolution.
    pub guard: PreloadedSpectrumGuard,
    /// Final quadrature interval count; zero identifies the exact symmetric path.
    pub intervals: usize,
    /// Maximum relative frequency change of the last two grid refinements.
    pub relative_frequency_change: f64,
    /// Maximum of the last two changes in normalized DC plus
    /// sum(h²*abs(delta a_h/A)), including all guard coefficients.
    /// Controls refinement of phase derivatives, not an infinite-tail estimate.
    pub normalized_kinematic_change: f64,
}

/// Explicit domain, coupling, numerical-range or work-budget refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PreloadedSpectrumError {
    /// Invalid axis, amplitude, harmonic count or quadrature controls.
    InvalidInput,
    /// An original stress channel drives transverse coordinates on this axis.
    NonInvariantAxis,
    /// The scalar potential is not globally convex and confining.
    NonConvex,
    /// Required source coefficients or physical outputs are not representable.
    NumericalRange,
    /// Two consecutive grid-refinement agreements were not obtained.
    QuadratureUnresolved {
        /// Last evaluated interval count.
        intervals: usize,
        /// Last relative frequency change; infinity before the second grid.
        relative_frequency_change: f64,
        /// Last normalized coefficient change; infinity before the second grid.
        normalized_kinematic_change: f64,
    },
}

impl core::fmt::Display for PreloadedSpectrumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "preloaded conservative spectrum: {self:?}")
    }
}
impl std::error::Error for PreloadedSpectrumError {}

type Error = PreloadedSpectrumError;

impl PreloadedModalStorage {
    /// Harmonics of an admitted, exactly invariant source coordinate under its
    /// held dead load. Motion starts at positive displacement A with zero
    /// velocity; time reversal gives a cosine displacement series, including
    /// its physical DC shift and even orders. No forcing or damping is added.
    ///
    /// Each active channel must preserve the axis individually, including
    /// preload products. No cross-channel cancellation or small-coupling cutoff
    /// creates a scalar model. The scalar force is k*q + alpha*q² + beta*q³,
    /// using the original cached tangent and channel contractions. Its potential
    /// must be globally convex: k,beta >= 0 and alpha² <= 3*k*beta, with a
    /// nonzero confining coefficient. Nonconvex wells and separatrices refuse.
    ///
    /// A turning-point substitution removes both energy-integral endpoint
    /// singularities. Composite Simpson quadrature computes the period, physical
    /// phase and all coefficients through order 2H; two consecutive refinements
    /// must meet the requested tolerance. H is in 1..=128. Work is bounded by
    /// O(max_intervals*H), scratch by O(max_intervals+H); all work is cold.
    /// The symmetric law reuses the exact Jacobi kinematic owner instead.
    ///
    /// Finite guard content is reported separately from quadrature convergence.
    /// There is no claimed infinite-tail bound or coupled stability certificate,
    /// and retained coefficients are never renormalized to hide omitted motion.
    ///
    /// # Errors
    /// Invalid input, non-invariant axis, nonconvex scalar law, unrepresentable
    /// coefficients/output, or insufficient quadrature budget.
    pub fn invariant_axis_spectrum(
        &self,
        mode: usize,
        positive_turning_displacement: f64,
        harmonic_order: usize,
        budget: PreloadedSpectrumBudget,
    ) -> Result<PreloadedSpectrum, Error> {
        let amplitude = positive_turning_displacement;
        if mode >= self.mode_count()
            || !amplitude.is_finite()
            || amplitude < 0.0
            || !(1..=128).contains(&harmonic_order)
            || !budget.max_intervals.is_power_of_two()
            || !(64..=65536).contains(&budget.max_intervals)
            || !budget.tolerance.is_finite()
            || budget.tolerance <= 0.0
            || budget.tolerance >= 1.0
        {
            return Err(Error::InvalidInput);
        }
        let (k, alpha, beta) = self.axis_coefficients(mode)?;
        let left = Scaled::product(&[alpha, alpha]);
        let right = Scaled::product(&[3.0, k, beta]);
        let exponent = left.exponent.max(right.exponent);
        if k < 0.0
            || beta < 0.0
            || (k == 0.0 && beta == 0.0)
            || (alpha != 0.0 && (k == 0.0 || beta == 0.0 || left.at(exponent) > right.at(exponent)))
        {
            return Err(Error::NonConvex);
        }
        if alpha == 0.0 || amplitude == 0.0 {
            return symmetric(k, beta, amplitude, harmonic_order);
        }
        let scaled = [
            Scaled::product(&[k]),
            Scaled::product(&[alpha, amplitude]),
            Scaled::product(&[beta, amplitude, amplitude]),
        ];
        let rate_exponent = scaled
            .iter()
            .filter(|s| !s.mantissa.is_zero())
            .map(|s| s.exponent)
            .max()
            .ok_or(Error::NonConvex)?
            .div_euclid(2);
        let coefficients: Vec<f64> = scaled
            .iter()
            .map(|s| s.at(2 * rate_exponent).to_f64())
            .collect();
        if coefficients
            .iter()
            .zip(&scaled)
            .any(|(x, s)| !x.is_finite() || (*x == 0.0 && !s.mantissa.is_zero()))
        {
            return Err(Error::NumericalRange);
        }
        let (k, alpha, beta) = (coefficients[0], coefficients[1], coefficients[2]);
        let energy = 0.5 * k + alpha / 3.0 + 0.25 * beta;
        let potential = |y: f64| y * y * (0.5 * k + y * (alpha / 3.0 + 0.25 * beta * y));
        let mut upper = 1.0;
        for _ in 0..32 {
            if potential(-upper) >= energy {
                break;
            }
            upper *= 2.0;
        }
        if !energy.is_finite() || energy <= 0.0 || potential(-upper) < energy {
            return Err(Error::NumericalRange);
        }
        let mut lower = 0.0;
        for _ in 0..64 {
            let middle = 0.5 * (lower + upper);
            if middle == lower || middle == upper {
                break;
            }
            if potential(-middle) < energy {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        let opposite = -0.5 * (lower + upper);
        // 2(E-V(y))=(1-y)(y-opposite)*(r0+r1*y+r2*y²).
        let r2 = 0.5 * beta;
        let r1 = 2.0 * alpha / 3.0 + r2 * (1.0 + opposite);
        let polynomial = [-2.0 * energy / opposite, r1, r2];
        let mut intervals = (8 * (2 * harmonic_order + 1)).max(64).next_power_of_two();
        if intervals > budget.max_intervals {
            return Err(Error::InvalidInput);
        }
        let mut previous: Option<Grid> = None;
        let mut previous_changes = (f64::INFINITY, f64::INFINITY);
        let mut changes = previous_changes;
        loop {
            let current = integrate(opposite, polynomial, 2 * harmonic_order, intervals)?;
            if let Some(old) = &previous {
                changes = (
                    (current.omega - old.omega).abs() / current.omega,
                    current
                        .coefficients
                        .iter()
                        .zip(&old.coefficients)
                        .enumerate()
                        .map(|(h, (a, b))| (h * h).max(1) as f64 * (a - b).abs())
                        .sum(),
                );
                if changes.0.max(previous_changes.0) <= budget.tolerance
                    && changes.1.max(previous_changes.1) <= budget.tolerance
                {
                    let omega = Scaled {
                        mantissa: Qd::from_f64(current.omega),
                        exponent: rate_exponent,
                    }
                    .value()?;
                    if omega <= 0.0 {
                        return Err(Error::NumericalRange);
                    }
                    let mut result = physical_spectrum(
                        amplitude,
                        opposite,
                        omega,
                        &current.coefficients,
                        harmonic_order,
                    )?;
                    result.intervals = intervals;
                    result.relative_frequency_change = changes.0.max(previous_changes.0);
                    result.normalized_kinematic_change = changes.1.max(previous_changes.1);
                    return Ok(result);
                }
                previous_changes = changes;
            }
            if intervals == budget.max_intervals {
                return Err(Error::QuadratureUnresolved {
                    intervals,
                    relative_frequency_change: changes.0,
                    normalized_kinematic_change: changes.1,
                });
            }
            previous = Some(current);
            intervals *= 2;
        }
    }

    fn axis_coefficients(&self, mode: usize) -> Result<(f64, f64, f64), Error> {
        let n = self.mode_count();
        let mut alpha = Qd::ZERO;
        let mut beta = Qd::ZERO;
        for (source, reference) in self.source.channels.iter().zip(&self.channels) {
            if source.coefficient == 0.0 {
                continue;
            }
            let d = source.coupling[mode * n + mode];
            let b = reference.product[mode];
            for i in 0..n {
                if i == mode {
                    continue;
                }
                let off_diagonal = source.coupling[i * n + mode];
                if (d != 0.0 || b != 0.0) && (off_diagonal != 0.0 || reference.product[i] != 0.0)
                    || (d == 0.0 && b == 0.0 && reference.form != 0.0 && off_diagonal != 0.0)
                {
                    return Err(Error::NonInvariantAxis);
                }
            }
            alpha = alpha + Scaled::product(&[3.0, source.coefficient, b, d]).coefficient()?;
            beta = beta + Scaled::product(&[source.coefficient, d, d]).coefficient()?;
        }
        let k = self.local_stiffness()[mode * n + mode];
        let a = alpha.to_f64();
        let b = beta.to_f64();
        if !a.is_finite()
            || !b.is_finite()
            || (a == 0.0 && !alpha.is_zero())
            || (b == 0.0 && !beta.is_zero())
        {
            return Err(Error::NumericalRange);
        }
        Ok((k, a, b))
    }
}

struct Grid {
    omega: f64,
    coefficients: Vec<f64>,
}

fn integrate(opposite: f64, r: [f64; 3], highest: usize, intervals: usize) -> Result<Grid, Error> {
    let pi = core::f64::consts::PI;
    let step = pi / intervals as f64;
    let center = 0.5 * (1.0 + opposite);
    let radius = 0.5 * (1.0 - opposite);
    let mut nodes = Vec::with_capacity(intervals + 1);
    for j in 0..=intervals {
        let y = center + radius * det::cos(j as f64 * step);
        let speed_squared = r[0] + y * (r[1] + r[2] * y);
        if !speed_squared.is_finite() || speed_squared <= 0.0 {
            return Err(Error::NumericalRange);
        }
        nodes.push([y, 1.0 / det::sqrt(speed_squared), 0.0]);
    }
    let mut time = Dd::ZERO;
    for j in (0..intervals).step_by(2) {
        let (a, b, c) = (nodes[j][1], nodes[j + 1][1], nodes[j + 2][1]);
        nodes[j + 1][2] = time.to_f64() + step / 12.0 * (5.0 * a + 8.0 * b - c);
        time = time + Dd::from_f64(step / 3.0 * (a + 4.0 * b + c));
        nodes[j + 2][2] = time.to_f64();
        if nodes[j + 1][2] <= nodes[j][2] || nodes[j + 1][2] >= nodes[j + 2][2] {
            return Err(Error::NumericalRange);
        }
    }
    let half_period = time.to_f64();
    let omega = pi / half_period;
    if !omega.is_finite() || omega <= 0.0 {
        return Err(Error::NumericalRange);
    }
    let mut sums = vec![Dd::ZERO; highest + 1];
    for (j, node) in nodes.iter().enumerate() {
        let weight = if j == 0 || j == intervals {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let factor = weight * node[0] * node[1];
        let phase = pi * (node[2] / half_period);
        sums[0] = sums[0] + Dd::from_f64(factor);
        for (h, sum) in sums.iter_mut().enumerate().skip(1) {
            *sum = *sum + Dd::from_f64(factor * det::cos(h as f64 * phase));
        }
    }
    let coefficients = sums
        .iter()
        .enumerate()
        .map(|(h, sum)| {
            sum.to_f64() * (step / (3.0 * half_period)) * if h == 0 { 1.0 } else { 2.0 }
        })
        .collect();
    Ok(Grid {
        omega,
        coefficients,
    })
}

fn physical_spectrum(
    amplitude: f64,
    opposite: f64,
    omega: f64,
    coefficients: &[f64],
    order: usize,
) -> Result<PreloadedSpectrum, Error> {
    let mut q = Vec::with_capacity(order);
    let mut v = Vec::with_capacity(order);
    let mut a = Vec::with_capacity(order);
    let mut guard = [Dd::ZERO; 3];
    for (h, &coefficient) in coefficients.iter().enumerate().skip(1) {
        let displacement = Scaled::product(&[amplitude, coefficient]).value()?;
        let velocity = -Scaled::product(&[amplitude, coefficient, omega, h as f64]).value()?;
        let acceleration =
            -Scaled::product(&[amplitude, coefficient, omega, omega, h as f64, h as f64])
                .value()?;
        if h <= order {
            q.push(displacement);
            v.push(velocity);
            a.push(acceleration);
        } else {
            for (sum, value) in guard.iter_mut().zip([displacement, velocity, acceleration]) {
                *sum = *sum + Dd::from_f64(value.abs());
            }
        }
    }
    if guard.iter().any(|s| !s.to_f64().is_finite()) {
        return Err(Error::NumericalRange);
    }
    Ok(PreloadedSpectrum {
        angular_frequency: omega,
        positive_turning_displacement: amplitude,
        negative_turning_displacement: Scaled::product(&[amplitude, opposite]).value()?,
        displacement_dc: Scaled::product(&[amplitude, coefficients[0]]).value()?,
        cosine_amplitudes: q,
        velocity_sine_amplitudes: v,
        acceleration_cosine_amplitudes: a,
        guard: PreloadedSpectrumGuard {
            last_harmonic: 2 * order,
            displacement_l1: guard[0].to_f64(),
            velocity_l1: guard[1].to_f64(),
            acceleration_l1: guard[2].to_f64(),
        },
        intervals: 0,
        relative_frequency_change: 0.0,
        normalized_kinematic_change: 0.0,
    })
}

fn symmetric(k: f64, beta: f64, amplitude: f64, order: usize) -> Result<PreloadedSpectrum, Error> {
    let exact = duffing_kinematic_spectrum_exact(det::sqrt(k), beta, amplitude, order)
        .map_err(|_| Error::NumericalRange)?;
    // Use the exact owner's physical derivatives directly: dividing by A would
    // lose representable derivatives when a displacement has rounded to zero.
    let mut result = PreloadedSpectrum {
        angular_frequency: exact.displacement.angular_frequency,
        positive_turning_displacement: amplitude,
        negative_turning_displacement: -amplitude,
        displacement_dc: 0.0,
        cosine_amplitudes: vec![0.0; order],
        velocity_sine_amplitudes: vec![0.0; order],
        acceleration_cosine_amplitudes: vec![0.0; order],
        guard: PreloadedSpectrumGuard {
            last_harmonic: 2 * order,
            displacement_l1: 0.0,
            velocity_l1: 0.0,
            acceleration_l1: 0.0,
        },
        intervals: 0,
        relative_frequency_change: 0.0,
        normalized_kinematic_change: 0.0,
    };
    for j in 0..order {
        let h = 2 * j + 1;
        if h <= order {
            result.cosine_amplitudes[h - 1] = exact.displacement.cosine_amplitudes[j];
            result.velocity_sine_amplitudes[h - 1] = exact.velocity_sine_amplitudes[j];
            result.acceleration_cosine_amplitudes[h - 1] = exact.acceleration_cosine_amplitudes[j];
        } else {
            result.guard.displacement_l1 += exact.displacement.cosine_amplitudes[j].abs();
            result.guard.velocity_l1 += exact.velocity_sine_amplitudes[j].abs();
            result.guard.acceleration_l1 += exact.acceleration_cosine_amplitudes[j].abs();
        }
    }
    if [
        result.guard.displacement_l1,
        result.guard.velocity_l1,
        result.guard.acceleration_l1,
    ]
    .iter()
    .any(|v| !v.is_finite())
    {
        return Err(Error::NumericalRange);
    }
    Ok(result)
}

/// Carry physical unit exponents through products before rounding. This also
/// compares the convexity margin without overflowing alpha² or k*beta.
struct Scaled {
    mantissa: Qd,
    exponent: i32,
}
impl Scaled {
    fn product(factors: &[f64]) -> Self {
        let mut result = Self {
            mantissa: Qd::from_f64(1.0),
            exponent: 0,
        };
        for &factor in factors {
            if factor == 0.0 {
                return Self {
                    mantissa: Qd::ZERO,
                    exponent: 0,
                };
            }
            let (value, correction) = if factor.abs() < f64::MIN_POSITIVE {
                (factor.abs() * 18_014_398_509_481_984.0, -54)
            } else {
                (factor.abs(), 0)
            };
            let bits = value.to_bits();
            result.exponent += ((bits >> 52) & 0x7ff) as i32 - 1023 + correction;
            let mantissa = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52));
            result.mantissa = result.mantissa * Qd::from_f64(mantissa.copysign(factor));
        }
        result
    }
    fn at(&self, exponent: i32) -> Qd {
        self.mantissa.scale_power_of_two(self.exponent - exponent)
    }
    fn value(&self) -> Result<f64, Error> {
        let value = self.at(0).to_f64();
        if value.is_finite() {
            Ok(value)
        } else {
            Err(Error::NumericalRange)
        }
    }
    fn coefficient(&self) -> Result<Qd, Error> {
        let value = self.at(0);
        if !value.is_finite() || (value.is_zero() && !self.mantissa.is_zero()) {
            Err(Error::NumericalRange)
        } else {
            Ok(value)
        }
    }
}
