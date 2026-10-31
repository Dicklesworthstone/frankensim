//! Homogeneous-sphere scattering with physical cross sections and angular data.
//!
//! The nonmagnetic, independent-particle Mie series uses a real, nonabsorbing
//! host index and the particle index `n + i k` with `exp(-i omega t)` time
//! convention. The logarithmic derivative is evaluated downwards; small real
//! arguments use a convergent Riccati--Bessel series rather than cancellation
//! in `sin(x)/x - cos(x)`. No Henyey--Greenstein fit replaces the angular law.
//!
//! Equations: Bohren & Huffman, *Absorption and Scattering of Light by Small
//! Particles*, §§4.4, 4.6; the independently published algorithm and reference
//! values are reproduced in Scott Prahl's [Mie algorithms notes]
//! (https://miepython.readthedocs.io/en/3.0.0/07_algorithm.html).
//! The explicit domain is `1e-4 <= x <= 128`, `|m x| <= 512`. A finite series
//! and a sampled angular-table discrepancy are numerical estimates, not
//! certified error bounds. Nonspherical particles, dependent scattering,
//! absorbing hosts, coherent interference between particles and polarization
//! transport are outside this model.

use super::{MAX_PHASE_KNOTS, TabulatedPhase, VolumeTransportError, checkpoint};
use fs_exec::Cx;
use fs_math::{c64::C64, det};

mod distribution;
mod population;
mod spectrum;
pub use distribution::{LognormalRadiusReport, ParticleSizeDistribution};
pub use population::{
    MAX_MIE_POPULATION_BINS, MiePopulation, MiePopulationBudget, MiePopulationSolution,
};
pub use spectrum::MieOpticalSpectrum;
pub(crate) use spectrum::prepare_visible_spectrum;

/// Maximum admitted size parameter `2 pi n_host r / lambda_vacuum`.
pub const MAX_MIE_SIZE_PARAMETER: f64 = 128.0;

/// Explicit preparation work limits, independent of a rendering sample count.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MieBudget {
    /// Maximum number of retained electric/magnetic multipole pairs.
    pub max_multipoles: usize,
    /// Maximum order from which the logarithmic derivative may descend.
    pub max_recurrence_order: usize,
}

impl Default for MieBudget {
    fn default() -> Self {
        Self {
            max_multipoles: 160,
            max_recurrence_order: 576,
        }
    }
}

/// Bounded model, arithmetic and preparation failures.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MieError {
    /// Radius, indices, density, wavelength or tolerance is invalid.
    Parameters,
    /// Size or internal size lies outside the supported numerical domain.
    Regime,
    /// Explicit recurrence, storage or multipole work budget was exhausted.
    Budget,
    /// The requested angular discrepancy could not fit in the knot budget.
    AngularResolution,
    /// Zero scattering has no normalized angular distribution.
    NoScattering,
    /// Arithmetic, recurrence convergence or passive-energy admission failed.
    Numerical,
    /// Cancellation or the surrounding execution budget refused preparation.
    Execution,
}

impl core::fmt::Display for MieError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sphere scattering refused: {self:?}")
    }
}
impl core::error::Error for MieError {}

/// One homogeneous sphere with explicit SI size and material optical inputs.
///
/// Indices supplied here are constants for this declaration. A dispersive
/// caller constructs the sphere at each wavelength using that material's
/// measured or constitutive optical index. Density is the particle material's
/// density, not the carrier-gas density; transport grids contain dispersed
/// particle mass per unit mixture volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MieSphere {
    radius_m: f64,
    particle_index_n: f64,
    extinction_k: f64,
    host_index_n: f64,
    particle_density_kg_m3: f64,
}

impl MieSphere {
    /// Admit positive finite radius, real indices and particle density, and a
    /// finite nonnegative extinction index (passive particles only).
    pub fn try_new(
        radius_m: f64,
        particle_index_n: f64,
        extinction_k: f64,
        host_index_n: f64,
        particle_density_kg_m3: f64,
    ) -> Result<Self, MieError> {
        if [
            radius_m,
            particle_index_n,
            host_index_n,
            particle_density_kg_m3,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
            || !extinction_k.is_finite()
            || extinction_k < 0.0
        {
            return Err(MieError::Parameters);
        }
        Ok(Self {
            radius_m,
            particle_index_n,
            extinction_k,
            host_index_n,
            particle_density_kg_m3,
        })
    }

    /// Physical particle radius, metres.
    #[must_use]
    pub const fn radius_m(self) -> f64 {
        self.radius_m
    }

    /// Real nonabsorbing host index used by this sphere.
    #[must_use]
    pub const fn host_index_n(self) -> f64 {
        self.host_index_n
    }

    /// Particle mass `4*pi*r^3*rho/3`, kg. Population reductions require
    /// this positive quantity to be representable in addition to cross sections.
    pub fn particle_mass_kg(self) -> Result<f64, MieError> {
        let mass = (4.0 * core::f64::consts::PI / 3.0)
            * (self.radius_m * self.radius_m)
            * (self.particle_density_kg_m3 * self.radius_m);
        if mass.is_finite() && mass > 0.0 {
            Ok(mass)
        } else {
            Err(MieError::Numerical)
        }
    }

    /// Prepare the complete retained Mie series at one vacuum wavelength, nm.
    ///
    /// The truncation `ceil(x + 4.05 cbrt(x) + 2)` follows the conventional
    /// Wiscombe size rule. Two downward starts sixteen orders apart must agree
    /// before coefficients are admitted. Preparation polls cancellation at
    /// least every 32 recurrence or multipole steps.
    pub fn at_wavelength(
        self,
        cx: &Cx<'_>,
        vacuum_nm: f64,
        budget: MieBudget,
    ) -> Result<MieSolution, MieError> {
        poll(cx)?;
        if !vacuum_nm.is_finite() || vacuum_nm <= 0.0 {
            return Err(MieError::Parameters);
        }
        let x = core::f64::consts::TAU * self.host_index_n * (self.radius_m / (vacuum_nm * 1e-9));
        let m = C64::new(
            self.particle_index_n / self.host_index_n,
            self.extinction_k / self.host_index_n,
        );
        let z = m.scale(x);
        if !(1e-4..=MAX_MIE_SIZE_PARAMETER).contains(&x) || !z.abs().is_finite() || z.abs() > 512.0
        {
            return Err(MieError::Regime);
        }
        let count = (x + 4.05 * det::exp(det::ln(x) / 3.0) + 2.0).ceil() as usize;
        let start = count.max(z.abs().ceil() as usize) + 48;
        if count > budget.max_multipoles || start + 16 > budget.max_recurrence_order {
            return Err(MieError::Budget);
        }
        let mut coefficients = Vec::new();
        coefficients
            .try_reserve_exact(count)
            .map_err(|_| MieError::Budget)?;
        // Exact index matching must not manufacture numerical haze.
        if m == C64::ONE {
            return self.finish(vacuum_nm, x, coefficients);
        }
        let derivatives = logarithmic_derivatives(cx, z, count, start)?;
        let check = logarithmic_derivatives(cx, z, count, start + 16)?;
        if derivatives
            .iter()
            .zip(&check)
            .any(|(&a, &b)| !finite(a) || !finite(b) || (a - b).abs() > 2e-12 * (1.0 + b.abs()))
        {
            return Err(MieError::Numerical);
        }
        let (sin_x, cos_x) = det::sin_cos(x);
        let mut psi_previous = sin_x;
        let mut psi = if x < 0.25 {
            small_psi(1, x)
        } else {
            sin_x / x - cos_x
        };
        // Outgoing h^(1) under exp(-i omega t): xi = psi + i x*y_n.
        let mut chi_previous = -cos_x;
        let mut chi = -cos_x / x - sin_x;
        for (index, derivative) in derivatives.into_iter().enumerate() {
            if index % 32 == 0 {
                poll(cx)?;
            }
            let n = (index + 1) as f64;
            let xi = C64::new(psi, chi);
            let xi_previous = C64::new(psi_previous, chi_previous);
            let electric = derivative / m + C64::from_re(n / x);
            let magnetic = derivative * m + C64::from_re(n / x);
            let a =
                (electric.scale(psi) - C64::from_re(psi_previous)) / (electric * xi - xi_previous);
            let b =
                (magnetic.scale(psi) - C64::from_re(psi_previous)) / (magnetic * xi - xi_previous);
            if !finite(a) || !finite(b) {
                return Err(MieError::Numerical);
            }
            coefficients.push((a, b));
            let next_psi = if x < 0.25 {
                small_psi(index + 2, x)
            } else {
                ((2.0 * n + 1.0) / x).mul_add(psi, -psi_previous)
            };
            let next_chi = ((2.0 * n + 1.0) / x).mul_add(chi, -chi_previous);
            (psi_previous, psi) = (psi, next_psi);
            (chi_previous, chi) = (chi, next_chi);
        }
        poll(cx)?;
        self.finish(vacuum_nm, x, coefficients)
    }

    fn finish(
        self,
        vacuum_nm: f64,
        x: f64,
        coefficients: Vec<(C64, C64)>,
    ) -> Result<MieSolution, MieError> {
        let (mut scattering_sum, mut absorption_sum, mut moment_sum) = (0.0, 0.0, 0.0);
        for (index, &(a, b)) in coefficients.iter().enumerate() {
            let n = (index + 1) as f64;
            let power = a.norm_sq() + b.norm_sq();
            let loss = (a.re - a.norm_sq()) + (b.re - b.norm_sq());
            // Roundoff near a passive lossless multipole is allowed; active
            // amplification or a material absorption deficit is not.
            if loss < -128.0 * f64::EPSILON * (a.abs() + b.abs()) {
                return Err(MieError::Numerical);
            }
            scattering_sum += (2.0 * n + 1.0) * power;
            if self.extinction_k > 0.0 {
                absorption_sum += (2.0 * n + 1.0) * loss.max(0.0);
            }
            moment_sum += (2.0 * n + 1.0) / (n * (n + 1.0)) * (a * b.conj()).re;
            if let Some(&(next_a, next_b)) = coefficients.get(index + 1) {
                moment_sum +=
                    n * (n + 2.0) / (n + 1.0) * ((a * next_a.conj()).re + (b * next_b.conj()).re);
            }
        }
        let q_scattering = 2.0 * scattering_sum / (x * x);
        let q_absorption = 2.0 * absorption_sum / (x * x);
        let area = core::f64::consts::PI * self.radius_m * self.radius_m;
        // C / mass = 3 Q / (4 rho r), without an r^3 intermediate.
        let mass_factor = 0.75 / (self.particle_density_kg_m3 * self.radius_m);
        let cross_sections = MieCrossSections {
            scattering_m2: q_scattering * area,
            absorption_m2: q_absorption * area,
            mass_scattering_m2_kg: q_scattering * mass_factor,
            mass_absorption_m2_kg: q_absorption * mass_factor,
            mean_cosine: if scattering_sum > 0.0 {
                2.0 * moment_sum / scattering_sum
            } else {
                0.0
            },
        };
        if [
            area,
            mass_factor,
            cross_sections.scattering_m2,
            cross_sections.absorption_m2,
            cross_sections.mass_scattering_m2_kg,
            cross_sections.mass_absorption_m2_kg,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
            || area == 0.0
            || mass_factor == 0.0
            || (q_scattering > 0.0
                && (cross_sections.scattering_m2 == 0.0
                    || cross_sections.mass_scattering_m2_kg == 0.0))
            || (q_absorption > 0.0
                && (cross_sections.absorption_m2 == 0.0
                    || cross_sections.mass_absorption_m2_kg == 0.0))
            || !(-1.0..=1.0).contains(&cross_sections.mean_cosine)
            || !cross_sections.extinction_m2().is_finite()
            || !cross_sections.transport_m2().is_finite()
        {
            return Err(MieError::Numerical);
        }
        Ok(MieSolution {
            vacuum_nm,
            size_parameter: x,
            coefficients,
            scattering_sum,
            cross_sections,
        })
    }
}

/// SI cross sections and particle-mass coefficients from the same series.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MieCrossSections {
    /// Elastic scattering cross section per particle, m².
    pub scattering_m2: f64,
    /// Absorption cross section per particle, m².
    pub absorption_m2: f64,
    /// Scattering per kilogram of dispersed particle material, m²/kg.
    pub mass_scattering_m2_kg: f64,
    /// Absorption per kilogram of dispersed particle material, m²/kg.
    pub mass_absorption_m2_kg: f64,
    /// Angular first moment of the retained Mie series.
    pub mean_cosine: f64,
}

impl MieCrossSections {
    /// Total beam removal, m²; absorption and redistribution remain distinct.
    #[must_use]
    pub fn extinction_m2(self) -> f64 {
        self.scattering_m2 + self.absorption_m2
    }

    /// Transport cross section, `C_abs + (1-g) C_sca`, m².
    /// This controls the independent-particle radiation-pressure efficiency.
    #[must_use]
    pub fn transport_m2(self) -> f64 {
        self.absorption_m2 + (1.0 - self.mean_cosine) * self.scattering_m2
    }
}

/// Prepared multipole coefficients; angle evaluations allocate no storage.
#[derive(Clone, Debug, PartialEq)]
pub struct MieSolution {
    vacuum_nm: f64,
    size_parameter: f64,
    coefficients: Vec<(C64, C64)>,
    scattering_sum: f64,
    cross_sections: MieCrossSections,
}

impl MieSolution {
    /// Vacuum wavelength used for this preparation, nm.
    #[must_use]
    pub const fn wavelength_nm(&self) -> f64 {
        self.vacuum_nm
    }
    /// Host-medium size parameter.
    #[must_use]
    pub const fn size_parameter(&self) -> f64 {
        self.size_parameter
    }
    /// Number of retained electric/magnetic multipole pairs.
    #[must_use]
    pub fn multipoles(&self) -> usize {
        self.coefficients.len()
    }
    /// Cross sections and asymmetry from the retained coefficients.
    #[must_use]
    pub const fn cross_sections(&self) -> MieCrossSections {
        self.cross_sections
    }

    /// Normalized unpolarized Mie phase density, sr⁻¹, evaluated directly from
    /// the retained series at the propagation-direction cosine.
    pub fn eval_cos(&self, cosine: f64) -> Result<f64, MieError> {
        if !cosine.is_finite() || cosine.abs() > 1.0 + 2e-10 {
            return Err(MieError::Parameters);
        }
        if self.scattering_sum == 0.0 {
            return Err(MieError::NoScattering);
        }
        let mu = cosine.clamp(-1.0, 1.0);
        let (mut pi_previous, mut pi) = (0.0, 1.0);
        let (mut s1, mut s2) = (C64::ZERO, C64::ZERO);
        for (index, &(a, b)) in self.coefficients.iter().enumerate() {
            let n = (index + 1) as f64;
            let tau = n * mu * pi - (n + 1.0) * pi_previous;
            let factor = (2.0 * n + 1.0) / (n * (n + 1.0));
            s1 = s1 + (a.scale(pi) + b.scale(tau)).scale(factor);
            s2 = s2 + (a.scale(tau) + b.scale(pi)).scale(factor);
            let next = ((2.0 * n + 1.0) * mu * pi - (n + 1.0) * pi_previous) / n;
            (pi_previous, pi) = (pi, next);
        }
        let density =
            (s1.norm_sq() + s2.norm_sq()) / (4.0 * core::f64::consts::PI * self.scattering_sum);
        if density.is_finite() && density >= 0.0 {
            Ok(density)
        } else {
            Err(MieError::Numerical)
        }
    }

    /// Prepare a normalized piecewise-linear phase and its matched sampler.
    ///
    /// An initial angular mesh resolves the retained polynomial degree; each
    /// interval is probed at its quarter points and midpoint. Refinement stops
    /// only when sampled pointwise discrepancy, normalization error and first
    /// moment error are all below `tolerance`. The discrepancy divides by the
    /// larger local density or `1/(4 pi)`. It is an observed interpolation
    /// error, not a rigorous supremum bound. No wavelength interpolation is
    /// performed: use one such report per spectral knot in a spectral law.
    pub fn phase_table(
        &self,
        cx: &Cx<'_>,
        tolerance: f64,
        max_knots: usize,
    ) -> Result<MiePhaseTable, MieError> {
        prepare_phase_table(self, cx, tolerance, max_knots)
    }
}

/// Shared preparation seam for single spheres, populations, and material
/// adapters. Only preparation evaluates the full angular series; transport
/// uses the native tabulated phase and its matched sampler.
pub(crate) trait MieAngularSolution {
    fn cross_sections(&self) -> MieCrossSections;
    fn multipoles(&self) -> usize;
    fn eval_cos(&self, cosine: f64) -> Result<f64, MieError>;

    fn phase_table(
        &self,
        cx: &Cx<'_>,
        tolerance: f64,
        max_knots: usize,
    ) -> Result<MiePhaseTable, MieError> {
        prepare_phase_table(self, cx, tolerance, max_knots)
    }
}

impl MieAngularSolution for MieSolution {
    fn cross_sections(&self) -> MieCrossSections {
        self.cross_sections()
    }
    fn multipoles(&self) -> usize {
        self.multipoles()
    }
    fn eval_cos(&self, cosine: f64) -> Result<f64, MieError> {
        self.eval_cos(cosine)
    }
}

fn prepare_phase_table<S: MieAngularSolution + ?Sized>(
    solution: &S,
    cx: &Cx<'_>,
    tolerance: f64,
    max_knots: usize,
) -> Result<MiePhaseTable, MieError> {
    poll(cx)?;
    if !tolerance.is_finite() || !(0.0..=0.1).contains(&tolerance) || tolerance == 0.0 {
        return Err(MieError::Parameters);
    }
    if solution.cross_sections().scattering_m2 == 0.0 {
        return Err(MieError::NoScattering);
    }
    let intervals = (4 * solution.multipoles()).max(16);
    if max_knots > MAX_PHASE_KNOTS || max_knots <= intervals {
        return Err(MieError::AngularResolution);
    }
    let mut knots = Vec::new();
    knots
        .try_reserve_exact(max_knots)
        .map_err(|_| MieError::Budget)?;
    for index in 0..=intervals {
        poll(cx)?;
        let cosine = if index == 0 {
            -1.0
        } else if index == intervals {
            1.0
        } else {
            -det::cos(core::f64::consts::PI * index as f64 / intervals as f64)
        };
        knots.push((cosine, solution.eval_cos(cosine)?));
    }
    loop {
        let mut next = Vec::new();
        next.try_reserve_exact(max_knots)
            .map_err(|_| MieError::Budget)?;
        next.push(knots[0]);
        let (mut maximum_error, mut normalization) = (0.0_f64, 0.0);
        for pair in knots.windows(2) {
            poll(cx)?;
            let [(lo, left), (hi, right)] = [pair[0], pair[1]];
            let mut interval_error = 0.0_f64;
            for fraction in [0.25, 0.5, 0.75] {
                let actual = solution.eval_cos((hi - lo).mul_add(fraction, lo))?;
                let interpolated = left * (1.0 - fraction) + right * fraction;
                interval_error = interval_error.max(
                    (actual - interpolated).abs()
                        / actual
                            .max(interpolated)
                            .max(1.0 / (4.0 * core::f64::consts::PI)),
                );
            }
            maximum_error = maximum_error.max(interval_error);
            normalization += core::f64::consts::PI * (hi - lo) * (left + right);
            if interval_error > tolerance * 0.5 {
                if next.len() + 1 >= max_knots {
                    return Err(MieError::AngularResolution);
                }
                let middle = (hi - lo).mul_add(0.5, lo);
                next.push((middle, solution.eval_cos(middle)?));
            }
            if next.len() == max_knots {
                return Err(MieError::AngularResolution);
            }
            next.push((hi, right));
        }
        if next.len() != knots.len() {
            knots = next;
            continue;
        }
        let phase = TabulatedPhase::try_new(&knots).map_err(|_| MieError::Numerical)?;
        let normalization_error = (normalization - 1.0).abs();
        let mean_cosine_error = (phase.mean_cosine() - solution.cross_sections().mean_cosine).abs();
        // Normalization changes every ordinate. Include that change in
        // the retained pointwise diagnostic rather than hiding it.
        let normalized_error = maximum_error + normalization_error / normalization;
        if normalized_error <= tolerance && mean_cosine_error <= tolerance {
            poll(cx)?;
            return Ok(MiePhaseTable {
                phase,
                normalization_error,
                mean_cosine_error,
                max_sampled_relative_error: normalized_error,
            });
        }
        // A global normalization/moment failure can survive local tests;
        // refine every interval instead of accepting a normalized alias.
        if 2 * knots.len() - 1 > max_knots {
            return Err(MieError::AngularResolution);
        }
        let mut refined = Vec::new();
        refined
            .try_reserve_exact(max_knots)
            .map_err(|_| MieError::Budget)?;
        refined.push(knots[0]);
        for pair in knots.windows(2) {
            poll(cx)?;
            let middle = (pair[1].0 - pair[0].0).mul_add(0.5, pair[0].0);
            refined.push((middle, solution.eval_cos(middle)?));
            refined.push(pair[1]);
        }
        knots = refined;
    }
}

/// Retained angular interpolant plus observed discrepancies against Mie.
#[derive(Clone, Debug, PartialEq)]
pub struct MiePhaseTable {
    /// Native normalized phase table; no HG approximation is involved.
    pub phase: TabulatedPhase,
    /// Absolute deviation of the raw angular integral from unity.
    pub normalization_error: f64,
    /// Absolute error in the normalized interpolant's mean cosine.
    pub mean_cosine_error: f64,
    /// Maximum sampled relative interpolation error, including normalization.
    pub max_sampled_relative_error: f64,
}

fn logarithmic_derivatives(
    cx: &Cx<'_>,
    z: C64,
    count: usize,
    start: usize,
) -> Result<Vec<C64>, MieError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| MieError::Budget)?;
    result.resize(count, C64::ZERO);
    let mut derivative = C64::ZERO;
    for n in (1..=start).rev() {
        if n % 32 == 0 {
            poll(cx)?;
        }
        let ratio = C64::from_re(n as f64) / z;
        derivative = ratio - (derivative + ratio).recip();
        if n > 1 && n - 1 <= count {
            result[n - 2] = derivative;
        }
    }
    Ok(result)
}

fn small_psi(n: usize, x: f64) -> f64 {
    let mut leading = x;
    for j in 1..=n {
        leading *= x / (2 * j + 1) as f64;
    }
    let (mut sum, mut term) = (1.0, 1.0);
    for k in 1..=16 {
        term *= -x * x / (2 * k * (2 * n + 2 * k + 1)) as f64;
        sum += term;
        if term.abs() < f64::EPSILON * sum.abs() {
            break;
        }
    }
    leading * sum
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn poll(cx: &Cx<'_>) -> Result<(), MieError> {
    checkpoint(cx).map_err(|_: VolumeTransportError| MieError::Execution)
}

#[cfg(test)]
mod tests;
