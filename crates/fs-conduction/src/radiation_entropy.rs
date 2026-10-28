//! Photon entropy of explicitly gray, diffuse, unpolarized thermal emission.
//!
//! Independent chaotic photon modes have entropy k_B*g(n), where
//! g(n)=(1+n)ln(1+n)-n ln(n). A gray emitter into vacuum has
//! n=epsilon/(exp(hc/(lambda k_B T))-1). Entropy is NOT epsilon times
//! blackbody entropy, nor emitted energy divided by material temperature.
//! No coherent beam, polarization, participating medium, reflected illumination
//! or radiative back-reaction is inferred from an energy-only observation.
//!
//! Spectra use W/(m² sr nm) and W/(m² sr K nm). Hemispherical flux is pi
//! times isotropic radiance. The full-spectrum integral uses an infinite-domain
//! change of variable, an explicit work budget and a numerical error ESTIMATE.
//! See Kabelac and Conrad, Entropy 14 (2012), 717–735, doi:10.3390/e14040717.

use fs_math::det;

const H: f64 = 6.626_070_15e-34;
const C: f64 = 299_792_458.0;
const K: f64 = 1.380_649e-23;
const C1_NM: f64 = 2.0 * H * C * C * 1.0e36;
const C2_NM_K: f64 = H * C / K * 1.0e9;
const C3_NM: f64 = 2.0 * C * K * 1.0e27;
const PI: f64 = core::f64::consts::PI;
const BLACK_INTEGRAL: f64 = 4.0 * PI * PI * PI * PI / 45.0;

/// Invalid physical coordinates, unrepresentable arithmetic or exhausted work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RadiationEntropyError {
    /// Temperature must be finite and strictly positive, in kelvin.
    Temperature,
    /// Gray emissivity must be finite and in [0,1].
    Emissivity,
    /// Vacuum wavelength must be finite and strictly positive, in nanometres.
    Wavelength,
    /// Relative tolerance or quadrature work bound is invalid.
    Quadrature,
    /// The requested integral did not converge within its evaluation budget.
    BudgetExceeded { evaluations: usize },
    /// Cancellation preserves the caller's original physical state.
    Cancelled,
    /// A finite physical result cannot be represented in binary64.
    Unrepresentable,
}
impl core::fmt::Display for RadiationEntropyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "photon entropy refused: {self:?}")
    }
}
impl std::error::Error for RadiationEntropyError {}

/// One emitted spectrum, before any exposure, color transform or geometry weight.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonEmissionSpectrum {
    /// Emitted energy radiance, W/(m² sr nm).
    pub energy_w_m2_sr_nm: f64,
    /// Emitted photon entropy radiance, W/(m² sr K nm).
    pub entropy_w_m2_sr_k_nm: f64,
}

/// Declared emission model, not a new material measurement or mutable heat store.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrayPhotonEmission {
    temperature_k: f64,
    emissivity: f64,
}

/// Relative accuracy for a positive dimensionless integral, separate from model error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonEntropyQuadrature {
    /// In [64*EPSILON, 0.1). This controls a Richardson error estimate, not a bound.
    pub relative_tolerance: f64,
    /// Maximum integrand evaluations, in 33..=1_048_577; no heap allocation.
    pub maximum_evaluations: usize,
}
impl Default for PhotonEntropyQuadrature {
    fn default() -> Self {
        Self { relative_tolerance: 1e-9, maximum_evaluations: 65_537 }
    }
}

/// Full-spectrum fluxes per area of ONE emitting side, with no added energy debit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonEmissionFlux {
    /// Analytic epsilon*sigma*T^4, W/m².
    pub energy_w_m2: f64,
    /// Integrated physical photon entropy, W/(m² K).
    pub entropy_w_m2_k: f64,
    /// Estimated absolute quadrature error in entropy, W/(m² K), not certified.
    pub estimated_entropy_error_w_m2_k: f64,
    /// Integrand evaluations. Exact black or nonemitting limits need none.
    pub evaluations: usize,
}

impl GrayPhotonEmission {
    /// Explicitly admit the wavelength-independent gray/diffuse approximation.
    /// A hemispherical-total emissivity claim alone does not establish it.
    pub fn new(temperature_k: f64, emissivity: f64) -> Result<Self, RadiationEntropyError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(RadiationEntropyError::Temperature);
        }
        if !emissivity.is_finite() || !(0.0..=1.0).contains(&emissivity) {
            return Err(RadiationEntropyError::Emissivity);
        }
        Ok(Self { temperature_k, emissivity })
    }

    /// Evaluate both spectra without allocating or first forming an underflowing
    /// photon occupation. Energy may underflow independently of entropy in a
    /// deep Wien tail. A zero spectrum means floating-point underflow, not an
    /// exact absence of photons, unless emissivity is exactly zero.
    pub fn spectrum_nm(self, wavelength_nm: f64)
        -> Result<PhotonEmissionSpectrum, RadiationEntropyError>
    {
        if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
            return Err(RadiationEntropyError::Wavelength);
        }
        if self.emissivity == 0.0 {
            return Ok(PhotonEmissionSpectrum { energy_w_m2_sr_nm: 0.0, entropy_w_m2_sr_k_nm: 0.0 });
        }
        let log_lambda = det::ln(wavelength_nm);
        let log_x = det::ln(C2_NM_K) - log_lambda - det::ln(self.temperature_k);
        let raw_x = (C2_NM_K / wavelength_nm) / self.temperature_k;
        let x = if raw_x.is_finite() && raw_x > 0.0 { raw_x } else { det::exp(log_x) };
        let log_n = if x == 0.0 {
            det::ln(self.emissivity) - log_x
        } else {
            log_occupation(x, det::ln(self.emissivity))
        };
        Ok(PhotonEmissionSpectrum {
            energy_w_m2_sr_nm: represented_exp(det::ln(C1_NM) - 5.0*log_lambda + log_n)?,
            entropy_w_m2_sr_k_nm: represented_exp(det::ln(C3_NM) - 4.0*log_lambda + log_mode_entropy(log_n))?,
        })
    }

    /// Integrate emitted entropy over ALL wavelengths, not just the render band.
    /// The integrand is normalized by epsilon before quadrature, keeping very
    /// weak gray emission resolvable. Refinement reuses previous evaluations and
    /// polls cancellation at most every 64 new nodes. No state is mutated.
    pub fn bolometric(self, controls: PhotonEntropyQuadrature, mut cancelled: impl FnMut() -> bool)
        -> Result<PhotonEmissionFlux, RadiationEntropyError>
    {
        poll(&mut cancelled)?;
        if !controls.relative_tolerance.is_finite()
            || !(64.0*f64::EPSILON..0.1).contains(&controls.relative_tolerance)
            || !(33..=1_048_577).contains(&controls.maximum_evaluations)
        { return Err(RadiationEntropyError::Quadrature); }
        if self.emissivity == 0.0 {
            return Ok(PhotonEmissionFlux { energy_w_m2: 0.0, entropy_w_m2_k: 0.0,
                estimated_entropy_error_w_m2_k: 0.0, evaluations: 0 });
        }
        let log_e = det::ln(self.emissivity);
        let log_t = det::ln(self.temperature_k);
        let log_sigma = det::ln(crate::radiation::STEFAN_BOLTZMANN_W_M2_K4);
        let energy = represented_exp(log_e + log_sigma + 4.0*log_t)?;
        let log_scale = log_e + log_sigma + det::ln(4.0/3.0) + 3.0*log_t - det::ln(BLACK_INTEGRAL);
        if self.emissivity == 1.0 {
            let entropy = represented_exp(log_sigma + det::ln(4.0/3.0) + 3.0*log_t)?;
            poll(&mut cancelled)?;
            return Ok(PhotonEmissionFlux { energy_w_m2: energy, entropy_w_m2_k: entropy,
                estimated_entropy_error_w_m2_k: 0.0, evaluations: 0 });
        }
        let mut n = 8usize;
        let mut sum = Sum::default();
        sum.add(0.5*(integrand(0.0, log_e)? + integrand(1.0, log_e)?))?;
        for i in 1..n { sum.add(integrand(i as f64/n as f64, log_e)?)?; }
        let mut trapezoid = sum.total()?/n as f64;
        let mut previous_simpson = None;
        loop {
            poll(&mut cancelled)?;
            let next_n = 2*n;
            if next_n + 1 > controls.maximum_evaluations {
                return Err(RadiationEntropyError::BudgetExceeded { evaluations: n + 1 });
            }
            let mut added = Sum::default();
            for i in 0..n {
                if i % 64 == 0 { poll(&mut cancelled)?; }
                added.add(integrand((2*i+1) as f64/next_n as f64, log_e)?)?;
            }
            let refined = 0.5*trapezoid + added.total()?/next_n as f64;
            let simpson = finite((4.0*refined-trapezoid)/3.0)?;
            if let Some(previous) = previous_simpson {
                let error = f64::abs(simpson-previous)/15.0;
                if simpson > 0.0 && error <= controls.relative_tolerance*simpson {
                    let entropy = represented_exp(log_scale + det::ln(simpson))?;
                    let estimate = if error == 0.0 { 0.0 }
                        else { represented_exp(log_scale + det::ln(error))? };
                    poll(&mut cancelled)?;
                    return Ok(PhotonEmissionFlux { energy_w_m2: energy, entropy_w_m2_k: entropy,
                        estimated_entropy_error_w_m2_k: estimate, evaluations: next_n+1 });
                }
            }
            previous_simpson = Some(simpson);
            trapezoid = refined;
            n = next_n;
        }
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RadiationEntropyError> {
    if cancelled() { Err(RadiationEntropyError::Cancelled) } else { Ok(()) }
}
fn finite(x: f64) -> Result<f64, RadiationEntropyError> {
    if x.is_finite() { Ok(x) } else { Err(RadiationEntropyError::Unrepresentable) }
}
fn represented_exp(x: f64) -> Result<f64, RadiationEntropyError> {
    let value = det::exp(x);
    if value.is_finite() && value >= 0.0 { Ok(value) }
    else { Err(RadiationEntropyError::Unrepresentable) }
}
// All callers have x in [0,1]; compensated log1p retains sub-ulp increments.
fn log1p(x: f64) -> f64 {
    let y = 1.0+x;
    if y == 1.0 { x } else { det::ln(y)*(x/(y-1.0)) }
}
fn log_occupation(x: f64, log_e: f64) -> f64 {
    if x >= 50.0 { log_e-x }
    else { log_e-det::ln(det::expm1(x)) }
}
// Evaluate log(g(n)) from log(n), never subtracting two O(n log n) terms.
fn log_mode_entropy(log_n: f64) -> f64 {
    if log_n == f64::NEG_INFINITY { return f64::NEG_INFINITY; }
    if log_n < -18.0 {
        let n = det::exp(log_n);
        log_n + det::ln(1.0-log_n + n*(0.5-n/6.0))
    } else if log_n < 0.0 {
        let n = det::exp(log_n);
        det::ln((1.0+n)*log1p(n)-n*log_n)
    } else {
        let inverse = det::exp(-log_n);
        let l = log1p(inverse);
        det::ln(log_n+l+if inverse == 0.0 { 1.0 } else { l/inverse })
    }
}
// x=(t/(1-t))^4 maps the COMPLETE Planck integral to [0,1]. The quartic
// map smooths the logarithmic low-frequency endpoint even as epsilon tends
// to zero; this is essential to the Simpson error estimate. Both limits vanish.
// Normalizing by emissivity avoids underflow of weak emission.
fn integrand(t: f64, log_e: f64) -> Result<f64, RadiationEntropyError> {
    if t == 0.0 || t == 1.0 { return Ok(0.0); }
    let ratio = t/(1.0-t);
    let square = ratio*ratio;
    let x = square*square;
    let log_g = log_mode_entropy(log_occupation(x, log_e));
    represented_exp(det::ln(4.0) + 11.0*det::ln(t) - 13.0*det::ln(1.0-t) + log_g-log_e)
}
#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, x: f64) -> Result<(), RadiationEntropyError> {
        let next = finite(self.value+x)?;
        self.correction = finite(self.correction + if self.value.abs() >= x.abs() {
            (self.value-next)+x
        } else { (x-next)+self.value })?;
        self.value = next;
        Ok(())
    }
    fn total(self) -> Result<f64, RadiationEntropyError> { finite(self.value+self.correction) }
}

#[cfg(test)]
#[path = "radiation_entropy/tests.rs"]
mod tests;
