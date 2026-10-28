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

mod ambient;
pub use ambient::MaterialRadiationEntropyExchange;
/// Combined reflected/emitted photon modes in an explicitly black isotropic bath.
pub mod bath;

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
    log_temperature: f64,
    log_emissivity: f64,
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
    /// Integrand evaluations used to prepare this factor, not repeated for a
    /// temperature-only evaluation. Exact black or nonemitting limits need none.
    pub evaluations: usize,
}

impl GrayPhotonEmission {
    /// Prepare the wavelength-independent gray/diffuse approximation, caching
    /// temperature/emissivity logarithms for repeated allocation-free spectra.
    /// A hemispherical-total emissivity claim alone does not establish it.
    pub fn new(temperature_k: f64, emissivity: f64) -> Result<Self, RadiationEntropyError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(RadiationEntropyError::Temperature);
        }
        if !emissivity.is_finite() || !(0.0..=1.0).contains(&emissivity) {
            return Err(RadiationEntropyError::Emissivity);
        }
        Ok(Self { temperature_k, emissivity,
            log_temperature: det::ln(temperature_k),
            log_emissivity: if emissivity == 0.0 { f64::NEG_INFINITY } else { det::ln(emissivity) },
        })
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
        let log_x = det::ln(C2_NM_K) - log_lambda - self.log_temperature;
        let raw_x = (C2_NM_K / wavelength_nm) / self.temperature_k;
        let x = if raw_x.is_finite() && raw_x > 0.0 { raw_x } else { det::exp(log_x) };
        let log_n = if x == 0.0 {
            self.log_emissivity - log_x
        } else {
            log_occupation(x, self.log_emissivity)
        };
        Ok(PhotonEmissionSpectrum {
            energy_w_m2_sr_nm: represented_exp(det::ln(C1_NM) - 5.0*log_lambda + log_n)?,
            entropy_w_m2_sr_k_nm: represented_exp(det::ln(C3_NM) - 4.0*log_lambda + log_mode_entropy(log_n))?,
        })
    }

    /// Integrate emitted entropy over ALL wavelengths, not just the render band.
    /// For many temperatures at unchanged emissivity, prepare once with
    /// `prepare_bolometric` and reuse its temperature-only evaluation.
    pub fn bolometric(self, controls: PhotonEntropyQuadrature, mut cancelled: impl FnMut() -> bool)
        -> Result<PhotonEmissionFlux, RadiationEntropyError>
    {
        let prepared = self.prepare_bolometric(controls, &mut cancelled)?;
        let result = prepared.flux_at_log_temperature(self.log_temperature)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Prepare the temperature-independent full-spectrum entropy factor once.
    /// Quadrature is normalized by emissivity to retain weak gray emission;
    /// refinement reuses prior evaluations, allocates nothing and checks
    /// cancellation at most every 64 new nodes. The result may be reused at
    /// arbitrary positive temperatures under the SAME explicitly gray model.
    /// Its reuse does not extend any material source's temperature validity.
    pub fn prepare_bolometric(self, controls: PhotonEntropyQuadrature,
        mut cancelled: impl FnMut() -> bool) -> Result<PreparedPhotonEntropy, RadiationEntropyError>
    {
        poll(&mut cancelled)?;
        if !controls.relative_tolerance.is_finite()
            || !(64.0*f64::EPSILON..0.1).contains(&controls.relative_tolerance)
            || !(33..=1_048_577).contains(&controls.maximum_evaluations)
        { return Err(RadiationEntropyError::Quadrature); }
        if self.emissivity == 0.0 {
            return Ok(PreparedPhotonEntropy { log_emissivity: f64::NEG_INFINITY,
                normalized_integral: 0.0, estimated_error: 0.0, evaluations: 0 });
        }
        let log_e = self.log_emissivity;
        if self.emissivity == 1.0 {
            poll(&mut cancelled)?;
            return Ok(PreparedPhotonEntropy { log_emissivity: 0.0,
                normalized_integral: BLACK_INTEGRAL, estimated_error: 0.0, evaluations: 0 });
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
                    poll(&mut cancelled)?;
                    return Ok(PreparedPhotonEntropy { log_emissivity: log_e,
                        normalized_integral: simpson, estimated_error: error, evaluations: next_n+1 });
                }
            }
            previous_simpson = Some(simpson);
            trapezoid = refined;
            n = next_n;
        }
    }
}

/// Prepared full-spectrum coefficient for one unchanged gray emissivity.
/// Stores no thermal history and allocates nothing. A temperature sweep needs
/// no further quadrature, material query or cancellation-controlled iteration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedPhotonEntropy {
    log_emissivity: f64,
    normalized_integral: f64,
    estimated_error: f64,
    evaluations: usize,
}
impl PreparedPhotonEntropy {
    /// Evaluate analytic T^4 energy and T^3 photon entropy scaling. The caller
    /// remains responsible for the original material's validity at this T.
    pub fn flux_at_temperature(self, temperature_k: f64)
        -> Result<PhotonEmissionFlux, RadiationEntropyError>
    {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(RadiationEntropyError::Temperature);
        }
        self.flux_at_log_temperature(det::ln(temperature_k))
    }
    fn flux_at_log_temperature(self, log_t: f64)
        -> Result<PhotonEmissionFlux, RadiationEntropyError>
    {
        if self.normalized_integral == 0.0 {
            return Ok(PhotonEmissionFlux { energy_w_m2: 0.0, entropy_w_m2_k: 0.0,
                estimated_entropy_error_w_m2_k: 0.0, evaluations: 0 });
        }
        let log_sigma = det::ln(crate::radiation::STEFAN_BOLTZMANN_W_M2_K4);
        let energy = represented_exp(self.log_emissivity+log_sigma+4.0*log_t)?;
        let log_base = self.log_emissivity+log_sigma+det::ln(4.0/3.0)+3.0*log_t;
        let entropy = if self.log_emissivity == 0.0 {
            represented_exp(log_base)?
        } else {
            represented_exp(log_base-det::ln(BLACK_INTEGRAL)+det::ln(self.normalized_integral))?
        };
        let error = if self.estimated_error == 0.0 { 0.0 } else {
            represented_exp(log_base-det::ln(BLACK_INTEGRAL)+det::ln(self.estimated_error))?
        };
        Ok(PhotonEmissionFlux { energy_w_m2: energy, entropy_w_m2_k: entropy,
            estimated_entropy_error_w_m2_k: error, evaluations: self.evaluations })
    }
}

impl crate::AmbientRadiationPatch {
    /// Explicit gray/diffuse photon observation using this original patch's
    /// selected emissivity, actual surface temperature and black bath.
    /// A total emissivity measurement alone does not establish spectral grayness.
    /// The existing heat law admits the surface's source domain; evaluating it
    /// here does not apply heat. The returned closure is not another cooling
    /// debit and does not carry new source authority. Retain this original patch
    /// and rebind after temperature/source changes. Spatially nonuniform traces
    /// require local observations, not entropy at an averaged temperature.
    pub fn gray_photon_exchange(&self, surface_temperature_k: f64)
        -> Result<bath::GrayBlackBathExchange, crate::ConductionError>
    {
        self.heat_flux_w_m2(surface_temperature_k)?;
        GrayPhotonEmission::new(surface_temperature_k,self.emissivity().value())
            .and_then(|emission| emission.in_black_bath(self.ambient_temperature_k()))
            .map_err(|error| crate::ConductionError::Radiation {
                surface: self.region().into(), what: error.to_string(),
                fix: "use a representable, explicitly gray/diffuse surface and black bath".into(),
            })
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

#[cfg(test)]
mod bath_source_tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
        PropertyKey,PropertyValue,Provenance,SelectionPolicy,UncertaintyModel};
    fn patch(ambient: f64) -> crate::AmbientRadiationPatch {
        let dims = crate::EMISSIVITY_DIMS;
        let mut claims = ClaimSet::new();
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(crate::SURFACE_EMISSIVITY_PROPERTY,dims),
            value: PropertyValue::Scalar { value: 0.4,dims },
            validity: ValidityDomain::unconstrained().with("T",250.0,500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: Vec::new(),
            provenance: Provenance { source: "analytic photon bath fixture".into(),
                license: "internal-test-use".into(), artifact: None },
        }).unwrap();
        let card = MaterialCard::assemble(MaterialStateId { chemistry: "gray-fixture".into(),
            phase: "solid".into(),process: "test".into(),revision: 0 },claims,Vec::new()).unwrap();
        let source = crate::SurfaceEmissivity::from_card("surface",&card,350.0,
            SelectionPolicy::SingleClaimOnly).unwrap();
        crate::AmbientRadiationPatch::new("surface",source,ambient).unwrap()
    }
    #[test]
    fn original_sources_and_surface_domains_survive_the_photon_binding() {
        let original = patch(1000.0);
        let saved = original.clone();
        let model = original.gray_photon_exchange(400.0).unwrap();
        assert_eq!(model.surface_temperature_k(),400.0);
        assert_eq!(model.bath_temperature_k(),1000.0);
        assert_eq!(model.emissivity(),original.emissivity().value());
        assert!(model.spectrum_nm(10_000.0).unwrap().heat_out_w_m2_nm<0.0);
        assert!(original.gray_photon_exchange(501.0).is_err());
        assert!(original.gray_photon_exchange(0.0).is_err());
        assert_eq!(original,saved);
    }
    #[test]
    fn equilibrium_patch_has_no_photon_or_material_production() {
        let original = patch(350.0);
        let model = original.gray_photon_exchange(350.0).unwrap();
        assert_eq!(model.spectrum_nm(10_000.0).unwrap().total_generation_w_m2_k_nm,0.0);
        assert_eq!(original.material_entropy_exchange(350.0).unwrap().generation_w_m2_k,0.0);
    }
}
