//! Original ordinary-water refraction from temperature, density and wavelength.
//!
//! Implements equation (1) and Table 1 of IAPWS R9-97, using ITS-90 temperature,
//! mass density and VACUUM wavelength. Attribution: International Association
//! for the Properties of Water and Steam, September 1997,
//! <https://iapws.org/documents/release/Rindex.download>.
//!
//! Only the endorsed 200--1100 nm, 261.15--773.15 K and 0--1060 kg/m3 domain
//! is admitted. The release's separately discussed infrared extrapolation is
//! not enabled. Optical uncertainty remains the condition-dependent estimate
//! in its Table 2; analytic derivatives do not establish smaller error bars.
//!
//! This is an absolute, scalar phase-index law for homogeneous ordinary water
//! or steam. A supplied T/rho point does not establish phase stability, purity
//! or a density closure. Neither two-phase bulk densities, ice, brine nor
//! arbitrary mixtures inherit this law. Absorption, scattering and birefringence
//! are separate properties and are not inferred from a real refractive index.
//! Original fluid solvers own density: no pressure-to-density EOS is hidden in
//! the optical evaluator. Temperature derivatives here hold density fixed.

use fs_blake3::{ContentHash, DomainHasher};
use fs_math::det;
use fs_thermochem::water::if97::If97State;

/// Original released coefficients and domain, not a fitted local approximation.
pub const SOURCE_RELEASE: &str = "IAPWS R9-97, equation (1) and Table 1, ITS-90";
pub const SOURCE_URL: &str = "https://iapws.org/documents/release/Rindex.download";
/// Inclusive original source support, in vacuum nanometres.
pub const WAVELENGTH_RANGE_NM: [f64; 2] = [200.0, 1100.0];
/// Inclusive original source support, in kelvin on ITS-90.
pub const TEMPERATURE_RANGE_K: [f64; 2] = [261.15, 773.15];
/// Inclusive original source support, in kg/m3. Zero is the vacuum limit.
pub const DENSITY_RANGE_KG_M3: [f64; 2] = [0.0, 1060.0];

const T_STAR_K: f64 = 273.15;
const RHO_STAR_KG_M3: f64 = 1000.0;
const LAMBDA_STAR_NM: f64 = 589.0;
const LAMBDA_UV: f64 = 0.229_202_0;
const LAMBDA_IR: f64 = 5.432_937;
const A: [f64; 8] = [
    0.244_257_733,
    9.746_344_76e-3,
    -3.732_349_96e-3,
    2.686_784_72e-4,
    1.589_205_70e-3,
    2.459_342_59e-3,
    0.900_704_920,
    -1.666_262_19e-2,
];

/// Source-domain or numerical refusal. No input is clipped to a source endpoint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaterRefractionError {
    WavelengthOutsideRange,
    TemperatureOutsideRange,
    DensityOutsideRange,
    NumericalRange,
}

impl core::fmt::Display for WaterRefractionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "IAPWS water refraction: {self:?}")
    }
}
impl core::error::Error for WaterRefractionError {}

/// The one immutable released law. Wavelength preparation and state sampling
/// use fixed storage and allocate nothing.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IapwsWaterRefraction;

impl IapwsWaterRefraction {
    /// Identity of the original law, its coefficient bits, domain and evaluator.
    /// This identity does not certify the state provider or its material purity.
    pub fn source_identity() -> ContentHash {
        let mut hash = DomainHasher::new("fs-material:iapws-water-refraction:v1");
        hash.update(SOURCE_RELEASE.as_bytes());
        for value in A
            .into_iter()
            .chain([
                T_STAR_K,
                RHO_STAR_KG_M3,
                LAMBDA_STAR_NM,
                LAMBDA_UV,
                LAMBDA_IR,
            ])
            .chain(WAVELENGTH_RANGE_NM)
            .chain(TEMPERATURE_RANGE_K)
            .chain(DENSITY_RANGE_KG_M3)
        {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.finalize()
    }

    /// Prepare the source's rational dispersion terms once for a vacuum
    /// wavelength. The wavelength remains the same for every sampled state.
    pub fn prepare(wavelength_nm: f64) -> Result<PreparedWaterRefraction, WaterRefractionError> {
        if !wavelength_nm.is_finite()
            || !(WAVELENGTH_RANGE_NM[0]..=WAVELENGTH_RANGE_NM[1]).contains(&wavelength_nm)
        {
            return Err(WaterRefractionError::WavelengthOutsideRange);
        }
        let lambda = wavelength_nm / LAMBDA_STAR_NM;
        let square = lambda * lambda;
        let uv = square - LAMBDA_UV * LAMBDA_UV;
        let ir = square - LAMBDA_IR * LAMBDA_IR;
        let constant = A[0] + A[4] / square + A[5] / uv + A[6] / ir;
        let temperature = (A[2] + A[3] * square) / T_STAR_K;
        let wavelength_constant = -2.0 * lambda / LAMBDA_STAR_NM
            * (A[4] / (square * square) + A[5] / (uv * uv) + A[6] / (ir * ir));
        let wavelength_temperature = 2.0 * A[3] * lambda / (LAMBDA_STAR_NM * T_STAR_K);
        Ok(PreparedWaterRefraction {
            wavelength_nm,
            constant,
            temperature,
            wavelength_constant,
            wavelength_temperature,
        })
    }
}

/// Original wavelength coefficients, reusable across a spatial or transient
/// ordinary-water state without repeated divisions by dispersion denominators.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedWaterRefraction {
    wavelength_nm: f64,
    constant: f64,
    temperature: f64,
    wavelength_constant: f64,
    wavelength_temperature: f64,
}

/// Absolute phase/group indices and exact partial derivatives of this same law.
/// Density is mass density in kg/m3; wavelength is in vacuum nm. Temperature
/// partials hold density fixed, and density partials hold temperature fixed.
/// A physical isobaric path must additionally differentiate its original EOS.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterRefractiveSample {
    pub phase_index: f64,
    /// Stable n-1, including dilute states where rounded n equals one.
    pub phase_index_minus_one: f64,
    /// `n - lambda_vac * dn/dlambda_vac`, from analytic source dispersion.
    pub group_index: f64,
    pub wavelength_derivative_per_nm: f64,
    pub temperature_derivative_per_k: f64,
    pub density_derivative_m3_per_kg: f64,
    pub temperature_second_derivative_per_k2: f64,
    pub density_second_derivative_m6_per_kg2: f64,
    pub temperature_density_derivative_m3_per_kg_k: f64,
}

impl PreparedWaterRefraction {
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }

    /// Evaluate original equation (1) and its state jet without allocation.
    /// The state owner is responsible for homogeneous ordinary-water identity
    /// and phase validity; a numerical point alone supplies neither.
    pub fn sample(
        &self,
        temperature_k: f64,
        density_kg_m3: f64,
    ) -> Result<WaterRefractiveSample, WaterRefractionError> {
        if !temperature_k.is_finite()
            || !(TEMPERATURE_RANGE_K[0]..=TEMPERATURE_RANGE_K[1]).contains(&temperature_k)
        {
            return Err(WaterRefractionError::TemperatureOutsideRange);
        }
        if !density_kg_m3.is_finite()
            || !(DENSITY_RANGE_KG_M3[0]..=DENSITY_RANGE_KG_M3[1]).contains(&density_kg_m3)
        {
            return Err(WaterRefractionError::DensityOutsideRange);
        }
        let r = density_kg_m3 / RHO_STAR_KG_M3;
        let base = self.temperature.mul_add(temperature_k, self.constant);
        // L=(n²-1)/(n²+2)=r*(base+a1*r+a7*r²), the released Lorentz-Lorenz law.
        let lorentz = r * (A[7].mul_add(r, A[1]).mul_add(r, base));
        let denominator = 1.0 - lorentz;
        let n2_minus_one = 3.0 * lorentz / denominator;
        let n = det::sqrt(1.0 + n2_minus_one);
        let refractivity = n2_minus_one / (n + 1.0);
        // Differentiate n(L). Retaining the analytic curvature is essential:
        // affine T/rho fields can still have a nonzero optical Hessian.
        let first = 1.5 / (n * denominator * denominator);
        let second = first * (2.0 / denominator - first / n);
        let l_t = r * self.temperature;
        let l_r = (3.0 * A[7] * r * r + 2.0 * A[1] * r + base) / RHO_STAR_KG_M3;
        let l_rr = (2.0 * A[1] + 6.0 * A[7] * r) / (RHO_STAR_KG_M3 * RHO_STAR_KG_M3);
        let l_tr = self.temperature / RHO_STAR_KG_M3;
        let l_lambda = r * self
            .wavelength_temperature
            .mul_add(temperature_k, self.wavelength_constant);
        let wavelength_derivative_per_nm = first * l_lambda;
        let sample = WaterRefractiveSample {
            phase_index: 1.0 + refractivity,
            phase_index_minus_one: refractivity,
            group_index: (-self.wavelength_nm)
                .mul_add(wavelength_derivative_per_nm, 1.0 + refractivity),
            wavelength_derivative_per_nm,
            temperature_derivative_per_k: first * l_t,
            density_derivative_m3_per_kg: first * l_r,
            temperature_second_derivative_per_k2: second * l_t * l_t,
            density_second_derivative_m6_per_kg2: second.mul_add(l_r * l_r, first * l_rr),
            temperature_density_derivative_m3_per_kg_k: second.mul_add(l_t * l_r, first * l_tr),
        };
        if denominator <= 0.0
            || sample.phase_index < 1.0
            || sample.group_index <= 0.0
            || ![
                sample.phase_index,
                sample.phase_index_minus_one,
                sample.group_index,
                sample.wavelength_derivative_per_nm,
                sample.temperature_derivative_per_k,
                sample.density_derivative_m3_per_kg,
                sample.temperature_second_derivative_per_k2,
                sample.density_second_derivative_m6_per_kg2,
                sample.temperature_density_derivative_m3_per_kg_k,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(WaterRefractionError::NumericalRange);
        }
        Ok(sample)
    }

    /// Observe one original admitted IF97 Gibbs fluid state, retaining its
    /// actual temperature and density. `If97State` represents one homogeneous
    /// liquid/vapor/supercritical branch, never ice or a bulk two-phase mixture.
    ///
    /// This explicitly combines the R9-97 T/rho optical law with IF97 density.
    /// It does not replace that density with IAPWS-95 or claim the IAPWS-95
    /// pressure-based verification precision from R9-97 Table 3. The optical
    /// temperature/density domain remains checked after the state extraction.
    pub fn sample_if97_state(
        &self,
        state: &If97State,
    ) -> Result<WaterRefractiveSample, WaterRefractionError> {
        self.sample(state.temperature().value(), state.density().value())
    }
}

#[cfg(test)]
mod tests;
