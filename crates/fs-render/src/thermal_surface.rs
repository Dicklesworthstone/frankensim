//! Physical spectral surface emission in local thermal equilibrium.
//!
//! The gray law is an explicit constitutive approximation: one emissivity applies at
//! every wavelength and direction, with diffuse reflectance `1 - emissivity`.
//! A hemispherical-total measurement alone does not establish that assumption.
//! The material-state adapter binds this declared approximation to the SAME
//! source sample and temperature used by the conduction radiation owner.
//!
//! Spectral radiance is W/(m² sr nm). No RGB lift, exposure or display mapping
//! enters the physical spectrum. The visible tracer integrates its own band;
//! cooling retains the complete Stefan--Boltzmann spectrum in its own owner.
//! The gray and spectral diffuse surfaces are two-sided: each face emits into
//! its outward hemisphere. Exitance is per emitting side, so a sheet exposing
//! both sides has twice its one-face area for thermal accounting. A coherent
//! stack instead selects one physical observation side explicitly.
//! This spectrum is for emission into vacuum (the tracer's exterior air
//! approximation). A thermal hit within an active dielectric medium refuses;
//! no refractive-index correction is invented for the original cooling model.
//!
//! [`SpectralDiffuseThermalSurface`] additionally accepts a declared finite-band
//! emissivity curve and uses it for both reciprocal reflection and emission.
//! That finite-band law makes no claim about bolometric cooling or out-of-band
//! continuation; its thermal owner remains responsible for total energy.
//!
//! [`RoughDiffuseThermalSurface`] retains explicitly declared EON facet albedo
//! and roughness. Its directional emissivity is the absorptance of that same
//! reciprocal reflection law, evaluated at the actual outgoing direction.
//! Macroscopic emissivity measurements are not converted into facet data.
//!
//! [`CoatingThermalEmission`] resolves reciprocal directional emission from
//! passive coherent films at separate temperatures. Its selected outgoing port
//! is vacuum, its spectrum remains inside the original optical support, and
//! its film order and polarization are explicit rather than diffuse or gray.

#[cfg(feature = "material-state")]
mod entropy;

mod planck;
mod planck_band;
pub(crate) use planck::{
    planck_radiance_difference_w_m2_sr_nm, planck_radiance_temperature_derivative_w_m2_sr_nm_k,
    planck_radiance_temperature_secant_w_m2_sr_nm_k, planck_radiative_exchange_w_and_w_k,
};

mod spectral;
pub use spectral::{DiffuseThermalSurface, SpectralDiffuseThermalSurface, SpectralEmissivity};

mod rough;
pub use rough::RoughDiffuseThermalSurface;

#[cfg(feature = "photometry")]
mod observables;
#[cfg(feature = "photometry")]
pub use observables::{
    HeldThermalEmissionInterval, ThermalEmissionError, ThermalEmissionPatch, ThermalVisibleEmission,
    ThermalVisibleExposure, ThermalVisibleObservables,
};

/// Original thermal spectra observed from moving inertial vacuum frames.
pub mod doppler;

/// Source-consistent finite-band exchange between diffuse thermal surfaces.
pub mod diffuse_exchange;

mod coating_emission;
pub use coating_emission::{CoatingEmissionError, CoatingEmissionSample, CoatingThermalEmission};

/// Full polarized Smith reflection with the same local Kirchhoff thermal source.
pub mod smith;
pub use smith::{SmithCoatingThermalEmission, SmithThermalError, SmithThermalSample};

#[cfg(feature = "crystal-emission")]
mod crystal_emission;
#[cfg(feature = "crystal-emission")]
pub use crystal_emission::{
    CrystalEmissionError, PreparedUniaxialSlabEmission, UniaxialSlabThermalEmission,
};

#[cfg(feature = "crystal-emission")]
mod crystal_stack_emission;
#[cfg(feature = "crystal-emission")]
pub use crystal_stack_emission::{
    CrystalStackEmissionError, PreparedUniaxialStackEmission, UniaxialStackThermalEmission,
};

mod coating_exchange;
pub use coating_exchange::{
    CoatingExchangeBudget, CoatingExchangeError, CoatingExchangeQuadrature, CoatingExchangeRates,
    CoatingExchangeResponse,
};

mod coating_gap;
pub use coating_gap::CoatingGapResponse;

mod coating_external;
pub use coating_external::{
    CoatingExternalBathBudget, CoatingExternalBathRates, CoatingExternalBathResponse,
};

mod coating_transmitting_pair;
pub use coating_transmitting_pair::{
    CoatingTransmittingPairBudget, CoatingTransmittingPairRates, CoatingTransmittingPairResponse,
};

mod coating_radiative_stack;
pub use coating_radiative_stack::{
    CoatingRadiativePairRate, CoatingRadiativeStackBudget, CoatingRadiativeStackRates,
    CoatingRadiativeStackResponse, MAX_RADIATIVE_STACK_BODIES,
};

/// Physical incident surface absorption on the original native thermal clock.
#[cfg(feature = "optical-heating")]
pub mod heating;

use fs_blake3::{ContentHash, DomainHasher};
use fs_math::det;

#[cfg(test)]
use planck_band::PLANCK_J_S;
use planck_band::{BOLTZMANN_J_K, C1_NM, C2_NM_K, LIGHT_M_S};
const SIGMA_W_M2_K4: f64 = 5.670_374_419e-8;

/// A gray/diffuse equilibrium model with a retained physical-source identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrayDiffuseThermalSurface {
    temperature_k: f64,
    emissivity: f64,
    source_identity: ContentHash,
}

impl GrayDiffuseThermalSurface {
    /// Declare a gray/diffuse approximation explicitly. This constructor does
    /// not turn the supplied identity into measurement authority. Material-card
    /// consumers should use `MaterialThermalBinding::try_gray_diffuse`.
    pub fn try_new(
        temperature_k: f64,
        emissivity: f64,
        source_identity: ContentHash,
    ) -> Result<Self, ThermalSurfaceError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(ThermalSurfaceError::Temperature);
        }
        if !emissivity.is_finite() || !(0.0..=1.0).contains(&emissivity) {
            return Err(ThermalSurfaceError::Emissivity);
        }
        if source_identity == ContentHash([0; 32]) {
            return Err(ThermalSurfaceError::SourceIdentity);
        }
        let surface = Self {
            temperature_k,
            emissivity,
            source_identity,
        };
        if !surface.bolometric_exitance_w_m2().is_finite() {
            return Err(ThermalSurfaceError::Unrepresentable);
        }
        Ok(surface)
    }

    /// Accepted thermodynamic temperature, K.
    #[must_use]
    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }

    /// The declared wavelength- and direction-independent emissivity.
    #[must_use]
    pub const fn emissivity(self) -> f64 {
        self.emissivity
    }

    /// Opaque reciprocal diffuse reflectance at every wavelength.
    #[must_use]
    pub fn reflectance(self) -> f64 {
        1.0 - self.emissivity
    }

    /// Hemispherical emitted power per area of one emitting side, W/m².
    /// This is an observation, never an additional debit to a thermal state.
    #[must_use]
    pub fn bolometric_exitance_w_m2(self) -> f64 {
        let t2 = self.temperature_k * self.temperature_k;
        self.emissivity * SIGMA_W_M2_K4 * t2 * t2
    }

    /// Emitted spectral radiance, W/(m² sr nm), at positive wavelength in nm.
    pub fn radiance_w_m2_sr_nm(self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
            return Err(ThermalSurfaceError::Wavelength);
        }
        if self.emissivity == 0.0 {
            return Ok(0.0);
        }
        let value =
            self.emissivity * planck_radiance_w_m2_sr_nm(self.temperature_k, wavelength_nm)?;
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(ThermalSurfaceError::Unrepresentable)
        }
    }

    /// Source, temperature and explicit gray/diffuse-model identity.
    #[must_use]
    pub fn identity(self) -> ContentHash {
        let mut h = DomainHasher::new("org.frankensim.render.gray-diffuse-thermal.v1");
        h.update(self.source_identity.as_bytes());
        h.update(&self.temperature_k.to_bits().to_le_bytes());
        h.update(&self.emissivity.to_bits().to_le_bytes());
        h.finalize()
    }
}

/// Shared SI blackbody source law for surfaces and absorbing/emitting media.
pub(crate) fn planck_radiance_w_m2_sr_nm(
    temperature_k: f64,
    wavelength_nm: f64,
) -> Result<f64, ThermalSurfaceError> {
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(ThermalSurfaceError::Temperature);
    }
    if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
        return Err(ThermalSurfaceError::Wavelength);
    }
    let x = (C2_NM_K / wavelength_nm) / temperature_k;
    checked_planck_exp(planck_log_radiance(temperature_k, wavelength_nm, x))
}

// Shared algebra of the original SI source law. Retaining the caller's original
// exponent/order preserves the established surface and volume radiance values.
fn planck_log_radiance(temperature_k: f64, wavelength_nm: f64, x: f64) -> f64 {
    // Log evaluation avoids lambda^-5 overflowing before the Wien tail
    // vanishes. The small-x branch is the representable Rayleigh--Jeans
    // limit; x>=50 makes log(expm1(x)) indistinguishable from x in f64.
    if x == 0.0 {
        det::ln(2.0 * LIGHT_M_S * BOLTZMANN_J_K * 1.0e27) + det::ln(temperature_k)
            - 4.0 * det::ln(wavelength_nm)
    } else {
        let log_denominator = if x >= 50.0 { x } else { det::ln(det::expm1(x)) };
        det::ln(C1_NM) - 5.0 * det::ln(wavelength_nm) - log_denominator
    }
}

fn checked_planck_exp(log_value: f64) -> Result<f64, ThermalSurfaceError> {
    let value = det::exp(log_value);
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(ThermalSurfaceError::Unrepresentable)
    }
}

/// Refusal before a thermal appearance enters a scene.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermalSurfaceError {
    /// Temperature must be finite and strictly positive in K.
    Temperature,
    /// Emissivity must be finite and in [0,1].
    Emissivity,
    /// The emissivity table must contain 2–16 strictly increasing positive
    /// finite wavelength knots. The legacy visible constructor requires exact
    /// 380 and 780 nm endpoints.
    SpectralEmissivity,
    /// The physical-source identity must not be zero.
    SourceIdentity,
    /// Wavelength must be finite and inside the model's declared domain.
    Wavelength,
    /// An emission or absorption cosine must be finite and in [0,1].
    Direction,
    /// A directional rough source requires its actual outgoing cosine.
    DirectionRequired,
    /// A Planck exchange coefficient or quadrature weight was negative or nonfinite.
    ExchangeCoefficient,
    /// The requested physical radiance or power cannot be represented.
    Unrepresentable,
}

impl core::fmt::Display for ThermalSurfaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SpectralEmissivity => {
                f.write_str("invalid source-band spectral emissivity table")
            }
            _ => write!(f, "invalid gray/diffuse thermal surface: {self:?}"),
        }
    }
}
impl core::error::Error for ThermalSurfaceError {}

#[cfg(test)]
mod tests;
