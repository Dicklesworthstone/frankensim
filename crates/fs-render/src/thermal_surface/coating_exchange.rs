//! Finite-band exchange between an isothermal coherent stack and a black bath.
//!
//! The original positive reciprocal absorption supplies the same directional
//! emissivity as [`CoatingThermalEmission`]. For one vacuum-facing port,
//! `P_out = area * integral_lambda integral_hemisphere
//! A_unpolarized(lambda,mu) * [B_lambda(T)-B_lambda(T_bath)] * mu dOmega dlambda`.
//! Axisymmetry supplies `dOmega = 2*pi*dmu`; Planck already includes both
//! polarizations, so `A_unpolarized=(A_s+A_p)/2`. No extra Lambertian pi or
//! photon-index factor enters. Optical constants remain frozen in each response.
//!
//! A common temperature for all films and the opaque substrate makes internal
//! radiative exchange cancel. This owner therefore does not approximate full
//! cooling of a nonisothermal stack by summing only its outward emission. A
//! transmitting response covers just its selected vacuum port; any second
//! exposed port is a separate physical exchange. Blackbody surroundings are
//! isotropic over that entire port's outgoing hemisphere.

use super::{
    CoatingEmissionError, CoatingThermalEmission, ThermalSurfaceError,
    planck_radiative_exchange_w_and_w_k,
};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};
use fs_math::det;
use std::sync::Arc;

/// Explicit composite-midpoint rule over vacuum wavelength and outward cosine.
///
/// Wavelength bins cover exactly 380--780 nm; cosine bins cover `(0,1)` with
/// uniform widths `1/cosine_bins`. The owner supplies the azimuthal factor 2*pi
/// and the projected-area cosine. Nodes avoid the singular exactly tangent mode.
/// These are approximation controls, not an error certificate: refine both
/// counts for sharp coherent spectral or angular features.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoatingExchangeQuadrature {
    pub wavelength_bins: usize,
    pub cosine_bins: usize,
}

/// Bounds on retained spectral storage and original coherent optical solves.
/// A single solve retains the optical owner's existing 64-film work bound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoatingExchangeBudget {
    pub max_spectral_nodes: usize,
    pub max_optical_evaluations: usize,
}

/// Signed outward watts and the positive frozen-optics thermal conductance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoatingExchangeRates {
    /// Positive means heat leaves the stack; negative means the bath heats it.
    pub outward_power_w: f64,
    /// P_out/(T-T_bath), W/K. Equal temperatures give the analytic tangent.
    /// For a vacuum bath this is P_out/T, never a secant through a rounded heat.
    pub secant_w_k: f64,
}

#[derive(Clone, Copy, Debug)]
struct SpectralNode {
    wavelength_nm: f64,
    // Numerical factorization of area * angular extraction * wavelength width.
    // These two positive scales enter the original log-scaled Planck owner
    // separately, even when their product would underflow or overflow.
    coefficient: f64,
    weight: f64,
}

/// Reusable hemispherical exchange response for a declared isothermal stack.
///
/// Preparation integrates the original positive optical work over angle once
/// at each wavelength. A temperature query then evaluates only the native
/// stable Planck differences and secants, with no optical solves or allocation.
/// The response retains fixed optics, area and quadrature, not the preparation
/// temperature. A temperature/enthalpy-dependent law must be resolved again
/// before reuse; [`Self::optical_identity`] identifies the retained optical data.
///
/// Only the original 380--780 nm band contributes. This is neither bolometric
/// cooling nor an infrared continuation of a visible optical table. It is a
/// physical power response; a native thermal owner applies its signed heat and
/// combines its secant with the native caloric slope and other stability limits.
#[derive(Clone, Debug)]
pub struct CoatingExchangeResponse {
    nodes: Arc<[SpectralNode]>,
    area_m2: f64,
    quadrature: CoatingExchangeQuadrature,
    optical_evaluations: usize,
    optical_identity: ContentHash,
    identity: ContentHash,
}

impl CoatingExchangeResponse {
    /// Prepare exactly one physical vacuum-facing port from the existing
    /// immutable emitter. Every film and opaque-substrate temperature must be
    /// exactly equal. An empty transmitting interface has no emitting material
    /// and is admitted as a zero exchange response without a fictitious body.
    ///
    /// The area is a positive physical emitting-side area, in square metres.
    /// Storage/work admission occurs before allocation or any optical solve.
    /// Checkpoints run before each bounded solve and before publication. Failed
    /// preparation publishes no partial response and never changes the emitter.
    pub fn from_emitter(
        emitter: &CoatingThermalEmission,
        area_m2: f64,
        quadrature: CoatingExchangeQuadrature,
        budget: CoatingExchangeBudget,
        mut checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<Self, CoatingExchangeError> {
        let optical_evaluations = admit(quadrature, budget, area_m2)?;
        let common_temperature = emitter
            .substrate_temperature_k()
            .or_else(|| emitter.film_temperature_k(0));
        for index in 0..emitter.film_count() {
            if emitter.film_temperature_k(index) != common_temperature {
                return Err(CoatingExchangeError::Nonisothermal);
            }
        }
        checkpoint()?;
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(quadrature.wavelength_bins)
            .map_err(|_| CoatingExchangeError::Allocation)?;
        let wavelength_width = (LAMBDA_MAX - LAMBDA_MIN) / quadrature.wavelength_bins as f64;
        let cosine_count = quadrature.cosine_bins as f64;
        let mut hash = DomainHasher::new("org.frankensim.render.coating-exchange-response.v1");
        hash.update(emitter.optical_identity().as_bytes());
        hash.update(&area_m2.to_bits().to_le_bytes());
        for count in [quadrature.wavelength_bins, quadrature.cosine_bins] {
            hash.update(&(count as u64).to_le_bytes());
        }
        for wavelength_index in 0..quadrature.wavelength_bins {
            let wavelength_nm = LAMBDA_MIN + (wavelength_index as f64 + 0.5) * wavelength_width;
            let mut angular = Sum::default();
            for cosine_index in 0..quadrature.cosine_bins {
                checkpoint()?;
                let cosine = (cosine_index as f64 + 0.5) / cosine_count;
                let [s, p] = emitter.reciprocal_absorbed_fraction(wavelength_nm, cosine)?;
                // Keep the polarization sum before halving: dividing tiny
                // positive absorptances separately can erase a real source.
                // 2*pi and 1/2 combine to pi below, after the positive sum.
                angular.add(positive_product(cosine, s + p)?)?;
            }
            let angular_value = angular.value()?;
            let mean_projected = angular_value / cosine_count;
            if angular_value > 0.0 && mean_projected == 0.0 {
                return Err(ThermalSurfaceError::Unrepresentable.into());
            }
            let extraction_sr = PI * mean_projected;
            let (coefficient, weight) = scaled_factors(area_m2, extraction_sr, wavelength_width)?;
            for value in [wavelength_nm, extraction_sr, coefficient, weight] {
                hash.update(&value.to_bits().to_le_bytes());
            }
            nodes.push(SpectralNode {
                wavelength_nm,
                coefficient,
                weight,
            });
        }
        checkpoint()?;
        Ok(Self {
            nodes: nodes.into(),
            area_m2,
            quadrature,
            optical_evaluations,
            optical_identity: emitter.optical_identity(),
            identity: hash.finalize(),
        })
    }

    /// Net outward power and positive temperature secant at the supplied
    /// common stack temperature and isotropic blackbody bath temperature.
    /// `None` means an explicitly empty incoming vacuum bath (`B_bath=0`).
    /// Every supplied temperature must be finite and strictly positive, even
    /// for a dark response. Equality returns exactly zero power and the native
    /// analytic tangent. Swapping two bath/body temperatures negates the heat
    /// and preserves the secant for these same frozen optical coefficients.
    ///
    /// The secant is a frozen-coefficient exchange slope, not the derivative
    /// of a temperature-dependent optical law and not a time-accuracy bound.
    /// Evaluate resolved state laws at accepted thermal steps and refine time.
    pub fn evaluate(
        &self,
        temperature_k: f64,
        bath_temperature_k: Option<f64>,
        mut checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<CoatingExchangeRates, CoatingExchangeError> {
        checkpoint()?;
        // Native source admission remains binding even if all optical weights
        // are exactly zero. This call performs no gross radiance conversion.
        planck_radiative_exchange_w_and_w_k(
            temperature_k,
            bath_temperature_k,
            LAMBDA_MIN,
            0.0,
            0.0,
        )?;
        let mut power = Sum::default();
        let mut secant = Sum::default();
        for node in self.nodes.iter() {
            checkpoint()?;
            let (heat, conductance) = planck_radiative_exchange_w_and_w_k(
                temperature_k,
                bath_temperature_k,
                node.wavelength_nm,
                node.coefficient,
                node.weight,
            )?;
            power.add(heat)?;
            secant.add(conductance)?;
        }
        let outward_power_w = power.value()?;
        let secant_w_k = secant.value()?;
        if outward_power_w != 0.0 && secant_w_k == 0.0 {
            return Err(ThermalSurfaceError::Unrepresentable.into());
        }
        checkpoint()?;
        Ok(CoatingExchangeRates {
            outward_power_w,
            secant_w_k,
        })
    }

    /// Positive physical area of the selected observed port, m².
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Explicit spectral and angular approximation controls.
    #[must_use]
    pub const fn quadrature(&self) -> CoatingExchangeQuadrature {
        self.quadrature
    }

    /// Number of original coherent solves performed during preparation.
    /// Queries of [`Self::evaluate`] perform zero optical solves.
    #[must_use]
    pub const fn optical_evaluations(&self) -> usize {
        self.optical_evaluations
    }

    /// Fixed optical laws, sources, original film geometry and selected side.
    /// Temperature is absent; state-dependent laws must be resolved before
    /// comparing this identity to a newer emitter's optical identity.
    #[must_use]
    pub const fn optical_identity(&self) -> ContentHash {
        self.optical_identity
    }

    /// Retained optical response, area, quadrature and scaled spectral weights.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Explicit admission, cancellation and original optical/source refusals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoatingExchangeError {
    /// Both midpoint counts must be positive and have a representable product.
    Quadrature,
    /// The declared spectral-storage or optical-solve budget is insufficient.
    Budget,
    /// Fallible storage reservation failed before a response was published.
    Allocation,
    /// A physical emitting-side area must be finite and strictly positive.
    Area,
    /// Independent film/substrate temperatures need an internal-exchange model.
    Nonisothermal,
    /// A caller checkpoint cancelled preparation or spectral evaluation.
    Cancelled,
    /// Original reciprocal optical or emitted-state admission.
    Emission(CoatingEmissionError),
    /// Original stable Planck or representability refusal.
    Thermal(ThermalSurfaceError),
}
impl From<CoatingEmissionError> for CoatingExchangeError {
    fn from(error: CoatingEmissionError) -> Self {
        Self::Emission(error)
    }
}
impl From<ThermalSurfaceError> for CoatingExchangeError {
    fn from(error: ThermalSurfaceError) -> Self {
        Self::Thermal(error)
    }
}
impl core::fmt::Display for CoatingExchangeError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "coherent coating exchange refused: {self:?}")
    }
}
impl core::error::Error for CoatingExchangeError {}

fn admit(
    quadrature: CoatingExchangeQuadrature,
    budget: CoatingExchangeBudget,
    area_m2: f64,
) -> Result<usize, CoatingExchangeError> {
    if !area_m2.is_finite() || area_m2 <= 0.0 {
        return Err(CoatingExchangeError::Area);
    }
    if quadrature.wavelength_bins == 0 || quadrature.cosine_bins == 0 {
        return Err(CoatingExchangeError::Quadrature);
    }
    let evaluations = quadrature
        .wavelength_bins
        .checked_mul(quadrature.cosine_bins)
        .ok_or(CoatingExchangeError::Quadrature)?;
    if quadrature.wavelength_bins > budget.max_spectral_nodes
        || evaluations > budget.max_optical_evaluations
    {
        return Err(CoatingExchangeError::Budget);
    }
    Ok(evaluations)
}

fn positive_product(first: f64, second: f64) -> Result<f64, ThermalSurfaceError> {
    let value = first * second;
    if !value.is_finite() || value < 0.0 || (first > 0.0 && second > 0.0 && value == 0.0) {
        Err(ThermalSurfaceError::Unrepresentable)
    } else {
        Ok(value)
    }
}

fn scaled_factors(
    area: f64,
    extraction: f64,
    wavelength_width: f64,
) -> Result<(f64, f64), ThermalSurfaceError> {
    if !extraction.is_finite() || extraction < 0.0 {
        return Err(ThermalSurfaceError::Unrepresentable);
    }
    if extraction == 0.0 {
        return Ok((0.0, wavelength_width));
    }
    let coefficient = area * extraction;
    if coefficient.is_normal() {
        return Ok((coefficient, wavelength_width));
    }
    // Keep the area separate when the first product overflows or enters the
    // subnormal range. The Planck owner applies both logarithmic scales before
    // rounding any power or conductance.
    let weight = extraction * wavelength_width;
    if weight.is_normal() {
        return Ok((area, weight));
    }
    // All three positive factors may individually be representable while both
    // straightforward products underflow. A balanced numerical factorization
    // preserves the same product without creating a zero exchange coefficient.
    let log_scale = det::ln(area) + det::ln(extraction) + det::ln(wavelength_width);
    let factor = det::exp(0.5 * log_scale);
    if !factor.is_normal() {
        return Err(ThermalSurfaceError::Unrepresentable);
    }
    Ok((factor, factor))
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), ThermalSurfaceError> {
        let next = self.value + value;
        if !value.is_finite() || !next.is_finite() {
            return Err(ThermalSurfaceError::Unrepresentable);
        }
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
        Ok(())
    }
    fn value(&self) -> Result<f64, ThermalSurfaceError> {
        let total = self.value + self.correction;
        if total.is_finite() {
            Ok(total)
        } else {
            Err(ThermalSurfaceError::Unrepresentable)
        }
    }
}

#[cfg(test)]
mod tests;
