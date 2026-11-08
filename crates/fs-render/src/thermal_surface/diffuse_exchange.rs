//! Reciprocal finite-band heat exchange from the original diffuse emissivity.
//!
//! One physical emitting side exchanges with an isotropic blackbody bath:
//! `P_out = area * pi * integral epsilon(lambda) * [B(T)-B(T_bath)] dlambda`.
//! The same piecewise-linear source used for observed light supplies absorption
//! and emission. No gray fit, RGB conversion or out-of-band continuation enters.

use super::coating_exchange::{Sum, scaled_factors};
use super::{SpectralEmissivity, ThermalSurfaceError, planck_radiative_exchange_w_and_w_k};
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};
use fs_math::det;
use std::sync::Arc;

/// Positive four-point Gauss rule in log wavelength, split at every source knot.
/// The original emissivity still interpolates linearly in wavelength. Subdivide
/// and compare results to assess integration error; this is not an error bound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DiffuseExchangeQuadrature {
    pub subdivisions_per_interval: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DiffuseExchangeBudget {
    pub max_spectral_nodes: usize,
}

/// Signed in-band watts, the frozen-emissivity secant, and unmodeled coverage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffuseExchangeRates {
    /// Positive means heat leaves the body; negative means bath heating.
    pub outward_power_w: f64,
    /// P_out/(T-T_bath), or the analytic tangent at equality; vacuum uses P/T.
    /// This excludes the derivative of a temperature/enthalpy-dependent law.
    pub secant_w_k: f64,
    /// Blackbody upper bound on the magnitude of omitted-wavelength exchange.
    /// This never supplies extra heat and is independent of quadrature error.
    pub omitted_exchange_upper_bound_w: f64,
}

#[derive(Clone, Copy, Debug)]
struct Node {
    wavelength_nm: f64,
    wavelength_weight_nm: f64,
}

/// Prepared finite-band response with reusable, source-knot-respecting nodes.
///
/// Emissivity is the original accepted spectral law. State-dependent consumers
/// resolve that same law and use [`Self::with_emissivity`] to share the grid.
/// Temperature evaluation allocates nothing and uses the existing stable SI
/// Planck difference/analytic-secant owner, including near thermal equilibrium.
#[derive(Clone, Debug)]
pub struct DiffuseExchangeResponse {
    nodes: Arc<[Node]>,
    emissivity: SpectralEmissivity,
    area_m2: f64,
    quadrature: DiffuseExchangeQuadrature,
    identity: ContentHash,
}

impl DiffuseExchangeResponse {
    pub fn from_emissivity(
        emissivity: SpectralEmissivity,
        area_m2: f64,
        quadrature: DiffuseExchangeQuadrature,
        budget: DiffuseExchangeBudget,
        mut checkpoint: impl FnMut() -> Result<(), DiffuseExchangeError>,
    ) -> Result<Self, DiffuseExchangeError> {
        checkpoint()?;
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(DiffuseExchangeError::Area);
        }
        let count = (emissivity.knots().len() - 1)
            .checked_mul(quadrature.subdivisions_per_interval)
            .and_then(|count| count.checked_mul(4))
            .filter(|&count| count > 0 && count <= budget.max_spectral_nodes && count <= 65_536)
            .ok_or(DiffuseExchangeError::Budget)?;
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(count)
            .map_err(|_| DiffuseExchangeError::Budget)?;
        const GAUSS: [(f64, f64); 4] = [
            (-0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
            (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
            (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
            (0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
        ];
        for pair in emissivity.knots().windows(2) {
            let log_left = det::ln(pair[0].0);
            let log_span = det::ln(pair[1].0) - log_left;
            let width = log_span / quadrature.subdivisions_per_interval as f64;
            if !width.is_finite() || width <= 0.0 {
                return Err(DiffuseExchangeError::Quadrature);
            }
            for index in 0..quadrature.subdivisions_per_interval {
                checkpoint()?;
                for (abscissa, weight) in GAUSS {
                    let fraction = index as f64 + 0.5 * (1.0 + abscissa);
                    let wavelength_nm = det::exp(fraction.mul_add(width, log_left));
                    let wavelength_weight_nm = wavelength_nm * (0.5 * width * weight);
                    if !wavelength_nm.is_finite()
                        || wavelength_nm <= pair[0].0
                        || wavelength_nm >= pair[1].0
                        || !wavelength_weight_nm.is_finite()
                        || wavelength_weight_nm <= 0.0
                    {
                        return Err(DiffuseExchangeError::Quadrature);
                    }
                    nodes.push(Node {
                        wavelength_nm,
                        wavelength_weight_nm,
                    });
                }
            }
        }
        let mut response = Self {
            nodes: nodes.into(),
            emissivity,
            area_m2,
            quadrature,
            identity: ContentHash([0; 32]),
        };
        response.identity = response.compute_identity();
        checkpoint()?;
        Ok(response)
    }

    /// Reuse this exact quadrature for a new accepted state of the SAME knot
    /// grid. This changes no existing response and never extrapolates support.
    pub fn with_emissivity(
        &self,
        emissivity: SpectralEmissivity,
    ) -> Result<Self, DiffuseExchangeError> {
        if emissivity.knots().len() != self.emissivity.knots().len()
            || emissivity
                .knots()
                .iter()
                .zip(self.emissivity.knots())
                .any(|(left, right)| left.0.to_bits() != right.0.to_bits())
        {
            return Err(DiffuseExchangeError::SourceGrid);
        }
        let mut next = self.clone();
        next.emissivity = emissivity;
        next.identity = next.compute_identity();
        Ok(next)
    }

    /// Only the declared finite band contributes to the applied heat. Unknown
    /// wavelengths are bounded using passivity (0 <= epsilon <= 1), never filled
    /// with guessed infrared data. The omitted bound is exactly zero at equal
    /// temperatures and otherwise uses the hotter blackbody's omitted emission.
    pub fn evaluate(
        &self,
        temperature_k: f64,
        bath_temperature_k: Option<f64>,
        mut checkpoint: impl FnMut() -> Result<(), DiffuseExchangeError>,
    ) -> Result<DiffuseExchangeRates, DiffuseExchangeError> {
        checkpoint()?;
        planck_radiative_exchange_w_and_w_k(
            temperature_k,
            bath_temperature_k,
            self.wavelength_bounds_nm().0,
            0.0,
            0.0,
        )?;
        let mut power = Sum::default();
        let mut secant = Sum::default();
        for node in self.nodes.iter() {
            checkpoint()?;
            let epsilon = self.emissivity.at(node.wavelength_nm)?;
            let (coefficient, weight) =
                scaled_factors(self.area_m2, PI * epsilon, node.wavelength_weight_nm)?;
            let (heat, conductance) = planck_radiative_exchange_w_and_w_k(
                temperature_k,
                bath_temperature_k,
                node.wavelength_nm,
                coefficient,
                weight,
            )?;
            power.add(heat)?;
            secant.add(conductance)?;
        }
        let outward_power_w = power.value()?;
        let secant_w_k = secant.value()?;
        if outward_power_w != 0.0 && secant_w_k == 0.0 {
            return Err(ThermalSurfaceError::Unrepresentable.into());
        }
        let omitted_exchange_upper_bound_w = if bath_temperature_k == Some(temperature_k) {
            0.0
        } else {
            super::planck_band::omitted_blackbody_power_upper_bound_w(
                temperature_k.max(bath_temperature_k.unwrap_or(temperature_k)),
                self.wavelength_bounds_nm(),
                self.area_m2,
            )
            .ok_or(ThermalSurfaceError::Unrepresentable)?
        };
        checkpoint()?;
        Ok(DiffuseExchangeRates {
            outward_power_w,
            secant_w_k,
            omitted_exchange_upper_bound_w,
        })
    }

    pub const fn emissivity(&self) -> SpectralEmissivity {
        self.emissivity
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn quadrature(&self) -> DiffuseExchangeQuadrature {
        self.quadrature
    }
    pub fn spectral_node_count(&self) -> usize {
        self.nodes.len()
    }
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.emissivity.wavelength_bounds_nm()
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    fn compute_identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("org.frankensim.render.diffuse-exchange-response.v1");
        hash.update(self.emissivity.identity().as_bytes());
        hash.update(&self.area_m2.to_bits().to_le_bytes());
        hash.update(&(self.quadrature.subdivisions_per_interval as u64).to_le_bytes());
        hash.finalize()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum DiffuseExchangeError {
    Area,
    Budget,
    Quadrature,
    SourceGrid,
    Cancelled,
    Surface(ThermalSurfaceError),
}
impl From<ThermalSurfaceError> for DiffuseExchangeError {
    fn from(value: ThermalSurfaceError) -> Self {
        Self::Surface(value)
    }
}
impl core::fmt::Display for DiffuseExchangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "diffuse spectral exchange refused: {self:?}")
    }
}
impl std::error::Error for DiffuseExchangeError {}

#[cfg(test)]
mod tests;
