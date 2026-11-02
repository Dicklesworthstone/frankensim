//! Layer-resolved exchange with one explicitly selected exterior vacuum bath.
//!
//! Each original coherent solve supplies the positive work absorbed in every
//! film under illumination from this port. Local Kirchhoff reciprocity gives
//! the same film's escaped thermal emission into the reversed vacuum mode:
//! `Q_i,out = area * integral pi*mu*(A_i,s+A_i,p)*[B(T_i)-B(T_bath)] dmu dlambda`.
//! See Greffet et al., *Light Emission by Nonequilibrium Bodies: Local Kirchhoff
//! Law*, Phys. Rev. X 8, 021008 (2018), equations (5)--(6),
//! <https://doi.org/10.1103/PhysRevX.8.021008>.
//!
//! This is EXTERIOR exchange, not a complete nonisothermal radiative model.
//! Radiation emitted by one layer and absorbed by another needs a separate
//! internal-transfer operator. Those transfers cancel in the total specimen
//! energy but can change individual layer temperatures. No internal coefficients
//! are inferred from exterior absorptances. Each layer is locally thermalized,
//! spatially uniform in temperature, passive, isotropic and reciprocal.
//!
//! A transmitting stack covers only its selected vacuum-facing port. A finite
//! vacuum/vacuum specimen can author both independent exterior baths; reverse
//! traversal retains authored film/temperature order. A material-side port is
//! not vacuum, and transmitted bath light is not film heating. A single absorbing
//! thermal subsystem, or a complete internally isothermal specimen, needs no
//! additional internal radiative redistribution to close its thermal account.

use super::coating_exchange::{
    CoatingExchangeBudget, CoatingExchangeError, CoatingExchangeQuadrature, CoatingExchangeRates,
    SpectralNode, Sum, admit, positive_product, scaled_factors,
};
use super::{CoatingThermalEmission, ThermalSurfaceError, planck_radiative_exchange_w_and_w_k};
use crate::conductor::MAX_COHERENT_STACK_LAYERS;
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};
use std::sync::Arc;

const MAX_CHANNELS: usize = MAX_COHERENT_STACK_LAYERS + 1;
const ZERO_RATE: CoatingExchangeRates = CoatingExchangeRates {
    outward_power_w: 0.0,
    secant_w_k: 0.0,
};

/// Bounds on retained layer/spectral slots and original coherent optical work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoatingExternalBathBudget {
    /// `(film_count + opaque_substrate_count) * wavelength_bins`, including
    /// explicitly lossless/zero-thickness films. Angles are not retained.
    pub max_channel_spectral_nodes: usize,
    /// One coherent solve per wavelength/cosine node, independent of film count.
    pub max_optical_evaluations: usize,
}

/// Exterior contribution to each original thermal channel, in authored order.
/// Positive watts leave that channel through this port. A positive secant is
/// relative to this bath at frozen optics, not an internal-transfer conductance
/// or a bound on derivatives of temperature/enthalpy-dependent optical laws.
#[derive(Clone, Debug, PartialEq)]
pub struct CoatingExternalBathRates {
    channels: [CoatingExchangeRates; MAX_CHANNELS],
    film_count: usize,
    has_substrate: bool,
    total_outward_power_w: f64,
}
impl CoatingExternalBathRates {
    #[must_use]
    pub fn film_rates(&self) -> &[CoatingExchangeRates] {
        &self.channels[..self.film_count]
    }
    #[must_use]
    pub const fn substrate_rate(&self) -> Option<CoatingExchangeRates> {
        if self.has_substrate {
            Some(self.channels[self.film_count])
        } else {
            None
        }
    }
    /// Signed sum of this port's layer/substrate contributions. Opposite layer
    /// signs can cancel; no artificial common-temperature secant is assigned.
    #[must_use]
    pub const fn total_outward_power_w(&self) -> f64 {
        self.total_outward_power_w
    }
}

/// Immutable finite-band angular integration of all exterior layer channels.
///
/// Preparation performs one original coherent solve per wavelength/cosine,
/// accumulates every layer simultaneously and retains only spectral weights.
/// Evaluation allocates nothing and performs no optical solve. Temperatures
/// do not enter this response's identity; state-dependent optical laws must be
/// resolved at the accepted state before comparing [`Self::optical_identity`].
#[derive(Clone, Debug)]
pub struct CoatingExternalBathResponse {
    nodes: Arc<[SpectralNode]>,
    film_count: usize,
    has_substrate: bool,
    area_m2: f64,
    quadrature: CoatingExchangeQuadrature,
    optical_evaluations: usize,
    wavelength_bounds_nm: (f64, f64),
    optical_identity: ContentHash,
    identity: ContentHash,
}
impl CoatingExternalBathResponse {
    /// Prepare one physical vacuum-facing port within its exact common optical
    /// source band. Every authored film retains a channel, even when its work is
    /// zero. Only an opaque stack has a substrate channel. A bare transmitting
    /// interface retains zero channel slots and exchanges zero heat.
    ///
    /// Storage/work limits are checked before allocation or optical solves.
    /// No isothermal restriction is imposed: each query supplies independent
    /// temperatures for this exterior contribution, not full internal cooling.
    pub fn from_emitter(
        emitter: &CoatingThermalEmission,
        area_m2: f64,
        quadrature: CoatingExchangeQuadrature,
        budget: CoatingExternalBathBudget,
        mut checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<Self, CoatingExchangeError> {
        let film_count = emitter.film_count();
        let has_substrate = emitter.substrate_temperature_k().is_some();
        let channel_count = film_count + usize::from(has_substrate);
        let slots = channel_count
            .checked_mul(quadrature.wavelength_bins)
            .ok_or(CoatingExchangeError::Budget)?;
        let optical_evaluations = admit(
            quadrature,
            CoatingExchangeBudget {
                max_spectral_nodes: quadrature.wavelength_bins,
                max_optical_evaluations: budget.max_optical_evaluations,
            },
            area_m2,
            1,
        )?;
        if slots > budget.max_channel_spectral_nodes {
            return Err(CoatingExchangeError::Budget);
        }
        checkpoint()?;
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(slots)
            .map_err(|_| CoatingExchangeError::Allocation)?;
        let wavelength_bounds_nm = emitter.wavelength_bounds_nm();
        let width =
            (wavelength_bounds_nm.1 - wavelength_bounds_nm.0) / quadrature.wavelength_bins as f64;
        if !width.is_finite() || width <= 0.0 {
            return Err(CoatingExchangeError::Quadrature);
        }
        let cosine_count = quadrature.cosine_bins as f64;
        let mut hash = DomainHasher::new("org.frankensim.render.coating-external-bath.v1");
        hash.update(emitter.optical_identity().as_bytes());
        hash.update(&area_m2.to_bits().to_le_bytes());
        for count in [
            film_count,
            usize::from(has_substrate),
            quadrature.wavelength_bins,
            quadrature.cosine_bins,
        ] {
            hash.update(&(count as u64).to_le_bytes());
        }
        for wavelength_index in 0..quadrature.wavelength_bins {
            let wavelength_nm = wavelength_bounds_nm.0 + (wavelength_index as f64 + 0.5) * width;
            let mut angular: [Sum; MAX_CHANNELS] = core::array::from_fn(|_| Sum::default());
            for cosine_index in 0..quadrature.cosine_bins {
                checkpoint()?;
                let cosine = (cosine_index as f64 + 0.5) / cosine_count;
                let absorbed = emitter.reciprocal_layer_absorption(wavelength_nm, cosine)?;
                for (sum, &[s, p]) in angular.iter_mut().zip(
                    absorbed
                        .film_absorbed()
                        .iter()
                        .chain(absorbed.substrate.iter()),
                ) {
                    // Sum polarizations before halving. Planck includes both;
                    // its 1/2 combines with the azimuthal 2*pi exactly once.
                    sum.add(positive_product(cosine, s + p)?)?;
                }
            }
            for sum in angular.iter().take(channel_count) {
                let angular_value = sum.value()?;
                let mean_projected = angular_value / cosine_count;
                if angular_value > 0.0 && mean_projected == 0.0 {
                    return Err(ThermalSurfaceError::Unrepresentable.into());
                }
                let extraction = PI * mean_projected;
                let (coefficient, weight) = scaled_factors(area_m2, extraction, width)?;
                for value in [wavelength_nm, extraction, coefficient, weight] {
                    hash.update(&value.to_bits().to_le_bytes());
                }
                nodes.push(SpectralNode {
                    wavelength_nm,
                    coefficient,
                    weight,
                });
            }
        }
        checkpoint()?;
        Ok(Self {
            nodes: nodes.into(),
            film_count,
            has_substrate,
            area_m2,
            quadrature,
            optical_evaluations,
            wavelength_bounds_nm,
            optical_identity: emitter.optical_identity(),
            identity: hash.finalize(),
        })
    }

    /// Signed exterior heat and frozen-optics secants at independent channel
    /// temperatures. `None` for the bath means zero incoming vacuum radiance;
    /// material temperatures always remain strictly positive. Substrate presence
    /// must exactly match preparation. Even a dark film keeps its temperature
    /// admission. Equality with the bath returns exact zero and the native
    /// analytic tangent, without subtracting gross emitted/absorbed powers.
    pub fn evaluate(
        &self,
        film_temperatures_k: &[f64],
        substrate_temperature_k: Option<f64>,
        bath_temperature_k: Option<f64>,
        mut checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<CoatingExternalBathRates, CoatingExchangeError> {
        checkpoint()?;
        if film_temperatures_k.len() != self.film_count
            || substrate_temperature_k.is_some() != self.has_substrate
        {
            return Err(CoatingExchangeError::TemperatureChannels);
        }
        if bath_temperature_k.is_some_and(|t| !t.is_finite() || t <= 0.0) {
            return Err(ThermalSurfaceError::Temperature.into());
        }
        let channel_count = self.film_count + usize::from(self.has_substrate);
        let mut temperatures = [0.0; MAX_CHANNELS];
        for (slot, &temperature) in temperatures.iter_mut().zip(
            film_temperatures_k
                .iter()
                .chain(substrate_temperature_k.iter()),
        ) {
            // Retain original Planck admission, including optically dark films.
            planck_radiative_exchange_w_and_w_k(
                temperature,
                bath_temperature_k,
                self.wavelength_bounds_nm.0,
                0.0,
                0.0,
            )?;
            *slot = temperature;
        }
        let mut powers: [Sum; MAX_CHANNELS] = core::array::from_fn(|_| Sum::default());
        let mut secants: [Sum; MAX_CHANNELS] = core::array::from_fn(|_| Sum::default());
        if channel_count != 0 {
            for spectral in self.nodes.chunks_exact(channel_count) {
                checkpoint()?;
                for (index, node) in spectral.iter().enumerate() {
                    let (power, secant) = planck_radiative_exchange_w_and_w_k(
                        temperatures[index],
                        bath_temperature_k,
                        node.wavelength_nm,
                        node.coefficient,
                        node.weight,
                    )?;
                    powers[index].add(power)?;
                    secants[index].add(secant)?;
                }
            }
        }
        let mut result = CoatingExternalBathRates {
            channels: [ZERO_RATE; MAX_CHANNELS],
            film_count: self.film_count,
            has_substrate: self.has_substrate,
            total_outward_power_w: 0.0,
        };
        let mut total = Sum::default();
        for (index, slot) in result.channels.iter_mut().take(channel_count).enumerate() {
            let outward_power_w = powers[index].value()?;
            let secant_w_k = secants[index].value()?;
            if outward_power_w != 0.0 && secant_w_k == 0.0 {
                return Err(ThermalSurfaceError::Unrepresentable.into());
            }
            *slot = CoatingExchangeRates {
                outward_power_w,
                secant_w_k,
            };
            total.add(outward_power_w)?;
        }
        result.total_outward_power_w = total.value()?;
        checkpoint()?;
        Ok(result)
    }

    #[must_use]
    pub const fn film_count(&self) -> usize {
        self.film_count
    }
    #[must_use]
    pub const fn has_substrate(&self) -> bool {
        self.has_substrate
    }
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    #[must_use]
    pub const fn quadrature(&self) -> CoatingExchangeQuadrature {
        self.quadrature
    }
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.wavelength_bounds_nm
    }
    #[must_use]
    pub const fn optical_evaluations(&self) -> usize {
        self.optical_evaluations
    }
    #[must_use]
    pub fn channel_spectral_nodes(&self) -> usize {
        self.nodes.len()
    }
    #[must_use]
    pub const fn optical_identity(&self) -> ContentHash {
        self.optical_identity
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

#[cfg(test)]
mod tests;
