//! Complete finite-band radiation exchange for two transmitting planar bodies.
//!
//! Each body is an isothermal, internally coherent reciprocal film stack, with
//! vacuum on both sides. The selected sides face a unit-view-factor vacuum gap;
//! the other sides face independent isotropic black baths or empty vacuum.
//! Gap crossings add powers after phase averaging, while all interference
//! inside each body remains in the original coherent optical solve. This is
//! the coherent-stack/incoherent-space construction of Byrnes, *Multilayer
//! optical calculations*, arXiv:1603.02720v5, section 6, equation (25):
//! <https://arxiv.org/abs/1603.02720>. Thermal sources use the original local
//! Kirchhoff owner in [`CoatingThermalEmission`].
//!
//! At a matched wavelength, cosine and polarization, the right/left gap powers
//! satisfy `I+ = T1*BL + R1*I- + A1*B1` and
//! `I- = T2*BR + R2*I+ + A2*B2`. Eliminating these intensities gives six positive
//! pair/bath exchange channels. Positive `A+T` supplies the gap escape, avoiding
//! a catastrophic `1-R1*R2` subtraction. No product of independently averaged
//! angular or polarization responses replaces these matched modes.
//!
//! The caller must justify phase averaging; no supplied distance establishes
//! it. Infinite/coextensive parallel plates, specular propagating modes and
//! fixed optical laws are assumed. This is not near-field exchange, coherent
//! gap interference, finite-edge escape, rough scattering or internal exchange
//! between nonisothermal films. Only the common original source band is used.

use super::coating_emission::ReciprocalPortPower;
use super::coating_exchange::{
    SpectralNode, Sum, admit, positive_product, scaled_factors, validate_isothermal,
};
use super::{
    CoatingExchangeBudget, CoatingExchangeError, CoatingExchangeQuadrature, CoatingExchangeRates,
    CoatingThermalEmission, ThermalSurfaceError, planck_radiative_exchange_w_and_w_k,
};
use crate::conductor::MAX_COHERENT_STACK_LAYERS;
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};
use std::sync::Arc;

const CHANNELS: usize = 6;
const ROUNDING_ALLOWANCE: f64 = 512.0 * f64::EPSILON * (MAX_COHERENT_STACK_LAYERS + 1) as f64;

/// Storage for six spectral channels and four coherent solves per angular node.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoatingTransmittingPairBudget {
    pub max_channel_spectral_nodes: usize,
    pub max_optical_evaluations: usize,
}

/// Signed rates on the six edges of the two-body/two-bath radiation network.
/// A native thermal owner applies `pair` conservatively between the bodies and
/// the four body/bath rates as external heat. Bath throughput heats neither body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoatingTransmittingPairRates {
    /// Positive watts leave the first body and enter the second body.
    pub pair: CoatingExchangeRates,
    /// Positive watts leave the first body for [left bath, right bath].
    pub first_baths: [CoatingExchangeRates; 2],
    /// Positive watts leave the second body for [left bath, right bath].
    pub second_baths: [CoatingExchangeRates; 2],
    /// Positive watts pass from left bath to right bath, bypassing body heat.
    /// Two vacuum baths have zero power and zero conductance.
    pub bath_transfer: CoatingExchangeRates,
    body_outward_power_w: [f64; 2],
    bath_outward_power_w: [f64; 2],
}
impl CoatingTransmittingPairRates {
    /// Net signed heat leaving [first body, second body].
    #[must_use]
    pub const fn body_outward_power_w(&self) -> [f64; 2] {
        self.body_outward_power_w
    }

    /// Net signed heat leaving [left bath, right bath], including throughput.
    /// Together with body powers these close the four-terminal energy balance.
    #[must_use]
    pub const fn bath_outward_power_w(&self) -> [f64; 2] {
        self.bath_outward_power_w
    }
}

/// Reusable complete response for two isothermal vacuum/vacuum film stacks.
///
/// Both input emitters' selected sides face the gap. The opposite face of each
/// same authored stack is solved internally, preserving asymmetric absorption.
/// Four optical solves at each wavelength/cosine serve all six channels. Each
/// temperature query reuses original stable Planck differences/secants with no
/// optical solves or allocations. Clones share immutable spectral storage.
#[derive(Clone, Debug)]
pub struct CoatingTransmittingPairResponse {
    nodes: Arc<[SpectralNode]>,
    area_m2: f64,
    quadrature: CoatingExchangeQuadrature,
    optical_evaluations: usize,
    wavelength_bounds_nm: (f64, f64),
    optical_identities: [ContentHash; 2],
    optical_identity: ContentHash,
    identity: ContentHash,
}
impl CoatingTransmittingPairResponse {
    /// Prepare the common supported band and hemispherical midpoint rule.
    /// Both half spaces of each stack must be vacuum across the entire band.
    /// The shared facing area is physical m², with no distance attenuation or
    /// invented view factor. Budgets are admitted before storage or solves;
    /// cancellation and any original optical refusal publish no partial owner.
    pub fn from_emitters(
        first: &CoatingThermalEmission,
        second: &CoatingThermalEmission,
        area_m2: f64,
        quadrature: CoatingExchangeQuadrature,
        budget: CoatingTransmittingPairBudget,
        mut checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<Self, CoatingExchangeError> {
        let optical_evaluations = admit(
            quadrature,
            CoatingExchangeBudget {
                max_spectral_nodes: budget.max_channel_spectral_nodes / CHANNELS,
                max_optical_evaluations: budget.max_optical_evaluations,
            },
            area_m2,
            4,
        )?;
        let slots = quadrature
            .wavelength_bins
            .checked_mul(CHANNELS)
            .ok_or(CoatingExchangeError::Budget)?;
        if !first.is_vacuum_transmitting() || !second.is_vacuum_transmitting() {
            return Err(CoatingExchangeError::VacuumTransmittingRequired);
        }
        validate_isothermal(first)?;
        validate_isothermal(second)?;
        let first_band = first.wavelength_bounds_nm();
        let second_band = second.wavelength_bounds_nm();
        let wavelength_bounds_nm = (
            first_band.0.max(second_band.0),
            first_band.1.min(second_band.1),
        );
        if wavelength_bounds_nm.0 >= wavelength_bounds_nm.1 {
            return Err(CoatingExchangeError::WavelengthDomain);
        }
        let width =
            (wavelength_bounds_nm.1 - wavelength_bounds_nm.0) / quadrature.wavelength_bins as f64;
        if !width.is_finite() || width <= 0.0 {
            return Err(CoatingExchangeError::Quadrature);
        }
        checkpoint()?;
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(slots)
            .map_err(|_| CoatingExchangeError::Allocation)?;
        let optical_identities = [first.optical_identity(), second.optical_identity()];
        let mut optical_hash =
            DomainHasher::new("org.frankensim.render.coating-transmitting-pair-optics.v1");
        for identity in optical_identities {
            optical_hash.update(identity.as_bytes());
        }
        let optical_identity = optical_hash.finalize();
        let mut hash = DomainHasher::new("org.frankensim.render.coating-transmitting-pair.v1");
        hash.update(optical_identity.as_bytes());
        for value in [area_m2, wavelength_bounds_nm.0, wavelength_bounds_nm.1] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        for count in [quadrature.wavelength_bins, quadrature.cosine_bins] {
            hash.update(&(count as u64).to_le_bytes());
        }
        let cosine_count = quadrature.cosine_bins as f64;
        for wavelength_index in 0..quadrature.wavelength_bins {
            let wavelength_nm = wavelength_bounds_nm.0 + (wavelength_index as f64 + 0.5) * width;
            let mut angular: [Sum; CHANNELS] = core::array::from_fn(|_| Sum::default());
            for cosine_index in 0..quadrature.cosine_bins {
                let cosine = (cosine_index as f64 + 0.5) / cosine_count;
                let mut ports = [ReciprocalPortPower {
                    reflected: [0.0; 2],
                    transmitted: [0.0; 2],
                    absorbed: [0.0; 2],
                }; 4];
                for (port, (emitter, opposite)) in ports.iter_mut().zip([
                    (first, false),
                    (first, true),
                    (second, false),
                    (second, true),
                ]) {
                    checkpoint()?;
                    *port = emitter
                        .reciprocal_vacuum_port_power(wavelength_nm, cosine, opposite)?
                        .ok_or(CoatingExchangeError::VacuumTransmittingRequired)?;
                }
                let fractions = exchange_fractions(ports)?;
                for (sum, [s, p]) in angular.iter_mut().zip(fractions) {
                    // Planck contains both polarizations: its 1/2 combines
                    // with the azimuthal 2*pi exactly once.
                    sum.add(positive_product(cosine, s + p)?)?;
                }
            }
            for sum in angular {
                let total = sum.value()?;
                let mean = total / cosine_count;
                if total > 0.0 && mean == 0.0 {
                    return Err(ThermalSurfaceError::Unrepresentable.into());
                }
                let extraction = PI * mean;
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
            area_m2,
            quadrature,
            optical_evaluations,
            wavelength_bounds_nm,
            optical_identities,
            optical_identity,
            identity: hash.finalize(),
        })
    }

    /// Evaluate signed watts and positive frozen-optics W/K on every channel.
    /// `None` is zero incoming bath radiance; every supplied temperature must
    /// be finite and strictly positive, even for a transparent body or port.
    /// Equal temperatures give exact zero heat and the analytic tangent. The
    /// bath-throughput diagnostic never contributes to either body's heat.
    pub fn evaluate(
        &self,
        first_temperature_k: f64,
        second_temperature_k: f64,
        left_bath_temperature_k: Option<f64>,
        right_bath_temperature_k: Option<f64>,
        mut checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<CoatingTransmittingPairRates, CoatingExchangeError> {
        let first = Some(first_temperature_k);
        let second = Some(second_temperature_k);
        let left = left_bath_temperature_k;
        let right = right_bath_temperature_k;
        if [first, second, left, right]
            .into_iter()
            .flatten()
            .any(|temperature| !temperature.is_finite() || temperature <= 0.0)
        {
            return Err(ThermalSurfaceError::Temperature.into());
        }
        let temperatures = [
            (first, second),
            (first, left),
            (first, right),
            (second, left),
            (second, right),
            (left, right),
        ];
        let mut powers: [Sum; CHANNELS] = core::array::from_fn(|_| Sum::default());
        let mut secants: [Sum; CHANNELS] = core::array::from_fn(|_| Sum::default());
        for nodes in self.nodes.chunks_exact(CHANNELS) {
            checkpoint()?;
            for (channel, node) in nodes.iter().enumerate() {
                let (a, b) = temperatures[channel];
                let (power, secant) = reservoir_exchange(a, b, *node)?;
                powers[channel].add(power)?;
                secants[channel].add(secant)?;
            }
        }
        let mut rates = [CoatingExchangeRates {
            outward_power_w: 0.0,
            secant_w_k: 0.0,
        }; CHANNELS];
        for (channel, rate) in rates.iter_mut().enumerate() {
            rate.outward_power_w = powers[channel].value()?;
            rate.secant_w_k = secants[channel].value()?;
            if rate.secant_w_k < 0.0 || (rate.outward_power_w != 0.0 && rate.secant_w_k == 0.0) {
                return Err(ThermalSurfaceError::Unrepresentable.into());
            }
        }
        let p = rates.map(|rate| rate.outward_power_w);
        let result = CoatingTransmittingPairRates {
            pair: rates[0],
            first_baths: [rates[1], rates[2]],
            second_baths: [rates[3], rates[4]],
            bath_transfer: rates[5],
            body_outward_power_w: [sum3(p[0], p[1], p[2])?, sum3(-p[0], p[3], p[4])?],
            bath_outward_power_w: [sum3(-p[1], -p[3], p[5])?, sum3(-p[2], -p[4], -p[5])?],
        };
        checkpoint()?;
        Ok(result)
    }

    #[must_use]
    pub const fn optical_identities(&self) -> [ContentHash; 2] {
        self.optical_identities
    }
    #[must_use]
    pub const fn optical_identity(&self) -> ContentHash {
        self.optical_identity
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.wavelength_bounds_nm
    }
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    #[must_use]
    pub const fn quadrature(&self) -> CoatingExchangeQuadrature {
        self.quadrature
    }
    /// Four original optical solves per wavelength/cosine node.
    #[must_use]
    pub const fn optical_evaluations(&self) -> usize {
        self.optical_evaluations
    }
    /// Six retained channels per wavelength node, including bath throughput.
    #[must_use]
    pub fn channel_spectral_nodes(&self) -> usize {
        self.nodes.len()
    }
}

fn sum3(a: f64, b: f64, c: f64) -> Result<f64, ThermalSurfaceError> {
    let mut sum = Sum::default();
    for value in [a, b, c] {
        sum.add(value)?;
    }
    sum.value()
}

fn reservoir_exchange(
    first: Option<f64>,
    second: Option<f64>,
    node: SpectralNode,
) -> Result<(f64, f64), ThermalSurfaceError> {
    let (temperature, other, sign) = match (first, second) {
        (Some(first), second) => (first, second, 1.0),
        (None, Some(second)) => (second, None, -1.0),
        (None, None) => return Ok((0.0, 0.0)),
    };
    let (power, secant) = planck_radiative_exchange_w_and_w_k(
        temperature,
        other,
        node.wavelength_nm,
        node.coefficient,
        node.weight,
    )?;
    Ok((sign * power, secant))
}

fn exchange_fractions(
    ports: [ReciprocalPortPower; 4],
) -> Result<[[f64; 2]; CHANNELS], ThermalSurfaceError> {
    for port in ports {
        for pol in 0..2 {
            let values = [
                port.reflected[pol],
                port.transmitted[pol],
                port.absorbed[pol],
            ];
            if values
                .into_iter()
                .any(|v| !v.is_finite() || v < 0.0 || v > 1.0 + ROUNDING_ALLOWANCE)
                || (values.into_iter().sum::<f64>() - 1.0).abs() > ROUNDING_ALLOWANCE
            {
                return Err(ThermalSurfaceError::Unrepresentable);
            }
        }
    }
    let [first, first_outer, second, second_outer] = ports;
    let mut channels = [[0.0; 2]; CHANNELS];
    for pol in 0..2 {
        // Both half spaces are identical vacuum and the constitutive owners
        // reciprocal. Keep the selected side's original positive T as the
        // common transmission; refuse a mismatch beyond solver roundoff.
        for (gap, outer) in [(first, first_outer), (second, second_outer)] {
            let a = gap.transmitted[pol];
            let b = outer.transmitted[pol];
            if (a - b).abs() > ROUNDING_ALLOWANCE * a.max(b) + 8.0 * f64::from_bits(1) {
                return Err(ThermalSurfaceError::Unrepresentable);
            }
        }
        let a = [first.absorbed[pol], second.absorbed[pol]];
        let t = [first.transmitted[pol], second.transmitted[pol]];
        let escape = [a[0] + t[0], a[1] + t[1]];
        let large = escape[0].max(escape[1]);
        channels[1][pol] = first_outer.absorbed[pol];
        channels[4][pol] = second_outer.absorbed[pol];
        if large == 0.0 {
            continue;
        }
        let small = escape[0].min(escape[1]);
        let scaled_denominator = 1.0 + small / large - small;
        if !scaled_denominator.is_finite() || scaled_denominator <= 0.0 {
            return Err(ThermalSurfaceError::Unrepresentable);
        }
        // Every numerator factor is <= large. Dividing before multiplying
        // retains representable work even when A1*A2 or T1*T2 underflows.
        let quotient = |x: f64, y: f64| -> Result<f64, ThermalSurfaceError> {
            if x == 0.0 || y == 0.0 {
                return Ok(0.0);
            }
            let scale = (x.max(y) / large) / scaled_denominator;
            if scale <= 0.0 || !scale.is_finite() {
                return Err(ThermalSurfaceError::Unrepresentable);
            }
            positive_product(x.min(y), scale)
        };
        channels[0][pol] = quotient(a[0], a[1])?;
        if second.reflected[pol] > 0.0 {
            channels[1][pol] += positive_product(quotient(a[0], t[0])?, second.reflected[pol])?;
        }
        channels[2][pol] = quotient(a[0], t[1])?;
        channels[3][pol] = quotient(a[1], t[0])?;
        if first.reflected[pol] > 0.0 {
            channels[4][pol] += positive_product(quotient(a[1], t[1])?, first.reflected[pol])?;
        }
        channels[5][pol] = quotient(t[0], t[1])?;
    }
    Ok(channels)
}

#[cfg(test)]
mod tests;
