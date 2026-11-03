//! Multiple transmitting bodies with mutual screening and thermal radiation.
//!
//! Bodies are ordered from left to right; EVERY selected emitter port faces
//! LEFT. Each is an internally coherent, isothermal, reciprocal film stack
//! between vacuum half spaces. Propagating powers, rather than amplitudes,
//! cross the phase-averaged gaps. This is the coherent/incoherent construction
//! of Byrnes, *Multilayer optical calculations*, arXiv:1603.02720v5, section 6:
//! <https://arxiv.org/abs/1603.02720>. Local thermal sources retain the original
//! Kirchhoff absorption owner; no angular/polarization averages are multiplied.
//!
//! A positive forward elimination and backward sweep solve both directional
//! intensities. Escape is carried independently of rounded unit reflection.
//! Prefix source radiance is normalized by its positive escape, so a weakly
//! coupled cavity does not require huge unnormalized internal source factors.
//! Every body/body and body/bath coefficient is retained; exterior throughput
//! is separate. Optical preparation is O(N²) per angular/spectral mode with
//! only 2N original coherent face solves. Temperature queries allocate nothing.
//!
//! This is finite-band, far-field transport through coextensive parallel
//! specular bodies. The caller must justify phase averaging. There is no gap
//! phase, near-field transfer, finite-edge escape, rough scattering or internal
//! radiative redistribution between nonisothermal films inside one body.

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

/// Fixed scratch/result bound; callers also supply their own smaller budget.
pub const MAX_RADIATIVE_STACK_BODIES: usize = 16;
const MAX_PAIRS: usize = MAX_RADIATIVE_STACK_BODIES * (MAX_RADIATIVE_STACK_BODIES - 1) / 2;
const MAX_TERMINALS: usize = MAX_RADIATIVE_STACK_BODIES + 2;
const MAX_CHANNELS: usize = MAX_TERMINALS * (MAX_TERMINALS - 1) / 2;
const ROUNDING_ALLOWANCE: f64 = 512.0 * f64::EPSILON * (MAX_COHERENT_STACK_LAYERS + 1) as f64;
const ZERO_RATE: CoatingExchangeRates = CoatingExchangeRates {
    outward_power_w: 0.0,
    secant_w_k: 0.0,
};
const ZERO_PORT: ReciprocalPortPower = ReciprocalPortPower {
    reflected: [0.0; 2],
    transmitted: [0.0; 2],
    absorbed: [0.0; 2],
};

/// Admission bounds before allocation or optical work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoatingRadiativeStackBudget {
    pub max_bodies: usize,
    /// `(N + 2) * (N + 1) / 2` retained channels per wavelength node.
    pub max_channel_spectral_nodes: usize,
    /// `2 * N * wavelength_bins * cosine_bins` coherent solves.
    pub max_optical_evaluations: usize,
}

/// One reciprocal edge; positive watts leave `first` and enter `second`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoatingRadiativePairRate {
    pub first: usize,
    pub second: usize,
    pub exchange: CoatingExchangeRates,
}

/// All signed exchange channels, in bounded allocation-free storage.
#[derive(Clone, Debug, PartialEq)]
pub struct CoatingRadiativeStackRates {
    body_count: usize,
    pairs: [CoatingRadiativePairRate; MAX_PAIRS],
    baths: [[CoatingExchangeRates; 2]; MAX_RADIATIVE_STACK_BODIES],
    /// Left-to-right exterior transmission, never heat deposited in a body.
    pub bath_transfer: CoatingExchangeRates,
    body_outward: [f64; MAX_RADIATIVE_STACK_BODIES],
    bath_outward: [f64; 2],
}
impl CoatingRadiativeStackRates {
    /// Lexicographic unordered body pairs, with `first < second`.
    #[must_use]
    pub fn body_pairs(&self) -> &[CoatingRadiativePairRate] {
        &self.pairs[..self.body_count * (self.body_count - 1) / 2]
    }
    /// Per-body outward power toward `[left bath, right bath]`.
    #[must_use]
    pub fn body_baths(&self) -> &[[CoatingExchangeRates; 2]] {
        &self.baths[..self.body_count]
    }
    #[must_use]
    pub fn body_outward_power_w(&self) -> &[f64] {
        &self.body_outward[..self.body_count]
    }
    /// Includes direct bath throughput. All body and bath powers sum to zero
    /// up to the explicitly ordered floating-point accumulation.
    #[must_use]
    pub const fn bath_outward_power_w(&self) -> [f64; 2] {
        self.bath_outward
    }
}

/// Shared immutable spectral response of a complete ordered radiative stack.
#[derive(Clone, Debug)]
pub struct CoatingRadiativeStackResponse {
    nodes: Arc<[SpectralNode]>,
    body_count: usize,
    optical_identities: [ContentHash; MAX_RADIATIVE_STACK_BODIES],
    area_m2: f64,
    quadrature: CoatingExchangeQuadrature,
    optical_evaluations: usize,
    wavelength_bounds_nm: (f64, f64),
    optical_identity: ContentHash,
    identity: ContentHash,
}
impl CoatingRadiativeStackResponse {
    /// Emitters are ordered left to right, each selected port facing LEFT.
    /// Each whole body must be internally isothermal and vacuum/vacuum across
    /// its entire original band. Only their common band is integrated.
    pub fn from_emitters(
        emitters: &[CoatingThermalEmission],
        area_m2: f64,
        quadrature: CoatingExchangeQuadrature,
        budget: CoatingRadiativeStackBudget,
        mut checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<Self, CoatingExchangeError> {
        let n = emitters.len();
        if n == 0 || n > MAX_RADIATIVE_STACK_BODIES || n > budget.max_bodies {
            return Err(CoatingExchangeError::Budget);
        }
        let channels = channel_count(n);
        let slots = channels
            .checked_mul(quadrature.wavelength_bins)
            .ok_or(CoatingExchangeError::Budget)?;
        let optical_evaluations = admit(
            quadrature,
            CoatingExchangeBudget {
                max_spectral_nodes: budget.max_channel_spectral_nodes / channels,
                max_optical_evaluations: budget.max_optical_evaluations,
            },
            area_m2,
            2 * n,
        )?;
        let mut band = (0.0_f64, f64::INFINITY);
        let mut optical_identities = [ContentHash([0; 32]); MAX_RADIATIVE_STACK_BODIES];
        let mut optical_hash =
            DomainHasher::new("org.frankensim.render.coating-radiative-stack-optics.v1");
        optical_hash.update(&(n as u64).to_le_bytes());
        for (index, emitter) in emitters.iter().enumerate() {
            if !emitter.is_vacuum_transmitting() {
                return Err(CoatingExchangeError::VacuumTransmittingRequired);
            }
            validate_isothermal(emitter)?;
            let bounds = emitter.wavelength_bounds_nm();
            band = (band.0.max(bounds.0), band.1.min(bounds.1));
            optical_identities[index] = emitter.optical_identity();
            optical_hash.update(emitter.optical_identity().as_bytes());
        }
        if band.0 >= band.1 {
            return Err(CoatingExchangeError::WavelengthDomain);
        }
        let width = (band.1 - band.0) / quadrature.wavelength_bins as f64;
        if !width.is_finite() || width <= 0.0 {
            return Err(CoatingExchangeError::Quadrature);
        }
        checkpoint()?;
        let optical_identity = optical_hash.finalize();
        let mut hash = DomainHasher::new("org.frankensim.render.coating-radiative-stack.v1");
        hash.update(optical_identity.as_bytes());
        for v in [area_m2, band.0, band.1] {
            hash.update(&v.to_bits().to_le_bytes());
        }
        for v in [quadrature.wavelength_bins, quadrature.cosine_bins] {
            hash.update(&(v as u64).to_le_bytes());
        }
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(slots)
            .map_err(|_| CoatingExchangeError::Allocation)?;
        let cosine_count = quadrature.cosine_bins as f64;
        for wavelength_index in 0..quadrature.wavelength_bins {
            let wavelength_nm = band.0 + (wavelength_index as f64 + 0.5) * width;
            let mut angular: [Sum; MAX_CHANNELS] = core::array::from_fn(|_| Sum::default());
            for cosine_index in 0..quadrature.cosine_bins {
                let cosine = (cosine_index as f64 + 0.5) / cosine_count;
                let mut ports = [[ZERO_PORT; 2]; MAX_RADIATIVE_STACK_BODIES];
                for (faces, emitter) in ports.iter_mut().zip(emitters) {
                    for (opposite, face) in faces.iter_mut().enumerate() {
                        checkpoint()?;
                        *face = emitter
                            .reciprocal_vacuum_port_power(wavelength_nm, cosine, opposite != 0)?
                            .ok_or(CoatingExchangeError::VacuumTransmittingRequired)?;
                    }
                }
                let fractions = exchange_fractions(&ports[..n], &mut checkpoint)?;
                for (sum, [s, p]) in angular[..channels].iter_mut().zip(fractions) {
                    sum.add(positive_product(cosine, s + p)?)?;
                }
            }
            for sum in &angular[..channels] {
                let total = sum.value()?;
                let mean = total / cosine_count;
                if total > 0.0 && mean == 0.0 {
                    return Err(ThermalSurfaceError::Unrepresentable.into());
                }
                // Planck counts both polarizations; its half combines with
                // the azimuthal 2*pi, giving pi*(s+p)*mu exactly once.
                let extraction = PI * mean;
                let (coefficient, weight) = scaled_factors(area_m2, extraction, width)?;
                for v in [wavelength_nm, extraction, coefficient, weight] {
                    hash.update(&v.to_bits().to_le_bytes());
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
            body_count: n,
            optical_identities,
            area_m2,
            quadrature,
            optical_evaluations,
            wavelength_bounds_nm: band,
            optical_identity,
            identity: hash.finalize(),
        })
    }

    /// Signed heat and positive frozen-optics secants. `None` means zero bath
    /// radiance. Equal reservoir temperatures produce exact zero exchange.
    /// No optical solves or allocation occur during this query.
    pub fn evaluate(
        &self,
        temperatures_k: &[f64],
        left_bath_temperature_k: Option<f64>,
        right_bath_temperature_k: Option<f64>,
        mut checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<CoatingRadiativeStackRates, CoatingExchangeError> {
        if temperatures_k.len() != self.body_count {
            return Err(CoatingExchangeError::TemperatureChannels);
        }
        if temperatures_k
            .iter()
            .copied()
            .chain(
                [left_bath_temperature_k, right_bath_temperature_k]
                    .into_iter()
                    .flatten(),
            )
            .any(|t| !t.is_finite() || t <= 0.0)
        {
            return Err(ThermalSurfaceError::Temperature.into());
        }
        let n = self.body_count;
        let channels = channel_count(n);
        let mut temperature_pairs = [(None, None); MAX_CHANNELS];
        let mut edge = 0;
        for first in 0..n {
            for second in first + 1..n {
                temperature_pairs[edge] =
                    (Some(temperatures_k[first]), Some(temperatures_k[second]));
                edge += 1;
            }
        }
        for &temperature in temperatures_k {
            for bath in [left_bath_temperature_k, right_bath_temperature_k] {
                temperature_pairs[edge] = (Some(temperature), bath);
                edge += 1;
            }
        }
        temperature_pairs[edge] = (left_bath_temperature_k, right_bath_temperature_k);
        let mut powers: [Sum; MAX_CHANNELS] = core::array::from_fn(|_| Sum::default());
        let mut secants: [Sum; MAX_CHANNELS] = core::array::from_fn(|_| Sum::default());
        for nodes in self.nodes.chunks_exact(channels) {
            checkpoint()?;
            for (edge, node) in nodes.iter().enumerate() {
                let (a, b) = temperature_pairs[edge];
                let (power, secant) = reservoir_exchange(a, b, *node)?;
                powers[edge].add(power)?;
                secants[edge].add(secant)?;
            }
        }
        let mut rates = [ZERO_RATE; MAX_CHANNELS];
        for edge in 0..channels {
            rates[edge] = CoatingExchangeRates {
                outward_power_w: powers[edge].value()?,
                secant_w_k: secants[edge].value()?,
            };
            if rates[edge].secant_w_k < 0.0
                || (rates[edge].outward_power_w != 0.0 && rates[edge].secant_w_k == 0.0)
            {
                return Err(ThermalSurfaceError::Unrepresentable.into());
            }
        }
        let mut result = CoatingRadiativeStackRates {
            body_count: n,
            pairs: [CoatingRadiativePairRate {
                first: 0,
                second: 0,
                exchange: ZERO_RATE,
            }; MAX_PAIRS],
            baths: [[ZERO_RATE; 2]; MAX_RADIATIVE_STACK_BODIES],
            bath_transfer: rates[channels - 1],
            body_outward: [0.0; MAX_RADIATIVE_STACK_BODIES],
            bath_outward: [0.0; 2],
        };
        let mut bodies: [Sum; MAX_RADIATIVE_STACK_BODIES] =
            core::array::from_fn(|_| Sum::default());
        let mut baths: [Sum; 2] = core::array::from_fn(|_| Sum::default());
        edge = 0;
        for first in 0..n {
            for second in first + 1..n {
                let rate = rates[edge];
                result.pairs[edge] = CoatingRadiativePairRate {
                    first,
                    second,
                    exchange: rate,
                };
                bodies[first].add(rate.outward_power_w)?;
                bodies[second].add(-rate.outward_power_w)?;
                edge += 1;
            }
        }
        for (body, sums) in bodies[..n].iter_mut().enumerate() {
            for bath in 0..2 {
                let rate = rates[edge];
                result.baths[body][bath] = rate;
                sums.add(rate.outward_power_w)?;
                baths[bath].add(-rate.outward_power_w)?;
                edge += 1;
            }
            result.body_outward[body] = sums.value()?;
        }
        baths[0].add(result.bath_transfer.outward_power_w)?;
        baths[1].add(-result.bath_transfer.outward_power_w)?;
        result.bath_outward = [baths[0].value()?, baths[1].value()?];
        checkpoint()?;
        Ok(result)
    }

    #[must_use]
    pub const fn body_count(&self) -> usize {
        self.body_count
    }
    #[must_use]
    pub fn optical_identities(&self) -> &[ContentHash] {
        &self.optical_identities[..self.body_count]
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
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    #[must_use]
    pub const fn quadrature(&self) -> CoatingExchangeQuadrature {
        self.quadrature
    }
    #[must_use]
    pub const fn optical_evaluations(&self) -> usize {
        self.optical_evaluations
    }
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.wavelength_bounds_nm
    }
    #[must_use]
    pub fn channel_spectral_nodes(&self) -> usize {
        self.nodes.len()
    }
}

const fn channel_count(n: usize) -> usize {
    (n + 2) * (n + 1) / 2
}

fn reservoir_exchange(
    a: Option<f64>,
    b: Option<f64>,
    node: SpectralNode,
) -> Result<(f64, f64), ThermalSurfaceError> {
    let (first, second, sign) = match (a, b) {
        (Some(first), second) => (first, second, 1.0),
        (None, Some(second)) => (second, None, -1.0),
        (None, None) => return Ok((0.0, 0.0)),
    };
    let (power, secant) = planck_radiative_exchange_w_and_w_k(
        first,
        second,
        node.wavelength_nm,
        node.coefficient,
        node.weight,
    )?;
    Ok((sign * power, secant))
}

#[derive(Clone, Copy, Default)]
struct Sweep {
    left: f64,
    own: f64,
    right: f64,
    prefix_source: f64,
    body_source: f64,
}

/// Positive product divided by a positive denominator when the final result
/// is bounded. Dividing the smaller numerator first avoids product underflow;
/// in this probability-valued use its intermediate cannot overflow.
fn probability_product(x: f64, y: f64, denominator: f64) -> Result<f64, ThermalSurfaceError> {
    if x == 0.0 || y == 0.0 {
        return Ok(0.0);
    }
    positive_product(x.min(y) / denominator, x.max(y))
}

fn exchange_fractions(
    ports: &[[ReciprocalPortPower; 2]],
    mut checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
) -> Result<[[f64; 2]; MAX_CHANNELS], CoatingExchangeError> {
    let n = ports.len();
    let mut result = [[0.0; 2]; MAX_CHANNELS];
    for faces in ports {
        for pol in 0..2 {
            for port in faces {
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
                    return Err(ThermalSurfaceError::Unrepresentable.into());
                }
            }
            let (a, b) = (faces[0].transmitted[pol], faces[1].transmitted[pol]);
            if (a - b).abs() > ROUNDING_ALLOWANCE * a.max(b) + 8.0 * f64::from_bits(1) {
                return Err(ThermalSurfaceError::Unrepresentable.into());
            }
        }
    }
    for pol in 0..2 {
        let mut sweep = [Sweep::default(); MAX_RADIATIVE_STACK_BODIES];
        let (mut prefix_reflection, mut prefix_escape) = (0.0_f64, 1.0_f64);
        for (i, [left, right]) in ports.iter().enumerate() {
            let (a, b, t) = (
                left.absorbed[pol],
                right.absorbed[pol],
                left.transmitted[pol],
            );
            let escape = a + t;
            let scale = prefix_escape.max(escape);
            if scale == 0.0 {
                // A wholly lossless trapped gap has no incoming thermal source.
                // The new body's other side can still emit or absorb normally.
                sweep[i].body_source = f64::from(b > 0.0);
                prefix_reflection = right.reflected[pol];
                prefix_escape = b;
                continue;
            }
            let denominator = prefix_escape / scale + prefix_reflection * (escape / scale);
            if !denominator.is_finite() || denominator <= 0.0 {
                return Err(ThermalSurfaceError::Unrepresentable.into());
            }
            let transmission_ratio = (t / scale) / denominator;
            let left_weight = (prefix_escape / scale) / denominator;
            let own_weight = (prefix_reflection * (a / scale)) / denominator;
            let right_weight = prefix_reflection * transmission_ratio;
            let escape_scale = b.max(prefix_escape).max(a);
            let (new_escape, prefix_source, body_source) = if escape_scale == 0.0 {
                (0.0, 0.0, 0.0)
            } else {
                let previous = prefix_escape / escape_scale;
                let left_absorption = prefix_reflection * (a / escape_scale);
                let direct = b / escape_scale;
                let scaled_escape = direct + transmission_ratio * (previous + left_absorption);
                let new_escape = positive_product(escape_scale, scaled_escape)?;
                if new_escape == 0.0 {
                    (0.0, 0.0, 0.0)
                } else {
                    (
                        new_escape,
                        probability_product(transmission_ratio, previous, scaled_escape)?,
                        direct / scaled_escape
                            + probability_product(
                                transmission_ratio,
                                left_absorption,
                                scaled_escape,
                            )?,
                    )
                }
            };
            sweep[i] = Sweep {
                left: left_weight,
                own: own_weight,
                right: right_weight,
                prefix_source,
                body_source,
            };
            // This bounded reflected contribution may round away under its
            // large direct term; the independent positive escape retains loss.
            prefix_reflection = right.reflected[pol] + right_weight * t;
            prefix_escape = new_escape;
        }
        let mut absorption = [[0.0; MAX_TERMINALS]; MAX_TERMINALS];
        for source in 0..n + 2 {
            checkpoint()?;
            let mut prefix = [0.0; MAX_RADIATIVE_STACK_BODIES];
            let mut radiance = f64::from(source == n);
            for i in 0..n {
                prefix[i] = radiance;
                radiance = positive_product(sweep[i].prefix_source, radiance)?
                    + if source == i {
                        sweep[i].body_source
                    } else {
                        0.0
                    };
            }
            let mut from_right = f64::from(source == n + 1);
            for i in (0..n).rev() {
                let [left, right] = ports[i];
                let own = f64::from(source == i);
                let from_left = positive_product(sweep[i].left, prefix[i])?
                    + positive_product(sweep[i].right, from_right)?
                    + sweep[i].own * own;
                absorption[i][source] = positive_product(left.absorbed[pol], from_left)?
                    + positive_product(right.absorbed[pol], from_right)?;
                if i == n - 1 {
                    absorption[n + 1][source] = positive_product(left.transmitted[pol], from_left)?
                        + positive_product(right.reflected[pol], from_right)?
                        + right.absorbed[pol] * own;
                }
                from_right = positive_product(left.reflected[pol], from_left)?
                    + positive_product(left.transmitted[pol], from_right)?
                    + left.absorbed[pol] * own;
            }
            absorption[n][source] = from_right;
        }
        let tolerance = ROUNDING_ALLOWANCE * (2 * n + 1) as f64;
        let reciprocal = |first: usize, second: usize| -> Result<f64, ThermalSurfaceError> {
            let (a, b) = (absorption[second][first], absorption[first][second]);
            if !a.is_finite()
                || !b.is_finite()
                || a < 0.0
                || b < 0.0
                || (a - b).abs() > tolerance * a.max(b) + 32.0 * f64::from_bits(1)
            {
                return Err(ThermalSurfaceError::Unrepresentable);
            }
            Ok(a)
        };
        let mut edge = 0;
        for first in 0..n {
            for second in first + 1..n {
                result[edge][pol] = reciprocal(first, second)?;
                edge += 1;
            }
        }
        for body in 0..n {
            for bath in n..n + 2 {
                result[edge][pol] = reciprocal(body, bath)?;
                edge += 1;
            }
        }
        result[edge][pol] = reciprocal(n, n + 1)?;
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
