//! Original positive film loss channels and their complete frequency tangent.

use super::super::PreparedLayer;
use super::*;

mod loss;
use loss::{SlabLossDirection, scaled_loss_derivatives};

/// Cycle-mean material heat per area for a fixed complex peak incident E.
/// This is monochromatic physical dissipation, not a thermal-state update.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbedPowerDensityDispersion {
    power_w_m2: f64,
    omega: f64,
    scaled_derivative_w_m2: f64,
}

impl AbsorbedPowerDensityDispersion {
    pub const fn power_w_m2(self) -> f64 {
        self.power_w_m2
    }
    /// omega dP/domega, holding the supplied physical incident E fixed.
    pub const fn scaled_derivative_w_m2(self) -> f64 {
        self.scaled_derivative_w_m2
    }
    /// dP/domega [W s/m²].
    pub fn derivative_w_s_m2(self) -> f64 {
        self.scaled_derivative_w_m2 / self.omega
    }

    pub(super) fn square_root_observation(self) -> Option<(C64, C64)> {
        if self.power_w_m2 == 0.0 {
            return None;
        }
        let a = det::sqrt(self.power_w_m2);
        let d = (self.scaled_derivative_w_m2 / a) * 0.5;
        (d.is_finite() && (d / self.omega).is_finite())
            .then_some((C64::from_re(a), C64::from_re(d)))
    }
}

/// Native scattering and absorption in every original film, in authored order.
/// The positive even/odd loss channels retain coherent interference before
/// powers are formed. Neither the carrier nor its slope uses 1-R-T.
#[derive(Clone, Debug)]
pub struct LayeredAbsorptionDispersion {
    response: LayeredFrequencyResponse,
    layers: Vec<PowerFractionDispersion>,
    total: PowerFractionDispersion,
    omega: f64,
    incident_admittance_s: f64,
    incident_admittance_scaled_derivative_s: f64,
    power_balance_relative_residual: f64,
}

impl LayeredAbsorptionDispersion {
    pub const fn response(&self) -> LayeredFrequencyResponse {
        self.response
    }
    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    pub fn layer_absorption(&self, layer_index: usize) -> Option<PowerFractionDispersion> {
        self.layers.get(layer_index).copied()
    }
    pub const fn total_absorption(&self) -> PowerFractionDispersion {
        self.total
    }
    /// Independent normalized defect of D(R+T+A)=0. A is computed from the
    /// original positive film loss, not defined by this identity.
    pub const fn power_balance_relative_residual(&self) -> f64 {
        self.power_balance_relative_residual
    }

    /// Heat and its complete slope for a FIXED physical complex peak incident E.
    /// Includes the original incident admittance slope, unlike absorptance alone.
    pub fn total_power_density(
        &self,
        incident_tangential_e_v_m: C64,
    ) -> Result<AbsorbedPowerDensityDispersion, ReflectionError> {
        self.power_density(self.total, incident_tangential_e_v_m)
    }

    pub fn layer_power_density(
        &self,
        layer_index: usize,
        incident_tangential_e_v_m: C64,
    ) -> Result<AbsorbedPowerDensityDispersion, ReflectionError> {
        let fraction = self
            .layer_absorption(layer_index)
            .ok_or_else(|| bad("absorption layer is outside the original stack"))?;
        self.power_density(fraction, incident_tangential_e_v_m)
    }

    fn power_density(
        &self,
        fraction: PowerFractionDispersion,
        source: C64,
    ) -> Result<AbsorbedPowerDensityDispersion, ReflectionError> {
        if !finite(source) {
            return Err(bad("incident absorption source must be finite"));
        }
        let power = source_power(self.incident_admittance_s, fraction.fraction(), source).to_f64();
        let slope = (source_power(
            self.incident_admittance_s,
            fraction.scaled_derivative(),
            source,
        ) + source_power(
            self.incident_admittance_scaled_derivative_s,
            fraction.fraction(),
            source,
        ))
        .to_f64();
        if !power.is_finite() || !slope.is_finite() || !(slope / self.omega).is_finite() {
            return Err(bad("physical absorbed-power tangent is unrepresentable"));
        }
        Ok(AbsorbedPowerDensityDispersion {
            power_w_m2: power,
            omega: self.omega,
            scaled_derivative_w_m2: slope,
        })
    }
}

impl PlanarDielectricLayers<'_> {
    /// Differentiate original material loss AND the coherent field in each film.
    ///
    /// The existing bounded native transfer/tangent sweep is reused, followed
    /// by one forward source sweep. Cost is O(layers), with fixed 64-film scratch
    /// storage and one bounded output vector. There is no depth mesh, frequency
    /// differencing, extra solve per film, or growing transfer exponential.
    ///
    /// Zero-thickness and lossless films retain exactly zero heat. A lossless
    /// carrier with a nonzero constitutive loss slope refuses. Propagation or
    /// composed-source underflow also refuses the interior tangent; the native
    /// carrier and front-reflection API remain independently available.
    pub fn frequency_absorption(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<LayeredAbsorptionDispersion, ReflectionError> {
        let mut work = AbsorptionWork::new();
        let response = self.frequency_dispersion_impl(
            incidence_cosine,
            polarization,
            constraint,
            Some(&mut work),
        )?;
        work.finish(response, self.layers.len())
    }
}

pub(super) fn search_observation(
    source: &PlanarDielectricLayers<'_>,
    cosine: f64,
    polarization: Polarization,
    constraint: IncidenceFrequencyConstraint,
    observation: LayeredResonanceObservation,
    need_derivative: bool,
) -> Result<Option<(C64, C64, f64)>, ReflectionError> {
    use LayeredResonanceObservation::*;
    if !matches!(
        observation,
        AbsorbedPowerFraction
            | LayerAbsorbedPowerFraction { .. }
            | AbsorbedPowerDensity
            | LayerAbsorbedPowerDensity { .. }
    ) {
        return Ok(None);
    }
    if let LayerAbsorbedPowerFraction { layer_index } | LayerAbsorbedPowerDensity { layer_index } =
        observation
    {
        if layer_index >= source.layers.len() {
            return Err(bad(
                "resonance observation layer is outside the original stack",
            ));
        }
    }
    // The single O(N) preparation supplies every selected channel. Keeping its
    // original positive carrier at half-height samples avoids a second field
    // owner or dividing a possibly underflowed physical heat by incident flux.
    let response = source.frequency_absorption(cosine, polarization, constraint)?;
    let (power, slope) = match observation {
        AbsorbedPowerFraction | LayerAbsorbedPowerFraction { .. } => {
            let jet = match observation {
                LayerAbsorbedPowerFraction { layer_index } => {
                    response.layer_absorption(layer_index).ok_or_else(|| {
                        bad("resonance observation layer is outside the original stack")
                    })?
                }
                _ => response.total_absorption(),
            };
            if !need_derivative {
                (C64::from_re(det::sqrt(jet.fraction())), C64::ZERO)
            } else {
                jet.square_root_observation().ok_or_else(|| {
                    bad("positive resolved absorption is required for a peak tangent")
                })?
            }
        }
        AbsorbedPowerDensity | LayerAbsorbedPowerDensity { .. } => {
            let jet = match observation {
                LayerAbsorbedPowerDensity { layer_index } => {
                    response.layer_power_density(layer_index, C64::ONE)?
                }
                _ => response.total_power_density(C64::ONE)?,
            };
            if !need_derivative {
                (C64::from_re(det::sqrt(jet.power_w_m2())), C64::ZERO)
            } else {
                jet.square_root_observation().ok_or_else(|| {
                    bad("positive resolved heating is required for a peak tangent")
                })?
            }
        }
        _ => return Ok(None),
    };
    Ok(Some((
        power,
        divided(slope, source.omega),
        response
            .response
            .relative_equation_residual()
            .max(response.power_balance_relative_residual()),
    )))
}

pub(super) struct AbsorptionStep {
    pub slab: Slab,
    pub tail_reflection: C64,
    pub tail_scaled_derivative: C64,
    pub gain: C64,
    pub gain_scaled_derivative: C64,
}

#[derive(Clone, Copy)]
struct LossRecord {
    gain: C64,
    gain_slope: C64,
    reflection: C64,
    reflection_slope: C64,
    loss: [f64; 2],
    loss_slope: [f64; 2],
}
impl LossRecord {
    const IDENTITY: Self = Self {
        gain: C64::ONE,
        gain_slope: C64::ZERO,
        reflection: C64::ZERO,
        reflection_slope: C64::ZERO,
        loss: [0.0; 2],
        loss_slope: [0.0; 2],
    };
}

pub(super) struct AbsorptionWork {
    records: [LossRecord; MAX_LAYERS],
    omega: f64,
    wave_number: f64,
    wave_slope: f64,
    epsilon: f64,
    epsilon_slope: f64,
    cosine: f64,
    cosine_slope: f64,
    admittance: f64,
    admittance_slope: f64,
    transverse_electric: bool,
}

impl AbsorptionWork {
    fn new() -> Self {
        Self {
            records: [LossRecord::IDENTITY; MAX_LAYERS],
            omega: 0.0,
            wave_number: 0.0,
            wave_slope: 0.0,
            epsilon: 0.0,
            epsilon_slope: 0.0,
            cosine: 0.0,
            cosine_slope: 0.0,
            admittance: 0.0,
            admittance_slope: 0.0,
            transverse_electric: true,
        }
    }

    pub(super) fn set_reference(
        &mut self,
        source: &PlanarDielectricLayers<'_>,
        cosine: f64,
        cosine_slope: f64,
        wave_slope: f64,
        epsilon_slope: f64,
        polarization: Polarization,
    ) -> Result<(), ReflectionError> {
        self.omega = source.omega;
        self.wave_number = source.wave_number;
        self.wave_slope = wave_slope;
        self.epsilon = source.incident_epsilon_f_m;
        self.epsilon_slope = epsilon_slope;
        self.cosine = cosine;
        self.cosine_slope = cosine_slope;
        self.transverse_electric = polarization == Polarization::Te;
        let normal = positive(det::sqrt(self.epsilon) / det::sqrt(source.incident_mu_h_m))?;
        let (y, log_slope) = match polarization {
            Polarization::Te => (normal * cosine, 0.5 * epsilon_slope + cosine_slope / cosine),
            Polarization::Tm => (normal / cosine, 0.5 * epsilon_slope - cosine_slope / cosine),
        };
        self.admittance = positive(y)?;
        self.admittance_slope = y * log_slope;
        if !self.admittance_slope.is_finite() {
            return Err(bad(
                "incident absorbed-power port tangent is unrepresentable",
            ));
        }
        Ok(())
    }

    pub(super) fn record(
        &mut self,
        index: usize,
        layer: PreparedLayer<'_>,
        material: BulkWaveDispersion<'_>,
        step: AbsorptionStep,
    ) -> Result<(), ReflectionError> {
        if step.gain == C64::ZERO {
            return Err(bad(
                "underflowed film gain has no resolved absorption tangent",
            ));
        }
        let material_slope = material.scaled_permittivity_derivatives_f_m()[0];
        let slope = if layer.wave.epsilon_f_m.im == 0.0 {
            if material_slope.im != 0.0 {
                return Err(bad(
                    "lossless carrier has a nonzero constitutive loss tangent",
                ));
            }
            [0.0; 2]
        } else {
            let ratio_slope = divided(material_slope, self.epsilon)
                - layer.wave.epsilon_ratio.scale(self.epsilon_slope);
            let q = layer.wave.normal(self.cosine)?;
            if q == C64::ZERO {
                return Err(bad("lossy film normal-wave tangent is unresolved"));
            }
            let dq = (ratio_slope.scale(layer.wave.mu_ratio)
                + C64::from_re(2.0 * self.cosine * self.cosine_slope))
                / q;
            let depth = positive(self.wave_number * layer.source.thickness_m)?;
            scaled_loss_derivatives(SlabLossDirection {
                normal: q,
                normal_scaled_derivative: dq.scale(0.5),
                permittivity_ratio: layer.wave.epsilon_ratio,
                permittivity_ratio_scaled_derivative: ratio_slope,
                permeability_ratio: layer.wave.mu_ratio,
                cosine: self.cosine,
                cosine_scaled_derivative: self.cosine_slope,
                depth,
                depth_scaled_derivative: depth * self.wave_slope,
                phase_parts: Dd::from_f64(q.re)
                    * (Dd::from_f64(self.wave_number) * Dd::from_f64(layer.source.thickness_m)),
                loss: [step.slab.even_loss, step.slab.odd_loss],
                transverse_electric: self.transverse_electric,
            })
            .ok_or_else(|| bad("positive film loss tangent is unrepresentable"))?
        };
        self.records[index] = LossRecord {
            gain: step.gain,
            gain_slope: step.gain_scaled_derivative,
            reflection: step.tail_reflection,
            reflection_slope: step.tail_scaled_derivative,
            loss: [step.slab.even_loss, step.slab.odd_loss],
            loss_slope: slope,
        };
        Ok(())
    }

    fn finish(
        self,
        response: LayeredFrequencyResponse,
        count: usize,
    ) -> Result<LayeredAbsorptionDispersion, ReflectionError> {
        let mut layers = Vec::with_capacity(count);
        let mut incoming = C64::ONE;
        let mut incoming_slope = C64::ZERO;
        let mut total_slope = Dd::ZERO;
        for record in &self.records[..count] {
            let outgoing = incoming * record.gain;
            if outgoing == C64::ZERO && incoming != C64::ZERO {
                return Err(bad(
                    "composed interior source underflow has no absorption tangent",
                ));
            }
            let outgoing_slope =
                products(&[(incoming_slope, record.gain), (incoming, record.gain_slope)])?;
            let opposite = record.reflection * outgoing;
            let opposite_slope = products(&[
                (record.reflection_slope, outgoing),
                (record.reflection, outgoing_slope),
            ])?;
            let channels = [incoming + opposite, incoming - opposite];
            let slopes = [
                incoming_slope + opposite_slope,
                incoming_slope - opposite_slope,
            ];
            let mut fraction = 0.0;
            let mut derivative = Dd::ZERO;
            for mode in 0..2 {
                let u = channels[mode];
                let du = slopes[mode];
                let lambda = record.loss[mode];
                let dlambda = record.loss_slope[mode];
                fraction += super::super::loss_amplitude(u, lambda).norm_sq();
                derivative = derivative
                    + product3(dlambda, u.re, u.re) * Dd::from_f64(0.5)
                    + product3(dlambda, u.im, u.im) * Dd::from_f64(0.5)
                    + product3(lambda, u.re, du.re)
                    + product3(lambda, u.im, du.im);
            }
            let jet =
                PowerFractionDispersion::from_scaled(fraction, self.omega, derivative.to_f64())
                    .ok_or_else(|| bad("coherent absorbed-power tangent is unrepresentable"))?;
            total_slope = total_slope + derivative;
            layers.push(jet);
            incoming = outgoing;
            incoming_slope = outgoing_slope;
        }
        let total = PowerFractionDispersion::from_scaled(
            response.carrier.absorbed_power_fraction,
            self.omega,
            total_slope.to_f64(),
        )
        .ok_or_else(|| bad("total absorbed-power tangent is unrepresentable"))?;
        let reflected = response.reflected_power_fraction()?.scaled_derivative();
        let transmitted = response
            .transmitted_power_fraction()?
            .ok_or_else(|| bad("transmitted flux tangent is unavailable for the absorption audit"))?
            .scaled_derivative();
        let absorbed = total.scaled_derivative();
        let scale = 1.0_f64
            .max(reflected.abs())
            .max(transmitted.abs())
            .max(absorbed.abs());
        let residual = (Dd::from_f64(reflected / scale)
            + Dd::from_f64(transmitted / scale)
            + Dd::from_f64(absorbed / scale))
        .to_f64()
        .abs()
            / (1.0 / scale
                + reflected.abs() / scale
                + transmitted.abs() / scale
                + absorbed.abs() / scale);
        if !residual.is_finite() || residual > 4096.0 * f64::EPSILON {
            return Err(bad("differentiated film Poynting balance failed"));
        }
        Ok(LayeredAbsorptionDispersion {
            response,
            layers,
            total,
            omega: self.omega,
            incident_admittance_s: self.admittance,
            incident_admittance_scaled_derivative_s: self.admittance_slope,
            power_balance_relative_residual: residual,
        })
    }
}

// Multiply the smallest and largest magnitudes first. This preserves tiny loss
// multiplied by large coherent fields without first squaring those fields.
fn product3(a: f64, b: f64, c: f64) -> Dd {
    let mut factors = [a, b, c];
    factors.sort_by(|a, b| a.abs().total_cmp(&b.abs()));
    (Dd::from_f64(factors[0]) * Dd::from_f64(factors[2])) * Dd::from_f64(factors[1])
}

// Include the physical source before any squaring or small-power product.
fn source_power(admittance: f64, fraction: f64, source: C64) -> Dd {
    let term = |component| {
        let mut factors = [0.5, admittance, fraction, component, component];
        factors.sort_by(|a, b| a.abs().total_cmp(&b.abs()));
        let (mut left, mut right) = (0, factors.len());
        let mut value = Dd::from_f64(1.0);
        while left < right {
            let factor = if value.hi.abs() >= 1.0 {
                let f = factors[left];
                left += 1;
                f
            } else {
                right -= 1;
                factors[right]
            };
            value = value * Dd::from_f64(factor);
        }
        value
    };
    term(source.re) + term(source.im)
}

#[cfg(test)]
mod scale_tests;
