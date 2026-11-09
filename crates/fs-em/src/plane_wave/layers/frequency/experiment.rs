//! One immutable physical experiment shared by peak and cancellation searches.

use super::*;

pub(super) struct FrequencyExperiment<'s, 'm> {
    source: &'s PlanarDielectricLayers<'m>,
    layers: [DielectricLayer<'m>; MAX_LAYERS],
    cosine: f64,
    polarization: Polarization,
    constraint: IncidenceFrequencyConstraint,
    observation: LayeredResonanceObservation,
    transverse: f64,
}

impl<'s, 'm> FrequencyExperiment<'s, 'm> {
    pub(super) fn new(
        source: &'s PlanarDielectricLayers<'m>,
        cosine: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
        observation: LayeredResonanceObservation,
    ) -> Result<Self, ReflectionError> {
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(bad("incidence cosine must lie in (0,1]"));
        }
        let mut layers = [DielectricLayer {
            response: source.exit,
            thickness_m: 0.0,
        }; MAX_LAYERS];
        for (target, layer) in layers.iter_mut().zip(&source.layers) {
            *target = layer.source;
        }
        Ok(Self {
            source,
            layers,
            cosine,
            polarization,
            constraint,
            observation,
            transverse: source.wave_number * det::sqrt((1.0 - cosine) * (1.0 + cosine)),
        })
    }

    pub(super) fn sample(
        &self,
        omega: f64,
        need_derivative: bool,
    ) -> Result<(C64, C64, f64), ReflectionError> {
        let sample = PlanarDielectricLayers::try_new_with_incident_response(
            self.source.incident_response,
            self.source.exit,
            omega,
            &self.layers[..self.source.layers.len()],
        )?;
        let cosine = match self.constraint {
            IncidenceFrequencyConstraint::FixedAngle => self.cosine,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber => {
                let sine = self.transverse / sample.wave_number;
                if !sine.is_finite() || sine >= 1.0 {
                    return Err(bad(
                        "fixed transverse wave number lost propagating incidence",
                    ));
                }
                positive(det::sqrt((1.0 - sine) * (1.0 + sine)))?
            }
        };
        if let Some(value) = absorption::search_observation(
            &sample,
            cosine,
            self.polarization,
            self.constraint,
            self.observation,
            need_derivative,
        )? {
            return Ok(value);
        }
        if !need_derivative {
            let response = sample.response(cosine, self.polarization)?;
            return Ok((
                select_carrier(response, self.observation)
                    .ok_or_else(|| bad("absorption observation was not resolved"))?,
                C64::ZERO,
                0.0,
            ));
        }
        let response = sample.frequency_dispersion(cosine, self.polarization, self.constraint)?;
        let jet = match self.observation {
            LayeredResonanceObservation::ReflectedElectricAmplitude => {
                Some(response.reflected_tangential_e())
            }
            LayeredResonanceObservation::TransmittedElectricAmplitude => {
                response.transmitted_tangential_e()
            }
            LayeredResonanceObservation::TransmittedMagneticAmplitude => {
                response.transmitted_conjugate_h()
            }
            LayeredResonanceObservation::TransmittedPowerFraction => {
                let power = response.transmitted_power_fraction()?.ok_or_else(|| {
                    bad("propagation underflow has no admitted transmitted-power slope")
                })?;
                let (amplitude, scaled_derivative) = power
                    .square_root_observation()
                    .ok_or_else(|| bad("positive resolved flux required for a power-peak query"))?;
                Some(LayerAmplitudeDerivative {
                    amplitude,
                    omega,
                    scaled_derivative,
                })
            }
            LayeredResonanceObservation::AbsorbedPowerFraction
            | LayeredResonanceObservation::LayerAbsorbedPowerFraction { .. }
            | LayeredResonanceObservation::AbsorbedPowerDensity
            | LayeredResonanceObservation::LayerAbsorbedPowerDensity { .. } => None,
        }
        .ok_or_else(|| bad("propagation underflow has no admitted transmission frequency jet"))?;
        Ok((
            jet.amplitude(),
            jet.derivative_s(),
            response.relative_equation_residual(),
        ))
    }
}
