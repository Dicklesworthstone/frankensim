//! Original coherent film field, power and heating frequency resonances.
//!
//! The native bounded E/H recursion owns every carrier value. Its tangent uses
//! the original dielectric jets, with no frequency differencing or growing
//! transfer exponential. Entire functions of the longitudinal phase squared
//! retain the analytic finite-film limit at critical refraction.
//!
//! exp(+i omega t) gives delay -Im(H'/H). This is a local phase slope, not
//! an energy or signal-front velocity. Transmitted electric/magnetic amplitude
//! maxima need not maximize transmitted power when the exit admittance varies.

use super::{
    DielectricLayer, LayeredReflection, MAX_LAYERS, PlanarDielectricLayers, Slab, admittance_state,
    bad, divided, maximum, positive, propagate_state,
};
use crate::plane_wave::{
    PlanarDielectricInterface, Polarization, ReflectionError, WaveMedium,
    dispersion::{BulkWaveDispersion, IncidenceFrequencyConstraint, PowerFractionDispersion},
    finite,
};
use fs_math::{c64::C64, dd::Dd, det};
use fs_phs::harmonic::resonance::{
    DrivenResonance, ResonanceControls, ResonanceSearchError, search_resonance,
};

mod absorption;
pub use absorption::{AbsorbedPowerDensityDispersion, LayeredAbsorptionDispersion};

/// One coherent field amplitude and its analytic first frequency derivative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerAmplitudeDerivative {
    amplitude: C64,
    omega: f64,
    scaled_derivative: C64,
}

impl LayerAmplitudeDerivative {
    pub const fn amplitude(self) -> C64 {
        self.amplitude
    }
    /// omega dH/domega, retaining optical frequencies without squaring them.
    pub const fn scaled_derivative(self) -> C64 {
        self.scaled_derivative
    }
    /// dH/domega, per rad/s.
    pub fn derivative_s(self) -> C64 {
        divided(self.scaled_derivative, self.omega)
    }
    /// Signed -d arg(H)/domega [s], absent at a zero or unrepresentable quotient.
    /// Reflection is referenced to the first face; transmission to the last.
    pub fn group_delay_s(self) -> Option<f64> {
        if self.amplitude == C64::ZERO {
            return None;
        }
        let delay = -(self.scaled_derivative / self.amplitude).im / self.omega;
        delay.is_finite().then_some(delay)
    }
}

/// Native passive carrier plus analytic derivatives of that same experiment.
#[derive(Clone, Copy, Debug)]
pub struct LayeredFrequencyResponse {
    carrier: LayeredReflection,
    constraint: IncidenceFrequencyConstraint,
    reflected: LayerAmplitudeDerivative,
    transmitted: Option<[LayerAmplitudeDerivative; 2]>,
    relative_equation_residual: f64,
}

impl LayeredFrequencyResponse {
    pub const fn carrier(self) -> LayeredReflection {
        self.carrier
    }
    pub const fn frequency_constraint(self) -> IncidenceFrequencyConstraint {
        self.constraint
    }
    pub const fn reflected_tangential_e(self) -> LayerAmplitudeDerivative {
        self.reflected
    }
    /// None when propagation underflowed. A zero rounded carrier does not prove
    /// that its derivative is zero, so no transmission jet is invented there.
    pub fn transmitted_tangential_e(self) -> Option<LayerAmplitudeDerivative> {
        self.transmitted.map(|value| value[0])
    }
    /// Exit H/(Y_inc E_inc), including the original incident admittance slope.
    /// The same propagation-underflow refusal applies to both transmitted fields.
    pub fn transmitted_conjugate_h(self) -> Option<LayerAmplitudeDerivative> {
        self.transmitted.map(|value| value[1])
    }

    /// Reflected normal power fraction and its complete analytic frequency slope.
    pub fn reflected_power_fraction(self) -> Result<PowerFractionDispersion, ReflectionError> {
        let r = self.reflected;
        PowerFractionDispersion::from_fields(
            self.carrier.reflected_power_fraction,
            r.omega,
            r.amplitude,
            r.amplitude,
            r.scaled_derivative,
            r.scaled_derivative,
        )
        .ok_or_else(|| bad("reflected film Poynting-flux tangent is unrepresentable"))
    }

    /// Transmitted normal energy-flux fraction. The slope uses both E and H,
    /// including the original dispersive incident/exit admittances. It is None
    /// when the native field tangent is unavailable after propagation underflow.
    pub fn transmitted_power_fraction(
        self,
    ) -> Result<Option<PowerFractionDispersion>, ReflectionError> {
        let Some([e, h]) = self.transmitted else {
            return Ok(None);
        };
        PowerFractionDispersion::from_fields(
            self.carrier.transmitted_power_fraction,
            e.omega,
            e.amplitude,
            h.amplitude,
            e.scaled_derivative,
            h.scaled_derivative,
        )
        .map(Some)
        .ok_or_else(|| bad("transmitted film Poynting-flux tangent is unrepresentable"))
    }
    /// Largest normalized backward defect of the bounded E/H transfer and
    /// tangent equations, audited with double-double product accumulation.
    /// This is not a forward amplitude error bound near a narrow resonance.
    pub const fn relative_equation_residual(self) -> f64 {
        self.relative_equation_residual
    }
}

/// Observation in the native one-sided experiment. Power objectives use a
/// phase-free square root, so their squared amplitude is the named physical
/// fraction or heat density, rather than an optical field intensity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LayeredResonanceObservation {
    /// |r| squared is the native reflected power fraction.
    ReflectedElectricAmplitude,
    /// |t_E| squared is an electric-amplitude ratio, not generally transmittance.
    TransmittedElectricAmplitude,
    /// |H_exit/(Y_inc E_inc)| squared is a magnetic-amplitude ratio.
    TransmittedMagneticAmplitude,
    /// Actual transmitted normal Poynting fraction. The shared search uses
    /// sqrt(T) as a phase-free observation, so observation.norm_sq() is T and
    /// the reported half-height crossings refer to transmitted POWER. Strictly
    /// positive resolved transmission is required at derivative samples.
    TransmittedPowerFraction,
    /// Power dissipated in all finite films, divided by incident normal power.
    /// This includes coherent standing waves and original material dispersion.
    AbsorbedPowerFraction,
    /// Power dissipated in one original film, divided by incident normal power.
    LayerAbsorbedPowerFraction { layer_index: usize },
    /// Actual total film heat per area for a fixed unit peak incident E [V/m].
    /// The squared observation has units W/m² and includes the incident-port
    /// admittance slope; it can peak at a different frequency from absorptance.
    AbsorbedPowerDensity,
    /// Actual heat per area in one original film, for fixed unit peak incident E.
    LayerAbsorbedPowerDensity { layer_index: usize },
}

impl<'a> PlanarDielectricLayers<'a> {
    /// Differentiate the original finite stack at its prepared carrier.
    ///
    /// The incident material jet must remain lossless, and the source's scalar
    /// permeability snapshot is fixed, as in BulkWaveDispersion. Both explicit
    /// incidence constraints are supported. Finite critical films use their
    /// analytic limit; an exactly critical exit half-space retains the native
    /// Fresnel-jet refusal. Zero-thickness films have exactly zero tangent.
    ///
    /// Source finite-band, derivative-domain and representability refusals
    /// propagate. Opaque films retain the native front-reflection limit, while
    /// transmitted derivatives become unavailable rather than silently zero.
    /// The stack recursion is O(layer_count) with no per-layer heap workspace;
    /// original constitutive evaluations retain their own cost. The native cap
    /// is 64 films.
    pub fn frequency_dispersion(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<LayeredFrequencyResponse, ReflectionError> {
        self.frequency_dispersion_impl(incidence_cosine, polarization, constraint, None)
    }

    fn frequency_dispersion_impl(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
        mut absorption: Option<&mut absorption::AbsorptionWork>,
    ) -> Result<LayeredFrequencyResponse, ReflectionError> {
        let carrier = self.response(incidence_cosine, polarization)?;
        let interface = PlanarDielectricInterface::try_new_with_incident_response(
            self.incident_response,
            self.exit,
        )?;
        let exit = interface.frequency_dispersion(
            self.omega,
            incidence_cosine,
            polarization,
            constraint,
        )?;
        let incident = exit.incident_dispersion();
        let incident_slope =
            incident.scaled_permittivity_derivatives_f_m()[0].re / self.incident_epsilon_f_m;
        let wave_slope = 1.0 + 0.5 * incident_slope;
        let sine_squared = (1.0 - incidence_cosine) * (1.0 + incidence_cosine);
        let cosine_slope = match constraint {
            IncidenceFrequencyConstraint::FixedAngle => 0.0,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber => {
                sine_squared * wave_slope / incidence_cosine
            }
        };
        if !wave_slope.is_finite() || !cosine_slope.is_finite() {
            return Err(bad("incident film frequency tangent is unrepresentable"));
        }
        if let Some(work) = absorption.as_deref_mut() {
            work.set_reference(
                self,
                incidence_cosine,
                cosine_slope,
                wave_slope,
                incident_slope,
                polarization,
            )?;
        }
        let mut state = admittance_state(self.exit_wave, incidence_cosine, polarization)?;
        let mut reflected_slope = exit.reflected_tangential_e().scaled_derivatives()[0];
        let mut transmitted = [
            exit.carrier().transmitted_tangential_e,
            state[1].scale(2.0) / (state[0] + state[1]),
        ];
        let mut transmitted_slope = [
            exit.transmitted_tangential_e().scaled_derivatives()[0],
            exit.transmitted_conjugate_h().scaled_derivatives()[0],
        ];
        let mut transmission_available = true;
        let mut residual: f64 = 0.0;
        for (index, layer) in self.layers.iter().enumerate().rev() {
            if layer.source.thickness_m == 0.0 {
                continue;
            }
            let slab = Slab::new(
                layer.wave,
                self.wave_number,
                incidence_cosine,
                polarization,
                layer.source.thickness_m,
            )?;
            if slab.scale == C64::ZERO {
                if absorption.is_some() {
                    return Err(bad(
                        "opaque propagation has no resolved interior absorption tangent",
                    ));
                }
                let front = PlanarDielectricInterface::try_new_with_incident_response(
                    self.incident_response,
                    layer.source.response,
                )?
                .frequency_dispersion(
                    self.omega,
                    incidence_cosine,
                    polarization,
                    constraint,
                )?;
                state = [slab.matrix[0][0], slab.matrix[1][0]];
                reflected_slope = front.reflected_tangential_e().scaled_derivatives()[0];
                transmitted = [C64::ZERO; 2];
                transmission_available = false;
                continue;
            }
            let material = layer.source.response.bulk_dispersion(self.omega)?;
            let tangent = slab_tangent(
                slab,
                layer.wave,
                material,
                incident,
                self.wave_number,
                incidence_cosine,
                cosine_slope,
                wave_slope,
                incident_slope,
                polarization,
                layer.source.thickness_m,
            )?;
            let tail = state;
            let tail_reflected_slope = reflected_slope;
            // A homogeneous state has a free common scalar. Choose its tangent
            // gauge so D(E+H)=0; carrier E/H still retain their native numerators.
            let tail_slope = reflected_slope * (tail[0] + tail[1]).scale(0.5);
            let tail_slope = [tail_slope, -tail_slope];
            let (next, gain) = propagate_state(slab, tail)?;
            let raw = [
                slab.matrix[0][0] * tail[0] + slab.matrix[0][1] * tail[1],
                slab.matrix[1][0] * tail[0] + slab.matrix[1][1] * tail[1],
            ];
            let normalization = positive(maximum(raw[0]).max(maximum(raw[1])))?;
            let mut slope = [C64::ZERO; 2];
            for row in 0..2 {
                let terms = [
                    (slab.matrix[row][0], tail_slope[0]),
                    (slab.matrix[row][1], tail_slope[1]),
                    (tangent[row][0], tail[0]),
                    (tangent[row][1], tail[1]),
                ];
                slope[row] = divided(products(&terms)?, normalization);
                if !finite(slope[row]) {
                    return Err(bad("film state tangent is unrepresentable"));
                }
                residual = residual.max(equation_residual(
                    next[row].scale(normalization),
                    &[
                        (slab.matrix[row][0], tail[0]),
                        (slab.matrix[row][1], tail[1]),
                    ],
                )?);
                residual =
                    residual.max(equation_residual(slope[row].scale(normalization), &terms)?);
            }
            let sum = next[0] + next[1];
            // Numerical slab scale and state normalization are frozen gauges.
            // The complete physical transfer tangent already includes all source
            // and reference-port variation, so neither gauge is differentiated.
            let gain_slope = -gain * ((slope[0] + slope[1]) / sum);
            reflected_slope =
                (products(&[(slope[0], next[1]), (-next[0], slope[1])])?.scale(2.0) / sum) / sum;
            if !finite(gain_slope) || !finite(reflected_slope) {
                return Err(bad("coherent film amplitude tangent is unrepresentable"));
            }
            if let Some(work) = absorption.as_deref_mut() {
                work.record(
                    index,
                    *layer,
                    material,
                    absorption::AbsorptionStep {
                        slab,
                        tail_reflection: (tail[0] - tail[1]) / (tail[0] + tail[1]),
                        tail_scaled_derivative: tail_reflected_slope,
                        gain,
                        gain_scaled_derivative: gain_slope,
                    },
                )?;
            }
            // A nonzero bounded slab scale can still yield an underflowed
            // composed gain; its lost transmission tangent is not zero.
            if gain == C64::ZERO {
                transmission_available = false;
            }
            if transmission_available {
                for component in 0..2 {
                    transmitted_slope[component] = products(&[
                        (gain_slope, transmitted[component]),
                        (gain, transmitted_slope[component]),
                    ])?;
                    let value = gain * transmitted[component];
                    if value == C64::ZERO
                        && gain != C64::ZERO
                        && transmitted[component] != C64::ZERO
                    {
                        transmission_available = false;
                    }
                    transmitted[component] = value;
                }
            }
            state = next;
        }
        let make = |amplitude, scaled_derivative| LayerAmplitudeDerivative {
            amplitude,
            omega: self.omega,
            scaled_derivative,
        };
        let reflected = make(carrier.reflected_tangential_e, reflected_slope);
        let transmitted = transmission_available.then(|| {
            [
                make(carrier.transmitted_tangential_e, transmitted_slope[0]),
                make(carrier.transmitted_conjugate_h, transmitted_slope[1]),
            ]
        });
        if !finite(reflected.derivative_s())
            || transmitted.is_some_and(|values| values.iter().any(|x| !finite(x.derivative_s())))
            || residual > 512.0 * f64::EPSILON
        {
            return Err(bad(
                "film frequency tangent or transfer residual is unrepresentable",
            ));
        }
        Ok(LayeredFrequencyResponse {
            carrier,
            constraint,
            reflected,
            transmitted,
            relative_equation_residual: residual,
        })
    }

    /// Refine a driven amplitude maximum and optional selected-amplitude
    /// half-height bandwidth, using the shared resonance-search owner.
    ///
    /// The caller interval must isolate an increasing-then-decreasing maximum;
    /// this is not a pole solve, mode count or guarantee of a unique peak.
    /// Every sample rebuilds the original material laws at its actual frequency.
    /// A frequency-prepared source that cannot admit another frequency refuses.
    /// FixedTransverseWavenumber retains k_parallel from this prepared carrier
    /// and incidence angle; loss of propagating incidence refuses.
    ///
    /// The selected field, normal-flux fraction or heat density remains the same
    /// experiment across the interval. Power fractions include port-admittance
    /// slopes; physical heat densities additionally hold unit peak incident E
    /// fixed. Absorption retains original material loss and coherent internal
    /// fields in every selected film. Power observations use a phase-free square
    /// root, with no optical phase assigned. Other choices measure field amplitudes.
    /// An observed half-height Q need not equal a pole-decay Q. The 64-film cap,
    /// shared total evaluation budget and before/after-query cancellation apply.
    pub fn driven_resonance(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
        observation: LayeredResonanceObservation,
        interval_rad_s: (f64, f64),
        controls: ResonanceControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenResonance, ResonanceSearchError<ReflectionError>> {
        if !incidence_cosine.is_finite() || incidence_cosine <= 0.0 || incidence_cosine > 1.0 {
            return Err(ResonanceSearchError::Source(bad(
                "incidence cosine must lie in (0,1]",
            )));
        }
        let mut layers = [DielectricLayer {
            response: self.exit,
            thickness_m: 0.0,
        }; MAX_LAYERS];
        for (target, source) in layers.iter_mut().zip(&self.layers) {
            *target = source.source;
        }
        let transverse =
            self.wave_number * det::sqrt((1.0 - incidence_cosine) * (1.0 + incidence_cosine));
        search_resonance(
            interval_rad_s,
            controls,
            cancelled,
            |omega, need_derivative, _| {
                let sample = Self::try_new_with_incident_response(
                    self.incident_response,
                    self.exit,
                    omega,
                    &layers[..self.layers.len()],
                )?;
                let cosine = match constraint {
                    IncidenceFrequencyConstraint::FixedAngle => incidence_cosine,
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber => {
                        let sine = transverse / sample.wave_number;
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
                    polarization,
                    constraint,
                    observation,
                    need_derivative,
                )? {
                    return Ok(value);
                }
                if !need_derivative {
                    let response = sample.response(cosine, polarization)?;
                    return Ok((
                        select_carrier(response, observation)
                            .ok_or_else(|| bad("absorption observation was not resolved"))?,
                        C64::ZERO,
                        0.0,
                    ));
                }
                let response = sample.frequency_dispersion(cosine, polarization, constraint)?;
                let jet = match observation {
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
                        let (amplitude, scaled_derivative) =
                            power.square_root_observation().ok_or_else(|| {
                                bad("positive resolved flux required for a power-peak query")
                            })?;
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
                .ok_or_else(|| {
                    bad("propagation underflow has no admitted transmission frequency jet")
                })?;
                Ok((
                    jet.amplitude(),
                    jet.derivative_s(),
                    response.relative_equation_residual(),
                ))
            },
        )
    }
}

fn select_carrier(
    value: LayeredReflection,
    observation: LayeredResonanceObservation,
) -> Option<C64> {
    Some(match observation {
        LayeredResonanceObservation::ReflectedElectricAmplitude => value.reflected_tangential_e,
        LayeredResonanceObservation::TransmittedElectricAmplitude => value.transmitted_tangential_e,
        LayeredResonanceObservation::TransmittedMagneticAmplitude => value.transmitted_conjugate_h,
        LayeredResonanceObservation::TransmittedPowerFraction => {
            C64::from_re(det::sqrt(value.transmitted_power_fraction))
        }
        LayeredResonanceObservation::AbsorbedPowerFraction
        | LayeredResonanceObservation::LayerAbsorbedPowerFraction { .. }
        | LayeredResonanceObservation::AbsorbedPowerDensity
        | LayeredResonanceObservation::LayerAbsorbedPowerDensity { .. } => return None,
    })
}

// D denotes omega*d/domega. Only the original epsilon varies with frequency;
// native scalar permeability is a fixed material-state snapshot.
fn slab_tangent(
    slab: Slab,
    medium: WaveMedium,
    material: BulkWaveDispersion<'_>,
    incident: BulkWaveDispersion<'_>,
    wave_number: f64,
    cosine: f64,
    cosine_slope: f64,
    wave_slope: f64,
    incident_slope: f64,
    polarization: Polarization,
    thickness: f64,
) -> Result<[[C64; 2]; 2], ReflectionError> {
    let q = medium.normal(cosine)?;
    let depth = positive(wave_number * thickness)?;
    let delta = q.scale(depth);
    let z = delta * delta;
    let epsilon_slope = divided(
        material.scaled_permittivity_derivatives_f_m()[0],
        incident.permittivity_f_m().re,
    ) - medium.epsilon_ratio.scale(incident_slope);
    let q_squared_slope =
        epsilon_slope.scale(medium.mu_ratio) + C64::from_re(2.0 * cosine * cosine_slope);
    let (s, s_slope) = match polarization {
        Polarization::Te => (
            C64::from_re(medium.mu_ratio * cosine),
            C64::from_re(medium.mu_ratio * cosine_slope),
        ),
        Polarization::Tm => (
            medium.epsilon_ratio.scale(cosine),
            epsilon_slope.scale(cosine) + medium.epsilon_ratio.scale(cosine_slope),
        ),
    };
    if s == C64::ZERO {
        return Err(bad("film frequency coupling is unrepresentable"));
    }
    let a = s.scale(depth);
    let b = (q / s) * delta;
    let da = a.scale(wave_slope) + s_slope.scale(depth);
    let db = b.scale(wave_slope) + (q_squared_slope / s).scale(depth) - b * (s_slope / s);
    let dz = z.scale(2.0 * wave_slope) + q_squared_slope.scale(depth).scale(depth);
    if [delta, z, a, b, da, db, dz].iter().any(|&x| !finite(x)) {
        return Err(bad("film longitudinal frequency jet is unrepresentable"));
    }
    let (p_sinc, p_sinc_z) = bounded_sinc_jet(q, wave_number, thickness, delta, z)?;
    let alpha = slab.transfer_scale_over_propagator;
    let diagonal = -(alpha * p_sinc).scale(0.5) * dz;
    let times_i = |value: C64| C64::new(-value.im, value.re);
    let sa = times_i((alpha * p_sinc) * da + ((alpha * a) * p_sinc_z) * dz);
    let sb = times_i((alpha * p_sinc) * db + ((alpha * b) * p_sinc_z) * dz);
    let (upper, lower) = match polarization {
        Polarization::Te => (sa, sb),
        Polarization::Tm => (sb, sa),
    };
    let result = [[diagonal, upper], [lower, diagonal]];
    if result.iter().flatten().any(|&x| !finite(x)) {
        return Err(bad("bounded film transfer tangent is unrepresentable"));
    }
    Ok(result)
}

// p=exp(-i delta) always decays. Compute p*sinc(delta) and
// p*d(sinc(delta))/d(delta²), without differentiating this numerical gauge.
fn bounded_sinc_jet(
    q: C64,
    wave_number: f64,
    thickness: f64,
    delta: C64,
    z: C64,
) -> Result<(C64, C64), ReflectionError> {
    let attenuation = det::exp(delta.im);
    let phase = Dd::from_f64(q.re) * (Dd::from_f64(wave_number) * Dd::from_f64(thickness));
    let (sine, cosine) = det::sin_cos_dd(phase);
    let (sine, cosine) = (sine.to_f64(), cosine.to_f64());
    let p = C64::new(attenuation * cosine, -attenuation * sine);
    let result = if z.abs() <= 1.0 {
        let [sinc, sinc_z, _] = super::dispersion::sinc_series(z);
        (p * sinc, p * sinc_z)
    } else {
        let power = det::exp(2.0 * delta.im);
        let loss = -det::expm1(2.0 * delta.im);
        let difference = C64::new(
            loss + 2.0 * power * sine * sine,
            2.0 * power * sine * cosine,
        );
        let p_cosine = C64::new(
            loss + 2.0 * power * cosine * cosine,
            -2.0 * power * sine * cosine,
        )
        .scale(0.5);
        let p_sinc = (difference / delta) * C64::new(0.0, -0.5);
        (p_sinc, ((p_cosine - p_sinc) / z).scale(0.5))
    };
    if !finite(result.0) || !finite(result.1) {
        return Err(bad("bounded film entire-function jet is unrepresentable"));
    }
    Ok(result)
}

fn products(terms: &[(C64, C64)]) -> Result<C64, ReflectionError> {
    let mut real = Dd::ZERO;
    let mut imaginary = Dd::ZERO;
    for &(a, b) in terms {
        real = real + Dd::from_f64(a.re) * Dd::from_f64(b.re)
            - Dd::from_f64(a.im) * Dd::from_f64(b.im);
        imaginary = imaginary
            + Dd::from_f64(a.re) * Dd::from_f64(b.im)
            + Dd::from_f64(a.im) * Dd::from_f64(b.re);
    }
    let value = C64::new(real.to_f64(), imaginary.to_f64());
    if finite(value) {
        Ok(value)
    } else {
        Err(bad("film tangent product is unrepresentable"))
    }
}

fn equation_residual(lhs: C64, terms: &[(C64, C64)]) -> Result<f64, ReflectionError> {
    let rhs = products(terms)?;
    let mut scale = maximum(lhs).max(maximum(rhs));
    for &(a, b) in terms {
        scale = scale.max(maximum(a) * maximum(b));
    }
    if scale == 0.0 {
        return Ok(0.0);
    }
    if !scale.is_finite() {
        return Err(bad("film tangent residual scale is unrepresentable"));
    }
    let mut denominator = divided(lhs, scale).abs();
    for &(a, b) in terms {
        denominator += (maximum(a) * maximum(b)) / scale;
    }
    let residual = (divided(lhs, scale) - divided(rhs, scale)).abs() / denominator;
    if residual.is_finite() {
        Ok(residual)
    } else {
        Err(bad("film tangent residual is unrepresentable"))
    }
}
