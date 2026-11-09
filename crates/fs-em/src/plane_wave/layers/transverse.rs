//! Fixed-frequency spatial dispersion of the original coherent planar stack.
//!
//! Derivatives follow signed k_parallel in one fixed incidence plane, using the
//! parent's fixed tangential TE/TM axes. The same scalar second-order algebra,
//! entire finite-film propagator and bounded E/H recurrence as frequency
//! dispersion retain internal critical layers without a square-root cusp.
//! Material values stay at the original prepared frequency; no material
//! frequency derivative or inferred broadband law participates.
//!
//! For exp(+i omega t - i k_parallel x), the local lateral phase shift is
//! +d arg(H)/dk_parallel. This is the narrow-spectrum Artmann displacement
//! along the interface, not a finite-beam centroid approximation at arbitrary
//! width, a transverse spin-Hall shift, or a displacement at an exterior plane.
//! Full finite beams remain with the original angular-spectrum field owner.
//! Spatial/angular beam-shift distinction: Töppel, Ornigotti and Aiello (2013),
//! <https://arxiv.org/abs/1307.6057>.

use super::{
    LayeredReflection, PlanarDielectricLayers, Slab, bad,
    dispersion::{Jet, propagate_jets, reflection_jet, slab_jets_from_normal},
    divided, maximum,
};
use crate::plane_wave::{
    Polarization, ReflectionError, WaveMedium, angular_spectrum::modes::incident_cosine, finite,
    positive,
};
use fs_math::c64::C64;

type Result<T> = core::result::Result<T, ReflectionError>;

/// Original complex field amplitude and its fixed-frequency k_parallel jet.
/// The parameter is a signed angular wave number [rad/m], not incidence angle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransverseAmplitudeDispersion {
    amplitude: C64,
    incident_wave_number: f64,
    scaled_derivatives: [C64; 2],
}

/// Local phase and amplitude changes along the declared incidence plane.
/// The phase gradient can represent a narrow reflected beam's lateral shift;
/// finite angular spread, angular filtering and singular zeros require the
/// original full angular-spectrum field, not this local interpretation alone.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransverseLogarithmicDispersion {
    /// +d arg(H)/dk_parallel [m], measured parallel to the interface.
    pub lateral_shift_m: f64,
    /// +d² arg(H)/dk_parallel² [m²], the spatial phase curvature.
    pub phase_curvature_m2: f64,
    /// d ln|H|/dk_parallel [m]; controls local angular-spectrum filtering.
    pub log_amplitude_derivative_m: f64,
    /// d² ln|H|/dk_parallel² [m²].
    pub log_amplitude_second_derivative_m2: f64,
}

impl TransverseAmplitudeDispersion {
    fn from_jet(amplitude: C64, incident_wave_number: f64, jet: Jet) -> Result<Self> {
        let value = Self {
            amplitude,
            incident_wave_number,
            scaled_derivatives: jet.derivatives(),
        };
        if !finite(amplitude)
            || !jet.is_finite()
            || !finite(value.derivative_m())
            || !finite(value.second_derivative_m2())
        {
            return Err(bad("transverse field derivatives are unrepresentable"));
        }
        Ok(value)
    }
    pub const fn amplitude(self) -> C64 {
        self.amplitude
    }
    /// [k_i dH/dk_parallel, k_i² d²H/dk_parallel²].
    /// Both entries are ordinary derivatives, not Taylor coefficients.
    pub const fn scaled_derivatives(self) -> [C64; 2] {
        self.scaled_derivatives
    }
    pub fn derivative_m(self) -> C64 {
        divided(self.scaled_derivatives[0], self.incident_wave_number)
    }
    pub fn second_derivative_m2(self) -> C64 {
        divided(
            divided(self.scaled_derivatives[1], self.incident_wave_number),
            self.incident_wave_number,
        )
    }
    /// No logarithmic phase/shift is assigned to an exact Brewster/Fresnel
    /// zero or an unrepresentable quotient. Complex derivatives remain valid.
    pub fn logarithmic_dispersion(self) -> Option<TransverseLogarithmicDispersion> {
        if self.amplitude == C64::ZERO {
            return None;
        }
        let first = self.scaled_derivatives[0] / self.amplitude;
        let second = self.scaled_derivatives[1] / self.amplitude - first * first;
        let first = divided(first, self.incident_wave_number);
        let second = divided(
            divided(second, self.incident_wave_number),
            self.incident_wave_number,
        );
        if !finite(first) || !finite(second) {
            return None;
        }
        Some(TransverseLogarithmicDispersion {
            lateral_shift_m: first.im,
            phase_curvature_m2: second.im,
            log_amplitude_derivative_m: first.re,
            log_amplitude_second_derivative_m2: second.re,
        })
    }
}

/// Source-bound coherent spatial response at one prepared frequency.
/// Reflection refers to the entrance; transmission includes the actual film
/// propagation to the exit face, including geometrical lateral displacement.
#[derive(Clone, Copy)]
pub struct LayeredTransverseDispersion<'s, 'm> {
    source: &'s PlanarDielectricLayers<'m>,
    transverse: f64,
    cosine: f64,
    polarization: Polarization,
    carrier: LayeredReflection,
    reflected: TransverseAmplitudeDispersion,
    transmitted: Option<[TransverseAmplitudeDispersion; 2]>,
}

impl<'s, 'm> LayeredTransverseDispersion<'s, 'm> {
    pub const fn source(self) -> &'s PlanarDielectricLayers<'m> {
        self.source
    }
    pub const fn transverse_wavenumber_rad_m(self) -> f64 {
        self.transverse
    }
    pub const fn incidence_cosine(self) -> f64 {
        self.cosine
    }
    pub const fn polarization(self) -> Polarization {
        self.polarization
    }
    pub const fn carrier(self) -> LayeredReflection {
        self.carrier
    }
    pub const fn reflected_tangential_e(self) -> TransverseAmplitudeDispersion {
        self.reflected
    }
    /// None after propagation underflow; a rounded zero does not prove a
    /// vanishing first or second spatial derivative.
    pub fn transmitted_tangential_e(self) -> Option<TransverseAmplitudeDispersion> {
        self.transmitted.map(|value| value[0])
    }
    pub fn transmitted_conjugate_h(self) -> Option<TransverseAmplitudeDispersion> {
        self.transmitted.map(|value| value[1])
    }
}

impl<'m> PlanarDielectricLayers<'m> {
    /// Differentiate the original coherent fields with respect to signed
    /// k_parallel at fixed frequency and material state. This plane contains
    /// the normal and one fixed tangential axis; k_parallel may be negative.
    /// Every supported monochromatic constitutive law remains supported.
    ///
    /// |k_parallel| must lie strictly inside the incident light cone. The
    /// original compensated angular-spectrum admission owns this check.
    /// Exact critical exits refuse a singular exterior derivative, while
    /// internal critical films retain their entire-function limit. Work is
    /// O(layers), allocation-free, with the original 64-film bound. There is
    /// no angular differencing, paraxial propagation or extra material fit.
    pub fn transverse_dispersion(
        &self,
        transverse_wavenumber_rad_m: f64,
        polarization: Polarization,
    ) -> Result<LayeredTransverseDispersion<'_, 'm>> {
        let cosine = incident_cosine(self.wave_number, transverse_wavenumber_rad_m, 0.0)
            .ok_or_else(|| bad("transverse carrier must lie inside the incident light cone"))?;
        let carrier = self.response(cosine, polarization)?;
        let p = transverse_wavenumber_rad_m / self.wave_number;
        let qi = Jet::new(
            C64::from_re(cosine),
            [
                C64::from_re(-p / cosine),
                C64::from_re(-((1.0 / cosine) / cosine) / cosine),
            ],
        );
        let normal_derivatives = [C64::from_re(-2.0 * p), C64::from_re(-2.0)];
        if !qi.is_finite() {
            return Err(bad("incident transverse jet is unrepresentable"));
        }
        let mut state = admittance(self.exit_wave, cosine, polarization, qi, normal_derivatives)?;
        let sum = state[0] + state[1];
        let mut transmitted = Some([state[0].scale(2.0) / sum, state[1].scale(2.0) / sum]);
        for layer in self.layers.iter().rev() {
            if layer.source.thickness_m == 0.0 {
                continue;
            }
            let slab = Slab::new(
                layer.wave,
                self.wave_number,
                cosine,
                polarization,
                layer.source.thickness_m,
            )?;
            if slab.scale == C64::ZERO {
                state = admittance(layer.wave, cosine, polarization, qi, normal_derivatives)?;
                transmitted = None;
                continue;
            }
            let q_squared = Jet::new(layer.wave.normal_squared(cosine), normal_derivatives);
            let a = match polarization {
                Polarization::Te => qi.scale(layer.wave.mu_ratio),
                Polarization::Tm => Jet::constant(layer.wave.epsilon_ratio) * qi,
            };
            let depth = Jet::constant(C64::from_re(self.wave_number * layer.source.thickness_m));
            let matrix = slab_jets_from_normal(slab, q_squared, a, depth, polarization)?;
            (state, transmitted) = propagate_jets(slab, matrix, state, transmitted)?;
        }
        let reflected = TransverseAmplitudeDispersion::from_jet(
            carrier.reflected_tangential_e,
            self.wave_number,
            reflection_jet(state),
        )?;
        let transmitted = transmitted
            .filter(|_| {
                carrier.transmitted_tangential_e != C64::ZERO
                    && carrier.transmitted_conjugate_h != C64::ZERO
            })
            .map(|jets| -> Result<_> {
                Ok([
                    TransverseAmplitudeDispersion::from_jet(
                        carrier.transmitted_tangential_e,
                        self.wave_number,
                        jets[0],
                    )?,
                    TransverseAmplitudeDispersion::from_jet(
                        carrier.transmitted_conjugate_h,
                        self.wave_number,
                        jets[1],
                    )?,
                ])
            })
            .transpose()?;
        Ok(LayeredTransverseDispersion {
            source: self,
            transverse: transverse_wavenumber_rad_m,
            cosine,
            polarization,
            carrier,
            reflected,
            transmitted,
        })
    }
}

fn admittance(
    medium: WaveMedium,
    cosine: f64,
    polarization: Polarization,
    qi: Jet,
    normal_derivatives: [C64; 2],
) -> Result<[Jet; 2]> {
    let normal = medium.normal(cosine)?;
    if normal == C64::ZERO {
        return Err(bad("critical exit has no finite transverse jet"));
    }
    let q = Jet::new(medium.normal_squared(cosine), normal_derivatives).root_at(normal);
    let (a, b) = match polarization {
        Polarization::Te => (qi.scale(medium.mu_ratio), q),
        Polarization::Tm => (q, Jet::constant(medium.epsilon_ratio) * qi),
    };
    let normalization = positive(maximum(a.value).max(maximum(b.value)))?;
    let state = [a.divide_real(normalization), b.divide_real(normalization)];
    if state.into_iter().any(|value| !value.is_finite()) {
        return Err(bad("exterior transverse jet is unrepresentable"));
    }
    Ok(state)
}
