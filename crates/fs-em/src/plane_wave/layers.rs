//! Coherent monochromatic reflection/transmission through finite isotropic films.
//!
//! Layers are ordered from the incident face to the exit half-space. Sources,
//! SI epsilon/mu and exp(+i omega t) phasors are the same as the parent owner.
//! Incoming light is lossless; each film and the exit may use its original Debye
//! law. Transmitted power enters the exit face; finite-film absorption excludes
//! subsequent dissipation in that semi-infinite medium.
//!
//! Composition retains bounded homogeneous E/H states in common real reference
//! ports and reconstructs forward scattering gains. It never subtracts a product
//! of nearly unit reflections at a film resonance. Only exp(-i k_z d), whose
//! magnitude is <=1, propagates. The analytic k_z=0 limit retains critical
//! films. No growing evanescent transfer matrix is formed.
//! Absorption uses dissipative even/odd slab excitations, avoiding subtraction
//! of nearly equal reflected and transmitted powers.
//!
//! Plane-wave continuity and Poynting conventions: S. J. Byrnes, "Multilayer
//! optical calculations", sections 4.1, 4.3-4.5, <https://arxiv.org/abs/1603.02720>.
//! That reference uses the opposite phasor and rotating optical p basis; this
//! owner retains the parent's fixed tangential-E axis and explicit sourced mu.
//! Finite coherence length, roughness, anisotropy and a mesh Maxwell solve are
//! outside this exact homogeneous planar model.

use super::{
    DielectricResponse, Polarization, ReflectionError, WaveMedium, finite, interface_response,
    positive,
};
use fs_material::{
    field_storage::{linear_electric_permittivity, linear_magnetic_permeability},
    state_point::ResolvedMaterialStatePoint,
};
use fs_math::{c64::C64, det};

/// Preparation storage and allocation-free angle work are bounded by this cap.
pub const MAX_LAYERS: usize = 64;

/// A finite film, ordered from the incident face toward the exit half-space.
#[derive(Clone, Copy, Debug)]
pub struct DielectricLayer<'a> {
    pub response: DielectricResponse<'a>,
    pub thickness_m: f64,
}

#[derive(Clone, Copy)]
struct PreparedLayer<'a> {
    source: DielectricLayer<'a>,
    wave: WaveMedium,
}

/// Frequency-resolved original material laws; angle queries allocate nothing.
pub struct PlanarDielectricLayers<'a> {
    incident: &'a ResolvedMaterialStatePoint,
    exit: DielectricResponse<'a>,
    exit_wave: WaveMedium,
    omega: f64,
    wave_number: f64,
    layers: Vec<PreparedLayer<'a>>,
}

/// Tangential-E amplitudes and normal Poynting fractions for one-sided incidence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayeredReflection {
    /// Referenced to the incident face, with the parent's fixed tangential axis.
    pub reflected_tangential_e: C64,
    /// Referenced to the last film's exit face, including coherent film phase.
    pub transmitted_tangential_e: C64,
    pub reflected_power_fraction: f64,
    /// Flux entering the exit half-space, before any subsequent attenuation.
    pub transmitted_power_fraction: f64,
    /// Physical loss in the finite films only; zero for wholly lossless films.
    pub absorbed_power_fraction: f64,
    pub incident_normal_wavenumber_rad_m: f64,
    pub exit_normal_wavenumber_rad_m: C64,
    pub exit_permittivity_f_m: C64,
}

impl<'a> PlanarDielectricLayers<'a> {
    /// Admit and cache every original law at one positive frequency.
    /// Zero-thickness films are identities but still require valid sources.
    pub fn try_new(
        incident: &'a ResolvedMaterialStatePoint,
        exit: DielectricResponse<'a>,
        omega_rad_s: f64,
        layers: &[DielectricLayer<'a>],
    ) -> Result<Self, ReflectionError> {
        if layers.len() > MAX_LAYERS
            || layers
                .iter()
                .any(|l| !l.thickness_m.is_finite() || l.thickness_m < 0.0)
        {
            return Err(bad(
                "finite layers require <=64 finite nonnegative thicknesses",
            ));
        }
        let omega = positive(omega_rad_s)?;
        let epsilon_i = linear_electric_permittivity(incident)?;
        let mu_i = linear_magnetic_permeability(incident)?;
        let wave_number = positive(omega * det::sqrt(epsilon_i) * det::sqrt(mu_i))?;
        let prepare = |source: DielectricResponse<'a>| {
            WaveMedium::new(
                source.permittivity(omega)?,
                linear_magnetic_permeability(source.material())?,
                epsilon_i,
                mu_i,
            )
        };
        let exit_wave = prepare(exit)?;
        let mut prepared = Vec::with_capacity(layers.len());
        for &source in layers {
            prepared.push(PreparedLayer {
                source,
                wave: prepare(source.response)?,
            });
        }
        Ok(Self {
            incident,
            exit,
            exit_wave,
            omega,
            wave_number,
            layers: prepared,
        })
    }

    pub const fn incident_material(&self) -> &'a ResolvedMaterialStatePoint {
        self.incident
    }
    pub fn exit_material(&self) -> &'a ResolvedMaterialStatePoint {
        self.exit.material()
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    pub fn layer(&self, index: usize) -> Option<DielectricLayer<'a>> {
        self.layers.get(index).map(|layer| layer.source)
    }
    pub fn layer_permittivity_f_m(&self, index: usize) -> Option<C64> {
        self.layers.get(index).map(|layer| layer.wave.epsilon_f_m)
    }

    /// Coherent TE/TM scattering, including evanescent tunnelling and film loss.
    /// Exact grazing has zero incoming normal power and is explicitly refused.
    pub fn response(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
    ) -> Result<LayeredReflection, ReflectionError> {
        let c = positive(incidence_cosine)?;
        if c > 1.0 {
            return Err(bad("incidence cosine must lie in (0,1]"));
        }
        let exit = interface_response(self.exit_wave, self.wave_number, c, polarization)?;
        let mut slabs = [Slab::IDENTITY; MAX_LAYERS];
        let mut states = [[C64::ZERO; 2]; MAX_LAYERS + 1];
        let mut gains = [C64::ZERO; MAX_LAYERS];
        let n = self.layers.len();
        for (i, layer) in self.layers.iter().enumerate() {
            slabs[i] = Slab::new(
                layer.wave,
                self.wave_number,
                c,
                polarization,
                layer.source.thickness_m,
            )?;
        }
        // Homogeneous [E, H/Y_inc] retains both 1+r and 1-r without
        // subtracting a reflection already rounded to +/-1. Initialize from
        // the original admittance pair, including an exactly critical exit.
        states[n] = admittance_state(self.exit_wave, c, polarization)?;
        let mut transmitted = exit.transmitted_tangential_e;
        let mut transmitted_power = exit.transmitted_power_fraction;
        for i in (0..n).rev() {
            let slab = slabs[i];
            if self.layers[i].source.thickness_m == 0.0 {
                gains[i] = C64::ONE;
                states[i] = states[i + 1];
                continue;
            }
            if slab.scale == C64::ZERO {
                // Attenuation underflow gives the opaque front-interface
                // limit. Its original admittance pair is retained directly.
                states[i] = [slab.matrix[0][0], slab.matrix[1][0]];
                transmitted = C64::ZERO;
                transmitted_power = 0.0;
            } else {
                let right = states[i + 1];
                let electric = slab.matrix[0][0] * right[0] + slab.matrix[0][1] * right[1];
                let magnetic = slab.matrix[1][0] * right[0] + slab.matrix[1][1] * right[1];
                let normalization = positive(maximum(electric).max(maximum(magnetic)))?;
                states[i] = [
                    divided(electric, normalization),
                    divided(magnetic, normalization),
                ];
                // For a passive tail E and H/Y_inc have nonnegative real
                // cross-power: their sum avoids the cancellation of
                // 1-r_left*r_right at a high-reflectivity film resonance.
                let incident = states[i][0] + states[i][1];
                if incident == C64::ZERO || !finite(incident) {
                    return Err(bad("film incident field is unrepresentable"));
                }
                let gain = divided(slab.scale, normalization) * ((right[0] + right[1]) / incident);
                gains[i] = gain;
                transmitted = gain * transmitted;
                transmitted_power = gain.scale(det::sqrt(transmitted_power)).norm_sq();
                if !finite(gain) || !finite(transmitted) || !transmitted_power.is_finite() {
                    return Err(bad("composed film scattering is unrepresentable"));
                }
            }
        }
        let reflected = (states[0][0] - states[0][1]) / (states[0][0] + states[0][1]);
        let mut incoming = C64::ONE;
        let mut absorbed = 0.0;
        for i in 0..n {
            let outgoing = incoming * gains[i];
            let right = states[i + 1];
            let tail_reflection = (right[0] - right[1]) / (right[0] + right[1]);
            let opposite = tail_reflection * outgoing;
            let even = (incoming + opposite).scale(det::sqrt(0.5 * slabs[i].even_loss));
            let odd = (incoming - opposite).scale(det::sqrt(0.5 * slabs[i].odd_loss));
            absorbed += even.norm_sq() + odd.norm_sq();
            incoming = outgoing;
        }
        let reflected_power = reflected.norm_sq();
        let tolerance = 256.0 * f64::EPSILON * (n + 1) as f64;
        if !finite(reflected)
            || !absorbed.is_finite()
            || absorbed < 0.0
            || !transmitted_power.is_finite()
            || transmitted_power < 0.0
            || (reflected_power + transmitted_power + absorbed - 1.0).abs() > tolerance
        {
            return Err(bad("finite-film Poynting balance failed"));
        }
        Ok(LayeredReflection {
            reflected_tangential_e: reflected,
            transmitted_tangential_e: transmitted,
            reflected_power_fraction: reflected_power,
            transmitted_power_fraction: transmitted_power,
            absorbed_power_fraction: absorbed,
            incident_normal_wavenumber_rad_m: exit.incident_normal_wavenumber_rad_m,
            exit_normal_wavenumber_rad_m: exit.transmitted_normal_wavenumber_rad_m,
            exit_permittivity_f_m: self.exit_wave.epsilon_f_m,
        })
    }
}

fn bad(message: &'static str) -> ReflectionError {
    ReflectionError::Invalid(message)
}
fn maximum(z: C64) -> f64 {
    z.re.abs().max(z.im.abs())
}
fn divided(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}

fn admittance_state(
    medium: WaveMedium,
    cosine: f64,
    polarization: Polarization,
) -> Result<[C64; 2], ReflectionError> {
    let normal = medium.normal(cosine)?;
    let (a, b) = match polarization {
        Polarization::Te => (C64::from_re(medium.mu_ratio * cosine), normal),
        Polarization::Tm => (normal, medium.epsilon_ratio.scale(cosine)),
    };
    let scale = positive(maximum(a).max(maximum(b)))?;
    Ok([divided(a, scale), divided(b, scale)])
}

#[derive(Clone, Copy)]
struct Slab {
    // Physical left state = matrix * right state / scale. Both are bounded;
    // scale=0 explicitly denotes the opaque front-interface limit.
    matrix: [[C64; 2]; 2],
    scale: C64,
    even_loss: f64,
    odd_loss: f64,
}
impl Slab {
    const IDENTITY: Self = Self {
        matrix: [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]],
        scale: C64::ONE,
        even_loss: 0.0,
        odd_loss: 0.0,
    };

    fn new(
        medium: WaveMedium,
        wave_number: f64,
        cosine: f64,
        polarization: Polarization,
        thickness: f64,
    ) -> Result<Self, ReflectionError> {
        if thickness == 0.0 {
            return Ok(Self::IDENTITY);
        }
        let q = medium.normal(cosine)?;
        let depth = positive(wave_number * thickness)?;
        let decay = q.im * depth;
        let attenuation = if decay == f64::NEG_INFINITY {
            0.0
        } else {
            det::exp(decay)
        };
        if attenuation == 0.0 {
            // An opaque slab is its front interface plus absorption: it
            // needs neither a growing reciprocal exponential nor its phase.
            let interface = interface_response(medium, wave_number, cosine, polarization)?;
            let state = admittance_state(medium, cosine, polarization)?;
            return Ok(Self {
                matrix: [[state[0], C64::ZERO], [state[1], C64::ZERO]],
                scale: C64::ZERO,
                even_loss: interface.transmitted_power_fraction,
                odd_loss: interface.transmitted_power_fraction,
            });
        }
        let phase = q.re * depth;
        if !phase.is_finite() {
            return Err(bad("film propagation phase is unrepresentable"));
        }
        let (sine, cos_phase) = det::sin_cos(phase);
        let p = C64::new(attenuation * cos_phase, -attenuation * sine);
        let power_decay = det::exp(2.0 * decay);
        let loss = -det::expm1(2.0 * decay);
        let difference = C64::new(
            loss + 2.0 * power_decay * sine * sine,
            2.0 * power_decay * sine * cos_phase,
        );
        let s = match polarization {
            Polarization::Te => C64::from_re(medium.mu_ratio * cosine),
            Polarization::Tm => medium.epsilon_ratio.scale(cosine),
        };
        let scale = positive(maximum(s).max(maximum(q)))?;
        let s = divided(s, scale);
        let normal = divided(q, scale);
        let s_squared = s * s;
        let q_squared = normal * normal;
        if s == C64::ZERO
            || s_squared == C64::ZERO
            || (normal != C64::ZERO && q_squared == C64::ZERO)
        {
            return Err(bad("film coupling is below representable precision"));
        }
        // Cancel the normal-wave factor before the critical limit. Choose
        // u=(1-p^2)/(q/h) or its reciprocal, whichever remains bounded.
        let (coefficient, reciprocal) = if q == C64::ZERO {
            let critical = depth * scale;
            if critical <= 0.5 {
                (C64::new(0.0, 2.0 * critical), false)
            } else {
                (C64::new(0.0, -(0.5 / depth) / scale), true)
            }
        } else if maximum(difference) > maximum(normal) {
            (normal / difference, true)
        } else {
            if normal == C64::ZERO {
                return Err(bad("film normal-wave contrast is unrepresentable"));
            }
            (difference / normal, false)
        };
        let phase_sum = C64::new(
            loss + 2.0 * power_decay * cos_phase * cos_phase,
            -2.0 * power_decay * sine * cos_phase,
        );
        let (diagonal, s_term, q_term, transfer_scale) = if reciprocal {
            (
                s * phase_sum * coefficient,
                s_squared,
                q_squared,
                s.scale(2.0) * p * coefficient,
            )
        } else {
            (
                s * phase_sum,
                s_squared * coefficient,
                q_squared * coefficient,
                s.scale(2.0) * p,
            )
        };
        let (upper, lower) = match polarization {
            Polarization::Te => (s_term, q_term),
            Polarization::Tm => (q_term, s_term),
        };
        // Common fixed scaling keeps every matrix/state product bounded.
        // The same factor MUST accompany the retained transmission scale.
        let matrix = [
            [diagonal.scale(1.0 / 16.0), upper.scale(1.0 / 16.0)],
            [lower.scale(1.0 / 16.0), diagonal.scale(1.0 / 16.0)],
        ];
        let transfer_scale = transfer_scale.scale(1.0 / 16.0);
        let (a, b) = match polarization {
            Polarization::Te => (s, normal),
            Polarization::Tm => (normal, s),
        };
        let (even_loss, odd_loss) = if medium.epsilon_f_m.im == 0.0 {
            (0.0, 0.0)
        } else {
            let cross = a * b.conj();
            let (even_work, odd_work) = if phase.abs().max(decay.abs()) < 0.1 {
                // Positive small-depth dissipation, including the O(d^3)
                // odd mode. Direct subtraction would erase very thin loss.
                let x = phase * phase;
                let y = decay * decay;
                let sinc_defect =
                    x * (1.0 / 6.0 - x * (1.0 / 120.0 - x * (1.0 / 5040.0 - x / 362880.0)));
                let sinhc_excess =
                    y * (1.0 / 6.0 + y * (1.0 / 120.0 + y * (1.0 / 5040.0 + y / 362880.0)));
                let plus = 2.0 + sinhc_excess - sinc_defect;
                let minus = sinhc_excess + sinc_defect;
                match polarization {
                    Polarization::Te => {
                        let factor =
                            attenuation * a.re * (-medium.index_squared_ratio.im / scale) * depth;
                        (factor * plus, factor * minus)
                    }
                    Polarization::Tm => {
                        let factor =
                            attenuation * ((-medium.epsilon_ratio.im * cosine) / scale) * depth;
                        let normal_weight = q.re * normal.re + q.im * normal.im;
                        let lateral_weight = (1.0 - cosine) * (1.0 + cosine) / scale;
                        (
                            factor * (normal_weight * plus + lateral_weight * minus),
                            factor * (normal_weight * minus + lateral_weight * plus),
                        )
                    }
                }
            } else {
                (
                    cross.re * loss - 2.0 * cross.im * p.im,
                    cross.re * loss + 2.0 * cross.im * p.im,
                )
            };
            let sum = C64::ONE + p;
            // Stable 1-p, including weak attenuation and small phase.
            let half_sine = det::sin(0.5 * phase);
            let difference = C64::new(
                -det::expm1(decay) + 2.0 * attenuation * half_sine * half_sine,
                attenuation * sine,
            );
            let even_den = positive((a * sum + b * difference).norm_sq())?;
            let odd_den = positive((a * difference + b * sum).norm_sq())?;
            (4.0 * even_work / even_den, 4.0 * odd_work / odd_den)
        };
        if matrix.iter().flatten().any(|&value| !finite(value))
            || !finite(transfer_scale)
            || transfer_scale == C64::ZERO
            || !even_loss.is_finite()
            || !odd_loss.is_finite()
            || even_loss < 0.0
            || odd_loss < 0.0
            || even_loss > 1.0 + 256.0 * f64::EPSILON
            || odd_loss > 1.0 + 256.0 * f64::EPSILON
        {
            return Err(bad(
                "finite film scattering or dissipative work is unrepresentable",
            ));
        }
        Ok(Self {
            matrix,
            scale: transfer_scale,
            even_loss,
            odd_loss,
        })
    }
}
