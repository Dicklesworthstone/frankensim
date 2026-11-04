//! Coherent porous/fluid multilayers, reflection and transmission.
//!
//! Ordered homogeneous isotropic equivalent-fluid layers share ONE saturating
//! fluid snapshot. Pore drag, density, bulk modulus and source-frequency support
//! remain owned by PorousFluidLaw. Gaps use the same fluid's lossless acoustics;
//! explicit fluid layers retain their own original acoustic source. No elastic
//! frame, lateral finite-size effect or fluid attenuation is invented.
//!
//! exp(+i omega t), forward exp(-i k_n z). Tangential wavenumber is conserved,
//! k_n^2=k^2-k_0^2 sin^2(theta), Z_n=Z_c k/k_n. Pressure and normal velocity
//! are continuous between layers. RigidBackedLayers terminates in v_n=0;
//! TransmittingLayers retains an outgoing receiving-fluid half-space instead.
//! A bounded attenuating-wave transfer replaces growing cos/sin matrices:
//! E=exp(-2i k_n d); [p; Z_ref v] propagates through
//! [[1+E, (Z_n/Z_ref)(1-E)],[(Z_ref/Z_n)(1-E),1+E]], up to a common scalar.
//! Each step is rescaled, retaining impedance poles without dividing by p or v.
//!
//! Reference: COMSOL Porous Absorber, Eq. (5), rigid-frame oblique response:
//! <https://doc.comsol.com/6.1/doc/com.comsol.help.models.aco.porous_absorber/models.aco.porous_absorber.pdf>.
//! Distinct-fluid layers use the same pressure/normal-velocity conditions as
//! Oelze, ECE/TAM 373, Chapter 6, sections 6.3–6.4:
//! <https://jontalle.web.engr.illinois.edu/uploads/473.F18/Lectures/Chapter_6b.pdf>.
//! Plane reflection here is extended reaction. A local BEM boundary must
//! explicitly select the NORMAL-incidence admittance instead.

use super::{PorousFluidLaw, PorousFluidResponse};
use crate::MaterialError;
use crate::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::{c64::C64, det};

/// Coherent pressure and power transmission into a sourced receiving fluid.
pub mod transmission;
pub use transmission::{LayeredTransmission, TransmittingLayers};

/// Coherent refraction and P/S conversion into an original isotropic solid.
pub mod elastic;
pub use elastic::{ElasticTransmittingLayers, LayeredElasticTransmission};
/// Incoherent incident-power weighting of the original angular responses.
pub mod angular;
pub use angular::{
    AngularPower, AngularPowerBudget, AngularPowerError, DiffuseIncidence, IncidentPower,
    MAX_ANGULAR_POWER_SAMPLES,
};

/// Bounded preparation storage and per-angle work.
pub const MAX_LAYERS: usize = 64;

/// Ordered from the exposed face toward the termination. Zero thickness is
/// an exact identity, but an included porous source still must admit frequency.
#[derive(Clone, Copy)]
pub enum Layer<'a> {
    Porous {
        law: &'a PorousFluidLaw,
        thickness_m: f64,
    },
    /// Lossless gap filled by the original coherent saturating/exterior fluid.
    FluidGap { thickness_m: f64 },
    /// Lossless homogeneous fluid with its own density and adiabatic sound
    /// speed. The original source is retained; no transport loss is inferred.
    FluidLayer {
        fluid: &'a dyn FluidAcoustics,
        thickness_m: f64,
    },
}
impl core::fmt::Debug for Layer<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Porous { law, thickness_m } => f
                .debug_struct("Porous")
                .field("law", law)
                .field("thickness_m", thickness_m)
                .finish(),
            Self::FluidGap { thickness_m } => f
                .debug_struct("FluidGap")
                .field("thickness_m", thickness_m)
                .finish(),
            Self::FluidLayer { thickness_m, .. } => f
                .debug_struct("FluidLayer")
                .field("thickness_m", thickness_m)
                .finish_non_exhaustive(),
        }
    }
}
impl Layer<'_> {
    fn thickness_m(self) -> f64 {
        match self {
            Self::Porous { thickness_m, .. }
            | Self::FluidGap { thickness_m }
            | Self::FluidLayer { thickness_m, .. } => thickness_m,
        }
    }
}

#[derive(Clone, Copy)]
struct PreparedLayer<'a> {
    source: Layer<'a>,
    porous: Option<PorousFluidResponse>,
    fluid: Option<FluidAcousticState>,
}

/// Prepared material responses at ONE frequency. The incident/pore fluid has
/// one coherent state; each distinct explicit fluid source has its own snapshot.
/// Original sources remain borrowed, not replaced by fitted coefficients.
/// Each angle query is allocation-free.
pub struct RigidBackedLayers<'a> {
    source: &'a dyn ThermoviscousFluid,
    fluid: ThermoviscousFluidState,
    omega: f64,
    k0: f64,
    z0: f64,
    layers: Vec<PreparedLayer<'a>>,
}

/// Complex pressure reflection and normal absorbed/incident power fraction.
/// Values are not clipped; zero absorption may have floating-point residue
/// within the declared arithmetic balance check.
#[derive(Clone, Copy, Debug)]
pub struct LayeredReflection {
    pub pressure_reflection: C64,
    pub absorbed_fraction: f64,
    pressure: C64,
    normalized_velocity: C64,
    cosine: f64,
}
impl LayeredReflection {
    /// beta=Z0*v_n/p in exp(+i omega t). For a BEM owner with exp(-i omega t),
    /// conjugate ONCE. A pressure node has no representable finite admittance;
    /// reflection itself remains valid there. There is no panel-area factor.
    pub fn normalized_surface_admittance(&self) -> Result<C64, MaterialError> {
        if self.pressure == C64::ZERO {
            return Err(bad(
                "layered boundary has a pressure node, not a finite admittance",
            ));
        }
        let beta = self.normalized_velocity.scale(self.cosine) / self.pressure;
        if !finite(beta) {
            return Err(bad("layered surface admittance is unrepresentable"));
        }
        Ok(beta)
    }
}

fn bad(message: &'static str) -> MaterialError {
    MaterialError::State {
        what: message.into(),
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn magnitude_scale(z: C64) -> f64 {
    z.re.abs().max(z.im.abs())
}
fn divided(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}

// One source may occur as the incident fluid, several finite layers and the
// receiving half-space. Reuse its snapshot by original identity in every role.
fn retained_fluid_snapshot(
    original: &dyn FluidAcoustics,
    incident: &dyn FluidAcoustics,
    exterior: FluidAcousticState,
    layers: &[PreparedLayer<'_>],
) -> Option<FluidAcousticState> {
    if core::ptr::eq(original, incident) {
        Some(exterior)
    } else {
        layers.iter().find_map(|layer| match layer.source {
            Layer::FluidLayer { fluid, .. } if core::ptr::eq(original, fluid) => layer.fluid,
            _ => None,
        })
    }
}

impl<'a> RigidBackedLayers<'a> {
    pub fn new(
        fluid: &'a dyn ThermoviscousFluid,
        omega: f64,
        layers: &[Layer<'a>],
    ) -> Result<Self, MaterialError> {
        if layers.len() > MAX_LAYERS
            || !omega.is_finite()
            || omega <= 0.0
            || layers
                .iter()
                .any(|l| !l.thickness_m().is_finite() || l.thickness_m() < 0.0)
        {
            return Err(bad(
                "layered reflection needs <=64 layers, positive frequency and finite nonnegative thickness",
            ));
        }
        // The actual source is evaluated exactly once. Re-querying it per
        // layer could accidentally compose different evolving fluid states.
        let snapshot = fluid.thermoviscous_state()?;
        let acoustic = snapshot.acoustic();
        let k0 = omega / acoustic.sound_speed_m_s();
        let z0 = acoustic.characteristic_impedance_pa_s_m();
        if !k0.is_finite() || k0 <= 0.0 || !z0.is_finite() || z0 <= 0.0 {
            return Err(bad("layered exterior wave data are unrepresentable"));
        }
        let mut prepared: Vec<PreparedLayer<'a>> = Vec::with_capacity(layers.len());
        for &source in layers {
            let porous = match source {
                Layer::Porous { law, .. } => Some(law.effective_fluid(&snapshot, omega)?),
                Layer::FluidGap { .. } | Layer::FluidLayer { .. } => None,
            };
            let layer_fluid = if let Layer::FluidLayer {
                fluid: original, ..
            } = source
            {
                let incident: &dyn FluidAcoustics = fluid;
                let prior = retained_fluid_snapshot(original, incident, acoustic, &prepared);
                Some(match prior {
                    Some(state) => state,
                    None => original.acoustic_state()?,
                })
            } else {
                None
            };
            prepared.push(PreparedLayer {
                source,
                porous,
                fluid: layer_fluid,
            });
        }
        Ok(Self {
            source: fluid,
            fluid: snapshot,
            omega,
            k0,
            z0,
            layers: prepared,
        })
    }

    pub const fn fluid_source(&self) -> &'a dyn ThermoviscousFluid {
        self.source
    }
    pub const fn fluid(&self) -> &ThermoviscousFluidState {
        &self.fluid
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn exterior_wavenumber_rad_m(&self) -> f64 {
        self.k0
    }
    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    pub fn layer(&self, index: usize) -> Option<Layer<'a>> {
        self.layers.get(index).map(|l| l.source)
    }
    pub fn material_response(&self, index: usize) -> Option<&PorousFluidResponse> {
        self.layers.get(index).and_then(|l| l.porous.as_ref())
    }
    /// Coherent snapshot of an explicitly sourced finite fluid layer.
    /// Shared-fluid gaps and porous layers retain their existing response views.
    pub fn fluid_response(&self, index: usize) -> Option<&FluidAcousticState> {
        self.layers.get(index).and_then(|l| l.fluid.as_ref())
    }

    /// Exact planar extended-reaction reflection within the admitted material
    /// model. Positive incidence cosine is measured from the interface normal.
    /// Empty layers give the rigid wall. Source laws are never re-evaluated here.
    pub fn reflection(&self, cosine: f64) -> Result<LayeredReflection, MaterialError> {
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(bad("layered reflection needs incidence cosine in (0,1]"));
        }
        // Reference is Z0/cos(theta), the incident normal wave impedance.
        // Represent its ratio algebraically; forming Z0/cos could overflow.
        let (mut pressure, mut velocity) = (C64::ONE, C64::ZERO);
        for layer in self.layers.iter().rev() {
            let d = layer.source.thickness_m();
            if d == 0.0 {
                continue;
            }
            (pressure, velocity, _, _) = self.propagate_layer(layer, cosine, pressure, velocity)?;
        }
        let (reflection, absorbed) = reflected_state(pressure, velocity)?;
        Ok(LayeredReflection {
            pressure_reflection: reflection,
            absorbed_fraction: absorbed,
            pressure,
            normalized_velocity: velocity,
            cosine,
        })
    }

    // Both terminations consume the same normal-wave and critical-angle law.
    fn propagate_layer(
        &self,
        layer: &PreparedLayer<'_>,
        cosine: f64,
        pressure: C64,
        velocity: C64,
    ) -> Result<(C64, C64, f64, C64), MaterialError> {
        let d = layer.source.thickness_m();
        let (p, v, scale, kn) = if let Some(fluid) = layer.fluid {
            let speed_ratio = self.fluid.acoustic().sound_speed_m_s() / fluid.sound_speed_m_s();
            if !speed_ratio.is_finite() || speed_ratio <= 0.0 {
                return Err(bad("finite fluid sound-speed contrast is unrepresentable"));
            }
            let q = normal_wavenumber_ratio(
                speed_ratio,
                cosine,
                self.fluid.acoustic().sound_speed_m_s(),
                fluid.sound_speed_m_s(),
            );
            let kn = q.scale(self.k0);
            // Z_n/Z_ref = a/q, q=k_n/k0. Form the finite matrix entries
            // directly: at q=0, sin(k_n*d)/q -> k0*d. Dividing by k_n
            // would falsely reject critical refraction or overflow beside it.
            let a = (fluid.density_kg_m3() / self.fluid.acoustic().density_kg_m3()) * cosine;
            if !finite(kn) || (q != C64::ZERO && kn == C64::ZERO) || !a.is_finite() || a <= 0.0 {
                return Err(bad("finite fluid normal-wave contrast is unrepresentable"));
            }
            let (diagonal, upper, lower) = if q.im == 0.0 {
                let phase = kn.re * d;
                if !phase.is_finite() {
                    return Err(bad("finite fluid layer phase is unrepresentable"));
                }
                let (sine, cos) = det::sin_cos(phase);
                let sine_over_q = if phase.abs() < 1.0e-4 {
                    (self.k0 * d) * if phase == 0.0 { 1.0 } else { sine / phase }
                } else {
                    sine / q.re
                };
                (
                    C64::from_re(cos),
                    C64::new(0.0, a * sine_over_q),
                    C64::new(0.0, (q.re / a) * sine),
                )
            } else {
                let decay = -kn.im * d;
                let difference = -det::expm1(-2.0 * decay);
                let difference_over_q = if decay < 1.0e-4 {
                    (self.k0 * d)
                        * if decay == 0.0 {
                            2.0
                        } else {
                            difference / decay
                        }
                } else {
                    difference / -q.im
                };
                (
                    C64::from_re(2.0 - difference),
                    C64::new(0.0, a * difference_over_q),
                    C64::new(0.0, (q.im / a) * difference),
                )
            };
            let (p, v, scale) = normalized_step(
                pressure,
                velocity,
                diagonal.scale(1.0 / 16.0),
                upper.scale(1.0 / 16.0),
                lower.scale(1.0 / 16.0),
            )?;
            (p, v, scale, kn)
        } else {
            let (kn, ratio) = self.normal_wave(layer, cosine)?;
            let (sum, difference) = bounded_propagator(kn, d)?;
            let (p, v, scale) = transfer_step(pressure, velocity, sum, difference, ratio)?;
            (p, v, scale, kn)
        };
        Ok((p, v, scale, kn))
    }

    // Shared by rigid and transmitting terminations: no second material law.
    fn normal_wave(
        &self,
        layer: &PreparedLayer<'_>,
        cosine: f64,
    ) -> Result<(C64, C64), MaterialError> {
        let (kn, ratio) = match layer.porous {
            None => (C64::from_re(self.k0 * cosine), C64::ONE),
            Some(response) => {
                let (zc, k) = response.characteristic();
                let kn = if cosine == 1.0 {
                    k
                } else {
                    let tangential = self.k0 * det::sqrt((1.0 - cosine) * (1.0 + cosine));
                    let mut value = (k * k - C64::from_re(tangential * tangential)).sqrt();
                    if value.im > 0.0 {
                        value = -value;
                    }
                    value
                };
                let zn = zc * (k / kn);
                (kn, divided(zn, self.z0).scale(cosine))
            }
        };
        if !finite(kn)
            || kn.re < 0.0
            || kn.im > 0.0
            || kn == C64::ZERO
            || !finite(ratio)
            || ratio.re <= 0.0
        {
            return Err(bad(
                "layer normal wave impedance is unrepresentable or nonpassive",
            ));
        }
        Ok((kn, ratio))
    }
}

// q=k_n/k_incident for a lossless fluid. Share this Snell branch between
// finite fluid layers and the outgoing receiving half-space. Keep matched
// grazing and small normal contrasts; ratio^2-(1-cosine^2) loses those limits.
// The cached ratio alone cannot classify critical incidence: its rounding can
// change the sign. Retain the original sound speeds for contrast and refinement.
fn normal_wavenumber_ratio(
    ratio: f64,
    cosine: f64,
    incident_speed: f64,
    transmitted_speed: f64,
) -> C64 {
    if cosine == 1.0 {
        return C64::from_re(ratio);
    }
    if incident_speed == transmitted_speed {
        return C64::from_re(cosine);
    }
    if ratio > 2.0 {
        let sine_over_ratio = det::sqrt((1.0 - cosine) * (1.0 + cosine)) / ratio;
        return C64::from_re(ratio * det::sqrt(1.0 - sine_over_ratio * sine_over_ratio));
    }
    let (squared, cancellation_scale) = if cosine > 0.5 {
        let tangential_squared = (1.0 - cosine) * (1.0 + cosine);
        (
            ratio * ratio - tangential_squared,
            ratio * ratio + tangential_squared,
        )
    } else {
        let contrast = ((incident_speed - transmitted_speed) / transmitted_speed) * (ratio + 1.0);
        (
            cosine.mul_add(cosine, contrast),
            cosine * cosine + contrast.abs(),
        )
    };
    let squared = if incident_speed < transmitted_speed
        && squared.abs() <= 64.0 * f64::EPSILON * cancellation_scale
    {
        use fs_math::dd::Dd;
        // Power-of-two scaling preserves the original speeds exactly and keeps
        // the precision ladder away from squared-speed overflow/underflow.
        let scale = f64::from_bits(transmitted_speed.to_bits() & 0x7ff0_0000_0000_0000)
            .max(f64::MIN_POSITIVE);
        let incident = Dd::from_f64(incident_speed / scale);
        let transmitted = Dd::from_f64(transmitted_speed / scale);
        let normal = transmitted * Dd::from_f64(cosine);
        (((incident - transmitted) * (incident + transmitted) + normal * normal)
            / (transmitted * transmitted))
            .to_f64()
    } else {
        squared
    };
    if squared >= 0.0 {
        C64::from_re(det::sqrt(squared))
    } else {
        C64::new(0.0, -det::sqrt(-squared))
    }
}

// Return the normalization as well, so a transmitting owner can retain the
// common pressure amplitude that cancels out of rigid reflection.
fn transfer_step(
    pressure: C64,
    velocity: C64,
    sum: C64,
    difference: C64,
    ratio: C64,
) -> Result<(C64, C64, f64), MaterialError> {
    let inverse = ratio.recip();
    if !finite(inverse) {
        return Err(bad("layer inverse wave impedance is unrepresentable"));
    }
    // A fixed common power-of-two factor leaves headroom for complex
    // products. Contrast-dependent scaling would erase O(cosine^2)
    // entries before state normalization at extreme grazing angles.
    let diagonal = sum.scale(1.0 / 16.0);
    let off_diagonal = difference.scale(1.0 / 16.0);
    let upper = ratio * off_diagonal;
    let lower = inverse * off_diagonal;
    normalized_step(pressure, velocity, diagonal, upper, lower)
}

fn normalized_step(
    pressure: C64,
    velocity: C64,
    diagonal: C64,
    upper: C64,
    lower: C64,
) -> Result<(C64, C64, f64), MaterialError> {
    let next_p = diagonal * pressure + upper * velocity;
    let next_v = lower * pressure + diagonal * velocity;
    let normalization = magnitude_scale(next_p).max(magnitude_scale(next_v));
    if !finite(next_p) || !finite(next_v) || normalization == 0.0 {
        return Err(bad("layer transfer state is unrepresentable"));
    }
    Ok((
        divided(next_p, normalization),
        divided(next_v, normalization),
        normalization,
    ))
}

fn reflected_state(pressure: C64, velocity: C64) -> Result<(C64, f64), MaterialError> {
    let denominator = pressure + velocity;
    let norm = denominator.norm_sq();
    if !norm.is_finite() || norm <= 0.0 {
        return Err(bad("layered incident pressure is unrepresentable"));
    }
    // Preserve the small complex correction when R approaches +/-1;
    // direct subtraction and division would lose grazing reflection phase.
    let reflection = if magnitude_scale(pressure) < 0.25 * magnitude_scale(velocity) {
        -C64::ONE + (pressure / denominator).scale(2.0)
    } else if magnitude_scale(velocity) < 0.25 * magnitude_scale(pressure) {
        C64::ONE - (velocity / denominator).scale(2.0)
    } else {
        (pressure - velocity) / denominator
    };
    // Evaluate actual surface work. 1-|R|^2 would erase weak absorption
    // whenever the reflected magnitude rounds to unity.
    let entered = 4.0 * (pressure.re * velocity.re + pressure.im * velocity.im) / norm;
    let tolerance = 128.0 * f64::EPSILON;
    if !finite(reflection)
        || !entered.is_finite()
        || entered < -tolerance
        || (reflection.norm_sq() + entered - 1.0).abs() > tolerance
    {
        return Err(bad(
            "layered reflection failed passive surface-work balance",
        ));
    }
    Ok((reflection, entered))
}

/// Bounded common-scaled cos/sin pair. The lossless branch is purely real /
/// imaginary, so a stack of fluid gaps retains exactly zero dissipative work.
fn bounded_propagator(kn: C64, thickness: f64) -> Result<(C64, C64), MaterialError> {
    let decay = 2.0 * (kn.im * thickness);
    let attenuation = if decay == f64::NEG_INFINITY {
        0.0
    } else {
        det::exp(decay)
    };
    if attenuation == 0.0 {
        return Ok((C64::ONE, C64::ONE));
    }
    let phase = kn.re * thickness;
    if !phase.is_finite() {
        return Err(bad("layer phase is unrepresentable"));
    }
    let (sine, cosine) = det::sin_cos(phase);
    if kn.im == 0.0 {
        return Ok((C64::from_re(cosine), C64::new(0.0, sine)));
    }
    let loss = -det::expm1(decay);
    let cross = 2.0 * attenuation * sine * cosine;
    // Stable 1-E and 1+E even at thin layers and near quarter-wave poles.
    Ok((
        C64::new(loss + 2.0 * attenuation * cosine * cosine, -cross),
        C64::new(loss + 2.0 * attenuation * sine * sine, cross),
    ))
}

#[cfg(test)]
mod normal_wave_tests {
    use super::*;

    fn normal(incident: f64, transmitted: f64, cosine: f64) -> C64 {
        normal_wavenumber_ratio(incident / transmitted, cosine, incident, transmitted)
    }

    #[test]
    fn g1_critical_side_uses_original_speeds_and_exact_input_oracles() {
        // Independent rational arithmetic on the exact binary64 inputs,
        // followed by an 85-decimal-digit square root. No rounded speed ratio
        // or floating sine/cosine reconstruction participates in these oracles.
        for (incident, transmitted, cosine, expected) in [
            (
                340.0,
                500.0,
                0.733_212_111_192_934_4,
                C64::new(0.0, -3.952_211_632_119_969_6e-9),
            ),
            (340.0, 425.0, 0.6, C64::new(0.0, -5.161_913_655_903_57e-9)),
            (
                340.0,
                425.0,
                0.599_999_999_999_999_9,
                C64::new(0.0, -1.264_405_455_326_820_9e-8),
            ),
            (
                1.0,
                1.0_f64.next_up(),
                2.107_342_425_544_701_4e-8,
                C64::from_re(2.630_083_138_697_363_2e-16),
            ),
            (
                340.0,
                340.0_f64.next_up(),
                2.107_342_425_544_701_4e-8,
                C64::from_re(1.047_454_809_181_105_1e-8),
            ),
        ] {
            let actual = normal(incident, transmitted, cosine);
            assert!(
                (actual - expected).abs() <= 8e-15 * expected.abs(),
                "{actual:?} != {expected:?}"
            );
            assert_eq!(actual.re == 0.0, expected.re == 0.0);
            assert_eq!(actual.im == 0.0, expected.im == 0.0);
        }
        let cosine = 0.733_212_111_192_934_4_f64;
        assert!(normal(340.0, 500.0, cosine.next_down()).im < 0.0);
        assert!(normal(340.0, 500.0, cosine.next_up()).re > 0.0);
    }

    #[test]
    fn g3_normal_wave_preserves_speed_scaling_and_exact_fast_paths() {
        let cosine = 0.733_212_111_192_934_4;
        let expected = normal(340.0, 500.0, cosine);
        for exponent in [-1000, -500, 0, 500, 1000] {
            let scale = 2.0_f64.powi(exponent);
            assert_eq!(normal(340.0 * scale, 500.0 * scale, cosine), expected);
        }
        for speed in [1e-150, 340.0, 1e150] {
            for cosine in [1.0, 0.25, 1e-200] {
                assert_eq!(normal(speed, speed, cosine), C64::from_re(cosine));
            }
        }
        for (incident, transmitted) in [(340.0, 500.0), (500.0, 340.0), (1.0, 1e-150)] {
            assert_eq!(
                normal(incident, transmitted, 1.0),
                C64::from_re(incident / transmitted)
            );
        }
    }
}
