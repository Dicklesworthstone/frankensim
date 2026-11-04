//! Coherent transmission through the existing prepared equivalent-fluid layers.
//!
//! The incident/saturating state and every porous law retain their original
//! owners. The receiving source supplies density and sound speed for an explicit
//! LOSSLESS homogeneous half-space. No attenuation is inferred from that trait.
//!
//! exp(+i omega t), outgoing exp(-i k_n z): pressure and normal velocity are
//! continuous, with k_n^2=(omega/c_receiver)^2-k_t^2. The outgoing branch has
//! nonnegative real k_n, or negative imaginary k_n under total internal reflection.
//! Pressure T is referenced to the exit face; it is not sqrt(transmitted power).
//! Power is weighted by the receiving/incident normal admittance.
//!
//! Reference: M. Oelze, ECE/TAM 373, Chapter 6, fluid interfaces and fluid layers:
//! <https://jontalle.web.engr.illinois.edu/uploads/473.F18/Lectures/Chapter_6a.pdf>.
//! The rigid-frame, infinite-planar and shared-fluid limits of the parent remain.

use super::{
    Layer, RigidBackedLayers, bad, divided, finite, magnitude_scale, normal_wavenumber_ratio,
    reflected_state, retained_fluid_snapshot,
};
use crate::MaterialError;
use crate::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::{c64::C64, det};

/// Prepared porous/gap layers terminated by an outgoing receiving fluid.
/// Both original fluid sources remain borrowed. Each is evaluated once during
/// preparation; angle queries allocate nothing and re-query no material source.
pub struct TransmittingLayers<'a> {
    prepared: RigidBackedLayers<'a>,
    receiver_source: &'a dyn FluidAcoustics,
    receiver: FluidAcousticState,
    speed_ratio: f64,
    density_ratio: f64,
    lossless: bool,
}

/// Pressure amplitudes at the entrance/exit planes and normal power fractions.
/// An evanescent receiver retains complex T but carries zero outgoing normal
/// power. Absorption is entrance surface work minus outgoing power; for lossy
/// stacks, values below floating-point subtraction resolution are not resolved.
/// No coefficient is clipped to force a balance.
#[derive(Clone, Copy, Debug)]
pub struct LayeredTransmission {
    pub pressure_reflection: C64,
    pub pressure_transmission: C64,
    pub receiving_normal_wavenumber_rad_m: C64,
    pub reflected_fraction: f64,
    pub transmitted_fraction: f64,
    pub absorbed_fraction: f64,
}

impl<'a> TransmittingLayers<'a> {
    /// The receiver is explicitly modeled as lossless using its original
    /// acoustic state. Gaps and pore fluid retain the incident fluid snapshot;
    /// explicit `FluidLayer` entries retain their own acoustic snapshots.
    /// A source reused in any of these roles shares its one prepared state.
    /// Material support and the parent's MAX_LAYERS limit are unchanged.
    pub fn new(
        incident: &'a dyn ThermoviscousFluid,
        receiver_source: &'a dyn FluidAcoustics,
        omega: f64,
        layers: &[Layer<'a>],
    ) -> Result<Self, MaterialError> {
        let prepared = RigidBackedLayers::new(incident, omega, layers)?;
        let exterior = prepared.fluid().acoustic();
        let original_incident: &dyn FluidAcoustics = incident;
        let receiver = match retained_fluid_snapshot(
            receiver_source,
            original_incident,
            exterior,
            &prepared.layers,
        ) {
            Some(state) => state,
            None => receiver_source.acoustic_state()?,
        };
        let speed_ratio = exterior.sound_speed_m_s() / receiver.sound_speed_m_s();
        let density_ratio = exterior.density_kg_m3() / receiver.density_kg_m3();
        if !speed_ratio.is_finite()
            || speed_ratio <= 0.0
            || !density_ratio.is_finite()
            || density_ratio <= 0.0
        {
            return Err(bad("receiving acoustic contrast is unrepresentable"));
        }
        let lossless = prepared
            .layers
            .iter()
            .all(|layer| layer.porous.is_none() || layer.source.thickness_m() == 0.0);
        Ok(Self {
            prepared,
            receiver_source,
            receiver,
            speed_ratio,
            density_ratio,
            lossless,
        })
    }

    /// The same original layer preparation, also available for a rigid-backed
    /// comparison without re-evaluating the incident fluid or any porous law.
    pub const fn prepared_layers(&self) -> &RigidBackedLayers<'a> {
        &self.prepared
    }
    pub const fn receiving_source(&self) -> &'a dyn FluidAcoustics {
        self.receiver_source
    }
    pub const fn receiving_state(&self) -> &FluidAcousticState {
        &self.receiver
    }

    /// Exact planar reflection/refraction/transmission within the admitted
    /// material model. Critical incidence (receiving k_n=0) remains finite.
    pub fn scattering(&self, cosine: f64) -> Result<LayeredTransmission, MaterialError> {
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(bad("layer transmission needs incidence cosine in (0,1]"));
        }
        let q = normal_wavenumber_ratio(
            self.speed_ratio,
            cosine,
            self.prepared.fluid().acoustic().sound_speed_m_s(),
            self.receiver.sound_speed_m_s(),
        );
        let receiving_kn = q.scale(self.prepared.k0);
        if !finite(receiving_kn) || (q != C64::ZERO && receiving_kn == C64::ZERO) {
            return Err(bad("receiving normal wavenumber is unrepresentable"));
        }
        // Y_receiver/Y_incident = (rho_incident/rho_receiver)*q/cos(theta).
        // Never divide by receiving k_n: it may be exactly zero at criticality.
        let load = if q == C64::ZERO {
            C64::ZERO
        } else {
            let value = positive_product_ratio(self.density_ratio, magnitude_scale(q), cosine)?;
            if q.re > 0.0 {
                C64::from_re(value)
            } else {
                C64::new(0.0, -value)
            }
        };
        let state = transmit_state(&self.prepared, cosine, C64::ONE, load)?;
        let reflection = state.reflection;
        let transmission = state.coefficient(C64::ONE)?;
        // Compute weighted power before underflowing |T|^2. A small pressure
        // amplitude can still carry representable power into a high admittance.
        let transmitted = if load.re == 0.0 {
            0.0
        } else {
            state.power(C64::ONE, det::ln(load.re))
        };
        let reflected = reflection.norm_sq();
        let absorbed = if self.lossless {
            0.0
        } else {
            state.entering - transmitted
        };
        let tolerance = 512.0 * f64::EPSILON;
        if !finite(transmission)
            || !transmitted.is_finite()
            || transmitted < 0.0
            || !absorbed.is_finite()
            || absorbed < -tolerance
            || (reflected + transmitted + absorbed - 1.0).abs() > tolerance
        {
            return Err(bad(
                "layer transmission failed passive normal-power balance",
            ));
        }
        Ok(LayeredTransmission {
            pressure_reflection: reflection,
            pressure_transmission: transmission,
            receiving_normal_wavenumber_rad_m: receiving_kn,
            reflected_fraction: reflected,
            transmitted_fraction: transmitted,
            absorbed_fraction: absorbed,
        })
    }
}

/// Propagate a homogeneous receiver boundary pair without dividing pressure by
/// velocity. Elastic mode conversion can have zero pressure or zero normal
/// velocity, so an admittance-only representation would introduce false poles.
pub(super) fn transmit_state(
    prepared: &RigidBackedLayers<'_>,
    cosine: f64,
    exit_pressure: C64,
    exit_velocity: C64,
) -> Result<TransmittedState, MaterialError> {
    let initial_scale = magnitude_scale(exit_pressure).max(magnitude_scale(exit_velocity));
    if !finite(exit_pressure) || !finite(exit_velocity) || initial_scale == 0.0 {
        return Err(bad("receiving pressure/velocity pair is unrepresentable"));
    }
    let mut pressure = divided(exit_pressure, initial_scale);
    let mut velocity = divided(exit_velocity, initial_scale);
    // [P,Q] = S*[physical pressure, Z_ref*physical velocity] for a unit
    // receiver modal amplitude. Retain the scalar canceled by reflection.
    let mut amplitude = LogAmplitude::new(-det::ln(initial_scale));
    for layer in prepared.layers.iter().rev() {
        let d = layer.source.thickness_m();
        if d == 0.0 {
            continue;
        }
        let kn;
        let normalization;
        (pressure, velocity, normalization, kn) =
            prepared.propagate_layer(layer, cosine, pressure, velocity)?;
        if kn.im == 0.0 {
            amplitude.add(-log_scaled(16.0, normalization));
        } else {
            amplitude.add(kn.im * d - log_scaled(8.0, normalization));
            amplitude.rotate(-kn.re * d);
        }
    }
    let (reflection, entering) = reflected_state(pressure, velocity)?;
    let denominator = pressure + velocity;
    let denominator_abs = denominator.abs();
    amplitude.add(core::f64::consts::LN_2 - det::ln(denominator_abs));
    amplitude.phase = amplitude.phase * divided(denominator.conj(), denominator_abs);
    Ok(TransmittedState {
        reflection,
        entering,
        amplitude,
    })
}

pub(super) struct TransmittedState {
    pub reflection: C64,
    pub entering: f64,
    amplitude: LogAmplitude,
}
impl TransmittedState {
    /// Receiver coefficient times its common transmitted amplitude. Multiply
    /// in log magnitude before rounding a deeply attenuated wave to zero.
    pub fn coefficient(&self, coefficient: C64) -> Result<C64, MaterialError> {
        self.scaled_coefficient(coefficient, 0.0)
    }

    pub fn scaled_coefficient(
        &self,
        coefficient: C64,
        log_scale: f64,
    ) -> Result<C64, MaterialError> {
        if coefficient == C64::ZERO {
            return Ok(C64::ZERO);
        }
        let scale = magnitude_scale(coefficient);
        if !finite(coefficient) || !scale.is_finite() {
            return Err(bad("layer transmission coefficient is unrepresentable"));
        }
        let direction = divided(coefficient, scale);
        let norm = direction.abs();
        let magnitude = det::exp(self.amplitude.log + log_scale + det::ln(scale) + det::ln(norm));
        if !magnitude.is_finite() {
            return Err(bad("layer transmission amplitude is unrepresentable"));
        }
        if magnitude == 0.0 {
            return Ok(C64::ZERO);
        }
        if !self.amplitude.phase_valid {
            return Err(bad("layer transmission phase is unrepresentable"));
        }
        Ok((self.amplitude.phase * divided(direction, norm)).scale(magnitude))
    }

    /// Normal mode power before squaring an underflowed field coefficient.
    pub fn power(&self, coefficient: C64, log_weight: f64) -> f64 {
        if coefficient == C64::ZERO || log_weight == f64::NEG_INFINITY {
            return 0.0;
        }
        let scale = magnitude_scale(coefficient);
        let log_coefficient = det::ln(scale) + det::ln(divided(coefficient, scale).abs());
        det::exp(log_weight + 2.0 * (self.amplitude.log + log_coefficient))
    }
}

fn positive_product_ratio(a: f64, b: f64, denominator: f64) -> Result<f64, MaterialError> {
    // Prefer removing the incidence cosine before multiplying a small density
    // contrast. A positive subnormal intermediate still loses relative accuracy.
    let ratio = b / denominator;
    let direct = a * ratio;
    let product = a * b;
    let alternate = product / denominator;
    let value = if ratio.is_normal() && direct.is_finite() && direct > 0.0 {
        direct
    } else if product.is_normal() && alternate.is_finite() && alternate > 0.0 {
        alternate
    } else {
        det::exp(det::ln(a) + det::ln(b) - det::ln(denominator))
    };
    if !value.is_finite() || value <= 0.0 {
        return Err(bad(
            "receiving normal admittance contrast is unrepresentable",
        ));
    }
    Ok(value)
}

fn log_scaled(factor: f64, value: f64) -> f64 {
    if value <= f64::MAX / factor {
        det::ln(factor * value)
    } else {
        det::ln(factor) + det::ln(value)
    }
}

// Compensated logarithmic magnitude prevents intermediate attenuation underflow;
// the separately normalized phase keeps the common complex factor explicit.
struct LogAmplitude {
    log: f64,
    correction: f64,
    phase: C64,
    phase_valid: bool,
}
impl LogAmplitude {
    fn new(log: f64) -> Self {
        Self {
            log,
            correction: 0.0,
            phase: C64::ONE,
            phase_valid: true,
        }
    }
    fn add(&mut self, value: f64) {
        if self.log == f64::NEG_INFINITY || value == f64::NEG_INFINITY {
            self.log = f64::NEG_INFINITY;
            self.correction = 0.0;
            return;
        }
        let corrected = value - self.correction;
        let sum = self.log + corrected;
        self.correction = if sum.is_finite() {
            (sum - self.log) - corrected
        } else {
            0.0
        };
        self.log = sum;
    }
    fn rotate(&mut self, radians: f64) {
        if self.log == f64::NEG_INFINITY {
            return;
        }
        if !radians.is_finite() {
            self.phase_valid = false;
            return;
        }
        let (sine, cosine) = det::sin_cos(radians);
        self.phase = self.phase * C64::new(cosine, sine);
        self.phase = divided(self.phase, self.phase.abs());
    }
}
