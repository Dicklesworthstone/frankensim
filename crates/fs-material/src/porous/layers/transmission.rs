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
    Layer, LayerTransfer, RigidBackedLayers, bad, divided, finite, magnitude_scale,
    normal_wavenumber_ratio, reflected_state, retained_fluid_snapshot,
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

    /// beta=Z0*v_n/p in exp(+i omega t), from the original scaled entrance
    /// pressure/velocity pair. Unlike reconstructing it from rounded R, this
    /// retains weak conductance when R rounds to unity. A pressure node has no
    /// finite Robin representation and refuses. No loss integration is done.
    pub fn normalized_surface_admittance(&self, cosine: f64) -> Result<C64, MaterialError> {
        let (state, _, _) =
            self.surface_state_with_layers(cosine, |_, _| Ok::<_, MaterialError>(()))?;
        state.normalized_surface_admittance(cosine)
    }

    /// Exact planar reflection/refraction/transmission within the admitted
    /// material model. Critical incidence (receiving k_n=0) remains finite.
    pub fn scattering(&self, cosine: f64) -> Result<LayeredTransmission, MaterialError> {
        self.scattering_with_layers(cosine, |_, _| Ok(()))
    }

    pub(super) fn scattering_with_layers<E: From<MaterialError>>(
        &self,
        cosine: f64,
        observe: impl FnMut(usize, LayerTransfer) -> Result<(), E>,
    ) -> Result<LayeredTransmission, E> {
        let (state, load, receiving_kn) = self.surface_state_with_layers(cosine, observe)?;
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
            return Err(bad("layer transmission failed passive normal-power balance").into());
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

    pub(super) fn surface_state_with_layers<E: From<MaterialError>>(
        &self,
        cosine: f64,
        observe: impl FnMut(usize, LayerTransfer) -> Result<(), E>,
    ) -> Result<(TransmittedState, C64, C64), E> {
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(bad("layer transmission needs incidence cosine in (0,1]").into());
        }
        let q = normal_wavenumber_ratio(
            self.speed_ratio,
            cosine,
            self.prepared.fluid().acoustic().sound_speed_m_s(),
            self.receiver.sound_speed_m_s(),
        );
        let receiving_kn = q.scale(self.prepared.k0);
        if !finite(receiving_kn) || (q != C64::ZERO && receiving_kn == C64::ZERO) {
            return Err(bad("receiving normal wavenumber is unrepresentable").into());
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
        let state = transmit_state_with_layers(&self.prepared, cosine, C64::ONE, load, observe)?;
        Ok((state, load, receiving_kn))
    }
}

/// Propagate a homogeneous receiver boundary pair without dividing pressure by
/// velocity. Elastic mode conversion can have zero pressure or zero normal
/// velocity, so an admittance-only representation would introduce false poles.
pub(super) fn transmit_state_with_layers<E: From<MaterialError>>(
    prepared: &RigidBackedLayers<'_>,
    cosine: f64,
    exit_pressure: C64,
    exit_velocity: C64,
    mut observe: impl FnMut(usize, LayerTransfer) -> Result<(), E>,
) -> Result<TransmittedState, E> {
    let initial_scale = magnitude_scale(exit_pressure).max(magnitude_scale(exit_velocity));
    if !finite(exit_pressure) || !finite(exit_velocity) || initial_scale == 0.0 {
        return Err(bad("receiving pressure/velocity pair is unrepresentable").into());
    }
    let mut pressure = divided(exit_pressure, initial_scale);
    let mut velocity = divided(exit_velocity, initial_scale);
    // [P,Q] = S*[physical pressure, Z_ref*physical velocity] for a unit
    // receiver modal amplitude. Retain the scalar canceled by reflection.
    let mut amplitude = LogAmplitude::new(-det::ln(initial_scale));
    for (index, layer) in prepared.layers.iter().enumerate().rev() {
        let d = layer.source.thickness_m();
        if d == 0.0 {
            observe(index, LayerTransfer::identity(pressure, velocity))?;
            continue;
        }
        let back_pressure = pressure;
        let back_velocity = velocity;
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
        observe(
            index,
            LayerTransfer {
                front_pressure: pressure,
                front_velocity: velocity,
                back_pressure,
                back_velocity,
                normalization,
                normal_wavenumber: kn,
            },
        )?;
    }
    TransmittedState::from_pair(pressure, velocity, amplitude).map_err(E::from)
}

pub(super) struct TransmittedState {
    pub reflection: C64,
    pub entering: f64,
    entrance_pressure: C64,
    entrance_velocity: C64,
    // Common scaled pressure/velocity factor before incident normalization.
    entrance_scale_log: f64,
    amplitude: LogAmplitude,
}
impl TransmittedState {
    /// Finish any bounded pressure/normal-velocity propagation in the same
    /// original incident normalization, retaining its accumulated scalar.
    pub(super) fn from_pair(
        pressure: C64,
        velocity: C64,
        mut amplitude: LogAmplitude,
    ) -> Result<Self, MaterialError> {
        let (reflection, entering) = reflected_state(pressure, velocity)?;
        let entrance_scale_log = amplitude.log;
        let denominator = pressure + velocity;
        let denominator_abs = denominator.abs();
        amplitude.add(core::f64::consts::LN_2 - det::ln(denominator_abs));
        amplitude.phase = amplitude.phase * divided(denominator.conj(), denominator_abs);
        Ok(Self {
            reflection,
            entering,
            entrance_pressure: pressure,
            entrance_velocity: velocity,
            entrance_scale_log,
            amplitude,
        })
    }

    pub(super) fn normalized_surface_admittance(&self, cosine: f64) -> Result<C64, MaterialError> {
        if self.entrance_pressure == C64::ZERO {
            return Err(bad("transmitting surface has a pressure node, not finite admittance"));
        }
        let beta = (self.entrance_velocity / self.entrance_pressure).scale(cosine);
        if !finite(beta) || beta.re < 0.0 {
            return Err(bad("transmitting surface admittance is nonfinite or active"));
        }
        Ok(beta)
    }

    /// Log amplitude for unit entrance pressure, before any squared field or
    /// incident-power fraction can underflow. The same normalization supplies
    /// independent local-loss integration and receiving-fluid power.
    pub(super) fn surface_pressure_log(&self) -> Result<f64, MaterialError> {
        if self.entrance_pressure == C64::ZERO {
            return Err(bad("transmitting surface power needs finite entrance pressure"));
        }
        Ok(-det::ln(self.entrance_pressure.abs()))
    }

    pub(super) fn transmitted_surface_conductance(
        &self,
        load: C64,
        cosine: f64,
    ) -> Result<f64, MaterialError> {
        if load.re == 0.0 {
            return Ok(0.0);
        }
        let logarithm = det::ln(cosine) + det::ln(load.re)
            + 2.0 * (self.entrance_scale_log + self.surface_pressure_log()?);
        let conductance = det::exp(logarithm);
        if !conductance.is_finite() || conductance < 0.0 {
            return Err(bad("transmitted surface conductance is unrepresentable"));
        }
        Ok(conductance)
    }

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
        self.amplitude.coefficient(coefficient, log_scale)
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

pub(super) fn positive_product_ratio(
    a: f64,
    b: f64,
    denominator: f64,
) -> Result<f64, MaterialError> {
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

pub(super) fn log_scaled(factor: f64, value: f64) -> f64 {
    if value <= f64::MAX / factor {
        det::ln(factor * value)
    } else {
        det::ln(factor) + det::ln(value)
    }
}

// Compensated logarithmic magnitude prevents intermediate attenuation underflow;
// the separately normalized phase keeps the common complex factor explicit.
#[derive(Clone)]
pub(super) struct LogAmplitude {
    pub(super) log: f64,
    correction: f64,
    phase: C64,
    phase_valid: bool,
}
impl LogAmplitude {
    pub(super) fn new(log: f64) -> Self {
        Self {
            log,
            correction: 0.0,
            phase: C64::ONE,
            phase_valid: true,
        }
    }
    pub(super) fn add(&mut self, value: f64) {
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
    pub(super) fn rotate(&mut self, radians: f64) {
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

    /// Multiply without rounding a representable downstream field to zero.
    pub(super) fn multiply(&mut self, coefficient: C64) -> Result<(), MaterialError> {
        if !finite(coefficient) {
            return Err(bad("layer field amplitude is unrepresentable"));
        }
        if coefficient == C64::ZERO {
            self.add(f64::NEG_INFINITY);
            return Ok(());
        }
        let scale = magnitude_scale(coefficient);
        let direction = divided(coefficient, scale);
        let norm = direction.abs();
        self.add(det::ln(scale) + det::ln(norm));
        self.phase = self.phase * divided(direction, norm);
        self.phase = divided(self.phase, self.phase.abs());
        Ok(())
    }

    /// Apply a retained complex scale only after including the field's units.
    pub(super) fn coefficient(
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
        let magnitude = det::exp(self.log + log_scale + det::ln(scale) + det::ln(norm));
        if !magnitude.is_finite() {
            return Err(bad("layer transmission amplitude is unrepresentable"));
        }
        if magnitude == 0.0 {
            return Ok(C64::ZERO);
        }
        if !self.phase_valid {
            return Err(bad("layer transmission phase is unrepresentable"));
        }
        Ok((self.phase * divided(direction, norm)).scale(magnitude))
    }
}

#[cfg(test)]
mod surface_tests {
    use super::*;
    use crate::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
    use crate::fluid_transport::FluidTransportState;
    use crate::porous::{PorousCard, PorousFluidLaw};
    use crate::porous::layers::absorption::{LayerAbsorptionBudget, LayerAbsorptionError};

    fn fluid() -> ThermoviscousFluidState {
        ThermoviscousFluidState::from_heat_capacity_ratio(
            FluidAcousticState::try_new(1.2, 340.0).unwrap(),
            FluidTransportState::try_new(1.2, 1005.0, 0.0257, 1.81e-5).unwrap(),
            1.4,
        ).unwrap()
    }

    #[test]
    fn surface_admittance_keeps_weak_receiving_conductance() {
        let gas = fluid();
        let receiver = FluidAcousticState::try_new(1.2e18, 340.0).unwrap();
        let stack = TransmittingLayers::new(&gas, &receiver, 2000.0, &[]).unwrap();
        let scattering = stack.scattering(1.0).unwrap();
        assert_eq!(scattering.pressure_reflection, C64::ONE);
        let beta = stack.normalized_surface_admittance(1.0).unwrap();
        assert!((beta.re / 1e-18 - 1.0).abs() < 1e-14);
        assert_eq!(beta.im, 0.0);
        let power = stack.normal_surface_power(
            LayerAbsorptionBudget { max_layers: 0 }, || false,
        ).unwrap();
        assert_eq!(power.absorbed_conductance, 0.0);
        assert!((power.transmitted_conductance / beta.re - 1.0).abs() < 2e-13);

        let node = TransmittedState::from_pair(
            C64::ZERO, C64::ONE, LogAmplitude::new(0.0),
        ).unwrap();
        assert!(node.normalized_surface_admittance(1.0).is_err());
    }

    #[test]
    fn transmitting_surface_separates_positive_heat_from_outgoing_power() {
        let gas = fluid();
        let receiver = gas.acoustic();
        let law = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
        let omega = core::f64::consts::TAU * 700.0;
        let lossless = TransmittingLayers::new(
            &gas, &receiver, omega, &[Layer::FluidGap { thickness_m: 0.03 }],
        ).unwrap();
        let budget = LayerAbsorptionBudget { max_layers: 1 };
        let power = lossless.normal_surface_power(budget, || false).unwrap();
        assert_eq!(power.absorbed_conductance, 0.0);
        assert!((power.transmitted_conductance - 1.0).abs() < 2e-13);
        for thickness_m in [0.03, 1e-20] {
            let stack = TransmittingLayers::new(
                &gas, &receiver, omega, &[Layer::Porous { law: &law, thickness_m }],
            ).unwrap();
            let power = stack.normal_surface_power(budget, || false).unwrap();
            let original = stack.layer_absorption(1.0, 1.0, budget, || false).unwrap();
            let scattering = stack.scattering(1.0).unwrap();
            let surface_pressure_squared = (C64::ONE + scattering.pressure_reflection).norm_sq();
            assert!(power.absorbed_conductance > 0.0, "even sub-roundoff heat must survive");
            assert!(power.transmitted_conductance > 0.0);
            assert!((power.absorbed_conductance * surface_pressure_squared
                / original.absorbed_fraction - 1.0).abs() < 2e-12);
            assert!((power.transmitted_conductance * surface_pressure_squared
                / original.transmitted_fraction - 1.0).abs() < 2e-12);
            assert!(power.balance_residual.abs() < 2e-12 * power.normalized_admittance.re);
        }
    }

    #[test]
    fn transmitting_surface_power_retains_native_work_and_cancellation_limits() {
        let gas = fluid();
        let receiver = gas.acoustic();
        let stack = TransmittingLayers::new(
            &gas, &receiver, 2000.0, &[Layer::FluidGap { thickness_m: 0.03 }],
        ).unwrap();
        assert!(matches!(
            stack.normal_surface_power(LayerAbsorptionBudget { max_layers: 0 }, || false),
            Err(LayerAbsorptionError::Invalid(_))
        ));
        let mut polls = 0;
        assert!(matches!(
            stack.normal_surface_power(LayerAbsorptionBudget { max_layers: 1 }, || {
                polls += 1;
                polls == 2
            }),
            Err(LayerAbsorptionError::Cancelled)
        ));
        assert!(stack.normalized_surface_admittance(0.0).is_err());
    }
}
