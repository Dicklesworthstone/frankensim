//! The complete transient field of the original prescribed Rayleigh surface.
//!
//! Pressure and the three particle-velocity components observe the SAME
//! accepted source trace and physical time. In particular, velocity retains
//! both Q'/(2 pi c R) and Q/(2 pi R²), with each source's own propagation
//! direction and delay. Inferring velocity from p/(rho c) would lose the
//! reactive near field and would make a pressure null a false velocity null.
//!
//! This composes the existing differentiating and first-order receivers; it
//! introduces no propagation, source, quadrature or antialiasing law. A field
//! frame can be read through any ideal first-order microphone orientation
//! without resetting its history. The observation position remains fixed.
//! Signed instantaneous intensity is p*u after coherent field superposition,
//! not a sum of per-source intensities and not a time-averaged sound power.
//!
//! Scope: prescribed small-amplitude motion in a flat infinite rigid baffle,
//! inviscid homogeneous stationary fluid, and zero prehistory. No finite-body
//! scattering, material absorption, moving receiver or acoustic backreaction.
//! Source/surface refinement and the declared output bandwidth remain required.

use super::{BaffledPressure, CircularOutletReceiver, RayleighMedium, SurfaceSample};
use fs_bem::{helmholtz::Medium, near_field::FirstOrder};
use fs_math::c64::C64;

/// Co-located real fields at one accepted output sample, including the common
/// reported numerical latency. Particle velocity is fluid motion in the
/// source frame, not the source's normal velocity or a diaphragm velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldFrame {
    pub pressure_pa: f64,
    pub particle_velocity_m_s: [f64; 3],
    medium: RayleighMedium,
}

impl FieldFrame {
    /// The original common medium, retained for physical microphone scaling.
    pub const fn medium(&self) -> RayleighMedium {
        self.medium
    }

    /// Signed instantaneous Cartesian energy flux [W/m²]. Near-field storage
    /// and interference can reverse its direction during a cycle. There is no
    /// peak-phasor factor of one half for real time-domain samples.
    pub fn intensity_w_m2(&self) -> Result<[f64; 3], String> {
        let intensity = self.particle_velocity_m_s.map(|u| self.pressure_pa * u);
        if intensity.iter().any(|value| !value.is_finite()) {
            return Err("transient Rayleigh intensity is unrepresentable".into());
        }
        Ok(intensity)
    }

    /// Read this field through the existing ideal pressure/velocity pattern.
    /// The front axis points from receiver toward its front source. Pattern
    /// and orientation may change on every frame: this is a read-only local
    /// projection and never rebuilds, retimes or resets a source history.
    pub fn microphone_pressure_pa(&self, microphone: FirstOrder) -> Result<f64, String> {
        microphone
            .observe(
                C64::from_re(self.pressure_pa),
                self.particle_velocity_m_s.map(C64::from_re),
                Medium {
                    density: self.medium.density,
                    sound_speed: self.medium.sound_speed,
                },
            )
            .map(|pressure| pressure.re)
            .map_err(|error| error.to_string())
    }
}

/// Allocation-free warm observation of all four real acoustic field channels.
/// Cold preparation retains the pressure bank plus three native first-order
/// banks (each with its reactive companion); every bank keeps the existing
/// bounded history. No scalar pressure-to-velocity approximation is used.
#[derive(Clone)]
pub struct BaffledField {
    pressure: BaffledPressure,
    velocity: [BaffledPressure; 3],
    impedance: f64,
    last: FieldFrame,
}

impl BaffledField {
    /// Prepare the full field on one source/observation clock.
    pub fn new(
        surface: &[SurfaceSample],
        inputs: usize,
        rate: u32,
        position_m: [f64; 3],
        medium: RayleighMedium,
    ) -> Result<Self, String> {
        Self::new_multirate(surface, inputs, rate, rate, position_m, medium)
    }

    /// Observe every accepted mechanics substep before antialiasing to the
    /// output clock. The native integer clock ratio (1..=16), surface, medium,
    /// bandlimited derivative, fractional delay and history limits apply.
    /// Radiating and reactive terms retain their native common-time alignment,
    /// including fractional alignment for non-dyadic rate ratios.
    pub fn new_multirate(
        surface: &[SurfaceSample],
        inputs: usize,
        mechanics_rate: u32,
        output_rate: u32,
        position_m: [f64; 3],
        medium: RayleighMedium,
    ) -> Result<Self, String> {
        let impedance = medium.density * medium.sound_speed;
        if !impedance.is_finite() || impedance <= 0.0 {
            return Err(
                "transient Rayleigh field needs representable positive fluid impedance".into(),
            );
        }
        let pressure = BaffledPressure::new_multirate_differentiating(
            surface,
            inputs,
            mechanics_rate,
            output_rate,
            position_m,
            medium,
        )?;
        let component = |axis| {
            let mut front = [0.0; 3];
            front[axis] = -1.0;
            // alpha=0, front=-e_i gives +rho*c*u_i in the native owner.
            let pattern = FirstOrder::new(0.0, front).map_err(|error| error.to_string())?;
            BaffledPressure::new_multirate_first_order(
                surface,
                inputs,
                mechanics_rate,
                output_rate,
                position_m,
                medium,
                pattern,
            )
        };
        let velocity = [component(0)?, component(1)?, component(2)?];
        Ok(Self {
            pressure,
            velocity,
            impedance,
            last: FieldFrame {
                pressure_pa: 0.0,
                particle_velocity_m_s: [0.0; 3],
                medium,
            },
        })
    }

    /// Full field of a uniform circular outlet. Input is the actual total
    /// outward volume flow [m³/s]. The original positive-area disk quadrature
    /// and its spatial phase-resolution admission are shared with pressure.
    pub fn circular_outlet(
        radius_m: f64,
        mechanics_rate: u32,
        output_rate: u32,
        receiver: CircularOutletReceiver,
        medium: RayleighMedium,
    ) -> Result<Self, String> {
        if receiver.maximum_frequency_hz > 0.45 * f64::from(output_rate) {
            return Err("Rayleigh field outlet band exceeds its output reconstruction".into());
        }
        let surface = BaffledPressure::circular_surface(radius_m, receiver, medium)?;
        Self::new_multirate(
            &surface,
            1,
            mechanics_rate,
            output_rate,
            receiver.position_m,
            medium,
        )
    }

    /// Numerical filter latency in output samples, separate from physical
    /// travel times. All four channels represent this same delayed time.
    pub fn filter_delay_samples(&self) -> f64 {
        self.pressure.filter_delay_samples()
    }

    /// Last accepted complete field; a refused trace leaves it unchanged.
    pub const fn last_frame(&self) -> FieldFrame {
        self.last
    }

    /// One-rate convenience; multirate sources must use the complete trace.
    pub fn step(&mut self, velocity: &[f64]) -> Result<FieldFrame, String> {
        self.step_trace(velocity)
    }

    /// Interleaved END velocities of EVERY accepted mechanics substep.
    /// All component previews and SI conversions succeed before any accepted
    /// history advances. The warm path allocates no memory and performs no
    /// geometry, material, exponential or filter-design work.
    pub fn step_trace(&mut self, velocity: &[f64]) -> Result<FieldFrame, String> {
        let pressure_pa = self.pressure.preview_trace(velocity)?;
        let mut particle_velocity_m_s = [0.0; 3];
        for (value, receiver) in particle_velocity_m_s.iter_mut().zip(&mut self.velocity) {
            *value = receiver.preview_trace(velocity)? / self.impedance;
            if !value.is_finite() {
                return Err("transient Rayleigh particle velocity is unrepresentable".into());
            }
        }
        self.pressure.accept_trace(velocity);
        for receiver in &mut self.velocity {
            receiver.accept_trace(velocity);
        }
        self.last = FieldFrame {
            pressure_pa,
            particle_velocity_m_s,
            medium: self.last.medium,
        };
        Ok(self.last)
    }
}

#[cfg(test)]
#[path = "baffled_pressure_field_tests.rs"]
mod tests;
