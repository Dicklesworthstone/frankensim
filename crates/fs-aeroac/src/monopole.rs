//! Causal, three-dimensional compact volume sources in a homogeneous fluid.
//!
//! The outgoing linear monopole is `p'(r,t+r/c) = rho V''(t)/(4 pi r)`.
//! The volume acceleration is a physical source input, not an incompressible
//! solver's pressure or an authored sound effect. In particular a spherical
//! bubble supplies `V'' = 4 pi R (R R'' + 2 R'^2)` from its dynamics.
//! See Brennen, *Cavitation and Bubble Dynamics*, chapters 2 and 3:
//! <https://media.library.caltech.edu/CaltechBOOK:1995.001/>.
//!
//! This is a one-way, free-field observation. It does not put acoustic radiation
//! damping back into the source, model shocks or boundaries, or certify measured
//! sound levels. The caller declares the source bandwidth; it must include
//! nonlinear harmonics, not just the excitation frequency. Each sample retains
//! its arrival time, so propagation never becomes an instantaneous coupling.

use crate::AeroacError;
use fs_material::fluid_acoustics::FluidAcousticState;

/// Explicit approximation limits for a compact, distant source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MonopoleLimits {
    /// Largest admitted `2 pi f_max R / c`; must be in `(0, 0.3]`.
    pub maximum_ka: f64,
    /// Smallest admitted observer distance/source radius; at least ten.
    pub minimum_distance_radii: f64,
    /// Largest admitted `|rho V''/(4 pi R)|/K_s` at the source surface,
    /// in `(0, 0.1]`. Moving the observer farther away cannot hide nonlinear
    /// source strength from this conservative linear-propagation check.
    pub maximum_pressure_ratio: f64,
}

/// One physical pressure sample and its causal observation time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MonopoleSample {
    /// Accepted source time [s].
    pub source_time_s: f64,
    /// Observer time `source_time + distance/c` [s].
    pub arrival_time_s: f64,
    /// Acoustic pressure perturbation [Pa], relative to the ambient state.
    pub pressure_pa: f64,
    /// The physical source volume acceleration [m³/s²].
    pub volume_acceleration_m3_s2: f64,
    /// Actual source compactness at the declared maximum frequency.
    pub ka: f64,
}

/// Prepared propagation from a fixed source center to a fixed observer.
///
/// Density and sound speed are the same checked material snapshot. There is no
/// atmospheric default and no gas-only equation of state. Preparation and
/// sampling do not allocate. Sources and any original material receipts remain
/// owned by the caller, as in the ordinary fluid-acoustic propagation APIs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactMonopoleObserver {
    medium: FluidAcousticState,
    distance_m: f64,
    maximum_frequency_hz: f64,
    delay_s: f64,
    pressure_per_volume_acceleration: f64,
    limits: MonopoleLimits,
}

fn invalid(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}

impl CompactMonopoleObserver {
    /// Prepare the physical `1/r` transfer under a caller-declared bandwidth.
    pub fn new(
        medium: FluidAcousticState,
        distance_m: f64,
        maximum_frequency_hz: f64,
        limits: MonopoleLimits,
    ) -> Result<Self, AeroacError> {
        if !distance_m.is_finite()
            || distance_m <= 0.0
            || !maximum_frequency_hz.is_finite()
            || maximum_frequency_hz <= 0.0
            || !limits.maximum_ka.is_finite()
            || limits.maximum_ka <= 0.0
            || limits.maximum_ka > 0.3
            || !limits.minimum_distance_radii.is_finite()
            || limits.minimum_distance_radii < 10.0
            || !limits.maximum_pressure_ratio.is_finite()
            || limits.maximum_pressure_ratio <= 0.0
            || limits.maximum_pressure_ratio > 0.1
        {
            return Err(invalid(
                "invalid compact-monopole distance, bandwidth or limits",
            ));
        }
        let delay_s = distance_m / medium.sound_speed_m_s();
        let pressure_per_volume_acceleration =
            (medium.density_kg_m3() / distance_m) / (4.0 * core::f64::consts::PI);
        if !delay_s.is_finite()
            || delay_s <= 0.0
            || !pressure_per_volume_acceleration.is_finite()
            || pressure_per_volume_acceleration <= 0.0
        {
            return Err(invalid("compact-monopole transfer is not representable"));
        }
        Ok(Self {
            medium,
            distance_m,
            maximum_frequency_hz,
            delay_s,
            pressure_per_volume_acceleration,
            limits,
        })
    }

    /// Original coherent local propagation coefficients.
    pub const fn medium(&self) -> FluidAcousticState {
        self.medium
    }

    /// Travel time from source center to observer [s].
    pub const fn propagation_delay_s(&self) -> f64 {
        self.delay_s
    }

    /// Observe a physical volume acceleration at its retarded arrival time.
    ///
    /// The radius bounds the entire compact source. An expanding bubble must
    /// therefore supply its current radius, not its initial equilibrium size.
    /// This routine does not interpolate, resample, or modify the source state.
    pub fn sample(
        &self,
        source_time_s: f64,
        source_radius_m: f64,
        volume_acceleration_m3_s2: f64,
    ) -> Result<MonopoleSample, AeroacError> {
        if !source_time_s.is_finite()
            || !source_radius_m.is_finite()
            || source_radius_m <= 0.0
            || !volume_acceleration_m3_s2.is_finite()
        {
            return Err(invalid("invalid compact-monopole source sample"));
        }
        let ka = (source_radius_m / self.medium.sound_speed_m_s())
            * self.maximum_frequency_hz
            * (2.0 * core::f64::consts::PI);
        if !ka.is_finite()
            || ka <= 0.0
            || ka > self.limits.maximum_ka
            || source_radius_m > self.distance_m / self.limits.minimum_distance_radii
        {
            return Err(invalid("source exceeds compact or distant-observer domain"));
        }
        let arrival_time_s = source_time_s + self.delay_s;
        let pressure_pa = self.pressure_per_volume_acceleration * volume_acceleration_m3_s2;
        let source_pressure_ratio = ((volume_acceleration_m3_s2 / source_radius_m).abs()
            / self.medium.sound_speed_m_s())
            / self.medium.sound_speed_m_s()
            / (4.0 * core::f64::consts::PI);
        if !arrival_time_s.is_finite()
            || arrival_time_s <= source_time_s
            || !pressure_pa.is_finite()
            || (volume_acceleration_m3_s2 != 0.0 && pressure_pa == 0.0)
            || !source_pressure_ratio.is_finite()
            || source_pressure_ratio > self.limits.maximum_pressure_ratio
        {
            return Err(invalid(
                "compact-monopole pressure or arrival time is not representable",
            ));
        }
        Ok(MonopoleSample {
            source_time_s,
            arrival_time_s,
            pressure_pa,
            volume_acceleration_m3_s2,
            ka,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn observer(distance: f64) -> CompactMonopoleObserver {
        CompactMonopoleObserver::new(
            FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
            distance,
            1000.0,
            MonopoleLimits {
                maximum_ka: 0.1,
                minimum_distance_radii: 10.0,
                maximum_pressure_ratio: 0.01,
            },
        )
        .unwrap()
    }

    #[test]
    fn harmonic_volume_has_analytic_pressure_and_retarded_time() {
        let radius: f64 = 0.001;
        let amplitude: f64 = 1e-8;
        let omega = 2.0 * core::f64::consts::PI * 1000.0;
        // At the maximum of V = V0 + A cos(omega t), V'' = -omega² A.
        let sample = observer(2.0)
            .sample(0.0, radius, -omega * omega * amplitude)
            .unwrap();
        let expected = -1000.0 * core::f64::consts::PI * 1e6 * amplitude / 2.0;
        assert!((sample.pressure_pa - expected).abs() < 2e-14);
        assert_eq!(sample.arrival_time_s, 2.0 / 1500.0);
        assert!(sample.pressure_pa < 0.0);
    }

    #[test]
    fn spherical_spreading_preserves_radiative_intensity_scaling() {
        let near = observer(1.0).sample(0.0, 1e-3, 0.03).unwrap();
        let far = observer(4.0).sample(0.0, 1e-3, 0.03).unwrap();
        assert_eq!(near.pressure_pa, 4.0 * far.pressure_pa);
        let power_near = 4.0 * core::f64::consts::PI * near.pressure_pa.powi(2) / 1.5e6;
        let power_far = 4.0 * core::f64::consts::PI * 16.0 * far.pressure_pa.powi(2) / 1.5e6;
        assert_eq!(power_near, power_far);
        assert_eq!(far.arrival_time_s, 4.0 * near.arrival_time_s);
    }

    #[test]
    fn physical_medium_controls_both_strength_and_delay() {
        let other = CompactMonopoleObserver::new(
            FluidAcousticState::try_new(500.0, 750.0).unwrap(),
            1.0,
            1000.0,
            MonopoleLimits {
                maximum_ka: 0.1,
                minimum_distance_radii: 10.0,
                maximum_pressure_ratio: 0.01,
            },
        )
        .unwrap();
        let a = observer(1.0).sample(0.0, 1e-3, 0.03).unwrap();
        let b = other.sample(0.0, 1e-3, 0.03).unwrap();
        assert_eq!(a.pressure_pa, 2.0 * b.pressure_pa);
        assert_eq!(b.arrival_time_s, 2.0 * a.arrival_time_s);
    }

    #[test]
    fn refuses_near_field_unresolved_source_and_unrepresentable_delay() {
        assert!(observer(0.005).sample(0.0, 0.001, 1.0).is_err());
        assert!(observer(2.0).sample(0.0, 0.1, 1.0).is_err());
        assert!(observer(1.0).sample(1e30, 1e-3, 1.0).is_err());
        assert!(observer(1.0).sample(0.0, 1e-3, f64::NAN).is_err());
        assert!(observer(1.0).sample(0.0, 1e-3, 1e12).is_err());
        assert!(observer(1e6).sample(0.0, 1e-3, 1e12).is_err());
        assert_eq!(
            observer(1.0).sample(0.0, 1e-3, 0.0).unwrap().pressure_pa,
            0.0
        );
    }
}
