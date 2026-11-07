//! Passive bound-charge polarization, with one shared exact oscillator owner.
//!
//! `P'' + gamma P' + omega_0^2 P = delta_epsilon omega_0^2 E` and
//! `U = (P^2 + (P'/omega_0)^2)/(2 delta_epsilon)` imply
//! `E P' = U' + gamma (P'/omega_0)^2/delta_epsilon`.
//! The coordinate `q=P/(omega_0 sqrt(delta_epsilon))` therefore maps exactly
//! to `fs-modal`'s unit-mass oscillator; no second integrator or linewidth
//! approximation is introduced. Here its work and energy are densities.
//!
//! Reference: COMSOL RF, "Drude-Lorentz Polarization",
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.rf/rf_ug_radio_frequency.07.068.html>.
//! This leaf admits explicit coefficients, not measured material authority.

use fs_math::{c64::C64, det};
use fs_modal::time::{LinearForceTransition, ModalAcousticState, ModalAcousticTimeError};

mod harmonic;
pub use harmonic::PreparedLorentzHarmonic;
mod transfer;
pub use transfer::LorentzHarmonicTransfer;
mod midpoint;
pub use midpoint::{LorentzMidpointPoleStep, PreparedLorentzMidpoint};
pub(super) use midpoint::{
    field_energy as midpoint_field_energy, field_increment as midpoint_field_increment,
};

/// Invalid physical data, an unresolved arithmetic range, or an actual pole.
#[derive(Clone, Debug, PartialEq)]
pub enum LorentzError {
    Invalid(&'static str),
    /// A lossless oscillator driven exactly at its natural frequency has no
    /// finite periodic response. Its transient still has a well-defined limit.
    SingularResonance,
    Evolution(ModalAcousticTimeError),
}
impl core::fmt::Display for LorentzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Lorentz polarization: {self:?}")
    }
}
impl std::error::Error for LorentzError {}
impl From<ModalAcousticTimeError> for LorentzError {
    fn from(value: ModalAcousticTimeError) -> Self {
        Self::Evolution(value)
    }
}
type Result<T> = core::result::Result<T, LorentzError>;

fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(LorentzError::Invalid(
            "polarization arithmetic is not representable",
        ))
    }
}
fn finite_complex(value: C64) -> Result<C64> {
    finite(value.re)?;
    finite(value.im)?;
    Ok(value)
}

/// One isotropic bound-charge process, before geometry reduction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LorentzPole {
    increment: f64,
    omega: f64,
    damping: f64,
    root_increment: f64,
    coupling: f64,
}

/// Original polarization and polarization current, in one to three dimensions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LorentzState<const N: usize> {
    /// Bound charge polarization [C/m2].
    pub polarization_c_m2: [f64; N],
    /// `dP/dt` [A/m2]; distinct from Ohmic current.
    pub current_a_m2: [f64; N],
}
impl<const N: usize> Default for LorentzState<N> {
    fn default() -> Self {
        Self {
            polarization_c_m2: [0.0; N],
            current_a_m2: [0.0; N],
        }
    }
}

/// Fixed-interval exact transition for a linearly interpolated electric field.
#[derive(Clone, Copy, Debug)]
pub struct PreparedLorentzStep {
    pole: LorentzPole,
    transition: LinearForceTransition,
}

/// Physical polarization energy and heat, without instantaneous field storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LorentzStep<const N: usize> {
    pub state: LorentzState<N>,
    pub stored_before_j_m3: f64,
    pub stored_after_j_m3: f64,
    /// Integral of `E(t) dot P'(t)` for the supplied electric-field history.
    pub supplied_work_j_m3: f64,
    /// Independent nonnegative integral of physical polarization friction.
    pub heat_j_m3: f64,
    pub heat_roundoff_tolerance_j_m3: f64,
    /// Stored after minus stored before plus heat minus supplied work.
    pub energy_residual_j_m3: f64,
}

impl LorentzPole {
    /// Positive absolute permittivity increment [F/m] and natural frequency
    /// [rad/s], with nonnegative physical damping rate [1/s]. Zero damping is
    /// an explicit lossless law. Free carriers require a separate Drude owner.
    pub fn try_new(
        permittivity_increment_f_m: f64,
        angular_frequency_rad_s: f64,
        damping_rate_per_s: f64,
    ) -> Result<Self> {
        if !permittivity_increment_f_m.is_finite()
            || permittivity_increment_f_m <= 0.0
            || !angular_frequency_rad_s.is_finite()
            || angular_frequency_rad_s <= 0.0
            || !damping_rate_per_s.is_finite()
            || damping_rate_per_s < 0.0
        {
            return Err(LorentzError::Invalid(
                "positive finite increment/frequency and nonnegative damping required",
            ));
        }
        let root_increment = det::sqrt(permittivity_increment_f_m);
        let coupling = finite(angular_frequency_rad_s * root_increment)?;
        if coupling == 0.0 {
            return Err(LorentzError::Invalid(
                "polarization normalization underflows",
            ));
        }
        Ok(Self {
            increment: permittivity_increment_f_m,
            omega: angular_frequency_rad_s,
            damping: damping_rate_per_s,
            root_increment,
            coupling,
        })
    }

    pub const fn permittivity_increment_f_m(self) -> f64 {
        self.increment
    }
    pub const fn angular_frequency_rad_s(self) -> f64 {
        self.omega
    }
    pub const fn damping_rate_per_s(self) -> f64 {
        self.damping
    }

    /// Absolute susceptibility [F/m] under `exp(+i omega t)`. Factored detuning
    /// retains high-Q near-resonance differences. Scaling avoids squaring a
    /// dimensional optical frequency or silently overflowing its denominator.
    pub fn susceptibility_f_m(self, angular_frequency_rad_s: f64) -> Result<C64> {
        let omega = angular_frequency_rad_s;
        if !omega.is_finite() || omega < 0.0 {
            return Err(LorentzError::Invalid(
                "finite nonnegative driving angular frequency required",
            ));
        }
        if omega == 0.0 {
            return Ok(C64::from_re(self.increment));
        }
        if omega == self.omega && self.damping == 0.0 {
            return Err(LorentzError::SingularResonance);
        }
        let scale = self.omega.max(omega).max(self.damping);
        let natural = self.omega / scale;
        let driven = omega / scale;
        let real = ((self.omega - omega) / scale) * (natural + driven);
        let imaginary = (self.damping / scale) * driven;
        let strength = (self.increment * natural) * natural;
        if strength == 0.0
            || (self.damping > 0.0 && imaginary == 0.0)
            || (real == 0.0 && omega != self.omega)
        {
            return Err(LorentzError::Invalid(
                "scaled polarization response loses a nonzero coefficient",
            ));
        }
        finite_complex(C64::from_re(strength) / C64::new(real, imaginary))
    }

    /// Full bound-charge kinetic plus restoring energy [J/m3], not epsilon
    /// times field squared. Inertial storage remains positive above resonance
    /// even when the real part of effective permittivity is negative.
    pub fn stored_energy_j_m3<const N: usize>(self, state: LorentzState<N>) -> Result<f64> {
        if !(1..=3).contains(&N) {
            return Err(LorentzError::Invalid(
                "one to three polarization components required",
            ));
        }
        let mut energy = 0.0;
        for i in 0..N {
            let position = finite(state.polarization_c_m2[i] / self.root_increment)?;
            let velocity = finite(state.current_a_m2[i] / self.coupling)?;
            energy = finite(energy + (0.5 * position) * position + (0.5 * velocity) * velocity)?;
        }
        Ok(energy)
    }

    /// Prepare once for repeated field samples. No Nyquist restriction is
    /// invented; the caller must resolve its actual field and observations.
    pub fn prepare_step(self, duration_s: f64) -> Result<PreparedLorentzStep> {
        Ok(PreparedLorentzStep {
            pole: self,
            transition: LinearForceTransition::try_new_with_damping_rate(
                self.omega,
                self.damping,
                duration_s,
            )?,
        })
    }
}

impl PreparedLorentzStep {
    /// Advance without allocation or transcendental calls. The input state is
    /// immutable and no partial vector is returned on a failed component.
    pub fn advance<const N: usize>(
        &self,
        before: LorentzState<N>,
        field_start_v_m: [f64; N],
        field_end_v_m: [f64; N],
    ) -> Result<LorentzStep<N>> {
        let stored_before = self.pole.stored_energy_j_m3(before)?;
        let mut state = LorentzState::default();
        let mut work = 0.0;
        let mut heat = 0.0;
        let mut tolerance = 0.0;
        for i in 0..N {
            let motion = self.transition.advance_with_energy(
                ModalAcousticState {
                    displacement_m_sqrt_kg: finite(
                        before.polarization_c_m2[i] / self.pole.coupling,
                    )?,
                    velocity_m_sqrt_kg_per_s: finite(before.current_a_m2[i] / self.pole.coupling)?,
                },
                finite(self.pole.coupling * field_start_v_m[i])?,
                finite(self.pole.coupling * field_end_v_m[i])?,
            )?;
            state.polarization_c_m2[i] =
                finite(self.pole.coupling * motion.state.displacement_m_sqrt_kg)?;
            state.current_a_m2[i] =
                finite(self.pole.coupling * motion.state.velocity_m_sqrt_kg_per_s)?;
            work = finite(work + motion.input_work_j)?;
            heat = finite(heat + motion.viscous_dissipation_j)?;
            tolerance = finite(tolerance + motion.dissipation_roundoff_tolerance_j)?;
        }
        let stored_after = self.pole.stored_energy_j_m3(state)?;
        let residual = finite((stored_after - stored_before) + (heat - work))?;
        Ok(LorentzStep {
            state,
            stored_before_j_m3: stored_before,
            stored_after_j_m3: stored_after,
            supplied_work_j_m3: work,
            heat_j_m3: heat,
            heat_roundoff_tolerance_j_m3: tolerance,
            energy_residual_j_m3: residual,
        })
    }
}

#[cfg(test)]
#[path = "oscillator_tests.rs"]
mod tests;
