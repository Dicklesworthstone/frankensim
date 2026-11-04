//! Scalar harmonic motion with independently integrated work and viscous heat.

use super::{HarmonicForceTransition, integrals::HarmonicIntegrals};
use crate::time::{ModalAcousticState, ModalAcousticTimeError};

/// Endpoint and physical power integrals under one continuous force tone.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicForceFrame {
    /// Motion at the interval endpoint in mass-normalized coordinates.
    pub state: ModalAcousticState,
    /// Signed integral of `F(t)*v(t)` [J], including quadrature interference.
    pub input_work_j: f64,
    /// Scale-aware floating-point allowance for force work [J].
    pub input_work_roundoff_tolerance_j: f64,
    /// Nonnegative integral of the physical viscous power `c*v(t)^2` [J].
    pub viscous_dissipation_j: f64,
    /// Scale-aware floating-point allowance for viscous heat [J].
    pub dissipation_roundoff_tolerance_j: f64,
}

/// Prepared exact scalar harmonic motion, signed force work and viscous heat.
///
/// The equation is `q'' + c*q' + omega^2*q = F(t)`, with
/// `F(t)=Fc*cos(Omega*t)+Fs*sin(Omega*t)` measured from the supplied initial
/// state. Forces and states use the mass-normalized units of
/// [`HarmonicForceTransition`]. A damping ratio supplies `c=2*zeta*omega`;
/// [`Self::try_new_with_damping_rate`] admits a physical `c` directly,
/// including viscous drag on an exactly zero-stiffness coordinate.
///
/// The motion and integral owners are shared with [`super::PreparedModalHarmonic`].
/// Preparation resolves a short interval before composing it to the requested
/// duration. Heat is a positive square-root Gramian; work is a separate signed
/// quadratic form. Neither is obtained by subtracting endpoint energies.
/// Warm advancement allocates nothing and evaluates no transcendental functions.
/// The smaller motion-only transition does not carry these power forms.
///
/// Exact refers to the continuous linear forcing law. Returned roundoff
/// allowances are numerical estimates, not rigorous enclosures. No acoustic
/// observation or state/energy budget is inferred. This is one frequency:
/// heats from independent tones cannot be summed to model simultaneous tones,
/// whose velocities have additional cross terms. Use [`super::PreparedModalTones`]
/// for simultaneous-tone work and heat in a ratio-parameterized modal bank.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicForce {
    motion: HarmonicForceTransition,
    integrals: HarmonicIntegrals,
}

impl PreparedHarmonicForce {
    /// Prepare one interval using a nonnegative damping ratio.
    /// A zero natural frequency requires a zero ratio, just as in
    /// [`HarmonicForceTransition::try_new`]. Call outside a real-time callback.
    ///
    /// # Errors
    /// Invalid inputs or unrepresentable motion/work/heat coefficients.
    pub fn try_new(
        angular_frequency_rad_s: f64,
        damping_ratio: f64,
        forcing_angular_frequency_rad_s: f64,
        duration_s: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        let motion = HarmonicForceTransition::try_new(
            angular_frequency_rad_s,
            damping_ratio,
            forcing_angular_frequency_rad_s,
            duration_s,
        )?;
        Self::prepare(motion, forcing_angular_frequency_rad_s)
    }

    /// Prepare with physical viscous rate `c` [1/s], including zero stiffness.
    /// For mass `m` [kg] and drag coefficient `b` [kg/s], supply `c=b/m` and
    /// mass-normalize state and forces. No damping ratio or stiffness floor is
    /// needed. Call outside a real-time callback.
    ///
    /// # Errors
    /// Invalid inputs or unrepresentable motion/work/heat coefficients.
    pub fn try_new_with_damping_rate(
        angular_frequency_rad_s: f64,
        damping_rate_per_s: f64,
        forcing_angular_frequency_rad_s: f64,
        duration_s: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        let motion = HarmonicForceTransition::try_new_with_damping_rate(
            angular_frequency_rad_s,
            damping_rate_per_s,
            forcing_angular_frequency_rad_s,
            duration_s,
        )?;
        Self::prepare(motion, forcing_angular_frequency_rad_s)
    }

    fn prepare(
        motion: HarmonicForceTransition,
        forcing_angular_frequency_rad_s: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        let integrals = HarmonicIntegrals::try_new(&motion, forcing_angular_frequency_rad_s)?;
        Ok(Self { motion, integrals })
    }

    /// Borrow the exact motion and phase transition used by the energy step.
    /// Its pure preview and phase rotation leave this preparation unchanged.
    #[must_use]
    pub fn motion(&self) -> &HarmonicForceTransition {
        &self.motion
    }

    /// Advance a supplied state and integrate signed work and positive heat.
    /// Nothing is retained or mutated, so rejected trials and checkpoint replay
    /// cannot advance accepted state or forcing phase. After accepting a step,
    /// rotate its forces with [`HarmonicForceTransition::advance_force_phase`]
    /// through [`Self::motion`].
    ///
    /// # Errors
    /// Nonfinite state/forces or unrepresentable motion or power arithmetic.
    pub fn advance_with_energy(
        &self,
        state: ModalAcousticState,
        cosine_force: f64,
        sine_force: f64,
    ) -> Result<HarmonicForceFrame, ModalAcousticTimeError> {
        let next = self.motion.advance(state, cosine_force, sine_force)?;
        let (work, work_allowance, heat, heat_allowance) =
            self.integrals.evaluate(state, cosine_force, sine_force)?;
        Ok(HarmonicForceFrame {
            state: next,
            input_work_j: work,
            input_work_roundoff_tolerance_j: work_allowance,
            viscous_dissipation_j: heat,
            dissipation_roundoff_tolerance_j: heat_allowance,
        })
    }
}
