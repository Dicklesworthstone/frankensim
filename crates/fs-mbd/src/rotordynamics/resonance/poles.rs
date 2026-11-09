//! Physical simple modes, force participation and shaft-speed sensitivities.
//!
//! The complete original pencil is T(s) = s^2 M + s(C + Omega G) + K.
//! The shared fs-modal owner supplies binary coordinate/rate equilibration,
//! physical unit-mass normalization, full left/right participation and original
//! source residual checks. No proportional-damping or whirl projection is used.
//!
//! Rotor harmonics use exp(+i omega t): s = i omega, omega_pole = -i s_pole,
//! and an angular-frequency Laurent residue is -i times its Laplace residue.
//! Stable poles therefore have NEGATIVE real s and POSITIVE imaginary omega.
//! Residues describe one simple pole's contribution, not a driven peak, a
//! complete modal expansion or a certificate of stability or observability.

use crate::rotordynamics::{MAX_ROTOR_DOFS, RotorOperatingPoint};
use fs_math::c64::C64;
use fs_modal::quadratic_modes::{
    PreparedQuadraticModes, PreparedQuadraticResponse, PreparedQuadraticSensitivities,
    QuadraticModeError, QuadraticParameterDirection, QuadraticResponseError,
};

pub use fs_modal::quadratic_modes::{
    PoleOptions, PoleResponseOptions, QuadraticForceResidue, QuadraticMode,
    QuadraticModeDerivative, QuadraticObservation, QuadraticResponseOptions,
};

/// One locally resolved simple mode at the declared signed shaft speed.
#[derive(Clone, Debug)]
pub struct RotorResonanceMode {
    /// Frozen shaft speed [rad/s] used in the original velocity operator.
    pub spin_rad_s: f64,
    /// Original displacement/velocity shapes and source-equation residual.
    /// The displacement satisfies phi^H M phi = 1 with a canonical phase.
    pub physical: QuadraticMode,
}

impl RotorResonanceMode {
    /// Complex angular frequency for the rotor's exp(+i omega t) convention.
    /// Both decay/growth and the signed oscillation frequency are retained.
    pub fn angular_frequency_rad_s(&self) -> C64 {
        let pole = self.physical.laplace_frequency_per_s;
        C64::new(pole.im, -pole.re)
    }
}

/// Refusal leaves the source and any previously prepared factor unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum RotorModeError {
    /// Mode metadata do not belong to the declared operating speed.
    Invalid(&'static str),
    /// Shared physical mode, source residual or resource admission refused.
    Mode(QuadraticModeError),
    /// Shared residue/sensitivity factor, batch or physical recovery refused.
    Response(QuadraticResponseError),
    /// Cancellation observed before any result is published.
    Cancelled,
}

impl core::fmt::Display for RotorModeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rotor resonance mode refused: {self:?}")
    }
}
impl std::error::Error for RotorModeError {}
impl From<QuadraticModeError> for RotorModeError {
    fn from(error: QuadraticModeError) -> Self {
        match error {
            QuadraticModeError::Cancelled => Self::Cancelled,
            error => Self::Mode(error),
        }
    }
}
impl From<QuadraticResponseError> for RotorModeError {
    fn from(error: QuadraticResponseError) -> Self {
        match error {
            QuadraticResponseError::Cancelled => Self::Cancelled,
            error => Self::Response(error),
        }
    }
}

/// Borrowed original pencil and reusable physical mass/coordinate preparation.
#[derive(Debug)]
pub struct PreparedRotorModes<'source> {
    spin_rad_s: f64,
    gyro: &'source [f64],
    numerical: PreparedQuadraticModes<'source>,
}

impl RotorOperatingPoint<'_> {
    /// Prepare local modes of the COMPLETE original M, C + Omega G and K.
    ///
    /// PoleOptions retains the shared owner's physical reference rate, full
    /// preparation/refinement payload and arithmetic admission. The original
    /// source arrays and process stack are excluded from its byte ceiling.
    /// Complex seeds are in the original displacement coordinates. Existing
    /// poles() supplies complete pole estimates; it does not supply shape seeds.
    pub fn prepare_resonance_modes(
        &self,
        options: PoleOptions,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRotorModes<'_>, RotorModeError> {
        let numerical = PreparedQuadraticModes::new(
            &self.model.mass,
            &self.velocity_operator,
            &self.model.stiffness,
            self.model.n,
            options,
            cancelled,
        )?;
        Ok(PreparedRotorModes {
            spin_rad_s: self.spin_rad_s,
            gyro: &self.model.gyro,
            numerical,
        })
    }
}

impl<'source> PreparedRotorModes<'source> {
    /// Refine one nearby physical Laplace pole and nonzero displacement seed.
    /// Repeated/defective roots may refuse; no mode identity across spin changes
    /// or complete spectrum is inferred. Active growth is never suppressed.
    pub fn refine(
        &self,
        laplace_seed_per_s: C64,
        displacement_seed: &[C64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<RotorResonanceMode, RotorModeError> {
        let physical = self.numerical.refine(
            laplace_seed_per_s,
            displacement_seed,
            cancelled,
        )?;
        Ok(RotorResonanceMode {
            spin_rad_s: self.spin_rad_s,
            physical,
        })
    }

    /// Recheck the mode against this original pencil and factor its full
    /// bordered linearization once for BOTH force residues and mode slopes.
    ///
    /// The shared response and sensitivity resource ceilings both apply;
    /// numerical.max_vectors bounds every subsequent force or direction batch.
    /// A mode from another speed refuses even if its eigenvalue coincides.
    pub fn prepare_response<'owner>(
        &'owner self,
        mode: &RotorResonanceMode,
        options: QuadraticResponseOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRotorPole<'owner, 'source>, RotorModeError> {
        if cancelled() {
            return Err(RotorModeError::Cancelled);
        }
        if mode.spin_rad_s != self.spin_rad_s {
            return Err(RotorModeError::Invalid(
                "mode and prepared source must have the same signed shaft speed",
            ));
        }
        let numerical = self.numerical.prepare_sensitivities(
            &mode.physical,
            options,
            &mut cancelled,
        )?;
        Ok(PreparedRotorPole {
            gyro: self.gyro,
            numerical,
        })
    }
}

/// One immutable shared bordered factor, borrowing its original source pencil.
/// Multiple force/sensor queries and shaft-speed slopes reuse that same factor.
#[derive(Debug)]
pub struct PreparedRotorPole<'owner, 'source> {
    gyro: &'source [f64],
    numerical: PreparedQuadraticSensitivities<'owner, 'source>,
}

impl PreparedRotorPole<'_, '_> {
    /// Canonical physical mode, source residual and underlying response factor.
    pub fn response(&self) -> &PreparedQuadraticResponse {
        self.numerical.response()
    }

    /// Residues for complete original-coordinate forces/torques and linear
    /// displacement/velocity observations, retaining true LEFT participation.
    ///
    /// Forces and sensor weights are their analytic values at the Laplace pole.
    /// Returned coefficients multiply 1/(s-s0), with velocity residue s0 times
    /// displacement residue. Under exp(+i omega t), multiply these coefficients
    /// by -i for residues in omega. These are neither amplitudes nor powers.
    /// Zero forces and dark receivers return zero without inventing coupling.
    pub fn forcing_residues(
        &self,
        forces: &[&[C64]],
        observations: &[QuadraticObservation<'_>],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<QuadraticForceResidue>, RotorModeError> {
        Ok(self.numerical.response().forcing_residues(
            forces,
            observations,
            cancelled,
        )?)
    }

    /// Analytic derivative with respect to signed shaft speed [rad/s], holding
    /// the supplied M, bearing C, K and per-speed G fixed. The full direction is
    /// dT/dOmega = s G; the pole, physical unit-mass shape and velocity all move.
    /// Speed-dependent bearings/preload need their own actual coefficient law.
    ///
    /// The shared factor and physical derivative audits are reused. The shape's
    /// baseline phase coordinate is held fixed; no crossing or moving-pivot
    /// derivative is claimed. Scratch is a bounded MAX_ROTOR_DOFS^2 real stack
    /// array, excluded from the shared heap-payload ceiling.
    pub fn spin_derivative(
        &self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<QuadraticModeDerivative, RotorModeError> {
        let zero = [0.0; MAX_ROTOR_DOFS * MAX_ROTOR_DOFS];
        let zero = &zero[..self.gyro.len()];
        let direction = QuadraticParameterDirection {
            mass: zero,
            damping: self.gyro,
            stiffness: zero,
        };
        let mut derivatives = self.numerical.parameter_derivatives(
            &[direction],
            cancelled,
        )?;
        derivatives.pop().ok_or(RotorModeError::Invalid(
            "shared owner returned no shaft-speed derivative",
        ))
    }
}
