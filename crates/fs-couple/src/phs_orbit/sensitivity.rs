//! Original pHS forcing/material directions of a complete driven HB waveform.

use super::{
    HbOrbit, HbScales, PeriodicPortInput, PhsOrbitError, PhsStorageHessianAction, bad,
    forced_problem, tangent::StorageTangent,
};
use fs_orbit::sensitivity::{
    ForcedOrbitSensitivity, ForcedSensitivityControls, ForcedSensitivityDirection,
    ForcedSensitivityError, PreparedForcedSensitivity,
};
use fs_phs::PortHamiltonian;

/// Source admission and exact shared harmonic derivative refusals.
#[derive(Clone, Debug, PartialEq)]
pub enum PhsForcedSensitivityError {
    /// Original input/storage/source failed its authored domain.
    Source(PhsOrbitError),
    /// Shared exact AFT derivative or numerical admission failed.
    Sensitivity(ForcedSensitivityError),
}
impl core::fmt::Display for PhsForcedSensitivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pHS forced harmonic sensitivity refused: {self:?}")
    }
}
impl std::error::Error for PhsForcedSensitivityError {}
impl From<PhsOrbitError> for PhsForcedSensitivityError {
    fn from(value: PhsOrbitError) -> Self {
        Self::Source(value)
    }
}
impl From<ForcedSensitivityError> for PhsForcedSensitivityError {
    fn from(value: ForcedSensitivityError) -> Self {
        Self::Sensitivity(value)
    }
}

/// One original pHS tangent/factor shared by arbitrary physical directions.
///
/// The source is borrowed unchanged. Every original nonlinear mode and explicit
/// memory coordinate participates. Storage-owned Hessians are discovered through
/// the same owner as forced HB; no linearized replacement Hamiltonian is built.
pub struct PreparedPhsForcedSensitivity<'a> {
    system: &'a PortHamiltonian,
    prepared: PreparedForcedSensitivity,
}

impl<'a> PreparedPhsForcedSensitivity<'a> {
    /// Prepare from the original source's native continuous storage tangent.
    /// The supplied orbit must use the same physical clock and cover the entire
    /// authored forcing bandwidth. All acceptance checks belong to fs-orbit.
    pub fn new(
        system: &'a PortHamiltonian,
        input: &PeriodicPortInput,
        orbit: &HbOrbit,
        controls: ForcedSensitivityControls,
        scales: &HbScales,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PhsForcedSensitivityError> {
        Self::prepare(system, input, orbit, controls, scales, None, cancelled)
    }

    /// Prepare with an explicitly authored COMPLETE original storage Hessian.
    /// Refusal or incomplete entries never fall back to numerical probes.
    pub fn new_analytic(
        system: &'a PortHamiltonian,
        input: &PeriodicPortInput,
        orbit: &HbOrbit,
        controls: ForcedSensitivityControls,
        scales: &HbScales,
        hessian: &PhsStorageHessianAction<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PhsForcedSensitivityError> {
        Self::prepare(
            system,
            input,
            orbit,
            controls,
            scales,
            Some(hessian),
            cancelled,
        )
    }

    fn prepare(
        system: &'a PortHamiltonian,
        input: &PeriodicPortInput,
        orbit: &HbOrbit,
        controls: ForcedSensitivityControls,
        scales: &HbScales,
        hessian: Option<&PhsStorageHessianAction<'_>>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PhsForcedSensitivityError> {
        let harmonics = orbit
            .coeffs
            .first()
            .map_or(0, |row| row.len().saturating_sub(1));
        if input.omega_rad_s() != orbit.omega || input.highest_nonzero_harmonic() > harmonics {
            return Err(bad(
                "sensitivity base must use the original clock and complete forcing bandwidth",
            )
            .into());
        }
        let mut problem = forced_problem(system, input)?;
        if let Some(action) = hessian {
            problem.tangent = StorageTangent::hessian(system.state_dim(), action);
        }
        let result = PreparedForcedSensitivity::new(&problem, orbit, controls, scales, cancelled);
        if problem.source_refused() {
            return Err(PhsOrbitError::SourceRefused.into());
        }
        Ok(Self {
            system,
            prepared: result?,
        })
    }

    /// Shared full-state preparation diagnostics and complete residual-direction
    /// entry point, including source-dependent frequency ports when authored.
    pub fn prepared(&self) -> &PreparedForcedSensitivity {
        &self.prepared
    }

    /// Differentiate all original physical port Fourier coefficients at once.
    ///
    /// `port_derivative` contains `du_n/dp`, with real DC and complex positive
    /// HALF-amplitudes. Its clock must equal the unchanged base clock; this entry
    /// requires zero `direction.angular_frequency_derivative`. Nonzero forcing
    /// above the retained order refuses, including arbitrarily small coefficients.
    /// The original pHS field owner forms `G du/dp`; the cached coupled tangent
    /// determines every state/harmonic response. No tangent/source re-fitting.
    pub fn solve_port_direction(
        &mut self,
        direction: ForcedSensitivityDirection,
        port_derivative: &PeriodicPortInput,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity, PhsForcedSensitivityError> {
        if port_derivative.omega_rad_s() != self.prepared.angular_frequency()
            || port_derivative.port_count() != self.system.port_dim()
            || port_derivative.highest_nonzero_harmonic() > self.prepared.harmonics()
        {
            return Err(bad(
                "port derivative must preserve the clock, all ports and retained bandwidth",
            )
            .into());
        }
        let zero_effort = vec![0.0; self.system.state_dim()];
        let mut ports = vec![0.0; self.system.port_dim()];
        let system = self.system;
        self.prepared
            .solve_island_direction(
                direction,
                |time, _state, out| {
                    port_derivative.sample_into(time, &mut ports)
                        && system
                            .vector_field_from_effort_into(&zero_effort, &ports, out)
                            .is_ok()
                },
                cancelled,
            )
            .map_err(Into::into)
    }

    /// Differentiate a physical storage/material parameter at fixed J, R, G,
    /// drive coefficients and clock. The callback writes the COMPLETE original
    /// `partial(gradient H)/partial p` at the unchanged physical state, including
    /// all nonlinear and explicit memory coordinates. It must not include the
    /// already-accounted-for Hessian times the unknown state derivative.
    ///
    /// The original field owner applies `(J-R)` to that effort direction. This
    /// admits source-authored stiffness, nonlinear and memory-material derivatives
    /// without duplicating their laws or reducing coupled modal coordinates.
    /// Changes of structure/port matrices use the complete residual entry instead.
    pub fn solve_storage_direction(
        &mut self,
        direction: ForcedSensitivityDirection,
        mut effort_derivative: impl FnMut(&[f64], &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity, PhsForcedSensitivityError> {
        let mut effort = vec![f64::NAN; self.system.state_dim()];
        let zero_ports = vec![0.0; self.system.port_dim()];
        let system = self.system;
        self.prepared
            .solve_island_direction(
                direction,
                |_time, state, out| {
                    effort.fill(f64::NAN);
                    effort_derivative(state, &mut effort)
                        && effort.iter().all(|value| value.is_finite())
                        && system
                            .vector_field_from_effort_into(&effort, &zero_ports, out)
                            .is_ok()
                },
                cancelled,
            )
            .map_err(Into::into)
    }
}

#[cfg(test)]
mod tests;
