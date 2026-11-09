//! Original pHS forcing/material directions of a complete driven HB waveform.

use super::{
    HbOrbit, HbScales, PeriodicPortInput, PhsOrbitError, PhsStorageHessianAction, bad,
    forced_problem, tangent::StorageTangent,
};
use fs_orbit::sensitivity::{
    ForcedMeanObjective, ForcedMeanObjectiveControls, ForcedOrbitSensitivity,
    ForcedSensitivityControls, ForcedSensitivityDirection, ForcedSensitivityError,
    PreparedForcedSensitivity,
};
use fs_phs::{PhsError, PortHamiltonian};

mod adjoint;
pub use adjoint::PreparedPhsForcedAdjoint;
mod parameter;
pub(super) mod power;
pub use parameter::PhsForcedParameterDirection;
pub(super) use parameter::ode_clock_column;
pub use power::{PhsPowerObjective, PreparedPhsForcedPowerAdjoint};

/// Source admission and exact shared harmonic derivative refusals.
#[derive(Clone, Debug, PartialEq)]
pub enum PhsForcedSensitivityError {
    /// Original input/storage/source failed its authored domain.
    Source(PhsOrbitError),
    /// Original physical coefficient direction or contraction refused.
    Physical(PhsError),
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
impl From<PhsError> for PhsForcedSensitivityError {
    fn from(value: PhsError) -> Self {
        Self::Physical(value)
    }
}

/// One original pHS tangent/factor shared by arbitrary physical directions.
///
/// The source is borrowed unchanged. Every original nonlinear mode and explicit
/// memory coordinate participates. Storage-owned Hessians are discovered through
/// the same owner as forced HB; no linearized replacement Hamiltonian is built.
pub struct PreparedPhsForcedSensitivity<'a> {
    system: &'a PortHamiltonian,
    // Immutable original Fourier drive, bound to the prepared trajectory.
    input: PeriodicPortInput,
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
        mut controls: ForcedSensitivityControls,
        scales: &HbScales,
        hessian: Option<&PhsStorageHessianAction<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PhsForcedSensitivityError> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
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
        let coefficients = input
            .coefficients
            .iter()
            .map(|row| row.len() as u128)
            .sum::<u128>();
        let source_bytes = core::mem::size_of::<PeriodicPortInput>() as u128
            + (input.port_count() as u128) * core::mem::size_of::<Vec<fs_math::c64::C64>>() as u128
            + coefficients * core::mem::size_of::<fs_math::c64::C64>() as u128
            + 16 * (system.state_dim() as u128 + system.port_dim() as u128);
        if source_bytes > controls.maximum_payload_bytes as u128 {
            return Err(ForcedSensitivityError::Budget {
                quantity: "physical forced drive and direction payload",
                requested: source_bytes,
                cap: controls.maximum_payload_bytes as u128,
            }
            .into());
        }
        controls.maximum_payload_bytes -= source_bytes as usize;
        let mut problem = forced_problem(system, input)?;
        if let Some(action) = hessian {
            problem.tangent = StorageTangent::hessian(system.state_dim(), action);
        }
        let result =
            PreparedForcedSensitivity::new(&problem, orbit, controls, scales, &mut cancelled);
        if problem.source_refused() {
            return Err(PhsOrbitError::SourceRefused.into());
        }
        let prepared = result?;
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        Ok(Self {
            system,
            input: input.clone(),
            prepared,
        })
    }

    /// Shared full-state preparation diagnostics and complete residual-direction
    /// entry point, including source-dependent frequency ports when authored.
    pub fn prepared(&self) -> &PreparedForcedSensitivity {
        &self.prepared
    }

    /// Sample a differentiable physical observation on the original full-state
    /// cycle. The shared AFT owner checks both its mean and complete Fourier
    /// gradient on independent grids. The callback writes its value and EVERY
    /// state partial, including explicit material/controller memory coordinates.
    /// Direct parameter and clock dependence remain separate adjoint partials.
    pub fn mean_objective(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        observation: impl FnMut(f64, &[f64], &mut f64, &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedMeanObjective, PhsForcedSensitivityError> {
        self.prepared
            .mean_objective(controls, observation, cancelled)
            .map_err(Into::into)
    }

    /// Mean ORIGINAL Hamiltonian and its complete state covector, including
    /// nonlinear storage and explicit memory. A storage parameter's direct
    /// `partial_p H` must still be supplied to the objective adjoint separately.
    pub fn mean_stored_energy(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedMeanObjective, PhsForcedSensitivityError> {
        let system = self.system;
        self.mean_objective(
            controls,
            |_time, state, value, gradient| {
                *value = system.hamiltonian(state);
                system.storage().gradient(state, gradient);
                value.is_finite() && gradient.iter().all(|v| v.is_finite())
            },
            cancelled,
        )
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
        let action = port_direction(
            self.system,
            port_derivative,
            self.prepared.angular_frequency(),
            self.prepared.harmonics(),
        )?;
        self.prepared
            .solve_island_direction(direction, action, cancelled)
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
        effort_derivative: impl FnMut(&[f64], &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity, PhsForcedSensitivityError> {
        self.prepared
            .solve_island_direction(
                direction,
                storage_direction(self.system, effort_derivative),
                cancelled,
            )
            .map_err(Into::into)
    }
}

// The waveform solve and objective adjoint use the SAME original physical
// direction. No parallel assembly of G or (J-R), and no numerical source probes.
fn port_direction<'a>(
    system: &'a PortHamiltonian,
    input: &'a PeriodicPortInput,
    omega: f64,
    harmonics: usize,
) -> Result<impl FnMut(f64, &[f64], &mut [f64]) -> bool + 'a, PhsForcedSensitivityError> {
    if input.omega_rad_s() != omega
        || input.port_count() != system.port_dim()
        || input.highest_nonzero_harmonic() > harmonics
    {
        return Err(bad(
            "port derivative must preserve the clock, all ports and retained bandwidth",
        )
        .into());
    }
    let zero_effort = vec![0.0; system.state_dim()];
    let mut ports = vec![0.0; system.port_dim()];
    Ok(move |time: f64, _state: &[f64], out: &mut [f64]| {
        input.sample_into(time, &mut ports)
            && system
                .vector_field_from_effort_into(&zero_effort, &ports, out)
                .is_ok()
    })
}

fn storage_direction<'a>(
    system: &'a PortHamiltonian,
    mut effort_derivative: impl FnMut(&[f64], &mut [f64]) -> bool + 'a,
) -> impl FnMut(f64, &[f64], &mut [f64]) -> bool + 'a {
    let mut effort = vec![f64::NAN; system.state_dim()];
    let zero_ports = vec![0.0; system.port_dim()];
    move |_time: f64, state: &[f64], out: &mut [f64]| {
        effort.fill(f64::NAN);
        effort_derivative(state, &mut effort)
            && effort.iter().all(|value| value.is_finite())
            && system
                .vector_field_from_effort_into(&effort, &zero_ports, out)
                .is_ok()
    }
}

#[cfg(test)]
mod tests;
