//! Complete physical directions through the original self-sustained loss law.

use super::{Adapter, DissipativeAutonomousHamiltonian, DissipativeError};
use crate::phs_orbit::PhsAutonomousParameterDirection;
use fs_orbit::sensitivity::{
    AutonomousMeanObjective, AutonomousObjectiveGradient, AutonomousOrbitSensitivity,
    ForcedMeanObjectiveControls, ForcedSensitivityControls, ForcedSensitivityError,
    PreparedAutonomousAdjoint, PreparedAutonomousSensitivity,
};
use fs_orbit::{HbOrbit, HbScales};
use fs_phs::harmonic::equilibrium::{DissipationParameterAction, LinearizationError};

mod parameter;
mod power;
pub use power::{DissipativeAutonomousPowerObjective, PreparedDissipativeAutonomousPowerAdjoint};

/// One original physical parameter direction, holding the state and frequency
/// fixed. The existing feedback derivative is TOTAL, including the physical
/// output's G/H chain. The separate loss callback is D_p at fixed state AND
/// effort; D_e * partial_p(grad H) is supplied automatically exactly once.
pub struct DissipativeAutonomousParameterDirection<'a> {
    pub physical: PhsAutonomousParameterDirection<'a>,
    /// None explicitly holds the original D law fixed. Its effort chain still
    /// participates when storage changes. A present callback fills every state.
    pub dissipation: Option<&'a DissipationParameterAction<'a>>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DissipativeAutonomousSensitivityError {
    Source(DissipativeError),
    Physical(fs_phs::PhsError),
    Sensitivity(ForcedSensitivityError),
}
impl core::fmt::Display for DissipativeAutonomousSensitivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "self-sustained nonlinear-loss sensitivity refused: {self:?}"
        )
    }
}
impl std::error::Error for DissipativeAutonomousSensitivityError {}
impl From<DissipativeError> for DissipativeAutonomousSensitivityError {
    fn from(value: DissipativeError) -> Self {
        Self::Source(value)
    }
}
impl From<fs_phs::PhsError> for DissipativeAutonomousSensitivityError {
    fn from(value: fs_phs::PhsError) -> Self {
        Self::Physical(value)
    }
}
impl From<ForcedSensitivityError> for DissipativeAutonomousSensitivityError {
    fn from(value: ForcedSensitivityError) -> Self {
        Self::Sensitivity(value)
    }
}
impl From<LinearizationError> for DissipativeAutonomousSensitivityError {
    fn from(value: LinearizationError) -> Self {
        Self::Source(DissipativeError::Dissipation(value))
    }
}
type Result<T> = core::result::Result<T, DissipativeAutonomousSensitivityError>;

/// One original bordered tangent/factor for complete waveform, frequency and
/// objective derivatives. Every source and authored tangent is bound unchanged
/// to the original prepared Hamiltonian. No neighboring model or refit is used.
pub struct PreparedDissipativeAutonomousSensitivity<'s, 'a> {
    source: &'s DissipativeAutonomousHamiltonian<'a>,
    prepared: PreparedAutonomousSensitivity,
}

impl<'a> DissipativeAutonomousHamiltonian<'a> {
    /// Prepare the unchanged full-state balance and full-waveform phase border.
    /// The original nonlinear-loss tangent reaches both independent AFT grids.
    /// Natural frequency remains an unknown; no neutral phase equation is
    /// discarded. Phase reference is a complete autonomous packed seed.
    ///
    /// Reserve 6*n+4*m scalar query scratch in addition to the generic owner's
    /// AFT/LU payload. Already-prepared source/tangent storage and callback work
    /// retain their source-defined limits. Singular/nonisolated cycles, bad
    /// source derivatives and unresolved finite-grid evidence refuse.
    pub fn prepare_sensitivity(
        &self,
        orbit: &HbOrbit,
        phase_reference: &[f64],
        mut controls: ForcedSensitivityControls,
        scales: &HbScales,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDissipativeAutonomousSensitivity<'_, 'a>> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let (n, m) = (self.system.state_dim(), self.system.port_dim());
        let bytes = 8_u128 * (6 * n as u128 + 4 * m as u128);
        if bytes > controls.maximum_payload_bytes as u128 {
            return Err(ForcedSensitivityError::Budget {
                quantity: "self-sustained nonlinear-loss parameter scratch",
                requested: bytes,
                cap: controls.maximum_payload_bytes as u128,
            }
            .into());
        }
        controls.maximum_payload_bytes -= bytes as usize;
        *self.failure.borrow_mut() = None;
        let prepared = PreparedAutonomousSensitivity::new(
            &Adapter(self),
            orbit,
            phase_reference,
            controls,
            scales,
            |point, out| {
                // This original full-state source has the ODE frequency port sI.
                // Its frequency column is i*k*X_k, with all memories retained.
                out.fill(0.0);
                for harmonic in 1..=point.harmonics {
                    for state in 0..n {
                        let re = n + (harmonic - 1) * 2 * n + state;
                        let im = re + n;
                        out[re] = -(harmonic as f64) * point.coefficients[im];
                        out[im] = harmonic as f64 * point.coefficients[re];
                    }
                }
                out.iter().all(|v| v.is_finite())
            },
            &mut cancelled,
        );
        if matches!(prepared, Err(ForcedSensitivityError::Cancelled)) {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        if let Some(error) = self.failure.borrow_mut().take() {
            return Err(error.into());
        }
        Ok(PreparedDissipativeAutonomousSensitivity {
            source: self,
            prepared: prepared?,
        })
    }
}

impl<'s, 'a> PreparedDissipativeAutonomousSensitivity<'s, 'a> {
    pub fn prepared(&self) -> &PreparedAutonomousSensitivity {
        &self.prepared
    }

    /// Complete fixed-phase Fourier and solved-frequency derivatives. Source
    /// changes give (J-R)e_p+(J_p-R_p)e+G_p*u+G*u_p-D_e*e_p-D_p; the full
    /// original state chain belongs to the cached bordered factor. The phase
    /// reference stays fixed and parameter_scale only conditions acceptance.
    pub fn solve_parameter(
        &mut self,
        parameter_scale: f64,
        direction: DissipativeAutonomousParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousOrbitSensitivity> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        parameter::with_samples(self.source, direction, |samples| {
            self.prepared.solve_island_direction(
                parameter_scale,
                |_time, state, out| samples.borrow_mut().field(state, out),
                &mut cancelled,
            )
        })
    }

    /// Original two-grid cycle mean and complete coefficient/frequency
    /// covector. The callback's final scalar is its explicit fixed-phase
    /// frequency partial; parameter partials remain separate from this chain.
    pub fn mean_objective(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        observation: impl FnMut(f64, &[f64], &mut f64, &mut [f64], &mut f64) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousMeanObjective> {
        self.prepared
            .mean_objective(controls, observation, cancelled)
            .map_err(Into::into)
    }

    /// Reuse the complete bordered factor for one objective transpose solve.
    /// The independent q_omega partial belongs here, not in a parameter q_p.
    pub fn prepare_adjoint(
        &mut self,
        gradient: &[f64],
        angular_frequency_partial: f64,
        objective_scale: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDissipativeAutonomousAdjoint<'_, 's, 'a>> {
        Ok(PreparedDissipativeAutonomousAdjoint {
            source: self.source,
            prepared: self.prepared.prepare_adjoint(
                gradient,
                angular_frequency_partial,
                objective_scale,
                cancelled,
            )?,
        })
    }
}

/// Complete waveform/frequency adjoint bound to the original loss-aware source.
pub struct PreparedDissipativeAutonomousAdjoint<'p, 's, 'a> {
    source: &'s DissipativeAutonomousHamiltonian<'a>,
    prepared: PreparedAutonomousAdjoint<'p>,
}
impl PreparedDissipativeAutonomousAdjoint<'_, '_, '_> {
    pub fn prepared(&self) -> &PreparedAutonomousAdjoint<'_> {
        &self.prepared
    }

    /// q_p is fixed-state/fixed-frequency. The complete original physical
    /// residual direction is sampled on each cached grid; the dual supplies
    /// all state and unknown-frequency chains without another state solve.
    pub fn gradient_parameter(
        &mut self,
        parameter_scale: f64,
        explicit_derivative: f64,
        direction: DissipativeAutonomousParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousObjectiveGradient> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        parameter::with_samples(self.source, direction, |samples| {
            self.prepared.gradient_island_direction(
                parameter_scale,
                explicit_derivative,
                |_time, state, out| samples.borrow_mut().field(state, out),
                &mut cancelled,
            )
        })
    }
}
