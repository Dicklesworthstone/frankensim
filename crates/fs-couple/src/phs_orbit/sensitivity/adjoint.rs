//! Source-bound objectives through the original forced pHS harmonic balance.

use super::{
    ForcedSensitivityDirection, PeriodicPortInput, PhsForcedSensitivityError,
    PreparedPhsForcedSensitivity, port_direction, storage_direction,
};
use fs_orbit::sensitivity::{ForcedObjectiveGradient, PreparedForcedAdjoint};
use fs_phs::PortHamiltonian;

/// One objective's original physical pHS adjoint. The shared fs-orbit owner
/// reuses the already prepared coupled factor; each source direction only
/// evaluates/project its original physical law and contracts the cached dual.
pub struct PreparedPhsForcedAdjoint<'a, 's> {
    system: &'s PortHamiltonian,
    input: &'a PeriodicPortInput,
    omega: f64,
    harmonics: usize,
    prepared: PreparedForcedAdjoint<'a>,
}

impl<'s> PreparedPhsForcedSensitivity<'s> {
    /// Prepare one complete ordinary real Fourier covector at the fixed forcing
    /// clock. A result from mean_objective/mean_stored_energy supplies this
    /// gradient directly. Objective units condition acceptance only. All
    /// original states, nonlinear couplings and material memories are retained.
    pub fn prepare_adjoint(
        &mut self,
        objective_gradient: &[f64],
        objective_scale: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsForcedAdjoint<'_, 's>, PhsForcedSensitivityError> {
        let omega = self.prepared.angular_frequency();
        let harmonics = self.prepared.harmonics();
        Ok(PreparedPhsForcedAdjoint {
            system: self.system,
            input: &self.input,
            omega,
            harmonics,
            prepared: self.prepared.prepare_adjoint(
                objective_gradient,
                objective_scale,
                cancelled,
            )?,
        })
    }
}

impl PreparedPhsForcedAdjoint<'_, '_> {
    /// Total objective derivative for simultaneous original J/R/G, storage,
    /// Fourier-drive and clock changes. The direct objective partial includes
    /// any explicit forcing-frequency dependence at fixed Fourier coordinates.
    /// Physical state and source directions are assembled by the same original
    /// owners as the waveform derivative; the cached dual performs no new solve.
    pub fn gradient_parameter(
        &mut self,
        chart: ForcedSensitivityDirection,
        explicit_derivative: f64,
        direction: super::PhsForcedParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedObjectiveGradient, PhsForcedSensitivityError> {
        if cancelled() {
            return Err(super::ForcedSensitivityError::Cancelled.into());
        }
        let n = self.system.state_dim();
        super::parameter::with_parameter_direction(
            self.system,
            self.input,
            self.harmonics,
            direction,
            |action| {
                self.prepared.gradient_split_direction(
                    chart,
                    explicit_derivative,
                    action,
                    |point, out| {
                        super::parameter::ode_clock_column(
                            n,
                            chart.angular_frequency_derivative,
                            point,
                            out,
                        )
                    },
                    &mut cancelled,
                )
            },
        )
    }

    /// Shared physical dual diagnostics and complete residual-direction entry.
    pub fn prepared(&self) -> &PreparedForcedAdjoint<'_> {
        &self.prepared
    }

    /// Total objective derivative for original forcing Fourier coefficients.
    /// Clock, port count and bandwidth must match the base. Nonzero prescribed
    /// frequency derivatives refuse on this local-island entry. The explicit
    /// partial is q_p at fixed state coefficients; the waveform chain is added
    /// exactly once by the adjoint. No per-parameter waveform solve is performed.
    pub fn gradient_port_direction(
        &mut self,
        direction: ForcedSensitivityDirection,
        explicit_derivative: f64,
        port_derivative: &PeriodicPortInput,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedObjectiveGradient, PhsForcedSensitivityError> {
        let action = port_direction(self.system, port_derivative, self.omega, self.harmonics)?;
        self.prepared
            .gradient_island_direction(direction, explicit_derivative, action, cancelled)
            .map_err(Into::into)
    }

    /// Total objective derivative for an original storage/material direction at
    /// fixed J/R/G and forcing clock. The callback writes complete partial_p
    /// grad H; it excludes the unknown state chain. Direct objective dependence
    /// on the material belongs in explicit_derivative, even for stored energy.
    pub fn gradient_storage_direction(
        &mut self,
        direction: ForcedSensitivityDirection,
        explicit_derivative: f64,
        effort_derivative: impl FnMut(&[f64], &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedObjectiveGradient, PhsForcedSensitivityError> {
        self.prepared
            .gradient_island_direction(
                direction,
                explicit_derivative,
                storage_direction(self.system, effort_derivative),
                cancelled,
            )
            .map_err(Into::into)
    }
}
