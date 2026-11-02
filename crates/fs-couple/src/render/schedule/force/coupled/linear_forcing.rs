//! Continuous piecewise-linear forcing of the complete bilateral network.
//!
//! All original modes, signed damping columns, material memories and radiation
//! states evolve together. Preparation is tied to the unchanged coefficients
//! and sample period; each warm step reuses the exact flow and positive heat
//! caches. The existing physical candidate and port-account allocations remain.
use super::*;

/// Explicit cold dimension, memory and multiply-term allowance for the full
/// trajectory and every positive loss channel. The native source-square owner
/// limits twice the extended state dimension to 256; no states are discarded.
pub use fs_phs::linear_exact::DenseExponentialBudget as LinearForceBudget;

/// A fixed-coefficient coupled system with cached continuous force ramps.
///
/// Read-only access preserves the original states and diagnostics. Consume the
/// wrapper to change coefficients or attachments, then prepare again. This
/// prevents an old transition or heat factor from surviving a material update.
pub struct LinearForceModalSystem {
    system: CoupledModalSystem,
    prepared: exact::linear_forcing::PreparedLinearForcing,
}

impl CoupledModalSystem {
    /// Prepare exact linear force records at this system's existing clock.
    /// Accepted mechanical, Maxwell and acoustic histories are preserved, even
    /// when this conversion follows earlier held-force steps. Harmonic-only
    /// constructions and over-budget complete operators refuse explicitly.
    pub fn into_linear_forcing(
        self,
        budget: LinearForceBudget,
        gate: &CancelGate,
    ) -> Result<LinearForceModalSystem, ModalCouplingError> {
        let prepared = exact::linear_forcing::PreparedLinearForcing::new(&self, budget, gate)?;
        Ok(LinearForceModalSystem {
            system: self,
            prepared,
        })
    }
}

impl LinearForceModalSystem {
    /// Original accepted state, physical accounts and port observations.
    #[must_use]
    pub const fn system(&self) -> &CoupledModalSystem {
        &self.system
    }

    /// Retain every accepted state and clock while releasing fixed coefficients.
    #[must_use]
    pub fn into_system(self) -> CoupledModalSystem {
        self.system
    }

    /// Advance `F(t) = start + (end-start)*t/dt` over one complete sample.
    /// Forces use the original flattened mass-normalized component order.
    ///
    /// Work is the signed integral of the actual force times velocity. Each
    /// component, connection, material and radiation heat channel independently
    /// integrates its original positive source squares on that same trajectory.
    /// Original energy, pressure, reaction, material-domain and cancellation
    /// gates run before any state or account publishes. This API does not change
    /// nonlinear contact's separately solved held reaction law.
    pub fn step(
        &mut self,
        start: &[f64],
        end: &[f64],
        gate: &CancelGate,
    ) -> Result<&CoupledModalFrame, ModalCouplingError> {
        self.prepared.stage(&mut self.system, start, end, gate)?;
        self.system.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        self.system.publish_staged();
        Ok(&self.system.frame)
    }
}

#[cfg(test)]
mod tests;
