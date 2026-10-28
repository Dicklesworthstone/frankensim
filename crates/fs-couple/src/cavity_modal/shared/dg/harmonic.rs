//! Full continuous cavity/body/neck/DG response, with every duct state retained.

use super::*;
use crate::cavity_modal::shared::harmonic::SharedCavityHarmonicResponse;
pub use crate::render::schedule::force::coupled::harmonic::dg::{
    DgHarmonicBranchResponse, DgHarmonicTerminalResponse,
};
use fs_duct::dg::DgHarmonicResponse;

/// Peak exp(-i omega t) response of the complete cavity and downstream duct.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedCavityDgHarmonicResponse {
    /// Physical cavity observations. The nested modal input power includes the
    /// opposing duct pressure; surface/incident source powers remain external.
    pub cavity: SharedCavityHarmonicResponse,
    /// Original DG energy coordinates, physical pressure/flow fields and losses.
    pub duct: DgHarmonicResponse,
    /// Complete finite side branches and their original junction exchanges.
    pub branches: Vec<DgHarmonicBranchResponse>,
    /// Original terminal/radiation state phasors and channel exchanges.
    pub terminal: Option<DgHarmonicTerminalResponse>,
    /// Signed mean power into the original terminal, outside known storage [W].
    pub terminal_power_w: f64,
    /// Physical inlet pressure from the coupled distributed field [Pa].
    pub duct_inlet_pressure_pa: C64,
    /// Physical flow from the neck INTO that same field [m3/s].
    pub duct_inlet_flow_m3_s: C64,
    /// Mean work of authored surface forces and incident neck pressure [W].
    pub authored_external_power_w: f64,
    /// Known cavity and duct energy [J]; generic terminal storage is unavailable.
    pub mean_stored_energy_j: f64,
    /// Independent balance of modal, duct and terminal interface powers [W].
    pub interface_power_residual_w: f64,
    /// Known mean loss plus terminal exchange minus authored external power [W].
    pub power_residual_w: f64,
    /// Infinity-norm backward residual of the complete continuous equations.
    pub solve_relative_residual: f64,
}
impl SharedCavityDgHarmonicResponse {
    /// Mechanical/material loss plus actual DG wall loss and matched absorption.
    pub fn total_dissipation_w(&self) -> f64 {
        self.cavity.network.total_dissipation_w()
            + self.duct.wall_dissipation_w
            + self.duct.matched_absorption_w
            + self.branches.iter().fold(0.0, |sum, branch| {
                sum + branch.duct.wall_dissipation_w + branch.duct.matched_absorption_w
            })
    }
}

impl SharedCavityDg {
    /// Solve the complete CONTINUOUS response of this cavity and original duct.
    ///
    /// This includes the duct's distributed resonances and actual viscothermal
    /// or compliant-wall memory. It does not invert an isolated duct impedance,
    /// so a bare duct pole can become a finite joint antiresonance. Surface
    /// forces and inward incident neck pressure retain the original conventions.
    /// The result is independent of accepted transients and of the time backend's
    /// coupling clock. It is not the finite-step transfer of the mean-port path.
    ///
    /// The original DG bore with its closed/open/matched distal boundary or its
    /// complete rational terminal/radiation matrix is supported, including every
    /// finite side duct and its zero-storage pressure-continuous junction. All
    /// original pole states, mutual channels and frequency limits are retained. Singular
    /// terminal feedthrough remains a valid retained-state frequency problem.
    /// Finite gap channels and driven mouths refuse until their states
    /// are retained too. The complete dense work limit applies. Nothing advances.
    pub fn harmonic_response(
        &self,
        omega_rad_s: f64,
        surface_force: &[C64],
        incident_neck_pressure_pa: C64,
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<SharedCavityDgHarmonicResponse, SharedCavityDgError> {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled.into());
        }
        let force = self.cavity.prepare_harmonic_force(
            omega_rad_s,
            surface_force,
            incident_neck_pressure_pa,
            gate,
        )?;
        let mut column = vec![0.0; self.cavity.system.mode_count()];
        column[self.cavity.surface_modes] = self.cavity.source.neck_volume_shape_m2_per_sqrt_kg();
        let loaded = self.cavity.system.harmonic_response_with_load(
            omega_rad_s,
            &force,
            &self.load,
            &column,
            maximum_work,
            gate,
        )?;
        let cavity = self.cavity.observe_harmonic_network(
            omega_rad_s,
            surface_force,
            incident_neck_pressure_pa,
            loaded.modal,
            gate,
        )?;
        Ok(SharedCavityDgHarmonicResponse {
            cavity,
            duct: loaded.duct,
            branches: loaded.branches,
            terminal: loaded.terminal,
            terminal_power_w: loaded.terminal_power_w,
            duct_inlet_pressure_pa: loaded.inlet_pressure_pa,
            duct_inlet_flow_m3_s: loaded.inlet_flow_m3_s,
            authored_external_power_w: loaded.authored_external_power_w,
            mean_stored_energy_j: loaded.mean_stored_energy_j,
            interface_power_residual_w: loaded.interface_power_residual_w,
            power_residual_w: loaded.power_residual_w,
            solve_relative_residual: loaded.solve_relative_residual,
        })
    }
}

#[cfg(test)]
mod tests;
