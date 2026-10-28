//! Full continuous cavity/body/neck/DG response, with every duct state retained.

use super::*;
use crate::cavity_modal::shared::harmonic::SharedCavityHarmonicResponse;
use fs_duct::dg::DgHarmonicResponse;

/// Peak exp(-i omega t) response of the complete cavity and downstream duct.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedCavityDgHarmonicResponse {
    /// Physical cavity observations. The nested modal input power includes the
    /// opposing duct pressure; surface/incident source powers remain external.
    pub cavity: SharedCavityHarmonicResponse,
    /// Original DG energy coordinates, physical pressure/flow fields and losses.
    pub duct: DgHarmonicResponse,
    /// Physical inlet pressure from the coupled distributed field [Pa].
    pub duct_inlet_pressure_pa: C64,
    /// Physical flow from the neck INTO that same field [m3/s].
    pub duct_inlet_flow_m3_s: C64,
    /// Mean work of authored surface forces and incident neck pressure [W].
    pub authored_external_power_w: f64,
    /// Mean energy of the complete retained physical system [J].
    pub mean_stored_energy_j: f64,
    /// Independent difference between opposing modal and DG input powers [W].
    pub interface_power_residual_w: f64,
    /// Complete mean loss minus authored external power [W].
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
    /// One direct DG duct with its original closed/open/matched distal boundary
    /// is supported. Branches, rational terminals and terminal channels refuse
    /// here until their complete continuous states are retained; their loads
    /// are never silently omitted. All original diffusion validity bands and
    /// the complete dense harmonic work limit apply. Nothing advances.
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
        if self.load.terminal().is_some()
            || !self.load.branches().is_empty()
            || !self.load.terminal_channels().is_empty()
        {
            return Err(network_invalid("continuous shared-cavity response requires one direct DG duct with no omitted terminal, branch or gap states").into());
        }
        let force = self.cavity.prepare_harmonic_force(
            omega_rad_s,
            surface_force,
            incident_neck_pressure_pa,
            gate,
        )?;
        let mut column = vec![0.0; self.cavity.system.mode_count()];
        column[self.cavity.surface_modes] = self.cavity.source.neck_volume_shape_m2_per_sqrt_kg();
        let loaded = self.cavity.system.harmonic_response_with_duct(
            omega_rad_s,
            &force,
            self.load.duct(),
            self.load.inlet(),
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
