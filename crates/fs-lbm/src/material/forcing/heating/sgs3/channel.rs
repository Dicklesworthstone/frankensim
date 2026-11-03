//! Dynamic plane-channel transport with the retained material's SI heat account.

use super::*;
use crate::thermal::sources::sgs3::channel::{
    ThermalChannelSgsConfig, ThermalChannelSgsReport, ThermalChannelSgsWorkspace3,
    step_channel_sgs_mapped,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialThermalChannelSgsStep {
    pub heat: HeatBalance,
    pub sgs: ThermalChannelSgsReport,
}

impl MaterialSlab<ThermalLbm3> {
    /// Advance prescribed power [W/m3] and external force [N/m3] with
    /// independently fitted momentum and heat on the existing halfway slab.
    /// The same frozen density, heat capacity, molecular transport, normalized
    /// temperature map and dx/dt conversion own both the source and wall heat.
    /// Filter spacing is one lattice cell; no turbulent Prandtl is prescribed.
    ///
    /// Original SI balance and endpoint checks precede publication. This is
    /// sensible-temperature Boussinesq transport with resolved wall layers,
    /// without viscous heating, pressure work or total-energy qualification.
    /// Authored local molecular rates use the lattice-unit entry directly.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_power_density_and_channel_sgs(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        config: ThermalChannelSgsConfig,
        workspace: &mut ThermalChannelSgsWorkspace3,
        gate: &CancelGate,
    ) -> Result<MaterialThermalChannelSgsStep, ThermalSgsError> {
        if gate.is_requested() {
            return Err(ThermalSgsError::Cancelled);
        }
        if self.lattice.t_bottom != 1.0 || self.lattice.t_top != 0.0 {
            return Err(invalid(
                "material heating refuses independently changed normalized wall temperatures",
            )
            .into());
        }
        let input = self.heat_inputs(
            power_w_m3,
            forces_n_m3,
            self.lattice.gpop.len(),
            self.force_units.cell_size_m,
            limits,
        )?;
        let (heat, sgs) = step_channel_sgs_mapped(
            &mut self.lattice,
            &input.heat,
            &input.force,
            config,
            None,
            workspace,
            input.limits,
            gate,
            |balance| input.balance(balance),
        )?;
        Ok(MaterialThermalChannelSgsStep { heat, sgs })
    }
}
