//! SGS transport on the retained-material three-dimensional thermal slab.

use super::{HeatBalance, HeatSourceLimits, MaterialSlab, ThermalLbm3, invalid};
use crate::thermal::sources::sgs3::{
    ThermalSgsConfig, ThermalSgsError, ThermalSgsReport, ThermalSgsWorkspace3, step_sgs_mapped,
};
use fs_exec::CancelGate;

/// Accepted sensible heat in SI units and the exact applied SGS diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialThermalSgsStep {
    pub heat: HeatBalance,
    pub sgs: ThermalSgsReport,
}

impl MaterialSlab<ThermalLbm3> {
    /// Advance prescribed power [W/m3] and external force [N/m3] with the
    /// shared SGS closure and explicit turbulent Prandtl number. Filter width
    /// remains explicitly in lattice cells; `nu_t * dx^2/dt` gives SI viscosity.
    /// Density, heat capacity, molecular transport, dx and dt retain the same
    /// material snapshot and normalized temperature map as `step_with_power_density`.
    ///
    /// The original SI source conversion and independently measured wall heat
    /// are checked before the existing thermal owner publishes its candidate.
    /// This frozen-property Boussinesq step adds no viscous/dissipation heating,
    /// pressure work, latent heat or kinetic-energy account. For authored local
    /// molecular fields, use `ThermalLbm3::step_with_sources_and_sgs` explicitly.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_power_density_and_sgs(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        config: ThermalSgsConfig,
        workspace: &mut ThermalSgsWorkspace3,
        gate: &CancelGate,
    ) -> Result<MaterialThermalSgsStep, ThermalSgsError> {
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
        let (heat, sgs) = step_sgs_mapped(
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
        Ok(MaterialThermalSgsStep { heat, sgs })
    }
}
