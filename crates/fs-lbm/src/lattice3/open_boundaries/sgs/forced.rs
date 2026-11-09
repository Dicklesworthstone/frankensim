//! Local turbulent heat on the original open-flow force-history owner.
//!
//! Incoming gradients use F_n/2, while the native transaction retains
//! F_(n+1) for endpoint velocity and prescribed heat-source half-kicks.

use super::{
    EquilibriumLattice3, GridError3, OpenThermalSgsConfig3, OpenThermalSgsError3,
    OpenThermalSgsStep3, OpenThermalSgsWorkspace3, SgsLattice3, SgsStepError,
};
use crate::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;
use crate::thermal::open::ThermalOpenField3;
use fs_exec::Cx;

impl<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize> ForcedOpenBoundaryGrid3<L, Q> {
    /// The original reusable local momentum/heat workspace budget.
    pub fn required_thermal_sgs_workspace_bytes(&self) -> Result<usize, SgsStepError<L::Error>> {
        self.inner.required_thermal_sgs_workspace_bytes()
    }

    /// Reuse the original six-f64-plus-collision storage at each cell.
    /// This workspace holds no separate force, thermal or boundary history.
    pub fn thermal_sgs_workspace(
        &self,
        maximum_workspace_bytes: usize,
    ) -> Result<OpenThermalSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        self.inner.thermal_sgs_workspace(maximum_workspace_bytes)
    }

    /// Hold the accepted spatial body force while transporting SGS momentum
    /// and heat through the original open-grid transaction.
    ///
    /// The FULL incoming physical velocity includes each cell's accepted
    /// F_n/2 shift. Using raw momentum would turn a spatial force gradient
    /// into fictitious strain. Reservoir and prescribed halfway-wall stencils
    /// retain their native geometry and actual wall velocities.
    ///
    /// One nu_t supplies both the native shear collision and
    /// alpha_t=nu_t/Pr_t. Molecular collision choices, thermal diffusivity,
    /// coefficient-ratio admissions and exact disabled controls are shared
    /// with the ordinary open thermal SGS path. Supply fresh molecular
    /// coefficients; previous augmented coefficients must not be fed back.
    ///
    /// Sources, reservoir heat, finite stationary walls, wall work, scalar
    /// flux history and counters publish only after the native transaction
    /// succeeds. Scratch may change on refusal. Mechanical work is not added
    /// to the fixed-reference-capacity scalar heat account; no total-energy,
    /// wall-function or turbulence-validation claim is supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_thermal_sgs(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        config: OpenThermalSgsConfig3,
        molecular_collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        molecular_thermal_relaxation: &[f64],
        workspace: &mut OpenThermalSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
    ) -> Result<OpenThermalSgsStep3, OpenThermalSgsError3<L::Error>> {
        let (sgs, maximum_eddy_diffusivity_lu, maximum_diffusivity_to_molecular) =
            self.inner.prepare_open_thermal_sgs(
                config,
                molecular_collisions,
                molecular_thermal_relaxation,
                workspace,
                Some(&self.forces),
                || cx.checkpoint().map_err(|_| GridError3::Cancelled),
            )?;
        let transport = self.step_serial_with_thermal(
            thermal,
            cx,
            Some(workspace.momentum_collisions()),
            heat_per_step,
            Some(workspace.thermal_relaxation()),
            origin,
        )?;
        Ok(OpenThermalSgsStep3 {
            transport,
            sgs,
            maximum_eddy_diffusivity_lu,
            maximum_diffusivity_to_molecular,
        })
    }

    /// Advance an authored force endpoint with the same incoming turbulent
    /// momentum/heat coefficients and original joint publication.
    ///
    /// Coefficients are fitted ONLY to the accepted populations and F_n.
    /// The proposed F_(n+1) changes neither this gradient nor its Prandtl
    /// heat closure. It enters the original endpoint velocity and final
    /// prescribed-source half-kick. No extra force reframe or impulse occurs.
    ///
    /// On success forces_at_end receives the old accepted force allocation.
    /// On any refusal/cancellation it and all accepted physical state retain
    /// their original values. The SGS workspace introduces no step allocation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_thermal_sgs(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        config: OpenThermalSgsConfig3,
        molecular_collisions: Option<&[L::Collision]>,
        forces_at_end: &mut Vec<[f64; 3]>,
        heat_per_step: &[f64],
        molecular_thermal_relaxation: &[f64],
        workspace: &mut OpenThermalSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
    ) -> Result<OpenThermalSgsStep3, OpenThermalSgsError3<L::Error>> {
        let (sgs, maximum_eddy_diffusivity_lu, maximum_diffusivity_to_molecular) =
            self.inner.prepare_open_thermal_sgs(
                config,
                molecular_collisions,
                molecular_thermal_relaxation,
                workspace,
                Some(&self.forces),
                || cx.checkpoint().map_err(|_| GridError3::Cancelled),
            )?;
        let transport = self.advance_serial_with_thermal(
            thermal,
            cx,
            Some(workspace.momentum_collisions()),
            forces_at_end,
            heat_per_step,
            Some(workspace.thermal_relaxation()),
            origin,
        )?;
        Ok(OpenThermalSgsStep3 {
            transport,
            sgs,
            maximum_eddy_diffusivity_lu,
            maximum_diffusivity_to_molecular,
        })
    }
}

#[cfg(test)]
mod tests;
