//! Local modeled stress and heat in one native buoyant-flow transaction.

use super::*;
use crate::lattice3::open_boundaries::sgs::{
    OpenThermalSgsConfig3, OpenThermalSgsError3, OpenThermalSgsWorkspace3,
};
use crate::sgs3::{SgsLattice3, SgsReport, SgsStepError};

/// Closure preparation or the original joint flow/force/heat owner refused.
#[derive(Clone, Debug)]
pub enum BoussinesqOpenSgsError3<E> {
    /// Shared local momentum/thermal coefficient admission.
    Sgs(OpenThermalSgsError3<E>),
    /// Native buoyancy, flow, source, thermal or finite-wall admission.
    Transport(BoussinesqOpenError3<E>),
}
impl<E> From<OpenThermalSgsError3<E>> for BoussinesqOpenSgsError3<E> {
    fn from(error: OpenThermalSgsError3<E>) -> Self {
        Self::Sgs(error)
    }
}
impl<E> From<BoussinesqOpenError3<E>> for BoussinesqOpenSgsError3<E> {
    fn from(error: BoussinesqOpenError3<E>) -> Self {
        Self::Transport(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for BoussinesqOpenSgsError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native buoyant open SGS: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for BoussinesqOpenSgsError3<E> {}

/// One jointly accepted native temperature, force and turbulent flow endpoint.
#[derive(Clone, Debug)]
pub struct BoussinesqOpenSgsStep3 {
    /// Original body/background/buoyancy impulses and reservoir/wall heat.
    pub transport: BoussinesqOpenStep3,
    /// Incoming modeled/resolved shear transfer, never an added heat source.
    pub sgs: SgsReport,
    /// Largest represented nu_t/Pr_t over fluid cells [lattice units].
    pub maximum_eddy_diffusivity_lu: f64,
    /// Largest admitted alpha_t/alpha_molecular over fluid cells.
    pub maximum_diffusivity_to_molecular: f64,
}

impl<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize> BoussinesqOpenGrid3<L, Q> {
    /// Required payload of the existing reusable native SGS/thermal workspace.
    pub fn required_thermal_sgs_workspace_bytes(&self) -> Result<usize, SgsStepError<L::Error>> {
        self.flow.inner.required_thermal_sgs_workspace_bytes()
    }

    /// Reserve six f64 values plus one native collision per original grid cell.
    /// The original population, thermal and force owners keep their own caps.
    pub fn thermal_sgs_workspace(
        &self,
        maximum_workspace_bytes: usize,
    ) -> Result<OpenThermalSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        self.flow.inner.thermal_sgs_workspace(maximum_workspace_bytes)
    }

    /// Advance heat and buoyant momentum with Smagorinsky, WALE or Vreman stress.
    ///
    /// The shared native gradient samples the accepted populations with their
    /// OWN spatial Guo half-force. Its nu_t augments the supplied fresh molecular
    /// collisions and alpha_t=nu_t/Pr_t augments the fresh molecular D3Q7 tau.
    /// Neither a previous augmented field nor a future buoyancy/background force
    /// may replace those incoming physical coefficients and velocity samples.
    ///
    /// The original owner then transports heat, computes endpoint buoyancy,
    /// collides with opening forces and reconstructs with endpoint forces. Its
    /// scalar sources, finite conducting walls, reservoir/backflow exchanges,
    /// force consistency and final publication gates remain unchanged. There
    /// is no second heat solve or deposition of modeled/force work as heat.
    ///
    /// All preparation is scratch only. Any late refusal preserves populations,
    /// thermal history, finite walls, force epochs, time and accepted receipts.
    /// This remains explicit fixed-reference-capacity Boussinesq transport with
    /// the original split and exploratory closure limits; it adds no variable-
    /// density EOS, total-energy equation or new temporal-order claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: OpenThermalSgsConfig3,
        molecular_collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        molecular_thermal_relaxation: &[f64],
        workspace: &mut OpenThermalSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
    ) -> Result<BoussinesqOpenSgsStep3, BoussinesqOpenSgsError3<L::Error>> {
        self.advance_sgs_with_poll(
            None,
            config,
            molecular_collisions,
            heat_per_step,
            molecular_thermal_relaxation,
            workspace,
            origin,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// The same turbulent buoyant transaction with a prescribed NEXT mechanical
    /// background. SGS still uses the retained opening total force. Only the
    /// original successful publication swaps `background_at_end` with the old
    /// background, so a refused attempt leaves the caller's buffer unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: OpenThermalSgsConfig3,
        molecular_collisions: Option<&[L::Collision]>,
        background_at_end: &mut Vec<[f64; 3]>,
        heat_per_step: &[f64],
        molecular_thermal_relaxation: &[f64],
        workspace: &mut OpenThermalSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
    ) -> Result<BoussinesqOpenSgsStep3, BoussinesqOpenSgsError3<L::Error>> {
        self.advance_sgs_with_poll(
            Some(background_at_end),
            config,
            molecular_collisions,
            heat_per_step,
            molecular_thermal_relaxation,
            workspace,
            origin,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_sgs_with_poll(
        &mut self,
        background_at_end: Option<&mut Vec<[f64; 3]>>,
        config: OpenThermalSgsConfig3,
        molecular_collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        molecular_thermal_relaxation: &[f64],
        workspace: &mut OpenThermalSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<BoussinesqOpenSgsStep3, BoussinesqOpenSgsError3<L::Error>> {
        let (sgs, maximum_eddy_diffusivity_lu, maximum_diffusivity_to_molecular) = self
            .flow
            .inner
            .prepare_open_thermal_sgs(
                config,
                molecular_collisions,
                molecular_thermal_relaxation,
                workspace,
                Some(&self.flow.forces),
                &mut *poll,
            )?;
        let transport = self.advance_with_diffusion_poll(
            background_at_end,
            Some(workspace.momentum_collisions()),
            heat_per_step,
            Some(workspace.thermal_relaxation()),
            None,
            origin,
            poll,
        )?;
        Ok(BoussinesqOpenSgsStep3 {
            transport,
            sgs,
            maximum_eddy_diffusivity_lu,
            maximum_diffusivity_to_molecular,
        })
    }
}

#[cfg(test)]
mod tests;
