//! Local SGS on the original native reservoir geometry and transaction.
//!
//! Momentum and heat consume one accepted force-centered velocity gradient.
//! An on-site inlet is differentiated with its original one-sided stencil;
//! a pressure face retains its native zero-normal-gradient donor convention.
//! Neither is interpreted as an impermeable halfway wall.

use super::*;
use crate::sgs3::{SgsConfig, SgsLattice3, SgsReport, SgsStepError, SgsWorkspace3};

mod forcing;
mod thermal;
pub use thermal::{
    OpenThermalSgsConfig3, OpenThermalSgsError3, OpenThermalSgsStep3, OpenThermalSgsWorkspace3,
};

/// Original macroscopic/collision buffers plus the actual prepared eddy
/// viscosity, in native x-fast order. All fields are scratch, valid only after
/// a successful preparation. Reuse never changes the molecular collision law.
pub struct OpenSgsWorkspace3<C> {
    local: SgsWorkspace3<C>,
    eddy: Vec<f64>,
}

impl<C> OpenSgsWorkspace3<C> {
    pub fn collisions(&self) -> &[C] {
        self.local.collisions()
    }
    pub fn eddy_viscosities(&self) -> &[f64] {
        &self.eddy
    }
    /// Actual retained numerical Vec capacities, excluding fixed metadata.
    pub fn workspace_bytes(&self) -> usize {
        self.local.velocity.capacity() * core::mem::size_of::<[f64; 3]>()
            + (self.local.density.capacity() + self.eddy.capacity()) * core::mem::size_of::<f64>()
            + self.local.collisions.capacity() * core::mem::size_of::<C>()
    }
}

/// Native reservoir exchange and the incoming modeled shear diagnostics.
#[derive(Clone, Debug)]
pub struct OpenSgsStep3 {
    pub momentum: OpenBoundaryMomentumTransport3,
    pub sgs: SgsReport,
}

impl<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    pub fn required_sgs_workspace_bytes(&self) -> Result<usize, SgsStepError<L::Error>> {
        let per_cell = 5usize
            .checked_mul(core::mem::size_of::<f64>())
            .and_then(|n| n.checked_add(core::mem::size_of::<L::Collision>()))
            .ok_or(GridError3::Input("open SGS workspace size overflow"))?;
        self.grid
            .f
            .len()
            .checked_mul(per_cell)
            .ok_or_else(|| GridError3::Input("open SGS workspace size overflow").into())
    }

    /// Allocate once under an explicit numerical payload budget. This adds
    /// five f64 values and one native collision value per grid cell. Existing
    /// flow/thermal populations retain their original independent budgets.
    pub fn sgs_workspace(
        &self,
        maximum_workspace_bytes: usize,
    ) -> Result<OpenSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        if self.required_sgs_workspace_bytes()? > maximum_workspace_bytes {
            return Err(GridError3::Input("open SGS workspace byte budget exceeded").into());
        }
        let local = self.grid.sgs_workspace()?;
        let mut eddy = Vec::new();
        eddy.try_reserve_exact(self.grid.f.len())
            .map_err(|_| GridError3::Allocation)?;
        eddy.resize(self.grid.f.len(), 0.0);
        let workspace = OpenSgsWorkspace3 { local, eddy };
        if workspace.workspace_bytes() > maximum_workspace_bytes {
            return Err(GridError3::Input("open SGS workspace capacity budget exceeded").into());
        }
        Ok(workspace)
    }

    /// Prepare Smagorinsky, WALE or Vreman using the original open-flow full
    /// gradient. Accepted velocities include the held force's half-step shift.
    /// Interior/periodic and stationary-wall stencils are unchanged. Velocity
    /// reservoir slopes are quadratic one-sided (linear with a solid second
    /// donor); pressure faces retain their native zero normal derivative.
    ///
    /// Optional collisions are fresh MOLECULAR transport, never a preceding
    /// augmented field. The original collision family, bulk and higher modes
    /// are retained; the exact zero model preserves collision bits. Model and
    /// molecular-ratio admissions remain the shared exploratory policy.
    pub fn prepare_sgs(
        &self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut OpenSgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        self.prepare_open_sgs(config, molecular, workspace, None, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })
    }

    pub(super) fn prepare_open_sgs(
        &self,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut OpenSgsWorkspace3<L::Collision>,
        forces: Option<&[[f64; 3]]>,
        checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        if workspace.eddy.len() != self.grid.f.len() {
            return Err(GridError3::Input("open SGS eddy field length mismatch").into());
        }
        self.grid.prepare_sgs_with_gradient(
            config,
            molecular,
            &mut workspace.local,
            forces,
            |cell, velocity| self.rheology_gradient(cell, velocity),
            |cell, eddy| workspace.eddy[cell] = eddy,
            checkpoint,
        )
    }

    /// One original native step with prepared turbulent shear and its actual
    /// reservoir mass/momentum/work receipt. No model persists into ordinary
    /// later steps. All failure/cancellation publication remains native.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut OpenSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
    ) -> Result<OpenSgsStep3, SgsStepError<L::Error>> {
        let sgs = self.prepare_sgs(cx, config, molecular, workspace)?;
        let momentum = self.step_serial(cx, Some(workspace.collisions()), origin)?;
        Ok(OpenSgsStep3 { momentum, sgs })
    }

    /// Pooled native collision/streaming with the same serial coefficient
    /// preparation, fixed reduction order, reservoir rules and cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled_with_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut OpenSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
    ) -> Result<OpenSgsStep3, SgsStepError<L::Error>> {
        let sgs = self.prepare_open_sgs(config, molecular, workspace, None, || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        })?;
        let momentum = self.step_pooled(pool, gate, Some(workspace.collisions()), origin)?;
        Ok(OpenSgsStep3 { momentum, sgs })
    }
}

#[cfg(test)]
mod tests;
