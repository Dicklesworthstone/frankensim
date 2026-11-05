//! Original source rheology composed with native local or dynamic SGS.
//!
//! Material samples describe molecular transport at the accepted incoming
//! resolved strain. SGS is then prepared against those exact molecular rates;
//! its modeled viscosity is never written into a material response. Both use
//! the same accepted Guo force frame. A next-time force enters only the native
//! transaction, never the incoming material query or turbulence fit.
//!
//! The existing source and SGS workspaces own all scratch. Source epochs and
//! the applied collision field publish only after the original step succeeds.
//! This explicit resolved-shear closure does not model correlations between
//! unresolved strain and nonlinear molecular viscosity, or add a heat equation.

use super::*;
use crate::sgs3::dynamic::DynamicSgsConfig;
use crate::sgs3::dynamic::channel::ChannelDynamicSgsConfig;
use crate::sgs3::{
    ChannelDynamicSgsWorkspace3, DynamicSgsStepReport, DynamicSgsWorkspace3, SgsConfig,
    SgsLattice3, SgsReport, SgsStepError, SgsWorkspace3,
};

mod moving_walls;
pub use moving_walls::ChannelRheologySgs3;

/// One explicit closure and its existing reusable scratch, borrowed for this
/// step only. Dynamic variants retain their original periodic/channel topology,
/// numerical capacity admission and exploration-only policy.
pub enum RheologySgs3<'a, C> {
    Local {
        config: SgsConfig,
        workspace: &'a mut SgsWorkspace3<C>,
    },
    PeriodicDynamic {
        config: DynamicSgsConfig,
        workspace: &'a mut DynamicSgsWorkspace3<C>,
    },
    ChannelDynamic {
        config: ChannelDynamicSgsConfig,
        workspace: &'a mut ChannelDynamicSgsWorkspace3<C>,
    },
}

/// Diagnostics from the same incoming source/velocity epoch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RheologySgsReport3 {
    Local(SgsReport),
    Dynamic(DynamicSgsStepReport),
}

/// Published source epoch and turbulence diagnostics from one accepted step.
#[derive(Debug)]
pub struct RheologySgsStep3<'a, C, F> {
    /// Samples and their relaxation times remain molecular. The collision
    /// vector is the actual applied molecular-plus-SGS field.
    pub rheology: &'a RheologyStep3<C, F>,
    pub sgs: RheologySgsReport3,
}

/// Original source/scaling refusals and SGS/native-step failures stay distinct.
#[derive(Debug)]
pub enum RheologySgsError3<E> {
    Rheology(RheologyError3<E>),
    Sgs(SgsStepError<E>),
}

impl<E> From<RheologyError3<E>> for RheologySgsError3<E> {
    fn from(error: RheologyError3<E>) -> Self {
        Self::Rheology(error)
    }
}
impl<E> From<SgsStepError<E>> for RheologySgsError3<E> {
    fn from(error: SgsStepError<E>) -> Self {
        Self::Sgs(error)
    }
}
impl<E> From<GridError3<E>> for RheologySgsError3<E> {
    fn from(error: GridError3<E>) -> Self {
        Self::Sgs(SgsStepError::Grid(error))
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RheologySgsError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-RHEOLOGY-SGS3: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for RheologySgsError3<E> {}

fn prepare_sgs<'a, L: SgsLattice3<Q>, const Q: usize>(
    grid: &Grid3<L, Q>,
    forces: Option<&[[f64; 3]]>,
    collisions: &mut [L::Collision],
    sgs: RheologySgs3<'a, L::Collision>,
    mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
) -> Result<(RheologySgsReport3, &'a [L::Collision]), SgsStepError<L::Error>> {
    let (report, applied) = match sgs {
        RheologySgs3::Local { config, workspace } => {
            let report = grid.prepare_sgs_with_forces(
                config,
                Some(&*collisions),
                workspace,
                forces,
                &mut checkpoint,
            )?;
            (RheologySgsReport3::Local(report), workspace.collisions())
        }
        RheologySgs3::PeriodicDynamic { config, workspace } => {
            let report = grid.prepare_dynamic_sgs_with_forces(
                config,
                Some(&*collisions),
                workspace,
                forces,
                &mut checkpoint,
            )?;
            (RheologySgsReport3::Dynamic(report), workspace.collisions())
        }
        RheologySgs3::ChannelDynamic { config, workspace } => {
            let report = grid.prepare_channel_dynamic_sgs_with_forces(
                config,
                Some(&*collisions),
                workspace,
                forces,
                &mut checkpoint,
            )?;
            (RheologySgsReport3::Dynamic(report), workspace.collisions())
        }
    };
    // Both preparations admit exactly the native cell count. Retain the
    // applied field in trial scratch BEFORE native publication, with bounded
    // cancellation. The source samples themselves keep molecular tau.
    for (trial, field) in collisions.chunks_mut(64).zip(applied.chunks(64)) {
        checkpoint()?;
        trial.copy_from_slice(field);
    }
    checkpoint()?;
    Ok((report, applied))
}

// No fallible or O(N) work follows native population/clock publication.
fn publish<'w, C, F>(
    workspace: &'w mut RheologyWorkspace3<C, F>,
    sgs: RheologySgsReport3,
) -> RheologySgsStep3<'w, C, F> {
    RheologySgsStep3 {
        rheology: workspace.publish(),
        sgs,
    }
}

impl<L: SgsLattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Resolve the original molecular source at the current physical strain,
    /// fit the selected SGS against that molecular field, and advance once.
    ///
    /// Source producers are queried once per fluid cell, serially. They must
    /// retain the reference density and must not advance external histories
    /// during a trial. The actual source mu supplies every SGS ratio bound and
    /// molecular dissipation term; constructor viscosity is never substituted.
    /// Source samples retain molecular tau, while the returned collision field
    /// includes SGS. Each following call resolves molecular transport afresh.
    ///
    /// Existing stationary/periodic walls, source-domain checks and independent
    /// nonshear collision rates remain authoritative. The original workspaces
    /// are reused without new adapter-owned allocation. Heat, unresolved
    /// rheological stress correlations and turbulence qualification are absent.
    pub fn step_serial_with_material_rheology_and_sgs<
        'w,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: RheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologySgsStep3<'w, L::Collision, F>, RheologySgsError3<L::Error>> {
        prepare_workspace(self, reference, None, workspace, resolve, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })?;
        let (report, applied) =
            prepare_sgs(self, None, &mut workspace.trial.collisions, sgs, || {
                cx.checkpoint().map_err(|_| GridError3::Cancelled)
            })?;
        self.step_serial_with_collisions(cx, applied)?;
        Ok(publish(workspace, report))
    }

    /// The same source and SGS preparation followed by the original pooled
    /// collision/stream transaction. Source objects need neither Clone nor Send.
    pub fn step_pooled_with_material_rheology_and_sgs<
        'w,
        P: KernelRunner,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        reference: &MaterialScaling<'_, R>,
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: RheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologySgsStep3<'w, L::Collision, F>, RheologySgsError3<L::Error>> {
        prepare_workspace(self, reference, None, workspace, resolve, || {
            gate_poll(gate)
        })?;
        let (report, applied) =
            prepare_sgs(self, None, &mut workspace.trial.collisions, sgs, || {
                gate_poll(gate)
            })?;
        self.step_pooled_with_collisions(pool, gate, applied)?;
        Ok(publish(workspace, report))
    }
}

impl<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Compose source rheology and SGS with the accepted spatial force F_n.
    /// The same force is held for the original native step.
    pub fn step_serial_with_material_rheology_and_sgs<
        'w,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: RheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologySgsStep3<'w, L::Collision, F>, RheologySgsError3<L::Error>> {
        prepare_workspace(
            &self.grid,
            reference,
            Some(self.forces()),
            workspace,
            resolve,
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )?;
        let (report, applied) = prepare_sgs(
            &self.grid,
            Some(self.forces()),
            &mut workspace.trial.collisions,
            sgs,
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )?;
        self.step_serial(cx, applied)?;
        Ok(publish(workspace, report))
    }

    /// Resolve both molecular and modeled viscosity at F_n, then pass F_(n+1)
    /// only to the original trapezoidal force advance. A failed endpoint keeps
    /// populations, force history, clock and the previous source report.
    pub fn advance_serial_with_material_rheology_and_sgs<
        'w,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        forces_at_end: Vec<[f64; 3]>,
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: RheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologySgsStep3<'w, L::Collision, F>, RheologySgsError3<L::Error>> {
        prepare_workspace(
            &self.grid,
            reference,
            Some(self.forces()),
            workspace,
            resolve,
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )?;
        let (report, applied) = prepare_sgs(
            &self.grid,
            Some(self.forces()),
            &mut workspace.trial.collisions,
            sgs,
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )?;
        self.advance_serial(cx, applied, forces_at_end)?;
        Ok(publish(workspace, report))
    }

    /// Pooled composition retaining the actual incoming spatial force frame.
    pub fn step_pooled_with_material_rheology_and_sgs<
        'w,
        P: KernelRunner,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        reference: &MaterialScaling<'_, R>,
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: RheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologySgsStep3<'w, L::Collision, F>, RheologySgsError3<L::Error>> {
        prepare_workspace(
            &self.grid,
            reference,
            Some(self.forces()),
            workspace,
            resolve,
            || gate_poll(gate),
        )?;
        let (report, applied) = prepare_sgs(
            &self.grid,
            Some(self.forces()),
            &mut workspace.trial.collisions,
            sgs,
            || gate_poll(gate),
        )?;
        self.step_pooled(pool, gate, applied)?;
        Ok(publish(workspace, report))
    }

    /// Pooled source/SGS composition with one original force-history advance.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_material_rheology_and_sgs<
        'w,
        P: KernelRunner,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        reference: &MaterialScaling<'_, R>,
        forces_at_end: Vec<[f64; 3]>,
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: RheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologySgsStep3<'w, L::Collision, F>, RheologySgsError3<L::Error>> {
        prepare_workspace(
            &self.grid,
            reference,
            Some(self.forces()),
            workspace,
            resolve,
            || gate_poll(gate),
        )?;
        let (report, applied) = prepare_sgs(
            &self.grid,
            Some(self.forces()),
            &mut workspace.trial.collisions,
            sgs,
            || gate_poll(gate),
        )?;
        self.advance_pooled(pool, gate, applied, forces_at_end)?;
        Ok(publish(workspace, report))
    }
}

#[cfg(test)]
mod tests;
