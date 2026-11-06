//! Lagrangian SGS memory on the native spatial/time-dependent force owner.
//!
//! Fits and pathline departures use the incoming Guo-centered velocity. The
//! original force transaction alone admits and publishes the endpoint force;
//! only after that succeeds may the matching turbulence history be committed.
//! A refused step leaves populations, force epoch, clock and history unchanged.

use super::{ForcedGrid3, ForcedLattice3, gate_checkpoint};
use crate::lattice3::GridError3;
use crate::sgs3::dynamic::lagrangian::LagrangianSgsConfig;
use crate::sgs3::{
    DynamicSgsStepReport, LagrangianSgsWorkspace3, SgsLattice3, SgsStepError,
};
use fs_exec::{CancelGate, Cx, KernelRunner};

impl<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Reuse the native periodic filter and history banks without detaching the
    /// spatial force owner. The byte gate and all-fluid topology checks remain
    /// those of the original Lagrangian workspace.
    pub fn lagrangian_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<LagrangianSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        self.grid.lagrangian_sgs_workspace(maximum_bytes)
    }

    /// Hold the current spatial force for one step with Lagrangian SGS.
    /// Fit and departure velocities include the current per-cell half-force.
    /// The constructor's molecular collision law, or the explicitly supplied
    /// molecular field, remains the source of every total shear rate.
    pub fn step_serial_with_lagrangian_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut LagrangianSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let (trial, dissipation) = self.grid.prepare_lagrangian_sgs_with_forces(
            config,
            molecular,
            &mut workspace.local,
            &mut workspace.filter,
            Some(self.forces()),
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )?;
        self.step_serial(cx, workspace.local.collisions())?;
        let fit = trial.report();
        trial.commit();
        Ok(DynamicSgsStepReport { fit, dissipation })
    }

    /// Fit at the current force epoch, then advance with the supplied NEXT
    /// force. The original Guo transaction supplies its trapezoidal source
    /// impulse; the endpoint force never enters the incoming SGS fit or
    /// pathline departure. History publishes only after the native step.
    pub fn advance_serial_with_lagrangian_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut LagrangianSgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let (trial, dissipation) = self.grid.prepare_lagrangian_sgs_with_forces(
            config,
            molecular,
            &mut workspace.local,
            &mut workspace.filter,
            Some(self.forces()),
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )?;
        self.advance_serial(cx, workspace.local.collisions(), forces_at_end)?;
        let fit = trial.report();
        trial.commit();
        Ok(DynamicSgsStepReport { fit, dissipation })
    }

    /// Deterministic fit followed by the native pooled held-force step.
    /// This reuses the original collision/filter arrays and adds no step
    /// preparation allocation. A native refusal discards the trial history.
    pub fn step_pooled_with_lagrangian_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut LagrangianSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let (trial, dissipation) = self.grid.prepare_lagrangian_sgs_with_forces(
            config,
            molecular,
            &mut workspace.local,
            &mut workspace.filter,
            Some(self.forces()),
            || gate_checkpoint(gate),
        )?;
        self.step_pooled(pool, gate, workspace.local.collisions())?;
        let fit = trial.report();
        trial.commit();
        Ok(DynamicSgsStepReport { fit, dissipation })
    }

    /// Pooled next-force advance with the same incoming fit and atomic
    /// population/force/clock/history publication as the serial entry.
    /// Periodic, all-fluid topology and exploratory-use restrictions remain.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_lagrangian_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut LagrangianSgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let (trial, dissipation) = self.grid.prepare_lagrangian_sgs_with_forces(
            config,
            molecular,
            &mut workspace.local,
            &mut workspace.filter,
            Some(self.forces()),
            || gate_checkpoint(gate),
        )?;
        self.advance_pooled(pool, gate, workspace.local.collisions(), forces_at_end)?;
        let fit = trial.report();
        trial.commit();
        Ok(DynamicSgsStepReport { fit, dissipation })
    }
}
