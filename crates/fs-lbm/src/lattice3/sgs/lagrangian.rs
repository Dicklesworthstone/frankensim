//! Accepted Lagrangian turbulence memory on the original lattice transaction.

use super::*;
use crate::sgs3::dynamic::DynamicSgsError;
use crate::sgs3::dynamic::lagrangian::{
    LagrangianDynamicSgs, LagrangianSgsConfig, LagrangianSgsTrial,
};
use crate::sgs3::{DynamicSgsStepReport, LagrangianSgsWorkspace3};

impl<L: SgsLattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Allocate a trajectory's local collision fields and two turbulence-memory
    /// banks. All three axes must be periodic, with no solids. The declared
    /// byte bound covers actual numerical capacities; stepping adds no arrays.
    pub fn lagrangian_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<LagrangianSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        self.admit_dynamic_topology(|| Ok(()))?;
        let cells = self.f.len();
        let local_bytes = cells
            .checked_mul(4 * core::mem::size_of::<f64>() + core::mem::size_of::<L::Collision>())
            .ok_or(DynamicSgsError::Invalid(
                "Lagrangian native workspace size overflow",
            ))?;
        let required = local_bytes
            .checked_add(LagrangianDynamicSgs::required_bytes(cells)?)
            .ok_or(DynamicSgsError::Invalid(
                "Lagrangian native workspace size overflow",
            ))?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            }
            .into());
        }
        let result = LagrangianSgsWorkspace3 {
            local: self.sgs_workspace()?,
            filter: LagrangianDynamicSgs::new(self.dims, [1.0; 3], maximum_bytes - local_bytes)?,
        };
        let required = result.workspace_bytes();
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            }
            .into());
        }
        Ok(result)
    }

    /// Advance the native D3Q19/D3Q27 flow and its Lagrangian dynamic SGS history
    /// in one transaction. Memory time is in lattice timesteps. Fits use the
    /// accepted Guo-centered physical velocity; the constructor's molecular
    /// collision or the explicitly supplied molecular field remains unchanged.
    /// D3Q27 retains its independent trace and higher collision rates.
    ///
    /// A closure, collision, Mach or cancellation refusal publishes neither
    /// populations, step count nor history. The fit clock is the incoming
    /// `completed_steps`, so failed trials can be retried with identical memory.
    /// This is exploratory periodic LES, without walls or a DNS claim.
    pub fn step_serial_with_lagrangian_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut LagrangianSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let (trial, dissipation) = self.prepare_lagrangian_sgs_with_forces(
            config,
            molecular,
            &mut workspace.local,
            &mut workspace.filter,
            None,
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )?;
        self.step_serial_with_collisions(cx, workspace.local.collisions())?;
        let fit = trial.report();
        trial.commit();
        Ok(DynamicSgsStepReport { fit, dissipation })
    }

    /// Identical deterministic fit followed by the original pooled step.
    /// Accepted history is swapped only after that native owner publishes.
    pub fn step_pooled_with_lagrangian_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut LagrangianSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let (trial, dissipation) = self.prepare_lagrangian_sgs_with_forces(
            config,
            molecular,
            &mut workspace.local,
            &mut workspace.filter,
            None,
            || {
                if gate.is_requested() {
                    Err(GridError3::Cancelled)
                } else {
                    Ok(())
                }
            },
        )?;
        self.step_pooled_with_collisions(pool, gate, workspace.local.collisions())?;
        let fit = trial.report();
        trial.commit();
        Ok(DynamicSgsStepReport { fit, dissipation })
    }

    // Keep the fit, departure velocity and dissipation in the same accepted
    // force frame. The endpoint force belongs to the native step alone.
    pub(super) fn prepare_lagrangian_sgs_with_forces<'a>(
        &self,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        local: &mut SgsWorkspace3<L::Collision>,
        filter: &'a mut LagrangianDynamicSgs,
        forces: Option<&[[f64; 3]]>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(LagrangianSgsTrial<'a>, SgsReport), SgsStepError<L::Error>> {
        checkpoint()?;
        self.admit_dynamic_topology(&mut checkpoint)?;
        config.validate()?;
        let cells = self.f.len();
        if local.velocity.len() != cells
            || local.density.len() != cells
            || local.collisions.len() != cells
            || filter.shape() != self.dims
            || molecular.is_some_and(|field| field.len() != cells)
            || forces.is_some_and(|field| field.len() != cells)
        {
            return Err(
                GridError3::Input("Lagrangian SGS workspace, molecular or force field mismatch").into(),
            );
        }
        // The public clock is integral. A rounded floating clock must never
        // make two distinct accepted steps appear to be the same fit time.
        if self.completed_steps() as u128 > (1u128 << 53) {
            return Err(DynamicSgsError::Invalid(
                "Lagrangian lattice clock exceeds exact f64 range",
            )
            .into());
        }
        for cell in 0..cells {
            if cell % 64 == 0 {
                checkpoint()?;
            }
            let force = forces.map_or(self.force, |field| field[cell]);
            let (density, velocity) =
                macro_cell::<L, Q>(&self.f[cell], force, cell, self.budget.maximum_mach)?;
            local.density[cell] = density;
            local.velocity[cell] = velocity;
        }
        let trial = filter.prepare(
            config,
            self.completed_steps() as f64,
            &local.density,
            &local.velocity,
            || checkpoint().is_err(),
        )?;
        let dissipation = self.apply_dynamic_sgs_viscosities(
            config.maximum_eddy_to_molecular,
            molecular,
            trial.viscosities(),
            local,
            checkpoint,
        )?;
        Ok((trial, dissipation))
    }
}
