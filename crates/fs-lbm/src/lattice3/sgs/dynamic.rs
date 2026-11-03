//! Dynamic SGS uses the original lattice collision and streaming transaction.

use super::*;
use crate::sgs3::dynamic::{DynamicSgsConfig, DynamicSgsError, PeriodicDynamicSgs};
use crate::sgs3::{DynamicSgsStepReport, DynamicSgsWorkspace3};

impl<L: SgsLattice3<Q>, const Q: usize> Grid3<L, Q> {
    fn admit_dynamic_topology(
        &self,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(), SgsStepError<L::Error>> {
        if self.periodic != [true; 3] {
            return Err(DynamicSgsError::Invalid(
                "dynamic test filtering requires a fully periodic domain without solids",
            )
            .into());
        }
        for tile in self.solid.chunks(64) {
            checkpoint()?;
            if tile.iter().any(|&solid| solid) {
                return Err(DynamicSgsError::Invalid(
                    "dynamic test filtering requires a fully periodic domain without solids",
                )
                .into());
            }
        }
        Ok(())
    }

    /// Allocate the complete macroscopic, collision and test-filter payload.
    /// The byte bound covers both workspaces, excluding allocator metadata.
    /// This does not change populations, collision rates or the accepted clock.
    pub fn dynamic_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<DynamicSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        self.admit_dynamic_topology(|| Ok(()))?;
        let cells = self.f.len();
        let local = cells
            .checked_mul(4 * core::mem::size_of::<f64>() + core::mem::size_of::<L::Collision>())
            .ok_or(DynamicSgsError::Invalid("workspace size overflow"))?;
        let required = local
            .checked_add(PeriodicDynamicSgs::required_bytes(cells)?)
            .ok_or(DynamicSgsError::Invalid("workspace size overflow"))?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            }
            .into());
        }
        Ok(DynamicSgsWorkspace3 {
            local: self.sgs_workspace()?,
            filter: PeriodicDynamicSgs::new(self.dims, 1.0, maximum_bytes - local)?,
        })
    }

    /// Prepare a density-weighted dynamic viscosity from same-time physical
    /// velocities. The base filter is one lattice cell, the test filter three.
    /// This retains each collision family's existing rate coupling: BGK has
    /// one rate and D3Q19 central moments couple all second-order rates;
    /// D3Q27 central/cumulant collisions retain independent trace and higher
    /// rates. No separate SGS kinetic-energy equation is added.
    pub fn prepare_dynamic_sgs(
        &self,
        cx: &Cx<'_>,
        config: DynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut DynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        self.prepare_dynamic_sgs_impl(config, molecular, workspace, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })
    }

    /// Step with an explicit dynamic closure. Failure publishes no populations
    /// or time, and the constructor's molecular collision is never overwritten.
    pub fn step_serial_with_dynamic_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: DynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut DynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.prepare_dynamic_sgs(cx, config, molecular, workspace)?;
        self.step_serial_with_collisions(cx, workspace.collisions())?;
        Ok(report)
    }

    /// Pooled stepping reuses the identical deterministic preparation and
    /// original transactional collision/streaming owner.
    pub fn step_pooled_with_dynamic_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: DynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut DynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.prepare_dynamic_sgs_impl(config, molecular, workspace, || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        })?;
        self.step_pooled_with_collisions(pool, gate, workspace.collisions())?;
        Ok(report)
    }

    fn prepare_dynamic_sgs_impl(
        &self,
        config: DynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut DynamicSgsWorkspace3<L::Collision>,
        checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        self.prepare_dynamic_sgs_with_forces(config, molecular, workspace, None, checkpoint)
    }

    pub(super) fn prepare_dynamic_sgs_with_forces(
        &self,
        config: DynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut DynamicSgsWorkspace3<L::Collision>,
        forces: Option<&[[f64; 3]]>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        checkpoint()?;
        self.admit_dynamic_topology(&mut checkpoint)?;
        config.validate()?;
        let cells = self.f.len();
        let local = &mut workspace.local;
        if local.velocity.len() != cells
            || local.density.len() != cells
            || local.collisions.len() != cells
            || workspace.filter.shape() != self.dims
            || molecular.is_some_and(|field| field.len() != cells)
            || forces.is_some_and(|field| field.len() != cells)
        {
            return Err(
                GridError3::Input("dynamic SGS workspace or molecular field mismatch").into(),
            );
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
        let fit = workspace
            .filter
            .prepare(config, &local.density, &local.velocity, || {
                checkpoint().is_err()
            })?;
        let dissipation = self.apply_dynamic_sgs_viscosities(
            config.maximum_eddy_to_molecular,
            molecular,
            workspace.filter.viscosities(),
            local,
            checkpoint,
        )?;
        Ok(DynamicSgsStepReport { fit, dissipation })
    }

    // Both nonlocal geometries retain the same native molecular law, shear
    // conversion and diagnostics. Prepared macro/viscosity lengths are admitted
    // by their caller; no filter or boundary is inferred in this conversion.
    pub(super) fn apply_dynamic_sgs_viscosities(
        &self,
        maximum_ratio: f64,
        molecular: Option<&[L::Collision]>,
        viscosities: &[f64],
        local: &mut SgsWorkspace3<L::Collision>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        let mut dissipation = SgsReport::default();
        for (cell, &eddy) in viscosities.iter().enumerate() {
            if cell % 64 == 0 {
                checkpoint()?;
            }
            let original = molecular.map_or(self.collision, |field| field[cell]);
            let nu = L::molecular_viscosity(original)
                .map_err(|source| GridError3::Collision { cell, source })?;
            if !(nu.is_finite() && nu > 0.0) {
                return Err(SgsError::InvalidMolecularViscosity { cell }.into());
            }
            let ratio = eddy / nu;
            if !ratio.is_finite() || ratio > maximum_ratio {
                return Err(SgsError::ViscosityRatioExceeded {
                    cell,
                    ratio,
                    maximum: maximum_ratio,
                }
                .into());
            }
            local.collisions[cell] = if eddy == 0.0 {
                original
            } else {
                L::with_shear_relaxation(original, 0.5 + (nu + eddy) / L::CS2)
                    .map_err(|source| GridError3::Collision { cell, source })?
            };
            let gradient = self.velocity_gradient(cell, &local.velocity);
            let third = (gradient[0][0] + gradient[1][1] + gradient[2][2]) / 3.0;
            let mut strain2 = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    let s =
                        0.5 * (gradient[i][j] + gradient[j][i]) - if i == j { third } else { 0.0 };
                    strain2 += s * s;
                }
            }
            let factor = 2.0 * local.density[cell] * strain2;
            dissipation.fluid_cells += 1;
            dissipation.resolved_dissipation_lu += factor * nu;
            dissipation.sgs_dissipation_lu += factor * eddy;
            dissipation.maximum_eddy_viscosity_lu = dissipation.maximum_eddy_viscosity_lu.max(eddy);
            dissipation.maximum_eddy_to_molecular =
                dissipation.maximum_eddy_to_molecular.max(ratio);
        }
        checkpoint()?;
        if !dissipation.resolved_dissipation_lu.is_finite()
            || !dissipation.sgs_dissipation_lu.is_finite()
        {
            return Err(SgsError::NonFiniteViscosity.into());
        }
        Ok(dissipation)
    }
}

#[cfg(test)]
mod tests;
