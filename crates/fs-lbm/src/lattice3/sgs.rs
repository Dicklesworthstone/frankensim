//! SGS integration through the existing per-cell collision owner.

use super::{Grid3, GridError3, Lattice3, coordinates, macro_cell};
use crate::kinematics::{FlowKinematics3, KinematicsError, KinematicsWorkspace3};
use crate::sgs3::{
    SgsConfig, SgsError, SgsLattice3, SgsReport, SgsStepError, SgsWorkspace3, eddy_viscosity,
};
use fs_exec::{CancelGate, Cx, KernelRunner};

mod channel;
mod dynamic;
mod forced;

impl<L: SgsLattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Allocate once for repeated SGS steps. Allocation failure is returned;
    /// no populations or molecular collision settings are changed.
    pub fn sgs_workspace(&self) -> Result<SgsWorkspace3<L::Collision>, GridError3<L::Error>> {
        let cells = self.f.len();
        let mut velocity = Vec::new();
        let mut density = Vec::new();
        let mut collisions = Vec::new();
        velocity
            .try_reserve_exact(cells)
            .map_err(|_| GridError3::Allocation)?;
        density
            .try_reserve_exact(cells)
            .map_err(|_| GridError3::Allocation)?;
        collisions
            .try_reserve_exact(cells)
            .map_err(|_| GridError3::Allocation)?;
        velocity.resize(cells, [0.0; 3]);
        density.resize(cells, 0.0);
        collisions.resize(cells, self.collision);
        Ok(SgsWorkspace3 {
            velocity,
            density,
            collisions,
        })
    }

    /// Prepare an explicit total-viscosity collision field from the current
    /// resolved velocity. `molecular` may supply the existing material-bound
    /// collision field; `None` uses the grid's constructor law. Never pass a
    /// previously augmented field as molecular input (that would compound SGS).
    ///
    /// Gradients are centered in the interior and periodic directions. At a
    /// stationary halfway wall a quadratic one-sided reconstruction uses the
    /// actual zero wall velocity, the cell centre, and its fluid neighbour.
    /// The normal derivative is exact for a quadratic profile at such walls;
    /// two opposed walls in a one-cell gap give the symmetric derivative zero.
    /// These are this grid's existing boundaries; open/moving walls are absent.
    ///
    /// The full velocity gradient is required by WALE. A cached macroscopic
    /// sweep avoids recomputing six neighbours' Q populations per cell. The
    /// scalar diagnostics and preparation use fixed cell order independent of
    /// the stepping worker count. Failed preparation leaves the grid unchanged;
    /// workspace contents must only be consumed after success.
    pub fn prepare_sgs(
        &self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        self.prepare_sgs_impl(config, molecular, workspace, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })
    }

    /// One opt-in exploratory SGS step, with failure-atomic publication and
    /// an instantaneous pre-step dissipation report. No model persists into a
    /// later ordinary step. The reusable workspace eliminates SGS allocations.
    pub fn step_serial_with_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        let report = self.prepare_sgs(cx, config, molecular, workspace)?;
        self.step_serial_with_collisions(cx, &workspace.collisions)?;
        Ok(report)
    }

    /// Pooled collision/streaming with exactly the same SGS arithmetic and
    /// bounded cancellation points as the serial owner. Preparation itself is
    /// a linear, allocation-free serial sweep; it does not spawn workers.
    pub fn step_pooled_with_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        let report = self.prepare_sgs_impl(config, molecular, workspace, || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        })?;
        self.step_pooled_with_collisions(pool, gate, &workspace.collisions)?;
        Ok(report)
    }

    fn prepare_sgs_impl(
        &self,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
        checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        self.prepare_sgs_with_forces(config, molecular, workspace, None, checkpoint)
    }

    // A spatial force slice always belongs to the SAME accepted epoch as f.
    // None preserves the original uniform-force path and arithmetic.
    fn prepare_sgs_with_forces(
        &self,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
        forces: Option<&[[f64; 3]]>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        checkpoint()?;
        config.validate()?;
        let cells = self.f.len();
        if workspace.velocity.len() != cells
            || workspace.density.len() != cells
            || workspace.collisions.len() != cells
            || molecular.is_some_and(|field| field.len() != cells)
            || forces.is_some_and(|field| field.len() != cells)
        {
            return Err(GridError3::Input("SGS workspace or molecular field length").into());
        }
        for cell in 0..cells {
            if cell % 64 == 0 {
                checkpoint()?;
            }
            if self.solid[cell] {
                workspace.velocity[cell] = [0.0; 3];
                workspace.density[cell] = 0.0;
            } else {
                let force = forces.map_or(self.force, |field| field[cell]);
                let (rho, u) =
                    macro_cell::<L, Q>(&self.f[cell], force, cell, self.budget.maximum_mach)?;
                workspace.velocity[cell] = u;
                workspace.density[cell] = rho;
            }
        }
        let mut report = SgsReport::default();
        for cell in 0..cells {
            if cell % 64 == 0 {
                checkpoint()?;
            }
            let original = molecular.map_or(self.collision, |field| field[cell]);
            workspace.collisions[cell] = original;
            let nu = L::molecular_viscosity(original)
                .map_err(|source| GridError3::Collision { cell, source })?;
            if self.solid[cell] {
                continue;
            }
            if !(nu.is_finite() && nu > 0.0) {
                return Err(SgsError::InvalidMolecularViscosity { cell }.into());
            }
            let gradient = self.velocity_gradient(cell, &workspace.velocity);
            let eddy = eddy_viscosity(config.model, config.filter_width_lu, gradient)?;
            let ratio = eddy / nu;
            if !ratio.is_finite() || ratio > config.maximum_eddy_to_molecular {
                return Err(SgsError::ViscosityRatioExceeded {
                    cell,
                    ratio,
                    maximum: config.maximum_eddy_to_molecular,
                }
                .into());
            }
            // Avoid a tau -> viscosity -> tau round trip for the zero control:
            // all collision families retain their original bits exactly.
            if eddy != 0.0 {
                let tau = 0.5 + (nu + eddy) / L::CS2;
                workspace.collisions[cell] = L::with_shear_relaxation(original, tau)
                    .map_err(|source| GridError3::Collision { cell, source })?;
            }
            let mut strain2 = 0.0;
            let dilation_third = (gradient[0][0] + gradient[1][1] + gradient[2][2]) / 3.0;
            for i in 0..3 {
                for j in 0..3 {
                    let s = 0.5 * (gradient[i][j] + gradient[j][i])
                        - if i == j { dilation_third } else { 0.0 };
                    strain2 += s * s;
                }
            }
            let factor = 2.0 * workspace.density[cell] * strain2;
            report.fluid_cells += 1;
            report.resolved_dissipation_lu += factor * nu;
            report.sgs_dissipation_lu += factor * eddy;
            report.maximum_eddy_viscosity_lu = report.maximum_eddy_viscosity_lu.max(eddy);
            report.maximum_eddy_to_molecular = report.maximum_eddy_to_molecular.max(ratio);
        }
        checkpoint()?;
        if !report.resolved_dissipation_lu.is_finite() || !report.sgs_dissipation_lu.is_finite() {
            return Err(SgsError::NonFiniteViscosity.into());
        }
        Ok(report)
    }
}

// One velocity-gradient stencil for molecular rheology and exploratory SGS.
// The caller supplies the same-time, force-centered macroscopic velocity.
impl<L: Lattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Resolved full gradients and vorticity of the accepted, force-centered
    /// state, without enabling an SGS closure. The returned lattice-unit
    /// field is x-fastest with `None` at solids. Reuses the rheology/SGS
    /// second-order stencil and its stationary halfway-wall geometry.
    pub fn kinematics<'a>(
        &self,
        gate: &CancelGate,
        workspace: &'a mut KinematicsWorkspace3,
    ) -> Result<&'a [Option<FlowKinematics3>], KinematicsError> {
        self.kinematics_with_forces(gate, workspace, None)
    }

    pub(super) fn kinematics_with_forces<'a>(
        &self,
        gate: &CancelGate,
        workspace: &'a mut KinematicsWorkspace3,
        forces: Option<&[[f64; 3]]>,
    ) -> Result<&'a [Option<FlowKinematics3>], KinematicsError> {
        workspace.admit(self.dims, gate)?;
        for cell in 0..self.f.len() {
            if cell % 64 == 0 {
                crate::kinematics::checkpoint(gate)?;
            }
            workspace.velocity[cell] = if self.solid[cell] {
                [0.0; 3]
            } else {
                let force = forces.map_or(self.force, |f| f[cell]);
                macro_cell::<L, Q>(&self.f[cell], force, cell, self.budget.maximum_mach)
                    .map_err(|_| KinematicsError::InvalidState {
                        cell: coordinates(cell, self.dims),
                    })?
                    .1
            };
        }
        for cell in 0..self.f.len() {
            if cell % 64 == 0 {
                crate::kinematics::checkpoint(gate)?;
            }
            workspace.field[cell] = if self.solid[cell] {
                None
            } else {
                Some(FlowKinematics3::from_gradient(
                    self.velocity_gradient(cell, &workspace.velocity),
                )?)
            };
        }
        crate::kinematics::checkpoint(gate)?;
        Ok(&workspace.field)
    }

    fn velocity_neighbour(
        &self,
        mut p: [usize; 3],
        axis: usize,
        positive: bool,
        velocity: &[[f64; 3]],
    ) -> Option<[f64; 3]> {
        if positive {
            if p[axis] + 1 == self.dims[axis] {
                if !self.periodic[axis] {
                    return None;
                }
                p[axis] = 0;
            } else {
                p[axis] += 1;
            }
        } else if p[axis] == 0 {
            if !self.periodic[axis] {
                return None;
            }
            p[axis] = self.dims[axis] - 1;
        } else {
            p[axis] -= 1;
        }
        let cell = (p[2] * self.dims[1] + p[1]) * self.dims[0] + p[0];
        if self.solid[cell] {
            None
        } else {
            Some(velocity[cell])
        }
    }

    pub(super) fn velocity_gradient(&self, cell: usize, velocity: &[[f64; 3]]) -> [[f64; 3]; 3] {
        self.velocity_gradient_with_walls(cell, velocity, |_, _, _| [0.0; 3])
    }

    // The caller must bind each missing neighbor to the SAME halfway-wall
    // velocity used by native streaming. Wall-relative differences preserve
    // uniform co-motion and the quadratic reconstruction at the fluid centre.
    pub(super) fn velocity_gradient_with_walls(
        &self,
        cell: usize,
        velocity: &[[f64; 3]],
        mut wall_velocity: impl FnMut(usize, usize, bool) -> [f64; 3],
    ) -> [[f64; 3]; 3] {
        let p = coordinates(cell, self.dims);
        let mut gradient = [[0.0; 3]; 3];
        for axis in 0..3 {
            let minus = self.velocity_neighbour(p, axis, false, velocity);
            let plus = self.velocity_neighbour(p, axis, true, velocity);
            let low_wall = if minus.is_none() {
                wall_velocity(cell, axis, false)
            } else {
                [0.0; 3]
            };
            let high_wall = if plus.is_none() {
                wall_velocity(cell, axis, true)
            } else {
                [0.0; 3]
            };
            for component in 0..3 {
                gradient[component][axis] = match (minus, plus) {
                    (Some(low), Some(high)) => 0.5 * (high[component] - low[component]),
                    // Keep the stationary arithmetic, including signed zero.
                    (None, Some(high)) if low_wall[component] == 0.0 => {
                        velocity[cell][component] + high[component] / 3.0
                    }
                    (Some(low), None) if high_wall[component] == 0.0 => {
                        -velocity[cell][component] - low[component] / 3.0
                    }
                    (None, Some(high)) => {
                        (velocity[cell][component] - low_wall[component])
                            + (high[component] - low_wall[component]) / 3.0
                    }
                    (Some(low), None) => {
                        -(velocity[cell][component] - high_wall[component])
                            - (low[component] - high_wall[component]) / 3.0
                    }
                    (None, None) if low_wall[component] == high_wall[component] => 0.0,
                    (None, None) => high_wall[component] - low_wall[component],
                };
            }
        }
        gradient
    }
}
