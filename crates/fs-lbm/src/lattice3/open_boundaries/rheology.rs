//! Original material rheology on native open reservoirs.

use super::*;
use crate::lattice3::ViscousLattice3;
use crate::lattice3::rheology::{
    RheologyError3, RheologyPoint3, RheologyStep3, RheologyWorkspace3,
    prepare_workspace_with_gradient,
};
use crate::material::MaterialScaling;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    // Preserve the closed-grid interior/periodic and halfway-wall stencil,
    // replacing ONLY the derivative normal to an actual open face. Otherwise
    // nonzero inlet velocity would be mistaken for wall shear. Pressure faces
    // inherit this owner's zero-normal-velocity-gradient donor condition,
    // except characteristic outlets whose normal velocity carries a wave.
    pub(super) fn rheology_gradient(&self, cell: usize, velocity: &[[f64; 3]]) -> [[f64; 3]; 3] {
        let mut gradient =
            self.grid
                .velocity_gradient_with_walls(cell, velocity, |cell, axis, upper| {
                    self.plan.gradient_wall_velocity(cell, axis, upper)
                });
        let Some((face, first)) = self.plan.locate(cell) else {
            return gradient;
        };
        let axis = face as usize / 2;
        if matches!(
            self.plan.boundaries.face(face),
            FaceBoundary3::Pressure { .. }
        ) && self.characteristic_outlets[face as usize].is_none()
        {
            for row in &mut gradient {
                row[axis] = 0.0;
            }
            return gradient;
        }
        let upper = face as usize % 2 == 1;
        let second = if upper {
            first - self.plan.stride
        } else {
            first + self.plan.stride
        };
        let sign = if upper { -1.0 } else { 1.0 };
        for component in 0..3 {
            let first_difference = velocity[first][component] - velocity[cell][component];
            // Quadratic one-sided slope at the on-site reservoir, exact for
            // affine/quadratic velocity. A solid second donor permits only the
            // linear one-sided estimate; never sample its synthetic zero.
            let derivative = if self.grid.solid[second] {
                first_difference
            } else {
                2.0 * first_difference
                    - 0.5 * (velocity[second][component] - velocity[cell][component])
            };
            gradient[component][axis] = sign * derivative;
        }
        gradient
    }
}

impl<L: EquilibriumLattice3<Q> + ViscousLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    /// Reserve the same two source epochs and velocity scratch used by native
    /// closed-grid material rheology. No material source is queried here.
    pub fn rheology_workspace<F: FluidTransport>(
        &self,
    ) -> Result<RheologyWorkspace3<L::Collision, F>, GridError3<L::Error>> {
        self.grid.rheology_workspace()
    }

    /// Resolve incoming physical strain and original material viscosity, then
    /// advance the original native open-reservoir transaction exactly once.
    ///
    /// Physical velocities include the current uniform half force. Interior,
    /// periodic and wall slopes retain the native grid stencil, with each
    /// exterior wall's actual prescribed velocity and stationary solid voxels.
    /// Velocity faces use an on-site quadratic one-sided slope, falling back
    /// to linear only when the second interior donor is solid. Pressure faces
    /// use the reservoir's zero-normal-velocity-gradient condition; attached
    /// characteristic outlets use the actual one-sided velocity slope. At seams,
    /// tangential derivatives still retain the native halfway-wall model.
    ///
    /// Sources are queried once per fluid cell at sqrt(2 dev(D):dev(D))/dt.
    /// The same frozen source mu/rho0 sets the actual shear relaxation; bulk,
    /// higher-order, TRT and equilibrium choices retain native semantics.
    /// Sources own temperature/composition and must not advance external
    /// history during tentative queries. No clipping, implicit constitutive
    /// solve, thermal evolution or SGS viscosity is added.
    ///
    /// Populations, reservoir exchange and source epoch publish together only
    /// after success. The workspace allocates no adapter-owned storage during
    /// steps; source-owned allocations and pooled scheduling are unchanged.
    pub fn step_serial_with_material_rheology_workspace<
        'w,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        origin: [f64; 3],
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        (
            &'w RheologyStep3<L::Collision, F>,
            OpenBoundaryMomentumTransport3,
        ),
        RheologyError3<L::Error>,
    > {
        prepare_workspace_with_gradient(
            &self.grid,
            reference,
            None,
            workspace,
            |cell, velocity| self.rheology_gradient(cell, velocity),
            resolve,
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )?;
        let exchange = self.step_serial(cx, Some(workspace.pending_collisions()), origin)?;
        Ok((workspace.publish(), exchange))
    }

    /// Pooled counterpart: original source queries stay serial and ordered;
    /// the existing native collision/stream kernels retain tile determinism.
    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled_with_material_rheology_workspace<
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
        origin: [f64; 3],
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        (
            &'w RheologyStep3<L::Collision, F>,
            OpenBoundaryMomentumTransport3,
        ),
        RheologyError3<L::Error>,
    > {
        prepare_workspace_with_gradient(
            &self.grid,
            reference,
            None,
            workspace,
            |cell, velocity| self.rheology_gradient(cell, velocity),
            resolve,
            || {
                if gate.is_requested() {
                    Err(GridError3::Cancelled)
                } else {
                    Ok(())
                }
            },
        )?;
        let exchange =
            self.step_pooled(pool, gate, Some(workspace.pending_collisions()), origin)?;
        Ok((workspace.publish(), exchange))
    }
}

#[cfg(test)]
mod tests;
