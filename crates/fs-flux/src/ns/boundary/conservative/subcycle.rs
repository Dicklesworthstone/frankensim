//! One native momentum/pressure solve after bounded conservative density work.
//!
//! The volume advector and exterior density/specific states remain frozen.
//! Density substeps supply their actual final inventories and ONE compatible
//! time-averaged mass flux to the original native conservative momentum owner.
//! Its old/new inertia, force, viscosity, pressure and solve use the caller's
//! macro interval. A larger admitted density interval does not imply improved
//! frozen-advector time accuracy or a fully implicit nonlinear coupling.

use super::{
    DensityMomentumAssembly, DensityMomentumControls, DensityMomentumStep, DensityTransportCell,
    FluxParams, FluxSystem, MaterialFluxNumerics, MaterialNavierStokesBoundary,
    PreparedDensityMomentum, Result, StokesAssembly, admit_transport_buffers, allocate, invalid,
    poll, positive,
};
use crate::ns::transport::density::subcycle::{DensitySubcycleControls, DensitySubcycleReport};

/// Explicit admission for density work within one original native macro solve.
#[derive(Clone, Copy, Debug)]
pub struct DensityMomentumSubcycleControls {
    /// Original physical normalization, geometry and per-owner storage caps.
    /// `maximum_courant` applies to EACH density substep; the macro interval
    /// can exceed that outgoing-volume Courant allowance.
    pub momentum: DensityMomentumControls,
    /// Maximum number of equal dyadic density steps before one momentum solve.
    pub maximum_substeps: usize,
}

/// One accepted coupled macrostep with its explicit density work/accounting.
pub struct DensityMomentumSubcycleStep<const N: usize> {
    /// Original native velocity, physical pressure, final density/specifics,
    /// macro inventory report and separately observed numerical kinetic losses.
    pub momentum: DensityMomentumStep<N>,
    /// Complete macro transport report plus actual density substep count/CFL.
    /// `transport` is identical to `momentum.transport`.
    pub subcycles: DensitySubcycleReport<N>,
}

/// Density work is complete; the one original macro momentum solve is pending.
pub struct DensityMomentumSubcycleAssembly<'a, 'm, const N: usize> {
    momentum: DensityMomentumAssembly<'a, 'm, N>,
    subcycles: DensitySubcycleReport<N>,
}

impl<'a, 'm, const N: usize> DensityMomentumSubcycleAssembly<'a, 'm, N> {
    /// The original native matrix/RHS, including pressure and boundary rows.
    pub fn momentum(&self) -> &StokesAssembly<'a, 'm> {
        self.momentum.momentum()
    }

    /// Positive final density and all actually transported specific inventories.
    pub fn transported_cells(&self) -> &[DensityTransportCell<N>] {
        self.momentum.transported_cells()
    }

    /// Complete macro physical inventory/boundary totals and bounded work.
    pub const fn subcycles(&self) -> &DensitySubcycleReport<N> {
        &self.subcycles
    }

    /// Perform exactly one original native momentum/pressure solve, then its
    /// complete physical observations before returning a coupled candidate.
    /// Cancellation semantics are the original synchronous solver's: checked
    /// before/after the solve and throughout observations, not inside it.
    pub fn solve(&self, cancelled: impl FnMut() -> bool) -> Result<DensityMomentumSubcycleStep<N>> {
        Ok(DensityMomentumSubcycleStep {
            momentum: self.momentum.solve(cancelled)?,
            subcycles: self.subcycles,
        })
    }
}

impl<'m> FluxSystem<'m> {
    /// Prepare one mass-consistent native BDM1/P0 macrostep with bounded density
    /// subcycling. Physical units, native boundary policies and endpoint force/
    /// viscosity laws match [`Self::assemble_density_momentum_bdf1`].
    ///
    /// The fixed original volume advector and exterior donors are used at EVERY
    /// density substep. Every passive inventory is actually transported; only
    /// PRE-step donor densities are averaged to reconstruct the single flux.
    /// Its macro divergence matches total density change, and its original
    /// quadratic face moments retain the native kinetic balance. Roundoff
    /// defects remain reported; density clipping or mass repair is not added.
    ///
    /// The original old/full-macro-new mass matrices and macro timestep enter
    /// one momentum assembly. Endpoint force/viscosity callbacks see only the
    /// final transported cells. The density loop performs no global solves.
    /// First-order frozen-advector accuracy still limits useful macro duration.
    ///
    /// Storage is bounded per original owner, with the ENTIRE reusable density
    /// subcycle workspace admitted together. Work and representability refusal
    /// or cancellation never publish a partially coupled state. As with the
    /// original API, the native linear-solver storage is separately owned.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn assemble_density_momentum_subcycled_bdf1<'a, const N: usize>(
        &'a self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        advector: &[f64],
        previous_velocity: &'a [f64],
        cells: &'a [DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumSubcycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensityMomentumSubcycleAssembly<'a, 'm, N>> {
        let (mut prepared, subcycles) = self.prepare_density_momentum_subcycled_bdf1(
            numerics,
            force_density,
            boundaries,
            advector,
            previous_velocity,
            cells,
            inflow,
            timestep_s,
            controls,
            &mut cancelled,
        )?;
        let weights = core::mem::take(&mut prepared.current_weights);
        let momentum = prepared.assemble(
            &|cell, point| dynamic_viscosity(cell, point, &prepared.cells[cell]),
            weights,
        )?;
        poll(&mut cancelled)?;
        Ok(DensityMomentumSubcycleAssembly {
            momentum: prepared.into_assembly(momentum),
            subcycles,
        })
    }

    // Both frozen-advector and endpoint Picard callers retain the original
    // substeps' accumulated mass flux, never reconstruct it from final density.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare_density_momentum_subcycled_bdf1<'a, const N: usize>(
        &'a self,
        numerics: MaterialFluxNumerics,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        advector: &[f64],
        previous_velocity: &'a [f64],
        cells: &'a [DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumSubcycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(PreparedDensityMomentum<'a, 'm, N>, DensitySubcycleReport<N>)> {
        poll(&mut cancelled)?;
        let maximum_substeps = controls.maximum_substeps;
        let controls = controls.momentum;
        if !positive(controls.reference_density_kg_m3)
            || !positive(controls.depth_m)
            || maximum_substeps == 0
        {
            return Err(invalid(
                "positive physical reference density/depth and nonzero density work budget required",
            ));
        }
        let inv_dt = self.admit_mixed_bdf1_history(previous_velocity, timestep_s)?;
        let params = FluxParams {
            nu: 1.0,
            sigma: numerics.penalty,
            tol: numerics.linear_tolerance,
            max_iters: numerics.max_linear_iterations,
        };
        let nc = self.mesh.tris.len();
        admit_transport_buffers::<N>(
            nc,
            self.mesh.edges.len(),
            controls.maximum_transport_storage_bytes,
        )?;
        let plan = self.prepare_density_transport(
            advector,
            controls.depth_m,
            controls.relative_volume_divergence_tolerance,
            controls.maximum_transport_storage_bytes,
            &mut cancelled,
        )?;
        let mut work = plan.subcycle_workspace::<N>(controls.maximum_transport_storage_bytes)?;
        let mut next = allocate(nc, DensityTransportCell::default())?;
        let subcycles = work.advance_into(
            cells,
            inflow,
            timestep_s,
            DensitySubcycleControls {
                maximum_courant: controls.maximum_courant,
                maximum_substeps,
            },
            &mut next,
            &mut cancelled,
        )?;
        let mass_flux =
            work.prepare_mass_flux(controls.maximum_transport_storage_bytes, &mut cancelled)?;
        drop(work);
        drop(plan);
        let prepared = self.prepare_transported_density_momentum_bdf1(
            params,
            force_density,
            boundaries,
            previous_velocity,
            cells,
            next,
            subcycles.transport,
            mass_flux,
            inv_dt,
            controls,
            &mut cancelled,
        )?;
        Ok((prepared, subcycles))
    }

    /// Prepare and solve [`Self::assemble_density_momentum_subcycled_bdf1`]
    /// with one native macro momentum/pressure solve and no caller mutation.
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_subcycled_bdf1_step<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        advector: &[f64],
        previous_velocity: &[f64],
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumSubcycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensityMomentumSubcycleStep<N>> {
        self.assemble_density_momentum_subcycled_bdf1(
            numerics,
            dynamic_viscosity,
            force_density,
            boundaries,
            advector,
            previous_velocity,
            cells,
            inflow,
            timestep_s,
            controls,
            &mut cancelled,
        )?
        .solve(&mut cancelled)
    }
}
