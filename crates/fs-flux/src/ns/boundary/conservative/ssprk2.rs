//! Conservative SSPRK2 density stages feeding one native BDF1 momentum solve.
//!
//! The density/inventory subsystem is second order in time for its frozen
//! original volume advector. The coupled macro momentum update remains the
//! original FIRST-ORDER BDF1 map. Its single assembly/solve uses original old
//! inertia, full-interval final density and the one RK-weighted compatible
//! mass flux. No momentum or pressure solve occurs inside the density stages.

use super::{
    DensityMomentumAssembly, DensityMomentumControls, DensityMomentumStep, DensityTransportCell,
    FluxParams, FluxSystem, MaterialFluxNumerics, MaterialNavierStokesBoundary, Result,
    StokesAssembly, admit_transport_buffers, allocate, invalid, poll, positive,
};
use crate::ns::transport::density::ssprk2::{DensitySspRk2Controls, DensitySspRk2Report};

/// Bounded density RK work inside one original native BDF1 macrostep.
#[derive(Clone, Copy, Debug)]
pub struct DensityMomentumSspRk2Controls {
    /// Original physical normalization, depth, divergence and per-owner storage
    /// admission. `maximum_courant` applies to EACH density Euler stage.
    pub momentum: DensityMomentumControls,
    /// Maximum dyadic physical density subintervals. Each has two original
    /// Euler evaluations; their actual count is recorded in `density_stages`.
    pub maximum_substeps: usize,
}

/// Accepted RK density/inventory state and one original native BDF1 flow state.
pub struct DensityMomentumSspRk2Step<const N: usize> {
    /// Original velocity, physical pressure, final cells, macro accounts and
    /// separated numerical kinetic losses. The momentum time order stays one.
    pub momentum: DensityMomentumStep<N>,
    /// Actual density RK work and complete weighted physical inventory report.
    /// `transport` is identical to `momentum.transport`.
    pub density_stages: DensitySspRk2Report<N>,
}

/// RK transport is complete; one original native macro flow solve is pending.
pub struct DensityMomentumSspRk2Assembly<'a, 'm, const N: usize> {
    momentum: DensityMomentumAssembly<'a, 'm, N>,
    density_stages: DensitySspRk2Report<N>,
}

impl<'a, 'm, const N: usize> DensityMomentumSspRk2Assembly<'a, 'm, N> {
    /// Original native matrix/RHS, pressure and boundary/reaction rows.
    pub fn momentum(&self) -> &StokesAssembly<'a, 'm> {
        self.momentum.momentum()
    }

    /// Final positive density and ALL actually transported specific inventories.
    pub fn transported_cells(&self) -> &[DensityTransportCell<N>] {
        self.momentum.transported_cells()
    }

    /// Accepted physical interval, actual density stage work and RK budgets.
    pub const fn density_stages(&self) -> &DensitySspRk2Report<N> {
        &self.density_stages
    }

    /// Execute exactly one original native BDF1 momentum/pressure solve and
    /// its physical observations before returning a coupled candidate. The
    /// existing synchronous solve is polled before/after, not inside it.
    pub fn solve(&self, cancelled: impl FnMut() -> bool) -> Result<DensityMomentumSspRk2Step<N>> {
        Ok(DensityMomentumSspRk2Step {
            momentum: self.momentum.solve(cancelled)?,
            density_stages: self.density_stages,
        })
    }
}

impl<'m> FluxSystem<'m> {
    /// Prepare conservative SSPRK2 density/inventory stages and ONE original
    /// native BDF1 macro momentum/pressure assembly. Physical units, boundary
    /// policies and endpoint force/viscosity laws follow
    /// [`Self::assemble_density_momentum_bdf1`].
    ///
    /// The volume advector, geometry and exterior donors remain frozen. Stage
    /// states are convex-combined in extensive coordinates, and actual face
    /// amounts retain the same RK weights as the averaged donor density. One
    /// original compatible mass flux then supplies the native quadratic face
    /// moments and macro density divergence. No averaged specific state is
    /// invented, and all N passive inventories use the actual RK transfers.
    ///
    /// Initial and full-interval final densities, the macro timestep and that
    /// one averaged flux enter the existing shared native assembly unchanged.
    /// Force and viscosity see final endpoint cells once per original sampling
    /// site. Pressure, reactions and numerical kinetic-loss separation retain
    /// their original meanings. There are no density-stage momentum solves.
    ///
    /// Only the frozen-advector density/inventory ODE has second-order time
    /// integration here. Density uses this system's configured donor-cell or
    /// bounded linear reconstruction. Specific-inventory spatial mixing, the
    /// native BDF1 momentum map and its frozen-advector nonlinear coupling
    /// remain first order. No new EOS,
    /// diffusion, phase source, moving geometry or nonlinear adjoint is supplied.
    ///
    /// Original per-owner storage caps include the ENTIRE density stage
    /// workspace together. The native matrix/linear solver has separate storage
    /// ownership. Refusal or cancellation returns no partial coupled candidate
    /// and leaves caller state unchanged.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn assemble_density_momentum_ssprk2_bdf1<'a, const N: usize>(
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
        controls: DensityMomentumSspRk2Controls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensityMomentumSspRk2Assembly<'a, 'm, N>> {
        poll(&mut cancelled)?;
        let maximum_substeps = controls.maximum_substeps;
        let controls = controls.momentum;
        if !positive(controls.reference_density_kg_m3)
            || !positive(controls.depth_m)
            || maximum_substeps == 0
        {
            return Err(invalid(
                "positive physical reference density/depth and nonzero density RK work budget required",
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
        let mut work = plan.ssprk2_workspace::<N>(controls.maximum_transport_storage_bytes)?;
        let mut next = allocate(nc, DensityTransportCell::default())?;
        let density_stages = work.advance_into(
            cells,
            inflow,
            timestep_s,
            DensitySspRk2Controls {
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
        let momentum = self.assemble_transported_density_momentum_bdf1(
            params,
            dynamic_viscosity,
            force_density,
            boundaries,
            previous_velocity,
            cells,
            next,
            density_stages.transport,
            mass_flux,
            inv_dt,
            controls,
            &mut cancelled,
        )?;
        Ok(DensityMomentumSspRk2Assembly {
            momentum,
            density_stages,
        })
    }

    /// Prepare and solve [`Self::assemble_density_momentum_ssprk2_bdf1`]
    /// with one original native BDF1 macro momentum/pressure solve.
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_ssprk2_bdf1_step<const N: usize>(
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
        controls: DensityMomentumSspRk2Controls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensityMomentumSspRk2Step<N>> {
        self.assemble_density_momentum_ssprk2_bdf1(
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
