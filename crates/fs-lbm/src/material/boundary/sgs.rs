//! Original material rheology and SGS on the same wall/open-boundary flow.
//!
//! The incoming resolved gradient supplies both molecular source shear and
//! eddy viscosity. Evaluating a nonlinear molecular law on the resolved strain
//! is an explicit LES approximation: it does not close the correlation between
//! unresolved strain and nonlinear molecular stress. No constitutive-memory,
//! thermal integration, immersed-body filtering or implicit source solve is
//! supplied here. Wall loads come from the actual total-viscosity collision.

use crate::d3q19::{
    BoundaryGrid3, BoundaryLink3, BoundarySgsError3, BoundarySgsReport3, BoundaryStepError3,
    CollisionModel3,
};
use crate::lattice3::{D3q19, ViscousLattice3};
use crate::material::immersed::{MaterialImmersedPoint3, MaterialImmersedSample3, strain_rate};
use crate::material::{MaterialLbmError, MaterialScaling, plan_momentum_scaling};
use crate::sgs3::{SgsConfig, SgsWorkspace3};
use fs_exec::{CancelGate, KernelRunner};
use fs_material::{MaterialError, fluid_transport::FluidTransport};

/// Original source errors remain distinct from closure and native-flow errors.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialBoundarySgsError3 {
    /// Native SGS, collision, streaming, wall-load or cancellation refusal.
    Sgs(BoundarySgsError3),
    /// Sampling budget or workspace binding was not admitted.
    Invalid(&'static str),
    /// Original material response, physical shear or momentum conversion failed.
    Source {
        cell: [usize; 3],
        source: MaterialLbmError,
    },
}

impl core::fmt::Display for MaterialBoundarySgsError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-BOUNDARY-MATERIAL-SGS3: {self:?}")
    }
}
impl std::error::Error for MaterialBoundarySgsError3 {}
impl From<BoundarySgsError3> for MaterialBoundarySgsError3 {
    fn from(error: BoundarySgsError3) -> Self {
        Self::Sgs(error)
    }
}

/// Reusable native SGS scratch and two preallocated original-source epochs.
/// The applied epoch and report change only after the native population swap.
/// Reconstruct after changing dimensions, physical scaling or collision law.
pub struct BoundaryRheologySgsWorkspace3<M> {
    flow: SgsWorkspace3<CollisionModel3>,
    dimensions: [usize; 3],
    scaling: [u64; 3],
    model: CollisionModel3,
    applied: Vec<MaterialImmersedSample3<M>>,
    pending: Vec<MaterialImmersedSample3<M>>,
    report: Option<BoundarySgsReport3>,
}

impl<M> BoundaryRheologySgsWorkspace3<M> {
    /// Original responses and incoming physical shear from the last successful
    /// step. Order is native x-fastest tiles, then x-fastest lanes, skipping
    /// permanent solids. `shear_relaxation_time` is MOLECULAR tau before SGS.
    /// Empty until the first accepted step; refusals retain this source epoch.
    pub fn applied(&self) -> &[MaterialImmersedSample3<M>] {
        &self.applied
    }

    /// SGS diagnostics and actual wall loads from the same accepted epoch.
    pub fn last_report(&self) -> Option<BoundarySgsReport3> {
        self.report
    }

    /// Retained numerical/source capacities, excluding allocator overhead and
    /// any additional heap storage owned by the original material responses.
    pub fn workspace_bytes(&self) -> usize {
        self.flow.velocity.capacity() * core::mem::size_of::<[f64; 3]>()
            + self.flow.density.capacity() * core::mem::size_of::<f64>()
            + self.flow.collisions.capacity() * core::mem::size_of::<CollisionModel3>()
            + (self.applied.capacity() + self.pending.capacity())
                * core::mem::size_of::<MaterialImmersedSample3<M>>()
    }

    fn accept(&mut self, report: BoundarySgsReport3) {
        core::mem::swap(&mut self.applied, &mut self.pending);
        self.report = Some(report);
    }
}

impl<R: FluidTransport + ?Sized> MaterialScaling<'_, R> {
    /// Allocate full-grid flow scratch and both source epochs once. The explicit
    /// cell budget includes solid padding; sources are sampled only in fluid.
    /// Borrowed concrete sources and `&dyn FluidTransport` remain supported.
    pub fn boundary_rheology_sgs_workspace<M>(
        &self,
        grid: &BoundaryGrid3,
        max_cells: usize,
        gate: &CancelGate,
    ) -> Result<BoundaryRheologySgsWorkspace3<M>, MaterialBoundarySgsError3> {
        checkpoint(gate)?;
        let dimensions = grid.dimensions();
        let count = dimensions
            .into_iter()
            .try_fold(1usize, usize::checked_mul)
            .filter(|&n| n <= max_cells)
            .ok_or(MaterialBoundarySgsError3::Invalid(
                "boundary rheology cell budget exceeded",
            ))?;
        let flow = grid.sgs_workspace()?;
        let mut applied = Vec::new();
        let mut pending = Vec::new();
        applied
            .try_reserve_exact(count)
            .map_err(|_| BoundarySgsError3::Allocation)?;
        pending
            .try_reserve_exact(count)
            .map_err(|_| BoundarySgsError3::Allocation)?;
        checkpoint(gate)?;
        Ok(BoundaryRheologySgsWorkspace3 {
            flow,
            dimensions,
            scaling: self.boundary_rheology_scaling(),
            model: grid.collision_model(),
            applied,
            pending,
            report: None,
        })
    }

    /// Advance the original wall/open-boundary owner with a fresh material
    /// response at each fluid cell's resolved shear and the selected SGS law.
    ///
    /// One force-centered gradient supplies BOTH `sqrt(2 D:D)/dt` [1/s] and
    /// eddy viscosity. The physical velocity and distance scales cancel dx.
    /// Native moving halfway-wall and on-site open-face derivatives are used;
    /// no second gradient, pressure boundary or wall model is introduced.
    /// Each source is queried once and its transport snapshot read once.
    /// Existing scaling converts that same mu/rho0 to molecular tau. The
    /// viscosity-ratio refusal and molecular dissipation use this LOCAL nu.
    /// Zero SGS retains the heterogeneous molecular field.
    ///
    /// This is incoming-time resolved-shear generalized-Newtonian rheology,
    /// with the nonlinear filtering approximation described by this module.
    /// Temperature/composition stay caller-owned and reference density must
    /// remain identical. Cp/k impose no thermal lattice bound or added heat.
    /// Retain the scaling's scenario speed bound. The original collision
    /// family, TRT magic product and non-shear controls stay authoritative.
    ///
    /// Optional canonical wall links measure impulse/work on the SAME augmented
    /// post-collision state. Native wall/open accuracy limits remain unchanged.
    /// Refusal leaves populations and `workspace.applied()` intact. Callbacks
    /// must not publish external material history: side effects cannot roll
    /// back. Source buffers allocate only when constructing the workspace;
    /// native stream scratch is reused after its first accepted allocation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_boundary_with_material_rheology_and_sgs<M: FluidTransport>(
        &self,
        grid: &mut BoundaryGrid3,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut BoundaryRheologySgsWorkspace3<M>,
        measurement: Option<(&[BoundaryLink3], [f64; 3])>,
        mut resolve: impl FnMut(MaterialImmersedPoint3) -> Result<M, MaterialError>,
    ) -> Result<BoundarySgsReport3, MaterialBoundarySgsError3> {
        self.admit_boundary_rheology_sgs(grid, gate, workspace)?;
        let pending = &mut workspace.pending;
        let model = workspace.model;
        let report = grid.step_with_sgs_sampled(
            gate,
            config,
            &mut workspace.flow,
            measurement,
            |cell, gradient| -> Result<_, MaterialBoundarySgsError3> {
                let (collision, sample) =
                    self.sample_boundary_rheology_sgs(model, cell, gradient, gate, &mut resolve)?;
                pending.push(sample);
                Ok(collision)
            },
        )?;
        workspace.accept(report);
        Ok(report)
    }

    /// Identical material/SGS arithmetic through the existing pooled collision
    /// and streaming owner. Source evaluation and wall-load reduction retain
    /// their canonical order independently of worker count. Pooled scheduling
    /// keeps its existing chunk-descriptor allocations and cancellation gates.
    #[allow(clippy::too_many_arguments)]
    pub fn step_boundary_pooled_with_material_rheology_and_sgs<
        M: FluidTransport,
        P: KernelRunner,
    >(
        &self,
        grid: &mut BoundaryGrid3,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut BoundaryRheologySgsWorkspace3<M>,
        measurement: Option<(&[BoundaryLink3], [f64; 3])>,
        mut resolve: impl FnMut(MaterialImmersedPoint3) -> Result<M, MaterialError>,
    ) -> Result<BoundarySgsReport3, MaterialBoundarySgsError3> {
        self.admit_boundary_rheology_sgs(grid, gate, workspace)?;
        let pending = &mut workspace.pending;
        let model = workspace.model;
        let report = grid.step_pooled_with_sgs_sampled(
            pool,
            gate,
            config,
            &mut workspace.flow,
            measurement,
            |cell, gradient| -> Result<_, MaterialBoundarySgsError3> {
                let (collision, sample) =
                    self.sample_boundary_rheology_sgs(model, cell, gradient, gate, &mut resolve)?;
                pending.push(sample);
                Ok(collision)
            },
        )?;
        workspace.accept(report);
        Ok(report)
    }

    fn boundary_rheology_scaling(&self) -> [u64; 3] {
        [
            self.cell_size_m().to_bits(),
            self.timestep_s().to_bits(),
            self.transport().density_kg_m3().to_bits(),
        ]
    }

    fn admit_boundary_rheology_sgs<M>(
        &self,
        grid: &BoundaryGrid3,
        gate: &CancelGate,
        workspace: &mut BoundaryRheologySgsWorkspace3<M>,
    ) -> Result<(), MaterialBoundarySgsError3> {
        checkpoint(gate)?;
        if workspace.dimensions != grid.dimensions()
            || workspace.scaling != self.boundary_rheology_scaling()
            || workspace.model != grid.collision_model()
        {
            return Err(MaterialBoundarySgsError3::Invalid(
                "boundary rheology dimensions, physical scaling or collision law changed",
            ));
        }
        workspace.pending.clear();
        Ok(())
    }

    fn sample_boundary_rheology_sgs<M: FluidTransport>(
        &self,
        model: CollisionModel3,
        cell: [usize; 3],
        gradient: [[f64; 3]; 3],
        gate: &CancelGate,
        resolve: &mut impl FnMut(MaterialImmersedPoint3) -> Result<M, MaterialError>,
    ) -> Result<(CollisionModel3, MaterialImmersedSample3<M>), MaterialBoundarySgsError3> {
        checkpoint(gate)?;
        let error = |source| MaterialBoundarySgsError3::Source { cell, source };
        let lattice_rate = strain_rate(gradient);
        let rate = lattice_rate / self.timestep_s();
        if !rate.is_finite() || (lattice_rate != 0.0 && rate == 0.0) {
            return Err(error(MaterialLbmError::Invalid(
                "boundary physical shear rate is unrepresentable",
            )));
        }
        let point = MaterialImmersedPoint3 {
            cell,
            shear_rate_per_s: rate,
        };
        let material = resolve(point);
        checkpoint(gate)?;
        let material = material.map_err(|source| error(source.into()))?;
        let transport = material.transport_state();
        checkpoint(gate)?;
        let transport = transport.map_err(|source| error(source.into()))?;
        if transport.density_kg_m3().to_bits() != self.transport().density_kg_m3().to_bits() {
            return Err(error(MaterialLbmError::Invalid(
                "boundary rheology requires unchanged reference density",
            )));
        }
        let tau = plan_momentum_scaling(
            transport.kinematic_viscosity_m2_s(),
            self.cell_size_m(),
            self.timestep_s(),
            0.0,
        )
        .map_err(error)?
        .0
        .tau;
        let collision = D3q19::with_shear_relaxation(model, tau).map_err(|source| {
            BoundarySgsError3::from(BoundaryStepError3::Collision { cell, source })
        })?;
        Ok((
            collision,
            MaterialImmersedSample3 {
                point,
                material,
                transport,
                shear_relaxation_time: tau,
            },
        ))
    }
}

fn checkpoint(gate: &CancelGate) -> Result<(), MaterialBoundarySgsError3> {
    if gate.is_requested() {
        Err(BoundarySgsError3::from(BoundaryStepError3::Cancelled).into())
    } else {
        Ok(())
    }
}
