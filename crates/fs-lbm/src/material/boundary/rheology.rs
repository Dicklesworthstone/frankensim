//! Source rheology in the original curved-wall D3Q19 step.
//!
//! The original Guo-frame stress estimator and safeguarded constitutive solver
//! determine local material viscosity. The existing material boundary adapter
//! still owns sampling admission, native collision/streaming and publication.

use super::{
    BoundaryGrid3, CurvedWall3, CurvedWallLink3, FluidTransport, FluidTransportState,
    MaterialBoundaryError3, MaterialCurvedWallStep3, MaterialError, MaterialLbmError,
    MaterialScaling, poll,
};
use crate::d3q19::{CollisionModel3, force_corrected_shear_rate3};
use crate::material::transport::rheology::{LocalRheologyError, solve_local};
use fs_exec::CancelGate;

pub use crate::material::transport::rheology::RheologySolve;

/// Actual incoming kinetic state supplied to the original material producer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialBoundaryRheologyPoint3 {
    pub cell: [usize; 3],
    /// sqrt(2 D:D) [1/s], in the physical held-force frame.
    pub shear_rate_per_s: f64,
}

/// Original source response retained at the accepted constitutive query.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialBoundaryRheology<M> {
    pub source: M,
    pub point: MaterialBoundaryRheologyPoint3,
    /// Actual source evaluations for this cell, including the accepted query.
    pub iterations: usize,
    /// Absolute tau residual of the retained source query.
    pub relaxation_defect: f64,
    transport: FluidTransportState,
}

impl<M> FluidTransport for MaterialBoundaryRheology<M> {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        // The source was read once in its local solve. The original boundary
        // adapter consumes this identical snapshot without another source read.
        Ok(self.transport)
    }
}

impl<R: FluidTransport + ?Sized> MaterialScaling<'_, R> {
    /// Resolve original material viscosity at self-consistent physical strain,
    /// then advance the same curved-wall momentum owner exactly once.
    ///
    /// Incoming populations and the ACTUAL held force supply A = gamma*tau.
    /// The shared bounded local solve closes tau = 1/2 + 3 nu(gamma) dt/dx^2
    /// together with gamma = A/(tau*dt). The constructor's shear relaxation is
    /// only an initial iterate. Every trial retains the reference density and
    /// uses the original source's supported rate domain; no viscosity floor,
    /// cap, constitutive fit or law-name dispatch is substituted.
    ///
    /// `resolve` owns temperature, pressure, composition and source selection.
    /// Each returned source's transport snapshot is read once per evaluation;
    /// the accepted source, physical query and convergence report are retained.
    /// Trials run in x-fastest fluid-cell order, skipping solids. A query must
    /// not advance or publish its caller's external thermal/material history.
    ///
    /// This is incoming-time generalized-Newtonian constitutive splitting.
    /// The native step preserves BGK/TRT/central-moment/reduced-cumulant controls,
    /// including D3Q19's bulk/shear tie, forcing and open-reservoir accounts.
    /// Geometry, wall velocity, origin and momentum/work retain lattice units.
    /// The physical scaling, low-Mach and scenario-speed obligations are those
    /// of `step_curved_walls_with_material`; next-time source admission belongs
    /// to the next call. No implicit time solve, thermal evolution, viscous heat,
    /// yield/viscoelastic history or general convergence guarantee is added.
    ///
    /// The complete stress field is admitted before material queries. Every
    /// local source evaluation polls cancellation, and every source/root/native
    /// refusal leaves accepted populations unchanged. Only the original native
    /// step publishes state; no fallible work follows that publication.
    #[allow(clippy::too_many_arguments)]
    pub fn step_curved_walls_with_material_rheology<M: FluidTransport>(
        &self,
        grid: &mut BoundaryGrid3,
        wall: &CurvedWall3,
        velocity: impl FnMut(CurvedWallLink3) -> [f64; 3],
        moment_origin: [f64; 3],
        max_cells: usize,
        solve: RheologySolve,
        gate: &CancelGate,
        mut resolve: impl FnMut(MaterialBoundaryRheologyPoint3) -> Result<M, MaterialError>,
    ) -> Result<MaterialCurvedWallStep3<MaterialBoundaryRheology<M>>, MaterialBoundaryError3> {
        poll(gate)?;
        solve.validate().map_err(|_| {
            MaterialBoundaryError3::Invalid("invalid local rheology solve controls")
        })?;
        let (numerators, [nx, ny, _]) = stress_field(grid, max_cells, gate)?;
        let initial_tau = match grid.collision_model() {
            CollisionModel3::Bgk { tau } | CollisionModel3::Trt { tau, .. } => tau,
            CollisionModel3::CentralMoment {
                second_order_rate, ..
            }
            | CollisionModel3::ReducedCumulant {
                second_order_rate, ..
            } => 1.0 / second_order_rate,
        };
        self.step_curved_walls_with_sampler(
            grid,
            wall,
            velocity,
            moment_origin,
            max_cells,
            gate,
            |cell| {
                let index = (cell[2] * ny + cell[1]) * nx + cell[0];
                let root = solve_local(
                    numerators[index],
                    self.timestep_s,
                    initial_tau,
                    solve,
                    || poll(gate),
                    |shear_rate_per_s| {
                        let point = MaterialBoundaryRheologyPoint3 {
                            cell,
                            shear_rate_per_s,
                        };
                        let error = |source| MaterialBoundaryError3::State { cell, source };
                        let source = resolve(point).map_err(|e| error(e.into()))?;
                        poll(gate)?;
                        let transport = source.transport_state().map_err(|e| error(e.into()))?;
                        poll(gate)?;
                        let tau = self.boundary_relaxation(cell, transport)?;
                        Ok(((source, transport), tau))
                    },
                )
                .map_err(|error| match error {
                    LocalRheologyError::Evaluation(error) => error,
                    LocalRheologyError::Invalid(what) => MaterialBoundaryError3::State {
                        cell,
                        source: MaterialLbmError::Invalid(what),
                    },
                    LocalRheologyError::Convergence {
                        iterations,
                        relaxation_defect,
                    } => MaterialBoundaryError3::RheologyConvergence {
                        cell,
                        iterations,
                        relaxation_defect,
                    },
                })?;
                Ok(MaterialBoundaryRheology {
                    source: root.value.0,
                    transport: root.value.1,
                    point: MaterialBoundaryRheologyPoint3 {
                        cell,
                        shear_rate_per_s: root.shear_rate_per_s,
                    },
                    iterations: root.iterations,
                    relaxation_defect: root.relaxation_defect,
                })
            },
        )
    }
}

fn stress_field(
    grid: &BoundaryGrid3,
    max_cells: usize,
    gate: &CancelGate,
) -> Result<(Vec<f64>, [usize; 3]), MaterialBoundaryError3> {
    let dims @ [nx, ny, nz] = grid.dimensions();
    let cells = nx
        .checked_mul(ny)
        .and_then(|plane| plane.checked_mul(nz))
        .filter(|&count| count <= max_cells)
        .ok_or(MaterialBoundaryError3::Invalid(
            "material boundary sampling exceeds its cell budget",
        ))?;
    let mut values = Vec::new();
    values.try_reserve_exact(cells).map_err(|_| {
        MaterialBoundaryError3::Invalid("material shear-stress allocation refused")
    })?;
    values.resize(cells, 0.0);
    for (index, numerator) in values.iter_mut().enumerate() {
        if index % 64 == 0 {
            poll(gate)?;
        }
        let cell @ [x, y, z] = [index % nx, (index / nx) % ny, index / (nx * ny)];
        if grid.is_solid(x, y, z) {
            continue;
        }
        let invalid = |what| MaterialBoundaryError3::State {
            cell,
            source: MaterialLbmError::Invalid(what),
        };
        let rho = grid.density(x, y, z);
        if !rho.is_finite() || rho <= 0.0 {
            return Err(invalid("rheology needs a positive finite native density"));
        }
        let u = grid.velocity(x, y, z);
        let force = grid.force_density(x, y, z);
        if !u.into_iter().chain(force).all(f64::is_finite) {
            return Err(invalid("rheology needs finite native velocity and held force"));
        }
        *numerator = force_corrected_shear_rate3(&grid.populations(x, y, z), rho, u, force, 1.0);
        if !numerator.is_finite() || *numerator < 0.0 {
            return Err(invalid("unrepresentable native non-equilibrium shear stress"));
        }
    }
    poll(gate)?;
    Ok((values, dims))
}
