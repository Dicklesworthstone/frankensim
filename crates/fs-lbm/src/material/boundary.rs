//! Source viscosity on the original curved-wall collision and momentum owner.
//!
//! Temperature, composition and source-law evaluation remain caller-owned.
//! One held transport snapshot per fluid cell supplies this step's viscosity;
//! no thermal state, material interpolation or second flow solver is added.

use fs_exec::CancelGate;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};

use super::{MaterialLbmError, MaterialScaling, plan_momentum_scaling};
use crate::d3q19::{
    BoundaryGrid3, BoundaryStepError3, CurvedWall3, CurvedWallLink3, CurvedWallStep3,
};

/// Original source value and its exact applied momentum conversion.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialBoundarySample3<M> {
    /// Original grid coordinate. Samples are ordered x fastest, omitting solids.
    pub cell: [usize; 3],
    /// Complete original material response, including its source receipts.
    /// Borrowed concrete sources and `&dyn FluidTransport` are also supported.
    pub material: M,
    /// The single admitted numerical snapshot extracted from `material`.
    pub transport: FluidTransportState,
    /// Applied tau = 1/2 + (mu/rho0) dt/(c_s^2 dx^2).
    pub shear_relaxation_time: f64,
}

/// Actual native momentum exchange and the material values that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialCurvedWallStep3<M> {
    /// Original wall, body-force and periodic accounts, in lattice units.
    pub momentum: CurvedWallStep3,
    /// Applied source epoch, not an inferred endpoint temperature evaluation.
    pub applied: Vec<MaterialBoundarySample3<M>>,
}

/// Source/conversion failures remain distinct from original kinetic refusals.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialBoundaryError3 {
    /// The existing grid exceeds the requested sampling/allocation budget.
    Invalid(&'static str),
    /// Original source, reference-density or momentum-scaling refusal.
    State {
        cell: [usize; 3],
        source: MaterialLbmError,
    },
    /// Original collision, geometry, momentum-exchange or cancellation refusal.
    Flow(BoundaryStepError3),
}

impl core::fmt::Display for MaterialBoundaryError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "material boundary: {reason}"),
            Self::State { cell, source } => {
                write!(f, "material boundary at {cell:?}: {source}")
            }
            Self::Flow(error) => write!(f, "material boundary: {error}"),
        }
    }
}
impl std::error::Error for MaterialBoundaryError3 {}
impl From<BoundaryStepError3> for MaterialBoundaryError3 {
    fn from(error: BoundaryStepError3) -> Self {
        Self::Flow(error)
    }
}

impl<R: FluidTransport + ?Sized> MaterialScaling<'_, R> {
    /// Apply current source viscosity to the original moving curved-wall step.
    ///
    /// `resolve` supplies each fluid cell's original material response at its
    /// current caller-owned temperature/composition/state. It is called once
    /// per fluid cell, in x-fastest order, and never for a solid. It can return
    /// a borrowed already-resolved source, avoiding new property interpolation
    /// for a deliberately frozen field. Source side effects cannot be undone;
    /// callbacks must not publish external thermal/material state.
    ///
    /// The retained reference rho0, dx and dt convert the SAME response's
    /// nu = mu/rho0 to local relaxation. Every source must retain rho0 exactly:
    /// this incompressible step cannot remap mass or change inertia. Cp and k
    /// are admitted by the original transport owner but neither applied nor
    /// required to stay constant. No temperature evolution, viscous heating,
    /// source pressure reconstruction or shear-rate constitutive solve is
    /// inferred. Reevaluate sources at the caller's next accepted thermal
    /// epoch; holding them is an explicit coefficient split.
    ///
    /// Geometry, surface velocities, moment origin, held forces and the
    /// returned momentum accounts use the ORIGINAL GRID'S LATTICE UNITS.
    /// Existing `velocity_lattice` converts prescribed SI velocity components.
    /// The grid's density one must represent this scaling's reference rho0;
    /// the caller retains that unit convention and its scenario speed bound.
    /// Interpolation drift, wall/periodic accounts and native validity limits
    /// are unchanged. No wall work is relabeled as a thermal source.
    ///
    /// The per-cell field is transient; the constructor's collision family and
    /// non-shear controls remain untouched. All source/conversion work precedes
    /// native stepping, and every native endpoint/load/cancellation gate still
    /// precedes its one population swap. Failure publishes neither a new flow
    /// state nor a successful source report. No fallible work follows the swap.
    #[allow(clippy::too_many_arguments)]
    pub fn step_curved_walls_with_material<M: FluidTransport>(
        &self,
        grid: &mut BoundaryGrid3,
        wall: &CurvedWall3,
        velocity: impl FnMut(CurvedWallLink3) -> [f64; 3],
        moment_origin: [f64; 3],
        max_cells: usize,
        gate: &CancelGate,
        mut resolve: impl FnMut([usize; 3]) -> Result<M, MaterialError>,
    ) -> Result<MaterialCurvedWallStep3<M>, MaterialBoundaryError3> {
        poll(gate)?;
        let [nx, ny, nz] = grid.dimensions();
        let cells = nx
            .checked_mul(ny)
            .and_then(|plane| plane.checked_mul(nz))
            .filter(|&count| count <= max_cells)
            .ok_or(MaterialBoundaryError3::Invalid(
                "material boundary sampling exceeds its cell budget",
            ))?;
        let mut relaxation = Vec::new();
        let mut applied = Vec::new();
        relaxation.try_reserve_exact(cells).map_err(|_| {
            MaterialBoundaryError3::Invalid("material relaxation allocation refused")
        })?;
        applied
            .try_reserve_exact(cells)
            .map_err(|_| MaterialBoundaryError3::Invalid("material sample allocation refused"))?;
        // Solid padding is never sampled or used by the native collision.
        relaxation.resize(cells, 1.0);
        for (index, tau) in relaxation.iter_mut().enumerate() {
            if index % 64 == 0 {
                poll(gate)?;
            }
            let cell = [index % nx, (index / nx) % ny, index / (nx * ny)];
            if grid.is_solid(cell[0], cell[1], cell[2]) {
                continue;
            }
            let error = |source| MaterialBoundaryError3::State { cell, source };
            let material = resolve(cell).map_err(|source| error(source.into()))?;
            let transport = material
                .transport_state()
                .map_err(|source| error(source.into()))?;
            if transport.density_kg_m3().to_bits() != self.transport.density_kg_m3().to_bits() {
                return Err(error(MaterialLbmError::Invalid(
                    "curved material viscosity requires unchanged reference density",
                )));
            }
            // Momentum conversion consumes this snapshot, not the source a
            // second time. Thermal diffusion imposes no extra lattice bound.
            *tau = plan_momentum_scaling(
                transport.kinematic_viscosity_m2_s(),
                self.cell_size_m,
                self.timestep_s,
                0.0,
            )
            .map_err(error)?
            .0
            .tau;
            applied.push(MaterialBoundarySample3 {
                cell,
                material,
                transport,
                shear_relaxation_time: *tau,
            });
        }
        poll(gate)?;
        let momentum = grid.step_with_curved_walls_and_relaxation(
            wall,
            velocity,
            moment_origin,
            &relaxation,
            gate,
        )?;
        Ok(MaterialCurvedWallStep3 { momentum, applied })
    }
}

fn poll(gate: &CancelGate) -> Result<(), MaterialBoundaryError3> {
    if gate.is_requested() {
        Err(BoundaryStepError3::Cancelled.into())
    } else {
        Ok(())
    }
}
