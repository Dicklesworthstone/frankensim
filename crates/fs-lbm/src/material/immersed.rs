//! Original molecular rheology on the native moving immersed-body owner.
//!
//! The accepted physical velocity gradient supplies the incoming shear rate.
//! The same source response controls fluid collision AND reciprocal body
//! coupling. This is explicit in the incoming constitutive coefficients; no
//! pure-fluid nonequilibrium-stress formula is imposed on a mixed solid cell.

/// Simultaneous moving bodies reuse the same incoming source preparation.
pub mod multiple;

use super::{MaterialLbmError, MaterialScaling, plan_momentum_scaling};
use crate::d3q19::{
    BodyFrameSampler3, BodyFrameShape3, BoundaryGrid3, MovingPartialSaturationBody3,
    MovingPartialSaturationExchange3, PartialSaturationError3,
};
use crate::kinematics::{KinematicsError, KinematicsWorkspace3};
use fs_exec::CancelGate;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState, viscous::DeformationRate},
};

/// Conditions at the accepted incoming fluid time level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialImmersedPoint3 {
    pub cell: [usize; 3],
    /// sqrt(2 dev(D):dev(D)) [1/s], including all deviatoric normal and shear terms.
    pub shear_rate_per_s: f64,
}

/// The original source and exact snapshot used for both momentum owners.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialImmersedSample3<M> {
    pub point: MaterialImmersedPoint3,
    pub material: M,
    pub transport: FluidTransportState,
    pub shear_relaxation_time: f64,
}

/// A successful native transaction and its applied, incoming material epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialImmersedStep3<M, T> {
    /// Original fluid/body/force/wall accounts, in native lattice units.
    pub momentum: T,
    /// X-fastest source responses, omitting permanent solid cells.
    pub applied: Vec<MaterialImmersedSample3<M>>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum MaterialImmersedError3 {
    Invalid(&'static str),
    State {
        cell: [usize; 3],
        source: MaterialLbmError,
    },
    Kinematics(KinematicsError),
    Flow(PartialSaturationError3),
}

impl core::fmt::Display for MaterialImmersedError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material immersed rheology: {self:?}")
    }
}
impl std::error::Error for MaterialImmersedError3 {}
impl From<KinematicsError> for MaterialImmersedError3 {
    fn from(error: KinematicsError) -> Self {
        Self::Kinematics(error)
    }
}
impl From<PartialSaturationError3> for MaterialImmersedError3 {
    fn from(error: PartialSaturationError3) -> Self {
        Self::Flow(error)
    }
}

pub(super) struct PreparedImmersedRheology3<M> {
    pub(super) relaxation: Vec<f64>,
    pub(super) applied: Vec<MaterialImmersedSample3<M>>,
}

impl<R: FluidTransport + ?Sized> MaterialScaling<'_, R> {
    /// Advance one actual rigid body and the fluid with sourced local viscosity.
    ///
    /// The native physical half-force velocity supplies the existing centered
    /// interior/periodic and wall-aware/open-face gradient. Its deviatoric
    /// symmetric strain determines the query, including in retained interior fluid.
    /// Rigid rotation and isotropic dilation contribute no material shear.
    /// No turbulence closure is enabled.
    /// Sampling the mixed-cell velocity gradient avoids treating its blended
    /// nonequilibrium stress as though it came from pure-fluid collision.
    ///
    /// `resolve` owns temperature, composition and original constitutive data.
    /// Each non-solid cell is queried once, in x-fastest order; its returned
    /// `FluidTransport` is read once. Sources retain this scale's reference
    /// density exactly. The same mu/rho0 gives the local native collision tau,
    /// saturation weight and implicit rigid-body response. Source validity
    /// failures propagate; no viscosity clipping, fit or replacement is added.
    /// Cp/k remain with the source but impose no thermal lattice constraint.
    ///
    /// This is a frozen incoming-coefficient time split: the next invocation
    /// resamples the then-accepted flow. It is not an implicit flow/rheology
    /// solve, endpoint source admission or a heat integrator. Callbacks must
    /// not publish external material history during tentative sampling.
    /// Geometry, body state, loads and returned accounts remain in native
    /// lattice units; density one represents the retained reference rho0.
    /// Body inertia is ADDITIONAL to retained interior fluid, as in the native
    /// owner. No subgrid film, contact or displaced-fluid correction is added.
    ///
    /// The caller-owned kinematics workspace is reusable and byte-budgeted.
    /// `max_cells` bounds source/rate buffers before allocation. Native rigid
    /// drift, selected higher collision rates, optional endpoint force history
    /// and all endpoint/cancellation checks retain one fluid/body publication.
    /// Failure changes scratch only; no fallible work follows the native swap.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_body_with_material_rheology<M: FluidTransport>(
        &self,
        grid: &mut BoundaryGrid3,
        sampler: &mut BodyFrameSampler3,
        shape: &impl BodyFrameShape3,
        body: &mut MovingPartialSaturationBody3,
        external_impulse: [f64; 6],
        moment_origin: [f64; 3],
        workspace: &mut KinematicsWorkspace3,
        max_cells: usize,
        endpoint_force_scale: Option<f64>,
        gate: &CancelGate,
        resolve: impl FnMut(MaterialImmersedPoint3) -> Result<M, MaterialError>,
    ) -> Result<MaterialImmersedStep3<M, MovingPartialSaturationExchange3>, MaterialImmersedError3>
    {
        let prepared = self.prepare_immersed_rheology(grid, workspace, max_cells, gate, resolve)?;
        let momentum = grid.step_with_moving_partial_saturation_body_and_relaxation(
            sampler,
            shape,
            body,
            external_impulse,
            moment_origin,
            &prepared.relaxation,
            endpoint_force_scale,
            gate,
        )?;
        Ok(MaterialImmersedStep3 {
            momentum,
            applied: prepared.applied,
        })
    }

    pub(super) fn prepare_immersed_rheology<M: FluidTransport>(
        &self,
        grid: &BoundaryGrid3,
        workspace: &mut KinematicsWorkspace3,
        max_cells: usize,
        gate: &CancelGate,
        mut resolve: impl FnMut(MaterialImmersedPoint3) -> Result<M, MaterialError>,
    ) -> Result<PreparedImmersedRheology3<M>, MaterialImmersedError3> {
        poll(gate)?;
        let [nx, ny, nz] = grid.dimensions();
        let cells = nx
            .checked_mul(ny)
            .and_then(|n| n.checked_mul(nz))
            .filter(|&n| n <= max_cells)
            .ok_or(MaterialImmersedError3::Invalid(
                "immersed material sampling exceeds its cell budget",
            ))?;
        let field = grid.kinematics(gate, workspace)?;
        let mut relaxation = Vec::new();
        let mut applied = Vec::new();
        relaxation.try_reserve_exact(cells).map_err(|_| {
            MaterialImmersedError3::Invalid("immersed relaxation allocation refused")
        })?;
        applied
            .try_reserve_exact(cells)
            .map_err(|_| MaterialImmersedError3::Invalid("immersed source allocation refused"))?;
        relaxation.resize(cells, 1.0);
        for (index, kinematics) in field.iter().enumerate() {
            poll(gate)?;
            let Some(kinematics) = kinematics else {
                continue;
            };
            let cell = [index % nx, index / nx % ny, index / (nx * ny)];
            let state_error = |source| MaterialImmersedError3::State { cell, source };
            let lattice_rate = strain_rate(kinematics.gradient())
                .map_err(|source| state_error(source.into()))?;
            let rate = lattice_rate / self.timestep_s;
            if !rate.is_finite() || (lattice_rate != 0.0 && rate == 0.0) {
                return Err(state_error(MaterialLbmError::Invalid(
                    "immersed physical shear rate is unrepresentable",
                )));
            }
            let point = MaterialImmersedPoint3 {
                cell,
                shear_rate_per_s: rate,
            };
            let material = resolve(point);
            poll(gate)?;
            let material = material.map_err(|source| state_error(source.into()))?;
            let transport = material.transport_state();
            poll(gate)?;
            let transport = transport.map_err(|source| state_error(source.into()))?;
            if transport.density_kg_m3().to_bits() != self.transport.density_kg_m3().to_bits() {
                return Err(state_error(MaterialLbmError::Invalid(
                    "immersed rheology requires unchanged reference density",
                )));
            }
            let tau = plan_momentum_scaling(
                transport.kinematic_viscosity_m2_s(),
                self.cell_size_m,
                self.timestep_s,
                0.0,
            )
            .map_err(state_error)?
            .0
            .tau;
            relaxation[index] = tau;
            applied.push(MaterialImmersedSample3 {
                point,
                material,
                transport,
                shear_relaxation_time: tau,
            });
        }
        poll(gate)?;
        Ok(PreparedImmersedRheology3 {
            relaxation,
            applied,
        })
    }
}

// Share the material law's objective shear/volume split, then let the owning
// physical adapter convert this lattice rate by its retained timestep.
pub(crate) fn strain_rate(g: [[f64; 3]; 3]) -> Result<f64, MaterialError> {
    DeformationRate::from_velocity_gradient(g).map(|rate| rate.shear_rate_per_s())
}

fn poll(gate: &CancelGate) -> Result<(), MaterialImmersedError3> {
    if gate.is_requested() {
        Err(PartialSaturationError3::Cancelled.into())
    } else {
        Ok(())
    }
}
