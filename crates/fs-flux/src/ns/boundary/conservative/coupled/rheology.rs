//! Cellwise molecular rheology in the coupled endpoint-advector owner.

use super::super::{
    DensityMomentumControls, DensityMomentumError, DensityTransportCell, FluxSystem,
    MaterialFluxNumerics, MaterialNavierStokesBoundary, Result, Sums, allocate, invalid, poll,
};
use super::{CoupledDensityMomentumControls, CoupledDensityMomentumStep, CoupledViscosity};
use crate::ns::boundary::{StokesBoundaryError, StokesRheologyControls};
use crate::ns::transport::density::DensityTransportError;
use fs_material::MaterialError;

#[path = "rheology/ports.rs"]
mod ports;
pub use ports::CoupledDensityMomentumRheologyFlowControlStep;

/// Pure molecular source at (cell, centroid [m], transported state, shear [1/s]).
/// Returns positive finite dynamic viscosity [Pa s] or its original source refusal.
pub type CellRheology<'a, const N: usize> = dyn Fn(usize, [f64; 2], &DensityTransportCell<N>, f64) -> core::result::Result<f64, MaterialError>
    + 'a;

/// Exact molecular source sample retained at the accepted endpoint.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CoupledDensityMomentumRheologyCell {
    /// Objective native symmetric-strain invariant sqrt(2 D:D) [1/s].
    pub strain_rate_per_s: f64,
    /// Source response [Pa s], used unchanged at all quadrature sites in this cell.
    pub dynamic_viscosity_pa_s: f64,
    /// Physical molecular mu * gamma² [W/m³], without numerical dissipation.
    pub viscous_heating_w_m3: f64,
}

/// One accepted transport/velocity/source endpoint and its physical dissipation.
pub struct CoupledDensityMomentumRheologyStep<const N: usize> {
    /// Accepted conservative inventories, native flow and nonlinear evidence.
    pub step: CoupledDensityMomentumStep<N>,
    /// Original triangle order; reporting never queries the source again.
    pub cells: Vec<CoupledDensityMomentumRheologyCell>,
    /// Volume integral of molecular viscous heating at the declared depth [W].
    pub viscous_heating_w: f64,
}

impl FluxSystem<'_> {
    /// Couple conservative transported density to native shear-dependent stress.
    ///
    /// Each nonlinear assembly evaluates the pure CELLWISE molecular source once
    /// per triangle, at its centroid, with that trial's transported density and
    /// specific inventories and the native candidate sqrt(2 D:D) shear rate.
    /// The positive finite dynamic viscosity [Pa s] is used unchanged at every
    /// volume and face quadrature site belonging to that cell. Finer spatial
    /// variation remains the original pointwise viscosity entry's responsibility.
    /// Existing material laws can be called directly; original MaterialError
    /// values are preserved with the failing cell. No law or cutoff is invented.
    ///
    /// The original shared Picard owner re-transports immutable old inventories
    /// with each trial velocity and requires the unrelaxed correction and both
    /// reassembled endpoint residuals to converge. Density, compatible mass flux,
    /// source samples, momentum, mechanical pressure and reactions therefore
    /// belong to the SAME accepted endpoint. All velocity/natural-traction
    /// boundaries, CFL, donor, cancellation and refusal semantics of
    /// [`Self::density_momentum_coupled_bdf1_step`] apply.
    ///
    /// Returned mu*gamma² is physical molecular heating, distinct from temporal,
    /// upwind and SIP numerical losses. It is not deposited into a carried
    /// inventory and does not close a thermal feedback loop. Reporting uses the
    /// retained samples without additional source calls. One reusable cell buffer
    /// is checked against the per-owner transport storage allowance; this does
    /// not cap the native linear solver or total process memory.
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_coupled_bdf1_rheology_step<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &CellRheology<'_, N>,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: impl Into<CoupledDensityMomentumControls>,
        nonlinear: StokesRheologyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoupledDensityMomentumRheologyStep<N>> {
        let controls = controls.into();
        let (solved, response) = self.coupled_density_momentum_bdf1_viscosity(
            numerics,
            CoupledViscosity::Rheology(dynamic_viscosity),
            force_density,
            boundaries,
            previous_velocity,
            cells,
            inflow,
            timestep_s,
            controls,
            nonlinear,
            None,
            &mut cancelled,
        )?;
        self.finish_coupled_rheology(solved.flow, response, controls.momentum(), &mut cancelled)
    }

    pub(super) fn finish_coupled_rheology<const N: usize>(
        &self,
        step: CoupledDensityMomentumStep<N>,
        mut response: Vec<CoupledDensityMomentumRheologyCell>,
        controls: DensityMomentumControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CoupledDensityMomentumRheologyStep<N>> {
        let mut power = Sums::<1>::new();
        for (cell, output) in response.iter_mut().enumerate() {
            poll(cancelled)?;
            let rate = output.strain_rate_per_s;
            let local = (output.dynamic_viscosity_pa_s * rate) * rate;
            let integrated = (local * self.mesh.areas[cell]) * controls.depth_m;
            if !local.is_finite()
                || (rate > 0.0 && local == 0.0)
                || !integrated.is_finite()
                || (local > 0.0 && integrated == 0.0)
            {
                return Err(invalid("unrepresentable molecular viscous heating"));
            }
            output.viscous_heating_w_m3 = local;
            power.add([integrated]);
        }
        let viscous_heating_w = power.values()[0];
        if !viscous_heating_w.is_finite() {
            return Err(invalid(
                "unrepresentable integrated molecular viscous heating",
            ));
        }
        poll(cancelled)?;
        Ok(CoupledDensityMomentumRheologyStep {
            step,
            cells: response,
            viscous_heating_w,
        })
    }
}

pub(super) fn response_buffer(
    system: &FluxSystem<'_>,
    controls: DensityMomentumControls,
) -> Result<Vec<CoupledDensityMomentumRheologyCell>> {
    let required = system
        .mesh
        .tris
        .len()
        .checked_mul(size_of::<CoupledDensityMomentumRheologyCell>())
        .ok_or_else(|| invalid("rheology buffer size overflow"))?;
    let maximum = controls.maximum_transport_storage_bytes;
    if required > maximum {
        return Err(DensityTransportError::Storage { required, maximum }.into());
    }
    allocate(
        system.mesh.tris.len(),
        CoupledDensityMomentumRheologyCell::default(),
    )
}

pub(super) fn sample_cells<const N: usize>(
    system: &FluxSystem<'_>,
    source: &CellRheology<'_, N>,
    velocity: &[f64],
    transported: &[DensityTransportCell<N>],
    response: &mut [CoupledDensityMomentumRheologyCell],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    for (cell, output) in response.iter_mut().enumerate() {
        poll(cancelled)?;
        let rate = system
            .try_cell_kinematics(velocity, cell)
            .ok_or_else(|| invalid("unrepresentable native rheology kinematics"))?
            .shear_rate();
        if !rate.is_finite() {
            return Err(invalid("unrepresentable native rheology shear rate"));
        }
        let point = system.mesh.centroids[cell];
        let viscosity = source(cell, point, &transported[cell], rate)
            .map_err(|error| DensityMomentumError::Material { cell, error })?;
        if !viscosity.is_finite() || viscosity <= 0.0 {
            return Err(StokesBoundaryError::InvalidViscosity {
                cell,
                point,
                value: viscosity,
            }
            .into());
        }
        *output = CoupledDensityMomentumRheologyCell {
            strain_rate_per_s: rate,
            dynamic_viscosity_pa_s: viscosity,
            viscous_heating_w_m3: 0.0,
        };
    }
    Ok(())
}
