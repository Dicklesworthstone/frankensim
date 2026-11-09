//! Existing bounded density integrators selected inside one endpoint Picard owner.

use super::super::ssprk2::DensityMomentumSspRk2Controls;
use super::super::subcycle::DensityMomentumSubcycleControls;
use super::super::{
    DensityMomentumControls, DensityTransportCell, FluxSystem, MaterialFluxNumerics,
    MaterialNavierStokesBoundary, PreparedDensityMomentum, Result,
};
use crate::ns::transport::density::DensityTransportReport;
use crate::ns::transport::density::ssprk2::DensitySspRk2Report;
use crate::ns::transport::density::subcycle::DensitySubcycleReport;

/// Physical momentum controls and explicit density work for EVERY Picard trial.
/// Existing coupled entry points accept all three original control types through
/// From conversions. Ordinary DensityMomentumControls keep strict single-step
/// admission; selecting bounded transport is an explicit caller decision.
#[derive(Clone, Copy, Debug)]
pub enum CoupledDensityMomentumControls {
    /// One original Euler density step; the macro interval must meet its CFL cap.
    SingleStep(DensityMomentumControls),
    /// Native equal dyadic Euler substeps, with their accumulated mass flux.
    SubcycledEuler(DensityMomentumSubcycleControls),
    /// Native equal dyadic SSPRK2 subintervals, with their RK-weighted mass flux.
    /// Only density at fixed trial advection has second-order integration;
    /// momentum, external storage and endpoint-advector coupling remain BDF1.
    SspRk2(DensityMomentumSspRk2Controls),
}

impl From<DensityMomentumControls> for CoupledDensityMomentumControls {
    fn from(controls: DensityMomentumControls) -> Self {
        Self::SingleStep(controls)
    }
}
impl From<DensityMomentumSubcycleControls> for CoupledDensityMomentumControls {
    fn from(controls: DensityMomentumSubcycleControls) -> Self {
        Self::SubcycledEuler(controls)
    }
}
impl From<DensityMomentumSspRk2Controls> for CoupledDensityMomentumControls {
    fn from(controls: DensityMomentumSspRk2Controls) -> Self {
        Self::SspRk2(controls)
    }
}

impl CoupledDensityMomentumControls {
    pub(super) const fn momentum(self) -> DensityMomentumControls {
        match self {
            Self::SingleStep(controls) => controls,
            Self::SubcycledEuler(controls) => controls.momentum,
            Self::SspRk2(controls) => controls.momentum,
        }
    }
}

/// Original physical transport report and actual work at the accepted endpoint.
/// These counts describe one accepted macro interval, not all discarded Picard
/// assemblies. Each trial obeys its own explicit substep/storage allowance;
/// the original nonlinear iteration budget additionally bounds repeated work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CoupledDensityTransportReport<const N: usize> {
    /// One physical Euler step, without hidden subdivisions.
    SingleStep(DensityTransportReport<N>),
    /// Original accepted dyadic Euler count, duration, CFL and macro accounts.
    SubcycledEuler(DensitySubcycleReport<N>),
    /// Original accepted dyadic RK count, stage work, CFL and weighted accounts.
    SspRk2(DensitySspRk2Report<N>),
}

impl<'m> FluxSystem<'m> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare_coupled_density_momentum<'a, const N: usize>(
        &'a self,
        numerics: MaterialFluxNumerics,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        advector: &[f64],
        previous_velocity: &'a [f64],
        cells: &'a [DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: CoupledDensityMomentumControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(
        PreparedDensityMomentum<'a, 'm, N>,
        CoupledDensityTransportReport<N>,
    )> {
        match controls {
            CoupledDensityMomentumControls::SingleStep(controls) => {
                let prepared = self.prepare_density_momentum_bdf1(
                    numerics,
                    force_density,
                    boundaries,
                    advector,
                    previous_velocity,
                    cells,
                    inflow,
                    timestep_s,
                    controls,
                    cancelled,
                )?;
                let report = CoupledDensityTransportReport::SingleStep(prepared.transport);
                Ok((prepared, report))
            }
            CoupledDensityMomentumControls::SubcycledEuler(controls) => {
                let (prepared, report) = self.prepare_density_momentum_subcycled_bdf1(
                    numerics,
                    force_density,
                    boundaries,
                    advector,
                    previous_velocity,
                    cells,
                    inflow,
                    timestep_s,
                    controls,
                    cancelled,
                )?;
                Ok((
                    prepared,
                    CoupledDensityTransportReport::SubcycledEuler(report),
                ))
            }
            CoupledDensityMomentumControls::SspRk2(controls) => {
                let (prepared, report) = self.prepare_density_momentum_ssprk2_bdf1(
                    numerics,
                    force_density,
                    boundaries,
                    advector,
                    previous_velocity,
                    cells,
                    inflow,
                    timestep_s,
                    controls,
                    cancelled,
                )?;
                Ok((prepared, CoupledDensityTransportReport::SspRk2(report)))
            }
        }
    }
}
