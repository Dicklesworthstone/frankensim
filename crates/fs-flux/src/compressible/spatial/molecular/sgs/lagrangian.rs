//! Pathline-averaged Favre--Germano stress on the physical gas trajectory.
//!
//! The shared history owner transports the signed Germano numerator and
//! nonnegative denominator separately, with caller-authored relaxation time
//! in seconds. One fit samples the ACCEPTED state at the START of each complete
//! molecular/Euler/molecular interval. Its nonnegative C_s squared is held
//! through both molecular halves, all SSP/RKL stages and all coupling retries.
//! Each conservative face still evaluates its actual stage density and strain.
//!
//! Only acceptance of the complete physical interval commits that start-time
//! sample. A failed or cancelled interval leaves both the accepted flow and
//! history unchanged. Thus the next fit's elapsed memory time is the previous
//! ACTUAL accepted duration, including an Euler/CFL-shortened interval; neither
//! a requested duration nor a molecular pseudo-time advances the history.
//!
//! This is first-order pathline/source integration with a coefficient held
//! over each physical step. A second-order flow integrator does not raise that
//! temporal order. The explicit departure-distance admission still applies;
//! the caller must resolve the memory time and trajectories with physical
//! timesteps. The initial sample uses the shared local 27-cell fit.
//!
//! Only fully periodic Cartesian grids with at least three cells per axis are
//! admitted. Heat uses the authored fixed turbulent Prandtl number, not an
//! independently fitted scalar memory. Wall/inflow memory, backscatter energy,
//! isotropic SGS energy and shock/turbulence modeling are not implemented.

use super::*;
use fs_lbm::sgs3::dynamic::lagrangian::{LagrangianDynamicSgs, LagrangianSgsTrial};

/// Physical-time Lagrangian stress controls and fixed gradient-diffusion heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LagrangianSubgridConfig {
    /// Memory time is in seconds; the departure bound is in grid cells.
    pub model: LagrangianSgsConfig,
    /// Positive finite turbulent Prandtl number; no scalar-memory fit is implied.
    pub turbulent_prandtl: f64,
}

impl LagrangianSubgridConfig {
    pub fn validate(self) -> Result<(), MolecularError> {
        self.model.validate()?;
        if !positive(self.turbulent_prandtl) {
            return Err(MolecularError::Invalid(
                "positive finite turbulent Prandtl number required",
            ));
        }
        Ok(())
    }
}

pub(super) struct LagrangianTransport {
    pub(super) config: LagrangianSubgridConfig,
    // The coupled owner takes the filter while a borrowed history trial is
    // alive. This also marks the only interval in which held faces may run.
    filter: Option<LagrangianDynamicSgs>,
    filter_width: f64,
    density: Vec<f64>,
    velocity: Vec<[f64; 3]>,
    last_report: Option<DynamicSgsReport>,
}

impl LagrangianTransport {
    pub(super) fn filter_width(&self) -> f64 {
        self.filter_width
    }

    pub(super) fn admit_prepared_step(&self) -> Result<(), MolecularError> {
        if self.filter.is_some() {
            return Err(MolecularError::Invalid(
                "Lagrangian SGS requires a prepared coupled physical step",
            ));
        }
        Ok(())
    }
}

impl SubgridTransport {
    pub(in crate::compressible::spatial::molecular) fn required_lagrangian_bytes(
        flow: &Euler3d,
    ) -> Result<usize, MolecularError> {
        let cells = flow.states.len();
        // Shared filter and two history banks: 168 bytes/cell. Density,
        // physical velocity and held face coefficient add 40 bytes/cell.
        let fields = cells
            .checked_mul(5 * core::mem::size_of::<f64>())
            .ok_or(MolecularError::Invalid("Lagrangian SGS workspace overflow"))?;
        LagrangianDynamicSgs::required_bytes(cells)?
            .checked_add(fields)
            .ok_or(MolecularError::Invalid("Lagrangian SGS workspace overflow"))
    }

    pub(in crate::compressible::spatial::molecular) fn new_lagrangian(
        flow: &Euler3d,
        config: LagrangianSubgridConfig,
        maximum_bytes: usize,
    ) -> Result<Self, MolecularError> {
        config.validate()?;
        dynamic::admit_domain(flow)?;
        dynamic::admit_bytes(Self::required_lagrangian_bytes(flow)?, maximum_bytes)?;
        let n = flow.states.len();
        let filter = LagrangianDynamicSgs::new(flow.grid.shape, flow.grid.spacing, maximum_bytes)?;
        let filter_width = filter.filter_width();
        let density = allocated(n, 0.0)?;
        let velocity = allocated(n, [0.0; 3])?;
        let dynamic_coefficient_squared = allocated(n, 0.0)?;
        let bytes = filter.workspace_bytes()
            + density.capacity() * core::mem::size_of::<f64>()
            + velocity.capacity() * core::mem::size_of::<[f64; 3]>()
            + dynamic_coefficient_squared.capacity() * core::mem::size_of::<f64>();
        dynamic::admit_bytes(bytes, maximum_bytes)?;
        Ok(Self {
            closure: SubgridClosure::Lagrangian(LagrangianTransport {
                config,
                filter: Some(filter),
                filter_width,
                density,
                velocity,
                last_report: None,
            }),
            wall_models: [None; 6],
            dynamic_coefficient_squared,
        })
    }

    pub(in crate::compressible::spatial::molecular) fn lagrangian_config(
        &self,
    ) -> Option<LagrangianSubgridConfig> {
        match &self.closure {
            SubgridClosure::Lagrangian(dynamic) => Some(dynamic.config),
            _ => None,
        }
    }

    pub(in crate::compressible::spatial::molecular) fn lagrangian_sample_time(
        &self,
    ) -> Option<Time> {
        match &self.closure {
            SubgridClosure::Lagrangian(dynamic) => dynamic
                .filter
                .as_ref()
                .and_then(LagrangianDynamicSgs::accepted_time)
                .map(Time::new),
            _ => None,
        }
    }

    pub(in crate::compressible::spatial::molecular) fn lagrangian_report(
        &self,
    ) -> Option<DynamicSgsReport> {
        match &self.closure {
            SubgridClosure::Lagrangian(dynamic) => dynamic.last_report,
            _ => None,
        }
    }

    pub(in crate::compressible::spatial::molecular) fn take_lagrangian_filter(
        &mut self,
    ) -> Option<LagrangianDynamicSgs> {
        match &mut self.closure {
            SubgridClosure::Lagrangian(dynamic) => dynamic.filter.take(),
            _ => None,
        }
    }

    // The returned borrow belongs only to the detached history filter. The
    // face cache and complete coupled solver remain available during its trial.
    pub(in crate::compressible::spatial::molecular) fn prepare_lagrangian<'a>(
        &mut self,
        filter: &'a mut LagrangianDynamicSgs,
        flow: &Euler3d,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LagrangianSgsTrial<'a>, MolecularError> {
        let SubgridClosure::Lagrangian(dynamic) = &mut self.closure else {
            return Err(MolecularError::Invalid("missing Lagrangian SGS owner"));
        };
        dynamic.admit_prepared_step()?;
        for (cell, state) in flow.states.iter().enumerate() {
            if cell % 64 == 0 {
                poll(cancelled)?;
            }
            let primitive = state.primitive(&flow.gas, &mut *cancelled)?;
            dynamic.density[cell] = primitive.density;
            dynamic.velocity[cell] = primitive.velocity;
        }
        let trial = filter.prepare(
            dynamic.config.model,
            flow.time,
            &dynamic.density,
            &dynamic.velocity,
            &mut *cancelled,
        )?;
        for (cell, (coefficient, fitted)) in self
            .dynamic_coefficient_squared
            .iter_mut()
            .zip(trial.coefficients_squared())
            .enumerate()
        {
            if cell % 64 == 0 {
                poll(cancelled)?;
            }
            *coefficient = fitted;
        }
        poll(cancelled)?;
        Ok(trial)
    }

    /// Restore the detached filter on BOTH success and error. Only a committed
    /// complete interval supplies a new public report.
    pub(in crate::compressible::spatial::molecular) fn restore_lagrangian_filter(
        &mut self,
        filter: LagrangianDynamicSgs,
        committed_report: Option<DynamicSgsReport>,
    ) {
        let SubgridClosure::Lagrangian(dynamic) = &mut self.closure else {
            unreachable!("the coupled step retains its original SGS owner");
        };
        debug_assert!(dynamic.filter.is_none());
        dynamic.filter = Some(filter);
        if let Some(report) = committed_report {
            dynamic.last_report = Some(report);
        }
    }
}
