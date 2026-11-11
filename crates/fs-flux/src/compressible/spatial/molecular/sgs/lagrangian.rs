//! Pathline-averaged Favre--Germano transport on the physical gas trajectory.
//!
//! One fit samples the ACCEPTED state at the START of each complete
//! molecular/Euler/molecular interval. Fitted nonnegative coefficients are held
//! through both molecular halves, all SSP/RKL stages and all coupling retries.
//! Conservative faces still use their actual stage density and velocity gradient.
//!
//! Only acceptance of the complete interval commits its start-time sample.
//! Failed/cancelled steps and diagnostic observations preserve accepted history.
//! Thus memory uses the previous ACTUAL accepted duration, including shortening,
//! never molecular pseudo-time or the requested next duration.
//!
//! This is first-order pathline/source integration with coefficients held over
//! a physical step. The explicit departure-distance admission remains active.
//! Initial history uses a local 27-cell periodic or nine-cell tangent fit.
//! Memory times are in seconds. Resolved channel histories never wrap through
//! the physical wall pair; their two tangent axes remain periodic.
//!
//! Cartesian grids require at least three cells per axis. Heat may use fixed
//! turbulent Prandtl or an independent scalar history. No inflow memory,
//! backscatter energy, isotropic SGS energy, pressure-dilatation, shock closure
//! or DNS qualification is supplied.

use super::*;
use fs_lbm::sgs3::dynamic::lagrangian::scalar::{LagrangianDynamicScalar, LagrangianScalarTrial};
use fs_lbm::sgs3::dynamic::lagrangian::{LagrangianDynamicSgs, LagrangianSgsTrial};
use fs_lbm::sgs3::dynamic::scalar::DynamicScalarReport;

mod thermal;
pub use thermal::LagrangianThermalSubgridConfig;

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

#[derive(Clone, Copy)]
enum LagrangianConfig {
    FixedPrandtl(LagrangianSubgridConfig),
    Thermal(LagrangianThermalSubgridConfig),
}

// Moving this enum detaches the complete trajectory, including both histories
// for a scalar fit. No second outer stepping implementation is needed.
pub(in crate::compressible::spatial::molecular) enum LagrangianFilter {
    Momentum(LagrangianDynamicSgs),
    Thermal(LagrangianDynamicScalar),
}

impl LagrangianFilter {
    fn accepted_time(&self) -> Option<f64> {
        match self {
            Self::Momentum(filter) => filter.accepted_time(),
            Self::Thermal(filter) => filter.accepted_time(),
        }
    }

    pub(in crate::compressible::spatial::molecular) fn filter_width(&self) -> f64 {
        match self {
            Self::Momentum(filter) => filter.filter_width(),
            Self::Thermal(filter) => filter.filter_width(),
        }
    }

    fn workspace_bytes(&self) -> usize {
        match self {
            Self::Momentum(filter) => filter.workspace_bytes(),
            Self::Thermal(filter) => filter.workspace_bytes(),
        }
    }
}

/// The lifetime belongs exclusively to the detached owner, not to face scratch.
pub(in crate::compressible::spatial::molecular) enum LagrangianTrial<'a> {
    Momentum(LagrangianSgsTrial<'a>),
    Thermal(LagrangianScalarTrial<'a>),
}

impl LagrangianTrial<'_> {
    pub(in crate::compressible::spatial::molecular) fn report(&self) -> DynamicSgsReport {
        match self {
            Self::Momentum(trial) => trial.report(),
            Self::Thermal(trial) => trial.momentum_report(),
        }
    }

    pub(in crate::compressible::spatial::molecular) fn reports(
        &self,
    ) -> (DynamicSgsReport, Option<DynamicScalarReport>) {
        let heat = match self {
            Self::Momentum(_) => None,
            Self::Thermal(trial) => Some(trial.scalar_report()),
        };
        (self.report(), heat)
    }

    pub(in crate::compressible::spatial::molecular) fn coefficients_squared(
        &self,
    ) -> impl ExactSizeIterator<Item = f64> + '_ {
        let n = match self {
            Self::Momentum(trial) => trial.viscosities().len(),
            Self::Thermal(trial) => trial.viscosities().len(),
        };
        (0..n).map(|cell| {
            match self {
                Self::Momentum(trial) => trial.coefficient_squared(cell),
                Self::Thermal(trial) => trial.coefficient_squared(cell),
            }
            .expect("an admitted trial retains its coefficient shape")
        })
    }

    fn scalar_coefficient(&self, cell: usize) -> Option<f64> {
        match self {
            Self::Momentum(_) => None,
            Self::Thermal(trial) => trial.scalar_coefficient(cell),
        }
    }

    pub(in crate::compressible::spatial::molecular) fn commit(self) {
        match self {
            Self::Momentum(trial) => trial.commit(),
            Self::Thermal(trial) => trial.commit(),
        }
    }
}

pub(super) struct LagrangianTransport {
    config: LagrangianConfig,
    // None only while the coupled step or observation holds an admitted trial.
    filter: Option<LagrangianFilter>,
    filter_width: f64,
    normal_axis: Option<usize>,
    density: Vec<f64>,
    velocity: Vec<[f64; 3]>,
    // Empty on the fixed-Prandtl path; no extra fields are allocated there.
    temperature: Vec<f64>,
    heat_coefficient: Vec<f64>,
    last_report: Option<DynamicSgsReport>,
    last_heat_report: Option<DynamicScalarReport>,
}

impl LagrangianTransport {
    /// Owned capacities while the filter is attached (it is `None` only while a
    /// coupled step or observation holds the admitted trial, which then owns its bytes).
    pub(super) fn workspace_bytes(&self) -> usize {
        self.filter.as_ref().map_or(0, LagrangianFilter::workspace_bytes)
            + (self.density.capacity() + self.temperature.capacity() + self.heat_coefficient.capacity())
                * core::mem::size_of::<f64>()
            + self.velocity.capacity() * core::mem::size_of::<[f64; 3]>()
    }

    pub(super) fn filter_width(&self) -> f64 {
        self.filter_width
    }

    pub(super) fn normal_axis(&self) -> Option<usize> {
        self.normal_axis
    }

    pub(super) fn fixed_prandtl_config(&self) -> Option<LagrangianSubgridConfig> {
        match self.config {
            LagrangianConfig::FixedPrandtl(config) => Some(config),
            LagrangianConfig::Thermal(_) => None,
        }
    }

    pub(super) fn admit_prepared_step(&self) -> Result<(), MolecularError> {
        if self.filter.is_some() {
            return Err(MolecularError::Invalid(
                "Lagrangian SGS requires a prepared coupled trajectory sample",
            ));
        }
        Ok(())
    }
}

impl SubgridTransport {
    pub(in crate::compressible::spatial::molecular) fn required_lagrangian_bytes(
        flow: &Euler3d,
    ) -> Result<usize, MolecularError> {
        lagrangian_bytes(flow.states.len(), false)
    }

    pub(in crate::compressible::spatial::molecular) fn new_lagrangian(
        flow: &Euler3d,
        config: LagrangianSubgridConfig,
        maximum_bytes: usize,
    ) -> Result<Self, MolecularError> {
        config.validate()?;
        Self::new_lagrangian_owner(
            flow,
            LagrangianConfig::FixedPrandtl(config),
            None,
            maximum_bytes,
        )
    }

    pub(in crate::compressible::spatial::molecular) fn new_lagrangian_channel(
        flow: &Euler3d,
        config: LagrangianSubgridConfig,
        maximum_bytes: usize,
    ) -> Result<Self, MolecularError> {
        config.validate()?;
        let normal = channel_axis(&flow.grid, flow.boundaries)?;
        Self::new_lagrangian_owner(
            flow,
            LagrangianConfig::FixedPrandtl(config),
            Some(normal),
            maximum_bytes,
        )
    }

    pub(in crate::compressible::spatial::molecular) fn new_lagrangian_channel_thermal(
        flow: &Euler3d,
        config: LagrangianThermalSubgridConfig,
        maximum_bytes: usize,
    ) -> Result<Self, MolecularError> {
        config.validate()?;
        let normal = channel_axis(&flow.grid, flow.boundaries)?;
        Self::new_lagrangian_owner(
            flow,
            LagrangianConfig::Thermal(config),
            Some(normal),
            maximum_bytes,
        )
    }

    pub(in crate::compressible::spatial::molecular) fn lagrangian_channel_normal_axis(
        &self,
    ) -> Option<usize> {
        match &self.closure {
            SubgridClosure::Lagrangian(dynamic) => dynamic.normal_axis(),
            _ => None,
        }
    }

    fn new_lagrangian_owner(
        flow: &Euler3d,
        config: LagrangianConfig,
        normal_axis: Option<usize>,
        maximum_bytes: usize,
    ) -> Result<Self, MolecularError> {
        if normal_axis.is_none() {
            dynamic::admit_domain(flow)?;
        }
        let n = flow.states.len();
        let thermal = matches!(config, LagrangianConfig::Thermal(_));
        dynamic::admit_bytes(lagrangian_bytes(n, thermal)?, maximum_bytes)?;
        let filter = if thermal {
            LagrangianFilter::Thermal(match normal_axis {
                Some(normal) => LagrangianDynamicScalar::new_channel(
                    flow.grid.shape,
                    flow.grid.spacing,
                    normal,
                    maximum_bytes,
                )?,
                None => {
                    LagrangianDynamicScalar::new(flow.grid.shape, flow.grid.spacing, maximum_bytes)?
                }
            })
        } else {
            LagrangianFilter::Momentum(match normal_axis {
                Some(normal) => LagrangianDynamicSgs::new_channel(
                    flow.grid.shape,
                    flow.grid.spacing,
                    normal,
                    maximum_bytes,
                )?,
                None => {
                    LagrangianDynamicSgs::new(flow.grid.shape, flow.grid.spacing, maximum_bytes)?
                }
            })
        };
        let filter_width = filter.filter_width();
        let density = allocated(n, 0.0)?;
        let velocity = allocated(n, [0.0; 3])?;
        let dynamic_coefficient_squared = allocated(n, 0.0)?;
        let temperature = if thermal {
            allocated(n, 0.0)?
        } else {
            Vec::new()
        };
        let heat_coefficient = if thermal {
            allocated(n, 0.0)?
        } else {
            Vec::new()
        };
        let transport = LagrangianTransport {
            config,
            filter: Some(filter),
            filter_width,
            normal_axis,
            density,
            velocity,
            temperature,
            heat_coefficient,
            last_report: None,
            last_heat_report: None,
        };
        let bytes = transport.workspace_bytes()
            + dynamic_coefficient_squared.capacity() * core::mem::size_of::<f64>();
        dynamic::admit_bytes(bytes, maximum_bytes)?;
        Ok(Self {
            closure: SubgridClosure::Lagrangian(transport),
            wall_models: [None; 6],
            dynamic_coefficient_squared,
        })
    }

    pub(in crate::compressible::spatial::molecular) fn lagrangian_config(
        &self,
    ) -> Option<LagrangianSubgridConfig> {
        match &self.closure {
            SubgridClosure::Lagrangian(dynamic) => dynamic.fixed_prandtl_config(),
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
                .and_then(LagrangianFilter::accepted_time)
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
    ) -> Option<LagrangianFilter> {
        match &mut self.closure {
            SubgridClosure::Lagrangian(dynamic) => dynamic.filter.take(),
            _ => None,
        }
    }

    // Periodic low-level regression seam; physical consumers supply actual walls.
    #[cfg(test)]
    pub(in crate::compressible::spatial::molecular) fn prepare_lagrangian<'a>(
        &mut self,
        filter: &'a mut LagrangianFilter,
        flow: &Euler3d,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LagrangianTrial<'a>, MolecularError> {
        self.prepare_lagrangian_with_walls(
            filter,
            flow,
            &[MolecularWall3d::SlipAdiabatic; 6],
            cancelled,
        )
    }

    // The returned borrow belongs only to the detached history filter.
    pub(in crate::compressible::spatial::molecular) fn prepare_lagrangian_with_walls<'a>(
        &mut self,
        filter: &'a mut LagrangianFilter,
        flow: &Euler3d,
        walls: &[MolecularWall3d; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LagrangianTrial<'a>, MolecularError> {
        let SubgridClosure::Lagrangian(dynamic) = &mut self.closure else {
            return Err(MolecularError::Invalid("missing Lagrangian SGS owner"));
        };
        dynamic.admit_prepared_step()?;
        let wall_data = dynamic
            .normal_axis
            .map(|normal| channel_wall_values(normal, walls))
            .transpose()?;
        for (cell, state) in flow.states.iter().enumerate() {
            if cell % 64 == 0 {
                poll(cancelled)?;
            }
            let primitive = state.primitive(&flow.gas, &mut *cancelled)?;
            dynamic.density[cell] = primitive.density;
            dynamic.velocity[cell] = primitive.velocity;
            if let Some(temperature) = dynamic.temperature.get_mut(cell) {
                *temperature = primitive.temperature;
            }
        }
        let trial = match (dynamic.config, filter) {
            (LagrangianConfig::FixedPrandtl(config), LagrangianFilter::Momentum(filter)) => {
                LagrangianTrial::Momentum(match wall_data {
                    Some((wall_velocity, _)) => filter.prepare_with_walls(
                        config.model,
                        flow.time,
                        &dynamic.density,
                        &dynamic.velocity,
                        wall_velocity,
                        &mut *cancelled,
                    )?,
                    None => filter.prepare(
                        config.model,
                        flow.time,
                        &dynamic.density,
                        &dynamic.velocity,
                        &mut *cancelled,
                    )?,
                })
            }
            (LagrangianConfig::Thermal(config), LagrangianFilter::Thermal(filter)) => {
                LagrangianTrial::Thermal(match wall_data {
                    Some((wall_velocity, scalar_walls)) => filter.prepare_with_walls(
                        config.momentum,
                        config.heat,
                        flow.time,
                        &dynamic.density,
                        &dynamic.velocity,
                        wall_velocity,
                        &dynamic.temperature,
                        scalar_walls,
                        &mut *cancelled,
                    )?,
                    None => filter.prepare(
                        config.momentum,
                        config.heat,
                        flow.time,
                        &dynamic.density,
                        &dynamic.velocity,
                        &dynamic.temperature,
                        &mut *cancelled,
                    )?,
                })
            }
            _ => {
                return Err(MolecularError::Invalid(
                    "Lagrangian history and controls differ",
                ));
            }
        };
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
        for (cell, coefficient) in dynamic.heat_coefficient.iter_mut().enumerate() {
            if cell % 64 == 0 {
                poll(cancelled)?;
            }
            *coefficient = trial
                .scalar_coefficient(cell)
                .ok_or(MolecularError::Invalid(
                    "Lagrangian scalar coefficient shape differs",
                ))?;
        }
        poll(cancelled)?;
        Ok(trial)
    }

    /// Restore on success and error; only whole-step acceptance publishes reports.
    pub(in crate::compressible::spatial::molecular) fn restore_lagrangian_filter(
        &mut self,
        filter: LagrangianFilter,
        committed_reports: Option<(DynamicSgsReport, Option<DynamicScalarReport>)>,
    ) {
        let SubgridClosure::Lagrangian(dynamic) = &mut self.closure else {
            unreachable!("the coupled owner retains its original SGS choice");
        };
        debug_assert!(dynamic.filter.is_none());
        dynamic.filter = Some(filter);
        if let Some((momentum, heat)) = committed_reports {
            dynamic.last_report = Some(momentum);
            dynamic.last_heat_report = heat;
        }
    }
}

fn lagrangian_bytes(cells: usize, thermal: bool) -> Result<usize, MolecularError> {
    let fields = cells
        .checked_mul(if thermal { 56 } else { 40 })
        .ok_or(MolecularError::Invalid("Lagrangian SGS workspace overflow"))?;
    let filter = if thermal {
        LagrangianDynamicScalar::required_bytes(cells)?
    } else {
        LagrangianDynamicSgs::required_bytes(cells)?
    };
    filter
        .checked_add(fields)
        .ok_or(MolecularError::Invalid("Lagrangian SGS workspace overflow"))
}
