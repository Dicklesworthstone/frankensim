//! Physical-time Lagrangian stress on the conservative mixture trajectory.
//!
//! One accepted START state supplies density and physical barycentric velocity
//! to the shared Favre--Germano history. The fitted nonnegative coefficient is
//! held through both molecular halves, SSP/RKL stages and whole-split retries.
//! The original face owner still evaluates current density, strain, composition,
//! source heat capacities and referenced partial enthalpies at every stage.
//! Both fixed-Schmidt and independently fitted species fluxes retain zero net
//! diffusive mass and carry their own referenced enthalpy; no separate viscous
//! or chemical heat source is added.
//!
//! Only the complete physical interval commits the history and its
//! report. Internal molecular clocks and rejected intervals cannot age it. The
//! next sample uses the previous ACTUAL accepted duration, including shortening.
//! This is first-order pathline/coefficient integration with an explicit memory
//! time in seconds. Fixed Prandtl/Schmidt controls remain an explicit option.
//! Independent heat/common-species histories instead use their own Germano
//! contractions and memory times, sharing the same momentum trajectory and
//! accepted sample clock. Cartesian periodic or resolved no-slip channel grids
//! use the same history owner; only tangent directions wrap at a wall. This
//! adds no inflow history, reaction-rate closure, SGS energy or shock model.

use super::*;
use fs_lbm::sgs3::dynamic::lagrangian::scalar::components::{
    LagrangianComponentsTrial, LagrangianDynamicComponents,
};
use fs_lbm::sgs3::dynamic::lagrangian::{LagrangianDynamicSgs, LagrangianSgsTrial};

mod channel;
mod species;
pub use species::{
    LagrangianComponentsConfig, LagrangianScalarConfig, MixtureLagrangianSpeciesSubgridConfig,
};

/// Pathline-averaged stress plus fixed gradient-diffusion heat/species controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureLagrangianSubgridConfig {
    /// Memory time is in seconds; physical fields use the mixture owner's SI units.
    pub transport: LagrangianSubgridConfig,
    /// Positive finite common turbulent Schmidt number for all species.
    pub turbulent_schmidt: f64,
}

impl MixtureLagrangianSubgridConfig {
    pub fn validate(self) -> Result<(), DissipationError> {
        self.transport.validate().map_err(|_| {
            DissipationError::Invalid("invalid Lagrangian mixture stress/heat controls")
        })?;
        if !positive(self.turbulent_schmidt) {
            return Err(DissipationError::Invalid(
                "positive finite turbulent Schmidt number required",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
enum LagrangianConfig {
    FixedTransport(MixtureLagrangianSubgridConfig),
    Species(MixtureLagrangianSpeciesSubgridConfig),
}

// The complete physical split detaches exactly one trajectory owner. The
// independent variant owns all three histories and cannot advance one alone.
pub(in crate::compressible::spatial::mixture) enum LagrangianFilter {
    Momentum(LagrangianDynamicSgs),
    Species(LagrangianDynamicComponents),
}

impl LagrangianFilter {
    fn accepted_time(&self) -> Option<f64> {
        match self {
            Self::Momentum(filter) => filter.accepted_time(),
            Self::Species(filter) => filter.accepted_time(),
        }
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn filter_width(&self) -> f64 {
        match self {
            Self::Momentum(filter) => filter.filter_width(),
            Self::Species(filter) => filter.filter_width(),
        }
    }

    fn workspace_bytes(&self) -> usize {
        match self {
            Self::Momentum(filter) => filter.workspace_bytes(),
            Self::Species(filter) => filter.workspace_bytes(),
        }
    }
}

pub(in crate::compressible::spatial::mixture) type LagrangianReports = (
    DynamicSgsReport,
    Option<DynamicScalarReport>,
    Option<DynamicScalarReport>,
);

/// The borrow belongs only to the detached history, leaving face scratch free.
pub(in crate::compressible::spatial::mixture) enum LagrangianTrial<'a> {
    Momentum(LagrangianSgsTrial<'a>),
    Species(LagrangianComponentsTrial<'a>),
}

impl LagrangianTrial<'_> {
    pub(in crate::compressible::spatial::mixture::molecular) fn report(&self) -> DynamicSgsReport {
        match self {
            Self::Momentum(trial) => trial.report(),
            Self::Species(trial) => trial.momentum_report(),
        }
    }

    pub(in crate::compressible::spatial::mixture) fn reports(&self) -> LagrangianReports {
        match self {
            Self::Momentum(_) => (self.report(), None, None),
            Self::Species(trial) => (
                self.report(),
                Some(trial.scalar_report()),
                Some(trial.components_report()),
            ),
        }
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn coefficients_squared(
        &self,
    ) -> impl ExactSizeIterator<Item = f64> + '_ {
        let cells = match self {
            Self::Momentum(trial) => trial.viscosities().len(),
            Self::Species(trial) => trial.viscosities().len(),
        };
        (0..cells).map(|cell| {
            match self {
                Self::Momentum(trial) => trial.coefficient_squared(cell),
                Self::Species(trial) => trial.coefficient_squared(cell),
            }
            .expect("an admitted Lagrangian trial retains its coefficient shape")
        })
    }

    fn scalar_coefficient(&self, cell: usize) -> Option<f64> {
        match self {
            Self::Momentum(_) => None,
            Self::Species(trial) => trial.scalar_coefficient(cell),
        }
    }

    fn component_coefficient(&self, cell: usize) -> Option<f64> {
        match self {
            Self::Momentum(_) => None,
            Self::Species(trial) => trial.component_coefficient(cell),
        }
    }

    pub(in crate::compressible::spatial::mixture) fn commit(self) {
        match self {
            Self::Momentum(trial) => trial.commit(),
            Self::Species(trial) => trial.commit(),
        }
    }
}

pub(super) struct LagrangianTransport {
    config: LagrangianConfig,
    // The complete coupled owner temporarily detaches this filter while its
    // borrowed history trial exists. Held face transport is allowed only then.
    filter: Option<LagrangianFilter>,
    // Detached storage still belongs to the transaction and its byte budget.
    filter_bytes: usize,
    filter_width: f64,
    normal_axis: Option<usize>,
    density: Vec<f64>,
    velocity: Vec<[f64; 3]>,
    // Empty for fixed Prandtl/Schmidt: its original capacity is unchanged.
    temperature: Vec<f64>,
    heat_coefficient: Vec<f64>,
    species_coefficient: Vec<f64>,
    held_nonzero: bool,
    last_report: Option<DynamicSgsReport>,
    last_heat_report: Option<DynamicScalarReport>,
    last_species_report: Option<DynamicScalarReport>,
}

impl LagrangianTransport {
    pub(super) fn workspace_bytes(&self) -> usize {
        self.filter_bytes
            + self.density.capacity() * size_of::<f64>()
            + self.velocity.capacity() * size_of::<[f64; 3]>()
            + (self.temperature.capacity()
                + self.heat_coefficient.capacity()
                + self.species_coefficient.capacity())
                * size_of::<f64>()
    }

    pub(super) fn filter_width(&self) -> f64 {
        self.filter_width
    }

    pub(super) fn fixed_transport_config(&self) -> Option<MixtureLagrangianSubgridConfig> {
        match self.config {
            LagrangianConfig::FixedTransport(config) => Some(config),
            LagrangianConfig::Species(_) => None,
        }
    }

    pub(super) fn transport_controls(&self) -> (Option<f64>, Option<f64>, f64) {
        match self.config {
            LagrangianConfig::FixedTransport(config) => (
                Some(config.transport.turbulent_prandtl),
                Some(config.turbulent_schmidt),
                config.transport.model.maximum_eddy_to_molecular,
            ),
            LagrangianConfig::Species(config) => {
                (None, None, config.momentum.maximum_eddy_to_molecular)
            }
        }
    }

    pub(super) fn admit_prepared_step(&self) -> Result<bool, DissipationError> {
        if self.filter.is_some() {
            return Err(DissipationError::Invalid(
                "Lagrangian mixture SGS requires a prepared coupled physical step",
            ));
        }
        Ok(self.held_nonzero)
    }
}

impl SubgridTransport {
    pub(in crate::compressible::spatial::mixture::molecular) fn required_lagrangian_bytes(
        flow: &MixtureEuler3d,
    ) -> Result<usize, DissipationError> {
        // Shared history/filter: 168 bytes/cell. Density/velocity staging,
        // held coefficient and the existing stage-local rho*Cv add 48.
        LagrangianDynamicSgs::required_bytes(flow.cell_count())
            .map_err(dynamic::dynamic_error)?
            .checked_add(flow.cell_count().checked_mul(6 * size_of::<f64>()).ok_or(
                DissipationError::Invalid("Lagrangian mixture workspace size overflow"),
            )?)
            .ok_or(DissipationError::Invalid(
                "Lagrangian mixture workspace size overflow",
            ))
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn new_lagrangian(
        flow: &MixtureEuler3d,
        config: MixtureLagrangianSubgridConfig,
    ) -> Result<Self, DissipationError> {
        config.validate()?;
        Self::new_lagrangian_owner(flow, LagrangianConfig::FixedTransport(config), None)
    }

    fn new_lagrangian_owner(
        flow: &MixtureEuler3d,
        config: LagrangianConfig,
        normal_axis: Option<usize>,
    ) -> Result<Self, DissipationError> {
        if normal_axis.is_none() {
            dynamic::admit_domain(flow)?;
        }
        let cells = flow.cell_count();
        let independent = matches!(config, LagrangianConfig::Species(_));
        let filter = if independent {
            LagrangianFilter::Species(
                match normal_axis {
                    Some(normal) => LagrangianDynamicComponents::new_channel(
                        flow.grid.shape,
                        flow.grid.spacing,
                        normal,
                        flow.gas.species_count(),
                        LagrangianDynamicComponents::required_bytes(cells)
                            .map_err(dynamic::dynamic_error)?,
                    ),
                    None => LagrangianDynamicComponents::new(
                        flow.grid.shape,
                        flow.grid.spacing,
                        flow.gas.species_count(),
                        LagrangianDynamicComponents::required_bytes(cells)
                            .map_err(dynamic::dynamic_error)?,
                    ),
                }
                .map_err(dynamic::dynamic_error)?,
            )
        } else {
            LagrangianFilter::Momentum(
                match normal_axis {
                    Some(normal) => LagrangianDynamicSgs::new_channel(
                        flow.grid.shape,
                        flow.grid.spacing,
                        normal,
                        LagrangianDynamicSgs::required_bytes(cells)
                            .map_err(dynamic::dynamic_error)?,
                    ),
                    None => LagrangianDynamicSgs::new(
                        flow.grid.shape,
                        flow.grid.spacing,
                        LagrangianDynamicSgs::required_bytes(cells)
                            .map_err(dynamic::dynamic_error)?,
                    ),
                }
                .map_err(dynamic::dynamic_error)?,
            )
        };
        let filter_width = filter.filter_width();
        let filter_bytes = filter.workspace_bytes();
        let scalar_cells = if independent { cells } else { 0 };
        Ok(Self {
            closure: ViscosityClosure::Lagrangian(LagrangianTransport {
                config,
                filter: Some(filter),
                filter_bytes,
                filter_width,
                normal_axis,
                density: allocated(cells, 0.0)?,
                velocity: allocated(cells, [0.0; 3])?,
                temperature: allocated(scalar_cells, 0.0)?,
                heat_coefficient: allocated(scalar_cells, 0.0)?,
                species_coefficient: allocated(scalar_cells, 0.0)?,
                held_nonzero: false,
                last_report: None,
                last_heat_report: None,
                last_species_report: None,
            }),
            wall_models: [None; 6],
            dynamic_coefficient_squared: allocated(cells, 0.0)?,
            heat_capacity_v_density: allocated(cells, 0.0)?,
        })
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn lagrangian_config(
        &self,
    ) -> Option<MixtureLagrangianSubgridConfig> {
        match &self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.fixed_transport_config(),
            _ => None,
        }
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn lagrangian_sample_time(
        &self,
    ) -> Option<Time> {
        match &self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic
                .filter
                .as_ref()
                .and_then(LagrangianFilter::accepted_time)
                .map(Time::new),
            _ => None,
        }
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn lagrangian_report(
        &self,
    ) -> Option<DynamicSgsReport> {
        match &self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.last_report,
            _ => None,
        }
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn take_lagrangian_filter(
        &mut self,
    ) -> Option<LagrangianFilter> {
        match &mut self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.filter.take(),
            _ => None,
        }
    }

    // Only the detached filter is borrowed by the result. The complete mixture
    // owner remains available for the original conservative step transaction.
    #[cfg(test)]
    pub(in crate::compressible::spatial::mixture::molecular) fn prepare_lagrangian<'a>(
        &mut self,
        filter: &'a mut LagrangianFilter,
        flow: &MixtureEuler3d,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LagrangianTrial<'a>, DissipationError> {
        self.prepare_lagrangian_with_walls(
            filter,
            flow,
            &[MolecularWall3d::SlipAdiabatic; 6],
            cancelled,
        )
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn prepare_lagrangian_with_walls<
        'a,
    >(
        &mut self,
        filter: &'a mut LagrangianFilter,
        flow: &MixtureEuler3d,
        walls: &[MolecularWall3d; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LagrangianTrial<'a>, DissipationError> {
        let ViscosityClosure::Lagrangian(dynamic) = &mut self.closure else {
            return Err(DissipationError::Invalid(
                "missing Lagrangian mixture owner",
            ));
        };
        dynamic.admit_prepared_step()?;
        let wall_data = dynamic
            .normal_axis
            .map(|normal| {
                crate::compressible::spatial::molecular::sgs::channel_wall_values(normal, walls)
                    .map_err(dynamic::dynamic_error)
            })
            .transpose()?;
        let species = flow.gas.species_count();
        for (index, row) in flow.states.chunks_exact(species + 4).enumerate() {
            if index % 64 == 0 {
                poll(cancelled)?;
            }
            let primitive = cell(row, species).primitive(&flow.gas, &mut *cancelled)?;
            dynamic.density[index] = primitive.density;
            dynamic.velocity[index] = primitive.velocity;
            if let Some(temperature) = dynamic.temperature.get_mut(index) {
                *temperature = primitive.temperature;
            }
        }
        let trial = match (dynamic.config, filter) {
            (LagrangianConfig::FixedTransport(config), LagrangianFilter::Momentum(filter)) => {
                let trial = match wall_data {
                    Some((wall_velocity, _)) => filter.prepare_with_walls(
                        config.transport.model,
                        flow.time,
                        &dynamic.density,
                        &dynamic.velocity,
                        wall_velocity,
                        &mut *cancelled,
                    ),
                    None => filter.prepare(
                        config.transport.model,
                        flow.time,
                        &dynamic.density,
                        &dynamic.velocity,
                        &mut *cancelled,
                    ),
                }
                .map_err(dynamic::dynamic_error)?;
                LagrangianTrial::Momentum(trial)
            }
            (LagrangianConfig::Species(config), LagrangianFilter::Species(filter)) => {
                let component = |component, index| {
                    flow.states[index * (species + 4) + component] / dynamic.density[index]
                };
                let trial = match wall_data {
                    Some((wall_velocity, scalar_walls)) => filter.prepare_with_walls(
                        config.momentum,
                        config.heat,
                        config.species,
                        flow.time,
                        &dynamic.density,
                        &dynamic.velocity,
                        wall_velocity,
                        &dynamic.temperature,
                        scalar_walls,
                        component,
                        &mut *cancelled,
                    ),
                    None => filter.prepare(
                        config.momentum,
                        config.heat,
                        config.species,
                        flow.time,
                        &dynamic.density,
                        &dynamic.velocity,
                        &dynamic.temperature,
                        component,
                        &mut *cancelled,
                    ),
                }
                .map_err(dynamic::dynamic_error)?;
                if trial.components_report().maximum_diffusivity
                    > config.maximum_species_diffusivity_m2_s
                {
                    return Err(DissipationError::Invalid(
                        "Lagrangian species diffusivity exceeds its physical admission bound",
                    ));
                }
                LagrangianTrial::Species(trial)
            }
            _ => {
                return Err(DissipationError::Invalid(
                    "Lagrangian mixture history and controls differ",
                ));
            }
        };
        for (index, (coefficient, fitted)) in self
            .dynamic_coefficient_squared
            .iter_mut()
            .zip(trial.coefficients_squared())
            .enumerate()
        {
            if index % 64 == 0 {
                poll(cancelled)?;
            }
            *coefficient = fitted;
        }
        for (index, (heat, species)) in dynamic
            .heat_coefficient
            .iter_mut()
            .zip(&mut dynamic.species_coefficient)
            .enumerate()
        {
            if index % 64 == 0 {
                poll(cancelled)?;
            }
            *heat = trial
                .scalar_coefficient(index)
                .ok_or(DissipationError::Invalid(
                    "Lagrangian mixture heat coefficient shape differs",
                ))?;
            *species = trial
                .component_coefficient(index)
                .ok_or(DissipationError::Invalid(
                    "Lagrangian mixture species coefficient shape differs",
                ))?;
        }
        poll(cancelled)?;
        let (momentum, heat, species) = trial.reports();
        dynamic.held_nonzero = momentum.maximum_coefficient_squared != 0.0
            || heat.is_some_and(|report| report.maximum_coefficient != 0.0)
            || species.is_some_and(|report| report.maximum_coefficient != 0.0);
        Ok(trial)
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn restore_lagrangian_filter(
        &mut self,
        filter: LagrangianFilter,
        committed_reports: Option<LagrangianReports>,
    ) {
        let ViscosityClosure::Lagrangian(dynamic) = &mut self.closure else {
            unreachable!("a coupled step retains its original mixture SGS owner");
        };
        debug_assert!(dynamic.filter.is_none());
        debug_assert_eq!(filter.workspace_bytes(), dynamic.filter_bytes);
        dynamic.filter = Some(filter);
        if let Some((momentum, heat, species)) = committed_reports {
            dynamic.last_report = Some(momentum);
            dynamic.last_heat_report = heat;
            dynamic.last_species_report = species;
        }
    }
}
