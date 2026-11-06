//! Physical-time Lagrangian stress on the conservative mixture trajectory.
//!
//! One accepted START state supplies density and physical barycentric velocity
//! to the shared Favre--Germano history. The fitted nonnegative coefficient is
//! held through both molecular halves, SSP/RKL stages and whole-split retries.
//! The original face owner still evaluates current density, strain, composition,
//! source heat capacities and referenced partial enthalpies at every stage.
//! Thus common-Schmidt species fluxes retain zero net diffusive mass and carry
//! their own enthalpy; no separate viscous or chemical heat source is added.
//!
//! Only the complete nonreacting physical interval commits the history and its
//! report. Internal molecular clocks and rejected intervals cannot age it. The
//! next sample uses the previous ACTUAL accepted duration, including shortening.
//! This is first-order pathline/coefficient integration with an explicit memory
//! time in seconds. Fixed Prandtl/Schmidt numbers are authored model parameters,
//! not dynamically fitted scalar histories. Fully periodic Cartesian grids only;
//! no wall/inflow history, reaction-rate closure, SGS energy or shock model.

use super::*;
use fs_lbm::sgs3::dynamic::lagrangian::{LagrangianDynamicSgs, LagrangianSgsTrial};

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

pub(super) struct LagrangianTransport {
    pub(super) config: MixtureLagrangianSubgridConfig,
    // The complete coupled owner temporarily detaches this filter while its
    // borrowed history trial exists. Held face transport is allowed only then.
    filter: Option<LagrangianDynamicSgs>,
    // Detached storage still belongs to the transaction and its byte budget.
    filter_bytes: usize,
    filter_width: f64,
    density: Vec<f64>,
    velocity: Vec<[f64; 3]>,
    held_nonzero: bool,
    last_report: Option<DynamicSgsReport>,
}

impl LagrangianTransport {
    pub(super) fn workspace_bytes(&self) -> usize {
        self.filter_bytes
            + self.density.capacity() * size_of::<f64>()
            + self.velocity.capacity() * size_of::<[f64; 3]>()
    }

    pub(super) fn filter_width(&self) -> f64 {
        self.filter_width
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
        dynamic::admit_domain(flow)?;
        let cells = flow.cell_count();
        let bytes = LagrangianDynamicSgs::required_bytes(cells).map_err(dynamic::dynamic_error)?;
        let filter = LagrangianDynamicSgs::new(flow.grid.shape, flow.grid.spacing, bytes)
            .map_err(dynamic::dynamic_error)?;
        let filter_width = filter.filter_width();
        let filter_bytes = filter.workspace_bytes();
        Ok(Self {
            closure: ViscosityClosure::Lagrangian(LagrangianTransport {
                config,
                filter: Some(filter),
                filter_bytes,
                filter_width,
                density: allocated(cells, 0.0)?,
                velocity: allocated(cells, [0.0; 3])?,
                held_nonzero: false,
                last_report: None,
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
            ViscosityClosure::Lagrangian(dynamic) => Some(dynamic.config),
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
                .and_then(LagrangianDynamicSgs::accepted_time)
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
    ) -> Option<LagrangianDynamicSgs> {
        match &mut self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.filter.take(),
            _ => None,
        }
    }

    // Only the detached filter is borrowed by the result. The complete mixture
    // owner remains available for the original conservative step transaction.
    pub(in crate::compressible::spatial::mixture::molecular) fn prepare_lagrangian<'a>(
        &mut self,
        filter: &'a mut LagrangianDynamicSgs,
        flow: &MixtureEuler3d,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LagrangianSgsTrial<'a>, DissipationError> {
        let ViscosityClosure::Lagrangian(dynamic) = &mut self.closure else {
            return Err(DissipationError::Invalid(
                "missing Lagrangian mixture owner",
            ));
        };
        dynamic.admit_prepared_step()?;
        let species = flow.gas.species_count();
        for (index, row) in flow.states.chunks_exact(species + 4).enumerate() {
            if index % 64 == 0 {
                poll(cancelled)?;
            }
            let primitive = cell(row, species).primitive(&flow.gas, &mut *cancelled)?;
            dynamic.density[index] = primitive.density;
            dynamic.velocity[index] = primitive.velocity;
        }
        let trial = filter
            .prepare(
                dynamic.config.transport.model,
                flow.time,
                &dynamic.density,
                &dynamic.velocity,
                &mut *cancelled,
            )
            .map_err(dynamic::dynamic_error)?;
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
        poll(cancelled)?;
        dynamic.held_nonzero = trial.report().maximum_coefficient_squared != 0.0;
        Ok(trial)
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn restore_lagrangian_filter(
        &mut self,
        filter: LagrangianDynamicSgs,
        committed_report: Option<DynamicSgsReport>,
    ) {
        let ViscosityClosure::Lagrangian(dynamic) = &mut self.closure else {
            unreachable!("a coupled step retains its original mixture SGS owner");
        };
        debug_assert!(dynamic.filter.is_none());
        debug_assert_eq!(filter.workspace_bytes(), dynamic.filter_bytes);
        dynamic.filter = Some(filter);
        if let Some(report) = committed_report {
            dynamic.last_report = Some(report);
        }
    }
}
