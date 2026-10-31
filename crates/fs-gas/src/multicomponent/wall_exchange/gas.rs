//! Borrowed original calorics for the common implicit gas/wall source.
//!
//! No source table, composition, energy reference or inversion is reconstructed.

use super::{Cell, MixtureGas, WallExchangeError};
use crate::dynamics::CaloricDomain;

pub(crate) struct WallGasState {
    pub(crate) density: f64,
    pub(crate) velocity: f64,
    pub(crate) temperature: f64,
}

/// Fixed constituent densities during one thermal source interval. Energy and
/// Cv below are volumetric in the original gas owner's sensible-energy gauge.
pub(crate) trait WallGasCaloric {
    fn domain(&self) -> CaloricDomain;
    fn primitive(
        &self,
        momentum: f64,
        energy: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<WallGasState, WallExchangeError>;
    fn caloric(&self, temperature: f64) -> Result<(f64, f64), WallExchangeError>;
}

pub(super) struct MixtureCaloric<'a> {
    pub(super) gas: &'a MixtureGas,
    pub(super) partial_densities: &'a [f64],
}

impl WallGasCaloric for MixtureCaloric<'_> {
    fn domain(&self) -> CaloricDomain {
        self.gas.domain()
    }

    fn primitive(
        &self,
        momentum: f64,
        energy: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<WallGasState, WallExchangeError> {
        let state = self.gas.primitive(
            Cell {
                partial_densities: self.partial_densities,
                momentum,
                energy,
            },
            cancelled,
        )?;
        Ok(WallGasState {
            density: state.density,
            velocity: state.velocity,
            temperature: state.temperature,
        })
    }

    fn caloric(&self, temperature: f64) -> Result<(f64, f64), WallExchangeError> {
        Ok(self.gas.caloric(self.partial_densities, temperature)?)
    }
}
